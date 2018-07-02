//
// This file is part of the libWetHair open source project
//
// The code is licensed solely for academic and non-commercial use under the
// terms of the Clear BSD License. The terms of the Clear BSD License are
// provided below. Other licenses may be obtained by contacting the faculty
// of the Columbia Computer Graphics Group or a Columbia University licensing officer.
//
// The Clear BSD License
//
// Copyright 2017 Yun (Raymond) Fei, Henrique Teles Maia, Christopher Batty,
// Changxi Zheng, and Eitan Grinspun
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted (subject to the limitations in the disclaimer
// below) provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice, this
//  list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright notice,
//  this list of conditions and the following disclaimer in the documentation
//  and/or other materials provided with the distribution.
//
// * Neither the name of the copyright holder nor the names of its contributors may be used
//  to endorse or promote products derived from this software without specific
//  prior written permission.
//
// NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE GRANTED BY THIS
// LICENSE. THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
// GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
// OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.

#include "WetHairCore.h"
#include "fluidsim2D.h"
#include "fluidsim3D.h"
#include "TimingUtilities.h"

#include <iomanip>

template<int DIM>
WetHairCore<DIM>::WetHairCore( TwoDScene<DIM>* scene, SceneStepper<DIM>* scene_stepper, std::function<bool(double)>&& script_callback )
: m_scene(scene)
, m_scene_stepper(scene_stepper)
, m_timing_statistics(17, 0)
, m_current_time(0)
, m_script_callback(std::move(script_callback))
{
}

template<int DIM>
WetHairCore<DIM>::~WetHairCore()
{
  if( m_scene != NULL )
  {
    delete m_scene;
    m_scene = NULL;
  }
  
  if( m_scene_stepper != NULL )
  {
    delete m_scene_stepper;
    m_scene_stepper = NULL;
  }
}

template<int DIM>
const scalar& WetHairCore<DIM>::getCurrentTime() const
{
  return m_current_time;
}

template<int DIM>
const std::vector<scalar>& WetHairCore<DIM>::getStepperTimingStatistics() const
{
  return m_scene_stepper->getTimingStatistics();
}
/////////////////////////////////////////////////////////////////////////////
// Simulation Control Functions

// #define USE_CFL

template<int DIM>
void WetHairCore<DIM>::stepSystem()
{
  assert( m_scene != NULL );
  assert( m_scene_stepper != NULL );
  
  const scalar dt = m_scene->getDt();
  
  const WetHairParameter& parameter = m_scene->getParameter();
  
  scalar t = 0;
  
  while(t < dt) {
#ifdef USE_CFL
    scalar substep = sim ? (sim->cfl() * 3.0) : dt;
#else
    scalar substep = dt;
#endif
    if(t + substep > dt)
      substep = dt - t;
    
    VectorXs oldpos = m_scene->getX();
    VectorXs oldvel = m_scene->getV();
    scalar t0 = timingutils::seconds();
    scalar t1;
    
    bool updateSDF = false;
    if(m_script_callback)
      updateSDF = m_script_callback(substep);
    
    m_scene->getFluidSim()->controlSources(m_current_time, dt);
    
    m_scene->applyScript(dt);
    
    std::cout << "[advect rigid bodies]" << std::endl;
    m_scene->advectRigidBodies(substep);
    
    if(updateSDF) {
      FluidSim* fluidsim = m_scene->getFluidSim();
      if(fluidsim) fluidsim->update_boundary();
    }
    
    // 0. advect the free-flow particles (Sec. 4.2).
    if(!parameter.no_fluids) {
      std::cout << "[advect fluid simulation]" << std::endl;
      m_scene->advectFluidSim(substep);
      
      t1 = timingutils::seconds();
      m_timing_statistics[0] += (t1 - t0); // fluids advection
      t0 = t1;
      
      std::cout << "[compute liquid phi]" << std::endl;
      m_scene->computeLiquidPhi();
      
      // update particle to grid with hair particles velocity and weight of bulk liquid,
      // the former is for tracing velocity of liquid on the space occupied by hairs, thus
      // we can compute drag later; the latter is for determine if current grid cell is inside
      // bulk liquid (used by cohesion force).
      std::cout << "[update particle flow to grid]" << std::endl;
      m_scene->updateParticleFlowsToGrid(true);
      
      t1 = timingutils::seconds();
      m_timing_statistics[1] += (t1 - t0); // particle to grid
      t0 = t1;
      
      std::cout << "[add gravity - fluid simulation]" << std::endl;
      m_scene->addGravityFluidSim(substep);
      
      t1 = timingutils::seconds();
      m_timing_statistics[2] += (t1 - t0); // add gravity to grid
      t0 = t1;
    }

    if(m_scene->getNumFlows() > 0) {
      
      // 1. handle the hair dynamics (Sec. 4.1), computing the drag force using velocities sampled from the grid (Sec. 4.8), as well as the adhesive/repulsive force between hairs (Sec. 4.5).
      scalar hairsubstep = substep / (scalar) parameter.hairsteps;
      
      m_scene->updateStrandParamsTimestep( hairsubstep );

      for(int i = 0; i < parameter.hairsteps; ++i) {
        std::cout << "[hair step: " << i << "]" << std::endl;
        // Step the simulated scene forward
        
        m_scene->updatePolygonalStructure(hairsubstep);
        
        t1 = timingutils::seconds();
        m_timing_statistics[3] += (t1 - t0); // search tree
        t0 = t1;
        
        m_scene->updateStrandStartStates();
        
        m_scene->interpolateMass();
        
        m_scene_stepper->stepScene( *m_scene, hairsubstep );
        
        m_scene_stepper->accept( *m_scene, hairsubstep );
        
        m_scene_stepper->PostStepScene( *m_scene, hairsubstep );
        
        t1 = timingutils::seconds();
        m_timing_statistics[4] += (t1 - t0); // hair integrator
        t0 = t1;
      }
      
      t0 = timingutils::seconds();
      
      if(!parameter.no_fluids) {
        // 2. advect the shallow water on hair, add force (Sec. 4.3 and Sec. 4.6), and compute inter-hair flow velocity (Sec. 4.4).
        std::cout << "[update interhair flow vars]" << std::endl;
        m_scene->computeInterHairDifferentialVars();
        
        t1 = timingutils::seconds();
        m_timing_statistics[5] += (t1 - t0); // compute interhair pressure etc.
        t0 = t1;
        
        std::cout << "[update flow geometric state]" << std::endl;
        m_scene->updateFlowGeometricState();
        
        t1 = timingutils::seconds();
        m_timing_statistics[6] += (t1 - t0); // update flow geometric vars.
        t0 = t1;
        
        if(!parameter.no_swe) {
          for(int i = 0; i < parameter.swesteps; ++i) {

            scalar swesubstep = substep / (scalar) parameter.swesteps;
            std::cout << "[advect hair flows]" << std::endl;
            m_scene->advectHairFlows(swesubstep);
            
            t1 = timingutils::seconds();
            m_timing_statistics[7] += (t1 - t0); // advect SWE.
            t0 = t1;
            
            std::cout << "[add force - hair flows]" << std::endl;
            m_scene->addForceHairFlows(m_scene_stepper->getAcceleration(), swesubstep);
            
            t1 = timingutils::seconds();
            m_timing_statistics[8] += (t1 - t0); // add force SWE.
            t0 = t1;
            
            std::cout << "[update hair flows height]" << std::endl;
            m_scene->updateHairFlowsHeight(m_scene_stepper->getAcceleration(), swesubstep);
            
            m_scene->updateReservoir(swesubstep);
            
            t1 = timingutils::seconds();
            m_timing_statistics[9] += (t1 - t0); // solve transfer SWE.
            t0 = t1;
          }
        }
        
        t0 = timingutils::seconds();
        
        // 3. transfer the velocity of on-hair liquid onto collocated grid (Sec. 4.7).
        std::cout << "[update hair flow to grid]" << std::endl;
        m_scene->updateHairFlowsToGrid(dt);
        
        std::cout << "[constrain hair particles]" << std::endl;
        m_scene->constrainHairParticles();
        
        t1 = timingutils::seconds();
        m_timing_statistics[10] += (t1 - t0); // update hair to grid
        t0 = t1;
      }
    }
    
    t0 = timingutils::seconds();
    
    // update particle to grid without hair particles (only bulk liquid)
    std::cout << "[update particle flow to grid - no hair particle]" << std::endl;
    m_scene->updateParticleFlowsToGrid(false);
    
    t1 = timingutils::seconds();
    m_timing_statistics[1] += (t1 - t0); // update particle to grid again
    t0 = t1;
    
    std::cout << "[add gravity - fluid simulation]" << std::endl;
    m_scene->addGravityFluidSim(substep);
    
    t1 = timingutils::seconds();
    m_timing_statistics[2] += (t1 - t0); // add gravity again
    t0 = t1;

    if(!parameter.no_fluids) {
      std::cout << "[add drag force - fluid simulation]" << std::endl;
      m_scene->addForceFluidSim(substep);
			
			// 6. solve the Viscous and/or Poisson equation and do pressure projection (Sec. 4.7).
			if(parameter.viscous_solve && parameter.viscosity > 0.0) {
				std::cout << "[viscous solve - fluid simulation]" << std::endl;
				m_scene->viscousSolveFluidSim(substep);
			}
			t1 = timingutils::seconds();
			m_timing_statistics[11] += (t1 - t0); // add drag force
			t0 = t1;
			
      // 7. combine the two velocity fields in a momentum-conservative style (Sec. 4.7).
      std::cout << "[combine velocity field]" << std::endl;
      
      m_scene->combineVelocityField();
      t1 = timingutils::seconds();
      m_timing_statistics[12] += (t1 - t0); // combine velocity
      t0 = t1;

      std::cout << "[pressure solve - fluid simulation]" << std::endl;
      m_scene->pressureSolveFluidSim(substep);

      t1 = timingutils::seconds();
      m_timing_statistics[13] += (t1 - t0); // solve poisson equation
      t0 = t1;
      
      if(m_scene->getNumFlows() > 0) {
        std::cout << "[update hair flows from grid]" << std::endl;
        // 8. transfer the velocity difference on grid back to hair (Sec. 4.7).
        m_scene->updateHairFlowsFromGrid(substep);

        t1 = timingutils::seconds();
        m_timing_statistics[14] += (t1 - t0); // update hair from grid
        t0 = t1;
      }
      
      std::cout << "[Ryoichi Ando's correction]" << std::endl;
      m_scene->getFluidSim()->correct(substep);
      
      t0 = timingutils::seconds();
      
      // 10. transfer the velocity on grid back to the particles with APIC (Sec. 4.7).
      std::cout << "[update particle flows from grid]" << std::endl;
      m_scene->updateParticleFlowsFromGrid();

      t1 = timingutils::seconds();
      m_timing_statistics[15] += (t1 - t0); // update particle from grid
      t0 = t1;
      
      // 11. absorb/release particles on hair.
      if(!(parameter.no_absorb && parameter.no_dripping) && m_scene->getNumFlows() > 0) {
        std::cout << "[absorb/release particle flows]" << std::endl;
        m_scene->absorbReleaseParticleFlows(substep, parameter.no_absorb, parameter.no_dripping);
      }
      
      t1 = timingutils::seconds();
      m_timing_statistics[16] += (t1 - t0); // absorb release liquid
      t0 = t1;
    }

    m_scene->globalVolumeAdjust();
#ifndef NDEBUG
    m_scene->checkConsistency();
#endif

    t += substep;
    
    //std::cout << "substep: " << substep << std::endl;
  }
  
  m_current_time += dt;
  
  m_scene->addVolSummary();
}

template<int DIM>
const std::vector<scalar>& WetHairCore<DIM>::getTimingStatistics() const
{
  return m_timing_statistics;
}

/////////////////////////////////////////////////////////////////////////////
// Status Functions

template<>
void WetHairCore<2>::getBoundingBox(Vectors<2>& bb_min, Vectors<2>& bb_max)
{
  const WetHairParameter& parameter = m_scene->getParameter();
  const VectorXs& x = getScene()->getX();
  // Compute the bounds on all particle positions
  scalar max_x = -std::numeric_limits<scalar>::infinity();
  scalar min_x =  std::numeric_limits<scalar>::infinity();
  scalar max_y = -std::numeric_limits<scalar>::infinity();
  scalar min_y =  std::numeric_limits<scalar>::infinity();
  for( int i = 0; i < getScene()->getNumParticles(); ++i )
  {
    if( x(2*i) > max_x )   max_x = x(2*i);
    if( x(2*i) < min_x )   min_x = x(2*i);
    if( x(2*i+1) > max_y ) max_y = x(2*i+1);
    if( x(2*i+1) < min_y ) min_y = x(2*i+1);
  }
  
  if(!parameter.no_fluids) {
    FluidSim* fluidsim = getScene()->getFluidSim();
    
    if(fluidsim) {
      Vector3s fsmin = fluidsim->getMinBBX();
      Vector3s fsmax = fluidsim->getMaxBBX();
      
      max_x = std::max(fsmax(0), max_x);
      max_y = std::max(fsmax(1), max_y);
      min_x = std::min(fsmin(0), min_x);
      min_y = std::min(fsmin(1), min_y);
    }
  }
  
  bb_min(0) = min_x;
  bb_min(1) = min_y;
  bb_max(0) = max_x;
  bb_max(1) = max_y;
}

template<>
void WetHairCore<3>::getBoundingBox(Vectors<3>& bb_min, Vectors<3>& bb_max)
{
  const WetHairParameter& parameter = m_scene->getParameter();
  const VectorXs& x = m_scene->getX();
  
  // Compute the bounds on all particle positions
  scalar max_x = -std::numeric_limits<scalar>::infinity();
  scalar min_x =  std::numeric_limits<scalar>::infinity();
  scalar max_y = -std::numeric_limits<scalar>::infinity();
  scalar min_y =  std::numeric_limits<scalar>::infinity();
  scalar max_z = -std::numeric_limits<scalar>::infinity();
  scalar min_z =  std::numeric_limits<scalar>::infinity();
  for( int i = 0; i < m_scene->getNumParticles(); ++i )
  {
    if( x( m_scene->getDof( i ) ) > max_x )   max_x = x(  m_scene->getDof( i ) );
    if( x( m_scene->getDof( i ) ) < min_x )   min_x = x( m_scene->getDof( i ) );
    if( x( m_scene->getDof( i ) + 1 ) > max_y ) max_y = x( m_scene->getDof( i ) + 1 );
    if( x( m_scene->getDof( i ) + 1 ) < min_y ) min_y = x( m_scene->getDof( i ) + 1 );
    if( x( m_scene->getDof( i ) + 2 ) > max_z ) max_z = x( m_scene->getDof( i ) + 2 );
    if( x( m_scene->getDof( i ) + 2 ) < min_z ) min_z = x( m_scene->getDof( i ) + 2 );
  }
  
  if(!parameter.no_fluids) {
    FluidSim* fluidsim = m_scene->getFluidSim();
    
    if(fluidsim) {
      Vector3s fsmin = fluidsim->getMinBBX();
      Vector3s fsmax = fluidsim->getMaxBBX();
      
      max_x = std::max(fsmax(0), max_x);
      max_y = std::max(fsmax(1), max_y);
      max_z = std::max(fsmax(2), max_z);
      min_x = std::min(fsmin(0), min_x);
      min_y = std::min(fsmin(1), min_y);
      min_z = std::min(fsmin(2), min_z);
    }
  }
  
  bb_min(0) = min_x;
  bb_min(1) = min_y;
  bb_min(2) = min_z;
  bb_max(0) = max_x;
  bb_max(1) = max_y;
  bb_max(2) = max_z;
}

template<int DIM>
TwoDScene<DIM>* WetHairCore<DIM>::getScene() const
{
  return m_scene;
}

// explicit instantiations at bottom
template class WetHairCore<2>;
template class WetHairCore<3>;
