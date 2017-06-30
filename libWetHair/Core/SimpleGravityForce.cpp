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

#include "SimpleGravityForce.h"
#include "fluidsim2D.h"
#include "fluidsim3D.h"
#include "TwoDScene.h"
#include <tbb/tbb.h>

const static char* simplegravityname = "simplegravity";

template<int DIM>
SimpleGravityForce<DIM>::SimpleGravityForce( const Vector3s& gravity, TwoDScene<DIM>* s )
: Force()
, m_gravity(gravity)
, m_scene(s)
{
  assert( (m_gravity.array()==m_gravity.array()).all() );
  assert( (m_gravity.array()!=std::numeric_limits<scalar>::infinity()).all() );
}

template<int DIM>
SimpleGravityForce<DIM>::~SimpleGravityForce()
{}

template<>
void SimpleGravityForce<2>::preCompute( const VectorXs& x, const VectorXs& v, const VectorXs& m, const scalar& dt )
{
  int np = m_scene->getNumParticles();
  
  m_buoyancy.resize(np * 2);
  
  FluidSim2D* fluid2d = (FluidSim2D*) m_scene->getFluidSim();
  
  tbb::parallel_for(0, np, 1, [&] (int i) {
    const Vector2s& pos = x.segment<2>( i * 2 ) - v.segment<2>( i * 2 ) * dt;
    scalar gamma = fluid2d->getClampedLiquidPhiValue(pos);
    m_buoyancy.segment<2>(i * 2) = fluid2d->get_pressure_gradient(pos) * gamma;
  });
}

template<>
void SimpleGravityForce<3>::preCompute( const VectorXs& x, const VectorXs& v, const VectorXs& m, const scalar& dt )
{
  int np = m_scene->getNumParticles();
  
  m_buoyancy.resize(np * 3);
  
  FluidSim3D* fluid3d = (FluidSim3D*) m_scene->getFluidSim();
  
  if(m_scene->viscositySolve()) {
    const scalar rho_L = m_scene->getLiquidDensity();
    tbb::parallel_for(0, np, 1, [&] (int i) {
      const Vector3s& pos = x.segment<3>( m_scene->getDof( i ) ) - v.segment<3>( m_scene->getDof( i ) ) * dt;
      const scalar gamma = fluid3d->getClampedLiquidPhiValue(pos);
      m_buoyancy.segment<3>(i * 3) = (fluid3d->get_pressure_gradient(pos) - fluid3d->get_visc_impulse(pos) * rho_L) * gamma;
    });
  } else {
    tbb::parallel_for(0, np, 1, [&] (int i) {
      const Vector3s& pos = x.segment<3>( m_scene->getDof( i ) ) - v.segment<3>( m_scene->getDof( i ) ) * dt;
      const scalar gamma = fluid3d->getClampedLiquidPhiValue(pos);
      m_buoyancy.segment<3>(i * 3) = fluid3d->get_pressure_gradient(pos) * gamma;
    });
  }
}

template<int DIM>
int SimpleGravityForce<DIM>::numConstraintPos()
{
  return 0;
}

template<int DIM>
int SimpleGravityForce<DIM>::numConstraintVel()
{
  return 0;
}

template<int DIM>
void SimpleGravityForce<DIM>::addEnergyToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, scalar& E )
{
  assert( x.size() == v.size() );
  assert( x.size() == m.size() );
  
  VectorXs g = m_gravity;
  
  // Assume 0 potential is at origin
  for( int i = 0; i < m_scene->getNumParticles(); ++i ){
    E -= m( m_scene->getDof(i) ) * g.segment<DIM>(0).dot( x.segment<DIM>( m_scene->getDof(i) ) );
  } 
}

template<int DIM>
void SimpleGravityForce<DIM>::addGradEToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, VectorXs& gradE )
{
  assert( x.size() == v.size() );
  assert( x.size() == m.size() );
  assert( x.size() == gradE.size() );
  
  const int np = m_scene->getNumParticles();
  
  VectorXs g = m_gravity;
  
  tbb::parallel_for(0, np, 1, [&] (int i) {
    scalar rho = m_scene->getHairDensity(i);
    gradE.segment<DIM>( m_scene->getDof(i) ) += m( m_scene->getDof(i) ) * (m_buoyancy.segment<DIM>(i * DIM) / rho - g.segment<DIM>(0));
  });
}

template<int DIM>
void SimpleGravityForce<DIM>::addHessXToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, TripletXs& hessE )
{
  assert( x.size() == v.size() );
  assert( x.size() == m.size() );
  // Nothing to do.
}

template<int DIM>
void SimpleGravityForce<DIM>::addHessVToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, TripletXs& hessE )
{
  assert( x.size() == v.size() );
  assert( x.size() == m.size() );
  // Nothing to do.
}

template<int DIM>
void SimpleGravityForce<DIM>::addGradEToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, VectorXs& gradE, int pidx )
{
  assert( x.size() == v.size() );
  assert( x.size() == m.size() );
  assert( x.size() == gradE.size() );

  int idir = m_scene->getComponent( pidx );
  if( idir == DIM ) return;

  scalar rho = m_scene->getHairDensity(pidx);

  gradE(pidx) = m(pidx) * (m_buoyancy(pidx) / rho - m_gravity(idir));
}

template<int DIM>
void SimpleGravityForce<DIM>::addHessXToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, VectorXs& hessE, int pidx )
{
  assert( x.size() == v.size() );
  assert( x.size() == m.size() );
  // Nothing to do.
}

template<int DIM>
void SimpleGravityForce<DIM>::addHessVToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, VectorXs& hessE, int pidx )
{
  assert( x.size() == v.size() );
  assert( x.size() == m.size() );
  // Nothing to do.
}

template<int DIM>
int SimpleGravityForce<DIM>::numJ()
{
  return 0;
}

template<int DIM>
int SimpleGravityForce<DIM>::numJv()
{
  return 0;
}

template<int DIM>
int SimpleGravityForce<DIM>::numJxv()
{
  return 0;
}

template<int DIM>
int SimpleGravityForce<DIM>::numTildeK()
{
  return 0;
}

template<int DIM>
bool SimpleGravityForce<DIM>::isParallelized()
{
  return true;
}

template<int DIM>
bool SimpleGravityForce<DIM>::isPrecomputationParallelized()
{
  return true;
}

template<int DIM>
const char* SimpleGravityForce<DIM>::name()
{
  return simplegravityname;
}

template<int DIM>
const char* SimpleGravityForce<DIM>::static_name()
{
  return simplegravityname;
}

template<int DIM>
Force* SimpleGravityForce<DIM>::createNewCopy()
{
  return new SimpleGravityForce(*this);
}

template<int DIM>
void SimpleGravityForce<DIM>::getAffectedVars( int pidx, std::unordered_set<int>& vars )
{
  if( m_scene->getComponent(pidx) < DIM ){
    vars.insert(pidx);
  }
}

template<int DIM>
bool SimpleGravityForce<DIM>::isContained( int pidx )
{
  if( m_scene->getComponent(pidx) == DIM ) return false;
  return true;
}

template<int DIM>
bool SimpleGravityForce<DIM>::isExternal()
{
  return true;
}
// explicit instantiations at bottom
template class SimpleGravityForce<2>;
template class SimpleGravityForce<3>;
