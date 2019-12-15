//
// This file is part of the libWetHair open source project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright 2017 Yun (Raymond) Fei, Henrique Teles Maia, Christopher Batty,
// Changxi Zheng, and Eitan Grinspun
//


#include "TwoDScene.h"
#include "SimpleGravityForce.h"
#include "MathUtilities.h"
#include "fluidsim2D.h"
#include "fluidsim3D.h"
#include "PolygonalCohesion.h"
#include "DER/StrandForce.h"
#include "ThreadUtils.h"
#include <iostream>
#include <set>
#include <stack>

WetHairParameter::WetHairParameter()
: dt(0.004),
rho(1.0),
sigma(72.0),
theta(M_PI / 3),
viscosity(8.9e-3),
friction(0.36),
max_limited_eta_prop(6.0),
latitude(0.71226395812),
earth_radius(6.37662216306e8),
earth_rotation(1.160576e-5 * 2.0 * M_PI),
height_smooth(1.0),
regularizer_shell(0.0),
capillary_accel_multiplier(0.0),
dripping_radius_multiplier(0.5),
bulk_threshold_multiplier(0.5),
absorptionRate(0.001),
airviscosity(0.0),
drag_radius_multiplier(1.0),
quadratic_dragging(1.0),
hair_hair_cohesion_multiplier(1.0),
hair_solid_cohesion_multiplier(1.0),
radius_multiplier(1.6),
collision_stiffness(10000.0),
radius_multiplier_planar(1.1),
collision_stiffness_planar(10000.0),
damping_multiplier(0.0),
damping_multiplier_planar(0.0),
friction_multiplier_planar(0.0),
hairsteps(1),
swesteps(1),
fluidcorrectionsteps(8),
drippingnear(true),
drippingfar(true),
drippingmiddle(true),
no_fluids(false),
no_swe(false),
no_absorb(false),
no_dripping(false),
no_fictitious(false),
use_ctcd(false),
global_volume_control(true),
individual_transfer(true),
volume_summary(false),
viscous_solve(false),
apply_coriolis(false),
mass_update_mode(MUM_MOMENTUM),
gravity(0.0, -981.0, 0.0)
{
}

template<int DIM>
TwoDScene<DIM>::TwoDScene( const bool& isMassSpring )
: m_x()
, m_base_x()
, m_v()
, m_m()
, m_interpolated_m()
, m_fixed()
, m_radii()
, m_edges()
, m_edge_radii()
, m_forces()
, m_particle_tags()
, m_fluid_sim(NULL)
, m_polygonal_cohesion(NULL)
, m_strout("log.txt")
, m_massSpringSim( isMassSpring )
{}

template<int DIM>
TwoDScene<DIM>::TwoDScene( const TwoDScene& otherscene, const bool& isMassSpring )
: m_x(otherscene.m_x)
, m_base_x(otherscene.m_base_x)
, m_v(otherscene.m_v)
, m_m(otherscene.m_m)
, m_interpolated_m(otherscene.m_interpolated_m)
, m_rest_m(otherscene.m_rest_m)
, m_fixed(otherscene.m_fixed)
, m_radii()
, m_edges()
, m_edge_radii()
, m_forces()
, m_particle_tags()
, m_fluid_sim(NULL)
, m_polygonal_cohesion(NULL)
, m_strout("log.txt")
, m_massSpringSim( isMassSpring )
{
  m_forces.resize(otherscene.m_forces.size());
  for( std::vector<Force*>::size_type i = 0; i < m_forces.size(); ++i ) m_forces[i] = otherscene.m_forces[i]->createNewCopy();
}

template<int DIM>
TwoDScene<DIM>::~TwoDScene()
{
  for( std::vector<Force*>::size_type i = 0; i < m_forces.size(); ++i )
  {
    assert( m_forces[i] != NULL );
    delete m_forces[i];
    m_forces[i] = NULL;
  }
  
  //annClose();
}

template<int DIM>
bool TwoDScene<DIM>::isMassSpring() const
{
  return m_massSpringSim;
}

template<int DIM>
bool TwoDScene<DIM>::isIndividualTransfer() const
{
  return m_parameters.individual_transfer;
}

template<int DIM>
bool TwoDScene<DIM>::applyCoriolis() const
{
  return m_parameters.apply_coriolis;
}


template<int DIM>
void TwoDScene<DIM>::updateHairConnectivity()
{
  int nf = m_flows.size();
  
  m_particle_to_hairs.resize(getNumParticles(), -1);
  m_particle_to_hair_local_indices.resize(getNumParticles(), -1);
  
  for(int i = 0; i < nf; ++i)
  {
    HairFlow<DIM>* flow = m_flows[i];
    auto& indices = flow->getParticleIndices();
    int nfp = indices.size();
    for(int j = 0; j < nfp; ++j)
    {
      int pidx = indices[j];
      m_particle_to_hairs[pidx] = i;
      m_particle_to_hair_local_indices[pidx] = j;
    }
  }
}

template<int DIM>
const std::vector<int>& TwoDScene<DIM>::getParticleToHairs() const
{
  return m_particle_to_hairs;
}

template<int DIM>
const std::vector<int>& TwoDScene<DIM>::getParticleToHairLocalIndices() const
{
  return m_particle_to_hair_local_indices;
}

template<int DIM>
const scalar& TwoDScene<DIM>::getLatitude() const
{
  return m_parameters.latitude;
}

template<int DIM>
const scalar& TwoDScene<DIM>::getEarthRadius() const
{
  return m_parameters.earth_radius;
}

template<int DIM>
const scalar& TwoDScene<DIM>::getEarthRotation() const
{
  return m_parameters.earth_rotation;
}

template<int DIM>
const scalar& TwoDScene<DIM>::getHeightSmooth() const
{
  return m_parameters.height_smooth;
}

template<int DIM>
const MASS_UPDATE_MODE TwoDScene<DIM>::getMassUpdateMode() const
{
  return m_parameters.mass_update_mode;
}

template<int DIM>
const scalar& TwoDScene<DIM>::getDrippingRadiusMultiplier() const
{
  return m_parameters.dripping_radius_multiplier;
}

template<int DIM>
const scalar& TwoDScene<DIM>::getDragRadiusMultiplier() const
{
  return m_parameters.drag_radius_multiplier;
}

template<int DIM>
const scalar& TwoDScene<DIM>::getBulkThresholdMultiplier() const
{
  return m_parameters.bulk_threshold_multiplier;
}

template<int DIM>
int TwoDScene<DIM>::getHairSteps() const
{
  return m_parameters.hairsteps;
}

template<int DIM>
int TwoDScene<DIM>::getSWESteps() const
{
  return m_parameters.swesteps;
}

template<int DIM>
const scalar& TwoDScene<DIM>::getDt() const
{
  return m_parameters.dt;
}

template<int DIM>
void TwoDScene<DIM>::setLiquidParameter(const WetHairParameter& parameter
                                        )
{
  m_parameters = parameter;
}

template<int DIM>
const scalar& TwoDScene<DIM>::getHairHairCohesionMultiplier() const
{
  return m_parameters.hair_hair_cohesion_multiplier;
}

template<int DIM>
const scalar& TwoDScene<DIM>::getHairSolidCohesionMultiplier() const
{
  return m_parameters.hair_solid_cohesion_multiplier;
}

template<int DIM>
const scalar& TwoDScene<DIM>::getCapillaryAccelMultiplier() const
{
  return m_parameters.capillary_accel_multiplier;
}

template<int DIM>
void TwoDScene<DIM>::setPolygonalCohesion(PolygonalCohesion<DIM>* cohesion)
{
  m_polygonal_cohesion = cohesion;
}

template<int DIM>
const PolygonalCohesion<DIM>* TwoDScene<DIM>::getPolygonalCohesion() const
{
  return m_polygonal_cohesion;
}

template<int DIM>
void TwoDScene<DIM>::setEdgeRestLength( int idx, const scalar& l0 )
{
  m_edge_rest_length(idx) = l0;
}

template<int DIM>
void TwoDScene<DIM>::setEdgePoissonRatio( int idx, const scalar& gamma )
{
  m_edge_poisson_ratio(idx) = gamma;
}

template<int DIM>
void TwoDScene<DIM>::updateFlowGeometricState()
{
  const int nflows = m_flows.size();
  threadutils::thread_pool::ParallelFor(0, nflows, [&] (int i) {
    m_flows[i]->updateGeometricState(m_x, m_v, m_fluid_sim);
  });
}

template<int DIM>
const scalar& TwoDScene<DIM>::getSearchRadius() const
{
  return m_search_radius;
}

template<int DIM>
const std::vector< std::vector<int> >& TwoDScene<DIM>::getParticleToEdge() const
{
  return m_particle_to_edge;
}

template<int DIM>
void TwoDScene<DIM>::initializeCenters()
{
  if(m_x.size() >= DIM) {
    VectorXs area(m_edges.size());
    
    mathutils::compute_edge_area<DIM>(m_x, m_edges, area, this);
    
    m_collision_keep_proximity = m_edge_radii.sum() / (scalar) (m_edge_radii.size());

    m_search_radius = m_collision_keep_proximity * m_parameters.max_limited_eta_prop * 2.0;
  }
}

template<int DIM>
void TwoDScene<DIM>::updateSearchRadius()
{
  int np = getNumParticles();
  scalar sum = 0.0;
  for(int i = 0; i < np; ++i)
  {
    sum += m_v.segment<DIM>(getDof(i)).norm();
  }
  if(np) sum /= (scalar) np;
  
  m_search_radius = std::min(m_collision_keep_proximity * m_parameters.max_limited_eta_prop * 6.0, m_collision_keep_proximity * m_parameters.max_limited_eta_prop * 2.0 + sum);
}

template<int DIM>
const FluidSim* TwoDScene<DIM>::getFluidSim() const
{
  return m_fluid_sim;
}

template<int DIM>
FluidSim* TwoDScene<DIM>::getFluidSim()
{
  return m_fluid_sim;
}

template<int DIM>
void TwoDScene<DIM>::setFluidSim(FluidSim* sim)
{
  m_fluid_sim = sim;
}

template<int DIM>
int TwoDScene<DIM>::getNumFlows() const
{
  return m_flows.size();
}

template<int DIM>
const scalar& TwoDScene<DIM>::getLiquidTheta() const
{
  return m_parameters.theta;
}

template<int DIM>
scalar TwoDScene<DIM>::getDStarPlanar(const scalar& radius) const
{
  return getPolygonalCohesion()->getDStarPlanar(radius);
}

template<int DIM>
scalar TwoDScene<DIM>::getStiffnessPlanar(const scalar& radius, const scalar& d0, const scalar& A_target, const scalar& pressure_weight) const
{
  return getPolygonalCohesion()->getStiffnessPlanar(radius, d0, A_target, pressure_weight);
}

template<int DIM>
scalar TwoDScene<DIM>::getDStar(const scalar& radius) const
{
  return getPolygonalCohesion()->getDStar(radius);
}

template<int DIM>
scalar TwoDScene<DIM>::getStiffness(const scalar& radius, const scalar& d0, const scalar& A_target, const scalar& pressure_weight) const
{
  return getPolygonalCohesion()->getStiffness(radius, d0, A_target, pressure_weight);
}

template<int DIM>
void TwoDScene<DIM>::updateBoundingBox()
{
  m_bb_min.setConstant(1e+20);
  m_bb_max.setConstant(-1e+20);
  
  int np = getNumParticles();
  for(int i = 0; i < np; ++i)
  {
    m_bb_min = m_bb_min.cwiseMin( m_x.segment<DIM>( getDof( i ) ) );
    m_bb_max = m_bb_max.cwiseMax( m_x.segment<DIM>( getDof( i ) ) );
  }
}

template<int DIM>
const Vectors<DIM>& TwoDScene<DIM>::getBoundingBoxMin() const
{
  return m_bb_min;
}

template<int DIM>
const Vectors<DIM>& TwoDScene<DIM>::getBoundingBoxMax() const
{
  return m_bb_max;
}

template<int DIM>
void TwoDScene<DIM>::initializeLiquids()
{
  initializeCenters();
  updateBoundingBox();
  if(m_polygonal_cohesion) {
    m_polygonal_cohesion->updateStructure(m_x);
    m_polygonal_cohesion->preCompute(m_x, m_v, m_m, 0.0);
    m_polygonal_cohesion->postStepScene(0.0);
  }
  computeInterHairDifferentialVars();
  updateFlowGeometricState();
  
  if(m_parameters.mass_update_mode != MUM_NONE || m_parameters.regularizer_shell > 0.0)
    updateHairMass();
  
  m_fluid_sim->init_hair_particles();
  m_fluid_sim->sort_particles();
  m_fluid_sim->compute_liquid_phi();
  
  updateVolumeRecord();
  m_volume_hair_old = m_volume_hair;
  m_volume_particle_old = m_volume_particle;
  m_volume_reservoir_old = m_volume_reservoir;
  
  m_volume_particle_removed = m_volume_particle_inserted = 0.0;
}

template<int DIM>
void TwoDScene<DIM>::updateHairMass()
{
  for (HairFlow<DIM>* flow : m_flows) {
    flow->updateHairMass();
  }
}

template<int DIM>
void TwoDScene<DIM>::advectHairFlows(const scalar& dt)
{
  const int nflows = m_flows.size();
  threadutils::thread_pool::ParallelFor(0, nflows, [&] (int i) {
    m_flows[i]->advance(m_x, dt);
  });
}

template<int DIM>
void TwoDScene<DIM>::computeInterHairDifferentialVars()
{
  if(!m_parameters.individual_transfer && m_polygonal_cohesion)
    m_polygonal_cohesion->computeInterHairVariables();
}

template<int DIM>
void TwoDScene<DIM>::addForceHairFlows(const VectorXs& accel, const scalar& dt)
{
  if(!m_parameters.individual_transfer && m_polygonal_cohesion)
    m_polygonal_cohesion->computeGlobalHairPressure();
  
  const int nflows = m_flows.size();
  threadutils::thread_pool::ParallelFor(0, nflows, [&] (int i) {
    m_flows[i]->add_force(m_x, accel, m_fluid_sim, dt);
  });
  
  if(!m_parameters.individual_transfer && m_polygonal_cohesion)
    m_polygonal_cohesion->computeInterHairVelocity(dt);
}

template<int DIM>
void TwoDScene<DIM>::updateHairFlowsToGrid(const scalar& dt)
{
  if(!m_fluid_sim) return;
  
  m_fluid_sim->prepare_update_from_hair();
  
  // update drag buffer
  m_fluid_drag_buffer.setZero();
  m_fluid_sim->preCompute(m_x, m_v, m_m, dt);
  m_fluid_sim->addGradEToTotal(m_x, m_v, m_m, m_fluid_drag_buffer);
  int ncol = m_m.size();
  for(int i = 0; i < ncol; ++i) {
    m_fluid_drag_buffer(i) /= m_m(i);
  }
  
  int nflows = m_flows.size();
  threadutils::thread_pool::ParallelFor(0, nflows, [&] (int iflow) {
    m_flows[iflow]->updateToFilteredGrid(m_x, m_v, m_fluid_sim, dt, iflow);
  });
  
  m_fluid_sim->done_update_from_hair();
}

template<int DIM>
const VectorXs& TwoDScene<DIM>::getFluidDragBuffer() const
{
  return m_fluid_drag_buffer;
}

template<int DIM>
void TwoDScene<DIM>::updateHairFlowsFromGrid(const scalar& dt)
{
  if(!m_fluid_sim) return;
  
  const int nflows = m_flows.size();
  threadutils::thread_pool::ParallelFor(0, nflows, [&] (int i) {
    m_flows[i]->updateFromFilteredGrid(m_x, m_v, m_fluid_sim, dt);
  });
}

template<int DIM>
void TwoDScene<DIM>::updatePolygonalStructure(const scalar& dt)
{
  m_polygonal_cohesion->updateStructure(m_x + m_v * dt);
}

template<int DIM>
void TwoDScene<DIM>::updateHairFlowsHeight(const VectorXs& accel, const scalar& dt)
{
  const int nflows = m_flows.size();
  threadutils::thread_pool::ParallelFor(0, nflows, [&] (int i) {
    m_flows[i]->preUpdateHairFlowHeight(dt);
  });
  
  if(!m_parameters.individual_transfer && m_polygonal_cohesion)
    m_polygonal_cohesion->computeInterHairRHS(dt);
  
  threadutils::thread_pool::ParallelFor(0, nflows, [&] (int i) {
    m_flows[i]->updateHairFlowHeight(dt);
    m_flows[i]->postUpdateHairFlowHeight(dt);
  });
}

template<int DIM>
void TwoDScene<DIM>::advectRigidBodies(const scalar& dt)
{
  if(!m_fluid_sim) return;
  
  m_fluid_sim->advect_boundary(dt);
}

template<int DIM>
const scalar& TwoDScene<DIM>::getAbsorptionRate() const
{
  return m_parameters.absorptionRate;
}

template<int DIM>
const scalar& TwoDScene<DIM>::getAirViscosity() const
{
  return m_parameters.airviscosity;
}

template<int DIM>
void TwoDScene<DIM>::uniformlyIncreaseHairLiquid(const scalar& dh, int group_idx)
{
  int np = getNumParticles();
  for(int i = 0; i < np; ++i)
  {
    if(m_script_group[i] == group_idx) {
      int hidx = getParticleToHairs()[i];
      int local_idx = getParticleToHairLocalIndices()[i];
      
      m_flows[hidx]->getEta()(local_idx) = std::max(0.0, m_flows[hidx]->getEta()(local_idx) + dh);
    }
  }
}

template<int DIM>
void TwoDScene<DIM>::constrainHairParticles()
{
  if(!m_fluid_sim) return;
  
  m_fluid_sim->constrain_hair_particles();
}

template<int DIM>
void TwoDScene<DIM>::advectFluidSim(const scalar& dt)
{
  if(!m_fluid_sim) return;
  
  std::cout << "[ADV: passively advect particles]" << std::endl;
  m_fluid_sim->advect_particles(dt);
}

template<int DIM>
void TwoDScene<DIM>::addForceFluidSim(const scalar& dt)
{
  if(!m_fluid_sim) return;
  
  m_fluid_sim->add_drag(dt);
}

template<int DIM>
void TwoDScene<DIM>::addGravityFluidSim(const scalar& dt)
{
  m_fluid_sim->add_gravity(dt);
}

template<int DIM>
void TwoDScene<DIM>::viscousSolveFluidSim(const scalar& dt)
{
  if(!m_fluid_sim) return;
  
  m_fluid_sim->apply_viscosity(dt);
}

template<int DIM>
void TwoDScene<DIM>::pressureSolveFluidSim(const scalar& dt)
{
  if(!m_fluid_sim) return;
  
  m_fluid_sim->project(dt);
  
  m_fluid_sim->constrain_velocity();
}

template<int DIM>
void TwoDScene<DIM>::globalVolumeAdjust()
{
  if(!m_parameters.global_volume_control) return;
  
  updateVolumeRecord();
  
  scalar new_vol_hair_lim = m_volume_hair_old + m_volume_particle_old - m_volume_particle + m_volume_particle_inserted - m_volume_particle_removed + m_volume_reservoir_old - m_volume_reservoir;
  
  std::cout << "[vol statistics: old hair: " << m_volume_hair_old << " new hair: " << m_volume_hair << " old particle: " << m_volume_particle_old << " new particle: " << m_volume_particle << " inserted particle: " << m_volume_particle_inserted << " removed particle: " << m_volume_particle_removed << " old reservoir: " << m_volume_reservoir_old << " new reservoir: " << m_volume_reservoir << "]" << std::endl;
  
  if(m_volume_hair > 0.0 && m_volume_hair > new_vol_hair_lim) {
    scalar prop = new_vol_hair_lim / m_volume_hair;
    std::cout << "[volume adj: vol: " << m_volume_hair << " limit: " << new_vol_hair_lim << " prop: " << prop << "]" << std::endl;
    
    for(auto& flow : m_flows)
    {
      flow->adjustVolumeGlobal(prop);
    }
    
    if(m_parameters.mass_update_mode != MUM_NONE)
    {
      updateHairMass();
    }
    
    updateVolumeRecord();
  }
}

template<int DIM>
void TwoDScene<DIM>::reportParticleRemoved(const scalar& vol_removed)
{
  m_volume_particle_removed = vol_removed;
}

template<int DIM>
void TwoDScene<DIM>::reportParticleAdded(const scalar& vol_added)
{
  m_volume_particle_inserted = vol_added;
}

template<int DIM>
void TwoDScene<DIM>::updateVolumeRecord()
{
  m_volume_hair_old = m_volume_hair;
  m_volume_particle_old = m_volume_particle;
  m_volume_reservoir_old = m_volume_reservoir;
  
  scalar sum_liq = 0.0;
  scalar sum_res = 0.0;
  for(auto& flow : m_flows)
  {
    sum_liq += flow->computeTotalLiquidVol();
    sum_res += flow->computeTotalReservoirVol();
  }
  
  scalar sum_free = m_fluid_sim->computeTotalLiquidVol();
  
  m_volume_hair = sum_liq;
  m_volume_reservoir = sum_res;
  m_volume_particle = sum_free;
}

template<int DIM>
void TwoDScene<DIM>::updateReservoir(const scalar& dt)
{
  for (HairFlow<DIM>* flow : m_flows) {
    flow->updateReservoir(m_fluid_sim, m_x, m_v, dt);
  }
}

template<int DIM>
void TwoDScene<DIM>::absorbReleaseParticleFlows(const scalar& dt, bool noabsorb, bool nodripping)
{
  if(m_fluid_sim) {
    if(!noabsorb) m_fluid_sim->shareParticleWithHairs(m_x, dt);
    if(!nodripping) m_fluid_sim->transferLiquidToGridParticle(dt);
  }
}

template<int DIM>
void TwoDScene<DIM>::computeLiquidPhi()
{
  if(!m_fluid_sim) return;
  
  m_fluid_sim->compute_liquid_phi();
}

template<int DIM>
void TwoDScene<DIM>::combineVelocityField()
{
  if(!m_fluid_sim) return;
  
  m_fluid_sim->combine_velocity_field();
}

template<int DIM>
void TwoDScene<DIM>::updateParticleFlowsToGrid(bool with_hair_particles)
{
  if(!m_fluid_sim) return;
  
  m_fluid_sim->map_p2g(with_hair_particles);
}

template<int DIM>
void TwoDScene<DIM>::updateParticleFlowsFromGrid()
{
  if(!m_fluid_sim) return;
  
  m_fluid_sim->map_g2p_apic();
}

template<int DIM>
const std::vector<HairFlow<DIM>*>& TwoDScene<DIM>::getFilmFlows() const
{
  return m_flows;
}

template<int DIM>
std::vector<HairFlow<DIM>*>& TwoDScene<DIM>::getFilmFlows()
{
  return m_flows;
}

template<int DIM>
void TwoDScene<DIM>::insertFilmFlow(HairFlow<DIM>* flow)
{
  m_flows.push_back(flow);
}

template<int DIM>
int TwoDScene<DIM>::getNumParticles() const
{
  if( m_massSpringSim ) return m_x.size()/DIM;
  else return ( m_x.size() + m_num_strands ) / 4;
}

template<int DIM>
int TwoDScene<DIM>::getNumEdges() const
{
  return m_edges.size();
}

template<int DIM>
int TwoDScene<DIM>::getNumDofs() const
{
  return m_x.size();
}

template<int DIM>
const VectorXs& TwoDScene<DIM>::getX() const
{
  return m_x;
}

template<int DIM>
VectorXs& TwoDScene<DIM>::getX()
{
  return m_x;
}

template<int DIM>
const VectorXs& TwoDScene<DIM>::getV() const
{
  return m_v;
}

template<int DIM>
VectorXs& TwoDScene<DIM>::getV()
{
  return m_v;
}

template<int DIM>
const VectorXs& TwoDScene<DIM>::getM() const
{
  return m_m;
}

template<int DIM>
VectorXs& TwoDScene<DIM>::getM()
{
  return m_m;
}

template<>
void TwoDScene<2>::interpolateMass()
{
  int np = getNumParticles();
  FluidSim2D* fluid2d = (FluidSim2D*) m_fluid_sim;
  
  threadutils::thread_pool::ParallelFor(0, np, [&] (int i) {
    const Vector2s& pos = m_x.segment<2>(i * 2);
    scalar clamped_particle_weight = fluid2d->getClampedLiquidPhiValue(pos);

    m_interpolated_m.segment<2>(i * 2) = m_rest_m.segment<2>(i * 2) * clamped_particle_weight + m_m.segment<2>(i * 2) * (1.0 - clamped_particle_weight);
  });
}

template<>
void TwoDScene<3>::interpolateMass()
{
  int np = getNumParticles();
  FluidSim3D* fluid3d = (FluidSim3D*) m_fluid_sim;
  
  threadutils::thread_pool::ParallelFor(0, np, [&] (int i) {
    const Vector3s& pos = m_x.segment<3>(getDof(i));
    scalar clamped_particle_weight = fluid3d->getClampedLiquidPhiValue(pos);
    
    if( m_massSpringSim || isTip( i ) ){
      m_interpolated_m.segment<3>(getDof(i)) = m_rest_m.segment<3>( getDof(i) ) * clamped_particle_weight + m_m.segment<3>( getDof(i) ) * (1.0 - clamped_particle_weight);
    } else {
      m_interpolated_m.segment<4>(getDof(i)) = m_rest_m.segment<4>( getDof(i) ) * clamped_particle_weight + m_m.segment<4>( getDof(i) ) * (1.0 - clamped_particle_weight);
    }
  });
}

template<int DIM>
const VectorXs& TwoDScene<DIM>::getInterpolatedM() const
{
  return m_interpolated_m;
}

template<int DIM>
VectorXs& TwoDScene<DIM>::getInterpolatedM()
{
  return m_interpolated_m;
}

template<int DIM>
const scalar& TwoDScene<DIM>::getMaxLimitEtaProp() const
{
  return m_parameters.max_limited_eta_prop;
}

template<int DIM>
const scalar& TwoDScene<DIM>::getQuadraticDragging() const
{
  return m_parameters.quadratic_dragging;
}

template<int DIM>
const scalar& TwoDScene<DIM>::getViscosity() const
{
  return m_parameters.viscosity;
}

template<int DIM>
const scalar& TwoDScene<DIM>::getHairFrictionCoeff() const
{
  return m_parameters.friction;
}

template<int DIM>
const VectorXs& TwoDScene<DIM>::getRadii() const
{
  return m_radii;
}

template<int DIM>
const Vector3s& TwoDScene<DIM>::getSimpleGravity() const
{
  return m_parameters.gravity;
}

template<int DIM>
const scalar& TwoDScene<DIM>::getLiquidDensity() const
{
  return m_parameters.rho;
}

template<int DIM>
scalar TwoDScene<DIM>::getHairDensity(int pidx) const
{
  if(m_strandParameters.size() > 0) {
    const StrandForce* sf = m_strands[ m_particle_to_hairs[pidx] ];
    return sf->m_strandParams->m_density;
  } else {
    return 1.32;
  }
}

template<int DIM>
const scalar& TwoDScene<DIM>::getLiquidTension() const
{
  return m_parameters.sigma;
}

template<int DIM>
const VectorXs& TwoDScene<DIM>::getHairRestMass() const
{
  return m_rest_m;
}

template<int DIM>
VectorXs& TwoDScene<DIM>::getHairRestMass()
{
  return m_rest_m;
}

template<int DIM>
const std::vector<unsigned char>& TwoDScene<DIM>::getFixed() const
{
  return m_fixed;
}

template<int DIM>
std::vector<unsigned char>& TwoDScene<DIM>::getFixed()
{
  return m_fixed;
}

template<int DIM>
void TwoDScene<DIM>::setScriptedGroup( int particle, int group_idx )
{
  m_script_group[particle] = group_idx;
}

template<int DIM>
int TwoDScene<DIM>::getScriptedGroup( int particle )
{
  return m_script_group[particle];
}

template<int DIM>
Vector3s& TwoDScene<DIM>::getScriptedTranslate(int sg_idx)
{
  return m_scripted_translate[ sg_idx ];
}

template<int DIM>
const Vector3s& TwoDScene<DIM>::getScriptedTranslate(int sg_idx) const
{
  return m_scripted_translate[ sg_idx ];
}

template<int DIM>
Eigen::Quaternion<scalar>& TwoDScene<DIM>::getScriptedRotation(int sg_idx)
{
  return m_scripted_rotation[ sg_idx ];
}

template<int DIM>
const Eigen::Quaternion<scalar>& TwoDScene<DIM>::getScriptedRotation(int sg_idx) const
{
  return m_scripted_rotation[ sg_idx ];
}


template<int DIM>
void TwoDScene<DIM>::initializeScriptedGroup( const std::vector<Script>& scripts )
{
  int maximal_group_idx = 0;
  for(int gidx : m_script_group)
  {
    maximal_group_idx = std::max(gidx, maximal_group_idx);
  }
  
  const int num_groups = maximal_group_idx + 1;
  m_scripted_rotation.resize( num_groups, Eigen::Quaternion<scalar>::Identity() );
  m_scripted_translate.resize( num_groups, Vector3s::Zero() );
  
  m_base_x = m_x;
  
  const int np = getNumParticles();
  
  for(const Script& s : scripts)
  {
    if(s.target == Script::ROOT && s.start < 0.0) {
      switch(s.type) {
        case Script::ROTATE:
        {
          m_scripted_rotation[std::max(0, s.index)] = Eigen::Quaternion<scalar>(Eigen::AngleAxis<scalar>(s.v(3), Vector3s(s.v(0), s.v(1), s.v(2)) ) );
          break;
        }
        case Script::TRANSLATE:
        {
          m_scripted_translate[std::max(0, s.index)] = Vector3s(s.v(0), s.v(1), s.v(2));
          break;
        }
        case Script::SCALE:
        {
          std::cerr << "Script::SCALE: NOT IMPLEMENTED!" << std::endl;
          break;
        }
        default:
          break;
      }
    } else if(s.target == Script::ALL && s.start < 0.0) {
      switch(s.type) {
        case Script::ROTATE:
        {
          Eigen::Quaternion<scalar> q = Eigen::Quaternion<scalar>(Eigen::AngleAxis<scalar>(s.v(3), Vector3s(s.v(0), s.v(1), s.v(2)) ) );
          m_scripted_rotation[std::max(0, s.index)] = q;
          for(int i = 0; i < np; ++i)
          {
            if(!isFixed(i) && (s.index < 0 || m_script_group[i] == s.index)) {
              const Vectors<DIM>& x0 = m_x.segment<DIM>( getDof(i) );
              Eigen::Quaternion<scalar> p(0.0, x0(0), x0(1), x0(2));
              VectorXs x1 = (q * p * q.inverse()).vec();
              m_x.segment<DIM>( getDof(i) ) = x1.segment<DIM>(0);
            }
          }
          break;
        }
        case Script::TRANSLATE:
        {
          m_scripted_translate[std::max(0, s.index)] = Vector3s(s.v(0), s.v(1), s.v(2));
          for(int i = 0; i < np; ++i)
          {
            if(!isFixed(i) && (s.index < 0 || m_script_group[i] == s.index)) {
              VectorXs xx = Vector3s(s.v(0), s.v(1), s.v(2));
              m_x.segment<DIM>( getDof(i) ) += xx.segment<DIM>(0);
            }
          }
          break;
        }
        case Script::SCALE:
        {
          std::cerr << "Script::SCALE: NOT IMPLEMENTED!" << std::endl;
          break;
        }
        default:
          break;
      }
    }
  }

  
  applyScript(0.0);
}

template<int DIM>
void TwoDScene<DIM>::applyScript(const scalar& dt)
{
  const int np = getNumParticles();
  for(int i = 0; i < np; ++i)
  {
    if(!isFixed(i)) continue;
    
    int sg_idx = m_script_group[i];
    if(sg_idx < 0 || sg_idx >= m_scripted_translate.size()) continue;
    
    const Eigen::Quaternion<scalar>& q = m_scripted_rotation[ sg_idx ];
    const Vector3s& t = m_scripted_translate[ sg_idx ];
    
    Eigen::Quaternion<scalar> p;
    const VectorXs& x0 = m_base_x.segment<DIM>( getDof( i ) );
    const Vectors<DIM>& xstar = m_x.segment<DIM>( getDof( i ) );
    Eigen::Quaternion<scalar> p0(0.0, x0(0), x0(1), DIM > 2 ? x0(2) : 0.0);
    VectorXs trans_x0 = (q * p0 * q.inverse()).vec() + t;
    
    if(dt == 0.0) {
      m_x.segment<DIM>( getDof( i ) ) = trans_x0.segment<DIM>(0);
    } else {
      m_v.segment<DIM>( getDof( i ) ) = (trans_x0.segment<DIM>(0) - xstar) / dt;
    }
  }
  
  if(DIM == 3 && !m_massSpringSim) {
    int nstrand = m_strandEquilibriumParameters.size();
    for(int i = 0; i < nstrand; ++i)
    {
      if(!m_strandEquilibriumParameters[i] || !m_strandEquilibriumParameters[i]->m_valid || !m_strandEquilibriumParameters[i]->m_dirty) continue;
      
      updateCurlyHair(m_strandEquilibriumParameters[i]->m_dL,
                      m_strandEquilibriumParameters[i]->m_vertices,
                      m_strandEquilibriumParameters[i]->m_curl_radius,
                      m_strandEquilibriumParameters[i]->m_curl_density,
                      m_strandEquilibriumParameters[i]->m_root_length);
      
      int nverts = (int) m_strandEquilibriumParameters[i]->m_vertices.size();
      VecX dof_restshape(nverts * 4 - 1);
      dof_restshape.setZero();
      
      for(int j = 0; j < nverts; ++j)
      {
        dof_restshape.segment<3>(j * 4) = m_strandEquilibriumParameters[i]->m_vertices[j];
      }
      m_strands[i]->updateRestShape( dof_restshape );
      
      m_strandEquilibriumParameters[i]->m_dirty = false;
    }
  }
}

template<int DIM>
void TwoDScene<DIM>::resizeSystem( int num_particles, int num_edges, int num_strands )
{
  if( !m_massSpringSim ){
    resizeSystemDER( num_particles, num_edges, num_strands );
    return;
  }
  assert( num_particles >= 0 );

  m_script_group.resize( num_particles, 0 );
  m_fluid_drag_buffer.resize(DIM*num_particles);
  m_x.resize(DIM*num_particles);
  m_base_x.resize(DIM*num_particles);
  m_v.resize(DIM*num_particles);
  m_m.resize(DIM*num_particles);
  m_interpolated_m.resize(DIM*num_particles);
  m_rest_m.resize(DIM*num_particles);
  m_fixed.resize(num_particles);
  m_radii.resize(num_particles);
  m_particle_tags.resize(num_particles);
  m_edges.resize(num_edges);
  m_edge_radii.resize(num_edges);
  m_edge_rest_radii.resize(num_edges);
  m_edge_rest_length.resize(num_edges);
  m_edge_poisson_ratio.resize(num_edges);
  m_particle_to_edge.resize(num_particles);
  m_num_strands = num_strands;
}

template<int DIM>
void TwoDScene<DIM>::resizeSystemDER( int num_particles, int num_edges, int num_strands )
{
  assert( num_particles >= 0 );

  int numDofs = (4 * num_particles) - num_strands; // [ x y z t ] for all strand vertices except tip,
  m_fluid_drag_buffer.resize( numDofs ); m_fluid_drag_buffer.setZero();
  m_x.resize( numDofs ); m_x.setZero();
  m_base_x.resize( numDofs ); m_base_x.setZero();
  m_v.resize( numDofs ); m_v.setZero();
  m_m.resize( numDofs ); m_m.setZero();
  m_interpolated_m.resize( numDofs ); m_interpolated_m.setZero();
  m_rest_m.resize( numDofs ); m_rest_m.setZero();

  m_script_group.resize( num_particles, 0 );
  m_particle_to_dofs.resize( num_particles );
  m_dofs_to_component.resize( numDofs );
  m_is_strand_tip.resize( num_particles );

  m_fixed.resize( num_particles );
  m_radii.resize( num_particles );
  m_particle_tags.resize( num_particles );
  m_edges.resize( num_edges );
  m_edge_radii.resize(num_edges);
  m_edge_rest_radii.resize(num_edges);
  m_edge_rest_length.resize(num_edges);
  m_edge_poisson_ratio.resize(num_edges);
  m_particle_to_edge.resize(num_particles);
  m_num_strands = num_strands;
  m_edge_to_hair.resize( num_edges );
  
}

template<int DIM>
void TwoDScene<DIM>::setPosition( int particle, const Vectors<DIM>& pos )
{
  assert( particle >= 0 );
  if( DIM == 2 ) assert( particle < getNumParticles() );

  m_x.segment<DIM>( getDof(particle) ) = pos;
}

template<int DIM>
Vectors<DIM> TwoDScene<DIM>::getPosition( int particle )
{
  assert( particle >= 0 );
  // assert( particle < getNumParticles() );
  assert( getDof(particle) < m_x.size() );

  return m_x.segment<DIM>( getDof(particle) );
}

template< int DIM >
int TwoDScene<DIM>::getDof( int particle ) const
{
  assert( particle >= 0 );

  if( m_massSpringSim ){
    assert( particle < getNumParticles() );
    return DIM * particle;
  }
  else{
    assert( particle < m_particle_to_dofs.size() );
    return m_particle_to_dofs[particle];
  }
}

template<int DIM>
int TwoDScene<DIM>::getComponent( int dof ) const
{
  assert( dof >= 0 );
  assert( dof < m_x.size() );

  if( m_massSpringSim ) return dof % DIM;
  else return m_dofs_to_component[dof];
}

template<int DIM>
int TwoDScene<DIM>::getVertFromDof( int dof ) const
{
  assert( dof >= 0 );
  assert( dof < m_x.size() );
  if( m_massSpringSim ) return dof / DIM;
  else return m_dofs_to_particle[dof];
}

template<int DIM>
bool TwoDScene<DIM>::isTip( int particle ) const
{
  assert( particle >= 0 );
  assert( particle < getNumParticles() );

  if( !m_massSpringSim ){
    return m_is_strand_tip[particle];
  }
  return false;
}

template<int DIM>
void TwoDScene<DIM>::setVelocity( int particle, const Vectors<DIM>& vel )
{
  assert( particle >= 0 );
  assert( getDof(particle) < m_v.size() );
  
  m_v.segment<DIM>( getDof(particle) ) = vel;
}

template<int DIM>
Vectors<DIM> TwoDScene<DIM>::getVelocity( int particle )
{
  assert( particle >= 0 );
  assert( particle < getNumParticles() );
  return m_v.segment<DIM>( getDof(particle) );
}

template<int DIM>
void TwoDScene<DIM>::setMass( int particle, const scalar& mass )
{
  assert( particle >= 0 );
  assert( particle < getNumParticles() );

  if( m_massSpringSim ){
    m_m.segment<DIM>(particle * DIM).setConstant(mass);
    m_rest_m.segment<DIM>(particle * DIM).setConstant(mass);
    m_interpolated_m.segment<DIM>(particle * DIM).setConstant(mass);
  }
  else{
     std::cerr << "mass computed via strand parameters in TwoDScene<3>::computeMassesAndRadiiFromStrands() " << std::endl;
     exit(1);
  }
}

template<int DIM>
void TwoDScene<DIM>::setFixed( int particle, bool fixed )
{
  assert( particle >= 0 );
  assert( particle < m_fixed.size() );

  m_fixed[particle] = fixed;
}

template<int DIM>
bool TwoDScene<DIM>::isFixed( int particle ) const
{
  assert( particle >= 0 );
  assert( particle < getNumParticles() );

  return m_fixed[particle];
}

template<int DIM>
const scalar& TwoDScene<DIM>::getRadius( int particle ) const
{
  assert( particle >= 0 );
  assert( particle < getNumParticles() );
  
  return m_radii[particle];  
}

template<int DIM>
void TwoDScene<DIM>::setRadius( int particle, scalar radius )
{
  assert( particle >= 0 );
  assert( particle < getNumParticles() );

  m_radii[particle] = radius;
}

template<int DIM>
void TwoDScene<DIM>::clearEdges()
{
  m_edges.clear();
}

template<int DIM>
void TwoDScene<DIM>::setEdge( int idx, const std::pair<int,int>& edge, scalar radius )
{
  m_edges[idx] = edge;
  m_particle_to_edge[edge.first].push_back(idx);
  m_particle_to_edge[edge.second].push_back(idx);

  if( m_massSpringSim ){
    m_edge_radii[idx] = radius;
    m_edge_rest_radii[idx] = radius;
  }
}

template<int DIM>
const std::vector<std::pair<int,int> >& TwoDScene<DIM>::getEdges() const
{
  return m_edges;
}

template<int DIM>
const VectorXs& TwoDScene<DIM>::getEdgeRadii() const
{
  return m_edge_radii;
}

template<int DIM>
const std::pair<int,int>& TwoDScene<DIM>::getEdge(int edg) const
{
  assert( edg >= 0 );
  assert( edg < (int) m_edges.size() );

  return m_edges[edg];
}

template<int DIM>
void TwoDScene<DIM>::insertForce( Force* newforce )
{
  m_forces.push_back(newforce);
}

template<int DIM>
void TwoDScene<DIM>::insertStrandParameters( StrandParameters* newparams )
{
  m_strandParameters.push_back(newparams);
}

template<int DIM>
void TwoDScene<DIM>::insertStrandEquilibriumParameters( StrandEquilibriumParameters* newparams )
{
  m_strandEquilibriumParameters.push_back(newparams);
}

template<int DIM>
scalar TwoDScene<DIM>::computeKineticEnergy() const
{
  return computeHairParticleKineticEnergy() + computeLiquidParticleKineticEnergy();
}

template<int DIM>
scalar TwoDScene<DIM>::computePotentialEnergy() const
{
  scalar U = 0.0;
  for( std::vector<Force*>::size_type i = 0; i < m_forces.size(); ++i ) m_forces[i]->addEnergyToTotal( m_x, m_v, m_m, U );
  return U;  
}

template<int DIM>
scalar TwoDScene<DIM>::computeTotalEnergy() const
{
  return computeKineticEnergy()+computePotentialEnergy();
}

template<int DIM>
std::vector< std::unordered_set<int> >& TwoDScene<DIM>::getBroadPhasePEPairs()
{
  return m_bp_particle_edge_pairs;
}

template<int DIM>
std::vector< std::unordered_set<int> >& TwoDScene<DIM>::getBroadPhasePPPairs()
{
  return m_bp_particle_particle_pairs;
}

template<int DIM>
std::vector< std::unordered_set<int> >& TwoDScene<DIM>::getBroadPhaseEEPairs()
{
  return m_bp_edge_edge_pairs;
}

template<int DIM>
const std::vector< std::unordered_set<int> >& TwoDScene<DIM>::getBroadPhasePEPairs() const
{
  return m_bp_particle_edge_pairs;
}

template<int DIM>
const std::vector< std::unordered_set<int> >& TwoDScene<DIM>::getBroadPhasePPPairs() const
{
  return m_bp_particle_particle_pairs;
}

template<int DIM>
const std::vector< std::unordered_set<int> >& TwoDScene<DIM>::getBroadPhaseEEPairs() const
{
  return m_bp_edge_edge_pairs;
}

template<int DIM>
StrandParameters* TwoDScene<DIM>::getStrandParameters( const int index )
{
  assert( 0 <= index );
  assert( index < m_strandParameters.size() );
  return m_strandParameters[index];
}

template<int DIM>
StrandEquilibriumParameters* TwoDScene<DIM>::getStrandEquilibriumParameters( const int index )
{
  assert( 0 <= index );
  assert( index < m_strandEquilibriumParameters.size() );
  return m_strandEquilibriumParameters[index];
}

template<int DIM>
void TwoDScene<DIM>::updateStrandParamsTimestep( const scalar& dt )
{
  for( unsigned p = 0; p < m_strandParameters.size(); ++p ){
    m_strandParameters[p]->computeViscousForceCoefficients( dt );
  }
}

template<int DIM>
void TwoDScene<DIM>::updateStrandStartStates()
{
  for( unsigned s = 0; s < m_num_strands; ++s ){
    m_strands[s]->updateStartDoFs( m_x );
  }
  if(m_polygonal_cohesion) {
    m_polygonal_cohesion->updateViscousStartPhi(m_x);
  }  
}

template<int DIM>
void TwoDScene<DIM>::accumulateExternalGradU( VectorXs& F, const VectorXs& dx, const VectorXs& dv )
{
  assert( F.size() == m_x.size() );
  assert( dx.size() == dv.size() );
  assert( dx.size() == 0 || dx.size() == F.size() );
  
  // Accumulate all energy gradients
  if( dx.size() == 0 ) for( std::vector<Force*>::size_type i = 0; i < m_external_forces.size(); ++i ) {
    m_external_forces[i]->addGradEToTotal( m_x, m_v, m_interpolated_m, F );
  }
  else for( std::vector<Force*>::size_type i = 0; i < m_external_forces.size(); ++i ) {
    m_external_forces[i]->addGradEToTotal( m_x+dx, m_v+dv, m_interpolated_m, F );
  }
}

template<int DIM>
void TwoDScene<DIM>::accumulateGradU( VectorXs& F, const VectorXs& dx, const VectorXs& dv )
{
  assert( F.size() == m_x.size() );
  assert( dx.size() == dv.size() );
  assert( dx.size() == 0 || dx.size() == F.size() );
  
  // Accumulate all energy gradients
  if( dx.size() == 0 ) for( std::vector<Force*>::size_type i = 0; i < m_forces.size(); ++i ) {
    m_forces[i]->addGradEToTotal( m_x, m_v, m_interpolated_m, F );
  }
  else for( std::vector<Force*>::size_type i = 0; i < m_forces.size(); ++i ) {
    m_forces[i]->addGradEToTotal( m_x+dx, m_v+dv, m_interpolated_m, F );
  }
}

template<int DIM>
void TwoDScene<DIM>::accumulateddUdxdx( TripletXs& A, const VectorXs& dx, const VectorXs& dv )
{
  if( dx.size() == 0 ) for( std::vector<Force*>::size_type i = 0; i < m_forces.size(); ++i ) m_forces[i]->addHessXToTotal( m_x, m_v, m_interpolated_m, A );
  else                 for( std::vector<Force*>::size_type i = 0; i < m_forces.size(); ++i ) m_forces[i]->addHessXToTotal( m_x+dx, m_v+dv, m_interpolated_m, A );
}

template<int DIM>
void TwoDScene<DIM>::accumulateddUdxdv( TripletXs& A, const VectorXs& dx, const VectorXs& dv )
{
  if( dx.size() == 0 ) for( std::vector<Force*>::size_type i = 0; i < m_forces.size(); ++i ) m_forces[i]->addHessVToTotal( m_x, m_v, m_interpolated_m, A );
  else                 for( std::vector<Force*>::size_type i = 0; i < m_forces.size(); ++i ) m_forces[i]->addHessVToTotal( m_x+dx, m_v+dv, m_interpolated_m, A );
}

template<int DIM>
void TwoDScene<DIM>::accumulateExternalddUdxdx( TripletXs& A, const VectorXs& dx, const VectorXs& dv )
{
  
  
  if( dx.size() == 0 ) for( std::vector<Force*>::size_type i = 0; i < m_external_forces.size(); ++i ) m_external_forces[i]->addHessXToTotal( m_x, m_v, m_interpolated_m, A );
  else                 for( std::vector<Force*>::size_type i = 0; i < m_external_forces.size(); ++i ) m_external_forces[i]->addHessXToTotal( m_x+dx, m_v+dv, m_interpolated_m, A );
}

template<int DIM>
void TwoDScene<DIM>::accumulateExternalddUdxdv( TripletXs& A, const VectorXs& dx, const VectorXs& dv )
{
  
  
  if( dx.size() == 0 ) for( std::vector<Force*>::size_type i = 0; i < m_external_forces.size(); ++i ) m_external_forces[i]->addHessVToTotal( m_x, m_v, m_interpolated_m, A );
  else                 for( std::vector<Force*>::size_type i = 0; i < m_external_forces.size(); ++i ) m_external_forces[i]->addHessVToTotal( m_x+dx, m_v+dv, m_interpolated_m, A );
}

template<int DIM>
void TwoDScene<DIM>::copyState( const TwoDScene& otherscene )
{
  m_x = otherscene.m_x;
  m_base_x = otherscene.m_base_x;
  m_v = otherscene.m_v;
  m_m = otherscene.m_m;
  m_interpolated_m = otherscene.m_interpolated_m;
  m_rest_m = otherscene.m_rest_m;
  m_fixed = otherscene.m_fixed;
  m_edges = otherscene.m_edges;
}

template<int DIM>
bool TwoDScene<DIM>::useCtcd() const
{
  return m_parameters.use_ctcd;
}

template<int DIM>
void TwoDScene<DIM>::checkConsistency()
{
  assert( m_x.size() == m_v.size() );
  assert( m_x.size() == m_m.size() );
  if(m_massSpringSim) assert( m_x.size() == (int) ( DIM * m_fixed.size() ) );
  else assert( m_x.size() == (int) (4*m_fixed.size() - getNumFlows()) );
  assert( (m_x.array()==m_x.array()).all() );
  assert( (m_v.array()==m_v.array()).all() );
  assert( (m_m.array()==m_m.array()).all() );
  
  for( std::vector<std::pair<int,int> >::size_type i = 0; i < m_edges.size(); ++i )
  {
    assert( m_edges[i].first >= 0 );
    assert( m_edges[i].first < getNumParticles() );
    assert( m_edges[i].second >= 0 );
    assert( m_edges[i].second < getNumParticles() );
  }

  // TODO: Add more checks
}

template<int DIM>
std::vector<std::string>& TwoDScene<DIM>::getParticleTags()
{
  return m_particle_tags;
}

template<int DIM>
const std::vector<std::string>& TwoDScene<DIM>::getParticleTags() const
{
  return m_particle_tags;
}

template<int DIM>
const std::vector<scalar>& TwoDScene<DIM>::getFreeVolSummary() const
{
  return m_liquid_free;
}

template<int DIM>
const std::vector<scalar>& TwoDScene<DIM>::getHairVolSummary() const
{
  return m_liquid_on_hair;
}

template<int DIM>
void TwoDScene<DIM>::addVolSummary()
{
  if( m_parameters.volume_summary ) {
    scalar sum_liq = 0.0;
    scalar sum_res = 0.0;
    for(auto& flow : m_flows)
    {
      sum_liq += flow->computeTotalLiquidVol();
      sum_res += flow->computeTotalReservoirVol();
    }
    
    scalar sum_free = m_fluid_sim->computeTotalLiquidVol();
    
    m_liquid_on_hair.push_back(sum_liq + sum_res);
    m_liquid_free.push_back(sum_free);
    
    std::cout << "[Volume Summary:\t" << sum_liq << " + " << sum_res << " + " << sum_free << " = " << (sum_liq + sum_res + sum_free) << "\t]" << std::endl;
  }
}

template<int DIM>
void TwoDScene<DIM>::preCompute(const VectorXs& dx, const VectorXs& dv, const scalar& dt)
{
  assert( dx.size() == dv.size() );

  int nf = m_forces.size();
  if( dx.size() == 0 ) {
    threadutils::thread_pool::ParallelFor(0, nf, [&] (int i) {
      if(!m_forces[i]->isPrecomputationParallelized()) m_forces[i]->preCompute( m_x, m_v, m_interpolated_m, dt );
    });
    
    for(int i = 0; i < nf; ++i)
    {
      if(m_forces[i]->isPrecomputationParallelized()) m_forces[i]->preCompute( m_x, m_v, m_interpolated_m, dt );
    }
  }
  else  {
    VectorXs nx = m_x + dx;
    VectorXs nv = m_v + dv;
    threadutils::thread_pool::ParallelFor(0, nf, [&] (int i) {
      if(!m_forces[i]->isPrecomputationParallelized()) m_forces[i]->preCompute( nx, nv, m_interpolated_m, dt );
    });
    
    for(int i = 0; i < nf; ++i)
    {
      if(m_forces[i]->isPrecomputationParallelized()) m_forces[i]->preCompute( nx, nv, m_interpolated_m, dt );
    }
  }
}

template<int DIM>
void TwoDScene<DIM>::postCompute(const scalar& dt)
{
  const int nforces = (int) m_forces.size();
  for(int i = 0; i < nforces; ++i) m_forces[i]->postStepScene(dt);
}

template<int DIM>
void TwoDScene<DIM>::preComputeLocal(const VectorXs& dx, const VectorXs& dv, const scalar& dt)
{
  assert( dx.size() == dv.size() );
  
  int nfl = m_hair_internal_forces.size();
  if( dx.size() == 0 ) {
    threadutils::thread_pool::ParallelFor(0, nfl, [&] (int i) {
      for(Force* f : m_hair_internal_forces[i]) {
        if(!f->isPrecomputationParallelized()) f->preCompute( m_x, m_v, m_interpolated_m, dt );
      }
    });
    
    for(int i = 0; i < nfl; ++i)
    {
      for(Force* f : m_hair_internal_forces[i]) {
        if(f->isPrecomputationParallelized()) f->preCompute( m_x, m_v, m_interpolated_m, dt );
      }
    }
    
    for(Force* f : m_external_forces)
    {
      f->preCompute( m_x, m_v, m_interpolated_m, dt );
    }
  }
  else  {
    VectorXs nx = m_x + dx;
    VectorXs nv = m_v + dv;
    threadutils::thread_pool::ParallelFor(0, nfl, [&] (int i) {
      for(Force* f : m_hair_internal_forces[i]) {
        if(!f->isPrecomputationParallelized()) f->preCompute( nx, nv, m_interpolated_m, dt );
      }
    });
    
    for(int i = 0; i < nfl; ++i)
    {
      for(Force* f : m_hair_internal_forces[i]) {
        if(f->isPrecomputationParallelized()) f->preCompute( nx, nv, m_interpolated_m, dt );
      }
    }
    
    for(Force* f : m_external_forces)
    {
      f->preCompute( nx, nv, m_interpolated_m, dt );
    }
  }
}

template<int DIM>
void TwoDScene<DIM>::preComputeInterhair(const VectorXs& dx, const VectorXs& dv, const scalar& dt)
{
  assert( dx.size() == dv.size() );
  
  int nf = m_inter_hair_forces.size();
  if( dx.size() == 0 ) {
    threadutils::thread_pool::ParallelFor(0, nf, [&] (int i) {
      if(!m_inter_hair_forces[i]->isPrecomputationParallelized()) m_inter_hair_forces[i]->preCompute( m_x, m_v, m_interpolated_m, dt );
    });
    
    for(int i = 0; i < nf; ++i)
    {
      if(m_inter_hair_forces[i]->isPrecomputationParallelized()) m_inter_hair_forces[i]->preCompute( m_x, m_v, m_interpolated_m, dt );
    }
  }
  else  {
    VectorXs nx = m_x + dx;
    VectorXs nv = m_v + dv;
    threadutils::thread_pool::ParallelFor(0, nf, [&] (int i) {
      if(!m_inter_hair_forces[i]->isPrecomputationParallelized()) m_inter_hair_forces[i]->preCompute( nx, nv, m_interpolated_m, dt );
    });
    
    for(int i = 0; i < nf; ++i)
    {
      if(m_inter_hair_forces[i]->isPrecomputationParallelized()) m_inter_hair_forces[i]->preCompute( nx, nv, m_interpolated_m, dt );
    }
  }
}

template<int DIM>
void TwoDScene<DIM>::postPreprocess(VectorXs& lambda, VectorXs& lambda_v,
                                    TripletXs& J, TripletXs& Jv, TripletXs& Jxv, TripletXs& tildeK,
                                    TripletXs& stiffness, TripletXs& damping, VectorXs& Phi, VectorXs& Phiv, const VectorXs& dx, const VectorXs& dv, const scalar& dt)
{
  const int nf = (int) m_internal_forces.size();
  if(dx.size() == 0) {
    // for unparallelized forces
    threadutils::thread_pool::ParallelFor(0, nf, [&] (int i) {
      if(!m_internal_forces[i]->isParallelized()) m_internal_forces[i]->computeIntegrationVars( m_x, m_v, m_interpolated_m, lambda, lambda_v, J, Jv, Jxv, tildeK, stiffness, damping, Phi, Phiv, dt);
    });
    
    // for parallelized forces
    for(int i = 0; i < nf; ++i)
    {
      if(m_internal_forces[i]->isParallelized()) m_internal_forces[i]->computeIntegrationVars( m_x, m_v, m_interpolated_m, lambda, lambda_v, J, Jv, Jxv, tildeK, stiffness, damping, Phi, Phiv, dt);
    }
  } else {
    VectorXs nx = m_x + dx;
    VectorXs nv = m_v + dv;
    // for unparallelized forces
    threadutils::thread_pool::ParallelFor(0, nf, [&] (int i) {
      if(!m_internal_forces[i]->isParallelized()) m_internal_forces[i]->computeIntegrationVars( nx, nv, m_interpolated_m, lambda, lambda_v, J, Jv, Jxv, tildeK, stiffness, damping, Phi, Phiv, dt);
    });
    
    // for parallelized forces
    for(int i = 0; i < nf; ++i)
    {
      if(m_internal_forces[i]->isParallelized()) m_internal_forces[i]->computeIntegrationVars( nx, nv, m_interpolated_m, lambda, lambda_v, J, Jv, Jxv, tildeK, stiffness, damping, Phi, Phiv, dt);
    }
  }
}

template<int DIM>
void TwoDScene<DIM>::interhairPostPreprocess(VectorXs& lambda, VectorXs& lambda_v,
                                             TripletXs& J, TripletXs& Jv, TripletXs& Jxv, TripletXs& tildeK,
                                             TripletXs& stiffness, TripletXs& damping, VectorXs& Phi, VectorXs& Phiv, const VectorXs& dx, const VectorXs& dv, const scalar& dt)
{
  const int nf = (int) m_inter_hair_forces.size();
  if(dx.size() == 0) {
    // for unparallelized forces
    threadutils::thread_pool::ParallelFor(0, nf, [&] (int i) {
      if(!m_inter_hair_forces[i]->isParallelized()) m_inter_hair_forces[i]->computeIntegrationVars( m_x, m_v, m_interpolated_m, lambda, lambda_v, J, Jv, Jxv, tildeK, stiffness, damping, Phi, Phiv, dt);
    });
    
    // for parallelized forces
    for(int i = 0; i < nf; ++i)
    {
      if(m_inter_hair_forces[i]->isParallelized()) m_inter_hair_forces[i]->computeIntegrationVars( m_x, m_v, m_interpolated_m, lambda, lambda_v, J, Jv, Jxv, tildeK, stiffness, damping, Phi, Phiv, dt);
    }
  } else {
    VectorXs nx = m_x + dx;
    VectorXs nv = m_v + dv;
    // for unparallelized forces
    threadutils::thread_pool::ParallelFor(0, nf, [&] (int i) {
      if(!m_inter_hair_forces[i]->isParallelized()) m_inter_hair_forces[i]->computeIntegrationVars( nx, nv, m_interpolated_m, lambda, lambda_v, J, Jv, Jxv, tildeK, stiffness, damping, Phi, Phiv, dt);
    });
    
    // for parallelized forces
    for(int i = 0; i < nf; ++i)
    {
      if(m_inter_hair_forces[i]->isParallelized()) m_inter_hair_forces[i]->computeIntegrationVars( nx, nv, m_interpolated_m, lambda, lambda_v, J, Jv, Jxv, tildeK, stiffness, damping, Phi, Phiv, dt);
    }
  }
}

template<int DIM>
void TwoDScene<DIM>::localPostPreprocess(VectorXs& lambda, VectorXs& lambda_v,
                                    TripletXs& J, TripletXs& Jv, TripletXs& Jxv, TripletXs& tildeK,
                                    TripletXs& stiffness, TripletXs& damping, VectorXs& Phi, VectorXs& Phiv, const VectorXs& dx, const VectorXs& dv, const scalar& dt)
{

  const int nf = (int) m_hair_internal_forces.size();
  
  if(dx.size() == 0) {
    // for unparallelized forces
    threadutils::thread_pool::ParallelFor(0, nf, [&] (int i) {
      for(Force* f : m_hair_internal_forces[i]) {
        if(!f->isParallelized()) f->computeIntegrationVars( m_x, m_v, m_interpolated_m, lambda, lambda_v, J, Jv, Jxv, tildeK, stiffness, damping, Phi, Phiv, dt);
      }
    });
    
    // for parallelized forces
    for(int i = 0; i < nf; ++i)
    {
      for(Force* f : m_hair_internal_forces[i]) {
        if(f->isParallelized()) f->computeIntegrationVars( m_x, m_v, m_interpolated_m, lambda, lambda_v, J, Jv, Jxv, tildeK, stiffness, damping, Phi, Phiv, dt);
      }
    }
  } else {
    VectorXs nx = m_x + dx;
    VectorXs nv = m_v + dv;
    // for unparallelized forces
    threadutils::thread_pool::ParallelFor(0, nf, [&] (int i) {
      for(Force* f : m_hair_internal_forces[i]) {
        if(!f->isParallelized()) f->computeIntegrationVars( nx, nv, m_interpolated_m, lambda, lambda_v, J, Jv, Jxv, tildeK, stiffness, damping, Phi, Phiv, dt);
      }
    });
    
    // for parallelized forces
    for(int i = 0; i < nf; ++i)
    {
      for(Force* f : m_hair_internal_forces[i]) {
        if(f->isParallelized()) f->computeIntegrationVars( nx, nv, m_interpolated_m, lambda, lambda_v, J, Jv, Jxv, tildeK, stiffness, damping, Phi, Phiv, dt);
      }
    }
  }
}

template<int DIM>
void TwoDScene<DIM>::getAffectedVars(int icol, std::unordered_set<int>& affected)
{
  for( std::vector<Force*>::size_type i = 0; i < m_forces.size(); ++i ) m_forces[i]->getAffectedVars(icol, affected);
}

template<int DIM>
void TwoDScene<DIM>::getAffectedForces(int icol, std::vector<Force*>& forces)
{
  for( std::vector<Force*>::size_type i = 0; i < m_forces.size(); ++i ) if(m_forces[i]->isContained(icol)) forces.push_back(m_forces[i]);
}

template<int DIM>
bool TwoDScene<DIM>::viscositySolve() const
{
  return m_parameters.viscous_solve;
}

template<int DIM>
bool TwoDScene<DIM>::drippingNear() const
{
  return m_parameters.drippingnear;
}

template<int DIM>
bool TwoDScene<DIM>::drippingFar() const
{
  return m_parameters.drippingfar;
}

template<int DIM>
bool TwoDScene<DIM>::drippingMiddle() const
{
  return m_parameters.drippingmiddle;
}


template<int DIM>
bool TwoDScene<DIM>::isFlowOnly() const
{
  return true;
}

template<int DIM>
bool TwoDScene<DIM>::doVolSummary() const
{
  return m_parameters.volume_summary;
}

template<int DIM>
void TwoDScene<DIM>::notifyGlobalIntegrator()
{
  if(m_polygonal_cohesion) {
    m_polygonal_cohesion->setUseDecoupledForce(false);
  }
}

template<int DIM>
void TwoDScene<DIM>::notifyFullIntegrator()
{
  if(m_polygonal_cohesion) {
    m_polygonal_cohesion->setUseParticlePOEMap(false);
  }
}

template<int DIM>
void TwoDScene<DIM>::notifyPartialIntegrator()
{
  if(m_polygonal_cohesion) {
    m_polygonal_cohesion->setUseParticlePOEMap(true);
  }
}

template<int DIM>
void TwoDScene<DIM>::notifyLocalIntegrator()
{
  if(m_polygonal_cohesion) {
    m_polygonal_cohesion->setUseDecoupledForce(true);
  }
}

template<int DIM>
void TwoDScene<DIM>::updateEdgeRadii()
{
  if( m_massSpringSim ){
    int ne = getNumEdges();
    for(int i = 0; i < ne; ++i)
    {
      auto& e = m_edges[i];
      scalar l = (m_x.segment<DIM>(e.first * DIM) - m_x.segment<DIM>(e.second * DIM)).norm();
      m_edge_radii(i) = std::max(m_edge_rest_radii(i) * 0.5, (1.0 - m_edge_poisson_ratio(i) * (l - m_edge_rest_length(i)) / m_edge_rest_length(i)) * m_edge_rest_radii(i));
    }
  }
  else{
    std::cerr<< "edge radii computed via strand parameters in TwoDScene<3>::computeMassesAndRadiiFromStrands(), exiting." << std::endl; exit(1);
  }
}

template<int DIM>
void TwoDScene<DIM>::updateMassWithLiquid(int idx, const scalar& liquid_mass, const scalar& liquid_radii)
{
  if( m_massSpringSim || isTip( idx ) ){
    for(int r = 0; r < DIM; ++r){
      m_m( getDof( idx ) + r) = m_rest_m( getDof( idx ) + r) + std::max( 0.0, liquid_mass ) ;
    }
  }
  else{
    for(int r = 0; r < 4; ++r){
      if(r == 3) {
        m_m( getDof( idx ) + r) = m_rest_m( getDof( idx ) + r) + std::max( 0.0, liquid_mass ) * (liquid_radii * liquid_radii + m_radii(idx) * m_radii(idx)) * 0.5;
      } else {
        m_m( getDof( idx ) + r) = m_rest_m( getDof( idx ) + r) + std::max( 0.0, liquid_mass );
      }
    }
  }
}

template<int DIM>
const scalar& TwoDScene<DIM>::getLiquidShell() const
{
  return m_parameters.regularizer_shell;
}

template<int DIM>
scalar& TwoDScene<DIM>::getLiquidShell()
{
  return m_parameters.regularizer_shell;
}

template<int DIM>
const VectorXs& TwoDScene<DIM>::getLiquidRestLength() const
{
  return m_edge_rest_length;
}

template<int DIM>
VectorXs& TwoDScene<DIM>::getLiquidRestLength()
{
  return m_edge_rest_length;
}

template<int DIM>
void TwoDScene<DIM>::saveFluidPressure( const std::string& szfn )
{
  if(!m_fluid_sim) return;
  
  m_fluid_sim->save_pressure(szfn);
}

template<int DIM>
Vectors<DIM> TwoDScene<DIM>::computeHairParticleMomentum() const
{
//[H]  todo, unclear if this has to change for twist
  int nf = getNumFlows();
  Vectors<DIM> p = Vectors<DIM>::Zero();
  
  for(int i = 0; i < nf; ++i)
  {
    p += m_flows[i]->computeHairLiquidMomentum(m_v);
  }
  
  return p;
}

template<int DIM>
Vector3s TwoDScene<DIM>::computeHairParticleAngularMomentum() const
{
//[H]  todo, unclear if this has to change for twist
  int nf = getNumFlows();
  Vector3s p = Vector3s::Zero();

  for(int i = 0; i < nf; ++i)
  {
    p += m_flows[i]->computeHairLiquidAngularMomentum(m_x, m_v, m_fluid_sim);
  }
  
  return p;
}

template<int DIM>
scalar TwoDScene<DIM>::computeHairParticleKineticEnergy() const
{
  int nf = getNumFlows();
  scalar E = 0.0;
  
  for(int i = 0; i < nf; ++i)
  {
    E += m_flows[i]->computeHairLiquidEnergy(m_v);
  }
  
  return E;
}

template<int DIM>
void TwoDScene<DIM>::categorizeForces()
{
  m_external_forces.clear();
  m_internal_forces.clear();
  m_inter_hair_forces.clear();
  m_hair_internal_forces.resize(m_flows.size());
  
  for( size_t i = 0; i < m_flows.size(); ++i)
  {
    m_hair_internal_forces[i].resize(0);
  }
  
  for( size_t i = 0; i < m_forces.size(); ++i )
  {
    if(m_forces[i]->isExternal()) {
      m_external_forces.push_back(m_forces[i]);
    } else {
      m_internal_forces.push_back(m_forces[i]);
      
      if(!m_forces[i]->isInterHair()) {
        int hidx = m_forces[i]->getAffectedHair(m_particle_to_hairs);
        if(hidx < 0) continue;
        m_hair_internal_forces[hidx].push_back(m_forces[i]);
      } else {
        m_inter_hair_forces.push_back(m_forces[i]);
      }
    }
  }
}

template<int DIM>
void TwoDScene<DIM>::updateNumConstraintsLocal(int& num_constraint_pos_, int& num_constraint_vel_, int& num_J_, int& num_Jv_, int& num_Jxv_, int& num_tildeK_)
{
  m_constraint_idx = Vector6i::Zero();
  
  int nfl = getNumFlows();
  for(int i = 0; i < nfl; ++i)
  {
    const std::vector< Force* >& hair_forces = m_hair_internal_forces[i];
    int nhair_forces = (int) hair_forces.size();
    
    Vector6i constraint_start = m_constraint_idx;
    
    for(int j = 0; j < nhair_forces; ++j)
    {
      int num_pos = hair_forces[j]->numConstraintPos();
      int num_vel = hair_forces[j]->numConstraintVel();
      int num_J = hair_forces[j]->numJ();
      int num_Jv = hair_forces[j]->numJv();
      int num_Jxv = hair_forces[j]->numJxv();
      int num_TildeK = hair_forces[j]->numTildeK();
      
      hair_forces[j]->setInternalIndex(m_constraint_idx(0), m_constraint_idx(1), m_constraint_idx(2), m_constraint_idx(3), m_constraint_idx(4), m_constraint_idx(5));
      m_constraint_idx(0) += num_pos;
      m_constraint_idx(1) += num_vel;
      m_constraint_idx(2) += num_J;
      m_constraint_idx(3) += num_Jv;
      m_constraint_idx(4) += num_Jxv;
      m_constraint_idx(5) += num_TildeK;
    }
    
    Vector6i num_constraints = m_constraint_idx - constraint_start;
    
    m_flows[i]->setConstraintParameters(constraint_start, num_constraints);
  }

  num_constraint_pos_ = m_constraint_idx(0);
  num_constraint_vel_ = m_constraint_idx(1);
  num_J_ = m_constraint_idx(2);
  num_Jv_ = m_constraint_idx(3);
  num_Jxv_ = m_constraint_idx(4);
  num_tildeK_ = m_constraint_idx(5);
}

template<int DIM>
void TwoDScene<DIM>::updateNumConstraintsInterHair(int& num_constraint_pos_, int& num_constraint_vel_, int& num_J_, int& num_Jv_, int& num_Jxv_, int& num_tildeK_,
                                                   Vector6i& interhair_param, Vector6i& interhair_num)
{
  interhair_param = m_constraint_idx;
  
  for( std::vector<Force*>::size_type i = 0; i < m_inter_hair_forces.size(); ++i )
  {
    int num_pos = m_inter_hair_forces[i]->numConstraintPos();
    int num_vel = m_inter_hair_forces[i]->numConstraintVel();
    int num_J = m_inter_hair_forces[i]->numJ();
    int num_Jv = m_inter_hair_forces[i]->numJv();
    int num_Jxv = m_inter_hair_forces[i]->numJxv();
    int num_TildeK = m_inter_hair_forces[i]->numTildeK();
    
    m_inter_hair_forces[i]->setInternalIndex(m_constraint_idx(0), m_constraint_idx(1), m_constraint_idx(2), m_constraint_idx(3), m_constraint_idx(4), m_constraint_idx(5));
    m_constraint_idx(0) += num_pos;
    m_constraint_idx(1) += num_vel;
    m_constraint_idx(2) += num_J;
    m_constraint_idx(3) += num_Jv;
    m_constraint_idx(4) += num_Jxv;
    m_constraint_idx(5) += num_TildeK;
  }
  
  num_constraint_pos_ = m_constraint_idx(0);
  num_constraint_vel_ = m_constraint_idx(1);
  num_J_ = m_constraint_idx(2);
  num_Jv_ = m_constraint_idx(3);
  num_Jxv_ = m_constraint_idx(4);
  num_tildeK_ = m_constraint_idx(5);
  
  interhair_num = m_constraint_idx - interhair_param;
}

template<int DIM>
void TwoDScene<DIM>::updateNumConstraints(int& num_constraint_pos_, int& num_constraint_vel_, int& num_J_, int& num_Jv_, int& num_Jxv_, int& num_tildeK_,
                                          Vector6i& interhair_param, Vector6i& interhair_num)
{
  m_constraint_idx = Vector6i::Zero();
  
  int nfl = getNumFlows();
  for(int i = 0; i < nfl; ++i)
  {
    const std::vector< Force* >& hair_forces = m_hair_internal_forces[i];
    int nhair_forces = (int) hair_forces.size();
    
    Vector6i constraint_start = m_constraint_idx;
    
    for(int j = 0; j < nhair_forces; ++j)
    {
      int num_pos = hair_forces[j]->numConstraintPos();
      int num_vel = hair_forces[j]->numConstraintVel();
      int num_J = hair_forces[j]->numJ();
      int num_Jv = hair_forces[j]->numJv();
      int num_Jxv = hair_forces[j]->numJxv();
      int num_TildeK = hair_forces[j]->numTildeK();
      
      hair_forces[j]->setInternalIndex(m_constraint_idx(0), m_constraint_idx(1), m_constraint_idx(2), m_constraint_idx(3), m_constraint_idx(4), m_constraint_idx(5));
      m_constraint_idx(0) += num_pos;
      m_constraint_idx(1) += num_vel;
      m_constraint_idx(2) += num_J;
      m_constraint_idx(3) += num_Jv;
      m_constraint_idx(4) += num_Jxv;
      m_constraint_idx(5) += num_TildeK;
    }
    
    Vector6i num_constraints = m_constraint_idx - constraint_start;
    
    m_flows[i]->setConstraintParameters(constraint_start, num_constraints);
  }
  
  interhair_param = m_constraint_idx;
  
  for( std::vector<Force*>::size_type i = 0; i < m_inter_hair_forces.size(); ++i )
  {
    int num_pos = m_inter_hair_forces[i]->numConstraintPos();
    int num_vel = m_inter_hair_forces[i]->numConstraintVel();
    int num_J = m_inter_hair_forces[i]->numJ();
    int num_Jv = m_inter_hair_forces[i]->numJv();
    int num_Jxv = m_inter_hair_forces[i]->numJxv();
    int num_TildeK = m_inter_hair_forces[i]->numTildeK();
    
    m_inter_hair_forces[i]->setInternalIndex(m_constraint_idx(0), m_constraint_idx(1), m_constraint_idx(2), m_constraint_idx(3), m_constraint_idx(4), m_constraint_idx(5));
    m_constraint_idx(0) += num_pos;
    m_constraint_idx(1) += num_vel;
    m_constraint_idx(2) += num_J;
    m_constraint_idx(3) += num_Jv;
    m_constraint_idx(4) += num_Jxv;
    m_constraint_idx(5) += num_TildeK;
  }
  
  num_constraint_pos_ = m_constraint_idx(0);
  num_constraint_vel_ = m_constraint_idx(1);
  num_J_ = m_constraint_idx(2);
  num_Jv_ = m_constraint_idx(3);
  num_Jxv_ = m_constraint_idx(4);
  num_tildeK_ = m_constraint_idx(5);
  
  interhair_num = m_constraint_idx - interhair_param;
}

template<int DIM>
void TwoDScene<DIM>::storeLambda(const VectorXs& lambda, const VectorXs& lambda_v)
{
  int nf = m_internal_forces.size();
  threadutils::thread_pool::ParallelFor(0, nf, [&] (int i) {
    m_internal_forces[i]->storeLambda(lambda, lambda_v);
  });
}

template<>
Vector2s TwoDScene<2>::computeLiquidParticleMomentum() const
{
  return ((FluidSim2D*) m_fluid_sim)->computeParticleMomentum();
}

template<>
Vector3s TwoDScene<2>::computeLiquidParticleAngularMomentum() const
{
  Vector3s v(0,0,0);
  v(2) = ((FluidSim2D*) m_fluid_sim)->computeParticleAngularMomentum();
  return v;
}

template<>
Vector3s TwoDScene<2>::computeLiquidGridAngularMomentum() const
{
  Vector3s v(0,0,0);
  v(2) = ((FluidSim2D*) m_fluid_sim)->computeParticleGridAngularMomentum();
  return v;
}

template<>
Vector3s TwoDScene<2>::computeHairGridAngularMomentum() const
{
  Vector3s v(0,0,0);
  v(2) = ((FluidSim2D*) m_fluid_sim)->computeHairGridAngularMomentum();
  return v;
}

template<>
Vector3s TwoDScene<2>::computeReweightedHairGridAngularMomentum() const
{
  Vector3s v(0,0,0);
  v(2) = ((FluidSim2D*) m_fluid_sim)->computeReweightedHairGridAngularMomentum();
  return v;
}

template<>
Vector3s TwoDScene<2>::computeReweightedLiquidGridAngularMomentum() const
{
  Vector3s v(0,0,0);
  v(2) = ((FluidSim2D*) m_fluid_sim)->computeReweightedParticleGridAngularMomentum();
  return v;
}

template<>
Vector3s TwoDScene<2>::computeCombinedGridAngularMomentum() const
{
  Vector3s v(0,0,0);
  v(2) = ((FluidSim2D*) m_fluid_sim)->computeCombinedGridAngularMomentum();
  return v;
}

template<>
scalar TwoDScene<2>::computeLiquidParticleKineticEnergy() const
{
  return ((FluidSim2D*) m_fluid_sim)->computeParticleKineticEnergy();
}

template<>
Vectors<2> TwoDScene<2>::computeLiquidGridMomentum() const
{
  return ((FluidSim2D*) m_fluid_sim)->computeParticleGridMomentum();
}

template<>
Vectors<2> TwoDScene<2>::computeReweightedLiquidGridMomentum() const
{
  return ((FluidSim2D*) m_fluid_sim)->computeReweightedParticleGridMomentum();
}

template<>
scalar TwoDScene<2>::computeLiquidGridKineticEnergy() const
{
  return ((FluidSim2D*) m_fluid_sim)->computeParticleKineticEnergy();
}

template<>
Vectors<2> TwoDScene<2>::computeHairGridMomentum() const
{
  return ((FluidSim2D*) m_fluid_sim)->computeHairGridMomentum();
}

template<>
Vectors<2> TwoDScene<2>::computeReweightedHairGridMomentum() const
{
  return ((FluidSim2D*) m_fluid_sim)->computeReweightedHairGridMomentum();
}

template<>
scalar TwoDScene<2>::computeHairGridKineticEnergy() const
{
  return ((FluidSim2D*) m_fluid_sim)->computeHairGridKineticEnergy();
}

template<>
Vectors<2> TwoDScene<2>::computeCombinedGridMomentum() const
{
  return ((FluidSim2D*) m_fluid_sim)->computeCombinedGridMomentum();
}

template<>
scalar TwoDScene<2>::computeCombinedGridKineticEnergy() const
{
  return ((FluidSim2D*) m_fluid_sim)->computeCombinedGridKineticEnergy();
}

template<>
Vectors<2> TwoDScene<2>::computeParticleWeightedCombinedGridMomentum() const
{
  return ((FluidSim2D*) m_fluid_sim)->computeParticleWeightedCombinedGridMomentum();
}

template<>
Vectors<2> TwoDScene<2>::computeHairWeightedCombinedGridMomentum() const
{
  return ((FluidSim2D*) m_fluid_sim)->computeHairWeightedCombinedGridMomentum();
}

template<>
Vector3s TwoDScene<2>::computeParticleWeightedCombinedGridAngularMomentum() const
{
  return Vector3s(0,0,((FluidSim2D*) m_fluid_sim)->computeParticleWeightedCombinedGridAngularMomentum());
}

template<>
Vector3s TwoDScene<2>::computeHairWeightedCombinedGridAngularMomentum() const
{
  return Vector3s(0,0,((FluidSim2D*) m_fluid_sim)->computeHairWeightedCombinedGridAngularMomentum());
}

template<>
scalar TwoDScene<2>::computeLiquidOverallDivergence() const
{
  return ((FluidSim2D*) m_fluid_sim)->computeOverallDivergence();
}

template<int DIM>
Vectors<DIM> TwoDScene<DIM>::computeHairDrag()
{
//[H]  todo, unclear if this has to change for twist

  int nf = getNumFlows();
  Vectors<DIM> p = Vectors<DIM>::Zero();

  for(int i = 0; i < nf; ++i)
  {
    p += m_flows[i]->computeHairDragForce(m_fluid_drag_buffer);
  }
  
  return p;
}

template<>
scalar TwoDScene<3>::computeLiquidOverallDivergence() const
{
  return ((FluidSim3D*) m_fluid_sim)->computeOverallDivergence();
}

template<>
Vectors<3> TwoDScene<3>::computeLiquidParticleMomentum() const
{
  return ((FluidSim3D*) m_fluid_sim)->computeParticleMomentum();
}

template<>
Vector3s TwoDScene<3>::computeLiquidParticleAngularMomentum() const
{
  return ((FluidSim3D*) m_fluid_sim)->computeParticleAngularMomentum();
}

template<>
Vector3s TwoDScene<3>::computeLiquidGridAngularMomentum() const
{
  return ((FluidSim3D*) m_fluid_sim)->computeParticleGridAngularMomentum();
}

template<>
scalar TwoDScene<3>::computeLiquidParticleKineticEnergy() const
{
  return ((FluidSim3D*) m_fluid_sim)->computeParticleKineticEnergy();
}

template<>
scalar TwoDScene<3>::computeLiquidGridKineticEnergy() const
{
  return ((FluidSim3D*) m_fluid_sim)->computeParticleKineticEnergy();
}

template<>
Vectors<3> TwoDScene<3>::computeCombinedGridMomentum() const
{
  return ((FluidSim3D*) m_fluid_sim)->computeCombinedGridMomentum();
}

template<>
scalar TwoDScene<3>::computeCombinedGridKineticEnergy() const
{
  return ((FluidSim3D*) m_fluid_sim)->computeCombinedGridKineticEnergy();
}

template<>
Vector3s TwoDScene<3>::computeReweightedLiquidGridAngularMomentum() const
{
  return ((FluidSim3D*) m_fluid_sim)->computeReweightedParticleGridAngularMomentum();
}

template<>
Vector3s TwoDScene<3>::computeCombinedGridAngularMomentum() const
{
  return ((FluidSim3D*) m_fluid_sim)->computeCombinedGridAngularMomentum();
}

template<int DIM>
bool TwoDScene<DIM>::confirmNotNeighborEdgesOnSameHair( const int& efirst, const int& esecond ) const
{
  if( DIM == 2 ) return false;

  if( m_edge_to_hair[efirst] != m_edge_to_hair[esecond] ) return true;
  else if( std::abs(efirst - esecond) < 3 ) return false;
  return true;
}

template<int DIM>
void TwoDScene<DIM>::setVertToDoFMap( const std::vector<int>& vert_to_dof, const VectorXi& dofs_to_vars, const std::vector<bool>& tipVerts, const VectorXi& dof_to_vert )
{
  if( m_massSpringSim ){
    std::cerr<< "setVertToDoFMap not defined for 2D scene, exiting." << std::endl; exit(1);
  }

  m_particle_to_dofs = vert_to_dof;
  m_dofs_to_component = dofs_to_vars;
  m_is_strand_tip = tipVerts;
  m_dofs_to_particle = dof_to_vert;
}

template<int DIM>
void TwoDScene<DIM>::computeMassesAndRadiiFromStrands()
{
  if( m_massSpringSim ){
    std::cerr<< "computeMassesAndRadiiFromStrands not defined for Mass-Spring scene, exiting." << std::endl; exit(1);
  }
  m_strands.resize( m_num_strands );

  int nStrand = 0;
  for( int f = 0; f < m_forces.size(); ++f )
  {
    StrandForce* const strand = dynamic_cast<StrandForce*>( m_forces[f] );
    if( strand == NULL ) continue;
    m_strands[nStrand] = strand;


    if( strand->m_strandParams->m_straightHairs != 1.0 ){
      Vec2Array& kappas = strand->alterRestKappas();
      for( int k = 0; k < kappas.size(); ++k ){
        kappas[k] *= strand->m_strandParams->m_straightHairs;
      }
    }

    for( int v = 0; v < strand->getNumVertices(); ++v ){
      const int globalVtx = strand->m_verts[v];
      const int globalEdx = globalVtx - nStrand;
      const int globalDof = getDof( globalVtx );
      const scalar r = strand->m_strandParams->getRadius( v, strand->getNumVertices() );

      m_radii[ globalVtx ] = r;
      m_m.segment<DIM>( globalDof ).setConstant( strand->m_vertexMasses[v] );

      if( v < strand->getNumEdges() ){
        m_edge_to_hair[globalEdx] = nStrand;
        // Edge radius, edge's should be indexed the same as 
        m_edge_radii[ globalEdx ] = r;
        m_edge_rest_radii[ globalEdx ] = r;

        // Twist Mass (Second moment of inertia * length)
        const scalar mass = strand->m_strandParams->m_density * M_PI * r * r * strand->m_restLengths[v];
        scalar vtm = 0.25 * mass * 2 * r * r;
        m_m[ globalDof + 3 ] = vtm;
      }
    }
    ++nStrand;
  }

  m_rest_m = m_m;
  m_interpolated_m = m_m;
  assert( nStrand == m_num_strands );
}

template<int DIM>
const WetHairParameter& TwoDScene<DIM>::getParameter() const
{
  return m_parameters;
}

template<int DIM>
WetHairParameter& TwoDScene<DIM>::getParameter()
{
  return m_parameters;
}

// explicit instantiations at bottom
template class TwoDScene<2>;
template class TwoDScene<3>;
