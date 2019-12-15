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


#include "HairFlow.h"
#include "TwoDScene.h"
#include "MathUtilities.h"
#include "fluidsim2D.h"
#include "array2_utils.h"

#include <iostream>
#include <fstream>

using namespace mathutils;

//#define NO_COHESION_FORCE

template<int DIM>
int HairFlow<DIM>::m_flow_counter = 0;

template<int DIM>
void HairFlow<DIM>::adjustVolumeGlobal(const scalar& prop)
{
  for(int i = 0; i < m_eta.size(); ++i)
  {
    m_eta(i) = sqrt(std::max(0.0, prop) * m_eta(i) * (m_eta(i) + m_rad_vec(i) * 2.0) + m_rad_vec(i) * m_rad_vec(i)) - m_rad_vec(i);
  }
}

template<int DIM>
int HairFlow<DIM>::index() const
{
  return m_flow_index;
}

template<int DIM>
int HairFlow<DIM>::size() const
{
  return m_global_to_local.size();
}

template<int DIM>
int HairFlow<DIM>::find(int idx_global) const
{
  auto itr = m_global_to_local.find(idx_global);
  if(itr == m_global_to_local.end()) {
    return -1;
  } else {
    return itr->second;
  }
}

template<int DIM>
HairFlow<DIM>::HairFlow(TwoDScene<DIM>* parent, const std::vector<int>& involved_particles, const VectorXs& eta, const std::vector<unsigned char>& particle_state)
: m_parent(parent), m_particle_indices(involved_particles), m_eta(eta), m_particle_state(particle_state), m_flow_index(m_flow_counter++)
{
  m_porosity.resize(m_eta.size());
  m_liquid_phi.resize(m_eta.size());
  
  auto& radii = m_parent->getRadii();
  
  int np = m_porosity.size();
  for(int i = 0; i < np; ++i)
  {
    scalar r = radii(m_particle_indices[i]);
    m_porosity(i) = 1.0 - r * r / ((m_eta(i) + r) * (m_eta(i) + r));
  }
  
  m_avg_eta = m_eta;
}

template<int DIM>
const VectorXs& HairFlow<DIM>::getVelocity() const
{
  return m_u;
}

template<int DIM>
VectorXs& HairFlow<DIM>::getVelocity()
{
  return m_u;
}

template<int DIM>
const std::vector<unsigned char>& HairFlow<DIM>::getState() const
{
  return m_particle_state;
}

template<int DIM>
const std::vector< std::vector<int> >& HairFlow<DIM>::getEdgeBridges() const
{
  return m_edge_bridges;
}

template<int DIM>
std::vector< std::vector<int> >& HairFlow<DIM>::getEdgeBridges()
{
  return m_edge_bridges;
}

template<int DIM>
const std::vector< int >& HairFlow<DIM>::getEdgeIndices() const
{
  return m_edge_indices;
}

template<int DIM>
std::unordered_map<int, int>& HairFlow<DIM>::getGlobalToLocal()
{
  return m_global_to_local;
}

template<int DIM>
const std::unordered_map<int, int>& HairFlow<DIM>::getGlobalToLocal() const
{
  return m_global_to_local;
}

template<int DIM>
const std::vector< int >& HairFlow<DIM>::getParticleIndices() const
{
  return m_particle_indices;
}

template<int DIM>
const std::vector< std::pair<int, int> >& HairFlow<DIM>::getLocalEdges() const
{
  return m_internal_edges;
}

template<int DIM>
const std::vector< std::pair<int, int> >& HairFlow<DIM>::getGlobalEdges() const
{
  return m_global_edges;
}

template<int DIM>
const scalar& HairFlow<DIM>::getAvgAreaE() const
{
  return m_avg_area_e;
}

template<int DIM>
const scalar& HairFlow<DIM>::getMinAreaE() const
{
  return m_min_area_e;
}

template<int DIM>
const scalar& HairFlow<DIM>::getMaxAreaE() const
{
  return m_max_area_e;
}

template<int DIM>
const scalar& HairFlow<DIM>::getMinAreaV() const
{
  return m_min_area_v;
}

template<int DIM>
const scalar& HairFlow<DIM>::getMaxAreaV() const
{
  return m_max_area_v;
}

template<int DIM>
const scalar& HairFlow<DIM>::getMinEta() const
{
  return m_min_eta;
}

template<int DIM>
const scalar& HairFlow<DIM>::getMaxEta() const
{
  return m_max_eta;
}

template<int DIM>
const VectorXs& HairFlow<DIM>::getEta() const
{
  return m_eta;
}

template<int DIM>
VectorXs& HairFlow<DIM>::getEta()
{
  return m_eta;
}

template<int DIM>
const VectorXs& HairFlow<DIM>::getAvgEta() const
{
  return m_avg_eta;
}

template<int DIM>
VectorXs& HairFlow<DIM>::getAvgEta()
{
  return m_avg_eta;
}

template<int DIM>
const VectorXs& HairFlow<DIM>::getPorosity() const
{
  return m_porosity;
}

template<int DIM>
VectorXs& HairFlow<DIM>::getPorosity()
{
  return m_porosity;
}

template<int DIM>
const VectorXs& HairFlow<DIM>::getAreaV() const
{
  return m_area_v;
}

template<int DIM>
VectorXs& HairFlow<DIM>::getAreaV()
{
  return m_area_v;
}

template<int DIM>
const VectorXs& HairFlow<DIM>::getAreaVHair() const
{
  return m_area_v_hair;
}

template<int DIM>
VectorXs& HairFlow<DIM>::getAreaVHair()
{
  return m_area_v_hair;
}

template<int DIM>
const VectorXs& HairFlow<DIM>::getAreaE() const
{
  return m_area_e;
}

template<int DIM>
VectorXs& HairFlow<DIM>::getAreaE()
{
  return m_area_e;
}

template<int DIM>
const VectorXs& HairFlow<DIM>::getNormalV() const
{
  return m_normal_v;
}

template<int DIM>
VectorXs& HairFlow<DIM>::getNormalV()
{
  return m_normal_v;
}

template<int DIM>
const VectorXs& HairFlow<DIM>::getNormalE() const
{
  return m_normal_e;
}

template<int DIM>
VectorXs& HairFlow<DIM>::getNormalE()
{
  return m_normal_e;
}

template<int DIM>
const VectorXs& HairFlow<DIM>::getRadiiE() const
{
  return m_edge_rad_vec;
}

template<int DIM>
VectorXs& HairFlow<DIM>::getRadiiE()
{
  return m_edge_rad_vec;
}

template<int DIM>
const VectorXs& HairFlow<DIM>::getRadiiV() const
{
  return m_rad_vec;
}

template<int DIM>
VectorXs& HairFlow<DIM>::getRadiiV()
{
  return m_rad_vec;
}

template<int DIM>
const MatrixXs& HairFlow<DIM>::getTangentE() const
{
  return m_dir_f;
}

template<int DIM>
MatrixXs& HairFlow<DIM>::getTangentE()
{
  return m_dir_f;
}


template<int DIM>
const MatrixXs& HairFlow<DIM>::getTangentV() const
{
  return m_dir_v;
}

template<int DIM>
MatrixXs& HairFlow<DIM>::getTangentV()
{
  return m_dir_v;
}

template<int DIM>
const MatrixXs& HairFlow<DIM>::getActualVertexVelocity() const
{
  return m_actual_u_v;
}

template<int DIM>
MatrixXs& HairFlow<DIM>::getActualVertexVelocity()
{
  return m_actual_u_v;
}

template<int DIM>
void HairFlow<DIM>::getAffectedVars( int colidx, std::unordered_set<int>& vars )
{
  int ip = m_parent->getVertFromDof(colidx);
  if(  m_parent->getComponent(colidx) == 0 && find(ip) != -1)
  {
    vars.insert(colidx);
  }
}

template<int DIM>
bool HairFlow<DIM>::isContained( int colidx )
{
  if( m_parent->getComponent(colidx) == DIM ) return false; //[H] todo, check if twist matters here or not
  int ip = m_parent->getVertFromDof(colidx);
  return find(ip) != -1;
}

template<int DIM>
void HairFlow<DIM>::read(const scalar* data)
{
  int neta = m_eta.size();
  int k = 0;
  for(int i = 0; i < neta; ++i)
  {
    m_eta(i) = data[k++];
  }
  int nu = m_u.size();
  for(int i = 0; i < nu; ++i)
  {
    m_u(i) = data[k++];
  }
}

template<int DIM>
void HairFlow<DIM>::writeReadable(std::ostream& o) const
{
  o << m_particle_indices.size() << "\n";
  
  const VectorXs& x = m_parent->getX();
  for(int pidx : m_particle_indices)
  {
    o << x.segment<DIM>( m_parent->getDof(pidx) ).transpose() << " " << m_parent->getScriptedGroup(pidx) << "\n";
  }
}

template<int DIM>
void HairFlow<DIM>::readReadable(std::istream& i)
{
  int np = 0;
  i >> np;
  VectorXs& x = m_parent->getX();
  int group = -1;
  for(int pidx = 0; pidx < np; ++pidx )
  {
    Vector3s pos;
    i >> pos[0] >> pos[1] >> pos[2] >> group;
    x.segment<3>( m_parent->getDof( m_particle_indices[ pidx ] ) ) = pos;
  }
}


template<int DIM>
void HairFlow<DIM>::write(std::vector<scalar>& data) const
{
  int neta = m_eta.size();
  for(int i = 0; i < neta; ++i)
  {
    data.push_back(m_eta(i));
  }
  int nu = m_u.size();
  for(int i = 0; i < nu; ++i)
  {
    data.push_back(m_u(i));
  }
}

template<int DIM>
size_t HairFlow<DIM>::serialized_size() const
{
  return (m_eta.size() + m_u.size()) * sizeof(scalar);
}

template<int DIM>
Vectors<DIM> HairFlow<DIM>::computeHairLiquidMomentum(const VectorXs& v) const
{
  //[H] todo, should twist affect computeHairLiquidMomentum?
  int ne = m_global_edges.size();
  
  scalar rho = m_parent->getLiquidDensity();
  
  Vectors<DIM> momentum = Vectors<DIM>::Zero();
  
  for(int i = 0; i < ne; ++i)
  {
    int i0 = m_internal_edges[i].first;
    int i1 = m_internal_edges[i].second;
    int global_i0 = m_global_edges[i].first;
    int global_i1 = m_global_edges[i].second;
    const Vectors<DIM>& v0 = v.segment<DIM>( m_parent->getDof(global_i0) );
    const Vectors<DIM>& v1 = v.segment<DIM>( m_parent->getDof(global_i1) );
    
    const scalar H0 = m_eta(i0) + m_rad_vec(i0);
    const scalar H1 = m_eta(i1) + m_rad_vec(i1);
    
    scalar vol = M_PI * (H0 * H0 + H1 * H1 - m_rad_vec(i0) * m_rad_vec(i0) - m_rad_vec(i1) * m_rad_vec(i1)) * 0.5 * m_area_e(i);
    scalar mass = rho * vol;
    
    Vectors<DIM> ve = (v0 + v1) * 0.5 + m_dir_f.row(i).transpose().normalized() * m_u(i);
    
    momentum += mass * ve;
  }
  
  return momentum;
}

template<int DIM>
const MatrixXs& HairFlow<DIM>::getAccelV() const
{
  return m_accel_v;
}

template<int DIM>
MatrixXs& HairFlow<DIM>::getAccelV()
{
  return m_accel_v;
}


template<int DIM>
Vectors<DIM> HairFlow<DIM>::computeHairDragForce(const VectorXs& v) const
{
  //[H] todo, should twist affect computeHairDragForce?
  int ne = m_global_edges.size();
  
  scalar rho = m_parent->getLiquidDensity();
  
  Vectors<DIM> momentum = Vectors<DIM>::Zero();
  
  for(int i = 0; i < ne; ++i)
  {
    int i0 = m_internal_edges[i].first;
    int i1 = m_internal_edges[i].second;
    int global_i0 = m_global_edges[i].first;
    int global_i1 = m_global_edges[i].second;
    const Vectors<DIM>& v0 = v.segment<DIM>( m_parent->getDof(global_i0) );
    const Vectors<DIM>& v1 = v.segment<DIM>( m_parent->getDof(global_i1) );
    
    const scalar H0 = m_eta(i0) + m_rad_vec(i0);
    const scalar H1 = m_eta(i1) + m_rad_vec(i1);
    
    scalar vol = M_PI * (H0 * H0 + H1 * H1 - m_rad_vec(i0) * m_rad_vec(i0) - m_rad_vec(i1) * m_rad_vec(i1)) * 0.5 * m_area_e(i);
    scalar mass = rho * vol;
    
    Vectors<DIM> ve = (v0 + v1) * 0.5;
    
    momentum += mass * ve;
  }
  
  return momentum;
}


template<int DIM>
scalar HairFlow<DIM>::computeHairLiquidEnergy(const VectorXs& v) const
{
  //[H] todo, should twist affect computeHairLiquidEnergy?
  int ne = m_global_edges.size();
  
  scalar rho = m_parent->getLiquidDensity();
  
  scalar E = 0;
  
  for(int i = 0; i < ne; ++i)
  {
    int i0 = m_internal_edges[i].first;
    int i1 = m_internal_edges[i].second;
    int global_i0 = m_global_edges[i].first;
    int global_i1 = m_global_edges[i].second;
    const Vectors<DIM>& v0 = v.segment<DIM>( m_parent->getDof(global_i0) );
    const Vectors<DIM>& v1 = v.segment<DIM>( m_parent->getDof(global_i1) );
    
    const scalar H0 = m_eta(i0) + m_rad_vec(i0);
    const scalar H1 = m_eta(i1) + m_rad_vec(i1);
    
    scalar vol = M_PI * (H0 * H0 + H1 * H1 - m_rad_vec(i0) * m_rad_vec(i0) - m_rad_vec(i1) * m_rad_vec(i1)) * 0.5 * m_area_e(i);
    scalar mass = rho * vol;
    
    Vectors<DIM> ve = (v0 + v1) * 0.5 + m_dir_f.row(i).transpose().normalized() * m_u(i);
    
    E += mass * ve.squaredNorm();
  }
  
  return E * 0.5;
}

template<int DIM>
void HairFlow<DIM>::setConstraintParameters(const Vector6i& start, const Vector6i& num)
{
  m_constraint_starts = start;
  m_constraint_length = num;
}

template<int DIM>
const Vector6i& HairFlow<DIM>::getConstraintIdx() const
{
  return m_constraint_starts;
}

template<int DIM>
const Vector6i& HairFlow<DIM>::getNumConstraints() const
{
  return m_constraint_length;
}

// explicit instantiations at bottom
template class HairFlow<2>;
template class HairFlow<3>;

template struct HairParticleBridge<2>;
template struct HairParticleBridge<3>;
