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

#include "PolygonalCohesion.h"
#include "TwoDScene.h"
#include "MathUtilities.h"
#include "HairFlow.h"
#include "CohesionTableGen.h"
#include "fluidsim2D.h"
#include "fluidsim3D.h"
#include "ThreadUtils.h"
#include "CTCD.h"

#include <stack>
#include <deque>
#include <nanoflann/nanoflann.hpp>
#include <set>
#include <numeric>
#include <algorithm>

//#define DEBUG_COHESION_TABLE
//#define DEBUG_PLANAR_COHESION_TABLE

template<int DIM>
PolygonalCohesion<DIM>::PolygonalCohesion(TwoDScene<DIM>* scene) :
m_parent(scene), m_tree(NULL), m_sorter(NULL), m_use_decoupled_force(false), m_compute_particle_poe_mapping(true), m_min_cohesion_table(NULL), m_max_cohesion_table(NULL)
{
  m_sorter = new Sorter(0, 0, 0);
  
  scalar gnorm = scene->getSimpleGravity().norm();
  gnorm = (gnorm == 0.0) ? 1000.0 : gnorm;
  const scalar max_eta = sqrt(scene->getLiquidTension() / (scene->getLiquidDensity() * gnorm)) * 0.5;
  const scalar sigma = scene->getLiquidTension();
  const scalar theta = scene->getLiquidTheta();

  const int disc = 256;
  
  if(m_parent->getNumParticles() > 0) {
    const WetHairParameter& parameter = scene->getParameter();
    
    const scalar max_rad = scene->getRadii().maxCoeff();
    const scalar min_rad = scene->getRadii().minCoeff();
    
    m_min_cohesion_table = new CohesionTable( parameter.radius_multiplier, parameter.collision_stiffness, parameter.radius_multiplier_planar, parameter.collision_stiffness_planar );
    
    m_max_cohesion_table = new CohesionTable( parameter.radius_multiplier, parameter.collision_stiffness, parameter.radius_multiplier_planar, parameter.collision_stiffness_planar );
    
    m_min_cohesion_table->setParameter(sigma, theta, min_rad, max_eta * 2.0, disc);
    
    m_max_cohesion_table->setParameter(sigma, theta, max_rad, max_eta * 2.0, disc);
    
    std::cout << "[construct adhesive/repulsive table]" << std::endl;
    m_min_cohesion_table->construct_alpha_table();
    m_max_cohesion_table->construct_alpha_table();
    
    m_min_cohesion_table->construct_planar_alpha_table();
    m_max_cohesion_table->construct_planar_alpha_table();
  }
  
#ifdef DEBUG_COHESION_TABLE
  scalar max_vol = M_PI * ((max_eta + avg_rad) * (max_eta + avg_rad) - avg_rad * avg_rad);
  std::cout << "[sweep adhesive/repulsive table]" << std::endl;
  
  MatrixXs stiffness_table(disc, disc);
  MatrixXs force_table(disc, disc);
  MatrixXs dEdd_table(disc, disc);
  const scalar d_star = m_cohesion_table.getDStar();
  for(int i = 0; i < disc; ++i)
  {
    for(int j = 0; j < disc; ++j)
    {
      scalar d = max_dist / (scalar) disc * (scalar) i;
      scalar A_L = max_vol / (scalar) disc * (scalar) j;
      scalar stiffness = m_cohesion_table.getStiffness(d, A_L, 1.0);
      scalar dEdd = m_cohesion_table.interpolate_dEdd(A_L, d);
      dEdd_table(j, i) = dEdd;
      stiffness_table(j, i) = stiffness;
      force_table(j, i) = stiffness * (d - d_star);
    }
  }
  
  std::cout << "[adhesive table]" << std::endl;
  std::cout << dEdd_table << std::endl;
  
  std::cout << "[stiffness table]" << std::endl;
  std::cout << stiffness_table << std::endl;
  
  std::cout << "[force table]" << std::endl;
  std::cout << force_table << std::endl;
#endif
  
#ifdef DEBUG_PLANAR_COHESION_TABLE
  scalar max_vol = M_PI * ((max_eta + avg_rad) * (max_eta + avg_rad) - avg_rad * avg_rad);
  std::cout << "[sweep planar adhesive/repulsive table]" << std::endl;
  
  MatrixXs stiffness_table(disc, disc);
  MatrixXs force_table(disc, disc);
  MatrixXs dEdd_table(disc, disc);
  const scalar d_star = m_cohesion_table.getDStarPlanar();
  for(int i = 0; i < disc; ++i)
  {
    for(int j = 0; j < disc; ++j)
    {
      scalar d = max_dist / (scalar) disc * (scalar) i;
      scalar A_L = max_vol / (scalar) disc * (scalar) j;
      scalar stiffness = m_cohesion_table.getStiffnessPlanar(d, A_L, 1.0);
      scalar dEdd = m_cohesion_table.interpolate_dEdd_planar(A_L, d);
      dEdd_table(j, i) = dEdd;
      stiffness_table(j, i) = stiffness;
      force_table(j, i) = stiffness * (d - d_star);
    }
  }
  
  std::cout << "[planar adhesive table]" << std::endl;
  std::cout << dEdd_table << std::endl;
  
  std::cout << "[planar stiffness table]" << std::endl;
  std::cout << stiffness_table << std::endl;
  
  std::cout << "[planar force table]" << std::endl;
  std::cout << force_table << std::endl;
#endif
}

template<int DIM>
PolygonalCohesion<DIM>::~PolygonalCohesion()
{
  if(m_sorter) delete m_sorter;
  if(m_max_cohesion_table) delete m_max_cohesion_table;
  if(m_min_cohesion_table) delete m_min_cohesion_table;
  
  for(auto& adj : m_adjacency_categorized)
  {
    for(auto& hp : adj) {
      if(hp.second) delete hp.second;
    }
  }
}

template<int DIM>
scalar PolygonalCohesion<DIM>::computeEdgeThickness(const scalar& vol_p, const int_Vectors_scalar<DIM>& ep)
{
  const scalar d0 = ep.v.norm();
  
  const scalar A = vol_p + M_PI * ep.eta * ep.eta;
  
  const scalar alpha_0 = m_min_cohesion_table->interpolate_alpha(A, d0);
  
  const scalar R_0 = m_min_cohesion_table->computeR(alpha_0, d0);
  
  const scalar H_0 = m_min_cohesion_table->computeH(R_0, alpha_0);
  
  const scalar alpha_1 = m_max_cohesion_table->interpolate_alpha(A, d0);
  
  const scalar R_1 = m_max_cohesion_table->computeR(alpha_1, d0);
  
  const scalar H_1 = m_max_cohesion_table->computeH(R_1, alpha_1);
  
  return (H_0 + H_1) * 0.5;
}

template<int DIM>
void PolygonalCohesion<DIM>::findEdgeEdgeContact( const VectorXs x, const VectorXs v, const scalar& dt, const int& base_eidx )
{
  if( DIM == 2 ) return;
  const std::vector<int>& particle_to_hairs = m_parent->getParticleToHairs();
  const std::vector<int>& particle_local_indices = m_parent->getParticleToHairLocalIndices();
  const std::vector< std::pair<int, int> >& scene_edges = m_parent->getEdges();
  
  const std::pair<int, int>& base_edge = scene_edges[ base_eidx ];
  
  const VectorXs& radii = m_parent->getRadii();
  const scalar& radius_0 = radii(base_edge.first);
  const scalar& radius_1 = radii(base_edge.second);
  
  auto& flows = m_parent->getFilmFlows();
  const scalar& theta = m_parent->getLiquidTheta();
  
  int pidx = base_edge.first;
  int hidx = particle_to_hairs[pidx];
  if(hidx < 0) return;
  
  int local_idx_0 = particle_local_indices[base_edge.first];
  int local_idx_1 = particle_local_indices[base_edge.second];
  const VectorXs& eta = flows[hidx]->getEta();
  const scalar peta_0 = eta(local_idx_0);
  const scalar peta_1 = eta(local_idx_1);
  const scalar H_0 = peta_0 + radius_0;
  const scalar H_1 = peta_1 + radius_1;
  const int pepair_neighbor_0 = std::max(1, (int) m_num_edge_connections[base_edge.first] );
  const int pepair_neighbor_1 = std::max(1, (int) m_num_edge_connections[base_edge.second] );
  
  const scalar vol_p_each = M_PI * ((H_0 * H_0 - radius_0 * radius_0) / (scalar) pepair_neighbor_0 + (H_1 * H_1 - radius_1 * radius_1) / (scalar) pepair_neighbor_1) * 0.5;
  const scalar radii_p = sqrt((radius_0 * radius_0 + radius_1 * radius_1) * 0.5);
  
  // Edge AB
  // end positions of base edge
  Vector3s xAe = x.segment<3>( m_parent->getDof( base_edge.first ) );
  Vector3s xBe = x.segment<3>( m_parent->getDof( base_edge.second ) );
  // start positions of edge
  Vector3s xAs = xAe - v.segment<3>( m_parent->getDof( base_edge.first ) ) * dt;
  Vector3s xBs = xBe - v.segment<3>( m_parent->getDof( base_edge.second) ) * dt;
  
  // grab center of Edge
  const Vector3s& pos = ( xAs + xBs ) * 0.5;
  const scalar& cellsize = m_parent->getSearchRadius();
  const Vectors<DIM>& bbx_min = m_parent->getBoundingBoxMin();
  int pindex[3] = {0, 0, 0};
  int range[3] = {0, 0, 0};
  const int nindex[3] = {m_sorter->ni, m_sorter->nj, m_sorter->nk};
  for(int r = 0; r < DIM; ++r)
  {
    pindex[r] = std::max(0, std::min(nindex[r] - 1, (int)((pos(r) - bbx_min(r)) / cellsize)));
    range[r] = 1;
  }
  
  auto& edge_pairs = m_edge_connections[base_eidx];
  for(auto it = std::begin(edge_pairs); it != std::end(edge_pairs); ++it)
  {
    it->second->updated = false;
  }
  
  const scalar H_each = sqrt(vol_p_each / M_PI + radii_p * radii_p);
  
  // get closest edges to each hair
  m_sorter->getNeigboringParticles_cell(pindex[0], pindex[1], pindex[2], -range[0], range[0], -range[1], range[1], -range[2], range[2], [&] (int eidx) {
    // ignore myself
    const std::pair<int, int>& e = scene_edges[eidx];
    if( e == base_edge ) return;
    
    // ignore edges on my strand
    int nhidx = particle_to_hairs[e.first];
    if(nhidx < 0 || nhidx == hidx) return;
    
    //if( eidx < base_eidx ) return; // avoid duplicates
    
    // Grab edge CD
    // end of timestep positions
    Vector3s xCe = x.segment<3>( m_parent->getDof( e.first ) );
    Vector3s xDe = x.segment<3>( m_parent->getDof( e.second ) );
    // start of timestep positions
    Vector3s xCs = xCe - v.segment<3>( m_parent->getDof( e.first ) ) * dt;
    Vector3s xDs = xDe - v.segment<3>( m_parent->getDof( e.second) ) * dt;
    
    // averaging radii_q since we don't have an alpha yet...
    const scalar radii_q = sqrt(mathutils::lerp(radii(e.first) * radii(e.first), radii(e.second) * radii(e.second), 0.5));
    HairFlow<DIM>* neighbor_flow = flows[nhidx];
    const VectorXs& neta = neighbor_flow->getEta();
    
    int nlidx_0 = particle_local_indices[e.first];
    int nlidx_1 = particle_local_indices[e.second];
    const int npepair_neighbor_0 = std::max(1, (int) m_num_edge_connections[e.first] );
    const int npepair_neighbor_1 = std::max(1, (int) m_num_edge_connections[e.second] );
    
    const scalar npeta_0 = neta(nlidx_0);
    const scalar npeta_1 = neta(nlidx_1);
    const scalar nH_0 = npeta_0 + radii(e.first);
    const scalar nH_1 = npeta_1 + radii(e.second);
    const scalar vol_q_each = M_PI * ((nH_0 * nH_0 - radii(e.first) * radii(e.first)) / (scalar) npepair_neighbor_0 + (nH_1 * nH_1 - radii(e.second) * radii(e.second)) / (scalar) npepair_neighbor_1) * 0.5;
    
    scalar nH_each = sqrt(vol_q_each / M_PI + radii_q * radii_q);
    
    auto itr = edge_pairs.find( eidx );
    EdgeEdgePairEEC* ref_pair = NULL;
    bool previous_linked = false;
    if( itr != edge_pairs.end() ){ // previously found
      ref_pair = itr->second;
      previous_linked = ref_pair->valid;
    }
    
    scalar max_dist;
    if(previous_linked) {
      max_dist = std::max(H_each + nH_each, (1.0 + 0.5 * theta) * sqrt(vol_p_each + vol_q_each) + radii_p + radii_q);
    } else {
      max_dist = std::min(H_each + nH_each, (1.0 + 0.5 * theta) * sqrt(vol_p_each + vol_q_each) + radii_p + radii_q);
    }
    
    double colTime = -1.0;
    if( CTCD::checkEEContact( xAs, xBs, xCs, xDs, xAe, xBe, xCe, xDe, max_dist, colTime ) ){
      if( !ref_pair ){
        ref_pair = new EdgeEdgePairEEC();
        
        edge_pairs[ eidx ] = ref_pair;
        m_num_edge_connections[ base_edge.first ] += 1;
        m_num_edge_connections[ base_edge.second ] += 1;
        m_num_edge_connections[ e.first ] += 1;
        m_num_edge_connections[ e.second ] += 1;
      }
      
      ref_pair->valid = true;
      ref_pair->base_eidx = base_eidx;
      ref_pair->neighbor_eidx = eidx;
      
      // Time of collision
      const Vector3s xAc = (1.0 - colTime) * xAs + colTime * xAe;
      const Vector3s xBc = (1.0 - colTime) * xBs + colTime * xBe;
      const Vector3s xCc = (1.0 - colTime) * xCs + colTime * xCe;
      const Vector3s xDc = (1.0 - colTime) * xDs + colTime * xDe;
      
      // points of contact on edges
      Vector3s xAB, xCD;
      scalar alpha, beta;
      double dist_squared = CTCD::ClosestPtSegmentSegment( xAc, xBc, xCc, xDc, alpha, beta, xAB, xCD );
      double d_contact = sqrt( dist_squared );
      
      EdgeEdgePairEEC* eep = ref_pair;
      eep->time = colTime;
      eep->alpha_contact = alpha;
      eep->neighbor_local_coord_contact = beta;
      
      // compute quadruture region from contact points
      // Edge d1 quad_region,
      scalar AC = (xAc - xCc).norm();
      scalar AD = (xAc - xDc).norm();
      scalar max_alpha = m_parent->isTip( base_edge.first ) ? 0.9 : 1.0;
      if( AC < AD ){
        if( AC == d_contact ){
          eep->alpha_0 = 0.0;
          eep->neighbor_local_coord_0 = 0.0;
        } else {
          eep->alpha_0 = mathutils::clamp( alpha * ( ( max_dist - AC) / (d_contact - AC) ) , 0.0, max_alpha );
          max_alpha = m_parent->isTip( e.first ) ? 0.9 : 1.0;
          eep->neighbor_local_coord_0 = mathutils::clamp( beta * ( ( max_dist - AC) / (d_contact - AC) ) , 0.0, max_alpha );
        }
      }
      else{
        if( AD == d_contact ){
          eep->alpha_0 = 0.0;
          eep->neighbor_local_coord_1 = 0.0;
        } else {
          eep->alpha_0 = mathutils::clamp( alpha * ( ( max_dist - AD) / (d_contact - AD) ) , 0.0, max_alpha );
          max_alpha = m_parent->isTip( e.second ) ? 0.9 : 1.0;
          eep->neighbor_local_coord_1 = mathutils::clamp( beta * ( ( max_dist - AD) / (d_contact - AD) ) , 0.0, max_alpha );
        }
      }
      
      // Edge d2 quad_region,
      scalar BC = (xBc - xCc).norm();
      scalar BD = (xBc - xDc).norm();
      max_alpha = m_parent->isTip( base_edge.second ) ? 0.9 : 1.0;
      if( BC < BD ){
        if( BC == d_contact ){
          eep->alpha_1 = 0.0;
          eep->neighbor_local_coord_0 = 0.0;
        } else {
          eep->alpha_1 = mathutils::clamp( alpha * ( ( max_dist - BC) / (d_contact - BC) ) , 0.0, max_alpha );
          max_alpha = m_parent->isTip( e.first ) ? 0.9 : 1.0;
          eep->neighbor_local_coord_0 = mathutils::clamp( beta * ( ( max_dist - BC) / (d_contact - BC) ) , 0.0, max_alpha );
        }
      }
      else{
        if( BD == d_contact ){
          eep->alpha_1 = 0.0;
          eep->neighbor_local_coord_1 = 0.0;
        } else {
          eep->alpha_1 = mathutils::clamp( alpha * ( ( max_dist - BD) / (d_contact - BD) ) , 0.0, max_alpha );
          max_alpha = m_parent->isTip( e.second ) ? 0.9 : 1.0;
          eep->neighbor_local_coord_1 = mathutils::clamp( beta * ( ( max_dist - BD) / (d_contact - BD) ) , 0.0, max_alpha );
        }
      }
      
      eep->avgpos = ( xAc + xBc + xCc + xDc ) * 0.25;
    } // CCD check
    else if( ref_pair ){
      ref_pair->valid = false;
    }
    
    if(ref_pair) ref_pair->updated = true;
  }); // edge neighbors
  
  //remove EdgeEdgePairs that are no longer valid/updated
  for(auto it = std::begin(edge_pairs); it != std::end(edge_pairs);)
  {
    if( !it->second->valid || !it->second->updated ){
      delete it->second;
      
      const std::pair<int, int>& e = scene_edges[it->first];
      m_num_edge_connections[ e.first ] = std::max( 0, m_num_edge_connections[ e.first ] - 1 );
      m_num_edge_connections[ e.second ] = std::max( 0, m_num_edge_connections[ e.second ] - 1 );
      m_num_edge_connections[ base_edge.first ] = std::max( 0, m_num_edge_connections[ base_edge.first ] - 1 );
      m_num_edge_connections[ base_edge.second ] = std::max( 0, m_num_edge_connections[ base_edge.second ] - 1 );
      
      it = edge_pairs.erase(it);
    }
    else ++it;
  }
}

template<int DIM>
void PolygonalCohesion<DIM>::findParticleEdgePairs(const VectorXs& x, int pidx)
{
  const std::vector<int>& particle_to_hairs = m_parent->getParticleToHairs();
  const std::vector<int>& particle_local_indices = m_parent->getParticleToHairLocalIndices();
  const std::vector< std::pair<int, int> >& scene_edges = m_parent->getEdges();
  const VectorXs& radii = m_parent->getRadii();
  const scalar& radius_i = radii(pidx);
  auto& flows = m_parent->getFilmFlows();
  const scalar& theta = m_parent->getLiquidTheta();
  const std::vector< std::vector<int> >& particle_edges = m_parent->getParticleToEdge();
  
  int hidx = particle_to_hairs[pidx];
  if(hidx < 0) return;
  
  int local_idx = particle_local_indices[pidx];
  
  const VectorXs& eta = flows[hidx]->getEta();
  
  const scalar peta = eta(local_idx);
  
  const scalar H_p = peta + radius_i;
  
  const scalar vol_p = M_PI * (H_p * H_p - radius_i * radius_i);
  
  const Vectors<DIM>& pos = x.segment<DIM>( m_parent->getDof( pidx ) );
  
  const scalar& cellsize = m_parent->getSearchRadius();
  const Vectors<DIM>& bbx_min = m_parent->getBoundingBoxMin();
  
  int pindex[3] = {0, 0, 0};
  int range[3] = {0, 0, 0};
  const int nindex[3] = {m_sorter->ni, m_sorter->nj, m_sorter->nk};
  for(int r = 0; r < DIM; ++r)
  {
    pindex[r] = std::max(0, std::min(nindex[r] - 1, (int)((pos(r) - bbx_min(r)) / cellsize)));
    range[r] = 1;
  }
  
  auto& hair_pairs = m_adjacency_categorized[pidx];
  const int npepair_local = std::max(1, m_num_adjacency_categorized[pidx]);
  const scalar vol_p_each = vol_p / (scalar) npepair_local;
  const scalar H_each = sqrt(vol_p / M_PI + radius_i * radius_i);
  
  for(auto it = std::begin(hair_pairs); it != std::end(hair_pairs);)
  {
    // delete all original pairs that is marked as should-be-deleted and not marked as a part of valid double-linked-edge
    // see the computation of point-edge pair for details
    if(!it->second->updated && it->second->should_be_deleted) {
      delete it->second;
      it = hair_pairs.erase(it);
    }
    else ++it;
  }
  
  // mark all the original pairs as not-updated & original
  for(auto hp : hair_pairs)
  {
    ParticleEdgePair* p = hp.second;
    p->updated = false;
    p->latest = false;
  }
  
  // get closest edge to each hair
  m_sorter->getNeigboringParticles_cell(pindex[0], pindex[1], pindex[2], -range[0], range[0], -range[1], range[1], -range[2], range[2], [&] (int eidx) {
    const std::pair<int, int>& e = scene_edges[eidx];
    if(e.first == pidx || e.second == pidx)
      return;
    
    int nhidx = particle_to_hairs[e.first];
    
    if(nhidx < 0 || nhidx == hidx) return;
    
    Vectors<DIM> gap_vec;
    scalar alpha;
    
    mathutils::pointedgevec<scalar, DIM>(pos, x.segment<DIM>( m_parent->getDof( e.first ) ), x.segment<DIM>( m_parent->getDof( e.second ) ), gap_vec, alpha);
    scalar dist = gap_vec.norm();
    
    auto itr = hair_pairs.find(nhidx);
    
    ParticleEdgePair* ref_pair = NULL;
    
    bool previous_linked = false;
    
    if(itr != hair_pairs.end()) {
      ref_pair = itr->second;
      if(ref_pair->updated && dist > ref_pair->dist) return;
      if(!ref_pair->latest) previous_linked = true;
    }
    
    const scalar radius_j = sqrt(mathutils::lerp(radii(e.first) * radii(e.first), radii(e.second) * radii(e.second), alpha));
    
    if(!ref_pair) {
      ref_pair = new ParticleEdgePair;
      ref_pair->latest = true;
      hair_pairs[nhidx] = ref_pair;
    }
    
    ref_pair->pidx = pidx;
    ref_pair->eidx = eidx;
    ref_pair->alpha = alpha;
    ref_pair->dist = dist;
    ref_pair->radii_j = radius_j;
    ref_pair->updated = true;
    
    HairFlow<DIM>* neighbor_flow = flows[nhidx];
    const VectorXs& neta = neighbor_flow->getEta();
    
    int nlidx_0 = particle_local_indices[e.first];
    int nlidx_1 = particle_local_indices[e.second];
    
    const int npepair_neighbor_0 = std::max(1, (int) m_num_adjacency_categorized[e.first]);
    const int npepair_neighbor_1 = std::max(1, (int) m_num_adjacency_categorized[e.second]);
    
    const scalar npeta_0 = neta(nlidx_0);
    const scalar npeta_1 = neta(nlidx_1);
    const scalar nH_0 = npeta_0 + radii(e.first);
    const scalar nH_1 = npeta_1 + radii(e.second);
    const scalar vol_q_each = M_PI * mathutils::lerp((nH_0 * nH_0 - radii(e.first) * radii(e.first)) / (scalar) npepair_neighbor_0
                                                     , (nH_1 * nH_1 - radii(e.second) * radii(e.second)) / (scalar) npepair_neighbor_1, alpha);
    
    scalar nH_each = sqrt(vol_q_each / M_PI + radius_j * radius_j);
    
    scalar max_dist;
    
    if(previous_linked) {
      max_dist = std::max(H_each + nH_each, (1.0 + 0.5 * theta) * sqrt(vol_p_each + vol_q_each) + radius_i + radius_j);
    } else {
      max_dist = std::min(H_each + nH_each, (1.0 + 0.5 * theta) * sqrt(vol_p_each + vol_q_each) + radius_i + radius_j);
    }
    
    ref_pair->valid = (dist <= max_dist);
    ref_pair->max_dist = max_dist;
  });
  
  // update all pairs that are not updated
  for(auto hp : hair_pairs)
  {
    ParticleEdgePair* pepair = hp.second;
    
    if(pepair->updated) continue;
    
    bool previous_linked = pepair->valid;
    
    const std::pair<int, int>& e = scene_edges[pepair->eidx];
    
    Vectors<DIM> gap_vec;
    scalar alpha;
    mathutils::pointedgevec<scalar, DIM>(pos, x.segment<DIM>( m_parent->getDof( e.first ) ), x.segment<DIM>( m_parent->getDof( e.second ) ), gap_vec, alpha);
    scalar dist = gap_vec.norm();
    
    pepair->alpha = alpha;
    pepair->dist = dist;
    
    int nhidx = particle_to_hairs[e.first];
    
    HairFlow<DIM>* neighbor_flow = flows[nhidx];
    const VectorXs& neta = neighbor_flow->getEta();
    
    int nlidx_0 = particle_local_indices[e.first];
    int nlidx_1 = particle_local_indices[e.second];
    
    const int npepair_neighbor_0 = std::max(1, (int) m_num_adjacency_categorized[e.first]);
    const int npepair_neighbor_1 = std::max(1, (int) m_num_adjacency_categorized[e.second]);
    
    const scalar npeta_0 = neta(nlidx_0);
    const scalar npeta_1 = neta(nlidx_1);
    const scalar nH_0 = npeta_0 + radii(e.first);
    const scalar nH_1 = npeta_1 + radii(e.second);
    const scalar vol_q_each = M_PI * mathutils::lerp((nH_0 * nH_0 - radii(e.first) * radii(e.first)) / (scalar) npepair_neighbor_0
                                                     , (nH_1 * nH_1 - radii(e.second) * radii(e.second)) / (scalar) npepair_neighbor_1, alpha);
    
    scalar nH_each = sqrt(vol_q_each / M_PI + pepair->radii_j * pepair->radii_j);
    
    scalar max_dist;
    
    if(previous_linked) {
      max_dist = std::max(H_each + nH_each, (1.0 + 0.5 * theta) * sqrt(vol_p_each + vol_q_each) + radius_i + pepair->radii_j);
    } else {
      max_dist = std::min(H_each + nH_each, (1.0 + 0.5 * theta) * sqrt(vol_p_each + vol_q_each) + radius_i + pepair->radii_j);
    }
    
    pepair->valid = (dist <= max_dist);
    pepair->max_dist = max_dist;
  }
  
  // mark all the as not-updated again for later use, also all as not to be deleted
  for(auto hp : hair_pairs)
  {
    ParticleEdgePair* p = hp.second;
    p->updated = false;
    p->should_be_deleted = false;
  }
}

template<int DIM>
void PolygonalCohesion<DIM>::findParticleParticlePairsEEC(const VectorXs& x)
{
  auto& edges = m_parent->getEdges();
  auto& radii = m_parent->getRadii();
  const int np = m_parent->getNumParticles();
  const int ne = m_parent->getNumEdges();
  if(!np) return;
  if(!ne) return;
  
  if(m_particle_to_pppairs.size() != np) m_particle_to_pppairs.resize(np);
  
  if(m_num_valid_edge_connections.size() != ne) m_num_valid_edge_connections.resize(ne);
  
  // foreach edge count valid eep
  threadutils::thread_pool::ParallelFor(0, ne, [&] (int eidx){
    int count_valid = 0;
    auto& adjs = m_edge_connections[eidx];
    for(auto itr = adjs.begin(); itr != adjs.end(); ++itr)
    {
      const EdgeEdgePairEEC* eep = itr->second;
      if(!eep->valid) continue;
      
      ++count_valid;
    }
    
    m_num_valid_edge_connections[eidx] = count_valid;
  });
  
  m_counting_valid_adjacency = m_num_valid_edge_connections;
  
  // prefix sum to get the locations of pp-pair
  std::partial_sum(m_counting_valid_adjacency.begin(), m_counting_valid_adjacency.end(), m_counting_valid_adjacency.begin());
  
  m_counting_pp_pairs.resize(m_counting_valid_adjacency[m_counting_valid_adjacency.size() - 1] * 2);
  if(m_counting_pp_pairs.size() == 0) {
    m_particle_particle_pairs.resize(0);
    return;
  }
  
  // foreach edge, store the hash into locations
  threadutils::thread_pool::ParallelFor(0, ne, [&] (int base_eidx){
    const std::pair<int, int>& base_e = edges[base_eidx];
    auto& adjs = m_edge_connections[base_eidx];
    int base_loc = (base_eidx == 0) ? 0 : m_counting_valid_adjacency[base_eidx - 1];
    
    int k = 0;
    for(auto itr = adjs.begin(); itr != adjs.end(); ++itr)
    {
      const EdgeEdgePairEEC* eep = itr->second;
      if(!eep->valid) continue;
      
      const int loc = base_loc + k;
      
      auto& neighbor_e = edges[itr->first];
      
      // check distance
      scalar d00 = (x.segment<DIM>( m_parent->getDof(base_e.first) ) - x.segment<DIM>( m_parent->getDof(neighbor_e.first) )).norm();
      scalar d01 = (x.segment<DIM>( m_parent->getDof(base_e.first) ) - x.segment<DIM>( m_parent->getDof(neighbor_e.second) )).norm();
      
      scalar d10 = (x.segment<DIM>( m_parent->getDof(base_e.second) ) - x.segment<DIM>( m_parent->getDof(neighbor_e.first) )).norm();
      scalar d11 = (x.segment<DIM>( m_parent->getDof(base_e.second) ) - x.segment<DIM>( m_parent->getDof(neighbor_e.second) )).norm();
      
      if(d00 < d01) {
        m_counting_pp_pairs[loc * 2 + 0] = makePairwisePPHash(base_e.first, neighbor_e.first);
      } else {
        m_counting_pp_pairs[loc * 2 + 0] = makePairwisePPHash(base_e.first, neighbor_e.second);
      }
      
      if(d10 < d11) {
        m_counting_pp_pairs[loc * 2 + 1] = makePairwisePPHash(base_e.second, neighbor_e.first);
      } else {
        m_counting_pp_pairs[loc * 2 + 1] = makePairwisePPHash(base_e.second, neighbor_e.second);
      }
      ++k;
    }
  });
  
  // sort the hash
  tbb::parallel_sort(m_counting_pp_pairs.begin(), m_counting_pp_pairs.end());
  
  const int n_pred_pairs = m_counting_pp_pairs.size();
  m_counting_pp_pair_location.resize(n_pred_pairs);
  
  // marking for the one where hash is different than previous
  threadutils::thread_pool::ParallelFor(0, n_pred_pairs, [&] (int i) {
    m_counting_pp_pair_location[i] = (i != 0 && m_counting_pp_pairs[i] != m_counting_pp_pairs[i - 1]);
  });
  
  // prefix sum to get relocated locations
  std::partial_sum(m_counting_pp_pair_location.begin(), m_counting_pp_pair_location.end(), m_counting_pp_pair_location.begin());
  
  // actually record the PP-pairs
  m_particle_particle_pairs.resize(m_counting_pp_pair_location[n_pred_pairs - 1] + 1);
  
  threadutils::thread_pool::ParallelFor(0, n_pred_pairs, [&] (int i) {
    if(i == 0 || m_counting_pp_pairs[i] != m_counting_pp_pairs[i - 1]) {
      int storing_loc = m_counting_pp_pair_location[i];
      uint64_t hash_code = m_counting_pp_pairs[i];
      ParticleParticlePair& ppp = m_particle_particle_pairs[storing_loc];
      ppp.pidx[0] = (int)((hash_code >> 32UL) & 0xFFFFFFFFUL);
      ppp.pidx[1] = (int)(hash_code & 0xFFFFFFFFUL);
      scalar r0 = radii(ppp.pidx[0]);
      scalar r1 = radii(ppp.pidx[1]);
      
      ppp.d = std::max(r0 + r1, (x.segment<DIM>( m_parent->getDof(ppp.pidx[0]) ) - x.segment<DIM>( m_parent->getDof(ppp.pidx[1]) )).norm());
      ppp.r = sqrt(mathutils::lerp(r0 * r0, r1 * r1, 0.5));
    }
  });
  
  // accumulate pidx -> other particles
  for(int i = 0; i < np; ++i)
  {
    m_particle_to_pppairs[i].resize(0);
  }
  const int npppairs = (int) m_particle_particle_pairs.size();
  for(int i = 0; i < npppairs; ++i)
  {
    ParticleParticlePair& ppp = m_particle_particle_pairs[i];
    m_particle_to_pppairs[ppp.pidx[0]].push_back(ppp.pidx[1]);
    m_particle_to_pppairs[ppp.pidx[1]].push_back(ppp.pidx[0]);
  }
}

template<int DIM>
void PolygonalCohesion<DIM>::findParticleParticlePairs(const VectorXs& x)
{
  auto& edges = m_parent->getEdges();
  auto& radii = m_parent->getRadii();
  const int np = m_parent->getNumParticles();
  auto& global_to_local = m_parent->getParticleToHairLocalIndices();
  auto& particle_hair = m_parent->getParticleToHairs();
  auto& flows = m_parent->getFilmFlows();
  
  if(!np) return;
  
  m_counting_valid_adjacency.resize(np);
  
  // for each particle, record num of valid PE-pairs
  threadutils::thread_pool::ParallelFor(0, np, [&] (int pidx) {
    auto& adjs = m_adjacency_categorized[pidx];
    int num_valid = 0;
    auto& flow = flows[particle_hair[pidx]];
    int local_idx = global_to_local[pidx];
    if(local_idx != 0 && local_idx != flow->getParticleIndices().size() - 1) {
      
      for(auto itr = adjs.begin(); itr != adjs.end(); ++itr)
      {
        if(itr->second->valid && itr->second->dist >= radii(pidx) * 4.0) ++num_valid;
      }
    }
    m_counting_valid_adjacency[pidx] = num_valid;
  });
  
  // prefix sum to get the locations of pp-pair
  std::partial_sum(m_counting_valid_adjacency.begin(), m_counting_valid_adjacency.end(), m_counting_valid_adjacency.begin());
  
  m_counting_pp_pairs.resize(m_counting_valid_adjacency[m_counting_valid_adjacency.size() - 1]);
  
  if(m_counting_pp_pairs.size() == 0) {
    m_particle_particle_pairs.resize(0);
    return;
  }
  
  // store the hash into locations
  threadutils::thread_pool::ParallelFor(0, np, [&] (int pidx) {
    auto& flow = flows[particle_hair[pidx]];
    int local_idx = global_to_local[pidx];
    if(local_idx != 0 && local_idx != flow->getParticleIndices().size() - 1) {
      int base_loc = (pidx == 0) ? 0 : m_counting_valid_adjacency[pidx - 1];
      int k = 0;
      auto& adjs = m_adjacency_categorized[pidx];
      for(auto itr = adjs.begin(); itr != adjs.end(); ++itr)
      {
        const ParticleEdgePair* pep = itr->second;
        if(!pep->valid || itr->second->dist < radii(pidx) * 4.0) continue;
        
        auto& e = edges[pep->eidx];
        scalar d0 = (x.segment<DIM>(m_parent->getDof(pep->pidx)) - x.segment<DIM>(m_parent->getDof(e.first))).norm();
        scalar d1 = (x.segment<DIM>(m_parent->getDof(pep->pidx)) - x.segment<DIM>(m_parent->getDof(e.second))).norm();
        
        if(d0 < d1) m_counting_pp_pairs[base_loc + k] = makePairwisePPHash(pep->pidx, e.first);
        else m_counting_pp_pairs[base_loc + k] = makePairwisePPHash(pep->pidx, e.second);
        ++k;
      }
    }
  });
  
  // sort the hash
  tbb::parallel_sort(m_counting_pp_pairs.begin(), m_counting_pp_pairs.end());
  
  const int n_pred_pairs = m_counting_pp_pairs.size();
  m_counting_pp_pair_location.resize(n_pred_pairs);
  
  // marking for the one where hash is different than previous
  threadutils::thread_pool::ParallelFor(0, n_pred_pairs, [&] (int i) {
    m_counting_pp_pair_location[i] = (i != 0 && m_counting_pp_pairs[i] != m_counting_pp_pairs[i - 1]);
  });
  
  // prefix sum to get relocated locations
  std::partial_sum(m_counting_pp_pair_location.begin(), m_counting_pp_pair_location.end(), m_counting_pp_pair_location.begin());
  
  // actually record the PP-pairs
  m_particle_particle_pairs.resize(m_counting_pp_pair_location[n_pred_pairs - 1] + 1);
  
  threadutils::thread_pool::ParallelFor(0, n_pred_pairs, [&] (int i) {
    if(i == 0 || m_counting_pp_pairs[i] != m_counting_pp_pairs[i - 1]) {
      int storing_loc = m_counting_pp_pair_location[i];
      uint64_t hash_code = m_counting_pp_pairs[i];
      ParticleParticlePair& ppp = m_particle_particle_pairs[storing_loc];
      ppp.pidx[0] = (int)((hash_code >> 32UL) & 0xFFFFFFFFUL);
      ppp.pidx[1] = (int)(hash_code & 0xFFFFFFFFUL);
      scalar r0 = radii(ppp.pidx[0]);
      scalar r1 = radii(ppp.pidx[1]);
      
      ppp.d = std::max(r0 + r1, (x.segment<DIM>( m_parent->getDof(ppp.pidx[0]) ) - x.segment<DIM>( m_parent->getDof(ppp.pidx[1]) )).norm());
      ppp.r = sqrt(mathutils::lerp(r0 * r0, r1 * r1, 0.5));
    }
  });
}

template<int DIM>
void PolygonalCohesion<DIM>::findPointEdgePairsEEC(const VectorXs& x, int base_eidx, std::vector<PointEdgePair>& poepairs)
{
  if( DIM == 2 ) return;
  
  auto& edges = m_parent->getEdges();
  auto& e = edges[base_eidx];
  auto& global_to_local = m_parent->getParticleToHairLocalIndices();
  auto& particle_hair = m_parent->getParticleToHairs();
  auto& flows = m_parent->getFilmFlows();
  
  for( auto& adjPair : m_edge_connections[base_eidx] )
  {
    const EdgeEdgePairEEC& eepair = *(adjPair.second);
    
    int num_linked_hairs_e0 = std::max(1, (int) m_num_edge_connections[e.first] );
    int num_linked_hairs_e1 = std::max(1, (int) m_num_edge_connections[e.second] );
    
    int local_base_e0 = global_to_local[e.first];
    int local_base_e1 = global_to_local[e.second];
    
    int hidx = particle_hair[e.first];
    HairFlow<DIM>* flow = flows[hidx];
    const VectorXs& eta = flow->getEta();
    const VectorXs& radii_v = flow->getRadiiV();
    
    int nhidx = particle_hair[ edges[eepair.neighbor_eidx].first ];
    HairFlow<DIM>* neighbor_flow = flows[nhidx];
    const VectorXs& neighbor_eta = neighbor_flow->getEta();
    const VectorXs& neighbor_radii_v = neighbor_flow->getRadiiV();
    
    const scalar peta_0 = eta(local_base_e0);
    const scalar peta_1 = eta(local_base_e1);
    const scalar radius_0 = radii_v(local_base_e0);
    const scalar radius_1 = radii_v(local_base_e1);
    const scalar H_0 = peta_0 + radius_0;
    const scalar H_1 = peta_1 + radius_1;
    
    auto& neighbor_e = edges[eepair.neighbor_eidx];
    int local_neighbor_e0 = global_to_local[neighbor_e.first];
    int local_neighbor_e1 = global_to_local[neighbor_e.second];
    
    const scalar vol_p_each = M_PI * mathutils::lerp((H_0 * H_0 - radius_0 * radius_0) / (scalar) num_linked_hairs_e0
                                                     , (H_1 * H_1 - radius_1 * radius_1) / (scalar) num_linked_hairs_e1, eepair.alpha_contact);
    
    // compute actual liquid area shared
    int neighbor_num_linked_hairs_e0 = std::max(1, (int) m_num_edge_connections[neighbor_e.first] );
    int neighbor_num_linked_hairs_e1 = std::max(1, (int) m_num_edge_connections[neighbor_e.second] );
    
    const scalar npeta_0 = neighbor_eta(local_neighbor_e0);
    const scalar npeta_1 = neighbor_eta(local_neighbor_e1);
    const scalar nradius_0 = neighbor_radii_v(local_neighbor_e0);
    const scalar nradius_1 = neighbor_radii_v(local_neighbor_e1);
    const scalar nH_0 = npeta_0 + nradius_0;
    const scalar nH_1 = npeta_1 + nradius_1;
    const scalar vol_q_each = M_PI * mathutils::lerp((nH_0 * nH_0 - nradius_0 * nradius_0) / (scalar) neighbor_num_linked_hairs_e0
                                                     , (nH_1 * nH_1 - nradius_1 * nradius_1) / (scalar) neighbor_num_linked_hairs_e1, eepair.neighbor_local_coord_contact);
    
    const scalar V = vol_p_each + vol_q_each;
    
    if(V < 1e-7) continue;
    
    const Vectors<DIM>& x0 = x.segment<DIM>( m_parent->getDof(e.first) );
    const Vectors<DIM>& x1 = x.segment<DIM>( m_parent->getDof(e.second) );
    const Vectors<DIM>& x2 = x.segment<DIM>( m_parent->getDof(neighbor_e.first) );
    const Vectors<DIM>& x3 = x.segment<DIM>( m_parent->getDof(neighbor_e.second) );
    
    const scalar area_base_e = (x1 - x0).norm();
    const scalar neighbor_area_base_e = (x2 - x3).norm();
    
    Vector3s xcc = eepair.avgpos;
    FluidSim3D* fluid3d = (FluidSim3D*) m_parent->getFluidSim();
    scalar cohesion_weight = 1.0 - fluid3d->get_clamped_particle_weight(xcc);
    scalar criterion = m_parent->getBulkThresholdMultiplier() * fluid3d->cellsize();
    
    PointEdgePair pep;
    pep.alpha_point = eepair.alpha_contact;
    pep.base_eidx = base_eidx;
    pep.quadrature_weight = ( mathutils::clamp(eepair.alpha_1 - eepair.alpha_0, 0.0, 1.0) * area_base_e / mathutils::lerp(num_linked_hairs_e0, num_linked_hairs_e1, pep.alpha_point) + mathutils::clamp(eepair.neighbor_local_coord_1 - eepair.neighbor_local_coord_0, 0.0, 1.0) * neighbor_area_base_e / mathutils::lerp(neighbor_num_linked_hairs_e0, neighbor_num_linked_hairs_e1, pep.neighbor_alpha_point) ) * 0.25;
    // FIX ME: another 0.5 is necessary since we duplicate the detection
    pep.neighbor_alpha_point = eepair.neighbor_local_coord_contact;
    pep.neighbor_eidx = eepair.neighbor_eidx;
    pep.V = V;
    pep.pressure_weight = cohesion_weight;
    pep.hash_code = ((uint64) pep.base_eidx << 32UL) | ((uint64) (pep.neighbor_eidx ) & 32UL);
    pep.time = eepair.time;
    
    poepairs.push_back(pep);
  }
}


template<int DIM>
void PolygonalCohesion<DIM>::findPointEdgePairs(const VectorXs& x, int base_eidx, std::vector<PointEdgePair>& poepairs)
{
  auto& edges = m_parent->getEdges();
  auto& e = edges[base_eidx];
  auto& radii = m_parent->getRadii();
  auto& global_to_local = m_parent->getParticleToHairLocalIndices();
  auto& particle_hair = m_parent->getParticleToHairs();
  auto& flows = m_parent->getFilmFlows();
  
  // get all the hairs linked with the edge
  std::unordered_map<int, std::pair<ParticleEdgePair*, ParticleEdgePair*> > hairs; // neighbor hair index -> Particle-Edge pair connected with that hair
  
  for(auto& hair_pair : m_adjacency_categorized[e.first])
  {
    std::pair<ParticleEdgePair*, ParticleEdgePair*>& pair = hairs[hair_pair.first];
    pair.first = NULL;
    pair.second = NULL;
  }
  
  for(auto& hair_pair : m_adjacency_categorized[e.second])
  {
    std::pair<ParticleEdgePair*, ParticleEdgePair*>& pair = hairs[hair_pair.first];
    pair.first = NULL;
    pair.second = NULL;
  }
  
  for(auto& hair_pair : m_adjacency_categorized[e.first])
  {
    hairs[hair_pair.first].first = hair_pair.second;
  }
  
  for(auto& hair_pair : m_adjacency_categorized[e.second])
  {
    hairs[hair_pair.first].second = hair_pair.second;
  }
  
  // for each hair, check the state for length
  for(auto& hair_pp : hairs)
  {
    const int nhidx = hair_pp.first;
    std::pair<ParticleEdgePair*, ParticleEdgePair*>& pairs = hair_pp.second;
    
    EdgeEdgePair eepair;
    eepair.base_eidx = base_eidx;
    eepair.neighbor_hair_idx = nhidx;
    eepair.alpha_0 = 0.0;
    eepair.alpha_1 = 1.0;
    
    if(pairs.first && pairs.second) // both defined, full length
    {
      // check if both pairs valid
      ParticleEdgePair* pair_0 = pairs.first;
      ParticleEdgePair* pair_1 = pairs.second;
      
      scalar raw_neighbor_coord_0 = computeLocalCoord(pair_0->eidx, pair_0->alpha);
      scalar raw_neighbor_coord_1 = computeLocalCoord(pair_1->eidx, pair_1->alpha);
      
      eepair.count_0 = pair_0->count;
      eepair.count_1 = pair_1->count;
      
      if(pair_0->valid && pair_1->valid) {
        // both valid, do nothing to the alphas, mark them as a part of the double-linked-edges
        pair_0->updated = true;
        pair_1->updated = true;
      } else if(!pair_0->valid && !pair_1->valid) {
        // both invalid, ignore this edge-edge pair, and mark both particle-edge pairs as "should-be-deleted"
        pair_0->should_be_deleted = true;
        pair_1->should_be_deleted = true;
        continue;
      } else if(pair_0->valid) {
        // alpha1 < 1.0
        // mark them as a part of the double-linked-edges
        pair_0->updated = true;
        pair_1->updated = true;
        
        const scalar& r0i = radii(pair_0->pidx);
        const scalar& r0j = pair_0->radii_j;
        const scalar& r1i = radii(pair_1->pidx);
        const scalar& r1j = pair_1->radii_j;
        scalar ddiv = (pair_0->max_dist - (pair_0->dist - (r0i + r0j))) - (pair_1->max_dist - (pair_1->dist - (r1i + r1j)));
        if(ddiv == 0.0) { // parallel
          // do nothing
        } else {
          scalar min_alpha = 0.0;
          if(global_to_local[pair_0->pidx] == 0){
            min_alpha = 0.1;
          }
          eepair.alpha_1 = mathutils::clamp((pair_0->max_dist - (pair_0->dist - (r0i + r0j))) / ddiv, min_alpha, 1.0);
          if(std::isnan(eepair.alpha_1)) {
            std::cerr << "Alpha1 is NAN!" << std::endl;
            std::cerr << pair_0->max_dist << std::endl;
            std::cerr << pair_0->dist << std::endl;
            std::cerr << r0i << std::endl;
            std::cerr << r0j << std::endl;
            std::cerr << ddiv << std::endl;
            exit(-1);
          }
        }
        
      } else if(pair_1->valid) {
        // alpha0 > 0.0
        // mark them as a part of the double-linked-edges
        pair_0->updated = true;
        pair_1->updated = true;
        
        const scalar& r0i = radii(pair_0->pidx);
        const scalar& r0j = pair_0->radii_j;
        const scalar& r1i = radii(pair_1->pidx);
        const scalar& r1j = pair_1->radii_j;
        scalar ddiv = (pair_0->max_dist - (pair_0->dist - (r0i + r0j))) - (pair_1->max_dist - (pair_1->dist - (r1i + r1j)));
        if(ddiv == 0.0) { // parallel
          // do nothing
        } else {
          scalar max_alpha = 1.0;
          if(global_to_local[pair_1->pidx] == (int) flows[particle_hair[pair_1->pidx]]->size() - 1){
            max_alpha = 0.9;
          }
          eepair.alpha_0 = mathutils::clamp((pair_0->max_dist - (pair_0->dist - (r0i + r0j))) / ddiv, 0.0, max_alpha);
          if(std::isnan(eepair.alpha_0)) {
            std::cerr << "Alpha0 is NAN!" << std::endl;
            std::cerr << pair_0->max_dist << std::endl;
            std::cerr << pair_0->dist << std::endl;
            std::cerr << r0i << std::endl;
            std::cerr << r0j << std::endl;
            std::cerr << ddiv << std::endl;
            exit(-1);
          }
        }
      }
      
      eepair.neighbor_local_coord_0 = mathutils::lerp(raw_neighbor_coord_0, raw_neighbor_coord_1, eepair.alpha_0);
      eepair.neighbor_local_coord_1 = mathutils::lerp(raw_neighbor_coord_0, raw_neighbor_coord_1, eepair.alpha_1);
    } else if(pairs.first && !pairs.second)
    {
      // alpha1 < 1, computed with extrapolation
      ParticleEdgePair* pair_0 = pairs.first;
      if(!pair_0->valid) {
        // mark pair_0 as "should-be-deleted"
        pair_0->should_be_deleted = true;
        continue;
      }
      auto& neighbor_e = edges[pair_0->eidx];
      Vectors<DIM> gap_vec;
      scalar neighbor_alpha;
      mathutils::pointedgevec<scalar, DIM>(x.segment<DIM>( m_parent->getDof( e.second ) ), x.segment<DIM>(m_parent->getDof( neighbor_e.first ) ), x.segment<DIM>( m_parent->getDof( neighbor_e.second ) ), gap_vec, neighbor_alpha);
      scalar pair_1_dist = gap_vec.norm();
      scalar ddiv = pair_1_dist - pair_0->dist;
      if(ddiv == 0.0) { // parallel
        // do nothing
      } else {
        const scalar& r0i = radii(pair_0->pidx);
        const scalar& r0j = pair_0->radii_j;
        scalar min_alpha = 0.0;
        if(global_to_local[pair_0->pidx] == 0){
          min_alpha = 0.1;
        }
        eepair.alpha_1 = mathutils::clamp((pair_0->max_dist - (pair_0->dist - (r0i + r0j))) / ddiv, min_alpha, 1.0);
        if(std::isnan(eepair.alpha_1)) {
          std::cerr << "Alpha1 is NAN*!" << std::endl;
          std::cerr << pair_0->max_dist << std::endl;
          std::cerr << pair_0->dist << std::endl;
          std::cerr << r0i << std::endl;
          std::cerr << r0j << std::endl;
          std::cerr << ddiv << std::endl;
          exit(-1);
        }
      }
      eepair.neighbor_local_coord_0 = computeLocalCoord(pair_0->eidx, pair_0->alpha);
      eepair.neighbor_local_coord_1 = computeLocalCoord(pair_0->eidx, neighbor_alpha);
      eepair.count_0 = eepair.count_1 = pair_0->count;
      
    } else if(!pairs.first && pairs.second)
    {
      // alpha0 > 0, computed with extrapolation
      ParticleEdgePair* pair_1 = pairs.second;
      if(!pair_1->valid) {
        // mark pair_1 as "should-be-deleted"
        pair_1->should_be_deleted = true;
        continue;
      }
      auto& neighbor_e = edges[pair_1->eidx];
      Vectors<DIM> gap_vec;
      scalar neighbor_alpha;
      mathutils::pointedgevec<scalar, DIM>(x.segment<DIM>( m_parent->getDof( e.first ) ), x.segment<DIM>( m_parent->getDof( neighbor_e.first ) ), x.segment<DIM>( m_parent->getDof( neighbor_e.second ) ), gap_vec, neighbor_alpha);
      scalar pair_0_dist = gap_vec.norm();
      scalar ddiv = pair_1->dist - pair_0_dist;
      if(ddiv == 0.0) { // parallel
        // do nothing
      } else {
        const scalar& r1i = radii(pair_1->pidx);
        const scalar& r1j = pair_1->radii_j;
        scalar max_alpha = 1.0;
        if(global_to_local[pair_1->pidx] == (int) flows[particle_hair[pair_1->pidx]]->size() - 1){
          max_alpha = 0.9;
        }
        eepair.alpha_0 = mathutils::clamp((pair_1->max_dist - (pair_0_dist - (r1i + r1j))) / ddiv, 0.0, max_alpha);
        if(std::isnan(eepair.alpha_0)) {
          std::cerr << "Alpha0 is NAN*!" << std::endl;
          std::cerr << pair_0_dist << std::endl;
          std::cerr << pair_1->max_dist << std::endl;
          std::cerr << pair_1->dist << std::endl;
          std::cerr << r1i << std::endl;
          std::cerr << r1j << std::endl;
          std::cerr << ddiv << std::endl;
          exit(-1);
        }
      }
      eepair.neighbor_local_coord_0 = computeLocalCoord(pair_1->eidx, neighbor_alpha);
      eepair.neighbor_local_coord_1 = computeLocalCoord(pair_1->eidx, pair_1->alpha);
      eepair.count_0 = eepair.count_1 = pair_1->count;
      
    } else {
      std::cerr << "ERROR: found undefined pairs!" << base_eidx << "[" << nhidx << "]->(" << pairs.first << ", " << pairs.second << ")" << std::endl;
      exit(EXIT_FAILURE);
    }
    
    int num_linked_hairs_e0 = std::max(1, (int) m_adjacency_categorized[e.first].size());
    int num_linked_hairs_e1 = std::max(1, (int) m_adjacency_categorized[e.second].size());
    
    int local_base_e0 = global_to_local[e.first];
    int local_base_e1 = global_to_local[e.second];
    
    int hidx = particle_hair[e.first];
    
    HairFlow<DIM>* flow = m_parent->getFilmFlows()[hidx];
    const VectorXs& eta = flow->getEta();
    
    HairFlow<DIM>* neighbor_flow = m_parent->getFilmFlows()[nhidx];
    const VectorXs& neighbor_eta = neighbor_flow->getEta();
    
    const scalar peta_0 = eta(local_base_e0);
    const scalar peta_1 = eta(local_base_e1);
    const scalar& radius_0 = radii(e.first);
    const scalar& radius_1 = radii(e.second);
    const scalar H_0 = peta_0 + radius_0;
    const scalar H_1 = peta_1 + radius_1;
    
    for(int k_gauss = 0; k_gauss < m_num_quadrature; ++k_gauss)
    {
      const scalar alpha_gauss = mathutils::lerp(eepair.alpha_0, eepair.alpha_1, mathutils::gauss_legendre_point[m_num_quadrature - 1][k_gauss]);
      
      const scalar neighbor_local_coord_gauss = mathutils::lerp(eepair.neighbor_local_coord_0, eepair.neighbor_local_coord_1, mathutils::gauss_legendre_point[m_num_quadrature - 1][k_gauss]);
      
      int neighbor_eidx;
      scalar neighbor_alpha;
      extractLocalCoord(neighbor_local_coord_gauss, eepair.neighbor_hair_idx, neighbor_eidx, neighbor_alpha);
      
      auto& neighbor_e = edges[neighbor_eidx];
      int local_neighbor_e0 = global_to_local[neighbor_e.first];
      int local_neighbor_e1 = global_to_local[neighbor_e.second];
      
      const scalar vol_p_each = M_PI * mathutils::lerp((H_0 * H_0 - radius_0 * radius_0) / (scalar) num_linked_hairs_e0,
                                                       (H_1 * H_1 - radius_1 * radius_1) / (scalar) num_linked_hairs_e0, alpha_gauss);
      
      // compute actual liquid area shared
      int neighbor_num_linked_hairs_e0 = std::max(1, (int) m_adjacency_categorized[neighbor_e.first].size());
      int neighbor_num_linked_hairs_e1 = std::max(1, (int) m_adjacency_categorized[neighbor_e.second].size());
      
      const scalar npeta_0 = neighbor_eta(local_neighbor_e0);
      const scalar npeta_1 = neighbor_eta(local_neighbor_e1);
      const scalar nradius_0 = radii(neighbor_e.first);
      const scalar nradius_1 = radii(neighbor_e.second);
      const scalar nH_0 = npeta_0 + nradius_0;
      const scalar nH_1 = npeta_1 + nradius_1;
      
      const scalar vol_q_each = M_PI * mathutils::lerp((nH_0 * nH_0 - nradius_0 * nradius_0) / (scalar) neighbor_num_linked_hairs_e0,
                                                       (nH_1 * nH_1 - nradius_1 * nradius_1) / (scalar) neighbor_num_linked_hairs_e1, neighbor_alpha);
      
      const scalar V = vol_p_each + vol_q_each;
      
      const Vectors<DIM>& x0 = x.segment<DIM>(m_parent->getDof( e.first ));
      const Vectors<DIM>& x1 = x.segment<DIM>(m_parent->getDof( e.second ));
      
      const scalar area_base_e = (x1 - x0).norm();
      
      scalar cohesion_weight = 1.0;
      if(DIM == 2) {
        Vector2s xcc = (x.segment<2>(e.first * 2) + x.segment<2>(e.second * 2) + x.segment<2>(neighbor_e.first * 2) + x.segment<2>(neighbor_e.second * 2)) * 0.25;
        FluidSim2D* fluid2d = (FluidSim2D*) m_parent->getFluidSim();
        cohesion_weight = 1.0 - fluid2d->get_clamped_particle_weight(xcc);
      } else if(DIM == 3) {
        Vector3s xcc = (x.segment<3>(m_parent->getDof( e.first )) + x.segment<3>(m_parent->getDof( e.second )) + x.segment<3>(m_parent->getDof( neighbor_e.first )) + x.segment<3>(m_parent->getDof( neighbor_e.second ))) * 0.25;
        FluidSim3D* fluid3d = (FluidSim3D*) m_parent->getFluidSim();
        cohesion_weight = 1.0 - fluid3d->get_clamped_particle_weight(xcc);
      }
      
      PointEdgePair pep;
      pep.alpha_point = alpha_gauss;
      pep.base_eidx = base_eidx;
      pep.quadrature_weight = mathutils::gauss_legendre_weight[m_num_quadrature - 1][k_gauss] * (eepair.alpha_1 - eepair.alpha_0) * mathutils::lerp(eepair.count_0, eepair.count_1, neighbor_alpha) * area_base_e / mathutils::lerp(num_linked_hairs_e0, num_linked_hairs_e1, alpha_gauss);
      pep.neighbor_alpha_point = neighbor_alpha;
      pep.neighbor_eidx = neighbor_eidx;
      pep.V = V;
      pep.pressure_weight = cohesion_weight;
      pep.hash_code = ((uint64) pep.base_eidx << 32UL) | ((uint64) (pep.neighbor_eidx * m_num_quadrature + k_gauss) & 32UL);
      pep.time = 0.0;
      
      poepairs.push_back(pep);
    }
  }
}

template<int DIM>
uint64_t PolygonalCohesion<DIM>::makePairwisePPHash(int i, int j)
{
  return ((uint64_t) std::min(i, j) << 32UL) | (uint64_t) std::max(i, j);
}

template<int DIM>
bool PolygonalCohesion<DIM>::checkPairwisePPHash(int i, int j)
{
  int np = m_parent->getNumParticles();
  
  return m_pp_pair_hash[j].find(i) != m_pp_pair_hash[j].end();
}

template<int DIM>
void PolygonalCohesion<DIM>::markPairwisePPHash(int i, int j)
{
  int np = m_parent->getNumParticles();
  
  m_pp_pair_hash[i].insert(j);
}

template<int DIM>
void PolygonalCohesion<DIM>::write(std::vector<scalar>& data) const
{
  
}

template<int DIM>
void PolygonalCohesion<DIM>::read(const scalar* data)
{
  
}

template<int DIM>
size_t PolygonalCohesion<DIM>::adjacency_size() const
{
  int np = m_adjacency_categorized.size();
  size_t sum = 0;
  for(int i = 0; i < np; ++i)
  {
    sum += m_adjacency_categorized[i].size();
  }
  
  return sum;
}

template<int DIM>
void PolygonalCohesion<DIM>::extractLocalCoord(const scalar& local_coord, int hidx, int& eidx, scalar& alpha) const
{
  HairFlow<DIM>* flow = m_parent->getFilmFlows()[hidx];
  auto& edge_indices = flow->getEdgeIndices();
  
  int local_eidx = std::min((int) edge_indices.size() - 1, (int) local_coord);
  eidx = edge_indices[local_eidx];
  alpha = local_coord - (scalar) local_eidx;
}

template<int DIM>
void PolygonalCohesion<DIM>::updatePointEdgePairs(const VectorXs& x, const VectorXs& v, const scalar& dt )
{ // X is end of timestep ( x + v*dt )
  const int ne = m_parent->getNumEdges();
  const int np = m_parent->getNumParticles();
  
  if(!np) return;
  
  if(m_compute_particle_poe_mapping) {
    for(int i = 0; i < np; ++i)
    {
      m_particle_to_point_edge_pairs[i].resize(0);
    }
  }
  
  if((int) m_point_edge_pairs_cache.size() != ne)
    m_point_edge_pairs_cache.resize(ne);
  
  if((int) m_counting_poe_pair_location.size() != ne)
    m_counting_poe_pair_location.resize(ne);
  
  if(m_parent->useCtcd()) {
    if((int) m_edge_connections.size() != ne)
      m_edge_connections.resize(ne);
    
    if((int) m_num_edge_connections.size() != np)
      m_num_edge_connections.resize(np);
  }
  
  if( DIM == 3 && m_parent->useCtcd() ){
    
    // Detect EdgeEdge CCD
    threadutils::thread_pool::ParallelFor(0, ne, [&] (int i){
      findEdgeEdgeContact(x, v, dt, i);
    });
    
    // this needed to get moved here since we did not have any contacts to work with in 3D earlier
    m_particle_particle_pairs.resize(0);
    
    findParticleParticlePairsEEC(x);
    
    m_point_edge_pairs.resize(0);
    // for each edge, sample Gauss quadrature
    threadutils::thread_pool::ParallelFor(0, ne, [&] (int i){
      std::vector<PointEdgePair>& poepairs = m_point_edge_pairs_cache[i];
      poepairs.resize(0);
      
      findPointEdgePairsEEC(x, i, poepairs);
      m_counting_poe_pair_location[i] = poepairs.size();
    });
  }
  else{
    m_point_edge_pairs.resize(0);
    // for each edge, sample Gauss quadrature
    threadutils::thread_pool::ParallelFor(0, ne, [&] (int i){
      std::vector<PointEdgePair>& poepairs = m_point_edge_pairs_cache[i];
      poepairs.resize(0);
      
      findPointEdgePairs(x, i, poepairs);
      m_counting_poe_pair_location[i] = poepairs.size();
    });
  }
  
  std::partial_sum(m_counting_poe_pair_location.begin(), m_counting_poe_pair_location.end(), m_counting_poe_pair_location.begin());
  
  const int npoe = m_counting_poe_pair_location[m_counting_poe_pair_location.size() - 1];
  m_point_edge_pairs.resize(npoe);
  
  threadutils::thread_pool::ParallelFor(0, ne, [&] (int i){
    int base_idx = (i == 0) ? 0 : m_counting_poe_pair_location[i - 1];
    std::vector<PointEdgePair>& poepairs = m_point_edge_pairs_cache[i];
    int k = 0;
    for(auto& pair : poepairs) {
      m_point_edge_pairs[base_idx + k] = pair;
      k++;
    }
    poepairs.resize(0);
  });
  
#ifdef VERBOSE_CONN
  for(int i = 0; i < npoe; ++i)
  {
    std::cout << "POE [" << i << "] = [" << m_point_edge_pairs[i] << "]" << std::endl;
  }
#endif
  
  if(m_compute_particle_poe_mapping) {
    // splat the POEP template onto Particles
    auto& edges = m_parent->getEdges();
    for(int i = 0; i < npoe; ++i)
    {
      const PointEdgePair& poepair = m_point_edge_pairs[i];
      const std::pair<int, int>& be = edges[poepair.base_eidx];
      const std::pair<int, int>& ne = edges[poepair.neighbor_eidx];
      m_particle_to_point_edge_pairs[be.first].push_back(i);
      m_particle_to_point_edge_pairs[be.second].push_back(i);
      m_particle_to_point_edge_pairs[ne.first].push_back(i);
      m_particle_to_point_edge_pairs[ne.second].push_back(i);
    }
  }
}

template<int DIM>
scalar PolygonalCohesion<DIM>::computeLocalCoord(int eidx, const scalar& alpha) const
{
  auto& global_local = m_parent->getParticleToHairLocalIndices();
  auto& edges = m_parent->getEdges();
  auto& e = edges[eidx];
  int base_local_idx = global_local[e.first];
  return (scalar) base_local_idx + alpha;
}

template<int DIM>
void PolygonalCohesion<DIM>::computeParticleEdgeCounts(int pidx)
{
  std::unordered_map<int, ParticleEdgePair*>& pepairs = m_adjacency_categorized[pidx];
  auto& edges = m_parent->getEdges();
  auto& flows = m_parent->getFilmFlows();
  auto& global_local = m_parent->getParticleToHairLocalIndices();
  auto& particle_hair = m_parent->getParticleToHairs();
  
  for (auto itr = pepairs.begin(); itr != pepairs.end(); ++itr) {
    ParticleEdgePair* pep = itr->second;
    int hidx = particle_hair[pep->pidx];
    HairFlow<DIM>* flow = flows[hidx];
    auto& hair_edges = flow->getEdgeIndices();
    int nlocale = (int) hair_edges.size();
    int local_pidx = global_local[pep->pidx];
    int local_eidx_low = mathutils::clamp(local_pidx - 1, 0, nlocale - 1);
    int local_eidx_high = mathutils::clamp(local_pidx, 0, nlocale - 1);
    int eidx_low = hair_edges[local_eidx_low];
    int eidx_high = hair_edges[local_eidx_high];
    
    
    auto& neighbor_e = edges[pep->eidx];
    int neighbor_pidxs[] = {neighbor_e.first, neighbor_e.second};
    
    bool find_inversed = false;
    
    for(int i = 0; i < 2; ++i) {
      const std::unordered_map<int, ParticleEdgePair*>& adj_pairs = m_adjacency_categorized[neighbor_pidxs[i]];
      auto itr = adj_pairs.find(hidx);
      if(itr != adj_pairs.end()) {
        // found current hair in neighbor hair
        const ParticleEdgePair* inv_pe_pair = itr->second;
        if((inv_pe_pair->eidx == eidx_low && inv_pe_pair->alpha >= 0.5) ||
           (inv_pe_pair->eidx == eidx_high && inv_pe_pair->alpha <= 0.5)) {
          find_inversed = true;
          break;
        }
      }
    }
    
    if(find_inversed) {
      pep->count = 0.5;
    } else {
      pep->count = 1.0;
    }
  }
}

//#define VERBOSE_CONN

template<int DIM>
void PolygonalCohesion<DIM>::computeConnections(const VectorXs& x)
{
  int np = m_parent->getNumParticles();
  int ne = m_parent->getNumEdges();
  
  if((int) m_pp_pair_hash.size() != np) m_pp_pair_hash.resize(np);
  
  if((int) m_adjacency_categorized.size() != np)
    m_adjacency_categorized.resize(np);
  
  if((int) m_num_adjacency_categorized.size() != np)
    m_num_adjacency_categorized.resize(np);
  
  if((int) m_particle_to_point_edge_pairs.size() != np)
    m_particle_to_point_edge_pairs.resize(np);
  
  threadutils::thread_pool::ParallelFor(0, np, [&] (int i) {
    auto& adjs = m_adjacency_categorized[i];
    // only count valid adjs
    int sum = 0;
    for(auto& adj : adjs)
    {
      if(adj.second->valid) ++sum;
    }
    
    m_num_adjacency_categorized[i] = sum;
  });
  
  // for each particle, get all the edges connected with, and build particle -> hair -> pairs
  threadutils::thread_pool::ParallelFor(0, np, [&] (int i)
                                        {
                                          findParticleEdgePairs(x, i);
                                        });
  
  threadutils::thread_pool::ParallelFor(0, np, [&] (int i)
                                        {
                                          computeParticleEdgeCounts(i);
                                        });
#ifdef VERBOSE_CONN
  for(int i = 0; i < npe; ++i)
  {
    std::cout << "PAE [" << i << "] = [" << m_particle_edge_pairs[i] << "]" << std::endl;
  }
#endif
  
  m_particle_particle_pairs.resize(0);
  
  if(!m_parent->isIndividualTransfer())
    findParticleParticlePairs(x);
}

template<int DIM>
void PolygonalCohesion<DIM>::buildSearchTree(const VectorXs& x)
{
  if( m_parent->useCtcd() ) m_parent->updateSearchRadius();
  
  const scalar cellsize = m_parent->getSearchRadius();
  
  m_parent->updateBoundingBox();
  
  const Vectors<DIM>& bbx_min = m_parent->getBoundingBoxMin();
  const Vectors<DIM>& bbx_max = m_parent->getBoundingBoxMax();
  
  int nindex[3] = {1, 1, 1};
  for(int r = 0; r < DIM; ++r)
  {
    nindex[r] = std::max(1, (int) ceil((bbx_max(r) - bbx_min(r)) / cellsize));
  }
  
  m_sorter->resize(nindex[0], nindex[1], nindex[2]);
  
  int ne = m_parent->getNumEdges();
  auto& edges = m_parent->getEdges();
  
  m_sorter->sort(ne, [&] (int eidx, int& i, int& j, int& k)
                 {
                   if(eidx < 0 || eidx >= edges.size()) return;
                   
                   auto& e = edges[eidx];
                   Vectors<DIM> xc = (x.segment<DIM>( m_parent->getDof( e.first ) ) + x.segment<DIM>( m_parent->getDof( e.second ) )) * 0.5;
                   
                   i = max(0, min((int)((xc(0) - bbx_min(0)) / cellsize), nindex[0]-1));
                   j = max(0, min((int)((xc(1) - bbx_min(1)) / cellsize), nindex[1]-1));
                   k = max(0, min((int)((xc(2) - bbx_min(2)) / cellsize), nindex[2]-1));
                 });
}


template<int DIM>
void PolygonalCohesion<DIM>::computeInterHairVariables()
{
  int M_inter = m_particle_particle_pairs.size();
  int N = m_parent->getNumParticles();
  const VectorXs& x = m_parent->getX();
  const VectorXs& radius = m_parent->getRadii();
  const std::vector< std::pair<int, int> >& edges = m_parent->getEdges();
  
  int M_local = edges.size();
  int M_global = M_inter + M_local;
  
  m_ce_inter_short_buffer.resize(M_inter);
  m_ce_inter_buffer.resize(M_inter * DIM);
  m_ce_global_buffer.resize(M_global * DIM);
  
  // compute area_e_global
  m_area_e_global.resize(M_global);
  m_area_e_global.setZero();
  
  tbb::parallel_for(0, M_local, 1, [&] (int i) {
    auto& e = edges[i];
    m_area_e_global(i) = (x.segment<DIM>( m_parent->getDof( e.first ) ) - x.segment<DIM>( m_parent->getDof( e.second ) )).norm();
  });
  stringutils::print(MAKE_REF(m_area_e_global));
  
  // compute G_f_inter & dirf_inter
  m_G_f_interhair.resize(M_inter * DIM);
  m_G_f_global.resize(M_global * DIM);
  m_dir_f_interhair.resize(M_inter, M_inter * DIM);
  TripletXs dir_f_interhair;
  dir_f_interhair.resize(M_inter * DIM);
  
  tbb::parallel_for(0, M_local, 1, [&] (int i) {
    for(int r = 0; r < DIM; ++r)
    {
      m_G_f_global(i * DIM + r) = m_area_e_global(i);
    }
  });
  
  tbb::parallel_for(0, M_inter, 1, [&] (int i) {
    const ParticleParticlePair& ppp = m_particle_particle_pairs[i];
    VectorXs dx = x.segment<DIM>( m_parent->getDof( ppp.pidx[1] ) ) - x.segment<DIM>( m_parent->getDof( ppp.pidx[0] ) );
    scalar ddx = ppp.d;
    m_area_e_global(M_local + i) = ddx;
    
    VectorXs ndx = dx / ddx;
    for(int r = 0; r < DIM; ++r)
    {
      m_G_f_interhair(i * DIM + r) = ddx;
      m_G_f_global[(i + M_local) * DIM + r] = ddx;
      dir_f_interhair[i * DIM + r] = Triplets(i, i * DIM + r, ndx(r));
    }
  });
  stringutils::print(MAKE_REF(m_area_e_global));
  m_dir_f_interhair.setFromTriplets(dir_f_interhair.begin(), dir_f_interhair.end()); stringutils::print(MAKE_REF(m_dir_f_interhair));
  
  // compute vertex area
  m_W_fv_interhair_T.resize(M_inter, N);
  m_pplink_count.resize(N);
  memset(&m_pplink_count[0], 0, N * sizeof(int));
  
  m_area_v_global.resize(N);
  m_area_v_global.setZero();
  for(int i = 0; i < M_inter; ++i)
  {
    const ParticleParticlePair& ppp = m_particle_particle_pairs[i];
    m_pplink_count[ppp.pidx[0]]++;
    m_pplink_count[ppp.pidx[1]]++;
    m_area_v_global(ppp.pidx[0]) += ppp.d;
    m_area_v_global(ppp.pidx[1]) += ppp.d;
  }
  const std::vector< std::vector<int> >& particle_edges = m_parent->getParticleToEdge();
  
  tbb::parallel_for(0, N, 1, [&] (int i) {
    m_pplink_count[i] += particle_edges[i].size();
    for(int eidx : particle_edges[i])
    {
      const scalar& d = m_area_e_global(eidx);
      m_area_v_global(i) += d;
    }
  });
  
  m_area_v_global *= 0.5;
  
  stringutils::print(MAKE_REF(m_area_v_global));
  
  // compute Gv^-1 global
  m_iG_v_global.resize(N);
  for(int i = 0; i < N; ++i)
  {
    m_iG_v_global[i] = 1.0 / m_area_v_global[i];
  }
  stringutils::print(MAKE_REF(m_iG_v_global));
  
  TripletXs W_fv_interhair_T;
  W_fv_interhair_T.resize(M_inter * 2);
  // compute W_fv_inter
  tbb::parallel_for(0, M_inter, 1, [&] (int i) {
    const ParticleParticlePair& ppp = m_particle_particle_pairs[i];
    const scalar& A_Fj = m_area_e_global(M_local + i);
    const scalar& A_Vi0 = m_area_v_global(ppp.pidx[0]);
    const scalar& A_Vi1 = m_area_v_global(ppp.pidx[1]);
    W_fv_interhair_T[i * 2 + 0] = Triplets(i, ppp.pidx[0], A_Fj / (A_Vi0 * 2));
    W_fv_interhair_T[i * 2 + 1] = Triplets(i, ppp.pidx[1], A_Fj / (A_Vi1 * 2));
  });
  m_W_fv_interhair_T.setFromTriplets(W_fv_interhair_T.begin(), W_fv_interhair_T.end()); stringutils::print(MAKE_REF(W_fv_interhair_T));
  
  // compute gradF_inter
  m_gradF_interhair.resize(M_inter * DIM, N);
  m_gradF_global.resize(M_global * DIM, N);
  TripletXs gradF_global;
  TripletXs gradF_interhair;
  gradF_global.resize(M_global * DIM * 2);
  gradF_interhair.resize(M_inter * DIM * 2);
  
  tbb::parallel_for(0, M_local, 1, [&] (int i) {
    auto& e = edges[i];
    const scalar& A_Fj = m_area_e_global(i);
    int i0 = e.first;
    int i1 = e.second;
    const VectorXs& x0 = x.segment<DIM>( m_parent->getDof( i0 ) );
    const VectorXs& x1 = x.segment<DIM>( m_parent->getDof( i1 ) );
    VectorXs dx = x1 - x0;
    
    scalar coeff = 1.0 / (A_Fj * A_Fj);
    
    for(int r = 0; r < DIM; ++r)
    {
      int k = i * DIM + r;
      gradF_global[k * 2 + 0] = Triplets(i * DIM + r, i0, -coeff * dx(r));
      gradF_global[k * 2 + 1] = Triplets(i * DIM + r, i1, coeff * dx(r));
    }
  });
  
  tbb::parallel_for(0, M_inter, 1, [&] (int i) {
    const ParticleParticlePair& ppp = m_particle_particle_pairs[i];
    const scalar& A_Fj = m_area_e_global(M_local + i);
    int i0 = ppp.pidx[0];
    int i1 = ppp.pidx[1];
    const VectorXs& x0 = x.segment<DIM>( m_parent->getDof( i0 ) );
    const VectorXs& x1 = x.segment<DIM>( m_parent->getDof( i1 ) );
    VectorXs dx = x1 - x0;
    
    scalar coeff = 1.0 / (A_Fj * A_Fj);
    
    for(int r = 0; r < DIM; ++r)
    {
      int k = i * DIM + r;
      gradF_interhair[k * 2 + 0] = Triplets(i * DIM + r, i0, -coeff * dx(r));
      gradF_interhair[k * 2 + 1] = Triplets(i * DIM + r, i1, coeff * dx(r));
      
      int kl = (M_local + i) * DIM + r;
      gradF_global[kl * 2 + 0] = Triplets((M_local + i) * DIM + r, i0, -coeff * dx(r));
      gradF_global[kl * 2 + 1] = Triplets((M_local + i) * DIM + r, i1, coeff * dx(r));
    }
  });
  m_gradF_interhair.setFromTriplets(gradF_interhair.begin(), gradF_interhair.end()); stringutils::print(MAKE_REF(m_gradF_interhair));
  m_gradF_global.setFromTriplets(gradF_global.begin(), gradF_global.end()); stringutils::print(MAKE_REF(m_gradF_global));
  
  m_gradF_interhair_T = m_gradF_interhair.transpose();
  m_gradF_global_T = m_gradF_global.transpose();
  m_dir_f_interhair_T = m_dir_f_interhair.transpose();
  
  m_gradF_global.makeCompressed();
  m_gradF_global_T.makeCompressed();
  m_gradF_interhair.makeCompressed();
  m_gradF_interhair_T.makeCompressed();
  m_dir_f_interhair.makeCompressed();
  m_dir_f_interhair_T.makeCompressed();
  m_W_fv_interhair_T.makeCompressed();
}

template<int DIM>
void PolygonalCohesion<DIM>::computeGlobalHairVariables()
{
}

template<int DIM>
void PolygonalCohesion<DIM>::computeGlobalHairPressure()
{
  // extract kappa from local to global
  int N = m_parent->getNumParticles();
  const VectorXs& v = m_parent->getV();
  
  auto& flows = m_parent->getFilmFlows();
  int nf = flows.size();
  
  m_cur_eta_v_global.resize(N);
  m_cv_buffer.resize(N);
  
  if(m_parent->isMassSpring()) {
    m_rhs_offset_v_global.resize(N, DIM + 1);
    m_cur_hhrr_v_global.resize(N, DIM + 1);
  } else {
    m_rhs_offset_v_global.resize(N, DIM + 2);
    m_cur_hhrr_v_global.resize(N, DIM + 2);
  }
  
  m_pressure_v_global.resize(N);
  
  m_cur_eta_v_global.setZero();
  m_cur_hhrr_v_global.setZero();
  
  tbb::parallel_for(0, nf, 1, [&] (int i) {
    const HairFlow<DIM>* flow = flows[i];
    const VectorXs& eta = flow->getEta();
    const VectorXs& radii_v = flow->getRadiiV();
    const std::vector<int>& indices = flow->getParticleIndices();
    const int nfp = indices.size();
    for(int j = 0; j < nfp; ++j)
    {
      int pidx = indices[j];
      int idof = m_parent->getDof(pidx);
      
      m_cur_eta_v_global[pidx] = eta[j];
      
      scalar A = (eta[j] + radii_v[j]) * (eta[j] + radii_v[j]) - radii_v[j] * radii_v[j];
      m_cur_hhrr_v_global(pidx, 0) = A;
      m_cur_hhrr_v_global(pidx, 1) = v(idof + 0) * A;
      m_cur_hhrr_v_global(pidx, 2) = v(idof + 1) * A;
      m_cur_hhrr_v_global(pidx, 3) = v(idof + 2) * A;
      
      if(!m_parent->isMassSpring()) {
        if(j == nfp - 1) {
          m_cur_hhrr_v_global(pidx, 4) = 0.0;
        } else {
          m_cur_hhrr_v_global(pidx, 4) = v(idof + 3) * A * (A / 2.0 + radii_v[j] * radii_v[j]);
        }
      }
    }
  });
  
  stringutils::print(MAKE_REF(m_cur_eta_v_global));
  stringutils::print(MAKE_REF(m_cur_hhrr_v_global));
  
  // compute pressure
  const scalar sigma = m_parent->getLiquidTension();
  //m_cur_eta_v_global * m_gradF_global^T * [m_G_f_global]
  mathutils::computeJTPhi_coeff(m_ce_global_buffer, m_G_f_global, m_cur_eta_v_global, m_gradF_global_T);
  //m_cur_eta_v_global * m_gradF_global^T * [m_G_f_global] * m_gradF_global * [m_iG_v_global] + 2.0 * m_kappa_v_global
  mathutils::computeJTPhi_coeff(m_pressure_v_global, m_iG_v_global, m_ce_global_buffer, m_gradF_global);
  m_pressure_v_global *= sigma;
  stringutils::print(MAKE_REF(m_pressure_v_global));
}

template<int DIM>
void PolygonalCohesion<DIM>::computeInterHairVelocity(const scalar& dt)
{
  int M_inter = m_particle_particle_pairs.size();
  const scalar& rho = m_parent->getLiquidDensity();
  const std::vector<int>& particle_to_hairs = m_parent->getParticleToHairs();
  const std::vector<HairFlow<DIM>*>& flows = m_parent->getFilmFlows();
  const std::vector<int>& global_local = m_parent->getParticleToHairLocalIndices();
  
  m_u_interhair.resize(M_inter);
  m_u_interhair.setZero();
  const scalar visc = m_parent->getViscosity();
  
  tbb::parallel_for(0, M_inter, 1, [&] (int i) {
    const ParticleParticlePair& ppp = m_particle_particle_pairs[i];
    int hidx0 = particle_to_hairs[ppp.pidx[0]];
    int hidx1 = particle_to_hairs[ppp.pidx[1]];
    const VectorXs& porosity0 = flows[hidx0]->getPorosity();
    const VectorXs& porosity1 = flows[hidx1]->getPorosity();
    int local0 = global_local[ppp.pidx[0]];
    int local1 = global_local[ppp.pidx[1]];
    
    scalar Sa = std::max(1e-63, mathutils::lerp(porosity0(local0), porosity1(local1), 0.5));
    scalar phi = std::max(1 - Sa, 1e-63);
    scalar beta = std::max(0.0, 4 * M_PI / (-log(phi)-1.476+2.*phi-0.5*phi*phi));
    
    scalar u = -(m_pressure_v_global(ppp.pidx[1]) - m_pressure_v_global(ppp.pidx[0])) * dt / rho;
    u *= 1.0 / (1.0 + beta * visc * dt + 0.14 * visc * rho * pow(Sa, -2.5) * sqrt(beta) * fabs(u * dt));
    
    m_u_interhair(i) = u;
  });
  
  stringutils::print(MAKE_REF(m_u_interhair));
}

template<int DIM>
void PolygonalCohesion<DIM>::computeInterHairRHS(const scalar& dt)
{
  for(int i = 0; i < m_cur_hhrr_v_global.cols(); ++i) {
    computeJTPhi(m_ce_inter_buffer, m_cur_hhrr_v_global.col(i), m_gradF_interhair_T);
    computeJTPhi_coeff(m_ce_inter_short_buffer, m_u_interhair, m_ce_inter_buffer, m_dir_f_interhair_T);
    computeJTPhi_coeff(m_rhs_offset_v_global, -dt, m_ce_inter_short_buffer, m_W_fv_interhair_T, i);
    
    computeJTPhi_coeff(m_ce_inter_buffer, m_G_f_interhair, m_u_interhair, m_dir_f_interhair);
    computeJTPhi_coeff(m_cv_buffer, m_iG_v_global, m_ce_inter_buffer, m_gradF_interhair);
    accumulate_cwiseProduct(m_cv_buffer, m_cur_hhrr_v_global.col(i));
    m_rhs_offset_v_global.col(i) += m_cv_buffer * dt;
  }
}


template<int DIM>
void PolygonalCohesion<DIM>::updatePorosity()
{
  const std::vector<HairFlow<DIM>*>& flows = m_parent->getFilmFlows();
  const int nf = flows.size();
  const VectorXs& radii = m_parent->getRadii();
  const std::vector<int>& particle_to_hairs = m_parent->getParticleToHairs();
  const std::vector<int>& particle_local_indices = m_parent->getParticleToHairLocalIndices();
  const std::vector< std::pair<int, int> >& edges = m_parent->getEdges();
  
  threadutils::thread_pool::ParallelFor(0, nf, [&] (int hidx) {
    HairFlow<DIM>* hair = flows[hidx];
    auto& particles = hair->getParticleIndices();
    VectorXs& eta = hair->getEta();
    VectorXs& porosity = hair->getPorosity();
    VectorXs& avg_eta = hair->getAvgEta();
    
    VectorXs weights(porosity.size());
    weights.setZero();
    
    porosity.setZero();
    
    int np = particles.size();
    
    for(int i = 0; i < np; ++i)
    {
      int pidx = particles[i];
      scalar accum_liquid = M_PI * eta(i) * eta(i);
      scalar accum_solid = M_PI * radii(pidx) * radii(pidx);
      
      scalar max_r = sqrt(accum_liquid / M_PI);
      int count = 1;
      
      const std::unordered_map<int, ParticleEdgePair*>& adjs = m_adjacency_categorized[pidx];
      for(auto& adj : adjs){
        const ParticleEdgePair& pair = *(adj.second);
        if(!pair.valid) continue;
        
        auto& neighbor_e = edges[pair.eidx];
        int nhidx = particle_to_hairs[neighbor_e.first];
        
        if(nhidx < 0) continue;
        
        const VectorXs& neta = flows[nhidx]->getEta();
        
        int ni_local_0 = particle_local_indices[neighbor_e.first];
        
        int ni_local_1 = particle_local_indices[neighbor_e.second];
        
        if(DIM == 2) {
          accum_solid += mathutils::lerp( radii(neighbor_e.first), radii(neighbor_e.second), pair.alpha );
          accum_liquid += mathutils::lerp( neta(ni_local_0) + radii(neighbor_e.first), neta(ni_local_1) + radii(neighbor_e.second), pair.alpha );
          count++;
        } else {
          accum_solid += M_PI * mathutils::lerp( radii(neighbor_e.first) * radii(neighbor_e.first), radii(neighbor_e.second) * radii(neighbor_e.second), pair.alpha );
          accum_liquid += M_PI * mathutils::lerp( (neta(ni_local_0) + radii(neighbor_e.first)) * (neta(ni_local_0) + radii(neighbor_e.first)),
                                                 (neta(ni_local_1) + radii(neighbor_e.second)) * (neta(ni_local_1) + radii(neighbor_e.second)), pair.alpha );
          count++;
        }
      }
      
      scalar a_eta = sqrt(accum_liquid / (scalar) count / M_PI);
      avg_eta(i) = a_eta;
      
      if(DIM == 2) {
        accum_liquid = std::min(accum_liquid, max_r);
      } else {
        accum_liquid = std::min(accum_liquid, M_PI * max_r * max_r);
      }
      
      
      scalar occupied = mathutils::clamp(accum_solid / accum_liquid, 0.0, 1.0);
      
      porosity(i) = 1.0 - occupied;
    }
  });
}


template<int DIM>
void PolygonalCohesion<DIM>::updatePorosityEEC()
{
  const std::vector<HairFlow<DIM>*>& flows = m_parent->getFilmFlows();
  const int nf = flows.size();
  const VectorXs& radii = m_parent->getRadii();
  const std::vector<int>& particle_to_hairs = m_parent->getParticleToHairs();
  const std::vector<int>& particle_local_indices = m_parent->getParticleToHairLocalIndices();
  
  threadutils::thread_pool::ParallelFor(0, nf, [&] (int hidx){
    HairFlow<DIM>* hair = flows[hidx];
    auto& particles = hair->getParticleIndices();
    VectorXs& eta = hair->getEta();
    VectorXs& porosity = hair->getPorosity();
    VectorXs& avg_eta = hair->getAvgEta();
    
    VectorXs weights(porosity.size());
    weights.setZero();
    
    porosity.setZero();
    
    int np = particles.size();
    
    for(int i = 0; i < np; ++i)
    {
      int pidx = particles[i];
      scalar accum_liquid = M_PI * eta(i) * eta(i);
      scalar accum_solid = M_PI * radii(pidx) * radii(pidx);
      
      scalar max_r = sqrt(accum_liquid / M_PI);
      int count = 1;
      
      //todo, update this based on ParticleParticlePairs...
      const std::vector<int>& adjs = m_particle_to_pppairs[pidx];
      
      for(int npidx : adjs){
        int nhidx = particle_to_hairs[npidx];
        
        if(nhidx < 0) continue;
        
        const VectorXs& neta = flows[nhidx]->getEta();
        
        int ni_local = particle_local_indices[npidx];
        
        accum_solid += M_PI * radii(npidx) * radii(npidx);
        
        accum_liquid += M_PI * (neta(ni_local) + radii(npidx)) * (neta(ni_local) + radii(npidx));
        
        count++;
      }
      
      scalar a_eta = sqrt(accum_liquid / (scalar) count / M_PI);
      avg_eta(i) = a_eta;
      
      accum_liquid = std::min(accum_liquid, M_PI * max_r * max_r);
      
      scalar occupied = mathutils::clamp(accum_solid / accum_liquid, 0.0, 1.0);
      
      porosity(i) = 1.0 - occupied;
    }
  });
}

template<int DIM>
void PolygonalCohesion<DIM>::addEnergyToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, scalar& E )
{
  const VectorXs& x0 = m_parent->getX();
  E += (m_gradE + m_hessE * (x - x0)).squaredNorm() * 0.5;
}

const scalar modulate_parameter = 1.0;

template<int DIM>
void PolygonalCohesion<DIM>::accumulateEFJPairwise(const VectorXs& x, const VectorXs& v, VectorXs& gradE, MatrixXs& hessE, MatrixXs& hessV)
{
}

template<int DIM>
void PolygonalCohesion<DIM>::addGradEToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, VectorXs& gradE )
{
}

template<int DIM>
void PolygonalCohesion<DIM>::setUseDecoupledForce(bool decoupled)
{
  m_use_decoupled_force = decoupled;
}

template<int DIM>
void PolygonalCohesion<DIM>::setUseParticlePOEMap(bool ppoemap)
{
  m_compute_particle_poe_mapping = ppoemap;
}

template<int DIM>
void PolygonalCohesion<DIM>::addHessXToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, TripletXs& hessE )
{
}

template<int DIM>
void PolygonalCohesion<DIM>::addHessVToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, TripletXs& hessE )
{
}

template<int DIM>
void PolygonalCohesion<DIM>::addGradEToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, VectorXs& gradE, int pidx )
{
}

template<int DIM>
void PolygonalCohesion<DIM>::addHessXToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, VectorXs& hessE, int pidx )
{
}

template<int DIM>
void PolygonalCohesion<DIM>::addHessVToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, VectorXs& hessE, int pidx )
{
}


template<int DIM>
void PolygonalCohesion<DIM>::updateViscousStartPhi( const VectorXs& x )
{
  if( m_parent->getParameter().damping_multiplier > 0. ){
    
    const std::vector< std::pair<int, int> >& edges = m_parent->getEdges();
    const int npep = (int) m_point_edge_pairs.size();
    m_viscous_start_phi.clear();
    
    const VectorXs& radii = m_parent->getRadii();
    for(int pair_idx = 0; pair_idx < npep; ++pair_idx) {
      
      Vectors<DIM> xc0;
      Vectors<DIM> xc1;
      
      const PointEdgePair& pair = m_point_edge_pairs[pair_idx];
      const std::pair<int, int>& base_e = edges[pair.base_eidx];
      const std::pair<int, int>& neighbor_e = edges[pair.neighbor_eidx];
      
      const Vectors<DIM> x0 = x.segment<DIM>( m_parent->getDof( base_e.first) );
      const Vectors<DIM> x1 = x.segment<DIM>( m_parent->getDof( base_e.second) );
      const Vectors<DIM> x2 = x.segment<DIM>( m_parent->getDof( neighbor_e.first) );
      const Vectors<DIM> x3 = x.segment<DIM>( m_parent->getDof( neighbor_e.second) );
      
      xc0 = x0 + pair.alpha_point * (x1 - x0);
      xc1 = x2 + pair.neighbor_alpha_point * (x3 - x2);
      scalar r0 = mathutils::lerp(radii( base_e.first ), radii(base_e.second), pair.alpha_point);
      
      scalar d0 = (xc1 - xc0).norm();
      scalar curr_phi = d0 - getDStar(r0);
      m_viscous_start_phi[ pair.hash_code ] = curr_phi;
      
    }
  }
}

template<int DIM>
void PolygonalCohesion<DIM>::computeIntegrationVars( const VectorXs& x, const VectorXs& v, const VectorXs& m,
                                                    VectorXs& lambda, VectorXs& lambda_v,
                                                    TripletXs& J, TripletXs& Jv, TripletXs& Jxv, TripletXs& tildeK,
                                                    TripletXs& stiffness, TripletXs& damping, VectorXs& Phi, VectorXs& Phiv, const scalar& dt)
{
  const std::vector< std::pair<int, int> >& edges = m_parent->getEdges();
  
  const int npep = (int) m_point_edge_pairs.size();
  
  const VectorXs& radii = m_parent->getRadii();
  
  const scalar& cohesion_multiplier = m_parent->getHairHairCohesionMultiplier();
  
  // in case it's missing
  if( m_parent->getParameter().damping_multiplier > 0. ){
    for(int pair_idx = 0; pair_idx < npep; ++pair_idx) {
      const PointEdgePair& pair = m_point_edge_pairs[pair_idx];
      
      std::unordered_map<uint64, scalar>::const_iterator poep = m_viscous_start_phi.find(pair.hash_code);
      if ( poep == m_viscous_start_phi.end() ){
        const std::pair<int, int>& base_e = edges[pair.base_eidx];
        const std::pair<int, int>& neighbor_e = edges[pair.neighbor_eidx];
        
        const Vectors<DIM> x0 = x.segment<DIM>( m_parent->getDof( base_e.first) ) + v.segment<DIM>( m_parent->getDof( base_e.first ) ) * dt * pair.time;
        const Vectors<DIM> x1 = x.segment<DIM>( m_parent->getDof( base_e.second) ) + v.segment<DIM>( m_parent->getDof( base_e.second ) ) * dt * pair.time;
        const Vectors<DIM> x2 = x.segment<DIM>( m_parent->getDof( neighbor_e.first) ) + v.segment<DIM>( m_parent->getDof( neighbor_e.first ) ) * dt * pair.time;
        const Vectors<DIM> x3 = x.segment<DIM>( m_parent->getDof( neighbor_e.second) ) + v.segment<DIM>( m_parent->getDof( neighbor_e.second ) ) * dt * pair.time;;
        
        Vectors<DIM> xc0;
        Vectors<DIM> xc1;
        
        xc0 = x0 + pair.alpha_point * (x1 - x0);
        xc1 = x2 + pair.neighbor_alpha_point * (x3 - x2);
        
        scalar r0 = mathutils::lerp( radii(base_e.first), radii(base_e.second), pair.alpha_point );
        
        scalar d0 = (xc1 - xc0).norm();
        scalar curr_phi = d0 - getDStar(r0);
        
        m_viscous_start_phi[ pair.hash_code ] = curr_phi;
      }
    }
  }
  

  threadutils::thread_pool::ParallelFor(0, npep, [&] (int pair_idx) {
    VectorXs dddx(DIM * 4);
    Vectors<DIM> xc0;
    Vectors<DIM> xc1;
    
    const PointEdgePair& pair = m_point_edge_pairs[pair_idx];
    const std::pair<int, int>& base_e = edges[pair.base_eidx];
    const std::pair<int, int>& neighbor_e = edges[pair.neighbor_eidx];
    
    const Vectors<DIM> x0 = x.segment<DIM>( m_parent->getDof( base_e.first) ) + v.segment<DIM>( m_parent->getDof( base_e.first ) ) * dt * pair.time;
    const Vectors<DIM> x1 = x.segment<DIM>( m_parent->getDof( base_e.second) ) + v.segment<DIM>( m_parent->getDof( base_e.second ) ) * dt * pair.time;
    const Vectors<DIM> x2 = x.segment<DIM>( m_parent->getDof( neighbor_e.first) ) + v.segment<DIM>( m_parent->getDof( neighbor_e.first ) ) * dt * pair.time;
    const Vectors<DIM> x3 = x.segment<DIM>( m_parent->getDof( neighbor_e.second) ) + v.segment<DIM>( m_parent->getDof( neighbor_e.second ) ) * dt * pair.time;;
    
    xc0 = x0 + pair.alpha_point * (x1 - x0);
    xc1 = x2 + pair.neighbor_alpha_point * (x3 - x2);
    scalar r0 = mathutils::lerp( radii(base_e.first), radii(base_e.second), pair.alpha_point );
    
    scalar d0 = (xc1 - xc0).norm();
    mathutils::grad_pointline_dist_fixed(x0, x1, x2, x3, pair.alpha_point, pair.neighbor_alpha_point, dddx);
    
    scalar k0 = getStiffness(r0, d0, pair.V, pair.pressure_weight * cohesion_multiplier) * pair.quadrature_weight;
    stiffness[m_internal_index_pos + pair_idx] = Triplets(m_internal_index_pos + pair_idx, m_internal_index_pos + pair_idx, k0 );
    // damping[m_internal_index_vel + pair_idx] = Triplets(m_internal_index_vel + pair_idx, m_internal_index_vel + pair_idx, k0 * dt);
    
    scalar curr_phi = d0 - getDStar(r0);
    Phi(m_internal_index_pos + pair_idx) = curr_phi;
    
    const scalar& damping_multiplier = m_parent->getParameter().damping_multiplier;
    
    if( damping_multiplier > 0. ){
      scalar viscous_phi = 0.;
      std::unordered_map<uint64, scalar>::const_iterator poep = m_viscous_start_phi.find(pair.hash_code);
      if ( poep != m_viscous_start_phi.end() ){
        viscous_phi = curr_phi - poep->second;
      }
      // std::cout << pair_idx << " " << viscous_phi << std::endl;
      stiffness[m_internal_index_pos + npep + pair_idx] = Triplets(m_internal_index_pos + npep + pair_idx, m_internal_index_pos + npep + pair_idx, k0 * damping_multiplier );
      Phi(m_internal_index_pos + npep + pair_idx) = viscous_phi;
    }
    
    // hess_pointline_dist_fixed(x0, x1, x2, x3, pair.alpha_point, pair.neighbor_alpha_point, d2ddx2);
    const int indices[] = {base_e.first, base_e.second, neighbor_e.first, neighbor_e.second};
    for(int i = 0; i < 4; ++i)
    {
      for(int r = 0; r < DIM; ++r)
      {
        J[m_internal_index_J + pair_idx * DIM * 4 + DIM * i + r] = Triplets(m_internal_index_pos + pair_idx, m_parent->getDof( indices[i] ) + r, dddx(DIM * i + r));
        if( damping_multiplier > 0. ){
          J[m_internal_index_J + (npep * DIM * 4) + pair_idx * DIM * 4 + DIM * i + r] = Triplets(m_internal_index_pos + npep + pair_idx, m_parent->getDof( indices[i] ) + r, dddx(DIM * i + r));
        }
        //Jv[m_internal_index_Jv + pair_idx * DIM * 4 + DIM * i + r] = Triplets(m_internal_index_vel + pair_idx, m_parent->getDof( indices[i] ) + r, dddx(DIM * i + r));
      }
    }
    
    /*    for(int i = 0; i < 4; ++i) {
     for(int j = 0; j < 4; ++j) {
     for(int r = 0; r < DIM; ++r) {
     for(int s = 0; s < DIM; ++s) {
     int idx_local = i * DIM + r;
     int jdx_local = j * DIM + s;
     int ijdx_local = jdx_local * DIM * 4 + idx_local;
     tildeK[m_internal_index_tildeK + pair_idx * (DIM * 4) * (DIM * 4) + ijdx_local] = Triplets( m_parent->getDof( indices[i] ) + r,  m_parent->getDof( indices[j] ) + s, d2ddx2(i * DIM + r, j * DIM + s) * (-ll));
     }
     }
     }
     }*/
  });
}

template<int DIM>
int PolygonalCohesion<DIM>::numJ()
{
  if( m_parent->getParameter().damping_multiplier > 0. ){
    return m_point_edge_pairs.size() * DIM * 4 * 2;
  }
  return m_point_edge_pairs.size() * DIM * 4;
}

template<int DIM>
int PolygonalCohesion<DIM>::numJv()
{
  return 0;//m_point_edge_pairs.size() * DIM * 4;
}

template<int DIM>
int PolygonalCohesion<DIM>::numJxv()
{
  return 0;
}

template<int DIM>
int PolygonalCohesion<DIM>::numTildeK()
{
  return 0; //m_point_edge_pairs.size() * (DIM * 4) * (DIM * 4);
}

template<int DIM>
bool PolygonalCohesion<DIM>::isParallelized()
{
  return true;
}

template<int DIM>
void PolygonalCohesion<DIM>::storeLambda(const VectorXs& lambda, const VectorXs& lambda_v)
{
}

template<int DIM>
int PolygonalCohesion<DIM>::numConstraintPos()
{
  if( m_parent->getParameter().damping_multiplier > 0. ){
    return (int) 2 * m_point_edge_pairs.size();
  }
  return (int) m_point_edge_pairs.size();
}

template<int DIM>
int PolygonalCohesion<DIM>::numConstraintVel()
{
  return 0;//(int) m_point_edge_pairs.size();
}

template<int DIM>
Force* PolygonalCohesion<DIM>::createNewCopy()
{
  return new PolygonalCohesion<DIM>(m_parent);
}

template<int DIM>
const char* PolygonalCohesion<DIM>::name()
{
  return "polygonalcohesion";
}

template<int DIM>
void PolygonalCohesion<DIM>::getAffectedVars( int pidx, std::unordered_set<int>& vars )
{
  if( m_parent->getComponent(pidx) == DIM ) return;
  if(!m_use_decoupled_force) {
    int ip = m_parent->getVertFromDof( pidx );
    vars.insert(pidx);
    
    const std::vector< std::pair<int, int> >& edges = m_parent->getEdges();
    
    std::vector<int>& pepindices = m_particle_to_point_edge_pairs[ip];
    
    for(int pepidx : pepindices) {
      PointEdgePair& pep = m_point_edge_pairs[pepidx];
      
      auto& e = edges[pep.neighbor_eidx];
      for(int r = 0; r < DIM; ++r)
      {
        vars.insert( m_parent->getDof( e.first ) + r);
        vars.insert( m_parent->getDof( e.second ) + r);
      }
    }
  } else {
    int ip = m_parent->getVertFromDof( pidx );
    for(int r = 0; r < DIM; ++r)
      vars.insert( m_parent->getDof( ip ) + r);
  }
}

template<int DIM>
bool PolygonalCohesion<DIM>::isContained( int pidx )
{
  if( m_parent->getComponent(pidx) == DIM ) return false;
  return true;
}

template<int DIM>
void PolygonalCohesion<DIM>::preCompute( const VectorXs& x, const VectorXs& v, const VectorXs& m, const scalar& dt )
{
  if(DIM == 2 || !m_parent->useCtcd() ) computeConnections(x);
  updatePointEdgePairs(x, v, dt); // assumes end of time position given for X
#ifndef USE_RAW_LINEAR
  if(m_use_decoupled_force)
#endif
    accumulateEFJPairwise(x, v, m_gradE, m_hessE, m_hessV);
}

template<int DIM>
void PolygonalCohesion<DIM>::postStepScene( const scalar& dt )
{
  if(m_parent->useCtcd()) {
    updatePorosityEEC();
  } else {
    updatePorosity();
  }
  
  m_point_edge_pairs.resize(0);
}

template<int DIM>
bool PolygonalCohesion<DIM>::isPrecomputationParallelized()
{
  return true;
}

template<int DIM>
void PolygonalCohesion<DIM>::updateStructure(const VectorXs& x)
{
  buildSearchTree(x);
}

template<int DIM>
bool PolygonalCohesion<DIM>::isInterHair() const
{
  return true;
}

template<int DIM>
void PolygonalCohesion<DIM>::writeReadable(std::ostream& oss_pepairs, std::ostream& oss_poepairs, std::ostream& oss_pppairs) const
{
  const VectorXs& x = m_parent->getX();
  
  auto& edges = m_parent->getEdges();
  
  int np = m_adjacency_categorized.size();
  for (int i = 0; i < np; ++i) {
    const std::unordered_map<int, ParticleEdgePair*> pairs = m_adjacency_categorized[i];
    for(auto itr = pairs.begin(); itr != pairs.end(); ++itr) {
      ParticleEdgePair* pep = itr->second;
      
      const Vectors<DIM>& x0 = x.segment<DIM>(m_parent->getDof( pep->pidx ));
      auto& e = edges[pep->eidx];
      
      const Vectors<DIM>& x1 = x.segment<DIM>(m_parent->getDof(e.first));
      const Vectors<DIM>& x2 = x.segment<DIM>(m_parent->getDof(e.second));
      Vectors<DIM> xc = mathutils::lerp(x1, x2, pep->alpha);
      
      oss_pepairs << x0.transpose() << " " << xc.transpose() << " " << pep->valid << std::endl;
    }
  }
  
  int npep = m_point_edge_pairs.size();
  for(int pair_idx = 0; pair_idx < npep; ++pair_idx)
  {
    const PointEdgePair& pair = m_point_edge_pairs[pair_idx];
    
    const std::pair<int, int>& base_e = edges[pair.base_eidx];
    const std::pair<int, int>& e = edges[pair.neighbor_eidx];
    
    const Vectors<DIM>& x0 = x.segment<DIM>( m_parent->getDof( base_e.first ) );
    const Vectors<DIM>& x1 = x.segment<DIM>( m_parent->getDof( base_e.second ) );
    const Vectors<DIM>& x2 = x.segment<DIM>( m_parent->getDof( e.first ) );
    const Vectors<DIM>& x3 = x.segment<DIM>( m_parent->getDof( e.second ) );
    
    Vectors<DIM> xc0 = x0 + pair.alpha_point * (x1 - x0);
    Vectors<DIM> xc1 = x2 + pair.neighbor_alpha_point * (x3 - x2);
    
    oss_poepairs << xc0.transpose() << " " << xc1.transpose() << std::endl;
  }
  
  int nppp = m_particle_particle_pairs.size();
  for(int pair_idx = 0; pair_idx < nppp; ++pair_idx)
  {
    const ParticleParticlePair& pair = m_particle_particle_pairs[pair_idx];
    
    const Vectors<DIM>& x0 = x.segment<DIM>(m_parent->getDof(pair.pidx[0]));
    const Vectors<DIM>& x1 = x.segment<DIM>(m_parent->getDof(pair.pidx[1]));
    
    oss_pppairs << x0.transpose() << " " << x1.transpose() << std::endl;
  }
}

template<int DIM>
const VectorXs& PolygonalCohesion<DIM>::getAreaVGlobal() const
{
  return m_area_v_global;
}

template<int DIM>
VectorXs& PolygonalCohesion<DIM>::getAreaVGlobal()
{
  return m_area_v_global;
}

template<int DIM>
const VectorXs& PolygonalCohesion<DIM>::getAreaEGlobal() const
{
  return m_area_e_global;
}

template<int DIM>
VectorXs& PolygonalCohesion<DIM>::getAreaEGlobal()
{
  return m_area_e_global;
}

template<int DIM>
const VectorXs& PolygonalCohesion<DIM>::getPressureVGlobal() const
{
  return m_pressure_v_global;
}

template<int DIM>
VectorXs& PolygonalCohesion<DIM>::getPressureVGlobal()
{
  return m_pressure_v_global;
}

template<int DIM>
const MatrixXs& PolygonalCohesion<DIM>::getRhsOffsetVGlobal() const
{
  return m_rhs_offset_v_global;
}

template<int DIM>
MatrixXs& PolygonalCohesion<DIM>::getRhsOffsetVGlobal()
{
  return m_rhs_offset_v_global;
}

template<int DIM>
const std::vector<int>& PolygonalCohesion<DIM>::getPPPCountV() const
{
  return m_pplink_count;
}

template<int DIM>
scalar PolygonalCohesion<DIM>::getDStarPlanar(const scalar& radius) const
{
  if(m_min_cohesion_table->getRadii() == m_max_cohesion_table->getRadii()) {
    return m_min_cohesion_table->getDStarPlanar();
  } else {
    if(radius <= m_min_cohesion_table->getRadii()) {
      return m_min_cohesion_table->getDStarPlanar();
    } else if(radius >= m_max_cohesion_table->getRadii()) {
      return m_max_cohesion_table->getDStarPlanar();
    } else {
      scalar frac = m_min_cohesion_table->getRadii() * (m_max_cohesion_table->getRadii() - radius) /
      (radius * (m_max_cohesion_table->getRadii() - m_min_cohesion_table->getRadii()));
      return mathutils::lerp(m_max_cohesion_table->getDStarPlanar(), m_min_cohesion_table->getDStarPlanar(), frac);
    }
  }
}

template<int DIM>
scalar PolygonalCohesion<DIM>::getStiffnessPlanar(const scalar& radius, const scalar& d0, const scalar& A_target, const scalar& pressure_weight) const
{
  if(m_min_cohesion_table->getRadii() == m_max_cohesion_table->getRadii()) {
    return m_min_cohesion_table->getStiffnessPlanar(d0, A_target, pressure_weight);
  } else {
    if(radius <= m_min_cohesion_table->getRadii()) {
      return m_min_cohesion_table->getStiffnessPlanar(d0, A_target, pressure_weight);
    } else if(radius >= m_max_cohesion_table->getRadii()) {
      return m_max_cohesion_table->getStiffnessPlanar(d0, A_target, pressure_weight);
    } else {
      scalar frac = m_min_cohesion_table->getRadii() * (m_max_cohesion_table->getRadii() - radius) /
      (radius * (m_max_cohesion_table->getRadii() - m_min_cohesion_table->getRadii()));
      return mathutils::lerp(m_max_cohesion_table->getStiffnessPlanar(d0, A_target, pressure_weight), m_min_cohesion_table->getStiffnessPlanar(d0, A_target, pressure_weight), frac);
    }
  }
}

template<int DIM>
scalar PolygonalCohesion<DIM>::getDStar(const scalar& radius) const
{
  if(m_min_cohesion_table->getRadii() == m_max_cohesion_table->getRadii()) {
    return m_min_cohesion_table->getDStar();
  } else {
    if(radius <= m_min_cohesion_table->getRadii()) {
      return m_min_cohesion_table->getDStar();
    } else if(radius >= m_max_cohesion_table->getRadii()) {
      return m_max_cohesion_table->getDStar();
    } else {
      scalar frac = m_min_cohesion_table->getRadii() * (m_max_cohesion_table->getRadii() - radius) /
      (radius * (m_max_cohesion_table->getRadii() - m_min_cohesion_table->getRadii()));
      return mathutils::lerp(m_max_cohesion_table->getDStar(), m_min_cohesion_table->getDStar(), frac);
    }
  }
}

template<int DIM>
scalar PolygonalCohesion<DIM>::getStiffness(const scalar& radius, const scalar& d0, const scalar& A_target, const scalar& pressure_weight) const
{
  if(m_min_cohesion_table->getRadii() == m_max_cohesion_table->getRadii()) {
    return m_min_cohesion_table->getStiffness(d0, A_target, pressure_weight);
  } else {
    if(radius <= m_min_cohesion_table->getRadii()) {
      return m_min_cohesion_table->getStiffness(d0, A_target, pressure_weight);
    } else if(radius >= m_max_cohesion_table->getRadii()) {
      return m_max_cohesion_table->getStiffness(d0, A_target, pressure_weight);
    } else {
      scalar frac = m_min_cohesion_table->getRadii() * (m_max_cohesion_table->getRadii() - radius) /
      (radius * (m_max_cohesion_table->getRadii() - m_min_cohesion_table->getRadii()));
      return mathutils::lerp(m_max_cohesion_table->getStiffness(d0, A_target, pressure_weight), m_min_cohesion_table->getStiffness(d0, A_target, pressure_weight), frac);
    }
  }
}

// explicit instantiations at bottom
template class PolygonalCohesion<2>;
template class PolygonalCohesion<3>;

std::ostream& operator<<(std::ostream &os, const EdgeEdgePair &info)
{
  os << "[" << info.toString() << "]";
  return os;
}

std::ostream& operator<<(std::ostream &os, const ParticleEdgePair &info)
{
  os << "[" << info.toString() << "]";
  return os;
}

std::ostream& operator<<(std::ostream &os, const PointEdgePair &info)
{
  os << "[" << info.toString() << "]";
  return os;
}

std::ostream& operator<<(std::ostream &os, const ParticleParticlePair &info)
{
  os << "[" << info.toString() << "]";
  return os;
}
