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

#include "LevelSetForce.h"
#include "TwoDimensionalDisplayController.h"
#include "StringUtilities.h"
#include "MathUtilities.h"
#include "TwoDScene.h"
#include "fluidsim2D.h"
#include "fluidsim3D.h"
#include "array2_utils.h"
#include "array3_utils.h"
#include "CohesionTableGen.h"

const static char* g_szLevelSetForceName = "Levelset Force";

template<int DIM>
LevelSetForce<DIM>::LevelSetForce( TwoDScene<DIM>* parent, FluidSim* fluidsim, int hidx )
: Force()
, m_hidx(hidx)
, m_fluidsim(fluidsim)
, m_parent(parent)
{
  const HairFlow<DIM>* flow = m_parent->getFilmFlows()[m_hidx];
  m_particle_ls_pairs.resize(flow->size());
  m_point_ls_pairs.resize(flow->getGlobalEdges().size() * m_num_quadrature);
}

template<int DIM>
LevelSetForce<DIM>::~LevelSetForce()
{}

template<>
void LevelSetForce<3>::preCompute( const VectorXs& x, const VectorXs& v, const VectorXs& m, const scalar& dt )
{
  const HairFlow<3>* flow = m_parent->getFilmFlows()[m_hidx];
  FluidSim3D* fluid3d = (FluidSim3D*) m_fluidsim;
  const scalar& theta = m_parent->getLiquidTheta();
  const VectorXs& radii = m_parent->getRadii();
  
  auto& flow_particles = flow->getParticleIndices();
  const VectorXs& eta = flow->getEta();
  int nfp = flow_particles.size();
  
  for(int j = 0; j < nfp; ++j)
  {
    int particle = flow_particles[j];
    
    // Compute the cohesion
    const scalar& peta = eta(j);
    
    scalar H = peta + radii(particle);
    scalar max_dist = (1.0 + 0.5 * theta) * sqrt(M_PI * (H * H - radii(particle) * radii(particle)));
    
    scalar fq = fluid3d->get_nodal_solid_phi(x.segment<3>( m_parent->getDof( particle ) ));
    
    m_particle_ls_pairs[j].dist = std::max(0.0, fq);
    m_particle_ls_pairs[j].max_dist = max_dist;
    m_particle_ls_pairs[j].valid = (m_parent->isFixed(particle)) ? false : (fq < max_dist);
    
    if(j > 0 && j < nfp - 1) {
      if(m_parent->isFixed(flow_particles[j + 1]) || m_parent->isFixed(flow_particles[j - 1])) {
        m_particle_ls_pairs[j].valid = false;
      }
    }
  }

  m_point_ls_pairs.resize(flow->getGlobalEdges().size() * m_num_quadrature);
  memset(&m_point_ls_pairs[0], 0, m_point_ls_pairs.size() * sizeof(PointLSPair));
  
  scalar bb = m_parent->getParameter().damping_multiplier_planar;
  scalar friction = m_parent->getParameter().friction_multiplier_planar;
  
  // for each edge
  int ne = flow->getGlobalEdges().size();
  const std::vector<std::pair<int, int> >& edges = flow->getGlobalEdges();
  const std::vector<std::pair<int, int> >& local_edges = flow->getLocalEdges();
  
  for(int i = 0; i < ne; ++i)
  {
    const std::pair<int, int>& e = edges[i];
    const std::pair<int, int>& local_e = local_edges[i];
    const ParticleLSPair& pair_0 = m_particle_ls_pairs[local_e.first];
    const ParticleLSPair& pair_1 = m_particle_ls_pairs[local_e.second];
    
    scalar alpha_0 = 0.0;
    scalar alpha_1 = 1.0;
    
    if(pair_0.valid && pair_1.valid) {
      // both valid, do nothing to the alphas
    } else if(!pair_0.valid && !pair_1.valid) {
      // both invalid, ignore this edge-ls pair
      continue;
    } else if(pair_0.valid) {
      // alpha_1 < 1.0
      const scalar& r0i = radii(e.first);
      const scalar& r1i = radii(e.second);
      scalar ddiv = (pair_0.max_dist - (pair_0.dist - r0i)) - (pair_1.max_dist - (pair_1.dist - r1i));
      if(ddiv == 0.0) {
        // parallel, do nothing
      } else {
        scalar min_alpha = 0.0;
        if(i == 0) {
          min_alpha = 0.1;
        }
        alpha_1 = mathutils::clamp((pair_0.max_dist - (pair_0.dist - r0i)) / ddiv, min_alpha, 1.0);
        if(alpha_1 == min_alpha) continue;
      }
    } else if(pair_1.valid) {
      const scalar& r0i = radii(e.first);
      const scalar& r1i = radii(e.second);
      scalar ddiv = (pair_0.max_dist - (pair_0.dist - r0i)) - (pair_1.max_dist - (pair_1.dist - r1i));
      if(ddiv == 0.0) { // parallel
        // do nothing
      } else {
        scalar max_alpha = 1.0;
        if(i == ne - 1){
          max_alpha = 0.9;
        }
        alpha_0 = mathutils::clamp((pair_0.max_dist - (pair_0.dist - r0i)) / ddiv, 0.0, max_alpha);
        if(alpha_0 == max_alpha) continue;
      }
    }
    
    // build point-ls pairs
    for(int k_gauss = 0; k_gauss < m_num_quadrature; ++k_gauss)
    {
      int local_i0 = local_e.first;
      int local_i1 = local_e.second;
      const scalar alpha_gauss = mathutils::lerp(alpha_0, alpha_1, mathutils::gauss_legendre_point[m_num_quadrature - 1][k_gauss]);
      const scalar H0 = eta(local_i0) + radii(e.first);
      const scalar H1 = eta(local_i1) + radii(e.second);
      
      const scalar V = M_PI * mathutils::lerp(H0 * H0 - radii(e.first) * radii(e.first), H1 * H1 - radii(e.second) * radii(e.second), alpha_gauss);
      
      const Vector3s& x0 = x.segment<3>( m_parent->getDof( e.first ) );
      const Vector3s& x1 = x.segment<3>( m_parent->getDof( e.second ) );
      
      Vector3s xc0 = mathutils::lerp(x0, x1, alpha_gauss);
      scalar r0 = mathutils::lerp(radii(e.first), radii(e.second), alpha_gauss);
      
      const scalar area_base_e = (x1 - x0).norm();
      
      PointLSPair& plp = m_point_ls_pairs[i * m_num_quadrature + k_gauss];
      plp.eidx = i;
      plp.alpha_point = alpha_gauss;
      plp.pressure_weight = 1.0 - fluid3d->get_clamped_particle_weight(xc0);
      plp.quadrature_weight = mathutils::gauss_legendre_weight[m_num_quadrature - 1][k_gauss] * (alpha_1 - alpha_0) * area_base_e;
      plp.V = V;
      plp.k_gauss = k_gauss;
      
      if( bb > 0. || friction > 0. ) {
        Vector3s fsv = fluid3d->get_solid_velocity(xc0);
        
        const Vector3s& xsa = x.segment<3>( m_parent->getDof( e.first ) ) - v.segment<3>( m_parent->getDof( e.first ) ) * dt;
        const Vector3s& xsb = x.segment<3>( m_parent->getDof( e.second ) ) - v.segment<3>( m_parent->getDof( e.second ) ) * dt;
        Vector3s xc = xsa + plp.alpha_point * (xsb - xsa) + fsv * dt;
        
        plp.x0 = xc;

        scalar d0 = fluid3d->get_nodal_solid_phi(xc);

        scalar d_star = m_parent->getDStarPlanar(r0);
        plp.viscous_phi = (d0 - d_star);
      }
    }
  }
  
  m_point_ls_pairs.erase(std::remove_if(m_point_ls_pairs.begin(), m_point_ls_pairs.end(), [] (const PointLSPair& p) {
    return p.quadrature_weight == 0.0;
  }), m_point_ls_pairs.end());
}

template<>
void LevelSetForce<2>::preCompute( const VectorXs& x, const VectorXs& v, const VectorXs& m, const scalar& dt )
{
  const HairFlow<2>* flow = m_parent->getFilmFlows()[m_hidx];
  FluidSim2D* fluid2d = (FluidSim2D*) m_fluidsim;
  const scalar& theta = m_parent->getLiquidTheta();
  const VectorXs& radii = m_parent->getRadii();
  
  auto& flow_particles = flow->getParticleIndices();
  const VectorXs& eta = flow->getEta();
  int nfp = flow_particles.size();
  
  for(int j = 0; j < nfp; ++j)
  {
    int particle = flow_particles[j];
    
    // Compute the cohesion
    const scalar& peta = eta(j);
    
    scalar H = peta + radii(particle);
    scalar max_dist = (1.0 + 0.5 * theta) * sqrt(M_PI * (H * H - radii(particle) * radii(particle)));
    
    scalar fq = fluid2d->get_nodal_solid_phi(x.segment<2>(particle * 2));
    
    m_particle_ls_pairs[j].dist = std::max(0.0, fq);
    m_particle_ls_pairs[j].max_dist = max_dist;
    m_particle_ls_pairs[j].valid = (m_parent->isFixed(particle)) ? false : (fq < max_dist);
    
    if(j > 0 && j < nfp - 1) {
      if(m_parent->isFixed(flow_particles[j + 1]) || m_parent->isFixed(flow_particles[j - 1])) {
        m_particle_ls_pairs[j].valid = false;
      }
    }
  }
  
  m_point_ls_pairs.resize(flow->getGlobalEdges().size() * m_num_quadrature);
  memset(&m_point_ls_pairs[0], 0, m_point_ls_pairs.size() * sizeof(PointLSPair));
  
  // for each edge
  int ne = flow->getGlobalEdges().size();
  const std::vector<std::pair<int, int> >& edges = flow->getGlobalEdges();
  const std::vector<std::pair<int, int> >& local_edges = flow->getLocalEdges();
  
  scalar bb = m_parent->getParameter().damping_multiplier_planar;
  scalar friction = m_parent->getParameter().friction_multiplier_planar;
  
  for(int i = 0; i < ne; ++i)
  {
    const std::pair<int, int>& e = edges[i];
    const std::pair<int, int>& local_e = local_edges[i];
    const ParticleLSPair& pair_0 = m_particle_ls_pairs[local_e.first];
    const ParticleLSPair& pair_1 = m_particle_ls_pairs[local_e.second];
    
    scalar alpha_0 = 0.0;
    scalar alpha_1 = 1.0;
    
    if(pair_0.valid && pair_1.valid) {
      // both valid, do nothing to the alphas
    } else if(!pair_0.valid && !pair_1.valid) {
      // both invalid, ignore this edge-ls pair
      continue;
    } else if(pair_0.valid) {
      // alpha_1 < 1.0
      const scalar& r0i = radii(e.first);
      const scalar& r1i = radii(e.second);
      scalar ddiv = (pair_0.max_dist - (pair_0.dist - r0i)) - (pair_1.max_dist - (pair_1.dist - r1i));
      if(ddiv == 0.0) {
        // parallel, do nothing
      } else {
        scalar min_alpha = 0.0;
        if(i == 0) {
          min_alpha = 0.1;
        }
        alpha_1 = mathutils::clamp((pair_0.max_dist - (pair_0.dist - r0i)) / ddiv, min_alpha, 1.0);
        if(min_alpha == alpha_1) continue;
      }
    } else if(pair_1.valid) {
      const scalar& r0i = radii(e.first);
      const scalar& r1i = radii(e.second);
      scalar ddiv = (pair_0.max_dist - (pair_0.dist - r0i)) - (pair_1.max_dist - (pair_1.dist - r1i));
      if(ddiv == 0.0) { // parallel
        // do nothing
      } else {
        scalar max_alpha = 1.0;
        if(i == ne - 1){
          max_alpha = 0.9;
        }
        alpha_0 = mathutils::clamp((pair_0.max_dist - (pair_0.dist - r0i)) / ddiv, 0.0, max_alpha);
        if(alpha_0 == max_alpha) continue;
      }
    }
    
    // build point-ls pairs
    for(int k_gauss = 0; k_gauss < m_num_quadrature; ++k_gauss)
    {
      int local_i0 = local_e.first;
      int local_i1 = local_e.second;
      const scalar alpha_gauss = mathutils::lerp(alpha_0, alpha_1, mathutils::gauss_legendre_point[m_num_quadrature - 1][k_gauss]);
      const scalar H0 = eta(local_i0) + radii(e.first);
      const scalar H1 = eta(local_i1) + radii(e.second);
      
      const scalar V = M_PI * mathutils::lerp(H0 * H0 - radii(e.first) * radii(e.first), H1 * H1 - radii(e.second) * radii(e.second), alpha_gauss);
      
      const Vector2s& x0 = x.segment<2>(e.first * 2);
      const Vector2s& x1 = x.segment<2>(e.second * 2);
      
      Vector2s xc0 = mathutils::lerp(x0, x1, alpha_gauss);
      const scalar area_base_e = (x1 - x0).norm();
      scalar r0 = mathutils::lerp(radii(e.first), radii(e.second), alpha_gauss);

      PointLSPair& plp = m_point_ls_pairs[i * m_num_quadrature + k_gauss];
      plp.eidx = i;
      plp.alpha_point = alpha_gauss;
      plp.pressure_weight = 1.0 - fluid2d->get_clamped_particle_weight(xc0);
      plp.quadrature_weight = mathutils::gauss_legendre_weight[m_num_quadrature - 1][k_gauss] * (alpha_1 - alpha_0) * area_base_e;
      plp.V = V;
      plp.k_gauss = k_gauss;
      
      if( bb > 0. || friction > 0. ) {
        Vector2s fsv = fluid2d->get_solid_velocity(xc0);
        
        const Vector2s& xsa = x.segment<2>(e.first * 2) - v.segment<2>(e.first * 2) * dt;
        const Vector2s& xsb = x.segment<2>(e.second * 2) - v.segment<2>(e.second * 2) * dt;
        Vector2s xc = xsa + plp.alpha_point * (xsb - xsa) + fsv * dt;
      
        plp.x0 = xc;

        scalar d0 = fluid2d->get_nodal_solid_phi(xc);

        scalar d_star = m_parent->getDStarPlanar(r0);
        plp.viscous_phi = (d0 - d_star);
      }
    }
  }
  
  m_point_ls_pairs.erase(std::remove_if(m_point_ls_pairs.begin(), m_point_ls_pairs.end(), [] (const PointLSPair& p) {
    return p.quadrature_weight == 0.0;
  }), m_point_ls_pairs.end());
}

template<>
void LevelSetForce<2>::addEnergyToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, scalar& E )
{
  assert( x.size() == v.size() );
  assert( x.size()%2 == 0 );
}

template<>
void LevelSetForce<3>::addEnergyToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, scalar& E )
{
  assert( x.size() == v.size() );
  // assert( x.size()%3 == 0 );
}

template<>
void LevelSetForce<2>::addGradEToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, VectorXs& gradE )
{
  assert( x.size() == v.size() );
  assert( x.size() == gradE.size() );
  assert( x.size()%2 == 0 );
}

template<>
void LevelSetForce<3>::addGradEToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, VectorXs& gradE )
{
  assert( x.size() == v.size() );
  assert( x.size() == gradE.size() );     
}

template<>
void LevelSetForce<2>::addHessXToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, TripletXs& hessE )
{
  assert( x.size() == v.size() );
  assert( x.size() == m.size() );
  assert( x.size()%2 == 0 );
}

template<>
void LevelSetForce<3>::addHessXToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, TripletXs& hessE )
{
  assert( x.size() == v.size() );
  assert( x.size() == m.size() );
}

template<>
void LevelSetForce<2>::addHessVToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, TripletXs& hessE )
{
  assert( x.size() == v.size() );
  assert( x.size() == m.size() );
  assert( x.size()%2 == 0 );
}

template<>
void LevelSetForce<3>::addHessVToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, TripletXs& hessE )
{
  assert( x.size() == v.size() );
  assert( x.size() == m.size() );
}

template<int DIM>
void LevelSetForce<DIM>::addGradEToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, VectorXs& gradE, int pidx )
{
  assert( x.size() == v.size() );
  assert( x.size() == gradE.size() );
}

template<int DIM>
void LevelSetForce<DIM>::addHessXToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, VectorXs& hessE, int pidx )
{
  assert( x.size() == v.size() );
  assert( x.size() == m.size() );
}

template<int DIM>
void LevelSetForce<DIM>::addHessVToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, VectorXs& hessE, int pidx )
{
  assert( x.size() == v.size() );
  assert( x.size() == m.size() );
}

template<int DIM>
void LevelSetForce<DIM>::getAffectedVars( int pidx, std::unordered_set<int>& vars )
{
  int ip = m_parent->getVertFromDof(pidx);
  for(int r = 0; r < DIM; ++r)
    vars.insert( m_parent->getDof(ip) + r);
}

template<int DIM>
int LevelSetForce<DIM>::getAffectedHair( const std::vector<int> particle_to_hairs )
{
  return m_hidx;
}

template<int DIM>
Force* LevelSetForce<DIM>::createNewCopy()
{
  return new LevelSetForce(*this);
}

template<int DIM>
const char* LevelSetForce<DIM>::name()
{
  return g_szLevelSetForceName;
}

template<int DIM>
bool LevelSetForce<DIM>::isContained( int pidx )
{
  if( m_parent->getComponent(pidx) < DIM ) return true;
  return false;
}

template<>
void LevelSetForce<2>::postStepScene(const scalar& dt )
{
  VectorXs& x = m_parent->getX();
  
  const HairFlow<2>* flow = m_parent->getFilmFlows()[m_hidx];
  FluidSim2D* fluid2d = (FluidSim2D*) m_fluidsim;
  const VectorXs& radii = m_parent->getRadii();
  
  auto& flow_particles = flow->getParticleIndices();
  int nfp = flow_particles.size();
  
  for(int j = 0; j < nfp; ++j)
  {
    if(m_parent->isFixed(flow_particles[j])) continue;
    
    Vector2s p = x.segment<2>( m_parent->getDof(flow_particles[j]) );
    
    scalar phi_value = fluid2d->get_nodal_solid_phi(p) - radii(flow_particles[j]) * 2.0;
    if(phi_value < 0) {
      Vector2s normal = fluid2d->get_nodal_solid_phi_gradient(p);
      normal.normalize();
      
      x.segment<2>( m_parent->getDof(flow_particles[j]) ) -= phi_value*normal;
    }
  }
}

template<>
void LevelSetForce<3>::postStepScene(const scalar& dt )
{
  VectorXs& x = m_parent->getX();
  
  const HairFlow<3>* flow = m_parent->getFilmFlows()[m_hidx];
  FluidSim3D* fluid3d = (FluidSim3D*) m_fluidsim;
  const VectorXs& radii = m_parent->getRadii();
  
  auto& flow_particles = flow->getParticleIndices();
  int nfp = flow_particles.size();
  
  for(int j = 0; j < nfp; ++j)
  {
    if(m_parent->isFixed(flow_particles[j])) continue;
    
    Vector3s p = x.segment<3>( m_parent->getDof(flow_particles[j]) );
    
    scalar phi_value = fluid3d->get_nodal_solid_phi(p) - radii(flow_particles[j]) * 2.0;
    if(phi_value < 0) {
      Vector3s normal = fluid3d->get_nodal_solid_phi_gradient(p);
      normal.normalize();
      
      x.segment<3>( m_parent->getDof(flow_particles[j]) ) -= phi_value*normal;
    }
  }
}

template<>
void LevelSetForce<2>::computeIntegrationVars( const VectorXs& x, const VectorXs& v, const VectorXs& m,
                                                VectorXs& lambda, VectorXs& lambda_v,
                                                TripletXs& J, TripletXs& Jv, TripletXs& Jxv, TripletXs& tildeK,
                                                TripletXs& stiffness, TripletXs& damping, VectorXs& Phi, VectorXs& Phiv, const scalar& dt)
{
  FluidSim2D* fluid2d = (FluidSim2D*) m_fluidsim;
  int nplp = m_point_ls_pairs.size();
  
  const HairFlow<2>* flow = m_parent->getFilmFlows()[m_hidx];
  const std::vector<std::pair<int, int> >& edges = flow->getGlobalEdges();
  const VectorXs& radii = m_parent->getRadii();
  const scalar& cohesion_multiplier = m_parent->getHairSolidCohesionMultiplier();
  
  scalar bb = m_parent->getParameter().damping_multiplier_planar;
  scalar friction = m_parent->getParameter().friction_multiplier_planar;
  
  for(int pair_idx = 0; pair_idx < nplp; ++pair_idx) {
    const PointLSPair& pair = m_point_ls_pairs[pair_idx];
    const std::pair<int, int>& e = edges[pair.eidx];

    const Vector2s& x0 = x.segment<2>(e.first * 2);
    const Vector2s& x1 = x.segment<2>(e.second * 2);

    Vector2s xc = x0 + pair.alpha_point * (x1 - x0);
    scalar d0 = fluid2d->get_nodal_solid_phi(xc);
    
    scalar r0 = mathutils::lerp(radii(e.first), radii(e.second), pair.alpha_point);
    
    Vector2s dddx = fluid2d->get_nodal_solid_phi_gradient(xc);
    dddx.normalize();
    
    scalar d_star = m_parent->getDStarPlanar(r0);
    
    int base = m_internal_index_pos + pair_idx;
    scalar curr_phi = (d0 - d_star);
    
    scalar k0 = m_parent->getStiffnessPlanar(r0, d0, pair.V, pair.pressure_weight * cohesion_multiplier) * pair.quadrature_weight;
    
    Phi(base) = curr_phi;
    stiffness[base] = Triplets(base, base, k0);
    base += nplp;
    
    if( bb > 0. ){
      Phi(base) = curr_phi - pair.viscous_phi;
      stiffness[base] = Triplets(base, base, k0 * bb);
      base += nplp;
    }
    
    if( friction > 0. ) {
      Phi(base) = sqrt((xc - pair.x0).transpose() * (Matrix2s::Identity() - dddx * dddx.transpose()) * (xc - pair.x0));
      if(Phi(base) < 1e-20)
        stiffness[base] = Triplets(base, base, 0.0);
      else
        stiffness[base] = Triplets(base, base, k0 * friction);
      base += nplp;
    }

    const int indices[] = {e.first, e.second};
    const scalar multipliers[] = {(1.0 - pair.alpha_point), pair.alpha_point};
    
    for(int i = 0; i < 2; ++i) {
      for(int r = 0; r < 2; ++r) {
        int base_J = m_internal_index_J + pair_idx * (2 * 2) + 2 * i + r;
        int base_idx = m_internal_index_pos + pair_idx;
        
        J[base_J] = Triplets(base_idx, indices[i] * 2 + r, dddx(r) * multipliers[i]);
        base_J += nplp * 2 * 2;
        base_idx += nplp;
        
        if( bb > 0. ){
          J[base_J] = Triplets(base_idx, indices[i] * 2 + r, dddx(r) * multipliers[i]);
          base_J += nplp * 2 * 2;
          base_idx += nplp;
        }
        
        if( friction > 0. ) {
          Vector2s dddx_tangent = (Matrix2s::Identity() - dddx * dddx.transpose()) * (xc - pair.x0);
          if(dddx_tangent.norm() < 1e-20) dddx_tangent = Vector2s(1, 0);
          else {
            dddx_tangent.normalize();
          }
          J[base_J] = Triplets(base_idx, indices[i] * 2 + r, dddx_tangent(r) * multipliers[i]);
          base_J += nplp * 2 * 2;
          base_idx += nplp;
        }
      }
    }
  }
}

template<>
void LevelSetForce<3>::computeIntegrationVars( const VectorXs& x, const VectorXs& v, const VectorXs& m,
                                              VectorXs& lambda, VectorXs& lambda_v,
                                              TripletXs& J, TripletXs& Jv, TripletXs& Jxv, TripletXs& tildeK,
                                              TripletXs& stiffness, TripletXs& damping, VectorXs& Phi, VectorXs& Phiv, const scalar& dt)
{
  FluidSim3D* fluid3d = (FluidSim3D*) m_fluidsim;
  int nplp = m_point_ls_pairs.size();
  
  const HairFlow<3>* flow = m_parent->getFilmFlows()[m_hidx];
  const std::vector<std::pair<int, int> >& edges = flow->getGlobalEdges();
  const VectorXs& radii = m_parent->getRadii();
  const scalar& cohesion_multiplier = m_parent->getHairSolidCohesionMultiplier();
  
  scalar bb = m_parent->getParameter().damping_multiplier_planar;
  scalar friction = m_parent->getParameter().friction_multiplier_planar;
  
  for(int pair_idx = 0; pair_idx < nplp; ++pair_idx) {
    const PointLSPair& pair = m_point_ls_pairs[pair_idx];
    const std::pair<int, int>& e = edges[pair.eidx];

    const Vector3s& x0 = x.segment<3>( m_parent->getDof( e.first ) );
    const Vector3s& x1 = x.segment<3>( m_parent->getDof( e.second ) );
    
    Vector3s xc = x0 + pair.alpha_point * (x1 - x0);
    scalar d0 = fluid3d->get_nodal_solid_phi(xc);
    
    scalar r0 = mathutils::lerp(radii(e.first), radii(e.second), pair.alpha_point);
    
    Vector3s dddx = fluid3d->get_nodal_solid_phi_gradient(xc);;
    dddx.normalize();
    
    scalar d_star = m_parent->getDStarPlanar(r0);
    
    int base = m_internal_index_pos + pair_idx;
    scalar curr_phi = d0 - d_star;
    
    scalar k0 = m_parent->getStiffnessPlanar(r0, d0, pair.V, pair.pressure_weight * cohesion_multiplier) * pair.quadrature_weight;
    
    Phi(base) = curr_phi;
    stiffness[base] = Triplets(base, base, k0);
    base += nplp;
    
    if( bb > 0. ){
      Phi(base) = curr_phi - pair.viscous_phi;
      stiffness[base] = Triplets(base, base, k0 * bb );
      base += nplp;
    }
    
    if( friction > 0. ) {
      // friction
      Phi(base) = sqrt((xc - pair.x0).transpose() * (Matrix3s::Identity() - dddx * dddx.transpose()) * (xc - pair.x0));
      if(Phi(base) < 1e-20) {
        stiffness[base] = Triplets(base, base, 0.0 );
      } else {
        stiffness[base] = Triplets(base, base, k0 * friction );
      }
      base += nplp;
    }

    const int indices[] = {e.first, e.second};
    const scalar multipliers[] = {(1.0 - pair.alpha_point), pair.alpha_point};
    
    for(int i = 0; i < 2; ++i) {
      for(int r = 0; r < 3; ++r) {
        int base_J = m_internal_index_J + pair_idx * (3 * 2) + 3 * i + r;
        int base_idx = m_internal_index_pos + pair_idx;
        
        J[base_J] = Triplets(base_idx, m_parent->getDof(indices[i]) + r, dddx(r) * multipliers[i]);
        base_J += nplp * 3 * 2;
        base_idx += nplp;
        
        if( bb > 0. ){
          J[base_J] = Triplets(base_idx, m_parent->getDof(indices[i]) + r, dddx(r) * multipliers[i] );
          base_J += nplp * 3 * 2;
          base_idx += nplp;
        }
        
        if( friction > 0. ) {
          Vector3s dddx_tangent = (Matrix3s::Identity() - dddx * dddx.transpose()) * (xc - pair.x0);
          if(dddx_tangent.norm() < 1e-20) dddx_tangent = Vector3s(1, 0, 0);
          else {
            dddx_tangent.normalize();
          }
          J[base_J] = Triplets(base_idx, m_parent->getDof(indices[i]) + r, dddx_tangent(r) * multipliers[i] );
          base_J += nplp * 3 * 2;
          base_idx += nplp;
        }
      }
    }
  }
}


template<int DIM>
int LevelSetForce<DIM>::numJ()
{
  int num_J = m_point_ls_pairs.size() * DIM * 2;
  scalar bb = m_parent->getParameter().damping_multiplier_planar;
  scalar friction = m_parent->getParameter().friction_multiplier_planar;
  
  if( bb > 0. ){
    num_J += m_point_ls_pairs.size() * DIM * 2;
  }
  
  if( friction > 0. ){
    num_J += m_point_ls_pairs.size() * DIM * 2;
  }
  
  return num_J;   // TODO: Implementation
}

template<int DIM>
int LevelSetForce<DIM>::numJv()
{
  return 0;//m_point_ls_pairs.size() * DIM * 2;   // TODO: Implementation
}

template<int DIM>
int LevelSetForce<DIM>::numJxv()
{
  return 0;   // TODO: Implementation
}

template<int DIM>
int LevelSetForce<DIM>::numTildeK()
{
  return 0;//m_point_ls_pairs.size() * (DIM * 2) * (DIM * 2);   // TODO: Implementation
}

template<int DIM>
bool LevelSetForce<DIM>::isParallelized()
{
  return false;   // TODO: Implementation
}

template<int DIM>
bool LevelSetForce<DIM>::isPrecomputationParallelized()
{
  return false;
}

template<int DIM>
void LevelSetForce<DIM>::storeLambda(const VectorXs& lambda, const VectorXs& lambda_v)
{
}

template<int DIM>
int LevelSetForce<DIM>::numConstraintPos()
{
  int num_cp = (int) m_point_ls_pairs.size();
  scalar bb = m_parent->getParameter().damping_multiplier_planar;
  scalar friction = m_parent->getParameter().friction_multiplier_planar;
  
  if( bb > 0. ){
    num_cp += (int) m_point_ls_pairs.size();
  }
  
  if( friction > 0.) {
    num_cp += (int) m_point_ls_pairs.size();
  }
  
  return num_cp;
}

template<int DIM>
int LevelSetForce<DIM>::numConstraintVel()
{
  return 0;//(int) m_point_ls_pairs.size();
}

// explicit instantiations at bottom
template class LevelSetForce<2>;
template class LevelSetForce<3>;

