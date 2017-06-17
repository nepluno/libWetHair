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

#include "fluidsim2D.h"
#include "MathUtilities.h"
#include "ThreadUtils.h"
#include "FluidDragForce.h"

#include "array2_utils.h"

#include "pcgsolver/sparse_matrix.h"
#include "pcgsolver/pcg_solver.h"

#include "TwoDScene.h"
#include "sorter.h"

#include <fstream>
#include <numeric>

//#define SURFACE_TENSION_ON_HAIRS

using namespace mathutils;

const static int MAX_BRIDGE_PER_EDGE = 64;

void extrapolate(Array2s& grid, Array2s& old_grid, const Array2s& grid_weight, const Array2s& grid_liquid_weight, Array2c& valid, Array2c old_valid, const Vector2i& offset);

FluidSim2D* g_fluid2d_iptr = NULL;

static void sorter_callback(int pidx, int& i, int& j, int& k)
{
  auto& particles = g_fluid2d_iptr->get_particles();
  
  auto& p = particles[pidx];
  
  int pi = (int)((p.x(0) - g_fluid2d_iptr->get_origin()(0)) / g_fluid2d_iptr->cellsize());
  int pj = (int)((p.x(1) - g_fluid2d_iptr->get_origin()(1)) / g_fluid2d_iptr->cellsize());
  
  i = max(0, min(g_fluid2d_iptr->get_ni()-1, pi));
  j = max(0, min(g_fluid2d_iptr->get_ni()-1, pj));
  k = 0;
}

FluidSim2D::~FluidSim2D()
{
  if(m_sorter) delete m_sorter;
}

FluidSim2D::FluidSim2D(const Vector2s& origin_,
                       scalar width,
                       int ni_, int nj_,
                       const scalar& dz_,
                       const std::vector< Boundary<2>* >& boundaries_,
                       const std::vector< SourceBoundary<2>* >& sources_,
                       TwoDScene<2>* scene)
: m_parent(scene)
{
  g_fluid2d_iptr = this;
  ryoichi_correction_counter = 0;
  origin = origin_;
  boundaries = boundaries_;
  sources = sources_;
  ni = ni_;
  nj = nj_;
  dx = width / (scalar) ni;
  dz = dz_;
  
  u.resize(ni+1,nj); temp_u.resize(ni+1,nj); u_weights.resize(ni+1,nj); u_weight_hair.resize(ni+1,nj);
  u_valid.resize(ni+1,nj); u_hair.resize(ni+1,nj); u_particle.resize(ni+1, nj); u_weight_particle.resize(ni+1,nj);
  u_drag.resize(ni+1, nj); u_weight_total.resize(ni+1, nj); u_pressure_grad.resize(ni+1,nj); u_solid.resize(ni+1,nj);
  v.resize(ni,nj+1); temp_v.resize(ni,nj+1); v_weights.resize(ni,nj+1); v_weight_hair.resize(ni,nj+1);
  v_valid.resize(ni,nj+1); v_hair.resize(ni,nj+1); v_particle.resize(ni, nj+1); v_weight_particle.resize(ni,nj+1);
  v_drag.resize(ni, nj+1); v_weight_total.resize(ni, nj+1); v_pressure_grad.resize(ni+1,nj); v_solid.resize(ni,nj+1);
  
  u.set_zero();
  v.set_zero();
  u_pressure_grad.set_zero();
  v_pressure_grad.set_zero();
  u_hair.set_zero();
  v_hair.set_zero();
  u_particle.set_zero();
  v_particle.set_zero();
  
  u_weight_hair.set_zero();
  u_weight_particle.set_zero();
  v_weight_hair.set_zero();
  v_weight_particle.set_zero();
  u_drag.set_zero();
  v_drag.set_zero();
  u_solid.set_zero();
  v_solid.set_zero();
  
  u_weight_total.set_zero();
  v_weight_total.set_zero();
  
  temp_u.set_zero();
  temp_v.set_zero();
  
  nodal_solid_phi.resize(ni+1,nj+1);
  valid.resize(ni+1, nj+1);
  old_valid.resize(ni+1, nj+1);
  liquid_phi.resize(ni,nj);
  liquid_phi.assign(3*dx);
  
  pressure.resize(ni * nj, 0.0);
 
  m_bridges.reserve(scene->getNumEdges() * MAX_BRIDGE_PER_EDGE);
  
  update_boundary();
  m_sorter = new Sorter(ni, nj, 1);
  particles.clear();
  

  int nflows = m_parent->getNumFlows();
  drag_forces.resize(nflows);
  for(int i = 0; i < nflows; ++i) {
    drag_forces[i] = new FluidDragForce<2>(*scene, i);
    m_parent->insertForce(drag_forces[i]);
  }

  m_pool_liquid_index_cache.resize(nj);
  m_pool_liquid_particle_cache.resize(nj);
  
  for(auto& s : sources)
  {
    s->sample(this);
  }
}

void FluidSim2D::advect_boundary(const scalar& dt)
{
  int nb = boundaries.size();
  threadutils::thread_pool::ParallelFor(0, nb, [&] (int i) {
    if(boundaries[i]->type == BT_UNION || boundaries[i]->type == BT_INTERSECT) return;
    boundaries[i]->advance(dt);
  });
  
  int ns = sources.size();
  threadutils::thread_pool::ParallelFor(0, ns, [&] (int i) {
    if(sources[i]->type == BT_UNION || sources[i]->type == BT_INTERSECT) return;
    sources[i]->advance(dt);
  });
}

//Initialize the grid-based signed distance field that dictates the position of the solid boundary
void FluidSim2D::update_boundary() {
  threadutils::thread_pool::ParallelFor(0, nj + 1, [&] (int j) {
    for(int i = 0; i < ni+1; ++i) {
      Vector2s pos(i*dx,j*dx);
      Vector2s vel;
      nodal_solid_phi(i,j) = compute_phi_vel(pos + origin, vel);
    }
  });
  
  threadutils::thread_pool::ParallelFor(0, nj, [&] (int j) {
    for(int i = 0; i < ni+1; ++i) {
      Vector2s pos(i*dx,(j+0.5)*dx);
      Vector2s vel;
      compute_phi_vel(pos + origin, vel);
      u_solid(i, j) = vel[0];
    }
  });
  
  threadutils::thread_pool::ParallelFor(0, nj + 1, [&] (int j) {
    for(int i = 0; i < ni; ++i) {
      Vector2s pos((i+0.5)*dx,j*dx);
      Vector2s vel;
      compute_phi_vel(pos + origin, vel);
      v_solid(i, j) = vel[1];
    }
  });
}

scalar FluidSim2D::computeOverallDivergence()
{
  scalar div_sum = 0;
  for(int j = 0; j < nj; ++j) for(int i = 0; i < ni; ++i) {
    if(liquid_phi(i,j) > 0) continue;
    scalar divi = (u_weights(i + 1, j) * u(i + 1, j) - u_weights(i, j) * u(i, j) + v_weights(i, j + 1) * v(i, j + 1) - v_weights(i, j) * v(i, j)) / dx;
    div_sum += divi * divi;
  }
  
  return sqrt(div_sum);
}

void FluidSim2D::transferLiquidToGridParticle(const scalar& dt)
{
  int np = particles.size();
  threadutils::thread_pool::ParallelFor(0, np, [&] (int i) {
    particles[i].fresh *= 0.99;
  });
  
  int nhairp = m_parent->getNumParticles();
  const VectorXs& x = m_parent->getX();
  const VectorXs& v = m_parent->getV();
  VectorXs& m = m_parent->getM();
  const VectorXs& rest_m = m_parent->getHairRestMass();
  const scalar& rho_liq = m_parent->getLiquidDensity();
  const scalar& gamma_liq = m_parent->getLiquidTension();
  
  m_sorter->sort(nhairp, [&] (int pidx, int& i, int& j, int& k) {
    const Vector2s& pos = x.segment<2>(pidx * 2);
    int pi = (int)((pos(0) - origin(0)) / dx);
    int pj = (int)((pos(1) - origin(1)) / dx);
    
    i = max(0, min(ni-1, pi));
    j = max(0, min(nj-1, pj));
    k = 0;
  });
  
  scalar release_radius = dx * default_radius_multiplier() * sqrt(3.0) * 0.5;
  
  scalar release_vol = dropvol(release_radius);
  
  const std::vector<int>& particle_to_hair = m_parent->getParticleToHairs();
  const std::vector<int>& global_to_local = m_parent->getParticleToHairLocalIndices();
  std::vector<HairFlow<2>*>& hairs = m_parent->getFilmFlows();
  
  const scalar epsilon = 1e-7 * dx / dt;
  
  const int nflows = hairs.size();
  // record liquid pool into cache
  if(m_pool_liquid_vol_cache.size() != nflows * 2) m_pool_liquid_vol_cache.resize(nflows * 2);
  m_pool_liquid_vol_cache.setZero();
  
  bool dripping_zero_end = m_parent->getParameter().drippingnear;
  bool dripping_far_end = m_parent->getParameter().drippingfar;
  bool dripping_middle = m_parent->getParameter().drippingmiddle;
  
  threadutils::thread_pool::ParallelFor(0, nflows, [&] (int i) {
    const HairFlow<2>* hair = hairs[i];
    const MatrixXs& actual_u_v = hair->getActualVertexVelocity();
    const MatrixXs& tan_v = hair->getTangentV();
    int iend_local = actual_u_v.rows() - 1;
    
    scalar u_v0 = actual_u_v.row(0).dot(tan_v.row(0));
    scalar u_vn = actual_u_v.row(iend_local).dot(tan_v.row(iend_local));
    
    int count_pool = 0;
    if(dripping_zero_end && u_v0 < -epsilon) ++count_pool;
    if(dripping_far_end && u_vn > epsilon) ++count_pool;
    
    if(!count_pool) return;
    scalar each_pool_liquid = hair->getPoolSize() / (scalar) count_pool;
    
    if(dripping_zero_end && u_v0 < -epsilon) m_pool_liquid_vol_cache(i * 2 + 0) = each_pool_liquid;
    if(dripping_far_end && u_vn > epsilon) m_pool_liquid_vol_cache(i * 2 + 1) = each_pool_liquid;
  });

  MASS_UPDATE_MODE mum = m_parent->getMassUpdateMode();
  const scalar& liquid_shell = m_parent->getLiquidShell();
  const scalar& dripping_radius_multiplier = m_parent->getDrippingRadiusMultiplier();
  
  threadutils::thread_pool::ParallelFor(0, nj, [&] (int j) {
    m_pool_liquid_index_cache[j].resize(0);
    m_pool_liquid_particle_cache[j].resize(0);
    m_regular_liquid_particle_cache[j].resize(0);
    
    for(int i = 0; i < ni; ++i)
    {
      scalar total_pool_liquid = 0.0;
      scalar total_regular_liquid = 0.0;
      Vector2s avg_vel = Vector2s::Zero();
      Vector2s avg_pos = Vector2s::Zero();
      Vector2s avg_dir = Vector2s::Zero();
      Vector2s avg_pdirpl = Vector2s::Zero();
      
      scalar equivalent_radius = 0.0;
      scalar avg_area_v = 0.0;
      Vector2s avg_accel_norm = Vector2s::Zero();
      
      std::vector<int> hair_indices;
      
      int count = 0;
      m_sorter->getCellAt(i, j, 0, [&] (int pidx) {
        int hidx = particle_to_hair[pidx];
        int local_idx = global_to_local[pidx];
        HairFlow<2>* hair = hairs[hidx];
        int nhair_particles = (int) (hair->getParticleIndices().size());
        const MatrixXs& actual_u_v = hair->getActualVertexVelocity();
        const MatrixXs& accel_v = hair->getAccelV();
        const MatrixXs& tangent_v = hair->getTangentV();
        const MatrixXs& tangent_e = hair->getTangentE();
        const VectorXs& eta_v = hair->getEta();
        const VectorXs& radii_v = hair->getRadiiV();
        const VectorXs& area_v = hair->getAreaV();
        const VectorXs& area_e = hair->getAreaE();
        
        if(local_idx == 0) {
          total_pool_liquid += m_pool_liquid_vol_cache(hidx * 2 + 0);
          if(m_pool_liquid_vol_cache(hidx * 2 + 0) > 0.0) {
            avg_pos += x.segment<2>( m_parent->getDof( pidx ) ) * m_pool_liquid_vol_cache(hidx * 2 + 0);
            avg_vel += (v.segment<2>( m_parent->getDof( pidx ) ) + actual_u_v.row(local_idx).transpose()) * m_pool_liquid_vol_cache(hidx * 2 + 0);
            avg_dir += tangent_v.row( local_idx ).transpose() * m_pool_liquid_vol_cache(hidx * 2 + 0);
            avg_pdirpl += ( tangent_v.row( local_idx + 1 ).transpose() - tangent_v.row( local_idx ).transpose() ) / area_e(0) * m_pool_liquid_vol_cache(hidx * 2 + 0);
          }
        } else if(local_idx == nhair_particles - 1) {
          total_pool_liquid += m_pool_liquid_vol_cache(hidx * 2 + 1);
          if(m_pool_liquid_vol_cache(hidx * 2 + 1) > 0.0)
          {
            avg_pos += x.segment<2>( m_parent->getDof( pidx ) ) * m_pool_liquid_vol_cache(hidx * 2 + 1);
            avg_vel += (v.segment<2>( m_parent->getDof( pidx ) ) + actual_u_v.row(local_idx).transpose()) * m_pool_liquid_vol_cache(hidx * 2 + 1);
            avg_dir += tangent_v.row( local_idx ).transpose() * m_pool_liquid_vol_cache(hidx * 2 + 1);
            avg_pdirpl += ( tangent_v.row( local_idx ).transpose() - tangent_v.row( local_idx - 1 ).transpose() ) / area_e(local_idx - 1) * m_pool_liquid_vol_cache(hidx * 2 + 0);
          } else if(dripping_middle) {
            // dripping from regular vertex
            // compute regular volume
            scalar reg_vol = M_PI * ((eta_v(local_idx) + radii_v(local_idx)) * (eta_v(local_idx) + radii_v(local_idx)) - radii_v(local_idx) * radii_v(local_idx)) * area_v(local_idx);
            total_regular_liquid += reg_vol;
            
            // update avg pos and vel
            avg_pos += x.segment<2>( m_parent->getDof( pidx ) ) * reg_vol;
            avg_vel += (v.segment<2>( m_parent->getDof( pidx ) ) + actual_u_v.row(local_idx).transpose()) * reg_vol;
            avg_dir += tangent_v.row( local_idx ).transpose() * reg_vol;
            avg_pdirpl += ( tangent_e.row( local_idx ).transpose() - tangent_e.row( local_idx - 1 ).transpose() ) / area_v(local_idx) * reg_vol;
            
            // update equivalent radius
            equivalent_radius += radii_v(local_idx) * radii_v(local_idx);
            
            // compute normal acceleration
            Vector2s tan_dir = tangent_v.row(local_idx).transpose().normalized();
            Matrix2s proj = Matrix2s::Identity() - tan_dir * tan_dir.transpose();
            
            // acceleration to normal dir
            Vector2s accel_vert = accel_v.row(local_idx).transpose();
            avg_accel_norm += proj * accel_vert * reg_vol;
            
            // average segment length
            avg_area_v += area_v(local_idx) * reg_vol;

            ++count;
          }
        }
      });
      
      scalar total_liquid = total_regular_liquid + total_pool_liquid;
      
      if(total_liquid == 0.0) continue;
      
      avg_pos /= total_liquid;
      avg_vel /= total_liquid;
      avg_dir /= total_liquid;
      avg_pdirpl /= total_liquid;
      
      if(total_pool_liquid > release_vol) {
        scalar actual_release_vol = total_pool_liquid;
        scalar radii_release = dropradius(actual_release_vol);
        
        m_pool_liquid_particle_cache[j].push_back(Particle<2>(avg_pos, avg_vel, radii_release, PT_LIQUID));
        
        // mark released liquid from hairs
        m_sorter->getCellAt(i, j, 0, [&] (int pidx) {
          int hidx = particle_to_hair[pidx];
          int local_idx = global_to_local[pidx];
          HairFlow<2>* hair = hairs[hidx];
          int nhair_particles = (int) (hair->getParticleIndices().size());
          
          if(local_idx == 0) {
            m_pool_liquid_index_cache[j].push_back(hidx * 2 + 0);
          }
          
          if(local_idx == nhair_particles - 1) {
            m_pool_liquid_index_cache[j].push_back(hidx * 2 + 1);
          }
        });
      } else if(dripping_middle && total_regular_liquid > 0.0) {
        // try release from regular vertices
        scalar accel_v_norm = (avg_accel_norm / total_regular_liquid).norm();
        scalar inv_kappa_sqr = gamma_liq / (rho_liq * accel_v_norm);
        equivalent_radius = sqrt(equivalent_radius);
        avg_area_v /= total_regular_liquid;
        
        // We follow Equ.(3) of [Lorenceau 2004] to determine the maximum radius of liquid attached to hair
        // Lorenceau, Élise, Christophe Clanet, and David Quéré. "Capturing drops with a thin fiber."
        // Journal of colloid and interface science 279.1 (2004): 192-197.
        scalar V_max = 4.0 * M_PI * equivalent_radius * inv_kappa_sqr * (dripping_radius_multiplier * dripping_radius_multiplier);
        if(total_regular_liquid > V_max + release_vol) {
          // really release particles
          scalar actual_release_vol = total_regular_liquid - V_max;
          int N_release = std::max(1, (int) floor(actual_release_vol / release_vol));
          scalar adjusted_release_vol = actual_release_vol / (scalar) N_release;
          
          scalar radii_release = dropradius(adjusted_release_vol);
          
          for(int p = 0; p < N_release; ++p) {
            scalar rand_dist = avg_area_v * ((scalar) rand() / (scalar) RAND_MAX * 2.0 - 1.0);
            Vector2s proj_dir = (avg_dir + avg_pdirpl * rand_dist).normalized();
            
            Vector2s rand_vec_projected = proj_dir * rand_dist;
            m_regular_liquid_particle_cache[j].push_back(Particle<2>(avg_pos + rand_vec_projected, avg_vel, radii_release, PT_LIQUID));
          }
          
          // the proportion to remove liquid
          scalar prop = (total_regular_liquid - actual_release_vol) / total_regular_liquid;
          
          // remove volume from hair
          m_sorter->getCellAt(i, j, 0, [&] (int pidx) {
            int hidx = particle_to_hair[pidx];
            int local_idx = global_to_local[pidx];
            HairFlow<2>* hair = hairs[hidx];
            VectorXs& eta_v = hair->getEta();
            const VectorXs& radii_v = hair->getRadiiV();
            const VectorXs& area_v = hair->getAreaV();
            
            scalar new_eta = sqrt(prop * eta_v(local_idx) * (eta_v(local_idx) + radii_v(local_idx) * 2.0) + radii_v(local_idx) * radii_v(local_idx)) - radii_v(local_idx);
            eta_v(local_idx) = new_eta;
            
            // update mass
            if(mum != MUM_NONE) {
              scalar new_H = std::max(liquid_shell * radii_v(local_idx), new_eta + radii_v(local_idx));
              scalar new_mass_liq = (new_H * new_H - radii_v(local_idx) * radii_v(local_idx)) * M_PI * area_v(local_idx) * rho_liq;
              
              int idof = m_parent->getDof(pidx);
              scalar new_total_mass = rest_m(idof) + new_mass_liq;
              m.segment<2>( idof ).setConstant(new_total_mass);
            }
          });
        }
      }
    }
  });
  
  for(int j = 0; j < nj; ++j)
  {
    // combine released particle into global particles
    particles.insert(particles.end(), m_pool_liquid_particle_cache[j].begin(), m_pool_liquid_particle_cache[j].end());
    
    // remove liquid from hair liquid pool
    const std::vector<int>& indices = m_pool_liquid_index_cache[j];
    for(int idx : indices)
    {
      int hidx = idx / 2;
      HairFlow<2>* hair = hairs[hidx];
      hair->getPoolSize() -= m_pool_liquid_vol_cache(idx);
    }
  }
}

void FluidSim2D::shareParticleWithHairs( VectorXs& x, scalar dt )
{
  m_sorter->sort(particles.size(), sorter_callback);
  
  std::vector<HairFlow<2>*>& hairs = m_parent->getFilmFlows();
  // build bridges
  int nhair = hairs.size();
  
  const scalar theta = m_parent->getLiquidTheta();
  const scalar maxetaprop = m_parent->getMaxLimitEtaProp();
  const scalar liq_rho = m_parent->getLiquidDensity();
  
  m_bridges.resize(0);
  
  int np = particles.size();
  for(int i = 0; i < np; ++i)
  {
    particles[i].bridges.resize(0);
  }
  m_hair_bridge_buffer.resize(nhair);
  
  VectorXs& v = m_parent->getV();
  VectorXs& m = m_parent->getM();
  
  threadutils::thread_pool::ParallelFor(0, nhair, [&] (int i) {
    HairFlow<2>* hair = hairs[i];
    auto& edges = hair->getGlobalEdges();
    const VectorXs& eta = hair->getEta();
    const VectorXs& area_e = hair->getAreaE();
    const VectorXs& radii_e = hair->getRadiiE();
    const VectorXs& radii_v = hair->getRadiiV();
    std::vector< std::vector<int> >& edge_bridges = hair->getEdgeBridges();
    edge_bridges.resize(edges.size());
    
    m_hair_bridge_buffer[i].resize(0);
    
    int ne = edge_bridges.size();
    for(int j = 1; j < ne - 1; ++j)
    {
      assert(j < (int) edge_bridges.size());
      
      auto& eb = edge_bridges[j];
      eb.resize(0);
      
      scalar H0 = eta(j) + radii_v(j);
      scalar H1 = eta(j + 1) + radii_v(j + 1);
      
      scalar Hj = sqrt((H0 * H0 + H1 * H1) * 0.5);
      
      if(Hj > radii_e(j) * maxetaprop) continue;
      
      auto& e = edges[j];
      // construct bounding box
      const Vector2s& p0 = x.segment<2>(e.first * 2);
      const Vector2s& p1 = x.segment<2>(e.second * 2);
      
      const scalar search_dist = eta(j) + dx / sqrt(2.0);
      
      scalar pmin_x = std::min(p0(0) - search_dist, p1(0) - search_dist);
      scalar pmin_y = std::min(p0(1) - search_dist, p1(1) - search_dist);
      scalar pmax_x = std::max(p0(0) + search_dist, p1(0) + search_dist);
      scalar pmax_y = std::max(p0(1) + search_dist, p1(1) + search_dist);
      
      int imin_x = std::max(0, std::min(ni - 1, (int)( ( pmin_x - origin(0)) / dx ) ) );
      int imin_y = std::max(0, std::min(nj - 1, (int)( ( pmin_y - origin(1)) / dx ) ) );
      int imax_x = std::max(0, std::min(ni - 1, (int)( ceil( pmax_x - origin(0)) / dx ) ) );
      int imax_y = std::max(0, std::min(nj - 1, (int)( ceil( pmax_y - origin(1)) / dx ) ) );
      
      scalar vol_hair = M_PI * (H0 * H0 + H1 * H1 - radii_v(j) * radii_v(j) - radii_v(j + 1) * radii_v(j + 1)) * 0.5 * area_e(j);
      
      for(int s = imin_y; s <= imax_y; ++s)
      {
        for(int r = imin_x; r <= imax_x; ++r)
        {
          m_sorter->getCellAt(r, s, 0, [&] (int pidx) {
            Particle<2>& p = particles[pidx];
            if(p.type != PT_LIQUID) return;
            
            scalar vol_particle = dropvol(p.radii);
            
            const scalar rupture_dist = std::max(p.radii + eta(j), (1.0 + 0.5 * theta) * pow(vol_particle + vol_hair, 1.0 / 3.0));
            
            scalar alpha;
            
            scalar dist = mathutils::pointedgedist(p.x, p0, p1, alpha);
            
            Vector2s pc = p0 * (1.0 - alpha) + p1 * alpha;
            
            Vector2s pp = (pc - origin) / dx - Vector2s(0.5, 0.5);
            
            scalar phi = interpolate_value(pp, liquid_phi);
            
            if(dist < rupture_dist && phi < 0.0) {
              HairParticleBridge<2> b;
              b.volume = 0.0;
              b.vel.setZero();
              b.alpha = alpha;
              b.pidx = pidx;
              b.eidx = j;
              
              m_hair_bridge_buffer[i].push_back(b);
            }
          });
        }
      }
    }
  });
  
  int kidx = 0;
  for(int i = 0; i < nhair; ++i)
  {
    HairFlow<2>* hair = hairs[i];
    std::vector< std::vector<int> >& edge_bridges = hair->getEdgeBridges();
    auto& buffer = m_hair_bridge_buffer[i];
    int nb = buffer.size();
    
    for(int j = 0; j < nb; ++j)
    {
      const HairParticleBridge<2>& b = buffer[j];
      particles[b.pidx].bridges.push_back(kidx + j);
      edge_bridges[b.eidx].push_back(kidx + j);
    }
    
    m_bridges.insert(m_bridges.end(), buffer.begin(), buffer.end());
    kidx += nb;
  }
  
  // distribute liquid to bridges
  threadutils::thread_pool::ParallelFor(0, np, [&] (int i) {
    auto& p = particles[i];
    if(p.type != PT_LIQUID) return;
    
    auto& bridge_indices = p.bridges;
    if(bridge_indices.size() == 0) return;
    
    scalar vol_particle = dropvol(p.radii);
    scalar vol_bridge = vol_particle / (scalar) bridge_indices.size();
    
    for(int bidx : bridge_indices)
    {
      assert(bidx >= 0 && bidx < (int) m_bridges.size());
      
      auto& bridge = m_bridges[bidx];
      bridge.volume = vol_bridge;
      bridge.vel = p.v;
    }
    
    p.radii = 0.0;
  });
  
  MASS_UPDATE_MODE mum = m_parent->getMassUpdateMode();
  const scalar& absorption_rate = m_parent->getAbsorptionRate();
  
  // hairs absorb liquid on bridges
  for(int i = 0; i < nhair; ++i) {
    HairFlow<2>* hair = hairs[i];
    auto& edges = hair->getGlobalEdges();
    VectorXs& eta = hair->getEta();
    const VectorXs& area_v = hair->getAreaVHair();
    const VectorXs& area_e = hair->getAreaE();
    std::vector< std::vector<int> >& edge_bridges = hair->getEdgeBridges();
    const VectorXs& radii_v = hair->getRadiiV();
    const std::vector<int>& particle_indices = hair->getParticleIndices();
    VectorXs& velocity_e = hair->getVelocity();
    
    int ne = edges.size();
    int np = eta.size();

    for(int j = 1; j < np - 1; ++j)
    {
      scalar max_H = radii_v(j) * maxetaprop;
      scalar H = eta(j) + radii_v(j);
      scalar modified_H = std::max(H, m_parent->getLiquidShell() * radii_v(j));
      
      scalar old_liq_vol = M_PI * (H * H - radii_v(j) * radii_v(j)) * area_v(j);
      scalar old_modified_vol = M_PI * (modified_H * modified_H - radii_v(j) * radii_v(j)) * area_v(j);
      scalar old_modified_liq_mass = old_modified_vol * liq_rho;
      
      scalar max_liq_vol = M_PI * (max_H * max_H - radii_v(j) * radii_v(j)) * area_v(j);
      
      scalar max_diff_vol = max_liq_vol - old_liq_vol;
      
      if(max_diff_vol <= 0.0) continue;
      
      scalar sum_liq_vol = 0.0;
      
      // a vertex links to one or two edges
      int local_eidxs[] = {j - 1, j}; //previous edge, next edge
      scalar sum_num_bridges = 0.0;
      
      for(int k = 0; k < 2; ++k) {
        int local_eidx = local_eidxs[k];
        if(local_eidx < 0 || local_eidx >= ne) continue;
        
        auto& bridges_indices = edge_bridges[local_eidx];
        for(int bidx : bridges_indices)
        {
          if(bidx < 0 || bidx >= (int) m_bridges.size()) continue;
          
          auto& bridge = m_bridges[bidx];
          
          if(bridge.volume <= 0.0) continue;
          
          if(k == 0) {
            sum_num_bridges += bridge.alpha;
          } else {
            sum_num_bridges += (1.0 - bridge.alpha);
          }
        }
      }
      
      if(sum_num_bridges == 0.0) continue;
      
      scalar avg_absorb_vol = max_diff_vol / sum_num_bridges;
      
      for(int k = 0; k < 2; ++k) {
        int local_eidx = local_eidxs[k];
        if(local_eidx < 0 || local_eidx >= ne) continue;
        
        auto& bridges_indices = edge_bridges[local_eidx];
        if(bridges_indices.size() == 0) continue;
        
        for(int bidx : bridges_indices)
        {
          if(bidx < 0 || bidx >= (int) m_bridges.size()) continue;
          
          auto& bridge = m_bridges[bidx];
          
          if(bridge.volume <= 0.0) continue;
          
          scalar actual_vol;
          if(k == 0) {
            actual_vol = bridge.volume * bridge.alpha;
          } else {
            actual_vol = bridge.volume * (1.0 - bridge.alpha);
          }
          actual_vol = std::min(actual_vol, avg_absorb_vol * absorption_rate);
          
          // accumulate volume on vertex
          sum_liq_vol += actual_vol;
          
          // decrease bridge volume
          bridge.volume = std::max(0.0, bridge.volume - actual_vol);
        }
      }
      
      // finished gathering, update values on hair
      // new liquid volume on vertex
      scalar new_liq_vol = old_liq_vol + sum_liq_vol;
      
      // hair vertex volume
      scalar hair_vol = M_PI * radii_v(j) * radii_v(j) * area_v(j);
      
      // hair vertex global idx
      int pidx = particle_indices[j];
      
      int idof = m_parent->getDof(pidx);
      
      // hair vertex density
      const scalar& hair_rho = m_parent->getHairDensity( pidx );
      
      // hair vertex mass
      const scalar hair_mass = hair_vol * hair_rho;
      
      const scalar old_modified_mass = hair_mass + old_modified_liq_mass;
      
      // update liquid height
      const scalar new_H = sqrt(new_liq_vol / M_PI / area_v(j) + radii_v(j) * radii_v(j));
      
      const scalar new_modified_H = std::max(new_H, m_parent->getLiquidShell() * radii_v(j));
      
      const scalar new_modified_liq_mass = M_PI * (new_modified_H * new_modified_H - radii_v(j) * radii_v(j)) * area_v(j) * liq_rho;
      
      const scalar new_modified_mass = hair_mass + new_modified_liq_mass;
      
      eta(j) = new_H - radii_v(j);
      
      // update hair velocity
      if((mum == MUM_DIRECT_DIV || mum == MUM_MOMENTUM) && !m_parent->isFixed(pidx)) {
        if( new_modified_mass > 0.0 ) v.segment<2>( idof ) *= old_modified_mass / new_modified_mass;
        else v.segment<2>( idof ).setZero();
      }
      
      if(mum != MUM_NONE)
        m.segment<2>( idof ).setConstant(new_modified_mass);
      
      // update hair vertex angular momentum
      const int local_neighbor_idx[] = {j - 1, j + 1};
      // tangent momentum
      for(int k = 0; k < 2; ++k)
      {
        int local_eidx = local_eidxs[k];
        if(local_eidx < 0 || local_eidx >= ne) continue;
        
        int local_nj = local_neighbor_idx[k];
        
        const scalar H_nj = eta(local_nj) + radii_v(local_nj);
        const scalar modified_H_nj = std::max(H_nj, m_parent->getLiquidShell() * radii_v(local_nj));
        
        const scalar old_liq_mass_edge = M_PI * std::max(0.0, new_modified_H * new_modified_H + modified_H_nj * modified_H_nj - radii_v(local_nj) * radii_v(local_nj) - radii_v(j) * radii_v(j)) * 0.5 * area_e(local_eidxs[k]) * liq_rho;
        const scalar new_liq_mass_edge = M_PI * std::max(0.0, new_modified_H * new_modified_H + modified_H_nj * modified_H_nj - radii_v(local_nj) * radii_v(local_nj) - radii_v(j) * radii_v(j)) * 0.5 * area_e(local_eidxs[k]) * liq_rho;
        
        // update tangent velocity
        if(mum == MUM_DIRECT_DIV || mum == MUM_MOMENTUM) {
          if(new_liq_mass_edge > 0.0) velocity_e(local_eidx) *= old_liq_mass_edge / new_liq_mass_edge;
          else velocity_e(local_eidx) = 0.0;
        }
      }
    }
  }
  
  // particles absorb remain liquid on bridges
  threadutils::thread_pool::ParallelFor(0, np, [&] (int i) {
    auto& p = particles[i];
    if(p.type != PT_LIQUID) return;
    
    auto& bridge_indices = p.bridges;
    if(bridge_indices.size() == 0) return;
    
    scalar sum_vol = 0.0;
    for(int bidx : bridge_indices)
    {
      assert(bidx >= 0 && bidx < (int) m_bridges.size());
      auto& bridge = m_bridges[bidx];
      sum_vol += std::max(0.0, bridge.volume);
    }
    
    p.radii = dropradius(sum_vol);
  });
  
  threadutils::thread_pool::ParallelFor(0, np, [&] (int i) {
    particles[i].bridges.resize(0);
  });
  
  particles.erase(std::remove_if(particles.begin(), particles.end(), [] (const Particle<2>& p) {
    return p.type == PT_LIQUID && p.radii <= 1e-7;
  }), particles.end());
  
  m_sorter->sort(particles.size(), sorter_callback);
}

void FluidSim2D::resample(Vector2s& p, Vector2s& u, Matrix2s& c)
{
  scalar wsum = 0.0;
  Vector2s save = u;
  Matrix2s csave = c;
  u = Vector2s::Zero();
  c = Matrix2s::Zero();
  
  const scalar rho = m_parent->getLiquidDensity();
  
  int ix = std::max(0, std::min(ni - 1, (int)((p(0) - origin(0)) / dx)));
  int iy = std::max(0, std::min(nj - 1, (int)((p(1) - origin(1)) / dx)));
  
  const scalar re = dx;
  
  m_sorter->getNeigboringParticles_cell(ix, iy, 0, -1, 1, -1, 1, 0, 0, [&] (int pidx) {
    Particle<2>& np = particles[pidx];
    Vector2s diff = np.x - p;
    scalar w = dropvol(np.radii) * rho * mathutils::linear_kernel(diff, re);
    u += w * np.v;
    c += w * np.c;
    wsum += w;
  });

  if(wsum) {
    u /= wsum;
    c /= wsum;
  } else {
    u = save;
    c = csave;
  }
}

void FluidSim2D::correct(scalar dt)
{
  int np = (int) particles.size();
  
  if(!np) return;
  
#ifdef USE_MERGE_PARTICLES
  
  const scalar standard_radius = dx / sqrt(2.0);
  const scalar maximal_radius = standard_radius * sqrt(6.0) / 2.0;
  const scalar maximal_vol = dropvol(maximal_radius);
  // merge particles
  for(Particle<2> p : particles)
  {
    p.deceased = false;
  }
  
  std::cout << "<merge particles: " << particles.size() << " -> ";
  // for each particle find its closest neighbor, ignore the ones has been deceased or with radius larger than standard
  threadutils::thread_pool::ParallelFor(0, nj, [&] (int j) {
    for(int i = 0; i < ni; ++i) {
      m_sorter->getCellAt(i, j, 0, [&] (int pidx_i) {
        if(particles[pidx_i].type != PT_LIQUID || particles[pidx_i].deceased || particles[pidx_i].pressure * dx < 10.0) return;
        scalar vol_i = dropvol(particles[pidx_i].radii);
        
        m_sorter->getCellAt(i, j, 0, [&] (int pidx_j) {
          if(particles[pidx_j].type != PT_LIQUID || pidx_j <= pidx_i || particles[pidx_j].deceased || particles[pidx_j].pressure * dx < 10.0) return;
          
          // for each pair, check their add-up volume
          scalar vol_j = dropvol(particles[pidx_j].radii);
          if(vol_i + vol_j > maximal_vol) return;
          
          scalar dist = (particles[pidx_j].x - particles[pidx_i].x).norm();
          if(dist > particles[pidx_i].radii + particles[pidx_j].radii) return;
          
          // combine if their distance < add-up radius and add-up volume < maximal volume
          particles[pidx_j].deceased = true;
          scalar combined_vol = vol_i + vol_j;
          scalar alpha = vol_j / combined_vol;
          particles[pidx_i].x = mathutils::lerp(particles[pidx_i].x, particles[pidx_j].x, alpha);
          particles[pidx_i].v = mathutils::lerp(particles[pidx_i].v, particles[pidx_j].v, alpha);
          particles[pidx_i].c = mathutils::lerp(particles[pidx_i].c, particles[pidx_j].c, alpha);
          particles[pidx_i].radii = dropradius(combined_vol);
        });
      });
    }
  });
  
  particles.erase(std::remove_if(particles.begin(), particles.end(), [] (const Particle<2>& p) {
    return p.deceased;
  }), particles.end());
  
  m_sorter->sort(particles.size(), sorter_callback);
  std::cout << particles.size() << ">\n";
#endif
  scalar coeff = sqrt(ni * nj);
  
  int ryoichi_correction_step = m_parent->getParameter().fluidcorrectionsteps;
  
  // Compute Pseudo Moved Point
  threadutils::thread_pool::ParallelFor(0, np, [&] (int n) {
    Particle<2> &p = particles[n];
    
    if(p.type != PT_LIQUID) return;
    
    if(n % ryoichi_correction_step != ryoichi_correction_counter % ryoichi_correction_step) {
      return;
    }
    
    Vector2s spring = Vector2s::Zero();

    int ix = std::max(0, std::min((int)((p.x(0) - origin(0)) / dx), ni));
    int iy = std::max(0, std::min((int)((p.x(1) - origin(1)) / dx), nj));
    
    m_sorter->getNeigboringParticles_cell(ix, iy, 0, -1, 1, -1, 1, 0, 0, [&] (int pidx) {
      const Particle<2>& np = particles[pidx];
      if(n != pidx) {
        scalar re = 1.5 * sqrt(p.radii * np.radii);
        scalar dist = (p.x - np.x).norm();
        scalar w = coeff * mathutils::smooth_kernel(dist * dist, re);
        if( dist > 1e-4 * re )
        {
          spring += w * (p.x - np.x) / dist * re;
        } else {
          spring(0) += re / dt * (rand() & 0xFF) / 255.0;
          spring(1) += re / dt * (rand() & 0xFF) / 255.0;
        }
      }
    });
    
    p.buf0 = p.x + dt * spring;
    
    Vector2s pp = (p.buf0 - origin)/dx;
    scalar phi_value = interpolate_value(pp, nodal_solid_phi);
    if(phi_value < 0) {
      Vector2s normal;
      interpolate_gradient(normal, pp, nodal_solid_phi);
      normal.normalize();
      p.buf0 -= phi_value*normal;
    }
  });
  
  
  // Update
  threadutils::thread_pool::ParallelFor(0, np, [&] (int n) {
    Particle<2>& p = particles[n];
    if(p.type != PT_LIQUID) return;
    
    if(n % ryoichi_correction_step != ryoichi_correction_counter % ryoichi_correction_step) {
      return;
    }
    
    p.x = p.buf0;
  });
  
  m_sorter->sort(particles.size(), sorter_callback);
  
  ryoichi_correction_counter++;
}

void FluidSim2D::add_drag(scalar dt)
{  // drag
  threadutils::thread_pool::ParallelFor(0, u.nj, [&] (int j){
    for(int i = 0; i < u.ni; ++i)
    {
      u_particle(i, j) += u_drag(i, j) * dt / m_parent->getLiquidDensity();
    }
  });
  
  threadutils::thread_pool::ParallelFor(0, v.nj, [&] (int j){
    for(int i = 0; i < v.ni; ++i)
    {
      v_particle(i, j) += v_drag(i, j) * dt / m_parent->getLiquidDensity();
    }
  });
}

void FluidSim2D::add_gravity(scalar dt)
{
  const Vector3s& gravity = m_parent->getSimpleGravity();
  
  // gravity
  threadutils::thread_pool::ParallelFor(0, u.nj, [&] (int j){
    for(int i = 0; i < u.ni; ++i)
    {
      u_particle(i, j) += gravity(0) * dt;
    }
  });
  
  threadutils::thread_pool::ParallelFor(0, v.nj, [&] (int j){
    for(int i = 0; i < v.ni; ++i)
    {
      v_particle(i, j) += gravity(1) * dt;
    }
  });
}

//For extrapolated points, replace the normal component
//of velocity with the object velocity (in this case zero).
void FluidSim2D::constrain_velocity() {
  temp_u = u;
  temp_v = v;
  
  //(At lower grid resolutions, the normal estimate from the signed
  //distance function is poor, so it doesn't work quite as well.
  //An exact normal would do better.)
  
  //constrain u
  threadutils::thread_pool::ParallelFor(0, u.nj, [&] (int j){
    for(int i = 0; i < u.ni; ++i) {
      if(u_weights(i,j) == 0) {
        //apply constraint
        Vector2s pos = Vector2s(i*dx, (j+0.5)*dx) + origin;
        Vector2s vel = get_velocity(pos);
        Vector2s normal(0,0);
        interpolate_gradient(normal, Vector2s(i, j + 0.5), nodal_solid_phi);
        normal.normalize();
        scalar perp_component = vel.dot(normal);
        vel -= perp_component*normal;
        Vector2s vel_sol = get_solid_velocity(pos);
        vel += vel_sol.dot(normal) * normal;
        temp_u(i,j) = vel[0];
      }
    }
  });
  
  //constrain v
  threadutils::thread_pool::ParallelFor(0, v.nj, [&] (int j){
    for(int i = 0; i < v.ni; ++i) {
      if(v_weights(i,j) == 0) {
        //apply constraint
        Vector2s pos = Vector2s((i+0.5)*dx, j*dx) + origin;
        Vector2s vel = get_velocity(pos);
        Vector2s normal(0,0);
        interpolate_gradient(normal, Vector2s(i + 0.5, j), nodal_solid_phi);
        normal.normalize();
        scalar perp_component = vel.dot(normal);
        vel -= perp_component*normal;
        Vector2s vel_sol = get_solid_velocity(pos);
        vel += vel_sol.dot(normal) * normal;
        temp_v(i,j) = vel[1];
      }
    }
  });
  
  //update
  u = temp_u;
  v = temp_v;
  
}

//Add a tracer particle for visualization
void FluidSim2D::add_particle(const Particle<2>& p) {
  particles.push_back(p);
}

//move the particles in the fluid
void FluidSim2D::advect_particles(scalar dt) {
  int np = particles.size();
  threadutils::thread_pool::ParallelFor(0, np, [&](int p){
    particles[p].x += particles[p].v * dt;
    Vector2s pp = (particles[p].x - origin)/dx;
    
    //Particles can still occasionally leave the domain due to truncation errors,
    //interpolation error, or large timesteps, so we project them back in for good measure.
    
    //Try commenting this section out to see the degree of accumulated error.
    scalar phi_value = interpolate_value(pp, nodal_solid_phi);
    if(phi_value < 0) {
      Vector2s normal;
      interpolate_gradient(normal, pp, nodal_solid_phi);
      normal.normalize();
      particles[p].x -= phi_value*normal;
    }
  });
  
  sort_particles();
  
  for(SourceBoundary<2>* s : sources)
  {
    if(!s->activated) continue;
    
    for(const Vector2s& p : s->detectors)
    {
      if(p(0) < origin(0) || p(0) > origin(0) + ((scalar) ni+1.) * dx
         || p(1) < origin(1) || p(1) > origin(1) + ((scalar) nj+1.) * dx) continue;
      
      Vector2s pp = (p - origin)/dx;
      
      int ix = (int) pp(0);
      int iy = (int) pp(1);
      
      int num_p_need = 2;
      m_sorter->getCellAt(ix, iy, 0, [&] (int i) {
        if(num_p_need <= 0) return;
        num_p_need--;
      });
      
      for(int k = 0; k < num_p_need; ++k)
      {
        scalar x = ((scalar) ix + 0.5 + (((scalar)rand() / (scalar)RAND_MAX) * 0.5 - 0.5) ) * dx;
        scalar y = ((scalar) iy + 0.5 + (((scalar)rand() / (scalar)RAND_MAX) * 0.5 - 0.5) ) * dx;
        Vector2s pt = Vector2s(x, y) + origin;
        
        scalar phi_value = interpolate_value(pp, nodal_solid_phi);
        if(phi_value > 0.0) {
          Vector2s vel = Vector2s::Zero();
          s->compute_phi_vel(pt, vel);
          vel += s->eject_vel;
          add_particle(Particle<2>(pt, vel, dx / sqrt(2.0), PT_LIQUID));
        }
      }
    }
  }
  
  particles.erase( std::remove_if(particles.begin(), particles.end(), [&] (const Particle<2>& p) {
    return
    p.x(0) < origin(0) || p.x(0) > origin(0) + ((scalar) ni+1.) * dx
    || p.x(1) < origin(1) || p.x(1) > origin(1) + ((scalar) nj+1.) * dx;
  }), particles.end());
  
  sort_particles();
  
  
}


void FluidSim2D::constrain_hair_particles()
{
  int np = particles.size();
  
  const VectorXs& x = m_parent->getX();
  const std::vector< std::pair<int, int> >& edges = m_parent->getEdges();
  const std::vector< HairFlow<2>* >& flows = m_parent->getFilmFlows();
  const std::vector< int >& particle_hairs = m_parent->getParticleToHairs();
  const std::vector< int >& local_indices = m_parent->getParticleToHairLocalIndices();
  
  threadutils::thread_pool::ParallelFor(0, np, [&](int p){
    if(particles[p].type == PT_HAIR) {
      auto& e = edges[particles[p].edge_idx];
      const HairFlow<2>* flow = flows[ particle_hairs[e.first] ];
      const VectorXs& eta = flow->getEta();
      const VectorXs& radii_v = flow->getRadiiV();
      
      int local_idx0 = local_indices[e.first];
      int local_idx1 = local_indices[e.second];
      
      const Vector2s& x0 = x.segment<2>(m_parent->getDof(e.first));
      const Vector2s& x1 = x.segment<2>(m_parent->getDof(e.second));
      particles[p].x = mathutils::lerp(x0, x1, particles[p].edge_alpha);
      
      scalar H0 = eta(local_idx0) + radii_v(local_idx0);
      scalar H1 = eta(local_idx1) + radii_v(local_idx1);
      
      particles[p].radii = sqrt(mathutils::lerp(H0 * H0, H1 * H1, particles[p].edge_alpha));
    }
  });
  
  m_sorter->sort(particles.size(), sorter_callback);
}


void FluidSim2D::compute_liquid_phi()
{
  liquid_phi.assign(3*dx);
  std::cout << "[Liquid-Phi: CVT]" << std::endl;
  
  threadutils::thread_pool::ParallelFor(0, nj, [&] (int j){
    for(int i = 0; i < ni; ++i) {
      Vector2s pos = Vector2s((i+0.5)*dx, (j+0.5)*dx) + origin;
      scalar phi_min = 1e+20;
      
      m_sorter->getNeigboringParticles_cell(i, j, 0, -2, 2, -2, 2, 0, 0, [&] (int pidx) {
        const Particle<2>& p = particles[pidx];
        phi_min = std::min(phi_min, (pos - p.x).norm() - 1.02 * std::max(dx / sqrt(3.0), p.radii));
      });
      liquid_phi(i,j) = std::min(liquid_phi(i,j), phi_min);
    }
  });
  
  std::cout << "[Liquid-Phi: extrapolate phi into solids]" << std::endl;
  threadutils::thread_pool::ParallelFor(0, nj, [&](int j){
    for(int i = 0; i < ni; ++i) {
      if(liquid_phi(i,j) < 0.5*dx) {
        float solid_phi_val = 0.25 * (
                                      nodal_solid_phi(i,j) +
                                      nodal_solid_phi(i+1,j) +
                                      nodal_solid_phi(i,j+1) +
                                      nodal_solid_phi(i+1,j+1));
        if(solid_phi_val < 0)
          liquid_phi(i,j) = -0.5*dx;
      }
    }
  });
  
  //write_matlab_array(std::cout, liquid_phi, "phi");
}

void FluidSim2D::combine_velocity_field()
{
  threadutils::thread_pool::ParallelFor(0, u.nj, [&](int j){
    for(int i = 0; i < u.ni; ++i) {
      scalar weight = u_weight_particle(i, j) + u_weight_hair(i, j);
      if(weight > 0) {
        u(i, j) = (u_weight_particle(i, j) * u_particle(i, j) + u_weight_hair(i, j) * u_hair(i, j)) / weight;
      } else {
        u(i, j) = 0.0;
      }
      u_weight_total(i, j) = weight;
    }
  });
  
  threadutils::thread_pool::ParallelFor(0, v.nj, [&](int j){
    for(int i = 0; i < v.ni; ++i) {
      scalar weight = v_weight_particle(i, j) + v_weight_hair(i, j);
      if(weight > 0) {
        v(i, j) = (v_weight_particle(i, j) * v_particle(i, j) + v_weight_hair(i, j) * v_hair(i, j)) / weight;
      } else {
        v(i, j) = 0.0;
      }
      v_weight_total(i, j) = weight;
    }
  });
}

void FluidSim2D::project(scalar dt) {
  compute_weights();
  
  //Set up and solve the variational pressure solve.
  
  solve_pressure(dt);
  
  temp_u = u;
  temp_v = v;
  
  old_valid = valid;
  
  extrapolate(u, temp_u, u_weights, liquid_phi, valid, old_valid, Vector2i(-1, 0));
  extrapolate(v, temp_v, v_weights, liquid_phi, valid, old_valid, Vector2i(0, -1));
}

scalar FluidSim2D::get_nodal_solid_phi(const Vector2s& position) const
{
  Vector2s p = (position - origin) / dx;
  return interpolate_value(p, nodal_solid_phi);
}

Vector2s FluidSim2D::get_nodal_solid_phi_gradient(const Vector2s& position) const
{
  Vector2s p = (position - origin) / dx;
  Vector2s normal;
  interpolate_gradient(normal, p, nodal_solid_phi);
  return normal;
}

Vector2s FluidSim2D::get_velocity(const Vector2s& position, const Array2s& u_, const Array2s& v_) const
{
  //Interpolate the velocity from the u and v grids
  Vector2s p = (position - origin) / dx;
  Vector2s p0 = p - Vector2s(0, 0.5);
  Vector2s p1 = p - Vector2s(0.5, 0);
  scalar u_value = interpolate_value(p0, u_);
  scalar v_value = interpolate_value(p1, v_);
  
  return Vector2s(u_value, v_value);
}

Vector2s FluidSim2D::get_pressure_gradient(const Vector2s& position) const
{
  return get_velocity(position, u_pressure_grad, v_pressure_grad);
}

scalar FluidSim2D::get_pressure(const Vector2s& position) const
{
  Vector2s p = (position - origin) / dx - Vector2s(0.5, 0.5);
  return interpolate_value(p, pressure, ni, nj);
}

//Interpolate velocity from the MAC grid.
Vector2s FluidSim2D::get_velocity(const Vector2s& position) const {
  return get_velocity(position, u, v);
}

Vector2s FluidSim2D::get_solid_velocity(const Vector2s& position) const {
  return get_velocity(position, u_solid, v_solid);
}

//Interpolate drag from the MAC grid.
Vector2s FluidSim2D::get_particle_drag(const Vector2s& position) const {
  return get_velocity(position, u_drag, v_drag);
}

//Interpolate hair velocity from the MAC grid.
Vector2s FluidSim2D::get_hair_velocity(const Vector2s& position) const {
  return get_velocity(position, u_hair, v_hair);
}

//Interpolate particle velocity from the MAC grid.
Vector2s FluidSim2D::get_particle_velocity(const Vector2s& position) const {
  return get_velocity(position, u_particle, v_particle);
}

scalar FluidSim2D::getLiquidPhiValue(const Vector2s& position) const
{
  Vector2s pp = (position - origin) / dx - Vector2s(0.5, 0.5);
  return interpolate_value(pp, liquid_phi);
}

scalar FluidSim2D::getClampedLiquidPhiValue(const Vector2s& position) const
{
  const scalar standard_radius = dx * default_radius_multiplier();
  scalar w = getLiquidPhiValue(position);
  scalar criterion = -m_parent->getBulkThresholdMultiplier() * standard_radius;
  return mathutils::clamp(w / criterion, 0.0, 1.0);
}

Matrix2s FluidSim2D::get_affine_matrix(const Vector2s& position) const
{
  Vector2s p = (position - origin) / dx;
  Vector2s p0 = p - Vector2s(0, 0.5);
  Vector2s p1 = p - Vector2s(0.5, 0);
  
  Matrix2s c;
  c.col(0) = affine_interpolate_value(p0, u) / dx;
  c.col(1) = affine_interpolate_value(p1, v) / dx;
  
  return c;
}

//Compute finite-volume style face-weights for fluid from nodal signed distances
void FluidSim2D::compute_weights() {
  threadutils::thread_pool::ParallelFor(0, u_weights.nj, [&] (int j) {
    for(int i = 0; i < u_weights.ni; ++i) {
      u_weights(i,j) = 1 - mathutils::fraction_inside(nodal_solid_phi(i,j+1), nodal_solid_phi(i,j));
      u_weights(i,j) = hardclamp(u_weights(i,j), 0.0, 1.0);
    }
  });
  threadutils::thread_pool::ParallelFor(0, v_weights.nj, [&] (int j) {
    for(int i = 0; i < v_weights.ni; ++i) {
      v_weights(i,j) = 1 - mathutils::fraction_inside(nodal_solid_phi(i+1,j), nodal_solid_phi(i,j));
      v_weights(i,j) = hardclamp(v_weights(i,j), 0.0, 1.0);
    }
  });
  
}

//An implementation of the variational pressure projection solve for static geometry
void FluidSim2D::solve_pressure(scalar dt) {
  //This linear system could be simplified, but I've left it as is for clarity
  //and consistency with the standard naive discretization
  
  int ni = v.ni;
  int nj = u.nj;
  int system_size = ni*nj;
  if(rhs.size() != system_size) {
    rhs.resize(system_size);
    pressure.resize(system_size);
    matrix.resize(system_size);
  }
  matrix.zero();
  
  for(int i = 0; i < rhs.size(); ++i)
  {
    rhs[i] = 0.0;
  }
  
  const scalar rho = m_parent->getLiquidDensity();
  
  //Build the linear system for pressure
  for(int j = 1; j < nj-1; ++j) {
    for(int i = 1; i < ni-1; ++i) {
      int index = i + ni*j;
      rhs[index] = 0;
      pressure[index] = 0;
      scalar centre_phi = liquid_phi(i,j);
      if(centre_phi < 0) {
        
        //right neighbour
        scalar term = u_weights(i+1,j) * dt / sqr(dx) / rho;
        scalar right_phi = liquid_phi(i+1,j);
        if(right_phi < 0) {
          matrix.add_to_element(index, index, term);
          matrix.add_to_element(index, index + 1, -term);
        }
        else {
          scalar theta = mathutils::fraction_inside(centre_phi, right_phi);
          if(theta < 0.01) theta = 0.01;
          matrix.add_to_element(index, index, term/theta);
        }
        rhs[index] -= u_weights(i+1,j)*u(i+1,j) / dx;
        
        //left neighbour
        term = u_weights(i,j) * dt / sqr(dx) / rho;
        scalar left_phi = liquid_phi(i-1,j);
        if(left_phi < 0) {
          matrix.add_to_element(index, index, term);
          matrix.add_to_element(index, index - 1, -term);
        }
        else {
          scalar theta = mathutils::fraction_inside(centre_phi, left_phi);
          if(theta < 0.01) theta = 0.01;
          matrix.add_to_element(index, index, term/theta);
        }
        rhs[index] += u_weights(i,j)*u(i,j) / dx;
        
        //top neighbour
        term = v_weights(i,j+1) * dt / sqr(dx) / rho;
        scalar top_phi = liquid_phi(i,j+1);
        if(top_phi < 0) {
          matrix.add_to_element(index, index, term);
          matrix.add_to_element(index, index + ni, -term);
        }
        else {
          scalar theta = mathutils::fraction_inside(centre_phi, top_phi);
          if(theta < 0.01) theta = 0.01;
          matrix.add_to_element(index, index, term/theta);
        }
        rhs[index] -= v_weights(i,j+1)*v(i,j+1) / dx;
        
        //bottom neighbour
        term = v_weights(i,j) * dt / sqr(dx) / rho;
        scalar bot_phi = liquid_phi(i,j-1);
        if(bot_phi < 0) {
          matrix.add_to_element(index, index, term);
          matrix.add_to_element(index, index - ni, -term);
        }
        else {
          scalar theta = mathutils::fraction_inside(centre_phi, bot_phi);
          if(theta < 0.01) theta = 0.01;
          matrix.add_to_element(index, index, term/theta);
        }
        rhs[index] += v_weights(i,j)*v(i,j) / dx;
      }
    }
  }
  
  //Solve the system using Robert Bridson's incomplete Cholesky PCG solver
  
  scalar tolerance;
  int iterations;
  bool success = solver.solve(matrix, rhs, pressure, tolerance, iterations);
  if(!success) {
    printf("WARNING: Pressure solve failed!************************************************\n");
    
    std::cout << "rhs=[";
    for(scalar s : rhs) {
      std::cout << s << "; ";
    }
    std::cout << "];" << std::endl;
    
    std::cout << "pressure=[";
    for(scalar s : pressure) {
      std::cout << s << "; ";
    }
    std::cout << "];" << std::endl;
    
    write_matlab_array(std::cout, u, "u");
    write_matlab_array(std::cout, v, "v");

    write_matlab_array(std::cout, u_hair, "u_hair");
    write_matlab_array(std::cout, v_hair, "v_hair");
    
    write_matlab_array(std::cout, u_particle, "u_particle");
    write_matlab_array(std::cout, v_particle, "v_particle");
    
    write_matlab_array(std::cout, u_drag, "u_drag");
    write_matlab_array(std::cout, v_drag, "v_drag");
    
    write_matlab_array(std::cout, u_weights, "u_weights");
    write_matlab_array(std::cout, v_weights, "v_weights");
   
    write_matlab_array(std::cout, u_weight_particle, "u_weight_particle");
    write_matlab_array(std::cout, v_weight_particle, "v_weight_particle");
    
    write_matlab_array(std::cout, u_weight_hair, "u_weight_hair");
    write_matlab_array(std::cout, v_weight_hair, "v_weight_hair");
    
    write_matlab_array(std::cout, u_weight_total, "u_weight_total");
    write_matlab_array(std::cout, v_weight_total, "v_weight_total");
    
    write_matlab_array(std::cout, liquid_phi, "liquid_phi");
    
    
    exit(0);
  }

  //Apply the velocity update
  u_valid.assign(0);
  u_pressure_grad.assign(0.0);
  threadutils::thread_pool::ParallelFor(0, u.nj, [&](int j) {
    for(int i = 1; i < u.ni-1; ++i) {
      int index = i + j*ni;
      if(u_weights(i,j) > 0) {
        if(liquid_phi(i,j) < 0 || liquid_phi(i-1,j) < 0) {
          float theta = 1;
          if(liquid_phi(i,j) >= 0 || liquid_phi(i-1,j) >= 0)
            theta = mathutils::fraction_inside(liquid_phi(i-1,j), liquid_phi(i,j));
          if(theta < 0.01) theta = 0.01;
          scalar pressure_grad = (pressure[index] - pressure[index-1]) / dx / theta;
          u(i,j) -= dt * pressure_grad / rho;
          u_pressure_grad(i, j) = pressure_grad;
        }
        u_valid(i,j) = 1;
      }
    }
  });
  v_valid.assign(0);
  v_pressure_grad.assign(0.0);
  threadutils::thread_pool::ParallelFor(1, v.nj-1, [&](int j) {
    for(int i = 0; i < v.ni; ++i) {
      int index = i + j*ni;
      if(v_weights(i,j) > 0) {
        if(liquid_phi(i,j) < 0 || liquid_phi(i,j-1) < 0) {
          float theta = 1;
          if(liquid_phi(i,j) >= 0 || liquid_phi(i,j-1) >= 0)
            theta = mathutils::fraction_inside(liquid_phi(i,j-1), liquid_phi(i,j));
          if(theta < 0.01) theta = 0.01;
          scalar pressure_grad = (pressure[index] - pressure[index-ni]) / dx / theta;
          v(i,j) -= dt * pressure_grad / rho;
          v_pressure_grad(i, j) = pressure_grad;
        }
        v_valid(i,j) = 1;
      }
    }
  });
  
  int num;
  
  num = u_valid.a.size();
  threadutils::thread_pool::ParallelFor(0, num, [&](int i) {
    if(u_valid.a[i] == 0)
      u.a[i] = 0;
  });
  
  num = v_valid.a.size();
  threadutils::thread_pool::ParallelFor(0, num, [&](int i){
    if(v_valid.a[i] == 0)
      v.a[i] = 0;
  });
}

scalar FluidSim2D::compute_phi_vel(const Vector2s& pos, Vector2s& vel) const
{
  scalar min_phi = std::numeric_limits<scalar>::max();
  
  vel.setZero();
  for(auto& b : boundaries)
  {
    if(!b->is_root()) continue;
    Vector2s temp_vel;
    scalar phi = b->compute_phi_vel(pos, temp_vel);
    if(phi < min_phi) {
      min_phi = phi;
      vel = temp_vel;
    }
  }
  
  return min_phi;
}

void FluidSim2D::init_hair_particles()
{
  const std::vector< HairFlow<2>* >& flows = m_parent->getFilmFlows();
  const VectorXs& x = m_parent->getX();
  const VectorXs& v = m_parent->getV();
  const std::vector<std::pair<int, int> > edges = m_parent->getEdges();
  
  const scalar default_radius = dx / sqrt(3.0) * 0.25;
  const scalar default_distance = default_radius / sqrt(2.0);
  
  for(const HairFlow<2>* flow : flows)
  {
    scalar accu_length = 0.0;
    scalar accu_next = 0.0;
    
    const std::vector< int >& edge_indices = flow->getEdgeIndices();
    const std::vector< std::pair<int, int> >& local_edges = flow->getLocalEdges();
    const VectorXs& eta = flow->getEta();
    const VectorXs& radii_v = flow->getRadiiV();
    int ne = edge_indices.size();
    
    for(int i = 0; i < ne; ++i)
    {
      int eidx = edge_indices[i];
      const std::pair<int, int>& local_e = local_edges[i];
      std::pair<int, int> e = edges[eidx];
      scalar prev_accu_next = accu_next;
      accu_next += (x.segment<2>( m_parent->getDof(e.first) ) - x.segment<2>( m_parent->getDof(e.second) )).norm();
      if(accu_next - accu_length > default_distance) {
        scalar alpha = (default_distance + accu_length - prev_accu_next) / (accu_next - prev_accu_next);
        const Vector2s& x0 = x.segment<2>( m_parent->getDof(e.first) );
        const Vector2s& x1 = x.segment<2>( m_parent->getDof(e.second) );
        const Vector2s& v0 = v.segment<2>( m_parent->getDof(e.first) );
        const Vector2s& v1 = v.segment<2>( m_parent->getDof(e.second) );
        
        Vector2s pos = mathutils::lerp(x0, x1, alpha);
        Vector2s vel = mathutils::lerp(v0, v1, alpha);
        
        int local_idx0 = local_e.first;
        int local_idx1 = local_e.second;
        
        scalar H0 = eta(local_idx0) + radii_v(local_idx0);
        scalar H1 = eta(local_idx1) + radii_v(local_idx1);
        
        scalar radius = sqrt(mathutils::lerp(H0 * H0, H1 * H1, alpha));
        
        particles.push_back(Particle<2>(pos, vel, radius, PT_HAIR, eidx, alpha));
        
        accu_length = accu_next;
      }
    }
  }
}

void FluidSim2D::init_random_particles(const scalar& rl, const scalar& rr, const scalar& rb, const scalar& rt)
{
  particles.clear();
  
  int num_particle = ni * nj;
  
  particles.reserve(num_particle);
  
  scalar srl = rl * ni * dx;
  scalar srr = rr * ni * dx;
  
  scalar srb = rb * nj * dx;
  scalar srt = rt * nj * dx;
  
  for(int i = 0; i < ni; ++i)
  {
    for(int j = 0; j < nj; ++j) {
      scalar x = ((scalar) i + 0.5 + (((scalar)rand() / (scalar)RAND_MAX) * 0.5 - 0.5) ) * dx;
      scalar y = ((scalar) j + 0.5 + (((scalar)rand() / (scalar)RAND_MAX) * 0.5 - 0.5) ) * dx;
      Vector2s pt = Vector2s(x,y) + origin;
      Vector2s vel;
      if(compute_phi_vel(pt, vel) > 0 && x >= srl && x <= srr && y >= srb && y <= srt )
        add_particle(Particle<2>(pt, Vector2s::Zero(), dx / sqrt(2.0), PT_LIQUID));
    }
  }
}

void FluidSim2D::sort_particles()
{
  m_sorter->sort(particles.size(), sorter_callback);
}


void FluidSim2D::map_p2g(bool with_hair_particles)
{
  const scalar rho = m_parent->getLiquidDensity();
  
  //u-component of velocity
  threadutils::thread_pool::ParallelFor(0, nj, [&] (int j){
    //for(int j = 0; j < nj; ++j)
    for(int i = 0; i < ni+1; ++i) {
      Vector2s pos = Vector2s(i*dx, (j+0.5)*dx) + origin;
      scalar sumw = 0.0;
      scalar sumu = 0.0;
      
      m_sorter->getNeigboringParticles_cell(i, j, 0, -1, 0, -1, 1, 0, 0, [&] (int pidx) {
        Particle<2>& p = particles[pidx];
        Vector2s diff = p.x - pos;
        
        scalar w = dropvol(p.radii) * rho * mathutils::linear_kernel(diff, dx);
        sumu += w * (p.v(0) - p.c.col(0).dot(diff));
        sumw += w;
      });

      u_particle(i, j) = sumw ? sumu / sumw : 0.0;
      u_weight_particle(i, j) = sumw;
    }
  });
                                        
  //v-component of velocity
  threadutils::thread_pool::ParallelFor(0, ni, [&] (int i){
    for(int j = 0; j < nj+1; ++j) {
      Vector2s pos = Vector2s((i+0.5)*dx, j*dx) + origin;

      scalar sumw = 0.0;
      scalar sumu = 0.0;
      m_sorter->getNeigboringParticles_cell(i, j, 0, -1, 1, -1, 0, 0, 0, [&] (int pidx) {
        Particle<2>& p = particles[pidx];
        Vector2s diff = p.x - pos;
        
        scalar w = dropvol(p.radii) * rho * linear_kernel(diff, dx);
        sumu += w * (p.v(1) - p.c.col(1).dot(diff));
        sumw += w;
      });
      
      v_particle(i, j) = sumw ? sumu / sumw : 0.0;
      v_weight_particle(i, j) = sumw;
    }
  });
}


void FluidSim2D::map_g2p_apic()
{
  int np = particles.size();
  threadutils::thread_pool::ParallelFor(0, np, [&] (int k){
    Particle<2>& p = particles[k];
    p.v = get_velocity(p.x);
    p.c = get_affine_matrix(p.x);
    p.pressure = get_pressure(p.x);
  });
}

void FluidSim2D::prepare_update_from_hair()
{
  int nflows = m_parent->getNumFlows();
  if(u_edge_vel_drag.size() != nflows) u_edge_vel_drag.resize(nflows);
  if(v_edge_vel_drag.size() != nflows) v_edge_vel_drag.resize(nflows);
  
  if(u_num_edge_voxel_intersections.size() != nflows) u_num_edge_voxel_intersections.resize(nflows);
  if(v_num_edge_voxel_intersections.size() != nflows) v_num_edge_voxel_intersections.resize(nflows);
}

void FluidSim2D::done_update_from_hair()
{
  int nflows = m_parent->getNumFlows();
  if(!nflows) return;
  
  threadutils::thread_pool::ParallelFor(0, nflows, [&] (int i) {
    u_num_edge_voxel_intersections[i] = u_edge_vel_drag[i].size();
    v_num_edge_voxel_intersections[i] = v_edge_vel_drag[i].size();
  });
  
  std::partial_sum(u_num_edge_voxel_intersections.begin(), u_num_edge_voxel_intersections.end(), u_num_edge_voxel_intersections.begin());
  std::partial_sum(v_num_edge_voxel_intersections.begin(), v_num_edge_voxel_intersections.end(), v_num_edge_voxel_intersections.begin());
  
  int usize = u_num_edge_voxel_intersections[nflows - 1];
  int vsize = v_num_edge_voxel_intersections[nflows - 1];
  
  if(usize == 0 && vsize == 0) return;
  
  u_vel_drag.resize(usize);
  v_vel_drag.resize(vsize);
  
  threadutils::thread_pool::ParallelFor(0, nflows, [&] (int i) {
    int u_base_idx = (i == 0) ? 0 : u_num_edge_voxel_intersections[i - 1];
    memcpy(&u_vel_drag[u_base_idx], &u_edge_vel_drag[i][0], u_edge_vel_drag[i].size() * sizeof(EdgeVelDragIntersection<2>));
    
    int v_base_idx = (i == 0) ? 0 : v_num_edge_voxel_intersections[i - 1];
    memcpy(&v_vel_drag[v_base_idx], &v_edge_vel_drag[i][0], v_edge_vel_drag[i].size() * sizeof(EdgeVelDragIntersection<2>));
  });
  
  m_sorter->sort(u_vel_drag.size(), [&] (int pidx, int& i, int& j, int& k) {
    i = std::max(0, std::min(m_sorter->ni - 1, u_vel_drag[pidx].coord(0)));
    j = std::max(0, std::min(m_sorter->nj - 1, u_vel_drag[pidx].coord(1)));
    k = 0;
  });
  
  const scalar rho_L = m_parent->getLiquidDensity();
  
  threadutils::thread_pool::ParallelFor(0, u.nj, [&] (int j){
    for(int i = 0; i < u.ni; ++i) {
      scalar sum_vel = 0.0;
      scalar sum_drag = 0.0;
      scalar sum_weight = 0.0;
      scalar sum_linear_weight = 0.0;
      scalar sum_vol = 0.0;
      
      m_sorter->getCellAt(i, j, 0, [&] (int idx) {
        const EdgeVelDragIntersection<2>& inter = u_vel_drag[idx];
        sum_vel += inter.vel_weighted;
        sum_drag += inter.drag_weighted;
        sum_weight += inter.weight;
        sum_linear_weight += inter.linear_weight;
        sum_vol += inter.vol_weighted;
      });
      
      if(sum_weight > 0) {
        u_hair(i, j) = sum_vel / sum_weight;
        scalar multiplier2 = mathutils::clamp(m_parent->getDragRadiusMultiplier() * m_parent->getDragRadiusMultiplier(), 0.0, sum_linear_weight * dx * dx * dz / sum_vol);
        u_drag(i, j) = sum_drag / (sum_linear_weight * dx * dx * dx * rho_L) * multiplier2;
      } else {
        u_hair(i, j) = 0.0;
        u_drag(i, j) = 0.0;
      }
      
      u_weight_hair(i, j) = sum_weight;
    }
  });
  
  m_sorter->sort(v_vel_drag.size(), [&] (int pidx, int& i, int& j, int& k) {
    i = std::max(0, std::min(m_sorter->ni - 1, v_vel_drag[pidx].coord(0)));
    j = std::max(0, std::min(m_sorter->nj - 1, v_vel_drag[pidx].coord(1)));
    k = 0;
  });
  
  threadutils::thread_pool::ParallelFor(0, v.nj, [&] (int j){
    for(int i = 0; i < v.ni; ++i) {
      scalar sum_vel = 0.0;
      scalar sum_drag = 0.0;
      scalar sum_weight = 0.0;
      scalar sum_linear_weight = 0.0;
      scalar sum_vol = 0.0;
      
      m_sorter->getCellAt(i, j, 0, [&] (int idx) {
        const EdgeVelDragIntersection<2>& inter = v_vel_drag[idx];
        sum_vel += inter.vel_weighted;
        sum_drag += inter.drag_weighted;
        sum_weight += inter.weight;
        sum_linear_weight += inter.linear_weight;
        sum_vol += inter.vol_weighted;
      });
      
      if(sum_weight > 0) {
        v_hair(i, j) = sum_vel / sum_weight;
        scalar multiplier2 = mathutils::clamp(m_parent->getDragRadiusMultiplier() * m_parent->getDragRadiusMultiplier(), 0.0, sum_linear_weight * dx * dx * dz / sum_vol);
        v_drag(i, j) = sum_drag / (sum_linear_weight * dx * dx * dx * rho_L) * multiplier2;
      } else {
        v_hair(i, j) = 0.0;
        v_drag(i, j) = 0.0;
      }
      
      v_weight_hair(i, j) = sum_weight;
    }
  });
}

const std::vector<Particle<2> >& FluidSim2D::get_particles() const
{
  return particles;
}

scalar FluidSim2D::cfl()
{
  scalar maxvel = 0;
  const int npu = u.a.size();
  const int npv = v.a.size();
  for(int i = 0; i < npu; ++i)
    maxvel = max(maxvel, fabs(u.a[i]));
  for(int i = 0; i < npv; ++i)
    maxvel = max(maxvel, fabs(v.a[i]));
  return dx / maxvel;
}

const std::vector< FluidSim2D::Boundary<2>* >& FluidSim2D::get_boundaries() const
{
  return boundaries;
}

const std::vector< FluidSim2D::SourceBoundary<2>* >& FluidSim2D::get_sources() const
{
  return sources;
}

const Vector2s& FluidSim2D::get_origin() const
{
  return origin;
}

int FluidSim2D::get_ni() const
{
  return ni;
}

int FluidSim2D::get_nj() const
{
  return nj;
}

int FluidSim2D::get_u_ni() const
{
  return u.ni;
}

int FluidSim2D::get_v_ni() const
{
  return v.ni;
}

int FluidSim2D::get_u_nj() const
{
  return u.nj;
}

int FluidSim2D::get_v_nj() const
{
  return v.nj;
}

std::vector< std::vector<EdgeVelDragIntersection<2> > >& FluidSim2D::get_u_edge_vel_drag()
{
  return u_edge_vel_drag;
}

std::vector< std::vector<EdgeVelDragIntersection<2> > >& FluidSim2D::get_v_edge_vel_drag()
{
  return v_edge_vel_drag;
}

scalar FluidSim2D::dropvol(const scalar& radii) const
{
  return M_PI * radii * radii * dz;
}

scalar FluidSim2D::dropradius(const scalar& vol) const
{
  return sqrt(vol / M_PI / dz);
}

scalar FluidSim2D::cellsize() const
{
  return dx;
}

void FluidSim2D::save_pressure(const std::string szfn)
{
  using namespace std;
  
  ofstream ofs(szfn.c_str());
  ofs << "pressure = [" << endl;
  for(int j = 0; j < nj; ++j) {
    for(int i = 0; i < ni; ++i) {
      ofs << pressure[j * ni + i] << "\t";
    }
    ofs << "\n";
  }
  
  ofs << "];" << endl;
  ofs.close();
}

void FluidSim2D::save_particles_off(const std::string szfn)
{
  using namespace std;
  
  ofstream ofs(szfn.c_str());
  ofs << "OFF" << endl;
  ofs << particles.size() << " 0 0" << endl;
  for(auto& p : particles)
  {
    ofs << p.x(0) << " " << p.x(1) << " " << 0.0 << " " << p.v(0) << " " << p.v(1) << " " << 0.0 << " " << p.radii << endl;
  }
  ofs.close();
}

void FluidSim2D::load_particles_off(const std::string szfn)
{
  particles.clear();
  
  using namespace std;
  
  ifstream ifs(szfn.c_str());
  
  string line;
  
  getline(ifs, line);
  if(line.find("OFF") == string::npos) {
    cerr << "Not an OFF file!" << endl;
    exit(0);
  }
  
  getline(ifs, line);
  vector<string> buf = stringutils::split(line, ' ');
  int nexpected = 0;
  if(buf.size() == 0) {
    cerr << "Not an OFF file!" << endl;
    exit(0);
  }
  
  stringutils::extractFromString(buf[0], nexpected);
  particles.reserve(nexpected);
  
  while(getline(ifs, line))
  {
    buf = stringutils::split(line, ' ');
    
    Particle<2> p;
    
    if(buf.size() >= 1) stringutils::extractFromString(buf[0], p.x(0));
    
    if(buf.size() >= 2) stringutils::extractFromString(buf[1], p.x(1));
    
    if(buf.size() >= 4) stringutils::extractFromString(buf[3], p.v(0));
    
    if(buf.size() >= 5) stringutils::extractFromString(buf[4], p.v(1));
    
    if(buf.size() >= 7) stringutils::extractFromString(buf[6], p.radii);
    
    particles.push_back(p);
  }
  
  ifs.close();
  
  m_sorter->sort(particles.size(), sorter_callback);
  
  compute_liquid_phi();
}

Vector2s FluidSim2D::computeParticleMomentum()
{
  int np = particles.size();
  
  Vector2s sum = Vector2s::Zero();
  
  scalar rho = m_parent->getLiquidDensity();
  
  for(int i = 0; i < np; ++i)
  {
    auto& p = particles[i];
    scalar m = M_PI * p.radii * p.radii * dz * rho;
    sum += p.v * m;
  }
  
  return sum;
}

scalar FluidSim2D::computeParticleAngularMomentum()
{
  int np = particles.size();
  
  scalar sum = 0.0;
  
  scalar rho = m_parent->getLiquidDensity();
  
  for(int i = 0; i < np; ++i)
  {
    auto& p = particles[i];
    scalar m = M_PI * p.radii * p.radii * dz * rho;
    
    Vector2s ppu = (p.x - origin) / dx - Vector2s(0.0, 0.5);
    Vector2s ppv = (p.x - origin) / dx - Vector2s(0.5, 0.0);
    
    int iu_low = std::max(0, std::min(u.ni-2, (int) ppu(0)));
    int iu_high = iu_low + 1;
    
    int ju_low = std::max(0, std::min(u.nj-2, (int) ppu(1)));
    int ju_high = ju_low + 1;
    
    scalar iu_frac = ppu(0) - floor(ppu(0));
    scalar ju_frac = ppu(1) - floor(ppu(1));
    
    Vector2s xu00 = Vector2s(iu_low*dx, (ju_low+0.5)*dx) + origin;
    Vector2s xu10 = Vector2s(iu_high*dx, (ju_low+0.5)*dx) + origin;
    Vector2s xu01 = Vector2s(iu_low*dx, (ju_high+0.5)*dx) + origin;
    Vector2s xu11 = Vector2s(iu_high*dx, (ju_high+0.5)*dx) + origin;
    
    Vector2s pu00 = m * (1.0 - iu_frac) * (1.0 - ju_frac) * Vector2s(p.v(0) + p.c.col(0).dot(xu00 - p.x), 0.0);
    Vector2s pu10 = m * iu_frac * (1.0 - ju_frac) * Vector2s(p.v(0) + p.c.col(0).dot(xu10 - p.x), 0.0);
    Vector2s pu01 = m * (1.0 - iu_frac) * ju_frac * Vector2s(p.v(0) + p.c.col(0).dot(xu01 - p.x), 0.0);
    Vector2s pu11 = m * iu_frac * ju_frac * Vector2s(p.v(0) + p.c.col(0).dot(xu11 - p.x), 0.0);
    
    sum += mathutils::cross2(xu00, pu00) + mathutils::cross2(xu01, pu01) + mathutils::cross2(xu10, pu10) + mathutils::cross2(xu11, pu11);
    
    int iv_low = std::max(0, std::min(v.ni-2, (int) ppv(0)));
    int iv_high = iv_low + 1;
    
    int jv_low = std::max(0, std::min(v.nj-2, (int) ppv(1)));
    int jv_high = jv_low + 1;
    
    scalar iv_frac = ppv(0) - floor(ppv(0));
    scalar jv_frac = ppv(1) - floor(ppv(1));
    
    Vector2s xv00 = Vector2s((iv_low+0.5)*dx, jv_low*dx) + origin;
    Vector2s xv10 = Vector2s((iv_high+0.5)*dx, jv_low*dx) + origin;
    Vector2s xv01 = Vector2s((iv_low+0.5)*dx, jv_high*dx) + origin;
    Vector2s xv11 = Vector2s((iv_high+0.5)*dx, jv_high*dx) + origin;
    
    Vector2s pv00 = m * (1.0 - iv_frac) * (1.0 - jv_frac) * Vector2s(0.0, p.v(1) + p.c.col(1).dot(xv00 - p.x));
    Vector2s pv10 = m * iv_frac * (1.0 - jv_frac) * Vector2s(0.0, p.v(1) + p.c.col(1).dot(xv10 - p.x));
    Vector2s pv01 = m * (1.0 - iv_frac) * jv_frac * Vector2s(0.0, p.v(1) + p.c.col(1).dot(xv01 - p.x));
    Vector2s pv11 = m * iv_frac * jv_frac * Vector2s(0.0, p.v(1) + p.c.col(1).dot(xv11 - p.x));
    
    sum += mathutils::cross2(xv00, pv00) + mathutils::cross2(xv01, pv01) + mathutils::cross2(xv10, pv10) + mathutils::cross2(xv11, pv11);
  }
  
  return sum;
}

int FluidSim2D::num_particles() const
{
  return particles.size();
}

scalar FluidSim2D::computeParticleGridAngularMomentum()
{
  scalar sum = 0.0;
  
  for(int j = 0; j < u_particle.nj; ++j) for(int i = 0; i < u_particle.ni; ++i)
  {
    Vector2s pos = Vector2s(i*dx, (j+0.5)*dx) + origin;
    Vector2s u_vec = Vector2s(u_weight_particle(i,j) * u_particle(i,j), 0.0);
    sum += mathutils::cross2(pos, u_vec);
  }
  
  for(int j = 0; j < v_particle.nj; ++j) for(int i = 0; i < v_particle.ni; ++i)
  {
    Vector2s pos = Vector2s((i+0.5)*dx, j*dx) + origin;
    Vector2s v_vec = Vector2s(0.0, v_weight_particle(i,j) * v_particle(i,j));
    sum += mathutils::cross2(pos, v_vec);
  }
  
  return sum;
}

scalar FluidSim2D::computeHairGridAngularMomentum()
{
  scalar sum = 0.0;
  
  for(int j = 0; j < u_hair.nj; ++j) for(int i = 0; i < u_hair.ni; ++i)
  {
    Vector2s pos = Vector2s(i*dx, (j+0.5)*dx) + origin;
    Vector2s u_vec = Vector2s(u_weight_hair(i,j) * u_hair(i,j), 0.0);
    sum += mathutils::cross2(pos, u_vec);
  }
  
  for(int j = 0; j < v_hair.nj; ++j) for(int i = 0; i < v_hair.ni; ++i)
  {
    Vector2s pos = Vector2s((i+0.5)*dx, j*dx) + origin;
    Vector2s v_vec = Vector2s(0.0, v_weight_hair(i,j) * v_hair(i,j));
    sum += mathutils::cross2(pos, v_vec);
  }
  
  return sum;
}

scalar FluidSim2D::computeCombinedGridAngularMomentum()
{
  scalar sum = 0.0;
  
  for(int j = 0; j < u.nj; ++j) for(int i = 0; i < u.ni; ++i)
  {
    Vector2s pos = Vector2s(i*dx, (j+0.5)*dx) + origin;
    Vector2s u_vec = Vector2s(u_weight_total(i,j) * u(i,j), 0.0);
    sum += mathutils::cross2(pos, u_vec);
  }
  
  for(int j = 0; j < v.nj; ++j) for(int i = 0; i < v.ni; ++i)
  {
    Vector2s pos = Vector2s((i+0.5)*dx, j*dx) + origin;
    Vector2s v_vec = Vector2s(0.0, v_weight_total(i,j) * v(i,j));
    sum += mathutils::cross2(pos, v_vec);
  }
  
  return sum;
}

scalar FluidSim2D::computeParticleWeightedCombinedGridAngularMomentum()
{
  scalar sum = 0.0;
  
  for(int j = 0; j < u.nj; ++j) for(int i = 0; i < u.ni; ++i)
  {
    Vector2s pos = Vector2s(i*dx, (j+0.5)*dx) + origin;
    Vector2s u_vec = Vector2s(u_weight_particle(i,j) * u(i,j), 0.0);
    sum += mathutils::cross2(pos, u_vec);
  }
  
  for(int j = 0; j < v.nj; ++j) for(int i = 0; i < v.ni; ++i)
  {
    Vector2s pos = Vector2s((i+0.5)*dx, j*dx) + origin;
    Vector2s v_vec = Vector2s(0.0, v_weight_particle(i,j) * v(i,j));
    sum += mathutils::cross2(pos, v_vec);
  }
  
  return sum;
}

scalar FluidSim2D::computeHairWeightedCombinedGridAngularMomentum()
{
  scalar sum = 0.0;
  
  for(int j = 0; j < u.nj; ++j) for(int i = 0; i < u.ni; ++i)
  {
    Vector2s pos = Vector2s(i*dx, (j+0.5)*dx) + origin;
    Vector2s u_vec = Vector2s(u_weight_hair(i,j) * u(i,j), 0.0);
    sum += mathutils::cross2(pos, u_vec);
  }
  
  for(int j = 0; j < v.nj; ++j) for(int i = 0; i < v.ni; ++i)
  {
    Vector2s pos = Vector2s((i+0.5)*dx, j*dx) + origin;
    Vector2s v_vec = Vector2s(0.0, v_weight_hair(i,j) * v(i,j));
    sum += mathutils::cross2(pos, v_vec);
  }
  
  return sum;
}

scalar FluidSim2D::computeReweightedHairGridAngularMomentum()
{
  scalar sum = 0.0;
  
  for(int j = 0; j < u_hair.nj; ++j) for(int i = 1; i < u_hair.ni-1; ++i)
  {
    scalar theta = 1.0;
    if(liquid_phi(i,j) >= 0 || liquid_phi(i-1,j) >= 0)
      theta = mathutils::fraction_inside(liquid_phi(i-1,j), liquid_phi(i,j));
    
    Vector2s pos = Vector2s(i*dx, (j+0.5)*dx) + origin;
    Vector2s u_vec = Vector2s(u_weight_hair(i,j) * u_hair(i,j) * (1.0 - theta), 0.0);
    sum += mathutils::cross2(pos, u_vec);
  }
  
  for(int j = 1; j < v_hair.nj-1; ++j) for(int i = 0; i < v_hair.ni; ++i)
  {
    scalar theta = 1.0;
    if(liquid_phi(i,j) >= 0 || liquid_phi(i,j-1) >= 0)
      theta = mathutils::fraction_inside(liquid_phi(i,j-1), liquid_phi(i,j));
    
    Vector2s pos = Vector2s((i+0.5)*dx, j*dx) + origin;
    Vector2s v_vec = Vector2s(0.0, v_weight_hair(i,j) * v_hair(i,j) * (1.0 - theta));
    sum += mathutils::cross2(pos, v_vec);
  }
  
  return sum;
}

scalar FluidSim2D::computeReweightedParticleGridAngularMomentum()
{
  scalar sum = 0.0;
  
  for(int j = 0; j < u_particle.nj; ++j) for(int i = 1; i < u_particle.ni-1; ++i)
  {
    scalar theta = 1.0;
    if(liquid_phi(i,j) >= 0 || liquid_phi(i-1,j) >= 0)
      theta = mathutils::fraction_inside(liquid_phi(i-1,j), liquid_phi(i,j));
    
    Vector2s pos = Vector2s(i*dx, (j+0.5)*dx) + origin;
    Vector2s u_vec = Vector2s(u_weight_particle(i,j) * u_particle(i,j) * theta, 0.0);
    sum += mathutils::cross2(pos, u_vec);
  }
  
  for(int j = 1; j < v_particle.nj-1; ++j) for(int i = 0; i < v_particle.ni; ++i)
  {
    scalar theta = 1.0;
    if(liquid_phi(i,j) >= 0 || liquid_phi(i,j-1) >= 0)
      theta = mathutils::fraction_inside(liquid_phi(i,j-1), liquid_phi(i,j));
    
    Vector2s pos = Vector2s((i+0.5)*dx, j*dx) + origin;
    Vector2s v_vec = Vector2s(0.0, v_weight_particle(i,j) * v_particle(i,j) * theta);
    sum += mathutils::cross2(pos, v_vec);
  }
  
  return sum;
}

scalar FluidSim2D::computeTotalLiquidVol() const
{
  scalar sum = 0.0;
  for(auto& p : particles)
  {
    sum += M_PI * p.radii * p.radii * dz;
  }
  return sum;
}

Vector2s FluidSim2D::computeHairGridMomentum()
{
  Vector2s sum = Vector2s::Zero();
  
  for(int j = 0; j < u_hair.nj; ++j) for(int i = 0; i < u_hair.ni; ++i)
  {
    sum(0) += u_weight_hair(i, j) * u_hair(i, j);
  }
  
  for(int j = 0; j < v_hair.nj; ++j) for(int i = 0; i < v_hair.ni; ++i)
  {
    sum(1) += v_weight_hair(i, j) * v_hair(i, j);
  }
  
  return sum;
}

Vector2s FluidSim2D::computeParticleGridMomentum()
{
  Vector2s sum = Vector2s::Zero();
  
  for(int j = 0; j < u_hair.nj; ++j) for(int i = 0; i < u_hair.ni; ++i)
  {
    sum(0) += u_weight_particle(i, j) * u_particle(i, j);
  }
  
  for(int j = 0; j < v_hair.nj; ++j) for(int i = 0; i < v_hair.ni; ++i)
  {
    sum(1) += v_weight_particle(i, j) * v_particle(i, j);
  }
  
  return sum;
}


Vector2s FluidSim2D::computeReweightedHairGridMomentum()
{
  Vector2s sum = Vector2s::Zero();
  
  for(int j = 0; j < u_hair.nj; ++j) for(int i = 1; i < u_hair.ni - 1; ++i)
  {
    scalar theta = 1.0;
    if(liquid_phi(i,j) >= 0 || liquid_phi(i-1,j) >= 0)
      theta = mathutils::fraction_inside(liquid_phi(i-1,j), liquid_phi(i,j));
    
    sum(0) += u_weight_hair(i, j) * u_hair(i, j) * (1. - theta);
  }
  
  for(int j = 1; j < v_hair.nj - 1; ++j) for(int i = 0; i < v_hair.ni; ++i)
  {
    scalar theta = 1.0;
    if(liquid_phi(i,j) >= 0 || liquid_phi(i,j-1) >= 0)
      theta = mathutils::fraction_inside(liquid_phi(i,j-1), liquid_phi(i,j));
    
    sum(1) += v_weight_hair(i, j) * v_hair(i, j) * (1. - theta);
  }
  
  return sum;
}

Vector2s FluidSim2D::computeReweightedParticleGridMomentum()
{
  Vector2s sum = Vector2s::Zero();
  
  for(int j = 0; j < u_particle.nj; ++j) for(int i = 1; i < u_particle.ni-1; ++i)
  {
    scalar theta = 1.0;
    if(liquid_phi(i,j) >= 0 || liquid_phi(i-1,j) >= 0)
      theta = mathutils::fraction_inside(liquid_phi(i-1,j), liquid_phi(i,j));
    
    sum(0) += u_weight_particle(i, j) * u_particle(i, j) * theta;
  }
  
  for(int j = 1; j < v_particle.nj-1; ++j) for(int i = 0; i < v_particle.ni; ++i)
  {
    scalar theta = 1.0;
    if(liquid_phi(i,j) >= 0 || liquid_phi(i,j-1) >= 0)
      theta = mathutils::fraction_inside(liquid_phi(i,j-1), liquid_phi(i,j));
    
    sum(1) += v_weight_particle(i, j) * v_particle(i, j) * theta;
  }
  
  return sum;
}

Vector2s FluidSim2D::computeParticleWeightedCombinedGridMomentum()
{
  Vector2s sum = Vector2s::Zero();
  
  for(int j = 0; j < u.nj; ++j) for(int i = 1; i < u.ni-1; ++i)
  {
    sum(0) += u_weight_particle(i, j) * u(i, j);
  }
  
  for(int j = 1; j < v.nj-1; ++j) for(int i = 0; i < v.ni; ++i)
  {
    sum(1) += v_weight_particle(i, j) * v(i, j);
  }
  
  return sum;
}

Vector2s FluidSim2D::computeHairWeightedCombinedGridMomentum()
{
  Vector2s sum = Vector2s::Zero();
  
  for(int j = 0; j < u.nj; ++j) for(int i = 1; i < u.ni-1; ++i)
  {
    sum(0) += u_weight_hair(i, j) * u(i, j);
  }
  
  for(int j = 1; j < v.nj-1; ++j) for(int i = 0; i < v.ni; ++i)
  {
    sum(1) += v_weight_hair(i, j) * v(i, j);
  }
  
  return sum;
}

Vector2s FluidSim2D::computeCombinedGridMomentum()
{
  Vector2s sum = Vector2s::Zero();
  
  for(int j = 0; j < u.nj; ++j) for(int i = 0; i < u.ni; ++i)
  {
    sum(0) += u_weight_total(i, j) * u(i, j);
  }
  
  for(int j = 0; j < v.nj; ++j) for(int i = 0; i < v.ni; ++i)
  {
    sum(1) += v_weight_total(i, j) * v(i, j);
  }
  
  return sum;
}

scalar FluidSim2D::computeParticleKineticEnergy()
{
  scalar T = 0.0;
  
  scalar rho = m_parent->getLiquidDensity();
  
  int np = particles.size();
  
  for( int i = 0; i < np; ++i ) {
    auto& p = particles[i];
    scalar m = M_PI * p.radii * p.radii * dz * rho;
    T += m * p.v.squaredNorm();
  }
  return 0.5*T;
}

scalar FluidSim2D::computeHairGridKineticEnergy()
{
  scalar sum = 0;
  
  for(int j = 0; j < nj; ++j) for(int i = 0; i < ni; ++i)
  {
    Vector2s pos((i+0.5) * dx, (j+0.5) * dx);
    Vector2s velocity = get_velocity(pos, u_hair, v_hair);
    scalar mass = get_velocity(pos, u_weight_hair, v_weight_hair).norm();
    
    sum += mass * velocity.squaredNorm();
  }
  
  return sum * 0.5;
}

scalar FluidSim2D::computeParticleGridKineticEnergy()
{
  scalar sum = 0;
  
  for(int j = 0; j < nj; ++j) for(int i = 0; i < ni; ++i)
  {
    Vector2s pos((i+0.5) * dx, (j+0.5) * dx);
    Vector2s velocity = get_velocity(pos, u_particle, v_particle);
    scalar mass = get_velocity(pos, u_weight_particle, v_weight_particle).norm();
    
    sum += mass * velocity.squaredNorm();
  }
  
  return sum * 0.5;
}

scalar FluidSim2D::computeCombinedGridKineticEnergy()
{
  scalar sum = 0;
  
  for(int j = 0; j < nj; ++j) for(int i = 0; i < ni; ++i)
  {
    Vector2s pos((i+0.5) * dx, (j+0.5) * dx);
    Vector2s velocity = get_velocity(pos);
    scalar mass = get_velocity(pos, u_weight_total, v_weight_total).norm();
    
    sum += mass * velocity.squaredNorm();
  }
  
  return sum * 0.5;
}

Vector2s FluidSim2D::computeHairGridDrag()
{
  Vector2s F = Vector2s::Zero();
  
  for(int j = 0; j < u.nj; ++j) for(int i = 0; i < u.ni; ++i)
  {
    F(0) += u_drag(i, j) * u_weight_hair(i, j);
  }
  
  for(int j = 0; j < v.nj; ++j) for(int i = 0; i < v.ni; ++i)
  {
    F(1) += v_drag(i, j) * v_weight_hair(i, j);
  }
  
  return F;
}

Vector3s FluidSim2D::getMinBBX() const
{
  return Vector3s(origin(0), origin(1), 0.0);
}

Vector3s FluidSim2D::getMaxBBX() const
{
  return Vector3s(origin(0) + ni * dx, origin(1) + nj * dx, 0.0);
}

void FluidSim2D::controlSources(const scalar& current_time, const scalar& dt)
{
  for(auto& s : sources)
  {
    s->activated = (current_time >= s->start && current_time + dt < s->end);
  }
}

void FluidSim2D::write(std::vector<scalar>& data) const
{
  for(auto& p : particles)
  {
    p.write(data);
  }
}

void FluidSim2D::read(const scalar* data, size_t size_particles, size_t size_boundaries)
{
  
  size_t np = size_particles / Particle<2>::size();
  particles.resize(np);
  
  size_t k = 0;
  for(size_t i = 0; i < np; ++i)
  {
    particles[i].read(data + k);
    k += Particle<2>::size() / sizeof(scalar);
  }
  
  int i = 0;
  int num = (size_particles + size_boundaries) / sizeof(scalar);
  while(k < num)
  {
    if(boundaries[i]->type == BT_UNION || boundaries[i]->type == BT_INTERSECT) continue;
    boundaries[i]->read(data + k);
    k += boundaries[i]->size() / sizeof(scalar);
    ++i;
  }
}

void FluidSim2D::writeReadable(std::vector<std::ostringstream>& oss) const
{
  const scalar default_radius = dx * default_radius_multiplier();
  const std::vector< std::pair<int, int> >& edges = m_parent->getEdges();
  const VectorXs& radius = m_parent->getRadii();
  
  for(auto& p : particles)
  {
    if(p.type == PT_HAIR) {
      const std::pair<int, int>& e = edges[ p.edge_idx ];
      const scalar radii_c = mathutils::lerp(radius(e.first), radius(e.second), p.edge_alpha);
      if(p.radii < 1.01 * radii_c) continue;
    }
    
    int level = std::max(0, std::min((int) oss.size() - 1, (int) (log( default_radius / p.radii ) / log(2.0))));
    
    oss[level] << p.x.transpose() << " " << p.radii << "\n";
  }
}

void FluidSim2D::readReadable( std::ifstream& file )
{
  std::cerr << "readReadable NOT IMPLEMENTED!" << std::endl;
}

size_t FluidSim2D::particle_size() const
{
  return particles.size() * Particle<2>::size();
}

size_t FluidSim2D::boundary_size() const
{
  int sum = 0;
  for(auto& b : boundaries)
  {
    if(b->type == BT_UNION || b->type == BT_INTERSECT) continue;
    sum += b->size();
  }
  
  return sum;
}

size_t FluidSim2D::crucial_grid_size() const
{
  return (u.a.size() + v.a.size() + u_pressure_grad.a.size() + v_pressure_grad.a.size() + liquid_phi.a.size()) * sizeof(scalar);
}

void FluidSim2D::add_particle(const VectorXs& pos, const VectorXs& vel, const scalar& radii, ParticleType type)
{
  add_particle(Particle<2>(pos.segment<2>(0), vel.segment<2>(0), radii, type));
}

void FluidSim2D::preCompute( const VectorXs& x, const VectorXs& v, const VectorXs& m, const scalar& dt )
{
  int nforces = drag_forces.size();
  threadutils::thread_pool::ParallelFor(0, nforces, [&] (int i) {
    drag_forces[i]->preCompute(x, v, m, dt);
  });
}

void FluidSim2D::addGradEToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, VectorXs& gradE )
{
  int nforces = drag_forces.size();
  threadutils::thread_pool::ParallelFor(0, nforces, [&] (int i) {
    drag_forces[i]->addGradEToTotal(x, v, m, gradE);
  });
}

scalar FluidSim2D::get_particle_weight(const Vector2s& position) const
{
  return get_velocity(position, u_weight_particle, v_weight_particle).norm();
}

scalar FluidSim2D::default_radius_multiplier() const
{
  return 1.0 / sqrt(2.0) / 2.0;
}

int FluidSim2D::default_particle_in_cell() const
{
  return 4;
}

void FluidSim2D::apply_viscosity(scalar dt)
{
  std::cerr << "FluidSim2D::apply_viscosity UNIMPLIMENTED!" << std::endl;
}

scalar FluidSim2D::get_clamped_particle_weight(const Vector2s& position) const
{
  const scalar standard_radius = dx * default_radius_multiplier();
  const scalar standard_vol = dropvol(standard_radius);
  scalar w = get_particle_weight(position);
  scalar criterion = m_parent->getBulkThresholdMultiplier() * standard_vol * m_parent->getLiquidDensity();
  return mathutils::clamp(w / criterion, 0.0, 1.0);
}

//Apply several iterations of a very simple "Jacobi"-style propagation of valid velocity data in all directions
void extrapolate(Array2s& grid, Array2s& old_grid, const Array2s& grid_weight, const Array2s& grid_liquid_weight, Array2c& valid, Array2c old_valid, const Vector2i& offset) {
  
  //Initialize the list of valid cells
  for(int j = 0; j < valid.nj; ++j) valid(0,j) = valid(valid.ni-1,j) = 0;
  for(int i = 0; i < valid.ni; ++i) valid(i,0) = valid(i,valid.nj-1) = 0;
  
  for(int j = 1; j < grid.nj - 1; ++j) for(int i = 1; i < grid.ni - 1; ++i)
    valid(i,j) = grid_weight(i,j) > 0;// && (grid_liquid_weight(i, j) < 0 || grid_liquid_weight(i + offset(0), j + offset(1)) < 0 );
  
  Array2s* pgrid[2] = {&grid, &old_grid};
  Array2c* pvalid[2] = {&valid, &old_valid};
  
  for(int layers = 0; layers < 4; ++layers) {
    Array2s* pgrid_source = pgrid[layers & 1];
    Array2s* pgrid_target = pgrid[!(layers & 1)];
    
    Array2c* pvalid_source = pvalid[layers & 1];
    Array2c* pvalid_target = pvalid[!(layers & 1)];
    threadutils::thread_pool::ParallelFor(1, grid.nj - 1, [&] (int j) {
      for(int i = 1; i < grid.ni-1; ++i) {
        scalar sum = 0;
        int count = 0;
        
        if(!(*pvalid_source)(i,j)) {
          
          if((*pvalid_source)(i+1,j)) {
            sum += (*pgrid_source)(i+1,j);
            ++count;
          }
          if((*pvalid_source)(i-1,j)) {
            sum += (*pgrid_source)(i-1,j);
            ++count;
          }
          if((*pvalid_source)(i,j+1)) {
            sum += (*pgrid_source)(i,j+1);
            ++count;
          }
          if((*pvalid_source)(i,j-1)) {
            sum += (*pgrid_source)(i,j-1);
            ++count;
          }
          
          //If any of neighbour cells were valid,
          //assign the cell their average value and tag it as valid
          if(count > 0) {
            (*pgrid_target)(i,j) = sum /(scalar)count;
            (*pvalid_target)(i,j) = 1;
          }
        }
      }
    });
    
    *pvalid_source = *pvalid_target;
    *pgrid_source = *pgrid_target;
  }
}



