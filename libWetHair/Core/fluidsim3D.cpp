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

#include "fluidsim3D.h"
#include "MathUtilities.h"

#include "array3_utils.h"

#include "pcgsolver/sparse_matrix.h"
#include "pcgsolver/pcg_solver.h"

#include "ThreadUtils.h"
#include "AlgebraicMultigrid.h"

#include "TwoDScene.h"
#include "sorter.h"
#include "FluidDragForce.h"

#include "volume_fractions.h"
#include "viscosity3d.h"

#include <fstream>
#include <numeric>
#include <iomanip>

//#define SURFACE_TENSION_ON_HAIRS

//#define USE_MERGE_PARTICLES

using namespace mathutils;

const static int MAX_BRIDGE_PER_EDGE = 64;

void extrapolate(Array3s& grid, Array3s& old_grid, const Array3s& grid_weight, const Array3s& grid_liquid_weight, Array3c& valid, Array3c old_valid, const Vector3i& offset);
void estimate_volume_fractions(Array3s& volumes,
                               const Vector3s& start_centre, const scalar dx,
                               const Array3s& phi, const Vector3s& phi_origin, const scalar phi_dx);

FluidSim3D* g_fluid3d_iptr = NULL;

static void sorter_callback(int pidx, int& i, int& j, int& k)
{
  auto& particles = g_fluid3d_iptr->get_particles();
  
  auto& p = particles[pidx];
  
  const Vector3s& origin = g_fluid3d_iptr->get_origin();
  scalar cellsize = g_fluid3d_iptr->cellsize();
  
  int pi = (int)((p.x(0) - origin(0)) / cellsize);
  int pj = (int)((p.x(1) - origin(1)) / cellsize);
  int pk = (int)((p.x(2) - origin(2)) / cellsize);
  
  i = max(0, min(g_fluid3d_iptr->get_ni()-1, pi));
  j = max(0, min(g_fluid3d_iptr->get_nj()-1, pj));
  k = max(0, min(g_fluid3d_iptr->get_nk()-1, pk));
}

scalar FluidSim3D::default_radius_multiplier() const
{
  return 1.0 / sqrt(2.0) / 2.0;
}

int FluidSim3D::default_particle_in_cell() const
{
  return 8;
}

FluidSim3D::~FluidSim3D()
{
  if(m_sorter) delete m_sorter;
}

FluidSim3D::FluidSim3D(const Vector3s& origin_, scalar width, int ni_, int nj_, int nk_,
                       const std::vector< Boundary<3>* >& boundaries_, const std::vector< SourceBoundary<3>* >& sources_, TwoDScene<3>* scene)
: m_parent(scene)
{
  g_fluid3d_iptr = this;
  ryoichi_correction_counter = 0;
  origin = origin_;
  boundaries = boundaries_;
  sources = sources_;
  ni = ni_;
  nj = nj_;
  nk = nk_;
  dx = width / (scalar) ni;
  
  u.resize(ni+1,nj,nk); temp_u.resize(ni+1,nj,nk); u_weights.resize(ni+1,nj,nk);
  u_weight_hair.resize(ni+1,nj,nk); u_valid.resize(ni+1,nj,nk); u_hair.resize(ni+1,nj,nk);
  u_particle.resize(ni+1, nj,nk); u_weight_particle.resize(ni+1,nj,nk);
  u_drag.resize(ni+1, nj,nk); u_weight_total.resize(ni+1, nj, nk);
  u_pressure_grad.resize(ni+1,nj,nk); u_solid.resize(ni+1, nj, nk); u_visc_impulse.resize(ni+1,nj,nk);
  
  v.resize(ni,nj+1,nk); temp_v.resize(ni,nj+1,nk); v_weights.resize(ni,nj+1,nk);
  v_weight_hair.resize(ni,nj+1,nk); v_valid.resize(ni,nj+1,nk); v_hair.resize(ni,nj+1,nk);
  v_particle.resize(ni, nj+1, nk); v_weight_particle.resize(ni,nj+1,nk);
  v_drag.resize(ni, nj+1, nk); v_weight_total.resize(ni, nj+1, nk);
  v_pressure_grad.resize(ni, nj+1, nk); v_solid.resize(ni, nj+1, nk); v_visc_impulse.resize(ni,nj+1,nk);
  
  w.resize(ni,nj,nk+1); temp_w.resize(ni,nj,nk+1); w_weights.resize(ni,nj,nk+1);
  w_weight_hair.resize(ni,nj,nk+1); w_valid.resize(ni,nj,nk+1); w_hair.resize(ni,nj,nk+1);
  w_particle.resize(ni, nj, nk+1); w_weight_particle.resize(ni,nj,nk+1);
  w_drag.resize(ni, nj, nk+1); w_weight_total.resize(ni, nj, nk+1);
  w_pressure_grad.resize(ni, nj, nk+1); w_solid.resize(ni, nj, nk+1); w_visc_impulse.resize(ni,nj,nk+1);
  
  u.set_zero();
  v.set_zero();
  w.set_zero();
  u_visc_impulse.set_zero();
  v_visc_impulse.set_zero();
  w_visc_impulse.set_zero();
  u_pressure_grad.set_zero();
  v_pressure_grad.set_zero();
  w_pressure_grad.set_zero();
  u_hair.set_zero();
  v_hair.set_zero();
  w_hair.set_zero();
  u_particle.set_zero();
  v_particle.set_zero();
  w_particle.set_zero();
  
  u_solid.set_zero();
  v_solid.set_zero();
  w_solid.set_zero();
  
  u_weight_hair.set_zero();
  v_weight_hair.set_zero();
  w_weight_hair.set_zero();
  u_weight_particle.set_zero();
  v_weight_particle.set_zero();
  w_weight_particle.set_zero();
  u_drag.set_zero();
  v_drag.set_zero();
  w_drag.set_zero();
  
  u_weight_total.set_zero();
  v_weight_total.set_zero();
  w_weight_total.set_zero();
  
  temp_u.set_zero();
  temp_v.set_zero();
  temp_w.set_zero();
  
  nodal_solid_phi.resize(ni+1,nj+1,nk+1);
  cell_solid_phi.resize(ni,nj,nk);
  valid.resize(ni+1, nj+1, nk+1);
  old_valid.resize(ni+1, nj+1, nk+1);
  liquid_phi.resize(ni, nj, nk);
  liquid_phi.assign(3*dx);
  
  c_vol_liquid.resize(ni, nj, nk, 0);
  u_vol_liquid.resize(ni+1, nj, nk, 0);
  v_vol_liquid.resize(ni, nj+1, nk, 0);
  w_vol_liquid.resize(ni, nj, nk+1, 0);
  ex_vol_liquid.resize(ni, nj+1, nk+1, 0);
  ey_vol_liquid.resize(ni+1, nj, nk+1, 0);
  ez_vol_liquid.resize(ni+1, nj+1, nk, 0);
  
  pressure.resize(ni * nj * nk, 0.0);
 
  m_bridges.reserve(scene->getNumEdges() * MAX_BRIDGE_PER_EDGE);
  
  update_boundary();
  m_sorter = new Sorter(ni, nj, nk);
  particles.clear();
  
  int nflows = m_parent->getNumFlows();
  drag_forces.resize(nflows);
  for(int i = 0; i < nflows; ++i) {
    drag_forces[i] = new FluidDragForce<3>(*scene, i);
    m_parent->insertForce(drag_forces[i]);
  }
  
  m_pool_liquid_index_cache.resize(nk);
  m_pool_liquid_particle_cache.resize(nk);
  m_regular_liquid_particle_cache.resize(nk);
  
  for(auto& s : sources)
  {
    s->sample(this);
  }
}

void FluidSim3D::advect_boundary(const scalar& dt)
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
void FluidSim3D::update_boundary() {
  threadutils::thread_pool::ParallelFor(0, nk+1, [&] (int k) {
    for(int j = 0; j < nj+1; ++j) for(int i = 0; i < ni+1; ++i) {
      Vector3s pos(i*dx,j*dx,k*dx);
      Vector3s vel;
      nodal_solid_phi(i,j,k) = compute_phi_vel(pos + origin, vel);
    }
  });
  
  threadutils::thread_pool::ParallelFor(0, nk, [&] (int k) {
    for(int j = 0; j < nj; ++j) for(int i = 0; i < ni; ++i) {
      Vector3s pos_cell((i + 0.5)*dx,(j + 0.5)*dx,(k + 0.5)*dx);
      Vector3s vel;
      
      cell_solid_phi(i,j,k) = compute_phi_vel(pos_cell + origin, vel);
    }
  });
  
  threadutils::thread_pool::ParallelFor(0, nk, [&] (int k) {
    for(int j = 0; j < nj; ++j) for(int i = 0; i < ni+1; ++i) {
      Vector3s pos(i*dx,(j+0.5)*dx,(k+0.5)*dx);
      Vector3s vel;
      compute_phi_vel(pos + origin, vel);
      u_solid(i, j, k) = vel[0];
    }
  });
  
  threadutils::thread_pool::ParallelFor(0, nk, [&] (int k) {
    for(int j = 0; j < nj+1; ++j) for(int i = 0; i < ni; ++i) {
      Vector3s pos((i+0.5)*dx,j*dx,(k+0.5)*dx);
      Vector3s vel;
      compute_phi_vel(pos + origin, vel);
      v_solid(i, j, k) = vel[1];
    }
  });
  
  threadutils::thread_pool::ParallelFor(0, nk+1, [&] (int k) {
    for(int j = 0; j < nj; ++j) for(int i = 0; i < ni; ++i) {
      Vector3s pos((i+0.5)*dx,(j+0.5)*dx,k*dx);
      Vector3s vel;
      compute_phi_vel(pos + origin, vel);
      w_solid(i, j, k) = vel[2];
    }
  });
}

int FluidSim3D::num_particles() const
{
  return particles.size();
}

scalar FluidSim3D::computeOverallDivergence()
{
  scalar div_sum = 0;
  for(int k = 0; k < nk; ++k) for(int j = 0; j < nj; ++j) for(int i = 0; i < ni; ++i) {
    if(liquid_phi(i,j,k) > 0) continue;
    scalar divi = (u_weights(i + 1, j, k) * u(i + 1, j, k) - u_weights(i, j, k) * u(i, j, k) + v_weights(i, j + 1, k) * v(i, j + 1, k) - v_weights(i, j, k) * v(i, j, k) + w_weights(i, j, k + 1) * u(i, j, k + 1) - w_weights(i, j, k) * w(i, j, k)) / dx;
    div_sum += divi * divi;
  }
  
  return sqrt(div_sum);
}

void FluidSim3D::transferLiquidToGridParticle(const scalar& dt)
{
  // set all particles as old
  int np = particles.size();
  threadutils::thread_pool::ParallelFor(0, np, [&] (int i) {
    particles[i].fresh *= 0.99;
  });
  
  int nhairp = m_parent->getNumParticles();
  const VectorXs& x = m_parent->getX();
  VectorXs& v = m_parent->getV();
  VectorXs& m = m_parent->getM();
  const VectorXs& rest_m = m_parent->getHairRestMass();
  const scalar& rho_liq = m_parent->getLiquidDensity();
  const scalar& gamma_liq = m_parent->getLiquidTension();
  
  m_sorter->sort(nhairp, [&] (int pidx, int& i, int& j, int& k) {
    const Vector3s& pos = x.segment<3>( m_parent->getDof( pidx ) );
    int pi = (int)((pos(0) - origin(0)) / dx);
    int pj = (int)((pos(1) - origin(1)) / dx);
    int pk = (int)((pos(2) - origin(2)) / dx);
    
    i = max(0, min(ni-1, pi));
    j = max(0, min(nj-1, pj));
    k = max(0, min(nk-1, pk));
  });
  
  scalar release_radius = dx * default_radius_multiplier() * sqrt(3.0) * 0.5;
  
  scalar release_vol = dropvol(release_radius);
  
  const std::vector<int>& particle_to_hair = m_parent->getParticleToHairs();
  const std::vector<int>& global_to_local = m_parent->getParticleToHairLocalIndices();
  std::vector<HairFlow<3>*>& hairs = m_parent->getFilmFlows();

  const scalar epsilon = 1e-7 * dx / dt;
  
  const int nflows = hairs.size();
  // record liquid pool into cache
  if(m_pool_liquid_vol_cache.size() != nflows * 2) m_pool_liquid_vol_cache.resize(nflows * 2);
  m_pool_liquid_vol_cache.setZero();
  
  bool dripping_zero_end = m_parent->getParameter().drippingnear;
  bool dripping_far_end = m_parent->getParameter().drippingfar;
  bool dripping_middle = m_parent->getParameter().drippingmiddle;
  
  threadutils::thread_pool::ParallelFor(0, nflows, [&] (int i) {
    const HairFlow<3>* hair = hairs[i];
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
  
  threadutils::thread_pool::ParallelFor(0, nk, [&] (int k) {
    m_pool_liquid_index_cache[k].resize(0);
    m_pool_liquid_particle_cache[k].resize(0);
    m_regular_liquid_particle_cache[k].resize(0);
    
    for(int j = 0; j < nj; ++j) for(int i = 0; i < ni; ++i)
    {
      scalar total_pool_liquid = 0.0;
      scalar total_regular_liquid = 0.0;
      Vector3s avg_vel = Vector3s::Zero();
      Vector3s avg_pos = Vector3s::Zero();
      Vector3s avg_dir = Vector3s::Zero();
      Vector3s avg_pdirpl = Vector3s::Zero();
      
      scalar equivalent_radius = 0.0;
      scalar avg_area_v = 0.0;
      Vector3s avg_accel_norm = Vector3s::Zero();
      
      std::vector<int> hair_indices;
      
      int count = 0;
      m_sorter->getCellAt(i, j, k, [&] (int pidx) {
        int hidx = particle_to_hair[pidx];
        int local_idx = global_to_local[pidx];
        HairFlow<3>* hair = hairs[hidx];
        int nhair_particles = (int) (hair->getParticleIndices().size());
        const MatrixXs& actual_u_v = hair->getActualVertexVelocity();
        const MatrixXs& accel_v = hair->getAccelV();
        const MatrixXs& tangent_v = hair->getTangentV();
        const MatrixXs& tangent_e = hair->getTangentE();
        const VectorXs& eta_v = hair->getEta();
        const VectorXs& radii_v = hair->getRadiiV();
        const VectorXs& area_v = hair->getAreaV();
        const VectorXs& area_e = hair->getAreaE();
        
        // dripping from tips
        if(local_idx == 0) {
          total_pool_liquid += m_pool_liquid_vol_cache(hidx * 2 + 0);
          if(m_pool_liquid_vol_cache(hidx * 2 + 0) > 0.0) {
            avg_pos += x.segment<3>( m_parent->getDof( pidx ) ) * m_pool_liquid_vol_cache(hidx * 2 + 0);
            avg_vel += (v.segment<3>( m_parent->getDof( pidx ) ) + actual_u_v.row(local_idx).transpose()) * m_pool_liquid_vol_cache(hidx * 2 + 0);
            avg_dir += tangent_v.row( local_idx ).transpose() * m_pool_liquid_vol_cache(hidx * 2 + 0);
            avg_pdirpl += ( tangent_v.row( local_idx + 1 ).transpose() - tangent_v.row( local_idx ).transpose() ) / area_e(0) * m_pool_liquid_vol_cache(hidx * 2 + 0);
          }
        } else if(local_idx == nhair_particles - 1) {
          total_pool_liquid += m_pool_liquid_vol_cache(hidx * 2 + 1);
          if(m_pool_liquid_vol_cache(hidx * 2 + 1) > 0.0)
          {
            avg_pos += x.segment<3>( m_parent->getDof( pidx ) ) * m_pool_liquid_vol_cache(hidx * 2 + 1);
            avg_vel += (v.segment<3>( m_parent->getDof( pidx ) ) + actual_u_v.row(local_idx).transpose()) * m_pool_liquid_vol_cache(hidx * 2 + 1);
            avg_dir += tangent_v.row( local_idx ).transpose() * m_pool_liquid_vol_cache(hidx * 2 + 1);
            avg_pdirpl += ( tangent_v.row( local_idx ).transpose() - tangent_v.row( local_idx - 1 ).transpose() ) / area_e(local_idx - 1) * m_pool_liquid_vol_cache(hidx * 2 + 0);
          }
        } else if(dripping_middle) {
          // dripping from regular vertex
          // compute regular volume
          scalar reg_vol = M_PI * ((eta_v(local_idx) + radii_v(local_idx)) * (eta_v(local_idx) + radii_v(local_idx)) - radii_v(local_idx) * radii_v(local_idx)) * area_v(local_idx);
          total_regular_liquid += reg_vol;
          
          // update avg pos and vel
          avg_pos += x.segment<3>( m_parent->getDof( pidx ) ) * reg_vol;
          avg_vel += (v.segment<3>( m_parent->getDof( pidx ) ) + actual_u_v.row(local_idx).transpose()) * reg_vol;
          avg_dir += tangent_v.row( local_idx ).transpose() * reg_vol;
          avg_pdirpl += ( tangent_e.row( local_idx ).transpose() - tangent_e.row( local_idx - 1 ).transpose() ) / area_v(local_idx) * reg_vol;
          
          // update equivalent radius
          equivalent_radius += radii_v(local_idx) * radii_v(local_idx);
          
          // compute normal acceleration
          Vector3s tan_dir = tangent_v.row(local_idx).transpose().normalized();
          Matrix3s proj = Matrix3s::Identity() - tan_dir * tan_dir.transpose();
          
          // acceleration to normal dir
          Vector3s accel_vert = accel_v.row(local_idx).transpose();
          avg_accel_norm += proj * accel_vert * reg_vol;
          
          // average segment length
          avg_area_v += area_v(local_idx) * reg_vol;
          

          
          ++count;
        }
      });
      
      scalar total_liquid = total_regular_liquid + total_pool_liquid;
      
      if(total_liquid == 0.0) continue;
      
      avg_pos /= total_liquid;
      avg_vel /= total_liquid;
      avg_dir /= total_liquid;
      avg_pdirpl /= total_liquid;
      
      if(total_pool_liquid > release_vol) { // release from pool
        scalar actual_release_vol = total_pool_liquid;
        scalar radii_release = dropradius(actual_release_vol);
        
        m_pool_liquid_particle_cache[k].push_back(Particle<3>(avg_pos, avg_vel, radii_release, PT_LIQUID));
        
        // mark released liquid from hairs
        m_sorter->getCellAt(i, j, k, [&] (int pidx) {
          int hidx = particle_to_hair[pidx];
          int local_idx = global_to_local[pidx];
          HairFlow<3>* hair = hairs[hidx];
          int nhair_particles = (int) (hair->getParticleIndices().size());
          
          if(local_idx == 0) {
            m_pool_liquid_index_cache[k].push_back(hidx * 2 + 0);
          }
          
          if(local_idx == nhair_particles - 1) {
            m_pool_liquid_index_cache[k].push_back(hidx * 2 + 1);
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
            Vector3s proj_dir = (avg_dir + avg_pdirpl * rand_dist).normalized();
            
            Vector3s rand_vec_projected = proj_dir * rand_dist;
            m_regular_liquid_particle_cache[k].push_back(Particle<3>(avg_pos + rand_vec_projected, avg_vel, radii_release, PT_LIQUID));
          }
          
          // the proportion to remove liquid
          scalar prop = (total_regular_liquid - actual_release_vol) / total_regular_liquid;
          
          // remove volume from hair
          m_sorter->getCellAt(i, j, k, [&] (int pidx) {
            int hidx = particle_to_hair[pidx];
            int local_idx = global_to_local[pidx];
            HairFlow<3>* hair = hairs[hidx];
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
              m.segment<3>( idof ).setConstant(new_total_mass);
              if(!m_parent->isMassSpring() && local_idx != eta_v.size()) {
                scalar new_inertia = rest_m( idof + 3 ) + new_mass_liq * (new_H * new_H + radii_v(local_idx) * radii_v(local_idx)) * 0.5;
                m( idof + 3 ) = new_inertia;
              }
            }
          });
        }
      }
    }
  });
  
  for(int k = 0; k < nk; ++k)
  {
    // combine released particle into global particles
    particles.insert(particles.end(), m_pool_liquid_particle_cache[k].begin(), m_pool_liquid_particle_cache[k].end());
    particles.insert(particles.end(), m_regular_liquid_particle_cache[k].begin(), m_regular_liquid_particle_cache[k].end());
    
    // remove liquid from hair liquid pool
    const std::vector<int>& indices = m_pool_liquid_index_cache[k];
    for(int idx : indices)
    {
      int hidx = idx / 2;
      HairFlow<3>* hair = hairs[hidx];
      hair->getPoolSize() -= m_pool_liquid_vol_cache(idx);
    }
  }
}


void FluidSim3D::shareParticleWithHairs( VectorXs& x, scalar dt )
{
  m_sorter->sort(particles.size(), sorter_callback);
  
  std::vector<HairFlow<3>*>& hairs = m_parent->getFilmFlows();
  // build bridges
  int nhair = hairs.size();
  
  const scalar theta = m_parent->getLiquidTheta();
  const scalar maxetaprop = m_parent->getMaxLimitEtaProp();
  const scalar liq_rho = m_parent->getLiquidDensity();
  
  m_bridges.resize(0);

  int np = particles.size();
  threadutils::thread_pool::ParallelFor(0, np, [&] (int i) {
    particles[i].bridges.resize(0);
  });
  m_hair_bridge_buffer.resize(nhair);
  
  VectorXs& v = m_parent->getV();
  VectorXs& m = m_parent->getM();
  
  threadutils::thread_pool::ParallelFor(0, nhair, [&] (int i) {
    HairFlow<3>* hair = hairs[i];
    auto& edges = hair->getGlobalEdges();
    const VectorXs& eta = hair->getEta();
    const VectorXs& area_e = hair->getAreaE();
    const VectorXs& radii_e = hair->getRadiiE();
    const VectorXs& radii_v = hair->getRadiiV();
    std::vector< std::vector<int> >& edge_bridges = hair->getEdgeBridges();
    edge_bridges.resize(edges.size());
    
    m_hair_bridge_buffer[i].resize(0);
    
    int ne = edge_bridges.size();

    // we rule out the first and the last edges to avoid absorbing newly generated particles
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
      const Vector3s& p0 = x.segment<3>( m_parent->getDof( e.first ) );
      const Vector3s& p1 = x.segment<3>( m_parent->getDof( e.second ) );
      
      const scalar search_dist = Hj + dx / sqrt(2.0);
      
      scalar pmin_x = std::min(p0(0) - search_dist, p1(0) - search_dist);
      scalar pmin_y = std::min(p0(1) - search_dist, p1(1) - search_dist);
      scalar pmin_z = std::min(p0(2) - search_dist, p1(2) - search_dist);
      scalar pmax_x = std::max(p0(0) + search_dist, p1(0) + search_dist);
      scalar pmax_y = std::max(p0(1) + search_dist, p1(1) + search_dist);
      scalar pmax_z = std::max(p0(2) + search_dist, p1(2) + search_dist);
      
      int imin_x = std::max(0, std::min(ni - 1, (int)( ( pmin_x - origin(0)) / dx ) ) );
      int imin_y = std::max(0, std::min(nj - 1, (int)( ( pmin_y - origin(1)) / dx ) ) );
      int imin_z = std::max(0, std::min(nk - 1, (int)( ( pmin_z - origin(2)) / dx ) ) );
      int imax_x = std::max(0, std::min(ni - 1, (int)( ceil( pmax_x - origin(0)) / dx ) ) );
      int imax_y = std::max(0, std::min(nj - 1, (int)( ceil( pmax_y - origin(1)) / dx ) ) );
      int imax_z = std::max(0, std::min(nk - 1, (int)( ceil( pmax_z - origin(2)) / dx ) ) );
      
      scalar vol_hair = M_PI * (H0 * H0 + H1 * H1 - radii_v(j) * radii_v(j) - radii_v(j + 1) * radii_v(j + 1)) * 0.5 * area_e(j);
      
      for(int t = imin_z; t <= imax_z; ++t)
      {
        for(int s = imin_y; s <= imax_y; ++s)
        {
          for(int r = imin_x; r <= imax_x; ++r)
          {
            m_sorter->getCellAt(r, s, t, [&] (int pidx) {
              Particle<3>& p = particles[pidx];
              if(p.type != PT_LIQUID) return;
              
              scalar vol_particle = dropvol(p.radii);
              
              const scalar rupture_dist = std::min(search_dist, (1.0 + 0.5 * theta) * pow(vol_particle + vol_hair, 1.0 / 3.0));
              
              scalar alpha;
              
              scalar dist = mathutils::pointedgedist(p.x, p0, p1, alpha);
              
              Vector3s pc = p0 * (1.0 - alpha) + p1 * alpha;
              
              Vector3s pp = (pc - origin) / dx - Vector3s(0.5, 0.5, 0.5);
              
              scalar phi = interpolate_value(pp, liquid_phi);
              
              if(dist < rupture_dist && phi < 0.0) {
                HairParticleBridge<3> b;
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
    }
  });

  int kidx = 0;
  for(int i = 0; i < nhair; ++i)
  {
    HairFlow<3>* hair = hairs[i];
    std::vector< std::vector<int> >& edge_bridges = hair->getEdgeBridges();
    auto& buffer = m_hair_bridge_buffer[i];
    int nb = buffer.size();
    
    for(int j = 0; j < nb; ++j)
    {
      const HairParticleBridge<3>& b = buffer[j];
      particles[b.pidx].bridges.push_back(kidx + j);
      edge_bridges[b.eidx].push_back(kidx + j);
    }
    
    m_bridges.insert(m_bridges.end(), buffer.begin(), buffer.end());
    kidx += nb;
  }

  if(m_bridges.size() == 0) return;
  
  // distribute liquid to bridges
  threadutils::thread_pool::ParallelFor(0, np, [&] (int i)
  {
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
  threadutils::thread_pool::ParallelFor(0, nhair, [&] (int i) {
    HairFlow<3>* hair = hairs[i];
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
    
    // here we use hair vertices for convenience, also avoid first and last ones
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
        if( new_modified_mass > 0.0 ) v.segment<3>( idof ) *= old_modified_mass / new_modified_mass;
        else v.segment<3>( idof ).setZero();
      }
      
      if(mum != MUM_NONE)
        m.segment<3>( idof ).setConstant(new_modified_mass);
      
      // update hair vertex angular momentum
      if(!m_parent->isMassSpring() && j != np - 1)
      {
        const scalar hair_moment_inertia = hair_mass * 0.5 * radii_v(j) * radii_v(j);
        const scalar old_liq_moment_inertia = old_modified_liq_mass * 0.5 * (modified_H * modified_H + radii_v(j) * radii_v(j));
        const scalar old_total_moment_inertia = hair_moment_inertia + old_liq_moment_inertia;
        
        const scalar new_liq_moment_inertia = new_modified_liq_mass * 0.5 * (new_modified_H * new_modified_H + radii_v(j) * radii_v(j));
        const scalar new_total_moment_inertia = hair_moment_inertia + new_liq_moment_inertia;
        
        if((mum == MUM_DIRECT_DIV || mum == MUM_MOMENTUM) && !m_parent->isFixed(pidx)) {
          if(new_total_moment_inertia > 0.0) v(idof + 3) *= old_total_moment_inertia / new_total_moment_inertia;
          else v(idof + 3) = 0.0;
        }
        
        if(mum != MUM_NONE)
          m(idof + 3) = new_total_moment_inertia;
      }
      
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
  });
  
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
  
  particles.erase(std::remove_if(particles.begin(), particles.end(), [] (const Particle<3>& p) {
    return p.type == PT_LIQUID && p.radii <= 1e-7;
  }), particles.end());
  
  m_sorter->sort(particles.size(), sorter_callback);
}

void FluidSim3D::resample(Vector3s& p, Vector3s& u, Matrix3s& c)
{
  scalar wsum = 0.0;
  Vector3s save = u;
  Matrix3s csave = c;
  u = Vector3s::Zero();
  c = Matrix3s::Zero();
  
  const scalar rho = m_parent->getLiquidDensity();
  
  int ix = std::max(0, std::min(ni - 1, (int)((p(0) - origin(0)) / dx)));
  int iy = std::max(0, std::min(nj - 1, (int)((p(1) - origin(1)) / dx)));
  int iz = std::max(0, std::min(nk - 1, (int)((p(2) - origin(2)) / dx)));
  
  const scalar re = dx;
  
  m_sorter->getNeigboringParticles_cell(ix, iy, iz, -1, 1, -1, 1, -1, 1, [&] (int pidx) {
    Particle<3>& np = particles[pidx];
    Vector3s diff = np.x - p;
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

void FluidSim3D::correct(scalar dt)
{
  int ryoichi_correction_step = m_parent->getParameter().fluidcorrectionsteps;
  
  if(!ryoichi_correction_step ) return;
  
  int np = (int) particles.size();
  
  if(!np) return;
  
#ifdef USE_MERGE_PARTICLES
  
  const scalar standard_radius = dx * default_radius_multiplier();
  const scalar maximal_radius = standard_radius * sqrt(6.0) / 2.0;
  const scalar maximal_vol = dropvol(maximal_radius);
  // merge particles
  for(Particle<3> p : particles)
  {
    p.deceased = false;
  }
  
  std::cout << "<merge particles: " << particles.size() << " -> ";
  // for each particle find its closest neighbor, ignore the ones has been deceased or with radius larger than standard
  threadutils::thread_pool::ParallelFor(0, nk, [&] (int k) {
    for(int j = 0; j < nj; ++j) for(int i = 0; i < ni; ++i) {
      m_sorter->getCellAt(i, j, k, [&] (int pidx_i) {
        if(particles[pidx_i].type != PT_LIQUID || particles[pidx_i].deceased || particles[pidx_i].pressure * dx < 10.0) return;
        scalar vol_i = dropvol(particles[pidx_i].radii);

        m_sorter->getCellAt(i, j, k, [&] (int pidx_j) {
          if(particles[pidx_j].type != PT_LIQUID || pidx_j <= pidx_i || particles[pidx_j].deceased || particles[pidx_j].pressure * dx < 10.0) return;

          // for each pair, check their add-up volume
          scalar vol_j = dropvol(particles[pidx_j].radii);
          scalar combined_vol = vol_i + vol_j;
          if(combined_vol > maximal_vol) return;
   
          scalar dist = (particles[pidx_j].x - particles[pidx_i].x).norm();
          if(dist > particles[pidx_i].radii + particles[pidx_j].radii) return;

          // combine if their distance < add-up radius and add-up volume < maximal volume
          particles[pidx_j].deceased = true;

          scalar alpha = vol_j / combined_vol;
          particles[pidx_i].x = mathutils::lerp(particles[pidx_i].x, particles[pidx_j].x, alpha);
          particles[pidx_i].v = mathutils::lerp(particles[pidx_i].v, particles[pidx_j].v, alpha);
          particles[pidx_i].c = mathutils::lerp(particles[pidx_i].c, particles[pidx_j].c, alpha);
          
          vol_i = combined_vol;
        });
        
        particles[pidx_i].radii = dropradius(vol_i);
      });
    }
  });
  
  particles.erase(std::remove_if(particles.begin(), particles.end(), [] (const Particle<3>& p) {
    return p.deceased;
  }), particles.end());
  
  m_sorter->sort(particles.size(), sorter_callback);
  std::cout << particles.size() << ">\n";
#endif
  scalar coeff = pow(ni * nj * nk, 1.0 / 3.0);
  
  threadutils::thread_pool::ParallelFor(0, np, [&] (int n)
  {
    Particle<3> &p = particles[n];
    
    if(p.type != PT_LIQUID) return;
    
    if(n % ryoichi_correction_step != ryoichi_correction_counter % ryoichi_correction_step) {
      return;
    }
    
    Vector3s spring = Vector3s::Zero();
    
    int ix = std::max(0, std::min((int)((p.x(0) - origin(0)) / dx), ni));
    int iy = std::max(0, std::min((int)((p.x(1) - origin(1)) / dx), nj));
    int iz = std::max(0, std::min((int)((p.x(2) - origin(2)) / dx), nk));
    
    m_sorter->getNeigboringParticles_cell(ix, iy, iz, -1, 1, -1, 1, -1, 1, [&] (int pidx) {
      const Particle<3> &np = particles[pidx];
      if(n != pidx)
      {
        scalar re = sqrt(p.radii * np.radii);
        scalar dist = (p.x - np.x).norm();
        scalar w = coeff * mathutils::smooth_kernel(dist * dist, re);
        if( dist > 1e-4 * re )
        {
          spring += w * (p.x - np.x) / dist * re;
        } else {
          spring(0) += re * (rand() & 0xFF) / 255.0;
          spring(1) += re * (rand() & 0xFF) / 255.0;
          spring(2) += re * (rand() & 0xFF) / 255.0;
        }
      }
    });
    
    p.buf0 = p.x + dt * spring;
    
    Vector3s pp = (p.buf0 - origin)/dx;
    scalar phi_value = interpolate_value(pp, nodal_solid_phi);
    if(phi_value < 0) {
      Vector3s normal;
      interpolate_gradient(normal, pp, nodal_solid_phi);
      normal.normalize();
      p.buf0 -= phi_value*normal;
    }
  });
  

  // Update
  threadutils::thread_pool::ParallelFor(0, np, [&] (int n)
  {
    Particle<3>& p = particles[n];
    if(p.type != PT_LIQUID) return;
    
    if(n % ryoichi_correction_step != ryoichi_correction_counter % ryoichi_correction_step) {
      return;
    }
    
    p.x = p.buf0;
  });

  m_sorter->sort(particles.size(), sorter_callback);
  
  ryoichi_correction_counter++;
}

void FluidSim3D::add_drag(scalar dt)
{
  // drag
  for(int k = 0; k < u.nk; ++k) for(int j = 0; j < u.nj; ++j) for(int i = 0; i < u.ni; ++i)
  {
    u_particle(i, j, k) += u_drag(i, j, k) * dt;
  }
  
  for(int k = 0; k < v.nk; ++k) for(int j = 0; j < v.nj; ++j) for(int i = 0; i < v.ni; ++i)
  {
    v_particle(i, j, k) += v_drag(i, j, k) * dt;
  }
  
  for(int k = 0; k < w.nk; ++k) for(int j = 0; j < w.nj; ++j) for(int i = 0; i < w.ni; ++i)
  {
    w_particle(i, j, k) += w_drag(i, j, k) * dt;
  }
}

void FluidSim3D::add_gravity(scalar dt)
{
  const Vector3s& gravity = m_parent->getSimpleGravity();
  
  if(m_parent->applyCoriolis()) {
    const scalar& lat = m_parent->getLatitude();
    const Vector3s omega = Vector3s(-sin(lat), cos(lat), 0.0) * m_parent->getEarthRotation();
    const Vector3s earth_R = Vector3s(0.0, m_parent->getEarthRadius(), 0.0);
    
    temp_u = u_particle;
    temp_v = v_particle;
    temp_w = w_particle;
    
    threadutils::thread_pool::ParallelFor(0, u.nk, [&] (int k){
      for(int j = 0; j < u.nj; ++j) for(int i = 0; i < u.ni; ++i)
      {
        Vector3s pos = Vector3s(i*dx, (j+0.5)*dx, (k+0.5)*dx) + origin;
        
        Vector3s vel = get_temp_velocity(pos);
        
        Vector3s r0 = earth_R + Vector3s(0.0, pos(1), 0.0);
        Vector3s acc_centri = omega.cross(omega.cross(r0));
        Vector3s acc_coriolis = 2.0 * omega.cross(vel);
        
        u_particle(i, j, k) -= (-gravity(0) + acc_coriolis(0) + acc_centri(0)) * dt;
      }
    });
    
    threadutils::thread_pool::ParallelFor(0, v.nk, [&] (int k){
      for(int j = 0; j < v.nj; ++j) for(int i = 0; i < v.ni; ++i)
      {
        Vector3s pos = Vector3s((i+0.5)*dx, j*dx, (k+0.5)*dx) + origin;
        
        Vector3s vel = get_temp_velocity(pos);
        
        Vector3s r0 = earth_R + Vector3s(0.0, pos(1), 0.0);
        Vector3s acc_centri = omega.cross(omega.cross(r0));
        Vector3s acc_coriolis = 2.0 * omega.cross(vel);
        
        v_particle(i, j, k) -= (-gravity(1) + acc_coriolis(1) + acc_centri(1)) * dt;
      }
    });
    
    threadutils::thread_pool::ParallelFor(0, w.nk, [&] (int k){
      for(int j = 0; j < w.nj; ++j) for(int i = 0; i < w.ni; ++i)
      {
        Vector3s pos = Vector3s((i+0.5)*dx, (j+0.5)*dx, k*dx) + origin;
        
        Vector3s vel = get_temp_velocity(pos);
        
        Vector3s r0 = earth_R + Vector3s(0.0, pos(1), 0.0);
        Vector3s acc_centri = omega.cross(omega.cross(r0));
        Vector3s acc_coriolis = 2.0 * omega.cross(vel);
        
        w_particle(i, j, k) -= (-gravity(2) + acc_coriolis(2) + acc_centri(2)) * dt;
      }
    });
  } else {
    threadutils::thread_pool::ParallelFor(0, u.nk, [&] (int k){
      for(int j = 0; j < u.nj; ++j) for(int i = 0; i < u.ni; ++i)
      {
        u_particle(i, j, k) -= (-gravity(0)) * dt;
      }
    });
    
    threadutils::thread_pool::ParallelFor(0, v.nk, [&] (int k){
      for(int j = 0; j < v.nj; ++j) for(int i = 0; i < v.ni; ++i)
      {
        v_particle(i, j, k) -= (-gravity(1)) * dt;
      }
    });
    
    threadutils::thread_pool::ParallelFor(0, w.nk, [&] (int k){
      for(int j = 0; j < w.nj; ++j) for(int i = 0; i < w.ni; ++i)
      {
        w_particle(i, j, k) -= (-gravity(2)) * dt;
      }
    });
  }
}

//For extrapolated points, replace the normal component
//of velocity with the object velocity (in this case zero).
void FluidSim3D::constrain_velocity() {
  temp_u = u;
  temp_v = v;
  temp_w = w;
  
  //(At lower grid resolutions, the normal estimate from the signed
  //distance function is poor, so it doesn't work quite as well.
  //An exact normal would do better.)
  
  //constrain u
  threadutils::thread_pool::ParallelFor(0, u.nk, [&] (int k){
    for(int j = 0; j < u.nj; ++j) for(int i = 0; i < u.ni; ++i) {
      if(u_weights(i,j,k) == 0) {
        //apply constraint
        Vector3s pos = Vector3s(i*dx, (j+0.5)*dx, (k+0.5)*dx) + origin;
        Vector3s vel = get_velocity(pos);
        Vector3s normal(0,0,0);
        interpolate_gradient(normal, Vector3s(i, j + 0.5, k + 0.5), nodal_solid_phi);
        normal.normalize();
        scalar perp_component = vel.dot(normal);
        vel -= perp_component*normal;
        Vector3s vel_sol = get_solid_velocity(pos);
        vel += vel_sol.dot(normal) * normal;
        temp_u(i,j,k) = vel[0];
      }
    }
  });
  
  //constrain v
  threadutils::thread_pool::ParallelFor(0, v.nk, [&] (int k){
    for(int j = 0; j < v.nj; ++j) for(int i = 0; i < v.ni; ++i) {
      if(v_weights(i,j,k) == 0) {
        //apply constraint
        Vector3s pos = Vector3s((i+0.5)*dx, j*dx, (k+0.5)*dx) + origin;
        Vector3s vel = get_velocity(pos);
        Vector3s normal(0,0,0);
        interpolate_gradient(normal, Vector3s(i + 0.5, j, k + 0.5), nodal_solid_phi);
        normal.normalize();
        scalar perp_component = vel.dot(normal);
        vel -= perp_component*normal;
        Vector3s vel_sol = get_solid_velocity(pos);
        vel += vel_sol.dot(normal) * normal;
        temp_v(i,j,k) = vel[1];
      }
    }
  });
  
  //constrain w
  threadutils::thread_pool::ParallelFor(0, w.nk, [&] (int k){
    for(int j = 0; j < w.nj; ++j) for(int i = 0; i < w.ni; ++i) {
      if(w_weights(i,j,k) == 0) {
        //apply constraint
        Vector3s pos = Vector3s((i+0.5)*dx, (j+0.5)*dx, k*dx) + origin;
        Vector3s vel = get_velocity(pos);
        Vector3s normal(0,0,0);
        interpolate_gradient(normal, Vector3s(i + 0.5, j + 0.5, k), nodal_solid_phi);
        normal.normalize();
        scalar perp_component = vel.dot(normal);
        vel -= perp_component*normal;
        Vector3s vel_sol = get_solid_velocity(pos);
        vel += vel_sol.dot(normal) * normal;
        temp_w(i,j,k) = vel[2];
      }
    }
  });
  
  //update
  u = temp_u;
  v = temp_v;
  w = temp_w;
}



//Add a tracer particle for visualization
void FluidSim3D::add_particle(const Particle<3>& p) {
  particles.push_back(p);
}

//move the particles in the fluid
void FluidSim3D::advect_particles(scalar dt) {
  int np = particles.size();

  const scalar air_viscosity = m_parent->getAirViscosity();
  const scalar default_radius = dx * default_radius_multiplier();
  if(air_viscosity > 0.0) {
    const scalar sirignano = 0.72;
    const scalar air_density = 0.0012;
    const scalar coeff = 9.0 / pow(2.0, sirignano) * air_density / m_parent->getLiquidDensity() * pow(air_viscosity, sirignano);

    threadutils::thread_pool::ParallelFor(0, np, [&](int p){
      // symplectic integrate air-drag
      if(particles[p].radii < 0.5 * default_radius) {
        const scalar beta_w = coeff * pow(particles[p].radii, -(1.0 + sirignano));
        const scalar airdrag = beta_w * pow(particles[p].v.norm(), 1.0 - sirignano);
        particles[p].v /= (1.0 + dt * airdrag);
        particles[p].c /= (1.0 + dt * airdrag);
      }
      
      particles[p].x += particles[p].v * dt;
      Vector3s pp = (particles[p].x - origin)/dx;
      
      //Particles can still occasionally leave the domain due to truncation errors,
      //interpolation error, or large timesteps, so we project them back in for good measure.
      
      //Try commenting this section out to see the degree of accumulated error.
      scalar phi_value = interpolate_value(pp, nodal_solid_phi);
      if(phi_value < 0) {
        Vector3s normal;
        interpolate_gradient(normal, pp, nodal_solid_phi);
        normal.normalize();
        particles[p].x -= phi_value*normal;
      }
    });
  } else {
    threadutils::thread_pool::ParallelFor(0, np, [&](int p){
      particles[p].x += particles[p].v * dt;
      Vector3s pp = (particles[p].x - origin)/dx;
      
      //Particles can still occasionally leave the domain due to truncation errors,
      //interpolation error, or large timesteps, so we project them back in for good measure.
      
      //Try commenting this section out to see the degree of accumulated error.
      scalar phi_value = interpolate_value(pp, nodal_solid_phi);
      if(phi_value < 0) {
        Vector3s normal;
        interpolate_gradient(normal, pp, nodal_solid_phi);
        normal.normalize();
        particles[p].x -= phi_value*normal;
      }
    });
  }
  
  
  sort_particles();
  
  scalar inserted_volume = 0.0;

  for(SourceBoundary<3>* s : sources)
  {
    if(!s->activated) continue;
    
    scalar generate_radius = s->drop_radius_prop * default_radius;
    scalar generate_vol = dropvol(generate_radius);
    
    
    for(const Vector3s& p : s->detectors)
    {
      if(p(0) < origin(0) || p(0) > origin(0) + ((scalar) ni+1.) * dx
         || p(1) < origin(1) || p(1) > origin(1) + ((scalar) nj+1.) * dx
         || p(2) < origin(2) || p(2) > origin(2) + ((scalar) nk+1.) * dx) continue;
         
      Vector3s pp = (p - origin)/dx;
      
      int ix = (int) pp(0);
      int iy = (int) pp(1);
      int iz = (int) pp(2);
      
      scalar sum_vol_existed = 0.0;
      m_sorter->getCellAt(ix, iy, iz, [&] (int i) {
        sum_vol_existed += dropvol(particles[i].radii);
      });
      
      int num_p_need = default_particle_in_cell();//(int)((default_vol - sum_vol_existed) / generate_vol);
      
      Vector3s cone_dir = s->eject_vel.normalized();
      Vector3s q = -cone_dir.cross(Vector3s(0, 0, 1));
      scalar angle = atan2(q.norm(), cone_dir.dot(Vector3s(0, 0, 1)));
      
      Matrix3s rot_mat = Eigen::AngleAxis<scalar>(angle, q).toRotationMatrix();
      
      for(int k = 0; k < num_p_need; ++k)
      {
        scalar x = ((scalar) ix + 0.5 + (((scalar)rand() / (scalar)RAND_MAX) * 0.5 - 0.5) ) * dx;
        scalar y = ((scalar) iy + 0.5 + (((scalar)rand() / (scalar)RAND_MAX) * 0.5 - 0.5) ) * dx;
        scalar z = ((scalar) iz + 0.5 + (((scalar)rand() / (scalar)RAND_MAX) * 0.5 - 0.5) ) * dx;
        Vector3s pt = Vector3s(x, y, z) + origin;
        
        scalar phi_value = interpolate_value(pp, nodal_solid_phi);
        if(phi_value > 0.0) {
          Vector3s vel = Vector3s::Zero();
          s->compute_phi_vel(pt, vel);
          Vector3s ej_vel = s->eject_vel;
          if(s->spray_angle > 0.0 && s->eject_vel.norm() > 0.0) {
            scalar z = (scalar) rand() / (scalar) RAND_MAX * (1.0 - cos(s->spray_angle)) + cos(s->spray_angle);
            scalar phi = (scalar) rand() / (scalar) RAND_MAX * 2.0 * M_PI;
            
            Vector3s vz = Vector3s(sqrt(1.0 - z * z) * cos(phi), sqrt(1.0 - z * z) * sin(phi), z).normalized();
            if(vz(0) == 0.0 && vz(1) == 0.0 && vz(2) == -1.0) {
              vz = -vz;
            } else if(!(vz(0) == 0.0 && vz(1) == 0.0 && vz(2) == 1.0)) {
              vz = rot_mat * vz;
            }
/*            Vector3s vz = Vector3s((scalar) rand() / (scalar) RAND_MAX * 2.0 - 1.0,
                                   (scalar) rand() / (scalar) RAND_MAX * 2.0 - 1.0,
                                   (scalar) rand() / (scalar) RAND_MAX * 2.0 - 1.0).normalized();*/
            
            ej_vel = vz * (s->eject_vel.norm() * ((scalar) rand() / (scalar) RAND_MAX * 0.5 + 0.5));
          }
          
          vel += ej_vel;
          add_particle(Particle<3>(pt + vel * dt, vel, generate_radius, PT_LIQUID));
          inserted_volume += generate_vol;
        }
      }
    }
  }
  
  scalar volume_removed = 0.0;
  
  particles.erase( std::remove_if(particles.begin(), particles.end(), [&] (const Particle<3>& p) {
    bool removed =
    p.x(0) < origin(0) + 0.5 * dx || p.x(0) > origin(0) + ((scalar) ni+0.5) * dx
    || p.x(1) < origin(1) + 0.5 * dx || p.x(1) > origin(1) + ((scalar) nj+0.5) * dx
    || p.x(2) < origin(2) + 0.5 * dx || p.x(2) > origin(2) + ((scalar) nk+0.5) * dx;
    
    if(removed) {
      volume_removed += dropvol(p.radii);
    }
    return removed;
  }), particles.end());
  
  m_parent->reportParticleAdded(inserted_volume);
  m_parent->reportParticleRemoved(volume_removed);
  
  sort_particles();
}

void FluidSim3D::constrain_hair_particles()
{
  int np = particles.size();
  
  const VectorXs& x = m_parent->getX();
  const std::vector< std::pair<int, int> >& edges = m_parent->getEdges();
  const std::vector< HairFlow<3>* >& flows = m_parent->getFilmFlows();
  const std::vector< int >& particle_hairs = m_parent->getParticleToHairs();
  const std::vector< int >& local_indices = m_parent->getParticleToHairLocalIndices();
  
  threadutils::thread_pool::ParallelFor(0, np, [&](int p){
    if(particles[p].type == PT_HAIR) {
      auto& e = edges[particles[p].edge_idx];
      const HairFlow<3>* flow = flows[ particle_hairs[e.first] ];
      const VectorXs& eta = flow->getEta();
      const VectorXs& radii_v = flow->getRadiiV();
      
      int local_idx0 = local_indices[e.first];
      int local_idx1 = local_indices[e.second];
      
      const Vector3s& x0 = x.segment<3>(m_parent->getDof(e.first));
      const Vector3s& x1 = x.segment<3>(m_parent->getDof(e.second));
      particles[p].x = mathutils::lerp(x0, x1, particles[p].edge_alpha);
      
      scalar H0 = eta(local_idx0) + radii_v(local_idx0);
      scalar H1 = eta(local_idx1) + radii_v(local_idx1);
      
      particles[p].radii = sqrt(mathutils::lerp(H0 * H0, H1 * H1, particles[p].edge_alpha));
    }
  });
  
  m_sorter->sort(particles.size(), sorter_callback);
}

void FluidSim3D::compute_liquid_phi()
{
  liquid_phi.assign(3*dx);
  std::cout << "[Liquid-Phi: CVT]" << std::endl;

  threadutils::thread_pool::ParallelFor(0, nk, [&] (int k){
    for(int j = 0; j < nj; ++j) for(int i = 0; i < ni; ++i) {
      Vector3s pos = Vector3s((i+0.5)*dx, (j+0.5)*dx, (k+0.5)*dx) + origin;
      scalar phi_min = 1e+20;
      
      m_sorter->getNeigboringParticles_cell(i, j, k, -2, 2, -2, 2, -2, 2, [&] (int pidx) {
        const Particle<3>& p = particles[pidx];
        phi_min = std::min(phi_min, (pos - p.x).norm() - std::max(dx * 0.883644, p.radii));
      });
      liquid_phi(i,j,k) = std::min(liquid_phi(i,j,k), phi_min);
    }
  });
  
  std::cout << "[Liquid-Phi: extrapolate phi into solids]" << std::endl;
  threadutils::thread_pool::ParallelFor(0, nk, [&](int k){
    for(int j = 0; j < nj; ++j) {
      for(int i = 0; i < ni; ++i) {
          Vector3s pos = Vector3s((i+0.5)*dx, (j+0.5)*dx, (k+0.5)*dx) + origin;
          Vector3s vel;
          scalar solid_phi_val = compute_phi_vel(pos, vel);
          liquid_phi(i,j,k) = std::min(liquid_phi(i,j,k), solid_phi_val);
      }
    }
  });
  
  //write_matlab_array(std::cout, liquid_phi, "phi");
}

void FluidSim3D::combine_velocity_field()
{
  // Combine the velocity field only for pressure computation later,
  // we replace the penalty method in [Bandara & Soga 2015] with a Poisson solver.
  
  threadutils::thread_pool::ParallelFor(0, u.nk, [&](int k){
    for(int j = 0; j < u.nj; ++j) for(int i = 0; i < u.ni; ++i) {
      scalar weight = u_weight_particle(i, j, k) + u_weight_hair(i, j, k);
      if(weight > 0) {
        u(i, j, k) = (u_weight_particle(i, j, k) * u_particle(i, j, k) + u_weight_hair(i, j, k) * u_hair(i, j, k)) / weight;
      } else {
        u(i, j, k) = 0.0;
      }
      u_weight_total(i, j, k) = weight;
    }
  });
  
  threadutils::thread_pool::ParallelFor(0, v.nk, [&](int k){
    for(int j = 0; j < v.nj; ++j) for(int i = 0; i < v.ni; ++i) {
      scalar weight = v_weight_particle(i, j, k) + v_weight_hair(i, j, k);
      if(weight > 0) {
        v(i, j, k) = (v_weight_particle(i, j, k) * v_particle(i, j, k) + v_weight_hair(i, j, k) * v_hair(i, j, k)) / weight;
      } else {
        v(i, j, k) = 0.0;
      }
      v_weight_total(i, j, k) = weight;
    }
  });
  
  threadutils::thread_pool::ParallelFor(0, w.nk, [&](int k){
    for(int j = 0; j < w.nj; ++j) for(int i = 0; i < w.ni; ++i) {
      scalar weight = w_weight_particle(i, j, k) + w_weight_hair(i, j, k);
      if(weight > 0) {
        w(i, j, k) = (w_weight_particle(i, j, k) * w_particle(i, j, k) + w_weight_hair(i, j, k) * w_hair(i, j, k)) / weight;
      } else {
        w(i, j, k) = 0.0;
      }
      w_weight_total(i, j, k) = weight;
    }
  });
}

void FluidSim3D::project(scalar dt) {
  compute_weights();
  
  //Set up and solve the variational pressure solve.
  
  solve_pressure(dt);
  
  temp_u = u;
  temp_v = v;
  temp_w = w;
  
  old_valid = valid;
  
  extrapolate(u, temp_u, u_weights, liquid_phi, valid, old_valid, Vector3i(-1, 0, 0));
  extrapolate(v, temp_v, v_weights, liquid_phi, valid, old_valid, Vector3i(0, -1, 0));
  extrapolate(w, temp_w, w_weights, liquid_phi, valid, old_valid, Vector3i(0, 0, -1));
}

scalar FluidSim3D::get_nodal_solid_phi(const Vector3s& position) const
{
  Vector3s p = (position - origin) / dx;
  return interpolate_value(p, nodal_solid_phi);
}

Vector3s FluidSim3D::get_nodal_solid_phi_gradient(const Vector3s& position) const
{
  Vector3s p = (position - origin) / dx;
  Vector3s normal;
  interpolate_gradient(normal, p, nodal_solid_phi);
  return normal;
}

Vector3s FluidSim3D::get_velocity(const Vector3s& position, const Array3s& u_, const Array3s& v_, const Array3s& w_) const
{
  //Interpolate the velocity from the u and v grids
  Vector3s p = (position - origin) / dx;
  Vector3s p0 = p - Vector3s(0, 0.5, 0.5);
  Vector3s p1 = p - Vector3s(0.5, 0, 0.5);
  Vector3s p2 = p - Vector3s(0.5, 0.5, 0);
  scalar u_value = interpolate_value(p0, u_);
  scalar v_value = interpolate_value(p1, v_);
  scalar w_value = interpolate_value(p2, w_);
  
  return Vector3s(u_value, v_value, w_value);
}

Vector3s FluidSim3D::get_visc_impulse(const Vector3s& position) const
{
  return get_velocity(position, u_visc_impulse, v_visc_impulse, w_visc_impulse);
}

Vector3s FluidSim3D::get_pressure_gradient(const Vector3s& position) const
{
  return get_velocity(position, u_pressure_grad, v_pressure_grad, w_pressure_grad);
}

scalar FluidSim3D::get_pressure(const Vector3s& position) const
{
  Vector3s p = (position - origin) / dx - Vector3s(0.5, 0.5, 0.5);
  return interpolate_value(p, pressure, ni, nj, nk);
}

//Interpolate velocity from the MAC grid.
Vector3s FluidSim3D::get_velocity(const Vector3s& position) const {
  return get_velocity(position, u, v, w);
}

//Interpolate solid velocity from the MAC grid.
Vector3s FluidSim3D::get_solid_velocity(const Vector3s& position) const {
  return get_velocity(position, u_solid, v_solid, w_solid);
}

//Interpolate drag from the MAC grid.
Vector3s FluidSim3D::get_particle_drag(const Vector3s& position) const {
  return get_velocity(position, u_drag, v_drag, w_drag);
}

//Interpolate hair velocity from the MAC grid.
Vector3s FluidSim3D::get_hair_velocity(const Vector3s& position) const {
  return get_velocity(position, u_hair, v_hair, w_hair);
}

//Interpolate particle velocity from the MAC grid.
Vector3s FluidSim3D::get_particle_velocity(const Vector3s& position) const {
  return get_velocity(position, u_particle, v_particle, w_particle);
}

//Interpolate particle velocity from the MAC grid.
Vector3s FluidSim3D::get_temp_velocity(const Vector3s& position) const {
  return get_velocity(position, temp_u, temp_v, temp_w);
}

scalar FluidSim3D::getLiquidPhiValue(const Vector3s& position) const
{
  Vector3s pp = (position - origin) / dx - Vector3s(0.5, 0.5, 0.5);
  return interpolate_value(pp, liquid_phi);
}

scalar FluidSim3D::getClampedLiquidPhiValue(const Vector3s& position) const
{
  const scalar standard_radius = dx * default_radius_multiplier();
  scalar w = getLiquidPhiValue(position);
  scalar criterion = -m_parent->getBulkThresholdMultiplier() * standard_radius;
  return mathutils::clamp(w / criterion, 0.0, 1.0);
}

Matrix3s FluidSim3D::get_affine_matrix(const Vector3s& position) const
{
  Vector3s p = (position - origin) / dx;
  Vector3s p0 = p - Vector3s(0, 0.5, 0.5);
  Vector3s p1 = p - Vector3s(0.5, 0, 0.5);
  Vector3s p2 = p - Vector3s(0.5, 0.5, 0);
  
  Matrix3s c;
  c.col(0) = affine_interpolate_value(p0, u) / dx;
  c.col(1) = affine_interpolate_value(p1, v) / dx;
  c.col(2) = affine_interpolate_value(p2, w) / dx;
  
  return c;
}

//Compute finite-volume style face-weights for fluid from nodal signed distances
void FluidSim3D::compute_weights() {
  threadutils::thread_pool::ParallelFor(0, u_weights.nk, [&] (int k) {
    for(int j = 0; j < u_weights.nj; ++j) for(int i = 0; i < u_weights.ni; ++i) {
      u_weights(i,j,k) = 1 - mathutils::fraction_inside(nodal_solid_phi(i,j,k), nodal_solid_phi(i,j+1,k), nodal_solid_phi(i,j,k+1), nodal_solid_phi(i,j+1,k+1));
      u_weights(i,j,k) = hardclamp(u_weights(i,j,k), 0.0, 1.0);
    }
  });
  threadutils::thread_pool::ParallelFor(0, v_weights.nk, [&] (int k) {
    for(int j = 0; j < v_weights.nj; ++j) for(int i = 0; i < v_weights.ni; ++i) {
      v_weights(i,j,k) = 1 - mathutils::fraction_inside(nodal_solid_phi(i,j,k), nodal_solid_phi(i+1,j,k), nodal_solid_phi(i,j,k+1), nodal_solid_phi(i+1,j,k+1));
      v_weights(i,j,k) = hardclamp(v_weights(i,j,k), 0.0, 1.0);
    }
  });
  threadutils::thread_pool::ParallelFor(0, w_weights.nk, [&] (int k) {
    for(int j = 0; j < w_weights.nj; ++j) for(int i = 0; i < w_weights.ni; ++i) {
      w_weights(i,j,k) = 1 - mathutils::fraction_inside(nodal_solid_phi(i,j,k), nodal_solid_phi(i,j+1,k), nodal_solid_phi(i+1,j,k), nodal_solid_phi(i+1,j+1,k));
      w_weights(i,j,k) = hardclamp(w_weights(i,j,k), 0.0, 1.0);
    }
  });
}

//#define USE_ROBERTS_SOLVER

//An implementation of the variational pressure projection solve for static geometry
void FluidSim3D::solve_pressure(scalar dt) {
  //This linear system could be simplified, but I've left it as is for clarity
  //and consistency with the standard naive discretization
  
  int ni = v.ni;
  int nj = u.nj;
  int nk = u.nk;
  
  int system_size = ni*nj*nk;
  int slice = ni*nj;
  if((int) pressure.size() != system_size) {
    pressure.resize(system_size);
  }
  
  pressure.assign(pressure.size(), 0);
  
  std::vector<double> x;
  std::vector<Vector3i> dof_ijk;

  dof_ijk.resize(0);
  dof_ijk.reserve(system_size);
  
  Array3ui dof_index;
  dof_index.resize(ni,nj,nk);
  dof_index.assign(0);
  
  const scalar rho = m_parent->getLiquidDensity();
  
  for (int k=1;k<nk-1;k++)for(int j=1;j<nj-1;j++)for(int i=1;i<ni-1;i++)
  {
    if (liquid_phi(i,j,k)<0 && (u_weights(i, j, k) > 1e-12 || u_weights(i+1, j, k) > 1e-12 ||
                                v_weights(i, j, k) > 1e-12 || v_weights(i, j+1, k) > 1e-12 ||
                                w_weights(i, j, k) > 1e-12 || w_weights(i, j, k+1) > 1e-12))
    {
      dof_index(i,j,k) = dof_ijk.size();
      dof_ijk.push_back(Vector3i(i,j,k));
    }
  }
  
  x.resize(dof_ijk.size());
  x.assign(dof_ijk.size(), 0);
    
  if(rhs.size() != x.size()) {
    rhs.resize(x.size());
    matrix.resize(x.size());
  }
  
  rhs.assign(rhs.size(), 0);
  matrix.zero();
  
  threadutils::thread_pool::ParallelFor(0, (int) dof_ijk.size(), [&] (int dof_idx) {
    const int i = dof_ijk[dof_idx](0);
    const int j = dof_ijk[dof_idx](1);
    const int k = dof_ijk[dof_idx](2);
    
    scalar centre_phi = liquid_phi(i,j,k);
    rhs[dof_idx] = 0;
    //right neighbour
    scalar term = u_weights(i+1,j,k) * dt / sqr(dx) / rho;
    scalar right_phi = liquid_phi(i+1,j,k);
    if(right_phi < 0) {
      matrix.add_to_element(dof_idx,
                            dof_idx, term);
      matrix.add_to_element(dof_idx,
                            dof_index(i+1,j,k), -term);
    }
    else {
      scalar theta = fraction_inside(centre_phi, right_phi);
      if(theta < 0.01) theta = 0.01;
      matrix.add_to_element(dof_idx,
                            dof_idx, term/theta);
    }
    rhs[dof_idx] -= (u_weights(i+1,j,k)*u(i+1,j,k) + (1.0-u_weights(i+1,j,k)) * u_solid(i+1,j,k)) / dx;
    
    //left neighbour
    term = u_weights(i,j,k) * dt / sqr(dx) / rho;
    scalar left_phi = liquid_phi(i-1,j,k);
    if(left_phi < 0) {
      matrix.add_to_element(dof_idx,
                            dof_idx, term);
      matrix.add_to_element(dof_idx,
                            dof_index(i-1,j,k), -term);
    }
    else {
      scalar theta = fraction_inside(centre_phi, left_phi);
      if(theta < 0.01) theta = 0.01;
      matrix.add_to_element(dof_idx,
                            dof_idx, term/theta);
    }
    rhs[dof_idx] += (u_weights(i,j,k)*u(i,j,k) + (1.0-u_weights(i,j,k)) * u_solid(i,j,k)) / dx;
    
    //top neighbour
    term = v_weights(i,j+1,k) * dt / sqr(dx) / rho;
    scalar top_phi = liquid_phi(i,j+1,k);
    if(top_phi < 0) {
      matrix.add_to_element(dof_idx,
                            dof_idx, term);
      matrix.add_to_element(dof_idx,
                            dof_index(i,j+1,k), -term);
    }
    else {
      scalar theta = fraction_inside(centre_phi, top_phi);
      if(theta < 0.01) theta = 0.01;
      matrix.add_to_element(dof_idx,
                            dof_idx, term/theta);
    }
    rhs[dof_idx] -= (v_weights(i,j+1,k)*v(i,j+1,k) + (1.0-v_weights(i,j+1,k)) * v_solid(i,j+1,k)) / dx;
    
    //bottom neighbour
    term = v_weights(i,j,k) * dt / sqr(dx) / rho;
    scalar bot_phi = liquid_phi(i,j-1,k);
    if(bot_phi < 0) {
      matrix.add_to_element(dof_idx,
                            dof_idx, term);
      matrix.add_to_element(dof_idx,
                            dof_index(i,j-1,k), -term);
    }
    else {
      scalar theta = fraction_inside(centre_phi, bot_phi);
      if(theta < 0.01) theta = 0.01;
      matrix.add_to_element(dof_idx,
                            dof_idx, term/theta);
    }
    rhs[dof_idx] += (v_weights(i,j,k)*v(i,j,k) + (1.0-v_weights(i,j,k)) * v_solid(i,j,k)) / dx;
    
    
    //far neighbour
    term = w_weights(i,j,k+1) * dt / sqr(dx) / rho;
    scalar far_phi = liquid_phi(i,j,k+1);
    if(far_phi < 0) {
      matrix.add_to_element(dof_idx,
                            dof_idx, term);
      matrix.add_to_element(dof_idx,
                            dof_index(i,j,k+1), -term);
    }
    else {
      scalar theta = fraction_inside(centre_phi, far_phi);
      if(theta < 0.01) theta = 0.01;
      matrix.add_to_element(dof_idx,
                            dof_idx, term/theta);
    }
    rhs[dof_idx] -= (w_weights(i,j,k+1)*w(i,j,k+1) + (1.0-w_weights(i,j,k+1)) * w_solid(i,j,k+1)) / dx;
    
    //near neighbour
    term = w_weights(i,j,k) * dt / sqr(dx) / rho;
    scalar near_phi = liquid_phi(i,j,k-1);
    if(near_phi < 0) {
      matrix.add_to_element(dof_idx,
                            dof_idx, term);
      matrix.add_to_element(dof_idx,
                            dof_index(i,j,k-1), -term);
    }
    else {
      scalar theta = fraction_inside(centre_phi, near_phi);
      if(theta < 0.01) theta = 0.01;
      matrix.add_to_element(dof_idx,
                            dof_idx, term/theta);
    }
    rhs[dof_idx] += (w_weights(i,j,k)*w(i,j,k) + (1.0-w_weights(i,j,k)) * w_solid(i,j,k)) / dx;
  });

  //Solve the system using Robert Bridson's incomplete Cholesky PCG solver
  
  scalar tolerance;
  int iterations;
  
  bool success = false;
#ifdef USE_ROBERTS_SOLVER
  solver.set_solver_parameters(1e-18, 1000);
  success = solver.solve(matrix, rhs, pressure, tolerance, iterations);
#else
  success = AMGPCGSolveSparse(matrix,rhs,x,dof_ijk,1e-6,500,tolerance,iterations,ni,nj,nk);
  
#endif

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
    write_matlab_array(std::cout, w, "w");

    write_matlab_array(std::cout, u_solid, "u_solid");
    write_matlab_array(std::cout, v_solid, "v_solid");
    write_matlab_array(std::cout, w_solid, "w_solid");
    
    write_matlab_array(std::cout, u_hair, "u_hair");
    write_matlab_array(std::cout, v_hair, "v_hair");
    write_matlab_array(std::cout, w_hair, "w_hair");
    
    write_matlab_array(std::cout, u_particle, "u_particle");
    write_matlab_array(std::cout, v_particle, "v_particle");
    write_matlab_array(std::cout, w_particle, "w_particle");
    
    write_matlab_array(std::cout, u_drag, "u_drag");
    write_matlab_array(std::cout, v_drag, "v_drag");
    write_matlab_array(std::cout, w_drag, "w_drag");
    
    write_matlab_array(std::cout, u_weights, "u_weights");
    write_matlab_array(std::cout, v_weights, "v_weights");
    write_matlab_array(std::cout, w_weights, "w_weights");
   
    write_matlab_array(std::cout, u_weight_particle, "u_weight_particle");
    write_matlab_array(std::cout, v_weight_particle, "v_weight_particle");
    write_matlab_array(std::cout, w_weight_particle, "w_weight_particle");
    
    write_matlab_array(std::cout, u_weight_hair, "u_weight_hair");
    write_matlab_array(std::cout, v_weight_hair, "v_weight_hair");
    write_matlab_array(std::cout, w_weight_hair, "w_weight_hair");
    
    write_matlab_array(std::cout, u_weight_total, "u_weight_total");
    write_matlab_array(std::cout, v_weight_total, "v_weight_total");
    write_matlab_array(std::cout, w_weight_total, "w_weight_total");
    
    write_matlab_array(std::cout, liquid_phi, "liquid_phi");
    
    
    exit(0);
  }
  
  threadutils::thread_pool::ParallelFor(0, (int) dof_ijk.size(), [&] (int dof_idx) {
    const int i = dof_ijk[dof_idx](0);
    const int j = dof_ijk[dof_idx](1);
    const int k = dof_ijk[dof_idx](2);
    
    pressure[k * ni * nj + j * ni + i] = x[dof_idx];
  });
  
  //Apply the velocity update
  u_valid.assign(0);
  u_pressure_grad.assign(0.0);
  int compute_num = u.ni*u.nj*u.nk;
  slice = u.ni*u.nj;
  threadutils::thread_pool::ParallelFor(0, compute_num, [&](int thread_idx)
  {
    int k = thread_idx/slice;
    int j = (thread_idx%slice)/u.ni;
    int i = thread_idx%u.ni;
    if(k<u.nk && j<u.nj && i<u.ni-1 && i>0)
    {
      int index = i + j*ni + k*ni*nj;
      if(u_weights(i,j,k) > 0) {
        if(liquid_phi(i,j,k) < 0 || liquid_phi(i-1,j,k) < 0) {
          float theta = 1;
          if(liquid_phi(i,j,k) >= 0 || liquid_phi(i-1,j,k) >= 0)
            theta = fraction_inside(liquid_phi(i-1,j,k), liquid_phi(i,j,k));
          if(theta < 0.01) theta = 0.01;
          scalar pressure_grad = (pressure[index] - pressure[index-1]) / dx / theta;
          u(i,j,k) = u_particle(i, j, k) - dt * pressure_grad / rho;
          u_pressure_grad(i, j, k) = pressure_grad;
          u_valid(i,j,k) = 1;
        }
      } else {
        u(i, j, k) = 0.0;
      }
    }
  });
  
  v_valid.assign(0);
  v_pressure_grad.assign(0.0);
  compute_num = v.ni*v.nj*v.nk;
  slice = v.ni*v.nj;
  threadutils::thread_pool::ParallelFor(0, compute_num, [&](int thread_idx)
  {
    int k = thread_idx/slice;
    int j = (thread_idx%slice)/v.ni;
    int i = thread_idx%v.ni;
    if(k<v.nk && j>0 && j<v.nj-1 && i<v.ni)
    {
      int index = i + j*ni + k*ni*nj;
      if(v_weights(i,j,k) > 0) {
        if (liquid_phi(i,j,k) < 0 || liquid_phi(i,j-1,k) < 0) {
          float theta = 1;
          if(liquid_phi(i,j,k) >= 0 || liquid_phi(i,j-1,k) >= 0)
            theta = fraction_inside(liquid_phi(i,j-1,k), liquid_phi(i,j,k));
          if(theta < 0.01) theta = 0.01;
          scalar pressure_grad = (pressure[index] - pressure[index-ni]) / dx / theta;
          v(i,j,k) = v_particle(i, j, k) - dt * pressure_grad / rho;
          v_pressure_grad(i, j, k) = pressure_grad;
          v_valid(i,j,k) = 1;
        }
      } else {
        v(i, j, k) = 0.0;
      }
    }
  });

  w_valid.assign(0);
  w_pressure_grad.assign(0.0);
  compute_num = w.ni*w.nj*w.nk;
  slice = w.ni*w.nj;
  threadutils::thread_pool::ParallelFor(0, compute_num, [&](int thread_idx)
  {
    int k = thread_idx/slice;
    int j = (thread_idx%slice)/w.ni;
    int i = thread_idx%w.ni;
    if(k>0 && k<w.nk-1 && j<w.nj && i<w.ni)
    {
      int index = i + j*ni + k*ni*nj;
      if(w_weights(i,j,k) > 0) {
        if(liquid_phi(i,j,k) < 0 || liquid_phi(i,j,k-1) < 0) {
          float theta = 1;
          if(liquid_phi(i,j,k) >= 0 || liquid_phi(i,j,k-1) >= 0)
            theta = fraction_inside(liquid_phi(i,j,k-1), liquid_phi(i,j,k));
          if(theta < 0.01) theta = 0.01;
          scalar pressure_grad = (pressure[index] - pressure[index-ni*nj]) / dx / theta;
          w(i,j,k) = w_particle(i, j, k) - dt * pressure_grad / rho;
          w_pressure_grad(i, j, k) = pressure_grad;
          w_valid(i,j,k) = 1;
        }
      } else {
        w(i, j, k) = 0.0;
      }
    }
  });
}

scalar FluidSim3D::compute_phi_vel(const Vector3s& pos, Vector3s& vel) const
{
  scalar min_phi = std::numeric_limits<scalar>::max();
  
  vel.setZero();
  for(auto& b : boundaries)
  {
    if(!b->is_root()) continue;
    Vector3s temp_vel;
    scalar phi = b->compute_phi_vel(pos, temp_vel);
    if(phi < min_phi) {
      min_phi = phi;
      vel = temp_vel;
    }
  }

  return min_phi;
}


void FluidSim3D::init_hair_particles()
{
  const std::vector< HairFlow<3>* >& flows = m_parent->getFilmFlows();
  const VectorXs& x = m_parent->getX();
  const VectorXs& v = m_parent->getV();
  const std::vector<std::pair<int, int> > edges = m_parent->getEdges();
  
  const scalar release_radius = dx * sqrt(3.0) * 0.25;
  // const scalar default_distance = release_radius / sqrt(3.0);
  
  for(const HairFlow<3>* flow : flows)
  {
    scalar accu_length = 0.0;
    
    const std::vector< int >& edge_indices = flow->getEdgeIndices();
    const std::vector< int >& indices = flow->getParticleIndices();
    const VectorXs& eta = flow->getEta();
    const VectorXs& radii_v = flow->getRadiiV();
    int np = eta.size();
    
    scalar total_length = flow->getTotalLength();
    
    while(accu_length < total_length) {
      int ipos;
      scalar alpha;
      flow->geodesic_to_local(accu_length, ipos, alpha);
      if(ipos >= np - 1) break;
      
      const Vector3s& x0 = x.segment<3>( m_parent->getDof(indices[ipos]) );
      const Vector3s& x1 = x.segment<3>( m_parent->getDof(indices[ipos + 1]) );
      const Vector3s& v0 = v.segment<3>( m_parent->getDof(indices[ipos]) );
      const Vector3s& v1 = v.segment<3>( m_parent->getDof(indices[ipos + 1]) );
      
      Vector3s pos = mathutils::lerp(x0, x1, alpha);
      Vector3s vel = mathutils::lerp(v0, v1, alpha);
      
      scalar H0 = eta(ipos) + radii_v(ipos);
      scalar H1 = eta(ipos + 1) + radii_v(ipos + 1);
      
      scalar radius = sqrt(mathutils::lerp(H0 * H0, H1 * H1, alpha));
      
      particles.push_back(Particle<3>(pos, vel, radius, PT_HAIR, edge_indices[ipos], alpha));
      
      accu_length += release_radius;
    }
  }
}

void FluidSim3D::init_random_particles(const scalar& rl, const scalar& rr, const scalar& rb, const scalar& rt, const scalar& rf, const scalar& rk)
{
  particles.clear();
  
  int num_particle = ni * nj * nk;
  
  particles.reserve(num_particle * default_particle_in_cell());
  
  scalar srl = rl * ni * dx;
  scalar srr = rr * ni * dx;
  
  scalar srb = rb * nj * dx;
  scalar srt = rt * nj * dx;
  
  scalar srf = rf * nk * dx;
  scalar srk = rk * nk * dx;
  
  for(int k = 0; k < nk; ++k)
  {
    for(int i = 0; i < ni; ++i)
    {
      for(int j = 0; j < nj; ++j) {
        for(int r = 0; r < default_particle_in_cell(); ++r) {
          scalar x = ((scalar) i + 0.5 + (((scalar)rand() / (scalar)RAND_MAX) * 2.0 - 1.0) ) * dx;
          scalar y = ((scalar) j + 0.5 + (((scalar)rand() / (scalar)RAND_MAX) * 2.0 - 1.0) ) * dx;
          scalar z = ((scalar) k + 0.5 + (((scalar)rand() / (scalar)RAND_MAX) * 2.0 - 1.0) ) * dx;
          Vector3s pt = Vector3s(x,y,z) + origin;
          
          Vector3s vel;
          if(compute_phi_vel(pt, vel) > 0 && x >= srl && x <= srr && y >= srb && y <= srt && z >= srf && z <= srk )
            add_particle(Particle<3>(pt, Vector3s::Zero(), dx * default_radius_multiplier(), PT_LIQUID));
        }
      }
    }
  }
}

void FluidSim3D::sort_particles()
{
  m_sorter->sort(particles.size(), sorter_callback);
}

scalar FluidSim3D::dropvol(const scalar& radii) const
{
  return 4.0 / 3.0 * M_PI * radii * radii * radii;
}

scalar FluidSim3D::dropradius(const scalar& vol) const
{
  return pow(vol / M_PI * 0.75, 1.0 / 3.0);
}

void FluidSim3D::map_p2g(bool with_hair_particles)
{
  //u-component of velocity
  threadutils::thread_pool::ParallelFor(0, nk, [&] (int k){
    for(int j = 0; j < nj; ++j) for(int i = 0; i < ni+1; ++i) {
      Vector3s pos = Vector3s(i*dx, (j+0.5)*dx, (k+0.5)*dx) + origin;
      scalar sumw = 0.0;
      scalar sumw_pure = 0.0;
      scalar sumu = 0.0;
      
      m_sorter->getNeigboringParticles_cell(i, j, k, -1, 0, -1, 1, -1, 1, [&] (int pidx) {
        const Particle<3>& p = particles[pidx];
        if(!with_hair_particles && p.type == PT_HAIR) return;
        
        Vector3s diff = p.x - pos;
        
        scalar w = dropvol(p.radii) * linear_kernel(diff, dx);
        sumu += w * (p.v(0) - p.c.col(0).dot(diff));
        sumw += w;
        if(p.type == PT_LIQUID) sumw_pure += w;
      });
      
      u_particle(i, j, k) = sumw ? sumu / sumw : 0.0;
      u_weight_particle(i, j, k) = sumw_pure;
    }
  });
  
  //v-component of velocity
  threadutils::thread_pool::ParallelFor(0, nk, [&] (int k){
    for(int j = 0; j < nj+1; ++j) for(int i = 0; i < ni; ++i) {
      Vector3s pos = Vector3s((i+0.5)*dx, j*dx, (k+0.5)*dx) + origin;

      scalar sumw = 0.0;
      scalar sumw_pure = 0.0;
      scalar sumu = 0.0;
      
      m_sorter->getNeigboringParticles_cell(i, j, k, -1, 1, -1, 0, -1, 1, [&] (int pidx) {
        const Particle<3>& p = particles[pidx];
        if(!with_hair_particles && p.type == PT_HAIR) return;
        Vector3s diff = p.x - pos;
        
        scalar w = dropvol(p.radii) * linear_kernel(diff, dx);
        sumu += w * (p.v(1) - p.c.col(1).dot(diff));
        sumw += w;
        if(p.type == PT_LIQUID) sumw_pure += w;
      });
      
      v_particle(i, j, k) = sumw ? sumu / sumw : 0.0;
      v_weight_particle(i, j, k) = sumw_pure;
    }
  });
  
  //w-component of velocity
  threadutils::thread_pool::ParallelFor(0, nk+1, [&] (int k){
    for(int j = 0; j < nj; ++j) for(int i = 0; i < ni; ++i) {
      Vector3s pos = Vector3s((i+0.5)*dx, (j+0.5)*dx, k*dx) + origin;
      
      scalar sumw = 0.0;
      scalar sumu = 0.0;
      scalar sumw_pure = 0.0;
      m_sorter->getNeigboringParticles_cell(i, j, k, -1, 1, -1, 1, -1, 0, [&] (int pidx) {
        const Particle<3>& p = particles[pidx];
        if(!with_hair_particles && p.type == PT_HAIR) return;
        Vector3s diff = p.x - pos;
        
        scalar w = dropvol(p.radii) * linear_kernel(diff, dx);
        sumu += w * (p.v(2) - p.c.col(2).dot(diff));
        sumw += w;
        if(p.type == PT_LIQUID) sumw_pure += w;
      });

      w_particle(i, j, k) = sumw ? sumu / sumw : 0.0;
      w_weight_particle(i, j, k) = sumw_pure;
    }
  });
}

void FluidSim3D::map_g2p_apic()
{
  int np = particles.size();
  threadutils::thread_pool::ParallelFor(0, np, [&] (int k){
    Particle<3>& p = particles[k];
    p.v = get_velocity(p.x);
    p.c = get_affine_matrix(p.x);
    p.pressure = get_pressure(p.x);
  });
}

void FluidSim3D::prepare_update_from_hair()
{
  int nflows = m_parent->getNumFlows();
  if(u_edge_vel_drag.size() != nflows) u_edge_vel_drag.resize(nflows);
  if(v_edge_vel_drag.size() != nflows) v_edge_vel_drag.resize(nflows);
  if(w_edge_vel_drag.size() != nflows) w_edge_vel_drag.resize(nflows);
  
  if(u_num_edge_voxel_intersections.size() != nflows) u_num_edge_voxel_intersections.resize(nflows);
  if(v_num_edge_voxel_intersections.size() != nflows) v_num_edge_voxel_intersections.resize(nflows);
  if(w_num_edge_voxel_intersections.size() != nflows) w_num_edge_voxel_intersections.resize(nflows);
}

void FluidSim3D::done_update_from_hair()
{
  int nflows = m_parent->getNumFlows();
  if(!nflows) return;
  
  threadutils::thread_pool::ParallelFor(0, nflows, [&] (int i) {
    u_num_edge_voxel_intersections[i] = u_edge_vel_drag[i].size();
    v_num_edge_voxel_intersections[i] = v_edge_vel_drag[i].size();
    w_num_edge_voxel_intersections[i] = w_edge_vel_drag[i].size();
  });
  
  std::partial_sum(u_num_edge_voxel_intersections.begin(), u_num_edge_voxel_intersections.end(), u_num_edge_voxel_intersections.begin());
  std::partial_sum(v_num_edge_voxel_intersections.begin(), v_num_edge_voxel_intersections.end(), v_num_edge_voxel_intersections.begin());
  std::partial_sum(w_num_edge_voxel_intersections.begin(), w_num_edge_voxel_intersections.end(), w_num_edge_voxel_intersections.begin());
  
  int usize = u_num_edge_voxel_intersections[nflows - 1];
  int vsize = v_num_edge_voxel_intersections[nflows - 1];
  int wsize = w_num_edge_voxel_intersections[nflows - 1];
  
  if(usize == 0 && vsize == 0 && wsize == 0) return;
  
  u_vel_drag.resize(usize);
  v_vel_drag.resize(vsize);
  w_vel_drag.resize(wsize);
  
  threadutils::thread_pool::ParallelFor(0, nflows, [&] (int i) {
    int u_base_idx = (i == 0) ? 0 : u_num_edge_voxel_intersections[i - 1];
    memcpy(&u_vel_drag[u_base_idx], &u_edge_vel_drag[i][0], u_edge_vel_drag[i].size() * sizeof(EdgeVelDragIntersection<3>));
    
    int v_base_idx = (i == 0) ? 0 : v_num_edge_voxel_intersections[i - 1];
    memcpy(&v_vel_drag[v_base_idx], &v_edge_vel_drag[i][0], v_edge_vel_drag[i].size() * sizeof(EdgeVelDragIntersection<3>));
    
    int w_base_idx = (i == 0) ? 0 : w_num_edge_voxel_intersections[i - 1];
    memcpy(&w_vel_drag[w_base_idx], &w_edge_vel_drag[i][0], w_edge_vel_drag[i].size() * sizeof(EdgeVelDragIntersection<3>));
  });

  m_sorter->sort(u_vel_drag.size(), [&] (int pidx, int& i, int& j, int& k) {
    i = std::max(0, std::min(m_sorter->ni - 1, u_vel_drag[pidx].coord(0)));
    j = std::max(0, std::min(m_sorter->nj - 1, u_vel_drag[pidx].coord(1)));
    k = std::max(0, std::min(m_sorter->nk - 1, u_vel_drag[pidx].coord(2)));
  });
  
  const scalar rho_L = m_parent->getLiquidDensity();
  
  threadutils::thread_pool::ParallelFor(0, u.nk, [&] (int k){
    for(int j = 0; j < u.nj; ++j) for(int i = 0; i < u.ni; ++i) {
      scalar sum_vel = 0.0;
      scalar sum_drag = 0.0;
      scalar sum_weight = 0.0;
      scalar sum_linear_weight = 0.0;
      scalar sum_vol = 0.0;
      
      m_sorter->getCellAt(i, j, k, [&] (int idx) {
        const EdgeVelDragIntersection<3>& inter = u_vel_drag[idx];
        sum_vel += inter.vel_weighted;
        sum_drag += inter.drag_weighted;
        sum_weight += inter.weight;
        sum_linear_weight += inter.linear_weight;
        sum_vol += inter.vol_weighted;
      });
      
      if(sum_weight > 0) {
        u_hair(i, j, k) = sum_vel / sum_weight;
        scalar multiplier2 = mathutils::clamp(m_parent->getDragRadiusMultiplier() * m_parent->getDragRadiusMultiplier(), 0.0, sum_linear_weight * dx * dx * dx / sum_vol);
        u_drag(i, j, k) = sum_drag / (sum_linear_weight * dx * dx * dx * rho_L) * multiplier2;
      } else {
        u_hair(i, j, k) = 0.0;
        u_drag(i, j, k) = 0.0;
      }
      
      u_weight_hair(i, j, k) = sum_vol;
    }
  });
  
  m_sorter->sort(v_vel_drag.size(), [&] (int pidx, int& i, int& j, int& k) {
    i = std::max(0, std::min(m_sorter->ni - 1, v_vel_drag[pidx].coord(0)));
    j = std::max(0, std::min(m_sorter->nj - 1, v_vel_drag[pidx].coord(1)));
    k = std::max(0, std::min(m_sorter->nk - 1, v_vel_drag[pidx].coord(2)));
  });
  
  threadutils::thread_pool::ParallelFor(0, v.nk, [&] (int k){
    for(int j = 0; j < v.nj; ++j) for(int i = 0; i < v.ni; ++i) {
      scalar sum_vel = 0.0;
      scalar sum_drag = 0.0;
      scalar sum_weight = 0.0;
      scalar sum_linear_weight = 0.0;
      scalar sum_vol = 0.0;
      
      m_sorter->getCellAt(i, j, k, [&] (int idx) {
        const EdgeVelDragIntersection<3>& inter = v_vel_drag[idx];
        sum_vel += inter.vel_weighted;
        sum_drag += inter.drag_weighted;
        sum_weight += inter.weight;
        sum_linear_weight += inter.linear_weight;
        sum_vol += inter.vol_weighted;
      });
      
      if(sum_weight > 0) {
        v_hair(i, j, k) = sum_vel / sum_weight;
        scalar multiplier2 = mathutils::clamp(m_parent->getDragRadiusMultiplier() * m_parent->getDragRadiusMultiplier(), 0.0, sum_linear_weight * dx * dx * dx / sum_vol);
        v_drag(i, j, k) = sum_drag / (sum_linear_weight * dx * dx * dx * rho_L) * multiplier2;
      } else {
        v_hair(i, j, k) = 0.0;
        v_drag(i, j, k) = 0.0;
      }
      
      v_weight_hair(i, j, k) = sum_vol;
    }
  });
  
  m_sorter->sort(w_vel_drag.size(), [&] (int pidx, int& i, int& j, int& k) {
    i = std::max(0, std::min(m_sorter->ni - 1, w_vel_drag[pidx].coord(0)));
    j = std::max(0, std::min(m_sorter->nj - 1, w_vel_drag[pidx].coord(1)));
    k = std::max(0, std::min(m_sorter->nk - 1, w_vel_drag[pidx].coord(2)));
  });
  
  threadutils::thread_pool::ParallelFor(0, w.nk, [&] (int k){
    for(int j = 0; j < w.nj; ++j) for(int i = 0; i < w.ni; ++i) {
      scalar sum_vel = 0.0;
      scalar sum_drag = 0.0;
      scalar sum_weight = 0.0;
      scalar sum_linear_weight = 0.0;
      scalar sum_vol = 0.0;
      
      m_sorter->getCellAt(i, j, k, [&] (int idx) {
        const EdgeVelDragIntersection<3>& inter = w_vel_drag[idx];
        sum_vel += inter.vel_weighted;
        sum_drag += inter.drag_weighted;
        sum_weight += inter.weight;
        sum_linear_weight += inter.linear_weight;
        sum_vol += inter.vol_weighted;
      });
      
      if(sum_weight > 0) {
        w_hair(i, j, k) = sum_vel / sum_weight;
        scalar multiplier2 = mathutils::clamp(m_parent->getDragRadiusMultiplier() * m_parent->getDragRadiusMultiplier(), 0.0, sum_linear_weight * dx * dx * dx / sum_vol);
        w_drag(i, j, k) = sum_drag / (sum_linear_weight * dx * dx * dx * rho_L) * multiplier2;
      } else {
        w_hair(i, j, k) = 0.0;
        w_drag(i, j, k) = 0.0;
      }
      
      w_weight_hair(i, j, k) = sum_vol;
    }
  });
}

scalar FluidSim3D::get_particle_weight(const Vector3s& position) const
{
  return get_velocity(position, u_weight_particle, v_weight_particle, w_weight_particle).norm();
}

scalar FluidSim3D::get_clamped_particle_weight(const Vector3s& position) const
{
  const scalar standard_radius = dx * default_radius_multiplier();
  const scalar standard_vol = dropvol(standard_radius);
  scalar w = get_particle_weight(position);
  scalar criterion = m_parent->getBulkThresholdMultiplier() * standard_vol * m_parent->getLiquidDensity();
  return mathutils::clamp(w / criterion, 0.0, 1.0);
}

scalar FluidSim3D::cfl()
{
  scalar maxvel = 0;
  const int npu = u.a.size();
  const int npv = v.a.size();
  const int npw = w.a.size();
  for(int i = 0; i < npu; ++i)
    maxvel = max(maxvel, fabs(u.a[i]));
  for(int i = 0; i < npv; ++i)
    maxvel = max(maxvel, fabs(v.a[i]));
  for(int i = 0; i < npw; ++i)
    maxvel = max(maxvel, fabs(w.a[i]));
  return dx / maxvel;
}

const std::vector< FluidSim3D::Boundary<3>* >& FluidSim3D::get_boundaries() const
{
  return boundaries;
}

const std::vector< FluidSim3D::SourceBoundary<3>* >& FluidSim3D::get_sources() const
{
  return sources;
}

const Vector3s& FluidSim3D::get_origin() const
{
  return origin;
}

int FluidSim3D::get_ni() const
{
  return ni;
}

int FluidSim3D::get_nj() const
{
  return nj;
}

int FluidSim3D::get_nk() const
{
  return nk;
}

int FluidSim3D::get_u_ni() const
{
  return u.ni;
}

int FluidSim3D::get_v_ni() const
{
  return v.ni;
}

int FluidSim3D::get_w_ni() const
{
  return w.ni;
}

int FluidSim3D::get_u_nj() const
{
  return u.nj;
}

int FluidSim3D::get_v_nj() const
{
  return v.nj;
}

int FluidSim3D::get_w_nj() const
{
  return w.nj;
}

int FluidSim3D::get_u_nk() const
{
  return u.nk;
}

int FluidSim3D::get_v_nk() const
{
  return v.nk;
}

int FluidSim3D::get_w_nk() const
{
  return w.nk;
}

std::vector< std::vector<EdgeVelDragIntersection<3> > >& FluidSim3D::get_u_edge_vel_drag()
{
  return u_edge_vel_drag;
}

std::vector< std::vector<EdgeVelDragIntersection<3> > >& FluidSim3D::get_v_edge_vel_drag()
{
  return v_edge_vel_drag;
}

std::vector< std::vector<EdgeVelDragIntersection<3> > >& FluidSim3D::get_w_edge_vel_drag()
{
  return w_edge_vel_drag;
}

const std::vector<Particle<3> >& FluidSim3D::get_particles() const
{
  return particles;
}

scalar FluidSim3D::cellsize() const
{
  return dx;
}

void FluidSim3D::save_pressure(const std::string szfn)
{
  using namespace std;
  
  ofstream ofs(szfn.c_str());
  ofs << "pressure = [" << endl;
  for(int k = 0; k < nk; ++k) {
    for(int j = 0; j < nj; ++j) {
      for(int i = 0; i < ni; ++i) {
        ofs << pressure[k * ni * nj + j * ni + i] << "\t";
      }
      ofs << "\n";
    }
  }
  
  ofs << "];" << endl;
  ofs.close();
}

Vector3s FluidSim3D::computeParticleMomentum()
{
  int np = particles.size();
  
  Vector3s sum = Vector3s::Zero();
  
  scalar rho = m_parent->getLiquidDensity();
  
  for(int i = 0; i < np; ++i)
  {
    auto& p = particles[i];
    scalar m = 4.0 / 3.0 * M_PI * p.radii * p.radii * p.radii * rho;
    sum += p.v * m;
  }
  
  return sum;
}

Vector3s FluidSim3D::computeParticleAngularMomentum()
{
  int np = particles.size();
  
  Vector3s sum = Vector3s::Zero();
  
  scalar rho = m_parent->getLiquidDensity();
  
  for(int i = 0; i < np; ++i)
  {
    auto& p = particles[i];
    scalar m = 4.0 / 3.0 * M_PI * p.radii * p.radii * p.radii * rho;
    
    Vector3s ppu = (p.x - origin) / dx - Vector3s(0.0, 0.5, 0.5);
    Vector3s ppv = (p.x - origin) / dx - Vector3s(0.5, 0.0, 0.5);
    Vector3s ppw = (p.x - origin) / dx - Vector3s(0.0, 0.5, 0.5);
    
    int iu_low = std::max(0, std::min(u.ni-2, (int) ppu(0)));
    int iu_high = iu_low + 1;
    
    int ju_low = std::max(0, std::min(u.nj-2, (int) ppu(1)));
    int ju_high = ju_low + 1;
    
    int ku_low = std::max(0, std::min(u.nk-2, (int) ppu(2)));
    int ku_high = ku_low + 1;
    
    scalar iu_frac = ppu(0) - floor(ppu(0));
    scalar ju_frac = ppu(1) - floor(ppu(1));
    scalar ku_frac = ppu(2) - floor(ppu(2));
    
    Vector3s xu000 = Vector3s(iu_low*dx, (ju_low+0.5)*dx, (ku_low+0.5)*dx) + origin;
    Vector3s xu100 = Vector3s(iu_high*dx, (ju_low+0.5)*dx, (ku_low+0.5)*dx) + origin;
    Vector3s xu010 = Vector3s(iu_low*dx, (ju_high+0.5)*dx, (ku_low+0.5)*dx) + origin;
    Vector3s xu110 = Vector3s(iu_high*dx, (ju_high+0.5)*dx, (ku_low+0.5)*dx) + origin;
    Vector3s xu001 = Vector3s(iu_low*dx, (ju_low+0.5)*dx, (ku_high+0.5)*dx) + origin;
    Vector3s xu101 = Vector3s(iu_high*dx, (ju_low+0.5)*dx, (ku_high+0.5)*dx) + origin;
    Vector3s xu011 = Vector3s(iu_low*dx, (ju_high+0.5)*dx, (ku_high+0.5)*dx) + origin;
    Vector3s xu111 = Vector3s(iu_high*dx, (ju_high+0.5)*dx, (ku_high+0.5)*dx) + origin;
    
    Vector3s pu000 = m * (1.0 - iu_frac) * (1.0 - ju_frac) * (1.0 - ku_frac) * Vector3s(p.v(0) + p.c.col(0).dot(xu000 - p.x), 0.0, 0.0);
    Vector3s pu100 = m * iu_frac * (1.0 - ju_frac) * (1.0 - ku_frac) * Vector3s(p.v(0) + p.c.col(0).dot(xu100 - p.x), 0.0, 0.0);
    Vector3s pu010 = m * (1.0 - iu_frac) * ju_frac * (1.0 - ku_frac) * Vector3s(p.v(0) + p.c.col(0).dot(xu010 - p.x), 0.0, 0.0);
    Vector3s pu110 = m * iu_frac * ju_frac * (1.0 - ku_frac) * Vector3s(p.v(0) + p.c.col(0).dot(xu110 - p.x), 0.0, 0.0);
    Vector3s pu001 = m * (1.0 - iu_frac) * (1.0 - ju_frac) * ku_frac * Vector3s(p.v(0) + p.c.col(0).dot(xu001 - p.x), 0.0, 0.0);
    Vector3s pu101 = m * iu_frac * (1.0 - ju_frac) * ku_frac * Vector3s(p.v(0) + p.c.col(0).dot(xu101 - p.x), 0.0, 0.0);
    Vector3s pu011 = m * (1.0 - iu_frac) * ju_frac * ku_frac * Vector3s(p.v(0) + p.c.col(0).dot(xu011 - p.x), 0.0, 0.0);
    Vector3s pu111 = m * iu_frac * ju_frac * ku_frac * Vector3s(p.v(0) + p.c.col(0).dot(xu111 - p.x), 0.0, 0.0);
    
    sum +=
    xu000.cross(pu000) +
    xu100.cross(pu100) +
    xu010.cross(pu010) +
    xu110.cross(pu110) +
    xu001.cross(pu001) +
    xu101.cross(pu101) +
    xu011.cross(pu011) +
    xu111.cross(pu111);
    
    int iv_low = std::max(0, std::min(v.ni-2, (int) ppv(0)));
    int iv_high = iv_low + 1;
    
    int jv_low = std::max(0, std::min(v.nj-2, (int) ppv(1)));
    int jv_high = jv_low + 1;
    
    int kv_low = std::max(0, std::min(v.nk-2, (int) ppv(2)));
    int kv_high = kv_low + 1;
    
    scalar iv_frac = ppv(0) - floor(ppv(0));
    scalar jv_frac = ppv(1) - floor(ppv(1));
    scalar kv_frac = ppv(2) - floor(ppv(2));
    
    Vector3s xv000 = Vector3s((iv_low+0.5)*dx, jv_low*dx, (kv_low+0.5)*dx) + origin;
    Vector3s xv100 = Vector3s((iv_high+0.5)*dx, jv_low*dx, (kv_low+0.5)*dx) + origin;
    Vector3s xv010 = Vector3s((iv_low+0.5)*dx, jv_high*dx, (kv_low+0.5)*dx) + origin;
    Vector3s xv110 = Vector3s((iv_high+0.5)*dx, jv_high*dx, (kv_low+0.5)*dx) + origin;
    Vector3s xv001 = Vector3s((iv_low+0.5)*dx, jv_low*dx, (kv_high+0.5)*dx) + origin;
    Vector3s xv101 = Vector3s((iv_high+0.5)*dx, jv_low*dx, (kv_high+0.5)*dx) + origin;
    Vector3s xv011 = Vector3s((iv_low+0.5)*dx, jv_high*dx, (kv_high+0.5)*dx) + origin;
    Vector3s xv111 = Vector3s((iv_high+0.5)*dx, jv_high*dx, (kv_high+0.5)*dx) + origin;
    
    Vector3s pv000 = m * (1.0 - iv_frac) * (1.0 - jv_frac) * (1.0 - kv_frac) * Vector3s(0.0, p.v(1) + p.c.col(1).dot(xv000 - p.x), 0.0);
    Vector3s pv100 = m * iv_frac * (1.0 - jv_frac) * (1.0 - kv_frac) * Vector3s(0.0, p.v(1) + p.c.col(1).dot(xv100 - p.x), 0.0);
    Vector3s pv010 = m * (1.0 - iv_frac) * jv_frac * (1.0 - kv_frac) * Vector3s(0.0, p.v(1) + p.c.col(1).dot(xv010 - p.x), 0.0);
    Vector3s pv110 = m * iv_frac * jv_frac * (1.0 - kv_frac) * Vector3s(0.0, p.v(1) + p.c.col(1).dot(xv110 - p.x), 0.0);
    Vector3s pv001 = m * (1.0 - iv_frac) * (1.0 - jv_frac) * kv_frac * Vector3s(0.0, p.v(1) + p.c.col(1).dot(xv001 - p.x), 0.0);
    Vector3s pv101 = m * iv_frac * (1.0 - jv_frac) * kv_frac * Vector3s(0.0, p.v(1) + p.c.col(1).dot(xv101 - p.x), 0.0);
    Vector3s pv011 = m * (1.0 - iv_frac) * jv_frac * kv_frac * Vector3s(0.0, p.v(1) + p.c.col(1).dot(xv011 - p.x), 0.0);
    Vector3s pv111 = m * iv_frac * jv_frac * kv_frac * Vector3s(0.0, p.v(1) + p.c.col(1).dot(xv111 - p.x), 0.0);

    sum +=
    xv000.cross(pv000) +
    xv100.cross(pv100) +
    xv010.cross(pv010) +
    xv110.cross(pv110) +
    xv001.cross(pv001) +
    xv101.cross(pv101) +
    xv011.cross(pv011) +
    xv111.cross(pv111);
    
    
    int iw_low = std::max(0, std::min(w.ni-2, (int) ppw(0)));
    int iw_high = iw_low + 1;
    
    int jw_low = std::max(0, std::min(w.nj-2, (int) ppw(1)));
    int jw_high = jw_low + 1;
    
    int kw_low = std::max(0, std::min(w.nk-2, (int) ppw(2)));
    int kw_high = kw_low + 1;
    
    scalar iw_frac = ppw(0) - floor(ppw(0));
    scalar jw_frac = ppw(1) - floor(ppw(1));
    scalar kw_frac = ppw(2) - floor(ppw(2));
    
    Vector3s xw000 = Vector3s((iw_low+0.5)*dx, (jw_low+0.5)*dx, kw_low*dx) + origin;
    Vector3s xw100 = Vector3s((iw_high+0.5)*dx, (jw_low+0.5)*dx, kw_low*dx) + origin;
    Vector3s xw010 = Vector3s((iw_low+0.5)*dx, (jw_high+0.5)*dx, kw_low*dx) + origin;
    Vector3s xw110 = Vector3s((iw_high+0.5)*dx, (jw_high+0.5)*dx, kw_low*dx) + origin;
    Vector3s xw001 = Vector3s((iw_low+0.5)*dx, (jw_low+0.5)*dx, kw_high*dx) + origin;
    Vector3s xw101 = Vector3s((iw_high+0.5)*dx, (jw_low+0.5)*dx, kw_high*dx) + origin;
    Vector3s xw011 = Vector3s((iw_low+0.5)*dx, (jw_high+0.5)*dx, kw_high*dx) + origin;
    Vector3s xw111 = Vector3s((iw_high+0.5)*dx, (jw_high+0.5)*dx, kw_high*dx) + origin;
    
    Vector3s pw000 = m * (1.0 - iw_frac) * (1.0 - jw_frac) * (1.0 - kw_frac) * Vector3s(0.0, 0.0, p.v(2) + p.c.col(2).dot(xv000 - p.x));
    Vector3s pw100 = m * iw_frac * (1.0 - jw_frac) * (1.0 - kw_frac) * Vector3s(0.0, 0.0, p.v(2) + p.c.col(2).dot(xv100 - p.x));
    Vector3s pw010 = m * (1.0 - iw_frac) * jw_frac * (1.0 - kw_frac) * Vector3s(0.0, 0.0, p.v(2) + p.c.col(2).dot(xv010 - p.x));
    Vector3s pw110 = m * iw_frac * jw_frac * (1.0 - kw_frac) * Vector3s(0.0, 0.0, p.v(2) + p.c.col(2).dot(xv110 - p.x));
    Vector3s pw001 = m * (1.0 - iw_frac) * (1.0 - jw_frac) * kw_frac * Vector3s(0.0, 0.0, p.v(2) + p.c.col(2).dot(xv001 - p.x));
    Vector3s pw101 = m * iw_frac * (1.0 - jw_frac) * kw_frac * Vector3s(0.0, 0.0, p.v(2) + p.c.col(2).dot(xv101 - p.x));
    Vector3s pw011 = m * (1.0 - iw_frac) * jw_frac * kw_frac * Vector3s(0.0, 0.0, p.v(2) + p.c.col(2).dot(xv011 - p.x));
    Vector3s pw111 = m * iw_frac * jw_frac * kw_frac * Vector3s(0.0, 0.0, p.v(2) + p.c.col(2).dot(xv111 - p.x));
    
    sum +=
    xw000.cross(pw000) +
    xw100.cross(pw100) +
    xw010.cross(pw010) +
    xw110.cross(pw110) +
    xw001.cross(pw001) +
    xw101.cross(pw101) +
    xw011.cross(pw011) +
    xw111.cross(pw111);
  }
  
  return sum;
}

Vector3s FluidSim3D::computeParticleGridAngularMomentum()
{
  Vector3s sum = Vector3s::Zero();
  
  for(int k = 0; k < u_particle.nk; ++k) for(int j = 0; j < u_particle.nj; ++j) for(int i = 0; i < u_particle.ni; ++i)
  {
    Vector3s pos = Vector3s(i*dx, (j+0.5)*dx, (k+0.5)*dx) + origin;
    Vector3s u_vec = Vector3s(u_weight_particle(i,j,k) * u_particle(i,j,k), 0.0, 0.0);
    sum += pos.cross(u_vec);
  }
  
  for(int k = 0; k < v_particle.nk; ++k) for(int j = 0; j < v_particle.nj; ++j) for(int i = 0; i < v_particle.ni; ++i)
  {
    Vector3s pos = Vector3s((i+0.5)*dx, j*dx, (k+0.5)*dx) + origin;
    Vector3s v_vec = Vector3s(0.0, v_weight_particle(i,j,k) * v_particle(i,j,k), 0.0);
    sum += pos.cross(v_vec);
  }
  
  for(int k = 0; k < w_particle.nk; ++k) for(int j = 0; j < w_particle.nj; ++j) for(int i = 0; i < w_particle.ni; ++i)
  {
    Vector3s pos = Vector3s((i+0.5)*dx, (j+0.5)*dx, k*dx) + origin;
    Vector3s w_vec = Vector3s(0.0, 0.0, w_weight_particle(i,j,k) * w_particle(i,j,k));
    sum += pos.cross(w_vec);
  }
  
  return sum;
}

Vector3s FluidSim3D::computeHairGridAngularMomentum()
{
  Vector3s sum = Vector3s::Zero();
  
  for(int k = 0; k < u_hair.nk; ++k) for(int j = 0; j < u_hair.nj; ++j) for(int i = 0; i < u_hair.ni; ++i)
  {
    Vector3s pos = Vector3s(i*dx, (j+0.5)*dx, (k+0.5)*dx) + origin;
    Vector3s u_vec = Vector3s(u_weight_hair(i,j,k) * u_hair(i,j,k), 0.0, 0.0);
    sum += pos.cross(u_vec);
  }
  
  for(int k = 0; k < v_hair.nk; ++k) for(int j = 0; j < v_hair.nj; ++j) for(int i = 0; i < v_hair.ni; ++i)
  {
    Vector3s pos = Vector3s((i+0.5)*dx, j*dx, (k+0.5)*dx) + origin;
    Vector3s v_vec = Vector3s(0.0, v_weight_hair(i,j,k) * v_hair(i,j,k), 0.0);
    sum += pos.cross(v_vec);
  }
  
  for(int k = 0; k < w_hair.nk; ++k) for(int j = 0; j < w_hair.nj; ++j) for(int i = 0; i < w_hair.ni; ++i)
  {
    Vector3s pos = Vector3s((i+0.5)*dx, (j+0.5)*dx, k*dx) + origin;
    Vector3s w_vec = Vector3s(0.0, 0.0, w_weight_hair(i,j,k) * w_hair(i,j,k));
    sum += pos.cross(w_vec);
  }
  
  return sum;
}

Vector3s FluidSim3D::computeCombinedGridAngularMomentum()
{
  Vector3s sum = Vector3s::Zero();
  
  for(int k = 0; k < u.nk; ++k) for(int j = 0; j < u.nj; ++j) for(int i = 0; i < u.ni; ++i)
  {
    Vector3s pos = Vector3s(i*dx, (j+0.5)*dx, (k+0.5)*dx) + origin;
    Vector3s u_vec = Vector3s(u_weight_total(i,j,k) * u(i,j,k), 0.0, 0.0);
    sum += pos.cross(u_vec);
  }
  
  for(int k = 0; k < v.nk; ++k) for(int j = 0; j < v.nj; ++j) for(int i = 0; i < v.ni; ++i)
  {
    Vector3s pos = Vector3s((i+0.5)*dx, j*dx, (k+0.5)*dx) + origin;
    Vector3s v_vec = Vector3s(0.0, v_weight_total(i,j,k) * v(i,j,k), 0.0);
    sum += pos.cross(v_vec);
  }
  
  for(int k = 0; k < w.nk; ++k) for(int j = 0; j < w.nj; ++j) for(int i = 0; i < w.ni; ++i)
  {
    Vector3s pos = Vector3s((i+0.5)*dx, (j+0.5)*dx, k*dx) + origin;
    Vector3s w_vec = Vector3s(0.0, 0.0, w_weight_total(i,j,k) * w(i,j,k));
    sum += pos.cross(w_vec);
  }
  
  return sum;
}

Vector3s FluidSim3D::computeParticleWeightedCombinedGridAngularMomentum()
{
  Vector3s sum = Vector3s::Zero();
  
  for(int k = 0; k < u.nk; ++k) for(int j = 0; j < u.nj; ++j) for(int i = 0; i < u.ni; ++i)
  {
    Vector3s pos = Vector3s(i*dx, (j+0.5)*dx, (k+0.5)*dx) + origin;
    Vector3s u_vec = Vector3s(u_weight_particle(i,j,k) * u(i,j,k), 0.0, 0.0);
    sum += pos.cross(u_vec);
  }
  
  for(int k = 0; k < v.nk; ++k) for(int j = 0; j < v.nj; ++j) for(int i = 0; i < v.ni; ++i)
  {
    Vector3s pos = Vector3s((i+0.5)*dx, j*dx, (k+0.5)*dx) + origin;
    Vector3s v_vec = Vector3s(0.0, v_weight_particle(i,j,k) * v(i,j,k), 0.0);
    sum += pos.cross(v_vec);
  }
  
  for(int k = 0; k < w.nk; ++k) for(int j = 0; j < w.nj; ++j) for(int i = 0; i < w.ni; ++i)
  {
    Vector3s pos = Vector3s((i+0.5)*dx, (j+0.5)*dx, k*dx) + origin;
    Vector3s w_vec = Vector3s(0.0, 0.0, w_weight_particle(i,j,k) * w(i,j,k));
    sum += pos.cross(w_vec);
  }
  
  return sum;
}

Vector3s FluidSim3D::computeHairWeightedCombinedGridAngularMomentum()
{
  Vector3s sum = Vector3s::Zero();
  
  for(int k = 0; k < u.nk; ++k) for(int j = 0; j < u.nj; ++j) for(int i = 0; i < u.ni; ++i)
  {
    Vector3s pos = Vector3s(i*dx, (j+0.5)*dx, (k+0.5)*dx) + origin;
    Vector3s u_vec = Vector3s(u_weight_hair(i,j,k) * u(i,j,k), 0.0, 0.0);
    sum += pos.cross(u_vec);
  }
  
  for(int k = 0; k < v.nk; ++k) for(int j = 0; j < v.nj; ++j) for(int i = 0; i < v.ni; ++i)
  {
    Vector3s pos = Vector3s((i+0.5)*dx, j*dx, (k+0.5)*dx) + origin;
    Vector3s v_vec = Vector3s(0.0, v_weight_hair(i,j,k) * v(i,j,k), 0.0);
    sum += pos.cross(v_vec);
  }
  
  for(int k = 0; k < w.nk; ++k) for(int j = 0; j < w.nj; ++j) for(int i = 0; i < w.ni; ++i)
  {
    Vector3s pos = Vector3s((i+0.5)*dx, (j+0.5)*dx, k*dx) + origin;
    Vector3s w_vec = Vector3s(0.0, 0.0, w_weight_hair(i,j,k) * w(i,j,k));
    sum += pos.cross(w_vec);
  }
  
  return sum;
}

Vector3s FluidSim3D::computeReweightedHairGridAngularMomentum()
{
  Vector3s sum = Vector3s::Zero();
  
  for(int k = 0; k < u_hair.nk; ++k) for(int j = 0; j < u_hair.nj; ++j) for(int i = 1; i < u_hair.ni-1; ++i)
  {
    scalar theta = 1.0;
    if(liquid_phi(i,j,k) >= 0 || liquid_phi(i-1,j,k) >= 0)
      theta = mathutils::fraction_inside(liquid_phi(i-1,j,k), liquid_phi(i,j,k));
    
    Vector3s pos = Vector3s(i*dx, (j+0.5)*dx, (k+0.5)*dx) + origin;
    Vector3s u_vec = Vector3s(u_weight_hair(i,j,k) * u_hair(i,j,k) * (1.0 - theta), 0.0, 0.0);
    sum += pos.cross(u_vec);
  }
  
  for(int k = 0; k < v_hair.nk; ++k) for(int j = 1; j < v_hair.nj-1; ++j) for(int i = 0; i < v_hair.ni; ++i)
  {
    scalar theta = 1.0;
    if(liquid_phi(i,j,k) >= 0 || liquid_phi(i,j-1,k) >= 0)
      theta = mathutils::fraction_inside(liquid_phi(i,j-1,k), liquid_phi(i,j,k));
    
    Vector3s pos = Vector3s((i+0.5)*dx, j*dx, (k+0.5)*dx) + origin;
    Vector3s v_vec = Vector3s(0.0, v_weight_hair(i,j,k) * v_hair(i,j,k) * (1.0 - theta), 0.0);
    sum += pos.cross(v_vec);
  }
  
  for(int k = 1; k < w_hair.nk-1; ++k) for(int j = 0; j < w_hair.nj; ++j) for(int i = 0; i < w_hair.ni; ++i)
  {
    scalar theta = 1.0;
    if(liquid_phi(i,j,k) >= 0 || liquid_phi(i,j,k-1) >= 0)
      theta = mathutils::fraction_inside(liquid_phi(i,j,k-1), liquid_phi(i,j,k));
    
    Vector3s pos = Vector3s((i+0.5)*dx, (j+0.5)*dx, k*dx) + origin;
    Vector3s w_vec = Vector3s(0.0, 0.0, w_weight_hair(i,j,k) * w_hair(i,j,k) * (1.0 - theta));
    sum += pos.cross(w_vec);
  }
  
  return sum;
}

Vector3s FluidSim3D::computeReweightedParticleGridAngularMomentum()
{
  Vector3s sum = Vector3s::Zero();
  
  for(int k = 0; k < u_particle.nk; ++k) for(int j = 0; j < u_particle.nj; ++j) for(int i = 1; i < u_particle.ni-1; ++i)
  {
    scalar theta = 1.0;
    if(liquid_phi(i,j,k) >= 0 || liquid_phi(i-1,j,k) >= 0)
      theta = mathutils::fraction_inside(liquid_phi(i-1,j,k), liquid_phi(i,j,k));
    
    Vector3s pos = Vector3s(i*dx, (j+0.5)*dx, (k+0.5)*dx) + origin;
    Vector3s u_vec = Vector3s(u_weight_particle(i,j,k) * u_particle(i,j,k) * theta, 0.0, 0.0);
    sum += pos.cross(u_vec);
  }
  
  for(int k = 0; k < v_particle.nk; ++k) for(int j = 1; j < v_particle.nj-1; ++j) for(int i = 0; i < v_particle.ni; ++i)
  {
    scalar theta = 1.0;
    if(liquid_phi(i,j,k) >= 0 || liquid_phi(i,j-1,k) >= 0)
      theta = mathutils::fraction_inside(liquid_phi(i,j-1,k), liquid_phi(i,j,k));
    
    Vector3s pos = Vector3s((i+0.5)*dx, j*dx, (k+0.5)*dx) + origin;
    Vector3s v_vec = Vector3s(0.0, v_weight_particle(i,j,k) * v_particle(i,j,k) * theta, 0.0);
    sum += pos.cross(v_vec);
  }
  
  for(int k = 1; k < w_particle.nk-1; ++k) for(int j = 0; j < w_particle.nj; ++j) for(int i = 0; i < w_particle.ni; ++i)
  {
    scalar theta = 1.0;
    if(liquid_phi(i,j,k) >= 0 || liquid_phi(i,j,k-1) >= 0)
      theta = mathutils::fraction_inside(liquid_phi(i,j,k-1), liquid_phi(i,j,k));
    
    Vector3s pos = Vector3s((i+0.5)*dx, (j+0.5)*dx, k*dx) + origin;
    Vector3s w_vec = Vector3s(0.0, 0.0, w_weight_particle(i,j,k) * w_particle(i,j,k) * theta);
    sum += pos.cross(w_vec);
  }
  
  return sum;
}

scalar FluidSim3D::computeTotalLiquidVol() const
{
  scalar sum = 0.0;
  for(auto& p : particles)
  {
    if(p.type != PT_LIQUID) continue;
    sum += 4.0 / 3.0 * M_PI * p.radii * p.radii * p.radii;
  }
  return sum;
}

Vector3s FluidSim3D::getMinBBX() const
{
  return Vector3s(origin(0), origin(1), origin(2));
}

Vector3s FluidSim3D::getMaxBBX() const
{
  return Vector3s(origin(0) + ni * dx, origin(1) + nj * dx, origin(2) + nk * dx);
}

Vector3s FluidSim3D::computeHairGridMomentum()
{
  Vector3s sum = Vector3s::Zero();
  
  for(int k = 0; k < u_hair.nk; ++k) for(int j = 0; j < u_hair.nj; ++j) for(int i = 0; i < u_hair.ni; ++i)
  {
    sum(0) += u_weight_hair(i, j, k) * u_hair(i, j, k);
  }
  
  for(int k = 0; k < v_hair.nk; ++k) for(int j = 0; j < v_hair.nj; ++j) for(int i = 0; i < v_hair.ni; ++i)
  {
    sum(1) += v_weight_hair(i, j, k) * v_hair(i, j, k);
  }
  
  for(int k = 0; k < w_hair.nk; ++k) for(int j = 0; j < w_hair.nj; ++j) for(int i = 0; i < w_hair.ni; ++i)
  {
    sum(2) += w_weight_hair(i, j, k) * w_hair(i, j, k);
  }
  
  return sum;
}

Vector3s FluidSim3D::computeParticleGridMomentum()
{
  Vector3s sum = Vector3s::Zero();
  
  for(int k = 0; k < u_hair.nk; ++k) for(int j = 0; j < u_hair.nj; ++j) for(int i = 0; i < u_hair.ni; ++i)
  {
    sum(0) += u_weight_particle(i, j, k) * u_particle(i, j, k);
  }
  
  for(int k = 0; k < v_hair.nk; ++k) for(int j = 0; j < v_hair.nj; ++j) for(int i = 0; i < v_hair.ni; ++i)
  {
    sum(1) += v_weight_particle(i, j, k) * v_particle(i, j, k);
  }
  
  for(int k = 0; k < w_hair.nk; ++k) for(int j = 0; j < w_hair.nj; ++j) for(int i = 0; i < w_hair.ni; ++i)
  {
    sum(2) += w_weight_particle(i, j, k) * w_particle(i, j, k);
  }
  return sum;
}

Vector3s FluidSim3D::computeReweightedHairGridMomentum()
{
  Vector3s sum = Vector3s::Zero();
  
  for(int k = 0; k < u_hair.nk; ++k) for(int j = 0; j < u_hair.nj; ++j) for(int i = 1; i < u_hair.ni - 1; ++i)
  {
    scalar theta = 1.0;
    if(liquid_phi(i,j,k) >= 0 || liquid_phi(i-1,j,k) >= 0)
      theta = mathutils::fraction_inside(liquid_phi(i-1,j,k), liquid_phi(i,j,k));
    
    sum(0) += u_weight_hair(i, j, k) * u_hair(i, j, k) * (1. - theta);
  }
  
  for(int k = 0; k < v_hair.nk; ++k) for(int j = 1; j < v_hair.nj - 1; ++j) for(int i = 0; i < v_hair.ni; ++i)
  {
    scalar theta = 1.0;
    if(liquid_phi(i,j,k) >= 0 || liquid_phi(i,j-1,k) >= 0)
      theta = mathutils::fraction_inside(liquid_phi(i,j-1,k), liquid_phi(i,j,k));
    
    sum(1) += v_weight_hair(i, j, k) * v_hair(i, j, k) * (1. - theta);
  }
  
  for(int k = 1; k < w_hair.nk - 1; ++k) for(int j = 0; j < w_hair.nj; ++j) for(int i = 0; i < w_hair.ni; ++i)
  {
    scalar theta = 1.0;
    if(liquid_phi(i,j,k) >= 0 || liquid_phi(i,j,k-1) >= 0)
      theta = mathutils::fraction_inside(liquid_phi(i,j,k-1), liquid_phi(i,j,k));
    
    sum(2) += w_weight_hair(i, j, k) * w_hair(i, j, k) * (1. - theta);
  }
  
  return sum;
}

Vector3s FluidSim3D::computeReweightedParticleGridMomentum()
{
  Vector3s sum = Vector3s::Zero();
  
  for(int k = 0; k < u_particle.nk; ++k) for(int j = 0; j < u_particle.nj; ++j) for(int i = 1; i < u_particle.ni-1; ++i)
  {
    scalar theta = 1.0;
    if(liquid_phi(i,j,k) >= 0 || liquid_phi(i-1,j,k) >= 0)
      theta = mathutils::fraction_inside(liquid_phi(i-1,j,k), liquid_phi(i,j,k));
    
    sum(0) += u_weight_particle(i, j, k) * u_particle(i, j, k) * theta;
  }
  
  for(int k = 0; k < u_particle.nk; ++k) for(int j = 1; j < v_particle.nj-1; ++j) for(int i = 0; i < v_particle.ni; ++i)
  {
    scalar theta = 1.0;
    if(liquid_phi(i,j,k) >= 0 || liquid_phi(i,j-1,k) >= 0)
      theta = mathutils::fraction_inside(liquid_phi(i,j-1,k), liquid_phi(i,j,k));
    
    sum(1) += v_weight_particle(i, j, k) * v_particle(i, j, k) * theta;
  }
  
  for(int k = 1; k < w_particle.nk-1; ++k) for(int j = 0; j < w_particle.nj; ++j) for(int i = 0; i < w_particle.ni; ++i)
  {
    scalar theta = 1.0;
    if(liquid_phi(i,j,k) >= 0 || liquid_phi(i,j,k-1) >= 0)
      theta = mathutils::fraction_inside(liquid_phi(i,j,k-1), liquid_phi(i,j,k));
    
    sum(2) += w_weight_particle(i, j, k) * w_particle(i, j, k) * theta;
  }
  
  return sum;
}

Vector3s FluidSim3D::computeParticleWeightedCombinedGridMomentum()
{
  Vector3s sum = Vector3s::Zero();
  
  for(int k = 0; k < u.nk; ++k) for(int j = 0; j < u.nj; ++j) for(int i = 1; i < u.ni-1; ++i)
  {
    sum(0) += u_weight_particle(i, j, k) * u(i, j, k);
  }
  
  for(int k = 0; k < v.nk; ++k) for(int j = 1; j < v.nj-1; ++j) for(int i = 0; i < v.ni; ++i)
  {
    sum(1) += v_weight_particle(i, j, k) * v(i, j, k);
  }
  
  for(int k = 1; k < w.nk-1; ++k) for(int j = 0; j < w.nj; ++j) for(int i = 0; i < w.ni; ++i)
  {
    sum(2) += w_weight_particle(i, j, k) * v(i, j, k);
  }
  
  return sum;
}

Vector3s FluidSim3D::computeHairWeightedCombinedGridMomentum()
{
  Vector3s sum = Vector3s::Zero();
  
  for(int k = 0; k < u.nk; ++k) for(int j = 0; j < u.nj; ++j) for(int i = 1; i < u.ni-1; ++i)
  {
    sum(0) += u_weight_hair(i, j, k) * u(i, j, k);
  }
  
  for(int k = 0; k < v.nk; ++k) for(int j = 1; j < v.nj-1; ++j) for(int i = 0; i < v.ni; ++i)
  {
    sum(1) += v_weight_hair(i, j, k) * v(i, j, k);
  }
  
  for(int k = 1; k < w.nk-1; ++k) for(int j = 0; j < w.nj; ++j) for(int i = 0; i < w.ni; ++i)
  {
    sum(2) += w_weight_hair(i, j, k) * w(i, j, k);
  }
  
  return sum;
}

Vector3s FluidSim3D::computeCombinedGridMomentum()
{
  Vector3s sum = Vector3s::Zero();
  
  for(int k = 0; k < u.nk; ++k) for(int j = 0; j < u.nj; ++j) for(int i = 0; i < u.ni; ++i)
  {
    sum(0) += u_weight_total(i, j, k) * u(i, j, k);
  }
  
  for(int k = 0; k < v.nk; ++k) for(int j = 0; j < v.nj; ++j) for(int i = 0; i < v.ni; ++i)
  {
    sum(1) += v_weight_total(i, j, k) * v(i, j, k);
  }
  
  for(int k = 0; k < w.nk; ++k) for(int j = 0; j < w.nj; ++j) for(int i = 0; i < w.ni; ++i)
  {
    sum(2) += w_weight_total(i, j, k) * w(i, j, k);
  }
  
  return sum;
}

scalar FluidSim3D::computeParticleKineticEnergy()
{
  scalar T = 0.0;
  
  scalar rho = m_parent->getLiquidDensity();
  
  int np = particles.size();
  
  for( int i = 0; i < np; ++i ) {
    auto& p = particles[i];
    scalar m = 4.0 / 3.0 * M_PI * p.radii * p.radii * p.radii * rho;
    T += m * p.v.squaredNorm();
  }
  return 0.5*T;
}

scalar FluidSim3D::computeHairGridKineticEnergy()
{
  scalar sum = 0;
  
  for(int k = 0; k < nk; ++k) for(int j = 0; j < nj; ++j) for(int i = 0; i < ni; ++i)
  {
    Vector3s pos((i+0.5) * dx, (j+0.5) * dx, (k+0.5) * dx);
    Vector3s velocity = get_velocity(pos, u_hair, v_hair, w_hair);
    scalar mass = get_velocity(pos, u_weight_hair, v_weight_hair, w_weight_hair).norm();
    
    sum += mass * velocity.squaredNorm();
  }
  
  return sum * 0.5;
}

scalar FluidSim3D::computeParticleGridKineticEnergy()
{
  scalar sum = 0;
  
  for(int k = 0; k < nk; ++k) for(int j = 0; j < nj; ++j) for(int i = 0; i < ni; ++i)
  {
    Vector3s pos((i+0.5) * dx, (j+0.5) * dx, (k+0.5) * dx);
    Vector3s velocity = get_velocity(pos, u_particle, v_particle, w_particle);
    scalar mass = get_velocity(pos, u_weight_particle, v_weight_particle, w_weight_particle).norm();
    
    sum += mass * velocity.squaredNorm();
  }
  
  return sum * 0.5;
}

scalar FluidSim3D::computeCombinedGridKineticEnergy()
{
  scalar sum = 0;
  
  for(int k = 0; k < nk; ++k) for(int j = 0; j < nj; ++j) for(int i = 0; i < ni; ++i)
  {
    Vector3s pos((i+0.5) * dx, (j+0.5) * dx, (k+0.5) * dx);
    Vector3s velocity = get_velocity(pos);
    scalar mass = get_velocity(pos, u_weight_total, v_weight_total, w_weight_total).norm();
    
    sum += mass * velocity.squaredNorm();
  }
  
  return sum * 0.5;
}

Vector3s FluidSim3D::computeHairGridDrag()
{
  Vector3s F = Vector3s::Zero();
  
  for(int k = 0; k < u.nk; ++k) for(int j = 0; j < u.nj; ++j) for(int i = 0; i < u.ni; ++i)
  {
    F(0) += u_drag(i, j, k) * u_weight_hair(i, j, k);
  }
  
  for(int k = 0; k < v.nk; ++k) for(int j = 0; j < v.nj; ++j) for(int i = 0; i < v.ni; ++i)
  {
    F(1) += v_drag(i, j, k) * v_weight_hair(i, j, k);
  }
  
  for(int k = 0; k < w.nk; ++k) for(int j = 0; j < w.nj; ++j) for(int i = 0; i < w.ni; ++i)
  {
    F(2) += w_drag(i, j, k) * w_weight_hair(i, j, k);
  }
  
  return F;
}

void FluidSim3D::add_particle(const VectorXs& pos, const VectorXs& vel, const scalar& radii, ParticleType type)
{
  add_particle(Particle<3>(pos.segment<3>(0), vel.segment<3>(0), radii, type));
}

void FluidSim3D::controlSources(const scalar& current_time, const scalar& dt)
{
  for(auto& s : sources)
  {
    scalar total_sub_length = s->sub_activate_length + s->sub_inactivate_length;
    s->activated = (current_time >= s->start && current_time + dt < s->end);
    if(total_sub_length != 0.0) {
      int iCycle = (int) (current_time / total_sub_length);
      scalar frac_sub_length = current_time - (scalar) iCycle * total_sub_length;
      s->activated = s->activated && (frac_sub_length + dt < s->sub_activate_length);
    }
  }
}

void FluidSim3D::writeReadable(std::vector<std::ostringstream>& oss) const
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
    
    oss[level] << std::setprecision(15) << p.x.transpose() << " " << p.radii << " " << (int)(p.type == PT_HAIR) << " " << p.fresh << "\n";
  }
}

void FluidSim3D::readReadable( std::ifstream& file )
{
  int np = 0;
  while (!file.eof()) {

    Vector3s pos;
    scalar radius;
    int type;

    file >> pos[0] >> pos[1] >> pos[2] >> radius >> type;
    particles[np].x = pos;
    particles[np].radii = radius;
    ++np;
  }

  particles.resize(np);
}


void FluidSim3D::write(std::vector<scalar>& data) const
{
  for(auto& p : particles)
  {
    p.write(data);
  }
  
  for(auto& b : boundaries)
  {
    if(b->type == BT_UNION || b->type == BT_INTERSECT) continue;
    b->write(data);
  }
  
  data.insert(data.end(), u.a.begin(), u.a.end());
  data.insert(data.end(), v.a.begin(), v.a.end());
  data.insert(data.end(), w.a.begin(), w.a.end());
  
  data.insert(data.end(), u_pressure_grad.a.begin(), u_pressure_grad.a.end());
  data.insert(data.end(), v_pressure_grad.a.begin(), v_pressure_grad.a.end());
  data.insert(data.end(), w_pressure_grad.a.begin(), w_pressure_grad.a.end());
  
  data.insert(data.end(), liquid_phi.a.begin(), liquid_phi.a.end());
}

void FluidSim3D::read(const scalar* data, size_t size_particles, size_t size_boundaries)
{
  
  size_t np = size_particles / Particle<3>::size();
  particles.resize(np);
  
  size_t k = 0;
  for(size_t i = 0; i < np; ++i)
  {
    particles[i].read(data + k);
    k += Particle<3>::size() / sizeof(scalar);
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
  
  memcpy(&u.a[0], data + k, u.a.size() * sizeof(scalar)); k += u.a.size();
  memcpy(&v.a[0], data + k, v.a.size() * sizeof(scalar)); k += v.a.size();
  memcpy(&w.a[0], data + k, w.a.size() * sizeof(scalar)); k += w.a.size();
  
  memcpy(&u_pressure_grad.a[0], data + k, u_pressure_grad.a.size() * sizeof(scalar)); k += u_pressure_grad.a.size();
  memcpy(&v_pressure_grad.a[0], data + k, v_pressure_grad.a.size() * sizeof(scalar)); k += v_pressure_grad.a.size();
  memcpy(&w_pressure_grad.a[0], data + k, w_pressure_grad.a.size() * sizeof(scalar)); k += w_pressure_grad.a.size();
  
  memcpy(&liquid_phi.a[0], data + k, liquid_phi.a.size() * sizeof(scalar)); k += liquid_phi.a.size();
}

size_t FluidSim3D::particle_size() const
{
  return particles.size() * Particle<3>::size();
}

size_t FluidSim3D::boundary_size() const
{
  int sum = 0;
  for(auto& b : boundaries)
  {
    if(b->type == BT_UNION || b->type == BT_INTERSECT) continue;
    sum += b->size();
  }
  
  return sum;
}

size_t FluidSim3D::crucial_grid_size() const
{
  return (u.a.size() + v.a.size() + w.a.size() + u_pressure_grad.a.size() + v_pressure_grad.a.size() + w_pressure_grad.a.size() + liquid_phi.a.size()) * sizeof(scalar);
}

void FluidSim3D::preCompute( const VectorXs& x, const VectorXs& v, const VectorXs& m, const scalar& dt )
{
  int nforces = drag_forces.size();
  threadutils::thread_pool::ParallelFor(0, nforces, [&] (int i) {
    drag_forces[i]->preCompute(x, v, m, dt);
  });
}

void FluidSim3D::addGradEToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, VectorXs& gradE )
{
  int nforces = drag_forces.size();
  threadutils::thread_pool::ParallelFor(0, nforces, [&] (int i) {
    drag_forces[i]->addGradEToTotal(x, v, m, gradE);
  });
}

void FluidSim3D::apply_viscosity(scalar dt) {
  compute_viscosity_weights();
  
  const scalar& visc = m_parent->getViscosity();
  
  advance_viscosity_implicit_weighted(u, v, w, u_visc_impulse, v_visc_impulse, w_visc_impulse,
                                      u_vol_liquid, v_vol_liquid, w_vol_liquid,
                                      c_vol_liquid, ex_vol_liquid, ey_vol_liquid, ez_vol_liquid, cell_solid_phi, visc, dt, dx);
}

void FluidSim3D::compute_viscosity_weights()
{
  estimate_volume_fractions(c_vol_liquid,  Vector3s(0.5*dx, 0.5*dx, 0.5*dx), dx, liquid_phi, Vector3s(0,0,0), dx);
  estimate_volume_fractions(u_vol_liquid,  Vector3s(0,       0.5*dx, 0.5*dx), dx, liquid_phi, Vector3s(0,0,0), dx);
  estimate_volume_fractions(v_vol_liquid,  Vector3s(0.5*dx, 0,       0.5*dx), dx, liquid_phi, Vector3s(0,0,0), dx);
  estimate_volume_fractions(w_vol_liquid,  Vector3s(0.5*dx, 0.5*dx, 0),       dx, liquid_phi, Vector3s(0,0,0), dx);
  estimate_volume_fractions(ex_vol_liquid, Vector3s(0.5*dx, 0,       0),       dx, liquid_phi, Vector3s(0,0,0), dx);
  estimate_volume_fractions(ey_vol_liquid, Vector3s(0,       0.5*dx, 0),       dx, liquid_phi, Vector3s(0,0,0), dx);
  estimate_volume_fractions(ez_vol_liquid, Vector3s(0,       0,       0.5*dx), dx, liquid_phi, Vector3s(0,0,0), dx);
}

void estimate_volume_fractions(Array3s& volumes,
                               const Vector3s& start_centre, const scalar dx,
                               const Array3s& phi, const Vector3s& phi_origin, const scalar phi_dx)
{
  threadutils::thread_pool::ParallelFor(0, volumes.nk, [&] (int k) {
    for(int j = 0; j < volumes.nj; ++j) for(int i = 0; i < volumes.ni; ++i)  {
      Vector3s centre = start_centre + Vector3s(i*dx, j*dx, k*dx);
      
      scalar offset = 0.5*dx;
      
      scalar phi000 = interpolate_value(Vector3s(centre + Vector3s(-offset,-offset,-offset)), phi, phi_origin, phi_dx);
      scalar phi001 = interpolate_value(Vector3s(centre + Vector3s(-offset,-offset,+offset)), phi, phi_origin, phi_dx);
      scalar phi010 = interpolate_value(Vector3s(centre + Vector3s(-offset,+offset,-offset)), phi, phi_origin, phi_dx);
      scalar phi011 = interpolate_value(Vector3s(centre + Vector3s(-offset,+offset,+offset)), phi, phi_origin, phi_dx);
      scalar phi100 = interpolate_value(Vector3s(centre + Vector3s(+offset,-offset,-offset)), phi, phi_origin, phi_dx);
      scalar phi101 = interpolate_value(Vector3s(centre + Vector3s(+offset,-offset,+offset)), phi, phi_origin, phi_dx);
      scalar phi110 = interpolate_value(Vector3s(centre + Vector3s(+offset,+offset,-offset)), phi, phi_origin, phi_dx);
      scalar phi111 = interpolate_value(Vector3s(centre + Vector3s(+offset,+offset,+offset)), phi, phi_origin, phi_dx);
      
      volumes(i,j,k) = volume_fraction(phi000, phi100, phi010, phi110, phi001, phi101, phi011, phi111);
    }
  });
}

//Apply several iterations of a very simple "Jacobi"-style propagation of valid velocity data in all directions
void extrapolate(Array3s& grid, Array3s& old_grid, const Array3s& grid_weight, const Array3s& grid_liquid_weight, Array3c& valid, Array3c old_valid, const Vector3i& offset) {
  
  //Initialize the list of valid cells
  
  for(int k = 0; k < valid.nk; ++k) for(int j = 0; j < valid.nj; ++j) valid(0,j,k) = valid(valid.ni-1,j,k) = 0;
  for(int k = 0; k < valid.nk; ++k) for(int i = 0; i < valid.ni; ++i) valid(i,0,k) = valid(i,valid.nj-1,k) = 0;
  for(int j = 0; j < valid.nj; ++j) for(int i = 0; i < valid.ni; ++i) valid(i,j,0) = valid(i,j,valid.nk-1) = 0;
  
  for(int k = 1; k < grid.nk - 1; ++k) for(int j = 1; j < grid.nj - 1; ++j) for(int i = 1; i < grid.ni - 1; ++i)
    valid(i,j,k) = grid_weight(i,j,k) > 0 && (grid_liquid_weight(i, j, k) < 0 || grid_liquid_weight(i + offset(0), j + offset(1), k + offset(2)) < 0 );
  
  Array3s* pgrid[2] = {&grid, &old_grid};
  Array3c* pvalid[2] = {&valid, &old_valid};
  
  for(int layers = 0; layers < 4; ++layers) {
    int num = grid.ni*grid.nj*grid.nk;
    int slice =  grid.ni*grid.nj;
    Array3s* pgrid_source = pgrid[layers & 1];
    Array3s* pgrid_target = pgrid[!(layers & 1)];
    
    Array3c* pvalid_source = pvalid[layers & 1];
    Array3c* pvalid_target = pvalid[!(layers & 1)];
    
    threadutils::thread_pool::ParallelFor(0,num, [&](int thread_idx)
    {
      int k = thread_idx/slice;
      int j = (thread_idx%slice)/grid.ni;
      int i = thread_idx%grid.ni;
      //for(int k = 1; k < grid.nk-1; ++k) for(int j = 1; j < grid.nj-1; ++j) for(int i = 1; i < grid.ni-1; ++i) {
      if(i>0&&i<grid.ni-1 &&j>0&&j<grid.nj-1 &&k>0&&k<grid.nk-1)
      {
        scalar sum = 0;
        int count = 0;
        
        if(!(*pvalid_source)(i,j,k)) {
          
          if((*pvalid_source)(i+1,j,k)) {
            sum += (*pgrid_source)(i+1,j,k);
            ++count;
          }
          if((*pvalid_source)(i-1,j,k)) {
            sum += (*pgrid_source)(i-1,j,k);
            ++count;
          }
          if((*pvalid_source)(i,j+1,k)) {
            sum += (*pgrid_source)(i,j+1,k);
            ++count;
          }
          if((*pvalid_source)(i,j-1,k)) {
            sum += (*pgrid_source)(i,j-1,k);
            ++count;
          }
          if((*pvalid_source)(i,j,k+1)) {
            sum += (*pgrid_source)(i,j,k+1);
            ++count;
          }
          if((*pvalid_source)(i,j,k-1)) {
            sum += (*pgrid_source)(i,j,k-1);
            ++count;
          }
          
          //If any of neighbour cells were valid,
          //assign the cell their average value and tag it as valid
          if(count > 0) {
            (*pgrid_target)(i,j,k) = sum / (scalar)count;
            (*pvalid_target)(i,j,k) = 1;
          }
        }
      }
    });
    
    *pvalid_source = *pvalid_target;
    *pgrid_source = *pgrid_target;
  }
}
                    

