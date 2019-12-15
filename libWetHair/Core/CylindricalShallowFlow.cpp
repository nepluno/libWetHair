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


#include "CylindricalShallowFlow.h"
#include "TwoDScene.h"
#include "MathUtilities.h"
#include "HairFlow.h"
#include "fluidsim2D.h"
#include "fluidsim3D.h"
#include "liangbarsky.h"
#include "bresemham.h"
#include "array2_utils.h"
#include "PolygonalCohesion.h"
#include <vector>

static const char* cylindrical_film_static_name = "cylindricalshallowflow";

using namespace mathutils;
using namespace mathutils::shallowflow;

const int N_GAUSS = 2;

template<int DIM>
CylindricalShallowFlow<DIM>::CylindricalShallowFlow(TwoDScene<DIM>* parent, const std::vector<int>& involved_particles, const VectorXs& eta, const std::vector<unsigned char>& particle_state) :
HairFlow<DIM>(parent, involved_particles, eta, particle_state)
{
  int np = involved_particles.size();
  for(int i = 0; i < np; ++i)
  {
    HairFlow<DIM>::m_global_to_local[involved_particles[i]] = i;
  }
  
  HairFlow<DIM>::m_particle_to_edges.resize(np);
  auto& edges = HairFlow<DIM>::m_parent->getEdges();
  int ne = edges.size();
  for(int i = 0; i < ne; ++i)
  {
    auto& e = edges[i];
    auto itr0 = HairFlow<DIM>::m_global_to_local.find(e.first);
    auto itr1 = HairFlow<DIM>::m_global_to_local.find(e.second);
    if(itr0 != HairFlow<DIM>::m_global_to_local.end() && itr1 != HairFlow<DIM>::m_global_to_local.end())
    {
      HairFlow<DIM>::m_edge_indices.push_back(i);
      HairFlow<DIM>::m_global_edges.push_back(e);
      HairFlow<DIM>::m_internal_edges.push_back(std::pair<int, int>(itr0->second, itr1->second));
      HairFlow<DIM>::m_particle_to_edges[itr0->second].push_back(i);
      HairFlow<DIM>::m_particle_to_edges[itr1->second].push_back(i);
    }
  }
  
  resizeSystem();
  HairFlow<DIM>::m_u.setZero();
  m_uvert.setZero();
  HairFlow<DIM>::m_rad_vec.setZero();
  HairFlow<DIM>::m_area_v.setZero();
  
  HairFlow<DIM>::m_porosity.setZero();
  HairFlow<DIM>::m_c_v.setZero();
  HairFlow<DIM>::m_c_star.setZero();
  
  compute_SP0_matrix(m_sp0);
  compute_SP1_matrix(m_sp1);

  copy_rest_mass<DIM>(HairFlow<DIM>::m_parent->getHairRestMass(), HairFlow<DIM>::m_particle_indices, m_mass_v, HairFlow<DIM>::m_parent);
  
  compute_edge_val(HairFlow<DIM>::m_eta, HairFlow<DIM>::m_internal_edges, HairFlow<DIM>::m_edge_eta);
  
  compute_rest_porosity(HairFlow<DIM>::m_eta, HairFlow<DIM>::m_parent->getRadii(), involved_particles, HairFlow<DIM>::m_porosity);
  m_pool_liquid = 0.0;
}

template<int DIM>
void CylindricalShallowFlow<DIM>::updateGeometricState(const VectorXs& x, const VectorXs& v, FluidSim* fluidsim)
{
  // compute geometry
  m_old_rad_vec = HairFlow<DIM>::m_rad_vec;
  m_old_area_v = HairFlow<DIM>::m_area_v; 
  const PolygonalCohesion<DIM>* cohesion = HairFlow<DIM>::m_parent->getPolygonalCohesion();

  const int ne = HairFlow<DIM>::m_edge_indices.size();
  const int np = HairFlow<DIM>::m_particle_indices.size();
  const std::vector<int>& ppp_count_global = cohesion->getPPPCountV();
  
  
  compute_radius_vector(HairFlow<DIM>::m_parent->getEdgeRadii(), HairFlow<DIM>::m_edge_indices, HairFlow<DIM>::m_edge_rad_vec);
  
  if(HairFlow<DIM>::m_parent->isIndividualTransfer())
  {
    compute_edge_area<DIM>(x, HairFlow<DIM>::m_global_edges, HairFlow<DIM>::m_area_e, HairFlow<DIM>::m_parent);
    
    compute_vertex_area<DIM>(HairFlow<DIM>::m_internal_edges, HairFlow<DIM>::m_area_e, HairFlow<DIM>::m_area_v_hair);
    
    HairFlow<DIM>::m_area_v = HairFlow<DIM>::m_area_v_hair;
  } else {
    const VectorXs& area_v_global = cohesion->getAreaVGlobal();
    const VectorXs& area_e_global = cohesion->getAreaEGlobal();
    
    for(int i = 0; i < ne; ++i) HairFlow<DIM>::m_area_e(i) = area_e_global(HairFlow<DIM>::m_edge_indices[i]);
    
    for(int i = 0; i < np; ++i) HairFlow<DIM>::m_area_v(i) = area_v_global(HairFlow<DIM>::m_particle_indices[i]);
    
    compute_vertex_area<DIM>(HairFlow<DIM>::m_internal_edges, HairFlow<DIM>::m_area_e, HairFlow<DIM>::m_area_v_hair);
  }
  

  
  // compute edge-vertex weight
  compute_bracket_matrix(HairFlow<DIM>::m_area_e, m_Gf, DIM);
  compute_bracket_matrix(HairFlow<DIM>::m_area_v, m_Gv);
  compute_inv_bracket_matrix(HairFlow<DIM>::m_area_v, m_iGv);
  
  compute_bracket_matrix(HairFlow<DIM>::m_area_v_hair, m_Gv_hair);
  compute_flat_bracket_matrix(HairFlow<DIM>::m_area_v_hair, m_area_v_hair_flat);
  compute_inv_bracket_matrix(HairFlow<DIM>::m_area_v_hair, m_iGv_hair);
  
  // compute local interpolation matrix
  compute_edge_to_vertex_matrix(HairFlow<DIM>::m_internal_edges, HairFlow<DIM>::m_area_v, HairFlow<DIM>::m_area_e, m_W_fv);
  compute_edge_to_vertex_matrix(HairFlow<DIM>::m_internal_edges, HairFlow<DIM>::m_area_v_hair, HairFlow<DIM>::m_area_e, m_W_fv_hair);
  
  // compute dirs
  compute_edge_dir<DIM>(x, HairFlow<DIM>::m_global_edges, HairFlow<DIM>::m_dir_f, HairFlow<DIM>::m_parent);
  expand_dir_f(HairFlow<DIM>::m_dir_f, m_dir_f_expand);
  compute_vertex_val(HairFlow<DIM>::m_dir_f, m_W_fv_hair, HairFlow<DIM>::m_dir_v);
  compute_normalized_vector<DIM>(HairFlow<DIM>::m_dir_v);
  
  compute_vertex_val(HairFlow<DIM>::m_edge_rad_vec, m_W_fv_hair, HairFlow<DIM>::m_rad_vec);
    // compute local gradient matrix
  compute_edge_gradient_matrix<DIM>(x, HairFlow<DIM>::m_global_edges, HairFlow<DIM>::m_internal_edges, HairFlow<DIM>::m_area_e, m_gradF, HairFlow<DIM>::m_parent);
  
  // compute local divergence matrix
  compute_vertex_divergence_matrix(m_iGv, m_Gf, m_gradF, m_divV);
  
  // compute local divergence matrix
  compute_vertex_divergence_matrix(m_iGv_hair, m_Gf, m_gradF, m_divV_hair);
  // compute local laplacian matrix
  compute_vertex_laplacian_matrix(m_gradF, m_divV_hair, m_L);
  
  // compute local coordinates
  m_sum_area_e =
  compute_accumulated_edge_length(HairFlow<DIM>::m_area_e, m_area_e_accu);
  compute_inverted_pos_mapping(m_area_e_accu, m_sum_area_e, m_area_e_inv_mapping);
  compute_edge_val(m_area_e_accu, HairFlow<DIM>::m_internal_edges, m_edge_x);
  
  // compute local porosity
  // compute_rest_porosity(HairFlow<DIM>::m_eta, HairFlow<DIM>::m_parent->getRadii(), HairFlow<DIM>::m_particle_indices, HairFlow<DIM>::m_porosity);
  compute_edge_val(HairFlow<DIM>::m_porosity, HairFlow<DIM>::m_internal_edges, m_porosity_e);
  HairFlow<DIM>::m_avg_area_e =
  compute_avg(HairFlow<DIM>::m_area_e);
  
  HairFlow<DIM>::m_max_area_e = HairFlow<DIM>::m_area_e.maxCoeff();
  HairFlow<DIM>::m_min_area_e = HairFlow<DIM>::m_area_e.minCoeff();
  
  HairFlow<DIM>::m_max_area_v = HairFlow<DIM>::m_area_v.maxCoeff();
  HairFlow<DIM>::m_min_area_v = HairFlow<DIM>::m_area_v.minCoeff();
  
}

template<int DIM>
void CylindricalShallowFlow<DIM>::updateReservoir(FluidSim* fluidsim, const VectorXs& x, const VectorXs& v, const scalar& dt)
{
  if(fluidsim) {
    if(HairFlow<DIM>::m_parent->isIndividualTransfer()) {
      int ne = HairFlow<DIM>::m_area_e.size();
      for(int i = 0; i < ne; ++i)
      {
        int p0 = i;
        int p1 = i + 1;
        
        scalar H0 = HairFlow<DIM>::m_eta(p0) + HairFlow<DIM>::m_rad_vec(p0);
        scalar H1 = HairFlow<DIM>::m_eta(p1) + HairFlow<DIM>::m_rad_vec(p1);
        
        scalar H0_old = m_old_eta(p0) + HairFlow<DIM>::m_rad_vec(p0);
        scalar H1_old = m_old_eta(p1) + HairFlow<DIM>::m_rad_vec(p1);
        
        scalar dA0 = H0 * H0 - H0_old * H0_old;
        scalar dA1 = H1 * H1 - H1_old * H1_old;
        
        m_pool_liquid += M_PI * (- dA0 - dA1) * HairFlow<DIM>::m_area_e(i) * 0.5;
        
        // std::cout << "[SGP: " << pidx0 << ", " << dA0 << ", " << Q0 << "]" << std::endl;
      }
    } else {
      const PolygonalCohesion<DIM>* cohesion = HairFlow<DIM>::m_parent->getPolygonalCohesion();
      const MatrixXs& inter_rhs = cohesion->getRhsOffsetVGlobal();
      
      int ne = HairFlow<DIM>::m_area_e.size();
      for(int i = 0; i < ne; ++i)
      {
        int p0 = i;
        int p1 = i + 1;
        
        int pidx0 = HairFlow<DIM>::m_particle_indices[p0];
        int pidx1 = HairFlow<DIM>::m_particle_indices[p1];
        
        scalar H0 = HairFlow<DIM>::m_eta(p0) + HairFlow<DIM>::m_rad_vec(p0);
        scalar H1 = HairFlow<DIM>::m_eta(p1) + HairFlow<DIM>::m_rad_vec(p1);
        
        scalar H0_old = m_old_eta(p0) + HairFlow<DIM>::m_rad_vec(p0);
        scalar H1_old = m_old_eta(p1) + HairFlow<DIM>::m_rad_vec(p1);
        
        scalar dA0 = H0 * H0 - H0_old * H0_old;
        scalar dA1 = H1 * H1 - H1_old * H1_old;
        
        const scalar& Q0 = inter_rhs(pidx0, 0);
        const scalar& Q1 = inter_rhs(pidx1, 0);
        
        m_pool_liquid += M_PI * (Q0 + Q1 - dA0 - dA1) * HairFlow<DIM>::m_area_e(i) * 0.5;
        
        // std::cout << "[SGP: " << pidx0 << ", " << dA0 << ", " << Q0 << "]" << std::endl;
      }
    }
    
    m_pool_liquid = std::max(0.0, m_pool_liquid);
  }
}


template<>
void CylindricalShallowFlow<2>::updateToFilteredGrid(const VectorXs& x, const VectorXs& v, FluidSim* fluidsim, const scalar& dt, int ibuffer)
{
  FluidSim2D* fluid2d = (FluidSim2D*) fluidsim;
  scalar cellsize = fluid2d->cellsize();
  
  const VectorXs& drag_buffer = m_parent->getFluidDragBuffer();
  std::vector<EdgeVelDragIntersection<2> >& u_edge_vel_drag = fluid2d->get_u_edge_vel_drag()[ibuffer];
  std::vector<EdgeVelDragIntersection<2> >& v_edge_vel_drag = fluid2d->get_v_edge_vel_drag()[ibuffer];
  
  u_edge_vel_drag.resize(0);
  v_edge_vel_drag.resize(0);

  const scalar& rho = m_parent->getLiquidDensity();
  
  // for each edge
  int ne = HairFlow<2>::m_global_edges.size();
  for(int eidx = 0; eidx < ne; ++eidx)
  {
    auto& e = HairFlow<2>::m_global_edges[eidx];
    int local_first = m_global_to_local.find(e.first)->second;
    int local_second = m_global_to_local.find(e.second)->second;
    
    const scalar eta_first = m_eta(local_first);
    const scalar eta_second = m_eta(local_second);
    const scalar radii_first = m_rad_vec(local_first);
    const scalar radii_second = m_rad_vec(local_second);
    const scalar er_first = eta_first + radii_first;
    const scalar er_second = eta_second + eta_second;
    const scalar rho_h_first = m_parent->getHairDensity(e.first);
    const scalar rho_h_second = m_parent->getHairDensity(e.second);
    
    const Vector2s& u_first = m_actual_u_v.row(local_first).transpose();
    const Vector2s& u_second = m_actual_u_v.row(local_second).transpose();
    
    const Vector2s& x0 = x.segment<2>(e.first * 2);
    const Vector2s& x1 = x.segment<2>(e.second * 2);
    const Vector2s& v0 = v.segment<2>(e.first * 2);
    const Vector2s& v1 = v.segment<2>(e.second * 2);
    const Vector2s& d0 = drag_buffer.segment<2>(e.first * 2);
    const Vector2s& d1 = drag_buffer.segment<2>(e.second * 2);
    const Vector4s& c0 = m_c_v.block<1,4>(local_first,0).transpose();
    const Vector4s& c1 = m_c_v.block<1,4>(local_second,0).transpose();
    
    scalar xmin;
    scalar xmax;
    scalar ymin;
    scalar ymax;
    
    if(x0(0) - er_first < x1(0) - er_second) {
      xmin = x0(0) - std::max(cellsize * 2.0, er_first);
    } else {
      xmin = x1(0) - std::max(cellsize * 2.0, er_second);
    }
    
    if(x0(0) + er_first > x1(0) + er_second) {
      xmax = x0(0) + std::max(cellsize * 2.0, er_first);
    } else {
      xmax = x1(0) + std::max(cellsize * 2.0, er_second);
    }
    
    if(x0(1) - er_first < x1(1) - er_second) {
      ymin = x0(1) - std::max(cellsize * 2.0, er_first);
    } else {
      ymin = x1(1) - std::max(cellsize * 2.0, er_second);
    }
    
    if(x0(1) + er_first > x1(1) + er_second) {
      ymax = x0(1) + std::max(cellsize * 2.0, er_first);
    } else {
      ymax = x1(1) + std::max(cellsize * 2.0, er_second);
    }
    
    int xmin_u = std::max(0, std::min(fluid2d->get_u_ni()-1, (int) ceil((xmin - fluid2d->get_origin()(0)) / cellsize)));
    int xmin_v = std::max(0, std::min(fluid2d->get_v_ni()-1, (int) ceil((xmin - fluid2d->get_origin()(0)) / cellsize - 0.5)));
    int xmax_u = std::max(0, std::min(fluid2d->get_u_ni()-1, (int) floor((xmax - fluid2d->get_origin()(0)) / cellsize)));
    int xmax_v = std::max(0, std::min(fluid2d->get_v_ni()-1, (int) floor((xmax - fluid2d->get_origin()(0)) / cellsize - 0.5)));
    
    int ymin_u = std::max(0, std::min(fluid2d->get_u_nj()-1, (int) ceil((ymin - fluid2d->get_origin()(1)) / cellsize - 0.5)));
    int ymin_v = std::max(0, std::min(fluid2d->get_v_nj()-1, (int) ceil((ymin - fluid2d->get_origin()(1)) / cellsize)));
    int ymax_u = std::max(0, std::min(fluid2d->get_u_nj()-1, (int) floor((ymax - fluid2d->get_origin()(1)) / cellsize - 0.5)));
    int ymax_v = std::max(0, std::min(fluid2d->get_v_nj()-1, (int) floor((ymax - fluid2d->get_origin()(1)) / cellsize)));
    
    std::vector<Vector2i> line_buffer;
    for(int j = ymin_u; j <= ymax_u; ++j) for(int i = xmin_u; i <= xmax_u; ++i)
    {
      Vector2s pos = Vector2s(i*cellsize,(j+0.5)*cellsize) + fluid2d->get_origin();
      Vector4s clipping_window(pos(0) - cellsize, pos(1) - cellsize,
                               pos(0) + cellsize, pos(1) + cellsize);
      
      Vector2s q0 = x0, q1 = x1;
      scalar alpha0 = 0.0, alpha1 = 1.0;
      if(!liangbarsky::clip_line(clipping_window, q0, q1, alpha0, alpha1)) continue;
      
      scalar length_segment = (q1 - q0).norm();
      
      for(int k = 0; k < N_GAUSS; ++k)
      {
        const scalar beta = mathutils::gauss_legendre_point[N_GAUSS-1][k];
        const scalar w = mathutils::gauss_legendre_weight[N_GAUSS-1][k];
        
        scalar alpha = mathutils::lerp(alpha0, alpha1, beta);
        Vector2s xc = x0 * (1.0 - alpha) + x1 * alpha;
        Vector2s vel_c = v0 * (1.0 - alpha) + v1 * alpha;
        Vector2s xcb = xc - vel_c * dt;
        Vector2i ixcb = Vector2i( (int)((xcb(0) - fluid2d->get_origin()(0)) / cellsize), (int)((xcb(1) - fluid2d->get_origin()(1)) / cellsize - 0.5));
        
        line_buffer.resize(0);
        bresemham::Bresenham2D(ixcb(0), ixcb(1), i, j, [&] (int ii, int jj) {
          line_buffer.push_back(Vector2i(ii, jj));
        });
        if(line_buffer.size() == 0) line_buffer.push_back(Vector2i(i, j));
        
        scalar ubase_c = vel_c(0);
        scalar uflow_c = u_first(0) * (1.0 - alpha) + u_second(0) * alpha;
        scalar vol_c = M_PI * mathutils::lerp(er_first * er_first - radii_first * radii_first, er_second * er_second - radii_second * radii_second, alpha) * length_segment;
        scalar mass_c = rho * vol_c;
        scalar rho_base_c = mathutils::lerp(rho_h_first, rho_h_second, alpha);
        scalar vol_base_c = M_PI * mathutils::lerp(radii_first * radii_first, radii_second * radii_second, alpha) * length_segment;
        scalar mass_base_c = rho_base_c * vol_base_c;
        scalar dragc = d0(0) * (1.0 - alpha) + d1(0) * alpha;
        Vector2s c_c = c0.segment<2>(0) * (1.0 - alpha) + c1.segment<2>(0) * alpha;
        
        Vector2s diff = pos - xc;
        
        scalar w_linear = mathutils::linear_kernel(diff, cellsize);
        scalar weight = (mass_c + mass_base_c) * w_linear;
        
        scalar u_base_world = ubase_c;
        scalar u_flow_world = u_base_world + uflow_c + c_c.dot(pos - xc);
        
        for(const Vector2i& cd : line_buffer) {
          EdgeVelDragIntersection<2> tmp_inter;
          tmp_inter.vel_weighted = (u_base_world * mass_base_c + u_flow_world * mass_c) * w_linear * w;
          tmp_inter.drag_weighted = dragc * weight * w / (scalar) line_buffer.size();
          tmp_inter.vol_weighted = (vol_c + vol_base_c) * w_linear * w;
          tmp_inter.weight = weight * w;
          tmp_inter.linear_weight = w_linear * w;
          tmp_inter.coord = cd;
          u_edge_vel_drag.push_back(tmp_inter);
        }
      }
    }
    
    for(int j = ymin_v; j <= ymax_v; ++j) for(int i = xmin_v; i <= xmax_v; ++i)
    {
      Vector2s pos = Vector2s((i+0.5)*cellsize,j*cellsize) + fluid2d->get_origin();
      Vector4s clipping_window(pos(0) - cellsize, pos(1) - cellsize,
                               pos(0) + cellsize, pos(1) + cellsize);
      Vector2s q0 = x0, q1 = x1;
      scalar alpha0 = 0.0, alpha1 = 1.0;
      if(!liangbarsky::clip_line(clipping_window, q0, q1, alpha0, alpha1)) continue;
      
      scalar length_segment = (q1 - q0).norm();
      
      for(int k = 0; k < N_GAUSS; ++k)
      {
        const scalar beta = mathutils::gauss_legendre_point[N_GAUSS-1][k];
        const scalar w = mathutils::gauss_legendre_weight[N_GAUSS-1][k];
        
        scalar alpha = mathutils::lerp(alpha0, alpha1, beta);
        Vector2s xc = x0 * (1.0 - alpha) + x1 * alpha;
        Vector2s vel_c = v0 * (1.0 - alpha) + v1 * alpha;
        Vector2s xcb = xc - vel_c * dt;
        Vector2i ixcb = Vector2i( (int)((xcb(0) - fluid2d->get_origin()(0)) / cellsize - 0.5), (int)((xcb(1) - fluid2d->get_origin()(1)) / cellsize));
        
        line_buffer.resize(0);
        bresemham::Bresenham2D(ixcb(0), ixcb(1), i, j, [&] (int ii, int jj) {
          line_buffer.push_back(Vector2i(ii, jj));
        });
        if(line_buffer.size() == 0) line_buffer.push_back(Vector2i(i, j));
        
        scalar vbase_c = vel_c(1);
        scalar vflow_c = u_first(1) * (1.0 - alpha) + u_second(1) * alpha;
        scalar vol_c = M_PI * mathutils::lerp(er_first * er_first - radii_first * radii_first, er_second * er_second - radii_second * radii_second, alpha) * length_segment;
        scalar mass_c = rho * vol_c;
        scalar rho_base_c = mathutils::lerp(rho_h_first, rho_h_second, alpha);
        scalar vol_base_c = M_PI * mathutils::lerp(radii_first * radii_first, radii_second * radii_second, alpha) * length_segment;
        scalar mass_base_c = rho_base_c * vol_base_c;
        scalar dragc = d0(1) * (1.0 - alpha) + d1(1) * alpha;
        Vector2s c_c = c0.segment<2>(2) * (1.0 - alpha) + c1.segment<2>(2) * alpha;
        
        Vector2s diff = pos - xc;
        
        scalar w_linear = mathutils::linear_kernel(diff, fluid2d->cellsize());
        scalar weight = (mass_c + mass_base_c) * w_linear;
        
        scalar v_base_world = vbase_c;
        scalar v_flow_world = v_base_world + vflow_c + c_c.dot(pos - xc);

        for(const Vector2i& cd : line_buffer) {
          EdgeVelDragIntersection<2> tmp_inter;
          tmp_inter.vel_weighted = (v_base_world * mass_base_c + v_flow_world * mass_c) * w_linear * w;
          tmp_inter.vol_weighted = (vol_c + vol_base_c) * w_linear * w;
          tmp_inter.weight = weight * w;
          tmp_inter.linear_weight = w_linear * w;
          tmp_inter.coord = cd;
          
          v_edge_vel_drag.push_back(tmp_inter);
        }
      }
    }
  }
}

template<>
Vector3s CylindricalShallowFlow<2>::computeHairLiquidAngularMomentum(const VectorXs& x, const VectorXs& v, FluidSim* fluidsim) const
{
  FluidSim2D* fluid2d = (FluidSim2D*) fluidsim;
  scalar cellsize = fluid2d->cellsize();
  
  VectorXs gradE;
  gradE.resize(x.size());
  
  scalar sum = 0.0;
  
  const scalar& rho = m_parent->getLiquidDensity();
  // for each edge
  for(auto& e : HairFlow<2>::m_global_edges)
  {
    int local_first = m_global_to_local.find(e.first)->second;
    int local_second = m_global_to_local.find(e.second)->second;
    
    const scalar eta_first = m_eta(local_first);
    const scalar eta_second = m_eta(local_second);
    const scalar radii_first = m_rad_vec(local_first);
    const scalar radii_second = m_rad_vec(local_second);
    const scalar er_first = eta_first + radii_first;
    const scalar er_second = eta_second + eta_second;
    const scalar rho_h_first = m_parent->getHairDensity(e.first);
    const scalar rho_h_second = m_parent->getHairDensity(e.second);
    
    const scalar u_first = m_uvert(local_first);
    const scalar u_second = m_uvert(local_second);
    
    const Vector2s& x0 = x.segment<2>(e.first * 2);
    const Vector2s& x1 = x.segment<2>(e.second * 2);
    const Vector2s& v0 = v.segment<2>(e.first * 2);
    const Vector2s& v1 = v.segment<2>(e.second * 2);
    const Vector4s& c0 = m_c_v.block<1,4>(local_first,0).transpose();
    const Vector4s& c1 = m_c_v.block<1,4>(local_second,0).transpose();
    
    Vector2s dir = (x1 - x0).normalized();
    
    scalar xmin;
    scalar xmax;
    scalar ymin;
    scalar ymax;
    
    if(x0(0) - eta_first < x1(0) - eta_second) {
      xmin = x0(0) - std::max(fluid2d->cellsize() * 2.0, eta_first);
    } else {
      xmin = x1(0) - std::max(fluid2d->cellsize() * 2.0, eta_second);
    }
    
    if(x0(0) + eta_first > x1(0) + eta_second) {
      xmax = x0(0) + std::max(fluid2d->cellsize() * 2.0, eta_first);
    } else {
      xmax = x1(0) + std::max(fluid2d->cellsize() * 2.0, eta_second);
    }
    
    if(x0(1) - eta_first < x1(1) - eta_second) {
      ymin = x0(1) - std::max(fluid2d->cellsize() * 2.0, eta_first);
    } else {
      ymin = x1(1) - std::max(fluid2d->cellsize() * 2.0, eta_second);
    }
    
    if(x0(1) + eta_first > x1(1) + eta_second) {
      ymax = x0(1) + std::max(fluid2d->cellsize() * 2.0, eta_first);
    } else {
      ymax = x1(1) + std::max(fluid2d->cellsize() * 2.0, eta_second);
    }
    
    int xmin_u = std::max(0, std::min(fluid2d->get_u_ni()-1, (int) ceil((xmin - fluid2d->get_origin()(0)) / cellsize)));
    int xmin_v = std::max(0, std::min(fluid2d->get_v_ni()-1, (int) ceil((xmin - fluid2d->get_origin()(0)) / cellsize - 0.5)));
    int xmax_u = std::max(0, std::min(fluid2d->get_u_ni()-1, (int) floor((xmax - fluid2d->get_origin()(0)) / cellsize)));
    int xmax_v = std::max(0, std::min(fluid2d->get_v_ni()-1, (int) floor((xmax - fluid2d->get_origin()(0)) / cellsize - 0.5)));
    
    int ymin_u = std::max(0, std::min(fluid2d->get_u_nj()-1, (int) ceil((ymin - fluid2d->get_origin()(1)) / cellsize - 0.5)));
    int ymin_v = std::max(0, std::min(fluid2d->get_v_nj()-1, (int) ceil((ymin - fluid2d->get_origin()(1)) / cellsize)));
    int ymax_u = std::max(0, std::min(fluid2d->get_u_nj()-1, (int) floor((ymax - fluid2d->get_origin()(1)) / cellsize - 0.5)));
    int ymax_v = std::max(0, std::min(fluid2d->get_v_nj()-1, (int) floor((ymax - fluid2d->get_origin()(1)) / cellsize)));
    
    for(int j = ymin_u; j <= ymax_u; ++j) for(int i = xmin_u; i <= xmax_u; ++i)
    {
      Vector2s pos = Vector2s(i*cellsize,(j+0.5)*cellsize) + fluid2d->get_origin();
      Vector4s clipping_window(pos(0) - cellsize, pos(1) - cellsize,
                               pos(0) + cellsize, pos(1) + cellsize);
      
      Vector2s q0 = x0, q1 = x1;
      scalar alpha0, alpha1;
      if(!liangbarsky::clip_line(clipping_window, q0, q1, alpha0, alpha1)) continue;
      
      scalar length_segment = (q1 - q0).norm();
      
      for(int k = 0; k < N_GAUSS; ++k)
      {
        const scalar beta = mathutils::gauss_legendre_point[N_GAUSS-1][k];
        const scalar w = mathutils::gauss_legendre_weight[N_GAUSS-1][k];
        
        scalar alpha = mathutils::lerp(alpha0, alpha1, beta);
        Vector2s xc = x0 * (1.0 - alpha) + x1 * alpha;
        scalar ubase_c = v0(0) * (1.0 - alpha) + v1(0) * alpha;
        scalar uflow_c = u_first * (1.0 - alpha) + u_second * alpha;
        scalar mass_c = rho * M_PI * mathutils::lerp(er_first * er_first - radii_first * radii_first, er_second * er_second - radii_second * radii_second, alpha) * length_segment;
        scalar rho_base_c = mathutils::lerp(rho_h_first, rho_h_second, alpha);
        scalar mass_base_c = rho_base_c * M_PI * mathutils::lerp(radii_first * radii_first, radii_second * radii_second, alpha) * length_segment;

        Vector2s c_c = c0.segment<2>(0) * (1.0 - alpha) + c1.segment<2>(0) * alpha;
        
        Vector2s diff = pos - xc;
        
        scalar w_linear = mathutils::linear_kernel(diff, cellsize);
        scalar u_base_world = ubase_c;
        scalar u_flow_world = u_base_world + uflow_c + c_c.dot(pos - xc);
        
        Vector2s pu = Vector2s((u_base_world * mass_base_c + u_flow_world * mass_c) * w_linear * w, 0.0);
        sum += cross2(pos, pu);
      }
    }
    
    for(int j = ymin_v; j <= ymax_v; ++j) for(int i = xmin_v; i <= xmax_v; ++i)
    {
      Vector2s pos = Vector2s((i+0.5)*cellsize,j*cellsize) + fluid2d->get_origin();
      Vector4s clipping_window(pos(0) - cellsize, pos(1) - cellsize,
                               pos(0) + cellsize, pos(1) + cellsize);
      Vector2s q0 = x0, q1 = x1;
      scalar alpha0, alpha1;
      if(!liangbarsky::clip_line(clipping_window, q0, q1, alpha0, alpha1)) continue;
      
      scalar length_segment = (q1 - q0).norm();
      
      for(int k = 0; k < N_GAUSS; ++k)
      {
        const scalar beta = mathutils::gauss_legendre_point[N_GAUSS-1][k];
        const scalar w = mathutils::gauss_legendre_weight[N_GAUSS-1][k];
        
        scalar alpha = mathutils::lerp(alpha0, alpha1, beta);
        Vector2s xc = x0 * (1.0 - alpha) + x1 * alpha;
        
        scalar vbase_c = v0(1) * (1.0 - alpha) + v1(1) * alpha;
        scalar vflow_c = u_first * (1.0 - alpha) + u_second * alpha;
        scalar mass_c = rho * M_PI * mathutils::lerp(er_first * er_first - radii_first * radii_first, er_second * er_second - radii_second * radii_second, alpha) * length_segment;
        scalar rho_base_c = mathutils::lerp(rho_h_first, rho_h_second, alpha);
        scalar mass_base_c = rho_base_c * M_PI * mathutils::lerp(radii_first * radii_first, radii_second * radii_second, alpha) * length_segment;
        Vector2s c_c = c0.segment<2>(2) * (1.0 - alpha) + c1.segment<2>(2) * alpha;
        
        Vector2s diff = pos - xc;
        
        scalar w_linear = mathutils::linear_kernel(diff, cellsize);
        scalar v_base_world = vbase_c;
        scalar v_flow_world = v_base_world + vflow_c + c_c.dot(pos - xc);
        
        Vector2s pv = Vector2s(0.0, (v_base_world * mass_base_c + v_flow_world * mass_c) * w_linear * w);
        sum += cross2(pos, pv);
      }
    }
  }
  
  return Vector3s(0,0,sum);
}

template<>
void CylindricalShallowFlow<2>::updateFromFilteredGrid(const VectorXs& x, VectorXs& v, FluidSim* fluidsim, const scalar& dt)
{
  FluidSim2D* fluid2d = (FluidSim2D*) fluidsim;
  
  // for each edge
  int ne = HairFlow<2>::m_global_edges.size();
  int np = HairFlow<2>::m_particle_indices.size();
  
  const scalar rho_L = m_parent->getLiquidDensity();
  
  for(int i = 0; i < ne; ++i)
  {
    auto& e = HairFlow<2>::m_global_edges[i];
    
    const Vector2s& x0 = x.segment<2>( m_parent->getDof( e.first ) );
    const Vector2s& x1 = x.segment<2>( m_parent->getDof( e.second ) );
    
    Vector2s dir = (x1 - x0).normalized();
    
    Vector2s edge_center = (x0 + x1) * 0.5;
    
    Vector2s pressure_grad = fluid2d->get_pressure_gradient( edge_center );
    HairFlow<2>::m_u(i) -= pressure_grad.dot(dir) * dt / rho_L;
  }
  
  for(int i = 0; i < np; ++i)
  {
    int pidx = HairFlow<2>::m_particle_indices[i];
    const Vector2s& xp = x.segment<2>( m_parent->getDof( pidx ) );
    Matrix2s c = fluid2d->get_affine_matrix(xp);
    m_c_v.block<1, 2>(i, 0) = c.col(0).transpose();
    m_c_v.block<1, 2>(i, 2) = c.col(1).transpose();
  }
}

template<>
void CylindricalShallowFlow<3>::updateToFilteredGrid(const VectorXs& x, const VectorXs& v, FluidSim* fluidsim, const scalar& dt, int ibuffer)
{
  FluidSim3D* fluid3d = (FluidSim3D*) fluidsim;
  scalar cellsize = fluid3d->cellsize();
  
  const VectorXs& drag_buffer = m_parent->getFluidDragBuffer();
  std::vector<EdgeVelDragIntersection<3> >& u_edge_vel_drag = fluid3d->get_u_edge_vel_drag()[ibuffer];
  std::vector<EdgeVelDragIntersection<3> >& v_edge_vel_drag = fluid3d->get_v_edge_vel_drag()[ibuffer];
  std::vector<EdgeVelDragIntersection<3> >& w_edge_vel_drag = fluid3d->get_w_edge_vel_drag()[ibuffer];
  
  u_edge_vel_drag.resize(0);
  v_edge_vel_drag.resize(0);
  w_edge_vel_drag.resize(0);
  
  const scalar& rho = m_parent->getLiquidDensity();

  // for each edge
  int ne = HairFlow<3>::m_global_edges.size();
  for(int eidx = 0; eidx < ne; ++eidx) {
    auto& e = HairFlow<3>::m_global_edges[eidx];
    int local_first = m_global_to_local.find(e.first)->second;
    int local_second = m_global_to_local.find(e.second)->second;
    
    const scalar eta_first = m_eta(local_first);
    const scalar eta_second = m_eta(local_second);
    const scalar radii_first = m_rad_vec(local_first);
    const scalar radii_second = m_rad_vec(local_second);
    const scalar er_first = eta_first + radii_first;
    const scalar er_second = eta_second + eta_second;
    const scalar rho_h_first = m_parent->getHairDensity(e.first);
    const scalar rho_h_second = m_parent->getHairDensity(e.second);
    
    const Vector3s& u_first = m_actual_u_v.row(local_first).transpose();
    const Vector3s& u_second = m_actual_u_v.row(local_second).transpose();
    
    const Vector3s& x0 = x.segment<3>( m_parent->getDof( e.first ) );
    const Vector3s& x1 = x.segment<3>( m_parent->getDof( e.second ) );
    const Vector3s& v0 = v.segment<3>( m_parent->getDof( e.first ) );
    const Vector3s& v1 = v.segment<3>( m_parent->getDof( e.second ) );
    const Vector3s& d0 = drag_buffer.segment<3>( m_parent->getDof( e.first ) );
    const Vector3s& d1 = drag_buffer.segment<3>( m_parent->getDof( e.second ) );
    const Vector9s& c0 = m_c_v.block<1,9>(local_first,0).transpose();
    const Vector9s& c1 = m_c_v.block<1,9>(local_second,0).transpose();
    
    scalar xmin;
    scalar xmax;
    scalar ymin;
    scalar ymax;
    scalar zmin;
    scalar zmax;
    
    if(x0(0) - eta_first < x1(0) - eta_second) {
      xmin = x0(0) - std::max(cellsize * 2.0, er_first);
    } else {
      xmin = x1(0) - std::max(cellsize * 2.0, er_second);
    }
    
    if(x0(0) + eta_first > x1(0) + eta_second) {
      xmax = x0(0) + std::max(cellsize * 2.0, er_first);
    } else {
      xmax = x1(0) + std::max(cellsize * 2.0, er_second);
    }
    
    if(x0(1) - eta_first < x1(1) - eta_second) {
      ymin = x0(1) - std::max(cellsize * 2.0, er_first);
    } else {
      ymin = x1(1) - std::max(cellsize * 2.0, er_second);
    }
    
    if(x0(1) + eta_first > x1(1) + eta_second) {
      ymax = x0(1) + std::max(cellsize * 2.0, er_first);
    } else {
      ymax = x1(1) + std::max(cellsize * 2.0, er_second);
    }
    
    if(x0(2) - eta_first < x1(2) - eta_second) {
      zmin = x0(2) - std::max(cellsize * 2.0, er_first);
    } else {
      zmin = x1(2) - std::max(cellsize * 2.0, er_second);
    }
    
    if(x0(2) + eta_first > x1(2) + eta_second) {
      zmax = x0(2) + std::max(cellsize * 2.0, er_first);
    } else {
      zmax = x1(2) + std::max(cellsize * 2.0, er_second);
    }
    
    int xmin_u = std::max(0, std::min(fluid3d->get_u_ni()-1, (int) ceil((xmin - fluid3d->get_origin()(0)) / cellsize)));
    int xmin_v = std::max(0, std::min(fluid3d->get_v_ni()-1, (int) ceil((xmin - fluid3d->get_origin()(0)) / cellsize - 0.5)));
    int xmin_w = std::max(0, std::min(fluid3d->get_w_ni()-1, (int) ceil((xmin - fluid3d->get_origin()(0)) / cellsize - 0.5)));
    int xmax_u = std::max(0, std::min(fluid3d->get_u_ni()-1, (int) floor((xmax - fluid3d->get_origin()(0)) / cellsize)));
    int xmax_v = std::max(0, std::min(fluid3d->get_v_ni()-1, (int) floor((xmax - fluid3d->get_origin()(0)) / cellsize - 0.5)));
    int xmax_w = std::max(0, std::min(fluid3d->get_w_ni()-1, (int) floor((xmax - fluid3d->get_origin()(0)) / cellsize - 0.5)));
    
    int ymin_u = std::max(0, std::min(fluid3d->get_u_nj()-1, (int) ceil((ymin - fluid3d->get_origin()(1)) / cellsize - 0.5)));
    int ymin_v = std::max(0, std::min(fluid3d->get_v_nj()-1, (int) ceil((ymin - fluid3d->get_origin()(1)) / cellsize)));
    int ymin_w = std::max(0, std::min(fluid3d->get_w_nj()-1, (int) ceil((ymin - fluid3d->get_origin()(1)) / cellsize - 0.5)));
    int ymax_u = std::max(0, std::min(fluid3d->get_u_nj()-1, (int) floor((ymax - fluid3d->get_origin()(1)) / cellsize - 0.5)));
    int ymax_v = std::max(0, std::min(fluid3d->get_v_nj()-1, (int) floor((ymax - fluid3d->get_origin()(1)) / cellsize)));
    int ymax_w = std::max(0, std::min(fluid3d->get_w_nj()-1, (int) floor((ymax - fluid3d->get_origin()(1)) / cellsize - 0.5)));
    
    int zmin_u = std::max(0, std::min(fluid3d->get_u_nk()-1, (int) ceil((zmin - fluid3d->get_origin()(2)) / cellsize - 0.5)));
    int zmin_v = std::max(0, std::min(fluid3d->get_v_nk()-1, (int) ceil((zmin - fluid3d->get_origin()(2)) / cellsize - 0.5)));
    int zmin_w = std::max(0, std::min(fluid3d->get_w_nk()-1, (int) ceil((zmin - fluid3d->get_origin()(2)) / cellsize)));
    int zmax_u = std::max(0, std::min(fluid3d->get_u_nk()-1, (int) floor((zmax - fluid3d->get_origin()(2)) / cellsize - 0.5)));
    int zmax_v = std::max(0, std::min(fluid3d->get_v_nk()-1, (int) floor((zmax - fluid3d->get_origin()(2)) / cellsize - 0.5)));
    int zmax_w = std::max(0, std::min(fluid3d->get_w_nk()-1, (int) floor((zmax - fluid3d->get_origin()(2)) / cellsize)));
    
    std::vector<Vector3i> line_buffer;
    for(int k = zmin_u; k <= zmax_u; ++k) for(int j = ymin_u; j <= ymax_u; ++j) for(int i = xmin_u; i <= xmax_u; ++i)
    {
      Vector3s pos = Vector3s(i*cellsize, (j+0.5)*cellsize, (k+0.5)*cellsize) + fluid3d->get_origin();
      Vector6s clipping_window;
      clipping_window << (pos(0) - cellsize), (pos(1) - cellsize), (pos(2) - cellsize),
      (pos(0) + cellsize), (pos(1) + cellsize), (pos(2) + cellsize);
      
      Vector3s q0 = x0, q1 = x1;
      scalar alpha0 = 0.0, alpha1 = 1.0;
      if(!liangbarsky::clip_line(clipping_window, q0, q1, alpha0, alpha1)) continue;
      
      scalar length_segment = (q1 - q0).norm();
      
      for(int r = 0; r < N_GAUSS; ++r)
      {
        const scalar beta = mathutils::gauss_legendre_point[N_GAUSS-1][r];
        const scalar w = mathutils::gauss_legendre_weight[N_GAUSS-1][r];
        
        scalar alpha = mathutils::lerp(alpha0, alpha1, beta);
        Vector3s xc = x0 * (1.0 - alpha) + x1 * alpha;
        Vector3s vel_c = v0 * (1.0 - alpha) + v1 * alpha;
        Vector3s xcb = xc - vel_c * dt;
        Vector3i ixcb = Vector3i( (int)((xcb(0) - fluid3d->get_origin()(0)) / cellsize), (int)((xcb(1) - fluid3d->get_origin()(1)) / cellsize - 0.5), (int)((xcb(2) - fluid3d->get_origin()(2)) / cellsize - 0.5 ));
        
        line_buffer.resize(0);
        bresemham::Bresenham3D(ixcb(0), ixcb(1), ixcb(2), i, j, k, [&] (int ii, int jj, int kk) {
          line_buffer.push_back(Vector3i(ii, jj, kk));
        });
        if(line_buffer.size() == 0) line_buffer.push_back(Vector3i(i, j, k));
        
        scalar ubase_c = vel_c(0);
        scalar uflow_c = u_first(0) * (1.0 - alpha) + u_second(0) * alpha;
        scalar vol_c = M_PI * mathutils::lerp(er_first * er_first - radii_first * radii_first, er_second * er_second - radii_second * radii_second, alpha) * length_segment;
        scalar mass_c = rho * vol_c;
        scalar rho_base_c = mathutils::lerp(rho_h_first, rho_h_second, alpha);
        scalar vol_base_c = M_PI * mathutils::lerp(radii_first * radii_first, radii_second * radii_second, alpha) * length_segment;
        scalar mass_base_c = rho_base_c * vol_base_c;
        scalar dragc = d0(0) * (1.0 - alpha) + d1(0) * alpha;
        Vector3s c_c = c0.segment<3>(0) * (1.0 - alpha) + c1.segment<3>(0) * alpha;
        
        Vector3s diff = pos - xc;
        
        scalar w_linear = mathutils::linear_kernel(diff, cellsize);
        scalar weight = (mass_c + mass_base_c) * w_linear;
        
        scalar u_base_world = ubase_c;
        scalar u_flow_world = u_base_world + uflow_c + c_c.dot(pos - xc);
        
        for(const Vector3i& cd : line_buffer) {
          EdgeVelDragIntersection<3> tmp_inter;
          tmp_inter.vel_weighted = (u_base_world * mass_base_c + u_flow_world * mass_c) * w_linear * w;
          tmp_inter.drag_weighted = dragc * weight * w / (scalar) line_buffer.size();
          tmp_inter.vol_weighted = (vol_c + vol_base_c) * w_linear * w;
          tmp_inter.weight = weight * w;
          tmp_inter.linear_weight = w_linear * w;
          tmp_inter.coord = cd;
          u_edge_vel_drag.push_back(tmp_inter);
        }
      }
    }
 
    for(int k = zmin_v; k <= zmax_v; ++k) for(int j = ymin_v; j <= ymax_v; ++j) for(int i = xmin_v; i <= xmax_v; ++i)
    {
      Vector3s pos = Vector3s((i+0.5)*cellsize,j*cellsize,(k+0.5)*cellsize) + fluid3d->get_origin();
      Vector6s clipping_window;
      clipping_window << (pos(0) - cellsize), (pos(1) - cellsize), (pos(2) - cellsize),
      (pos(0) + cellsize), (pos(1) + cellsize), (pos(2) + cellsize);
      
      Vector3s q0 = x0, q1 = x1;
      scalar alpha0 = 0.0, alpha1 = 1.0;
      if(!liangbarsky::clip_line(clipping_window, q0, q1, alpha0, alpha1)) continue;
      
      scalar length_segment = (q1 - q0).norm();
      
      for(int r = 0; r < N_GAUSS; ++r)
      {
        const scalar beta = mathutils::gauss_legendre_point[N_GAUSS-1][r];
        const scalar w = mathutils::gauss_legendre_weight[N_GAUSS-1][r];
        
        scalar alpha = mathutils::lerp(alpha0, alpha1, beta);
        Vector3s xc = x0 * (1.0 - alpha) + x1 * alpha;
        Vector3s vel_c = v0 * (1.0 - alpha) + v1 * alpha;
        Vector3s xcb = xc - vel_c * dt;
        Vector3i ixcb = Vector3i( (int)((xcb(0) - fluid3d->get_origin()(0)) / cellsize - 0.5), (int)((xcb(1) - fluid3d->get_origin()(1)) / cellsize), (int)((xcb(2) - fluid3d->get_origin()(2)) / cellsize - 0.5 ));
        
        line_buffer.resize(0);
        bresemham::Bresenham3D(ixcb(0), ixcb(1), ixcb(2), i, j, k, [&] (int ii, int jj, int kk) {
          line_buffer.push_back(Vector3i(ii, jj, kk));
        });
        if(line_buffer.size() == 0) line_buffer.push_back(Vector3i(i, j, k));
        
        scalar vbase_c = vel_c(1);
        scalar vflow_c = u_first(1) * (1.0 - alpha) + u_second(1) * alpha;
        scalar vol_c = M_PI * mathutils::lerp(er_first * er_first - radii_first * radii_first, er_second * er_second - radii_second * radii_second, alpha) * length_segment;
        scalar mass_c = rho * vol_c;
        scalar rho_base_c = mathutils::lerp(rho_h_first, rho_h_second, alpha);
        scalar vol_base_c = M_PI * mathutils::lerp(radii_first * radii_first, radii_second * radii_second, alpha) * length_segment;
        scalar mass_base_c = rho_base_c * vol_base_c;
        scalar dragc = d0(1) * (1.0 - alpha) + d1(1) * alpha;
        Vector3s c_c = c0.segment<3>(3) * (1.0 - alpha) + c1.segment<3>(3) * alpha;
        
        Vector3s diff = pos - xc;
        
        scalar w_linear = mathutils::linear_kernel(diff, cellsize);
        scalar weight = (mass_c + mass_base_c) * w_linear;
        
        scalar v_base_world = vbase_c;
        scalar v_flow_world = v_base_world + vflow_c + c_c.dot(pos - xc);
        
        for(const Vector3i& cd : line_buffer) {
          EdgeVelDragIntersection<3> tmp_inter;
          tmp_inter.vel_weighted = (v_base_world * mass_base_c + v_flow_world * mass_c) * w_linear * w;
          tmp_inter.drag_weighted = dragc * weight * w / (scalar) line_buffer.size();
          tmp_inter.vol_weighted = (vol_c + vol_base_c) * w_linear * w;
          tmp_inter.weight = weight * w;
          tmp_inter.linear_weight = w_linear * w;
          tmp_inter.coord = cd;
          
          v_edge_vel_drag.push_back(tmp_inter);
        }
      }
    }
    
    for(int k = zmin_w; k <= zmax_w; ++k) for(int j = ymin_w; j <= ymax_w; ++j) for(int i = xmin_w; i <= xmax_w; ++i)
    {
      Vector3s pos = Vector3s((i+0.5)*cellsize,(j+0.5)*cellsize,k*cellsize) + fluid3d->get_origin();
      Vector6s clipping_window;
      clipping_window << (pos(0) - cellsize), (pos(1) - cellsize), (pos(2) - cellsize),
      (pos(0) + cellsize), (pos(1) + cellsize), (pos(2) + cellsize);
      
      Vector3s q0 = x0, q1 = x1;
      scalar alpha0 = 0.0, alpha1 = 1.0;
      if(!liangbarsky::clip_line(clipping_window, q0, q1, alpha0, alpha1)) continue;
      
      scalar length_segment = (q1 - q0).norm();
      
      for(int r = 0; r < N_GAUSS; ++r)
      {
        const scalar beta = mathutils::gauss_legendre_point[N_GAUSS-1][r];
        const scalar w = mathutils::gauss_legendre_weight[N_GAUSS-1][r];
        
        scalar alpha = mathutils::lerp(alpha0, alpha1, beta);
        Vector3s xc = x0 * (1.0 - alpha) + x1 * alpha;
        Vector3s vel_c = v0 * (1.0 - alpha) + v1 * alpha;
        Vector3s xcb = xc - vel_c * dt;
        Vector3i ixcb = Vector3i( (int)((xcb(0) - fluid3d->get_origin()(0)) / cellsize - 0.5), (int)((xcb(1) - fluid3d->get_origin()(1)) / cellsize - 0.5), (int)((xcb(2) - fluid3d->get_origin()(2)) / cellsize));
        
        line_buffer.resize(0);
        bresemham::Bresenham3D(ixcb(0), ixcb(1), ixcb(2), i, j, k, [&] (int ii, int jj, int kk) {
          line_buffer.push_back(Vector3i(ii, jj, kk));
        });
        if(line_buffer.size() == 0) line_buffer.push_back(Vector3i(i, j, k));

        scalar wbase_c = vel_c(2);
        scalar wflow_c = u_first(2) * (1.0 - alpha) + u_second(2) * alpha;
        scalar vol_c = M_PI * mathutils::lerp(er_first * er_first - radii_first * radii_first, er_second * er_second - radii_second * radii_second, alpha) * length_segment;
        scalar mass_c = rho * vol_c;
        scalar rho_base_c = mathutils::lerp(rho_h_first, rho_h_second, alpha);
        scalar vol_base_c = M_PI * mathutils::lerp(radii_first * radii_first, radii_second * radii_second, alpha) * length_segment;
        scalar mass_base_c = rho_base_c * vol_base_c;
        scalar dragc = d0(2) * (1.0 - alpha) + d1(2) * alpha;
        Vector3s c_c = c0.segment<3>(6) * (1.0 - alpha) + c1.segment<3>(6) * alpha;
        
        Vector3s diff = pos - xc;
        
        scalar w_linear = mathutils::linear_kernel(diff, cellsize);
        scalar weight = (mass_c + mass_base_c) * w_linear;
        
        scalar w_base_world = wbase_c;
        scalar w_flow_world = w_base_world + wflow_c + c_c.dot(pos - xc);
        
        for(const Vector3i& cd : line_buffer) {
          EdgeVelDragIntersection<3> tmp_inter;
          tmp_inter.vel_weighted = (w_base_world * mass_base_c + w_flow_world * mass_c) * w_linear * w;
          tmp_inter.drag_weighted = dragc * weight * w / (scalar) line_buffer.size();
          tmp_inter.vol_weighted = (vol_c + vol_base_c) * w_linear * w;
          tmp_inter.weight = weight * w;
          tmp_inter.linear_weight = w_linear * w;
          tmp_inter.coord = cd;
          
          w_edge_vel_drag.push_back(tmp_inter);
        }
      }
    }
  }
}

template<>
void CylindricalShallowFlow<3>::updateFromFilteredGrid(const VectorXs& x, VectorXs& v, FluidSim* fluidsim, const scalar& dt)
{
  FluidSim3D* fluid3d = (FluidSim3D*) fluidsim;
  
  // for each edge
  int ne = HairFlow<3>::m_global_edges.size();
  int np = HairFlow<3>::m_particle_indices.size();
  
  const scalar rho_L = m_parent->getLiquidDensity();
  bool visc_solve = m_parent->viscositySolve();
  
  for(int i = 0; i < ne; ++i)
  {
    auto& e = HairFlow<3>::m_global_edges[i];
    
    const Vector3s& x0 = x.segment<3>( m_parent->getDof( e.first ) );
    const Vector3s& x1 = x.segment<3>( m_parent->getDof( e.second ) );
    
    Vector3s dir = (x1 - x0).normalized();
    
    Vector3s edge_center = (x0 + x1) * 0.5;

    Vector3s pressure_grad = fluid3d->get_pressure_gradient( edge_center );
    
    HairFlow<3>::m_u(i) -= pressure_grad.dot(dir) * dt / rho_L;
    
    if(visc_solve) {
      HairFlow<3>::m_u(i) += fluid3d->get_visc_impulse( edge_center ).dot(dir) * dt;
    }
  }
  
  for(int i = 0; i < np; ++i)
  {
    int pidx = HairFlow<3>::m_particle_indices[i];
    const Vector3s& xp = x.segment<3>( m_parent->getDof( pidx ) );
    Matrix3s c = fluid3d->get_affine_matrix(xp);
    m_c_v.block<1, 3>(i, 0) = c.col(0).transpose();
    m_c_v.block<1, 3>(i, 3) = c.col(1).transpose();
    m_c_v.block<1, 3>(i, 6) = c.col(2).transpose();
  }
}

template<>
Vector3s CylindricalShallowFlow<3>::computeHairLiquidAngularMomentum(const VectorXs& x, const VectorXs& v, FluidSim* fluidsim) const
{
  FluidSim3D* fluid3d = (FluidSim3D*) fluidsim;
  scalar cellsize = fluid3d->cellsize();
  
  VectorXs gradE;
  gradE.resize(x.size());
  
  Vector3s sum = Vector3s::Zero();
  
  const scalar& rho = m_parent->getLiquidDensity();
  // for each edge
  for(auto& e : HairFlow<3>::m_global_edges)
  {
    int local_first = m_global_to_local.find(e.first)->second;
    int local_second = m_global_to_local.find(e.second)->second;
    
    const scalar eta_first = m_eta(local_first);
    const scalar eta_second = m_eta(local_second);
    const scalar radii_first = m_rad_vec(local_first);
    const scalar radii_second = m_rad_vec(local_second);
    const scalar er_first = eta_first + radii_first;
    const scalar er_second = eta_second + eta_second;
    const scalar rho_h_first = m_parent->getHairDensity(e.first);
    const scalar rho_h_second = m_parent->getHairDensity(e.second);
    
    const scalar u_first = m_uvert(local_first);
    const scalar u_second = m_uvert(local_second);
    
    const Vector3s& x0 = x.segment<3>( m_parent->getDof( e.first ) );
    const Vector3s& x1 = x.segment<3>( m_parent->getDof( e.second ) );
    const Vector3s& v0 = v.segment<3>( m_parent->getDof( e.first ) );
    const Vector3s& v1 = v.segment<3>( m_parent->getDof( e.second ) );
    const Vector9s& c0 = m_c_v.block<1,9>(local_first,0).transpose();
    const Vector9s& c1 = m_c_v.block<1,9>(local_second,0).transpose();
    
    Vector3s dir = (x1 - x0).normalized();
    
    scalar xmin;
    scalar xmax;
    scalar ymin;
    scalar ymax;
    scalar zmin;
    scalar zmax;
    
    if(x0(0) - eta_first < x1(0) - eta_second) {
      xmin = x0(0) - std::max(cellsize * 2.0, eta_first);
    } else {
      xmin = x1(0) - std::max(cellsize * 2.0, eta_second);
    }
    
    if(x0(0) + eta_first > x1(0) + eta_second) {
      xmax = x0(0) + std::max(cellsize * 2.0, eta_first);
    } else {
      xmax = x1(0) + std::max(cellsize * 2.0, eta_second);
    }
    
    if(x0(1) - eta_first < x1(1) - eta_second) {
      ymin = x0(1) - std::max(cellsize * 2.0, eta_first);
    } else {
      ymin = x1(1) - std::max(cellsize * 2.0, eta_second);
    }
    
    if(x0(1) + eta_first > x1(1) + eta_second) {
      ymax = x0(1) + std::max(cellsize * 2.0, eta_first);
    } else {
      ymax = x1(1) + std::max(cellsize * 2.0, eta_second);
    }
    
    if(x0(2) - eta_first < x1(2) - eta_second) {
      zmin = x0(2) - std::max(cellsize * 2.0, eta_first);
    } else {
      zmin = x1(2) - std::max(cellsize * 2.0, eta_second);
    }
    
    if(x0(2) + eta_first > x1(2) + eta_second) {
      zmax = x0(2) + std::max(cellsize * 2.0, eta_first);
    } else {
      zmax = x1(2) + std::max(cellsize * 2.0, eta_second);
    }
    
    int xmin_u = std::max(0, std::min(fluid3d->get_u_ni()-1, (int) ceil((xmin - fluid3d->get_origin()(0)) / cellsize)));
    int xmin_v = std::max(0, std::min(fluid3d->get_v_ni()-1, (int) ceil((xmin - fluid3d->get_origin()(0)) / cellsize - 0.5)));
    int xmin_w = std::max(0, std::min(fluid3d->get_w_ni()-1, (int) ceil((xmin - fluid3d->get_origin()(0)) / cellsize - 0.5)));
    int xmax_u = std::max(0, std::min(fluid3d->get_u_ni()-1, (int) floor((xmax - fluid3d->get_origin()(0)) / cellsize)));
    int xmax_v = std::max(0, std::min(fluid3d->get_v_ni()-1, (int) floor((xmax - fluid3d->get_origin()(0)) / cellsize - 0.5)));
    int xmax_w = std::max(0, std::min(fluid3d->get_w_ni()-1, (int) floor((xmax - fluid3d->get_origin()(0)) / cellsize - 0.5)));
    
    int ymin_u = std::max(0, std::min(fluid3d->get_u_nj()-1, (int) ceil((ymin - fluid3d->get_origin()(1)) / cellsize - 0.5)));
    int ymin_v = std::max(0, std::min(fluid3d->get_v_nj()-1, (int) ceil((ymin - fluid3d->get_origin()(1)) / cellsize)));
    int ymin_w = std::max(0, std::min(fluid3d->get_w_nj()-1, (int) ceil((ymin - fluid3d->get_origin()(1)) / cellsize - 0.5)));
    int ymax_u = std::max(0, std::min(fluid3d->get_u_nj()-1, (int) floor((ymax - fluid3d->get_origin()(1)) / cellsize - 0.5)));
    int ymax_v = std::max(0, std::min(fluid3d->get_v_nj()-1, (int) floor((ymax - fluid3d->get_origin()(1)) / cellsize)));
    int ymax_w = std::max(0, std::min(fluid3d->get_w_nj()-1, (int) floor((ymax - fluid3d->get_origin()(1)) / cellsize - 0.5)));
    
    int zmin_u = std::max(0, std::min(fluid3d->get_u_nk()-1, (int) ceil((zmin - fluid3d->get_origin()(2)) / cellsize - 0.5)));
    int zmin_v = std::max(0, std::min(fluid3d->get_v_nk()-1, (int) ceil((zmin - fluid3d->get_origin()(2)) / cellsize - 0.5)));
    int zmin_w = std::max(0, std::min(fluid3d->get_w_nk()-1, (int) ceil((zmin - fluid3d->get_origin()(2)) / cellsize)));
    int zmax_u = std::max(0, std::min(fluid3d->get_u_nk()-1, (int) floor((zmax - fluid3d->get_origin()(2)) / cellsize - 0.5)));
    int zmax_v = std::max(0, std::min(fluid3d->get_v_nk()-1, (int) floor((zmax - fluid3d->get_origin()(2)) / cellsize - 0.5)));
    int zmax_w = std::max(0, std::min(fluid3d->get_w_nk()-1, (int) floor((zmax - fluid3d->get_origin()(2)) / cellsize)));
    
    for(int k = zmin_u; k <= zmax_u; ++k) for(int j = ymin_u; j <= ymax_u; ++j) for(int i = xmin_u; i <= xmax_u; ++i)
    {
      Vector3s pos = Vector3s(i*cellsize, (j+0.5)*cellsize, (k+0.5)*cellsize) + fluid3d->get_origin();
      Vector6s clipping_window;
      clipping_window << (pos(0) - cellsize), (pos(1) - cellsize), (pos(2) - cellsize),
      (pos(0) + cellsize), (pos(1) + cellsize), (pos(2) + cellsize);
      
      Vector3s q0 = x0, q1 = x1;
      scalar alpha0, alpha1;
      if(!liangbarsky::clip_line(clipping_window, q0, q1, alpha0, alpha1)) continue;
      
      scalar length_segment = (q1 - q0).norm();
      
      for(int r = 0; r < N_GAUSS; ++r)
      {
        const scalar beta = mathutils::gauss_legendre_point[N_GAUSS-1][r];
        const scalar w = mathutils::gauss_legendre_weight[N_GAUSS-1][r];
        
        scalar alpha = mathutils::lerp(alpha0, alpha1, beta);
        Vector3s xc = x0 * (1.0 - alpha) + x1 * alpha;
        scalar ubase_c = v0(0) * (1.0 - alpha) + v1(0) * alpha;
        scalar uflow_c = u_first * (1.0 - alpha) + u_second * alpha;
        scalar mass_c = rho * M_PI * mathutils::lerp(er_first * er_first - radii_first * radii_first, er_second * er_second - radii_second * radii_second, alpha) * length_segment;
        scalar rho_base_c = mathutils::lerp(rho_h_first, rho_h_second, alpha);
        scalar mass_base_c = rho_base_c * M_PI * mathutils::lerp(radii_first * radii_first, radii_second * radii_second, alpha) * length_segment;
        Vector3s c_c = c0.segment<3>(0) * (1.0 - alpha) + c1.segment<3>(0) * alpha;
        
        Vector3s diff = pos - xc;
        
        scalar w_linear = mathutils::linear_kernel(diff, cellsize);
        
        scalar u_base_world = ubase_c;
        scalar u_flow_world = u_base_world + uflow_c + c_c.dot(pos - xc);
        
        Vector3s pu = Vector3s((u_base_world * mass_base_c + u_flow_world * mass_c) * w_linear * w, 0.0, 0.0);
        sum += pos.cross(pu);
      }
    }
    
    for(int k = zmin_v; k <= zmax_v; ++k) for(int j = ymin_v; j <= ymax_v; ++j) for(int i = xmin_v; i <= xmax_v; ++i)
    {
      Vector3s pos = Vector3s((i+0.5)*cellsize,j*cellsize,(k+0.5)*cellsize) + fluid3d->get_origin();
      Vector6s clipping_window;
      clipping_window << (pos(0) - cellsize), (pos(1) - cellsize), (pos(2) - cellsize),
      (pos(0) + cellsize), (pos(1) + cellsize), (pos(2) + cellsize);
      Vector3s q0 = x0, q1 = x1;
      scalar alpha0, alpha1;
      if(!liangbarsky::clip_line(clipping_window, q0, q1, alpha0, alpha1)) continue;
      
      scalar length_segment = (q1 - q0).norm();
      
      for(int r = 0; r < N_GAUSS; ++r)
      {
        const scalar beta = mathutils::gauss_legendre_point[N_GAUSS-1][r];
        const scalar w = mathutils::gauss_legendre_weight[N_GAUSS-1][r];
        
        scalar alpha = mathutils::lerp(alpha0, alpha1, beta);
        Vector3s xc = x0 * (1.0 - alpha) + x1 * alpha;
        scalar vbase_c = v0(1) * (1.0 - alpha) + v1(1) * alpha;
        scalar vflow_c = u_first * (1.0 - alpha) + u_second * alpha;
        
        scalar mass_c = rho * M_PI * mathutils::lerp(er_first * er_first - radii_first * radii_first, er_second * er_second - radii_second * radii_second, alpha) * length_segment;
        scalar rho_base_c = mathutils::lerp(rho_h_first, rho_h_second, alpha);
        scalar mass_base_c = rho_base_c * M_PI * mathutils::lerp(radii_first * radii_first, radii_second * radii_second, alpha) * length_segment;
        Vector3s c_c = c0.segment<3>(3) * (1.0 - alpha) + c1.segment<3>(3) * alpha;
        
        Vector3s diff = pos - xc;
        
        scalar w_linear = mathutils::linear_kernel(diff, cellsize);
        scalar v_base_world = vbase_c;
        scalar v_flow_world = v_base_world + vflow_c + c_c.dot(pos - xc);
        
        Vector3s pv = Vector3s(0.0, (v_base_world * mass_base_c + v_flow_world * mass_c) * w_linear * w, 0.0);
        sum += pos.cross(pv);
      }
    }
    
    for(int k = zmin_w; k <= zmax_w; ++k) for(int j = ymin_w; j <= ymax_w; ++j) for(int i = xmin_w; i <= xmax_w; ++i)
    {
      Vector3s pos = Vector3s(i*cellsize,(j+0.5)*cellsize,(k+0.5)*cellsize) + fluid3d->get_origin();
      Vector6s clipping_window;
      clipping_window << (pos(0) - cellsize), (pos(1) - cellsize), (pos(2) - cellsize),
      (pos(0) + cellsize), (pos(1) + cellsize), (pos(2) + cellsize);
      Vector3s q0 = x0, q1 = x1;
      scalar alpha0, alpha1;
      if(!liangbarsky::clip_line(clipping_window, q0, q1, alpha0, alpha1)) continue;
      
      scalar length_segment = (q1 - q0).norm();
      
      for(int r = 0; r < N_GAUSS; ++r)
      {
        const scalar beta = mathutils::gauss_legendre_point[N_GAUSS-1][k];
        const scalar w = mathutils::gauss_legendre_weight[N_GAUSS-1][k];
        
        scalar alpha = mathutils::lerp(alpha0, alpha1, beta);
        Vector3s xc = x0 * (1.0 - alpha) + x1 * alpha;
        scalar wbase_c = v0(2) * (1.0 - alpha) + v1(2) * alpha;
        scalar wflow_c = u_first * (1.0 - alpha) + u_second * alpha;
        scalar mass_c = rho * M_PI * mathutils::lerp(er_first * er_first - radii_first * radii_first, er_second * er_second - radii_second * radii_second, alpha) * length_segment;
        scalar rho_base_c = mathutils::lerp(rho_h_first, rho_h_second, alpha);
        scalar mass_base_c = rho_base_c * M_PI * mathutils::lerp(radii_first * radii_first, radii_second * radii_second, alpha) * length_segment;
        Vector3s c_c = c0.segment<3>(6) * (1.0 - alpha) + c1.segment<3>(6) * alpha;
        
        Vector3s diff = pos - xc;
        
        scalar w_linear = mathutils::linear_kernel(diff, cellsize);
        scalar w_base_world = wbase_c;
        scalar w_flow_world = w_base_world + wflow_c + c_c.dot(pos - xc);
        
        Vector3s pv = Vector3s(0.0, 0.0, (w_base_world * mass_base_c + w_flow_world * mass_c) * w_linear * w);
        sum += pos.cross(pv);
      }
    }
  }
  
  return sum;
}

template<int DIM>
void CylindricalShallowFlow<DIM>::advance(const VectorXs& x, const scalar& dt)
{
  // update height field due to manifold change
  //compute_stretched_liquid_height(m_old_rad_vec, HairFlow<DIM>::m_rad_vec, m_old_area_v, HairFlow<DIM>::m_area_v, HairFlow<DIM>::m_eta); stringutils::print(HairFlow<DIM>::m_eta);
  
  // advection
  compute_vertex_val(HairFlow<DIM>::m_u, m_W_fv_hair, m_uvert);
  
  compute_backtracing<DIM>(m_edge_x, dt, m_sum_area_e, m_area_e_accu, m_area_e_inv_mapping,
                           m_uvert, HairFlow<DIM>::m_c_v, m_ustar, HairFlow<DIM>::m_c_star);
  
  compute_vertex_val(HairFlow<DIM>::m_c_star, m_W_fv_hair, HairFlow<DIM>::m_c_v);
}

template<int DIM>
void CylindricalShallowFlow<DIM>::geodesic_to_local(const scalar& geopos, int& pidx_low, scalar& alpha) const
{
  scalar pos_hair = find_pos_on_hair(geopos, m_sum_area_e, m_area_e_accu, m_area_e_inv_mapping);
  pidx_low = (int) pos_hair;
  alpha = pos_hair - floor(pos_hair);
}

template<int DIM>
scalar CylindricalShallowFlow<DIM>::getTotalLength() const
{
  return m_sum_area_e;
}

template<int DIM>
scalar CylindricalShallowFlow<DIM>::computeTotalLiquidVol() const
{
  scalar sum = 0.0;
  const int ne = HairFlow<DIM>::m_area_e.size();
  
  for(int i = 0; i < ne; ++i)
  {
    int p0 = i;
    int p1 = i + 1;
    scalar H0 = HairFlow<DIM>::m_eta(p0) + HairFlow<DIM>::m_rad_vec(p0);
    scalar H1 = HairFlow<DIM>::m_eta(p1) + HairFlow<DIM>::m_rad_vec(p1);
    
    sum += (H1 * H1 - HairFlow<DIM>::m_rad_vec(p1) * HairFlow<DIM>::m_rad_vec(p1) + H0 * H0 - HairFlow<DIM>::m_rad_vec(p0) * HairFlow<DIM>::m_rad_vec(p0)) * M_PI * HairFlow<DIM>::m_area_e(i) * 0.5;
  }
  return sum;
}

template<int DIM>
scalar CylindricalShallowFlow<DIM>::computeTotalReservoirVol() const
{
  return m_pool_liquid;
}

template<int DIM>
const scalar& CylindricalShallowFlow<DIM>::getPoolSize() const
{
  return m_pool_liquid;
}

template<int DIM>
scalar& CylindricalShallowFlow<DIM>::getPoolSize()
{
  return m_pool_liquid;
}

template<int DIM>
void CylindricalShallowFlow<DIM>::add_force(const VectorXs& x, const VectorXs& accel, FluidSim* fluidsim, const scalar& dt)
{
  const int np = HairFlow<DIM>::m_particle_indices.size();
  const PolygonalCohesion<DIM>* cohesion = HairFlow<DIM>::m_parent->getPolygonalCohesion();
  
  const scalar& cap_multiplier = HairFlow<DIM>::m_parent->getCapillaryAccelMultiplier();
  const scalar& gamma = HairFlow<DIM>::m_parent->getLiquidTension();
  
  if(HairFlow<DIM>::m_parent->isIndividualTransfer())
  {
    compute_next_pressure(m_L, HairFlow<DIM>::m_eta, HairFlow<DIM>::m_parent->getLiquidTension(), m_pressure);
  } else {
    const VectorXs& global_pressure = cohesion->getPressureVGlobal();
    for(int i = 0; i < np; ++i)
      m_pressure(i) = global_pressure(HairFlow<DIM>::m_particle_indices[i]);
  }
  

  auto& gv = HairFlow<DIM>::m_parent->getSimpleGravity();
  VectorXs g(DIM);
  for(int i = 0; i < DIM; ++i) g(i) = gv(i);
  // add force
  compute_accel_v(accel, HairFlow<DIM>::m_particle_indices, g, HairFlow<DIM>::m_accel_v, HairFlow<DIM>::m_parent);
  
  
  compute_edge_val(HairFlow<DIM>::m_accel_v, HairFlow<DIM>::m_internal_edges, m_accel_f);
  compute_accel_e(m_accel_f, HairFlow<DIM>::m_dir_f, m_accel_e);
  if(cap_multiplier > 0.0) {
    compute_capillary_accel_e(HairFlow<DIM>::m_rad_vec, HairFlow<DIM>::m_eta, m_area_e_accu, gamma, m_sum_area_e, cap_multiplier, m_cap_e);
    
    m_accel_e += m_cap_e;
  }
  
  compute_next_velocity(m_ustar,
                        m_gradF,
                        m_dir_f_expand,
                        m_pressure,
                        m_porosity_e,
                        HairFlow<DIM>::m_edge_eta,
                        HairFlow<DIM>::m_edge_rad_vec,
                        dt, m_accel_e,
                        HairFlow<DIM>::m_parent->getLiquidDensity(),
                        HairFlow<DIM>::m_parent->getViscosity(),
                        HairFlow<DIM>::m_parent->getHairFrictionCoeff(),
                        HairFlow<DIM>::m_u);
  
  
  
  if(HairFlow<DIM>::m_parent->drippingNear())
    HairFlow<DIM>::m_u(0) = std::min(0.0, HairFlow<DIM>::m_u(0));
  else
    HairFlow<DIM>::m_u(0) = 0.0;
  
  if(HairFlow<DIM>::m_parent->drippingFar())
    HairFlow<DIM>::m_u(HairFlow<DIM>::m_u.size() - 1) = std::max(0.0, HairFlow<DIM>::m_u(HairFlow<DIM>::m_u.size() - 1));
  else
    HairFlow<DIM>::m_u(HairFlow<DIM>::m_u.size() - 1) = 0.0;
  
  for(int i = 0; i < HairFlow<DIM>::m_u.size(); ++i) {
    if(HairFlow<DIM>::m_edge_eta(i) - HairFlow<DIM>::m_edge_rad_vec(i) < 1e-6) HairFlow<DIM>::m_u(i) = 0.0;
  }
  
  compute_actual_u(HairFlow<DIM>::m_u, HairFlow<DIM>::m_dir_f, HairFlow<DIM>::m_actual_u_f);
  compute_vertex_val(HairFlow<DIM>::m_actual_u_f, m_W_fv_hair, HairFlow<DIM>::m_actual_u_v);
}

template<int DIM>
void CylindricalShallowFlow<DIM>::preUpdateHairFlowHeight(const scalar& dt)
{
  m_old_eta = HairFlow<DIM>::m_eta;

  compute_D_matrix(m_W_fv, HairFlow<DIM>::m_u, m_gradF, m_dir_f_expand, m_D);
  compute_A_matrix(m_D, m_divV, HairFlow<DIM>::m_u, m_dir_f_expand, m_L, HairFlow<DIM>::m_parent->getHeightSmooth(), m_A);
  
  compute_LHS_matrix(m_sp0, m_sp1, m_A, dt, m_LHS);
  
  const int np = HairFlow<DIM>::m_particle_indices.size();
  const VectorXs& v = HairFlow<DIM>::m_parent->getV();
  const scalar& liquid_shell = HairFlow<DIM>::m_parent->getLiquidShell();
  
  // compute rhs vector
  for(int i = 0; i < np; ++i)
  {
    int pidx = HairFlow<DIM>::m_particle_indices[i];
    const scalar& r = HairFlow<DIM>::m_rad_vec(i);
    const scalar& eta = HairFlow<DIM>::m_eta(i);
    const int idof = HairFlow<DIM>::m_parent->getDof(pidx);
    
    scalar A = (eta + r) * (eta + r) - r * r;
    m_rhs(i, 0) = A;
    
    scalar H = std::max(eta + r, liquid_shell * r);
    scalar modified_A = H * H - r * r;
    m_rhs(i, 1) = v(idof + 0) * modified_A;
    m_rhs(i, 2) = v(idof + 1) * modified_A;
    m_rhs(i, 3) = v(idof + 2) * modified_A;
    if(!HairFlow<DIM>::m_parent->isMassSpring()) {
      if(i == np - 1) {
        m_rhs(i, 4) = 0.0;
      } else {
        m_rhs(i, 4) = v(idof + 3) * modified_A * (modified_A / 2.0 + r * r);
      }
    }
  }
  
  m_LHS.makeCompressed();
  
  m_solver.compute(m_LHS);
}

template<int DIM>
void CylindricalShallowFlow<DIM>::postUpdateHairFlowHeight(const scalar& dt)
{
  compute_edge_val(HairFlow<DIM>::m_eta, HairFlow<DIM>::m_internal_edges, HairFlow<DIM>::m_edge_eta);
  HairFlow<DIM>::m_max_eta = HairFlow<DIM>::m_eta.maxCoeff();
  HairFlow<DIM>::m_min_eta = HairFlow<DIM>::m_eta.minCoeff();
}

template<int DIM>
void CylindricalShallowFlow<DIM>::updateHairFlowHeight(const scalar& dt)
{
  const int np = HairFlow<DIM>::m_eta.size();
  const PolygonalCohesion<DIM>* cohesion = HairFlow<DIM>::m_parent->getPolygonalCohesion();
  VectorXs& v = HairFlow<DIM>::m_parent->getV();
  VectorXs& m = HairFlow<DIM>::m_parent->getM();
  
  if(HairFlow<DIM>::m_parent->isIndividualTransfer()) {
    for(int i = 1; i < np - 1; ++i) m_rhs_plus.row(i) = m_rhs.row(i);
  } else {
    const MatrixXs& inter_rhs = cohesion->getRhsOffsetVGlobal();
    for(int i = 1; i < np - 1; ++i) m_rhs_plus.row(i) = m_rhs.row(i) + inter_rhs.row(HairFlow<DIM>::m_particle_indices[i]);
  }

  m_rhs_plus.row(0).setZero();
  m_rhs_plus.row(m_rhs_plus.rows() - 1).setZero();
  
  const scalar& rho_liq = HairFlow<DIM>::m_parent->getLiquidDensity();
  
  if(m_solver.info() == Eigen::Success)
  {
    m_sol = m_solver.solve(m_rhs_plus);
    
    MASS_UPDATE_MODE mum = HairFlow<DIM>::m_parent->getMassUpdateMode();
    const scalar& liquid_shell = HairFlow<DIM>::m_parent->getLiquidShell();
    
    if(mum == MUM_NONE) {
      for(int i = 0; i < np; ++i)
      {
        m_sol(i, 0) = std::max(0.0, m_sol(i, 0));
        
        HairFlow<DIM>::m_eta(i) = sqrt(std::max(m_sol(i, 0) + HairFlow<DIM>::m_rad_vec(i) * HairFlow<DIM>::m_rad_vec(i), HairFlow<DIM>::m_rad_vec(i) * HairFlow<DIM>::m_rad_vec(i))) - HairFlow<DIM>::m_rad_vec(i);
      }
    } else if(mum == MUM_MASS_ONLY) {
      for(int i = 0; i < np; ++i)
      {
        int pidx = HairFlow<DIM>::m_particle_indices[i];
        int idof = HairFlow<DIM>::m_parent->getDof(pidx);
        m_sol(i, 0) = std::max(0.0, m_sol(i, 0));
        
        HairFlow<DIM>::m_eta(i) = sqrt(m_sol(i, 0) + HairFlow<DIM>::m_rad_vec(i) * HairFlow<DIM>::m_rad_vec(i)) - HairFlow<DIM>::m_rad_vec(i);
        
        scalar H = std::max(HairFlow<DIM>::m_eta(i) + HairFlow<DIM>::m_rad_vec(i), liquid_shell * HairFlow<DIM>::m_rad_vec(i));
        scalar modified_A = H * H - HairFlow<DIM>::m_rad_vec(i) * HairFlow<DIM>::m_rad_vec(i);
        
        const scalar& rho_hair = HairFlow<DIM>::m_parent->getHairDensity(pidx);
        const scalar w = HairFlow<DIM>::m_rad_vec(i) * HairFlow<DIM>::m_rad_vec(i) * rho_hair + modified_A * rho_liq;
        if(w != 0.0) {
          m.segment<DIM>(idof).setConstant(w * M_PI * HairFlow<DIM>::m_area_v_hair(i));
        }
        
        if(i != np - 1 && !HairFlow<DIM>::m_parent->isMassSpring()) {
          const scalar m_hair = HairFlow<DIM>::m_rad_vec(i) * HairFlow<DIM>::m_rad_vec(i) * rho_hair;
          const scalar I_hair = m_hair * HairFlow<DIM>::m_rad_vec(i) * HairFlow<DIM>::m_rad_vec(i) * 0.5;
          const scalar m_liq = modified_A * rho_liq;
          const scalar I_liq = m_liq * (modified_A * 0.5 + HairFlow<DIM>::m_rad_vec(i) * HairFlow<DIM>::m_rad_vec(i));
          
          const scalar I = I_hair + I_liq;
          if(I != 0.0) {
            m(idof + 3) = I * M_PI * HairFlow<DIM>::m_area_v_hair(i);
          }
        }
      }
    } else if(mum == MUM_DIRECT_DIV) {
      for(int i = 0; i < np; ++i)
      {
        int pidx = HairFlow<DIM>::m_particle_indices[i];
        int idof = HairFlow<DIM>::m_parent->getDof(pidx);
        m_sol(i, 0) = std::max(0.0, m_sol(i, 0));
        
        HairFlow<DIM>::m_eta(i) = sqrt(m_sol(i, 0) + HairFlow<DIM>::m_rad_vec(i) * HairFlow<DIM>::m_rad_vec(i)) - HairFlow<DIM>::m_rad_vec(i);
        
        scalar H = std::max(HairFlow<DIM>::m_eta(i) + HairFlow<DIM>::m_rad_vec(i), liquid_shell * HairFlow<DIM>::m_rad_vec(i));
        scalar modified_A = H * H - HairFlow<DIM>::m_rad_vec(i) * HairFlow<DIM>::m_rad_vec(i);
        
        const scalar& rho_hair = HairFlow<DIM>::m_parent->getHairDensity(pidx);
        const scalar w = HairFlow<DIM>::m_rad_vec(i) * HairFlow<DIM>::m_rad_vec(i) * rho_hair + modified_A * rho_liq;
        if(w != 0.0) {
          if(!HairFlow<DIM>::m_parent->isFixed(pidx))
            v.segment<DIM>(idof) *= m(idof) / (w * M_PI * HairFlow<DIM>::m_area_v_hair(i));
          m.segment<DIM>(idof).setConstant(w * M_PI * HairFlow<DIM>::m_area_v_hair(i));
        }
        
        if(i != np - 1 && !HairFlow<DIM>::m_parent->isMassSpring()) {
          const scalar m_hair = HairFlow<DIM>::m_rad_vec(i) * HairFlow<DIM>::m_rad_vec(i) * rho_hair;
          const scalar I_hair = m_hair * HairFlow<DIM>::m_rad_vec(i) * HairFlow<DIM>::m_rad_vec(i) * 0.5;
          const scalar m_liq = modified_A * rho_liq;
          const scalar I_liq = m_liq * (modified_A * 0.5 + HairFlow<DIM>::m_rad_vec(i) * HairFlow<DIM>::m_rad_vec(i));
          
          const scalar I = I_hair + I_liq;
          if(I != 0.0) {
            if(!HairFlow<DIM>::m_parent->isFixed(pidx))
              v(idof + 3) *= m(idof + 3) / (I * M_PI * HairFlow<DIM>::m_area_v_hair(i));
            m(idof + 3) = I * M_PI * HairFlow<DIM>::m_area_v_hair(i);
          }
        }
      }
    } else if(mum == MUM_MOMENTUM) {
      for(int i = 0; i < np; ++i)
      {
        int pidx = HairFlow<DIM>::m_particle_indices[i];
        int idof = HairFlow<DIM>::m_parent->getDof(pidx);
        m_sol(i, 0) = std::max(0.0, m_sol(i, 0));
        
        HairFlow<DIM>::m_eta(i) = sqrt(m_sol(i, 0) + HairFlow<DIM>::m_rad_vec(i) * HairFlow<DIM>::m_rad_vec(i)) - HairFlow<DIM>::m_rad_vec(i);
        
        scalar H = std::max(HairFlow<DIM>::m_eta(i) + HairFlow<DIM>::m_rad_vec(i), liquid_shell * HairFlow<DIM>::m_rad_vec(i));
        scalar modified_A = H * H - HairFlow<DIM>::m_rad_vec(i) * HairFlow<DIM>::m_rad_vec(i);
        
        const scalar& rho_hair = HairFlow<DIM>::m_parent->getHairDensity(pidx);
        const Vectors<DIM> p_hair = v.segment<DIM>(idof) * HairFlow<DIM>::m_rad_vec(i) * HairFlow<DIM>::m_rad_vec(i) * rho_hair;
        const Vectors<DIM> p_liq = m_sol.block<1, DIM>(i, 1) * rho_liq;
        const scalar w = HairFlow<DIM>::m_rad_vec(i) * HairFlow<DIM>::m_rad_vec(i) * rho_hair + modified_A * rho_liq;
        if(w != 0.0) {
          if(!HairFlow<DIM>::m_parent->isFixed(pidx))
            v.segment<DIM>(idof) = (p_hair + p_liq) / w;
          m.segment<DIM>(idof).setConstant(w * M_PI * HairFlow<DIM>::m_area_v_hair(i));
        }
        
        if(i != np - 1 && !HairFlow<DIM>::m_parent->isMassSpring()) {
          const scalar m_hair = HairFlow<DIM>::m_rad_vec(i) * HairFlow<DIM>::m_rad_vec(i) * rho_hair;
          const scalar I_hair = m_hair * HairFlow<DIM>::m_rad_vec(i) * HairFlow<DIM>::m_rad_vec(i) * 0.5;
          const scalar L_hair = v(idof + 3) * I_hair;
          const scalar m_liq = modified_A * rho_liq;
          const scalar I_liq = m_liq * (modified_A * 0.5 + HairFlow<DIM>::m_rad_vec(i) * HairFlow<DIM>::m_rad_vec(i));
          const scalar L_liq = m_sol(i, 4) * rho_liq;
          
          const scalar I = I_hair + I_liq;
          if(I != 0.0) {
            if(!HairFlow<DIM>::m_parent->isFixed(pidx))
              v(idof + 3) = (L_hair + L_liq) / I;
            m(idof + 3) = I * M_PI * HairFlow<DIM>::m_area_v_hair(i);
          }
        }
      }
    }
  } else {
    std::cout << "Solver Error: " << m_solver.info() << std::endl;
    std::cout << "LHS: " << std::endl;
    std::cout << m_LHS << std::endl;
    std::cout << "rhs: " << std::endl;
    std::cout << m_rhs << std::endl;
    std::cout << "u: " << std::endl;
    std::cout << HairFlow<DIM>::m_u << std::endl;
    
    std::cout << "po: " << std::endl;
    std::cout << m_porosity_e << std::endl;
    exit(1);
  }
  
  
}

template<int DIM>
void CylindricalShallowFlow<DIM>::updateHairMass()
{
  const int np = HairFlow<DIM>::m_eta.size();
  
  compute_liquid_mass(HairFlow<DIM>::m_eta, HairFlow<DIM>::m_area_v_hair, HairFlow<DIM>::m_rad_vec, HairFlow<DIM>::m_parent->getLiquidDensity(), m_liquid_mass, HairFlow<DIM>::m_parent->getLiquidShell());
  
  for(int i = 0; i < np; ++i)
  {
    HairFlow<DIM>::m_parent->updateMassWithLiquid(HairFlow<DIM>::m_particle_indices[i], m_liquid_mass(i), HairFlow<DIM>::m_eta(i) + HairFlow<DIM>::m_rad_vec(i));
  }
}


template<int DIM>
void CylindricalShallowFlow<DIM>::resizeSystem()
{
  const int num_particles = HairFlow<DIM>::m_particle_indices.size();
  const int num_edges = HairFlow<DIM>::m_internal_edges.size();
  
  HairFlow<DIM>::m_area_v.resize(num_particles);
  HairFlow<DIM>::m_area_v_hair.resize(num_particles);
  HairFlow<DIM>::m_area_e.resize(num_edges);
  
  HairFlow<DIM>::m_dir_f.resize(num_edges, DIM);
  HairFlow<DIM>::m_dir_v.resize(num_particles, DIM);
  
  HairFlow<DIM>::m_actual_u_f.resize(num_edges, DIM);
  HairFlow<DIM>::m_actual_u_v.resize(num_particles, DIM);
  
  HairFlow<DIM>::m_edge_rad_vec.resize(num_edges);
  HairFlow<DIM>::m_rad_vec.resize(num_particles);
  HairFlow<DIM>::m_stored_friction_coeff.resize(num_edges);
  
  m_area_v_hair_flat.resize(1, num_particles);
  m_Gv.resize(num_particles, num_particles);
  m_iGv.resize(num_particles, num_particles);
  m_Gv_hair.resize(num_particles, num_particles);
  m_iGv_hair.resize(num_particles, num_particles);
  m_Gf.resize(num_edges * DIM, num_edges * DIM);
  m_W_fv.resize(num_particles, num_edges);
  m_W_fv_hair.resize(num_particles, num_edges);
  m_gradF.resize(num_edges * DIM, num_particles);
  m_divV.resize(num_particles, num_edges * DIM);
  m_dir_f_expand.resize(num_edges, num_edges * DIM);
  m_D.resize(num_particles, num_particles);
  m_A.resize(num_particles, num_particles);
  m_LHS.resize(num_particles, num_particles);
  m_L.resize(num_particles, num_particles);
  m_sp0.resize(num_particles, num_particles);
  m_sp1.resize(num_particles, num_particles);
  
  HairFlow<DIM>::m_accel_v.resize(num_particles, DIM);
  m_accel_f.resize(num_edges, DIM);
  m_porosity_e.resize(num_edges);
  m_pressure.resize(num_particles);
  
  HairFlow<DIM>::m_c_v.resize(num_particles, DIM * DIM);
  HairFlow<DIM>::m_c_star.resize(num_edges, DIM * DIM);
  
  m_edge_x.resize(num_edges);
  HairFlow<DIM>::m_u.resize(num_edges);
  m_ustar.resize(num_edges);
  m_accel_e.resize(num_edges);
  m_cap_e.resize(num_edges);
  m_cap_e.setZero();
  HairFlow<DIM>::m_edge_eta.resize(num_edges);
  
  m_area_e_accu.resize(num_particles);
  m_area_e_inv_mapping.resize(num_particles);
  m_uvert.resize(num_particles);
  
  m_old_rad_vec.resize(num_particles);
  m_old_area_v.resize(num_particles);
  m_delta_sol.resize(num_particles);
  
  if(HairFlow<DIM>::m_parent->isMassSpring()) {
    m_rhs.resize(num_particles, DIM + 1);
    m_rhs_plus.resize(num_particles, DIM + 1);
    m_sol.resize(num_particles, DIM + 1);
  } else {
    m_rhs.resize(num_particles, DIM + 2);
    m_rhs_plus.resize(num_particles, DIM + 2);
    m_sol.resize(num_particles, DIM + 2);
  }

  m_liquid_mass.resize(num_particles);
  m_mass_v.resize(num_particles);
  
  m_old_eta.resize(num_particles);
  
  HairFlow<DIM>::m_edge_bridges.resize(num_edges);
}

template<int DIM>
Force* CylindricalShallowFlow<DIM>::createNewCopy()
{
	return NULL;
}

template<int DIM>
const char* CylindricalShallowFlow<DIM>::name()
{
  return cylindrical_film_static_name;
}

// explicit instantiations at bottom
template class CylindricalShallowFlow<2>;
template class CylindricalShallowFlow<3>;
