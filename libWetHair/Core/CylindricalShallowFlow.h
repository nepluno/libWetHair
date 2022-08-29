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

#ifndef LIBWETHAIR_CORE_CYLINDRICAL_SHALLOW_FLOW_H_
#define LIBWETHAIR_CORE_CYLINDRICAL_SHALLOW_FLOW_H_

#include <Eigen/Core>
#include <cassert>
#include <iostream>
#include <unordered_map>
#include <vector>

#include "Force.h"
#include "HairFlow.h"
#include "MathDefs.h"
#include "MathUtilities.h"

namespace libwethair {
namespace mathutils {
namespace shallowflow {

template <int DIM>
inline void compute_edge_gradient_matrix(
    const VectorXs& x, const std::vector<std::pair<int, int> >& edges,
    const std::vector<std::pair<int, int> >& edge_local, const VectorXs& area_e,
    SparseXs& gradF, TwoDScene<DIM>* parent) {
  const int ne = edges.size();

  TripletXs tri;

  for (int i = 0; i < ne; ++i) {
    int i0 = edges[i].first;
    int i1 = edges[i].second;
    const VectorXs& x0 = x.segment<DIM>(parent->getDof(i0));
    const VectorXs& x1 = x.segment<DIM>(parent->getDof(i1));
    VectorXs dx = x1 - x0;

    scalar coeff = 1.0 / (area_e(i) * area_e(i));

    int local_i0 = edge_local[i].first;
    int local_i1 = edge_local[i].second;

    for (int r = 0; r < DIM; ++r) {
      tri.push_back(Triplets(i * DIM + r, local_i0, -coeff * dx(r)));
      tri.push_back(Triplets(i * DIM + r, local_i1, coeff * dx(r)));
    }
  }

  gradF.setFromTriplets(tri.begin(), tri.end());
}

template <int DIM>
inline void compute_vertex_curvature(const VectorXs& x,
                                     const std::vector<int>& particles,
                                     VectorXs& kappa, TwoDScene<DIM>* parent) {
  int np = particles.size();
  for (int i = 0; i < np; ++i) {
    int idx = particles[i];
    if (i == 0 || i == np - 1)
      kappa(i) = 0.0;
    else {
      int idx_prev = particles[i - 1];
      int idx_next = particles[i + 1];

      Vectors<DIM> ei = x.segment<DIM>(parent->getDof(idx)) -
                        x.segment<DIM>(parent->getDof(idx_prev));
      Vectors<DIM> ej = x.segment<DIM>(parent->getDof(idx)) -
                        x.segment<DIM>(parent->getDof(idx_next));

      scalar curvature =
          2.0 * unicross(ei, ej) / (ei.norm() * ej.norm() + fabs(ei.dot(ej)));

      kappa(i) = curvature;
    }
  }
}

inline void compute_edge_val(const VectorXs& vert_val, const SparseXs& W_fv,
                             VectorXs& edge_val) {
  edge_val = W_fv.transpose() * vert_val;
}

inline void compute_edge_val(const MatrixXs& vert_val, const SparseXs& W_fv,
                             MatrixXs& edge_val) {
  edge_val = W_fv.transpose() * vert_val;
}

inline void compute_edge_val(const VectorXs& vert_val,
                             const std::vector<std::pair<int, int> >& edges,
                             VectorXs& edge_val) {
  int ne = edges.size();

  for (int i = 0; i < ne; ++i) {
    int i0 = edges[i].first;
    int i1 = edges[i].second;
    edge_val(i) = (vert_val(i0) + vert_val(i1)) * 0.5;
  }
}

inline void compute_edge_val(const MatrixXs& vert_val,
                             const std::vector<std::pair<int, int> >& edges,
                             MatrixXs& edge_val) {
  int ne = edges.size();

  for (int i = 0; i < ne; ++i) {
    int i0 = edges[i].first;
    int i1 = edges[i].second;
    edge_val.row(i) = (vert_val.row(i0) + vert_val.row(i1)) * 0.5;
  }
}

template <int DIM>
inline void compute_accel_v(const VectorXs& accel,
                            const std::vector<int>& particles,
                            const VectorXs& g, MatrixXs& accel_f,
                            TwoDScene<DIM>* parent) {
  int np = particles.size();
  for (int i = 0; i < np; ++i) {
    int idx = particles[i];
    accel_f.row(i) =
        (-accel.segment<DIM>(parent->getDof(idx)) + g.segment<DIM>(0))
            .transpose();
  }
}

inline void compute_capillary_accel_e(const VectorXs& radii_v,
                                      const VectorXs& eta_v,
                                      const VectorXs& accu_v,
                                      const scalar& gamma, const scalar& Lmax,
                                      const scalar& multiplier,
                                      VectorXs& cap_e) {
  int ne = cap_e.size();
  for (int i = 0; i < ne; ++i) {
    scalar invL = (accu_v(i) + accu_v(i + 1)) * 0.5;
    scalar H0 = eta_v(i) + radii_v(i);
    scalar H1 = eta_v(i + 1) + radii_v(i + 1);
    scalar R0 = sqrt((H0 * H0 + H1 * H1) * 0.5);
    scalar rf =
        sqrt((radii_v(i) * radii_v(i) + radii_v(i + 1) * radii_v(i + 1)) * 0.5);
    scalar L = Lmax - invL;
    if (L <= 0.0 || L >= 1.0 || R0 + rf == 0.0)
      cap_e(i) = 0.0;
    else
      cap_e(i) = -4.0 * gamma * rf / L / ((R0 + rf) * (R0 + rf)) * multiplier;
  }
}

inline scalar compute_accumulated_edge_length(const VectorXs& area_e,
                                              VectorXs& area_e_accu) {
  int np = area_e_accu.size();
  int ne = area_e.size();
  scalar sum = 0.0;
  for (int i = 0; i < np; ++i) {
    area_e_accu(i) = sum;
    if (i < ne) sum += area_e(i);
  }

  return sum;
}

inline void compute_inverted_pos_mapping(const VectorXs& area_e_accu,
                                         const scalar& sum_area,
                                         VectorXs& area_e_inv_mapping) {
  const int N_disc = area_e_inv_mapping.size();
  const scalar incr = sum_area / N_disc;

  int last_update = 0;

  const int np = area_e_accu.size();

  for (int i = 1; i <= np; ++i) {
    scalar grid_loc = ((i == np) ? sum_area : area_e_accu(i)) / incr;
    scalar prev_grid_loc = area_e_accu(i - 1) / incr;

    int i_grid_loc = (int)ceil(grid_loc);

    for (int j = last_update; j < i_grid_loc && j < N_disc; ++j) {
      scalar frac = ((scalar)j - prev_grid_loc) / (grid_loc - prev_grid_loc);
      area_e_inv_mapping(j) = (scalar)(i - 1) + frac;
    }

    last_update = i_grid_loc;
  }
}

inline scalar find_pos_on_hair(const scalar& val, const scalar& sum_area,
                               const VectorXs& area_e_accu,
                               const VectorXs& area_e_inv_mapping) {
  const int N_disc = area_e_inv_mapping.size();
  const scalar incr = sum_area / N_disc;

  scalar grid_pos = val / incr;

  const int np = area_e_accu.size();

  int lowb =
      (int)area_e_inv_mapping(std::max(0, std::min(np - 1, (int)grid_pos)));
  int highb = (int)ceil(
      area_e_inv_mapping(std::max(0, std::min(np - 1, (int)ceil(grid_pos)))));

  int k = highb;
  while (k > lowb && area_e_accu(k) > val) {
    --k;
  }

  scalar lowv = area_e_accu(k);
  scalar highv = (k == np - 1) ? sum_area : area_e_accu(k + 1);
  if (lowv == highv)
    return (scalar)k;
  else
    return (scalar)k + (val - lowv) / (highv - lowv);
}

inline scalar interpolate(const scalar& pos, const scalar& sum_area,
                          const VectorXs& area_e_accu,
                          const VectorXs& area_e_inv_mapping,
                          const VectorXs& source) {
  scalar pos_hair =
      find_pos_on_hair(pos, sum_area, area_e_accu, area_e_inv_mapping);
  int np = source.size();
  int ip = std::max(0, std::min(np - 2, (int)pos_hair));
  scalar frac = pos_hair - floor(pos_hair);

  return lerp(source(ip), source(ip + 1), frac);
}

inline void interpolate(const scalar& pos, const scalar& sum_area,
                        const VectorXs& area_e_accu,
                        const VectorXs& area_e_inv_mapping,
                        const MatrixXs& source, VectorXsT& target) {
  scalar pos_hair =
      find_pos_on_hair(pos, sum_area, area_e_accu, area_e_inv_mapping);
  int np = source.rows();  // size();
  int ip = std::max(0, std::min(np - 2, (int)pos_hair));
  scalar frac = pos_hair - floor(pos_hair);

  target = source.row(ip) * (1.0 - frac) + source.row(ip + 1) * frac;
}

inline scalar trace_rk2(const scalar& pos, const scalar& dt,
                        const scalar& sum_area, const VectorXs& area_e_accu,
                        const VectorXs& area_e_inv_mapping,
                        const VectorXs& u_vert) {
  scalar input = pos;
  scalar velocity =
      interpolate(input, sum_area, area_e_accu, area_e_inv_mapping, u_vert);
  velocity = interpolate(input + 0.5 * dt * velocity, sum_area, area_e_accu,
                         area_e_inv_mapping, u_vert);
  input += dt * velocity;
  return input;
}

template <int DIM>
inline void compute_backtracing(const VectorXs& edge_x, const scalar& dt,
                                const scalar& sum_area,
                                const VectorXs& area_e_accu,
                                const VectorXs& area_e_inv_mapping,
                                const VectorXs& u_vert, const MatrixXs& c_vert,
                                VectorXs& u_star, MatrixXs& c_star) {
  int ne = u_star.size();
  VectorXsT buf;

  for (int i = 0; i < ne; ++i) {
    scalar pos = edge_x(i);
    pos =
        trace_rk2(pos, -dt, sum_area, area_e_accu, area_e_inv_mapping, u_vert);
    u_star(i) = interpolate(pos, sum_area, area_e_accu, area_e_inv_mapping,
                            u_vert);  // * porosity_e(i);
    interpolate(pos, sum_area, area_e_accu, area_e_inv_mapping, c_vert, buf);
    c_star.row(i) = buf;
  }
}

inline void compute_accel_e(const MatrixXs& accel, const MatrixXs& dir,
                            VectorXs& accel_e) {
  int np = accel_e.size();

  for (int i = 0; i < np; ++i) {
    accel_e(i) = accel.row(i).dot(dir.row(i));
  }
}

inline void compute_next_velocity(
    const VectorXs& u_star, const SparseXs& gradF, const SparseXs& dir_f_exp,
    const VectorXs& pressure, const VectorXs& porosity_e,
    const VectorXs& eta_edge, const VectorXs& radii_edge, const scalar& dt,
    const VectorXs& accel_e, const scalar& rho, const scalar& visc,
    const scalar& friction, VectorXs& u_next) {
  u_next = u_star + (accel_e - (1.0 / rho) * dir_f_exp * gradF * pressure) * dt;

  int np = u_next.size();

  if (friction == 0.0 || visc == 0.0) return;

  for (int i = 0; i < np; ++i) {
    // Here we generalize [Gerbeau & Perthame 2000] for porous flow
    scalar H = eta_edge(i) + radii_edge(i);
    scalar A = M_PI * H * H * porosity_e(i);
    scalar vb = 3.0 * visc * A + friction * A * A;
    scalar vc = vb / (vb + 3.0 * friction * dt * visc);
    u_next(i) *= vc;
  }
}

inline void compute_LHS_matrix(const SparseXs& sp0, const SparseXs& sp1,
                               const SparseXs& A, const scalar& dt,
                               SparseXs& lhs) {
  lhs = sp1 + sp0 * A * dt;
}

inline void compute_SP0_matrix(SparseXs& sp0) {
  TripletXs tri;
  int np = sp0.rows();
  tri.reserve(np);
  for (int i = 1; i < np - 1; ++i) {
    tri.push_back(Triplets(i, i, 1.0));
  }

  sp0.setFromTriplets(tri.begin(), tri.end());
}

inline void compute_SP1_matrix(SparseXs& sp1) {
  TripletXs tri;
  int np = sp1.rows();
  tri.reserve(np);
  for (int i = 0; i < np; ++i) {
    tri.push_back(Triplets(i, i, 1.0));
  }
  tri.push_back(Triplets(0, 1, -1.0));
  tri.push_back(Triplets(np - 1, np - 2, -1.0));

  sp1.setFromTriplets(tri.begin(), tri.end());
}

inline void compute_rhs_vector(const std::vector<int>& particle_indices,
                               const VectorXs& hair_vel, const VectorXs& h,
                               const VectorXs& r, VectorXs& rhs) {
  for (int i = 0; i < rhs.size(); ++i) {
    rhs(i, 0) = (h(i) + r(i)) * (h(i) + r(i));

    int pidx = particle_indices[i];
  }
}

inline void compute_stretched_liquid_height(const VectorXs& old_rad,
                                            const VectorXs& new_rad,
                                            const VectorXs& old_area,
                                            const VectorXs& new_area,
                                            VectorXs& eta) {
  int np = eta.size();
  for (int i = 0; i < np; ++i) {
    scalar na = new_area(i);
    if (na == 0.0) continue;

    eta(i) = sqrt(
        std::max(0.0, new_rad(i) * new_rad(i) +
                          old_area(i) / new_area(i) *
                              (eta(i) * eta(i) - old_rad(i) * old_rad(i))));
  }
}

inline void compute_liquid_mass(const VectorXs& eta, const VectorXs& area,
                                const VectorXs& rad, const scalar& rho,
                                VectorXs& mass, const scalar& shell) {
  int np = eta.size();
  for (int i = 0; i < np; ++i) {
    scalar H = std::max(eta(i) + rad(i), shell * rad(i));
    mass(i) = rho * M_PI * area(i) * (H * H - rad(i) * rad(i));
  }
}

inline void compute_liquid_mass(const scalar max_eta_prop, const VectorXs& area,
                                const VectorXs& rad, const scalar& rho,
                                VectorXs& mass) {
  int np = rad.size();
  for (int i = 0; i < np; ++i) {
    mass(i) = rho * M_PI * area(i) * (max_eta_prop * max_eta_prop - 1.0) *
              rad(i) * rad(i);
  }
}

template <int DIM>
inline void copy_rest_mass(const VectorXs& rest_mass,
                           const std::vector<int>& indices, VectorXs& mass_v,
                           TwoDScene<DIM>* parent) {
  int np = mass_v.size();

  for (int i = 0; i < np; ++i) {
    mass_v(i) = rest_mass(parent->getDof(indices[i]));
  }
}

inline void compute_rest_length(const VectorXs& global_rest_e,
                                const std::vector<int>& edge_idx,
                                VectorXs& rest_e) {
  int np = rest_e.size();
  for (int i = 0; i < np; ++i) {
    rest_e(i) = global_rest_e(edge_idx[i]);
  }
}

inline scalar integrate_volume_difference(
    const VectorXs& eta, const VectorXs& old_eta, const VectorXs& area_e,
    const std::vector<std::pair<int, int> >& edges) {
  int ne = edges.size();

  scalar sum = 0.0;
  for (int i = 0; i < ne; ++i) {
    auto& e = edges[i];
    scalar eta0 = eta(e.first);
    scalar oldeta0 = old_eta(e.first);
    scalar area_diff0 = (eta0 * eta0 - oldeta0 * oldeta0);

    scalar eta1 = eta(e.second);
    scalar oldeta1 = old_eta(e.second);
    scalar area_diff1 = (eta1 * eta1 - oldeta1 * oldeta1);

    sum += (area_diff0 + area_diff1) * area_e(i) * 0.5;
  }

  return sum;
}

inline void compute_actual_u(const VectorXs& u, const MatrixXs& dir,
                             MatrixXs& actual_u) {
  int np = dir.rows();
  for (int i = 0; i < np; ++i) {
    actual_u.row(i) = u(i) * dir.row(i).normalized();
  }
}

template <int DIM>
inline int get_max_normal_velocity_index(const MatrixXs& actual_u_v,
                                         const MatrixXs& dir_v) {
  int np = actual_u_v.rows();
  scalar max_norm_v = 0.0;
  int index = np - 1;
  for (int i = 0; i < np; ++i) {
    scalar vv;
    if (i == 0) {
      vv = std::max(0.0, -dir_v.row(0).dot(actual_u_v.row(0)));
    } else if (i == np - 1) {
      vv = std::max(0.0, dir_v.row(i).dot(actual_u_v.row(i)));
    } else {
      vv = (actual_u_v.row(i) * (Matrixs<DIM>::Identity() -
                                 dir_v.row(i).transpose() * dir_v.row(i)))
               .norm();
    }
    if (vv > max_norm_v) {
      index = i;
      max_norm_v = vv;
    }
  }

  return index;
}

inline void compute_rest_porosity(const VectorXs& eta,
                                  const VectorXs& global_radii,
                                  const std::vector<int>& indices,
                                  VectorXs& porosity) {
  int np = eta.size();
  for (int i = 0; i < np; ++i) {
    scalar ra = global_radii(indices[i]);
    porosity(i) = 1.0 - mathutils::clamp(ra * ra / (eta(i) * eta(i)), 0.0, 1.0);
  }
}

inline void expand_dir_f(const MatrixXs& dir_f, SparseXs& dir_f_exp) {
  TripletXs tris;
  tris.reserve(dir_f.size());

  int ne = dir_f.rows();
  int dim = dir_f.cols();

  for (int i = 0; i < ne; ++i) {
    for (int j = 0; j < dim; ++j) {
      tris.push_back(Triplets(i, i * dim + j, dir_f(i, j)));
    }
  }

  dir_f_exp.setFromTriplets(tris.begin(), tris.end());
}
};  // namespace shallowflow
};  // namespace mathutils

template <int DIM>
class CylindricalShallowFlow : public HairFlow<DIM> {
 public:
  CylindricalShallowFlow(TwoDScene<DIM>* parent,
                         const std::vector<int>& involved_particles,
                         const VectorXs& eta,
                         const std::vector<unsigned char>& particle_state);

  virtual void updateGeometricState(const VectorXs& x, const VectorXs& v,
                                    FluidSim* fluidsim);

  virtual void updateToFilteredGrid(const VectorXs& x, const VectorXs& v,
                                    FluidSim* fluidsim, const scalar& dt,
                                    int ibuffer);

  virtual void updateFromFilteredGrid(const VectorXs& x, VectorXs& v,
                                      FluidSim* fluidsim, const scalar& dt);

  virtual void advance(const VectorXs& x, const scalar& dt);

  virtual void add_force(const VectorXs& x, const VectorXs& accel,
                         FluidSim* fluidsim, const scalar& dt);

  virtual void updateHairMass();

  virtual void updateHairFlowHeight(const scalar& dt);

  virtual void preUpdateHairFlowHeight(const scalar& dt);

  virtual void postUpdateHairFlowHeight(const scalar& dt);

  virtual void resizeSystem();

  virtual Force* createNewCopy();

  virtual const char* name();

  virtual void updateReservoir(FluidSim* fluidsim, const VectorXs& x,
                               const VectorXs& v, const scalar& dt);

  virtual Vector3s computeHairLiquidAngularMomentum(const VectorXs& x,
                                                    const VectorXs& v,
                                                    FluidSim* fluidsim) const;

  virtual scalar computeTotalLiquidVol() const;

  virtual scalar computeTotalReservoirVol() const;

  virtual const scalar& getPoolSize() const;

  virtual scalar& getPoolSize();

  virtual void geodesic_to_local(const scalar& geopos, int& pidx_low,
                                 scalar& alpha) const;

  virtual scalar getTotalLength() const;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
 protected:
  SparseXs m_Gv;
  SparseXs m_iGv;
  SparseXs m_Gv_hair;
  SparseXs m_iGv_hair;
  SparseXs m_Gf;
  SparseXs m_W_fv;
  SparseXs m_W_fv_hair;
  SparseXs m_gradF;
  SparseXs m_divV;
  SparseXs m_divV_hair;
  SparseXs m_L;
  SparseXs m_dir_f_expand;
  SparseXs m_D;
  SparseXs m_A;
  SparseXs m_LHS;
  SparseXs m_sp0;
  SparseXs m_sp1;

  SparseXs
      m_area_v_hair_flat;  // used by igl::active_set for volume computation

  MatrixXs m_accel_f;

  VectorXs m_delta_sol;  // used by igl::active_set for result
  VectorXs m_edge_x;
  VectorXs m_area_e_accu;
  VectorXs m_area_e_inv_mapping;
  VectorXs m_uvert;
  VectorXs m_ustar;
  VectorXs m_accel_e;
  VectorXs m_cap_e;
  MatrixXs m_rhs;
  MatrixXs m_rhs_plus;
  MatrixXs m_sol;
  VectorXs m_mass_v;

  VectorXs m_porosity_e;
  VectorXs m_pressure;

  VectorXs m_old_eta;
  VectorXs m_old_area_v;
  VectorXs m_old_rad_vec;

  VectorXs m_liquid_mass;
  Eigen::SparseLU<SparseXs> m_solver;

  scalar m_sum_area_e;
  scalar m_pool_liquid;
};

}  // namespace libwethair

#endif  // LIBWETHAIR_CORE_CYLINDRICAL_SHALLOW_FLOW_H_
