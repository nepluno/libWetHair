//
// This file is part of the libWetHair open source project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright 2008 Christopher Batty, Robert Bridson
//

#include "viscosity3d.h"

#include <tbb/tbb.h>

#include <cmath>
#include <fstream>

#include "MathDefs.h"
#include "array3.h"
#include "pcgsolver/pcg_solver.h"
#include "pcgsolver/sparse_matrix.h"

namespace libwethair {
using namespace robertbridson;

Array3c u_state;  //(nx+1,ny,nz,0);
Array3c v_state;  //(nx,ny+1,nz,0);
Array3c w_state;  //(nx,ny,nz+1,0);

SparseMatrixd matrix;      //(dim,dim,15);
std::vector<double> rhs;   //(dim);
std::vector<double> soln;  //(dim);
std::vector<Vector3i> dof_ijk;
Array3ui dof_index;

inline unsigned u_ind(int i, int j, int k, int nx, int ny) {
  return dof_index(i, j, k);
}

inline unsigned v_ind(int i, int j, int k, int nx, int ny, int nz) {
  return dof_index(i, j, k + nz);
}

inline unsigned w_ind(int i, int j, int k, int nx, int ny, int nz) {
  return dof_index(i, j, k + nz * 2);
}

void advance_viscosity_implicit_weighted(
    Array3s& u, Array3s& v, Array3s& w, Array3s& u_visc_impulse,
    Array3s& v_visc_impulse, Array3s& w_visc_impulse, const Array3s& vol_u,
    const Array3s& vol_v, const Array3s& vol_w, const Array3s& vol_c,
    const Array3s& vol_ex, const Array3s& vol_ey, const Array3s& vol_ez,
    const Array3s& solid_phi, const scalar viscosity, scalar dt, scalar dx) {
  scalar over_dx = 1.0 / dx;
  int nx = solid_phi.ni;
  int ny = solid_phi.nj;
  int nz = solid_phi.nk;

  int dim = (nx + 1) * ny * nz + nx * (ny + 1) * nz + nx * ny * (nz + 1);
  if (u_state.ni != u.ni) {
    u_state.resize(nx + 1, ny, nz);
    v_state.resize(nx, ny + 1, nz);
    w_state.resize(nx, ny, nz + 1);
  }

  dof_ijk.resize(0);
  dof_ijk.reserve(dim);

  dof_index.resize(nx, ny, nz * 3);
  dof_index.assign(0);

  u_state.assign(0);
  v_state.assign(0);
  w_state.assign(0);

  const int SOLID = 3;
  const int FLUID = 2;
  // const int AIR = 1;

  // check if interpolated velocity positions are inside solid
  for (int k = 0; k < nz; ++k)
    for (int j = 0; j < ny; ++j)
      for (int i = 0; i < nx + 1; ++i) {
        if (i - 1 < 0 || i >= nx ||
            solid_phi(i - 1, j, k) + solid_phi(i, j, k) <= 0)
          u_state(i, j, k) = SOLID;
        else
          u_state(i, j, k) = FLUID;
      }

  for (int k = 0; k < nz; ++k)
    for (int j = 0; j < ny + 1; ++j)
      for (int i = 0; i < nx; ++i) {
        if (j - 1 < 0 || j >= ny ||
            solid_phi(i, j - 1, k) + solid_phi(i, j, k) <= 0)
          v_state(i, j, k) = SOLID;
        else
          v_state(i, j, k) = FLUID;
      }

  for (int k = 0; k < nz + 1; ++k)
    for (int j = 0; j < ny; ++j)
      for (int i = 0; i < nx; ++i) {
        if (k - 1 < 0 || k >= nz ||
            solid_phi(i, j, k - 1) + solid_phi(i, j, k) <= 0)
          w_state(i, j, k) = SOLID;
        else
          w_state(i, j, k) = FLUID;
      }

  const int u_base_idx = 0;

  for (int k = 1; k < nz; k++)
    for (int j = 1; j < ny; j++)
      for (int i = 1; i < nx; i++) {
        if (u_state(i, j, k) == FLUID &&
            (vol_u(i, j, k) > 0.0 || vol_c(i, j, k) > 0.0 ||
             vol_c(i - 1, j, k) > 0.0 || vol_ez(i, j + 1, k) > 0.0 ||
             vol_ez(i, j, k) > 0.0 || vol_ey(i, j, k + 1) > 0.0 ||
             vol_ey(i, j, k) > 0.0)) {
          dof_index(i, j, k) = dof_ijk.size();
          dof_ijk.push_back(Vector3i(i, j, k));
        }
      }

  const int v_base_idx = dof_ijk.size();

  for (int k = 1; k < nz; k++)
    for (int j = 1; j < ny; j++)
      for (int i = 1; i < nx; i++) {
        if (v_state(i, j, k) == FLUID &&
            (vol_v(i, j, k) > 0.0 || vol_ez(i + 1, j, k) > 0. ||
             vol_ez(i, j, k) > 0. || vol_c(i, j, k) > 0. ||
             vol_c(i, j - 1, k) > 0. || vol_ex(i, j, k + 1) > 0. ||
             vol_ex(i, j, k) > 0.)) {
          dof_index(i, j, k + nz) = dof_ijk.size();
          dof_ijk.push_back(Vector3i(i, j, k + nz));
        }
      }

  const int w_base_idx = dof_ijk.size();

  for (int k = 1; k < nz; k++)
    for (int j = 1; j < ny; j++)
      for (int i = 1; i < nx; i++) {
        if (w_state(i, j, k) == FLUID &&
            (vol_w(i, j, k) > 0.0 || vol_ey(i + 1, j, k) > 0. ||
             vol_ey(i, j, k) > 0. || vol_ex(i, j + 1, k) > 0. ||
             vol_ex(i, j, k) > 0. || vol_c(i, j, k) > 0. ||
             vol_c(i, j, k - 1) > 0.)) {
          dof_index(i, j, k + nz * 2) = dof_ijk.size();
          dof_ijk.push_back(Vector3i(i, j, k + nz * 2));
        }
      }

  const int total_ijk_size = dof_ijk.size();

  soln.resize(dof_ijk.size());
  soln.assign(dof_ijk.size(), 0);

  if (rhs.size() != soln.size()) {
    rhs.resize(soln.size());
    matrix.resize(soln.size());
  }

  rhs.assign(rhs.size(), 0);
  matrix.zero();

  scalar factor = dt * mathutils::sqr(over_dx);

  threadutils::thread_pool::ParallelFor(u_base_idx, v_base_idx, [&](int index) {
    const int i = dof_ijk[index](0);
    const int j = dof_ijk[index](1);
    const int k = dof_ijk[index](2);

    rhs[index] = vol_u(i, j, k) * u(i, j, k);
    matrix.set_element(index, index, vol_u(i, j, k));

    scalar vol_right = vol_c(i, j, k);
    scalar vol_left = vol_c(i - 1, j, k);

    scalar vol_top = vol_ez(i, j + 1, k);
    scalar vol_bottom = vol_ez(i, j, k);

    scalar vol_front = vol_ey(i, j, k + 1);
    scalar vol_back = vol_ey(i, j, k);

    // u_x_right
    if (vol_right > 0.) {
      matrix.add_to_element(index, index, 2 * factor * viscosity * vol_right);
      if (u_state(i + 1, j, k) == FLUID)
        matrix.add_to_element(index, u_ind(i + 1, j, k, nx, ny),
                              -2 * factor * viscosity * vol_right);
      else if (u_state(i + 1, j, k) == SOLID)
        rhs[index] -= -2 * factor * viscosity * vol_right * u(i + 1, j, k);
    }

    // u_x_left
    if (vol_left > 0.) {
      matrix.add_to_element(index, index, 2 * factor * viscosity * vol_left);
      if (u_state(i - 1, j, k) == FLUID)
        matrix.add_to_element(index, u_ind(i - 1, j, k, nx, ny),
                              -2 * factor * viscosity * vol_left);
      else if (u_state(i - 1, j, k) == SOLID)
        rhs[index] -= -2 * factor * viscosity * vol_left * u(i - 1, j, k);
    }

    // u_y_top
    if (vol_top > 0.) {
      matrix.add_to_element(index, index, +factor * viscosity * vol_top);
      if (u_state(i, j + 1, k) == FLUID)
        matrix.add_to_element(index, u_ind(i, j + 1, k, nx, ny),
                              -factor * viscosity * vol_top);
      else if (u_state(i, j + 1, k) == SOLID)
        rhs[index] -= -u(i, j + 1, k) * factor * viscosity * vol_top;
    }

    // u_y_bottom
    if (vol_bottom > 0.) {
      matrix.add_to_element(index, index, +factor * viscosity * vol_bottom);
      if (u_state(i, j - 1, k) == FLUID)
        matrix.add_to_element(index, u_ind(i, j - 1, k, nx, ny),
                              -factor * viscosity * vol_bottom);
      else if (u_state(i, j - 1, k) == SOLID)
        rhs[index] -= -u(i, j - 1, k) * factor * viscosity * vol_bottom;
    }

    // u_z_front
    if (vol_front > 0.) {
      matrix.add_to_element(index, index, +factor * viscosity * vol_front);
      if (u_state(i, j, k + 1) == FLUID)
        matrix.add_to_element(index, u_ind(i, j, k + 1, nx, ny),
                              -factor * viscosity * vol_front);
      else if (u_state(i, j, k + 1) == SOLID)
        rhs[index] -= -u(i, j, k + 1) * factor * viscosity * vol_front;
    }

    // u_z_back
    if (vol_back > 0.) {
      matrix.add_to_element(index, index, +factor * viscosity * vol_back);
      if (u_state(i, j, k - 1) == FLUID)
        matrix.add_to_element(index, u_ind(i, j, k - 1, nx, ny),
                              -factor * viscosity * vol_back);
      else if (u_state(i, j, k - 1) == SOLID)
        rhs[index] -= -u(i, j, k - 1) * factor * viscosity * vol_back;
    }

    // v_x_top
    if (vol_top > 0.) {
      if (v_state(i, j + 1, k) == FLUID)
        matrix.add_to_element(index, v_ind(i, j + 1, k, nx, ny, nz),
                              -factor * viscosity * vol_top);
      else if (v_state(i, j + 1, k) == SOLID)
        rhs[index] -= -v(i, j + 1, k) * factor * viscosity * vol_top;

      if (v_state(i - 1, j + 1, k) == FLUID)
        matrix.add_to_element(index, v_ind(i - 1, j + 1, k, nx, ny, nz),
                              factor * viscosity * vol_top);
      else if (v_state(i - 1, j + 1, k) == SOLID)
        rhs[index] -= v(i - 1, j + 1, k) * factor * viscosity * vol_top;
    }

    // v_x_bottom
    if (vol_bottom > 0.) {
      if (v_state(i, j, k) == FLUID)
        matrix.add_to_element(index, v_ind(i, j, k, nx, ny, nz),
                              +factor * viscosity * vol_bottom);
      else if (v_state(i, j, k) == SOLID)
        rhs[index] -= v(i, j, k) * factor * viscosity * vol_bottom;

      if (v_state(i - 1, j, k) == FLUID)
        matrix.add_to_element(index, v_ind(i - 1, j, k, nx, ny, nz),
                              -factor * viscosity * vol_bottom);
      else if (v_state(i - 1, j, k) == SOLID)
        rhs[index] -= -v(i - 1, j, k) * factor * viscosity * vol_bottom;
    }

    // w_x_front
    if (vol_front > 0.) {
      if (w_state(i, j, k + 1) == FLUID)
        matrix.add_to_element(index, w_ind(i, j, k + 1, nx, ny, nz),
                              -factor * viscosity * vol_front);
      else if (w_state(i, j, k + 1) == SOLID)
        rhs[index] -= -w(i, j, k + 1) * factor * viscosity * vol_front;

      if (w_state(i - 1, j, k + 1) == FLUID)
        matrix.add_to_element(index, w_ind(i - 1, j, k + 1, nx, ny, nz),
                              factor * viscosity * vol_front);
      else if (w_state(i - 1, j, k + 1) == SOLID)
        rhs[index] -= w(i - 1, j, k + 1) * factor * viscosity * vol_front;
    }

    // w_x_back
    if (vol_back > 0.) {
      if (w_state(i, j, k) == FLUID)
        matrix.add_to_element(index, w_ind(i, j, k, nx, ny, nz),
                              +factor * viscosity * vol_back);
      else if (w_state(i, j, k) == SOLID)
        rhs[index] -= w(i, j, k) * factor * viscosity * vol_back;

      if (w_state(i - 1, j, k) == FLUID)
        matrix.add_to_element(index, w_ind(i - 1, j, k, nx, ny, nz),
                              -factor * viscosity * vol_back);
      else if (w_state(i - 1, j, k) == SOLID)
        rhs[index] -= -w(i - 1, j, k) * factor * viscosity * vol_back;
    }
  });

  threadutils::thread_pool::ParallelFor(v_base_idx, w_base_idx, [&](int index) {
    const int i = dof_ijk[index](0);
    const int j = dof_ijk[index](1);
    const int k = dof_ijk[index](2) - nz;

    rhs[index] = vol_v(i, j, k) * v(i, j, k);
    matrix.set_element(index, index, vol_v(i, j, k));

    scalar vol_right = vol_ez(i + 1, j, k);
    scalar vol_left = vol_ez(i, j, k);

    scalar vol_top = vol_c(i, j, k);
    scalar vol_bottom = vol_c(i, j - 1, k);

    scalar vol_front = vol_ex(i, j, k + 1);
    scalar vol_back = vol_ex(i, j, k);

    // v_x_right
    if (vol_right > 0.0) {
      matrix.add_to_element(index, index, +factor * viscosity * vol_right);
      if (v_state(i + 1, j, k) == FLUID)
        matrix.add_to_element(index, v_ind(i + 1, j, k, nx, ny, nz),
                              -factor * viscosity * vol_right);
      else if (v_state(i + 1, j, k) == SOLID)
        rhs[index] -= -v(i + 1, j, k) * factor * viscosity * vol_right;
    }

    // v_x_left
    if (vol_left > 0.0) {
      matrix.add_to_element(index, index, +factor * viscosity * vol_left);
      if (v_state(i - 1, j, k) == FLUID)
        matrix.add_to_element(index, v_ind(i - 1, j, k, nx, ny, nz),
                              -factor * viscosity * vol_left);
      else if (v_state(i - 1, j, k) == SOLID)
        rhs[index] -= -v(i - 1, j, k) * factor * viscosity * vol_left;
    }

    // vy_top
    if (vol_top > 0.0) {
      matrix.add_to_element(index, index, +2 * factor * viscosity * vol_top);
      if (v_state(i, j + 1, k) == FLUID)
        matrix.add_to_element(index, v_ind(i, j + 1, k, nx, ny, nz),
                              -2 * factor * viscosity * vol_top);
      else if (v_state(i, j + 1, k) == SOLID)
        rhs[index] -= -2 * factor * viscosity * vol_top * v(i, j + 1, k);
    }

    // vy_bottom
    if (vol_bottom > 0.0) {
      matrix.add_to_element(index, index, +2 * factor * viscosity * vol_bottom);
      if (v_state(i, j - 1, k) == FLUID)
        matrix.add_to_element(index, v_ind(i, j - 1, k, nx, ny, nz),
                              -2 * factor * viscosity * vol_bottom);
      else if (v_state(i, j - 1, k) == SOLID)
        rhs[index] -= -2 * factor * viscosity * vol_bottom * v(i, j - 1, k);
    }

    // v_z_front
    if (vol_front > 0.0) {
      matrix.add_to_element(index, index, +factor * viscosity * vol_front);
      if (v_state(i, j, k + 1) == FLUID)
        matrix.add_to_element(index, v_ind(i, j, k + 1, nx, ny, nz),
                              -factor * viscosity * vol_front);
      else if (v_state(i, j, k + 1) == SOLID)
        rhs[index] -= -v(i, j, k + 1) * factor * viscosity * vol_front;
    }

    // v_z_back
    if (vol_back > 0.0) {
      matrix.add_to_element(index, index, +factor * viscosity * vol_back);
      if (v_state(i, j, k - 1) == FLUID)
        matrix.add_to_element(index, v_ind(i, j, k - 1, nx, ny, nz),
                              -factor * viscosity * vol_back);
      else if (v_state(i, j, k - 1) == SOLID)
        rhs[index] -= -v(i, j, k - 1) * factor * viscosity * vol_back;
    }

    // u_y_right
    if (vol_right > 0.0) {
      if (u_state(i + 1, j, k) == FLUID)
        matrix.add_to_element(index, u_ind(i + 1, j, k, nx, ny),
                              -factor * viscosity * vol_right);
      else if (u_state(i + 1, j, k) == SOLID)
        rhs[index] -= -u(i + 1, j, k) * factor * viscosity * vol_right;

      if (u_state(i + 1, j - 1, k) == FLUID)
        matrix.add_to_element(index, u_ind(i + 1, j - 1, k, nx, ny),
                              factor * viscosity * vol_right);
      else if (u_state(i + 1, j - 1, k) == SOLID)
        rhs[index] -= u(i + 1, j - 1, k) * factor * viscosity * vol_right;
    }

    // u_y_left
    if (vol_left > 0.0) {
      if (u_state(i, j, k) == FLUID)
        matrix.add_to_element(index, u_ind(i, j, k, nx, ny),
                              factor * viscosity * vol_left);
      else if (u_state(i, j, k) == SOLID)
        rhs[index] -= u(i, j, k) * factor * viscosity * vol_left;

      if (u_state(i, j - 1, k) == FLUID)
        matrix.add_to_element(index, u_ind(i, j - 1, k, nx, ny),
                              -factor * viscosity * vol_left);
      else if (u_state(i, j - 1, k) == SOLID)
        rhs[index] -= -u(i, j - 1, k) * factor * viscosity * vol_left;
    }

    // w_y_front
    if (vol_front > 0.0) {
      if (w_state(i, j, k + 1) == FLUID)
        matrix.add_to_element(index, w_ind(i, j, k + 1, nx, ny, nz),
                              -factor * viscosity * vol_front);
      else if (w_state(i, j, k + 1) == SOLID)
        rhs[index] -= -w(i, j, k + 1) * factor * viscosity * vol_front;

      if (w_state(i, j - 1, k + 1) == FLUID)
        matrix.add_to_element(index, w_ind(i, j - 1, k + 1, nx, ny, nz),
                              factor * viscosity * vol_front);
      else if (w_state(i, j - 1, k + 1) == SOLID)
        rhs[index] -= w(i, j - 1, k + 1) * factor * viscosity * vol_front;
    }

    // w_y_back
    if (vol_back > 0.0) {
      if (w_state(i, j, k) == FLUID)
        matrix.add_to_element(index, w_ind(i, j, k, nx, ny, nz),
                              factor * viscosity * vol_back);
      else if (w_state(i, j, k) == SOLID)
        rhs[index] -= w(i, j, k) * factor * viscosity * vol_back;

      if (w_state(i, j - 1, k) == FLUID)
        matrix.add_to_element(index, w_ind(i, j - 1, k, nx, ny, nz),
                              -factor * viscosity * vol_back);
      else if (w_state(i, j - 1, k) == SOLID)
        rhs[index] -= -w(i, j - 1, k) * factor * viscosity * vol_back;
    }
  });

  threadutils::thread_pool::ParallelFor(
      w_base_idx, total_ijk_size, [&](int index) {
        const int i = dof_ijk[index](0);
        const int j = dof_ijk[index](1);
        const int k = dof_ijk[index](2) - (nz * 2);

        rhs[index] = vol_w(i, j, k) * w(i, j, k);
        matrix.set_element(index, index, vol_w(i, j, k));

        scalar vol_right = vol_ey(i + 1, j, k);
        scalar vol_left = vol_ey(i, j, k);

        scalar vol_top = vol_ex(i, j + 1, k);
        scalar vol_bottom = vol_ex(i, j, k);

        scalar vol_front = vol_c(i, j, k);
        scalar vol_back = vol_c(i, j, k - 1);

        // w_x_right
        if (vol_right > 0.) {
          matrix.add_to_element(index, index, +factor * viscosity * vol_right);
          if (w_state(i + 1, j, k) == FLUID)
            matrix.add_to_element(index, w_ind(i + 1, j, k, nx, ny, nz),
                                  -factor * viscosity * vol_right);
          else if (w_state(i + 1, j, k) == SOLID)
            rhs[index] -= -factor * viscosity * vol_right * w(i + 1, j, k);
        }

        // w_x_left
        if (vol_left > 0.) {
          matrix.add_to_element(index, index, factor * viscosity * vol_left);
          if (w_state(i - 1, j, k) == FLUID)
            matrix.add_to_element(index, w_ind(i - 1, j, k, nx, ny, nz),
                                  -factor * viscosity * vol_left);
          else if (w_state(i - 1, j, k) == SOLID)
            rhs[index] -= -factor * viscosity * vol_left * w(i - 1, j, k);
        }

        // w_y_top
        if (vol_top > 0.) {
          matrix.add_to_element(index, index, +factor * viscosity * vol_top);
          if (w_state(i, j + 1, k) == FLUID)
            matrix.add_to_element(index, w_ind(i, j + 1, k, nx, ny, nz),
                                  -factor * viscosity * vol_top);
          else if (w_state(i, j + 1, k) == SOLID)
            rhs[index] -= -factor * viscosity * vol_top * w(i, j + 1, k);
        }

        // w_y_bottom
        if (vol_bottom > 0.) {
          matrix.add_to_element(index, index, factor * viscosity * vol_bottom);
          if (w_state(i, j - 1, k) == FLUID)
            matrix.add_to_element(index, w_ind(i, j - 1, k, nx, ny, nz),
                                  -factor * viscosity * vol_bottom);
          else if (w_state(i, j - 1, k) == SOLID)
            rhs[index] -= -factor * viscosity * vol_bottom * w(i, j - 1, k);
        }

        // w_z_front
        if (vol_front > 0.) {
          matrix.add_to_element(index, index,
                                +2 * factor * viscosity * vol_front);
          if (w_state(i, j, k + 1) == FLUID)
            matrix.add_to_element(index, w_ind(i, j, k + 1, nx, ny, nz),
                                  -2 * factor * viscosity * vol_front);
          else if (w_state(i, j, k + 1) == SOLID)
            rhs[index] -= -2 * factor * viscosity * vol_front * w(i, j, k + 1);
        }

        // w_z_back
        if (vol_back > 0.) {
          matrix.add_to_element(index, index,
                                +2 * factor * viscosity * vol_back);
          if (w_state(i, j, k - 1) == FLUID)
            matrix.add_to_element(index, w_ind(i, j, k - 1, nx, ny, nz),
                                  -2 * factor * viscosity * vol_back);
          else if (w_state(i, j, k - 1) == SOLID)
            rhs[index] -= -2 * factor * viscosity * vol_back * w(i, j, k - 1);
        }

        // u_z_right
        if (vol_right > 0.) {
          if (u_state(i + 1, j, k) == FLUID)
            matrix.add_to_element(index, u_ind(i + 1, j, k, nx, ny),
                                  -factor * viscosity * vol_right);
          else if (u_state(i + 1, j, k) == SOLID)
            rhs[index] -= -u(i + 1, j, k) * factor * viscosity * vol_right;

          if (u_state(i + 1, j, k - 1) == FLUID)
            matrix.add_to_element(index, u_ind(i + 1, j, k - 1, nx, ny),
                                  factor * viscosity * vol_right);
          else if (u_state(i + 1, j, k - 1) == SOLID)
            rhs[index] -= u(i + 1, j, k - 1) * factor * viscosity * vol_right;
        }

        // u_z_left
        if (vol_left > 0.) {
          if (u_state(i, j, k) == FLUID)
            matrix.add_to_element(index, u_ind(i, j, k, nx, ny),
                                  factor * viscosity * vol_left);
          else if (u_state(i, j, k) == SOLID)
            rhs[index] -= u(i, j, k) * factor * viscosity * vol_left;

          if (u_state(i, j, k - 1) == FLUID)
            matrix.add_to_element(index, u_ind(i, j, k - 1, nx, ny),
                                  -factor * viscosity * vol_left);
          else if (u_state(i, j, k - 1) == SOLID)
            rhs[index] -= -u(i, j, k - 1) * factor * viscosity * vol_left;
        }

        // v_z_top
        if (vol_top > 0.) {
          if (v_state(i, j + 1, k) == FLUID)
            matrix.add_to_element(index, v_ind(i, j + 1, k, nx, ny, nz),
                                  -factor * viscosity * vol_top);
          else if (v_state(i, j + 1, k) == SOLID)
            rhs[index] -= -v(i, j + 1, k) * factor * viscosity * vol_top;

          if (v_state(i, j + 1, k - 1) == FLUID)
            matrix.add_to_element(index, v_ind(i, j + 1, k - 1, nx, ny, nz),
                                  factor * viscosity * vol_top);
          else if (v_state(i, j + 1, k - 1) == SOLID)
            rhs[index] -= v(i, j + 1, k - 1) * factor * viscosity * vol_top;
        }

        // v_z_bottom
        if (vol_bottom > 0.) {
          if (v_state(i, j, k) == FLUID)
            matrix.add_to_element(index, v_ind(i, j, k, nx, ny, nz),
                                  +factor * viscosity * vol_bottom);
          else if (v_state(i, j, k) == SOLID)
            rhs[index] -= v(i, j, k) * factor * viscosity * vol_bottom;

          if (v_state(i, j, k - 1) == FLUID)
            matrix.add_to_element(index, v_ind(i, j, k - 1, nx, ny, nz),
                                  -factor * viscosity * vol_bottom);
          else if (v_state(i, j, k - 1) == SOLID)
            rhs[index] -= -v(i, j, k - 1) * factor * viscosity * vol_bottom;
        }
      });

  PCGSolver<double> solver;
  double res_out;
  int iter_out;
  solver.set_solver_parameters(1e-9, 10000, 0.97, 0.1);
  bool success = false;

  success = solver.solve(matrix, rhs, soln, res_out, iter_out);

  if (iter_out >= 1000) {
    std::cerr << "\n\n\n**********FAILED**************\n\n\n" << std::endl;
    exit(0);
  }

  int parallel_size = (nx + 1) * (ny + 1) * (nz + 1);
  int slice = (nx + 1) * (ny + 1);
  tbb::parallel_for(0, parallel_size, 1, [&](int thread_idx) {
    int k = thread_idx / slice;
    int j = (thread_idx % slice) / (nx + 1);
    int i = thread_idx % (nx + 1);

    if (i >= 0 && i < nx + 1 && j >= 0 && j < ny && k >= 0 && k < nz) {
      if (u_state(i, j, k) == FLUID) {
        scalar new_u = (scalar)soln[u_ind(i, j, k, nx, ny)];
        u_visc_impulse(i, j, k) = (new_u - u(i, j, k)) / dt;
        u(i, j, k) = new_u;
      } else {
        u_visc_impulse(i, j, k) = 0.0;
        u(i, j, k) = 0;
      }
    }

    if (i >= 0 && i < nx && j >= 0 && j < ny + 1 && k >= 0 && k < nz) {
      if (v_state(i, j, k) == FLUID) {
        scalar new_v = (scalar)soln[v_ind(i, j, k, nx, ny, nz)];
        v_visc_impulse(i, j, k) = (new_v - v(i, j, k)) / dt;
        v(i, j, k) = new_v;
      } else {
        v_visc_impulse(i, j, k) = 0.0;
        v(i, j, k) = 0;
      }
    }

    if (i >= 0 && i < nx && j >= 0 && j < ny && k >= 0 && k < nz + 1) {
      if (w_state(i, j, k) == FLUID) {
        scalar new_w = (scalar)soln[w_ind(i, j, k, nx, ny, nz)];
        w_visc_impulse(i, j, k) = (new_w - w(i, j, k)) / dt;
        w(i, j, k) = new_w;
      } else {
        w_visc_impulse(i, j, k) = 0.0;
        w(i, j, k) = 0;
      }
    }
  });
}
}  // namespace libwethair
