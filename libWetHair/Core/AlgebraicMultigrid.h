//
// This file is part of the libWetHair open source project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright 2015 Xinxin Zhang
//

#ifndef LIBWETHAIR_CORE_ALGEBRAICMULTIGRID_H_
#define LIBWETHAIR_CORE_ALGEBRAICMULTIGRID_H_

#include <tbb/tbb.h>

#include <cmath>
#include <iostream>
#include <memory>
#include <thread>
#include <vector>

#include "GeometricLevelGen.h"
#include "MathDefs.h"
#include "pcgsolver/blas_wrapper.h"
#include "pcgsolver/sparse_matrix.h"
/*
given A_L, R_L, P_L, b,compute x using
Multigrid Cycles.

*/
//#define AMG_VERBOSE

template <class T>
void RBGS(const robertbridson::FixedSparseMatrix<T> &A, const std::vector<T> &b, std::vector<T> &x,
          int ni, int nj, int nk, int iternum) {
  for (int iter = 0; iter < iternum; iter++) {
    size_t num = ni * nj * nk;
    size_t slice = ni * nj;
    tbb::parallel_for((size_t)0, num, (size_t)1, [&](size_t thread_idx) {
      int k = thread_idx / slice;
      int j = (thread_idx % slice) / ni;
      int i = thread_idx % ni;
      if (k < nk && j < nj && i < ni) {
        if ((i + j + k) % 2 == 1) {
          unsigned int index = (unsigned int)thread_idx;
          T sum = 0;
          T diag = 0;
          for (int ii = A.rowstart[index]; ii < A.rowstart[index + 1]; ii++) {
            if (A.colindex[ii] != index)  // none diagonal terms
            {
              sum += A.value[ii] * x[A.colindex[ii]];
            } else  // record diagonal value A(i,i)
            {
              diag = A.value[ii];
            }
          }  // A(i,:)*x for off-diag terms
          if (diag != 0) {
            x[index] = (b[index] - sum) / diag;
          } else {
            x[index] = 0;
          }
        }
      }
    });

    tbb::parallel_for((size_t)0, num, (size_t)1, [&](size_t thread_idx) {
      int k = thread_idx / slice;
      int j = (thread_idx % slice) / ni;
      int i = thread_idx % ni;
      if (k < nk && j < nj && i < ni) {
        if ((i + j + k) % 2 == 0) {
          unsigned int index = (unsigned int)thread_idx;
          T sum = 0;
          T diag = 0;
          for (int ii = A.rowstart[index]; ii < A.rowstart[index + 1]; ii++) {
            if (A.colindex[ii] != index)  // none diagonal terms
            {
              sum += A.value[ii] * x[A.colindex[ii]];
            } else  // record diagonal value A(i,i)
            {
              diag = A.value[ii];
            }
          }  // A(i,:)*x for off-diag terms
          if (diag != 0) {
            x[index] = (b[index] - sum) / diag;
          } else {
            x[index] = 0;
          }
        }
      }
    });
  }
}

template <class T>
void restriction(const robertbridson::FixedSparseMatrix<T> &R, const robertbridson::FixedSparseMatrix<T> &A,
                 const std::vector<T> &x, const std::vector<T> &b_curr,
                 std::vector<T> &b_next) {
  b_next.assign(b_next.size(), 0);
  std::vector<T> r = b_curr;
  multiply_and_subtract(A, x, r);
  multiply(R, r, b_next);
  r.resize(0);
}
template <class T>
void prolongatoin(const robertbridson::FixedSparseMatrix<T> &P, const std::vector<T> &x_curr,
                  std::vector<T> &x_next) {
  static std::vector<T> xx;
  xx.resize(x_next.size());
  xx.assign(xx.size(), 0);
  multiply(P, x_curr, xx);  // xx = P*x_curr;
  robertbridson::BLAS::add_scaled(1.0, xx, x_next);
  xx.resize(0);
}

template <class T>
void RBGS_with_pattern(const robertbridson::FixedSparseMatrix<T> &A, const std::vector<T> &b,
                       std::vector<T> &x, std::vector<bool> &pattern, int iternum) {
  for (int iter = 0; iter < iternum; iter++) {
    size_t num = x.size();

    tbb::parallel_for((size_t)0, num, (size_t)1, [&](size_t thread_idx) {
      if (thread_idx >= A.rowstart.size() - 1 || thread_idx >= pattern.size())
        return;
      if (pattern[thread_idx] == true) {
        T sum = 0;
        T diag = 0;
        for (int ii = A.rowstart[thread_idx]; ii < A.rowstart[thread_idx + 1];
             ii++) {
          if (ii >= A.colindex.size() || ii >= A.value.size()) continue;

          if (A.colindex[ii] != thread_idx &&
              A.colindex[ii] < num)  // none diagonal terms
          {
            sum += A.value[ii] * x[A.colindex[ii]];
          } else  // record diagonal value A(i,i)
          {
            diag = A.value[ii];
          }
        }  // A(i,:)*x for off-diag terms
        if (diag != 0) {
          x[thread_idx] = (b[thread_idx] - sum) / diag;
        } else {
          x[thread_idx] = 0;
        }
      }
    });

    tbb::parallel_for((size_t)0, num, (size_t)1, [&](size_t thread_idx) {
      if (thread_idx >= A.rowstart.size() - 1 || thread_idx >= pattern.size())
        return;
      if (pattern[thread_idx] == false) {
        T sum = 0;
        T diag = 0;
        for (int ii = A.rowstart[thread_idx]; ii < A.rowstart[thread_idx + 1];
             ii++) {
          if (ii >= A.colindex.size() || ii >= A.value.size()) continue;

          if (A.colindex[ii] != thread_idx &&
              A.colindex[ii] < num)  // none diagonal terms
          {
            sum += A.value[ii] * x[A.colindex[ii]];
          } else  // record diagonal value A(i,i)
          {
            diag = A.value[ii];
          }
        }  // A(i,:)*x for off-diag terms
        if (diag != 0) {
          x[thread_idx] = (b[thread_idx] - sum) / diag;
        } else {
          x[thread_idx] = 0;
        }
      }
    });
  }
}

template <class T>
void amgVCycleCompressed(std::vector<std::shared_ptr<robertbridson::FixedSparseMatrix<T>>> &A_L,
                         std::vector<robertbridson::FixedSparseMatrix<T>> &R_L,
                         std::vector<robertbridson::FixedSparseMatrix<T>> &P_L,
                         std::vector<std::vector<bool>> &p_L, std::vector<T> &x,
                         const std::vector<T> &b) {
  int total_level = A_L.size();
  std::vector<std::vector<T>> x_L;
  std::vector<std::vector<T>> b_L;
  x_L.resize(total_level);
  b_L.resize(total_level);
  b_L[0] = b;
  x_L[0] = x;
  for (int i = 1; i < total_level; i++) {
    int unknowns = A_L[i]->n;
    x_L[i].resize(unknowns);
    x_L[i].assign(x_L[i].size(), 0);
    b_L[i].resize(unknowns);
    b_L[i].assign(b_L[i].size(), 0);
  }

  for (int i = 0; i < total_level - 1; i++) {
#ifdef AMG_VERBOSE
    printf("level: %d, RBGS\n", i);
#endif
    RBGS_with_pattern(*(A_L[i]), b_L[i], x_L[i], (p_L[i]), 4);
#ifdef AMG_VERBOSE
    printf("level: %d, restriction\n", i);
#endif
    restriction((R_L[i]), *(A_L[i]), x_L[i], b_L[i], b_L[i + 1]);
  }
  int i = total_level - 1;
#ifdef AMG_VERBOSE
  printf("level: %d, top solve\n", i);
#endif
  RBGS_with_pattern(*(A_L[i]), b_L[i], x_L[i], (p_L[i]), 200);
  for (int i = total_level - 2; i >= 0; i--) {
#ifdef AMG_VERBOSE
    printf("level: %d, prolongation\n", i);
#endif
    prolongatoin((P_L[i]), x_L[i + 1], x_L[i]);
#ifdef AMG_VERBOSE
    printf("level: %d, RBGS\n", i);
#endif
    RBGS_with_pattern(*(A_L[i]), b_L[i], x_L[i], (p_L[i]), 4);
  }
  x = x_L[0];

  for (int i = 0; i < total_level; i++) {
    x_L[i].resize(0);
    x_L[i].shrink_to_fit();
    b_L[i].resize(0);
    b_L[i].shrink_to_fit();
  }
}
template <class T>
void amgPrecondCompressed(std::vector<std::shared_ptr<robertbridson::FixedSparseMatrix<T>>> &A_L,
                          std::vector<robertbridson::FixedSparseMatrix<T>> &R_L,
                          std::vector<robertbridson::FixedSparseMatrix<T>> &P_L,
                          std::vector<std::vector<bool>> &p_L, std::vector<T> &x,
                          const std::vector<T> &b) {
  // printf("preconditioning begin\n");
  x.resize(b.size());
  x.assign(x.size(), 0);
  amgVCycleCompressed(A_L, R_L, P_L, p_L, x, b);
  // printf("preconditioning finished\n");
}

template <class T>
bool AMGPCGSolveSparse(const robertbridson::SparseMatrix<T> &matrix, const std::vector<T> &rhs,
                       std::vector<T> &result, std::vector<Vector3i> &Dof_ijk,
                       T tolerance_factor, int max_iterations, T &residual_out,
                       int &iterations_out, int ni, int nj, int nk) {
  static std::shared_ptr<robertbridson::FixedSparseMatrix<T>> fixed_matrix =
      std::make_shared<robertbridson::FixedSparseMatrix<T>>();
  fixed_matrix->construct_from_matrix(matrix);
  static std::vector<std::shared_ptr<robertbridson::FixedSparseMatrix<T>>> A_L;
  static std::vector<robertbridson::FixedSparseMatrix<T>> R_L;
  static std::vector<robertbridson::FixedSparseMatrix<T>> P_L;
  static std::vector<std::vector<bool>> p_L;
  std::vector<T> m, z, s, r;
  int total_level;
  levelGen<T> amg_levelGen;
#ifdef AMG_VERBOSE
  std::cout << "[AMG: generate levels]" << std::endl;
#endif
  amg_levelGen.generateLevelsGalerkinCoarseningSparse(
      A_L, R_L, P_L, p_L, total_level, fixed_matrix, Dof_ijk, ni, nj, nk);

  unsigned int n = matrix.n;
  if (m.size() != n) {
    m.resize(n);
    s.resize(n);
    z.resize(n);
    r.resize(n);
  }
  mathutils::zero(result);
  r = rhs;
  residual_out = robertbridson::BLAS::abs_max(r);
  if (residual_out == 0) {
    iterations_out = 0;

    for (int i = 0; i < total_level; i++) {
      A_L[i]->clear();
    }
    for (int i = 0; i < total_level - 1; i++) {
      R_L[i].clear();
      P_L[i].clear();
    }

    return true;
  }
  double tol = tolerance_factor * residual_out;
#ifdef AMG_VERBOSE
  std::cout << "[AMG: preconditioning]" << std::endl;
#endif
  amgPrecondCompressed(A_L, R_L, P_L, p_L, z, r);
#ifdef AMG_VERBOSE
  std::cout << "[AMG: first precond done]" << std::endl;
#endif
  double rho = robertbridson::BLAS::dot(z, r);
  if (rho == 0 || rho != rho) {
    for (int i = 0; i < total_level; i++) {
      A_L[i]->clear();
    }
    for (int i = 0; i < total_level - 1; i++) {
      R_L[i].clear();
      P_L[i].clear();
    }
    iterations_out = 0;
    return false;
  }

  s = z;
#ifdef AMG_VERBOSE
  std::cout << "[AMG: iterative solve]" << std::endl;
#endif
  int iteration;
  for (iteration = 0; iteration < max_iterations; ++iteration) {
    multiply(*fixed_matrix, s, z);
    // printf("multiply done\n");
    double alpha = rho / robertbridson::BLAS::dot(s, z);
    // printf("%d,%d,%d,%d\n",s.size(),z.size(),r.size(),result.size());
    robertbridson::BLAS::add_scaled(alpha, s, result);
    robertbridson::BLAS::add_scaled(-alpha, z, r);
    residual_out = robertbridson::BLAS::abs_max(r);

    if (residual_out <= tol) {
      iterations_out = iteration + 1;

      for (int i = 0; i < total_level; i++) {
        A_L[i]->clear();
      }
      for (int i = 0; i < total_level - 1; i++) {
        R_L[i].clear();
        P_L[i].clear();
      }

      return true;
    }
#ifdef AMG_VERBOSE
    std::cout << "[AMG: iterative preconditioning]" << std::endl;
#endif
    amgPrecondCompressed(A_L, R_L, P_L, p_L, z, r);
#ifdef AMG_VERBOSE
    std::cout << "[AMG: second precond done]" << std::endl;
#endif
    double rho_new = robertbridson::BLAS::dot(z, r);
    double beta = rho_new / rho;
    robertbridson::BLAS::add_scaled(beta, s, z);
    s.swap(z);  // s=beta*s+z
    rho = rho_new;
  }
  iterations_out = iteration;
  for (int i = 0; i < total_level; i++) {
    A_L[i]->clear();
  }
  for (int i = 0; i < total_level - 1; i++) {
    R_L[i].clear();
    P_L[i].clear();
  }
  return false;
}

#endif  // LIBWETHAIR_CORE_ALGEBRAICMULTIGRID_H_
