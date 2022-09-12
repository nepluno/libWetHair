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

#ifndef LIBWETHAIR_CORE_STRAND_COMPLIANT_EULER_
#define LIBWETHAIR_CORE_STRAND_COMPLIANT_EULER_

#include <iostream>
//#include <ceres/ceres.h>

#include "MathUtilities.h"
#include "SceneStepper.h"

namespace libwethair {

template <int DIM>
class StrandCompliantManager;

template <int DIM>
class StrandCompliantEuler {
 public:
  StrandCompliantEuler(StrandCompliantManager<DIM>* parent, int hidx);

  virtual ~StrandCompliantEuler();

  virtual bool stepScene(TwoDScene<DIM>& scene, scalar dt);

  virtual bool stepScene(TwoDScene<DIM>& scene, scalar dt, VectorXs& r,
                         const VectorXs& b);

  virtual bool PreconditionScene(TwoDScene<DIM>& scene, scalar dt, VectorXs& r,
                                 const VectorXs& b);

  virtual void preIterate(TwoDScene<DIM>& scene, scalar dt);

  virtual void computeRHS(TwoDScene<DIM>& scene, scalar dt, VectorXs& b);

  virtual void computeRHSIncremental(TwoDScene<DIM>& scene, scalar dt,
                                     VectorXs& b, const VectorXs& vplus);

  virtual void computeAp(const VectorXs& p, VectorXs& b);

  virtual void updateLambda(TwoDScene<DIM>& scene, const VectorXs& dx,
                            const VectorXs& dv, scalar dt);

  virtual void updateNextV(TwoDScene<DIM>& scene, const VectorXs& vplus);

 private:
  TripletXs m_A_nz;
  TripletXs m_J_nz;
  TripletXs m_Jv_nz;
  TripletXs m_Jxv_nz;
  TripletXs m_invC_nz;
  TripletXs m_invCv_nz;

  TripletXs m_J_inter_nz;
  TripletXs m_Jv_inter_nz;
  TripletXs m_invC_inter_nz;
  TripletXs m_invCv_inter_nz;

  SparseXs m_A;
  SparseXs m_J;
  SparseXs m_Jv;
  SparseXs m_Jxv;
  SparseXs m_invC;
  SparseXs m_invCv;

  SparseXs m_J_inter;
  SparseXs m_Jv_inter;
  SparseXs m_invC_inter;
  SparseXs m_invCv_inter;

  VectorXs m_A_inv_diag;

  Eigen::SimplicialLDLT<SparseXs> m_solver;

  int m_hidx;

  StrandCompliantManager<DIM>* m_parent;

  int m_start_global_dof;
  int m_num_global_dof;
};

}  // namespace libwethair

#endif  // LIBWETHAIR_CORE_STRAND_COMPLIANT_EULER_
