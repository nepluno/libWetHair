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

#ifndef LIBWETHAIR_CORE_COMPLIANT_IMPLICIT_EULER_H_
#define LIBWETHAIR_CORE_COMPLIANT_IMPLICIT_EULER_H_

#include <Eigen/Dense>
#include <iostream>

#include "MathUtilities.h"
#include "SceneStepper.h"

namespace libwethair {

template <int DIM>
class TwoDScene;

template <int DIM>
class CompliantImplicitEuler : public SceneStepper<DIM> {
 public:
  CompliantImplicitEuler(TwoDScene<DIM>* scene, int max_iters, scalar criterion,
                         bool autoUpdateNextInit = true);

  virtual ~CompliantImplicitEuler();

  virtual bool stepScene(TwoDScene<DIM>& scene, scalar dt,
                         bool updatePreCompute = true);

  virtual std::string getName() const;

 private:
  void zeroFixedDoFs(const TwoDScene<DIM>& scene, VectorXs& vec);
  void updateNumConstraints(const VectorXs& dx, const VectorXs& dv, scalar dt);

  bool m_bAutoUpdateNextInit;
  int m_max_iters;
  scalar m_criterion;

  TwoDScene<DIM>* m_scene;

  TripletXs m_Kext_nz;
  TripletXs m_A_nz;
  TripletXs m_J_nz;
  TripletXs m_Jv_nz;
  TripletXs m_Jxv_nz;
  TripletXs m_Fix_nz;
  TripletXs m_M_nz;
  TripletXs m_invC_nz;
  TripletXs m_invCv_nz;

  SparseXs m_Kext;
  SparseXs m_A;
  SparseXs m_J;
  SparseXs m_JC;
  SparseXs m_Jv;
  SparseXs m_JvC;
  SparseXs m_Jxv;
  SparseXs m_Fix;
  SparseXs m_M;
  SparseXs m_invC;
  SparseXs m_invCv;

  VectorXs m_lambda;
  VectorXs m_lambda_v;
  VectorXs m_gradU;
  VectorXs m_Phi;
  VectorXs m_Phi_v;
  VectorXs m_b;
  VectorXs m_vplus;

  Vector6i m_interhair_idx;
  Vector6i m_interhair_num;

  Eigen::SimplicialLDLT<SparseXs> m_solver;
  Eigen::ConjugateGradient<SparseXs> m_iterative_solver;
};

}  // namespace libwethair

#endif  // LIBWETHAIR_CORE_COMPLIANT_IMPLICIT_EULER_H_
