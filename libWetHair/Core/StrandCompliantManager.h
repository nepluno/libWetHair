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

#ifndef LIBWETHAIR_CORE_STRAND_COMPLIANT_MANAGER_
#define LIBWETHAIR_CORE_STRAND_COMPLIANT_MANAGER_

#include <future>
#include <vector>

#include "MathDefs.h"
#include "SceneStepper.h"
#include "TwoDScene.h"

template <int DIM>
class StrandCompliantEuler;

template <int DIM>
class StrandCompliantManager : public SceneStepper<DIM> {
 public:
  StrandCompliantManager(TwoDScene<DIM>* scene, int max_newton, int max_iters,
                         scalar criterion, bool compute_interhair,
                         bool use_preconditioner);

  virtual ~StrandCompliantManager();

  virtual bool stepScene(TwoDScene<DIM>& scene, scalar dt,
                         bool updatePreCompute = true);

  virtual bool stepSceneNonlinear(TwoDScene<DIM>& scene, scalar dt,
                                  bool updatePreCompute = true);

  virtual bool stepSceneLinear(TwoDScene<DIM>& scene, scalar dt,
                               bool updatePreCompute = true);

  virtual std::string getName() const;

  virtual void zeroFixedDoFs(const TwoDScene<DIM>& scene, VectorXs& vec);

  virtual void localPreIterate(TwoDScene<DIM>& scene, scalar dt);

  virtual void localStepScene(TwoDScene<DIM>& scene, scalar dt, VectorXs& r,
                              const VectorXs& b);

  virtual void localPreconditionScene(TwoDScene<DIM>& scene, scalar dt,
                                      VectorXs& r, const VectorXs& b);

  virtual void localUpdateLambda(TwoDScene<DIM>& scene, const VectorXs& dx,
                                 const VectorXs& dv, scalar dt);

  virtual void updateNextV(TwoDScene<DIM>& scene, const VectorXs& vplus);

  virtual void localUpdateNumConstraints(const VectorXs& dx, const VectorXs& dv,
                                         scalar dt);

  virtual void interhairUpdateNumConstraints(const VectorXs& dx,
                                             const VectorXs& dv, scalar dt);

  virtual void computeRHSLocal(TwoDScene<DIM>& scene, scalar dt, VectorXs& b);

  virtual void computeRHSIncrementalLocal(TwoDScene<DIM>& scene, scalar dt,
                                          VectorXs& b, const VectorXs& vplus);

  virtual void computeRHSIncrementalInterhair(TwoDScene<DIM>& scene, scalar dt,
                                              VectorXs& b,
                                              const VectorXs& vplus);

  virtual void computeAp(TwoDScene<DIM>& scene, const VectorXs& p, VectorXs& b,
                         scalar dt);

  virtual void write(std::vector<scalar>&) const;

  virtual void read(const scalar* data);

  virtual size_t size();

 protected:
  std::vector<StrandCompliantEuler<DIM>*> m_integrators;

  scalar m_dt;
  TwoDScene<DIM>* m_scene;

  TripletXs m_A_nz;
  TripletXs m_J_nz;
  TripletXs m_Jv_nz;
  TripletXs m_Jxv_nz;
  TripletXs m_Fix_nz;
  TripletXs m_M_nz;
  TripletXs m_invC_nz;
  TripletXs m_invCv_nz;

  TripletXs m_J_interhair_nz;
  TripletXs m_Jv_interhair_nz;
  SparseXs m_J_interhair;
  SparseXs m_Jv_interhair;
  SparseXs m_JT_interhair;
  SparseXs m_JvT_interhair;
  VectorXs m_invC_interhair;
  VectorXs m_invCv_interhair;

  VectorXs m_lambda;
  VectorXs m_lambda_v;
  VectorXs m_gradU;
  VectorXs m_Phi;
  VectorXs m_Phi_interhair;
  VectorXs m_Phi_v;
  VectorXs m_Phiv_interhair;
  VectorXs m_vplus;

  VectorXs m_dvi;
  VectorXs m_dxi;
  VectorXs m_dv;
  VectorXs m_dx;
  VectorXs m_dx_scripted;

  VectorXs m_v;
  VectorXs m_r;
  VectorXs m_p;
  VectorXs m_q;
  VectorXs m_t;
  VectorXs m_z;
  VectorXs m_rhs;

  VectorXs m_cv_buffer;
  VectorXs m_c_buffer;

  Vector6i m_interhair_idx;
  Vector6i m_interhair_num;

  int m_max_num_newton;
  int m_max_num_iters;
  scalar m_criterion;

  bool m_compute_interhair;
  bool m_use_preconditioner;

  friend class StrandCompliantEuler<DIM>;
};

#endif  // LIBWETHAIR_CORE_STRAND_COMPLIANT_MANAGER_
