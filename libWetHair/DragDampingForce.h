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

#ifndef LIBWETHAIR_CORE_DRAG_DAMPING_FORCE_H_
#define LIBWETHAIR_CORE_DRAG_DAMPING_FORCE_H_

#include <Eigen/Core>
#include <iostream>

#include "Force.h"
#include "MathDefs.h"

namespace libwethair {

template <int DIM>
class TwoDScene;

template <int DIM>
class DragDampingForce : public Force {
 public:
  DragDampingForce(const TwoDScene<DIM>& scene, const scalar& b, int hidx);

  virtual ~DragDampingForce();

  virtual void computeIntegrationVars(const VectorXs& x, const VectorXs& v,
                                      const VectorXs& m, VectorXs& lambda,
                                      VectorXs& lambda_v, TripletXs& J,
                                      TripletXs& Jv, TripletXs& Jxv,
                                      TripletXs& tildeK, TripletXs& stiffness,
                                      TripletXs& damping, VectorXs& Phi,
                                      VectorXs& Phiv, const scalar& dt);

  virtual int numConstraintPos();

  virtual int numConstraintVel();

  virtual int numJ();

  virtual int numJv();

  virtual int numJxv();

  virtual int numTildeK();

  virtual bool isParallelized();

  virtual bool isPrecomputationParallelized();

  virtual void storeLambda(const VectorXs& lambda, const VectorXs& lambda_v);

  virtual Force* createNewCopy();

  virtual void getAffectedVars(int pidx, std::unordered_set<int>& vars);

  virtual int getAffectedHair(const std::vector<int> particle_to_hairs);

  virtual bool isContained(int pidx);

  virtual const char* name();

 private:
  scalar m_b;
  int m_dofs;
  int m_hidx;
  VectorXs m_lambda_v;
  const TwoDScene<DIM>& m_scene;
};

}  // namespace libwethair

#endif  // LIBWETHAIR_CORE_DRAG_DAMPING_FORCE_H_
