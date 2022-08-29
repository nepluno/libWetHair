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

#ifndef LIBWETHAIR_CORE_LINEAR_SPRING_FORCE_H_
#define LIBWETHAIR_CORE_LINEAR_SPRING_FORCE_H_

#include <Eigen/Core>
#include <iostream>

#include "Force.h"

namespace libwethair {

template <int DIM>
class TwoDScene;

template <int DIM>
class LinearSpringForce : public Force {
 public:
  LinearSpringForce(TwoDScene<DIM>* parent,
                    const std::pair<int, int>& endpoints, const scalar& k,
                    const scalar& l0, const scalar& b = 0.0);

  virtual ~LinearSpringForce();

  virtual Force* createNewCopy();

  virtual void preCompute(const VectorXs& x, const VectorXs& v,
                          const VectorXs& m, const scalar& dt);

  virtual void getAffectedVars(int pidx, std::unordered_set<int>& vars);

  virtual int getAffectedHair(const std::vector<int> particle_to_hairs);

  virtual bool isContained(int pidx);

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

  virtual const char* name();

 private:
  std::pair<int, int> m_endpoints;
  scalar m_k;
  scalar m_l0;
  scalar m_b;

  Vectors<DIM> m_D;
  Vectors<DIM> m_lambda_v;
  Vectors<DIM> m_lambda;

  TwoDScene<DIM>* m_scene;
};

}  // namespace libwethair

#endif  // LIBWETHAIR_CORE_LINEAR_SPRING_FORCE_H_
