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

#ifndef LIBWETHAIR_CORE_SPRING_FORCE_H_
#define LIBWETHAIR_CORE_SPRING_FORCE_H_

#include <Eigen/Core>
#include <iostream>

#include "Force.h"
#include "TwoDScene.h"

template <int DIM>
class SpringForce : public Force {
 public:
  SpringForce(const std::pair<int, int>& endpoints, const scalar& k,
              const scalar& l0, TwoDScene<DIM>* scene, const scalar& b = 0.0);

  virtual ~SpringForce();

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

  virtual Force* createNewCopy();

  virtual void getAffectedVars(int pidx, std::unordered_set<int>& vars);

  virtual int getAffectedHair(const std::vector<int> particle_to_hairs);

  virtual bool isContained(int pidx);

  virtual void storeLambda(const VectorXs& lambda, const VectorXs& lambda_v);

  virtual const char* name();

 private:
  std::pair<int, int> m_endpoints;
  scalar m_k;
  scalar m_l0;
  scalar m_b;

  TwoDScene<DIM>* m_scene;

  scalar m_lambda;
  scalar m_lambda_v;
};

#endif  // LIBWETHAIR_CORE_SPRING_FORCE_H_
