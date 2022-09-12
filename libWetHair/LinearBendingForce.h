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

#ifndef LIBWETHAIR_CORE_LINEAR_BENDING_FORCE_H_
#define LIBWETHAIR_CORE_LINEAR_BENDING_FORCE_H_

#include <Eigen/Core>
#include <iostream>

#include "Force.h"

namespace libwethair {

template <int DIM>
class TwoDScene;

template <int DIM>
class LinearBendingForce : public Force {
 public:
  LinearBendingForce(TwoDScene<DIM>* parent, int idx1, int idx2, int idx3,
                     const scalar& alpha, const scalar& beta,
                     const Vectors<DIM - 1>& theta0, const scalar& eb1n,
                     const scalar& eb2n);

  virtual ~LinearBendingForce();

  virtual Force* createNewCopy();

  virtual void preCompute(const VectorXs& x, const VectorXs& v,
                          const VectorXs& m, const scalar& dt);

  virtual void getAffectedVars(int pidx, std::unordered_set<int>& vars);

  virtual int getAffectedHair(const std::vector<int> particle_to_hairs);

  virtual bool isContained(int pidx);

  virtual const char* name();

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

 private:
  int m_idx1;
  int m_idx2;
  int m_idx3;

  Matrixs<DIM> m_R;  // rotation matrix

  scalar m_alpha;             // stiffness coefficient
  scalar m_beta;              // damping coefficient
  Vectors<DIM - 1> m_theta0;  // rest angle
  scalar m_eb1n;              // norm of e1 bar
  scalar m_eb2n;              // norm of e2 bar

  Vectors<DIM> m_x1;
  Vectors<DIM> m_x2;
  Vectors<DIM> m_x3;

  Vectors<DIM> m_L0;

  Vectors<DIM> m_RL0;

  Vectors<DIM> m_lambda_v;
  Vectors<DIM> m_lambda;

  scalar m_c1;
  scalar m_c2;

  TwoDScene<DIM>* m_scene;
};

}  // namespace libwethair

#endif  // LIBWETHAIR_CORE_LINEAR_BENDING_FORCE_H_
