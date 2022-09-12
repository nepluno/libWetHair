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

#ifndef LIBWETHAIR_CORE_SIMPLE_GRAVITY_FORCE_H_
#define LIBWETHAIR_CORE_SIMPLE_GRAVITY_FORCE_H_

#include <Eigen/Core>
#include <iostream>

#include "Force.h"

namespace libwethair {

template <int DIM>
class TwoDScene;

template <int DIM>
class SimpleGravityForce : public Force {
  VectorXs m_buoyancy;
  TwoDScene<DIM>* m_scene;

 public:
  SimpleGravityForce(const Vector3s& gravity, TwoDScene<DIM>* scene);

  virtual ~SimpleGravityForce();

  virtual void addEnergyToTotal(const VectorXs& x, const VectorXs& v,
                                const VectorXs& m, scalar& E);

  virtual void addGradEToTotal(const VectorXs& x, const VectorXs& v,
                               const VectorXs& m, VectorXs& gradE);

  virtual void addHessXToTotal(const VectorXs& x, const VectorXs& v,
                               const VectorXs& m, TripletXs& hessE);

  virtual void addHessVToTotal(const VectorXs& x, const VectorXs& v,
                               const VectorXs& m, TripletXs& hessE);

  virtual void addGradEToTotal(const VectorXs& x, const VectorXs& v,
                               const VectorXs& m, VectorXs& gradE, int pidx);

  virtual void addHessXToTotal(const VectorXs& x, const VectorXs& v,
                               const VectorXs& m, VectorXs& hessE, int pidx);

  virtual void addHessVToTotal(const VectorXs& x, const VectorXs& v,
                               const VectorXs& m, VectorXs& hessE, int pidx);

  virtual void preCompute(const VectorXs& x, const VectorXs& v,
                          const VectorXs& m, const scalar& dt);

  virtual int numConstraintPos();

  virtual int numConstraintVel();

  virtual int numJ();

  virtual int numJv();

  virtual int numJxv();

  virtual int numTildeK();

  virtual bool isParallelized();

  virtual bool isPrecomputationParallelized();

  virtual Force* createNewCopy();

  virtual const char* name();

  static const char* static_name();

  virtual void getAffectedVars(int pidx, std::unordered_set<int>& vars);

  virtual bool isContained(int pidx);

  virtual bool isExternal();

  const Vector3s& m_gravity;
};

}  // namespace libwethair

#endif  // LIBWETHAIR_CORE_SIMPLE_GRAVITY_FORCE_H_
