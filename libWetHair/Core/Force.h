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

#ifndef LIBWETHAIR_CORE_FORCE_H_
#define LIBWETHAIR_CORE_FORCE_H_

#include <Eigen/Core>
#include <unordered_set>

#include "MathDefs.h"

class Force {
 protected:
  int m_internal_index_pos;
  int m_internal_index_vel;

  int m_internal_index_J;
  int m_internal_index_Jv;
  int m_internal_index_Jxv;
  int m_internal_index_tildeK;

 public:
  virtual ~Force();

  virtual void addEnergyToTotal(const VectorXs& x, const VectorXs& v,
                                const VectorXs& m, scalar& E) {}
  virtual void addGradEToTotal(const VectorXs& x, const VectorXs& v,
                               const VectorXs& m, VectorXs& gradE) {}

  virtual void addHessXToTotal(const VectorXs& x, const VectorXs& v,
                               const VectorXs& m, TripletXs& hessE) {}

  virtual void addHessVToTotal(const VectorXs& x, const VectorXs& v,
                               const VectorXs& m, TripletXs& hessE) {}

  virtual void addGradEToTotal(const VectorXs& x, const VectorXs& v,
                               const VectorXs& m, VectorXs& gradE, int pidx) {}

  virtual void addHessXToTotal(const VectorXs& x, const VectorXs& v,
                               const VectorXs& m, VectorXs& hessE, int pidx) {}

  virtual void addHessVToTotal(const VectorXs& x, const VectorXs& v,
                               const VectorXs& m, VectorXs& hessE, int pidx) {}

  virtual void computeIntegrationVars(const VectorXs& x, const VectorXs& v,
                                      const VectorXs& m, VectorXs& lambda,
                                      VectorXs& lambda_v, TripletXs& J,
                                      TripletXs& Jv, TripletXs& Jxv,
                                      TripletXs& tildeK, TripletXs& stiffness,
                                      TripletXs& damping, VectorXs& Phi,
                                      VectorXs& Phiv, const scalar& dt) {}

  virtual int numConstraintPos() = 0;

  virtual int numConstraintVel() = 0;

  virtual int numJ() = 0;

  virtual int numJv() = 0;

  virtual int numJxv() = 0;

  virtual int numTildeK() = 0;

  virtual bool isParallelized() = 0;

  virtual bool isPrecomputationParallelized() = 0;

  virtual void storeLambda(const VectorXs& lambda, const VectorXs& lambda_v);

  virtual void setInternalIndex(int index_pos, int index_vel, int index_J,
                                int index_Jv, int index_Jxv, int index_tildeK);

  virtual Force* createNewCopy() = 0;

  virtual const char* name() = 0;

  virtual void getAffectedVars(int pidx, std::unordered_set<int>& vars) = 0;

  virtual int getAffectedHair(const std::vector<int> particle_to_hairs);

  virtual bool isContained(int pidx) = 0;

  virtual bool isExternal();

  virtual void preCompute(const VectorXs& x, const VectorXs& v,
                          const VectorXs& m, const scalar& dt);

  virtual bool isInterHair() const;

  virtual void postStepScene(const scalar& dt);
};

#endif  // LIBWETHAIR_CORE_FORCE_H_
