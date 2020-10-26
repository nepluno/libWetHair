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

#ifndef LIBWETHAIR_CORE_LEVEL_SET_FORCE_H_
#define LIBWETHAIR_CORE_LEVEL_SET_FORCE_H_

#include <Eigen/Core>
#include <iostream>
#include <stack>

#include "Force.h"

class FluidSim;

template <int DIM>
class TwoDScene;

template <int DIM>
class LevelSetForce : public Force {
  struct ParticleLSPair {
    scalar dist;
    scalar max_dist;
    bool valid;
  };

  struct PointLSPair {
    int eidx;
    int k_gauss;
    scalar alpha_point;
    scalar V;
    scalar quadrature_weight;
    scalar pressure_weight;
    scalar viscous_phi;
    Vectors<DIM> x0;
  };

 public:
  const int m_num_quadrature = 1;

  LevelSetForce(TwoDScene<DIM>* parent, FluidSim* fluidsim, int hidx);

  virtual ~LevelSetForce();

  virtual void preCompute(const VectorXs& x, const VectorXs& v,
                          const VectorXs& m, const scalar& dt);

  virtual void computeIntegrationVars(const VectorXs& x, const VectorXs& v,
                                      const VectorXs& m, VectorXs& lambda,
                                      VectorXs& lambda_v, TripletXs& J,
                                      TripletXs& Jv, TripletXs& Jxv,
                                      TripletXs& tildeK, TripletXs& stiffness,
                                      TripletXs& damping, VectorXs& Phi,
                                      VectorXs& Phiv, const scalar& dt);

  virtual void postStepScene(const scalar& dt);

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

  virtual const char* name();

  virtual void getAffectedVars(int pidx, std::unordered_set<int>& vars);

  virtual int getAffectedHair(const std::vector<int> particle_to_hairs);

  virtual bool isContained(int pidx);

 private:
  int m_hidx;
  FluidSim* m_fluidsim;
  TwoDScene<DIM>* m_parent;

  std::vector<ParticleLSPair> m_particle_ls_pairs;
  std::vector<PointLSPair> m_point_ls_pairs;
};

#endif  // LIBWETHAIR_CORE_LEVEL_SET_FORCE_H_
