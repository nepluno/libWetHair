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

#include "DragDampingForce.h"

#include <tbb/tbb.h>

#include "HairFlow.h"
#include "TwoDScene.h"

namespace libwethair {

template <int DIM>
DragDampingForce<DIM>::DragDampingForce(const TwoDScene<DIM>& scene,
                                        const scalar& b, int hidx)
    : Force(),
      m_b(b),
      m_dofs(scene.getFilmFlows()[hidx]->size() * DIM),
      m_hidx(hidx),
      m_scene(scene) {
  assert(m_b >= 0.0);

  m_lambda_v.resize(m_dofs);
  m_lambda_v.setZero();
}

template <int DIM>
DragDampingForce<DIM>::~DragDampingForce() {}

template <int DIM>
void DragDampingForce<DIM>::computeIntegrationVars(
    const VectorXs& x, const VectorXs& v, const VectorXs& m, VectorXs& lambda,
    VectorXs& lambda_v, TripletXs& J, TripletXs& Jv, TripletXs& Jxv,
    TripletXs& tildeK, TripletXs& stiffness, TripletXs& damping, VectorXs& Phi,
    VectorXs& Phiv, const scalar& dt) {
  const HairFlow<DIM>* flow = m_scene.getFilmFlows()[m_hidx];

  int np = flow->size();
  const std::vector<int>& particles = flow->getParticleIndices();

  tbb::parallel_for(0, np, 1, [&](int k) {
    int pidx = particles[k];
    for (int r = 0; r < DIM; ++r) {
      int i = m_scene.getDof(pidx) + r;
      Jv[m_internal_index_Jv + k * DIM + r] =
          Triplets(m_internal_index_vel + k * DIM + r, i, 1.0);
      damping[m_internal_index_vel + k * DIM + r] =
          Triplets(m_internal_index_vel + k * DIM + r,
                   m_internal_index_vel + k * DIM + r, m_b);
      Phiv[m_internal_index_vel + k * DIM + r] = v(i);
      lambda_v[m_internal_index_vel + k * DIM + r] = m_lambda_v(k * DIM + r);
    }
  });
}

template <int DIM>
int DragDampingForce<DIM>::numJ() {
  return 0;
}

template <int DIM>
int DragDampingForce<DIM>::numJv() {
  const HairFlow<DIM>* flow = m_scene.getFilmFlows()[m_hidx];

  int np = flow->size();

  return np * DIM;
}

template <int DIM>
int DragDampingForce<DIM>::numJxv() {
  return 0;
}

template <int DIM>
int DragDampingForce<DIM>::numTildeK() {
  return 0;
}

template <int DIM>
bool DragDampingForce<DIM>::isParallelized() {
  return false;
}

template <int DIM>
bool DragDampingForce<DIM>::isPrecomputationParallelized() {
  return true;
}

template <int DIM>
int DragDampingForce<DIM>::numConstraintPos() {
  return 0;
}

template <int DIM>
int DragDampingForce<DIM>::numConstraintVel() {
  const HairFlow<DIM>* flow = m_scene.getFilmFlows()[m_hidx];

  int np = flow->size();

  return np * DIM;
}

template <int DIM>
const char* DragDampingForce<DIM>::name() {
  return "dragdamping";
}

template <int DIM>
void DragDampingForce<DIM>::storeLambda(const VectorXs& lambda,
                                        const VectorXs& lambda_v) {
  m_lambda_v = lambda_v.segment(m_internal_index_vel, m_dofs);
}

template <int DIM>
Force* DragDampingForce<DIM>::createNewCopy() {
  return new DragDampingForce(*this);
}

template <int DIM>
void DragDampingForce<DIM>::getAffectedVars(int colidx,
                                            std::unordered_set<int>& vars) {
  if (m_scene.getComponent(colidx) == DIM) return;
  vars.insert(colidx);
}

template <int DIM>
int DragDampingForce<DIM>::getAffectedHair(
    const std::vector<int> particle_to_hairs) {
  return m_hidx;
}

template <int DIM>
bool DragDampingForce<DIM>::isContained(int colidx) {
  int idir = m_scene.getComponent(colidx);
  if (idir == DIM) return false;
  int pidx = m_scene.getVertFromDof(colidx);
  const std::vector<int>& particle_hairs = m_scene.getParticleToHairs();
  return particle_hairs[pidx] == m_hidx;
}

// explicit instantiations at bottom
template class DragDampingForce<2>;
template class DragDampingForce<3>;

}  // namespace libwethair
