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

#include "StrandCompliantEuler.h"

#include "HairFlow.h"
#include "MathUtilities.h"
#include "StrandCompliantManager.h"

namespace libwethair {

template <int DIM>
StrandCompliantEuler<DIM>::StrandCompliantEuler(
    StrandCompliantManager<DIM>* parent, int hidx)
    : m_parent(parent), m_hidx(hidx) {
  HairFlow<DIM>* flow = m_parent->m_scene->getFilmFlows()[hidx];
  int ndof = m_parent->m_scene->isMassSpring() ? flow->size() * DIM
                                               : flow->size() * 4 - 1;
  m_A.resize(ndof, ndof);

  m_num_global_dof = ndof;
  m_start_global_dof = m_parent->m_scene->getDof(flow->getParticleIndices()[0]);
}

template <int DIM>
StrandCompliantEuler<DIM>::~StrandCompliantEuler() {}

template <int DIM>
void StrandCompliantEuler<DIM>::computeRHSIncremental(TwoDScene<DIM>& scene,
                                                      scalar dt, VectorXs& b,
                                                      const VectorXs& vplus) {
  VectorXs dv = vplus.segment(m_start_global_dof, m_num_global_dof) -
                scene.getV().segment(m_start_global_dof, m_num_global_dof);
  HairFlow<DIM>* flow = m_parent->m_scene->getFilmFlows()[m_hidx];
  const VectorXs& m =
      scene.getInterpolatedM().segment(m_start_global_dof, m_num_global_dof);
  const VectorXs& gradu =
      m_parent->m_gradU.segment(m_start_global_dof, m_num_global_dof);
  const Vector6i& constraint_idx = flow->getConstraintIdx();
  const Vector6i& constraint_num = flow->getNumConstraints();

  const VectorXs& phi =
      m_parent->m_Phi.segment(constraint_idx(0), constraint_num(0));

  b.segment(m_start_global_dof, m_num_global_dof) =
      -dv.cwiseProduct(m) - gradu * dt -
      dt * (m_J.transpose() * (m_invC * phi));
  if (constraint_num(1) > 0) {
    const VectorXs& phi_v =
        m_parent->m_Phi_v.segment(constraint_idx(1), constraint_num(1));
    b.segment(m_start_global_dof, m_num_global_dof) +=
        -dt * (m_Jv.transpose() * (m_invCv * phi_v));
  }
}

template <int DIM>
void StrandCompliantEuler<DIM>::computeRHS(TwoDScene<DIM>& scene, scalar dt,
                                           VectorXs& b) {
  const VectorXs& v =
      scene.getV().segment(m_start_global_dof, m_num_global_dof);
  HairFlow<DIM>* flow = m_parent->m_scene->getFilmFlows()[m_hidx];
  const VectorXs& m =
      scene.getInterpolatedM().segment(m_start_global_dof, m_num_global_dof);
  const VectorXs& gradu =
      m_parent->m_gradU.segment(m_start_global_dof, m_num_global_dof);
  const Vector6i& constraint_idx = flow->getConstraintIdx();
  const Vector6i& constraint_num = flow->getNumConstraints();

  const VectorXs& phi =
      m_parent->m_Phi.segment(constraint_idx(0), constraint_num(0));

  b.segment(m_start_global_dof, m_num_global_dof) =
      v.cwiseProduct(m) - gradu * dt - dt * (m_J.transpose() * (m_invC * phi));

  if (constraint_num(1) > 0) {
    const VectorXs& phi_v =
        m_parent->m_Phi_v.segment(constraint_idx(1), constraint_num(1));
    b.segment(m_start_global_dof, m_num_global_dof) +=
        -dt * (m_Jv.transpose() * (m_invCv * (phi_v - m_Jv * v)));
  }
}

template <int DIM>
void StrandCompliantEuler<DIM>::computeAp(const VectorXs& p, VectorXs& b) {
  b.segment(m_start_global_dof, m_num_global_dof) =
      m_A * p.segment(m_start_global_dof, m_num_global_dof);
}

template <int DIM>
bool StrandCompliantEuler<DIM>::PreconditionScene(TwoDScene<DIM>& scene,
                                                  scalar dt, VectorXs& r,
                                                  const VectorXs& b) {
  if (m_parent->m_use_preconditioner) {
    return stepScene(scene, dt, r, b);
  } else {
    const VectorXs& rr = b.segment(m_start_global_dof, m_num_global_dof);
    r.segment(m_start_global_dof, m_num_global_dof) =
        m_A_inv_diag.cwiseProduct(rr);
    return true;
  }
}

template <int DIM>
bool StrandCompliantEuler<DIM>::stepScene(TwoDScene<DIM>& scene, scalar dt,
                                          VectorXs& r, const VectorXs& b) {
  const VectorXs& rr = b.segment(m_start_global_dof, m_num_global_dof);
  r.segment(m_start_global_dof, m_num_global_dof) =
      m_solver.solve(rr);  // m_A_diag.cwiseProduct(rr);
  return true;
}

template <int DIM>
bool StrandCompliantEuler<DIM>::stepScene(TwoDScene<DIM>& scene, scalar dt) {
  return false;
}

template <int DIM>
void StrandCompliantEuler<DIM>::updateLambda(TwoDScene<DIM>& scene,
                                             const VectorXs& dx,
                                             const VectorXs& dv, scalar dt) {
  HairFlow<DIM>* flow = scene.getFilmFlows()[m_hidx];
  const std::vector<int>& particles = flow->getParticleIndices();

  const Vector6i& constraint_idx = flow->getConstraintIdx();
  const Vector6i& constraint_num = flow->getNumConstraints();

  const VectorXs& phi =
      m_parent->m_Phi.segment(constraint_idx(0), constraint_num(0));

  const VectorXs& ddx = dx.segment(m_start_global_dof, m_num_global_dof);
  const VectorXs& ddv = dv.segment(m_start_global_dof, m_num_global_dof);

  m_parent->m_lambda.segment(constraint_idx(0), constraint_num(0)) =
      -m_invC * (m_J * ddx + phi);
  if (constraint_num(1) > 0) {
    const VectorXs& phi_v =
        m_parent->m_Phi_v.segment(constraint_idx(1), constraint_num(1));
    m_parent->m_lambda_v.segment(constraint_idx(1), constraint_num(1)) =
        -m_invCv * (m_Jxv * ddx + m_Jv * ddv + phi_v);
  }
}

template <int DIM>
void StrandCompliantEuler<DIM>::updateNextV(TwoDScene<DIM>& scene,
                                            const VectorXs& vplus) {
  VectorXs& v_next = scene.getV();
  HairFlow<DIM>* flow = scene.getFilmFlows()[m_hidx];
  const std::vector<int>& particles = flow->getParticleIndices();
  for (int pidx : particles) {
    if (!scene.isFixed(pidx)) {
      int numdofs = scene.isMassSpring() || scene.isTip(pidx) ? DIM : 4;
      v_next.segment(scene.getDof(pidx), numdofs) =
          vplus.segment(scene.getDof(pidx), numdofs);
    }
  }
}

template <int DIM>
void StrandCompliantEuler<DIM>::preIterate(TwoDScene<DIM>& scene, scalar dt) {
  // Foreach local forces:
  //	compute an A matrix (localized version of compute Integration vars),
  HairFlow<DIM>* flow = m_parent->m_scene->getFilmFlows()[m_hidx];
  int ndof = scene.isMassSpring() ? flow->size() * DIM : flow->size() * 4 - 1;
  const VectorXs& m = scene.getInterpolatedM();
  const Vector6i& constraint_idx = flow->getConstraintIdx();
  const Vector6i& constraint_num = flow->getNumConstraints();
  const std::vector<int>& global_local =
      m_parent->m_scene->getParticleToHairLocalIndices();
  const std::vector<int>& particle_hair =
      m_parent->m_scene->getParticleToHairs();

  m_A_nz.resize(0);
  m_A_nz.reserve(constraint_num(5) + ndof);
  const TripletXs& m_A_nz_ref = m_parent->m_A_nz;
  int nK = constraint_num(5);
  int base_K = constraint_idx(5);
  for (int i = 0; i < nK; ++i) {
    const Triplets& t = m_A_nz_ref[base_K + i];
    int ip = scene.getVertFromDof(t.row());
    int jp = scene.getVertFromDof(t.col());
    int ip_local = global_local[ip];
    int jp_local = global_local[jp];
    int idir = scene.getComponent(t.row());
    int jdir = scene.getComponent(t.col());

    scalar val = 0;
    if (!scene.isFixed(ip) && !scene.isFixed(jp)) {
      val += t.value() * dt * dt;
    }

    if (val != 0.0) {
      if (scene.isMassSpring())
        m_A_nz.push_back(
            Triplets(ip_local * DIM + idir, jp_local * DIM + jdir, val));
      else
        m_A_nz.push_back(
            Triplets(ip_local * 4 + idir, jp_local * 4 + jdir, val));
    }
  }

  for (int i = 0; i < ndof; ++i) {
    if (scene.isFixed(scene.getVertFromDof(m_start_global_dof + i)))
      m_A_nz.push_back(Triplets(i, i, 1.0));
    else
      m_A_nz.push_back(Triplets(i, i, m(m_start_global_dof + i)));
  }

  m_A.setFromTriplets(m_A_nz.begin(), m_A_nz.end());

  int nconstraint = constraint_num(0);
  int nconstraint_v = constraint_num(1);
  int base_constraint = constraint_idx(0);
  int base_constraint_v = constraint_idx(1);

  m_J_nz.resize(0);
  m_J_nz.reserve(constraint_num(2));
  int nJ = constraint_num(2);
  int base_J = constraint_idx(2);
  const TripletXs& m_J_nz_ref = m_parent->m_J_nz;
  for (int i = 0; i < nJ; ++i) {
    const Triplets& t = m_J_nz_ref[base_J + i];
    int jp = scene.getVertFromDof(t.col());
    if (scene.isFixed(jp)) continue;
    int jp_local = global_local[jp];
    int jdir = scene.getComponent(t.col());
    if (scene.isMassSpring())
      m_J_nz.push_back(Triplets(t.row() - base_constraint,
                                jp_local * DIM + jdir, t.value()));
    else
      m_J_nz.push_back(
          Triplets(t.row() - base_constraint, jp_local * 4 + jdir, t.value()));
  }
  m_J.resize(nconstraint, ndof);
  m_J.setFromTriplets(m_J_nz.begin(), m_J_nz.end());

  m_Jv_nz.resize(0);
  m_Jv_nz.reserve(constraint_num(3));
  int nJv = constraint_num(3);
  int base_Jv = constraint_idx(3);
  const TripletXs& m_Jv_nz_ref = m_parent->m_Jv_nz;
  for (int i = 0; i < nJv; ++i) {
    const Triplets& t = m_Jv_nz_ref[base_Jv + i];
    int jp = scene.getVertFromDof(t.col());
    if (scene.isFixed(jp)) continue;
    int jp_local = global_local[jp];
    int jdir = scene.getComponent(t.col());
    if (scene.isMassSpring())
      m_Jv_nz.push_back(Triplets(t.row() - base_constraint_v,
                                 jp_local * DIM + jdir, t.value()));
    else
      m_Jv_nz.push_back(Triplets(t.row() - base_constraint_v,
                                 jp_local * 4 + jdir, t.value()));
  }
  m_Jv.resize(nconstraint_v, ndof);
  m_Jv.setFromTriplets(m_Jv_nz.begin(), m_Jv_nz.end());

  m_Jxv_nz.resize(0);
  m_Jxv_nz.reserve(constraint_num(4));
  int nJxv = constraint_num(4);
  int base_Jxv = constraint_idx(4);
  const TripletXs& m_Jxv_nz_ref = m_parent->m_Jxv_nz;
  for (int i = 0; i < nJxv; ++i) {
    const Triplets& t = m_Jxv_nz_ref[base_Jxv + i];
    int jp = scene.getVertFromDof(t.col());
    if (scene.isFixed(jp)) continue;
    int jp_local = global_local[jp];
    int jdir = scene.getComponent(t.col());
    if (scene.isMassSpring())
      m_Jxv_nz.push_back(Triplets(t.row() - base_constraint_v,
                                  jp_local * DIM + jdir, t.value()));
    else
      m_Jxv_nz.push_back(Triplets(t.row() - base_constraint_v,
                                  jp_local * 4 + jdir, t.value()));
  }
  m_Jxv.resize(nconstraint_v, ndof);
  m_Jxv.setFromTriplets(m_Jxv_nz.begin(), m_Jxv_nz.end());

  m_invC_nz.resize(nconstraint);
  const TripletXs& m_invC_nz_ref = m_parent->m_invC_nz;
  for (int i = 0; i < nconstraint; ++i) {
    const Triplets& t = m_invC_nz_ref[base_constraint + i];
    m_invC_nz[i] = Triplets(i, i, t.value());
  }
  m_invC.resize(nconstraint, nconstraint);
  m_invC.setFromTriplets(m_invC_nz.begin(), m_invC_nz.end());

  m_invCv_nz.resize(nconstraint_v);
  const TripletXs& m_invCv_nz_ref = m_parent->m_invCv_nz;
  for (int i = 0; i < nconstraint_v; ++i) {
    const Triplets& t = m_invCv_nz_ref[base_constraint_v + i];
    m_invCv_nz[i] = Triplets(i, i, t.value());
  }
  m_invCv.resize(nconstraint_v, nconstraint_v);
  m_invCv.setFromTriplets(m_invCv_nz.begin(), m_invCv_nz.end());

  // pre-factor it with simplicial LDLT
  m_A += (SparseXs(m_J.transpose()) * ((m_invC * dt * dt) * m_J));

  if (nconstraint_v > 0) {
    m_A +=
        (SparseXs(m_Jv.transpose()) * ((m_invCv * dt) * (m_Jv + m_Jxv * dt)));
  }

  m_solver.compute(m_A);

  if (!m_parent->m_use_preconditioner) {
    m_A_inv_diag = m_A.diagonal().cwiseInverse();
  }
}

// explicit instantiations at bottom
template class StrandCompliantEuler<2>;
template class StrandCompliantEuler<3>;

}  // namespace libwethair
