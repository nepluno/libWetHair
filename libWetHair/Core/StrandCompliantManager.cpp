//
// This file is part of the libWetHair open source project
//
// The code is licensed solely for academic and non-commercial use under the
// terms of the Clear BSD License. The terms of the Clear BSD License are
// provided below. Other licenses may be obtained by contacting the faculty
// of the Columbia Computer Graphics Group or a Columbia University licensing officer.
//
// The Clear BSD License
//
// Copyright 2017 Yun (Raymond) Fei, Henrique Teles Maia, Christopher Batty,
// Changxi Zheng, and Eitan Grinspun
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted (subject to the limitations in the disclaimer
// below) provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice, this
//  list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright notice,
//  this list of conditions and the following disclaimer in the documentation
//  and/or other materials provided with the distribution.
//
// * Neither the name of the copyright holder nor the names of its contributors may be used
//  to endorse or promote products derived from this software without specific
//  prior written permission.
//
// NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE GRANTED BY THIS
// LICENSE. THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
// GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
// OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.

#include "StrandCompliantManager.h"
#include "StrandCompliantEuler.h"
#include "TwoDScene.h"
#include "ThreadUtils.h"
#include "TimingUtilities.h"
#include <unordered_map>

using namespace threadutils;

template<int DIM>
StrandCompliantManager<DIM>::~StrandCompliantManager()
{
  int nintegrators = m_integrators.size();
  for(int i = 0; i < nintegrators; ++i)
  {
    delete m_integrators[i];
  }
  m_integrators.clear();
}

template<int DIM>
StrandCompliantManager<DIM>::StrandCompliantManager(TwoDScene<DIM>* scene, int max_newton, int max_iters, scalar criterion, bool compute_interhair, bool use_preconditioner)
: m_max_num_newton(max_newton), m_max_num_iters(max_iters), m_criterion(criterion), m_scene(scene), m_compute_interhair(compute_interhair), m_use_preconditioner(use_preconditioner)
{
  int numhairs = scene->getNumFlows();

  m_integrators.resize(numhairs);
  for(int i = 0; i < numhairs; ++i)
  {
    m_integrators[i] = new StrandCompliantEuler<DIM>(this, i);
  }
  
  m_dt = 0.0;
  
  int ndof = scene->getNumDofs();
  m_gradU.resize(ndof);
  m_r.resize(ndof);
  m_v.resize(ndof);
  m_p.resize(ndof);
  m_q.resize(ndof);
  m_t.resize(ndof);
  m_z.resize(ndof);
  m_rhs.resize(ndof);
  m_vplus.resize(ndof);
  
  m_dv.resize(ndof);
  m_dx.resize(ndof);
  m_dvi.resize(ndof);
  m_dxi.resize(ndof);
  m_dx_scripted.resize(ndof);
  
  SceneStepper<DIM>::m_timing_statistics.resize(6, 0.0);
}

template<int DIM>
void StrandCompliantManager<DIM>::localPreIterate(TwoDScene<DIM> & scene, scalar dt)
{
  int nhairs = m_integrators.size();
  
  thread_pool::ParallelFor(0, nhairs, [&] (int hidx) {
    m_integrators[hidx]->preIterate(scene, dt);
  });
}

template<int DIM>
void StrandCompliantManager<DIM>::localStepScene(TwoDScene<DIM> & scene, scalar dt, VectorXs& r, const VectorXs& b)
{
  scalar t0 = timingutils::seconds();
  int nhairs = m_integrators.size();
  thread_pool::ParallelFor(0, nhairs, [&] (int hidx) {
    m_integrators[hidx]->stepScene(scene, dt, r, b);
  });
  scalar t1 = timingutils::seconds();
  SceneStepper<DIM>::m_timing_statistics[4] += (t1 - t0); // local Solve
  t0 = t1;
}

template<int DIM>
void StrandCompliantManager<DIM>::localPreconditionScene(TwoDScene<DIM> & scene, scalar dt, VectorXs& r, const VectorXs& b)
{
  scalar t0 = timingutils::seconds();
  int nhairs = m_integrators.size();
  thread_pool::ParallelFor(0, nhairs, [&] (int hidx) {
    m_integrators[hidx]->PreconditionScene(scene, dt, r, b);
  });
  scalar t1 = timingutils::seconds();
  SceneStepper<DIM>::m_timing_statistics[4] += (t1 - t0); // local Solve
  t0 = t1;
}

template<int DIM>
void StrandCompliantManager<DIM>::localUpdateLambda(TwoDScene<DIM> & scene, const VectorXs& dx, const VectorXs& dv, scalar dt)
{
  int nhairs = m_integrators.size();
  thread_pool::ParallelFor(0, nhairs, [&] (int hidx) {
    m_integrators[hidx]->updateLambda(scene, dx, dv, dt);
  });
}

template<int DIM>
void StrandCompliantManager<DIM>::updateNextV(TwoDScene<DIM> & scene, const VectorXs& vplus)
{
  int nhairs = m_integrators.size();
  thread_pool::ParallelFor(0, nhairs, [&] (int hidx) {
    m_integrators[hidx]->updateNextV(scene, vplus);
  });
}


template<int DIM>
void StrandCompliantManager<DIM>::zeroFixedDoFs( const TwoDScene<DIM>& scene, VectorXs& vec )
{
  int nprts = scene.getNumParticles();
  for( int i = 0; i < nprts; ++i ){
    if( scene.isFixed(i) ){
      int numdofs = scene.isMassSpring() || scene.isTip(i) ? DIM : 4;
      vec.segment(scene.getDof(i), numdofs).setZero();
    }
  }
}

template<int DIM>
void StrandCompliantManager<DIM>::computeRHSIncrementalLocal( TwoDScene<DIM> & scene, scalar dt, VectorXs& b, const VectorXs& vplus )
{
  m_gradU.setZero();
  scene.accumulateExternalGradU(m_gradU);
  zeroFixedDoFs(scene, m_gradU);
  
  // compute local b
  int nhairs = m_integrators.size();
  thread_pool::ParallelFor(0, nhairs, [&] (int hidx) {
    m_integrators[hidx]->computeRHSIncremental(scene, dt, b, vplus);
  });
}

template<int DIM>
void StrandCompliantManager<DIM>::computeRHSIncrementalInterhair( TwoDScene<DIM> & scene, scalar dt, VectorXs& b, const VectorXs& vplus )
{
  if(m_interhair_num(0) > 0) {
    mathutils::compute_cwiseProduct(m_c_buffer, m_invC_interhair, m_Phi_interhair);
    mathutils::accumulateJTPhi_coeff(b, -dt, m_c_buffer, m_J_interhair);
  }
  
  if(m_interhair_num(1) > 0) {
    mathutils::compute_cwiseProduct(m_cv_buffer, m_invCv_interhair, m_Phiv_interhair);
    mathutils::accumulateJTPhi_coeff(b, -dt, m_cv_buffer, m_Jv_interhair);
  }
}

template<int DIM>
void StrandCompliantManager<DIM>::computeRHSLocal( TwoDScene<DIM> & scene, scalar dt, VectorXs& b )
{
  m_gradU.setZero();
  scene.accumulateExternalGradU(m_gradU);
  zeroFixedDoFs(scene, m_gradU);
  
  // compute local b
  int nhairs = m_integrators.size();
  thread_pool::ParallelFor(0, nhairs, [&] (int hidx) {
    m_integrators[hidx]->computeRHS(scene, dt, b);
  });
}

template<int DIM>
void StrandCompliantManager<DIM>::computeAp( TwoDScene<DIM> & scene, const VectorXs& p, VectorXs& b, scalar dt )
{
  scalar t0 = timingutils::seconds();
  
  // compute local Ap
  int nhairs = m_integrators.size();
  thread_pool::ParallelFor(0, nhairs, [&] (int hidx) {
    m_integrators[hidx]->computeAp(p, b);
  });
  
  if(m_interhair_num(2) > 0) {
    mathutils::computeJTPhi_coeff(m_c_buffer, m_invC_interhair, p, m_JT_interhair);
    mathutils::accumulateJTPhi_coeff(b, dt * dt, m_c_buffer, m_J_interhair);
  }
  
  if(m_interhair_num(3) > 0) {
    mathutils::computeJTPhi_coeff(m_cv_buffer, m_invCv_interhair, p, m_JvT_interhair);
    mathutils::accumulateJTPhi_coeff(b, dt, m_cv_buffer, m_Jv_interhair);
  }
  
  scalar t1 = timingutils::seconds();
  SceneStepper<DIM>::m_timing_statistics[5] += (t1 - t0); // global CG
  t0 = t1;
}

template<int DIM>
void StrandCompliantManager<DIM>::localUpdateNumConstraints(const VectorXs& dx, const VectorXs& dv, scalar dt)
{
  int num_pos, num_vel, num_J, num_Jv, num_Jxv, num_tildeK;
  m_scene->updateNumConstraintsLocal(num_pos, num_vel, num_J, num_Jv, num_Jxv, num_tildeK);
  
  m_lambda.resize(num_pos);
  m_lambda_v.resize(num_vel);
  m_A_nz.resize(num_tildeK);
  m_J_nz.resize(num_J);
  m_Jv_nz.resize(num_Jv);
  m_Jxv_nz.resize(num_Jxv);
  m_invC_nz.resize(num_pos);
  m_invCv_nz.resize(num_vel);
  m_Phi.resize(num_pos);
  m_Phi_v.resize(num_vel);
  
  m_scene->localPostPreprocess(m_lambda, m_lambda_v, m_J_nz, m_Jv_nz, m_Jxv_nz, m_A_nz, m_invC_nz, m_invCv_nz, m_Phi, m_Phi_v, dx, dv, dt);
}

template<int DIM>
void StrandCompliantManager<DIM>::interhairUpdateNumConstraints(const VectorXs& dx, const VectorXs& dv, scalar dt)
{
  int num_pos, num_vel, num_J, num_Jv, num_Jxv, num_tildeK;
  m_scene->updateNumConstraintsInterHair(num_pos, num_vel, num_J, num_Jv, num_Jxv, num_tildeK, m_interhair_idx, m_interhair_num);
  
  m_lambda.conservativeResize(num_pos);
  m_lambda_v.conservativeResize(num_vel);
  m_A_nz.resize(num_tildeK);
  m_J_nz.resize(num_J);
  m_Jv_nz.resize(num_Jv);
  m_Jxv_nz.resize(num_Jxv);
  m_invC_nz.resize(num_pos);
  m_invCv_nz.resize(num_vel);
  m_Phi.conservativeResize(num_pos);
  m_Phi_v.conservativeResize(num_vel);
  
  m_scene->interhairPostPreprocess(m_lambda, m_lambda_v, m_J_nz, m_Jv_nz, m_Jxv_nz, m_A_nz, m_invC_nz, m_invCv_nz, m_Phi, m_Phi_v, dx, dv, dt);
  
  m_c_buffer.resize(m_interhair_num(0));
  m_cv_buffer.resize(m_interhair_num(1));
  
  int ndof = m_scene->getNumDofs();
  
  int num_J_inter = m_interhair_num(2);
  m_J_interhair_nz.resize(num_J_inter);
  thread_pool::ParallelFor(0, num_J_inter, [&] (int i) {
    const Triplets& t = m_J_nz[i + m_interhair_idx(2)];
    m_J_interhair_nz[i] = Triplets(t.row() - m_interhair_idx(0), t.col(), t.value());
  });
  m_J_interhair.resize(m_interhair_num(0), ndof);
  m_J_interhair.setFromTriplets(m_J_interhair_nz.begin(), m_J_interhair_nz.end());
  m_JT_interhair = m_J_interhair.transpose();
  m_J_interhair.makeCompressed();
  m_JT_interhair.makeCompressed();
  
  int num_Jv_inter = m_interhair_num(3);
  m_Jv_interhair_nz.resize(num_Jv_inter);
  thread_pool::ParallelFor(0, num_Jv_inter, [&] (int i) {
    const Triplets& t = m_Jv_nz[i + m_interhair_idx(3)];
    m_Jv_interhair_nz[i] = Triplets(t.row() - m_interhair_idx(1), t.col(), t.value());
  });
  m_Jv_interhair.resize(m_interhair_num(1), ndof);
  m_Jv_interhair.setFromTriplets(m_Jv_interhair_nz.begin(), m_Jv_interhair_nz.end());
  m_JvT_interhair = m_Jv_interhair.transpose();
  m_Jv_interhair.makeCompressed();
  m_JvT_interhair.makeCompressed();
  
  m_invC_interhair.resize(m_interhair_num(0));
  thread_pool::ParallelFor(0, m_interhair_num(0), [&] (int i) {
    m_invC_interhair(i) = m_invC_nz[i + m_interhair_idx(0)].value();
  });
  
  m_invCv_interhair.resize(m_interhair_num(1));
  thread_pool::ParallelFor(0, m_interhair_num(1), [&] (int i) {
    m_invCv_interhair(i) = m_invCv_nz[i + m_interhair_idx(1)].value();
  });
  
  m_Phi_interhair = m_Phi.segment(m_interhair_idx(0), m_interhair_num(0));
  m_Phiv_interhair = m_Phi_v.segment(m_interhair_idx(1), m_interhair_num(1));
}

#define PCG_VERBOSE
#define NEWTON_VERBOSE

template<int DIM>
bool StrandCompliantManager<DIM>::stepSceneNonlinear( TwoDScene<DIM> & scene, scalar dt, bool updatePreCompute )
{
  m_dt = dt;
  m_scene = &scene;
  
  VectorXs& x = scene.getX();
  VectorXs& v = scene.getV();
  
  m_dv.setZero();
  m_dx = v * dt;
  
  std::cout << "[pre-compute-local]" << std::endl;
  scene.preComputeLocal(m_dx, m_dv, dt);
  
  std::cout << "[compute-assist-vars-local]" << std::endl;
  localUpdateNumConstraints(VectorXs(), VectorXs(), dt);
  
  localPreIterate(scene, dt);
  computeRHSLocal(scene, dt, m_rhs);
  
  // solve Mv_+=b
  localStepScene(scene, dt, m_vplus, m_rhs);

  m_dvi = m_vplus - v;
  m_dxi = m_dvi * dt;
  
  m_dv += m_dvi;
  m_dx += m_dxi;
  
  std::cout << "[update-lambda-local]" << std::endl;
  localUpdateLambda(scene, m_dxi, m_dvi, dt);
  m_scene->storeLambda(m_lambda, m_lambda_v);
  
  if(m_compute_interhair) {
    std::cout << "[pre-compute-interhair]" << std::endl;
    scene.preComputeInterhair(m_dx, m_dv, dt);
  }
  
  // Start Newton-Krylov Iteration
  int nkiters = 0;
  scalar newton_res_norm = 0.0;
  
  for(; nkiters < m_max_num_newton; ++nkiters) {
    std::cout << "[compute-assist-vars-local]" << std::endl;
    localUpdateNumConstraints(m_dx, m_dv, dt);
    
    if(m_compute_interhair) {
      std::cout << "[compute-assist-vars-interhair]" << std::endl;
      interhairUpdateNumConstraints(m_dx, m_dv, dt);
    }
    
    // update Pre-conditioner
    localPreIterate(scene, dt);
    computeRHSIncrementalLocal(scene, dt, m_rhs, m_vplus);
    
    if(m_compute_interhair) {
      computeRHSIncrementalInterhair(scene, dt, m_rhs, m_vplus);
    }
    // r = b - A * 0
    m_r = m_rhs;

    newton_res_norm = m_r.norm();
    
#ifdef NEWTON_VERBOSE
    std::cout << "[Newton total iter: " << nkiters << ", res: " << newton_res_norm << "]" << std::endl;
#endif
    
    int iter = 0;
    if(newton_res_norm < m_criterion) {
      break;
    }
 
    m_dvi.setZero();
    // solve Mz = r
    localPreconditionScene(scene, dt, m_z, m_r);
    
    // p = z
    m_p = m_z;
    
    scalar rho = m_r.dot(m_z);
    
    // q = Ap
    computeAp(scene, m_p, m_q, dt);
    
    // alpha = rho / (p, q)
    scalar alpha = rho / m_p.dot(m_q);
    
    // x = x + alpha*p
    m_dvi += m_p * alpha;
    
    // r = r - alpha*q
    m_r -= m_q * alpha;
    
    scalar cg_res_norm = m_r.norm();
    
    scalar rho_old, beta;
    for(; iter < m_max_num_iters && cg_res_norm > m_criterion; ++iter)
    {
      rho_old = rho;
      
      // solve Mz = r
      localPreconditionScene(scene, dt, m_z, m_r);
      
      rho = m_r.dot(m_z);
      
      beta = rho / rho_old;
      
      // p = beta * p + z
      m_p = m_z + m_p * beta;
      
      // q = Ap
      computeAp(scene, m_p, m_q, dt);
      
      // alpha = rho / (p, q)
      alpha = rho / m_p.dot(m_q);
      
      // x = x + alpha*p
      m_dvi += m_p * alpha;
      
      // r = r - alpha*q
      m_r -= m_q * alpha;
      
      cg_res_norm = m_r.norm();
      
#ifdef PCG_VERBOSE
      std::cout << "[PCG iter: " << iter << ", res: " << cg_res_norm << "]" << std::endl;
#endif
    }
    
    std::cout << "[PCG total iter: " << iter << ", res: " << cg_res_norm << "]" << std::endl;
    
    m_dxi = m_dvi * dt;
    
    m_dv += m_dvi;
    m_vplus = v + m_dv;
    m_dx = m_vplus * dt;

    std::cout << "[update-lambda-local]" << std::endl;
    localUpdateLambda(scene, m_dxi, m_dvi, dt);
    m_scene->storeLambda(m_lambda, m_lambda_v);
  }
  
  std::cout << "[Newton total iter: " << nkiters << ", res: " << newton_res_norm << "]" << std::endl;
  
  updateNextV(scene, m_vplus);
  
  SceneStepper<DIM>::m_next_x = x + dt * v;
  
  return true;
}

template<int DIM>
void StrandCompliantManager<DIM>::write(std::vector<scalar>& buf) const
{
  int n_lambda = m_lambda.size();
  for(int i = 0; i < n_lambda; ++i)
  {
    buf.push_back(m_lambda(i));
  }
  
  int n_lambda_v = m_lambda_v.size();
  for(int i = 0; i < n_lambda_v; ++i)
  {
    buf.push_back(m_lambda_v(i));
  }
}

template<int DIM>
void StrandCompliantManager<DIM>::read(const scalar* data)
{
  m_scene->preComputeLocal(VectorXs(), VectorXs(), 0.0);
  int num_pos, num_vel, num_J, num_Jv, num_Jxv, num_tildeK;
  m_scene->updateNumConstraintsLocal(num_pos, num_vel, num_J, num_Jv, num_Jxv, num_tildeK);
  
  m_lambda.resize(num_pos);
  m_lambda_v.resize(num_vel);
  
  memcpy(m_lambda.data(), data, num_pos * sizeof(scalar) );
  memcpy(m_lambda_v.data(), data + num_pos, num_vel * sizeof(scalar) );
  
  m_scene->storeLambda(m_lambda, m_lambda_v);
}

template<int DIM>
size_t StrandCompliantManager<DIM>::size()
{
  return (m_lambda.size() + m_lambda_v.size()) * sizeof(scalar);
}


template<int DIM>
bool StrandCompliantManager<DIM>::stepSceneLinear( TwoDScene<DIM> & scene, scalar dt, bool updatePreCompute )
{
  std::cout << "[pre-compute]" << std::endl;
  scalar t0 = timingutils::seconds();
  scalar t1;
  
  m_dt = dt;
  m_scene = &scene;
  
  VectorXs& x = scene.getX();
  VectorXs& v = scene.getV();
  
  m_dx = v * dt;
  m_dv.setZero();
  m_dxi = m_dv;
  
  m_dx_scripted.setZero();
  
  int np = m_scene->getNumParticles();
  
  for(int i = 0; i < np; ++i)
  {
    if(scene.isFixed(i)) {
      int numdofs = scene.isMassSpring() || scene.isTip(i) ? DIM : 4;
      m_dx_scripted.segment( scene.getDof(i), numdofs ) = v.segment( scene.getDof(i), numdofs ) * dt;
    }
  }
  
  scene.preComputeLocal(m_dx_scripted, m_dv, dt);
  
  t1 = timingutils::seconds();
  SceneStepper<DIM>::m_timing_statistics[0] += (t1 - t0); // local precomputation
  t0 = t1;
  
  std::cout << "[compute-assist-vars]" << std::endl;
  localUpdateNumConstraints(m_dx_scripted, m_dv, dt);
  
  t1 = timingutils::seconds();
  SceneStepper<DIM>::m_timing_statistics[1] += (t1 - t0); // local Jacobian
  t0 = t1;
  
  localPreIterate(scene, dt);
  computeRHSLocal(scene, dt, m_rhs);
  
  t1 = timingutils::seconds();
  SceneStepper<DIM>::m_timing_statistics[2] += (t1 - t0); // local Matrix Composition
  t0 = t1;
  
  std::cout << "[solve-equations]" << std::endl;
  // solve Mv_+=b
  localStepScene(scene, dt, m_vplus, m_rhs);
  
  m_dv = m_vplus - v;
  m_dx = m_vplus * dt;
  
  if(!m_compute_interhair) {
    localUpdateLambda(scene, m_dx, m_dv, dt);
    
    m_scene->storeLambda(m_lambda, m_lambda_v);
    
    updateNextV(scene, m_vplus);
    
    SceneStepper<DIM>::m_next_x = x + dt * v;
    
    return true;
  }
  
  t0 = timingutils::seconds();
  
  scene.preComputeInterhair(m_dx, m_dv, dt);
  interhairUpdateNumConstraints(m_dxi, m_dv, dt);
  
  computeRHSIncrementalInterhair(scene, dt, m_rhs, m_vplus);
  
  t1 = timingutils::seconds();
  SceneStepper<DIM>::m_timing_statistics[3] += (t1 - t0); // compute Interhair Variables
  t0 = t1;
  
  // Ap = Ax0
  computeAp(scene, m_vplus, m_r, dt);
  // r = b - Ax0
  m_r = m_rhs - m_r;
  
  // solve Mz = r
  localPreconditionScene(scene, dt, m_z, m_r);
  
  // p = z
  m_p = m_z;
  
  scalar rho = m_r.dot(m_z);
  
  scalar res_norm = m_r.norm();
  
  int iter = 0;
  if(res_norm < m_criterion) {
    std::cout << "[pcg total iter: " << iter << ", res: " << res_norm << "]" << std::endl;
    
    localUpdateLambda(scene, m_dx, m_dv, dt);
    
    m_scene->storeLambda(m_lambda, m_lambda_v);
    
    updateNextV(scene, m_vplus);
    
    SceneStepper<DIM>::m_next_x = x + dt * v;
    
    return true;
  }
  
  // q = Ap
  computeAp(scene, m_p, m_q, dt);
  
  // alpha = rho / (p, q)
  scalar alpha = rho / m_p.dot(m_q);
  
  // x = x + alpha*p
  m_vplus += m_p * alpha;
  
  // r = r - alpha*q
  m_r -= m_q * alpha;
  
  res_norm = m_r.norm();
  
  scalar rho_old, beta;
  for(; iter < m_max_num_iters && res_norm > m_criterion; ++iter)
  {
    rho_old = rho;
    
    // solve Mz = r
    localPreconditionScene(scene, dt, m_z, m_r);
    
    rho = m_r.dot(m_z);
    
    beta = rho / rho_old;
    
    // p = beta * p + z
    m_p = m_z + m_p * beta;
    
    // q = Ap
    computeAp(scene, m_p, m_q, dt);
    
    // alpha = rho / (p, q)
    alpha = rho / m_p.dot(m_q);
    
    // x = x + alpha*p
    m_vplus += m_p * alpha;
    
    // r = r - alpha*q
    m_r -= m_q * alpha;
    
    res_norm = m_r.norm();
    
#ifdef PCG_VERBOSE
    std::cout << "[pcg iter: " << iter << ", res: " << res_norm << "]" << std::endl;
#endif
  }
  
  std::cout << "[pcg total iter: " << iter << ", res: " << res_norm << "]" << std::endl;
  
  m_dv = m_vplus - v;
  m_dx = m_vplus * dt;
  
  localUpdateLambda(scene, m_dx, m_dv, dt);
  
  m_scene->storeLambda(m_lambda, m_lambda_v);
  
  updateNextV(scene, m_vplus);
  
  SceneStepper<DIM>::m_next_x = x + dt * v;
  
  return true;
}

template<int DIM>
bool StrandCompliantManager<DIM>::stepScene( TwoDScene<DIM> & scene, scalar dt, bool updatePreCompute )
{
  if(m_max_num_newton == 0) return stepSceneLinear(scene, dt, updatePreCompute);
  else return stepSceneNonlinear(scene, dt, updatePreCompute);
}

template<int DIM>
std::string StrandCompliantManager<DIM>::getName() const
{
  return "Linear Compliant Implicit Euler With Preconditioner";
}

// explicit instantiations at bottom
template class StrandCompliantManager<2>;
template class StrandCompliantManager<3>;
