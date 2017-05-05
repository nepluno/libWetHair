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


#ifndef __STRAND_COMPLIANT_MANAGER__
#define __STRAND_COMPLIANT_MANAGER__

#include "TwoDScene.h"

#include "MathDefs.h"

#include "SceneStepper.h"

#include <vector>

#include <future>

template<int DIM>
class StrandCompliantEuler;

template<int DIM>
class StrandCompliantManager : public SceneStepper<DIM>
{
public:
  StrandCompliantManager(TwoDScene<DIM>* scene, int max_newton, int max_iters, scalar criterion, bool compute_interhair, bool use_preconditioner);
  
  virtual ~StrandCompliantManager();
  
  virtual bool stepScene( TwoDScene<DIM> & scene, scalar dt, bool updatePreCompute = true );
  
  virtual bool stepSceneNonlinear( TwoDScene<DIM> & scene, scalar dt, bool updatePreCompute = true );
  
  virtual bool stepSceneLinear( TwoDScene<DIM> & scene, scalar dt, bool updatePreCompute = true );
  
  virtual std::string getName() const;
  
  virtual void zeroFixedDoFs( const TwoDScene<DIM>& scene, VectorXs& vec );
  
  virtual void localPreIterate(TwoDScene<DIM> & scene, scalar dt);
  
  virtual void localStepScene(TwoDScene<DIM> & scene, scalar dt, VectorXs& r, const VectorXs& b);
  
  virtual void localPreconditionScene(TwoDScene<DIM> & scene, scalar dt, VectorXs& r, const VectorXs& b);
  
  virtual void localUpdateLambda(TwoDScene<DIM> & scene, const VectorXs& dx, const VectorXs& dv, scalar dt);
  
  virtual void updateNextV(TwoDScene<DIM> & scene, const VectorXs& vplus);
  
  virtual void localUpdateNumConstraints(const VectorXs& dx, const VectorXs& dv, scalar dt);
  
  virtual void interhairUpdateNumConstraints(const VectorXs& dx, const VectorXs& dv, scalar dt);
  
  virtual void computeRHSLocal( TwoDScene<DIM> & scene, scalar dt, VectorXs& b );
  
  virtual void computeRHSIncrementalLocal( TwoDScene<DIM> & scene, scalar dt, VectorXs& b, const VectorXs& vplus );
  
  virtual void computeRHSIncrementalInterhair( TwoDScene<DIM> & scene, scalar dt, VectorXs& b, const VectorXs& vplus );
  
  virtual void computeAp( TwoDScene<DIM> & scene, const VectorXs& p, VectorXs& b, scalar dt );
  
  virtual void write(std::vector<scalar>&) const;
  
  virtual void read(const scalar* data);
  
  virtual size_t size();
  
protected:
  std::vector<StrandCompliantEuler<DIM>*> m_integrators;
  
  scalar m_dt;
  TwoDScene<DIM>* m_scene;
  
  TripletXs m_A_nz;
  TripletXs m_J_nz;
  TripletXs m_Jv_nz;
  TripletXs m_Jxv_nz;
  TripletXs m_Fix_nz;
  TripletXs m_M_nz;
  TripletXs m_invC_nz;
  TripletXs m_invCv_nz;
  
  TripletXs m_J_interhair_nz;
  TripletXs m_Jv_interhair_nz;
  SparseXs m_J_interhair;
  SparseXs m_Jv_interhair;
  SparseXs m_JT_interhair;
  SparseXs m_JvT_interhair;
  VectorXs m_invC_interhair;
  VectorXs m_invCv_interhair;
  
  VectorXs m_lambda;
  VectorXs m_lambda_v;
  VectorXs m_gradU;
  VectorXs m_Phi;
  VectorXs m_Phi_interhair;
  VectorXs m_Phi_v;
  VectorXs m_Phiv_interhair;
  VectorXs m_vplus;
  
  VectorXs m_dvi;
  VectorXs m_dxi;
  VectorXs m_dv;
  VectorXs m_dx;
  VectorXs m_dx_scripted;
  
  VectorXs m_v;
  VectorXs m_r;
  VectorXs m_p;
  VectorXs m_q;
  VectorXs m_t;
  VectorXs m_z;
  VectorXs m_rhs;
  
  VectorXs m_cv_buffer;
  VectorXs m_c_buffer;
  
  Vector6i m_interhair_idx;
  Vector6i m_interhair_num;
  
  int m_max_num_newton;
  int m_max_num_iters;
  scalar m_criterion;
  
  bool m_compute_interhair;
  bool m_use_preconditioner;
  
  friend class StrandCompliantEuler<DIM>;
};

#endif
