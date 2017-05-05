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

#ifndef __COMPLIANT_IMPLICIT_EULER__
#define __COMPLIANT_IMPLICIT_EULER__

#include <Eigen/Dense>
#include <iostream>

#include "SceneStepper.h"
#include "MathUtilities.h"
#include "StringUtilities.h"

template<int DIM>
class TwoDScene;

template<int DIM>
class CompliantImplicitEuler : public SceneStepper<DIM>
{
public:
  CompliantImplicitEuler(TwoDScene<DIM>* scene, int max_iters, scalar criterion, bool autoUpdateNextInit = true);
  
  virtual ~CompliantImplicitEuler();
  
  virtual bool stepScene( TwoDScene<DIM>& scene, scalar dt, bool updatePreCompute = true );
  
  virtual std::string getName() const;
private:
  void zeroFixedDoFs( const TwoDScene<DIM>& scene, VectorXs& vec );
  void updateNumConstraints(const VectorXs& dx, const VectorXs& dv, scalar dt);
  
  bool m_bAutoUpdateNextInit;
  int m_max_iters;
  scalar m_criterion;
  
  TwoDScene<DIM>* m_scene;
  
  TripletXs m_Kext_nz;
  TripletXs m_A_nz;
  TripletXs m_J_nz;
  TripletXs m_Jv_nz;
  TripletXs m_Jxv_nz;
  TripletXs m_Fix_nz;
  TripletXs m_M_nz;
  TripletXs m_invC_nz;
  TripletXs m_invCv_nz;
  
  SparseXs m_Kext;
  SparseXs m_A;
  SparseXs m_J;
  SparseXs m_JC;
  SparseXs m_Jv;
  SparseXs m_JvC;
  SparseXs m_Jxv;
  SparseXs m_Fix;
  SparseXs m_M;
  SparseXs m_invC;
  SparseXs m_invCv;
  
  VectorXs m_lambda;
  VectorXs m_lambda_v;
  VectorXs m_gradU;
  VectorXs m_Phi;
  VectorXs m_Phi_v;
  VectorXs m_b;
  VectorXs m_vplus;
  
  Vector6i m_interhair_idx;
  Vector6i m_interhair_num;

  Eigen::SimplicialLDLT< SparseXs > m_solver;
  Eigen::ConjugateGradient< SparseXs > m_iterative_solver;
};

#endif
