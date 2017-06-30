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

#ifndef __STRAND_COMPLIANT_EULER__
#define __STRAND_COMPLIANT_EULER__

#include <iostream>
//#include <ceres/ceres.h>

#include "SceneStepper.h"
#include "MathUtilities.h"

template<int DIM>
class StrandCompliantManager;

template<int DIM>
class StrandCompliantEuler
{
public:
  StrandCompliantEuler(StrandCompliantManager<DIM>* parent, int hidx);
  
  virtual ~StrandCompliantEuler();
  
  virtual bool stepScene( TwoDScene<DIM> & scene, scalar dt );
  
  virtual bool stepScene( TwoDScene<DIM> & scene, scalar dt, VectorXs& r, const VectorXs& b );
  
  virtual bool PreconditionScene( TwoDScene<DIM> & scene, scalar dt, VectorXs& r, const VectorXs& b );
  
  virtual void preIterate( TwoDScene<DIM>& scene, scalar dt );
  
  virtual void computeRHS( TwoDScene<DIM> & scene, scalar dt, VectorXs& b );
  
  virtual void computeRHSIncremental( TwoDScene<DIM> & scene, scalar dt, VectorXs& b, const VectorXs& vplus );
  
  virtual void computeAp( const VectorXs& p, VectorXs& b );
  
  virtual void updateLambda( TwoDScene<DIM>& scene, const VectorXs& dx, const VectorXs& dv, scalar dt );
  
  virtual void updateNextV( TwoDScene<DIM>& scene, const VectorXs& vplus );
private:
  TripletXs m_A_nz;
  TripletXs m_J_nz;
  TripletXs m_Jv_nz;
  TripletXs m_Jxv_nz;
  TripletXs m_invC_nz;
  TripletXs m_invCv_nz;
  
  TripletXs m_J_inter_nz;
  TripletXs m_Jv_inter_nz;
  TripletXs m_invC_inter_nz;
  TripletXs m_invCv_inter_nz;
  
  SparseXs m_A;
  SparseXs m_J;
  SparseXs m_Jv;
  SparseXs m_Jxv;
  SparseXs m_invC;
  SparseXs m_invCv;
  
  SparseXs m_J_inter;
  SparseXs m_Jv_inter;
  SparseXs m_invC_inter;
  SparseXs m_invCv_inter;
  
  VectorXs m_A_inv_diag;
  
  Eigen::SimplicialLDLT<SparseXs> m_solver;
  
  int m_hidx;

  
  StrandCompliantManager<DIM>* m_parent;
  
  int m_start_global_dof;
  int m_num_global_dof;
};

#endif
