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

#ifndef __LINEAR_BENDING_FORCE_H__
#define __LINEAR_BENDING_FORCE_H__

#include <Eigen/Core>
#include "Force.h"
#include <iostream>

template<int DIM>
class TwoDScene;

template<int DIM>
class LinearBendingForce : public Force
{
public:
  
  LinearBendingForce( TwoDScene<DIM>* parent, int idx1, int idx2, int idx3, const scalar& alpha, const scalar& beta, const Vectors<DIM-1>& theta0, const scalar& eb1n, const scalar& eb2n );
  
  virtual ~LinearBendingForce();
    
  virtual Force* createNewCopy();
  
  virtual void preCompute( const VectorXs& x, const VectorXs& v, const VectorXs& m, const scalar& dt );

  virtual void getAffectedVars( int pidx, std::unordered_set<int>& vars );
  
  virtual int getAffectedHair( const std::vector<int> particle_to_hairs );
  
  virtual bool isContained( int pidx );
  
  virtual const char* name();
  
  virtual void computeIntegrationVars( const VectorXs& x, const VectorXs& v, const VectorXs& m,
                                      VectorXs& lambda, VectorXs& lambda_v,
                                      TripletXs& J, TripletXs& Jv, TripletXs& Jxv, TripletXs& tildeK,
                                      TripletXs& stiffness, TripletXs& damping, VectorXs& Phi, VectorXs& Phiv, const scalar& dt);
  
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
  
  Matrixs<DIM> m_R;       // rotation matrix
  
  scalar m_alpha;     // stiffness coefficient
  scalar m_beta;      // damping coefficient
  Vectors<DIM-1> m_theta0;    // rest angle
  scalar m_eb1n;      // norm of e1 bar
  scalar m_eb2n;      // norm of e2 bar
  
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

#endif
