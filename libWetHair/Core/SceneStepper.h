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

#ifndef __SCENE_STEPPER__
#define __SCENE_STEPPER__

#include "TwoDScene.h"

#include "MathDefs.h"

template<int DIM>
class SceneStepper
{
protected:
  std::vector<scalar> m_timing_statistics;
  
  VectorXs m_old_v;
  VectorXs m_a;
  VectorXs m_next_x;
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  
  virtual ~SceneStepper();
  
  virtual bool stepScene( TwoDScene<DIM> & scene, scalar dt, bool updatePreCompute = true ) = 0;
  
  virtual void setNextX( const VectorXs& nextx );
  
  virtual const VectorXs& getNextX() const;
  
  virtual void accept( TwoDScene<DIM> & scene, scalar dt );
  
  virtual std::string getName() const = 0;
  
  virtual const VectorXs& getAcceleration() const;
  
  virtual void PostStepScene( TwoDScene<DIM> & scene, scalar dt );
  
  virtual const std::vector<scalar>& getTimingStatistics() const;
  
  virtual void write(std::vector<scalar>&) const;
  
  virtual void read(const scalar* data);
  
  virtual size_t size();
  
  virtual void init( TwoDScene<DIM> & scene );
};

#endif
