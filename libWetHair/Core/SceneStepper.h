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
