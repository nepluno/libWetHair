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


#ifndef __WET_HAIR_CORE_H__
#define __WET_HAIR_CORE_H__

#include "TwoDScene.h"
#include "SceneStepper.h"

template<int DIM>
class WetHairCore
{
public:
  
  WetHairCore( TwoDScene<DIM>* scene, SceneStepper<DIM>* scene_stepper, std::function<bool(double)>&& script_callback);
  ~WetHairCore();
  /////////////////////////////////////////////////////////////////////////////
  // Simulation Control Functions
  
  virtual void stepSystem();

  /////////////////////////////////////////////////////////////////////////////
  // Status Functions
  
  virtual void getBoundingBox(Vectors<DIM>& bb_min, Vectors<DIM>& bb_max);

  virtual const std::vector<scalar>& getTimingStatistics() const;
  virtual const std::vector<scalar>& getStepperTimingStatistics() const;
  
  virtual TwoDScene<DIM>* getScene() const;
  
  virtual const scalar& getCurrentTime() const;
  
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
private:
  TwoDScene<DIM>* m_scene;
  
  SceneStepper<DIM>* m_scene_stepper;
  
  std::vector<scalar> m_timing_statistics;
  std::function<bool(double)> m_script_callback;
  
  scalar m_current_time;
};

#endif
