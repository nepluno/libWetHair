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


#ifndef __TWO_D_SCENE_SERIALIZER_H__
#define __TWO_D_SCENE_SERIALIZER_H__

#include <fstream>
#include <iostream>

#include "StringUtilities.h"

template<int DIM>
class TwoDScene;

template<int DIM>
class SceneStepper;

template<int DIM>
class TwoDSceneRenderer;

template<int DIM>
class TwoDSceneSerializer
{
public:
  void serializeFluidReadable( TwoDScene<DIM>& scene, std::vector< std::ostringstream >& outputstream ) const;
  
  void serializeHairReadable( TwoDScene<DIM>& scene, std::ostream& outputstream ) const;
  
  void serializeShallowFlowReadable( const TwoDSceneRenderer<DIM>* renderer, TwoDScene<DIM>& scene, std::ostream& outputstream ) const;
  
  void serializeBoundariesReadable( const TwoDSceneRenderer<DIM>* renderer, TwoDScene<DIM>& scene, std::ostream& os_boundary_single, std::ostream& os_boundary_double);
  
  void serializePolygonalCohesionReadable( const TwoDSceneRenderer<DIM>* renderer, TwoDScene<DIM>& scene, std::ostream& os_pe, std::ostream& os_poe, std::ostream& os_ppp ) const;
  
  bool deSerializeFluidReadable( TwoDScene<DIM>& scene, const std::vector< std::string > &filename_fluids );

  bool deSerializeHairReadable( TwoDScene<DIM>& scene, const std::string &filename_hairs );

  bool deSerializeShallowFlowReadable( const TwoDSceneRenderer<DIM>* renderer, TwoDScene<DIM>& scene, const std::string &filename_flows );
  
  void serializeScene( TwoDScene<DIM>& scene, SceneStepper<DIM>* stepper, std::ostream& outputstream ) const;

  void loadScene( TwoDScene<DIM>& scene, SceneStepper<DIM>* stepper, std::istream& inputstream ) const;
};

#endif
