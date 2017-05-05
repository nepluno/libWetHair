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

#ifndef __EXECUTABLE_SIMULATION_H__
#define __EXECUTABLE_SIMULATION_H__

#include <string>
#include <Eigen/Core>
#include "MathDefs.h"
#include "RenderingUtilities.h"
#include "Camera.h"

template<int DIM>
class TwoDSceneRenderer;

class ExecutableSimulation
{
public:

  virtual ~ExecutableSimulation()
  {}

  /////////////////////////////////////////////////////////////////////////////
  // Simulation Control Functions
  
  virtual void stepSystem() = 0;
  
  /////////////////////////////////////////////////////////////////////////////
  // OpenGL Rendering Functions
  
  // TODO: Combine both of these
  virtual void renderSceneOpenGL() = 0;
  
  // TODO: Scrap these two functions, they are no longer needed
  virtual void renderSceneDifferencesOpenGL() = 0;
  virtual void initializeOpenGLRenderer() = 0;
  
  virtual void updateOpenGLRendererState() = 0;

  virtual void computeCameraCenter( renderingutils::Viewport& view ) = 0;

  /////////////////////////////////////////////////////////////////////////////
  // Serialization Functions
  
  virtual void serializeScene( std::ostream& outputstream ) = 0;
  virtual void serializeSceneReadable( std::vector< std::ostringstream >& osfluid, std::ostream& oshair, std::ostream& osflow, std::ostream& os_boundary_single, std::ostream& os_boundary_double, std::ostream& os_pe, std::ostream& os_poe, std::ostream& os_ppp ) = 0;
  virtual bool deSerializeSceneReadable( const std::vector< std::string > &filename_fluids, const std::string &filename_hairs, const std::string &filename_flows, const std::string &filename_bd_single, const std::string &filename_bd_double ) = 0;
  virtual void loadSerializedScene( std::istream& inputstream ) = 0;
  
  /////////////////////////////////////////////////////////////////////////////
  // Status Functions

  virtual std::string getSolverName() = 0;

  virtual void centerCamera(bool b_reshape = true) = 0;
  virtual void keyboard( unsigned char key, int x, int y ) = 0;
  virtual void reshape( int w, int h ) = 0;
  
  virtual void special( int key, int x, int y ) = 0;
  
  virtual void mouse( int button, int state, int x, int y ) = 0;
  
  virtual void translateView( double dx, double dy ) = 0;
  
  virtual void zoomView( double dx, double dy ) = 0;
  
  virtual void motion( int x, int y ) = 0;
  
  virtual const scalar& getDt() const = 0;
  virtual const scalar& getCurrentTime() const = 0;
  virtual scalar getCurrentFrame() const = 0;
  
  virtual int getWindowWidth() const = 0;
  virtual int getWindowHeight() const = 0;
  
  virtual void setWindowWidth(int w) = 0;
  virtual void setWindowHeight(int h) = 0;
  
  virtual void setCenterX( double x ) = 0;
  virtual void setCenterY( double y ) = 0;
  virtual void setScaleFactor( double scale ) = 0;
  
  virtual void setCamera( const Camera& cam ) = 0;
  virtual void setView( const renderingutils::Viewport& view ) = 0;
  
  virtual const std::vector<scalar>& getStepperTimingStatistics() const = 0;
  
  virtual const std::vector<scalar>& getTimingStatistics() const = 0;
};

#endif
