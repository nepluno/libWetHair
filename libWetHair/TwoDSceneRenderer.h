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


#ifndef __TWO_D_SCENE_RENDERER_H__
#define __TWO_D_SCENE_RENDERER_H__

#include <Eigen/StdVector>
#include <iostream>
#include "TwoDScene.h"
#include "MathUtilities.h"
#include "RenderingUtilities.h"
#include "HairFlow.h"

#include "Icosphere.h"
#include "Capsule.h"
#include "RoundCornerBox.h"

// TODO: Get display controller out of here
// TODO: Make particle system and rigid body renderers that inherit from this

template<int DIM>
class TwoDimensionalDisplayController;

template<int DIM>
class TwoDSceneRenderer
{
public:
  enum PARTICLE_VIS_MODE
  {
    PVM_NONE,
    PVM_ETA,
    PVM_AREA,
    PVM_KAPPA,
    
    PVM_COUNT
  };
  
  enum EDGE_VIS_MODE
  {
    EVM_NONE,
    EVM_AREA,
    EVM_KAPPA,
    
    EVM_COUNT
  };
  
  // TODO: Gut this method
  TwoDSceneRenderer( const TwoDScene<DIM>& scene,
                     const std::vector<renderingutils::Color>& particle_colors, const std::vector<renderingutils::Color>& edge_colors);
  
  TwoDSceneRenderer();
  
  
  void initializeOpenGLRenderer(const TwoDScene<DIM>& scene);
  void updateOpenGLRenderer(const TwoDScene<DIM>& scene, bool updateDevice);
  void renderParticleSimulation( const TwoDScene<DIM>& scene ) const;
  
  void renderVolumeGraph( const TwoDScene<DIM>& scene ) const;

  // Returns a reference to the vector containing particle colors
  std::vector<renderingutils::Color>& getParticleColors();
  const std::vector<renderingutils::Color>& getParticleColors() const;
  
  void selectNextParticleVisMode();
  
  void selectNextEdgeVisMode();
  
  void switchShowParticleNormal();
  
  void switchShowEdgeNormal();
  
  void switchShowLiquidPolygon();
  
  void setDC(const TwoDimensionalDisplayController<DIM>*);
  
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  
  void initializeCircleRenderer( int num_points );
  void initializeSemiCircleRenderer( int num_points );
  void initializeCylinderRenderer( int num_points, const TwoDScene<DIM>& scene );
  void initializeBoundaryRenderer( const TwoDScene<DIM>& scene );
  
  void renderSolidCircle( const Vectors<DIM>& center, double radius ) const;

  void writeTransformedHairFlow( std::ostream& o, const TwoDScene<DIM>& scene ) const;
  
  void writeBoundaries(std::ostream& os_single, std::ostream& os_double) const;
  
  int renderDebuggingInfo( const HairFlow<DIM> &flow ) const;
  
  void renderAnalyticalShape( const TwoDScene<DIM>& scene, int layer ) const;
  
  const TwoDimensionalDisplayController<DIM>* m_dc;
  
  // TODO: Move this out of here and into some subclass
  // Particle System rendering state
  std::vector<renderingutils::Color> m_particle_colors;
  std::vector<renderingutils::Color> m_edge_colors;
  
  // Precomputed points for a circle
  std::vector<std::pair<double,double> > m_circle_points;
  std::vector<std::pair<double,double> > m_semi_circle_points;
  
  IcosphereCreator m_icosphere;
  std::vector< CapsuleCreator > m_capsules;
  std::vector< RoundCornerBox > m_roundcornerboxes;
  
  std::vector<Vector3s> m_cylinder_points;
  std::vector<int> m_cylinder_elements;
  
  bool m_draw_grid;
  bool m_draw_particles;
  bool m_draw_velocities;
  bool m_draw_boundaries;
  
  bool m_show_particle_normal;
  bool m_show_edge_normal;
  bool m_show_liquid_polygon;
  int m_pvm;
  int m_evm;
  
  GLuint m_element_hairs;
  GLuint m_vertex_hair_core;
  GLuint m_vertex_hair_flow;
  GLuint m_vertex_fluids;
  
  const TwoDScene<DIM>* m_scene;
};

#endif
