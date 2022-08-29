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

#ifndef LIBWETHAIR_APP_TWO_D_SCENE_RENDERER_H_
#define LIBWETHAIR_APP_TWO_D_SCENE_RENDERER_H_

#include <Eigen/StdVector>
#include <iostream>

#include "Capsule.h"
#include "HairFlow.h"
#include "Icosphere.h"
#include "MathUtilities.h"
#include "RenderingUtilities.h"
#include "RoundCornerBox.h"
#include "TwoDScene.h"

// TODO: Get display controller out of here
// TODO: Make particle system and rigid body renderers that inherit from this

template <int DIM>
class TwoDimensionalDisplayController;

template <int DIM>
class TwoDSceneRenderer {
 public:
  enum PARTICLE_VIS_MODE {
    PVM_NONE,
    PVM_ETA,
    PVM_AREA,
    PVM_KAPPA,

    PVM_COUNT
  };

  enum EDGE_VIS_MODE {
    EVM_NONE,
    EVM_AREA,
    EVM_KAPPA,

    EVM_COUNT
  };

  // TODO: Gut this method
  TwoDSceneRenderer(const libwethair::TwoDScene<DIM>& scene,
                    const std::vector<libwethair::renderingutils::Color>& particle_colors,
                    const std::vector<libwethair::renderingutils::Color>& edge_colors);

  TwoDSceneRenderer();

  void initializeOpenGLRenderer(const libwethair::TwoDScene<DIM>& scene);
  void updateOpenGLRenderer(const libwethair::TwoDScene<DIM>& scene, bool updateDevice);
  void renderParticleSimulation(const libwethair::TwoDScene<DIM>& scene) const;

  void renderVolumeGraph(const libwethair::TwoDScene<DIM>& scene) const;

  // Returns a reference to the vector containing particle colors
  std::vector<libwethair::renderingutils::Color>& getParticleColors();
  const std::vector<libwethair::renderingutils::Color>& getParticleColors() const;

  void selectNextParticleVisMode();

  void selectNextEdgeVisMode();

  void switchShowParticleNormal();

  void switchShowEdgeNormal();

  void switchShowLiquidPolygon();

  void setDC(const TwoDimensionalDisplayController<DIM>*);

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  void initializeCircleRenderer(int num_points);
  void initializeSemiCircleRenderer(int num_points);
  void initializeCylinderRenderer(int num_points, const libwethair::TwoDScene<DIM>& scene);
  void initializeBoundaryRenderer(const libwethair::TwoDScene<DIM>& scene);

  void renderSolidCircle(const libwethair::Vectors<DIM>& center, double radius) const;

  void writeTransformedHairFlow(std::ostream& o,
                                const libwethair::TwoDScene<DIM>& scene) const;

  void writeBoundaries(std::ostream& os_single, std::ostream& os_double) const;

  int renderDebuggingInfo(const libwethair::HairFlow<DIM>& flow) const;

  void renderAnalyticalShape(const libwethair::TwoDScene<DIM>& scene, int layer) const;

  const TwoDimensionalDisplayController<DIM>* m_dc;

  // TODO: Move this out of here and into some subclass
  // Particle System rendering state
  std::vector<libwethair::renderingutils::Color> m_particle_colors;
  std::vector<libwethair::renderingutils::Color> m_edge_colors;

  // Precomputed points for a circle
  std::vector<std::pair<double, double> > m_circle_points;
  std::vector<std::pair<double, double> > m_semi_circle_points;

  libwethair::IcosphereCreator m_icosphere;
  std::vector<libwethair::CapsuleCreator> m_capsules;
  std::vector<libwethair::RoundCornerBox> m_roundcornerboxes;

  std::vector<libwethair::Vector3s> m_cylinder_points;
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

  const libwethair::TwoDScene<DIM>* m_scene;
};

#endif  // LIBWETHAIR_APP_TWO_D_SCENE_RENDERER_H_
