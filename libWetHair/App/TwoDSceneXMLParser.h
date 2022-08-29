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

#ifndef LIBWETHAIR_APP_TWO_D_SCENE_XML_PARSER_H_
#define LIBWETHAIR_APP_TWO_D_SCENE_XML_PARSER_H_

#include <Eigen/StdVector>
#include <fstream>
#include <iostream>
#include <limits>

#include "Camera.h"
#include "CompliantImplicitEuler.h"
#include "CylindricalShallowFlow.h"
#include "DER/StrandForce.h"
#include "DER/StrandParameters.h"
#include "DragDampingForce.h"
#include "ExecutableSimulation.h"
#include "LevelSetForce.h"
#include "LinearBendingForce.h"
#include "LinearSpringForce.h"
#include "ParticleSimulation.h"
#include "PolygonalCohesion.h"
#include "RenderingUtilities.h"
#include "SimpleGravityForce.h"
#include "SpringForce.h"
#include "StrandCompliantManager.h"
#include "StringUtilities.h"
#include "TwoDScene.h"
#include "TwoDSceneRenderer.h"
#include "TwoDimensionalDisplayController.h"
#include "fluidsim2D.h"
#include "fluidsim3D.h"
#include "rapidxml.hpp"

// REALLY USEFULL TODOs
//   TODO: Improve error messages to display all valid options, etc. Could
//   define an option class that knows its valid options and bounds on values.

// LESS USEFULL TODOs
//   TODO: Write method for computing number of a given property
//   TODO: Add some additional error checking for repeat properties, etc
//   TODO: Abstract out common code
//   TODO: Check for invalid properties

class TwoDSceneXMLParser {
 public:
  void loadExecutableSimulation(
      const std::string& file_name,
      const char* memory_str,
      const std::string& init_filename,
      bool simulate_comparison,
      bool rendering_enabled,
      ExecutableSimulation** execsim,
      libwethair::renderingutils::Viewport& view, Camera& cam,
      libwethair::scalar& max_time,
      libwethair::scalar& steps_per_sec_cap,
      libwethair::renderingutils::Color& bgcolor,
      std::string& description,
      std::string& scenetag,
      bool& cam_inited,
      bool& view_inited);

  // TODO: NEED AN EIGEN_ALIGNED_THING_HERE ?
 private:
  template <int DIM>
  void loadParticleSimulation(bool simulate_comparison, bool rendering_enabled,
                              ExecutableSimulation** execsim,
                              libwethair::renderingutils::Viewport& view,
                              libwethair::renderingutils::Color& bgcolor,
                              rapidxml::xml_node<>* node,
                              const std::string& init_filename);

  void loadDERSimulation(bool simulate_comparison, bool rendering_enabled,
                         ExecutableSimulation** execsim,
                         libwethair::renderingutils::Viewport& view,
                         libwethair::renderingutils::Color& bgcolor,
                         rapidxml::xml_node<>* node,
                         const std::string& init_filename);

  void loadXMLFile(const std::string& filename, std::vector<char>& xmlchars,
                   rapidxml::xml_document<>& doc);

  void loadXMLFileFromMemory(const char* memory_str,
                             std::vector<char>& xmlchars,
                             rapidxml::xml_document<>& doc);

  bool loadTextFileIntoString(const std::string& filename,
                              std::string& filecontents);

  void loadSimulationType(rapidxml::xml_node<>* node, std::string& simtype);

  template <int DIM>
  void loadLiquidParameter(rapidxml::xml_node<>* node,
                           libwethair::TwoDScene<DIM>& twoscene);

  template <int DIM>
  void loadFluidSim(rapidxml::xml_node<>* node, libwethair::TwoDScene<DIM>& twoscene);

  template <int DIM>
  void loadLinearBendingForces(rapidxml::xml_node<>* node,
                               libwethair::TwoDScene<DIM>& twodscene);

  template <int DIM>
  void loadParticleEdges(rapidxml::xml_node<>* node, libwethair::TwoDScene<DIM>& twodscene);

  void loadStrandParticleEdges(rapidxml::xml_node<>* node, libwethair::TwoDScene<3>& scene);

  template <int DIM>
  void loadFlow(rapidxml::xml_node<>* node, libwethair::TwoDScene<DIM>& twodscene);

  void loadSceneTag(rapidxml::xml_node<>* node, std::string& scenetag);

  template <int DIM>
  void loadSpringForces(rapidxml::xml_node<>* node, libwethair::TwoDScene<DIM>& twodscene);

  void loadScripts(rapidxml::xml_node<>* node, std::vector<libwethair::Script>& scripts);

  template <int DIM>
  void loadLinearSpringForces(rapidxml::xml_node<>* node,
                              libwethair::TwoDScene<DIM>& twodscene);

  template <int DIM>
  void loadDragDampingForces(rapidxml::xml_node<>* node,
                             libwethair::TwoDScene<DIM>& twodscene);

  template <int DIM>
  void loadIntegrator(rapidxml::xml_node<>* node, libwethair::TwoDScene<DIM>& twodscene,
                      libwethair::SceneStepper<DIM>** scenestepper);

  template <int DIM>
  void loadStrandParameters(rapidxml::xml_node<>* node, libwethair::TwoDScene<DIM>& scene);

  template <int DIM>
  void loadStrandForces(rapidxml::xml_node<>* node, libwethair::TwoDScene<DIM>& scene);

  void loadMaxTime(rapidxml::xml_node<>* node, libwethair::scalar& max_t);

  void loadMaxSimFrequency(rapidxml::xml_node<>* node, libwethair::scalar& max_freq);

  bool loadViewport(rapidxml::xml_node<>* node, libwethair::renderingutils::Viewport& view);

  void loadBackgroundColor(rapidxml::xml_node<>* node,
                           libwethair::renderingutils::Color& color);

  void loadParticleColors(rapidxml::xml_node<>* node,
                          std::vector<libwethair::renderingutils::Color>& particle_colors);

  void loadEdgeColors(rapidxml::xml_node<>* node,
                      std::vector<libwethair::renderingutils::Color>& edge_colors);

  void loadSceneDescriptionString(rapidxml::xml_node<>* node,
                                  std::string& description_string);

  bool loadCamera(rapidxml::xml_node<>* node, Camera& camera);
};

#endif  // LIBWETHAIR_APP_TWO_D_SCENE_XML_PARSER_H_
