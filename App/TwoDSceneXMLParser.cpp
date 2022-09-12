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

#include "TwoDSceneXMLParser.h"

#include <libWetHair/DER/Dependencies/ElasticStrandUtils.h>
#include <libWetHair/ThreadUtils.h>

using namespace libwethair;

void TwoDSceneXMLParser::loadExecutableSimulation(
    GLFWwindow* window,
    const std::string& file_name, const char* memory_str,
    const std::string& init_filename, bool simulate_comparison,
    bool rendering_enabled, ExecutableSimulation** execsim,
    renderingutils::Viewport& view, Camera& cam, scalar& max_time,
    scalar& steps_per_sec_cap, renderingutils::Color& bgcolor,
    std::string& description, std::string& scenetag, bool& cam_inited,
    bool& view_inited) {
  // Load the xml document
  std::vector<char> xmlchars;
  rapidxml::xml_document<> doc;

  if (memory_str) {
    loadXMLFileFromMemory(memory_str, xmlchars, doc);
  } else {
    loadXMLFile(file_name, xmlchars, doc);
  }

  // Attempt to locate the root node
  rapidxml::xml_node<>* node = doc.first_node("scene");
  if (node == NULL) {
    std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse xml "
                 "scene file. Failed to locate root <scene> node. Exiting."
              << std::endl;
    exit(1);
  }

  int dimension = 2;
  if (node->first_attribute("dim")) {
    std::string attribute(node->first_attribute("dim")->value());
    if (!stringutils::extractFromString(attribute, dimension) ||
        (dimension != 2 && dimension != 3)) {
      std::cerr
          << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse value "
             "of dimension attribute for scene. Value must be 2 or 3. Exiting."
          << std::endl;
      exit(1);
    }
  }
  // Determine what simulation type this is (particle, rigid body, etc)
  std::string simtype;
  loadSimulationType(node, simtype);

  // Parse common state
  loadMaxTime(node, max_time);
  loadMaxSimFrequency(node, steps_per_sec_cap);
  view_inited = loadViewport(node, view);
  loadBackgroundColor(node, bgcolor);
  loadSceneDescriptionString(node, description);
  loadSceneTag(node, scenetag);

  cam_inited = loadCamera(node, cam);

  // Parse the user-requested simulation type. The default is a particle
  // simulation.
  if (simtype == "" || simtype == "particle-system" ||
      simtype == "mass-spring") {
    switch (dimension) {
      case 2:
        loadParticleSimulation<2>(window, simulate_comparison, rendering_enabled,
                                  execsim, view, bgcolor, node, init_filename);
        break;
      case 3:
        loadParticleSimulation<3>(window, simulate_comparison, rendering_enabled,
                                  execsim, view, bgcolor, node, init_filename);
        break;
      default:
        std::cerr
            << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Invalid dimension '"
            << dimension
            << "' specified. Valid options are '2' and '3'. Exiting."
            << std::endl;
        exit(1);
        break;
    }
  } else if (simtype == "DiscreteElasticRods" && dimension == 3) {
    loadDERSimulation(window, simulate_comparison, rendering_enabled, execsim, view,
                      bgcolor, node, init_filename);
  } else {
    std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Invalid simtype '"
              << simtype
              << "' specified. Valid options is 'particle-system'. Exiting."
              << std::endl;
    exit(1);
  }
}

void TwoDSceneXMLParser::loadScripts(rapidxml::xml_node<>* node,
                                     std::vector<Script>& scripts) {
  for (rapidxml::xml_node<>* nd = node->first_node("script"); nd;
       nd = nd->next_sibling("script")) {
    Script scr;
    rapidxml::xml_attribute<>* typend = nd->first_attribute("target");
    if (typend) {
      std::string handlertype(typend->value());
      if (handlertype == "camera")
        scr.target = Script::CAMERA;
      else if (handlertype == "root")
        scr.target = Script::ROOT;
      else if (handlertype == "solid")
        scr.target = Script::SOLID;
      else if (handlertype == "source")
        scr.target = Script::SOURCE;
      else if (handlertype == "curlradius")
        scr.target = Script::CURLRADIUS;
      else if (handlertype == "curldensity")
        scr.target = Script::CURLDENSITY;
      else if (handlertype == "all")
        scr.target = Script::ALL;
      else {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Invalid script "
                     "'target' attribute specified. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    typend = nd->first_attribute("type");
    if (typend) {
      std::string handlertype(typend->value());
      if (handlertype == "rotate")
        scr.type = Script::ROTATE;
      else if (handlertype == "translate")
        scr.type = Script::TRANSLATE;
      else if (handlertype == "scale")
        scr.type = Script::SCALE;
      else if (handlertype == "absorption")
        scr.type = Script::ABSORB;
      else {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Invalid script "
                     "'type' attribute specified. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    scr.func = Script::CUBIC;
    typend = nd->first_attribute("func");
    if (typend) {
      std::string handlertype(typend->value());
      if (handlertype == "cubic")
        scr.func = Script::CUBIC;
      else if (handlertype == "cosine")
        scr.func = Script::COSINE;
      else if (handlertype == "weno")
        scr.func = Script::WENO;
      else {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Invalid script "
                     "'func' attribute specified. Exiting."
                  << std::endl;
        exit(1);
      }
    }
    scr.base_pos = 0.0;

    if (nd->first_attribute("x")) {
      std::string attribute(nd->first_attribute("x")->value());
      if (!stringutils::extractFromString(attribute, scr.v(0))) {
        std::cerr
            << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse value "
               "of x attribute for script. Value must be scalar. Exiting."
            << std::endl;
        exit(1);
      }
    }
    if (nd->first_attribute("y")) {
      std::string attribute(nd->first_attribute("y")->value());
      if (!stringutils::extractFromString(attribute, scr.v(1))) {
        std::cerr
            << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse value "
               "of y attribute for script. Value must be scalar. Exiting."
            << std::endl;
        exit(1);
      }
    }
    if (nd->first_attribute("z")) {
      std::string attribute(nd->first_attribute("z")->value());
      if (!stringutils::extractFromString(attribute, scr.v(2))) {
        std::cerr
            << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse value "
               "of z attribute for script. Value must be scalar. Exiting."
            << std::endl;
        exit(1);
      }
    }
    if (nd->first_attribute("w")) {
      std::string attribute(nd->first_attribute("w")->value());
      if (!stringutils::extractFromString(attribute, scr.v(3))) {
        std::cerr
            << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse value "
               "of w attribute for script. Value must be scalar. Exiting."
            << std::endl;
        exit(1);
      }
    }

    if (nd->first_attribute("ox")) {
      std::string attribute(nd->first_attribute("ox")->value());
      if (!stringutils::extractFromString(attribute, scr.origin(0))) {
        std::cerr
            << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse value "
               "of ox attribute for script. Value must be scalar. Exiting."
            << std::endl;
        exit(1);
      }
    }
    if (nd->first_attribute("oy")) {
      std::string attribute(nd->first_attribute("oy")->value());
      if (!stringutils::extractFromString(attribute, scr.origin(1))) {
        std::cerr
            << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse value "
               "of oy attribute for script. Value must be scalar. Exiting."
            << std::endl;
        exit(1);
      }
    }
    if (nd->first_attribute("oz")) {
      std::string attribute(nd->first_attribute("oz")->value());
      if (!stringutils::extractFromString(attribute, scr.origin(2))) {
        std::cerr
            << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse value "
               "of oz attribute for script. Value must be scalar. Exiting."
            << std::endl;
        exit(1);
      }
    }

    if (nd->first_attribute("start")) {
      std::string attribute(nd->first_attribute("start")->value());
      if (!stringutils::extractFromString(attribute, scr.start)) {
        std::cerr
            << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse value "
               "of start attribute for script. Value must be scalar. Exiting."
            << std::endl;
        exit(1);
      }
    }
    if (nd->first_attribute("end")) {
      std::string attribute(nd->first_attribute("end")->value());
      if (!stringutils::extractFromString(attribute, scr.end)) {
        std::cerr
            << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse value "
               "of end attribute for script. Value must be scalar. Exiting."
            << std::endl;
        exit(1);
      }
    }

    scr.ease_start = scr.ease_end = (scr.end - scr.start) / 3.0;

    if (nd->first_attribute("easestart")) {
      std::string attribute(nd->first_attribute("easestart")->value());
      if (!stringutils::extractFromString(attribute, scr.ease_start)) {
        std::cerr
            << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse value "
               "of start attribute for script. Value must be scalar. Exiting."
            << std::endl;
        exit(1);
      }
    }
    if (nd->first_attribute("easeend")) {
      std::string attribute(nd->first_attribute("easeend")->value());
      if (!stringutils::extractFromString(attribute, scr.ease_end)) {
        std::cerr
            << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse value "
               "of start attribute for script. Value must be scalar. Exiting."
            << std::endl;
        exit(1);
      }
    }

    scr.amplitude = 1.0;
    if (nd->first_attribute("amplitude")) {
      std::string attribute(nd->first_attribute("amplitude")->value());
      if (!stringutils::extractFromString(attribute, scr.amplitude)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of amplitude attribute for script. Value must be "
                     "scalar. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    scr.base_dt = 1.0;
    if (nd->first_attribute("dt")) {
      std::string attribute(nd->first_attribute("dt")->value());
      if (!stringutils::extractFromString(attribute, scr.base_dt)) {
        std::cerr
            << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse value "
               "of dt attribute for script. Value must be scalar. Exiting."
            << std::endl;
        exit(1);
      }
    }

    if (nd->first_attribute("base")) {
      std::string attribute(nd->first_attribute("base")->value());
      std::vector<std::string> bases = stringutils::split(attribute, ' ');

      scr.base_vertices.reserve(bases.size());
      for (const std::string& str : bases) {
        scalar y = 0;
        stringutils::extractFromString(str, y);
        scr.base_vertices.push_back(y);
      }
    }

    scr.frequency = 1.0;
    if (nd->first_attribute("frequency")) {
      std::string attribute(nd->first_attribute("frequency")->value());
      if (!stringutils::extractFromString(attribute, scr.frequency)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of frequency attribute for script. Value must be "
                     "scalar. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    scr.index = -1;
    if (nd->first_attribute("i")) {
      std::string attribute(nd->first_attribute("i")->value());
      if (!stringutils::extractFromString(attribute, scr.index)) {
        std::cerr
            << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse value "
               "of i attribute for script. Value must be integer. Exiting."
            << std::endl;
        exit(1);
      }
    }

    scr.updateSDF = false;
    if (nd->first_attribute("updatesdf")) {
      std::string attribute(nd->first_attribute("updatesdf")->value());
      if (!stringutils::extractFromString(attribute, scr.updateSDF)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of updatesdf attribute for script. Value must be "
                     "boolean. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    scr.transform_global = false;
    if (nd->first_attribute("global")) {
      std::string attribute(nd->first_attribute("global")->value());
      if (!stringutils::extractFromString(attribute, scr.transform_global)) {
        std::cerr
            << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse value "
               "of global attribute for script. Value must be boolean. Exiting."
            << std::endl;
        exit(1);
      }
    }

    scripts.push_back(scr);
  }
}

template <>
void TwoDSceneXMLParser::loadFluidSim(rapidxml::xml_node<>* node,
                                      TwoDScene<2>& twoscene) {
  rapidxml::xml_node<>* nd = node->first_node("fluidsim");

  if (!nd) return;

  Vector2s origin;
  origin.setZero();
  scalar grid_width;
  int grid_resolution_x = 1;
  int grid_resolution_y = 1;

  std::vector<FluidSim::Boundary<2>*> boundaries;
  std::vector<FluidSim::SourceBoundary<2>*> sources;

  if (nd->first_attribute("ox")) {
    std::string attribute(nd->first_attribute("ox")->value());
    if (!stringutils::extractFromString(attribute, origin(0))) {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                   "value of ox attribute for fluidsim parameters. Value must "
                   "be numeric. Exiting."
                << std::endl;
      exit(1);
    }
  }

  if (nd->first_attribute("oy")) {
    std::string attribute(nd->first_attribute("oy")->value());
    if (!stringutils::extractFromString(attribute, origin(1))) {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                   "value of oy attribute for fluidsim parameters. Value must "
                   "be numeric. Exiting."
                << std::endl;
      exit(1);
    }
  }

  if (nd->first_attribute("width")) {
    std::string attribute(nd->first_attribute("width")->value());
    if (!stringutils::extractFromString(attribute, grid_width)) {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                   "value of width attribute for fluidsim parameters. Value "
                   "must be numeric. Exiting."
                << std::endl;
      exit(1);
    }
  }

  if (nd->first_attribute("gx")) {
    std::string attribute(nd->first_attribute("gx")->value());
    if (!stringutils::extractFromString(attribute, grid_resolution_x)) {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                   "value of gx attribute for fluidsim parameters. Value must "
                   "be numeric. Exiting."
                << std::endl;
      exit(1);
    }
  }

  if (nd->first_attribute("gy")) {
    std::string attribute(nd->first_attribute("gy")->value());
    if (!stringutils::extractFromString(attribute, grid_resolution_y)) {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                   "value of gy attribute for fluidsim parameters. Value must "
                   "be numeric. Exiting."
                << std::endl;
      exit(1);
    }
  }

  scalar grid_dz = grid_width / (scalar)grid_resolution_x;
  if (nd->first_attribute("dz")) {
    std::string attribute(nd->first_attribute("dz")->value());
    if (!stringutils::extractFromString(attribute, grid_dz)) {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                   "value of dz attribute for fluidsim parameters. Value must "
                   "be numeric. Exiting."
                << std::endl;
      exit(1);
    }
  }

  for (rapidxml::xml_node<>* subnd = nd->first_node("boundary"); subnd;
       subnd = subnd->next_sibling("boundary")) {
    FluidSim::BOUNDARY_TYPE bt = FluidSim::BT_COUNT;
    if (subnd->first_attribute("type")) {
      std::string handlertype(subnd->first_attribute("type")->value());
      if (handlertype == "circle")
        bt = FluidSim::BT_CIRCLE;
      else if (handlertype == "box")
        bt = FluidSim::BT_BOX;
      else if (handlertype == "capsule")
        bt = FluidSim::BT_CAPSULE;
      else if (handlertype == "union")
        bt = FluidSim::BT_UNION;
      else if (handlertype == "intersect")
        bt = FluidSim::BT_INTERSECT;
      else {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of type attribute for boundary parameters. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    if (bt == FluidSim::BT_BOX || bt == FluidSim::BT_CIRCLE ||
        bt == FluidSim::BT_CAPSULE) {
      Vector2s eject_vel = Vector2s::Zero();
      Vector2s center;
      scalar rangle = 0.0;
      bool inside = false;
      bool is_source = false;
      scalar spray = 0.0;
      scalar dropradiusprop = 1.0;

      if (subnd->first_attribute("cx")) {
        std::string attribute(subnd->first_attribute("cx")->value());
        if (!stringutils::extractFromString(attribute, center(0))) {
          std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to "
                       "parse value of center(0) attribute for boundary "
                       "parameters. Value must be numeric. Exiting."
                    << std::endl;
          exit(1);
        }
      }

      if (subnd->first_attribute("cy")) {
        std::string attribute(subnd->first_attribute("cy")->value());
        if (!stringutils::extractFromString(attribute, center(1))) {
          std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to "
                       "parse value of center(1) attribute for boundary "
                       "parameters. Value must be numeric. Exiting."
                    << std::endl;
          exit(1);
        }
      }

      if (subnd->first_attribute("rw")) {
        std::string attribute(subnd->first_attribute("rw")->value());
        if (!stringutils::extractFromString(attribute, rangle)) {
          std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to "
                       "parse value of rw attribute for boundary parameters. "
                       "Value must be numeric. Exiting."
                    << std::endl;
          exit(1);
        }
      }

      if (subnd->first_attribute("inside")) {
        std::string attribute(subnd->first_attribute("inside")->value());
        if (!stringutils::extractFromString(attribute, inside)) {
          std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to "
                       "parse value of inside attribute for fluidsim "
                       "parameters. Value must be boolean. Exiting."
                    << std::endl;
          exit(1);
        }
      }

      if (subnd->first_attribute("evx")) {
        std::string attribute(subnd->first_attribute("evx")->value());
        if (!stringutils::extractFromString(attribute, eject_vel(0))) {
          std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to "
                       "parse value of evx attribute for fluidsim parameters. "
                       "Value must be boolean. Exiting."
                    << std::endl;
          exit(1);
        }

        is_source = true;
      }

      if (subnd->first_attribute("evy")) {
        std::string attribute(subnd->first_attribute("evy")->value());
        if (!stringutils::extractFromString(attribute, eject_vel(1))) {
          std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to "
                       "parse value of evy attribute for fluidsim parameters. "
                       "Value must be boolean. Exiting."
                    << std::endl;
          exit(1);
        }

        is_source = true;
      }

      scalar start_time = 0.0;
      scalar end_time = 1e+20;

      if (subnd->first_attribute("start")) {
        std::string attribute(subnd->first_attribute("start")->value());
        if (!stringutils::extractFromString(attribute, start_time)) {
          std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to "
                       "parse value of start attribute for fluidsim "
                       "parameters. Value must be boolean. Exiting."
                    << std::endl;
          exit(1);
        }

        is_source = true;
      }

      if (subnd->first_attribute("end")) {
        std::string attribute(subnd->first_attribute("end")->value());
        if (!stringutils::extractFromString(attribute, end_time)) {
          std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to "
                       "parse value of end attribute for fluidsim parameters. "
                       "Value must be boolean. Exiting."
                    << std::endl;
          exit(1);
        }

        is_source = true;
      }

      if (subnd->first_attribute("spray")) {
        std::string attribute(subnd->first_attribute("spray")->value());
        if (!stringutils::extractFromString(attribute, spray)) {
          std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to "
                       "parse value of spray attribute for fluidsim "
                       "parameters. Value must be numeric. Exiting."
                    << std::endl;
          exit(1);
        }

        is_source = true;
      }

      if (subnd->first_attribute("dropradiusprop")) {
        std::string attribute(
            subnd->first_attribute("dropradiusprop")->value());
        if (!stringutils::extractFromString(attribute, dropradiusprop)) {
          std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to "
                       "parse value of dropradius attribute for fluidsim "
                       "parameters. Value must be numeric. Exiting."
                    << std::endl;
          exit(1);
        }

        is_source = true;
      }

      scalar activate = 0.0;
      scalar inactivate = 0.0;

      if (subnd->first_attribute("activate")) {
        std::string attribute(subnd->first_attribute("activate")->value());
        if (!stringutils::extractFromString(attribute, activate)) {
          std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to "
                       "parse value of activate attribute for fluidsim "
                       "parameters. Value must be numeric. Exiting."
                    << std::endl;
          exit(1);
        }

        is_source = true;
      }

      if (subnd->first_attribute("inactivate")) {
        std::string attribute(subnd->first_attribute("inactivate")->value());
        if (!stringutils::extractFromString(attribute, inactivate)) {
          std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to "
                       "parse value of inactivate attribute for fluidsim "
                       "parameters. Value must be numeric. Exiting."
                    << std::endl;
          exit(1);
        }

        is_source = true;
      }

      Vector4s parameter(0, 0, 0, 0);

      switch (bt) {
        case FluidSim::BT_CIRCLE:
          if (subnd->first_attribute("radius")) {
            std::string attribute(subnd->first_attribute("radius")->value());
            if (!stringutils::extractFromString(attribute, parameter(0))) {
              std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to "
                           "parse value of radius attribute for fluidsim "
                           "parameters. Value must be numeric. Exiting."
                        << std::endl;
              exit(1);
            }
          }
          break;
        case FluidSim::BT_BOX:
          if (subnd->first_attribute("ex")) {
            std::string attribute(subnd->first_attribute("ex")->value());
            if (!stringutils::extractFromString(attribute, parameter(0))) {
              std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to "
                           "parse value of ex attribute for fluidsim "
                           "parameters. Value must be numeric. Exiting."
                        << std::endl;
              exit(1);
            }
          }

          if (subnd->first_attribute("ey")) {
            std::string attribute(subnd->first_attribute("ey")->value());
            if (!stringutils::extractFromString(attribute, parameter(1))) {
              std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to "
                           "parse value of ey attribute for fluidsim "
                           "parameters. Value must be numeric. Exiting."
                        << std::endl;
              exit(1);
            }
          }

          if (subnd->first_attribute("radius")) {
            std::string attribute(subnd->first_attribute("radius")->value());
            if (!stringutils::extractFromString(attribute, parameter(3))) {
              std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to "
                           "parse value of radius attribute for fluidsim "
                           "parameters. Value must be numeric. Exiting."
                        << std::endl;
              exit(1);
            }
          }
          break;
        case FluidSim::BT_CAPSULE:
          if (subnd->first_attribute("radius")) {
            std::string attribute(subnd->first_attribute("radius")->value());
            if (!stringutils::extractFromString(attribute, parameter(0))) {
              std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to "
                           "parse value of radius attribute for fluidsim "
                           "parameters. Value must be numeric. Exiting."
                        << std::endl;
              exit(1);
            }
          }
          if (subnd->first_attribute("halflength")) {
            std::string attribute(
                subnd->first_attribute("halflength")->value());
            if (!stringutils::extractFromString(attribute, parameter(1))) {
              std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to "
                           "parse value of halflength attribute for fluidsim "
                           "parameters. Value must be numeric. Exiting."
                        << std::endl;
              exit(1);
            }
          }
          break;
        default:
          break;
      }

      if (is_source) {
        sources.push_back(new FluidSim::SourceBoundary<2>(
            center, parameter, bt, inside, Vector3s(0, 0, 1.0), rangle,
            eject_vel, start_time, end_time, spray, dropradiusprop, activate,
            inactivate));
      } else {
        boundaries.push_back(new FluidSim::SolidBoundary<2>(
            center, parameter, bt, inside, Vector3s(0, 0, 1.0), rangle));
      }
    } else if (bt == FluidSim::BT_UNION || bt == FluidSim::BT_INTERSECT) {
      FluidSim::OperatorBoundary<2>* ob = new FluidSim::OperatorBoundary<2>(bt);

      if (subnd->first_attribute("i")) {
        std::string attributes(subnd->first_attribute("i")->value());
        std::vector<std::string> split_attr =
            stringutils::split(attributes, ' ');
        for (const std::string& si : split_attr) {
          int i;
          if (!stringutils::extractFromString(si, i)) {
            std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to "
                         "parse value of i attribute for boundary parameters. "
                         "Value must be integer sequences. Exiting."
                      << std::endl;
            delete ob;
            exit(1);
          } else if (i >= boundaries.size()) {
            std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to "
                         "parse value of i attribute for boundary parameters. "
                         "Value must refers to existing Boundaries. Exiting."
                      << std::endl;
            delete ob;
            exit(1);
          }
          ob->children.push_back(boundaries[i]);
          boundaries[i]->parent = (int)boundaries.size();
        }
      } else {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of i attribute for boundary parameters. Value must "
                     "be integer sequences. Exiting."
                  << std::endl;
        delete ob;
        exit(1);
      }

      boundaries.push_back(ob);
    }
  }

  FluidSim2D* fluidsim =
      new FluidSim2D(origin, grid_width, grid_resolution_x, grid_resolution_y,
                     grid_dz, boundaries, sources, &twoscene);

  if (nd->first_attribute("init")) {
    std::string attribute(nd->first_attribute("init")->value());
    if (attribute == "random") {
      scalar rl = 0.0;
      scalar rr = 1.0;
      scalar rb = 0.0;
      scalar rt = 1.0;
      if (nd->first_attribute("rl")) {
        std::string sub_attr(nd->first_attribute("rl")->value());
        stringutils::extractFromString(sub_attr, rl);
      }
      if (nd->first_attribute("rr")) {
        std::string sub_attr(nd->first_attribute("rr")->value());
        stringutils::extractFromString(sub_attr, rr);
      }
      if (nd->first_attribute("rb")) {
        std::string sub_attr(nd->first_attribute("rb")->value());
        stringutils::extractFromString(sub_attr, rb);
      }
      if (nd->first_attribute("rt")) {
        std::string sub_attr(nd->first_attribute("rt")->value());
        stringutils::extractFromString(sub_attr, rt);
      }

      fluidsim->init_random_particles(rl, rr, rb, rt);
    } else if (attribute == "none") {
      // do nothing
    }
  }

  fluidsim->init_hair_particles();
  fluidsim->sort_particles();
  fluidsim->compute_liquid_phi();

  twoscene.setFluidSim(fluidsim);

  int nfl = twoscene.getNumFlows();
  for (int i = 0; i < nfl; ++i)
    twoscene.insertForce(new LevelSetForce<2>(&twoscene, fluidsim, i));
}

template <>
void TwoDSceneXMLParser::loadFluidSim(rapidxml::xml_node<>* node,
                                      TwoDScene<3>& twoscene) {
  rapidxml::xml_node<>* nd = node->first_node("fluidsim");

  if (!nd) return;

  Vector3s origin;
  origin.setZero();
  scalar grid_width;
  int grid_resolution_x = 1;
  int grid_resolution_y = 1;
  int grid_resolution_z = 1;

  std::vector<FluidSim::Boundary<3>*> boundaries;
  std::vector<FluidSim::SourceBoundary<3>*> sources;

  if (nd->first_attribute("ox")) {
    std::string attribute(nd->first_attribute("ox")->value());
    if (!stringutils::extractFromString(attribute, origin(0))) {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                   "value of ox attribute for fluidsim parameters. Value must "
                   "be numeric. Exiting."
                << std::endl;
      exit(1);
    }
  }

  if (nd->first_attribute("oy")) {
    std::string attribute(nd->first_attribute("oy")->value());
    if (!stringutils::extractFromString(attribute, origin(1))) {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                   "value of oy attribute for fluidsim parameters. Value must "
                   "be numeric. Exiting."
                << std::endl;
      exit(1);
    }
  }

  if (nd->first_attribute("oz")) {
    std::string attribute(nd->first_attribute("oz")->value());
    if (!stringutils::extractFromString(attribute, origin(2))) {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                   "value of oy attribute for fluidsim parameters. Value must "
                   "be numeric. Exiting."
                << std::endl;
      exit(1);
    }
  }

  if (nd->first_attribute("width")) {
    std::string attribute(nd->first_attribute("width")->value());
    if (!stringutils::extractFromString(attribute, grid_width)) {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                   "value of width attribute for fluidsim parameters. Value "
                   "must be numeric. Exiting."
                << std::endl;
      exit(1);
    }
  }

  if (nd->first_attribute("gx")) {
    std::string attribute(nd->first_attribute("gx")->value());
    if (!stringutils::extractFromString(attribute, grid_resolution_x)) {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                   "value of gx attribute for fluidsim parameters. Value must "
                   "be numeric. Exiting."
                << std::endl;
      exit(1);
    }
  }

  if (nd->first_attribute("gy")) {
    std::string attribute(nd->first_attribute("gy")->value());
    if (!stringutils::extractFromString(attribute, grid_resolution_y)) {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                   "value of gy attribute for fluidsim parameters. Value must "
                   "be numeric. Exiting."
                << std::endl;
      exit(1);
    }
  }

  if (nd->first_attribute("gz")) {
    std::string attribute(nd->first_attribute("gz")->value());
    if (!stringutils::extractFromString(attribute, grid_resolution_z)) {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                   "value of gz attribute for fluidsim parameters. Value must "
                   "be numeric. Exiting."
                << std::endl;
      exit(1);
    }
  }

  for (rapidxml::xml_node<>* subnd = nd->first_node("boundary"); subnd;
       subnd = subnd->next_sibling("boundary")) {
    FluidSim::BOUNDARY_TYPE bt = FluidSim::BT_COUNT;
    if (subnd->first_attribute("type")) {
      std::string handlertype(subnd->first_attribute("type")->value());
      if (handlertype == "circle")
        bt = FluidSim::BT_CIRCLE;
      else if (handlertype == "box")
        bt = FluidSim::BT_BOX;
      else if (handlertype == "capsule")
        bt = FluidSim::BT_CAPSULE;
      else if (handlertype == "union")
        bt = FluidSim::BT_UNION;
      else if (handlertype == "intersect")
        bt = FluidSim::BT_INTERSECT;
      else {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of type attribute for boundary parameters. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    if (bt == FluidSim::BT_BOX || bt == FluidSim::BT_CIRCLE ||
        bt == FluidSim::BT_CAPSULE) {
      Vector3s eject_vel = Vector3s::Zero();
      Vector3s center;
      Vector3s raxis = Vector3s(0, 1, 0);
      scalar spray = 0.0;
      scalar dropradiusprop = 1.0;
      scalar rangle = 0.0;
      bool is_source = false;

      if (subnd->first_attribute("cx")) {
        std::string attribute(subnd->first_attribute("cx")->value());
        if (!stringutils::extractFromString(attribute, center(0))) {
          std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to "
                       "parse value of center(0) attribute for fluidsim "
                       "parameters. Value must be numeric. Exiting."
                    << std::endl;
          exit(1);
        }
      }

      if (subnd->first_attribute("cy")) {
        std::string attribute(subnd->first_attribute("cy")->value());
        if (!stringutils::extractFromString(attribute, center(1))) {
          std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to "
                       "parse value of center(1) attribute for fluidsim "
                       "parameters. Value must be numeric. Exiting."
                    << std::endl;
          exit(1);
        }
      }

      if (subnd->first_attribute("cz")) {
        std::string attribute(subnd->first_attribute("cz")->value());
        if (!stringutils::extractFromString(attribute, center(2))) {
          std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to "
                       "parse value of center(2) attribute for fluidsim "
                       "parameters. Value must be numeric. Exiting."
                    << std::endl;
          exit(1);
        }
      }

      if (subnd->first_attribute("rx")) {
        std::string attribute(subnd->first_attribute("rx")->value());
        if (!stringutils::extractFromString(attribute, raxis(0))) {
          std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to "
                       "parse value of center(0) attribute for fluidsim "
                       "parameters. Value must be numeric. Exiting."
                    << std::endl;
          exit(1);
        }
      }

      if (subnd->first_attribute("ry")) {
        std::string attribute(subnd->first_attribute("ry")->value());
        if (!stringutils::extractFromString(attribute, raxis(1))) {
          std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to "
                       "parse value of center(1) attribute for fluidsim "
                       "parameters. Value must be numeric. Exiting."
                    << std::endl;
          exit(1);
        }
      }

      if (subnd->first_attribute("rz")) {
        std::string attribute(subnd->first_attribute("rz")->value());
        if (!stringutils::extractFromString(attribute, raxis(2))) {
          std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to "
                       "parse value of center(2) attribute for fluidsim "
                       "parameters. Value must be numeric. Exiting."
                    << std::endl;
          exit(1);
        }
      }

      if (subnd->first_attribute("rw")) {
        std::string attribute(subnd->first_attribute("rw")->value());
        if (!stringutils::extractFromString(attribute, rangle)) {
          std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to "
                       "parse value of center(2) attribute for fluidsim "
                       "parameters. Value must be numeric. Exiting."
                    << std::endl;
          exit(1);
        }
      }

      bool inside = false;
      if (subnd->first_attribute("inside")) {
        std::string attribute(subnd->first_attribute("inside")->value());
        if (!stringutils::extractFromString(attribute, inside)) {
          std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to "
                       "parse value of inside attribute for fluidsim "
                       "parameters. Value must be boolean. Exiting."
                    << std::endl;
          exit(1);
        }
      }

      if (subnd->first_attribute("evx")) {
        std::string attribute(subnd->first_attribute("evx")->value());
        if (!stringutils::extractFromString(attribute, eject_vel(0))) {
          std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to "
                       "parse value of evx attribute for fluidsim parameters. "
                       "Value must be boolean. Exiting."
                    << std::endl;
          exit(1);
        }

        is_source = true;
      }

      if (subnd->first_attribute("evy")) {
        std::string attribute(subnd->first_attribute("evy")->value());
        if (!stringutils::extractFromString(attribute, eject_vel(1))) {
          std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to "
                       "parse value of evy attribute for fluidsim parameters. "
                       "Value must be boolean. Exiting."
                    << std::endl;
          exit(1);
        }

        is_source = true;
      }

      if (subnd->first_attribute("evz")) {
        std::string attribute(subnd->first_attribute("evz")->value());
        if (!stringutils::extractFromString(attribute, eject_vel(2))) {
          std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to "
                       "parse value of evz attribute for fluidsim parameters. "
                       "Value must be boolean. Exiting."
                    << std::endl;
          exit(1);
        }

        is_source = true;
      }

      scalar start_time = 0.0;
      scalar end_time = 1e+20;

      if (subnd->first_attribute("start")) {
        std::string attribute(subnd->first_attribute("start")->value());
        if (!stringutils::extractFromString(attribute, start_time)) {
          std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to "
                       "parse value of start attribute for fluidsim "
                       "parameters. Value must be boolean. Exiting."
                    << std::endl;
          exit(1);
        }

        is_source = true;
      }

      if (subnd->first_attribute("end")) {
        std::string attribute(subnd->first_attribute("end")->value());
        if (!stringutils::extractFromString(attribute, end_time)) {
          std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to "
                       "parse value of end attribute for fluidsim parameters. "
                       "Value must be boolean. Exiting."
                    << std::endl;
          exit(1);
        }

        is_source = true;
      }

      if (subnd->first_attribute("spray")) {
        std::string attribute(subnd->first_attribute("spray")->value());
        if (!stringutils::extractFromString(attribute, spray)) {
          std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to "
                       "parse value of spray attribute for fluidsim "
                       "parameters. Value must be numeric. Exiting."
                    << std::endl;
          exit(1);
        }

        is_source = true;
      }

      if (subnd->first_attribute("dropradiusprop")) {
        std::string attribute(
            subnd->first_attribute("dropradiusprop")->value());
        if (!stringutils::extractFromString(attribute, dropradiusprop)) {
          std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to "
                       "parse value of dropradius attribute for fluidsim "
                       "parameters. Value must be numeric. Exiting."
                    << std::endl;
          exit(1);
        }

        is_source = true;
      }

      scalar activate = 0.0;
      scalar inactivate = 0.0;

      if (subnd->first_attribute("activate")) {
        std::string attribute(subnd->first_attribute("activate")->value());
        if (!stringutils::extractFromString(attribute, activate)) {
          std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to "
                       "parse value of activate attribute for fluidsim "
                       "parameters. Value must be numeric. Exiting."
                    << std::endl;
          exit(1);
        }

        is_source = true;
      }

      if (subnd->first_attribute("inactivate")) {
        std::string attribute(subnd->first_attribute("inactivate")->value());
        if (!stringutils::extractFromString(attribute, inactivate)) {
          std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to "
                       "parse value of inactivate attribute for fluidsim "
                       "parameters. Value must be numeric. Exiting."
                    << std::endl;
          exit(1);
        }

        is_source = true;
      }

      Vector4s parameter(0, 0, 0, 0);

      switch (bt) {
        case FluidSim::BT_CIRCLE:
          if (subnd->first_attribute("radius")) {
            std::string attribute(subnd->first_attribute("radius")->value());
            if (!stringutils::extractFromString(attribute, parameter(0))) {
              std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to "
                           "parse value of radius attribute for fluidsim "
                           "parameters. Value must be numeric. Exiting."
                        << std::endl;
              exit(1);
            }
          }
          break;
        case FluidSim::BT_BOX:
          if (subnd->first_attribute("ex")) {
            std::string attribute(subnd->first_attribute("ex")->value());
            if (!stringutils::extractFromString(attribute, parameter(0))) {
              std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to "
                           "parse value of ex attribute for fluidsim "
                           "parameters. Value must be numeric. Exiting."
                        << std::endl;
              exit(1);
            }
          }

          if (subnd->first_attribute("ey")) {
            std::string attribute(subnd->first_attribute("ey")->value());
            if (!stringutils::extractFromString(attribute, parameter(1))) {
              std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to "
                           "parse value of ey attribute for fluidsim "
                           "parameters. Value must be numeric. Exiting."
                        << std::endl;
              exit(1);
            }
          }

          if (subnd->first_attribute("ez")) {
            std::string attribute(subnd->first_attribute("ez")->value());
            if (!stringutils::extractFromString(attribute, parameter(2))) {
              std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to "
                           "parse value of ez attribute for fluidsim "
                           "parameters. Value must be numeric. Exiting."
                        << std::endl;
              exit(1);
            }
          }

          if (subnd->first_attribute("radius")) {
            std::string attribute(subnd->first_attribute("radius")->value());
            if (!stringutils::extractFromString(attribute, parameter(3))) {
              std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to "
                           "parse value of radius attribute for fluidsim "
                           "parameters. Value must be numeric. Exiting."
                        << std::endl;
              exit(1);
            }
          }
          break;
        case FluidSim::BT_CAPSULE:
          if (subnd->first_attribute("radius")) {
            std::string attribute(subnd->first_attribute("radius")->value());
            if (!stringutils::extractFromString(attribute, parameter(0))) {
              std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to "
                           "parse value of radius attribute for fluidsim "
                           "parameters. Value must be numeric. Exiting."
                        << std::endl;
              exit(1);
            }
          }
          if (subnd->first_attribute("halflength")) {
            std::string attribute(
                subnd->first_attribute("halflength")->value());
            if (!stringutils::extractFromString(attribute, parameter(1))) {
              std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to "
                           "parse value of halflength attribute for fluidsim "
                           "parameters. Value must be numeric. Exiting."
                        << std::endl;
              exit(1);
            }
          }
          break;
        default:
          break;
      }

      if (is_source) {
        sources.push_back(new FluidSim::SourceBoundary<3>(
            center, parameter, bt, inside, raxis, rangle, eject_vel, start_time,
            end_time, spray, dropradiusprop, activate, inactivate));
      } else {
        boundaries.push_back(new FluidSim::SolidBoundary<3>(
            center, parameter, bt, inside, raxis, rangle));
      }
    } else if (bt == FluidSim::BT_UNION || bt == FluidSim::BT_INTERSECT) {
      FluidSim::OperatorBoundary<3>* ob = new FluidSim::OperatorBoundary<3>(bt);

      if (subnd->first_attribute("i")) {
        std::string attributes(subnd->first_attribute("i")->value());
        std::vector<std::string> split_attr =
            stringutils::split(attributes, ' ');
        for (const std::string& si : split_attr) {
          int i;
          if (!stringutils::extractFromString(si, i)) {
            std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to "
                         "parse value of i attribute for boundary parameters. "
                         "Value must be integer sequences. Exiting."
                      << std::endl;
            delete ob;
            exit(1);
          } else if (i >= boundaries.size()) {
            std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to "
                         "parse value of i attribute for boundary parameters. "
                         "Value must refers to existing Boundaries. Exiting."
                      << std::endl;
            delete ob;
            exit(1);
          }
          ob->children.push_back(boundaries[i]);
          boundaries[i]->parent = (int)boundaries.size();
        }
      } else {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of i attribute for boundary parameters. Value must "
                     "be integer sequences. Exiting."
                  << std::endl;
        delete ob;
        exit(1);
      }

      boundaries.push_back(ob);
    }
  }
  FluidSim3D* fluidsim =
      new FluidSim3D(origin, grid_width, grid_resolution_x, grid_resolution_y,
                     grid_resolution_z, boundaries, sources, &twoscene);

  if (nd->first_attribute("init")) {
    std::string attribute(nd->first_attribute("init")->value());
    if (attribute == "random") {
      scalar rl = 0.0;
      scalar rr = 1.0;
      scalar rb = 0.0;
      scalar rt = 1.0;
      scalar rback = 0.0;
      scalar rforw = 1.0;
      if (nd->first_attribute("rl")) {
        std::string sub_attr(nd->first_attribute("rl")->value());
        stringutils::extractFromString(sub_attr, rl);
      }
      if (nd->first_attribute("rr")) {
        std::string sub_attr(nd->first_attribute("rr")->value());
        stringutils::extractFromString(sub_attr, rr);
      }
      if (nd->first_attribute("rb")) {
        std::string sub_attr(nd->first_attribute("rb")->value());
        stringutils::extractFromString(sub_attr, rb);
      }
      if (nd->first_attribute("rt")) {
        std::string sub_attr(nd->first_attribute("rt")->value());
        stringutils::extractFromString(sub_attr, rt);
      }
      if (nd->first_attribute("rback")) {
        std::string sub_attr(nd->first_attribute("rback")->value());
        stringutils::extractFromString(sub_attr, rback);
      }
      if (nd->first_attribute("rforw")) {
        std::string sub_attr(nd->first_attribute("rforw")->value());
        stringutils::extractFromString(sub_attr, rforw);
      }

      fluidsim->init_random_particles(rl, rr, rb, rt, rback, rforw);
    } else if (attribute == "none") {
      // do nothing
    }
  }

  twoscene.setFluidSim(fluidsim);

  int nfl = twoscene.getNumFlows();
  for (int i = 0; i < nfl; ++i)
    twoscene.insertForce(new LevelSetForce<3>(&twoscene, fluidsim, i));
}

template <int DIM>
void TwoDSceneXMLParser::loadLiquidParameter(rapidxml::xml_node<>* node,
                                             TwoDScene<DIM>& twoscene) {
  rapidxml::xml_node<>* nd = node->first_node("liquid");

  WetHairParameter parameter;

  if (nd) {
    rapidxml::xml_attribute<>* dtnd = nd->first_attribute("dt");
    if (dtnd) {
      if (!stringutils::extractFromString(std::string(dtnd->value()),
                                          parameter.dt)) {
        std::cerr
            << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse 'dt' "
               "attribute for integrator. Value must be numeric. Exiting."
            << std::endl;
        exit(1);
      }
    }

    rapidxml::xml_attribute<>* hnsnd = nd->first_attribute("hairstep");
    if (hnsnd) {
      if (!stringutils::extractFromString(std::string(hnsnd->value()),
                                          parameter.hairsteps)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "'hairstep' attribute for integrator. Value must be "
                     "integer. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    rapidxml::xml_attribute<>* swesnd = nd->first_attribute("swestep");
    if (swesnd) {
      if (!stringutils::extractFromString(std::string(swesnd->value()),
                                          parameter.swesteps)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "'swestep' attribute for integrator. Value must be "
                     "integer. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    if (nd->first_attribute("rho")) {
      std::string attribute(nd->first_attribute("rho")->value());
      if (!stringutils::extractFromString(attribute, parameter.rho)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of rho attribute for liquid parameters. Value must "
                     "be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    if (nd->first_attribute("sigma")) {
      std::string attribute(nd->first_attribute("sigma")->value());
      if (!stringutils::extractFromString(attribute, parameter.sigma)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of sigma attribute for liquid parameters. Value "
                     "must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    // Attempt to load the duration value
    rapidxml::xml_attribute<>* timend = nd->first_attribute("theta");
    if (timend != NULL) {
      if (!stringutils::extractFromString(std::string(timend->value()),
                                          parameter.theta)) {
        std::cerr
            << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
               "'theta' attribute for liquid. Value must be numeric. Exiting."
            << std::endl;
        exit(1);
      }
    }

    timend = nd->first_attribute("viscosity");
    if (timend != NULL) {
      if (!stringutils::extractFromString(std::string(timend->value()),
                                          parameter.viscosity)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "'viscosity' attribute for liquid. Value must be numeric. "
                     "Exiting."
                  << std::endl;
        exit(1);
      }
    }

    timend = nd->first_attribute("quadraticdragging");
    if (timend != NULL) {
      if (!stringutils::extractFromString(std::string(timend->value()),
                                          parameter.quadratic_dragging)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "'quadraticdragging' attribute for liquid. Value must be "
                     "numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    timend = nd->first_attribute("friction");
    if (timend != NULL) {
      if (!stringutils::extractFromString(std::string(timend->value()),
                                          parameter.friction)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "'friction' attribute for liquid. Value must be numeric. "
                     "Exiting."
                  << std::endl;
        exit(1);
      }
    }

    timend = nd->first_attribute("radiusmultiplier");
    if (timend != NULL) {
      if (!stringutils::extractFromString(std::string(timend->value()),
                                          parameter.radius_multiplier)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "'radiusmultiplier' attribute for liquid. Value must be "
                     "numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    timend = nd->first_attribute("collisionstiffness");
    if (timend != NULL) {
      if (!stringutils::extractFromString(std::string(timend->value()),
                                          parameter.collision_stiffness)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "'collisionstiffness' attribute for liquid. Value must be "
                     "numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    timend = nd->first_attribute("radiusmultiplierplanar");
    if (timend != NULL) {
      if (!stringutils::extractFromString(std::string(timend->value()),
                                          parameter.radius_multiplier_planar)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "'radiusmultiplierplanar' attribute for liquid. Value "
                     "must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    timend = nd->first_attribute("dragradiusmultiplier");
    if (timend != NULL) {
      if (!stringutils::extractFromString(std::string(timend->value()),
                                          parameter.drag_radius_multiplier)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "'dragradiusmultiplier' attribute for liquid. Value must "
                     "be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    timend = nd->first_attribute("collisionstiffnessplanar");
    if (timend != NULL) {
      if (!stringutils::extractFromString(
              std::string(timend->value()),
              parameter.collision_stiffness_planar)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "'collisionstiffnessplanar' attribute for liquid. Value "
                     "must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    timend = nd->first_attribute("dampingmultiplier");
    if (timend != NULL) {
      if (!stringutils::extractFromString(std::string(timend->value()),
                                          parameter.damping_multiplier)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "'dampingmultiplier' attribute for liquid. Value must be "
                     "numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    timend = nd->first_attribute("dampingmultiplierplanar");
    if (timend != NULL) {
      if (!stringutils::extractFromString(
              std::string(timend->value()),
              parameter.damping_multiplier_planar)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "'dampingmultiplierplanar' attribute for liquid. Value "
                     "must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    timend = nd->first_attribute("bulkthresholdmultiplier");
    if (timend != NULL) {
      if (!stringutils::extractFromString(
              std::string(timend->value()),
              parameter.bulk_threshold_multiplier)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "'bulkthresholdmultiplier' attribute for liquid. Value "
                     "must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    timend = nd->first_attribute("hairhaircohesion");
    if (timend != NULL) {
      if (!stringutils::extractFromString(
              std::string(timend->value()),
              parameter.hair_hair_cohesion_multiplier)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "'hairhaircohesion' attribute for liquid. Value must be "
                     "numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    timend = nd->first_attribute("hairsolidcohesion");
    if (timend != NULL) {
      if (!stringutils::extractFromString(
              std::string(timend->value()),
              parameter.hair_solid_cohesion_multiplier)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "'hairsolidcohesion' attribute for liquid. Value must be "
                     "numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    timend = nd->first_attribute("maxetaprop");
    if (timend != NULL) {
      if (!stringutils::extractFromString(std::string(timend->value()),
                                          parameter.max_limited_eta_prop)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "'maxetaprop' attribute for liquid. Value must be "
                     "numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    timend = nd->first_attribute("regularizershell");
    if (timend != NULL) {
      if (!stringutils::extractFromString(std::string(timend->value()),
                                          parameter.regularizer_shell)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "'regularizershell' attribute for liquid. Value must be "
                     "numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    timend = nd->first_attribute("latitude");
    if (timend != NULL) {
      if (!stringutils::extractFromString(std::string(timend->value()),
                                          parameter.latitude)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "'latitude' attribute for liquid. Value must be numeric. "
                     "Exiting."
                  << std::endl;
        exit(1);
      }
    }

    timend = nd->first_attribute("earthradius");
    if (timend != NULL) {
      if (!stringutils::extractFromString(std::string(timend->value()),
                                          parameter.earth_radius)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "'earthradius' attribute for liquid. Value must be "
                     "numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    timend = nd->first_attribute("earthrotation");
    if (timend != NULL) {
      if (!stringutils::extractFromString(std::string(timend->value()),
                                          parameter.earth_rotation)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "'earthrotation' attribute for liquid. Value must be "
                     "numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    timend = nd->first_attribute("heightsmooth");
    if (timend != NULL) {
      if (!stringutils::extractFromString(std::string(timend->value()),
                                          parameter.height_smooth)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "'heightsmooth' attribute for liquid. Value must be "
                     "numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    timend = nd->first_attribute("capillaryaccel");
    if (timend != NULL) {
      if (!stringutils::extractFromString(
              std::string(timend->value()),
              parameter.capillary_accel_multiplier)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "'capillaryaccel' attribute for liquid. Value must be "
                     "numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    timend = nd->first_attribute("drippingradius");
    if (timend != NULL) {
      if (!stringutils::extractFromString(
              std::string(timend->value()),
              parameter.dripping_radius_multiplier)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "'drippingradius' attribute for liquid. Value must be "
                     "numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    timend = nd->first_attribute("airviscosity");
    if (timend != NULL) {
      if (!stringutils::extractFromString(std::string(timend->value()),
                                          parameter.airviscosity)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "'airviscosity' attribute for liquid. Value must be "
                     "numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    timend = nd->first_attribute("massupdate");
    if (timend != NULL) {
      std::string mum = std::string(timend->value());
      if (mum == "none")
        parameter.mass_update_mode = MUM_NONE;
      else if (mum == "only")
        parameter.mass_update_mode = MUM_MASS_ONLY;
      else if (mum == "direct")
        parameter.mass_update_mode = MUM_DIRECT_DIV;
      else if (mum == "momentum")
        parameter.mass_update_mode = MUM_MOMENTUM;
    }

    timend = nd->first_attribute("globalvolumecontrol");
    if (timend != NULL) {
      if (!stringutils::extractFromString(std::string(timend->value()),
                                          parameter.global_volume_control)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "'globalvolumecontrol' attribute for liquid. Value must "
                     "be boolean. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    timend = nd->first_attribute("viscoussolve");
    if (timend != NULL) {
      if (!stringutils::extractFromString(std::string(timend->value()),
                                          parameter.viscous_solve)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "'viscoussolve' attribute for liquid. Value must be "
                     "boolean. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    if (nd->first_attribute("absorption")) {
      std::string attribute(nd->first_attribute("absorption")->value());
      if (!stringutils::extractFromString(attribute,
                                          parameter.absorptionRate)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of absorption attribute for fluidsim parameters. "
                     "Value must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    if (node->first_attribute("ctcd")) {
      std::string attribute(node->first_attribute("ctcd")->value());
      if (!stringutils::extractFromString(attribute, parameter.use_ctcd)) {
        std::cerr
            << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse value "
               "of ctcd attribute for scene. Value must be boolean. Exiting."
            << std::endl;
        exit(1);
      }
    }

    if (node->first_attribute("individualtransfer")) {
      std::string attribute(
          node->first_attribute("individualtransfer")->value());
      if (!stringutils::extractFromString(attribute,
                                          parameter.individual_transfer)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of individualtransfer attribute for scene. Value "
                     "must be boolean. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    if (node->first_attribute("nofluids")) {
      std::string attribute(node->first_attribute("nofluids")->value());
      if (!stringutils::extractFromString(attribute, parameter.no_fluids)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of nofluids attribute for scene. Value must be "
                     "boolean. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    if (node->first_attribute("noswe")) {
      std::string attribute(node->first_attribute("noswe")->value());
      if (!stringutils::extractFromString(attribute, parameter.no_swe)) {
        std::cerr
            << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse value "
               "of noswe attribute for scene. Value must be boolean. Exiting."
            << std::endl;
        exit(1);
      }
    }

    if (node->first_attribute("noabsorb")) {
      std::string attribute(node->first_attribute("noabsorb")->value());
      if (!stringutils::extractFromString(attribute, parameter.no_absorb)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of noabsorb attribute for scene. Value must be "
                     "boolean. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    if (node->first_attribute("nodripping")) {
      std::string attribute(node->first_attribute("nodripping")->value());
      if (!stringutils::extractFromString(attribute, parameter.no_dripping)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of nodripping attribute for scene. Value must be "
                     "boolean. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    if (node->first_attribute("nofictitious")) {
      std::string attribute(node->first_attribute("nofictitious")->value());
      if (!stringutils::extractFromString(attribute, parameter.no_fictitious)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of nofictitious attribute for scene. Value must be "
                     "boolean. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    if (node->first_attribute("gravityx")) {
      std::string attribute(nd->first_attribute("gravityx")->value());
      if (!stringutils::extractFromString(attribute, parameter.gravity(0))) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of gravityx attribute for constantforce. Value "
                     "must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    // Extract the y component of the force
    if (nd->first_attribute("gravityy")) {
      std::string attribute(nd->first_attribute("gravityy")->value());
      if (!stringutils::extractFromString(attribute, parameter.gravity(1))) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of fy attribute for constantforce. Value must be "
                     "numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    if (DIM == 3) {
      // Extract the y component of the force
      if (nd->first_attribute("gravityz")) {
        std::string attribute(nd->first_attribute("gravityz")->value());
        if (!stringutils::extractFromString(attribute, parameter.gravity(2))) {
          std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to "
                       "parse value of gravityz attribute for constantforce. "
                       "Value must be numeric. Exiting."
                    << std::endl;
          exit(1);
        }
      }
    }

    if (nd->first_attribute("frictionmultiplierplanar")) {
      std::string attribute(
          nd->first_attribute("frictionmultiplierplanar")->value());
      if (!stringutils::extractFromString(
              attribute, parameter.friction_multiplier_planar)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of frictionmultiplierplanar attribute for fluidsim "
                     "parameters. Value must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    if (nd->first_attribute("correction")) {
      std::string attribute(nd->first_attribute("correction")->value());
      if (!stringutils::extractFromString(attribute,
                                          parameter.fluidcorrectionsteps)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of correction attribute for fluidsim parameters. "
                     "Value must be boolean. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    if (nd->first_attribute("coriolis")) {
      std::string attribute(nd->first_attribute("coriolis")->value());
      if (!stringutils::extractFromString(attribute,
                                          parameter.apply_coriolis)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of coriolis attribute for fluidsim parameters. "
                     "Value must be boolean. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    if (nd->first_attribute("drippingmiddle")) {
      std::string attribute(nd->first_attribute("drippingmiddle")->value());
      if (!stringutils::extractFromString(attribute,
                                          parameter.drippingmiddle)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of drippingmiddle attribute for fluidsim "
                     "parameters. Value must be boolean. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    if (nd->first_attribute("drippingnear")) {
      std::string attribute(nd->first_attribute("drippingnear")->value());
      if (!stringutils::extractFromString(attribute, parameter.drippingnear)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of drippingnear attribute for fluidsim parameters. "
                     "Value must be boolean. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    if (nd->first_attribute("drippingfar")) {
      std::string attribute(nd->first_attribute("drippingfar")->value());
      if (!stringutils::extractFromString(attribute, parameter.drippingfar)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of drippingfar attribute for fluidsim parameters. "
                     "Value must be boolean. Exiting."
                  << std::endl;
        exit(1);
      }
    }
  }

  twoscene.setLiquidParameter(parameter);
  twoscene.insertForce(
      new SimpleGravityForce<DIM>(twoscene.getParameter().gravity, &twoscene));
}

template <int DIM>
void TwoDSceneXMLParser::loadLinearBendingForces(rapidxml::xml_node<>* node,
                                                 TwoDScene<DIM>& twodscene) {
  assert(node != NULL);

  int forcenum = 0;
  for (rapidxml::xml_node<>* nd = node->first_node("linearbendingforce"); nd;
       nd = nd->next_sibling("linearbendingforce")) {
    // Extract the vertices between which the force acts
    int idx1 = -1;

    if (nd->first_attribute("i1")) {
      std::string attribute(nd->first_attribute("i1")->value());
      if (!stringutils::extractFromString(attribute, idx1)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of i1 attribute for elasticbodybendingforce "
                  << forcenum << ". Value must be integer. Exiting."
                  << std::endl;
        exit(1);
      }
    } else {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                   "value of i1 attribute for elasticbodybendingforce "
                << forcenum << ". Exiting." << std::endl;
      exit(1);
    }

    int idx2 = -1;

    if (nd->first_attribute("i2")) {
      std::string attribute(nd->first_attribute("i2")->value());
      if (!stringutils::extractFromString(attribute, idx2)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of i2 attribute for elasticbodybendingforce "
                  << forcenum << ". Value must be integer. Exiting."
                  << std::endl;
        exit(1);
      }
    } else {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                   "value of i2 attribute for elasticbodybendingforce "
                << forcenum << ". Exiting." << std::endl;
      exit(1);
    }

    int idx3 = -1;

    if (nd->first_attribute("i3")) {
      std::string attribute(nd->first_attribute("i3")->value());
      if (!stringutils::extractFromString(attribute, idx3)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of i3 attribute for elasticbodybendingforce "
                  << forcenum << ". Value must be integer. Exiting."
                  << std::endl;
        exit(1);
      }
    } else {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                   "value of i3 attribute for elasticbodybendingforce "
                << forcenum << ". Exiting." << std::endl;
      exit(1);
    }

    // Extract the bending force stiffness
    scalar alpha = std::numeric_limits<scalar>::signaling_NaN();
    if (nd->first_attribute("alpha")) {
      std::string attribute(nd->first_attribute("alpha")->value());
      if (!stringutils::extractFromString(attribute, alpha)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of alpha attribute for elasticbodybendingforce "
                  << forcenum << ". Value must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    } else {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                   "alpha attribute for elasticbodybendingforce "
                << forcenum << ". Exiting." << std::endl;
      exit(1);
    }

    scalar beta = 0.0;
    if (nd->first_attribute("beta")) {
      std::string attribute(nd->first_attribute("beta")->value());
      if (!stringutils::extractFromString(attribute, beta)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of beta attribute for elasticbodybendingforce "
                  << forcenum << ". Value must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    // Extract the rest angle
    Vectors<DIM - 1> theta0 = Vectors<DIM - 1>::Zero();
    if (nd->first_attribute("theta0")) {
      std::string attribute(nd->first_attribute("theta0")->value());
      if (!stringutils::extractFromString(attribute, theta0(0))) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of theta0 attribute for elasticbodybendingforce "
                  << forcenum << ". Value must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    } else {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                   "theta0 attribute for elasticbodybendingforce "
                << forcenum << ". Exiting." << std::endl;
      exit(1);
    }

    if (DIM == 3) {
      if (nd->first_attribute("phi0")) {
        std::string attribute(nd->first_attribute("phi0")->value());
        if (!stringutils::extractFromString(attribute, theta0(1))) {
          std::cerr
              << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                 "value of phi0 attribute for elasticbodybendingforce "
              << forcenum << ". Value must be numeric. Exiting." << std::endl;
          exit(1);
        }
      }
    }

    assert(twodscene.getNumParticles() > idx1 &&
           "Please put the bending force tag below all nodes it references!");
    assert(twodscene.getNumParticles() > idx2 &&
           "Please put the bending force tag below all nodes it references!");
    assert(twodscene.getNumParticles() > idx3 &&
           "Please put the bending force tag below all nodes it references!");

    const VectorXs& x = twodscene.getX();
    Vectors<DIM> xb1 = x.segment<DIM>(twodscene.getDof(idx1));
    Vectors<DIM> xb2 = x.segment<DIM>(twodscene.getDof(idx2));
    Vectors<DIM> xb3 = x.segment<DIM>(twodscene.getDof(idx3));

    twodscene.insertForce(new LinearBendingForce<DIM>(
        &twodscene, idx1, idx2, idx3, alpha, beta, theta0, (xb2 - xb1).norm(),
        (xb3 - xb2).norm()));

    ++forcenum;
  }
}

template <int DIM>
void TwoDSceneXMLParser::loadParticleSimulation(
    GLFWwindow* window,
    bool simulate_comparison, bool rendering_enabled,
    ExecutableSimulation** execsim, renderingutils::Viewport& view,
    renderingutils::Color& bgcolor, rapidxml::xml_node<>* node,
    const std::string& init_filename) {
  TwoDScene<DIM>* scene = new TwoDScene<DIM>(true);

  // Scripts
  std::vector<Script> scripts;
  loadScripts(node, scripts);

  // Scene
  loadParticleEdges(node, *scene);

  loadLiquidParameter(node, *scene);
  loadFlow(node, *scene);
  loadFluidSim(node, *scene);

  PolygonalCohesion<DIM>* cohesion = new PolygonalCohesion<DIM>(scene);
  scene->insertForce(cohesion);
  scene->setPolygonalCohesion(cohesion);

  // Forces
  loadSpringForces(node, *scene);
  loadDragDampingForces(node, *scene);
  loadLinearSpringForces(node, *scene);
  loadLinearBendingForces(node, *scene);

  // Integrator/solver
  SceneStepper<DIM>* scene_stepper = NULL;
  loadIntegrator(node, *scene, &scene_stepper);
  scene_stepper->init(*scene);
  assert(scene_stepper != NULL);

  scene->updateHairConnectivity();
  scene->categorizeForces();
  scene->initializeScriptedGroup(scripts);

  if (!init_filename.empty()) {
    std::ifstream binary_input;
    // Attempt to open the binary
    binary_input.open(init_filename);
    if (binary_input.fail()) {
      std::cerr << outputmod::startred
                << "ERROR IN INITIALIZATION: " << outputmod::endred
                << "Failed to open binary input file: "
                << " `" << init_filename << "`   Exiting." << std::endl;
      exit(1);
    }
    assert(scene);

    TwoDSceneSerializer<DIM> serializer;
    serializer.loadScene(*scene, scene_stepper, binary_input);
  }

  scene->initializeLiquids();

  // Rendering state
  std::vector<renderingutils::Color> particle_colors;
  particle_colors.resize(
      scene->getNumParticles(),
      renderingutils::Color(0.650980392156863, 0.294117647058824, 0.0));
  loadParticleColors(node, particle_colors);

  std::vector<renderingutils::Color> edge_colors;
  edge_colors.resize(scene->getNumEdges(),
                     renderingutils::Color(0.0, 0.0, 0.0));
  loadEdgeColors(node, edge_colors);

  TwoDSceneRenderer<DIM>* scene_renderer = NULL;
  if (rendering_enabled) {
    scene_renderer =
        new TwoDSceneRenderer<DIM>(*scene, particle_colors, edge_colors);
  }

  ParticleSimulation<DIM>* ps = new ParticleSimulation<DIM>(
      window, scene, scene_stepper, scene_renderer, scripts);
  *execsim = ps;

  if (rendering_enabled) {
    scene_renderer->setDC(ps->getDC());
    ps->getDC()->setRender(scene_renderer);
  }
}

void TwoDSceneXMLParser::loadDERSimulation(GLFWwindow* window,
                                           bool simulate_comparison,
                                           bool rendering_enabled,
                                           ExecutableSimulation** execsim,
                                           renderingutils::Viewport& view,
                                           renderingutils::Color& bgcolor,
                                           rapidxml::xml_node<>* node,
                                           const std::string& init_filename) {
  const int DIM = 3;
  TwoDScene<DIM>* scene = new TwoDScene<DIM>(false);

  // Scripts
  std::vector<Script> scripts;
  loadScripts(node, scripts);

  // Scene
  // Load Strands & vertices/edges
  loadLiquidParameter(node, *scene);
  loadStrandParameters(node, *scene);
  loadStrandParticleEdges(node, *scene);

  // // Forces make it so that this code runs as a mass spring system as well.
  // flag to turn off TWIST force (bending too since uses twist dofs?), and all
  // other forces below simulatenously supported
  loadSpringForces(node, *scene);
  loadDragDampingForces(node, *scene);
  loadLinearSpringForces(node, *scene);
  loadLinearBendingForces(node, *scene);

  loadFluidSim(node, *scene);

  PolygonalCohesion<DIM>* cohesion = new PolygonalCohesion<DIM>(scene);
  scene->insertForce(cohesion);
  scene->setPolygonalCohesion(cohesion);

  // Integrator/solver
  SceneStepper<DIM>* scene_stepper = NULL;
  loadIntegrator(node, *scene, &scene_stepper);
  scene_stepper->init(*scene);
  assert(scene_stepper != NULL);

  scene->updateHairConnectivity();
  scene->categorizeForces();
  scene->initializeScriptedGroup(scripts);

  if (!init_filename.empty()) {
    std::ifstream binary_input;
    // Attempt to open the binary
    binary_input.open(init_filename);
    if (binary_input.fail()) {
      std::cerr << outputmod::startred
                << "ERROR IN INITIALIZATION: " << outputmod::endred
                << "Failed to open binary input file: "
                << " `" << init_filename << "`   Exiting." << std::endl;
      exit(1);
    }
    assert(scene);

    TwoDSceneSerializer<DIM> serializer;
    serializer.loadScene(*scene, scene_stepper, binary_input);
  }

  scene->initializeLiquids();

  // Rendering state
  std::vector<renderingutils::Color> particle_colors;
  particle_colors.resize(
      scene->getNumParticles(),
      renderingutils::Color(0.650980392156863, 0.294117647058824, 0.0));
  loadParticleColors(node, particle_colors);

  std::vector<renderingutils::Color> edge_colors;
  edge_colors.resize(
      scene->getNumEdges(),
      renderingutils::Color(scene->getStrandParameters(0)->getColor()));
  loadEdgeColors(node, edge_colors);

  TwoDSceneRenderer<DIM>* scene_renderer = NULL;
  // if( rendering_enabled )
  {
    scene_renderer =
        new TwoDSceneRenderer<DIM>(*scene, particle_colors, edge_colors);
  }

  ParticleSimulation<DIM>* ps = new ParticleSimulation<DIM>(
      window, scene, scene_stepper, scene_renderer, scripts);
  *execsim = ps;

  if (rendering_enabled) {
    scene_renderer->setDC(ps->getDC());
    ps->getDC()->setRender(scene_renderer);
  }
}

void TwoDSceneXMLParser::loadXMLFileFromMemory(const char* memory_str,
                                               std::vector<char>& xmlchars,
                                               rapidxml::xml_document<>& doc) {
  xmlchars.resize(strlen(memory_str));
  strcpy(&xmlchars[0], memory_str);
  doc.parse<0>(&xmlchars[0]);
}

void TwoDSceneXMLParser::loadXMLFile(const std::string& filename,
                                     std::vector<char>& xmlchars,
                                     rapidxml::xml_document<>& doc) {
  // Attempt to read the text from the user-specified xml file
  std::string filecontents;
  if (!loadTextFileIntoString(filename, filecontents)) {
    std::cerr << "\033[31;1mERROR IN TWODSCENEXMLPARSER:\033[m XML scene file "
              << filename << ". Failed to read file." << std::endl;
    exit(1);
  }

  // Copy string into an array of characters for the xml parser
  for (int i = 0; i < (int)filecontents.size(); ++i)
    xmlchars.push_back(filecontents[i]);
  xmlchars.push_back('\0');

  // Initialize the xml parser with the character vector
  doc.parse<0>(&xmlchars[0]);
}

bool TwoDSceneXMLParser::loadTextFileIntoString(const std::string& filename,
                                                std::string& filecontents) {
  // Attempt to open the text file for reading
  std::ifstream textfile(filename.c_str(), std::ifstream::in);
  if (!textfile) return false;

  // Read the entire file into a single string
  std::string line;
  while (getline(textfile, line)) filecontents.append(line);

  textfile.close();

  return true;
}

void TwoDSceneXMLParser::loadSimulationType(rapidxml::xml_node<>* node,
                                            std::string& simtype) {
  assert(node != NULL);
  rapidxml::xml_node<>* nd = node->first_node("simtype");

  if (node->first_node("simtype"))
    if (node->first_node("simtype")->first_attribute("type"))
      simtype = node->first_node("simtype")->first_attribute("type")->value();
}

template <int DIM>
void TwoDSceneXMLParser::loadFlow(rapidxml::xml_node<>* node,
                                  TwoDScene<DIM>& twodscene) {
  for (rapidxml::xml_node<>* nd = node->first_node("flow"); nd;
       nd = nd->next_sibling("flow")) {
    std::vector<int> particle_indices;
    std::vector<scalar> particle_eta;
    std::vector<unsigned char> particle_state;

    LIQUID_SIM_TYPE lst = LST_SHALLOW;

    rapidxml::xml_attribute<>* typend = nd->first_attribute("type");
    if (typend != NULL) {
      std::string handlertype(typend->value());
      if (handlertype == "shallow")
        lst = LST_SHALLOW;
      else {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Invalid flow "
                     "'type' attribute specified. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    for (rapidxml::xml_node<>* subnd = nd->first_node("node"); subnd;
         subnd = subnd->next_sibling("node")) {
      int idx = -1;
      if (subnd->first_attribute("i")) {
        std::string attribute(subnd->first_attribute("i")->value());
        if (!stringutils::extractFromString(attribute, idx)) {
          std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to "
                       "parse value of i attribute for film flow. Value must "
                       "be numeric. Exiting."
                    << std::endl;
          exit(1);
        }
      }

      if (idx == -1) continue;

      scalar eta = -1e+20;
      if (subnd->first_attribute("eta")) {
        std::string attribute(subnd->first_attribute("eta")->value());
        if (!stringutils::extractFromString(attribute, eta)) {
          std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to "
                       "parse value of eta attribute for film flow. Value must "
                       "be numeric. Exiting."
                    << std::endl;
          exit(1);
        }
      }

      if (eta < 0) continue;

      bool isSource = false;
      if (subnd->first_attribute("source")) {
        std::string attribute(subnd->first_attribute("source")->value());
        if (!stringutils::extractFromString(attribute, isSource)) {
          std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to "
                       "parse value of source attribute for film flow. Value "
                       "must be boolean. Exiting."
                    << std::endl;
          exit(1);
        }
      }

      particle_indices.push_back(idx);
      particle_eta.push_back(eta);
      particle_state.push_back(isSource);
    }

    HairFlow<DIM>* flow = NULL;

    switch (lst) {
      case LST_SHALLOW:
        flow = new CylindricalShallowFlow<DIM>(
            &twodscene, particle_indices,
            Eigen::Map<VectorXs>(&particle_eta[0], particle_eta.size()),
            particle_state);
        break;
      default:
        break;
    }

    if (flow) {
      twodscene.insertFilmFlow(flow);
    }
  }
}

template <int DIM>
void TwoDSceneXMLParser::loadParticleEdges(rapidxml::xml_node<>* node,
                                           TwoDScene<DIM>& twodscene) {
  assert(node != NULL);

  // Count the number of particles
  int numparticles = 0;
  for (rapidxml::xml_node<>* nd = node->first_node("particle"); nd;
       nd = nd->next_sibling("particle"))
    ++numparticles;

  int numedges = 0;
  for (rapidxml::xml_node<>* nd = node->first_node("edge"); nd;
       nd = nd->next_sibling("edge"))
    ++numedges;

  twodscene.resizeSystem(numparticles, numedges);

  // std::cout << "Num particles " << numparticles << std::endl;

  std::vector<std::string>& tags = twodscene.getParticleTags();

  int particle = 0;
  for (rapidxml::xml_node<>* nd = node->first_node("particle"); nd;
       nd = nd->next_sibling("particle")) {
    // Extract the particle's initial position
    Vectors<DIM> pos;
    if (nd->first_attribute("px")) {
      std::string attribute(nd->first_attribute("px")->value());
      if (!stringutils::extractFromString(attribute, pos(0))) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of px attribute for particle "
                  << particle << ". Value must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    } else {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                   "px attribute for particle "
                << particle << ". Exiting." << std::endl;
      exit(1);
    }

    if (nd->first_attribute("py")) {
      std::string attribute(nd->first_attribute("py")->value());
      if (!stringutils::extractFromString(attribute, pos(1))) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of py attribute for particle "
                  << particle << ". Value must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    } else {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                   "py attribute for particle "
                << particle << ". Exiting." << std::endl;
      exit(1);
    }

    if (DIM == 3) {
      if (nd->first_attribute("pz")) {
        std::string attribute(nd->first_attribute("pz")->value());
        if (!stringutils::extractFromString(attribute, pos(2))) {
          std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to "
                       "parse value of pz attribute for particle "
                    << particle << ". Value must be numeric. Exiting."
                    << std::endl;
          exit(1);
        }
      } else {
        pos(2) = 0.0;
      }
    }

    twodscene.setPosition(particle, pos);

    // Extract the particle's initial velocity
    Vectors<DIM> vel;
    if (nd->first_attribute("vx")) {
      std::string attribute(nd->first_attribute("vx")->value());
      if (!stringutils::extractFromString(attribute, vel(0))) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of vx attribute for particle "
                  << particle << ". Value must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    } else {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                   "vx attribute for particle "
                << particle << ". Exiting." << std::endl;
      exit(1);
    }

    if (nd->first_attribute("vy")) {
      std::string attribute(nd->first_attribute("vy")->value());
      if (!stringutils::extractFromString(attribute, vel(1))) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of vy attribute for particle "
                  << particle << ". Value must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    } else {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                   "vy attribute for particle "
                << particle << ". Exiting." << std::endl;
      exit(1);
    }

    if (DIM == 3) {
      if (nd->first_attribute("vz")) {
        std::string attribute(nd->first_attribute("vz")->value());
        if (!stringutils::extractFromString(attribute, vel(2))) {
          std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to "
                       "parse value of vz attribute for particle "
                    << particle << ". Value must be numeric. Exiting."
                    << std::endl;
          exit(1);
        }
      } else {
        vel(2) = 0.0;
      }
    }

    twodscene.setVelocity(particle, vel);

    // Extract the particle's mass
    scalar mass = std::numeric_limits<scalar>::signaling_NaN();
    if (nd->first_attribute("m")) {
      std::string attribute(nd->first_attribute("m")->value());
      if (!stringutils::extractFromString(attribute, mass)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of m attribute for particle "
                  << particle << ". Value must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    } else {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse m "
                   "attribute for particle "
                << particle << ". Exiting." << std::endl;
      exit(1);
    }
    twodscene.setMass(particle, mass);

    // Determine if the particle is fixed
    bool fixed;
    if (nd->first_attribute("fixed")) {
      std::string attribute(nd->first_attribute("fixed")->value());
      if (!stringutils::extractFromString(attribute, fixed)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of fixed attribute for particle "
                  << particle << ". Value must be boolean. Exiting."
                  << std::endl;
        exit(1);
      }
    } else {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                   "fixed attribute for particle "
                << particle << ". Exiting." << std::endl;
      exit(1);
    }
    twodscene.setFixed(particle, fixed);

    // Extract the particle's radius, if present
    scalar radius = 0.1;
    if (nd->first_attribute("radius")) {
      std::string attribute(nd->first_attribute("radius")->value());
      if (!stringutils::extractFromString(attribute, radius)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "radius attribute for particle "
                  << particle << ". Value must be scalar. Exiting."
                  << std::endl;
        exit(1);
      }
    }
    twodscene.setRadius(particle, radius);

    // Extract the particle's tag, if present
    if (nd->first_attribute("tag")) {
      std::string tag(nd->first_attribute("tag")->value());
      tags[particle] = tag;
    }

    int group = 0;
    if (nd->first_attribute("group")) {
      std::string attribute(nd->first_attribute("group")->value());
      if (!stringutils::extractFromString(attribute, group)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of group attribute for particle "
                  << particle << ". Value must be integer. Exiting."
                  << std::endl;
        exit(1);
      }
    }
    twodscene.setScriptedGroup(particle, group);

    // std::cout << "Particle: " << particle << "    x: " << pos.transpose() <<
    // "   v: " << vel.transpose() << "   m: " << mass << "   fixed: " << fixed
    // << std::endl; std::cout << tags[particle] << std::endl;

    ++particle;
  }

  const VectorXs& x = twodscene.getX();

  int edge = 0;
  for (rapidxml::xml_node<>* nd = node->first_node("edge"); nd;
       nd = nd->next_sibling("edge")) {
    std::pair<int, int> newedge;
    if (nd->first_attribute("i")) {
      std::string attribute(nd->first_attribute("i")->value());
      if (!stringutils::extractFromString(attribute, newedge.first)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of i attribute for edge "
                  << edge << ". Value must be integer. Exiting." << std::endl;
        exit(1);
      }
    } else {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                   "value of i attribute for edge "
                << edge << ". Exiting." << std::endl;
      exit(1);
    }

    if (nd->first_attribute("j")) {
      std::string attribute(nd->first_attribute("j")->value());
      if (!stringutils::extractFromString(attribute, newedge.second)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of j attribute for edge "
                  << edge << ". Value must be integer. Exiting." << std::endl;
        exit(1);
      }
    } else {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                   "value of j attribute for edge "
                << edge << ". Exiting." << std::endl;
      exit(1);
    }

    scalar radius = 0.015;
    if (nd->first_attribute("radius")) {
      std::string attribute(nd->first_attribute("radius")->value());
      if (!stringutils::extractFromString(attribute, radius)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "radius attribute for edge "
                  << edge << ". Value must be scalar. Exiting." << std::endl;
        exit(1);
      }
    }

    scalar poisson = 0.377;
    if (nd->first_attribute("poisson")) {
      std::string attribute(nd->first_attribute("poisson")->value());
      if (!stringutils::extractFromString(attribute, poisson)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "poisson attribute for edge "
                  << edge << ". Value must be scalar. Exiting." << std::endl;
        exit(1);
      }
    }

    // std::cout << "Edge: " << edge << "    i: " << newedge.first << "   j: "
    // << newedge.second << std::endl;

    twodscene.setEdge(edge, newedge, radius);
    twodscene.setEdgePoissonRatio(edge, poisson);
    twodscene.setEdgeRestLength(
        edge, (x.segment<DIM>(twodscene.getDof(newedge.first)) -
               x.segment<DIM>(twodscene.getDof(newedge.second)))
                  .norm());

    ++edge;
  }
}

void TwoDSceneXMLParser::loadStrandParticleEdges(rapidxml::xml_node<>* node,
                                                 TwoDScene<3>& scene) {
  assert(node != NULL);

  // Count the number of particles, edges, and strands
  int numstrands = 0;
  int numparticles = 0;
  int numedges = 0;
  for (rapidxml::xml_node<>* nd = node->first_node("StrandMaterialForces"); nd;
       nd = nd->next_sibling("StrandMaterialForces")) {
    ++numstrands;
    int edgeCount = 0;
    for (rapidxml::xml_node<>* subnd = nd->first_node("particle"); subnd;
         subnd = subnd->next_sibling("particle")) {
      ++numparticles;
      ++edgeCount;
    }
    numedges += edgeCount - 1;
  }

  for (rapidxml::xml_node<>* nd = node->first_node("ParameterizedStrand"); nd;
       nd = nd->next_sibling("ParameterizedStrand")) {
    int nv = 0;
    if (nd->first_attribute("nv")) {
      std::string attribute(nd->first_attribute("nv")->value());
      if (!stringutils::extractFromString(attribute, nv)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "nv for strand "
                  << numstrands << ". Value must be scalar. Exiting."
                  << std::endl;
        exit(1);
      }
    }
    if (nv < 2) continue;
    numparticles += nv;
    numedges += nv - 1;
    ++numstrands;
  }

  scene.resizeSystem(numparticles, numedges, numstrands);

  // compute vertex ID to DoF mapping
  std::vector<int> vert_to_dof(numparticles);
  VectorXi dofVars(numparticles * 4 - numstrands + 1);
  VectorXi dofVerts(numparticles * 4 - numstrands + 1);
  Vector4i dofs = Vector4i(0, 1, 2, 3);

  std::vector<bool> tipVerts(numparticles);
  int dof = 0;
  numparticles = 0;
  for (rapidxml::xml_node<>* nd = node->first_node("StrandMaterialForces"); nd;
       nd = nd->next_sibling("StrandMaterialForces")) {
    for (rapidxml::xml_node<>* subnd = nd->first_node("particle"); subnd;
         subnd = subnd->next_sibling("particle")) {
      dofVars.segment<4>(dof) = dofs;
      dofVerts.segment<4>(dof).setConstant(numparticles);
      vert_to_dof[numparticles++] = dof;
      dof += 4;
    }
    --dof;
    tipVerts[numparticles - 1] = true;
  }
  for (rapidxml::xml_node<>* nd = node->first_node("ParameterizedStrand"); nd;
       nd = nd->next_sibling("ParameterizedStrand")) {
    int nv = 0;
    if (nd->first_attribute("nv")) {
      std::string attribute(nd->first_attribute("nv")->value());
      if (!stringutils::extractFromString(attribute, nv)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "nv for strand "
                  << numstrands << ". Value must be scalar. Exiting."
                  << std::endl;
        exit(1);
      }
    }
    if (nv < 2) continue;
    for (int i = 0; i < nv; ++i) {
      dofVars.segment<4>(dof) = dofs;
      dofVerts.segment<4>(dof).setConstant(numparticles);
      vert_to_dof[numparticles++] = dof;
      dof += 4;
    }
    --dof;
    tipVerts[numparticles - 1] = true;
  }
  dofVars.conservativeResize(numparticles * 4 - numstrands);
  dofVerts.conservativeResize(numparticles * 4 - numstrands);
  scene.setVertToDoFMap(vert_to_dof, dofVars, tipVerts, dofVerts);

  std::vector<std::string>& tags = scene.getParticleTags();
  int globalStrandID = 0;
  std::vector<std::vector<int> > particle_indices_vec;
  std::vector<std::vector<scalar> > particle_eta_vec;
  std::vector<std::vector<unsigned char> > particle_state_vec;
  std::vector<LIQUID_SIM_TYPE> liquidtype_vec;

  std::vector<Vector3s> particle_pos;
  Vector3s pos;
  Vector3s vel;
  int vtx = 0;
  int edx = 0;
  for (rapidxml::xml_node<>* nd = node->first_node("StrandMaterialForces"); nd;
       nd = nd->next_sibling("StrandMaterialForces")) {
    int paramsIndex = -1;
    if (nd->first_attribute("params")) {
      std::string attribute(nd->first_attribute("params")->value());
      if (!stringutils::extractFromString(attribute, paramsIndex)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of params (StrandParameters index) for StrandForce "
                  << globalStrandID << ". Value must be integer. Exiting."
                  << std::endl;
        exit(1);
      }
    } else {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                   "value of params (StrandParameters index) for StrandForce "
                << globalStrandID << ". Exiting." << std::endl;
      exit(1);
    }

    LIQUID_SIM_TYPE lst = LST_SHALLOW;
    rapidxml::xml_attribute<>* typend = nd->first_attribute("type");
    if (typend != NULL) {
      std::string handlertype(typend->value());
      if (handlertype == "shallow")
        lst = LST_SHALLOW;
      else {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Invalid flow "
                     "'type' attribute specified. Exiting."
                  << std::endl;
        exit(1);
      }
    }
    liquidtype_vec.push_back(lst);

    std::vector<scalar> particle_eta;
    std::vector<unsigned char> particle_state;
    std::vector<int> particle_indices;
    particle_indices.clear();
    particle_pos.clear();
    int globalvtx = vtx;
    for (rapidxml::xml_node<>* subnd = nd->first_node("particle"); subnd;
         subnd = subnd->next_sibling("particle")) {
      particle_indices.push_back(vtx);

      scalar eta = 0.0;
      if (subnd->first_attribute("eta")) {
        std::string attribute(subnd->first_attribute("eta")->value());
        if (!stringutils::extractFromString(attribute, eta)) {
          std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to "
                       "parse value of eta attribute for film flow. Value must "
                       "be numeric. Exiting."
                    << std::endl;
          exit(1);
        }
      }
      particle_eta.push_back(eta);

      bool isSource = false;
      if (subnd->first_attribute("source")) {
        std::string attribute(subnd->first_attribute("source")->value());
        if (!stringutils::extractFromString(attribute, isSource)) {
          std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to "
                       "parse value of source attribute for film flow. Value "
                       "must be boolean. Exiting."
                    << std::endl;
          exit(1);
        }
      }
      particle_state.push_back(isSource);

      // Extract the particle's initial position
      pos.setZero();
      if (subnd->first_attribute("x")) {
        std::string position(subnd->first_attribute("x")->value());
        if (!stringutils::readList(position, ' ', pos)) {
          std::cerr << "Failed to load x, y, and z positions for particle "
                    << vtx << std::endl;
          exit(1);
        }
      } else {
        std::cerr
            << "Failed to find x, y, and z position attributes for particle "
            << vtx << std::endl;
        exit(1);
      }
      scene.setPosition(vtx, pos);
      particle_pos.push_back(pos);

      // Extract the particle's initial velocity
      vel.setZero();
      if (subnd->first_attribute("v")) {
        std::string velocity(subnd->first_attribute("v")->value());
        if (!stringutils::readList(velocity, ' ', vel)) {
          std::cerr << "Failed to load x, y, and z positions for particle "
                    << vtx << std::endl;
          exit(1);
        }
      }
      scene.setVelocity(vtx, vel);

      // Determine if the particle is fixed
      bool fixed = false;
      if (subnd->first_attribute("fixed")) {
        std::string attribute(subnd->first_attribute("fixed")->value());
        if (!stringutils::extractFromString(attribute, fixed)) {
          std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to "
                       "parse value of fixed attribute for particle "
                    << vtx << ". Value must be boolean. Exiting." << std::endl;
          exit(1);
        }
      }
      scene.setFixed(vtx, fixed);

      // Extract the particle's tag, if present
      if (subnd->first_attribute("tag")) {
        std::string tag(subnd->first_attribute("tag")->value());
        tags[vtx] = tag;
      }

      int group = 0;
      if (subnd->first_attribute("group")) {
        std::string attribute(subnd->first_attribute("group")->value());
        if (!stringutils::extractFromString(attribute, group)) {
          std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to "
                       "parse value of group attribute for particle "
                    << vtx << ". Value must be integer. Exiting." << std::endl;
          exit(1);
        }
      }
      scene.setScriptedGroup(vtx, group);

      if (vtx != globalvtx) {
        std::pair<int, int> newedge(vtx - 1, vtx);
        scene.setEdge(edx, newedge);
        scene.setEdgeRestLength(edx, (scene.getPosition(newedge.first) -
                                      scene.getPosition(newedge.second))
                                         .norm());
        ++edx;
      }

      ++vtx;
    }

    particle_indices_vec.push_back(particle_indices);
    particle_eta_vec.push_back(particle_eta);
    particle_state_vec.push_back(particle_state);

    scene.insertForce(new StrandForce(&scene, particle_indices, paramsIndex,
                                      globalStrandID++));
    scene.insertStrandEquilibriumParameters(
        new StrandEquilibriumParameters(particle_pos, 0., 0., 0., 0., false));
  }

  for (rapidxml::xml_node<>* nd = node->first_node("ParameterizedStrand"); nd;
       nd = nd->next_sibling("ParameterizedStrand")) {
    int paramsIndex = -1;
    if (nd->first_attribute("params")) {
      std::string attribute(nd->first_attribute("params")->value());
      if (!stringutils::extractFromString(attribute, paramsIndex)) {
        std::cerr
            << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse value "
               "of params (StrandParameters index) for ParameterizedStrand "
            << globalStrandID << ". Value must be integer. Exiting."
            << std::endl;
        exit(1);
      }
    } else {
      std::cerr
          << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse value "
             "of params (StrandParameters index) for ParameterizedStrand "
          << globalStrandID << ". Exiting." << std::endl;
      exit(1);
    }

    LIQUID_SIM_TYPE lst = LST_SHALLOW;
    rapidxml::xml_attribute<>* typend = nd->first_attribute("type");
    if (typend != NULL) {
      std::string handlertype(typend->value());
      if (handlertype == "shallow")
        lst = LST_SHALLOW;
      else {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Invalid flow "
                     "'type' attribute specified. Exiting."
                  << std::endl;
        exit(1);
      }
    }
    liquidtype_vec.push_back(lst);

    Vec3 initnorm;
    Vec3 startpoint;
    scalar length;
    scalar curl_radius = 0.0;
    scalar curl_density = 1e+63;
    int nv;

    rapidxml::xml_attribute<>* attr = nd->first_attribute("nx");
    if (attr != NULL) {
      std::string attribute(attr->value());
      if (!stringutils::extractFromString(attribute, initnorm(0))) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of nx for ParameterizedStrand "
                  << globalStrandID << ". Value must be real. Exiting."
                  << std::endl;
        exit(1);
      }
    } else {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                   "value of nx for ParameterizedStrand "
                << globalStrandID << ". Value must be real. Exiting."
                << std::endl;
      exit(1);
    }

    attr = nd->first_attribute("ny");
    if (attr != NULL) {
      std::string attribute(attr->value());
      if (!stringutils::extractFromString(attribute, initnorm(1))) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of ny for ParameterizedStrand "
                  << globalStrandID << ". Value must be real. Exiting."
                  << std::endl;
        exit(1);
      }
    } else {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                   "value of ny for ParameterizedStrand "
                << globalStrandID << ". Value must be real. Exiting."
                << std::endl;
      exit(1);
    }

    attr = nd->first_attribute("nz");
    if (attr != NULL) {
      std::string attribute(attr->value());
      if (!stringutils::extractFromString(attribute, initnorm(2))) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of nz for ParameterizedStrand "
                  << globalStrandID << ". Value must be real. Exiting."
                  << std::endl;
        exit(1);
      }
    } else {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                   "value of nz for ParameterizedStrand "
                << globalStrandID << ". Value must be real. Exiting."
                << std::endl;
      exit(1);
    }

    attr = nd->first_attribute("px");
    if (attr != NULL) {
      std::string attribute(attr->value());
      if (!stringutils::extractFromString(attribute, startpoint(0))) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of px for ParameterizedStrand "
                  << globalStrandID << ". Value must be real. Exiting."
                  << std::endl;
        exit(1);
      }
    } else {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                   "value of px for ParameterizedStrand "
                << globalStrandID << ". Value must be real. Exiting."
                << std::endl;
      exit(1);
    }

    attr = nd->first_attribute("py");
    if (attr != NULL) {
      std::string attribute(attr->value());
      if (!stringutils::extractFromString(attribute, startpoint(1))) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of py for ParameterizedStrand "
                  << globalStrandID << ". Value must be real. Exiting."
                  << std::endl;
        exit(1);
      }
    } else {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                   "value of py for ParameterizedStrand "
                << globalStrandID << ". Value must be real. Exiting."
                << std::endl;
      exit(1);
    }

    attr = nd->first_attribute("pz");
    if (attr != NULL) {
      std::string attribute(attr->value());
      if (!stringutils::extractFromString(attribute, startpoint(2))) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of pz for ParameterizedStrand "
                  << globalStrandID << ". Value must be real. Exiting."
                  << std::endl;
        exit(1);
      }
    } else {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                   "value of pz for ParameterizedStrand "
                << globalStrandID << ". Value must be real. Exiting."
                << std::endl;
      exit(1);
    }

    attr = nd->first_attribute("length");
    if (attr != NULL) {
      std::string attribute(attr->value());
      if (!stringutils::extractFromString(attribute, length)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of length for ParameterizedStrand "
                  << globalStrandID << ". Value must be real. Exiting."
                  << std::endl;
        exit(1);
      }
    } else {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                   "value of length for ParameterizedStrand "
                << globalStrandID << ". Value must be real. Exiting."
                << std::endl;
      exit(1);
    }

    attr = nd->first_attribute("nv");
    if (attr != NULL) {
      std::string attribute(attr->value());
      if (!stringutils::extractFromString(attribute, nv)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of nv for ParameterizedStrand "
                  << globalStrandID << ". Value must be integer. Exiting."
                  << std::endl;
        exit(1);
      }
    } else {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                   "value of nv for ParameterizedStrand "
                << globalStrandID << ". Value must be integer. Exiting."
                << std::endl;
      exit(1);
    }

    attr = nd->first_attribute("curlradius");
    if (attr != NULL) {
      std::string attribute(attr->value());
      if (!stringutils::extractFromString(attribute, curl_radius)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of curlradius for ParameterizedStrand "
                  << globalStrandID << ". Value must be real. Exiting."
                  << std::endl;
        exit(1);
      }
    } else {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                   "value of curlradius for ParameterizedStrand "
                << globalStrandID << ". Value must be real. Exiting."
                << std::endl;
      exit(1);
    }

    attr = nd->first_attribute("curldensity");
    if (attr != NULL) {
      std::string attribute(attr->value());
      if (!stringutils::extractFromString(attribute, curl_density)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of curldensity for ParameterizedStrand "
                  << globalStrandID << ". Value must be real. Exiting."
                  << std::endl;
        exit(1);
      }
    } else {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                   "value of curldensity for ParameterizedStrand "
                << globalStrandID << ". Value must be real. Exiting."
                << std::endl;
      exit(1);
    }

    scalar dL = length / (double)(nv - 1);
    std::vector<Vec3> strand_vertices;

    genCurlyHair(initnorm, startpoint, dL, nv, strand_vertices, curl_radius,
                 curl_density, dL);

    int global_vtx = vtx;
    for (int i = 0; i < nv; ++i) {
      scene.setPosition(global_vtx + i, strand_vertices[i]);
    }

    std::vector<scalar> particle_eta;
    std::vector<unsigned char> particle_state;
    std::vector<int> particle_indices;
    particle_indices.clear();
    particle_pos.clear();
    for (rapidxml::xml_node<>* subnd = nd->first_node("particle"); subnd;
         subnd = subnd->next_sibling("particle")) {
      particle_indices.push_back(vtx);

      scalar eta = 0.0;
      if (subnd->first_attribute("eta")) {
        std::string attribute(subnd->first_attribute("eta")->value());
        if (!stringutils::extractFromString(attribute, eta)) {
          std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to "
                       "parse value of eta attribute for film flow. Value must "
                       "be numeric. Exiting."
                    << std::endl;
          exit(1);
        }
      }
      particle_eta.push_back(eta);

      bool isSource = false;
      if (subnd->first_attribute("source")) {
        std::string attribute(subnd->first_attribute("source")->value());
        if (!stringutils::extractFromString(attribute, isSource)) {
          std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to "
                       "parse value of source attribute for film flow. Value "
                       "must be boolean. Exiting."
                    << std::endl;
          exit(1);
        }
      }
      particle_state.push_back(isSource);

      // Extract the particle's initial velocity
      vel.setZero();
      if (subnd->first_attribute("v")) {
        std::string velocity(subnd->first_attribute("v")->value());
        if (!stringutils::readList(velocity, ' ', vel)) {
          std::cerr << "Failed to load x, y, and z positions for particle "
                    << vtx << std::endl;
          exit(1);
        }
      }
      scene.setVelocity(vtx, vel);

      // Determine if the particle is fixed
      bool fixed = false;
      if (subnd->first_attribute("fixed")) {
        std::string attribute(subnd->first_attribute("fixed")->value());
        if (!stringutils::extractFromString(attribute, fixed)) {
          std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to "
                       "parse value of fixed attribute for particle "
                    << vtx << ". Value must be boolean. Exiting." << std::endl;
          exit(1);
        }
      }
      scene.setFixed(vtx, fixed);

      // Extract the particle's tag, if present
      if (subnd->first_attribute("tag")) {
        std::string tag(subnd->first_attribute("tag")->value());
        tags[vtx] = tag;
      }

      int group = 0;
      if (subnd->first_attribute("group")) {
        std::string attribute(subnd->first_attribute("group")->value());
        if (!stringutils::extractFromString(attribute, group)) {
          std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to "
                       "parse value of group attribute for particle "
                    << vtx << ". Value must be integer. Exiting." << std::endl;
          exit(1);
        }
      }
      scene.setScriptedGroup(vtx, group);

      if (vtx != global_vtx) {
        std::pair<int, int> newedge(vtx - 1, vtx);
        scene.setEdge(edx, newedge);
        scene.setEdgeRestLength(edx, (scene.getPosition(newedge.first) -
                                      scene.getPosition(newedge.second))
                                         .norm());
        ++edx;
      }

      ++vtx;
    }

    while (vtx < global_vtx + nv) {
      particle_indices.push_back(vtx);
      particle_eta.push_back(0.0);
      particle_state.push_back(false);
      scene.setVelocity(vtx, Vec3::Zero());
      scene.setFixed(vtx, false);
      scene.setScriptedGroup(vtx, 0);
      if (vtx != global_vtx) {
        std::pair<int, int> newedge(vtx - 1, vtx);
        scene.setEdge(edx, newedge);
        scene.setEdgeRestLength(edx, (scene.getPosition(newedge.first) -
                                      scene.getPosition(newedge.second))
                                         .norm());
        ++edx;
      }
      ++vtx;
    }

    particle_indices_vec.push_back(particle_indices);
    particle_eta_vec.push_back(particle_eta);
    particle_state_vec.push_back(particle_state);

    scene.insertForce(new StrandForce(&scene, particle_indices, paramsIndex,
                                      globalStrandID++));
    scene.insertStrandEquilibriumParameters(new StrandEquilibriumParameters(
        strand_vertices, curl_radius, curl_density, dL, dL, true));
  }

  scene.computeMassesAndRadiiFromStrands();

  scene.getFilmFlows().resize(particle_indices_vec.size(), NULL);

  threadutils::thread_pool::ParallelFor(
      0, (int)particle_indices_vec.size(), [&](int p) {
        scene.getFilmFlows()[p] = new CylindricalShallowFlow<3>(
            &scene, particle_indices_vec[p],
            Eigen::Map<VectorXs>(&particle_eta_vec[p][0],
                                 particle_eta_vec[p].size()),
            particle_state_vec[p]);
      });
}

void TwoDSceneXMLParser::loadSceneTag(rapidxml::xml_node<>* node,
                                      std::string& scenetag) {
  assert(node != NULL);

  if (node->first_node("scenetag")) {
    if (node->first_node("scenetag")->first_attribute("tag")) {
      scenetag = node->first_node("scenetag")->first_attribute("tag")->value();
    } else {
      std::cerr
          << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse value "
             "of tag attribute for scenetag. Value must be string. Exiting."
          << std::endl;
      exit(1);
    }
  }
}

template <int DIM>
void TwoDSceneXMLParser::loadLinearSpringForces(rapidxml::xml_node<>* node,
                                                TwoDScene<DIM>& twodscene) {
  assert(node != NULL);

  // Extract the edge the force acts across
  int forcenum = 0;
  for (rapidxml::xml_node<>* nd = node->first_node("linearspringforce"); nd;
       nd = nd->next_sibling("linearspringforce")) {
    int edge = -1;

    if (nd->first_attribute("edge")) {
      std::string attribute(nd->first_attribute("edge")->value());
      if (!stringutils::extractFromString(attribute, edge)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of edge attribute for springforce "
                  << forcenum << ". Value must be integer. Exiting."
                  << std::endl;
        exit(1);
      }
    } else {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                   "value of edge attribute for springforce "
                << forcenum << ". Exiting." << std::endl;
      exit(1);
    }

    std::pair<int, int> newedge(twodscene.getEdge(edge));

    // Extract the spring stiffness
    scalar k = std::numeric_limits<scalar>::signaling_NaN();
    if (nd->first_attribute("k")) {
      std::string attribute(nd->first_attribute("k")->value());
      if (!stringutils::extractFromString(attribute, k)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of k attribute for springforce "
                  << forcenum << ". Value must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    } else {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse k "
                   "attribute for springforce "
                << forcenum << ". Exiting." << std::endl;
      exit(1);
    }

    // Extract the spring rest length
    scalar l0 = std::numeric_limits<scalar>::signaling_NaN();
    if (nd->first_attribute("l0")) {
      std::string attribute(nd->first_attribute("l0")->value());
      if (!stringutils::extractFromString(attribute, l0)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of l0 attribute for springforce "
                  << forcenum << ". Value must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    } else {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                   "l0 attribute for springforce "
                << forcenum << ". Exiting." << std::endl;
      exit(1);
    }

    // Extract the optional damping coefficient
    scalar b = 0.0;
    if (nd->first_attribute("b")) {
      std::string attribute(nd->first_attribute("b")->value());
      if (!stringutils::extractFromString(attribute, b)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of b attribute for springforce "
                  << forcenum << ". Value must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    // std::cout << "Springforce: " << forcenum << "    i: " << newedge.first <<
    // "   j: " << newedge.second << "   k: " << k << "   l0: " << l0 <<
    // std::endl;

    twodscene.insertForce(
        new LinearSpringForce<DIM>(&twodscene, newedge, k, l0, b));
    twodscene.setEdgeRestLength(edge, l0);

    ++forcenum;
  }
}

template <int DIM>
void TwoDSceneXMLParser::loadSpringForces(rapidxml::xml_node<>* node,
                                          TwoDScene<DIM>& twodscene) {
  assert(node != NULL);

  // Extract the edge the force acts across
  int forcenum = 0;
  for (rapidxml::xml_node<>* nd = node->first_node("springforce"); nd;
       nd = nd->next_sibling("springforce")) {
    int edge = -1;

    if (nd->first_attribute("edge")) {
      std::string attribute(nd->first_attribute("edge")->value());
      if (!stringutils::extractFromString(attribute, edge)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of edge attribute for springforce "
                  << forcenum << ". Value must be integer. Exiting."
                  << std::endl;
        exit(1);
      }
    } else {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                   "value of edge attribute for springforce "
                << forcenum << ". Exiting." << std::endl;
      exit(1);
    }

    std::pair<int, int> newedge(twodscene.getEdge(edge));

    // Extract the spring stiffness
    scalar k = std::numeric_limits<scalar>::signaling_NaN();
    if (nd->first_attribute("k")) {
      std::string attribute(nd->first_attribute("k")->value());
      if (!stringutils::extractFromString(attribute, k)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of k attribute for springforce "
                  << forcenum << ". Value must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    } else {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse k "
                   "attribute for springforce "
                << forcenum << ". Exiting." << std::endl;
      exit(1);
    }

    // Extract the spring rest length
    scalar l0 = std::numeric_limits<scalar>::signaling_NaN();
    if (nd->first_attribute("l0")) {
      std::string attribute(nd->first_attribute("l0")->value());
      if (!stringutils::extractFromString(attribute, l0)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of l0 attribute for springforce "
                  << forcenum << ". Value must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    } else {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                   "l0 attribute for springforce "
                << forcenum << ". Exiting." << std::endl;
      exit(1);
    }

    // Extract the optional damping coefficient
    scalar b = 0.0;
    if (nd->first_attribute("b")) {
      std::string attribute(nd->first_attribute("b")->value());
      if (!stringutils::extractFromString(attribute, b)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of b attribute for springforce "
                  << forcenum << ". Value must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    // std::cout << "Springforce: " << forcenum << "    i: " << newedge.first <<
    // "   j: " << newedge.second << "   k: " << k << "   l0: " << l0 <<
    // std::endl;

    twodscene.insertForce(new SpringForce<DIM>(newedge, k, l0, &twodscene, b));

    if (DIM == 2) twodscene.setEdgeRestLength(edge, l0);

    ++forcenum;
  }

  // SpringForce( const std::pair<int,int>& endpoints, const scalar& k, const
  // scalar& l0 )
}

template <int DIM>
void TwoDSceneXMLParser::loadDragDampingForces(rapidxml::xml_node<>* node,
                                               TwoDScene<DIM>& twodscene) {
  assert(node != NULL);

  int forcenum = 0;
  for (rapidxml::xml_node<>* nd = node->first_node("dragdamping"); nd;
       nd = nd->next_sibling("dragdamping")) {
    Vector2s constforce;
    constforce.setConstant(std::numeric_limits<scalar>::signaling_NaN());

    // Extract the linear damping coefficient
    scalar b = std::numeric_limits<scalar>::signaling_NaN();
    if (nd->first_attribute("b")) {
      std::string attribute(nd->first_attribute("b")->value());
      if (!stringutils::extractFromString(attribute, b)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of b attribute for dragdamping "
                  << forcenum << ". Value must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    } else {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse b "
                   "attribute for dragdamping "
                << forcenum << ". Exiting." << std::endl;
      exit(1);
    }

    int i = std::numeric_limits<int>::signaling_NaN();
    if (nd->first_attribute("i")) {
      std::string attribute(nd->first_attribute("i")->value());
      if (!stringutils::extractFromString(attribute, i)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of i attribute for dragdamping "
                  << forcenum << ". Value must be integer. Exiting."
                  << std::endl;
        exit(1);
      }
    } else {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse i "
                   "attribute for dragdamping "
                << forcenum << ". Exiting." << std::endl;
      exit(1);
    }

    // std::cout << "x: " << constforce.transpose() << std::endl;

    twodscene.insertForce(new DragDampingForce<DIM>(twodscene, b, i));

    ++forcenum;
  }
}

template <int DIM>
void TwoDSceneXMLParser::loadIntegrator(rapidxml::xml_node<>* node,
                                        TwoDScene<DIM>& twodscene,
                                        SceneStepper<DIM>** scenestepper) {
  assert(node != NULL);

  // Attempt to locate the integrator node
  rapidxml::xml_node<>* nd = node->first_node("integrator");
  if (nd == NULL) {
    std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m No integrator "
                 "specified. Exiting."
              << std::endl;
    exit(1);
  }

  // Attempt to load the integrator type
  rapidxml::xml_attribute<>* typend = nd->first_attribute("type");
  if (typend == NULL) {
    std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m No integrator "
                 "'type' attribute specified. Exiting."
              << std::endl;
    exit(1);
  }
  std::string integratortype(typend->value());

  if (integratortype == "direct-solver") {
    rapidxml::xml_attribute<>* dtnd = nd->first_attribute("maxiters");
    int max_iters = 0;
    if (dtnd) {
      if (!stringutils::extractFromString(std::string(dtnd->value()),
                                          max_iters)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "'maxiters' attribute for integrator. Value must be "
                     "integer. Exiting."
                  << std::endl;
        exit(1);
      }
    }
    dtnd = nd->first_attribute("criterion");
    scalar criterion = 1e-8;
    if (dtnd) {
      if (!stringutils::extractFromString(std::string(dtnd->value()),
                                          criterion)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "'criterion' attribute for integrator. Value must be "
                     "numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }
    *scenestepper =
        new CompliantImplicitEuler<DIM>(&twodscene, max_iters, criterion);
    twodscene.notifyGlobalIntegrator();
    twodscene.notifyFullIntegrator();
  } else {
    rapidxml::xml_attribute<>* dtnd = nd->first_attribute("maxiters");
    int max_iters = 50;
    if (dtnd) {
      if (!stringutils::extractFromString(std::string(dtnd->value()),
                                          max_iters)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "'maxiters' attribute for integrator. Value must be "
                     "integer. Exiting."
                  << std::endl;
        exit(1);
      }
    }
    dtnd = nd->first_attribute("maxnewton");
    int max_newton = 0;
    if (dtnd) {
      if (!stringutils::extractFromString(std::string(dtnd->value()),
                                          max_newton)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "'max_newton' attribute for integrator. Value must be "
                     "numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }
    dtnd = nd->first_attribute("criterion");
    scalar criterion = 1e-8;
    if (dtnd) {
      if (!stringutils::extractFromString(std::string(dtnd->value()),
                                          criterion)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "'criterion' attribute for integrator. Value must be "
                     "numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }
    dtnd = nd->first_attribute("interhair");
    bool compute_interhair = true;
    if (dtnd) {
      if (!stringutils::extractFromString(std::string(dtnd->value()),
                                          compute_interhair)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "'interhair' attribute for integrator. Value must be "
                     "numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }
    dtnd = nd->first_attribute("preconditioner");
    bool use_preconditioner = true;
    if (dtnd) {
      if (!stringutils::extractFromString(std::string(dtnd->value()),
                                          use_preconditioner)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "'interhair' attribute for integrator. Value must be "
                     "numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }
    *scenestepper = new StrandCompliantManager<DIM>(
        &twodscene, max_newton, max_iters, criterion, compute_interhair,
        use_preconditioner);
    twodscene.notifyGlobalIntegrator();
    twodscene.notifyFullIntegrator();
  }
}

template <int DIM>
void TwoDSceneXMLParser::loadStrandParameters(rapidxml::xml_node<>* node,
                                              TwoDScene<DIM>& scene) {
  /* Available options, if not defined take default below
    <StrandParameters>
      <radius value="" />
      <youngsModulus value="" />
      <shearModulus value="" />
      <density value="" />
      <viscosity value="" />
      <baseRotation value=""/>
      <accumulateWithViscous value="1" />
      <accumulateViscousOnlyForBendingModes value="1" />
      <variableRadiusHair value="1" />
      <straightHairs value="0" />
    </StrandParameters>
  */

  rapidxml::xml_node<>* nd;

  int paramsCount = 0;
  for (nd = node->first_node("StrandParameters"); nd;
       nd = nd->next_sibling("StrandParameters")) {
    // default values:
    scalar radius = 0.0025;
    scalar YoungsModulus = 1e10;
    scalar shearModulus = 3.4e9;
    scalar density = 1.3;
    scalar viscosity = 1e3;
    scalar stretchingMultiplier = 1.0;
    scalar baseRotation = 0.;
    bool accumulateWithViscous = true;
    bool accumulateViscousOnlyForBendingModes = true;
    bool variableRadiusHair = false;
    scalar straightHairs = 1.;
    Vec3 haircolor = Vec3(0, 0, 0);
    rapidxml::xml_node<>* subnd;

    if ((subnd = nd->first_node("haircolor"))) {
      std::string attributer(subnd->first_attribute("r")->value());
      if (!stringutils::extractFromString(attributer, haircolor(0))) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of haircolor attribute for StrandParameters "
                  << paramsCount << ". Value must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
      std::string attributeg(subnd->first_attribute("g")->value());
      if (!stringutils::extractFromString(attributeg, haircolor(1))) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of haircolor attribute for StrandParameters "
                  << paramsCount << ". Value must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
      std::string attributeb(subnd->first_attribute("b")->value());
      if (!stringutils::extractFromString(attributeb, haircolor(2))) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of haircolor attribute for StrandParameters "
                  << paramsCount << ". Value must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    if ((subnd = nd->first_node("radius"))) {
      std::string attribute(subnd->first_attribute("value")->value());
      if (!stringutils::extractFromString(attribute, radius)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of radius attribute for StrandParameters "
                  << paramsCount << ". Value must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    if ((subnd = nd->first_node("youngsModulus"))) {
      std::string attribute(subnd->first_attribute("value")->value());
      if (!stringutils::extractFromString(attribute, YoungsModulus)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of youngsModulus attribute for StrandParameters "
                  << paramsCount << ". Value must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    if ((subnd = nd->first_node("stretchingMultiplier"))) {
      std::string attribute(subnd->first_attribute("value")->value());
      if (!stringutils::extractFromString(attribute, stretchingMultiplier)) {
        std::cerr
            << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse value "
               "of stretchingMultiplier attribute for StrandParameters "
            << paramsCount << ". Value must be numeric. Exiting." << std::endl;
        exit(1);
      }
    }

    if ((subnd = nd->first_node("shearModulus"))) {
      std::string attribute(subnd->first_attribute("value")->value());
      if (!stringutils::extractFromString(attribute, shearModulus)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of shearModulus attribute for StrandParameters "
                  << paramsCount << ". Value must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    if ((subnd = nd->first_node("density"))) {
      std::string attribute(subnd->first_attribute("value")->value());
      if (!stringutils::extractFromString(attribute, density)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of density attribute for StrandParameters "
                  << paramsCount << ". Value must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    if ((subnd = nd->first_node("viscosity"))) {
      std::string attribute(subnd->first_attribute("value")->value());
      if (!stringutils::extractFromString(attribute, viscosity)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of viscosity attribute for StrandParameters "
                  << paramsCount << ". Value must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    if ((subnd = nd->first_node("baseRotation"))) {
      std::string attribute(subnd->first_attribute("value")->value());
      if (!stringutils::extractFromString(attribute, baseRotation)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of baseRotation attribute for StrandParameters "
                  << paramsCount << ". Value must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    if ((subnd = nd->first_node("accumulateWithViscous"))) {
      std::string attribute(subnd->first_attribute("value")->value());
      if (!stringutils::extractFromString(attribute, accumulateWithViscous)) {
        std::cerr
            << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse value "
               "of accumulateWithViscous attribute for StrandParameters "
            << paramsCount << ". Value must be numeric. Exiting." << std::endl;
        exit(1);
      }
    }

    if ((subnd = nd->first_node("accumulateViscousOnlyForBendingModes"))) {
      std::string attribute(subnd->first_attribute("value")->value());
      if (!stringutils::extractFromString(
              attribute, accumulateViscousOnlyForBendingModes)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of accumulateViscousOnlyForBendingModes attribute "
                     "for StrandParameters "
                  << paramsCount << ". Value must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    if ((subnd = nd->first_node("variableRadiusHair"))) {
      std::string attribute(subnd->first_attribute("value")->value());
      if (!stringutils::extractFromString(attribute, variableRadiusHair)) {
        std::cerr
            << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse value "
               "of variableRadiusHair attribute for StrandParameters "
            << paramsCount << ". Value must be numeric. Exiting." << std::endl;
        exit(1);
      }
    }

    if ((subnd = nd->first_node("straightHairs"))) {
      std::string attribute(subnd->first_attribute("value")->value());
      if (!stringutils::extractFromString(attribute, straightHairs)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of straightHairs attribute for StrandParameters "
                  << paramsCount << ". Value must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    scene.insertStrandParameters(new StrandParameters(
        radius, YoungsModulus, shearModulus, stretchingMultiplier, density,
        viscosity, baseRotation, scene.getDt() / scene.getHairSteps(),
        accumulateWithViscous, accumulateViscousOnlyForBendingModes,
        variableRadiusHair, straightHairs, haircolor));
    ++paramsCount;
  }
}

void TwoDSceneXMLParser::loadMaxTime(rapidxml::xml_node<>* node,
                                     scalar& max_t) {
  assert(node != NULL);

  // Attempt to locate the duraiton node
  rapidxml::xml_node<>* nd = node->first_node("duration");
  if (nd == NULL) {
    std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m No duration "
                 "specified. Exiting."
              << std::endl;
    exit(1);
  }

  // Attempt to load the duration value
  rapidxml::xml_attribute<>* timend = nd->first_attribute("time");
  if (timend == NULL) {
    std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m No duration 'time' "
                 "attribute specified. Exiting."
              << std::endl;
    exit(1);
  }

  max_t = std::numeric_limits<scalar>::signaling_NaN();
  if (!stringutils::extractFromString(std::string(timend->value()), max_t)) {
    std::cerr
        << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse 'time' "
           "attribute for duration. Value must be numeric. Exiting."
        << std::endl;
    exit(1);
  }
}

bool TwoDSceneXMLParser::loadViewport(rapidxml::xml_node<>* node,
                                      renderingutils::Viewport& view) {
  assert(node != NULL);

  if (node->first_node("viewport")) {
    rapidxml::xml_attribute<>* cx =
        node->first_node("viewport")->first_attribute("cx");
    if (cx == NULL) {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m No viewport 'cx' "
                   "attribute specified. Exiting."
                << std::endl;
      exit(1);
    }
    if (!stringutils::extractFromString(std::string(cx->value()), view.cx)) {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                   "'cx' attribute for viewport. Value must be scalar. Exiting."
                << std::endl;
      exit(1);
    }
    rapidxml::xml_attribute<>* cy =
        node->first_node("viewport")->first_attribute("cy");
    if (cy == NULL) {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m No viewport 'cy' "
                   "attribute specified. Exiting."
                << std::endl;
      exit(1);
    }
    if (!stringutils::extractFromString(std::string(cy->value()), view.cy)) {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                   "'cy' attribute for viewport. Value must be scalar. Exiting."
                << std::endl;
      exit(1);
    }
    rapidxml::xml_attribute<>* cz =
        node->first_node("viewport")->first_attribute("cz");
    if (cz) {
      if (!stringutils::extractFromString(std::string(cz->value()), view.cz)) {
        std::cerr
            << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse 'cz' "
               "attribute for viewport. Value must be scalar. Exiting."
            << std::endl;
        exit(1);
      }
    }
    rapidxml::xml_attribute<>* size =
        node->first_node("viewport")->first_attribute("size");
    if (size == NULL) {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m No viewport "
                   "'size' attribute specified. Exiting."
                << std::endl;
      exit(1);
    }
    if (!stringutils::extractFromString(std::string(size->value()),
                                        view.size)) {
      std::cerr
          << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse 'size' "
             "attribute for viewport. Value must be scalar. Exiting."
          << std::endl;
      exit(1);
    }

    return true;
  }

  return false;
}

void TwoDSceneXMLParser::loadMaxSimFrequency(rapidxml::xml_node<>* node,
                                             scalar& max_freq) {
  assert(node != NULL);

  // Attempt to locate the duraiton node
  if (node->first_node("maxsimfreq")) {
    // Attempt to load the duration value
    rapidxml::xml_attribute<>* atrbnde =
        node->first_node("maxsimfreq")->first_attribute("max");
    if (atrbnde == NULL) {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m No maxsimfreq "
                   "'max' attribute specified. Exiting."
                << std::endl;
      exit(1);
    }

    if (!stringutils::extractFromString(std::string(atrbnde->value()),
                                        max_freq)) {
      std::cerr
          << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse 'max' "
             "attribute for maxsimfreq. Value must be scalar. Exiting."
          << std::endl;
      exit(1);
    }
  }
}

void TwoDSceneXMLParser::loadBackgroundColor(rapidxml::xml_node<>* node,
                                             renderingutils::Color& color) {
  if (rapidxml::xml_node<>* nd = node->first_node("backgroundcolor")) {
    // Read in the red color channel
    double red = -1.0;
    if (nd->first_attribute("r")) {
      std::string attribute(nd->first_attribute("r")->value());
      if (!stringutils::extractFromString(attribute, red)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of r attribute for backgroundcolor. Value must be "
                     "scalar. Exiting."
                  << std::endl;
        exit(1);
      }
    } else {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                   "value of r attribute for backgroundcolor. Exiting."
                << std::endl;
      exit(1);
    }

    if (red < 0.0 || red > 1.0) {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                   "value of r attribute for backgroundcolor. Invalid color "
                   "specified. Valid range is "
                << 0.0 << "..." << 1.0 << std::endl;
      exit(1);
    }

    // Read in the green color channel
    double green = -1.0;
    if (nd->first_attribute("g")) {
      std::string attribute(nd->first_attribute("g")->value());
      if (!stringutils::extractFromString(attribute, green)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of g attribute for backgroundcolor. Value must be "
                     "scalar. Exiting."
                  << std::endl;
        exit(1);
      }
    } else {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                   "value of g attribute for backgroundcolor. Exiting."
                << std::endl;
      exit(1);
    }

    if (green < 0.0 || green > 1.0) {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                   "value of g attribute for backgroundcolor. Invalid color "
                   "specified. Valid range is "
                << 0.0 << "..." << 1.0 << std::endl;
      exit(1);
    }

    // Read in the blue color channel
    double blue = -1.0;
    if (nd->first_attribute("b")) {
      std::string attribute(nd->first_attribute("b")->value());
      if (!stringutils::extractFromString(attribute, blue)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of b attribute for backgroundcolor. Value must be "
                     "scalar. Exiting."
                  << std::endl;
        exit(1);
      }
    } else {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                   "value of b attribute for backgroundcolor. Exiting."
                << std::endl;
      exit(1);
    }

    if (blue < 0.0 || blue > 1.0) {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                   "value of b attribute for backgroundcolor. Invalid color "
                   "specified. Valid range is "
                << 0.0 << "..." << 1.0 << std::endl;
      exit(1);
    }

    // std::cout << red << "   " << green << "   " << blue << std::endl;

    color.r = red;
    color.g = green;
    color.b = blue;
  }
}

void TwoDSceneXMLParser::loadParticleColors(
    rapidxml::xml_node<>* node,
    std::vector<renderingutils::Color>& particle_colors) {
  int particlecolornum = 0;
  for (rapidxml::xml_node<>* nd = node->first_node("particlecolor"); nd;
       nd = nd->next_sibling("particlecolor")) {
    // Determine which particle this color corresponds to
    int particle = -1;
    if (nd->first_attribute("i")) {
      std::string attribute(nd->first_attribute("i")->value());
      if (!stringutils::extractFromString(attribute, particle)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of i attribute for particlecolor "
                  << particlecolornum << ". Value must be integer. Exiting."
                  << std::endl;
        exit(1);
      }
    } else {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                   "value of i attribute for particlecolor "
                << particlecolornum << ". Exiting." << std::endl;
      exit(1);
    }

    if (particle < 0 || particle >= (int)particle_colors.size()) {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                   "value of i attribute for particlecolor "
                << particlecolornum
                << ". Invalid particle specified. Valid range is " << 0 << "..."
                << particle_colors.size() - 1 << std::endl;
      exit(1);
    }

    // Read in the red color channel
    double red = -1.0;
    if (nd->first_attribute("r")) {
      std::string attribute(nd->first_attribute("r")->value());
      if (!stringutils::extractFromString(attribute, red)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of r attribute for particlecolor "
                  << particlecolornum << ". Value must be scalar. Exiting."
                  << std::endl;
        exit(1);
      }
    } else {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                   "value of r attribute for particlecolor "
                << particlecolornum << ". Exiting." << std::endl;
      exit(1);
    }

    if (red < 0.0 || red > 1.0) {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                   "value of r attribute for particlecolor "
                << particlecolornum
                << ". Invalid color specified. Valid range is " << 0.0 << "..."
                << 1.0 << std::endl;
      exit(1);
    }

    // Read in the green color channel
    double green = -1.0;
    if (nd->first_attribute("g")) {
      std::string attribute(nd->first_attribute("g")->value());
      if (!stringutils::extractFromString(attribute, green)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of g attribute for particlecolor "
                  << particlecolornum << ". Value must be scalar. Exiting."
                  << std::endl;
        exit(1);
      }
    } else {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                   "value of g attribute for particlecolor "
                << particlecolornum << ". Exiting." << std::endl;
      exit(1);
    }

    if (green < 0.0 || green > 1.0) {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                   "value of g attribute for particlecolor "
                << particlecolornum
                << ". Invalid color specified. Valid range is " << 0.0 << "..."
                << 1.0 << std::endl;
      exit(1);
    }

    // Read in the blue color channel
    double blue = -1.0;
    if (nd->first_attribute("b")) {
      std::string attribute(nd->first_attribute("b")->value());
      if (!stringutils::extractFromString(attribute, blue)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of b attribute for particlecolor "
                  << particlecolornum << ". Value must be scalar. Exiting."
                  << std::endl;
        exit(1);
      }
    } else {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                   "value of b attribute for particlecolor "
                << particlecolornum << ". Exiting." << std::endl;
      exit(1);
    }

    if (blue < 0.0 || blue > 1.0) {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                   "value of b attribute for particlecolor "
                << particlecolornum
                << ". Invalid color specified. Valid range is " << 0.0 << "..."
                << 1.0 << std::endl;
      exit(1);
    }

    particle_colors[particle] = renderingutils::Color(red, green, blue);

    ++particlecolornum;
  }
}

void TwoDSceneXMLParser::loadEdgeColors(
    rapidxml::xml_node<>* node,
    std::vector<renderingutils::Color>& edge_colors) {
  int edgecolornum = 0;
  for (rapidxml::xml_node<>* nd = node->first_node("edgecolor"); nd;
       nd = nd->next_sibling("edgecolor")) {
    // Determine which particle this color corresponds to
    int edge = -1;
    if (nd->first_attribute("i")) {
      std::string attribute(nd->first_attribute("i")->value());
      if (!stringutils::extractFromString(attribute, edge)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of i attribute for edgecolor "
                  << edgecolornum << ". Value must be integer. Exiting."
                  << std::endl;
        exit(1);
      }
    } else {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                   "value of i attribute for edgecolor "
                << edgecolornum << ". Exiting." << std::endl;
      exit(1);
    }

    if (edge < 0 || edge > (int)edge_colors.size()) {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                   "value of i attribute for edgecolor "
                << edgecolornum << ". Invalid edge specified. Valid range is "
                << 0 << "..." << edge_colors.size() - 1 << std::endl;
      exit(1);
    }

    // Read in the red color channel
    double red = -1.0;
    if (nd->first_attribute("r")) {
      std::string attribute(nd->first_attribute("r")->value());
      if (!stringutils::extractFromString(attribute, red)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of r attribute for edgecolor "
                  << edgecolornum << ". Value must be scalar. Exiting."
                  << std::endl;
        exit(1);
      }
    } else {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                   "value of r attribute for edgecolor "
                << edgecolornum << ". Exiting." << std::endl;
      exit(1);
    }

    if (red < 0.0 || red > 1.0) {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                   "value of r attribute for edgecolor "
                << edgecolornum << ". Invalid color specified. Valid range is "
                << 0.0 << "..." << 1.0 << std::endl;
      exit(1);
    }

    // Read in the green color channel
    double green = -1.0;
    if (nd->first_attribute("g")) {
      std::string attribute(nd->first_attribute("g")->value());
      if (!stringutils::extractFromString(attribute, green)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of g attribute for edgecolor "
                  << edgecolornum << ". Value must be scalar. Exiting."
                  << std::endl;
        exit(1);
      }
    } else {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                   "value of g attribute for edgecolor "
                << edgecolornum << ". Exiting." << std::endl;
      exit(1);
    }

    if (green < 0.0 || green > 1.0) {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                   "value of g attribute for edgecolor "
                << edgecolornum << ". Invalid color specified. Valid range is "
                << 0.0 << "..." << 1.0 << std::endl;
      exit(1);
    }

    // Read in the blue color channel
    double blue = -1.0;
    if (nd->first_attribute("b")) {
      std::string attribute(nd->first_attribute("b")->value());
      if (!stringutils::extractFromString(attribute, blue)) {
        std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                     "value of b attribute for edgecolor "
                  << edgecolornum << ". Value must be scalar. Exiting."
                  << std::endl;
        exit(1);
      }
    } else {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                   "value of b attribute for edgecolor "
                << edgecolornum << ". Exiting." << std::endl;
      exit(1);
    }

    if (blue < 0.0 || blue > 1.0) {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse "
                   "value of b attribute for edgecolor "
                << edgecolornum << ". Invalid color specified. Valid range is "
                << 0.0 << "..." << 1.0 << std::endl;
      exit(1);
    }

    // std::cout << "edge: " << edge << " r: " << red << " g: " << green << " b:
    // " << blue << std::endl;

    edge_colors[edge] = renderingutils::Color(red, green, blue);

    ++edgecolornum;
  }
}

void TwoDSceneXMLParser::loadSceneDescriptionString(
    rapidxml::xml_node<>* node, std::string& description_string) {
  assert(node != NULL);

  description_string = "No description specified.";

  // Attempt to locate the integrator node
  rapidxml::xml_node<>* nd = node->first_node("description");
  if (nd != NULL) {
    rapidxml::xml_attribute<>* typend = nd->first_attribute("text");
    if (typend == NULL) {
      std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m No text attribute "
                   "specified for description. Exiting."
                << std::endl;
      exit(1);
    }
    description_string = typend->value();
  }
}

bool TwoDSceneXMLParser::loadCamera(rapidxml::xml_node<>* node,
                                    Camera& camera) {
  rapidxml::xml_node<>* nd = node->first_node("camera");
  if (nd) {
    Eigen::Quaterniond rotation(1, 0, 0, 0);
    rapidxml::xml_node<>* nd_rot = nd->first_node("rotation");

    if (nd_rot) {
      if (nd_rot->first_attribute("x")) {
        std::string attribute(nd_rot->first_attribute("x")->value());
        if (!stringutils::extractFromString(attribute, rotation.x())) {
          std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to "
                       "parse value of x attribute for rotation. Value must be "
                       "scalar. Exiting."
                    << std::endl;
          exit(1);
        }
      }
      if (nd_rot->first_attribute("y")) {
        std::string attribute(nd_rot->first_attribute("y")->value());
        if (!stringutils::extractFromString(attribute, rotation.y())) {
          std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to "
                       "parse value of y attribute for rotation. Value must be "
                       "scalar. Exiting."
                    << std::endl;
          exit(1);
        }
      }
      if (nd_rot->first_attribute("z")) {
        std::string attribute(nd_rot->first_attribute("z")->value());
        if (!stringutils::extractFromString(attribute, rotation.z())) {
          std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to "
                       "parse value of z attribute for rotation. Value must be "
                       "scalar. Exiting."
                    << std::endl;
          exit(1);
        }
      }
      if (nd_rot->first_attribute("w")) {
        std::string attribute(nd_rot->first_attribute("w")->value());
        if (!stringutils::extractFromString(attribute, rotation.w())) {
          std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to "
                       "parse value of w attribute for rotation. Value must be "
                       "scalar. Exiting."
                    << std::endl;
          exit(1);
        }
      }
    }

    Eigen::Vector3d center(0, 0, 0);
    rapidxml::xml_node<>* nd_center = nd->first_node("center");

    if (nd_center) {
      if (nd_center->first_attribute("x")) {
        std::string attribute(nd_center->first_attribute("x")->value());
        if (!stringutils::extractFromString(attribute, center.x())) {
          std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to "
                       "parse value of x attribute for center. Value must be "
                       "scalar. Exiting."
                    << std::endl;
          exit(1);
        }
      }
      if (nd_center->first_attribute("y")) {
        std::string attribute(nd_center->first_attribute("y")->value());
        if (!stringutils::extractFromString(attribute, center.y())) {
          std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to "
                       "parse value of y attribute for center. Value must be "
                       "scalar. Exiting."
                    << std::endl;
          exit(1);
        }
      }
      if (nd_center->first_attribute("z")) {
        std::string attribute(nd_center->first_attribute("z")->value());
        if (!stringutils::extractFromString(attribute, center.z())) {
          std::cerr << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to "
                       "parse value of z attribute for center. Value must be "
                       "scalar. Exiting."
                    << std::endl;
          exit(1);
        }
      }
    }

    scalar dist = 0.0;
    scalar radius = 100.0;
    scalar fov = 40.0;

    if (nd->first_attribute("dist")) {
      std::string attribute(nd->first_attribute("dist")->value());
      if (!stringutils::extractFromString(attribute, dist)) {
        std::cerr
            << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse value "
               "of dist attribute for camera. Value must be scalar. Exiting."
            << std::endl;
        exit(1);
      }
    }

    if (nd->first_attribute("radius")) {
      std::string attribute(nd->first_attribute("radius")->value());
      if (!stringutils::extractFromString(attribute, radius)) {
        std::cerr
            << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse value "
               "of radius attribute for camera. Value must be scalar. Exiting."
            << std::endl;
        exit(1);
      }
    }

    if (nd->first_attribute("fov")) {
      std::string attribute(nd->first_attribute("fov")->value());
      if (!stringutils::extractFromString(attribute, fov)) {
        std::cerr
            << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse value "
               "of fov attribute for camera. Value must be scalar. Exiting."
            << std::endl;
        exit(1);
      }
    }

    camera.rotation_ = rotation;
    camera.center_ = center;
    camera.dist_ = dist;
    camera.radius_ = radius;
    camera.fov_ = fov;

    return true;
  }

  return false;
}
