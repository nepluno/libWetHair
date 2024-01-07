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

#ifndef LIBWETHAIR_APP_EXECUTABLE_SIMULATION_H_
#define LIBWETHAIR_APP_EXECUTABLE_SIMULATION_H_

#include <Eigen/Core>
#include <string>

#include <libWetHair/MathDefs.h>

#include "Camera.h"
#include "RenderingUtilities.h"

template <int DIM>
class TwoDSceneRenderer;

class ExecutableSimulation {
 public:
  virtual ~ExecutableSimulation() {}

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

  virtual void computeCameraCenter(
      libwethair::renderingutils::Viewport& view) = 0;

  /////////////////////////////////////////////////////////////////////////////
  // Serialization Functions

  virtual void serializeScene(std::ostream& outputstream) = 0;
  virtual void serializeSceneReadable(std::vector<std::ostringstream>& osfluid,
                                      std::ostream& oshair,
                                      std::ostream& osflow,
                                      std::ostream& os_boundary_single,
                                      std::ostream& os_boundary_double,
                                      std::ostream& os_pe, std::ostream& os_poe,
                                      std::ostream& os_ppp) = 0;
  virtual bool deSerializeSceneReadable(
      const std::vector<std::string>& filename_fluids,
      const std::string& filename_hairs, const std::string& filename_flows,
      const std::string& filename_bd_single,
      const std::string& filename_bd_double) = 0;
  virtual void loadSerializedScene(std::istream& inputstream) = 0;

  /////////////////////////////////////////////////////////////////////////////
  // Status Functions

  virtual std::string getSolverName() = 0;

  virtual void centerCamera(bool b_reshape = true) = 0;
  virtual void keyboard(int key) = 0;
  virtual void reshape(int w, int h) = 0;

  virtual void mouse(int button, int state, int res_scale) = 0;

  virtual void translateView(double dx, double dy) = 0;

  virtual void zoomView(double dx, double dy) = 0;

  virtual void motion(int x, int y) = 0;

  virtual const libwethair::scalar& getDt() const = 0;
  virtual const libwethair::scalar& getCurrentTime() const = 0;
  virtual libwethair::scalar getCurrentFrame() const = 0;

  virtual int getWindowWidth() const = 0;
  virtual int getWindowHeight() const = 0;

  virtual void setWindowWidth(int w) = 0;
  virtual void setWindowHeight(int h) = 0;

  virtual void setCenterX(double x) = 0;
  virtual void setCenterY(double y) = 0;
  virtual void setScaleFactor(double scale) = 0;

  virtual void setCamera(const Camera& cam) = 0;
  virtual void setView(const libwethair::renderingutils::Viewport& view) = 0;

  virtual const std::vector<libwethair::scalar>& getStepperTimingStatistics()
      const = 0;

  virtual const std::vector<libwethair::scalar>& getTimingStatistics()
      const = 0;
};

#endif  // LIBWETHAIR_APP_EXECUTABLE_SIMULATION_H_
