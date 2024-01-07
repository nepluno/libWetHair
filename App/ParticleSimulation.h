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

#ifndef LIBWETHAIR_APP_PARTICLE_SIMULATION_H_
#define LIBWETHAIR_APP_PARTICLE_SIMULATION_H_

#include <libWetHair/SceneStepper.h>
#include <libWetHair/TwoDScene.h>
#include <libWetHair/WetHairCore.h>

#include "ExecutableSimulation.h"
#include "TwoDSceneRenderer.h"
#include "TwoDSceneSerializer.h"
#include "TwoDimensionalDisplayController.h"

extern bool g_bPrintSummary;
extern bool g_rendering_enabled;

// TODO: Move code out of header!
template <int DIM>
class ParticleSimulation : public ExecutableSimulation {
 public:
  ParticleSimulation(GLFWwindow* window, libwethair::TwoDScene<DIM>* scene,
                     libwethair::SceneStepper<DIM>* scene_stepper,
                     TwoDSceneRenderer<DIM>* scene_renderer,
                     const std::vector<libwethair::Script>& scripts);
  ~ParticleSimulation();

  virtual void stepSystem();
  /////////////////////////////////////////////////////////////////////////////
  // Rendering Functions

  virtual void initializeOpenGLRenderer();
  virtual void renderSceneOpenGL();
  virtual void renderSceneDifferencesOpenGL();
  virtual void updateOpenGLRendererState();
  virtual void computeCameraCenter(libwethair::renderingutils::Viewport& view);
  /////////////////////////////////////////////////////////////////////////////
  // Serialization Functions

  virtual void serializeScene(std::ostream& outputstream);
  virtual void serializeSceneReadable(std::vector<std::ostringstream>& osfluid,
                                      std::ostream& oshair,
                                      std::ostream& osflow,
                                      std::ostream& os_boundary_single,
                                      std::ostream& os_boundary_double,
                                      std::ostream& os_pe, std::ostream& os_poe,
                                      std::ostream& os_ppp);
  virtual bool deSerializeSceneReadable(
      const std::vector<std::string>& filename_fluids,
      const std::string& filename_hairs, const std::string& filename_flows,
      const std::string& filename_bd_single,
      const std::string& filename_bd_double);
  virtual void loadSerializedScene(std::istream& inputstream);
  /////////////////////////////////////////////////////////////////////////////
  // Status Functions

  virtual void centerCamera(bool b_reshape = true);
  virtual void keyboard(int key);
  virtual void reshape(int w, int h);

  virtual void mouse(int button, int state, int res_scale);

  virtual void translateView(double dx, double dy);

  virtual void zoomView(double dx, double dy);

  virtual void motion(int x, int y);

  virtual int getWindowWidth() const;
  virtual int getWindowHeight() const;

  virtual void setWindowWidth(int w);
  virtual void setWindowHeight(int h);

  virtual const libwethair::scalar& getDt() const;
  virtual const libwethair::scalar& getCurrentTime() const;
  virtual libwethair::scalar getCurrentFrame() const;

  virtual TwoDimensionalDisplayController<DIM>* getDC();

  virtual void setCenterX(double x);
  virtual void setCenterY(double y);
  virtual void setCenterZ(double z);
  virtual void setScaleFactor(double scale);

  virtual void setCamera(const Camera& cam);
  virtual void setView(const libwethair::renderingutils::Viewport& view);

  virtual bool stepScript(const libwethair::scalar& dt);

  virtual std::string getSolverName();

  virtual const std::vector<libwethair::scalar>& getTimingStatistics() const;
  virtual const std::vector<libwethair::scalar>& getStepperTimingStatistics()
      const;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
 private:
  libwethair::WetHairCore<DIM>* m_wet_hair_core;

  libwethair::TwoDScene<DIM>* m_scene;

  libwethair::SceneStepper<DIM>* m_scene_stepper;

  TwoDSceneRenderer<DIM>* m_scene_renderer;

  TwoDimensionalDisplayController<DIM> m_display_controller;

  TwoDSceneSerializer<DIM> m_scene_serializer;

  std::vector<libwethair::Script> m_scripts;
};

#endif  // LIBWETHAIR_APP_PARTICLE_SIMULATION_H_
