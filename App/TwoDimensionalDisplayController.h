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

#ifndef LIBWETHAIR_APP_TWO_DIMENSIONAL_DISPLAY_CONTROLLER_H_
#define LIBWETHAIR_APP_TWO_DIMENSIONAL_DISPLAY_CONTROLLER_H_

#include <cassert>
#include <cmath>
#include <iostream>
#include <stack>

#include "Camera.h"
#include "RenderingUtilities.h"

extern bool g_rendering_enabled;

template <int DIM>
class TwoDSceneRenderer;

template <int DIM>
class TwoDimensionalDisplayController {
 public:
  TwoDimensionalDisplayController(GLFWwindow* window, int width, int height);

  void reshape(int w, int h);

  void keyboard(int key);

  void mouse(int button, int state, int res_scale);

  void translateView(double dx, double dy);

  void translateViewZ(double dx, double dz);

  void zoomView(double dx, double dy);

  void motion(int x, int y);

  int getWindowWidth() const;
  int getWindowHeight() const;

  int getWorldWidth() const;
  int getWorldHeight() const;

  double getCenterX() const;
  double getCenterY() const;
  void setCenterX(double x);
  void setCenterY(double y);
  void setCenterZ(double z);

  void setScaleFactor(double scale);
  double getScaleFactor() const;
  double getRatio() const;

  void applyProjection() const;

  void applyOrtho(double ratio, double add_scale = 1.0) const;

  void applyOrthoZ(double ratio, double add_scale = 1.0) const;

  void setRender(TwoDSceneRenderer<DIM>* render);

  void setWindowWidth(int w);
  void setWindowHeight(int h);
  void initCamera(const libwethair::renderingutils::Viewport& view);
  void getEye(libwethair::Vector3s& eye) const;

  const Camera& getCamera() const;

  Camera& getCamera();

  void initDefaultCamera();

 private:
  // Width of the window in pixels
  int m_window_width;
  // Height of the window in pixels
  int m_window_height;
  // Factor to 'zoom' in or out by
  double m_scale_factor;
  // Center of the display, x coord
  double m_center_x;
  // Center of the display, y coord
  double m_center_y;
  // Center of the display, y coord
  double m_center_z;

  // True if the user is dragging the display left
  bool m_left_drag;
  // True if the user is dragging the display right
  bool m_right_drag;
  // Last position of the cursor in a drag, x coord
  double m_last_x;
  // Last position of the cursor in a drag, y coord
  double m_last_y;

  bool m_right_part_click;

  TwoDSceneRenderer<DIM>* m_render;

  Camera m_camera;

  std::stack<Camera> m_cam_stack;

  GLFWwindow* m_window;
};

template <int DIM>
TwoDimensionalDisplayController<DIM>* GetDisplayController();

#endif  // LIBWETHAIR_APP_TWO_DIMENSIONAL_DISPLAY_CONTROLLER_H_
