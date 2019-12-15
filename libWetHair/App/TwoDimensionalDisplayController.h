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


#ifndef __TWO_DIMENSIONAL_DISPLAY_CONTROLLER_H__
#define __TWO_DIMENSIONAL_DISPLAY_CONTROLLER_H__

#include <cmath>
#include <iostream>
#include <cassert>
#include <stack>

#include "RenderingUtilities.h"
#include "Camera.h"

extern bool g_rendering_enabled;

template<int DIM>
class TwoDSceneRenderer;

template<int DIM>
class TwoDimensionalDisplayController
{
public:
  TwoDimensionalDisplayController( int width, int height );
  
  void reshape( int w, int h );
  
  void keyboard( unsigned char key, int x, int y );
  
  void special( int key, int x, int y );

  void mouse( int button, int state, int x, int y );
  
  void translateView( double dx, double dy );
  
  void translateViewZ( double dx, double dz );
  
  void zoomView( double dx, double dy );
  
  void motion( int x, int y );
  
  int getWindowWidth() const;  
  int getWindowHeight() const;

  int getWorldWidth() const;
  int getWorldHeight() const;
  
  double getCenterX() const;
  double getCenterY() const;
  void setCenterX( double x );
  void setCenterY( double y );
  void setCenterZ( double z );
  
  void setScaleFactor( double scale );
  double getScaleFactor() const;
  double getRatio() const;
  
  void applyProjection() const;
  
  void applyOrtho( double ratio, double add_scale = 1.0 ) const;
  
  void applyOrthoZ( double ratio, double add_scale = 1.0 ) const;
  
  void setRender(TwoDSceneRenderer<DIM>* render);
  
  void setWindowWidth(int w);
  void setWindowHeight(int h);
  void initCamera(const renderingutils::Viewport& view);
  void getEye(Vector3s& eye) const;
  
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
  
  unsigned m_modifiers;
  
  TwoDSceneRenderer<DIM>* m_render;
  
  Camera m_camera;
  
  std::stack<Camera> m_cam_stack;
};

template<int DIM>
TwoDimensionalDisplayController<DIM>* GetDisplayController();

#endif
