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

#ifndef LIBWETHAIR_APP_RENDERING_UTILITIES_H_
#define LIBWETHAIR_APP_RENDERING_UTILITIES_H_

#ifdef WIN32
#define NOMINMAX
#include <Windows.h>
#endif
#include <GL/glew.h>
#include <GLFW/glfw3.h>

#ifdef __APPLE__
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include <cstdio>
#include <iostream>
#include <list>

#include <libWetHair/JET.h>
#include <libWetHair/MathDefs.h>
#include <libWetHair/MathUtilities.h>
#include "StringUtilities.h"

namespace libwethair {
namespace renderingutils {
class Color {
 public:
  Color();

  Color(double r, double g, double b);

  Color(const Vector3s&);

  Vector3s toVector() const;

  double r;
  double g;
  double b;
};
}  // namespace renderingutils
}  // namespace libwethair

extern int getDCWindowWidth();
extern int getDCWindowHeight();
extern libwethair::renderingutils::Color& getDCBackgroundColor();

namespace libwethair {
namespace renderingutils {
// False => error
bool checkGLErrors();

// TODO: Move these functions to scene renderer?
inline void setOrthographicProjection() {
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();

  gluOrtho2D(0, getDCWindowWidth(), 0, getDCWindowHeight());

  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();

  assert(renderingutils::checkGLErrors());
}

inline Vector3s interpolateColor(const scalar& x, const scalar& xmin,
                                 const scalar& xmax) {
  scalar dm = (xmax - xmin);

  scalar a;
  if (dm == 0.0)
    a = x;
  else
    a = (x - xmin) / dm * (scalar)(jetmapping_size - 1);

  int isel = std::max(std::min((int)a, jetmapping_size - 1), 0);
  int inext = (isel + 1) % (jetmapping_size);
  scalar fraca = std::max(std::min(a - (scalar)isel, 1.0), 0.0);

  return mathutils::lerp(jetmapping_real[isel], jetmapping_real[inext], fraca);
}

struct Viewport {
 public:
  Viewport()
      : cx(0.0), cy(0.0), cz(0.0), rx(1.0), ry(1.0), rz(1.0), size(1.2) {}
  double cx;
  double cy;
  double cz;
  double rx;
  double ry;
  double rz;
  double size;
};

}  // namespace renderingutils
}  // namespace libwethair

#endif  // LIBWETHAIR_APP_RENDERING_UTILITIES_H_
