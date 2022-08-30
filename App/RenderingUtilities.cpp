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

#include "RenderingUtilities.h"

namespace libwethair {
namespace renderingutils {

bool checkGLErrors() {
  /*
  GLenum errCode;
  const GLubyte *errString;

  if ((errCode = glGetError()) != GL_NO_ERROR)
  {
    errString = gluErrorString(errCode);
    std::cout << outputmod::startred << "OpenGL Error:" << outputmod::endred <<
  std::flush; fprintf(stderr, " %s\n", errString); return false;
  }*/
  return true;
}

Color::Color(const Vector3s& c) : r(c(0)), g(c(1)), b(c(2)) {}

Color::Color() : r(0.0), g(0.0), b(0.0) {}

Color::Color(double r, double g, double b) : r(r), g(g), b(b) {
  assert(r >= 0.0);
  assert(r <= 1.0);
  assert(g >= 0.0);
  assert(g <= 1.0);
  assert(b >= 0.0);
  assert(b <= 1.0);
}

Vector3s Color::toVector() const { return Vector3s(r, g, b); }

}  // namespace renderingutils
}  // namespace libwethair
