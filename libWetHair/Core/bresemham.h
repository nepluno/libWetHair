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

#ifndef LIBWETHAIR_CORE_BRESEMHAM_H_
#define LIBWETHAIR_CORE_BRESEMHAM_H_

#include "MathDefs.h"

namespace bresemham {
template <typename Callable>
void Bresenham3D(int x1, int y1, int z1, const int x2, const int y2,
                 const int z2, Callable func) {
  int i, dx, dy, dz, l, m, n, x_inc, y_inc, z_inc, err_1, err_2, dx2, dy2, dz2;
  int point[3];

  point[0] = x1;
  point[1] = y1;
  point[2] = z1;
  dx = x2 - x1;
  dy = y2 - y1;
  dz = z2 - z1;
  x_inc = (dx < 0) ? -1 : 1;
  l = abs(dx);
  y_inc = (dy < 0) ? -1 : 1;
  m = abs(dy);
  z_inc = (dz < 0) ? -1 : 1;
  n = abs(dz);
  dx2 = l << 1;
  dy2 = m << 1;
  dz2 = n << 1;

  if ((l >= m) && (l >= n)) {
    err_1 = dy2 - l;
    err_2 = dz2 - l;
    for (i = 0; i < l; i++) {
      func(point[0], point[1], point[2]);
      if (err_1 > 0) {
        point[1] += y_inc;
        err_1 -= dx2;
      }
      if (err_2 > 0) {
        point[2] += z_inc;
        err_2 -= dx2;
      }
      err_1 += dy2;
      err_2 += dz2;
      point[0] += x_inc;
    }
  } else if ((m >= l) && (m >= n)) {
    err_1 = dx2 - m;
    err_2 = dz2 - m;
    for (i = 0; i < m; i++) {
      func(point[0], point[1], point[2]);
      if (err_1 > 0) {
        point[0] += x_inc;
        err_1 -= dy2;
      }
      if (err_2 > 0) {
        point[2] += z_inc;
        err_2 -= dy2;
      }
      err_1 += dx2;
      err_2 += dz2;
      point[1] += y_inc;
    }
  } else {
    err_1 = dy2 - n;
    err_2 = dx2 - n;
    for (i = 0; i < n; i++) {
      func(point[0], point[1], point[2]);
      if (err_1 > 0) {
        point[1] += y_inc;
        err_1 -= dz2;
      }
      if (err_2 > 0) {
        point[0] += x_inc;
        err_2 -= dz2;
      }
      err_1 += dy2;
      err_2 += dx2;
      point[2] += z_inc;
    }
  }
  func(point[0], point[1], point[2]);
}

template <typename Callable>
void Bresenham2D(int x1, int y1, const int x2, const int y2, Callable func) {
  int i, dx, dy, l, m, x_inc, y_inc, err, dx2, dy2;
  int point[2];

  point[0] = x1;
  point[1] = y1;
  dx = x2 - x1;
  dy = y2 - y1;
  x_inc = (dx < 0) ? -1 : 1;
  l = abs(dx);
  y_inc = (dy < 0) ? -1 : 1;
  m = abs(dy);
  dx2 = l << 1;
  dy2 = m << 1;

  if (l >= m) {
    err = dy2 - l;
    for (i = 0; i < l; i++) {
      func(point[0], point[1]);
      if (err > 0) {
        point[1] += y_inc;
        err -= dx2;
      }
      err += dy2;
      point[0] += x_inc;
    }
  } else {
    err = dx2 - m;
    for (i = 0; i < m; i++) {
      func(point[0], point[1]);
      if (err > 0) {
        point[0] += x_inc;
        err -= dy2;
      }
      err += dx2;
      point[1] += y_inc;
    }
  }
  func(point[0], point[1]);
}

};  // namespace bresemham

#endif  // LIBWETHAIR_CORE_BRESEMHAM_H_
