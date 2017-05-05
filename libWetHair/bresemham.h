//
// This file is part of the libWetHair open source project
//
// The code is licensed solely for academic and non-commercial use under the
// terms of the Clear BSD License. The terms of the Clear BSD License are
// provided below. Other licenses may be obtained by contacting the faculty
// of the Columbia Computer Graphics Group or a Columbia University licensing officer.
//
// The Clear BSD License
//
// Copyright 2017 Yun (Raymond) Fei, Henrique Teles Maia, Christopher Batty,
// Changxi Zheng, and Eitan Grinspun
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted (subject to the limitations in the disclaimer
// below) provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice, this
//  list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright notice,
//  this list of conditions and the following disclaimer in the documentation
//  and/or other materials provided with the distribution.
//
// * Neither the name of the copyright holder nor the names of its contributors may be used
//  to endorse or promote products derived from this software without specific
//  prior written permission.
//
// NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE GRANTED BY THIS
// LICENSE. THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
// GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
// OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.


#ifndef BRESEMHAM_H__
#define BRESEMHAM_H__

#include "MathDefs.h"

namespace bresemham {
  template<typename Callable>
  void Bresenham3D(int x1, int y1, int z1, const int x2, const int y2, const int z2, Callable func){
    
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
  
  
  template<typename Callable>
  void Bresenham2D(int x1, int y1, const int x2, const int y2, Callable func){
    
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

};

#endif
