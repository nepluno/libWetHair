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

#include "liangbarsky.h"
#include "MathDefs.h"

namespace liangbarsky
{
  int clip_line(const Vector4s& c, Vector2s& q1, Vector2s& q2, scalar& t0, scalar& t1)
  {
    t0 = 0.0;
    t1 = 1.0;
    double xdelta = q2(0)-q1(0);
    double ydelta = q2(1)-q1(1);
    double p=1.0,q=0.0,r;
    
    for(int edge=0; edge<4; edge++) {   // Traverse through left, right, bottom, top edges.
      if (edge==0) {  p = -xdelta;    q = -(c(0)-q1(0));  }
      if (edge==1) {  p = xdelta;     q =  (c(2)-q1(0)); }
      if (edge==2) {  p = -ydelta;    q = -(c(1)-q1(1));}
      if (edge==3) {  p = ydelta;     q =  (c(3)-q1(1));   }
      r = q/p;
      if(p==0 && q<0) return false;   // Don't draw line at all. (parallel line outside)
      
      if(p<0) {
        if(r>t1) return false;         // Don't draw line at all.
        else if(r>t0) t0=r;            // Line is clipped!
      } else if(p>0) {
        if(r<t0) return false;      // Don't draw line at all.
        else if(r<t1) t1=r;         // Line is clipped!
      }
    }
    
    Vector2s tq0 = q1;
    
    q1 = tq0 + t0 * Vector2s(xdelta, ydelta);
    q2 = tq0 + t1 * Vector2s(xdelta, ydelta);

    return true;        // (clipped) line is drawn
  }
  
  int clip_line(const Vector6s& c, Vector3s& q1, Vector3s& q2, scalar& t0, scalar& t1)
  {
    t0 = 0.0;
    t1 = 1.0;
    double xdelta = q2(0)-q1(0);
    double ydelta = q2(1)-q1(1);
    double zdelta = q2(2)-q1(2);
    double p=1.0,q=0.0,r;
    
    for(int edge=0; edge<6; edge++) {   // Traverse through left, right, bottom, top edges.
      if (edge==0) {  p = -xdelta;    q = -(c(0)-q1(0));  }
      if (edge==1) {  p = xdelta;     q =  (c(3)-q1(0)); }
      if (edge==2) {  p = -ydelta;    q = -(c(1)-q1(1));}
      if (edge==3) {  p = ydelta;     q =  (c(4)-q1(1));   }
      if (edge==4) {  p = -zdelta;    q = -(c(2)-q1(2));}
      if (edge==5) {  p = zdelta;     q =  (c(5)-q1(2));   }
      
      r = q/p;
      if(p==0 && q<0) return false;   // Don't draw line at all. (parallel line outside)
      
      if(p<0) {
        if(r>t1) return false;         // Don't draw line at all.
        else if(r>t0) t0=r;            // Line is clipped!
      } else if(p>0) {
        if(r<t0) return false;      // Don't draw line at all.
        else if(r<t1) t1=r;         // Line is clipped!
      }
    }
    
    Vector3s tq0 = q1;
    
    q1 = tq0 + t0 * Vector3s(xdelta, ydelta, zdelta);
    q2 = tq0 + t1 * Vector3s(xdelta, ydelta, zdelta);
    
    return true;        // (clipped) line is drawn
  }
}

