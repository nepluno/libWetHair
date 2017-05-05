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

#include "openglutils.h"

#ifdef __APPLE__
#include <GLUT/glut.h> // why does Apple have to put glut.h here...
#else
#ifdef WIN32
#define NOMINMAX
#include <Windows.h>
#endif
#include <GL/glut.h> // ...when everyone else puts it here?
#endif

#include "MathDefs.h"
#include "MathUtilities.h"
#include <cfloat>

#ifdef WIN32
#define _USE_MATH_DEFINES
#endif
#include <cmath>

void draw_circle2d(const Vector2s& centre, scalar rad, int segs)
{
  glBegin(GL_POLYGON);
  for(int i=0;i<segs;i++){
    scalar cosine=rad*cos(i*2*PI/(scalar)(segs));
    scalar sine=rad* sin(i*2*PI/(scalar)(segs));
    Vector2s tmp = (Vector2s(cosine,sine) + centre);
    glVertex2dv(tmp.data());
  }
  glEnd();
}

void draw_grid2d(const Vector2s& origin, scalar dx, int nx, int ny) {
  scalar width = nx*dx;
  scalar height = ny*dx;
  
  glBegin(GL_LINES);
  for(int i = 0; i <= nx; i++) {
    Vector2s a(i*dx, 0);
    Vector2s b(i*dx, height);
    Vector2s oa = origin + a;
    Vector2s ob = origin + b;
    glVertex2dv(oa.data());
    glVertex2dv(ob.data());
  }
  for(int j = 0; j <= ny; ++j) {
    Vector2s a(0,j*dx);
    Vector2s b(width,j*dx);
    Vector2s oa = origin + a;
    Vector2s ob = origin + b;
    glVertex2dv(oa.data());
    glVertex2dv(ob.data());
  }
  glEnd();
}

void draw_grid3d(const Vector3s& origin, scalar dx, int nx, int ny, int nz) {
  scalar width = nx*dx;
  scalar height = ny*dx;
  scalar depth = nz*dx;
  
  glBegin(GL_LINES);
  for(int k = 0; k <= nz; ++k) for(int i = 0; i <= nx; ++i) {
    Vector3s a(i*dx, 0, k*dx);
    Vector3s b(i*dx, height, k*dx);
    Vector3s oa = origin + a;
    Vector3s ob = origin + b;
    glVertex3dv(oa.data());
    glVertex3dv(ob.data());
  }
  for(int k = 0; k <= nz; ++k) for(int j = 0; j <= ny; ++j) {
    Vector3s a(0,j*dx,k*dx);
    Vector3s b(width,j*dx,k*dx);
    Vector3s oa = origin + a;
    Vector3s ob = origin + b;
    glVertex3dv(oa.data());
    glVertex3dv(ob.data());
  }
  for(int j = 0; j <= ny; ++j) for(int i = 0; i <= nx; ++i) {
    Vector3s a(i*dx,j*dx,0);
    Vector3s b(i*dx,j*dx,depth);
    Vector3s oa = origin + a;
    Vector3s ob = origin + b;
    glVertex3dv(oa.data());
    glVertex3dv(ob.data());
  }
  glEnd();
}

void draw_box2d(const Vector2s& origin, scalar width, scalar height) {
  Vector2s o1 = origin + Vector2s(0, height);
  Vector2s o2 = origin + Vector2s(width, height);
  Vector2s o3 = origin + Vector2s(width, 0);
  
  glBegin(GL_POLYGON);
  glVertex2dv(origin.data());
  glVertex2dv(o1.data());
  glVertex2dv(o2.data());
  glVertex2dv(o3.data());
  glEnd();
}

void draw_segmentset2d(const std::vector<Vector2s>& vertices, const std::vector<Vector2i>& edges) {
  glBegin(GL_LINES);
  for(unsigned int i = 0; i < edges.size(); ++i) {
    glVertex2dv(vertices[edges[i][0]].data());
    glVertex2dv(vertices[edges[i][1]].data());
  }
  glEnd();
}

void draw_points2d(const std::vector<Vector2s>& points) {
  glBegin(GL_POINTS);
  for(unsigned int i = 0; i < points.size(); ++i) {
    glVertex2dv(points[i].data());
  }
  glEnd();
}

void draw_polygon2d(const std::vector<Vector2s>& vertices) {
  glBegin(GL_POLYGON);
  for(unsigned int i = 0; i < vertices.size(); ++i)
    glVertex2dv(vertices[i].data());
  glEnd();
}

void draw_polygon2d(const std::vector<Vector2s>& vertices, const std::vector<int>& order) {
  glBegin(GL_POLYGON);
  for(unsigned int i = 0; i < order.size(); ++i)
    glVertex2dv(vertices[order[i]].data());
  glEnd();
  
}
void draw_segment2d(const Vector2s& start, const Vector2s& end) {
  glBegin(GL_LINES);
  glVertex2dv(start.data());
  glVertex2dv(end.data());
  glEnd();
}

void draw_arrow2d(const Vector2s& start, const Vector2s& end, scalar arrow_head_len)
{
  Vector2s direction = end - start;
  
  Vector2s dir_norm = direction;
  
  //TODO Possibly automatically scale arrowhead length based on vector magnitude
  if(dir_norm.norm() < 1e-14)
    return;
  
  dir_norm.normalize();
  Vector2s perp(dir_norm[1],-dir_norm[0]);
  
  Vector2s tip_left = end + arrow_head_len/(scalar)sqrt(2.0)*(-dir_norm + perp);
  Vector2s tip_right = end + arrow_head_len/(scalar)sqrt(2.0)*(-dir_norm - perp);
  
  glBegin(GL_LINES);
  glVertex2dv(start.data());
  glVertex2dv(end.data());
  glVertex2dv(end.data());
  glVertex2dv(tip_left.data());
  glVertex2dv(end.data());
  glVertex2dv(tip_right.data());
  glEnd();
  
}

void draw_arrow3d(const Vector3s& start, const Vector3s& end, const Vector3s& view_dir, scalar arrow_head_len)
{
  Vector3s direction = end - start;
  
  Vector3s dir_norm = direction;
  
  if(dir_norm.norm() < 1e-14)
    return;
  
  dir_norm.normalize();
  
  Vector3s perp = dir_norm.cross(view_dir);
  if(perp.norm() < 1e-14) perp = dir_norm.cross(Vector3s(0, 0, 1));
  
  perp.normalize();
  
  Vector3s tip_left = end + arrow_head_len/(scalar)sqrt(2.0)*(-dir_norm + perp);
  Vector3s tip_right = end + arrow_head_len/(scalar)sqrt(2.0)*(-dir_norm - perp);
  
  glBegin(GL_LINES);
  glVertex3dv(start.data());
  glVertex3dv(end.data());
  glVertex3dv(end.data());
  glVertex3dv(tip_left.data());
  glVertex3dv(end.data());
  glVertex3dv(tip_right.data());
  glEnd();
}

void draw_trimesh2d(const std::vector<Vector2s>& vertices, const std::vector<Vector3i>& tris) {
  glBegin(GL_TRIANGLES);
  for(unsigned int i = 0; i < tris.size(); ++i) {
    glVertex2dv(vertices[tris[i][0]].data());
    glVertex2dv(vertices[tris[i][1]].data());
    glVertex2dv(vertices[tris[i][2]].data());
  }
  glEnd();
}


void hueToRGB(scalar hue, scalar sat, scalar val, scalar &r, scalar &g, scalar &b) {
  //compute hue (adapted from an older Wikipedia article)
  int Hi = (int)(floor(hue / 60.0f)) % 6;
  scalar f = hue / 60 - Hi;
  scalar p = val * (1 - sat);
  scalar q = val * (1- f * sat);
  scalar t = val * (1 - (1 - f) * sat);
  
  switch(Hi) {
    case 0:
      r=val;
      g=t;
      b=p;
      break;
    case 1:
      r=q;
      g=val;
      b=p;
      break;
    case 2:
      r=p;
      g=val;
      b=t;
      break;
    case 3:
      r=p;
      g=q;
      b=val;
      break;
    case 4:
      r=t;
      g=p;
      b=val;
      break;
    case 5:
      r=val;
      g=p;
      b=q;
      break;
  }
}

void draw_grid_data2d(Array2s& data, Vector2s origin, scalar dx, bool color) {
  scalar max_val = FLT_MIN;
  scalar min_val = FLT_MAX;
  for(int j = 0; j < data.nj; ++j) for(int i = 0; i < data.ni; ++i) {
    max_val = std::max(data(i,j), max_val);
    min_val = std::min(data(i,j), min_val);
  }
  
  for(int j = 0; j < data.nj; ++j) {
    for(int i = 0; i < data.ni; ++i) {
      glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
      Vector2s bl = origin + Vector2s(i*dx,j*dx);
      scalar r,g,b;
      if(color) {
        hueToRGB(240*(data(i,j) - min_val)/(max_val-min_val), 1, 1, r,g,b);
      }
      else {
        scalar gray = (data(i,j) - min_val)/(max_val-min_val);
        r = g = b = gray;
      }
      //TODO Black body colormap, if I can find it.
      glColor3f(r,g,b);
      draw_box2d(bl, dx, dx);
    }
  }
  
}

void draw_trimesh3d(const std::vector<Vector3s>& vertices, const std::vector<Vector3i>& tris) {
  glBegin(GL_TRIANGLES);
  for(unsigned int i = 0; i < tris.size(); ++i) {
    glVertex3dv(vertices[tris[i][0]].data());
    glVertex3dv(vertices[tris[i][1]].data());
    glVertex3dv(vertices[tris[i][2]].data());
  }
  glEnd();
}

void write_trimesh3d(std::ostream& o, const std::vector<Vector3s>& vertices, const std::vector<Vector3i>& tris, const Vector3s& center, const Vector3s& scaling, const Eigen::Quaternion<scalar>& rot, bool double_sided)
{
  for(size_t i = 0; i < vertices.size(); ++i)
  {
    Eigen::Quaternion<scalar> rotp(0.0, vertices[i](0), vertices[i](1), vertices[i](2));
    Vector3s rot_x = (rot * rotp * rot.inverse()).vec().cwiseProduct(scaling) + center;
    
    o << "v " << rot_x.transpose() << "\n";
  }
  
  if(double_sided) {
    for(size_t i = 0; i < vertices.size(); ++i)
    {
      Eigen::Quaternion<scalar> rotp(0.0, vertices[i](0), vertices[i](1), vertices[i](2));
      Vector3s rot_x = (rot * rotp * rot.inverse()).vec().cwiseProduct(scaling * 1.01) + center;
      
      o << "v " << rot_x.transpose() << "\n";
    }
    
    for(unsigned int i = 0; i < tris.size(); ++i)
    {
      const Vector3i& tri = tris[i];
      Vector3i tri_inversed = Vector3i(tri(0), tri(2), tri(1));
      o << "f " << tri_inversed.transpose() << "\n";
    }
    
    for(unsigned int i = 0; i < tris.size(); ++i)
    {
      Vector3i tri = tris[i] + Vector3i::Constant(vertices.size());
      o << "f " << tri.transpose() << "\n";
    }
  } else {
    for(unsigned int i = 0; i < tris.size(); ++i)
    {
      o << "f " << tris[i].transpose() << "\n";
    }
  }
}

void draw_trimesh3d(const std::vector<Vector3s>& vertices, const std::vector<Vector3i>& tris, const Vector3s& center, const Vector3s& scaling, const Eigen::Quaternion<scalar>& rot)
{
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  //glLoadIdentity();
  Eigen::AngleAxis<scalar> rotp(rot);
  
  glTranslated(center(0), center(1), center(2));
  glScaled(scaling(0), scaling(1), scaling(2));
  glRotated(rotp.angle() * 180.0 / M_PI, rotp.axis()(0), rotp.axis()(1), rotp.axis()(2));
  //std::cout << rotp.axis() << " [" << rotp.angle() << "]" << std::endl;
  draw_trimesh3d(vertices, tris);
  
  glPopMatrix();
}

void draw_trimesh3d(const std::vector<Vector3s>& vertices, const std::vector<Vector3i>& tris, const std::vector<Vector3s> & normals, const Eigen::Quaternion<scalar>& rot) {
  glBegin(GL_TRIANGLES);
  for(unsigned int i = 0; i < tris.size(); ++i) {
    glNormal3dv(normals[tris[i][0]].data());
    glVertex3dv(vertices[tris[i][0]].data());
    glNormal3dv(normals[tris[i][1]].data());
    glVertex3dv(vertices[tris[i][1]].data());
    glNormal3dv(normals[tris[i][2]].data());
    glVertex3dv(vertices[tris[i][2]].data());
  }
  glEnd();
}

void write_box3d(std::vector<Vector3s>& vertices, std::vector<Vector3i>& tris)
{
  vertices.reserve(8);
  vertices.push_back(Vector3s(1.0, -1.0, -1.0));
  vertices.push_back(Vector3s(1.0, -1.0, 1.0));
  vertices.push_back(Vector3s(-1.0, -1.0, 1.0));
  vertices.push_back(Vector3s(-1.0, -1.0, -1.0));
  
  vertices.push_back(Vector3s(1.0, 1.0, -1.0));
  vertices.push_back(Vector3s(1.0, 1.0, 1.0));
  vertices.push_back(Vector3s(-1.0, 1.0, 1.0));
  vertices.push_back(Vector3s(-1.0, 1.0, -1.0));
  
  tris.reserve(12);
  tris.push_back(Vector3i(0, 1, 2));
  tris.push_back(Vector3i(0, 2, 3));
  tris.push_back(Vector3i(4, 7, 6));
  tris.push_back(Vector3i(4, 6, 5));
  
  tris.push_back(Vector3i(0, 4, 5));
  tris.push_back(Vector3i(0, 5, 1));
  tris.push_back(Vector3i(1, 5, 6));
  tris.push_back(Vector3i(1, 6, 2));
  
  tris.push_back(Vector3i(2, 6, 7));
  tris.push_back(Vector3i(2, 7, 3));
  tris.push_back(Vector3i(4, 0, 3));
  tris.push_back(Vector3i(4, 3, 7));
}

void write_box3d(std::ostream& o, const Vector3s& dimensions, const Vector3s& center, const Eigen::Quaternion<scalar>& rot, bool double_sided)
{
  scalar width = dimensions[0];
  scalar height = dimensions[1];
  scalar depth = dimensions[2];
  
  std::vector<Vector3s> verts;
  std::vector<Vector3i> tris;
  
  write_box3d(verts, tris);
  write_trimesh3d(o, verts, tris, center, Vector3s(width * 0.5, height * 0.5, depth * 0.5), rot, double_sided);
}


void draw_box3d(const Vector3s& dimensions, const Vector3s& center, const Eigen::Quaternion<scalar>& rot)
{
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  //glLoadIdentity();
  
  Eigen::AngleAxis<scalar> rotp(rot);
  
  glTranslated(center(0), center(1), center(2));
  glRotated(rotp.angle() * 180.0 / M_PI, rotp.axis()(0), rotp.axis()(1), rotp.axis()(2));
  draw_box3d(dimensions);
  glPopMatrix();
}

void draw_box3d(const Vector3s& dimensions) {
  
  //Draw an axis-aligned box with specified dimensions,
  //where the midpoint of the box is at the origin
  
  scalar width = dimensions[0];
  scalar height = dimensions[1];
  scalar depth = dimensions[2];
  
  glBegin(GL_POLYGON);
  glNormal3d(-1,0,0);
  glVertex3d(-0.5*width, -0.5*height, 0.5*depth);
  glVertex3d(-0.5*width, 0.5*height, 0.5*depth);
  glVertex3d(-0.5*width, 0.5*height, -0.5*depth);
  glVertex3d(-0.5*width, -0.5*height, -0.5*depth);
  glEnd();
  
  glBegin(GL_POLYGON);
  glNormal3d(1,0,0);
  glVertex3d(0.5*width, -0.5*height, 0.5*depth);
  glVertex3d(0.5*width, 0.5*height, 0.5*depth);
  glVertex3d(0.5*width, 0.5*height, -0.5*depth);
  glVertex3d(0.5*width, -0.5*height, -0.5*depth);
  glEnd();
  
  glBegin(GL_POLYGON);
  glNormal3d(0,0,-1);
  glVertex3d(-0.5*width, -0.5*height, -0.5*depth);
  glVertex3d(0.5*width, -0.5*height, -0.5*depth);
  glVertex3d(0.5*width, 0.5*height, -0.5*depth);
  glVertex3d(-0.5*width, 0.5*height, -0.5*depth);
  glEnd();
  
  glBegin(GL_POLYGON);
  glNormal3d(0,0,1);
  glVertex3d(-0.5*width, -0.5*height, 0.5*depth);
  glVertex3d(0.5*width, -0.5*height, 0.5*depth);
  glVertex3d(0.5*width, 0.5*height, 0.5*depth);
  glVertex3d(-0.5*width, 0.5*height, 0.5*depth);
  glEnd();
  
  glBegin(GL_POLYGON);
  glNormal3d(0,-1,0);
  glVertex3d(-0.5*width, -0.5*height, 0.5*depth);
  glVertex3d(0.5*width, -0.5*height, 0.5*depth);
  glVertex3d(0.5*width, -0.5*height, -0.5*depth);
  glVertex3d(-0.5*width, -0.5*height, -0.5*depth);
  glEnd();
  
  glBegin(GL_POLYGON);
  glNormal3d(0,1,0);
  glVertex3d(-0.5*width, 0.5*height, 0.5*depth);
  glVertex3d(0.5*width, 0.5*height, 0.5*depth);
  glVertex3d(0.5*width, 0.5*height, -0.5*depth);
  glVertex3d(-0.5*width, 0.5*height, -0.5*depth);
  glEnd();
}
