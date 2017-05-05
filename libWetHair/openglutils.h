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


#ifndef OPENGL_UTILS_H
#define OPENGL_UTILS_H

#include <vector>
#include "MathDefs.h"
#include "array2.h"

void draw_circle2d(const Vector2s& centre, scalar rad, int segs);
void draw_grid2d(const Vector2s& origin, scalar dx, int nx, int ny);
void draw_grid3d(const Vector3s& origin, scalar dx, int nx, int ny, int nz);
void draw_box2d(const Vector2s& origin, scalar width, scalar height);
void draw_segmentset2d(const std::vector<Vector2s>& vertices, const std::vector<Vector2i>& edges);
void draw_points2d(const std::vector<Vector2s>& points);
void draw_polygon2d(const std::vector<Vector2s>& vertices);
void draw_polygon2d(const std::vector<Vector2s>& vertices, const std::vector<int>& order);
void draw_segment2d(const Vector2s& start, const Vector2s& end);
void draw_arrow2d(const Vector2s& start, const Vector2s& end, scalar arrow_head_len);
void draw_arrow3d(const Vector3s& start, const Vector3s& end, const Vector3s& view_dir, scalar arrow_head_len);
void draw_grid_data2d(Array2s& data, Vector2s origin, scalar dx, bool color = false);
void draw_trimesh2d(const std::vector<Vector2s>& vertices, const std::vector<Vector3i>& tris);    
   
void draw_trimesh3d(const std::vector<Vector3s>& vertices, const std::vector<Vector3i>& tris);
void draw_trimesh3d(const std::vector<Vector3s>& vertices, const std::vector<Vector3i>& tris, const Vector3s& center, const Vector3s& scaling, const Eigen::Quaternion<scalar>& rot);
void write_trimesh3d(std::ostream& o, const std::vector<Vector3s>& vertices, const std::vector<Vector3i>& tris, const Vector3s& center, const Vector3s& scaling, const Eigen::Quaternion<scalar>& rot, bool double_sided);
void draw_trimesh3d(const std::vector<Vector3s>& vertices, const std::vector<Vector3i>& tris, const std::vector<Vector3s>& normals, const Eigen::Quaternion<scalar>& rot);
void draw_box3d(const Vector3s& dimensions);
void draw_box3d(const Vector3s& dimensions, const Vector3s& center, const Eigen::Quaternion<scalar>& rot);
void write_box3d(std::ostream& o, const Vector3s& dimensions, const Vector3s& center, const Eigen::Quaternion<scalar>& rot, bool double_sided);
void write_box3d(std::vector<Vector3s>& vertices, std::vector<Vector3i>& tris);

#endif
