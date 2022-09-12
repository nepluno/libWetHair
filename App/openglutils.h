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

#ifndef LIBWETHAIR_APP_OPENGL_UTILS_H_
#define LIBWETHAIR_APP_OPENGL_UTILS_H_

#include <vector>

#include <libWetHair/MathDefs.h>
#include <libWetHair/array2.h>

void draw_circle2d(const libwethair::Vector2s& centre, libwethair::scalar rad, int segs);
void draw_grid2d(const libwethair::Vector2s& origin, libwethair::scalar dx, int nx, int ny);
void draw_grid3d(const libwethair::Vector3s& origin, libwethair::scalar dx, int nx, int ny, int nz);
void draw_box2d(const libwethair::Vector2s& origin, libwethair::scalar width, libwethair::scalar height);
void draw_segmentset2d(const std::vector<libwethair::Vector2s>& vertices,
                       const std::vector<libwethair::Vector2i>& edges);
void draw_points2d(const std::vector<libwethair::Vector2s>& points);
void draw_polygon2d(const std::vector<libwethair::Vector2s>& vertices);
void draw_polygon2d(const std::vector<libwethair::Vector2s>& vertices,
                    const std::vector<int>& order);
void draw_segment2d(const libwethair::Vector2s& start, const libwethair::Vector2s& end);
void draw_arrow2d(const libwethair::Vector2s& start, const libwethair::Vector2s& end,
                  libwethair::scalar arrow_head_len);
void draw_arrow3d(const libwethair::Vector3s& start, const libwethair::Vector3s& end,
                  const libwethair::Vector3s& view_dir, libwethair::scalar arrow_head_len);
void draw_grid_data2d(libwethair::Array2s& data, libwethair::Vector2s origin, libwethair::scalar dx,
                      bool color = false);
void draw_trimesh2d(const std::vector<libwethair::Vector2s>& vertices,
                    const std::vector<libwethair::Vector3i>& tris);

void draw_trimesh3d(const std::vector<libwethair::Vector3s>& vertices,
                    const std::vector<libwethair::Vector3i>& tris);
void draw_trimesh3d(const std::vector<libwethair::Vector3s>& vertices,
                    const std::vector<libwethair::Vector3i>& tris, const libwethair::Vector3s& center,
                    const libwethair::Vector3s& scaling,
                    const Eigen::Quaternion<libwethair::scalar>& rot);
void write_trimesh3d(std::ostream& o, const std::vector<libwethair::Vector3s>& vertices,
                     const std::vector<libwethair::Vector3i>& tris, const libwethair::Vector3s& center,
                     const libwethair::Vector3s& scaling,
                     const Eigen::Quaternion<libwethair::scalar>& rot, bool double_sided);
void draw_trimesh3d(const std::vector<libwethair::Vector3s>& vertices,
                    const std::vector<libwethair::Vector3i>& tris,
                    const std::vector<libwethair::Vector3s>& normals,
                    const Eigen::Quaternion<libwethair::scalar>& rot);
void draw_box3d(const libwethair::Vector3s& dimensions);
void draw_box3d(const libwethair::Vector3s& dimensions, const libwethair::Vector3s& center,
                const Eigen::Quaternion<libwethair::scalar>& rot);
void write_box3d(std::ostream& o, const libwethair::Vector3s& dimensions,
                 const libwethair::Vector3s& center, const Eigen::Quaternion<libwethair::scalar>& rot,
                 bool double_sided);
void write_box3d(std::vector<libwethair::Vector3s>& vertices, std::vector<libwethair::Vector3i>& tris);

#endif  // LIBWETHAIR_APP_OPENGL_UTILS_H_
