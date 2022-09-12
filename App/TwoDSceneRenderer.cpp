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

#include "TwoDSceneRenderer.h"

#include <sstream>

#include <libWetHair/HairFlow.h>
#include <libWetHair/fluidsim2D.h>
#include <libWetHair/fluidsim3D.h>
#include <libWetHair/PolygonalCohesion.h>

#include "TwoDimensionalDisplayController.h"
#include "openglutils.h"

#define RADPERDEG 0.0174533

using namespace libwethair;

const static Vector3s liquid_color =
    Vector3s(64.0 / 255.0, 164.0 / 255.0, 223.0 / 255.0);
const static int disc_arc = 5;
const static scalar eta_epsilon = 5e-2;

std::vector<Vector3f> m_dynamic_hair_core;
std::vector<Vector3f> m_dynamic_hair_flow;
std::vector<Vector3f> m_dynamic_fluid_particles;

template <int DIM>
TwoDSceneRenderer<DIM>::TwoDSceneRenderer(
    const TwoDScene<DIM>& scene,
    const std::vector<renderingutils::Color>& particle_colors,
    const std::vector<renderingutils::Color>& edge_colors)
    : m_particle_colors(particle_colors),
      m_edge_colors(edge_colors)
      // Precomputed rendering state
      ,
      m_circle_points(),
      m_semi_circle_points(),
      m_pvm(PVM_ETA),
      m_evm(EVM_NONE),
      m_show_edge_normal(false),
      m_show_particle_normal(false),
      m_show_liquid_polygon(true),
      m_draw_grid(false),
      m_draw_boundaries(true),
      m_draw_particles(true),
      m_draw_velocities(false),
      m_scene(&scene) {
  initializeCircleRenderer(disc_arc);
  initializeSemiCircleRenderer(disc_arc);
  initializeCylinderRenderer(disc_arc, scene);
  initializeBoundaryRenderer(scene);
}

template <int DIM>
void TwoDSceneRenderer<DIM>::setDC(
    const TwoDimensionalDisplayController<DIM>* dc) {
  m_dc = dc;
}

template <int DIM>
TwoDSceneRenderer<DIM>::TwoDSceneRenderer()
    : m_pvm(PVM_NONE),
      m_evm(EVM_NONE),
      m_show_edge_normal(false),
      m_show_particle_normal(false),
      m_show_liquid_polygon(true),
      m_draw_grid(false),
      m_draw_boundaries(true),
      m_draw_particles(true),
      m_draw_velocities(false),
      m_scene(NULL) {}

template <>
void TwoDSceneRenderer<3>::initializeBoundaryRenderer(
    const TwoDScene<3>& scene) {
  m_icosphere.Create(disc_arc);

  const FluidSim3D* sim = (const FluidSim3D*)scene.getFluidSim();
  auto& boundaries = sim->get_boundaries();

  int nb = boundaries.size();
  for (int i = 0; i < nb; ++i) {
    FluidSim3D::Boundary<3>* b = boundaries[i];
    if (b->type == FluidSim3D::BT_CAPSULE) {
      FluidSim3D::SolidBoundary<3>* sb = (FluidSim3D::SolidBoundary<3>*)b;
      m_capsules.push_back(CapsuleCreator());
      CapsuleCreator& cc = m_capsules[m_capsules.size() - 1];
      cc.Create(32, sb->parameter(0), sb->parameter(1));
    } else if (b->type == FluidSim3D::BT_BOX) {
      FluidSim3D::SolidBoundary<3>* sb = (FluidSim3D::SolidBoundary<3>*)b;
      if (sb->parameter(3) > 0.0) {
        m_roundcornerboxes.push_back(RoundCornerBox());
        RoundCornerBox& rcb = m_roundcornerboxes[m_roundcornerboxes.size() - 1];
        rcb.Create(
            4, Vector3s(sb->parameter(0), sb->parameter(1), sb->parameter(2)),
            sb->parameter(3));
      }
    }
  }
}

template <>
void TwoDSceneRenderer<2>::initializeBoundaryRenderer(
    const TwoDScene<2>& scene) {}

template <int DIM>
void TwoDSceneRenderer<DIM>::initializeCircleRenderer(int num_points) {
  m_circle_points.resize(num_points);
  double dtheta = 2.0 * PI / ((double)num_points);
  for (int i = 0; i < num_points; ++i) {
    m_circle_points[i].first = cos(((double)i) * dtheta);
    m_circle_points[i].second = sin(((double)i) * dtheta);
  }
}

template <int DIM>
void TwoDSceneRenderer<DIM>::initializeSemiCircleRenderer(int num_points) {
  double dtheta = PI / ((double)(num_points - 1));
  m_semi_circle_points.resize(num_points);
  for (int i = 0; i < num_points; ++i) {
    m_semi_circle_points[i].first = -sin(((double)i) * dtheta);
    m_semi_circle_points[i].second = cos(((double)i) * dtheta);
  }
}

template <int DIM>
void TwoDSceneRenderer<DIM>::initializeCylinderRenderer(
    int num_points, const TwoDScene<DIM>& scene) {
  double dtheta = 2.0 * PI / ((double)num_points);
  m_cylinder_points.resize(num_points);
  for (int i = 0; i < num_points; ++i) {
    m_cylinder_points[i] =
        Vector3s(0.0, cos((scalar)i * dtheta), sin((scalar)i * dtheta));
  }

  auto& edges = scene.getEdges();
  int ne = scene.getNumEdges();
  int np = scene.getNumParticles();
  int nfluidp = scene.getFluidSim()->num_particles();

  m_cylinder_elements.resize(ne * num_points * 4);

  for (int k = 0; k < ne; ++k) {
    int base = k * num_points * 4;
    int vidx0 = edges[k].first;
    int vidx1 = edges[k].second;
    for (int i = 0; i < num_points; ++i) {
      int inext = (i + 1) % num_points;
      m_cylinder_elements[base + i * 4 + 0] = vidx0 * num_points + i;
      m_cylinder_elements[base + i * 4 + 1] = vidx0 * num_points + inext;
      m_cylinder_elements[base + i * 4 + 2] = vidx1 * num_points + inext;
      m_cylinder_elements[base + i * 4 + 3] = vidx1 * num_points + i;
    }
  }

  m_dynamic_hair_core.resize(num_points * np);
  m_dynamic_hair_flow.resize(num_points * np);
  m_dynamic_fluid_particles.resize(nfluidp);
}

template <int DIM>
void TwoDSceneRenderer<DIM>::initializeOpenGLRenderer(
    const TwoDScene<DIM>& scene) {
  int np = scene.getNumParticles();
  int nfluids_particles = scene.getFluidSim()->num_particles();

  glGenBuffersARB(1, &m_element_hairs);
  glGenBuffersARB(1, &m_vertex_hair_core);
  glGenBuffersARB(1, &m_vertex_hair_flow);
  glGenBuffersARB(1, &m_vertex_fluids);

  glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER, m_element_hairs);
  glBufferDataARB(GL_ELEMENT_ARRAY_BUFFER,
                  m_cylinder_elements.size() * sizeof(int),
                  &m_cylinder_elements[0], GL_STATIC_DRAW);
  glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER, 0);

  glBindBufferARB(GL_ARRAY_BUFFER, m_vertex_hair_core);
  glBufferDataARB(GL_ARRAY_BUFFER, disc_arc * np * sizeof(Vector3f), NULL,
                  GL_DYNAMIC_DRAW);

  glBindBufferARB(GL_ARRAY_BUFFER, m_vertex_hair_flow);
  glBufferDataARB(GL_ARRAY_BUFFER, disc_arc * np * sizeof(Vector3f), NULL,
                  GL_DYNAMIC_DRAW);

  glBindBufferARB(GL_ARRAY_BUFFER, m_vertex_fluids);
  glBufferDataARB(GL_ARRAY_BUFFER, nfluids_particles * sizeof(Vector3f), NULL,
                  GL_DYNAMIC_DRAW);

  glBindBufferARB(GL_ARRAY_BUFFER, 0);
}

template <>
void TwoDSceneRenderer<2>::renderSolidCircle(const Vector2s& center,
                                             double radius) const {
  glBegin(GL_TRIANGLE_FAN);
  glVertex2d(center.x(), center.y());

  for (std::vector<std::pair<double, double> >::size_type i = 0;
       i < m_circle_points.size(); ++i)
    glVertex2d(radius * m_circle_points[i].first + center.x(),
               radius * m_circle_points[i].second + center.y());

  glVertex2d(radius * m_circle_points.front().first + center.x(),
             radius * m_circle_points.front().second + center.y());
  glEnd();
}

template <>
void TwoDSceneRenderer<3>::renderSolidCircle(const Vector3s& center,
                                             double radius) const {
  glBegin(GL_TRIANGLE_FAN);
  glVertex3dv(center.data());

  for (std::vector<std::pair<double, double> >::size_type i = 0;
       i < m_circle_points.size(); ++i)
    glVertex3d(radius * m_circle_points[i].first + center.x(),
               radius * m_circle_points[i].second + center.y(), center.z());

  glVertex3d(radius * m_circle_points.front().first + center.x(),
             radius * m_circle_points.front().second + center.y(), center.z());
  glEnd();
}

template <int DIM>
void TwoDSceneRenderer<DIM>::writeTransformedHairFlow(
    std::ostream& o, const TwoDScene<DIM>& scene) const {
  // vertex
  for (auto& v : m_dynamic_hair_flow) {
    o << "v " << v(0) << " " << v(1) << " " << v(2) << std::endl;
  }

  int num_quads = m_cylinder_elements.size() / 4;

  for (int i = 0; i < num_quads; ++i) {
    o << "f " << m_cylinder_elements[i * 4 + 0] << " "
      << m_cylinder_elements[i * 4 + 1] << " " << m_cylinder_elements[i * 4 + 2]
      << " " << m_cylinder_elements[i * 4 + 3] << std::endl;
  }
}

template <>
void TwoDSceneRenderer<3>::writeBoundaries(std::ostream& os_single,
                                           std::ostream& os_double) const {
  auto& boundaries =
      ((const FluidSim3D*)m_scene->getFluidSim())->get_boundaries();

  int iCapsule = 0;
  int iRoundBox = 0;
  for (auto& b : boundaries) {
    if (b->type == FluidSim3D::BT_UNION || b->type == FluidSim3D::BT_INTERSECT)
      continue;
    FluidSim3D::SolidBoundary<3>* sb = (FluidSim3D::SolidBoundary<3>*)b;

    scalar inner_scale = 1.0;

    switch (sb->type) {
      case FluidSim3D::BT_CIRCLE:
        write_trimesh3d(sb->sign <= 0.0 ? os_double : os_single,
                        m_icosphere.vertices, m_icosphere.indices, sb->center,
                        Vector3s::Constant(sb->parameter(0) * inner_scale),
                        sb->rot, sb->sign <= 0.0);
        break;
      case FluidSim3D::BT_BOX:
        if (sb->parameter(3) > 0.0) {
          write_trimesh3d(sb->sign <= 0.0 ? os_double : os_single,
                          m_roundcornerboxes[iRoundBox].vertices,
                          m_roundcornerboxes[iRoundBox].indices, sb->center,
                          Vector3s::Constant(inner_scale), sb->rot,
                          sb->sign <= 0.0);
        } else {
          write_box3d(sb->sign <= 0.0 ? os_double : os_single,
                      sb->parameter.segment<3>(0) * 2.0 * inner_scale,
                      sb->center, sb->rot, sb->sign <= 0.0);
        }
        ++iRoundBox;
        break;
      case FluidSim3D::BT_CAPSULE:
        write_trimesh3d(sb->sign <= 0.0 ? os_double : os_single,
                        m_capsules[iCapsule].vertices,
                        m_capsules[iCapsule].indices, sb->center,
                        Vector3s::Constant(inner_scale), sb->rot,
                        sb->sign <= 0.0);
        ++iCapsule;
        break;
      default:
        break;
    }
  }
}

template <>
void TwoDSceneRenderer<2>::writeBoundaries(std::ostream& os_single,
                                           std::ostream& os_double) const {
  std::cerr << "writeReadableBoundary NOT IMPLEMENTED!" << std::endl;
}

inline Vector3f toFloatVec(const Vector3s& v) {
  return Vector3f((float)v(0), (float)v(1), (float)v(2));
}

template <int DIM>
void TwoDSceneRenderer<DIM>::updateOpenGLRenderer(const TwoDScene<DIM>& scene,
                                                  bool updateDevice) {
  const VectorXs& x = scene.getX();
  const int np = scene.getNumParticles();
  const std::vector<int>& hair_idx = scene.getParticleToHairs();
  const std::vector<int>& local_idx = scene.getParticleToHairLocalIndices();
  const auto& flows = scene.getFilmFlows();

  if (DIM == 3) {
    for (int i = 0; i < np; ++i) {
      int ihair = hair_idx[i];
      int ilocal = local_idx[i];

      const MatrixXs& vertex_dirs = flows[ihair]->getTangentV();
      const VectorXs& eta = flows[ihair]->getEta();

      const Vector3sT& ve = vertex_dirs.row(ilocal);
      const Vector3s& x0 = x.segment<3>(scene.getDof(i));
      const scalar& radius = scene.getRadius(i);
      const scalar H = radius + eta(ilocal);

      double phi0 = atan2(ve(1), ve(0));
      double length0 = ve.norm();
      double gamma0 = -asin(ve(2) / length0);

      Eigen::Affine3d trh =
          Eigen::Affine3d(Eigen::Translation3d(Vector3s(x0(0), x0(1), x0(2)))) *
          Eigen::Affine3d(Eigen::AngleAxisd(phi0, Vector3s(0, 0, 1)) *
                          Eigen::AngleAxisd(gamma0, Vector3s(0, 1, 0))) *
          Eigen::Affine3d(Eigen::Scaling(0.0, radius, radius));

      Eigen::Affine3d trf =
          Eigen::Affine3d(Eigen::Translation3d(Vector3s(x0(0), x0(1), x0(2)))) *
          Eigen::Affine3d(Eigen::AngleAxisd(phi0, Vector3s(0, 0, 1)) *
                          Eigen::AngleAxisd(gamma0, Vector3s(0, 1, 0))) *
          Eigen::Affine3d(Eigen::Scaling(0.0, H, H));

      for (int j = 0; j < disc_arc; ++j) {
        // render hairs
        m_dynamic_hair_core[i * disc_arc + j] =
            toFloatVec((trh * m_cylinder_points[j]).eval());
        m_dynamic_hair_flow[i * disc_arc + j] =
            toFloatVec((trf * m_cylinder_points[j]).eval());
      }
    }

    const auto& fluid_particles =
        ((FluidSim3D*)scene.getFluidSim())->get_particles();
    int nfp = fluid_particles.size();
    m_dynamic_fluid_particles.resize(0);
    m_dynamic_fluid_particles.reserve(nfp);
    for (int i = 0; i < nfp; ++i) {
      if (fluid_particles[i].type == PT_LIQUID) {
        m_dynamic_fluid_particles.push_back(Vector3f(fluid_particles[i].x(0),
                                                     fluid_particles[i].x(1),
                                                     fluid_particles[i].x(2)));
      }
    }

    if (updateDevice && m_dynamic_fluid_particles.size() > 0) {
      glBindBufferARB(GL_ARRAY_BUFFER, m_vertex_fluids);
      glBufferDataARB(GL_ARRAY_BUFFER,
                      m_dynamic_fluid_particles.size() * sizeof(Vector3f), NULL,
                      GL_DYNAMIC_DRAW);
      glBufferSubDataARB(GL_ARRAY_BUFFER, 0,
                         m_dynamic_fluid_particles.size() * sizeof(Vector3f),
                         &m_dynamic_fluid_particles[0]);
    }
  } else {
    for (int i = 0; i < np; ++i) {
      int ihair = hair_idx[i];
      int ilocal = local_idx[i];

      const MatrixXs& vertex_dirs = flows[ihair]->getTangentV();
      const VectorXs& eta = flows[ihair]->getEta();

      const Vector2sT& ve = vertex_dirs.row(ilocal);
      const Vector2s& x0 = x.segment<2>(scene.getDof(i));
      const scalar& radius = scene.getRadius(i);
      const scalar H = radius + eta(ilocal);

      double phi0 = atan2(ve(1), ve(0));

      Eigen::Affine3d trh =
          Eigen::Affine3d(Eigen::Translation3d(Vector3s(x0(0), x0(1), 0.0))) *
          Eigen::Affine3d(Eigen::AngleAxisd(phi0, Vector3s(0, 0, 1))) *
          Eigen::Affine3d(Eigen::Scaling(0.0, radius, radius));

      Eigen::Affine3d trf =
          Eigen::Affine3d(Eigen::Translation3d(Vector3s(x0(0), x0(1), 0.0))) *
          Eigen::Affine3d(Eigen::AngleAxisd(phi0, Vector3s(0, 0, 1))) *
          Eigen::Affine3d(Eigen::Scaling(0.0, H, H));

      for (int j = 0; j < disc_arc; ++j) {
        // render hairs
        m_dynamic_hair_core[i * disc_arc + j] =
            toFloatVec((trh * m_cylinder_points[j]).eval());
        m_dynamic_hair_flow[i * disc_arc + j] =
            toFloatVec((trf * m_cylinder_points[j]).eval());
      }
    }

    const auto& fluid_particles =
        ((FluidSim2D*)scene.getFluidSim())->get_particles();
    int nfp = fluid_particles.size();
    m_dynamic_fluid_particles.resize(0);
    m_dynamic_fluid_particles.reserve(nfp);
    for (int i = 0; i < nfp; ++i) {
      if (fluid_particles[i].type == PT_LIQUID) {
        m_dynamic_fluid_particles.push_back(
            Vector3f(fluid_particles[i].x(0), fluid_particles[i].x(1), 0.0));
      }
    }

    if (updateDevice && m_dynamic_fluid_particles.size() > 0) {
      glBindBufferARB(GL_ARRAY_BUFFER, m_vertex_fluids);
      glBufferDataARB(GL_ARRAY_BUFFER,
                      m_dynamic_fluid_particles.size() * sizeof(Vector3f), NULL,
                      GL_DYNAMIC_DRAW);
      glBufferSubDataARB(GL_ARRAY_BUFFER, 0,
                         m_dynamic_fluid_particles.size() * sizeof(Vector3f),
                         &m_dynamic_fluid_particles[0]);
    }
  }

  if (updateDevice) {
    glBindBufferARB(GL_ARRAY_BUFFER, m_vertex_hair_core);
    glBufferSubDataARB(GL_ARRAY_BUFFER, 0,
                       m_dynamic_hair_core.size() * sizeof(Vector3f),
                       &m_dynamic_hair_core[0]);

    glBindBufferARB(GL_ARRAY_BUFFER, m_vertex_hair_flow);
    glBufferSubDataARB(GL_ARRAY_BUFFER, 0,
                       m_dynamic_hair_flow.size() * sizeof(Vector3f),
                       &m_dynamic_hair_flow[0]);

    glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER, 0);
  }
}

template <>
void TwoDSceneRenderer<3>::renderParticleSimulation(
    const TwoDScene<3>& scene) const {
  assert(renderingutils::checkGLErrors());

  glEnable(GL_DEPTH_TEST);

  glViewport(0, 0, m_dc->getWindowWidth(), m_dc->getWindowHeight());

  assert(renderingutils::checkGLErrors());
  const FluidSim3D* sim = (const FluidSim3D*)scene.getFluidSim();

  if (sim) {
    if (m_draw_particles) {
      glEnable(GL_BLEND);
      glColor4d(0, 0, 1, 0.25);
      glPointSize(3);

      glEnableClientState(GL_VERTEX_ARRAY);
      glBindBufferARB(GL_ARRAY_BUFFER, m_vertex_fluids);
      glVertexPointer(3, GL_FLOAT, 0, NULL);
      glDrawArrays(GL_POINTS, 0, m_dynamic_fluid_particles.size());
      glBindBufferARB(GL_ARRAY_BUFFER, 0);
      glDisableClientState(GL_VERTEX_ARRAY);

      glDisable(GL_BLEND);
    }

    if (m_draw_grid) {
      glEnable(GL_BLEND);
      glColor4d(0.5, 0.5, 0.5, 0.08);
      glLineWidth(1);
      draw_grid3d(sim->get_origin(), sim->cellsize(), sim->get_ni(),
                  sim->get_nj(), sim->get_nk());
      glDisable(GL_BLEND);
    }

    if (m_draw_velocities) {
      glEnable(GL_BLEND);
      glColor4d(1, 0, 0, 0.25);
      Vector3s vd;
      m_dc->getCamera().getViewDir(vd);
      for (int k = 0; k < sim->get_nk(); ++k)
        for (int j = 0; j < sim->get_nj(); ++j)
          for (int i = 0; i < sim->get_ni(); ++i) {
            Vector3s pos = Vector3s((i + 0.5) * sim->cellsize(),
                                    (j + 0.5) * sim->cellsize(),
                                    (k + 0.5) * sim->cellsize()) +
                           sim->get_origin();
            draw_arrow3d(pos, pos + 0.01 * sim->get_velocity(pos), vd,
                         0.1 * sim->cellsize());
          }
      glDisable(GL_BLEND);
      assert(renderingutils::checkGLErrors());
    }

    if (m_draw_boundaries) {
      int iCapsule = 0;
      int iRoundBox = 0;
      auto& boundaries = sim->get_boundaries();

      for (auto& b : boundaries) {
        if (b->type == FluidSim3D::BT_UNION ||
            b->type == FluidSim3D::BT_INTERSECT)
          continue;
        FluidSim3D::SolidBoundary<3>* sb = (FluidSim3D::SolidBoundary<3>*)b;

        if (sb->sign < 0.0) {
          glEnable(GL_BLEND);
          glColor4d(0.650980392156863, 0.294117647058824, 0.0, 0.05);
        } else {
          glDisable(GL_BLEND);
          glColor3d(1.0, 1.0, 1.0);
        }

        scalar inner_scale = 1.0;  //(sb->sign > 0.0 ? 0.95 : 1.01);
        scalar outer_scale = 1.0;  // 0.96;

        switch (sb->type) {
          case FluidSim3D::BT_CIRCLE:
            if (sb->sign < 0.0) {
              glDisable(GL_DEPTH_TEST);
            } else {
              glEnable(GL_POLYGON_STIPPLE);
            }
            draw_trimesh3d(m_icosphere.vertices, m_icosphere.indices,
                           sb->center,
                           Vector3s(sb->parameter(0) * inner_scale,
                                    sb->parameter(0) * inner_scale,
                                    sb->parameter(0) * inner_scale),
                           sb->rot);
            glDisable(GL_POLYGON_STIPPLE);
            glEnable(GL_DEPTH_TEST);
            if (sb->sign > 0.0) {
              glColor3d(0.85, 0.85, 0.85);
              glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
              draw_trimesh3d(m_icosphere.vertices, m_icosphere.indices,
                             sb->center,
                             Vector3s(sb->parameter(0) * outer_scale,
                                      sb->parameter(0) * outer_scale,
                                      sb->parameter(0) * outer_scale),
                             sb->rot);
              glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
            }
            break;
          case FluidSim3D::BT_BOX:
            if (sb->sign < 0.0) {
              glDisable(GL_DEPTH_TEST);
            } else {
              glEnable(GL_POLYGON_STIPPLE);
            }
            if (sb->parameter(3) > 0.0) {
              draw_trimesh3d(m_roundcornerboxes[iRoundBox].vertices,
                             m_roundcornerboxes[iRoundBox].indices, sb->center,
                             Vector3s(inner_scale, inner_scale, inner_scale),
                             sb->rot);
            } else {
              draw_box3d(Vector3s(sb->parameter(0), sb->parameter(1),
                                  sb->parameter(2)) *
                             2.0 * inner_scale,
                         sb->center, sb->rot);
            }
            glDisable(GL_POLYGON_STIPPLE);
            glEnable(GL_DEPTH_TEST);
            if (sb->sign > 0.0) {
              glColor3d(0.85, 0.85, 0.85);
              glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
              if (sb->parameter(3) > 0.0) {
                draw_trimesh3d(
                    m_roundcornerboxes[iRoundBox].vertices,
                    m_roundcornerboxes[iRoundBox].indices, sb->center,
                    Vector3s(outer_scale, outer_scale, outer_scale), sb->rot);
              } else {
                draw_box3d(Vector3s(sb->parameter(0), sb->parameter(1),
                                    sb->parameter(2)) *
                               2.0 * outer_scale,
                           sb->center, sb->rot);
              }
              glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
            }
            if (sb->parameter(3) > 0.0) {
              ++iRoundBox;
            }
            break;
          case FluidSim3D::BT_CAPSULE:
            if (sb->sign < 0.0) {
              glDisable(GL_DEPTH_TEST);
            } else {
              glEnable(GL_POLYGON_STIPPLE);
            }
            draw_trimesh3d(m_capsules[iCapsule].vertices,
                           m_capsules[iCapsule].indices, sb->center,
                           Vector3s(inner_scale, inner_scale, inner_scale),
                           sb->rot);
            glDisable(GL_POLYGON_STIPPLE);
            glEnable(GL_DEPTH_TEST);
            if (sb->sign > 0.0) {
              glColor3d(0.85, 0.85, 0.85);
              glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
              draw_trimesh3d(m_capsules[iCapsule].vertices,
                             m_capsules[iCapsule].indices, sb->center,
                             Vector3s(outer_scale, outer_scale, outer_scale),
                             sb->rot);
              glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
            }
            ++iCapsule;
            break;
          default:
            break;
        }
      }
      glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
      glDisable(GL_BLEND);
    }
  }

  assert(renderingutils::checkGLErrors());

  glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER, m_element_hairs);

  glColor4d(0, 0, 0, 1.0);

  glEnableClientState(GL_VERTEX_ARRAY);
  glBindBufferARB(GL_ARRAY_BUFFER, m_vertex_hair_core);
  glVertexPointer(3, GL_FLOAT, 0, NULL);

  glDrawElements(GL_QUADS, (int)m_cylinder_elements.size(), GL_UNSIGNED_INT, 0);

  // render hair flow
  glColor4d(liquid_color(0), liquid_color(1), liquid_color(2), 0.2);
  glEnable(GL_POLYGON_STIPPLE);

  glBindBufferARB(GL_ARRAY_BUFFER, m_vertex_hair_flow);
  glVertexPointer(3, GL_FLOAT, 0, NULL);

  glDrawElements(GL_QUADS, (int)m_cylinder_elements.size(), GL_UNSIGNED_INT, 0);

  glBindBufferARB(GL_ARRAY_BUFFER, 0);
  glDisableClientState(GL_VERTEX_ARRAY);

  glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER, 0);

  glDisable(GL_POLYGON_STIPPLE);

  if (scene.doVolSummary()) {
    renderVolumeGraph(scene);
  }

  glViewport(0, 0, m_dc->getWindowWidth(), m_dc->getWindowHeight());
}

template <>
void TwoDSceneRenderer<2>::renderParticleSimulation(
    const TwoDScene<2>& scene) const {
  assert(renderingutils::checkGLErrors());

  glViewport(0, 0, m_dc->getWindowWidth(), m_dc->getWindowHeight());

  assert(renderingutils::checkGLErrors());
  const FluidSim2D* sim = (const FluidSim2D*)scene.getFluidSim();

  if (sim) {
    assert(renderingutils::checkGLErrors());

    if (m_draw_grid) {
      glColor3d(0.96, 0.96, 0.96);
      glLineWidth(1);
      draw_grid2d(sim->get_origin(), sim->cellsize(), sim->get_ni(),
                  sim->get_nj());
      assert(renderingutils::checkGLErrors());
    }

    glColor3d(0.650980392156863, 0.294117647058824, 0.0);

    if (m_draw_boundaries) {
      glLineWidth(3);
      glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
      auto& boundaries = sim->get_boundaries();

      for (auto& b : boundaries) {
        if (b->type == FluidSim2D::BT_UNION ||
            b->type == FluidSim2D::BT_INTERSECT)
          continue;
        FluidSim2D::SolidBoundary<2>* sb = (FluidSim2D::SolidBoundary<2>*)b;

        switch (sb->type) {
          case FluidSim2D::BT_CIRCLE:
            draw_circle2d(sb->center, sb->parameter(0), 128);
            break;
          case FluidSim2D::BT_BOX:
            draw_box2d(
                sb->center - Vector2s(sb->parameter(0), sb->parameter(1)),
                sb->parameter(0) * 2., sb->parameter(1) * 2.);
            break;
          default:
            break;
        }
      }
      glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
      glLineWidth(1);
      assert(renderingutils::checkGLErrors());
    }

    if (m_draw_particles) {
      auto& particles = sim->get_particles();

      glColor3d(0, 0, 1);
      glPointSize(5);

      glEnableClientState(GL_VERTEX_ARRAY);
      glBindBufferARB(GL_ARRAY_BUFFER, m_vertex_fluids);
      glVertexPointer(3, GL_FLOAT, 0, NULL);
      glDrawArrays(GL_POINTS, 0, m_dynamic_fluid_particles.size());

      glBindBufferARB(GL_ARRAY_BUFFER, 0);
      glDisableClientState(GL_VERTEX_ARRAY);

      assert(renderingutils::checkGLErrors());
    }

    if (m_draw_velocities) {
      glEnable(GL_BLEND);
      glColor4d(1, 0, 0, 0.25);
      for (int j = 0; j < sim->get_nj(); ++j)
        for (int i = 0; i < sim->get_ni(); ++i) {
          Vector2s pos = Vector2s((i + 0.5) * sim->cellsize(),
                                  (j + 0.5) * sim->cellsize()) +
                         sim->get_origin();
          draw_arrow2d(pos, pos + 0.01 * sim->get_velocity(pos),
                       0.1 * sim->cellsize());
        }
      glDisable(GL_BLEND);
      assert(renderingutils::checkGLErrors());
    }
  }

  assert(renderingutils::checkGLErrors());

  glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER, m_element_hairs);

  glColor4d(0, 0, 0, 1.0);

  glEnableClientState(GL_VERTEX_ARRAY);
  glBindBufferARB(GL_ARRAY_BUFFER, m_vertex_hair_core);
  glVertexPointer(3, GL_FLOAT, 0, NULL);

  glDrawElements(GL_QUADS, (int)m_cylinder_elements.size(), GL_UNSIGNED_INT, 0);

  // render hair flow
  glColor4d(liquid_color(0), liquid_color(1), liquid_color(2), 0.2);
  glEnable(GL_POLYGON_STIPPLE);

  glBindBufferARB(GL_ARRAY_BUFFER, m_vertex_hair_flow);
  glVertexPointer(3, GL_FLOAT, 0, NULL);

  glDrawElements(GL_QUADS, (int)m_cylinder_elements.size(), GL_UNSIGNED_INT, 0);

  glBindBufferARB(GL_ARRAY_BUFFER, 0);
  glDisableClientState(GL_VERTEX_ARRAY);

  glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER, 0);

  glDisable(GL_POLYGON_STIPPLE);

  if (scene.doVolSummary()) {
    renderVolumeGraph(scene);
  }

  glViewport(0, 0, m_dc->getWindowWidth(), m_dc->getWindowHeight());
}

template <int DIM>
void TwoDSceneRenderer<DIM>::renderVolumeGraph(
    const TwoDScene<DIM>& scene) const {
  renderingutils::setOrthographicProjection();
  const std::vector<scalar>& hair_vol = scene.getHairVolSummary();
  const std::vector<scalar>& free_vol = scene.getFreeVolSummary();

  int window = std::min((int)(1.0 / scene.getDt()), (int)hair_vol.size());
  int start_idx = hair_vol.size() - window;

  scalar stretch_width = (scalar)m_dc->getWindowWidth() * 0.15;
  scalar stretch_height = (scalar)m_dc->getWindowHeight() * 0.15;

  glDisable(GL_DEPTH_TEST);

  glColor3d(0.0, 0.0, 0.0);
  glLineWidth(3.0);
  glBegin(GL_LINES);
  for (int i = 0; i < window - 1; ++i) {
    scalar sx = (scalar)i / (scalar)window * stretch_width;
    scalar fvc = (hair_vol[start_idx + i] + free_vol[start_idx + i]) /
                 (free_vol[0] + hair_vol[0]);
    scalar sy = (fvc + 0.25) * stretch_height;
    glVertex2d(sx, sy);

    scalar sx1 = (scalar)(i + 1) / (scalar)window * stretch_width;
    scalar sy1 = ((hair_vol[start_idx + i + 1] + free_vol[start_idx + i + 1]) /
                      (free_vol[0] + hair_vol[0]) +
                  0.25) *
                 stretch_height;

    glVertex2d(sx1, sy1);
  }
  glEnd();
  glLineWidth(1.0);

  glColor3d(1.0, 0.0, 0.0);
  glBegin(GL_LINES);
  for (int i = 0; i < window - 1; ++i) {
    scalar sx = (scalar)i / (scalar)window * stretch_width;
    scalar fvc = (free_vol[start_idx + i]) / (free_vol[0] + hair_vol[0]);
    scalar sy = (fvc + 0.25) * stretch_height;
    glVertex2d(sx, sy);

    scalar sx1 = (scalar)(i + 1) / (scalar)window * stretch_width;
    scalar fvc1 = (free_vol[start_idx + i + 1]) / (free_vol[0] + hair_vol[0]);
    scalar sy1 = (fvc1 + 0.25) * stretch_height;
    glVertex2d(sx1, sy1);
  }
  glEnd();

  glColor3d(0.0, 0.0, 1.0);
  glBegin(GL_LINES);
  for (int i = 0; i < window - 1; ++i) {
    scalar sx = (scalar)i / (scalar)window * stretch_width;
    scalar sy = (hair_vol[start_idx + i] / (free_vol[0] + hair_vol[0]) + 0.25) *
                stretch_height;
    glVertex2d(sx, sy);

    scalar sx1 = (scalar)(i + 1) / (scalar)window * stretch_width;
    scalar sy1 =
        (hair_vol[start_idx + i + 1] / (free_vol[0] + hair_vol[0]) + 0.25) *
        stretch_height;
    glVertex2d(sx1, sy1);
  }
  glEnd();

  glEnable(GL_DEPTH_TEST);

  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
}

template <int DIM>
int TwoDSceneRenderer<DIM>::renderDebuggingInfo(
    const HairFlow<DIM>& flow) const {
  int line = 2 * flow.index() + 1;
  if (m_pvm != PVM_NONE) {
    std::ostringstream oss;
    switch (m_pvm) {
      case PVM_ETA:
        // oss << "Flow " << flow.index() << " - vertex eta (Min: " <<
        // flow.getMinEta() << ", Max: " << flow.getMaxEta() << ", Pool: " <<
        // flow.getPoolSize() << ")";
        oss << "Flow " << flow.index()
            << " - vertex eta (Min: " << flow.getMinEta()
            << ", Max: " << flow.getMaxEta() << ")";
        break;

      default:
        break;
    }

    renderingutils::drawHUDString(line++, oss.str());
  }

  if (m_evm != EVM_NONE) {
    std::ostringstream oss;
    switch (m_evm) {
      case EVM_AREA:
        oss << "Flow " << flow.index()
            << " - edge area (Min: " << flow.getMinAreaE()
            << ", Max: " << flow.getMaxAreaE() << ")";
        break;
      default:
        break;
    }

    renderingutils::drawHUDString(line++, oss.str());
  }

  return line;
}

template <int DIM>
void TwoDSceneRenderer<DIM>::selectNextParticleVisMode() {
  m_pvm = (m_pvm + 1) % PVM_COUNT;
}

template <int DIM>
void TwoDSceneRenderer<DIM>::selectNextEdgeVisMode() {
  m_evm = (m_evm + 1) % EVM_COUNT;
}

template <int DIM>
void TwoDSceneRenderer<DIM>::switchShowParticleNormal() {
  m_show_particle_normal = !m_show_particle_normal;
}

template <int DIM>
void TwoDSceneRenderer<DIM>::switchShowEdgeNormal() {
  m_show_edge_normal = !m_show_edge_normal;
}

template <int DIM>
void TwoDSceneRenderer<DIM>::switchShowLiquidPolygon() {
  m_show_liquid_polygon = !m_show_liquid_polygon;
}

template <int DIM>
std::vector<renderingutils::Color>&
TwoDSceneRenderer<DIM>::getParticleColors() {
  return m_particle_colors;
}

template <int DIM>
const std::vector<renderingutils::Color>&
TwoDSceneRenderer<DIM>::getParticleColors() const {
  return m_particle_colors;
}

// explicit instantiations at bottom
template class TwoDSceneRenderer<2>;
template class TwoDSceneRenderer<3>;
