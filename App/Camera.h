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

#ifndef LIBWETHAIR_CORE_CAMERA_H_
#define LIBWETHAIR_CORE_CAMERA_H_
#include <Eigen/Dense>

#include <libWetHair/MathDefs.h>

class Camera {
 public:
  Eigen::Quaterniond rotation_;  // rotation
  Eigen::Vector3d center_;       // center point
  double dist_;                  // point of view to center
  double radius_;                // bounding sphere
  double fov_;                   // angle
  Camera(const Camera& that);
  void operator=(const Camera& that);

  explicit Camera(const double fov = 40);
  explicit Camera(Eigen::Quaterniond& rot, Eigen::Vector3d& center,
                  const double& dist, const double& radius, const double& fov);
  void init(const Eigen::Vector3d& bmin, const Eigen::Vector3d& bmax);
  void clone(const Camera& that);
  void getViewDir(Eigen::Vector3d& viewdir) const;
  void getLookAt(Eigen::Vector3d& eye, Eigen::Vector3d& center,
                 Eigen::Vector3d& up) const;
  void getEye(Eigen::Vector3d& eye) const;
  void getPerspective(double& fov, double& zNear, double& zFar) const;
  void rotate(const double oldx, const double oldy, const double newx,
              const double newy);
  void zoom(const double oldx, const double oldy, const double newx,
            const double newy);
  void pan(const double oldx, const double oldy, const double newx,
           const double newy);
  void rotate(const libwethair::Vector3s& axis, const libwethair::scalar& angle, bool global);
  void project_to_sphere(const double& radius, Eigen::Vector3d& p) const;
  friend std::ostream& operator<<(std::ostream& output, const Camera& cam);
};
#endif  // LIBWETHAIR_CORE_CAMERA_H_
