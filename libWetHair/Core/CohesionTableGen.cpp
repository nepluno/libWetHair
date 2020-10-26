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

#include "CohesionTableGen.h"

#include <iostream>
#include <numeric>

#include "MathUtilities.h"

#define COLLISION_COHESION_TABLE

inline scalar spring_func(const scalar& x, const scalar& dmin, const scalar& k0,
                          const scalar& k1, const scalar& x1, const scalar& h) {
  const scalar t2 = 1.0 / h;
  const scalar t3 = x * x;
  const scalar t4 = dmin + h;
  const scalar t5 = h * h;
  const scalar t6 = k1 * t5;
  const scalar t7 = dmin * h * k1;
  return -t2 - t2 * t3 * 1.0 / (t4 * t4) * (t6 + t7 - h * x1 * 3.0 - 3.0) +
         t2 * t3 * 1.0 / (t4 * t4 * t4) * x * (t6 + t7 - h * x1 * 2.0 - 2.0);
}

inline scalar grad_spring_func(const scalar& x, const scalar& dmin,
                               const scalar& k0, const scalar& k1,
                               const scalar& x1, const scalar& h) {
  const scalar t2 = 1.0 / h;
  const scalar t3 = dmin + h;
  const scalar t4 = h * h;
  const scalar t5 = k1 * t4;
  const scalar t6 = dmin * h * k1;
  return t2 * 1.0 / (t3 * t3 * t3) * (x * x) * (t5 + t6 - h * x1 * 2.0 - 2.0) *
             3.0 -
         t2 * 1.0 / (t3 * t3) * x * (t5 + t6 - h * x1 * 3.0 - 3.0) * 2.0;
}

CohesionTable::CohesionTable(const scalar& radius_multiplier,
                             const scalar& collision_stiffness,
                             const scalar& radius_multiplier_planar,
                             const scalar& collision_stiffness_planar)
    : m_sigma(72.0),
      m_theta(M_PI / 4),
      m_radii(0.004),
      m_max_alpha(M_PI - m_theta),
      m_max_d0(0.2),
      m_min_d0(2.0 * m_radii),
      m_min_d0_planar(m_radii),
      m_radius_multiplier(radius_multiplier),
      m_collision_stiffness(collision_stiffness),
      m_collision_stiffness_planar(collision_stiffness_planar),
      m_radius_multiplier_planar(radius_multiplier_planar),
      m_discretization(256) {}

void CohesionTable::print_energy_data(std::ostream& oss,
                                      bool first_time) const {
  scalar d_inc = m_max_d0 / (scalar)m_discretization;
  scalar A_max = m_max_As.maxCoeff();
  scalar A_inc = A_max / (scalar)m_discretization;

  MatrixXs E(m_discretization, m_discretization);

  for (int i = 0; i < m_discretization; ++i) {
    scalar A = (scalar)i * A_inc;
    scalar d_breaking;
    if (first_time) {
      d_breaking = sqrt(2.0 * A / M_PI + 4.0 * m_radii * m_radii);
    } else {
      d_breaking = (1.0 + 0.5 * m_theta) * sqrt(A);
    }

    for (int j = 0; j < m_discretization; ++j) {
      scalar d0 = (scalar)j * d_inc;

      if (A == 0.0 || d0 <= d_breaking) {
        scalar k = getStiffness(d0, A, 1.0);

        E(i, j) = 0.5 * k * (d0 - getDStar()) * (d0 - getDStar());
      } else {
        scalar A0 = A * 0.5;
        scalar H0 = sqrt(A0 / M_PI);
        scalar S0 = 2.0 * M_PI * H0;
        E(i, j) = 2.0 * S0 * m_sigma;
      }
    }
  }

  oss << "Max D: " << m_max_d0 << " Max A: " << A_max << std::endl;
  oss << E << std::endl;
}

void CohesionTable::print_table(std::ostream& oss, const MatrixXs& mat,
                                const scalar& dmin) const {
  scalar d_inc = m_max_d0 / (scalar)m_discretization;

  for (int j = 0; j < m_discretization; ++j) {
    for (int i = 0; i < m_discretization; ++i) {
      scalar d0 = (scalar)i * d_inc + dmin;
      scalar a_inc = (m_max_As(i) - m_min_As(i)) / m_discretization;

      scalar A = (scalar)j * a_inc + m_min_As(i);

      scalar value = mat(j, i);

      oss << A << " " << d0 << " " << value << std::endl;
    }
  }
}

void CohesionTable::print_dEdd_table(std::ostream& oss) const {
  print_table(oss, m_dEdd_table, 0.0);
}

scalar CohesionTable::interpolate_table(const scalar& A, const scalar& d0,
                                        const MatrixXs& mat,
                                        const scalar& dmin) const {
  scalar d_inc = m_max_d0 / (scalar)m_discretization;
  scalar p = std::max(0., d0 - dmin) / d_inc;
  scalar fp = std::min(p - floor(p), 1.0);
  int ip0 = std::max(0, std::min(m_discretization - 2, (int)p));
  int ip1 = ip0 + 1;

  scalar a_inc0 = (m_max_As(ip0) - m_min_As(ip0)) / m_discretization;
  scalar a_inc1 = (m_max_As(ip1) - m_min_As(ip1)) / m_discretization;

  if (a_inc0 == 0.0 || a_inc1 == 0.0) return 0.0;

  scalar q0 =
      (std::max(m_min_As(ip0), std::min(m_max_As(ip0), A)) - m_min_As(ip0)) /
      a_inc0;
  scalar q1 =
      (std::max(m_min_As(ip1), std::min(m_max_As(ip1), A)) - m_min_As(ip1)) /
      a_inc1;

  scalar fq0 = std::min(q0 - floor(q0), 1.0);
  scalar fq1 = std::min(q1 - floor(q1), 1.0);

  int iq00 = std::max(0, std::min(m_discretization - 2, (int)q0));
  int iq10 = std::max(0, std::min(m_discretization - 2, (int)q1));

  int iq01 = iq00 + 1;
  int iq11 = iq10 + 1;

  scalar dEdd0;
  scalar dEdd1;

  scalar v00 = mat(iq00, ip0);
  scalar v01 = mat(iq01, ip0);
  scalar v10 = mat(iq10, ip1);
  scalar v11 = mat(iq11, ip1);

  dEdd0 = mathutils::lerp(v00, v01, fq0);
  dEdd1 = mathutils::lerp(v10, v11, fq1);

  if (std::isnan(dEdd0)) {
    std::cout << "dEdd0 is NAN!" << std::endl;
    std::cout << "A: " << A << std::endl;
    std::cout << "d0: " << d0 << std::endl;
    std::cout << "mat: " << mat << std::endl;
    std::cout << "dmin: " << dmin << std::endl;
    std::cout << "fp: " << fp << std::endl;

    exit(0);
  }

  if (std::isnan(dEdd1)) {
    std::cout << "dEdd1 is NAN!" << std::endl;
    std::cout << "A: " << A << std::endl;
    std::cout << "d0: " << d0 << std::endl;
    std::cout << "mat: " << mat << std::endl;
    std::cout << "dmin: " << dmin << std::endl;
    std::cout << "fp: " << fp << std::endl;

    exit(0);
  }

  return mathutils::lerp(dEdd0, dEdd1, fp);
}

scalar CohesionTable::interpolate_table_grad(const scalar& A, const scalar& d0,
                                             const MatrixXs& mat,
                                             const scalar& dmin) const {
  scalar d_inc = m_max_d0 / (scalar)m_discretization;
  scalar p = std::max(0., d0 - dmin) / d_inc;
  int ip0 = std::max(0, std::min(m_discretization - 2, (int)p));
  int ip1 = ip0 + 1;

  scalar a_inc0 = (m_max_As(ip0) - m_min_As(ip0)) / m_discretization;
  scalar a_inc1 = (m_max_As(ip1) - m_min_As(ip1)) / m_discretization;

  if (a_inc0 == 0.0 || a_inc1 == 0.0) return 0.0;

  scalar q0 =
      (std::max(m_min_As(ip0), std::min(m_max_As(ip0), A)) - m_min_As(ip0)) /
      a_inc0;
  scalar q1 =
      (std::max(m_min_As(ip1), std::min(m_max_As(ip1), A)) - m_min_As(ip1)) /
      a_inc1;

  scalar fq0 = std::min(q0 - floor(q0), 1.0);
  scalar fq1 = std::min(q1 - floor(q1), 1.0);

  int iq00 = std::max(0, std::min(m_discretization - 2, (int)q0));
  int iq10 = std::max(0, std::min(m_discretization - 2, (int)q1));

  int iq01 = iq00 + 1;
  int iq11 = iq10 + 1;

  scalar dEdd0;
  scalar dEdd1;

  scalar v00 = mat(iq00, ip0);
  scalar v01 = mat(iq01, ip0);
  scalar v10 = mat(iq10, ip1);
  scalar v11 = mat(iq11, ip1);

  dEdd0 = mathutils::lerp(v00, v01, fq0);
  dEdd1 = mathutils::lerp(v10, v11, fq1);

  return (dEdd1 - dEdd0) / d_inc;
}

scalar CohesionTable::interpolate_dEdd(const scalar& A,
                                       const scalar& d0) const {
  return interpolate_table(A, d0, m_dEdd_table, 0.0);
}

scalar CohesionTable::interpolate_d2Edd2(const scalar& A,
                                         const scalar& d0) const {
  return interpolate_table_grad(A, d0, m_dEdd_table, 0.0);
}

scalar CohesionTable::interpolate_alpha(const scalar& A,
                                        const scalar& d0) const {
  return interpolate_table(A, d0, m_alpha_table, m_min_d0);
}

scalar CohesionTable::interpolate_dEdd_planar(const scalar& A,
                                              const scalar& d0) const {
  return interpolate_table(A, d0, m_dEdd_planar_table, 0.0);
}

scalar CohesionTable::interpolate_d2Edd2_planar(const scalar& A,
                                                const scalar& d0) const {
  return interpolate_table_grad(A, d0, m_dEdd_planar_table, 0.0);
}

scalar CohesionTable::interpolate_alpha_planar(const scalar& A,
                                               const scalar& d0) const {
  return interpolate_table(A, d0, m_alpha_planar_table, m_min_d0_planar);
}

void CohesionTable::setParameter(const scalar& sigma, const scalar& theta,
                                 const scalar& radii, const scalar& max_d0,
                                 const int disc) {
  m_sigma = sigma;
  m_theta = theta;
  m_radii = radii;
  m_max_d0 = max_d0;
  m_discretization = disc;
  m_min_d0 = radii * 2.0;
  m_max_alpha = M_PI - m_theta;
}

scalar CohesionTable::getRadii() const { return m_radii; }

scalar CohesionTable::computeH(const scalar& R, const scalar& alpha) const {
  return m_radii * sin(alpha) - R * (1.0 - sin(m_theta + alpha));
}

scalar CohesionTable::computeR(const scalar& alpha, const scalar& d0) const {
  return (d0 - 2.0 * m_radii * cos(alpha)) / (2.0 * cos(m_theta + alpha));
}

scalar CohesionTable::computeA(const scalar& R, const scalar& alpha) const {
  return 2.0 * R * R *
             (alpha + m_theta - M_PI / 2 + 0.5 * sin(2.0 * (alpha + m_theta))) +
         2.0 * m_radii * R * (sin(2.0 * alpha + m_theta) - sin(m_theta)) -
         m_radii * m_radii * (2.0 * alpha - sin(2.0 * alpha));
}

scalar CohesionTable::computeApproxA(const scalar& alpha,
                                     const scalar& d0) const {
  const scalar gamma = alpha + m_theta;
  const scalar t2 = m_radii * m_radii;
  const scalar t3 = d0 * d0;
  const scalar t4 = cos(m_theta);
  const scalar t5 = t4 * t4;
  const scalar t6 = sin(m_theta);
  return -t2 * sin(m_theta * 2.0) + t2 * M_PI * (1.0 / 3.0) -
         t3 * M_PI * (1.0 / 6.0) - gamma * t2 * (8.0 / 3.0) +
         gamma * t3 * (1.0 / 3.0) + t2 * m_theta * 2.0 -
         t2 * t5 * M_PI * (4.0 / 3.0) + d0 * m_radii * t4 * 2.0 +
         gamma * t2 * t5 * (8.0 / 3.0) +
         d0 * gamma * m_radii * t6 * (2.0 / 3.0) -
         d0 * m_radii * t6 * M_PI * (1.0 / 3.0);
}

scalar CohesionTable::computeApproxdEdd(const scalar& alpha,
                                        const scalar& d0) const {
  const scalar gamma = alpha + m_theta;

  const scalar t2 = sin(m_theta);
  const scalar t3 = m_radii * m_radii;
  const scalar t4 = cos(m_theta);
  const scalar t5 = d0 * d0;
  const scalar t6 = d0 * m_radii * t2 * 2.0;
  const scalar t7 = m_theta * 2.0;
  const scalar t8 = sin(t7);
  return (m_sigma *
          (t3 * -8.0 + t5 + t6 + t3 * (t4 * t4) * 8.0 + t3 * t8 * M_PI * 2.0 -
           gamma * t3 * t8 * 4.0 - d0 * gamma * m_radii * t4 * 2.0 +
           d0 * m_radii * t4 * M_PI) *
          2.0) /
         (t5 + t6 - (t2 * t2) * t3 * 8.0);
}

scalar CohesionTable::computedEddPlanar(const scalar& R,
                                        const scalar& alpha) const {
  if (R == 0.0) {
    return 0.0;
  } else {
    const scalar t2 = m_theta * 2.0;
    const scalar t3 = sin(alpha);
    const scalar t4 = alpha + m_theta;
    const scalar t5 = sin(t4);
    const scalar t6 = alpha + t2;
    const scalar t7 = m_radii * m_radii;
    const scalar t8 = sin(t6);
    const scalar t9 = R * R;
    const scalar t10 = alpha * 2.0;
    const scalar t11 = sin(t2);
    const scalar t12 = t2 + t10;
    const scalar t13 = sin(t12);
    const scalar t14 = m_theta * 3.0;
    const scalar t15 = alpha + t14;
    const scalar t16 = cos(t10);
    const scalar t17 = cos(t12);
    const scalar t18 = cos(m_theta);
    const scalar t19 = t10 + m_theta;
    const scalar t20 = cos(t19);
    return (m_sigma *
            (t7 * M_PI * -2.0 - t9 * M_PI * 2.0 + alpha * t7 * 2.0 +
             alpha * t9 * 2.0 + t3 * t7 * 4.0 + t3 * t9 * 4.0 + t7 * t8 * 2.0 +
             t8 * t9 * 4.0 - t7 * t11 * 2.0 + t7 * t13 * 2.0 - t9 * t11 +
             t9 * t13 * 3.0 + t7 * m_theta * 4.0 + t9 * m_theta * 4.0 +
             t7 * sin(t10) * 2.0 + t7 * sin(alpha - t2) * 2.0 +
             t9 * sin(t10 + m_theta * 4.0) - R * m_radii * sin(t14) +
             R * m_radii * sin(t15) * 2.0 + R * m_radii * sin(t19) * 5.0 -
             R * m_radii * sin(m_theta) * 3.0 +
             R * m_radii * sin(alpha - m_theta) * 6.0 + t7 * t16 * M_PI * 2.0 +
             t9 * t17 * M_PI * 2.0 + R * m_radii * t5 * 8.0 -
             alpha * t7 * t16 * 2.0 - alpha * t9 * t17 * 2.0 -
             t7 * t16 * m_theta * 4.0 - t9 * t17 * m_theta * 4.0 +
             R * m_radii * sin(t10 + t14) * 3.0 +
             R * m_radii * t18 * m_theta * 8.0 -
             R * m_radii * t20 * m_theta * 8.0 -
             R * m_radii * t18 * M_PI * 4.0 + R * m_radii * t20 * M_PI * 4.0 +
             R * alpha * m_radii * t18 * 4.0 -
             R * alpha * m_radii * t20 * 4.0)) /
           (R * (m_radii * 2.0 + R * t18 * 4.0 + R * cos(t4) * 3.0 +
                 R * cos(t15) + m_radii * cos(alpha) * 2.0 +
                 m_radii * cos(t2) * 2.0 + m_radii * cos(t6) * 2.0 -
                 R * t5 * M_PI * 2.0 + R * alpha * t5 * 2.0 -
                 m_radii * t3 * M_PI * 2.0 + R * t5 * m_theta * 4.0 +
                 alpha * m_radii * t3 * 2.0 + m_radii * t3 * m_theta * 4.0));
  }
}

scalar CohesionTable::computedEdd(const scalar& R, const scalar& alpha) const {
  if (R == 0.0) {
    return 0.0;
  } else {
    const scalar t2 = sin(alpha);
    const scalar t3 = alpha + m_theta;
    const scalar t4 = sin(t3);
    const scalar t5 = R * R;
    const scalar t6 = m_radii * m_radii;
    const scalar t7 = m_theta * 2.0;
    const scalar t8 = alpha * 2.0;
    const scalar t9 = t7 + t8;
    const scalar t10 = sin(t9);
    const scalar t11 = cos(t8);
    const scalar t12 = cos(t9);
    const scalar t13 = cos(m_theta);
    const scalar t14 = t8 + m_theta;
    const scalar t15 = cos(t14);
    return (m_sigma *
            (-t5 * M_PI - t6 * M_PI + alpha * t5 * 2.0 + alpha * t6 * 2.0 +
             t5 * t10 * 2.0 + t6 * t10 + t5 * m_theta * 2.0 +
             t6 * m_theta * 2.0 - t6 * sin(t7) + t6 * sin(t8) +
             R * m_radii * sin(t14) * 3.0 - R * m_radii * sin(m_theta) * 2.0 +
             R * m_radii * sin(t8 + m_theta * 3.0) + t5 * t12 * M_PI +
             t6 * t11 * M_PI - alpha * t5 * t12 * 2.0 - alpha * t6 * t11 * 2.0 -
             t5 * t12 * m_theta * 2.0 - t6 * t11 * m_theta * 2.0 +
             R * m_radii * t13 * m_theta * 4.0 -
             R * m_radii * t15 * m_theta * 4.0 -
             R * m_radii * t13 * M_PI * 2.0 + R * m_radii * t15 * M_PI * 2.0 +
             R * alpha * m_radii * t13 * 4.0 -
             R * alpha * m_radii * t15 * 4.0)) /
           (R * (m_radii * cos(alpha + t7) + R * cos(t3) * 2.0 +
                 m_radii * cos(alpha) - R * t4 * M_PI + R * alpha * t4 * 2.0 -
                 m_radii * t2 * M_PI + R * t4 * m_theta * 2.0 +
                 alpha * m_radii * t2 * 2.0 + m_radii * t2 * m_theta * 2.0));
  }
}

scalar CohesionTable::computeRPlanar(const scalar& alpha,
                                     const scalar& d0) const {
  return (d0 - m_radii * cos(alpha)) / (cos(m_theta + alpha) + cos(m_theta));
}

scalar CohesionTable::computeAPlanar(const scalar& R,
                                     const scalar& alpha) const {
  return 2.0 * (0.5 * m_radii * m_radii * sin(alpha) * cos(alpha) +
                m_radii * sin(alpha) * R * cos(m_theta + alpha) +
                0.5 * R * R * sin(m_theta + alpha) * cos(m_theta + alpha)) +
         2.0 *
             (R * sin(m_theta + alpha) - R * sin(m_theta) +
              m_radii * sin(alpha)) *
             R * cos(m_theta) +
         R * R * sin(m_theta) * cos(m_theta) -
         (alpha * m_radii * m_radii + R * R * (M_PI - 2.0 * m_theta - alpha));
}

scalar CohesionTable::computeHPlanar(const scalar& R,
                                     const scalar& alpha) const {
  return computeH(R, alpha);
}

scalar CohesionTable::computeApproxAPlanar(const scalar& alpha,
                                           const scalar& d0) const {
  if (m_theta == 0.0)
    return 0.0;
  else {
    const scalar gamma = alpha + m_theta * 2.0;
    const scalar t2 = m_theta * 3.0;
    const scalar t3 = cos(t2);
    const scalar t4 = m_radii * m_radii;
    const scalar t5 = d0 * d0;
    const scalar t6 = cos(m_theta);
    const scalar t7 = M_PI * M_PI;
    const scalar t8 = m_theta * 5.0;
    const scalar t9 = cos(t8);
    const scalar t10 = gamma * gamma;
    const scalar t11 = sin(m_theta);
    const scalar t12 = sin(t2);
    const scalar t13 = sin(t8);
    return 1.0 / (t11 * t11 * t11) *
           (t3 * t4 * 6.0 - t3 * t5 * 1.2E1 + t5 * t6 * 1.2E1 - t4 * t9 * 6.0 -
            t4 * t11 * M_PI * 4.0 - t4 * t12 * M_PI * 1.6E1 -
            t5 * t11 * M_PI * 3.2E1 + t4 * t13 * M_PI * 4.0 -
            d0 * m_radii * t3 * 2.4E1 + d0 * m_radii * t6 * 2.4E1 -
            gamma * t4 * t11 * 3.2E1 + gamma * t4 * t12 * 2.8E1 +
            gamma * t5 * t11 * 3.2E1 - gamma * t4 * t13 * 4.0 -
            t3 * t4 * t7 * 3.0 - t3 * t4 * t10 * 3.0 + t4 * t6 * t7 * 2.2E1 +
            t5 * t6 * t7 * 2.0E1 + t4 * t6 * t10 * 2.2E1 + t4 * t7 * t9 +
            t5 * t6 * t10 * 2.0E1 + t4 * t9 * t10 + t4 * t11 * m_theta * 7.2E1 -
            t4 * t12 * m_theta * 2.4E1 + d0 * gamma * m_radii * t11 * 4.0E1 +
            d0 * gamma * m_radii * t12 * 8.0 + d0 * m_radii * t6 * t7 * 4.0E1 +
            d0 * m_radii * t6 * t10 * 4.0E1 -
            d0 * m_radii * t11 * M_PI * 4.0E1 -
            d0 * m_radii * t12 * M_PI * 8.0 + gamma * t3 * t4 * M_PI * 6.0 -
            gamma * t4 * t6 * M_PI * 4.4E1 - gamma * t5 * t6 * M_PI * 4.0E1 -
            gamma * t4 * t9 * M_PI * 2.0 -
            d0 * gamma * m_radii * t6 * M_PI * 8.0E1) *
           (1.0 / 4.8E1);
  }
}

scalar CohesionTable::computeApproxdEddPlanar(const scalar& alpha,
                                              const scalar& d0) const {
  const scalar gamma = alpha + m_theta * 2.0;

  if (m_theta == 0.0) {
    const scalar t2 = d0 * d0;
    const scalar t3 = m_radii * m_radii;
    return m_sigma * 1.0 / pow(d0 + m_radii, 2.0) *
           (t2 * M_PI * 4.0 + t3 * M_PI * 4.0 - gamma * t2 * 4.0 -
            gamma * t3 * 4.0 + d0 * m_radii * M_PI * 8.0 -
            d0 * gamma * m_radii * 8.0) *
           (1.0 / 2.0);
  } else {
    const scalar t2 = m_theta * 2.0;
    const scalar t3 = cos(t2);
    const scalar t4 = d0 * d0;
    const scalar t5 = m_radii * m_radii;
    const scalar t6 = t3 * t3;
    const scalar t7 = M_PI * M_PI;
    const scalar t8 = gamma * gamma;
    const scalar t9 = sin(t2);
    const scalar t10 = m_theta * 4.0;
    const scalar t11 = sin(t10);
    return (m_sigma * 1.0 / pow(d0 + m_radii * t3, 2.0) *
            (t4 * -2.0 + t3 * t4 * 2.0 + t4 * t7 - t5 * t6 * 2.0 + t4 * t8 -
             t5 * t7 * 2.0 - t5 * t8 * 2.0 - gamma * t4 * M_PI * 2.0 +
             gamma * t5 * M_PI * 4.0 - t4 * t9 * M_PI * 2.0 - t5 * t11 * M_PI -
             d0 * m_radii * t3 * 4.0 + d0 * m_radii * t6 * 4.0 -
             d0 * m_radii * t7 - d0 * m_radii * t8 + gamma * t4 * t9 * 2.0 +
             gamma * t5 * t11 - t3 * t4 * t7 + t3 * t5 * t6 * 2.0 -
             t3 * t4 * t8 + t3 * t5 * t7 + t3 * t5 * t8 + t5 * t6 * t7 +
             t5 * t6 * t8 + d0 * gamma * m_radii * t9 * 2.0 +
             d0 * gamma * m_radii * t11 + d0 * m_radii * t6 * t7 +
             d0 * m_radii * t6 * t8 + d0 * gamma * m_radii * M_PI * 2.0 -
             d0 * m_radii * t9 * M_PI * 2.0 - d0 * m_radii * t11 * M_PI +
             gamma * t3 * t4 * M_PI * 2.0 - gamma * t3 * t5 * M_PI * 2.0 -
             gamma * t5 * t6 * M_PI * 2.0 -
             d0 * gamma * m_radii * t6 * M_PI * 2.0) *
            (-1.0 / 2.0)) /
           sin(m_theta);
  }
}

scalar CohesionTable::computedEddAreaDist(const scalar& A_target,
                                          const scalar& d0) const {
  if (d0 < getDStar()) return computedEddAreaDist(A_target, getDStar());

  if (d0 <
      2.0 * sqrt(A_target / M_PI + 2.0 * m_radii * m_radii) - m_radii * 2.0)
    return 0.0;

  scalar alpha = interpolate_alpha(A_target, d0);
  scalar gamma = alpha + m_theta;

  scalar dEdd;

  if (gamma < PI / 2. + m_ang_epsilon && gamma > PI / 2. - m_ang_epsilon) {
    dEdd = computeApproxdEdd(alpha, d0);
  } else {
    scalar R_target = computeR(alpha, d0);
    dEdd = computedEdd(R_target, alpha);
  }

  return std::max(0.0, dEdd);
}

scalar CohesionTable::computedEddAreaDistPlanar(const scalar& A_target,
                                                const scalar& d0) const {
  if (d0 < getDStarPlanar())
    return computedEddAreaDistPlanar(A_target, getDStarPlanar());

  if (m_theta > 0.0 && d0 < sqrt((1.0 - cos(m_theta)) * (1.0 - cos(m_theta)) *
                                 (A_target + M_PI * m_radii * m_radii) /
                                 (m_theta - 0.5 * sin(2.0 * m_theta))) -
                                m_radii)
    return 0.0;

  scalar alpha = interpolate_alpha_planar(A_target, d0);

  scalar gamma = alpha + m_theta * 2.0;

  scalar dEdd;

  if (gamma < M_PI + m_ang_epsilon && gamma > M_PI - m_ang_epsilon) {
    dEdd = computeApproxdEddPlanar(alpha, d0);
  } else {
    scalar R_target = computeRPlanar(alpha, d0);
    dEdd = computedEddPlanar(R_target, alpha);
  }

  return std::max(0.0, dEdd);
}

void CohesionTable::construct_alpha_table() {
  scalar alpha_inc = m_max_alpha / (scalar)m_discretization;
  scalar d_inc = m_max_d0 / (scalar)m_discretization;

  m_alpha_table.resize(m_discretization, m_discretization);
  m_A_table.resize(m_discretization, m_discretization);
  m_dEdd_table.resize(m_discretization, m_discretization);
  m_max_As.resize(m_discretization);
  m_min_As.resize(m_discretization);

  const scalar dmin = getDMin();

  for (int i = 0; i < m_discretization; ++i) {
    scalar d0 = d_inc * (scalar)i + dmin;

    scalar maxA = -1e+20;
    scalar minA = 1e+20;
    for (int j = 0; j < m_discretization; ++j) {
      scalar alpha = alpha_inc * (scalar)j;

      scalar gamma = alpha + m_theta;

      scalar A;

      if (gamma < PI / 2. + m_ang_epsilon && gamma > PI / 2. - m_ang_epsilon) {
        A = computeApproxA(alpha, d0);
      } else {
        scalar R = computeR(alpha, d0);
        A = computeA(R, alpha);
      }

      m_A_table(j, i) = A;
      maxA = std::max(maxA, A);
      minA = std::min(minA, A);
    }

    m_max_As(i) = maxA;
    m_min_As(i) = minA;
  }

  // std::cout << m_A_table << std::endl;

  for (int i = 0; i < m_discretization; ++i) {
    const VectorXs& v = m_A_table.col(i);

    scalar A_inc = (m_max_As(i) - m_min_As(i)) / m_discretization;

    for (int j = 0; j < m_discretization; ++j) {
      scalar A_target = A_inc * (scalar)j + m_min_As(i);

      int retidx = std::min(mathutils::bipart_closest(v, A_target),
                            m_discretization - 1);
      int otheridx = std::min(m_discretization - 1, retidx + 1);

      scalar A0 = m_A_table(retidx, i);
      scalar A1 = m_A_table(otheridx, i);

      scalar alpha_target;

      if (A1 == A0)
        alpha_target = retidx * alpha_inc;
      else
        alpha_target =
            ((A_target - A0) / (A1 - A0) * (scalar)(otheridx - retidx) +
             retidx) *
            alpha_inc;

      if (alpha_target > m_max_alpha) {
        alpha_target = m_max_alpha;
      } else if (alpha_target < 0) {
        alpha_target = 0.0;
      }

      m_alpha_table(j, i) = alpha_target;
    }
  }

  for (int j = 0; j < m_discretization; ++j) {
    for (int i = 0; i < m_discretization; ++i) {
      scalar A_inc = (m_max_As(i) - m_min_As(i)) / m_discretization;
      scalar A_target = std::max(m_min_As(i), A_inc * (scalar)j + m_min_As(i));

      scalar d0 = d_inc * (scalar)i;
      m_dEdd_table(j, i) = computedEddAreaDist(A_target, d0);
    }
  }

  // std::cout << m_dEdd_table << std::endl;
  // std::cout << m_d2Edd2_table << std::endl;
}

void CohesionTable::construct_planar_alpha_table() {
  scalar alpha_inc = m_max_alpha / (scalar)m_discretization;
  scalar d_inc = m_max_d0 / (scalar)m_discretization;

  m_alpha_planar_table.resize(m_discretization, m_discretization);
  m_A_planar_table.resize(m_discretization, m_discretization);
  m_dEdd_planar_table.resize(m_discretization, m_discretization);
  m_max_A_planars.resize(m_discretization);
  m_min_A_planars.resize(m_discretization);

  const scalar dmin = getDMinPlanar();

  for (int i = 0; i < m_discretization; ++i) {
    scalar d0 = d_inc * (scalar)i + dmin;

    scalar maxA = -1e+20;
    scalar minA = 1e+20;
    for (int j = 0; j < m_discretization; ++j) {
      scalar alpha = alpha_inc * (scalar)j;

      scalar gamma = alpha + m_theta * 2.0;

      scalar A;

      if (gamma < M_PI + m_ang_epsilon && gamma > M_PI - m_ang_epsilon) {
        A = computeApproxAPlanar(alpha, d0);
      } else {
        scalar R = computeRPlanar(alpha, d0);
        A = computeAPlanar(R, alpha);
      }

      m_A_planar_table(j, i) = A;
      maxA = std::max(maxA, A);
      minA = std::min(minA, A);

      // std::cout << maxA << ", " << minA << ", " << alpha << ", " << j <<
      // std::endl;
    }

    m_max_A_planars(i) = maxA;
    m_min_A_planars(i) = minA;
  }

  // std::cout << m_A_planar_table << std::endl;

  for (int i = 0; i < m_discretization; ++i) {
    const VectorXs& v = m_A_planar_table.col(i);
    const size_t num_results = 1;

    scalar A_inc = (m_max_A_planars(i) - m_min_A_planars(i)) / m_discretization;

    for (int j = 0; j < m_discretization; ++j) {
      scalar A_target = A_inc * (scalar)j + m_min_A_planars(i);

      int retidx = std::min(mathutils::bipart_closest(v, A_target),
                            m_discretization - 1);
      int otheridx = std::min(m_discretization - 1, retidx + 1);

      scalar A0 = m_A_planar_table(retidx, i);
      scalar A1 = m_A_planar_table(otheridx, i);

      scalar alpha_target;

      if (A1 == A0)
        alpha_target = retidx * alpha_inc;
      else
        alpha_target =
            ((A_target - A0) / (A1 - A0) * (scalar)(otheridx - retidx) +
             retidx) *
            alpha_inc;

      if (alpha_target > m_max_alpha) {
        alpha_target = m_max_alpha;
      } else if (alpha_target < 0) {
        alpha_target = 0.0;
      }

      m_alpha_planar_table(j, i) = alpha_target;
    }
  }

  for (int j = 0; j < m_discretization; ++j) {
    for (int i = 0; i < m_discretization; ++i) {
      scalar A_inc =
          (m_max_A_planars(i) - m_min_A_planars(i)) / m_discretization;
      scalar A_target =
          std::max(m_min_A_planars(i), A_inc * (scalar)j + m_min_A_planars(i));

      scalar d0 = d_inc * (scalar)i;
      m_dEdd_planar_table(j, i) = computedEddAreaDistPlanar(A_target, d0);
    }
  }

  // std::cout << m_dEdd_table << std::endl;
  // std::cout << m_d2Edd2_table << std::endl;
}

scalar CohesionTable::getStiffness(const scalar& d0, const scalar& A_target,
                                   const scalar& pressure_weight) const {
  const scalar d_star = getDStar();
  const scalar d_hat = getDHat();
  const scalar dmin = getDMin();
  const scalar dry_criterion = 1e-12;

  if (d0 < d_hat) {
    if (d0 < dmin) {
      return m_collision_stiffness;
    } else {
      scalar dEdd_d_hat =
          A_target < dry_criterion
              ? 0.0
              : (interpolate_dEdd(A_target, d_hat) * pressure_weight);

      scalar stiffness_d_hat = dEdd_d_hat / (d_hat - d_star);

      if (d0 < d_star) {
        // linear interpolate
        scalar stiffness =
            ((stiffness_d_hat - m_collision_stiffness) * d0 +
             m_collision_stiffness * d_star - stiffness_d_hat * dmin) /
            (d_star - dmin);
        return stiffness;
      } else {
        // extrapolate the stiffness at d_hat
        return stiffness_d_hat;
      }
    }
  } else {
    scalar dEdd_d0 = A_target < dry_criterion
                         ? 0.0
                         : (interpolate_dEdd(A_target, d0) * pressure_weight);
    scalar stiffness_d0 = dEdd_d0 / (d0 - d_star);
    return stiffness_d0;
  }
}

scalar CohesionTable::getStiffnessPlanar(const scalar& d0,
                                         const scalar& A_target,
                                         const scalar& pressure_weight) const {
  const scalar d_star = getDStarPlanar();
  const scalar d_hat = getDHatPlanar();
  const scalar dmin = getDMinPlanar();
  const scalar dry_criterion = 1e-12;

  if (d0 < d_hat) {
    if (d0 < dmin) {
      return m_collision_stiffness_planar;
    } else {
      scalar dEdd_d_hat =
          A_target < dry_criterion
              ? 0.0
              : (interpolate_dEdd_planar(A_target, d_hat) * pressure_weight);

      scalar stiffness_d_hat = dEdd_d_hat / (d_hat - d_star);

      if (d0 < d_star) {
        // linear interpolate
        scalar stiffness =
            ((stiffness_d_hat - m_collision_stiffness_planar) * d0 +
             m_collision_stiffness_planar * d_star - stiffness_d_hat * dmin) /
            (d_star - dmin);
        return stiffness;
      } else {
        // extrapolate the stiffness at d_hat
        return stiffness_d_hat;
      }
    }
  } else {
    scalar dEdd_d0 =
        A_target < dry_criterion
            ? 0.0
            : (interpolate_dEdd_planar(A_target, d0) * pressure_weight);
    scalar stiffness_d0 = dEdd_d0 / (d0 - d_star);
    return stiffness_d0;
  }
}

scalar CohesionTable::getDMinPlanar() const { return m_min_d0_planar; }

scalar CohesionTable::getDHatPlanar() const {
  return getDMinPlanar() * (2.0 * m_radius_multiplier_planar - 1.0);
}

scalar CohesionTable::getDStarPlanar() const {
  return getDMinPlanar() * m_radius_multiplier_planar;
}

scalar CohesionTable::getDMin() const { return m_min_d0; }

scalar CohesionTable::getDHat() const {
  return getDMin() * (2.0 * m_radius_multiplier - 1.0);
}

scalar CohesionTable::getDStar() const {
  return getDMin() * m_radius_multiplier;
}

scalar CohesionTable::getRadiusMultiplierPlanar() const {
  return m_radius_multiplier_planar;
}

scalar CohesionTable::getCollisionStiffnessPlanar() const {
  return m_collision_stiffness_planar;
}

scalar CohesionTable::getRadiusMultiplier() const {
  return m_radius_multiplier;
}

scalar CohesionTable::getCollisionStiffness() const {
  return m_collision_stiffness;
}
