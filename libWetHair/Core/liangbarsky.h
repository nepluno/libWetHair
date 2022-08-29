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

#ifndef LIBWETHAIR_CORE_LIANG_BARSKY_H_
#define LIBWETHAIR_CORE_LIANG_BARSKY_H_

/* Liang-Barsky clipping algorithm.
 */

#include "MathDefs.h"

namespace libwethair {
namespace liangbarsky {
/* clip_line()
 * modifies parameters in place to clip the line,
 * returns 0 if line is totally outside clip window
 * returns 1 if line is not totally outside clip window
 */
int clip_line(const Vector4s& c, Vector2s& q1, Vector2s& q2, scalar& alpha0,
              scalar& alpha1);
int clip_line(const Vector6s& c, Vector3s& q1, Vector3s& q2, scalar& t0,
              scalar& t1);
}  // namespace liangbarsky
}  // namespace libwethair

#endif  // LIBWETHAIR_CORE_LIANG_BARSKY_H_
