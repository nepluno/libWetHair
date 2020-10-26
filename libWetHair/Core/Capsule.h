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

#ifndef LIBWETHAIR_CORE_CAPSULE_H_
#define LIBWETHAIR_CORE_CAPSULE_H_

/* Capsule generation algorithm.
 * Adapted from Paul Bourke's C implementation found here:
 * http://paulbourke.net/geometry/capsule/
 */

#include <unordered_map>
#include <vector>

#include "MathDefs.h"

class CapsuleCreator {
 public:
  std::vector<Vector3s> vertices;
  std::vector<Vector3i> indices;

  CapsuleCreator();
  void Create(int N, const scalar& radius, const scalar& halfheight);
};

#endif  // LIBWETHAIR_CORE_CAPSULE_H_
