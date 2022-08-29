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

#ifndef LIBWETHAIR_CORE_ICOSPHERE_H_
#define LIBWETHAIR_CORE_ICOSPHERE_H_

/* Icosphere generation algorithm.
 * Adapted from Andreas Kahler's C# implementation found here:
 * http://blog.andreaskahler.com/2009/06/creating-icosphere-mesh-in-code.html
 */

#include <unordered_map>
#include <vector>

#include "MathDefs.h"

namespace libwethair {

class IcosphereCreator {
 public:
  std::vector<Vector3s> vertices;
  std::vector<Vector3i> indices;
  std::unordered_map<uint64, int> middlePointIndexCache;
  int index;

  int addVertex(const Vector3s& p);
  int getMiddlePoint(int p1, int p2);

  IcosphereCreator();
  void Create(int recursionLevel);
};

}  // namespace libwethair

#endif  // LIBWETHAIR_CORE_ICOSPHERE_H_
