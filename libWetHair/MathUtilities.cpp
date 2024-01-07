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

#include "MathUtilities.h"

namespace libwethair {
namespace mathutils {

bool approxSymmetric(const MatrixXs& A, const scalar& eps) {
  for (int i = 0; i < A.rows(); ++i)
    for (int j = i + 1; j < A.cols(); ++j)
      if (fabs(A(i, j) - A(j, i)) >= eps)
        return false;
  return true;
}

}  // namespace mathutils

// explicit instantiations at bottom
template struct int_Vectors<2>;
template struct int_Vectors<3>;

}  // namespace libwethair
