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

#ifndef LIBWETHAIR_CORE_DER_FORCE_FORCE_ACCUMULATOR_H_
#define LIBWETHAIR_CORE_DER_FORCE_FORCE_ACCUMULATOR_H_

#include <libWetHair/DER/Definitions.h>
#include <libWetHair/DER/StrandForce.h>

namespace libwethair {

template <typename ForceT>
class ForceAccumulator {
 public:
  template <typename AccumulatedT>
  static void accumulateCurrent(AccumulatedT& accumulated,
                                StrandForce& strand) {
    // Yes, this statement does nothing... It is there to guarantee that
    // accumulateCurrent will never be called on dissipative forces, as they
    // require evolution to be estimated. So if you get a compilation error
    // here, that's why. If thas is a problem just comment out the line below,
    // but make sure that dissipative forces still make sense when called on the
    // current state (most likely they should return zero).
    ForceT::NonDissipativeForce;

    accumulate(accumulated, strand);
  }

  template <typename AccumulatedT>
  static void accumulateFuture(AccumulatedT& accumulated, StrandForce& strand) {
    accumulate(accumulated, strand);
  }

  static void accumulate(scalar& energy, const StrandForce& strand) {
    for (IndexType vtx = ForceT::s_first;
         vtx < strand.getNumVertices() - ForceT::s_last; ++vtx) {
      energy += ForceT::localEnergy(strand, vtx);
    }
  }

  static void accumulate(VecX& force, const StrandForce& strand) {
    typename ForceT::LocalForceType localF;
    for (IndexType vtx = ForceT::s_first;
         vtx < strand.getNumVertices() - ForceT::s_last; ++vtx) {
      ForceT::computeLocal(localF, strand, vtx);
      ForceT::addInPosition(force, vtx, localF);
    }
  }

  // Jacobian of the Force <==>  - Hessian of the Energy
  static void accumulate(TripletXs& hessianOfEnergy,
                         const StrandForce& strand) {
    typename ForceT::LocalJacobianType localJ;
    for (IndexType vtx = ForceT::s_first;
         vtx < strand.getNumVertices() - ForceT::s_last; ++vtx) {
      ForceT::computeLocal(localJ, strand, vtx);

      if (localJ.rows() > 6) {  // (Bending & Twisting)
        for (IndexType r = 0; r < localJ.rows(); ++r) {
          for (IndexType c = 0; c < localJ.cols(); ++c) {
            if (isSmall(localJ(r, c))) continue;
            hessianOfEnergy.push_back(
                Triplets((vtx - 1) * 4 + r, (vtx - 1) * 4 + c, localJ(r, c)));
          }
        }
      } else {  // Stretch
        int trCount = 0;
        for (IndexType r = 0; r < localJ.rows(); ++r) {
          if (r == 3) {  // skip twist dof
            ++trCount;
          }
          int tcCount = 0;
          for (IndexType c = 0; c < localJ.cols(); ++c) {
            if (c == 3) {  // skip twist dof
              ++tcCount;
            }
            if (isSmall(localJ(r, c))) continue;
            hessianOfEnergy.push_back(Triplets(
                vtx * 4 + r + trCount, vtx * 4 + c + tcCount, localJ(r, c)));
          }
        }
      }
    }
  }
};

}  // namespace libwethair

#endif  // LIBWETHAIR_CORE_DER_FORCE_FORCE_ACCUMULATOR_H_
