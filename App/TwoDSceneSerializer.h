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

#ifndef LIBWETHAIR_APP_TWO_D_SCENE_SERIALIZER_H_
#define LIBWETHAIR_APP_TWO_D_SCENE_SERIALIZER_H_

#include <fstream>
#include <iostream>

#include "StringUtilities.h"

namespace libwethair {

template <int DIM>
class TwoDScene;

template <int DIM>
class SceneStepper;

}  // namespace libwethair

template <int DIM>
class TwoDSceneRenderer;

template <int DIM>
class TwoDSceneSerializer {
 public:
  void serializeFluidReadable(
      libwethair::TwoDScene<DIM>& scene,
      std::vector<std::ostringstream>& outputstream) const;

  void serializeHairReadable(libwethair::TwoDScene<DIM>& scene,
                             std::ostream& outputstream) const;

  void serializeShallowFlowReadable(const TwoDSceneRenderer<DIM>* renderer,
                                    libwethair::TwoDScene<DIM>& scene,
                                    std::ostream& outputstream) const;

  void serializeBoundariesReadable(const TwoDSceneRenderer<DIM>* renderer,
                                   libwethair::TwoDScene<DIM>& scene,
                                   std::ostream& os_boundary_single,
                                   std::ostream& os_boundary_double);

  void serializePolygonalCohesionReadable(
      const TwoDSceneRenderer<DIM>* renderer, libwethair::TwoDScene<DIM>& scene,
      std::ostream& os_pe, std::ostream& os_poe, std::ostream& os_ppp) const;

  bool deSerializeFluidReadable(
      libwethair::TwoDScene<DIM>& scene,
      const std::vector<std::string>& filename_fluids);

  bool deSerializeHairReadable(libwethair::TwoDScene<DIM>& scene,
                               const std::string& filename_hairs);

  bool deSerializeShallowFlowReadable(const TwoDSceneRenderer<DIM>* renderer,
                                      libwethair::TwoDScene<DIM>& scene,
                                      const std::string& filename_flows);

  void serializeScene(libwethair::TwoDScene<DIM>& scene,
                      libwethair::SceneStepper<DIM>* stepper,
                      std::ostream& outputstream) const;

  void loadScene(libwethair::TwoDScene<DIM>& scene,
                 libwethair::SceneStepper<DIM>* stepper,
                 std::istream& inputstream) const;
};

#endif  // LIBWETHAIR_APP_TWO_D_SCENE_SERIALIZER_H_
