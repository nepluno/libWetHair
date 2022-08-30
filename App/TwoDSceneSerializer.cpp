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

#include "TwoDSceneSerializer.h"

#include <fstream>
#include <iostream>
#include <string>

#include <libWetHair/HairFlow.h>
#include <libWetHair/PolygonalCohesion.h>
#include <libWetHair/SceneStepper.h>
#include <libWetHair/TwoDScene.h>
#include <libWetHair/fluidsim2D.h>
#include <libWetHair/fluidsim3D.h>

#include "TwoDSceneRenderer.h"

#define LOAD_HAIR_ONLY

using namespace libwethair;

template <int DIM>
void TwoDSceneSerializer<DIM>::serializeFluidReadable(
    TwoDScene<DIM>& scene,
    std::vector<std::ostringstream>& outputstreams) const {
  const FluidSim* sim = scene.getFluidSim();
  sim->writeReadable(outputstreams);
}

template <int DIM>
void TwoDSceneSerializer<DIM>::serializeBoundariesReadable(
    const TwoDSceneRenderer<DIM>* renderer,
    TwoDScene<DIM>& scene,
    std::ostream& os_boundary_single,
    std::ostream& os_boundary_double) {
  renderer->writeBoundaries(os_boundary_single, os_boundary_double);
}

template <int DIM>
void TwoDSceneSerializer<DIM>::serializeHairReadable(
    TwoDScene<DIM>& scene, std::ostream& outputstream) const {
  const std::vector<HairFlow<DIM>*> flows = scene.getFilmFlows();
  for (const HairFlow<DIM>* flow : flows) {
    flow->writeReadable(outputstream);
  }
}

template <int DIM>
void TwoDSceneSerializer<DIM>::serializePolygonalCohesionReadable(
    const TwoDSceneRenderer<DIM>* renderer,
    TwoDScene<DIM>& scene,
    std::ostream& os_pe, std::ostream& os_poe, std::ostream& os_ppp) const {
  const PolygonalCohesion<DIM>* cohesion = scene.getPolygonalCohesion();
  cohesion->writeReadable(os_pe, os_poe, os_ppp);
}

template <int DIM>
void TwoDSceneSerializer<DIM>::serializeShallowFlowReadable(
    const TwoDSceneRenderer<DIM>* renderer,
    TwoDScene<DIM>& scene,
    std::ostream& outputstream) const {
  renderer->writeTransformedHairFlow(outputstream, scene);
}

//////////////////

template <int DIM>
bool TwoDSceneSerializer<DIM>::deSerializeFluidReadable(
    TwoDScene<DIM>& scene, const std::vector<std::string>& filename_fluids) {
  bool f = false;
  std::ifstream file(filename_fluids[0].c_str());
  if (file.is_open()) {
    f = true;
    FluidSim* sim = scene.getFluidSim();
    sim->readReadable(file);
  }
  file.close();

  // sim->readReadableBoundary(os_boundary_single, os_boundary_double); // for
  // now
  return f;
}

template <int DIM>
bool TwoDSceneSerializer<DIM>::deSerializeHairReadable(
    TwoDScene<DIM>& scene, const std::string& filename_hairs) {
  bool h = false;
  std::ifstream file(filename_hairs.c_str());
  if (file.is_open()) {
    h = true;
    std::vector<HairFlow<DIM>*> flows = scene.getFilmFlows();
    for (HairFlow<DIM>* flow : flows) {
      flow->readReadable(file);
      flow->updateGeometricState(scene.getX(), scene.getV(),
                                 scene.getFluidSim());
    }
  }
  file.close();
  return h;
}

template <int DIM>
bool TwoDSceneSerializer<DIM>::deSerializeShallowFlowReadable(
    const TwoDSceneRenderer<DIM>* renderer, TwoDScene<DIM>& scene,
    const std::string& filename_flows) {
  // renderer->writeTransformedHairFlow(outputstream, scene); // for now
  return false;
}

//////////////////

template <int DIM>
void TwoDSceneSerializer<DIM>::serializeScene(
    TwoDScene<DIM>& scene, SceneStepper<DIM>* stepper,
    std::ostream& outputstream) const {
  int ndof = scene.getNumDofs();
  scalar* xdata = scene.getX().data();
  outputstream.write((char*)xdata, ndof * sizeof(scalar));
  scalar* vdata = scene.getV().data();
  outputstream.write((char*)vdata, ndof * sizeof(scalar));

  size_t flowdatasize = 0;
  size_t fluiddatasize = 0;
  size_t fluidcrucialgridsize = 0;
  size_t boundarydatasize = 0;
  size_t integratorsize = 0;

  std::vector<scalar> buf;

  const std::vector<HairFlow<DIM>*> flows = scene.getFilmFlows();
  for (const HairFlow<DIM>* flow : flows) {
    flow->write(buf);
    flowdatasize += flow->serialized_size();
  }

  if (DIM == 2) {
    const FluidSim2D* fluid2d = (const FluidSim2D*)scene.getFluidSim();
    fluid2d->write(buf);
    fluiddatasize += fluid2d->particle_size();
    boundarydatasize += fluid2d->boundary_size();
    fluidcrucialgridsize += fluid2d->crucial_grid_size();
  } else {
    const FluidSim3D* fluid3d = (const FluidSim3D*)scene.getFluidSim();
    fluid3d->write(buf);
    fluiddatasize += fluid3d->particle_size();
    boundarydatasize += fluid3d->boundary_size();
    fluidcrucialgridsize += fluid3d->crucial_grid_size();
  }

  stepper->write(buf);
  integratorsize += stepper->size();

  outputstream.write((char*)&flowdatasize, sizeof(size_t));
  outputstream.write((char*)&fluiddatasize, sizeof(size_t));
  outputstream.write((char*)&boundarydatasize, sizeof(size_t));
  outputstream.write((char*)&fluidcrucialgridsize, sizeof(size_t));
  outputstream.write((char*)&integratorsize, sizeof(size_t));

  outputstream.write((char*)&buf[0], buf.size() * sizeof(scalar));
  outputstream.flush();
}

template <int DIM>
void TwoDSceneSerializer<DIM>::loadScene(TwoDScene<DIM>& scene,
                                         SceneStepper<DIM>* stepper,
                                         std::istream& inputstream) const {
  int ndof = scene.getNumDofs();
  scalar* xdata = scene.getX().data();
  inputstream.read((char*)xdata, ndof * sizeof(scalar));
#ifndef LOAD_HAIR_ONLY
  scalar* vdata = scene.getV().data();
  inputstream.read((char*)vdata, ndof * sizeof(scalar));

  size_t flowdatasize = 0;
  size_t fluiddatasize = 0;
  size_t boundarydatasize = 0;
  size_t integratorsize = 0;
  size_t fluidcrucialgridsize = 0;

  inputstream.read((char*)&flowdatasize, sizeof(size_t));
  inputstream.read((char*)&fluiddatasize, sizeof(size_t));
  inputstream.read((char*)&boundarydatasize, sizeof(size_t));
  inputstream.read((char*)&fluidcrucialgridsize, sizeof(size_t));
  inputstream.read((char*)&integratorsize, sizeof(size_t));

  std::vector<scalar> buf;
  size_t datasize = flowdatasize + fluiddatasize + boundarydatasize +
                    fluidcrucialgridsize + integratorsize;

  buf.resize(datasize / sizeof(scalar));
  inputstream.read((char*)&buf[0], datasize);

  if (inputstream.fail()) {
    std::cout << outputmod::startred
              << "Error in TwoDSceneSerialized: " << outputmod::endred
              << "Failed to load timestep. Exiting." << std::endl;
    exit(1);
  }

  int k = 0;
  const std::vector<HairFlow<DIM>*> flows = scene.getFilmFlows();
  for (HairFlow<DIM>* flow : flows) {
    const char* data_ptr = ((const char*)&buf[0]) + k;
    flow->read((const scalar*)data_ptr);
    k += flow->serialized_size();
  }

  if (DIM == 2) {
    FluidSim2D* fluid2d = (FluidSim2D*)scene.getFluidSim();
    const char* data_ptr = ((const char*)&buf[0]) + k;
    fluid2d->read((const scalar*)data_ptr, fluiddatasize, boundarydatasize);
    k += (fluiddatasize + boundarydatasize + fluidcrucialgridsize);
  } else {
    FluidSim3D* fluid3d = (FluidSim3D*)scene.getFluidSim();
    const char* data_ptr = ((const char*)&buf[0]) + k;
    fluid3d->read((const scalar*)data_ptr, fluiddatasize, boundarydatasize);
    k += (fluiddatasize + boundarydatasize + fluidcrucialgridsize);
  }

  {
    const char* data_ptr = ((const char*)&buf[0]) + k;
    stepper->read((const scalar*)data_ptr);
  }
#endif
}

// explicit instantiations at bottom
template class TwoDSceneSerializer<2>;
template class TwoDSceneSerializer<3>;
