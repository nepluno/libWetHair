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

#ifndef LIBWETHAIR_CORE_TWO_D_SCENE_H_
#define LIBWETHAIR_CORE_TWO_D_SCENE_H_

#include <Eigen/Core>
#include <Eigen/StdVector>
#include <fstream>
#include <iostream>

#include "CohesionTableGen.h"
#include "DER/StrandParameters.h"
#include "Force.h"
#include "HairFlow.h"

namespace libwethair {

template <int DIM>
class PolygonalCohesion;

template <int DIM>
class HairFlow;

class FluidSim;

class StrandForce;

struct Script {
  enum TARGET {
    CAMERA,
    ROOT,
    SOLID,
    SOURCE,
    CURLRADIUS,
    CURLDENSITY,
    ALL,

    TARGET_COUNT
  };

  enum TYPE {
    ROTATE,
    TRANSLATE,
    SCALE,
    ABSORB,

    TYPE_COUNT
  };

  enum FUNC {
    CUBIC,
    COSINE,
    WENO,

    FUNC_COUNT
  };

  TARGET target;
  TYPE type;
  FUNC func;
  int index;
  Vector4s v;
  Vector3s origin;
  scalar start;
  scalar end;
  scalar ease_start;
  scalar ease_end;
  scalar amplitude;
  scalar frequency;
  scalar base_dt;
  scalar base_pos;
  bool updateSDF;
  bool transform_global;

  std::vector<scalar> base_vertices;
};

enum MASS_UPDATE_MODE {
  MUM_NONE,
  MUM_MASS_ONLY,
  MUM_DIRECT_DIV,
  MUM_MOMENTUM,

  MUM_COUNT
};

struct WetHairParameter {
  scalar dt;
  scalar rho;
  scalar sigma;
  scalar theta;
  scalar viscosity;
  scalar quadratic_dragging;
  scalar friction;
  scalar max_limited_eta_prop;
  scalar latitude;
  scalar earth_radius;
  scalar earth_rotation;
  scalar height_smooth;
  scalar regularizer_shell;
  scalar capillary_accel_multiplier;
  scalar dripping_radius_multiplier;
  scalar bulk_threshold_multiplier;
  scalar absorptionRate;
  scalar airviscosity;
  scalar drag_radius_multiplier;
  scalar hair_hair_cohesion_multiplier;
  scalar hair_solid_cohesion_multiplier;
  scalar radius_multiplier;
  scalar collision_stiffness;
  scalar radius_multiplier_planar;
  scalar collision_stiffness_planar;
  scalar damping_multiplier;
  scalar damping_multiplier_planar;
  scalar friction_multiplier_planar;

  int hairsteps;
  int swesteps;
  int fluidcorrectionsteps;

  bool no_fluids;
  bool no_swe;
  bool no_absorb;
  bool no_dripping;
  bool no_fictitious;

  bool drippingnear;
  bool drippingfar;
  bool drippingmiddle;
  bool use_ctcd;
  bool apply_coriolis;
  bool global_volume_control;
  bool individual_transfer;
  bool volume_summary;
  bool viscous_solve;

  MASS_UPDATE_MODE mass_update_mode;
  Vector3s gravity;

  // The max ratio allowed between the magnitude of the new and old
  // velocity.  Going over this will cause a runtime exception.
  scalar max_velocity_ratio;

  WetHairParameter();
};

template <int DIM>
class TwoDScene {
 public:
  TwoDScene(const bool& isMassSpring);

  TwoDScene(const TwoDScene& otherscene, const bool& isMassSpring);

  ~TwoDScene();

  int getNumParticles() const;
  int getNumEdges() const;
  int getNumHalfplanes() const;
  int getNumFlows() const;
  int getNumDofs() const;

  bool isMassSpring() const;

  bool useCtcd() const;

  const VectorXs& getX() const;

  VectorXs& getX();

  const VectorXs& getV() const;

  VectorXs& getV();

  const VectorXs& getM() const;

  VectorXs& getM();

  const VectorXs& getInterpolatedM() const;

  VectorXs& getInterpolatedM();

  const VectorXs& getLiquidRestLength() const;

  VectorXs& getLiquidRestLength();

  const scalar& getLiquidShell() const;

  scalar& getLiquidShell();

  const scalar getMaxVelocityRatio() const;

  const VectorXs& getHairRestMass() const;

  VectorXs& getHairRestMass();

  const VectorXs& getRadii() const;

  void resizeSystem(int num_particles, int num_edges, int num_strands = 0);

  void setPosition(int particle, const Vectors<DIM>& pos);
  Vectors<DIM> getPosition(int particle);

  void setScriptedGroup(int particle, int group_idx);
  int getScriptedGroup(int particle);

  void initializeScriptedGroup(const std::vector<Script>& scripts);

  int getDof(int particle) const;

  int getComponent(int dof) const;

  int getVertFromDof(int dof) const;

  void setVelocity(int particle, const Vectors<DIM>& vel);

  Vectors<DIM> getVelocity(int particle);

  void setMass(int particle, const scalar& mass);

  void setFixed(int particle, bool fixed);

  bool isFixed(int particle) const;

  bool isTip(int particle) const;

  const scalar& getRadius(int particle) const;
  void setRadius(int particle, scalar radius);

  void clearEdges();

  void clearHalfplanes();

  void setEdge(int idx, const std::pair<int, int>& edge, scalar radius = 0.055);

  void setEdgeRestLength(int idx, const scalar& l0);

  void setEdgePoissonRatio(int idx, const scalar& gamma);

  void insertHalfplane(const std::pair<VectorXs, VectorXs>& halfplane);

  const std::vector<std::pair<int, int>>& getEdges() const;

  // Each halfplane is a pair of two 2D vectors, a position and a normal.
  // The the Theme 02 Milestone 01 PDF for a description of these vectors.
  const std::vector<std::pair<VectorXs, VectorXs>>& getHalfplanes() const;

  const std::pair<VectorXs, VectorXs>& getHalfplane(int idx) const;

  const VectorXs& getEdgeRadii() const;

  const std::vector<int>& getParticleToHairs() const;

  const std::vector<int>& getParticleToHairLocalIndices() const;

  const std::pair<int, int>& getEdge(int edg) const;

  void insertForce(Force* newforce);

  void insertStrandParameters(StrandParameters* newparams);

  void insertStrandEquilibriumParameters(
      StrandEquilibriumParameters* newparams);

  void insertFilmFlow(HairFlow<DIM>* flow);

  void accumulateExternalGradU(VectorXs& F, const VectorXs& dx = VectorXs(),
                               const VectorXs& dv = VectorXs());

  void accumulateGradU(VectorXs& F, const VectorXs& dx = VectorXs(),
                       const VectorXs& dv = VectorXs());

  void accumulateExternalddUdxdx(TripletXs& A, const VectorXs& dx = VectorXs(),
                                 const VectorXs& dv = VectorXs());

  // Kind of a misnomer.
  void accumulateExternalddUdxdv(TripletXs& A, const VectorXs& dx = VectorXs(),
                                 const VectorXs& dv = VectorXs());

  void accumulateddUdxdx(TripletXs& A, const VectorXs& dx = VectorXs(),
                         const VectorXs& dv = VectorXs());

  // Kind of a misnomer.
  void accumulateddUdxdv(TripletXs& A, const VectorXs& dx = VectorXs(),
                         const VectorXs& dv = VectorXs());

  scalar computeKineticEnergy() const;
  scalar computePotentialEnergy() const;
  scalar computeTotalEnergy() const;

  Vectors<DIM> computeHairParticleMomentum() const;
  Vector3s computeHairParticleAngularMomentum() const;
  scalar computeHairParticleKineticEnergy() const;

  Vectors<DIM> computeLiquidParticleMomentum() const;
  Vector3s computeLiquidParticleAngularMomentum() const;
  scalar computeLiquidParticleKineticEnergy() const;

  Vectors<DIM> computeLiquidGridMomentum() const;
  Vector3s computeLiquidGridAngularMomentum() const;
  Vectors<DIM> computeReweightedLiquidGridMomentum() const;
  Vector3s computeReweightedLiquidGridAngularMomentum() const;
  scalar computeLiquidGridKineticEnergy() const;

  Vectors<DIM> computeHairGridMomentum() const;
  Vector3s computeHairGridAngularMomentum() const;
  Vectors<DIM> computeReweightedHairGridMomentum() const;
  Vector3s computeReweightedHairGridAngularMomentum() const;
  scalar computeHairGridKineticEnergy() const;

  scalar computeLiquidOverallDivergence() const;

  Vectors<DIM> computeCombinedGridMomentum() const;
  Vector3s computeCombinedGridAngularMomentum() const;
  Vectors<DIM> computeParticleWeightedCombinedGridMomentum() const;

  Vectors<DIM> computeHairWeightedCombinedGridMomentum() const;

  Vector3s computeParticleWeightedCombinedGridAngularMomentum() const;

  Vector3s computeHairWeightedCombinedGridAngularMomentum() const;

  scalar computeCombinedGridKineticEnergy() const;

  Vectors<DIM> computeHairDrag();

  void copyState(const TwoDScene& otherscene);

  void checkConsistency();

  std::vector<std::string>& getParticleTags();
  const std::vector<std::string>& getParticleTags() const;

  void gatherSimpleGravity();

  const scalar& getLatitude() const;

  const scalar& getEarthRadius() const;

  const scalar& getEarthRotation() const;

  const scalar& getHeightSmooth() const;

  bool applyCoriolis() const;

  const MASS_UPDATE_MODE getMassUpdateMode() const;

  void setLiquidParameter(const WetHairParameter& parameters);

  const std::vector<unsigned char>& getFixed() const;

  std::vector<unsigned char>& getFixed();

  const scalar& getMaxLimitEtaProp() const;

  const scalar& getViscosity() const;

  const scalar& getQuadraticDragging() const;

  const scalar& getHairFrictionCoeff() const;

  const scalar& getCapillaryAccelMultiplier() const;

  const scalar& getDrippingRadiusMultiplier() const;

  const scalar& getDragRadiusMultiplier() const;

  const scalar& getBulkThresholdMultiplier() const;

  const scalar& getHairHairCohesionMultiplier() const;

  const scalar& getHairSolidCohesionMultiplier() const;

  const scalar& getEpsilon() const;

  const scalar& getDt() const;

  int getHairSteps() const;

  int getSWESteps() const;

  const FluidSim* getFluidSim() const;

  FluidSim* getFluidSim();

  void setFluidSim(FluidSim*);

  void setPolygonalCohesion(PolygonalCohesion<DIM>*);

  const PolygonalCohesion<DIM>* getPolygonalCohesion() const;

  const scalar& getLiquidDensity() const;

  scalar getHairDensity(int pidx) const;

  const scalar& getLiquidTension() const;

  const scalar& getLiquidTheta() const;

  const Vector3s& getSimpleGravity() const;

  const std::vector<HairFlow<DIM>*>& getFilmFlows() const;

  std::vector<HairFlow<DIM>*>& getFilmFlows();

  void updateFlowGeometricState();

  void updateHairFlowsHeight(const VectorXs& accel, const scalar& dt);

  void advectHairFlows(const scalar& dt);

  void addForceHairFlows(const VectorXs& accel, const scalar& dt);

  void updateHairFlowsToGrid(const scalar& dt);

  void updateHairFlowsFromGrid(const scalar& dt);

  void updatePolygonalStructure(const scalar& dt);

  void updateHairMass();

  const scalar& getAbsorptionRate() const;

  const scalar& getAirViscosity() const;

  void uniformlyIncreaseHairLiquid(const scalar& dh, int group_idx);

  void viscousSolveFluidSim(const scalar& dt);

  void pressureSolveFluidSim(const scalar& dt);

  void advectRigidBodies(const scalar& dt);

  void advectFluidSim(const scalar& dt);

  void addForceFluidSim(const scalar& dt);

  void addGravityFluidSim(const scalar& dt);

  void updateBoundingBox();

  void initializeLiquids();

  void initializeCenters();

  void initializeMinimizers();

  void createMinimizers(const scalar& dt);

  void updateHairConnectivity();

  void saveCurrentX();

  scalar getDStarPlanar(const scalar& radius) const;

  scalar getStiffnessPlanar(const scalar& radius, const scalar& d0,
                            const scalar& A_target,
                            const scalar& pressure_weight) const;

  scalar getDStar(const scalar& radius) const;

  scalar getStiffness(const scalar& radius, const scalar& d0,
                      const scalar& A_target,
                      const scalar& pressure_weight) const;

  const scalar& getSearchRadius() const;

  const std::vector<std::vector<int>>& getParticleToEdge() const;

  void preCompute(const VectorXs& dx, const VectorXs& dv, const scalar& dt);

  void preComputeLocal(const VectorXs& dx, const VectorXs& dv,
                       const scalar& dt);

  void preComputeInterhair(const VectorXs& dx, const VectorXs& dv,
                           const scalar& dt);

  void postCompute(const scalar& dt);

  void postPreprocess(VectorXs& lambda, VectorXs& lambda_v, TripletXs& J,
                      TripletXs& Jv, TripletXs& Jxv, TripletXs& tildeK,
                      TripletXs& stiffness, TripletXs& damping, VectorXs& Phi,
                      VectorXs& Phiv, const VectorXs& dx, const VectorXs& dv,
                      const scalar& dt);

  void localPostPreprocess(VectorXs& lambda, VectorXs& lambda_v, TripletXs& J,
                           TripletXs& Jv, TripletXs& Jxv, TripletXs& tildeK,
                           TripletXs& stiffness, TripletXs& damping,
                           VectorXs& Phi, VectorXs& Phiv, const VectorXs& dx,
                           const VectorXs& dv, const scalar& dt);

  void interhairPostPreprocess(VectorXs& lambda, VectorXs& lambda_v,
                               TripletXs& J, TripletXs& Jv, TripletXs& Jxv,
                               TripletXs& tildeK, TripletXs& stiffness,
                               TripletXs& damping, VectorXs& Phi,
                               VectorXs& Phiv, const VectorXs& dx,
                               const VectorXs& dv, const scalar& dt);

  void getAffectedVars(int icol, std::unordered_set<int>& affected);

  void getAffectedForces(int icol, std::vector<Force*>& forces);

  void updateNumConstraints(int& num_constraint_pos, int& num_constraint_vel,
                            int& num_J, int& num_Jv, int& num_Jxv,
                            int& num_tildeK, Vector6i& interhair_param,
                            Vector6i& interhair_num);

  void updateNumConstraintsLocal(int& num_constraint_pos,
                                 int& num_constraint_vel, int& num_J,
                                 int& num_Jv, int& num_Jxv, int& num_tildeK);

  void updateNumConstraintsInterHair(int& num_constraint_pos,
                                     int& num_constraint_vel, int& num_J,
                                     int& num_Jv, int& num_Jxv, int& num_tildeK,
                                     Vector6i& interhair_param,
                                     Vector6i& interhair_num);

  void categorizeForces();

  void storeLambda(const VectorXs& lambda, const VectorXs& lambda_v);

  void computeLiquidPhi();

  void combineVelocityField();

  void interpolateMass();

  void updateParticleFlowsToGrid(bool with_hair_particles);

  void updateParticleFlowsFromGrid();

  void constrainHairParticles();

  void absorbReleaseParticleFlows(const scalar& dt, bool noabsorb,
                                  bool nodripping);

  void updateReservoir(const scalar& dt);

  void globalVolumeAdjust();

  void updateVolumeRecord();

  void reportParticleRemoved(const scalar& vol_removed);

  void reportParticleAdded(const scalar& vol_added);

  void notifyFullIntegrator();

  void notifyPartialIntegrator();

  void notifyGlobalIntegrator();

  void notifyLocalIntegrator();

  bool isIndividualTransfer() const;

  bool isFlowOnly() const;

  void updateEdgeRadii();

  void updateMassWithLiquid(int idx, const scalar& liquid_mass,
                            const scalar& liquid_radii);

  void saveFluidPressure(const std::string& szfn);

  const VectorXs& getFluidDragBuffer() const;

  std::vector<std::unordered_set<int>>& getBroadPhasePEPairs();

  std::vector<std::unordered_set<int>>& getBroadPhasePPPairs();

  std::vector<std::unordered_set<int>>& getBroadPhaseEEPairs();

  const std::vector<std::unordered_set<int>>& getBroadPhasePEPairs() const;

  const std::vector<std::unordered_set<int>>& getBroadPhasePPPairs() const;

  const std::vector<std::unordered_set<int>>& getBroadPhaseEEPairs() const;

  const Vectors<DIM>& getBoundingBoxMin() const;

  const Vectors<DIM>& getBoundingBoxMax() const;

  StrandParameters* getStrandParameters(const int index);

  StrandEquilibriumParameters* getStrandEquilibriumParameters(const int index);

  void updateStrandParamsTimestep(const scalar& dt);
  void updateStrandStartStates();

  void computeInterHairDifferentialVars();

  void setVertToDoFMap(const std::vector<int>& vert_to_dof,
                       const VectorXi& dofs_to_vars,
                       const std::vector<bool>& tipVerts,
                       const VectorXi& dof_to_vert);

  bool confirmNotNeighborEdgesOnSameHair(const int& efirst,
                                         const int& esecond) const;

  void computeMassesAndRadiiFromStrands();

  void addVolSummary();

  const std::vector<scalar>& getFreeVolSummary() const;

  const std::vector<scalar>& getHairVolSummary() const;

  bool doVolSummary() const;

  bool drippingNear() const;

  bool drippingFar() const;

  bool drippingMiddle() const;

  bool viscositySolve() const;

  void applyScript(const scalar& dt);

  Vector3s& getScriptedTranslate(int sg_idx);

  const Vector3s& getScriptedTranslate(int sg_idx) const;

  Eigen::Quaternion<scalar>& getScriptedRotation(int sg_idx);

  const Eigen::Quaternion<scalar>& getScriptedRotation(int sg_idx) const;

  void updateSearchRadius();

  const WetHairParameter& getParameter() const;

  WetHairParameter& getParameter();

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

 private:
  void resizeSystemDER(int num_particles, int num_edges, int num_strands);

  VectorXs m_base_x;
  VectorXs m_x;
  VectorXs m_v;
  VectorXs m_m;
  VectorXs m_rest_m;
  VectorXs m_radii;
  VectorXs m_interpolated_m;

  VectorXs m_fluid_drag_buffer;

  std::vector<int> m_particle_to_dofs;
  VectorXi m_dofs_to_component;
  std::vector<bool> m_is_strand_tip;
  VectorXi m_dofs_to_particle;
  std::vector<StrandForce*> m_strands;
  int m_num_strands;
  std::vector<int> m_edge_to_hair;

  std::vector<scalar> m_liquid_on_hair;
  std::vector<scalar> m_liquid_free;

  scalar m_search_radius;
  scalar m_collision_keep_proximity;

  Vectors<DIM> m_bb_min;
  Vectors<DIM> m_bb_max;

  FluidSim* m_fluid_sim;

  PolygonalCohesion<DIM>* m_polygonal_cohesion;

  std::vector<unsigned char> m_fixed;

  std::ofstream m_strout;
  // Vertex radii
  std::vector<std::pair<int, int>> m_edges;
  VectorXs m_edge_radii;
  VectorXs m_edge_rest_radii;
  VectorXs m_edge_rest_length;
  VectorXs m_edge_poisson_ratio;
  // Forces. Note that the scene inherits responsibility for deleting forces.
  std::vector<Force*> m_forces;
  std::vector<Force*> m_external_forces;
  std::vector<Force*> m_internal_forces;
  std::vector<std::vector<Force*>>
      m_hair_internal_forces;  // hair idx -> forces only affecting that hair
  std::vector<Force*> m_inter_hair_forces;
  // String 'tags' assigned to particles. Can be used to identify and single out
  // particles for special treatment.
  std::vector<std::string> m_particle_tags;

  std::vector<int> m_particle_to_hairs;

  std::vector<int> m_script_group;

  std::vector<Vector3s> m_scripted_translate;

  std::vector<Eigen::Quaternion<scalar>> m_scripted_rotation;

  std::vector<int> m_particle_to_hair_local_indices;

  std::vector<HairFlow<DIM>*> m_flows;

  std::vector<std::vector<int>> m_particle_to_edge;

  std::vector<std::unordered_set<int>> m_bp_edge_edge_pairs;
  std::vector<std::unordered_set<int>> m_bp_particle_edge_pairs;
  std::vector<std::unordered_set<int>> m_bp_particle_particle_pairs;

  std::vector<StrandParameters*> m_strandParameters;
  std::vector<StrandEquilibriumParameters*> m_strandEquilibriumParameters;

  const bool m_massSpringSim;

  scalar m_volume_hair;
  scalar m_volume_particle;
  scalar m_volume_particle_inserted;
  scalar m_volume_particle_removed;
  scalar m_volume_reservoir;

  scalar m_volume_hair_old;
  scalar m_volume_particle_old;
  scalar m_volume_reservoir_old;

  Vector6i m_constraint_idx;

  WetHairParameter m_parameters;
};

}  // namespace libwethair

#endif  // LIBWETHAIR_CORE_TWO_D_SCENE_H_
