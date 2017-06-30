//
// This file is part of the libWetHair open source project
//
// The code is licensed solely for academic and non-commercial use under the
// terms of the Clear BSD License. The terms of the Clear BSD License are
// provided below. Other licenses may be obtained by contacting the faculty
// of the Columbia Computer Graphics Group or a Columbia University licensing officer.
//
// The Clear BSD License
//
// Copyright 2017 Yun (Raymond) Fei, Henrique Teles Maia, Christopher Batty,
// Changxi Zheng, and Eitan Grinspun
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted (subject to the limitations in the disclaimer
// below) provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice, this
//  list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright notice,
//  this list of conditions and the following disclaimer in the documentation
//  and/or other materials provided with the distribution.
//
// * Neither the name of the copyright holder nor the names of its contributors may be used
//  to endorse or promote products derived from this software without specific
//  prior written permission.
//
// NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE GRANTED BY THIS
// LICENSE. THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
// GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
// OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.


#ifndef __HAIR_FLOW_H__
#define __HAIR_FLOW_H__

#include <Eigen/Core>
#include <iostream>
#include <cassert>
#include <vector>
#include <unordered_map>

#include "Force.h"

#include "MathDefs.h"
#include "MathUtilities.h"

template<int DIM>
class TwoDScene;

class FluidSim;

enum LIQUID_SIM_TYPE
{
  LST_SHALLOW,
  
  LST_COUNT
};

template<int DIM>
struct HairParticleBridge
{
  Vectors<DIM> vel;
  scalar volume;
  scalar alpha;
  int pidx;
  int eidx;
};

struct HairPlaneIntersection
{
  int lidx;
  int pidx_0;
  int pidx_1;
  int flow_idx;
  int hair_idx;
  int prev_inter_idx;
  int next_inter_idx;
  scalar dist_to_prev;
  scalar dist_to_next;
  scalar flow_height;
  scalar alpha;
  scalar pos_y;
  scalar porosity;
  Vector2s jac;
  Matrix2s hess;
};

template<int DIM>
class HairFlow
{
protected:
  static int m_flow_counter;
  
public:
  
  enum STATE
  {
    NORMAL,
    ENDPOINT,
    STAR
  };
  
  HairFlow(TwoDScene<DIM>* parent, const std::vector<int>& involved_particles, const VectorXs& eta, const std::vector<unsigned char>& particle_state);
  
  virtual const scalar& getAvgAreaE() const;
  virtual const scalar& getMinAreaE() const;
  virtual const scalar& getMaxAreaE() const;
  virtual const scalar& getMinAreaV() const;
  virtual const scalar& getMaxAreaV() const;
  virtual const scalar& getMinEta() const;
  virtual const scalar& getMaxEta() const;
  
  virtual const VectorXs& getEta() const;
  
  virtual VectorXs& getEta();
  
  virtual const VectorXs& getAvgEta() const;
  
  virtual VectorXs& getAvgEta();
  
  virtual const VectorXs& getVelocity() const;
  
  virtual VectorXs& getVelocity();
  
  virtual const MatrixXs& getActualVertexVelocity() const;
  
  virtual MatrixXs& getActualVertexVelocity();
  
  virtual const VectorXs& getPorosity() const;
  
  virtual VectorXs& getPorosity();
  
  virtual const VectorXs& getAreaV() const;
  
  virtual VectorXs& getAreaV();
  
  virtual const VectorXs& getAreaVHair() const;
  
  virtual VectorXs& getAreaVHair();
  
  virtual const MatrixXs& getAccelV() const;
  
  virtual MatrixXs& getAccelV();
  
  virtual const VectorXs& getAreaE() const;
  
  virtual VectorXs& getAreaE();
  
  virtual const VectorXs& getNormalV() const;
  
  virtual VectorXs& getNormalV();
  
  virtual const VectorXs& getRadiiV() const;
  
  virtual VectorXs& getRadiiV();
  
  virtual const VectorXs& getRadiiE() const;
  
  virtual VectorXs& getRadiiE();
  
  virtual const VectorXs& getNormalE() const;
  
  virtual VectorXs& getNormalE();
  
  virtual const MatrixXs& getTangentV() const;
  
  virtual MatrixXs& getTangentV();
  
  virtual const MatrixXs& getTangentE() const;
  
  virtual MatrixXs& getTangentE();
  
  virtual scalar getTotalLength() const = 0;
  
  virtual void adjustVolumeGlobal(const scalar& prop);
  
  virtual std::unordered_map<int, int>& getGlobalToLocal();
  
  virtual const std::unordered_map<int, int>& getGlobalToLocal() const;
  
  virtual const std::vector< std::vector<int> >& getEdgeBridges() const;
  
  virtual std::vector< std::vector<int> >& getEdgeBridges();
  
  virtual const std::vector< std::pair<int, int> >& getLocalEdges() const;
  
  virtual const std::vector< std::pair<int, int> >& getGlobalEdges() const;
  
  virtual const std::vector< int >& getEdgeIndices() const;
  
  virtual const std::vector< int >& getParticleIndices() const;
  
  virtual void updateGeometricState(const VectorXs& x, const VectorXs& v, FluidSim* fluidsim) = 0;
  
  virtual void updateToFilteredGrid(const VectorXs& x, const VectorXs& v, FluidSim* fluidsim, const scalar& dt, int ibuffer) = 0;
  
  virtual void updateFromFilteredGrid(const VectorXs& x, VectorXs& v, FluidSim* fluidsim, const scalar& dt) = 0;
  
  virtual void advance(const VectorXs& x, const scalar& dt) = 0;
  
  virtual void add_force(const VectorXs& x, const VectorXs& accel, FluidSim* fluidsim, const scalar& dt) = 0;
  
  virtual void updateHairMass() = 0;
  
  virtual void updateHairFlowHeight(const scalar& dt) = 0;
  
  virtual void preUpdateHairFlowHeight(const scalar& dt) = 0;
  
  virtual void postUpdateHairFlowHeight(const scalar& dt) = 0;
  
  virtual void updateReservoir(FluidSim* fluidsim, const VectorXs& x, const VectorXs& v, const scalar& dt) = 0;
  
  virtual const std::vector<unsigned char>& getState() const;
  
  virtual void resizeSystem() = 0;
  
  virtual int index() const;
  
  virtual int find(int idx_global) const;
  
  virtual int size() const;
  
  virtual void getAffectedVars( int pidx, std::unordered_set<int>& vars );
  
  virtual Vectors<DIM> computeHairLiquidMomentum(const VectorXs& v) const;
  
  virtual Vector3s computeHairLiquidAngularMomentum(const VectorXs& x, const VectorXs& v, FluidSim* fluidsim) const = 0;
  
  virtual Vectors<DIM> computeHairDragForce(const VectorXs& v) const;
  
  virtual scalar computeHairLiquidEnergy(const VectorXs& v) const;
  
  virtual bool isContained( int pidx );
  
  virtual scalar computeTotalLiquidVol() const = 0;
  
  virtual scalar computeTotalReservoirVol() const = 0;
  
  virtual const scalar& getPoolSize() const = 0;
  
  virtual scalar& getPoolSize() = 0;
  
  virtual void read(const scalar* data);
  
  virtual void write(std::vector<scalar>& data) const;
  
  virtual void writeReadable(std::ostream&) const;
  virtual void readReadable(std::istream& i);
  
  virtual size_t serialized_size() const;

  virtual void setConstraintParameters(const Vector6i& start, const Vector6i& num);
  
  virtual const Vector6i& getConstraintIdx() const;
  
  virtual const Vector6i& getNumConstraints() const;
  
  virtual void geodesic_to_local(const scalar& geopos, int& pidx_low, scalar& alpha) const = 0;
  
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
protected:
  TwoDScene<DIM>* m_parent;
  
  std::vector<int> m_particle_indices;
  std::vector<unsigned char> m_particle_state;
  std::vector<int> m_edge_indices;
  std::vector< std::pair<int, int> > m_internal_edges;
  std::vector< std::pair<int, int> > m_global_edges;
  std::vector< std::vector<int> > m_particle_to_edges;
  std::unordered_map<int, int> m_global_to_local;
  
  // liquid state
  VectorXs m_liquid_phi;
  VectorXs m_porosity;
  VectorXs m_eta;
  VectorXs m_avg_eta;
  VectorXs m_edge_eta;
  VectorXs m_u;
  VectorXs m_stored_friction_coeff;
  
  // geometric cache
  MatrixXs m_accel_v;
  VectorXs m_area_v;
  VectorXs m_area_v_hair;
  VectorXs m_area_e;
  VectorXs m_normal_e;
  VectorXs m_normal_v;
  VectorXs m_rad_vec;
  VectorXs m_edge_rad_vec;
  MatrixXs m_dir_v;
  MatrixXs m_dir_f;
  
  MatrixXs m_actual_u_f;
  MatrixXs m_actual_u_v;
  
  // apic
  MatrixXs m_c_v;
  MatrixXs m_c_star;
  
  scalar m_avg_area_e;
  scalar m_min_area_e;
  scalar m_max_area_e;
  
  scalar m_min_area_v;
  scalar m_max_area_v;
  scalar m_max_eta;
  scalar m_min_eta;
  
  Vector6i m_constraint_starts;
  Vector6i m_constraint_length;
  
  std::vector< std::vector<int> > m_edge_bridges;
  
  int m_flow_index;
};

#endif
