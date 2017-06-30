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

#ifndef __POLYGONAL_COHESION_H__
#define __POLYGONAL_COHESION_H__

#include <Eigen/Core>

#include "MathDefs.h"
#include "CohesionTableGen.h"
#include "HairFlow.h"
#include "Force.h"
#include <unordered_set>
#include <set>
#include <unordered_map>
#include <fstream>
#include <sstream>
#include <deque>
#include <nanoflann/nanoflann.hpp>
#include <thread>
#include "sorter.h"

template<int DIM>
class TwoDScene;

struct PolygonalRegion
{
  int hair_idx_0;
  int hair_idx_1;
  std::vector<int> local_indices_0;
  std::vector<int> local_indices_1;
};

struct EdgeEdgePair
{
  int base_eidx;
  int neighbor_hair_idx;
  scalar alpha_0;
  scalar alpha_1;
  scalar neighbor_local_coord_0;
  scalar neighbor_local_coord_1;
  scalar count_0;
  scalar count_1;
  
  std::string toString() const
  {
    std::stringstream ss;
    ss << "base: " << base_eidx << " neighborhair: " << neighbor_hair_idx << " alpha: (" << alpha_0 << ", " << alpha_1 << ") ncoord: (" << neighbor_local_coord_0 << ", " << neighbor_local_coord_1 << ") count: (" << count_0 << ", " << count_1 << ")";
    return ss.str();
  }
};

struct EdgeEdgePairEEC
{
  int base_eidx;
  int neighbor_eidx;
  scalar alpha_0;
  scalar alpha_1;
  scalar neighbor_local_coord_0;
  scalar neighbor_local_coord_1;

  scalar time;
  scalar alpha_contact;
  scalar neighbor_local_coord_contact;
  Vector3s avgpos;
  
  scalar count_0; //outdated
  scalar count_1; //outdated

  bool valid;
  bool updated;
  
  std::string toString() const
  {
    std::stringstream ss;
    ss << "base: " << base_eidx << " neighboredge: " << neighbor_eidx << " alpha: (" << alpha_0 << ", " << alpha_1 << ") ncoord: (" << neighbor_local_coord_0 << ", " << neighbor_local_coord_1 << ")";// count: (" << count_0 << ", " << count_1 << ")";
    return ss.str();
  }
};

struct ParticleEdgePair
{
  int pidx;
  int eidx;
  scalar alpha;
  scalar dist;
  scalar radii_j;
  scalar max_dist;
  scalar count;
  bool valid;
  bool updated;
  bool should_be_deleted;
  bool latest;
  
  std::string toString() const
  {
    std::stringstream ss;
    ss << "pidx: " << pidx << " eidx: " << eidx << " alpha: " << alpha << " dist: " << dist << " edgeradii: " << radii_j << " maxd: " << max_dist << " count: " << count << " valid: " << (valid ? "true" : "false") << " updated: " << (updated ? "true" : "false") << ")";
    return ss.str();
  }
  
  static size_t size()
  {
    return ((sizeof(ParticleEdgePair) + sizeof(scalar) - 1) / sizeof(scalar)) * sizeof(scalar);
  }
  
  void write(std::vector<scalar>& data) const
  {
    size_t N = size() / sizeof(scalar);
    scalar* buf = new scalar[N];
    memset(buf, 0, N * sizeof(scalar));
    memcpy(buf, this, sizeof(ParticleEdgePair));
    for(int i = 0; i < N; ++i)
    {
      data.push_back(buf[i]);
    }
    
    delete[] buf;
  }
  
  void read(const scalar* data)
  {
    const ParticleEdgePair* b = (ParticleEdgePair*) data;
    pidx = b->pidx;
    eidx = b->eidx;
    alpha = b->alpha;
    dist = b->dist;
    radii_j = b->radii_j;
    max_dist = b->max_dist;
    count = b->count;
    valid = b->valid;
    updated = b->updated;
    latest = b->latest;
  }
};

struct PointEdgePair
{
  int base_eidx;
  scalar alpha_point;
  int neighbor_eidx;
  scalar neighbor_alpha_point;
  scalar V;
  scalar quadrature_weight;
  scalar pressure_weight;
  uint64 hash_code;
  
  scalar time;
  std::string toString() const
  {
    std::stringstream ss;
    ss << "base: (" << base_eidx << ", " << alpha_point << ") neighbor: (" << neighbor_eidx << ", " << neighbor_alpha_point << ") V: " << V << " q-weight: " << quadrature_weight << " p-weight: " << pressure_weight;
    return ss.str();
  }
};

struct ParticleParticlePair
{
  int pidx[2];
  scalar d;
  scalar r;
  
  std::string toString() const
  {
    std::stringstream ss;
    ss << "pidx: (" << pidx[0] << ", " << pidx[1] << ") d: " << d;
    return ss.str();
  }
};

std::ostream& operator<<(std::ostream &os, const EdgeEdgePair &info);
std::ostream& operator<<(std::ostream &os, const ParticleEdgePair &info);
std::ostream& operator<<(std::ostream &os, const PointEdgePair &info);
std::ostream& operator<<(std::ostream &os, const ParticleParticlePair &info);

template<int DIM>
class PolygonalCohesion: public Force
{
  const int m_meddling_stencil = 2;
  const int m_num_quadrature = 1;
  
  TwoDScene<DIM>* m_parent;
  std::vector< std::unordered_map<int, ParticleEdgePair*> > m_adjacency_categorized; // pidx -> (hair_idx, pair_idx)
  std::vector< std::unordered_map<int, EdgeEdgePairEEC*> > m_edge_connections; // eidx -> (neighbor_eidx_idx -> EE_pair)
  std::vector< int > m_num_valid_edge_connections; // eidx -> count of valid EEP
  std::vector< int > m_num_edge_connections; // pidx -> count
  std::vector< std::vector<int> > m_particle_to_pppairs; // pidx -> other particles connected
  // std::vector< std::vector< EdgeEdgePair > > m_edge_edge_contacts;
  std::vector< int > m_num_adjacency_categorized;
  std::vector< int > m_counting_valid_adjacency;
  std::vector< int > m_counting_pp_pair_location;
  std::vector< uint64_t > m_counting_pp_pairs;
  std::vector< ParticleParticlePair > m_particle_particle_pairs;
  std::vector< PointEdgePair > m_point_edge_pairs;
  
  std::vector< std::vector<PointEdgePair> > m_point_edge_pairs_cache;
  std::vector< int > m_counting_poe_pair_location;
  
  VectorXs m_particle_adjacency_hair_size;
  
  std::unordered_map<int, std::unordered_set<int> > m_adjacency_hair_edges_buffer;
  
  std::vector< std::vector<int> > m_particle_to_point_edge_pairs;
  
  std::vector< std::unordered_set<int> > m_pp_pair_hash;
  
  VectorXs m_particle_length;
  
  nanoflann::KDTreeEigenMatrixAdaptor< MatrixXs, DIM>* m_tree;
  Sorter* m_sorter;
  
  MatrixXs m_edge_buffer;
  
  // buffers
  VectorXs m_gradE;
  MatrixXs m_hessE;
  MatrixXs m_hessV;
  VectorXs m_pair_counts;
  
  TripletXs m_hess_buffer;
  bool m_use_decoupled_force;
  bool m_compute_particle_poe_mapping;
  
  // variables for inter-hair shallow water dynamics

  SparseXs m_W_fv_interhair_T;
  SparseXs m_gradF_global;
  SparseXs m_gradF_global_T;
  SparseXs m_gradF_interhair;
  SparseXs m_gradF_interhair_T;
  SparseXs m_dir_f_interhair;
  SparseXs m_dir_f_interhair_T;
  
  VectorXs m_G_f_global;
  VectorXs m_iG_v_global;
  VectorXs m_G_f_interhair;

  VectorXs m_u_interhair;
  VectorXs m_divu_interhair;
  
 
  VectorXs m_area_v_global;
  VectorXs m_area_e_global;
  VectorXs m_cur_eta_v_global;
  VectorXs m_pressure_v_global;
  
  MatrixXs m_rhs_offset_v_global;
  MatrixXs m_cur_hhrr_v_global;
  
  VectorXs m_ce_inter_buffer;
  VectorXs m_ce_inter_short_buffer;
  VectorXs m_ce_global_buffer;
  VectorXs m_cv_buffer;
  
  std::vector<int> m_pplink_count;
  
  std::unordered_map<uint64, scalar> m_poep_lambda;
  
  std::mutex m_pep_mutex;
  std::mutex m_eep_mutex;
  std::mutex m_poep_mutex;
  
  CohesionTable* m_min_cohesion_table;
  CohesionTable* m_max_cohesion_table;

  std::unordered_map<uint64, scalar> m_viscous_start_phi;

public:
  PolygonalCohesion(TwoDScene<DIM>* scene);
  
  virtual ~PolygonalCohesion();
  
  virtual void addEnergyToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, scalar& E );
  
  virtual void addGradEToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, VectorXs& gradE );
  
  virtual void addHessXToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, TripletXs& hessE );
  
  virtual void addHessVToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, TripletXs& hessE );
  
  virtual void addGradEToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, VectorXs& gradE, int pidx );
  
  virtual void addHessXToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, VectorXs& hessE, int pidx );
  
  virtual void addHessVToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, VectorXs& hessE, int pidx );
  
  virtual void computeIntegrationVars( const VectorXs& x, const VectorXs& v, const VectorXs& m,
                                      VectorXs& lambda, VectorXs& lambda_v,
                                      TripletXs& J, TripletXs& Jv, TripletXs& Jxv, TripletXs& tildeK,
                                      TripletXs& stiffness, TripletXs& damping, VectorXs& Phi, VectorXs& Phiv, const scalar& dt);
  
  virtual int numConstraintPos();
  
  virtual int numConstraintVel();
  
  virtual int numJ();
  
  virtual int numJv();
  
  virtual int numJxv();
  
  virtual int numTildeK();
  
  virtual bool isParallelized();
  
  virtual bool isPrecomputationParallelized();
  
  virtual void storeLambda(const VectorXs& lambda, const VectorXs& lambda_v);
  
  virtual Force* createNewCopy();
  
  virtual const char* name();
  
  virtual void getAffectedVars( int pidx, std::unordered_set<int>& vars );
  
  virtual bool isContained( int pidx );
  
  virtual void preCompute( const VectorXs& x, const VectorXs& v, const VectorXs& m, const scalar& dt );

  virtual void setUseDecoupledForce(bool decoupled);
  
  virtual void setUseParticlePOEMap(bool ppoemap);
  
  virtual void computeInterHairVariables();
  
  virtual void computeGlobalHairVariables();
  
  virtual void computeGlobalHairPressure();
  
  virtual void computeInterHairVelocity(const scalar& dt);
  
  virtual void computeInterHairRHS(const scalar& dt);
  
  // get-set for film flow
  const VectorXs& getAreaVGlobal() const;
  
  VectorXs& getAreaVGlobal();
  
  const VectorXs& getAreaEGlobal() const;
  
  VectorXs& getAreaEGlobal();
  
  const VectorXs& getPressureVGlobal() const;
  
  VectorXs& getPressureVGlobal();
  
  const MatrixXs& getRhsOffsetVGlobal() const;
  
  MatrixXs& getRhsOffsetVGlobal();
  
  const std::vector<int>& getPPPCountV() const;
  
  virtual void updateStructure(const VectorXs& x);
  
  virtual void updatePointEdgePairs(const VectorXs& x, const VectorXs& v, const scalar& dt );
  
  virtual bool isInterHair() const;
  
  virtual void postStepScene(const scalar& dt );
  
  virtual void write(std::vector<scalar>& data) const;
  
  virtual void read(const scalar* data);
  
  virtual void writeReadable(std::ostream& oss_pepairs, std::ostream& oss_poepairs, std::ostream& oss_pppairs) const;
  
  virtual size_t adjacency_size() const;
  
  virtual scalar getDStarPlanar(const scalar& radius) const;
  
  virtual scalar getStiffnessPlanar(const scalar& radius, const scalar& d0, const scalar& A_target, const scalar& pressure_weight) const;
  
  virtual scalar getDStar(const scalar& radius) const;
  
  virtual scalar getStiffness(const scalar& radius, const scalar& d0, const scalar& A_target, const scalar& pressure_weight) const;
  
  void updateViscousStartPhi( const VectorXs& x );

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
protected:
  scalar computeLocalCoord(int eidx, const scalar& alpha) const;
  void extractLocalCoord(const scalar& local_coord, int hidx, int& eidx, scalar& alpha) const;
  
  void computeConnections(const VectorXs& hair_vertices);
  void computeParticleEdgeCounts(int pe_idx);

  void buildSearchTree(const VectorXs& );
 
  void updatePorosityEEC();
  
  void updatePorosity();
  
  scalar computeEdgeThickness(const scalar& vol_p, const int_Vectors_scalar<DIM>& ep);
  
  bool checkPairwisePPHash(int i, int j);
  void markPairwisePPHash(int i, int j);
  inline uint64_t makePairwisePPHash(int i, int j);

  void findParticleEdgePairs(const VectorXs& x, int pidx);
  
  void findParticleParticlePairs(const VectorXs& x);

  void findParticleParticlePairsEEC(const VectorXs& x);
  void findEdgeEdgeContact( const VectorXs x, const VectorXs v, const scalar& dt, const int& base_eidx );
  void findPointEdgePairsEEC(const VectorXs& x, int base_eidx, std::vector<PointEdgePair>& poepairs);

  void findPointEdgePairs(const VectorXs& x, int base_eidx, std::vector<PointEdgePair>& poepairs);
  
  void accumulateEFJPairwise(const VectorXs& x, const VectorXs& v, VectorXs& gradE, MatrixXs& hessE, MatrixXs& hessV);
};


#endif
