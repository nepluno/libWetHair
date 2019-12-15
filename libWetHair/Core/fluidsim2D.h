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


#ifndef FLUIDSIM2D_H
#define FLUIDSIM2D_H

#include "fluidsim.h"
#include "MathUtilities.h"
#include "array2.h"
#include "pcgsolver/pcg_solver.h"

//#define USE_SURFACE_TENSION

#include <vector>

template<int DIM>
class TwoDScene;

class Sorter;

template<int DIM>
class FluidDragForce;

class FluidSim2D : public FluidSim
{
  
public:
  
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  
  FluidSim2D(const Vector2s& origin_, scalar width, int ni_, int nj_, const scalar& dz_,
             const std::vector< Boundary<2>* >& boundaries_, const std::vector< SourceBoundary<2>* >& sources_, TwoDScene<2>* scene);

  virtual ~FluidSim2D();
  
  virtual void advect_boundary(const scalar& dt);
  virtual void update_boundary();
  virtual void init_random_particles(const scalar& rl, const scalar& rr, const scalar& rb, const scalar& rt);
  virtual void init_hair_particles();
  virtual void controlSources(const scalar& current_time, const scalar& dt);

  const std::vector<Particle<2> >& get_particles() const;
  const std::vector< Boundary<2>* >& get_boundaries() const;
  const std::vector< SourceBoundary<2>* >& get_sources() const;
  const Vector2s& get_origin() const;
  int get_ni() const;
  int get_nj() const;
  
  int get_u_ni() const;
  int get_v_ni() const;
  
  int get_u_nj() const;
  int get_v_nj() const;
  
  scalar compute_phi_vel(const Vector2s& pos, Vector2s& vel) const;
  
  Vector2s get_velocity(const Vector2s& position, const Array2s& u, const Array2s& v) const;
  Vector2s get_velocity(const Vector2s& position) const;
  scalar get_nodal_solid_phi(const Vector2s& position) const;
  Vector2s get_nodal_solid_phi_gradient(const Vector2s& position) const;
  Vector2s get_solid_velocity(const Vector2s& position) const;
  Vector2s get_pressure_gradient(const Vector2s& position) const;
  scalar get_pressure(const Vector2s& position) const;
  Vector2s get_hair_velocity(const Vector2s& position) const;
  Vector2s get_particle_velocity(const Vector2s& position) const;
  Vector2s get_particle_drag(const Vector2s& position) const;
  Matrix2s get_affine_matrix(const Vector2s& position) const;
  
  std::vector< std::vector<EdgeVelDragIntersection<2> > >& get_u_edge_vel_drag();
  std::vector< std::vector<EdgeVelDragIntersection<2> > >& get_v_edge_vel_drag();
  
  scalar getLiquidPhiValue(const Vector2s& position) const;
  scalar getClampedLiquidPhiValue(const Vector2s& position) const;
  
  virtual int num_particles() const;
  
  virtual void advect_particles(scalar dt);
  virtual void constrain_hair_particles();
  
  //fluid velocity operations
  virtual void add_drag(scalar dt);
  virtual void add_gravity(scalar dt);
  
  virtual void project(scalar dt);
  virtual void compute_weights();
  virtual void solve_pressure(scalar dt);
  virtual scalar get_particle_weight(const Vector2s& position) const;
  virtual scalar get_clamped_particle_weight(const Vector2s& position) const;
  
  virtual void constrain_velocity();
  virtual void prepare_update_from_hair();
  virtual void done_update_from_hair();
  
  virtual scalar cfl();
  virtual scalar cellsize() const;
  
  virtual scalar default_radius_multiplier() const;
  
  virtual int default_particle_in_cell() const;
  
  virtual void add_particle(const Particle<2>& p);
  
  virtual void map_p2g(bool with_hair_particles);
  virtual void map_g2p_apic();
  
  virtual void combine_velocity_field();

  virtual void correct(scalar dt);
  virtual void resample(Vector2s& p, Vector2s& u, Matrix2s& c);
  
  virtual void shareParticleWithHairs( VectorXs& x, scalar dt );
  
  virtual void transferLiquidToGridParticle(const scalar& dt);
  
  virtual void compute_liquid_phi();
  
  virtual void save_pressure(const std::string szfn);

  virtual void sort_particles();
  
  virtual Vector2s computeParticleMomentum();
  
  virtual scalar computeParticleAngularMomentum();
  
  virtual Vector2s computeHairGridMomentum();
  
  virtual scalar computeHairGridAngularMomentum();
  
  virtual Vector2s computeParticleGridMomentum();
  
  virtual scalar computeParticleGridAngularMomentum();
  
  virtual Vector2s computeReweightedHairGridMomentum();
  
  virtual Vector2s computeReweightedParticleGridMomentum();
  
  virtual scalar computeReweightedHairGridAngularMomentum();
  
  virtual scalar computeReweightedParticleGridAngularMomentum();
  
  virtual Vector2s computeParticleWeightedCombinedGridMomentum();
  
  virtual Vector2s computeHairWeightedCombinedGridMomentum();
  
  virtual scalar computeParticleWeightedCombinedGridAngularMomentum();
  
  virtual scalar computeHairWeightedCombinedGridAngularMomentum();
  
  virtual Vector2s computeCombinedGridMomentum();
  
  virtual scalar computeCombinedGridAngularMomentum();
  
  virtual scalar computeTotalLiquidVol() const;
  
  virtual scalar computeParticleKineticEnergy();
  
  virtual scalar computeHairGridKineticEnergy();
  
  virtual scalar computeParticleGridKineticEnergy();
  
  virtual scalar computeCombinedGridKineticEnergy();
  
  virtual scalar computeOverallDivergence();

  virtual Vector2s computeHairGridDrag();
  
  virtual void preCompute( const VectorXs& x, const VectorXs& v, const VectorXs& m, const scalar& dt );
  
  virtual void addGradEToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, VectorXs& gradE );

  virtual Vector3s getMinBBX() const;
  
  virtual Vector3s getMaxBBX() const;
  
  virtual scalar dropvol(const scalar& radii) const;
  
  virtual scalar dropradius(const scalar& vol) const;
  
  virtual void add_particle(const VectorXs& pos, const VectorXs& vel, const scalar& radii, ParticleType type);
  
  virtual void write(std::vector<scalar>& data) const;
  
  virtual void writeReadable(std::vector<std::ostringstream>&) const;

  virtual void readReadable( std::ifstream& file );
  
  virtual void read(const scalar* data, size_t size_particles, size_t size_boundaries);
  
  virtual size_t particle_size() const;
  
  virtual size_t boundary_size() const;
  
  virtual size_t crucial_grid_size() const;
  
  virtual void apply_viscosity( scalar dt );
  
protected:
  
  int ryoichi_correction_counter;
  
  // Boundaries
  std::vector< Boundary<2>* > boundaries;
  std::vector< SourceBoundary<2>* > sources;
  
  // Grid Origin
  Vector2s origin;
  
  // Grid dimensions
  int ni,nj;
  scalar dx;
  scalar dz;
  
  // Fluid velocity
  Array2s u, v;
  Array2s temp_u, temp_v;
  
  Array2s u_weight_hair;
  Array2s v_weight_hair;
  
  Array2s u_weight_particle;
  Array2s v_weight_particle;
  
  Array2s u_weight_total;
  Array2s v_weight_total;
  
  Array2s u_hair, v_hair;
  Array2s u_pressure_grad, v_pressure_grad;
  Array2s u_particle, v_particle;
  Array2s u_solid, v_solid;
  
  // Hair -> Voxel Intersections
  std::vector< std::vector<EdgeVelDragIntersection<2> > > u_edge_vel_drag;
  std::vector< std::vector<EdgeVelDragIntersection<2> > > v_edge_vel_drag;
  
  std::vector< int > u_num_edge_voxel_intersections;
  std::vector< int > v_num_edge_voxel_intersections;
  
  std::vector< EdgeVelDragIntersection<2> > u_vel_drag;
  std::vector< EdgeVelDragIntersection<2> > v_vel_drag;
  
  // Tracer particles
  std::vector<Particle<2> > particles;
  
  // Static geometry representation
  Array2s nodal_solid_phi;
  Array2s u_weights, v_weights;
  Array2c u_valid, v_valid;
  
  Array2s liquid_phi;
  
  Array2s u_drag;
  Array2s v_drag;
  
  // Data arrays for extrapolation
  Array2c valid, old_valid;
  
  std::vector<FluidDragForce<2>*> drag_forces;
  
  // Solver data
  robertbridson::PCGSolver<scalar> solver;
  robertbridson::SparseMatrix<scalar> matrix;
  std::vector<double> rhs;
  std::vector<double> pressure;
  
  TwoDScene<2>* m_parent;
  Sorter* m_sorter;
  
  std::vector< std::vector<int> > m_pool_liquid_index_cache;
  std::vector< std::vector<Particle<2> > > m_pool_liquid_particle_cache;
  std::vector< std::vector<Particle<2> > > m_regular_liquid_particle_cache;
  VectorXs m_pool_liquid_vol_cache;
  
  std::vector<HairParticleBridge<2> > m_bridges;
  std::vector< std::vector<HairParticleBridge<2> > > m_hair_bridge_buffer;
  
  std::vector<int> m_hair_particle_affected;
};

#endif
