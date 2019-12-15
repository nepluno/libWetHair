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


#ifndef FLUIDSIM3D_H
#define FLUIDSIM3D_H

#include "fluidsim.h"
#include "MathUtilities.h"
#include "array3.h"
#include "pcgsolver/pcg_solver.h"

//#define USE_SURFACE_TENSION

#include <vector>

template<int DIM>
class TwoDScene;

class Sorter;

template<int DIM>
class FluidDragForce;

class FluidSim3D : public FluidSim
{
  
public:
  
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  
  FluidSim3D(const Vector3s& origin_, scalar width, int ni_, int nj_, int nk_,
             const std::vector< Boundary<3>* >& boundaries_, const std::vector< SourceBoundary<3>* >& sources_, TwoDScene<3>* scene);

  virtual ~FluidSim3D();
  
  virtual void advect_boundary(const scalar& dt);
  virtual void update_boundary();
  virtual void init_random_particles(const scalar& rl, const scalar& rr, const scalar& rb, const scalar& rt, const scalar& rf, const scalar& rk);
  virtual void init_hair_particles();
  virtual void controlSources(const scalar& current_time, const scalar& dt);
  
  const std::vector<Particle<3> >& get_particles() const;
  const std::vector< Boundary<3>* >& get_boundaries() const;
  const std::vector< SourceBoundary<3>* >& get_sources() const;
  const Vector3s& get_origin() const;
  int get_ni() const;
  int get_nj() const;
  int get_nk() const;
  
  int get_u_ni() const;
  int get_v_ni() const;
  int get_w_ni() const;
  
  int get_u_nj() const;
  int get_v_nj() const;
  int get_w_nj() const;
  
  int get_u_nk() const;
  int get_v_nk() const;
  int get_w_nk() const;
  
  scalar compute_phi_vel(const Vector3s& pos, Vector3s& vel) const;
  
  Vector3s get_velocity(const Vector3s& position, const Array3s& u, const Array3s& v, const Array3s& w) const;
  Vector3s get_velocity(const Vector3s& position) const;
  scalar get_nodal_solid_phi(const Vector3s& position) const;
  Vector3s get_nodal_solid_phi_gradient(const Vector3s& position) const;
  Vector3s get_solid_velocity(const Vector3s& position) const;
  Vector3s get_pressure_gradient(const Vector3s& position) const;
  Vector3s get_visc_impulse(const Vector3s& position) const;
  scalar get_pressure(const Vector3s& position) const;
  Vector3s get_hair_velocity(const Vector3s& position) const;
  Vector3s get_particle_velocity(const Vector3s& position) const;
  Vector3s get_temp_velocity(const Vector3s& position) const;
  Vector3s get_particle_drag(const Vector3s& position) const;
  Matrix3s get_affine_matrix(const Vector3s& position) const;
  
  std::vector< std::vector<EdgeVelDragIntersection<3> > >& get_u_edge_vel_drag();
  std::vector< std::vector<EdgeVelDragIntersection<3> > >& get_v_edge_vel_drag();
  std::vector< std::vector<EdgeVelDragIntersection<3> > >& get_w_edge_vel_drag();
  
  scalar getLiquidPhiValue(const Vector3s& position) const;
  scalar getClampedLiquidPhiValue(const Vector3s& position) const;
  
  //tracer particle operations
  virtual int num_particles() const;
  
  virtual void advect_particles(scalar dt);
  virtual void constrain_hair_particles();
  
  //fluid velocity operations
  virtual void add_drag(scalar dt);
  virtual void add_gravity(scalar dt);
  
  virtual void project(scalar dt);
  virtual void compute_weights();
  virtual void solve_pressure(scalar dt);
  virtual scalar get_particle_weight(const Vector3s& position) const;
  virtual scalar get_clamped_particle_weight(const Vector3s& position) const;
  
  virtual void constrain_velocity();
  virtual void prepare_update_from_hair();
  virtual void done_update_from_hair();
  
  virtual scalar cfl();
  virtual scalar cellsize() const;
  
  virtual scalar default_radius_multiplier() const;
  
  virtual int default_particle_in_cell() const;
  
  virtual void add_particle(const Particle<3>& p);
  
  virtual void map_p2g(bool with_hair_particles);
  virtual void map_g2p_apic();
  
  virtual void combine_velocity_field();

  virtual void correct(scalar dt);
  virtual void resample(Vector3s& p, Vector3s& u, Matrix3s& c);
  
  virtual void apply_viscosity( scalar dt );
  
  virtual void compute_viscosity_weights();
  
  virtual void shareParticleWithHairs( VectorXs& x, scalar dt );
  
  virtual void transferLiquidToGridParticle(const scalar& dt);
  
  virtual void compute_liquid_phi();
  
  virtual void save_pressure(const std::string szfn);

  virtual void sort_particles();
  
  virtual Vector3s computeParticleMomentum();
  
  virtual Vector3s computeParticleAngularMomentum();
  
  virtual Vector3s computeParticleGridAngularMomentum();
  
  virtual Vector3s computeReweightedParticleGridAngularMomentum();
  
  virtual Vector3s computeCombinedGridMomentum();
  
  virtual Vector3s computeCombinedGridAngularMomentum();
  
  virtual scalar computeTotalLiquidVol() const;
  
  virtual scalar computeParticleKineticEnergy();
  
  virtual scalar computeParticleGridKineticEnergy();
  
  virtual scalar computeCombinedGridKineticEnergy();
  
  virtual scalar computeOverallDivergence();

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
  
protected:
  // Boundaries
  std::vector< Boundary<3>* > boundaries;
  std::vector< SourceBoundary<3>* > sources;
  
  // Grid Origin
  Vector3s origin;
  
  // Grid dimensions
  int ni,nj,nk;
  scalar dx;
  
  // Fluid velocity
  Array3s u, v, w;
  Array3s temp_u, temp_v, temp_w;
  
  Array3s u_weight_particle;
  Array3s v_weight_particle;
  Array3s w_weight_particle;
  
  Array3s u_weight_total;
  Array3s v_weight_total;
  Array3s w_weight_total;
  
  Array3s u_pressure_grad, v_pressure_grad, w_pressure_grad;
  Array3s u_particle, v_particle, w_particle;
  Array3s u_solid, v_solid, w_solid;
  
  // Hair -> Voxel Intersections
  std::vector< std::vector<EdgeVelDragIntersection<3> > > u_edge_vel_drag;
  std::vector< std::vector<EdgeVelDragIntersection<3> > > v_edge_vel_drag;
  std::vector< std::vector<EdgeVelDragIntersection<3> > > w_edge_vel_drag;
  
  std::vector< int > u_num_edge_voxel_intersections;
  std::vector< int > v_num_edge_voxel_intersections;
  std::vector< int > w_num_edge_voxel_intersections;
  
  std::vector< EdgeVelDragIntersection<3> > u_vel_drag;
  std::vector< EdgeVelDragIntersection<3> > v_vel_drag;
  std::vector< EdgeVelDragIntersection<3> > w_vel_drag;
  
  // Tracer particles
  std::vector<Particle<3> > particles;
  
  // Static geometry representation
  Array3s nodal_solid_phi;
  Array3s cell_solid_phi;
  Array3s u_weights, v_weights, w_weights;
  Array3c u_valid, v_valid, w_valid;
  
  Array3s liquid_phi;
  
  Array3s u_drag;
  Array3s v_drag;
  Array3s w_drag;
  
  // Data arrays for extrapolation
  Array3c valid, old_valid;
  
  std::vector<FluidDragForce<3>*> drag_forces;
  
  //Data for viscosity solve
  Array3s u_vol_liquid, v_vol_liquid, w_vol_liquid,
  ex_vol_liquid, ey_vol_liquid, ez_vol_liquid, c_vol_liquid;
  
  Array3s u_visc_impulse, v_visc_impulse, w_visc_impulse;
  
  // Solver data
  robertbridson::PCGSolver<scalar> solver;
  robertbridson::SparseMatrix<scalar> matrix;
  std::vector<double> rhs;
  std::vector<double> pressure;
  
  TwoDScene<3>* m_parent;
  Sorter* m_sorter;
  
  std::vector< std::vector<int> > m_pool_liquid_index_cache;
  std::vector< std::vector<Particle<3> > > m_pool_liquid_particle_cache;
  std::vector< std::vector<Particle<3> > > m_regular_liquid_particle_cache;
  VectorXs m_pool_liquid_vol_cache;
  
  std::vector<HairParticleBridge<3> > m_bridges;
  std::vector< std::vector<HairParticleBridge<3> > > m_hair_bridge_buffer;
  
  std::vector<int> m_hair_particle_affected;
  
  int ryoichi_correction_counter;
  
};

#endif
