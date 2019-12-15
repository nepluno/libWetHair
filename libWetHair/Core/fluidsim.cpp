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


#include "fluidsim.h"
#include "TwoDScene.h"
#include "fluidsim2D.h"
#include "fluidsim3D.h"

inline scalar circle_phi(const Vector2s& position, const Vector2s& centre, scalar radius) {
  return ((position-centre).norm() - radius);
}

inline scalar box_phi(const Vector2s& position, const Vector2s& centre, const Vector2s& expand, const scalar& radius) {
  scalar dx = fabs(position[0] - centre[0]) - expand[0];
  scalar dy = fabs(position[1] - centre[1]) - expand[1];
  scalar dax = max(dx, 0.0);
  scalar day = max(dy, 0.0);
  return min(max(dx, dy), 0.0) + sqrt(dax * dax + day * day) - radius;
}

inline scalar circle_phi(const Vector3s& position, const Vector3s& centre, scalar radius) {
  return ((position-centre).norm() - radius);
}

inline scalar capsule_phi(const Vector2s& position, const Vector2s& centre, const scalar& radius, const scalar& halflength)
{
  Vector2s a = centre - Vector2s(halflength, 0);
  Vector2s pa = position - a;
  Vector2s ba = Vector2s(2.0 * halflength, 0);
  scalar h = mathutils::clamp(pa.dot(ba) / (4.0 * halflength * halflength), 0.0, 1.0);
  return (pa - ba * h).norm() - radius;
}


inline scalar box_phi(const Vector3s& position, const Vector3s& centre, const Vector3s& expand, const scalar& radius) {
  scalar dx = fabs(position[0] - centre[0]) - expand[0];
  scalar dy = fabs(position[1] - centre[1]) - expand[1];
  scalar dz = fabs(position[2] - centre[2]) - expand[2];
  scalar dax = max(dx, 0.0);
  scalar day = max(dy, 0.0);
  scalar daz = max(dz, 0.0);
  return min(max(max(dx, dy), dz), 0.0) + sqrt(dax * dax + day * day + daz * daz) - radius;
}

inline scalar capsule_phi(const Vector3s& position, const Vector3s& centre, const scalar& radius, const scalar& halflength)
{
  Vector3s a = centre - Vector3s(halflength, 0, 0);
  Vector3s pa = position - a;
  Vector3s ba = Vector3s(2.0 * halflength, 0, 0);
  scalar h = mathutils::clamp(pa.dot(ba) / (4.0 * halflength * halflength), 0.0, 1.0);
  return (pa - ba * h).norm() - radius;
}

template<int DIM>
FluidSim::Boundary<DIM>::Boundary(BOUNDARY_TYPE type_, int parent_)
: type(type_), parent(parent_)
{
}

template<int DIM>
FluidSim::SolidBoundary<DIM>::SolidBoundary(const Vectors<DIM>& center_, const VectorXs& parameter_, BOUNDARY_TYPE type_, bool inside, const Vector3s& raxis, const scalar& rangle)
: Boundary<DIM>(type_), center(center_), parameter(parameter_), sign(inside ? -1.0 : 1.0), rot(Eigen::AngleAxis<scalar>(rangle, raxis))
{
  V.setZero();
  omega.setZero();
  
  future_center = center;
  future_rot = rot;
}

template<int DIM>
FluidSim::SourceBoundary<DIM>::SourceBoundary(const Vectors<DIM>& center_, const VectorXs& parameter_, BOUNDARY_TYPE type_, bool inside, const Vector3s& raxis, const scalar& rangle, const Vectors<DIM>& eject_vel_, const scalar& start_, const scalar& end_, const scalar& spray_angle_, const scalar& drop_radius_prop_, const scalar& sub_activate_length_, const scalar& sub_inactivate_length_)
: SolidBoundary<DIM>(center_, parameter_, type_, inside, raxis, rangle), eject_vel(eject_vel_), last_parent(NULL), start(start_), end(end_), activated(false), spray_angle(spray_angle_), drop_radius_prop(drop_radius_prop_), sub_activate_length(sub_activate_length_), sub_inactivate_length(sub_inactivate_length_)
{}


template<>
scalar FluidSim::SolidBoundary<3>::compute_phi_vel(const Vectors<3>& pos, Vectors<3>& vel) const
{
  scalar phi = 0.0;
  
  Vector3s dx = pos - center;
  
  switch (type) {
    case BT_BOX:
    {
      Eigen::Quaternion<scalar> p0(0.0, dx(0), dx(1), dx(2));
      Eigen::Quaternion<scalar> irot = rot.conjugate();
      Vector3s rotp = (irot * p0 * irot.inverse()).vec() + center;
      phi = sign * box_phi(rotp, center, Vector3s(parameter(0), parameter(1), parameter(2)), parameter(3));
      break;
    }
      
    case BT_CIRCLE:
    {
      phi = sign * circle_phi(pos, center, parameter(0));
      break;
    }
      
    case BT_CAPSULE:
    {
      Eigen::Quaternion<scalar> p0(0.0, dx(0), dx(1), dx(2));
      Eigen::Quaternion<scalar> irot = rot.conjugate();
      Vector3s rotp = (irot * p0 * irot.inverse()).vec() + center;
      phi = sign * capsule_phi(rotp, center, parameter(0), parameter(1));
      break;
    }
    default:
      break;
  }
  
  vel = V + omega.cross(dx);
  
  return phi;
}

template<>
scalar FluidSim::SolidBoundary<2>::compute_phi_vel(const Vectors<2>& pos, Vectors<2>& vel) const
{
  scalar phi = 0.0;
  
  Vector2s dx = pos - center;
  
  switch (type) {
    case BT_BOX:
    {
      Eigen::Quaternion<scalar> p0(0.0, dx(0), dx(1), 0.0);
      Eigen::Quaternion<scalar> irot = rot.conjugate();
      Vector3s rotp = (irot * p0 * irot.inverse()).vec();
      Vector2s rotp2 = Vector2s(rotp(0), rotp(1)) + center;
      phi = sign * box_phi(rotp2, center, Vector2s(parameter(0), parameter(1)), parameter(3));
      break;
    }
    case BT_CIRCLE:
    {
      phi = sign * circle_phi(pos, center, parameter(0));
      break;
    }
    case BT_CAPSULE:
    {
      Eigen::Quaternion<scalar> p0(0.0, dx(0), dx(1), 0.0);
      Eigen::Quaternion<scalar> irot = rot.conjugate();
      Vector3s rotp = (irot * p0 * irot.inverse()).vec();
      Vector2s rotp2 = Vector2s(rotp(0), rotp(1)) + center;
      phi = sign * capsule_phi(rotp2, center, parameter(0), parameter(1));
      break;
    }
    default:
      break;
  }
  
  
  vel = V + omega(2) * Vector2s(-dx[1], dx[0]);
  
  return phi;
}

template<>
void FluidSim::SourceBoundary<2>::sample(FluidSim* parent)
{
  if(!parent) return;
  
  FluidSim2D* fluid2d = (FluidSim2D*) parent;
  
  detectors.clear();
  
  int ni = fluid2d->get_ni();
  int nj = fluid2d->get_nj();
  scalar cellsize = fluid2d->cellsize();
  
  for(int j = 0; j < nj; ++j) for(int i = 0; i < ni; ++i)
  {
    Vector2s pos = Vector2s((i + 0.5) * cellsize, (j + 0.5) * cellsize) + fluid2d->get_origin();
    Vector2s v;
    
    if(compute_phi_vel(pos, v) < 0.0) {
      detectors.push_back(pos);
    }
  }
  
  last_parent = parent;
}

template<>
void FluidSim::SourceBoundary<3>::sample(FluidSim* parent)
{
  if(!parent) return;
  
  FluidSim3D* fluid3d = (FluidSim3D*) parent;
  
  int ni = fluid3d->get_ni();
  int nj = fluid3d->get_nj();
  int nk = fluid3d->get_nk();
  scalar cellsize = fluid3d->cellsize();
  
  detectors.clear();
  for(int k = 0; k < nk; ++k) for(int j = 0; j < nj; ++j) for(int i = 0; i < ni; ++i)
  {
    Vector3s pos = Vector3s((i + 0.5) * cellsize, (j + 0.5) * cellsize, (k + 0.5) * cellsize) + fluid3d->get_origin();
    Vector3s v;
    
    if(compute_phi_vel(pos, v) < 0.0) {
      detectors.push_back(pos);
    }
  }
  
  last_parent = parent;
}




template<>
void FluidSim::SourceBoundary<2>::advance(const scalar& dt)
{
  Eigen::Quaternion<scalar> p = Eigen::Quaternion<scalar>(0.0, eject_vel(0), eject_vel(1), 0);
  Vector3s ev = (FluidSim::SolidBoundary<2>::future_rot * FluidSim::SolidBoundary<2>::rot.inverse() * p * FluidSim::SolidBoundary<2>::rot * FluidSim::SolidBoundary<2>::future_rot.inverse()).vec();
  eject_vel(0) = ev(0); eject_vel(1) = ev(1);
  
  FluidSim::SolidBoundary<2>::advance(dt);

  sample(last_parent);
}

template<>
void FluidSim::SourceBoundary<3>::advance(const scalar& dt)
{
  Eigen::Quaternion<scalar> p = Eigen::Quaternion<scalar>(0.0, eject_vel(0), eject_vel(1), eject_vel(2));
  eject_vel = (FluidSim::SolidBoundary<3>::future_rot * FluidSim::SolidBoundary<3>::rot.inverse() * p * FluidSim::SolidBoundary<3>::rot * FluidSim::SolidBoundary<3>::future_rot.inverse()).vec();
  
  FluidSim::SolidBoundary<3>::advance(dt);

  sample(last_parent);
}

template<int DIM>
void FluidSim::SolidBoundary<DIM>::advance(const scalar& dt)
{
  V = (future_center - center) / dt;
  Eigen::Quaternion<scalar> q = future_rot * rot.conjugate();
  scalar len = q.vec().norm();
  if(len > 0.0) {
    scalar angle = 2.0 * atan2(len, q.w());
    omega = q.vec() / len * angle / dt;
  } else {
    omega.setZero();
  }
  
  center = future_center;
  rot = future_rot;
}

template<int DIM>
void FluidSim::SolidBoundary<DIM>::write(std::vector<scalar>& buf) const
{
  for(int r = 0; r < DIM; ++r) {
    buf.push_back(center(r));
  }
  for(int r = 0; r < DIM; ++r) {
    buf.push_back(parameter(r));
  }
  
  buf.push_back(rot.x());
  buf.push_back(rot.y());
  buf.push_back(rot.z());
  buf.push_back(rot.w());
  
  for(int r = 0; r < 3; ++r) {
    buf.push_back(omega(r));
  }
  
  for(int r = 0; r < DIM; ++r) {
    buf.push_back(V(r));
  }
}

template<int DIM>
void FluidSim::SolidBoundary<DIM>::read(const scalar* data)
{
  int k = 0;
  for(int r = 0; r < DIM; ++r)
  {
    center(r) = data[k++];
  }
  for(int r = 0; r < DIM; ++r) {
    parameter(r) = data[k++];
  }
  rot.x() = data[k++];
  rot.y() = data[k++];
  rot.z() = data[k++];
  rot.w() = data[k++];
  
  for(int r = 0; r < 3; ++r) {
    omega(r) = data[k++];
  }
  
  for(int r = 0; r < DIM; ++r) {
    V(r) = data[k++];
  }
  
  future_rot = rot;
  future_center = center;
}

template<int DIM>
size_t FluidSim::SolidBoundary<DIM>::size() const
{
  return (DIM + DIM + 4 + 3 + DIM) * sizeof(scalar);
}
template<int DIM>
FluidSim::OperatorBoundary<DIM>::OperatorBoundary(BOUNDARY_TYPE type_)
: Boundary<DIM>(type_)
{
}

template<int DIM>
void FluidSim::OperatorBoundary<DIM>::advance(const scalar& dt)
{
  int nb = children.size();
  threadutils::thread_pool::ParallelFor(0, nb, [&] (int i) {
    children[i]->advance(dt);
  });
}

template<int DIM>
scalar FluidSim::OperatorBoundary<DIM>::compute_phi_vel(const Vectors<DIM>& pos, Vectors<DIM>& vel) const
{
  switch (Boundary<DIM>::type) {
    case BT_UNION:
    {
      scalar min_phi = 1e+20;
      
      for(const Boundary<DIM>* child : children)
      {
        Vectors<DIM> sub_vel;
        scalar phi = child->compute_phi_vel(pos, sub_vel);
        if(phi < min_phi)
        {
          min_phi = phi;
          vel = sub_vel;
        }
      }
      
      return min_phi;
    }
    case BT_INTERSECT:
    {
      scalar max_phi = -1e+20;
      
      for(const Boundary<DIM>* child : children)
      {
        Vectors<DIM> sub_vel;
        scalar phi = child->compute_phi_vel(pos, sub_vel);
        if(phi > max_phi)
        {
          max_phi = phi;
          vel = sub_vel;
        }
      }
      
      return max_phi;
    }
    default:
      vel = Vectors<DIM>::Zero();
      return 1e+20;
  }
}

template<int DIM>
void FluidSim::OperatorBoundary<DIM>::write(std::vector<scalar>& buf) const
{
  for(const Boundary<DIM>* child : children)
    child->write(buf);
}

template<int DIM>
void FluidSim::OperatorBoundary<DIM>::read(const scalar* data)
{
  int offset = 0;
  for(Boundary<DIM>* child : children) {
    child->read(data + offset);
    offset += child->size() / sizeof(scalar);
  }
}

template<int DIM>
size_t FluidSim::OperatorBoundary<DIM>::size() const
{
  int sum = 0;
  for(const Boundary<DIM>* child : children)
    sum += child->size();

  return sum;
}

template<int DIM>
Particle<DIM>::Particle(const Vectors<DIM>& x_, const Vectors<DIM>& v_, const scalar& radii_, ParticleType type_)
: x(x_), v(v_), radii(radii_), type(type_), deceased(false), fresh(1.0), edge_idx(-1), edge_alpha(0.0), pressure(0.0)
{
  c.setZero();
  buf0.setZero();
  
  bridges.clear();

}

template<int DIM>
Particle<DIM>::Particle(const Vectors<DIM>& x_, const Vectors<DIM>& v_, const scalar& radii_, ParticleType type_, int edge_idx_, const scalar& edge_alpha_)
: x(x_), v(v_), radii(radii_), type(type_), deceased(false), fresh(1.0), edge_idx(edge_idx_), edge_alpha(edge_alpha_), pressure(0.0)
{
  c.setZero();
  buf0.setZero();
  
  bridges.clear();
  
}

template<int DIM>
Particle<DIM>::Particle()
: x(Vectors<DIM>::Zero()), v(Vectors<DIM>::Zero()), radii(0.0), fresh(1.0), deceased(false), edge_idx(-1), edge_alpha(0.0), pressure(0.0), type(PT_LIQUID)
{
  c.setZero();
  buf0.setZero();
  
  bridges.clear();

}

template<int DIM>
Particle<DIM>::Particle(const Particle<DIM>& p)
: x(p.x), v(p.v), buf0(p.buf0), c(p.c), bridges(p.bridges), radii(p.radii), type(p.type), deceased(p.deceased), edge_idx(p.edge_idx), edge_alpha(p.edge_alpha), pressure(p.pressure), fresh(1.0)
{
}

template<int DIM>
void Particle<DIM>::write(std::vector<scalar>& buf) const
{
  buf.push_back(radii);
  for(int i = 0; i < DIM; ++i)
  {
    buf.push_back(x(i));
  }
  for(int i = 0; i < DIM; ++i)
  {
    buf.push_back(v(i));
  }
  for(int i = 0; i < DIM; ++i)
  {
    for(int j = 0; j < DIM; ++j) {
      buf.push_back(c(i,j));
    }
  }
}

template<int DIM>
void Particle<DIM>::read(const scalar* data)
{
  int k = 0;
  radii = data[k++];
  for(int i = 0; i < DIM; ++i)
  {
    x(i) = data[k++];
  }
  for(int i = 0; i < DIM; ++i)
  {
    v(i) = data[k++];
  }
  for(int i = 0; i < DIM; ++i)
  {
    for(int j = 0; j < DIM; ++j) {
      c(i, j) = data[k++];
    }
  }
}

template<int DIM>
size_t Particle<DIM>::size()
{
  size_t s = 0;
  s += sizeof(radii);
  s += DIM * sizeof(scalar);
  s += DIM * sizeof(scalar);
  s += (DIM * DIM) * sizeof(scalar);
  
  return s;
}

// explicit instantiations at bottom
template struct FluidSim::Boundary<2>;
template struct FluidSim::Boundary<3>;

template struct FluidSim::SolidBoundary<2>;
template struct FluidSim::SolidBoundary<3>;

template struct FluidSim::OperatorBoundary<2>;
template struct FluidSim::OperatorBoundary<3>;

template struct FluidSim::SourceBoundary<2>;
template struct FluidSim::SourceBoundary<3>;

template struct Particle<2>;
template struct Particle<3>;

template struct EdgeVelDragIntersection<2>;
template struct EdgeVelDragIntersection<3>;

template struct EdgePhiIntersection<2>;
template struct EdgePhiIntersection<3>;


