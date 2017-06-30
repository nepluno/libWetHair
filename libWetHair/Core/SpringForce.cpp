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

#include "SpringForce.h"
#include "TwoDScene.h"

template<int DIM>
SpringForce<DIM>::SpringForce( const std::pair<int,int>& endpoints, const scalar& k, const scalar& l0, TwoDScene<DIM>* scene, const scalar& b )
: Force()
, m_endpoints(endpoints)
, m_k(k)
, m_l0(l0)
, m_b(b)
, m_scene( scene )
, m_lambda(0.0)
, m_lambda_v(0.0)
{
  assert( m_endpoints.first >= 0 );
  assert( m_endpoints.second >= 0 );
  assert( m_endpoints.first != m_endpoints.second );
  assert( m_k >= 0.0 );
  assert( m_l0 >= 0.0 );
  assert( m_b >= 0.0 );

}

template<int DIM>
SpringForce<DIM>::~SpringForce()
{}

template<int DIM>
void SpringForce<DIM>::computeIntegrationVars( const VectorXs& x, const VectorXs& v, const VectorXs& m,
                                         VectorXs& lambda, VectorXs& lambda_v,
                                         TripletXs& J, TripletXs& Jv, TripletXs& Jxv, TripletXs& tildeK,
                                         TripletXs& stiffness, TripletXs& damping, VectorXs& Phi, VectorXs& Phiv, const scalar& dt)
{
  Vectors<DIM> nhat = x.segment<DIM>( m_scene->getDof( m_endpoints.second ))-x.segment<DIM>( m_scene->getDof( m_endpoints.first ));
  scalar l = nhat.norm();
  nhat /= l;
  
  Phi(m_internal_index_pos) = l - m_l0;
  stiffness[m_internal_index_pos] = Triplets(m_internal_index_pos, m_internal_index_pos, m_k);
  lambda[m_internal_index_pos] = m_lambda;
  
  for(int r = 0; r < DIM; ++r)
  {
    J[m_internal_index_J + r] = Triplets(m_internal_index_pos, m_scene->getDof(m_endpoints.second) + r, nhat(r));
    J[m_internal_index_J + DIM + r] = Triplets(m_internal_index_pos, m_scene->getDof(m_endpoints.first) + r, -nhat(r));
  }
  
  Matrixs<DIM> dJdx = (Matrixs<DIM>::Identity() - nhat * nhat.transpose()) / l;
  
  scalar weight = -m_lambda;
  
  if(m_b > 0.0) weight += -m_lambda_v;
  
  for(int r = 0; r < DIM; ++r) {
    for(int s = 0; s < DIM; ++s) {

      tildeK[m_internal_index_tildeK + r * (DIM * 2) + s] = Triplets( m_scene->getDof(m_endpoints.first) + r, m_scene->getDof(m_endpoints.first) + s, dJdx(r, s) * weight);
      tildeK[m_internal_index_tildeK + (r + DIM) * (DIM * 2) + DIM + s] = Triplets( m_scene->getDof(m_endpoints.second) + r, m_scene->getDof(m_endpoints.second) + s, dJdx(r, s) * weight);
      
      tildeK[m_internal_index_tildeK + r * (DIM * 2) + DIM + s] = Triplets( m_scene->getDof(m_endpoints.first) + r, m_scene->getDof(m_endpoints.second) + s, -dJdx(r, s) * weight);
      tildeK[m_internal_index_tildeK + (r + DIM) * (DIM * 2) + s] = Triplets( m_scene->getDof(m_endpoints.second) + r, m_scene->getDof(m_endpoints.first) + s, -dJdx(r, s) * weight);
    }
  }
  
  if(m_b > 0.0) {
    Vectors<DIM> dv = v.segment<DIM>( m_scene->getDof(m_endpoints.second) ) - v.segment<DIM>( m_scene->getDof(m_endpoints.first) );
    
    Phiv(m_internal_index_vel) = nhat.dot(dv);
    lambda_v[m_internal_index_vel] = m_lambda_v;
    damping[m_internal_index_vel] = Triplets(m_internal_index_vel, m_internal_index_vel, m_b);
    
    Vectors<DIM> Jxv_base = dJdx * dv;
    
    for(int r = 0; r < DIM; ++r)
    {
      Jv[m_internal_index_Jv + r] = Triplets(m_internal_index_vel, m_scene->getDof( m_endpoints.second ) + r, nhat(r));
      Jv[m_internal_index_Jv + DIM + r] = Triplets(m_internal_index_vel, m_scene->getDof(m_endpoints.first) + r, -nhat(r));
      
      Jxv[m_internal_index_Jxv + r] = Triplets(m_internal_index_vel, m_scene->getDof(m_endpoints.second) + r, Jxv_base(r));
      Jxv[m_internal_index_Jxv + DIM + r] = Triplets(m_internal_index_vel, m_scene->getDof(m_endpoints.first) + r, -Jxv_base(r));
    }
  }
}

template<int DIM>
int SpringForce<DIM>::numJ()
{
  return DIM * 2;   // TODO: Implementation
}

template<int DIM>
int SpringForce<DIM>::numJv()
{
  return m_b > 0.0 ? DIM * 2 : 0;   // TODO: Implementation
}

template<int DIM>
int SpringForce<DIM>::numJxv()
{
  return m_b > 0.0 ? DIM * 2 : 0;   // TODO: Implementation
}

template<int DIM>
int SpringForce<DIM>::numTildeK()
{
  return (DIM * 2) * (DIM * 2);   // TODO: Implementation
}

template<int DIM>
bool SpringForce<DIM>::isParallelized()
{
  return false;   // TODO: Implementation
}

template<int DIM>
bool SpringForce<DIM>::isPrecomputationParallelized()
{
  return false;
}

template<int DIM>
int SpringForce<DIM>::numConstraintPos()
{
  return 1;
}

template<int DIM>
int SpringForce<DIM>::numConstraintVel()
{
  return m_b > 0.0 ? 1 : 0;
}

template<int DIM>
void SpringForce<DIM>::storeLambda(const VectorXs& lambda, const VectorXs& lambda_v)
{
  m_lambda = lambda(m_internal_index_pos);
  
  if(m_b > 0.0) m_lambda_v = lambda_v(m_internal_index_vel);
}

template<int DIM>
const char* SpringForce<DIM>::name()
{
  return "springforce";
}

template<int DIM>
Force* SpringForce<DIM>::createNewCopy()
{
  return new SpringForce(*this);
}

template<int DIM>
void SpringForce<DIM>::getAffectedVars( int pidx, std::unordered_set<int>& vars )
{
  int idir = m_scene->getComponent(pidx);
  if( idir == DIM ) return;
  int ip = m_scene->getVertFromDof(pidx);

  if(ip == m_endpoints.first || ip == m_endpoints.second) {
    for(int r = 0; r < DIM; ++r) {
      vars.insert( m_scene->getDof( m_endpoints.first ) + r);
      vars.insert( m_scene->getDof( m_endpoints.second ) + r);
    }
  }
}

template<int DIM>
int SpringForce<DIM>::getAffectedHair( const std::vector<int> particle_to_hairs )
{
  return particle_to_hairs[m_endpoints.first];
}

template<int DIM>
bool SpringForce<DIM>::isContained( int pidx )
{
  int idir = m_scene->getComponent(pidx);
  if( idir == DIM ) return false;
  int ip = m_scene->getVertFromDof(pidx);

  return (ip == m_endpoints.first || ip == m_endpoints.second);
}

// explicit instantiations at bottom
template class SpringForce<2>;
template class SpringForce<3>;
