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

#include "LinearSpringForce.h"
#include "TwoDScene.h"

template<int DIM>
LinearSpringForce<DIM>::LinearSpringForce( TwoDScene<DIM>* parent, const std::pair<int,int>& endpoints, const scalar& k, const scalar& l0, const scalar& b )
: Force()
, m_endpoints(endpoints)
, m_k(k)
, m_l0(l0)
, m_b(b)
, m_scene( parent )
{
  m_lambda.setZero();
  m_lambda_v.setZero();
  
  assert( m_endpoints.first >= 0 );
  assert( m_endpoints.second >= 0 );
  assert( m_endpoints.first != m_endpoints.second );
  assert( m_k >= 0.0 );
  assert( m_l0 >= 0.0 );
  assert( m_b >= 0.0 );
}

template<int DIM>
LinearSpringForce<DIM>::~LinearSpringForce()
{}

template<int DIM>
void LinearSpringForce<DIM>::preCompute( const VectorXs& x, const VectorXs& v, const VectorXs& m, const scalar& dt )
{
  m_D = (x.segment<DIM>( m_scene->getDof(m_endpoints.first) ) - x.segment<DIM>( m_scene->getDof(m_endpoints.second) ) ).normalized() * m_l0;
}

template<int DIM>
void LinearSpringForce<DIM>::computeIntegrationVars( const VectorXs& x, const VectorXs& v, const VectorXs& m,
                                                    VectorXs& lambda, VectorXs& lambda_v,
                                                    TripletXs& J, TripletXs& Jv, TripletXs& Jxv, TripletXs& tildeK,
                                                    TripletXs& stiffness, TripletXs& damping, VectorXs& Phi, VectorXs& Phiv, const scalar& dt)
{
  Phi.segment<DIM>(m_internal_index_pos) = x.segment<DIM>( m_scene->getDof(m_endpoints.first) ) - x.segment<DIM>( m_scene->getDof(m_endpoints.second) ) - m_D;
  lambda.segment<DIM>(m_internal_index_pos) = m_lambda;
  
  if(m_b > 0.0) {
    Phiv.segment<DIM>(m_internal_index_vel) = v.segment<DIM>( m_scene->getDof(m_endpoints.first) ) - v.segment<DIM>( m_scene->getDof(m_endpoints.second) );
    lambda_v.segment<DIM>(m_internal_index_vel) = m_lambda_v;
  }
  
  for(int r = 0; r < DIM; ++r)
  {
    stiffness[m_internal_index_pos + r] = Triplets(m_internal_index_pos + r, m_internal_index_pos + r, m_k);
    
    J[m_internal_index_J + r] = Triplets(m_internal_index_pos + r, m_scene->getDof(m_endpoints.first) + r, 1.0);
    J[m_internal_index_J + DIM + r] = Triplets(m_internal_index_pos + r, m_scene->getDof(m_endpoints.second) + r, -1.0);
    
    if(m_b > 0.0) {
      damping[m_internal_index_vel + r] = Triplets(m_internal_index_vel + r, m_internal_index_vel + r, m_b);
      
      Jv[m_internal_index_Jv + r] = Triplets(m_internal_index_vel + r, m_scene->getDof(m_endpoints.first) + r, 1.0);
      Jv[m_internal_index_Jv + DIM + r] = Triplets(m_internal_index_vel + r, m_scene->getDof(m_endpoints.second) + r, -1.0);
    }
  }
}

template<int DIM>
int LinearSpringForce<DIM>::numJ()
{
  return DIM * 2;
}

template<int DIM>
int LinearSpringForce<DIM>::numJv()
{
  return m_b > 0.0 ? DIM * 2 : 0;
}

template<int DIM>
int LinearSpringForce<DIM>::numJxv()
{
  return 0;
}

template<int DIM>
int LinearSpringForce<DIM>::numTildeK()
{
  return 0;
}

template<int DIM>
bool LinearSpringForce<DIM>::isParallelized()
{
  return false;
}

template<int DIM>
bool LinearSpringForce<DIM>::isPrecomputationParallelized()
{
  return false;
}

template<int DIM>
int LinearSpringForce<DIM>::numConstraintPos()
{
  return DIM;
}

template<int DIM>
int LinearSpringForce<DIM>::numConstraintVel()
{
  return m_b > 0.0 ? DIM : 0;
}

template<int DIM>
Force* LinearSpringForce<DIM>::createNewCopy()
{
  return new LinearSpringForce(*this);
}

template<int DIM>
const char* LinearSpringForce<DIM>::name()
{
  return "linearspringforce";
}

template<int DIM>
void LinearSpringForce<DIM>::getAffectedVars( int colidx, std::unordered_set<int>& vars )
{
  int idir = m_scene->getComponent( colidx );
  if( idir == DIM ) return;
  int ip = m_scene->getVertFromDof( colidx );

  if(ip == m_endpoints.first || ip == m_endpoints.second) {
    for(int r = 0; r < DIM; ++r) {
      vars.insert( m_scene->getDof(m_endpoints.first) + r);
      vars.insert( m_scene->getDof(m_endpoints.second) + r);
    }
  }
}

template<int DIM>
int LinearSpringForce<DIM>::getAffectedHair( const std::vector<int> particle_to_hairs )
{
  return particle_to_hairs[m_endpoints.first];
}

template<int DIM>
bool LinearSpringForce<DIM>::isContained( int colidx )
{
  int idir = m_scene->getComponent( colidx );
  if( idir == DIM ) return false;
  int ip = m_scene->getVertFromDof( colidx );

  if(ip == m_endpoints.first || ip == m_endpoints.second) return true;
  else return false;
}


template<int DIM>
void LinearSpringForce<DIM>::storeLambda(const VectorXs& lambda, const VectorXs& lambda_v)
{
  m_lambda = lambda.segment<DIM>(m_internal_index_pos);
  if(m_b > 0.0) m_lambda_v = lambda_v.segment<DIM>(m_internal_index_vel);
}

// explicit instantiations at bottom
template class LinearSpringForce<2>;
template class LinearSpringForce<3>;
