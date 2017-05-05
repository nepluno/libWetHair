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

#include "DragDampingForce.h"
#include "TwoDScene.h"
#include "HairFlow.h"
#include <tbb/tbb.h>

template<int DIM>
DragDampingForce<DIM>::DragDampingForce( const TwoDScene<DIM>& scene, const scalar& b, int hidx )
: Force()
, m_b(b)
, m_dofs(scene.getFilmFlows()[hidx]->size() * DIM)
, m_hidx(hidx)
, m_scene(scene)
{
  assert( m_b >= 0.0 );

  m_lambda_v.resize(m_dofs);
  m_lambda_v.setZero();
}

template<int DIM>
DragDampingForce<DIM>::~DragDampingForce()
{}

template<int DIM>
void DragDampingForce<DIM>::addEnergyToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, scalar& E )
{
  assert( x.size() == v.size() );
  assert( x.size() == m.size() );

  std::cerr << outputmod::startred << "WARNING IN DRAGDAMPINGFORCE: " << outputmod::endred << "No energy defined for DragDampingForce." << std::endl;
}

template<int DIM>
void DragDampingForce<DIM>::addGradEToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, VectorXs& gradE )
{
  assert( x.size() == v.size() );
  assert( x.size() == m.size() );
  assert( x.size() == gradE.size() );
  
  const HairFlow<DIM>* flow = m_scene.getFilmFlows()[m_hidx];
  
  int np = flow->size();
  const std::vector<int>& particles = flow->getParticleIndices();
  
  for( int i = 0; i < np; ++i ) {
    gradE.segment<DIM>( m_scene.getDof( particles[i] ) ) += m_b * v.segment<DIM>( m_scene.getDof( particles[i] ) );
  }
}

template<int DIM>
void DragDampingForce<DIM>::addHessXToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, TripletXs& hessE )
{
  assert( x.size() == v.size() );
  assert( x.size() == m.size() );
  // Nothing to do.
}

template<int DIM>
void DragDampingForce<DIM>::addHessVToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, TripletXs& hessE )
{
  assert( x.size() == v.size() );
  assert( x.size() == m.size() );
  
  const HairFlow<DIM>* flow = m_scene.getFilmFlows()[m_hidx];
  
  int np = flow->size();
  const std::vector<int>& particles = flow->getParticleIndices();
  
  for(int i = 0; i < np; ++i) {
    int pidx = particles[i];
    for(int r = 0; r < DIM; ++r)
      hessE.push_back(Triplets(m_scene.getDof( pidx ) + r, m_scene.getDof( pidx ) + r, m_b));
  }
  
}

template<int DIM>
void DragDampingForce<DIM>::addGradEToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, VectorXs& gradE, int pidx )
{
  assert( x.size() == v.size() );
  assert( x.size() == m.size() );
  assert( x.size() == gradE.size() );

  if( m_scene.getComponent(pidx) == DIM ) return;
  gradE(pidx) += m_b * v(pidx);
}

template<int DIM>
void DragDampingForce<DIM>::addHessXToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, VectorXs& hessE, int pidx )
{
  assert( x.size() == v.size() );
  assert( x.size() == m.size() );
  assert( x.size() == hessE.rows() );
  assert( x.size() == hessE.cols() );
  // Nothing to do.
}

template<int DIM>
void DragDampingForce<DIM>::addHessVToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, VectorXs& hessE, int pidx )
{
  assert( x.size() == v.size() );
  assert( x.size() == m.size() );
  assert( x.size() == hessE.rows() );
  assert( x.size() == hessE.cols() );

  if( m_scene.getComponent(pidx) == DIM ) return;
  hessE(pidx) += m_b;
}

template<int DIM>
void DragDampingForce<DIM>::computeIntegrationVars( const VectorXs& x, const VectorXs& v, const VectorXs& m,
                                                   VectorXs& lambda, VectorXs& lambda_v,
                                                   TripletXs& J, TripletXs& Jv, TripletXs& Jxv, TripletXs& tildeK,
                                                   TripletXs& stiffness, TripletXs& damping, VectorXs& Phi, VectorXs& Phiv, const scalar& dt)
{
  const HairFlow<DIM>* flow = m_scene.getFilmFlows()[m_hidx];
  
  int np = flow->size();
  const std::vector<int>& particles = flow->getParticleIndices();
  
  tbb::parallel_for(0, np, 1, [&] (int k) {
    int pidx = particles[k];
    for(int r = 0; r < DIM; ++r) {
      int i = m_scene.getDof(pidx) + r;
      Jv[m_internal_index_Jv + k * DIM + r] = Triplets(m_internal_index_vel + k * DIM + r, i, 1.0);
      damping[m_internal_index_vel + k * DIM + r] = Triplets(m_internal_index_vel + k * DIM + r, m_internal_index_vel + k * DIM + r, m_b);
      Phiv[m_internal_index_vel + k * DIM + r] = v(i);
      lambda_v[m_internal_index_vel + k * DIM + r] = m_lambda_v(k * DIM + r);
    }
  });
}

template<int DIM>
int DragDampingForce<DIM>::numJ()
{
  return 0;
}

template<int DIM>
int DragDampingForce<DIM>::numJv()
{
  const HairFlow<DIM>* flow = m_scene.getFilmFlows()[m_hidx];
  
  int np = flow->size();
  
  return np * DIM;
}

template<int DIM>
int DragDampingForce<DIM>::numJxv()
{
  return 0;
}

template<int DIM>
int DragDampingForce<DIM>::numTildeK()
{
  return 0;
}

template<int DIM>
bool DragDampingForce<DIM>::isParallelized()
{
  return false;
}

template<int DIM>
bool DragDampingForce<DIM>::isPrecomputationParallelized()
{
  return true;
}

template<int DIM>
int DragDampingForce<DIM>::numConstraintPos()
{
  return 0;
}

template<int DIM>
int DragDampingForce<DIM>::numConstraintVel()
{
  const HairFlow<DIM>* flow = m_scene.getFilmFlows()[m_hidx];
  
  int np = flow->size();
  
  return np * DIM;
}

template<int DIM>
const char* DragDampingForce<DIM>::name()
{
  return "dragdamping";
}

template<int DIM>
void DragDampingForce<DIM>::storeLambda(const VectorXs& lambda, const VectorXs& lambda_v)
{
  m_lambda_v = lambda_v.segment(m_internal_index_vel, m_dofs);
}

template<int DIM>
Force* DragDampingForce<DIM>::createNewCopy()
{
  return new DragDampingForce(*this);
}

template<int DIM>
void DragDampingForce<DIM>::getAffectedVars( int colidx, std::unordered_set<int>& vars )
{
  if( m_scene.getComponent(colidx) == DIM ) return;
  vars.insert(colidx);
}

template<int DIM>
int DragDampingForce<DIM>::getAffectedHair( const std::vector<int> particle_to_hairs )
{
  return m_hidx;
}

template<int DIM>
bool DragDampingForce<DIM>::isContained( int colidx )
{
  int idir = m_scene.getComponent(colidx);
  if( idir == DIM ) return false;
  int pidx = m_scene.getVertFromDof( colidx );
  const std::vector<int>& particle_hairs = m_scene.getParticleToHairs();
  return particle_hairs[pidx] == m_hidx;
}

// explicit instantiations at bottom
template class DragDampingForce<2>;
template class DragDampingForce<3>;
