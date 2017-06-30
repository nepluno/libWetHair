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

#include "StrandForce.h"

#include "Forces/ForceAccumulator.h"
#include "Forces/BendingForce.h"
#include "Forces/StretchingForce.h"
#include "Forces/TwistingForce.h"
#include "Forces/ViscousOrNotViscous.h"
#include "Dependencies/BendingProducts.h"

// To match with rest of FilmFlow framework sign convention
//  (we compute Forces and Force Jacobians, FilmFlow expects Energy gradients and Hessians)
#define FORCE_SIGN -1.0
#define HESS_SIGN -1.0

#define STRETCH
#define TWIST
#define BEND

StrandForce::StrandForce( 
	TwoDScene<3>* scene,
	const std::vector<int>& consecutiveStrandVertices, 
	const int& parameterIndex, 
	int globalIndex ) :
		m_verts( consecutiveStrandVertices ),
        m_globalIndex( globalIndex ),
        m_strandParams( NULL ),
        m_scene( scene ),
        m_requiresExactForceJacobian( true ),
        m_strandEnergyUpdate( 0. ),
        // m_strandForceUpdate( getNumVertices() * 4 - 1 ),
        m_strandHessianUpdate(),
        m_strandState( NULL ),
        m_startState( NULL )
{
    m_strandParams = m_scene->getStrandParameters( parameterIndex );

    VecX initDoFs( getNumVertices() * 4 - 1 );
    for( int i = 0; i < getNumVertices(); ++i ){
        if( m_scene->isTip( m_verts[i] ) ) initDoFs.segment<3>( i * 4 ) = m_scene->getX().segment<3>( m_scene->getDof( m_verts[i] ) );
        else initDoFs.segment<4>( i * 4 ) = m_scene->getX().segment<4>( m_scene->getDof( m_verts[i] ) );
    }
	m_strandState = new StrandState( initDoFs, m_strandParams->getBendingMatrixBase() );
    m_startState = new StartState( initDoFs );

    resizeInternals();
    freezeRestShape( 0, getNumEdges() ); // for now the rest shape is the shape in which the strand is created, unless this is called later on.
    
// m_hess.resize( getNumVertices() * 4 - 1, getNumVertices() * 4 - 1 );

    // wont need to update first step's DoFs,
    // therefore must initialize the stored quantities as well
// recomputeGlobal();

    m_lambda.resize( numConstraintNonViscous() ); m_lambda.setZero();
    m_lambda_v.resize( numConstraintViscous() ); m_lambda_v.setZero();
}

StrandForce::~StrandForce()
{}

StrandState::StrandState( const VecX& initDoFs, BendingMatrixBase& bendingMatrixBase ):
    m_dofs( initDoFs ),
    m_edges( m_dofs ),
    m_lengths( m_edges ),
    m_tangents( m_edges, m_lengths ),
    m_referenceFrames1( m_tangents ),
    m_referenceFrames2( m_tangents, m_referenceFrames1 ),
    m_referenceTwists( m_tangents, m_referenceFrames1 ),
    m_twists( m_referenceTwists, m_dofs ),
    m_curvatureBinormals( m_tangents ),
    m_trigThetas( m_dofs ),
    m_materialFrames1( m_trigThetas, m_referenceFrames1, m_referenceFrames2 ),
    m_materialFrames2( m_trigThetas, m_referenceFrames1, m_referenceFrames2 ),
    m_kappas( m_curvatureBinormals, m_materialFrames1, m_materialFrames2 ),
    m_gradKappas( m_lengths, m_tangents, m_curvatureBinormals, m_materialFrames1, m_materialFrames2, m_kappas ),
    m_gradTwists( m_lengths, m_curvatureBinormals ),
    m_gradTwistsSquared( m_gradTwists ),
    m_hessKappas( m_lengths, m_tangents, m_curvatureBinormals, m_materialFrames1, m_materialFrames2, m_kappas ),
    m_hessTwists( m_tangents, m_lengths, m_curvatureBinormals ),
    m_bendingProducts( bendingMatrixBase, m_gradKappas )
{}

StartState::StartState( const VecX& initDoFs ):
    m_dofs( initDoFs ),
    m_edges( m_dofs ),
    m_lengths( m_edges ),
    m_tangents( m_edges, m_lengths ),
    m_referenceFrames1( m_tangents ),
    m_referenceFrames2( m_tangents, m_referenceFrames1 ),
    m_referenceTwists( m_tangents, m_referenceFrames1 ),
    m_twists( m_referenceTwists, m_dofs ),
    m_curvatureBinormals( m_tangents ),
    m_trigThetas( m_dofs ),
    m_materialFrames1( m_trigThetas, m_referenceFrames1, m_referenceFrames2 ),
    m_materialFrames2( m_trigThetas, m_referenceFrames1, m_referenceFrames2 ),
    m_kappas( m_curvatureBinormals, m_materialFrames1, m_materialFrames2 )
{}

void StrandForce::updateStartDoFs( const VecX& x_startOfStep )
{
    const VecX& currentStrandDoFs = x_startOfStep.segment( m_scene->getDof( m_verts[0] ), getNumVertices() * 4 - 1 );
    m_startState->m_dofs.set( currentStrandDoFs );
}

void StrandForce::updateRestShape( const VecX& dof_restshape, scalar damping )
{
  StartState restshape_state( dof_restshape );
  
  int nedges = getNumEdges();
  for( IndexType vtx = 0; vtx < nedges; ++vtx )
  { // Fix rest lengths
    m_restLengths[vtx] = ( 1. - damping ) * restshape_state.m_lengths[vtx] + damping * m_restLengths[vtx];
  }
  updateEverythingThatDependsOnRestLengths();
  
  for( IndexType vtx = 0; vtx < nedges; ++vtx )
  {
    m_restKappas[vtx] = ( 1. - damping ) * restshape_state.m_kappas[vtx] + damping * m_restKappas[vtx];
    m_restTwists[vtx] = ( 1. - damping ) * restshape_state.m_twists[vtx] + damping * m_restTwists[vtx];
  }
}

void StrandForce::resizeInternals()
{ // To be called on creation
    m_restLengths.resize( getNumEdges() );
    m_restKappas.resize( getNumEdges() );
    m_restTwists.resize( getNumEdges() );
    m_vertexMasses.resize( getNumVertices() );
    m_VoronoiLengths.resize( getNumVertices() );
    m_invVoronoiLengths.resize( getNumVertices() );
}

void StrandForce::freezeRestShape( unsigned begin, unsigned end, scalar damping )
{ // Take the current configuration as rest shape

    for( IndexType vtx = begin; vtx < end; ++vtx )
    { // Fix rest lengths
        m_restLengths[vtx] = ( 1. - damping ) * m_strandState->m_lengths[vtx] + damping * m_restLengths[vtx];
    }
    updateEverythingThatDependsOnRestLengths();

    for( IndexType vtx = begin; vtx < end; ++vtx )
    {
        m_restKappas[vtx] = ( 1. - damping ) * m_strandState->m_kappas[vtx] + damping * m_restKappas[vtx];
        m_restTwists[vtx] = ( 1. - damping ) * m_strandState->m_twists[vtx] + damping * m_restTwists[vtx];
    }
}

void StrandForce::updateEverythingThatDependsOnRestLengths()
{
    // Total rest length
    m_totalRestLength = 0.0;
    for( IndexType vtx = 0; vtx < getNumEdges(); ++vtx ){
        m_totalRestLength += m_restLengths[vtx];
    }

    // Compute Voronoi lengths
    m_VoronoiLengths[0] = 0.5 * m_restLengths[0];
    for( IndexType vtx = 1; vtx < getNumEdges(); ++vtx ){
        m_VoronoiLengths[vtx] = 0.5 * ( m_restLengths[vtx - 1] + m_restLengths[vtx] );
    }
    m_VoronoiLengths[getNumEdges()] = 0.5 * m_restLengths[getNumVertices() - 2];

    // Compute masses and inverse of Voronoi lengths
    for( IndexType vtx = 0; vtx < getNumVertices(); ++vtx ){
        m_vertexMasses[vtx] = m_strandParams->m_density * m_VoronoiLengths[vtx] * 
                                M_PI * m_strandParams->getRadius( vtx, getNumVertices() ) * m_strandParams->getRadius( vtx, getNumVertices() );
        m_invVoronoiLengths[vtx] = 1.0 / m_VoronoiLengths[vtx];
    }
}

///////////////////////////////////////////// Force functions //////////////////////////////////////////////////////////

void StrandForce::clearStored()
{
    m_strandEnergyUpdate = 0.;
    m_strandForceUpdate.setZero();
    m_strandHessianUpdate.clear();
    m_hess.setZero();
}

void StrandForce::recomputeGlobal()
{
    clearStored();
    accumulateQuantity( m_strandEnergyUpdate );
    accumulateQuantity( m_strandForceUpdate );
    accumulateQuantity( m_strandHessianUpdate );
    m_hess.setFromTriplets( m_strandHessianUpdate.begin(), m_strandHessianUpdate.end() );

    // Free some memory
    m_strandState->m_hessTwists.free();
    m_strandState->m_hessKappas.free();
}

template<typename AccumulatedT>
void StrandForce::accumulateQuantity( AccumulatedT& accumulated )
{
    ForceAccumulator< StretchingForce< NonViscous > >::accumulate( accumulated, *this );
    ForceAccumulator< TwistingForce< NonViscous > >::accumulate( accumulated, *this );
    ForceAccumulator< BendingForce< NonViscous > >::accumulate( accumulated, *this );

    if( m_strandParams->m_accumulateWithViscous ){
        if( !m_strandParams->m_accumulateViscousOnlyForBendingModes )
        {
            ForceAccumulator< StretchingForce< Viscous > >::accumulate( accumulated, *this );
        }
        ForceAccumulator< TwistingForce<Viscous> >::accumulate( accumulated, *this );
        ForceAccumulator< BendingForce<Viscous> >::accumulate( accumulated, *this );
    }
}

Force* StrandForce::createNewCopy()
{
    return new StrandForce(*this);
}

void StrandForce::preCompute( const VecX& x, const VecX& v, const VecX& m, const scalar& dt )
{  
    /* nothing to do here, updateStartDoFs called separately and otherwise need to update every time we compute (in case nonlinear) */
}

void StrandForce::computeIntegrationVars( const VectorXs& x, const VectorXs& v, const VectorXs& m,
                                         VectorXs& lambda, VectorXs& lambda_v,
                                         TripletXs& J, TripletXs& Jv, TripletXs& Jxv, TripletXs& tildeK,
                                         TripletXs& stiffness, TripletXs& damping, VectorXs& Phi, VectorXs& Phiv, const scalar& dt)
{
    // update DoFs: 'future' if Nonlinear, these are current if Linear
    const VecX& futureStrandDoFs = x.segment( m_scene->getDof( m_verts[0] ), getNumVertices() * 4 - 1 );
    if( futureStrandDoFs != m_strandState->m_dofs.get() ){
        m_strandState->m_dofs.set( futureStrandDoFs );
    }  

    int ncnv = numConstraintNonViscous();
    lambda.segment( m_internal_index_pos, ncnv ) = m_lambda;
    if( m_strandParams->m_accumulateWithViscous ) lambda.segment( m_internal_index_pos + ncnv, numConstraintViscous() ) = m_lambda_v;

    VectorXs combinedLambda = m_lambda;
    if( m_strandParams->m_accumulateWithViscous ){
        if( !m_strandParams->m_accumulateViscousOnlyForBendingModes ){
            combinedLambda += m_lambda_v;
        }
        else{
            combinedLambda.segment( getNumEdges(), m_lambda_v.size() ) += m_lambda_v;
        }
    }

    const unsigned global_start_dof = m_scene->getDof( m_verts[0] );

    unsigned poffset = 0;
    unsigned jOffset = 0;
    unsigned koffset = 0;

#ifdef STRETCH
    StretchingForce< NonViscous >::accumulateIntegrationVars( 
                                        m_internal_index_pos, m_internal_index_J, m_internal_index_tildeK,
                                        global_start_dof, *this, combinedLambda, J, tildeK, stiffness, Phi, poffset );
    poffset += getNumEdges();
    jOffset += getNumEdges() * 6;
    koffset += getNumEdges() * 36;
#endif

#ifdef TWIST
    TwistingForce< NonViscous >::accumulateIntegrationVars( 
                                        m_internal_index_pos + poffset,
                                        m_internal_index_J + jOffset,
                                        m_internal_index_tildeK + koffset,
                                        global_start_dof, *this, combinedLambda, J, tildeK, stiffness, Phi, poffset );
    poffset += (getNumVertices() - 2);
    jOffset += (getNumVertices() - 2) * 11;
    koffset += (getNumVertices() - 2) * 121;
#endif

#ifdef BEND
    BendingForce< NonViscous >::accumulateIntegrationVars( 
                                        m_internal_index_pos + poffset,
                                        m_internal_index_J + jOffset,
                                        m_internal_index_tildeK + koffset,
                                        global_start_dof, *this, combinedLambda, J, tildeK, stiffness, Phi, poffset );

    poffset += 2 * (getNumVertices() - 2);
    jOffset += 2 * (getNumVertices() - 2) * 11;
    koffset += 2 * (getNumVertices() - 2) * 121;
#endif

    if( m_strandParams->m_accumulateWithViscous ){
        if( !m_strandParams->m_accumulateViscousOnlyForBendingModes ){
#ifdef STRETCH            
            StretchingForce< Viscous >::accumulateIntegrationVars( 
                                        m_internal_index_pos + poffset, 
                                        m_internal_index_J + jOffset, 
                                        m_internal_index_tildeK + koffset,
                                        global_start_dof, *this, combinedLambda, J, tildeK, stiffness, Phi, poffset );
    poffset += getNumEdges();
    jOffset += getNumEdges() * 6;
    koffset += getNumEdges() * 36;
#endif            
        }
#ifdef TWIST        
        TwistingForce< Viscous >::accumulateIntegrationVars( 
                                        m_internal_index_pos + poffset, 
                                        m_internal_index_J + jOffset, 
                                        m_internal_index_tildeK + koffset,
                                        global_start_dof, *this, combinedLambda, J, tildeK, stiffness, Phi, poffset );
    poffset += (getNumVertices() - 2);
    jOffset += (getNumVertices() - 2) * 11;
    koffset += (getNumVertices() - 2) * 121;
#endif        

#ifdef BEND
        BendingForce< Viscous >::accumulateIntegrationVars( 
                                        m_internal_index_pos + poffset, 
                                        m_internal_index_J + jOffset, 
                                        m_internal_index_tildeK + koffset,
                                        global_start_dof, *this, combinedLambda, J, tildeK, stiffness, Phi, poffset );
#endif
    }
}

int StrandForce::numJ()
{
    int numJ = 0;

#ifdef STRETCH
    numJ += (getNumEdges() * 6); // Stretch
#endif
#ifdef TWIST
    numJ += ((getNumVertices() - 2 ) * 11); // Twist
#endif
#ifdef BEND
    numJ += 2 * ((getNumVertices() - 2 ) * 11); // Bend
#endif

    if( m_strandParams->m_accumulateWithViscous )
    { // double the active forces for viscous
        if( !m_strandParams->m_accumulateViscousOnlyForBendingModes ){
#ifdef STRETCH
            numJ += (getNumEdges() * 6); // Stretch
#endif
        }
#ifdef TWIST
        numJ += ((getNumVertices() - 2 ) * 11); // Twist
#endif
#ifdef BEND
        numJ += 2 * ((getNumVertices() - 2 ) * 11); // Bend
#endif
    }

    return numJ;
}

int StrandForce::numJv()
{ // viscous now treated via positions
    return 0;
}

int StrandForce::numJxv()
{
    return numJv();
}

int StrandForce::numTildeK()
{ //number of nonzeros in hessian, split up by force
    int numTildeK = 0;

#ifdef STRETCH
    numTildeK += getNumEdges() * 36; // Stretch
#endif

#ifdef TWIST
    numTildeK += (getNumVertices() - 2) * 121; // Twist
#endif

#ifdef BEND
    numTildeK += 2 * (getNumVertices() - 2) * 121; // Bend
#endif

    return numTildeK;
}

bool StrandForce::isParallelized()
{
    return false; // TODO, could PARALLELIZE each force calls
}

bool StrandForce::isPrecomputationParallelized()
{
  return false;
}

int StrandForce::numConstraintPos()
{
    int numConstraintPos = numConstraintNonViscous();
    numConstraintPos += numConstraintViscous();
    return numConstraintPos;
}

int StrandForce::numConstraintVel()
{ // now treated via pos
    return 0;
}

int StrandForce::numConstraintNonViscous()
{ //Spring = numEdges  //Twist = NumVertices - 2  //Bending = 2*(NumVertices - 2)
    int numConstraintNonViscous = 0;
#ifdef STRETCH
    numConstraintNonViscous += getNumEdges(); // Stretch
#endif

#ifdef TWIST
    numConstraintNonViscous += (getNumVertices() - 2); // Twist
#endif

#ifdef BEND
    numConstraintNonViscous += 2 * (getNumVertices() - 2); // Bend
#endif
    return numConstraintNonViscous;
}

int StrandForce::numConstraintViscous()
{
    int numConstraintViscous = 0;
    if( m_strandParams->m_accumulateWithViscous ){
        if( !m_strandParams->m_accumulateViscousOnlyForBendingModes ){
#ifdef STRETCH
            numConstraintViscous += getNumEdges(); // Stretch
#endif
        }
#ifdef TWIST
        numConstraintViscous += (getNumVertices() - 2); // Twist
#endif
#ifdef BEND
        numConstraintViscous += 2 * (getNumVertices() - 2); // Bend
#endif
    }
    return numConstraintViscous;
}

void StrandForce::storeLambda( const VectorXs& lambda, const VectorXs& lambda_v )
{
    int ncnv = numConstraintNonViscous();
    m_lambda = lambda.segment( m_internal_index_pos, ncnv );
    if( m_strandParams->m_accumulateWithViscous ){
        m_lambda_v = lambda.segment( m_internal_index_pos + ncnv, numConstraintViscous() );
    }
}

void StrandForce::getAffectedVars( int colidx, std::unordered_set<int>& vars )
{
    int ip = m_scene->getVertFromDof( colidx );
    for( int v = 0; v < getNumVertices(); ++v ){
        if( m_verts[v] == ip )
        {
            if( m_scene->isTip( ip ) ){
                for( int r = 0; r < 3; ++r ) vars.insert( m_scene->getDof( ip ) + r );
            }
            else{
                for( int r = 0; r < 4; ++r ) vars.insert( m_scene->getDof( ip ) + r );
            }

            // include previous and next vertices affected by this vert's stretch/bend/twist
            if( v != 0 ){
                for( int r = 0; r < 4; ++r ) vars.insert( m_scene->getDof( m_verts[v - 1] ) + r );
            }
            if( v != getNumVertices() - 1 ){
                for( int r = 0; r < 3; ++r ) vars.insert( m_scene->getDof( m_verts[v + 1] ) + r );
            }
            break;
        }
    }
}

void StrandForce::getLocalAffectedVars( int colidx, std::vector< std::pair<int,int> >& vars )
{ // local, global
    int ip = m_scene->getVertFromDof( colidx );
    for( int v = 0; v < getNumVertices(); ++v ){
        if( m_verts[v] == ip )
        {
            if( m_scene->isTip( ip ) ){
                for( int r = 0; r < 3; ++r ) vars.push_back( std::pair<int,int>( 4 * v + r, m_scene->getDof( ip ) + r ) );
            }
            else{
                for( int r = 0; r < 4; ++r ) vars.push_back( std::pair<int,int>( 4 * v + r, m_scene->getDof( ip ) + r ) );
            }

            // include previous and next vertices affected by this vert's stretch/bend/twist
            if( v != 0 ){
                for( int r = 0; r < 4; ++r ) vars.push_back( std::pair<int,int>( 4 * (v - 1) + r, m_scene->getDof( m_verts[v - 1] ) + r ) );
            }
            if( v != getNumVertices() - 1 ){
                for( int r = 0; r < 3; ++r ) vars.push_back( std::pair<int,int>( 4 * (v + 1) + r, m_scene->getDof( m_verts[v + 1] ) + r ) );
            }
            break;
        }
    }
}

int StrandForce::getAffectedHair( const std::vector<int> particle_to_hairs )
{
  // TODO:: return the hair index
    return particle_to_hairs[ m_verts[0] ];
  // return -1;
}

bool StrandForce::isContained( int pidx )
{
    int ip = m_scene->getVertFromDof( pidx );
    for( int v = 0; v < getNumVertices(); ++v ){
        if( m_verts[v] == ip ){
            return true;
        }
    }
    return false;
}

const char* StrandForce::name(){ return "Strand Material Forces"; }

