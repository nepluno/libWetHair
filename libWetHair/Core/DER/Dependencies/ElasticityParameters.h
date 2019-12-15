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


#ifndef ELASTICPARAMETERS_H_
#define ELASTICPARAMETERS_H_

#include "DependencyNode.h"
#define PI_4 0.785398163397448309616


/**
 * Unit: cm
 */
class PhysicalRadius: public DependencyNode< scalar >
{
public:
    PhysicalRadius( scalar radius ) :
            DependencyNode< scalar >( radius )
    {
#ifdef VERBOSE_DEPENDENCY_NODE
        std::cout << "Creating " << name() << ' ' << this << '\n';
#endif

        setClean();
    }

    virtual const char* name() const
    {
        return "PhysicalRadius";
    }

protected    :
    virtual void compute()
    {}
};

/**
 * Unit: no dimension
 */
class BaseRotation: public DependencyNode<scalar>
{
public:
    BaseRotation( scalar baseRotation ) :
            DependencyNode<scalar>( baseRotation )
    {
#ifdef VERBOSE_DEPENDENCY_NODE
        std::cout << "Creating " << name() << ' ' << this << '\n';
#endif

        setClean();
    }

    virtual const char* name() const
    {
        return "BaseRotation";
    }

protected:
    virtual void compute()
    {}
};

/**
 * \brief This contains the bending matrix base, must be multiplied by the appropriate viscous or non-viscous
 * coefficient (with optional interpolation factor).
 *
 * Unit: cm^4
 */
class BendingMatrixBase: public DependencyNode<Mat2>
{
public:
    BendingMatrixBase( PhysicalRadius& rad, BaseRotation& baseRotation ) :
            DependencyNode<Mat2>( Mat2() ), //
            m_physicalRadius( rad ), //
            m_baseRotation( baseRotation )
    {
#ifdef VERBOSE_DEPENDENCY_NODE
        std::cout << "Creating " << name() << ' ' << this << '\n';
#endif
        m_value.setZero();
        m_physicalRadius.addDependent( this );
        m_baseRotation.addDependent( this );
    }

    virtual const char* name() const
    {
        return "BendingMatrixBase";
    }

protected:
    virtual void compute()
    {
        const scalar& radius = m_physicalRadius.get();
        const scalar baseRotation = m_baseRotation.get();

        Mat2& B = m_value;
        B( 0, 0 ) = PI_4 * radius * cube( radius );
        B( 1, 1 ) = PI_4 * radius * cube( radius );
        // rotate cross section by a constant angle
        const Mat2& rot = Eigen::Rotation2D<scalar>( baseRotation ).toRotationMatrix();
        B = rot * B * rot.transpose();
        B( 0, 1 ) = B( 1, 0 ) = 0.5 * ( B( 0, 1 ) + B( 1, 0 ) ); // For perfect numerical symmetry

        setDependentsDirty();
    }

    PhysicalRadius& m_physicalRadius;
    BaseRotation& m_baseRotation;
};

/**
 * Unit: dPa = g cm^-1 s^-2
 */
class YoungsModulus: public DependencyNode<scalar>
{
public:
    YoungsModulus( scalar youngsModulus ) :
            DependencyNode<scalar>( youngsModulus )
    {
        setClean();
    }

    virtual const char* name() const
    {
        return "YoungsModulus";
    }

protected:
    virtual void compute()
    {}
};

/**
 * Unit: dPa = g cm^-1 s^-2
 */
class ShearModulus: public DependencyNode<scalar>
{
public:
    ShearModulus( scalar shearModulus ) :
            DependencyNode<scalar>( shearModulus )
    {
        setClean();
    }

    virtual const char* name() const
    {
        return "ShearModulus";
    }

protected:
    virtual void compute()
    {}
};

/**
 * Unit: 10^-5 N = g cm s^-2
 */
class ElasticKs: public DependencyNode<scalar>
{
public:
    ElasticKs( PhysicalRadius& rad, YoungsModulus& ym ) :
            DependencyNode<scalar>( std::numeric_limits<scalar>::signaling_NaN() ), //
            m_physicalRadius( rad ), //
            m_youngsModulus( ym ) //
    {
        m_physicalRadius.addDependent( this );
        m_youngsModulus.addDependent( this );
    }

    virtual const char* name() const
    {
        return "ElasticKs";
    }

protected:
    virtual void compute()
    {
        const scalar& radius = m_physicalRadius.get();
        const scalar youngsModulus = m_youngsModulus.get();

        m_value = M_PI * radius * radius * youngsModulus;

        setDependentsDirty();
    }

    PhysicalRadius& m_physicalRadius;
    YoungsModulus& m_youngsModulus;
};

/**
 * Unit: 10^-5 cm^2 N = g cm^3 s^-2
 */
class ElasticKt: public DependencyNode<scalar>
{
public:
    ElasticKt( PhysicalRadius& rad, ShearModulus& sm ) :
            DependencyNode<scalar>( std::numeric_limits<scalar>::signaling_NaN() ), //
            m_physicalRadius( rad ), //
            m_shearModulus( sm )
    {
        m_physicalRadius.addDependent( this );
        m_shearModulus.addDependent( this );
    }

    virtual const char* name() const
    {
        return "ElasticKt";
    }

protected:
    virtual void compute()
    {
        const scalar& radius = m_physicalRadius.get();
        const scalar shearModulus = m_shearModulus.get();

        m_value = PI_4 * radius * radius
                * ( radius * radius + radius * radius ) * shearModulus;

        setDependentsDirty();
    }

    PhysicalRadius& m_physicalRadius;
    ShearModulus& m_shearModulus;
};

#endif
