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


#include "Force.h"

Force::~Force()
{}

void Force::preCompute( const VectorXs& x, const VectorXs& v, const VectorXs& m, const scalar& dt )
{}

bool Force::isInterHair() const
{
  return false;
}

void Force::storeLambda(const VectorXs& lambda, const VectorXs& lambda_v)
{}

void Force::postStepScene(const scalar& dt )
{}

int Force::getAffectedHair( const std::vector<int> particle_to_hairs )
{
  return -1;
}

bool Force::isExternal()
{
  return false;
}

void Force::setInternalIndex(int index_pos,
                             int index_vel,
                             int index_J,
                             int index_Jv,
                             int index_Jxv,
                             int index_tildeK)
{
  m_internal_index_pos = index_pos;
  m_internal_index_vel = index_vel;
  m_internal_index_J = index_J;
  m_internal_index_Jv = index_Jv;
  m_internal_index_Jxv = index_Jxv;
  m_internal_index_tildeK = index_tildeK;
}
