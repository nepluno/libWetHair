//
// This file is part of the libWetHair open source project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright 2008 Christopher Batty, Robert Bridson
//

#ifndef LIBWETHAIR_CORE_VISCOSITY_3D_H_
#define LIBWETHAIR_CORE_VISCOSITY_3D_H_

#include "array3.h"

namespace libwethair {

void advance_viscosity_implicit_weighted(
    Array3s& u, Array3s& v, Array3s& w, Array3s& u_visc_impulse,
    Array3s& v_visc_impulse, Array3s& w_visc_impulse, const Array3s& vol_u,
    const Array3s& vol_v, const Array3s& vol_w, const Array3s& vol_c,
    const Array3s& vol_ex, const Array3s& vol_ey, const Array3s& vol_ez,
    const Array3s& solid_phi, const scalar viscosity, scalar dt, scalar dx);

}  // namespace libwethair

#endif  // LIBWETHAIR_CORE_VISCOSITY_3D_H_
