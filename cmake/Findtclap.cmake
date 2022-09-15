#
# This file is part of the libWetHair open source project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
# Copyright 2017 Yun (Raymond) Fei
# Copyright 2022 Fredrik Salomonsson
#
find_package(PkgConfig REQUIRED)
pkg_check_modules(tclap REQUIRED IMPORTED_TARGET tclap)
add_library(tclap::tclap ALIAS PkgConfig::tclap)
