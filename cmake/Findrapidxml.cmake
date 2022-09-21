#
# This file is part of the libWetHair open source project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
# Copyright 2017 Breannan Smith (smith@cs.columbia.edu)
# Copyright 2017 Yun (Raymond) Fei
# Copyright 2022 Fredrik Salomonsson
#

find_path(RAPIDXML_INCLUDE_DIR rapidxml.hpp
    PATHS
    ${RAPIDXML_ROOT}
    PATH_SUFFIXES
    include
    include/rapidxml)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(rapidxml DEFAULT_MSG RAPIDXML_INCLUDE_DIR)
mark_as_advanced(RAPIDXML_INCLUDE_DIR)

add_library(rapidxml::rapidxml INTERFACE IMPORTED)
target_include_directories(rapidxml::rapidxml INTERFACE "${RAPIDXML_INCLUDE_DIR}")
