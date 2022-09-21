#
# This file is part of the libWetHair open source project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
# Copyright 2022 Fredrik Salomonsson
#
find_path(STB_IMAGE_WRITE_INCLUDE_DIR stb_image_write.h
    PATHS
    ${STB_IMAGE_WRITE_ROOT}
    PATH_SUFFIXES
    include)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(stb DEFAULT_MSG STB_IMAGE_WRITE_INCLUDE_DIR)
mark_as_advanced(STB_IMAGE_WRITE_INCLUDE_DIR)

add_library(stb::stb INTERFACE IMPORTED)
target_include_directories(stb::stb INTERFACE "${STB_IMAGE_WRITE_INCLUDE_DIR}")
