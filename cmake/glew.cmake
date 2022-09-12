# 
# This file is part of the libWetHair open source project
# 
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
# 
# Copyright 2022 Yun (Raymond) Fei, Henrique Teles Maia, Christopher Batty,
# Changxi Zheng, and Eitan Grinspun
# 

if(TARGET glew::glew)
    return()
endif()

message(STATUS "Third-party (external): creating target 'glew::glew'")

include(FetchContent)
FetchContent_Declare(
    glew
    GIT_REPOSITORY https://github.com/Perlmint/glew-cmake.git
    GIT_TAG        f456deace7b408655109aaeff71421ef2d3858c6
)

set(CMAKE_INSTALL_DEFAULT_COMPONENT_NAME "glew")
set(glew-cmake_BUILD_SHARED OFF CACHE BOOL "")
FetchContent_MakeAvailable(glew)

add_library(glew::glew ALIAS libglew_static)

set_target_properties(libglew_static PROPERTIES FOLDER third_party)