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

if(TARGET anttweakbar::anttweakbar)
    return()
endif()

message(STATUS "Third-party (external): creating target 'anttweakbar::anttweakbar'")

include(FetchContent)
FetchContent_Declare(
    anttweakbar
    GIT_REPOSITORY https://github.com/tschw/AntTweakBar.git
    GIT_TAG        4ff73b96515164cbe66b484022e56e771ea6c224
)

set(CMAKE_INSTALL_DEFAULT_COMPONENT_NAME "anttweakbar")
set(ATB_BUILD_EXAMPLES OFF CACHE BOOL "")
FetchContent_MakeAvailable(anttweakbar)

add_library(AntTweakBar::AntTweakBar ALIAS AntTweakBar)

set_target_properties(AntTweakBar PROPERTIES FOLDER third_party)
