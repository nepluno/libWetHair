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

if(TARGET tclap::tclap)
    return()
endif()

message(STATUS "Third-party (external): creating target 'tclap::tclap'")

include(FetchContent)
FetchContent_Declare(
    tclap
    GIT_REPOSITORY https://github.com/mirror/tclap.git
    GIT_TAG v1.2.5
    GIT_SHALLOW TRUE
)
FetchContent_GetProperties(tclap)
if(NOT tclap_POPULATED)
    FetchContent_Populate(tclap)
endif()
set(TCLAP_INCLUDE_DIRS ${tclap_SOURCE_DIR}/include)

install(DIRECTORY ${TCLAP_INCLUDE_DIRS}/tclap
    DESTINATION include
)

add_library(tclap_tclap INTERFACE)
add_library(tclap::tclap ALIAS tclap_tclap)

include(GNUInstallDirs)
target_include_directories(tclap_tclap SYSTEM INTERFACE
    $<BUILD_INTERFACE:${TCLAP_INCLUDE_DIRS}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# Install rules
set(CMAKE_INSTALL_DEFAULT_COMPONENT_NAME tclap)
set_target_properties(tclap_tclap PROPERTIES EXPORT_NAME tclap)
install(DIRECTORY ${TCLAP_INCLUDE_DIRS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(TARGETS tclap_tclap EXPORT tclap_Targets)
install(EXPORT tclap_Targets DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/tclap NAMESPACE tclap::)
