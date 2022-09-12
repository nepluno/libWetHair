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

if(TARGET rapidxml::rapidxml)
    return()
endif()

message(STATUS "Third-party (external): creating target 'rapidxml::rapidxml'")

include(FetchContent)
FetchContent_Declare(
    rapidxml
    GIT_REPOSITORY https://github.com/discord/rapidxml.git
    GIT_TAG 2ae4b2888165a393dfb6382168825fddf00c27b9
)
FetchContent_GetProperties(rapidxml)
if(NOT rapidxml_POPULATED)
    FetchContent_Populate(rapidxml)
endif()
set(RAPIDXML_INCLUDE_DIRS ${rapidxml_SOURCE_DIR})

install(FILES 
    ${RAPIDXML_INCLUDE_DIRS}/rapidxml.hpp
    ${RAPIDXML_INCLUDE_DIRS}/rapidxml_iterators.hpp
    ${RAPIDXML_INCLUDE_DIRS}/rapidxml_print.hpp
    ${RAPIDXML_INCLUDE_DIRS}/rapidxml_utils.hpp
    DESTINATION include
)

add_library(rapidxml_rapidxml INTERFACE)
add_library(rapidxml::rapidxml ALIAS rapidxml_rapidxml)

include(GNUInstallDirs)
target_include_directories(rapidxml_rapidxml SYSTEM INTERFACE
    $<BUILD_INTERFACE:${RAPIDXML_INCLUDE_DIRS}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# Install rules
set(CMAKE_INSTALL_DEFAULT_COMPONENT_NAME rapidxml)
set_target_properties(rapidxml_rapidxml PROPERTIES EXPORT_NAME rapidxml)
install(DIRECTORY ${RAPIDXML_INCLUDE_DIRS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(TARGETS rapidxml_rapidxml EXPORT rapidxml_Targets)
install(EXPORT rapidxml_Targets DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/rapidxml NAMESPACE rapidxml::)
