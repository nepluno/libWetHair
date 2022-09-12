#
# Copyright 2020 Adobe. All rights reserved.
# This file is licensed to you under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License. You may obtain a copy
# of the License at http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under
# the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR REPRESENTATIONS
# OF ANY KIND, either express or implied. See the License for the specific language
# governing permissions and limitations under the License.
#
if(TARGET glfw::glfw)
    return()
endif()

if(EMSCRIPTEN)
    # Use the glfw implementation provided by the Emscripten SDK.
    message(STATUS "Third-party (built-in): creating target 'glfw::glfw'")
    add_library(glfw INTERFACE)
    add_library(glfw::glfw ALIAS glfw)
    target_link_options(glfw INTERFACE
        "SHELL:-s MAX_WEBGL_VERSION=2"
        "SHELL:-s USE_GLFW=3"
        -lglfw
    )
    return()
endif()

message(STATUS "Third-party (external): creating target 'glfw::glfw'")

include(FetchContent)
FetchContent_Declare(
    glfw
    GIT_REPOSITORY https://github.com/glfw/glfw.git
    GIT_TAG tags/3.3.6
    GIT_SHALLOW TRUE
)

option(GLFW_BUILD_EXAMPLES "Build the GLFW example programs" OFF)
option(GLFW_BUILD_TESTS "Build the GLFW test programs" OFF)
option(GLFW_BUILD_DOCS "Build the GLFW documentation" OFF)
option(GLFW_INSTALL "Generate installation target" OFF)
option(GLFW_VULKAN_STATIC "Use the Vulkan loader statically linked into application" OFF)
FetchContent_MakeAvailable(glfw)

add_library(glfw::glfw ALIAS glfw)

foreach(name IN ITEMS glfw update_mappings)
    if(TARGET ${name})
        set_target_properties(${name} PROPERTIES FOLDER third_party)
    endif()
endforeach()

# Warning config
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
    target_compile_options(glfw PRIVATE
        "-Wno-missing-field-initializers"
        "-Wno-objc-multiple-method-names"
    )
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    target_compile_options(glfw PRIVATE
        "-Wno-missing-field-initializers"
        "-Wno-objc-multiple-method-names"
    )
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    target_compile_options(glfw PRIVATE
        "-Wno-missing-field-initializers"
        "-Wno-sign-compare"
        "-Wno-unused-parameter"
    )
endif()
