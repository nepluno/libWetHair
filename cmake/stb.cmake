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
if(TARGET stb::stb)
    return()
endif()

message(STATUS "Third-party (external): creating target 'stb::stb'")

include(FetchContent)
FetchContent_Declare(
    stb
    GIT_REPOSITORY https://github.com/nothings/stb.git
    GIT_TAG f67165c2bb2af3060ecae7d20d6f731173485ad0
)
FetchContent_MakeAvailable(stb)

# Generate implementation file
file(WRITE "${stb_BINARY_DIR}/stb_image.cpp.in" [[
    #define STB_IMAGE_IMPLEMENTATION
    #include <stb_image.h>

    #define STB_IMAGE_WRITE_IMPLEMENTATION
    #include <stb_image_write.h>
]])

configure_file(${stb_BINARY_DIR}/stb_image.cpp.in ${stb_BINARY_DIR}/stb_image.cpp COPYONLY)

# Define stb library
add_library(stb ${stb_BINARY_DIR}/stb_image.cpp)
add_library(stb::stb ALIAS stb)

target_include_directories(stb PUBLIC "${stb_SOURCE_DIR}")

set_target_properties(stb PROPERTIES FOLDER third_party)
