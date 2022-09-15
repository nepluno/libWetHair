# Require a fairly modern cmake version
cmake_minimum_required (VERSION 3.14.0)

project (libWetHair)

set(CMAKE_CXX_STANDARD          11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS        OFF)

option(LIBWETHAIR_BUILD_CORE "Enable the core library, it is ON by default. Use BUILD_SHARED_LIBS to control the type of the library" ON)
option(LIBWETHAIR_BUILD_APP "Enable the libWetHair executable, it is ON by default." ON)
option(LIBWETHAIR_INSTALL_ASSETS "Install the assests, it is ON by default." ON)
option(LIBWETHAIR_FIND_DEPENDENCIES "Find dependencies instead of downloading them" OFF)

include(GNUInstallDirs)

if (NOT CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/install")
endif (NOT CMAKE_INSTALL_PREFIX)

# Initialize the build type (Release, Debug, etc)
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE Release CACHE STRING
    "Choose the type of build, options are: Debug Release."
    FORCE)
endif (NOT CMAKE_BUILD_TYPE)

add_definitions (-DCMAKE_BUILD_TYPE="${CMAKE_BUILD_TYPE}")

if (CMAKE_BUILD_TYPE MATCHES Debug)
  add_definitions (-DDEBUG)
endif (CMAKE_BUILD_TYPE MATCHES Debug)

# Add directory with macros
list (APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

# Don't build in the source directory
if ("${CMAKE_BINARY_DIR}" STREQUAL "${CMAKE_SOURCE_DIR}")
  message (SEND_ERROR "Do not build in the source directory.")
  message (FATAL_ERROR "Remove the created \"CMakeCache.txt\" file and the \"CMakeFiles\" directory, then create a build directory and call \"${CMAKE_COMMAND} <path to the sources>\".")
endif ("${CMAKE_BINARY_DIR}" STREQUAL "${CMAKE_SOURCE_DIR}")

if (LIBWETHAIR_BUILD_CORE)
  add_subdirectory (libWetHair)
endif (LIBWETHAIR_BUILD_CORE)

if (LIBWETHAIR_BUILD_APP)
  add_subdirectory (App)
endif (LIBWETHAIR_BUILD_APP)

if (LIBWETHAIR_INSTALL_ASSETS)
  install(DIRECTORY ${CMAKE_SOURCE_DIR}/assets
    DESTINATION ${CMAKE_INSTALL_DATADIR}/libWetHair)
endif (LIBWETHAIR_INSTALL_ASSETS)
