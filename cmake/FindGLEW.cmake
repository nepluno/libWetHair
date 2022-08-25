# Copyright (c) 2009 Boudewijn Rempt <boud@valdyas.org>
# Copyright (c) 2022 Fredrik Salomonsson <fredriks@d2.com>
#
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying COPYING-CMAKE-SCRIPTS file.
#
# - try to find glew library and include files
#  GLEW_INCLUDE_DIR, where to find GL/glew.h, etc.
#  GLEW_LIBRARIES, the libraries to link against
#  GLEW_FOUND, If false, do not try to use GLEW.
# Also defined, but not for general use are:
#  GLEW_GLEW_LIBRARY = the full path to the glew library.

IF (WIN32)

  IF(CYGWIN)

    FIND_PATH( GLEW_INCLUDE_DIR GL/glew.h)

    FIND_LIBRARY( GLEW_GLEW_LIBRARY glew32
      ${OPENGL_LIBRARY_DIR}
      /usr/lib/w32api
      /usr/X11R6/lib
    )


  ELSE(CYGWIN)

    FIND_PATH( GLEW_INCLUDE_DIR GL/glew.h
      $ENV{GLEW_ROOT_PATH}/include
    )

    FIND_LIBRARY( GLEW_GLEW_LIBRARY
      NAMES glew glew32
      PATHS
      $ENV{GLEW_ROOT_PATH}/lib
      ${OPENGL_LIBRARY_DIR}
    )

  ENDIF(CYGWIN)

ELSE (WIN32)

  IF (APPLE)
# These values for Apple could probably do with improvement.
    FIND_PATH( GLEW_INCLUDE_DIR GL/glew.h
      /opt/homebrew/opt/glew/include
    )

    FIND_LIBRARY( GLEW_GLEW_LIBRARY GLEW
      /opt/homebrew/opt/glew/lib
    )
  ELSE (APPLE)
    FIND_PACKAGE(PkgConfig REQUIRED)
    PKG_CHECK_MODULES(GLEW REQUIRED IMPORTED_TARGET glew)
    ADD_LIBRARY(GLEW::glew ALIAS PkgConfig::GLEW)
  ENDIF (APPLE)

ENDIF (WIN32)

IF(NOT TARGET GLEW::glew)
  SET( GLEW_FOUND "NO" )
  IF(GLEW_INCLUDE_DIR)
    IF(GLEW_GLEW_LIBRARY)
      # Is -lXi and -lXmu required on all platforms that have it?
      # If not, we need some way to figure out what platform we are on.
      SET( GLEW_LIBRARIES
        ${GLEW_GLEW_LIBRARY}
        ${GLEW_cocoa_LIBRARY}
        )
      SET( GLEW_FOUND "YES" )

      #The following deprecated settings are for backwards compatibility with CMake1.4
      SET (GLEW_LIBRARY ${GLEW_LIBRARIES})
      SET (GLEW_INCLUDE_PATH ${GLEW_INCLUDE_DIR})

    ENDIF(GLEW_GLEW_LIBRARY)
  ENDIF(GLEW_INCLUDE_DIR)

  IF(GLEW_FOUND)
    IF(NOT GLEW_FIND_QUIETLY)
      MESSAGE(STATUS "Found Glew: ${GLEW_LIBRARIES}")
    ENDIF(NOT GLEW_FIND_QUIETLY)
  ELSE(GLEW_FOUND)
    IF(GLEW_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR "Could not find Glew")
    ENDIF(GLEW_FIND_REQUIRED)
  ENDIF(GLEW_FOUND)

  ADD_LIBRARY(GLEW::glew INTERFACE IMPORTED)
  TARGET_INCLUDE_DIRECTORIES(GLEW::glew INTERFACE "${GLEW_INCLUDE_DIR}")
  TARGET_LINK_LIBRARIES(GLEW::glew INTERFACE "${GLEW_LIBRARIES}")
ENDIF(NOT TARGET GLEW::glew)

MARK_AS_ADVANCED(
  GLEW_INCLUDE_DIR
  GLEW_GLEW_LIBRARY
  GLEW_Xmu_LIBRARY
  GLEW_Xi_LIBRARY
)
