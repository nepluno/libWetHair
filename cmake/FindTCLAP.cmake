# - Find TCLAP
# Find the TCLAP headers
#
# TCLAP_INCLUDE_DIR - where to find the TCLAP headers
# TCLAP_FOUND       - True if TCLAP is found

FIND_PACKAGE(PkgConfig QUIET)

if (PKG_CONFIG_FOUND)
  pkg_check_modules (TCLAP REQUIRED IMPORTED_TARGET tclap)
  add_library (TCLAP::tclap ALIAS PkgConfig::TCLAP)
else (PKG_CONFIG_FOUND)
  if (TCLAP_INCLUDE_DIR)
    # already in cache, be silent
    set (TCLAP_FIND_QUIETLY TRUE)
  endif (TCLAP_INCLUDE_DIR)

  # find the headers
  find_path (TCLAP_INCLUDE_PATH tclap/CmdLine.h
    PATHS
    ${TCLAP_ROOT}/include
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_INSTALL_PREFIX}/include
    )

  # handle the QUIETLY and REQUIRED arguments and set TCLAP_FOUND to
  # TRUE if all listed variables are TRUE
  include (FindPackageHandleStandardArgs)
  find_package_handle_standard_args (TCLAP "TCLAP (http://tclap.sourceforge.net/) could not be found. Set TCLAP_INC_DIR to point to the headers." TCLAP_INCLUDE_PATH)

  if (TCLAP_FOUND)
    set (TCLAP_INCLUDE_DIR ${TCLAP_INCLUDE_PATH})
  endif (TCLAP_FOUND)

  mark_as_advanced(TCLAP_INCLUDE_PATH)

  add_library(TCLAP::tclap INTERFACE IMPORTED)
  target_include_directories(TCLAP::tclap INTERFACE "${TCLAP_INCLUDE_PATH}")

endif (PKG_CONFIG_FOUND)
