

if (NOT TARGET libWetHair::WetHairCore)
  # Find dependencies for WetHairCore
  include (CMakeFindDependencyMacro)
  find_dependency (Eigen3)
  find_dependency (TBB)
  find_dependency (Threads)

  include (${CMAKE_CURRENT_LIST_DIR}/libWetHair_WetHairCore.cmake)
endif ()
