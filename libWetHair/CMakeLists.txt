# Find dependencies
find_package (Eigen3 REQUIRED CONFIG)
find_package (TBB REQUIRED CONFIG)
set (THREADS_PREFER_PTHREAD_FLAG ON)
find_package (Threads REQUIRED)

add_library (WetHairCore
  "DER/Dependencies/DegreesOfFreedom.cpp"
  "DER/Dependencies/ElasticStrandUtils.cpp"
  "DER/Dependencies/Kappas.cpp"
  "DER/Dependencies/MaterialFrames.cpp"
  "DER/Dependencies/ReferenceFrames.cpp"
  "DER/Dependencies/Twists.cpp"
  "DER/Forces/BendingForce.cpp"
  "DER/Forces/StretchingForce.cpp"
  "DER/Forces/TwistingForce.cpp"
  "DER/StrandForce.cpp"
  "CTCD.cpp"
  "Capsule.cpp"
  "CohesionTableGen.cpp"
  "CompliantImplicitEuler.cpp"
  "CylindricalShallowFlow.cpp"
  "DragDampingForce.cpp"
  "FluidDragForce.cpp"
  "Force.cpp"
  "GeometricLevelGen.cpp"
  "HairFlow.cpp"
  "Icosphere.cpp"
  "LevelSetForce.cpp"
  "LinearBendingForce.cpp"
  "LinearSpringForce.cpp"
  "MathUtilities.cpp"
  "PolygonalCohesion.cpp"
  "SceneStepper.cpp"
  "SimpleGravityForce.cpp"
  "SpringForce.cpp"
  "StrandCompliantEuler.cpp"
  "StrandCompliantManager.cpp"
  "TimingUtilities.cpp"
  "TwoDScene.cpp"
  "WetHairCore.cpp"
  "fluidsim.cpp"
  "fluidsim2D.cpp"
  "fluidsim3D.cpp"
  "liangbarsky.cpp"
  "sorter.cpp"
  "viscosity3d.cpp"
  "volume_fractions.cpp")
add_library(libWetHair::WetHairCore ALIAS WetHairCore)

target_link_libraries(WetHairCore
  PUBLIC
  Eigen3::Eigen
  TBB::tbb
  PRIVATE
  Threads::Threads)

target_include_directories(WetHairCore PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

install(TARGETS WetHairCore
  EXPORT WetHairCore)

install(EXPORT WetHairCore
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/libWetHair"
  NAMESPACE libWetHair::
  FILE libWetHair_WetHairCore.cmake
  COMPONENT WetHairCore)

install(FILES libWetHairConfig.cmake
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/libWetHair")

# Install headers
install(DIRECTORY .
  DESTINATION include/libWetHair
  FILES_MATCHING
  PATTERN *.h)

