# libWetHair Executable

# Find dependencies
find_package (ANTTWEAKBAR REQUIRED)
find_package (Eigen3 REQUIRED CONFIG)
find_package (FreeGLUT REQUIRED CONFIG)
find_package (GLEW REQUIRED)
find_package (OpenGL REQUIRED)
find_package (PNG REQUIRED)
find_package (RapidXML REQUIRED)
find_package (TBB REQUIRED CONFIG)
find_package (TCLAP REQUIRED)
set (THREADS_PREFER_PTHREAD_FLAG ON)
find_package (Threads REQUIRED)

add_library (WetHairCore
  "Core/DER/Dependencies/DegreesOfFreedom.cpp"
  "Core/DER/Dependencies/ElasticStrandUtils.cpp"
  "Core/DER/Dependencies/Kappas.cpp"
  "Core/DER/Dependencies/MaterialFrames.cpp"
  "Core/DER/Dependencies/ReferenceFrames.cpp"
  "Core/DER/Dependencies/Twists.cpp"
  "Core/DER/Forces/BendingForce.cpp"
  "Core/DER/Forces/StretchingForce.cpp"
  "Core/DER/Forces/TwistingForce.cpp"
  "Core/DER/StrandForce.cpp"
  "Core/CTCD.cpp"
  "Core/Capsule.cpp"
  "Core/CohesionTableGen.cpp"
  "Core/CompliantImplicitEuler.cpp"
  "Core/CylindricalShallowFlow.cpp"
  "Core/DragDampingForce.cpp"
  "Core/FluidDragForce.cpp"
  "Core/Force.cpp"
  "Core/GeometricLevelGen.cpp"
  "Core/HairFlow.cpp"
  "Core/Icosphere.cpp"
  "Core/LevelSetForce.cpp"
  "Core/LinearBendingForce.cpp"
  "Core/LinearSpringForce.cpp"
  "Core/MathUtilities.cpp"
  "Core/PolygonalCohesion.cpp"
  "Core/SceneStepper.cpp"
  "Core/SimpleGravityForce.cpp"
  "Core/SpringForce.cpp"
  "Core/StrandCompliantEuler.cpp"
  "Core/StrandCompliantManager.cpp"
  "Core/TimingUtilities.cpp"
  "Core/TwoDScene.cpp"
  "Core/WetHairCore.cpp"
  "Core/fluidsim.cpp"
  "Core/fluidsim2D.cpp"
  "Core/fluidsim3D.cpp"
  "Core/liangbarsky.cpp"
  "Core/sorter.cpp"
  "Core/viscosity3d.cpp"
  "Core/volume_fractions.cpp")
add_library(libWetHair::WetHairCore ALIAS WetHairCore)

target_link_libraries(WetHairCore
  PUBLIC
  Eigen3::Eigen
  TBB::tbb
  PRIVATE
  Threads::Threads)

target_include_directories(WetHairCore INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Core>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/libWetHair>)

add_executable (libWetHair
  "App/Camera.cpp"
  "App/ParticleSimulation.cpp"
  "App/RenderingUtilities.cpp"
  "App/StringUtilities.cpp"
  "App/TwoDSceneRenderer.cpp"
  "App/TwoDSceneSerializer.cpp"
  "App/TwoDSceneXMLParser.cpp"
  "App/TwoDimensionalDisplayController.cpp"
  "App/YImage.cpp"
  "App/main.cpp"
  "App/openglutils.cpp")

target_link_libraries(libWetHair PRIVATE
  libWetHair::WetHairCore
  ANTTWEAKBAR::AntTweakBar
  Eigen3::Eigen
  FreeGLUT::freeglut
  GLEW::glew
  OpenGL::GLU
  PNG::PNG
  RapidXML::RapidXML
  TCLAP::tclap)

install(TARGETS libWetHair)

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
install(DIRECTORY Core/
  DESTINATION include/libWetHair
  FILES_MATCHING
  PATTERN *.h)

