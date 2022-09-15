# Find dependencies
if (NOT TARGET libWetHair::WetHairCore)
  find_package (libWetHair REQUIRED CONFIG)
endif (NOT TARGET libWetHair::WetHairCore)
include (glfw)
include (glew)
include (anttweakbar)
include (tclap)
include (rapidxml)
include (stb)

find_package (OpenGL REQUIRED)

add_executable (libWetHair
  "Camera.h"
  "ExecutableSimulation.h"
  "ParticleSimulation.h"
  "RenderingUtilities.h"
  "StringUtilities.h"
  "TwoDSceneRenderer.h"
  "TwoDSceneSerializer.h"
  "TwoDSceneXMLParser.h"
  "TwoDimensionalDisplayController.h"
  "openglutils.h"
  "Camera.cpp"
  "ParticleSimulation.cpp"
  "RenderingUtilities.cpp"
  "StringUtilities.cpp"
  "TwoDSceneRenderer.cpp"
  "TwoDSceneSerializer.cpp"
  "TwoDSceneXMLParser.cpp"
  "TwoDimensionalDisplayController.cpp"
  "main.cpp"
  "openglutils.cpp")

target_link_libraries(libWetHair PRIVATE
  libWetHair::WetHairCore
  glfw::glfw
  anttweakbar::anttweakbar
  Eigen3::Eigen
  GLEW::glew
  OpenGL::GLU
  rapidxml::rapidxml
  tclap::tclap
  stb::stb)

install(TARGETS libWetHair)
