# Find dependencies
if (NOT TARGET libWetHair::WetHairCore)
  find_package (libWetHair REQUIRED CONFIG)
endif (NOT TARGET libWetHair::WetHairCore)
find_package (ANTTWEAKBAR REQUIRED)
find_package (Eigen3 REQUIRED CONFIG)
find_package (FreeGLUT REQUIRED CONFIG)
find_package (GLEW REQUIRED)
find_package (OpenGL REQUIRED)
find_package (PNG REQUIRED)
find_package (RapidXML REQUIRED)
find_package (TCLAP REQUIRED)

add_executable (libWetHair
  "Camera.cpp"
  "ParticleSimulation.cpp"
  "RenderingUtilities.cpp"
  "StringUtilities.cpp"
  "TwoDSceneRenderer.cpp"
  "TwoDSceneSerializer.cpp"
  "TwoDSceneXMLParser.cpp"
  "TwoDimensionalDisplayController.cpp"
  "YImage.cpp"
  "main.cpp"
  "openglutils.cpp")

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
