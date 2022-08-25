# libWetHair Executable

append_files (Headers "h" . App Core Core/DER Core/DER/Forces Core/DER/Dependencies Core/pcgsolver)
append_files (Sources "cpp" . App Core Core/DER Core/DER/Forces Core/DER/Dependencies Core/pcgsolver)
append_files (Sources "c" . App Core Core/DER Core/DER/Forces Core/DER/Dependencies Core/pcgsolver)

include_directories (Core)

#find_package (TCLAP REQUIRED)
#include_directories (${TCLAP_INCLUDE_DIR})

# Eigen3 library is required
find_package (Eigen3 REQUIRED CONFIG)

# Locate OpenGL
find_package (OpenGL REQUIRED)

find_package (GLEW REQUIRED)
find_package (TBB REQUIRED CONFIG)

# Locate GLUT
find_package (FreeGLUT REQUIRED CONFIG)

# RapidXML library is required
find_package (RapidXML REQUIRED)

# TCLAP library is required
find_package (TCLAP REQUIRED)

# Locate AntTweakBar
find_package (ANTTWEAKBAR REQUIRED)

find_package (PNG REQUIRED)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

add_executable (libWetHair ${Headers} ${Templates} ${Sources})
target_link_libraries(libWetHair PRIVATE ${LIBWETHAIR_LIBRARIES}
  ANTTWEAKBAR::AntTweakBar
  Eigen3::Eigen
  FreeGLUT::freeglut
  GLEW::glew
  OpenGL::GLU
  PNG::PNG
  RapidXML::RapidXML
  TBB::tbb
  TCLAP::tclap
  Threads::Threads)

install(TARGETS libWetHair)
