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
if (OPENGL_FOUND)
  include_directories (${OPENGL_INCLUDE_DIR})
  set (LIBWETHAIR_LIBRARIES ${LIBWETHAIR_LIBRARIES} ${OPENGL_LIBRARIES})
else (OPENGL_FOUND)
  message (SEND_ERROR "Unable to locate OpenGL")
endif (OPENGL_FOUND)

find_package (GLEW REQUIRED)
find_package (TBB REQUIRED CONFIG)

# Locate GLUT
find_package (FreeGLUT REQUIRED CONFIG)

# RapidXML library is required
find_package (RapidXML REQUIRED)
if (RAPIDXML_FOUND)
  include_directories (${RAPIDXML_INCLUDE_DIR})
else (RAPIDXML_FOUND)
  message (SEND_ERROR "Unable to locate RapidXML")
endif (RAPIDXML_FOUND)

# TCLAP library is required
find_package (TCLAP REQUIRED)

# Locate AntTweakBar
find_package (ANTTWEAKBAR REQUIRED)
if (ANTTWEAKBAR_FOUND)
  include_directories( ${ANT_TWEAK_BAR_INCLUDE_DIR} )
  set (LIBWETHAIR_LIBRARIES ${LIBWETHAIR_LIBRARIES} ${ANT_TWEAK_BAR_LIBRARY})  
else (ANTTWEAKBAR_FOUND)
  message (SEND_ERROR "Unable to locate ANTTWEAKBAR")
endif (ANTTWEAKBAR_FOUND)  

find_package (PNG REQUIRED)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

add_executable (libWetHair ${Headers} ${Templates} ${Sources})
target_link_libraries(libWetHair PRIVATE ${LIBWETHAIR_LIBRARIES}
  Eigen3::Eigen
  FreeGLUT::freeglut
  GLEW::glew
  PNG::PNG
  TBB::tbb
  TCLAP::tclap
  Threads::Threads)

install(TARGETS libWetHair)
