# libWetHair Executable

append_files (Headers "h" . App Core Core/DER Core/DER/Forces Core/DER/Dependencies Core/pcgsolver)
append_files (Sources "cpp" . App Core Core/DER Core/DER/Forces Core/DER/Dependencies Core/pcgsolver)
append_files (Sources "c" . App Core Core/DER Core/DER/Forces Core/DER/Dependencies Core/pcgsolver)

include_directories (Core)

#find_package (TCLAP REQUIRED)
#include_directories (${TCLAP_INCLUDE_DIR})

# Locate OpenGL
find_package (OpenGL REQUIRED)
if (OPENGL_FOUND)
  include_directories (${OPENGL_INCLUDE_DIR})
  set (LIBWETHAIR_LIBRARIES ${LIBWETHAIR_LIBRARIES} ${OPENGL_LIBRARIES})
else (OPENGL_FOUND)
  message (SEND_ERROR "Unable to locate OpenGL")
endif (OPENGL_FOUND)

find_package(GLEW REQUIRED)
if (GLEW_FOUND)
  include_directories(${GLEW_INCLUDE_DIRS})
  set (LIBWETHAIR_LIBRARIES ${LIBWETHAIR_LIBRARIES} ${GLEW_LIBRARIES})
else (GLEW_FOUND)
  message (SEND_ERROR "Unable to locate GLEW")
endif()

find_package (TBB REQUIRED)
if (TBB_FOUND)
  include_directories (${TBB_INCLUDE_DIRS})
  set (LIBWETHAIR_LIBRARIES ${LIBWETHAIR_LIBRARIES} ${TBB_LIBRARIES})
else (TBB_FOUND)
  message (SEND_ERROR "Unable to locate TBB")
endif (TBB_FOUND)

# Locate GLUT
find_package (FreeGLUT REQUIRED)
if (FREEGLUT_FOUND)
  include_directories (${FREEGLUT_INCLUDE_DIRS})
  set (LIBWETHAIR_LIBRARIES ${LIBWETHAIR_LIBRARIES} ${FREEGLUT_LIBRARIES})
else (FREEGLUT_FOUND)
  message (SEND_ERROR "Unable to locate FREEGLUT")
endif (FREEGLUT_FOUND)

# RapidXML library is required
find_package (RapidXML REQUIRED)
if (RAPIDXML_FOUND)
  include_directories (${RAPIDXML_INCLUDE_DIR})
else (RAPIDXML_FOUND)
  message (SEND_ERROR "Unable to locate RapidXML")
endif (RAPIDXML_FOUND)

# TCLAP library is required
find_package (TCLAP REQUIRED)
if (TCLAP_FOUND)
  include_directories (${TCLAP_INCLUDE_PATH})
else (TCLAP_FOUND)
  message (SEND_ERROR "Unable to locate TCLAP")
endif (TCLAP_FOUND)

# Locate AntTweakBar
find_package (ANTTWEAKBAR REQUIRED)
if (ANTTWEAKBAR_FOUND)
  include_directories( ${ANT_TWEAK_BAR_INCLUDE_DIR} )
  set (LIBWETHAIR_LIBRARIES ${LIBWETHAIR_LIBRARIES} ${ANT_TWEAK_BAR_LIBRARY})  
else (ANTTWEAKBAR_FOUND)
  message (SEND_ERROR "Unable to locate ANTTWEAKBAR")
endif (ANTTWEAKBAR_FOUND)  

find_package (PNG REQUIRED)
if (PNG_FOUND)
  add_definitions (${PNG_DEFINITIONS})
  include_directories (${PNG_INCLUDE_DIR})
  set (LIBWETHAIR_LIBRARIES ${LIBWETHAIR_LIBRARIES} ${PNG_LIBRARIES})
endif (PNG_FOUND)

find_package(Threads REQUIRED)
if(THREADS_HAVE_PTHREAD_ARG)
  set_property(TARGET libWetHair PROPERTY COMPILE_OPTIONS "-pthread")
  set_property(TARGET libWetHair PROPERTY INTERFACE_COMPILE_OPTIONS "-pthread")
endif()
if(CMAKE_THREAD_LIBS_INIT)
  set(LIBWETHAIR_LIBRARIES ${LIBWETHAIR_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
endif()

add_executable (libWetHair ${Headers} ${Templates} ${Sources})
target_link_libraries (libWetHair ${LIBWETHAIR_LIBRARIES} "-lpthread" )

INSTALL_TARGETS(/bin libWetHair)
