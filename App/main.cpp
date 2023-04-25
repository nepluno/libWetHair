//
// This file is part of the libWetHair open source project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright 2017 Yun (Raymond) Fei, Henrique Teles Maia, Christopher Batty,
// Changxi Zheng, and Eitan Grinspun
//

#include <AntTweakBar.h>
#include <sys/stat.h>
#include <tclap/CmdLine.h>

#include <Eigen/StdVector>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>

#ifdef WIN32
#include <direct.h>
#endif

// TODO: Clean these up, we don't need them all anymore.
#include <libWetHair/Force.h>
#include <libWetHair/MathDefs.h>
#include <libWetHair/TwoDScene.h>
#include <libWetHair/SceneStepper.h>
#include <libWetHair/SpringForce.h>
#include <libWetHair/TimingUtilities.h>

#include "DefaultXMLScene.h"
#include "ExecutableSimulation.h"
#include "ParticleSimulation.h"
#include "RenderingUtilities.h"
#include "StringUtilities.h"
#include "TwoDSceneRenderer.h"
#include "TwoDSceneSerializer.h"
#include "TwoDSceneXMLParser.h"
#include "TwoDimensionalDisplayController.h"

using namespace libwethair;

///////////////////////////////////////////////////////////////////////////////
// Contains the actual simulation, renderer, parser, and serializer
ExecutableSimulation* g_executable_simulation;

extern int getDCWindowWidth() {
  return g_executable_simulation->getWindowWidth();
}

extern int getDCWindowHeight() {
  return g_executable_simulation->getWindowHeight();
}

///////////////////////////////////////////////////////////////////////////////
// Rendering State
GLFWwindow* g_window = NULL;
bool g_rendering_enabled = true;
int g_dump_png = 0;
int g_dump_binary = 0;
int g_dump_readable = 0;
int g_res_scale = 1;
double g_sec_per_frame;
double g_last_time = timingutils::seconds();

bool g_should_dump_state = false;
bool g_should_dump_state_readable = false;
bool g_should_load_state = false;

int g_fluid_stream_level = 1;

renderingutils::Color g_bgcolor(0.93, 0.93, 0.93);

extern renderingutils::Color& getDCBackgroundColor() { return g_bgcolor; };
///////////////////////////////////////////////////////////////////////////////
// SVG Rendering State
bool g_svg_enabled = false;
std::string g_movie_dir;

///////////////////////////////////////////////////////////////////////////////
// Parser state
std::string g_xml_scene_file;
std::string g_description;
std::string g_scene_tag = "";

///////////////////////////////////////////////////////////////////////////////
// Scene input/output state
std::string g_binary_file_name;
std::string g_short_file_name;

bool g_simulate_initstate = false;
std::string g_initstate_file_name;
std::ifstream g_binary_input;

///////////////////////////////////////////////////////////////////////////////
// Simulation state
bool g_paused = true;
int g_num_steps = 0;
bool g_simulation_ran_to_completion = false;

int g_prev_step = 0;

///////////////////////////////////////////////////////////////////////////////
// Simulation functions

void dumpPNG(const std::string& filename);
void dumpstate(const std::string& filename);
void dumpstate_readable(const std::vector<std::string>& filename_fluids,
                        const std::string& filename_hairs,
                        const std::string& filename_flows,
                        const std::string& filename_bd_single,
                        const std::string& filename_bd_double,
                        const std::string& fn_pe, const std::string& fn_poe,
                        const std::string& fn_ppp);

void stepSystem() {
  // Step the system forward in time
  g_executable_simulation->stepSystem();

  int current_step = (int)g_executable_simulation->getCurrentFrame();
  scalar current_time = g_executable_simulation->getCurrentTime();

  // Determine if the simulation is complete
  if (!g_simulation_ran_to_completion && current_step >= g_num_steps) {
    std::cout << outputmod::startpink << "message: " << outputmod::endpink
              << "Simulation complete at time " << current_time << std::endl;
    g_simulation_ran_to_completion = true;
    g_paused = true;
    return;
  }

  auto& statistics = g_executable_simulation->getTimingStatistics();
  std::cout << "************************" << std::endl;
  scalar time_total = 0.0;
  for (size_t i = 0; i < statistics.size(); ++i) {
    time_total += statistics[i];
  }

  auto& stepper_stat = g_executable_simulation->getStepperTimingStatistics();
  scalar time_total_stepper = 0.0;
  for (size_t i = 0; i < stepper_stat.size(); ++i) {
    time_total_stepper += stepper_stat[i];
  }

  std::cout << "Stats: [Anim Time: " << current_time
            << " Compute Time: " << time_total
            << " FPS: " << (current_step / time_total) << "]" << std::endl;

  for (size_t i = 0; i < statistics.size(); ++i) {
    std::cout << i << ": " << (statistics[i] / (scalar)current_step) << "sec - "
              << (statistics[i] / time_total * 100.0) << "%" << std::endl;
  }

  std::cout << "Stepper Stats:" << std::endl;

  for (size_t i = 0; i < stepper_stat.size(); ++i) {
    std::cout << i << ": " << (stepper_stat[i] / (scalar)current_step)
              << "sec - " << (stepper_stat[i] / time_total_stepper * 100.0)
              << "%" << std::endl;
  }
  std::cout << "************************" << std::endl;

  // If the user wants to save output to a binary
  if (g_dump_binary && !(current_step % g_dump_binary)) {
    std::stringstream oss;
    oss << g_short_file_name << "/frame" << std::setw(5) << std::setfill('0')
        << (current_step / g_dump_binary) << ".bin";
    dumpstate(oss.str());
  }

  if (g_dump_readable && !(current_step % g_dump_readable)) {
    std::vector<std::stringstream> oss_fluids(g_fluid_stream_level);
    std::vector<std::string> str_fluids(g_fluid_stream_level);
    for (int i = 0; i < g_fluid_stream_level; ++i) {
      oss_fluids[i] << g_short_file_name << "/frame" << std::setw(5)
                    << std::setfill('0')
                    << (current_step / std::max(g_dump_readable, 1))
                    << "_fluids_" << i << ".txt";
      str_fluids[i] = oss_fluids[i].str();
    }

    std::stringstream oss_hairs;
    oss_hairs << g_short_file_name << "/frame" << std::setw(5)
              << std::setfill('0')
              << (current_step / std::max(g_dump_readable, 1)) << "_hairs.txt";

    std::stringstream oss_flows;
    oss_flows << g_short_file_name << "/frame" << std::setw(5)
              << std::setfill('0')
              << (current_step / std::max(g_dump_readable, 1)) << "_flows.txt";

    std::stringstream oss_bd_single;
    oss_bd_single << g_short_file_name << "/frame" << std::setw(5)
                  << std::setfill('0')
                  << (current_step / std::max(g_dump_readable, 1))
                  << "_single.txt";

    std::stringstream oss_bd_double;
    oss_bd_double << g_short_file_name << "/frame" << std::setw(5)
                  << std::setfill('0')
                  << (current_step / std::max(g_dump_readable, 1))
                  << "_double.txt";

    std::stringstream oss_pe;
    oss_pe << g_short_file_name << "/frame" << std::setw(5) << std::setfill('0')
           << (current_step / std::max(g_dump_readable, 1)) << "_pe.txt";

    std::stringstream oss_poe;
    oss_poe << g_short_file_name << "/frame" << std::setw(5)
            << std::setfill('0')
            << (current_step / std::max(g_dump_readable, 1)) << "_poe.txt";

    std::stringstream oss_ppp;
    oss_ppp << g_short_file_name << "/frame" << std::setw(5)
            << std::setfill('0')
            << (current_step / std::max(g_dump_readable, 1)) << "_ppp.txt";

    dumpstate_readable(str_fluids, oss_hairs.str(), oss_flows.str(),
                       oss_bd_single.str(), oss_bd_double.str(), oss_pe.str(),
                       oss_poe.str(), oss_ppp.str());

    for (int i = 0; i < g_fluid_stream_level; ++i) {
      oss_fluids[i].clear();
      str_fluids[i].clear();
    }

    oss_hairs.clear();
    oss_flows.clear();
    oss_bd_single.clear();
    oss_bd_double.clear();
    oss_pe.clear();
    oss_poe.clear();
    oss_ppp.clear();
  }

  // Update the state of the renderers
  if (g_rendering_enabled) g_executable_simulation->updateOpenGLRendererState();

  // If the user wants to generate a PNG movie
  if (g_dump_png && !(current_step % g_dump_png)) {
    std::stringstream oss;
    oss << g_short_file_name << "/frame" << std::setw(5) << std::setfill('0')
        << (current_step / g_dump_png) << ".png";
    dumpPNG(oss.str());
    oss.clear();
  }
}

void syncScene() {}

void headlessSimLoop() {
  while (true) {
    stepSystem();
    int current_step = (int)g_executable_simulation->getCurrentFrame();

    if (current_step > g_num_steps) break;
  }
}

void dumpPNGsubprog(char* image, char* fnstr, int w, int h) {
  unsigned* uimage = (unsigned*) image;
  for (int i = 0; i < h / 2; ++i) {
    for (int j = 0; j < w; ++j) {
      std::swap(uimage[i * w + j], uimage[(h - i - 1) * w + j]);
    }
  }
  stbi_write_png(fnstr, w, h, 4, image, w * 4);
  delete image;
  delete fnstr;
}

///////////////////////////////////////////////////////////////////////////////
// Rendering and UI functions

void dumpPNG(const std::string& filename) {
  char* fnstr = new char[filename.length() + 1];
  strcpy(fnstr, filename.data());
  int w, h;
  glfwGetFramebufferSize(g_window, &w, &h);

  char* img_data = new char[w * h * 4];

  glFinish();

  glPixelStorei(GL_PACK_ALIGNMENT, 4);
  glPixelStorei(GL_PACK_ROW_LENGTH, 0);
  glPixelStorei(GL_PACK_SKIP_ROWS, 0);
  glPixelStorei(GL_PACK_SKIP_PIXELS, 0);
  glReadBuffer(GL_BACK);

  glFinish();
  glReadPixels(0, 0, w, h, GL_RGBA,
               GL_UNSIGNED_BYTE, img_data);

  std::thread t(std::bind(dumpPNGsubprog, img_data, fnstr, w, h));

  t.detach();
}

void dumpStatesubprog(char* data, int length, char* fnstr) {
  std::ofstream ofs(fnstr);
  ofs.write(data, length);
  ofs.flush();
  ofs.close();

  delete[] data;
  delete[] fnstr;
}

void dumpstate(const std::string& filename) {
  char* fnstr = new char[filename.length() + 1];
  strcpy(fnstr, filename.data());
  // Attempt to open the binary
  std::ostringstream binary_output;
  if (binary_output.fail()) {
    std::cerr << outputmod::startred
              << "ERROR IN INITIALIZATION: " << outputmod::endred
              << "Failed to open stream. Exiting." << std::endl;
    exit(1);
  }
  // Save the initial conditions
  g_executable_simulation->serializeScene(binary_output);

  const std::string& str = binary_output.str();
  char* data = new char[str.size()];
  memcpy(data, str.c_str(), str.size());

  std::thread t(std::bind(dumpStatesubprog, data, (int)str.size(), fnstr));

  t.detach();
}

void dumpstate_readable(const std::vector<std::string>& filename_fluids,
                        const std::string& filename_hairs,
                        const std::string& filename_flows,
                        const std::string& filename_bd_single,
                        const std::string& filename_bd_double,
                        const std::string& filename_pe,
                        const std::string& filename_poe,
                        const std::string& filename_ppp) {
  std::vector<char*> fnstr_fluids(g_fluid_stream_level);

  for (int i = 0; i < g_fluid_stream_level; ++i) {
    fnstr_fluids[i] = new char[filename_fluids[i].length() + 1];
    strcpy(fnstr_fluids[i], filename_fluids[i].data());
  }

  char* fnstr_hairs = new char[filename_hairs.length() + 1];
  strcpy(fnstr_hairs, filename_hairs.data());

  char* fnstr_flows = new char[filename_flows.length() + 1];
  strcpy(fnstr_flows, filename_flows.data());

  char* fnstr_bd_single = new char[filename_bd_single.length() + 1];
  strcpy(fnstr_bd_single, filename_bd_single.data());

  char* fnstr_bd_double = new char[filename_bd_double.length() + 1];
  strcpy(fnstr_bd_double, filename_bd_double.data());

  char* fnstr_pe = new char[filename_pe.length() + 1];
  strcpy(fnstr_pe, filename_pe.data());

  char* fnstr_poe = new char[filename_poe.length() + 1];
  strcpy(fnstr_poe, filename_poe.data());

  char* fnstr_ppp = new char[filename_ppp.length() + 1];
  strcpy(fnstr_ppp, filename_ppp.data());

  std::vector<std::ostringstream> os_fluids(g_fluid_stream_level);
  std::ostringstream os_hairs;
  std::ostringstream os_flows;
  std::ostringstream os_bd_single;
  std::ostringstream os_bd_double;
  std::ostringstream os_pe;
  std::ostringstream os_poe;
  std::ostringstream os_ppp;

  if (os_hairs.fail() || os_flows.fail() || os_bd_single.fail() ||
      os_bd_double.fail() || os_pe.fail() || os_poe.fail() || os_ppp.fail()) {
    std::cerr << outputmod::startred
              << "ERROR IN INITIALIZATION: " << outputmod::endred
              << "Failed to open stream. Exiting." << std::endl;
    exit(1);
  }

  g_executable_simulation->serializeSceneReadable(os_fluids, os_hairs, os_flows,
                                                  os_bd_single, os_bd_double,
                                                  os_pe, os_poe, os_ppp);

  for (int i = 0; i < g_fluid_stream_level; ++i) os_fluids[i].flush();

  os_hairs.flush();
  os_flows.flush();
  os_bd_single.flush();
  os_bd_double.flush();
  os_pe.flush();
  os_poe.flush();
  os_ppp.flush();

  for (int i = 0; i < g_fluid_stream_level; ++i) {
    const std::string& str_fluids = os_fluids[i].str();
    char* data_fluids = new char[str_fluids.size()];
    memcpy(data_fluids, str_fluids.c_str(), str_fluids.size());

    std::thread(std::bind(dumpStatesubprog, data_fluids, (int)str_fluids.size(),
                          fnstr_fluids[i]))
        .detach();
  }

  const std::string& str_hairs = os_hairs.str();
  char* data_hairs = new char[str_hairs.size()];
  memcpy(data_hairs, str_hairs.c_str(), str_hairs.size());

  const std::string& str_flows = os_flows.str();
  char* data_flows = new char[str_flows.size()];
  memcpy(data_flows, str_flows.c_str(), str_flows.size());

  const std::string& str_bd_single = os_bd_single.str();
  char* data_bd_single = new char[str_bd_single.size()];
  memcpy(data_bd_single, str_bd_single.c_str(), str_bd_single.size());

  const std::string& str_bd_double = os_bd_double.str();
  char* data_bd_double = new char[str_bd_double.size()];
  memcpy(data_bd_double, str_bd_double.c_str(), str_bd_double.size());

  const std::string& str_pe = os_pe.str();
  char* data_pe = new char[str_pe.size()];
  memcpy(data_pe, str_pe.c_str(), str_pe.size());

  const std::string& str_poe = os_poe.str();
  char* data_poe = new char[str_poe.size()];
  memcpy(data_poe, str_poe.c_str(), str_poe.size());

  const std::string& str_ppp = os_ppp.str();
  char* data_ppp = new char[str_ppp.size()];
  memcpy(data_ppp, str_ppp.c_str(), str_ppp.size());

  std::thread(std::bind(dumpStatesubprog, data_hairs, (int)str_hairs.size(),
                        fnstr_hairs))
      .detach();
  std::thread(std::bind(dumpStatesubprog, data_flows, (int)str_flows.size(),
                        fnstr_flows))
      .detach();
  std::thread(std::bind(dumpStatesubprog, data_bd_single,
                        (int)str_bd_single.size(), fnstr_bd_single))
      .detach();
  std::thread(std::bind(dumpStatesubprog, data_bd_double,
                        (int)str_bd_double.size(), fnstr_bd_double))
      .detach();
  std::thread(
      std::bind(dumpStatesubprog, data_pe, (int)str_pe.size(), fnstr_pe))
      .detach();
  std::thread(
      std::bind(dumpStatesubprog, data_poe, (int)str_poe.size(), fnstr_poe))
      .detach();
  std::thread(
      std::bind(dumpStatesubprog, data_ppp, (int)str_ppp.size(), fnstr_ppp))
      .detach();

  for (int i = 0; i < g_fluid_stream_level; ++i) {
    os_fluids[i].clear();
  }

  os_hairs.clear();
  os_flows.clear();
  os_bd_single.clear();
  os_bd_double.clear();
  os_pe.clear();
  os_poe.clear();
  os_ppp.clear();
}

void reshape(GLFWwindow* window, int w, int h) {
  TwWindowSize(w, h);

  g_executable_simulation->reshape(w, h);

  assert(renderingutils::checkGLErrors());
}

void drawHUD() { TwDraw(); }

void display() {
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glMatrixMode(GL_MODELVIEW);

  assert(g_executable_simulation != NULL);
  g_executable_simulation->renderSceneOpenGL();

  drawHUD();

  assert(renderingutils::checkGLErrors());
}

void keyboard(GLFWwindow* window, int key, int scancode, int action, int mods) {
  if (TwEventCharGLFW(key, action) || TwEventKeyGLFW3(window, key, scancode, action, mods)) {
    return;
  }

  if (action == GLFW_PRESS) {
    g_executable_simulation->keyboard(key);

    if (key == GLFW_KEY_ESCAPE) {
      TwTerminate();
      glfwDestroyWindow(g_window);
      glfwTerminate();

      exit(0);
    } else if (key == GLFW_KEY_SPACE) {
      g_paused = !g_paused;
    } else if (key == GLFW_KEY_C) {
      g_executable_simulation->centerCamera();
    } else if (key == GLFW_KEY_S) {
      stepSystem();
    } else if (key == GLFW_KEY_B) {
      g_should_dump_state = true;
    } else if (key == GLFW_KEY_G) {
      int current_step = (int)g_executable_simulation->getCurrentFrame();
      int sdiv = g_dump_png ? g_dump_png : 1;
      std::stringstream oss;
      oss << g_short_file_name << "/frame" << std::setw(5) << std::setfill('0')
          << (current_step / sdiv) << ".png";
      dumpPNG(oss.str());
      oss.clear();
    }
  }

  assert(renderingutils::checkGLErrors());
}

void mouse(GLFWwindow *window, int button, int action, int mods) {
  if (TwEventMouseButtonGLFW3(window, button, action, mods)) {
    return;
  }

  g_executable_simulation->mouse(button, action, g_res_scale);

  assert(renderingutils::checkGLErrors());
}

void motion(GLFWwindow *window, double x, double y) {
  if (TwEventCursorPosGLFW3(window, x * g_res_scale, y * g_res_scale)) {
    return;
  }

  g_executable_simulation->motion((int) x * g_res_scale, (int) y * g_res_scale);

  assert(renderingutils::checkGLErrors());
}

void idle() {
  // std::cout << "g_last_time: " << g_last_time << std::endl;
  // Trigger the next timestep
  double current_time = timingutils::seconds();
  // std::cout << "current_time: " << current_time << std::endl;
  // std::cout << "g_sec_per_frame: " << g_sec_per_frame << std::endl;
  if (!g_paused && current_time - g_last_time >= g_sec_per_frame) {
    g_last_time = current_time;
    stepSystem();

  }

  int current_step = (int)g_executable_simulation->getCurrentFrame();

  if (g_should_dump_state) {
    std::stringstream oss;
    oss << g_short_file_name << "/frame" << std::setw(5) << std::setfill('0')
        << (current_step / std::max(g_dump_binary, 1)) << ".bin";
    dumpstate(oss.str());
    g_should_dump_state = false;
  }

  if (g_should_dump_state_readable) {
    std::vector<std::stringstream> oss_fluids;
    std::vector<std::string> str_fluids;
    for (int i = 0; i < g_fluid_stream_level; ++i) {
      oss_fluids[i] << g_short_file_name << "/frame" << std::setw(5)
                    << std::setfill('0')
                    << (current_step / std::max(g_dump_readable, 1))
                    << "_fluids_" << i << ".txt";
      str_fluids[i] = oss_fluids[i].str();
    }

    std::stringstream oss_hairs;
    oss_hairs << g_short_file_name << "/frame" << std::setw(5)
              << std::setfill('0')
              << (current_step / std::max(g_dump_readable, 1)) << "_hairs.txt";

    std::stringstream oss_flows;
    oss_flows << g_short_file_name << "/frame" << std::setw(5)
              << std::setfill('0')
              << (current_step / std::max(g_dump_readable, 1)) << "_flows.txt";

    std::stringstream oss_bd_single;
    oss_bd_single << g_short_file_name << "/frame" << std::setw(5)
                  << std::setfill('0')
                  << (current_step / std::max(g_dump_readable, 1))
                  << "_single.txt";

    std::stringstream oss_bd_double;
    oss_bd_double << g_short_file_name << "/frame" << std::setw(5)
                  << std::setfill('0')
                  << (current_step / std::max(g_dump_readable, 1))
                  << "_double.txt";

    std::stringstream oss_pe;
    oss_bd_single << g_short_file_name << "/frame" << std::setw(5)
                  << std::setfill('0')
                  << (current_step / std::max(g_dump_readable, 1)) << "_pe.txt";

    std::stringstream oss_poe;
    oss_bd_double << g_short_file_name << "/frame" << std::setw(5)
                  << std::setfill('0')
                  << (current_step / std::max(g_dump_readable, 1))
                  << "_poe.txt";

    std::stringstream oss_ppp;
    oss_bd_double << g_short_file_name << "/frame" << std::setw(5)
                  << std::setfill('0')
                  << (current_step / std::max(g_dump_readable, 1))
                  << "_ppp.txt";

    dumpstate_readable(str_fluids, oss_hairs.str(), oss_flows.str(),
                       oss_bd_single.str(), oss_bd_double.str(), oss_pe.str(),
                       oss_poe.str(), oss_ppp.str());
    g_should_dump_state_readable = false;
  }

  assert(renderingutils::checkGLErrors());
}

void initializeOpenGLandGLFW(int argc, char** argv) {
  // Initialize GLFW
  if (!glfwInit())
    exit(-1);
  
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);
  
  g_window = glfwCreateWindow(1920, 1080, "HairSim", NULL, NULL);
  
  glfwSetKeyCallback(g_window, keyboard);
  glfwSetMouseButtonCallback(g_window, mouse);
  glfwSetCursorPosCallback(g_window, motion);
  glfwSetFramebufferSizeCallback(g_window, reshape);
  
  glfwMakeContextCurrent(g_window);
  glfwSwapInterval(1);

  GLenum err = glewInit();
  if (GLEW_OK != err) {
    std::cerr << "Error: " << glewGetErrorString(err) << std::endl;
    exit(-1);
  }

  glEnable(GL_MULTISAMPLE_ARB);
  glHint(GL_MULTISAMPLE_FILTER_HINT_NV, GL_NICEST);

  TwInit(TW_OPENGL, NULL);

  // Initialize OpenGL
  int w, h;
  glfwGetFramebufferSize(g_window, &w, &h);
  int ww, wh;
  glfwGetWindowSize(g_window, &ww, &wh);
  g_res_scale = w / ww;
  glViewport(0, 0, w, h);
  glClearColor(g_bgcolor.r, g_bgcolor.g, g_bgcolor.b, 1.0);

  glEnable(GL_BLEND);
  glEnable(GL_POLYGON_SMOOTH);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);

  GLubyte halftone[] = {
      0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55, 0xAA, 0xAA, 0xAA, 0xAA,
      0x55, 0x55, 0x55, 0x55, 0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55,
      0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55, 0xAA, 0xAA, 0xAA, 0xAA,
      0x55, 0x55, 0x55, 0x55, 0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55,
      0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55, 0xAA, 0xAA, 0xAA, 0xAA,
      0x55, 0x55, 0x55, 0x55, 0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55,
      0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55, 0xAA, 0xAA, 0xAA, 0xAA,
      0x55, 0x55, 0x55, 0x55, 0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55,
      0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55, 0xAA, 0xAA, 0xAA, 0xAA,
      0x55, 0x55, 0x55, 0x55, 0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55,
      0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55};

  glPolygonStipple(halftone);

  assert(renderingutils::checkGLErrors());
}

///////////////////////////////////////////////////////////////////////////////
// Parser functions

void loadScene(const std::string& file_name, const char* memory_str) {
  // Maximum time in the simulation to run for. This has nothing to do with run
  // time, cpu time, etc. This is time in the 'virtual world'.
  scalar max_time;
  // Maximum frequency, in wall clock time, to execute the simulation for. This
  // serves as a cap for simulations that run too fast to see a solution.
  scalar steps_per_sec_cap = 100.0;
  // Contains the center and 'scale factor' of the view
  renderingutils::Viewport view;

  Camera cam;
  // Load the simulation and pieces of rendring and UI state
  assert(g_executable_simulation == NULL);
  TwoDSceneXMLParser xml_scene_parser;

  bool cam_init = false;
  bool view_init = false;
  xml_scene_parser.loadExecutableSimulation(
      g_window,
      file_name, memory_str, g_initstate_file_name, false, g_rendering_enabled,
      &g_executable_simulation, view, cam, max_time, steps_per_sec_cap,
      g_bgcolor, g_description, g_scene_tag, cam_init, view_init);
  assert(g_executable_simulation != NULL);

  if (cam_init) {
    g_executable_simulation->setCamera(cam);
  } else if (view_init) {
    g_executable_simulation->setView(view);
  } else {
    // If the user did not request a custom viewport, try to compute a
    // reasonable default.
    g_executable_simulation->centerCamera(false);
  }

  // To cap the framerate, compute the minimum time a single timestep should
  // take
  g_sec_per_frame = 1.0 / steps_per_sec_cap;
  // Integer number of timesteps to take
  g_num_steps = (int)ceil(max_time / g_executable_simulation->getDt());
}

void parseCommandLine(int argc, char** argv) {
  try {
    TCLAP::CmdLine cmd(
        "A Multi-Scale Model for Simulating Liquid-Hair Interactions");

    // XML scene file to load
    TCLAP::ValueArg<std::string> scene("s", "scene",
                                       "Simulation to run; an xml scene file",
                                       false, "", "string", cmd);

    // Begin the scene paused or running
    TCLAP::ValueArg<bool> paused(
        "p", "paused", "Begin the simulation paused if 1, running if 0", false,
        true, "boolean", cmd);

    // Begin the scene paused or running
    TCLAP::ValueArg<int> dumppng("g", "generate", "Generate PNG if 1, not if 0",
                                 false, 0, "integer", cmd);

    // Run the simulation with rendering enabled or disabled
    TCLAP::ValueArg<bool> display(
        "d", "display",
        "Run the simulation with display enabled if 1, without if 0", false,
        true, "boolean", cmd);

    // Readable-file (used by Houdini) to save output to
    TCLAP::ValueArg<int> dumpreadable(
        "o", "readableoutput", "readable file to save simulation state to",
        false, 0, "integer", cmd);

    // File to load for init
    TCLAP::ValueArg<std::string> init(
        "j", "initfile",
        "Binary file to load simulation state from for initialization", false,
        "", "string", cmd);

    cmd.parse(argc, argv);

    if (!scene.isSet()) {
      cmd.getOutput()->usage(cmd);
    }

    g_xml_scene_file = scene.getValue();
    g_paused = paused.getValue();
    g_rendering_enabled = display.getValue();
    g_dump_png = dumppng.getValue();
    g_dump_readable = dumpreadable.getValue();

    if (init.isSet()) {
      g_simulate_initstate = true;
      g_initstate_file_name = init.getValue();
    }

  } catch (TCLAP::ArgException& e) {
    std::cerr << "error: " << e.what() << std::endl;
    exit(1);
  }
}

///////////////////////////////////////////////////////////////////////////////
// Various support functions

void cleanupAtExit() {
  if (g_executable_simulation != NULL) {
    delete g_executable_simulation;
    g_executable_simulation = NULL;
  }
}

std::ostream& libwethair_header(std::ostream& stream) {
  stream << outputmod::startgreen <<

      ".__  ._____.   __      __        __     ___ ___        .__        \n"
         << "|  | |__\\_ |__/  \\    /  \\ _____/  |_  /   |   \\_____  "
            "|__|______ \n"
         << "|  | |  || __ \\   \\/\\/   // __ \\   __\\/    ~    \\__  \\ |  "
            "\\_  __ \\\n"
         << "|  |_|  || \\_\\ \\        /\\  ___/|  |  \\    Y    // __ \\|  "
            "||  | \\/\n"
         << "|____/__||___  /\\__/\\  /  \\___  >__|   \\___|_  /(____  "
            "/__||__|   \n"
         << "             \\/      \\/       \\/             \\/      \\/      "
            "     \n"
         << std::endl;

  return stream;
}

int main(int argc, char** argv) {
  // srand('SIGG' ^ 'RAPH');
  srand('S' ^ 'R');

  std::ios::sync_with_stdio(false);

  // Parse command line arguments
  parseCommandLine(argc, argv);

  bool loadDefaultScene = false;

  if (g_xml_scene_file.empty()) {
    g_xml_scene_file = default_xml_scene_name;
    g_short_file_name = default_xml_scene_short_name;
    loadDefaultScene = true;
  } else {
    std::vector<std::string> pathes;

    stringutils::split(g_xml_scene_file, '/', pathes);

    std::vector<std::string> path_first;

    stringutils::split(pathes[pathes.size() - 1], '.', path_first);

    g_short_file_name = path_first[0];
  }

#ifdef WIN32
  _mkdir(g_short_file_name.c_str());
#else
  mkdir(g_short_file_name.c_str(), 0777);
#endif

  // Function to cleanup at progarm exit
  atexit(cleanupAtExit);

  // Initialization for OpenGL and GLFW
  if (g_rendering_enabled) initializeOpenGLandGLFW(argc, argv);
  
  // Load the user-specified scene
  if (loadDefaultScene) {
    loadScene(g_xml_scene_file, default_xml_str);
  } else {
    loadScene(g_xml_scene_file, NULL);
  }

  if (g_rendering_enabled) {
    int w, h;
    glfwGetFramebufferSize(g_window, &w, &h);
    reshape(g_window, w, h);
    g_executable_simulation->initializeOpenGLRenderer();
  }

  // Print a header
  std::cout << libwethair_header << std::endl;

  // Print some status info about this HairFlow build
#ifdef NDEBUG
  std::cout << outputmod::startblue << "Build type: " << outputmod::endblue
            << "Release" << std::endl;
#else
  std::cout << outputmod::startblue << "Build type: " << outputmod::endblue
            << "Debug" << std::endl;
#endif
#ifdef EIGEN_VECTORIZE
  std::cout << outputmod::startblue << "Vectorization: " << outputmod::endblue
            << "Enabled" << std::endl;
#else
  std::cout << outputmod::startblue << "Vectorization: " << outputmod::endblue
            << "Disabled" << std::endl;
#endif

  std::cout << outputmod::startblue << "Scene: " << outputmod::endblue
            << g_xml_scene_file << std::endl;
  std::cout << outputmod::startblue << "Integrator: " << outputmod::endblue
            << g_executable_simulation->getSolverName() << std::endl;
  std::cout << outputmod::startblue << "Description: " << outputmod::endblue
            << g_description << std::endl;

  if (g_rendering_enabled) {
    while (!glfwWindowShouldClose(g_window)) {
      idle();
      display();
      
      glfwSwapBuffers(g_window);
      glfwPollEvents();
    }
    
    glfwDestroyWindow(g_window);
    glfwTerminate();
  } else {
    headlessSimLoop();
  }
  return EXIT_SUCCESS;
}
