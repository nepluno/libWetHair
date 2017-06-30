//
// This file is part of the libWetHair open source project
//
// The code is licensed solely for academic and non-commercial use under the
// terms of the Clear BSD License. The terms of the Clear BSD License are
// provided below. Other licenses may be obtained by contacting the faculty
// of the Columbia Computer Graphics Group or a Columbia University licensing officer.
//
// The Clear BSD License
//
// Copyright 2017 Yun (Raymond) Fei, Henrique Teles Maia, Christopher Batty,
// Changxi Zheng, and Eitan Grinspun
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted (subject to the limitations in the disclaimer
// below) provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice, this
//  list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright notice,
//  this list of conditions and the following disclaimer in the documentation
//  and/or other materials provided with the distribution.
//
// * Neither the name of the copyright holder nor the names of its contributors may be used
//  to endorse or promote products derived from this software without specific
//  prior written permission.
//
// NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE GRANTED BY THIS
// LICENSE. THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
// GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
// OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.

#include "TwoDimensionalDisplayController.h"
#include "TwoDSceneRenderer.h"

void* g_display_controller = NULL;

template<int DIM>
TwoDimensionalDisplayController<DIM>::TwoDimensionalDisplayController( int width, int height )
: m_window_width(width)
, m_window_height(height)
, m_scale_factor(1.0)
, m_center_x(0.0)
, m_center_y(0.0)
, m_left_drag(false)
, m_right_drag(false)
, m_last_x(0)
, m_last_y(0)
, m_render(NULL)
, m_modifiers(0)
, m_right_part_click(false)
{
  g_display_controller = (void*) this;
}

template<int DIM>
void TwoDimensionalDisplayController<DIM>::initDefaultCamera()
{
  if(m_cam_stack.size() == 0) {
    m_cam_stack.push(m_camera);
  }
}

template<int DIM>
void TwoDimensionalDisplayController<DIM>::setWindowWidth(int w)
{
  m_window_width = w;
}

template<int DIM>
void TwoDimensionalDisplayController<DIM>::setWindowHeight(int h)
{
  m_window_height = h;
}

template<int DIM>
void TwoDimensionalDisplayController<DIM>::reshape( int w, int h )
{
  assert( renderingutils::checkGLErrors() );
  // Record the new width and height
  m_window_width = w;
  m_window_height = h;

  applyProjection();
  
  assert( renderingutils::checkGLErrors() );
}

template<>
void TwoDimensionalDisplayController<3>::applyProjection() const
{
  // Reset the coordinate system before modifying
  if(g_rendering_enabled) {
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    
    double fov, zNear, zFar;
    m_camera.getPerspective( fov, zNear, zFar );
    gluPerspective( fov, getRatio(),  zNear, zFar );
    
    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity();
    Eigen::Vector3d eye, center, up;
    m_camera.getLookAt( eye, center, up );
    gluLookAt( eye.x(), eye.y(), eye.z(), center.x(), center.y(), center.z(), up.x(), up.y(), up.z() );
    
    glutPostRedisplay();
  }
}

template<>
void TwoDimensionalDisplayController<2>::applyProjection() const
{
  // Reset the coordinate system before modifying
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  double ratio = getRatio();
  applyOrtho(ratio);
  // Render the scene
  glutPostRedisplay();
}

template<int DIM>
void TwoDimensionalDisplayController<DIM>::applyOrtho( double ratio, double add_scale ) const
{
  // Set the coordinate system to achieve the desired zoom level, center
  gluOrtho2D(m_center_x-m_scale_factor*ratio*add_scale,m_center_x+m_scale_factor*ratio*add_scale,m_center_y-m_scale_factor*add_scale,m_center_y+m_scale_factor*add_scale);
}

template<int DIM>
void TwoDimensionalDisplayController<DIM>::applyOrthoZ( double ratio, double add_scale ) const
{
  // Set the coordinate system to achieve the desired zoom level, center
  gluOrtho2D(m_center_x-m_scale_factor*ratio*add_scale,m_center_x+m_scale_factor*ratio*add_scale,m_center_z-m_scale_factor*add_scale,m_center_z+m_scale_factor*add_scale);
}

template<int DIM>
double TwoDimensionalDisplayController<DIM>::getRatio() const
{
  return (double) m_window_width / (double) m_window_height;
}

template<int DIM>
void TwoDimensionalDisplayController<DIM>::setRender(TwoDSceneRenderer<DIM>* render)
{
  m_render = render;
}

template<int DIM>
int TwoDimensionalDisplayController<DIM>::getWorldWidth() const
{
  double ratio = getRatio();
  return 2*m_scale_factor/ratio;
}

template<int DIM>
int TwoDimensionalDisplayController<DIM>::getWorldHeight() const
{
  return 2*m_scale_factor;
}

template<int DIM>
void TwoDimensionalDisplayController<DIM>::keyboard( unsigned char key, int x, int y )
{ 
  if( key == '-' || key == '_' )
  {
    m_scale_factor += 0.1;
    reshape(m_window_width,m_window_height);
  }
  else if( key == '=' || key == '+' )
  {
    m_scale_factor = std::max(0.1,m_scale_factor-0.1);
    reshape(m_window_width,m_window_height);
  }
  else if( m_render && (key == 'P' || key == 'p') )
  {
    m_render->selectNextParticleVisMode();
    glutPostRedisplay();
  }
  else if( m_render && (key == 'K' || key == 'k') )
  {
    m_cam_stack.push(m_camera);
  }
  else if( m_render && (key == 'L' || key == 'l') )
  {
    if(m_cam_stack.size() > 0) {
      m_camera = m_cam_stack.top();
      if(m_cam_stack.size() > 1)
        m_cam_stack.pop();
      
      applyProjection();
    }
  }
  else if( m_render && (key == 'E' || key == 'e') )
  {
    m_render->selectNextEdgeVisMode();
    glutPostRedisplay();
  }
  else if( m_render && (key == 'M' || key == 'm') )
  {
    m_render->switchShowEdgeNormal();
    glutPostRedisplay();
  }
  else if( m_render && (key == 'N' || key == 'n') )
  {
    m_render->switchShowParticleNormal();
    glutPostRedisplay();
  }
  else if( key == 'I' || key == 'i' )
  {
    if(DIM == 2) {
      std::cout << "<viewport cx=\"" << m_center_x << "\" cy=\"" << m_center_y << "\" size=\"" << m_scale_factor << "\" cz=\"0\"/>" << std::endl;
    } else {
      std::cout << m_camera << std::endl;
    }
  }
}

template<int DIM>
void TwoDimensionalDisplayController<DIM>::special( int key, int x, int y )
{
  if( GLUT_KEY_UP == key ) 
  {
    m_center_y += 0.1;
    reshape(m_window_width,m_window_height);
  }
  else if( GLUT_KEY_DOWN == key ) 
  {
    m_center_y -= 0.1;
    reshape(m_window_width,m_window_height);
  }
  else if( GLUT_KEY_LEFT == key ) 
  {
    m_center_x -= 0.1;
    reshape(m_window_width,m_window_height);
  }
  else if( GLUT_KEY_RIGHT == key ) 
  {
    m_center_x += 0.1;
    reshape(m_window_width,m_window_height);
  }
}  

template<>
void TwoDimensionalDisplayController<2>::mouse( int button, int state, int x, int y )
{
  if( !m_right_drag && button == GLUT_LEFT_BUTTON && state == GLUT_DOWN )
  {
    m_left_drag = true;
    m_last_x = x;
    m_last_y = y;
    m_modifiers = glutGetModifiers();
    m_right_part_click = false;
  }
  if( button == GLUT_LEFT_BUTTON && state == GLUT_UP )
  {
    m_left_drag = false;
    m_modifiers = 0;
    m_right_part_click = false;
  }
  
  if( !m_left_drag && button == GLUT_RIGHT_BUTTON && state == GLUT_DOWN )
  {
    m_right_drag = true;
    m_last_x = x;
    m_last_y = y;
    m_modifiers = glutGetModifiers();
    m_right_part_click = false;
  }
  if( button == GLUT_RIGHT_BUTTON && state == GLUT_UP )
  {
    m_right_drag = false;
    m_modifiers = 0;
    m_right_part_click = false;
  }
}

template<>
void TwoDimensionalDisplayController<3>::mouse( int button, int state, int x, int y )
{
  if( !m_right_drag && button == GLUT_LEFT_BUTTON && state == GLUT_DOWN )
  {
    m_right_part_click = (x > m_window_width / 2);

    m_left_drag = true;
    int iPart = 1;
    int mx = x % (m_window_width / iPart);
    int my = y;
    
    double r = (m_window_width / iPart) < m_window_height ? (m_window_width / iPart) : m_window_height;
    
    double nx = (2.0 * mx - m_window_width / iPart) / r;
    double ny = (m_window_height - 2.0 * my) / r;
    
    m_last_x = nx;
    m_last_y = ny;
    
    m_modifiers = glutGetModifiers();
  }
  if( button == GLUT_LEFT_BUTTON && state == GLUT_UP )
  {
    m_right_part_click = false;
    m_left_drag = false;
    m_modifiers = 0;
  }
  
  if( !m_left_drag && button == GLUT_RIGHT_BUTTON && state == GLUT_DOWN )
  {
    m_right_part_click = (x > m_window_width / 2);

    int iPart = 1;
    int mx = x % (m_window_width / iPart);
    int my = y;
    
    double r = (m_window_width / iPart) < m_window_height ? (m_window_width / iPart) : m_window_height;
    
    double nx = (2.0 * mx - m_window_width / iPart) / r;
    double ny = (m_window_height - 2.0 * my) / r;
    
    m_last_x = nx;
    m_last_y = ny;
    m_right_drag = true;
    m_modifiers = glutGetModifiers();

  }
  if( button == GLUT_RIGHT_BUTTON && state == GLUT_UP )
  {
    m_right_part_click = false;
    m_right_drag = false;
    m_modifiers = 0;
  }
}

template<int DIM>
void TwoDimensionalDisplayController<DIM>::translateViewZ( double dx, double dz )
{
  double percent_x = dx/((double)m_window_width);
  double percent_z = dz/((double)m_window_height);
  double translate_x = percent_x*2.0*m_scale_factor*((double)m_window_width)/((double)m_window_height);
  double translate_z = percent_z*2.0*m_scale_factor;
  m_center_x -= translate_x;
  m_center_z += translate_z;
  reshape(m_window_width,m_window_height);
}

template<int DIM>
void TwoDimensionalDisplayController<DIM>::translateView( double dx, double dy )
{
  double percent_x = dx/((double)m_window_width);
  double percent_y = dy/((double)m_window_height);
  double translate_x = percent_x*2.0*m_scale_factor*((double)m_window_width)/((double)m_window_height);
  double translate_y = percent_y*2.0*m_scale_factor;
  m_center_x -= translate_x;
  m_center_y += translate_y;
  reshape(m_window_width,m_window_height);
}

template<int DIM>
void TwoDimensionalDisplayController<DIM>::zoomView( double dx, double dy )
{
  double percent_x = dx/((double)m_window_width);
  double percent_y = dy/((double)m_window_height);
  
  double scale;
  if( std::fabs(percent_x) > std::fabs(percent_y) ) scale = -percent_x;
  else scale = percent_y;
  
  m_scale_factor += 2.0*scale;
  
  reshape(m_window_width,m_window_height);
}

template<int DIM>
void TwoDimensionalDisplayController<DIM>::getEye(Vector3s& eye) const
{
  m_camera.getEye(eye);
}

template<>
void TwoDimensionalDisplayController<2>::motion( int x, int y )
{
  if( m_left_drag ) 
  {
    double dx = x - m_last_x;
    double dy = y - m_last_y;
    m_last_x = x;
    m_last_y = y;
    translateView( dx, dy );
  }
  if( m_right_drag ) 
  {
    double dx = x - m_last_x;
    double dy = y - m_last_y;
    m_last_x = x;
    m_last_y = y;
    zoomView( dx, dy );
  }
}

template<>
void TwoDimensionalDisplayController<3>::motion( int x, int y )
{
  if( m_left_drag )
  {

    double ox = m_last_x;
    double oy = m_last_y;
    
    int iPart = 1;
    
    int mx = x % (m_window_width / iPart);
    int my = y;
    
    double r = (m_window_width / iPart) < m_window_height ? (m_window_width / iPart) : m_window_height;
    
    double nx = (2.0 * mx - m_window_width / iPart) / r;
    double ny = (m_window_height - 2.0 * my) / r;
    
    m_last_x = nx;
    m_last_y = ny;
    if(m_modifiers & GLUT_ACTIVE_SHIFT) {
      m_camera.pan(ox, oy, nx, ny);
    } else {
      m_camera.rotate(ox, oy, nx, ny);
    }
    applyProjection();


  }
  if( m_right_drag )
  {
    double ox = m_last_x;
    double oy = m_last_y;
    
    int iPart = 1;
    
    int mx = x % (m_window_width / iPart);
    int my = y;
    
    double r = (m_window_width / iPart) < m_window_height ? (m_window_width / iPart) : m_window_height;
    
    double nx = (2.0 * mx - m_window_width / iPart) / r;
    double ny = (m_window_height - 2.0 * my) / r;
    
    m_last_x = nx;
    m_last_y = ny;
    
    m_camera.zoom(ox, oy, nx, ny);
    
    applyProjection();
  }
}

template<int DIM>
int TwoDimensionalDisplayController<DIM>::getWindowWidth() const
{
  return m_window_width;
}

template<int DIM>
int TwoDimensionalDisplayController<DIM>::getWindowHeight() const
{
  return m_window_height;
}

template<int DIM>
double TwoDimensionalDisplayController<DIM>::getCenterX() const
{
  return m_center_x;
}

template<int DIM>
double TwoDimensionalDisplayController<DIM>::getCenterY() const
{
  return m_center_y;
}


template<int DIM>
void TwoDimensionalDisplayController<DIM>::setCenterX( double x )
{
  m_center_x = x;
}

template<int DIM>
void TwoDimensionalDisplayController<DIM>::setCenterY( double y )
{
  m_center_y = y;
}

template<int DIM>
void TwoDimensionalDisplayController<DIM>::setCenterZ( double z )
{
  m_center_z = z;
}

template<int DIM>
void TwoDimensionalDisplayController<DIM>::setScaleFactor( double scale )
{
  m_scale_factor = scale;
}

template<int DIM>
double TwoDimensionalDisplayController<DIM>::getScaleFactor() const
{
  return m_scale_factor;
}

template<int DIM>
void TwoDimensionalDisplayController<DIM>::initCamera(const renderingutils::Viewport& view)
{
  Vector3s bmin(view.cx - view.rx, view.cy - view.ry, view.cz - view.rz);
  Vector3s bmax(view.cx + view.rx, view.cy + view.ry, view.cz + view.rz);
  m_camera.init(bmin, bmax);
}

template<int DIM>
const Camera& TwoDimensionalDisplayController<DIM>::getCamera() const
{
  return m_camera;
}

template<int DIM>
Camera& TwoDimensionalDisplayController<DIM>::getCamera()
{
  return m_camera;
}

// explicit instantiations at bottom
template class TwoDimensionalDisplayController<2>;
template class TwoDimensionalDisplayController<3>;

template<> TwoDimensionalDisplayController<2>* GetDisplayController()
{
  return (TwoDimensionalDisplayController<2>*) g_display_controller;
}

template<> TwoDimensionalDisplayController<3>* GetDisplayController()
{
  return (TwoDimensionalDisplayController<3>*) g_display_controller;
}

