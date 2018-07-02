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

#include "ParticleSimulation.h"
#include "TwoDSceneRenderer.h"
#include "fluidsim2D.h"
#include "fluidsim3D.h"
#include "TimingUtilities.h"

#include <iomanip>
#include <AntTweakBar.h>

//#define WRITE_PC_VARS

template<int DIM>
ParticleSimulation<DIM>::ParticleSimulation( TwoDScene<DIM>* scene, SceneStepper<DIM>* scene_stepper, TwoDSceneRenderer<DIM>* scene_renderer, const std::vector<Script>& scripts )
: m_scene(scene)
, m_scene_stepper(scene_stepper)
, m_scene_renderer(scene_renderer)
, m_scripts(scripts)
, m_display_controller(1920, 1080)
{
  m_wet_hair_core = new WetHairCore<DIM>( scene, scene_stepper, [&] (double dt) {
    return stepScript(dt);
  });
}

template<int DIM>
ParticleSimulation<DIM>::~ParticleSimulation()
{
  if( m_wet_hair_core != NULL )
  {
    delete m_wet_hair_core;
    m_wet_hair_core = NULL;
  }
  
  if( m_scene_renderer != NULL )
  {
    delete m_scene_renderer;
    m_scene_renderer = NULL;
  }
}

template<int DIM>
const std::vector<scalar>& ParticleSimulation<DIM>::getStepperTimingStatistics() const
{
  return m_wet_hair_core->getStepperTimingStatistics();
}
/////////////////////////////////////////////////////////////////////////////
// Simulation Control Functions

// #define USE_CFL

template<int DIM>
void ParticleSimulation<DIM>::stepSystem( )
{
  m_wet_hair_core->stepSystem();
  
  if( m_scene_renderer ) m_scene_renderer->updateOpenGLRenderer( *(m_wet_hair_core->getScene()), g_rendering_enabled );
}

/////////////////////////////////////////////////////////////////////////////
// Rendering Functions

template<int DIM>
void ParticleSimulation<DIM>::initializeOpenGLRenderer()
{
  m_scene_renderer->initializeOpenGLRenderer( *(m_wet_hair_core->getScene()) );
  m_scene_renderer->updateOpenGLRenderer( *(m_wet_hair_core->getScene()), true );
  
  TwBar* bar = TwNewBar("Parameters");
  TwDefine(" GLOBAL help='A Multi-Scale Model for Simulating Liquid-Hair Interactions' ");
  TwDefine(" Parameters size='400 800' color='166 75 0' ");
  
  WetHairParameter& parameters = m_scene->getParameter();
  
  TwAddVarRO(bar, "t", TW_TYPE_DOUBLE, &m_wet_hair_core->getCurrentTime(), " help='Current Time (s)' ");
  
  TwAddVarRW(bar, "dt", TW_TYPE_DOUBLE, &parameters.dt, " min=0.0001 max=0.01 step=0.0001 help='Time step (s)' ");
  TwAddVarRW(bar, "liquid density", TW_TYPE_DOUBLE, &parameters.rho, " min=0.1 max=10.0 step=0.1 help='Liquid density  (g/cm^3)' ");
  TwAddVarRW(bar, "surface tension", TW_TYPE_DOUBLE, &parameters.sigma, " min=0.0 max=400.0 step=0.1 help='Surface tension coefficient (dyn/cm)' ");
  TwAddVarRW(bar, "viscosity", TW_TYPE_DOUBLE, &parameters.viscosity, " min=0.0 max=10000.0 step=0.0001 help='Liquid viscosity (dyn·s/cm^2)' ");
  TwAddVarRW(bar, "quadratic dragging", TW_TYPE_DOUBLE, &parameters.quadratic_dragging, " min=0.0 max=100.0 step=0.0001 help='Liquid quadratic dragging coefficient' ");
  TwAddVarRW(bar, "hair friction", TW_TYPE_DOUBLE, &parameters.friction, " min=0.0 max=1.0 step=0.0001 help='Hair friction coefficient' ");
  TwAddVarRW(bar, "max flow/radius ratio", TW_TYPE_DOUBLE, &parameters.max_limited_eta_prop, " min=0.0 max=30.0 step=0.1 help='Maximal ratio between height of reduced-flow and radius of hair' ");
  TwAddVarRW(bar, "min flow/radius ratio", TW_TYPE_DOUBLE, &parameters.regularizer_shell, " min=0.0 max=5.0 step=0.1 help='Minimal ratio between height of reduced-flow and radius of hair, used to alleviating singularity in stiffness matrix' ");
  TwAddVarRW(bar, "latitude", TW_TYPE_DOUBLE, &parameters.latitude, " min=-1.5707 max=1.5707 step=0.001 help='Geographical latitude (for Coriolis effect, radian)' ");
  TwAddVarRW(bar, "planet radius", TW_TYPE_DOUBLE, &parameters.earth_radius, " min=0.0 max=1e10 step=1e6 help='Planet radius (for Coriolis effect, cm)' ");
  TwAddVarRW(bar, "planet rotation", TW_TYPE_DOUBLE, &parameters.earth_rotation, " min=0.0 max=1e-3 step=1e-6 help='Angular velocity of planet rotation (for Coriolis effect, rad/s)' ");
  TwAddVarRW(bar, "reduced-liquid smoothness", TW_TYPE_DOUBLE, &parameters.height_smooth, " min=0.0 max=5.0 step=0.1 help='Coefficient on the diffusive term of reduced-liquid motion' ");
  TwAddVarRW(bar, "capillary sucking multiplier", TW_TYPE_DOUBLE, &parameters.capillary_accel_multiplier, " min=0.0 max=10.0 step=0.1 help='Multiplier on additional capillary sucking effect for hairs with different radius' ");
  TwAddVarRW(bar, "dripping height multiplier", TW_TYPE_DOUBLE, &parameters.dripping_radius_multiplier, " min=0.0 max=1.0 step=0.01 help='Multiplier on the maximal height of reduced liquid held by hairs' ");
  TwAddVarRW(bar, "bulk threshold multiplier", TW_TYPE_DOUBLE, &parameters.bulk_threshold_multiplier, " min=0.0 max=1.0 step=0.01 help='Multiplier on the threshold to differentiate bulk liquid and reduced one' ");
  TwAddVarRW(bar, "capturing rate", TW_TYPE_DOUBLE, &parameters.absorptionRate, " min=0.0 max=0.1 step=0.001 help='Rate for capturing liquid from hairs (s)' ");
  TwAddVarRW(bar, "air viscosity", TW_TYPE_DOUBLE, &parameters.airviscosity, " min=0.0 max=0.01 step=0.0001 help='Air viscosity (dyn·s/cm^2)' ");
  TwAddVarRW(bar, "porosity multiplier", TW_TYPE_DOUBLE, &parameters.drag_radius_multiplier, " min=0.0 max=30.0 step=0.01 help='Multiplier on the porosity term of drag force applied onto liquid (clamped into [0, 1]).' ");
  TwAddVarRW(bar, "hair-hair cohesion multiplier", TW_TYPE_DOUBLE, &parameters.hair_hair_cohesion_multiplier, " min=0.0 max=10.0 step=0.001 help='Multiplier on the hair-hair cohesion effect' ");
  TwAddVarRW(bar, "hair-solid cohesion multiplier", TW_TYPE_DOUBLE, &parameters.hair_solid_cohesion_multiplier, " min=0.0 max=10.0 step=0.001 help='Multiplier on the hair-solid cohesion effect' ");
  TwAddVarRW(bar, "hair-hair radius multiplier", TW_TYPE_DOUBLE, &parameters.radius_multiplier, " min=0.0 max=10.0 step=0.01 help='Multiplier for hair radius used in hair-hair collision' ");
  TwAddVarRW(bar, "hair-solid radius multiplier", TW_TYPE_DOUBLE, &parameters.radius_multiplier_planar, " min=0.0 max=10.0 step=0.01 help='Multiplier for hair radius used in hair-solid collision' ");
  TwAddVarRW(bar, "hair-hair collision", TW_TYPE_DOUBLE, &parameters.collision_stiffness, " min=0.0 max=50000.0 step=100 help='Stiffness for hair-hair collision' ");
  TwAddVarRW(bar, "hair-solid collision", TW_TYPE_DOUBLE, &parameters.collision_stiffness_planar, " min=0.0 max=50000.0 step=100 help='Stiffness for hair-solid collision' ");
  TwAddVarRW(bar, "hair-hair damping", TW_TYPE_DOUBLE, &parameters.damping_multiplier, " min=0.0 max=1.0 step=0.01 help='Damping for hair-hair collision' ");
  TwAddVarRW(bar, "hair-solid damping", TW_TYPE_DOUBLE, &parameters.damping_multiplier_planar, " min=0.0 max=1.0 step=0.01 help='Damping for hair-solid collision' ");
  TwAddVarRW(bar, "hair-solid friction", TW_TYPE_DOUBLE, &parameters.friction_multiplier_planar, " min=0.0 max=1.0 step=0.01 help='Friction for hair-solid collision' ");
  
  TwAddVarRW(bar, "gravity x", TW_TYPE_DOUBLE, &(parameters.gravity.data()[0]), " min=-10000.0 max=10000.0 step=1.0 help='Gravitational Acceleration in X-Axis' ");
  
  TwAddVarRW(bar, "gravity y", TW_TYPE_DOUBLE, &(parameters.gravity.data()[1]), " min=-10000.0 max=10000.0 step=1.0 help='Gravitational Acceleration in Y-Axis' ");
  
  TwAddVarRW(bar, "gravity z", TW_TYPE_DOUBLE, &(parameters.gravity.data()[2]), " min=-10000.0 max=10000.0 step=1.0 help='Gravitational Acceleration in Z-Axis' ");
  
  TwAddVarRW(bar, "# hair substeps", TW_TYPE_UINT32, &parameters.hairsteps, " min=1 max=100 step=1 help='Number of Substeps for Hairs' ");
  TwAddVarRW(bar, "# reduced-liquid substeps", TW_TYPE_UINT32, &parameters.swesteps, " min=1 max=100 step=1 help='Number of Substeps for Reduced Liquid' ");
  TwAddVarRW(bar, "# fluid-correction steps", TW_TYPE_UINT32, &parameters.fluidcorrectionsteps, " min=0 max=100 step=1 help='Number of Steps to Perform Correction to Fluid Particles' ");
  
  TwAddVarRW(bar, "NO fluids", TW_TYPE_BOOL8, &parameters.no_fluids, " help='Turn On to Ignore Fluids' ");
  TwAddVarRW(bar, "NO reduced-liquid", TW_TYPE_BOOL8, &parameters.no_swe, " help='Turn On to Ignore Reduced Liquid' ");
  TwAddVarRW(bar, "NO capturing", TW_TYPE_BOOL8, &parameters.no_absorb, " help='Turn On to Forbid Liquid Capturing' ");
  TwAddVarRW(bar, "NO dripping", TW_TYPE_BOOL8, &parameters.no_dripping, " help='Turn On to Forbid Liquid Dripping' ");
  TwAddVarRW(bar, "NO fictitious", TW_TYPE_BOOL8, &parameters.no_fictitious, " help='Turn On to Forbid Fictitious Force' ");
  
  TwAddVarRW(bar, "dripping @ root", TW_TYPE_BOOL8, &parameters.drippingnear, " help='Turn On to Drip from Roots' ");
  TwAddVarRW(bar, "dripping @ tip", TW_TYPE_BOOL8, &parameters.drippingfar, " help='Turn On to Drip from Tips' ");
  TwAddVarRW(bar, "dripping @ middle", TW_TYPE_BOOL8, &parameters.drippingmiddle, " help='Turn On to Drip in the Middle' ");
  TwAddVarRW(bar, "viscous solve", TW_TYPE_BOOL8, &parameters.viscous_solve, " help='Turn On to Solve Viscous Effect for Bulk Liquid' ");
  TwAddVarRW(bar, "CTCD", TW_TYPE_BOOL8, &parameters.use_ctcd, " help='Turn On for Continuous-Time Collision Detection on Hairs' ");
  TwAddVarRW(bar, "apply coriolis force", TW_TYPE_BOOL8, &parameters.apply_coriolis, " help='Turn On for Coriolis force on 3D Liquid' ");
  TwAddVarRW(bar, "reduce volume loss", TW_TYPE_BOOL8, &parameters.global_volume_control, " help='Turn On to Reduce Volume Loss for Bulk Liquid' ");
  
  TwAddVarRW(bar, "individual transfer", TW_TYPE_BOOL8, &parameters.individual_transfer, " help='Turn On to Forbid Liquid Transfer between Hairs' ");
  TwAddVarRW(bar, "volume summary", TW_TYPE_BOOL8, &parameters.volume_summary, " help='Turn On to Summarize for Volume Info' ");
  
  TwEnumVal massmodeEV[ MUM_COUNT ] = {
    {MUM_MOMENTUM, "Update Momentum"}, {MUM_MASS_ONLY, "Mass-Only"}, {MUM_DIRECT_DIV, "Velocity-Scaling"}, {MUM_NONE, "Constant Mass"}
  };
  
  TwType massmodeType = TwDefineEnum("MassModeType", massmodeEV, MUM_COUNT);
  
  TwAddVarRW(bar, "Reduced-Flow Update Mode", massmodeType, &parameters.mass_update_mode, " help='Modes with different variables updated during the advection of reduced-liquid. Update Momentum: update both liquid mass and momentum of hairs; Mass-Only: only update liquid mass; Velocity-Scaling: directly divide hair velocity with new mass; Constant Mass: never update liquid mass on hairs or momentum of hairs' ");
  
  TwAddVarRW(bar, "Render Grid", TW_TYPE_BOOL8, &m_scene_renderer->m_draw_grid, " help='Turn On to Render Grid' ");
  TwAddVarRW(bar, "Render Particles", TW_TYPE_BOOL8, &m_scene_renderer->m_draw_particles, " help='Turn On to Render Particles' ");
  TwAddVarRW(bar, "Render Fluid Velocities", TW_TYPE_BOOL8, &m_scene_renderer->m_draw_velocities, " help='Turn On to Render Fluid Velocities' ");
  TwAddVarRW(bar, "Render Solid Objects", TW_TYPE_BOOL8, &m_scene_renderer->m_draw_boundaries, " help='Turn On to Render Solid Objects' ");  
}

template<int DIM>
void ParticleSimulation<DIM>::renderSceneOpenGL()
{
  assert( m_scene_renderer != NULL );
  m_scene_renderer->renderParticleSimulation(*(m_wet_hair_core->getScene()));
}

template<int DIM>
void ParticleSimulation<DIM>::renderSceneDifferencesOpenGL()
{
}

template<int DIM>
void ParticleSimulation<DIM>::updateOpenGLRendererState()
{
}

template<>
void ParticleSimulation<2>::computeCameraCenter( renderingutils::Viewport& view )
{
  Vector2s bb_min, bb_max;
  m_wet_hair_core->getBoundingBox(bb_min, bb_max);
  
  // Set center of view to center of bounding box
  view.cx = 0.5*(bb_max(0)+bb_min(0));
  view.cy = 0.5*(bb_max(1)+bb_min(1));
  
  // Set the zoom such that all particles are in view
  view.rx = 0.5*(bb_max(0)-bb_min(0));
  if( view.rx == 0.0 ) view.rx = 1.0;
  view.ry = 0.5*(bb_max(1)-bb_min(1));
  if( view.ry == 0.0 ) view.ry = 1.0;
}

template<>
void ParticleSimulation<3>::computeCameraCenter( renderingutils::Viewport& view )
{
  Vector3s bb_min, bb_max;
  m_wet_hair_core->getBoundingBox(bb_min, bb_max);
  
  // Set center of view to center of bounding box
  view.cx = 0.5*(bb_max(0)+bb_min(0));
  view.cy = 0.5*(bb_max(1)+bb_min(1));
  view.cz = 0.5*(bb_max(2)+bb_min(2));
  
  // Set the zoom such that all particles are in view
  view.rx = 0.5*(bb_max(0)-bb_min(0));
  if( view.rx == 0.0 ) view.rx = 1.0;
  view.ry = 0.5*(bb_max(1)-bb_min(1));
  if( view.ry == 0.0 ) view.ry = 1.0;
  view.rz = 0.5*(bb_max(2)-bb_min(2));
  if( view.rz == 0.0 ) view.rz = 1.0;
}

/////////////////////////////////////////////////////////////////////////////
// Serialization Functions
template<int DIM>
const std::vector<scalar>& ParticleSimulation<DIM>::getTimingStatistics() const
{
  return m_wet_hair_core->getTimingStatistics();
}

template<int DIM>
void ParticleSimulation<DIM>::serializeScene( std::ostream& outputstream )
{
  assert( m_scene != NULL );
  m_scene_serializer.serializeScene( *m_scene, m_scene_stepper, outputstream );
}

template<int DIM>
void ParticleSimulation<DIM>::serializeSceneReadable( std::vector< std::ostringstream >& osfluid, std::ostream& oshair, std::ostream& osflow, std::ostream& os_boundary_single, std::ostream& os_boundary_double, std::ostream& os_pe, std::ostream& os_poe, std::ostream& os_ppp )
{
  m_scene_serializer.serializeFluidReadable( *m_scene, osfluid );
  m_scene_serializer.serializeHairReadable( *m_scene, oshair );
  m_scene_serializer.serializeShallowFlowReadable( m_scene_renderer, *m_scene, osflow );
  m_scene_serializer.serializeBoundariesReadable( m_scene_renderer, *m_scene, os_boundary_single, os_boundary_double );
#ifdef WRITE_PC_VARS
  m_scene_serializer.serializePolygonalCohesionReadable(m_scene_renderer, *m_scene, os_pe, os_poe, os_ppp);
#endif
}

template<int DIM>
bool ParticleSimulation<DIM>::deSerializeSceneReadable( const std::vector< std::string > &filename_fluids, const std::string &filename_hairs, const std::string &filename_flows, const std::string &filename_bd_single, const std::string &filename_bd_double )
{
  bool f = m_scene_serializer.deSerializeFluidReadable( *m_scene, filename_fluids );
  bool h = m_scene_serializer.deSerializeHairReadable( *m_scene, filename_hairs );
  bool r = m_scene_serializer.deSerializeShallowFlowReadable( m_scene_renderer, *m_scene, filename_flows );

  return f || h || r;
}

template<int DIM>
void ParticleSimulation<DIM>::loadSerializedScene( std::istream& inputstream )
{
  assert( m_scene != NULL );
  m_scene_serializer.loadScene( *m_scene, m_scene_stepper, inputstream );
}
/////////////////////////////////////////////////////////////////////////////
// Status Functions

template<int DIM>
std::string ParticleSimulation<DIM>::getSolverName()
{
  assert( m_scene_stepper != NULL );
  return m_scene_stepper->getName();
}

template<int DIM>
const scalar& ParticleSimulation<DIM>::getDt() const
{
  return m_scene->getDt();
}

template<int DIM>
const scalar& ParticleSimulation<DIM>::getCurrentTime() const
{
  return m_wet_hair_core->getCurrentTime();
}

template<int DIM>
scalar ParticleSimulation<DIM>::getCurrentFrame() const
{
  return m_wet_hair_core->getCurrentTime() / m_scene->getDt();
}

template<int DIM>
void ParticleSimulation<DIM>::centerCamera(bool b_reshape)
{
  renderingutils::Viewport view;
  
  computeCameraCenter(view);
  scalar ratio;
  
  ratio = ((scalar)m_display_controller.getWindowHeight())/((scalar)m_display_controller.getWindowWidth());
  
  view.size = 1.2*std::max(ratio*view.rx,view.ry);
  
  m_display_controller.setCenterX(view.cx);
  m_display_controller.setCenterY(view.cy);
  m_display_controller.setCenterZ(view.cz);
  m_display_controller.setScaleFactor(view.size);
  m_display_controller.initCamera(view);
  
  if(b_reshape) m_display_controller.reshape(m_display_controller.getWindowWidth(),m_display_controller.getWindowHeight());
  
  m_display_controller.initDefaultCamera();
}

template<int DIM>
void ParticleSimulation<DIM>::keyboard( unsigned char key, int x, int y )
{
  m_display_controller.keyboard(key,x,y);
}

template<int DIM>
void ParticleSimulation<DIM>::reshape( int w, int h )
{
  m_display_controller.reshape(w,h);
}

template<int DIM>
void ParticleSimulation<DIM>::special( int key, int x, int y )
{
  m_display_controller.special(key, x, y);
}

template<int DIM>
void ParticleSimulation<DIM>::mouse( int button, int state, int x, int y )
{
  m_display_controller.mouse(button,state,x,y);
}

template<int DIM>
void ParticleSimulation<DIM>::translateView( double dx, double dy )
{
  m_display_controller.translateView(dx,dy);
}

template<int DIM>
void ParticleSimulation<DIM>::zoomView( double dx, double dy )
{
  m_display_controller.zoomView(dx,dy);
}

template<int DIM>
void ParticleSimulation<DIM>::motion( int x, int y )
{
  m_display_controller.motion(x,y);
}

template<int DIM>
int ParticleSimulation<DIM>::getWindowWidth() const
{
  return m_display_controller.getWindowWidth();
}

template<int DIM>
int ParticleSimulation<DIM>::getWindowHeight() const
{
  return m_display_controller.getWindowHeight();
}

template<int DIM>
void ParticleSimulation<DIM>::setWindowWidth(int w)
{
  m_display_controller.setWindowWidth(w);
}

template<int DIM>
void ParticleSimulation<DIM>::setWindowHeight(int h)
{
  m_display_controller.setWindowHeight(h);
}

template<int DIM>
TwoDimensionalDisplayController<DIM>* ParticleSimulation<DIM>::getDC()
{
  return &m_display_controller;
}

template<int DIM>
void ParticleSimulation<DIM>::setCenterX( double x )
{
  m_display_controller.setCenterX(x);
}

template<int DIM>
void ParticleSimulation<DIM>::setCenterY( double y )
{
  m_display_controller.setCenterY(y);
}

template<int DIM>
void ParticleSimulation<DIM>::setCenterZ( double z )
{
  m_display_controller.setCenterZ(z);
}

template<int DIM>
void ParticleSimulation<DIM>::setScaleFactor( double scale )
{
  m_display_controller.setScaleFactor(scale);
}

template<int DIM>
void ParticleSimulation<DIM>::setCamera( const Camera& cam )
{
  m_display_controller.getCamera() = cam;
  
  m_display_controller.setCenterX(cam.center_(0));
  m_display_controller.setCenterY(cam.center_(1));
  m_display_controller.setCenterZ(cam.center_(2));
  
  scalar size = 1.2 * cam.radius_;
  m_display_controller.setScaleFactor(size);
  m_display_controller.initDefaultCamera();
}

template<int DIM>
void ParticleSimulation<DIM>::setView( const renderingutils::Viewport& view )
{
  m_display_controller.setCenterX(view.cx);
  m_display_controller.setCenterY(view.cy);
  m_display_controller.setCenterZ(view.cz);
  m_display_controller.setScaleFactor(view.size);
  m_display_controller.initCamera(view);
  m_display_controller.initDefaultCamera();
}

inline scalar cosine_ease_function(const scalar& t, const scalar& t0, const scalar& t1, const scalar& ta, const scalar& tb, const scalar& amp,  const scalar& freq)
{
  const scalar Ta = ta - t0;
  const scalar Tt = t - t0;
  const scalar Tb = t1 - tb;
  const scalar Te = t1 - t0;
  const scalar w = 2.0 * M_PI * freq;
  
  if(t < t0 || t > t1) return 0.0;
  else if(t < ta) {
    const scalar t2 = Ta*w;
    const scalar t3 = cos(t2);
    const scalar t4 = sin(t2);
    return 1.0/(Ta*Ta*Ta)*(Tt*Tt)*(amp*t3*2.0+amp*Ta*t4*w)*-3.0+amp*1.0/(Ta*Ta)*Tt*(t3*3.0+Ta*t4*w)*2.0;
  } else if(t > tb) {
    const scalar t2 = Tb-Te;
    const scalar t3 = t2*w;
    const scalar t4 = Te-Tt;
    const scalar t5 = 1.0/(Tb*Tb*Tb);
    const scalar t6 = cos(t3);
    const scalar t7 = sin(t3);
    return -amp*t5*(Te*2.0-Tt*2.0)*(Tb*t6*3.0-Te*t6*2.0+t6*Tt*2.0+(Tb*Tb)*t7*w+Tb*t7*w*Tt-Tb*Te*t7*w)+amp*(t4*t4)*t5*(t6*2.0+Tb*t7*w);
  } else {
    return -amp * w * sin(w * Tt);
  }
}

inline scalar weno_ease_function(const scalar& t, const scalar& dt, const scalar& t0, const scalar& t1, const scalar& base_dt, const scalar& cur_pos, const std::vector<scalar>& bases)
{
  if(t < t0 || t > t1) return 0.0;
  
  scalar spos = (t - t0) / base_dt;
  int ipos = (int) spos;
  scalar fpos = spos - (scalar) ipos;
  
  int nb = (int) bases.size();
  
  scalar extracted_bases[] = {
    bases[mathutils::clamp(ipos - 2, 0, nb - 1)],
    bases[mathutils::clamp(ipos - 1, 0, nb - 1)],
    bases[mathutils::clamp(ipos - 0, 0, nb - 1)],
    bases[mathutils::clamp(ipos + 1, 0, nb - 1)],
    bases[mathutils::clamp(ipos + 2, 0, nb - 1)],
    bases[mathutils::clamp(ipos + 3, 0, nb - 1)]
  };
  
  scalar target_pos = mathutils::lerp_weno(extracted_bases, fpos);
  return (target_pos - cur_pos) / dt;
}

inline scalar cubic_ease_function(const scalar& t, const scalar& t0, const scalar& t1, const scalar& ta, const scalar& tb, const scalar& L)
{
  scalar yh = (L * 2.0) / (t1 - t0 + tb - ta);
  if(t < t0 || t > t1) return 0.0;
  else {
    if(t < ta) return (yh * (t0 - t) * (t0 - t) * (t0 - 3.0 * ta + 2.0 * t)) / ((t0 - ta) * (t0 - ta) * (t0 - ta));
    else if (t > tb) return (yh * (t1 - t) * (t1 - t) * (t1 - 3.0 * tb + 2.0 * t)) / ((t1 - tb) * (t1 - tb) * (t1 - tb));
    else return yh;
  }
}

template<int DIM>
bool ParticleSimulation<DIM>::stepScript( const scalar& dt )
{
  const scalar& current_time = m_wet_hair_core->getCurrentTime();
  
  bool updateSDF = false;
  for(auto& script : m_scripts)
  {
    if(current_time >= script.start && current_time + dt <= script.end)
    {
      switch (script.target) {
        case Script::CURLRADIUS:
        {
          scalar dx = script.v(0);
          scalar vel = 0;
          if(script.func == Script::CUBIC) {
            vel = cubic_ease_function(current_time + dt, script.start, script.end, script.start + script.ease_start, script.end - script.ease_end, dx);
          } else {
            std::cout << "UNIMPLEMENTED SCRIPT TYPE [" << script.type << "]!" << std::endl;
          }
          StrandEquilibriumParameters* parameter = m_scene->getStrandEquilibriumParameters( script.index );
          if(!parameter || !parameter->m_valid) break;
          parameter->m_curl_radius += vel * dt;
          parameter->m_dirty = true;
        }
          break;
        case Script::CURLDENSITY:
        {
          scalar dx = script.v(0);
          scalar vel = 0;
          if(script.func == Script::CUBIC) {
            vel = cubic_ease_function(current_time + dt, script.start, script.end, script.start + script.ease_start, script.end - script.ease_end, dx);
          } else {
            std::cout << "UNIMPLEMENTED SCRIPT TYPE [" << script.type << "]!" << std::endl;
          }
          StrandEquilibriumParameters* parameter = m_scene->getStrandEquilibriumParameters( script.index );
          if(!parameter || !parameter->m_valid) break;
          parameter->m_curl_density += vel * dt;
          parameter->m_dirty = true;
        }
          break;
        case Script::CAMERA:
        {
          switch (script.type) {
            case Script::ROTATE:
            {
              static scalar last_pos = 0.0;
              scalar vel = 0;
              if(script.func == Script::COSINE) {
                vel = cosine_ease_function(current_time + dt, script.start, script.end, script.start + script.ease_start, script.end - script.ease_end, script.amplitude, script.frequency);
              } else if(script.func == Script::CUBIC) {
                vel = cubic_ease_function(current_time + dt, script.start, script.end, script.start + script.ease_start, script.end - script.ease_end, script.v(3));
              } else if(script.func == Script::WENO) {
                vel = weno_ease_function(current_time + dt, dt, script.start, script.end, script.base_dt, script.base_pos, script.base_vertices);
                script.base_pos += vel * dt;
              }
              scalar da = vel * dt;
              Vector3s axis(script.v(0),script.v(1),script.v(2));
              m_display_controller.getCamera().rotate(axis, da, script.transform_global );
              m_display_controller.applyProjection();
            }
              break;
              
            default:
              std::cout << "UNIMPLEMENTED SCRIPT TYPE [" << script.type << "]!" << std::endl;
              break;
          }
        }
          break;
        case Script::ROOT:
        {
          switch (script.type) {
            case Script::TRANSLATE:
            {
              Vector3s vec(script.v(0),script.v(1),script.v(2));
              scalar dx = vec.norm();
              
              int group = std::max(0, script.index);
              Vector3s& trans = m_scene->getScriptedTranslate( group );
              if(dx > 0.0) {
                scalar vel = 0;
                VectorXs nv = vec / dx;

                if(script.func == Script::COSINE) {
                  vel = cosine_ease_function(current_time + dt, script.start, script.end, script.start + script.ease_start, script.end - script.ease_end, script.amplitude, script.frequency);
                } else if(script.func == Script::CUBIC) {
                  vel = cubic_ease_function(current_time + dt, script.start, script.end, script.start + script.ease_start, script.end - script.ease_end, dx);
                } else if(script.func == Script::WENO) {
                  vel = weno_ease_function(current_time + dt, dt, script.start, script.end, script.base_dt, script.base_pos, script.base_vertices);
                  script.base_pos += vel * dt;
                }
                
                trans += nv * vel * dt;
              }
              
              break;
            }
            case Script::ROTATE:
            {
              scalar vel = 0;
              if(script.func == Script::COSINE) {
                vel = cosine_ease_function(current_time + dt, script.start, script.end, script.start + script.ease_start, script.end - script.ease_end, script.amplitude, script.frequency);
              } else if(script.func == Script::CUBIC) {
                vel = cubic_ease_function(current_time + dt, script.start, script.end, script.start + script.ease_start, script.end - script.ease_end, script.v(3));
              } else if(script.func == Script::WENO) {
                vel = weno_ease_function(current_time + dt, dt, script.start, script.end, script.base_dt, script.base_pos, script.base_vertices);
                script.base_pos += vel * dt;
              }
              scalar da = vel * dt;
              Vector3s axis(script.v(0),script.v(1),script.v(2));
              
              int group = std::max(0, script.index);
              Eigen::AngleAxis<scalar> rot(da, axis);
              Eigen::Quaternion<scalar> qrot(rot);
              
              Eigen::Quaternion<scalar>& prot = m_scene->getScriptedRotation( group );
              if(script.transform_global) {
                prot = qrot * prot;
              } else {
                prot *= qrot;
              }
              break;
            }
            case Script::ABSORB:
            {
              scalar vel = 0;
              if(script.func == Script::COSINE) {
                vel = cosine_ease_function(current_time + dt, script.start, script.end, script.start + script.ease_start, script.end - script.ease_end, script.amplitude, script.frequency);
              } else if(script.func == Script::CUBIC) {
                vel = cubic_ease_function(current_time + dt, script.start, script.end, script.start + script.ease_start, script.end - script.ease_end, script.v(3));
              } else if(script.func == Script::WENO) {
                vel = weno_ease_function(current_time + dt, dt, script.start, script.end, script.base_dt, script.base_pos, script.base_vertices);
                script.base_pos += vel * dt;
              }
              scalar dh = vel * dt;
              
              int group = std::max(0, script.index);
              m_scene->uniformlyIncreaseHairLiquid(dh, group);
              m_scene->updateVolumeRecord();
              break;
            }
            default:
              std::cout << "UNIMPLEMENTED SCRIPT TYPE [" << script.type << "]!" << std::endl;
              break;
          }
          break;
        }
        case Script::SOLID:
        {
          FluidSim* fluidsim = m_scene->getFluidSim();
          
          switch (script.type) {
            case Script::TRANSLATE:
            {
              if(DIM == 2) {
                FluidSim2D* fluid2d = (FluidSim2D*) fluidsim;
                const std::vector< FluidSim::Boundary<2>* >& boundaries = fluid2d->get_boundaries();
                Vector2s dx = Vector2s(script.v(0), script.v(1));
                scalar ldx = dx.norm();
                if(script.index >= 0) {
                  if(boundaries[script.index]->type == FluidSim::BT_UNION || boundaries[script.index]->type == FluidSim::BT_INTERSECT) break;
                  FluidSim::SolidBoundary<2>* boundary = (FluidSim::SolidBoundary<2>*) boundaries[script.index];
                  if(ldx > 0.0) {
                    Vector2s ndx = dx / ldx;
                    scalar vel = 0;
                    if(script.func == Script::COSINE) {
                      vel = cosine_ease_function(current_time + dt, script.start, script.end, script.start + script.ease_start, script.end - script.ease_end, script.amplitude, script.frequency);
                    } else if(script.func == Script::CUBIC) {
                      vel = cubic_ease_function(current_time + dt, script.start, script.end, script.start + script.ease_start, script.end - script.ease_end, ldx);
                    } else if(script.func == Script::WENO) {
                      vel = weno_ease_function(current_time + dt, dt, script.start, script.end, script.base_dt, script.base_pos, script.base_vertices);
                      script.base_pos += vel * dt;
                    }
                    boundary->future_center(0) += ndx(0) * vel * dt;
                    boundary->future_center(1) += ndx(1) * vel * dt;
                  }
                } else {
                  for(auto& b : boundaries) {
                    if(b->type == FluidSim::BT_UNION || b->type == FluidSim::BT_INTERSECT) continue;
                    FluidSim::SolidBoundary<2>* boundary = (FluidSim::SolidBoundary<2>*) b;

                    if(ldx > 0.0) {
                      Vector2s ndx = dx / ldx;
                      scalar vel = 0;
                      if(script.func == Script::COSINE) {
                        vel = cosine_ease_function(current_time + dt, script.start, script.end, script.start + script.ease_start, script.end - script.ease_end, script.amplitude, script.frequency);
                      } else if(script.func == Script::CUBIC) {
                        vel = cubic_ease_function(current_time + dt, script.start, script.end, script.start + script.ease_start, script.end - script.ease_end, ldx);
                      } else if(script.func == Script::WENO) {
                        vel = weno_ease_function(current_time + dt, dt, script.start, script.end, script.base_dt, script.base_pos, script.base_vertices);
                        script.base_pos += vel * dt;
                      }
                      boundary->future_center(0) += ndx(0) * vel * dt;
                      boundary->future_center(1) += ndx(1) * vel * dt;
                    }
                  }
                }
              } else {
                FluidSim3D* fluid3d = (FluidSim3D*) fluidsim;
                const std::vector< FluidSim::Boundary<3>* >& boundaries = fluid3d->get_boundaries();
                Vector3s dx = Vector3s(script.v(0), script.v(1), script.v(2));
                scalar ldx = dx.norm();
                if(script.index >= 0) {
                  if(boundaries[script.index]->type == FluidSim::BT_UNION || boundaries[script.index]->type == FluidSim::BT_INTERSECT) break;
                  FluidSim::SolidBoundary<3>* boundary = (FluidSim::SolidBoundary<3>*) boundaries[script.index];
                  if(ldx > 0.0) {
                    Vector3s ndx = dx / ldx;
                    scalar vel = 0;
                    if(script.func == Script::COSINE) {
                      vel = cosine_ease_function(current_time + dt, script.start, script.end, script.start + script.ease_start, script.end - script.ease_end, script.amplitude, script.frequency);
                    } else if(script.func == Script::CUBIC) {
                      vel = cubic_ease_function(current_time + dt, script.start, script.end, script.start + script.ease_start, script.end - script.ease_end, ldx);
                    } else if(script.func == Script::WENO) {
                      vel = weno_ease_function(current_time + dt, dt, script.start, script.end, script.base_dt, script.base_pos, script.base_vertices);
                      script.base_pos += vel * dt;
                    }
                    boundary->future_center(0) += ndx(0) * vel * dt;
                    boundary->future_center(1) += ndx(1) * vel * dt;
                    boundary->future_center(2) += ndx(2) * vel * dt;
                  }
                } else {
                  for(auto& b : boundaries) {
                    if(b->type == FluidSim::BT_UNION || b->type == FluidSim::BT_INTERSECT) continue;
                    FluidSim::SolidBoundary<3>* boundary = (FluidSim::SolidBoundary<3>*) b;
                    
                    if(ldx > 0.0) {
                      Vector3s ndx = dx / ldx;
                      scalar vel = 0;
                      if(script.func == Script::COSINE) {
                        vel = cosine_ease_function(current_time + dt, script.start, script.end, script.start + script.ease_start, script.end - script.ease_end, script.amplitude, script.frequency);
                      } else if(script.func == Script::CUBIC) {
                        vel = cubic_ease_function(current_time + dt, script.start, script.end, script.start + script.ease_start, script.end - script.ease_end, ldx);
                      } else if(script.func == Script::WENO) {
                        vel = weno_ease_function(current_time + dt, dt, script.start, script.end, script.base_dt, script.base_pos, script.base_vertices);
                        script.base_pos += vel * dt;
                      }
                      boundary->future_center(0) += ndx(0) * vel * dt;
                      boundary->future_center(1) += ndx(1) * vel * dt;
                      boundary->future_center(2) += ndx(2) * vel * dt;
                    }
                  }
                }
              }
              updateSDF = updateSDF | script.updateSDF;
              
              break;
            }
            case Script::ROTATE:
            {
              if(DIM == 2) {
                FluidSim2D* fluid2d = (FluidSim2D*) fluidsim;
                const std::vector< FluidSim::Boundary<2>* >& boundaries = fluid2d->get_boundaries();
                scalar vel = 0;
                if(script.func == Script::COSINE) {
                  vel = cosine_ease_function(current_time + dt, script.start, script.end, script.start + script.ease_start, script.end - script.ease_end, script.amplitude, script.frequency);
                } else if(script.func == Script::CUBIC) {
                  vel = cubic_ease_function(current_time + dt, script.start, script.end, script.start + script.ease_start, script.end - script.ease_end, script.v(3));
                } else if(script.func == Script::WENO) {
                  vel = weno_ease_function(current_time + dt, dt, script.start, script.end, script.base_dt, script.base_pos, script.base_vertices);
                  script.base_pos += vel * dt;
                }
                scalar da = vel * dt;
                Vector3s axis(0.0,0.0,1.0);
                
                Eigen::AngleAxis<scalar> rot(da, axis);
                Eigen::Quaternion<scalar> qrot(rot);
                
                if(script.index >= 0) {
                  if(boundaries[script.index]->type == FluidSim::BT_UNION || boundaries[script.index]->type == FluidSim::BT_INTERSECT) break;
                  FluidSim::SolidBoundary<2>* boundary = (FluidSim::SolidBoundary<2>*) boundaries[script.index];
                  
                  Eigen::Quaternion<scalar>& prot = boundary->future_rot;
                  if(script.transform_global) {
                    prot = qrot * prot;
                  } else {
                    prot *= qrot;
                  }
                } else {
                  for(auto& b : boundaries)
                  {
                    if(b->type == FluidSim::BT_UNION || b->type == FluidSim::BT_INTERSECT) continue;
                    FluidSim::SolidBoundary<2>* boundary = (FluidSim::SolidBoundary<2>*) b;
                    
                    Eigen::Quaternion<scalar>& prot = boundary->future_rot;
                    if(script.transform_global) {
                      prot = qrot * prot;
                    } else {
                      prot *= qrot;
                    }
                  }
                }
              } else {
                FluidSim3D* fluid3d = (FluidSim3D*) fluidsim;
                const std::vector< FluidSim::Boundary<3>* >& boundaries = fluid3d->get_boundaries();
                
                scalar vel = 0;
                if(script.func == Script::COSINE) {
                  vel = cosine_ease_function(current_time + dt, script.start, script.end, script.start + script.ease_start, script.end - script.ease_end, script.amplitude, script.frequency);
                } else if(script.func == Script::CUBIC) {
                  vel = cubic_ease_function(current_time + dt, script.start, script.end, script.start + script.ease_start, script.end - script.ease_end, script.v(3));
                } else if(script.func == Script::WENO) {
                  vel = weno_ease_function(current_time + dt, dt, script.start, script.end, script.base_dt, script.base_pos, script.base_vertices);
                  script.base_pos += vel * dt;
                }
                scalar da = vel * dt;
                Vector3s axis(script.v(0),script.v(1),script.v(2));
                
                Eigen::AngleAxis<scalar> rot(da, axis);
                Eigen::Quaternion<scalar> qrot(rot);
                
                if(script.index >= 0) {
                  if(boundaries[script.index]->type == FluidSim::BT_UNION || boundaries[script.index]->type == FluidSim::BT_INTERSECT) break;
                  FluidSim::SolidBoundary<3>* boundary = (FluidSim::SolidBoundary<3>*) boundaries[script.index];
                  Eigen::Quaternion<scalar>& prot = boundary->future_rot;
                  if(script.transform_global) {
                    prot = qrot * prot;
                  } else {
                    prot *= qrot;
                  }
                } else {
                  for(auto& b : boundaries)
                  {
                    if(b->type == FluidSim::BT_UNION || b->type == FluidSim::BT_INTERSECT) continue;
                    FluidSim::SolidBoundary<3>* boundary = (FluidSim::SolidBoundary<3>*) b;
                    
                    Eigen::Quaternion<scalar>& prot = boundary->future_rot;
                    if(script.transform_global) {
                      prot = qrot * prot;
                    } else {
                      prot *= qrot;
                    }
                    
                  }
                }
              }
              
              updateSDF = updateSDF | script.updateSDF;
              break;
            }
            case Script::SCALE:
            {
              break;
            }
            default:
              std::cout << "UNIMPLEMENTED SCRIPT TYPE [" << script.type << "]!" << std::endl;
              break;
          }
          break;
        }
        case Script::SOURCE:
        {
          FluidSim* fluidsim = m_scene->getFluidSim();
          
          switch (script.type) {
            case Script::TRANSLATE:
            {
              if(DIM == 2) {
                FluidSim2D* fluid2d = (FluidSim2D*) fluidsim;
                const std::vector< FluidSim::SourceBoundary<2>* >& boundaries = fluid2d->get_sources();
                Vector2s dx = Vector2s(script.v(0), script.v(1));
                scalar ldx = dx.norm();
                if(script.index >= 0) {
                  if(boundaries[script.index]->type == FluidSim::BT_UNION || boundaries[script.index]->type == FluidSim::BT_INTERSECT) break;
                  FluidSim::SourceBoundary<2>* boundary = boundaries[script.index];
                  if(ldx > 0.0) {
                    Vector2s ndx = dx / ldx;
                    scalar vel = 0;
                    if(script.func == Script::COSINE) {
                      vel = cosine_ease_function(current_time + dt, script.start, script.end, script.start + script.ease_start, script.end - script.ease_end, script.amplitude, script.frequency);
                    } else if(script.func == Script::CUBIC) {
                      vel = cubic_ease_function(current_time + dt, script.start, script.end, script.start + script.ease_start, script.end - script.ease_end, ldx);
                    } else if(script.func == Script::WENO) {
                      vel = weno_ease_function(current_time + dt, dt, script.start, script.end, script.base_dt, script.base_pos, script.base_vertices);
                      script.base_pos += vel * dt;
                    }
                    boundary->future_center(0) += ndx(0) * vel * dt;
                    boundary->future_center(1) += ndx(1) * vel * dt;
                  }
                } else {
                  for(auto& b : boundaries) {
                    if(b->type == FluidSim::BT_UNION || b->type == FluidSim::BT_INTERSECT) continue;
                    
                    if(ldx > 0.0) {
                      Vector2s ndx = dx / ldx;
                      scalar vel = 0;
                      if(script.func == Script::COSINE) {
                        vel = cosine_ease_function(current_time + dt, script.start, script.end, script.start + script.ease_start, script.end - script.ease_end, script.amplitude, script.frequency);
                      } else if(script.func == Script::CUBIC) {
                        vel = cubic_ease_function(current_time + dt, script.start, script.end, script.start + script.ease_start, script.end - script.ease_end, ldx);
                      } else if(script.func == Script::WENO) {
                        vel = weno_ease_function(current_time + dt, dt, script.start, script.end, script.base_dt, script.base_pos, script.base_vertices);
                        script.base_pos += vel * dt;
                      }
                      b->future_center(0) += ndx(0) * vel * dt;
                      b->future_center(1) += ndx(1) * vel * dt;
                    }
                  }
                }
              } else {
                FluidSim3D* fluid3d = (FluidSim3D*) fluidsim;
                const std::vector< FluidSim::SourceBoundary<3>* >& boundaries = fluid3d->get_sources();
                Vector3s dx = Vector3s(script.v(0), script.v(1), script.v(2));
                scalar ldx = dx.norm();
                if(script.index >= 0) {
                  if(boundaries[script.index]->type == FluidSim::BT_UNION || boundaries[script.index]->type == FluidSim::BT_INTERSECT) break;
                  FluidSim::SourceBoundary<3>* boundary = boundaries[script.index];
                  if(ldx > 0.0) {
                    Vector3s ndx = dx / ldx;
                    scalar vel = 0;
                    if(script.func == Script::COSINE) {
                      vel = cosine_ease_function(current_time + dt, script.start, script.end, script.start + script.ease_start, script.end - script.ease_end, script.amplitude, script.frequency);
                    } else if(script.func == Script::CUBIC) {
                      vel = cubic_ease_function(current_time + dt, script.start, script.end, script.start + script.ease_start, script.end - script.ease_end, ldx);
                    } else if(script.func == Script::WENO) {
                      vel = weno_ease_function(current_time + dt, dt, script.start, script.end, script.base_dt, script.base_pos, script.base_vertices);
                      script.base_pos += vel * dt;
                    }
                    boundary->future_center(0) += ndx(0) * vel * dt;
                    boundary->future_center(1) += ndx(1) * vel * dt;
                    boundary->future_center(2) += ndx(2) * vel * dt;
                  }
                } else {
                  for(auto& b : boundaries) {
                    if(b->type == FluidSim::BT_UNION || b->type == FluidSim::BT_INTERSECT) continue;
                    
                    if(ldx > 0.0) {
                      Vector3s ndx = dx / ldx;
                      scalar vel = 0;
                      if(script.func == Script::COSINE) {
                        vel = cosine_ease_function(current_time + dt, script.start, script.end, script.start + script.ease_start, script.end - script.ease_end, script.amplitude, script.frequency);
                      } else if(script.func == Script::CUBIC) {
                        vel = cubic_ease_function(current_time + dt, script.start, script.end, script.start + script.ease_start, script.end - script.ease_end, ldx);
                      } else if(script.func == Script::WENO) {
                        vel = weno_ease_function(current_time + dt, dt, script.start, script.end, script.base_dt, script.base_pos, script.base_vertices);
                        script.base_pos += vel * dt;
                      }
                      b->future_center(0) += ndx(0) * vel * dt;
                      b->future_center(1) += ndx(1) * vel * dt;
                      b->future_center(2) += ndx(2) * vel * dt;
                    }
                  }
                }
              }
              updateSDF = updateSDF | script.updateSDF;
              
              break;
            }
            case Script::ROTATE:
            {
              if(DIM == 2) {
                FluidSim2D* fluid2d = (FluidSim2D*) fluidsim;
                const std::vector< FluidSim::SourceBoundary<2>* >& boundaries = fluid2d->get_sources();
                scalar vel = 0;
                if(script.func == Script::COSINE) {
                  vel = cosine_ease_function(current_time + dt, script.start, script.end, script.start + script.ease_start, script.end - script.ease_end, script.amplitude, script.frequency);
                } else if(script.func == Script::CUBIC) {
                  vel = cubic_ease_function(current_time + dt, script.start, script.end, script.start + script.ease_start, script.end - script.ease_end, script.v(3));
                } else if(script.func == Script::WENO) {
                  vel = weno_ease_function(current_time + dt, dt, script.start, script.end, script.base_dt, script.base_pos, script.base_vertices);
                  script.base_pos += vel * dt;
                }
                scalar da = vel * dt;
                Vector3s axis(0.0,0.0,1.0);
                
                Eigen::AngleAxis<scalar> rot(da, axis);
                Eigen::Quaternion<scalar> qrot(rot);
                
                if(script.index >= 0) {
                  if(boundaries[script.index]->type == FluidSim::BT_UNION || boundaries[script.index]->type == FluidSim::BT_INTERSECT) break;
                  FluidSim::SourceBoundary<2>* boundary = boundaries[script.index];
                  
                  Eigen::Quaternion<scalar>& prot = boundary->future_rot;
                  if(script.transform_global) {
                    prot = qrot * prot;
                  } else {
                    prot *= qrot;
                  }
                } else {
                  for(auto& b : boundaries)
                  {
                    if(b->type == FluidSim::BT_UNION || b->type == FluidSim::BT_INTERSECT) continue;
                    FluidSim::SourceBoundary<2>* boundary = b;
                    
                    Eigen::Quaternion<scalar>& prot = boundary->future_rot;
                    if(script.transform_global) {
                      prot = qrot * prot;
                    } else {
                      prot *= qrot;
                    }
                  }
                }
              } else {
                FluidSim3D* fluid3d = (FluidSim3D*) fluidsim;
                const std::vector< FluidSim::SourceBoundary<3>* >& boundaries = fluid3d->get_sources();
                
                scalar vel = 0;
                if(script.func == Script::COSINE) {
                  vel = cosine_ease_function(current_time + dt, script.start, script.end, script.start + script.ease_start, script.end - script.ease_end, script.amplitude, script.frequency);
                } else if(script.func == Script::CUBIC) {
                  vel = cubic_ease_function(current_time + dt, script.start, script.end, script.start + script.ease_start, script.end - script.ease_end, script.v(3));
                } else if(script.func == Script::WENO) {
                  vel = weno_ease_function(current_time + dt, dt, script.start, script.end, script.base_dt, script.base_pos, script.base_vertices);
                  script.base_pos += vel * dt;
                }
                scalar da = vel * dt;
                Vector3s axis(script.v(0),script.v(1),script.v(2));
                
                Eigen::AngleAxis<scalar> rot(da, axis);
                Eigen::Quaternion<scalar> qrot(rot);
                
                if(script.index >= 0) {
                  if(boundaries[script.index]->type == FluidSim::BT_UNION || boundaries[script.index]->type == FluidSim::BT_INTERSECT) break;
                  FluidSim::SourceBoundary<3>* boundary = boundaries[script.index];
                  Eigen::Quaternion<scalar>& prot = boundary->future_rot;
                  if(script.transform_global) {
                    prot = qrot * prot;
                  } else {
                    prot *= qrot;
                  }
                } else {
                  for(auto& b : boundaries)
                  {
                    if(b->type == FluidSim::BT_UNION || b->type == FluidSim::BT_INTERSECT) continue;
                    FluidSim::SourceBoundary<3>* boundary = b;
                    
                    Eigen::Quaternion<scalar>& prot = boundary->future_rot;
                    if(script.transform_global) {
                      prot = qrot * prot;
                    } else {
                      prot *= qrot;
                    }
                    
                  }
                }
              }
              
              updateSDF = updateSDF | script.updateSDF;
              break;
            }
            case Script::SCALE:
            {
              break;
            }
            default:
              std::cout << "UNIMPLEMENTED SCRIPT TYPE [" << script.type << "]!" << std::endl;
              break;
          }
          break;
        }
        default:
          std::cout << "UNIMPLEMENTED SCRIPT TARGET [" << script.target << "]!" << std::endl;
          break;
      }
    }
  }
  
  return updateSDF;
}

// explicit instantiations at bottom
template class ParticleSimulation<2>;
template class ParticleSimulation<3>;
