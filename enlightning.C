/*
   ===========================================================================
   Copyright (C) 2017 Jon Rood.

   This file is part of Enlightning source code.

   Enlightning source code is free software; you can redistribute it
   and/or modify it under the terms of the GNU General Public License as
   published by the Free Software Foundation; either version 3 of the License,
   or (at your option) any later version.

   Enlightning source code is distributed in the hope that it will be
   useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with Enlightning; if not, see <http://www.gnu.org/licenses/>.
   ===========================================================================
 */

#include "enlightning.h"

#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cfloat>

#include "SAMRAI/hier/BoundaryBox.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellIndex.h"
#include "SAMRAI/pdat/CellIterator.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/RestartManager.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/Timer.h"

#include "enlightningfort.h"

// Fast define for squaring a number.
#define SQUARE(x) ((x) * (x))

// Source geometry defines.
#define POINT (10)
#define SINE (20)
#define LIGHTNING (30)

// Source action defines.
#define PULSE (10)
#define DRIVEN (20)

// Boundary condition defines.
#define REFLECT (10)
#define ABSORB (20)

// Number of ghost cells.
#define CELLG (3)

#define PI 3.14159265358979323846264338

// Timers.
boost::shared_ptr<tbox::Timer> enlightning::t_init;
boost::shared_ptr<tbox::Timer> enlightning::t_rhs;
boost::shared_ptr<tbox::Timer> enlightning::t_rk;
boost::shared_ptr<tbox::Timer> enlightning::t_state;
boost::shared_ptr<tbox::Timer> enlightning::t_absorb;
boost::shared_ptr<tbox::Timer> enlightning::t_bc;
boost::shared_ptr<tbox::Timer> enlightning::t_tag;
boost::shared_ptr<tbox::Timer> enlightning::t_addsource;
boost::shared_ptr<tbox::Timer> enlightning::t_record;

// Constructor to intialize simulation variables.
enlightning::enlightning(
  const string                                & object_name,
  const tbox::Dimension                       & dim,
  boost::shared_ptr<tbox::Database>             input_db,
  boost::shared_ptr<geom::CartesianGridGeometry>grid_geom)
  :
  algs::MethodOfLinesPatchStrategy::MethodOfLinesPatchStrategy(),
  p_object_name(object_name),
  p_grid_geometry(grid_geom),
  p_dim(dim),
  p_nghosts(dim, CELLG),

// p_nghosts(hier::IntVector(dim, CELLG)),
  p_zero_ghosts(p_dim, 0),

// Allocate grid cell variables.
  p_w_n(new pdat::CellVariable<double>(dim, "w_n", NEQU)),
  p_K(new pdat::CellVariable<double>(dim, "K", NEQU)),
  p_pressure(new pdat::CellVariable<double>(dim, "pressure", 1)),
  p_temperature(new pdat::CellVariable<double>(dim, "temperature", 1)),
  p_source(new pdat::CellVariable<double>(dim, "source", 1)),

// Initialize variables to NaN, etc.
  p_src_type(tbox::MathUtilities<string>::getMax()),
  p_src_type_int(tbox::MathUtilities<int>::getSignalingNaN()),
  p_src_mode(tbox::MathUtilities<string>::getMax()),
  p_src_mode_int(tbox::MathUtilities<int>::getSignalingNaN()),
  p_cfl(tbox::MathUtilities<double>::getSignalingNaN()),
  p_specific_dt(tbox::MathUtilities<double>::getSignalingNaN()),
  p_start_time(tbox::MathUtilities<double>::getSignalingNaN()),
  p_end_time(tbox::MathUtilities<double>::getSignalingNaN()),
  p_loop_time(tbox::MathUtilities<double>::getSignalingNaN()),
  p_max_timesteps(tbox::MathUtilities<int>::getSignalingNaN()),
  p_regrid_step(tbox::MathUtilities<int>::getSignalingNaN()),
  p_tag_buffer(tbox::MathUtilities<int>::getSignalingNaN()),
  p_viz_pressure(tbox::MathUtilities<int>::getSignalingNaN()),
  p_viz_source(tbox::MathUtilities<int>::getSignalingNaN()),
  p_viz_temperature(tbox::MathUtilities<int>::getSignalingNaN()),
  p_viz_w_n(tbox::MathUtilities<int>::getSignalingNaN()),
  p_iteration_number(tbox::MathUtilities<int>::getSignalingNaN()),
  p_bc_type_L(tbox::MathUtilities<string>::getMax()),
  p_bc_type_R(tbox::MathUtilities<string>::getMax()),
  p_bc_type_B(tbox::MathUtilities<string>::getMax()),
  p_bc_type_T(tbox::MathUtilities<string>::getMax()),
  p_bc_type_L_int(tbox::MathUtilities<int>::getSignalingNaN()),
  p_bc_type_R_int(tbox::MathUtilities<int>::getSignalingNaN()),
  p_bc_type_B_int(tbox::MathUtilities<int>::getSignalingNaN()),
  p_bc_type_T_int(tbox::MathUtilities<int>::getSignalingNaN()),
  p_bc_absorb_width(tbox::MathUtilities<double>::getSignalingNaN()),
  p_src_amp(tbox::MathUtilities<double>::getSignalingNaN()),
  p_src_spread(tbox::MathUtilities<double>::getSignalingNaN()),
  p_src_sine_amp(tbox::MathUtilities<double>::getSignalingNaN()),
  p_src_sine_freq(tbox::MathUtilities<double>::getSignalingNaN()),
  p_src_omega(tbox::MathUtilities<double>::getSignalingNaN()),
  p_src_num_pulses(tbox::MathUtilities<int>::getSignalingNaN()),
  p_p0(tbox::MathUtilities<double>::getSignalingNaN()),
  p_sampled_p0(tbox::MathUtilities<double>::getSignalingNaN()),
  p_T0(tbox::MathUtilities<double>::getSignalingNaN()),
  p_c_lr(tbox::MathUtilities<double>::getSignalingNaN()),
  p_rho0(tbox::MathUtilities<double>::getSignalingNaN()),
  p_c(tbox::MathUtilities<double>::getSignalingNaN()),
  p_gamma1(tbox::MathUtilities<double>::getSignalingNaN()),
  p_RH(tbox::MathUtilities<double>::getSignalingNaN()),
  p_R(tbox::MathUtilities<double>::getSignalingNaN()),
  p_c_p(tbox::MathUtilities<double>::getSignalingNaN()),
  p_T_star_n(tbox::MathUtilities<double>::getSignalingNaN()),
  p_T_star_o(tbox::MathUtilities<double>::getSignalingNaN()),
  p_mu(tbox::MathUtilities<double>::getSignalingNaN()),
  p_muB(tbox::MathUtilities<double>::getSignalingNaN()),
  p_kappa(tbox::MathUtilities<double>::getSignalingNaN()),
  p_weno_alpha(tbox::MathUtilities<double>::getSignalingNaN()),
  p_tau_n(tbox::MathUtilities<double>::getSignalingNaN()),
  p_tau_o(tbox::MathUtilities<double>::getSignalingNaN()),
  p_R_tilde(tbox::MathUtilities<double>::getSignalingNaN()),
  p_c_v(tbox::MathUtilities<double>::getSignalingNaN()),
  p_M(tbox::MathUtilities<double>::getSignalingNaN()),
  p_record_audio(tbox::MathUtilities<int>::getSignalingNaN()),
  p_num_mics(tbox::MathUtilities<int>::getSignalingNaN()),
  p_use_sample_rate_dt(tbox::MathUtilities<int>::getSignalingNaN()),
  p_wav_sample_rate(tbox::MathUtilities<int>::getSignalingNaN()),
  p_lightning_input_line_count(tbox::MathUtilities<int>::getSignalingNaN())
{
  TBOX_ASSERT(!object_name.empty());
  TBOX_ASSERT(input_db);
  TBOX_ASSERT(grid_geom);

  tbox::MathUtilities<double>::setArrayToSignalingNaN(p_tag_tolerance, 3),
  tbox::MathUtilities<double>::setArrayToSignalingNaN(p_src_pos, 2),
  p_src_pulse_times.resize(0);
  p_mic_pos_x.resize(0);
  p_mic_pos_y.resize(0);
  p_lightning_coords.resize(0);

  // Register timers.
  if (!t_init) {
    t_init = tbox::TimerManager::getManager()->
                  getTimer("apps::enlightning::initializeDataOnPatch()");
    t_rhs = tbox::TimerManager::getManager()->
                  getTimer("apps::enlightning::singleStep()_rhs");
    t_rk = tbox::TimerManager::getManager()->
                  getTimer("apps::enlightning::singleStep()_rk");
    t_state = tbox::TimerManager::getManager()->
                  getTimer("apps::enlightning::singleStep()_state");
    t_absorb = tbox::TimerManager::getManager()->
                  getTimer("apps::enlightning::setAbsorbingBoundaryConditions()");
    t_bc = tbox::TimerManager::getManager()->
                  getTimer("apps::enlightning::setPhysicalBoundaryConditions()");
    t_tag = tbox::TimerManager::getManager()->
                  getTimer("apps::enlightning::tagGradientDetectorCells()");
    t_addsource = tbox::TimerManager::getManager()->
                  getTimer("apps::enlightning::addSource()");
    t_record = tbox::TimerManager::getManager()->
                  getTimer("apps::enlightning::recordAudio()");
  }

  // Register the enlightning object with the restart manager.
  tbox::RestartManager::getManager()->registerRestartItem(p_object_name, this);

  // If restarting load data from checkpoint, otherwise load initial conditions.
  bool is_from_restart = tbox::RestartManager::getManager()->isFromRestart();

  if (is_from_restart) {
    getFromRestart();
  } else {
    getFromInput(input_db);

    // The ambient pressure deviates a bit after the first timestep,
    // so p_sampled_p0 will use a sampled pressure value as the ambient
    // value when recording audio with virtual microphones.
    p_sampled_p0       = p_p0;
    p_loop_time        = p_start_time;
    p_iteration_number = 0;

    // Set simulation constants.
    const double T_ref = p_T0;
    const double p_ref = p_p0;
    const double zeta  = 10.79586 * (1 - (273.16 / p_T0)) - 5.02808
                         * log10(p_T0 / 273.16) + (1.50474e-4)
                         * (1 - pow(10, (-8.29692 * ((p_T0 / 273.16) - 1))))
                         - (4.2873e-4) *
                         (1 - pow(10, (-4.76955 * ((273.16 / p_T0) - 1))))
                         - 2.2195983;
    const double p_vp = p_ref * pow(10, zeta);
    const double h    = (p_RH * p_vp) / p_p0;
    const double eta  = 4.17 * (pow(T_ref / p_T0, 1.0 / 3.0) - 1);
    const double f_n  = (p_p0 / p_ref) * (pow((T_ref / p_T0), 0.5)
                                          * (9 + 280 * h * exp(-eta)));
    const double f_o = (p_p0 / p_ref) * (24 + (4.04e4) * h
                                         * ((0.02 + h) / (0.391 + h)));
    p_tau_n = 1 / (2 * PI * f_n);
    p_tau_o = 1 / (2 * PI * f_o);

    // Read in lightning channel geometry if necessary.
    if (p_src_type_int == LIGHTNING) readLightningInput("lightning.txt");
  }
}

// Destructor only sets timers to NULL.
enlightning::~enlightning()
{
  t_init.reset();
  t_rhs.reset();
  t_rk.reset();
  t_state.reset();
  t_absorb.reset();
  t_bc.reset();
  t_tag.reset();
  t_addsource.reset();
  t_record.reset();
}

// Register cell variables with SAMRAI.
void enlightning::registerModelVariables(algs::MethodOfLinesIntegrator *integrator)
{
  // Cell variables that need interpolation.
  integrator->registerVariable(p_w_n, p_nghosts,
                               algs::MethodOfLinesIntegrator::SOLN,
                               p_grid_geometry,
                               "CONSERVATIVE_COARSEN",
                               "LINEAR_REFINE");

  integrator->registerVariable(p_pressure, p_nghosts,
                               algs::MethodOfLinesIntegrator::SOLN,
                               p_grid_geometry,
                               "CONSERVATIVE_COARSEN",
                               "LINEAR_REFINE");

  integrator->registerVariable(p_temperature, p_nghosts,
                               algs::MethodOfLinesIntegrator::SOLN,
                               p_grid_geometry,
                               "CONSERVATIVE_COARSEN",
                               "LINEAR_REFINE");

  integrator->registerVariable(p_source, p_nghosts,
                               algs::MethodOfLinesIntegrator::SOLN,
                               p_grid_geometry,
                               "CONSERVATIVE_COARSEN",
                               "LINEAR_REFINE");

  // Scratch variable that doesn't need interpolation.
  integrator->registerVariable(p_K, p_nghosts,
                               algs::MethodOfLinesIntegrator::RHS,
                               p_grid_geometry,
                               "NO_COARSEN",
                               "NO_REFINE");

  // Object to get data from the cell variables.
  hier::VariableDatabase *vardb = hier::VariableDatabase::getDatabase();

  int grid_id;

  // Output visualization for grids in the w_n vector.
  if (p_viz_w_n == 1) {
    grid_id = vardb->mapVariableAndContextToIndex(p_w_n, getInteriorContext());
    string dump_name = "w_n #";
    const int size   = static_cast<int>(dump_name.length()) + 16;
    char     *buffer = new char[size];

    for (int n = 0; n < NEQU; n++) {
      sprintf(buffer, "%s%01d", dump_name.c_str(), n);
      string variable_name(buffer);

      if (p_visit_writer) {
        p_visit_writer->registerPlotQuantity(variable_name,
                                             "SCALAR",
                                             grid_id, n);
      }

      if (!p_visit_writer) {
        TBOX_WARNING(p_object_name
                     << ": registerModelVariables()\n"
                     << "Visit data writer was not registered.\n"
                     << "Consequently, no plot data will\n"
                     << "be written."
                     << endl);
      }
    }
    delete[] buffer;
  }

  // Output visualization for the pressure grid.
  if (p_viz_pressure == 1) {
    grid_id =
      vardb->mapVariableAndContextToIndex(p_pressure, getInteriorContext());

    if (p_visit_writer) {
      p_visit_writer->registerPlotQuantity("pressure",
                                           "SCALAR",
                                           grid_id,
                                           0);
    } else {
      TBOX_WARNING(p_object_name
                   << ": registerModelVariables()\n"
                   << "Visit data writer was not registered.\n"
                   << "Consequently, no plot data will\n"
                   << "be written."
                   << endl);
    }
  }

  // Output visualization for the temperature grid.
  if (p_viz_temperature == 1) {
    grid_id = vardb->mapVariableAndContextToIndex(p_temperature,
                                                  getInteriorContext());

    if (p_visit_writer) {
      p_visit_writer->registerPlotQuantity("temperature",
                                           "SCALAR",
                                           grid_id, 0);
    } else {
      TBOX_WARNING(p_object_name
                   << ": registerModelVariables()\n"
                   << "Visit data writer was not registered.\n"
                   << "Consequently, no plot data will\n"
                   << "be written."
                   << endl);
    }
  }

  // Output visualization for the source grid.
  if (p_viz_source == 1) {
    grid_id = vardb->mapVariableAndContextToIndex(p_source, getInteriorContext());

    if (p_visit_writer) {
      p_visit_writer->registerPlotQuantity("source", "SCALAR",
                                           grid_id, 0);
    } else {
      TBOX_WARNING(p_object_name
                   << ": registerModelVariables()\n"
                   << "Visit data writer was not registered.\n"
                   << "Consequently, no plot data will\n"
                   << "be written."
                   << endl);
    }
  }
}

// Initialize each patch each timestep.
void enlightning::initializeDataOnPatch(hier::Patch& patch,
                                        const double time,
                                        const bool   initial_time) const
{
  (void)time;
  (void)initial_time;

  // Start timer for initialization.
  t_init->start();

  // Get source grid object.
  boost::shared_ptr<pdat::CellData<double> > source_grid(
    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
      patch.getPatchData(p_source, getInteriorContext())));

  // Set source to 0 every timestep.
  source_grid->fillAll(0);

  // Use initial conditions at first timestep.
  if (initial_time) {
    const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
      BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
        patch.getPatchGeometry()));

    boost::shared_ptr<pdat::CellData<double> > w_n(
      BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
        patch.getPatchData(p_w_n, getInteriorContext())));

    boost::shared_ptr<pdat::CellData<double> > pressure_grid(
      BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
        patch.getPatchData(p_pressure, getInteriorContext())));

    boost::shared_ptr<pdat::CellData<double> > temperature_grid(
      BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
        patch.getPatchData(p_temperature, getInteriorContext())));

    pressure_grid->fillAll(p_p0);
    temperature_grid->fillAll(p_T0);
    w_n->fill(p_rho0,        0);
    w_n->fill(0,             1);
    w_n->fill(0,             2);
    w_n->fill(0,             3);
    w_n->fill(p_rho0 * p_T0, 4);
    // w_n->fill(p_rho0*p_T0,5);
  }

  // Stop timer for initialization.
  t_init->stop();
}

// Decide what timestep to use on each patch.
double enlightning::computeStableDtOnPatch(hier::Patch& patch,
                                           const double time) const
{
  (void)time;
  (void)patch;

  const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
      patch.getPatchGeometry()));
  const double *dx = patch_geom->getDx();
  double dt        = 0.0;

  // Compute dt.
  if (dx[0] <= dx[1]) dt = p_cfl * (dx[0] / p_c);
  else dt = p_cfl * (dx[1] / p_c);

  // Use a specific dt.
  if (p_specific_dt > 0) dt = p_specific_dt;

  // Use a .wav sample rate dt.
  if (p_use_sample_rate_dt == 1) dt = 1.0 / p_wav_sample_rate;

  return dt;
}

// Perform the calculation on a patch over a timestep.
void enlightning::singleStep(hier::Patch& patch,
                             const double dt,
                             const double alpha_1,
                             const double alpha_2,
                             const double beta) const
{
  // Get pointer to w_n grid for next timestep.
  boost::shared_ptr<pdat::CellData<double> > w_n(
    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
      patch.getPatchData(p_w_n, getInteriorWithGhostsContext())));

  // Get pointer to w grid for this timestep.
  boost::shared_ptr<pdat::CellData<double> > w(
    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
      patch.getPatchData(p_w_n, getInteriorContext())));

  boost::shared_ptr<pdat::CellData<double> > source_grid(
    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
      patch.getPatchData(p_source, getInteriorContext())));

  // Scratch grid for calculating the RHS.
  boost::shared_ptr<pdat::CellData<double> > K(
    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
      patch.getPatchData(p_K, getInteriorContext())));

  boost::shared_ptr<pdat::CellData<double> > pressure_grid(
    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
      patch.getPatchData(p_pressure, getInteriorWithGhostsContext())));

  boost::shared_ptr<pdat::CellData<double> > temperature_grid(
    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
      patch.getPatchData(p_temperature, getInteriorWithGhostsContext())));

  // Get information for this patch for passing to fortran.
  const hier::Index ifirst = patch.getBox().lower();
  const hier::Index ilast  = patch.getBox().upper();

  const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
      patch.getPatchGeometry()));

  const double *dx       = patch_geom->getDx();
  const double *patchXLo = patch_geom->getXLower();

  const int level = patch.getPatchLevelNumber();

  // Call fortran function for calculating the right-hand side of the equation.
  t_rhs->start();
  SAMRAI_F77_FUNC(rhs, RHS) (ifirst(0),
                             ilast(0),
                             ifirst(1),
                             ilast(1),
                             p_nghosts(0),
                             p_nghosts(1),
                             patchXLo,
                             dx,
                             w_n->getPointer(),
                             source_grid->getPointer(),
                             K->getPointer(),
                             pressure_grid->getPointer(),
                             temperature_grid->getPointer(),
                             p_weno_alpha,
                             p_kappa,
                             p_mu,
                             p_muB,
                             p_R,
                             p_T_star_n,
                             p_T_star_o,
                             p_tau_n,
                             p_tau_o,
                             p_tag_tolerance,
                             level,
                             p_c,
                             p_c_lr,
                             NEQU);
  t_rhs->stop();

  // Call fortran function to perform the runge kutta.
  t_rk->start();
  SAMRAI_F77_FUNC(runge_kutta, RUNGE_KUTTA) (ifirst(0),
                                             ilast(0),
                                             ifirst(1),
                                             ilast(1),
                                             p_nghosts(0),
                                             p_nghosts(1),
                                             dt,
                                             alpha_1,
                                             alpha_2,
                                             beta,
                                             w_n->getPointer(),
                                             w->getPointer(),
                                             K->getPointer(),
                                             NEQU);
  t_rk->stop();

  // Call fortran function to update the state variables.
  t_state->start();
  SAMRAI_F77_FUNC(update_state, UPDATE_STATE) (dx,
                                               patchXLo,
                                               ifirst(0),
                                               ilast(0),
                                               ifirst(1),
                                               ilast(1),
                                               p_nghosts(0),
                                               p_nghosts(1),
                                               w_n->getPointer(),
                                               pressure_grid->getPointer(),
                                               temperature_grid->getPointer(),
                                               p_rho0,
                                               p_c,
                                               p_gamma1,
                                               p_c_p,
                                               p_p0,
                                               p_T0,
                                               p_c_lr,
                                               p_c_v,
                                               p_R_tilde,
                                               p_M,
                                               p_R,
                                               NEQU);
  t_state->stop();

  // Set absorbing boundary conditions if necessary.
  setAbsorbingBoundaryConditions(patch);
}

// Tag cells on patch that are not smooth.
void enlightning::tagGradientDetectorCells(hier::Patch& patch,
                                           const double regrid_time,
                                           const bool   initial_error,
                                           const int    tag_index,
                                           const bool   uses_richardson_extrapolation_too)
{
  (void)regrid_time;
  (void)initial_error;
  (void)uses_richardson_extrapolation_too;

  // Get pointers to patch data to send to fortran functions.
  boost::shared_ptr<pdat::CellData<int> > tags(
    BOOST_CAST<pdat::CellData<int>, hier::PatchData>(
      patch.getPatchData(tag_index)));

  // Only use w_n and source grids for smoothness detection.
  boost::shared_ptr<pdat::CellData<double> > w_n(
    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
      patch.getPatchData(p_w_n, getInteriorWithGhostsContext())));

  boost::shared_ptr<pdat::CellData<double> > source_grid(
    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
      patch.getPatchData(p_source, getInteriorWithGhostsContext())));

  const hier::Index ifirst = patch.getBox().lower();
  const hier::Index ilast  = patch.getBox().upper();

  // Call fortran function to tag cells that are not smooth.
  t_tag->start();
  SAMRAI_F77_FUNC(tag_cells, TAG_CELLS) (ifirst(0),
                                         ilast(0),
                                         ifirst(1),
                                         ilast(1),
                                         p_nghosts(0),
                                         p_nghosts(1),
                                         tags->getPointer(),
                                         w_n->getPointer(),
                                         source_grid->getPointer(),
                                         true,
                                         p_tag_tolerance,
                                         NEQU);
  t_tag->stop();
}

// Handle absorbing boundary conditions. The setPhysicalBoundaryConditions()
// function is only called on patches that touch the physical boundary,
// therefore this function is always called to apply the absorbing
// envelopes to all patches if absorbing boundary conditions are used.
void enlightning::setAbsorbingBoundaryConditions(hier::Patch& patch) const
{
  t_absorb->start();

  // Any of the 4 sides of the domain can be set to use an absorbing boundary.
  if ((p_bc_type_L_int == ABSORB) || (p_bc_type_R_int == ABSORB)
      || (p_bc_type_B_int == ABSORB) || (p_bc_type_T_int == ABSORB)) {
    boost::shared_ptr<pdat::CellData<double> > w_n(
      BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
        patch.getPatchData(p_w_n, getInteriorWithGhostsContext())));

    boost::shared_ptr<pdat::CellData<double> > pressure_grid(
      BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
        patch.getPatchData(p_pressure, getInteriorWithGhostsContext())));
    boost::shared_ptr<pdat::CellData<double> > temperature_grid(
      BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
        patch.getPatchData(p_temperature, getInteriorWithGhostsContext())));
    const hier::Index ifirst = patch.getBox().lower();
    const hier::Index ilast  = patch.getBox().upper();
    const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
      BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
        patch.getPatchGeometry()));
    const double *dx                  = patch_geom->getDx();
    const double *gridXLo             = p_grid_geometry->getXLower();
    const double *gridXHi             = p_grid_geometry->getXUpper();
    double *w_n_data                  = w_n->getPointer();
    double *pressure_data             = pressure_grid->getPointer();
    double *temperature_data          = temperature_grid->getPointer();
    const hier::IntVector ghost_cells = w_n->getGhostCellWidth();
    const int patch_width             =
      (ilast(0) + 1) - ifirst(0) + 2 * ghost_cells(0);
    const int patch_area              =
      ((ilast(1) + 1) - ifirst(1) + 2 * ghost_cells(1)) * patch_width;
    double xc[2];
    const double *patchXLo = patch_geom->getXLower();
    double tmp;
    int    idx2d;
    double ep = 1.0e-8;

    // Apply envelope to absorb waves at the left boundary.
    if (p_bc_type_L_int == ABSORB) {
      xc[0] = patchXLo[0] + dx[0] * (0.5);
      tmp = -exp(-log(2.0) / SQUARE(p_bc_absorb_width) * SQUARE(xc[0] - gridXLo[0])) + 1;

      for (int j = ifirst(1); j <= ilast(1); j++) {
        xc[1] = patchXLo[1] + dx[1] * ((double)(j - ifirst(1)) + 0.5);

        for (int i = ifirst(0); (fabs(1.0 - tmp) > ep) && (i <= ilast(0)); i++) {
          xc[0] = patchXLo[0] + dx[0] * ((double)(i - ifirst(0)) + 0.5);
          tmp = -exp(-log(2.0) / SQUARE(p_bc_absorb_width) * SQUARE(xc[0] - gridXLo[0])) + 1;
          idx2d = (j - ifirst(1) + ghost_cells(1)) * patch_width +
                  (i - ifirst(0) + ghost_cells(0));
          pressure_data[idx2d] = (pressure_data[idx2d] - p_p0) * tmp + p_p0;
          temperature_data[idx2d] = (temperature_data[idx2d] - p_T0) * tmp + p_T0;
          w_n_data[0 * patch_area + idx2d] =
            (w_n_data[0 * patch_area + idx2d] - p_rho0) * tmp + p_rho0;
          w_n_data[1 * patch_area + idx2d] = w_n_data[1 * patch_area + idx2d] * tmp;
          w_n_data[2 * patch_area + idx2d] = w_n_data[2 * patch_area + idx2d] * tmp;
          w_n_data[3 * patch_area + idx2d] = w_n_data[3 * patch_area + idx2d] * tmp;
          w_n_data[4 * patch_area + idx2d] = (w_n_data[4 * patch_area + idx2d] 
			                   - p_T0 * p_rho0) * tmp + p_T0 * p_rho0;
          // w_n_data[5*patch_area + idx2d] =
          // (w_n_data[5*patch_area + idx2d]-p_T0*p_rho0)*tmp+p_T0*p_rho0;
        }
        tmp = 0;
      }
    }

    // Apply envelope to absorb waves at the right boundary.
    if (p_bc_type_R_int == ABSORB) {
      xc[0] = patchXLo[0] + dx[0] * ((double)(ilast(0) - ifirst(0)) + 0.5);
      tmp = -exp(-log(2.0) / SQUARE(p_bc_absorb_width) * SQUARE(gridXHi[0] - xc[0])) + 1;

      for (int j = ifirst(1); j <= ilast(1); j++) {
        xc[1] = patchXLo[1] + dx[1] * ((double)(j - ifirst(1)) + 0.5);

        for (int i = ilast(0); (fabs(1.0 - tmp) > ep) && (i >= ifirst(0)); i--) {
          xc[0] = patchXLo[0] + dx[0] * ((double)(i - ifirst(0)) + 0.5);
          tmp = -exp(-log(2.0) / SQUARE(p_bc_absorb_width) * SQUARE(gridXHi[0] - xc[0])) + 1;
          idx2d = (j - ifirst(1) + ghost_cells(1)) * patch_width +
                  (i - ifirst(0) + ghost_cells(0));
          pressure_data[idx2d] = (pressure_data[idx2d] - p_p0) * tmp + p_p0;
          temperature_data[idx2d] = (temperature_data[idx2d] - p_T0) * tmp + p_T0;
          w_n_data[0 * patch_area + idx2d] =
            (w_n_data[0 * patch_area + idx2d] - p_rho0) * tmp + p_rho0;
          w_n_data[1 * patch_area + idx2d] = w_n_data[1 * patch_area + idx2d] * tmp;
          w_n_data[2 * patch_area + idx2d] = w_n_data[2 * patch_area + idx2d] * tmp;
          w_n_data[3 * patch_area + idx2d] = w_n_data[3 * patch_area + idx2d] * tmp;
          w_n_data[4 * patch_area + idx2d] = (w_n_data[4 * patch_area + idx2d] 
			                   - p_T0 * p_rho0) * tmp + p_T0 * p_rho0;
          // w_n_data[5*patch_area + idx2d] =
          // (w_n_data[5*patch_area + idx2d]-p_T0*p_rho0)*tmp+p_T0*p_rho0;
        }
        tmp = 0;
      }
    }

    // Apply envelope to absorb waves at the bottom boundary.
    if (p_bc_type_B_int == ABSORB) {
      xc[1] = patchXLo[1] + dx[1] * (0.5);
      tmp = -exp(-log(2.0) / SQUARE(p_bc_absorb_width) * SQUARE(xc[1] - gridXLo[1])) + 1;

      for (int j = ifirst(1); (fabs(1.0 - tmp) > ep) && (j <= ilast(1)); j++) {
        xc[1] = patchXLo[1] + dx[1] * ((double)(j - ifirst(1)) + 0.5);
        tmp = -exp(-log(2.0) / SQUARE(p_bc_absorb_width) * SQUARE(xc[1] - gridXLo[1])) + 1;

        for (int i = ifirst(0); i <= ilast(0); i++) {
          idx2d = (j - ifirst(1) + ghost_cells(1)) * patch_width +
                  (i - ifirst(0) + ghost_cells(0));
          pressure_data[idx2d] = (pressure_data[idx2d] - p_p0) * tmp + p_p0;
          temperature_data[idx2d] = (temperature_data[idx2d] - p_T0) * tmp + p_T0;
          w_n_data[0 * patch_area + idx2d] =
            (w_n_data[0 * patch_area + idx2d] - p_rho0) * tmp + p_rho0;
          w_n_data[1 * patch_area + idx2d] = w_n_data[1 * patch_area + idx2d] * tmp;
          w_n_data[2 * patch_area + idx2d] = w_n_data[2 * patch_area + idx2d] * tmp;
          w_n_data[3 * patch_area + idx2d] = w_n_data[3 * patch_area + idx2d] * tmp;
          w_n_data[4 * patch_area + idx2d] = (w_n_data[4 * patch_area + idx2d]
			                   - p_T0 * p_rho0) * tmp + p_T0 * p_rho0;
          // w_n_data[5*patch_area + idx2d] =
          // (w_n_data[5*patch_area + idx2d]-p_T0*p_rho0)*tmp+p_T0*p_rho0;
        }
      }
    }

    // Apply envelope to absorb waves at the top boundary.
    if (p_bc_type_T_int == ABSORB) {
      xc[1] = patchXLo[1] + dx[1] * ((double)(ilast(1) - ifirst(1)) + 0.5);
      tmp = -exp(-log(2.0) / SQUARE(p_bc_absorb_width) * SQUARE(gridXHi[1] - xc[1])) + 1;

      for (int j = ilast(1); (fabs(1.0 - tmp) > ep) && (j >= ifirst(1)); j--) {
        xc[1] = patchXLo[1] + dx[1] * ((double)(j - ifirst(1)) + 0.5);
        tmp = -exp(-log(2.0) / SQUARE(p_bc_absorb_width) * SQUARE(gridXHi[1] - xc[1])) + 1;

        for (int i = ifirst(0); i <= ilast(0); i++) {
          idx2d = (j - ifirst(1) + ghost_cells(1)) * patch_width +
                  (i - ifirst(0) + ghost_cells(0));
          pressure_data[idx2d] = (pressure_data[idx2d] - p_p0) * tmp + p_p0;
          temperature_data[idx2d] = (temperature_data[idx2d] - p_T0) * tmp + p_T0;
          w_n_data[0 * patch_area + idx2d] =
            (w_n_data[0 * patch_area + idx2d] - p_rho0) * tmp + p_rho0;
          w_n_data[1 * patch_area + idx2d] = w_n_data[1 * patch_area + idx2d] * tmp;
          w_n_data[2 * patch_area + idx2d] = w_n_data[2 * patch_area + idx2d] * tmp;
          w_n_data[3 * patch_area + idx2d] = w_n_data[3 * patch_area + idx2d] * tmp;
          w_n_data[4 * patch_area + idx2d] = (w_n_data[4 * patch_area + idx2d]
			                   - p_T0 * p_rho0) * tmp + p_T0 * p_rho0;
          // w_n_data[5*patch_area + idx2d] =
          // (w_n_data[5*patch_area + idx2d]-p_T0*p_rho0)*tmp+p_T0*p_rho0;
        }
      }
    }
  }
  t_absorb->stop();
}

// Handle reflecting boundary conditions. Only called when a patch
// touches the physical boundary of the domain.
void enlightning::setPhysicalBoundaryConditions(hier::Patch          & patch,
                                                const double           fill_time,
                                                const hier::IntVector& ghost_width_to_fill)
{
  (void)fill_time;

  t_bc->start();

  const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
      patch.getPatchGeometry()));

  // if (patch_geom->getTouchesRegularBoundary() == true) {
  const hier::Index ifirstP = patch.getBox().lower();
  const hier::Index ilastP  = patch.getBox().upper();
  boost::shared_ptr<pdat::CellData<double> > w_n(
    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
      patch.getPatchData(p_w_n, getInteriorWithGhostsContext())));
  boost::shared_ptr<pdat::CellData<double> > pressure_grid(
    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
      patch.getPatchData(p_pressure, getInteriorWithGhostsContext())));
  boost::shared_ptr<pdat::CellData<double> > temperature_grid(
    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
      patch.getPatchData(p_temperature, getInteriorWithGhostsContext())));
  boost::shared_ptr<pdat::CellData<double> > source_grid(
    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
      patch.getPatchData(p_source, getInteriorWithGhostsContext())));

  // Get pointers to cell data and information about the patch.
  double *w_n_data                  = w_n->getPointer();
  double *pressure_data             = pressure_grid->getPointer();
  double *temperature_data          = temperature_grid->getPointer();
  double *source_data               = source_grid->getPointer();
  const hier::IntVector ghost_cells = w_n->getGhostCellWidth();
  const int patch_width             = (ilastP(0) + 1) - ifirstP(0) + 2 *
                                      ghost_cells(0);
  const int patch_area =
    ((ilastP(1) + 1) - ifirstP(1) + 2 * ghost_cells(1)) * patch_width;

  // Iterate over patch edges that touch the boundary.
  const std::vector<hier::BoundaryBox> edge_bdry =
    patch_geom->getEdgeBoundaries();

  for (int q = 0; q < edge_bdry.size(); q++) {
    int bdry_loc = edge_bdry[q].getLocationIndex();
    hier::Box interior(patch.getBox());
    hier::Box fill_box =
      patch_geom->getBoundaryFillBox(edge_bdry[q], interior, ghost_width_to_fill);
    const hier::Index ifirst = fill_box.lower();
    const hier::Index ilast  = fill_box.upper();

    // Apply reflection to left boundary.
    if (bdry_loc == 0) {
      for (int j = ifirst(1); j <= ilast(1); j++) {
        for (int i = ilast(0); i >= ifirst(0); i--) {
          int idx2d =
            (j - ifirstP(1) + ghost_cells(1)) * patch_width +
            (i  - ifirstP(0) + ghost_cells(0));
          int idx2dIP1 =
            (j - ifirstP(1) + ghost_cells(1)) * patch_width +
            (i + 1 - ifirstP(0) + ghost_cells(0));
          pressure_data[idx2d]    = pressure_data[idx2dIP1];
          temperature_data[idx2d] = temperature_data[idx2dIP1];
          source_data[idx2d]      = source_data[idx2dIP1];

          for (int k = 0; k < NEQU; k++) {
            int idx3d    = k * patch_area + idx2d;
            int idx3dIP1 = k * patch_area + idx2dIP1;
            w_n_data[idx3d] = w_n_data[idx3dIP1];

            if ((k == 1) && (i == ilast(0)) &&
                (p_bc_type_L_int == REFLECT)) w_n_data[idx3d] = -w_n_data[idx3d];
          }
        }
      }
    }

    // Apply reflection to right boundary.
    if (bdry_loc == 1) {
      for (int j = ifirst(1); j <= ilast(1); j++) {
        for (int i = ifirst(0); i <= ilast(0); i++) {
          int idx2d =
            (j - ifirstP(1) + ghost_cells(1)) * patch_width +
            (i  - ifirstP(0) + ghost_cells(0));
          int idx2dIM1 =
            (j - ifirstP(1) + ghost_cells(1)) * patch_width +
            (i - 1 - ifirstP(0) + ghost_cells(0));
          pressure_data[idx2d]    = pressure_data[idx2dIM1];
          temperature_data[idx2d] = temperature_data[idx2dIM1];
          source_data[idx2d]      = source_data[idx2dIM1];

          for (int k = 0; k < NEQU; k++) {
            int idx3d    = k * patch_area + idx2d;
            int idx3dIM1 = k * patch_area + idx2dIM1;
            w_n_data[idx3d] = w_n_data[idx3dIM1];

            if ((k == 1) && (i == ifirst(0)) &&
                (p_bc_type_R_int == REFLECT)) w_n_data[idx3d] = -w_n_data[idx3d];
          }
        }
      }
    }

    // Apply reflection to bottom boundary.
    if (bdry_loc == 2) {
      for (int j = ilast(1); j >= ifirst(1); j--) {
        for (int i = ilast(0); i >= ifirst(0); i--) {
          int idx2d =
            (j  - ifirstP(1) + ghost_cells(1)) * patch_width +
            (i - ifirstP(0) + ghost_cells(0));
          int idx2dJP1 =
            (j + 1 - ifirstP(1) + ghost_cells(1)) * patch_width +
            (i - ifirstP(0) + ghost_cells(0));
          pressure_data[idx2d]    = pressure_data[idx2dJP1];
          temperature_data[idx2d] = temperature_data[idx2dJP1];
          source_data[idx2d]      = source_data[idx2dJP1];

          for (int k = 0; k < NEQU; k++) {
            int idx3d    = k * patch_area + idx2d;
            int idx3dJP1 = k * patch_area + idx2dJP1;
            w_n_data[idx3d] = w_n_data[idx3dJP1];

            if ((k == 2) && (j == ilast(1)) &&
                (p_bc_type_B_int == REFLECT)) w_n_data[idx3d] = -w_n_data[idx3d];
          }
        }
      }
    }

    // Apply reflection to top boundary.
    if (bdry_loc == 3) {
      for (int j = ifirst(1); j <= ilast(1); j++) {
        for (int i = ifirst(0); i <= ilast(0); i++) {
          int idx2d =
            (j  - ifirstP(1) + ghost_cells(1)) * patch_width +
            (i - ifirstP(0) + ghost_cells(0));
          int idx2dJM1 =
            (j - 1 - ifirstP(1) + ghost_cells(1)) * patch_width +
            (i - ifirstP(0) + ghost_cells(0));
          pressure_data[idx2d]    = pressure_data[idx2dJM1];
          temperature_data[idx2d] = temperature_data[idx2dJM1];
          source_data[idx2d]      = source_data[idx2dJM1];

          for (int k = 0; k < NEQU; k++) {
            int idx3d    = k * patch_area + idx2d;
            int idx3dJM1 = k * patch_area + idx2dJM1;
            w_n_data[idx3d] = w_n_data[idx3dJM1];

            if ((k == 2) && (j == ifirst(1)) &&
                (p_bc_type_T_int == REFLECT)) w_n_data[idx3d] = -w_n_data[idx3d];
          }
        }
      }
    }
  }

  // I am not certain at the moment why I needed this section.
  // I will need to investigate it to give a more appropriate comment.
  const std::vector<hier::BoundaryBox> node_bdry =
    patch_geom->getNodeBoundaries();

  for (int r = 0; r < node_bdry.size(); r++) {
    int bdry_loc = node_bdry[r].getLocationIndex();
    hier::Box interior(patch.getBox());
    hier::Box fill_box =
      patch_geom->getBoundaryFillBox(node_bdry[r], interior, ghost_width_to_fill);
    const hier::Index ifirst = fill_box.lower();
    const hier::Index ilast  = fill_box.upper();

    if ((bdry_loc == 0) || (bdry_loc == 2)) {
      for (int j = ifirst(1); j <= ilast(1); j++) {
        for (int i = ilast(0); i >= ifirst(0); i--) {
          int idx2d =
            (j - ifirstP(1) + ghost_cells(1)) * patch_width +
            (i  - ifirstP(0) + ghost_cells(0));
          int idx2dIP1 =
            (j - ifirstP(1) + ghost_cells(1)) * patch_width +
            (i + 1 - ifirstP(0) + ghost_cells(0));
          pressure_data[idx2d]    = pressure_data[idx2dIP1];
          temperature_data[idx2d] = temperature_data[idx2dIP1];
          source_data[idx2d]      = source_data[idx2dIP1];

          for (int k = 0; k < NEQU; k++) {
            int idx3d    = k * patch_area + idx2d;
            int idx3dIP1 = k * patch_area + idx2dIP1;
            w_n_data[idx3d] = w_n_data[idx3dIP1];

            if ((k == 1) && (i == ilast(0)) &&
                (p_bc_type_L_int == REFLECT)) w_n_data[idx3d] = -w_n_data[idx3d];
          }
        }
      }
    }

    if ((bdry_loc == 1) || (bdry_loc == 3)) {
      for (int j = ifirst(1); j <= ilast(1); j++) {
        for (int i = ifirst(0); i <= ilast(0); i++) {
          int idx2d =
            (j - ifirstP(1) + ghost_cells(1)) * patch_width +
            (i  - ifirstP(0) + ghost_cells(0));
          int idx2dIM1 =
            (j - ifirstP(1) + ghost_cells(1)) * patch_width +
            (i - 1 - ifirstP(0) + ghost_cells(0));
          pressure_data[idx2d]    = pressure_data[idx2dIM1];
          temperature_data[idx2d] = temperature_data[idx2dIM1];
          source_data[idx2d]      = source_data[idx2dIM1];

          for (int k = 0; k < NEQU; k++) {
            int idx3d    = k * patch_area + idx2d;
            int idx3dIM1 = k * patch_area + idx2dIM1;
            w_n_data[idx3d] = w_n_data[idx3dIM1];

            if ((k == 1) && (i == ifirst(0)) &&
                (p_bc_type_R_int == REFLECT)) w_n_data[idx3d] = -w_n_data[idx3d];
          }
        }
      }
    }
  }

  // }

  t_bc->stop();
}

// Store visit writer handle in a private variable.
void enlightning::registerVisItDataWriter(
  boost::shared_ptr<appu::VisItDataWriter>viz_writer)
{
  p_visit_writer = viz_writer;
}

// Print information about all variables and constants, etc. used in the
// simulation.
void enlightning::printClassData(ostream& os) const
{
  fflush(stdout);

  os << "ptr enlightning = " << (enlightning *)this << endl;

  os << "ptr grid geometry = " << p_grid_geometry.get() << endl;

  os << "cfl: " << p_cfl << endl;
  os << "specific_dt: " << p_specific_dt << endl;

  os << "tag_tolerance: " << p_tag_tolerance[0] << ", " << p_tag_tolerance[1] <<
    endl;

  os << "max_timesteps: " << p_max_timesteps << endl;
  os << "end_time: " << p_end_time << endl;
  os << "regrid_step: " << p_regrid_step << endl;
  os << "tag_buffer: " << p_tag_buffer << endl;
  os << "viz_pressure: " << p_viz_pressure << endl;
  os << "viz_source: " << p_viz_source << endl;
  os << "viz_temperature: " << p_viz_temperature << endl;
  os << "viz_w_n: " << p_viz_w_n << endl;

  os << endl;
  os << "Boundary:" << endl;

  os << "bc_type_L: " << p_bc_type_L << endl;
  os << "bc_type_R: " << p_bc_type_R << endl;
  os << "bc_type_B: " << p_bc_type_B << endl;
  os << "bc_type_T: " << p_bc_type_T << endl;
  os << "bc_absorb_width: " << p_bc_absorb_width << endl;

  os << endl;
  os << "Recording:" << endl;

  os << "record_audio: " << p_record_audio << endl;
  os << "num_mics: " << p_num_mics << endl;

  for (int k = 0; k < p_num_mics; k++) {
    os << "mic_pos_x[" << k << "]: " << p_mic_pos_x[k] << endl;
    os << "mic_pos_y[" << k << "]: " << p_mic_pos_y[k] << endl;
  }
  os << "use_sample_rate_dt: " << p_use_sample_rate_dt << endl;
  os << "wav_sample_rate: " << p_wav_sample_rate << endl;

  os << endl;
  os << "Source:" << endl;

  os << "src_type: " << p_src_type << endl;
  os << "src_mode: " << p_src_mode << endl;
  os << "src_omega: " << p_src_omega << endl;
  os << "src_amp: " << p_src_amp << endl;
  os << "src_spread: " << p_src_spread << endl;
  os << "src_pos: " << p_src_pos[0] << ", " << p_src_pos[1] << endl;
  os << "src_sine_amp: " << p_src_sine_amp << endl;
  os << "src_sine_freq: " << p_src_sine_freq << endl;
  os << "src_num_pulses: " << p_src_num_pulses << endl;

  for (int k = 0; k < p_src_num_pulses; k++) {
    os << "src_pulse_times[" << k << "]: " << p_src_pulse_times[k] << endl;
  }

  os << endl;
  os << "Constants:" << endl;

  os << "weno_alpha: " << p_weno_alpha << endl;
  os << "c: " << p_c << endl;
  os << "gamma1: " << p_gamma1 << endl;
  os << "p0: " << p_p0 << endl;
  os << "RH: " << p_RH << endl;
  os << "T0: " << p_T0 << endl;
  os << "c_lr: " << p_c_lr << endl;
  os << "rho0: " << p_rho0 << endl;
  os << "R: " << p_R << endl;
  os << "c_p: " << p_c_p << endl;
  os << "T_star_n: " << p_T_star_n << endl;
  os << "T_star_o: " << p_T_star_o << endl;
  os << "mu: " << p_mu << endl;
  os << "muB: " << p_muB << endl;
  os << "kappa: " << p_kappa << endl;
  os << "c_v: " << p_c_v << endl;
  os << "R_tilde: " << p_R_tilde << endl;
  os << "M: " << p_M << endl;
  os << endl << endl;
}

// Use SAMRAI to parse enlightning section of the input file.
// Also checks for correctness of input parameters.
void enlightning::getFromInput(boost::shared_ptr<tbox::Database>db)
{
  if (db->keyExists("cfl")) {
    p_cfl = db->getDouble("cfl");

    if ((p_cfl <= 0) || (p_cfl >= 1)) {
      TBOX_ERROR(p_object_name << ":  "
                               << "Input `cfl' bad value."
                               << endl);
    }
  } else {
    TBOX_ERROR(p_object_name << ":  "
                             << "Input `cfl' not found."
                             << endl);
  }

  if (db->keyExists("specific_dt")) {
    p_specific_dt = db->getDouble("specific_dt");

    if (p_specific_dt < 0) {
      TBOX_ERROR(p_object_name << ":  "
                               << "Input `specific_dt' bad value."
                               << endl);
    }
  } else {
    TBOX_ERROR(p_object_name << ":  "
                             << "Input `specific_dt' not found."
                             << endl);
  }

  if (db->keyExists("tag_tolerance")) {
    db->getDoubleArray("tag_tolerance",
                       p_tag_tolerance, 3);

    if ((p_tag_tolerance[0] <= 0) || (p_tag_tolerance[1] <= 0)) {
      TBOX_ERROR(p_object_name << ":  "
                               << "Input `tag_tolerance' bad value."
                               << endl);
    }
  } else {
    TBOX_ERROR(p_object_name << ":  "
                             << "Input `tag_tolerance' not found."
                             << endl);
  }

  if (db->keyExists("max_timesteps")) {
    p_max_timesteps = db->getInteger("max_timesteps");

    if (p_max_timesteps <= 0) {
      TBOX_ERROR(p_object_name << ":  "
                               << "Input `max_timesteps' bad value."
                               << endl);
    }
  } else {
    TBOX_ERROR(p_object_name << ":  "
                             << "Input `max_timesteps' not found."
                             << endl);
  }

  p_start_time = 0.0;

  if (db->keyExists("end_time")) {
    p_end_time = db->getDouble("end_time");

    if (p_end_time <= 0) {
      TBOX_ERROR(p_object_name << ":  "
                               << "Input `end_time' bad value."
                               << endl);
    }
  } else {
    TBOX_ERROR(p_object_name << ":  "
                             << "Input `end_time' not found."
                             << endl);
  }

  if (db->keyExists("regrid_step")) {
    p_regrid_step = db->getInteger("regrid_step");

    if (p_regrid_step <= 0) {
      TBOX_ERROR(p_object_name << ":  "
                               << "Input `regrid_step' bad value."
                               << endl);
    }
  } else {
    TBOX_ERROR(p_object_name << ":  "
                             << "Input `regrid_step' not found."
                             << endl);
  }

  if (db->keyExists("tag_buffer")) {
    p_tag_buffer = db->getInteger("tag_buffer");

    if (p_tag_buffer <= 0) {
      TBOX_ERROR(p_object_name << ":  "
                               << "Input `tag_buffer' bad value."
                               << endl);
    }
  } else {
    TBOX_ERROR(p_object_name << ":  "
                             << "Input `tag_buffer' not found."
                             << endl);
  }

  if (db->keyExists("viz_pressure")) {
    p_viz_pressure = db->getInteger("viz_pressure");

    if ((p_viz_pressure != 0) && (p_viz_pressure != 1)) {
      TBOX_ERROR(p_object_name << ":  "
                               << "Input `viz_pressure' bad value."
                               << endl);
    }
  } else {
    TBOX_ERROR(p_object_name << ":  "
                             << "Input `viz_pressure' not found."
                             << endl);
  }

  if (db->keyExists("viz_source")) {
    p_viz_source = db->getInteger("viz_source");

    if ((p_viz_source != 0) && (p_viz_source != 1)) {
      TBOX_ERROR(p_object_name << ":  "
                               << "Input `viz_source' bad value."
                               << endl);
    }
  } else {
    TBOX_ERROR(p_object_name << ":  "
                             << "Input `viz_source' not found."
                             << endl);
  }

  if (db->keyExists("viz_temperature")) {
    p_viz_temperature = db->getInteger("viz_temperature");

    if ((p_viz_temperature != 0) && (p_viz_temperature != 1)) {
      TBOX_ERROR(p_object_name << ":  "
                               << "Input `viz_temperature' bad value."
                               << endl);
    }
  } else {
    TBOX_ERROR(p_object_name << ":  "
                             << "Input `viz_temperature' not found."
                             << endl);
  }

  if (db->keyExists("viz_w_n")) {
    p_viz_w_n = db->getInteger("viz_w_n");

    if ((p_viz_w_n != 0) && (p_viz_w_n != 1)) {
      TBOX_ERROR(p_object_name << ":  "
                               << "Input `viz_w_n' bad value."
                               << endl);
    }
  } else {
    TBOX_ERROR(p_object_name << ":  "
                             << "Input `viz_w_n' not found."
                             << endl);
  }

  boost::shared_ptr<tbox::Database> boundary_db;

  if (db->keyExists("boundary")) {
    boundary_db = db->getDatabase("boundary");
  } else {
    TBOX_ERROR(p_object_name << ": "
                             << "Input `boundary' database not found."
                             << endl);
  }

  if (boundary_db->keyExists("bc_type_L")) {
    p_bc_type_L = boundary_db->getString("bc_type_L");

    if (p_bc_type_L == "REFLECT") {
      p_bc_type_L_int = REFLECT;
    } else if (p_bc_type_L == "ABSORB") {
      p_bc_type_L_int = ABSORB;
    } else {
      TBOX_ERROR(p_object_name << ": "
                               << "Unknown bc_type_L string = "
                               << p_bc_type_L
                               << endl);
    }
  } else {
    TBOX_ERROR(p_object_name << ": "
                             << "Input `bc_type_L' value not found."
                             << endl);
  }

  if (boundary_db->keyExists("bc_type_R")) {
    p_bc_type_R = boundary_db->getString("bc_type_R");

    if (p_bc_type_R == "REFLECT") {
      p_bc_type_R_int = REFLECT;
    } else if (p_bc_type_R == "ABSORB") {
      p_bc_type_R_int = ABSORB;
    } else {
      TBOX_ERROR(p_object_name << ": "
                               << "Unknown bc_type_R string = "
                               << p_bc_type_R
                               << endl);
    }
  } else {
    TBOX_ERROR(p_object_name << ": "
                             << "Input `bc_type_R' value not found."
                             << endl);
  }

  if (boundary_db->keyExists("bc_type_B")) {
    p_bc_type_B = boundary_db->getString("bc_type_B");

    if (p_bc_type_B == "REFLECT") {
      p_bc_type_B_int = REFLECT;
    } else if (p_bc_type_B == "ABSORB") {
      p_bc_type_B_int = ABSORB;
    } else {
      TBOX_ERROR(p_object_name << ": "
                               << "Unknown bc_type_B string = "
                               << p_bc_type_B
                               << endl);
    }
  } else {
    TBOX_ERROR(p_object_name << ": "
                             << "Input `bc_type_B' value not found."
                             << endl);
  }

  if (boundary_db->keyExists("bc_type_T")) {
    p_bc_type_T = boundary_db->getString("bc_type_T");

    if (p_bc_type_T == "REFLECT") {
      p_bc_type_T_int = REFLECT;
    } else if (p_bc_type_T == "ABSORB") {
      p_bc_type_T_int = ABSORB;
    } else {
      TBOX_ERROR(p_object_name << ": "
                               << "Unknown bc_type_T string = "
                               << p_bc_type_T
                               << endl);
    }
  } else {
    TBOX_ERROR(p_object_name << ": "
                             << "Input `bc_type_T' value not found."
                             << endl);
  }

  if (boundary_db->keyExists("bc_absorb_width")) {
    p_bc_absorb_width = boundary_db->getDouble("bc_absorb_width");

    if (p_bc_absorb_width < 0) {
      TBOX_ERROR(p_object_name << ":  "
                               << "Input `bc_absorb_width' bad value."
                               << endl);
    }
  } else {
    TBOX_ERROR(p_object_name << ":  "
                             << "Input `bc_absorb_width' not found."
                             << endl);
  }

  boost::shared_ptr<tbox::Database> recording_db;

  if (db->keyExists("recording")) {
    recording_db = db->getDatabase("recording");
  } else {
    TBOX_ERROR(p_object_name << ": "
                             << "Input `recording' database not found."
                             << endl);
  }

  if (recording_db->keyExists("record_audio")) {
    p_record_audio = recording_db->getInteger("record_audio");

    if ((p_record_audio != 0) && (p_record_audio != 1)) {
      TBOX_ERROR(p_object_name << ":  "
                               << "Input `record_audio' bad value."
                               << endl);
    }
  } else {
    TBOX_ERROR(p_object_name << ":  "
                             << "Input `record_audio' not found."
                             << endl);
  }

  if (recording_db->keyExists("num_mics")) {
    p_num_mics = recording_db->getInteger("num_mics");

    if (p_num_mics < 1) {
      TBOX_ERROR(p_object_name << ":  "
                               << "Input `num_mics' bad value."
                               << endl);
    }
  } else {
    TBOX_ERROR(p_object_name << ":  "
                             << "Input `num_mics' not found."
                             << endl);
  }

  if (p_record_audio == 1) {
    if (recording_db->keyExists("mic_pos_x")) {
      p_mic_pos_x = recording_db->getDoubleVector("mic_pos_x");

      for (int i = 0; i < p_mic_pos_x.size(); i++) {
        if ((p_mic_pos_x[i] < 0) /*|| (p_mic_pos_x[i] > hi)*/) {
          TBOX_ERROR(p_object_name << ":  "
                                   << "Input `mic_pos_x['"
                                   << i << "] bad value."
                                   << endl);
        }
      }
    } else {
      TBOX_ERROR(p_object_name << ":  "
                               << "Input `mic_pos_x' not found."
                               << endl);
    }

    if (recording_db->keyExists("mic_pos_y")) {
      p_mic_pos_y = recording_db->getDoubleVector("mic_pos_y");

      for (int i = 0; i < p_mic_pos_y.size(); i++) {
        if ((p_mic_pos_y[i] < 0) /*|| (p_mic_pos_y[i] > hi)*/) {
          TBOX_ERROR(p_object_name << ":  "
                                   << "Input `mic_pos_y['"
                                   << i << "] bad value."
                                   << endl);
        }
      }
    } else {
      TBOX_ERROR(p_object_name << ":  "
                               << "Input `mic_pos_y' not found."
                               << endl);
    }

    if ((p_mic_pos_x.size() != p_num_mics)
        || (p_mic_pos_y.size() != p_num_mics)) {
      TBOX_ERROR(p_object_name << ":  "
                               << "'num_mics' does not match 'mic_pos' arrays"
                               << endl);
    }

    if (recording_db->keyExists("use_sample_rate_dt")) {
      p_use_sample_rate_dt = recording_db->getInteger("use_sample_rate_dt");

      if ((p_use_sample_rate_dt != 0) && (p_use_sample_rate_dt != 1)) {
        TBOX_ERROR(p_object_name << ":  "
                                 << "Input `use_sample_rate_dt' bad value."
                                 << endl);
      }
    } else {
      TBOX_ERROR(p_object_name << ":  "
                               << "Input `use_sample_rate_dt' not found."
                               << endl);
    }

    if (recording_db->keyExists("wav_sample_rate")) {
      p_wav_sample_rate = recording_db->getInteger("wav_sample_rate");

      if (p_wav_sample_rate <= 0) {
        TBOX_ERROR(p_object_name << ":  "
                                 << "Input `wav_sample_rate' bad value."
                                 << endl);
      }
    } else {
      TBOX_ERROR(p_object_name << ":  "
                               << "Input `wav_sample_rate' not found."
                               << endl);
    }
  } else {
    p_num_mics = 0;
    p_mic_pos_x.resize(0);
    p_mic_pos_y.resize(0);
    p_wav_sample_rate    = 22050;
    p_use_sample_rate_dt = 0;
  }

  boost::shared_ptr<tbox::Database> src_data_db;

  if (db->keyExists("src_data")) {
    src_data_db = db->getDatabase("src_data");
  } else {
    TBOX_ERROR(p_object_name << ": "
                             << "Input `src_data' database not found."
                             << endl);
  }

  if (src_data_db->keyExists("src_type")) {
    p_src_type = src_data_db->getString("src_type");

    if (p_src_type == "POINT") {
      p_src_type_int = POINT;
    } else if (p_src_type == "SINE") {
      p_src_type_int = SINE;
    } else if (p_src_type == "LIGHTNING") {
      p_src_type_int = LIGHTNING;
    } else {
      TBOX_ERROR(p_object_name << ": "
                               << "Unknown src_type string = "
                               << p_src_type
                               << endl);
    }
  } else {
    TBOX_ERROR(p_object_name << ": "
                             << "Input `src_type' value not found."
                             << endl);
  }

  if (src_data_db->keyExists("src_mode")) {
    p_src_mode = src_data_db->getString("src_mode");

    if (p_src_mode == "PULSE") {
      p_src_mode_int = PULSE;
    } else if (p_src_mode == "DRIVEN") {
      p_src_mode_int = DRIVEN;
    } else {
      TBOX_ERROR(p_object_name << ": "
                               << "Unknown src_mode string = "
                               << p_src_mode
                               << endl);
    }
  } else {
    TBOX_ERROR(p_object_name << ": "
                             << "Input `src_mode' value not found."
                             << endl);
  }

  if (src_data_db->keyExists("src_omega")) {
    p_src_omega = src_data_db->getDouble("src_omega");

    if (p_src_omega < 0) {
      TBOX_ERROR(p_object_name << ":  "
                               << "Input `src_omega' bad value."
                               << endl);
    }
  } else {
    TBOX_ERROR(p_object_name << ": "
                             << "Input `src_omega' not found."
                             << endl);
  }

  if (src_data_db->keyExists("src_amp")) {
    p_src_amp = src_data_db->getDouble("src_amp");

    if (p_src_amp < 0) {
      TBOX_ERROR(p_object_name << ":  "
                               << "Input `src_amp' bad value."
                               << endl);
    }
  } else {
    TBOX_ERROR(p_object_name << ": "
                             << "Input `src_amp' not found."
                             << endl);
  }

  if (src_data_db->keyExists("src_spread")) {
    p_src_spread = src_data_db->getDouble("src_spread");

    if (p_src_spread < 0) {
      TBOX_ERROR(p_object_name << ":  "
                               << "Input `src_spread' bad value."
                               << endl);
    }
  } else {
    TBOX_ERROR(p_object_name << ": "
                             << "Input `src_spread' not found."
                             << endl);
  }

  if (src_data_db->keyExists("src_pos")) {
    src_data_db->getDoubleArray("src_pos", p_src_pos, p_dim.getValue());

    if ((p_src_pos[0] < 0) || (p_src_pos[1] < 0) /*|| (p_src_pos[0] >= hi) 
						   || (p_src_pos[1] >= hi)*/) {
      TBOX_ERROR(p_object_name << ":  "
                               << "Input `src_pos' bad value."
                               << endl);
    }
  } else {
    TBOX_ERROR(p_object_name << ": "
                             << "Input `src_pos' not found."
                             << endl);
  }

  if (src_data_db->keyExists("src_sine_amp")) {
    p_src_sine_amp = src_data_db->getDouble("src_sine_amp");

    if (p_src_sine_amp < 0) {
      TBOX_ERROR(p_object_name << ":  "
                               << "Input `src_sine_amp' bad value."
                               << endl);
    }
  } else {
    p_src_sine_amp = 1.0;
  }

  if (src_data_db->keyExists("src_sine_freq")) {
    p_src_sine_freq = src_data_db->getDouble("src_sine_freq");

    if (p_src_sine_freq < 0) {
      TBOX_ERROR(p_object_name << ":  "
                               << "Input `src_sine_freq' bad value."
                               << endl);
    }
  } else {
    p_src_sine_freq = 1.0;
  }

  if (src_data_db->keyExists("src_num_pulses")) {
    p_src_num_pulses = src_data_db->getInteger("src_num_pulses");

    if (p_src_num_pulses < 1) {
      TBOX_ERROR(p_object_name << ":  "
                               << "Input `src_num_pulses' bad value."
                               << endl);
    }
  } else {
    TBOX_ERROR(p_object_name << ":  "
                             << "Input `src_num_pulses' not found."
                             << endl);
  }

  if (src_data_db->keyExists("src_pulse_times")) {
    p_src_pulse_times = src_data_db->getIntegerVector("src_pulse_times");

    for (int i = 0; i < p_src_pulse_times.size(); i++) {
      if (p_src_pulse_times[i] < 0) {
        TBOX_ERROR(p_object_name << ":  "
                                 << "Input `src_pulse_times['"
                                 << i << "] bad value."
                                 << endl);
      }
    }
  } else {
    TBOX_ERROR(p_object_name << ":  "
                             << "Input `src_pulse_times' not found."
                             << endl);
  }

  if (p_src_pulse_times.size() != p_src_num_pulses) {
    TBOX_ERROR(p_object_name << ":  "
                             << "'src_num_pulses' does not match 'src_pulse_times' array"
                             << endl);
  }

  boost::shared_ptr<tbox::Database> constants_db;

  if (db->keyExists("constants")) {
    constants_db = db->getDatabase("constants");
  } else {
    TBOX_ERROR(p_object_name << ": "
                             << "Input `constants' database not found."
                             << endl);
  }

  if (constants_db->keyExists("weno_alpha")) {
    p_weno_alpha = constants_db->getDouble("weno_alpha");

    if ((p_weno_alpha < 300) || (p_weno_alpha > 400)) {
      TBOX_ERROR(p_object_name << ":  "
                               << "Input `weno_alpha' bad value."
                               << endl);
    }
  } else {
    TBOX_ERROR(p_object_name << ": "
                             << "Input `weno_alpha' not found."
                             << endl);
  }

  if (constants_db->keyExists("c")) {
    p_c = constants_db->getDouble("c");

    if (p_c <= 0) {
      TBOX_ERROR(p_object_name << ":  "
                               << "Input `c' bad value."
                               << endl);
    }
  } else {
    TBOX_ERROR(p_object_name << ": "
                             << "Input `c' not found."
                             << endl);
  }

  if (constants_db->keyExists("gamma1")) {
    p_gamma1 = constants_db->getDouble("gamma1");

    if (p_gamma1 <= 0) {
      TBOX_ERROR(p_object_name << ":  "
                               << "Input `gamma1' bad value."
                               << endl);
    }
  } else {
    TBOX_ERROR(p_object_name << ": "
                             << "Input `gamma1' not found."
                             << endl);
  }

  if (constants_db->keyExists("p0")) {
    p_p0 = constants_db->getDouble("p0");

    if (p_p0 <= 0) {
      TBOX_ERROR(p_object_name << ":  "
                               << "Input `p0' bad value."
                               << endl);
    }
  } else {
    TBOX_ERROR(p_object_name << ": "
                             << "Input `p0' not found."
                             << endl);
  }

  if (constants_db->keyExists("RH")) {
    p_RH = constants_db->getDouble("RH");

    if ((p_RH < 0) || (p_RH > 100)) {
      TBOX_ERROR(p_object_name << ":  "
                               << "Input `RH' bad value."
                               << endl);
    }
  } else {
    TBOX_ERROR(p_object_name << ": "
                             << "Input `RH' not found."
                             << endl);
  }

  if (constants_db->keyExists("T0")) {
    p_T0 = constants_db->getDouble("T0");

    if (p_T0 <= 0) {
      TBOX_ERROR(p_object_name << ":  "
                               << "Input `T0' bad value."
                               << endl);
    }
  } else {
    TBOX_ERROR(p_object_name << ": "
                             << "Input `T0' not found."
                             << endl);
  }

  if (constants_db->keyExists("c_lr")) {
    p_c_lr = constants_db->getDouble("c_lr");
  } else {
    TBOX_ERROR(p_object_name << ": "
                             << "Input `c_lr' not found."
                             << endl);
  }

  if (constants_db->keyExists("rho0")) {
    p_rho0 = constants_db->getDouble("rho0");

    if (p_rho0 <= 0) {
      TBOX_ERROR(p_object_name << ":  "
                               << "Input `rho0' bad value."
                               << endl);
    }
  } else {
    TBOX_ERROR(p_object_name << ": "
                             << "Input `rho0' not found."
                             << endl);
  }

  if (constants_db->keyExists("R")) {
    p_R = constants_db->getDouble("R");

    if (p_R <= 0) {
      TBOX_ERROR(p_object_name << ":  "
                               << "Input `R' bad value."
                               << endl);
    }
  } else {
    TBOX_ERROR(p_object_name << ": "
                             << "Input `R' not found."
                             << endl);
  }

  if (constants_db->keyExists("c_p")) {
    p_c_p = constants_db->getDouble("c_p");
  } else {
    TBOX_ERROR(p_object_name << ": "
                             << "Input `c_p' not found."
                             << endl);
  }

  if (constants_db->keyExists("T_star_n")) {
    p_T_star_n = constants_db->getDouble("T_star_n");
  } else {
    TBOX_ERROR(p_object_name << ": "
                             << "Input `T_star_n' not found."
                             << endl);
  }

  if (constants_db->keyExists("T_star_o")) {
    p_T_star_o = constants_db->getDouble("T_star_o");
  } else {
    TBOX_ERROR(p_object_name << ": "
                             << "Input `T_star_o' not found."
                             << endl);
  }

  if (constants_db->keyExists("mu")) {
    p_mu = constants_db->getDouble("mu");
  } else {
    TBOX_ERROR(p_object_name << ": "
                             << "Input `mu' not found."
                             << endl);
  }

  if (constants_db->keyExists("muB")) {
    p_muB = constants_db->getDouble("muB");
  } else {
    TBOX_ERROR(p_object_name << ": "
                             << "Input `muB' not found."
                             << endl);
  }

  if (constants_db->keyExists("kappa")) {
    p_kappa = constants_db->getDouble("kappa");
  } else {
    TBOX_ERROR(p_object_name << ": "
                             << "Input `kappa' not found."
                             << endl);
  }

  if (constants_db->keyExists("c_v")) {
    p_c_v = constants_db->getDouble("c_v");
  } else {
    TBOX_ERROR(p_object_name << ": "
                             << "Input `c_v' not found."
                             << endl);
  }

  if (constants_db->keyExists("R_tilde")) {
    p_R_tilde = constants_db->getDouble("R_tilde");
  } else {
    TBOX_ERROR(p_object_name << ": "
                             << "Input `R_tilde' not found."
                             << endl);
  }

  if (constants_db->keyExists("M")) {
    p_M = constants_db->getDouble("M");
  } else {
    TBOX_ERROR(p_object_name << ": "
                             << "Input `M' not found."
                             << endl);
  }
}

// Save all variable and constants, etc. to database for checkpointing
// simulation.
void enlightning::putToRestart(const boost::shared_ptr<tbox::Database>& db) const
{
  db->putIntegerArray("p_nghosts",     &p_nghosts[0],     p_dim.getValue());
  db->putIntegerArray("p_zero_ghosts", &p_zero_ghosts[0], p_dim.getValue());

  db->putInteger("p_max_timesteps", p_max_timesteps);
  db->putDouble("p_start_time", p_start_time);
  db->putDouble("p_end_time",   p_end_time);
  db->putInteger("p_regrid_step", p_regrid_step);
  db->putInteger("p_tag_buffer",  p_tag_buffer);
  db->putDouble("p_loop_time", p_loop_time);
  db->putInteger("p_iteration_number", p_iteration_number);

  db->putDouble("p_cfl",         p_cfl);
  db->putDouble("p_specific_dt", p_specific_dt);
  db->putDoubleArray("p_tag_tolerance", p_tag_tolerance, 3);
  db->putInteger("p_viz_pressure",    p_viz_pressure);
  db->putInteger("p_viz_source",      p_viz_source);
  db->putInteger("p_viz_temperature", p_viz_temperature);
  db->putInteger("p_viz_w_n",         p_viz_w_n);

  db->putString("p_bc_type_L", p_bc_type_L);
  db->putInteger("p_bc_type_L_int", p_bc_type_L_int);
  db->putString("p_bc_type_R", p_bc_type_R);
  db->putInteger("p_bc_type_R_int", p_bc_type_R_int);
  db->putString("p_bc_type_B", p_bc_type_B);
  db->putInteger("p_bc_type_B_int", p_bc_type_B_int);
  db->putString("p_bc_type_T", p_bc_type_T);
  db->putInteger("p_bc_type_T_int", p_bc_type_T_int);
  db->putDouble("p_bc_absorb_width", p_bc_absorb_width);

  db->putInteger("p_record_audio", p_record_audio);
  db->putInteger("p_num_mics",     p_num_mics);

  if (p_num_mics > 0) {
    db->putDoubleVector("p_mic_pos_x", p_mic_pos_x);
    db->putDoubleVector("p_mic_pos_y", p_mic_pos_y);
  }
  db->putInteger("p_use_sample_rate_dt", p_use_sample_rate_dt);
  db->putInteger("p_wav_sample_rate",    p_wav_sample_rate);

  db->putString("p_src_type", p_src_type);
  db->putInteger("p_src_type_int", p_src_type_int);
  db->putString("p_src_mode", p_src_mode);
  db->putInteger("p_src_mode_int", p_src_mode_int);
  db->putDouble("p_src_omega",  p_src_omega);
  db->putDouble("p_src_amp",    p_src_amp);
  db->putDouble("p_src_spread", p_src_spread);
  db->putDoubleArray("p_src_pos", p_src_pos, 2);
  db->putDouble("p_src_sine_amp",  p_src_sine_amp);
  db->putDouble("p_src_sine_freq", p_src_sine_freq);
  db->putInteger("p_src_num_pulses", p_src_num_pulses);
  db->putIntegerVector("p_src_pulse_times", p_src_pulse_times);

  db->putDouble("p_p0",         p_p0);
  db->putDouble("p_sampled_p0", p_sampled_p0);
  db->putDouble("p_T0",         p_T0);
  db->putDouble("p_c_lr",       p_c_lr);
  db->putDouble("p_rho0",       p_rho0);
  db->putDouble("p_c",          p_c);
  db->putDouble("p_gamma1",     p_gamma1);
  db->putDouble("p_RH",         p_RH);
  db->putDouble("p_R",          p_R);
  db->putDouble("p_c_p",        p_c_p);
  db->putDouble("p_T_star_n",   p_T_star_n);
  db->putDouble("p_T_star_o",   p_T_star_o);
  db->putDouble("p_mu",         p_mu);
  db->putDouble("p_muB",        p_muB);
  db->putDouble("p_kappa",      p_kappa);
  db->putDouble("p_weno_alpha", p_weno_alpha);
  db->putDouble("p_tau_n",      p_tau_n);
  db->putDouble("p_tau_o",      p_tau_o);
  db->putDouble("p_R_tilde",    p_R_tilde);
  db->putDouble("p_c_v",        p_c_v);
  db->putDouble("p_M",          p_M);

  if (p_src_type_int == LIGHTNING) {
    db->putInteger("p_lightning_input_line_count",
                   p_lightning_input_line_count);
    db->putDoubleVector("p_lightning_coords",
                        p_lightning_coords);
  }
}

// Retrieve all variables and constants, etc. from restart database for
// restarting a simulation.
void enlightning::getFromRestart()
{
  boost::shared_ptr<tbox::Database> root_db =
    tbox::RestartManager::getManager()->getRootDatabase();

  boost::shared_ptr<tbox::Database> restart_db;

  if (root_db->isDatabase(p_object_name)) {
    restart_db = root_db->getDatabase(p_object_name);
  } else {
    TBOX_ERROR("Restart database corresponding to "
               << p_object_name << " not found in the restart file.");
  }

  int *tmp_nghosts = &p_nghosts[0];
  restart_db->getIntegerArray("p_nghosts", tmp_nghosts, p_dim.getValue());

  for (int i = 0; i < p_dim.getValue(); i++) {
    if (p_nghosts(i) != CELLG) {
      TBOX_ERROR(p_object_name << ": "
                               << "Key data `p_nghosts' in restart file != CELLG." <<
                 endl);
    }
  }
  int *tmp_zero_ghosts = &p_zero_ghosts[0];
  restart_db->getIntegerArray("p_zero_ghosts", tmp_zero_ghosts, p_dim.getValue());

  for (int i = 0; i < p_dim.getValue(); i++) {
    if (p_zero_ghosts(i) != 0) {
      TBOX_ERROR(p_object_name << ": "
                               << "Key data `p_zero_ghosts' in restart file != 0." <<
                 endl);
    }
  }

  p_max_timesteps    = restart_db->getInteger("p_max_timesteps");
  p_start_time       = restart_db->getDouble("p_start_time");
  p_end_time         = restart_db->getDouble("p_end_time");
  p_regrid_step      = restart_db->getInteger("p_regrid_step");
  p_tag_buffer       = restart_db->getInteger("p_tag_buffer");
  p_loop_time        = restart_db->getDouble("p_loop_time");
  p_iteration_number = restart_db->getInteger("p_iteration_number");

  p_cfl              = restart_db->getDouble("p_cfl");
  p_specific_dt      = restart_db->getDouble("p_specific_dt");
  restart_db->getDoubleArray("p_tag_tolerance", p_tag_tolerance, 3);
  p_viz_pressure     = restart_db->getInteger("p_viz_pressure");
  p_viz_source       = restart_db->getInteger("p_viz_source");
  p_viz_temperature  = restart_db->getInteger("p_viz_temperature");
  p_viz_w_n          = restart_db->getInteger("p_viz_w_n");

  p_bc_type_L        = restart_db->getString("p_bc_type_L");
  p_bc_type_L_int    = restart_db->getInteger("p_bc_type_L_int");
  p_bc_type_R        = restart_db->getString("p_bc_type_R");
  p_bc_type_R_int    = restart_db->getInteger("p_bc_type_R_int");
  p_bc_type_B        = restart_db->getString("p_bc_type_B");
  p_bc_type_B_int    = restart_db->getInteger("p_bc_type_B_int");
  p_bc_type_T        = restart_db->getString("p_bc_type_T");
  p_bc_type_T_int    = restart_db->getInteger("p_bc_type_T_int");
  p_bc_absorb_width  = restart_db->getDouble("p_bc_absorb_width");

  p_record_audio = restart_db->getInteger("p_record_audio");
  p_num_mics     = restart_db->getInteger("p_num_mics");

  if (p_num_mics > 0) {
    p_mic_pos_x = restart_db->getDoubleVector("p_mic_pos_x");
    p_mic_pos_y = restart_db->getDoubleVector("p_mic_pos_y");
  } else {
    p_mic_pos_x.resize(0);
    p_mic_pos_y.resize(0);
  }
  p_use_sample_rate_dt = restart_db->getInteger("p_use_sample_rate_dt");
  p_wav_sample_rate    = restart_db->getInteger("p_wav_sample_rate");

  p_src_type     = restart_db->getString("p_src_type");
  p_src_type_int = restart_db->getInteger("p_src_type_int");
  p_src_mode     = restart_db->getString("p_src_mode");
  p_src_mode_int = restart_db->getInteger("p_src_mode_int");
  p_src_omega    = restart_db->getDouble("p_src_omega");
  p_src_amp      = restart_db->getDouble("p_src_amp");
  p_src_spread   = restart_db->getDouble("p_src_spread");
  restart_db->getDoubleArray("p_src_pos", p_src_pos, 2);
  p_src_sine_amp    = restart_db->getDouble("p_src_sine_amp");
  p_src_sine_freq   = restart_db->getDouble("p_src_sine_freq");
  p_src_num_pulses  = restart_db->getInteger("p_src_num_pulses");
  p_src_pulse_times = restart_db->getIntegerVector("p_src_pulse_times");

  p_p0         = restart_db->getDouble("p_p0");
  p_sampled_p0 = restart_db->getDouble("p_sampled_p0");
  p_T0         = restart_db->getDouble("p_T0");
  p_c_lr       = restart_db->getDouble("p_c_lr");
  p_rho0       = restart_db->getDouble("p_rho0");
  p_c          = restart_db->getDouble("p_c");
  p_gamma1     = restart_db->getDouble("p_gamma1");
  p_RH         = restart_db->getDouble("p_RH");
  p_R          = restart_db->getDouble("p_R");
  p_c_p        = restart_db->getDouble("p_c_p");
  p_T_star_n   = restart_db->getDouble("p_T_star_n");
  p_T_star_o   = restart_db->getDouble("p_T_star_o");
  p_mu         = restart_db->getDouble("p_mu");
  p_muB        = restart_db->getDouble("p_muB");
  p_kappa      = restart_db->getDouble("p_kappa");
  p_weno_alpha = restart_db->getDouble("p_weno_alpha");
  p_tau_n      = restart_db->getDouble("p_tau_n");
  p_tau_o      = restart_db->getDouble("p_tau_o");
  p_R_tilde    = restart_db->getDouble("p_R_tilde");
  p_c_v        = restart_db->getDouble("p_c_v");
  p_M          = restart_db->getDouble("p_M");

  if (p_src_type_int == LIGHTNING) {
    p_lightning_input_line_count =
      restart_db->getInteger("p_lightning_input_line_count");
    p_lightning_coords =
      restart_db->getDoubleVector("p_lightning_coords");
  }
}

int enlightning::getMaxTimesteps()
{
  return p_max_timesteps;
}

double enlightning::getStartTime()
{
  return p_start_time;
}

double enlightning::getEndTime()
{
  return p_end_time;
}

int enlightning::getRegridStep()
{
  return p_regrid_step;
}

int enlightning::getTagBuffer()
{
  return p_tag_buffer;
}

double enlightning::getLoopTime()
{
  return p_loop_time;
}

void enlightning::setLoopTime(const double loop_time)
{
  p_loop_time = loop_time;
}

int enlightning::getIterationNumber()
{
  return p_iteration_number;
}

int enlightning::shouldRecordAudio()
{
  return p_record_audio;
}

void enlightning::setIterationNumber(const int iter_num)
{
  p_iteration_number = iter_num;
}

int enlightning::getSourceMode()
{
  return p_src_mode_int;
}

// Return true if the source should be drawn, based on its type.
int enlightning::shouldAddSource()
{
  int add_source = 0;

  if ((p_src_mode_int == DRIVEN) &&
      (p_iteration_number >= p_src_pulse_times[0])) {
    add_source = 1;
  } else if (p_src_mode_int == PULSE) {
    for (int i = 0; i < p_src_num_pulses; i++) {
      // p_iteration_number+1 because this is called before main() increments
      // it.
      if ((p_iteration_number + 1) == p_src_pulse_times[i]) {
        add_source = 1;
      }
    }
  }

  return add_source;
}

// Draw source on patch (point, sine, or lightning channel pulse, or driven
// point source).
void enlightning::addSource(boost::shared_ptr<hier::PatchHierarchy>patch_hier,
                            double dt)
{
  t_addsource->start();
  const int num_levels = patch_hier->getNumberOfLevels();

  for (int q = 0; q < num_levels; q++) {
    boost::shared_ptr<hier::PatchLevel> level(patch_hier->getPatchLevel(q));

    for (hier::PatchLevel::iterator patch_iterator(level->begin());
         patch_iterator != level->end();
         ++patch_iterator) {
      boost::shared_ptr<hier::Patch> patch(
        level->getPatch(patch_iterator->getBox().getBoxId()));

      boost::shared_ptr<pdat::CellData<double> > source_grid(
        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
          patch->getPatchData(p_source, getInteriorContext())));

      const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
        BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
          patch->getPatchGeometry()));

      const double *dx       = patch_geom->getDx();
      const double *patchXLo = patch_geom->getXLower();

      const hier::IntVector ghost_cells = source_grid->getGhostCellWidth();

      const hier::Index ifirst = patch->getBox().lower();
      const hier::Index ilast  = patch->getBox().upper();

      double   *source_data = source_grid->getPointer();
      const int patch_width = (ilast(0) + 1) - ifirst(0) + 2 * ghost_cells(0);
      double    xc[2], x0, x1;
      int idx;

      if (p_src_mode_int == PULSE) {
        // Draw a point source pulse.
        if (p_src_type_int == POINT) {
          for (int j = ifirst(1); j <= ilast(1); j++) {
            xc[1] = patchXLo[1] + dx[1] * ((double)(j - ifirst(1)) + 0.5);
            x1    = xc[1] - p_src_pos[1];

            for (int i = ifirst(0); i <= ilast(0); i++) {
              xc[0] = patchXLo[0] + dx[0] * ((double)(i - ifirst(0)) + 0.5);
              x0    = xc[0] - p_src_pos[0];

              idx              = (j - ifirst(1)) * patch_width + (i - ifirst(0));
              source_data[idx] = p_src_amp * exp(-log(2.0) / SQUARE(p_src_spread) *
                                 (x0 * x0 + x1 * x1));
            }
          }
        }

        // Draw a sine wave shaped pulse (used before lightning channel
        // capability was added).
        if (p_src_type_int == SINE) {
          for (int j = ifirst(1); j <= ilast(1); j++) {
            xc[1] = patchXLo[1] + dx[1] * ((double)(j - ifirst(1)) + 0.5);
            x1    = xc[1] - p_src_pos[1];

            for (int i = ifirst(0); i <= ilast(0); i++) {
              xc[0] = patchXLo[0] + dx[0] * ((double)(i - ifirst(0)) + 0.5);
              x0    = xc[0] - p_src_pos[0];

              idx              = (j - ifirst(1)) * patch_width + (i - ifirst(0));
              source_data[idx] = p_src_amp * exp(-log(2.0) / SQUARE(p_src_spread) *
                                 SQUARE((x0) - p_src_sine_amp * sin(p_src_sine_freq * (x1))));
            }
          }
        }

        // Draw a lightning channel source using coordinates already read from
        // file.
        if (p_src_type_int == LIGHTNING) {
          for (int k = 0; k < p_lightning_input_line_count; k++) {
            int orientation = 0;
            double theta    = -atan((p_lightning_coords[k * 4 + 3]
		    - p_lightning_coords[k * 4 + 1])
                    / (p_lightning_coords[k * 4 + 2] -
                       p_lightning_coords[k * 4 + 0]));

            if (((p_lightning_coords[k * 4 + 0] > p_lightning_coords[k * 4 + 2])
                 && (p_lightning_coords[k * 4 + 1] >
                     p_lightning_coords[k * 4 + 3]))
                || ((p_lightning_coords[k * 4 + 0] <=
                     p_lightning_coords[k * 4 + 2])
                    && (p_lightning_coords[k * 4 + 1] >=
                        p_lightning_coords[k * 4 + 3]))) orientation = 1;
            else if (((p_lightning_coords[k * 4 + 0] >=
                       p_lightning_coords[k * 4 + 2])
                      && (p_lightning_coords[k * 4 + 1] <=
                          p_lightning_coords[k * 4 + 3]))
                     || ((p_lightning_coords[k * 4 + 0] <
                          p_lightning_coords[k * 4 + 2])
                         && (p_lightning_coords[k * 4 + 1] <
                             p_lightning_coords[k * 4 + 3]))) orientation = 2;

            for (int j = ifirst(1); j <= ilast(1); j++) {
              xc[1] = patchXLo[1] + dx[1] * ((double)(j - ifirst(1)) + 0.5);

              for (int i = ifirst(0); i <= ilast(0); i++) {
                xc[0] = patchXLo[0] + dx[0] * ((double)(i - ifirst(0)) + 0.5);
                double line0 = tan(PI / 2.0 - theta)
                               * (xc[0] - p_lightning_coords[k * 4 + 0]) +
                               p_lightning_coords[k * 4 + 1];
                double line1 = tan(PI / 2.0 - theta)
                               * (xc[0] - p_lightning_coords[k * 4 + 2]) +
                               p_lightning_coords[k * 4 + 3];
                idx = (j - ifirst(1)) * patch_width + (i - ifirst(0));

                if (((orientation == 1) && (xc[1] >= line0))
                    || ((orientation == 2) && (xc[1] <= line0))) {
                  double temp = p_src_amp * exp(-log(2.0) / SQUARE(p_src_spread)
                         * (SQUARE(xc[0] - p_lightning_coords[k * 4 + 0])
                         + SQUARE(xc[1] - p_lightning_coords[k * 4 + 1])));
                  source_data[idx] = source_data[idx] + temp -
                                     (source_data[idx] / p_src_amp) * temp;
                }

                if (((orientation == 1) && (xc[1] <= line1))
                    || ((orientation == 2) && (xc[1] >= line1))) {
                  double temp = p_src_amp * exp(-log(2.0) / SQUARE(p_src_spread)
                         * (SQUARE(xc[0] - p_lightning_coords[k * 4 + 2])
                         + SQUARE(xc[1] - p_lightning_coords[k * 4 + 3])));
                  source_data[idx] = source_data[idx] + temp -
                                     (source_data[idx] / p_src_amp) * temp;
                }

                if (((orientation == 1) && (xc[1] < line0) && (xc[1] > line1))
                    || ((orientation == 2) && (xc[1] > line0) &&
                        (xc[1] < line1))) {
                  double temp = p_src_amp * exp(-log(2.0) / SQUARE(p_src_spread)
                         * (SQUARE(sin(theta) * (xc[0] - p_lightning_coords[k * 4 + 0])
                         + cos( theta) * (xc[1] - p_lightning_coords[k * 4 + 1]))));
                  source_data[idx] = source_data[idx] + temp -
                                     (source_data[idx] / p_src_amp) * temp;
                }
              }
            }
          }
        }
      } else if (p_src_mode_int == DRIVEN) {
        // If a driven point source, need to oscillate it.
        if (p_src_type_int == POINT) {
          for (int j = ifirst(1); j <= ilast(1); j++) {
            xc[1] = patchXLo[1] + dx[1] * ((double)(j - ifirst(1)) + 0.5);
            x1    = xc[1] - p_src_pos[1];

            for (int i = ifirst(0); i <= ilast(0); i++) {
              xc[0] = patchXLo[0] + dx[0] * ((double)(i - ifirst(0)) + 0.5);
              x0    = xc[0] - p_src_pos[0];

              idx              = (j - ifirst(1)) * patch_width + (i - ifirst(0));
              source_data[idx] = p_src_amp * exp(-log(2.0) / SQUARE(p_src_spread) *
                                 (x0 * x0 + x1 * x1)) * sin(p_src_omega *
                                 (p_iteration_number + 1 - p_src_pulse_times[0]) * dt);
            }
          }
        }

        if ((p_src_type_int == SINE) || (p_src_type_int == LIGHTNING)) TBOX_ERROR(
            "DRIVEN SINE or LIGHTNING source mode not implemented" << endl);
      }
    }
  }
  t_addsource->stop();
}

// Used to zero out the entire source grid.
void enlightning::setSourceToZero(
  boost::shared_ptr<hier::PatchHierarchy>patch_hier)
{
  const int num_levels = patch_hier->getNumberOfLevels();

  for (int q = 0; q < num_levels; q++) {
    boost::shared_ptr<hier::PatchLevel> level(patch_hier->getPatchLevel(q));

    for (hier::PatchLevel::iterator patch_iterator(level->begin());
         patch_iterator != level->end();
         ++patch_iterator) {
      boost::shared_ptr<hier::Patch> mypatch(
        level->getPatch(patch_iterator->getBox().getBoxId()));
      boost::shared_ptr<pdat::CellData<double> > source_grid(
        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
          mypatch->getPatchData(p_source, getInteriorContext())));
      source_grid->fillAll(0);
    }
  }
}

// If source is a lightning channel, read in the geometry of the channel from
// file.
void enlightning::readLightningInput(const char *input_lightning_file)
{
  FILE *fr;
  float x1, y1, x2, y2;
  char  line[80];

  p_lightning_input_line_count = 0;

  fr = fopen(input_lightning_file, "rt");

  if (fr == NULL) {
    TBOX_ERROR(p_object_name << ": "
                             << "Can't open '" << input_lightning_file << "' for reading."
                             << endl);
  }

  // Count the number of line segments first.
  while (fgets(line, 80, fr) != NULL) {
    sscanf(line, "%f %f %f %f", &x1, &y1, &x2, &y2);
    p_lightning_input_line_count++;
  }

  if (p_lightning_input_line_count < 1) {
    TBOX_ERROR(p_object_name << ": "
                             << "File '" << input_lightning_file << "' is empty."
                             << endl);
  }

  p_lightning_coords.resize(p_lightning_input_line_count * 4);

  // Initialize array that will store segment coordinates in memory.
  for (int i = 0; i < p_lightning_input_line_count; i++) {
    p_lightning_coords[i * 4 + 0] = 0;
    p_lightning_coords[i * 4 + 1] = 0;
    p_lightning_coords[i * 4 + 2] = 0;
    p_lightning_coords[i * 4 + 3] = 0;
  }

  // Read segment coordinates in from file.
  rewind(fr);
  int j = 0;

  while (fgets(line, 80, fr) != NULL) {
    sscanf(line, "%f %f %f %f", &x1, &y1, &x2, &y2);
    p_lightning_coords[j * 4 + 0] = x1;
    p_lightning_coords[j * 4 + 1] = y1;
    p_lightning_coords[j * 4 + 2] = x2;
    p_lightning_coords[j * 4 + 3] = y2;

    // printf("%f %f %f %f\n", p_lightning_coords[j*4+0],
    // p_lightning_coords[j*4+1],p_lightning_coords[j*4+2],
    // p_lightning_coords[j*4+3]);
    j++;
  }

  fclose(fr);
}

// Record the pressure values at each timestep for each requested microphone.
void enlightning::recordAudio(boost::shared_ptr<hier::PatchHierarchy>patch_hier)
{
  t_record->start();
  const int num_levels = patch_hier->getNumberOfLevels();

  // Iterate over grid levels.
  for (int q = 0; q < num_levels; q++) {
    boost::shared_ptr<hier::PatchLevel> level(patch_hier->getPatchLevel(q));

    // Iterate over patches on grid level.
    for (hier::PatchLevel::iterator patch_iterator(level->begin());
         patch_iterator != level->end();
         ++patch_iterator) {
      boost::shared_ptr<hier::Patch> mypatch(
        level->getPatch(patch_iterator->getBox().getBoxId()));
      const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
        BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
          mypatch->getPatchGeometry()));
      const double *patchXLo = patch_geom->getXLower();
      const double *patchXHi = patch_geom->getXUpper();

      // Iterate over number of microphones.
      for (int k = 0; k < p_num_mics; k++) {
        // If the microphone location is in the patch, record its value.
        if ((p_mic_pos_x[k] >= patchXLo[0]) && (p_mic_pos_x[k] <= patchXHi[0])
            && (p_mic_pos_y[k] >= patchXLo[1]) &&
            (p_mic_pos_y[k] <= patchXHi[1])) {
          FILE *fw1, *fw2;
          char  filename1[64];
          char  filename2[64];

          // Create two file per microphone (one for info about the microphone,
          // the other is the actual recorded pressure data).
          sprintf(filename1, "enlightning.record/mic-%d-info.txt", k + 1);
          sprintf(filename2, "enlightning.record/mic-%d.txt",      k + 1);
          const double *dx         = patch_geom->getDx();
          const hier::Index ifirst = mypatch->getBox().lower();
          const hier::Index ilast  = mypatch->getBox().upper();

          // Get pointer to pressure data.
          boost::shared_ptr<pdat::CellData<double> > pressure_grid(
            BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
              mypatch->getPatchData(p_pressure, getInteriorContext())));
          double *p_data = pressure_grid->getPointer();
          int     i      = (int)((p_mic_pos_x[k] - patchXLo[0]) / dx[0]);
          int     j      = (int)((p_mic_pos_y[k] - patchXLo[1]) / dx[1]);
          int     idx    = ((j) * (ilast(0) + 1 - ifirst(0))) + (i);

          // Need to get some pressure samples initially to get a baseline
          // for converting to .wav files after the simulation is done.
          if (p_iteration_number < p_src_pulse_times[0]) {
            p_sampled_p0 = p_data[idx];

            // Open file for appending.
            fw1 = fopen(filename1, "a");

            if (fw1 != NULL) {
              // Output information about microphone at beginning of file.
              if (p_iteration_number == 1) {
                fprintf(fw1, "%d\n",      k + 1);
                fprintf(fw1, "%d\n",      p_wav_sample_rate);
                fprintf(fw1, "(%g,%g)\n", p_mic_pos_x[k], p_mic_pos_y[k]);
              }

              // Storing q shows us the grid level the value is recorded from.
              fprintf(fw1, "%d %d %e %e\n",
                      p_iteration_number, q, p_p0, p_sampled_p0);
              fclose(fw1);
            } else {
              TBOX_ERROR(p_object_name << ": "
                                       << "Can't open '" << filename1 << "' for recording."
                                       << endl);
            }
          }

          // Open file to record the actual pressure data.
          fw2 = fopen(filename2, "a");

          if (fw2 != NULL) {
            fprintf(fw2, "%d %d %e\n", p_iteration_number, q, p_data[idx]);
            fclose(fw2);
          } else {
            TBOX_ERROR(p_object_name << ": "
                                     << "Can't open '" << filename2 << "' for recording."
                                     << endl);
          }

          // mic_recorder[k*(p_max_timesteps+2)+p_iteration_number]=p_data[idx];
          // cout << mic_recorder[k*(p_max_timesteps+2)+p_iteration_number]
          // << endl;//-p_sampled_p0 << endl;
          // cout << p_sampled_p0 << endl;
          // cout << k*(p_max_timesteps+2)+p_iteration_number << endl;
          // hier::LocalId::LocalId mylocal=patch_iterator->getLocalId();
          // cout << "patch: " << mylocal.getValue() << endl;
          // cout << "i: " << i << " j: " << j << endl;
          // cout << "ifirst[0]: " << ifirst[0] << " ifirst[1]: " << ifirst[1]
          // << endl;
          // cout << "ilast[0]: " << ilast[0] << " ilast[1]: " << ilast[1] <<
          // endl;
          // cout << "patchXLo[0]: " << patchXLo[0] << " patchXLo[1]: " <<
          // patchXLo[1] << endl;
          // cout << "patchXHi[0]: " << patchXHi[0] << " patchXHi[1]: " <<
          // patchXHi[1] << endl;
          // cout << "worldx: " << patchXLo[0]+dx[0]*(i+0.5) << " worldy: "
          // << patchXLo[1]+dx[1]*(j+0.5) << endl;
          // cout << "p_mic_pos_x[ " << k << "]: " << p_mic_pos_x[k] << endl;
          // cout << "p_mic_pos_y[ " << k << "]: " << p_mic_pos_y[k] << endl;
          // cout << mic_recorder[p_iteration_number] << endl;
          // fflush(stdout);
        }
      }
    }
  }
  t_record->stop();
}
