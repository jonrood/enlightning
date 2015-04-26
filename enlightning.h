/*
   ===========================================================================
   Copyright (C) 2015 Jon Rood.

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

#ifndef included_enlightningXD
#define included_enlightningXD

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/tbox/MessageStream.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/appu/BoundaryUtilityStrategy.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/algs/MethodOfLinesIntegrator.h"
#include "SAMRAI/algs/MethodOfLinesPatchStrategy.h"
#include "SAMRAI/hier/Patch.h"

#include "SAMRAI/tbox/Serializable.h"
#include "SAMRAI/hier/VariableContext.h"
#include "SAMRAI/appu/VisItDataWriter.h"

#define included_String
#include <string>
#include <vector>
using namespace std;

// Define the number for the set of Navier-Stokes equations.
#define NEQU (5)

using namespace SAMRAI;

class enlightning :
  public tbox::Serializable,
  public algs::MethodOfLinesPatchStrategy,
  public appu::BoundaryUtilityStrategy {
public:

  enlightning(const string                                & object_name,
              const tbox::Dimension                       & dim,
              boost::shared_ptr<tbox::Database>             input_db,
              boost::shared_ptr<geom::CartesianGridGeometry>grid_geom);

  ~enlightning();

  void
  registerModelVariables(algs::MethodOfLinesIntegrator *integrator);

  void
  initializeDataOnPatch(hier::Patch& patch,
                        const double time,
                        const bool   initial_time) const;

  double
  computeStableDtOnPatch(hier::Patch& patch,
                         const double time) const;

  void putToRestart(const boost::shared_ptr<tbox::Database>& db) const;

  hier::IntVector getRefineOpStencilWidth(const tbox::Dimension& dim) const {
    return hier::IntVector::getOne(dim);
  }

  void
  preprocessRefine(
    hier::Patch          & fine,
    const hier::Patch    & coarse,
    const hier::Box      & fine_box,
    const hier::IntVector& ratio) {
    NULL_USE(fine);
    NULL_USE(coarse);
    NULL_USE(fine_box);
    NULL_USE(ratio);
  }

  void
  postprocessRefine(
    hier::Patch          & fine,
    const hier::Patch    & coarse,
    const hier::Box      & fine_box,
    const hier::IntVector& ratio) {
    NULL_USE(fine);
    NULL_USE(coarse);
    NULL_USE(fine_box);
    NULL_USE(ratio);
  }

  hier::IntVector
  getCoarsenOpStencilWidth(const tbox::Dimension& dim) const {
    return hier::IntVector::getZero(dim);
  }

  void
  preprocessCoarsen(
    hier::Patch          & coarse,
    const hier::Patch    & fine,
    const hier::Box      & coarse_box,
    const hier::IntVector& ratio) {
    NULL_USE(coarse);
    NULL_USE(fine);
    NULL_USE(coarse_box);
    NULL_USE(ratio);
  }

  void
  postprocessCoarsen(
    hier::Patch          & coarse,
    const hier::Patch    & fine,
    const hier::Box      & coarse_box,
    const hier::IntVector& ratio) {
    NULL_USE(coarse);
    NULL_USE(fine);
    NULL_USE(coarse_box);
    NULL_USE(ratio);
  }

  void
  singleStep(hier::Patch& patch,
             const double dt,
             const double alpha_1,
             const double alpha_2,
             const double beta) const;

  void
  tagGradientDetectorCells(hier::Patch& patch,
                           const double regrid_time,
                           const bool   initial_error,
                           const int    tag_index,
                           const bool   uses_richardson_extrapolation_too);

  void
  setPhysicalBoundaryConditions(hier::Patch& patch,
                                const double fill_time,
                                const hier::IntVector&
                                ghost_width_to_fill);

  void
  readDirichletBoundaryDataEntry(
    const boost::shared_ptr<tbox::Database>& db,
    string                                 & db_name,
    int                                      bdry_location_index) {
    NULL_USE(db);
    NULL_USE(db_name);
    NULL_USE(bdry_location_index);
  }

  void
  readNeumannBoundaryDataEntry(
    const boost::shared_ptr<tbox::Database>& db,
    string                                 & db_name,
    int                                      bdry_location_index) {
    NULL_USE(db);
    NULL_USE(db_name);
    NULL_USE(bdry_location_index);
  }

  void
  registerVisItDataWriter(boost::shared_ptr<appu::VisItDataWriter>viz_writer);

  void
  printClassData(ostream& os) const;

  int
  getMaxTimesteps();

  double
  getStartTime();

  double
  getEndTime();

  int
  getRegridStep();

  int
  getTagBuffer();

  double
  getLoopTime();

  void
  setLoopTime(const double loop_time);

  int
  getIterationNumber();

  void
  setIterationNumber(const int iter_num);

  int
  shouldRecordAudio();

  int
  getSourceMode();

  int
  shouldAddSource();

  void
  addSource(boost::shared_ptr<hier::PatchHierarchy>patch_hier,
            double dt);

  void
  setSourceToZero(boost::shared_ptr<hier::PatchHierarchy>patch_hier);

  void
  recordAudio(boost::shared_ptr<hier::PatchHierarchy>patch_hier);

private:

  void
  getFromInput(boost::shared_ptr<tbox::Database>db);

  void
  getFromRestart();

  void
  readLightningInput(const char *input_lightning_file);

  void
  setAbsorbingBoundaryConditions(hier::Patch& patch) const;

  // p_ denotes variables private to the enlightning class
  string p_object_name;
  tbox::Dimension p_dim;
  boost::shared_ptr<geom::CartesianGridGeometry> p_grid_geometry;
  boost::shared_ptr<appu::VisItDataWriter> p_visit_writer;
  boost::shared_ptr<pdat::CellVariable<double> > p_w_n;
  boost::shared_ptr<pdat::CellVariable<double> > p_K;
  boost::shared_ptr<pdat::CellVariable<double> > p_pressure;
  boost::shared_ptr<pdat::CellVariable<double> > p_temperature;
  boost::shared_ptr<pdat::CellVariable<double> > p_source;
  hier::IntVector p_nghosts;
  hier::IntVector p_zero_ghosts;
  string p_src_type;
  int    p_src_type_int;
  string p_src_mode;
  int    p_src_mode_int;
  double p_cfl;
  double p_specific_dt;
  double p_tag_tolerance[3];
  double p_start_time;
  double p_end_time;
  double p_loop_time;
  int    p_max_timesteps;
  int    p_regrid_step;
  int    p_tag_buffer;
  int    p_viz_pressure;
  int    p_viz_source;
  int    p_viz_temperature;
  int    p_viz_w_n;
  int    p_iteration_number;
  string p_bc_type_L;
  string p_bc_type_R;
  string p_bc_type_B;
  string p_bc_type_T;
  int    p_bc_type_L_int;
  int    p_bc_type_R_int;
  int    p_bc_type_B_int;
  int    p_bc_type_T_int;
  double p_bc_absorb_width;
  double p_src_amp;
  double p_src_spread;
  double p_src_pos[2];
  double p_src_sine_amp;
  double p_src_sine_freq;
  double p_src_omega;
  double p_p0;
  double p_sampled_p0;
  double p_T0;
  double p_c_lr;
  double p_rho0;
  double p_c;
  double p_gamma1;
  double p_RH;
  double p_R;
  double p_c_p;
  double p_T_star_n;
  double p_T_star_o;
  double p_mu;
  double p_muB;
  double p_kappa;
  double p_weno_alpha;
  double p_tau_n;
  double p_tau_o;
  double p_R_tilde;
  double p_c_v;
  double p_M;
  int    p_record_audio;
  int    p_num_mics;
  std::vector<double> p_mic_pos_x;
  std::vector<double> p_mic_pos_y;
  int p_src_num_pulses;
  std::vector<int> p_src_pulse_times;
  int p_use_sample_rate_dt;
  int p_wav_sample_rate;
  std::vector<double> p_lightning_coords;
  int p_lightning_input_line_count;

  static boost::shared_ptr<tbox::Timer> t_init;
  static boost::shared_ptr<tbox::Timer> t_rhs;
  static boost::shared_ptr<tbox::Timer> t_rk;
  static boost::shared_ptr<tbox::Timer> t_state;
  static boost::shared_ptr<tbox::Timer> t_absorb;
  static boost::shared_ptr<tbox::Timer> t_bc;
  static boost::shared_ptr<tbox::Timer> t_tag;
  static boost::shared_ptr<tbox::Timer> t_addsource;
  static boost::shared_ptr<tbox::Timer> t_record;
};

#endif // ifndef included_enlightningXD
