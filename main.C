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

#include "SAMRAI/SAMRAI_config.h"

#include "enlightning.h"

#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/mesh/BergerRigoutsos.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/mesh/TreeLoadBalancer.h"
#include "SAMRAI/mesh/ChopAndPackLoadBalancer.h"
#include "SAMRAI/tbox/InputDatabase.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/RestartManager.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/tbox/BalancedDepthFirstTree.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/mesh/StandardTagAndInitialize.h"
#include "SAMRAI/algs/MethodOfLinesIntegrator.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/tbox/Timer.h"
#include "SAMRAI/tbox/TimerManager.h"

#include "boost/shared_ptr.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <time.h>
using namespace std;

#include <sys/stat.h>

using namespace SAMRAI;

// Defines for source type.
#define PULSE (10)
#define DRIVEN (20)

// Number of ghost cells.
#define CELLG (3)

// Size of moving window used for estimating run time left.
#define TIME_WINDOW_SIZE (20)

int main(int argc, char *argv[])
{
  // SAMRAI initializations.
  tbox::SAMRAI_MPI::init(&argc, &argv);
  tbox::SAMRAIManager::initialize();
  tbox::SAMRAIManager::startup();
  const tbox::SAMRAI_MPI&
  mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());

  {
    // Program specific initializations.
    int time_window[TIME_WINDOW_SIZE] = { 0 };
    time_t begin, check_t1;
    string input_filename;
    string restart_read_dirname;
    int    restore_num     = 0;
    bool   is_from_restart = false;

    // Check for correct number of program arguments.
    if ((argc != 2) && (argc != 4)) {
      tbox::pout << "USAGE:  " << argv[0]
                 << " <input filename> "
                 << "<restart dir> <restore number>\n"
                 << endl;
      tbox::SAMRAI_MPI::abort();
      return -1;
    } else {
      input_filename = argv[1];

      // Check for provided restart information.
      if (argc == 4) {
        restart_read_dirname = argv[2];
        restore_num          = atoi(argv[3]);
        is_from_restart      = true;
      }
    }

    tbox::plog << "input_filename = "
               << input_filename << endl;
    tbox::plog << "restart_read_dirname = "
               << restart_read_dirname << endl;
    tbox::plog << "restore_num = "
               << restore_num << endl;

    // Create object for the input database.
    boost::shared_ptr<tbox::InputDatabase>input_db(new tbox::InputDatabase("input_db"));

    // Parse and store the input file named "input_filename" for referencing.
    tbox::InputManager::getManager()->parseInputFile(input_filename,
                                                     input_db);

    // Get object for parsing "Main" section of input file.
    boost::shared_ptr<tbox::Database> main_db =
      input_db->getDatabase("Main");

    // Set a two-dimensional simulation.
    const tbox::Dimension dim(2);

    // Get the base name of the simulation from the input file.
    const std::string base_name =
      main_db->getStringWithDefault("base_name",
                                    "unnamed");

    // Use base name in input file to name the log file.
    const std::string log_filename =
      main_db->getStringWithDefault("log_filename",
                                    base_name + ".log");

    // Generate logs from all nodes if the input file demands it, otherwise only
    // log node 0.
    bool log_all_nodes = false;

    if (main_db->keyExists("log_all_nodes")) {
      log_all_nodes = main_db->getBool("log_all_nodes");
    }

    if (log_all_nodes) {
      tbox::PIO::logAllNodes(log_filename);
    } else {
      tbox::PIO::logOnlyNodeZero(log_filename);
    }

    // Get the interval of timesteps in which visual plots will be written to
    // disk.
    int viz_dump_interval = 0;

    if (main_db->keyExists("viz_dump_interval")) {
      viz_dump_interval = main_db->getInteger("viz_dump_interval");
    }

    // Set the output directory in which visualization data output will be
    // written to.
    const std::string visit_dump_dirname =
      main_db->getStringWithDefault("viz_dump_dirname",
                                    base_name + ".visit");

    int visit_number_procs_per_file = 1;

    if (viz_dump_interval > 0) {
      if (main_db->keyExists("visit_number_procs_per_file")) {
        visit_number_procs_per_file =
          main_db->getInteger("visit_number_procs_per_file");
      }
    }

    // Used to know if viz data should be output at any certain timestep.
    const bool viz_dump_data = (viz_dump_interval > 0);

    // Set interval of timesteps in which restart (or checkpoint) data is
    // output.
    int restart_interval = 0;

    if (main_db->keyExists("restart_interval")) {
      restart_interval = main_db->getInteger("restart_interval");
    }

    // Set output directory for restart data.
    const std::string restart_write_dirname =
      main_db->getStringWithDefault("restart_write_dirname",
                                    base_name + ".restart");

    // Only write restart data if interval is set from input file and if the
    // restart directory is initially empty.
    const bool write_restart = (restart_interval > 0)
                               && !(restart_write_dirname.empty());

    // Create object for the restart manager.
    tbox::RestartManager *restart_manager =
      tbox::RestartManager::getManager();

    // If we want to restart, open the restart file for reading from the
    // specified timestep.
    if (is_from_restart) {
      restart_manager->
      openRestartFile(restart_read_dirname,
                      restore_num,
                      mpi.getSize());
    }

    // Haven't added any source to the source grid yet.
    bool added_source = false;

    // Create timer manager object and read and set options for the timer
    // manager from the input file.
    tbox::TimerManager::createManager(input_db->getDatabase("TimerManager"));

    // Create a grid object and read and set options for a base cartesian grid
    // from the input file.
    boost::shared_ptr<geom::CartesianGridGeometry>
    grid_geometry(new geom::CartesianGridGeometry(dim,
                                                  "CartesianGeometry",
                                                  input_db->getDatabase(
                                                    "CartesianGeometry")));

    // Create a patch hierarchy object and read and set options for how higher
    // grid levels are defined from the input file.
    boost::shared_ptr<hier::PatchHierarchy>
    patch_hierarchy(new hier::PatchHierarchy("PatchHierarchy",
                                             grid_geometry,
                                             input_db->getDatabase(
                                               "PatchHierarchy")));

    // Create an object specific to the enlightning simulation and read and
    // store data from "enlightning" section of the input file.
    enlightning *enlightning_model =
      new enlightning("enlightning",
                      dim,
                      input_db->getDatabase("Enlightning"),
                      grid_geometry);

    // Create an integrator object and read and store data for the integrator
    // section of the input file.
    boost::shared_ptr<algs::MethodOfLinesIntegrator>
    mol_integrator(new algs::MethodOfLinesIntegrator("MethodOfLinesIntegrator",
                                                     input_db->getDatabase(
                                                       "MethodOfLinesIntegrator"),
                                                     enlightning_model));

    // Create a cell tagging object and read and store data for the cell tagger
    // section of the input file.
    boost::shared_ptr<mesh::StandardTagAndInitialize>
    error_detector(new mesh::StandardTagAndInitialize("StandardTagAndInitialize",
                                                      mol_integrator.get(),
                                                      input_db->getDatabase(
                                                        "StandardTagAndInitialize")));

    // Create a box generator object and read and store data for the box
    // generator section of the input file.
    boost::shared_ptr<mesh::BergerRigoutsos>
    box_generator(new mesh::BergerRigoutsos(dim,
                                            input_db->getDatabaseWithDefault(
                                              "BergerRigoutsos",
                                              boost::    shared_ptr<tbox::Database>())));

    // Set chosen load balancer type from the input file.
    const std::string lb = main_db->getStringWithDefault("load_balancer", "tlb");

    // Create a tree load balancer object and read settings for the tree load
    // balancer section of the input file.
    boost::shared_ptr<mesh::TreeLoadBalancer>
    tree_load_balancer(new mesh::TreeLoadBalancer(dim,
                                                  "TreeLoadBalancer",
                                                  input_db->getDatabase(
                                                    "TreeLoadBalancer"),
                                                  boost::shared_ptr<tbox::RankTreeStrategy>(
                                                    new tbox::BalancedDepthFirstTree)));

    // Initialize MPI for the tree load balancer.
    tree_load_balancer->setSAMRAI_MPI(tbox::SAMRAI_MPI::getSAMRAIWorld());

    // Create chop and pack load balancer object and read setting for the chop
    // and pack load balancer section of the input file.
    boost::shared_ptr<mesh::ChopAndPackLoadBalancer>
    chop_load_balancer(new mesh::ChopAndPackLoadBalancer(dim,
                                                         "ChopLoadBalancer",
                                                         input_db->getDatabase(
                                                           "ChopLoadBalancer")));

    // Create object to set as the chosen load balancer.
    boost::shared_ptr<mesh::LoadBalanceStrategy> my_load_balancer;

    // Set chosen load balancer.
    if (lb == "tlb") {
      my_load_balancer = tree_load_balancer;
    } else if (lb == "caplb") {
      my_load_balancer = chop_load_balancer;
    } else {
      TBOX_ERROR("Unknown load balancer" << endl);
    }

    // Create gridding algorithm object and read settings for the gridding
    // algorithm section of the input file.
    boost::shared_ptr<mesh::GriddingAlgorithm>
    gridding_algorithm(new mesh::GriddingAlgorithm(patch_hierarchy,
                                                   "GriddingAlgorithm",
                                                   input_db->getDatabase(
                                                     "GriddingAlgorithm"),
                                                   error_detector,
                                                   box_generator,
                                                   my_load_balancer));

    // Create visit data writer object for outputing hdf5 plots.
    boost::shared_ptr<appu::VisItDataWriter>
    visit_data_writer(new appu::VisItDataWriter(dim,
                                                "Enlightning VisIt Writer",
                                                visit_dump_dirname,
                                                visit_number_procs_per_file));

    // Register the visit hdf5 writer with the simulation.
    enlightning_model->registerVisItDataWriter(visit_data_writer);

    // Log inputs used in the simulation.
    tbox::plog << "Check input data and variables before simulation:"
	       << endl << endl;

    tbox::plog << "Input database..." << endl << endl;
    input_db->printClassData(tbox::plog);

    tbox::plog << "\nVariable database..." << endl;
    hier::VariableDatabase::getDatabase()->printClassData(tbox::plog);
    tbox::plog << endl << endl;

    tbox::plog << "\nCheck Enlightning model data... " << endl;
    enlightning_model->printClassData(tbox::plog);

    // Initialize integrator.
    mol_integrator->initializeIntegrator(gridding_algorithm);

    // Create array of cell tagging objects for each level of the grid.
    std::vector<int>
    tag_buffer_array(patch_hierarchy->getMaxNumberOfLevels());

    for (int il = 0; il < patch_hierarchy->getMaxNumberOfLevels(); il++) {
      tag_buffer_array[il] = enlightning_model->getTagBuffer();
    }

    // Create object for deciding when to regrid.
    std::vector<double>
    regrid_start_time(patch_hierarchy->getMaxNumberOfLevels());

    // Initialize variable that holds the simulation time elapsed, i.e. not the
    // run time elapsed.
    double loop_time  = enlightning_model->getLoopTime();
    int    loop_cycle = enlightning_model->getIterationNumber();

    // If restarting, load the saved grid hierarchy, otherwise initialize a grid
    // hierarchy.
    if (tbox::RestartManager::getManager()->isFromRestart()) {
      patch_hierarchy->initializeHierarchy();

      // patch_hierarchy->getFromRestart();
      gridding_algorithm->getTagAndInitializeStrategy()->
      resetHierarchyConfiguration(patch_hierarchy,
                                  0,
                                  patch_hierarchy->getFinestLevelNumber());
    } else {
      gridding_algorithm->makeCoarsestLevel(loop_time);
      bool done          = false;
      bool initial_cycle = true;

      for (int ln = 0;
           patch_hierarchy->levelCanBeRefined(ln) && !done;
           ++ln) {
        gridding_algorithm->makeFinerLevel(
          tag_buffer_array[ln],
          initial_cycle,
          loop_cycle,
          loop_time);
        done = !(patch_hierarchy->finerLevelExists(ln));
      }
    }

    // Close the restart database.
    tbox::RestartManager::getManager()->closeRestartFile();

    // Set the iteration number, and save the iteration number the current run
    // started at.
    int iteration_num       = enlightning_model->getIterationNumber();
    int iteration_num_start = iteration_num;

    // Output the initial hdf5 plot before the first timestep.
    if (viz_dump_data) {
      visit_data_writer->writePlotData(patch_hierarchy,
                                       iteration_num,
                                       loop_time);
    }

    // Save the time when the run started to calculate the elapsed run time.
    time(&begin);

    if (!is_from_restart) {
      tbox::pout << "Running simulation..." << endl;
    } else {
      tbox::pout << "Restarting simulation..." << endl;
    }

    // **************** MAIN LOOP START ****************
    while ((loop_time < enlightning_model->getEndTime()) &&
           (iteration_num < enlightning_model->getMaxTimesteps())) {
      iteration_num = enlightning_model->getIterationNumber();
      iteration_num++;

      tbox::pout << "------------------------" << endl;
      tbox::pout << "step: " << iteration_num - 1 << endl;
      tbox::pout << setiosflags(ios::scientific) << setprecision(6);
      tbox::pout << "time: " << loop_time << endl;

      double dt = mol_integrator->getTimestep(patch_hierarchy, loop_time);
      tbox::pout << "dt:   " << dt << endl;
      tbox::pout << resetiosflags(ios::scientific) << setprecision(-1);

      // Add source at this timestep if necessary.
      if (enlightning_model->shouldAddSource()) {
        bool done          = false;
        bool initial_cycle = false;
        int  loop_cycle    = enlightning_model->getIterationNumber();
        added_source = true;
        tbox::pout << "adding source..." << endl;
        enlightning_model->addSource(patch_hierarchy, dt);

        // Refine the source all the way up the grid hierarchy.
        for (int ln = 0;
             patch_hierarchy->levelCanBeRefined(ln) && !done;
             ln++) {
          gridding_algorithm->makeFinerLevel(
            tag_buffer_array[ln],
            initial_cycle,
            loop_cycle,
            loop_time);
          done = !(patch_hierarchy->finerLevelExists(ln));
          enlightning_model->addSource(patch_hierarchy, dt);
        }
      }

      // Calculate the next timestep in the simulation.
      mol_integrator->advanceHierarchy(patch_hierarchy, loop_time, dt);

      loop_time += dt;

      // Tell the model what step and time the simulation is at.
      enlightning_model->setLoopTime(loop_time);
      enlightning_model->setIterationNumber(iteration_num);

      // If at a restart interval, write restart data.
      if (write_restart) {
        if ((iteration_num % restart_interval) == 0) {
          tbox::RestartManager::getManager()->
          writeRestartFile(restart_write_dirname,
                           iteration_num);
        }
      }

      // If at a visualization interval, write visualization data.
      if (viz_dump_data) {
        if (((iteration_num % viz_dump_interval) == 0)
            || ((added_source == true)
                && (enlightning_model->getSourceMode() == PULSE))) {
          visit_data_writer->writePlotData(patch_hierarchy,
                                           iteration_num,
                                           loop_time);
        }
      }

      // If a source was added, zero the source grid out.
      if (added_source) {
        enlightning_model->setSourceToZero(patch_hierarchy);
        added_source = false;
      }

      // Output microphone data if recording is requested.
      if (enlightning_model->shouldRecordAudio()) {
        enlightning_model->recordAudio(
          patch_hierarchy);
      }

      // Regrid if necessary.
      if (((iteration_num % enlightning_model->getRegridStep()) == 0) &&
          (patch_hierarchy->getMaxNumberOfLevels() > 1)) {
        tbox::plog << endl << "REGRID: #levels ";
        tbox::plog << patch_hierarchy->getFinestLevelNumber() << "->";

        gridding_algorithm->regridAllFinerLevels(0,
                                                 tag_buffer_array,
                                                 iteration_num,
                                                 loop_time,
                                                 regrid_start_time);

        tbox::plog << patch_hierarchy->getFinestLevelNumber() << endl;
        tbox::plog << endl;
      }

      // Save the time for calculating elapsed run time.
      time(&check_t1);

      // Print elapsed run time.
      int elapsed = (int)(check_t1 - begin);

      if (elapsed <= 60) {
        tbox::pout << "elapsed: " << elapsed << "s"
                   << endl;
      } else if ((elapsed > 60) && (elapsed <= 3600)) {
        tbox::pout << "elapsed: " << elapsed / 60 << "m"
                   << elapsed % 60 << "s" << endl;
      } else if ((elapsed > 3600) && (elapsed <= 86400)) {
        tbox::pout << "elapsed: " << elapsed / 3600 << "h"
                   << (elapsed % 3600) / 60 << "m"
                   << (elapsed % 3600) % 60 << "s" << endl;
      } else if (elapsed > 86400) {
        tbox::pout << "elapsed: " << elapsed / 86400 << "d"
                   << (elapsed % 86400) / 3600 << "h"
                   << ((elapsed % 86400) % 3600) / 60 << "m"
                   << ((elapsed % 86400) % 3600) % 60 << "s" << endl;
      }

      // Guess the estimated run time left of the simulation.
      // Since the time to calculate each timestep varies widely,
      // a "moving window" of TIME_WINDOW_SIZE times to calculate
      // previous timesteps is used in the estimation which
      // makes it much more accurate than an estimation using
      // the entire elapsed simulation run time. This algorithm
      // should be improved since the estimation jumps back and
      // forth due to time(&check_t1) being elapsed seconds
      // and not higher resolution. A larger TIME_WINDOW_SIZE
      // could also be used.
      double progress, time_window_percent;
      double progress_step = 100.0 * ((double)iteration_num)
                             / ((double)(enlightning_model->getMaxTimesteps()));
      double progress_time =
        100.0 * (loop_time / (enlightning_model->getEndTime()));

      if (progress_step > progress_time) {
        progress            = progress_step;
        time_window_percent = 100.0 * ((double)(TIME_WINDOW_SIZE))
                              / ((double)(enlightning_model->getMaxTimesteps()));
      } else {
        progress            = progress_time;
        time_window_percent = 100.0 * ((double)(TIME_WINDOW_SIZE)*dt)
                              / (enlightning_model->getEndTime());
      }

      // First in, first out array of stored times.
      for (int i = 0; i < (TIME_WINDOW_SIZE - 1);
           i++) time_window[i] = time_window[i + 1];

      time_window[TIME_WINDOW_SIZE - 1] = (int)(check_t1);

      tbox::pout << setiosflags(ios::fixed | ios::showpoint)
                 << setprecision(1);

      if ((iteration_num - iteration_num_start) > TIME_WINDOW_SIZE) {
        int diff_t       = time_window[TIME_WINDOW_SIZE - 1] - time_window[0];
        int seconds_left = (int)(((100.0 - progress)
                                  / time_window_percent) * ((double)(diff_t)));

        if (seconds_left <= 60) {
          tbox::pout << progress << "% done, "
                     << seconds_left << "s left" << endl;
        } else if ((seconds_left > 60) && (seconds_left <= 3600)) {
          tbox::pout << progress << "% done, "
                     << seconds_left / 60 << "m"
                     << seconds_left % 60 << "s left" << endl;
        } else if ((seconds_left > 3600) && (seconds_left <= 86400)) {
          tbox::pout << progress << "% done, "
                     << seconds_left / 3600 << "h"
                     << (seconds_left % 3600) / 60 << "m"
                     << (seconds_left % 3600) % 60 << "s left" << endl;
        } else if (seconds_left > 86400) {
          tbox::pout << progress << "% done, "
                     << seconds_left / 86400 << "d"
                     << (seconds_left % 86400) / 3600 << "h"
                     << ((seconds_left % 86400) % 3600) / 60 << "m"
                     << ((seconds_left % 86400) % 3600) % 60
                     << "s" << endl;
        }
      } else {
        tbox::pout << progress << "% done, estimating" << endl;
      }
      tbox::pout << resetiosflags(ios::fixed | ios::showpoint)
                 << setprecision(-1);
    }

    // ******************* MAIN LOOP END *******************

    // Output timer information after run is done.
    tbox::TimerManager::getManager()->print(tbox::plog);

    // Throw away all the objects.
    box_generator.reset();
    tree_load_balancer.reset();
    chop_load_balancer.reset();
    my_load_balancer.reset();
    gridding_algorithm.reset();
    visit_data_writer.reset();
    error_detector.reset();
    mol_integrator.reset();
    patch_hierarchy.reset();
    grid_geometry.reset();
    main_db.reset();
    input_db.reset();

    if (enlightning_model) {
      delete enlightning_model;
    }
  }

  // Shut down SAMRAI.
  tbox::SAMRAIManager::shutdown();
  tbox::SAMRAIManager::finalize();
  tbox::SAMRAI_MPI::finalize();

  return 0;
}
