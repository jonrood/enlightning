# Program Description

Enlightning is a program which uses the SAMRAI adaptive mesh refinement library (v3.1.0-beta at the time of writing, but recently updated to use v3.10.0) developed by Lawrence Livermore National Laboratory for simulation of the nonlinear acoustics of thunder created by lightning strikes. This software accompanies the dissertation "Numerical Simulation of the Acoustical Propagation of Thunder" by Jon Rood, 2012.

# File Descriptions:

`LICENSE.txt` - GPL license information.

`enlightning.C` - Implementations of SAMRAI routines.

`enlightning.h` - Header accompanying enlightning.C.

`enlightningFort.h` - Header for interface to the program functions that were written in Fortran.

`input.txt` - Enlightning input file that describes the parameters used in the simulation.

`lightning.txt` - File describing the vertices of the geometry of the lightning channel used in the simulation.

`main.C` - Main code that orchestrates the setup and driving of the simulation using SAMRAI.

`makefile` - Makefile for building (can pick and choose the method for solving the RHS by specifying the file for it in the makefile, for example if you wanted to solve using the DRP method by specifying the rhs_drp.m4 file).

`MOREINFO.md` - This file.

`enlightning.record/`

   * `convert_all_mics.sh` - Script for running the mic converter program on all microphone data files in the directory.

   * `makefile` - Makefile for building the mic converter program.

   * `mic.c` - Code for converting the pressure data captured by the virtual microphones, and converting it to .wav files or plot files.

`enlightning.restart/`

   * (Houses data for checkpointing simulation during runtime.)

`enlightning.viz/`

   * (Houses visualization data output during simulation runtime.)

`scripts/`

   * `enlightning.sh` - Script showing usage of executable.

   * `enlightningr.sh` - Script showing using of restarting a checkpointed simulation.

   * `saverun.sh` - Script for compressing and archiving the relevant simulation data from a run.

   * `saverunr.sh` - Same as the saverun.sh script, but includes the data for restarting from a checkpoint as well.

`fortran/`

   * `rhs_drp.m4` - This function calculates the RHS using a DRP scheme.

   * `rhs_hybrid_low.m4` - This function calculates the RHS using a hybrid scheme with a third-order WENO subscheme using the alpha from the input file for flux splitting and a DRP subscheme.

   * `rhs_weno_high_llf.m4` - This function calculates the RHS using a fifth-order WENO scheme with local lax-friedrichs flux splitting.

   * `runge_kutta.m4` - This function performs the runge kutta.

   * `rhs_hybrid_high.m4` - This function calculates the RHS using a hybrid scheme with a fifth-order WENO subscheme using the alpha from the input file for flux-splitting and a DRP subscheme.

   * `rhs_hybrid_low_llf.m4` - This function calculates the RHS using the hybrid scheme with a third-order WENO subscheme using local lax-friedrichs flux splitting and a DRP subscheme.

   * `rhs_weno_low.m4` - This function calculates the RHS using a third-order WENO scheme with the alpha from the input file for flux-splitting.

   * `tag_cells.m4` - This function calculates the smoothness or gradient in order to tag cells that require refinement.

   * `rhs_hybrid_high_llf.m4` - This function calculates the RHS using the hybrid scheme with a fifth-order WENO subscheme using local lax-friedrichs flux splitting and a DRP subscheme.

   * `rhs_weno_high.m4` - This function calculates the RHS using a fifth-order WENO scheme using the alpha value from the input file for the flux-splitting.

   * `rhs_weno_low_llf.m4` - This function calculates the RHS using a third-order WENO scheme with local lax-friedrichs flux splitting.

   * `update_state.m4` - This function updates the pressure and temperature state variables.

`matlab/`

   * `LightningChannel.m` - Matlab program for generating lightning channel geometries. Try running with “LightningChannel(2000,2000,250,2000,700,2,16,1,1,20,0.006,8,4);” for example.

   * `hybrid.m` - Matlab program that is a prototype of the hybrid method solver. Try running with “hybrid(400,51,51,26,26,1e5,1,1000,10,0.1,1,0.001);” for example.

   * `DrawSource.m` - Matlab program that is a prototype for drawing source segments on a grid. Try running with “DrawSource(1000,1000,2000,2000,6,1);” for example.

`opencl/`

   * `main.cpp` - Main program for OpenCL based prototype of the program.

   * `kernel.cl` - Kernel for calculating the RHS.

   * `makefile` - Makefile for compiling on OSX using Xcode on the command line.

`media/`

   * `mach_stem.mov` - Movie showing an example mach stem formation.

   * `mach_stem_amr.mov` - Movie showing mach stem with AMR levels and patches shown.

   * `mach_stem1.png` - Screenshot of mach stem formation.

   * `mach_stem2.png` - Zoomed in capture of AMR mesh.

   * `lightning_strike.png` - Example of a lightning channel simulation.

`paper/

   * `jr_dissertation` - Dissertation detailing the simulation.
