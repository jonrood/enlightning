# Installing Enlightning and its Dependencies

Enlightning relies upon the SAMRAI software library so you will need to install it. On Linux I would recommend installing SAMRAI using [Spack](https://github.com/llnl/spack), and you should basically be able to do `spack install samrai`, and all the dependencies should install as well. Then you can follow the Mac instructions at step 6 to install Enlightning itself. On the Mac I would recommend installing with Homebrew.

## OSX 10.12 general install notes (NOT NECESSARILY EXACT):

1. Install gcc and gfortran using Homebrew (tested with GCC 6.3.0):

   `brew install gcc`

2. Install MPI (tested with mpich 3.2.2) using Homebrew:

   `brew install mpich`

   You may have to add your `hostname -f` output to /etc/hosts as '127.0.0.1'

3. Install hdf5 (tested with 1.10.0) and boost (tested with 1.63.0) using Homebrew:

   `brew install homebrew/science/hdf5 --with-fortran --with-mpi --with-tests`

   `brew install boost`

4. Install SAMRAI:

   `mkdir ~/SAMRAI-v3.11.2`

   `mv ~/Downloads/SAMRAI-v3.11.2.tar.gz ~/SAMRAI-v3.11.2/`

   `cd ~/SAMRAI-v3.11.2/`

   `tar xf SAMRAI-v3.11.2.tar.gz`

   `mkdir build`

   `cd build`

   `../SAMRAI/configure --prefix=${HOME}/SAMRAI-v3.11.2/install --with-hdf5=/usr/local --with-CC=mpicc --with-CXX=mpicxx --with-F77=mpif90 --with-zlib --without-petsc --enable-opt=-O3 --disable-debug --with-boost=/usr/local`

   `make library`

   `make install`

5. You may need to fix an error in SAMRAI's Fortran location:

   Change the line `INCLUDE_SAM   = $(SAMRAI)/source` in `~/SAMRAI-v3.11.2/install/config/Makefile.config` to be `INCLUDE_SAM   = $(SAMRAI)/include`.

6. Install Enlightning:

   Change the SAMRAI variable in the Enlightning makefile to the location of SAMRAI source, i.e. `~/SAMRAI-v3.11.2/install`

   Then issue the command `make` in the enlightning directory.
