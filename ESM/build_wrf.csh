#!/bin/csh -f
#
# git $Id$
# svn $Id: build_wrf.csh 1202 2023-10-24 15:36:07Z arango $
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Copyright (c) 2002-2023 The ROMS/TOMS Group                           :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::: Hernan G. Arango :::
#                                                                       :::
# WRF Compiling CSH Script: WRF Versions 4.5.2 and higher               :::
#                                                                       :::
# Script to compile WRF where source code files are kept separate       :::
# from the application configuration and build objects.                 :::
#                                                                       :::
# Q: How/why does this script work?                                     :::
#                                                                       :::
# A:                                                                    :::
#                                                                       :::
# Usage:                                                                :::
#                                                                       :::
#    ./build_wrf.csh [options]                                          :::
#                                                                       :::
# Options:                                                              :::
#                                                                       :::
#    -j [N]         Compile in parallel using N CPUs                    :::
#                     omit argument for all available CPUs              :::
#                                                                       :::
#    -b [branch]    Compile myroms forked WRF GitHub branch             :::
#                                                                       :::
#                     build_wrf.csh -j 5 -b                             :::
#                     build_wrf.csh -j 5 -b feature/coupling            :::
#                                                                       :::
#    -noclean       Do not run clean -a script                          :::
#                                                                       :::
#    -noconfig      Do not run configure compilation script             :::
#                                                                       :::
# The branch option -b will download the WRF code from the myroms fork  :::
# of WRF (https://github.com/myroms/WRF) and checkout the given branch  :::
# or the 'Coupling' branch if none is provided.                         :::
#                                                                       :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

setenv which_MPI openmpi                        # default, overwriten below

# Initialize.

set separator = `perl -e "print ':' x 100;"`

set parallel = 0
set clean = 1
set config = 1
set branch = 0

setenv CPLFLAG ''
setenv MY_CPP_FLAGS ''

while ( ($#argv) > 0 )
  switch ($1)
    case "-j"
      shift
      set parallel = 1
      if (`echo $1 | grep '^[0-9]\+$'` != "" ) then
        set NCPUS = $1
        shift
      else
        set NCPUS = 2
      endif
    breaksw

    case "-b"
      shift
      set branch = 1
      set move = 1
      set branch_name = `echo $1 | grep -v '^-'`
      if ( "$branch_name" == "" ) then
        set branch_name = 'Coupling'
      else
        shift
      endif
      echo "Using myroms WRF $branch_name branch"
    breaksw

    case "-noclean"
      shift
      set clean = 0
    breaksw

    case "-noconfig"
      shift
      set config = 0
    breaksw

    case "-*":
      echo ""
      echo "$0 : Unknown option [ $1 ]"
      echo ""
      echo "Available Options:"
      echo ""
      echo "-j [N]            Compile in parallel using N CPUs"
      echo ""
      echo "-b [branch_name]  Compile specific myroms forked WRF GitHub branch"
      echo "                    For example:  build_wrf.csh -b Coupling"
      echo ""
      echo "-noclean          Do not run clean script"
      echo ""
      echo "-nocconfig        Do not run configure compilation script"
      echo ""
      exit 1
    breaksw

  endsw
end

#--------------------------------------------------------------------------
# Set a local environmental variable to define the root path to the
# directories where ROMS and WRF source files are located.  The ROMS
# source directory is needed for replacing several WRF files for
# ESMF/NUOPC Coupling (see below).
#--------------------------------------------------------------------------

 setenv ROMS_SRC_DIR         ${HOME}/ocean/repository/git/roms

 setenv WRF_ROOT_DIR         ${HOME}/ocean/repository/git/WRF

#--------------------------------------------------------------------------
# Set a local environmental variable to define the path of the working
# application directory where all this project's files are kept.
#--------------------------------------------------------------------------

 setenv MY_PROJECT_DIR       ${PWD}
 setenv MY_PROJECT_DATA      `dirname ${PWD}`/Data

#--------------------------------------------------------------------------
# Set Fortran compiler and MPI library to use.
#--------------------------------------------------------------------------

# Set path of the directory containing makefile configuration (*.mk) files.
# The user has the option to specify a customized version of these files
# in a different directory than the one distributed with the source code,
# ${COAMPS_ROOT_DIR}/compilers. If this is the case, the you need to keep
# these configurations files up-to-date.

 setenv USE_MPI              on          # distributed-memory parallelism
 setenv USE_MPIF90           on          # compile with mpif90 script
#setenv which_MPI            intel       # compile with mpiifort library
#setenv which_MPI            mpich       # compile with MPICH library
#setenv which_MPI            mpich2      # compile with MPICH2 library
#setenv which_MPI            mvapich2    # compile with MVAPICH2 library
 setenv which_MPI            openmpi     # compile with OpenMPI library

 setenv FORT                 ifort
#setenv FORT                 gfortran
#setenv FORT                 pgi

#setenv USE_REAL_DOUBLE      on          # use real double precision (-r8)
#setenv USE_DEBUG            on          # use Fortran debugging flags
 setenv USE_NETCDF           on          # compile with NetCDF
 setenv USE_NETCDF4          on          # compile with NetCDF-4 library
                                         # (Must also set USE_NETCDF)

 unsetenv PNETCDF                        # disable compiling with pNetCDF

#--------------------------------------------------------------------------
# Use my specified library paths. It is not needed but it is added for
# for checking in the future.
#--------------------------------------------------------------------------

 setenv USE_MY_LIBS no           # use system default library paths
#setenv USE_MY_LIBS yes          # use my customized library paths

 set MY_PATHS = ${ROMS_SRC_DIR}/Compilers/my_build_paths.csh
#set MY_PATHS = ${HOME}/Compilers/ROMS/my_build_paths.csh

if ($USE_MY_LIBS == 'yes') then
  source ${MY_PATHS} ${MY_PATHS}
endif

#--------------------------------------------------------------------------
# WRF build and executable directory.
#--------------------------------------------------------------------------

# Put the *.a, .f and .f90 files in a project specific Build directory to
# avoid conflict with other projects.

if ($?USE_DEBUG) then
  setenv  WRF_BUILD_DIR   ${MY_PROJECT_DIR}/Build_wrfG
else
  setenv  WRF_BUILD_DIR   ${MY_PROJECT_DIR}/Build_wrf
endif

# Put WRF executables in the following directory.

setenv  WRF_BIN_DIR       ${WRF_BUILD_DIR}/Bin

#--------------------------------------------------------------------------
# Download WRF code if requested
#--------------------------------------------------------------------------

if ( $branch == 1 ) then
  if ( ! -d ${MY_PROJECT_DIR}/wrf ) then
    echo ""
    echo "Downloading WRF source code from myroms GitHub: https://www.github.com/myroms/WRF"
    echo ""
    git clone https://www.github.com/myroms/WRF.git wrf
  endif
  echo ""
  echo "Checking out myroms WRF GitHub branch: $branch_name"
  echo ""
  cd wrf
  git checkout $branch_name
  setenv WRF_ROOT_DIR ${PWD}
else
  echo ""
  echo "Using ROMS source code from: ${MY_ROMS_SRC}"
  echo ""
endif

# Decode WRF version from its README file to decide the appropriate data
# file links needed.

 set wrf_ver = `grep 'WRF Model Version ' ${WRF_ROOT_DIR}/README | sed -e 's/[^0-9]*\([0-9]\.[0-9]\).*/\1/'`
 setenv WRF_VERSION          ${wrf_ver}

# Go to the users source directory to compile. The options set above will
# pick up the application-specific code from the appropriate place.

 cd ${WRF_ROOT_DIR}

#--------------------------------------------------------------------------
# Configure. It creates configure.wrf script used for compilation.
#
#   -d    build with debugging information and no optimization
#   -D    build with -d AND floating traps, traceback, uninitialized variables
#   -r8   build with 8 byte reals
#
# During configuration the WRF/arch/Config.pl perl script is executed and
# we need to interactively select the combination of compiler and parallel
# comunications option. For example, for Darwin operating system we get:
#
# Please select from among the following Darwin ARCH options:
#
#  1. (serial)   2. (smpar)   3. (dmpar)   4. (dm+sm)   PGI (pgf90/pgcc)
#  5. (serial)   6. (smpar)   7. (dmpar)   8. (dm+sm)   INTEL (ifort/icc)
#  9. (serial)  10. (smpar)  11. (dmpar)  12. (dm+sm)   INTEL (ifort/clang)
# 13. (serial)               14. (dmpar)                GNU (g95/gcc)
# 15. (serial)  16. (smpar)  17. (dmpar)  18. (dm+sm)   GNU (gfortran/gcc)
# 19. (serial)  20. (smpar)  21. (dmpar)  22. (dm+sm)   GNU (gfortran/clang)
# 23. (serial)               24. (dmpar)                IBM (xlf90_r/cc)
# 25. (serial)  26. (smpar)  27. (dmpar)  28. (dm+sm)   PGI (pgf90/pgcc): -f90=pgf90
# 29. (serial)  30. (smpar)  31. (dmpar)  32. (dm+sm)   INTEL (ifort/icc): Open MPI
# 33. (serial)  34. (smpar)  35. (dmpar)  36. (dm+sm)   INTEL/GNU (ifort/gcc): Open MPI
# 37. (serial)  38. (smpar)  39. (dmpar)  40. (dm+sm)   GNU (gfortran/gcc): Open MPI
#
# Enter selection [1-40] : ??
#
# For coupling with ESMF/NUOPC, we need select an option from the (dmpar)
# column for distributed-memory configuration.
#--------------------------------------------------------------------------

# Clean source code.

if ( $clean == 1 ) then
  echo ""
  echo "${separator}"
  echo "Cleaning WRF source code:  ${WRF_ROOT_DIR}/clean -a"
  echo "${separator}"
  echo ""
  ${WRF_ROOT_DIR}/clean -a
endif

if ($?USE_DEBUG) then
#  set DEBUG_FLAG = -d
   set DEBUG_FLAG = -D
endif

setenv CONFIG_FLAGS ''

if ( $config == 1 ) then
  if ($?USE_DEBUG && $?USE_REAL_DOUBLE) then
    setenv CONFIG_FLAGS   "${DEBUG_FLAG} -r8"
  else if ($?USE_DEBUG) then
    setenv CONFIG_FLAGS   ${DEBUG_FLAG}
  else if ($?USE_REAL_DOUBLE) then
    setenv CONFIG_FLAGS   -r8
  endif

  echo ""
  echo "${separator}"
  echo "Configuring WRF code:  ${WRF_ROOT_DIR}/configure ${CONFIG_FLAGS}"
  echo "${separator}"
  echo ""

  ${WRF_ROOT_DIR}/configure ${CONFIG_FLAGS}

# If which_MPI is "intel" then we need to replace DM_FC and DM_CC in configure.wrf

  if ( "${which_MPI}" == "intel" ) then
    perl -i -pe 's/^DM_FC(\s*)=(\s*)mpif90/DM_FC$1=$2mpiifort/' ${WRF_ROOT_DIR}/configure.wrf
    perl -i -pe 's/^DM_CC(\s*)=(\s*)mpicc/DM_CC$1=$2mpiicc/' ${WRF_ROOT_DIR}/configure.wrf
  endif
endif

#--------------------------------------------------------------------------
# WRF Compile script:
#
# Usage:
#
#     compile [-j n] wrf   compile wrf in run dir (NOTE: no real.exe,
#                          ndown.exe, or ideal.exe generated)
#
#   or choose a test case (see README_test_cases for details) :
#
#     compile [-j n] em_b_wave
#     compile [-j n] em_convrad
#     compile [-j n] em_esmf_exp
#     compile [-j n] em_fire
#     compile [-j n] em_grav2d_x
#     compile [-j n] em_heldsuarez
#     compile [-j n] em_hill2d_x
#     compile [-j n] em_les
#     compile [-j n] em_quarter_ss
#     compile [-j n] em_real
#     compile [-j n] em_scm_xy
#     compile [-j n] em_seabreeze2d_x
#     compile [-j n] em_squall2d_x
#     compile [-j n] em_squall2d_y
#     compile [-j n] em_tropical_cyclone
#     compile [-j n] nmm_real
#     compile [-j n] nmm_tropical_cyclone
#
#     compile -j n         parallel make using n tasks if supported
#                          (default 2)
#
#     compile -h           help message
#--------------------------------------------------------------------------

setenv WRF_DA_CORE 0             # no Data Assimilation core
setenv WRF_EM_CORE 1             # Eurelian Mass-coordinate core
setenv WRF_NMM_CORE 0            # Nonhydrostatic Mesoscale Model core

# Remove existing build directory.

if ( $clean == 1 ) then
  echo ""
  echo "${separator}"
  echo "Removing WRF build directory:  ${WRF_BUILD_DIR}"
  echo "${separator}"
  echo ""
  /bin/rm -rf ${WRF_BUILD_DIR}
endif

# Compile (if -move is set, the binaries will go to WRF_BIN_DIR set above).

#set WRF_CASE = wrf
 set WRF_CASE = em_real

if ( $parallel == 1 ) then
  setenv J "-j ${NCPUS}"
else
  setenv J "-j 2"
endif

echo ""
echo "${separator}"
echo "Compiling WRF using  ${MY_PROJECT_DIR}/${0}:"
echo ""
echo "   ${WRF_ROOT_DIR}/compile ${WRF_CASE}"
echo "        WRF_DA_CORE = ${WRF_DA_CORE},    Data Assimilation core"
echo "        WRF_EM_CORE = ${WRF_EM_CORE},    Eurelian Mass-coordinate core"
echo "        WRF_NMM_CORE = ${WRF_NMM_CORE},   Nonhydrostatic Mesoscale Model core"
echo "        J = ${J},          number of compiling CPUs"
echo "${separator}"
echo ""

${WRF_ROOT_DIR}/compile ${WRF_CASE}

#--------------------------------------------------------------------------
# Move WRF objects and executables to working project directory.
#--------------------------------------------------------------------------

${ROMS_SRC_DIR}/ESM/wrf_move.csh

#--------------------------------------------------------------------------
# Create WRF data links to working project directory.
#--------------------------------------------------------------------------

if ( $WRF_CASE == "em_real" ) then
  ${ROMS_SRC_DIR}/ESM/wrf_links.csh
endif
