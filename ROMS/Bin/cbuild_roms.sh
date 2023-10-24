#!/bin/bash
#
# git $Id$
# svn $Id: cbuild_roms.sh 1202 2023-10-24 15:36:07Z arango $
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Copyright (c) 2002-2023 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.md                                                 :::
#:::::::::::::::::::::::::::::::::::::::::::::::::::::: David Robertson :::
#                                                                       :::
# ROMS CMake Compiling BASH Script                                      :::
#                                                                       :::
# Script to configure and compile a user application where the          :::
# application-specific files are kept separate from the ROMS            :::
# source code.                                                          :::
#                                                                       :::
# Q: How/why does this script work?                                     :::
#                                                                       :::
# A: The ROMS makefile configures user-defined options with a set of    :::
#    flags such as ROMS_APPLICATION. Browse the makefile to see these.  :::
#    If an option in the makefile uses the syntax ?= in setting the     :::
#    default, this means that make will check whether an environment    :::
#    variable by that name is set in the shell that calls make. If so   :::
#    the environment variable value overrides the default (and the      :::
#    user need not maintain separate makefiles, or frequently edit      :::
#    the makefile, to run separate applications).                       :::
#                                                                       :::
# Usage:                                                                :::
#                                                                       :::
#    ./cbuild_roms.sh [options]                                         :::
#                                                                       :::
# Options:                                                              :::
#                                                                       :::
#    -j [N]         Compile in parallel using N CPUs                    :::
#                     omit argument for all available CPUs              :::
#                                                                       :::
#    -b             Compile a specific ROMS GitHub branch               :::
#                                                                       :::
#                     cbuild_roms.sh -j 5 -b feature/kernel             :::
#                                                                       :::
#    -p macro       Prints any Makefile macro value. For example,       :::
#                                                                       :::
#                     cbuild_roms.sh -p MY_CPP_FLAGS                    :::
#                                                                       :::
#    -noclean       Do not clean already compiled objects               :::
#                                                                       :::
#    -v             Compile in verbose mode (VERBOSE=1)                 :::
#                                                                       :::
# The branch option -b is only possible for ROMS source code from       :::
# https://github.com/myroms. Such versions are under development        :::
# and targeted to advanced users, superusers, and beta testers.         :::
# Regular and novice users must use the default 'develop' branch.       :::
#                                                                       :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

export which_MPI=openmpi                       # default, overwritten below

parallel=0
clean=1
dprint=0
Verbose=0
branch=0

command="cbuild_roms.sh $@"

separator=`perl -e "print '<>' x 50;"`

export MY_CPP_FLAGS=

while [ $# -gt 0 ]
do
  case "$1" in
    -j )
      shift
      parallel=1
      test=`echo $1 | grep '^[0-9]\+$'`
      if [ "$test" != "" ]; then
        NCPUS="-j $1"
        shift
      else
        NCPUS="-j"
      fi
      ;;

    -p )
      shift
      clean=0
      dprint=1
      debug="$1"
      shift
      ;;

    -v )
      shift
      Verbose=1
      ;;

    -noclean )
      shift
      clean=0
      ;;

    -b )
      shift
      branch=1
      branch_name=`echo $1 | grep -v '^-'`
      if [ "$branch_name" == "" ]; then
        echo "Please enter a ROMS GitHub branch name."
        exit 1
      fi
      shift
      ;;

    * )
      echo ""
      echo "${separator}"
      echo "$0 : Unknown option [ $1 ]"
      echo ""
      echo "Available Options:"
      echo ""
      echo "-j [N]          Compile in parallel using N CPUs"
      echo "                  omit argument for all avaliable CPUs"
      echo ""
      echo "-b branch_name  Compile specific ROMS GitHub branch name"
      echo "                  For example:  cbuild_roms.sh -b feature/kernel"
      echo ""
      echo "-p macro        Prints any Makefile macro value"
      echo "                  For example:  cbuild_roms.sh -p FFLAGS"
      echo ""
      echo "-noclean        Do not clean already compiled objects"
      echo ""
      echo "-v              Compile in verbose mode"
      echo "${separator}"
      echo ""
      exit 1
      ;;
  esac
done

# Set the CPP option defining the particular application. This will
# determine the name of the ".h" header file with the application
# CPP definitions. REQUIRED

export   ROMS_APPLICATION=UPWELLING

# Set a local environmental variable to define the path to the directories
# where the ROMS source code is located (MY_ROOT_DIR), and this project's
# configuration and files are kept (MY_PROJECT_DIR). Notice that if the
# User sets the ROMS_ROOT_DIR environment variable in their computer logging
# script describing the location from where the ROMS source code was cloned
# or downloaded, it uses that value.

if [ -n "${ROMS_ROOT_DIR:+1}" ]; then
  export      MY_ROOT_DIR=${ROMS_ROOT_DIR}
else
  export      MY_ROOT_DIR=${HOME}/ocean/repository/git
fi

export     MY_PROJECT_DIR=${PWD}

# The path to the user's local current ROMS source code.
#
# If downloading ROMS locally, this would be the user's Working Copy Path.
# One advantage of maintaining your source code copy is that when working
# simultaneously on multiple machines (e.g., a local workstation, a local
# cluster, and a remote supercomputer), you can update with the latest ROMS
# release and always get an up-to-date customized source on each machine.
# This script allows for differing paths to the code and inputs on other
# computers.

 export       MY_ROMS_SRC=${MY_ROOT_DIR}/roms

# Which type(s) of libraries would you like?
#
# NOTE: If you choose both and also choose to build the ROMS executable,
#       it will be linked to the shared version of the library.
#
# Valid options are SHARED, STATIC, and BOTH.

 export           LIBTYPE=STATIC

# Do you want to build the ROMS executable?
#
# Valid values are: ON (build the executable) and OFF (do NOT build the
# executable). If you comment this out the executable WILL be built.

 export   ROMS_EXECUTABLE=ON

# Set path of the directory containing "my_build_paths.sh".
#
# The user has the option to specify a customized version of this file
# in a different directory than the one distributed with the source code,
# ${MY_ROMS_SRC}/Compilers. If this is the case, you need to keep these
# configurations files up-to-date.

 export         COMPILERS=${MY_ROMS_SRC}/Compilers
#export         COMPILERS=${HOME}/Compilers/ROMS

#--------------------------------------------------------------------------
# Set tunable CPP options.
#--------------------------------------------------------------------------
#
# Sometimes it is desirable to activate one or more CPP options to run
# different variants of the same application without modifying its header
# file. If this is the case, specify each option here.
#
# Notice also that you need to use shell's quoting syntax to enclose the
# definition.
#
#export      MY_CPP_FLAGS="${MY_CPP_FLAGS} -DAVERAGES"
#export      MY_CPP_FLAGS="${MY_CPP_FLAGS} -DDEBUGGING"
#
# can be used to write time-averaged fields. Notice that you can have as
# many definitions as you want by appending values.

#--------------------------------------------------------------------------
# Compilation options.
#--------------------------------------------------------------------------

 export           USE_MPI=on               # distributed-memory parallelism
 export        USE_MPIF90=on               # compile with mpif90 script
#export         which_MPI=intel            # compile with mpiifort library
#export         which_MPI=mpich            # compile with MPICH library
#export         which_MPI=mpich2           # compile with MPICH2 library
#export         which_MPI=mvapich2         # compile with MVAPICH2 library
 export         which_MPI=openmpi          # compile with OpenMPI library

 export              FORT=ifort
#export              FORT=gfortran
#export              FORT=pgi

#export         USE_DEBUG=on               # use Fortran debugging flags

# ROMS I/O choices and combinations. A more complete description of the
# available options can be found in the wiki (https://myroms.org/wiki/IO).
# Most users will want to enable at least USE_NETCDF4 because that will
# instruct the ROMS build system to use nf-config to determine the
# necessary libraries and paths to link into the ROMS executable.

 export       USE_NETCDF4=on               # compile with NetCDF-4 library
#export   USE_PARALLEL_IO=on               # Parallel I/O with NetCDF-4/HDF5
#export           USE_PIO=on               # Parallel I/O with PIO library
#export       USE_SCORPIO=on               # Parallel I/O with SCORPIO library

# If any of the coupling component use the HDF5 Fortran API for primary
# I/O, we need to compile the main driver with the HDF5 library.

#export          USE_HDF5=on               # compile with HDF5 library

#--------------------------------------------------------------------------
# If coupling Earth System Models (ESM), set the location of the ESM
# component libraries and modules.
#--------------------------------------------------------------------------

source ${MY_ROMS_SRC}/ESM/esm_libs.sh ${MY_ROMS_SRC}/ESM/esm_libs.sh

#--------------------------------------------------------------------------
# If applicable, use my specified library paths.
#--------------------------------------------------------------------------

 export USE_MY_LIBS=no            # use system default library paths
#export USE_MY_LIBS=yes           # use my customized library paths

MY_PATHS=${COMPILERS}/my_build_paths.sh

if [ "${USE_MY_LIBS}" = "yes" ]; then
  source ${MY_PATHS} ${MY_PATHS}
fi

# Set location of the application header file.

 export     MY_HEADER_DIR=${MY_PROJECT_DIR}

# If you have custom analytical functions to include, enter the path here.

 export MY_ANALYTICAL_DIR=${MY_PROJECT_DIR}

 echo ""
 echo "${separator}"

# Put the CMake files in a project specific Build directory to avoid conflict
# with other projects.

if [ "${USE_DEBUG:-x}" = "on" ]; then
  export        BUILD_DIR=${MY_PROJECT_DIR}/CBuild_romsG
else
  if [ "${USE_MPI:-x}" = "on" ]; then
    export      BUILD_DIR=${MY_PROJECT_DIR}/CBuild_romsM
  else
    export      BUILD_DIR=${MY_PROJECT_DIR}/CBuild_roms
  fi
fi

# For backward compatibility, set deprecated SCRATCH_DIR to compile
# older released versions of ROMS.

export SCRATCH_DIR=${BUILD_DIR}

# If necessary, create ROMS build directory.

if [ $dprint -eq 0 ]; then
  if [ -d ${BUILD_DIR} ]; then
    if [ $clean -eq 1 ]; then
      echo ""
      echo "Removing ROMS build directory: ${BUILD_DIR}"
      echo ""
      rm -rf ${BUILD_DIR}
      echo ""
      echo "Creating ROMS build directory: ${BUILD_DIR}"
      echo ""
      mkdir ${BUILD_DIR}
    fi
  else
    if [ $clean -eq 1 ]; then
      mkdir ${BUILD_DIR}
      cd ${BUILD_DIR}
    else
      echo ""
      echo "Option -noclean activated when the ROMS build directory did not exist"
      echo "Creating ROMS build directory and disabling -noclean"
      echo ""
      clean=1
      mkdir ${BUILD_DIR}
      cd ${BUILD_DIR}
    fi
  fi

  # If requested, check out requested branch from ROMS GitHub

  if [ $branch -eq 1 ]; then
    if [ ! -d ${MY_PROJECT_DIR}/src ]; then
      echo ""
      echo "Downloading ROMS source code from GitHub: https://www.github.com/myroms"
      echo ""
      git clone https://www.github.com/myroms/roms.git src
    fi
    echo ""
    echo "Checking out ROMS GitHub branch: $branch_name"
    echo ""
    cd src
    git checkout $branch_name

    # If we are using the COMPILERS from the ROMS source code
    # overide the value set above

    if [[ ${COMPILERS} == ${MY_ROMS_SRC}* ]]; then
      export COMPILERS=${MY_PROJECT_DIR}/src/Compilers
    fi
    export MY_ROMS_SRC=${MY_PROJECT_DIR}/src

  else
    echo ""
    echo "Using ROMS source code from: ${MY_ROMS_SRC}"
    echo ""
    cd ${MY_ROMS_SRC}
  fi
fi

#--------------------------------------------------------------------------
# Add environmental variables constructed in 'makefile' to MY_CPP_FLAGS
# so can be passed to ROMS.
#--------------------------------------------------------------------------

ANALYTICAL_DIR="ANALYTICAL_DIR='${MY_ANALYTICAL_DIR}'"
HEADER=`echo ${ROMS_APPLICATION} | tr '[:upper:]' '[:lower:]'`.h
HEADER_DIR="HEADER_DIR='${MY_HEADER_DIR}'"
ROOT_DIR="ROOT_DIR='${MY_ROMS_SRC}'"

mycppflags="${MY_CPP_FLAGS}"

export       MY_CPP_FLAGS="${MY_CPP_FLAGS} -D${ANALYTICAL_DIR}"
export       MY_CPP_FLAGS="${MY_CPP_FLAGS} -D${HEADER_DIR}"
export       MY_CPP_FLAGS="${MY_CPP_FLAGS} -D${ROOT_DIR}"

cd ${BUILD_DIR}

#--------------------------------------------------------------------------
# Configure.
#--------------------------------------------------------------------------

# Construct the "cmake" command.

if [ ! -z "${LIBTYPE}" ]; then
  ltype="-DLIBTYPE=${LIBTYPE}"
else
  ltype=""
fi

if [ ! -z "${FORT}" ]; then
  if [ ${FORT} == "ifort" ]; then
    compiler="-DCMAKE_Fortran_COMPILER=ifort"
  elif [ ${FORT} == "gfortran" ]; then
    compiler="-DCMAKE_Fortran_COMPILER=gfortran"
  else
    compiler=""
  fi
fi

if [ ! -z "${MY_CPP_FLAGS}" ]; then
  tmp=`echo ${MY_CPP_FLAGS} | sed 's/^ *-D//' | sed 's/ *-D/;/g'`
  extra_flags="-DMY_CPP_FLAGS=${tmp}"
else
  extra_flags=""
fi

if [ ! -z "${PARPACK_LIBDIR}" ]; then
  parpack_ldir="-DPARPACK_LIBDIR=${PARPACK_LIBDIR}"
else
  parpack_ldir=""
fi

if [ ! -z "${ARPACK_LIBDIR}" ]; then
  arpack_ldir="-DARPACK_LIBDIR=${ARPACK_LIBDIR}"
else
  arpack_ldir=""
fi

if [ ! -z "${USE_SCORPIO}" ]; then
  if [[ ! -z "${PIO_LIBDIR}" && ! -z "${PIO_INCDIR}" ]]; then
    pio_ldir="-DPIO_LIBDIR=${PIO_LIBDIR}"
    pio_idir="-DPIO_INCDIR=${PIO_INCDIR}"
    if [[ ! -z "${PNETCDF_LIBDIR}" && ! -z "${PNETCDF_INCDIR}" ]]; then
      pnetcdf_ldir="-DPNETCDF_LIBDIR=${PNETCDF_LIBDIR}"
      pnetcdf_idir="-DPNETCDF_INCDIR=${PNETCDF_INCDIR}"
    else
      pnetcdf_ldir=""
      pnetcdf_idir=""
    fi
  else
    pio_ldir=""
    pio_idir=""
    pnetcdf_ldir=""
    pnetcdf_idir=""
  fi
fi

if [[ ! -z "${USE_MPI}" && "${USE_MPI}" == "on" ]]; then
  mpi="-DMPI=ON"
else
  mpi=""
fi

if [[ ! -z "${USE_MPIF90}" && "${USE_MPIF90}" == "on" ]]; then
  comm="-DCOMM=${which_MPI}"
else
  comm=""
fi

if [ ! -z "${ROMS_EXECUTABLE}" ]; then
  if [[ "${ROMS_EXECUTABLE}" == "ON" ]]; then
    roms_exec="-DROMS_EXECUTABLE=ON"
  else
    roms_exec="-DROMS_EXECUTABLE=OFF"
  fi
else
  roms_exec=""
fi

if [[ ! -z "${USE_DEBUG}" && "${USE_DEBUG}" == "on" ]]; then
  dbg="-DCMAKE_BUILD_TYPE=Debug"
else
  dbg="-DCMAKE_BUILD_TYPE=Release"
fi

#--------------------------------------------------------------------------
# Run the chosen build command.
#--------------------------------------------------------------------------

my_hdir="-DMY_HEADER_DIR=${MY_HEADER_DIR}"

if [[ $dprint -eq 0 && $clean -eq 1 ]]; then
  echo ""
  echo "Configuring CMake for ROMS application:"
  echo ""
  cmake -DROMS_APP=${ROMS_APPLICATION} \
                   ${my_hdir} \
                   ${ltype} \
                   ${compiler} \
                   ${extra_flags} \
                   ${parpack_ldir} \
                   ${arpack_ldir} \
                   ${pio_ldir} \
                   ${pio_idir} \
                   ${pnetcdf_ldir} \
                   ${pnetcdf_idir} \
                   ${mpi} \
                   ${comm} \
                   ${roms_exec} \
                   ${dbg} \
                   ${MY_ROMS_SRC}
fi

if [ $? -ne 0 ]; then
  echo "cmake did not complete successfully"
  exit 1
fi

#--------------------------------------------------------------------------
# Compile.
#--------------------------------------------------------------------------

if [ $dprint -eq 1 ]; then
  echo $debug:"${!debug}"
else

  echo ""
  echo "Compiling ROMS source code:"
  echo ""

  if [ $parallel -eq 1 ]; then
    if [ $Verbose -eq 1 ]; then
      make VERBOSE=1 $NCPUS
    else
      make $NCPUS
    fi
  else
    if [ $Verbose -eq 1 ]; then
      make VERBOSE=1
    else
      make
    fi
  fi
  make install

  echo ""
  echo "${separator}"
  echo "CMake Build script command:    ${command}"
  echo "ROMS source directory:         ${MY_ROMS_SRC}"
  echo "ROMS build  directory:         ${BUILD_DIR}"
  if [ $branch -eq 1 ]; then
    echo "ROMS downloaded from:          https://github.com/myroms/roms.git"
    echo "ROMS compiled branch:          $branch_name"
  fi
  echo "ROMS Application:              ${ROMS_APPLICATION}"
  FFLAGS=`cat fortran_flags`
  echo "Fortran compiler:              ${FORT}"
  echo "Fortran flags:                ${FFLAGS}"
  if [ -n "${mycppflags:+1}" ]; then
    echo "Added CPP Options:            ${mycppflags}"
  fi
  echo "${separator}"
  echo ""
fi

cd ${MY_PROJECT_DIR}

# If ROMS_EXECUTABLE is set to OFF remove the symlink from
# previous builds if present.

if [[ ! -z "${ROMS_EXECUTABLE}" && "${ROMS_EXECUTABLE}" == "OFF" ]]; then
  if [[ ! -z "${USE_DEBUG}" && "${USE_DEBUG}" == "on" ]]; then
    if [ -L romsG ]; then
      rm -f romsG
    fi
  elif [[ ! -z "${USE_MPI}" && "${USE_MPI}" == "on" ]]; then
    if [ -L romsM ]; then
      rm -f romsM
    fi
  else
    if [ -L romsS ]; then
      rm -f romsS
    fi
  fi
fi

# Copy executable to project directory. This should work even
# if ROMS was linked to the shared library (libROMS.{so|dylib})
# because CMAKE_BUILD_WITH_INSTALL_RPATH is set to FALSE so that
# RPATH/RUNPATH are set correctly for both the build tree and
# installed locations of the ROMS executable.

if [[ -z "${ROMS_EXECUTABLE}" || "${ROMS_EXECUTABLE}" == "ON" ]]; then
  if [ $dprint -eq 0 ]; then
    if [[ ! -z "${USE_DEBUG}" && "${USE_DEBUG}" == "on" ]]; then
      cp -pfv ${BUILD_DIR}/romsG .
    elif [[ ! -z "${USE_MPI}" && "${USE_MPI}" == "on" ]]; then
      cp -pfv ${BUILD_DIR}/romsM .
    else
      cp -pfv ${BUILD_DIR}/romsS .
    fi
  fi
fi
