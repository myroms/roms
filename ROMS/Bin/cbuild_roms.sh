#!/bin/bash
#
# git $Id$
# svn $Id: cbuild_roms.sh 1051 2020-12-04 23:09:05Z arango $
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Copyright (c) 2002-2020 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#:::::::::::::::::::::::::::::::::::::::::::::::::::::: David Robertson :::
#                                                                       :::
# ROMS ecbuild (CMake) Compiling BASH Script                            :::
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
#    -j [N]      Compile in parallel using N CPUs                       :::
#                  omit argument for all available CPUs                 :::
#                                                                       :::
#    -p macro    Prints any Makefile macro value. For example,          :::
#                                                                       :::
#                  cbuild_roms.sh -p MY_CPP_FLAGS                       :::
#                                                                       :::
#    -noclean    Do not clean already compiled objects                  :::
#                                                                       :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

export which_MPI=openmpi                       # default, overwritten below

parallel=0
clean=1
dprint=0

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
      dprint=1
      debug="$1"
      shift
      ;;

    -noclean )
      shift
      clean=0
      ;;

    * )
      echo ""
      echo "$0 : Unknown option [ $1 ]"
      echo ""
      echo "Available Options:"
      echo ""
      echo "-j [N]      Compile in parallel using N CPUs"
      echo "              omit argument for all avaliable CPUs"
      echo ""
      echo "-p macro    Prints any Makefile macro value"
      echo "              For example:  cbuild_roms.sh -p FFLAGS"
      echo ""
      echo "-noclean    Do not clean already compiled objects"
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
# where all this project's files are kept.

export        MY_ROOT_DIR=${HOME}/ocean/repository
export     MY_PROJECT_DIR=${PWD}

# The path to the user's local current ROMS source code.

 export       MY_ROMS_SRC=${MY_ROOT_DIR}/trunk

# Which type(s) of libraries would you like?
# NOTE: If you choose both and also choose to build the ROMS executable,
#       it will be linked to the static version of the library.
#
# Valid options are SHARED, STATIC, and BOTH.

 export           LIBTYPE=BOTH

# Set the path to the users header (REQUIRED).

 export     MY_HEADER_DIR=${MY_PROJECT_DIR}

# If you have custom analytical functions to include, enter the path here.

 export MY_ANALYTICAL_DIR=${MY_PROJECT_DIR}

# Does your application include a biological or sediment model?
# The logic here needs adjusting. For now, leave commented unless you want
# a biological model

#export           BIOLOGY=ON
#export          SEDIMENT=ON

# Do you want to build the ROMS executable?
# Valid values are: ON (build the executable) and OFF (do NOT build the
# executable). If you comment this out the executable WILL be built.

 export   ROMS_EXECUTABLE=ON

# Set path of the directory containing my_build_paths.sh.
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

#--------------------------------------------------------------------------
# Compilation options and paths.
#--------------------------------------------------------------------------

 export           USE_MPI=on                 # distributed-memory parallelism
 export        USE_MPIF90=on                 # compile with mpif90 script
#export         which_MPI=mpich              # compile with MPICH library
#export         which_MPI=mpich2             # compile with MPICH2 library
#export         which_MPI=mvapich2           # compile with MVAPICH2 library
 export         which_MPI=openmpi            # compile with OpenMPI library

#export              FORT=ifort
 export              FORT=gfortran
#export              FORT=pgi

 export                FC=mpifort            # ENV var that CMake looks for

#--------------------------------------------------------------------------
# If applicable, use my specified library paths.
#--------------------------------------------------------------------------

 export USE_MY_LIBS=no            # use system default library paths
#export USE_MY_LIBS=yes           # use my customized library paths

MY_PATHS=${COMPILERS}/my_build_paths.sh

if [ "${USE_MY_LIBS}" = "yes" ]; then
  source ${MY_PATHS} ${MY_PATHS}
fi

# Put the CMake files in a project specific Build directory to avoid conflict
# with other projects.

if [ -n "${USE_DEBUG:+1}" ]; then
  export      SCRATCH_DIR=${MY_PROJECT_DIR}/CBuild_romsG
else
  export      SCRATCH_DIR=${MY_PROJECT_DIR}/CBuild_roms
fi

# Create the build directory specified above and change into it.

if [ -d ${SCRATCH_DIR} ]; then
  if [ $clean -eq 1 ]; then
    rm -rf ${SCRATCH_DIR}
    mkdir ${SCRATCH_DIR}
    cd ${SCRATCH_DIR}
  else
    cd ${SCRATCH_DIR}
  fi
else
  if [ $clean -eq 1 ]; then
    mkdir ${SCRATCH_DIR}
    cd ${SCRATCH_DIR}
  else
    echo "-noclean option activated when the build directory didn't exist"
    echo "creating the directory and disabling -noclean"
    clean=1
    mkdir ${SCRATCH_DIR}
    cd ${SCRATCH_DIR}
  fi
fi

#--------------------------------------------------------------------------
# Configure.
#--------------------------------------------------------------------------

# Construct the ecbuild command.

if [ ! -z "${LIBTYPE}" ]; then
  ltype="-DLIBTYPE=${LIBTYPE}"
else
  ltype=""
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

if [[ ! -z "${USE_MPI}" && "${USE_MPI}" == "on" ]]; then
  mpi="-DMPI=ON"
else
  mpi=""
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

if [[ ! -z "${BIOLOGY}" && "${BIOLOGY}" == "ON" ]]; then
 bio="-DBIOLOGY=ON"
else
 bio=""
fi

if [[ ! -z "${SEDIMENT}" && "${SEDIMENT}" == "ON" ]]; then
 sediment="-DSEDIMENT=ON"
else
 sediment=""
fi

#--------------------------------------------------------------------------
# Run the ecbuild command.
#--------------------------------------------------------------------------

my_hdir="-DMY_HEADER_DIR=${MY_HEADER_DIR}"

if [ $dprint -eq 0 ]; then
  ecbuild -DAPP=${ROMS_APPLICATION} \
                ${my_hdir} \
                ${ltype} \
                ${extra_flags} \
                ${parpack_ldir} \
                ${arpack_ldir} \
                ${mpi} \
                ${roms_exec} \
                ${bio} \
                ${sediment} \
                ${MY_ROMS_SRC}
fi

if [ $? -ne 0 ]; then
  echo "ecbuild did not complete successfully"
  exit 1
fi

#--------------------------------------------------------------------------
# Compile
#--------------------------------------------------------------------------

if [ $dprint -eq 1 ]; then
  echo $debug:"${!debug}"
else
  if [ $parallel -eq 1 ]; then
    make $NCPUS
  else
    make
  fi
fi

cd ${MY_PROJECT_DIR}

# This needs more logic and a way to chose the right executable name.

if [ $dprint -eq 0 ]; then
  if [ ! -h romsM ]; then
    ln -s ${SCRATCH_DIR}/romsM
  fi
fi

# The alternative to this is to set your LD_LIBRARY_PATH appropriately
# outside of this script.

if [ $dprint -eq 0 ]; then
  if [[ "${LIBTYPE}" == "SHARED" || "${LIBTYPE}" == "BOTH" ]]; then
    if [ ! -h libROMS.so ]; then
      ln -s ${SCRATCH_DIR}/libROMS.so
    fi
  fi
fi
