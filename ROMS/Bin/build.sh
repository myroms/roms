#!/bin/csh -f
# 
# svn $Id$
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::: John Wilkin :::
# Copyright (c) 2002-2008 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::: Hernan G. Arango :::
#                                                                       :::
# ROMS/TOMS Compiling Script                                            :::
#                                                                       :::
# Script to compile an user application where the application-specific  :::
# files are kept separate from the ROMS source code.                    :::
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
#    ./build.sh [options]                                               :::
#                                                                       :::
# Options:                                                              :::
#                                                                       :::
#    -j [N]      Compile in parallel using N CPUs                       :::
#                  omit argument for all available CPUs                 :::
#    -noclean    Do not clean already compiled objects                  :::
#                                                                       :::
# Notice that sometimes the parallel compilation fail to find MPI       :::
# include file "mpif.h".                                                :::
#                                                                       :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

set parallel = 0
set clean = 1

while ( ($#argv) > 0 )
  switch ($1)
    case "-noclean"
      shift
      set clean = 0
    breaksw

    case "-j"
      shift
      set parallel = 1
      if (`echo $1 | grep -P '^\d+$'` != "" ) then
        set NCPUS = "-j $1"
        shift
      else
        set NCPUS = "-j"
      endif
    breaksw

    case "-*":
      echo ""
      echo "$0 : Unknonw option [ $1 ]"
      echo ""
      echo "Available Options:"
      echo ""
      echo "-j [N]      Compile in parallel using N CPUs"
      echo "              omit argument for all avaliable CPUs"
      echo "-noclean    Do not clean already compiled objects"
      echo ""
      exit 1
    breaksw

  endsw
end

# Set the CPP option defining the particular application. This will
# determine the name of the ".h" header file with the application 
# CPP definitions.

setenv ROMS_APPLICATION     CBLAST

# Set a local environmental variable to define the path to the directories
# where all this project's files are kept.

setenv MY_ROOT_DIR          /home/arango/ocean/toms/repository
setenv MY_PROJECT_DIR       ${MY_ROOT_DIR}/Projects/cblast

# The path to the user's local current ROMS source code. 
#
# If using svn locally, this would be the user's Working Copy Path (WCPATH). 
# Note that one advantage of maintaining your source code locally with svn 
# is that when working simultaneously on multiple machines (e.g. a local 
# workstation, a local cluster and a remote supercomputer) you can checkout 
# the latest release and always get an up-to-date customized source on each 
# machine. This script is designed to more easily allow for differing paths 
# to the code and inputs on differing machines. 

setenv MY_ROMS_SRC          ${MY_ROOT_DIR}/branches/arango

# Set tunable CPP options.
#
# Sometimes it is desirable to activate one or more CPP options to run
# different variants of the same application without modifying its header
# file. If this is the case, specify each options here using the -D syntax.
# Notice also that you need to use shell's quoting syntax to enclose the
# definition.  Both single or double quotes works. For example, to write
# time-averaged fields set:
#
#    setenv MY_CPP_FLAGS "-DAVERAGES"

#setenv MY_CPP_FLAGS "-D"

# Other user defined environmental variables. See the ROMS makefile for
# details on other options the user might want to set here. 

 setenv USE_MPI             on
 setenv USE_MPIF90          on
 setenv FORT                pgi

#setenv USE_OpenMP          on

#setenv USE_DEBUG           on
 setenv USE_LARGE           on
#setenv USE_NETCDF4         on

# The rest of this script sets the path to the users header file and
# analytical source files, if any. See the templates in User/Functionals.

 setenv MY_HEADER_DIR       ${MY_PROJECT_DIR}/Forward

#setenv MY_ANALYTICAL_DIR   ${MY_PROJECT_DIR}/Functionals

# Put the binary to execute in the following directory.

 setenv BINDIR              ${MY_PROJECT_DIR}/Forward

# Put the f90 files in a project specific Build directory to avoid conflict
# with other projects. 

 setenv SCRATCH_DIR         ${MY_PROJECT_DIR}/Forward/Build

# Go to the users source directory to compile. The options set above will
# pick up the application-specific code from the appropriate place.

 cd ${MY_ROMS_SRC}

# Stop if activating both MPI and OpenMP at the same time.

if ( ${?USE_MPI} & ${?USE_OpenMP} ) then
  echo "You cannot activate USE_MPI and USE_OpenMP at the same time!"
  exit 1
endif

# Remove build directory. 

if ( $clean == 1 ) then
  make clean
endif

# Compile (the binary will go to BINDIR set above).  

if ( $parallel == 1 ) then
  make $NCPUS
else
  make
endif
