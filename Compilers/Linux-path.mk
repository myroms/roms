# svn $Id$
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Copyright (c) 2002-2007 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
# Include file for PathScale compiler on Linux
# -------------------------------------------------------------------------
#
# ARPACK_LIBDIR  ARPACK libary directory
# FC             Name of the fortran compiler to use
# FFLAGS         Flags to the fortran compiler
# CPP            Name of the C-preprocessor
# CPPFLAGS       Flags to the C-preprocessor
# CLEAN          Name of cleaning executable after C-preprocessing
# NETCDF_INCDIR  NetCDF include directory
# NETCDF_LIBDIR  NetCDF libary directory
# LD             Program to load the objects into an executable
# LDFLAGS        Flags to the loader
# RANLIB         Name of ranlib command
# MDEPFLAGS      Flags for sfmakedepend  (-s if you keep .f files)
#
# First the defaults
#
               FC := pathf90
           FFLAGS := -fno-second-underscore -I/usr/include
              CPP := /usr/bin/cpp
         CPPFLAGS := -P -traditional -I/usr/include -DpgiFortran
               LD := $(FC)
          LDFLAGS := 
               AR := ar
          ARFLAGS := r
               RM := rm -f
           RANLIB := ranlib
	     PERL := perl

        MDEPFLAGS := --cpp --fext=f90 --file=-

#
# Library locations, can be overridden by environment variables.
#

       MCT_LIBDIR ?= /usr/local/mct/lib
    NETCDF_INCDIR ?= /opt/pathscalesoft/netcdf-3.6.0-p1/include
    NETCDF_LIBDIR ?= /opt/pathscalesoft/netcdf-3.6.0-p1/lib

         CPPFLAGS += -I$(NETCDF_INCDIR)
             LIBS := -L$(NETCDF_LIBDIR) -lnetcdf

ifdef ARPACK
 ifdef MPI
   PARPACK_LIBDIR ?= /opt/pathscalesoft/PARPACK
             LIBS += -L$(PARPACK_LIBDIR) -lparpack
 endif
    ARPACK_LIBDIR ?= /opt/pathscalesoft/PARPACK
             LIBS += -L$(ARPACK_LIBDIR) -larpack
endif

ifdef MPI
         CPPFLAGS += -DMPI
 ifdef MPIF90
               FC := /opt/pathscalesoft/mpich/bin/mpif90
               LD := $(FC)
 else
             LIBS += -lfmpi-pgi -lmpi-pgi 
 endif
endif

ifdef OpenMP
         CPPFLAGS += -D_OPENMP
endif

ifdef DEBUG
           FFLAGS += -g -C
else
           FFLAGS += -Ofast
endif

ifdef SWAN_COUPLE
           FFLAGS += -fixedform -I/usr/local/mct/include
             LIBS += -L$(MCT_LIBDIR) -lmct -lmpeu
endif
