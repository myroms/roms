# svn $Id$
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Copyright (c) 2002-2007 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
# Include file for GNU G95 compiler on Linux
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
               FC := gfortran
           FFLAGS := -frepack-arrays
              CPP := /usr/bin/cpp
         CPPFLAGS := -P -traditional
               LD := $(FC)
          LDFLAGS :=
               AR := ar
          ARFLAGS := -r
               RM := rm -f
           RANLIB := ranlib
             PERL := perl

        MDEPFLAGS := --cpp --fext=f90 --file=-

#
# Library locations, can be overridden by environment variables.
#

       MCT_LIBDIR ?= /usr/local/mct/lib
    NETCDF_INCDIR ?= /usr/local/include
    NETCDF_LIBDIR ?= /usr/local/lib

         CPPFLAGS += -I$(NETCDF_INCDIR)
             LIBS := -L$(NETCDF_LIBDIR) -lnetcdf

ifdef ARPACK
    ARPACK_LIBDIR ?= /usr/local/lib
             LIBS += -L$(ARPACK_LIBDIR) -larpack
endif

ifdef MPI
         CPPFLAGS += -DMPI
 ifdef MPIF90
               FC := mpif90
               LD := $(FC)
 else
             LIBS += -lfmpi -lmpi
 endif
endif

ifdef OpenMP
         CPPFLAGS += -D_OPENMP
           FFLAGS += -fopenmp
endif

ifdef DEBUG
           FFLAGS += -g -fbounds-check -Wall -Wno-unused-variable -Wno-unused-labels
else
           FFLAGS += -O3 -ffast-math
endif

ifdef SWAN_COUPLE
           FFLAGS += -ffixed-form -I/usr/local/mct/include
             LIBS += -L$(MCT_LIBDIR) -lmct -lmpeu
endif

# Turn off bounds checking for function def_var, as "dimension(*)"
# declarations confuse Gnu Fortran 95 bounds-checking code.

def_var.o: FFLAGS += -fno-bounds-check
