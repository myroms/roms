# svn $Id$
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Copyright (c) 2002-2007 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
# Include file for IBM xlf95_r Fortran Compiler on the Macintosh
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
               FC := xlf95_r
           FFLAGS := -qsuffix=f=f90 -qmaxmem=-1 -qextname
              CPP := /usr/bin/cpp
         CPPFLAGS := -P -DMAC -traditional
               LD = $(FC)
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

#    NETCDF_INCDIR ?= /usr/local/include
#    NETCDF_LIBDIR ?= /usr/local/lib
    NETCDF_INCDIR ?= ${HOME}/include
    NETCDF_LIBDIR ?= ${HOME}/lib
         CPPFLAGS += -I$(NETCDF_INCDIR)
             LIBS := -L$(NETCDF_LIBDIR) -lnetcdf

ifdef ARPACK
 ifdef MPI
   PARPACK_LIBDIR ?= /usr/local/lib
             LIBS += -L$(PARPACK_LIBDIR) -lparpack
 endif
    ARPACK_LIBDIR ?= /usr/local/lib
             LIBS += -L$(ARPACK_LIBDIR) -larpack
endif

ifdef MPI
         CPPFLAGS += -DMPI
               FC := mpxlf95_r
endif

ifdef OpenMP
         CPPFLAGS += -D_OPENMP
           FFLAGS += -qsmp=omp
endif

ifdef DEBUG
           FFLAGS += -g -qfullpath
else
           FFLAGS += -O3 -qstrict
endif
