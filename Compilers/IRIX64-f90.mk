# svn $Id$
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Copyright (c) 2002-2007 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
# Include file for IRIX F90 compiler on SGI
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
# LDR            Program to load the objects into an executable
# LDFLAGS        Flags to the loader
# RANLIB         Name of ranlib command
# MDEPFLAGS      Flags for sfmakedepend  (-s if you keep .f files)
#
# First the defaults
#
               FC := f90
           FFLAGS := -mips4 -u -TENV:X=3
              CPP := /usr/lib32/cmplrs/cpp
         CPPFLAGS := -P -DSGI
               LD := $(FC)
          LDFLAGS :=
               AR := ar
          ARFLAGS := -r
               RM := rm -f
           RANLIB := touch
	     PERL := perl

        MDEPFLAGS := --cpp --fext=f90 --file=-

#
# Library locations, can be overridden by environment variables.
#

ifdef LARGE
           FFLAGS += -64
       MCT_LIBDIR ?= /usr/local/mct/lib
    NETCDF_INCDIR ?= $(HOME)/netcdf/include
    NETCDF_LIBDIR ?= $(HOME)/netcdf/lib64
else
           FFLAGS += -n32
       MCT_LIBDIR ?= /usr/local/mct/lib
    NETCDF_INCDIR ?= /usr/local/include
    NETCDF_LIBDIR ?= /usr/local/lib
endif
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
           FFLAGS += -LANG:recursive=on
          LDFLAGS += -mp -mp_schedtype=simple
             LIBS += -lmpi
endif

ifdef OpenMP
         CPPFLAGS += -D_OPENMP
           FFLAGS += -mp -MP:open_mp=ON
endif

ifdef DEBUG
           FFLAGS += -g -C
else
           FFLAGS += -O3 -OPT:Olimit=4038
endif

ifdef SWAN_COUPLE
           FFLAGS += -fixed -I/usr/local/mct/include
             LIBS += -L$(MCT_LIBDIR) -lmct -lmpeu
endif
