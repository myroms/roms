# svn $Id$
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Copyright (c) 2002-2007 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
# Include file for PGI Fortran compiler on Linux
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
               FC := pgf90
           FFLAGS :=
              CPP := /usr/bin/cpp
         CPPFLAGS := -P -traditional
               LD := $(FC)
          LDFLAGS := 
               AR := ar
          ARFLAGS := r
            MKDIR := mkdir -p
               RM := rm -f
           RANLIB := ranlib
	     PERL := perl
             TEST := test

        MDEPFLAGS := --cpp --fext=f90 --file=- --objdir=$(SCRATCH_DIR)

#
# Library locations, can be overridden by environment variables.
#

       MCT_LIBDIR ?= /usr/local/mct/lib
    NETCDF_INCDIR ?= /opt/pgisoft/netcdf/include
    NETCDF_LIBDIR ?= /opt/pgisoft/netcdf/lib

         CPPFLAGS += -I$(NETCDF_INCDIR)
             LIBS := -L$(NETCDF_LIBDIR) -lnetcdf
#            LIBS := -L$(NETCDF_LIBDIR) -lnetcdf -lgmalloc

ifdef ARPACK
 ifdef MPI
#  PARPACK_LIBDIR ?= /opt/pgisoft/PARPACK
   PARPACK_LIBDIR ?= $(ROMSHOME)/Lib
  ifdef MPIF90
   ifdef DEBUG
             LIBS += -L$(PARPACK_LIBDIR) -lparpack_dbg_daggoo
   else
             LIBS += -L$(PARPACK_LIBDIR) -lparpack_daggoo
   endif
  else
   ifdef DEBUG
             LIBS += -L$(PARPACK_LIBDIR) -lparpack_dbg_moby
   else
             LIBS += -L$(PARPACK_LIBDIR) -lparpack_moby
#            LIBS += -L$(PARPACK_LIBDIR) -lparpack
   endif
  endif
 endif
#   ARPACK_LIBDIR ?= /opt/pgisoft/PARPACK
    ARPACK_LIBDIR ?= $(ROMSHOME)/Lib
 ifdef MPIF90
  ifdef DEBUG
             LIBS += -L$(ARPACK_LIBDIR) -larpack_dbg_daggoo
  else
             LIBS += -L$(ARPACK_LIBDIR) -larpack_daggoo
  endif
 else
  ifdef DEBUG
             LIBS += -L$(ARPACK_LIBDIR) -larpack_dbg_moby
  else
             LIBS += -L$(ARPACK_LIBDIR) -larpack_moby
#            LIBS += -L$(ARPACK_LIBDIR) -larpack
  endif
 endif
endif

ifdef MPI
         CPPFLAGS += -DMPI
 ifdef MPIF90
               FC := mpif90
               LD := $(FC)
 else
             LIBS += -lfmpi-pgi -lmpi-pgi 
 endif
endif

ifdef OpenMP
         CPPFLAGS += -D_OPENMP
endif

ifdef DEBUG
#          FFLAGS += -g -C -Mchkstk -Mchkfpstk
           FFLAGS += -g -C
#          FFLAGS += -g
else
           FFLAGS += -u -Bstatic -fastsse -Mipa=fast -tp k8-64
endif

ifdef SWAN_COUPLE
           FFLAGS += -Mfixed -I/usr/local/mct/include
             LIBS += -L$(MCT_LIBDIR) -lmct -lmpeu
endif

       clean_list += ifc* work.pc*

#
# Set free form format in source files to allow long string for
# local directory and compilation flags inside the code.
#

$(SCRATCH_DIR)/mod_ncparam.o: FFLAGS += -Mfree
$(SCRATCH_DIR)/mod_strings.o: FFLAGS += -Mfree
