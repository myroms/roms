#
# Include file for PathScale compiler on Linux
# -----------------------------------------------------------------
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
           FFLAGS := -fno-second-underscore
              CPP := /usr/bin/cpp
         CPPFLAGS := -P -traditional -DLINUX -I/usr/include -DpgiFortran
            CLEAN := Bin/cpp_clean
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

    NETCDF_INCDIR ?= /opt/netcdf-3.6.0-pathscale/include
    NETCDF_LIBDIR ?= /opt/netcdf-3.6.0-pathscale/include

         CPPFLAGS += -I$(NETCDF_INCDIR)
             LIBS := -L$(NETCDF_LIBDIR) -lnetcdf

ifdef ARPACK
    ARPACK_LIBDIR ?= /opt/pathscalesoft/ARPACK
             LIBS += -L$(ARPACK_LIBDIR) -larpack_LINUX
endif

ifdef MPI
         CPPFLAGS += -DMPI
             LIBS +=  -lfmpi-pgi -lmpi-pgi 
endif

ifdef OpenMP
         CPPFLAGS += -D_OPENMP
endif

ifdef DEBUG
           FFLAGS += -g -C
else
           FFLAGS += -Ofast
endif
