#
# Include file for IBM xlf95_r Fortran Compiler
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
               FC := xlf95_r
           FFLAGS := -qsuffix=f=f90 -qmaxmem=-1 -qarch=pwr4 -qtune=pwr4
              CPP := /usr/lib/cpp
         CPPFLAGS := -P -DAIX
            CLEAN := Bin/cpp_clean
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

ifdef LARGE
           FFLAGS += -q64
          ARFLAGS += -X 64
          LDFLAGS += -bmaxdata:0x200000000
    NETCDF_INCDIR ?= /usr/local/pkg/netcdf/netcdf-3.5.0_64/include
    NETCDF_LIBDIR ?= /usr/local/pkg/netcdf/netcdf-3.5.0_64/lib
else
          LDFLAGS += -bmaxdata:0x70000000
    NETCDF_INCDIR ?= /usr/local/include
    NETCDF_LIBDIR ?= /usr/local/lib
endif
         CPPFLAGS += -I$(NETCDF_INCDIR)
             LIBS := -L$(NETCDF_LIBDIR) -lnetcdf

ifdef ARPACK
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
