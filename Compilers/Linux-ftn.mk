#
# Include file for CRAY FTN cross-compiler with Linux
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
               FC := ftn
           FFLAGS := -e I -e m
              CPP := /usr/bin/cpp
         CPPFLAGS := -P -traditional -DCRAYX1 -DCRAY
            CLEAN := Bin/cpp_clean
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

    NETCDF_INCDIR ?= $(HOME)/include
    NETCDF_LIBDIR ?= $(HOME)/lib

         CPPFLAGS += -I$(NETCDF_INCDIR)
             LIBS := -L$(NETCDF_LIBDIR) -lnetcdf

ifdef ARPACK
    ARPACK_LIBDIR ?= /usr/local/lib
             LIBS += -L$(ARPACK_LIBDIR) -larpack
endif

ifdef MPI
         CPPFLAGS += -DMPI
endif

ifdef OpenMP
         CPPFLAGS += -D_OPENMP
endif

ifdef DEBUG
           FFLAGS += -G 0
else
           FFLAGS += -O 3,aggress
endif

# Cray specials

lmd_skpp.o: FFLAGS += -O inlinefrom=lmd_wscale.f90
lmd_bkpp.o: FFLAGS += -O inlinefrom=lmd_wscale.f90

lmd_skpp.o: lmd_wscale.f90
lmd_bkpp.o: lmd_wscale.f90
