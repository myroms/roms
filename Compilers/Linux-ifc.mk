#
# Include file for Intel IFC (version 7.x) compiler on Linux
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
               FC := ifc
           FFLAGS := 
              CPP := /usr/bin/cpp
         CPPFLAGS := -P -traditional -DLINUX
            CLEAN := Bin/cpp_clean
               LD := $(FC)
          LDFLAGS := -Vaxlib
               AR := ar
          ARFLAGS := r
               RM := rm -f
           RANLIB := ranlib
	     PERL := perl

        MDEPFLAGS := --cpp --fext=f90 --file=-

#
# Library locations, can be overridden by environment variables.
#

    NETCDF_INCDIR ?= /opt/intel/compiler70/ia32/include
    NETCDF_LIBDIR ?= /opt/intel/compiler70/ia32/lib

         CPPFLAGS += -I$(NETCDF_INCDIR)
             LIBS := -L$(NETCDF_LIBDIR) -lnetcdf

ifdef ARPACK
    ARPACK_LIBDIR ?= /opt/intelsoft//ARPACK
             LIBS += -L$(ARPACK_LIBDIR) -larpack_LINUX
endif

ifdef MPI
         CPPFLAGS += -DMPI
endif

ifdef OpenMP
         CPPFLAGS += -D_OPENMP
           FFLAGS += -openmp
endif

ifdef DEBUG
           FFLAGS += -g -CA -CB -CS -CU
else
           FFLAGS += -ip -O3
 ifeq ($(CPU),i686)
           FFLAGS += -pc80 -xW
 endif
 ifeq ($(CPU),x86_64)
           FFLAGS += -xW
 endif
endif

       clean_list += ifc* work.pc*
