# git $Id$
# svn $Id: Darwin-gfortran.mk 1202 2023-10-24 15:36:07Z arango $
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Copyright (c) 2002-2023 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.md                                                 :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
# Include file for GNU Fortran compiler on Darwin
# -------------------------------------------------------------------------
#
# ARPACK_LIBDIR  ARPACK libary directory
# FC             Name of the fortran compiler to use
# FFLAGS         Flags to the fortran compiler
# CPP            Name of the C-preprocessor
# CPPFLAGS       Flags to the C-preprocessor
# HDF5_INCDIR    HDF5 include directory
# HDF5_LIBDIR    HDF5 library directory
# HDF5_LIBS      HDF5 library switches
# LIBS           Required libraries during linking
# ROMS_LIB       Directory and name for ROMS library
# NF_CONFIG      NetCDF Fortran configuration script
# NETCDF_INCDIR  NetCDF include directory
# NETCDF_LIBDIR  NetCDF library directory
# NETCDF_LIBS    NetCDF library switches
# PIO_INCDIR     Parallel-IO (PIO) from SCORPIO library include directory
# PIO_LIBDIR     Parallel-IO (PIO) from SCORPIO libary directory
# PIO_LIBS       Parallel-IO (PIO) from SCORPIO library switches
# PNETCDF_INCDIR PNetCDF include directory
# PNETCDF_LIBDIR PNetCDF libary directory
# PNETCDF_LIBS   PNetCDF library switches

# LD             Program to load the objects into an executable or shared library
# LDFLAGS        Flags to the loader
# RANLIB         Name of ranlib command
# MDEPFLAGS      Flags for sfmakedepend  (-s if you keep .f files)
#
# First the defaults
#
               FC := gfortran
           FFLAGS := -frepack-arrays
       FIXEDFLAGS := -ffixed-form
        FREEFLAGS := -ffree-form -ffree-line-length-none
              CPP := /usr/bin/cpp
         CPPFLAGS := -P -traditional-cpp -w          # -w turns off warnings
           INCDIR := /usr/include /usr/local/bin
            SLIBS := -L/usr/local/lib -L/usr/lib
            ULIBS :=
             LIBS :=
         ROMS_LIB := -L$(BUILD_DIR) -lROMS
       MOD_SUFFIX := mod
               LD := $(FC)
          LDFLAGS :=
               AR := ar
          ARFLAGS := -r
            MKDIR := mkdir -p
               CP := cp -p -v
               RM := rm -f
           RANLIB := ranlib -no_warning_for_no_symbols
             PERL := perl
             TEST := test
      ST_LIB_NAME := libROMS.a
      SH_LIB_NAME := libROMS.dylib

#--------------------------------------------------------------------------
# Checking for minimum version and versions requiring extra compile flags
#--------------------------------------------------------------------------

# Because of the recursive derived types in Utility/yaml_parser.F, gfortran
# version 7 or greater is required.

          MIN_VER := $(shell expr `$(FC) -dumpversion | cut -f1 -d.` \<= 6)
ifeq "$(MIN_VER)" "1"
      $(error gfortran version 7 or greater is required)
endif

# Starting with gfortran version 10, GNU takes a much stricter approach to argument
# type and size matching. This conflicts with ROMS high level MPI calls so this flag
# is added for compilation to complete successfully.

       STRICT_VER := $(shell expr `$(FC) -dumpversion | cut -f1 -d.` \>= 10)
ifeq "$(STRICT_VER)" "1"
           FFLAGS += -fallow-argument-mismatch       # needed for gfortran 10 and higher
endif

#--------------------------------------------------------------------------
# Compiling flags for ROMS Applications.
#--------------------------------------------------------------------------

ifdef USE_ROMS
 ifdef USE_DEBUG
           FFLAGS += -g
           FFLAGS += -fbounds-check
           FFLAGS += -fbacktrace
           FFLAGS += -fcheck=all
#          FFLAGS += -fsanitize=address -fsanitize=undefined
           FFLAGS += -finit-real=nan -ffpe-trap=invalid,zero,overflow
           FFLAGS += -fmax-stack-var-size=64000000
 else
           FFLAGS += -O3
#          FFLAGS += -ffast-math
 endif
 ifdef SHARED
          LDFLAGS += -Wl,-rpath,$(BUILD_DIR)

       SH_LDFLAGS += -dynamiclib
       SH_LDFLAGS += -Wl,-headerpad_max_install_names
       SH_LDFLAGS += -Wl,-undefined,dynamic_lookup

       SH_LDFLAGS += -install_name @rpath/$(SH_LIB_NAME)
 endif

        MDEPFLAGS := --cpp --fext=f90 --file=- --objdir=$(BUILD_DIR)
endif

#--------------------------------------------------------------------------
# Compiling flags for CICE Applications.
#--------------------------------------------------------------------------

ifdef CICE_APPLICATION
          CPPDEFS := -DLINUS $(MY_CPP_FLAGS)
 ifdef USE_DEBUG
           FFLAGS += -g
#          FFLAGS += -O3
           FFLAGS += -fbounds-check
#          FFLAGS += -fcheck=all
#          FFLAGS += -fsanitize=address -fsanitize=undefined
 else
           FFLAGS := -O3 -w
 endif
endif

#--------------------------------------------------------------------------
# Coupled models libraries.
#--------------------------------------------------------------------------

ifdef USE_COAMPS
             LIBS += $(COAMPS_LIB_DIR)/libashare.a
             LIBS += $(COAMPS_LIB_DIR)/libcoamps.a
             LIBS += $(COAMPS_LIB_DIR)/libaa.a
             LIBS += $(COAMPS_LIB_DIR)/libam.a
             LIBS += $(COAMPS_LIB_DIR)/libfishpak.a
             LIBS += $(COAMPS_LIB_DIR)/libfnoc.a
             LIBS += $(COAMPS_LIB_DIR)/libtracer.a
#            LIBS += $(COAMPS_LIB_DIR)/libnl_beq.a
endif

ifdef CICE_APPLICATION
            SLIBS += $(SLIBS) $(LIBS)
endif

ifdef USE_WRF
             LIBS += $(WRF_LIB_DIR)/module_wrf_top.o
             LIBS += $(WRF_LIB_DIR)/libwrflib.a
             LIBS += $(WRF_LIB_DIR)/libfftpack.a
             LIBS += $(WRF_LIB_DIR)/libio_grib1.a
             LIBS += $(WRF_LIB_DIR)/libio_grib_share.a
             LIBS += $(WRF_LIB_DIR)/libwrfio_int.a
             LIBS += $(WRF_LIB_DIR)/libesmf_time.a
             LIBS += $(WRF_LIB_DIR)/librsl_lite.a
             LIBS += $(WRF_LIB_DIR)/module_internal_header_util.o
             LIBS += $(WRF_LIB_DIR)/pack_utils.o
             LIBS += $(WRF_LIB_DIR)/libwrfio_nf.a
endif

#--------------------------------------------------------------------------
# Library locations, can be overridden by environment variables.
#--------------------------------------------------------------------------


ifdef USE_PIO
       PIO_INCDIR ?= /opt/gfortransoft/openmpi/pio/include
       PIO_LIBDIR ?= /opt/gfortransoft/openmpi/pio/lib
           FFLAGS += -I$(PIO_INCDIR)
             LIBS += -L$(PIO_LIBDIR) -lpiof -lpioc

   PNETCDF_INCDIR ?= /opt/gfortransoft/openmpi/pnetcdf/include
   PNETCDF_LIBDIR ?= /opt/gfortransoft/openmpi/pnetcdf/lib
           FFLAGS += -I$(PNETCDF_INCDIR)
             LIBS += -L$(PNETCDF_LIBDIR) -lpnetcdf
endif

ifdef USE_SCORPIO
       PIO_INCDIR ?= /opt/gfortransoft/openmpi/scorpio/include
       PIO_LIBDIR ?= /opt/gfortransoft/openmpi/scorpio/lib
           FFLAGS += -I$(PIO_INCDIR)
             LIBS += -L$(PIO_LIBDIR) -lpiof -lpioc

   PNETCDF_INCDIR ?= /opt/intelsoft/openmpi/pnetcdf/include
   PNETCDF_LIBDIR ?= /opt/intelsoft/openmpi/pnetcdf/lib
           FFLAGS += -I$(PNETCDF_INCDIR)
             LIBS += -L$(PNETCDF_LIBDIR) -lpnetcdf
endif

ifdef USE_NETCDF4
        NF_CONFIG ?= nf-config
    NETCDF_INCDIR ?= $(shell $(NF_CONFIG) --prefix)/include
             LIBS += $(shell $(NF_CONFIG) --flibs)
           INCDIR += $(NETCDF_INCDIR) $(INCDIR)
else
    NETCDF_INCDIR ?= /opt/gfortransoft/serial/netcdf3/include
    NETCDF_LIBDIR ?= /opt/gfortransoft/serial/netcdf3/lib
      NETCDF_LIBS ?= -lnetcdf
             LIBS += -L$(NETCDF_LIBDIR) $(NETCDF_LIBS)
           INCDIR += $(NETCDF_INCDIR) $(INCDIR)
endif

ifdef USE_HDF5
      HDF5_INCDIR ?= /opt/gfortransoft/serial/hdf5/include
      HDF5_LIBDIR ?= /opt/gfortransoft/serial/hdf5/lib
        HDF5_LIBS ?= -lhdf5_fortran -lhdf5hl_fortran -lhdf5 -lz
             LIBS += -L$(HDF5_LIBDIR) $(HDF5_LIBS)
           INCDIR += $(HDF5_INCDIR)
endif

ifdef USE_ARPACK
 ifdef USE_MPI
   PARPACK_LIBDIR ?= /opt/gfortransoft/PARPACK
  ifdef USE_DEBUG
             LIBS += -L$(PARPACK_LIBDIR) -lparpack-debug
  else
             LIBS += -L$(PARPACK_LIBDIR) -lparpack
  endif
 endif
    ARPACK_LIBDIR ?= /opt/gfortransoft/ARPACK
 ifdef USE_DEBUG
             LIBS += -L$(ARPACK_LIBDIR) -larpack-debug
 else
             LIBS += -L$(ARPACK_LIBDIR) -larpack
 endif
endif

ifdef USE_MPI
         CPPFLAGS += -DMPI
 ifdef USE_MPIF90
               FC := mpif90
 else
             LIBS += -lfmpi -lmpi
 endif
endif

ifdef USE_OpenMP
         CPPFLAGS += -D_OPENMP
           FFLAGS += -fopenmp -static-libgcc
#            LIBS += -lgomp
endif

ifdef USE_MCT
       MCT_INCDIR ?= /usr/local/mct/include
       MCT_LIBDIR ?= /usr/local/mct/lib
           FFLAGS += -I$(MCT_INCDIR)
             LIBS += -L$(MCT_LIBDIR) -lmct -lmpeu
           INCDIR += $(MCT_INCDIR) $(INCDIR)
endif

ifdef USE_ESMF
                     include $(ESMFMKFILE)
          ESMF_OS ?= $(OS)
      ESMF_SUBDIR := $(ESMF_OS).$(ESMF_COMPILER).$(ESMF_ABI).$(ESMF_COMM).$(ESMF_SITE)
      ESMF_MK_DIR ?= $(ESMF_DIR)/lib/lib$(ESMF_BOPT)/$(ESMF_SUBDIR)
           FFLAGS += $(ESMF_F90COMPILEPATHS)
             LIBS += $(ESMF_F90LINKPATHS) $(ESMF_F90ESMFLINKLIBS)
endif

# Use full path of compiler.

               FC := $(shell which ${FC})
               LD := $(FC)

#--------------------------------------------------------------------------
# ROMS specific rules.
#--------------------------------------------------------------------------

# Turn off bounds checking for function def_var, as "dimension(*)"
# declarations confuse Gnu Fortran 95 bounds-checking code.

ifdef USE_ROMS
 $(BUILD_DIR)/def_var.o: FFLAGS += -fno-bounds-check
endif

# Allow integer overflow in ran_state.F.  This is not allowed
# during -O3 optimization. This option should be applied only for
# Gfortran versions >= 4.2.

ifdef USE_ROMS
 FC_TEST := $(findstring $(shell ${FC} --version | head -1 | cut -d " " -f 5 | \
                           cut -d "." -f 1-2),4.0 4.1)

 ifeq "${FC_TEST}" ""
  $(BUILD_DIR)/ran_state.o: FFLAGS += -fno-strict-overflow
 endif
endif

# Set free form format in some ROMS source files to allow long string for
# local directory and compilation flags inside the code.

ifdef USE_ROMS
 $(BUILD_DIR)/mod_ncparam.o: FFLAGS += $(FREEFLAGS)
 $(BUILD_DIR)/mod_strings.o: FFLAGS += $(FREEFLAGS)
 $(BUILD_DIR)/analytical.o: FFLAGS += $(FREEFLAGS)
 $(BUILD_DIR)/biology.o: FFLAGS += $(FREEFLAGS)

 ifdef USE_ADJOINT
  $(BUILD_DIR)/ad_biology.o: FFLAGS += $(FREEFLAGS)
 endif
 ifdef USE_REPRESENTER
  $(BUILD_DIR)/rp_biology.o: FFLAGS += $(FREEFLAGS)
 endif
 ifdef USE_TANGENT
  $(BUILD_DIR)/tl_biology.o: FFLAGS += $(FREEFLAGS)
 endif
endif

#--------------------------------------------------------------------------
# Model coupling specific rules.
#--------------------------------------------------------------------------

# Add COAMPS library directory to include path of ESMF coupling files.

ifdef USE_COAMPS
 $(BUILD_DIR)/esmf_atm.o: FFLAGS += -I$(COAMPS_LIB_DIR)
 $(BUILD_DIR)/esmf_esm.o: FFLAGS += -I$(COAMPS_LIB_DIR)
endif

# Add WRF library directory to include path of ESMF coupling files.

ifdef USE_WRF
 $(BUILD_DIR)/esmf_atm.o: FFLAGS += -I$(WRF_LIB_DIR)
endif

# Supress free format in SWAN source files since there are comments
# beyond column 72.

ifdef USE_SWAN
 $(BUILD_DIR)/ocpcre.o: FFLAGS += $(FIXEDFLAGS)
 $(BUILD_DIR)/ocpids.o: FFLAGS += $(FIXEDFLAGS)
 $(BUILD_DIR)/ocpmix.o: FFLAGS += $(FIXEDFLAGS)
 $(BUILD_DIR)/swancom1.o: FFLAGS += $(FIXEDFLAGS)
 $(BUILD_DIR)/swancom2.o: FFLAGS += $(FIXEDFLAGS)
 $(BUILD_DIR)/swancom3.o: FFLAGS += $(FIXEDFLAGS)
 $(BUILD_DIR)/swancom4.o: FFLAGS += $(FIXEDFLAGS)
 $(BUILD_DIR)/swancom5.o: FFLAGS += $(FIXEDFLAGS)
 $(BUILD_DIR)/swanmain.o: FFLAGS += $(FIXEDFLAGS)
 $(BUILD_DIR)/swanout1.o: FFLAGS += $(FIXEDFLAGS)
 $(BUILD_DIR)/swanout2.o: FFLAGS += $(FIXEDFLAGS)
 $(BUILD_DIR)/swanparll.o: FFLAGS += $(FIXEDFLAGS)
 $(BUILD_DIR)/swanpre1.o: FFLAGS += $(FIXEDFLAGS)
 $(BUILD_DIR)/swanpre2.o: FFLAGS += $(FIXEDFLAGS)
 $(BUILD_DIR)/swanser.o: FFLAGS += $(FIXEDFLAGS)
 $(BUILD_DIR)/swmod1.o: FFLAGS += $(FIXEDFLAGS)
 $(BUILD_DIR)/swmod2.o: FFLAGS += $(FIXEDFLAGS)
 $(BUILD_DIR)/m_constants.o: FFLAGS += $(FREEFLAGS)
 $(BUILD_DIR)/m_fileio.o: FFLAGS += $(FREEFLAGS)
 $(BUILD_DIR)/mod_xnl4v5.o: FFLAGS += $(FREEFLAGS)
 $(BUILD_DIR)/serv_xnl4v5.o: FFLAGS += $(FREEFLAGS)
endif
