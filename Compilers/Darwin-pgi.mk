# git $Id$
# svn $Id: Darwin-pgi.mk 1120 2022-04-08 19:14:36Z arango $
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Copyright (c) 2002-2022 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
# Include file for PGI compiler on Darwin
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
               FC := pgf90
           FFLAGS :=
       FIXEDFLAGS := -Mnofree
        FREEFLAGS := -Mfree
              CPP := /usr/bin/cpp
         CPPFLAGS := -P -traditional-cpp -w          # -w turns off warnings
           INCDIR := /usr/include /usr/local/bin
            SLIBS := -L/usr/local/lib -L/usr/lib
            ULIBS :=
             LIBS :=
         ROMS_LIB := -L$(SCRATCH_DIR) -lROMS
       MOD_SUFFIX := mod
               LD := $(FC)
          LDFLAGS :=
               AR := ar
          ARFLAGS := r
            MKDIR := mkdir -p
               CP := cp -p -v
               RM := rm -f
           RANLIB := ranlib -no_warning_for_no_symbols
             PERL := perl
             TEST := test
      ST_LIB_NAME := libROMS.a
      SH_LIB_NAME := libROMS.dylib

#--------------------------------------------------------------------------
# Compiling flags for ROMS Applications.
#--------------------------------------------------------------------------
#
# Perform floating-point operations in strict conformance with the
# IEEE standard by using -Kieee. This may slow down computations
# because some optimizations are disabled.  However, we noticed a
# speed-up. The user may want to have this option to allow similar,
# if not identical solutions between different of the PGI compiler.

ifdef USE_ROMS
 ifdef USE_DEBUG
           FFLAGS += -g
           FFLAGS += -C
           FFLAGS += -Mchkstk -Mchkfpstk
 else
           FFLAGS += -fastsse -Mipa=fast
 endif
           FFLAGS += -Kieee
 ifdef SHARED
          LDFLAGS += -Wl,-rpath,$(SCRATCH_DIR)

       SH_LDFLAGS += -dynamiclib
       SH_LDFLAGS += -Wl,-headerpad_max_install_names
       SH_LDFLAGS += -Wl,-undefined,dynamic_lookup

       SH_LDFLAGS += -install_name @rpath/$(SH_LIB_NAME)
 endif

        MDEPFLAGS := --cpp --fext=f90 --file=- --objdir=$(SCRATCH_DIR)
endif

#--------------------------------------------------------------------------
# Compiling flags for CICE Applications.
#--------------------------------------------------------------------------

ifdef CICE_APPLICATION
          CPPDEFS := -DLINUS $(MY_CPP_FLAGS)
 ifdef USE_DEBUG
           FFLAGS += -g
           FFLAGS += -C
           FFLAGS += -Mchkstk -Mchkfpstk
 else
           FFLAGS += -fastsse -Mipa=fast
 endif
           FFLAGS += -Kieee
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
 ifeq "$(strip $(WRF_LIB_DIR))" "$(WRF_SRC_DIR)"
             LIBS += $(WRF_LIB_DIR)/main/module_wrf_top.o
             LIBS += $(WRF_LIB_DIR)/main/libwrflib.a
             LIBS += $(WRF_LIB_DIR)/external/fftpack/fftpack5/libfftpack.a
             LIBS += $(WRF_LIB_DIR)/external/io_grib1/libio_grib1.a
             LIBS += $(WRF_LIB_DIR)/external/io_grib_share/libio_grib_share.a
             LIBS += $(WRF_LIB_DIR)/external/io_int/libwrfio_int.a
             LIBS += $(WRF_LIB_DIR)/external/esmf_time_f90/libesmf_time.a
             LIBS += $(WRF_LIB_DIR)/external/RSL_LITE/librsl_lite.a
             LIBS += $(WRF_LIB_DIR)/frame/module_internal_header_util.o
             LIBS += $(WRF_LIB_DIR)/frame/pack_utils.o
             LIBS += $(WRF_LIB_DIR)/external/io_netcdf/libwrfio_nf.a
     WRF_MOD_DIRS  = main frame phys share external/esmf_time_f90
 else
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
endif

#--------------------------------------------------------------------------
# Library locations, can be overridden by environment variables.
#--------------------------------------------------------------------------
#
# According to the PGI manual, the -u -Bstatic flags initializes
# the symbol table with -Bstatic, which is undefined for the linker.
# An undefined symbol triggers loading of the first member of an
# archive library. The -u flag fails with version 7.x of the compiler
# because it expects an argument.


ifdef USE_NETCDF4
        NF_CONFIG ?= nf-config
    NETCDF_INCDIR ?= $(shell $(NF_CONFIG) --prefix)/include
             LIBS += $(shell $(NF_CONFIG) --flibs)
           INCDIR += $(NETCDF_INCDIR) $(INCDIR)
else
    NETCDF_INCDIR ?= /opt/pgisoft/serial/netcdf3/include
    NETCDF_LIBDIR ?= /opt/pgisoft/serial/netcdf3/lib
      NETCDF_LIBS ?= -lnetcdf
             LIBS += -L$(NETCDF_LIBDIR) $(NETCDF_LIBS)
           INCDIR += $(NETCDF_INCDIR) $(INCDIR)
endif

ifdef USE_HDF5
      HDF5_INCDIR ?= /opt/pgisoft/serial/hdf5/include
      HDF5_LIBDIR ?= /opt/pgisoft/serial/hdf5/lib
        HDF5_LIBS ?= -lhdf5_fortran -lhdf5hl_fortran -lhdf5 -lz
             LIBS += -L$(HDF5_LIBDIR) $(HDF5_LIBS)
           INCDIR += $(HDF5_INCDIR)
endif

ifdef USE_ARPACK
 ifdef USE_MPI
   PARPACK_LIBDIR ?= /opt/pgisoft/PARPACK
             LIBS += -L$(PARPACK_LIBDIR) -lparpack
 endif
    ARPACK_LIBDIR ?= /opt/pgisoft/ARPACK
             LIBS += -L$(ARPACK_LIBDIR) -larpack
endif

ifdef USE_MPI
         CPPFLAGS += -DMPI
 ifdef USE_MPIF90
               FC := mpif90
 else
             LIBS += -Bdynamic -lfmpi-pgi -lmpi-pgi -Bstatic
 endif
endif

ifdef USE_OpenMP
         CPPFLAGS += -D_OPENMP
           FFLAGS += -mp
endif

ifdef USE_MCT
       MCT_INCDIR ?= /opt/pgisoft/mct/include
       MCT_LIBDIR ?= /opt/pgisoft/mct/lib
           FFLAGS += -I$(MCT_INCDIR)
             LIBS += -L$(MCT_LIBDIR) -lmct -lmpeu
           INCDIR += $(MCT_INCDIR) $(INCDIR)
endif

ifdef USE_ESMF
          ESMF_OS ?= $(OS)
      ESMF_SUBDIR := $(ESMF_OS).$(ESMF_COMPILER).$(ESMF_ABI).$(ESMF_COMM).$(ESMF_SITE)
      ESMF_MK_DIR ?= $(ESMF_DIR)/lib/lib$(ESMF_BOPT)/$(ESMF_SUBDIR)
                     include $(ESMF_MK_DIR)/esmf.mk
           FFLAGS += $(ESMF_F90COMPILEPATHS)
             LIBS += $(ESMF_F90LINKPATHS) $(ESMF_F90ESMFLINKLIBS)
endif

# Use full path of compiler.

               FC := $(shell which ${FC})
               LD := $(FC)

#--------------------------------------------------------------------------
# ROMS specific rules.
#--------------------------------------------------------------------------

# Set free form format in some ROMS source files to allow long string for
# local directory and compilation flags inside the code.

ifdef USE_ROMS
 $(SCRATCH_DIR)/mod_ncparam.o: FFLAGS += $(FREEFLAGS)
 $(SCRATCH_DIR)/mod_strings.o: FFLAGS += $(FREEFLAGS)
 $(SCRATCH_DIR)/analytical.o: FFLAGS += $(FREEFLAGS)
 $(SCRATCH_DIR)/biology.o: FFLAGS += $(FREEFLAGS)

 ifdef USE_ADJOINT
  $(SCRATCH_DIR)/ad_biology.o: FFLAGS += $(FREEFLAGS)
 endif
 ifdef USE_REPRESENTER
  $(SCRATCH_DIR)/rp_biology.o: FFLAGS += $(FREEFLAGS)
 endif
 ifdef USE_TANGENT
  $(SCRATCH_DIR)/tl_biology.o: FFLAGS += $(FREEFLAGS)
 endif
endif

#--------------------------------------------------------------------------
# Model coupling specific rules.
#--------------------------------------------------------------------------

# Add COAMPS library directory to include path of ESMF coupling files.

ifdef USE_COAMPS
 $(SCRATCH_DIR)/esmf_atm.o: FFLAGS += -I$(COAMPS_LIB_DIR)
 $(SCRATCH_DIR)/esmf_esm.o: FFLAGS += -I$(COAMPS_LIB_DIR)
endif

# Add WRF library directory to include path of ESMF coupling files.

ifdef USE_WRF
 ifeq "$(strip $(WRF_LIB_DIR))" "$(WRF_SRC_DIR)"
  $(SCRATCH_DIR)/esmf_atm.o: FFLAGS += $(addprefix -I$(WRF_LIB_DIR)/,$(WRF_MOD_DIRS))
 else
  $(SCRATCH_DIR)/esmf_atm.o: FFLAGS += -I$(WRF_LIB_DIR)
 endif
endif

# Supress free format in SWAN source files since there are comments
# beyond column 72.

ifdef USE_SWAN
 $(SCRATCH_DIR)/ocpcre.o: FFLAGS += $(FIXEDFLAGS)
 $(SCRATCH_DIR)/ocpids.o: FFLAGS += $(FIXEDFLAGS)
 $(SCRATCH_DIR)/ocpmix.o: FFLAGS += $(FIXEDFLAGS)
 $(SCRATCH_DIR)/swancom1.o: FFLAGS += $(FIXEDFLAGS)
 $(SCRATCH_DIR)/swancom2.o: FFLAGS += $(FIXEDFLAGS)
 $(SCRATCH_DIR)/swancom3.o: FFLAGS += $(FIXEDFLAGS)
 $(SCRATCH_DIR)/swancom4.o: FFLAGS += $(FIXEDFLAGS)
 $(SCRATCH_DIR)/swancom5.o: FFLAGS += $(FIXEDFLAGS)
 $(SCRATCH_DIR)/swanmain.o: FFLAGS += $(FIXEDFLAGS)
 $(SCRATCH_DIR)/swanout1.o: FFLAGS += $(FIXEDFLAGS)
 $(SCRATCH_DIR)/swanout2.o: FFLAGS += $(FIXEDFLAGS)
 $(SCRATCH_DIR)/swanparll.o: FFLAGS += $(FIXEDFLAGS)
 $(SCRATCH_DIR)/swanpre1.o: FFLAGS += $(FIXEDFLAGS)
 $(SCRATCH_DIR)/swanpre2.o: FFLAGS += $(FIXEDFLAGS)
 $(SCRATCH_DIR)/swanser.o: FFLAGS += $(FIXEDFLAGS)
 $(SCRATCH_DIR)/swmod1.o: FFLAGS += $(FIXEDFLAGS)
 $(SCRATCH_DIR)/swmod2.o: FFLAGS += $(FIXEDFLAGS)
 $(SCRATCH_DIR)/m_constants.o: FFLAGS += $(FREEFLAGS)
 $(SCRATCH_DIR)/m_fileio.o: FFLAGS += $(FREEFLAGS)
 $(SCRATCH_DIR)/mod_xnl4v5.o: FFLAGS += $(FREEFLAGS)
 $(SCRATCH_DIR)/serv_xnl4v5.o: FFLAGS += $(FREEFLAGS)
endif
