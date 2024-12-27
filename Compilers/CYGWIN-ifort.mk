# git $Id$
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Copyright (c) 2002-2025 The ROMS Group                                :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.md                                                 :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
#
# Include file for Intel compiler on Cygwin
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

# LD             Program to load the objects into an executable or shared library
# LDFLAGS        Flags to the loader
# RANLIB         Name of ranlib command
# MDEPFLAGS      Flags for sfmakedepend  (-s if you keep .f files)
#
# First the defaults
#

# This needs an if block to prevent building the executable
# when none is requested.
ifdef BIN
              BIN := $(BIN).exe
endif

               FC := ifort
           FFLAGS := /align /G7 /MD
       FIXEDFLAGS := /nofree
        FREEFLAGS := /free
              CPP := /usr/bin/cpp
         CPPFLAGS := -P -DCYGWIN -DCYGWIN_ifort -traditional
           INCDIR := /usr/include /usr/local/bin
            SLIBS := -L/usr/local/lib -L/usr/lib -L/usr/lib64
            ULIBS :=
             LIBS :=
         ROMS_LIB := -L$(BUILD_DIR) -lROMS
       MOD_SUFFIX := mod
               LD := $(FC)
          LDFLAGS := /link /stack:67108864
               AR := ar
          ARFLAGS := r
            MKDIR := mkdir -p
               CP := cp -p -v
               RM := rm -f
           RANLIB := ranlib
             PERL := perl
             TEST := test
      ST_LIB_NAME := libROMS.a
      SH_LIB_NAME := cygROMS.dll

#--------------------------------------------------------------------------
# Compiling flags for ROMS Applications.
#--------------------------------------------------------------------------

ifdef USE_ROMS
 ifdef USE_DEBUG
           FFLAGS += /debug:full
           FFLAGS += /traceback
           FFLAGS += /Od
           FFLAGS += /Zi
 else
           FFLAGS += /Ox
 endif
 ifdef SHARED
           FFLAGS += -fPIC
       SH_LDFLAGS += -shared
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
           FFLAGS += -fbounds-check
#          FFLAGS += -fcheck=all
 else
           FFLAGS := -O3 -w
 endif
endif

#--------------------------------------------------------------------------
# Coupled models.  Notice Linux needs the libraries repeated for
# dependencies for some of the coupled components.
#--------------------------------------------------------------------------

ifdef USE_COAMPS
             LIBS += $(COAMPS_LIB_DIR)/coamps_driver.a
             LIBS += $(COAMPS_LIB_DIR)/libaa.a
             LIBS += $(COAMPS_LIB_DIR)/libam.a
             LIBS += $(COAMPS_LIB_DIR)/libashare.a
             LIBS += $(COAMPS_LIB_DIR)/libcoamps.a
             LIBS += $(COAMPS_LIB_DIR)/libfnoc.a
             LIBS += $(COAMPS_LIB_DIR)/libaa.a
             LIBS += $(COAMPS_LIB_DIR)/libam.a
             LIBS += $(COAMPS_LIB_DIR)/libashare.a
             LIBS += $(COAMPS_LIB_DIR)/libcoamps.a
             LIBS += $(COAMPS_LIB_DIR)/libfnoc.a
             LIBS += $(COAMPS_LIB_DIR)/libfishpak.a
             LIBS += $(COAMPS_LIB_DIR)/libtracer.a
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
ifneq ($(NETCDFPAR),)
             LIBS += $(WRF_LIB_DIR)/libwrfio_nfpar.a
endif
             LIBS += $(WRF_LIB_DIR)/libwrfio_nf.a
endif

#--------------------------------------------------------------------------
# Library locations, can be overridden by environment variables.
# These are be specified in Unix form and will be converted as
# necessary to Windows form for Windows-native commands. The default
# values below assume that Cygwin mounts have been defined pointing to
# the NETCDF and ARPACK library locations.
#

ifdef USE_NETCDF4
        NC_CONFIG ?= nc-config
   TEST_NC_CONFIG := $(shell which $(NC_CONFIG))
  ifneq ($(TEST_NC_CONFIG),)
             LIBS += $(shell $(NC_CONFIG) --libs)
  endif
        NF_CONFIG ?= nf-config
    NETCDF_INCDIR ?= $(shell $(NF_CONFIG) --prefix)/include
             LIBS += $(shell $(NF_CONFIG) --flibs)
           INCDIR += $(NETCDF_INCDIR) $(INCDIR)
else
    NETCDF_INCDIR ?= /usr/include
    NETCDF_LIBDIR ?= /usr/lib
      NETCDF_LIBS ?= -lnetcdf
             LIBS += -L$(NETCDF_LIBDIR) $(NETCDF_LIBS)
           INCDIR += $(NETCDF_INCDIR) $(INCDIR)
endif

ifdef USE_HDF5
      HDF5_INCDIR ?= /usr/include
      HDF5_LIBDIR ?= /usr/lib
        HDF5_LIBS ?= -lhdf5_fortran -lhdf5hl_fortran -lhdf5 -lz
             LIBS += -L$(HDF5_LIBDIR) $(HDF5_LIBS)
           INCDIR += $(HDF5_INCDIR)
endif

ifdef USE_ARPACK
 ifdef USE_MPI
   PARPACK_LIBDIR ?= /usr/local/lib
             LIBS += -L$(PARPACK_LIBDIR) -lparpack
 endif
    ARPACK_LIBDIR ?= /usr/local/lib
             LIBS += -L$(ARPACK_LIBDIR) -larpack
endif

ifdef USE_MPI
       MCT_LIBDIR ?= c:\\work\\models\\MCT_v2.2\\mct
      MPEU_LIBDIR ?= c:\\work\\models\\MCT_v2.2\\mpeu
       MPI_INCDIR ?= c:\\work\\models\\MPICH2\\include
       MPI_LIBDIR ?= c:\\work\\models\\MPICH2\\lib
       LIBS_WIN32 += "$(MPI_LIBDIR)\fmpich2.lib "
         CPPFLAGS += -DMPI -I$(MPI_INCDIR)
           FFLAGS += -I$(MPI_INCDIR)
endif

ifdef USE_OpenMP
         CPPFLAGS += -D_OPENMP
           FFLAGS += /Qopenmp /Qopenmp_report1
endif

ifdef USE_MCT
         CPPFLAGS += -traditional-cpp
           FFLAGS += -I$(MCT_LIBDIR) -I$(MPEU_LIBDIR)
           FFLAGS += /noextend_source -assume:byterecl
       LIBS_WIN32 += "$(MCT_LIBDIR)\libmct.a" "$(MPEU_LIBDIR)\libmpeu.a"
           INCDIR += $(MCT_INCDIR) $(INCDIR)
endif

ifdef USE_ESMF
                     include $(ESMFMKFILE)
          ESMF_OS ?= $(OS)
      ESMF_SUBDIR := $(ESMF_OS).$(ESMF_COMPILER).$(ESMF_ABI).$(ESMF_COMM).$(ESMF_SITE)
      ESMF_MK_DIR ?= $(ESMF_DIR)/lib/lib$(ESMF_BOPT)/$(ESMF_SUBDIR)
           FFLAGS += $(ESMF_F90COMPILEPATHS)
       LIBS_WIN32 += $(ESMF_F90LINKPATHS) $(ESMF_F90ESMFLINKLIBS)
endif

#
# For a Windows compiler, create variables pointing to the Windows
# file names needed when linking. Use of the "=" sign means that
# variables will be evaluated only when needed.
#
         BIN_WIN32 = "$$(cygpath --windows $(BIN))"
        LIBS_WIN32 += "$$(cygpath --windows $(NETCDF_LIB))"
        LIBS_WIN32 += "c:\cygwin\lib\gcc\i686-pc-mingw32\3.4.4\libgcc.a"
ifdef USE_ARPACK
        LIBS_WIN32 += "$$(cygpath --windows $(ARPACK_LIB))"
endif

        LD_WINDOWS := on

#
# Use full path of compiler. Notice that a very special editing is
# done for the FC defintion to allow blank spaces and parenthesis in the
# compiler path. For example, we can have:
#
#      c:\\Software\\My Compilers (64 bit)
# or
#      /cygdrive/c/Program Files (x86)/Intel/Compiler/11.1/051/bin/ia32/ifort
# or
#      /cygdrive/c/Program Files/Intel/Compiler/11.1/051/bin/ia32/ifort
#
               FC := $(shell which ${FC} | sed 's|\([ |(|)]\)|\\\1|g')
               LD := $(FC)

#--------------------------------------------------------------------------
# ROMS specific rules.
#--------------------------------------------------------------------------

# CYGWIN can only load user compiled .dll files located in the same
# directory as the executable. This rule will copy the cygROMS.dll
# to the $(BINDIR) so the executable can run. This is only needed when
# EXEC and SHARED are set and STATIC is NOT. If STATIC is set then
# CYGWIN will automatically link with the static library.

# Blank it out to be sure
       CYG_DLL_CP :=

ifdef SHARED
 ifdef EXEC
  ifndef STATIC
       CYG_DLL_CP := cyg_dll_cp

.PHONY: cyg_dll_cp
cyg_dll_cp: $(BIN)
            $(CP) $(BUILD_DIR)/$(SH_LIB_NAME) $(BINDIR)

  endif
 endif
endif

#
# For a Windows compiler, override the compilation rule
#

%.o: %.f90
	cd $(BUILD_DIR); $(FC) -c $(FFLAGS) $(notdir $<) /object:$(notdir $@)

#
# Override rule for object files. Use .o instead of .obj
#

define one-compile-rule
  $1: $2 $3
	cd $$(BUILD_DIR); $$(FC) -c $$(FFLAGS) $(notdir $2) /object:$(notdir $1)

  $2: $3
	$$(CPP) $$(CPPFLAGS) $$(MY_CPP_FLAGS) $$< > $$@
	$$(CLEAN) $$@

endef

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
