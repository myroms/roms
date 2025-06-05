# git $Id$
#::::::::::::::::::::::::::::::::::::::::::::::::::::: Hernan G. Arango :::
# Copyright (c) 2002-2025 The ROMS Group                  Kate Hedstrom :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.md                                                 :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#                                                                       :::
#  ROMS Framework Master Makefile                                       :::
#                                                                       :::
#  This makefile is designed to work only with GNU Make version 3.80 or :::
#  higher. It can be used in any architecture provided that there is a  :::
#  machine/compiler rules file in the  "Compilers"  subdirectory.  You  :::
#  may need to modify the rules file to specify the  correct path  for  :::
#  the NetCDF and ARPACK libraries. The ARPACK library is only used in  :::
#  the Generalized Stability Theory analysis and Laczos algorithm.      :::
#                                                                       :::
#  If appropriate,  the USER needs to modify the  macro definitions in  :::
#  in user-defined section below.  To activate an option set the macro  :::
#  to "on". For example, if you want to compile with debugging options  :::
#  set:                                                                 :::
#                                                                       :::
#      USE_DEBUG := on                                                  :::
#                                                                       :::
#  Otherwise, leave macro definition blank.                             :::
#                                                                       :::
#  The USER needs to provide a value for the  macro FORT.  Choose  the  :::
#  appropriate value from the list below.                               :::
#                                                                       :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

ifneq (3.80,$(firstword $(sort $(MAKE_VERSION) 3.80)))
 $(error This makefile requires GNU make version 3.80 or higher. \
		Your current version is: $(MAKE_VERSION))
endif

#--------------------------------------------------------------------------
#  Initialize some things.
#--------------------------------------------------------------------------

sources :=

#--------------------------------------------------------------------------
#  Check that at least one of SHARED, STATIC, or EXEC are set. If none are
#  set, then ROMS defaults to creating a statically linked executable.
#  This is to safeguard against old build scripts and misconfigured
#  new ones (e.g. EXEC defined but neither SHARED nor STATIC defined).
#--------------------------------------------------------------------------

ifndef SHARED
 ifndef STATIC
    STATIC := on
  ifndef EXEC
    EXEC   := on
  endif
 endif
endif

#==========================================================================
#  Start of user-defined options. In some macro definitions below: "on" or
#  any other string means TRUE while blank (or spaces) is FALSE.
#==========================================================================
#
#  The CPP option defining a particular application is specified below.
#  See header file "ROMS/Include/cppdefs.h" for all available idealized
#  and realistic applications CPP flags. For example, to activate the
#  upwelling test case (UPWELLING) set:
#
#    ROMS_APPLICATION ?= UPWELLING
#
#  Notice that this makefile will include the associated application header
#  file, which is located either in the "ROMS/Include" or MY_HEADER_DIR
#  directory.  This makefile is designed to search in both directories.
#  The only constrain is that the application CPP option must be unique
#  and header file name is the lowercase value of ROMS_APPLICATION with
#  the .h extension. For example, the upwelling application includes the
#  "upwelling.h" header file.

ROMS_APPLICATION ?= UPWELLING

#  If application header files is not located in "ROMS/Include",
#  provide an alternate directory FULL PATH.

MY_HEADER_DIR ?=

#  If your application requires analytical expressions and they are
#  not located in "ROMS/Functionals", provide an alternate directory.
#  Notice that a set analytical expressions templates can be found in
#  "User/Functionals".
#
#  If applicable, also used this directory to place your customized
#  biology model header file (like fennel.h, nemuro.h, ecosim.h, etc).

MY_ANALYTICAL_DIR ?=

#  Sometimes it is desirable to activate one or more CPP options to
#  run different variants of the same application without modifying
#  its header file. If this is the case, specify such options here
#  using the -D syntax.  For example, to write time-averaged fields
#  set:
#
#    MY_CPP_FLAGS ?= -DAVERAGES
#

MY_CPP_FLAGS ?=

#  Activate debugging compiler options:

   USE_DEBUG ?=

#  If parallel applications, use at most one of these definitions
#  (leave both definitions blank in serial applications):

     USE_MPI ?=
  USE_OpenMP ?=

#  If distributed-memory, turn on compilation via the script "mpif90".
#  This is needed in some Linux operating systems. In some systems with
#  native MPI libraries the compilation does not require MPICH type
#  scripts. This macro is also convient when there are several fortran
#  compiliers (ifort, pgf90, pathf90) in the system that use mpif90.
#  In this, case the user need to select the desired compiler below and
#  turn on both USE_MPI and USE_MPIF90 macros.

  USE_MPIF90 ?=

#  If applicable, activate 64-bit compilation:

   USE_LARGE ?= on

#  If applicable, link with NetCDF-4 library. Notice that the NetCDF-4
#  library needs both the HDF5 and MPI libraries.

 USE_NETCDF4 ?=

#--------------------------------------------------------------------------
#  We are going to include a file with all the settings that depend on
#  the system and the compiler. We are going to build up the name of the
#  include file using information on both. Set your compiler here from
#  the following list:
#
#  Operating System        Compiler(s)
#
#     AIX:                    xlf
#     ALPHA:                  f90
#     CYGWIN:                 g95, df, ifort
#     Darwin:                 f90, xlf
#     IRIX:                   f90
#     Linux:                  ftn, ifc, ifort, pgi, path, g95, gfortran
#     SunOS:                  f95
#     UNICOS-mp:              ftn
#     SunOS/Linux:            ftn (Cray cross-compiler)
#
#  Feel free to send us additional rule files to include! Also, be sure
#  to check the appropriate file to make sure it has the right paths to
#  NetCDF and so on.
#--------------------------------------------------------------------------

        FORT ?= ifort

#--------------------------------------------------------------------------
#  Set directory for executable.
#--------------------------------------------------------------------------

      BINDIR ?= .

#==========================================================================
#  End of user-defined options. See also the machine-dependent include
#  file being used above.
#==========================================================================

#--------------------------------------------------------------------------
#  Set ROMS Build directory for processing and compiling files.
#--------------------------------------------------------------------------

BUILD_DIR ?= Build_roms

#  Backward compatability with old build scripts and make configuration
#  files (*.mk). The BUILD_DIR macro is preferred.

ifdef SCRATCH_DIR
  BUILD_DIR := $(SCRATCH_DIR)
else
  SCRATCH_DIR := $(BUILD_DIR)
endif

#  Define cleaning macros for compiling.

clean_list := core *.ipo $(BUILD_DIR)

ifeq "$(strip $(BUILD_DIR))" "."
  clean_list := core *.o *.oo *.mod *.f90 lib*.a *.bak
  clean_list += $(CURDIR)/*.ipo
endif
ifeq "$(strip $(BUILD_DIR))" "./"
  clean_list := core *.o *.oo *.ipo *.mod *.f90 lib*.a *.bak
  clean_list += $(CURDIR)/*.ipo
endif

#--------------------------------------------------------------------------
#  Notice that the token "libraries" is initialized with the ROMS/Utility
#  library to account for calls to objects in other ROMS libraries or
#  cycling dependencies. These types of dependencies are problematic in
#  some compilers during linking. Such libraries appear twice at linking
#  step (beginning and almost the end of ROMS library list).
#--------------------------------------------------------------------------

libraries :=

#--------------------------------------------------------------------------
#  Set Pattern rules.
#--------------------------------------------------------------------------

%.o: %.F

%.o: %.f90
	cd $(BUILD_DIR); $(FC) -c $(FFLAGS) $(notdir $<)

%.f90: %.F
	$(CPP) $(CPPFLAGS) $(MY_CPP_FLAGS) $< > $*.f90
	$(CLEAN) $*.f90

CLEAN := ROMS/Bin/cpp_clean

#--------------------------------------------------------------------------
#  Set C-preprocessing flags associated with ROMS application. They are
#  used in "ROMS/Include/cppdefs.h" to include the appropriate application
#  header file.
#--------------------------------------------------------------------------

ifdef ROMS_APPLICATION
        HEADER := $(addsuffix .h, \
			$(shell echo ${ROMS_APPLICATION} | tr [A-Z] [a-z]))
 ROMS_CPPFLAGS := -D$(ROMS_APPLICATION)
 ROMS_CPPFLAGS += -D'HEADER="$(HEADER)"'
 ifdef MY_HEADER_DIR
  ROMS_CPPFLAGS += -D'ROMS_HEADER="$(MY_HEADER_DIR)/$(HEADER)"'
 else
  ROMS_CPPFLAGS += -D'ROMS_HEADER="$(HEADER)"'
 endif
 ifdef MY_CPP_FLAGS
  ROMS_CPPFLAGS += $(MY_CPP_FLAGS)
 endif
endif

#--------------------------------------------------------------------------
#  Internal macro definitions used to select the code to compile and
#  additional libraries to link. It uses the CPP activated in the
#  header file ROMS/Include/cppdefs.h to determine macro definitions.
#--------------------------------------------------------------------------

  COMPILERS ?= $(CURDIR)/Compilers

MAKE_MACROS := $(shell echo ${HOME} | sed 's| |\\ |g')/make_macros.mk

ifneq ($(MAKECMDGOALS),clean)
  ifneq ($(MAKECMDGOALS),tarfile)
    MACROS := $(shell cpp -P $(ROMS_CPPFLAGS) Compilers/make_macros.h > \
                $(MAKE_MACROS); $(CLEAN) $(MAKE_MACROS))

    GET_MACROS := $(wildcard $(BUILD_DIR)/make_macros.*)

    ifdef GET_MACROS
      include $(BUILD_DIR)/make_macros.mk
    else
      include $(MAKE_MACROS)
    endif
  endif
endif

clean_list += $(MAKE_MACROS)

#--------------------------------------------------------------------------
#  Make functions for putting the temporary files in $(BUILD_DIR)
#  DO NOT modify this section; spaces and blank lines are needed.
#--------------------------------------------------------------------------

# $(call source-dir-to-binary-dir, directory-list)
source-dir-to-binary-dir = $(addprefix $(BUILD_DIR)/, $(notdir $1))

# $(call source-to-object, source-file-list)
source-to-object = $(call source-dir-to-binary-dir,   \
                   $(subst .F,.o,$1))

# $(call make-static-library, library-name, source-file-list)
define make-static-library
   sources   += $2

   $(BUILD_DIR)/$1: $(call source-dir-to-binary-dir,    \
                      $(subst .F,.o,$2))
	$(AR) $(ARFLAGS) $$@ $$^
	$(RANLIB) $$@
endef

# $(call make-shared-library, library-name, source-file-list)
define make-shared-library
   $(BUILD_DIR)/$1: $(call source-dir-to-binary-dir,    \
                      $(subst .F,.o,$2))
	$(LD) $(FFLAGS) $(SH_LDFLAGS) -o $$@ $$^ $(LIBS)
endef

# $(call f90-source, source-file-list)
f90-source = $(call source-dir-to-binary-dir,     \
                   $(subst .F,.f90,$1))

# $(compile-rules)
define compile-rules
  $(foreach f, $(local_src),       \
    $(call one-compile-rule,$(call source-to-object,$f), \
    $(call f90-source,$f),$f))
endef

# $(call one-compile-rule, binary-file, f90-file, source-files)
define one-compile-rule
  $1: $2 $3
	cd $$(BUILD_DIR); $$(FC) -c $$(FFLAGS) $(notdir $2)

  $2: $3
	$$(CPP) $$(CPPFLAGS) $$(MY_CPP_FLAGS) $$< > $$@
	$$(CLEAN) $$@

endef

#--------------------------------------------------------------------------
#  Set ROMS executable file name.
#--------------------------------------------------------------------------

ifdef EXEC
  ifdef USE_DEBUG
    BIN ?= $(BINDIR)/romsG
  else
   ifdef USE_MPI
     BIN ?= $(BINDIR)/romsM
   else
    ifdef USE_OpenMP
      BIN ?= $(BINDIR)/romsO
    else
      BIN ?= $(BINDIR)/romsS
    endif
   endif
  endif
endif

#--------------------------------------------------------------------------
#  Set name of module files for netCDF F90 interface. On some platforms
#  these will need to be overridden in the machine-dependent include file.
#--------------------------------------------------------------------------

   NETCDF_MODFILE := netcdf.mod
TYPESIZES_MODFILE := typesizes.mod

#--------------------------------------------------------------------------
#  "uname -s" should return the OS or kernel name and "uname -m" should
#  return the CPU or hardware name. In practice the results can be pretty
#  flaky. Run the results through sed to convert "/" and " " to "-",
#  then apply platform-specific conversions.
#--------------------------------------------------------------------------

OS := $(shell uname -s | sed 's/[\/ ]/-/g')
OS := $(patsubst CYGWIN_%,CYGWIN,$(OS))
OS := $(patsubst MINGW%,MINGW,$(OS))
OS := $(patsubst sn%,UNICOS-sn,$(OS))

CPU := $(shell uname -m | sed 's/[\/ ]/-/g')

GITURL := $(shell git config remote.origin.url)
GITREV := $(shell git log -n 1 --format=%H)

ROOTDIR := $(shell pwd)

ifndef FORT
  $(error Variable FORT not set)
endif

ifneq ($(MAKECMDGOALS),clean)
  ifneq ($(MAKECMDGOALS),tarfile)
    MKFILE := $(COMPILERS)/$(OS)-$(strip $(FORT)).mk
    include $(MKFILE)
  endif
endif

ifdef USE_MPI
 ifdef USE_OpenMP
  $(error You cannot activate USE_MPI and USE_OpenMP at the same time!)
 endif
endif

ifdef STATIC
  libraries += $(BUILD_DIR)/$(ST_LIB_NAME)
endif

ifdef SHARED
  libraries += $(BUILD_DIR)/$(SH_LIB_NAME)
endif

#--------------------------------------------------------------------------
#  Pass the platform variables to the preprocessor as macros. Convert to
#  valid, upper-case identifiers. Attach ROMS application  CPP options.
#--------------------------------------------------------------------------

CPPFLAGS += -D$(shell echo ${OS} | tr "-" "_" | tr [a-z] [A-Z])
CPPFLAGS += -D$(shell echo ${CPU} | tr "-" "_" | tr [a-z] [A-Z])
CPPFLAGS += -D$(shell echo ${FORT} | tr "-" "_" | tr [a-z] [A-Z])

CPPFLAGS += -D'ROOT_DIR="$(ROOTDIR)"'
ifdef ROMS_APPLICATION
  CPPFLAGS  += $(ROMS_CPPFLAGS)
  MDEPFLAGS += -DROMS_HEADER="$(HEADER)"
endif

ifndef MY_ANALYTICAL_DIR
  MY_ANALYTICAL_DIR := $(ROOTDIR)/ROMS/Functionals
endif
ifeq (,$(findstring ROMS/Functionals,$(MY_ANALYTICAL_DIR)))
  MY_ANALYTICAL := on
endif
CPPFLAGS += -D'ANALYTICAL_DIR="$(MY_ANALYTICAL_DIR)"'

ifdef MY_ANALYTICAL
  CPPFLAGS += -D'MY_ANALYTICAL="$(MY_ANALYTICAL)"'
endif

CPPFLAGS += -D'GIT_URL="$(GITURL)"'
CPPFLAGS += -D'GIT_REV="$(GITREV)"'

#--------------------------------------------------------------------------
#  Build target directories.
#--------------------------------------------------------------------------

.PHONY: all

all: $(BUILD_DIR) $(BUILD_DIR)/MakeDepend $(libraries) $(BIN) rm_macros $(CYG_DLL_CP)

 modules  :=
ifdef USE_ADJOINT
 modules  +=	ROMS/Adjoint \
		ROMS/Adjoint/Biology
endif
ifdef USE_REPRESENTER
 modules  +=	ROMS/Representer \
		ROMS/Representer/Biology
endif
ifdef USE_SEAICE
 modules  +=	ROMS/Nonlinear/SeaIce
endif
ifdef USE_TANGENT
 modules  +=	ROMS/Tangent \
		ROMS/Tangent/Biology
endif
 modules  +=	ROMS/Nonlinear \
		ROMS/Nonlinear/BBL \
		ROMS/Nonlinear/Biology \
		ROMS/Nonlinear/Sediment \
		ROMS/Nonlinear/Vegetation \
		ROMS/Nonlinear/WEC \
		ROMS/Functionals \
		ROMS/Utility \
		ROMS/Drivers \
		ROMS/Modules

 includes :=	ROMS/Include
ifdef MY_ANALYTICAL
 includes +=	$(MY_ANALYTICAL_DIR)
endif
ifdef USE_ADJOINT
 includes +=	ROMS/Adjoint \
		ROMS/Adjoint/Biology
endif
ifdef USE_REPRESENTER
 includes +=	ROMS/Representer \
		ROMS/Representer/Biology
endif
ifdef USE_SEAICE
 includes +=	ROMS/Nonlinear/SeaIce
endif
ifdef USE_TANGENT
 includes +=	ROMS/Tangent \
		ROMS/Tangent/Biology
endif
 includes +=	ROMS/Nonlinear \
		ROMS/Nonlinear/BBL \
		ROMS/Nonlinear/Biology \
		ROMS/Nonlinear/Sediment \
		ROMS/Nonlinear/Vegetation \
		ROMS/Utility \
		ROMS/Drivers \
                ROMS/Functionals
ifdef MY_HEADER_DIR
 includes +=	$(MY_HEADER_DIR)
endif

ifdef USE_PIO
 includes +=	$(PIO_INCDIR)
endif

ifdef USE_COAMPS
 includes +=	$(COAMPS_LIB_DIR)
endif

ifdef USE_WRF
 ifeq "$(strip $(WRF_LIB_DIR))" "$(WRF_SRC_DIR)"
  includes +=	$(addprefix $(WRF_LIB_DIR)/,$(WRF_MOD_DIRS))
 else
  includes +=	$(WRF_LIB_DIR)
 endif
endif

modules  +=	Master
includes +=	Master Compilers

vpath %.F $(modules)
vpath %.h $(includes)
vpath %.f90 $(BUILD_DIR)
vpath %.o $(BUILD_DIR)

include $(addsuffix /Module.mk,$(modules))

MDEPFLAGS += $(patsubst %,-I %,$(includes)) --silent --moddir $(BUILD_DIR)

CPPFLAGS  += $(patsubst %,-I%,$(includes))

ifdef MY_HEADER_DIR
  CPPFLAGS += -D'HEADER_DIR="$(MY_HEADER_DIR)"'
else
  CPPFLAGS += -D'HEADER_DIR="$(ROOTDIR)/ROMS/Include"'
endif

$(BUILD_DIR):
	$(shell $(TEST) -d $(BUILD_DIR) || $(MKDIR) $(BUILD_DIR) )

#--------------------------------------------------------------------------
#  Special CPP macros for mod_strings.F
#--------------------------------------------------------------------------

$(BUILD_DIR)/mod_strings.f90: CPPFLAGS += -DMY_OS='"$(OS)"' \
              -DMY_CPU='"$(CPU)"' -DMY_FORT='"$(FORT)"' \
              -DMY_FC='"$(FC)"' -DMY_FFLAGS='"$(FFLAGS)"'

#--------------------------------------------------------------------------
#  ROMS libraries.
#--------------------------------------------------------------------------

ifdef SHARED
  $(eval $(call make-shared-library,$(SH_LIB_NAME),$(sources)))
endif

ifdef STATIC
  $(eval $(call make-static-library,$(ST_LIB_NAME),$(sources)))
endif

MYLIB := libroms.a

.PHONY: libraries

libraries: $(libraries)

#--------------------------------------------------------------------------
#  Target to create ROMS dependecies.
#--------------------------------------------------------------------------

ifneq ($(MAKECMDGOALS),tarfile)
$(BUILD_DIR)/$(NETCDF_MODFILE): | $(BUILD_DIR)
	cp -f $(NETCDF_INCDIR)/$(NETCDF_MODFILE) $(BUILD_DIR)

$(BUILD_DIR)/$(TYPESIZES_MODFILE): | $(BUILD_DIR)
	cp -f $(NETCDF_INCDIR)/$(TYPESIZES_MODFILE) $(BUILD_DIR)

$(BUILD_DIR)/MakeDepend: makefile \
                           $(BUILD_DIR)/$(NETCDF_MODFILE) \
                           $(BUILD_DIR)/$(TYPESIZES_MODFILE) \
                           | $(BUILD_DIR)
	@ $(SFMAKEDEPEND) $(MDEPFLAGS) $(sources) > $(BUILD_DIR)/MakeDepend
	cp -p $(MAKE_MACROS) $(BUILD_DIR)

.PHONY: depend

SFMAKEDEPEND := ./ROMS/Bin/sfmakedepend

depend: $(BUILD_DIR)
	$(SFMAKEDEPEND) $(MDEPFLAGS) $(sources) > $(BUILD_DIR)/MakeDepend
endif

ifneq ($(MAKECMDGOALS),clean)
  -include $(BUILD_DIR)/MakeDepend
endif

#--------------------------------------------------------------------------
#  Target to create ROMS tar file.
#--------------------------------------------------------------------------

.PHONY: tarfile

tarfile:
		tar --exclude=".git" -cvf roms-4_2.tar *

.PHONY: zipfile

zipfile:
		zip -r roms-4_2.zip *

.PHONY: gzipfile

gzipfile:
		gzip -v roms-4_2.gzip *

#--------------------------------------------------------------------------
#  Cleaning targets.
#--------------------------------------------------------------------------

.PHONY: clean

clean:
	$(RM) -r $(clean_list)

.PHONY: rm_macros

rm_macros:
	$(RM) -r $(MAKE_MACROS)

#--------------------------------------------------------------------------
#  A handy debugging target. This will allow to print the value of any
#  makefile defined macro (see http://tinyurl.com/8ax3j). For example,
#  to find the value of CPPFLAGS execute:
#
#        gmake print-CPPFLAGS
#  or
#        make print-CPPFLAGS
#--------------------------------------------------------------------------

.PHONY: print-%

print-%:
	@echo $* = $($*)
# DO NOT DELETE THIS LINE - used by make depend
coupler.o: mct_coupler.h cppdefs.h ducknc.h globaldefs.h mct_roms_swan.h tile.h
coupler.o: mct_roms_wrf.h esmf_coupler.h
coupler.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/distribute.o
coupler.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_coupler.o
coupler.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_esmf_esm.o
coupler.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
coupler.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
coupler.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
coupler.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_kinds.o
coupler.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
coupler.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
coupler.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
coupler.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
coupler.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sedbed.o
coupler.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sediment.o
coupler.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
coupler.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

esmf_atm.o: cppdefs.h ducknc.h globaldefs.h esmf_atm_wrf.h esmf_atm_coamps.h
esmf_atm.o: esmf_atm_regcm.h esmf_atm_void.h
esmf_atm.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_esmf_esm.o
esmf_atm.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

esmf_data.o: cppdefs.h ducknc.h globaldefs.h
esmf_data.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/dateclock.o
esmf_data.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_esmf_esm.o
esmf_data.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
esmf_data.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
esmf_data.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
esmf_data.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
esmf_data.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
esmf_data.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_strings.o
esmf_data.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

esmf_esm.o: cppdefs.h ducknc.h globaldefs.h
esmf_esm.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/coupler.o
esmf_esm.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/esmf_atm.o
esmf_esm.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/esmf_atm.o
esmf_esm.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/esmf_atm.o
esmf_esm.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/esmf_atm.o
esmf_esm.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/esmf_data.o
esmf_esm.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/esmf_ice.o
esmf_esm.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/esmf_roms.o
esmf_esm.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/esmf_roms.o
esmf_esm.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/esmf_wav.o
esmf_esm.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/esmf_wav.o
esmf_esm.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_esmf_esm.o

esmf_ice.o: esmf_ice_void.h cppdefs.h ducknc.h globaldefs.h esmf_ice_cice.h
esmf_ice.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_esmf_esm.o

esmf_roms.o: cmeps_roms.h cppdefs.h ducknc.h globaldefs.h esmf_roms.h
esmf_roms.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/bc_2d.o
esmf_roms.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/dateclock.o
esmf_roms.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_2d.o
esmf_roms.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/get_metadata.o
esmf_roms.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_esmf_esm.o
esmf_roms.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
esmf_roms.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
esmf_roms.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
esmf_roms.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_kinds.o
esmf_roms.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_mixing.o
esmf_roms.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
esmf_roms.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
esmf_roms.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
esmf_roms.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
esmf_roms.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
esmf_roms.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
esmf_roms.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
esmf_roms.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_strings.o
esmf_roms.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o
esmf_roms.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/roms_kernel.o
esmf_roms.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/stdinp_mod.o
esmf_roms.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/stdout_mod.o
esmf_roms.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o
esmf_roms.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/yaml_parser.o

esmf_wav.o: esmf_wav_wam.h esmf_wav_void.h cppdefs.h ducknc.h globaldefs.h
esmf_wav.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_esmf_esm.o

master.o: mct_driver.h esmf_driver.h cppdefs.h ducknc.h globaldefs.h roms.h
master.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/coupler.o
master.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/esmf_esm.o
master.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_arrays.o
master.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_coupler.o
master.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_esmf_esm.o
master.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
master.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
master.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
master.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
master.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/roms_kernel.o
master.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/set_pio.o

mod_esmf_esm.o: cppdefs.h ducknc.h globaldefs.h
mod_esmf_esm.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/dateclock.o
mod_esmf_esm.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_dim.o
mod_esmf_esm.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_var.o
mod_esmf_esm.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/distribute.o
mod_esmf_esm.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/get_metadata.o
mod_esmf_esm.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/inp_decode.o
mod_esmf_esm.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
mod_esmf_esm.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_kinds.o
mod_esmf_esm.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
mod_esmf_esm.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
mod_esmf_esm.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
mod_esmf_esm.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
mod_esmf_esm.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
mod_esmf_esm.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_strings.o
mod_esmf_esm.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/stdout_mod.o
mod_esmf_esm.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

propagator.o: propagator_fte.h propagator_so_semi.h propagator_afte.h
propagator.o: propagator_hso.h propagator_so.h propagator_hop.h cppdefs.h
propagator.o: ducknc.h globaldefs.h propagator_fsv.h propagator_op.h
propagator.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/close_io.o
propagator.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/dotproduct.o
propagator.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/ini_adjust.o
propagator.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/inner2state.o
propagator.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_coupling.o
propagator.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
propagator.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
propagator.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_kinds.o
propagator.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
propagator.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
propagator.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
propagator.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
propagator.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
propagator.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
propagator.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_storage.o
propagator.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/packing.o
propagator.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/set_depth.o
propagator.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

roms_kernel.o: afte_roms.h hessian_op_roms.h rbl4dvar_roms.h tlcheck_roms.h
roms_kernel.o: obs_sen_rbl4dvar_analysis.h i4dvar_roms.h hessian_so_roms.h
roms_kernel.o: pert_roms.h tl_rbl4dvar_roms.h adsen_roms.h tl_r4dvar_roms.h
roms_kernel.o: obs_sen_r4dvar_analysis.h so_semi_roms.h tl_roms.h nl_roms.h
roms_kernel.o: rp_roms.h fsv_roms.h ad_roms.h split_r4dvar_roms.h so_roms.h
roms_kernel.o: obs_sen_rbl4dvar_forecast.h optobs_roms.h r4dvar_roms.h
roms_kernel.o: picard_roms.h fte_roms.h correlation.h obs_sen_i4dvar_analysis.h
roms_kernel.o: op_roms.h jedi_roms.h cppdefs.h ducknc.h globaldefs.h
roms_kernel.o: split_rbl4dvar_roms.h split_i4dvar_roms.h array_modes.h
roms_kernel.o: symmetry.h
roms_kernel.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/analytical.o
roms_kernel.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/array_modes.o
roms_kernel.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/close_io.o
roms_kernel.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/congrad.o
roms_kernel.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/convolve.o
roms_kernel.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/coupler.o
roms_kernel.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/dateclock.o
roms_kernel.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_dai.o
roms_kernel.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_gst.o
roms_kernel.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_impulse.o
roms_kernel.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_mod.o
roms_kernel.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_norm.o
roms_kernel.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/distribute.o
roms_kernel.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/dotproduct.o
roms_kernel.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/get_gst.o
roms_kernel.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/get_state.o
roms_kernel.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/get_wetdry.o
roms_kernel.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/i4dvar.o
roms_kernel.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/ini_adjust.o
roms_kernel.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/ini_hmixcoef.o
roms_kernel.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/inp_par.o
roms_kernel.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_arrays.o
roms_kernel.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_boundary.o
roms_kernel.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
roms_kernel.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_fourdvar.o
roms_kernel.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
roms_kernel.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
roms_kernel.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_nesting.o
roms_kernel.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
roms_kernel.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
roms_kernel.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
roms_kernel.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
roms_kernel.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
roms_kernel.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
roms_kernel.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
roms_kernel.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_storage.o
roms_kernel.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nesting.o
roms_kernel.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/normalization.o
roms_kernel.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/omega.o
roms_kernel.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/packing.o
roms_kernel.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/post_initial.o
roms_kernel.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/propagator.o
roms_kernel.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/r4dvar.o
roms_kernel.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/rbl4dvar.o
roms_kernel.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/rho_eos.o
roms_kernel.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/rpcg_lanczos.o
roms_kernel.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/set_depth.o
roms_kernel.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/set_masks.o
roms_kernel.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/set_massflux.o
roms_kernel.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/stats_modobs.o
roms_kernel.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/stdinp_mod.o
roms_kernel.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/stdout_mod.o
roms_kernel.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/stiffness.o
roms_kernel.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o
roms_kernel.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wetdry.o
roms_kernel.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wrt_dai.o
roms_kernel.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wrt_gst.o
roms_kernel.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wrt_impulse.o
roms_kernel.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wrt_ini.o
roms_kernel.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wrt_rst.o
roms_kernel.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/zeta_balance.o

i4dvar.o: cppdefs.h ducknc.h globaldefs.h
i4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/back_cost.o
i4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/background_std.o
i4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/cgradient.o
i4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/close_io.o
i4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/cost_grad.o
i4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_hessian.o
i4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_ini.o
i4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_mod.o
i4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_norm.o
i4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_std.o
i4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/get_state.o
i4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/ini_adjust.o
i4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_boundary.o
i4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
i4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_fourdvar.o
i4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
i4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_kinds.o
i4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_mixing.o
i4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
i4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
i4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
i4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
i4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
i4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
i4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
i4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
i4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/normalization.o
i4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/set_masks.o
i4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o
i4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/sum_grad.o
i4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wrt_evolved.o
i4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wrt_ini.o
i4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wrt_std.o
i4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/zeta_balance.o

r4dvar.o: cppdefs.h ducknc.h globaldefs.h
r4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/background_std.o
r4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/close_io.o
r4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/congrad.o
r4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/convolve.o
r4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_dai.o
r4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_error.o
r4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_hessian.o
r4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_impulse.o
r4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_ini.o
r4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_mod.o
r4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_norm.o
r4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_std.o
r4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/get_state.o
r4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_boundary.o
r4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
r4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_fourdvar.o
r4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
r4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
r4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
r4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
r4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
r4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
r4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
r4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
r4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
r4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/normalization.o
r4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/posterior.o
r4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/posterior_var.o
r4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/random_ic.o
r4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o
r4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wrt_hessian.o
r4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wrt_impulse.o
r4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wrt_ini.o
r4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wrt_std.o
r4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/zeta_balance.o

rbl4dvar.o: cppdefs.h ducknc.h globaldefs.h
rbl4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/background_std.o
rbl4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/comp_Jb0.o
rbl4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/congrad.o
rbl4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/convolve.o
rbl4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_error.o
rbl4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_hessian.o
rbl4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_impulse.o
rbl4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_ini.o
rbl4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_mod.o
rbl4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_norm.o
rbl4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_std.o
rbl4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/frc_iau.o
rbl4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/get_state.o
rbl4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/ini_adjust.o
rbl4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/inp_decode.o
rbl4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_boundary.o
rbl4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
rbl4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_fourdvar.o
rbl4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
rbl4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_kinds.o
rbl4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_mixing.o
rbl4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
rbl4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
rbl4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
rbl4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
rbl4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
rbl4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
rbl4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
rbl4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
rbl4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/normalization.o
rbl4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/posterior.o
rbl4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/posterior_var.o
rbl4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/random_ic.o
rbl4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/rpcg_lanczos.o
rbl4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/state_regrid.o
rbl4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o
rbl4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/sum_grad.o
rbl4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/sum_imp.o
rbl4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wrt_aug_imp.o
rbl4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wrt_error.o
rbl4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wrt_hessian.o
rbl4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wrt_impulse.o
rbl4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wrt_ini.o
rbl4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wrt_std.o
rbl4dvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/zeta_balance.o

analytical.o: ana_mask.h set_bounds.h tile.h ana_diag.h ana_btflux.h
analytical.o: ana_humid.h ana_passive.h ana_smflux.h ana_perturb.h ana_rain.h
analytical.o: ana_pair.h ana_spinning.h ana_fsobc.h ana_m2clima.h ana_tclima.h
analytical.o: ana_m3obc.h ana_nudgcoef.h ana_srflux.h ana_sediment.h ana_sss.h
analytical.o: ana_psource.h ana_cloud.h ana_drag.h ana_biology.h ana_winds.h
analytical.o: ana_m3clima.h ana_tobc.h ana_dqdsst.h cppdefs.h ducknc.h
analytical.o: globaldefs.h ana_grid.h ana_m2obc.h ana_wtype.h ana_specir.h
analytical.o: ana_respiration.h ana_vmix.h ana_sst.h ana_initial.h ana_wwave.h
analytical.o: ana_sponge.h ana_ssh.h ana_tair.h ana_stflux.h ana_scope.h
analytical.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/dateclock.o
analytical.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/distribute.o
analytical.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/erf.o
analytical.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_2d.o
analytical.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_3d.o
analytical.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_biology.o
analytical.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_boundary.o
analytical.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_clima.o
analytical.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_eclight.o
analytical.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
analytical.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
analytical.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
analytical.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_mixing.o
analytical.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
analytical.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
analytical.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
analytical.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
analytical.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
analytical.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sedbed.o
analytical.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sediment.o
analytical.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sources.o
analytical.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
analytical.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o
analytical.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/stats.o

mod_arrays.o: cppdefs.h ducknc.h globaldefs.h
mod_arrays.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_average.o
mod_arrays.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_bbl.o
mod_arrays.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_boundary.o
mod_arrays.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_clima.o
mod_arrays.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_coupling.o
mod_arrays.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_diags.o
mod_arrays.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_extract.o
mod_arrays.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
mod_arrays.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_fourdvar.o
mod_arrays.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
mod_arrays.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ice.o
mod_arrays.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
mod_arrays.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_mixing.o
mod_arrays.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_nesting.o
mod_arrays.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
mod_arrays.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
mod_arrays.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
mod_arrays.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
mod_arrays.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sedbed.o
mod_arrays.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sedflocs.o
mod_arrays.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sources.o
mod_arrays.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_tides.o
mod_arrays.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_vegetation.o
mod_arrays.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/set_pio.o

mod_average.o: set_bounds.h cppdefs.h ducknc.h globaldefs.h
mod_average.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/destroy.o
mod_average.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_kinds.o
mod_average.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
mod_average.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
mod_average.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o

mod_bbl.o: cppdefs.h ducknc.h globaldefs.h set_bounds.h
mod_bbl.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/destroy.o
mod_bbl.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_kinds.o
mod_bbl.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o

mod_behavior.o: oyster_floats_mod.h cppdefs.h ducknc.h globaldefs.h
mod_behavior.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o

mod_biology.o: red_tide_mod.h npzd_Franks_mod.h cppdefs.h ducknc.h globaldefs.h
mod_biology.o: ecosim_mod.h fennel_mod.h npzd_Powell_mod.h hypoxia_srm_mod.h
mod_biology.o: nemuro_mod.h npzd_iron_mod.h
mod_biology.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
mod_biology.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o

mod_boundary.o: set_bounds.h cppdefs.h ducknc.h globaldefs.h
mod_boundary.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/destroy.o
mod_boundary.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_kinds.o
mod_boundary.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
mod_boundary.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
mod_boundary.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o

mod_clima.o: set_bounds.h cppdefs.h ducknc.h globaldefs.h
mod_clima.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/destroy.o
mod_clima.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_kinds.o
mod_clima.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
mod_clima.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o

mod_coupler.o: cppdefs.h ducknc.h globaldefs.h
mod_coupler.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
mod_coupler.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
mod_coupler.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
mod_coupler.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
mod_coupler.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o

mod_coupling.o: cppdefs.h ducknc.h globaldefs.h set_bounds.h
mod_coupling.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/destroy.o
mod_coupling.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_kinds.o
mod_coupling.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o

mod_diags.o: cppdefs.h ducknc.h globaldefs.h set_bounds.h
mod_diags.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/destroy.o
mod_diags.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_biology.o
mod_diags.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_kinds.o
mod_diags.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o

mod_eclight.o: cppdefs.h ducknc.h globaldefs.h
mod_eclight.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_biology.o
mod_eclight.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o

mod_eoscoef.o: cppdefs.h ducknc.h globaldefs.h
mod_eoscoef.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_kinds.o

mod_extract.o: set_bounds.h cppdefs.h ducknc.h globaldefs.h
mod_extract.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/destroy.o
mod_extract.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_kinds.o
mod_extract.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
mod_extract.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o

mod_floats.o: cppdefs.h ducknc.h globaldefs.h
mod_floats.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/destroy.o
mod_floats.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
mod_floats.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o

mod_forces.o: cppdefs.h ducknc.h globaldefs.h set_bounds.h
mod_forces.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/destroy.o
mod_forces.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_biology.o
mod_forces.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_kinds.o
mod_forces.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
mod_forces.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o

mod_fourdvar.o: cppdefs.h ducknc.h globaldefs.h
mod_fourdvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
mod_fourdvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
mod_fourdvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
mod_fourdvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
mod_fourdvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
mod_fourdvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
mod_fourdvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
mod_fourdvar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

mod_grid.o: cppdefs.h ducknc.h globaldefs.h set_bounds.h
mod_grid.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/destroy.o
mod_grid.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_kinds.o
mod_grid.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
mod_grid.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o

mod_ice.o: cppdefs.h ducknc.h globaldefs.h

mod_iounits.o: cppdefs.h ducknc.h globaldefs.h
mod_iounits.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_kinds.o

mod_kinds.o: cppdefs.h ducknc.h globaldefs.h

mod_mixing.o: cppdefs.h ducknc.h globaldefs.h set_bounds.h
mod_mixing.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/destroy.o
mod_mixing.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_kinds.o
mod_mixing.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
mod_mixing.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o

mod_ncparam.o: red_tide_var.h fennel_var.h npzd_Franks_var.h hypoxia_srm_var.h
mod_ncparam.o: npzd_Powell_var.h cppdefs.h ducknc.h globaldefs.h ecosim_var.h
mod_ncparam.o: vegetation_var.h nemuro_var.h npzd_iron_var.h sediment_var.h
mod_ncparam.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/get_metadata.o
mod_ncparam.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_biology.o
mod_ncparam.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ice.o
mod_ncparam.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
mod_ncparam.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
mod_ncparam.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
mod_ncparam.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
mod_ncparam.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sediment.o
mod_ncparam.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_vegetation.o
mod_ncparam.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

mod_nesting.o: cppdefs.h ducknc.h globaldefs.h
mod_nesting.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/destroy.o
mod_nesting.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_boundary.o
mod_nesting.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_kinds.o
mod_nesting.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
mod_nesting.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o

mod_netcdf.o: cppdefs.h ducknc.h globaldefs.h
mod_netcdf.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/dateclock.o
mod_netcdf.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/distribute.o
mod_netcdf.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
mod_netcdf.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_kinds.o
mod_netcdf.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
mod_netcdf.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
mod_netcdf.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
mod_netcdf.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
mod_netcdf.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

mod_ocean.o: cppdefs.h ducknc.h globaldefs.h set_bounds.h
mod_ocean.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/destroy.o
mod_ocean.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_kinds.o
mod_ocean.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o

mod_parallel.o: cppdefs.h ducknc.h globaldefs.h
mod_parallel.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
mod_parallel.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_kinds.o
mod_parallel.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
mod_parallel.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_strings.o

mod_param.o: cppdefs.h ducknc.h globaldefs.h
mod_param.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_kinds.o

mod_pio_netcdf.o: cppdefs.h ducknc.h globaldefs.h
mod_pio_netcdf.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/dateclock.o
mod_pio_netcdf.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
mod_pio_netcdf.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_kinds.o
mod_pio_netcdf.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
mod_pio_netcdf.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
mod_pio_netcdf.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
mod_pio_netcdf.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
mod_pio_netcdf.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
mod_pio_netcdf.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_vegetation.o
mod_pio_netcdf.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

mod_scalars.o: cppdefs.h ducknc.h globaldefs.h
mod_scalars.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_kinds.o
mod_scalars.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o

mod_sedbed.o: sedbed_mod.h set_bounds.h cppdefs.h ducknc.h globaldefs.h
mod_sedbed.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/destroy.o
mod_sedbed.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_kinds.o
mod_sedbed.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
mod_sedbed.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
mod_sedbed.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sediment.o

mod_sedflocs.o: cppdefs.h ducknc.h globaldefs.h sedflocs_mod.h set_bounds.h
mod_sedflocs.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_kinds.o
mod_sedflocs.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
mod_sedflocs.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
mod_sedflocs.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sediment.o

mod_sediment.o: sediment_mod.h cppdefs.h ducknc.h globaldefs.h
mod_sediment.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o

mod_sources.o: cppdefs.h ducknc.h globaldefs.h
mod_sources.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/destroy.o
mod_sources.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
mod_sources.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_kinds.o
mod_sources.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
mod_sources.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
mod_sources.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
mod_sources.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
mod_sources.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
mod_sources.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
mod_sources.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

mod_stepping.o: cppdefs.h ducknc.h globaldefs.h

mod_storage.o: cppdefs.h ducknc.h globaldefs.h
mod_storage.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
mod_storage.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o

mod_strings.o: cppdefs.h ducknc.h globaldefs.h

mod_tides.o: cppdefs.h ducknc.h globaldefs.h set_bounds.h
mod_tides.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/destroy.o
mod_tides.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
mod_tides.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_kinds.o
mod_tides.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
mod_tides.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
mod_tides.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
mod_tides.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
mod_tides.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
mod_tides.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
mod_tides.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
mod_tides.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

mod_vegetation.o: vegetation_mod.h set_bounds.h cppdefs.h ducknc.h globaldefs.h
mod_vegetation.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_kinds.o
mod_vegetation.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o

bbl.o: mb_bbl.h tile.h set_bounds.h cppdefs.h ducknc.h globaldefs.h ssw_bbl.h
bbl.o: sg_bbl.h
bbl.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/bc_2d.o
bbl.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_bbl.o
bbl.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
bbl.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
bbl.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
bbl.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
bbl.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
bbl.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
bbl.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
bbl.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sedbed.o
bbl.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sediment.o
bbl.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
bbl.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o

bbl_output.o: cppdefs.h ducknc.h globaldefs.h
bbl_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_var.o
bbl_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/extract_sta.o
bbl_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_average.o
bbl_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_bbl.o
bbl_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
bbl_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
bbl_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
bbl_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_mixing.o
bbl_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
bbl_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
bbl_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
bbl_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
bbl_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
bbl_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
bbl_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
bbl_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sedbed.o
bbl_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sediment.o
bbl_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
bbl_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fwrite2d.o
bbl_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fwrite3d.o
bbl_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/omega.o
bbl_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

biology.o: npzd_Franks.h set_bounds.h tile.h nemuro.h npzd_Powell.h npzd_iron.h
biology.o: hypoxia_srm.h cppdefs.h ducknc.h globaldefs.h ecosim.h fennel.h
biology.o: red_tide.h
biology.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/dateclock.o
biology.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_biology.o
biology.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_diags.o
biology.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_eclight.o
biology.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
biology.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
biology.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
biology.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_kinds.o
biology.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
biology.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
biology.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
biology.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
biology.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o

biology_floats.o: cppdefs.h ducknc.h globaldefs.h oyster_floats.h
biology_floats.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_behavior.o
biology_floats.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_biology.o
biology_floats.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_floats.o
biology_floats.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
biology_floats.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
biology_floats.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
biology_floats.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
biology_floats.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
biology_floats.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o

mod_vandera_funcs.o: cppdefs.h ducknc.h globaldefs.h
mod_vandera_funcs.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_kinds.o
mod_vandera_funcs.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o

sed_bed.o: tile.h set_bounds.h cppdefs.h ducknc.h globaldefs.h
sed_bed.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/bc_3d.o
sed_bed.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_2d.o
sed_bed.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_bbl.o
sed_bed.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
sed_bed.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
sed_bed.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
sed_bed.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
sed_bed.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
sed_bed.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sedbed.o
sed_bed.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sediment.o
sed_bed.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
sed_bed.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o

sed_bed2.o: tile.h cppdefs.h ducknc.h globaldefs.h set_bounds.h
sed_bed2.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/bc_3d.o
sed_bed2.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_2d.o
sed_bed2.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_bbl.o
sed_bed2.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_biology.o
sed_bed2.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_diags.o
sed_bed2.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
sed_bed2.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
sed_bed2.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
sed_bed2.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
sed_bed2.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
sed_bed2.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sedbed.o
sed_bed2.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sediment.o
sed_bed2.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
sed_bed2.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o

sed_bed_cohesive.o: tile.h set_bounds.h cppdefs.h ducknc.h globaldefs.h
sed_bed_cohesive.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/bc_3d.o
sed_bed_cohesive.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_2d.o
sed_bed_cohesive.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_bbl.o
sed_bed_cohesive.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
sed_bed_cohesive.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
sed_bed_cohesive.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
sed_bed_cohesive.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
sed_bed_cohesive.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
sed_bed_cohesive.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sedbed.o
sed_bed_cohesive.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sediment.o
sed_bed_cohesive.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
sed_bed_cohesive.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o

sed_bedload.o: set_bounds.h cppdefs.h ducknc.h globaldefs.h tile.h
sed_bedload.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/bc_3d.o
sed_bedload.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_2d.o
sed_bedload.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_bbl.o
sed_bedload.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
sed_bedload.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
sed_bedload.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
sed_bedload.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
sed_bedload.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
sed_bedload.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
sed_bedload.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sedbed.o
sed_bedload.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sediment.o
sed_bedload.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
sed_bedload.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o

sed_bedload_vandera.o: set_bounds.h cppdefs.h ducknc.h globaldefs.h tile.h
sed_bedload_vandera.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/bc_2d.o
sed_bedload_vandera.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/bc_3d.o
sed_bedload_vandera.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_2d.o
sed_bedload_vandera.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_bbl.o
sed_bedload_vandera.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
sed_bedload_vandera.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
sed_bedload_vandera.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_kinds.o
sed_bedload_vandera.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
sed_bedload_vandera.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
sed_bedload_vandera.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
sed_bedload_vandera.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
sed_bedload_vandera.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sedbed.o
sed_bedload_vandera.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sediment.o
sed_bedload_vandera.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
sed_bedload_vandera.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_vandera_funcs.o
sed_bedload_vandera.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o

sed_biodiff.o: set_bounds.h cppdefs.h ducknc.h globaldefs.h tile.h
sed_biodiff.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/bc_3d.o
sed_biodiff.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_bbl.o
sed_biodiff.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
sed_biodiff.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
sed_biodiff.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
sed_biodiff.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
sed_biodiff.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
sed_biodiff.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sedbed.o
sed_biodiff.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sediment.o
sed_biodiff.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
sed_biodiff.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o

sed_flocs.o: set_bounds.h cppdefs.h ducknc.h globaldefs.h tile.h
sed_flocs.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_bbl.o
sed_flocs.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
sed_flocs.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
sed_flocs.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_mixing.o
sed_flocs.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
sed_flocs.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
sed_flocs.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
sed_flocs.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sedflocs.o
sed_flocs.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sediment.o
sed_flocs.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o

sed_fluxes.o: tile.h cppdefs.h ducknc.h globaldefs.h set_bounds.h
sed_fluxes.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_bbl.o
sed_fluxes.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
sed_fluxes.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
sed_fluxes.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
sed_fluxes.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
sed_fluxes.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
sed_fluxes.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sedbed.o
sed_fluxes.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sediment.o
sed_fluxes.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o

sed_settling.o: cppdefs.h ducknc.h globaldefs.h set_bounds.h tile.h
sed_settling.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
sed_settling.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
sed_settling.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
sed_settling.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
sed_settling.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
sed_settling.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sedbed.o
sed_settling.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sediment.o
sed_settling.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o

sed_surface.o: cppdefs.h ducknc.h globaldefs.h set_bounds.h tile.h
sed_surface.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/bc_3d.o
sed_surface.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
sed_surface.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
sed_surface.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
sed_surface.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sedbed.o
sed_surface.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sediment.o
sed_surface.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
sed_surface.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o

sediment.o: cppdefs.h ducknc.h globaldefs.h
sediment.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_vandera_funcs.o
sediment.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/sed_bed.o
sediment.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/sed_bed2.o
sediment.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/sed_bed_cohesive.o
sediment.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/sed_bedload.o
sediment.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/sed_bedload_vandera.o
sediment.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/sed_biodiff.o
sediment.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/sed_flocs.o
sediment.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/sed_fluxes.o
sediment.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/sed_settling.o
sediment.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/sed_surface.o
sediment.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/sedtr_decay.o
sediment.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/sedtr_reactions_pom.o

sediment_output.o: cppdefs.h ducknc.h globaldefs.h
sediment_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_var.o
sediment_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/extract_sta.o
sediment_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_average.o
sediment_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_bbl.o
sediment_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
sediment_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
sediment_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
sediment_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_mixing.o
sediment_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
sediment_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
sediment_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
sediment_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
sediment_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
sediment_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
sediment_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
sediment_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sedbed.o
sediment_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sediment.o
sediment_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
sediment_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fwrite2d.o
sediment_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fwrite3d.o
sediment_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/omega.o
sediment_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

sedtr_decay.o: tile.h cppdefs.h ducknc.h globaldefs.h set_bounds.h
sedtr_decay.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/bc_3d.o
sedtr_decay.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_2d.o
sedtr_decay.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
sedtr_decay.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
sedtr_decay.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
sedtr_decay.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
sedtr_decay.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
sedtr_decay.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sedbed.o
sedtr_decay.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sediment.o
sedtr_decay.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
sedtr_decay.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o

sedtr_reactions_pom.o: set_bounds.h cppdefs.h ducknc.h globaldefs.h tile.h
sedtr_reactions_pom.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/bc_3d.o
sedtr_reactions_pom.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_2d.o
sedtr_reactions_pom.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_biology.o
sedtr_reactions_pom.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
sedtr_reactions_pom.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
sedtr_reactions_pom.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
sedtr_reactions_pom.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
sedtr_reactions_pom.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
sedtr_reactions_pom.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sedbed.o
sedtr_reactions_pom.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sediment.o
sedtr_reactions_pom.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
sedtr_reactions_pom.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o

sedtr_reactions_sed_decay.o: cppdefs.h ducknc.h globaldefs.h set_bounds.h
sedtr_reactions_sed_decay.o: tile.h
sedtr_reactions_sed_decay.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/bc_3d.o
sedtr_reactions_sed_decay.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_2d.o
sedtr_reactions_sed_decay.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
sedtr_reactions_sed_decay.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
sedtr_reactions_sed_decay.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
sedtr_reactions_sed_decay.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
sedtr_reactions_sed_decay.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
sedtr_reactions_sed_decay.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sedbed.o
sedtr_reactions_sed_decay.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sediment.o
sedtr_reactions_sed_decay.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
sedtr_reactions_sed_decay.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o

marsh_dynamics.o: cppdefs.h ducknc.h globaldefs.h
marsh_dynamics.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/marsh_sed_erosion.o
marsh_dynamics.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/marsh_tidal_range.o
marsh_dynamics.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/marsh_vert_growth.o
marsh_dynamics.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/marsh_wave_thrust.o
marsh_dynamics.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
marsh_dynamics.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
marsh_dynamics.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o

marsh_sed_erosion.o: set_bounds.h cppdefs.h ducknc.h globaldefs.h tile.h
marsh_sed_erosion.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/bc_2d.o
marsh_sed_erosion.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/bc_3d.o
marsh_sed_erosion.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/marsh_wave_thrust.o
marsh_sed_erosion.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
marsh_sed_erosion.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
marsh_sed_erosion.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
marsh_sed_erosion.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
marsh_sed_erosion.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
marsh_sed_erosion.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
marsh_sed_erosion.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sedbed.o
marsh_sed_erosion.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sediment.o
marsh_sed_erosion.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
marsh_sed_erosion.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_vegetation.o
marsh_sed_erosion.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o

marsh_tidal_range.o: tile.h set_bounds.h cppdefs.h ducknc.h globaldefs.h
marsh_tidal_range.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/bc_2d.o
marsh_tidal_range.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/bc_3d.o
marsh_tidal_range.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
marsh_tidal_range.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
marsh_tidal_range.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
marsh_tidal_range.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
marsh_tidal_range.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
marsh_tidal_range.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
marsh_tidal_range.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_vegetation.o
marsh_tidal_range.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o

marsh_vert_growth.o: set_bounds.h cppdefs.h ducknc.h globaldefs.h tile.h
marsh_vert_growth.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/bc_2d.o
marsh_vert_growth.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/bc_3d.o
marsh_vert_growth.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
marsh_vert_growth.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
marsh_vert_growth.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
marsh_vert_growth.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
marsh_vert_growth.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
marsh_vert_growth.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sedbed.o
marsh_vert_growth.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sediment.o
marsh_vert_growth.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
marsh_vert_growth.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_vegetation.o
marsh_vert_growth.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o

marsh_wave_thrust.o: tile.h set_bounds.h cppdefs.h ducknc.h globaldefs.h
marsh_wave_thrust.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/bc_2d.o
marsh_wave_thrust.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_2d.o
marsh_wave_thrust.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
marsh_wave_thrust.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
marsh_wave_thrust.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
marsh_wave_thrust.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
marsh_wave_thrust.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
marsh_wave_thrust.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
marsh_wave_thrust.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_vegetation.o
marsh_wave_thrust.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o

vegetation_biomass.o: tile.h cppdefs.h ducknc.h globaldefs.h set_bounds.h
vegetation_biomass.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
vegetation_biomass.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
vegetation_biomass.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
vegetation_biomass.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sedbed.o
vegetation_biomass.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sediment.o
vegetation_biomass.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o

vegetation_drag.o: tile.h set_bounds.h cppdefs.h ducknc.h globaldefs.h
vegetation_drag.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/bc_3d.o
vegetation_drag.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_2d.o
vegetation_drag.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_3d.o
vegetation_drag.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
vegetation_drag.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
vegetation_drag.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
vegetation_drag.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
vegetation_drag.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
vegetation_drag.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_vegetation.o
vegetation_drag.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o

vegetation_output.o: cppdefs.h ducknc.h globaldefs.h
vegetation_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_var.o
vegetation_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
vegetation_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
vegetation_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
vegetation_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
vegetation_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
vegetation_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
vegetation_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
vegetation_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
vegetation_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_vegetation.o
vegetation_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fwrite2d.o
vegetation_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fwrite3d.o
vegetation_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

vegetation_stream.o: set_bounds.h cppdefs.h ducknc.h globaldefs.h tile.h
vegetation_stream.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
vegetation_stream.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
vegetation_stream.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
vegetation_stream.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
vegetation_stream.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_vegetation.o

vegetation_turb.o: cppdefs.h ducknc.h globaldefs.h set_bounds.h tile.h
vegetation_turb.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
vegetation_turb.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_mixing.o
vegetation_turb.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
vegetation_turb.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
vegetation_turb.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
vegetation_turb.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
vegetation_turb.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_vegetation.o

wec_dissip.o: tile.h set_bounds.h cppdefs.h ducknc.h globaldefs.h
wec_dissip.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/bc_2d.o
wec_dissip.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_2d.o
wec_dissip.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_diags.o
wec_dissip.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
wec_dissip.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
wec_dissip.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
wec_dissip.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
wec_dissip.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
wec_dissip.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o

wec_output.o: cppdefs.h ducknc.h globaldefs.h
wec_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_var.o
wec_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/extract_sta.o
wec_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_average.o
wec_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
wec_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
wec_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
wec_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_mixing.o
wec_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
wec_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
wec_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
wec_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
wec_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
wec_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
wec_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
wec_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
wec_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fwrite2d.o
wec_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fwrite3d.o
wec_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/omega.o
wec_output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

wec_roller.o: tile.h cppdefs.h ducknc.h globaldefs.h set_bounds.h
wec_roller.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/bc_2d.o
wec_roller.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_diags.o
wec_roller.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
wec_roller.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
wec_roller.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
wec_roller.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
wec_roller.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
wec_roller.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
wec_roller.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o

wec_stokes.o: tile.h set_bounds.h cppdefs.h ducknc.h globaldefs.h
wec_stokes.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/bc_3d.o
wec_stokes.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_2d.o
wec_stokes.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_3d.o
wec_stokes.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_diags.o
wec_stokes.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
wec_stokes.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
wec_stokes.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
wec_stokes.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
wec_stokes.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
wec_stokes.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o
wec_stokes.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wec_u2dbc_im.o
wec_stokes.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wec_u3dbc_im.o
wec_stokes.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wec_v2dbc_im.o
wec_stokes.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wec_v3dbc_im.o

wec_streaming.o: tile.h set_bounds.h cppdefs.h ducknc.h globaldefs.h
wec_streaming.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/bc_2d.o
wec_streaming.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/bc_3d.o
wec_streaming.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_diags.o
wec_streaming.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
wec_streaming.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
wec_streaming.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_mixing.o
wec_streaming.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
wec_streaming.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
wec_streaming.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
wec_streaming.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_vegetation.o
wec_streaming.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o
wec_streaming.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/vegetation_stream.o

wec_u2dbc_im.o: cppdefs.h ducknc.h globaldefs.h set_bounds.h tile.h
wec_u2dbc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_boundary.o
wec_u2dbc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
wec_u2dbc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
wec_u2dbc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
wec_u2dbc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
wec_u2dbc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
wec_u2dbc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o

wec_u3dbc_im.o: set_bounds.h cppdefs.h ducknc.h globaldefs.h tile.h
wec_u3dbc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_boundary.o
wec_u3dbc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
wec_u3dbc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
wec_u3dbc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
wec_u3dbc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
wec_u3dbc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
wec_u3dbc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o

wec_v2dbc_im.o: cppdefs.h ducknc.h globaldefs.h set_bounds.h tile.h
wec_v2dbc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_boundary.o
wec_v2dbc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
wec_v2dbc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
wec_v2dbc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
wec_v2dbc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
wec_v2dbc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
wec_v2dbc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o

wec_v3dbc_im.o: tile.h cppdefs.h ducknc.h globaldefs.h set_bounds.h
wec_v3dbc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_boundary.o
wec_v3dbc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
wec_v3dbc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
wec_v3dbc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
wec_v3dbc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
wec_v3dbc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
wec_v3dbc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o

wec_vf.o: tile.h cppdefs.h ducknc.h globaldefs.h set_bounds.h
wec_vf.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/bc_2d.o
wec_vf.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/bc_3d.o
wec_vf.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_diags.o
wec_vf.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
wec_vf.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
wec_vf.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_mixing.o
wec_vf.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
wec_vf.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
wec_vf.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
wec_vf.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
wec_vf.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o

wec_wave_mix.o: set_bounds.h cppdefs.h ducknc.h globaldefs.h tile.h
wec_wave_mix.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_3d.o
wec_wave_mix.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_coupling.o
wec_wave_mix.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_diags.o
wec_wave_mix.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
wec_wave_mix.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
wec_wave_mix.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_mixing.o
wec_wave_mix.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
wec_wave_mix.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
wec_wave_mix.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
wec_wave_mix.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
wec_wave_mix.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o

wec_wvelocity.o: tile.h cppdefs.h ducknc.h globaldefs.h set_bounds.h
wec_wvelocity.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/bc_3d.o
wec_wvelocity.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_coupling.o
wec_wvelocity.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
wec_wvelocity.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
wec_wvelocity.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
wec_wvelocity.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
wec_wvelocity.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
wec_wvelocity.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o

bc_2d.o: cppdefs.h ducknc.h globaldefs.h set_bounds.h
bc_2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_2d.o
bc_2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_boundary.o
bc_2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
bc_2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
bc_2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
bc_2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o

bc_3d.o: set_bounds.h cppdefs.h ducknc.h globaldefs.h
bc_3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_3d.o
bc_3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_boundary.o
bc_3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
bc_3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
bc_3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
bc_3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o

bc_4d.o: set_bounds.h cppdefs.h ducknc.h globaldefs.h
bc_4d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_4d.o
bc_4d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_boundary.o
bc_4d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
bc_4d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
bc_4d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
bc_4d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o

bc_bry2d.o: set_bounds.h cppdefs.h ducknc.h globaldefs.h
bc_bry2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
bc_bry2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o

bc_bry3d.o: set_bounds.h cppdefs.h ducknc.h globaldefs.h
bc_bry3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
bc_bry3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o

bulk_flux.o: tile.h set_bounds.h cppdefs.h ducknc.h globaldefs.h
bulk_flux.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_2d.o
bulk_flux.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
bulk_flux.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
bulk_flux.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ice.o
bulk_flux.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_kinds.o
bulk_flux.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_mixing.o
bulk_flux.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
bulk_flux.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
bulk_flux.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
bulk_flux.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
bulk_flux.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o

bvf_mix.o: set_bounds.h cppdefs.h ducknc.h globaldefs.h tile.h
bvf_mix.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_3d.o
bvf_mix.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_mixing.o
bvf_mix.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
bvf_mix.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
bvf_mix.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o

conv_2d.o: cppdefs.h ducknc.h globaldefs.h set_bounds.h
conv_2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/bc_2d.o
conv_2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
conv_2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
conv_2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o

conv_3d.o: cppdefs.h ducknc.h globaldefs.h set_bounds.h
conv_3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/bc_3d.o
conv_3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
conv_3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
conv_3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o

conv_bry2d.o: cppdefs.h ducknc.h globaldefs.h set_bounds.h
conv_bry2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/bc_bry2d.o
conv_bry2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
conv_bry2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
conv_bry2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o

conv_bry3d.o: cppdefs.h ducknc.h globaldefs.h set_bounds.h
conv_bry3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/bc_bry3d.o
conv_bry3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
conv_bry3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
conv_bry3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o

diag.o: set_bounds.h cppdefs.h ducknc.h globaldefs.h tile.h
diag.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/analytical.o
diag.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/distribute.o
diag.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
diag.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
diag.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
diag.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
diag.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
diag.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
diag.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o

exchange_2d.o: set_bounds.h cppdefs.h ducknc.h globaldefs.h
exchange_2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
exchange_2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o

exchange_2d_xtr.o: cppdefs.h ducknc.h globaldefs.h set_bounds_xtr.h
exchange_2d_xtr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
exchange_2d_xtr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o

exchange_3d.o: cppdefs.h ducknc.h globaldefs.h set_bounds.h
exchange_3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
exchange_3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o

exchange_3d_xtr.o: set_bounds_xtr.h cppdefs.h ducknc.h globaldefs.h
exchange_3d_xtr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
exchange_3d_xtr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o

exchange_4d.o: set_bounds.h cppdefs.h ducknc.h globaldefs.h
exchange_4d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
exchange_4d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o

forcing.o: tile.h set_bounds.h cppdefs.h ducknc.h globaldefs.h
forcing.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_coupling.o
forcing.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
forcing.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
forcing.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
forcing.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
forcing.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o

frc_adjust.o: tile.h set_bounds.h cppdefs.h ducknc.h globaldefs.h
frc_adjust.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
frc_adjust.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
frc_adjust.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
frc_adjust.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o

get_data.o: cppdefs.h ducknc.h globaldefs.h
get_data.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_biology.o
get_data.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_boundary.o
get_data.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_clima.o
get_data.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
get_data.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
get_data.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
get_data.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
get_data.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
get_data.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
get_data.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
get_data.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sources.o
get_data.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
get_data.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

get_idata.o: cppdefs.h ducknc.h globaldefs.h
get_idata.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_tides.o
get_idata.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_biology.o
get_idata.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
get_idata.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
get_idata.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
get_idata.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
get_idata.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
get_idata.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
get_idata.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
get_idata.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
get_idata.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
get_idata.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sources.o
get_idata.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
get_idata.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_tides.o
get_idata.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fread3d.o
get_idata.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fread4d.o
get_idata.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

gls_corstep.o: tile.h set_bounds.h cppdefs.h ducknc.h globaldefs.h
gls_corstep.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_3d.o
gls_corstep.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
gls_corstep.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
gls_corstep.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_mixing.o
gls_corstep.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
gls_corstep.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
gls_corstep.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
gls_corstep.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
gls_corstep.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
gls_corstep.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_vegetation.o
gls_corstep.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o
gls_corstep.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/tkebc_im.o
gls_corstep.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/vegetation_turb.o

gls_prestep.o: tile.h set_bounds.h cppdefs.h ducknc.h globaldefs.h
gls_prestep.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_3d.o
gls_prestep.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
gls_prestep.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_mixing.o
gls_prestep.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
gls_prestep.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
gls_prestep.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
gls_prestep.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
gls_prestep.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o
gls_prestep.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/tkebc_im.o

hmixing.o: tile.h cppdefs.h ducknc.h globaldefs.h set_bounds.h
hmixing.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_3d.o
hmixing.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
hmixing.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_mixing.o
hmixing.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
hmixing.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
hmixing.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
hmixing.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
hmixing.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
hmixing.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_vegetation.o
hmixing.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o

ini_fields.o: cppdefs.h ducknc.h globaldefs.h set_bounds.h tile.h
ini_fields.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_2d.o
ini_fields.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_3d.o
ini_fields.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_coupling.o
ini_fields.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
ini_fields.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
ini_fields.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
ini_fields.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
ini_fields.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
ini_fields.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sedbed.o
ini_fields.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sediment.o
ini_fields.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
ini_fields.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o
ini_fields.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/t3dbc_im.o
ini_fields.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/u2dbc_im.o
ini_fields.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/u3dbc_im.o
ini_fields.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/uv_var_change.o
ini_fields.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/v2dbc_im.o
ini_fields.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/v3dbc_im.o
ini_fields.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/zetabc.o

initial.o: cppdefs.h ducknc.h globaldefs.h
initial.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/analytical.o
initial.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/close_io.o
initial.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/coupler.o
initial.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/dateclock.o
initial.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_ini.o
initial.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/distribute.o
initial.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/get_state.o
initial.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/get_wetdry.o
initial.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/ini_adjust.o
initial.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/ini_fields.o
initial.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/ini_hmixcoef.o
initial.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_bbl.o
initial.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_fourdvar.o
initial.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
initial.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
initial.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
initial.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_nesting.o
initial.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
initial.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
initial.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
initial.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
initial.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
initial.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nesting.o
initial.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/obs_initial.o
initial.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/omega.o
initial.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/rho_eos.o
initial.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/set_depth.o
initial.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/set_masks.o
initial.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/set_massflux.o
initial.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/stiffness.o
initial.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o
initial.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wetdry.o
initial.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wpoints.o

interp_floats.o: cppdefs.h ducknc.h globaldefs.h
interp_floats.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_floats.o
interp_floats.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
interp_floats.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
interp_floats.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o

lmd_bkpp.o: cppdefs.h ducknc.h globaldefs.h set_bounds.h tile.h
lmd_bkpp.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/bc_2d.o
lmd_bkpp.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
lmd_bkpp.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
lmd_bkpp.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_mixing.o
lmd_bkpp.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
lmd_bkpp.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
lmd_bkpp.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
lmd_bkpp.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
lmd_bkpp.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o
lmd_bkpp.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/shapiro.o

lmd_skpp.o: set_bounds.h cppdefs.h ducknc.h globaldefs.h tile.h
lmd_skpp.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/bc_2d.o
lmd_skpp.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
lmd_skpp.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
lmd_skpp.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_mixing.o
lmd_skpp.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
lmd_skpp.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
lmd_skpp.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
lmd_skpp.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
lmd_skpp.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o
lmd_skpp.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/shapiro.o

lmd_swfrac.o: set_bounds.h cppdefs.h ducknc.h globaldefs.h
lmd_swfrac.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_mixing.o
lmd_swfrac.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
lmd_swfrac.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o

lmd_vmix.o: tile.h cppdefs.h ducknc.h globaldefs.h set_bounds.h
lmd_vmix.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/bc_3d.o
lmd_vmix.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/lmd_bkpp.o
lmd_vmix.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/lmd_skpp.o
lmd_vmix.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
lmd_vmix.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_mixing.o
lmd_vmix.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
lmd_vmix.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
lmd_vmix.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
lmd_vmix.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
lmd_vmix.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o

main2d.o: cppdefs.h ducknc.h globaldefs.h
main2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/coupler.o
main2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/dateclock.o
main2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/diag.o
main2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/dotproduct.o
main2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/equilibrium_tide.o
main2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/forcing.o
main2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/frc_adjust.o
main2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/ini_fields.o
main2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_coupler.o
main2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
main2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_nesting.o
main2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
main2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
main2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
main2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
main2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nesting.o
main2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/obc_adjust.o
main2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/set_avg.o
main2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/set_tides.o
main2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/set_vbc.o
main2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/step2d.o
main2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/step_floats.o
main2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

main3d.o: cppdefs.h ducknc.h globaldefs.h
main3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/analytical.o
main3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/bbl.o
main3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/biology.o
main3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/bulk_flux.o
main3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/bvf_mix.o
main3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/coupler.o
main3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/dateclock.o
main3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/diag.o
main3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/dotproduct.o
main3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/equilibrium_tide.o
main3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/forcing.o
main3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/frc_adjust.o
main3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/frc_iau.o
main3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/gls_corstep.o
main3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/gls_prestep.o
main3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/hmixing.o
main3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/lmd_vmix.o
main3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/marsh_dynamics.o
main3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_coupler.o
main3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
main3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_nesting.o
main3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
main3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
main3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
main3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
main3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_vegetation.o
main3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/my25_corstep.o
main3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/my25_prestep.o
main3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nesting.o
main3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/obc_adjust.o
main3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/omega.o
main3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/post_initial.o
main3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/rho_eos.o
main3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/rhs3d.o
main3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/sediment.o
main3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/set_avg.o
main3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/set_depth.o
main3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/set_massflux.o
main3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/set_tides.o
main3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/set_vbc.o
main3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/set_zeta.o
main3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/step2d.o
main3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/step3d_t.o
main3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/step3d_uv.o
main3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/step_floats.o
main3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o
main3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wec_dissip.o
main3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wec_roller.o
main3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wec_stokes.o
main3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wec_vf.o
main3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wec_wave_mix.o
main3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wec_wvelocity.o
main3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wvelocity.o

mpdata_adiff.o: cppdefs.h ducknc.h globaldefs.h set_bounds.h
mpdata_adiff.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
mpdata_adiff.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
mpdata_adiff.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o

my25_corstep.o: cppdefs.h ducknc.h globaldefs.h set_bounds.h tile.h
my25_corstep.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_3d.o
my25_corstep.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
my25_corstep.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
my25_corstep.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_mixing.o
my25_corstep.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
my25_corstep.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
my25_corstep.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
my25_corstep.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
my25_corstep.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
my25_corstep.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o
my25_corstep.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/tkebc_im.o

my25_prestep.o: tile.h set_bounds.h cppdefs.h ducknc.h globaldefs.h
my25_prestep.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_3d.o
my25_prestep.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
my25_prestep.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_mixing.o
my25_prestep.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
my25_prestep.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
my25_prestep.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
my25_prestep.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
my25_prestep.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
my25_prestep.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o
my25_prestep.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/tkebc_im.o

nesting.o: tile.h cppdefs.h ducknc.h globaldefs.h set_bounds.h
nesting.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/distribute.o
nesting.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_2d.o
nesting.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_3d.o
nesting.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_clima.o
nesting.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_coupling.o
nesting.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
nesting.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
nesting.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
nesting.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
nesting.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_nesting.o
nesting.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
nesting.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
nesting.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
nesting.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
nesting.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
nesting.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o
nesting.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/set_depth.o
nesting.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

obc_adjust.o: set_bounds.h cppdefs.h ducknc.h globaldefs.h tile.h
obc_adjust.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_boundary.o
obc_adjust.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
obc_adjust.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
obc_adjust.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o

obc_volcons.o: cppdefs.h ducknc.h globaldefs.h set_bounds.h tile.h
obc_volcons.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/distribute.o
obc_volcons.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
obc_volcons.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
obc_volcons.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
obc_volcons.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
obc_volcons.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
obc_volcons.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o

omega.o: cppdefs.h ducknc.h globaldefs.h set_bounds.h tile.h
omega.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/bc_3d.o
omega.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_3d.o
omega.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
omega.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
omega.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
omega.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
omega.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
omega.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sedbed.o
omega.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sources.o
omega.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
omega.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o

output.o: cppdefs.h ducknc.h globaldefs.h
output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/close_io.o
output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_avg.o
output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_diags.o
output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_extract.o
output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_floats.o
output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_his.o
output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_quick.o
output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_rst.o
output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_station.o
output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/distribute.o
output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_floats.o
output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_fourdvar.o
output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/obs_read.o
output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/obs_write.o
output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o
output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wrt_avg.o
output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wrt_diags.o
output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wrt_extract.o
output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wrt_floats.o
output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wrt_his.o
output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wrt_quick.o
output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wrt_rst.o
output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wrt_station.o
output.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wrt_tides.o

post_initial.o: cppdefs.h ducknc.h globaldefs.h
post_initial.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/ini_fields.o
post_initial.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_nesting.o
post_initial.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
post_initial.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
post_initial.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
post_initial.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nesting.o
post_initial.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/set_depth.o

pre_step3d.o: cppdefs.h ducknc.h globaldefs.h set_bounds.h tile.h
pre_step3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_3d.o
pre_step3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_diags.o
pre_step3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
pre_step3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
pre_step3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_mixing.o
pre_step3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
pre_step3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
pre_step3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
pre_step3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sources.o
pre_step3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
pre_step3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o
pre_step3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/t3dbc_im.o

prsgrd.o: prsgrd40.h set_bounds.h tile.h prsgrd31.h prsgrd32.h cppdefs.h
prsgrd.o: ducknc.h globaldefs.h prsgrd42.h prsgrd44.h
prsgrd.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_diags.o
prsgrd.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
prsgrd.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
prsgrd.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
prsgrd.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
prsgrd.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
prsgrd.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o

rho_eos.o: cppdefs.h ducknc.h globaldefs.h set_bounds.h tile.h
rho_eos.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_2d.o
rho_eos.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_3d.o
rho_eos.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_coupling.o
rho_eos.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_eoscoef.o
rho_eos.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
rho_eos.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_mixing.o
rho_eos.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
rho_eos.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
rho_eos.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
rho_eos.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sediment.o
rho_eos.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
rho_eos.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o

rhs3d.o: tile.h set_bounds.h cppdefs.h ducknc.h globaldefs.h
rhs3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_clima.o
rhs3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_coupling.o
rhs3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_diags.o
rhs3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
rhs3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
rhs3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_mixing.o
rhs3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
rhs3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
rhs3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
rhs3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
rhs3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_vegetation.o
rhs3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/pre_step3d.o
rhs3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/prsgrd.o
rhs3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/t3dmix.o
rhs3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/t3dmix.o
rhs3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/uv3dmix.o
rhs3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/uv3dmix.o
rhs3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/vegetation_drag.o
rhs3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wec_streaming.o
rhs3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wec_vf.o

set_avg.o: tile.h set_bounds.h cppdefs.h ducknc.h globaldefs.h
set_avg.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_2d.o
set_avg.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_3d.o
set_avg.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_average.o
set_avg.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_bbl.o
set_avg.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_coupling.o
set_avg.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
set_avg.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
set_avg.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_mixing.o
set_avg.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
set_avg.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
set_avg.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
set_avg.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
set_avg.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sedbed.o
set_avg.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sediment.o
set_avg.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
set_avg.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_tides.o
set_avg.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o
set_avg.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/set_masks.o
set_avg.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/uv_rotate.o
set_avg.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/vorticity.o

set_data.o: tile.h set_bounds.h cppdefs.h ducknc.h globaldefs.h
set_data.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/analytical.o
set_data.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/distribute.o
set_data.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_2d.o
set_data.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_biology.o
set_data.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_boundary.o
set_data.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_clima.o
set_data.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
set_data.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
set_data.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_mixing.o
set_data.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
set_data.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
set_data.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
set_data.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
set_data.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sources.o
set_data.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
set_data.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o
set_data.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/set_2dfld.o
set_data.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/set_3dfld.o
set_data.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

set_depth.o: tile.h cppdefs.h ducknc.h globaldefs.h set_bounds.h
set_depth.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_2d.o
set_depth.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_3d.o
set_depth.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_boundary.o
set_depth.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_coupling.o
set_depth.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
set_depth.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
set_depth.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
set_depth.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
set_depth.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
set_depth.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
set_depth.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o

set_massflux.o: tile.h set_bounds.h cppdefs.h ducknc.h globaldefs.h
set_massflux.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_3d.o
set_massflux.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_coupling.o
set_massflux.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
set_massflux.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
set_massflux.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
set_massflux.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
set_massflux.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
set_massflux.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o

set_tides.o: cppdefs.h ducknc.h globaldefs.h set_bounds.h tile.h
set_tides.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/distribute.o
set_tides.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_2d.o
set_tides.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_boundary.o
set_tides.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_clima.o
set_tides.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
set_tides.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
set_tides.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
set_tides.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
set_tides.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
set_tides.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_tides.o
set_tides.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o

set_vbc.o: cppdefs.h ducknc.h globaldefs.h set_bounds.h tile.h
set_vbc.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/bc_2d.o
set_vbc.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
set_vbc.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
set_vbc.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
set_vbc.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
set_vbc.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
set_vbc.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
set_vbc.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o

set_zeta.o: tile.h set_bounds.h cppdefs.h ducknc.h globaldefs.h
set_zeta.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_2d.o
set_zeta.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_coupling.o
set_zeta.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
set_zeta.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
set_zeta.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
set_zeta.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o

step2d.o: step2d_LF_AM3.h set_bounds.h tile.h cppdefs.h ducknc.h globaldefs.h
step2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_2d.o
step2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_clima.o
step2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_coupling.o
step2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_diags.o
step2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
step2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
step2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_mixing.o
step2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
step2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
step2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
step2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
step2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sedbed.o
step2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sources.o
step2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
step2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_vegetation.o
step2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o
step2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/obc_volcons.o
step2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/u2dbc_im.o
step2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/v2dbc_im.o
step2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/vegetation_drag.o
step2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wetdry.o
step2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/zetabc.o

step3d_t.o: tile.h cppdefs.h ducknc.h globaldefs.h set_bounds.h
step3d_t.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_3d.o
step3d_t.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_clima.o
step3d_t.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_diags.o
step3d_t.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
step3d_t.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_mixing.o
step3d_t.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
step3d_t.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_nesting.o
step3d_t.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
step3d_t.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
step3d_t.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
step3d_t.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sedbed.o
step3d_t.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sources.o
step3d_t.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
step3d_t.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o
step3d_t.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mpdata_adiff.o
step3d_t.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nesting.o
step3d_t.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/t3dbc_im.o

step3d_uv.o: tile.h set_bounds.h cppdefs.h ducknc.h globaldefs.h
step3d_uv.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_2d.o
step3d_uv.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_3d.o
step3d_uv.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_coupling.o
step3d_uv.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_diags.o
step3d_uv.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
step3d_uv.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
step3d_uv.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_mixing.o
step3d_uv.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
step3d_uv.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
step3d_uv.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
step3d_uv.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
step3d_uv.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sources.o
step3d_uv.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
step3d_uv.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o
step3d_uv.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/u3dbc_im.o
step3d_uv.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/uv_var_change.o
step3d_uv.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/v3dbc_im.o

step_floats.o: cppdefs.h ducknc.h globaldefs.h
step_floats.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/biology_floats.o
step_floats.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/distribute.o
step_floats.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/interp_floats.o
step_floats.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_floats.o
step_floats.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
step_floats.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
step_floats.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
step_floats.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
step_floats.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
step_floats.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
step_floats.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
step_floats.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
step_floats.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/vwalk_floats.o

t3dbc_im.o: set_bounds.h cppdefs.h ducknc.h globaldefs.h tile.h
t3dbc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_boundary.o
t3dbc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_clima.o
t3dbc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
t3dbc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
t3dbc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
t3dbc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
t3dbc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
t3dbc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o

t3dmix.o: t3dmix4_s.h tile.h set_bounds.h t3dmix2_geo.h t3dmix2_iso.h cppdefs.h
t3dmix.o: ducknc.h globaldefs.h t3dmix4_iso.h t3dmix2_s.h t3dmix4_geo.h
t3dmix.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_clima.o
t3dmix.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_diags.o
t3dmix.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
t3dmix.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_mixing.o
t3dmix.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
t3dmix.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
t3dmix.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
t3dmix.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
t3dmix.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o

tkebc_im.o: tile.h cppdefs.h ducknc.h globaldefs.h set_bounds.h
tkebc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_boundary.o
tkebc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
tkebc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_mixing.o
tkebc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
tkebc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
tkebc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
tkebc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o

u2dbc_im.o: cppdefs.h ducknc.h globaldefs.h set_bounds.h tile.h
u2dbc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_boundary.o
u2dbc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_clima.o
u2dbc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
u2dbc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
u2dbc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
u2dbc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_nesting.o
u2dbc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
u2dbc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
u2dbc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
u2dbc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o

u3dbc_im.o: set_bounds.h cppdefs.h ducknc.h globaldefs.h tile.h
u3dbc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_boundary.o
u3dbc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_clima.o
u3dbc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
u3dbc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
u3dbc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
u3dbc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
u3dbc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
u3dbc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o

uv3dmix.o: cppdefs.h ducknc.h globaldefs.h uv3dmix2_geo.h tile.h set_bounds.h
uv3dmix.o: uv3dmix4_geo.h uv3dmix4_s.h uv3dmix2_s.h
uv3dmix.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_coupling.o
uv3dmix.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_diags.o
uv3dmix.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
uv3dmix.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_mixing.o
uv3dmix.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
uv3dmix.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
uv3dmix.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
uv3dmix.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
uv3dmix.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o

v2dbc_im.o: cppdefs.h ducknc.h globaldefs.h set_bounds.h tile.h
v2dbc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_boundary.o
v2dbc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_clima.o
v2dbc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
v2dbc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
v2dbc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
v2dbc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_nesting.o
v2dbc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
v2dbc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
v2dbc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
v2dbc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o

v3dbc_im.o: cppdefs.h ducknc.h globaldefs.h set_bounds.h tile.h
v3dbc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_boundary.o
v3dbc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_clima.o
v3dbc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
v3dbc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
v3dbc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
v3dbc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
v3dbc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
v3dbc_im.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o

vwalk_floats.o: cppdefs.h ducknc.h globaldefs.h
vwalk_floats.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/distribute.o
vwalk_floats.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/interp_floats.o
vwalk_floats.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_floats.o
vwalk_floats.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
vwalk_floats.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_mixing.o
vwalk_floats.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
vwalk_floats.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
vwalk_floats.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
vwalk_floats.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
vwalk_floats.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
vwalk_floats.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
vwalk_floats.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nrutil.o

wetdry.o: cppdefs.h ducknc.h globaldefs.h set_bounds.h tile.h
wetdry.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_2d.o
wetdry.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_coupling.o
wetdry.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
wetdry.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
wetdry.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
wetdry.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
wetdry.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
wetdry.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sources.o
wetdry.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o

wvelocity.o: tile.h set_bounds.h cppdefs.h ducknc.h globaldefs.h
wvelocity.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/bc_3d.o
wvelocity.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_2d.o
wvelocity.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_coupling.o
wvelocity.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
wvelocity.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
wvelocity.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
wvelocity.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
wvelocity.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
wvelocity.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
wvelocity.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o

zetabc.o: cppdefs.h ducknc.h globaldefs.h set_bounds.h tile.h
zetabc.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_boundary.o
zetabc.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
zetabc.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
zetabc.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
zetabc.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
zetabc.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
zetabc.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o

ADfromTL.o: set_bounds.h cppdefs.h ducknc.h globaldefs.h tile.h
ADfromTL.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_clima.o
ADfromTL.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
ADfromTL.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
ADfromTL.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
ADfromTL.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
ADfromTL.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
ADfromTL.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o

array_modes.o: cppdefs.h ducknc.h globaldefs.h
array_modes.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/distribute.o
array_modes.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_fourdvar.o
array_modes.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
array_modes.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
array_modes.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
array_modes.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o

back_cost.o: tile.h set_bounds.h cppdefs.h ducknc.h globaldefs.h
back_cost.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/distribute.o
back_cost.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_boundary.o
back_cost.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
back_cost.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_fourdvar.o
back_cost.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
back_cost.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
back_cost.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
back_cost.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
back_cost.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
back_cost.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o

background_std.o: tile.h cppdefs.h ducknc.h globaldefs.h set_bounds.h
background_std.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_2d.o
background_std.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_3d.o
background_std.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_coupling.o
background_std.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
background_std.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
background_std.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_kinds.o
background_std.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
background_std.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
background_std.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
background_std.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
background_std.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
background_std.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o
background_std.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/set_depth.o

cgradient.o: set_bounds.h cppdefs.h ducknc.h globaldefs.h tile.h
cgradient.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/distribute.o
cgradient.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/lapack_mod.o
cgradient.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_boundary.o
cgradient.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_coupling.o
cgradient.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
cgradient.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_fourdvar.o
cgradient.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
cgradient.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
cgradient.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
cgradient.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
cgradient.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
cgradient.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
cgradient.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
cgradient.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
cgradient.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
cgradient.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
cgradient.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fread2d.o
cgradient.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fread2d_bry.o
cgradient.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fread3d.o
cgradient.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fread3d_bry.o
cgradient.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/state_addition.o
cgradient.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/state_copy.o
cgradient.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/state_dotprod.o
cgradient.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/state_initialize.o
cgradient.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/state_read.o
cgradient.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/state_scale.o
cgradient.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o
cgradient.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wrt_hessian.o

check_multifile.o: cppdefs.h ducknc.h globaldefs.h
check_multifile.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/dateclock.o
check_multifile.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
check_multifile.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
check_multifile.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
check_multifile.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
check_multifile.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
check_multifile.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
check_multifile.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
check_multifile.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

checkadj.o: cppdefs.h ducknc.h globaldefs.h
checkadj.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
checkadj.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
checkadj.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
checkadj.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
checkadj.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
checkadj.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_strings.o
checkadj.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

checkdefs.o: cppdefs.h ducknc.h globaldefs.h
checkdefs.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
checkdefs.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
checkdefs.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
checkdefs.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
checkdefs.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_strings.o
checkdefs.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

checkerror.o: cppdefs.h ducknc.h globaldefs.h
checkerror.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
checkerror.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
checkerror.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
checkerror.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o

checkvars.o: cppdefs.h ducknc.h globaldefs.h
checkvars.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_biology.o
checkvars.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ice.o
checkvars.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
checkvars.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
checkvars.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
checkvars.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
checkvars.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
checkvars.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
checkvars.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
checkvars.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sediment.o
checkvars.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_vegetation.o
checkvars.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

close_io.o: cppdefs.h ducknc.h globaldefs.h
close_io.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/dateclock.o
close_io.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
close_io.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
close_io.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
close_io.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
close_io.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
close_io.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
close_io.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
close_io.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

comp_Jb0.o: tile.h set_bounds.h cppdefs.h ducknc.h globaldefs.h
comp_Jb0.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_boundary.o
comp_Jb0.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
comp_Jb0.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_fourdvar.o
comp_Jb0.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
comp_Jb0.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
comp_Jb0.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
comp_Jb0.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
comp_Jb0.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
comp_Jb0.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
comp_Jb0.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
comp_Jb0.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
comp_Jb0.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/state_dotprod.o
comp_Jb0.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

congrad.o: cppdefs.h ducknc.h globaldefs.h
congrad.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/distribute.o
congrad.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/lapack_mod.o
congrad.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_fourdvar.o
congrad.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
congrad.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
congrad.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
congrad.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
congrad.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
congrad.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
congrad.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
congrad.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

convolve.o: cppdefs.h ducknc.h globaldefs.h
convolve.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/comp_Jb0.o
convolve.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/get_state.o
convolve.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/ini_adjust.o
convolve.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_boundary.o
convolve.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
convolve.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_fourdvar.o
convolve.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
convolve.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
convolve.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
convolve.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
convolve.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
convolve.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
convolve.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o
convolve.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/sum_grad.o
convolve.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/time_corr.o
convolve.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wrt_hessian.o
convolve.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wrt_ini.o

cost_grad.o: tile.h cppdefs.h ducknc.h globaldefs.h set_bounds.h
cost_grad.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_boundary.o
cost_grad.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
cost_grad.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
cost_grad.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
cost_grad.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
cost_grad.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o

dateclock.o: cppdefs.h ducknc.h globaldefs.h
dateclock.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_kinds.o
dateclock.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
dateclock.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
dateclock.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/round.o

def_avg.o: cppdefs.h ducknc.h globaldefs.h
def_avg.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/bbl_output.o
def_avg.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_dim.o
def_avg.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_info.o
def_avg.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_var.o
def_avg.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_biology.o
def_avg.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_fourdvar.o
def_avg.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
def_avg.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
def_avg.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
def_avg.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
def_avg.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
def_avg.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
def_avg.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
def_avg.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sediment.o
def_avg.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/sediment_output.o
def_avg.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o
def_avg.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wec_output.o
def_avg.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wrt_info.o

def_dai.o: cppdefs.h ducknc.h globaldefs.h
def_dai.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_dim.o
def_dai.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_info.o
def_dai.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_var.o
def_dai.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_biology.o
def_dai.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_fourdvar.o
def_dai.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
def_dai.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
def_dai.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
def_dai.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
def_dai.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
def_dai.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
def_dai.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
def_dai.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sediment.o
def_dai.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o
def_dai.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wrt_info.o

def_diags.o: cppdefs.h ducknc.h globaldefs.h
def_diags.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_dim.o
def_diags.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_info.o
def_diags.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_var.o
def_diags.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_biology.o
def_diags.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_fourdvar.o
def_diags.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
def_diags.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
def_diags.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
def_diags.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
def_diags.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
def_diags.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
def_diags.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
def_diags.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sediment.o
def_diags.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o
def_diags.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wrt_info.o

def_dim.o: cppdefs.h ducknc.h globaldefs.h
def_dim.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/distribute.o
def_dim.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
def_dim.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
def_dim.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
def_dim.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
def_dim.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
def_dim.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

def_error.o: cppdefs.h ducknc.h globaldefs.h
def_error.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_dim.o
def_error.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_info.o
def_error.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_var.o
def_error.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_biology.o
def_error.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_fourdvar.o
def_error.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
def_error.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
def_error.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
def_error.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
def_error.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
def_error.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
def_error.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
def_error.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sediment.o
def_error.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o
def_error.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wrt_info.o

def_extract.o: cppdefs.h ducknc.h globaldefs.h
def_extract.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/bbl_output.o
def_extract.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_dim.o
def_extract.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_info.o
def_extract.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_var.o
def_extract.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_biology.o
def_extract.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_fourdvar.o
def_extract.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
def_extract.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
def_extract.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
def_extract.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
def_extract.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
def_extract.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
def_extract.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
def_extract.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sediment.o
def_extract.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/sediment_output.o
def_extract.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o
def_extract.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wec_output.o
def_extract.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wrt_info.o

def_floats.o: cppdefs.h ducknc.h globaldefs.h
def_floats.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_dim.o
def_floats.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_info.o
def_floats.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_var.o
def_floats.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_biology.o
def_floats.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_floats.o
def_floats.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_fourdvar.o
def_floats.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
def_floats.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
def_floats.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
def_floats.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
def_floats.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
def_floats.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
def_floats.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
def_floats.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
def_floats.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sediment.o
def_floats.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o
def_floats.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wrt_info.o

def_gst.o: cppdefs.h ducknc.h globaldefs.h
def_gst.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_dim.o
def_gst.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_info.o
def_gst.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_var.o
def_gst.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
def_gst.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
def_gst.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
def_gst.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
def_gst.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
def_gst.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
def_gst.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
def_gst.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_storage.o
def_gst.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

def_hessian.o: cppdefs.h ducknc.h globaldefs.h
def_hessian.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_dim.o
def_hessian.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_info.o
def_hessian.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_var.o
def_hessian.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_biology.o
def_hessian.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_fourdvar.o
def_hessian.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
def_hessian.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
def_hessian.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
def_hessian.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
def_hessian.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
def_hessian.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
def_hessian.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
def_hessian.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sediment.o
def_hessian.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o
def_hessian.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wrt_info.o

def_his.o: cppdefs.h ducknc.h globaldefs.h
def_his.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/bbl_output.o
def_his.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_dim.o
def_his.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_info.o
def_his.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_var.o
def_his.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_biology.o
def_his.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_fourdvar.o
def_his.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
def_his.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
def_his.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
def_his.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
def_his.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
def_his.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
def_his.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
def_his.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sediment.o
def_his.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_vegetation.o
def_his.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/sediment_output.o
def_his.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o
def_his.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/vegetation_output.o
def_his.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wec_output.o
def_his.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wrt_info.o

def_impulse.o: cppdefs.h ducknc.h globaldefs.h
def_impulse.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_dim.o
def_impulse.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_info.o
def_impulse.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_var.o
def_impulse.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_biology.o
def_impulse.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_fourdvar.o
def_impulse.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
def_impulse.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
def_impulse.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
def_impulse.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
def_impulse.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
def_impulse.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
def_impulse.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
def_impulse.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sediment.o
def_impulse.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o
def_impulse.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wrt_info.o

def_info.o: npzd_Powell_def.h fennel_def.h npzd_iron_def_pio.h ecosim_def_pio.h
def_info.o: cppdefs.h ducknc.h globaldefs.h hypoxia_srm_def_pio.h
def_info.o: sediment_def.h nemuro_def_pio.h sediment_def_pio.h
def_info.o: npzd_Franks_def_pio.h oyster_floats_def.h npzd_iron_def.h
def_info.o: oyster_floats_def_pio.h red_tide_def.h fennel_def_pio.h
def_info.o: red_tide_def_pio.h ecosim_def.h npzd_Franks_def.h
def_info.o: npzd_Powell_def_pio.h hypoxia_srm_def.h nemuro_def.h
def_info.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_dim.o
def_info.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_var.o
def_info.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/distribute.o
def_info.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/lbc.o
def_info.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_fourdvar.o
def_info.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
def_info.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
def_info.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
def_info.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
def_info.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
def_info.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
def_info.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
def_info.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
def_info.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_strings.o
def_info.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o
def_info.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/tadv.o

def_ini.o: cppdefs.h ducknc.h globaldefs.h
def_ini.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_dim.o
def_ini.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_var.o
def_ini.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
def_ini.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
def_ini.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
def_ini.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
def_ini.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
def_ini.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
def_ini.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
def_ini.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

def_lanczos.o: cppdefs.h ducknc.h globaldefs.h
def_lanczos.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_dim.o
def_lanczos.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_info.o
def_lanczos.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_var.o
def_lanczos.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_biology.o
def_lanczos.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_fourdvar.o
def_lanczos.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
def_lanczos.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
def_lanczos.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
def_lanczos.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
def_lanczos.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
def_lanczos.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
def_lanczos.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
def_lanczos.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sediment.o
def_lanczos.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o
def_lanczos.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wrt_info.o

def_mod.o: cppdefs.h ducknc.h globaldefs.h
def_mod.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/dateclock.o
def_mod.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_dim.o
def_mod.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_var.o
def_mod.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/distribute.o
def_mod.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_fourdvar.o
def_mod.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
def_mod.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
def_mod.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
def_mod.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
def_mod.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
def_mod.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
def_mod.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
def_mod.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_strings.o
def_mod.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

def_norm.o: cppdefs.h ducknc.h globaldefs.h
def_norm.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_dim.o
def_norm.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_info.o
def_norm.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_var.o
def_norm.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_biology.o
def_norm.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_fourdvar.o
def_norm.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
def_norm.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
def_norm.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
def_norm.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
def_norm.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
def_norm.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
def_norm.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
def_norm.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sediment.o
def_norm.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o
def_norm.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wrt_info.o

def_quick.o: cppdefs.h ducknc.h globaldefs.h
def_quick.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/bbl_output.o
def_quick.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_dim.o
def_quick.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_info.o
def_quick.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_var.o
def_quick.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_biology.o
def_quick.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_fourdvar.o
def_quick.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
def_quick.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
def_quick.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
def_quick.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
def_quick.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
def_quick.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
def_quick.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
def_quick.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sediment.o
def_quick.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_vegetation.o
def_quick.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/sediment_output.o
def_quick.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o
def_quick.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/vegetation_output.o
def_quick.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wec_output.o
def_quick.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wrt_info.o

def_rst.o: cppdefs.h ducknc.h globaldefs.h
def_rst.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_dim.o
def_rst.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_info.o
def_rst.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_var.o
def_rst.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_biology.o
def_rst.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_fourdvar.o
def_rst.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
def_rst.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
def_rst.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
def_rst.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
def_rst.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
def_rst.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
def_rst.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
def_rst.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sediment.o
def_rst.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_vegetation.o
def_rst.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o
def_rst.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/vegetation_output.o
def_rst.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wrt_info.o

def_state.o: cppdefs.h ducknc.h globaldefs.h
def_state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_dim.o
def_state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_info.o
def_state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_var.o
def_state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_biology.o
def_state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_fourdvar.o
def_state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
def_state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
def_state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
def_state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
def_state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
def_state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
def_state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
def_state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sediment.o
def_state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o
def_state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wrt_info.o

def_station.o: cppdefs.h ducknc.h globaldefs.h
def_station.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/bbl_output.o
def_station.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_dim.o
def_station.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_info.o
def_station.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_var.o
def_station.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_biology.o
def_station.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_fourdvar.o
def_station.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
def_station.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
def_station.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
def_station.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
def_station.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
def_station.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
def_station.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
def_station.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sediment.o
def_station.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/sediment_output.o
def_station.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o
def_station.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wec_output.o
def_station.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wrt_info.o

def_std.o: cppdefs.h ducknc.h globaldefs.h
def_std.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_dim.o
def_std.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_info.o
def_std.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_var.o
def_std.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_biology.o
def_std.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_fourdvar.o
def_std.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
def_std.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
def_std.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
def_std.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
def_std.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
def_std.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
def_std.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
def_std.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sediment.o
def_std.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o
def_std.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wrt_info.o

def_tides.o: cppdefs.h ducknc.h globaldefs.h
def_tides.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_dim.o
def_tides.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_info.o
def_tides.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_var.o
def_tides.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
def_tides.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
def_tides.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
def_tides.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
def_tides.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
def_tides.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
def_tides.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
def_tides.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
def_tides.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_tides.o
def_tides.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o
def_tides.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wrt_info.o

def_var.o: cppdefs.h ducknc.h globaldefs.h
def_var.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/distribute.o
def_var.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
def_var.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
def_var.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
def_var.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
def_var.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
def_var.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
def_var.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

destroy.o: cppdefs.h ducknc.h globaldefs.h
destroy.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
destroy.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_kinds.o
destroy.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o

distribute.o: cppdefs.h ducknc.h globaldefs.h set_bounds.h
distribute.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
distribute.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
distribute.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
distribute.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
distribute.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o

dotproduct.o: tile.h cppdefs.h ducknc.h globaldefs.h set_bounds.h
dotproduct.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/distribute.o
dotproduct.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_fourdvar.o
dotproduct.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
dotproduct.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
dotproduct.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
dotproduct.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
dotproduct.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
dotproduct.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
dotproduct.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
dotproduct.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o

edit_multifile.o: cppdefs.h ducknc.h globaldefs.h
edit_multifile.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/close_io.o
edit_multifile.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
edit_multifile.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
edit_multifile.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
edit_multifile.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o

equilibrium_tide.o: tile.h cppdefs.h ducknc.h globaldefs.h set_bounds.h
equilibrium_tide.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/dateclock.o
equilibrium_tide.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_2d.o
equilibrium_tide.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
equilibrium_tide.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_kinds.o
equilibrium_tide.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
equilibrium_tide.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
equilibrium_tide.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
equilibrium_tide.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o

erf.o: cppdefs.h ducknc.h globaldefs.h
erf.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
erf.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_kinds.o
erf.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
erf.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o

extract_field.o: set_bounds.h cppdefs.h ducknc.h globaldefs.h
extract_field.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/distribute.o
extract_field.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/interpolate.o
extract_field.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_extract.o
extract_field.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
extract_field.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
extract_field.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
extract_field.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
extract_field.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
extract_field.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
extract_field.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

extract_obs.o: cppdefs.h ducknc.h globaldefs.h
extract_obs.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_fourdvar.o
extract_obs.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_kinds.o
extract_obs.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
extract_obs.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o

extract_slice.o: cppdefs.h ducknc.h globaldefs.h
extract_slice.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
extract_slice.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o

extract_sta.o: cppdefs.h ducknc.h globaldefs.h
extract_sta.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/distribute.o
extract_sta.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
extract_sta.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
extract_sta.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
extract_sta.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
extract_sta.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o

frc_iau.o: cppdefs.h ducknc.h globaldefs.h set_bounds.h tile.h
frc_iau.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
frc_iau.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_fourdvar.o
frc_iau.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
frc_iau.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
frc_iau.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
frc_iau.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
frc_iau.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o

frc_weak.o: set_bounds.h cppdefs.h ducknc.h globaldefs.h tile.h
frc_weak.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
frc_weak.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_fourdvar.o
frc_weak.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
frc_weak.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
frc_weak.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
frc_weak.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o

gasdev.o: cppdefs.h ducknc.h globaldefs.h
gasdev.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_kinds.o
gasdev.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nrutil.o

get_2dfld.o: cppdefs.h ducknc.h globaldefs.h
get_2dfld.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/dateclock.o
get_2dfld.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/inquiry.o
get_2dfld.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
get_2dfld.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
get_2dfld.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
get_2dfld.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
get_2dfld.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
get_2dfld.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
get_2dfld.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
get_2dfld.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fread2d.o
get_2dfld.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fread3d.o
get_2dfld.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

get_2dfldr.o: cppdefs.h ducknc.h globaldefs.h
get_2dfldr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/dateclock.o
get_2dfldr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/inquiry.o
get_2dfldr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
get_2dfldr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
get_2dfldr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
get_2dfldr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
get_2dfldr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
get_2dfldr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
get_2dfldr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
get_2dfldr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fread2d.o
get_2dfldr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fread3d.o
get_2dfldr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

get_3dfld.o: cppdefs.h ducknc.h globaldefs.h
get_3dfld.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/dateclock.o
get_3dfld.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/inquiry.o
get_3dfld.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
get_3dfld.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
get_3dfld.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
get_3dfld.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
get_3dfld.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
get_3dfld.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
get_3dfld.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
get_3dfld.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fread3d.o
get_3dfld.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

get_3dfldr.o: cppdefs.h ducknc.h globaldefs.h
get_3dfldr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/dateclock.o
get_3dfldr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/inquiry.o
get_3dfldr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
get_3dfldr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
get_3dfldr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
get_3dfldr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
get_3dfldr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
get_3dfldr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
get_3dfldr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
get_3dfldr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fread3d.o
get_3dfldr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

get_bounds.o: cppdefs.h ducknc.h globaldefs.h
get_bounds.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
get_bounds.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_nesting.o
get_bounds.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
get_bounds.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
get_bounds.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o

get_cycle.o: cppdefs.h ducknc.h globaldefs.h
get_cycle.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
get_cycle.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
get_cycle.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
get_cycle.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
get_cycle.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
get_cycle.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
get_cycle.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
get_cycle.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

get_env.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
get_env.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o

get_extract.o: cppdefs.h ducknc.h globaldefs.h
get_extract.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_2d_xtr.o
get_extract.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/extract_field.o
get_extract.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_extract.o
get_extract.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
get_extract.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
get_extract.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
get_extract.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
get_extract.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
get_extract.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
get_extract.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
get_extract.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o
get_extract.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fread2d_xtr.o
get_extract.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

get_grid.o: cppdefs.h ducknc.h globaldefs.h
get_grid.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_2d.o
get_grid.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
get_grid.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
get_grid.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_mixing.o
get_grid.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
get_grid.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_nesting.o
get_grid.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
get_grid.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
get_grid.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
get_grid.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
get_grid.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
get_grid.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o
get_grid.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nesting.o
get_grid.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fread2d.o
get_grid.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

get_gst.o: cppdefs.h ducknc.h globaldefs.h
get_gst.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/distribute.o
get_gst.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
get_gst.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
get_gst.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
get_gst.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
get_gst.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
get_gst.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
get_gst.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
get_gst.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_storage.o
get_gst.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

get_hash.o: cppdefs.h ducknc.h globaldefs.h
get_hash.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/distribute.o
get_hash.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
get_hash.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_kinds.o
get_hash.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
get_hash.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

get_metadata.o: cppdefs.h ducknc.h globaldefs.h
get_metadata.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
get_metadata.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_kinds.o
get_metadata.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
get_metadata.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
get_metadata.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o
get_metadata.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/yaml_parser.o

get_ngfld.o: cppdefs.h ducknc.h globaldefs.h
get_ngfld.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/dateclock.o
get_ngfld.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/get_hash.o
get_ngfld.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/inquiry.o
get_ngfld.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
get_ngfld.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
get_ngfld.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
get_ngfld.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
get_ngfld.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
get_ngfld.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
get_ngfld.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
get_ngfld.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

get_ngfldr.o: cppdefs.h ducknc.h globaldefs.h
get_ngfldr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/dateclock.o
get_ngfldr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/get_hash.o
get_ngfldr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/inquiry.o
get_ngfldr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
get_ngfldr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
get_ngfldr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
get_ngfldr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
get_ngfldr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
get_ngfldr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
get_ngfldr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
get_ngfldr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

get_nudgcoef.o: cppdefs.h ducknc.h globaldefs.h
get_nudgcoef.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_2d.o
get_nudgcoef.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_3d.o
get_nudgcoef.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_clima.o
get_nudgcoef.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
get_nudgcoef.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
get_nudgcoef.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
get_nudgcoef.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
get_nudgcoef.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
get_nudgcoef.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
get_nudgcoef.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
get_nudgcoef.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
get_nudgcoef.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o
get_nudgcoef.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fread2d.o
get_nudgcoef.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fread3d.o
get_nudgcoef.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

get_state.o: cppdefs.h ducknc.h globaldefs.h
get_state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/checkvars.o
get_state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/dateclock.o
get_state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/lbc.o
get_state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_boundary.o
get_state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
get_state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_fourdvar.o
get_state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
get_state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ice.o
get_state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
get_state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_mixing.o
get_state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
get_state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
get_state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
get_state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
get_state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
get_state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
get_state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
get_state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sedbed.o
get_state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sediment.o
get_state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
get_state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_strings.o
get_state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_vegetation.o
get_state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o
get_state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fread2d.o
get_state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fread2d_bry.o
get_state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fread3d.o
get_state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fread3d_bry.o
get_state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fread4d.o
get_state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

get_varcoords.o: cppdefs.h ducknc.h globaldefs.h
get_varcoords.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
get_varcoords.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
get_varcoords.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
get_varcoords.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
get_varcoords.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
get_varcoords.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
get_varcoords.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
get_varcoords.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

get_wetdry.o: cppdefs.h ducknc.h globaldefs.h
get_wetdry.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_2d.o
get_wetdry.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
get_wetdry.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
get_wetdry.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
get_wetdry.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
get_wetdry.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
get_wetdry.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
get_wetdry.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
get_wetdry.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
get_wetdry.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o
get_wetdry.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fread2d.o
get_wetdry.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

grid_coords.o: cppdefs.h ducknc.h globaldefs.h
grid_coords.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/distribute.o
grid_coords.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/interpolate.o
grid_coords.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_floats.o
grid_coords.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
grid_coords.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
grid_coords.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
grid_coords.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o

ini_adjust.o: tile.h cppdefs.h ducknc.h globaldefs.h set_bounds.h
ini_adjust.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_2d.o
ini_adjust.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_3d.o
ini_adjust.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_boundary.o
ini_adjust.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_coupling.o
ini_adjust.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
ini_adjust.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_fourdvar.o
ini_adjust.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
ini_adjust.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
ini_adjust.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
ini_adjust.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
ini_adjust.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
ini_adjust.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
ini_adjust.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
ini_adjust.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sedbed.o
ini_adjust.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
ini_adjust.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o
ini_adjust.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/set_depth.o
ini_adjust.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/state_addition.o
ini_adjust.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/state_assign.o
ini_adjust.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/state_copy.o
ini_adjust.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/t3dbc_im.o
ini_adjust.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/u2dbc_im.o
ini_adjust.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/u3dbc_im.o
ini_adjust.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/v2dbc_im.o
ini_adjust.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/v3dbc_im.o
ini_adjust.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/zetabc.o

ini_hmixcoef.o: cppdefs.h ducknc.h globaldefs.h set_bounds.h tile.h
ini_hmixcoef.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_2d.o
ini_hmixcoef.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
ini_hmixcoef.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_mixing.o
ini_hmixcoef.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
ini_hmixcoef.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
ini_hmixcoef.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
ini_hmixcoef.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o

ini_lanczos.o: set_bounds.h cppdefs.h ducknc.h globaldefs.h tile.h
ini_lanczos.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_boundary.o
ini_lanczos.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
ini_lanczos.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_fourdvar.o
ini_lanczos.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
ini_lanczos.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
ini_lanczos.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
ini_lanczos.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
ini_lanczos.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
ini_lanczos.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
ini_lanczos.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
ini_lanczos.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
ini_lanczos.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
ini_lanczos.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/state_addition.o
ini_lanczos.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/state_dotprod.o
ini_lanczos.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/state_initialize.o
ini_lanczos.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/state_read.o
ini_lanczos.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/state_scale.o
ini_lanczos.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

inner2state.o: tile.h set_bounds.h cppdefs.h ducknc.h globaldefs.h
inner2state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/distribute.o
inner2state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_boundary.o
inner2state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
inner2state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_fourdvar.o
inner2state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
inner2state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
inner2state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
inner2state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
inner2state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
inner2state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
inner2state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
inner2state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
inner2state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
inner2state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
inner2state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/state_addition.o
inner2state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/state_copy.o
inner2state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/state_dotprod.o
inner2state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/state_initialize.o
inner2state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/state_read.o
inner2state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/state_scale.o
inner2state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

inp_decode.o: cppdefs.h ducknc.h globaldefs.h
inp_decode.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
inp_decode.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_kinds.o
inp_decode.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
inp_decode.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
inp_decode.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
inp_decode.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
inp_decode.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
inp_decode.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

inp_par.o: cppdefs.h ducknc.h globaldefs.h
inp_par.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/dateclock.o
inp_par.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/distribute.o
inp_par.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/lbc.o
inp_par.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
inp_par.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_kinds.o
inp_par.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
inp_par.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
inp_par.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
inp_par.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
inp_par.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
inp_par.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
inp_par.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_strings.o
inp_par.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/ran_state.o
inp_par.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/set_contact.o
inp_par.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/stdinp_mod.o
inp_par.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o
inp_par.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/tadv.o
inp_par.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/tile_indices.o

inquiry.o: cppdefs.h ducknc.h globaldefs.h
inquiry.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/get_cycle.o
inquiry.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
inquiry.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
inquiry.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
inquiry.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
inquiry.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
inquiry.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
inquiry.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
inquiry.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

interpolate.o: cppdefs.h ducknc.h globaldefs.h
interpolate.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/distribute.o
interpolate.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
interpolate.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
interpolate.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
interpolate.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
interpolate.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
interpolate.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
interpolate.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

lanc_resid.o: cppdefs.h ducknc.h globaldefs.h set_bounds.h tile.h
lanc_resid.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_boundary.o
lanc_resid.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
lanc_resid.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
lanc_resid.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
lanc_resid.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
lanc_resid.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o

lapack_mod.o: cppdefs.h ducknc.h globaldefs.h
lapack_mod.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_kinds.o

lbc.o: cppdefs.h ducknc.h globaldefs.h
lbc.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/distribute.o
lbc.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
lbc.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
lbc.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
lbc.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
lbc.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
lbc.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
lbc.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

lubksb.o: cppdefs.h ducknc.h globaldefs.h
lubksb.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_kinds.o

ludcmp.o: cppdefs.h ducknc.h globaldefs.h
ludcmp.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_kinds.o

memory.o: cppdefs.h ducknc.h globaldefs.h
memory.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/distribute.o
memory.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
memory.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
memory.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
memory.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
memory.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o

metrics.o: tile.h cppdefs.h ducknc.h globaldefs.h set_bounds.h
metrics.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/distribute.o
metrics.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_2d.o
metrics.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_fourdvar.o
metrics.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
metrics.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
metrics.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_mixing.o
metrics.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
metrics.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_nesting.o
metrics.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
metrics.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
metrics.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
metrics.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
metrics.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
metrics.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o
metrics.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/set_depth.o

mp_exchange.o: set_bounds.h set_bounds_xtr.h cppdefs.h ducknc.h globaldefs.h
mp_exchange.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
mp_exchange.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
mp_exchange.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
mp_exchange.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o

mp_routines.o: cppdefs.h ducknc.h globaldefs.h
mp_routines.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_kinds.o

nf_fread2d.o: cppdefs.h ducknc.h globaldefs.h
nf_fread2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/distribute.o
nf_fread2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/get_bounds.o
nf_fread2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/get_hash.o
nf_fread2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
nf_fread2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
nf_fread2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
nf_fread2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
nf_fread2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
nf_fread2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
nf_fread2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
nf_fread2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
nf_fread2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/regrid.o
nf_fread2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

nf_fread2d_bry.o: cppdefs.h ducknc.h globaldefs.h
nf_fread2d_bry.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/distribute.o
nf_fread2d_bry.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/get_hash.o
nf_fread2d_bry.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
nf_fread2d_bry.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
nf_fread2d_bry.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
nf_fread2d_bry.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
nf_fread2d_bry.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
nf_fread2d_bry.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
nf_fread2d_bry.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
nf_fread2d_bry.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

nf_fread2d_xtr.o: cppdefs.h ducknc.h globaldefs.h
nf_fread2d_xtr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/distribute.o
nf_fread2d_xtr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/get_bounds.o
nf_fread2d_xtr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/get_hash.o
nf_fread2d_xtr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
nf_fread2d_xtr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
nf_fread2d_xtr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
nf_fread2d_xtr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
nf_fread2d_xtr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
nf_fread2d_xtr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
nf_fread2d_xtr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
nf_fread2d_xtr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
nf_fread2d_xtr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/regrid.o
nf_fread2d_xtr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

nf_fread3d.o: cppdefs.h ducknc.h globaldefs.h
nf_fread3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/distribute.o
nf_fread3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/get_bounds.o
nf_fread3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/get_hash.o
nf_fread3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
nf_fread3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
nf_fread3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
nf_fread3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
nf_fread3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
nf_fread3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
nf_fread3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
nf_fread3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

nf_fread3d_bry.o: cppdefs.h ducknc.h globaldefs.h
nf_fread3d_bry.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/distribute.o
nf_fread3d_bry.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/get_hash.o
nf_fread3d_bry.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
nf_fread3d_bry.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
nf_fread3d_bry.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
nf_fread3d_bry.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
nf_fread3d_bry.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
nf_fread3d_bry.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
nf_fread3d_bry.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
nf_fread3d_bry.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

nf_fread3d_xtr.o: cppdefs.h ducknc.h globaldefs.h
nf_fread3d_xtr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/distribute.o
nf_fread3d_xtr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/get_bounds.o
nf_fread3d_xtr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/get_hash.o
nf_fread3d_xtr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
nf_fread3d_xtr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
nf_fread3d_xtr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
nf_fread3d_xtr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
nf_fread3d_xtr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
nf_fread3d_xtr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
nf_fread3d_xtr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
nf_fread3d_xtr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

nf_fread4d.o: cppdefs.h ducknc.h globaldefs.h
nf_fread4d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/distribute.o
nf_fread4d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/get_bounds.o
nf_fread4d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/get_hash.o
nf_fread4d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
nf_fread4d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
nf_fread4d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
nf_fread4d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
nf_fread4d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
nf_fread4d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
nf_fread4d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
nf_fread4d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

nf_fwrite2d.o: cppdefs.h ducknc.h globaldefs.h
nf_fwrite2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/distribute.o
nf_fwrite2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/get_bounds.o
nf_fwrite2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
nf_fwrite2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
nf_fwrite2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
nf_fwrite2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
nf_fwrite2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
nf_fwrite2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
nf_fwrite2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/pack_field.o
nf_fwrite2d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/stats.o

nf_fwrite2d_bry.o: cppdefs.h ducknc.h globaldefs.h
nf_fwrite2d_bry.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/distribute.o
nf_fwrite2d_bry.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
nf_fwrite2d_bry.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
nf_fwrite2d_bry.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
nf_fwrite2d_bry.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
nf_fwrite2d_bry.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
nf_fwrite2d_bry.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
nf_fwrite2d_bry.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/pack_field.o

nf_fwrite3d.o: cppdefs.h ducknc.h globaldefs.h
nf_fwrite3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/distribute.o
nf_fwrite3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/get_bounds.o
nf_fwrite3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
nf_fwrite3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
nf_fwrite3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
nf_fwrite3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
nf_fwrite3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
nf_fwrite3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
nf_fwrite3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/pack_field.o
nf_fwrite3d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/stats.o

nf_fwrite3d_bry.o: cppdefs.h ducknc.h globaldefs.h
nf_fwrite3d_bry.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/distribute.o
nf_fwrite3d_bry.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
nf_fwrite3d_bry.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
nf_fwrite3d_bry.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
nf_fwrite3d_bry.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
nf_fwrite3d_bry.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
nf_fwrite3d_bry.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
nf_fwrite3d_bry.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/pack_field.o

nf_fwrite4d.o: cppdefs.h ducknc.h globaldefs.h
nf_fwrite4d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/distribute.o
nf_fwrite4d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/get_bounds.o
nf_fwrite4d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
nf_fwrite4d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
nf_fwrite4d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
nf_fwrite4d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
nf_fwrite4d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
nf_fwrite4d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
nf_fwrite4d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/pack_field.o
nf_fwrite4d.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/stats.o

normalization.o: tile.h cppdefs.h ducknc.h globaldefs.h set_bounds.h
normalization.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/bc_2d.o
normalization.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/bc_3d.o
normalization.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/bc_bry2d.o
normalization.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/bc_bry3d.o
normalization.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/distribute.o
normalization.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_boundary.o
normalization.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
normalization.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_fourdvar.o
normalization.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
normalization.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
normalization.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_kinds.o
normalization.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_mixing.o
normalization.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
normalization.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
normalization.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
normalization.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
normalization.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
normalization.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
normalization.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
normalization.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sedbed.o
normalization.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
normalization.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o
normalization.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fwrite2d.o
normalization.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fwrite3d.o
normalization.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/set_depth.o
normalization.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o
normalization.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/white_noise.o

nrutil.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
nrutil.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_kinds.o

ntimestep.o: cppdefs.h ducknc.h globaldefs.h
ntimestep.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
ntimestep.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_nesting.o
ntimestep.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
ntimestep.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
ntimestep.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o

obs_cost.o: cppdefs.h ducknc.h globaldefs.h
obs_cost.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_fourdvar.o
obs_cost.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
obs_cost.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
obs_cost.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o

obs_depth.o: cppdefs.h ducknc.h globaldefs.h set_bounds.h
obs_depth.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_fourdvar.o
obs_depth.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
obs_depth.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
obs_depth.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o

obs_initial.o: cppdefs.h ducknc.h globaldefs.h
obs_initial.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_fourdvar.o
obs_initial.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
obs_initial.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
obs_initial.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
obs_initial.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
obs_initial.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
obs_initial.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
obs_initial.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
obs_initial.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

obs_k2z.o: cppdefs.h ducknc.h globaldefs.h
obs_k2z.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_kinds.o
obs_k2z.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o

obs_read.o: cppdefs.h ducknc.h globaldefs.h
obs_read.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/dateclock.o
obs_read.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_fourdvar.o
obs_read.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
obs_read.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
obs_read.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
obs_read.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
obs_read.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
obs_read.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
obs_read.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
obs_read.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

obs_write.o: cppdefs.h ducknc.h globaldefs.h set_bounds.h
obs_write.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/distribute.o
obs_write.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/extract_obs.o
obs_write.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_fourdvar.o
obs_write.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
obs_write.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
obs_write.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
obs_write.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
obs_write.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
obs_write.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
obs_write.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
obs_write.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
obs_write.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
obs_write.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
obs_write.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fwrite2d.o
obs_write.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

pack_field.o: cppdefs.h ducknc.h globaldefs.h
pack_field.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/distribute.o
pack_field.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/extract_field.o
pack_field.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
pack_field.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o

packing.o: tile.h set_bounds.h cppdefs.h ducknc.h globaldefs.h
packing.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/distribute.o
packing.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_2d.o
packing.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_3d.o
packing.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
packing.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
packing.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
packing.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
packing.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
packing.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
packing.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
packing.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
packing.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
packing.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
packing.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
packing.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_storage.o
packing.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o
packing.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fread2d.o
packing.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fread3d.o
packing.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

posterior.o: set_bounds.h cppdefs.h ducknc.h globaldefs.h tile.h
posterior.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/distribute.o
posterior.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/lapack_mod.o
posterior.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_boundary.o
posterior.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_coupling.o
posterior.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
posterior.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_fourdvar.o
posterior.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
posterior.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
posterior.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
posterior.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
posterior.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
posterior.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
posterior.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
posterior.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
posterior.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
posterior.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
posterior.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/state_addition.o
posterior.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/state_copy.o
posterior.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/state_dotprod.o
posterior.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/state_initialize.o
posterior.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/state_read.o
posterior.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/state_scale.o
posterior.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o
posterior.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wrt_hessian.o

posterior_var.o: cppdefs.h ducknc.h globaldefs.h set_bounds.h tile.h
posterior_var.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/distribute.o
posterior_var.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_boundary.o
posterior_var.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_coupling.o
posterior_var.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
posterior_var.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_fourdvar.o
posterior_var.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
posterior_var.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
posterior_var.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
posterior_var.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
posterior_var.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
posterior_var.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
posterior_var.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
posterior_var.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
posterior_var.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/state_addition.o
posterior_var.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/state_copy.o
posterior_var.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/state_initialize.o
posterior_var.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/state_product.o
posterior_var.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/state_read.o

ran1.o: cppdefs.h ducknc.h globaldefs.h
ran1.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_kinds.o
ran1.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/ran_state.o

ran_state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_kinds.o
ran_state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nrutil.o

random_ic.o: set_bounds.h cppdefs.h ducknc.h globaldefs.h tile.h
random_ic.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/distribute.o
random_ic.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_boundary.o
random_ic.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_coupling.o
random_ic.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
random_ic.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_fourdvar.o
random_ic.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
random_ic.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
random_ic.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
random_ic.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
random_ic.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
random_ic.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
random_ic.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
random_ic.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
random_ic.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
random_ic.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o
random_ic.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/white_noise.o

read_asspar.o: cppdefs.h ducknc.h globaldefs.h
read_asspar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/inp_decode.o
read_asspar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_fourdvar.o
read_asspar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
read_asspar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
read_asspar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
read_asspar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
read_asspar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
read_asspar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

read_biopar.o: nemuro_inp.h red_tide_inp.h cppdefs.h ducknc.h globaldefs.h
read_biopar.o: hypoxia_srm_inp.h npzd_Franks_inp.h ecosim_inp.h fennel_inp.h
read_biopar.o: npzd_iron_inp.h npzd_Powell_inp.h
read_biopar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/inp_decode.o
read_biopar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_biology.o
read_biopar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_eclight.o
read_biopar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
read_biopar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
read_biopar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
read_biopar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o

read_couplepar.o: cppdefs.h ducknc.h globaldefs.h
read_couplepar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/distribute.o
read_couplepar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/inp_decode.o
read_couplepar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_coupler.o
read_couplepar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
read_couplepar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
read_couplepar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
read_couplepar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o

read_fltbiopar.o: cppdefs.h ducknc.h globaldefs.h oyster_floats_inp.h
read_fltbiopar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/inp_decode.o
read_fltbiopar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_behavior.o
read_fltbiopar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
read_fltbiopar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
read_fltbiopar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
read_fltbiopar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
read_fltbiopar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o

read_fltpar.o: cppdefs.h ducknc.h globaldefs.h
read_fltpar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/inp_decode.o
read_fltpar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_floats.o
read_fltpar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
read_fltpar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
read_fltpar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
read_fltpar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
read_fltpar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
read_fltpar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

read_icepar.o: cppdefs.h ducknc.h globaldefs.h

read_phypar.o: cppdefs.h ducknc.h globaldefs.h
read_phypar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/dateclock.o
read_phypar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/inp_decode.o
read_phypar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_biology.o
read_phypar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_coupler.o
read_phypar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_fourdvar.o
read_phypar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ice.o
read_phypar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
read_phypar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
read_phypar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_nesting.o
read_phypar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
read_phypar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
read_phypar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
read_phypar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
read_phypar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
read_phypar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sediment.o
read_phypar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
read_phypar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_storage.o
read_phypar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_strings.o
read_phypar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_vegetation.o
read_phypar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/set_pio.o
read_phypar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

read_sedpar.o: cppdefs.h ducknc.h globaldefs.h sediment_inp.h
read_sedpar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/inp_decode.o
read_sedpar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
read_sedpar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
read_sedpar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
read_sedpar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
read_sedpar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sedflocs.o
read_sedpar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sediment.o

read_stapar.o: cppdefs.h ducknc.h globaldefs.h
read_stapar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/inp_decode.o
read_stapar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ice.o
read_stapar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
read_stapar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
read_stapar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
read_stapar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
read_stapar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
read_stapar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sediment.o

read_vegpar.o: vegetation_inp.h cppdefs.h ducknc.h globaldefs.h
read_vegpar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/inp_decode.o
read_vegpar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
read_vegpar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
read_vegpar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
read_vegpar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
read_vegpar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sediment.o
read_vegpar.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_vegetation.o

regrid.o: cppdefs.h ducknc.h globaldefs.h
regrid.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/distribute.o
regrid.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/get_varcoords.o
regrid.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/interpolate.o
regrid.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
regrid.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
regrid.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
regrid.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
regrid.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
regrid.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
regrid.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o
regrid.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/shapiro.o
regrid.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

rep_matrix.o: cppdefs.h ducknc.h globaldefs.h
rep_matrix.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/distribute.o
rep_matrix.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_fourdvar.o
rep_matrix.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
rep_matrix.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
rep_matrix.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
rep_matrix.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o

roms_interp.o: cppdefs.h ducknc.h globaldefs.h set_bounds.h
roms_interp.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/distribute.o
roms_interp.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/interpolate.o
roms_interp.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
roms_interp.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
roms_interp.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
roms_interp.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
roms_interp.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
roms_interp.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
roms_interp.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
roms_interp.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
roms_interp.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/stats.o
roms_interp.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

round.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_kinds.o

rpcg_lanczos.o: cppdefs.h ducknc.h globaldefs.h
rpcg_lanczos.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/distribute.o
rpcg_lanczos.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/lapack_mod.o
rpcg_lanczos.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_fourdvar.o
rpcg_lanczos.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
rpcg_lanczos.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
rpcg_lanczos.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
rpcg_lanczos.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
rpcg_lanczos.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
rpcg_lanczos.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
rpcg_lanczos.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
rpcg_lanczos.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

set_2dfld.o: set_bounds.h cppdefs.h ducknc.h globaldefs.h
set_2dfld.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_2d.o
set_2dfld.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
set_2dfld.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
set_2dfld.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
set_2dfld.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
set_2dfld.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
set_2dfld.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o

set_2dfldr.o: set_bounds.h cppdefs.h ducknc.h globaldefs.h
set_2dfldr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_2d.o
set_2dfldr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
set_2dfldr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
set_2dfldr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
set_2dfldr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
set_2dfldr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
set_2dfldr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o

set_3dfld.o: cppdefs.h ducknc.h globaldefs.h set_bounds.h
set_3dfld.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/dateclock.o
set_3dfld.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_3d.o
set_3dfld.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
set_3dfld.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
set_3dfld.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
set_3dfld.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
set_3dfld.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
set_3dfld.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o

set_3dfldr.o: set_bounds.h cppdefs.h ducknc.h globaldefs.h
set_3dfldr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/dateclock.o
set_3dfldr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_3d.o
set_3dfldr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
set_3dfldr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
set_3dfldr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
set_3dfldr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
set_3dfldr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
set_3dfldr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o

set_contact.o: cppdefs.h ducknc.h globaldefs.h
set_contact.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
set_contact.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_nesting.o
set_contact.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
set_contact.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
set_contact.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
set_contact.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
set_contact.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
set_contact.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

set_diags.o: tile.h cppdefs.h ducknc.h globaldefs.h set_bounds.h
set_diags.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/bc_2d.o
set_diags.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/bc_3d.o
set_diags.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/bc_4d.o
set_diags.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_2d.o
set_diags.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_biology.o
set_diags.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_diags.o
set_diags.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
set_diags.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
set_diags.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
set_diags.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
set_diags.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
set_diags.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o

set_grid.o: cppdefs.h ducknc.h globaldefs.h
set_grid.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/analytical.o
set_grid.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/distribute.o
set_grid.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/equilibrium_tide.o
set_grid.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/get_extract.o
set_grid.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/get_grid.o
set_grid.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/get_nudgcoef.o
set_grid.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/metrics.o
set_grid.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_nesting.o
set_grid.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
set_grid.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
set_grid.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
set_grid.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nesting.o
set_grid.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

set_masks.o: tile.h cppdefs.h ducknc.h globaldefs.h set_bounds.h
set_masks.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_2d.o
set_masks.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
set_masks.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
set_masks.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
set_masks.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sources.o
set_masks.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o

set_ngfld.o: cppdefs.h ducknc.h globaldefs.h
set_ngfld.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
set_ngfld.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
set_ngfld.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
set_ngfld.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
set_ngfld.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o

set_ngfldr.o: cppdefs.h ducknc.h globaldefs.h
set_ngfldr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
set_ngfldr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
set_ngfldr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
set_ngfldr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
set_ngfldr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o

set_pio.o: cppdefs.h ducknc.h globaldefs.h
set_pio.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
set_pio.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_kinds.o
set_pio.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
set_pio.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
set_pio.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
set_pio.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
set_pio.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_storage.o
set_pio.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_vegetation.o
set_pio.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

set_scoord.o: cppdefs.h ducknc.h globaldefs.h
set_scoord.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
set_scoord.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
set_scoord.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
set_scoord.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
set_scoord.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o

set_weights.o: cppdefs.h ducknc.h globaldefs.h
set_weights.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
set_weights.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
set_weights.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
set_weights.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o

shapiro.o: set_bounds.h cppdefs.h ducknc.h globaldefs.h
shapiro.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o

sqlq.o: cppdefs.h ducknc.h globaldefs.h
sqlq.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_kinds.o

state_addition.o: cppdefs.h ducknc.h globaldefs.h set_bounds.h
state_addition.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
state_addition.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
state_addition.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o

state_assign.o: cppdefs.h ducknc.h globaldefs.h
state_assign.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
state_assign.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
state_assign.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
state_assign.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o

state_copy.o: set_bounds.h cppdefs.h ducknc.h globaldefs.h
state_copy.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
state_copy.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
state_copy.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o

state_dotprod.o: cppdefs.h ducknc.h globaldefs.h set_bounds.h
state_dotprod.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/distribute.o
state_dotprod.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
state_dotprod.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
state_dotprod.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
state_dotprod.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o

state_initialize.o: cppdefs.h ducknc.h globaldefs.h set_bounds.h
state_initialize.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
state_initialize.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
state_initialize.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o

state_join.o: cppdefs.h ducknc.h globaldefs.h
state_join.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/dateclock.o
state_join.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_coupling.o
state_join.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
state_join.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
state_join.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_mixing.o
state_join.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
state_join.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
state_join.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
state_join.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
state_join.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
state_join.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
state_join.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
state_join.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fread2d.o
state_join.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fread3d.o
state_join.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fwrite2d.o
state_join.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fwrite3d.o
state_join.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

state_product.o: cppdefs.h ducknc.h globaldefs.h set_bounds.h
state_product.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/distribute.o
state_product.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
state_product.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
state_product.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
state_product.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o

state_read.o: cppdefs.h ducknc.h globaldefs.h set_bounds.h
state_read.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/dateclock.o
state_read.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/distribute.o
state_read.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
state_read.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
state_read.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
state_read.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
state_read.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
state_read.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
state_read.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
state_read.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fread2d.o
state_read.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fread2d_bry.o
state_read.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fread3d.o
state_read.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fread3d_bry.o
state_read.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

state_regrid.o: tile.h set_bounds.h cppdefs.h ducknc.h globaldefs.h
state_regrid.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
state_regrid.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
state_regrid.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
state_regrid.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
state_regrid.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
state_regrid.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
state_regrid.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
state_regrid.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
state_regrid.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/roms_interp.o
state_regrid.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/state_assign.o
state_regrid.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

state_scale.o: cppdefs.h ducknc.h globaldefs.h set_bounds.h
state_scale.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
state_scale.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
state_scale.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o

stats.o: cppdefs.h ducknc.h globaldefs.h
stats.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/distribute.o
stats.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/get_hash.o
stats.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
stats.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
stats.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
stats.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o

stats_modobs.o: cppdefs.h ducknc.h globaldefs.h
stats_modobs.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/distribute.o
stats_modobs.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_fourdvar.o
stats_modobs.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
stats_modobs.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
stats_modobs.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
stats_modobs.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
stats_modobs.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
stats_modobs.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
stats_modobs.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
stats_modobs.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
stats_modobs.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/obs_k2z.o
stats_modobs.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

stdinp_mod.o: cppdefs.h ducknc.h globaldefs.h
stdinp_mod.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/distribute.o
stdinp_mod.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/inp_decode.o
stdinp_mod.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
stdinp_mod.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_kinds.o
stdinp_mod.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
stdinp_mod.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

stdout_mod.o: cppdefs.h ducknc.h globaldefs.h
stdout_mod.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
stdout_mod.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
stdout_mod.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
stdout_mod.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

stiffness.o: set_bounds.h cppdefs.h ducknc.h globaldefs.h tile.h
stiffness.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/distribute.o
stiffness.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
stiffness.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
stiffness.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
stiffness.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
stiffness.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
stiffness.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o

strings.o: cppdefs.h ducknc.h globaldefs.h
strings.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/distribute.o
strings.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
strings.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o

sum_grad.o: tile.h set_bounds.h cppdefs.h ducknc.h globaldefs.h
sum_grad.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_boundary.o
sum_grad.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
sum_grad.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
sum_grad.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
sum_grad.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
sum_grad.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o

sum_imp.o: tile.h set_bounds.h cppdefs.h ducknc.h globaldefs.h
sum_imp.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_boundary.o
sum_imp.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
sum_imp.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
sum_imp.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
sum_imp.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o

tadv.o: cppdefs.h ducknc.h globaldefs.h
tadv.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
tadv.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
tadv.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
tadv.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
tadv.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
tadv.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
tadv.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

tides_date.o: cppdefs.h ducknc.h globaldefs.h
tides_date.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/dateclock.o
tides_date.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
tides_date.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
tides_date.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
tides_date.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
tides_date.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
tides_date.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
tides_date.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
tides_date.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

tile_indices.o: cppdefs.h ducknc.h globaldefs.h
tile_indices.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/get_bounds.o
tile_indices.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
tile_indices.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
tile_indices.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
tile_indices.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o

time_corr.o: set_bounds.h cppdefs.h ducknc.h globaldefs.h
time_corr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
time_corr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
time_corr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
time_corr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
time_corr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
time_corr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
time_corr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
time_corr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
time_corr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
time_corr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fread2d.o
time_corr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fread3d.o
time_corr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fwrite2d.o
time_corr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fwrite3d.o
time_corr.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

timers.o: cppdefs.h ducknc.h globaldefs.h
timers.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/distribute.o
timers.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
timers.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
timers.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
timers.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_strings.o
timers.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

uv_rotate.o: cppdefs.h ducknc.h globaldefs.h set_bounds.h
uv_rotate.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_2d.o
uv_rotate.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_3d.o
uv_rotate.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
uv_rotate.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
uv_rotate.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o

uv_var_change.o: tile.h set_bounds.h cppdefs.h ducknc.h globaldefs.h
uv_var_change.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
uv_var_change.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
uv_var_change.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
uv_var_change.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
uv_var_change.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
uv_var_change.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o

vorticity.o: tile.h set_bounds.h cppdefs.h ducknc.h globaldefs.h
vorticity.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_2d.o
vorticity.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_3d.o
vorticity.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_average.o
vorticity.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
vorticity.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
vorticity.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
vorticity.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
vorticity.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
vorticity.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
vorticity.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o

white_noise.o: cppdefs.h ducknc.h globaldefs.h
white_noise.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/distribute.o
white_noise.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_kinds.o
white_noise.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
white_noise.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
white_noise.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
white_noise.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o
white_noise.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nrutil.o

wpoints.o: tile.h cppdefs.h ducknc.h globaldefs.h set_bounds.h
wpoints.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/distribute.o
wpoints.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
wpoints.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
wpoints.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
wpoints.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
wpoints.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
wpoints.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
wpoints.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
wpoints.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_storage.o
wpoints.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/set_pio.o

wrt_aug_imp.o: set_bounds.h cppdefs.h ducknc.h globaldefs.h
wrt_aug_imp.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
wrt_aug_imp.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
wrt_aug_imp.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
wrt_aug_imp.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
wrt_aug_imp.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
wrt_aug_imp.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
wrt_aug_imp.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
wrt_aug_imp.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
wrt_aug_imp.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
wrt_aug_imp.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fwrite2d.o
wrt_aug_imp.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fwrite3d.o
wrt_aug_imp.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

wrt_avg.o: cppdefs.h ducknc.h globaldefs.h
wrt_avg.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/bbl_output.o
wrt_avg.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_average.o
wrt_avg.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
wrt_avg.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
wrt_avg.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
wrt_avg.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
wrt_avg.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
wrt_avg.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
wrt_avg.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
wrt_avg.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
wrt_avg.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
wrt_avg.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sedbed.o
wrt_avg.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sediment.o
wrt_avg.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_tides.o
wrt_avg.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fwrite2d.o
wrt_avg.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fwrite3d.o
wrt_avg.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/sediment_output.o
wrt_avg.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o
wrt_avg.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wec_output.o

wrt_dai.o: set_bounds.h cppdefs.h ducknc.h globaldefs.h
wrt_dai.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
wrt_dai.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
wrt_dai.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_mixing.o
wrt_dai.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
wrt_dai.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
wrt_dai.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
wrt_dai.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
wrt_dai.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
wrt_dai.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
wrt_dai.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
wrt_dai.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
wrt_dai.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fwrite2d.o
wrt_dai.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fwrite3d.o
wrt_dai.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

wrt_diags.o: cppdefs.h ducknc.h globaldefs.h
wrt_diags.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_biology.o
wrt_diags.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_diags.o
wrt_diags.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
wrt_diags.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
wrt_diags.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
wrt_diags.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
wrt_diags.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
wrt_diags.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
wrt_diags.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
wrt_diags.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
wrt_diags.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fwrite2d.o
wrt_diags.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fwrite3d.o
wrt_diags.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fwrite4d.o
wrt_diags.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

wrt_error.o: cppdefs.h ducknc.h globaldefs.h
wrt_error.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_boundary.o
wrt_error.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
wrt_error.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_fourdvar.o
wrt_error.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
wrt_error.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
wrt_error.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_mixing.o
wrt_error.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
wrt_error.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
wrt_error.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
wrt_error.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
wrt_error.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
wrt_error.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
wrt_error.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
wrt_error.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
wrt_error.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fwrite2d.o
wrt_error.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fwrite2d_bry.o
wrt_error.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fwrite3d.o
wrt_error.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fwrite3d_bry.o
wrt_error.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

wrt_evolved.o: cppdefs.h ducknc.h globaldefs.h
wrt_evolved.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_boundary.o
wrt_evolved.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
wrt_evolved.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
wrt_evolved.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
wrt_evolved.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_mixing.o
wrt_evolved.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
wrt_evolved.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
wrt_evolved.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
wrt_evolved.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
wrt_evolved.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
wrt_evolved.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
wrt_evolved.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
wrt_evolved.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
wrt_evolved.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fwrite2d.o
wrt_evolved.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fwrite2d_bry.o
wrt_evolved.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fwrite3d.o
wrt_evolved.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fwrite3d_bry.o
wrt_evolved.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

wrt_extract.o: set_bounds.h cppdefs.h ducknc.h globaldefs.h
wrt_extract.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/bbl_output.o
wrt_extract.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/extract_field.o
wrt_extract.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_bbl.o
wrt_extract.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_boundary.o
wrt_extract.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_coupling.o
wrt_extract.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_extract.o
wrt_extract.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
wrt_extract.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
wrt_extract.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
wrt_extract.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_mixing.o
wrt_extract.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
wrt_extract.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
wrt_extract.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
wrt_extract.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
wrt_extract.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
wrt_extract.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
wrt_extract.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
wrt_extract.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sedbed.o
wrt_extract.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sediment.o
wrt_extract.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
wrt_extract.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fwrite2d.o
wrt_extract.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fwrite2d_bry.o
wrt_extract.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fwrite3d.o
wrt_extract.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fwrite3d_bry.o
wrt_extract.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/omega.o
wrt_extract.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/sediment_output.o
wrt_extract.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o
wrt_extract.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/uv_rotate.o
wrt_extract.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wec_output.o

wrt_floats.o: cppdefs.h ducknc.h globaldefs.h
wrt_floats.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_floats.o
wrt_floats.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
wrt_floats.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
wrt_floats.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
wrt_floats.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
wrt_floats.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
wrt_floats.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
wrt_floats.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
wrt_floats.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
wrt_floats.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

wrt_gst.o: cppdefs.h ducknc.h globaldefs.h
wrt_gst.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/distribute.o
wrt_gst.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
wrt_gst.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
wrt_gst.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
wrt_gst.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
wrt_gst.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
wrt_gst.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
wrt_gst.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
wrt_gst.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_storage.o
wrt_gst.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

wrt_hessian.o: cppdefs.h ducknc.h globaldefs.h
wrt_hessian.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_boundary.o
wrt_hessian.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
wrt_hessian.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
wrt_hessian.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
wrt_hessian.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_mixing.o
wrt_hessian.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
wrt_hessian.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
wrt_hessian.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
wrt_hessian.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
wrt_hessian.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
wrt_hessian.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
wrt_hessian.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
wrt_hessian.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
wrt_hessian.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fwrite2d.o
wrt_hessian.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fwrite2d_bry.o
wrt_hessian.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fwrite3d.o
wrt_hessian.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fwrite3d_bry.o
wrt_hessian.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

wrt_his.o: cppdefs.h ducknc.h globaldefs.h set_bounds.h
wrt_his.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/bbl_output.o
wrt_his.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/extract_slice.o
wrt_his.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_bbl.o
wrt_his.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_boundary.o
wrt_his.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_coupling.o
wrt_his.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
wrt_his.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
wrt_his.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
wrt_his.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_mixing.o
wrt_his.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
wrt_his.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
wrt_his.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
wrt_his.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
wrt_his.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
wrt_his.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
wrt_his.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
wrt_his.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sedbed.o
wrt_his.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sediment.o
wrt_his.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
wrt_his.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fwrite2d.o
wrt_his.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fwrite2d_bry.o
wrt_his.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fwrite3d.o
wrt_his.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fwrite3d_bry.o
wrt_his.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/omega.o
wrt_his.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/sediment_output.o
wrt_his.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o
wrt_his.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/uv_rotate.o
wrt_his.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/vegetation_output.o
wrt_his.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/vorticity.o
wrt_his.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wec_output.o

wrt_impulse.o: cppdefs.h ducknc.h globaldefs.h set_bounds.h
wrt_impulse.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
wrt_impulse.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
wrt_impulse.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
wrt_impulse.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
wrt_impulse.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
wrt_impulse.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
wrt_impulse.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
wrt_impulse.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
wrt_impulse.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
wrt_impulse.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fread2d.o
wrt_impulse.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fread3d.o
wrt_impulse.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fwrite2d.o
wrt_impulse.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fwrite3d.o
wrt_impulse.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

wrt_info.o: hypoxia_srm_wrt_pio.h cppdefs.h ducknc.h globaldefs.h
wrt_info.o: oyster_floats_wrt_pio.h oyster_floats_wrt.h nemuro_wrt_pio.h
wrt_info.o: npzd_Powell_wrt.h red_tide_wrt_pio.h npzd_Franks_wrt_pio.h
wrt_info.o: npzd_iron_wrt.h red_tide_wrt.h ecosim_wrt_pio.h fennel_wrt.h
wrt_info.o: npzd_Franks_wrt.h nemuro_wrt.h sediment_wrt_pio.h sediment_wrt.h
wrt_info.o: ecosim_wrt.h fennel_wrt_pio.h npzd_iron_wrt_pio.h hypoxia_srm_wrt.h
wrt_info.o: npzd_Powell_wrt_pio.h
wrt_info.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/distribute.o
wrt_info.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/extract_sta.o
wrt_info.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_behavior.o
wrt_info.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_biology.o
wrt_info.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_eclight.o
wrt_info.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_extract.o
wrt_info.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_fourdvar.o
wrt_info.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
wrt_info.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
wrt_info.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
wrt_info.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
wrt_info.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
wrt_info.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
wrt_info.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
wrt_info.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
wrt_info.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sediment.o
wrt_info.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sources.o
wrt_info.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_storage.o
wrt_info.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fwrite2d.o
wrt_info.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

wrt_ini.o: cppdefs.h ducknc.h globaldefs.h
wrt_ini.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/dateclock.o
wrt_ini.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_boundary.o
wrt_ini.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
wrt_ini.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_fourdvar.o
wrt_ini.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
wrt_ini.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
wrt_ini.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_mixing.o
wrt_ini.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
wrt_ini.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
wrt_ini.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
wrt_ini.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
wrt_ini.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
wrt_ini.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
wrt_ini.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
wrt_ini.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sediment.o
wrt_ini.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
wrt_ini.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fwrite2d.o
wrt_ini.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fwrite2d_bry.o
wrt_ini.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fwrite3d.o
wrt_ini.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fwrite3d_bry.o
wrt_ini.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

wrt_quick.o: set_bounds.h cppdefs.h ducknc.h globaldefs.h
wrt_quick.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/bbl_output.o
wrt_quick.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/extract_slice.o
wrt_quick.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_bbl.o
wrt_quick.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_coupling.o
wrt_quick.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
wrt_quick.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
wrt_quick.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
wrt_quick.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_mixing.o
wrt_quick.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
wrt_quick.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
wrt_quick.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
wrt_quick.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
wrt_quick.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
wrt_quick.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
wrt_quick.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
wrt_quick.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sedbed.o
wrt_quick.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sediment.o
wrt_quick.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
wrt_quick.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fwrite2d.o
wrt_quick.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fwrite3d.o
wrt_quick.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/omega.o
wrt_quick.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/sediment_output.o
wrt_quick.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o
wrt_quick.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/uv_rotate.o
wrt_quick.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/vorticity.o
wrt_quick.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wec_output.o

wrt_rst.o: cppdefs.h ducknc.h globaldefs.h
wrt_rst.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
wrt_rst.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
wrt_rst.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_mixing.o
wrt_rst.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
wrt_rst.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
wrt_rst.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
wrt_rst.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
wrt_rst.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
wrt_rst.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
wrt_rst.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
wrt_rst.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sedbed.o
wrt_rst.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sediment.o
wrt_rst.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
wrt_rst.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fwrite2d.o
wrt_rst.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fwrite3d.o
wrt_rst.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fwrite4d.o
wrt_rst.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

wrt_state.o: cppdefs.h ducknc.h globaldefs.h
wrt_state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/dateclock.o
wrt_state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_boundary.o
wrt_state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_coupling.o
wrt_state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
wrt_state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
wrt_state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
wrt_state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_mixing.o
wrt_state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
wrt_state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
wrt_state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
wrt_state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
wrt_state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
wrt_state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
wrt_state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
wrt_state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
wrt_state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fwrite2d.o
wrt_state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fwrite2d_bry.o
wrt_state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fwrite3d.o
wrt_state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fwrite3d_bry.o
wrt_state.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

wrt_station.o: cppdefs.h ducknc.h globaldefs.h
wrt_station.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/bbl_output.o
wrt_station.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/extract_sta.o
wrt_station.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_bbl.o
wrt_station.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
wrt_station.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
wrt_station.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
wrt_station.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_mixing.o
wrt_station.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
wrt_station.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
wrt_station.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
wrt_station.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
wrt_station.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
wrt_station.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
wrt_station.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
wrt_station.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sedbed.o
wrt_station.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_sediment.o
wrt_station.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
wrt_station.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/sediment_output.o
wrt_station.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o
wrt_station.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/uv_rotate.o
wrt_station.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wec_output.o

wrt_std.o: cppdefs.h ducknc.h globaldefs.h
wrt_std.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_forces.o
wrt_std.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_fourdvar.o
wrt_std.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
wrt_std.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
wrt_std.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_mixing.o
wrt_std.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
wrt_std.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
wrt_std.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
wrt_std.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
wrt_std.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
wrt_std.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
wrt_std.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
wrt_std.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
wrt_std.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fwrite2d.o
wrt_std.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fwrite3d.o
wrt_std.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

wrt_tides.o: cppdefs.h ducknc.h globaldefs.h
wrt_tides.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
wrt_tides.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
wrt_tides.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
wrt_tides.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_netcdf.o
wrt_tides.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
wrt_tides.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
wrt_tides.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_pio_netcdf.o
wrt_tides.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
wrt_tides.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
wrt_tides.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_tides.o
wrt_tides.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fwrite3d.o
wrt_tides.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fwrite4d.o
wrt_tides.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings.o

yaml_parser.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
yaml_parser.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_kinds.o
yaml_parser.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
yaml_parser.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o

zeta_balance.o: tile.h set_bounds.h cppdefs.h ducknc.h globaldefs.h
zeta_balance.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/distribute.o
zeta_balance.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_2d.o
zeta_balance.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_3d.o
zeta_balance.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_coupling.o
zeta_balance.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_fourdvar.o
zeta_balance.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_grid.o
zeta_balance.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_iounits.o
zeta_balance.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_mixing.o
zeta_balance.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ncparam.o
zeta_balance.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_ocean.o
zeta_balance.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_parallel.o
zeta_balance.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_param.o
zeta_balance.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_scalars.o
zeta_balance.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mod_stepping.o
zeta_balance.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange.o
zeta_balance.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/rho_eos.o
zeta_balance.o: /storage/arango/ROMS/Projects/DUCKNC/Build_romsM/set_depth.o

/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/adfromtl_mod.mod: ADfromTL.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/analytical_mod.mod: analytical.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/array_modes_mod.mod: array_modes.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/back_cost_mod.mod: back_cost.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/background_std_mod.mod: background_std.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/bbl_mod.mod: bbl.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/bbl_output_mod.mod: bbl_output.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/bc_2d_mod.mod: bc_2d.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/bc_3d_mod.mod: bc_3d.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/bc_4d_mod.mod: bc_4d.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/bc_bry2d_mod.mod: bc_bry2d.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/bc_bry3d_mod.mod: bc_bry3d.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/biology_floats_mod.mod: biology_floats.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/biology_mod.mod: biology.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/bulk_flux_mod.mod: bulk_flux.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/bvf_mix_mod.mod: bvf_mix.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/cgradient_mod.mod: cgradient.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/checkvars_mod.mod: checkvars.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/close_io_mod.mod: close_io.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/cmeps_roms_mod.mod: esmf_roms.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/comp_jb0_mod.mod: comp_Jb0.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/congrad_mod.mod: congrad.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/conv_2d_mod.mod: conv_2d.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/conv_3d_bry_mod.mod: conv_bry3d.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/conv_3d_mod.mod: conv_3d.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/conv_bry2d_mod.mod: conv_bry2d.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/convolve_mod.mod: convolve.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/cost_grad_mod.mod: cost_grad.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/coupler_mod.mod: coupler.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/dateclock_mod.mod: dateclock.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_avg_mod.mod: def_avg.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_dai_mod.mod: def_dai.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_diags_mod.mod: def_diags.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_dim_mod.mod: def_dim.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_error_mod.mod: def_error.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_extract_mod.mod: def_extract.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_floats_mod.mod: def_floats.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_gst_mod.mod: def_gst.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_hessian_mod.mod: def_hessian.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_his_mod.mod: def_his.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_impulse_mod.mod: def_impulse.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_info_mod.mod: def_info.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_ini_mod.mod: def_ini.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_lanczos_mod.mod: def_lanczos.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_mod_mod.mod: def_mod.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_norm_mod.mod: def_norm.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_quick_mod.mod: def_quick.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_rst_mod.mod: def_rst.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_state_mod.mod: def_state.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_station_mod.mod: def_station.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_std_mod.mod: def_std.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_tides_mod.mod: def_tides.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/def_var_mod.mod: def_var.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/destroy_mod.mod: destroy.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/diag_mod.mod: diag.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/distribute_mod.mod: distribute.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/dotproduct_mod.mod: dotproduct.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/equilibrium_tide_mod.mod: equilibrium_tide.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/erf_mod.mod: erf.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/esmf_atm_mod.mod: esmf_atm.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/esmf_coamps_mod.mod: esmf_atm.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/esmf_coupler_mod.mod: coupler.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/esmf_data_mod.mod: esmf_data.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/esmf_esm_mod.mod: esmf_esm.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/esmf_ice_mod.mod: esmf_ice.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/esmf_regcm_mod.mod: esmf_atm.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/esmf_roms_mod.mod: esmf_roms.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/esmf_wam_mod.mod: esmf_wav.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/esmf_wav_mod.mod: esmf_wav.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/esmf_wrf_mod.mod: esmf_atm.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_2d_mod.mod: exchange_2d.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_2d_xtr_mod.mod: exchange_2d_xtr.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_3d_mod.mod: exchange_3d.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_3d_xtr_mod.mod: exchange_3d_xtr.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/exchange_4d_mod.mod: exchange_4d.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/extract_field_mod.mod: extract_field.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/extract_obs_mod.mod: extract_obs.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/extract_slice_mod.mod: extract_slice.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/extract_sta_mod.mod: extract_sta.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/forcing_mod.mod: forcing.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/frc_adjust_mod.mod: frc_adjust.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/frc_iau_mod.mod: frc_iau.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/frc_weak_mod.mod: frc_weak.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/get_bounds_mod.mod: get_bounds.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/get_cycle_mod.mod: get_cycle.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/get_env_mod.mod: get_env.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/get_extract_mod.mod: get_extract.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/get_grid_mod.mod: get_grid.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/get_gst_mod.mod: get_gst.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/get_hash_mod.mod: get_hash.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/get_metadata_mod.mod: get_metadata.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/get_nudgcoef_mod.mod: get_nudgcoef.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/get_state_mod.mod: get_state.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/get_varcoords_mod.mod: get_varcoords.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/get_wetdry_mod.mod: get_wetdry.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/gls_corstep_mod.mod: gls_corstep.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/gls_prestep_mod.mod: gls_prestep.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/hmixing_mod.mod: hmixing.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/i4dvar_mod.mod: i4dvar.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/ini_adjust_mod.mod: ini_adjust.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/ini_fields_mod.mod: ini_fields.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/ini_hmixcoef_mod.mod: ini_hmixcoef.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/ini_lanczos_mod.mod: ini_lanczos.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/inner2state_mod.mod: inner2state.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/inp_decode_mod.mod: inp_decode.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/inp_par_mod.mod: inp_par.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/inquiry_mod.mod: inquiry.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/interp_floats_mod.mod: interp_floats.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/lanc_resid_mod.mod: lanc_resid.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/lbc_mod.mod: lbc.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/lmd_bkpp_mod.mod: lmd_bkpp.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/lmd_skpp_mod.mod: lmd_skpp.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/lmd_vmix_mod.mod: lmd_vmix.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/marsh_dynamics_mod.mod: marsh_dynamics.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/marsh_sed_erosion_mod.mod: marsh_sed_erosion.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/marsh_tidal_range_mod.mod: marsh_tidal_range.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/marsh_vert_growth_mod.mod: marsh_vert_growth.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/marsh_wave_thrust_mod.mod: marsh_wave_thrust.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mct_coupler_mod.mod: coupler.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/metrics_mod.mod: metrics.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mp_exchange_mod.mod: mp_exchange.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/mpdata_adiff_mod.mod: mpdata_adiff.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/my25_corstep_mod.mod: my25_corstep.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/my25_prestep_mod.mod: my25_prestep.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nesting_mod.mod: nesting.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fread2d_bry_mod.mod: nf_fread2d_bry.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fread2d_mod.mod: nf_fread2d.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fread2d_xtr_mod.mod: nf_fread2d_xtr.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fread3d_bry_mod.mod: nf_fread3d_bry.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fread3d_mod.mod: nf_fread3d.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fread3d_xtr_mod.mod: nf_fread3d_xtr.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fread4d_mod.mod: nf_fread4d.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fwrite2d_bry_mod.mod: nf_fwrite2d_bry.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fwrite2d_mod.mod: nf_fwrite2d.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fwrite3d_bry_mod.mod: nf_fwrite3d_bry.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fwrite3d_mod.mod: nf_fwrite3d.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/nf_fwrite4d_mod.mod: nf_fwrite4d.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/normalization_mod.mod: normalization.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/obc_adjust_mod.mod: obc_adjust.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/obc_volcons_mod.mod: obc_volcons.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/obs_initial_mod.mod: obs_initial.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/obs_k2z_mod.mod: obs_k2z.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/obs_read_mod.mod: obs_read.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/obs_write_mod.mod: obs_write.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/omega_mod.mod: omega.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/pack_field_mod.mod: pack_field.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/packing_mod.mod: packing.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/post_initial_mod.mod: post_initial.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/posterior_mod.mod: posterior.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/posterior_var_mod.mod: posterior_var.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/pre_step3d_mod.mod: pre_step3d.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/propagator_mod.mod: propagator.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/prsgrd_mod.mod: prsgrd.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/r4dvar_mod.mod: r4dvar.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/random_ic_mod.mod: random_ic.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/rbl4dvar_mod.mod: rbl4dvar.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/regrid_mod.mod: regrid.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/rho_eos_mod.mod: rho_eos.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/rhs3d_mod.mod: rhs3d.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/roms_interp_mod.mod: roms_interp.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/roms_interpolate_mod.mod: interpolate.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/roms_kernel_mod.mod: roms_kernel.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/round_mod.mod: round.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/rpcg_lanczos_mod.mod: rpcg_lanczos.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/sed_bed_cohesive_mod.mod: sed_bed_cohesive.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/sed_bed_mod.mod: sed_bed.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/sed_bed_mod2.mod: sed_bed2.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/sed_bedload_mod.mod: sed_bedload.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/sed_bedload_vandera_mod.mod: sed_bedload_vandera.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/sed_biodiff_mod.mod: sed_biodiff.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/sed_biomass_mod.mod: vegetation_biomass.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/sed_flocs_mod.mod: sed_flocs.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/sed_fluxes_mod.mod: sed_fluxes.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/sed_settling_mod.mod: sed_settling.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/sed_surface_mod.mod: sed_surface.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/sediment_mod.mod: sediment.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/sediment_output_mod.mod: sediment_output.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/sedtr_decay_mod.mod: sedtr_decay.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/sedtr_reactions_pom_mod.mod: sedtr_reactions_pom.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/sedtr_reactions_sed_decay_mod.mod: sedtr_reactions_sed_decay.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/set_2dfld_mod.mod: set_2dfld.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/set_2dfldr_mod.mod: set_2dfldr.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/set_3dfld_mod.mod: set_3dfld.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/set_3dfldr_mod.mod: set_3dfldr.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/set_avg_mod.mod: set_avg.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/set_contact_mod.mod: set_contact.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/set_depth_mod.mod: set_depth.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/set_masks_mod.mod: set_masks.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/set_massflux_mod.mod: set_massflux.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/set_pio_mod.mod: set_pio.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/set_tides_mod.mod: set_tides.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/set_vbc_mod.mod: set_vbc.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/set_zeta_mod.mod: set_zeta.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/shapiro_mod.mod: shapiro.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/state_addition_mod.mod: state_addition.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/state_assign_mod.mod: state_assign.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/state_copy_mod.mod: state_copy.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/state_dotprod_mod.mod: state_dotprod.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/state_initialize_mod.mod: state_initialize.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/state_join_mod.mod: state_join.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/state_product_mod.mod: state_product.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/state_read_mod.mod: state_read.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/state_regrid_mod.mod: state_regrid.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/state_scale_mod.mod: state_scale.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/stats_mod.mod: stats.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/stats_modobs_mod.mod: stats_modobs.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/step2d_mod.mod: step2d.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/step3d_t_mod.mod: step3d_t.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/step3d_uv_mod.mod: step3d_uv.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/step_floats_mod.mod: step_floats.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/stiffness_mod.mod: stiffness.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/strings_mod.mod: strings.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/sum_grad_mod.mod: sum_grad.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/sum_imp_mod.mod: sum_imp.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/t3dbc_mod.mod: t3dbc_im.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/t3dmix2_mod.mod: t3dmix.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/t3dmix4_mod.mod: t3dmix.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/tadv_mod.mod: tadv.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/tile_indices_mod.mod: tile_indices.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/time_corr_mod.mod: time_corr.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/tkebc_mod.mod: tkebc_im.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/u2dbc_mod.mod: u2dbc_im.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/u3dbc_mod.mod: u3dbc_im.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/uv3dmix2_mod.mod: uv3dmix.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/uv3dmix4_mod.mod: uv3dmix.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/uv_rotate_mod.mod: uv_rotate.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/uv_var_change_mod.mod: uv_var_change.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/v2dbc_mod.mod: v2dbc_im.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/v3dbc_mod.mod: v3dbc_im.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/vegetation_drag_mod.mod: vegetation_drag.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/vegetation_output_mod.mod: vegetation_output.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/vegetation_stream_mod.mod: vegetation_stream.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/vegetation_turb_mod.mod: vegetation_turb.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/vorticity_mod.mod: vorticity.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/vwalk_floats_mod.mod: vwalk_floats.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wec_dissip_mod.mod: wec_dissip.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wec_output_mod.mod: wec_output.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wec_roller_mod.mod: wec_roller.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wec_stokes_mod.mod: wec_stokes.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wec_streaming_mod.mod: wec_streaming.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wec_u2dbc_mod.mod: wec_u2dbc_im.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wec_u3dbc_mod.mod: wec_u3dbc_im.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wec_v2dbc_mod.mod: wec_v2dbc_im.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wec_v3dbc_mod.mod: wec_v3dbc_im.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wec_vf_mod.mod: wec_vf.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wec_wave_mix_mod.mod: wec_wave_mix.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wec_wvelocity_mod.mod: wec_wvelocity.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wetdry_mod.mod: wetdry.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/white_noise_mod.mod: white_noise.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wpoints_mod.mod: wpoints.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wrt_aug_imp_mod.mod: wrt_aug_imp.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wrt_avg_mod.mod: wrt_avg.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wrt_dai_mod.mod: wrt_dai.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wrt_diags_mod.mod: wrt_diags.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wrt_error_mod.mod: wrt_error.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wrt_evolved_mod.mod: wrt_evolved.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wrt_extract_mod.mod: wrt_extract.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wrt_floats_mod.mod: wrt_floats.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wrt_gst_mod.mod: wrt_gst.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wrt_hessian_mod.mod: wrt_hessian.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wrt_his_mod.mod: wrt_his.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wrt_impulse_mod.mod: wrt_impulse.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wrt_info_mod.mod: wrt_info.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wrt_ini_mod.mod: wrt_ini.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wrt_quick_mod.mod: wrt_quick.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wrt_rst_mod.mod: wrt_rst.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wrt_state_mod.mod: wrt_state.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wrt_station_mod.mod: wrt_station.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wrt_std_mod.mod: wrt_std.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wrt_tides_mod.mod: wrt_tides.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/wvelocity_mod.mod: wvelocity.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/yaml_parser_mod.mod: yaml_parser.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/zeta_balance_mod.mod: zeta_balance.o
/storage/arango/ROMS/Projects/DUCKNC/Build_romsM/zetabc_mod.mod: zetabc.o
