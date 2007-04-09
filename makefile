# $Id$
#::::::::::::::::::::::::::::::::::::::::::::::::::::: Hernan G. Arango :::
# Copyright (c) 2002-2007 The ROMS/TOMS Group             Kate Hedstrom :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#                                                                       :::
#  ROMS/TOMS Framework Master Makefile                                  :::
#                                                                       :::
#  This makefile is designed to work only with GNU Make version 3.77 or :::
#  higher. It can be used in any architecture provided that there is a  :::
#  machine/compiler rules file in the  "Compilers"  subdirectory.  You  :::
#  may need to modify the rules file to specify the  correct path  for  :::
#  the NetCDF and ARPACK libraries. The ARPACK library is only used in  :::
#  the Generalized Stability Theory analysis.                           :::
#                                                                       :::
#  If appropriate,  the USER needs to modify the  macro definitions in  :::
#  in user-defined section below.  To activate an option set the macro  :::
#  to "on". For example, if you want to compile with debugging options  :::
#  set:                                                                 :::
#                                                                       :::
#      DEBUG := on                                                      :::
#                                                                       :::
#  Otherwise, leave macro definition blank.                             :::
#                                                                       :::
#  The USER needs to provide a value for the  macro FORT.  Choose  the  :::
#  appropriate value from the list below.                               :::
#                                                                       :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


#--------------------------------------------------------------------------
#  Initialize some things.
#--------------------------------------------------------------------------

  clean_list := core *.o *.oo *.mod *.f90 lib*.a *.bak
  sources    := 
  path_srcs  := 
  libraries  :=

     objects  = $(subst .F,.o,$(sources))

#==========================================================================
#  Start of user-defined options. Modify macro variables: on is TRUE while
#  blank is FALSE.
#==========================================================================
#
#  The CPP option defining a particular application is specified below.
#  See header file "ROMS/Include/cppdefs.h" for all available idealized
#  and realistic applications CPP flags. For example, to activate the
#  upwelling test case (UPWELLING) set:
#
#    ROMS_APPLICATION := UPWELLING
#
#  Notice that this makefile will include the associated application header
#  file, which is located either in the "ROMS/Include" or MY_HEADER_DIR
#  directory.  This makefile is designed to search in both directories.
#  The only constrain is that the application CPP option must be unique
#  and header file name is the lowercase value of ROMS_APPLICATION with
#  the .h extension. For example, the upwelling application includes the
#  "upwelling.h" header file.  

ROMS_APPLICATION := UPWELLING

#  If application header files is not located in the "ROMS/Include",
#  provide alternate directory.

MY_HEADER_DIR :=

#  Activate option below to include user analytical expressions from
#  the provided templates in "User/Functionals" or MY_ANALYTICAL_DIR
#  directory instead of official distributed expressions in
#  "ROMS/Functionals". Notice that the user have the choice of 
#  specifying any other directory for these expressions in macro
#  MY_ANALYTICAL_DIR.

MY_ANALYTICAL := 

#  Provide directory path for your own analytical expressions. Notice
#  that a template is provided in "User/Functionals".

ifdef MY_ANALYTICAL
  MY_ANALYTICAL_DIR :=
endif

#  Set number of ROMS nested and/or composed grid.  Currently, only
#  one grid is supported.  This option will be available in the near
#  future.

 NestedGrids := 1

#  Activate debugging compiler options:

       DEBUG := on

#  If parallel applications, use at most one of these definitions
#  (leave both definitions blank in serial applications):

         MPI := 
      OpenMP :=

#  If distributed-memory, turn on compilation via the script "mpif90".
#  This is needed in some Linux operating systems. In some systems with
#  native MPI libraries the compilation does not require MPICH type
#  scripts. This macro is also convient when there are several fortran
#  compiliers (ifort, pgf90, pathf90) in the system that use mpif90.
#  In this, case the user need to select the desired compiler below and
#  turn on both MPI and MPIF90 macros.

      MPIF90 :=

#  If applicable, compile with the ARPACK library (GST analysis):

      ARPACK :=

#  If applicable, activate 64-bit compilation:

       LARGE :=

#  If applicable, activate Coupling to SWAN wave model.

      SWAN_COUPLE := 

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
#     Darwin:                 f90, xlf
#     Linux:                  ifc, ifort, pgi, path, g95
#     SunOS:                  f95
#     UNICOS-mp:              ftn
#     SunOS/Linux:            ftn (Cray cross-compiler)
#
#  Feel free to send us additional rule files to include! Also, be sure
#  to check the appropriate file to make sure it has the right paths to
#  NetCDF and so on.
#--------------------------------------------------------------------------

        FORT ?= pgi

#--------------------------------------------------------------------------
#  Set directory for executable.
#--------------------------------------------------------------------------

      BINDIR := .

#==========================================================================
#  End of user-defined options. See also the machine-dependent include
#  file being used above.
#==========================================================================

#--------------------------------------------------------------------------
#  Set Pattern rules.
#--------------------------------------------------------------------------

%.o: %.F

%.o: %.f90
	$(FC) -c $(FFLAGS) $<

%.f90: %.F
	$(CPP) $(CPPFLAGS) $< > $*.f90
	$(CLEAN) $*.f90

CLEAN := ROMS/Bin/cpp_clean

#--------------------------------------------------------------------------
#  Set ROMS/TOMS executable file name.
#--------------------------------------------------------------------------

BIN := $(BINDIR)/oceanS
ifdef DEBUG
  BIN := $(BINDIR)/oceanG
else
 ifdef MPI
   BIN := $(BINDIR)/oceanM
 endif
 ifdef OpenMP
   BIN := $(BINDIR)/oceanO
 endif
endif

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

SVNREV := $(shell svnversion -n .)

ROOTDIR := $(shell pwd)

COMPILERS := ./Compilers

ifndef FORT
  $(error Variable FORT not set)
endif

ifneq "$(MAKECMDGOALS)" "clean"
  include $(COMPILERS)/$(OS)-$(strip $(FORT)).mk
endif

#--------------------------------------------------------------------------
#  Pass the platform variables to the preprocessor as macros. Convert to
#  valid, upper-case identifiers. If applicable, attach ROMS application
#  CPP option.
#--------------------------------------------------------------------------

CPPFLAGS += -D$(shell echo ${OS} | tr "-" "_" | tr [a-z] [A-Z])
CPPFLAGS += -D$(shell echo ${CPU} | tr "-" "_" | tr [a-z] [A-Z])
CPPFLAGS += -D$(shell echo ${FORT} | tr "-" "_" | tr [a-z] [A-Z])


CPPFLAGS += -D'ROOT_DIR="$(ROOTDIR)"'
ifdef ROMS_APPLICATION
  HEADER := $(addsuffix .h,$(shell echo ${ROMS_APPLICATION} | tr [A-Z] [a-z]))
  CPPFLAGS += -D$(ROMS_APPLICATION)
  CPPFLAGS += -D'ROMS_HEADER="$(HEADER)"'
  CPPFLAGS += -DNestedGrids=$(NestedGrids)
endif

ifndef MY_ANALYTICAL
  MY_ANALYTICAL_DIR := ./ROMS/Functionals
endif
CPPFLAGS += -D'ANALYTICAL_DIR="$(MY_ANALYTICAL_DIR)"'

ifdef MY_ANALYTICAL
  CPPFLAGS += -D'MY_ANALYTICAL="$(MY_ANALYTICAL)"'
endif

ifdef SVNREV
  CPPFLAGS += -D'SVN_REV="$(SVNREV)"'
else
  SVNREV := $(shell grep Revision ./ROMS/Version | sed 's/.* \([0-9]*\) .*/\1/')
  CPPFLAGS += -D'SVN_REV="$(SVNREV)"'
endif  

#--------------------------------------------------------------------------
#  Build target directories.
#--------------------------------------------------------------------------

.PHONY: all

all: $(BIN)

modules   :=	ROMS/Adjoint \
		ROMS/Representer \
		ROMS/Tangent \
		ROMS/Nonlinear \
		ROMS/Functionals \
		ROMS/SeaIce \
		ROMS/Utility \
		ROMS/Modules

includes  :=	ROMS/Include \
		ROMS/Functionals \
		ROMS/Adjoint \
		ROMS/Nonlinear \
		ROMS/Representer \
		ROMS/Tangent \
		ROMS/SeaIce \
		ROMS/Utility \
		ROMS/Drivers

ifdef MY_ANALYTICAL
 includes +=	$(MY_ANALYTICAL_DIR)
endif

ifdef MY_HEADER_DIR
 includes +=	$(MY_HEADER_DIR)
endif

ifdef SWAN_COUPLE
 modules  +=	SWAN/Src
 includes +=	SWAN/Src
endif

modules   +=	Master
includes  +=	Master

vpath %.F $(modules)
vpath %.h $(includes)

include $(addsuffix /Module.mk,$(modules))

MDEPFLAGS += $(patsubst %,-I %,$(includes))

CPPFLAGS += $(patsubst %,-I%,$(includes))

ifdef MY_HEADER_DIR
  CPPFLAGS += -D'HEADER_DIR="$(MY_HEADER_DIR)"'
else
  CPPFLAGS += -D'HEADER_DIR="./ROMS/Include"'
endif

#--------------------------------------------------------------------------
#  Add profiling.
#--------------------------------------------------------------------------

# FFLAGS += -check bounds                 # ifort
# FFLAGS += -C                            # pgi
# FFLAGS += -xpg                          # Sun
# FFLAGS += -pg                           # g95
# FFLAGS += -qp                           # ifort
# FFLAGS += -Mprof=func,lines             # pgi
# FFLAGS += -Mprof=mpi,lines              # pgi
# FFLAGS += -Mprof=mpi,hwcts              # pgi
# FFLAGS += -Mprof=func                   # pgi

#--------------------------------------------------------------------------
#  Special CPP macros for mod_strings.F
#--------------------------------------------------------------------------

mod_strings.f90: CPPFLAGS += -DMY_OS="'$(OS)'" -DMY_CPU="'$(CPU)'" \
                             -DMY_FORT="'$(FORT)'" -DMY_FC="'$(FC)'" \
                             -DMY_FFLAGS="'$(FFLAGS)'"

#--------------------------------------------------------------------------
#  ROMS/TOMS libraries.
#--------------------------------------------------------------------------

MYLIB := libocean.a

.PHONY: libraries

libraries: $(libraries)

#--------------------------------------------------------------------------
#  Target to create ROMS/TOMS dependecies.
#--------------------------------------------------------------------------

.PHONY: depend

SFMAKEDEPEND := ./ROMS/Bin/sfmakedepend

depend:
	mv $(COMPILERS)/MakeDepend $(COMPILERS)/MakeDepend.orig
	$(SFMAKEDEPEND) $(MDEPFLAGS) $(path_srcs) > $(COMPILERS)/MakeDepend 

ifneq "$(MAKECMDGOALS)" "clean"
  -include $(COMPILERS)/MakeDepend
endif

#--------------------------------------------------------------------------
#  Target to create ROMS/TOMS tar file.
#--------------------------------------------------------------------------

.PHONY: tarfile

tarfile:
		tar -cvf roms-3_0.tar *

.PHONY: zipfile

zipfile:
		zip -r roms-3_0.zip *

.PHONY: gzipfile

gzipfile:
		gzip -v roms-3_0.gzip *

#--------------------------------------------------------------------------
#  Cleaning targets.
#--------------------------------------------------------------------------

.PHONY: clean

clean:
	$(RM) $(clean_list)

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

# DO NOT DELETE
