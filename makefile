#:::::::::::::::::::::::::::::::::::::::::::::::::::::::: Kate Hedstrom :::
#                                                                       :::
#  ROMS/TOMS Version 2.2 Master Makefile                                :::
#                                                                       :::
#  This makefile is designed to work only with GNU Make version 3.79 or :::
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
#::::::::::::::::::::::::::::::::::::::::::::::::::::: Hernan G. Arango :::


#--------------------------------------------------------------------------
#  Initialize some things.
#--------------------------------------------------------------------------

  clean_list := core *.o *.mod *.f90 lib*.a *.bak
  sources    := 
  path_srcs  := 
  libraries  :=

     objects  = $(subst .F,.o,$(sources)

#==========================================================================
#  Start of user-defined options. Modify macro variables: on is TRUE while
#  blank is FALSE.
#==========================================================================
#
#  Activate debugging compiler options:

       DEBUG :=

#  If parallel applications, use at most one of these definitions
#  (leave both definitions blank in serial applications):

         MPI := 
      OpenMP :=

#  If applicable, compile with the ARPACK library (GST analysis):

      ARPACK :=

#  If applicable, activate 64-bit compilation:

       LARGE := on

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
#     CYGWIN:                 g95, df
#     Darwin:                 f90
#     IRIX:                   f90
#     Linux:                  ifc, ifort, pgi, path, g95, mpif90
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

#--------------------------------------------------------------------------
#  Set ROMS/TOMS executable file name.
#--------------------------------------------------------------------------

BIN := $(BINDIR)/oceanS
ifdef MPI
  BIN := $(BINDIR)/oceanM
endif
ifdef OpenMP
  BIN := $(BINDIR)/oceanO
endif

#--------------------------------------------------------------------------
#  "uname -s" should return the OS or kernel name and "uname -m" should
#  return the CPU or hardware name. In practice the results can be pretty
#  flaky. Run the results through sed to convert "/" and " " to "-",
#  then apply platform-specific conversions.
#--------------------------------------------------------------------------

OS := $(shell uname -s | sed 's/[\/ ]/-/g')
OS := $(patsubst CYGWIN_%,CYGWIN,$(OS))
OS := $(patsubst sn%,UNICOS-sn,$(OS))

CPU := $(shell uname -m | sed 's/[\/ ]/-/g')

ifndef FORT
  $(error Variable FORT not set)
endif

ifneq "$(MAKECMDGOALS)" "clean"
  include Compilers/$(OS)-$(strip $(FORT)).mk
endif

#--------------------------------------------------------------------------
#  Pass the platform variables to the preprocessor as macros. Convert to
#  valid, upper-case identifiers.
#--------------------------------------------------------------------------

CPPFLAGS += -D$(shell echo ${OS} | tr "-" "_" | tr [a-z] [A-Z])
CPPFLAGS += -D$(shell echo ${CPU} | tr "-" "_" | tr [a-z] [A-Z])
CPPFLAGS += -D$(shell echo ${FORT} | tr "-" "_" | tr [a-z] [A-Z])

#--------------------------------------------------------------------------
#  Build target directories.
#--------------------------------------------------------------------------

.PHONY: all

all: $(BIN)

modules  := Nonlinear Utility Modules Drivers

includes := Include Nonlinear Drivers

vpath %.F $(modules)
vpath %.h $(includes)

include $(addsuffix /Module.mk,$(modules))

MDEPFLAGS += $(patsubst %,-I %,$(includes))

CPPFLAGS += $(patsubst %,-I%,$(includes))

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

depend:
	mv Compilers/MakeDepend Compilers/MakeDepend.orig
	Bin/sfmakedepend $(MDEPFLAGS) $(path_srcs) > Compilers/MakeDepend 

ifneq "$(MAKECMDGOALS)" "clean"
  -include Compilers/MakeDepend
endif

#--------------------------------------------------------------------------
#  Target to create ROMS/TOMS tar file.
#--------------------------------------------------------------------------

.PHONY: tarfile

tarfile:
		tar -cvf ocean-2_2.tar *

.PHONY: zipfile

zipfile:
		zip -r ocean-2_2.zip *


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
