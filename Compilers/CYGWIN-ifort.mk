# svn $Id$
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Copyright (c) 2002-2007 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
#
# Include file for Compaq Visual Fortran compiler on Cygwin
# -------------------------------------------------------------------------
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

              BIN := $(BIN).exe

               FC := ifort
#          FFLAGS := /align /G7 /MD /iface:cvf
#          FFLAGS := /align /G7 /MD /assume:underscore /names:lowercase
           FFLAGS := /align /G7 /MD
              CPP := /usr/bin/cpp
         CPPFLAGS := -P -DCYGWIN -DCYGWIN_ifort -traditional
               LD := $(FC)
#         LDFLAGS := /link /nodefaultlib:libcmt  /nodefaultlib:msvcrt /stack:67108864
          LDFLAGS := /link /stack:67108864
               AR := ar
          ARFLAGS := r
               RM := rm -f
           RANLIB := ranlib
             PERL := perl

        MDEPFLAGS := --cpp --fext=f90 --file=-

#
# Library locations, can be overridden by environment variables.
# These are be specified in Unix form and will be converted as
# necessary to Windows form for Windows-native commands. The default
# values below assume that Cygwin mounts have been defined pointing to
# the NETCDF and ARPACK library locations.
#

   NETCDF_INCDIR ?= /netcdf-win32/include
   NETCDF_LIBDIR ?= /netcdf-win32/lib

         CPPFLAGS += -I$(NETCDF_INCDIR)
       NETCDF_LIB := $(NETCDF_LIBDIR)/libnetcdf.a

ifdef ARPACK
    ARPACK_LIBDIR ?= /arpack-win32/lib
       ARPACK_LIB := $(ARPACK_LIBDIR)/arpack.lib
endif

#
# Compiler flags
#

ifdef OpenMP
         CPPFLAGS += -D_OPENMP
           FFLAGS += /Qopenmp /Qopenmp_report1
endif

ifdef DEBUG
           FFLAGS += /debug:full /traceback /Od /Zi
else
           FFLAGS += /Ox
endif

ifdef MPI
       MPI_INCDIR ?= c:\\work\\models\\MPICH2\\include
       MPI_LIBDIR ?= c:\\work\\models\\MPICH2\\lib
       LIBS_WIN32 += "$(MPI_LIBDIR)\fmpich2.lib "
         CPPFLAGS += -DMPI -I$(MPI_INCDIR)
           FFLAGS += -I$(MPI_INCDIR)
endif

ifdef SWAN_COUPLE
       MCT_LIBDIR ?= c:\\work\\models\\MCT_v2.2\\mct
      MPEU_LIBDIR ?= c:\\work\\models\\MCT_v2.2\\mpeu
         CPPFLAGS += -traditional-cpp
           FFLAGS += -I$(MCT_LIBDIR) -I$(MPEU_LIBDIR) 
           FFLAGS += /fixed /noextend_source -assume:byterecl
       LIBS_WIN32 += "$(MCT_LIBDIR)\libmct.a" "$(MPEU_LIBDIR)\libmpeu.a"
endif

#
# For a Windows compiler, create variables pointing to the Windows
# file names needed when linking. Use of the "=" sign means that
# variables will be evaluated only when needed.
#

         BIN_WIN32 = "$$(cygpath --windows $(BIN))"
        LIBS_WIN32 += "$$(cygpath --windows $(NETCDF_LIB))"
        LIBS_WIN32 += "c:\cygwin\lib\gcc\i686-pc-mingw32\3.4.4\libgcc.a"
ifdef ARPACK
        LIBS_WIN32 += "$$(cygpath --windows $(ARPACK_LIB))"
endif

        LD_WINDOWS := on

#
# For a Windows compiler, override the compilation rule
#

%.o: %.f90
	$(FC) $(FFLAGS) -c $< /object:$@
