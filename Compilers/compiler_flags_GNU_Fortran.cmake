# git $Id$
# svn $Id: compiler_flags_GNU_Fortran.cmake 1051 2020-12-04 23:09:05Z arango $
#:::::::::::::::::::::::::::::::::::::::::::::::::::::: David Robertson :::
# Copyright (c) 2002-2020 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
# CMake Flags for the GNU Fortran Compiler.

###########################################################################
# FLAGS COMMON TO ALL BUILD TYPES
###########################################################################

set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -frepack-arrays")

###########################################################################
# RELEASE FLAGS
###########################################################################

set( CMAKE_Fortran_FLAGS_RELEASE "-O3 -ffast-math" )

###########################################################################
# DEBUG FLAGS
###########################################################################

set( CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -fbounds-check -fbacktrace -fcheck=all -fsanitize=address -fsanitize=undefined -finit-real=nan -ffpe-trap=invalid,zero,overflow" )

###########################################################################
# BIT REPRODUCIBLE FLAGS
###########################################################################

set( CMAKE_Fortran_FLAGS_BIT     "-O2 -fbounds-check" )

###########################################################################
# LINK FLAGS
###########################################################################

set( CMAKE_Fortran_LINK_FLAGS    "" )

###########################################################################
# ROMS Definitions
###########################################################################

set( my_fort   "gfortran" )
set( my_fc     "${CMAKE_Fortran_COMPILER}" )

if( ${CMAKE_BUILD_TYPE} EQUAL "Debug" )
  set( my_fflags "${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_DEBUG}" )
else()
  set( my_fflags "${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_RELEASE}" )
endif()

set( CPPFLAGS  "-P" "--traditional-cpp" "-w" )

###########################################################################


