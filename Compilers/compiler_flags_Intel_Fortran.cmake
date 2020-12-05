# git $Id$
# svn $Id: compiler_flags_Intel_Fortran.cmake 1051 2020-12-04 23:09:05Z arango $
#:::::::::::::::::::::::::::::::::::::::::::::::::::::: David Robertson :::
# Copyright (c) 2002-2020 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
# CMake Flags for the Intel Fortran Compiler.

###########################################################################
# FLAGS COMMON TO ALL BUILD TYPES
###########################################################################

set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fp-model precise" )

###########################################################################
# RELEASE FLAGS
###########################################################################

set( CMAKE_Fortran_FLAGS_RELEASE "-ip -O3 -traceback -check uninit" )

###########################################################################
# DEBUG FLAGS
###########################################################################

set( CMAKE_Fortran_FLAGS_DEBUG   "-g -check all -check bounds -traceback -warn interfaces,nouncalled -gen-interfaces" )

###########################################################################
# BIT REPRODUCIBLE FLAGS
###########################################################################

set( CMAKE_Fortran_FLAGS_BIT     "-O2" )

###########################################################################
# LINK FLAGS
###########################################################################

set( CMAKE_Fortran_LINK_FLAGS    "" )

###########################################################################
# ROMS Definitions
###########################################################################

set( my_fort   "ifort" )
set( my_fc      "${CMAKE_Fortran_COMPILER}" )

if( ${CMAKE_BUILD_TYPE} EQUAL "Debug" )
  set( my_fflags "${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_DEBUG}" )
else()
  set( my_fflags "${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_RELEASE}" )
endif()

set( CPPFLAGS  "-P" "--traditional-cpp" "-w" )

###########################################################################
