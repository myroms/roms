# git $Id$
# svn $Id: roms_compiler_flags.cmake 1210 2024-01-03 22:03:03Z arango $
#:::::::::::::::::::::::::::::::::::::::::::::::::::::: David Robertson :::
# Copyright (c) 2002-2024 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.md                                                 :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
# CMake Fortran Compiler Flags.

set( my_os     "${CMAKE_SYSTEM_NAME}" )
set( my_cpu    "${CMAKE_SYSTEM_PROCESSOR}" )

###########################################################################
# Fortran
###########################################################################

if( CMAKE_Fortran_COMPILER_ID MATCHES "GNU" )
  include( compiler_flags_GNU_Fortran )
elseif( CMAKE_Fortran_COMPILER_ID MATCHES "Intel" )
  include( compiler_flags_Intel_Fortran )
else()
  message( STATUS "Fortran compiler with ID ${CMAKE_Fortran_COMPILER_ID} will be used with CMake default options")
endif()

string( TOUPPER ${CMAKE_BUILD_TYPE} build_flags )
message( STATUS "ROMS BUILD TYPE = ${CMAKE_BUILD_TYPE}" )
message( STATUS "ROMS COMPILER FLAGS = ${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_${build_flags}}" )

###########################################################################
# C-Preprocessor Definitions
###########################################################################

# Set "my_fflags" used in ROMS for compiling information.

if( ${CMAKE_BUILD_TYPE} MATCHES "Debug" )
  set( my_fflags "${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_DEBUG}" )
elseif( ${CMAKE_BUILD_TYPE} MATCHES "RelWithDebInfo" )
  set( my_fflags "${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_RELWITHDEBINFO}" )
else()
  set( my_fflags "${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_RELEASE}" )
endif()

string( TOUPPER ${my_os} OS )
string( TOUPPER ${my_cpu} CPU )
string( TOUPPER ${my_fort} FORT )
add_compile_definitions ( ${OS} ${CPU} ${FORT} )

