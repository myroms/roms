# git $Id$
# svn $Id: CMakeLists.txt 1051 2020-12-04 23:09:05Z arango $
#:::::::::::::::::::::::::::::::::::::::::::::::::::::: David Robertson :::
# Copyright (c) 2002-2020 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
# Top-Level CMake File Definitions.

cmake_minimum_required( VERSION 3.12.0 FATAL_ERROR )

set( CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/Compilers;${CMAKE_MODULE_PATH})
Message( STATUS "module path = ${CMAKE_MODULE_PATH}" )

project(roms VERSION 3.9 LANGUAGES Fortran)

set( ECBUILD_DEFAULT_BUILD_TYPE Release )
set( ENABLE_OS_TESTS           OFF CACHE BOOL "Disable OS tests" FORCE )
set( ENABLE_LARGE_FILE_SUPPORT OFF CACHE BOOL "Disable testing of large file support" FORCE )
set( ENABLE_MPI ON CACHE BOOL "Compile with MPI" )

include( ecbuild_system NO_POLICY_SCOPE )

ecbuild_requires_macro_version( 2.5 )

###########################################################################
# Project
###########################################################################

ecbuild_declare_project()

ecbuild_enable_fortran( REQUIRED )

set( ROMS_LINKER_LANGUAGE Fortran )

###########################################################################
# Dependencies
###########################################################################

# NetCDF Library.

find_package( NetCDF REQUIRED COMPONENTS Fortran )
include_directories( ${NETCDF_INCLUDE_DIRS} )

# MPI Library.

ecbuild_find_mpi( COMPONENTS Fortran REQUIRED )
ecbuild_include_mpi()
link_libraries(${MPI_Fortran_LIBRARIES})

###########################################################################
# ROMS
###########################################################################

include( roms_compiler_flags )
include( roms_config )
include( roms_functions )

add_subdirectory( Master )
add_subdirectory( ROMS/Drivers )
add_subdirectory( ROMS/Functionals )
add_subdirectory( ROMS/Modules )
add_subdirectory( ROMS/Nonlinear )
add_subdirectory( ROMS/Utility )

include_directories(
  ${CMAKE_CURRENT_BINARY_DIR}
  "/usr/local/include"
  "$ENV{NETCDF}/include"
  "Master"
  "ROMS/Drivers"
  "ROMS/Functionals"
  "ROMS/Include"
  "ROMS/Modules"
  "ROMS/Nonlinear"
  "ROMS/Utility"
)

link_directories(
  ${CMAKE_CURRENT_BINARY_DIR}
  "/usr/local/lib"
  "$ENV{NETCDF}/lib"
)

list(APPEND srcs
     ${Master_files}
     ${ROMS_Drivers_files}
     ${ROMS_Functionals_files}
     ${ROMS_Modules_files}
     ${ROMS_Nonlinear_files}
     ${ROMS_Utility_files}
)

if( BIOLOGY )
  add_subdirectory( ROMS/Nonlinear/Biology )
  include_directories( "ROMS/Nonlinear/Biology" )
  list( APPEND srcs
        ${ROMS_Nonlinear_Biology_files}
  )
endif()

if( SEDIMENT )
  add_subdirectory( ROMS/Nonlinear/Sediment )
  include_directories( "ROMS/Nonlinear/Sediment" )
  list( APPEND srcs
        ${ROMS_Nonlinear_Sediment_files}
  )
endif()

if( ADJOINT )
  if( BIOLOGY )
    add_subdirectory( ROMS/Adjoint/Biology )
    include_directories( "ROMS/Adjoint/Biology" )
    list( APPEND srcs
          ${ROMS_Adjoint_Biology_files}
    )
  endif()

  add_subdirectory( ROMS/Adjoint )
  include_directories( "ROMS/Adjoint" )
  list( APPEND srcs
        ${ROMS_Adjoint_files}
  )
endif()

if( REPRESENTER )
  if( BIOLOGY )
    add_subdirectory( ROMS/Representer/Biology )
    include_directories( "ROMS/Representer/Biology" )
    list( APPEND srcs
          ${ROMS_Representer_Biology_files}
    )
  endif()

  add_subdirectory( ROMS/Representer )
  include_directories( "ROMS/Representer" )
  list( APPEND srcs
        ${ROMS_Representer_files}
  )
endif()

if( TANGENT )
  if( BIOLOGY )
    add_subdirectory( ROMS/Tangent/Biology )
    include_directories( "ROMS/Tangent/Biology" )
    list( APPEND srcs
          ${ROMS_Tangent_Biology_files}
    )
  endif()

  add_subdirectory( ROMS/Tangent )
  include_directories( "ROMS/Tangent" )
  list( APPEND srcs
        ${ROMS_Tangent_files}
  )
endif()

# Special CPP definitions for "mod_strings.F".

set_property(
  SOURCE ROMS/Modules/mod_strings.F
  APPEND_STRING PROPERTY COMPILE_DEFINITIONS
    MY_OS='${my_os}' MY_CPU='${my_cpu}' MY_FORT='${my_fort}'
    MY_FC='${my_fc}'
    MY_FFLAGS='${my_fflags}'
)

###########################################################################
#
###########################################################################


# Remove the .F files from "DependInfo.cmake" for the Objects target.

add_custom_command(
  OUTPUT fixdep
  COMMAND ${PERL} ${CMAKE_CURRENT_SOURCE_DIR}/ROMS/Bin/FixDependInfo.pl ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/Objects.dir/DependInfo.cmake
  COMMAND touch fixdep
  VERBATIM
)

# Custom target needed to make the custom command above run

add_custom_target(fix ALL
  DEPENDS fixdep
)

# The "preprocess_fortran" function tells CMake how to create .f90 files.
# It allows CMake crippled Fortran depency tracker to look for the right
# modules using the pre-processed .f90 files instead of the .F files.
#
# The dependency generator is easily confused by CPP if-directives.
# This also allows us to feed .f90 files to CMake as the sources to all
# the ROMS code and eliminate compiler assumptions related to .F files.
#
# The function "preprocess_fortran" is defined in the "roms_functions.cmake"
# file in the "Compilers" directory.

preprocess_fortran( ${srcs} )
set( All_f90s "${f90srcs}" )

# Adding all the sources as an OBJECT library prevents CMake from trying
# to link an executable or collect the object files into a library.
#
# The libXXX.(a|so) libraries will be created below.

add_library( Objects OBJECT ${All_f90s} )

# Add dependency on custom target "fix" so "DependInfo.cmake" is cleaned up
# before the the dependencies are calculated.

add_dependencies( Objects fix )

# Shared libraries need PIC. There is perfomance degradation to position
# independent code so we only set this if shared libraries are requested.

set( deplib )

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_BINARY_DIR} CACHE PATH "..." FORCE)
endif()

if( LIBSHARED )
  set_property(TARGET Objects PROPERTY POSITION_INDEPENDENT_CODE 1)
  add_library( ROMS_shared SHARED $<TARGET_OBJECTS:Objects> )
  set_target_properties(ROMS_shared PROPERTIES OUTPUT_NAME ROMS)
  install( TARGETS ROMS_shared LIBRARY )

  # This is needed to add to the depency list of the ROMS executable
  # since it is linking with libROMS.so instead of libROMS_shared.so
  # The executable doesn't see "ROMS" as a target generated by this
  # project so it starts linking before the library is fully assembled.

  list( APPEND deplib ROMS_shared )
endif()

if( LIBSTATIC )
  add_library( ROMS_static STATIC $<TARGET_OBJECTS:Objects> )
  set_target_properties(ROMS_static PROPERTIES OUTPUT_NAME ROMS)
  install( TARGETS ROMS_static ARCHIVE )

  # This is needed to add to the depency list of the ROMS executable
  # since it is linking with libROMS.a instead of libROMS_static.a
  # The executable doesn't see "ROMS" as a target generated by this
  # project so it starts linking before the library is fully assembled.

  list( APPEND deplib ROMS_static )
endif()

# If PARPACK/ARPACK is needed add directory for the linker

if( ARPACK )
  link_directories( ${PARPACK_LIBDIR} ${ARPACK_LIBDIR} )
endif()

# There needs to be at least one source file for the "add_executable" line
# so "Master/master.F" is removed from "Master/CMakeLists.txt" and added
# manually here.
#
# I would have done all the files in "Master" here to simplify things
# but "ocean_control.F", "coupler.F", and "propagator.F" include
# modules needed in other parts of the code outside "Master".

if( ROMS_EXECUTABLE )
  preprocess_fortran( Master/master.F )
  set( master_f90 "${f90srcs}" )
  add_executable( "${BIN}" ${master_f90} )
  add_dependencies( "${BIN}" Objects ${deplib} )

  if( LIBSHARED AND LIBSTATIC )
    Message( STATUS "Both LIBSHARED and LIBSTATIC are enabled, the linker will probably default to the shared ROMS library." )
    target_link_libraries( "${BIN}" ROMS netcdff netcdf )
  elseif( LIBSHARED )
    Message( STATUS "Linking ROMS using shared ROMS library." )
    target_link_libraries( "${BIN}" ROMS netcdff netcdf )
  elseif( LIBSTATIC )
    Message( STATUS "Linking ROMS using static ROMS library." )
    target_link_libraries( "${BIN}" ROMS netcdff netcdf )
  else()
    Message( FATAL_ERROR "Unknown configuration!" )
  endif()

  # If ARPACK is needed, add to linking library list. The path to the library is set above.

  if( ARPACK )
    target_link_libraries( "${BIN}" parpack arpack )
  endif()

  install( TARGETS "${BIN}" RUNTIME )

endif()
