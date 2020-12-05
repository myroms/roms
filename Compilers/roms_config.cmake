# git $Id$
# svn $Id: roms_config.cmake 1051 2020-12-04 23:09:05Z arango $
#:::::::::::::::::::::::::::::::::::::::::::::::::::::: David Robertson :::
# Copyright (c) 2002-2020 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
# CMake configuration for any ROMS application.

string( TOLOWER "${APP}.h" HEADER )
string( TOLOWER "${APP}" application )
add_definitions( -DROOT_DIR="${CMAKE_CURRENT_SOURCE_DIR}" -D${APP} -DHEADER="${HEADER}" )

# Do you want a shared or static library?

if( NOT DEFINED LIBTYPE )
  option( LIBSHARED "Build ROMS as a Shared Library" OFF )
  option( LIBSTATIC "Build ROMS as a Static Library" ON )
  Message( STATUS "No LIBTYPE chosen, defaulting to static. Valid types are SHARED, STATIC, and BOTH" )
else()
  string( TOLOWER "${LIBTYPE}" libtype )
  if( ${libtype} STREQUAL "shared" )
    option( LIBSHARED "Build ROMS as a Shared Library" ON )
    option( LIBSTATIC "Build ROMS as a Static Library" OFF )
    Message( STATUS "Building with shared ROMS library." )
  elseif( ${libtype} STREQUAL "static" )
    option( LIBSHARED "Build ROMS as a Shared Library" OFF )
    option( LIBSTATIC "Build ROMS as a Static Library" ON )
    Message( STATUS "Building with static ROMS library." )
  elseif( ${libtype} STREQUAL "both" )
    option( LIBSHARED "Build ROMS as a Shared Library" ON )
    option( LIBSTATIC "Build ROMS as a Static Library" ON )
    Message( STATUS "Building both static and shared ROMS libraries." )
  else()
    Message( FATAL_ERROR "Invalid value for LIBTYPE. Valid types are SHARED, STATIC, and BOTH" )
  endif()
endif()

# Do you want a ROMS executable?

if( NOT DEFINED ROMS_EXECUTABLE )
  option( ROMS_EXECUTABLE "Turn on/off building ROMS executable" ON )
  Message( STATUS "ROMS_EXECUTABLE not set, defaulting to ON. The ROMS executable will be built." )
elseif( ROMS_EXECUTABLE )
  option( ROMS_EXECUTABLE"Turn on/off building ROMS executable" ON )
  Message( STATUS "The ROMS executable will be built." )
elseif( NOT ROMS_EXECUTABLE )
  option( ROMS_EXECUTABLE "Turn on/off building ROMS executable" OFF )
  Message( STATUS "The ROMS executable will NOT be built." )
else()
  Message( FATAL_ERROR "Invalid value for ROMS_EXECUTABLE. Valid values are ON and OFF" )
endif()

# Does the application have a biological model?

if( DEFINED BIOLOGY )
  option( BIOLOGY "Turn on/off Biological Models" ON )
else()
  option( BIOLOGY "Turn on/off Biological Models" OFF )
endif()

#Does the application have a sediment model?

if( DEFINED SEDIMENT )
  option( SEDIMENT "Turn on/off Sediment Transport Model" ON )
else()
  option( SEDIMENT "Turn on/off Sediment Transport Model" OFF )
endif()

# Where is your <application>.h file?

set( HEADER_DIR  ${MY_HEADER_DIR} )

# Any custom analytical files?

if( NOT DEFINED MY_ANALYTICAL_DIR )
  set( ANALYTICAL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/ROMS/Functionals" )
  add_definitions( -DANALYTICAL_DIR="${ANALYTICAL_DIR}" )
  # no need for "include_directories" since it is already included
else()
  add_definitions( -DANALYTICAL_DIR="${MY_ANALYTICAL_DIR}" )
  include_directories( ${MY_ANALYTICAL_DIR} )
endif()

# Location(s) of ARPACK/PARPACK libraries. If both libarpack.a and libparpack.a
# both live in the same directory, leave ARPACK_LIBDIR blank.
#
# The decision about whether to use them in linking is computed below.
# This CMake setup will NOT build ARPACK/PARPACK for you.

if(DEFINED PARPACK_LIBDIR)
  set( PARPACK_LIBDIR "${PARPACK_LIBDIR}" )
else()
  set( PARPACK_LIBDIR "" )
endif()
if(DEFINED ARPACK_LIBDIR)
  set( ARPACK_LIBDIR "${ARPACK_LIBDIR}" )
else()
  set( ARPACK_LIBDIR "" )
endif()

Message( STATUS "PARPACK_LIBDIR = ${PARPACK_LIBDIR}" )
Message( STATUS "ARPACK_LIBDIR  = ${ARPACK_LIBDIR}" )

# Set ROMS SVN repository information.

set( SVN_URL "${MY_SVN_URL}" )
set( SVN_REV "${MY_SVN_REV}" )

set( ROMS_HEADER ${HEADER_DIR}/${HEADER} )

add_definitions(
  -DROMS_HEADER="${ROMS_HEADER}"
  -DSVN_URL="${SVN_URL}"
  -DSVN_REV="${SVN_REV}"
)

# Set ROMS Executable Name.

if( ${CMAKE_BUILD_TYPE} EQUAL "Debug" )
  set( BIN "romsG" )
elseif( MPI )
  set( BIN "romsM" )
  add_definitions( -DMPI )
else()
  set( BIN "romsS" )
endif()

if( MY_CPP_FLAGS )
  foreach( flag ${MY_CPP_FLAGS} )
    add_definitions( -D${flag} )
  endforeach()
  list( APPEND check_flags ${MY_CPP_FLAGS} )
endif()

list( APPEND USE_ADJOINT
      "ARRAY_MODES"
      "I4DVAR"
      "I4DVAR_ANA_SENSITIVITY"
      "NORMALIZATION"
      "RBL4DVAR"
      "R4DVAR"
      "SPLIT_I4DVAR"
      "SPLIT_RBL4DVAR"
      "SPLIT_R4DVAR"
      "SPLIT_SP4DVAR"
      "SP4DVAR"
)

list( APPEND USE_REPRESENTER
      "ARRAY_MODES"
      "R4DVAR"
      "R4DVAR_ANA_SENSITIVITY"
)

list( APPEND USE_TANGENT
      "ARRAY_MODES"
      "I4DVAR"
      "I4DVAR_ANA_SENSITIVITY"
      "NORMALIZATION"
      "RBL4DVAR"
      "R4DVAR"
      "SPLIT_I4DVAR"
      "SPLIT_RBL4DVAR"
      "SPLIT_R4DVAR"
      "SPLIT_SP4DVAR"
      "SP4DVAR"
)

# ARPACK/PARPACK

list( APPEND USE_ARPACK
      "ARRAY_MODES"
      "CLIPPING"
      "I4DVAR"
      "PROPAGATOR"
      "RBL4DVAR"
      "RBL4DVAR_ANA_SENSITIVITY"
      "R4DVAR"
      "R4DVAR_ANA_SENSITIVITY"
      "SP4DVAR"
)

foreach( flag ${check_flags} )

  # Be very careful with the IN_LIST syntax it must be no quotes around
  # the variable containing the string you are looking for and just the
  # name of the list (listename) not the list variable (${listname})

  if( NOT ADJOINT AND ( ${flag} IN_LIST USE_ADJOINT ) )
    option( ADJOINT "Turn on/off Adjoint Model" ON )
    Message( STATUS "Adjoint Model ENABLED" )
  endif()

  # If running ARRAY_MODES R4DVAR R4DVAR_ANA_SENSITIVITY enable
  # the Representer model.

  if( NOT REPRESENTER AND ( ${flag} IN_LIST USE_REPRESENTER ) )
    option( ADJOINT "Turn on/off Representer Model" ON )
    Message( STATUS "Representer Model ENABLED" )
  endif()

  if( NOT TANGENT AND ( ${flag} IN_LIST USE_TANGENT ) )
    option( TANGENT "Turn on/off Tangent Linear Model" ON )
    message( STATUS "Tangent Linear Model ENABLED" )
  endif()

  if( NOT ARPACK AND ( ${flag} IN_LIST USE_ARPACK ) )
    option( ARPACK "ARPACK" ON )
  endif()

endforeach()

