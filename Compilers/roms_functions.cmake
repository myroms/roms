# git $Id$
# svn $Id: roms_functions.cmake 1051 2020-12-04 23:09:05Z arango $
#:::::::::::::::::::::::::::::::::::::::::::::::::::::: David Robertson :::
# Copyright (c) 2002-2020 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
# Functions used in CMake to overcome its rudimentary capabilities.

find_program(
  PERL
  NAMES "perl"
  DOC "Perl command"
)

find_program(
  CPP_CLEAN
  NAMES "cpp_clean"
  HINTS "${CMAKE_SOURCE_DIR}/ROMS/Bin/"
  DOC "ROMS CPP Clean command"
)

find_program(
  CPP_EXECUTABLE
  NAMES "cpp"
  DOC "C-preprocessor command"
)

# This function teaches CMake how to use CPP to create .f90 files from
# ROMS .F source files. It is imperative that we generate .f90 files
# BEFORE CMake attempts to determine Fortran dependencies because
# its Fortran dependency tracker is inconsistent and easily confused
# by C-preprocessor if-directives.
#
# First, make the directory for the .f90 files.

file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/f90")

function(preprocess_fortran f90srcs)

   # Determine include directories.

   get_directory_property(incdirs INCLUDE_DIRECTORIES)
   set(incflags)
   foreach(i ${incdirs})
     list(APPEND incflags "-I${i}")
   endforeach()

   # Set-up "defines" (not configuration specific).

   get_directory_property(defines COMPILE_DEFINITIONS)
   set(defflags)
   foreach(d ${defines})
     list(APPEND defflags "-D${d}")
   endforeach()

   # Loop over all arguments.

   set(fsrcs)
   foreach(f ${ARGV})

     # Is it a Fortran file?

     if(f MATCHES "\\.[Ff](9[05])?")
       message(STATUS "Got fortran file: ${f}")

       # Get special definitions for file, if any.

       set(fdefflags)
       if(${f} MATCHES ".*mod_strings.*")
         get_property(fdefines
                      SOURCE ${f}
                      PROPERTY COMPILE_DEFINITIONS
         )
         foreach(fd ${fdefines})
           list(APPEND fdefflags "-D${fd}")
         endforeach()
       endif()

       # Construct output filename.

       if(NOT IS_ABSOLUTE "${f}")
         get_filename_component(f "${f}" ABSOLUTE)
       endif()
       file(RELATIVE_PATH r "${CMAKE_CURRENT_SOURCE_DIR}" "${f}")
       get_filename_component(e "${r}" EXT)
       get_filename_component(n "${r}" NAME_WE)
       get_filename_component(p "${r}" PATH)
       set(of "${CMAKE_CURRENT_BINARY_DIR}/f90/${n}.f90")
       message(STATUS "Output name: ${of}")

       # Preprocess the thing.

       add_custom_command(OUTPUT "${of}"
         COMMAND ${CPP_EXECUTABLE} ${CPPFLAGS} ${defflags} ${fdefflags} ${incflags} "${f}" > "${of}"
         COMMAND ${CPP_CLEAN} "${of}"
         IMPLICIT_DEPENDS Fortran "${f}"
         COMMENT "Preprocessing ${f}"
         COMMAND_EXPAND_LISTS
         VERBATIM
       )
       list(APPEND fsrcs "${of}")
     else()
       list(APPEND fsrcs "${f}")
     endif()
   endforeach()

   # Return the (preprocessed) sources.

   set(f90srcs "${fsrcs}" PARENT_SCOPE)
endfunction()

# UNDER CONSTRUCTION. This function is not working quite right yet.

function(link_netcdf)
  set(NF-CONFIG $ENV{NETCDF}/bin/nf-config)
  execute_process(COMMAND ${NF-CONFIG} --flibs
                  OUTPUT_VARIABLE tmp
  )

  # Strip out massive trailing whitespace.

  string(STRIP "${tmp}" linkline)
  string(REGEX MATCHALL "-L[^ \t]*" ldirs ${linkline})
  string(REGEX MATCHALL "-l[^ \t]*" libs ${linkline})
  set(netcdf_ldirs "${ldirs}" PARENT_SCOPE)
  set(netcdf_libs "${libs}" PARENT_SCOPE)
  Message(STATUS "netcdf_libs = ${libs} netcdf_ldirs = ${ldirs}")
endfunction()
