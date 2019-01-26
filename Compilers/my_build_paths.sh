# svn $Id$
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Copyright (c) 2002-2019 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::: Hernan G. Arango :::
#                                                                       :::
# ROMS/TOMS Customized Compiling Libraries Script                       :::
#                                                                       :::
# This C-shell script sets the customized library paths needed by the   :::
# build script when the enviromental variable USE_MY_LIBS has a 'yes'   :::
# value.                                                                :::
#                                                                       :::
# For example, in build_roms.sh we have:                                :::
#                                                                       :::
#       if ($USE_MY_LIBS == 'yes') then                                 :::
#          source ${COMPILERS}/my_build_paths.sh                        :::
#       endif                                                           :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#--------------------------------------------------------------------------
# Add MPI library to compile.
#--------------------------------------------------------------------------
#
# Recall also that the MPI library comes in several flavors:
# MPICH, MPICH2, OpenMPI, etc.

# There are several MPI libraries available (MPICH, MPICH2, OpenMPI,
# etc.). Here, we set the desired  "mpif90" script to use during
# compilation. This only works if the make configuration file (say,
# Linux-ifort.mk) in the "Compilers" directory has the following
# definition for FC (Fortran Compiler) in the USE_MPI section:
#
#              FC := mpif90
#
# that is, "mpif90" defined without any path. Notice that the path
# where the MPI library is installed is computer dependent. Recall
# that you still need to use the appropriate "mpirun" to execute.

if ($?USE_MPIF90) then
  switch ($FORT)

    case "ifort"
      if ($which_MPI == "mpich" ) then
        setenv PATH /opt/intelsoft/mpich/bin:$PATH
      else if ($which_MPI == "mpich2" ) then
        setenv PATH /opt/intelsoft/mpich2/bin:$PATH
      else if ($which_MPI == "openmpi" ) then
        setenv PATH /opt/intelsoft/openmpi/bin:$PATH
      endif
    breaksw

    case "pgi"
      if ($which_MPI == "mpich" ) then
        setenv PATH /opt/pgisoft/mpich/bin:$PATH
      else if ($which_MPI == "mpich2" ) then
        setenv PATH /opt/pgisoft/mpich2/bin:$PATH
      else if ($which_MPI == "openmpi" ) then
        setenv PATH /opt/pgisoft/openmpi/bin:$PATH
      endif
    breaksw

    case "gfortran"
      if ($which_MPI == "mpich2" ) then
        setenv PATH /opt/gfortransoft/mpich2/bin:$PATH
      else if ($which_MPI == "openmpi" ) then
        setenv PATH /opt/gfortransoft/openmpi/bin:$PATH
      endif
    breaksw

  endsw
endif

#--------------------------------------------------------------------------
# Set libraries to compile and link.
#--------------------------------------------------------------------------
#
# The path of the libraries required by ROMS can be set here using
# environmental variables which take precedence to the values
# specified in the make macro definitions file (Compilers/*.mk).
# For most applications, only the location of the NetCDF library
# is needed during compilation.
#
# Notice that when the USE_NETCDF4 macro is activated, we need the
# serial or parallel version of the NetCDF-4/HDF5 library. The
# configuration script NF_CONFIG (available since NetCDF 4.0.1)
# is used to set up all the required libraries according to the
# installed options (openDAP, netCDF4/HDF5 file format). The
# parallel library uses the MPI-I/O layer (usually available
# in MPICH2 and OpenMPI) requiring compiling with the selected
# MPI library.
#
# In ROMS distributed-memory applications, you may use either the
# serial or parallel version of the NetCDF-4/HDF5 library. The
# parallel version is required when parallel I/O is activated
# (ROMS cpp option PARALLEL_IO and HDF5).
#
# However, in serial or shared-memory ROMS applications, we need
# to use the serial version of the NetCDF-4/HDF5 to avoid conflicts
# with the compiler. We cannot activate MPI constructs in serial
# or shared-memory ROMS code. Hybrid parallelism is not possible.

switch ($FORT)

# Intel Compiler:

  case "ifort"
    setenv ESMF_COMPILER        intelgcc
    if ($?USE_DEBUG) then
      setenv ESMF_BOPT          g
    else
      setenv ESMF_BOPT          O
    endif
    setenv ESMF_ABI             64
    setenv ESMF_COMM            ${which_MPI}
    setenv ESMF_SITE            default

    setenv ARPACK_LIBDIR        /opt/intelsoft/serial/ARPACK
    if ($?USE_MPI) then
      if ($which_MPI == "mpich" ) then
        setenv MCT_INCDIR       /opt/intelsoft/mpich/mct/include
        setenv MCT_LIBDIR       /opt/intelsoft/mpich/mct/lib
        setenv PARPACK_LIBDIR   /opt/intelsoft/mpich/PARPACK
      else if ($which_MPI == "mpich2" ) then
        setenv MCT_INCDIR       /opt/intelsoft/mpich2/mct/include
        setenv MCT_LIBDIR       /opt/intelsoft/mpich2/mct/lib
        setenv PARPACK_LIBDIR   /opt/intelsoft/mpich2/PARPACK
      else if ($which_MPI == "openmpi" ) then
        setenv MCT_INCDIR       /opt/intelsoft/openmpi/mct/include
        setenv MCT_LIBDIR       /opt/intelsoft/openmpi/mct/lib
        setenv PARPACK_LIBDIR   /opt/intelsoft/openmpi/PARPACK
      endif
    endif

    if ($?USE_NETCDF4) then
      if ($?USE_PARALLEL_IO && $?USE_MPI) then
        if ($which_MPI == "mpich" ) then
          setenv ESMF_DIR       /opt/intelsoft/mpich/esmf_nc4
          setenv NF_CONFIG      /opt/intelsoft/mpich/netcdf4/bin/nf-config
          setenv NETCDF_INCDIR  /opt/intelsoft/mpich/netcdf4/include
        else if ($which_MPI == "mpich2" ) then
          setenv ESMF_DIR       /opt/intelsoft/mpich2/esmf_nc4
          setenv NF_CONFIG      /opt/intelsoft/mpich2/netcdf4/bin/nf-config
          setenv NETCDF_INCDIR  /opt/intelsoft/mpich2/netcdf4/include
        else if ($which_MPI == "openmpi" ) then
          setenv ESMF_DIR       /opt/intelsoft/openmpi/esmf_nc4
          setenv NF_CONFIG      /opt/intelsoft/openmpi/netcdf4/bin/nf-config
          setenv NETCDF_INCDIR  /opt/intelsoft/openmpi/netcdf4/include
        endif
      else
        if ($which_MPI == "mpich" ) then
          setenv ESMF_DIR       /opt/intelsoft/mpich/esmf_nc4
        else if ($which_MPI == "mpich2" ) then
          setenv ESMF_DIR       /opt/intelsoft/mpich2/esmf_nc4
        else if ($which_MPI == "openmpi" ) then
          setenv ESMF_DIR       /opt/intelsoft/openmpi/esmf_nc4
        endif
        setenv NF_CONFIG        /opt/intelsoft/serial/netcdf4/bin/nf-config
        setenv NETCDF_INCDIR    /opt/intelsoft/serial/netcdf4/include
      endif
    else
      if ($which_MPI == "mpich" ) then
        setenv ESMF_DIR         /opt/intelsoft/mpich/esmf_nc3
      else if ($which_MPI == "mpich2" ) then
        setenv ESMF_DIR         /opt/intelsoft/mpich2/esmf_nc3
      else if ($which_MPI == "openmpi" ) then
        setenv ESMF_DIR         /opt/intelsoft/openmpi/esmf_nc3
      endif
      setenv NETCDF_INCDIR      /opt/intelsoft/serial/netcdf3/include
      setenv NETCDF_LIBDIR      /opt/intelsoft/serial/netcdf3/lib
    endif

    if ($?USE_HDF5) then
      if ($?USE_PARALLEL_IO && $?USE_MPI) then
        if ( $which_MPI == "mpich" ) then
          setenv HDF5_LIBDIR     /opt/intelsoft/mpich/hdf5/lib
          setenv HDF5_INCDIR     /opt/intelsoft/mpich/hdf5/include
        else if ( $which_MPI == "mpich2" ) then
          setenv HDF5_LIBDIR     /opt/intelsoft/mpich2/hdf5/lib
          setenv HDF5_INCDIR     /opt/intelsoft/mpich2/hdf5/include
        else if ( $which_MPI == "openmpi" ) then
          setenv HDF5_LIBDIR     /opt/intelsoft/openmpi/hdf5/lib
          setenv HDF5_INCDIR     /opt/intelsoft/openmpi/hdf5/include
        endif
      else
        setenv HDF5_LIBDIR       /opt/intelsoft/serial/hdf5/lib
        setenv HDF5_INCDIR       /opt/intelsoft/serial/hdf5/include
      endif
    endif

  breaksw

# PGI Compiler:

  case "pgi"
    setenv ESMF_COMPILER        pgi
    if ($?USE_DEBUG) then
      setenv ESMF_BOPT          g
    else
      setenv ESMF_BOPT          O
    endif
    setenv ESMF_ABI             64
    setenv ESMF_COMM            ${which_MPI}
    setenv ESMF_SITE            default

    setenv ARPACK_LIBDIR        /opt/pgisoft/serial/ARPACK
    if ($?USE_MPI) then
      if ($which_MPI == "mpich" ) then
        setenv MCT_INCDIR       /opt/pgisoft/mpich/mct/include
        setenv MCT_LIBDIR       /opt/pgisoft/mpich/mct/lib
        setenv PARPACK_LIBDIR   /opt/pgisoft/mpich/PARPACK
      else if ($which_MPI == "mpich2" ) then
        setenv MCT_INCDIR       /opt/pgisoft/mpich2/mct/include
        setenv MCT_LIBDIR       /opt/pgisoft/mpich2/mct/lib
        setenv PARPACK_LIBDIR   /opt/pgisoft/mpich2/PARPACK
      else if ($which_MPI == "openmpi" ) then
        setenv MCT_INCDIR       /opt/pgisoft/openmpi/mct/include
        setenv MCT_LIBDIR       /opt/pgisoft/openmpi/mct/lib
        setenv PARPACK_LIBDIR   /opt/pgisoft/openmpi/PARPACK
      endif
    endif

    if ($?USE_NETCDF4) then
      if ($?USE_PARALLEL_IO && $?USE_MPI) then
        if ($which_MPI == "mpich" ) then
          setenv ESMF_DIR       /opt/pgisoft/mpich/esmf_nc4
          setenv NF_CONFIG      /opt/pgisoft/mpich/netcdf4/bin/nf-config
          setenv NETCDF_INCDIR  /opt/pgisoft/mpich/netcdf4/include
        else if ($which_MPI == "mpich2" ) then
          setenv ESMF_DIR       /opt/pgisoft/mpich2/esmf_nc4
          setenv NF_CONFIG      /opt/pgisoft/mpich2/netcdf4/bin/nf-config
          setenv NETCDF_INCDIR  /opt/pgisoft/mpich2/netcdf4/include
        else if ($which_MPI == "openmpi" ) then
          setenv ESMF_DIR       /opt/pgisoft/openmpi/esmf_nc4
          setenv NF_CONFIG      /opt/pgisoft/openmpi/netcdf4/bin/nf-config
          setenv NETCDF_INCDIR  /opt/pgisoft/openmpi/netcdf4/include
        endif
      else
        if ($which_MPI == "mpich" ) then
          setenv ESMF_DIR       /opt/pgisoft/mpich/esmf_nc4
        else if ($which_MPI == "mpich2" ) then
          setenv ESMF_DIR       /opt/pgisoft/mpich2/esmf_nc4
        else if ($which_MPI == "openmpi" ) then
          setenv ESMF_DIR       /opt/pgisoft/openmpi/esmf_nc4
        endif
        setenv NF_CONFIG        /opt/pgisoft/serial/netcdf4/bin/nf-config
        setenv NETCDF_INCDIR    /opt/pgisoft/serial/netcdf4/include
      endif
    else
      if ($which_MPI == "mpich" ) then
        setenv ESMF_DIR         /opt/pgisoft/mpich/esmf_nc3
      else if ($which_MPI == "mpich2" ) then
        setenv ESMF_DIR         /opt/pgisoft/mpich2/esmf_nc3
      else if ($which_MPI == "openmpi" ) then
        setenv ESMF_DIR         /opt/pgisoft/openmpi/esmf_nc3
      endif
      setenv NETCDF_INCDIR      /opt/pgisoft/serial/netcdf3/include
      setenv NETCDF_LIBDIR      /opt/pgisoft/serial/netcdf3/lib
    endif

    if ($?USE_HDF5) then
      if ($?USE_PARALLEL_IO && $?USE_MPI) then
        if ( $which_MPI == "mpich" ) then
          setenv HDF5_LIBDIR     /opt/pgisoft/mpich/hdf5/lib
          setenv HDF5_INCDIR     /opt/pgisoft/mpich/hdf5/include
        else if ( $which_MPI == "mpich2" ) then
          setenv HDF5_LIBDIR     /opt/pgisoft/mpich2/hdf5/lib
          setenv HDF5_INCDIR     /opt/pgisoft/mpich2/hdf5/include
        else if ( $which_MPI == "openmpi" ) then
          setenv HDF5_LIBDIR     /opt/pgisoft/openmpi/hdf5/lib
          setenv HDF5_INCDIR     /opt/pgisoft/openmpi/hdf5/include
        endif
      else
        setenv HDF5_LIBDIR       /opt/pgisoft/serial/hdf5/lib
        setenv HDF5_INCDIR       /opt/pgisoft/serial/hdf5/include
      endif
    endif

  breaksw

# GNU Compiler:

  case "gfortran"
    setenv ESMF_COMPILER        gfortran
    if ($?USE_DEBUG) then
      setenv ESMF_BOPT          g
    else
      setenv ESMF_BOPT          O
    endif
    setenv ESMF_ABI             64
    setenv ESMF_COMM            ${which_MPI}
    setenv ESMF_SITE            default

    setenv ARPACK_LIBDIR        /opt/gfortransoft/serial/ARPACK
    if ($?USE_MPI) then
      if ($which_MPI == "mpich2" ) then
        setenv MCT_INCDIR       /opt/gfortransoft/mpich2/mct/include
        setenv MCT_LIBDIR       /opt/gfortransoft/mpich2/mct/lib
        setenv PARPACK_LIBDIR   /opt/gfortransoft/mpich2/PARPACK
      else if ($which_MPI == "openmpi" ) then
        setenv MCT_INCDIR       /opt/gfortransoft/openmpi/mct/include
        setenv MCT_LIBDIR       /opt/gfortransoft/openmpi/mct/lib
        setenv PARPACK_LIBDIR   /opt/gfortransoft/openmpi/PARPACK
      endif
    endif

    if ($?USE_NETCDF4) then
      if ($?USE_PARALLEL_IO && $?USE_MPI) then
        if ($which_MPI == "mpich2" ) then
          setenv ESMF_DIR       /opt/gfortransoft/mpich2/esmf_nc4
          setenv NF_CONFIG      /opt/gfortransoft/mpich2/netcdf4/bin/nf-config
          setenv NETCDF_INCDIR  /opt/gfortransoft/mpich2/netcdf4/include
        else if ($which_MPI == "openmpi" ) then
          setenv ESMF_DIR       /opt/gfortransoft/openmpi/esmf_nc4
          setenv NF_CONFIG      /opt/gfortransoft/openmpi/netcdf4/bin/nf-config
          setenv NETCDF_INCDIR  /opt/gfortransoft/openmpi/netcdf4/include
        endif
      else
        if ($which_MPI == "mpich2" ) then
          setenv ESMF_DIR       /opt/gfortransoft/mpich2/esmf_nc4
        else if ($which_MPI == "openmpi" ) then
            setenv ESMF_DIR       /opt/gfortransoft/openmpi/esmf_nc4
        endif
        setenv NF_CONFIG        /opt/gfortransoft/serial/netcdf4/bin/nf-config
        setenv NETCDF_INCDIR    /opt/gfortransoft/serial/netcdf4/include
      endif
    else
      if ($which_MPI == "mpich2" ) then
        setenv ESMF_DIR         /opt/gfortransoft/mpich2/esmf_nc3
      else if ($which_MPI == "openmpi" ) then
        setenv ESMF_DIR         /opt/gfortransoft/openmpi/esmf_nc3
      endif
      setenv NETCDF_INCDIR      /opt/gfortransoft/serial/netcdf3/include
      setenv NETCDF_LIBDIR      /opt/gfortransoft/serial/netcdf3/lib
    endif

    if ($?USE_HDF5) then
      if ($?USE_PARALLEL_IO && $?USE_MPI) then
        if ( $which_MPI == "mpich2" ) then
          setenv HDF5_LIBDIR     /opt/gfortransoft/mpich2/hdf5/lib
          setenv HDF5_INCDIR     /opt/gfortransoft/mpich2/hdf5/include
        else if ( $which_MPI == "openmpi" ) then
          setenv HDF5_LIBDIR     /opt/gfortransoft/openmpi/hdf5/lib
          setenv HDF5_INCDIR     /opt/gfortransoft/openmpi/hdf5/include
        endif
      else
        setenv HDF5_LIBDIR       /opt/gfortransoft/serial/hdf5/lib
        setenv HDF5_INCDIR       /opt/gfortransoft/serial/hdf5/include
      endif
    endif

  breaksw

endsw
