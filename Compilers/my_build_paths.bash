# svn $Id$
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Copyright (c) 2002-2019 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::: Hernan G. Arango :::
#                                                                       :::
# ROMS/TOMS Customized Compiling Libraries Script                       :::
#                                                                       :::
# This bash script sets the customized library paths needed by the      :::
# build script when the enviromental variable USE_MY_LIBS has a 'yes'   :::
# value.                                                                :::
#                                                                       :::
# For example, in build_roms.bash we have:                              :::
#                                                                       :::
#       if [ "${USE_MY_LIBS}" = "yes" ]; then                           :::
#         source ${COMPILERS}/my_build_paths.bash                       :::
#       fi                                                              :::
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

if [ -n "${USE_MPIF90:+1}" ]; then
  case "$FORT" in
    ifort )
      if [ "${which_MPI}" = "mpich" ]; then
        export PATH=/opt/intelsoft/mpich/bin:$PATH
      elif [ "${which_MPI}" = "mpich2" ]; then
        export PATH=/opt/intelsoft/mpich2/bin:$PATH
      elif [ "${which_MPI}" = "openmpi" ]; then
        export PATH=/opt/intelsoft/openmpi/bin:$PATH
      fi
      ;;

    pgi )
      if [ "${which_MPI}" = "mpich" ]; then
        export PATH=/opt/pgisoft/mpich/bin:$PATH
      elif [ "${which_MPI}" = "mpich2" ]; then
        export PATH=/opt/pgisoft/mpich2/bin:$PATH
      elif [ "${which_MPI}" = "openmpi" ]; then
        export PATH=/opt/pgisoft/openmpi/bin:$PATH
      fi
      ;;

    gfortran )
      if [ "${which_MPI}" = "mpich2" ]; then
        export PATH=/opt/gfortransoft/mpich2/bin:$PATH
      elif [ "${which_MPI}" = "openmpi" ]; then
        export PATH=/opt/gfortransoft/openmpi/bin:$PATH
      fi
      ;;

  esac
fi

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

case "$FORT" in
  ifort )
    export           MPI_ROOT=""
    export      ESMF_COMPILER=intelgcc
    if [ -n "${USE_DEBUG:+1}" ]; then
      export        ESMF_BOPT=g
    else
      export        ESMF_BOPT=O
    fi
    export           ESMF_ABI=64
    export          ESMF_COMM=${which_MPI}
    export          ESMF_SITE=default

    export       ARPACK_LIBDIR=/opt/intelsoft/serial/ARPACK
    if [ -n "${USE_MPI:+1}" ]; then
      if [ "${which_MPI}" = "mpich" ]; then
        export       MPI_ROOT=/opt/intelsoft/mpich
      elif [ "${which_MPI}" = "mpich2" ]; then
        export       MPI_ROOT=/opt/intelsoft/mpich2
      elif [ "${which_MPI}" = "openmpi" ]; then
        export       MPI_ROOT=/opt/intelsoft/openmpi
      fi
      export       MCT_INCDIR=${MPI_ROOT}/mct/include
      export       MCT_LIBDIR=${MPI_ROOT}/mct/lib
      export   PARPACK_LIBDIR=${MPI_ROOT}/PARPACK
    fi

    if [ -n "${USE_NETCDF4:+1}" ]; then
      if [ -n "${USE_PARALLEL_IO:+1}" ] && [ -n "${USE_MPI:+1}" ]; then
        export       ESMF_DIR=${MPI_ROOT}/esmf_nc4
        export         NETCDF=${MPI_ROOT}/netcdf4
        export      NF_CONFIG=${NETCDF}/bin/nf-config
        export  NETCDF_INCDIR=${NETCDF}/include
        export        NETCDF4=1
      else
        export       ESMF_DIR=${MPI_ROOT}/esmf_nc4
        export         NETCDF=/opt/intelsoft/serial/netcdf4
        export      NF_CONFIG=${NETCDF}/bin/nf-config
        export  NETCDF_INCDIR=${NETCDF}/include
        export        NETCDF4=1
      fi
    else
      export         ESMF_DIR=${MPI_ROOT}/esmf_nc3
      export           NETCDF=/opt/intelsoft/serial/netcdf3
      export    NETCDF_INCDIR=${NETCDF}/include
      export    NETCDF_LIBDIR=${NETCDF}/lib
      export   NETCDF_classic=1
    fi

    if [ -n "${USE_HDF5:+1}" ]; then
      if [ -n "${USE_PARALLEL_IO:+1}" ] && [ -n "${USE_MPI:+1}" ]; then
        export           HDF5=${MPI_ROOT}/hdf5
        export    HDF5_LIBDIR=${HDF5}/lib
        export    HDF5_INCDIR=${HDF5}/include
      else
        export           HDF5=/opt/intelsoft/serial/hdf5
        export    HDF5_LIBDIR=${HDF5}/lib
        export    HDF5_INCDIR=${HDF5}/include
      fi
    fi
    ;;

  pgi )
    export           MPI_ROOT=""
    export      ESMF_COMPILER=pgi
    if [ -n "${USE_DEBUG:+1}" ]; then
      export        ESMF_BOPT=g
    else
      export        ESMF_BOPT=O
    fi
    export           ESMF_ABI=64
    export          ESMF_COMM=${which_MPI}
    export          ESMF_SITE=default

    export      ARPACK_LIBDIR=/opt/pgisoft/serial/ARPACK
    if [ -n "${USE_MPI:+1}" ]; then
      if [ "${which_MPI}" = "mpich" ]; then
        export       MPI_ROOT=/opt/pgisoft/mpich
      elif [ "${which_MPI}" = "mpich2" ]; then
        export       MPI_ROOT=/opt/pgisoft/mpich2
      elif [ "${which_MPI}" = "openmpi" ]; then
        export       MPI_ROOT=/opt/pgisoft/openmpi
      fi
      export       MCT_INCDIR=${MPI_ROOT}/mct/include
      export       MCT_LIBDIR=${MPI_ROOT}/mct/lib
      export   PARPACK_LIBDIR=${MPI_ROOT}/PARPACK
    fi

    if [ -n "${USE_NETCDF4:+1}" ]; then
      if [ -n "${USE_PARALLEL_IO:+1}" ] && [ -n "${USE_MPI:+1}" ]; then
        export       ESMF_DIR=${MPI_ROOT}/esmf_nc4
        export         NETCDF=${MPI_ROOT}/netcdf4
        export      NF_CONFIG=${NETCDF}/bin/nf-config
        export  NETCDF_INCDIR=${NETCDF}/include
        export        NETCDF4=1
      else
        export       ESMF_DIR=${MPI_ROOT}/esmf_nc4
        export         NETCDF=/opt/pgisoft/serial/netcdf4
        export      NF_CONFIG=${NETCDF}/bin/nf-config
        export  NETCDF_INCDIR=${NETCDF}/include
        export        NETCDF4=1
      fi
    else
      export         ESMF_DIR=${MPI_ROOT}/esmf_nc3
      export           NETCDF=/opt/pgisoft/serial/netcdf3
      export    NETCDF_INCDIR=${NETCDF}/include
      export    NETCDF_LIBDIR=${NETCDF}/lib
      export   NETCDF_classic=1
    fi

    if [ -n "${USE_HDF5:+1}" ]; then
      if [ -n "${USE_PARALLEL_IO:+1}" ] && [ -n "${USE_MPI:+1}" ]; then
        export           HDF5=${MPI_ROOT}/hdf5
        export    HDF5_LIBDIR=${HDF5}/lib
        export    HDF5_INCDIR=${HDF5}/include
      else
        export           HDF5=/opt/pgisoft/serial/hdf5
        export    HDF5_LIBDIR=${HDF5}/lib
        export    HDF5_INCDIR=${HDF5}/include
      fi
    fi
    ;;

  gfortran )
    export           MPI_ROOT=""
    export      ESMF_COMPILER=gfortran
    if [ -n "${USE_DEBUG:+1}" ]; then
      export        ESMF_BOPT=g
    else
      export        ESMF_BOPT=O
    fi
    export           ESMF_ABI=64
    export          ESMF_COMM=${which_MPI}
    export          ESMF_SITE=default

    export      ARPACK_LIBDIR=/opt/gfortransoft/serial/ARPACK
    if [ -n "${USE_MPI:+1}" ]; then
      if [ "${which_MPI}" = "mpich2" ]; then
        export       MPI_ROOT=/opt/gfortransoft/mpich2
      elif [ "${which_MPI}" = "openmpi" ]; then
        export       MPI_ROOT=/opt/gfortransoft/openmpi
      fi
      export       MCT_INCDIR=${MPI_ROOT}/mct/include
      export       MCT_LIBDIR=${MPI_ROOT}/mct/lib
      export   PARPACK_LIBDIR=${MPI_ROOT}/PARPACK
    fi

    if [ -n "${USE_NETCDF4:+1}" ]; then
      if [ -n "${USE_PARALLEL_IO:+1}" ] && [ -n "${USE_MPI:+1}" ]; then
        export       ESMF_DIR=${MPI_ROOT}/esmf_nc4
        export         NETCDF=${MPI_ROOT}/netcdf4
        export      NF_CONFIG=${NETCDF}/bin/nf-config
        export  NETCDF_INCDIR=${NETCDF}/include
        export        NETCDF4=1
      else
        export       ESMF_DIR=${MPI_ROOT}/esmf_nc4
        export         NETCDF=/opt/gfortransoft/serial/netcdf4
        export      NF_CONFIG=${NETCDF}/bin/nf-config
        export  NETCDF_INCDIR=${NETCDF}/include
        export        NETCDF4=1
      fi
    else
      export         ESMF_DIR=${MPI_ROOT}/esmf_nc3
      export           NETCDF=/opt/gfortransoft/serial/netcdf3
      export    NETCDF_INCDIR=${NETCDF}/include
      export    NETCDF_LIBDIR=${NETCDF}/lib
      export   NETCDF_classic=1
    fi

    if [ -n "${USE_HDF5:+1}" ]; then
      if [ -n "${USE_PARALLEL_IO:+1}" ] && [ -n "${USE_MPI:+1}" ]; then
        export           HDF5=${MPI_ROOT}/hdf5
        export    HDF5_LIBDIR=${HDF5}/lib
        export    HDF5_INCDIR=${HDF5}/include
      else
        export           HDF5=/opt/gfortransoft/serial/hdf5
        export    HDF5_LIBDIR=${HDF5}/lib
        export    HDF5_INCDIR=${HDF5}/include
      fi
    fi
    ;;

esac
