#!/bin/csh -f
#
# git $Id$
# svn $Id: wrf_links.csh 1202 2023-10-24 15:36:07Z arango $
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Copyright (c) 2002-2023 The ROMS/TOMS Group                           :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::: Hernan G. Arango :::
#                                                                       :::
# WRF em_real symlinks CSH Script: WRF Versions 4.1 and up              :::
#                                                                       :::
# Script called from "build_wrf.csh" if ${WRF_CASE} is set to           :::
# "em_real". This script creates the data links for running the         :::
# 'em_real' executable.                                                 :::
#                                                                       :::
# Usage:                                                                :::
#                                                                       :::
#    This script is called from "build_wrf.csh".                        :::
#                                                                       :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Initialize.

# ${MY_PROJECT_DIR}, ${WRF_ROOT_DIR}, ${USE_REAL_DOUBLE}, ${WRF_BUILD_DIR}
# are environment variables set in "build_wrf.csh".

set separator = `perl -e "print ':' x 100;"`

echo ""
echo "${separator}"
echo "Creating WRF data links:  Case em_real"
echo "${separator}"
echo ""

cd ${MY_PROJECT_DIR}

# If WRF source code is in the project directory (i.e. the -b options
# was used), not setting 'maxdepth' will result in symlinks in the
# source code being removed.

find ./ -maxdepth 1 -type l -exec /bin/rm -f {} \;

foreach file ( ${WRF_ROOT_DIR}/run/* )
  if ( $file =~ *namelist.input* || $file =~ *README* ) then
    continue
  endif
  if ( $file =~ *ETAMPNEW* || $file =~ *RRTM* ) then
    if ( $?USE_REAL_DOUBLE ) then
      if ( $file =~ *_DBL ) then
        set fout = `basename $file _DBL`
        ln -sfv $file $fout
      endif
    else
      if ( $file !~ *_DBL ) then
        ln -sfv $file .
      endif
    endif
  else
    ln -sfv $file .
  endif
end

# Remove all symlinks inside the WRF_BUILD_DIR

find ${WRF_BUILD_DIR} -type l -exec /bin/rm -f {} \;

# Set directory structure needed to compile WPS.

echo ""
echo "${separator}"
echo "Creating WPS directory layout."
echo "${separator}"
echo ""

/bin/mkdir -vp ${WRF_BUILD_DIR}/external/io_grib1
/bin/mkdir -vp ${WRF_BUILD_DIR}/external/io_grib_share
/bin/mkdir -vp ${WRF_BUILD_DIR}/external/io_int
/bin/mkdir -vp ${WRF_BUILD_DIR}/external/io_netcdf
/bin/mkdir -vp ${WRF_BUILD_DIR}/frame

# We do not copy files out of the inc directory so we create a symlink
# do that directory.

echo "cd ${WRF_BUILD_DIR}"
cd ${WRF_BUILD_DIR}
ln -sfv ${WRF_ROOT_DIR}/inc .

# io_grib1

echo "cd ${WRF_BUILD_DIR}/external/io_grib1"
cd ${WRF_BUILD_DIR}/external/io_grib1
ln -sfv ../../libio_grib1.a .

# io_grib_share

echo "cd ${WRF_BUILD_DIR}/external/io_grib_share"
cd ${WRF_BUILD_DIR}/external/io_grib_share
ln -sfv ../../libio_grib_share.a .
ln -sfv ../../wrf_io_flags.h     .
ln -sfv ../../wrf_status_codes.h .

# io_int

echo "cd ${WRF_BUILD_DIR}/external/io_int"
cd ${WRF_BUILD_DIR}/external/io_int
ln -sfv ../../libwrfio_int.a                  .
ln -sfv ../../module_internal_header_util.mod .

# io_netcdf

echo "cd ${WRF_BUILD_DIR}/external/io_netcdf"
cd ${WRF_BUILD_DIR}/external/io_netcdf
ln -sfv ../../libwrfio_nf.a .

# frame

echo "cd ${WRF_BUILD_DIR}/frame"
cd ${WRF_BUILD_DIR}/frame
ln -sfv ../pack_utils.o                  .
ln -sfv ../module_machine.o              .
ln -sfv ../module_internal_header_util.o .
ln -sfv ../module_driver_constants.o     .
