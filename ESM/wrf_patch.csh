#!/bin/csh -f
#
# git $Id$
# svn $Id: wrf_patch.csh 1071 2021-06-14 04:22:26Z arango $
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Copyright (c) 2002-2021 The ROMS/TOMS Group                           :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::: Hernan G. Arango :::
#                                                                       :::
# WRF Patching CSH Script: WRF Versions 4.1 and up                      :::
#                                                                       :::
# Script called from "build_wrf.csh" to check whether WRF source code   :::
# has been patched for NetCDF4 library depencies, added configure       :::
# options, creating clean .f90 files for debugging, rename modules to   :::
# WRF_ESMF_*, and correct optional argment from defaultCalendar to      :::
# defaultCalKind in ESMF_Initialize call.                               :::
#                                                                       :::
# Usage:                                                                :::
#                                                                       :::
#    This script is called from "build_wrf.csh".                        :::
#                                                                       :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Initialize.

# ${WRF_ROOT_DIR} and ${ROMS_SRC_DIR} are environment variables set in
# "build_wrf.csh".

set separator = `perl -e "print ':' x 100;"`
set CHECK_STRING = 'WRF-ROMS ESMF-NUOPC Coupling'

echo ""
echo "${separator}"
echo "If applicable, replacing several WRF files for WRF ESMF/NUOPC Coupling"
echo "${separator}"
echo ""

# Reworking linking NetCDF4 library dependencies

if (! `grep -c "${CHECK_STRING}" ${WRF_ROOT_DIR}/configure`) then
  mv -v ${WRF_ROOT_DIR}/configure ${WRF_ROOT_DIR}/configure.orig
  cp -fv ${ROMS_SRC_DIR}/ESM/wrf_configure ${WRF_ROOT_DIR}/configure
else
  echo "   No need to replace: ${WRF_ROOT_DIR}/configure"
endif

if (! `grep -c "${CHECK_STRING}" ${WRF_ROOT_DIR}/Makefile`) then
  mv -v  ${WRF_ROOT_DIR}/Makefile ${WRF_ROOT_DIR}/Makefile.orig
  cp -fv ${ROMS_SRC_DIR}/ESM/wrf_Makefile  ${WRF_ROOT_DIR}/Makefile
else
  echo "   No need to replace: ${WRF_ROOT_DIR}/Makefile"
endif

if (! `grep -c "${CHECK_STRING}" ${WRF_ROOT_DIR}/arch/postamble`) then
  mv -v  ${WRF_ROOT_DIR}/arch/postamble ${WRF_ROOT_DIR}/arch/postamble.orig
  cp -fv ${ROMS_SRC_DIR}/ESM/wrf_postamble ${WRF_ROOT_DIR}/arch/postamble
else
  echo "   No need to replace: ${WRF_ROOT_DIR}/arch/postamble"
endif

# Adding MacOS Intel/GNU with OpenMPI

if (! `grep -c "${CHECK_STRING}" ${WRF_ROOT_DIR}/arch/configure.defaults`) then
  mv -v  ${WRF_ROOT_DIR}/arch/configure.defaults ${WRF_ROOT_DIR}/arch/configure.defaults.orig
  cp -fv ${ROMS_SRC_DIR}/ESM/wrf_configure.defaults ${WRF_ROOT_DIR}/arch/configure.defaults
else
  echo "   No need to replace: ${WRF_ROOT_DIR}/arch/configure.defaults"
endif

# Create clean .f90 files for debugging and rename modules to WRF_ESMF_*

if (! `grep -c "${CHECK_STRING}" ${WRF_ROOT_DIR}/external/esmf_time_f90/Makefile`) then
  mv -v  ${WRF_ROOT_DIR}/external/esmf_time_f90/Makefile ${WRF_ROOT_DIR}/external/esmf_time_f90/Makefile.orig
  cp -fv ${ROMS_SRC_DIR}/ESM/wrf_Makefile.esmf ${WRF_ROOT_DIR}/external/esmf_time_f90/Makefile
else
  echo "   No need to replace: ${WRF_ROOT_DIR}/external/esmf_time_f90/Makefile"
endif

# Correcting optional argument from defaultCalendar to defaultCalKind in
# ESMF_Initialize call

if (! `grep -c "${CHECK_STRING}" ${WRF_ROOT_DIR}/external/esmf_time_f90/Test1.F90`) then
  mv -v  ${WRF_ROOT_DIR}/external/esmf_time_f90/Test1.F90 ${WRF_ROOT_DIR}/external/esmf_time_f90/Test1.F90.orig
  cp -fv ${ROMS_SRC_DIR}/ESM/wrf_Test1.F90 ${WRF_ROOT_DIR}/external/esmf_time_f90/Test1.F90
else
  echo "   No need to replace: ${WRF_ROOT_DIR}/external/esmf_time_f90/Test1.F90"
endif
