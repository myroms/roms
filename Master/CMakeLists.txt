# git $Id$
# svn $Id: CMakeLists.txt 1202 2023-10-24 15:36:07Z arango $
#:::::::::::::::::::::::::::::::::::::::::::::::::::::: David Robertson :::
# Copyright (c) 2002-2024 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.md                                                 :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
# Source code list for sub-directory "Master". We cannot add "master.F" to
# list because it is needed to define the ROMS executable in the main
# CMakeList.txt.

list( APPEND _files
      Master/coupler.F
      Master/esmf_atm.F
      Master/esmf_data.F
      Master/esmf_esm.F
      Master/esmf_ice.F
      Master/esmf_roms.F
      Master/esmf_wav.F
      Master/mod_esmf_esm.F
      Master/propagator.F
      Master/roms_kernel.F
)

set( Master_files
     ${_files}
     PARENT_SCOPE
)
