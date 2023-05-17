# git $Id$
# svn $Id: CMakeLists.txt 1166 2023-05-17 20:11:58Z arango $
#:::::::::::::::::::::::::::::::::::::::::::::::::::::: David Robertson :::
# Copyright (c) 2002-2023 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
# Source code list for sub-directory "Master"

list( APPEND _files
      Master/cmeps_roms.F
      Master/coupler.F
      Master/esmf_roms.F
      Master/propagator.F
      Master/roms_kernel.F
)

set( Master_files
     ${_files}
     PARENT_SCOPE
)
