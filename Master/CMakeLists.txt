# git $Id$
# svn $Id: CMakeLists.txt 1092 2021-09-29 04:30:19Z arango $
#:::::::::::::::::::::::::::::::::::::::::::::::::::::: David Robertson :::
# Copyright (c) 2002-2021 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
# Source code list for sub-directory "Master"

list( APPEND _files
      Master/coupler.F
      Master/esmf_roms.F
      Master/propagator.F
      Master/roms_kernel.F
)

set( Master_files
     ${_files}
     PARENT_SCOPE
)
