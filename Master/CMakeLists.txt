# git $Id$
# svn $Id: CMakeLists.txt 1099 2022-01-06 21:01:01Z arango $
#:::::::::::::::::::::::::::::::::::::::::::::::::::::: David Robertson :::
# Copyright (c) 2002-2022 The ROMS/TOMS Group                           :::
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
