# git $Id$
# svn $Id: CMakeLists.txt 1210 2024-01-03 22:03:03Z arango $
#:::::::::::::::::::::::::::::::::::::::::::::::::::::: David Robertson :::
# Copyright (c) 2002-2024 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.md                                                 :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
# Source code list for sub-directory "ROMS/Nonlinear"

list( APPEND _files
      ROMS/Nonlinear/BBL/bbl.F
      ROMS/Nonlinear/BBL/bbl_output.F
)

set ( ROMS_Nonlinear_BBL_files
      ${_files}
      PARENT_SCOPE
)
