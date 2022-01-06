# git $Id$
# svn $Id: CMakeLists.txt 1099 2022-01-06 21:01:01Z arango $
#:::::::::::::::::::::::::::::::::::::::::::::::::::::: David Robertson :::
# Copyright (c) 2002-2022 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
# Source code list for sub-directory "ROMS/Tangent/Biology"

list( APPEND _files
      ROMS/Tangent/Biology/tl_biology.F
)

set( ROMS_Tangent_Biology_files
     ${_files}
     PARENT_SCOPE
)
