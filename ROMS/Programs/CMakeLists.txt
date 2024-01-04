# git $Id$
# svn $Id: CMakeLists.txt 1195 2023-08-31 01:52:55Z arango $
#:::::::::::::::::::::::::::::::::::::::::::::::::::::: David Robertson :::
# Copyright (c) 2002-2024 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.md                                                 :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
# Source code list for sub-directory "ROMS/Programs"

list( APPEND _files
      ROMS/Programs/types.F
      ROMS/Programs/yaml_parser_test.F
)

set ( ROMS_Programs_files
      ${_files}
      PARENT_SCOPE
)
