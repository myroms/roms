# git $Id$
#::::::::::::::::::::::::::::::::::::::::::::::::::::: Hernan G. Arango :::
# Copyright (c) 2002-2024 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.md                                                 :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

local_sub  := ROMS/Nonlinear/BBL

local_src  := $(wildcard $(local_sub)/*.F)

sources    += $(local_src)

$(eval $(compile-rules))
