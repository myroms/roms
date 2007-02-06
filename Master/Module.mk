# svn $Id$
#::::::::::::::::::::::::::::::::::::::::::::::::::::: Hernan G. Arango :::
# Copyright (c) 2002-2007 The ROMS/TOMS Group             Kate Hedstrom :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

local_sub  := Master

local_src  := $(wildcard $(local_sub)/*.F)
path_srcs  += $(local_src)

local_src  := $(patsubst $(local_sub)/%.F,%.F,$(local_src))
local_objs := $(subst .F,.o,$(local_src))

sources    += $(local_src)

ifeq ($(OS)-$(strip $(FORT)),CYGWIN-df)
$(BIN):	$(libraries) $(local_objs)
	$(LD) $(FFLAGS) $(local_objs) $(libraries) $(LIBS_WIN32) /exe:$(BIN_WIN32) /link $(LDFLAGS)
else
$(BIN):	$(libraries) $(local_objs)
	$(LD) $(FFLAGS) $(LDFLAGS) $(local_objs) -o $@ $(libraries) $(LIBS)
endif
