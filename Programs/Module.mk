local_sub  := Programs

local_src  := $(wildcard $(local_sub)/*.F)
path_srcs  += $(local_src)

local_src  := $(patsubst $(local_sub)/%.F,%.F,$(local_src))
local_objs := $(subst .F,$(OBJ_EXT),$(local_src))

sources   += $(local_src)
