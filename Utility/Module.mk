local_sub  := Utility

local_lib  := libUTIL.a
local_src  := $(wildcard $(local_sub)/*.F)
path_srcs  += $(local_src)

local_src  := $(patsubst $(local_sub)/%.F,%.F,$(local_src))
local_objs := $(subst .F,.o,$(local_src))

libraries += $(local_lib)
sources   += $(local_src)

$(local_lib): $(local_objs)
	$(AR) $(ARFLAGS) $@ $^
