local_sub  := Nonlinear

local_lib  := libNLM.a
local_src  := $(wildcard $(local_sub)/*.F)
path_srcs  += $(local_src)

local_src  := $(patsubst $(local_sub)/%.F,%.F,$(local_src))
local_objs := $(subst .F,.o,$(local_src))

libraries += $(local_lib)
sources   += $(local_src)

$(local_lib): $(local_objs)
	$(AR) $(ARFLAGS) $@ $^
