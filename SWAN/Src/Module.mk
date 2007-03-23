local_sub  := SWAN/Src

local_lib  := libSWAN.a
local_src  := $(wildcard $(local_sub)/*.F)
path_srcs  += $(local_src)

local_src  := $(patsubst $(local_sub)/%.F,%.F,$(local_src))
local_objs := $(subst .F,.o,$(local_src))

libraries += $(local_lib)
sources   += $(local_src)

swch = -unix -f95 -mpi

#	@perl SWAN/switch.pl $(swch) SWAN/*.ftn
$(local_lib): $(local_objs)
	$(AR) $(ARFLAGS) $@ $^
	$(RANLIB) $@


