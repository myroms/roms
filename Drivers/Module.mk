local_sub  := Drivers

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
