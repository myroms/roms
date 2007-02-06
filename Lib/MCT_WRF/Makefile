#
CPP=/lib/cpp
CPPFLAGS=-traditional -I. -P
#
F90FREE=mpif90 -Mfree
F90FIXED=mpif90 -Mnofree
FFLAGS=-fast -IMCT/mct -IMCT/mpeu -I$(NETCDF_ROOT)/include

INCLUDES= \
	arrays.xpnd.code \
	pipe_recv.code \
	recv_data.code \
	get_key_value.code \
	pipe_send.code \
	send_data.code \
	macros.h \
	wrf_io_flags.h \

OBJS = \
	mct_wrf_io_utilities.o \
	mct_wrf_io_interpolator.o \
	mct_wrf_io_pipe.o \
	mct_wrf_io_module.o \
	mct_wrf_io.o \

LIBRARY=libwrf_mct_io_api.a

$(LIBRARY): $(OBJS)
	ar r $(LIBRARY) $(OBJS)

$(OBJS) : %.o: %.F
	$(CPP) $(CPPFLAGS) $*.F > $*_cpp.f
	$(F90FREE) $(FFLAGS) -c $*_cpp.f
	@ mv $*_cpp.o $@

SRCS = ${OBJS:.o=.F}

tar :
	tar cvf mct_io_api.tar $(INCLUDES) $(SRCS)

clean:
	@ rm -f *.o
