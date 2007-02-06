#!/bin/make

MODULE		= twocon

makefiles	= configure	Makefile	Makefile.rules	\
		  Makefile.conf 

SRCS_C		= 

SRCS_F		= 

SRCS_f		= 

SRCS_F90	= twocomponent.concurrent.F90 twocomponent.sequential.F90 \
		  twocomponentUneven.sequential.F90

SRCS_f90	= 

SRCS_ALL	= $(SRCS_f) $(SRCS_F) $(SRCS_f90) $(SRCS_F90) $(SRCS_C)

OBJS_ALL	= $(SRCS_F:.F=.o)	$(SRCS_F90:.F90=.o)	\
		  $(SRCS_f90:.f90=.o)	$(SRCS_f:.f=.o) \
		  $(SRCS_C:.c=.o)

ARCH_TARGETS	= $(makefiles) $(SRCS_ALL) 

MPEU_IDIR       = ../../mpeu
MPEU_LDIR       = ../../mpeu
MCT_IDIR        = ../../mct
MCT_LDIR        = ../../mct

include Makefile.rules
include Makefile.conf

# targets

pt:
	protex $(SRCS_ALL)

twocon:  twocomponent.concurrent.o
	$(LD) -o $@ twocomponent.concurrent.o $(LDFLAGS) $(_L) $(_l)

twoseq:  twocomponent.sequential.o
	$(LD) -o $@ twocomponent.sequential.o $(LDFLAGS) $(_L) $(_l)

twosequn:  twocomponentUneven.sequential.o
	$(LD) -o $@ twocomponentUneven.sequential.o $(LDFLAGS) $(_L) $(_l)

all:	twocon twoseq twosequn



# dependencies:





