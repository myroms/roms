#!/bin/make

MODULE		= utmct

makefiles	= configure	Makefile	Makefile.rules	\
		  Makefile.conf

SRCS_C		= 

SRCS_F		= 

SRCS_f		= 

SRCS_F90	= mph.F90 master.F90 convertgauss.F90 convertPOPT.F90 \
	          cpl.F90 ccm.F90 pop.F90 mpe_logging_init.F90 \
                  ReadSparseMatrixAsc.F90

SRCS_f90	= 

SRCS_ALL	= $(SRCS_f) $(SRCS_F) $(SRCS_f90) $(SRCS_F90) $(SRCS_C)

OBJS_ALL	= $(SRCS_F:.F=.o)	$(SRCS_F90:.F90=.o)	\
		  $(SRCS_f90:.f90=.o)	$(SRCS_f:.f=.o) \
		  $(SRCS_C:.c=.o)

ARCH_TARGETS	= $(makefiles) $(SRCS_ALL) ut_mp_sym_Cxpy.F90 \
		  ut_mp_rec_Cxpy.F90

include Makefile.rules
include Makefile.conf

# targets

pt:
	protex $(SRCS_ALL)

utmct:	$(OBJS_ALL) 
	$(LD) -o $@ $(OBJS_ALL) $(LDFLAGS) $(_L) $(_l)

all:	utmct



# dependencies:





