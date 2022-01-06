#include "cppdefs.h"
#if defined NONLINEAR && defined FLOATS && defined FLOAT_BIOLOGY
!
!git $Id$
!svn $Id: biology_floats.F 1099 2022-01-06 21:01:01Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2022 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  These routine computes the behavior for biological floats.          !
!                                                                      !
!=======================================================================
!
# ifdef FLOAT_OYSTER
#  include <oyster_floats.h>
# endif
#else
      MODULE biology_floats_mod
      END MODULE biology_floats_mod
#endif
