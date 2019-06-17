#include "cppdefs.h"
      MODULE biology_floats_mod
#if defined NONLINEAR && defined FLOATS && defined FLOAT_BIOLOGY
!
!git $Id: 560fb376ff8a4576170ebcd4b459de6bcce908f6 $
!svn $Id: biology_floats.F 937 2019-01-28 06:13:04Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2019 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  These routine computes the behavior for biological floats.          !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: biology_floats

      CONTAINS

# ifdef FLOAT_OYSTER
#  include <oyster_floats.h>
# endif

#endif

      END MODULE biology_floats_mod
