#include "cppdefs.h"
#if defined ADJOINT && defined SOLVE3D
!!
!!git $Id$
!!svn $Id: ad_prsgrd.F 1210 2024-01-03 22:03:03Z arango $
!!================================================= Hernan G. Arango ===
!! Copyright (c) 2002-2024 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.md                                               !
!!======================================================================
!!                                                                     !
!! This routine computes the adjoint baroclinic hydrostatic pressure   !
!! gradient term.                                                      !
!!                                                                     !
!!======================================================================
!!
# if defined PJ_GRADPQ4
!!#  include "ad_prsgrd44.h"
# elif defined PJ_GRADPQ2
!!#  include "ad_prsgrd42.h"
# elif defined PJ_GRADP
#  include "ad_prsgrd40.h"
# elif defined DJ_GRADPS
#  include "ad_prsgrd32.h"
# else
#  include "ad_prsgrd31.h"
# endif
#else
      MODULE ad_prsgrd_mod
      END MODULE ad_prsgrd_mod
#endif
