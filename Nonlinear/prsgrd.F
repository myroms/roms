#include "cppdefs.h"
      MODULE prsgrd_mod
#ifdef SOLVE3D
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine computes the baroclinic hydrostatic pressure gradient  !
!  term.                                                               !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC prsgrd

      CONTAINS

# if defined PJ_GRADPQ4
#  include "prsgrd44.h"
# elif defined PJ_GRADQ2
#  include "prsgrd42.h"
# elif defined PJ_GRADP
#  include "prsgrd40.h"
# elif defined DJ_GRADPS
#  include "prsgrd32.h"
# else
#  include "prsgrd31.h"
# endif

#endif

      END MODULE prsgrd_mod
