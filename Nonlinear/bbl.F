#include "cppdefs.h"
      MODULE bbl_mod
#ifdef BBL_MODEL
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine computes bottom momentun stress via a bottom boundary  !
!  layer formulation.                                                  !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC bblm

      CONTAINS

# ifdef SG_BBL
#  include "sg_bbl.h"
# endif

#endif

      END MODULE bbl_mod
