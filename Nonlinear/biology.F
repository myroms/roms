#include "cppdefs.h"
      MODULE biology_mod
#ifdef BIOLOGY
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine computes the source and sink terms for selected        !
!   biology model.                                                     !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC biology_tile

      CONTAINS

# ifdef BIO_FASHAM
#  include "fasham.h"
# endif

#endif
      END MODULE biology_mod
