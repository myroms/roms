#include "cppdefs.h"
      MODULE mod_stepping
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This MODULE contains time stepping indices.                         !
!                                                                      !
!  Lnew      New descent algorithm state solution index.               !
!  Lold      Previous descent algorithm state solution index.          !
!  knew      Fast time index corresponding to the newest time step     !
!              for variables associated with 2D primitive equations    !
!  krhs      Fast time index corresponding to the time step 2D         !
!              primitive variables used to calculate the right-hand    !
!              terms.                                                  !
!  kstp      Fast time index corresponding to the time step 2D         !
!              primitive variables to which the current changes        !
!              are added to compute the new time step values.          !
!                                                                      !
!  nfm3      Float index for time level "n-3".                         !
!  nfm2      Float index for time level "n-2".                         !
!  nfm1      Float index for time level "n-1".                         !
!  nf        Float index for time level "n".                           !
!  nfp1      Float index for time level "n+1".                         !
!                                                                      !
!  nnew                                                                !
!  nrhs                                                                !
!  nstp                                                                !
!  NTC       Number of tidal components to consider.                   !
!                                                                      !
!=======================================================================
!
        USE mod_param

        implicit none

        integer, private :: ig

        integer, dimension(Ngrids) :: knew = (/ (1, ig=1,Ngrids) /)
        integer, dimension(Ngrids) :: krhs = (/ (1, ig=1,Ngrids) /)
        integer, dimension(Ngrids) :: kstp = (/ (1, ig=1,Ngrids) /)

        integer, dimension(Ngrids) :: nnew = (/ (1, ig=1,Ngrids) /)
        integer, dimension(Ngrids) :: nrhs = (/ (1, ig=1,Ngrids) /)
        integer, dimension(Ngrids) :: nstp = (/ (1, ig=1,Ngrids) /)

#ifdef FLOATS
        integer, dimension(Ngrids) :: nf   = (/ (0, ig=1,Ngrids) /)
        integer, dimension(Ngrids) :: nfp1 = (/ (1, ig=1,Ngrids) /)
        integer, dimension(Ngrids) :: nfm3 = (/ (2, ig=1,Ngrids) /)
        integer, dimension(Ngrids) :: nfm2 = (/ (3, ig=1,Ngrids) /)
        integer, dimension(Ngrids) :: nfm1 = (/ (4, ig=1,Ngrids) /)
#endif

#if defined SSH_TIDES || defined UV_TIDES
        integer, dimension(Ngrids) :: NTC
#endif
        integer, dimension(Ngrids) :: Lnew = (/ (1, ig=1,Ngrids) /)
        integer, dimension(Ngrids) :: Lold = (/ (1, ig=1,Ngrids) /)

      END MODULE mod_stepping

