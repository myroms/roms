#include "cppdefs.h"
      MODULE mod_stepping
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2011 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This MODULE contains time stepping indices.                         !
!                                                                      !
#ifdef BOUNDARY
!  Lbinp     Open boundary adjustment input  fields index.             !
!  Lbout     Open boundary adjustment output fields index.             !
#endif
#if defined ADJUST_BOUNDARY || defined ADJUST_STFLUX || \
    defined ADJUST_WSTRESS
!  Lfinp     Surface forcing adjustment input  fields index.           !
!  Lfout     Surface forcing adjustment output fields index.           !
#endif
!  Lnew      New descent algorithm state solution index.               !
!  Lold      Previous descent algorithm state solution index.          !
!                                                                      !
!  knew      Barotropic (fast) time-step index corresponding to the    !
!              newest values for 2D primitive equation variables.      !
!  krhs      Barotropic (fast) time-step index used to compute the     !
!              right-hand-terms of 2D primitive equation variables.    !
!  kstp      Barotropic (fast) time-step index to which the current    !
!              changes are added to compute new 2D primitive equation  !
!              variables.                                              !
!                                                                      !
!  nfm3      Float index for time level "n-3".                         !
!  nfm2      Float index for time level "n-2".                         !
!  nfm1      Float index for time level "n-1".                         !
!  nf        Float index for time level "n".                           !
!  nfp1      Float index for time level "n+1".                         !
!                                                                      !
!  nnew      Baroclinic (slow) time-step index corresponding to the    !
!              newest values for 3D primitive equation variables.      !
!  nrhs      Baroclinic (slow) time-step index used to compute the     !
!              right-hand-terms of 3D primitive equation variables.    !
!  nstp      Baroclinic (slow) time-step index to which the current    !
!              changes are added to compute new 3D primitive equation  !
!              variables.                                              !
#if defined SSH_TIDES || defined UV_TIDES
!                                                                      !
!  NTC       Number of tidal components to consider.                   !
#endif
!                                                                      !
!=======================================================================
!
        USE mod_param
!
        implicit none
!
        integer, allocatable :: knew(:)
        integer, allocatable :: krhs(:)
        integer, allocatable :: kstp(:)

        integer, allocatable :: nnew(:)
        integer, allocatable :: nrhs(:)
        integer, allocatable :: nstp(:)

#ifdef FLOATS
        integer, allocatable :: nf(:)
        integer, allocatable :: nfp1(:)
        integer, allocatable :: nfm3(:)
        integer, allocatable :: nfm2(:)
        integer, allocatable :: nfm1(:)
#endif

#ifdef ADJUST_BOUNDARY
        integer, allocatable :: Lbinp(:)
        integer, allocatable :: Lbout(:)
#endif

#if defined ADJUST_BOUNDARY || defined ADJUST_STFLUX || \
    defined ADJUST_WSTRESS
        integer, allocatable :: Lfinp(:)
        integer, allocatable :: Lfout(:)
#endif

        integer, allocatable :: Lnew(:)
        integer, allocatable :: Lold(:)

#if defined SSH_TIDES || defined UV_TIDES
        integer, allocatable :: NTC(:)
#endif
!
      CONTAINS
!
      SUBROUTINE allocate_stepping
!
!=======================================================================
!                                                                      !
!  This routine allocates several variables in the module that depend  !
!  on the number of nested grids.                                      !
!                                                                      !
!=======================================================================
!
!  Local variable declarations.
!
      integer :: ng
!
!-----------------------------------------------------------------------
!  Allocate variables.
!-----------------------------------------------------------------------
!
      allocate ( knew(Ngrids) )
      allocate ( krhs(Ngrids) )
      allocate ( kstp(Ngrids) )

      allocate ( nnew(Ngrids) )
      allocate ( nrhs(Ngrids) )
      allocate ( nstp(Ngrids) )

#ifdef FLOATS
      allocate ( nf(Ngrids) )
      allocate ( nfp1(Ngrids) )
      allocate ( nfm3(Ngrids) )
      allocate ( nfm2(Ngrids) )
      allocate ( nfm1(Ngrids) )
#endif

#ifdef ADJUST_BOUNDARY
      allocate ( Lbinp(Ngrids) )
      allocate ( Lbout(Ngrids) )
#endif

#if defined ADJUST_BOUNDARY || defined ADJUST_STFLUX || \
    defined ADJUST_WSTRESS
      allocate ( Lfinp(Ngrids) )
      allocate ( Lfout(Ngrids) )
#endif

      allocate ( Lnew(Ngrids) )
      allocate ( Lold(Ngrids) )

#if defined SSH_TIDES || defined UV_TIDES
      allocate ( NTC(Ngrids) )
#endif
!
!-----------------------------------------------------------------------
!  Initialize time-stepping indices.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        knew(ng)=1
        krhs(ng)=1
        kstp(ng)=1

        nnew(ng)=1
        nrhs(ng)=1
        nstp(ng)=1

#ifdef FLOATS
        nf  (ng)=0
        nfp1(ng)=1
        nfm3(ng)=2
        nfm2(ng)=3
        nfm1(ng)=4
#endif

#ifdef ADJUST_BOUNDARY
        Lbinp(ng)=1
        Lbout(ng)=1
#endif

#if defined ADJUST_BOUNDARY || defined ADJUST_STFLUX || \
    defined ADJUST_WSTRESS
        Lfinp(ng)=1
        Lfout(ng)=1
#endif

        Lnew(ng)=1
        Lold(ng)=1
      END DO

      RETURN
      END SUBROUTINE allocate_stepping

      END MODULE mod_stepping

