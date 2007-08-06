#include "cppdefs.h"
      SUBROUTINE get_bounds (ng, tile, gtype, Nghost, Itile, Jtile,     &
     &                       LBi, UBi, LBj, UBj)
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2007 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine compute grid bounds in the I- and J-directions.        !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     tile       Domain partition.                                     !
!     gtype      C-grid type. If zero, compute array allocation bounds.!
!                  Otherwise, compute bounds for IO processing.        !
!     Nghost     Number of ghost-points in the halo region:            !
!                  Nghost = 0,  compute non-overlaping bounds.         !
!                  Nghost > 0,  compute overlaping bounds.             !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Itile      Tile coordinate in the I-direction, use only in       !
!                  distributed-memory applications.                    !
!     Jtile      Tile coordinate in the J-direction, use only in       !
!                  distributed-memory applications.                    !
!     LBi        I-dimension Lower bound.                              !
!     UBi        I-dimension Upper bound.                              !
!     LBj        J-dimension Lower bound.                              !
!     UBj        J-dimension Upper bound.                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_ncparam
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, gtype, Nghost
      integer, intent(out) :: Itile, Jtile, LBi, UBi, LBj, UBj
#ifdef DISTRIBUTE
!
!  Local variable declarations.
!
      integer :: Iend, Istr, Jend, Jstr
      integer :: IstrR, IstrU, IendR, JstrR, JstrV, JendR
      integer :: Imin, Imax, Jmin, Jmax
      integer :: MyType
!
!-----------------------------------------------------------------------
!  Set array bounds in the I- and J-direction for distributed-memory
!  configurations.
!-----------------------------------------------------------------------
!
!  Set first and last grid-points according to staggered C-grid
!  classification.  If gtype = 0, it returns the values needed for
!  array allocation. Otherwise, it returns the values needed for IO
!  processing.
!
      MyType=ABS(gtype)
      IF (MyType.eq.0) THEN
        Imin=LOWER_BOUND_I
        Imax=UPPER_BOUND_I
        Jmin=LOWER_BOUND_J
        Jmax=UPPER_BOUND_J
      ELSE
        Imin=0
        Imax=Lm(ng)+1
        Jmin=0
        Jmax=Mm(ng)+1
        IF ((MyType.eq.p2dvar).or.(MyType.eq.u2dvar).or.                &
     &      (MyType.eq.p3dvar).or.(MyType.eq.u3dvar)) Imin=1
        IF ((MyType.eq.p2dvar).or.(MyType.eq.v2dvar).or.                &
     &      (MyType.eq.p3dvar).or.(MyType.eq.v3dvar)) Jmin=1
      END IF
!
!  Set physical, overlaping (Nghost>0) or non-overlaping (Nghost=0) 
!  grid bounds according to tile rank.
!
      CALL get_tile (ng, tile, Itile, Jtile, Istr, Iend, Jstr, Jend,    &
     &               IstrR, IstrU, IendR, JstrR, JstrV, JendR)

!
      IF (Itile.eq.0) THEN
        LBi=Imin
      ELSE
        LBi=Istr-Nghost
      END IF
      IF (Itile.eq.(NtileI(ng)-1)) THEN
        UBi=Imax
      ELSE
        UBi=Iend+Nghost
      END IF
      IF (Jtile.eq.0) THEN
        LBj=Jmin
      ELSE
        LBj=Jstr-Nghost
      END IF
      IF (Jtile.eq.(NtileJ(ng)-1)) THEN
        UBj=Jmax
      ELSE
        UBj=Jend+Nghost
      END IF
#else
!
!-----------------------------------------------------------------------
!  Set array allocation bounds in the I- and J-direction for serial and
!  shared-memory configurations.
!-----------------------------------------------------------------------
!
      Itile=-1
      Jtile=-1
      LBi=LOWER_BOUND_I
      UBi=UPPER_BOUND_I
      LBj=LOWER_BOUND_J
      UBj=UPPER_BOUND_J
#endif
      RETURN
      END SUBROUTINE get_bounds

      SUBROUTINE get_tile (ng, tile, Itile, Jtile,                      &
     &                     Istr, Iend, Jstr, Jend,                      &
     &                     IstrR, IstrU, IendR,                         &
     &                     JstrR, JstrV, JendR)
!
!=======================================================================
!                                                                      !
!  This routine computes the starting and ending horizontal indices    !
!  for each sub-domain partition or tile.                              !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number (integer).                         !
!     tile       Sub-domain partition.                                 !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Itile      I-tile coordinate (a value from 0 to NtileI(ng)).     !
!     Jtile      J-tile coordinate (a value from 0 to NtileJ(ng)).     !
!     Istr       Starting tile index in the I-direction.               !
!     Iend       Ending   tile index in the I-direction.               !
!     Jstr       Starting tile index in the J-direction.               !
!     Jend       Ending   tile index in the J-direction.               !
!                                                                      !
!     IstrR      Starting tile index in the I-direction (RHO-points).  !
!     IstrU      Starting tile index in the I-direction (U-points).    !
!     IendR      Ending   tile index in the I-direction (RHO_points).  !
!                                                                      !
!     JstrR      Starting tile index in the J-direction (RHO-points).  !
!     JstrV      Starting tile index in the J-direction (V-points).    !
!     JendR      Ending   tile index in the J-direction (RHO_points).  !
!                                                                      !
!======================================================================!
!
      USE mod_param
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(out) :: Itile, Jtile
      integer, intent(out) :: Iend, Istr, Jend, Jstr
      integer, intent(out) :: IstrR, IstrU, IendR
      integer, intent(out) :: JstrR, JstrV, JendR
!
!  Local variable declarations.
!
      integer :: ChunkSizeI, ChunkSizeJ, MarginI, MarginJ
!
!-----------------------------------------------------------------------
!  Set physical non-overlaping grid bounds according to tile rank.
!-----------------------------------------------------------------------
!
      ChunkSizeI=(Lm(ng)+NtileI(ng)-1)/NtileI(ng)
      ChunkSizeJ=(Mm(ng)+NtileJ(ng)-1)/NtileJ(ng)
      MarginI=(NtileI(ng)*ChunkSizeI-Lm(ng))/2
      MarginJ=(NtileJ(ng)*ChunkSizeJ-Mm(ng))/2
      Jtile=tile/NtileI(ng)
      Itile=tile-Jtile*NtileI(ng)
!
!  Grid bounds in the I-direction.
!
      Istr=1+Itile*ChunkSizeI-MarginI
      Iend=Istr+ChunkSizeI-1
      Istr=MAX(Istr,1)
      Iend=MIN(Iend,Lm(ng))
!
!  Grid bounds in the J-direction.
!
      Jstr=1+Jtile*ChunkSizeJ-MarginJ
      Jend=Jstr+ChunkSizeJ-1
      Jstr=MAX(Jstr,1)
      Jend=MIN(Jend,Mm(ng))

#include "set_bounds.h"

      RETURN
      END SUBROUTINE get_tile
