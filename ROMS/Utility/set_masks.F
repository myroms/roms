#include "cppdefs.h"
      MODULE set_masks_mod
#ifdef MASKING
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2011 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  These routines set internal Land/Sea masking arrays that are used   !
!  to process fields into output NetCDF files.  The Land grid points   !
!  are replaced by the _FillValue in the output files to  facilitate   !
!  post-processing with generic tools.                                 !
!                                                                      !
# ifdef UV_PSOURCE
!  If point sources, insure that masks at point source locations are   !
!  set to water to avoid masking with _FillValue at those locations.   !
# endif
# ifdef WET_DRY
!
!  If wetting and drying, masks array are time dependent and changed   !
!  at every time-step in routine "wetdry".  Notice that time-average   !
!  masks are needed for the selected time window.
# endif
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC :: set_masks
# if defined WET_DRY     && \
     defined AVERAGES    || \
    (defined AD_AVERAGES && defined ADJOINT) || \
    (defined RP_AVERAGES && defined TL_IOMS) || \
    (defined TL_AVERAGES && defined TANGENT)
      PUBLIC :: set_avg_masks
# endif
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE set_masks (ng, tile, model)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
# ifdef UV_PSOURCE
      USE mod_sources
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, model, 2)
# endif
      CALL set_masks_tile (ng, tile, model,                             &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
# ifdef UV_PSOURCE
     &                     Msrc(ng), Nsrc(ng),                          &
     &                     SOURCES(ng) % Isrc,                          &
     &                     SOURCES(ng) % Jsrc,                          &
     &                     SOURCES(ng) % Dsrc,                          &
# endif
     &                     GRID(ng) % pmask,                            &
     &                     GRID(ng) % rmask,                            &
     &                     GRID(ng) % umask,                            &
     &                     GRID(ng) % vmask,                            &
# if defined AVERAGES    || \
    (defined AD_AVERAGES && defined ADJOINT) || \
    (defined RP_AVERAGES && defined TL_IOMS) || \
    (defined TL_AVERAGES && defined TANGENT)
     &                     GRID(ng) % pmask_avg,                        &
     &                     GRID(ng) % rmask_avg,                        &
     &                     GRID(ng) % umask_avg,                        &
     &                     GRID(ng) % vmask_avg,                        &
# endif
# ifdef DIAGNOSTICS
     &                     GRID(ng) % pmask_dia,                        &
     &                     GRID(ng) % rmask_dia,                        &
     &                     GRID(ng) % umask_dia,                        &
     &                     GRID(ng) % vmask_dia,                        &
# endif
     &                     GRID(ng) % pmask_io,                         &
     &                     GRID(ng) % rmask_io,                         &
     &                     GRID(ng) % umask_io,                         &
     &                     GRID(ng) % vmask_io)
# ifdef PROFILE
      CALL wclock_off (ng, model, 2)
# endif

      RETURN
      END SUBROUTINE set_masks
!
!***********************************************************************
      SUBROUTINE set_masks_tile (ng, tile, model,                       &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
# ifdef UV_PSOURCE
     &                           Msrc, Nsrc,                            &
     &                           Isrc, Jsrc, Dsrc,                      &
# endif
     &                           pmask, rmask,                          &
     &                           umask, vmask,                          &
# if defined AVERAGES    || \
    (defined AD_AVERAGES && defined ADJOINT) || \
    (defined RP_AVERAGES && defined TL_IOMS) || \
    (defined TL_AVERAGES && defined TANGENT)
     &                           pmask_avg, rmask_avg,                  &
     &                           umask_avg, vmask_avg,                  &
# endif
# ifdef DIAGNOSTICS
     &                           pmask_dia, rmask_dia,                  &
     &                           umask_dia, vmask_dia,                  &
# endif
     &                           pmask_io, rmask_io,                    &
     &                           umask_io, vmask_io)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE exchange_2d_mod
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
# ifdef UV_PSOURCE
      integer, intent(in) :: Msrc, Nsrc
# endif
!
# ifdef ASSUMED_SHAPE
#  ifdef UV_PSOURCE
      integer, intent(in) :: Isrc(:)
      integer, intent(in) :: Jsrc(:)

      real(r8), intent(in) :: Dsrc(:)
#  endif
      real(r8), intent(in) :: pmask(LBi:,LBj:)
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  if defined AVERAGES    || \
     (defined AD_AVERAGES && defined ADJOINT) || \
     (defined RP_AVERAGES && defined TL_IOMS) || \
     (defined TL_AVERAGES && defined TANGENT)
      real(r8), intent(inout) :: pmask_avg(LBi:,LBj:)
      real(r8), intent(inout) :: rmask_avg(LBi:,LBj:)
      real(r8), intent(inout) :: umask_avg(LBi:,LBj:)
      real(r8), intent(inout) :: vmask_avg(LBi:,LBj:)
#  endif
#  ifdef DIAGNOSTICS
      real(r8), intent(inout) :: pmask_dia(LBi:,LBj:)
      real(r8), intent(inout) :: rmask_dia(LBi:,LBj:)
      real(r8), intent(inout) :: umask_dia(LBi:,LBj:)
      real(r8), intent(inout) :: vmask_dia(LBi:,LBj:)
#  endif
      real(r8), intent(inout) :: pmask_io(LBi:,LBj:)
      real(r8), intent(inout) :: rmask_io(LBi:,LBj:)
      real(r8), intent(inout) :: umask_io(LBi:,LBj:)
      real(r8), intent(inout) :: vmask_io(LBi:,LBj:)
# else
#  ifdef UV_PSOURCE
      integer, intent(in) :: Isrc(Msrc)
      integer, intent(in) :: Jsrc(Msrc)

      real(r8), intent(in) :: Dsrc(Msrc)
#  endif
      real(r8), intent(in) :: pmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  if defined AVERAGES    || \
     (defined AD_AVERAGES && defined ADJOINT) || \
     (defined RP_AVERAGES && defined TL_IOMS) || \
     (defined TL_AVERAGES && defined TANGENT)
      real(r8), intent(inout) :: pmask_avg(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: rmask_avg(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: umask_avg(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: vmask_avg(LBi:UBi,LBj:UBj)
#  endif
#  ifdef DIAGNOSTICS
      real(r8), intent(inout) :: pmask_dia(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: rmask_dia(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: umask_dia(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: vmask_dia(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(inout) :: pmask_io(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: rmask_io(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: umask_io(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: vmask_io(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: i, j
# ifdef UV_PSOURCE
      integer :: is
# endif

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Initialize internal history files Land/Sea masks with its respective
!  application grid mask.
!-----------------------------------------------------------------------
!
      DO j=Jstr,JendR
        DO i=Istr,IendR
          pmask_io(i,j)=pmask(i,j)
        END DO
      END DO
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          rmask_io(i,j)=rmask(i,j)
        END DO
      END DO
      DO j=JstrR,JendR
        DO i=Istr,IendR
          umask_io(i,j)=umask(i,j)
        END DO
      END DO
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          vmask_io(i,j)=vmask(i,j)
        END DO
      END DO

# ifdef UV_PSOURCE
!
!  Insure that masks at mass point source locations are set to water
!  to avoid masking with _FillValue at those locations.
!
      DO is=1,Nsrc
        i=Isrc(is)
        j=Jsrc(is)
        IF (((IstrR.le.i).and.(i.le.IendR)).and.                        &
     &      ((JstrR.le.j).and.(j.le.JendR))) THEN
          IF (INT(Dsrc(is)).eq.0) THEN
            umask_io(i,j)=1.0_r8
          ELSE
            vmask_io(i,j)=1.0_r8
          END IF
        END IF
      END DO
# endif
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_p2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          pmask_io)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          rmask_io)
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          umask_io)
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          vmask_io)
      END IF

# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 4,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    pmask_io, rmask_io, umask_io, vmask_io)
# endif

# if defined AVERAGES    || \
    (defined AD_AVERAGES && defined ADJOINT) || \
    (defined RP_AVERAGES && defined TL_IOMS) || \
    (defined TL_AVERAGES && defined TANGENT)
!
!-----------------------------------------------------------------------
!  Initialize average file Land/Sea masks for time-averaged fields.
!-----------------------------------------------------------------------
!
      DO j=Jstr,JendR
        DO i=Istr,IendR
#  ifdef WET_DRY
          pmask_avg(i,j)=0.0_r8
#  else
          pmask_avg(i,j)=pmask_io(i,j)
#  endif
        END DO
      END DO

      DO j=JstrR,JendR
        DO i=IstrR,IendR
#  ifdef WET_DRY
          rmask_avg(i,j)=0.0_r8
#  else
          rmask_avg(i,j)=rmask_io(i,j)
#  endif
        END DO
      END DO

      DO j=JstrR,JendR
        DO i=Istr,IendR
#  ifdef WET_DRY
          umask_avg(i,j)=0.0_r8
#  else
          umask_avg(i,j)=umask_io(i,j)
#  endif
        END DO
      END DO

      DO j=Jstr,JendR
        DO i=IstrR,IendR
#  ifdef WET_DRY
          vmask_avg(i,j)=0.0_r8
#  else
          vmask_avg(i,j)=vmask_io(i,j)
#  endif
        END DO
      END DO
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_p2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          pmask_avg)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          rmask_avg)
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          umask_avg)
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          vmask_avg)
      END IF

#  ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 4,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    pmask_avg, rmask_avg, umask_avg, vmask_avg)
#  endif
# endif

# ifdef DIAGNOSTICS
!
!-----------------------------------------------------------------------
!  Initialize diagnostic file Land/Sea masks for time-averaged fields.
!-----------------------------------------------------------------------
!
      DO j=Jstr,JendR
        DO i=Istr,IendR
#  ifdef WET_DRY
          pmask_dia(i,j)=0.0_r8
#  else
          pmask_dia(i,j)=pmask_io(i,j)
#  endif
        END DO
      END DO

      DO j=JstrR,JendR
        DO i=IstrR,IendR
#  ifdef WET_DRY
          rmask_dia(i,j)=0.0_r8
#  else
          rmask_dia(i,j)=rmask_io(i,j)
#  endif
        END DO
      END DO

      DO j=JstrR,JendR
        DO i=Istr,IendR
#  ifdef WET_DRY
          umask_dia(i,j)=0.0_r8
#  else
          umask_dia(i,j)=umask_io(i,j)
#  endif
        END DO
      END DO

      DO j=Jstr,JendR
        DO i=IstrR,IendR
#  ifdef WET_DRY
          vmask_dia(i,j)=0.0_r8
#  else
          vmask_dia(i,j)=vmask_io(i,j)
#  endif
        END DO
      END DO
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_p2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          pmask_dia)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          rmask_dia)
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          umask_dia)
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          vmask_dia)
      END IF

#  ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 4,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    pmask_dia, rmask_dia, umask_dia, vmask_dia)
#  endif

# endif

      RETURN
      END SUBROUTINE set_masks_tile

# if defined WET_DRY     && \
     defined AVERAGES    || \
    (defined AD_AVERAGES && defined ADJOINT) || \
    (defined RP_AVERAGES && defined TL_IOMS) || \
    (defined TL_AVERAGES && defined TANGENT)
!
!***********************************************************************
      SUBROUTINE set_avg_masks (ng, tile, model,                        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          pmask_avg, rmask_avg,                   &
     &                          umask_avg, vmask_avg)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE exchange_2d_mod
#  ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
#  endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
#  ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: pmask_avg(LBi:,LBj:)
      real(r8), intent(inout) :: rmask_avg(LBi:,LBj:)
      real(r8), intent(inout) :: umask_avg(LBi:,LBj:)
      real(r8), intent(inout) :: vmask_avg(LBi:,LBj:)
#  else
      real(r8), intent(inout) :: pmask_avg(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: rmask_avg(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: umask_avg(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: vmask_avg(LBi:UBi,LBj:UBj)
#  endif
!
!
!  Local variable declarations.
!
      integer :: i, j

#  include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Return if time-averaging window is zero.
!-----------------------------------------------------------------------
!
      IF (nAVG(ng).eq.0) RETURN
!
!-----------------------------------------------------------------------
!  If last time-step of average window, convert time dependent counters
!  for wet points to time-averaged Land/Sea masks (dry=0, wet=1) for
!  the current average window period. Notice that a grid point is wet
!  if the count is greater than zero for the current time average
!  window.
!-----------------------------------------------------------------------
!
      IF ((iic(ng).gt.ntsAVG(ng)).and.                                  &
     &    (MOD(iic(ng)-1,nAVG(ng)).eq.0).and.                           &
     &    ((iic(ng).ne.ntstart(ng)).or.(nrrec(ng).eq.0))) THEN

        DO j=Jstr,JendR
          DO i=Istr,IendR
            pmask_avg(i,j)=MIN(1.0_r8, pmask_avg(i,j))
          END DO
        END DO
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            rmask_avg(i,j)=MIN(1.0_r8, rmask_avg(i,j))
          END DO
        END DO
        DO j=JstrR,JendR
          DO i=Istr,IendR
            umask_avg(i,j)=MIN(1.0_r8, umask_avg(i,j))
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            vmask_avg(i,j)=MIN(1.0_r8, vmask_avg(i,j))
          END DO
        END DO
!
        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
          CALL exchange_p2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            pmask_avg)
          CALL exchange_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            rmask_avg)
          CALL exchange_u2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            umask_avg)
          CALL exchange_v2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            vmask_avg)
        END IF

#  ifdef DISTRIBUTE
        CALL mp_exchange2d (ng, tile, model, 4,                         &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      pmask_avg, rmask_avg, umask_avg, vmask_avg)
#  endif

      END IF

      RETURN

      END SUBROUTINE set_avg_masks
# endif
#endif
      END MODULE set_masks_mod
