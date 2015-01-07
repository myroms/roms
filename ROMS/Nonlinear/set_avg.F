#include "cppdefs.h"
      MODULE set_avg_mod
#ifdef AVERAGES
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2015 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine accumulates and computes output time-averaged       !
!  fields.  Due to synchronization, the time-averaged fields are       !
!  computed in delayed mode. All averages are accumulated at the       !
!  beggining of the next time-step.                                    !
!                                                                      !
# if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
!  It computes leas-squares coefficients to detide time-averaged       !
!  fields. Notice that "set_detide" is called last since we need       !
!  the regular time-averages for those fields to detide.               !
!                                                                      !
# endif
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC :: set_avg

      CONTAINS
!
!***********************************************************************
      SUBROUTINE set_avg (ng, tile)
!***********************************************************************
!
      USE mod_param
# ifdef WET_DRY
      USE mod_grid
# endif
      USE mod_stepping
# if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
      USE mod_tides
# endif
# ifdef WET_DRY
!
      USE set_masks_mod, ONLY : set_avg_masks
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 5)
# endif
      CALL set_avg_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
# ifdef SOLVE3D
     &                   NOUT,                                          &
# endif
     &                   KOUT)

# if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
      CALL set_detide_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      NTC(ng), KOUT,                              &
#  ifdef SOLVE3D
     &                      NOUT,                                       &
#  endif
     &                      TIDES(ng) % CosOmega,                       &
     &                      TIDES(ng) % SinOmega,                       &
     &                      TIDES(ng) % CosW_avg,                       &
     &                      TIDES(ng) % CosW_sum,                       &
     &                      TIDES(ng) % SinW_avg,                       &
     &                      TIDES(ng) % SinW_sum,                       &
     &                      TIDES(ng) % CosWCosW,                       &
     &                      TIDES(ng) % SinWSinW,                       &
     &                      TIDES(ng) % SinWCosW)
# endif

# ifdef WET_DRY
      CALL set_avg_masks (ng, tile, iNLM,                               &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    GRID(ng) % pmask_avg,                         &
     &                    GRID(ng) % rmask_avg,                         &
     &                    GRID(ng) % umask_avg,                         &
     &                    GRID(ng) % vmask_avg)
# endif

# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 5)
# endif
      RETURN
      END SUBROUTINE set_avg
!
!***********************************************************************
      SUBROUTINE set_avg_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
# ifdef SOLVE3D
     &                         Nout,                                    &
# endif
     &                         Kout)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_average
# if defined FORWARD_WRITE && defined SOLVE3D
      USE mod_coupling
# endif
      USE mod_forces
      USE mod_grid
# ifdef SOLVE3D
      USE mod_mixing
# endif
      USE mod_ocean
      USE mod_scalars
# if defined SEDIMENT && defined BEDLOAD
      USE mod_sedbed
      USE mod_sediment
# endif
!
      USE exchange_2d_mod
# ifdef SOLVE3D
      USE exchange_3d_mod
# endif
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
#  ifdef SOLVE3D
      USE mp_exchange_mod, ONLY : mp_exchange3d
#  endif
# endif
      USE uv_rotate_mod, ONLY : uv_rotate2d
# ifdef SOLVE3D
      USE uv_rotate_mod, ONLY : uv_rotate3d
# endif
      USE vorticity_mod, ONLY : vorticity_tile
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Kout
# ifdef SOLVE3D
      integer, intent(in) :: Nout
# endif
!
!
!  Local variable declarations.
!
      integer :: i, it, j, k

      real(r8) :: fac

      real(r8) :: pfac(IminS:ImaxS,JminS:JmaxS)
      real(r8) :: rfac(IminS:ImaxS,JminS:JmaxS)
      real(r8) :: ufac(IminS:ImaxS,JminS:JmaxS)
      real(r8) :: vfac(IminS:ImaxS,JminS:JmaxS)

# ifdef SOLVE3D
      real(r8) :: potvor(LBi:UBi,LBj:UBj,N(ng))
      real(r8) :: relvor(LBi:UBi,LBj:UBj,N(ng))
# endif
      real(r8) :: potvor_bar(LBi:UBi,LBj:UBj)
      real(r8) :: relvor_bar(LBi:UBi,LBj:UBj)

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Return if time-averaging window is zero.
!-----------------------------------------------------------------------
!
      IF (nAVG(ng).eq.0) RETURN
!
!-----------------------------------------------------------------------
!  Compute vorticity fields.
!-----------------------------------------------------------------------
!
      IF (Aout(id2dPV,ng).or.Aout(id2dRV,ng).or.                        &
     &    Aout(id3dPV,ng).or.Aout(id3dRV,ng)) THEN
        CALL vorticity_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
# ifdef SOLVE3D
     &                       Kout, Nout,                                &
# else
     &                       Kout,                                      &
# endif
# ifdef MASKING
     &                       GRID(ng) % pmask,                          &
     &                       GRID(ng) % umask,                          &
     &                       GRID(ng) % vmask,                          &
# endif
     &                       GRID(ng) % fomn,                           &
     &                       GRID(ng) % h,                              &
     &                       GRID(ng) % om_u,                           &
     &                       GRID(ng) % on_v,                           &
     &                       GRID(ng) % pm,                             &
     &                       GRID(ng) % pn,                             &
# ifdef SOLVE3D
     &                       GRID(ng) % z_r,                            &
     &                       OCEAN(ng) % pden,                          &
     &                       OCEAN(ng) % u,                             &
     &                       OCEAN(ng) % v,                             &
# endif
     &                       OCEAN(ng) % ubar,                          &
     &                       OCEAN(ng) % vbar,                          &
     &                       OCEAN(ng) % zeta,                          &
# ifdef SOLVE3D
     &                       potvor, relvor,                            &
# endif
     &                       potvor_bar, relvor_bar)
      END IF
!
!-----------------------------------------------------------------------
!  Initialize time-averaged arrays when appropriate.  Notice that
!  fields are initilized twice during re-start.  However, the time-
!  averaged fields are computed correctly.
!-----------------------------------------------------------------------
!
      IF (((iic(ng).gt.ntsAVG(ng)).and.                                 &
     &     (MOD(iic(ng)-1,nAVG(ng)).eq.1)).or.                          &
     &    ((iic(ng).ge.ntsAVG(ng)).and.(nAVG(ng).eq.1)).or.             &
     &    ((nrrec(ng).gt.0).and.(iic(ng).eq.ntstart(ng)))) THEN

# ifdef WET_DRY
!
!  If wetting and drying, initialize time dependent counters for wet
!  points. The time averaged field at each point is only accumulated
!  over wet points since its multiplied by the appropriate mask.
!
        DO j=Jstr,JendR
          DO i=Istr,IendR
            GRID(ng)%pmask_avg(i,j)=MAX(0.0_r8,                         &
     &                                  MIN(GRID(ng)%pmask_full(i,j),   &
     &                                      1.0_r8))
          END DO
        END DO
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            GRID(ng)%rmask_avg(i,j)=MAX(0.0_r8,                         &
     &                                  MIN(GRID(ng)%rmask_full(i,j),   &
     &                                      1.0_r8))
          END DO
        END DO
        DO j=JstrR,JendR
          DO i=Istr,IendR
            GRID(ng)%umask_avg(i,j)=MAX(0.0_r8,                         &
     &                                  MIN(GRID(ng)%umask_full(i,j),   &
     &                                      1.0_r8))
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            GRID(ng)%vmask_avg(i,j)=MAX(0.0_r8,                         &
     &                                  MIN(GRID(ng)%vmask_full(i,j),   &
     &                                      1.0_r8))
          END DO
        END DO
# endif
!
!  Initialize state variables.
!
        IF (Aout(idFsur,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgzeta(i,j)=OCEAN(ng)%zeta(i,j,Kout)
# ifdef WET_DRY
              AVERAGE(ng)%avgzeta(i,j)=AVERAGE(ng)%avgzeta(i,j)*        &
     &                                 GRID(ng)%rmask_full(i,j)
# endif
            END DO
          END DO
        END IF

        IF (Aout(idUbar,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgu2d(i,j)=OCEAN(ng)%ubar(i,j,Kout)
# ifdef WET_DRY
              AVERAGE(ng)%avgu2d(i,j)=AVERAGE(ng)%avgu2d(i,j)*          &
     &                                GRID(ng)%umask_full(i,j)
# endif
            END DO
          END DO
        END IF
        IF (Aout(idVbar,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgv2d(i,j)=OCEAN(ng)%vbar(i,j,Kout)
# ifdef WET_DRY
              AVERAGE(ng)%avgv2d(i,j)=AVERAGE(ng)%avgv2d(i,j)*          &
     &                                GRID(ng)%vmask_full(i,j)
# endif
            END DO
          END DO
        END IF

        IF (Aout(idu2dE,ng).and.Aout(idv2dN,ng)) THEN
          CALL uv_rotate2d (ng, tile, .FALSE., .FALSE.,                 &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      GRID(ng) % CosAngler,                       &
     &                      GRID(ng) % SinAngler,                       &
# ifdef MASKING
     &                      GRID(ng)%rmask_full,                        &
# endif
     &                      OCEAN(ng) % ubar(:,:,Kout),                 &
     &                      OCEAN(ng) % vbar(:,:,Kout),                 &
     &                      AVERAGE(ng)%avgu2dE,                        &
     &                      AVERAGE(ng)%avgv2dN)
        END IF

# ifdef SOLVE3D
        IF (Aout(idUvel,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgu3d(i,j,k)=OCEAN(ng)%u(i,j,k,Nout)
#  ifdef WET_DRY
                AVERAGE(ng)%avgu3d(i,j,k)=AVERAGE(ng)%avgu3d(i,j,k)*    &
     &                                    GRID(ng)%umask_full(i,j)
#  endif
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idVvel,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgv3d(i,j,k)=OCEAN(ng)%v(i,j,k,Nout)
#  ifdef WET_DRY
                AVERAGE(ng)%avgv3d(i,j,k)=AVERAGE(ng)%avgv3d(i,j,k)*    &
     &                                    GRID(ng)%vmask_full(i,j)
#  endif
              END DO
            END DO
          END DO
        END IF

        IF (Aout(idu3dE,ng).and.Aout(idv3dN,ng)) THEN
          CALL uv_rotate3d (ng, tile, .FALSE., .FALSE.,                 &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      GRID(ng) % CosAngler,                       &
     &                      GRID(ng) % SinAngler,                       &
#  ifdef MASKING
     &                      GRID(ng)%rmask_full,                        &
#  endif
     &                      OCEAN(ng) % u(:,:,:,Nout),                  &
     &                      OCEAN(ng) % v(:,:,:,Nout),                  &
     &                      AVERAGE(ng)%avgu3dE,                        &
     &                      AVERAGE(ng)%avgv3dN)
        END IF

        IF (Aout(idOvel,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgw3d(i,j,k)=OCEAN(ng)%W(i,j,k)*           &
     &                                    GRID(ng)%pm(i,j)*             &
     &                                    GRID(ng)%pn(i,j)
#  ifdef WET_DRY
                AVERAGE(ng)%avgw3d(i,j,k)=AVERAGE(ng)%avgw3d(i,j,k)*    &
     &                                    GRID(ng)%rmask_full(i,j)
#  endif
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idWvel,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgwvel(i,j,k)=OCEAN(ng)%wvel(i,j,k)
#  ifdef WET_DRY
                AVERAGE(ng)%avgwvel(i,j,k)=AVERAGE(ng)%avgwvel(i,j,k)*  &
     &                                     GRID(ng)%rmask_full(i,j)
#  endif
              END DO
            END DO
          END DO
        END IF

        IF (Aout(idDano,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgrho(i,j,k)=OCEAN(ng)%rho(i,j,k)
#  ifdef WET_DRY
                AVERAGE(ng)%avgrho(i,j,k)=AVERAGE(ng)%avgrho(i,j,k)*    &
     &                                    GRID(ng)%rmask_full(i,j)
#  endif
              END DO
            END DO
          END DO
        END IF
        DO it=1,NT(ng)
          IF (Aout(idTvar(it),ng)) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgt(i,j,k,it)=OCEAN(ng)%t(i,j,k,Nout,it)
#  ifdef WET_DRY
                  AVERAGE(ng)%avgt(i,j,k,it)=AVERAGE(ng)%avgt(i,j,k,it)*&
     &                                       GRID(ng)%rmask_full(i,j)
#  endif
                END DO
              END DO
            END DO
          END IF
        END DO

#  if defined SEDIMENT && defined BEDLOAD
        DO it=1,NST
          IF (Aout(idUbld(it),ng)) THEN
            DO j=JstrR,JendR
              DO i=Istr,IendR
                SEDBED(ng)%avgbedldu(i,j,it)=SEDBED(ng)%bedldu(i,j,it)
#   ifdef WET_DRY
                SEDBED(ng)%avgbedldu(i,j,it)=SEDBED(ng)%bedldu(i,j,it)* &
     &                                       GRID(ng)%umask_full(i,j)
#   endif
              END DO
            END DO
          END IF

          IF (Aout(idVbld(it),ng)) THEN
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                SEDBED(ng)%avgbedldv(i,j,it)=SEDBED(ng)%bedldv(i,j,it)
#   ifdef WET_DRY
                SEDBED(ng)%avgbedldv(i,j,it)=SEDBED(ng)%bedldv(i,j,it)* &
     &                                       GRID(ng)%vmask_full(i,j)
#   endif
              END DO
            END DO
          END IF
        END DO
#  endif

#  if defined LMD_MIXING || defined MY25_MIXING || defined GLS_MIXING
        IF (Aout(idVvis,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgAKv(i,j,k)=MIXING(ng)%Akv(i,j,k)
#   ifdef WET_DRY
                AVERAGE(ng)%avgAKv(i,j,k)=AVERAGE(ng)%avgAKv(i,j,k)*    &
     &                                    GRID(ng)%rmask_full(i,j)
#   endif
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idTdif,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgAKt(i,j,k)=MIXING(ng)%Akt(i,j,k,itemp)
#   ifdef WET_DRY
                AVERAGE(ng)%avgAKt(i,j,k)=AVERAGE(ng)%avgAKt(i,j,k)*    &
     &                                    GRID(ng)%rmask_full(i,j)
#   endif
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idSdif,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgAKs(i,j,k)=MIXING(ng)%Akt(i,j,k,isalt)
#   ifdef WET_DRY
                AVERAGE(ng)%avgAKs(i,j,k)=AVERAGE(ng)%avgAKs(i,j,k)*    &
     &                                    GRID(ng)%rmask_full(i,j)
#   endif
              END DO
            END DO
          END DO
        END IF
#  endif
#  ifdef LMD_SKPP
        IF (Aout(idHsbl,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avghsbl(i,j)=MIXING(ng)%hsbl(i,j)
#   ifdef WET_DRY
              AVERAGE(ng)%avghsbl(i,j)=AVERAGE(ng)%avghsbl(i,j)*        &
     &                                 GRID(ng)%rmask_full(i,j)
#   endif
            END DO
          END DO
        END IF
#  endif
#  ifdef LMD_BKPP
        IF (Aout(idHbbl,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avghbbl(i,j)=MIXING(ng)%hbbl(i,j)
#   ifdef WET_DRY
              AVERAGE(ng)%avghbbl(i,j)=AVERAGE(ng)%avghbbl(i,j)*        &
     &                                 GRID(ng)%rmask_full(i,j)
#   endif
            END DO
          END DO
        END IF
#  endif
# endif

# if defined FORWARD_WRITE && defined SOLVE3D
!
!  Initialize 2D/3D coupling terms.
!
        IF (Aout(idUfx1,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgDU_avg1(i,j)=COUPLING(ng)%DU_avg1(i,j)
#  ifdef WET_DRY
              AVERAGE(ng)%avgDU_avg1(i,j)=AVERAGE(ng)%avgDU_avg1(i,j)*  &
     &                                    GRID(ng)%umask_full(i,j)
#  endif
            END DO
          END DO
        END IF
        IF (Aout(idUfx2,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgDU_avg2(i,j)=COUPLING(ng)%DU_avg2(i,j)
#  ifdef WET_DRY
              AVERAGE(ng)%avgDU_avg2(i,j)=AVERAGE(ng)%avgDU_avg2(i,j)*  &
     &                                    GRID(ng)%umask_full(i,j)
#  endif
            END DO
          END DO
        END IF

        IF (Aout(idVfx1,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgDV_avg1(i,j)=COUPLING(ng)%DV_avg1(i,j)
#  ifdef WET_DRY
              AVERAGE(ng)%avgDV_avg1(i,j)=AVERAGE(ng)%avgDV_avg1(i,j)*  &
     &                                    GRID(ng)%vmask_full(i,j)
#  endif
            END DO
          END DO
        END IF
        IF (Aout(idVfx2,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgDV_avg2(i,j)=COUPLING(ng)%DV_avg2(i,j)
#  ifdef WET_DRY
              AVERAGE(ng)%avgDV_avg2(i,j)=AVERAGE(ng)%avgDV_avg2(i,j)*  &
     &                                    GRID(ng)%vmask_full(i,j)
#  endif
            END DO
          END DO
        END IF
# endif
!
!  Initialize surface and bottom fluxes.
!
        IF (Aout(idUsms,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgsus(i,j)=FORCES(ng)%sustr(i,j)
# ifdef WET_DRY
              AVERAGE(ng)%avgsus(i,j)=AVERAGE(ng)%avgsus(i,j)*          &
     &                                GRID(ng)%umask_full(i,j)
# endif
            END DO
          END DO
        END IF
        IF (Aout(idVsms,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgsvs(i,j)=FORCES(ng)%svstr(i,j)
# ifdef WET_DRY
              AVERAGE(ng)%avgsvs(i,j)=AVERAGE(ng)%avgsvs(i,j)*          &
     &                                GRID(ng)%vmask_full(i,j)
# endif
            END DO
          END DO
        END IF

        IF (Aout(idUbms,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgbus(i,j)=FORCES(ng)%bustr(i,j)
# ifdef WET_DRY
              AVERAGE(ng)%avgbus(i,j)=AVERAGE(ng)%avgbus(i,j)*          &
     &                                GRID(ng)%umask_full(i,j)
# endif
            END DO
          END DO
        END IF
        IF (Aout(idVbms,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgbvs(i,j)=FORCES(ng)%bvstr(i,j)
# ifdef WET_DRY
              AVERAGE(ng)%avgbvs(i,j)=AVERAGE(ng)%avgbvs(i,j)*          &
     &                                GRID(ng)%vmask_full(i,j)
# endif
            END DO
          END DO
        END IF

# ifdef SOLVE3D
#  if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
        IF (Aout(idPair,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgPair(i,j)=FORCES(ng)%Pair(i,j)
#   ifdef WET_DRY
              AVERAGE(ng)%avgPair(i,j)=AVERAGE(ng)%avgPair(i,j)*        &
     &                                 GRID(ng)%rmask_full(i,j)
#   endif
            END DO
          END DO
        END IF
#  endif
#  if defined BULK_FLUXES || defined ECOSIM
        IF (Aout(idUair,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgUwind(i,j)=FORCES(ng)%Uwind(i,j)
#   ifdef WET_DRY
              AVERAGE(ng)%avgUwind(i,j)=AVERAGE(ng)%avgUwind(i,j)*      &
     &                                  GRID(ng)%rmask_full(i,j)
#   endif
            END DO
          END DO
        END IF

        IF (Aout(idVair,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgVwind(i,j)=FORCES(ng)%Vwind(i,j)
#   ifdef WET_DRY
              AVERAGE(ng)%avgVwind(i,j)=AVERAGE(ng)%avgVwind(i,j)*      &
     &                                  GRID(ng)%rmask_full(i,j)
#   endif
            END DO
          END DO
        END IF
#  endif

        IF (Aout(idTsur(itemp),ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgstf(i,j)=FORCES(ng)%stflx(i,j,itemp)
#  ifdef WET_DRY
              AVERAGE(ng)%avgstf(i,j)=AVERAGE(ng)%avgstf(i,j)*          &
     &                                GRID(ng)%rmask_full(i,j)
#  endif
            END DO
          END DO
        END IF
        IF (Aout(idTsur(isalt),ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgswf(i,j)=FORCES(ng)%stflx(i,j,isalt)
#  ifdef WET_DRY
              AVERAGE(ng)%avgswf(i,j)=AVERAGE(ng)%avgswf(i,j)*          &
     &                                GRID(ng)%rmask_full(i,j)
#  endif
            END DO
          END DO
        END IF
#  ifdef SHORTWAVE
        IF (Aout(idSrad,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgsrf(i,j)=FORCES(ng)%srflx(i,j)
#   ifdef WET_DRY
              AVERAGE(ng)%avgsrf(i,j)=AVERAGE(ng)%avgsrf(i,j)*          &
     &                                GRID(ng)%rmask_full(i,j)
#   endif
            END DO
          END DO
        END IF
#  endif
#  ifdef BULK_FLUXES
        IF (Aout(idLhea,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avglhf(i,j)=FORCES(ng)%lhflx(i,j)
#   ifdef WET_DRY
              AVERAGE(ng)%avglhf(i,j)=AVERAGE(ng)%avglhf(i,j)*          &
     &                                GRID(ng)%rmask_full(i,j)
#   endif
            END DO
          END DO
        END IF
        IF (Aout(idLrad,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avglrf(i,j)=FORCES(ng)%lrflx(i,j)
#   ifdef WET_DRY
              AVERAGE(ng)%avglrf(i,j)=AVERAGE(ng)%avglrf(i,j)*          &
     &                                GRID(ng)%rmask_full(i,j)
#   endif
            END DO
          END DO
        END IF
        IF (Aout(idShea,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgshf(i,j)=FORCES(ng)%shflx(i,j)
#   ifdef WET_DRY
              AVERAGE(ng)%avgshf(i,j)=AVERAGE(ng)%avgshf(i,j)*          &
     &                                GRID(ng)%rmask_full(i,j)
#   endif
            END DO
          END DO
        END IF
#   ifdef EMINUSP
        IF (Aout(idevap,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgevap(i,j)=FORCES(ng)%evap(i,j)
#    ifdef WET_DRY
              AVERAGE(ng)%avgevap(i,j)=AVERAGE(ng)%avgevap(i,j)*        &
     &                                 GRID(ng)%rmask_full(i,j)
#    endif
            END DO
          END DO
        END IF
        IF (Aout(idrain,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgrain(i,j)=FORCES(ng)%rain(i,j)
#    ifdef WET_DRY
              AVERAGE(ng)%avgrain(i,j)=AVERAGE(ng)%avgrain(i,j)*        &
     &                                 GRID(ng)%rmask_full(i,j)
#    endif
            END DO
          END DO
        END IF
#   endif
#  endif
# endif

# ifdef NEARSHORE_MELLOR
!
!  Initialize radiation stresses.
!
        IF (Aout(idU2Sd,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgu2Sd(i,j)=OCEAN(ng)%ubar_stokes(i,j)
#  ifdef WET_DRY
              AVERAGE(ng)%avgu2Sd(i,j)=AVERAGE(ng)%avgu2Sd(i,j)*        &
     &                                 GRID(ng)%umask_full(i,j)
#  endif
            END DO
          END DO
        END IF
        IF (Aout(idV2Sd,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgv2Sd(i,j)=OCEAN(ng)%vbar_stokes(i,j)
#  ifdef WET_DRY
              AVERAGE(ng)%avgv2Sd(i,j)=AVERAGE(ng)%avgv2Sd(i,j)*        &
     &                                 GRID(ng)%vmask_full(i,j)
#  endif
            END DO
          END DO
        END IF

        IF (Aout(idU2rs,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgu2RS(i,j)=MIXING(ng)%rustr2d(i,j)
#  ifdef WET_DRY
              AVERAGE(ng)%avgu2RS(i,j)=AVERAGE(ng)%avgu2RS(i,j)*        &
     &                                 GRID(ng)%umask_full(i,j)
#  endif
            END DO
          END DO
        END IF
        IF (Aout(idV2rs,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgv2RS(i,j)=MIXING(ng)%rvstr2d(i,j)
#  ifdef WET_DRY
              AVERAGE(ng)%avgv2RS(i,j)=AVERAGE(ng)%avgv2RS(i,j)*        &
     &                                 GRID(ng)%vmask_full(i,j)
#  endif
            END DO
          END DO
        END IF

        IF (Aout(idW2xx,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgSxx2d(i,j)=MIXING(ng)%Sxx_bar(i,j)
#  ifdef WET_DRY
              AVERAGE(ng)%avgSxx2d(i,j)=AVERAGE(ng)%avgSxx2d(i,j)*      &
     &                                  GRID(ng)%rmask_full(i,j)
#  endif
            END DO
          END DO
        END IF
        IF (Aout(idW2xy,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgSxy2d(i,j)=MIXING(ng)%Sxy_bar(i,j)
#  ifdef WET_DRY
              AVERAGE(ng)%avgSxy2d(i,j)=AVERAGE(ng)%avgSxy2d(i,j)*      &
     &                                  GRID(ng)%rmask_full(i,j)
#  endif
            END DO
          END DO
        END IF
        IF (Aout(idW2yy,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgSyy2d(i,j)=MIXING(ng)%Syy_bar(i,j)
#  ifdef WET_DRY
              AVERAGE(ng)%avgSyy2d(i,j)=AVERAGE(ng)%avgSyy2d(i,j)*      &
     &                                  GRID(ng)%rmask_full(i,j)
#  endif
            END DO
          END DO
        END IF

#  ifdef SOLVE3D
        IF (Aout(idU3Sd,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgu3Sd(i,j,k)=OCEAN(ng)%u_stokes(i,j,k)
#   ifdef WET_DRY
                AVERAGE(ng)%avgu3Sd(i,j,k)=AVERAGE(ng)%avgu3Sd(i,j,k)*  &
     &                                     GRID(ng)%umask_full(i,j)
#   endif
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idV3Sd,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgv3Sd(i,j,k)=OCEAN(ng)%v_stokes(i,j,k)
#   ifdef WET_DRY
                AVERAGE(ng)%avgv3Sd(i,j,k)=AVERAGE(ng)%avgv3Sd(i,j,k)*  &
     &                                     GRID(ng)%vmask_full(i,j)
#   endif
              END DO
            END DO
          END DO
        END IF

        IF (Aout(idU3rs,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgu3RS(i,j,k)=MIXING(ng)%rustr3d(i,j,k)
#   ifdef WET_DRY
                AVERAGE(ng)%avgu3RS(i,j,k)=AVERAGE(ng)%avgu3RS(i,j,k)*  &
     &                                     GRID(ng)%umask_full(i,j)
#   endif
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idV3rs,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgv3RS(i,j,k)=MIXING(ng)%rvstr3d(i,j,k)
#   ifdef WET_DRY
                AVERAGE(ng)%avgv3RS(i,j,k)=AVERAGE(ng)%avgv3RS(i,j,k)*  &
     &                                     GRID(ng)%vmask_full(i,j)
#   endif
              END DO
            END DO
          END DO
        END IF

        IF (Aout(idW3xx,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgSxx3d(i,j,k)=MIXING(ng)%Sxx(i,j,k)
#   ifdef WET_DRY
                AVERAGE(ng)%avgSxx3d(i,j,k)=AVERAGE(ng)%avgSxx3d(i,j,k)*&
     &                                      GRID(ng)%rmask_full(i,j)
#   endif
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idW3xy,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgSxy3d(i,j,k)=MIXING(ng)%Sxy(i,j,k)
#   ifdef WET_DRY
                AVERAGE(ng)%avgSxy3d(i,j,k)=AVERAGE(ng)%avgSxy3d(i,j,k)*&
     &                                      GRID(ng)%rmask_full(i,j)
#   endif
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idW3yy,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgSyy3d(i,j,k)=MIXING(ng)%Syy(i,j,k)
#   ifdef WET_DRY
                AVERAGE(ng)%avgSyy3d(i,j,k)=AVERAGE(ng)%avgSyy3d(i,j,k)*&
     &                                      GRID(ng)%rmask_full(i,j)
#   endif
              END DO
            END DO
          END DO
        END IF
#   ifdef NEARSHORE_MELLOR05
        IF (Aout(idW3zx,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgSzx3d(i,j,k)=MIXING(ng)%Szx(i,j,k)
#    ifdef WET_DRY
                AVERAGE(ng)%avgSzx3d(i,j,k)=AVERAGE(ng)%avgSzx3d(i,j,k)*&
     &                                      GRID(ng)%rmask_full(i,j)
#    endif
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idW3zy,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgSzy3d(i,j,k)=MIXING(ng)%Szy(i,j,k)
#    ifdef WET_DRY
                AVERAGE(ng)%avgSzy3d(i,j,k)=AVERAGE(ng)%avgSzy3d(i,j,k)*&
     &                                      GRID(ng)%rmask_full(i,j)
#    endif
              END DO
            END DO
          END DO
        END IF
#   endif
#  endif
# endif
!
!  Initialize vorticity fields.
!
        IF (Aout(id2dPV,ng)) THEN
          DO j=Jstr,Jend
            DO i=Istr,Iend
              AVERAGE(ng)%avgpvor2d(i,j)=potvor_bar(i,j)
# ifdef WET_DRY
              AVERAGE(ng)%avgpvor2d(i,j)=AVERAGE(ng)%avgpvor2d(i,j)*    &
     &                                   GRID(ng)%pmask_full(i,j)
# endif
            END DO
          END DO
        END IF
        IF (Aout(id2dRV,ng)) THEN
          DO j=Jstr,Jend
            DO i=Istr,Iend
              AVERAGE(ng)%avgrvor2d(i,j)=relvor_bar(i,j)
# ifdef WET_DRY
              AVERAGE(ng)%avgrvor2d(i,j)=AVERAGE(ng)%avgrvor2d(i,j)*    &
     &                                   GRID(ng)%pmask_full(i,j)
# endif
            END DO
          END DO
        END IF
# ifdef SOLVE3D
        IF (Aout(id3dPV,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              DO i=Istr,Iend
                AVERAGE(ng)%avgpvor3d(i,j,k)=potvor(i,j,k)
#  ifdef WET_DRY
                AVERAGE(ng)%avgpvor3d(i,j,k)=AVERAGE(ng)%avgpvor3d(i,j, &
     &                                                             k)*  &
     &                                       GRID(ng)%pmask_full(i,j)
#  endif
              END DO
            END DO
          END DO
        END IF
        IF (Aout(id3dRV,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              DO i=Istr,Iend
                AVERAGE(ng)%avgrvor3d(i,j,k)=relvor(i,j,k)
#  ifdef WET_DRY
                AVERAGE(ng)%avgrvor3d(i,j,k)=AVERAGE(ng)%avgrvor3d(i,j, &
     &                                                             k)*  &
     &                                       GRID(ng)%pmask_full(i,j)
#  endif
              END DO
            END DO
          END DO
        END IF
# endif
!
!  Initialize quadratic fields.
!
        IF (Aout(idZZav,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgZZ(i,j)=OCEAN(ng)%zeta(i,j,Kout)*          &
     &                               OCEAN(ng)%zeta(i,j,Kout)
# ifdef WET_DRY
              AVERAGE(ng)%avgZZ(i,j)=AVERAGE(ng)%avgZZ(i,j)*            &
     &                               GRID(ng)%rmask_full(i,j)
# endif
            END DO
          END DO
        END IF
        IF (Aout(idU2av,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgU2(i,j)=OCEAN(ng)%ubar(i,j,Kout)*          &
     &                               OCEAN(ng)%ubar(i,j,Kout)
# ifdef WET_DRY
              AVERAGE(ng)%avgU2(i,j)=AVERAGE(ng)%avgU2(i,j)*            &
     &                               GRID(ng)%umask_full(i,j)
# endif
            END DO
          END DO
        END IF
        IF (Aout(idV2av,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgV2(i,j)=OCEAN(ng)%vbar(i,j,Kout)*          &
     &                               OCEAN(ng)%vbar(i,j,Kout)
# ifdef WET_DRY
              AVERAGE(ng)%avgV2(i,j)=AVERAGE(ng)%avgV2(i,j)*            &
     &                               GRID(ng)%vmask_full(i,j)
# endif
            END DO
          END DO
        END IF

# ifdef SOLVE3D
        IF (Aout(idUUav,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgUU(i,j,k)=OCEAN(ng)%u(i,j,k,Nout)*       &
     &                                   OCEAN(ng)%u(i,j,k,Nout)
#  ifdef WET_DRY
                AVERAGE(ng)%avgUU(i,j,k)=AVERAGE(ng)%avgUU(i,j,k)*      &
     &                                   GRID(ng)%umask_full(i,j)
#  endif
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idVVav,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgVV(i,j,k)=OCEAN(ng)%v(i,j,k,Nout)*       &
     &                                   OCEAN(ng)%v(i,j,k,Nout)
#  ifdef WET_DRY
                AVERAGE(ng)%avgVV(i,j,k)=AVERAGE(ng)%avgVV(i,j,k)*      &
     &                                   GRID(ng)%vmask_full(i,j)
#  endif
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idUVav,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              DO i=Istr,Iend
                AVERAGE(ng)%avgUV(i,j,k)=0.25_r8*                       &
     &                                   (OCEAN(ng)%u(i  ,j  ,k,Nout)+  &
     &                                    OCEAN(ng)%u(i+1,j  ,k,Nout))* &
     &                                   (OCEAN(ng)%v(i  ,j  ,k,Nout)+  &
     &                                    OCEAN(ng)%v(i  ,j+1,k,Nout))
#  ifdef WET_DRY
                AVERAGE(ng)%avgUV(i,j,k)=AVERAGE(ng)%avgUV(i,j,k)*      &
     &                                   GRID(ng)%rmask_full(i,j)
#  endif
              END DO
            END DO
          END DO
        END IF

        IF (Aout(idHUav,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgHuon(i,j,k)=GRID(ng)%Huon(i,j,k)
#  ifdef WET_DRY
                AVERAGE(ng)%avgHuon(i,j,k)=AVERAGE(ng)%avgHuon(i,j,k)*  &
     &                                     GRID(ng)%umask_full(i,j)
#  endif
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idHVav,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgHvom(i,j,k)=GRID(ng)%Hvom(i,j,k)
#  ifdef WET_DRY
                AVERAGE(ng)%avgHvom(i,j,k)=AVERAGE(ng)%avgHvom(i,j,k)*  &
     &                                     GRID(ng)%vmask_full(i,j)
#  endif
              END DO
            END DO
          END DO
        END IF

        DO it=1,NT(ng)
          IF (Aout(idTTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgTT(i,j,k,it)=OCEAN(ng)%t(i,j,k,        &
     &                                                    Nout,it)*     &
     &                                        OCEAN(ng)%t(i,j,k,        &
     &                                                    Nout,it)
#  ifdef WET_DRY
                  AVERAGE(ng)%avgTT(i,j,k,it)=AVERAGE(ng)%avgTT(i,j,k,  &
     &                                                          it)*    &
     &                                        GRID(ng)%rmask_full(i,j)
#  endif
                END DO
              END DO
            END DO
          END IF
          IF (Aout(idUTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=Istr,Iend
                  AVERAGE(ng)%avgUT(i,j,k,it)=0.5_r8*                   &
     &                                        OCEAN(ng)%u(i,j,k,Nout)*  &
     &                                        (OCEAN(ng)%t(i-1,j,k,     &
     &                                                     Nout,it)+    &
     &                                         OCEAN(ng)%t(i  ,j,k,     &
     &                                                     Nout,it))
#  ifdef WET_DRY
                  AVERAGE(ng)%avgUT(i,j,k,it)=AVERAGE(ng)%avgUT(i,j,k,  &
     &                                                          it)*    &
     &                                        GRID(ng)%umask_full(i,j)
#  endif
                END DO
              END DO
            END DO
          END IF
          IF (Aout(idVTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=Jstr,Jend
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgVT(i,j,k,it)=0.5_r8*                   &
     &                                        OCEAN(ng)%v(i,j,k,Nout)*  &
     &                                        (OCEAN(ng)%t(i,j-1,k,     &
     &                                                     Nout,it)+    &
     &                                         OCEAN(ng)%t(i,j  ,k,     &
     &                                                     Nout,it))
#  ifdef WET_DRY
                  AVERAGE(ng)%avgVT(i,j,k,it)=AVERAGE(ng)%avgVT(i,j,k,  &
     &                                                          it)*    &
     &                                        GRID(ng)%vmask_full(i,j)
#  endif
                END DO
              END DO
            END DO
          END IF

          IF (Aout(iHUTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=Istr,Iend
                  AVERAGE(ng)%avgHuonT(i,j,k,it)=0.5_r8*                &
     &                                           GRID(ng)%Huon(i,j,k)*  &
     &                                           (OCEAN(ng)%t(i-1,j,k,  &
     &                                                        Nout,it)+ &
     &                                            OCEAN(ng)%t(i  ,j,k,  &
     &                                                        Nout,it))
#  ifdef WET_DRY
                  AVERAGE(ng)%avgHuonT(i,j,k,it)=AVERAGE(ng)%avgHuonT   &
     &                                                      (i,j,k,it)* &
     &                                           GRID(ng)%umask_full(i, &
     &                                                               j)
#  endif
                END DO
              END DO
            END DO
          END IF
          IF (Aout(iHVTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=Jstr,Jend
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgHvomT(i,j,k,it)=0.5_r8*                &
     &                                           GRID(ng)%Hvom(i,j,k)*  &
     &                                           (OCEAN(ng)%t(i,j-1,k,  &
     &                                                        Nout,it)+ &
     &                                            OCEAN(ng)%t(i,j  ,k,  &
     &                                                        Nout,it))
#  ifdef WET_DRY
                  AVERAGE(ng)%avgHvomT(i,j,k,it)=AVERAGE(ng)%avgHvomT   &
     &                                                      (i,j,k,it)* &
     &                                           GRID(ng)%vmask_full(i, &
     &                                                               j)
#  endif
                END DO
              END DO
            END DO
          END IF
        END DO
# endif
!
!-----------------------------------------------------------------------
!  Accumulate time-averaged fields.
!-----------------------------------------------------------------------
!
      ELSE IF (iic(ng).gt.ntsAVG(ng)) THEN

# ifdef WET_DRY
!
!  If wetting and drying, accumulate wet points counters.
!  points. The time averaged field at each point is only accumulated
!  over wet points since its multiplied by the appropriate mask.
!
        DO j=Jstr,JendR
          DO i=Istr,IendR
            GRID(ng)%pmask_avg(i,j)=GRID(ng)%pmask_avg(i,j)+            &
     &                              MAX(0.0_r8,                         &
     &                                  MIN(GRID(ng)%pmask_full(i,j),   &
     &                                      1.0_r8))
          END DO
        END DO
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            GRID(ng)%rmask_avg(i,j)=GRID(ng)%rmask_avg(i,j)+            &
     &                              MAX(0.0_r8,                         &
     &                                  MIN(GRID(ng)%rmask_full(i,j),   &
     &                                      1.0_r8))
          END DO
        END DO
        DO j=JstrR,JendR
          DO i=Istr,IendR
            GRID(ng)%umask_avg(i,j)=GRID(ng)%umask_avg(i,j)+            &
     &                              MAX(0.0_r8,                         &
     &                                  MIN(GRID(ng)%umask_full(i,j),   &
     &                                      1.0_r8))
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            GRID(ng)%vmask_avg(i,j)=GRID(ng)%vmask_avg(i,j)+            &
     &                              MAX(0.0_r8,                         &
     &                                  MIN(GRID(ng)%vmask_full(i,j),   &
     &                                      1.0_r8))
          END DO
        END DO
# endif
!
!  Accumulate state variables.
!
        IF (Aout(idFsur,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgzeta(i,j)=AVERAGE(ng)%avgzeta(i,j)+        &
# ifdef WET_DRY
     &                                 GRID(ng)%rmask_full(i,j)*        &
# endif
     &                                 OCEAN(ng)%zeta(i,j,Kout)
            END DO
          END DO
        END IF

        IF (Aout(idUbar,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgu2d(i,j)=AVERAGE(ng)%avgu2d(i,j)+          &
# ifdef WET_DRY
     &                                GRID(ng)%umask_full(i,j)*         &
# endif
     &                                OCEAN(ng)%ubar(i,j,Kout)
            END DO
          END DO
        END IF
        IF (Aout(idVbar,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgv2d(i,j)=AVERAGE(ng)%avgv2d(i,j)+          &
# ifdef WET_DRY
     &                                GRID(ng)%vmask_full(i,j)*         &
# endif
     &                                OCEAN(ng)%vbar(i,j,Kout)
            END DO
          END DO
        END IF

        IF (Aout(idu2dE,ng).and.Aout(idv2dN,ng)) THEN
          CALL uv_rotate2d (ng, tile, .TRUE., .FALSE.,                  &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      GRID(ng) % CosAngler,                       &
     &                      GRID(ng) % SinAngler,                       &
# ifdef MASKING
     &                      GRID(ng)%rmask_full,                        &
# endif
     &                      OCEAN(ng) % ubar(:,:,Kout),                 &
     &                      OCEAN(ng) % vbar(:,:,Kout),                 &
     &                      AVERAGE(ng)%avgu2dE,                        &
     &                      AVERAGE(ng)%avgv2dN)
        END IF

# ifdef SOLVE3D
        IF (Aout(idUvel,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgu3d(i,j,k)=AVERAGE(ng)%avgu3d(i,j,k)+    &
#  ifdef WET_DRY
     &                                    GRID(ng)%umask_full(i,j)*     &
#  endif
     &                                    OCEAN(ng)%u(i,j,k,Nout)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idVvel,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgv3d(i,j,k)=AVERAGE(ng)%avgv3d(i,j,k)+    &
#  ifdef WET_DRY
     &                                    GRID(ng)%vmask_full(i,j)*     &
#  endif
     &                                    OCEAN(ng)%v(i,j,k,Nout)
              END DO
            END DO
          END DO
        END IF

        IF (Aout(idu3dE,ng).and.Aout(idv3dN,ng)) THEN
          CALL uv_rotate3d (ng, tile, .TRUE., .FALSE.,                  &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      GRID(ng) % CosAngler,                       &
     &                      GRID(ng) % SinAngler,                       &
#  ifdef MASKING
     &                      GRID(ng)%rmask_full,                        &
#  endif
     &                      OCEAN(ng) % u(:,:,:,Nout),                  &
     &                      OCEAN(ng) % v(:,:,:,Nout),                  &
     &                      AVERAGE(ng)%avgu3dE,                        &
     &                      AVERAGE(ng)%avgv3dN)
        END IF

        IF (Aout(idOvel,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgw3d(i,j,k)=AVERAGE(ng)%avgw3d(i,j,k)+    &
#  ifdef WET_DRY
     &                                    GRID(ng)%rmask_full(i,j)*     &
#  endif
     &                                    OCEAN(ng)%W(i,j,k)*           &
     &                                    GRID(ng)%pm(i,j)*             &
     &                                    GRID(ng)%pn(i,j)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idWvel,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgwvel(i,j,k)=AVERAGE(ng)%avgwvel(i,j,k)+  &
#  ifdef WET_DRY
     &                                     GRID(ng)%rmask_full(i,j)*    &
#  endif
     &                                     OCEAN(ng)%wvel(i,j,k)
              END DO
            END DO
          END DO
        END IF

        IF (Aout(idDano,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgrho(i,j,k)=AVERAGE(ng)%avgrho(i,j,k)+    &
#  ifdef WET_DRY
     &                                    GRID(ng)%rmask_full(i,j)*     &
#  endif
     &                                    OCEAN(ng)%rho(i,j,k)
              END DO
            END DO
          END DO
        END IF
        DO it=1,NT(ng)
          IF (Aout(idTvar(it),ng)) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgt(i,j,k,it)=AVERAGE(ng)%avgt(i,j,k,it)+&
#  ifdef WET_DRY
     &                                       GRID(ng)%rmask_full(i,j)*  &
#  endif
     &                                       OCEAN(ng)%t(i,j,k,Nout,it)
                END DO
              END DO
            END DO
          END IF
        END DO

#  if defined SEDIMENT && defined BEDLOAD
        DO it=1,NST
          IF (Aout(idUbld(it),ng)) THEN
            DO j=JstrR,JendR
              DO i=Istr,IendR
                SEDBED(ng)%avgbedldu(i,j,it)=SEDBED(ng)%avgbedldu(i,j,  &
     &                                                            it)+  &
#   ifdef WET_DRY
     &                                       GRID(ng)%umask_full(i,j)*  &
#   endif
     &                                       SEDBED(ng)%bedldu(i,j,it)
              END DO
            END DO
          END IF

          IF (Aout(idVbld(it),ng)) THEN
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                SEDBED(ng)%avgbedldv(i,j,it)=SEDBED(ng)%avgbedldv(i,j,  &
     &                                                            it)+  &
#   ifdef WET_DRY
     &                                       GRID(ng)%vmask_full(i,j)*  &
#   endif
     &                                       SEDBED(ng)%bedldv(i,j,it)
              END DO
            END DO
          END IF
        END DO
#  endif

#  if defined LMD_MIXING || defined MY25_MIXING || defined GLS_MIXING
        IF (Aout(idVvis,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgAKv(i,j,k)=AVERAGE(ng)%avgAKv(i,j,k)+    &
#   ifdef WET_DRY
     &                                    GRID(ng)%rmask_full(i,j)*     &
#   endif
     &                                    MIXING(ng)%Akv(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idTdif,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgAKt(i,j,k)=AVERAGE(ng)%avgAKt(i,j,k)+    &
#   ifdef WET_DRY
     &                                    GRID(ng)%rmask_full(i,j)*     &
#   endif
     &                                    MIXING(ng)%Akt(i,j,k,itemp)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idSdif,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgAKs(i,j,k)=AVERAGE(ng)%avgAKs(i,j,k)+    &
#   ifdef WET_DRY
     &                                    GRID(ng)%rmask_full(i,j)*     &
#   endif
     &                                    MIXING(ng)%Akt(i,j,k,isalt)
              END DO
            END DO
          END DO
        END IF
#  endif
#  ifdef LMD_SKPP
        IF (Aout(idHsbl,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avghsbl(i,j)=AVERAGE(ng)%avghsbl(i,j)+        &
#   ifdef WET_DRY
     &                                 GRID(ng)%rmask_full(i,j)*        &
#   endif
     &                                 MIXING(ng)%hsbl(i,j)
            END DO
          END DO
        END IF
#  endif
#  ifdef LMD_BKPP
        IF (Aout(idHbbl,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avghbbl(i,j)=AVERAGE(ng)%avghbbl(i,j)+        &
#   ifdef WET_DRY
     &                                 GRID(ng)%rmask_full(i,j)*        &
#   endif
     &                                 MIXING(ng)%hbbl(i,j)
            END DO
          END DO
        END IF
#  endif
# endif

# if defined FORWARD_WRITE && defined SOLVE3D
!
!  Accumulate 2D/3D coupling terms.
!
        IF (Aout(idUfx1,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgDU_avg1(i,j)=AVERAGE(ng)%avgDU_avg1(i,j)+  &
#  ifdef WET_DRY
     &                                    GRID(ng)%umask_full(i,j)*     &
#  endif
     &                                    COUPLING(ng)%DU_avg1(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idUfx2,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgDU_avg2(i,j)=AVERAGE(ng)%avgDU_avg2(i,j)+  &
#  ifdef WET_DRY
     &                                    GRID(ng)%umask_full(i,j)*     &
#  endif
     &                                    COUPLING(ng)%DU_avg2(i,j)
            END DO
          END DO
        END IF

        IF (Aout(idVfx1,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgDV_avg1(i,j)=AVERAGE(ng)%avgDV_avg1(i,j)+  &
#  ifdef WET_DRY
     &                                    GRID(ng)%vmask_full(i,j)*     &
#  endif
     &                                    COUPLING(ng)%DV_avg1(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idVfx2,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgDV_avg2(i,j)=AVERAGE(ng)%avgDV_avg2(i,j)+  &
#  ifdef WET_DRY
     &                                    GRID(ng)%vmask_full(i,j)*     &
#  endif
     &                                    COUPLING(ng)%DV_avg2(i,j)
            END DO
          END DO
        END IF
# endif
!
!  Accumulate surface and bottom fluxes.
!
        IF (Aout(idUsms,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgsus(i,j)=AVERAGE(ng)%avgsus(i,j)+          &
# ifdef WET_DRY
     &                                GRID(ng)%umask_full(i,j)*         &
# endif
     &                                FORCES(ng)%sustr(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idVsms,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgsvs(i,j)=AVERAGE(ng)%avgsvs(i,j)+          &
# ifdef WET_DRY
     &                                GRID(ng)%vmask_full(i,j)*         &
# endif
     &                                FORCES(ng)%svstr(i,j)
            END DO
          END DO
        END IF

        IF (Aout(idUbms,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgbus(i,j)=AVERAGE(ng)%avgbus(i,j)+          &
# ifdef WET_DRY
     &                                GRID(ng)%umask_full(i,j)*         &
# endif
     &                                FORCES(ng)%bustr(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idVbms,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgbvs(i,j)=AVERAGE(ng)%avgbvs(i,j)+          &
# ifdef WET_DRY
     &                                GRID(ng)%vmask_full(i,j)*         &
# endif
     &                                FORCES(ng)%bvstr(i,j)
            END DO
          END DO
        END IF

# ifdef SOLVE3D

#  if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
        IF (Aout(idPair,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgPair(i,j)=AVERAGE(ng)%avgPair(i,j)+        &
#   ifdef WET_DRY
     &                                 GRID(ng)%rmask_full(i,j)*        &
#   endif
     &                                 FORCES(ng)%Pair(i,j)
            END DO
          END DO
        END IF
#  endif
#  if defined BULK_FLUXES || defined ECOSIM
        IF (Aout(idUair,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgUwind(i,j)=AVERAGE(ng)%avgUwind(i,j)+      &
#   ifdef WET_DRY
     &                                  GRID(ng)%rmask_full(i,j)*       &
#   endif
     &                                  FORCES(ng)%Uwind(i,j)
            END DO
          END DO
        END IF

        IF (Aout(idVair,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgVwind(i,j)=AVERAGE(ng)%avgVwind(i,j)+      &
#   ifdef WET_DRY
     &                                  GRID(ng)%rmask_full(i,j)*       &
#   endif
     &                                  FORCES(ng)%Vwind(i,j)
            END DO
          END DO
        END IF
#  endif

        IF (Aout(idTsur(itemp),ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgstf(i,j)=AVERAGE(ng)%avgstf(i,j)+          &
#  ifdef WET_DRY
     &                                GRID(ng)%rmask_full(i,j)*         &
#  endif
     &                                FORCES(ng)%stflx(i,j,itemp)
            END DO
          END DO
        END IF
        IF (Aout(idTsur(isalt),ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgswf(i,j)=AVERAGE(ng)%avgswf(i,j)+          &
#  ifdef WET_DRY
     &                                GRID(ng)%rmask_full(i,j)*         &
#  endif
     &                                FORCES(ng)%stflx(i,j,isalt)
            END DO
          END DO
        END IF
#  ifdef SHORTWAVE
        IF (Aout(idSrad,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgsrf(i,j)=AVERAGE(ng)%avgsrf(i,j)+          &
#   ifdef WET_DRY
     &                                GRID(ng)%rmask_full(i,j)*         &
#   endif
     &                                FORCES(ng)%srflx(i,j)
            END DO
          END DO
        END IF
#  endif
#  ifdef BULK_FLUXES
        IF (Aout(idLhea,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avglhf(i,j)=AVERAGE(ng)%avglhf(i,j)+          &
#   ifdef WET_DRY
     &                                GRID(ng)%rmask_full(i,j)*         &
#   endif
     &                                FORCES(ng)%lhflx(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idShea,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgshf(i,j)=AVERAGE(ng)%avgshf(i,j)+          &
#   ifdef WET_DRY
     &                                GRID(ng)%rmask_full(i,j)*         &
#   endif
     &                                FORCES(ng)%shflx(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idLrad,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avglrf(i,j)=AVERAGE(ng)%avglrf(i,j)+          &
#   ifdef WET_DRY
     &                                GRID(ng)%rmask_full(i,j)*         &
#   endif
     &                                FORCES(ng)%lrflx(i,j)
            END DO
          END DO
        END IF
#   ifdef EMINUSP
        IF (Aout(idevap,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgevap(i,j)=AVERAGE(ng)%avgevap(i,j)+        &
#    ifdef WET_DRY
     &                                 GRID(ng)%rmask_full(i,j)*        &
#    endif
     &                                 FORCES(ng)%evap(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idrain,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgrain(i,j)=AVERAGE(ng)%avgrain(i,j)+        &
#    ifdef WET_DRY
     &                                 GRID(ng)%rmask_full(i,j)*        &
#    endif
     &                                 FORCES(ng)%rain(i,j)
            END DO
          END DO
        END IF
#   endif
#  endif
# endif

# ifdef NEARSHORE_MELLOR
!
!  Accumulate radiation stresses.
!
        IF (Aout(idU2Sd,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgu2Sd(i,j)=AVERAGE(ng)%avgu2Sd(i,j)+        &
#  ifdef WET_DRY
     &                                 GRID(ng)%umask_full(i,j)*        &
#  endif
     &                                 OCEAN(ng)%ubar_stokes(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idV2Sd,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgv2Sd(i,j)=AVERAGE(ng)%avgv2Sd(i,j)+        &
#  ifdef WET_DRY
     &                                 GRID(ng)%vmask_full(i,j)*        &
#  endif
     &                                 OCEAN(ng)%vbar_stokes(i,j)
            END DO
          END DO
        END IF

        IF (Aout(idU2rs,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgu2RS(i,j)=AVERAGE(ng)%avgu2RS(i,j)+        &
#  ifdef WET_DRY
     &                                 GRID(ng)%umask_full(i,j)*        &
#  endif
     &                                 MIXING(ng)%rustr2d(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idV2rs,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgv2RS(i,j)=AVERAGE(ng)%avgv2RS(i,j)+        &
#  ifdef WET_DRY
     &                                 GRID(ng)%vmask_full(i,j)*        &
#  endif
     &                                 MIXING(ng)%rvstr2d(i,j)
            END DO
          END DO
        END IF

        IF (Aout(idW2xx,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgSxx2d(i,j)=AVERAGE(ng)%avgSxx2d(i,j)+      &
#  ifdef WET_DRY
     &                                  GRID(ng)%rmask_full(i,j)*       &
#  endif
     &                                  MIXING(ng)%Sxx_bar(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idW2xy,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgSxy2d(i,j)=AVERAGE(ng)%avgSxy2d(i,j)+      &
#  ifdef WET_DRY
     &                                  GRID(ng)%rmask_full(i,j)*       &
#  endif
     &                                  MIXING(ng)%Sxy_bar(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idW2yy,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgSyy2d(i,j)=AVERAGE(ng)%avgSyy2d(i,j)+      &
#  ifdef WET_DRY
     &                                  GRID(ng)%rmask_full(i,j)*       &
#  endif
     &                                  MIXING(ng)%Syy_bar(i,j)
            END DO
          END DO
        END IF

#  ifdef SOLVE3D
        IF (Aout(idU3Sd,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgu3Sd(i,j,k)=AVERAGE(ng)%avgu3Sd(i,j,k)+  &
#   ifdef WET_DRY
     &                                     GRID(ng)%umask_full(i,j)*    &
#   endif
     &                                     OCEAN(ng)%u_stokes(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idV3Sd,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgv3Sd(i,j,k)=AVERAGE(ng)%avgv3Sd(i,j,k)+  &
#   ifdef WET_DRY
     &                                     GRID(ng)%vmask_full(i,j)*    &
#   endif
     &                                     OCEAN(ng)%v_stokes(i,j,k)
              END DO
            END DO
          END DO
        END IF

        IF (Aout(idU3rs,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgu3RS(i,j,k)=AVERAGE(ng)%avgu3RS(i,j,k)+  &
#   ifdef WET_DRY
     &                                     GRID(ng)%umask_full(i,j)*    &
#   endif
     &                                     MIXING(ng)%rustr3d(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idV3rs,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgv3RS(i,j,k)=AVERAGE(ng)%avgv3RS(i,j,k)+  &
#   ifdef WET_DRY
     &                                     GRID(ng)%vmask_full(i,j)*    &
#   endif
     &                                     MIXING(ng)%rvstr3d(i,j,k)
              END DO
            END DO
          END DO
        END IF

        IF (Aout(idW3xx,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgSxx3d(i,j,k)=AVERAGE(ng)%avgSxx3d(i,j,k)+&
#   ifdef WET_DRY
     &                                      GRID(ng)%rmask_full(i,j)*   &
#   endif
     &                                      MIXING(ng)%Sxx(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idW3xy,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgSxy3d(i,j,k)=AVERAGE(ng)%avgSxy3d(i,j,k)+&
#   ifdef WET_DRY
     &                                      GRID(ng)%rmask_full(i,j)*   &
#   endif
     &                                      MIXING(ng)%Sxy(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idW3yy,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgSyy3d(i,j,k)=AVERAGE(ng)%avgSyy3d(i,j,k)+&
#   ifdef WET_DRY
     &                                      GRID(ng)%rmask_full(i,j)*   &
#   endif
     &                                      MIXING(ng)%Syy(i,j,k)
              END DO
            END DO
          END DO
        END IF
#   ifdef NEARSHORE_MELLOR05
        IF (Aout(idW3zx,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgSzx3d(i,j,k)=AVERAGE(ng)%avgSzx3d(i,j,k)+&
#    ifdef WET_DRY
     &                                      GRID(ng)%rmask_full(i,j)*   &
#    endif
     &                                      MIXING(ng)%Szx(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idW3zy,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgSzy3d(i,j,k)=AVERAGE(ng)%avgSzy3d(i,j,k)+&
#    ifdef WET_DRY
     &                                      GRID(ng)%rmask_full(i,j)*   &
#    endif
     &                                      MIXING(ng)%Szy(i,j,k)
              END DO
            END DO
          END DO
        END IF
#   endif
#  endif
# endif
!
!  Accumulate vorticity fields.
!
        IF (Aout(id2dPV,ng)) THEN
          DO j=Jstr,Jend
            DO i=Istr,Iend
              AVERAGE(ng)%avgpvor2d(i,j)=AVERAGE(ng)%avgpvor2d(i,j)+    &
# ifdef WET_DRY
     &                                   GRID(ng)%pmask_full(i,j)*      &
# endif
     &                                   potvor_bar(i,j)
            END DO
          END DO
        END IF
        IF (Aout(id2dRV,ng)) THEN
          DO j=Jstr,Jend
            DO i=Istr,Iend
              AVERAGE(ng)%avgrvor2d(i,j)=AVERAGE(ng)%avgrvor2d(i,j)+    &
# ifdef WET_DRY
     &                                   GRID(ng)%pmask_full(i,j)*      &
# endif
     &                                   relvor_bar(i,j)
            END DO
          END DO
        END IF
# ifdef SOLVE3D
        IF (Aout(id3dPV,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              DO i=Istr,Iend
                AVERAGE(ng)%avgpvor3d(i,j,k)=AVERAGE(ng)%avgpvor3d(i,j, &
     &                                                             k)+  &
#  ifdef WET_DRY
     &                                       GRID(ng)%pmask_full(i,j)*  &
#  endif
     &                                       potvor(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(id3dRV,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              DO i=Istr,Iend
                AVERAGE(ng)%avgrvor3d(i,j,k)=AVERAGE(ng)%avgrvor3d(i,j, &
     &                                                             k)+  &
#  ifdef WET_DRY
     &                                       GRID(ng)%pmask_full(i,j)*  &
#  endif
     &                                       relvor(i,j,k)
              END DO
            END DO
          END DO
        END IF
# endif
!
!  Accumulate quadratic fields.
!
        IF (Aout(idZZav,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgZZ(i,j)=AVERAGE(ng)%avgZZ(i,j)+            &
# ifdef WET_DRY
     &                               GRID(ng)%rmask_full(i,j)*          &
# endif
     &                               OCEAN(ng)%zeta(i,j,Kout)*          &
     &                               OCEAN(ng)%zeta(i,j,Kout)
            END DO
          END DO
        END IF
        IF (Aout(idU2av,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgU2(i,j)=AVERAGE(ng)%avgU2(i,j)+            &
# ifdef WET_DRY
     &                               GRID(ng)%umask_full(i,j)*          &
# endif
     &                               OCEAN(ng)%ubar(i,j,Kout)*          &
     &                               OCEAN(ng)%ubar(i,j,Kout)
            END DO
          END DO
        END IF
        IF (Aout(idV2av,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgV2(i,j)=AVERAGE(ng)%avgV2(i,j)+            &
# ifdef WET_DRY
     &                               GRID(ng)%vmask_full(i,j)*          &
# endif
     &                               OCEAN(ng)%vbar(i,j,Kout)*          &
     &                               OCEAN(ng)%vbar(i,j,Kout)
            END DO
          END DO
        END IF

# ifdef SOLVE3D
        IF (Aout(idUUav,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgUU(i,j,k)=AVERAGE(ng)%avgUU(i,j,k)+      &
#  ifdef WET_DRY
     &                                   GRID(ng)%umask_full(i,j)*      &
#  endif
     &                                   OCEAN(ng)%u(i,j,k,Nout)*       &
     &                                   OCEAN(ng)%u(i,j,k,Nout)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idVVav,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgVV(i,j,k)=AVERAGE(ng)%avgVV(i,j,k)+      &
#  ifdef WET_DRY
     &                                   GRID(ng)%vmask_full(i,j)*      &
#  endif
     &                                   OCEAN(ng)%v(i,j,k,Nout)*       &
     &                                   OCEAN(ng)%v(i,j,k,Nout)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idUVav,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              DO i=Istr,Iend
                AVERAGE(ng)%avgUV(i,j,k)=AVERAGE(ng)%avgUV(i,j,k)+      &
#  ifdef WET_DRY
     &                                   GRID(ng)%rmask_full(i,j)*      &
#  endif
     &                                   0.25_r8*                       &
     &                                   (OCEAN(ng)%u(i  ,j  ,k,Nout)+  &
     &                                    OCEAN(ng)%u(i+1,j  ,k,Nout))* &
     &                                   (OCEAN(ng)%v(i  ,j  ,k,Nout)+  &
     &                                    OCEAN(ng)%v(i  ,j+1,k,Nout))
              END DO
            END DO
          END DO
        END IF

        IF (Aout(idHUav,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgHuon(i,j,k)=AVERAGE(ng)%avgHuon(i,j,k)+  &
#  ifdef WET_DRY
     &                                     GRID(ng)%umask_full(i,j)*    &
#  endif
     &                                     GRID(ng)%Huon(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idHVav,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgHvom(i,j,k)=AVERAGE(ng)%avgHvom(i,j,k)+  &
#  ifdef WET_DRY
     &                                     GRID(ng)%vmask_full(i,j)*    &
#  endif
     &                                     GRID(ng)%Hvom(i,j,k)
              END DO
            END DO
          END DO
        END IF

        DO it=1,NT(ng)
          IF (Aout(idTTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgTT(i,j,k,it)=AVERAGE(ng)%avgTT(i,j,k,  &
     &                                                          it)+    &
#  ifdef WET_DRY
     &                                        GRID(ng)%rmask_full(i,j)* &
#  endif
     &                                        OCEAN(ng)%t(i,j,k,        &
     &                                                    Nout,it)*     &
     &                                        OCEAN(ng)%t(i,j,k,        &
     &                                                    Nout,it)
                END DO
              END DO
            END DO
          END IF
          IF (Aout(idUTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=Istr,Iend
                  AVERAGE(ng)%avgUT(i,j,k,it)=AVERAGE(ng)%avgUT(i,j,k,  &
     &                                                          it)+    &
#  ifdef WET_DRY
     &                                        GRID(ng)%umask_full(i,j)* &
#  endif
     &                                        0.5_r8*                   &
     &                                        OCEAN(ng)%u(i,j,k,Nout)*  &
     &                                        (OCEAN(ng)%t(i-1,j,k,     &
     &                                                     Nout,it)+    &
     &                                         OCEAN(ng)%t(i  ,j,k,     &
     &                                                     Nout,it))
                END DO
              END DO
            END DO
          END IF
          IF (Aout(idVTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=Jstr,Jend
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgVT(i,j,k,it)=AVERAGE(ng)%avgVT(i,j,k,  &
     &                                                          it)+    &
#  ifdef WET_DRY
     &                                        GRID(ng)%vmask_full(i,j)* &
#  endif
     &                                        0.5_r8*                   &
     &                                        OCEAN(ng)%v(i,j,k,Nout)*  &
     &                                        (OCEAN(ng)%t(i,j-1,k,     &
     &                                                     Nout,it)+    &
     &                                         OCEAN(ng)%t(i,j  ,k,     &
     &                                                     Nout,it))
                END DO
              END DO
            END DO
          END IF

          IF (Aout(iHUTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=Istr,Iend
                  AVERAGE(ng)%avgHuonT(i,j,k,it)=AVERAGE(ng)%avgHuonT(i,&
     &                                                       j,k,it)+   &
#  ifdef WET_DRY
     &                                           GRID(ng)%umask_full(i, &
     &                                                            j)*   &
#  endif
     &                                           0.5_r8*                &
     &                                           GRID(ng)%Huon(i,j,k)*  &
     &                                           (OCEAN(ng)%t(i-1,j,k,  &
     &                                                        Nout,it)+ &
     &                                            OCEAN(ng)%t(i  ,j,k,  &
     &                                                        Nout,it))
                END DO
              END DO
            END DO
          END IF
          IF (Aout(iHVTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=Jstr,Jend
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgHvomT(i,j,k,it)=AVERAGE(ng)%avgHvomT(i,&
     &                                                       j,k,it)+   &
#  ifdef WET_DRY
     &                                           GRID(ng)%vmask_full(i, &
     &                                                            j)*   &
#  endif
     &                                           0.5_r8*                &
     &                                           GRID(ng)%Hvom(i,j,k)*  &
     &                                           (OCEAN(ng)%t(i,j-1,k,  &
     &                                                        Nout,it)+ &
     &                                            OCEAN(ng)%t(i,j  ,k,  &
     &                                                        Nout,it))
                END DO
              END DO
            END DO
          END IF
        END DO
# endif
      END IF
!
!-----------------------------------------------------------------------
!  Convert accumulated sums into time-averages, if appropriate.
!  Notice that we need to apply periodic conditions, if any, since
!  the full I- and J-ranges are different.
!-----------------------------------------------------------------------
!
      IF (((iic(ng).gt.ntsAVG(ng)).and.                                 &
     &     (MOD(iic(ng)-1,nAVG(ng)).eq.0).and.                          &
     &     ((iic(ng).ne.ntstart(ng)).or.(nrrec(ng).eq.0))).or.          &
     &    ((iic(ng).ge.ntsAVG(ng)).and.(nAVG(ng).eq.1))) THEN
        IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
          IF (nAVG(ng).eq.1) THEN
            AVGtime(ng)=time(ng)
          ELSE
            AVGtime(ng)=AVGtime(ng)+REAL(nAVG(ng),r8)*dt(ng)
          END IF
        END IF
!
!  Set time-averaged factors for each C-grid variable type. Notice that
!  the I- and J-ranges are all grid types are the same for convinience.
# ifdef WET_DRY
!  In wetting and drying, the sums are devided by the number of times
!  that each qrid point is wet.
# endif
!
# ifdef WET_DRY
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            pfac(i,j)=1.0_r8/MAX(1.0_r8, GRID(ng)%pmask_avg(i,j))
            rfac(i,j)=1.0_r8/MAX(1.0_r8, GRID(ng)%rmask_avg(i,j))
            ufac(i,j)=1.0_r8/MAX(1.0_r8, GRID(ng)%umask_avg(i,j))
            vfac(i,j)=1.0_r8/MAX(1.0_r8, GRID(ng)%vmask_avg(i,j))
          END DO
        END DO
# else
        fac=1.0_r8/REAL(nAVG(ng),r8)
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            pfac(i,j)=fac
            rfac(i,j)=fac
            ufac(i,j)=fac
            vfac(i,j)=fac
          END DO
        END DO
# endif
!
!  Process state variables.
!
        IF (Aout(idFsur,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgzeta(i,j)=rfac(i,j)*                       &
     &                                 AVERAGE(ng)%avgzeta(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgzeta)
# ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgzeta)
# endif
          END IF
        END IF

        IF (Aout(idUbar,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgu2d(i,j)=ufac(i,j)*                        &
     &                                AVERAGE(ng)%avgu2d(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_u2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgu2d)
# ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgu2d)
# endif
          END IF
        END IF

        IF (Aout(idVbar,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgv2d(i,j)=vfac(i,j)*                        &
     &                                AVERAGE(ng)%avgv2d(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_v2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgv2d)
# ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgv2d)
# endif
          END IF
        END IF

        IF (Aout(idu2dE,ng).and.Aout(idv2dN,ng)) THEN
          DO j=Jstr,Jend
            DO i=Istr,Iend
              AVERAGE(ng)%avgu2dE(i,j)=rfac(i,j)*                       &
     &                                 AVERAGE(ng)%avgu2dE(i,j)
              AVERAGE(ng)%avgv2dN(i,j)=rfac(i,j)*                       &
     &                                 AVERAGE(ng)%avgv2dN(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgu2dE)
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgv2dN)
# ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 2,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgu2dE,                    &
     &                          AVERAGE(ng)%avgv2dN)
# endif
          END IF
        END IF

# ifdef SOLVE3D
        IF (Aout(idUvel,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgu3d(i,j,k)=ufac(i,j)*                    &
     &                                    AVERAGE(ng)%avgu3d(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_u3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              AVERAGE(ng)%avgu3d)
#  ifdef DISTRIBUTE
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgu3d)
#  endif
          END IF
        END IF

        IF (Aout(idVvel,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgv3d(i,j,k)=vfac(i,j)*                    &
     &                                    AVERAGE(ng)%avgv3d(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_v3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              AVERAGE(ng)%avgv3d)
#  ifdef DISTRIBUTE
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgv3d)
#  endif
          END IF
        END IF

        IF (Aout(idu3dE,ng).and.Aout(idv3dN,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              DO i=Istr,Iend
                AVERAGE(ng)%avgu3dE(i,j,k)=rfac(i,j)*                   &
     &                                     AVERAGE(ng)%avgu3dE(i,j,k)
                AVERAGE(ng)%avgv3dN(i,j,k)=rfac(i,j)*                   &
     &                                     AVERAGE(ng)%avgv3dN(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              AVERAGE(ng)%avgu3dE)
            CALL exchange_r3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              AVERAGE(ng)%avgv3dN)
#  ifdef DISTRIBUTE
            CALL mp_exchange3d (ng, tile, iNLM, 2,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgu3dE,                    &
     &                          AVERAGE(ng)%avgv3dN)
#  endif
          END IF
        END IF

        IF (Aout(idOvel,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgw3d(i,j,k)=rfac(i,j)*                    &
     &                                    AVERAGE(ng)%avgw3d(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_w3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 0, N(ng),       &
     &                              AVERAGE(ng)%avgw3d)
#  ifdef DISTRIBUTE
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 0, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgw3d)
#  endif
          END IF
        END IF

        IF (Aout(idWvel,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgwvel(i,j,k)=rfac(i,j)*                   &
     &                                     AVERAGE(ng)%avgwvel(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_w3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 0, N(ng),       &
     &                              AVERAGE(ng)%avgwvel)
#  ifdef DISTRIBUTE
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 0, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgwvel)
#  endif
          END IF
        END IF

        IF (Aout(idDano,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgrho(i,j,k)=rfac(i,j)*                    &
     &                                    AVERAGE(ng)%avgrho(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              AVERAGE(ng)%avgrho)
#  ifdef DISTRIBUTE
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgrho)
#  endif
          END IF
        END IF

        DO it=1,NT(ng)
          IF (Aout(idTvar(it),ng)) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgt(i,j,k,it)=rfac(i,j)*                 &
     &                                       AVERAGE(ng)%avgt(i,j,k,it)
                END DO
              END DO
            END DO
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_r3d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj, 1, N(ng),     &
     &                                AVERAGE(ng)%avgt(:,:,:,it))
#  ifdef DISTRIBUTE
              CALL mp_exchange3d (ng, tile, iNLM, 1,                    &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            NghostPoints,                         &
     &                            EWperiodic(ng), NSperiodic(ng),       &
     &                            AVERAGE(ng)%avgt(:,:,:,it))
#  endif
            END IF
          END IF
        END DO

#  if defined SEDIMENT && defined BEDLOAD
        DO it=1,NST
          IF (Aout(idUbld(it),ng)) THEN
            DO j=JstrR,JendR
              DO i=Istr,IendR
                SEDBED(ng)%avgbedldu(i,j,it)=ufac(i,j)*                 &
     &                                       SEDBED(ng)%avgbedldu(i,j,  &
     &                                                            it)
              END DO
            END DO
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_u2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                SEDBED(ng)%avgbedldu(:,:,it))
#   ifdef DISTRIBUTE
              CALL mp_exchange2d (ng, tile, iNLM, 1,                    &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            NghostPoints,                         &
     &                            EWperiodic(ng), NSperiodic(ng),       &
     &                            SEDBED(ng)%avgbedldu(:,:,it))
#   endif
            END IF
          END IF

          IF (Aout(idVbld(it),ng)) THEN
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                SEDBED(ng)%avgbedldv(i,j,it)=vfac(i,j)*                 &
     &                                       SEDBED(ng)%avgbedldv(i,j,  &
     &                                                            it)
              END DO
            END DO
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_v2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                SEDBED(ng)%avgbedldv(:,:,it))
#   ifdef DISTRIBUTE
              CALL mp_exchange2d (ng, tile, iNLM, 1,                    &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            NghostPoints,                         &
     &                            EWperiodic(ng), NSperiodic(ng),       &
     &                            SEDBED(ng)%avgbedldv(:,:,it))
#   endif
            END IF
          END IF
        END DO
#  endif

#  if defined LMD_MIXING || defined MY25_MIXING || defined GLS_MIXING
        IF (Aout(idVvis,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgAKv(i,j,k)=rfac(i,j)*                    &
     &                                    AVERAGE(ng)%avgAKv(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_w3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 0, N(ng),       &
     &                              AVERAGE(ng)%avgAKv)
#   ifdef DISTRIBUTE
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 0, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgAKv)
#   endif
          END IF
        END IF

        IF (Aout(idTdif,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgAKt(i,j,k)=rfac(i,j)*                    &
     &                                    AVERAGE(ng)%avgAKt(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_w3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 0, N(ng),       &
     &                              AVERAGE(ng)%avgAKt)
#   ifdef DISTRIBUTE
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 0, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgAKt)
#   endif
          END IF
        END IF

        IF (Aout(idSdif,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgAKs(i,j,k)=rfac(i,j)*                    &
     &                                    AVERAGE(ng)%avgAKs(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_w3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 0, N(ng),       &
     &                              AVERAGE(ng)%avgAKs)
#   ifdef DISTRIBUTE
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 0, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgAKs)
#   endif
          END IF
        END IF
#  endif

#  ifdef LMD_SKPP
        IF (Aout(idHsbl,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avghsbl(i,j)=rfac(i,j)*                       &
     &                                 AVERAGE(ng)%avghsbl(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avghsbl)
#   ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avghsbl)
#   endif
          END IF
        END IF
#  endif

#  ifdef LMD_BKPP
        IF (Aout(idHbbl,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avghbbl(i,j)=rfac(i,j)*                       &
     &                                 AVERAGE(ng)%avghbbl(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avghbbl)
#   ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avghbbl)
#   endif
          END IF
        END IF
#  endif
# endif

# if defined FORWARD_WRITE && defined SOLVE3D
!
!  Process 2D/3D coupling terms.
!
        IF (Aout(idUfx1,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgDU_avg1(i,j)=ufac(i,j)*                    &
     &                                    AVERAGE(ng)%avgDU_avg1(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_u2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgDU_avg1)
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgDU_avg1)
#  endif
          END IF
        END IF

        IF (Aout(idUfx2,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgDU_avg2(i,j)=ufac(i,j)*                    &
     &                                    AVERAGE(ng)%avgDU_avg2(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_u2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgDU_avg2)
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgDU_avg2)
#  endif
          END IF
        END IF

        IF (Aout(idVfx1,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgDV_avg1(i,j)=vfac(i,j)*                    &
     &                                    AVERAGE(ng)%avgDV_avg1(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_v2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgDV_avg1)
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgDV_avg1)
#  endif
          END IF
        END IF

        IF (Aout(idVfx2,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgDV_avg2(i,j)=vfac(i,j)*                    &
     &                                    AVERAGE(ng)%avgDV_avg2(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_v2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgDV_avg2)
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgDV_avg2)
#  endif
          END IF
        END IF
# endif
!
!  Process surface and bottom fluxes.
!
        IF (Aout(idUsms,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgsus(i,j)=ufac(i,j)*                        &
     &                                AVERAGE(ng)%avgsus(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_u2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgsus)
# ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgsus)
# endif
          END IF
        END IF

        IF (Aout(idVsms,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgsvs(i,j)=vfac(i,j)*                        &
     &                                AVERAGE(ng)%avgsvs(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_v2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgsvs)
# ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgsvs)
# endif
          END IF
        END IF

        IF (Aout(idUbms,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgbus(i,j)=ufac(i,j)*                        &
     &                                AVERAGE(ng)%avgbus(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_u2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgbus)
# ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgbus)
# endif
          END IF
        END IF

        IF (Aout(idVbms,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgbvs(i,j)=vfac(i,j)*                        &
     &                                AVERAGE(ng)%avgbvs(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_v2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgbvs)
# ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgbvs)
# endif
          END IF
        END IF

# ifdef SOLVE3D
#  if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
        IF (Aout(idPair,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgPair(i,j)=rfac(i,j)*                       &
     &                                 AVERAGE(ng)%avgPair(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgPair)
#   ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgPair)
#   endif
          END IF
        END IF
#  endif

#  if defined BULK_FLUXES || defined ECOSIM
        IF (Aout(idUair,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgUwind(i,j)=rfac(i,j)*                      &
     &                                  AVERAGE(ng)%avgUwind(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgUwind)
#   ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgUwind)
#   endif
          END IF
        END IF

        IF (Aout(idVair,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgVwind(i,j)=rfac(i,j)*                      &
     &                                  AVERAGE(ng)%avgVwind(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgVwind)
#   ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgVwind)
#   endif
          END IF
        END IF
#  endif

        IF (Aout(idTsur(itemp),ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgstf(i,j)=rfac(i,j)*                        &
     &                                AVERAGE(ng)%avgstf(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgstf)
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgstf)
#  endif
          END IF
        END IF

        IF (Aout(idTsur(isalt),ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgswf(i,j)=rfac(i,j)*                        &
     &                                AVERAGE(ng)%avgswf(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgswf)
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgswf)
#  endif
          END IF
        END IF

#  ifdef SHORTWAVE
        IF (Aout(idSrad,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgsrf(i,j)=rfac(i,j)*                        &
     &                                AVERAGE(ng)%avgsrf(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgsrf)
#   ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgsrf)
#   endif
          END IF
        END IF
#  endif

#  ifdef BULK_FLUXES
        IF (Aout(idLhea,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avglhf(i,j)=rfac(i,j)*                        &
     &                                AVERAGE(ng)%avglhf(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avglhf)
#   ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avglhf)
#   endif
          END IF
        END IF

        IF (Aout(idShea,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgshf(i,j)=rfac(i,j)*                        &
     &                                AVERAGE(ng)%avgshf(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgshf)
#   ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgshf)
#   endif
          END IF
        END IF

        IF (Aout(idLrad,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avglrf(i,j)=rfac(i,j)*                        &
     &                                AVERAGE(ng)%avglrf(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avglrf)
#   ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avglrf)
#   endif
          END IF
        END IF

#   ifdef EMINUSP
        IF (Aout(idevap,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgevap(i,j)=rfac(i,j)*                       &
     &                                 AVERAGE(ng)%avgevap(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgevap)
#    ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgevap)
#    endif
          END IF
        END IF

        IF (Aout(idrain,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgrain(i,j)=rfac(i,j)*                       &
     &                                 AVERAGE(ng)%avgrain(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgrain)
#    ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgrain)
#    endif
          END IF
        END IF
#   endif
#  endif
# endif

# ifdef NEARSHORE_MELLOR
!
!  Process radiation stresses.
!
        IF (Aout(idU2Sd,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgu2Sd(i,j)=ufac(i,j)*                       &
     &                                 AVERAGE(ng)%avgu2Sd(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_u2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgu2Sd)
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgu2Sd)
#  endif
          END IF
        END IF

        IF (Aout(idV2Sd,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgv2Sd(i,j)=vfac(i,j)*                       &
     &                                 AVERAGE(ng)%avgv2Sd(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_v2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgv2Sd)
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgv2Sd)
#  endif
          END IF
        END IF

        IF (Aout(idU2rs,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgu2RS(i,j)=ufac(i,j)*                       &
     &                                 AVERAGE(ng)%avgu2RS(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_u2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgu2RS)
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgu2RS)
#  endif
          END IF
        END IF

        IF (Aout(idV2rs,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgv2RS(i,j)=vfac(i,j)*                       &
     &                                 AVERAGE(ng)%avgv2RS(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_v2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgv2RS)
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgv2RS)
#  endif
          END IF
        END IF

        IF (Aout(idW2xx,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgSxx2d(i,j)=rfac(i,j)*                      &
     &                                  AVERAGE(ng)%avgSxx2d(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgSxx2d)
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgSxx2d)
#  endif
          END IF
        END IF

        IF (Aout(idW2xy,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgSxy2d(i,j)=rfac(i,j)*                      &
     &                                  AVERAGE(ng)%avgSxy2d(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgSxy2d)
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgSxy2d)
#  endif
          END IF
        END IF

        IF (Aout(idW2yy,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgSyy2d(i,j)=rfac(i,j)*                      &
     &                                  AVERAGE(ng)%avgSyy2d(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgSyy2d)
#  ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgSyy2d)
#  endif
          END IF
        END IF

#  ifdef SOLVE3D
        IF (Aout(idU3Sd,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgu3Sd(i,j,k)=ufac(i,j)*                   &
     &                                     AVERAGE(ng)%avgu3Sd(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_u3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              AVERAGE(ng)%avgu3Sd)
#   ifdef DISTRIBUTE
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgu3Sd)
#   endif
          END IF
        END IF

        IF (Aout(idV3Sd,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgv3Sd(i,j,k)=vfac(i,j)*                   &
     &                                     AVERAGE(ng)%avgv3Sd(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_v3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              AVERAGE(ng)%avgv3Sd)
#   ifdef DISTRIBUTE
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgv3Sd)
#   endif
          END IF
        END IF

        IF (Aout(idU3rs,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgu3RS(i,j,k)=ufac(i,j)*                   &
     &                                     AVERAGE(ng)%avgu3RS(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_u3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              AVERAGE(ng)%avgu3RS)
#   ifdef DISTRIBUTE
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgu3RS)
#   endif
          END IF
        END IF

        IF (Aout(idV3rs,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgv3RS(i,j,k)=vfac(i,j)*                   &
     &                                     AVERAGE(ng)%avgv3RS(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_v3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              AVERAGE(ng)%avgv3RS)
#   ifdef DISTRIBUTE
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgv3RS)
#   endif
          END IF
        END IF

        IF (Aout(idW3xx,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgSxx3d(i,j,k)=rfac(i,j)*                  &
     &                                      AVERAGE(ng)%avgSxx3d(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              AVERAGE(ng)%avgSxx3d)
#   ifdef DISTRIBUTE
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgSxx3d)
#   endif
          END IF
        END IF

        IF (Aout(idW3xy,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgSxy3d(i,j,k)=rfac(i,j)*                  &
     &                                      AVERAGE(ng)%avgSxy3d(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              AVERAGE(ng)%avgSxy3d)
#   ifdef DISTRIBUTE
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgSxy3d)
#   endif
          END IF
        END IF

        IF (Aout(idW3yy,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgSyy3d(i,j,k)=rfac(i,j)*                  &
     &                                      AVERAGE(ng)%avgSyy3d(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              AVERAGE(ng)%avgSyy3d)
#   ifdef DISTRIBUTE
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgSyy3d)
#   endif
          END IF
        END IF

#   ifdef NEARSHORE_MELLOR05
        IF (Aout(idW3zx,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgSzx3d(i,j,k)=rfac(i,j)*                  &
     &                                      AVERAGE(ng)%avgSzx3d(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              AVERAGE(ng)%avgSzx3d)
#    ifdef DISTRIBUTE
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgSzx3d)
#    endif
          END IF
        END IF

        IF (Aout(idW3zy,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgSzy3d(i,j,k)=rfac(i,j)*                  &
     &                                      AVERAGE(ng)%avgSzy3d(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              AVERAGE(ng)%avgSzy3d)
#    ifdef DISTRIBUTE
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgSzy3d)
#    endif
          END IF
        END IF
#   endif
#  endif
# endif
!
!  Process vorticity fields.
!
        IF (Aout(id2dPV,ng)) THEN
          DO j=Jstr,Jend
            DO i=Istr,Iend
              AVERAGE(ng)%avgpvor2d(i,j)=pfac(i,j)*                     &
     &                                   AVERAGE(ng)%avgpvor2d(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_p2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgpvor2d)
# ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgpvor2d)
# endif
          END IF
        END IF

        IF (Aout(id2dRV,ng)) THEN
          DO j=Jstr,Jend
            DO i=Istr,Iend
              AVERAGE(ng)%avgrvor2d(i,j)=pfac(i,j)*                     &
     &                                   AVERAGE(ng)%avgrvor2d(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_p2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgrvor2d)
# ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgrvor2d)
# endif
          END IF
        END IF

# ifdef SOLVE3D
        IF (Aout(id3dPV,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              DO i=Istr,Iend
                AVERAGE(ng)%avgpvor3d(i,j,k)=pfac(i,j)*                 &
     &                                      AVERAGE(ng)%avgpvor3d(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_p3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              AVERAGE(ng)%avgpvor3d)
#  ifdef DISTRIBUTE
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgpvor3d)
#  endif
          END IF
        END IF

        IF (Aout(id3dRV,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              DO i=Istr,Iend
                AVERAGE(ng)%avgrvor3d(i,j,k)=pfac(i,j)*                 &
     &                                      AVERAGE(ng)%avgrvor3d(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_p3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              AVERAGE(ng)%avgrvor3d)
#  ifdef DISTRIBUTE
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgrvor3d)
#  endif
          END IF
        END IF
# endif
!
!  Process quadratic fields.
!
        IF (Aout(idZZav,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgZZ(i,j)=rfac(i,j)*                         &
     &                               AVERAGE(ng)%avgZZ(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgZZ)
# ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgZZ)
# endif
          END IF
        END IF

        IF (Aout(idU2av,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgU2(i,j)=ufac(i,j)*                         &
     &                               AVERAGE(ng)%avgU2(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_u2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgU2)
# ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgU2)
# endif
          END IF
        END IF

        IF (Aout(idV2av,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgV2(i,j)=vfac(i,j)*                         &
     &                               AVERAGE(ng)%avgV2(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_v2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgV2)
# ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgV2)
# endif
          END IF
        END IF

# ifdef SOLVE3D
        IF (Aout(idUUav,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgUU(i,j,k)=ufac(i,j)*                     &
     &                                   AVERAGE(ng)%avgUU(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_u3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              AVERAGE(ng)%avgUU)
#  ifdef DISTRIBUTE
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgUU)
#  endif
          END IF
        END IF

        IF (Aout(idVVav,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgVV(i,j,k)=vfac(i,j)*                     &
     &                                   AVERAGE(ng)%avgVV(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_v3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              AVERAGE(ng)%avgVV)
#  ifdef DISTRIBUTE
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgVV)
#  endif
          END IF
        END IF

        IF (Aout(idUVav,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              DO i=Istr,Iend
                AVERAGE(ng)%avgUV(i,j,k)=rfac(i,j)*                     &
     &                                   AVERAGE(ng)%avgUV(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              AVERAGE(ng)%avgUV)
#  ifdef DISTRIBUTE
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgUV)
#  endif
          END IF
        END IF

        IF (Aout(idHUav,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgHuon(i,j,k)=ufac(i,j)*                   &
     &                                     AVERAGE(ng)%avgHuon(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_u3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              AVERAGE(ng)%avgHuon)
#  ifdef DISTRIBUTE
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgHuon)
#  endif
          END IF
        END IF

        IF (Aout(idHVav,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgHvom(i,j,k)=vfac(i,j)*                   &
     &                                     AVERAGE(ng)%avgHvom(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_v3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              AVERAGE(ng)%avgHvom)
#  ifdef DISTRIBUTE
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgHvom)
#  endif
          END IF
        END IF

        DO it=1,NT(ng)
          IF (Aout(idTTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgTT(i,j,k,it)=rfac(i,j)*                &
     &                                       AVERAGE(ng)%avgTT(i,j,k,it)
                END DO
              END DO
            END DO
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_r3d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj, 1, N(ng),     &
     &                                AVERAGE(ng)%avgTT(:,:,:,it))
#  ifdef DISTRIBUTE
              CALL mp_exchange3d (ng, tile, iNLM, 1,                    &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            NghostPoints,                         &
     &                            EWperiodic(ng), NSperiodic(ng),       &
     &                            AVERAGE(ng)%avgTT(:,:,:,it))
#  endif
            END IF
          END IF

          IF (Aout(idUTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=Istr,Iend
                  AVERAGE(ng)%avgUT(i,j,k,it)=ufac(i,j)*                &
     &                                       AVERAGE(ng)%avgUT(i,j,k,it)
                END DO
              END DO
            END DO
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_u3d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj, 1, N(ng),     &
     &                                AVERAGE(ng)%avgUT(:,:,:,it))
#  ifdef DISTRIBUTE
              CALL mp_exchange3d (ng, tile, iNLM, 1,                    &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            NghostPoints,                         &
     &                            EWperiodic(ng), NSperiodic(ng),       &
     &                            AVERAGE(ng)%avgUT(:,:,:,it))
#  endif
            END IF
          END IF

          IF (Aout(idVTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=Jstr,Jend
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgVT(i,j,k,it)=vfac(i,j)*                &
     &                                       AVERAGE(ng)%avgVT(i,j,k,it)
                END DO
              END DO
            END DO
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_v3d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj, 1, N(ng),     &
     &                                AVERAGE(ng)%avgVT(:,:,:,it))
#  ifdef DISTRIBUTE
              CALL mp_exchange3d (ng, tile, iNLM, 1,                    &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            NghostPoints,                         &
     &                            EWperiodic(ng), NSperiodic(ng),       &
     &                            AVERAGE(ng)%avgVT(:,:,:,it))
#  endif
            END IF
          END IF

          IF (Aout(iHUTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=Istr,Iend
                  AVERAGE(ng)%avgHuonT(i,j,k,it)=ufac(i,j)*             &
     &                                    AVERAGE(ng)%avgHuonT(i,j,k,it)
                END DO
              END DO
            END DO
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_u3d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj, 1, N(ng),     &
     &                                AVERAGE(ng)%avgHuonT(:,:,:,it))
#  ifdef DISTRIBUTE
              CALL mp_exchange3d (ng, tile, iNLM, 1,                    &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            NghostPoints,                         &
     &                            EWperiodic(ng), NSperiodic(ng),       &
     &                            AVERAGE(ng)%avgHuonT(:,:,:,it))
#  endif
            END IF
          END IF

          IF (Aout(iHVTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=Jstr,Jend
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgHvomT(i,j,k,it)=vfac(i,j)*             &
     &                                    AVERAGE(ng)%avgHvomT(i,j,k,it)
                END DO
              END DO
            END DO
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_v3d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj, 1, N(ng),     &
     &                                AVERAGE(ng)%avgHvomT(:,:,:,it))
#  ifdef DISTRIBUTE
              CALL mp_exchange3d (ng, tile, iNLM, 1,                    &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            NghostPoints,                         &
     &                            EWperiodic(ng), NSperiodic(ng),       &
     &                            AVERAGE(ng)%avgHvomT(:,:,:,it))
#  endif
            END IF
          END IF
        END DO
# endif
      END IF

      RETURN
      END SUBROUTINE set_avg_tile

# if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
!
!***********************************************************************
      SUBROUTINE set_detide_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            NTC, Kout,                            &
#  ifdef SOLVE3D
     &                            Nout,                                 &
#  endif
     &                            CosOmega, SinOmega,                   &
     &                            CosW_avg, CosW_sum,                   &
     &                            SinW_avg, SinW_sum,                   &
     &                            CosWCosW, SinWSinW, SinWCosW)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_average
#  ifdef WET_DRY
      USE mod_grid
#  endif
      USE mod_ocean
      USE mod_scalars
      USE mod_tides
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Kout
#  ifdef SOLVE3D
      integer, intent(in) :: Nout
#  endif
      integer, intent(in) :: NTC
!
#  ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: CosOmega(:)
      real(r8), intent(in) :: SinOmega(:)
      real(r8), intent(inout) :: CosW_avg(:)
      real(r8), intent(inout) :: CosW_sum(:)
      real(r8), intent(inout) :: SinW_avg(:)
      real(r8), intent(inout) :: SinW_sum(:)
      real(r8), intent(inout) :: CosWCosW(:,:)
      real(r8), intent(inout) :: SinWSinW(:,:)
      real(r8), intent(inout) :: SinWCosW(:,:)
#  else
      real(r8), intent(in) :: CosOmega(NTC)
      real(r8), intent(in) :: SinOmega(NTC)
      real(r8), intent(inout) :: CosW_avg(NTC)
      real(r8), intent(inout) :: CosW_sum(NTC)
      real(r8), intent(inout) :: SinW_avg(NTC)
      real(r8), intent(inout) :: SinW_sum(NTC)
      real(r8), intent(inout) :: CosWCosW(NTC,NTC)
      real(r8), intent(inout) :: SinWSinW(NTC,NTC)
      real(r8), intent(inout) :: SinWCosW(NTC,NTC)
#  endif
!
!  Local variable declarations.
!
      integer :: i, it, j, k
      integer :: NTC2, mk, nk

      integer, dimension(2*NTC+1) :: indx

      real(r8) :: fac, fac1
      real(r8) :: Hsum, d

      real(r8), dimension(0:2*NTC) :: Ak
      real(r8), dimension(0:2*NTC) :: tide_harmonics
      real(r8), dimension(0:2*NTC,0:2*NTC) :: C, Y

      real(r8) :: rfac(IminS:ImaxS,JminS:JmaxS)
      real(r8) :: ufac(IminS:ImaxS,JminS:JmaxS)
      real(r8) :: vfac(IminS:ImaxS,JminS:JmaxS)

#  include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Return if time-averaging window is zero.
!-----------------------------------------------------------------------
!
      IF (nAVG(ng).eq.0) RETURN
!
!-----------------------------------------------------------------------
!  Initialize time-averaged arrays when appropriate.  Notice that
!  fields are initilized twice during re-start.  However, the time-
!  averaged fields are computed correctly.
!-----------------------------------------------------------------------
!
      NTC2=2*NTC

      IF (((iic(ng).gt.ntsAVG(ng)).and.                                 &
     &     (MOD(iic(ng)-1,nAVG(ng)).eq.1)).or.                          &
     &    ((nrrec(ng).gt.0).and.(iic(ng).eq.ntstart(ng)))) THEN
!
!  Compute least-squares coefficients to detide time-averaged fields.
!  Notice that the coefficients are always accumulated and not
!  re-initialized.  This allows better tidal fit as the simulation
!  progresses.
!
        IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
          Hcount(ng)=Hcount(ng)+1
          DO nk=1,NTC
            SinW_avg(nk)=SinOmega(nk)
            CosW_avg(nk)=CosOmega(nk)
            SinW_sum(nk)=SinW_sum(nk)+SinOmega(nk)
            CosW_sum(nk)=CosW_sum(nk)+CosOmega(nk)
            DO mk=1,NTC
              SinWSinW(mk,nk)=SinWSinW(mk,nk)+SinOmega(mk)*SinOmega(nk)
              CosWCosW(mk,nk)=CosWCosW(mk,nk)+CosOmega(mk)*CosOmega(nk)
              SinWCosW(mk,nk)=SinWCosW(mk,nk)+SinOmega(mk)*CosOmega(nk)
            END DO
          END DO
          tide_harmonics(0)=1.0_r8
          DO nk=1,NTC
            tide_harmonics(nk    )=SinOmega(nk)
            tide_harmonics(nk+NTC)=CosOmega(nk)
          END DO
        END IF
!
!  Initialize.
!
        IF (Aout(idFsuD,ng)) THEN
          DO nk=0,NTC2
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                TIDES(ng)%zeta_tide(i,j,nk)=TIDES(ng)%zeta_tide(i,j,    &
     &                                                          nk)+    &
#  ifdef WET_DRY
     &                                      GRID(ng)%rmask_full(i,j)*   &
#  endif
     &                                      OCEAN(ng)%zeta(i,j,Kout)*   &
     &                                      tide_harmonics(nk)
              END DO
            END DO
          END DO
        END IF

        IF (Aout(idu2dD,ng)) THEN
          DO nk=0,NTC2
            DO j=JstrR,JendR
              DO i=Istr,IendR
                TIDES(ng)%ubar_tide(i,j,nk)=TIDES(ng)%ubar_tide(i,j,    &
     &                                                          nk)+    &
#  ifdef WET_DRY
     &                                      GRID(ng)%umask_full(i,j)*   &
#  endif
     &                                      OCEAN(ng)%ubar(i,j,Kout)*   &
     &                                      tide_harmonics(nk)
              END DO
            END DO
          END DO
        END IF

        IF (Aout(idv2dD,ng)) THEN
          DO nk=0,NTC2
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                TIDES(ng)%vbar_tide(i,j,nk)=TIDES(ng)%vbar_tide(i,j,    &
     &                                                          nk)+    &
#  ifdef WET_DRY
     &                                      GRID(ng)%vmask_full(i,j)*   &
#  endif
     &                                      OCEAN(ng)%vbar(i,j,Kout)*   &
     &                                      tide_harmonics(nk)
              END DO
            END DO
          END DO
        END IF

#  ifdef SOLVE3D
        IF (Aout(idu3dD,ng)) THEN
          DO nk=0,NTC2
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=Istr,IendR
                  TIDES(ng)%u_tide(i,j,k,nk)=TIDES(ng)%u_tide(i,j,k,    &
     &                                                        nk)+      &
#   ifdef WET_DRY
     &                                       GRID(ng)%umask_full(i,j)*  &
#   endif
     &                                       OCEAN(ng)%u(i,j,k,Nout)*   &
     &                                       tide_harmonics(nk)
                END DO
              END DO
            END DO
          END DO
        END IF

        IF (Aout(idv3dD,ng)) THEN
          DO nk=0,NTC2
            DO k=1,N(ng)
              DO j=Jstr,JendR
                DO i=IstrR,IendR
                  TIDES(ng)%v_tide(i,j,k,nk)=TIDES(ng)%v_tide(i,j,k,    &
     &                                                        nk)+      &
#   ifdef WET_DRY
     &                                       GRID(ng)%vmask_full(i,j)*  &
#   endif
     &                                       OCEAN(ng)%v(i,j,k,Nout)*   &
     &                                       tide_harmonics(nk)
                END DO
              END DO
            END DO
          END DO
        END IF

        DO it=1,NAT
          IF (Aout(idTrcD(it),ng)) THEN
            DO nk=0,NTC2
              DO k=1,N(ng)
                DO j=Jstr,JendR
                  DO i=IstrR,IendR
                    TIDES(ng)%t_tide(i,j,k,nk,it)=TIDES(ng)%t_tide(i,j, &
     &                                                        k,nk,it)+ &
#   ifdef WET_DRY
     &                                        GRID(ng)%rmask_full(i,j)* &
#   endif
     &                                            OCEAN(ng)%t(i,j,k,    &
     &                                                        Nout,it)* &
     &                                            tide_harmonics(nk)
                  END DO
                END DO
              END DO
            END DO
          END IF
        END DO
#  endif
!
!-----------------------------------------------------------------------
!  Accumulate time-averaged fields.
!-----------------------------------------------------------------------
!
      ELSE IF (iic(ng).gt.ntsAVG(ng)) THEN
!
!  Accumukate Detide least-squares coefficients. They only vary in time
!  since omega (as computed in set_tides) uses model time coordinate.
!
        IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
          Hcount(ng)=Hcount(ng)+1
          DO nk=1,NTC
            SinW_avg(nk)=SinW_avg(nk)+SinOmega(nk)
            CosW_avg(nk)=CosW_avg(nk)+CosOmega(nk)
            SinW_sum(nk)=SinW_sum(nk)+SinOmega(nk)
            CosW_sum(nk)=CosW_sum(nk)+CosOmega(nk)
            DO mk=1,NTC
              SinWSinW(mk,nk)=SinWSinW(mk,nk)+SinOmega(mk)*SinOmega(nk)
              CosWCosW(mk,nk)=CosWCosW(mk,nk)+CosOmega(mk)*CosOmega(nk)
              SinWCosW(mk,nk)=SinWCosW(mk,nk)+SinOmega(mk)*CosOmega(nk)
            END DO
          END DO
          tide_harmonics(0)=1.0_r8
          DO nk=1,NTC
            tide_harmonics(nk    )=SinOmega(nk)
            tide_harmonics(nk+NTC)=CosOmega(nk)
          END DO
        END IF
!
!  Accumulate.
!
        IF (Aout(idFsuD,ng)) THEN
          DO nk=0,NTC2
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                TIDES(ng)%zeta_tide(i,j,nk)=TIDES(ng)%zeta_tide(i,j,    &
     &                                                          nk)+    &
#  ifdef WET_DRY
     &                                      GRID(ng)%rmask_full(i,j)*   &
#  endif
     &                                      OCEAN(ng)%zeta(i,j,Kout)*   &
     &                                      tide_harmonics(nk)
              END DO
            END DO
          END DO
        END IF

        IF (Aout(idu2dD,ng)) THEN
          DO nk=0,NTC2
            DO j=JstrR,JendR
              DO i=Istr,IendR
                TIDES(ng)%ubar_tide(i,j,nk)=TIDES(ng)%ubar_tide(i,j,    &
     &                                                          nk)+    &
#  ifdef WET_DRY
     &                                      GRID(ng)%umask_full(i,j)*   &
#  endif
     &                                      OCEAN(ng)%ubar(i,j,Kout)*   &
     &                                      tide_harmonics(nk)
              END DO
            END DO
          END DO
        END IF

        IF (Aout(idv2dD,ng)) THEN
          DO nk=0,NTC2
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                TIDES(ng)%vbar_tide(i,j,nk)=TIDES(ng)%vbar_tide(i,j,    &
     &                                                          nk)+    &
#  ifdef WET_DRY
     &                                      GRID(ng)%vmask_full(i,j)*   &
#  endif
     &                                      OCEAN(ng)%vbar(i,j,Kout)*   &
     &                                      tide_harmonics(nk)
              END DO
            END DO
          END DO
        END IF

#  ifdef SOLVE3D
        IF (Aout(idu3dD,ng)) THEN
          DO nk=0,NTC2
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=Istr,IendR
                  TIDES(ng)%u_tide(i,j,k,nk)=TIDES(ng)%u_tide(i,j,k,    &
     &                                                        nk)+      &
#   ifdef WET_DRY
     &                                       GRID(ng)%umask_full(i,j)*  &
#   endif
     &                                       OCEAN(ng)%u(i,j,k,Nout)*   &
     &                                       tide_harmonics(nk)
                END DO
              END DO
            END DO
          END DO
        END IF

        IF (Aout(idv3dD,ng)) THEN
          DO nk=0,NTC2
            DO k=1,N(ng)
              DO j=Jstr,JendR
                DO i=IstrR,IendR
                  TIDES(ng)%v_tide(i,j,k,nk)=TIDES(ng)%v_tide(i,j,k,    &
     &                                                        nk)+      &
#   ifdef WET_DRY
     &                                       GRID(ng)%vmask_full(i,j)*  &
#   endif
     &                                       OCEAN(ng)%v(i,j,k,Nout)*   &
     &                                       tide_harmonics(nk)
                END DO
              END DO
            END DO
          END DO
        END IF

        DO it=1,NAT
          IF (Aout(idTrcD(it),ng)) THEN
            DO nk=0,NTC2
              DO k=1,N(ng)
                DO j=Jstr,JendR
                  DO i=IstrR,IendR
                    TIDES(ng)%t_tide(i,j,k,nk,it)=TIDES(ng)%t_tide(i,j, &
     &                                                        k,nk,it)+ &
#   ifdef WET_DRY
     &                                        GRID(ng)%umask_full(i,j)* &
#   endif
     &                                            OCEAN(ng)%t(i,j,k,    &
     &                                                        Nout,it)* &
     &                                            tide_harmonics(nk)
                  END DO
                END DO
              END DO
            END DO
          END IF
        END DO
#  endif
      END IF
!
!-----------------------------------------------------------------------
!  Convert accumulated sums into time-averages, if appropriate.
!-----------------------------------------------------------------------
!
      IF ((iic(ng).gt.ntsAVG(ng)).and.                                  &
     &    (MOD(iic(ng)-1,nAVG(ng)).eq.0).and.                           &
     &    ((iic(ng).ne.ntstart(ng)).or.(nrrec(ng).eq.0))) THEN
!
!  Compute detide least-squares coefficients.  Build coefficient squared
!  matrix C(0:2*NTC,0:2*NTC) to invert. It is 2*NTC because we are
!  solving for real components Ak and Bk. The zero rows and column is
!  for the coefficients associated with the time mean.
!
!        F(t) = Fmean + SUM [ Ak sin(omega(k)*t) ]
!                     + SUM [ Bk cos(omega(k)*t) ]   for k=1:NTC
!
!  In the code below, all the arrays are collapsed into a single
!  dimension index such that:
!
!            k=0               mean term
!            k=1:NTC           sine terms
!            k=NTC+1:2*NTC     cosine terms
!
        IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
          C(0,0)=1.0_r8                     ! time-averaged coefficient
          fac1=1.0_r8/REAL(Hcount(ng),r8)   ! global summation factor
          DO nk=1,NTC
            C(0,nk    )=fac1*SinW_sum(nk)
            C(0,nk+NTC)=fac1*CosW_sum(nk)
            C(nk,0    )=C(0,nk)             ! symmetric
            C(nk+NTC,0)=C(0,nk+NTC)         ! symmetric
            DO mk=1,NTC
              C(mk,nk)=fac1*SinWSinW(mk,nk)
              C(mk,nk+NTC)=fac1*SinWCosW(mk,nk)
              C(mk+NTC,nk)=fac1*SinWCosW(nk,mk)
              C(mk+NTC,nk+NTC)=fac1*CosWCosW(mk,nk)
            END DO
          END DO
          DO nk=0,NTC2
            DO mk=0,NTC2
              C(nk,mk)=C(mk,nk)
            END DO
          END DO
!
!  Invert least-squares coefficient matrix by LU decomposition.
!
          DO mk=0,NTC2
            DO nk=0,NTC2
              Y(mk,nk)=0.0_r8
            END DO
            Y(mk,mk)=1.0_r8                ! identity matrix
          END DO
          CALL ludcmp (C(0,0), NTC2+1, NTC2+1, indx, d)
!
!  Find inverse by columns. The matrix Y will now contain the inverse
!  of the least-squares coefficient matrix C, which will have been
!  destroyed.
!
          DO nk=0,NTC2
            CALL lubksb (C(0,0), NTC2+1, NTC2+1, indx, Y(0,nk))
          END DO
!
!  Compute time-averaged harmonics for current field average window.
!
          tide_harmonics(0)=1.0_r8
          DO nk=1,NTC
            tide_harmonics(nk    )=SinW_avg(nk)
            tide_harmonics(nk+NTC)=cosW_avg(nk)
          END DO
!
!  Scale inverse by the global summation factor.
!
          DO nk=0,NTC2
            DO mk=0,NTC2
              Y(nk,mk)=fac1*Y(nk,mk)
            END DO
          END DO
        END IF
!
!  Set time-averaged factors for each C-grid variable type. Notice that
!  the I- and J-ranges are all grid types are the same for convinience.
# ifdef WET_DRY
!  In wetting and drying, the sums are devided by the number of times
!  that each qrid point is wet. This was already computed in above in
!  routine "set_avg_tile" for nAVG steps time-average window.
# endif
!
# ifdef WET_DRY
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            rfac(i,j)=1.0_r8/MAX(1.0_r8, GRID(ng)%rmask_avg(i,j))
            ufac(i,j)=1.0_r8/MAX(1.0_r8, GRID(ng)%umask_avg(i,j))
            vfac(i,j)=1.0_r8/MAX(1.0_r8, GRID(ng)%vmask_avg(i,j))
          END DO
        END DO
# else
        fac=1.0_r8/REAL(nAVG(ng),r8)
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            rfac(i,j)=fac
            ufac(i,j)=fac
            vfac(i,j)=fac
          END DO
        END DO
# endif
!
!  Process accumulated detided averages. Notice that the regular
!  time-averaged field values used here are the ones computed in
!  routine "set_avg_tile".
!
        IF (Aout(idFsuD,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              Hsum=0.0_r8
              DO nk=0,NTC2
                Ak(nk)=0.0_r8
                DO mk=0,NTC2
!!                Ak(nk)=Ak(nk)+Y(nk,mk)*TIDES(ng)%zeta_tide(i,j,mk)
                  Ak(nk)=Ak(nk)+Y(mk,nk)*TIDES(ng)%zeta_tide(i,j,mk)
                END DO
                Hsum=Hsum+Ak(nk)*tide_harmonics(nk)*rfac(i,j)
              END DO
              TIDES(ng)%zeta_detided(i,j)=AVERAGE(ng)%avgzeta(i,j)-     &
     &                                    Hsum
            END DO
          END DO
        END IF

        IF (Aout(idu2dD,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              Hsum=0.0_r8
              DO nk=0,NTC2
                Ak(nk)=0.0_r8
                DO mk=0,NTC2
!!                Ak(nk)=Ak(nk)+Y(nk,mk)*TIDES(ng)%ubar_tide(i,j,mk)
                  Ak(nk)=Ak(nk)+Y(mk,nk)*TIDES(ng)%ubar_tide(i,j,mk)
                END DO
                Hsum=Hsum+Ak(nk)*tide_harmonics(nk)*ufac(i,j)
              END DO
              TIDES(ng)%ubar_detided(i,j)=AVERAGE(ng)%avgu2d(i,j)-      &
     &                                    Hsum
            END DO
          END DO
        END IF

        IF (Aout(idv2dD,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              Hsum=0.0_r8
              DO nk=0,NTC2
                Ak(nk)=0.0_r8
                DO mk=0,NTC2
!!                Ak(nk)=Ak(nk)+Y(nk,mk)*TIDES(ng)%vbar_tide(i,j,mk)
                  Ak(nk)=Ak(nk)+Y(mk,nk)*TIDES(ng)%vbar_tide(i,j,mk)
                END DO
                Hsum=Hsum+Ak(nk)*tide_harmonics(nk)*vfac(i,j)
              END DO
              TIDES(ng)%vbar_detided(i,j)=AVERAGE(ng)%avgv2d(i,j)-      &
     &                                    Hsum
            END DO
          END DO
        END IF

#  ifdef SOLVE3D
        IF (Aout(idu3dD,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                Hsum=0.0_r8
                DO nk=0,NTC2
                  Ak(nk)=0.0_r8
                  DO mk=0,NTC2
!!                  Ak(nk)=Ak(nk)+Y(nk,mk)*TIDES(ng)%u_tide(i,j,k,mk)
                    Ak(nk)=Ak(nk)+Y(mk,nk)*TIDES(ng)%u_tide(i,j,k,mk)
                  END DO
                  Hsum=Hsum+Ak(nk)*tide_harmonics(nk)*ufac(i,j)
                END DO
                TIDES(ng)%u_detided(i,j,k)=AVERAGE(ng)%avgu3d(i,j,k)-   &
     &                                     Hsum
              END DO
            END DO
          END DO
        END IF

        IF (Aout(idv3dD,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                Hsum=0.0_r8
                DO nk=0,NTC2
                  Ak(nk)=0.0_r8
                  DO mk=0,NTC2
!!                  Ak(nk)=Ak(nk)+Y(nk,mk)*TIDES(ng)%v_tide(i,j,k,mk)
                    Ak(nk)=Ak(nk)+Y(mk,nk)*TIDES(ng)%v_tide(i,j,k,mk)
                  END DO
                  Hsum=Hsum+Ak(nk)*tide_harmonics(nk)*vfac(i,j)
                END DO
                TIDES(ng)%v_detided(i,j,k)=AVERAGE(ng)%avgv3d(i,j,k)-   &
     &                                     Hsum
              END DO
            END DO
          END DO
        END IF

        DO it=1,NAT
          IF (Aout(idTrcD(it),ng)) THEN
            DO k=1,N(ng)
              DO j=Jstr,JendR
                DO i=IstrR,IendR
                  Hsum=0.0_r8
                  DO nk=0,NTC2
                    Ak(nk)=0.0_r8
                    DO mk=0,NTC2
!!                    Ak(nk)=Ak(nk)+                                    &
!!   &                       Y(nk,mk)*TIDES(ng)%t_tide(i,j,k,mk,it)
                      Ak(nk)=Ak(nk)+                                    &
     &                       Y(mk,nk)*TIDES(ng)%t_tide(i,j,k,mk,it)
                    END DO
                    Hsum=Hsum+Ak(nk)*tide_harmonics(nk)*rfac(i,j)
                  END DO
                  TIDES(ng)%t_detided(i,j,k,it)=AVERAGE(ng)%avgt(i,j,k, &
     &                                                           it)-   &
     &                                          Hsum
                END DO
              END DO
            END DO
          END IF
        END DO
#  endif
      END IF

      RETURN
      END SUBROUTINE set_detide_tile
# endif

#endif
      END MODULE set_avg_mod
