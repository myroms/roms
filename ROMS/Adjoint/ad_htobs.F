#include "cppdefs.h"
      MODULE ad_htobs_mod
#if defined ADJOINT && (defined WEAK_CONSTRAINT || defined IOM)
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2007 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes the adjoint observation operator,             !
!                                                                      !
!       transpose(H) * X                                               !
!                                                                      !
!  which loads the observation vector into the adjoint forcing arrays. !
!                                                                      !
!=======================================================================
!
      implicit none

      CONTAINS
!
!***********************************************************************
      SUBROUTINE ad_htobs (ng, tile, model)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_ocean
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
# include "tile.h"
!
      CALL ad_htobs_tile (ng, model, Istr, Iend, Jstr, Jend,            &
     &                    LBi, UBi, LBj, UBj,                           &
# ifdef MASKING
     &                    GRID(ng) % rmask,                             &
     &                    GRID(ng) % umask,                             &
     &                    GRID(ng) % vmask,                             &
# endif
# ifdef SOLVE3D
     &                    GRID(ng) % z_r,                               &
     &                    OCEAN(ng) % f_u,                              &
     &                    OCEAN(ng) % f_v,                              &
     &                    OCEAN(ng) % f_t,                              &
# endif
     &                    OCEAN(ng) % f_ubar,                           &
     &                    OCEAN(ng) % f_vbar,                           &
     &                    OCEAN(ng) % f_zeta)
      RETURN
      END SUBROUTINE ad_htobs
!
!***********************************************************************
      SUBROUTINE ad_htobs_tile (ng, model, Istr, Iend, Jstr, Jend,      &
     &                          LBi, UBi, LBj, UBj,                     &
# ifdef MASKING
     &                          rmask, umask, vmask,                    &
# endif
# ifdef SOLVE3D
     &                          z_r,                                    &
     &                          f_u, f_v, f_t,                          &
# endif
     &                          f_ubar, f_vbar, f_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
!
# ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_collect
# endif
      USE ad_extract_obs_mod, ONLY : ad_extract_obs2d
# ifdef SOLVE3D
      USE ad_extract_obs_mod, ONLY : ad_extract_obs3d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
#  ifdef SOLVE3D
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(inout) :: f_u(LBi:,LBj:,:)
      real(r8), intent(inout) :: f_v(LBi:,LBj:,:)
      real(r8), intent(inout) :: f_t(LBi:,LBj:,:,:)
#  endif
      real(r8), intent(inout) :: f_ubar(LBi:,LBj:)
      real(r8), intent(inout) :: f_vbar(LBi:,LBj:)
      real(r8), intent(inout) :: f_zeta(LBi:,LBj:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef SOLVE3D
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: f_u(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: f_v(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: f_t(LBi:UBi,LBj:UBj,N(ng),NT(ng))
#  endif
      real(r8), intent(inout) :: f_ubar(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: f_vbar(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: f_zeta(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: Mstr, Mend, ObsSum, ObsVoid
# ifdef DISTRIBUTE
      integer :: Ncollect
# endif
      integer :: i, ie, iobs, is, j

# ifdef SOLVE3D
      integer :: itrc, k

      real(r8), allocatable :: z_vel(:,:,:)
# endif
      real(r8), parameter :: IniVal = 0.0_r8

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute model minus observations adjoint misfit forcing. The
!  representer coefficients (or its approximation PSI) has been
!  already loaded into vector ADmodVal in the conjugate gradient
!  or read in.
!-----------------------------------------------------------------------
!
      IF (ProcessObs(ng)) THEN
!
!  Set starting and ending indices of representer coefficient vector to
!  proccess. The adjoint forcing is only computed for current time
!  survey observations.
!
        Mstr=NstrObs(ng)
        Mend=NendObs(ng)
!
!  Free-surface.
!
        DO i=LBi,UBi
          DO j=LBj,UBj
            f_zeta(i,j)=0.0_r8
          END DO
        END DO
        IF (FOURDVAR(ng)%ObsCount(isFsur).gt.0) THEN
          CALL ad_extract_obs2d (ng, model, Istr, Iend, Istr, Jend,     &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           rILB(ng), rIUB(ng),                    &
     &                           rJLB(ng), rJUB(ng),                    &
     &                           isFsur,                                &
     &                           Mobs, Mstr, Mend,                      &
     &                           time(ng), dt(ng),                      &
     &                           ObsType, ObsScale,                     &
     &                           Tobs, Xobs, Yobs,                      &
     &                           f_zeta,                                &
# ifdef MASKING
     &                           rmask,                                 &
# endif
     &                           ADmodVal)
        END IF
!
!  2D u-momentum component.
!
        DO i=LBi,UBi
          DO j=LBj,UBj
            f_ubar(i,j)=0.0_r8
          END DO
        END DO
        IF (FOURDVAR(ng)%ObsCount(isUbar).gt.0) THEN
          CALL ad_extract_obs2d (ng, model, Istr, Iend, Jstr, Jend,     &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           uILB(ng), uIUB(ng),                    &
     &                           uJLB(ng), uJUB(ng),                    &
     &                           isUbar,                                &
     &                           Mobs, Mstr, Mend,                      &
     &                           time(ng), dt(ng),                      &
     &                           ObsType, ObsScale,                     &
     &                           Tobs, Xobs, Yobs,                      &
     &                           f_ubar,                                &
# ifdef MASKING
     &                           umask,                                 &
# endif
     &                           ADmodVal)
        END IF
!
!  2D v-momentum component.
!
        DO i=LBi,UBi
          DO j=LBj,UBj
            f_vbar(i,j)=0.0_r8
          END DO
        END DO
        IF (FOURDVAR(ng)%ObsCount(isVbar).gt.0) THEN
          CALL ad_extract_obs2d (ng, model, Istr, Iend, Jstr, Jend,     &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           vILB(ng), vIUB(ng),                    &
     &                           vJLB(ng), vJUB(ng),                    &
     &                           isVbar,                                &
     &                           Mobs, Mstr, Mend,                      &
     &                           time(ng), dt(ng),                      &
     &                           ObsType, ObsScale,                     &
     &                           Tobs, Xobs, Yobs,                      &
     &                           f_vbar,                                &
# ifdef MASKING
     &                           vmask,                                 &
# endif
     &                           ADmodVal)
        END IF
# ifdef SOLVE3D
!
!  3D u-momentum component.
!
        DO k=1,N(ng)
          DO i=LBi,UBi
            DO j=LBj,UBj
              f_u(i,j,k)=0.0_r8
            END DO
          END DO
        END DO
        IF (FOURDVAR(ng)%ObsCount(isUvel).gt.0) THEN
          IF (Load_Zobs(ng)) THEN
            IF (.not.allocated(z_vel)) THEN
              allocate (z_vel(LBi:UBi,LBj:UBj,1:N(ng)))
              z_vel=0.0_r8
            END IF
            DO k=1,N(ng)
              DO j=Jstr-1,Jend+1
                DO i=IstrU-1,Iend+1
                  z_vel(i,j,k)=0.5_r8*(z_r(i-1,j,k)+                    &
     &                                 z_r(i  ,j,k))
                END DO
              END DO
            END DO
          END IF
          CALL ad_extract_obs3d (ng, model, Istr, Iend, Jstr, Jend,     &
     &                           LBi, UBi, LBj, UBj, 1, N(ng),          &
     &                           uILB(ng), uIUB(ng),                    &
     &                           uJLB(ng), uJUB(ng),                    &
     &                           isUvel,                                &
     &                           Mobs, Mstr, Mend,                      &
     &                           time(ng), dt(ng),                      &
     &                           ObsType, ObsScale,                     &
     &                           Tobs, Xobs, Yobs, Zobs,                &
     &                           f_u, z_vel,                            &
#  ifdef MASKING
     &                           umask,                                 &
#  endif
     &                           ADmodVal)
        END IF
!
!  3D v-momentum component.
!
        DO k=1,N(ng)
          DO i=LBi,UBi
            DO j=LBj,UBj
              f_v(i,j,k)=0.0_r8
            END DO
          END DO
        END DO
        IF (FOURDVAR(ng)%ObsCount(isVvel).gt.0) THEN
          IF (Load_Zobs(ng)) THEN
            IF (.not.allocated(z_vel)) THEN
              allocate (z_vel(LBi:UBi,LBj:UBj,1:N(ng)))
              z_vel=0.0_r8
            END IF
            DO k=1,N(ng)
              DO j=JstrV-1,Jend+1
                DO i=Istr-1,Iend+1
                  z_vel(i,j,k)=0.5_r8*(z_r(i,j-1,k)+                    &
     &                                 z_r(i,j  ,k))
                END DO
              END DO
            END DO
          END IF
          CALL ad_extract_obs3d (ng, model, Istr, Iend, Jstr, Jend,     &
     &                           LBi, UBi, LBj, UBj, 1, N(ng),          &
     &                           vILB(ng), vIUB(ng),                    &
     &                           vJLB(ng), vJUB(ng),                    &
     &                           isVvel,                                &
     &                           Mobs, Mstr, Mend,                      &
     &                           time(ng), dt(ng),                      &
     &                           ObsType, ObsScale,                     &
     &                           Tobs, Xobs, Yobs, Zobs,                &
     &                           f_v, z_r,                              &
#  ifdef MASKING
     &                           vmask,                                 &
#  endif
     &                           ADmodVal)
        END IF
!
!  Tracer type variables.
!
        DO itrc=1,NT(ng)
          DO k=1,N(ng)
            DO i=LBi,UBi
              DO j=LBj,UBj
                f_t(i,j,k,itrc)=0.0_r8
              END DO
            END DO
          END DO
          IF (FOURDVAR(ng)%ObsCount(isTvar(itrc)).gt.0) THEN
            CALL ad_extract_obs3d (ng, model, Istr, Iend, Jstr, Jend,   &
     &                             LBi, UBi, LBj, UBj, 1, N(ng),        &
     &                             rILB(ng), rIUB(ng),                  &
     &                             rJLB(ng), rJUB(ng),                  &
     &                             isTvar(itrc),                        &
     &                             Mobs, Mstr, Mend,                    &
     &                             time(ng), dt(ng),                    &
     &                             ObsType, ObsScale,                   &
     &                             Tobs, Xobs, Yobs, Zobs,              &
     &                             f_t(:,:,:,itrc), z_r,                &
#  ifdef MASKING
     &                             rmask,                               &
#  endif
     &                             ADmodVal)
          END IF
        END DO
# endif
# ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Collect all observations scaling factors.
!-----------------------------------------------------------------------
!
        Ncollect=Mend-Mstr+1
        CALL mp_collect (ng, model, Ncollect, IniVal,                   &
     &                   ObsScale(Mstr:))
# endif
!
!-----------------------------------------------------------------------
!  Set counters for number of rejected observations for each state
!  variable.
!-----------------------------------------------------------------------
!
        DO iobs=Mstr,Mend
          IF (ObsScale(iobs).lt.1.0) THEN
            IF  (ObsType(iobs).eq.isFsur) THEN
              FOURDVAR(ng)%ObsReject(isFsur)=                           &
     &                              FOURDVAR(ng)%ObsReject(isFsur)+1
            ELSE IF (ObsType(iobs).eq.isUbar) THEN
              FOURDVAR(ng)%ObsReject(isUbar)=                           &
     &                              FOURDVAR(ng)%ObsReject(isUbar)+1
            ELSE IF (ObsType(iobs).eq.isVbar) THEN
              FOURDVAR(ng)%ObsReject(isVbar)=                           &
     &                              FOURDVAR(ng)%ObsReject(isVbar)+1
# ifdef SOLVE3D
            ELSE IF (ObsType(iobs).eq.isUvel) THEN
              FOURDVAR(ng)%ObsReject(isUvel)=                           &
     &                              FOURDVAR(ng)%ObsReject(isUvel)+1
            ELSE IF (ObsType(iobs).eq.isVvel) THEN
              FOURDVAR(ng)%ObsReject(isVvel)=                           &
     &                              FOURDVAR(ng)%ObsReject(isVvel)+1
            ELSE
              DO itrc=1,NT(ng)
                IF (ObsType(iobs).eq.isTvar(itrc)) THEN
                  i=isTvar(itrc)
                  FOURDVAR(ng)%ObsReject(i)=FOURDVAR(ng)%ObsReject(i)+1
                END IF
              END DO
# endif
            END IF
          END IF
        END DO
!
!  Load total available and rejected observations into structure
!  array.
!
        DO i=1,NstateVar(ng)
          FOURDVAR(ng)%ObsCount(0)=FOURDVAR(ng)%ObsCount(0)+            &
     &                             FOURDVAR(ng)%ObsCount(i)
          FOURDVAR(ng)%ObsReject(0)=FOURDVAR(ng)%ObsReject(0)+          &
     &                              FOURDVAR(ng)%ObsReject(i)
        END DO
!
!  Report.
!
        IF (SOUTH_WEST_TEST) THEN
          IF (Master) THEN
            ObsSum=0
            ObsVoid=0
            is=NstrObs(ng)
            DO i=1,NstateVar(ng)
              IF (FOURDVAR(ng)%ObsCount(i).gt.0) THEN
                ie=is+FOURDVAR(ng)%ObsCount(i)-1
                WRITE (stdout,10) TRIM(Vname(1,idSvar(i))), is, ie,     &
     &                            ie-is+1, FOURDVAR(ng)%ObsReject(i)
                is=ie+1
                ObsSum=ObsSum+FOURDVAR(ng)%ObsCount(i)
                ObsVoid=ObsVoid+FOURDVAR(ng)%ObsReject(i)
              END IF
            END DO
            WRITE (stdout,20) ObsSum, ObsVoid,                          &
     &                        FOURDVAR(ng)%ObsCount(0),                 &
     &                        FOURDVAR(ng)%ObsReject(0)
            WRITE (stdout,30) tdays(ng), NstrObs(ng), NendObs(ng),      &
     &                        iic(ng)
  10        FORMAT (10x,a,t25,4(1x,i10))
  20        FORMAT (/,10x,'Total',t47,2(1x,i10),                        &
     &              /,10x,'Obs Tally',t47,2(1x,i10),/)
  30        FORMAT (3x,' AD_HTOBS    - Computed adjoint observations ', &
     &              'forcing,',t64,'t = ',f12.4,/,19x,'(Observation ',  &
     &              'records = ',i7.7,' - ',i7.7,', iic = ',i7.7,')')
          END IF
        END IF
      END IF        
      RETURN
      END SUBROUTINE ad_htobs_tile
#endif
      END MODULE ad_htobs_mod