#include "cppdefs.h"
      MODULE tl_omega_mod
#if defined TANGENT && defined SOLVE3D
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2024 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  This routine computes S-coordinate vertical velocity (m^3/s),       !
!                                                                      !
!                  W=[Hz/(m*n)]*omega,                                 !
!                                                                      !
!  diagnostically at horizontal RHO-points and vertical W-points.      !
!                                                                      !
# ifdef OMEGA_IMPLICIT
!  Added implicit vertical advection following Shchepetkin (2015).     !
!  The vertical velocity is split into explicit (W) and implicit (Wi)  !
!  parts that adjust automatically to the local flow conditions based  !
!  on the Courant number for stability, allowing larger time steps.    !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    Shchepetkin, A.F., 2015: An adaptive, Courant-number-dependent    !
!      implicit scheme for vertical advection in oceanic modeling,     !
!      Ocean Modelling, 91, 38-69, doi: 10.1016/j.ocemod.2015.03.006   !
!                                                                      !
# endif
!  BASIC STATE variables needed: W, z_w.                               !
!                                                                      !
!  NOTE: We need to recompute  basic state W in this routine since     !
!  ----  intermediate values of W are needed by the tangent linear     !
!        and adjoint routines.                                         !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC  :: tl_omega
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE tl_omega (ng, tile, model)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_ocean
# if defined SEDIMENT && defined SED_MORPH
      USE mod_sedbed
      USE mod_stepping
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, model, 13, __LINE__, MyFile)
# endif
      CALL tl_omega_tile (ng, tile, model,                              &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
# if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
     &                    nstp(ng), nnew(ng),                           &
     &                    GRID(ng) % omn,                               &
     &                    SEDBED(ng) % bed_thick,                       &
     &                    SEDBED(ng) % tl_bed_thick,                    &
# endif
     &                    GRID(ng) % Huon,                              &
     &                    GRID(ng) % Hvom,                              &
# ifdef OMEGA_IMPLICIT
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
# endif
     &                    GRID(ng) % z_w,                               &
     &                    GRID(ng) % tl_Huon,                           &
     &                    GRID(ng) % tl_Hvom,                           &
     &                    GRID(ng) % tl_z_w,                            &
# ifdef OMEGA_IMPLICIT
     &                    OCEAN(ng) % Wi,                               &
     &                    OCEAN(ng) % tl_Wi,                            &
# endif
     &                    OCEAN(ng) % W,                                &
     &                    OCEAN(ng) % tl_W)
# ifdef PROFILE
      CALL wclock_off (ng, model, 13, __LINE__, MyFile)
# endif
!
      RETURN
      END SUBROUTINE tl_omega
!
!***********************************************************************
      SUBROUTINE tl_omega_tile (ng, tile, model,                        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
# if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
     &                          nstp, nnew,                             &
     &                          omn,
     &                          bed_thick, tl_bed_thick,                &
# endif
     &                          Huon, Hvom,                             &
# ifdef OMEGA_IMPLICIT
     &                          pm, pn,                                 &
# endif
     &                          z_w,                                    &
     &                          tl_Huon, tl_Hvom, tl_z_w,               &
# ifdef OMEGA_IMPLICIT
     &                          Wi, tl_Wi,                              &
# endif
     &                          W, tl_W)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_sources
!
      USE bc_3d_mod, ONLY : bc_w3d_tile
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange3d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
# if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
      integer, intent(in) :: nstp, nnew
# endif
!
# ifdef ASSUMED_SHAPE
#  if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
      real(r8), intent(in) :: omn(LBi:,LBj:)
      real(r8), intent(in) :: bed_thick(LBi:,LBj:,:)
      real(r8), intent(in) :: tl_bed_thick(LBi:,LBj:,:)
#  endif
      real(r8), intent(in) :: Huon(LBi:,LBj:,:)
      real(r8), intent(in) :: Hvom(LBi:,LBj:,:)
#  ifdef OMEGA_IMPLICIT
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: tl_Huon(LBi:,LBj:,:)
      real(r8), intent(in) :: tl_Hvom(LBi:,LBj:,:)
      real(r8), intent(in) :: tl_z_w(LBi:,LBj:,0:)
#  ifdef OMEGA_IMPLICIT
      real(r8), intent(out) :: Wi(LBi:,LBj:,0:)
      real(r8), intent(out) :: tl_Wi(LBi:,LBj:,0:)
#  endif
      real(r8), intent(out) :: W(LBi:,LBj:,0:)
      real(r8), intent(out) :: tl_W(LBi:,LBj:,0:)

# else

#  if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
      real(r8), intent(in) :: omn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: bed_thick(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: tl_bed_thick(LBi:UBi,LBj:UBj,3)
#  endif
      real(r8), intent(in) :: Huon(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: Hvom(LBi:UBi,LBj:UBj,N(ng))
#  ifdef OMEGA_IMPLICIT
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(in) :: tl_Huon(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: tl_Hvom(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: tl_z_w(LBi:UBi,LBj:UBj,0:N(ng))
#  ifdef OMEGA_IMPLICIT
      real(r8), intent(out) :: Wi(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(out) :: tl_Wi(LBi:UBi,LBj:UBj,0:N(ng))
#  endif
      real(r8), intent(out) :: W(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(out) :: tl_W(LBi:UBi,LBj:UBj,0:N(ng))
# endif
!
!  Local variable declarations.
!
      integer :: i, ii, is, j, jj, k
      real(r8) :: cff, tl_cff
# if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
      real(r8) :: cff1
# endif
      real(r8), dimension(IminS:ImaxS) :: wrk
      real(r8), dimension(IminS:ImaxS) :: tl_wrk
# ifdef OMEGA_IMPLICIT
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: Cu_adv
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: tl_Cu_adv
      real(r8) :: cw, c2d, dh, cw_max, cw_max2, cw_min 
      real(r8) :: tl_dh
!
      real(r8), parameter :: amax = 0.75_r8    ! Maximum Courant number
      real(r8), parameter :: amin = 0.60_r8    ! Minimum Courant number
      real(r8), parameter :: cmnx_ratio = amin/amax
      real(r8), parameter :: cutoff = 2.0_r8-amin/amax
      real(r8), parameter :: r4cmx = 1.0_r8/(4.0_r8-4.0_r8*amin/amax)
# endif

# include "set_bounds.h"
!
!------------------------------------------------------------------------
!  Vertically integrate horizontal mass flux divergence.
!
!  Starting with zero vertical velocity at the bottom, integrate
!  from the bottom (k=0) to the free-surface (k=N).  The w(:,:,N(ng))
!  contains the vertical velocity at the free-surface, d(zeta)/d(t).
!  Notice that barotropic mass flux divergence is not used directly.
!------------------------------------------------------------------------
!
# if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
!  For sediment bed change, we need to include the mass change of
!  water volume due to change of the sea floor. This is similar to
!  the LwSrc point source approach.
!
      cff1=1.0_r8/dt(ng)
!
# endif
      DO j=Jstr,Jend
        DO i=Istr,Iend
          W(i,j,0)=0.0_r8
          tl_W(i,j,0)=0.0_r8
# if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
          wrk(i)=cff1*(bed_thick(i,j,nstp)-                             &
     &                 bed_thick(i,j,nnew))*omn(i,j)
          tl_wrk(i)=cff1*(tl_bed_thick(i,j,nstp)-                       &
     &                    tl_bed_thick(i,j,nnew))*omn(i,j)
# endif
        END DO
!
!  Code added to clear tl_W to be consistent with adjoint.
!
        DO k=1,N(ng)
          DO i=Istr,Iend
           tl_W(i,j,k)=0.0_r8
          END DO
        END DO
        DO k=1,N(ng)
          DO i=Istr,Iend
            W(i,j,k)=W(i,j,k-1)-                                        &
# if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
     &               wrk(i)-                                            &
# endif
     &               (Huon(i+1,j,k)-Huon(i,j,k)+                        &
     &                Hvom(i,j+1,k)-Hvom(i,j,k))
            tl_W(i,j,k)=tl_W(i,j,k-1)-                                  &
# if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
     &                  tl_wrk(i)-                                      &
# endif
     &                  (tl_Huon(i+1,j,k)-tl_Huon(i,j,k)+               &
     &                   tl_Hvom(i,j+1,k)-tl_Hvom(i,j,k))
# ifdef OMEGA_IMPLICIT
!
!  Compute the horizontal Courant number.
!
            Cu_adv(i,k)=MAX(Huon(i+1,j  ,k),0.0_r8)-                    &
     &                  MIN(Huon(i  ,j  ,k),0.0_r8)+                    &
     &                  MAX(Hvom(i  ,j+1,k),0.0_r8)-                    &
     &                  MIN(Hvom(i  ,j  ,k),0.0_r8)
            tl_Cu_adv(i,k)=(0.5_r8+SIGN(0.5_r8, Huon(i+1,j  ,k)))*      &
     &                     tl_Huon(i+1,j  ,k)-                          &
                           (0.5_r8+SIGN(0.5_r8,-Huon(i  ,j  ,k)))*      &
     &                     tl_Huon(i  ,j  ,k)+                          &
                           (0.5_r8+SIGN(0.5_r8, Hvom(i  ,j+1,k)))*      &
     &                     tl_Hvom(i  ,j+1,k)-                          &
                           (0.5_r8+SIGN(0.5_r8,-Hvom(i  ,j  ,k)))*      &
     &                     tl_Hvom(i  ,j  ,k)
# endif
          END DO
        END DO
!
!  Apply mass point sources (volume vertical influx), if any.
!
!  Overwrite W(Isrc,Jsrc,k) with the same divergence of Huon,Hvom as
!  above but add in point source Qsrc(k) and reaccumulate the vertical
!  sum to obtain the correct net Qbar given in user input - J. Levin
!  (Jupiter Intelligence Inc.) and J. Wilkin
!
!    Dsrc(is) = 2,  flow across grid cell w-face (positive or negative)
!
        IF (LwSrc(ng)) THEN
          DO is=1,Nsrc(ng)
            IF (INT(SOURCES(ng)%Dsrc(is)).eq.2) THEN
              ii=SOURCES(ng)%Isrc(is)
              jj=SOURCES(ng)%Jsrc(is)
              IF (((IstrR.le.ii).and.(ii.le.IendR)).and.                &
     &            ((JstrR.le.jj).and.(jj.le.JendR)).and.                &
     &            (j.eq.jj)) THEN
# if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
                wrk(ii)=cff1*(bed_thick(ii,jj,nstp)-                    &
     &                        bed_thick(ii,jj,nnew))*omn(ii,jj)
                tl_wrk(ii)=cff1*(tl_bed_thick(ii,jj,nstp)-              &
     &                           tl_bed_thick(ii,jj,nnew))*omn(ii,jj)
# endif
                DO k=1,N(ng)
                  W(ii,jj,k)=W(ii,jj,k-1)-                              &
# if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
     &                       wrk(ii)-                                   &
# endif
     &                       (Huon(ii+1,jj,k)-Huon(ii,jj,k)+            &
     &                        Hvom(ii,jj+1,k)-Hvom(ii,jj,k))+           &
     &                       SOURCES(ng)%Qsrc(is,k)
                  tl_W(ii,jj,k)=tl_W(ii,jj,k-1)-                        &
# if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
     &                          tl_wrk(ii)-                             &
# endif
     &                          (tl_Huon(ii+1,jj,k)-tl_Huon(ii,jj,k)+   &
     &                           tl_Hvom(ii,jj+1,k)-tl_Hvom(ii,jj,k))+  &
     &                          SOURCES(ng)%tl_Qsrc(is,k)
                END DO
              END IF
            END IF
          END DO
        END IF
!
        DO i=Istr,Iend
          cff=1.0_r8/(z_w(i,j,N(ng))-z_w(i,j,0))
          tl_cff=-cff*cff*(tl_z_w(i,j,N(ng))-tl_z_w(i,j,0))
          wrk(i)=cff*W(i,j,N(ng))
          tl_wrk(i)=tl_cff*W(i,j,N(ng))+cff*tl_W(i,j,N(ng))
# ifdef OMEGA_IMPLICIT
          Cu_adv(i,0)=dt(ng)*pm(i,j)*pn(i,j)
# endif
        END DO
!
!  In order to insure zero vertical velocity at the free-surface,
!  subtract the vertical velocities of the moving S-coordinates
!  isosurfaces. These isosurfaces are proportional to d(zeta)/d(t).
!  The proportionally coefficients are a linear function of the
!  S-coordinate with zero value at the bottom (k=0) and unity at
!  the free-surface (k=N).
!
        DO k=N(ng)-1,1,-1
          DO i=Istr,Iend
            W(i,j,k)=W(i,j,k)-                                          &
     &               wrk(i)*(z_w(i,j,k)-z_w(i,j,0))
            tl_W(i,j,k)=tl_W(i,j,k)-                                    &
     &                  tl_wrk(i)*(z_w(i,j,k)-z_w(i,j,0))-              &
     &                  wrk(i)*(tl_z_w(i,j,k)-tl_z_w(i,j,0))
# ifdef OMEGA_IMPLICIT
!
!  Determine implicit part "Wi" of vertical advection: "W"  becomes the
!  explicit part "We".
!
!  HGA: 'dh' cannot be linearized because it is the local vertical grid
!       spacing. It has to be positive!  We need to take the nonlinear
!       value.
!       Andy what should we do here for the linearization?
!
            Wi(i,j,k)=W(i,j,k)
            tl_Wi(i,j,k)=tl_W(i,j,k)
            IF (Wi(i,j,k).ge.0.0_r8) THEN     ! Three different variants
              c2d=Cu_adv(i,k)                 ! for computing 2D Courant
              dh=z_w(i,j,k)-z_w(i,j,k-1)      ! number at the interface:
            ELSE                              ! (1) use value from the 
              c2d=Cu_adv(i,k+1)               !     grid box upstream in
              dh=z_w(i,j,k+1)-z_w(i,j,k)      !     vertical direction;
            END IF
!
!!          c2d=0.5_r8*(Cu_adv(i,k  )+                                  &
!!   &                  Cu_adv(i,k+1))        ! (2) average the two; or
!!          dh=0.5_r8*(z_w(i,j,k+1)-                                    &
!!   &                 z_w(i,j,k-1))
!!
!!          c2d=MAX(Cu_adv(i,k  ),                                      &
!!   &              Cu_adv(i,k+1))            ! (3) pick its maximum
!!          dh=MIN(z_w(i,j,k+1)-z_w(i,j,k),                             &
!!                 z_w(i,j,k)-z_w(i,j,k-1))
!!
            cw_max=amax*dh-c2d*Cu_adv(i,0)    ! compare the vertical
            IF (cw_max.ge.0.0_r8) THEN        ! displacement to dz*amax.
              cw_max2=cw_max*cw_max           ! Partition W into Wi and
              cw_min=cw_max*cmnx_ratio        ! We.
              cw=ABS(Wi(i,j,k))*Cu_adv(i,0)
              IF (cw.le.cw_min) THEN
                cff=cw_max2
              ELSE IF (cw.le.cutoff*cw_max) THEN
                cff=cw_max2+r4cmx*(cw-cw_min)**2
              ELSE
                cff=cw_max*cw
              END IF
!
              W(i,j,k)=cw_max2*Wi(i,j,k)/cff
              tl_W(i,j,k)=cw_max2*tl_Wi(i,j,k)/cff
              Wi(i,j,k)=Wi(i,j,k)-W(i,j,k)
              tl_Wi(i,j,k)=tl_Wi(i,j,k)-tl_W(i,j,k)
            ELSE                              ! All the displacement is
              W(i,j,k)=0.0_r8                 ! greater than amax*dz, so
              tl_W(i,j,k)=0.0_r8              ! keep it all into Wi.
            END IF
# endif
          END DO
        END DO
        DO i=Istr,Iend
          W(i,j,N(ng))=0.0_r8
          tl_W(i,j,N(ng))=0.0_r8
        END DO
      END DO
!
!  Set lateral boundary conditions.
!
      CALL bc_w3d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, 0, N(ng),                   &
     &                  W)
      CALL bc_w3d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, 0, N(ng),                   &
     &                  tl_W)
# ifdef OMEGA_IMPLICIT
      CALL bc_w3d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, 0, N(ng),                   &
     &                  Wi)
      CALL bc_w3d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, 0, N(ng),                   &
     &                  tl_Wi)
# endif
# ifdef DISTRIBUTE
      CALL mp_exchange3d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj, 0, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    W)
      CALL mp_exchange3d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj, 0, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_W)
#  ifdef OMEGA_IMPLICIT
      CALL mp_exchange3d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj, 0, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    Wi)
      CALL mp_exchange3d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj, 0, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_Wi)
#  endif
# endif
!
      RETURN
      END SUBROUTINE tl_omega_tile
#endif
      END MODULE tl_omega_mod
