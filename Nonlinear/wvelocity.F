#include "cppdefs.h"
      MODULE wvelocity_mod
#ifdef SOLVE3D
!
!========================================== Alexander F. Shchepetkin ===
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This subroutines computes vertical velocity (m/s) at W-points       !
!  from the vertical mass flux (omega*hz/m*n).  This computation       !
!  is done solely for output purposes.                                 !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC wvelocity

      CONTAINS
!
!***********************************************************************
      SUBROUTINE wvelocity (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_coupling
      USE mod_grid
      USE mod_ocean
      USE mod_stepping
!
      integer, intent(in) :: ng, tile

# include "tile.h"

      CALL wvelocity_tile (ng, Istr, Iend, Jstr, Jend,                  &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     nstp(ng),                                    &
     &                     GRID(ng) % pm,                               &
     &                     GRID(ng) % pn,                               &
     &                     GRID(ng) % z_r,                              &
     &                     GRID(ng) % z_w,                              &
     &                     COUPLING(ng) % DU_avg1,                      &
     &                     COUPLING(ng) % DV_avg1,                      &
     &                     OCEAN(ng) % u,                               &
     &                     OCEAN(ng) % v,                               &
     &                     OCEAN(ng) % W,                               &
     &                     OCEAN(ng) % wvel)
      RETURN
      END SUBROUTINE wvelocity
!
!***********************************************************************
      SUBROUTINE wvelocity_tile (ng, Istr, Iend, Jstr,Jend,             &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           nstp,                                  &
     &                           pm, pn, z_r, z_w,                      &
     &                           DU_avg1, DV_avg1,                      &
     &                           u, v, W,                               &
     &                           wvel)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars

# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
!
      USE exchange_2d_mod
      USE exchange_3d_mod, ONLY : exchange_w3d_tile
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: nstp

# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(inout) :: DU_avg1(LBi:,LBj:)
      real(r8), intent(inout) :: DV_avg1(LBi:,LBj:)
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
      real(r8), intent(in) :: W(LBi:,LBj:,0:)
      real(r8), intent(out) :: wvel(LBi:,LBj:,0:)
# else
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(inout) :: DU_avg1(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: DV_avg1(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: W(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(out) :: wvel(LBi:UBi,LBj:UBj,0:N(ng))
# endif

!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j, k

      real(r8) :: cff1, cff2, cff3

      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY,N(ng)) :: vert

      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: wrk

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute "true" vertical velocity.
!-----------------------------------------------------------------------
    
# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
!
!  Exchange time-averaged fields.
!
      CALL exchange_u2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, DU_avg1)
      CALL exchange_v2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, DV_avg1)
# endif
!
!  Compute contribution due to quasi-horizontal motions along
!  S-coordinate surfaces:  (Ui +Vj)*GRADs(z).
!
      DO k=1,N(ng)
        DO j=Jstr,Jend
          DO i=Istr,Iend+1
            wrk(i,j)=u(i,j,k,nstp)*(z_r(i,j,k)-z_r(i-1,j,k))*           &
     &                             (pm(i-1,j)+pm(i,j))
          END DO
          DO i=Istr,Iend
            vert(i,j,k)=0.25_r8*(wrk(i,j)+wrk(i+1,j))
          END DO
        END DO
        DO j=Jstr,Jend+1
          DO i=Istr,Iend
            wrk(i,j)=v(i,j,k,nstp)*(z_r(i,j,k)-z_r(i,j-1,k))*           &
     &                             (pn(i,j-1)+pn(i,j))
          END DO
        END DO
        DO j=Jstr,Jend
          DO i=Istr,Iend
            vert(i,j,k)=vert(i,j,k)+0.25_r8*(wrk(i,j)+wrk(i,j+1))
          END DO
        END DO
      END DO
!
!  Compute contribution due to time tendency of the free-surface,
!  d(zeta)/d(t), which is the vertical velocity at the free-surface
!  and it is expressed in terms of barotropic mass flux divergence.
!  Notice that it is divided by the total depth of the water column.
!  This is needed because this contribution is linearly distributed
!  throughout the water column by multiplying it by the distance from
!  the bottom to the depth at which the vertical velocity is computed.
!
      J_LOOP : DO j=Jstr,Jend
        DO i=Istr,Iend
          wrk(i,j)=(DU_avg1(i,j)-DU_avg1(i+1,j)+                        &
     &              DV_avg1(i,j)-DV_avg1(i,j+1))/                       &
     &             (z_w(i,j,N(ng))-z_w(i,j,0))
        END DO
!
!  Notice that a cubic interpolation is used to shift the "vert"
!  contribution from vertical RHO- to W-points.
!
        cff1=3.0_r8/8.0_r8
        cff2=3.0_r8/4.0_r8
        cff3=1.0_r8/8.0_r8
        DO i=Istr,Iend
!!        wvel(i,j,0)=0.0_r8
          wvel(i,j,0)=cff1*(vert(i,j,1)-                                &
     &                      (vert(i,j,2)-vert(i,j,1))*                  &
     &                      2.0_r8*(z_r(i,j,1)-z_w(i,j,0)))+            &
     &                cff2*vert(i,j,1)-                                 &
     &                cff3*vert(i,j,2)
          wvel(i,j,1)=pm(i,j)*pn(i,j)*                                  &
     &                (W(i,j,1)+                                        &
     &                 wrk(i,j)*(z_w(i,j,1)-z_w(i,j,0)))+               &
     &                cff1*vert(i,j,1)+                                 &
     &                cff2*vert(i,j,2)-                                 &
     &                cff3*vert(i,j,3)
        END DO
        cff1=9.0_r8/16.0_r8
        cff2=1.0_r8/16.0_r8
        DO k=2,N(ng)-2
          DO i=Istr,Iend
            wvel(i,j,k)=pm(i,j)*pn(i,j)*                                &
     &                  (W(i,j,k)+                                      &
     &                   wrk(i,j)*(z_w(i,j,k)-z_w(i,j,0)))+             &
     &                  cff1*(vert(i,j,k  )+vert(i,j,k+1))-             &
     &                  cff2*(vert(i,j,k-1)+vert(i,j,k+2))
          END DO
        END DO
        cff1=3.0_r8/8.0_r8
        cff2=3.0_r8/4.0_r8
        cff3=1.0_r8/8.0_r8
        DO i=Istr,Iend
          wvel(i,j,N(ng))=pm(i,j)*pn(i,j)*                              &
     &                    wrk(i,j)*(z_w(i,j,N(ng))-z_w(i,j,0))+         &
     &                    cff1*(vert(i,j,N(ng))+                        &
     &                          (vert(i,j,N(ng))-vert(i,j,N(ng)-1))*    &
     &                          2.0_r8*(z_w(i,j,N(ng))-                 &
     &                                  z_r(i,j,N(ng))))+               &
     &                    cff2*vert(i,j,N(ng)  )-                       &
     &                    cff3*vert(i,j,N(ng)-1)
          wvel(i,j,N(ng)-1)=pm(i,j)*pn(i,j)*                            &
     &                      (W(i,j,N(ng)-1)+                            &
     &                       wrk(i,j)*(z_w(i,j,N(ng)-1)-z_w(i,j,0)))+   &
     &                      cff1*vert(i,j,N(ng)  )+                     &
     &                      cff2*vert(i,j,N(ng)-1)-                     &
     &                      cff3*vert(i,j,N(ng)-2)
        END DO
      END DO J_LOOP
!
!  Set lateral boundary conditions.
!
# ifndef EW_PERIODIC
      IF (WESTERN_EDGE) THEN
        DO k=0,N(ng)
          DO j=Jstr,Jend
            wvel(0,j,k)=wvel(1,j,k)
          END DO
        END DO
      END IF
      IF (EASTERN_EDGE) THEN
        DO k=0,N(ng)
          DO j=Jstr,Jend
            wvel(Lm(ng)+1,j,k)=wvel(Lm(ng),j,k)
          END DO
        END DO
      END IF
# endif
# ifndef NS_PERIODIC
#  ifdef EW_PERIODIC
#   define I_RANGE Istr,Iend
#  else
#   define I_RANGE IstrR,IendR
#  endif
      IF (SOUTHERN_EDGE) THEN
        DO k=0,N(ng)
          DO i=I_RANGE
            wvel(i,0,k)=wvel(i,1,k)
          END DO
        END DO
      END IF
      IF (NORTHERN_EDGE) THEN
        DO k=0,N(ng)
          DO i=I_RANGE
            wvel(i,Mm(ng)+1,k)=wvel(i,Mm(ng),k)
          END DO
        END DO
      END IF
#  undef I_RANGE
# endif
# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      CALL exchange_w3d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, wvel)
# endif
      RETURN
      END SUBROUTINE wvelocity_tile
#endif
      END MODULE wvelocity_mod
