#include "cppdefs.h"
      MODULE step3d_uv_mod
#ifdef SOLVE3D
!
!========================================== Alexander F. Shchepetkin ===
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine  performs the timestep for the horizontal momentum     !
!  equations.  The vertical viscosity terms are time-stepped using     !
!  an implicit algorithm.                                              !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC step3d_uv

      CONTAINS
!
!***********************************************************************
      SUBROUTINE step3d_uv (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_coupling
      USE mod_forces
      USE mod_grid
      USE mod_mixing
      USE mod_ocean
# ifdef UV_PSOURCE
      USE mod_sources
# endif
      USE mod_stepping
!
      integer, intent(in) :: ng, tile

# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, 34)
# endif
      CALL step3d_uv_tile (ng, Istr, Iend, Jstr, Jend,                  &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     nrhs(ng), nstp(ng), nnew(ng),                &
# ifdef UV_PSOURCE
     &                     Nsrc(ng),                                    &
     &                     SOURCES(ng) % Isrc,                          &
     &                     SOURCES(ng) % Jsrc,                          &
     &                     SOURCES(ng) % Dsrc,                          &
     &                     SOURCES(ng) % Qsrc,                          &
# endif
# ifdef MASKING
     &                     GRID(ng) % umask,                            &
     &                     GRID(ng) % vmask,                            &
# endif
     &                     GRID(ng) % Hz,                               &
     &                     GRID(ng) % om_v,                             &
     &                     GRID(ng) % on_u,                             &
     &                     GRID(ng) % pm,                               &
     &                     GRID(ng) % pn,                               &
     &                     GRID(ng) % z_r,                              &
     &                     GRID(ng) % z_w,                              &
     &                     MIXING(ng) % Akv,                            &
     &                     COUPLING(ng) % DU_avg1,                      &
     &                     COUPLING(ng) % DV_avg1,                      &
     &                     COUPLING(ng) % DU_avg2,                      &
     &                     COUPLING(ng) % DV_avg2,                      &
     &                     OCEAN(ng) % ru,                              &
     &                     OCEAN(ng) % rv,                              &
     &                     OCEAN(ng) % u,                               &
     &                     OCEAN(ng) % v,                               &
     &                     OCEAN(ng) % ubar,                            &
     &                     OCEAN(ng) % vbar,                            &
     &                     GRID(ng) % Huon,                             &
     &                     GRID(ng) % Hvom)
# ifdef PROFILE
      CALL wclock_off (ng, 34)
# endif
      RETURN
      END SUBROUTINE step3d_uv
!
!***********************************************************************
      SUBROUTINE step3d_uv_tile (ng, Istr, Iend, Jstr, Jend,            &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           nrhs, nstp, nnew,                      &
# ifdef UV_PSOURCE
     &                           Nsrc,                                  &
     &                           Isrc, Jsrc, Dsrc, Qsrc,                &
# endif
# ifdef MASKING
     &                           umask, vmask,                          &
# endif
     &                           Hz, om_v, on_u, pm, pn, z_r, z_w,      &
     &                           Akv,                                   &
     &                           DU_avg1, DV_avg1, DU_avg2, DV_avg2,    &
     &                           ru, rv,                                &
     &                           u, v, ubar, vbar, Huon, Hvom)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      USE exchange_2d_mod
      USE exchange_3d_mod
# endif
      USE u3dbc_mod, ONLY : u3dbc_tile
      USE v3dbc_mod, ONLY : v3dbc_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: nrhs, nstp, nnew

# ifdef ASSUMED_SHAPE
#  ifdef UV_PSOURCE
      integer, intent(in) :: Nsrc
      integer, intent(in) :: Isrc(:)
      integer, intent(in) :: Jsrc(:)
      real(r8), intent(in) :: Dsrc(:)
      real(r8), intent(in) :: Qsrc(:,:)
#  endif
#  ifdef MASKING
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: Akv(LBi:,LBj:,0:)
      real(r8), intent(in) :: DU_avg1(LBi:,LBj:)
      real(r8), intent(in) :: DV_avg1(LBi:,LBj:)
      real(r8), intent(in) :: DU_avg2(LBi:,LBj:)
      real(r8), intent(in) :: DV_avg2(LBi:,LBj:)
      real(r8), intent(in) :: ru(LBi:,LBj:,0:,:)
      real(r8), intent(in) :: rv(LBi:,LBj:,0:,:)
      real(r8), intent(inout) :: u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: v(LBi:,LBj:,:,:)
      real(r8), intent(out) :: ubar(LBi:,LBj:,:)
      real(r8), intent(out) :: vbar(LBi:,LBj:,:)
      real(r8), intent(out) :: Huon(LBi:,LBj:,:)
      real(r8), intent(out) :: Hvom(LBi:,LBj:,:)
# else
#  ifdef UV_PSOURCE
      integer, intent(in) :: Nsrc
      integer, intent(in) :: Isrc(Msrc)
      integer, intent(in) :: Jsrc(Msrc)
      real(r8), intent(in) :: Dsrc(Msrc)
      real(r8), intent(in) :: Qsrc(Msrc,N(ng))
#  endif
#  ifdef MASKING
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: om_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(in) :: Akv(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(in) :: DU_avg1(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: DV_avg1(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: DU_avg2(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: DV_avg2(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: ru(LBi:UBi,LBj:UBj,0:N(ng),2)
      real(r8), intent(in) :: rv(LBi:UBi,LBj:UBj,0:N(ng),2)
      real(r8), intent(inout) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: v(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(out) :: ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(out) :: vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(out) :: Huon(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(out) :: Hvom(LBi:UBi,LBj:UBj,N(ng))
# endif
!
!  Imported variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, is, j, k

      real(r8) :: cff, cff1, cff2

      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N(ng)) :: AK
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N(ng)) :: BC
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N(ng)) :: CF
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N(ng)) :: DC
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N(ng)) :: FC
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N(ng)) :: Hzk
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N(ng)) :: oHz

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Time step momentum equation in the XI-direction.
!-----------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO i=IstrU,Iend
          AK(i,0)=0.5_r8*(Akv(i-1,j,0)+Akv(i,j,0))
          DO k=1,N(ng)
            AK(i,k)=0.5_r8*(Akv(i-1,j,k)+Akv(i,j,k))
            Hzk(i,k)=0.5_r8*(Hz(i-1,j,k)+Hz(i,j,k))
# ifdef SPLINES
            oHz(i,k)=1.0_r8/Hzk(i,k)
# endif
          END DO
          AK(i,N(ng))=0.5_r8*(Akv(i-1,j,N(ng))+Akv(i,j,N(ng)))
        END DO
!
!  Time step right-hand-side terms.
!
        IF (iic(ng).eq.ntfirst) THEN
          cff=0.25_r8*dt(ng)
        ELSE IF (iic(ng).eq.(ntfirst+1)) THEN
          cff=0.25_r8*dt(ng)*3.0_r8/2.0_r8
        ELSE
          cff=0.25_r8*dt(ng)*23.0_r8/12.0_r8
        END IF
        DO i=IstrU,Iend
          DC(i,0)=cff*(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j))
        END DO
        DO k=1,N(ng)
          DO i=IstrU,Iend
            u(i,j,k,nnew)=(u(i,j,k,nnew)+DC(i,0)*ru(i,j,k,nrhs))
# ifdef SPLINES
            u(i,j,k,nnew)=u(i,j,k,nnew)*oHz(i,k)
# endif
          END DO
        END DO
# ifdef SPLINES
!
!  Use conservative, parabolic spline reconstruction of vertical
!  viscosity derivatives.  Then, time step vertical viscosity term
!  implicitly.
!
        cff1=1.0_r8/6.0_r8
        DO k=1,N(ng)-1
          DO i=IstrU,Iend
            FC(i,k)=cff1*Hzk(i,k  )-dt(ng)*AK(i,k-1)*oHz(i,k  )
            CF(i,k)=cff1*Hzk(i,k+1)-dt(ng)*AK(i,k+1)*oHz(i,k+1)
          END DO
        END DO
        DO i=IstrU,Iend
          CF(i,0)=0.0_r8
          DC(i,0)=0.0_r8
        END DO
!
        cff1=1.0_r8/3.0_r8
        DO k=1,N(ng)-1
          DO i=IstrU,Iend
            BC(i,k)=cff1*(Hzk(i,k)+Hzk(i,k+1))+                         &
     &              dt(ng)*AK(i,k)*(oHz(i,k)+oHz(i,k+1))
            cff=1.0_r8/(BC(i,k)-FC(i,k)*CF(i,k-1))
            CF(i,k)=cff*CF(i,k)
            DC(i,k)=cff*(u(i,j,k+1,nnew)-u(i,j,k,nnew)-                 &
     &                   FC(i,k)*DC(i,k-1))
          END DO
        END DO
!
        DO i=IstrU,Iend
          DC(i,N(ng))=0.0_r8
        END DO
        DO k=N(ng)-1,1,-1
          DO i=IstrU,Iend
            DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1)
          END DO
        END DO
!
        DO k=1,N(ng)
          DO i=IstrU,Iend
            DC(i,k)=DC(i,k)*AK(i,k)
            u(i,j,k,nnew)=u(i,j,k,nnew)+                                &
     &                    dt(ng)*oHz(i,k)*(DC(i,k)-DC(i,k-1))
          END DO
        END DO
# else
!
!  Compute off-diagonal coefficients [lambda*dt*Akv/Hz] for the
!  implicit vertical viscosity term at horizontal U-points and
!  vertical W-points.
!
        cff=-lambda*dt(ng)/0.5_r8
        DO k=1,N(ng)-1
          DO i=IstrU,Iend
            FC(i,k)=cff*AK(i,k)/(z_r(i,j,k+1)+z_r(i-1,j,k+1)-           &
     &                           z_r(i,j,k  )-z_r(i-1,j,k  ))
          END DO
        END DO
        DO i=IstrU,Iend
          FC(i,0)=0.0_r8
          FC(i,N(ng))=0.0_r8
        END DO
!
!  Solve the tridiagonal system.
!
        DO k=1,N(ng)
          DO i=IstrU,Iend
            DC(i,k)=u(i,j,k,nnew)
            BC(i,k)=Hzk(i,k)-FC(i,k)-FC(i,k-1)
          END DO
        END DO
        DO i=IstrU,Iend
          cff=1.0_r8/BC(i,1)
          CF(i,1)=cff*FC(i,1)
          DC(i,1)=cff*DC(i,1)
        END DO
        DO k=2,N(ng)-1
          DO i=IstrU,Iend
            cff=1.0_r8/(BC(i,k)-FC(i,k-1)*CF(i,k-1))
            CF(i,k)=cff*FC(i,k)
            DC(i,k)=cff*(DC(i,k)-FC(i,k-1)*DC(i,k-1))
          END DO
        END DO
!
!  Compute new solution by back subtitution.
!
        DO i=IstrU,Iend
          DC(i,N(ng))=(DC(i,N(ng))-FC(i,N(ng)-1)*DC(i,N(ng)-1))/        &
     &                (BC(i,N(ng))-FC(i,N(ng)-1)*CF(i,N(ng)-1))
          u(i,j,N(ng),nnew)=DC(i,N(ng))
        END DO
        DO k=N(ng)-1,1,-1
          DO i=IstrU,Iend
            DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1)
            u(i,j,k,nnew)=DC(i,k)
          END DO
        END DO
# endif
!
!  Replace incorrect vertical mean with more accurate barotropic
!  component, ubar=DU_avg1/(D*on_u).  Recall that, D=CF(:,0).
!
        DO i=IstrU,Iend
          CF(i,0)=Hzk(i,1)
          DC(i,0)=u(i,j,1,nnew)*Hzk(i,1)
        END DO
        DO k=2,N(ng)
          DO i=IstrU,Iend
            CF(i,0)=CF(i,0)+Hzk(i,k)
            DC(i,0)=DC(i,0)+u(i,j,k,nnew)*Hzk(i,k)
          END DO
        END DO
        DO i=IstrU,Iend
          DC(i,0)=(DC(i,0)*on_u(i,j)-DU_avg1(i,j))/                     &
     &            (CF(i,0)*on_u(i,j))
        END DO
!
!  Couple and update new solution.
!
        DO k=1,N(ng)
          DO i=IstrU,Iend
            u(i,j,k,nnew)=u(i,j,k,nnew)-DC(i,0)
# ifdef MASKING
            u(i,j,k,nnew)=u(i,j,k,nnew)*umask(i,j)
# endif
          END DO
        END DO
!
!-----------------------------------------------------------------------
!  Time step momentum equation in the ETA-direction.
!-----------------------------------------------------------------------
!
        IF (j.ge.JstrV) THEN
          DO i=Istr,Iend
            AK(i,0)=0.5_r8*(Akv(i,j-1,0)+Akv(i,j,0))
            DO k=1,N(ng)
              AK(i,k)=0.5_r8*(Akv(i,j-1,k)+Akv(i,j,k))
              Hzk(i,k)=0.5_r8*(Hz(i,j-1,k)+Hz(i,j,k))
# ifdef SPLINES
              oHz(i,k)=1.0_r8/Hzk(i,k)
# endif
            END DO
            AK(i,N(ng))=0.5_r8*(Akv(i,j-1,N(ng))+Akv(i,j,N(ng)))
          END DO
!
!  Time step right-hand-side terms.
!
          IF (iic(ng).eq.ntfirst) THEN
            cff=0.25_r8*dt(ng)
          ELSE IF (iic(ng).eq.(ntfirst+1)) THEN
            cff=0.25_r8*dt(ng)*3.0_r8/2.0_r8
          ELSE
            cff=0.25_r8*dt(ng)*23.0_r8/12.0_r8
          END IF
          DO i=Istr,Iend
            DC(i,0)=cff*(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))
          END DO
          DO k=1,N(ng)
            DO i=Istr,Iend
              v(i,j,k,nnew)=(v(i,j,k,nnew)+DC(i,0)*rv(i,j,k,nrhs))
# ifdef SPLINES
              v(i,j,k,nnew)=v(i,j,k,nnew)*oHz(i,k)
# endif
            END DO
          END DO
# ifdef SPLINES
!
!  Use conservative, parabolic spline reconstruction of vertical
!  viscosity derivatives.  Then, time step vertical viscosity term
!  implicitly.
!
          cff1=1.0_r8/6.0_r8
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              FC(i,k)=cff1*Hzk(i,k  )-dt(ng)*AK(i,k-1)*oHz(i,k  )
              CF(i,k)=cff1*Hzk(i,k+1)-dt(ng)*AK(i,k+1)*oHz(i,k+1)
            END DO
          END DO
          DO i=Istr,Iend
            CF(i,0)=0.0_r8
            DC(i,0)=0.0_r8
          END DO
!
          cff1=1.0_r8/3.0_r8
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              BC(i,k)=cff1*(Hzk(i,k)+Hzk(i,k+1))+                       &
     &                dt(ng)*AK(i,k)*(oHz(i,k)+oHz(i,k+1))
              cff=1.0_r8/(BC(i,k)-FC(i,k)*CF(i,k-1))
              CF(i,k)=cff*CF(i,k)
              DC(i,k)=cff*(v(i,j,k+1,nnew)-v(i,j,k,nnew)-               &
     &                     FC(i,k)*DC(i,k-1))
            END DO
          END DO
!
          DO i=Istr,Iend
            DC(i,N(ng))=0.0_r8
          END DO
          DO k=N(ng)-1,1,-1
            DO i=Istr,Iend
              DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1)
            END DO
          END DO
!
          DO k=1,N(ng)
            DO i=Istr,Iend
              DC(i,k)=DC(i,k)*AK(i,k)
              v(i,j,k,nnew)=v(i,j,k,nnew)+                              &
     &                      dt(ng)*oHz(i,k)*(DC(i,k)-DC(i,k-1))
            END DO
          END DO
# else
!
!  Compute off-diagonal coefficients [lambda*dt*Akv/Hz] for the
!  implicit vertical viscosity term at horizontal V-points and
!  vertical W-points.
!
          cff=-lambda*dt(ng)/0.5_r8
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              FC(i,k)=cff*AK(i,k)/(z_r(i,j,k+1)+z_r(i,j-1,k+1)-         &
     &                             z_r(i,j,k  )-z_r(i,j-1,k  ))
            END DO
          END DO
          DO i=Istr,Iend
            FC(i,0)=0.0_r8
            FC(i,N(ng))=0.0_r8
          END DO
!
!  Solve the tridiagonal system.
!
          DO k=1,N(ng)
            DO i=Istr,Iend
              DC(i,k)=v(i,j,k,nnew)
              BC(i,k)=Hzk(i,k)-FC(i,k)-FC(i,k-1)
            END DO
          END DO
          DO i=Istr,Iend
            cff=1.0_r8/BC(i,1)
            CF(i,1)=cff*FC(i,1)
            DC(i,1)=cff*DC(i,1)
          END DO
          DO k=2,N(ng)-1
            DO i=Istr,Iend
              cff=1.0_r8/(BC(i,k)-FC(i,k-1)*CF(i,k-1))
              CF(i,k)=cff*FC(i,k)
              DC(i,k)=cff*(DC(i,k)-FC(i,k-1)*DC(i,k-1))
            END DO
          END DO
!
!  Compute new solution by back subtitution.
!
          DO i=Istr,Iend
            DC(i,N(ng))=(DC(i,N(ng))-FC(i,N(ng)-1)*DC(i,N(ng)-1))/      &
     &                  (BC(i,N(ng))-FC(i,N(ng)-1)*CF(i,N(ng)-1))
            v(i,j,N(ng),nnew)=DC(i,N(ng))
          END DO
          DO k=N(ng)-1,1,-1
            DO i=Istr,Iend
              DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1)
              v(i,j,k,nnew)=DC(i,k)
            END DO
          END DO
# endif
!
!  Replace incorrect vertical mean with more accurate barotropic
!  component, vbar=DV_avg1/(D*om_v).  Recall that, D=CF(:,0).
!
          DO i=Istr,Iend
            CF(i,0)=Hzk(i,1)
            DC(i,0)=v(i,j,1,nnew)*Hzk(i,1)
          END DO
          DO k=2,N(ng)
            DO i=Istr,Iend
              CF(i,0)=CF(i,0)+Hzk(i,k)
              DC(i,0)=DC(i,0)+v(i,j,k,nnew)*Hzk(i,k)
            END DO
          END DO
          DO i=Istr,Iend
            DC(i,0)=(DC(i,0)*om_v(i,j)-DV_avg1(i,j))/                   &
     &              (CF(i,0)*om_v(i,j))
          END DO
!
!  Couple and update new solution.
!
          DO k=1,N(ng)
            DO i=Istr,Iend
              v(i,j,k,nnew)=v(i,j,k,nnew)-DC(i,0)
# ifdef MASKING
              v(i,j,k,nnew)=v(i,j,k,nnew)*vmask(i,j)
# endif
            END DO
          END DO
        END IF
      END DO
!
!-----------------------------------------------------------------------
! Set lateral boundary conditions.
!-----------------------------------------------------------------------
!
      CALL u3dbc_tile (ng, Istr, Iend, Jstr, Jend,                      &
     &                 LBi, UBi, LBj, UBj, N(ng),                       &
     &                 nstp, nnew,                                      &
     &                 u)
      CALL v3dbc_tile (ng, Istr, Iend, Jstr, Jend,                      &
     &                 LBi, UBi, LBj, UBj, N(ng),                       &
     &                 nstp, nnew,                                      &
     &                 v)
# ifdef UV_PSOURCE
!
!-----------------------------------------------------------------------
! Apply mass point sources.
!-----------------------------------------------------------------------
!
      DO is=1,Nsrc
        i=Isrc(is)
        j=Jsrc(is)
        IF (INT(Dsrc(is)).eq.0) THEN
          IF (((Istr.le.i).and.(i.le.Iend+1)).and.                      &
     &        ((Jstr.le.j).and.(j.le.Jend))) THEN
            DO k=1,N(ng)
              u(i,j,k,nnew)=Qsrc(is,k)/                                 &
     &                      (on_u(i,j)*                                 &
     &                       0.5_r8*(z_w(i-1,j,k)-z_w(i-1,j,k-1)+       &
     &                               z_w(i  ,j,k)-z_w(i  ,j,k-1)))
            END DO
          END IF
        ELSE
          IF (((Istr.le.i).and.(i.le.Iend)).and.                        &
     &        ((Jstr.le.j).and.(j.le.Jend+1))) THEN
            DO k=1,N(ng)
                v(i,j,k,nnew)=Qsrc(is,k)/                               &
     &                        (om_v(i,j)*                               &
     &                         0.5_r8*(z_w(i,j-1,k)-z_w(i,j-1,k-1)+     &
     &                                 z_w(i,j  ,k)-z_w(i,j  ,k-1)))
            END DO
          END IF
        END IF
      END DO
# endif
!
!-----------------------------------------------------------------------
!  Couple 2D and 3D momentum equations.
!-----------------------------------------------------------------------
!
# ifdef EW_PERIODIC
#  define IU_RANGE Istr,Iend
#  define IV_RANGE Istr,Iend
# else
#  define IU_RANGE Istr,IendR
#  define IV_RANGE IstrR,IendR
# endif
# ifdef NS_PERIODIC
#  define JU_RANGE Jstr,Jend
#  define JV_RANGE Jstr,Jend
# else
#  define JU_RANGE JstrR,JendR
#  define JV_RANGE Jstr,JendR
# endif
!
!  Couple velocity component in the XI-direction.
!
      DO j=JU_RANGE
        DO i=IU_RANGE
          DC(i,0)=0.0_r8
          CF(i,0)=0.0_r8
          FC(i,0)=0.0_r8
        END DO
!
!  Compute thicknesses of U-boxes DC(i,1:N), total depth of the water
!  column DC(i,0), and incorrect vertical mean CF(i,0).  Notice that
!  barotropic component is replaced with its fast-time averaged
!  values.
!
        DO k=1,N(ng)
          DO i=IU_RANGE
            DC(i,k)=0.5_r8*(Hz(i,j,k)+Hz(i-1,j,k))*on_u(i,j)
            DC(i,0)=DC(i,0)+DC(i,k)
            CF(i,0)=CF(i,0)+DC(i,k)*u(i,j,k,nnew)
          END DO
        END DO
        DO i=IU_RANGE
          DC(i,0)=1.0_r8/DC(i,0)
          CF(i,0)=DC(i,0)*(CF(i,0)-DU_avg1(i,j))
          ubar(i,j,1)=DC(i,0)*DU_avg1(i,j)
          ubar(i,j,2)=ubar(i,j,1)
        END DO
!
!  Replace incorrect vertical mean with more accurate barotropic
!  component, ubar=DU_avg1/(D*on_u).  Recall that, D=CF(:,0).
!
        DO k=N(ng),1,-1
          DO i=IU_RANGE
            u(i,j,k,nnew)=u(i,j,k,nnew)-CF(i,0)
!!# ifdef MASKING
!!          u(i,j,k,nnew)=u(i,j,k,nnew)*umask(i,j)
!!# endif
            Huon(i,j,k)=0.5_r8*(Huon(i,j,k)+u(i,j,k,nnew)*DC(i,k))
            FC(i,0)=FC(i,0)+Huon(i,j,k)
          END DO
        END DO
!
!  Compute correct mass flux, Hz*u/n.
!
        DO i=IU_RANGE
          FC(i,0)=DC(i,0)*(FC(i,0)-DU_avg2(i,j))
        END DO
        DO k=1,N(ng)
          DO i=IU_RANGE
            Huon(i,j,k)=Huon(i,j,k)-DC(i,k)*FC(i,0)
          END DO
        END DO
!
!  Couple velocity component in the ETA-direction.
!
        IF (j.ge.Jstr) THEN
          DO i=IV_RANGE
            DC(i,0)=0.0_r8
            CF(i,0)=0.0_r8
            FC(i,0)=0.0_r8
          END DO
!
!  Compute thicknesses of V-boxes DC(i,1:N), total depth of the water
!  column DC(i,0), and incorrect vertical mean CF(i,0).  Notice that
!  barotropic component is replaced with its fast-time averaged
!  values.
!
          DO k=1,N(ng)
            DO i=IV_RANGE
              DC(i,k)=0.5_r8*(Hz(i,j,k)+Hz(i,j-1,k))*om_v(i,j)
              DC(i,0)=DC(i,0)+DC(i,k)
              CF(i,0)=CF(i,0)+DC(i,k)*v(i,j,k,nnew)
            END DO
          END DO
          DO i=IV_RANGE
            DC(i,0)=1.0_r8/DC(i,0)
            CF(i,0)=DC(i,0)*(CF(i,0)-DV_avg1(i,j))
            vbar(i,j,1)=DC(i,0)*DV_avg1(i,j)
            vbar(i,j,2)=vbar(i,j,1)
          END DO
!
!  Replace incorrect vertical mean with more accurate barotropic
!  component, vbar=DV_avg1/(D*om_v).  Recall that, D=CF(:,0).
!
          DO k=N(ng),1,-1
            DO i=IV_RANGE
              v(i,j,k,nnew)=v(i,j,k,nnew)-CF(i,0)
!!# ifdef MASKING
!!            v(i,j,k,nnew)=v(i,j,k,nnew)*vmask(i,j)
!!# endif
              Hvom(i,j,k)=0.5_r8*(Hvom(i,j,k)+v(i,j,k,nnew)*DC(i,k))
              FC(i,0)=FC(i,0)+Hvom(i,j,k)
            END DO
          END DO
!
!  Compute correct mass flux, Hz*v/m.
!
          DO i=IV_RANGE
            FC(i,0)=DC(i,0)*(FC(i,0)-DV_avg2(i,j))
          END DO
          DO k=1,N(ng)
            DO i=IV_RANGE
              Hvom(i,j,k)=Hvom(i,j,k)-DC(i,k)*FC(i,0)
            END DO
          END DO
        END IF
      END DO
# undef IU_RANGE
# undef JU_RANGE
# undef IV_RANGE
# undef JV_RANGE
!
# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
!
!  Apply periodic boundaries conditions to mass fluxes and barotropic
!  velocity components, if so prescribed.
!
      CALL exchange_u3d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, u(:,:,:,nnew))
      CALL exchange_v3d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, v(:,:,:,nnew))
      CALL exchange_u3d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, Huon)
      CALL exchange_v3d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, Hvom)
!
      DO k=1,2
        CALL exchange_u2d_tile (ng, Istr, Iend, Jstr, Jend,             &
     &                          LBi, UBi, LBj, UBj, ubar(:,:,k))
        CALL exchange_v2d_tile (ng, Istr, Iend, Jstr, Jend,             &
     &                          LBi, UBi, LBj, UBj, vbar(:,:,k))
      END DO
# endif
      RETURN
      END SUBROUTINE step3d_uv_tile
#endif
      END MODULE step3d_uv_mod
