#include "cppdefs.h"
#define SASHA
      MODULE lmd_bkpp_mod
#ifdef LMD_BKPP
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine determines the depth of  bottom  oceanic boundary      !
!  layer,  hbbl,  as the deepest depth  where the bulk Richardson      !
!  number is equal to the critical value, Ric.                         !
!                                                                      !
!  Then,  it computes the vertical mixing coefficients  within the     !
!  boundary layer. They depend on surface forcing and the magnitude    !
!  and gradient of interior mixing below  the boundary layer.  The     !
!  ocean interior is allowed to force the boundary layer through a     !
!  dependence of the nondimensional vertical shape function G(sigma)   !
!  and its vertical derivative at  sigma=1  on the interior  mixing    !
!  coefficients, and it vertical derivative at d=hsbl. The boundary    !
!  layer mixing coefficients are computed by matching these values.    !
!                                                                      !
! Reference:                                                           !
!                                                                      !
!  Large, W.G., J.C. McWilliams, and S.C. Doney, 1994: A Review        !
!    and model with a nonlocal boundary layer parameterization,        !
!    Reviews of Geophysics, 32,363-403.                                !
!                                                                      !
!  This routine was adapted from Bill Large 1995 code.                 !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: lmd_bkpp

      CONTAINS
!
!***********************************************************************
      SUBROUTINE lmd_bkpp (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_forces
      USE mod_grid
      USE mod_mixing
      USE mod_ocean
      USE mod_stepping
!
      implicit none

      integer, intent(in) :: ng, tile

# include "tile.h"
!
      CALL lmd_bkpp_tile (ng, Istr, Iend, Jstr, Jend,                   &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    nstp(ng),                                     &
# ifdef MASKING
     &                    GRID(ng) % rmask,                             &
# endif
     &                    GRID(ng) % f,                                 &
     &                    GRID(ng) % h,                                 &
     &                    GRID(ng) % Hz,                                &
     &                    GRID(ng) % z_r,                               &
     &                    GRID(ng) % z_w,                               &
     &                    OCEAN(ng) % u,                                &
     &                    OCEAN(ng) % v,                                &
     &                    OCEAN(ng) % pden,                             &
     &                    FORCES(ng) % srflx,                           &
     &                    FORCES(ng) % stflx,                           &
     &                    FORCES(ng) % bustr,                           &
     &                    FORCES(ng) % bvstr,                           &
     &                    MIXING(ng) % alpha,                           &
# ifdef SALINITY
     &                    MIXING(ng) % beta,                            &
# endif
     &                    MIXING(ng) % bvf,                             &
     &                    MIXING(ng) % ksbl,                            &
     &                    MIXING(ng) % Akt,                             &
     &                    MIXING(ng) % Akv,                             &
     &                    MIXING(ng) % kbbl,                            &
     &                    MIXING(ng) % hbbl)
      RETURN
      END SUBROUTINE lmd_bkpp
!
!***********************************************************************
      SUBROUTINE lmd_bkpp_tile (ng, Istr, Iend, Jstr, Jend,             &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          nstp,                                   &
# ifdef MASKING
     &                          rmask,                                  &
# endif
     &                          f, h, Hz, z_r, z_w,                     &
     &                          u, v, pden,                             &
     &                          srflx, stflx,  bustr, bvstr,            &
     &                          alpha,                                  &
# ifdef SALINITY
     &                          beta,                                   &
# endif
     &                          bvf,                                    &
     &                          ksbl, Akt, Akv, kbbl, hbbl)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE bc_2d_mod, ONLY : bc_r2d_tile
      USE bc_3d_mod, ONLY : bc_w3d_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, Intent(in) :: nstp

# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: f(LBi:,LBj:)
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
      real(r8), intent(in) :: pden(LBi:,LBj:,:)
      real(r8), intent(in) :: srflx(LBi:,LBj:)
      real(r8), intent(in) :: stflx(LBi:,LBj:,:)
      real(r8), intent(in) :: bustr(LBi:,LBj:)
      real(r8), intent(in) :: bvstr(LBi:,LBj:)
      real(r8), intent(in) :: alpha(LBi:,LBj:)
#  ifdef SALINITY
      real(r8), intent(in) :: beta(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: bvf(LBi:,LBj:,0:)
      integer,  intent(in) :: ksbl(LBi:,LBj:)
      real(r8), intent(inout) :: Akt(LBi:,LBj:,0:,:)
      real(r8), intent(inout) :: Akv(LBi:,LBj:,0:)
      real(r8), intent(inout) :: hbbl(LBi:,LBj:)
      integer,  intent(out) :: kbbl(LBi:,LBj:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: f(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),3)
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),3)
      real(r8), intent(in) :: pden(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: srflx(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: stflx(LBi:UBi,LBj:UBj,NT(ng))
      real(r8), intent(in) :: bustr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: bvstr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: alpha(LBi:UBi,LBj:UBj)
#  ifdef SALINITY
      real(r8), intent(in) :: beta(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: bvf(LBi:UBi,LBj:UBj,0:N(ng))
      integer,  intent(in) :: ksbl(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: Akt(LBi:UBi,LBj:UBj,0:N(ng),NAT)
      real(r8), intent(inout) :: Akv(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(inout) :: hbbl(LBi:UBi,LBj:UBj)
      integer,  intent(out) :: kbbl(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, itrc, j, k

      real(r8), parameter :: eps = 1.0E-10_r8

      real(r8) :: Gm, Gt, Gs, K_bl, Ribot, Ritop, Rk, Rref
      real(r8) :: Uk, Uref, Vk, Vref, Vtc
      real(r8) :: a1, a2, a3, cff, cff1, cff2, cff_up, cff_dn
      real(r8) :: depth, dK_bl, hekman, sigma, zbl

      real(r8), dimension (PRIVATE_2D_SCRATCH_ARRAY,0:N(ng)) :: Bflux
      real(r8), dimension (PRIVATE_2D_SCRATCH_ARRAY,0:N(ng)) :: Kv
      real(r8), dimension (PRIVATE_2D_SCRATCH_ARRAY,0:N(ng)) :: Kt
# ifdef SALINITY
      real(r8), dimension (PRIVATE_2D_SCRATCH_ARRAY,0:N(ng)) :: Ks
# endif

      real(r8), dimension (PRIVATE_1D_SCRATCH_ARRAY,0:N(ng)) :: FC
      real(r8), dimension (PRIVATE_1D_SCRATCH_ARRAY,0:N(ng)) :: dR
      real(r8), dimension (PRIVATE_1D_SCRATCH_ARRAY,0:N(ng)) :: dU
      real(r8), dimension (PRIVATE_1D_SCRATCH_ARRAY,0:N(ng)) :: dV

      real(r8), dimension (PRIVATE_2D_SCRATCH_ARRAY) :: Bo
      real(r8), dimension (PRIVATE_2D_SCRATCH_ARRAY) :: Bosol
      real(r8), dimension (PRIVATE_2D_SCRATCH_ARRAY) :: Bfbot
      real(r8), dimension (PRIVATE_2D_SCRATCH_ARRAY) :: Gm1
      real(r8), dimension (PRIVATE_2D_SCRATCH_ARRAY) :: Gt1
      real(r8), dimension (PRIVATE_2D_SCRATCH_ARRAY) :: Gs1
      real(r8), dimension (PRIVATE_2D_SCRATCH_ARRAY) :: Ustar
      real(r8), dimension (PRIVATE_2D_SCRATCH_ARRAY) :: bl_dpth
      real(r8), dimension (PRIVATE_2D_SCRATCH_ARRAY) :: dGm1dS
      real(r8), dimension (PRIVATE_2D_SCRATCH_ARRAY) :: dGt1dS
      real(r8), dimension (PRIVATE_2D_SCRATCH_ARRAY) :: dGs1dS
      real(r8), dimension (PRIVATE_2D_SCRATCH_ARRAY) :: f1
      real(r8), dimension (PRIVATE_2D_SCRATCH_ARRAY) :: swdk
      real(r8), dimension (PRIVATE_2D_SCRATCH_ARRAY) :: wm
      real(r8), dimension (PRIVATE_2D_SCRATCH_ARRAY) :: ws
      real(r8), dimension (PRIVATE_2D_SCRATCH_ARRAY) :: zgrid

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Initialize relevant parameters.
!-----------------------------------------------------------------------
!
      Vtc=lmd_Cv*SQRT(-lmd_betaT)/(SQRT(lmd_cs*lmd_epsilon)*            &
     &                             lmd_Ric*vonKar*vonKar)
!
!-----------------------------------------------------------------------
!  Load current vertical viscosity and diffusivities into work
!  arrays.
!-----------------------------------------------------------------------
!
      DO k=0,N(ng)
        DO j=Jstr,Jend
          DO i=Istr,Iend
            Kv(i,j,k)=Akv(i,j,k)
            Kt(i,j,k)=Akt(i,j,k,itemp)
# ifdef SALINITY
            Ks(i,j,k)=Akt(i,j,k,isalt)
# endif
          END DO
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Get approximation of bottom layer depth using "lmd_eps" and
!  boundary layer depth from previous time step.
!-----------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          bl_dpth(i,j)=lmd_epsilon*(hbbl(i,j)-z_w(i,j,0))
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Compute turbulent friction velocity (m/s) "Ustar" from bottom
!  stress at RHO-points.
!-----------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          Ustar(i,j)=SQRT(SQRT((0.5_r8*(bustr(i,j)+bustr(i+1,j)))**2+   &
     &                         (0.5_r8*(bvstr(i,j)+bvstr(i,j+1)))**2))
# ifdef MASKING
          Ustar(i,j)=Ustar(i,j)*rmask(i,j)
# endif
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Compute bottom turbulent buoyancy forcing "Bo" (m2/s3). Compute
!  surface radiative buoyancy forcing "Bosol" (m2/s3) that can be
!  used in shallow areas when it can penetrate all the way to the
!  bottom.
!-----------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
# ifdef SALINITY
          Bo(i,j)=g*(alpha(i,j)*(stflx(i,j,itemp)-srflx(i,j))-          &
     &               beta (i,j)*stflx(i,j,isalt))
# else
          Bo(i,j)=g*alpha(i,j)*(stflx(i,j,itemp)-srflx(i,j))
# endif
!!        Bosol(i,j)=g*alpha(i,j)*srflx(i,j)
          Bosol(i,j)=0.0_r8
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Compute total buoyancy flux (m2/s3) at W-points.
!-----------------------------------------------------------------------
!
      DO k=0,N(ng)
        DO j=Jstr,Jend
          DO i=Istr,Iend
            zgrid(i,j)=z_w(i,j,N(ng))-z_w(i,j,k)
          END DO
        END DO
        CALL lmd_swfrac_tile (ng, Istr, Iend, Jstr, Jend, -1.0_r8,      &
     &                        zgrid, swdk)
        DO j=Jstr,Jend
          DO i=Istr,Iend
            Bflux(i,j,k)=(Bo(i,j)+Bosol(i,j)*(1.0_r8-swdk(i,j)))
# ifdef MASKING
            Bflux(i,j,k)=Bflux(i,j,k)*rmask(i,j)
# endif
          END DO
        END DO
      END DO
!
!=======================================================================
!  Compute bulk Richardson number "Rib" and then find depth of the
!  oceanic bottom boundary layer "hbbl", such that Rib(hbbl)=Ric.
!=======================================================================
!
      DO j=Jstr,Jend
# ifdef SPLINES
!
! Construct parabolic splines for vertical derivatives of potential
! density and velocity components at W-points.  FC is a scratch array.
!
        DO i=Istr,Iend
          FC(i,0)=0.0_r8
          dR(i,0)=0.0_r8
          dU(i,0)=0.0_r8
          dV(i,0)=0.0_r8
        END DO
        DO k=1,N(ng)-1
          DO i=Istr,Iend
            cff=1.0_r8/(2.0_r8*Hz(i,j,k+1)+                             &
     &                  Hz(i,j,k)*(2.0_r8-FC(i,k-1)))
            FC(i,k)=cff*Hz(i,j,k+1)
            dR(i,k)=cff*(6.0_r8*(pden(i,j,k+1)-pden(i,j,k))-            &
     &                   Hz(i,j,k)*dR(i,k-1))
            dU(i,k)=cff*(3.0_r8*(u(i  ,j,k+1,nstp)-u(i,  j,k,nstp)+     &
     &                           u(i+1,j,k+1,nstp)-u(i+1,j,k,nstp))-    &
     &                   Hz(i,j,k)*dU(i,k-1))
            dV(i,k)=cff*(3.0_r8*(v(i,j  ,k+1,nstp)-v(i,j  ,k,nstp)+     &
     &                           v(i,j+1,k+1,nstp)-v(i,j+1,k,nstp))-    &
     &                   Hz(i,j,k)*dV(i,k-1))
          END DO
        END DO
        DO i=Istr,Iend
          dR(i,N(ng))=0.0_r8
          dU(i,N(ng))=0.0_r8
          dV(i,N(ng))=0.0_r8
        END DO
        DO k=N(ng)-1,1,-1
          DO i=Istr,Iend
            dR(i,k)=dR(i,k)-FC(i,k)*dR(i,k+1)
            dU(i,k)=dU(i,k)-FC(i,k)*dU(i,k+1)
            dV(i,k)=dV(i,k)-FC(i,k)*dV(i,k+1)
          END DO
        END DO
# else
!
! Compute vertical derivatives of potential density and velocity
! components at W-points.
!
        DO k=1,N(ng)-1
          DO i=Istr,Iend
            cff=1.0_r8/(z_r(i,j,k+1)-z_r(i,j,k))
            dR(i,k)=cff*(pden(i,j,k+1)-pden(i,j,k))
            cff=0.5_r8*cff
            dU(i,k)=cff*(u(i  ,j,k+1,nstp)-u(i,  j,k,nstp)+             &
     &                   u(i+1,j,k+1,nstp)-u(i+1,j,k,nstp))
            dV(i,k)=cff*(v(i,j  ,k+1,nstp)-v(i,j  ,k,nstp)+             &
     &                   v(i,j+1,k+1,nstp)-v(i,j+1,k,nstp))
          END DO
        END DO
        DO i=Istr,Iend
          dR(i,0)=0.0_r8
          dR(i,N(ng))=0.0_r8
          dU(i,0)=0.0_r8
          dU(i,N(ng))=0.0_r8
          dV(i,0)=0.0_r8
          dV(i,N(ng))=0.0_r8
        END DO
# endif
!
!-----------------------------------------------------------------------
!  Compute bulk Richardson number "Rib" and then find depth of oceanic
!  bottom boundary layer "hbbl".
!
!                  [Br - B(d)] * d
!     Rib(d) = ----------------------- ;     Rib(hbbl)=Ric     (1)
!              |Vr - V(d)|^2 + Vt(d)^2
!
!  where "Br" and "Vr" are the bottom reference buoyancy and velocity
!  while "B(d)" and "V(d)" are the buoyancy and velocity at depth "d".
!
!  In the code below, the criterion "Rib(hbbl)=Ric" is reformulated
!  as follows:
!
!     Rib(d)       Ritop(d)
!     ------ = --------------- = 1                             (2)
!      Ric      Ric * Ribot(d)
!
!  where "Ritop" and "Ribot" are numerator and denominator in Eq. (1).
!  In its turn, Eq. (2) is rewritten in the following form:
!
!     FC(d) = Ritop(d) - Ric * Ribot(d) = 0                    (3)
!
!  That is, the planetary boundary layer extends to the depth where
!  the critical function "FC(d)" changes its sign.
!-----------------------------------------------------------------------
!
!  Compute potential density and velocity components bottom reference
!  values.
!
        cff1=1.0_r8/3.0_r8
        cff2=1.0_r8/6.0_r8
        DO i=Istr,Iend
          Rref=pden(i,j,1)-                                             &
     &         Hz(i,j,1)*(cff1*dR(i,0)+cff2*dR(i,1))
          Uref=0.5_r8*(u(i,j,1,nstp)+u(i+1,j,1,nstp))-                  &
     &         Hz(i,j,1)*(cff1*dU(i,0)+cff2*dU(i,1))
          Vref=0.5_r8*(v(i,j,1,nstp)+v(i,j+1,1,nstp))-                  &
     &         Hz(i,j,1)*(cff1*dV(i,0)+cff2*dV(i,1))
!
!  Compute critical function, FC, for bulk Richardson number.
!
          FC(i,0)=0.0_r8
          DO k=1,N(ng)
            depth=z_w(i,j,k)-z_w(i,j,0)
            IF (Bflux(i,j,k).lt.0.0_r8) THEN
              sigma=MIN(bl_dpth(i,j),depth)
            ELSE
              sigma=depth
            END IF
            CALL lmd_wscale (Bflux(i,j,k), Ustar(i,j), sigma,           &
     &                       wm(i,j), ws(i,j))
            Rk=pden(i,j,k)+                                             &
     &         Hz(i,j,k)*(cff1*dR(i,k)+cff2*dR(i,k-1))
            Uk=0.5_r8*(u(i,j,k,nstp)+u(i+1,j,k,nstp))+                  &
     &         Hz(i,j,k)*(cff1*dU(i,k)+cff2*dU(i,k-1))
            Vk=0.5_r8*(v(i,j,k,nstp)+v(i,j+1,k,nstp))+                  &
     &         Hz(i,j,k)*(cff1*dV(i,k)+cff2*dV(i,k-1))
!
            Ritop=-gorho0*(Rk-Rref)*depth
            Ribot=(Uk-Uref)**2+(Vk-Vref)**2+                            &
     &            Vtc*depth*ws(i,j)*SQRT(ABS(bvf(i,j,k)))
# ifdef SASHA
            FC(i,k)=Ritop-lmd_Ric*Ribot
# else
            FC(i,k)=Ritop/(Ribot+eps)
# endif
          END DO
        END DO
!
! Linearly interpolate to find "hbbl" where Rib/Ric=1.
!
        DO i=Istr,Iend
          kbbl(i,j)=N(ng)
          hbbl(i,j)=z_w(i,j,N(ng))
        END DO
# ifdef SASHA
        DO k=1,N(ng)-1
          DO i=Istr,Iend
            IF ((kbbl(i,j).eq.N(ng)).and.(FC(i,k).gt.0.0_r8)) THEN
              hbbl(i,j)=(z_w(i,j,k)*FC(i,k-1)-z_w(i,j,k-1)*FC(i,k))/    &
     &                  (FC(i,k-1)-FC(i,k))
              kbbl(i,j)=k
            END IF
          END DO
        END DO
# else
        DO k=1,N(ng)
          DO i=Istr,Iend
            IF ((kbbl(i,j).eq.N(ng)).and.((FC(i,k-1).lt.lmd_Ric).and.   &
     &                                    (FC(i,k  ).ge.lmd_Ric))) THEN
              hbbl(i,j)=((lmd_Ric-FC(i,k-1))*z_w(i,j,k  )+              &
     &                   (FC(i,k  )-lmd_Ric)*z_w(i,j,k-1))/             &
     &                  (FC(i,k)-FC(i,k-1))
              kbbl(i,j)=k
            END IF
          END DO
        END DO
# endif
      END DO
!
!  Compute total buoyancy flux at bottom boundary layer depth,
!  "Bfbot".
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          zgrid(i,j)=z_w(i,j,N(ng))-hbbl(i,j)
        END DO
      END DO
      CALL lmd_swfrac_tile (ng, Istr, Iend, Jstr, Jend, -1.0_r8,        &
     &                      zgrid, swdk)
      DO j=Jstr,Jend
        DO i=Istr,Iend
          Bfbot(i,j)=(Bo(i,j)+Bosol(i,j)*(1.0_r8-swdk(i,j)))
# ifdef MASKING
          Bfbot(i,j)=Bfbot(i,j)*rmask(i,j)
# endif
        END DO
      END DO
!
!  Under neutral and stable conditions, the depth of the bottom
!  boundary layer is required to be less than Ekman and Monin-Obukov
!  depths.
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          IF (Ustar(i,j).ge.0.0_r8) THEN
            hekman=lmd_cekman*Ustar(i,j)/MAX(ABS(f(i,j)),eps)-h(i,j)
            hbbl(i,j)= MIN(hekman,hbbl(i,j))
          END IF
          hbbl(i,j)=MIN(hbbl(i,j),z_w(i,j,N(ng)))
          hbbl(i,j)=MAX(hbbl(i,j),z_w(i,j,0))
# ifdef MASKING
          hbbl(i,j)=hbbl(i,j)*rmask(i,j)
# endif
        END DO
      END DO
!
!  Apply gradient or periodic boundary conditions.
!
      CALL bc_r2d_tile (ng, Istr, Iend, Jstr, Jend,                     &
     &                  LBi, UBi, LBj, UBj, hbbl)
!
!  Find new boundary layer index "kbbl".
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          kbbl(i,j)=N(ng)
          DO k=1,N(ng)
            IF ((kbbl(i,j).eq.N(ng)).and.(z_w(i,j,k).gt.hbbl(i,j))) THEN
              kbbl(i,j)=k
            END IF
          END DO
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Compute total buoyancy flux at final bottom boundary layer depth.
!-----------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          zgrid(i,j)=z_w(i,j,N(ng))-hbbl(i,j)
        END DO
      END DO
      CALL lmd_swfrac_tile (ng, Istr, Iend, Jstr, Jend, -1.0_r8,        &
     &                      zgrid, swdk)
      DO j=Jstr,Jend
        DO i=Istr,Iend
          Bfbot(i,j)=(Bo(i,j)+Bosol(i,j)*(1.0_r8-swdk(i,j)))
# ifdef MASKING
          Bfbot(i,j)=Bfbot(i,j)*rmask(i,j)
# endif
        END DO
      END DO
!
!=======================================================================
!  Compute vertical mixing coefficients within the planetary boundary
!  layer.
!=======================================================================
!
!  Compute tubulent velocity scales (wm,ws) at "hbbl".
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          bl_dpth(i,j)=lmd_epsilon*(hbbl(i,j)-z_w(i,j,0))
          IF (Bfbot(i,j).gt.0.0_r8) THEN
            cff=1.0_r8
          ELSE
            cff=lmd_epsilon
          END IF
          sigma=cff*(hbbl(i,j)-z_w(i,j,0))
          CALL lmd_wscale (Bfbot(i,j), Ustar(i,j), sigma,               &
     &                     wm(i,j), ws(i,j))
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Compute nondimensional shape function Gx(sigma) in terms of the
!  interior diffusivities at sigma=1 (Gm1, Gs1, Gt1) and its vertical
!  derivative evaluated "hbbl" via interpolation.
!-----------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          f1(i,j)=5.0_r8*MAX(0.0_r8,Bfbot(i,j))*vonKar/                 &
     &            (Ustar(i,j)*Ustar(i,j)*Ustar(i,j)*Ustar(i,j)+eps)
        END DO
      END DO
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          zbl=hbbl(i,j)-z_w(i,j,0)
          k=kbbl(i,j)
          cff=1.0_r8/(z_w(i,j,k)-z_w(i,j,k-1))
          cff_dn=cff*(hbbl(i,j)-z_w(i,j,k-1))
          cff_up=cff*(z_w(i,j,k)-hbbl(i,j))
!
!  Compute nondimensional shape function for viscosity "Gm1" and its
!  vertical derivative "dGm1dS" evaluated at "hbbl".
!
          K_bl=cff_dn*Kv(i,j,k)+cff_up*Kv(i,j,k-1)
          dK_bl=-cff*(Kv(i,j,k)-Kv(i,j,k-1))
          Gm1(i,j)=K_bl/(zbl*wm(i,j)+eps)
# ifdef MASKING
          Gm1(i,j)=Gm1(i,j)*rmask(i,j)
# endif
          dGm1dS(i,j)=MIN(0.0_r8, K_bl*f1(i,j)-dK_bl/(wm(i,j)+eps))
!
!  Compute nondimensional shape function for diffusion of temperature
!  "Gt1" and its vertical derivative "dGt1dS" evaluated at "hbbl".
!
          K_bl=cff_dn*Kt(i,j,k)+cff_up*Kt(i,j,k-1)
          dK_bl=-cff*(Kt(i,j,k)-Kt(i,j,k-1))
          Gt1(i,j)=K_bl/(zbl*ws(i,j)+eps)
# ifdef MASKING
          Gt1(i,j)=Gt1(i,j)*rmask(i,j)
# endif
          dGt1dS(i,j)=MIN(0.0_r8, K_bl*f1(i,j)-dK_bl/(ws(i,j)+eps))
# ifdef SALINITY
!
!  Compute nondimensional shape function for diffusion of salinity
!  "Gs1" and its vertical derivative "dGs1dS" evaluated at "hbbl".
!
          K_bl=cff_dn*Ks(i,j,k)+cff_up*Ks(i,j,k-1)
          dK_bl=-cff*(Ks(i,j,k)-Ks(i,j,k-1))
          Gs1(i,j)=K_bl/(zbl*ws(i,j)+eps)
#  ifdef MASKING
          Gs1(i,j)=Gs1(i,j)*rmask(i,j)
#  endif
          dGs1dS(i,j)=MIN(0.0_r8, K_bl*f1(i,j)-dK_bl/(ws(i,j)+eps))
# endif
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Compute bottom boundary layer mixing coefficients.
!-----------------------------------------------------------------------
!
      DO k=1,N(ng)-1
        DO j=Jstr,Jend
          DO i=Istr,Iend
            IF (z_w(i,j,k).lt.hbbl(i,j)) THEN
!
!  Compute turbulent velocity scales at vertical W-points.
!
              depth=z_w(i,j,k)-z_w(i,j,0)
              IF (Bflux(i,j,k).lt.0.0_r8) THEN
                sigma=MIN(bl_dpth(i,j),depth)
              ELSE
                sigma=depth
              END IF
              CALL lmd_wscale (Bflux(i,j,k), Ustar(i,j), sigma,         &
     &                         wm(i,j), ws(i,j))
!
!  Set polynomial coefficients for shape function.
!
              sigma=depth/(hbbl(i,j)-z_w(i,j,0)+eps)
# ifdef MASKING
              sigma=sigma*rmask(i,j)
# endif
              a1=sigma-2.0_r8
              a2=3.0_r8-2.0_r8*sigma
              a3=sigma-1.0_r8
!
!  Compute nondimesional shape functions.
!
              Gm=a1+a2*Gm1(i,j)+a3*dGm1dS(i,j)
              Gt=a1+a2*Gt1(i,j)+a3*dGt1dS(i,j)
# ifdef SALINITY
              Gs=a1+a2*Gs1(i,j)+a3*dGs1dS(i,j)
# endif
!
!  Compute boundary layer mixing coefficients, combine them
!  with interior mixing coefficients.
!
              Kv(i,j,k)=depth*wm(i,j)*(1.0_r8+sigma*Gm)
              Kt(i,j,k)=depth*ws(i,j)*(1.0_r8+sigma*Gt)
# ifdef SALINITY
              Ks(i,j,k)=depth*ws(i,j)*(1.0_r8+sigma*Gs)
# endif
            END IF
          END DO
        END DO
!
!-----------------------------------------------------------------------
!  Finally, set vertical mixing coefficients "Akv" and "Akt".
!-----------------------------------------------------------------------
!
!  Take the maximum estimate of vertical mixing only where the surface
!  and bottom boundary layers overlap. (Do not let the interior
!  overwrite the boundary layer estimate).
!
        DO j=Jstr,Jend
          DO i=Istr,Iend
            IF (z_w(i,j,k).lt.hbbl(i,j)) THEN
               if(k.gt.ksbl(i,j)) THEN
                  Akv(i,j,k)=MAX(Akv(i,j,k),Kv(i,j,k))
                  Akt(i,j,k,itemp)=MAX(Akt(i,j,k,itemp),Kt(i,j,k))
# ifdef SALINITY
                  Akt(i,j,k,isalt)=MAX(Akt(i,j,k,isalt),Ks(i,j,k))
# endif
               ELSE
                 Akv(i,j,k)=Kv(i,j,k)
                 Akt(i,j,k,itemp)=Kt(i,j,k)
# ifdef SALINITY
                 Akt(i,j,k,isalt)=Ks(i,j,k)
# endif
               END IF
            END IF
          END DO
        END DO
      END DO
!
!  Apply boundary conditions.
!
      CALL bc_w3d_tile (ng, Istr, Iend, Jstr, Jend,                     &
     &                  LBi, UBi, LBj, UBj, 0, N(ng), Akv)
      DO itrc=1,NAT
        CALL bc_w3d_tile (ng, Istr, Iend, Jstr, Jend,                   &
     &                    LBi, UBi, LBj, UBj, 0, N(ng), Akt(:,:,:,itrc))
      END DO
      RETURN
      END SUBROUTINE lmd_bkpp_tile
#endif
      END MODULE lmd_bkpp_mod
