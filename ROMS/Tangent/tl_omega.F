#include "cppdefs.h"
      MODULE tl_omega_mod
#if defined TANGENT && defined SOLVE3D
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2007 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!======================================================================= 
!                                                                      !
!  This routine computes S-coordinate vertical velocity (m^3/s),       !
!                                                                      !
!                  W=[Hz/(m*n)]*omega,                                 !
!                                                                      !
!  diagnostically at horizontal RHO-points and vertical W-points.      !
!                                                                      !
!  BASIC STATE variables needed: W, z_w.                               !
!                                                                      !
!  NOTE: We need to recompute  basic state W in this routine since     !
!  ----  intermediate values of W are needed by the tangent linear     !
!        and adjoint routines.                                         !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: tl_omega

      CONTAINS
!
!***********************************************************************
      SUBROUTINE tl_omega (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_ocean
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
      CALL wclock_on (ng, iTLM, 13)
# endif
      CALL tl_omega_tile (ng, Istr, Iend, Jstr, Jend,                   &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    GRID(ng) % Huon,                              &
     &                    GRID(ng) % Hvom,                              &
     &                    GRID(ng) % z_w,                               &
     &                    GRID(ng) % tl_Huon,                           &
     &                    GRID(ng) % tl_Hvom,                           &
     &                    GRID(ng) % tl_z_w,                            &
     &                    OCEAN(ng) % W,                                &
     &                    OCEAN(ng) % tl_W)
# ifdef PROFILE
      CALL wclock_off (ng, iTLM, 13)
# endif
      RETURN
      END SUBROUTINE tl_omega
!
!***********************************************************************
      SUBROUTINE tl_omega_tile (ng, Istr, Iend, Jstr, Jend,             &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Huon, Hvom, z_w,                        &
     &                          tl_Huon, tl_Hvom, tl_z_w,               &
     &                          W, tl_W)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE bc_3d_mod, ONLY : bc_w3d_tile
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange3d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: Huon(LBi:,LBj:,:)
      real(r8), intent(in) :: Hvom(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: tl_Huon(LBi:,LBj:,:)
      real(r8), intent(in) :: tl_Hvom(LBi:,LBj:,:)
      real(r8), intent(in) :: tl_z_w(LBi:,LBj:,0:)

      real(r8), intent(out) :: W(LBi:,LBj:,0:)
      real(r8), intent(out) :: tl_W(LBi:,LBj:,0:)
# else
      real(r8), intent(in) :: Huon(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: Hvom(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(in) :: tl_Huon(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: tl_Hvom(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: tl_z_w(LBi:UBi,LBj:UBj,0:N(ng))

      real(r8), intent(out) :: W(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(out) :: tl_W(LBi:UBi,LBj:UBj,0:N(ng))
# endif
!
!  Local variable declarations.
!
# ifdef DISTRIBUTE
#  ifdef EW_PERIODIC
      logical :: EWperiodic=.TRUE.
#  else
      logical :: EWperiodic=.FALSE.
#  endif
#  ifdef NS_PERIODIC
      logical :: NSperiodic=.TRUE.
#  else
      logical :: NSperiodic=.FALSE.
#  endif
# endif
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j, k

      real(r8) :: cff, tl_cff

      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY) :: wrk
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY) :: tl_wrk

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Vertically integrate horizontal mass flux divergence.
!-----------------------------------------------------------------------
!
!  Starting with zero vertical velocity at the bottom, integrate
!  from the bottom (k=0) to the free-surface (k=N).  The w(:,:,N)
!  contains the vertical velocity at the free-surface, d(zeta)/d(t).
!  Notice that barotropic mass flux divergence is not used directly.
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          W(i,j,0)=0.0_r8
          tl_W(i,j,0)=0.0_r8
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
     &               (Huon(i+1,j,k)-Huon(i,j,k)+                        &
     &                Hvom(i,j+1,k)-Hvom(i,j,k))
            tl_W(i,j,k)=tl_W(i,j,k-1)-                                  &
     &                  (tl_Huon(i+1,j,k)-tl_Huon(i,j,k)+               &
     &                   tl_Hvom(i,j+1,k)-tl_Hvom(i,j,k))
          END DO
        END DO
        DO i=Istr,Iend
          cff=1.0_r8/(z_w(i,j,N(ng))-z_w(i,j,0))
          tl_cff=-cff*cff*(tl_z_w(i,j,N(ng))-tl_z_w(i,j,0))
          wrk(i)=cff*W(i,j,N(ng))
          tl_wrk(i)=tl_cff*W(i,j,N(ng))+cff*tl_W(i,j,N(ng))
        END DO
!
!  In order to insure zero vertical velocity at the free-surface,
!  subtract the vertical velocities of the moving S-coordinates
!  isosurfaces. These isosurfaces are proportional to d(zeta)/d(t).
!  The proportionaly coefficients are a linear function of the
!  S-coordinate with zero value at the bottom (k=0) and unity at
!  the free-surface (k=N).
!
        DO k=N(ng)-1,1,-1
          DO i=Istr,Iend
            W(i,j,k)=W(i,j,k)-                                          &
     &               wrk(i)*(z_w(i,j,k)-z_w(i,j,0))
            tl_W(i,j,k)=tl_W(i,j,k)-                                    &
     &                  tl_wrk(i)*(z_w(i,j,k)-z_w(i,j,0))-              &
     &	                wrk(i)*(tl_z_w(i,j,k)-tl_z_w(i,j,0))
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
      CALL bc_w3d_tile (ng, Istr, Iend, Jstr, Jend,                     &
     &                  LBi, UBi, LBj, UBj, 0, N(ng),                   &
     &                  W)
      CALL bc_w3d_tile (ng, Istr, Iend, Jstr, Jend,                     &
     &                  LBi, UBi, LBj, UBj, 0, N(ng),                   &
     &                  tl_W)
# ifdef DISTRIBUTE
      CALL mp_exchange3d (ng, iTLM, 1, Istr, Iend, Jstr, Jend,          &
     &                    LBi, UBi, LBj, UBj, 0, N(ng),                 &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    W)
      CALL mp_exchange3d (ng, iTLM, 1, Istr, Iend, Jstr, Jend,          &
     &                    LBi, UBi, LBj, UBj, 0, N(ng),                 &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    tl_W)
# endif

      RETURN
      END SUBROUTINE tl_omega_tile
#endif
      END MODULE tl_omega_mod
