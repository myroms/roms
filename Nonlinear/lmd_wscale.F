#include "cppdefs.h"
#if defined NONLINEAR && (defined LMD_SKPP || defined LMD_BKPP)
      SUBROUTINE lmd_wscale (Bfsfc, Ustar, sigma, wm, ws)
!
!=======================================================================
!  Copyright (c) 2005 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine computes the turbulent velocity scale for momentum     !
!  and tracer as a function of turbulent friction velocity "Ustar"     !
!  and stability parameter "zetahat".                                  !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     Bfsfc    surface buoyancy forcing (m2/s2).                       !
!     Ustar    Tubulent friction velocity (m/s).                       !
!     sigma    boundary layer depth (m).                               !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     wm       turbulent velocity scale (m/s) at sigma for momentum.   !
!     ws       turbulent velocity scale (m/s) at sigma for tracer.     !
!                                                                      !
!  This routine was adapted from Bill Large 1995 code.                 !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_scalars

      implicit none
!
!  Imported variable declarations.
!
      real(r8), intent(in) :: Bfsfc, Ustar, sigma

      real(r8), intent(out) :: wm, ws
!
!  Local variable declarations.
!
      real(r8), parameter :: eps = 1.0E-20_r8
      real(r8), parameter :: r3 = 1.0_r8/3.0_r8

      real(r8) :: Ustar3, zetahat, zetapar
!
!-----------------------------------------------------------------------
!  Compute turbulent velocity scale for momentum and tracers.
!-----------------------------------------------------------------------
!
      Ustar3=Ustar*Ustar*Ustar
      zetahat=vonKar*sigma*Bfsfc
      zetapar=zetahat/(Ustar3+eps)
!
!  Stable regime.
!
      IF (zetahat.ge.0.0_r8) THEN
        wm=vonKar*Ustar/(1.0_r8+5.0_r8*zetapar)
        ws=wm
!
!  Unstable regime.
!
      ELSE
        IF (zetapar.gt.lmd_zetam) THEN
          wm=vonKar*Ustar*(1.0_r8-16.0_r8*zetapar)**0.25_r8
        ELSE
          wm=vonKar*(lmd_am*Ustar3-lmd_cm*zetahat)**r3
        END IF
        IF (zetapar.gt.lmd_zetas) THEN
          ws=vonKar*Ustar*(1.0_r8-16.0_r8*zetapar)**0.5_r8
        ELSE
          ws=vonKar*(lmd_as*Ustar3-lmd_cs*zetahat)**r3
        END IF
      END IF
#else
      SUBROUTINE lmd_wscale
#endif
      RETURN
      END SUBROUTINE lmd_wscale
