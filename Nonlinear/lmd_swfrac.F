#include "cppdefs.h"
#if (defined LMD_SKPP || defined SOLAR_SOURCE) && defined SOLVE3D
      SUBROUTINE lmd_swfrac_tile (ng, Istr, Iend, Jstr, Jend, Zscale,   &
     &                            Z, swdk)
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This subroutine computes the  fraction  of  solar shortwave flux    !
!  penetrating to specified depth (times Zscale) due to exponential    !
!  decay in Jerlov water type.                                         !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     Zscale   scale factor to apply to depth array.                   !
!     Z        vertical height (meters, negative) for                  !
!              desired solar short-wave fraction.                      !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     swdk     shortwave (radiation) fractional decay.                 !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!  Paulson, C.A., and J.J. Simpson, 1977: Irradiance meassurements     !
!     in the upper ocean, J. Phys. Oceanogr., 7, 952-956.              !
!                                                                      !
!  This routine was adapted from Bill Large 1995 code.                 !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_mixing
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr

      real(r8), intent(in) :: Zscale

      real(r8), intent(in) :: Z(PRIVATE_2D_SCRATCH_ARRAY)
      real(r8), intent(out) :: swdk(PRIVATE_2D_SCRATCH_ARRAY)
!
!  Local variable declarations.
!
      integer :: Jindex, i, j
!
!-----------------------------------------------------------------------
!  Use Paulson and Simpson (1977) two wavelength bands solar
!  absorption model.
!-----------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          Jindex=MIXING(ng)%Jwtype(i,j)
          swdk(i,j)=lmd_r1(Jindex)*                                     &
     &              EXP(Z(i,j)*Zscale/lmd_mu1(Jindex))+                 &
     &              (1.0_r8-lmd_r1(Jindex))*                            &
     &              EXP(Z(i,j)*Zscale/lmd_mu2(Jindex))
        END DO
      END DO
      END SUBROUTINE lmd_swfrac_tile
#else
      SUBROUTINE lmd_swfrac
      RETURN
      END SUBROUTINE lmd_swfrac
#endif

