#include "cppdefs.h"
      MODULE mod_sediment
#ifdef SEDIMENT
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!   bed      Sediment bed layer properties.                            !
!   bed_frac Sediment class fraction composition in each layer         !
!              (non-dimensional: 0.0-1.0).                             !
!   bottom   Bottom (exposed sediment layer) properties.               !
!                                                                      !
!  Parameters for sediment model:                                      !
!                                                                      !
!   Csed     Sediment concentration (mg/l), used during analytical     !
!              initialization.                                         !
!   Erate    Surface erosion rate (kg/m2/s).                           !
!   Sd50     Median sediment grain diameter (m).                       !
!   Srho     Sediment grain density (kg/m3).                           !
!   SedIter  Maximum number of iterations.                             !
!   Wsed     Particle settling velocity (m/s).                         !
!   poros    Porosity (non-dimensional: 0.0-1.0):                      !
!              Vwater/(Vwater+Vsed).                                   !
!   tau_ce   Kinematic critical shear for erosion (m2/s2).             !
!   tau_cd   Kinematic critical shear for deposition (m2/s2).          !
!                                                                      !
!  BED properties indices:                                             !
!                                                                      !
!   MBEDP    Number of bed properties (array dimension).               !
!   idSbed   IO indices for bed properties variables.                  !
!   idfrac   sediment class fraction (non-dimensional).                !
!   ithck    Sediment layer thickness (m).                             !
!   iaged    Sediment layer age (days).                                !
!   iporo    Sediment layer porosity (non-dimensional).                !
!   idiff    Bio-diffusivity at the bottom of each layer.              !
!                                                                      !
!  BOTTOM properties indices:                                          !
!                                                                      !
!   MBOTP    Number of bottom properties (array dimension).            !
!   idBott   IO indices for bottom properties variables.               !
!   irlen    Sediment ripple length (m).                               !
!   irhgt    Sediment ripple height (m).                               !
!   ishgt    Sediment saltation height (m).                            !
!   izNik    Nikuradse bottom roughness (m).                           !
!   izbio    Biological bottom roughness (m).                          !
!   izbld    Bed load bottom roughness (m).                            !
!   izbfm    Bed form bottom roughness (m).                            !
!   izapp    Apparent bottom roughness (m).                            !
!                                                                      !
!=======================================================================
!
        USE mod_param

        implicit none

        integer, parameter :: MBEDP = 4          ! Bed Properties:
        integer, parameter :: ithck = 1          ! layer thickness
        integer, parameter :: iaged = 2          ! layer age
        integer, parameter :: iporo = 3          ! layer porosity
        integer, parameter :: idiff = 4          ! bio-diffusivity

        integer, parameter :: MBOTP = 8          ! Bottom Properties:
        integer, parameter :: irlen = 1          ! ripple length
        integer, parameter :: irhgt = 2          ! ripple height
        integer, parameter :: ishgt = 3          ! saltation height
        integer, parameter :: izNik = 4          ! Nikuradse Zo
        integer, parameter :: izbio = 5          ! biological Zo
        integer, parameter :: izbld = 6          ! bed load Zo
        integer, parameter :: izbfm = 7          ! bed form Zo
        integer, parameter :: izapp = 8          ! apparent Zo

        integer  :: SedIter = 1                  ! number of iterations

        integer  :: idBott(MBOTP)                ! bottom properties IDs
        integer  :: idSbed(MBEDP)                ! bed properties IDs
        integer  :: idfrac(NST)                  ! class fraction indices

        real(r8) :: Csed(NST,Ngrids)             ! initial concentration
        real(r8) :: Erate(NST,Ngrids)            ! erosion rate
        real(r8) :: Sd50(NST,Ngrids)             ! mediam grain diameter
        real(r8) :: Srho(NST,Ngrids)             ! grain density
        real(r8) :: Wsed(NST,Ngrids)             ! settling velocity
        real(r8) :: poros(NST,Ngrids)            ! porosity
        real(r8) :: tau_ce(NST,Ngrids)           ! shear for erosion
        real(r8) :: tau_cd(NST,Ngrids)           ! shear for deposition

#endif
 
      END MODULE mod_sediment
