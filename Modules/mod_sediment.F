#include "cppdefs.h"
      MODULE mod_sediment
#if defined SEDIMENT || defined BBL_MODEL
!
!================================================== John C. Warner   ===
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  Parameters for sediment model:                                      !
!                                                                      !
!   Csed     Sediment concentration (kg/m3), used during analytical    !
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
!   iaged    Sediment layer age (s).                                   !
!   iporo    Sediment layer porosity (non-dimensional).                !
!   idiff    Sediment layer bio-diffusivity (m2/s).                    !
!                                                                      !
!  BOTTOM properties indices:                                          !
!                                                                      !
!   MBOTP    Number of bottom properties (array dimension).            !
!   idBott   IO indices for bottom properties variables.               !
!   isd50    Median sediment grain diameter (m).                       !
!   idens    Median sediment grain density (kg/m3).                    !
!   iwsed    Mean settling velocity (m/s).                             !
!   itauc    Mean critical erosion stress (m2/s2).                     !
!   irlen    Sediment ripple length (m).                               !
!   irhgt    Sediment ripple height (m).                               !
!   ibwav    Bed wave excursion amplitude (m).                         !
!   izNik    Nikuradse bottom roughness (m).                           !
!   izbio    Biological bottom roughness (m).                          !
!   izbfm    Bed form bottom roughness (m).                            !
!   izbld    Bed load bottom roughness (m).                            !
!   izapp    Apparent bottom roughness (m).                            !
!   izwbl    Bottom roughness used wave BBL (m).                       ! 
!   izdef    Default bottom roughness (m).                             !
!   iactv    Active layer thickness for erosive potential (m).         !
!   ishgt    Sediment saltation height (m).                            !
!                                                                      !
!=======================================================================
!
        USE mod_param

        implicit none

        integer, parameter :: MBEDP = 4    ! Bed Properties:
        integer, parameter :: ithck = 1    ! layer thickness
        integer, parameter :: iaged = 2    ! layer age
        integer, parameter :: iporo = 3    ! layer porosity
        integer, parameter :: idiff = 4    ! layer bio-diffusivity
 
        integer, parameter :: MBOTP = 16   ! Bottom Properties:
        integer, parameter :: isd50 = 1    ! mean grain diameter
        integer, parameter :: idens = 2    ! mean grain density
        integer, parameter :: iwsed = 3    ! mean settle velocity
        integer, parameter :: itauc = 4    ! critical erosion stress
        integer, parameter :: irlen = 5    ! ripple length
        integer, parameter :: irhgt = 6    ! ripple height
        integer, parameter :: ibwav = 7    ! wave excursion amplitude
        integer, parameter :: izdef = 8    ! default bottom roughness
        integer, parameter :: izapp = 9    ! apparent bottom roughness
        integer, parameter :: izNik = 10    ! Nikuradse bottom roughness
        integer, parameter :: izbio = 11    ! biological bottom roughness
        integer, parameter :: izbfm = 12   ! bed form bottom roughness
        integer, parameter :: izbld = 13   ! bed load bottom roughness
        integer, parameter :: izwbl = 14   ! wave bottom roughness
        integer, parameter :: iactv = 15   ! active layer thickness
        integer, parameter :: ishgt = 16   ! saltation height

        integer  :: idBmas(NST)            ! class mass indices
        integer  :: idBott(MBOTP)          ! bottom properties IDs
        integer  :: idfrac(NST)            ! class fraction indices
        integer  :: idSbed(MBEDP)          ! bed properties IDs

        real(r8) :: Csed(NST,Ngrids)       ! initial concentration
        real(r8) :: Erate(NST,Ngrids)      ! erosion rate
        real(r8) :: Sd50(NST,Ngrids)       ! mediam grain diameter
        real(r8) :: Srho(NST,Ngrids)       ! grain density
        real(r8) :: Wsed(NST,Ngrids)       ! settling velocity
        real(r8) :: poros(NST,Ngrids)      ! porosity
        real(r8) :: tau_ce(NST,Ngrids)     ! shear for erosion
        real(r8) :: tau_cd(NST,Ngrids)     ! shear for deposition

#endif
 
      END MODULE mod_sediment
