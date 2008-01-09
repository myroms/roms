#include "cppdefs.h"
      MODULE mod_sediment
#if defined SEDIMENT || defined BBL_MODEL
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group        John C. Warner   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
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
!   morph_fac  Morphological scale factor (nondimensional).             !
!   newlayer_thick    New layer deposit thickness criteria (m).        !
!   bedload_coeff     Bedload rate coefficient (nondimensional).       !
!                                                                      !
!  BED properties indices:                                             !
!                                                                      !
!   MBEDP    Number of bed properties (array dimension).               !
!   idBmas   Sedimen mass index.                                       !
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
!   izdef    Default bottom roughness (m).                             !
!   izapp    Apparent bottom roughness (m).                            !
!   izNik    Nikuradse bottom roughness (m).                           !
!   izbio    Biological bottom roughness (m).                          !
!   izbfm    Bed form bottom roughness (m).                            !
!   izbld    Bed load bottom roughness (m).                            !
!   izwbl    Bottom roughness used wave BBL (m).                       ! 
!   iactv    Active layer thickness for erosive potential (m).         !
!   ishgt    Sediment saltation height (m).                            !
!   idoff    Offset for calculation of dmix erodibility profile (m)    !
!   idslp    Slope for calculation of dmix or erodibility profile      !
!   idtim    Time scale for restoring erodibility profile (s)          !
!   idbmx    bed biodifusivity max                                     !
!   idbmm    bed biodifusivity m                                       !
!   idbzs    bed biodifusivity zs                                      !
!   idbzm    bed biodifusivity zm                                      !
!   idbzp    bed biodifusivity phi                                     !
!                                                                      !
!=======================================================================
!
        USE mod_param

        implicit none

# ifdef SED_BIODIFF
        integer, parameter :: MBEDP = 5    ! Bed Properties:
# else
        integer, parameter :: MBEDP = 4    ! Bed Properties:
# endif
        integer, parameter :: ithck = 1    ! layer thickness
        integer, parameter :: iaged = 2    ! layer age
        integer, parameter :: iporo = 3    ! layer porosity
        integer, parameter :: idiff = 4    ! layer bio-diffusivity
# ifdef SED_BIODIFF
        integer, parameter :: ibtcr = 5    ! layer critical stress
# endif 

# ifdef SED_BIODIFF
        integer, parameter :: MBOTP = 24   ! Bottom Properties:
# else
        integer, parameter :: MBOTP = 16   ! Bottom Properties:
# endif
        integer, parameter :: isd50 = 1    ! mean grain diameter
        integer, parameter :: idens = 2    ! mean grain density
        integer, parameter :: iwsed = 3    ! mean settle velocity
        integer, parameter :: itauc = 4    ! critical erosion stress
        integer, parameter :: irlen = 5    ! ripple length
        integer, parameter :: irhgt = 6    ! ripple height
        integer, parameter :: ibwav = 7    ! wave excursion amplitude
        integer, parameter :: izdef = 8    ! default bottom roughness
        integer, parameter :: izapp = 9    ! apparent bottom roughness
        integer, parameter :: izNik = 10   ! Nikuradse bottom roughness
        integer, parameter :: izbio = 11   ! biological bottom roughness
        integer, parameter :: izbfm = 12   ! bed form bottom roughness
        integer, parameter :: izbld = 13   ! bed load bottom roughness
        integer, parameter :: izwbl = 14   ! wave bottom roughness
        integer, parameter :: iactv = 15   ! active layer thickness
        integer, parameter :: ishgt = 16   ! saltation height
#ifdef SED_BIODIFF
        integer, parameter :: idoff = 17   ! mixing depth offset
        integer, parameter :: idslp = 18   ! mixing depth slope
        integer, parameter :: idtim = 19   ! erodibility time scale
        integer, parameter :: idbmx = 20   ! diffusivity db_max
        integer, parameter :: idbmm = 21   ! diffusivity db_m
        integer, parameter :: idbzs = 22   ! diffusivity db_zs
        integer, parameter :: idbzm = 23   ! diffusivity db_zm
        integer, parameter :: idbzp = 24   ! diffusivity db_zphi
#endif

        integer  :: idBott(MBOTP)          ! bottom properties IDs
        integer  :: idSbed(MBEDP)          ! bed properties IDs

        integer, allocatable :: idBmas(:)  ! class mass indices
        integer, allocatable :: idfrac(:)  ! class fraction indices
        integer, allocatable :: idUbld(:)  ! bed load u-points
        integer, allocatable :: idVbld(:)  ! bed load v-points

        real(r8) :: newlayer_thick(Ngrids) ! deposit thickness criteria
        real(r8) :: bedload_coeff(Ngrids)  ! bedload rate coefficient

        real(r8), allocatable :: Csed(:,:)       ! initial concentration
        real(r8), allocatable :: Erate(:,:)      ! erosion rate
        real(r8), allocatable :: Sd50(:,:)       ! mediam grain diameter
        real(r8), allocatable :: Srho(:,:)       ! grain density
        real(r8), allocatable :: Wsed(:,:)       ! settling velocity
        real(r8), allocatable :: poros(:,:)      ! porosity
        real(r8), allocatable :: tau_ce(:,:)     ! shear for erosion
        real(r8), allocatable :: tau_cd(:,:)     ! shear for deposition
        real(r8), allocatable :: morph_fac(:,:)  ! morphological factor
#endif
 
      END MODULE mod_sediment
