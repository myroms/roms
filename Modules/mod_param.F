#include "cppdefs.h"
      MODULE mod_param
!
!=======================================================================
!  Copyright (c) 2005 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  Grid parameters:                                                    !
!                                                                      !
!  Im         Number of global grid points in the XI-direction         !
!               for each nested grid.                                  !
!  Jm         Number of global grid points in the ETA-direction        !
!               for each nested grid.                                  !
!  Lm         Number of interior grid points in the XI-direction       !
!               for each nested grid.                                  !
!  Mm         Number of internal grid points in the ETA-direction.     !
!               for each nested grid.                                  !
!  N          Number of vertical levels for each nested grid.          !
!  Ngrids     Number of nested and/or connected grids to solve.        !
!  NtileI     Number of XI-direction tiles or domain partitions for    !
!               each nested grid. Values used to compute tile ranges.  !
!  NtileJ     Number of ETA-direction tiles or domain partitions for   !
!               each nested grid. Values used to compute tile ranges.  !
!  NtileX     Number of XI-direction tiles or domain partitions for    !
!               each nested grid. Values used in parallel loops.       !
!  NtileE     Number of ETA-direction tiles or domain partitions for   !
!               each nested grid. Values used in parallel loops.       !
!                                                                      !
!  Configuration parameters:                                           !
!                                                                      !
!  Nfloats    Number of floats tracjectories.                          !
!  Nstation   Number of output stations.                               !
!  Msev       Maximum number of state error eigenvectors.              !
!  MTC        Maximum number of tidal components.                      !
!  Msrc       Maximum number of point Sources/Sinks.                   !
!  NFT        Number of float time levels (minus one).                 !
!  NFV        Number of float variables.                               !
!  NSV        Number of model state variables.                         !
!                                                                      !
!  Tracer parameters:                                                  !
!                                                                      !
!  NAT        Number of active tracer type variables (usually,         !
!               NAT=2 for potential temperature and salinity).         !
!  NBT        Number of biological tracer type variables.              !
!  NST        Number of sediment tracer type variables (NCS+NNS).      !
!  NPT        Number of extra passive tracer type variables to         !
!               advect and diffuse only (dyes, etc).                   !
!  NT         Total number of tracer type variables.                   !
!  MT         Maximum number of tracer type variables.                 !
!                                                                      !
!  Nbed       Number of sediment bed layers.                           !
!  NCS        Number of cohesive (mud) sediment tracers.               !
!  NNS        Number of non-cohesive (sand) sediment tracers.          !
!                                                                      !
!  NBands     Number of spectral irradiance bands in EcoSim.           !
!  Nbac       Number of bacteria constituents in EcoSim.               !
!  Ndom       Number of dissolved matter constituents in EcoSim.       !
!  Nfec       Number of fecal matter constituents in EcoSim.           !
!  Nphy       Number of phytoplankton constituents in EcoSim.          !
!  Npig       Number of pigment constituents in EcoSim.                !
!  PHY        EcoSim indices of phytoplankton species considered.      !
!  PIG        EcoSim phytoplankton-pigment matrix.                     !
!                                                                      !
!  Diagnostic fields parameters:                                       !
!                                                                      !
!  NDbio2d    Number of diagnostic 2D biology fields.                  !
!  NDbio3d    Number of diagnostic 3D biology fields.                  !
!  NDT        Number of diagnostic tracer fields.                      !
!  NDM2d      Number of diagnostic 2D momentum fields.                 !
!  NDM3d      Number of diagnostic 3D momentum fields.                 !
!  NDrhs      Number of diagnostic 3D right-hand-side fields.          !
!                                                                      !
!=======================================================================
!
        USE mod_kinds

        implicit none

        integer, private :: ig

        integer, parameter :: Msev = 300
        integer, parameter :: Msrc = 200
#ifdef ECOSIM
        integer, parameter :: NBands = 60
#endif
#if defined BIO_FASHAM && defined DIAGNOSTICS_BIO
# ifdef CARBON
#  ifdef OXYGEN
        integer, parameter :: NDbio2d = 4
        integer, parameter :: NDbio3d = 2
#  else
        integer, parameter :: NDbio2d = 3
        integer, parameter :: NDbio3d = 2
#  endif
# else
#  ifdef OXYGEN
        integer, parameter :: NDbio2d = 2
        integer, parameter :: NDbio3d = 2
#  else
        integer, parameter :: NDbio2d = 1
        integer, parameter :: NDbio3d = 2
#  endif
# endif
#endif
!!      integer, parameter :: MTC = 7
        integer, parameter :: MTC = 20
        integer, parameter :: NFT = 4

#if defined ADRIA02
        integer, parameter :: Ngrids = 1
        integer, parameter :: NAT = 2
        integer, parameter, dimension(Ngrids) :: Lm = (/ 158 /)
        integer, parameter, dimension(Ngrids) :: Mm = (/ 58  /)
        integer, parameter, dimension(Ngrids) :: N  = (/ 20  /)
#elif defined BASIN
        integer, parameter :: Ngrids = 1
        integer, parameter :: NAT = 1
        integer, parameter, dimension(Ngrids) :: Lm = (/ 180 /)
        integer, parameter, dimension(Ngrids) :: Mm = (/ 140 /)
        integer, parameter, dimension(Ngrids) :: N  = (/ 12  /)
#elif defined BENCHMARK1
        integer, parameter :: Ngrids = 1
        integer, parameter :: NAT = 2
        integer, parameter, dimension(Ngrids) :: Lm = (/ 512 /)
        integer, parameter, dimension(Ngrids) :: Mm = (/ 64  /)
        integer, parameter, dimension(Ngrids) :: N  = (/ 30  /)
#elif defined BENCHMARK2
        integer, parameter :: Ngrids = 1
        integer, parameter :: NAT = 2
        integer, parameter, dimension(Ngrids) :: Lm = (/ 1024 /)
        integer, parameter, dimension(Ngrids) :: Mm = (/ 128  /)
        integer, parameter, dimension(Ngrids) :: N  = (/ 30   /)
#elif defined BENCHMARK3
        integer, parameter :: Ngrids = 1
        integer, parameter :: NAT = 2
        integer, parameter, dimension(Ngrids) :: Lm = (/ 2048 /)
        integer, parameter, dimension(Ngrids) :: Mm = (/ 256  /)
        integer, parameter, dimension(Ngrids) :: N  = (/ 30   /)
#elif defined BIO_TOY
        integer, parameter :: Ngrids = 1
        integer, parameter :: NAT = 2
        integer, parameter, dimension(Ngrids) :: Lm = (/ 4   /)
        integer, parameter, dimension(Ngrids) :: Mm = (/ 4   /)
        integer, parameter, dimension(Ngrids) :: N  = (/ 30  /)
#elif defined BISCAY
        integer, parameter :: Ngrids = 1
        integer, parameter :: NAT = 2
        integer, parameter, dimension(Ngrids) :: Lm = (/ 96  /)
        integer, parameter, dimension(Ngrids) :: Mm = (/ 128 /)
        integer, parameter, dimension(Ngrids) :: N  = (/ 20  /)
#elif defined BL_TEST
        integer, parameter :: Ngrids = 1
        integer, parameter :: NAT = 2
        integer, parameter, dimension(Ngrids) :: Lm = (/ 39  /)
        integer, parameter, dimension(Ngrids) :: Mm = (/ 4   /)
        integer, parameter, dimension(Ngrids) :: N  = (/ 30  /)
#elif defined CANYON_A
        integer, parameter :: Ngrids = 1
        integer, parameter :: NAT = 1
        integer, parameter, dimension(Ngrids) :: Lm = (/ 65  /)
        integer, parameter, dimension(Ngrids) :: Mm = (/ 48  /)
        integer, parameter, dimension(Ngrids) :: N  = (/ 10  /)
#elif defined CANYON_B
        integer, parameter :: Ngrids = 1
        integer, parameter :: NAT = 1
        integer, parameter, dimension(Ngrids) :: Lm = (/ 66  /)
        integer, parameter, dimension(Ngrids) :: Mm = (/ 48  /)
        integer, parameter, dimension(Ngrids) :: N  = (/ 12  /)
#elif defined CBLAST
        integer, parameter :: Ngrids = 1
        integer, parameter :: NAT = 2
        integer, parameter, dimension(Ngrids) :: Lm = (/ 160 /)
        integer, parameter, dimension(Ngrids) :: Mm = (/ 128 /)
        integer, parameter, dimension(Ngrids) :: N  = (/ 20  /)
#elif defined COUPLING_TEST
        integer, parameter :: Ngrids = 1
        integer, parameter :: NAT = 2
        integer, parameter, dimension(Ngrids) :: Lm = (/ 42  /)
        integer, parameter, dimension(Ngrids) :: Mm = (/ 42  /)
        integer, parameter, dimension(Ngrids) :: N  = (/ 15  /)
#elif defined DAMEE_4
        integer, parameter :: Ngrids = 1
        integer, parameter :: NAT = 2
        integer, parameter, dimension(Ngrids) :: Lm = (/ 128 /)
        integer, parameter, dimension(Ngrids) :: Mm = (/ 128 /)
        integer, parameter, dimension(Ngrids) :: N  = (/ 20  /)
#elif defined DOUBLE_GYRE
        integer, parameter :: Ngrids = 1
        integer, parameter ::  NAT = 1
        integer, parameter, dimension(Ngrids) :: Lm = (/ 54  /)
        integer, parameter, dimension(Ngrids) :: Mm = (/ 108 /)
        integer, parameter, dimension(Ngrids) :: N  = (/ 4   /)
!!      integer, parameter, dimension(Ngrids) :: N  = (/ 30  /)
#elif defined EAC
        integer, parameter :: Ngrids = 1
        integer, parameter ::  NAT = 2
        integer, parameter, dimension(Ngrids) :: Lm = (/ 128 /)
        integer, parameter, dimension(Ngrids) :: Mm = (/ 160 /)
        integer, parameter, dimension(Ngrids) :: N  = (/ 30  /)
#elif defined ESTUARY_TEST
        integer, parameter :: Ngrids = 1
        integer, parameter ::  NAT = 2
!!      integer, parameter, dimension(Ngrids) :: Lm = (/ 200 /)
        integer, parameter, dimension(Ngrids) :: Lm = (/ 50 /)
        integer, parameter, dimension(Ngrids) :: Mm = (/ 3   /)
        integer, parameter, dimension(Ngrids) :: N  = (/ 10  /)
#elif defined FLT_TEST
        integer, parameter :: Ngrids = 1
        integer, parameter ::  NAT = 1
        integer, parameter, dimension(Ngrids) :: Lm = (/ 14  /)
        integer, parameter, dimension(Ngrids) :: Mm = (/ 12  /)
        integer, parameter, dimension(Ngrids) :: N  = (/ 4   /)
#elif defined GRAV_ADJ
        integer, parameter :: Ngrids = 1
        integer, parameter :: NAT = 2
        integer, parameter, dimension(Ngrids) :: Lm = (/ 128 /)
        integer, parameter, dimension(Ngrids) :: Mm = (/ 4   /)
        integer, parameter, dimension(Ngrids) :: N  = (/ 40  /)
#elif defined IAS
        integer, parameter :: Ngrids = 1
        integer, parameter :: NAT = 2
        integer, parameter, dimension(Ngrids) :: Lm = (/ 256 /)
        integer, parameter, dimension(Ngrids) :: Mm = (/ 128 /)
        integer, parameter, dimension(Ngrids) :: N  = (/ 30  /)
#elif defined KELVIN
        integer, parameter :: Ngrids = 1
        integer, parameter :: NAT = 1
        integer, parameter, dimension(Ngrids) :: Lm = (/ 50  /)
        integer, parameter, dimension(Ngrids) :: Mm = (/ 30  /)
        integer, parameter, dimension(Ngrids) :: N  = (/ 10  /)
#elif defined LAB_CANYON
        integer, parameter :: Ngrids = 1
        integer, parameter :: NAT = 1
        integer, parameter, dimension(Ngrids) :: Lm = (/ 110 /)
        integer, parameter, dimension(Ngrids) :: Mm = (/ 250 /)
        integer, parameter, dimension(Ngrids) :: N  = (/ 25  /)
#elif defined LAKE_SIGNELL
        integer, parameter :: Ngrids = 1
        integer, parameter :: NAT = 2
        integer, parameter, dimension(Ngrids) :: Lm = (/ 100 /)
        integer, parameter, dimension(Ngrids) :: Mm = (/ 20  /)
        integer, parameter, dimension(Ngrids) :: N  = (/ 8   /)
#elif defined LMD_TEST
        integer, parameter :: Ngrids = 1
        integer, parameter :: NAT = 2
        integer, parameter, dimension(Ngrids) :: Lm = (/ 10  /)
        integer, parameter, dimension(Ngrids) :: Mm = (/ 10  /)
        integer, parameter, dimension(Ngrids) :: N  = (/ 50  /)
#elif defined NATL
        integer, parameter :: Ngrids = 1
        integer, parameter :: NAT = 2
        integer, parameter, dimension(Ngrids) :: Lm = (/ 1000 /)
        integer, parameter, dimension(Ngrids) :: Mm = (/ 1024 /)
        integer, parameter, dimension(Ngrids) :: N  = (/ 30   /)
#elif defined NENA
        integer, parameter :: Ngrids = 1
        integer, parameter :: NAT = 2
        integer, parameter, dimension(Ngrids) :: Lm = (/ 384 /)
        integer, parameter, dimension(Ngrids) :: Mm = (/ 128 /)
        integer, parameter, dimension(Ngrids) :: N  = (/ 30  /)
#elif defined NJ_BIGHT
        integer, parameter :: Ngrids = 1
        integer, parameter :: NAT = 2
        integer, parameter, dimension(Ngrids) :: Lm = (/ 109 /)
        integer, parameter, dimension(Ngrids) :: Mm = (/ 239 /)
        integer, parameter, dimension(Ngrids) :: N  = (/ 30  /)
!!      integer, parameter, dimension(Ngrids) :: N  = (/ 25  /)
#elif defined NPACIFIC
        integer, parameter :: Ngrids = 1
        integer, parameter :: NAT = 2
        integer, parameter, dimension(Ngrids) :: Lm = (/ 174 /)
        integer, parameter, dimension(Ngrids) :: Mm = (/ 94  /)
        integer, parameter, dimension(Ngrids) :: N  = (/ 20  /)
#elif defined OVERFLOW
        integer, parameter :: Ngrids = 1
        integer, parameter :: NAT = 1
        integer, parameter, dimension(Ngrids) :: Lm = (/ 4   /)
        integer, parameter, dimension(Ngrids) :: Mm = (/ 128 /)
        integer, parameter, dimension(Ngrids) :: N  = (/ 20  /)
#elif defined RIVERPLUME1
        integer, parameter :: Ngrids = 1
        integer, parameter :: NAT = 2
        integer, parameter, dimension(Ngrids) :: Lm = (/ 39  /)
        integer, parameter, dimension(Ngrids) :: Mm = (/ 67  /)
        integer, parameter, dimension(Ngrids) :: N  = (/ 13  /)
#elif defined RIVERPLUME2
        integer, parameter :: Ngrids = 1
        integer, parameter :: NAT = 2
        integer, parameter, dimension(Ngrids) :: Lm = (/ 50  /)
        integer, parameter, dimension(Ngrids) :: Mm = (/ 70  /)
        integer, parameter, dimension(Ngrids) :: N  = (/ 13  /)
#elif defined SCB
        integer, parameter :: Ngrids = 1
        integer, parameter :: NAT = 2
        integer, parameter, dimension(Ngrids) :: Lm = (/ 72  /)
        integer, parameter, dimension(Ngrids) :: Mm = (/ 50 /)
        integer, parameter, dimension(Ngrids) :: N  = (/ 20  /)
#elif defined SEAMOUNT
        integer, parameter :: Ngrids = 1
        integer, parameter :: NAT = 1
        integer, parameter, dimension(Ngrids) :: Lm = (/ 49  /)
        integer, parameter, dimension(Ngrids) :: Mm = (/ 48  /)
        integer, parameter, dimension(Ngrids) :: N  = (/ 13  /)
#elif defined SHELFRONT
        integer, parameter :: Ngrids = 1
        integer, parameter :: NAT = 2
        integer, parameter, dimension(Ngrids) :: Lm = (/ 4   /)
        integer, parameter, dimension(Ngrids) :: Mm = (/ 40  /)
        integer, parameter, dimension(Ngrids) :: N  = (/ 20  /)
#elif defined SOLITON
        integer, parameter :: Ngrids = 1
        integer, parameter :: NAT = 1
        integer, parameter, dimension(Ngrids) :: Lm = (/ 96  /)
        integer, parameter, dimension(Ngrids) :: Mm = (/ 32  /)
        integer, parameter, dimension(Ngrids) :: N  = (/ 3   /)
#elif defined SED_TEST1
        integer, parameter :: Ngrids = 1
        integer, parameter :: NAT = 2
        integer, parameter, dimension(Ngrids) :: Lm = (/ 100 /)
        integer, parameter, dimension(Ngrids) :: Mm = (/ 12  /)
        integer, parameter, dimension(Ngrids) :: N  = (/ 10  /)
!!      integer, parameter, dimension(Ngrids) :: Mm = (/ 3   /)
!!      integer, parameter, dimension(Ngrids) :: N  = (/ 40  /)
#elif defined SED_TOY
        integer, parameter :: Ngrids = 1
        integer, parameter :: NAT = 2
!!      integer, parameter, dimension(Ngrids) :: Lm = (/ 5 /)
        integer, parameter, dimension(Ngrids) :: Lm = (/ 100 /)
        integer, parameter, dimension(Ngrids) :: Mm = (/ 4  /)
        integer, parameter, dimension(Ngrids) :: N  = (/ 10  /)
!!      integer, parameter, dimension(Ngrids) :: N  = (/ 100  /)
#elif defined UPWELLING
        integer, parameter :: Ngrids = 1
        integer, parameter :: NAT = 2
        integer, parameter, dimension(Ngrids) :: Lm = (/ 41  /)
        integer, parameter, dimension(Ngrids) :: Mm = (/ 80  /)
        integer, parameter, dimension(Ngrids) :: N  = (/ 16  /)
#elif defined USWEST
        integer, parameter :: Ngrids = 1
        integer, parameter :: NAT = 2
        integer, parameter, dimension(Ngrids) :: Lm = (/ 126 /)
        integer, parameter, dimension(Ngrids) :: Mm = (/ 30  /)
        integer, parameter, dimension(Ngrids) :: N  = (/ 20  /)
#elif defined WEDDELL
        integer, parameter :: Ngrids = 1
        integer, parameter :: NAT = 2
        integer, parameter, dimension(Ngrids) :: Lm = (/ 300 /)
        integer, parameter, dimension(Ngrids) :: Mm = (/ 5   /)
        integer, parameter, dimension(Ngrids) :: N  = (/ 30  /)
#elif defined WINDBASIN
        integer, parameter :: Ngrids = 1
        integer, parameter :: NAT = 1
        integer, parameter, dimension(Ngrids) :: Lm = (/ 101  /)
        integer, parameter, dimension(Ngrids) :: Mm = (/ 51  /)
        integer, parameter, dimension(Ngrids) :: N  = (/ 25  /)
#else
        Illegal choice of  parameters Im, Jm, N, NAT
#endif
!
!  Set nonlinear, tangent linear, and adjoint models identifiers.
!
        integer, parameter :: iNLM = 1
        integer, parameter :: iTLM = 2
        integer, parameter :: iRPM = 3
        integer, parameter :: iADM = 4
!
!  Set number of extra passive tracers to advect and diffuse only
!  (like dyes, etc). This parameter is independent of the number of
!  biological and/or sediment tracers. It can be used 
!
#ifdef T_PASSIVE
        integer, parameter :: NPT = 1
#else
        integer, parameter :: NPT = 0
#endif

#ifdef SEDIMENT
!
!  Sediment tracers parameters.
!
# if defined ADRIA02
        integer, parameter :: Nbed = 8
        integer, parameter :: NCS = 8
        integer, parameter :: NNS = 0
        integer, parameter :: NST = NCS + NNS
# elif defined ESTUARY_TEST
        integer, parameter :: Nbed = 1
        integer, parameter :: NCS = 1
        integer, parameter :: NNS = 0
        integer, parameter :: NST = NCS + NNS
# elif defined LAKE_SIGNELL
        integer, parameter :: Nbed = 5
        integer, parameter :: NCS = 2
        integer, parameter :: NNS = 1
        integer, parameter :: NST = NCS + NNS
# elif defined SED_TEST1
        integer, parameter :: Nbed = 1
        integer, parameter :: NCS = 1
        integer, parameter :: NNS = 1
        integer, parameter :: NST = NCS + NNS
# elif defined SED_TOY
        integer, parameter :: Nbed = 1
        integer, parameter :: NCS = 1
        integer, parameter :: NNS = 0
        integer, parameter :: NST = NCS + NNS
# else
        Illegal choice of parameters Nbed, NCS, NNS, NST
# endif
#else
        integer, parameter :: Nbed = 0
        integer, parameter :: NST = 0
#endif

#ifdef BIOLOGY
!
!  Biological tracers parameters.
!
# if defined BIO_FASHAM
#  ifdef CARBON
#    ifdef OXYGEN
        integer, parameter :: NBT = 12
#    else
        integer, parameter :: NBT = 11
#    endif
#  else
#    ifdef OXYGEN
        integer, parameter :: NBT = 8
#    else
        integer, parameter :: NBT = 7
#    endif
#  endif
# elif defined ECOSIM
        integer, parameter :: Nbac = 1
        integer, parameter :: Ndom = 2
        integer, parameter :: Nfec = 2
        integer, parameter :: Nphy = 4
        integer, parameter :: Npig = 7
!
!  Determine number of EcoSim biological tracer. Currently, there is a
!  maximum of seven phytoplankton species and seven different pigments:
!
! [1] small diatom           [1] chlorophyll-a
! [2] large diatom           [2] chlorophyll-b
! [3] small dinoflagellate   [3] chlorophyll-c
! [4] large dinoflagellate   [4] photosythetic carotenoids
! [5] synechococcus          [5] photoprotective carotenoids
! [6] small prochlorococcus  [6] low  urobilin phycoeurythin carotenoids
! [7] large prochlorococcus  [7] high urobilin phycoeurythin carotenoids
!
!  The phytoplankton/pigment matrix is as follows:
!
!               P h y t o p l a n k t o n
!              [1]   [2]   [3]   [4]   [5]   [6]   [7]
!
!       t [7]   0     0     0     0     1     0     0
!       n [6]   0     0     0     0     0     0     0
!       e [5]   1     1     1     1     1     1     1
!       m [4]   1     1     1     1     0     0     0
!       g [3]   1     1     1     1     0     0     0
!       i [2]   0     0     0     0     0     1     1
!       P [1]   1     1     1     1     1     1     1
!
        integer, parameter, dimension(7,7) :: PIG = reshape (           &
     &                                      (/ 1, 1, 1, 1, 1, 1, 1,     &
     &                                         0, 0, 0, 0, 0, 1, 1,     &
     &                                         1, 1, 1, 1, 0, 0, 0,     &
     &                                         1, 1, 1, 1, 0, 0, 0,     &
     &                                         1, 1, 1, 1, 1, 1, 1,     &
     &                                         0, 0, 0, 0, 0, 0, 0,     &
     &                                         0, 0, 0, 0, 1, 0, 0 /),  &
     &                                      (/ 7, 7 /) )
!
!  Set phytoplankton species to consider (see above classification):
!
        integer, parameter, dimension(Nphy) :: PHY = (/ 1, 2, 4, 5 /)
# endif
#else
        integer, parameter :: NBT = 0
#endif
!
!  The following parameters are derived or only known during execution.
!
        common /param/ Im, Jm, MT, NT, NFV, NSV, Nfloats, Nstation,     &
#ifdef ECOSIM
     &                 NBT,                                             &
#endif
     &                 NtileI, NtileJ, NtileX, NtileE

        integer :: MT
#ifdef ECOSIM
        integer :: NBT
#endif
        integer, dimension(Ngrids) :: Im, Jm, NT, NFV, NSV
        integer, dimension(Ngrids) :: Nfloats, Nstation
        integer, dimension(Ngrids) :: NtileI, NtileJ
        integer, dimension(Ngrids) :: NtileX, NtileE
!
!  Set number of ghost-points in the halo region.  It is only used
!  in distributed-memory applications.
!
        integer :: NghostPoints = GHOST_POINTS
!
!  Diagnostic fields parameters.
!
        integer :: NDT, NDM2d, NDM3d, NDrhs

      CONTAINS

        SUBROUTINE initialize_param
!
!=======================================================================
!  Copyright (c) 2005 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine initializes several parameters in module "mod_param"   !
!  for all nested grids.                                               !
!                                                                      !
!=======================================================================
!
!  Local variable declarations
!
        integer :: I_padd, J_padd, ng
#ifdef ECOSIM
        integer :: i, j
!
!-----------------------------------------------------------------------
!  Determine number of EcoSim total bio-optical constituents:
!-----------------------------------------------------------------------
!
!       Nutrients: NO3, NO4, PO4, FeO, SiO, DIC  (6)
!        Bacteria: C, Fe, N, P                   (Nbac*4)
!             DOM: CDM, C, N, P                  (Ndom*4)
!           Fecal: C, Fe, N, P, Si               (Nfec*5)
!    Phytoplakton: C, Fe, N, P                   (Nfec*4 + Si)
!        Pigments: look table
!
        NBT=6+(Nbac*4)+(Ndom*4)+(Nfec*5)+(Nphy*4)
!
!  Add phytoplankton silica constituents.
!
        DO i=1,Nphy
          IF (PHY(i).le.2) NBT=NBT+1
        END DO
!
!  Add pigments.  Check phytoplankton-pigment table for values greater
!  than zero.
!
        DO j=1,Npig
          DO i=1,Nphy
            IF (PIG(PHY(i),j).eq.1) NBT=NBT+1
          END DO
        END DO
#endif

#ifdef DIAGNOSTICS
!
!-----------------------------------------------------------------------
!  Determine number of diagnostic variables.
!-----------------------------------------------------------------------

# ifdef DIAGNOSTICS_TS
!
!  Tracer diagnostics.
!
        NDT=4          ! Acceleration, advection, vertical diffusion
#  if defined TS_DIF2 || defined TS_DIF4
        NDT=NDT+1      ! Horizontal diffusion
#  endif
# else
        NDT=0          ! No tracer diagnostics
# endif
# ifdef DIAGNOSTICS_UV
!
!  2D Momentum diagnostics.
!
        NDM2d=3        ! Acceleration, 2D P-Gradient, stresses
#   ifdef UV_ADV
        NDM2d=NDM2d+1  ! Horizontal advection
#   endif
#   ifdef UV_COR
        NDM2d=NDM2d+1  ! Coriolis
#   endif
#   if defined UV_VIS2 || defined UV_VIS4
        NDM2d=NDM2d+1  ! Horizontal viscosity
#   endif
#  ifdef SOLVE3D
!
!  3D Momentum diagnostics and right-hand-side terms.
!
        NDM3d=3        ! Acceleration, 3D P-Gradient, vertical viscosity
        NDrhs=1        ! 3D P-Gradient
#   ifdef UV_ADV
        NDM3d=NDM3d+2  ! Horizontal and vertical advection
        NDrhs=NDrhs+2  
#   endif
#   ifdef UV_COR
        NDM3d=NDM3d+1  ! Coriolis
        NDrhs=NDrhs+1
#   endif
#   if defined UV_VIS2 || defined UV_VIS4
        NDM3d=NDM3d+1  ! Horizontal viscosity
#   endif
#   ifdef BODYFORCE
        NDrhs=NDrhs+1  ! Vertical viscosity
#   endif
#  else
        NDM3d=0        ! No 3D momentum diagnostics
        NDrhs=0
#  endif
# endif
#endif
!
!-----------------------------------------------------------------------
!  Derived dimension parameters.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          I_padd=(Lm(ng)+2)/2-(Lm(ng)+1)/2
          J_padd=(Mm(ng)+2)/2-(Mm(ng)+1)/2
          Im(ng)=Lm(ng)+I_padd
          Jm(ng)=Mm(ng)+J_padd
          NT(ng)=NAT+NBT+NST+NPT
#ifdef FLOAT_VWALK
          NFV(ng)=NT(ng)+12
#else
          NFV(ng)=NT(ng)+10
#endif
          NSV(ng)=NT(ng)+5
        END DO
!
!  Set maximum number of tracer between all nested grids.
!
        MT=MAX(2,MAXVAL(NT))

        RETURN
        END SUBROUTINE initialize_param

      END MODULE mod_param
