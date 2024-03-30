      MODULE mod_sediment
!
!git $Id$
!svn $Id: sediment_mod.h 1210 2024-01-03 22:03:03Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2024 The ROMS/TOMS Group        John C. Warner   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  Parameters for sediment model:                                      !
!  =============================                                       !
!                                                                      !
!   Csed            Sediment concentration (kg/m3), used during        !
!                     analytical initialization.                       !
!   Erate           Surface erosion rate (kg/m2/s).                    !
!   Sd50            Median sediment grain diameter (m).                !
!   Srho            Sediment grain density (kg/m3).                    !
!   SedIter         Maximum number of iterations.                      !
!   Wsed            Particle settling velocity (m/s).                  !
!   poros           Porosity (non-dimensional: 0.0-1.0):               !
!                     Vwater/(Vwater+Vsed).                            !
!   tau_ce          Kinematic critical shear for erosion (m2/s2).      !
!   tau_cd          Kinematic critical shear for deposition (m2/s2).   !
!                                                                      !
!   bedload_coeff   Bedload rate coefficient (nondimensional).         !
!   minlayer_thick  Minimum thickness for 2nd layer (m).               !
!   newlayer_thick  New layer deposit thickness criteria (m).          !
!   morph_fac       Morphological scale factor (nondimensional).       !
!                                                                      !
!  BED properties indices:                                             !
!  ======================                                              !
!                                                                      !
!   MBEDP           Number of bed properties (array dimension).        !
!   idBmas(:)       Sediment mass index.                               !
!   idSbed(:)       IO indices for bed properties variables.           !
!   idfrac(:)       Sediment class fraction (non-dimensional).         !
!   ithck           Sediment layer thickness (m).                      !
!   iaged           Sediment layer age (s).                            !
!   iporo           Sediment layer porosity (non-dimensional).         !
!   idiff           Sediment layer bio-diffusivity (m2/s).             !
!   ibtcr           Sediment critical stress for erosion (Pa).         !
!                                                                      !
!  BOTTOM properties indices:                                          !
!  =========================                                           !
!                                                                      !
!   MBOTP           Number of bottom properties (array dimension).     !
!   idBott(:)       IO indices for bottom properties variables.        !
!   isd50           Median sediment grain diameter (m).                !
!   idens           Median sediment grain density (kg/m3).             !
!   iwsed           Mean settling velocity (m/s).                      !
!   itauc           Mean critical erosion stress (m2/s2).              !
!   irlen           Sediment ripple length (m).                        !
!   irhgt           Sediment ripple height (m).                        !
!   ibwav           Bed wave excursion amplitude (m).                  !
!   izdef           Default bottom roughness (m).                      !
!   izapp           Apparent bottom roughness (m).                     !
!   izNik           Nikuradse bottom roughness (m).                    !
!   izbio           Biological bottom roughness (m).                   !
!   izbfm           Bed form bottom roughness (m).                     !
!   izbld           Bed load bottom roughness (m).                     !
!   izwbl           Bottom roughness used wave BBL (m).                !
!   iactv           Active layer thickness for erosive potential (m).  !
!   ishgt           Sediment saltation height (m).                     !
!   idefx           Erosion flux.                                      !
!   idnet           Erosion or deposition.                             !
!   idoff           Offset for calculation of dmix erodibility         !
!                     profile (m).                                     !
!   idslp           Slope  for calculation of dmix or erodibility      !
!                     profile.                                         !
!   idtim           Time scale for restoring erodibility profile (s).  !
!   idbmx           Bed biodifusivity maximum.                         !
!   idbmm           Bed biodifusivity minimum.                         !
!   idbzs           Bed biodifusivity zs.                              !
!   idbzm           Bed biodifusivity zm.                              !
!   idbzp           Bed biodifusivity phi.                             !
!   idprp           Cohesive behavior.                                 !
!                                                                      !
!=======================================================================
!
      USE mod_param
!
      implicit none
!
!-----------------------------------------------------------------------
!  Tracer identification indices.
!-----------------------------------------------------------------------
!
      integer, allocatable :: idsed(:)    ! Cohesive and non-cohesive
      integer, allocatable :: idmud(:)    ! Cohesive sediment
      integer, allocatable :: isand(:)    ! Non-cohesive sediment
!
!-----------------------------------------------------------------------
!  Bed and bottom properties indices.
!-----------------------------------------------------------------------
!
!  Set size of properties arrays.
!
#if defined COHESIVE_BED || defined SED_BIODIFF || defined MIXED_BED
      integer, parameter :: MBEDP = 5      ! Bed properties
#else
      integer, parameter :: MBEDP = 4      ! Bed properties
#endif
      integer  :: idSbed(MBEDP)            ! bed properties IDs
!
#if defined MIXED_BED
      integer, parameter :: MBOTP = 27     ! Bottom properties
#elif defined COHESIVE_BED || defined SED_BIODIFF
      integer, parameter :: MBOTP = 26     ! Bottom properties
#else
      integer, parameter :: MBOTP = 18     ! Bottom properties
#endif
      integer  :: idBott(MBOTP)            ! bottom properties IDs
!
!  Set properties indices.
!
      integer, parameter :: ithck = 1      ! layer thickness
      integer, parameter :: iaged = 2      ! layer age
      integer, parameter :: iporo = 3      ! layer porosity
      integer, parameter :: idiff = 4      ! layer bio-diffusivity
#if defined COHESIVE_BED || defined SED_BIODIFF || defined MIXED_BED
      integer, parameter :: ibtcr = 5      ! layer critical stress
#endif
!
      integer, parameter :: isd50 = 1      ! mean grain diameter
      integer, parameter :: idens = 2      ! mean grain density
      integer, parameter :: iwsed = 3      ! mean settle velocity
      integer, parameter :: itauc = 4      ! critical erosion stress
      integer, parameter :: irlen = 5      ! ripple length
      integer, parameter :: irhgt = 6      ! ripple height
      integer, parameter :: ibwav = 7      ! wave excursion amplitude
      integer, parameter :: izdef = 8      ! default bottom roughness
      integer, parameter :: izapp = 9      ! apparent bottom roughness
      integer, parameter :: izNik = 10     ! Nikuradse bottom roughness
      integer, parameter :: izbio = 11     ! biological bottom roughness
      integer, parameter :: izbfm = 12     ! bed form bottom roughness
      integer, parameter :: izbld = 13     ! bed load bottom roughness
      integer, parameter :: izwbl = 14     ! wave bottom roughness
      integer, parameter :: iactv = 15     ! active layer thickness
      integer, parameter :: ishgt = 16     ! saltation height
      integer, parameter :: idefx = 17     ! erosion flux
      integer, parameter :: idnet = 18     ! erosion or deposition
#if defined COHESIVE_BED || defined SED_BIODIFF || defined MIXED_BED
      integer, parameter :: idoff = 19     ! tau critical offset
      integer, parameter :: idslp = 20     ! tau critical slope
      integer, parameter :: idtim = 21     ! erodibility time scale
      integer, parameter :: idbmx = 22     ! diffusivity db_max
      integer, parameter :: idbmm = 23     ! diffusivity db_m
      integer, parameter :: idbzs = 24     ! diffusivity db_zs
      integer, parameter :: idbzm = 25     ! diffusivity db_zm
      integer, parameter :: idbzp = 26     ! diffusivity db_zphi
#endif
#if defined MIXED_BED
      integer, parameter :: idprp = 27     ! cohesive behavior
#endif
!
!  Sediment metadata indices vectors.
!
      integer, allocatable :: idBmas(:)    ! class mass indices
      integer, allocatable :: idfrac(:)    ! class fraction indices
      integer, allocatable :: idUbld(:)    ! bed load u-points
      integer, allocatable :: idVbld(:)    ! bed load v-points
!
!-----------------------------------------------------------------------
!  Input sediment parameters.
!-----------------------------------------------------------------------
!
      real(r8), allocatable :: newlayer_thick(:)   ! deposit thickness
      real(r8), allocatable :: minlayer_thick(:)   ! 2nd layer thickness
      real(r8), allocatable :: bedload_coeff(:)    ! bedload rate

      real(r8), allocatable :: Csed(:,:)       ! initial concentration
      real(r8), allocatable :: Erate(:,:)      ! erosion rate
      real(r8), allocatable :: Sd50(:,:)       ! mediam grain diameter
      real(r8), allocatable :: Srho(:,:)       ! grain density
      real(r8), allocatable :: Wsed(:,:)       ! settling velocity
      real(r8), allocatable :: poros(:,:)      ! porosity
      real(r8), allocatable :: tau_ce(:,:)     ! shear for erosion
      real(r8), allocatable :: tau_cd(:,:)     ! shear for deposition
      real(r8), allocatable :: morph_fac(:,:)  ! morphological factor

#if defined COHESIVE_BED || defined MIXED_BED
      real(r8), allocatable :: tcr_min(:)      ! minimum shear erosion
      real(r8), allocatable :: tcr_max(:)      ! maximum shear erosion
      real(r8), allocatable :: tcr_slp(:)      ! Tau_crit profile slope
      real(r8), allocatable :: tcr_off(:)      ! Tau_crit profile offset
      real(r8), allocatable :: tcr_tim(:)      ! Tau_crit consolidation
#endif

#if defined MIXED_BED
      real(r8), allocatable :: transC(:)       ! cohesive transition
      real(r8), allocatable :: transN(:)       ! noncohesive transition
#endif
!
      CONTAINS
!
      SUBROUTINE initialize_sediment
!
!=======================================================================
!                                                                      !
!  This routine sets several variables needed by the sediment model.   !
!  It allocates and assigns sediment tracers indices.                  !
!                                                                      !
!=======================================================================
!
!  Local variable declarations
!
      integer :: i, ic

      real(r8), parameter :: IniVal = 0.0_r8
!
!-----------------------------------------------------------------------
!  Initialize tracer identification indices.
!-----------------------------------------------------------------------
!
!  Allocate various nested grid depended parameters
!
      IF (.not.allocated(newlayer_thick)) THEN
        allocate ( newlayer_thick(Ngrids) )
        newlayer_thick = IniVal
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(minlayer_thick)) THEN
        allocate ( minlayer_thick(Ngrids) )
        minlayer_thick = IniVal
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(bedload_coeff)) THEN
        allocate ( bedload_coeff(Ngrids) )
        bedload_coeff = IniVal
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

#if defined COHESIVE_BED || defined MIXED_BED
      IF (.not.allocated(tcr_min)) THEN
        allocate ( tcr_min(Ngrids) )
        tcr_min = IniVal
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(tcr_max)) THEN
        allocate ( tcr_max(Ngrids) )
        tcr_max = IniVal
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(tcr_slp)) THEN
        allocate ( tcr_slp(Ngrids) )
        tcr_slp = IniVal
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(tcr_off)) THEN
        allocate ( tcr_off(Ngrids) )
        tcr_off = IniVal
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(tcr_tim)) THEN
        allocate ( tcr_tim(Ngrids) )
        tcr_tim = IniVal
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
#endif

#if defined MIXED_BED
      IF (.not.allocated(transC)) THEN
        allocate ( transC(Ngrids) )
        transC = IniVal
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(transN)) THEN
        allocate ( transN(Ngrids) )
        transN = IniVal
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
#endif
!
!  Allocate sediment tracers indices vectors.
!
      IF (.not.allocated(idsed)) THEN
        allocate ( idsed(MAX(1,NST)) )
        Dmem(1)=Dmem(1)+REAL(MAX(1,NST),r8)
      END IF

      IF (.not.allocated(idmud)) THEN
        allocate ( idmud(MAX(1,NCS)) )
        Dmem(1)=Dmem(1)+REAL(MAX(1,NCS),r8)
      END IF

      IF (.not.allocated(isand)) THEN
        allocate ( isand(MAX(1,NNS)) )
        Dmem(1)=Dmem(1)+REAL(MAX(1,NNS),r8)
      END IF

      IF (.not.allocated(idBmas)) THEN
        allocate ( idBmas(NST) )
        Dmem(1)=Dmem(1)+REAL(NST,r8)
      END IF

      IF (.not.allocated(idfrac)) THEN
        allocate ( idfrac(NST) )
        Dmem(1)=Dmem(1)+REAL(NST,r8)
      END IF

      IF (.not.allocated(idUbld)) THEN
        allocate ( idUbld(NST) )
        Dmem(1)=Dmem(1)+REAL(NST,r8)
      END IF

      IF (.not.allocated(idVbld)) THEN
        allocate ( idVbld(NST) )
        Dmem(1)=Dmem(1)+REAL(NST,r8)
      END IF
!
!  Set cohesive and noncohesive suspended sediment tracers
!  identification indices.
!
      ic=NAT+NPT
      DO i=1,NCS
        ic=ic+1
        idmud(i)=ic
        idsed(i)=idmud(i)
      END DO
      DO i=1,NNS
        ic=ic+1
        isand(i)=ic
        idsed(NCS+i)=isand(i)
      END DO
!
      RETURN
      END SUBROUTINE initialize_sediment

      END MODULE mod_sediment
