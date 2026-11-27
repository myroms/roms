!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2024 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  Parameters for IOP-based, NPZD ecosystem model:                     !
!                                                                      !
!  AttSW          Light attenuation due to sea water (1/m).            !
!  BioIter        Maximum number of iterations to achieve convergence  !
!                   of the nonlinear solution.                         !
!  BioIni         Initial concentration for analytical initial         !
!                   (uniform) conditions.                              !
!  AttSW          Light attenuation due to seawater (1/m).             !
!  BphyMap        Mapping from phytoplankton backscatter to Nitrogen   !
!                   biomass at 440 nm (millimole_N/m2).                !
!  BdetMap        Mapping from detritus backscatter to Nitrogen        !
!                   biomass at 440 nm (millimole_N/m2).                !
!  CDOM_LightAtt  Light attenuation factor due to CDOM absorption      !
!                   (nondimensional).                                  !
!  CDOM_sigma     Light-dependent degradation rate for CDOM absorption !
!                   (1/(Watt m-2 day)).                                !
!  DetRR          Detritus remineraliztion rate (1/day).               !
!  K_DIN          Half-saturation for phytoplankton nitrogen uptake    !
!                   (millimole_N m-3).                                 !
!  PARfrac        Fraction of shortwave radiation that is              !
!                   photosynthetically active (nondimensional).        !
!  PhotoRmax      Maximum photosynthetic rate (1/day).                 !
!  PhyIS          Phytoplankton, initial slope of the P-I curve        !
!                   (m2/Watt).                                         !
!  PhyMRD         Phytoplankton mortality rate to the Detritus pool    !
!                   (1/day).                                           !
!  PhyMRN         Phytoplankton loss rate to the Nitrogen pool         !
!                   (1/day).                                           !
!  ThetaM         Maximum ratio of phytoplankton backscatter to        !
!                   absorption (nondimensional).                       !
!  Vm_DIN         Nitrogen uptake rate (1/day).                        !
!  wDet           Detrital sinking rate (m/day).                       !
!  wPhy           Phytoplankton sinking rate (m/day).                  !
!                                                                      !
!=======================================================================
!
      USE mod_param
!
      implicit none
!
!  Set number of spectral bands to consider.
!
      integer, parameter :: NBands = 2
!
!  Set biological tracer identification indices.
!
      integer, allocatable :: idbio(:)   ! Biological tracers
      integer :: iAphy(NBands)           ! Phytoplankton absorption
      integer :: iBphy(NBands)           ! Phytoplankton backscattering
      integer :: iBdet(NBands)           ! Detritus backscattering
      integer :: aCDOM(NBands)           ! CDOM absorption
      integer :: iDIN_                   ! Dissolved Inorganic Nitrogen
      integer :: i440n                   ! 440nm spectral band
      integer :: i510n                   ! 510nm spectral band
!
!  Biological parameters.
!
      integer, allocatable :: BioIter(:)

      real(r8), allocatable :: BioIni(:,:)

      real(r8), allocatable :: AttSW(:)             ! 1/m
      real(r8), allocatable :: BphyMap(:)           ! mmole_N/m2
      real(r8), allocatable :: BdetMap(:)           ! mmole_N/m2
      real(r8), allocatable :: CDOM_LightAtt(:)     ! nondimensional
      real(r8), allocatable :: CDOM_sigma(:)        ! 1/(W m-2 day)
      real(r8), allocatable :: DetRR(:)             ! 1/day
      real(r8), allocatable :: K_DIN(:)             ! millimol_N m-3
      real(r8), allocatable :: PARfrac(:)           ! nondimensional
      real(r8), allocatable :: PhotoRmax(:)         ! 1/day
      real(r8), allocatable :: PhyIS(:)             ! m2/W
      real(r8), allocatable :: PhyMRD(:)            ! 1/day
      real(r8), allocatable :: PhyMRN(:)            ! 1/day
      real(r8), allocatable :: ThetaM(:)            ! nondimensional
      real(r8), allocatable :: Vm_DIN(:)            ! 1/day
      real(r8), allocatable :: wDet(:)              ! m/day
      real(r8), allocatable :: wPhy(:)              ! m/day

#ifdef TANGENT
      real(r8), allocatable :: tl_PARfrac(:)        ! nondimensional
#endif
#ifdef ADJOINT
      real(r8), allocatable :: ad_PARfrac(:)        ! nondimensional
#endif

      CONTAINS

      SUBROUTINE initialize_biology
!
!=======================================================================
!                                                                      !
!  This routine sets several variables needed by the biology model.    !
!  It allocates and assigns biological tracers indices.                !
!                                                                      !
!=======================================================================
!
!  Local variable declarations
!
      integer :: i, ic
!
!-----------------------------------------------------------------------
!  Set number of biological tracers.
!-----------------------------------------------------------------------
!
      NBT=9
!
!-----------------------------------------------------------------------
!  Allocate various module variables.
!-----------------------------------------------------------------------
!
      IF (.not.allocated(BioIter)) THEN
        allocate ( BioIter(Ngrids) )
      END IF
      IF (.not.allocated(AttSW)) THEN
        allocate ( AttSW(Ngrids) )
      END IF
      IF (.not.allocated(BphyMap)) THEN
        allocate ( BphyMap(Ngrids) )
      END IF
      IF (.not.allocated(BdetMap)) THEN
        allocate ( BdetMap(Ngrids) )
      END IF
      IF (.not.allocated(CDOM_LightAtt)) THEN
        allocate ( CDOM_LightAtt(Ngrids) )
      END IF
      IF (.not.allocated(CDOM_sigma)) THEN
        allocate ( CDOM_sigma(Ngrids) )
      END IF
      IF (.not.allocated(DetRR)) THEN
        allocate ( DetRR(Ngrids) )
      END IF
      IF (.not.allocated(K_DIN)) THEN
        allocate ( K_DIN(Ngrids) )
      END IF
      IF (.not.allocated(PARfrac)) THEN
        allocate ( PARfrac(Ngrids) )
      END IF
      IF (.not.allocated(PhotoRmax)) THEN
        allocate ( PhotoRmax(Ngrids) )
      END IF
      IF (.not.allocated(PhyIS)) THEN
        allocate ( PhyIS(Ngrids) )
      END IF
      IF (.not.allocated(PhyMRD)) THEN
        allocate ( PhyMRD(Ngrids) )
      END IF
      IF (.not.allocated(PhyMRN)) THEN
        allocate ( PhyMRN(Ngrids) )
      END IF
      IF (.not.allocated(ThetaM)) THEN
        allocate ( ThetaM(Ngrids) )
      END IF
      IF (.not.allocated(Vm_DIN)) THEN
        allocate ( Vm_DIN(Ngrids) )
      END IF
      IF (.not.allocated(wDet)) THEN
        allocate ( wDet(Ngrids) )
      END IF
      IF (.not.allocated(wPhy)) THEN
        allocate ( wPhy(Ngrids) )
      END IF

#ifdef TANGENT
      IF (.not.allocated(tl_PARfrac)) THEN
        allocate ( tl_PARfrac(Ngrids) )
      END IF
#endif
#ifdef ADJOINT
      IF (.not.allocated(ad_PARfrac)) THEN
        allocate ( ad_PARfrac(Ngrids) )
      END IF
#endif
!
!  Allocate biological tracer vector.
!
      IF (.not.allocated(idbio)) THEN
        allocate ( idbio(NBT) )
      END IF
!
!-----------------------------------------------------------------------
!  Initialize tracer identification indices.
!-----------------------------------------------------------------------
!
      ic=NAT+NPT+NCS+NNS
      DO i=1,NBT
        idbio(i)=ic+i
      END DO
      ic=ic+1
      iDIN_=ic
      DO i=1,NBands
        iAphy(i)=ic+1
        ic=ic+1
        iBphy(i)=ic+1
        ic=ic+1
        iBdet(i)=ic+1
        ic=ic+1
        aCDOM(i)=ic+1
        ic=ic+1
      END DO
      i440n=1
      i510n=2

      RETURN
      END SUBROUTINE initialize_biology
