!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2024 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  Parameters for IOP-based, CDOM ecosystem model:                     !
!                                                                      !
!  AttSW          Light attenuation due to sea water (1/m).            !
!  BioIter        Maximum number of iterations to achieve convergence  !
!                   of the nonlinear solution.                         !
!  BioIni         Initial concentration for analytical initial         !
!                   (uniform) conditions.                              !
!  CDOM_LightAtt  Factor for light attenuation due to CDOM             !
!                   (nondimensional).                                  !
!  CDOM_sigma     Light-dependent degradation term for CDOM            !
!                   (1/(Watt m-2 day)).                                !
!  PARfrac        Fraction of shortwave radiation that is available    !
!                   for photosyntesis (nondimensional).                !
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
      integer, allocatable :: idbio(:)    ! Biological tracers
      integer :: aCDOM(NBands)            ! CDOM absorption
      integer :: i440n                    ! 440nm spectral band
      integer :: i510n                    ! 510nm spectral band
!
!  Biological parameters.
!
      integer, allocatable :: BioIter(:)

      real(r8), allocatable :: BioIni(:,:)

      real(r8), allocatable :: AttSW(:)             ! 1/m
      real(r8), allocatable :: CDOM_LightAtt(:)     ! nondimensonal
      real(r8), allocatable :: CDOM_sigma(:)        ! 1/(W m-2 day)
      real(r8), allocatable :: PARfrac(:)           ! nondimensional

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
      NBT=2
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
      IF (.not.allocated(CDOM_LightAtt)) THEN
        allocate ( CDOM_LightAtt(Ngrids) )
      END IF
      IF (.not.allocated(CDOM_sigma)) THEN
        allocate ( CDOM_sigma(Ngrids) )
      END IF
      IF (.not.allocated(PARfrac)) THEN
        allocate ( PARfrac(Ngrids) )
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
      DO i=1,NBands
        aCDOM(i)=ic+1
        ic=ic+1
      END DO
      i440n=1
      i510n=2

      RETURN
      END SUBROUTINE initialize_biology
