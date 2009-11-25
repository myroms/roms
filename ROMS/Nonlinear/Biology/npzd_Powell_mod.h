!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Parameters for Powell et al. (2006) ecosystem model:                !
!                                                                      !
!  AttPhy    Light attenuation due to phytoplankton (self-shading      !
!              coefficient), [m2/millimole_N].                         !
!  AttSW     Light attenuation due to sea water, [1/m].                !
!  BioIter   Maximum number of iterations to achieve convergence of    !
!              the nonlinear solution.                                 !
!  BioIni    Initial concentration for analytical initial (uniform)    !
!              conditions.                                             !
!  DetRR     Detritus remineraliztion rate, [1/day].                   !
!  K_NO3     Half-saturation for phytoplankton nitrate uptake          !
!              [millimole_N m-3].                                      !
!  Ivlev     Ivlev constant for zooplankton grazing parameterization,  !
!              [nondimensional].                                       !
!  PARfrac   Fraction of shortwave radiation that is available for     !
!              photosyntesis [nondimensional].                         !
!  PhyIS     Phytoplankton, initial slope of the P-I curve [m2/W].     !
!  PhyMRD    Phytoplankton mortality rate to the Detritus pool,        !
!              [1/day].                                                !
!  PhyMRN    Phytoplankton mortality rate to the Nitrogen pool,        !
!              [1/day].                                                !
!  Vm_NO3    Nitrate uptake rate, [1/day].                             !
!  wDet      Detrital sinking rate, [m/day].                           !
!  wPhy      Phytoplankton sinking rate, [m/day].                      !
!  ZooEED    Zooplankton excretion efficiency to Detritus pool,        !
!              {nondimensional].                                       !
!  ZooEEN    Zooplankton excretion efficiency to Nitrogen pool,        !
!              {nondimensional].                                       !
!  ZooGR     Zooplankton grazing rate, [1/day].                        !
!  ZooMRD    Zooplankton mortality rate to Detritus pool, [1/day].     !
!  ZooMRN    Zooplankton mortality rate to Nitrogen pool, [1/day].     !
!                                                                      !
!=======================================================================
!
        USE mod_param
!
        implicit none
!
        integer, dimension(Ngrids) :: BioIter

#ifdef ANA_BIOLOGY
        real(r8), allocatable :: BioIni(:,:)
#endif
        real(r8), dimension(Ngrids) :: AttPhy        ! m2/mmole
        real(r8), dimension(Ngrids) :: AttSW         ! 1/m
        real(r8), dimension(Ngrids) :: DetRR         ! 1/day
        real(r8), dimension(Ngrids) :: K_NO3         ! mmol/m3
        real(r8), dimension(Ngrids) :: Ivlev         ! nondimensional
        real(r8), dimension(Ngrids) :: PARfrac       ! nondimensional
#ifdef TANGENT
        real(r8), dimension(Ngrids) :: tl_PARfrac
#endif
#ifdef ADJOINT
        real(r8), dimension(Ngrids) :: ad_PARfrac
#endif
        real(r8), dimension(Ngrids) :: PhyIS         ! m2/W
        real(r8), dimension(Ngrids) :: PhyMRD        ! 1/day
        real(r8), dimension(Ngrids) :: PhyMRN        ! 1/day
        real(r8), dimension(Ngrids) :: Vm_NO3        ! 1/day
        real(r8), dimension(Ngrids) :: wDet          ! m/day
#ifdef TANGENT
        real(r8), dimension(Ngrids) :: tl_wDet
#endif
#ifdef ADJOINT
        real(r8), dimension(Ngrids) :: ad_wDet
#endif
        real(r8), dimension(Ngrids) :: wPhy          ! m/day
#ifdef TANGENT
        real(r8), dimension(Ngrids) :: tl_wPhy
#endif
#ifdef ADJOINT
        real(r8), dimension(Ngrids) :: ad_wPhy
#endif
        real(r8), dimension(Ngrids) :: ZooEED        ! nondimensional
        real(r8), dimension(Ngrids) :: ZooEEN        ! nondimensional
        real(r8), dimension(Ngrids) :: ZooGR         ! 1/day
        real(r8), dimension(Ngrids) :: ZooMRD        ! 1/day
        real(r8), dimension(Ngrids) :: ZooMRN        ! 1/day
