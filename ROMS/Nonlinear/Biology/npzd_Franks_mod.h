!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Parameters for Franks et al. (1986) type model:                     !
!                                                                      !
!  BioIter   Maximum number of iterations to achieve convergence of    !
!              the nonlinear solution.                                 !
!  BioIni    Initial concentration for analytical initial (uniform)    !
!              conditions.                                             !
!  DetRR     Detritus remineraliztion rate, [1/day].                   !
!  K_ext     Light extinction coefficient, [1/m].                      !
!  K_NO3     Inverse half-saturation for phytoplankton nitrate uptake  !
!              [1/(millimole_N m-3)].                                  !
!  K_Phy     Phytoplankton saturation coefficient, [millimole_N m-3].  !
!  PhyMR     Phytoplankton senescence/mortality rate, [1/day].         !
!  Vm_NO3    Nitrate uptake rate, [1/day].                             !
!  wDet      Detrital sinking rate, [m/day].                           !
!  ZooGR     Zooplankton maximum growth rate, [1/day].                 !
!  ZooMR     Zooplankton mortality rate, [1/day].                      !
!  ZooMD     Zooplankton death bits rate, [1/day].                     !
!  ZooGA     Zooplankton grazing inefficiency, [nondimensional].       !
!  ZooEC     Zooplankton excreted fraction, [nondimensional].          !
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
        real(r8), dimension(Ngrids) :: DetRR         ! 1/day
        real(r8), dimension(Ngrids) :: K_ext         ! 1/m
        real(r8), dimension(Ngrids) :: K_NO3         ! 1/(mmol/m3)
        real(r8), dimension(Ngrids) :: K_Phy         ! mmol/m3
        real(r8), dimension(Ngrids) :: PhyMR         ! 1/day
        real(r8), dimension(Ngrids) :: Vm_NO3        ! 1/day
        real(r8), dimension(Ngrids) :: wDet          ! m/day
#ifdef TANGENT
        real(r8), dimension(Ngrids) :: tl_wDet
#endif
#ifdef ADJOINT
        real(r8), dimension(Ngrids) :: ad_wDet
#endif
        real(r8), dimension(Ngrids) :: ZooGR         ! 1/day
        real(r8), dimension(Ngrids) :: ZooMR         ! 1/day
        real(r8), dimension(Ngrids) :: ZooMD         ! 1/day
        real(r8), dimension(Ngrids) :: ZooGA         ! nondimensional
        real(r8), dimension(Ngrids) :: ZooEC         ! nondimensional
