#include "cppdefs.h"
      MODULE mod_coupling
#ifdef SOLVE3D
!
!=======================================================================
!  Copyright (c) 2005 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  DU_avg1   Time averaged U-flux for 2D equations (m3/s).             !
!  DU_avg2   Time averaged U-flux for 3D equations coupling (m3/s).    !
!  DV_avg1   Time averaged V-flux for 2D equations (m3/s).             !
!  DV_avg2   Time averaged V-flux for 3D equations coupling (m3/s).    !
!  Zt_avg1   Free-surface averaged over all short time-steps (m).      !
!  rhoA      Normalized vertical averaged density.                     !
!  rhoS      Normalized vertical averaged density perturbation.        !
!  rufrc     Right-hand-side forcing term for 2D U-momentum (m4/s2)    !
!  rvfrc     Right-hand-side forcing term for 2D V-momentum (m4/s2)    !
!                                                                      !
!=======================================================================
!
        USE mod_kinds

        implicit none

        TYPE T_COUPLING
!
!  Nonlinear model state.
!
          real(r8), pointer :: DU_avg1(:,:)
          real(r8), pointer :: DU_avg2(:,:)
          real(r8), pointer :: DV_avg1(:,:)
          real(r8), pointer :: DV_avg2(:,:)
          real(r8), pointer :: Zt_avg1(:,:)
          real(r8), pointer :: rufrc(:,:)
          real(r8), pointer :: rvfrc(:,:)
# ifdef VAR_RHO_2D
          real(r8), pointer :: rhoA(:,:)
          real(r8), pointer :: rhoS(:,:)
# endif

# if defined TANGENT || defined TL_IOMS
!
!  Tangent linear model state.
!
          real(r8), pointer :: tl_DU_avg1(:,:)
          real(r8), pointer :: tl_DU_avg2(:,:)
          real(r8), pointer :: tl_DV_avg1(:,:)
          real(r8), pointer :: tl_DV_avg2(:,:)
          real(r8), pointer :: tl_Zt_avg1(:,:)
          real(r8), pointer :: tl_rufrc(:,:)
          real(r8), pointer :: tl_rvfrc(:,:)
#  ifdef VAR_RHO_2D
          real(r8), pointer :: tl_rhoA(:,:)
          real(r8), pointer :: tl_rhoS(:,:)
#  endif
# endif

# ifdef ADJOINT
!
!  Adjoint model state.
!
          real(r8), pointer :: ad_DU_avg1(:,:)
          real(r8), pointer :: ad_DU_avg2(:,:)
          real(r8), pointer :: ad_DV_avg1(:,:)
          real(r8), pointer :: ad_DV_avg2(:,:)
          real(r8), pointer :: ad_Zt_avg1(:,:)
          real(r8), pointer :: ad_rufrc(:,:)
          real(r8), pointer :: ad_rvfrc(:,:)
#  ifdef VAR_RHO_2D
          real(r8), pointer :: ad_rhoA(:,:)
          real(r8), pointer :: ad_rhoS(:,:)
#  endif
# endif

# if defined FORWARD_READ && \
    (defined TANGENT || defined TL_IOMS || defined ADJOINT)
!
!  Latest two records of the nonlinear trajectory used to interpolate
!  the background state in the tangent linear and adjoint models.
!
          real(r8), pointer :: DU_avg1G(:,:,:)
          real(r8), pointer :: DU_avg2G(:,:,:)
          real(r8), pointer :: DV_avg1G(:,:,:)
          real(r8), pointer :: DV_avg2G(:,:,:)
          real(r8), pointer :: rufrcG(:,:,:)
          real(r8), pointer :: rvfrcG(:,:,:)
# endif

        END TYPE T_COUPLING

        TYPE (T_COUPLING), allocatable :: COUPLING(:)

        CONTAINS

      SUBROUTINE allocate_coupling (ng, LBi, UBi, LBj, UBj)
!
!=======================================================================
!  Copyright (c) 2005 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine allocates all variables in the module for all nested   !
!   grids.                                                             !
!                                                                      !
!=======================================================================
!
      USE mod_param
!
!  Local variable declarations.
!
      integer, intent(in) :: ng, LBi, UBi, LBj, UBj
!
!-----------------------------------------------------------------------
!  Initialize module variables.
!-----------------------------------------------------------------------
!
      IF (ng.eq.1) allocate ( COUPLING(Ngrids) )
!
!  Nonlinear model state.
!
      allocate ( COUPLING(ng) % DU_avg1(LBi:UBi,LBj:UBj) )
      allocate ( COUPLING(ng) % DU_avg2(LBi:UBi,LBj:UBj) )
      allocate ( COUPLING(ng) % DV_avg1(LBi:UBi,LBj:UBj) )
      allocate ( COUPLING(ng) % DV_avg2(LBi:UBi,LBj:UBj) )
      allocate ( COUPLING(ng) % Zt_avg1(LBi:UBi,LBj:UBj) )
      allocate ( COUPLING(ng) % rufrc(LBi:UBi,LBj:UBj) )
      allocate ( COUPLING(ng) % rvfrc(LBi:UBi,LBj:UBj) )

# ifdef VAR_RHO_2D
      allocate ( COUPLING(ng) % rhoA(LBi:UBi,LBj:UBj) )
      allocate ( COUPLING(ng) % rhoS(LBi:UBi,LBj:UBj) )
# endif

# if defined TANGENT || defined TL_IOMS
!
!  Tangent linear model state.
!
      allocate ( COUPLING(ng) % tl_DU_avg1(LBi:UBi,LBj:UBj) )
      allocate ( COUPLING(ng) % tl_DU_avg2(LBi:UBi,LBj:UBj) )
      allocate ( COUPLING(ng) % tl_DV_avg1(LBi:UBi,LBj:UBj) )
      allocate ( COUPLING(ng) % tl_DV_avg2(LBi:UBi,LBj:UBj) )
      allocate ( COUPLING(ng) % tl_Zt_avg1(LBi:UBi,LBj:UBj) )
      allocate ( COUPLING(ng) % tl_rufrc(LBi:UBi,LBj:UBj) )
      allocate ( COUPLING(ng) % tl_rvfrc(LBi:UBi,LBj:UBj) )

#  ifdef VAR_RHO_2D
      allocate ( COUPLING(ng) % tl_rhoA(LBi:UBi,LBj:UBj) )
      allocate ( COUPLING(ng) % tl_rhoS(LBi:UBi,LBj:UBj) )
#  endif
# endif

# ifdef ADJOINT
!
!  Adjoint model state.
!
      allocate ( COUPLING(ng) % ad_DU_avg1(LBi:UBi,LBj:UBj) )
      allocate ( COUPLING(ng) % ad_DU_avg2(LBi:UBi,LBj:UBj) )
      allocate ( COUPLING(ng) % ad_DV_avg1(LBi:UBi,LBj:UBj) )
      allocate ( COUPLING(ng) % ad_DV_avg2(LBi:UBi,LBj:UBj) )
      allocate ( COUPLING(ng) % ad_Zt_avg1(LBi:UBi,LBj:UBj) )
      allocate ( COUPLING(ng) % ad_rufrc(LBi:UBi,LBj:UBj) )
      allocate ( COUPLING(ng) % ad_rvfrc(LBi:UBi,LBj:UBj) )

#  ifdef VAR_RHO_2D
      allocate ( COUPLING(ng) % ad_rhoA(LBi:UBi,LBj:UBj) )
      allocate ( COUPLING(ng) % ad_rhoS(LBi:UBi,LBj:UBj) )
#  endif
# endif

# if defined FORWARD_READ && \
    (defined TANGENT || defined TL_IOMS || defined ADJOINT)
!
!  Latest two records of the nonlinear trajectory used to interpolate
!  the background state in the tangent linear and adjoint models.
!
      allocate ( COUPLING(ng) % DU_avg1G(LBi:UBi,LBj:UBj,2) )
      allocate ( COUPLING(ng) % DU_avg2G(LBi:UBi,LBj:UBj,2) )
      allocate ( COUPLING(ng) % DV_avg1G(LBi:UBi,LBj:UBj,2) )
      allocate ( COUPLING(ng) % DV_avg2G(LBi:UBi,LBj:UBj,2) )
      allocate ( COUPLING(ng) % rufrcG(LBi:UBi,LBj:UBj,2) )
      allocate ( COUPLING(ng) % rvfrcG(LBi:UBi,LBj:UBj,2) )
# endif

      RETURN
      END SUBROUTINE allocate_coupling

      SUBROUTINE initialize_coupling (ng, tile)
!
!=======================================================================
!  Copyright (c) 2005 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine initialize all variables in the module using first     !
!  touch distribution policy. In shared-memory configuration, this     !
!  operation actually performs propagation of the  "shared arrays"     !
!  across the cluster, unless another policy is specified to           !
!  override the default.                                               !
!                                                                      !
!=======================================================================
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
      integer :: Imin, Imax, Jmin, Jmax

      real(r8), parameter :: IniVal = 0.0_r8
!
# include "tile.h"
!
!  Set array initialization range.
!
#ifdef _OPENMP
      IF (WESTERN_EDGE) THEN
        Imin=LBi
      ELSE
        Imin=Istr
      END IF
      IF (EASTERN_EDGE) THEN
        Imax=UBi
      ELSE
        Imax=Iend
      END IF
      IF (SOUTHERN_EDGE) THEN
        Jmin=LBj
      ELSE
        Jmin=Jstr
      END IF
      IF (NORTHERN_EDGE) THEN
        Jmax=UBj
      ELSE
        Jmax=Jend
      END IF
#else
      Imin=LBi
      Imax=UBi
      Jmin=LBj
      Jmax=UBj
#endif
!
!-----------------------------------------------------------------------
!  Initialize module variables.
!-----------------------------------------------------------------------
!
!  Nonlinear model state.
!
      COUPLING(ng) % DU_avg1(Imin:Imax,Jmin:Jmax) = IniVal
      COUPLING(ng) % DU_avg2(Imin:Imax,Jmin:Jmax) = IniVal

      COUPLING(ng) % DV_avg1(Imin:Imax,Jmin:Jmax) = IniVal
      COUPLING(ng) % DV_avg2(Imin:Imax,Jmin:Jmax) = IniVal

      COUPLING(ng) % Zt_avg1(Imin:Imax,Jmin:Jmax) = IniVal

      COUPLING(ng) % rufrc(Imin:Imax,Jmin:Jmax) = IniVal
      COUPLING(ng) % rvfrc(Imin:Imax,Jmin:Jmax) = IniVal

# ifdef VAR_RHO_2D
      COUPLING(ng) % rhoA(Imin:Imax,Jmin:Jmax) = IniVal
      COUPLING(ng) % rhoS(Imin:Imax,Jmin:Jmax) = IniVal
# endif

# if defined TANGENT || defined TL_IOMS
!
!  Tangent linear model state.
!
      COUPLING(ng) % tl_DU_avg1(Imin:Imax,Jmin:Jmax) = IniVal
      COUPLING(ng) % tl_DU_avg2(Imin:Imax,Jmin:Jmax) = IniVal

      COUPLING(ng) % tl_DV_avg1(Imin:Imax,Jmin:Jmax) = IniVal
      COUPLING(ng) % tl_DV_avg2(Imin:Imax,Jmin:Jmax) = IniVal

      COUPLING(ng) % tl_Zt_avg1(Imin:Imax,Jmin:Jmax) = IniVal

      COUPLING(ng) % tl_rufrc(Imin:Imax,Jmin:Jmax) = IniVal
      COUPLING(ng) % tl_rvfrc(Imin:Imax,Jmin:Jmax) = IniVal

#  ifdef VAR_RHO_2D
      COUPLING(ng) % tl_rhoA(Imin:Imax,Jmin:Jmax) = IniVal
      COUPLING(ng) % tl_rhoS(Imin:Imax,Jmin:Jmax) = IniVal
#  endif
# endif

# ifdef ADJOINT
!
!  Adjoint model state.
!
      COUPLING(ng) % ad_DU_avg1(Imin:Imax,Jmin:Jmax) = IniVal
      COUPLING(ng) % ad_DU_avg2(Imin:Imax,Jmin:Jmax) = IniVal

      COUPLING(ng) % ad_DV_avg1(Imin:Imax,Jmin:Jmax) = IniVal
      COUPLING(ng) % ad_DV_avg2(Imin:Imax,Jmin:Jmax) = IniVal

      COUPLING(ng) % ad_Zt_avg1(Imin:Imax,Jmin:Jmax) = IniVal

      COUPLING(ng) % ad_rufrc(Imin:Imax,Jmin:Jmax) = IniVal
      COUPLING(ng) % ad_rvfrc(Imin:Imax,Jmin:Jmax) = IniVal

#  ifdef VAR_RHO_2D
      COUPLING(ng) % ad_rhoA(Imin:Imax,Jmin:Jmax) = IniVal
      COUPLING(ng) % ad_rhoS(Imin:Imax,Jmin:Jmax) = IniVal
#  endif
# endif

# if defined FORWARD_READ && \
    (defined TANGENT || defined TL_IOMS || defined ADJOINT)
!
!  Latest two records of the nonlinear trajectory used to interpolate
!  the background state in the tangent linear and adjoint models.
!
      COUPLING(ng) % DU_avg1G(Imin:Imax,Jmin:Jmax,1:2) = IniVal
      COUPLING(ng) % DU_avg2G(Imin:Imax,Jmin:Jmax,1:2) = IniVal

      COUPLING(ng) % DV_avg1G(Imin:Imax,Jmin:Jmax,1:2) = IniVal
      COUPLING(ng) % DV_avg2G(Imin:Imax,Jmin:Jmax,1:2) = IniVal

      COUPLING(ng) % rufrcG(Imin:Imax,Jmin:Jmax,1:2) = IniVal
      COUPLING(ng) % rvfrcG(Imin:Imax,Jmin:Jmax,1:2) = IniVal
# endif

      RETURN
      END SUBROUTINE initialize_coupling
#endif
      END MODULE mod_coupling
