#include "cppdefs.h"
      MODULE mod_coupling
#ifdef SOLVE3D
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
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

        END TYPE T_COUPLING

        TYPE (T_COUPLING), allocatable :: COUPLING(:)

        CONTAINS

      SUBROUTINE allocate_coupling (ng, LBi, UBi, LBj, UBj)
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
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

      RETURN
      END SUBROUTINE allocate_coupling

      SUBROUTINE initialize_coupling (ng, tile)
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
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

      RETURN
      END SUBROUTINE initialize_coupling
#endif
      END MODULE mod_coupling
