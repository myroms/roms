#include "cppdefs.h"
      MODULE mod_bbl
#ifdef BBL_MODEL
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  Ubed         Wind-induced, bed wave orbital U-velocity (m/s).       !
!  Ubot         Bottom U-momentum (m/s) above bed.                     !
!  Vbed         Wind-induced, bed wave orbital V-velocity (m/s).       !
!  Vbot         Bottom V-momentum (m/s) above bed.                     !
!  bustrw       Kinematic bottom stress (m2/s2) due to wind-induced    !
!                 waves the XI-direction at horizontal U-points.       !
!  bustrcwmax   Kinematic bottom stress (m2/s2) due to maximum wind    !
!                 and currents in the XI-direction at U-points.        !
!  bvstrw       Kinematic bottom stress (m2/s2) due to wind-induced    !
!                 waves the ETA-direction at horizontal V-points.      !
!  bvstrcwmax   Kinematic bottom stress (m2/s2) due to maximum wind    !
!                 and currents in the ETA-direction at V-points.       !
!                                                                      !
!=======================================================================
!
        USE mod_kinds

        implicit none

        TYPE T_BBL

          integer,  pointer :: Iconv(:,:)

          real(r8), pointer :: Ubed(:,:)
          real(r8), pointer :: Ubot(:,:)
          real(r8), pointer :: Vbed(:,:)
          real(r8), pointer :: Vbot(:,:)
          real(r8), pointer :: bustrw(:,:)
          real(r8), pointer :: bvstrw(:,:)
          real(r8), pointer :: bustrcwmax(:,:)
          real(r8), pointer :: bvstrcwmax(:,:)

        END TYPE T_BBL

        TYPE (T_BBL), allocatable :: BBL(:)

      CONTAINS

      SUBROUTINE allocate_bbl (ng, LBi, UBi, LBj, UBj)
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine allocates all variables in the module for all nested   !
!  grids.                                                              !
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
!  Allocate module variables.
!-----------------------------------------------------------------------
!
      IF (ng.eq.1) allocate ( BBL(Ngrids) )
!
      allocate ( BBL(ng) % Iconv(LBi:UBi,LBj:UBj) )
      allocate ( BBL(ng) % Ubed(LBi:UBi,LBj:UBj) )
      allocate ( BBL(ng) % Ubot(LBi:UBi,LBj:UBj) )
      allocate ( BBL(ng) % Vbed(LBi:UBi,LBj:UBj) )
      allocate ( BBL(ng) % Vbot(LBi:UBi,LBj:UBj) )
      allocate ( BBL(ng) % bustrw(LBi:UBi,LBj:UBj) )
      allocate ( BBL(ng) % bvstrw(LBi:UBi,LBj:UBj) )
      allocate ( BBL(ng) % bustrcwmax(LBi:UBi,LBj:UBj) )
      allocate ( BBL(ng) % bvstrcwmax(LBi:UBi,LBj:UBj) )

      RETURN
      END SUBROUTINE allocate_bbl

      SUBROUTINE initialize_bbl (ng, tile)
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
      BBL(ng) % Iconv(Imin:Imax,Jmin:Jmax) = 0

      BBL(ng) % Ubed(Imin:Imax,Jmin:Jmax) = IniVal
      BBL(ng) % Ubot(Imin:Imax,Jmin:Jmax) = IniVal

      BBL(ng) % Vbed(Imin:Imax,Jmin:Jmax) = IniVal
      BBL(ng) % Vbot(Imin:Imax,Jmin:Jmax) = IniVal

      BBL(ng) % bustrw(Imin:Imax,Jmin:Jmax) = IniVal
      BBL(ng) % bvstrw(Imin:Imax,Jmin:Jmax) = IniVal

      BBL(ng) % bustrcwmax(Imin:Imax,Jmin:Jmax) = IniVal
      BBL(ng) % bvstrcwmax(Imin:Imax,Jmin:Jmax) = IniVal

      RETURN
      END SUBROUTINE initialize_bbl
#endif
      END MODULE mod_bbl
