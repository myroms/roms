#include "cppdefs.h"
      MODULE mod_bbl
#ifdef BBL_MODEL
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  Abed         wind-induced, bed wave excursion amplitude (m).        !
!  Hripple      Bed ripple height (m).                                 !
!  Lripple      Bed ripple length (m).                                 !
!  Sdens        Sediment grain density (kg/m3) at RHO-points.          !
!  Ssize        Sediment grain diameter size (m) at RHO-points.        !
!  Ubed         Wind-induced, bed wave orbital U-velocity (m/s).       !
!  Ubot         Bottom U-momentum (m/s) above bed.                     !
!  Vbed         Wind-induced, bed wave orbital V-velocity (m/s).       !
!  Vbot         Bottom V-momentum (m/s) above bed.                     !
!  Zbnot        Hydraulic bottom roughness (m).                        !
!  Zbnotc       Apparent hydraulic bottom roughness (m).               !
!  bustrw       Kinematic bottom stress (m2/s2) due to wind-induced    !
!                 waves the XI-direction at horizontal U-points.       !
!  bvstrw       Kinematic bottom stress (m2/s2) due to wind-induced    !
!                 waves the ETA-direction at horizontal V-points.      !
!                                                                      !
!=======================================================================
!
        USE mod_kinds

        implicit none

        TYPE T_BBL

          integer,  pointer :: Iconv(:,:)

          real(r8), pointer :: Abed(:,:)
          real(r8), pointer :: Hripple(:,:)
          real(r8), pointer :: Lripple(:,:)
          real(r8), pointer :: Sdens(:,:)
          real(r8), pointer :: Ssize(:,:)
          real(r8), pointer :: Ubed(:,:)
          real(r8), pointer :: Ubot(:,:)
          real(r8), pointer :: Vbed(:,:)
          real(r8), pointer :: Vbot(:,:)
          real(r8), pointer :: Zbnot(:,:)
          real(r8), pointer :: Zbnotc(:,:)
          real(r8), pointer :: bustrw(:,:)
          real(r8), pointer :: bvstrw(:,:)

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
      allocate ( BBL(ng) % Abed(LBi:UBi,LBj:UBj) )
      allocate ( BBL(ng) % Hripple(LBi:UBi,LBj:UBj) )
      allocate ( BBL(ng) % Lripple(LBi:UBi,LBj:UBj) )
      allocate ( BBL(ng) % Sdens(LBi:UBi,LBj:UBj) )
      allocate ( BBL(ng) % Ssize(LBi:UBi,LBj:UBj) )
      allocate ( BBL(ng) % Ubed(LBi:UBi,LBj:UBj) )
      allocate ( BBL(ng) % Ubot(LBi:UBi,LBj:UBj) )
      allocate ( BBL(ng) % Vbed(LBi:UBi,LBj:UBj) )
      allocate ( BBL(ng) % Vbot(LBi:UBi,LBj:UBj) )
      allocate ( BBL(ng) % Zbnot(LBi:UBi,LBj:UBj) )
      allocate ( BBL(ng) % Zbnotc(LBi:UBi,LBj:UBj) )
      allocate ( BBL(ng) % bustrw(LBi:UBi,LBj:UBj) )
      allocate ( BBL(ng) % bvstrw(LBi:UBi,LBj:UBj) )

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
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV

      real(r8), parameter :: IniVal = 0.0_r8
!
# include "tile.h"
# ifdef DISTRIBUTE
      IstrR=LBi
      IendR=UBi
      JstrR=LBj
      JendR=UBj
# else
#  include "set_bounds.h"
# endif
!
!-----------------------------------------------------------------------
!  Initialize module variables.
!-----------------------------------------------------------------------
!
      BBL(ng) % Iconv(IstrR:IendR,JstrR:JendR) = 0

      BBL(ng) % Abed(IstrR:IendR,JstrR:JendR) = IniVal

      BBL(ng) % Hripple(IstrR:IendR,JstrR:JendR) = IniVal
      BBL(ng) % Lripple(IstrR:IendR,JstrR:JendR) = IniVal

      BBL(ng) % Sdens(IstrR:IendR,JstrR:JendR) = IniVal
      BBL(ng) % Ssize(IstrR:IendR,JstrR:JendR) = IniVal

      BBL(ng) % Ubed(IstrR:IendR,JstrR:JendR) = IniVal
      BBL(ng) % Ubot(IstrR:IendR,JstrR:JendR) = IniVal

      BBL(ng) % Vbed(IstrR:IendR,JstrR:JendR) = IniVal
      BBL(ng) % Vbot(IstrR:IendR,JstrR:JendR) = IniVal

      BBL(ng) % Zbnot(IstrR:IendR,JstrR:JendR) = IniVal
      BBL(ng) % Zbnotc(IstrR:IendR,JstrR:JendR) = IniVal

      BBL(ng) % bustrw(IstrR:IendR,JstrR:JendR) = IniVal
      BBL(ng) % bvstrw(IstrR:IendR,JstrR:JendR) = IniVal

      RETURN
      END SUBROUTINE initialize_bbl
#endif
      END MODULE mod_bbl
