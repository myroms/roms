#include "cppdefs.h"
      MODULE mod_average
#ifdef AVERAGES
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  2D Time-averaged fields for output purposes.                        !
!                                                                      !
!  avgu2d     2D velocity component (m/s) in the XI-direction.         !
!  avgv2d     2D velocity component (m/s) in the ETA-direction.        !
!  avgzeta    Free surface (m).                                        !
!                                                                      !
!  3D Time-averaged fields for output purposes.                        !
!                                                                      !
!  avgAKs     Vertical diffusion of Salinity (m2/s).                   !
!  avgAKt     Vertical diffusion of temperature (m2/s).                !
!  avgAKv     Vertical viscosity (m2/s).                               !
!  avghbbl    Depth of oceanic bottom boundary layer (m).              !
!  avghsbl    Depth of oceanic surface boundary layer (m).             !
!  avgrho     Density anomaly (kg/m3).                                 !
!  avgt       Tracer type variables (usually, potential temperature    !
!               and salinity).                                         !
!  avgu3d     3D velocity component (m/s) in the XI-direction.         !
!  avgv3d     3D velocity component (m/s) in the ETA-direction.        !
!  avgw3d     S-coordinate [omega*Hz/mn] vertical velocity (m3/s).     !
!                                                                      !
!=======================================================================
!
        USE mod_kinds

        implicit none

        TYPE T_AVERAGE

          real(r8), pointer :: avgu2d(:,:)
          real(r8), pointer :: avgv2d(:,:)
          real(r8), pointer :: avgzeta(:,:)

# ifdef SOLVE3D
          real(r8), pointer :: avgrho(:,:,:)
          real(r8), pointer :: avgt(:,:,:,:)
          real(r8), pointer :: avgu3d(:,:,:)
          real(r8), pointer :: avgv3d(:,:,:)
          real(r8), pointer :: avgw3d(:,:,:)
#  ifdef AVERAGES_AKS
          real(r8), pointer :: avgAKs(:,:,:)
#  endif
#  ifdef AVERAGES_AKT
          real(r8), pointer :: avgAKt(:,:,:)
#  endif
#  ifdef AVERAGES_AKV
          real(r8), pointer :: avgAKv(:,:,:)
#  endif
#  ifdef LMD_BKPP
          real(r8), pointer :: avghbbl(:,:)
#  endif
#  ifdef LMD_SKPP
          real(r8), pointer :: avghsbl(:,:)
#  endif
# endif

        END TYPE T_AVERAGE

        TYPE (T_AVERAGE), allocatable :: AVERAGE(:)

      CONTAINS

      SUBROUTINE allocate_average (ng, LBi, UBi, LBj, UBj)
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
      IF (ng.eq.1 ) allocate ( AVERAGE(Ngrids) )
!
      allocate ( AVERAGE(ng) % avgu2d(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE(ng) % avgv2d(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE(ng) % avgzeta(LBi:UBi,LBj:UBj) )
# ifdef SOLVE3D
      allocate ( AVERAGE(ng) % avgrho(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( AVERAGE(ng) % avgt(LBi:UBi,LBj:UBj,N(ng),NT(ng)) )
      allocate ( AVERAGE(ng) % avgu3d(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( AVERAGE(ng) % avgv3d(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( AVERAGE(ng) % avgw3d(LBi:UBi,LBj:UBj,0:N(ng)) )
#  ifdef AVERAGES_AKS
      allocate ( AVERAGE(ng) % avgAKs(LBi:UBi,LBj:UBj,0:N(ng)) )
#  endif
#  ifdef AVERAGES_AKT
      allocate ( AVERAGE(ng) % avgAKt(LBi:UBi,LBj:UBj,0:N(ng)) )
#  endif
#  ifdef AVERAGES_AKV
      allocate ( AVERAGE(ng) % avgAKv(LBi:UBi,LBj:UBj,0:N(ng)) )
#  endif

#  ifdef LMD_BKPP
      allocate ( AVERAGE(ng) % avghbbl(LBi:UBi,LBj:UBj) )
#  endif
#  ifdef LMD_SKPP
      allocate ( AVERAGE(ng) % avghsbl(LBi:UBi,LBj:UBj) )
#  endif
# endif

      RETURN
      END SUBROUTINE allocate_average

      SUBROUTINE initialize_average (ng, tile)
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
      AVERAGE(ng) % avgu2d(IstrR:IendR,JstrR:JendR) = IniVal
      AVERAGE(ng) % avgv2d(IstrR:IendR,JstrR:JendR) = IniVal
      AVERAGE(ng) % avgzeta(IstrR:IendR,JstrR:JendR) = IniVal

# ifdef SOLVE3D

      AVERAGE(ng) % avgrho(IstrR:IendR,JstrR:JendR,1:N(ng)) = IniVal
      AVERAGE(ng) % avgu3d(IstrR:IendR,JstrR:JendR,1:N(ng)) = IniVal
      AVERAGE(ng) % avgv3d(IstrR:IendR,JstrR:JendR,1:N(ng)) = IniVal
      AVERAGE(ng) % avgw3d(IstrR:IendR,JstrR:JendR,0:N(ng)) = IniVal
      AVERAGE(ng) % avgt(IstrR:IendR,JstrR:JendR,1:N(ng),1:NT(ng)) =    &
     &                  IniVal
#   ifdef AVERAGES_AKS
      AVERAGE(ng) % avgAKs(IstrR:IendR,JstrR:JendR,0:N(ng)) = IniVal
#   endif
#   ifdef AVERAGES_AKT
      AVERAGE(ng) % avgAKt(IstrR:IendR,JstrR:JendR,0:N(ng)) = IniVal
#   endif
#   ifdef AVERAGES_AKV
      AVERAGE(ng) % avgAKv(IstrR:IendR,JstrR:JendR,0:N(ng)) = IniVal
#   endif
#  ifdef LMD_BKPP
      AVERAGE(ng) % avghbbl(IstrR:IendR,JstrR:JendR) = IniVal
#  endif
#  ifdef LMD_SKPP
      AVERAGE(ng) % avghsbl(IstrR:IendR,JstrR:JendR) = IniVal
#  endif

# endif
      RETURN
      END SUBROUTINE initialize_average
#endif
      END MODULE mod_average

