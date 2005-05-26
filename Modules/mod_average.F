#include "cppdefs.h"
      MODULE mod_average
#ifdef AVERAGES
!
!=======================================================================
!  Copyright (c) 2005 ROMS/TOMS Group                                  !
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
!  avgbedldu  Bed load flux u-direction (kg/m2/s).                     !
!  avgbedldv  Bed load flux v-direction (kg/m2/s).                     !
!  avglhf     Latent heat flux (W/m2).                                 !
!  avglrf     Longwave radiation flux (W/m2).                          !
!  avghbbl    Depth of oceanic bottom boundary layer (m).              !
!  avghsbl    Depth of oceanic surface boundary layer (m).             !
!  avgrho     Density anomaly (kg/m3).                                 !
!  avgshf     Sensible heat flux (W/m2).                               !
!  avgsrf     Shortwave radiation flux (W/m2).                         !
!  avgstf     Surface net heat flux (W/m2).                            !
!  avgswf     Surface net salt flux (kg/m2/s).                         !
!  avgevap    Surface net evaporation (kg/m2/s).                       !
!  avgrain    Surface net rain fall (kg/m2/s).                         !
!  avgsus     Surface u-momentum stress (N/m2).                        !
!  avgsvs     Surface v-momentum stress (N/m2).                        !
!  avgt       Tracer type variables (usually, potential temperature    !
!               and salinity).                                         !
!  avgUT      Quadratic term <u*t> for potential temperature and       !
!               salinity at U-points.                                  !
!  avgVT      Quadratic term <v*t> for potential temperature and       !
!               salinity at V-points.                                  !
!  avgTT      Quadratic term <t*t> for tracers.                        !
!  avgUU      Quadratic term <u*u> for 3D momentum at U-points.        !
!  avgUV      Quadratic term <u*v> for 3D momentum at RHO-points.      !
!  avgVV      Quadratic term <v*v> for 3D momentum at V-points.        !
!  avgU2      Quadratic term <ubar*ubar> for 2D momentum at U-points.  !
!  avgV2      Quadratic term <vbar*vbar> for 2D momentum at V-points.  !
!  avgZZ      Quadratic term <zeta*zeta> for free-surface.             !
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
# ifdef AVERAGES_QUADRATIC
          real(r8), pointer :: avgU2(:,:)
          real(r8), pointer :: avgV2(:,:)
          real(r8), pointer :: avgZZ(:,:)
# endif
# ifdef SOLVE3D
          real(r8), pointer :: avgrho(:,:,:)
          real(r8), pointer :: avgt(:,:,:,:)
          real(r8), pointer :: avgu3d(:,:,:)
          real(r8), pointer :: avgv3d(:,:,:)
          real(r8), pointer :: avgw3d(:,:,:)
#  ifdef AVERAGES_QUADRATIC
#   ifdef SALINITY
          real(r8), pointer :: avgSS(:,:,:)
#   endif
          real(r8), pointer :: avgTT(:,:,:,:)
          real(r8), pointer :: avgUT(:,:,:,:)
          real(r8), pointer :: avgVT(:,:,:,:)
          real(r8), pointer :: avgUU(:,:,:)
          real(r8), pointer :: avgUV(:,:,:)
          real(r8), pointer :: avgVV(:,:,:)
#  endif
#  ifdef AVERAGES_AKS
          real(r8), pointer :: avgAKs(:,:,:)
#  endif
#  ifdef AVERAGES_AKT
          real(r8), pointer :: avgAKt(:,:,:)
#  endif
#  ifdef AVERAGES_AKV
          real(r8), pointer :: avgAKv(:,:,:)
#  endif
#  ifdef AVERAGES_FLUXES
          real(r8), pointer :: avgstf(:,:)
          real(r8), pointer :: avgswf(:,:)
#   ifdef BULK_FLUXES
          real(r8), pointer :: avglhf(:,:)
          real(r8), pointer :: avglrf(:,:)
          real(r8), pointer :: avgshf(:,:)
#    ifdef EMINUSP
          real(r8), pointer :: avgevap(:,:)
          real(r8), pointer :: avgrain(:,:)
#    endif
#   endif
#   ifdef SHORTWAVE
          real(r8), pointer :: avgsrf(:,:)
#   endif
#  endif
#  ifdef LMD_BKPP
          real(r8), pointer :: avghbbl(:,:)
#  endif
#  ifdef LMD_SKPP
          real(r8), pointer :: avghsbl(:,:)
#  endif
# endif
# ifdef AVERAGES_FLUXES
          real(r8), pointer :: avgsus(:,:)
          real(r8), pointer :: avgsvs(:,:)
# endif
# if defined SEDIMENT && defined BEDLOAD
          real(r8), pointer :: avgbedldu(:,:,:)
          real(r8), pointer :: avgbedldv(:,:,:)
# endif

        END TYPE T_AVERAGE

        TYPE (T_AVERAGE), allocatable :: AVERAGE(:)

      CONTAINS

      SUBROUTINE allocate_average (ng, LBi, UBi, LBj, UBj)
!
!=======================================================================
!  Copyright (c) 2005 ROMS/TOMS Group                                  !
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
# ifdef AVERAGES_QUADRATIC
      allocate ( AVERAGE(ng) % avgU2(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE(ng) % avgV2(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE(ng) % avgZZ(LBi:UBi,LBj:UBj) )
# endif
# ifdef SOLVE3D
      allocate ( AVERAGE(ng) % avgrho(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( AVERAGE(ng) % avgt(LBi:UBi,LBj:UBj,N(ng),NT(ng)) )
      allocate ( AVERAGE(ng) % avgu3d(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( AVERAGE(ng) % avgv3d(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( AVERAGE(ng) % avgw3d(LBi:UBi,LBj:UBj,0:N(ng)) )
#  ifdef AVERAGES_QUADRATIC
      allocate ( AVERAGE(ng) % avgTT(LBi:UBi,LBj:UBj,N(ng),NAT) )
      allocate ( AVERAGE(ng) % avgUT(LBi:UBi,LBj:UBj,N(ng),NAT) )
      allocate ( AVERAGE(ng) % avgVT(LBi:UBi,LBj:UBj,N(ng),NAT) )
      allocate ( AVERAGE(ng) % avgUU(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( AVERAGE(ng) % avgUV(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( AVERAGE(ng) % avgVV(LBi:UBi,LBj:UBj,N(ng)) )
#  endif
#  ifdef AVERAGES_AKS
      allocate ( AVERAGE(ng) % avgAKs(LBi:UBi,LBj:UBj,0:N(ng)) )
#  endif
#  ifdef AVERAGES_AKT
      allocate ( AVERAGE(ng) % avgAKt(LBi:UBi,LBj:UBj,0:N(ng)) )
#  endif
#  ifdef AVERAGES_AKV
      allocate ( AVERAGE(ng) % avgAKv(LBi:UBi,LBj:UBj,0:N(ng)) )
#  endif
#  ifdef AVERAGES_FLUXES
      allocate ( AVERAGE(ng) % avgstf(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE(ng) % avgswf(LBi:UBi,LBj:UBj) )
#   ifdef BULK_FLUXES
      allocate ( AVERAGE(ng) % avglhf(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE(ng) % avglrf(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE(ng) % avgshf(LBi:UBi,LBj:UBj) )
#    ifdef EMINUSP
      allocate ( AVERAGE(ng) % avgevap(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE(ng) % avgrain(LBi:UBi,LBj:UBj) )
#    endif
#   endif
#   ifdef SHORTWAVE
      allocate ( AVERAGE(ng) % avgsrf(LBi:UBi,LBj:UBj) )
#   endif
#  endif

#  ifdef LMD_BKPP
      allocate ( AVERAGE(ng) % avghbbl(LBi:UBi,LBj:UBj) )
#  endif
#  ifdef LMD_SKPP
      allocate ( AVERAGE(ng) % avghsbl(LBi:UBi,LBj:UBj) )
#  endif
# endif
# ifdef AVERAGES_FLUXES
      allocate ( AVERAGE(ng) % avgsus(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE(ng) % avgsvs(LBi:UBi,LBj:UBj) )
# endif

# if defined SEDIMENT && defined BEDLOAD
      allocate ( AVERAGE(ng) % avgbedldu(LBi:UBi,LBj:UBj,NST) )
      allocate ( AVERAGE(ng) % avgbedldv(LBi:UBi,LBj:UBj,NST) )
# endif


      RETURN
      END SUBROUTINE allocate_average

      SUBROUTINE initialize_average (ng, tile)
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
#if defined SEDIMENT || defined BBL_MODEL
      USE mod_sediment
#endif
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
      AVERAGE(ng) % avgu2d(Imin:Imax,Jmin:Jmax) = IniVal
      AVERAGE(ng) % avgv2d(Imin:Imax,Jmin:Jmax) = IniVal
      AVERAGE(ng) % avgzeta(Imin:Imax,Jmin:Jmax) = IniVal
# ifdef AVERAGES_QUADRATIC
      AVERAGE(ng) % avgU2(Imin:Imax,Jmin:Jmax) = IniVal
      AVERAGE(ng) % avgV2(Imin:Imax,Jmin:Jmax) = IniVal
      AVERAGE(ng) % avgZZ(Imin:Imax,Jmin:Jmax) = IniVal
# endif

# ifdef SOLVE3D
      AVERAGE(ng) % avgrho(Imin:Imax,Jmin:Jmax,1:N(ng)) = IniVal
      AVERAGE(ng) % avgu3d(Imin:Imax,Jmin:Jmax,1:N(ng)) = IniVal
      AVERAGE(ng) % avgv3d(Imin:Imax,Jmin:Jmax,1:N(ng)) = IniVal
      AVERAGE(ng) % avgw3d(Imin:Imax,Jmin:Jmax,0:N(ng)) = IniVal
      AVERAGE(ng) % avgt(Imin:Imax,Jmin:Jmax,1:N(ng),1:NT(ng)) = IniVal
#  ifdef AVERAGES_QUADRATIC
      AVERAGE(ng) % avgTT(Imin:Imax,Jmin:Jmax,1:N(ng),1:NAT) = IniVal
      AVERAGE(ng) % avgUT(Imin:Imax,Jmin:Jmax,1:N(ng),1:NAT) = IniVal
      AVERAGE(ng) % avgVT(Imin:Imax,Jmin:Jmax,1:N(ng),1:NAT) = IniVal
      AVERAGE(ng) % avgUU(Imin:Imax,Jmin:Jmax,1:N(ng)) = IniVal
      AVERAGE(ng) % avgUV(Imin:Imax,Jmin:Jmax,1:N(ng)) = IniVal
      AVERAGE(ng) % avgVV(Imin:Imax,Jmin:Jmax,1:N(ng)) = IniVal
#  endif
#  ifdef AVERAGES_AKS
      AVERAGE(ng) % avgAKs(Imin:Imax,Jmin:Jmax,0:N(ng)) = IniVal
#  endif
#  ifdef AVERAGES_AKT
      AVERAGE(ng) % avgAKt(Imin:Imax,Jmin:Jmax,0:N(ng)) = IniVal
#  endif
#  ifdef AVERAGES_AKV
      AVERAGE(ng) % avgAKv(Imin:Imax,Jmin:Jmax,0:N(ng)) = IniVal
#  endif
#  ifdef AVERAGES_FLUXES
      AVERAGE(ng) % avgstf(Imin:Imax,Jmin:Jmax) = IniVal
      AVERAGE(ng) % avgswf(Imin:Imax,Jmin:Jmax) = IniVal
#   ifdef BULK_FLUXES
      AVERAGE(ng) % avglhf(Imin:Imax,Jmin:Jmax) = IniVal
      AVERAGE(ng) % avglrf(Imin:Imax,Jmin:Jmax) = IniVal
      AVERAGE(ng) % avgshf(Imin:Imax,Jmin:Jmax) = IniVal
#    ifdef EMINUSP
      AVERAGE(ng) % avgevap(Imin:Imax,Jmin:Jmax) = IniVal
      AVERAGE(ng) % avgrain(Imin:Imax,Jmin:Jmax) = IniVal
#    endif
#   endif
#   ifdef SHORTWAVE
      AVERAGE(ng) % avgsrf(Imin:Imax,Jmin:Jmax) = IniVal
#   endif
#  endif
#  ifdef LMD_BKPP
      AVERAGE(ng) % avghbbl(Imin:Imax,Jmin:Jmax) = IniVal
#  endif
#  ifdef LMD_SKPP
      AVERAGE(ng) % avghsbl(Imin:Imax,Jmin:Jmax) = IniVal
#  endif

# endif

# ifdef AVERAGES_FLUXES
      AVERAGE(ng) % avgsus(Imin:Imax,Jmin:Jmax) = IniVal
      AVERAGE(ng) % avgsvs(Imin:Imax,Jmin:Jmax) = IniVal
# endif

# if defined SEDIMENT && defined BEDLOAD
      AVERAGE(ng) % avgbedldu(LBi:UBi,LBj:UBj,NST) = IniVal
      AVERAGE(ng) % avgbedldv(LBi:UBi,LBj:UBj,NST) = IniVal
# endif

      RETURN
      END SUBROUTINE initialize_average
#endif
      END MODULE mod_average
