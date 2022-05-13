#include "cppdefs.h"
      MODULE mod_average
#if defined AVERAGES    || \
   (defined AD_AVERAGES && defined ADJOINT) || \
   (defined RP_AVERAGES && defined TL_IOMS) || \
   (defined TL_AVERAGES && defined TANGENT)
!
!git $Id$
!svn $Id: mod_average.F 1126 2022-05-13 03:17:28Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2022 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  The strategy here is  to define all possible pointers in the        !
!  time-averaged structure and allocate only those requested by        !
!  the user. This will facilitate a better management of memory.       !
!                                                                      !
!  Time-averaged state variables for output purposes.                  !
!                                                                      !
!  avgu2d     2D velocity  component (m/s) in the XI-direction.        !
!  avgv2d     2D velocity  component (m/s) in the ETA-direction.       !
!  avgu2dE    2D Eastward  component (m/s) at RHO-points.              !
!  avgv2dN    2D Northward component (m/s) at RHO-points.              !
!  avgzeta    Free surface (m).                                        !

# ifdef SOLVE3D
!  avgu3d     3D velocity  component (m/s) in the XI-direction.        !
!  avgv3d     3D velocity  component (m/s) in the ETA-direction.       !
!  avgu3dE    3D Eastward  component (m/s) at RHO-points.              !
!  avgv3dN    3D Northward component (m/s) at RHO-points.              !
!  avgw3d     S-coordinate [omega*Hz/mn] vertical velocity (m3/s).     !
!  avgwvel    3D "true" vertical velocity (m/s).                       !
!  avgrho     Density anomaly (kg/m3).                                 !
!  avgt       Tracer type variables (usually, potential temperature    !
!               and salinity).                                         !

#  if defined LMD_MIXING || defined MY25_MIXING || defined GLS_MIXING
!  avgAKt     Vertical diffusion of temperature (m2/s).                !
!  avgAKv     Vertical viscosity (m2/s).                               !
!  avgAKs     Vertical diffusion of Salinity (m2/s).                   !
#  endif
#  ifdef LMD_SKPP
!  avghsbl    Depth of oceanic surface boundary layer (m).             !
#  endif
#  ifdef LMD_BKPP
!  avghbbl    Depth of oceanic bottom boundary layer (m).              !
#  endif
# endif
# if defined FORWARD_WRITE && defined SOLVE3D
!                                                                      !
!  Time-averaged 2D/3D coupling terms.                                 !
!                                                                      !
!  avgDU_avg1 time-averaged u-flux for 3D momentum coupling.           !
!  avgDU_avg2 time-averaged u-flux for 3D momentum coupling.           !
!  avgDV_avg1 time-averaged v-flux for 3D momentum coupling.           !
!  avgDV_avg2 time-averaged v-flux for 3D momentum coupling.           !
# endif
!                                                                      !
!  Time-averaged surface and bottom fluxes.                            !
!                                                                      !
!  avgsus     Surface u-momentum stress (N/m2).                        !
!  avgsvs     Surface v-momentum stress (N/m2).                        !
!  avgbus     Bottom u-momentum stress (N/m2).                         !
!  avgbvs     Bottom v-momentum stress (N/m2).                         !

# ifdef SOLVE3D
!  avgstf     Surface net heat flux (W/m2).                            !
!  avgswf     Surface net freshwater flux (kg/m2/s).                   !
#  ifdef SHORTWAVE
!  avgsrf     Shortwave radiation flux (W/m2).                         !
#  endif
#  if defined BULK_FLUXES || defined FRC_COUPLING
!  avglhf     Latent heat flux (W/m2).                                 !
!  avglrf     Longwave radiation flux (W/m2).                          !
!  avgshf     Sensible heat flux (W/m2).                               !
#  endif
#  if defined BULK_FLUXES && defined EMINUSP
!  avgevap    Surface net evaporation (kg/m2/s).                       !
!  avgrain    Surface net rain fall (kg/m2/s).                         !
#  endif
# endif

# ifdef NEARSHORE_MELLOR
!  Time-averaged radiation stresses.                                   !
!                                                                      !
!  avgu2Sd    2D stokes velocity component (m/s) in the XI-direction.  !
!  avgv2Sd    2D stokes velocity component (m/s) in the ETA-direction. !
!  avgu2RS    2D radiation stress tensor in the XI-direction.          !
!  avgv2RS    2D radiation stress tensor in the ETA-direction.         !
!  avgSxx2d   2D radiation stress, xx-component.                       !
!  avgSxy2d   2D radiation stress, xy-component.                       !
!  avgSyy2d   2D radiation stress, yy-component.                       !
#  ifdef SOLVE3D
!  avgu3Sd    3D stokes velocity component (m/s) in the XI-direction.  !
!  avgv3Sd    3D stokes velocity component (m/s) in the ETA-direction. !
!  avgu3RS    3D radiation stress tensor in the XI-direction.          !
!  avgv3RS    3D radiation stress tensor in the ETA-direction.         !
!  avgSxx3d   3D radiation stress, xx-component.                       !
!  avgSxy3d   3D radiation stress, xy-component.                       !
!  avgSyy3d   3D radiation stress, yy-component.                       !
#   ifdef NEARSHORE_MELLOR05
!  avgSzx3d   3D radiation stress, zx-component.                       !
!  avgSzy3d   3D radiation stress, zy-component.                       !
#   endif
#  endif
# endif
!                                                                      !
!  Time-averaged quadratic fields.                                     !
!                                                                      !
!  avgZZ      Quadratic term <zeta*zeta> for free-surface.             !
!  avgU2      Quadratic term <ubar*ubar> for 2D momentum at U-points.  !
!  avgV2      Quadratic term <vbar*vbar> for 2D momentum at V-points.  !

# ifdef SOLVE3D
!  avgUU      Quadratic term <u*u> for 3D momentum at U-points.        !
!  avgVV      Quadratic term <v*v> for 3D momentum at V-points.        !
!  avgUV      Quadratic term <u*v> for 3D momentum at RHO-points.      !
!  avgHuon    U-momentum flux, Hz*u/pn (m3/s).                         !
!  avgHvom    V-momentum flux, Hz*v/pm (m3/s).                         !
!  avgTT      Quadratic term <t*t> for tracers.                        !
!  avgUT      Quadratic term <u*t> for potential temperature and       !
!               salinity at U-points.                                  !
!  avgVT      Quadratic term <v*t> for potential temperature and       !
!               salinity at V-points.                                  !
!  avgHuonT   Tracer u-transport, Hz*u*t/pn (Tunits m3/s).             !
!  avgHvomT   Tracer v-transport, Hz*v*t/pn (Tunits m3/s).             !
# endif
!                                                                      !
!  Time-averages vorticity fields.                                     !
!                                                                      !
!  avgpvor2d  2D, vertically integrated, potential vorticity.          !
!  avgrvor2d  2D, vertically integrated, relative vorticity.           !
!  avgpvor3d  3D potential vorticity.                                  !
!  rvorvor2d  3D relative vorticity.                                   !
!                                                                      !
!=======================================================================
!
        USE mod_kinds
!
        implicit none
!
        PUBLIC :: allocate_average
        PUBLIC :: deallocate_average
        PUBLIC :: initialize_average
!
!-----------------------------------------------------------------------
!  Define T_AVERAGE structure.
!-----------------------------------------------------------------------
!
        TYPE T_AVERAGE
!
!  Time-averaged state variables.
!
          real(r8), pointer :: avgzeta(:,:)

          real(r8), pointer :: avgu2d(:,:)
          real(r8), pointer :: avgv2d(:,:)

          real(r8), pointer :: avgu2dE(:,:)
          real(r8), pointer :: avgv2dN(:,:)

# ifdef SOLVE3D
          real(r8), pointer :: avgu3d(:,:,:)
          real(r8), pointer :: avgv3d(:,:,:)

          real(r8), pointer :: avgu3dE(:,:,:)
          real(r8), pointer :: avgv3dN(:,:,:)

          real(r8), pointer :: avgw3d(:,:,:)
          real(r8), pointer :: avgwvel(:,:,:)

          real(r8), pointer :: avgrho(:,:,:)
          real(r8), pointer :: avgt(:,:,:,:)

#  if defined LMD_MIXING || defined MY25_MIXING || defined GLS_MIXING
          real(r8), pointer :: avgAKv(:,:,:)
          real(r8), pointer :: avgAKt(:,:,:)
          real(r8), pointer :: avgAKs(:,:,:)
#  endif
#  ifdef LMD_SKPP
          real(r8), pointer :: avghsbl(:,:)
#  endif
#  ifdef LMD_BKPP
          real(r8), pointer :: avghbbl(:,:)
#  endif
# endif

# if defined FORWARD_WRITE && defined SOLVE3D
!
!  Time-averaged 2D/3D coupling terms.
!
          real(r8), pointer :: avgDU_avg1(:,:)
          real(r8), pointer :: avgDU_avg2(:,:)

          real(r8), pointer :: avgDV_avg1(:,:)
          real(r8), pointer :: avgDV_avg2(:,:)
# endif
!
!  Time-averaged surface and bottom fluxes.
!
          real(r8), pointer :: avgsus(:,:)
          real(r8), pointer :: avgsvs(:,:)

          real(r8), pointer :: avgbus(:,:)
          real(r8), pointer :: avgbvs(:,:)

# ifdef SOLVE3D
#  if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
          real(r8), pointer :: avgPair(:,:)
#  endif
#  if defined BULK_FLUXES
          real(r8), pointer :: avgTair(:,:)
#  endif
#  if defined BULK_FLUXES || defined ECOSIM
          real(r8), pointer :: avgUwind(:,:)
          real(r8), pointer :: avgVwind(:,:)
#  endif
          real(r8), pointer :: avgstf(:,:)
          real(r8), pointer :: avgswf(:,:)
#  ifdef SHORTWAVE
          real(r8), pointer :: avgsrf(:,:)
#  endif
#  if defined BULK_FLUXES || defined FRC_COUPLING
          real(r8), pointer :: avglhf(:,:)
          real(r8), pointer :: avglrf(:,:)
          real(r8), pointer :: avgshf(:,:)
#  endif
#  if defined BULK_FLUXES && defined EMINUSP
          real(r8), pointer :: avgevap(:,:)
          real(r8), pointer :: avgrain(:,:)
#  endif
# endif

# ifdef NEARSHORE_MELLOR
!
!  Time-averaged radiation stresses.
!
          real(r8), pointer :: avgu2Sd(:,:)
          real(r8), pointer :: avgv2Sd(:,:)

          real(r8), pointer :: avgu2RS(:,:)
          real(r8), pointer :: avgv2RS(:,:)

          real(r8), pointer :: avgSxx2d(:,:)
          real(r8), pointer :: avgSxy2d(:,:)
          real(r8), pointer :: avgSyy2d(:,:)

#  ifdef SOLVE3D
          real(r8), pointer :: avgu3Sd(:,:,:)
          real(r8), pointer :: avgv3Sd(:,:,:)

          real(r8), pointer :: avgu3RS(:,:,:)
          real(r8), pointer :: avgv3RS(:,:,:)

          real(r8), pointer :: avgSxx3d(:,:,:)
          real(r8), pointer :: avgSxy3d(:,:,:)
          real(r8), pointer :: avgSyy3d(:,:,:)
#   ifdef NEARSHORE_MELLOR05
          real(r8), pointer :: avgSzx3d(:,:,:)
          real(r8), pointer :: avgSzy3d(:,:,:)
#   endif
#  endif
# endif
!
!  Time-averaged quadratic fields.
!
          real(r8), pointer :: avgZZ(:,:)
          real(r8), pointer :: avgU2(:,:)
          real(r8), pointer :: avgV2(:,:)
# ifdef SOLVE3D
          real(r8), pointer :: avgUU(:,:,:)
          real(r8), pointer :: avgUV(:,:,:)
          real(r8), pointer :: avgVV(:,:,:)

          real(r8), pointer :: avgHuon(:,:,:)
          real(r8), pointer :: avgHvom(:,:,:)

          real(r8), pointer :: avgTT(:,:,:,:)
          real(r8), pointer :: avgUT(:,:,:,:)
          real(r8), pointer :: avgVT(:,:,:,:)

          real(r8), pointer :: avgHuonT(:,:,:,:)
          real(r8), pointer :: avgHvomT(:,:,:,:)
# endif
!
!  Time-averaged vorticity fields.
!
          real(r8), pointer :: avgpvor2d(:,:)
          real(r8), pointer :: avgrvor2d(:,:)
# ifdef SOLVE3D
          real(r8), pointer :: avgpvor3d(:,:,:)
          real(r8), pointer :: avgrvor3d(:,:,:)
# endif

        END TYPE T_AVERAGE
!
        TYPE (T_AVERAGE), allocatable :: AVERAGE(:)
!
      CONTAINS
!
      SUBROUTINE allocate_average (ng, LBi, UBi, LBj, UBj)
!
!=======================================================================
!                                                                      !
!  This routine allocates all variables in the module for all nested   !
!  grids.                                                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, LBi, UBi, LBj, UBj
!
!  Local variable declarations.
!
      real(r8) :: size2d
!
!-----------------------------------------------------------------------
!  Allocate module variables.
!-----------------------------------------------------------------------
!
      IF (ng.eq.1 ) allocate ( AVERAGE(Ngrids) )
!
!  Set horizontal array size.
!
      size2d=REAL((UBi-LBi+1)*(UBj-LBj+1),r8)
!
!  Time-averaged state variables.
!
      IF (Aout(idFsur,ng)) THEN
        allocate ( AVERAGE(ng) % avgzeta(LBi:UBi,LBj:UBj) )
        Dmem(ng)=Dmem(ng)+size2d
      END IF

      IF (Aout(idUbar,ng)) THEN
        allocate ( AVERAGE(ng) % avgu2d(LBi:UBi,LBj:UBj) )
        Dmem(ng)=Dmem(ng)+size2d
      END IF
      IF (Aout(idVbar,ng)) THEN
        allocate ( AVERAGE(ng) % avgv2d(LBi:UBi,LBj:UBj) )
        Dmem(ng)=Dmem(ng)+size2d
      END IF

      IF (Aout(idu2dE,ng)) THEN
        allocate ( AVERAGE(ng) % avgu2dE(LBi:UBi,LBj:UBj) )
        Dmem(ng)=Dmem(ng)+size2d
      END IF
      IF (Aout(idv2dN,ng)) THEN
        allocate ( AVERAGE(ng) % avgv2dN(LBi:UBi,LBj:UBj) )
        Dmem(ng)=Dmem(ng)+size2d
      END IF

# ifdef SOLVE3D
      IF (Aout(idUvel,ng)) THEN
        allocate ( AVERAGE(ng) % avgu3d(LBi:UBi,LBj:UBj,N(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(N(ng),r8)*size2d
      END IF
      IF (Aout(idVvel,ng)) THEN
        allocate ( AVERAGE(ng) % avgv3d(LBi:UBi,LBj:UBj,N(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(N(ng),r8)*size2d
      END IF

      IF (Aout(idu3dE,ng)) THEN
        allocate ( AVERAGE(ng) % avgu3dE(LBi:UBi,LBj:UBj,N(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(N(ng),r8)*size2d
      END IF
      IF (Aout(idv3dN,ng)) THEN
        allocate ( AVERAGE(ng) % avgv3dN(LBi:UBi,LBj:UBj,N(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(N(ng),r8)*size2d
      END IF

      IF (Aout(idOvel,ng)) THEN
        allocate ( AVERAGE(ng) % avgw3d(LBi:UBi,LBj:UBj,0:N(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(N(ng)+1,r8)*size2d
      END IF
      IF (Aout(idWvel,ng)) THEN
        allocate ( AVERAGE(ng) % avgwvel(LBi:UBi,LBj:UBj,0:N(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(N(ng)+1,r8)*size2d
      END IF

      IF (Aout(idDano,ng)) THEN
        allocate ( AVERAGE(ng) % avgrho(LBi:UBi,LBj:UBj,N(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(N(ng),r8)*size2d
      END IF
      IF (ANY(Aout(idTvar(:),ng))) THEN
        allocate ( AVERAGE(ng) % avgt(LBi:UBi,LBj:UBj,N(ng),NT(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(N(ng)*NT(ng),r8)*size2d
      END IF

#  if defined LMD_MIXING || defined MY25_MIXING || defined GLS_MIXING
      IF (Aout(idVvis,ng)) THEN
        allocate ( AVERAGE(ng) % avgAKv(LBi:UBi,LBj:UBj,0:N(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(N(ng)+1,r8)*size2d
      END IF
      IF (Aout(idTdif,ng)) THEN
        allocate ( AVERAGE(ng) % avgAKt(LBi:UBi,LBj:UBj,0:N(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(N(ng)+1,r8)*size2d
      END IF
      IF (Aout(idSdif,ng)) THEN
        allocate ( AVERAGE(ng) % avgAKs(LBi:UBi,LBj:UBj,0:N(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(N(ng)+1,r8)*size2d
      END IF
#  endif
#  ifdef LMD_SKPP
      IF (Aout(idHsbl,ng)) THEN
        allocate ( AVERAGE(ng) % avghsbl(LBi:UBi,LBj:UBj) )
        Dmem(ng)=Dmem(ng)+size2d
      END IF
#  endif
#  ifdef LMD_BKPP
      IF (Aout(idHbbl,ng)) THEN
        allocate ( AVERAGE(ng) % avghbbl(LBi:UBi,LBj:UBj) )
        Dmem(ng)=Dmem(ng)+size2d
      END IF
#  endif
# endif

# if defined FORWARD_WRITE && defined SOLVE3D
!
!  Time-averaged 2D/3D coupling terms.
!
      IF (Aout(idUfx1,ng)) THEN
        allocate ( AVERAGE(ng) % avgDU_avg1(LBi:UBi,LBj:UBj) )
        Dmem(ng)=Dmem(ng)+size2d
      END IF
      IF (Aout(idUfx2,ng)) THEN
        allocate ( AVERAGE(ng) % avgDU_avg2(LBi:UBi,LBj:UBj) )
        Dmem(ng)=Dmem(ng)+size2d
      END IF

      IF (Aout(idVfx1,ng)) THEN
        allocate ( AVERAGE(ng) % avgDV_avg1(LBi:UBi,LBj:UBj) )
        Dmem(ng)=Dmem(ng)+size2d
      END IF
      IF (Aout(idVfx2,ng)) THEN
        allocate ( AVERAGE(ng) % avgDV_avg2(LBi:UBi,LBj:UBj) )
        Dmem(ng)=Dmem(ng)+size2d
      END IF
# endif
!
!  Time-averaged surface and bottom fluxes.
!
      IF (Aout(idUsms,ng)) THEN
        allocate ( AVERAGE(ng) % avgsus(LBi:UBi,LBj:UBj) )
        Dmem(ng)=Dmem(ng)+size2d
      END IF
      IF (Aout(idVsms,ng)) THEN
        allocate ( AVERAGE(ng) % avgsvs(LBi:UBi,LBj:UBj) )
        Dmem(ng)=Dmem(ng)+size2d
      END IF

      IF (Aout(idUbms,ng)) THEN
        allocate ( AVERAGE(ng) % avgbus(LBi:UBi,LBj:UBj) )
        Dmem(ng)=Dmem(ng)+size2d
      END IF
      IF (Aout(idVbms,ng)) THEN
        allocate ( AVERAGE(ng) % avgbvs(LBi:UBi,LBj:UBj) )
        Dmem(ng)=Dmem(ng)+size2d
      END IF

# ifdef SOLVE3D
#  if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
      IF (Aout(idPair,ng)) THEN
        allocate ( AVERAGE(ng) % avgPair(LBi:UBi,LBj:UBj) )
        Dmem(ng)=Dmem(ng)+size2d
      END IF
#  endif
#  if defined BULK_FLUXES
      IF (Aout(idTair,ng)) THEN
        allocate ( AVERAGE(ng) % avgTair(LBi:UBi,LBj:UBj) )
        Dmem(ng)=Dmem(ng)+size2d
      END IF
#  endif
#  if defined BULK_FLUXES || defined ECOSIM
      IF (Aout(idUair,ng)) THEN
        allocate ( AVERAGE(ng) % avgUwind(LBi:UBi,LBj:UBj) )
        Dmem(ng)=Dmem(ng)+size2d
      END IF
      IF (Aout(idVair,ng)) THEN
        allocate ( AVERAGE(ng) % avgVwind(LBi:UBi,LBj:UBj) )
        Dmem(ng)=Dmem(ng)+size2d
      END IF
#  endif
      IF (Aout(idTsur(itemp),ng)) THEN
        allocate ( AVERAGE(ng) % avgstf(LBi:UBi,LBj:UBj) )
        Dmem(ng)=Dmem(ng)+size2d
      END IF
#  ifdef SALINITY
      IF (Aout(idTsur(isalt),ng)) THEN
        allocate ( AVERAGE(ng) % avgswf(LBi:UBi,LBj:UBj) )
        Dmem(ng)=Dmem(ng)+size2d
      END IF
#  endif
#  ifdef SHORTWAVE
      IF (Aout(idSrad,ng)) THEN
        allocate ( AVERAGE(ng) % avgsrf(LBi:UBi,LBj:UBj) )
        Dmem(ng)=Dmem(ng)+size2d
      END IF
#  endif
#  if defined BULK_FLUXES || defined FRC_COUPLING
      IF (Aout(idLhea,ng)) THEN
        allocate ( AVERAGE(ng) % avglhf(LBi:UBi,LBj:UBj) )
        Dmem(ng)=Dmem(ng)+size2d
      END IF
      IF (Aout(idLrad,ng)) THEN
        allocate ( AVERAGE(ng) % avglrf(LBi:UBi,LBj:UBj) )
        Dmem(ng)=Dmem(ng)+size2d
      END IF
      IF (Aout(idShea,ng)) THEN
        allocate ( AVERAGE(ng) % avgshf(LBi:UBi,LBj:UBj) )
        Dmem(ng)=Dmem(ng)+size2d
      END IF
#  endif
#  if defined BULK_FLUXES && defined EMINUSP
      IF (Aout(idevap,ng)) THEN
        allocate ( AVERAGE(ng) % avgevap(LBi:UBi,LBj:UBj) )
        Dmem(ng)=Dmem(ng)+size2d
      END IF
      IF (Aout(idrain,ng)) THEN
        allocate ( AVERAGE(ng) % avgrain(LBi:UBi,LBj:UBj) )
        Dmem(ng)=Dmem(ng)+size2d
      END IF
#  endif
# endif

# ifdef NEARSHORE_MELLOR
!
!  Time-averaged radiation stresses.
!
      IF (Aout(idU2Sd,ng)) THEN
        allocate ( AVERAGE(ng) % avgu2Sd(LBi:UBi,LBj:UBj) )
        Dmem(ng)=Dmem(ng)+size2d
      END IF
      IF (Aout(idV2Sd,ng)) THEN
        allocate ( AVERAGE(ng) % avgv2Sd(LBi:UBi,LBj:UBj) )
        Dmem(ng)=Dmem(ng)+size2d
      END IF

      IF (Aout(idU2rs,ng)) THEN
        allocate ( AVERAGE(ng) % avgu2RS(LBi:UBi,LBj:UBj) )
        Dmem(ng)=Dmem(ng)+size2d
      END IF
      IF (Aout(idV2rs,ng)) THEN
        allocate ( AVERAGE(ng) % avgv2RS(LBi:UBi,LBj:UBj) )
        Dmem(ng)=Dmem(ng)+size2d
      END IF

      IF (Aout(idW2xx,ng)) THEN
        allocate ( AVERAGE(ng) % avgSxx2d(LBi:UBi,LBj:UBj) )
        Dmem(ng)=Dmem(ng)+size2d
      END IF
      IF (Aout(idW2xy,ng)) THEN
        allocate ( AVERAGE(ng) % avgSxy2d(LBi:UBi,LBj:UBj) )
        Dmem(ng)=Dmem(ng)+size2d
      END IF
      IF (Aout(idW2yy,ng)) THEN
        allocate ( AVERAGE(ng) % avgSyy2d(LBi:UBi,LBj:UBj) )
        Dmem(ng)=Dmem(ng)+size2d
      END IF

#  ifdef SOLVE3D
      IF (Aout(idU3Sd,ng)) THEN
        allocate ( AVERAGE(ng) % avgu3Sd(LBi:UBi,LBj:UBj,N(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(N(ng),r8)*size2d
      END IF
      IF (Aout(idV3Sd,ng)) THEN
        allocate ( AVERAGE(ng) % avgv3Sd(LBi:UBi,LBj:UBj,N(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(N(ng),r8)*size2d
      END IF

      IF (Aout(idU3rs,ng)) THEN
        allocate ( AVERAGE(ng) % avgu3RS(LBi:UBi,LBj:UBj,N(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(N(ng),r8)*size2d
      END IF
      IF (Aout(idV3rs,ng)) THEN
        allocate ( AVERAGE(ng) % avgv3RS(LBi:UBi,LBj:UBj,N(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(N(ng),r8)*size2d
      END IF

      IF (Aout(idW3xx,ng)) THEN
        allocate ( AVERAGE(ng) % avgSxx3d(LBi:UBi,LBj:UBj,N(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(N(ng),r8)*size2d
      END IF
      IF (Aout(idW3xy,ng)) THEN
        allocate ( AVERAGE(ng) % avgSxy3d(LBi:UBi,LBj:UBj,N(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(N(ng),r8)*size2d
      END IF
      IF (Aout(idW3yy,ng)) THEN
        allocate ( AVERAGE(ng) % avgSyy3d(LBi:UBi,LBj:UBj,N(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(N(ng),r8)*size2d
      END IF
#   ifdef NEARSHORE_MELLOR05
      IF (Aout(idW3zx,ng)) THEN
        allocate ( AVERAGE(ng) % avgSzx3d(LBi:UBi,LBj:UBj,N(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(N(ng),r8)*size2d
      END IF
      IF (Aout(idW3zy,ng)) THEN
        allocate ( AVERAGE(ng) % avgSzy3d(LBi:UBi,LBj:UBj,N(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(N(ng),r8)*size2d
      END IF
#   endif
#  endif
# endif
!
!  Time-averaged quadratic fields.
!
      IF (Aout(idZZav,ng)) THEN
        allocate ( AVERAGE(ng) % avgU2(LBi:UBi,LBj:UBj) )
        Dmem(ng)=Dmem(ng)+size2d
      END IF
      IF (Aout(idU2av,ng)) THEN
        allocate ( AVERAGE(ng) % avgV2(LBi:UBi,LBj:UBj) )
        Dmem(ng)=Dmem(ng)+size2d
      END IF
      IF (Aout(idV2av,ng)) THEN
        allocate ( AVERAGE(ng) % avgZZ(LBi:UBi,LBj:UBj) )
        Dmem(ng)=Dmem(ng)+size2d
      END IF

# ifdef SOLVE3D
      IF (Aout(idUUav,ng)) THEN
        allocate ( AVERAGE(ng) % avgUU(LBi:UBi,LBj:UBj,N(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(N(ng),r8)*size2d
      END IF
      IF (Aout(idVVav,ng)) THEN
        allocate ( AVERAGE(ng) % avgVV(LBi:UBi,LBj:UBj,N(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(N(ng),r8)*size2d
      END IF
      IF (Aout(idUVav,ng)) THEN
        allocate ( AVERAGE(ng) % avgUV(LBi:UBi,LBj:UBj,N(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(N(ng),r8)*size2d
      END IF

      IF (Aout(idHUav,ng)) THEN
        allocate ( AVERAGE(ng) % avgHuon(LBi:UBi,LBj:UBj,N(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(N(ng),r8)*size2d
      END IF
      IF (Aout(idHVav,ng)) THEN
        allocate ( AVERAGE(ng) % avgHvom(LBi:UBi,LBj:UBj,N(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(N(ng),r8)*size2d
      END IF

      IF (ANY(Aout(idTTav(:),ng))) THEN
        allocate ( AVERAGE(ng) % avgTT(LBi:UBi,LBj:UBj,N(ng),NT(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(N(ng)*NT(ng),r8)*size2d
      END IF
      IF (ANY(Aout(idUTav(:),ng))) THEN
        allocate ( AVERAGE(ng) % avgUT(LBi:UBi,LBj:UBj,N(ng),NT(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(N(ng)*NT(ng),r8)*size2d
      END IF
      IF (ANY(Aout(idVTav(:),ng))) THEN
        allocate ( AVERAGE(ng) % avgVT(LBi:UBi,LBj:UBj,N(ng),NT(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(N(ng)*NT(ng),r8)*size2d
      END IF

      IF (ANY(Aout(iHUTav(:),ng))) THEN
        allocate ( AVERAGE(ng)% avgHuonT(LBi:UBi,LBj:UBj,N(ng),NT(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(N(ng)*NT(ng),r8)*size2d
      END IF
      IF (ANY(Aout(iHVTav(:),ng))) THEN
        allocate ( AVERAGE(ng)% avgHvomT(LBi:UBi,LBj:UBj,N(ng),NT(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(N(ng)*NT(ng),r8)*size2d
      END IF
# endif
!
!  Time-averaged vorticity fields.
!
      IF (Aout(id2dPV,ng)) THEN
        allocate ( AVERAGE(ng) % avgpvor2d(LBi:UBi,LBj:UBj) )
        Dmem(ng)=Dmem(ng)+size2d
      END IF
      IF (Aout(id2dRV,ng)) THEN
        allocate ( AVERAGE(ng) % avgrvor2d(LBi:UBi,LBj:UBj) )
        Dmem(ng)=Dmem(ng)+size2d
      END IF
# ifdef SOLVE3D
      IF (Aout(id3dPV,ng)) THEN
        allocate ( AVERAGE(ng) % avgpvor3d(LBi:UBi,LBj:UBj,N(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(N(ng),r8)*size2d
      END IF
      IF (Aout(id3dRV,ng)) THEN
        allocate ( AVERAGE(ng) % avgrvor3d(LBi:UBi,LBj:UBj,N(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(N(ng),r8)*size2d
      END IF
# endif
!
      RETURN
      END SUBROUTINE allocate_average
!
      SUBROUTINE deallocate_average (ng)
!
!=======================================================================
!                                                                      !
!  This routine deallocates all variables in the module for all nested !
!  grids.                                                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars

# ifdef SUBOBJECT_DEALLOCATION
!
      USE destroy_mod, ONLY : destroy
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", deallocate_average"


# ifdef SUBOBJECT_DEALLOCATION
!
!-----------------------------------------------------------------------
!  Deallocate each variable in the derived-type T_AVERAGE structure
!  separately.
!-----------------------------------------------------------------------
!
!  Time-averaged state variables.
!
      IF (Aout(idFsur,ng)) THEN
        IF (.not.destroy(ng, AVERAGE(ng)%avgzeta, MyFile,               &
     &                   __LINE__, 'AVERAGE(ng)%avgzeta')) RETURN
      END IF

      IF (Aout(idUbar,ng)) THEN
        IF (.not.destroy(ng, AVERAGE(ng)%avgu2d, MyFile,                &
     &                   __LINE__, 'AVERAGE(ng)%avgu2d')) RETURN
      END IF

      IF (Aout(idVbar,ng)) THEN
        IF (.not.destroy(ng, AVERAGE(ng)%avgv2d, MyFile,                &
     &                   __LINE__, 'AVERAGE(ng)%avgv2d')) RETURN
      END IF

      IF (Aout(idu2dE,ng)) THEN
        IF (.not.destroy(ng, AVERAGE(ng)%avgu2dE, MyFile,               &
     &                   __LINE__, 'AVERAGE(ng)%avgu2dE')) RETURN
      END IF

      IF (Aout(idv2dN,ng)) THEN
        IF (.not.destroy(ng, AVERAGE(ng)%avgv2dN, MyFile,               &
     &                   __LINE__, 'AVERAGE(ng)%avgv2dN')) RETURN
      END IF

#  ifdef SOLVE3D
      IF (Aout(idUvel,ng)) THEN
        IF (.not.destroy(ng, AVERAGE(ng)%avgu3d, MyFile,                &
     &                   __LINE__, 'AVERAGE(ng)%avgu3d')) RETURN
      END IF

      IF (Aout(idVvel,ng)) THEN
        IF (.not.destroy(ng, AVERAGE(ng)%avgv3d, MyFile,                &
     &                   __LINE__, 'AVERAGE(ng)%avgv3d')) RETURN
      END IF

      IF (Aout(idu3dE,ng)) THEN
        IF (.not.destroy(ng, AVERAGE(ng)%avgu3dE, MyFile,               &
     &                   __LINE__, 'AVERAGE(ng)%avgu3dE')) RETURN
      END IF

      IF (Aout(idv3dN,ng)) THEN
        IF (.not.destroy(ng, AVERAGE(ng)%avgv3dN, MyFile,               &
     &                   __LINE__, 'AVERAGE(ng)%avgv3dN')) RETURN
      END IF

      IF (Aout(idOvel,ng)) THEN
        IF (.not.destroy(ng, AVERAGE(ng)%avgw3d, MyFile,                &
     &                   __LINE__, 'AVERAGE(ng)%avgw3d')) RETURN
      END IF

      IF (Aout(idWvel,ng)) THEN
        IF (.not.destroy(ng, AVERAGE(ng)%avgwvel, MyFile,               &
     &                   __LINE__, 'AVERAGE(ng)%avgwvel')) RETURN
      END IF

      IF (Aout(idDano,ng)) THEN
        IF (.not.destroy(ng, AVERAGE(ng)%avgrho, MyFile,                &
     &                   __LINE__, 'AVERAGE(ng)%avgrho')) RETURN
      END IF

      IF (ANY(Aout(idTvar(:),ng))) THEN
        IF (.not.destroy(ng, AVERAGE(ng)%avgt, MyFile,                  &
     &                   __LINE__, 'AVERAGE(ng)%avgt')) RETURN
      END IF

#   if defined LMD_MIXING || defined MY25_MIXING || defined GLS_MIXING
      IF (Aout(idVvis,ng)) THEN
        IF (.not.destroy(ng, AVERAGE(ng)%avgAKv, MyFile,                &
     &                   __LINE__, 'AVERAGE(ng)%avgAKv')) RETURN
      END IF

      IF (Aout(idTdif,ng)) THEN
        IF (.not.destroy(ng, AVERAGE(ng)%avgAKt, MyFile,                &
     &                   __LINE__, 'AVERAGE(ng)%avgAKt')) RETURN
      END IF

      IF (Aout(idSdif,ng)) THEN
        IF (.not.destroy(ng, AVERAGE(ng)%avgAKs, MyFile,                &
     &                   __LINE__, 'AVERAGE(ng)%avgAKs')) RETURN
      END IF
#   endif

#   ifdef LMD_SKPP
      IF (Aout(idHsbl,ng)) THEN
        IF (.not.destroy(ng, AVERAGE(ng)%avghsbl, MyFile,               &
     &                   __LINE__, 'AVERAGE(ng)%avghsbl')) RETURN
      END IF
#   endif

#   ifdef LMD_BKPP
      IF (Aout(idHbbl,ng)) THEN
        IF (.not.destroy(ng, AVERAGE(ng)%avghbbl, MyFile,               &
     &                   __LINE__, 'AVERAGE(ng)%avghbbl')) RETURN
      END IF
#   endif
#  endif

#  if defined FORWARD_WRITE && defined SOLVE3D
!
!  Time-averaged 2D/3D coupling terms.
!
      IF (Aout(idUfx1,ng)) THEN
        IF (.not.destroy(ng, AVERAGE(ng)%avgDU_avg1, MyFile,            &
     &                   __LINE__, 'AVERAGE(ng)%avgDU_avg1')) RETURN
      END IF

      IF (Aout(idUfx2,ng)) THEN
        IF (.not.destroy(ng, AVERAGE(ng)%avgDU_avg2, MyFile,            &
     &                   __LINE__, 'AVERAGE(ng)%avgDU_avg2')) RETURN
      END IF

      IF (Aout(idVfx1,ng)) THEN
        IF (.not.destroy(ng, AVERAGE(ng)%avgDV_avg1, MyFile,            &
     &                   __LINE__, 'AVERAGE(ng)%avgDV_avg1')) RETURN
      END IF

      IF (Aout(idVfx2,ng)) THEN
        IF (.not.destroy(ng, AVERAGE(ng)%avgDV_avg2, MyFile,            &
     &                   __LINE__, 'AVERAGE(ng)%avgDV_avg2')) RETURN
      END IF
#  endif
!
!  Time-averaged surface and bottom fluxes.
!
      IF (Aout(idUsms,ng)) THEN
        IF (.not.destroy(ng, AVERAGE(ng)%avgsus, MyFile,                &
     &                   __LINE__, 'AVERAGE(ng)%avgsus')) RETURN
      END IF

      IF (Aout(idVsms,ng)) THEN
        IF (.not.destroy(ng, AVERAGE(ng)%avgsvs, MyFile,                &
     &                   __LINE__, 'AVERAGE(ng)%avgsvs')) RETURN
      END IF

      IF (Aout(idUbms,ng)) THEN
        IF (.not.destroy(ng, AVERAGE(ng)%avgbus, MyFile,                &
     &                   __LINE__, 'AVERAGE(ng)%avgbus')) RETURN
      END IF

      IF (Aout(idVbms,ng)) THEN
        IF (.not.destroy(ng, AVERAGE(ng)%avgbvs, MyFile,                &
     &                   __LINE__, 'AVERAGE(ng)%avgbvs')) RETURN
      END IF

#  ifdef SOLVE3D

#   if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
      IF (Aout(idPair,ng)) THEN
        IF (.not.destroy(ng, AVERAGE(ng)%avgPair, MyFile,               &
     &                   __LINE__, 'AVERAGE(ng)%avgPair')) RETURN
      END IF
#   endif

#   if defined BULK_FLUXES
      IF (Aout(idTair,ng)) THEN
        IF (.not.destroy(ng, AVERAGE(ng)%avgTair, MyFile,               &
     &                   __LINE__, 'AVERAGE(ng)%avgTair')) RETURN
      END IF
#   endif

#   if defined BULK_FLUXES || defined ECOSIM
      IF (Aout(idUair,ng)) THEN
        IF (.not.destroy(ng, AVERAGE(ng)%avgUwind, MyFile,              &
     &                   __LINE__, 'AVERAGE(ng)%avgUwind')) RETURN
      END IF

      IF (Aout(idVair,ng)) THEN
        IF (.not.destroy(ng, AVERAGE(ng)%avgVwind, MyFile,              &
     &                   __LINE__, 'AVERAGE(ng)%avgVwind')) RETURN
      END IF
#   endif

      IF (Aout(idTsur(itemp),ng)) THEN
        IF (.not.destroy(ng, AVERAGE(ng)%avgstf, MyFile,                &
     &                   __LINE__, 'AVERAGE(ng)%avgstf')) RETURN
      END IF

#   ifdef SALINITY
      IF (Aout(idTsur(isalt),ng)) THEN
        IF (.not.destroy(ng, AVERAGE(ng)%avgswf, MyFile,                &
     &                   __LINE__, 'AVERAGE(ng)%avgswf')) RETURN
      END IF
#   endif

#   ifdef SHORTWAVE
      IF (Aout(idSrad,ng)) THEN
        IF (.not.destroy(ng, AVERAGE(ng)%avgsrf, MyFile,                &
     &                   __LINE__, 'AVERAGE(ng)%avgsrf')) RETURN
      END IF
#   endif

#   if defined BULK_FLUXES || defined FRC_COUPLING
      IF (Aout(idLhea,ng)) THEN
        IF (.not.destroy(ng, AVERAGE(ng)%avglhf, MyFile,                &
     &                   __LINE__, 'AVERAGE(ng)%avglhf')) RETURN
      END IF

      IF (Aout(idLrad,ng)) THEN
        IF (.not.destroy(ng, AVERAGE(ng)%avglrf, MyFile,                &
     &                   __LINE__, 'AVERAGE(ng)%avglrf')) RETURN
      END IF

      IF (Aout(idShea,ng)) THEN
        IF (.not.destroy(ng, AVERAGE(ng)%avgshf, MyFile,                &
     &                   __LINE__, 'AVERAGE(ng)%avgshf')) RETURN
      END IF
#   endif

#   if defined BULK_FLUXES && defined EMINUSP
      IF (Aout(idevap,ng)) THEN
        IF (.not.destroy(ng, AVERAGE(ng)%avgevap, MyFile,               &
     &                   __LINE__, 'AVERAGE(ng)%avgevap')) RETURN
      END IF

      IF (Aout(idrain,ng)) THEN
        IF (.not.destroy(ng, AVERAGE(ng)%avgrain, MyFile,               &
     &                   __LINE__, 'AVERAGE(ng)%avgrain')) RETURN
      END IF
#   endif
#  endif

#  ifdef NEARSHORE_MELLOR
!
!  Time-averaged radiation stresses.
!
      IF (Aout(idU2Sd,ng)) THEN
        IF (.not.destroy(ng, AVERAGE(ng)%avgu2Sd, MyFile,               &
     &                   __LINE__, 'AVERAGE(ng)%avgu2Sd')) RETURN
      END IF

      IF (Aout(idV2Sd,ng)) THEN
        IF (.not.destroy(ng, AVERAGE(ng)%avgv2Sd, MyFile,               &
     &                   __LINE__, 'AVERAGE(ng)%avgv2Sd')) RETURN
      END IF

      IF (Aout(idU2rs,ng)) THEN
        IF (.not.destroy(ng, AVERAGE(ng)%avgu2RS, MyFile,               &
     &                   __LINE__, 'AVERAGE(ng)%avgu2RS')) RETURN
      END IF

      IF (Aout(idV2rs,ng)) THEN
        IF (.not.destroy(ng, AVERAGE(ng)%avgv2RS, MyFile,               &
     &                   __LINE__, 'AVERAGE(ng)%avgv2RS')) RETURN
      END IF

      IF (Aout(idW2xx,ng)) THEN
        IF (.not.destroy(ng, AVERAGE(ng)%avgSxx2d, MyFile,              &
     &                   __LINE__, 'AVERAGE(ng)%avgSxx2d')) RETURN
      END IF

      IF (Aout(idW2xy,ng)) THEN
        IF (.not.destroy(ng, AVERAGE(ng)%avgSxy2d, MyFile,              &
     &                   __LINE__, 'AVERAGE(ng)%avgSxy2d')) RETURN
      END IF

      IF (Aout(idW2yy,ng)) THEN
        IF (.not.destroy(ng, AVERAGE(ng)%avgSyy2d, MyFile,              &
     &                   __LINE__, 'AVERAGE(ng)%avgSyy2d')) RETURN
      END IF

#   ifdef SOLVE3D
      IF (Aout(idU3Sd,ng)) THEN
        IF (.not.destroy(ng, AVERAGE(ng)%avgu3Sd, MyFile,               &
     &                   __LINE__, 'AVERAGE(ng)%avgu3Sd')) RETURN
      END IF

      IF (Aout(idV3Sd,ng)) THEN
        IF (.not.destroy(ng, AVERAGE(ng)%avgv3Sd, MyFile,               &
     &                   __LINE__, 'AVERAGE(ng)%avgv3Sd')) RETURN
      END IF

      IF (Aout(idU3rs,ng)) THEN
        IF (.not.destroy(ng, AVERAGE(ng)%avgu3RS, MyFile,               &
     &                   __LINE__, 'AVERAGE(ng)%avgu3RS')) RETURN
      END IF

      IF (Aout(idV3rs,ng)) THEN
        IF (.not.destroy(ng, AVERAGE(ng)%avgv3RS, MyFile,               &
     &                   __LINE__, 'AVERAGE(ng)%avgv3RS')) RETURN
      END IF

      IF (Aout(idW3xx,ng)) THEN
        IF (.not.destroy(ng, AVERAGE(ng)%avgSxx3d, MyFile,              &
     &                   __LINE__, 'AVERAGE(ng)%avgSxx3d')) RETURN
      END IF

      IF (Aout(idW3xy,ng)) THEN
        IF (.not.destroy(ng, AVERAGE(ng)%avgSxy3d, MyFile,              &
     &                   __LINE__, 'AVERAGE(ng)%avgSxy3d')) RETURN
      END IF

      IF (Aout(idW3yy,ng)) THEN
        IF (.not.destroy(ng, AVERAGE(ng)%avgSyy3d, MyFile,              &
     &                   __LINE__, 'AVERAGE(ng)%avgSyy3d')) RETURN
      END IF

#    ifdef NEARSHORE_MELLOR05
      IF (Aout(idW3zx,ng)) THEN
        IF (.not.destroy(ng, AVERAGE(ng)%avgSzx3d, MyFile,              &
     &                   __LINE__, 'AVERAGE(ng)%avgSzx3d')) RETURN
      END IF

      IF (Aout(idW3zy,ng)) THEN
        IF (.not.destroy(ng, AVERAGE(ng)%avgSzy3d, MyFile,              &
     &                   __LINE__, 'AVERAGE(ng)%avgSzy3d')) RETURN
      END IF
#    endif
#   endif
#  endif
!
!  Time-averaged quadratic fields.
!
      IF (Aout(idZZav,ng)) THEN
        IF (.not.destroy(ng, AVERAGE(ng)%avgU2, MyFile,                 &
     &                   __LINE__, 'AVERAGE(ng)%avgU2')) RETURN
      END IF

      IF (Aout(idU2av,ng)) THEN
        IF (.not.destroy(ng, AVERAGE(ng)%avgV2, MyFile,                 &
     &                   __LINE__, 'AVERAGE(ng)%avgV2')) RETURN
      END IF

      IF (Aout(idV2av,ng)) THEN
        IF (.not.destroy(ng, AVERAGE(ng)%avgZZ, MyFile,                 &
     &                   __LINE__, 'AVERAGE(ng)%avgZZ')) RETURN
      END IF

#  ifdef SOLVE3D
      IF (Aout(idUUav,ng)) THEN
        IF (.not.destroy(ng, AVERAGE(ng)%avgUU, MyFile,                 &
     &                   __LINE__, 'AVERAGE(ng)%avgUU')) RETURN
      END IF

      IF (Aout(idVVav,ng)) THEN
        IF (.not.destroy(ng, AVERAGE(ng)%avgVV, MyFile,                 &
     &                   __LINE__, 'AVERAGE(ng)%avgVV')) RETURN
      END IF

      IF (Aout(idUVav,ng)) THEN
        IF (.not.destroy(ng, AVERAGE(ng)%avgUV, MyFile,                 &
     &                   __LINE__, 'AVERAGE(ng)%avgUV')) RETURN
      END IF

      IF (Aout(idHUav,ng)) THEN
        IF (.not.destroy(ng, AVERAGE(ng)%avgHuon, MyFile,               &
     &                   __LINE__, 'AVERAGE(ng)%avgHuon')) RETURN
      END IF

      IF (Aout(idHVav,ng)) THEN
        IF (.not.destroy(ng, AVERAGE(ng)%avgHvom, MyFile,               &
     &                   __LINE__, 'AVERAGE(ng)%avgHvom')) RETURN
      END IF

      IF (ANY(Aout(idTTav(:),ng))) THEN
        IF (.not.destroy(ng, AVERAGE(ng)%avgTT, MyFile,                 &
     &                   __LINE__, 'AVERAGE(ng)%avgTT')) RETURN
      END IF

      IF (ANY(Aout(idUTav(:),ng))) THEN
        IF (.not.destroy(ng, AVERAGE(ng)%avgUT, MyFile,                 &
     &                   __LINE__, 'AVERAGE(ng)%avgUT')) RETURN
      END IF

      IF (ANY(Aout(idVTav(:),ng))) THEN
        IF (.not.destroy(ng, AVERAGE(ng)%avgVT, MyFile,                 &
     &                   __LINE__, 'AVERAGE(ng)%avgVT')) RETURN
      END IF

      IF (ANY(Aout(iHUTav(:),ng))) THEN
        IF (.not.destroy(ng, AVERAGE(ng)%avgHuonT, MyFile,              &
     &                   __LINE__, 'AVERAGE(ng)%avgHuonT')) RETURN
      END IF

      IF (ANY(Aout(iHVTav(:),ng))) THEN
        IF (.not.destroy(ng, AVERAGE(ng)%avgHvomT, MyFile,              &
     &                   __LINE__, 'AVERAGE(ng)%avgHvomT')) RETURN
      END IF
#  endif
!
!  Time-averaged vorticity fields.
!
      IF (Aout(id2dPV,ng)) THEN
        IF (.not.destroy(ng, AVERAGE(ng)%avgpvor2d, MyFile,             &
     &                   __LINE__, 'AVERAGE(ng)%avgpvor2d')) RETURN
      END IF

      IF (Aout(id2dRV,ng)) THEN
        IF (.not.destroy(ng, AVERAGE(ng)%avgrvor2d, MyFile,             &
     &                   __LINE__, 'AVERAGE(ng)%avgrvor2d')) RETURN
      END IF

#  ifdef SOLVE3D
      IF (Aout(id3dPV,ng)) THEN
        IF (.not.destroy(ng, AVERAGE(ng)%avgpvor3d, MyFile,             &
     &                   __LINE__, 'AVERAGE(ng)%avgpvor3d')) RETURN
      END IF

      IF (Aout(id3dRV,ng)) THEN
        IF (.not.destroy(ng, AVERAGE(ng)%avgrvor3d, MyFile,             &
     &                   __LINE__, 'AVERAGE(ng)%avgrvor3d')) RETURN
      END IF
#  endif
# endif
!
!-----------------------------------------------------------------------
!  Deallocate derived-type AVERAGE structure.
!-----------------------------------------------------------------------
!
      IF (ng.eq.Ngrids) THEN
        IF (allocated(AVERAGE)) deallocate ( AVERAGE )
      END IF
!
      RETURN
      END SUBROUTINE deallocate_average
!
      SUBROUTINE initialize_average (ng, tile)
!
!=======================================================================
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
      USE mod_ncparam
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, j
# ifdef SOLVE3D
      integer :: itrc, k
# endif

      real(r8), parameter :: IniVal = 0.0_r8

# include "set_bounds.h"
!
!  Set array initialization range.
!
# ifdef DISTRIBUTE
      Imin=BOUNDS(ng)%LBi(tile)
      Imax=BOUNDS(ng)%UBi(tile)
      Jmin=BOUNDS(ng)%LBj(tile)
      Jmax=BOUNDS(ng)%UBj(tile)
# else
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
        Imin=BOUNDS(ng)%LBi(tile)
      ELSE
        Imin=Istr
      END IF
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
        Imax=BOUNDS(ng)%UBi(tile)
      ELSE
        Imax=Iend
      END IF
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
        Jmin=BOUNDS(ng)%LBj(tile)
      ELSE
        Jmin=Jstr
      END IF
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
        Jmax=BOUNDS(ng)%UBj(tile)
      ELSE
        Jmax=Jend
      END IF
# endif
!
!-----------------------------------------------------------------------
!  Initialize module variables.
!-----------------------------------------------------------------------
!
!  Time-averaged state variables.
!
      IF (Aout(idFsur,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgzeta(i,j) = IniVal
          END DO
        END DO
      END IF

      IF (Aout(idUbar,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgu2d(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idVbar,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgv2d(i,j) = IniVal
          END DO
        END DO
      END IF

      IF (Aout(idu2dE,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgu2dE(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idv2dN,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgv2dN(i,j) = IniVal
          END DO
        END DO
      END IF

# ifdef SOLVE3D
      IF (Aout(idUvel,ng)) THEN
        DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              AVERAGE(ng) % avgu3d(i,j,k) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idVvel,ng)) THEN
        DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              AVERAGE(ng) % avgv3d(i,j,k) = IniVal
            END DO
          END DO
        END DO
      END IF

      IF (Aout(idu3dE,ng)) THEN
        DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              AVERAGE(ng) % avgu3dE(i,j,k) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idv3dN,ng)) THEN
        DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              AVERAGE(ng) % avgv3dN(i,j,k) = IniVal
            END DO
          END DO
        END DO
      END IF

      IF (Aout(idOvel,ng)) THEN
        DO k=0,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              AVERAGE(ng) % avgw3d(i,j,k) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idWvel,ng)) THEN
        DO k=0,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              AVERAGE(ng) % avgwvel(i,j,k) = IniVal
            END DO
          END DO
        END DO
      END IF

      IF (Aout(idDano,ng)) THEN
        DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              AVERAGE(ng) % avgrho(i,j,k) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (ANY(Aout(idTvar(:),ng))) THEN
        DO itrc=1,NT(ng)
          DO k=1,N(ng)
            DO j=Jmin,Jmax
              DO i=Imin,Imax
                AVERAGE(ng) % avgt(i,j,k,itrc) = IniVal
              END DO
            END DO
          END DO
        END DO
      END IF

#  if defined LMD_MIXING || defined MY25_MIXING || defined GLS_MIXING
      IF (Aout(idVvis,ng)) THEN
        DO k=0,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              AVERAGE(ng) % avgAKv(i,j,k) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idTdif,ng)) THEN
        DO k=0,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              AVERAGE(ng) % avgAKt(i,j,k) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idSdif,ng)) THEN
        DO k=0,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              AVERAGE(ng) % avgAKs(i,j,k) = IniVal
            END DO
          END DO
        END DO
      END IF
#  endif

#  ifdef LMD_SKPP
      IF (Aout(idHsbl,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avghsbl(i,j) = IniVal
          END DO
        END DO
      END IF
#  endif

#  ifdef LMD_BKPP
      IF (Aout(idHbbl,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avghbbl(i,j) = IniVal
          END DO
        END DO
      END IF
#  endif
# endif

# if defined FORWARD_WRITE && defined SOLVE3D
!
!  Time-averaged 2D/3D coupling terms.
!
      IF (Aout(idUfx1,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgDU_avg1(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idUfx2,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgDU_avg2(i,j) = IniVal
          END DO
        END DO
      END IF

      IF (Aout(idVfx1,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgDV_avg1(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idVfx2,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgDV_avg2(i,j) = IniVal
          END DO
        END DO
      END IF
# endif
!
!  Time-averaged surface and bottom fluxes.
!
      IF (Aout(idUsms,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgsus(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idVsms,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgsvs(i,j) = IniVal
          END DO
        END DO
      END IF

      IF (Aout(idUbms,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgbus(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idVbms,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgbvs(i,j) = IniVal
          END DO
        END DO
      END IF

# ifdef SOLVE3D
#  if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
      IF (Aout(idPair,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgPair(i,j) = IniVal
          END DO
        END DO
      END IF
#  endif
#  if defined BULK_FLUXES
      IF (Aout(idTair,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgTair(i,j) = IniVal
          END DO
        END DO
      END IF
#  endif
#  if defined BULK_FLUXES || defined ECOSIM
      IF (Aout(idUair,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgUwind(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idVair,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgVwind(i,j) = IniVal
          END DO
        END DO
      END IF
#  endif
      IF (Aout(idTsur(itemp),ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgstf(i,j) = IniVal
          END DO
        END DO
      END IF
#  ifdef SALINITY
      IF (Aout(idTsur(isalt),ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgswf(i,j) = IniVal
          END DO
        END DO
      END IF
#  endif
#  ifdef SHORTWAVE
      IF (Aout(idSrad,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgsrf(i,j) = IniVal
          END DO
        END DO
      END IF
#  endif

#  if defined BULK_FLUXES || defined FRC_COUPLING
      IF (Aout(idLhea,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avglhf(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idLrad,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avglrf(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idShea,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgshf(i,j) = IniVal
          END DO
        END DO
      END IF
#  endif

#  if defined BULK_FLUXES && defined EMINUSP
      IF (Aout(idevap,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgevap(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idrain,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgrain(i,j) = IniVal
          END DO
        END DO
      END IF
#  endif
# endif

# ifdef NEARSHORE_MELLOR
!
!  Time-averaged radiation stresses.
!
      IF (Aout(idU2Sd,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgu2Sd(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idV2Sd,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgv2Sd(i,j) = IniVal
          END DO
        END DO
      END IF

      IF (Aout(idU2rs,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgu2RS(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idV2rs,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgv2RS(i,j) = IniVal
          END DO
        END DO
      END IF

      IF (Aout(idW2xx,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgSxx2d(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idW2xy,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgSxy2d(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idW2yy,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgSyy2d(i,j) = IniVal
          END DO
        END DO
      END IF

#  ifdef SOLVE3D
      IF (Aout(idU3Sd,ng)) THEN
        DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              AVERAGE(ng) % avgu3Sd(i,j,k) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idV3Sd,ng)) THEN
        DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              AVERAGE(ng) % avgv3Sd(i,j,k) = IniVal
            END DO
          END DO
        END DO
      END IF

      IF (Aout(idU3rs,ng)) THEN
        DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              AVERAGE(ng) % avgu3RS(i,j,k) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idV3rs,ng)) THEN
        DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              AVERAGE(ng) % avgv3RS(i,j,k) = IniVal
            END DO
          END DO
        END DO
      END IF

      IF (Aout(idW3xx,ng)) THEN
        DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              AVERAGE(ng) % avgSxx3d(i,j,k) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idW3xy,ng)) THEN
        DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              AVERAGE(ng) % avgSxy3d(i,j,k) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idW3yy,ng)) THEN
        DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              AVERAGE(ng) % avgSyy3d(i,j,k) = IniVal
            END DO
          END DO
        END DO
      END IF
#   ifdef NEARSHORE_MELLOR05
      IF (Aout(idW3zx,ng)) THEN
        DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              AVERAGE(ng) % avgSzx3d(i,j,k) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idW3zy,ng)) THEN
        DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              AVERAGE(ng) % avgSzy3d(i,j,k) = IniVal
            END DO
          END DO
        END DO
      END IF
#   endif
#  endif
# endif
!
!  Time-averaged quadratic fields.
!
      IF (Aout(idZZav,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgU2(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idU2av,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgV2(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(idV2av,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgZZ(i,j) = IniVal
          END DO
        END DO
      END IF

# ifdef SOLVE3D
      IF (Aout(idUUav,ng)) THEN
        DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              AVERAGE(ng) % avgUU(i,j,k) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idVVav,ng)) THEN
        DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              AVERAGE(ng) % avgVV(i,j,k) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idUVav,ng)) THEN
        DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              AVERAGE(ng) % avgUV(i,j,k) = IniVal
            END DO
          END DO
        END DO
      END IF

      IF (Aout(idHUav,ng)) THEN
        DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              AVERAGE(ng) % avgHuon(i,j,k) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (Aout(idHVav,ng)) THEN
        DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              AVERAGE(ng) % avgHvom(i,j,k) = IniVal
            END DO
          END DO
        END DO
      END IF

      IF (ANY(Aout(idTTav(:),ng))) THEN
        DO itrc=1,NAT
          DO k=1,N(ng)
            DO j=Jmin,Jmax
              DO i=Imin,Imax
                AVERAGE(ng) % avgTT(i,j,k,itrc) = IniVal
              END DO
            END DO
          END DO
        END DO
      END IF
      IF (ANY(Aout(idUTav(:),ng))) THEN
        DO itrc=1,NAT
          DO k=1,N(ng)
            DO j=Jmin,Jmax
              DO i=Imin,Imax
                AVERAGE(ng) % avgUT(i,j,k,itrc) = IniVal
              END DO
            END DO
          END DO
        END DO
      END IF
      IF (ANY(Aout(idVTav(:),ng))) THEN
        DO itrc=1,NAT
          DO k=1,N(ng)
            DO j=Jmin,Jmax
              DO i=Imin,Imax
                AVERAGE(ng) % avgVT(i,j,k,itrc) = IniVal
              END DO
            END DO
          END DO
        END DO
      END IF

      IF (ANY(Aout(iHUTav(:),ng))) THEN
        DO itrc=1,NAT
          DO k=1,N(ng)
            DO j=Jmin,Jmax
              DO i=Imin,Imax
                AVERAGE(ng) % avgHuonT(i,j,k,itrc) = IniVal
              END DO
            END DO
          END DO
        END DO
      END IF
      IF (ANY(Aout(iHVTav(:),ng))) THEN
        DO itrc=1,NAT
          DO k=1,N(ng)
            DO j=Jmin,Jmax
              DO i=Imin,Imax
                AVERAGE(ng) % avgHvomT(i,j,k,itrc) = IniVal
              END DO
            END DO
          END DO
        END DO
      END IF
# endif
!
!  Time-averaged vorticity fields.
!
      IF (Aout(id2dPV,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgpvor2d(i,j) = IniVal
          END DO
        END DO
      END IF
      IF (Aout(id2dRV,ng)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            AVERAGE(ng) % avgrvor2d(i,j) = IniVal
          END DO
        END DO
      END IF
# ifdef SOLVE3D
      IF (Aout(id3dPV,ng)) THEN
        DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              AVERAGE(ng) % avgpvor3d(i,j,k) = IniVal
            END DO
          END DO
        END DO
      END IF
      IF (Aout(id3dRV,ng)) THEN
        DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              AVERAGE(ng) % avgrvor3d(i,j,k) = IniVal
            END DO
          END DO
        END DO
      END IF
# endif

      RETURN
      END SUBROUTINE initialize_average
#endif
      END MODULE mod_average
