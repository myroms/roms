#include "cppdefs.h"
      MODULE mod_forces
/*
************************************************************************
** Copyright (c) 2002 ROMS/TOMS Group                                 **
*************************************************** Hernan G. Arango ***
**                                                                    **
** Surface momentum stresses.                                         **
**                                                                    **
** sustr        Kinematic surface momentum flux (wind stress) in      **
**                the XI-direction (m2/s2) at horizontal U-points.    **
** sustrG       Latest two-time snapshots of input "sustr" grided     **
**                data used for interpolation.                        **
** svstr        Kinematic surface momentum flux (wind stress) in      **
**                the ETA-direction (m2/s2) at horizontal V-points.   **
** svstrG       Latest two-time snapshots of input "svstr" grided     **
**                data used for interpolation.                        **
**                                                                    **
** Bottom momentum stresses.                                          **
**                                                                    **
** bustr        Kinematic bottom momentum flux (bottom stress) in     **
**                the XI-direction (m2/s2) at horizontal U-points.    **
** bvstr        Kinematic bottom momentum flux (bottom stress) in     **
**                ETA-direction (m2/s2) at horizontal V-points.       **
**                                                                    **
** Surface wind induced waves.                                        **
**                                                                    **
** Awave        Surface wind induced wave amplitude (m).              **
** AwaveG       Latest two-time snapshots of input "Awave" grided     **
**                data used for interpolation.                        **
** Dwave        Surface wind induced wave direction (radians).        **
** DwaveG       Latest two-time snapshots of input "Dwave" grided     **
**                data used for interpolation.                        **
** Pwave        Surface wind induced wave period (s).                 **
** PwaveG       Latest two-time snapshots of input "Pwave" grided     **
**                data used for interpolation.                        **
** Ub_swan      Bottom orbital velocity read in from swan output      **
** Ub_swanG     Latest two-time snapshots of input "Ub_swan" grided   **
**                data used for interpolation.                        **
**                                                                    **
** Solar shortwave radiation flux.                                    **
**                                                                    **
** srflx        Kinematic surface shortwave solar radiation flux      **
**                (Celsius m/s) at horizontal RHO-points              **
** srflxG       Latest two-time snapshots of input "srflx" grided     **
**                data used for interpolation.                        **
**                                                                    **
** Cloud fraction.                                                    **
**                                                                    **
** cloud        Cloud fraction (percentage/100).                      **
** cloudG       Latest two-time snapshots of input "cloud" grided     **
**                data used for interpolation.                        **
**                                                                    **
** Surface heat fluxes, Atmosphere-Ocean bulk parameterization.       **
**                                                                    **
** lhflx        Kinematic net sensible heat flux (degC m/s).          **
** lrflx        Kinematic net longwave radiation (degC m/s).          **
** shflx        Kinematic net sensible heat flux (degC m/s).          **
**                                                                    **
** Surface air humidity.                                              **
**                                                                    **
** Hair         Surface air specific (g/kg) or relative humidity      **
**                (percentage).                                       **
** HairG        Latest two-time snapshots of input "Hair" grided      **
**                data used for interpolation.                        **
**                                                                    **
** Surface air pressure.                                              **
**                                                                    **
** Pair         Surface air pressure (mb).                            **
** PairG        Latest two-time snapshots of input "Pair" grided      **
**                data used for interpolation.                        **
**                                                                    **
** Surface air temperature.                                           **
**                                                                    **
** Tair         Surface air temperature (Celsius)                     **
** TairG        Latest two-time snapshots of input "Tair" grided      **
**                data used for interpolation.                        **
** Surface Winds.                                                     **
**                                                                    **
** Uwind        Surface wind in the XI-direction (m/s) at             **
**                horizontal RHO-points.                              **
** UwindG       Latest two-time snapshots of input "Uwind" grided     **
**                data used for interpolation.                        **
** Vwind        Surface wind in the ETA-direction (m/s) at            **
**                horizontal RHO-points.                              **
** VwindG       Latest two-time snapshots of input "Vwind" grided     **
**                data used for interpolation.                        **
**                                                                    **
** Rain fall rate.                                                    **
**                                                                    **
** rain         Rain fall rate (kg/m2/s).                             **
** rainG        Latest two-time snapshots of input "rain" grided      **
**                data used for interpolation.                        **
**                                                                    **
** Surface tracer fluxes.                                             **
**                                                                    **
** stflx        Kinematic surface flux of tracer type variables       **
**                (temperature: degC m/s; salinity: PSU m/s) at       **
**                horizontal RHO-points.                              **
** stflxG       Latest two-time snapshots of input "stflx" grided     **
**                data used for interpolation.                        **
**                                                                    **
** Bottom tracer fluxes.                                              **
**                                                                    **
** btflx        Kinematic bottom flux of tracer type variables        **
**                (temperature: degC m/s; salinity: PSU m/s) at       **
**                horizontal RHO-points.                              **
** btflxG       Latest two-time snapshots of input "btflx" grided     **
**                data used for interpolation.                        **
**                                                                    **
** Surface heat flux correction.                                      **
**                                                                    **
** dqdt         Kinematic surface net heat flux sensitivity to SST,   **
**                d(Q)/d(SST), (m/s).                                 **
** dqdtG        Latest two-time snapshots of input "dqdt" grided      **
**                data used for interpolation.                        **
** sst          Sea surface temperature (Celsius).                    **
** sstG         Latest two-time snapshots of input "sst" grided       **
**                data used for interpolation.                        **
**                                                                    **
** Surface freshwater flux correction.                                **
**                                                                    **
** sss          Sea surface salinity (PSU).                           **
** sssG         Latest two-time snapshots of input "sss" grided       **
**                data used for interpolation.                        **
**                                                                    **
** Surface spectral downwelling irradiance.                           **
**                                                                    **
** SpecIr       Spectral irradiance (NBands) from 400-700 nm at       **
**                5 nm bandwidth.                                     **
** avcos        Cosine of average zenith angle of downwelling         **
**                spectral photons.                                   **
**                                                                    **
************************************************************************
*/

        USE mod_kinds

        implicit none

        TYPE T_FORCES

          real(r8), pointer :: sustr(:,:)
          real(r8), pointer :: svstr(:,:)
#if !defined ANA_SMFLUX || !defined BULK_FLUXES
          real(r8), pointer :: sustrG(:,:,:)
          real(r8), pointer :: svstrG(:,:,:)
#endif
          real(r8), pointer :: bustr(:,:)
          real(r8), pointer :: bvstr(:,:)

#ifdef BBL_MODEL
# ifdef SWAN
          real(r8), pointer :: UB_swan(:,:)
# else
          real(r8), pointer :: Awave(:,:)
# endif
          real(r8), pointer :: Dwave(:,:)
          real(r8), pointer :: Pwave(:,:)
# ifndef ANA_WWAVE
#  ifdef SWAN
          real(r8), pointer :: UB_swanG(:,:,:)
#  else
          real(r8), pointer :: AwaveG(:,:,:)
#  endif
          real(r8), pointer :: DwaveG(:,:,:)
          real(r8), pointer :: PwaveG(:,:,:)
# endif
#endif

#ifdef SOLVE3D

# ifdef SHORTWAVE
          real(r8), pointer :: srflx(:,:)
#  ifndef ANA_SRFLUX
          real(r8), pointer :: srflxG(:,:,:)
#  endif
# endif

# ifdef CLOUDS
          real(r8), pointer :: cloud(:,:)
#  ifndef ANA_CLOUD
          real(r8), pointer :: cloudG(:,:,:)
#  endif
# endif
# ifdef BULK_FLUXES
          real(r8), pointer :: lhflx(:,:)
          real(r8), pointer :: lrflx(:,:)
#  ifndef LONGWAVE
          real(r8), pointer :: lrflxG(:,:,:)
#  endif
          real(r8), pointer :: shflx(:,:)
# endif

# if defined BULK_FLUXES || defined ECOSIM || \
    (defined SHORTWAVE && defined ANA_SRFLUX)
          real(r8), pointer :: Hair(:,:)
#  ifndef ANA_HUMIDITY
          real(r8), pointer :: HairG(:,:,:)
#  endif
          real(r8), pointer :: Tair(:,:)
#  ifndef ANA_TAIR
          real(r8), pointer :: TairG(:,:,:)
#  endif
# endif

# if defined BULK_FLUXES || defined ECOSIM
          real(r8), pointer :: Pair(:,:)
#  ifndef ANA_PAIR
          real(r8), pointer :: PairG(:,:,:)
#  endif
          real(r8), pointer :: Uwind(:,:)
          real(r8), pointer :: Vwind(:,:)
#  ifndef ANA_WINDS
          real(r8), pointer :: UwindG(:,:,:)
          real(r8), pointer :: VwindG(:,:,:)
#  endif
# endif

# ifdef BULK_FLUXES
          real(r8), pointer :: rain(:,:)
#  ifndef ANA_RAIN
          real(r8), pointer :: rainG(:,:,:)
#  endif
# endif

          real(r8), pointer :: stflx(:,:,:)
# if !defined ANA_STFLUX || !defined ANA_SSFLUX
          real(r8), pointer :: stflxG(:,:,:,:)
# endif
          real(r8), pointer :: btflx(:,:,:)
# if !defined ANA_BTFLUX || !defined ANA_BSFLUX
          real(r8), pointer :: btflxG(:,:,:,:)
# endif

# ifdef QCORRECTION
          real(r8), pointer :: dqdt(:,:)
          real(r8), pointer :: sst(:,:)
#  ifndef ANA_SST
          real(r8), pointer :: dqdtG(:,:,:)
          real(r8), pointer :: sstG(:,:,:)
#  endif
# endif

# if defined SALINITY && (defined SCORRECTION || defined SRELAXATION)
          real(r8), pointer :: sss(:,:)
#  ifndef ANA_SSS
          real(r8), pointer :: sssG(:,:,:)
#  endif
# endif

# ifdef ECOSIM
          real(r8), pointer :: SpecIr(:,:,:)
          real(r8), pointer :: avcos(:,:,:)
# endif
#endif

        END TYPE T_FORCES

        TYPE (T_FORCES), allocatable :: FORCES(:)

      CONTAINS

      SUBROUTINE allocate_forces (ng, LBi, UBi, LBj, UBj)
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
      IF (ng.eq.1) allocate ( FORCES(Ngrids) )
!
      allocate ( FORCES(ng) % sustr(LBi:UBi,LBj:UBj) )
      allocate ( FORCES(ng) % svstr(LBi:UBi,LBj:UBj) )

#if !defined ANA_SMFLUX || !defined BULK_FLUXES
      allocate ( FORCES(ng) % sustrG(LBi:UBi,LBj:UBj,2) )
      allocate ( FORCES(ng) % svstrG(LBi:UBi,LBj:UBj,2) )
#endif

      allocate ( FORCES(ng) % bustr(LBi:UBi,LBj:UBj) )
      allocate ( FORCES(ng) % bvstr(LBi:UBi,LBj:UBj) )

#ifdef BBL_MODEL
# ifdef SWAN
      allocate ( FORCES(ng) % UB_swan(LBi:UBi,LBj:UBj) )
# else
      allocate ( FORCES(ng) % Awave(LBi:UBi,LBj:UBj) )
# endif
      allocate ( FORCES(ng) % Dwave(LBi:UBi,LBj:UBj) )
      allocate ( FORCES(ng) % Pwave(LBi:UBi,LBj:UBj) )

# ifndef ANA_WWAVE
#  ifdef SWAN
      allocate ( FORCES(ng) % UB_swanG(LBi:UBi,LBj:UBj,2) )
#  else
      allocate ( FORCES(ng) % AwaveG(LBi:UBi,LBj:UBj,2) )
#  endif
      allocate ( FORCES(ng) % DwaveG(LBi:UBi,LBj:UBj,2) )
      allocate ( FORCES(ng) % PwaveG(LBi:UBi,LBj:UBj,2) )
# endif
#endif

#ifdef SOLVE3D

# ifdef SHORTWAVE
      allocate ( FORCES(ng) % srflx(LBi:UBi,LBj:UBj) )

#  ifndef ANA_SRFLUX
      allocate ( FORCES(ng) % srflxG(LBi:UBi,LBj:UBj,2) )
#  endif
# endif

# ifdef CLOUDS
      allocate ( FORCES(ng) % cloud(LBi:UBi,LBj:UBj) )

#  ifndef ANA_CLOUD
      allocate ( FORCES(ng) % cloudG(LBi:UBi,LBj:UBj,2) )
#  endif
# endif

# ifdef BULK_FLUXES
      allocate ( FORCES(ng) % lhflx(LBi:UBi,LBj:UBj) )
      allocate ( FORCES(ng) % lrflx(LBi:UBi,LBj:UBj) )

#  ifndef LONGWAVE
      allocate ( FORCES(ng) % lrflxG(LBi:UBi,LBj:UBj,2) )
#  endif

      allocate ( FORCES(ng) % shflx(LBi:UBi,LBj:UBj) )
# endif

# if defined BULK_FLUXES || defined ECOSIM || \
    (defined SHORTWAVE && defined ANA_SRFLUX)
      allocate ( FORCES(ng) % Hair(LBi:UBi,LBj:UBj) )

#  ifndef ANA_HUMIDITY
      allocate ( FORCES(ng) % HairG(LBi:UBi,LBj:UBj,2) )
#  endif

      allocate ( FORCES(ng) % Tair(LBi:UBi,LBj:UBj) )

#  ifndef ANA_TAIR
      allocate ( FORCES(ng) % TairG(LBi:UBi,LBj:UBj,2) )
#  endif
# endif

# if defined BULK_FLUXES || defined ECOSIM
      allocate ( FORCES(ng) % Pair(LBi:UBi,LBj:UBj) )

#  ifndef ANA_PAIR
      allocate ( FORCES(ng) % PairG(LBi:UBi,LBj:UBj,2) )
#  endif

      allocate ( FORCES(ng) % Uwind(LBi:UBi,LBj:UBj) )
      allocate ( FORCES(ng) % Vwind(LBi:UBi,LBj:UBj) )

#  ifndef ANA_WINDS
      allocate ( FORCES(ng) % UwindG(LBi:UBi,LBj:UBj,2) )
      allocate ( FORCES(ng) % VwindG(LBi:UBi,LBj:UBj,2) )
#  endif
# endif

# ifdef BULK_FLUXES
      allocate ( FORCES(ng) % rain(LBi:UBi,LBj:UBj) )

#  ifndef ANA_RAIN
      allocate ( FORCES(ng) % rainG(LBi:UBi,LBj:UBj,2) )
#  endif
# endif

      allocate ( FORCES(ng) % stflx(LBi:UBi,LBj:UBj,NT(ng)) )

# if !defined ANA_STFLUX || !defined ANA_SSFLUX
      allocate ( FORCES(ng) % stflxG(LBi:UBi,LBj:UBj,2,NT(ng)) )
# endif

      allocate ( FORCES(ng) % btflx(LBi:UBi,LBj:UBj,NT(ng)) )

# if !defined ANA_BTFLUX || !defined ANA_BSFLUX
      allocate ( FORCES(ng) % btflxG(LBi:UBi,LBj:UBj,2,NT(ng)) )
# endif

# ifdef QCORRECTION
      allocate ( FORCES(ng) % dqdt(LBi:UBi,LBj:UBj) )
      allocate ( FORCES(ng) % sst(LBi:UBi,LBj:UBj) )

#  ifndef ANA_SST
      allocate ( FORCES(ng) % dqdtG(LBi:UBi,LBj:UBj,2) )
      allocate ( FORCES(ng) % sstG(LBi:UBi,LBj:UBj,2) )
#  endif
# endif

# if defined SALINITY && (defined SCORRECTION || defined SRELAXATION)
      allocate ( FORCES(ng) % sss(LBi:UBi,LBj:UBj) )

#  ifndef ANA_SSS
      allocate ( FORCES(ng) % sssG(LBi:UBi,LBj:UBj,2) )
#  endif
# endif

# ifdef ECOSIM
      allocate ( FORCES(ng) % SpecIr(LBi:UBi,LBj:UBj,NBands) )
      allocate ( FORCES(ng) % avcos(LBi:UBi,LBj:UBj,NBands) )
# endif

#endif

      RETURN
      END SUBROUTINE allocate_forces

      SUBROUTINE initialize_forces (ng, tile)
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
#include "tile.h"
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
      FORCES(ng) % sustr(Imin:Imax,Jmin:Jmax) = IniVal
      FORCES(ng) % svstr(Imin:Imax,Jmin:Jmax) = IniVal

#if !defined ANA_SMFLUX || !defined BULK_FLUXES
      FORCES(ng) % sustrG(Imin:Imax,Jmin:Jmax,1:2) = IniVal
      FORCES(ng) % svstrG(Imin:Imax,Jmin:Jmax,1:2) = IniVal
#endif

      FORCES(ng) % bustr(Imin:Imax,Jmin:Jmax) = IniVal
      FORCES(ng) % bvstr(Imin:Imax,Jmin:Jmax) = IniVal

#ifdef BBL_MODEL
# ifdef SWAN
      FORCES(ng) % UB_swan(Imin:Imax,Jmin:Jmax) = IniVal
# else
      FORCES(ng) % Awave(Imin:Imax,Jmin:Jmax) = IniVal
# endif
      FORCES(ng) % Dwave(Imin:Imax,Jmin:Jmax) = IniVal
      FORCES(ng) % Pwave(Imin:Imax,Jmin:Jmax) = IniVal

# ifndef ANA_WWAVE
#  ifdef SWAN
      FORCES(ng) % Ub_swanG(Imin:Imax,Jmin:Jmax,1:2) = IniVal
#  else
      FORCES(ng) % AwaveG(Imin:Imax,Jmin:Jmax,1:2) = IniVal
#  endif
      FORCES(ng) % DwaveG(Imin:Imax,Jmin:Jmax,1:2) = IniVal
      FORCES(ng) % PwaveG(Imin:Imax,Jmin:Jmax,1:2) = IniVal
# endif
#endif

#ifdef SOLVE3D

# ifdef SHORTWAVE
      FORCES(ng) % srflx(Imin:Imax,Jmin:Jmax) = IniVal

#  ifndef ANA_SRFLUX
      FORCES(ng) % srflxG(Imin:Imax,Jmin:Jmax,1:2) = IniVal
#  endif
# endif

# ifdef CLOUDS
      FORCES(ng) % cloud(Imin:Imax,Jmin:Jmax) = IniVal

#  ifndef ANA_CLOUD
      FORCES(ng) % cloudG(Imin:Imax,Jmin:Jmax,1:2) = IniVal
#  endif
# endif

# ifdef BULK_FLUXES
      FORCES(ng) % lhflx(Imin:Imax,Jmin:Jmax) = IniVal
      FORCES(ng) % lrflx(Imin:Imax,Jmin:Jmax) = IniVal
      FORCES(ng) % shflx(Imin:Imax,Jmin:Jmax) = IniVal
# endif

# if defined BULK_FLUXES || defined ECOSIM || \
    (defined SHORTWAVE && defined ANA_SRFLUX)
      FORCES(ng) % Hair(Imin:Imax,Jmin:Jmax) = IniVal
      FORCES(ng) % Tair(Imin:Imax,Jmin:Jmax) = IniVal
# endif

# if defined BULK_FLUXES || defined ECOSIM
      FORCES(ng) % Pair(Imin:Imax,Jmin:Jmax) = IniVal

      FORCES(ng) % Uwind(Imin:Imax,Jmin:Jmax) = IniVal
      FORCES(ng) % Vwind(Imin:Imax,Jmin:Jmax) = IniVal
# endif

# ifdef BULK_FLUXES
      FORCES(ng) % rain(Imin:Imax,Jmin:Jmax) = IniVal
# endif

# if !defined LONGWAVE && defined BULK_FLUXES
      FORCES(ng) % lrflxG(Imin:Imax,Jmin:Jmax,1:2) = IniVal
# endif

# if defined BULK_FLUXES || defined ECOSIM
#  ifndef ANA_HUMIDITY
      FORCES(ng) % HairG(Imin:Imax,Jmin:Jmax,1:2) = IniVal
#  endif

#  ifndef ANA_PAIR
      FORCES(ng) % PairG(Imin:Imax,Jmin:Jmax,1:2) = IniVal
#  endif

#  ifndef ANA_TAIR
      FORCES(ng) % TairG(Imin:Imax,Jmin:Jmax,1:2) = IniVal
#  endif

#  ifndef ANA_WINDS
      FORCES(ng) % UwindG(Imin:Imax,Jmin:Jmax,1:2) = IniVal
      FORCES(ng) % VwindG(Imin:Imax,Jmin:Jmax,1:2) = IniVal
#  endif
# endif

# if !defined ANA_RAIN && defined BULK_FLUXES
      FORCES(ng) % rainG(Imin:Imax,Jmin:Jmax,1:2) = IniVal
# endif

      FORCES(ng) % stflx(Imin:Imax,Jmin:Jmax,1:NT(ng)) = IniVal
      FORCES(ng) % btflx(Imin:Imax,Jmin:Jmax,1:NT(ng)) = IniVal

# if !defined ANA_STFLUX || !defined ANA_SSFLUX
      FORCES(ng) % stflxG(Imin:Imax,Jmin:Jmax,1:2,1:NT(ng)) = IniVal
# endif

# if !defined ANA_BTFLUX || !defined ANA_BSFLUX
      FORCES(ng) % btflxG(Imin:Imax,Jmin:Jmax,1:2,1:NT(ng)) = IniVal
# endif

# ifdef QCORRECTION
      FORCES(ng) % dqdt(Imin:Imax,Jmin:Jmax) = IniVal
      FORCES(ng) % sst(Imin:Imax,Jmin:Jmax) = IniVal

#  ifndef ANA_SST
      FORCES(ng) % dqdtG(Imin:Imax,Jmin:Jmax,1:2) = IniVal
      FORCES(ng) % sstG(Imin:Imax,Jmin:Jmax,1:2) = IniVal
#  endif
# endif

# if defined SALINITY && (defined SCORRECTION || defined SRELAXATION)
      FORCES(ng) % sss(Imin:Imax,Jmin:Jmax) = IniVal

#  ifndef ANA_SSS
      FORCES(ng) % sssG(Imin:Imax,Jmin:Jmax,1:2) = IniVal
#  endif
# endif

# ifdef ECOSIM
      FORCES(ng) % SpecIr(Imin:Imax,Jmin:Jmax,1:NBands) = IniVal
      FORCES(ng) % avcos(Imin:Imax,Jmin:Jmax,1:NBands) = IniVal
# endif

#endif

      RETURN
      END SUBROUTINE initialize_forces

      END MODULE mod_forces
