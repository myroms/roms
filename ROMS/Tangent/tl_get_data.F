#include "cppdefs.h"
# ifdef TANGENT
      SUBROUTINE tl_get_data (ng)
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2011 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine reads in forcing, climatology and other data from      !
!  NetCDF files.  If there is more than one time-record,  data is      !
!  loaded into global  two-time  record arrays. The interpolation      !
!  is carried elsewhere.                                               !
!                                                                      !
!  Currently, this routine is only executed in serial mode by the      !
!  main thread.                                                        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_boundary
# ifdef FORWARD_READ
      USE mod_coupling
# endif
# ifdef CLIMATOLOGY
      USE mod_clima
# endif
      USE mod_forces
      USE mod_grid
      USE mod_iounits
      USE mod_mixing
      USE mod_ncparam
# ifdef FORWARD_READ
      USE mod_ocean
# endif
      USE mod_scalars
# if defined UV_PSOURCE || defined TS_PSOURCE || defined Q_PSOURCE
      USE mod_sources
# endif
      USE mod_stepping
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      logical, dimension(3) :: update =                                 &
     &         (/ .FALSE., .FALSE., .FALSE. /)

      integer :: ILB, IUB, JLB, JUB
      integer :: LBi, UBi, LBj, UBj
      integer :: i, my_tile

# ifdef FORWARD_MIXING
      real(r8) :: scale
# endif
!
!  Lower and upper bounds for nontiled (global values) boundary arrays.
!
      my_tile=-1                           ! for global values
      ILB=BOUNDS(ng)%LBi(my_tile)
      IUB=BOUNDS(ng)%UBi(my_tile)
      JLB=BOUNDS(ng)%LBj(my_tile)
      JUB=BOUNDS(ng)%UBj(my_tile)
!
!  Lower and upper bounds for tiled arrays.
!
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on input data time wall clock.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, iTLM, 3)
# endif
# ifdef FRC_FILE
!
!=======================================================================
!  Read in forcing data from FORCING NetCDF file.
!=======================================================================
# endif

# if !defined ANA_PSOURCE && (defined UV_PSOURCE || \
                              defined TS_PSOURCE || defined Q_PSOURCE)
!
!-----------------------------------------------------------------------
!  Point Sources/Sinks time dependent data.
!-----------------------------------------------------------------------
!
#  if defined UV_PSOURCE || defined Q_PSOURCE
!
!  Point Source/Sink vertically integrated mass transport.
!
      CALL get_ngfld (ng, iTLM, idRtra, ncFRCid(idRtra,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                1, Nsrc(ng), 1, 2, 1, Nsrc(ng), 1,                &
     &                SOURCES(ng) % QbarG(1,1))
      IF (exit_flag.ne.NoError) RETURN
#  endif

#  if defined TS_PSOURCE && defined SOLVE3D
!
!  Tracer Sources/Sinks.
!
      DO i=1,NT(ng)
        IF (LtracerSrc(i,ng)) THEN
          CALL get_ngfld (ng, iTLM, idRtrc(i), ncFRCid(idRtrc(i),ng),   &
     &                    nFfiles(ng), FRC(1,ng), update(1),            &
     &                    1, Nsrc(ng), N(ng), 2, 1, Nsrc(ng), N(ng),    &
     &                    SOURCES(ng) % TsrcG(1,1,1,i))
          IF (exit_flag.ne.NoError) RETURN
        END IF
      END DO
#  endif
# endif

# if !defined ANA_WINDS    && \
     ((defined BULK_FLUXES && !defined NL_BULK_FLUXES) || \
      defined ECOSIM)
!
!-----------------------------------------------------------------------
!  Surface wind components.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iTLM, idUair, ncFRCid(idUair,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#  endif
     &                FORCES(ng) % UwindG(LBi,LBj,1))
      IF (exit_flag.ne.NoError) RETURN

      CALL get_2dfld (ng , iTLM, idVair, ncFRCid(idVair,ng),            &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#  endif
     &                FORCES(ng) % VwindG(LBi,LBj,1))
      IF (exit_flag.ne.NoError) RETURN
# endif

# ifndef AIR_OCEAN
#  if !defined ANA_SMFLUX  && \
      !defined BULK_FLUXES && !defined NL_BULK_FLUXES
!
!-----------------------------------------------------------------------
!  Surface wind stress components.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iTLM, idUsms, ncFRCid(idUsms,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % umask(LBi,LBj),                        &
#   endif
     &                FORCES(ng) % sustrG(LBi,LBj,1))
      IF (exit_flag.ne.NoError) RETURN

      CALL get_2dfld (ng, iTLM, idVsms, ncFRCid(idVsms,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % vmask(LBi,LBj),                        &
#   endif
     &                FORCES(ng) % svstrG(LBi,LBj,1))
      IF (exit_flag.ne.NoError) RETURN
#  endif

#  ifdef NL_BULK_FLUXES
!
!-----------------------------------------------------------------------
!  Surface wind stress components from NLM bulk flux computation.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iTLM, idUsms, BLK(ng)%ncid,                   &
     &                1, BLK(ng), update(1),                            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % umask(LBi,LBj),                        &
#   endif
     &                FORCES(ng) % sustrG(LBi,LBj,1))
      IF (exit_flag.ne.NoError) RETURN

      CALL get_2dfld (ng, iTLM, idVsms, BLK(ng)%ncid,                   &
     &                1, BLK(ng), update(1),                            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % vmask(LBi,LBj),                        &
#   endif
     &                FORCES(ng) % svstrG(LBi,LBj,1))
      IF (exit_flag.ne.NoError) RETURN
#  endif
# endif

# if !defined ANA_PAIR && \
     ((defined BULK_FLUXES && !defined NL_BULK_FLUXES) || \
      defined ECOSIM       || defined ATM_PRESS)
!
!-----------------------------------------------------------------------
!  Surface air pressure.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iTLM, idPair, ncFRCid(idPair,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#  endif
     &                FORCES(ng) % PairG(LBi,LBj,1))
      IF (exit_flag.ne.NoError) RETURN
# endif

# if !defined ANA_WWAVE && defined WAVE_DATA
!
!-----------------------------------------------------------------------
!  Surface wind induced wave amplitude, direction and period.
!-----------------------------------------------------------------------
!
#  ifdef WAVES_DIR
      CALL get_2dfld (ng, iTLM, idWdir, ncFRCid(idWdir,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#   endif
     &                FORCES(ng) % DwaveG(LBi,LBj,1))
      IF (exit_flag.ne.NoError) RETURN
#  endif

#  ifdef WAVES_HEIGHT
      CALL get_2dfld (ng, iTLM, idWamp, ncFRCid(idWamp,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#   endif
     &                FORCES(ng) % HwaveG(LBi,LBj,1))
      IF (exit_flag.ne.NoError) RETURN
#  endif

#  ifdef WAVES_LENGTH
      CALL get_2dfld (ng, iTLM, idWlen, ncFRCid(idWlen,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#   endif
     &                FORCES(ng) % LwaveG(LBi,LBj,1))
      IF (exit_flag.ne.NoError) RETURN
#  endif

#  ifdef WAVES_TOP_PERIOD
      CALL get_2dfld (ng, iTLM, idWptp, ncFRCid(idWptp,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#   endif
     &                FORCES(ng) % Pwave_topG(LBi,LBj,1))
      IF (exit_flag.ne.NoError) RETURN
#  endif

#  ifdef WAVES_BOT_PERIOD
      CALL get_2dfld (ng, iTLM, idWpbt, ncFRCid(idWpbt,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#   endif
     &                FORCES(ng) % Pwave_botG(LBi,LBj,1))
      IF (exit_flag.ne.NoError) RETURN
#  endif

#  if defined WAVES_UB
      CALL get_2dfld (ng, iTLM, idWorb, ncFRCid(idWorb,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#   endif
     &                FORCES(ng) % Ub_swanG(LBi,LBj,1))
      IF (exit_flag.ne.NoError) RETURN
#  endif

#  if defined TKE_WAVEDISS
      CALL get_2dfld (ng, iTLM, idWdis, ncFRCid(idWdis,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#   endif
     &                FORCES(ng) % Wave_dissipG(LBi,LBj,1))
      IF (exit_flag.ne.NoError) RETURN
#  endif

#  if defined SVENDSEN_ROLLER
      CALL get_2dfld (ng, iTLM, idWbrk, ncFRCid(idWbrk,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#   endif
     &                FORCES(ng) % Wave_breakG(LBi,LBj,1))
      IF (exit_flag.ne.NoError) RETURN
#  endif
# endif

# ifdef SOLVE3D

#  if !defined ANA_CLOUD && defined CLOUDS
!
!-----------------------------------------------------------------------
!  Cloud fraction.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iTLM, idCfra, ncFRCid(idCfra,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#   endif
     &                FORCES(ng) % cloudG(LBi,LBj,1))
      IF (exit_flag.ne.NoError) RETURN
#  endif

#  if !defined ANA_SRFLUX && defined SHORTWAVE
!
!-----------------------------------------------------------------------
!  Surface solar shortwave radiation.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iTLM, idSrad, ncFRCid(idSrad,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#   endif
     &                FORCES(ng) % srflxG(LBi,LBj,1))
      IF (exit_flag.ne.NoError) RETURN
#  endif

#  if (defined BULK_FLUXES && !defined NL_BULK_FLUXES) && \
      !defined LONGWAVE    && !defined LONGWAVE_OUT
!
!-----------------------------------------------------------------------
!  Surface net longwave radiation.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iTLM, idLrad, ncFRCid(idLrad,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#   endif
     &                FORCES(ng) % lrflxG(LBi,LBj,1))
      IF (exit_flag.ne.NoError) RETURN
#  endif

#  if (defined BULK_FLUXES && !defined NL_BULK_FLUXES) && \
      defined LONGWAVE_OUT
!
!-----------------------------------------------------------------------
!  Surface downwelling longwave radiation.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iTLM, idLdwn, ncFRCid(idLdwn,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#   endif
     &                FORCES(ng) % lrflxG(LBi,LBj,1))
      IF (exit_flag.ne.NoError) RETURN
#  endif

#  if !defined ANA_TAIR && \
     ((defined BULK_FLUXES && !defined NL_BULK_FLUXES) || \
      defined ECOSIM       || \
      (defined SHORTWAVE   && defined ANA_SRFLUX && defined ALBEDO))
!
!-----------------------------------------------------------------------
!  Surface air temperature.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iTLM, idTair, ncFRCid(idTair,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#   endif
     &                FORCES(ng) % TairG(LBi,LBj,1))
      IF (exit_flag.ne.NoError) RETURN
#  endif

#  if !defined ANA_HUMIDITY && \
      ((defined BULK_FLUXES && !defined NL_BULK_FLUXES) || \
       defined ECOSIM)
!
!-----------------------------------------------------------------------
!  Surface air humidity.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iTLM, idQair, ncFRCid(idQair,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#   endif
     &                FORCES(ng) % HairG(LBi,LBj,1))
      IF (exit_flag.ne.NoError) RETURN
#  endif

#  if !defined ANA_RAIN    && \
      (defined BULK_FLUXES && !defined NL_BULK_FLUXES)
!
!-----------------------------------------------------------------------
!  Rain fall rate.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iTLM, idrain, ncFRCid(idrain,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#   endif
     &                FORCES(ng) % rainG(LBi,LBj,1))
      IF (exit_flag.ne.NoError) RETURN
#  endif

#  if !defined ANA_STFLUX  && \
      !defined BULK_FLUXES && !defined NL_BULK_FLUXES
!
!-----------------------------------------------------------------------
!  Surface net heat flux.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iTLM, idTsur(itemp),                          &
     &                ncFRCid(idTsur(itemp),ng),                        &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#   endif
     &                FORCES(ng) % stflxG(LBi,LBj,1,itemp))
      IF (exit_flag.ne.NoError) RETURN
#  endif

#  ifdef NL_BULK_FLUXES
!
!-----------------------------------------------------------------------
!  Surface net heat flux from NLM bulk flux computation.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iTLM, idTsur(itemp), BLK(ng)%ncid,            &
     &                1, BLK(ng), update(1),                            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#   endif
     &                FORCES(ng) % stflxG(LBi,LBj,1,itemp))
      IF (exit_flag.ne.NoError) RETURN
#  endif

#  if !defined ANA_SST && defined QCORRECTION
!
!-----------------------------------------------------------------------
!  Surface net heat flux correction fields: sea surface temperature
!  (SST) and heat flux sensitivity to SST (dQdSST).
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iTLM, idSSTc, ncFRCid(idSSTc,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#   endif
     &                FORCES(ng) % sstG(LBi,LBj,1))
      IF (exit_flag.ne.NoError) RETURN

      CALL get_2dfld (ng, iTLM, iddQdT, ncFRCid(iddQdT,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#   endif
     &                FORCES(ng) % dqdtG(LBi,LBj,1))
      IF (exit_flag.ne.NoError) RETURN
#  endif

#  ifndef ANA_BTFLUX
!
!-----------------------------------------------------------------------
!  Bottom net heat flux.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iTLM, idTbot(itemp),                          &
     &                ncFRCid(idTbot(itemp),ng),                        &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#   endif
     &                FORCES(ng) % btflxG(LBi,LBj,1,itemp))
      IF (exit_flag.ne.NoError) RETURN
#  endif

#  ifdef SALINITY
#   if !(defined ANA_SSFLUX || defined EMINUSP || defined SRELAXATION)
!
!-----------------------------------------------------------------------
!  Surface net freshwater flux: E-P.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iTLM, idsfwf, ncFRCid(idsfwf,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#    ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#    endif
     &                FORCES(ng) % stflxG(LBi,LBj,1,isalt))
      IF (exit_flag.ne.NoError) RETURN
#   endif

#   if defined EMINUSP && defined NL_BULK_FLUXES
!
!-----------------------------------------------------------------------
!  Surface net freshwater flux (E-P) from NLM bulk flux computation.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iTLM, idEmPf, BLK(ng)%ncid,                   &
     &                1, BLK(ng), update(1),                            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#    ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#    endif
     &                FORCES(ng) % stflxG(LBi,LBj,1,isalt))
      IF (exit_flag.ne.NoError) RETURN
#   endif

#   if !defined ANA_SSS && (defined SCORRECTION || defined SRELAXATION)
!
!-----------------------------------------------------------------------
!  Surface net freshwater flux correction field: sea surface salinity.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iTLM, idSSSc, ncFRCid(idSSSc,ng),             &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#    ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#    endif
     &                FORCES(ng) % sssG(LBi,LBj,1))
      IF (exit_flag.ne.NoError) RETURN
#   endif

#   ifndef ANA_BSFLUX
!
!-----------------------------------------------------------------------
!  Bottom net freshwater flux.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iTLM, idTbot(isalt),                          &
     &                ncFRCid(idTbot(isalt),ng),                        &
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#    ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#    endif
     &                FORCES(ng) % btflxG(LBi,LBj,1,isalt))
      IF (exit_flag.ne.NoError) RETURN
#   endif
#  endif

#  if defined BIOLOGY || defined SEDIMENT || defined T_PASSIVE
#   ifndef ANA_SPFLUX
!
!-----------------------------------------------------------------------
!  Passive tracers surface fluxes.
!-----------------------------------------------------------------------
!
      DO i=NAT+1,NT(ng)
        CALL get_2dfld (ng, iTLM, idTsur(i), ncFRCid(idTsur(i),ng),     &
     &                  nFfiles(ng), FRC(1,ng), update(1),              &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
#    ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
#    endif
     &                  FORCES(ng) % stflxG(LBi,LBj,1,i))
        IF (exit_flag.ne.NoError) RETURN
      END DO
#   endif

#   ifndef ANA_BPFLUX
!
!-----------------------------------------------------------------------
!  Passive tracers bottom fluxes.
!-----------------------------------------------------------------------
!
      DO i=NAT+1,NT(ng)
        CALL get_2dfld (ng, iTLM, idTbot(i), ncFRCid(idTbot(i),ng),     &
     &                  nFfiles(ng), FRC(1,ng), update(1),              &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
#    ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
#    endif
     &                  FORCES(ng) % btflxG(LBi,LBj,1,i))
        IF (exit_flag.ne.NoError) RETURN
      END DO
#   endif
#  endif
# endif
!
!=======================================================================
!  Read in open boundary conditions from BOUNDARY NetCDF file.
!=======================================================================

# ifndef ANA_FSOBC
!
      IF (tl_LBC(iwest,isFsur,ng)%acquire) THEN
        CALL get_ngfld (ng, iTLM, idZbry(iwest), BRY(ng)%ncid,          &
     &                  1, BRY(ng), update(1),                          &
     &                  JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,                 &
     &                  BOUNDARY(ng) % zetaG_west(JLB,1))
        IF (exit_flag.ne.NoError) RETURN
      END IF

      IF (tl_LBC(ieast,isFsur,ng)%acquire) THEN
        CALL get_ngfld (ng, iTLM, idZbry(ieast), BRY(ng)%ncid,          &
     &                  1, BRY(ng), update(1),                          &
     &                  JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,                 &
     &                  BOUNDARY(ng) % zetaG_east(JLB,1))
        IF (exit_flag.ne.NoError) RETURN
      END IF

      IF (tl_LBC(isouth,isFsur,ng)%acquire) THEN
        CALL get_ngfld (ng, iTLM, idZbry(isouth), BRY(ng)%ncid,         &
     &                  1, BRY(ng), update(1),                          &
     &                  ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,                 &
     &                  BOUNDARY(ng) % zetaG_south(ILB,1))
        IF (exit_flag.ne.NoError) RETURN
      END IF

      IF (tl_LBC(inorth,isFsur,ng)%acquire) THEN
        CALL get_ngfld (ng, iTLM, idZbry(inorth), BRY(ng)%ncid,         &
     &                  1, BRY(ng), update(1),                          &
     &                  ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,                 &
     &                  BOUNDARY(ng) % zetaG_north(ILB,1))
        IF (exit_flag.ne.NoError) RETURN
      END IF
# endif

# ifndef ANA_M2OBC
!
      IF (tl_LBC(iwest,isUbar,ng)%acquire) THEN
        CALL get_ngfld (ng, iTLM, idU2bc(iwest), BRY(ng)%ncid,          &
     &                  1, BRY(ng), update(1),                          &
     &                  JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,                 &
     &                  BOUNDARY(ng) % ubarG_west(JLB,1))
        IF (exit_flag.ne.NoError) RETURN
      END IF

      IF (tl_LBC(iwest,isVbar,ng)%acquire) THEN
        CALL get_ngfld (ng, iTLM, idV2bc(iwest), BRY(ng)%ncid,          &
     &                  1, BRY(ng), update(1),                          &
     &                  JLB, JUB, 1, 2, 1, Mm(ng)+1, 1,                 &
     &                  BOUNDARY(ng) % vbarG_west(JLB,1))
        IF (exit_flag.ne.NoError) RETURN
      END IF

      IF (tl_LBC(ieast,isUbar,ng)%acquire) THEN
        CALL get_ngfld (ng, iTLM, idU2bc(ieast), BRY(ng)%ncid,          &
     &                  1, BRY(ng), update(1),                          &
     &                  JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,                 &
     &                  BOUNDARY(ng) % ubarG_east(JLB,1))
        IF (exit_flag.ne.NoError) RETURN
      END IF

      IF (tl_LBC(ieast,isVbar,ng)%acquire) THEN
        CALL get_ngfld (ng, iTLM, idV2bc(ieast), BRY(ng)%ncid,          &
     &                  1, BRY(ng), update(1),                          &
     &                  JLB, JUB, 1, 2, 1, Mm(ng)+1, 1,                 &
     &                  BOUNDARY(ng) % vbarG_east(JLB,1))
        IF (exit_flag.ne.NoError) RETURN
      END IF

      IF (tl_LBC(isouth,isUbar,ng)%acquire) THEN
        CALL get_ngfld (ng, iTLM, idU2bc(isouth), BRY(ng)%ncid,         &
     &                  1, BRY(ng), update(1),                          &
     &                  ILB, IUB, 1, 2, 1, Lm(ng)+1, 1,                 &
     &                  BOUNDARY(ng) % ubarG_south(ILB,1))
        IF (exit_flag.ne.NoError) RETURN
      END IF

      IF (tl_LBC(isouth,isVbar,ng)%acquire) THEN
        CALL get_ngfld (ng, iTLM, idV2bc(isouth), BRY(ng)%ncid,         &
     &                  1, BRY(ng), update(1),                          &
     &                  ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,                 &
     &                  BOUNDARY(ng) % vbarG_south(ILB,1))
        IF (exit_flag.ne.NoError) RETURN
      END IF

      IF (tl_LBC(inorth,isUbar,ng)%acquire) THEN
        CALL get_ngfld (ng, iTLM, idU2bc(inorth), BRY(ng)%ncid,         &
     &                  1, BRY(ng), update(1),                          &
     &                  ILB, IUB, 1, 2, 1, Lm(ng)+1, 1,                 &
     &                  BOUNDARY(ng) % ubarG_north(ILB,1))
        IF (exit_flag.ne.NoError) RETURN
      END IF

      IF (tl_LBC(inorth,isVbar,ng)%acquire) THEN
        CALL get_ngfld (ng, iTLM, idV2bc(inorth), BRY(ng)%ncid,         &
     &                  1, BRY(ng), update(1),                          &
     &                  ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,                 &
     &                  BOUNDARY(ng) % vbarG_north(ILB,1))
        IF (exit_flag.ne.NoError) RETURN
      END IF
# endif

# ifdef SOLVE3D
#  ifndef ANA_M3OBC
!
      IF (tl_LBC(iwest,isUvel,ng)%acquire) THEN
        CALL get_ngfld (ng, iTLM, idU3bc(iwest), BRY(ng)%ncid,          &
     &                  1, BRY(ng), update(1),                          &
     &                  JLB, JUB, N(ng), 2, 0, Mm(ng)+1, N(ng),         &
     &                  BOUNDARY(ng) % uG_west(JLB,1,1))
        IF (exit_flag.ne.NoError) RETURN
      END IF

      IF (tl_LBC(iwest,isVvel,ng)%acquire) THEN
        CALL get_ngfld (ng, iTLM, idV3bc(iwest), BRY(ng)%ncid,          &
     &                  1, BRY(ng), update(1),                          &
     &                  JLB, JUB, N(ng), 2, 1, Mm(ng)+1, N(ng),         &
     &                  BOUNDARY(ng) % vG_west(JLB,1,1))
        IF (exit_flag.ne.NoError) RETURN
      END IF

      IF (tl_LBC(ieast,isUvel,ng)%acquire) THEN
        CALL get_ngfld (ng, iTLM, idU3bc(ieast), BRY(ng)%ncid,          &
     &                  1, BRY(ng), update(1),                          &
     &                  JLB, JUB, N(ng), 2, 0, Mm(ng)+1, N(ng),         &
     &                  BOUNDARY(ng) % uG_east(JLB,1,1))
        IF (exit_flag.ne.NoError) RETURN
      END IF

      IF (tl_LBC(ieast,isVvel,ng)%acquire) THEN
        CALL get_ngfld (ng, iTLM, idV3bc(ieast), BRY(ng)%ncid,          &
     &                  1, BRY(ng), update(1),                          &
     &                  JLB, JUB, N(ng), 2, 1, Mm(ng)+1, N(ng),         &
     &                  BOUNDARY(ng) % vG_east(JLB,1,1))
        IF (exit_flag.ne.NoError) RETURN
      END IF

      IF (tl_LBC(isouth,isUvel,ng)%acquire) THEN
        CALL get_ngfld (ng, iTLM, idU3bc(isouth), BRY(ng)%ncid,         &
     &                  1, BRY(ng), update(1),                          &
     &                  ILB, IUB, N(ng), 2, 1, Lm(ng)+1, N(ng),         &
     &                  BOUNDARY(ng) % uG_south(ILB,1,1))
        IF (exit_flag.ne.NoError) RETURN
      END IF

      IF (tl_LBC(isouth,isVvel,ng)%acquire) THEN
        CALL get_ngfld (ng, iTLM, idV3bc(isouth), BRY(ng)%ncid,         &
     &                  1, BRY(ng), update(1),                          &
     &                  ILB, IUB, N(ng), 2, 0, Lm(ng)+1, N(ng),         &
     &                  BOUNDARY(ng) % vG_south(ILB,1,1))
        IF (exit_flag.ne.NoError) RETURN
      END IF

      IF (tl_LBC(inorth,isUvel,ng)%acquire) THEN
        CALL get_ngfld (ng, iTLM, idU3bc(inorth), BRY(ng)%ncid,         &
     &                  1, BRY(ng), update(1),                          &
     &                  ILB, IUB, N(ng), 2, 1, Lm(ng)+1, N(ng),         &
     &                  BOUNDARY(ng) % uG_north(ILB,1,1))
        IF (exit_flag.ne.NoError) RETURN
      END IF

      IF (tl_LBC(inorth,isVvel,ng)%acquire) THEN
        CALL get_ngfld (ng, iTLM, idV3bc(inorth), BRY(ng)%ncid,         &
     &                  1, BRY(ng), update(1),                          &
     &                  ILB, IUB, N(ng), 2, 0, Lm(ng)+1, N(ng),         &
     &                  BOUNDARY(ng) % vG_north(ILB,1,1))
        IF (exit_flag.ne.NoError) RETURN
      END IF
#  endif

#  ifndef ANA_TOBC
!
      DO i=1,NT(ng)
        IF (tl_LBC(iwest,isTvar(i),ng)%acquire) THEN
          CALL get_ngfld (ng, iTLM, idTbry(iwest,i), BRY(ng)%ncid,      &
     &                    1, BRY(ng), update(1),                        &
     &                    JLB, JUB, N(ng), 2, 0, Mm(ng)+1, N(ng),       &
     &                    BOUNDARY(ng) % tG_west(JLB,1,1,i))
          IF (exit_flag.ne.NoError) RETURN
        END IF
      END DO

      DO i=1,NT(ng)
        IF (tl_LBC(ieast,isTvar(i),ng)%acquire) THEN
          CALL get_ngfld (ng, iTLM, idTbry(ieast,i), BRY(ng)%ncid,      &
     &                    1, BRY(ng), update(1),                        &
     &                    JLB, JUB, N(ng), 2, 0, Mm(ng)+1, N(ng),       &
     &                    BOUNDARY(ng) % tG_east(JLB,1,1,i))
          IF (exit_flag.ne.NoError) RETURN
        END IF
      END DO

      DO i=1,NT(ng)
        IF (tl_LBC(isouth,isTvar(i),ng)%acquire) THEN
          CALL get_ngfld (ng, iTLM, idTbry(isouth,i), BRY(ng)%ncid,     &
     &                    1, BRY(ng), update(1),                        &
     &                    ILB, IUB, N(ng), 2, 0, Lm(ng)+1, N(ng),       &
     &                    BOUNDARY(ng) % tG_south(ILB,1,1,i))
          IF (exit_flag.ne.NoError) RETURN
        END IF
      END DO

      DO i=1,NT(ng)
        IF (tl_LBC(inorth,isTvar(i),ng)%acquire) THEN
          CALL get_ngfld (ng, iTLM, idTbry(inorth,i), BRY(ng)%ncid,     &
     &                    1, BRY(ng), update(1),                        &
     &                    ILB, IUB, N(ng), 2, 0, Lm(ng)+1, N(ng),       &
     &                    BOUNDARY(ng) % tG_north(ILB,1,1,i))
          IF (exit_flag.ne.NoError) RETURN
        END IF
      END DO
#  endif
# endif

# ifdef CLM_FILE
!
!=======================================================================
!  Read in climatology data from CLIMATOLOGY NetCDF file.
!=======================================================================
!
#  if !defined ANA_SSH && defined ZCLIMATOLOGY
      CALL get_2dfld (ng, iTLM, idSSHc, CLM(ng)%ncid,                   &
     &                1, CLM(ng), update(1),                            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#   endif
     &                CLIMA(ng) % sshG(LBi,LBj,1))
      IF (exit_flag.ne.NoError) RETURN
#  endif

#  if !defined ANA_M2CLIMA && defined M2CLIMATOLOGY
      CALL get_2dfld (ng, iTLM, idUbcl, CLM(ng)%ncid,                   &
     &                1, CLM(ng), update(1),                            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % umask(LBi,LBj),                        &
#   endif
     &                CLIMA(ng) % ubarclmG(LBi,LBj,1))
      IF (exit_flag.ne.NoError) RETURN

      CALL get_2dfld (ng, iTLM, idVbcl, CLM(ng)%ncid,                   &
     &                1, CLM(ng), update(1),                            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % vmask(LBi,LBj),                        &
#   endif
     &                CLIMA(ng) % vbarclmG(LBi,LBj,1))
      IF (exit_flag.ne.NoError) RETURN
#  endif

#  ifdef SOLVE3D
#   if !defined ANA_TCLIMA && defined TCLIMATOLOGY
      DO i=1,NAT
        CALL get_3dfld (ng, iTLM, idTclm(i), CLM(ng)%ncid,              &
     &                  1, CLM(ng), update(1),                          &
     &                  LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,             &
#    ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
#    endif
     &                  CLIMA(ng) % tclmG(LBi,LBj,1,1,i))
        IF (exit_flag.ne.NoError) RETURN
      END DO
#   endif

#   if !defined ANA_M3CLIMA && defined M3CLIMATOLOGY
      CALL get_3dfld (ng, iTLM, idUclm, CLM(ng)%ncid,                   &
     &                1, CLM(ng), update(1),                            &
     &                LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,               &
#    ifdef MASKING
     &                GRID(ng) % umask(LBi,LBj),                        &
#    endif
     &                CLIMA(ng) % uclmG(LBi,LBj,1,1))
      IF (exit_flag.ne.NoError) RETURN

      CALL get_3dfld (ng, iTLM, idVclm, CLM(ng)%ncid,                   &
     &                1, CLM(ng), update(1),                            &
     &                LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,               &
#    ifdef MASKING
     &                GRID(ng) % vmask(LBi,LBj),                        &
#    endif
     &                CLIMA(ng) % vclmG(LBi,LBj,1,1))
      IF (exit_flag.ne.NoError) RETURN
#   endif
#  endif
# endif

# ifdef FORWARD_READ
!
!-----------------------------------------------------------------------
!  Read in forward state solution.
!-----------------------------------------------------------------------
!
!  Read in free-surface.
!
      CALL get_2dfld (ng, iTLM, idFsur, FWD(ng)%ncid,                   &
     &                1, FWD(ng), update(1),                            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#  endif
     &                OCEAN(ng) % zetaG(LBi,LBj,1))
      IF (exit_flag.ne.NoError) RETURN
!
!  Read 2D momentum.
!
      CALL get_2dfld (ng, iTLM, idUbar, FWD(ng)%ncid,                   &
     &                1, FWD(ng), update(1),                            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % umask(LBi,LBj),                        &
#  endif
     &                OCEAN(ng) % ubarG(LBi,LBj,1))
      IF (exit_flag.ne.NoError) RETURN

      CALL get_2dfld (ng, iTLM, idVbar, FWD(ng)%ncid,                   &
     &                1, FWD(ng), update(1),                            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#  ifdef MASKING
     &                GRID(ng) % vmask(LBi,LBj),                        &
#  endif
     &                OCEAN(ng) % vbarG(LBi,LBj,1))
      IF (exit_flag.ne.NoError) RETURN

#  ifdef FORWARD_RHS
!
!  Read in variables associated with 2D right-hand-side terms.
!
      CALL get_2dfld (ng, iTLM, idRzet, FWD(ng)%ncid,                   &
     &                1, FWD(ng), update(1),                            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#   endif
     &                OCEAN(ng) % rzetaG(LBi,LBj,1))
      IF (exit_flag.ne.NoError) RETURN

      CALL get_2dfld (ng, iTLM, idRu2d, FWD(ng)%ncid,                   &
     &                1, FWD(ng), update(1),                            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % umask(LBi,LBj),                        &
#   endif
     &                OCEAN(ng) % rubarG(LBi,LBj,1))
      IF (exit_flag.ne.NoError) RETURN

      CALL get_2dfld (ng, iTLM, idRv2d, FWD(ng)%ncid,                   &
     &                1, FWD(ng), update(1),                            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % vmask(LBi,LBj),                        &
#   endif
     &                OCEAN(ng) % rvbarG(LBi,LBj,1))
      IF (exit_flag.ne.NoError) RETURN
#  endif

#  ifdef SOLVE3D
!
!  Read in variables associated with time-averaged 2D momentum terms.
!
      CALL get_2dfld (ng, iTLM, idUfx1, FWD(ng)%ncid,                   &
     &                1, FWD(ng), update(1),                            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % umask(LBi,LBj),                        &
#   endif
     &                COUPLING(ng) % DU_avg1G(LBi,LBj,1))
      IF (exit_flag.ne.NoError) RETURN

      CALL get_2dfld (ng, iTLM, idUfx2, FWD(ng)%ncid,                   &
     &                1, FWD(ng), update(1),                            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % umask(LBi,LBj),                        &
#   endif
     &                COUPLING(ng) % DU_avg2G(LBi,LBj,1))
      IF (exit_flag.ne.NoError) RETURN

      CALL get_2dfld (ng, iTLM, idVfx1, FWD(ng)%ncid,                   &
     &                1, FWD(ng), update(1),                            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % vmask(LBi,LBj),                        &
#   endif
     &                COUPLING(ng) % DV_avg1G(LBi,LBj,1))
      IF (exit_flag.ne.NoError) RETURN

      CALL get_2dfld (ng, iTLM, idVfx2, FWD(ng)%ncid,                   &
     &                1, FWD(ng), update(1),                            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#   ifdef MASKING
     &                GRID(ng) % vmask(LBi,LBj),                        &
#   endif
     &                COUPLING(ng) % DV_avg2G(LBi,LBj,1))
      IF (exit_flag.ne.NoError) RETURN
!
!  Read in 3D momentum.
!
      CALL get_3dfld (ng, iTLM, idUvel, FWD(ng)%ncid,                   &
     &                1, FWD(ng), update(1),                            &
     &                LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,               &
#   ifdef MASKING
     &                GRID(ng) % umask(LBi,LBj),                        &
#   endif
     &                OCEAN(ng) % uG(LBi,LBj,1,1))
      IF (exit_flag.ne.NoError) RETURN

      CALL get_3dfld (ng, iTLM, idVvel, FWD(ng)%ncid,                   &
     &                1, FWD(ng), update(1),                            &
     &                LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,               &
#   ifdef MASKING
     &                GRID(ng) % vmask(LBi,LBj),                        &
#   endif
     &                OCEAN(ng) % vG(LBi,LBj,1,1))
      IF (exit_flag.ne.NoError) RETURN

#   ifdef FORWARD_RHS
!
!  Read in variables associated with 3D momentum right-hand-side terms.
!
      CALL get_2dfld (ng, iTLM, idRuct, FWD(ng)%ncid,                   &
     &                1, FWD(ng), update(1),                            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#    ifdef MASKING
     &                GRID(ng) % umask(LBi,LBj),                        &
#    endif
     &                COUPLING(ng) % rufrcG(LBi,LBj,1))
      IF (exit_flag.ne.NoError) RETURN

      CALL get_2dfld (ng, iTLM, idRvct, FWD(ng)%ncid,                   &
     &                1, FWD(ng), update(1),                            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#    ifdef MASKING
     &                GRID(ng) % vmask(LBi,LBj),                        &
#    endif
     &                COUPLING(ng) % rvfrcG(LBi,LBj,1))
      IF (exit_flag.ne.NoError) RETURN

      CALL get_3dfld (ng, iTLM, idRu3d, FWD(ng)%ncid,                   &
     &                1, FWD(ng), update(1),                            &
     &                LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,               &
#    ifdef MASKING
     &                GRID(ng) % umask(LBi,LBj),                        &
#    endif
     &                OCEAN(ng) % ruG(LBi,LBj,1,1))
      IF (exit_flag.ne.NoError) RETURN

      CALL get_3dfld (ng, iTLM, idRv3d, FWD(ng)%ncid,                   &
     &                1, FWD(ng), update(1),                            &
     &                LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,               &
#    ifdef MASKING
     &                GRID(ng) % vmask(LBi,LBj),                        &
#    endif
     &                OCEAN(ng) % rvG(LBi,LBj,1,1))
      IF (exit_flag.ne.NoError) RETURN
#   endif
!
!  Read in 3D tracers.
!
      DO i=1,NT(ng)
        CALL get_3dfld (ng, iTLM, idTvar(i), FWD(ng)%ncid,              &
     &                  1, FWD(ng), update(1),                          &
     &                  LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,             &
#   ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
#   endif
     &                  OCEAN(ng) % tG(LBi,LBj,1,1,i))
        IF (exit_flag.ne.NoError) RETURN
      END DO

#   ifdef FORWARD_MIXING
!
!  Read in vertical mixing variables.
!
      DO i=1,NAT
        scale=Fscale(idDiff(i),ng)                    ! save and rescale
        Fscale(idDiff(i),ng)=tl_Akt_fac(i,ng)
        CALL get_3dfld (ng, iTLM, idDiff(i), FWD(ng)%ncid,              &
     &                  1, FWD(ng), update(1),                          &
     &                  LBi, UBi, LBj, UBj, 0, N(ng), 2, 1,             &
#    ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
#    endif
     &                  MIXING(ng) % AktG(LBi,LBj,0,1,i))
        Fscale(idDiff(i),ng)=scale
        IF (exit_flag.ne.NoError) RETURN
      END DO

      scale=Fscale(idVvis,ng)                         ! save and rescale
      Fscale(idVvis,ng)=tl_Akv_fac(ng)
      CALL get_3dfld (ng, iTLM, idVvis, FWD(ng)%ncid,                   &
     &                1, FWD(ng), update(1),                            &
     &                LBi, UBi, LBj, UBj, 0, N(ng), 2, 1,               &
#    ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#    endif
     &                MIXING(ng) % AkvG(LBi,LBj,0,1))
      Fscale(idVvis,ng)=scale
      IF (exit_flag.ne.NoError) RETURN
#   endif

#   if defined MY25_MIXING_NOT_YET || defined GLS_MIXING_NOT_YET
!
!  Read in turbulent kinetic energy.
!
      CALL get_3dfld (ng, iTLM, idMtke, FWD(ng)%ncid,                   &
     &                1, FWD(ng), update(1),                            &
     &                LBi, UBi, LBj, UBj, 0, N(ng), 2, 1,               &
#    ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#    endif
     &                MIXING(ng) % tkeG(LBi,LBj,0,1))
      IF (exit_flag.ne.NoError) RETURN
!
!  Read in turbulent kinetic energy times length scale.
!
      CALL get_3dfld (ng, iTLM, idMtls, FWD(ng)%ncid,                   &
     &                1, FWD(ng), update(1),                            &
     &                LBi, UBi, LBj, UBj, 0, N(ng), 2, 1,               &
#    ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#    endif
     &                MIXING(ng) % glsG(LBi,LBj,0,1))
      IF (exit_flag.ne.NoError) RETURN
!
!  Read in vertical mixing length scale.
!
      CALL get_3dfld (ng, iTLM, idVmLS, FWD(ng)%ncid,                   &
     &                1, FWD(ng), update(1),                            &
     &                LBi, UBi, LBj, UBj, 0, N(ng), 2, 1,               &
#    ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#    endif
     &                MIXING(ng) % LscaleG(LBi,LBj,0,1))
      IF (exit_flag.ne.NoError) RETURN
!
!  Read in vertical mixing coefficient for turbulent kinetic energy.
!
      CALL get_3dfld (ng, iTLM, idVmKK, FWD(ng)%ncid,                   &
     &                1, FWD(ng), update(1),                            &
     &                LBi, UBi, LBj, UBj, 0, N(ng), 2, 1,               &
#    ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#    endif
     &                MIXING(ng) % AkkG(LBi,LBj,0,1))
      IF (exit_flag.ne.NoError) RETURN

#    ifdef GLS_MIXING_NOT_YET
!
!  Read in vertical mixing coefficient for turbulent length scale.
!
      CALL get_3dfld (ng, iTLM, idVmKP, FWD(ng)%ncid,                   &
     &                1, FWD(ng), update(1),                            &
     &                LBi, UBi, LBj, UBj, 0, N(ng), 2, 1,               &
#     ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#     endif
     &                MIXING(ng) % AkpG(LBi,LBj,0,1))
      IF (exit_flag.ne.NoError) RETURN
#    endif
#   endif

#   ifdef LMD_MIXING_NOT_YET
!
!  Read in depth of surface oceanic boundary layer.
!
      CALL get_2dfld (ng, iTLM, idHsbl, FWD(ng)%ncid,                   &
     &                1, FWD(ng), update(1),                            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#    ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#    endif
     &                MIXING(ng) % hsblG(LBi,LBj,1))
      IF (exit_flag.ne.NoError) RETURN
#   endif

#   ifdef LMD_BKPP_NOT_YET
!
!  Read in depth of bottom oceanic boundary layer.
!
      CALL get_2dfld (ng, iTLM, idHbbl, FWD(ng)%ncid,                   &
     &                1, FWD(ng), update(1),                            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
#    ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#    endif
     &                MIXING(ng) % hbblG(LBi,LBj,1))
      IF (exit_flag.ne.NoError) RETURN
#   endif

#   ifdef LMD_NONLOCAL_NOT_YET
!
!  Read in boundary layer nonlocal transport.
!
      DO i=1,NAT
        CALL get_3dfld (ng, iTLM, idGhat(i), FWD(ng)%ncid,              &
     &                  1, FWD(ng), update(1),                          &
     &                  LBi, UBi, LBj, UBj, 0, N(ng), 2, 1,             &
#    ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
#    endif
     &                  MIXING(ng) % ghatsG(LBi,LBj,0,i))
        IF (exit_flag.ne.NoError) RETURN
      END DO
#   endif
#  endif

#  ifdef WEAK_CONSTRAINT
!
!-----------------------------------------------------------------------
!  Read frequent impulse forcing for weak constraint.
!-----------------------------------------------------------------------
!
      IF (FrequentImpulse(ng)) THEN
        CALL get_2dfld (ng, iTLM, idZtlf, TLF(ng)%ncid,                 &
     &                  1, TLF(ng), update(1),                          &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
#   ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
#   endif
     &                  OCEAN(ng) % f_zetaG(LBi,LBj,1))
        IF (exit_flag.ne.NoError) RETURN

#   ifndef SOLVE3D
!
!  Read 2D momentum forcing.
!
        CALL get_2dfld (ng, iTLM, idUbtf, TLF(ng)%ncid,                 &
     &                  1, TLF(ng), update(1),                          &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
#    ifdef MASKING
     &                  GRID(ng) % umask(LBi,LBj),                      &
#    endif
     &                  OCEAN(ng) % f_ubarG(LBi,LBj,1))
        IF (exit_flag.ne.NoError) RETURN

        CALL get_2dfld (ng, iTLM, idVbtf, TLF(ng)%ncid,                 &
     &                  1, TLF(ng), update(1),                          &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
#    ifdef MASKING
     &                  GRID(ng) % vmask(LBi,LBj),                      &
#    endif
     &                  OCEAN(ng) % f_vbarG(LBi,LBj,1))
        IF (exit_flag.ne.NoError) RETURN
#   endif

#   ifdef SOLVE3D
!
!  Read in 3D momentum forcing.
!
        CALL get_3dfld (ng, iTLM, idUtlf, TLF(ng)%ncid,                 &
     &                  1, TLF(ng), update(1),                          &
     &                  LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,             &
#    ifdef MASKING
     &                  GRID(ng) % umask(LBi,LBj),                      &
#    endif
     &                  OCEAN(ng) % f_uG(LBi,LBj,1,1))
        IF (exit_flag.ne.NoError) RETURN

        CALL get_3dfld (ng, iTLM, idVtlf, TLF(ng)%ncid,                 &
     &                  1, TLF(ng), update(1),                          &
     &                  LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,             &
#    ifdef MASKING
     &                  GRID(ng) % vmask(LBi,LBj),                      &
#    endif
     &                  OCEAN(ng) % f_vG(LBi,LBj,1,1))
        IF (exit_flag.ne.NoError) RETURN
!
!  Read in 3D tracers forcing.
!
        DO i=1,NT(ng)
          CALL get_3dfld (ng, iTLM, idTtlf(i), TLF(ng)%ncid,            &
     &                    1, TLF(ng), update(1),                        &
     &                    LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,           &
#    ifdef MASKING
     &                    GRID(ng) % rmask(LBi,LBj),                    &
#    endif
     &                    OCEAN(ng) % f_tG(LBi,LBj,1,1,i))
          IF (exit_flag.ne.NoError) RETURN
        END DO
#   endif
      END IF
#  endif
# endif

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off input data time wall clock.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, iTLM, 3)
# endif

      RETURN
      END SUBROUTINE tl_get_data
#else
      SUBROUTINE tl_get_data
      RETURN
      END SUBROUTINE tl_get_data
#endif
