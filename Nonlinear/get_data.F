#include "cppdefs.h"
      SUBROUTINE get_data (ng)
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine reads in forcing, climatology and assimilation data    !
!  from input NetCDF files.  If there is more than one time-record,    !
!  data  is loaded  into global two-time record arrays.  The actual    !
!  interpolation is carried elsewhere.                                 !
!                                                                      !
!  Currently,  this routine is only executed in serial mode by the     !
!  main thread.                                                        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_boundary
#ifdef CLIMATOLOGY
      USE mod_clima
#endif
      USE mod_forces
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
#if defined ASSIMILATION || defined NUDGING
      USE mod_obs
#endif
      USE mod_scalars
#if defined UV_PSOURCE || defined TS_PSOURCE
      USE mod_sources
#endif
      USE mod_stepping
#if defined SSH_TIDES || defined UV_TIDES
      USE mod_tides
#endif
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
     &         (/ .false., .false., .false. /)

#ifdef OBC
      integer :: ILB, IUB, JLB, JUB
#endif
      integer :: LBi, UBi, LBj, UBj
      integer :: i, is

      real(r8) :: time_save = 0.0_r8

#ifdef OBC
!
!  Lower and upper bounds for nontiled boundary arrays.
!
      ILB=LOWER_BOUND_I
      IUB=UPPER_BOUND_I
      JLB=LOWER_BOUND_J
      JUB=UPPER_BOUND_J
#endif
!
!  Lower and upper bounds for tiled arrays.
!
      LBi=lbound(GRID(ng)%h,DIM=1)
      UBi=ubound(GRID(ng)%h,DIM=1)
      LBj=lbound(GRID(ng)%h,DIM=2)
      UBj=ubound(GRID(ng)%h,DIM=2)

#ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on input data time wall clock.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, 3)
#endif
!
!=======================================================================
!  Read in forcing data from FORCING NetCDF file.
!=======================================================================
!
#if defined SSH_TIDES || defined UV_TIDES
!
!-----------------------------------------------------------------------
!  Tide period, amplitude, phase, and currents.
!-----------------------------------------------------------------------
!
!  Tidal Period.
!
      IF (iic(ng).eq.0) THEN
        NTC(ng)=0
        CALL get_ngfld (ng, idTper, ncfrcid(idTper,ng),                 &
     &                  nFfiles(ng), frcname(1,ng),                     &
     &                  update(1), 1, MTC, 1, 1, 1, NTC(ng), 1,         &
     &                  TIDES(ng) % Tperiod(1))
      END IF
#endif
#ifdef SSH_TIDES
!
!  Tidal elevation amplitude and phase. In order to read data as a
!  function of tidal period, we need to reset the model time variables
!  temporarily.
!
      IF (iic(ng).eq.0) THEN
        time_save=time(ng)
        time(ng)=8640000.0_r8
        tdays(ng)=time(ng)*sec2day
        CALL get_2dfld (ng, idTzam, ncfrcid(idTzam,ng),                 &
     &                  nFfiles(ng), frcname(1,ng),                     &
     &                  update(1), LBi, UBi, LBj, UBj, MTC, NTC(ng),    &
# ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
# endif
     &                  TIDES(ng) % SSH_Tamp(LBi,LBj,1))
        CALL get_2dfld (ng, idTzph, ncfrcid(idTzph,ng),                 &
     &                  nFfiles(ng), frcname(1,ng),                     &
     &                  update(1), LBi, UBi, LBj, UBj, MTC, NTC(ng),    &
# ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
# endif
     &                  TIDES(ng) % SSH_Tphase(LBi,LBj,1))
        time(ng)=time_save
        tdays(ng)=time(ng)*sec2day
      END IF
#endif
#ifdef UV_TIDES
!
!  Tidal currents angle, phase, major and minor ellipse axis.
!
      IF (iic(ng).eq.0) THEN
        time_save=time(ng)
        time(ng)=8640000.0_r8
        tdays(ng)=time(ng)*sec2day
        CALL get_2dfld (ng, idTvan, ncfrcid(idTvan,ng),                 &
     &                  nFfiles(ng), frcname(1,ng),                     &
     &                  update(1), LBi, UBi, LBj, UBj, MTC, NTC(ng),    &
# ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
# endif
     &                  TIDES(ng) % UV_Tangle(LBi,LBj,1))
        CALL get_2dfld (ng, idTvph, ncfrcid(idTvph,ng),                 &
     &                  nFfiles(ng), frcname(1,ng),                     &
     &                  update(1), LBi, UBi, LBj, UBj, MTC, NTC(ng),    &
# ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
# endif
    &                  TIDES(ng) % UV_Tphase(LBi,LBj,1))
        CALL get_2dfld (ng, idTvma, ncfrcid(idTvma,ng),                 &
     &                  nFfiles(ng), frcname(1,ng),                     &
     &                  update(1), LBi, UBi, LBj, UBj, MTC, NTC(ng),    &
# ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
# endif
     &                  TIDES(ng) % UV_Tmajor(LBi,LBj,1))
        CALL get_2dfld (ng, idTvmi, ncfrcid(idTvmi,ng),                 &
     &                  nFfiles(ng), frcname(1,ng),                     &
     &                  update(1), LBi, UBi, LBj, UBj, MTC, NTC(ng),    &
# ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
# endif
     &                  TIDES(ng) % UV_Tminor(LBi,LBj,1))
        time(ng)=time_save
        tdays(ng)=time(ng)*sec2day
      END IF
#endif

#if !defined ANA_PSOURCE && (defined UV_PSOURCE || \
                             defined TS_PSOURCE)
!
!-----------------------------------------------------------------------
!  Point Sources/Sinks position, direction, special flag, and mass
!  transport nondimensional shape profile.  Point sources are at U-
!  and V-points.
!-----------------------------------------------------------------------
!
      IF (iic(ng).eq.0) THEN
        CALL get_ngfld (ng, idRxpo, ncfrcid(idRxpo,ng),                 &
     &                  nFfiles(ng), frcname(1,ng),                     &
     &                  update(1), 1, Msrc, 1, 1, 1, Nsrc(ng), 1,       &
     &                  SOURCES(ng) % Xsrc(1))
        CALL get_ngfld (ng, idRepo, ncfrcid(idRepo,ng),                 &
     &                  nFfiles(ng), frcname(1,ng),                     &
     &                  update(1), 1, Msrc, 1, 1, 1, Nsrc(ng), 1,       &
     &                  SOURCES(ng) % Ysrc(1))
        CALL get_ngfld (ng, idRdir, ncfrcid(idRdir,ng),                 &
     &                  nFfiles(ng), frcname(1,ng),                     &
     &                  update(1), 1, Msrc, 1, 1, 1, Nsrc(ng), 1,       &
     &                  SOURCES(ng) % Dsrc(1))
        CALL get_ngfld (ng, idRvsh, ncfrcid(idRvsh,ng),                 &
     &                  nFfiles(ng), frcname(1,ng),                     &
     &                  update(1), 1, Msrc, N, 1, 1, Nsrc(ng), N,       &
     &                  SOURCES(ng) % Qshape(1,1))
# ifdef TS_PSOURCE
        CALL get_ngfld (ng, idRflg, ncfrcid(idRflg,ng),                 &
     &                  nFfiles(ng), frcname(1,ng),                     &
     &                  update(1), 1, Msrc, 1, 1, 1, Nsrc(ng), 1,       &
     &                  SOURCES(ng) % Fsrc(1))
        IF (update(1)) THEN
          DO i=1,NT(ng)
            SOURCES(ng)%Ltracer(i)=.false.
            DO is=1,Nsrc(ng)
              SOURCES(ng)%Lsrc(is,i)=.false.
            END DO
          END DO
          DO is=1,Nsrc(ng)
            IF (SOURCES(ng)%Fsrc(is).eq.1.0_r8) THEN
              SOURCES(ng)%Lsrc(is,itemp)=.true.
              SOURCES(ng)%Ltracer(itemp)=.true.
            ELSE IF (SOURCES(ng)%Fsrc(is).eq.2.0_r8) THEN
              SOURCES(ng)%Lsrc(is,isalt)=.true.
              SOURCES(ng)%Ltracer(isalt)=.true.
            ELSE IF (SOURCES(ng)%Fsrc(is).eq.3.0_r8) THEN
              SOURCES(ng)%Lsrc(is,itemp)=.true.
              SOURCES(ng)%Lsrc(is,isalt)=.true.
              SOURCES(ng)%Ltracer(itemp)=.true.
              SOURCES(ng)%Ltracer(isalt)=.true.
            END IF
          END DO
        END IF
# endif
        DO is=1,Nsrc(ng)
          SOURCES(ng)%Isrc(is)=                                         &
     &                MAX(1,MIN(NINT(SOURCES(ng)%Xsrc(is)),Lm(ng)+1))
          SOURCES(ng)%Jsrc(is)=                                         &
     &                MAX(1,MIN(NINT(SOURCES(ng)%Ysrc(is)),Mm(ng)+1))
        END DO
      END IF
# ifdef UV_PSOURCE
!
!  Point Source/Sink vertically integrated mass transport.
!
      CALL get_ngfld (ng, idRtra, ncfrcid(idRtra,ng),                   &
     &                nFfiles(ng), frcname(1,ng),                       &
     &                update(1), 1, Msrc, 1, 2, 1, Nsrc(ng), 1,         &
     &                SOURCES(ng) % QbarG(1,1))
# endif
# if defined TS_PSOURCE && defined SOLVE3D
!
!  Tracer Sources/Sinks.
!
      DO i=1,NT(ng)
        IF (SOURCES(ng)%Ltracer(i)) THEN
          CALL get_ngfld (ng, idRtrc(i), ncfrcid(idRtrc(i),ng),         &
     &                    nFfiles(ng), frcname(1,ng),                   &
     &                    update(1), 1, Msrc, N, 2, 1, Nsrc(ng), N,     &
     &                    SOURCES(ng) % TsrcG(1,1,1,i))
        END IF
      END DO
# endif
#endif

#if !defined ANA_WINDS && (defined BULK_FLUXES || defined ECOSIM)
!
!-----------------------------------------------------------------------
!  Surface wind components.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, idUair, ncfrcid(idUair,ng),                   &
     &                nFfiles(ng), frcname(1,ng),                       &
     &                update(1), LBi, UBi, LBj, UBj, 2, 1,              &
# ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
# endif
     &                FORCES(ng) % UwindG(LBi,LBj,1))
      CALL get_2dfld (ng , idVair, ncfrcid(idVair,ng),                  &
     &                nFfiles(ng), frcname(1,ng),                       &
     &                update(1), LBi, UBi, LBj, UBj, 2, 1,              &
# ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
# endif
     &                FORCES(ng) % VwindG(LBi,LBj,1))
#endif

#ifndef AIR_OCEAN
# if !defined ANA_SMFLUX && !defined BULK_FLUXES
!
!-----------------------------------------------------------------------
!  Surface wind stress components.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, idUsms, ncfrcid(idUsms,ng),                   &
     &                nFfiles(ng), frcname(1,ng),                       &
     &                update(1), LBi, UBi, LBj, UBj, 2, 1,              &
#  ifdef MASKING
     &                GRID(ng) % umask(LBi,LBj),                        &
#  endif
     &                FORCES(ng) % sustrG(LBi,LBj,1))
      CALL get_2dfld (ng, idVsms, ncfrcid(idVsms,ng),                   &
     &                nFfiles(ng), frcname(1,ng),                       &
     &                update(1), LBi, UBi, LBj, UBj, 2, 1,              &
#  ifdef MASKING
     &                GRID(ng) % vmask(LBi,LBj),                        &
#  endif
     &                FORCES(ng) % svstrG(LBi,LBj,1))
# endif
#endif

#ifdef SOLVE3D

# if !defined ANA_WWAVE && defined BBL_MODEL
!
!-----------------------------------------------------------------------
!  Surface wind induced wave amplitude, direction and period.
!-----------------------------------------------------------------------
!
# ifdef SWAN
      CALL get_2dfld (ng, idWorb, ncfrcid(idWorb,ng),                   &
     &                nFfiles(ng), frcname(1,ng),                       &
     &                update(1), LBi, UBi, LBj, UBj, 2, 1,              &
#  ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#  endif
     &                FORCES(ng) % UB_swanG(LBi,LBj,1))
# else
      CALL get_2dfld (ng, idWamp, ncfrcid(idWamp,ng),                   &
     &                nFfiles(ng), frcname(1,ng),                       &
     &                update(1), LBi, UBi, LBj, UBj, 2, 1,              &
#  ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#  endif
     &                FORCES(ng) % AwaveG(LBi,LBj,1))
# endif
      CALL get_2dfld (ng, idWdir, ncfrcid(idWdir,ng),                   &
     &                nFfiles(ng), frcname(1,ng),                       &
     &                update(1), LBi, UBi, LBj, UBj, 2, 1,              &
#  ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#  endif
     &                FORCES(ng) % DwaveG(LBi,LBj,1))
      CALL get_2dfld (ng, idWper, ncfrcid(idWper,ng),                   &
     &                nFfiles(ng), frcname(1,ng),                       &
     &                update(1), LBi, UBi, LBj, UBj, 2, 1,              &
#  ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#  endif
     &                FORCES(ng) % PwaveG(LBi,LBj,1))
# endif

# if !defined ANA_CLOUD && defined CLOUDS
!
!-----------------------------------------------------------------------
!  Cloud fraction.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, idCfra, ncfrcid(idCfra,ng),                   &
     &                nFfiles(ng), frcname(1,ng),                       &
     &                update(1), LBi, UBi, LBj, UBj, 2, 1,              &
#  ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#  endif
     &                FORCES(ng) % cloudG(LBi,LBj,1))
# endif

# if !defined ANA_SRFLUX && defined SHORTWAVE
!
!-----------------------------------------------------------------------
!  Surface solar shortwave radiation.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, idSrad, ncfrcid(idSrad,ng),                   &
     &                nFfiles(ng), frcname(1,ng),                       &
     &                update(1), LBi, UBi, LBj, UBj, 2, 1,              &
#  ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#  endif
     &                FORCES(ng) % srflxG(LBi,LBj,1))
# endif

# if defined BULK_FLUXES && !defined LONGWAVE
!
!-----------------------------------------------------------------------
!  Surface net longwave radiation.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, idLrad, ncfrcid(idLrad,ng),                   &
     &                nFfiles(ng), frcname(1,ng),                       &
     &                update(1), LBi, UBi, LBj, UBj, 2, 1,              &
#  ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#  endif
     &                FORCES(ng) % lrflxG(LBi,LBj,1))
# endif

# if defined BULK_FLUXES && defined LONGWAVE_OUT
!
!-----------------------------------------------------------------------
!  Surface downwelling longwave radiation.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, idLdwn, ncfrcid(idLdwn,ng),                   &
     &                nFfiles(ng), frcname(1,ng),                       &
     &                update(1), LBi, UBi, LBj, UBj, 2, 1,              &
#  ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#  endif
     &                FORCES(ng) % lrflxG(LBi,LBj,1))
# endif

# if !defined ANA_PAIR && (defined BULK_FLUXES || defined ECOSIM)
!
!-----------------------------------------------------------------------
!  Surface air pressure.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, idPair, ncfrcid(idPair,ng),                   &
     &                nFfiles(ng), frcname(1,ng),                       &
     &                update(1), LBi, UBi, LBj, UBj, 2, 1,              &
#  ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#  endif
     &                FORCES(ng) % PairG(LBi,LBj,1))
# endif

# if !defined ANA_TAIR && \
    ( defined BULK_FLUXES || defined ECOSIM || \
     (defined SHORTWAVE && defined ANA_SRFLUX && defined ALBEDO) )
!
!-----------------------------------------------------------------------
!  Surface air temperature.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, idTair, ncfrcid(idTair,ng),                   &
     &                nFfiles(ng), frcname(1,ng),                       &
     &                update(1), LBi, UBi, LBj, UBj, 2, 1,              &
#  ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#  endif
     &                FORCES(ng) % TairG(LBi,LBj,1))
# endif

# if !defined ANA_HUMIDITY && (defined BULK_FLUXES || defined ECOSIM)
!
!-----------------------------------------------------------------------
!  Surface air humidity.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, idQair, ncfrcid(idQair,ng),                   &
     &                nFfiles(ng), frcname(1,ng),                       &
     &                update(1), LBi, UBi, LBj, UBj, 2, 1,              &
#  ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#  endif
     &                FORCES(ng) % HairG(LBi,LBj,1))
# endif

# if !defined ANA_RAIN && defined BULK_FLUXES
!
!-----------------------------------------------------------------------
!  Rain fall rate.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, idrain, ncfrcid(idrain,ng),                   &
     &                nFfiles(ng), frcname(1,ng),                       &
     &                update(1), LBi, UBi, LBj, UBj, 2, 1,              &
#  ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#  endif
     &                FORCES(ng) % rainG(LBi,LBj,1))
# endif

# if !defined ANA_STFLUX && !defined BULK_FLUXES
!
!-----------------------------------------------------------------------
!  Surface net heat flux.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, idTsur(itemp), ncfrcid(idTsur(itemp),ng),     &
     &                nFfiles(ng), frcname(1,ng),                       &
     &                update(1), LBi, UBi, LBj, UBj, 2, 1,              &
#  ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#  endif
     &                FORCES(ng) % stflxG(LBi,LBj,1,itemp))
# endif

# if !defined ANA_SST && defined QCORRECTION
!
!-----------------------------------------------------------------------
!  Surface net heat flux correction fields: sea surface temperature
!  (SST) and heat flux sensitivity to SST (dQdSST).
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, idSSTc, ncfrcid(idSSTc,ng),                   &
     &                nFfiles(ng), frcname(1,ng),                       &
     &                update(1), LBi, UBi, LBj, UBj, 2, 1,              &
#  ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#  endif
     &                FORCES(ng) % sstG(LBi,LBj,1))
      CALL get_2dfld (ng, iddQdT, ncfrcid(iddQdT,ng),                   &
     &                nFfiles(ng), frcname(1,ng),                       &
     &                update(1), LBi, UBi, LBj, UBj, 2, 1,              &
#  ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#  endif
     &                FORCES(ng) % dqdtG(LBi,LBj,1))
# endif

# ifndef ANA_BTFLUX
!
!-----------------------------------------------------------------------
!  Bottom net heat flux.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, idTbot(itemp), ncfrcid(idTbot(itemp),ng),     &
     &                nFfiles(ng), frcname(1,ng),                       &
     &                update(1), LBi, UBi, LBj, UBj, 2, 1,              &
#  ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#  endif
     &                FORCES(ng) % btflxG(LBi,LBj,1,itemp))
# endif

# ifdef SALINITY
#  ifndef ANA_SSFLUX
!
!-----------------------------------------------------------------------
!  Surface net freshwater flux: E-P.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, idTsur(isalt), ncfrcid(idTsur(isalt),ng),     &
     &                nFfiles(ng), frcname(1,ng),                       &
     &                update(1), LBi, UBi, LBj, UBj, 2, 1,              &
#   ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#   endif
     &                FORCES(ng) % stflxG(LBi,LBj,1,isalt))
#  endif

#  if !defined ANA_SSS && (defined SCORRECTION || defined SRELAXATION)
!
!-----------------------------------------------------------------------
!  Surface net freshwater flux correction field: sea surface salinity.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, idSSSc, ncfrcid(idSSSc,ng),                   &
     &                nFfiles(ng), frcname(1,ng),                       &
     &                update(1), LBi, UBi, LBj, UBj, 2, 1,              &
#   ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#   endif
     &                FORCES(ng) % sssG(LBi,LBj,1))
#  endif

#  ifndef ANA_BSFLUX
!
!-----------------------------------------------------------------------
!  Bottom net freshwater flux.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, idTbot(isalt), ncfrcid(idTbot(isalt),ng),     &
     &                nFfiles(ng), frcname(1,ng),                       &
     &                update(1), LBi, UBi, LBj, UBj, 2, 1,              &
#   ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#   endif
     &                FORCES(ng) % btflxG(LBi,LBj,1,isalt))
#  endif
# endif

# if defined SEDIMENT || defined BIOLOGY
#  ifndef ANA_SPFLUX
!
!-----------------------------------------------------------------------
!  Passive tracers surface fluxes.
!-----------------------------------------------------------------------
!
      DO i=NAT+1,NT(ng)
        CALL get_2dfld (ng, idTsur(i), ncfrcid(idTsur(i),ng),           &
     &                  nFfiles(ng), frcname(1,ng),                     &
     &                  update(1), LBi, UBi, LBj, UBj, 2, 1,            &
#   ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
#   endif
     &                  FORCES(ng) % stflxG(LBi,LBj,1,i))
      END DO
#  endif

#  ifndef ANA_BPFLUX
!
!-----------------------------------------------------------------------
!  Passive tracers bottom fluxes.
!-----------------------------------------------------------------------
!
      DO i=NAT+1,NT(ng)
        CALL get_2dfld (ng, idTbot(i), ncfrcid(idTbot(i),ng),           &
     &                  nFfiles(ng), frcname(1,ng),                     &
     &                  update(1), LBi, UBi, LBj, UBj, 2, 1,            &
#   ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
#   endif
     &                  FORCES(ng) % btflxG(LBi,LBj,1,i))
      END DO
#  endif
# endif
#endif

#ifdef OBC_DATA
!
!=======================================================================
!  Read in open boundary conditions from BOUNDARY NetCDF file.
!=======================================================================
!
# ifndef ANA_FSOBC
#  ifdef WEST_FSOBC
      CALL get_ngfld (ng, idZbry(iwest), ncbryid(ng), 1, bryname(ng),   &
     &                update(1), JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,        &
     &                BOUNDARY(ng) % zetaG_west(JLB,1))
#  endif
#  ifdef EAST_FSOBC
      CALL get_ngfld (ng, idZbry(ieast), ncbryid(ng), 1, bryname(ng),   &
     &                update(1), JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,        &
     &                BOUNDARY(ng) % zetaG_east(JLB,1))
#  endif
#  ifdef SOUTH_FSOBC
      CALL get_ngfld (ng, idZbry(isouth), ncbryid(ng), 1, bryname(ng),  &
     &                update(1), ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,        &
     &                BOUNDARY(ng) % zetaG_south(ILB,1))
#  endif
#  ifdef NORTH_FSOBC
      CALL get_ngfld (ng, idZbry(inorth), ncbryid(ng), 1, bryname(ng),  &
     &                update(1), ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,        &
     &                BOUNDARY(ng) % zetaG_north(ILB,1))
#  endif
# endif
# ifndef ANA_M2OBC
#  ifdef WEST_M2OBC
      CALL get_ngfld (ng, idU2bc(iwest), ncbryid(ng), 1, bryname(ng),   &
     &                update(1), JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,        &
     &                BOUNDARY(ng) % ubarG_west(JLB,1))
      CALL get_ngfld (ng, idV2bc(iwest), ncbryid(ng), 1, bryname(ng),   &
     &                update(1), JLB, JUB, 1, 2, 1, Mm(ng)+1, 1,        &
     &                BOUNDARY(ng) % vbarG_west(JLB,1))
#  endif
#  ifdef EAST_M2OBC
      CALL get_ngfld (ng, idU2bc(ieast), ncbryid(ng), 1, bryname(ng),   &
     &                update(1), JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,        &
     &                BOUNDARY(ng) % ubarG_east(JLB,1))
      CALL get_ngfld (ng, idV2bc(ieast), ncbryid(ng), 1, bryname(ng),   &
     &                update(1), JLB, JUB, 1, 2, 1, Mm(ng)+1, 1,        &
     &                BOUNDARY(ng) % vbarG_east(JLB,1))
#  endif
#  ifdef SOUTH_M2OBC
      CALL get_ngfld (ng, idU2bc(isouth), ncbryid(ng), 1, bryname(ng),  &
     &                update(1), ILB, IUB, 1, 2, 1, Lm(ng)+1, 1,        &
     &                BOUNDARY(ng) % ubarG_south(ILB,1))
      CALL get_ngfld (ng, idV2bc(isouth), ncbryid(ng), 1, bryname(ng),  &
     &                update(1), ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,        &
     &                BOUNDARY(ng) % vbarG_south(ILB,1))
#  endif
#  ifdef NORTH_M2OBC
      CALL get_ngfld (ng, idU2bc(inorth), ncbryid(ng), 1, bryname(ng),  &
     &                update(1), ILB, IUB, 1, 2, 1, Lm(ng)+1, 1,        &
     &                BOUNDARY(ng) % ubarG_north(ILB,1))
      CALL get_ngfld (ng, idV2bc(inorth), ncbryid(ng), 1, bryname(ng),  &
     &                update(1), ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,        &
     &                BOUNDARY(ng) % vbarG_north(ILB,1))
#  endif
# endif
# ifdef SOLVE3D
#  ifndef ANA_M3OBC
#   ifdef WEST_M3OBC
      CALL get_ngfld (ng, idU3bc(iwest), ncbryid(ng), 1, bryname(ng),   &
     &                update(1), JLB, JUB, N(ng), 2, 0, Mm(ng)+1,       &
     &                N(ng), BOUNDARY(ng) % uG_west(JLB,1,1))
      CALL get_ngfld (ng, idV3bc(iwest), ncbryid(ng), 1, bryname(ng),   &
     &                update(1), JLB, JUB, N(ng), 2, 1, Mm(ng)+1,       &
     &                N(ng), BOUNDARY(ng) % vG_west(JLB,1,1))
#   endif
#   ifdef EAST_M3OBC
      CALL get_ngfld (ng, idU3bc(ieast), ncbryid(ng), 1, bryname(ng),   &
     &                update(1), JLB, JUB, N(ng), 2, 0, Mm(ng)+1,       &
     &                N(ng), BOUNDARY(ng) % uG_east(JLB,1,1))
      CALL get_ngfld (ng, idV3bc(ieast), ncbryid(ng), 1, bryname(ng),   &
     &                update(1), JLB, JUB, N(ng), 2, 1, Mm(ng)+1,       &
     &                N(ng), BOUNDARY(ng) % vG_east(JLB,1,1))
#   endif
#   ifdef SOUTH_M3OBC
      CALL get_ngfld (ng, idU3bc(isouth), ncbryid(ng), 1, bryname(ng),  &
     &                update(1), ILB, IUB, N(ng), 2, 1, Lm(ng)+1,       &
     &                N(ng), BOUNDARY(ng) % uG_south(ILB,1,1))
      CALL get_ngfld (ng, idV3bc(isouth), ncbryid(ng), 1, bryname(ng),  &
     &                update(1), ILB, IUB, N(ng), 2, 0, Lm(ng)+1,       &
     &                N(ng), BOUNDARY(ng) % vG_south(ILB,1,1))
#   endif
#   ifdef NORTH_M3OBC
      CALL get_ngfld (ng, idU3bc(inorth), ncbryid(ng), 1, bryname(ng),  &
     &                update(1), ILB, IUB, N(ng), 2, 1, Lm(ng)+1,       &
     &                N(ng), BOUNDARY(ng) % uG_north(ILB,1,1))
      CALL get_ngfld (ng, idV3bc(inorth), ncbryid(ng), 1, bryname(ng),  &
     &                update(1), ILB, IUB, N(ng), 2, 0, Lm(ng)+1,       &
     &                N(ng), BOUNDARY(ng) % vG_north(ILB,1,1))
#   endif
#  endif
#  ifndef ANA_TOBC
#   ifdef WEST_TOBC
      DO i=1,NT(ng)
        CALL get_ngfld (ng, idTbry(iwest,i), ncbryid(ng),               &
     &                  1, bryname(ng),                                 &
     &                  update(1), JLB, JUB, N(ng), 2, 0, Mm(ng)+1,     &
     &                  N(ng), BOUNDARY(ng) % tG_west(JLB,1,1,i))
      END DO
#   endif
#   ifdef EAST_TOBC
      DO i=1,NT(ng)
        CALL get_ngfld (ng, idTbry(ieast,i), ncbryid(ng),               &
     &                  1, bryname(ng),                                 &
     &                  update(1), JLB, JUB, N(ng), 2, 0, Mm(ng)+1,     &
     &                  N(ng), BOUNDARY(ng) % tG_east(JLB,1,1,i))
      END DO
#   endif
#   ifdef SOUTH_TOBC
      DO i=1,NT(ng)
        CALL get_ngfld (ng, idTbry(isouth,i), ncbryid(ng),              &
     &                  1, bryname(ng),                                 &
     &                  update(1), ILB, IUB, N(ng), 2, 0, Lm(ng)+1,     &
     &                  N(ng), BOUNDARY(ng) % tG_south(ILB,1,1,i))
      END DO
#   endif
#   ifdef NORTH_TOBC
      DO i=1,NT(ng)
        CALL get_ngfld (ng, idTbry(inorth,i), ncbryid(ng),              &
     &                  1, bryname(ng),                                 &
     &                  update(1), ILB, IUB, N(ng), 2, 0, Lm(ng)+1,     &
     &                  N(ng), BOUNDARY(ng) % tG_north(ILB,1,1,i))
      END DO
#   endif
#  endif
# endif
#endif
#if defined ZCLIMATOLOGY || defined M2CLIMATOLOGY || \
    defined TCLIMATOLOGY || defined M3CLIMATOLOGY
!
!=======================================================================
!  Read in climatology data from CLIMATOLOGY NetCDF file.
!=======================================================================
!
# if !defined ANA_SSH && defined ZCLIMATOLOGY
      CALL get_2dfld (ng, idSSHc, ncclmid(ng), 1, clmname(ng),          &
     &                update(1), LBi, UBi, LBj, UBj, 2, 1,              &
#  ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#  endif
     &                CLIMA(ng) % sshG(LBi,LBj,1))
# endif
# if !defined ANA_M2CLIMA && defined M2CLIMATOLOGY
      CALL get_2dfld (ng, idUbcl, ncclmid(ng), 1, clmname(ng),          &
     &                update(1), LBi, UBi, LBj, UBj, 2, 1,              &
#  ifdef MASKING
     &                GRID(ng) % umask(LBi,LBj),                        &
#  endif
     &                CLIMA(ng) % ubarclmG(LBi,LBj,1))
      CALL get_2dfld (ng, idVbcl, ncclmid(ng), 1, clmname(ng),          &
     &                update(1), LBi, UBi, LBj, UBj, 2, 1,              &
#  ifdef MASKING
     &                GRID(ng) % vmask(LBi,LBj),                        &
#  endif
     &                CLIMA(ng) % vbarclmG(LBi,LBj,1))
# endif
# ifdef SOLVE3D
#  if !defined ANA_TCLIMA && defined TCLIMATOLOGY
      DO i=1,NAT
        CALL get_3dfld (ng, idTclm(i), ncclmid(ng), 1, clmname(ng),     &
     &                  update(1), LBi, UBi, LBj, UBj, 2, 1,            &
#   ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
#   endif
     &                  CLIMA(ng) % tclmG(LBi,LBj,1,1,i))
      END DO
#  endif
#  if !defined ANA_M3CLIMA && defined M3CLIMATOLOGY
      CALL get_3dfld (ng, idUclm, ncclmid(ng), 1, clmname(ng),          &
     &                update(1), LBi, UBi, LBj, UBj, 2, 1,              &
#   ifdef MASKING
     &                GRID(ng) % umask(LBi,LBj),                        &
#   endif
     &                CLIMA(ng) % uclmG(LBi,LBj,1,1))
      CALL get_3dfld (ng, idVclm, ncclmid(ng), 1, clmname(ng),          &
     &                update(1), LBi, UBi, LBj, UBj, 2, 1,              &
#   ifdef MASKING
     &                GRID(ng) % vmask(LBi,LBj),                        &
#   endif
     &                CLIMA(ng) % vclmG(LBi,LBj,1,1))
#  endif
# endif
#endif
#if defined ASSIMILATION || defined NUDGING
!
!=======================================================================
!  Read in assimilation data from ASSIMILATION NetCDF files.
!=======================================================================
!
# if defined NUDGING_SSH || defined ASSIMILATION_SSH
!
!-----------------------------------------------------------------------
!  Read in sea surface height observations and error variance.
!-----------------------------------------------------------------------
!
      IF (assi_SSH(ng)) THEN
#  ifdef NUDGING_SSH
        CALL get_2dfld (ng, idSSHo, ncSSHid(ng), 1, SSHname(ng),        &
     &                  update(1), LBi, UBi, LBj, UBj, 2, 1,            &
#   ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
#   endif
     &                  OBS(ng) % SSHdat(LBi,LBj,1))
        CALL get_2dfld (ng, idSSHe, ncSSHid(ng), 1, SSHname(ng),        &
     &                  update(2), LBi, UBi, LBj, UBj, 2, 1,            &
#   ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
#   endif
     &                  OBS(ng) % EdatSSH(LBi,LBj,1))
        IF (update(1).and.update(2)) update_SSH(ng)=.true.
#  else
        CALL get_2dfld (ng, idSSHo, ncSSHid(ng), 1, SSHname(ng),        &
     &                  update(1), LBi, UBi, LBj, UBj, 1, 1,            &
#   ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
#   endif
     &                  OBS(ng) % SSHobs(LBi,LBj))
        CALL get_2dfld (ng, idSSHe, ncSSHid(ng), 1, SSHname(ng),        &
     &                  update(2), LBi, UBi, LBj, UBj, 1, 1,            &
#   ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
#   endif
     &                  OBS(ng) % EobsSSH(LBi,LBj))
        IF (update(1).and.update(2)) THEN
          update_SSH(ng)=.true.
          tSSHobs(1,ng)=Vtime(1,idSSHo,ng)
          tsSSHobs(ng)=Vtime(1,idSSHo,ng)*day2sec
          Finfo(5,idSSHo,ng)=tsSSHobs(ng)
          Finfo(5,idSSHe,ng)=tsSSHobs(ng)
          EobsSSHmin(ng)=Finfo(4,idSSHe,ng)
          EobsSSHmax(ng)=Finfo(5,idSSHe,ng)
        END IF
#  endif
      END IF
# endif
# if defined NUDGING_SST || defined ASSIMILATION_SST
!
!-----------------------------------------------------------------------
!  Read in sea surface temperature observations and error variance.
!-----------------------------------------------------------------------
!
      IF (assi_SST(ng)) THEN
#  ifdef NUDGING_SST
        CALL get_2dfld (ng, idSSTo, ncSSTid(ng), 1, SSTname(ng),        &
     &                  update(1), LBi, UBi, LBj, UBj, 2, 1,            &
#   ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
#   endif
     &                  OBS(ng) % SSTdat(LBi,LBj,1))
        CALL get_2dfld (ng, idSSTe, ncSSTid(ng), 1, SSTname(ng),        &
     &                  update(2), LBi, UBi, LBj, UBj, 2, 1,            &
#   ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
#   endif
     &                  OBS(ng) % EdatSST(LBi,LBj,1))
#  else
        CALL get_2dfld (ng, idSSTo, ncSSTid(ng), 1, SSTname(ng),        &
     &                  update(1), LBi, UBi, LBj, UBj, 1, 1,            &
#   ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
#   endif
     &                  OBS(ng) % SSTobs(LBi,LBj))
        CALL get_2dfld (ng, idSSTe, ncSSTid(ng), 1, SSTname(ng),        &
     &                  update(2), LBi, UBi, LBj, UBj, 1, 1,            &
#   ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
#   endif
     &                  OBS(ng) % EobsSST(LBi,LBj))
#  endif
        IF (update(1).and.update(2)) update_SST(ng)=.true.
      END IF
# endif
# if defined NUDGING_T || defined ASSIMILATION_T
!
!-----------------------------------------------------------------------
!  Read in tracers observations and error variance.
!-----------------------------------------------------------------------
!
      DO i=1,NAT
        IF (assi_T(i,ng)) THEN
#  ifdef NUDGING_T
          CALL get_3dfld (ng, idTobs(i), ncTOBSid(ng), 1, TOBSname(ng), &
     &                    update(1), LBi, UBi, LBj, UBj, 2, 1,          &
#   ifdef MASKING
     &                    GRID(ng) % rmask(LBi,LBj),                    &
#   endif
     &                    OBS(ng) % Tdat(LBi,LBj,1,1,i))
          CALL get_3dfld (ng, idTerr(i), ncTOBSid(ng), 1, TOBSname(ng), &
     &                    update(2), LBi, UBi, LBj, UBj, 2, 1,          &
#   ifdef MASKING
     &                    GRID(ng) % rmask(LBi,LBj),                    &
#   endif
     &                    OBS(ng) % EdatT(LBi,LBj,1,1,i))
          IF (update(1).and.update(2)) update_T(i,ng)=.true.
#  else
          CALL get_3dfld (ng, idTobs(i), ncTOBSid(ng), 1, TOBSname(ng), &
     &                    update(1), LBi, UBi, LBj, UBj, 1, 1,          &
#   ifdef MASKING
     &                    GRID(ng) % rmask(LBi,LBj),                    &
#   endif
     &                    OBS(ng) % Tobs(LBi,LBj,1,i))
          CALL get_3dfld (ng, idTerr(i), ncTOBSid(ng), 1, TOBSname(ng), &
     &                    update(2), LBi, UBi, LBj, UBj, 1, 1,          &
#   ifdef MASKING
     &                    GRID(ng) % rmask(LBi,LBj),                    &
#   endif
     &                    OBS(ng) % EobsT(LBi,LBj,1,i))
          IF (update(1).and.update(2)) THEN
            update_T(i,ng)=.true.
            tTobs(1,i,ng)=Vtime(1,idTobs(i),ng)
            tsTobs(i,ng)=Vtime(1,idTobs(i),ng)*day2sec
            Finfo(5,idTobs(i),ng)=tsTobs(i,ng)
            Finfo(5,idTerr(i),ng)=tsTobs(i,ng)
            EobsTmin(i,ng)=Finfo(4,idTerr(i),ng)
            EobsTmax(i,ng)=Finfo(5,idTerr(i),ng)
          END IF
#  endif
        END IF
      END DO
# endif
# if defined NUDGING_UVsur || defined ASSIMILATION_UVsur
!
!-----------------------------------------------------------------------
!  Read in surface current observations and error variance.
!-----------------------------------------------------------------------
!
      IF (assi_UVsur(ng)) THEN
#  ifdef NUDGING_UVsur
        CALL get_2dfld (ng, idUsur, ncVSURid(ng), 1, VSURname(ng),      &
     &                  update(1), LBi, UBi, LBj, UBj, 2, 1,            &
#   ifdef MASKING
     &                  GRID(ng) % umask(LBi,LBj),                      &
#   endif
     &                  OBS(ng) % Usurdat(LBi,LBj,1))
        CALL get_2dfld (ng, idVsur, ncVSURid(ng), 1, VSURname(ng),      &
     &                  update(2), LBi, UBi, LBj, UBj, 2, 1,            &
#   ifdef MASKING
     &                  GRID(ng) % vmask(LBi,LBj),                      &
#   endif
     &                  OBS(ng) % Vsurdat(LBi,LBj,1))
        CALL get_2dfld (ng, idUVse, ncVSURid(ng), 1, VSURname(ng),      &
     &                  update(3), LBi, UBi, LBj, UBj, 2, 1,            &
#   ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
#   endif
     &                  OBS(ng) % EdatVsur(LBi,LBj,1))
#  else
        CALL get_2dfld (ng, idUsur, ncVSURid(ng), 1, VSURname(ng),      &
     &                  update(1), LBi, UBi, LBj, UBj, 1, 1,            &
#   ifdef MASKING
     &                  GRID(ng) % umask(LBi,LBj),                      &
#   endif
     &                  OBS(ng) % Usur(LBi,LBj))
        CALL get_2dfld (ng, idVsur, ncVSURid(ng), 1, VSURname(ng),      &
     &                  update(2), LBi, UBi, LBj, UBj, 1 ,1,            &
#   ifdef MASKING
     &                  GRID(ng) % vmask(LBi,LBj),                      &
#   endif
     &                  OBS(ng) % Vsur(LBi,LBj))
        CALL get_2dfld (ng, idUVse, ncVSURid(ng), 1, VSURname(ng),      &
     &                  update(3), LBi, UBi, LBj, UBj, 1, 1,            &
#   ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
#   endif
     &                  OBS(ng) % EobsVsur(LBi,LBj))
#  endif
        IF (update(1).and.update(2).and.update(3)) update_UVsur=.true.
      END IF
# endif
# if defined NUDGING_UV || defined ASSIMILATION_UV
!
!-----------------------------------------------------------------------
!  Read in horizontal current observations and error variance.
!-----------------------------------------------------------------------
!
      IF (assi_UV(ng)) THEN
#  ifdef NUDGING_UV
        CALL get_3dfld (ng, idUobs, ncVOBSid(ng), 1, VOBSname(ng),      &
     &                  update(1), LBi, UBi, LBj, UBj, 2, 1,            &
#   ifdef MASKING
     &                  GRID(ng) % umask(LBi,LBj),                      &
#   endif
     &                  OBS(ng) % Udat(LBi,LBj,1,1))
        CALL get_3dfld (ng, idVobs, ncVOBSid(ng), 1, VOBSname(ng),      &
     &                  update(2), LBi, UBi, LBj, UBj, 2, 1,            &
#   ifdef MASKING
     &                  GRID(ng) % vmask(LBi,LBj),                      &
#   endif
     &                  OBS(ng) % Vdat(LBi,LBj,1,1))
        CALL get_3dfld (ng, idUVer, ncVOBSid(ng), 1, VOBSname(ng),      &
     &                  update(3), LBi, UBi, LBj, UBj, 2, 1,            &
#   ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
#   endif
     &                  OBS(ng) % EdatUV(LBi,LBj,1,1))
        IF (update(1).and.update(2).and.update(3)) update_UV(ng)=.true.
#  else
        CALL get_3dfld (ng, idUobs, ncVOBSid(ng), 1, VOBSname(ng),      &
     &                  update(1), LBi, UBi, LBj, UBj, 1, 1,            &
#   ifdef MASKING
     &                  GRID(ng) % umask(LBi,LBj),                      &
#   endif
     &                  OBS(ng) % Uobs(LBi,LBj,1))
        CALL get_3dfld (ng, idVobs, ncVOBSid(ng), 1, VOBSname(ng),      &
     &                  update(2), LBi, UBi, LBj, UBj, 1, 1,            &
#   ifdef MASKING
     &                  GRID(ng) % vmask(LBi,LBj),                      &
#   endif
     &                  OBS(ng) % Vobs(LBi,LBj,1))
        CALL get_3dfld (ng, idUVer, ncVOBSid(ng), 1, VOBSname(ng),      &
     &                  update(3), LBi, UBi, LBj, UBj, 1, 1,            &
#   ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
#   endif
     &                  OBS(ng) % EobsUV(LBi,LBj,1))
        IF (update(1).and.update(2).and.update(3)) THEN
          update_UV(ng)=.true.
          tVobs(1,ng)=Vtime(1,idVobs,ng)
          tsVobs(ng)=Vtime(1,idVobs,ng)*day2sec
          Finfo(5,idUobs,ng)=tsVobs(ng)
          Finfo(5,idVobs,ng)=tsVobs(ng)
          Finfo(5,idUVer,ng)=tsVobs(ng)
          EobsUVmin=Finfo(4,idUVer,ng)
          EobsUVmax=Finfo(5,idUVer,ng)
        END IF
#  endif
      END IF
# endif
#endif
#ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off input data time wall clock.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, 3)
#endif
      RETURN
      END SUBROUTINE get_data
