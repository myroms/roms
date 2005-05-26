#include "cppdefs.h"
#ifdef NONLINEAR
      SUBROUTINE set_data (ng, tile)
!
!=======================================================================
!  Copyright (c) 2005 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This subroutine processes forcing, boundary, climatology, and       !
!  assimilation input data. It time-interpolates between snapshots.    !
!                                                                      !
!=======================================================================
!
      USE mod_param

      integer, intent(in) :: ng, tile

# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 4)
# endif
      CALL set_data_tile (ng, tile, Istr, Iend, Jstr, Jend,             &
     &                    LBi, UBi, LBj, UBj)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 4)
# endif
      RETURN
      END SUBROUTINE set_data
!
!***********************************************************************
      SUBROUTINE set_data_tile (ng, tile, Istr, Iend, Jstr, Jend,       &
     &                          LBi, UBi, LBj, UBj)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
# ifdef CLIMATOLOGY
      USE mod_clima
# endif
      USE mod_forces
      USE mod_grid
      USE mod_mixing
      USE mod_ncparam
      USE mod_ocean
      USE mod_stepping
# if defined ASSIMILATION || defined NUDGING
      USE mod_obs
# endif
      USE mod_scalars
# if defined TS_PSOURCE || defined UV_PSOURCE
      USE mod_sources
# endif
!
      USE analytical_mod
# if defined NUDGING_UVsur || defined ASSIMILATION_UVsur
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      USE exchange_3d_mod, ONLY : exchange_r3d_tile
#  endif
# endif
      USE set_2dfld_mod
# ifdef SOLVE3D
      USE set_3dfld_mod
# endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, Istr, Iend, Jstr, Jend
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
!  Local variable declarations.
!
      logical :: update = .FALSE.

# ifdef OBC
      integer :: ILB, IUB, JLB, JUB
# endif
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, itrc, j, k, order

      real(r8) :: Zr, cff, cff1, cff2

      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) ::  work1
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) ::  work2

# include "set_bounds.h"

# ifdef SOLVE3D

#  ifdef CLOUDS
!
!-----------------------------------------------------------------------
!  Set cloud fraction (nondimensional).  Notice that clouds are
!  processed first in case that they are used to adjust shortwave
!  radiation.
!-----------------------------------------------------------------------
!
#   ifdef ANA_CLOUD
      CALL ana_cloud (ng, tile, iNLM)
#   else
      CALL set_2dfld_tile (ng, iNLM, Istr, Iend, Jstr, Jend, idCfra,    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%cloudG,                           &
     &                     FORCES(ng)%cloud,                            &
     &                     update)
#   endif
#  endif

#  if defined BULK_FLUXES || defined ECOSIM || \
    (defined SHORTWAVE && defined ANA_SRFLUX && defined ALBEDO)
!
!-----------------------------------------------------------------------
!  Set surface air temperature (degC).
!-----------------------------------------------------------------------
!
#   ifdef ANA_TAIR
      CALL ana_tair (ng, tile, iNLM)
#   else
      CALL set_2dfld_tile (ng, iNLM, Istr, Iend, Jstr, Jend, idTair,    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%TairG,                            &
     &                     FORCES(ng)%Tair,                             &
     &                     update)
#   endif
#  endif

#  if defined BULK_FLUXES || defined ECOSIM || \
     (defined SHORTWAVE && defined ANA_SRFLUX && defined ALBEDO)
!
!-----------------------------------------------------------------------
!  Set surface air relative or specific humidity.
!-----------------------------------------------------------------------
!
#   ifdef ANA_HUMIDITY
      CALL ana_humid (ng, tile, iNLM)
#   else
      CALL set_2dfld_tile (ng, iNLM, Istr, Iend, Jstr, Jend, idQair,    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%HairG,                            &
     &                     FORCES(ng)%Hair,                             &
     &                     update)
#   endif
#  endif

#  ifdef SHORTWAVE
!
!-----------------------------------------------------------------------
!  Set kinematic surface solar shortwave radiation flux (degC m/s).
!-----------------------------------------------------------------------
!
#   ifdef ANA_SRFLUX
      CALL ana_srflux (ng, tile, iNLM)
#   else
      CALL set_2dfld_tile (ng, iNLM, Istr, Iend, Jstr, Jend, idSrad,    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%srflxG,                           &
     &                     FORCES(ng)%srflx,                            &
     &                     update)
#   endif
#   ifdef DIURNAL_SRFLUX
!
!  Modulate the averaged shortwave radiation flux by the local diurnal
!  cycle.
!
      CALL ana_srflux (ng, tile, iNLM)
#   endif
#  endif

#  if defined BULK_FLUXES && !defined LONGWAVE && !defined LONGWAVE_OUT
!
!-----------------------------------------------------------------------
!  Surface net longwave radiation (degC m/s).
!-----------------------------------------------------------------------
!
      CALL set_2dfld_tile (ng, iNLM, Istr, Iend, Jstr, Jend, idLrad,    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%lrflxG,                           &
     &                     FORCES(ng)%lrflx,                            &
     &                     update)
#  endif

#  if defined LONGWAVE_OUT && defined BULK_FLUXES
!
!-----------------------------------------------------------------------
!  Surface downwelling longwave radiation (degC m/s).
!-----------------------------------------------------------------------
!
      CALL set_2dfld_tile (ng, iNLM, Istr, Iend, Jstr, Jend, idLdwn,    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%lrflxG,                           &
     &                     FORCES(ng)%lrflx,                            &
     &                     update)
#  endif

#  if defined BULK_FLUXES || defined ECOSIM
!
!-----------------------------------------------------------------------
!  Set surface air pressure (mb).
!-----------------------------------------------------------------------
!
#   ifdef ANA_PAIR
      CALL ana_pair (ng, tile, iNLM)
#   else
      CALL set_2dfld_tile (ng, iNLM, Istr, Iend, Jstr, Jend, idPair,    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%PairG,                            &
     &                     FORCES(ng)%Pair,                             &
     &                     update)
#   endif
!
!-----------------------------------------------------------------------
!  Set surface winds (m/s).
!-----------------------------------------------------------------------
!
#   ifdef ANA_WINDS
      CALL ana_winds (ng, tile, iNLM)
#   else
      CALL set_2dfld_tile (ng, iNLM, Istr, Iend, Jstr, Jend, idUair,    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%UwindG,                           &
     &                     FORCES(ng)%Uwind,                            &
     &                     update)
      CALL set_2dfld_tile (ng, iNLM, Istr, Iend, Jstr, Jend, idVair,    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%VwindG,                           &
     &                     FORCES(ng)%Vwind,                            &
     &                     update)
#    ifdef CURVGRID
!
!  If input point surface winds, rotate to curvilinear grid.
!
      IF (.not.Linfo(1,idUair,ng)) THEN
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            cff1=FORCES(ng)%Uwind(i,j)*                                 &
     &           COS(GRID(ng)%angler(i,j))+                             &
     &           FORCES(ng)%Vwind(i,j)*                                 &
     &           SIN(GRID(ng)%angler(i,j))
            cff2=FORCES(ng)%Vwind(i,j)*                                 &
     &           COS(GRID(ng)%angler(i,j))-                             &
     &           FORCES(ng)%Uwind(i,j)*                                 &
     &           SIN(GRID(ng)%angler(i,j))
            FORCES(ng)%Uwind(i,j)=cff1
            FORCES(ng)%Vwind(i,j)=cff2
          END DO
        END DO
      END IF
#    endif
#   endif
#  endif

#  ifdef BULK_FLUXES
!
!-----------------------------------------------------------------------
!  Set rain fall rate (kg/m2/s).
!-----------------------------------------------------------------------
!
#   ifdef ANA_RAIN
      CALL ana_rain (ng, tile, iNLM)
#   else
      CALL set_2dfld_tile (ng, iNLM, Istr, Iend, Jstr, Jend, idrain,    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%rainG,                            &
     &                     FORCES(ng)%rain,                             &
     &                     update)
#   endif
#  endif

#  ifndef BULK_FLUXES
!
!-----------------------------------------------------------------------
!  Set kinematic surface net heat flux (degC m/s).
!-----------------------------------------------------------------------
!
#   ifdef ANA_STFLUX
      CALL ana_stflux (ng, tile, iNLM, itemp)
#   else
      CALL set_2dfld_tile (ng, iNLM, Istr, Iend, Jstr, Jend,            &
     &                     idTsur(itemp),                               &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%stflxG(:,:,:,itemp),              &
     &                     FORCES(ng)%stflx (:,:,itemp),                &
     &                     update)
#   endif
#  endif

#  ifdef QCORRECTION
!
!-----------------------------------------------------------------------
!  Set sea surface temperature (SST) and heat flux sensitivity to
!  SST (dQdSST) which are used for surface heat flux correction.
!-----------------------------------------------------------------------
!
#   ifdef ANA_SST
      CALL ana_sst (ng, tile, iNLM)
#   else
      CALL set_2dfld_tile (ng, iNLM, Istr, Iend, Jstr, Jend, idSSTc,    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%sstG,                             &
     &                     FORCES(ng)%sst,                              &
     &                     update)
      CALL set_2dfld_tile (ng, iNLM, Istr, Iend, Jstr, Jend, iddQdT,    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%dqdtG,                            &
     &                     FORCES(ng)%dqdt,                             &
     &                     update)
#   endif
#  endif
!
!-----------------------------------------------------------------------
!  Set kinematic bottom net heat flux (degC m/s).
!-----------------------------------------------------------------------
!
#  ifdef ANA_BTFLUX
      CALL ana_btflux (ng, tile, iNLM, itemp)
#  else
      CALL set_2dfld_tile (ng, iNLM, Istr, Iend, Jstr, Jend,            &
     &                     idTbot(itemp),                               &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%btflxG(:,:,:,itemp),              &
     &                     FORCES(ng)%btflx (:,:,itemp),                &
     &                     update)
#  endif

#  ifdef SALINITY
!
!-----------------------------------------------------------------------
!  Set kinematic surface freshwater (E-P) flux (m/s).
!-----------------------------------------------------------------------
!
#   ifdef ANA_SSFLUX
      CALL ana_stflux (ng, tile, iNLM, isalt)
#   else
#   if !(defined EMINUSP || defined SRELAXATION)
      CALL set_2dfld_tile (ng, iNLM, Istr, Iend, Jstr, Jend,            &
     &                     idsfwf,                                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%stflxG(:,:,:,isalt),              &
     &                     FORCES(ng)%stflx (:,:,isalt),                &
     &                     update)
#    endif
#   endif

#   if defined SCORRECTION || defined SRELAXATION
!
!-----------------------------------------------------------------------
!  Set surface salinity for freshwater flux correction.
!-----------------------------------------------------------------------
!
#    ifdef ANA_SSS
      CALL ana_sss (ng, tile, iNLM)
#    else
      CALL set_2dfld_tile (ng, iNLM, Istr, Iend, Jstr, Jend, idSSSc,    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%sssG,                             &
     &                     FORCES(ng)%sss,                              &
     &                     update)
#    endif
#   endif
!
!-----------------------------------------------------------------------
!  Set kinematic bottom salt flux (m/s).
!-----------------------------------------------------------------------
!
#   ifdef ANA_BSFLUX
      CALL ana_btflux (ng, tile, iNLM, isalt)
#   else
      CALL set_2dfld_tile (ng, iNLM, Istr, Iend, Jstr, Jend,            &
     &                     idTbot(isalt),                               &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%btflxG(:,:,:,isalt),              &
     &                     FORCES(ng)%btflx (:,:,isalt),                &
     &                     update)
#   endif
#  endif

#  if defined SEDIMENT || defined BIOLOGY
!
!-----------------------------------------------------------------------
!  Set kinematic surface and bottom pasive tracer fluxes (T m/s).
!-----------------------------------------------------------------------
!
      DO itrc=NAT+1,NT(ng)
#   ifdef ANA_SPFLUX
        CALL ana_stflux (ng, tile, iNLM, itrc)
#   else
        CALL set_2dfld_tile (ng, iNLM, Istr, Iend, Jstr, Jend,          &
     &                       idTsur(itrc),                              &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       FORCES(ng)%stflxG(:,:,:,itrc),             &
     &                       FORCES(ng)%stflx (:,:,itrc),               &
     &                       update)
#   endif
#   ifdef ANA_BPFLUX
        CALL ana_btflux (ng, tile, iNLM, itrc)
#   else
        CALL set_2dfld_tile (ng, iNLM, Istr, Iend, Jstr, Jend,          &
     &                       idTbot(itrc),                              &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       FORCES(ng)%btflxG(:,:,:,itrc),             &
     &                       FORCES(ng)%btflx (:,:,itrc),               &
     &                       update)
#   endif
      END DO
#  endif
# endif

# ifndef AIR_OCEAN
#  ifndef BULK_FLUXES
!
!-----------------------------------------------------------------------
!  Set kinematic surface momentum flux (m2/s2).
!-----------------------------------------------------------------------
!
#   ifdef ANA_SMFLUX
      CALL ana_smflux (ng, tile, iNLM)
#   else
      CALL set_2dfld_tile (ng, iNLM, Istr, Iend, Jstr, Jend, idUsms,    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%sustrG,                           &
     &                     FORCES(ng)%sustr,                            &
     &                     update)
      CALL set_2dfld_tile (ng, iNLM, Istr, Iend, Jstr, Jend, idVsms,    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%svstrG,                           &
     &                     FORCES(ng)%svstr,                            &
     &                     update)
#    ifdef CURVGRID
!
!  If input point wind stress, rotate to curvilinear grid.  Notice
!  that rotation is done at RHO-points. It does not matter.
!
      IF (.not.Linfo(1,idUsms,ng)) THEN
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            cff1=FORCES(ng)%sustr(i,j)*                                 &
     &           COS(GRID(ng)%angler(i,j))+                             &
     &           FORCES(ng)%svstr(i,j)*                                 &
     &           SIN(GRID(ng)%angler(i,j))
            cff2=FORCES(ng)%svstr(i,j)*                                 &
     &           COS(GRID(ng)%angler(i,j))-                             &
     &           FORCES(ng)%sustr(i,j)*                                 &
     &           SIN(GRID(ng)%angler(i,j))
            FORCES(ng)%sustr(i,j)=cff1
            FORCES(ng)%svstr(i,j)=cff2
          END DO
        END DO
      END IF
#    endif
#   endif
#  endif
# endif

# ifdef BBL_MODEL
!
!-----------------------------------------------------------------------
!  Set surface wind-induced wave amplitude, direction and period.
!-----------------------------------------------------------------------
!
#  ifdef ANA_WWAVE
      CALL ana_wwave (ng, tile, iNLM)
#  else
#   ifdef SWAN
      CALL set_2dfld_tile (ng, iNLM, Istr, Iend, Jstr, Jend, idWorb,    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%Ub_swanG,                         &
     &                     FORCES(ng)%Ub_swan,                          &
     &                     update)
#   else
      CALL set_2dfld_tile (ng, iNLM, Istr, Iend, Jstr, Jend, idWamp,    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%AwaveG,                           &
     &                     FORCES(ng)%Awave,                            &
     &                     update)
#   endif
      CALL set_2dfld_tile (ng, iNLM, Istr, Iend, Jstr, Jend, idWdir,    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%DwaveG,                           &
     &                     FORCES(ng)%Dwave,                            &
     &                     update)
      CALL set_2dfld_tile (ng, iNLM, Istr, Iend, Jstr, Jend, idWper,    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%PwaveG,                           &
     &                     FORCES(ng)%Pwave,                            &
     &                     update)
#   ifdef CURVGRID
!
!  If input point-data, rotate direction to curvilinear coordinates.
!
      IF (.not.Linfo(1,idWdir,ng)) THEN
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            FORCES(ng)%Dwave(i,j)=FORCES(ng)%Dwave(i,j)-                &
     &                            GRID(ng)%angler(i,j)
          END DO
        END DO
      END IF
#   endif
#  endif
# endif

# if defined ECOSIM && defined SOLVE3D
!
!-----------------------------------------------------------------------
!  Compute spectral irradiance and cosine of average zenith angle of
!  downwelling spectral photons.
!-----------------------------------------------------------------------
!
      CALL ana_specir (ng, tile, iNLM)
# endif     

# ifdef ANA_SPINNING
!
!-----------------------------------------------------------------------
!  Set time-varying rotation force (centripetal accelerations) for
!  polar coordinate grids.
!-----------------------------------------------------------------------
!
      CALL ana_spinning (ng, tile, iNLM)
# endif

# if defined UV_PSOURCE || defined TS_PSOURCE
!
!-----------------------------------------------------------------------
!  Set point Sources/Sinks (river runoff).
!-----------------------------------------------------------------------
!
      IF (SOUTH_WEST_TEST) THEN
#  ifdef ANA_PSOURCE
        CALL ana_psource (ng, tile, iNLM)
#  else
#   ifdef UV_PSOURCE
        CALL set_ngfld (ng, iNLM, idRtra, 1, Msrc, 1,                   &
     &                  1, Nsrc(ng), 1,                                 &
     &                  SOURCES(ng) % QbarG,                            &
     &                  SOURCES(ng) % Qbar,                             &
     &                  update)
#    ifdef SOLVE3D
        DO k=1,N(ng)
          DO i=1,Nsrc(ng)
            SOURCES(ng)%Qsrc(i,k)=SOURCES(ng)%Qbar(i)*                  &
     &                            SOURCES(ng)%Qshape(i,k)
          END DO
        END DO
#    endif
#   endif
#   if defined TS_PSOURCE && defined SOLVE3D
        DO itrc=1,NT(ng)
          IF (SOURCES(ng)%Ltracer(itrc)) THEN
            CALL set_ngfld (ng, iNLM, idRtrc(itrc), 1, Msrc, N(ng),     &
     &                      1, Nsrc(ng), N(ng),                         &
     &                      SOURCES(ng) % TsrcG(1,1,1,itrc),            &
     &                      SOURCES(ng) % Tsrc(1,1,itrc),               &
     &                      update)
            IF (.not.update) THEN
              DO i=1,Nsrc(ng)
                SOURCES(ng) % Lsrc(i,itrc)=.FALSE.
              END DO
              exit_flag=0
            END IF
          END IF
        END DO
#   endif
#  endif
      END IF
# endif

# ifdef OBC
!
!-----------------------------------------------------------------------
!  Set open boundary conditions fields.
!-----------------------------------------------------------------------
!
!  Lower and upper bounds for nontiled arrays.
!
      ILB=LOWER_BOUND_I
      IUB=UPPER_BOUND_I
      JLB=LOWER_BOUND_J
      JUB=UPPER_BOUND_J

#  ifdef ANA_FSOBC
      CALL ana_fsobc (ng, tile, iNLM)
#  else
      IF (SOUTH_WEST_TEST) THEN
#   ifdef WEST_FSOBC
        CALL set_ngfld (ng, iNLM, idZbry(iwest), JLB, JUB, 1,           &
     &                  0, Mm(ng)+1, 1,                                 &
     &                  BOUNDARY(ng) % zetaG_west(JLB,1),               &
     &                  BOUNDARY(ng) % zeta_west(JLB),                  &
     &                  update)
#   endif
#   ifdef EAST_FSOBC
        CALL set_ngfld (ng, iNLM, idZbry(ieast), JLB, JUB, 1,           &
     &                  0, Mm(ng)+1, 1,                                 &
     &                  BOUNDARY(ng) % zetaG_east(JLB,1),               &
     &                  BOUNDARY(ng) % zeta_east(JLB),                  &
     &                  update)
#   endif
#   ifdef SOUTH_FSOBC
        CALL set_ngfld (ng, iNLM, idZbry(isouth), ILB, IUB, 1,          &
     &                  0, Lm(ng)+1 ,1,                                 &
     &                  BOUNDARY(ng) % zetaG_south(ILB,1),              &
     &                  BOUNDARY(ng) % zeta_south(ILB),                 &
     &                  update)
#   endif
#   ifdef NORTH_FSOBC
        CALL set_ngfld (ng, iNLM, idZbry(inorth), ILB, IUB, 1,          &
     &                  0, Lm(ng)+1, 1,                                 &
     &                  BOUNDARY(ng) % zetaG_north(ILB,1),              &
     &                  BOUNDARY(ng) % zeta_north(ILB),                 &
     &                  update)
#   endif
      END IF
#  endif
#  ifdef ANA_M2OBC
      CALL ana_m2obc (ng, tile, iNLM)
#  else
      IF (SOUTH_WEST_TEST) THEN
#   ifdef WEST_M2OBC
        CALL set_ngfld (ng, iNLM, idU2bc(iwest), JLB, JUB, 1,           &
     &                  0, Mm(ng)+1, 1,                                 &
     &                  BOUNDARY(ng) % ubarG_west(JLB,1),               &
     &                  BOUNDARY(ng) % ubar_west(JLB),                  &
     &                  update)
        CALL set_ngfld (ng, iNLM, idV2bc(iwest), JLB, JUB, 1,           &
     &                  1, Mm(ng)+1, 1,                                 &
     &                  BOUNDARY(ng) % vbarG_west(JLB,1),               &
     &                  BOUNDARY(ng) % vbar_west(JLB),                  &
     &                  update)
#   endif
#   ifdef EAST_M2OBC
        CALL set_ngfld (ng, iNLM, idU2bc(ieast), JLB, JUB, 1,           &
     &                  0, Mm(ng)+1, 1,                                 &
     &                  BOUNDARY(ng) % ubarG_east(JLB,1),               &
     &                  BOUNDARY(ng) % ubar_east(JLB),                  &
     &                  update)
        CALL set_ngfld (ng, iNLM, idV2bc(ieast), JLB, JUB, 1,           &
     &                  1, Mm(ng)+1, 1,                                 &
     &                  BOUNDARY(ng) % vbarG_east(JLB,1),               &
     &                  BOUNDARY(ng) % vbar_east(JLB),                  &
     &                  update)
#   endif
#   ifdef SOUTH_M2OBC
        CALL set_ngfld (ng, iNLM, idU2bc(isouth), ILB, IUB, 1,          &
     &                  1, Lm(ng)+1, 1,                                 &
     &                  BOUNDARY(ng) % ubarG_south(ILB,1),              &
     &                  BOUNDARY(ng) % ubar_south(ILB),                 &
     &                  update)
        CALL set_ngfld (ng, iNLM, idV2bc(isouth), ILB, IUB, 1,          &
     &                  0, Lm(ng)+1, 1,                                 &
     &                  BOUNDARY(ng) % vbarG_south(ILB,1),              &
     &                  BOUNDARY(ng) % vbar_south(ILB),                 &
     &                  update)
#   endif
#   ifdef NORTH_M2OBC
        CALL set_ngfld (ng, iNLM, idU2bc(inorth), ILB, IUB, 1,          &
     &                  1, Lm(ng)+1, 1,                                 &
     &                  BOUNDARY(ng) % ubarG_north(ILB,1),              &
     &                  BOUNDARY(ng) % ubar_north(ILB),                 &
     &                  update)
        CALL set_ngfld (ng, iNLM, idV2bc(inorth), ILB, IUB, 1,          &
     &                  0, Lm(ng)+1, 1,                                 &
     &                  BOUNDARY(ng) % vbarG_north(ILB,1),              &
     &                  BOUNDARY(ng) % vbar_north(ILB),                 &
     &                  update)
#   endif
      END IF
#  endif
#  ifdef SOLVE3D
#   ifdef ANA_M3OBC
      CALL ana_m3obc (ng, tile, iNLM)
#   else
      IF (SOUTH_WEST_TEST) THEN
#    ifdef WEST_M3OBC
        CALL set_ngfld (ng, iNLM, idU3bc(iwest), JLB, JUB, N(ng),       &
     &                  0, Mm(ng)+1, N(ng),                             &
     &                  BOUNDARY(ng) % uG_west(JLB,1,1),                &
     &                  BOUNDARY(ng) % u_west(JLB,1),                   &
     &                  update)
        CALL set_ngfld (ng, iNLM, idV3bc(iwest), JLB, JUB, N(ng),       &
     &                  1, Mm(ng)+1, N(ng),                             &
     &                  BOUNDARY(ng) % vG_west(JLB,1,1),                &
     &                  BOUNDARY(ng) % v_west(JLB,1),                   &
     &                  update)
#    endif
#    ifdef EAST_M3OBC
        CALL set_ngfld (ng, iNLM, idU3bc(ieast), JLB, JUB, N(ng),       &
     &                  0, Mm(ng)+1, N(ng),                             &
     &                  BOUNDARY(ng) % uG_east(JLB,1,1),                &
     &                  BOUNDARY(ng) % u_east(JLB,1),                   &
     &                  update)
        CALL set_ngfld (ng, iNLM, idV3bc(ieast), JLB, JUB, N(ng),       &
     &                  1, Mm(ng)+1, N(ng),                             &
     &                  BOUNDARY(ng) % vG_east(JLB,1,1),                &
     &                  BOUNDARY(ng) % v_east(JLB,1),                   &
     &                  update)
#    endif
#    ifdef SOUTH_M3OBC
        CALL set_ngfld (ng, iNLM, idU3bc(isouth), ILB, IUB, N(ng),      &
     &                  1, Lm(ng)+1, N(ng),                             &
     &                  BOUNDARY(ng) % uG_south(ILB,1,1),               &
     &                  BOUNDARY(ng) % u_south(ILB,1),                  &
     &                  update)
        CALL set_ngfld (ng, iNLM, idV3bc(isouth), ILB, IUB, N(ng),      &
     &                  0, Lm(ng)+1, N(ng),                             &
     &                  BOUNDARY(ng) % vG_south(ILB,1,1),               &
     &                  BOUNDARY(ng) % v_south(ILB,1),                  &
     &                  update)
#    endif
#    ifdef NORTH_M3OBC
        CALL set_ngfld (ng, iNLM, idU3bc(inorth), ILB, IUB, N(ng),      &
     &                  1, Lm(ng)+1, N(ng),                             &
     &                  BOUNDARY(ng) % uG_north(ILB,1,1),               &
     &                  BOUNDARY(ng) % u_north(ILB,1),                  &
     &                  update)
        CALL set_ngfld (ng, iNLM, idV3bc(inorth), ILB, IUB, N(ng),      &
     &                  0, Lm(ng)+1, N(ng),                             &
     &                  BOUNDARY(ng) % vG_north(ILB,1,1),               &
     &                  BOUNDARY(ng) % v_north(ILB,1),                  &
     &                  update)
#    endif
      END IF
#   endif
#   ifdef ANA_TOBC
      CALL ana_tobc (ng, tile, iNLM)
#   else
      IF (SOUTH_WEST_TEST) THEN
        DO itrc=1,NT(ng)
#    ifdef WEST_TOBC
          CALL set_ngfld (ng, iNLM, idTbry(iwest,itrc),                 &
     &                    JLB, JUB, N(ng), 0, Mm(ng)+1, N(ng),          &
     &                    BOUNDARY(ng) % tG_west(JLB,1,1,itrc),         &
     &                    BOUNDARY(ng) % t_west(JLB,1,itrc),            &
     &                    update)
#    endif
#    ifdef EAST_TOBC
          CALL set_ngfld (ng, iNLM, idTbry(ieast,itrc),                 &
     &                    JLB, JUB, N(ng), 0, Mm(ng)+1, N(ng),          &
     &                    BOUNDARY(ng) % tG_east(JLB,1,1,itrc),         &
     &                    BOUNDARY(ng) % t_east(JLB,1,itrc),            &
     &                    update)
#    endif
#    ifdef SOUTH_TOBC
          CALL set_ngfld (ng, iNLM, idTbry(isouth,itrc),                &
     &                    ILB, IUB, N(ng), 0, Lm(ng)+1, N(ng),          &
     &                    BOUNDARY(ng) % tG_south(ILB,1,1,itrc),        &
     &                    BOUNDARY(ng) % t_south(ILB,1,itrc),           &
     &                    update)
#    endif
#    ifdef NORTH_TOBC
          CALL set_ngfld (ng, iNLM, idTbry(inorth,itrc),                &
     &                    ILB, IUB, N(ng), 0, Lm(ng)+1, N(ng),          &
     &                    BOUNDARY(ng) % tG_north(ILB,1,1,itrc),        &
     &                    BOUNDARY(ng) % t_north(ILB,1,itrc),           &
     &                    update)
#    endif
        END DO
      END IF
#   endif
#  endif
# endif

# ifdef ZCLIMATOLOGY
!
!-----------------------------------------------------------------------
!  Set sea surface height climatology (m).
!-----------------------------------------------------------------------
!
#  ifdef ANA_SSH
      CALL ana_ssh (ng, tile, iNLM)
#  else
      CALL set_2dfld_tile (ng, iNLM, Istr, Iend, Jstr, Jend, idSSHc,    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     CLIMA(ng)%sshG,                              &
     &                     CLIMA(ng)%ssh,                               &
     &                     update)
#  endif
# endif

# ifdef M2CLIMATOLOGY
!
!-----------------------------------------------------------------------
!  Set 2D momentum climatology (m/s).
!-----------------------------------------------------------------------
!
#  ifdef ANA_M2CLIMA
      CALL ana_m2clima (ng, tile, iNLM)
#  else
      CALL set_2dfld_tile (ng, iNLM, Istr, Iend, Jstr, Jend, idUbcl,    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     CLIMA(ng)%ubarclmG,                          &
     &                     CLIMA(ng)%ubarclm,                           &
     &                     update)
      CALL set_2dfld_tile (ng, iNLM, Istr, Iend, Jstr, Jend, idVbcl,    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     CLIMA(ng)%vbarclmG,                          &
     &                     CLIMA(ng)%vbarclm,                           &
     &                     update)
#  endif
# endif

# if defined SOLVE3D && defined TCLIMATOLOGY
!
!-----------------------------------------------------------------------
!  Set tracer climatology.
!-----------------------------------------------------------------------
!
#  ifdef ANA_TCLIMA
      CALL ana_tclima (ng, tile, iNLM)
#  else
      DO itrc=1,NAT
        CALL set_3dfld_tile (ng, iNLM, Istr, Iend, Jstr, Jend,          &
     &                       idTclm(itrc),                              &
     &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                       CLIMA(ng)%tclmG(:,:,:,:,itrc),             &
     &                       CLIMA(ng)%tclm (:,:,:,itrc),               &
     &                       update)
      END DO
#  endif
# endif

# if defined SOLVE3D && defined M3CLIMATOLOGY
!
!-----------------------------------------------------------------------
!  Set 3D momentum climatology (m/s).
!-----------------------------------------------------------------------
!
#  ifdef ANA_M3CLIMA
      CALL ana_m3clima (ng, tile, iNLM)
#  else
      CALL set_3dfld_tile (ng, iNLM, Istr, Iend, Jstr, Jend, idUclm,    &
     &                     LBi, UBi, LBj, UBj, 1, N(ng),                &
     &                     CLIMA(ng)%uclmG,                             &
     &                     CLIMA(ng)%uclm,                              &
     &                     update)
      CALL set_3dfld_tile (ng, iNLM, Istr, Iend, Jstr, Jend, idVclm,    &
     &                     LBi, UBi, LBj, UBj, 1, N(ng),                &
     &                     CLIMA(ng)%vclmG,                             &
     &                     CLIMA(ng)%vclm,                              &
     &                     update)
#  endif
# endif

# if defined NUDGING_SSH
!
!-----------------------------------------------------------------------
!  Set sea surface height observations and error variance.
!-----------------------------------------------------------------------
!
      IF (assi_SSH(ng)) THEN
        CALL set_2dfld_tile (ng, iNLM, Istr, Iend, Jstr, Jend, idSSHo,  &
     &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                       OBS(ng)%SSHdat,                            &
     &                       OBS(ng)%SSHobs,                            &
     &                       update)
        CALL set_2dfld_tile (ng, iNLM, Istr, Iend, Jstr, Jend, idSSHe,  &
     &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                       OBS(ng)%EdatSSH,                           &
     &                       OBS(ng)%EobsSSH,                           &
     &                       update)
        IF (.not.update.and.SOUTH_WEST_TEST) THEN
          update_SSH(ng)=.FALSE.
          exit_flag=0
        END IF
      END IF
# endif

# ifdef SOLVE3D
#  if defined NUDGING_SST || defined ASSIMILATION_SST
!
!-----------------------------------------------------------------------
!  Set sea surface temperature observations and error variance.
!-----------------------------------------------------------------------
!
      IF (assi_SST(ng)) THEN
#   ifdef NUDGING_SST
        CALL set_2dfld_tile (ng, iNLM, Istr, Iend, Jstr, Jend, idSSTo,  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       OBS(ng)%SSTdat,                            &
     &                       OBS(ng)%SSTobs,                            &
     &                       update)
        CALL set_2dfld_tile (ng, iNLM, Istr, Iend, Jstr, Jend, idSSTe,  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       OBS(ng)%EdatSST,                           &
     &                       OBS(ng)%EobsSST,                           &
     &                       update)
        IF (.not.update.and.SOUTH_WEST_TEST) THEN
          update_SST(ng)=.FALSE.
          update_T(itemp,ng)=.FALSE.
          exit_flag=0
        END IF
#   endif
!
!  Extend sea surface temperature and associated error variance using
!  provided basis function polynomials.
!
        IF (extend_SST(ng).and.update_SST(ng)) THEN
          IF (SOUTH_WEST_TEST) THEN
            update_SST(ng)=.FALSE.
            update_T(itemp,ng)=.TRUE.
#   ifdef ASSIMILATION_SST
            tTobs(1,itemp,ng)=Vtime(1,idSSTo,ng)
            tsTobs(itemp,ng)=Vtime(1,idSSTo,ng)*day2sec
            Finfo(7,idSSTo,ng)=tsTobs(itemp,ng)
            Finfo(7,idSSTe,ng)=tsTobs(itemp,ng)
            EobsTmin(itemp,ng)=Finfo(8,idSSTe,ng)
            EobsTmax(itemp,ng)=Finfo(9,idSSTe,ng)
#   endif
          END IF
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                Zr=GRID(ng)%z_r(i,j,k)/GRID(ng)%h(i,j)
                cff=perr_SST(npSST(ng),ng)
                cff1=pcoef_SST(npSST(ng),ng)
                DO order=npSST(ng)-1,0,-1
                  cff=Zr*cff+perr_SST(order,ng)
                  cff1=Zr*cff+pcoef_SST(order,ng)
                END DO
                OBS(ng)%EobsT(i,j,k,itemp)=MIN(1.0_r8,cff*              &
     &                                         OBS(ng)%EobsSST(i,j))
                OBS(ng)%Tobs(i,j,k,itemp)=cff1*OBS(ng)%SSTobs(i,j)
              END DO
            END DO
          END DO
        END IF
      END IF
#  endif

#  if defined NUDGING_T
!
!-----------------------------------------------------------------------
!  Set tracers observations and error variance.
!-----------------------------------------------------------------------
!
      DO itrc=1,NAT
        IF (assi_T(itrc,ng)) THEN
          CALL set_3dfld_tile (ng, iNLM, Istr, Iend, Jstr, Jend,        &
     &                         idTobs(itrc),                            &
     &                         LBi, UBi, LBj, UBj, 1, N(ng),            &
     &                         OBS(ng)%Tdat(:,:,:,:,itrc),              &
     &                         OBS(ng)%Tobs(:,:,:,itrc),                &
     &                         update)
          CALL set_3dfld_tile (ng, iNLM, Istr, Iend, Jstr, Jend,        &
     &                         idTerr(itrc),                            &
     &                         LBi, UBi, LBj, UBj, 1, N(ng),            &
     &                         OBS(ng)%EdatT(:,:,:,:,itrc),             &
     &                         OBS(ng)%EobsT(:,:,:,itrc),               &
     &                         update)
          IF (.not.update.and.SOUTH_WEST_TEST) THEN
            update_T(itrc,ng)=.FALSE.
            exit_flag=0
          END IF
        END IF
      END DO
#  endif

#  if defined NUDGING_UVsur || defined ASSIMILATION_UVsur
!
!-----------------------------------------------------------------------
!  Set surface current observations and error variance.
!-----------------------------------------------------------------------
!
      IF (assi_UVsur(ng)) THEN
#   ifdef NUDGING_UVsur
        CALL set_2dfld_tile (ng, iNLM, Istr, Iend, Jstr, Jend, idUsur,  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       OBS(ng)%Usurdat,                           &
     &                       OBS(ng)%Usur,                              &
     &                       update)
        CALL set_2dfld_tile (ng, iNLM, Istr, Iend, Jstr, Jend, idVsur,  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       OBS(ng)%Vsurdat,                           &
     &                       OBS(ng)%Vsur,                              &
     &                       update)
        CALL set_2dfld_tile (ng, iNLM, Istr, Iend, Jstr, Jend, idUVse,  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       OBS(ng)%EdatVsur,                          &
     &                       OBS(ng)%EobsVsur,                          &
     &                       update)
        IF (.not.update.and.SOUTH_WEST_TEST) THEN
          update_UVsur(ng)=.FALSE.
          update_UV(ng)=.FALSE.
          exit_flag=0
        END IF
#   endif
!
!  Extend surface currents observations and associated error variance
!  using provided basis function polynomials.
!
        IF (extend_UV(ng).and.update_UVsur(ng)) THEN
          IF (SOUTH_WEST_TEST) THEN
            update_UVsur(ng)=.FALSE.
            update_UV(ng)=.TRUE.
#  ifdef ASSIMILATION_UVsur
            tVobs(1,ng)=Vtime(1,idVsur,ng)
            tsVobs(ng)=Vtime(1,idVsur,ng)*day2sec
            Finfo(7,idUsur,ng)=tsVobs(ng)
            Finfo(7,idVsur,ng)=tsVobs(ng)
            Finfo(7,idUVse,ng)=tsVobs(ng)
            EobsUVmin(ng)=Finfo(8,idUVse,ng)
            EobsUVmax(ng)=Finfo(9,idUVse,ng)
#  endif
          END IF
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                Zr=GRID(ng)%z_r(i,j,k)/GRID(ng)%h(i,j)
                cff=perr_V(npUV(ng),ng)
                DO order=npUV(ng)-1,0,-1
                  cff=Zr*cff+perr_V(order,ng)
                END DO
                OBS(ng)%EobsUV(i,j,k)=MIN(1.0_r8,cff+                   &
     &                                    OBS(ng)%EobsVsur(i,j))
              END DO
            END DO
            DO j=JstrV-1,Jend
              DO i=IstrU-1,Iend
                Zr=GRID(ng)%z_r(i,j,k)/GRID(ng)%h(i,j)
                cff1=pcoef_U(npUV(ng),ng)
                cff2=pcoef_V(npUV(ng),ng)
                DO order=npUV(ng)-1,0,-1
                  cff1=Zr*cff1+pcoef_U(order,ng)
                  cff2=Zr*cff2+pcoef_V(order,ng)
                END DO
                work1(i,j)=cff1*OBS(ng)%Usur(i,j)-                      &
     &                     cff2*OBS(ng)%Vsur(i,j)
                work2(i,j)=cff2*OBS(ng)%Usur(i,j)+                      &
     &                     cff1*OBS(ng)%Vsur(i,j)
              END DO
            END DO
            DO j=Jstr,Jend
              DO i=IstrU,Iend
                OBS(ng)%Uobs(i,j,k)=0.5_r8*(work1(i-1,j)+work1(i,j))
              END DO
              IF (j.ge.JstrV) THEN
                DO i=Istr,Iend
                  OBS(ng)%Vobs(i,j,k)=0.5_r8*(work2(i,j-1)+work2(i,j))
                END DO
              END IF
            END DO
          END DO
#   if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
          CALL exchange_r3d_tile (ng, iNLM, Istr, Iend, Jstr, Jend,     &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            NghostPoints,                         &
     &                            OBS(ng)%EobsUV)
#  endif
        END IF
      END IF
#  endif

#  ifdef NUDGING_UV
!
!-----------------------------------------------------------------------
!  Set horizontal current observations and error variance.
!-----------------------------------------------------------------------
!
      IF (assi_UV(ng)) THEN
        CALL set_3dfld_tile (ng, iNLM, Istr, Iend, Jstr, Jend, idUobs,  &
     &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                       OBS(ng)%Udat,                              &
     &                       OBS(ng)%Uobs,                              &
     &                       update)
        CALL set_3dfld_tile (ng, iNLM, Istr, Iend, Jstr, Jend, idVobs,  &
     &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                       OBS(ng)%Vdat,                              &
     &                       OBS(ng)%Vobs,                              &
     &                       update)
        CALL set_3dfld_tile (ng, iNLM, Istr, Iend, Jstr, Jend, idUVer,  &
     &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                       OBS(ng)%EdatUV,                            &
     &                       OBS(ng)%EobsUV,                            &
     &                       update)
        IF (.not.update.and.SOUTH_WEST_TEST) THEN
          update_UV(ng)=.FALSE.
          exit_flag=0
        END IF
      END IF
#  endif
# endif
      RETURN
      END SUBROUTINE set_data_tile
#else
      SUBROUTINE set_data
      RETURN
      END SUBROUTINE set_data
#endif
