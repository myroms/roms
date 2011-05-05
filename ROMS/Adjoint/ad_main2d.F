#include "cppdefs.h"
#if defined ADJOINT && !defined SOLVE3D
      SUBROUTINE ad_main2d (RunInterval)
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2011 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine is the main driver for adjoint  ROMS/TOMS  when     !
!  configurated as shallow water (barotropic) ocean model only. It     !
!  advances backward  the  adjoint model  for all nested grids, if     !
!  any, by the specified time interval (seconds), RunInterval.         !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
# ifdef MODEL_COUPLING
      USE mod_coupler
# endif
# ifdef FOUR_DVAR
      USE mod_fourdvar
# endif
      USE mod_iounits
      USE mod_scalars
      USE mod_stepping
# ifdef SO_SEMI
      USE mod_storage
# endif
!
# if defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
     defined OPT_OBSERVATIONS || defined SENSITIVITY_4DVAR
      USE adsen_force_mod, ONLY : adsen_force
# endif
      USE ad_diag_mod, ONLY : ad_diag
# ifdef WEAK_CONSTRAINT
      USE ad_forcing_mod, ONLY : ad_forcing
# endif
# ifdef ADJUST_WSTRESS
      USE ad_frc_adjust_mod, ONLY : ad_frc_adjust
# endif
      USE ad_ini_fields_mod, ONLY : ad_ini_fields, ad_ini_zeta
# if defined FOUR_DVAR && defined OBSERVATIONS
#  ifdef WEAK_CONSTRAINT
      USE ad_htobs_mod, ONLY : ad_htobs
#  else
      USE ad_misfit_mod, ONLY : ad_misfit
#  endif
# endif
# ifdef ADJUST_BOUNDARY
      USE ad_obc_adjust_mod, ONLY : ad_obc_adjust
# endif
# ifdef NEARSHORE_MELLOR_NOT_YET
!!    USE ad_radiation_stress_mod, ONLY : ad_radiation_stress
# endif
# ifdef AD_AVERAGES
      USE ad_set_avg_mod, ONLY : ad_set_avg
# endif
# if defined SSH_TIDES_NOT_YET || defined UV_TIDES_NOT_YET
!!    USE ad_set_tides_mod, ONLY : ad_set_tides
# endif
      USE ad_set_vbc_mod, ONLY: ad_set_vbc
      USE ad_step2d_mod, ONLY : ad_step2d
# ifdef FLOATS_NOT_YET
!!    USE ad_step_floats_mod, ONLY : tl_step_floats
# endif
# ifdef WEAK_CONSTRAINT
      USE frc_weak_mod, ONLY : frc_ADgather, frc_clear
# endif
# ifdef AIR_OCEAN_NOT_YET
      USE ocean_coupler_mod, ONLY : atmos_coupling
# endif
# ifdef WAVES_OCEAN_NOT_YET
      USE ocean_coupler_mod, ONLY : waves_coupling
# endif
# ifdef SO_SEMI
      USE packing_mod, ONLY : ad_pack
# endif
!
      implicit none
!
!  Imported variable declarations.
!
      real(r8), intent(in) :: RunInterval
!
!  Local variable declarations.
!
      logical:: backward = .TRUE.

      integer :: next_indx1, ng, subs, tile, thread
# ifdef FLOATS_NOT_YET
      integer :: Lend, Lstr, chunk_size
# endif
      integer :: ksav, ktmp
      real(r8) :: my_StepTime
!
!=======================================================================
!  Time-step adjoint vertically integrated equations.
!=======================================================================
!
      my_StepTime=0.0_r8

      STEP_LOOP : DO WHILE (my_StepTime.le.RunInterval)

        my_StepTime=my_StepTime+MAXVAL(dt)
!
!  Set time clock.
!
        DO ng=1,Ngrids
          iic(ng)=iic(ng)-1
          time(ng)=time(ng)-dt(ng)
          tdays(ng)=time(ng)*sec2day
          CALL time_string (time(ng), time_code(ng))
        END DO
!
!-----------------------------------------------------------------------
!  Read in required data, if any, from input NetCDF files.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          CALL ad_get_data (ng)
          IF (exit_flag.ne.NoError) RETURN
        END DO
!
!-----------------------------------------------------------------------
!  Process input data, if any: time interpolate between snapshots.
!  If appropriate, compute and report diagnostics.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(numthreads)
          DO thread=0,numthreads-1
            subs=NtileX(ng)*NtileE(ng)/numthreads
            DO tile=subs*thread,subs*(thread+1)-1,+1      ! irreversible
              CALL ad_set_data (ng, TILE)
# ifdef AD_AVERAGES
              CALL ad_set_avg (ng, TILE)
# endif
# ifdef DIAGNOSTICS
!!            CALL ad_set_diags (ng, TILE)
# endif
              CALL ad_diag (ng, TILE)
            END DO
          END DO
!$OMP END PARALLEL DO
        END DO
        IF (exit_flag.ne.NoError) RETURN

# if (defined FOUR_DVAR    && !defined IS4DVAR_SENSITIVITY) && \
      defined OBSERVATIONS
!
!-----------------------------------------------------------------------
!  If appropriate, read observation and model state at observation
!  locations.  Then, compute adjoint misfit forcing terms.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
#  ifdef SENSITIVITY_4DVAR
          IF (.not.LsenPSAS(ng)) THEN
#  endif
            HalfDT=0.5_r8*dt(ng)
            IF (((time(ng)-HalfDT).le.ObsTime(ng)).and.                 &
     &          (ObsTime(ng).lt.(time(ng)+HalfDT))) THEN
              ProcessObs(ng)=.TRUE.
              CALL obs_read (ng, iADM, backward)
              IF (exit_flag.ne.NoError) RETURN
            ELSE
              ProcessObs(ng)=.FALSE.
            END IF
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(numthreads)
            DO thread=0,numthreads-1
              subs=NtileX(ng)*NtileE(ng)/numthreads
              DO tile=subs*thread,subs*(thread+1)-1,+1
#  ifdef WEAK_CONSTRAINT
                CALL ad_htobs (ng, TILE, iADM)
#  else
                CALL ad_misfit (ng, TILE, iADM)
#  endif
              END DO
            END DO
!$OMP END PARALLEL DO
            IF (exit_flag.ne.NoError) RETURN
#  ifdef SENSITIVITY_4DVAR
          END IF
#  endif
        END DO
# endif

# ifdef WEAK_CONSTRAINT
!
!-----------------------------------------------------------------------
!  If appropriate, add representer coefficients (Beta hat) impulse
!  forcing to adjoint solution. Read next impulse record, if available.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          IF (ProcessObs(ng)) THEN
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(numthreads)
            DO thread=0,numthreads-1
              subs=NtileX(ng)*NtileE(ng)/numthreads
              DO tile=subs*thread,subs*(thread+1)-1,+1
                CALL ad_forcing (ng, TILE, knew(ng), nnew(ng))
              END DO
            END DO
!$OMP END PARALLEL DO
          END IF
        END DO
# endif
!
!  Avoid time-stepping if additional delayed IO time-step.
!
        TIME_STEP: IF (MINVAL(iic).ne.MINVAL(ntstart)) THEN

# ifdef FLOATS_NOT_YET
!
!-----------------------------------------------------------------------
!  Compute Lagrangian drifters trajectories.
!-----------------------------------------------------------------------
!
!  Shift floats time indices.
!
          DO ng=1,Ngrids
            IF (Lfloats(Ng)) THEN
              nfp1(ng)=MOD(nfp1(ng)+1,NFT+1)
              nf  (ng)=MOD(nf  (ng)+1,NFT+1)
              nfm1(ng)=MOD(nfm1(ng)+1,NFT+1)
              nfm2(ng)=MOD(nfm2(ng)+1,NFT+1)
              nfm3(ng)=MOD(nfm3(ng)+1,NFT+1)
!
!$OMP PARALLEL DO PRIVATE(thread,chunk_size,Lstr,Lend)                  &
!$OMP&            SHARED(numthreads,Nfloats)
              DO thread=0,numthreads-1
                chunk_size=(Nfloats(ng)+numthreads-1)/numthreads
                Lstr=1+thread*chunk_size
                Lend=MIN(Nfloats(ng),Lstr+chunk_size-1)
                CALL ad_step_floats (ng, Lstr, Lend)
              END DO
!$OMP END PARALLEL DO
            END IF
          END DO
# endif
!
!-----------------------------------------------------------------------
!  Solve the vertically integrated primitive equations for the
!  free-surface and momentum components.
!-----------------------------------------------------------------------
!
!  Corrector step - Apply 2D time-step corrector scheme.  Notice that
!  ==============   there is not need for a corrector step during the
!  auxiliary (nfast+1) time-step.
!
          DO ng=1,Ngrids
            iif(ng)=1
            nfast(ng)=1
            IF (iif(ng).lt.(nfast(ng)+1)) THEN
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(numthreads)
              DO thread=0,numthreads-1
                subs=NtileX(ng)*NtileE(ng)/numthreads
                DO tile=subs*thread,subs*(thread+1)-1,+1
                  CALL ad_step2d (ng, TILE)
                END DO
              END DO
!$OMP END PARALLEL DO
            END IF
!
!  Set time indices for corrector step.
!
            next_indx1=3-indx1(ng)
            IF (.not.PREDICTOR_2D_STEP(ng)) THEN
              PREDICTOR_2D_STEP(ng)=.TRUE.
              ktmp=knew(ng)
              ksav=kstp(ng)
              knew(ng)=krhs(ng)
              kstp(ng)=ktmp
              krhs(ng)=ksav
            END IF
          END DO
!
!  Predictor step - Advance barotropic equations using 2D time-step
!  ==============   predictor scheme.
!
          DO ng=1,Ngrids
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(numthreads)
            DO thread=0,numthreads-1
              subs=NtileX(ng)*NtileE(ng)/numthreads
              DO tile=subs*(thread+1)-1,subs*thread,-1
                CALL ad_step2d (ng, TILE)
              END DO
            END DO
!$OMP END PARALLEL DO
!
!  Set time indices for predictor step. The PREDICTOR_2D_STEP switch
!  it is assumed to be false before the first time-step.
!
            IF (PREDICTOR_2D_STEP(ng).and.                              &
     &          (iic(ng).ne.ntend(ng))) THEN
              PREDICTOR_2D_STEP(ng)=.FALSE.
              ksav=knew(ng)
              knew(ng)=krhs(ng)
              krhs(ng)=ksav
            END IF
          END DO

        END IF TIME_STEP

# ifdef SO_SEMI
!
!-----------------------------------------------------------------------
!  If stochastic optimals with respect the seminorm of chosen
!  functional, pack adjoint state surface forcing needed by the
!  dynamical propagator.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          IF (MOD(iic(ng)-1,nADJ(ng)).eq.0) THEN
            SOrec(ng)=SOrec(ng)+1
!$OMP PARALLEL DO PRIVATE(thread,subs,tile,SOrec)                       &
!$OMP&            SHARED(numthreads,Nstr,Nend)
            DO thread=0,numthreads-1
              subs=NtileX(ng)*NtileE(ng)/numthreads
              DO tile=subs*thread,subs*(thread+1)-1,+1
                CALL ad_pack (ng, TILE, Nstr(ng), Nend(ng),             &
     &                        STORAGE(ng)%so_state(:,SOrec(ng)))
              END DO
            END DO
!$OMP END PARALLEL DO
          END IF
        END DO
# endif
!
!-----------------------------------------------------------------------
!  Set vertical boundary conditions. Process tidal forcing.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(numthreads)
          DO thread=0,numthreads-1
            subs=NtileX(ng)*NtileE(ng)/numthreads
            DO tile=subs*thread,subs*(thread+1)-1,+1
# if defined SSH_TIDES_NOT_YET || defined UV_TIDES_NOT_YET
              CALL ad_set_tides (ng, TILE)
# endif
              CALL ad_set_vbc (ng, TILE)
            END DO
          END DO
!$OMP END PARALLEL DO
        END DO

# ifdef NEARSHORE_MELLOR_NOT_YET
!
!-----------------------------------------------------------------------
!  Compute radiation stress terms.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(numthreads)
          DO thread=0,numthreads-1
            subs=NtileX(ng)*NtileE(ng)/numthreads
            DO tile=subs*(thread+1)-1,subs*thread,-1
              CALL ad_radiation_stress (ng, TILE)
            END DO
          END DO
!$OMP END PARALLEL DO
        END DO
# endif

# ifdef WAVES_OCEAN_NOT_YET
!
!-----------------------------------------------------------------------
!  Couple to waves model every CoupleSteps(Iwaves,ng) timesteps: get
!  waves/sea fluxes.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          IF ((iic(ng).ne.ntstart(ng)).and.                             &
     &        MOD(iic(ng)-1,CoupleSteps(Iwaves,ng)).eq.0) THEN
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(numthreads)
            DO thread=0,numthreads-1
              subs=NtileX(ng)*NtileE(ng)/numthreads
              DO tile=subs*thread,subs*(thread+1)-1,+1
                CALL waves_coupling (ng, TILE)
              END DO
            END DO
!$OMP END PARALLEL DO
          END IF
        END DO
# endif

# ifdef AIR_OCEAN_NOT_YET
!
!-----------------------------------------------------------------------
!  Couple to atmospheric model every CoupleSteps(Iatmos) timesteps: get
!  air/sea fluxes.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          IF ((iic(ng).ne.ntstart(ng)).and.                             &
     &        MOD(iic(ng)-1,CoupleSteps(Iatmos,ng)).eq.0) THEN
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(numthreads)
            DO thread=0,numthreads-1
              subs=NtileX(ng)*NtileE(ng)/numthreads
              DO tile=subs*(thread+1)-1,subs*thread,-1
                CALL atmos_coupling (ng, TILE)
              END DO
            END DO
!$OMP END PARALLEL DO
          END IF
        END DO
# endif
!
!-----------------------------------------------------------------------
!  If not a restart, initialize all time levels and compute other
!  initial fields.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          IF (iic(ng).eq.ntend(ng)) THEN
!
!  Initialize other state variables.
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(numthreads)
            DO thread=0,numthreads-1
              subs=NtileX(ng)*NtileE(ng)/numthreads
              DO tile=subs*(thread+1)-1,subs*thread,-1
                CALL ad_ini_fields (ng, TILE, iADM)
              END DO
            END DO
!$OMP END PARALLEL DO
!
!  Initialize free-surface.
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(numthreads)
            DO thread=0,numthreads-1
              subs=NtileX(ng)*NtileE(ng)/numthreads
              DO tile=subs*thread,subs*(thread+1)-1,+1
                CALL ad_ini_zeta (ng, TILE, iADM)
              END DO
            END DO
!$OMP END PARALLEL DO
          END IF
        END DO

# ifdef ADJUST_WSTRESS
!
!-----------------------------------------------------------------------
!  Interpolate surface forcing increments and adjust surface forcing.
!  Skip first timestep.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          IF (iic(ng).ne.ntstart(ng)) THEN
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(numthreads)
            DO thread=0,numthreads-1
              subs=NtileX(ng)*NtileE(ng)/numthreads
              DO tile=subs*thread,subs*(thread+1)-1,+1
                CALL ad_frc_adjust (ng, TILE, Lfout(ng))
              END DO
            END DO
!$OMP END PARALLEL DO
          END IF
        END DO
# endif

# ifdef ADJUST_BOUNDARY
!
!-----------------------------------------------------------------------
!  Interpolate open boundary increments and adjust open boundaries.
!  Skip first timestep.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          IF (iic(ng).ne.ntstart(ng)) THEN
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(numthreads)
            DO thread=0,numthreads-1
              subs=NtileX(ng)*NtileE(ng)/numthreads
              DO tile=subs*thread,subs*(thread+1)-1,+1
                CALL ad_obc_adjust (ng, TILE, Lbout(ng))
              END DO
            END DO
!$OMP END PARALLEL DO
          END IF
        END DO
# endif

# ifdef WEAK_CONSTRAINT
!
!-----------------------------------------------------------------------
!  Gather weak constraint forcing to storage arrays.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          IF (iic(ng).ne.ntstart(ng)) THEN
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(numthreads)
            DO thread=0,numthreads-1
              subs=NtileX(ng)*NtileE(ng)/numthreads
              DO tile=subs*thread,subs*(thread+1)-1,+1
                CALL frc_ADgather (ng, TILE)
              END DO
            END DO
!$OMP END PARALLEL DO
          END IF
        END DO
# endif
!
!-----------------------------------------------------------------------
!  If appropriate, write out fields into output NetCDF files.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          CALL ad_output (ng)
          IF (exit_flag.ne.NoError) RETURN
        END DO

# ifdef WEAK_CONSTRAINT
!
!-----------------------------------------------------------------------
!  Copy storage arrays index 1 into index 2, and clear index 1.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          IF (MOD(iic(ng)-1,nADJ(ng)).eq.0) THEN
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(numthreads)
            DO thread=0,numthreads-1
              subs=NtileX(ng)*NtileE(ng)/numthreads
              DO tile=subs*thread,subs*(thread+1)-1,+1
                CALL frc_clear (ng, TILE)
              END DO
            END DO
!$OMP END PARALLEL DO
          END IF
        END DO
# endif

# if defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
     defined OPT_OBSERVATIONS || defined SENSITIVITY_4DVAR
!
!-----------------------------------------------------------------------
!  Add appropriate forcing terms to the adjoint model. The form of the
!  forcing depends on the functional whose sensitivity is required.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
#  ifdef SENSITIVITY_4DVAR
          IF (LsenPSAS(ng)) THEN
#  endif
#  if !defined AD_IMPULSE
            IF ((DendS(ng).ge.tdays(ng)).and.                           &
     &          (tdays(ng).ge.DstrS(ng))) THEN
#  endif
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(numthreads)
              DO thread=0,numthreads-1
                subs=NtileX(ng)*NtileE(ng)/numthreads
                DO tile=subs*thread,subs*(thread+1)-1,+1
                  CALL adsen_force (ng, TILE)
                END DO
              END DO
!$OMP END PARALLEL DO
#  if !defined AD_IMPULSE
            END IF
#  endif
#  ifdef SENSITIVITY_4DVAR
          END IF
#  endif
        END DO
# endif
      END DO STEP_LOOP

      RETURN
      END SUBROUTINE ad_main2d
#else
      SUBROUTINE ad_main2d
      RETURN
      END SUBROUTINE ad_main2d
#endif
