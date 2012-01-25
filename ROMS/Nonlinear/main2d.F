#include "cppdefs.h"
#if defined NONLINEAR && !defined SOLVE3D
      SUBROUTINE main2d (RunInterval)
!
!svn $Id$
!=======================================================================
!  Copyright (c) 2002-2012 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This subroutine is the main driver for nonlinear ROMS/TOMS when     !
!  configurated as shallow water (barotropic) ocean model only. It     !
!  advances forward  the vertically integrated primitive equations     !
!  for all nested grids,  if any,  by the specified  time interval     !
!  (seconds), RunInterval.                                             !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
# ifdef MODEL_COUPLING
      USE mod_coupler
# endif
      USE mod_iounits
      USE mod_scalars
      USE mod_stepping
!
      USE diag_mod, ONLY : diag
# ifdef TLM_CHECK
      USE dotproduct_mod, ONLY : nl_dotproduct
# endif
# if defined W4DPSAS || defined NLM_OUTER || \
     defined W4DPSAS_SENSITIVITY
      USE forcing_mod, ONLY : forcing
# endif
# ifdef ADJUST_WSTRESS
      USE frc_adjust_mod, ONLY : frc_adjust
# endif
      USE ini_fields_mod, ONLY : ini_fields, ini_zeta
# ifdef NESTING
      USE nesting_mod, ONLY : nesting
# endif
# if defined ADJUST_BOUNDARY
      USE obc_adjust_mod, ONLY : obc_adjust, load_obc
# endif
# ifdef AIR_OCEAN
      USE ocean_coupler_mod, ONLY : ocn2atm_coupling
# endif
# ifdef WAVES_OCEAN
      USE ocean_coupler_mod, ONLY : ocn2wav_coupling
# endif
# ifdef NEARSHORE_MELLOR
      USE radiation_stress_mod, ONLY : radiation_stress
# endif
# ifdef AVERAGES
      USE set_avg_mod, ONLY : set_avg
# endif
# if defined SSH_TIDES || defined UV_TIDES
      USE set_tides_mod, ONLY : set_tides
# endif
      USE set_vbc_mod, ONLY: set_vbc
      USE step2d_mod, ONLY : step2d
# ifdef FLOATS
      USE step_floats_mod, ONLY : step_floats
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
      integer :: next_indx1, ng, subs, tile, thread
# ifdef FLOATS
      integer :: Lend, Lstr, chunk_size
# endif
      real(r8) :: MaxDT, my_StepTime
!
!=======================================================================
!  Time-step vertically integrated equations.
!=======================================================================
!
      my_StepTime=0.0_r8
      MaxDT=MAXVAL(dt)

      STEP_LOOP : DO WHILE (my_StepTime.le.(RunInterval+0.5_r8*MaxDT))

        my_StepTime=my_StepTime+MaxDT
!
!  Set time clock.
!
        DO ng=1,Ngrids
          iic(ng)=iic(ng)+1
          time(ng)=time(ng)+dt(ng)
          tdays(ng)=time(ng)*sec2day
          CALL time_string (time(ng), time_code(ng))
        END DO

# if defined W4DPSAS || defined NLM_OUTER || \
     defined W4DPSAS_SENSITIVITY
!
!-----------------------------------------------------------------------
!  If appropriate, add convolved adjoint solution impulse forcing to
!  the nonlinear model solution. Notice that the forcing is only needed
!  after finishing all the inner loops. The forcing is continuous.
!  That is, it is time interpolated at every time-step from available
!  snapshots (FrequentImpulse=TRUE).
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          IF (FrequentImpulse(ng)) THEN
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(numthreads)
            DO thread=0,numthreads-1
              subs=NtileX(ng)*NtileE(ng)/numthreads
              DO tile=subs*thread,subs*(thread+1)-1,+1
                CALL forcing (ng, TILE, kstp(ng), nstp(ng))
              END DO
            END DO
!$OMP END PARALLEL DO
          END IF
        END DO
# endif
!
!-----------------------------------------------------------------------
!  Read in required data, if any, from input NetCDF files.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          CALL get_data (ng)
          IF (exit_flag.ne.NoError) RETURN
        END DO
!
!-----------------------------------------------------------------------
!  If applicable, process input data: time interpolate between data
!  snapshots.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(numthreads)
          DO thread=0,numthreads-1
            subs=NtileX(ng)*NtileE(ng)/numthreads
            DO tile=subs*thread,subs*(thread+1)-1,+1
              CALL set_data (ng, TILE)
            END DO
          END DO
!$OMP END PARALLEL DO
        END DO
        IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Initialize all time levels and compute other initial fields.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          IF (iic(ng).eq.ntstart(ng)) THEN
!
!  Initialize free-surface.
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
            DO thread=0,numthreads-1
              subs=NtileX(ng)*NtileE(ng)/numthreads
              DO tile=subs*thread,subs*(thread+1)-1,+1
                CALL ini_zeta (ng, TILE, iNLM)
              END DO
            END DO
!$OMP END PARALLEL DO
!
!  Initialize other state variables.
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
            DO thread=0,numthreads-1
              subs=NtileX(ng)*NtileE(ng)/numthreads
              DO tile=subs*(thread+1)-1,subs*thread,-1
                CALL ini_fields (ng, TILE, iNLM)
              END DO
            END DO
!$OMP END PARALLEL DO
          END IF
        END DO
!
!-----------------------------------------------------------------------
!  Compute and report diagnostics. If appropriate, accumulate time-
!  averaged output data which needs a irreversible loop in shared-memory
!  jobs.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(Lnew,numthreads)
          DO thread=0,numthreads-1
            subs=NtileX(ng)*NtileE(ng)/numthreads
            DO tile=subs*thread,subs*(thread+1)-1,+1      ! irreversible
# ifdef AVERAGES
              CALL set_avg (ng, TILE)
# endif
# ifdef DIAGNOSTICS
              CALL set_diags (ng, TILE)
# endif
              CALL diag (ng, TILE)
# ifdef TLM_CHECK
              CALL nl_dotproduct (ng, TILE, Lnew(ng))
# endif
            END DO
          END DO
!$OMP END PARALLEL DO
        END DO
        IF (exit_flag.ne.NoError) RETURN

# ifdef AIR_OCEAN
!
!-----------------------------------------------------------------------
!  Couple ocean to atmosphere model every "CoupleSteps(Iatmos)"
!  timesteps: get air/sea fluxes.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          IF ((iic(ng).ne.ntstart(ng)).and.                             &
     &        MOD(iic(ng),CoupleSteps(Iatmos,ng)).eq.0) THEN
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
            DO thread=0,numthreads-1
              subs=NtileX(ng)*NtileE(ng)/numthreads
              DO tile=subs*(thread+1)-1,subs*thread,-1
                CALL ocn2atm_coupling (ng, TILE)
              END DO
            END DO
!$OMP END PARALLEL DO
          END IF
        END DO
# endif

# ifdef ADJUST_BOUNDARY
!
!-----------------------------------------------------------------------
!  Interpolate open boundary increments and adjust open boundary.
!  Load open boundary into storage arrays. Skip the last output
!  timestep.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          IF (iic(ng).lt.(ntend(ng)+1)) THEN
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
            DO thread=0,numthreads-1
              subs=NtileX(ng)*NtileE(ng)/numthreads
              DO tile=subs*thread,subs*(thread+1)-1,+1
                CALL obc_adjust (ng, TILE, Lbinp(ng))
                CALL load_obc (ng, TILE, Lbout(ng))
              END DO
            END DO
!$OMP END PARALLEL DO
          END IF
        END DO
# endif

# ifdef ADJUST_WSTRESS
!
!-----------------------------------------------------------------------
!  Interpolate surface forcing increments and adjust surface forcing.
!  Load surface forcing to storage arrays.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          IF (iic(ng).lt.(ntend(ng)+1)) THEN
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
            DO thread=0,numthreads-1
              subs=NtileX(ng)*NtileE(ng)/numthreads
              DO tile=subs*thread,subs*(thread+1)-1,+1
                CALL frc_adjust (ng, TILE, Lfinp(ng))
              END DO
            END DO
!$OMP END PARALLEL DO
          END IF
        END DO
# endif

# ifdef WAVES_OCEAN
!
!-----------------------------------------------------------------------
!  Couple ocean to waves model every "CoupleSteps(Iwaves)"
!  timesteps: get waves/sea fluxes.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          IF ((iic(ng).ne.ntstart(ng)).and.                             &
     &        MOD(iic(ng)-1,CoupleSteps(Iwaves,ng)).eq.0) THEN
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
            DO thread=0,numthreads-1
              subs=NtileX(ng)*NtileE(ng)/numthreads
              DO tile=subs*thread,subs*(thread+1)-1,+1
                CALL ocn2wav_coupling (ng, TILE)
              END DO
            END DO
!$OMP END PARALLEL DO
          END IF
        END DO
# endif

# ifdef NEARSHORE_MELLOR
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
              CALL radiation_stress (ng, TILE)
            END DO
          END DO
!$OMP END PARALLEL DO
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
              CALL set_vbc (ng, TILE)
# if defined SSH_TIDES || defined UV_TIDES
              CALL set_tides (ng, TILE)
# endif
            END DO
          END DO
!$OMP END PARALLEL DO
        END DO
# ifdef NESTING
        CALL nesting (100)
# endif
!
!-----------------------------------------------------------------------
!  If appropriate, write out fields into output NetCDF files.  Notice
!  that IO data is written in delayed and serial mode.  Exit if last
!  time step.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          CALL output (ng)
          IF ((exit_flag.ne.NoError).or.                                &
     &        ((iic(ng).eq.(ntend(ng)+1)).and.(ng.eq.Ngrids))) RETURN
        END DO
!
!-----------------------------------------------------------------------
!  Solve the vertically integrated primitive equations for the
!  free-surface and momentum components.
!-----------------------------------------------------------------------
!
!  Set time indices for predictor step. The PREDICTOR_2D_STEP switch
!  it is assumed to be false before the first time-step.
!
        DO ng=1,Ngrids
          iif(ng)=1
          nfast(ng)=1
          next_indx1=3-indx1(ng)
          IF (.not.PREDICTOR_2D_STEP(ng)) THEN
            PREDICTOR_2D_STEP(ng)=.TRUE.
            IF (FIRST_2D_STEP) THEN
              kstp(ng)=indx1(ng)
            ELSE
              kstp(ng)=3-indx1(ng)
            END IF
            knew(ng)=3
            krhs(ng)=indx1(ng)
          END IF
!
!  Predictor step - Advance barotropic equations using 2D time-step
!  ==============   predictor scheme.
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
          DO thread=0,numthreads-1
            subs=NtileX(ng)*NtileE(ng)/numthreads
            DO tile=subs*(thread+1)-1,subs*thread,-1
              CALL step2d (ng, TILE)
            END DO
          END DO
!$OMP END PARALLEL DO
        END DO
# ifdef NESTING
        CALL nesting (100)
# endif
!
!  Set time indices for corrector step.
!
        DO ng=1,Ngrids
          IF (PREDICTOR_2D_STEP(ng)) THEN
            PREDICTOR_2D_STEP(ng)=.FALSE.
            knew(ng)=next_indx1
            kstp(ng)=3-knew(ng)
            krhs(ng)=3
            IF (iif(ng).lt.(nfast(ng)+1)) indx1(ng)=next_indx1
          END IF
!
!  Corrector step - Apply 2D time-step corrector scheme.  Notice that
!  ==============   there is not need for a corrector step during the
!  auxiliary (nfast+1) time-step.
!
          IF (iif(ng).lt.(nfast(ng)+1)) THEN
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
            DO thread=0,numthreads-1
              subs=NtileX(ng)*NtileE(ng)/numthreads
              DO tile=subs*thread,subs*(thread+1)-1,+1
                CALL step2d (ng, TILE)
              END DO
            END DO
!$OMP END PARALLEL DO
          END IF
        END DO
# ifdef NESTING
        CALL nesting (100)
# endif

# ifdef FLOATS
!
!-----------------------------------------------------------------------
!  Compute Lagrangian drifters trajectories.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          IF (Lfloats(Ng)) THEN
!$OMP PARALLEL DO PRIVATE(thread,chunk_size,Lstr,Lend)                  &
!$OMP&            SHARED(numthreads,Nfloats)
            DO thread=0,numthreads-1
              chunk_size=(Nfloats(ng)+numthreads-1)/numthreads
              Lstr=1+thread*chunk_size
              Lend=MIN(Nfloats(ng),Lstr+chunk_size-1)
              CALL step_floats (ng, Lstr, Lend)
            END DO
!$OMP END PARALLEL DO
!
!  Shift floats time indices.
!
            nfp1(ng)=MOD(nfp1(ng)+1,NFT+1)
            nf  (ng)=MOD(nf  (ng)+1,NFT+1)
            nfm1(ng)=MOD(nfm1(ng)+1,NFT+1)
            nfm2(ng)=MOD(nfm2(ng)+1,NFT+1)
            nfm3(ng)=MOD(nfm3(ng)+1,NFT+1)
          END IF
        END DO
# endif
      END DO STEP_LOOP

#else
      SUBROUTINE main2d
#endif
      RETURN
      END SUBROUTINE main2d
