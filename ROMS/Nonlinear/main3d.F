#include "cppdefs.h"
#if defined NONLINEAR && defined SOLVE3D
      SUBROUTINE main3d (RunInterval)
!
!svn $Id$
!=======================================================================
!  Copyright (c) 2002-2012 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This subroutine is the main driver for nonlinear ROMS/TOMS when     !
!  configurated as a full 3D baroclinic ocean model.  It  advances     !
!  forward the primitive equations for all  nested  grids, if any,     !
!  for the specified time interval (seconds), RunInterval.             !
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
# ifdef ANA_VMIX
      USE analytical_mod, ONLY : ana_vmix
# endif
# ifdef BIOLOGY
      USE biology_mod, ONLY : biology
# endif
# ifdef BBL_MODEL
      USE bbl_mod, ONLY : bblm
# endif
# ifdef BULK_FLUXES
      USE bulk_flux_mod, ONLY : bulk_flux
# endif
# ifdef BVF_MIXING
      USE bvf_mix_mod, ONLY : bvf_mix
# endif
      USE diag_mod, ONLY : diag
# ifdef TLM_CHECK
      USE dotproduct_mod, ONLY : nl_dotproduct
# endif
# if defined W4DPSAS || defined NLM_OUTER || \
     defined W4DPSAS_SENSITIVITY
      USE forcing_mod, ONLY : forcing
# endif
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      USE frc_adjust_mod, ONLY : frc_adjust
# endif
# ifdef GLS_MIXING
      USE gls_corstep_mod, ONLY : gls_corstep
      USE gls_prestep_mod, ONLY : gls_prestep
# endif
# if defined DIFF_3DCOEF || defined VISC_3DCOEF
      USE hmixing_mod, ONLY : hmixing
# endif
      USE ini_fields_mod, ONLY : ini_fields, ini_zeta
# ifdef LMD_MIXING
      USE lmd_vmix_mod, ONLY : lmd_vmix
# endif
# ifdef MY25_MIXING
      USE my25_corstep_mod, ONLY : my25_corstep
      USE my25_prestep_mod, ONLY : my25_prestep
# endif
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
      USE omega_mod, ONLY : omega
# ifdef NEARSHORE_MELLOR
      USE radiation_stress_mod, ONLY : radiation_stress
# endif
# ifndef TS_FIXED
      USE rho_eos_mod, ONLY : rho_eos
# endif
      USE rhs3d_mod, ONLY : rhs3d
# ifdef SEDIMENT
      USE sediment_mod, ONLY : sediment
# endif
# ifdef AVERAGES
      USE set_avg_mod, ONLY : set_avg
# endif
      USE set_depth_mod, ONLY : set_depth
      USE set_massflux_mod, ONLY : set_massflux
# if defined SSH_TIDES || defined UV_TIDES
      USE set_tides_mod, ONLY : set_tides
# endif
      USE set_vbc_mod, ONLY : set_vbc
      USE set_zeta_mod, ONLY : set_zeta
      USE step2d_mod, ONLY : step2d
# ifndef TS_FIXED
      USE step3d_t_mod, ONLY : step3d_t
# endif
      USE step3d_uv_mod, ONLY : step3d_uv
# ifdef FLOATS
      USE step_floats_mod, ONLY : step_floats
# endif
      USE wvelocity_mod, ONLY : wvelocity
!
      implicit none
!
!  Imported variable declarations.
!
      real(r8), intent(in) :: RunInterval
!
!  Local variable declarations.
!
      integer :: ng, tile
      integer :: my_iif, next_indx1
# ifdef FLOATS
      integer :: Lend, Lstr, chunk_size
# endif
      real(r8) :: MaxDT, my_StepTime
!
!=======================================================================
!  Time-step nonlinear 3D primitive equations by the specified time.
!=======================================================================
!
      my_StepTime=0.0_r8
      MaxDT=MAXVAL(dt)

      STEP_LOOP : DO WHILE (my_StepTime.le.(RunInterval+0.5_r8*MaxDT))

        my_StepTime=my_StepTime+MaxDT
!
!  Set time indices and time clock.
!
        DO ng=1,Ngrids
          iic(ng)=iic(ng)+1
          nstp(ng)=1+MOD(iic(ng)-ntstart(ng),2)
          nnew(ng)=3-nstp(ng)
          nrhs(ng)=nstp(ng)
!$OMP MASTER
          time(ng)=time(ng)+dt(ng)
          tdays(ng)=time(ng)*sec2day
          CALL time_string (time(ng), time_code(ng))
!$OMP END MASTER
        END DO
!$OMP BARRIER
!
!-----------------------------------------------------------------------
!  Read in required data, if any, from input NetCDF files.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
!$OMP MASTER
          CALL get_data (ng)
!$OMP END MASTER
!$OMP BARRIER
          IF (exit_flag.ne.NoError) RETURN
        END DO
!
!-----------------------------------------------------------------------
!  If applicable, process input data: time interpolate between data
!  snapshots.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL set_data (ng, tile)
          END DO
!$OMP BARRIER
        END DO
        IF (exit_flag.ne.NoError) RETURN

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
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL forcing (ng, tile, kstp(ng), nstp(ng))
            END DO
!$OMP BARRIER
          END IF
        END DO
# endif
!
!-----------------------------------------------------------------------
!  Initialize all time levels and compute other initial fields.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          IF (iic(ng).eq.ntstart(ng)) THEN
!
!  Initialize free-surface and compute initial level thicknesses and
!  depths.
!
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL ini_zeta (ng, tile, iNLM)
              CALL set_depth (ng, tile)
            END DO
!$OMP BARRIER
!
!  Initialize other state variables.
!
            DO tile=last_tile(ng),first_tile(ng),-1
              CALL ini_fields (ng, tile, iNLM)
            END DO
!$OMP BARRIER
          END IF
        END DO
!
!-----------------------------------------------------------------------
!  Compute horizontal mass fluxes (Hz*u/n and Hz*v/m), density related
!  quatities and report global diagnostics.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL set_massflux (ng, tile)
# ifndef TS_FIXED
            CALL rho_eos (ng, tile)
# endif
            CALL diag (ng, tile)
# ifdef TLM_CHECK
            CALL nl_dotproduct (ng, tile, Lnew(ng))
# endif
          END DO
!$OMP BARRIER
        END DO
        IF (exit_flag.ne.NoError) RETURN
# ifdef NESTING
        CALL nesting (5)
# endif

# ifdef AIR_OCEAN
!
!-----------------------------------------------------------------------
!  Couple ocean to atmosphere model every "CoupleSteps(Iatmos)"
!  timesteps: get air/sea fluxes.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          IF ((iic(ng).ne.ntstart(ng)).and.                             &
     &        MOD(iic(ng)-1,CoupleSteps(Iatmos,ng)).eq.0) THEN
            DO tile=last_tile(ng),first_tile(ng),-1
              CALL ocn2atm_coupling (ng, tile)
            END DO
!$OMP BARRIER
          END IF
        END DO
# endif

# ifdef WAVES_OCEAN
!
!-----------------------------------------------------------------------
!  Couple to ocean to waves model every "CoupleSteps(Iwaves)"
!  timesteps: get waves/ocean fluxes.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          IF ((iic(ng).ne.ntstart(ng)).and.                             &
     &        MOD(iic(ng)-1,CoupleSteps(Iwaves,ng)).eq.0) THEN
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL ocn2wav_coupling (ng, tile)
            END DO
!$OMP BARRIER
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
          DO tile=last_tile(ng),first_tile(ng),-1
            CALL radiation_stress (ng, tile)
          END DO
!$OMP BARRIER
        END DO
#  ifdef NESTING
        CALL nesting (8)
#  endif
# endif
!
!-----------------------------------------------------------------------
!  Set fields for vertical boundary conditions. Process tidal forcing,
!  if any.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1
# ifdef BULK_FLUXES
#  if defined FOUR_DVAR && defined NL_BULK_FLUXES
            IF (Nrun.eq.1) CALL bulk_flux (ng, tile)
#  else
            CALL bulk_flux (ng, tile)
#  endif
# endif
# ifdef BBL_MODEL
            CALL bblm (ng, tile)
# endif
            CALL set_vbc (ng, tile)
# if defined SSH_TIDES || defined UV_TIDES
            CALL set_tides (ng, tile)
# endif
          END DO
!$OMP BARRIER
        END DO
# ifdef NESTING
        CALL nesting (9)
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
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL obc_adjust (ng, tile, Lbinp(ng))
              CALL load_obc (ng, tile, Lbout(ng))
            END DO
!$OMP BARRIER
          END IF
        END DO
# endif

# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
!
!-----------------------------------------------------------------------
!  Interpolate surface forcing increments and adjust surface forcing.
!  Load surface forcing into storage arrays. Skip the last output
!  timestep.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          IF (iic(ng).lt.(ntend(ng)+1)) THEN
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL frc_adjust (ng, tile, Lfinp(ng))
            END DO
!$OMP BARRIER
          END IF
        END DO
# endif
!
!-----------------------------------------------------------------------
!  Compute time-dependent vertical/horizontal mixing coefficients for
!  momentum and tracers. Compute S-coordinate vertical velocity,
!  diagnostically from horizontal mass divergence.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          DO tile=last_tile(ng),first_tile(ng),-1
# if defined ANA_VMIX
            CALL ana_vmix (ng, tile, iNLM)
# elif defined LMD_MIXING
            CALL lmd_vmix (ng, tile)
# elif defined BVF_MIXING
            CALL bvf_mix (ng, tile)
# endif
# if defined DIFF_3DCOEF || defined VISC_3DCOEF
            CALL hmixing (ng, tile)
# endif
            CALL omega (ng, tile)
            CALL wvelocity (ng, tile, nstp(ng))
          END DO
!$OMP BARRIER
        END DO
# ifdef NESTING
        CALL nesting (10)
# endif
!
!-----------------------------------------------------------------------
!  Set free-surface to it time-averaged value.  If applicable,
!  accumulate time-averaged output data which needs a irreversible
!  loop in shared-memory jobs.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1         ! irreversible
            CALL set_zeta (ng, tile)
#  ifdef DIAGNOSTICS
            CALL set_diags (ng, tile)
#  endif
#  ifdef AVERAGES
            CALL set_avg (ng, tile)
#  endif
          END DO
!$OMP BARRIER
        END DO
#  ifdef NESTING
        CALL nesting (11)
#  endif
!
!-----------------------------------------------------------------------
!  If appropriate, write out fields into output NetCDF files.  Notice
!  that IO data is written in delayed and serial mode.  Exit if last
!  time step.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
!$OMP MASTER
          CALL output (ng)
!$OMP END MASTER
!$OMP BARRIER
          IF ((exit_flag.ne.NoError).or.                                &
     &        ((iic(ng).eq.(ntend(ng)+1)).and.(ng.eq.Ngrids))) RETURN
        END DO
!
!-----------------------------------------------------------------------
!  Compute right-hand-side terms for 3D equations.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          DO tile=last_tile(ng),first_tile(ng),-1
            CALL rhs3d (ng, tile)
#  ifdef MY25_MIXING
            CALL my25_prestep (ng, tile)
#  elif defined GLS_MIXING
            CALL gls_prestep (ng, tile)
#  endif
          END DO
!$OMP BARRIER
        END DO
#  ifdef NESTING
        CALL nesting (12)
#  endif
!
!-----------------------------------------------------------------------
!  Solve the vertically integrated primitive equations for the
!  free-surface and barotropic momentum components.
!-----------------------------------------------------------------------
!
        LOOP_2D : DO my_iif=1,MAXVAL(nfast)+1
!
!  Set time indices for predictor step. The PREDICTOR_2D_STEP switch
!  it is assumed to be false before the first time-step.
!
          DO ng=1,Ngrids
            next_indx1=3-indx1(ng)
            IF (.not.PREDICTOR_2D_STEP(ng).and.                         &
     &          my_iif.le.(nfast(ng)+1)) THEN
              PREDICTOR_2D_STEP(ng)=.TRUE.
              iif(ng)=my_iif
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
!  ==============   predictor scheme.  No actual time-stepping is
!  performed during the auxiliary (nfast+1) time-step. It is needed
!  to finalize the fast-time averaging of 2D fields, if any, and
!  compute the new time-evolving depths.
!
            IF (my_iif.le.(nfast(ng)+1)) THEN
              DO tile=last_tile(ng),first_tile(ng),-1
                CALL step2d (ng, tile)
              END DO
!$OMP BARRIER
            END IF
          END DO
# ifdef NESTING
          CALL nesting (14)
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
              DO tile=first_tile(ng),last_tile(ng),+1
                CALL step2d (ng, tile)
              END DO
!$OMP BARRIER
            END IF
          END DO
# ifdef NESTING
          CALL nesting (15)
# endif
        END DO LOOP_2D
# ifdef NESTING
        CALL nesting (21)
# endif
!
!-----------------------------------------------------------------------
!  Recompute depths and thicknesses using the new time filtered
!  free-surface.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          DO tile=last_tile(ng),first_tile(ng),-1
            CALL set_depth (ng, tile)
          END DO
!$OMP BARRIER
        END DO
!
!-----------------------------------------------------------------------
!  Time-step 3D momentum equations.
!-----------------------------------------------------------------------
!
!  Time-step 3D momentum equations and couple with vertically
!  integrated equations.
!
        DO ng=1,Ngrids
          DO tile=last_tile(ng),first_tile(ng),-1
            CALL step3d_uv (ng, tile)
          END DO
!$OMP BARRIER
        END DO
# ifdef NESTING
        CALL nesting (16)
# endif
!
!-----------------------------------------------------------------------
!  Time-step vertical mixing turbulent equations and passive tracer
!  source and sink terms, if applicable.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL omega (ng, tile)
# ifdef MY25_MIXING
            CALL my25_corstep (ng, tile)
# elif defined GLS_MIXING
            CALL gls_corstep (ng, tile)
# endif
# ifdef BIOLOGY
            CALL biology (ng, tile)
# endif
# ifdef SEDIMENT
            CALL sediment (ng, tile)
# endif
          END DO
!$OMP BARRIER
        END DO
# ifdef NESTING
        CALL nesting (17)
# endif

# ifndef TS_FIXED
!
!-----------------------------------------------------------------------
!  Time-step tracer equations.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          DO tile=last_tile(ng),first_tile(ng),-1
            CALL step3d_t (ng, tile)
          END DO
!$OMP BARRIER
        END DO
#  ifdef NESTING
        CALL nesting (18)
#  endif
# endif

# ifdef FLOATS
!
!-----------------------------------------------------------------------
!  Compute Lagrangian drifters trajectories: Split all the drifters
!  between all the computational threads, except in distributed-memory
!  and serial configurations. In distributed-memory, the parallel node
!  containing the drifter is selected internally since the state
!  variables do not have a global scope.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          IF (Lfloats(Ng)) THEN
#  ifdef _OPENMP
            chunk_size=(Nfloats(ng)+numthreads-1)/numthreads
            Lstr=1+MyThread*chunk_size
            Lend=MIN(Nfloats(ng),Lstr+chunk_size-1)
#  else
            Lstr=1
            Lend=Nfloats(ng)
#  endif
            CALL step_floats (ng, Lstr, Lend)
!$OMP BARRIER
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

      RETURN
      END SUBROUTINE main3d
#else
      SUBROUTINE main3d
      RETURN
      END SUBROUTINE main3d
#endif
