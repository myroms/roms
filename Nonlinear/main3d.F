#include "cppdefs.h"
#ifdef SOLVE3D
      SUBROUTINE main3d (ng)
!
!========================================== Alexander F. Shchepetkin ===
!  Copyright (c) ROMS/TOMS Group                                       !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine is the main driver for ROMS/TOMS when configurated     !
!  as a full 3D baroclinic ocean model.  It advances the primitive     !
!  equations for a single time step.                                   !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
      USE mod_stepping
!
# ifdef ANA_VMIX
      USE analytical_mod, ONLY : ana_vmix
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
# if defined SSH_TIDES || defined UV_TIDES
      USE clm_tides_mod, ONLY : clm_tides
# endif
      USE diag_mod, ONLY : diag
# ifdef GLS_MIXING
      USE gls_corstep_mod, ONLY : gls_corstep
      USE gls_prestep_mod, ONLY : gls_prestep
# endif
# ifdef LMD_MIXING
      USE lmd_vmix_mod, ONLY : lmd_vmix
# endif
# ifdef MY25_MIXING
      USE my25_corstep_mod, ONLY : my25_corstep
      USE my25_prestep_mod, ONLY : my25_prestep
# endif
# ifdef ASSIMILATION
      USE oi_update_mod, ONLY : oi_update
# endif
      USE omega_mod, ONLY : omega
# ifndef DIAGNOSTIC
      USE rho_eos_mod, ONLY : rho_eos
# endif
      USE rhs3d_mod, ONLY : rhs3d
# ifdef SEDIMENT
      USE sediment_mod, ONLY : sediment
# endif
      USE set_massflux_mod, ONLY : set_massflux
      USE set_vbc_mod, ONLY : set_vbc
      USE step2d_mod, ONLY : step2d
      USE step3d_t_mod, ONLY : step3d_t
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
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      integer :: my_iif, next_indx1, subs, tile, thread
# ifdef FLOATS
      integer :: Lend, Lstr, chunk_size
# endif
!
!=======================================================================
!  Time-step 3D primitive equations.
!=======================================================================
!
!  Set time indices and time clock.
!
      nstp(ng)=1+MOD(iic(ng)-ntstart,2)
      nnew(ng)=3-nstp(ng)
      nrhs(ng)=nstp(ng)
      time(ng)=time(ng)+dt(ng)
      tdays(ng)=time(ng)*sec2day
!
!-----------------------------------------------------------------------
!  Read in forcing, climatology and assimilation data from input
!  NetCDF files, if any.
!-----------------------------------------------------------------------
!
      CALL get_data (ng)
      if (exit_flag.ne.0) RETURN
!
!-----------------------------------------------------------------------
!  Compute density related fields and horizontal mass fluxes
!  (Hz*u/n and Hz*v/m). Process forcing, climatology and assimilation
!  input data: time interpolate between snapshots.  If applicable,
!  compute and report diagnostics.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*thread,subs*(thread+1)-1,+1
          CALL set_massflux (ng, TILE)
# ifndef DIAGNOSTIC
          CALL rho_eos (ng, TILE)
# endif
          CALL set_data (ng, TILE)
          CALL diag (ng, TILE)
        END DO
      END DO
!$OMP END PARALLEL DO
      if (exit_flag.ne.0) RETURN
!
!-----------------------------------------------------------------------
!  Set fields for vertical boundary conditions. Process tidal forcing,
!  if any.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*thread,subs*(thread+1)-1,+1
# ifdef BULK_FLUXES
          CALL bulk_flux (ng, TILE)
# endif
# ifdef BBL_MODEL
          CALL bblm (ng, TILE)
# endif
          CALL set_vbc (ng, TILE)
# if defined SSH_TIDES || defined UV_TIDES
          CALL clm_tides (ng, TILE)
# endif
        END DO
      END DO
!$OMP END PARALLEL DO
!
!-----------------------------------------------------------------------
!  Compute vertical mixing coefficients for momentum and tracers.
!  Compute S-coordinate vertical velocity, diagnostically from
!  horizontal mass divergence.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*(thread+1)-1,subs*thread,-1
# if defined ANA_VMIX
          CALL ana_vmix (ng, TILE)
# elif defined LMD_MIXING
          CALL lmd_vmix (ng, TILE)
# elif defined BVF_MIXING
          CALL bvf_mix (ng, TILE)
# endif
          CALL omega (ng, TILE)
        END DO
      END DO
!$OMP END PARALLEL DO
!
!-----------------------------------------------------------------------
!  Compute right-hand-side terms for 3D equations.  If applicable,
!  accumulate time-averaged output data.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*thread,subs*(thread+1)-1,+1
          CALL rhs3d (ng, TILE)
# ifdef MY25_MIXING
          CALL my25_prestep (ng, TILE)
# elif defined GLS_MIXING
          CALL gls_prestep (ng, TILE)
# endif
# ifdef AVERAGES
          CALL set_avg (ng, TILE)
# endif
          CALL wvelocity (ng, TILE)
        END DO
      END DO
!$OMP END PARALLEL DO
!
!-----------------------------------------------------------------------
!  If appropriate, write out fields into output NetCDF files.  Notice
!  that IO data is written in delayed and serial mode.  If last time
!  step, stop all clocks.
!-----------------------------------------------------------------------
!
      CALL output (ng)
      if (exit_flag.ne.0) RETURN
      IF (iic(ng).eq.(ntend+1)) THEN
        IF (Master) write(stdout,10) ' Elapsed CPU time (seconds):'
 10     format(/,a,/)
!$OMP PARALLEL DO PRIVATE(thread) SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          CALL wclock_off (ng, 0)
        END DO
!$OMP END PARALLEL DO
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Solve the vertically integrated primitive equations for the
!  free-surface and barotropic momentum components.
!-----------------------------------------------------------------------
!
      DO my_iif=1,nfast(ng)+1
!
!  Set time indices for predictor step. The PREDICTOR_2D_STEP switch
!  it is assumed to be false before the first time-step.
!
        next_indx1=3-indx1(ng)
        IF (.not.PREDICTOR_2D_STEP(ng)) THEN
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
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*(thread+1)-1,subs*thread,-1
            CALL step2d (ng, TILE)
          END DO
        END DO
!$OMP END PARALLEL DO
!
!  Set time indices for corrector step.
!
        IF (PREDICTOR_2D_STEP(ng)) THEN
          PREDICTOR_2D_STEP(ng)=.false.
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
!
!-----------------------------------------------------------------------
!  Time-step 3D equations.
!-----------------------------------------------------------------------
!
!  Time-step 3D momentum equations and couple with vertically
!  integrated equations.
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*(thread+1)-1,subs*thread,-1
          CALL step3d_uv (ng, TILE)
        END DO
      END DO
!$OMP END PARALLEL DO
!
!-----------------------------------------------------------------------
!  Time-step tracer equations.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*thread,subs*(thread+1)-1,+1
          CALL omega (ng, TILE)
# ifdef MY25_MIXING
          CALL my25_corstep (ng, TILE)
# elif defined GLS_MIXING
          CALL gls_corstep (ng, TILE)
# endif
# ifdef SEDIMENT
          CALL sediment (ng, TILE)
# endif
# ifndef DIAGNOSTIC
          CALL step3d_t (ng, TILE)
# endif
# ifdef ASSIMILATION
          CALL oi_update (ng, TILE)
# endif
        END DO
      END DO
!$OMP END PARALLEL DO
# ifdef FLOATS
!
!-----------------------------------------------------------------------
!  Compute Lagrangian drifters trajectories.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,chunk_size,Lstr,Lend)                  &
!$OMP&            SHARED(ng,numthreads,Nfloats)
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
      nf(ng)  =MOD(nf(ng)  +1,NFT+1)
      nfm1(ng)=MOD(nfm1(ng)+1,NFT+1)
      nfm2(ng)=MOD(nfm2(ng)+1,NFT+1)
      nfm3(ng)=MOD(nfm3(ng)+1,NFT+1)
# endif
#else
      SUBROUTINE main3d
#endif
      RETURN
      END SUBROUTINE main3d
