#include "cppdefs.h"
#if defined NONLINEAR && !defined SOLVE3D
      SUBROUTINE main2d (ng)
!
!========================================== Alexander F. Shchepetkin ===
!  Copyright (c) 2005 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This subroutine is the main driver for nonlinear ROMS/TOMS when     !
!  configurated as shallow water (barotropic) ocean model only. It     !
!  advances forward the  vertically integrated primitive equations     !
!  for a single time step.                                             !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
      USE mod_stepping
!
# ifdef AIR_OCEAN
      USE atm_coupler_mod, ONLY : atmos_coupling
# endif
      USE diag_mod, ONLY : diag
# ifdef TLM_CHECK
      USE dotproduct_mod, ONLY : nl_dotproduct
# endif
# ifdef ASSIMILATION
      USE oi_update_mod, ONLY : oi_update
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
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      integer :: next_indx1, subs, tile, thread
# ifdef FLOATS
      integer :: Lend, Lstr, chunk_size
# endif
!
!=======================================================================
!  Time-step vertically integrated equations.
!=======================================================================
!
!  Set time clock.
!
      time(ng)=time(ng)+dt(ng)
      tdays(ng)=time(ng)*sec2day
!
!-----------------------------------------------------------------------
!  Read in required data, if any, from input NetCDF files.
!-----------------------------------------------------------------------
!
      CALL get_data (ng)
      IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Process input data, if any: time interpolate between snapshots.
!  If appropriate, compute and report diagnostics.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile)                             &
!$OMP&            SHARED(ng,Lnew,numthreads)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*thread,subs*(thread+1)-1,+1
          CALL set_data (ng, TILE)
# ifdef AVERAGES
          CALL set_avg (ng, TILE)
# endif
# ifdef DIAGNOSTICS_UV
          CALL set_diags (ng, TILE)
# endif
          CALL diag (ng, TILE)
# ifdef TLM_CHECK
          CALL nl_dotproduct (ng, TILE, Lnew(ng))
# endif
        END DO
      END DO
!$OMP END PARALLEL DO
      IF (exit_flag.ne.NoError) RETURN

# ifdef AIR_OCEAN
!
!-----------------------------------------------------------------------
!  Couple to atmospheric model every "nOCN_steps" timesteps: get
!  air/sea fluxes.
!-----------------------------------------------------------------------
!
      IF (MOD(iic(ng),nOCN_steps).eq.0) THEN
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*(thread+1)-1,subs*thread,-1
            CALL atmos_coupling (ng, TILE)
          END DO
        END DO
!$OMP END PARALLEL DO
      END IF
# endif

!
!-----------------------------------------------------------------------
!  Set vertical boundary conditions. Process tidal forcing.
!-----------------------------------------------------------------------
!
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
!
!-----------------------------------------------------------------------
!  If appropriate, write out fields into output NetCDF files.  Notice
!  that IO data is written in delayed and serial mode.  Exit if last
!  time step.
!-----------------------------------------------------------------------
!
      CALL output (ng)
      IF ((exit_flag.ne.NoError).or.(iic(ng).eq.(ntend+1))) RETURN
!
!-----------------------------------------------------------------------
!  Solve the vertically integrated primitive equations for the
!  free-surface and momentum components.
!-----------------------------------------------------------------------
!
!  Set time indices for predictor step. The PREDICTOR_2D_STEP switch
!  it is assumed to be false before the first time-step.
!
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
!
!  Set time indices for corrector step.
!
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

# ifdef ASSIMILATION
!
!-----------------------------------------------------------------------
!  Assimilation of observations via optimal interpolation.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*thread,subs*(thread+1)-1,+1
          CALL oi_update (ng, TILE)
        END DO
      END DO
!$OMP END PARALLEL DO
# endif

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
      SUBROUTINE main2d
#endif
      RETURN
      END SUBROUTINE main2d
