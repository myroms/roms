#include "cppdefs.h"
      SUBROUTINE initial (ng)
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine initializes all model variables.                       !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
#ifdef BBL_MODEL
      USE mod_bbl
#endif
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
      USE mod_ocean
      USE mod_scalars
      USE mod_stepping
!
      USE analytical_mod
      USE horz_mix_mod, ONLY : horz_mix
      USE metrics_mod, ONLY : metrics
#ifdef SOLVE3D
      USE set_depth_mod, ONLY : set_depth
      USE omega_mod, ONLY : omega
      USE rho_eos_mod, ONLY : rho_eos
      USE set_massflux_mod, ONLY : set_massflux
      USE stiffness_mod, ONLY : stiffness
#endif
#if defined MASKING && (defined READ_WATER || defined WRITE_WATER)
      USE wpoints_mod, ONLY : wpoints
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
      logical :: update = .false.

      integer :: LBi, UBi, LBj, UBj
      integer :: subs, tile, thread

      integer :: my_numthreads
!
!=======================================================================
!   Initialize model variables.
!=======================================================================
!
      IF (Master) THEN
#ifdef PERTURBATION
        WRITE (stdout,10) Nrun
 10     FORMAT (/,' <<<< Ensemble/Perturbation Run: ',i5.5,' >>>>',/)
#endif
        WRITE (stdout,20) ' INITIAL: Configurating and initializing ...'
 20     FORMAT (/,a,/)
      END IF
!
!-----------------------------------------------------------------------
!  Initialize time stepping indices and counters.
!-----------------------------------------------------------------------
!
      iif(ng)=1
      indx1(ng)=1
      kstp(ng)=1
      krhs(ng)=1
      knew(ng)=1
      PREDICTOR_2D_STEP=.false.
!
      iic(ng)=0
      ntend=ntimes
      nstp(ng)=1
      nrhs(ng)=1
      nnew(ng)=1
#ifdef FLOATS
      nf(ng)=0
      nfp1(ng)=1
      nfm1(ng)=4
      nfm2(ng)=3
      nfm3(ng)=2
#endif
      tdays(ng)=dstart
      time(ng)=tdays(ng)*day2sec
      LBi=lbound(GRID(ng)%h,DIM=1)
      UBi=ubound(GRID(ng)%h,DIM=1)
      LBj=lbound(GRID(ng)%h,DIM=2)
      UBj=ubound(GRID(ng)%h,DIM=2)
!
!-----------------------------------------------------------------------
!  Start time wall clocks.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread) SHARED(ng,numthreads)
      DO thread=0,numthreads-1
        CALL wclock_on (ng, 0)
#ifdef PROFILE
        CALL wclock_on (ng, 1)
#endif
      END DO
!$OMP END PARALLEL DO
!
!=======================================================================
!  On first pass of ensemble run loop, initialize model configuration.
!=======================================================================
!
      IF (Nrun.eq.1) THEN
!
!-----------------------------------------------------------------------
!  Set horizontal grid, bathymetry, and Land/Sea masking (if any).
!  Use analytical functions or read in from a grid NetCDF.
!-----------------------------------------------------------------------
!
#ifdef ANA_GRID
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1
            CALL ana_grid (ng, TILE)
# ifdef MASKING
            CALL ana_mask (ng, TILE)
# endif
          END DO
        END DO
!$OMP END PARALLEL DO
#else
        CALL get_grid (ng)
        if (exit_flag.ne.0) RETURN
#endif

#ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Set vertical S-coordinate transformation function.
!-----------------------------------------------------------------------
!
        CALL set_scoord (ng)
#endif

#ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Set barotropic time-steps average weighting function.
!-----------------------------------------------------------------------
!
        CALL set_weights (ng)
#endif
!
!-----------------------------------------------------------------------
!  Compute various metric term combinations.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1
            CALL metrics (ng, TILE)
#if defined MASKING && (defined READ_WATER || defined WRITE_WATER)
            CALL wpoints (ng, TILE)
#endif
          END DO
        END DO
!$OMP END PARALLEL DO

#if defined VISC_GRID || defined DIFF_GRID || defined SPONGE
!
!-----------------------------------------------------------------------
!  Set horizontal mixing coefficients. Rescale according to the local
!  grid size. If applicable, increases horizontal mixing in sponge
!  areas.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1
            CALL horz_mix (ng, TILE)
          END DO
        END DO
!$OMP END PARALLEL DO
#endif

#ifdef NUDGING_COFF
!
!-----------------------------------------------------------------------
!  If appropriate, set nudging coefficiests time scales.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1
            CALL set_nudgcof (ng, TILE)
          END DO
        END DO
!$OMP END PARALLEL DO
#endif
      END IF
!
!=======================================================================
!  Initialize model state variables and forcing.  This part is
!  executed for each ensemble run.
!=======================================================================

#if defined SOLVE3D && defined ANA_INITIAL
!
!-----------------------------------------------------------------------
!  Compute initial time-evolving depths with zero free-surface.  The
!  grid depths are sometimes needed during analytical initialization.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*thread,subs*(thread+1)-1
          CALL set_depth (ng, TILE)
        END DO
      END DO
!$OMP END PARALLEL DO
#endif
!
!-----------------------------------------------------------------------
!  Set primitive variables initial conditions.  Use analytical
!  functions or read from an initial or restart NetCDF file.
!-----------------------------------------------------------------------
!
#ifdef ANA_INITIAL
      IF (nrrec.eq.0) THEN
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1
            CALL ana_initial (ng, TILE)
          END DO
        END DO
!$OMP END PARALLEL DO
      ELSE
        CALL get_initial (ng)
      END IF
#else
      CALL get_initial (ng)
#endif

#if defined BIOLOGY && defined ANA_BIOLOGY && defined SOLVE3D
!
!  Analytical initial conditions for biology.
!
      IF (nrrec.eq.0) THEN
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1
            CALL ana_biology (ng, TILE)
          END DO
        END DO
!$OMP END PARALLEL DO
      END IF
#endif

#if (defined SEDIMENT || defined BBL_MODEL) && defined ANA_SEDIMENT && \
     defined SOLVE3D
!
!  Analytical initial conditions for sediment.
!
      IF (nrrec.eq.0) THEN
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1
            CALL ana_sediment (ng, TILE)
          END DO
        END DO
!$OMP END PARALLEL DO
      END IF
#endif
      IF (exit_flag.ne.0) RETURN

#ifdef PERTURBATION
!
!-----------------------------------------------------------------------
!  Perturb initial conditions with a ramdon linear combination of
!  error subspaces.
!-----------------------------------------------------------------------
!
# ifdef RANDOM_ESPERT
      CALL res_pert (ng)
# endif
#endif

#ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Compute initial time-evolving depths.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*thread,subs*(thread+1)-1
          CALL set_depth (ng, TILE)
        END DO
      END DO
!$OMP END PARALLEL DO
!
!-----------------------------------------------------------------------
!  Compute initial horizontal mass fluxes, Hz*u/n and Hz*v/m.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*thread,subs*(thread+1)-1
          CALL set_massflux (ng, TILE)
        END DO
      END DO
!$OMP END PARALLEL DO
!
!-----------------------------------------------------------------------
!  Compute initial S-coordinates vertical velocity. Compute initial
!  density anomaly from potential temperature and salinity via equation
!  of state for seawater.  Also compute other equation of state related
!  quatities.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*thread,subs*(thread+1)-1
          CALL omega (ng, TILE)
          CALL rho_eos (ng, TILE)
        END DO
      END DO
!$OMP END PARALLEL DO
#endif
!
!-----------------------------------------------------------------------
!  Read in initial forcing, climatology and assimilation data from
!  input NetCDF files.  It loads the first relevant data record for
!  the time-interpolation between snapshots.
!-----------------------------------------------------------------------
!
      CALL get_data (ng)
      IF (exit_flag.ne.0) RETURN

#ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Compute grid stiffness.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*thread,subs*(thread+1)-1
          CALL stiffness (ng, TILE)
        END DO
      END DO
!$OMP END PARALLEL DO
#endif

#if defined FLOATS || defined STATIONS
!
!-----------------------------------------------------------------------
!  If applicable, convert initial locations to fractional grid
!  coordinates.
!-----------------------------------------------------------------------
!
      CALL grid_coords (ng)
#endif

#ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off initiialization time wall clock.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread) SHARED(ng,numthreads)
      DO thread=0,numthreads-1
        CALL wclock_off (ng, 1)
      END DO
!$OMP END PARALLEL DO
#endif
      RETURN
      END SUBROUTINE initial
