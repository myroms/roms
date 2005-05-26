#include "cppdefs.h"
      SUBROUTINE initial (ng)
!
!=======================================================================
!  Copyright (c) 2005 ROMS/TOMS Group                                  !
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
#ifdef FOUR_DVAR
      USE mod_fourdvar
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
#ifdef IS4DVAR
      USE ini_adjust_mod, ONLY : ini_adjust
#endif
#ifdef TLM_CHECK
      USE ini_adjust_mod, ONLY : ini_perturb
#endif
      USE metrics_mod, ONLY : metrics
#ifdef SOLVE3D
      USE set_depth_mod, ONLY : set_depth
      USE omega_mod, ONLY : omega
      USE rho_eos_mod, ONLY : rho_eos
      USE set_massflux_mod, ONLY : set_massflux
#endif
      USE stiffness_mod, ONLY : stiffness
#if defined MASKING && \
   (defined READ_WATER || defined WRITE_WATER || defined PROPAGATOR)
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
      logical :: update = .FALSE.
      logical :: FirstPass = .TRUE.

      integer :: LBi, UBi, LBj, UBj
      integer :: IniRec, Tindex, subs, tile, thread

      integer :: my_numthreads
!
!=======================================================================
!   Initialize model variables.
!=======================================================================
!
      IF (Master) THEN
#if defined PERTURBATION
        WRITE (stdout,10) Nrun
 10     FORMAT (/,' <<<< Ensemble/Perturbation Run: ',i5.5,' >>>>',/)
#elif defined FOUR_DVAR
        WRITE (stdout,10) Nrun, Ipass
 10     FORMAT (/,' <<<< 4D Variational Data Assimilation, ',           &
     &          'Iteration = ',i5.5, ', Ipass = ',i1,' >>>>',/)
#endif
        WRITE (stdout,20) 'INITIAL: Configurating and initializing ',   &
     &                    'forward nonlinear model ...'
  20    FORMAT (/,1x,a,a,/)
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
      PREDICTOR_2D_STEP(ng)=.FALSE.
      synchro_flag(ng)=.TRUE.
!
      iic(ng)=0
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
      ntstart=INT((time(ng)-dstart*day2sec)/dt(ng))+1
      ntend=ntimes
      ntfirst=ntstart

      IniRec=nrrec
      Tindex=1

      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)
!
!-----------------------------------------------------------------------
!  Start time wall clocks.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread) SHARED(ng,numthreads)
      DO thread=0,numthreads-1
#ifdef PROFILE
        CALL wclock_on (ng, iNLM, 1)
#endif
      END DO
!$OMP END PARALLEL DO

#ifdef FOUR_DVAR
!
!-----------------------------------------------------------------------
!  If variational data assimilation, reset several IO switches and
!  variables.
!-----------------------------------------------------------------------
!
!  Set switch to create (true) nonlinear model initial conditions and
!  histroy NetCDF files or append (false) to an existing file.  Set
!  record to read from initial NetCDF file.
!
      IF ((Nrun.eq.ERstr).and.(Ipass.le.1)) THEN
# ifdef ANA_INITIAL
        LdefINI(ng)=.TRUE.
# endif
        LdefHIS(ng)=.TRUE.    
        CALL def_ini (ng)
        IF (exit_flag.ne.NoError) RETURN
        IniRec=nrrec
        tINIindx(ng)=IniRec
      ELSE
        IniRec=tINIindx(ng)
      END IF
!
!  Reset nonlinear history time record counters. These counters are
!  reset on every iteration pass. This file is created on the first
!  iteration pass.
!
      tHISindx(ng)=0
      NrecHIS(ng)=0

# if defined S4DVAR
!  
!  Avoid writing data into average, history, and restart files in
!  the second iteration pass (Ipass=2).
!
      IF (Ipass.eq.1) THEN
        LwrtAVG(ng)=.TRUE.
        LwrtHIS(ng)=.TRUE.
        LwrtRST(ng)=.TRUE.
      ELSEif (Ipass.eq.2) THEN
        LwrtAVG(ng)=.FALSE.
        LwrtHIS(ng)=.FALSE.
        LwrtRST(ng)=.FALSE.
      END IF

# elif defined IS4DVAR
!
!  Activate switches to writting data into average, history and
!  restart files.
!
      LwrtAVG(ng)=.TRUE.
      LwrtHIS(ng)=.TRUE.
      LwrtRST(ng)=.TRUE.
# endif
!
!  Open observations NetCDF file and initialize various variables
!  needed for processing the nonlinear state solution at observation
!  locations.
!
      CALL obs_initial (ng, iNLM, .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
#endif
!
!=======================================================================
!  On first pass of ensemble/perturbation/iteration loop, initialize
!  model configuration.
!=======================================================================
!
      IF (Nrun.eq.ERstr) THEN
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
            CALL ana_grid (ng, TILE, iNLM)
# ifdef MASKING
            CALL ana_mask (ng, TILE, iNLM)
# endif
          END DO
        END DO
!$OMP END PARALLEL DO
#else
        CALL get_grid (ng, iNLM)
        if (exit_flag.ne.NoError) RETURN
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
            CALL metrics (ng, TILE, iNLM)
#if defined MASKING && \
   (defined READ_WATER || defined WRITE_WATER || defined PROPAGATOR)
            CALL wpoints (ng, TILE, iNLM)
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
            CALL horz_mix (ng, TILE, iNLM)
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
            CALL set_nudgcof (ng, TILE, iNLM)
          END DO
        END DO
!$OMP END PARALLEL DO
#endif
      END IF
!
!=======================================================================
!  Initialize model state variables and forcing.  This part is
!  executed for each ensemble/perturbation/iteration run.
!=======================================================================

#ifdef TLM_CHECK
!
!  Clear state variables.
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1
            CALL initialize_ocean (ng, TILE, iNLM)
          END DO
        END DO
!$OMP END PARALLEL DO
#endif

#if defined SOLVE3D && !defined INI_FILE
!
!-----------------------------------------------------------------------
!  If analytical initial conditions, compute initial time-evolving
!  depths with zero free-surface.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile)                             &
!$OMP&            SHARED(ng,numthreads,Tindex,FirstPass)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*thread,subs*(thread+1)-1
          CALL set_depth (ng, TILE, Tindex, FirstPass)
        END DO
      END DO
!$OMP END PARALLEL DO
#endif
!
!-----------------------------------------------------------------------
!  Set primitive variables initial conditions.  Use analytical
!  functions or read from an initial or restart NetCDF file.
!-----------------------------------------------------------------------

#ifdef ANA_INITIAL
!
!  Analytical initial conditions for momentum and active tracers.
!
      IF (nrrec.eq.0) THEN
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1
            CALL ana_initial (ng, TILE, iNLM)
          END DO
        END DO
!$OMP END PARALLEL DO
      END IF
#endif

#if defined ANA_PASSIVE && defined SOLVE3D
!
!  Analytical initial conditions for inert passive tracers.
!
      IF (nrrec.eq.0) THEN
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1
            CALL ana_passive (ng, TILE, iNLM)
          END DO
        END DO
!$OMP END PARALLEL DO
      END IF
#endif

#if defined ANA_BIOLOGY && defined SOLVE3D
!
!  Analytical initial conditions for biology tracers.
!
      IF (nrrec.eq.0) THEN
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1
            CALL ana_biology (ng, TILE, iNLM)
          END DO
        END DO
!$OMP END PARALLEL DO
      END IF
#endif

#if defined ANA_SEDIMENT && defined SOLVE3D
!
!  Analytical initial conditions for sediment tracers.
!
      IF (nrrec.eq.0) THEN
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1
            CALL ana_sediment (ng, TILE, iNLM)
          END DO
        END DO
!$OMP END PARALLEL DO
      END IF
#endif
!
!  Read in initial conditions for initial or restart NetCDF file.
!
#ifdef INI_FILE
      IF (nrrec.eq.0) THEN
        CALL get_state (ng, iNLM, 1, INIname(ng), IniRec, Tindex)
      END IF
#else
      IF (nrrec.ne.0) THEN
        CALL get_state (ng, iNLM, 1, INIname(ng), IniRec, Tindex)
      END IF
#endif
      IF (exit_flag.ne.NoError) RETURN

#ifdef ANA_PERTURB
!
!  Perturb initial conditions.
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*thread,subs*(thread+1)-1
          CALL ana_perturb (ng, TILE, iNLM)
        END DO
      END DO
!$OMP END PARALLEL DO
#endif

#ifdef TLM_CHECK
!
!-----------------------------------------------------------------------
!  Add a perturbation to nonlinear state variable according to the outer
!  loop iteration with the steepest descent direction of the gradient
!  (adjoint state).
!-----------------------------------------------------------------------
!
      IF (outer.ge.1) THEN
!$OMP PARALLEL DO PRIVATE(thread,subs,tile,Tindex)                      &
!$OMP             SHARED(ng,numthreads,Lnew)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1
            CALL ini_perturb (ng, TILE, Lnew(ng), Tindex)
          END DO
        END DO
      END IF
!$OMP END PARALLEL DO
#endif

#ifdef IS4DVAR
!
!-----------------------------------------------------------------------
!  Set record counter for writing initial conditions into nonlinear
!  initial NetCDF file.
!-----------------------------------------------------------------------
!
!  Notice that initial conditions are written always into "tINIindx"
!  greater than one and the original initial conditions (usually, time
!  record one in NetCDF) are left unmodified for the backgound state
!  computations in the cost function. If cycling initial conditions
!  file, the adjusted initial conditions are written always into
!  record two and the counters "tINIindx" and "NrecINI" are reset to
!  unity. Otherwise, all the intermediate adjuted nonlinear solutions
!  are stored in the initial nonlinear NetCDF file.
!
      IF (LcycleINI(ng)) THEN
        tINIindx(ng)=1
        NrecINI(ng)=1
      END IF
      IF (Nrun.eq.ERstr) tINIindx(ng)=Tindex
!
!-----------------------------------------------------------------------
!  Add tangent linear increments to previous nonlinear initial
!  conditions.
!-----------------------------------------------------------------------
!
!  Add tangent linear increments to previous nonlinear initial
!  conditions.
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile,Tindex)                      &
!$OMP             SHARED(ng,numthreads,Lold)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*thread,subs*(thread+1)-1
          CALL ini_adjust (ng, TILE, Lold(ng), Tindex)
        END DO
      END DO
!$OMP END PARALLEL DO
!
!  Write out initial condition into nonlinear initial NetCDF file.
!
      CALL wrt_ini (ng, Tindex)
#endif

#ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Compute initial time-evolving depths.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile)                             &
!$OMP&            SHARED(ng,numthreads,Tindex,FirstPass)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*thread,subs*(thread+1)-1
          CALL set_depth (ng, TILE, Tindex, FirstPass)
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

#if defined FOUR_DVAR || !defined TANGENT || !defined ADJOINT
!
!-----------------------------------------------------------------------
!  Read in initial forcing, climatology and assimilation data from
!  input NetCDF files.  It loads the first relevant data record for
!  the time-interpolation between snapshots.
!-----------------------------------------------------------------------
!
      CALL get_data (ng)
      IF (exit_flag.ne.NoError) RETURN
#endif
!
!-----------------------------------------------------------------------
!  Compute grid stiffness.
!-----------------------------------------------------------------------
!
      IF (Lstiffness) THEN
        Lstiffness=.FALSE.
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1
            CALL stiffness (ng, TILE, iNLM)
          END DO
        END DO
!$OMP END PARALLEL DO
      END IF

#if defined FLOATS || defined STATIONS
!
!-----------------------------------------------------------------------
!  If applicable, convert initial locations to fractional grid
!  coordinates.
!-----------------------------------------------------------------------
!
      CALL grid_coords (ng, iNLM)
#endif

#ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off initiialization time wall clock.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread) SHARED(ng,numthreads)
      DO thread=0,numthreads-1
        CALL wclock_off (ng, iNLM, 1)
      END DO
!$OMP END PARALLEL DO
#endif
      RETURN
      END SUBROUTINE initial
