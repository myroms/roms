#include "cppdefs.h"
      SUBROUTINE initial (ng)
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2012 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
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
# ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti
# endif
#ifdef TLM_CHECK
      USE ini_adjust_mod, ONLY : ini_perturb
#endif
      USE ini_hmixcoef_mod, ONLY : ini_hmixcoef
      USE metrics_mod, ONLY : metrics
#ifdef SOLVE3D
      USE set_depth_mod, ONLY : set_depth
      USE omega_mod, ONLY : omega
      USE rho_eos_mod, ONLY : rho_eos
      USE set_massflux_mod, ONLY : set_massflux
#endif
#ifdef MASKING
      USE set_masks_mod, ONLY : set_masks
#endif
      USE stiffness_mod, ONLY : stiffness
#ifdef WAVES_OCEAN
      USE ocean_coupler_mod, ONLY : ocn2wav_coupling
#endif
#if defined PROPAGATOR || \
    (defined MASKING    && (defined READ_WATER || defined WRITE_WATER))
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
      logical, save :: First = .TRUE.
      logical :: update = .FALSE.

      integer :: LBi, UBi, LBj, UBj
      integer :: Fcount, IniRec, Tindex
      integer :: thread, tile
#if defined ADJUST_BOUNDARY || \
    defined ADJUST_STFLUX   || defined ADJUST_WSTRESS
      integer :: irec
#endif
!
!=======================================================================
!   Initialize model variables.
!=======================================================================
!
!$OMP MASTER
      IF (Master) THEN
#if defined PERTURBATION
        WRITE (stdout,10) Nrun
 10     FORMAT (/,' <<<< Ensemble/Perturbation Run: ',i5.5,' >>>>',/)
#elif defined IS4DVAR    || defined SENSITIVITY_4DVAR || \
      defined TL_W4DPSAS || defined TL_W4DVAR         || \
      defined W4DPSAS    || defined W4DVAR
        WRITE (stdout,10) outer, inner
 10     FORMAT (/,' <<<< 4D Variational Data Assimilation, ',           &
     &          'Outer = ',i3.3, ', Inner = ',i3.3,' >>>>',/)
#endif
        WRITE (stdout,20) 'INITIAL: Configuring and initializing ',     &
     &                    'forward nonlinear model ...'
 20     FORMAT (/,1x,a,a,/)
      END IF
!$OMP END MASTER
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
!
!$OMP MASTER
      synchro_flag(ng)=.TRUE.
      first_time(ng)=0
      tdays(ng)=dstart
      time(ng)=tdays(ng)*day2sec
      ntstart(ng)=INT((time(ng)-dstart*day2sec)/dt(ng))+1
      ntend(ng)=ntimes(ng)
      ntfirst(ng)=ntstart(ng)
      CALL time_string (time(ng), time_code(ng))
!$OMP END MASTER
!$OMP BARRIER

      IniRec=nrrec(ng)
      Tindex=1

      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)

#ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Start time wall clocks.
!-----------------------------------------------------------------------
!
      DO thread=THREAD_RANGE
        CALL wclock_on (ng, iNLM, 2)
      END DO
!$OMP BARRIER
#endif

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
!$OMP MASTER
      IF (First) THEN
        First=.FALSE.
# ifdef ANA_INITIAL
        LdefINI(ng)=.TRUE.
# endif
        LdefHIS(ng)=.TRUE.
        CALL def_ini (ng)
# ifdef DISTRIBUTE
        CALL mp_bcasti (ng, iNLM, exit_flag)
# endif
        IF (exit_flag.ne.NoError) RETURN
        IniRec=nrrec(ng)
        INI(ng)%Rindex=IniRec
      ELSE
        IniRec=INI(ng)%Rindex
      END IF

# ifdef ADJUST_BOUNDARY
!
!  Initialize open boundary counter for storage arrays.
!
      OBCcount(ng)=0
# endif
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
!
!  Initialize surface forcing counter for storage arrays.
!
      SFcount(ng)=0
# endif
!
!  Reset nonlinear history time record counters. These counters are
!  reset on every iteration pass. This file is created on the first
!  iteration pass.
!
      HIS(ng)%Rindex=0
      Fcount=HIS(ng)%Fcount
      HIS(ng)%Nrec(Fcount)=0

# ifdef IS4DVAR
!
!  Activate switches to writting data into average, history and
!  restart files.
!
      LwrtAVG(ng)=.TRUE.
      LwrtHIS(ng)=.TRUE.
      LwrtRST(ng)=.TRUE.
# endif
!$OMP END MASTER
!$OMP BARRIER
#endif
#ifdef OBSERVATIONS
!
!-----------------------------------------------------------------------
!  Open observations NetCDF file and initialize various variables
!  needed for processing the nonlinear state solution at observation
!  locations.
!-----------------------------------------------------------------------
!
!$OMP MASTER
      CALL obs_initial (ng, iNLM, .FALSE.)
!$OMP END MASTER
!$OMP BARRIER
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
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL ana_grid (ng, tile, iNLM)
# ifdef MASKING
          CALL ana_mask (ng, tile, iNLM)
# endif
        END DO
!$OMP BARRIER
#else
!$OMP MASTER
        CALL get_grid (ng, iNLM)
!$OMP END MASTER
!$OMP BARRIER
# ifdef DISTRIBUTE
        CALL mp_bcasti (ng, iNLM, exit_flag)
# endif
        if (exit_flag.ne.NoError) RETURN
#endif

#ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Set vertical S-coordinate transformation function.
!-----------------------------------------------------------------------
!
!$OMP MASTER
        CALL set_scoord (ng)
!$OMP END MASTER
!$OMP BARRIER
#endif

#ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Set barotropic time-steps average weighting function.
!-----------------------------------------------------------------------
!
!$OMP MASTER
        CALL set_weights (ng)
!$OMP END MASTER
!$OMP BARRIER
#endif
!
!-----------------------------------------------------------------------
!  Compute various metric term combinations.
!-----------------------------------------------------------------------
!
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL metrics (ng, tile, iNLM)
        END DO
!$OMP BARRIER
!
!-----------------------------------------------------------------------
!  If appropriate, set nudging coefficiests time scales.
!-----------------------------------------------------------------------
!
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL ana_nudgcoef (ng, tile, iNLM)
        END DO
!$OMP BARRIER
      END IF
!
!-----------------------------------------------------------------------
!  Initialize horizontal mixing coefficients.
!-----------------------------------------------------------------------
!
      DO tile=first_tile(ng),last_tile(ng),+1
        CALL ini_hmixcoef (ng, tile, iNLM)
      END DO
!$OMP BARRIER

#if defined VISC_GRID || defined DIFF_GRID || defined SPONGE
!
!-----------------------------------------------------------------------
!  Rescale horizontal mixing coefficients according to the local
!  grid size. If applicable, increases horizontal mixing in sponge
!  areas.
!-----------------------------------------------------------------------
!
      DO tile=first_tile(ng),last_tile(ng),+1
        CALL ana_hmixcoef (ng, tile, iNLM)
      END DO
!$OMP BARRIER
#endif
!
!=======================================================================
!  Initialize model state variables and forcing.  This part is
!  executed for each ensemble/perturbation/iteration run.
!=======================================================================

#ifdef TLM_CHECK
!
!  Clear state variables.
!
      DO tile=first_tile(ng),last_tile(ng),+1
        CALL initialize_ocean (ng, tile, iNLM)
      END DO
!$OMP BARRIER
#endif

#if defined SOLVE3D && !defined INI_FILE
!
!-----------------------------------------------------------------------
!  If analytical initial conditions, compute initial time-evolving
!  depths with zero free-surface.
!-----------------------------------------------------------------------
!
      DO tile=first_tile(ng),last_tile(ng),+1
        CALL set_depth (ng, tile)
      END DO
!$OMP BARRIER
#endif
!
!-----------------------------------------------------------------------
!  Set primitive variables initial conditions.
!-----------------------------------------------------------------------

#ifdef ANA_INITIAL
!
!  Analytical initial conditions for momentum and active tracers.
!
      IF (nrrec(ng).eq.0) THEN
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL ana_initial (ng, tile, iNLM)
        END DO
!$OMP BARRIER
      END IF
#endif

#if defined ANA_PASSIVE && defined SOLVE3D
!
!  Analytical initial conditions for inert passive tracers.
!
      IF (nrrec(ng).eq.0) THEN
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL ana_passive (ng, tile, iNLM)
        END DO
!$OMP BARRIER
      END IF
#endif

#if defined ANA_BIOLOGY && defined SOLVE3D
!
!  Analytical initial conditions for biology tracers.
!
      IF (nrrec(ng).eq.0) THEN
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL ana_biology (ng, tile, iNLM)
        END DO
!$OMP BARRIER
      END IF
#endif

#if defined ANA_SEDIMENT && defined SOLVE3D
!
!  Analytical initial conditions for sediment tracers.
!
      IF (nrrec(ng).eq.0) THEN
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL ana_sediment (ng, tile, iNLM)
        END DO
!$OMP BARRIER
      END IF
#endif

#ifdef INI_FILE
!
!  Read in initial conditions from initial NetCDF file.
!
!$OMP MASTER
      CALL get_state (ng, iNLM, 1, INI(ng)%name, IniRec, Tindex)
!$OMP END MASTER
!$OMP BARRIER
# ifdef DISTRIBUTE
      CALL mp_bcasti (ng, iNLM, exit_flag)
# endif
      IF (exit_flag.ne.NoError) RETURN
#else
!
!  If restart, read in initial conditions restart NetCDF file.
!
      IF (nrrec(ng).ne.0) THEN
!$OMP MASTER
        CALL get_state (ng, 0, 1, INI(ng)%name, IniRec, Tindex)
!$OMP END MASTER
!$OMP BARRIER
# ifdef DISTRIBUTE
        CALL mp_bcasti (ng, iNLM, exit_flag)
# endif
        IF (exit_flag.ne.NoError) RETURN
      END IF
#endif

#if (defined W4DPSAS || defined TL_W4DPSAS || \
     defined W4DPSAS_SENSITIVITY)          && \
    (defined ADJUST_BOUNDARY || defined ADJUST_WSTRESS ||\
     defined ADJUST_STFLUX)
!
!-----------------------------------------------------------------------
!  Read in the surface forcing and/or open boundary conditions
!  increments for PSAS from record IniRec of the NLM initial NetCDF
!  file.
!-----------------------------------------------------------------------
!
      IF (Nrun.gt.1) THEN
!$OMP MASTER
        CALL get_state (ng, 12, 12, INI(ng)%name, IniRec, Tindex)
!$OMP END MASTER
!$OMP BARRIER
# ifdef DISTRIBUTE
        CALL mp_bcasti (ng, iNLM, exit_flag)
# endif
        IF (exit_flag.ne.NoError) RETURN
      END IF
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
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL ini_perturb (ng, tile, Lnew(ng), Tindex)
        END DO
!$OMP BARRIER
      END IF
#endif

#ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Compute initial time-evolving depths.
!-----------------------------------------------------------------------
!
      DO tile=first_tile(ng),last_tile(ng),+1
        CALL set_depth (ng, tile)
      END DO
!$OMP BARRIER
!
!-----------------------------------------------------------------------
!  Compute initial horizontal mass fluxes, Hz*u/n and Hz*v/m.
!-----------------------------------------------------------------------
!
      DO tile=first_tile(ng),last_tile(ng),+1
        CALL set_massflux (ng, tile)
      END DO
!$OMP BARRIER
!
!-----------------------------------------------------------------------
!  Compute initial S-coordinates vertical velocity. Compute initial
!  density anomaly from potential temperature and salinity via equation
!  of state for seawater.  Also compute other equation of state related
!  quatities.
!-----------------------------------------------------------------------
!
      DO tile=first_tile(ng),last_tile(ng),+1
        CALL omega (ng, tile)
        CALL rho_eos (ng, tile)
      END DO
!$OMP BARRIER
#endif

#if defined ANA_PSOURCE && \
   (defined UV_PSOURCE  || defined TS_PSOURCE || defined Q_PSOURCE)
!
!-----------------------------------------------------------------------
!  Set point Sources/Sinks position, direction, special flag, and mass
!  transport nondimensional shape profile with analytcal expressions.
!  Point sources are at U- and V-points. We need to get their positions
!  to process internal Land/Sea masking arrays during initialization.
!-----------------------------------------------------------------------
!
      DO tile=first_tile(ng),last_tile(ng),+1
        CALL ana_psource (ng, tile, iNLM)
      END DO
!$OMP BARRIER
#endif

#if defined FOUR_DVAR || !defined TANGENT || !defined ADJOINT
!
!-----------------------------------------------------------------------
!  Read in initial forcing, climatology and assimilation data from
!  input NetCDF files.  It loads the first relevant data record for
!  the time-interpolation between snapshots.
!-----------------------------------------------------------------------
!
!$OMP MASTER

# ifdef ADJUST_BOUNDARY
!
!  If first pass of iteration loop, set time of open boundary
!  adjustment.
!
      IF (Nrun.eq.ERstr) THEN
        OBC_time(1,ng)=time(ng)
        DO irec=2,Nbrec(ng)
          OBC_time(irec,ng)=OBC_time(irec-1,ng)+nOBC(ng)*dt(ng)
        END DO
      END IF
# endif
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
!
!  If first pass of iteration loop, set time of surface forcing
!  adjustment.
!
      IF (Nrun.eq.ERstr) THEN
        SF_time(1,ng)=time(ng)
        DO irec=2,Nfrec(ng)
          SF_time(irec,ng)=SF_time(irec-1,ng)+nSFF(ng)*dt(ng)
        END DO
      END IF
# endif
!
!  If applicable, read in input data.
!
      CALL close_inp (ng, iNLM)
# ifdef TIMELESS_DATA
      CALL get_idata (ng)
# endif
      CALL get_data (ng)
# ifdef DISTRIBUTE
      CALL mp_bcasti (ng, iNLM, exit_flag)
# endif
!$OMP END MASTER
!$OMP BARRIER
      IF (exit_flag.ne.NoError) RETURN
#endif

#ifdef MASKING
!
!-----------------------------------------------------------------------
!  Set internal I/O mask arrays.
!-----------------------------------------------------------------------
!
      DO tile=first_tile(ng),last_tile(ng),+1
        CALL set_masks (ng, tile, iNLM)
      END DO
!$OMP BARRIER
#endif

#if defined PROPAGATOR || \
   (defined MASKING    && (defined READ_WATER || defined WRITE_WATER ))
!
!-----------------------------------------------------------------------
!  Set variables associated with the processing water points and/or
!  size of packed state arrays.
!-----------------------------------------------------------------------
!
      DO tile=first_tile(ng),last_tile(ng),+1
        CALL wpoints (ng, tile, iNLM)
      END DO
!$OMP BARRIER
#endif

#if defined NLM_OUTER || defined TL_W4DPSAS          || \
    defined W4DPSAS   || defined W4DPSAS_SENSITIVITY
!
!-----------------------------------------------------------------------
!  Read in convolved adjoint impulse forcing (first record) and its
!  application time.
!-----------------------------------------------------------------------
!
      IF (SporadicImpulse(ng)) THEN
        FrcRec(ng)=1
!$OMP MASTER
        CALL get_state (ng, 7, 7, TLF(ng)%name, FrcRec(ng), 1)
!$OMP END MASTER
!$OMP BARRIER
# ifdef DISTRIBUTE
        CALL mp_bcasti (ng, iTLM, exit_flag)
# endif
        IF (exit_flag.ne.NoError) RETURN
      END IF
#endif

#if defined ANA_DRAG && defined UV_DRAG_GRID
!
!-----------------------------------------------------------------------
!  Set analytical spatially varying bottom friction parameter.
!-----------------------------------------------------------------------
!
      IF (Nrun.eq.ERstr) THEN
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL ana_drag (ng, tile, iNLM)
        END DO
!$OMP BARRIER
      END IF
#endif
!
!-----------------------------------------------------------------------
!  Compute grid stiffness.
!-----------------------------------------------------------------------
!
      IF (Lstiffness) THEN
        Lstiffness=.FALSE.
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL stiffness (ng, tile, iNLM)
        END DO
!$OMP BARRIER
      END IF

#if defined FLOATS || defined STATIONS
!
!-----------------------------------------------------------------------
!  If applicable, convert initial locations to fractional grid
!  coordinates.
!-----------------------------------------------------------------------
!
!$OMP MASTER
      CALL grid_coords (ng, iNLM)
!$OMP END MASTER
!$OMP BARRIER
#endif

# ifdef WAVES_OCEAN
!
!-----------------------------------------------------------------------
!  Read in initial forcing from coupled wave model.
!-----------------------------------------------------------------------
!
      DO tile=first_tile(ng),last_tile(ng),+1
        CALL ocn2wav_coupling (ng, tile)
      END DO
!$OMP BARRIER
      IF (Master) WRITE (stdout,'(/)')
# endif
!
!-----------------------------------------------------------------------
!  Initialize time-stepping counter and clock.
!-----------------------------------------------------------------------
!
!  Subsract one time unit to avoid special case due to initialization
!  in the main time-stepping routine.
!
      iic(ng)=ntstart(ng)-1
!$OMP MASTER
      time(ng)=time(ng)-dt(ng)
!$OMP END MASTER
!$OMP BARRIER

#ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off initialization time wall clock.
!-----------------------------------------------------------------------
!
      DO thread=THREAD_RANGE
        CALL wclock_off (ng, iNLM, 2)
      END DO
!$OMP BARRIER
#endif
      RETURN
      END SUBROUTINE initial
