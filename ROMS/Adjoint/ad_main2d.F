#include "cppdefs.h"
#if defined ADJOINT && !defined SOLVE3D
      SUBROUTINE ad_main2d (RunInterval)
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2024 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  This routine is the main driver for ROMS/TOMS adjoint model (ADM)   !
!  when configured as a 2D barotropic shallow water ocean model. It    !
!  advances backward the ADM for all nested grids, if any, by the      !
!  specified time interval (seconds), RunInterval.                     !
!                                                                      !
# if defined STEP2D_FB_LF_AM3
!  Numerical 2D time-stepping kernel: FB AB3-AM4                       !
# elif defined STEP2D_FB_LF_AM3
!  Numerical 2D time-stepping kernel: FB LF-AM3                        !
# else
!  Numerical 2D time-stepping kernel: LF-AM3 (Legacy scheme)           !
# endif
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
# if defined MODEL_COUPLING && defined MCT_LIB
      USE mod_coupler
# endif
# ifdef FOUR_DVAR
      USE mod_fourdvar
# endif
      USE mod_iounits
# ifdef SENSITIVITY_4DVAR
      USE mod_ncparam
# endif
# ifdef NESTING
      USE mod_nesting
# endif
      USE mod_scalars
      USE mod_stepping
# ifdef SO_SEMI
      USE mod_storage
# endif
!
# if defined AD_SENSITIVITY   || defined I4DVAR_ANA_SENSITIVITY || \
     defined OPT_OBSERVATIONS || defined SENSITIVITY_4DVAR
      USE adsen_force_mod,         ONLY : adsen_force
# endif
      USE ad_diag_mod,             ONLY : ad_diag
# ifdef WEAK_CONSTRAINT
      USE ad_force_dual_mod,       ONLY : ad_force_dual
# endif
# ifdef FORCING_SV
      USE ad_forcing_mod,          ONLY : ad_forcing
# endif
# ifdef ADJUST_WSTRESS
      USE ad_frc_adjust_mod,       ONLY : ad_frc_adjust
# endif
      USE ad_ini_fields_mod,       ONLY : ad_ini_fields, ad_ini_zeta
# ifdef AD_OUTPUT_STATE
      USE ad_ini_fields_mod,       ONLY : ad_out_fields, ad_out_zeta
# endif
# if defined FOUR_DVAR && defined OBSERVATIONS
#  ifdef WEAK_CONSTRAINT
      USE ad_htobs_mod,            ONLY : ad_htobs
#  else
      USE ad_misfit_mod,           ONLY : ad_misfit
#  endif
# endif
# ifdef NESTING
      USE nesting_mod,             ONLY : nesting
      USE ad_nesting_mod,          ONLY : ad_nesting
#  ifndef ONE_WAY
      USE nesting_mod,             ONLY : do_twoway
#  endif
# endif
# ifdef ADJUST_BOUNDARY
      USE ad_obc_adjust_mod,       ONLY : ad_obc_adjust
# endif
# ifdef WEC_MELLOR_NOT_YET
!!    USE ad_radiation_stress_mod, ONLY : ad_radiation_stress
# endif
# ifdef AD_AVERAGES
      USE ad_set_avg_mod,          ONLY : ad_set_avg
# endif
# if defined SSH_TIDES_NOT_YET || defined UV_TIDES_NOT_YET
!!    USE ad_set_tides_mod,        ONLY : ad_set_tides
# endif
      USE ad_set_vbc_mod,          ONLY : ad_set_vbc
      USE ad_step2d_mod,           ONLY : ad_step2d
# ifdef FLOATS_NOT_YET
!!    USE ad_step_floats_mod,      ONLY : ad_step_floats
# endif
      USE dateclock_mod,           ONLY : time_string
# ifdef TIDE_GENERATING_FORCES
      USE equilibrium_tide_mod,    ONLY : equilibrium_tide
# endif
# ifdef WEAK_CONSTRAINT
      USE frc_weak_mod,            ONLY : frc_ADgather, frc_clear
# endif
# if defined ATM_COUPLING_NOT_YET && defined MCT_LIB
      USE mct_coupler_mod,         ONLY : ocn2atm_coupling
# endif
# if defined WAV_COUPLING_NOT_YET && defined MCT_LIB
      USE mct_coupler_mod,         ONLY : ocn2wav_coupling
# endif
# if (defined FOUR_DVAR    && !defined I4DVAR_ANA_SENSITIVITY) && \
      defined OBSERVATIONS
      USE obs_read_mod,            ONLY : obs_read
# endif
# ifdef SO_SEMI
      USE packing_mod,             ONLY : ad_pack
# endif
      USE strings_mod,             ONLY : FoundError
!
      implicit none
!
!  Imported variable declarations.
!
      real(dp), intent(in) :: RunInterval
!
!  Local variable declarations.
!
      logical :: backward = .TRUE.
      logical :: ad_advance
      logical :: DoNestLayer, Time_Step
      integer :: Nsteps, Rsteps
      integer :: ig, il, is, istep, ng, nl, tile
      integer :: my_iif, next_indx1
# ifdef FLOATS_NOT_YET
      integer :: Lend, Lstr, chunk_size
# endif
      integer :: ks, kt
# ifdef NESTING
      integer :: icount, itcount
# endif
!
      real(r8) :: HalfDT, MaxDT, my_StepTime
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
!=======================================================================
!  Time-step adjoint 2D primitive equations backwards.
!=======================================================================

# ifdef NESTING
!
!  Compute nonlinear model timestepping sequence and load it into
!  StepInfo array.  The sequence of values are reversed to timestep
!  the adjoint model backwards.
!
      CALL nlm_step_sequence (RunInterval, icount)
      IF (FoundError(exit_flag, NoError, __LINE__,                      &
     &               __FILE__)) RETURN
# endif
!
!  Initialize.
!
      Time_Step=.TRUE.
      DoNestLayer=.TRUE.
# ifdef NESTING
      itcount=icount+1
# endif
!
      KERNEL_LOOP : DO WHILE (Time_Step)
!
!  In nesting applications, the number of nesting layers (NestLayers) is
!  used to facilitate refinement grids and composite/refinament grids
!  combinations. Otherwise, the solution it is looped once for a single
!  grid application (NestLayers = 1).
!
        nl=0
# ifdef NESTING
        DO ng=1,Ngrids
          TwoWayInterval(ng)=dt(1)
        END DO
# endif
!
        NEST_LAYER : DO WHILE (DoNestLayer)

# ifdef NESTING
!
!  Determine current nested layer (nl), number of timesteps (Nsteps)
!  to execute for grids in current nested layer, number of timesteps
!  (Rsteps) to complete the RunInterval time window, and timestep
!  counter (step_counter) for each grid. Their values are reversed
!  from the saved nonlinear model time stepping sequence.
!
          itcount=itcount-1             ! count backwards for adjoint
          IF (itcount.gt.0) THEN
            nl=StepInfo(itcount,1)
            Nsteps=StepInfo(itcount,2)
            IF (itcount.eq.icount) Nsteps=1
            Rsteps=StepInfo(itcount,3)
            DO ng=1,Ngrids
              step_counter(ng)=StepInfo(itcount,ng+3)
            END DO
          ELSE
            nl=0
          END IF
# else
!
!  Determine number of time steps to compute in each nested grid layer
!  based on the specified time interval (seconds), RunInterval. Non
!  nesting applications have NestLayers=1. Notice that RunInterval is
!  set in the calling driver. Its value may span the full period of the
!  simulation, a multi-model coupling interval (RunInterval > ifac*dt),
!  or just a single step (RunInterval=0).
!
          CALL ntimesteps (iADM, RunInterval, nl, Nsteps, Rsteps)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
!
          IF ((nl.le.0).or.(nl.gt.NestLayers)) EXIT
!
!  Time-step governing equations for Nsteps.
!
          STEP_LOOP : DO istep=Nsteps,1,-1
!
!  Set time indices and time clock.
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
# ifdef JEDI
              jic(ng)=jic(ng)-1
              time4jedi(ng)=time4jedi(ng)-dt(ng)
# endif
              tdays(ng)=time(ng)*sec2day
              IF (step_counter(ng).eq.1) Time_Step=.FALSE.
              IF ((iic(ng).eq.ntstart(ng)).and.(ng.eq.Ngrids)) THEN
                ad_advance=.FALSE.
# ifdef NESTING
                DO ng=1,Ngrids
                  TwoWayInterval(ng)=0.0_r8
                END DO
# endif
              ELSE
                ad_advance=.TRUE.
              END IF
            END DO
!
!-----------------------------------------------------------------------
!  Read in required data, if any, from input NetCDF files.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              CALL ad_get_data (ng)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END DO
!
!-----------------------------------------------------------------------
!  Process input data, if any: time interpolate between snapshots.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              DO tile=first_tile(ng),last_tile(ng),+1
                CALL ad_set_data (ng, tile)
# ifdef AD_AVERAGES
                CALL ad_set_avg (ng, tile)
# endif
# ifdef TIDE_GENERATING_FORCES
                CALL equilibrium_tide (ng, tile, iADM)
# endif
# ifdef DIAGNOSTICS
!!              CALL ad_set_diags (ng, tile)
# endif
              END DO
            END DO
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Avoid time-stepping if additional delayed IO time-step.
!
            ADVANCE : IF (ad_advance) THEN

# ifdef FLOATS_NOT_YET
!
!-----------------------------------------------------------------------
!  Compute Lagrangian drifters trajectories: Split all the drifters
!  between all the computational threads, except in distributed-memory
!  and serial configurations. In distributed-memory, the parallel node
!  containing the drifter is selected internally since the state
!  variables do not have a global scope.
!-----------------------------------------------------------------------
!
!  Shift floats time indices.
!
              DO ig=1,GridsInLayer(nl)
                ng=GridNumber(ig,nl)
                IF (Lfloats(Ng)) THEN
                  nfp1(ng)=MOD(nfp1(ng)+1,NFT+1)
                  nf  (ng)=MOD(nf  (ng)+1,NFT+1)
                  nfm1(ng)=MOD(nfm1(ng)+1,NFT+1)
                  nfm2(ng)=MOD(nfm2(ng)+1,NFT+1)
                  nfm3(ng)=MOD(nfm3(ng)+1,NFT+1)
!
#  ifdef _OPENMP
                  chunk_size=(Nfloats(ng)+numthreads-1)/numthreads
                  Lstr=1+MyThread*chunk_size
                  Lend=MIN(Nfloats(ng),Lstr+chunk_size-1)
#  else
                  Lstr=1
                  Lend=Nfloats(ng)
#  endif
                  CALL ad_step_floats (ng, Lstr, Lend)
                END IF
              END DO
# endif

# ifdef NESTING
!
!-----------------------------------------------------------------------
!  If donor to a finer grid, extract data for the external contact
!  points. This is the latest solution for the coarser grid.
!
!  It is stored in the REFINED structure so it can be used for the
!  space-time interpolation when "nputD" argument is used above.
!-----------------------------------------------------------------------
!
              DO ig=1,GridsInLayer(nl)
                ng=GridNumber(ig,nl)
                IF (DonorToFiner(ng)) THEN
                  CALL ad_nesting (ng, iADM, ngetD)
                END IF
              END DO

#  ifndef ONE_WAY
!
!-----------------------------------------------------------------------
!  If refinement grids, perform two-way coupling between fine and
!  coarse grids. Correct coarse grid tracers values at the refinement
!  grid with refined accumulated fluxes.  Then, replace coarse grid
!  state variable with averaged refined grid values (two-way nesting).
!  Update coarse grid depth variables.
!
!  The two-way exchange of infomation between nested grids needs to be
!  done at the correct time-step and in the right sequence.
!-----------------------------------------------------------------------
!
              DO il=1,NestLayers
                DO ig=1,GridsInLayer(il)
                  ng=GridNumber(ig,il)
                  IF (do_twoway(iADM,nl,il,ng,istep)) THEN
                    CALL ad_nesting (ng, iADM, n2way)
                  END IF
                END DO
              END DO
#  endif
# endif

# ifdef NESTING
!
!
!  If composite or mosaic grids, process additional points in the
!  contact zone between connected grids for the state variables
!  associated with the 2D engine CORRECTOR STEP section (KNEW INDEX).
!
              DO ig=1,GridsInLayer(nl)
                ng=GridNumber(ig,nl)
                IF (ANY(CompositeGrid(:,ng))) THEN
                  CALL ad_nesting (ng, iADM, n2dCS)
                END IF
              END DO

#  if defined MASKING && defined WET_DRY
!
!  If nesting and wetting and drying, scale horizontal interpolation
!  weights to account for land/sea masking in contact areas. This needs
!  to be done at very time-step since the Land/Sea masking is time
!  dependent.
!
              DO ig=1,GridsInLayer(nl)
                ng=GridNumber(ig,nl)
                CALL ad_nesting (ng, iADM, nmask)
              END DO
#  endif
# endif

# ifdef STEP2D_FB_AB3_AM4
!
!-----------------------------------------------------------------------
!  Solve adjoint vertically integrated primitive equations for the
!  free-surface and barotropic momentum components using a generalized
!  Forward-Backward, 3rd-order Adams-Bashforth / 4th-order Adams-Moulton
!  (FB AB3-AM4) time stepping scheme (Shchepetkin and McWilliams, 2009).
!-----------------------------------------------------------------------
!
              DO ig=1,GridsInLayer(nl)
                ng=GridNumber(ig,nl)
                iif(ng)=1
                nfast(ng)=1

                IF (MOD(knew(ng),2).eq.0) THEN              ! zig-zag
                  DO tile=first_tile(ng),last_tile(ng),+1   ! processing
                    CALL ad_step2d (ng, tile)               ! sequence
                  END DO
                ELSE
                  DO tile=last_tile(ng),first_tile(ng),-1
                    CALL ad_step2d (ng, tile)
                  END DO
                END IF

!^              kstp(ng)=knew(ng)
!^              knew(ng)=kstp(ng)+1
!^              IF (knew(ng).gt.4) knew(ng)=1
!^
                knew(ng)=kstp(ng)
                kstp(ng)=knew(ng)-1
                IF (kstp(ng).lt.1) kstp(ng)=4
              END DO

# else

#  ifdef STEP2D_FB_LF_AM3
!
!-----------------------------------------------------------------------
!  Solve adjoint vertically integrated primitive equations for
!  free-surface and momentum components using a predictor-corrector
!  LeapFrog / 3rd-order Adams-Moulton with a Forward-Backward
!  feeback (FB LF-AM3) stepping scheme (Shchepetkin and McWilliams,
!  2009).
!-----------------------------------------------------------------------
!
!  Corrector AM3 substep with FB-feedback.
!
              DO ig=1,GridsInLayer(nl)
                ng=GridNumber(ig,nl)
                iif(ng)=1

                DO tile=first_tile(ng),last_tile(ng),+1
                  CALL tl_step2d (ng, tile)
                END DO

!^              knew(ng)=3-kstp(ng)
!^              next_kstp(ng)=knew(ng)
!^
                knew(ng)=next_kstp(ng)
                kstp(ng)=3-knew(ng)
              END DO
!
!  Predictor LF substep with FB-feedback.
!
              DO ig=1,GridsInLayer(nl)
                ng=GridNumber(ig,nl)
                iif(ng)=1

                DO tile=last_tile(ng),first_tile(ng),-1
                  CALL ad_step2d (ng, tile)
                END DO

!^              kstp(ng)=next_kstp(ng)
!^              knew(ng)=3
!^
!^ HGA/AMM: Need to compute the adjoint of the indices

              END DO

#  else
!
!-----------------------------------------------------------------------
!  Solve adjoint vertically integrated primitive equations for
!  free-surface and momentum components using a predictor-corrector
!  LeapFrog with 3rd-order Adams-Moulton (LF-AM3) time stepping scheme
!  (ROMS legacy 2D kernel).
!-----------------------------------------------------------------------
!
!  Corrector step - Apply 2D adjoint time-step corrector scheme.  Notice
!  ==============   that there is not need for a corrector step during
!  the auxiliary (nfast+1) time-step.
!
              DO ig=1,GridsInLayer(nl)
                ng=GridNumber(ig,nl)
                iif(ng)=1
                nfast(ng)=1

                IF (iif(ng).lt.(nfast(ng)+1)) THEN
                  DO tile=first_tile(ng),last_tile(ng),+1
                    CALL ad_step2d (ng, tile)
                  END DO
                END IF
!
!  Set time indices for adjoint predictor step.
!
                next_indx1=3-indx1(ng)
                IF (.not.PREDICTOR_2D_STEP(ng)) THEN
                  PREDICTOR_2D_STEP(ng)=.TRUE.
!^                knew(ng)=next_indx1
!^                kstp(ng)=3-knew(ng)
!^                krhs(ng)=3
!^
                  kt=knew(ng)
                  ks=kstp(ng)
                  knew(ng)=krhs(ng)
                  kstp(ng)=kt
                  krhs(ng)=ks
!^                IF (my_iif.lt.(nfast(ng)+1)) indx1(ng)=next_indx1
                END IF
              END DO

#   ifdef NESTING
!
!  If composite or mosaic grids, process additional points in the
!  contact zone between connected grids for the state variables
!  associated with the 2D engine PREDICTOR STEP section.
!
              DO ig=1,GridsInLayer(nl)
                ng=GridNumber(ig,nl)
                IF (ANY(CompositeGrid(:,ng))) THEN
                  CALL ad_nesting (ng, iADM, n2dPS)
                END IF
              END DO
#   endif
!
!  Predictor step - Advance adjoint barotropic equations using 2D
!  ==============   time-step predictor scheme.  No actual time-
!  stepping is performed during the auxiliary (nfast+1) time-step.
!  It is needed to finalize the fast-time averaging of 2D fields,
!  if any, and compute the new time-evolving depths.
!
              DO ig=1,GridsInLayer(nl)
                ng=GridNumber(ig,nl)
                DO tile=last_tile(ng),first_tile(ng),-1
                  CALL ad_step2d (ng, tile)
                END DO
!
!  Set time indices for next adjoint corrector step. The
!  PREDICTOR_2D_STEP switch it is assumed to be false before the
!  first time-step.
!
                IF (PREDICTOR_2D_STEP(ng).and.                          &
     &              my_iif.le.(nfast(ng)+1)) THEN
                  PREDICTOR_2D_STEP(ng)=.FALSE.
!^                IF (FIRST_2D_STEP) THEN
!^                  kstp(ng)=indx1(ng)
!^                ELSE
!^                  kstp(ng)=3-indx1(ng)
!^                END IF
!^                knew(ng)=3
!^                krhs(ng)=indx1(ng)
!^
                  ks=knew(ng)
                  knew(ng)=krhs(ng)
                  krhs(ng)=ks
                END IF
              END DO
#  endif

# endif
            END IF ADVANCE

# if (defined FOUR_DVAR    && !defined I4DVAR_ANA_SENSITIVITY) && \
      defined OBSERVATIONS
!
!-----------------------------------------------------------------------
!  If appropriate, read observation and model state at observation
!  locations.  Then, compute adjoint forcing terms due to observations.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
#  ifdef SENSITIVITY_4DVAR
#   ifdef RBL4DVAR_FCT_SENSITIVITY
              IF (.not.Lsen4DVAR(ng).and.LsenFCT(ng)) THEN
#   else
              IF (.not.Lsen4DVAR(ng)) THEN
#   endif
#  endif
              HalfDT=0.5_r8*dt(ng)
              IF (((time(ng)-HalfDT).le.ObsTime(ng)).and.               &
     &             (ObsTime(ng).lt.(time(ng)+HalfDT))) THEN
                ProcessObs(ng)=.TRUE.
                CALL obs_read (ng, iADM, backward)
              ELSE
                ProcessObs(ng)=.FALSE.
              END IF
!
#  ifdef SP4DVAR
!
!  Skip assimilation of obs on first timestep unless inter=Nsaddle.
!
              IF ((iic(ng).ne.ntstart(ng)).or.Lsadd(ng)) THEN
#  endif

              DO tile=first_tile(ng),last_tile(ng),+1
#  ifdef WEAK_CONSTRAINT
                CALL ad_htobs (ng, tile, iADM)
#  else
                CALL ad_misfit (ng, tile, iADM)
#  endif
              END DO
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
#  ifdef SENSITIVITY_4DVAR
              END IF
#  endif
#  ifdef SP4DVAR
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
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              IF (ProcessObs(ng)) THEN
                DO tile=first_tile(ng),last_tile(ng),+1
                  CALL ad_force_dual (ng, tile, kstp(ng), nstp(ng))
                END DO
              END IF
            END DO
# endif

# ifdef NESTING
!
!-----------------------------------------------------------------------
!  If refinement grid, interpolate (space, time) state variables
!  contact points from donor grid extracted data.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              IF (RefinedGrid(ng).and.(RefineScale(ng).gt.0)) THEN
                CALL ad_nesting (ng, iADM, nputD)
              END IF
            END DO
# endif

# ifdef NESTING
!
!  If composite or mosaic grids, process additional points in the
!  contact zone between connected grids for 3D kernel free-surface.
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              IF (ANY(CompositeGrid(:,ng))) THEN
                CALL ad_nesting (ng, iADM, nzeta)
              END IF
            END DO
# endif

# ifdef SO_SEMI
!
!-----------------------------------------------------------------------
!  If stochastic optimals with respect the seminorm of chosen
!  functional, pack adjoint state surface forcing needed by the
!  dynamical propagator.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              IF (MOD(iic(ng)-1,nADJ(ng)).eq.0) THEN
                SOrec(ng)=SOrec(ng)+1
                DO tile=first_tile(ng),last_tile(ng),+1
                  CALL ad_pack (ng, tile, Nstr(ng), Nend(ng),           &
     &                          STORAGE(ng)%so_state(:,SOrec(ng)))
                END DO
              END IF
            END DO
# endif

# ifdef NESTING
!
!  If composite or mosaic grids, process additional points in the
!  contact zone between connected grids for bottom stress variables.
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              IF (ANY(CompositeGrid(:,ng))) THEN
                CALL ad_nesting (ng, iADM, nbstr)
              END IF
            END DO
# endif
!
!-----------------------------------------------------------------------
!  Set adjoint fields for vertical boundary conditions. Process tidal
!  forcing, if any.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              DO tile=first_tile(ng),last_tile(ng),+1
# if defined SSH_TIDES_NOT_YET || defined UV_TIDES_NOT_YET
                CALL ad_set_tides (ng, tile)
# endif
                CALL ad_set_vbc (ng, tile)
              END DO
            END DO

# ifdef WEC_MELLOR_NOT_YET
!
!-----------------------------------------------------------------------
!  Compute radiation stress terms.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              DO tile=last_tile(ng),first_tile(ng),-1
                CALL ad_radiation_stress (ng, tile)
              END DO
            END DO
# endif

# if defined WAV_COUPLING_NOT_YET && defined MCT_LIB
!
!-----------------------------------------------------------------------
!  Couple to waves model every CoupleSteps(Iwaves) timesteps: get
!  waves/sea fluxes.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              IF ((iic(ng).ne.ntstart(ng)).and.                         &
     &            MOD(iic(ng)-1,CoupleSteps(Iwaves,ng)).eq.0) THEN
                DO tile=first_tile(ng),last_tile(ng),+1
                  CALL ocn2wav_coupling (ng, tile)
                END DO
              END IF
            END DO
# endif

# ifdef ADJUST_WSTRESS
!
!-----------------------------------------------------------------------
!  Interpolate surface forcing increments and adjust surface forcing.
!  Skip first timestep.
!-----------------------------------------------------------------------
!
#  ifdef RBL4DVAR_FCT_SENSITIVITY
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              IF (.not.Lsen4DVAR(ng)) THEN          ! ignore in forecast
                IF (iic(ng).ne.ntstart(ng)) THEN    !           interval
                  DO tile=first_tile(ng),last_tile(ng),+1
                    CALL ad_frc_adjust (ng, tile, Lfout(ng))
                  END DO
                END IF
              END IF
            END DO
#  else
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              IF (iic(ng).ne.ntstart(ng)) THEN
                DO tile=first_tile(ng),last_tile(ng),+1
                  CALL ad_frc_adjust (ng, tile, Lfout(ng))
                END DO
              END IF
            END DO
#  endif
# endif

# ifdef ADJUST_BOUNDARY
!
!-----------------------------------------------------------------------
!  Interpolate open boundary increments and adjust open boundaries.
!  Skip first timestep.
!-----------------------------------------------------------------------
!
#  ifdef RBL4DVAR_FCT_SENSITIVITY
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              IF (.not.Lsen4DVAR(ng)) THEN          ! ignore in forecast
                IF (iic(ng).ne.ntstart(ng)) THEN    !           interval
                  DO tile=first_tile(ng),last_tile(ng),+1
                    CALL ad_obc_adjust (ng, tile, Lbout(ng))
                  END DO
                END IF
              END IF
            END DO
#  else
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              IF (iic(ng).ne.ntstart(ng)) THEN
                DO tile=first_tile(ng),last_tile(ng),+1
                  CALL ad_obc_adjust (ng, tile, Lbout(ng))
                END DO
              END IF
            END DO
#  endif
# endif

# if defined ATM_COUPLING_NOT_YET && defined MCT_LIB
!
!-----------------------------------------------------------------------
!  Couple to atmospheric model every CoupleSteps(Iatmos) timesteps: get
!  air/sea fluxes.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              IF ((iic(ng).ne.ntstart(ng)).and.                         &
     &            MOD(iic(ng)-1,CoupleSteps(Iatmos,ng)).eq.0) THEN
                DO tile=last_tile(ng),first_tile(ng),-1
                  CALL ocn2atm_coupling (ng, tile)
                END DO
              END IF
            END DO
# endif
!
!-----------------------------------------------------------------------
!  Compute and report diagnostics. If appropriate, accumulate time-
!  averaged output data which needs a irreversible loop in shared-memory
!  jobs.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              DO tile=first_tile(ng),last_tile(ng),+1     ! irreversible
                CALL ad_diag (ng, tile)
# ifdef AD_AVERAGES
                CALL ad_set_avg (ng, tile)
# endif
              END DO
            END DO
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef FORCING_SV
!
!-----------------------------------------------------------------------
!  Compute the adjoint forcing for the forcing singular vectors.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              DO tile=first_tile(ng),last_tile(ng),+1
                CALL ad_forcing (ng, tile, kstp(ng), nstp(ng))
              END DO
            END DO
# endif

# ifdef AD_OUTPUT_STATE
!
!-----------------------------------------------------------------------
!  Set full adjoint output arrays. Due to the exact discrete adjoint,
!  the predictor/corrector time-stepping scheme with multiple time
!  levels, pieces of the adjoint solution are in two-time levels and
!  need to be added in the "_sol" arrays for output purposes.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              IF (iic(ng).ne.ntend(ng)) THEN
                DO tile=last_tile(ng),first_tile(ng),-1
                  CALL ad_out_fields (ng, tile, iADM)
                END DO
                DO tile=first_tile(ng),last_tile(ng),+1
                  CALL ad_out_zeta (ng, tile, iADM)
                END DO
              END IF
            END DO
# endif

# ifndef FORCING_SV
!
!-----------------------------------------------------------------------
!  If not a restart, initialize all time levels and compute other
!  initial fields.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              IF (iic(ng).eq.ntend(ng)) THEN
#  ifdef NESTING
!
!  Extract donor grid initial data at contact points and store it in
!  REFINED structure so it can be used for the space-time interpolation.
!
                IF (RefinedGrid(ng)) THEN
                  CALL ad_nesting (ng, iADM, ngetD)
                END IF
#  endif
!
!  Adjoint of initialize other state variables.
!
                DO tile=last_tile(ng),first_tile(ng),-1
                  CALL ad_ini_fields (ng, tile, iADM)
                END DO
!
!  Adjoint of initialize free-surface and compute initial level
!  thicknesses and depths.
!
                DO tile=first_tile(ng),last_tile(ng),+1
                  CALL ad_ini_zeta (ng, tile, iADM)
                END DO
              END IF
            END DO
# endif

# if defined WEAK_CONSTRAINT && !defined SP4DVAR
!
!-----------------------------------------------------------------------
!  Gather weak constraint forcing to storage arrays.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              IF (iic(ng).ne.ntstart(ng)) THEN
                DO tile=first_tile(ng),last_tile(ng),+1
                  CALL frc_ADgather (ng, tile)
                END DO
              END IF
            END DO
# endif
!
!-----------------------------------------------------------------------
!  If appropriate, write out fields into output NetCDF files.  Notice
!  that IO data is written in delayed and serial mode.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              CALL ad_output (ng)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END DO

# if defined WEAK_CONSTRAINT && !defined SP4DVAR
!
!-----------------------------------------------------------------------
!  Copy storage arrays index 1 into index 2, and clear index 1.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              IF (MOD(iic(ng)-1,nADJ(ng)).eq.0) THEN
                DO tile=first_tile(ng),last_tile(ng),+1
                  CALL frc_clear (ng, tile)
                END DO
              END IF
            END DO
# endif

# if (defined AD_SENSITIVITY   || defined I4DVAR_ANA_SENSITIVITY || \
      defined OPT_OBSERVATIONS || defined SENSITIVITY_4DVAR)     && \
     !defined OBS_SPACE
!
!-----------------------------------------------------------------------
!  Add appropriate forcing terms to the adjoint model. The form of the
!  forcing depends on the functional whose sensitivity is required.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
#  ifdef SENSITIVITY_4DVAR
              IF (Lsen4DVAR(ng)) THEN
#  endif
#  if !defined AD_IMPULSE
                IF ((DendS(ng).ge.tdays(ng)).and.                       &
     &              (tdays(ng).ge.DstrS(ng))) THEN
#  endif
                  DO tile=first_tile(ng),last_tile(ng),+1
                    CALL adsen_force (ng, tile)
                  END DO
#  if !defined AD_IMPULSE
                END IF
#  endif
#  ifdef SENSITIVITY_4DVAR
              END IF
#  endif
            END DO
# endif
!
!-----------------------------------------------------------------------
!  Set adjoint time indices and time clock.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              iic(ng)=iic(ng)-1
# ifndef NESTING
              step_counter(ng)=step_counter(ng)-1
# endif
              time(ng)=time(ng)-dt(ng)
              CALL time_string (time(ng), time_code(ng))
            END DO

          END DO STEP_LOOP

        END DO NEST_LAYER

      END DO KERNEL_LOOP
!
      RETURN
      END SUBROUTINE ad_main2d
#else
      SUBROUTINE ad_main2d
      RETURN
      END SUBROUTINE ad_main2d
#endif
