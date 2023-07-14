#include "cppdefs.h"
#if defined TANGENT && !defined SOLVE3D
      SUBROUTINE tl_main2d (RunInterval)
!
!svn $Id: tl_main2d.F 285 2023-05-29 15:50:59Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine is the main driver for ROMS/TOMS perturbation tangent  !
!  linear model (TLM) when configured as a 2D barotropic shallow       !
!  water ocean model. It advances advances forward the TLM for all     !
!  nested grids, if any, by the specified time interval (seconds),     !
!  RunInterval.                                                        !
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
      USE mod_iounits
# ifdef NESTING
      USE mod_nesting
# endif
      USE mod_scalars
      USE mod_stepping
!
      USE dateclock_mod,           ONLY : time_string
# ifdef TLM_CHECK
      USE dotproduct_mod,          ONLY : tl_dotproduct
# endif
# ifdef TIDE_GENERATING_FORCES
      USE equilibrium_tide_mod,    ONLY : equilibrium_tide
# endif
# if defined ATM_COUPLING_NOT_YET && defined MCT_LIB
      USE mct_coupler_mod,         ONLY : ocn2atm_coupling
# endif
# if defined WAV_COUPLING_NOT_YET && defined MCT_LIB
      USE mct_coupler_mod,         ONLY : ocn2wav_coupling
# endif
      USE strings_mod,             ONLY : FoundError
      USE tl_diag_mod,             ONLY : tl_diag
# if defined WEAK_CONSTRAINT || defined FORCING_SV
      USE tl_forcing_mod,          ONLY : tl_forcing
# endif
# ifdef ADJUST_WSTRESS
      USE tl_frc_adjust_mod,       ONLY : tl_frc_adjust
# endif
      USE tl_ini_fields_mod,       ONLY : tl_ini_fields, tl_ini_zeta
# ifdef ADJUST_BOUNDARY
      USE tl_obc_adjust_mod,       ONLY : tl_obc_adjust
# endif
# ifdef WEC_MELLOR_NOT_YET
!!    USE tl_radiation_stress_mod, ONLY : tl_radiation_stress
# endif
# ifdef TL_AVERAGES
      USE tl_set_avg_mod,          ONLY : tl_set_avg
# endif
# if defined SSH_TIDES_NOT_YET || defined UV_TIDES_NOT_YET
!!    USE tl_set_tides_mod,        ONLY : tl_set_tides
# endif
      USE tl_set_vbc_mod,          ONLY : tl_set_vbc
      USE tl_step2d_mod,           ONLY : tl_step2d
# ifdef FLOATS_NOT_YET
!!    USE tl_step_floats_mod,      ONLY : tl_step_floats
# endif
!
      implicit none
!
!  Imported variable declarations.
!
      real(dp), intent(in) :: RunInterval
!
!  Local variable declarations.
!
      logical :: DoNestLayer, Time_Step
!
      integer :: Nsteps, Rsteps
      integer :: ig, il, istep, ng, nl, tile
      integer :: next_indx1
# ifdef FLOATS_NOT_YET
      integer :: Lend, Lstr, chunk_size
# endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
!=======================================================================
!  Time-step tangent linear vertically integrated equations.
!=======================================================================
!
!  Time-step the 3D kernel for the specified time interval (seconds),
!  RunInterval.
!
      Time_Step=.TRUE.
      DoNestLayer=.TRUE.
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
        TwoWayInterval(1:Ngrids)=0.0_r8
# endif
!
        NEST_LAYER : DO WHILE (DoNestLayer)
!
!  Determine number of time steps to compute in each nested grid layer
!  based on the specified time interval (seconds), RunInterval. Non
!  nesting applications have NestLayers=1. Notice that RunInterval is
!  set in the calling driver. Its value may span the full period of the
!  simulation, a multi-model coupling interval (RunInterval > ifac*dt),
!  or just a single step (RunInterval=0).
!
          CALL ntimesteps (iTLM, RunInterval, nl, Nsteps, Rsteps)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          IF ((nl.le.0).or.(nl.gt.NestLayers)) EXIT
!
!  Time-step governing equations for Nsteps.
!
          STEP_LOOP : DO istep=1,Nsteps
!
!  Set time indices and time clock.
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              tdays(ng)=time(ng)*sec2day
              IF (step_counter(ng).eq.Rsteps) Time_Step=.FALSE.
            END DO
!
!-----------------------------------------------------------------------
!  Read in required data, if any, from input NetCDF files.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              CALL tl_get_data (ng)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END DO
!
!-----------------------------------------------------------------------
!  If applicable, process input data: time interpolate between data
!  snapshots.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              DO tile=first_tile(ng),last_tile(ng),+1
                CALL tl_set_data (ng, tile)
              END DO
            END DO
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# if (defined WEAK_CONSTRAINT || defined FORCING_SV) && \
     !defined SP4DVAR
!
!-----------------------------------------------------------------------
!  If appropriate, add convolved adjoint solution impulse forcing to
!  the tangent linear model solution. Notice that the forcing is only
!  needed after finishing all inner loops. The forcing is continuous.
!  That is, it is time interpolated at every time-step from available
!  snapshots (FrequentImpulse=TRUE).
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
#  ifdef WEAK_CONSTRAINT
#   ifdef WEAK_NOINTERP
              IF ((iic(ng).gt.1).and.(iic(ng).ne.ntend(ng)+1).and.      &
     &            (MOD(iic(ng)-1,nadj(ng)).eq.0)) THEN
                IF (Master) THEN
                  WRITE (stdout,*) ' FORCING TLM at iic = ',iic(ng)
                END IF
#   endif
                IF (FrequentImpulse(ng)) THEN
                  DO tile=first_tile(ng),last_tile(ng),+1
                    CALL tl_forcing (ng, tile, kstp(ng), nstp(ng))
                  END DO
                END IF
#   ifdef WEAK_NOINTERP
              END IF
#   endif
#  else
              DO tile=first_tile(ng),last_tile(ng),+1
                CALL tl_forcing (ng, tile, kstp(ng), nstp(ng))
              END DO
#  endif
            END DO
# endif
!
!-----------------------------------------------------------------------
!  Initialize all time levels and compute other initial fields.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              IF (iic(ng).eq.ntstart(ng)) THEN
!
!  Initialize free-surface.
!
                DO tile=first_tile(ng),last_tile(ng),+1
                  CALL tl_ini_zeta (ng, tile, iTLM)
                END DO
!
!  Initialize other state variables.
!
                DO tile=last_tile(ng),first_tile(ng),-1
                  CALL tl_ini_fields (ng, tile, iTLM)
                END DO

# ifdef NESTING
!
!  Extract donor grid initial data at contact points and store it in
!  REFINED structure so it can be used for the space-time interpolation.
!
                IF (RefinedGrid(ng)) THEN
                  CALL tl_nesting (ng, iTLM, ngetD)
                END IF
# endif
              END IF
            END DO
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
# ifdef TL_AVERAGES
                CALL tl_set_avg (ng, tile)
# endif
# ifdef DIAGNOSTICS
!!              CALL tl_set_diags (ng, tile)
# endif
# ifdef TIDE_GENERATING_FORCES
                CALL equilibrium_tide (ng, tile, iTLM)
# endif
                CALL tl_diag (ng, tile)
# ifdef TLM_CHECK
                CALL tl_dotproduct (ng, tile, Lnew(ng))
# endif
              END DO
            END DO
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# if defined ATM_COUPLING_NOT_YET && defined MCT_LIB
!
!-----------------------------------------------------------------------
!  Couple ocean to atmosphere model every "CoupleSteps(Iatmos)"
!  timesteps: get air/sea fluxes.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              IF ((iic(ng).ne.ntstart(ng)).and.                         &
     &            MOD(iic(ng),CoupleSteps(Iatmos,ng)).eq.0) THEN
                DO tile=last_tile(ng),first_tile(ng),-1
                  CALL ocn2atm_coupling (ng, tile)
                END DO
              END IF
            END DO
# endif

# ifdef ADJUST_BOUNDARY
!
!-----------------------------------------------------------------------
!  Interpolate open boundary increments and adjust open boundaries.
!  Skip the last output timestep.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              IF (iic(ng).lt.(ntend(ng)+1)) THEN
                DO tile=first_tile(ng),last_tile(ng),+1
                  CALL tl_obc_adjust (ng, tile, Lbinp(ng))
                END DO
              END IF
            END DO
# endif

# ifdef ADJUST_WSTRESS
!
!-----------------------------------------------------------------------
!  Interpolate surface forcing increments and adjust surface forcing.
!  Skip the last output timestep.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              IF (iic(ng).lt.(ntend(ng)+1)) THEN
                DO tile=first_tile(ng),last_tile(ng),+1
                  CALL tl_frc_adjust (ng, tile, Lfinp(ng))
                END DO
              END IF
            END DO
# endif

# if defined WAV_COUPLING_NOT_YET && defined MCT_LIB
!
!-----------------------------------------------------------------------
!  Couple ocean to waves model every "CoupleSteps(Iwaves)"
!  timesteps: get waves/sea fluxes.
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

# ifdef WEC_MELLOR_NOT_YET
!
!-----------------------------------------------------------------------
!  Compute radiation stress terms.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              DO tile=last_tile(ng),first_tile(ng),-1
                CALL tl_radiation_stress (ng, tile)
              END DO
            END DO
# endif
!
!-----------------------------------------------------------------------
!  Set vertical boundary conditions. Process tidal forcing.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              DO tile=first_tile(ng),last_tile(ng),+1
                CALL tl_set_vbc (ng, tile)
# if defined SSH_TIDES_NOT_YET || defined UV_TIDES_NOT_YET
                CALL tl_set_tides (ng, tile)
# endif
              END DO
            END DO

# ifdef NESTING
!
!  If composite or mosaic grids, process additional points in the
!  contact zone between connected grids for bottom stress variables.
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              IF (ANY(CompositeGrid(:,ng))) THEN
                CALL tl_nesting (ng, iTLM, nbstr)
              END IF
            END DO
# endif
!
!-----------------------------------------------------------------------
!  If appropriate, write out fields into output NetCDF files.  Notice
!  that IO data is written in delayed and serial mode.  Exit if last
!  time step.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              CALL tl_output (ng)
              IF ((FoundError(exit_flag, NoError, __LINE__, MyFile)).or.&
     &            ((iic(ng).eq.(ntend(ng)+1)).and.(ng.eq.Ngrids))) THEN
                RETURN
              END IF
            END DO

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
                CALL tl_nesting (ng, iTLM, nputD)
              END IF
            END DO
# endif

# ifdef STEP2D_FB_AB3_AM4
!
!-----------------------------------------------------------------------
!  Solve tangent linear vertically integrated primitive equations for
!  free-surface and barotropic momentum components using a generalized
!  Forward-Backward, 3rd-order Adams-Bashforth / 4th-order Adams-Moulton
!  (FB AB3-AM4) time stepping scheme (Shchepetkin and McWilliams,
!  2009).
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              iif(ng)=1
              nfast(ng)=1
              kstp(ng)=knew(ng)
              knew(ng)=kstp(ng)+1
              IF (knew(ng).gt.4) knew(ng)=1
              IF (MOD(knew(ng),2).eq.0) THEN                ! zig-zag
                DO tile=first_tile(ng),last_tile(ng),+1     ! processing
                  CALL tl_step2d (ng, tile)                 ! sequence
                END DO
              ELSE
                DO tile=last_tile(ng),first_tile(ng),-1
                  CALL tl_step2d (ng, tile)
                END DO
              END IF
            END DO

# else

#  ifdef STEP2D_FB_LF_AM3
!
!-----------------------------------------------------------------------
!  Solve tangent linear vertically integrated primitive equations for
!  free-surface and momentum components using a predictor-corrector
!  LeapFrog / 3rd-order Adams-Moulton with a Forward-Backward
!  feeback (FB LF-AM3) stepping scheme (Shchepetkin and McWilliams,
!  2009).
!-----------------------------------------------------------------------
!
!  Predictor LF substep with FB-feedback.
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              iif(ng)=1
              nfast(ng)=1
              kstp(ng)=next_kstp(ng)
              knew(ng)=3

              DO tile=last_tile(ng),first_tile(ng),-1
                CALL tl_step2d (ng, tile)
              END DO
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
                CALL tl_nesting (ng, iTLM, n2dPS)
              END IF
            END DO
#   endif
!
!  Corrector AM3 substep with FB-feedback.
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              knew(ng)=3-kstp(ng)
              next_kstp(ng)=knew(ng)

              DO tile=first_tile(ng),last_tile(ng),+1
                CALL tl_step2d (ng, tile)
              END DO
            END DO

#  else
!
!-----------------------------------------------------------------------
!  Solve tangent linear vertically integrated primitive equations for
!  free-surface and momentum components using a predictor-corrector
!  LeapFrog with 3rd-order Adams-Moulton (LF-AM3) time stepping scheme
!  (ROMS legacy 2D kernel).
!-----------------------------------------------------------------------
!
!  Set time indices for predictor step. The PREDICTOR_2D_STEP switch
!  it is assumed to be false before the first time-step.
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
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
              DO tile=last_tile(ng),first_tile(ng),-1
                CALL tl_step2d (ng, tile)
              END DO
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
                CALL tl_nesting (ng, iTLM, n2dPS)
              END IF
            END DO
#   endif
!
!  Set time indices for corrector step.
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
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
                  CALL tl_step2d (ng, tile)
                END DO
              END IF
            END DO
#  endif

# endif

# ifdef NESTING
#  if defined MASKING && defined WET_DRY
!
!-----------------------------------------------------------------------
!  If nesting and wetting and drying, scale horizontal interpolation
!  weights to account for land/sea masking in contact areas. This needs
!  to be done at very time-step since the Land/Sea masking is time
!  dependent.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              CALL tl_nesting (ng, iTLM, nmask)
            END DO
#  endif
!
!-----------------------------------------------------------------------
!  If composite or mosaic grids, process additional points in the
!  contact zone between connected grids for the state variables
!  associated with the 2D engine CORRECTOR STEP section (KNEW INDEX).
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              IF (ANY(CompositeGrid(:,ng))) THEN
                CALL tl_nesting (ng, iTLM, n2dCS)
              END IF
            END DO
# endif

# ifdef NESTING
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
            DO il=NestLayers,1,-1
              DO ig=1,GridsInLayer(il)
                ng=GridNumber(ig,il)
                IF (do_twoway(iNLM, nl, il, ng, istep)) THEN
                  CALL tl_nesting (ng, iTLM, n2way)
                END IF
              END DO
            END DO
#  endif
!
!-----------------------------------------------------------------------
!  If donor to a finer grid, extract data for the external contact
!  points. This is the latest solution for the coarser grid.
!
!  It is stored in the REFINED structure so it can be used for the
!  space-time interpolation.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              IF (DonorToFiner(ng)) THEN
                CALL tl_nesting (ng, iTLM, ngetD)
              END IF
            END DO
# endif

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
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              IF (Lfloats(ng)) THEN
#  ifdef _OPENMP
                chunk_size=(Nfloats(ng)+numthreads-1)/numthreads
                Lstr=1+MyThread*chunk_size
                Lend=MIN(Nfloats(ng),Lstr+chunk_size-1)
#  else
                Lstr=1
                Lend=Nfloats(ng)
#  endif
                CALL tl_step_floats (ng, Lstr, Lend)
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
!
!-----------------------------------------------------------------------
!  Advance time index and time clock.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              iic(ng)=iic(ng)+1
              time(ng)=time(ng)+dt(ng)
              step_counter(ng)=step_counter(ng)-1
              CALL time_string (time(ng), time_code(ng))
            END DO

          END DO STEP_LOOP

        END DO NEST_LAYER

      END DO KERNEL_LOOP
!
      RETURN
      END SUBROUTINE tl_main2d
#else
      SUBROUTINE tl_main2d
      RETURN
      END SUBROUTINE tl_main2d
#endif
