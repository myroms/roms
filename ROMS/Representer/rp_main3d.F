#include "cppdefs.h"
#if defined TL_IOMS && defined SOLVE3D
      SUBROUTINE rp_main3d (RunInterval)
!
!svn $Id: rp_main3d.F 285 2023-05-29 15:50:59Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine is the main driver for ROMS/TOMS finite amplitude      !
!  tangent linear (Representers) model (RPM) when configure as a       !
!  full 3D baroclinic ocean model. It advances forward the RPM for     !
!  all nested grids, if any, for the specified time interval           !
!  (seconds), RunInterval.                                             !
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
# ifdef ANA_VMIX
      USE analytical_mod,          ONLY : ana_vmix
# endif
      USE dateclock_mod,           ONLY : time_string
# ifdef TIDE_GENERATING_FORCES
      USE equilibrium_tide_mod,    ONLY : equilibrium_tide
# endif
# if defined ATM_COUPLING_NOT_YET && defined MCT_LIB
      USE mct_coupler_mod,         ONLY : ocn2atm_coupling
# endif
# if defined WAV_COUPLING_NOT_YET && defined MCT_LIB
      USE mct_coupler_mod,         ONLY : ocn2wav_coupling
# endif
# ifdef FORWARD_READ
      USE omega_mod,               ONLY : omega
      USE set_depth_mod,           ONLY : set_depth
      USE set_massflux_mod,        ONLY : set_massflux
# endif
      USE strings_mod,             ONLY : FoundError
# ifdef BIOLOGY
      USE rp_biology_mod,          ONLY : rp_biology
# endif
# ifdef BBL_MODEL_NOT_YET
!!    USE rp_bbl_mod,              ONLY : rp_bblm
# endif
# if defined BULK_FLUXES_NOT_YET && !defined PRIOR_BULK_FLUXES
      USE rp_bulk_flux_mod,        ONLY : rp_bulk_flux
# endif
# ifdef BVF_MIXING_NOT_YET
!!    USE rp_bvf_mix_mod,          ONLY : rp_bvf_mix
# endif
      USE rp_diag_mod,             ONLY : rp_diag
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      USE rp_frc_adjust_mod,       ONLY : rp_frc_adjust
# endif
# ifdef GLS_MIXING_NOT_YET
!!    USE rp_gls_corstep_mod,      ONLY : rp_gls_corstep
!!    USE rp_gls_prestep_mod,      ONLY : rp_gls_prestep
# endif
      USE rp_ini_fields_mod,       ONLY : rp_ini_fields, rp_ini_zeta
# ifdef LMD_MIXING_NOT_YET
!!    USE rp_lmd_vmix_mod,         ONLY : rp_lmd_vmix
# endif
# ifdef MY25_MIXING
!!    USE rp_my25_corstep_mod,     ONLY : rp_my25_corstep
!!    USE rp_my25_prestep_mod,     ONLY : rp_my25_prestep
# endif
# ifdef ADJUST_BOUNDARY
      USE rp_obc_adjust_mod,       ONLY : rp_obc_adjust
      USE rp_obc_adjust_mod,       ONLY : rp_obc2d_adjust
      USE rp_set_depth_mod,        ONLY : rp_set_depth_bry
# endif
      USE rp_omega_mod,            ONLY : rp_omega
# ifdef WEC_MELLOR_NOT_YET
!!    USE rp_radiation_stress_mod, ONLY : rp_radiation_stress
# endif
# ifndef TS_FIXED
      USE rp_rho_eos_mod,          ONLY : rp_rho_eos
# endif
      USE rp_rhs3d_mod,            ONLY : rp_rhs3d
# ifdef SEDIMENT_NOT_YET
!!    USE rp_sediment_mod,         ONLY : rp_sediment
# endif
      USE rp_set_depth_mod,        ONLY : rp_set_depth
      USE rp_set_massflux_mod,     ONLY : rp_set_massflux
# if defined SSH_TIDES_NOT_YET || defined UV_TIDES_NOT_YET
!!    USE rp_set_tides_mod,        ONLY : rp_set_tides
# endif
      USE rp_set_vbc_mod,          ONLY : rp_set_vbc
# if !(defined STEP2D_FB_AB3_AM4 || defined STEP2D_FB_LF_AM3)
      USE rp_set_zeta_mod,         ONLY : rp_set_zeta
# endif
      USE rp_step2d_mod,           ONLY : rp_step2d
# ifndef TS_FIXED
      USE rp_step3d_t_mod,         ONLY : rp_step3d_t
# endif
      USE rp_step3d_uv_mod,        ONLY : rp_step3d_uv
# ifdef FLOATS_NOT_YET
!!    USE rp_step_floats_mod,      ONLY : rp_step_floats
# endif
# ifdef WEAK_CONSTRAINT
      USE tl_forcing_mod,          ONLY : tl_forcing
# endif
# ifdef RP_AVERAGES
      USE tl_set_avg_mod,          ONLY : tl_set_avg
# endif
!!    USE wvelocity_mod,           ONLY : wvelocity
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

      integer :: Nsteps, Rsteps
      integer :: ig, il, istep, ng, nl, tile
      integer :: my_iif, next_indx1
# ifdef FLOATS_NOT_YET
      integer :: Lend, Lstr, chunk_size
# endif
!
      real(r8) :: MaxDT, my_StepTime
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
!=======================================================================
!  Time-step representers tangent linear 3D primitive equations.
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
#ifdef NESTING
        TwoWayInterval(1:Ngrids)=0.0_r8
#endif
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
          CALL ntimesteps (iRPM, RunInterval, nl, Nsteps, Rsteps)
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
              nstp(ng)=1+MOD(iic(ng)-ntstart(ng),2)
              nnew(ng)=3-nstp(ng)
              nrhs(ng)=nstp(ng)
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
              CALL rp_get_data (ng)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END DO
!
!-----------------------------------------------------------------------
!  If applicable, process input data: time interpolate between data
!  snapshots. Compute BASIC STATE depths and thickness.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              DO tile=first_tile(ng),last_tile(ng),+1
                CALL rp_set_data (ng, tile)
# ifdef FORWARD_READ
            CALL set_depth (ng, tile, iRPM)
# endif
              END DO
            END DO
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef FORWARD_READ
!
!-----------------------------------------------------------------------
!  Compute BASIC STATE horizontal mass fluxes (Hz*u/n and Hz*v/m).
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              DO tile=last_tile(ng),first_tile(ng),-1
                CALL set_massflux (ng, tile, iRPM)
              END DO
            END DO
# endif

# ifdef WEAK_CONSTRAINT
!
!-----------------------------------------------------------------------
!  If appropriate, add convolved adjoint solution impulse forcing to
!  the representer model solution. Notice that the forcing is only
!  needed after finishing all inner loops. The forcing is continuous.
!  That is, it is time interpolated at every time-step from available
!  snapshots (FrequentImpulse=TRUE).
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              IF (FrequentImpulse(ng)) THEN
                DO tile=first_tile(ng),last_tile(ng),+1
                  CALL tl_forcing (ng, tile, kstp(ng), nstp(ng))
                END DO
              END IF
            END DO
# endif
!
!-----------------------------------------------------------------------
!  If not a restart, initialize all time levels and compute other
!  initial fields.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              IF (iic(ng).eq.ntstart(ng)) THEN
!
!  Initialize free-surface and compute initial level thicknesses and
!  depths.
!
                DO tile=first_tile(ng),last_tile(ng),+1
                  CALL rp_ini_zeta (ng, tile, iRPM)
                  CALL rp_set_depth (ng, tile, iRPM)
                END DO
!
!  Initialize other state variables.
!
                DO tile=last_tile(ng),first_tile(ng),-1
                  CALL rp_ini_fields (ng, tile, iRPM)
                END DO

# ifdef NESTING
!
!  Extract donor grid initial data at contact points and store it in
!  REFINED structure so it can be used for the space-time interpolation.
!
                IF (RefinedGrid(ng)) THEN
                  CALL rp_nesting (ng, iRPM, ngetD)
                END IF
# endif
              END IF
            END DO
!
!-----------------------------------------------------------------------
!  Compute horizontal mass fluxes (Hz*u/n and Hz*v/m), density related
!  quatities and report global diagnostics. Compute BASIC STATE omega
!  vertical velocity.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              DO tile=first_tile(ng),last_tile(ng),+1
                CALL rp_set_massflux (ng, tile, iRPM)
# ifndef TS_FIXED
                CALL rp_rho_eos (ng, tile, iRPM)
# endif
# ifdef TIDE_GENERATING_FORCES
                CALL equilibrium_tide (ng, tile, iRPM)
# endif
                CALL rp_diag (ng, tile)
# ifdef FORWARD_READ
                CALL omega (ng, tile, iRPM)
# endif
              END DO
            END DO
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

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

# ifdef WEC_MELLOR_NOT_YET
!
!-----------------------------------------------------------------------
!  Compute radiation stress terms.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              DO tile=last_tile(ng),first_tile(ng),-1
                CALL rp_radiation_stress (ng, tile)
              END DO
            END DO
# endif
!
!-----------------------------------------------------------------------
!  Set fields for vertical boundary conditions. Process tidal forcing,
!  if any.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              DO tile=first_tile(ng),last_tile(ng),+1
# if defined BULK_FLUXES_NOT_YET && !defined PRIOR_BULK_FLUXES
                CALL rp_bulk_flux (ng, tile)
# endif
# ifdef BBL_MODEL_NOT_YET
                CALL rp_bblm (ng, tile)
# endif
                CALL rp_set_vbc (ng, tile)
# if defined SSH_TIDES_NOT_YET || defined UV_TIDES_NOT_YET
                CALL rp_set_tides (ng, tile)
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
                CALL rp_nesting (ng, iRPM, nbstr)
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
                  CALL rp_obc_adjust (ng, tile, Lbinp(ng))
                  CALL rp_set_depth_bry (ng, tile, iRPM)
                  CALL rp_obc2d_adjust (ng, tile, Lbinp(ng))
                END DO
              END IF
            END DO
# endif

# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
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
                  CALL rp_frc_adjust (ng, tile, Lfinp(ng))
                END DO
              END IF
            END DO
# endif
!
!-----------------------------------------------------------------------
!  Compute tangent linear vertical mixing coefficients for momentum and
!  tracers. Compute S-coordinate vertical velocity, diagnostically from
!  horizontal mass divergence.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              DO tile=last_tile(ng),first_tile(ng),-1
# if defined ANA_VMIX_NOT_YET
                CALL rp_ana_vmix (ng, tile)
# elif defined LMD_MIXING_NOT_YET
                CALL rp_lmd_vmix (ng, tile)
# elif defined BVF_MIXING_NOT_YET
                CALL rp_bvf_mix (ng, tile)
# endif
                CALL rp_omega (ng, tile, iRPM)
!!              CALL wvelocity (ng, tile, nstp(ng))
              END DO
            END DO

# if !(defined STEP2D_FB_AB3_AM4 || defined STEP2D_FB_LF_AM3) || \
       defined DIAGNOSTICS       || defined AVERAGES
!
!-----------------------------------------------------------------------
!  Set free-surface to it time-averaged value.  If applicable,
!  accumulate time-averaged output data which needs a irreversible
!  loop in shared-memory jobs.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              DO tile=first_tile(ng),last_tile(ng),+1     ! irreversible
#  if !(defined STEP2D_FB_AB3_AM4 || defined STEP2D_FB_LF_AM3)
                CALL rp_set_zeta (ng, tile)
#  endif
#  ifdef DIAGNOSTICS
!!              CALL rp_set_diags (ng, tile)
#  endif
#  ifdef RP_AVERAGES
                CALL tl_set_avg (ng, tile)
#  endif
              END DO
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
                CALL rp_nesting (ng, iRPM, nzeta)
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
              CALL rp_output (ng)
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
#  ifdef NESTING_DEBUG
!
!  Also, fill BRY_CONTACT(:,:)%Mflux to check for mass conservation
!  between coarse and fine grids.  This is only done for diagnostic
!  purposes. Also, debugging is possible with very verbose output
!  to fort.300 is allowed by activating uppercase(nesting_debug).
#  endif
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              IF (RefinedGrid(ng).and.(RefineScale(ng).gt.0)) THEN
                CALL rp_nesting (ng, iRPM, nputD)
#  ifdef NESTING_DEBUG
                CALL rp_nesting (ng, iRPM, nmflx)
#  endif
              END IF
            END DO
# endif
!
!-----------------------------------------------------------------------
!  Compute right-hand-side terms for 3D equations.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              DO tile=last_tile(ng),first_tile(ng),-1
                CALL rp_rhs3d (ng, tile)
#  ifdef MY25_MIXING_NOT_YET
                CALL rp_my25_prestep (ng, tile)
#  elif defined GLS_MIXING_NOT_YET
                CALL rp_gls_prestep (ng, tile)
#  endif
              END DO
            END DO

# ifdef NESTING
!
!  If composite or mosaic grids, process additional points in the
!  contact zone between connected grids for right-hand-side terms
!  (tracers).
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              IF (ANY(CompositeGrid(:,ng))) THEN
                CALL rp_nesting (ng, iRPM, nrhst)
              END IF
            END DO
# endif

# ifdef STEP2D_FB_AB3_AM4
!
!-----------------------------------------------------------------------
!  Solve representers vertically integrated primitive equations for
!  free-surface and barotropic momentum components using a generalized
!  Forward-Backward, 3rd-order Adams-Bashforth / 4th-order Adams-Moulton
!  (FB AB3-AM4) time stepping scheme (Shchepetkin and McWilliams,
!  2009).
!-----------------------------------------------------------------------
!
            LOOP_2D : DO my_iif=1,MAXVAL(nfast)

              DO ig=1,GridsInLayer(nl)
                ng=GridNumber(ig,nl)
                IF (my_iif.le.nfast(ng)) THEN
                  iif(ng)=my_iif
                  kstp(ng)=knew(ng)
                  knew(ng)=kstp(ng)+1
                  IF (knew(ng).gt.4) knew(ng)=1

                  IF (MOD(knew(ng),2).eq.0) THEN            ! zig-zag
                    DO tile=first_tile(ng),last_tile(ng),+1 ! processing
                      CALL rp_step2d (ng, tile)             ! sequence
                    END DO
                  ELSE
                    DO tile=last_tile(ng),first_tile(ng),-1
                      CALL rp_step2d (ng, tile)
                    END DO
                  END IF
                END IF
              END DO

#  ifdef NESTING
!
!  If composite or mosaic grids, process additional points in the
!  contact zone between connected grids for the state variables
!  associated with the 2D engine CORRECTOR STEP section (KNEW INDEX).
!
              DO ig=1,GridsInLayer(nl)
                ng=GridNumber(ig,nl)
                IF (ANY(CompositeGrid(:,ng))) THEN
                  CALL rp_nesting (ng, iRPM, n2dCS)
                END IF
              END DO
#  endif
            END DO LOOP_2D

# else

#  ifdef STEP2D_FB_LF_AM3
!
!-----------------------------------------------------------------------
!  Solve the vertically integrated primitive equations for the
!  free-surface and barotropic momentum components using a predictor-
!  corrector LeapFrog / 3rd-order Adams-Moulton with a Forward-Backward
!  feeback (FB LF-AM3) stepping scheme (Shchepetkin and McWilliams,
!  2009).
!-----------------------------------------------------------------------
!
            LOOP_2D : DO my_iif=1,MAXVAL(nfast)
!
!  Predictor LF substep with FB-feedback.
!
              DO ig=1,GridsInLayer(nl)
                ng=GridNumber(ig,nl)
                IF (my_iif.le.nfast(ng)) THEN
                  iif(ng)=my_iif
                  kstp(ng)=next_kstp(ng)
                  knew(ng)=3

                  DO tile=last_tile(ng),first_tile(ng),-1
                    CALL rp_step2d (ng, tile)
                  END DO
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
                  CALL rp_nesting (ng, iRPM, n2dPS)
                END IF
              END DO
#   endif
!
!  Corrector AM3 substep with FB-feedback.
!
              DO ig=1,GridsInLayer(nl)
                ng=GridNumber(ig,nl)
                IF (my_iif.le.nfast(ng)) THEN
                  knew(ng)=3-kstp(ng)
                  next_kstp(ng)=knew(ng)

                  DO tile=first_tile(ng),last_tile(ng),+1
                    CALL rp_step2d (ng, tile)
                  END DO
                END IF
              END DO

#   ifdef NESTING
!
!  If composite or mosaic grids, process additional points in the
!  contact zone between connected grids for the state variables
!  associated with the 2D engine CORRECTOR STEP section (KNEW INDEX).
!
              DO ig=1,GridsInLayer(nl)
                ng=GridNumber(ig,nl)
                IF (ANY(CompositeGrid(:,ng))) THEN
                  CALL rp_nesting (ng, iRPM, n2dCS)
                END IF
              END DO
#   endif
            END DO LOOP_2D

#  else
!
!-----------------------------------------------------------------------
!  Solve the vertically integrated primitive equations for the
!  free-surface and barotropic momentum components using a predictor-
!  corrector LeapFrog with 3rd-order Adams-Moulton (LF-AM3) time
!  stepping scheme.
!-----------------------------------------------------------------------
!
            LOOP_2D : DO my_iif=1,MAXVAL(nfast)+1
!
!  Set time indices for predictor step. The PREDICTOR_2D_STEP switch
!  it is assumed to be false before the first time-step.
!
              DO ig=1,GridsInLayer(nl)
                ng=GridNumber(ig,nl)
                next_indx1=3-indx1(ng)
                IF (.not.PREDICTOR_2D_STEP(ng).and.                     &
     &              my_iif.le.(nfast(ng)+1)) THEN
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
                    CALL rp_step2d (ng, tile)
                  END DO
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
                  CALL rp_nesting (ng, iRPM, n2dPS)
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
                    CALL rp_step2d (ng, tile)
                  END DO
                END IF
              END DO

#   ifdef NESTING
!
!  If composite or mosaic grids, process additional points in the
!  contact zone between connected grids for the state variables
!  associated with the 2D engine CORRECTOR STEP section (KNEW INDEX).
!
              DO ig=1,GridsInLayer(nl)
                ng=GridNumber(ig,nl)
                IF (ANY(CompositeGrid(:,ng))) THEN
                  CALL rp_nesting (ng, iRPM, n2dCS)
                END IF
              END DO
#   endif
            END DO LOOP_2D
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
              CALL rp_nesting (ng, iRPM, nmask)
            END DO
#  endif
!
!-----------------------------------------------------------------------
!  If composite or mosaic grids, process additional points in the
!  contact zone between connected grids for the time-averaged
!  momentum fluxes (DU_avg1, DV_avg1) and free-surface (Zt_avg).
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              IF (ANY(CompositeGrid(:,ng))) THEN
                CALL rp_nesting (ng, iRPM, n2dfx)
              END IF
            END DO
# endif

# if !(defined STEP2D_FB_AB3_AM4 || defined STEP2D_FB_LF_AM3)
!
!-----------------------------------------------------------------------
!  Recompute depths and thicknesses using the new time filtered
!  free-surface.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              DO tile=last_tile(ng),first_tile(ng),-1
                CALL rp_set_depth (ng, tile, iRPM)
              END DO
            END DO
# endif

# ifdef NESTING
!
!  If nesting, determine vertical indices and vertical interpolation
!  weights in the contact zone using new depth arrays.
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              CALL rp_nesting (ng, iRPM, nzwgt)
            END DO
# endif
!
!-----------------------------------------------------------------------
!  Time-step 3D momentum equations.
!-----------------------------------------------------------------------
!
!  Time-step 3D momentum equations and couple with vertically
!  integrated equations.
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              DO tile=last_tile(ng),first_tile(ng),-1
                CALL rp_step3d_uv (ng, tile)
              END DO
            END DO

# ifdef NESTING
!
!  If composite or mosaic grids, process additional points in the
!  contact zone between connected grids for 3D momentum (u,v),
!  adjusted 2D momentum (ubar,vbar), and fluxes (Huon, Hvom).
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              IF (ANY(CompositeGrid(:,ng))) THEN
                CALL rp_nesting (ng, iRPM, n3duv)
              END IF
            END DO
# endif
!
!-----------------------------------------------------------------------
!  Time-step vertical mixing turbulent equations and passive tracer
!  source and sink terms, if applicable.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              DO tile=first_tile(ng),last_tile(ng),+1
                CALL rp_omega (ng, tile, iRPM)
# ifdef MY25_MIXING_NOT_YET
                CALL rp_my25_corstep (ng, tile)
# elif defined GLS_MIXING_NOT_YET
                CALL rp_gls_corstep (ng, tile)
# endif
# ifdef BIOLOGY
                CALL rp_biology (ng, tile)
# endif
# ifdef SEDIMENT_NOT_YET
                CALL rp_sediment (ng, tile)
# endif
              END DO
            END DO

# ifndef TS_FIXED
!
!-----------------------------------------------------------------------
!  Time-step tracer equations.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              DO tile=last_tile(ng),first_tile(ng),-1
                CALL rp_step3d_t (ng, tile)
              END DO
            END DO

#  ifdef NESTING
!
!  If composite or mosaic grids, process additional points in the
!  contact zone between connected grids for Tracer Variables.
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              IF (ANY(CompositeGrid(:,ng))) THEN
                CALL rp_nesting (ng, iRPM, n3dTV)
              END IF
            END DO
#  endif
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
                IF (do_twoway(iRPM, nl, il, ng, istep)) THEN
                  CALL rp_nesting (ng, iRPM, n2way)
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
!  space-time interpolation when "nputD" argument is used above.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              IF (DonorToFiner(ng)) THEN
                CALL rp_nesting (ng, iRPM, ngetD)
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
                CALL rp_step_floats (ng, Lstr, Lend)
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
      END SUBROUTINE rp_main3d
#else
      SUBROUTINE rp_main3d
      RETURN
      END SUBROUTINE rp_main3d
#endif
