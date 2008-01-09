      MODULE ocean_control_mod
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  ROMS/TOMS 4DVAR Background-error Correlation Model                  !
!                                                                      !
!  This driver is used to build and test the 4DVAR background-error    !
!  correlation model using a generalized difussion operator:           !
!                                                                      !
!            B = S C S                                                 !
!                                                                      !
!            C = C^(1/2) C^(T/2)                                       !
!                                                                      !
!      C^(1/2) = G L^(1/2) W^(-1/2)                                    !
!                                                                      !
!      C^(T/2) = W^(T/2) L^(T/2) G                                     !
!                                                                      !
!  where                                                               !
!                                                                      !
!         B : background-error covariance                              !
!         C : background-error correlation                             !
!         G : normalization coefficient matrix                         !
!         L : self-adjoint diffusion filter                            !
!         S : background-error standard deviation                      !
!         W : Grid cell area or volume diagonal matrix                 !
!                                                                      !
!  The routines in this driver control the initialization,  time-      !
!  stepping, and finalization of  ROMS/TOMS  model following ESMF      !
!  conventions:                                                        !
!                                                                      !
!     initialize                                                       !
!     run                                                              !
!     finalize                                                         !
!                                                                      !
!  Reference                                                           !
!                                                                      !
!     Weaver, A. and P. Courtier, 2001: Correlation modelling on       !
!       the sphere using a generalized diffusion equation, Q. J.       !
!       R. Meteorol. Soc., 127, 1815-1845.                             !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: initialize, run, finalize

      CONTAINS

      SUBROUTINE initialize (first, MyCOMM)
!
!=======================================================================
!                                                                      !
!  This routine allocates and initializes ROMS/TOMS state variables    !
!  and internal and external parameters.                               !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_iounits
      USE mod_scalars
!
#ifdef AIR_OCEAN 
      USE ocean_coupler_mod, ONLY : initialize_atmos_coupling
#endif
#ifdef WAVES_OCEAN
      USE ocean_coupler_mod, ONLY : initialize_waves_coupling
#endif
#ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti
#endif
!
!  Imported variable declarations.
!
      logical, intent(inout) :: first

      integer, intent(in), optional :: MyCOMM
!
!  Local variable declarations.
!
      logical :: allocate_vars = .TRUE.

      integer :: STDrec, ng, thread

#ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Set distribute-memory (MPI) world communictor.
!-----------------------------------------------------------------------
!
      IF (PRESENT(MyCOMM)) THEN
        OCN_COMM_WORLD=MyCOMM
      ELSE
        OCN_COMM_WORLD=MPI_COMM_WORLD
      END IF
#endif
!
!-----------------------------------------------------------------------
!  On first pass, initialize model parameters a variables for all
!  nested/composed grids.  Notice that the logical switch "first"
!  is used to allow multiple calls to this routine during ensemble
!  configurations.
!-----------------------------------------------------------------------
!
      IF (first) THEN
        first=.FALSE.
!
!  Initialize parallel parameters.
!
        CALL initialize_parallel
!
!  Initialize wall clocks.
!
        IF (Master) THEN
          WRITE (stdout,10)
 10       FORMAT (' Process Information:',/)
        END IF
        DO ng=1,Ngrids
!$OMP PARALLEL DO PRIVATE(thread) SHARED(ng,numthreads)
          DO thread=0,numthreads-1
            CALL wclock_on (ng, iNLM, 0)
          END DO
!$OMP END PARALLEL DO
        END DO

#if defined AIR_OCEAN || defined WAVES_OCEAN
!
!  Initialize coupling streams between model(s).
!
        DO ng=1,Ngrids
# ifdef AIR_OCEAN
          CALL initialize_atmos_coupling (ng, MyRank)
# endif
# ifdef WAVES_OCEAN
          CALL initialize_waves_coupling (ng, MyRank)
# endif
        END DO
#endif
!
!  Read in model tunable parameters from standard input. Initialize
!  "mod_param", "mod_ncparam" and "mod_scalar" modules.
!
        CALL inp_par (iNLM)
        IF (exit_flag.ne.NoError) THEN
          IF (Master) THEN
            WRITE (stdout,'(/,a,i3,/)') Rerror(exit_flag), exit_flag
          END IF
          RETURN
        END IF
!
!  Allocate and initialize modules variables.
!
        CALL mod_arrays (allocate_vars)
!
!  Allocate and initialize observation arrays.
!
        CALL initialize_fourdvar
!
!  Read in background/model error standard deviation factors and
!  spatial convolution diffusion coefficients.
!  
        STDrec=1
        DO ng=1,Ngrids
          CALL get_state (ng, 6, 6, STDname(ng), STDrec, 1)
#ifdef DISTRIBUTE
          CALL mp_bcasti (ng, iNLM, exit_flag, 1)
#endif
          IF (exit_flag.ne.NoError) RETURN
        END DO

      END IF

      RETURN
      END SUBROUTINE initialize

      SUBROUTINE run
!
!=======================================================================
!                                                                      !
!  This routine computes background-error correlations.                !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
      USE mod_stepping
!
#ifdef BALANCE_OPERATOR
!!    USE ad_balance_mod, ONLY: ad_balance
#endif
      USE ad_convolution_mod, ONLY : ad_convolution
      USE ad_variability_mod, ONLY : ad_variability
      USE analytical_mod, ONLY : ana_perturb
#ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti
#endif
      USE ini_adjust_mod, ONLY : load_ADtoTL
      USE ini_adjust_mod, ONLY : load_TLtoAD
      USE normalization_mod, ONLY : normalization
#ifdef BALANCE_OPERATOR
!!    USE tl_balance_mod, ONLY: tl_balance
#endif
      USE tl_convolution_mod, ONLY : tl_convolution
      USE tl_variability_mod, ONLY : tl_variability
!
!  Local variable declarations.
!
      logical :: add
      integer :: i, ng, subs, tile, thread
!
!=======================================================================
!  Run model for all nested grids, if any.
!=======================================================================
!
      NEST_LOOP : DO ng=1,Ngrids
!
!-----------------------------------------------------------------------
!  Initialize metrics.
!-----------------------------------------------------------------------
!
        CALL initial (ng)
!
!-----------------------------------------------------------------------
!  Get background-error covariance normalization matrix.
!-----------------------------------------------------------------------
!
!  Compute or read in background-error covariance normalization factors.
!  If computing, write out factors to NetCDF. This is an expensive
!  computation and needs to be computed once for an application grid.
!  
        IF (LwrtNRM(ng)) THEN
          CALL def_norm (ng)
!$OMP PARALLEL DO PRIVATE(ng,thread,subs,tile) SHARED(numthreads)
          DO thread=0,numthreads-1
            subs=NtileX(ng)*NtileE(ng)/numthreads
            DO tile=subs*thread,subs*(thread+1)-1
              CALL normalization (ng, TILE, 2)
            END DO
          END DO
!$OMP END PARALLEL DO
          LdefNRM(ng)=.FALSE.
          LwrtNRM(ng)=.FALSE.
        ELSE
          tNRMindx(ng)=1
          CALL get_state (ng, 5, 5, NRMname(ng), tNRMindx(ng), 1)
#ifdef DISTRIBUTE
          CALL mp_bcasti (ng, iNLM, exit_flag, 1)
#endif
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!-----------------------------------------------------------------------
!  Test correlation model.
!-----------------------------------------------------------------------
!
!  Initialize adjoint model state with a delta function at specified
!  point. Use USER parameters from standard input to perturb solution
!  in routine "ana_perturb". Then, convolve solution with the adjoint
!  diffusion operator.
! 
        ADmodel=.TRUE.
        Lnew(ng)=1
!$OMP PARALLEL DO PRIVATE(ng,thread,subs,tile) SHARED(numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1,+1
#ifdef BALANCE_OPERATOR
            CALL ad_balance (ng, TILE, Lnew(ng))
#endif
            CALL ana_perturb (ng, TILE, iADM)
            CALL ad_convolution (ng, TILE, Lnew(ng), 2)
          END DO
        END DO
!$OMP END PARALLEL DO
        ADmodel=.FALSE.
!
!  Initialize tangent linear model with convolved adjoint solution.
!  Then, apply tangent linear convolution.
!
        add=.FALSE.
!$OMP PARALLEL DO PRIVATE(ng,add,thread,subs,tile) SHARED(numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1
            CALL load_ADtoTL (ng, TILE, Lnew(ng), Lnew(ng), add)
            CALL tl_convolution (ng, TILE, Lnew(ng), 2)
#ifdef BALANCE_OPERATOR
            CALL tl_balance (ng, TILE, Lnew(ng))
#endif
            CALL load_TLtoAD (ng, TILE, Lnew(ng), Lnew(ng), add)
          END DO
        END DO
!$OMP END PARALLEL DO
!
!  Write out background error correlation in adjoint history NetCDF
!  file.
!
        kstp(ng)=Lnew(ng)
#ifdef SOLVE3D
        nstp(ng)=Lnew(ng)
#endif
        LdefADJ(ng)=.TRUE.
        LwrtADJ(ng)=.TRUE.
        LwrtState2d(ng)=.TRUE.
        CALL ad_def_his (ng, LdefADJ(ng))
#if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
        Ladjusted(ng)=.TRUE.
#endif
        CALL ad_wrt_his (ng)    
#if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
        Ladjusted(ng)=.FALSE.
#endif

      END DO NEST_LOOP

      RETURN
      END SUBROUTINE run

      SUBROUTINE finalize
!
!=======================================================================
!                                                                      !
!  This routine terminates ROMS/TOMS nonlinear, tangent linear, and    !
!  adjoint models execution.                                           !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
!
!  Local variable declarations.
!
      integer :: ng, thread
!
!-----------------------------------------------------------------------
!  If blowing-up, save latest model state into RESTART NetCDF file.
!-----------------------------------------------------------------------
!
!  If cycling restart records, write solution into record 3.
!
      DO ng=1,Ngrids
        IF (LwrtRST(ng).and.(exit_flag.eq.1)) THEN
          IF (Master) WRITE (stdout,10)
 10       FORMAT (/,' Blowing-up: Saving latest model state into ',     &
     &              ' RESTART file',/)
          IF (LcycleRST(ng).and.(NrecRST(ng).ge.2)) THEN
            tRSTindx(ng)=2
            LcycleRST(ng)=.FALSE.
          END IF
          blowup=exit_flag
          exit_flag=NoError
          CALL wrt_rst (ng)
        END IF
      END DO
!
!-----------------------------------------------------------------------
!  Stop model and time profiling clocks.  Close output NetCDF files.
!-----------------------------------------------------------------------
!
!  Stop time clocks.
!
      IF (Master) THEN
        WRITE (stdout,20)
 20     FORMAT (/,'Elapsed CPU time (seconds):',/)
      END IF

      DO ng=1,Ngrids
!$OMP PARALLEL DO PRIVATE(ng,thread) SHARED(numthreads)
        DO thread=0,numthreads-1
          CALL wclock_off (ng, iNLM, 0)
        END DO
!$OMP END PARALLEL DO
      END DO
!
!  Close IO files.
!
      CALL close_io

      RETURN
      END SUBROUTINE finalize

      END MODULE ocean_control_mod
