      MODULE ocean_control_mod
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2007 The ROMS/TOMS Group           W. G. Zhang   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  ROMS/TOMS Optimal Observation Driver:                               !
!                                                                      !
!  This driver computes the  adjoint sensitivity  of a function or     !
!  index, J,  in terms of space and/or time integrals of the model     !
!  state, S(zeta,u,v,T,...).  Small changes, dS, in S will lead to     !
!  changes dJ in J:                                                    !
!                                                                      !
!  dJ = (dJ/dzeta) dzeta + (dJ/du) du + (dJ/dv) dv + (dJ/dt) dT ...    !
!                                                                      !
!  and                                                                 !
!                                                                      !
!  dJ/dS = transpose(R) S                                              !
!                                                                      !
!  where  transpose(R) is the adjoint propagator.  It implies that     !
!  the sensitivity for ALL variables,  parameters,  and space-time     !
!  points can be computed from a single integration of the adjoint     !
!  model.                                                              !
!                                                                      !
!  These routines control the initialization,  time-stepping,  and     !
!  finalization of  ROMS/TOMS  model following ESMF conventions:       !
!                                                                      !
!     initialize                                                       !
!     run                                                              !
!     finalize                                                         !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
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
      USE atm_coupler_mod, ONLY : initialize_coupling
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
            CALL wclock_on (ng, iADM, 0)
          END DO
!$OMP END PARALLEL DO
        END DO

#ifdef AIR_OCEAN 
!
!  Initialize coupling streams between atmosphere and ocean using the
!  Model Coupling Toolkit (MCT).
!
        CALL initialize_coupling (MyRank)
#endif
!
!  Read in model tunable parameters from standard input. Initialize
!  "mod_param", "mod_ncparam" and "mod_scalar" modules.
!
        CALL inp_par (iADM)
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
!  Read in background-error standard deviation factors and spatial
!  convolution diffusion coefficients.
!
        STDrec=1
        DO ng=1,Ngrids
          CALL get_state (ng, 6, 6, STDname(ng), STDrec, 1)
#ifdef DISTRIBUTE
          CALL mp_bcasti (ng, iADM, exit_flag, 1)
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
!  This routine computes the adjoint sensitivity analysis, dJ/dS,      !
!  to the specified functional J.  The sensitivity masking arrays      !
!  Rscope, Uscope, and Vscope are used to evaluate the functional      !
!  in the desired spatial area.                                        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
      USE mod_ocean
      USE mod_scalars
      USE mod_stepping
!
#ifdef BALANCE_OPERATOR
      USE ad_balance_mod, ONLY: ad_balance
#endif
      USE ad_convolution_mod, ONLY : ad_convolution
      USE ad_variability_mod, ONLY : ad_variability
#ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti
#endif
      USE ini_adjust_mod, ONLY : load_ADtoTL
      USE ini_adjust_mod, ONLY : load_TLtoAD
      USE normalization_mod, ONLY : normalization
#ifdef BALANCE_OPERATOR
      USE tl_balance_mod, ONLY: tl_balance
#endif
      USE tl_convolution_mod, ONLY : tl_convolution
      USE tl_variability_mod, ONLY : tl_variability
!
!  Local variable declarations.
!
      logical :: add
      
      integer :: Nrec, i, my_iic, ng, nvd, thread, subs, tile

      integer, dimension(4) :: Vsize

      real (r8) :: str_day, end_day

      character (len=80) :: fname
!
!=======================================================================
!  Run model for all nested grids, if any.
!=======================================================================
!
      NEST_LOOP : DO ng=1,Ngrids
!
!-----------------------------------------------------------------------
!  Initialize adjoint model and define sensitivity functional.
!-----------------------------------------------------------------------
!
        Lstiffness=.FALSE.
        CALL ad_initial (ng)
        IF (exit_flag.ne.NoError) THEN
          IF (Master) THEN
            WRITE (stdout,10) Rerror(exit_flag), exit_flag
          END IF
          RETURN
        END IF
!
!  Activate adjoint output.
!
        LdefADJ(ng)=.TRUE.
        LwrtADJ(ng)=.TRUE.
        LcycleADJ(ng)=.FALSE.
!
!-----------------------------------------------------------------------
!  Compute or read in background-error correlations normalization
!  factors.
!-----------------------------------------------------------------------
!
!  If computing, write out factors to NetCDF. This is an expensive
!  computation and needs to be computed once for a particular
!  application grid.
!
        IF (LdefNRM(ng)) THEN
          CALL def_norm (ng)
!$OMP PARALLEL DO PRIVATE(ng,thread,subs,tile)                          &
!$OMP&            SHARED(inner,numthreads)
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
!------------------------------------------------------------------------
!  Time-step adjoint model: Compute the gradient or index, dJ/dS, of
!  the sensitivity functional.
!------------------------------------------------------------------------
!
        str_day=time(ng)*sec2day
        end_day=str_day-ntimes*dt(ng)*sec2day
        IF ((DstrS(ng).eq.0.0_r8).and.(DendS(ng).eq.0.0_r8)) THEN
          DstrS(ng)=end_day
          DendS(ng)=str_day
        END IF
        IF (Master) THEN
          WRITE (stdout,20) ntstart, ntend, DendS(ng), DstrS(ng)
        END IF
        IF ((DstrS(ng).gt.str_day).or.(DstrS(ng).lt.end_day)) THEN
          IF (Master)  WRITE (stdout,30) 'DstrS = ', DstrS(ng),         &
     &                                   end_day, str_day
          exit_flag=7
          RETURN
        END IF
        IF ((DendS(ng).gt.str_day).or.(DendS(ng).lt.end_day)) THEN
          IF (Master)  WRITE (stdout,30) 'DendS = ', DendS(ng),         &
     &                                   end_day, str_day
          exit_flag=7
          RETURN
        END IF

        time(ng)=time(ng)+dt(ng)

        AD_LOOP : DO my_iic=ntstart,ntend,-1

          iic(ng)=my_iic
#ifdef SOLVE3D
          CALL ad_main3d (ng)
#else
          CALL ad_main2d (ng)
#endif
          IF (Master.and.(exit_flag.ne.NoError)) THEN
            WRITE (stdout,10) Rerror(exit_flag), exit_flag
            RETURN
          END IF

        END DO AD_LOOP
!
!-----------------------------------------------------------------------
!  Read in adjoint solution from adjoint history file (ADJname), apply 
!  spatial convolution, and then write output NetCDF file (ADJname).
!-----------------------------------------------------------------------
!
        IF (Master) THEN
          WRITE (stdout,40)
        END IF
!
!  Inquire about the number of records in input NetCDF.
!
        IF (InpThread) THEN
          CALL opencdf (ng, 1, ADJname(ng), fname, N(ng), 0, Nrec, nvd, &
     &                  Vsize)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Process last available adjoint record.
!
        Lnew(ng)=1
        tADJindx(ng)=0
        NrecADJ(ng)=0
!
!  Read adjoint solution.
!
        CALL get_state (ng, iADM, 4, ADJname(ng), Nrec, Lnew(ng))
!
!  First, multiply adjoint solution by the background-error standard
!  deviations.  Second, convolve resulting adjoint solution with the
!  adjoint diffusion operator which embeds background-error spatial
!  correlations. Notice that the spatial convolution is only done
!  for half of the diffusion steps.
!
!$OMP PARALLEL DO PRIVATE(ng,thread,subs,tile)                          &
!$OMP&            SHARED(inner,CGstepF,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1
#ifdef BALANCE_OPERATOR
            CALL ad_balance (ng, TILE, Lnew(ng))
#endif
            CALL ad_variability (ng, TILE, Lnew(ng))
            CALL ad_convolution (ng, TILE, Lnew(ng), 2)
          END DO
        END DO
!$OMP END PARALLEL DO
!
!  To insure symmetry, convolve resulting filtered adjoint solution
!  from above with the tangent linear diffusion operator for the
!  other half of steps. Then, multiply result with its corresponding
!  background-error standard deviations.
!
        add=.FALSE.
!$OMP PARALLEL DO PRIVATE(ng,thread,subs,tile)                          &
!$OMP&            SHARED(inner,add,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1,+1
            CALL load_ADtoTL (ng, TILE, Lnew(ng), Lnew(ng), add)
            CALL tl_convolution (ng, TILE, Lnew(ng), 2)
            CALL tl_variability (ng, TILE, Lnew(ng))
#ifdef BALANCE_OPERATOR
            CALL tl_balance (ng, TILE, Lnew(ng))
#endif
          END DO
        END DO
!$OMP END PARALLEL DO
!
!  Write out convolved solution to tangent linear initial NetCDF file.
!
        LdefITL(ng)=.TRUE.
        kstp(ng)=Lnew(ng)
#ifdef SOLVE3D
        nstp(ng)=Lnew(ng)
#endif
        CALL tl_def_ini (ng)
        CALL tl_wrt_ini (ng, Lnew(ng), 1)
        IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Run tangent linear model for all nested grids, if any.
!-----------------------------------------------------------------------
!
!  Initialize tangent linear model.
!
        CALL tl_initial (ng)
        IF (exit_flag.ne.NoError) THEN
          IF (Master) THEN
            WRITE (stdout,10) Rerror(exit_flag), exit_flag
          END IF
          RETURN
        END IF
!
!  Activate tangent linear output.
!
        LdefTLM(ng)=.TRUE.
        LwrtTLM(ng)=.TRUE.
        LcycleTLM(ng)=.FALSE.
!
!  Time-step tangent linear model.
!
        IF (Master) THEN
          WRITE (stdout,50) ntstart, ntend
        END IF

        time(ng)=time(ng)-dt(ng)

        TL_LOOP : DO my_iic=ntstart,ntend+1

          iic(ng)=my_iic
#ifdef SOLVE3D
          CALL tl_main3d (ng)
#else
          CALL tl_main2d (ng)
#endif
          IF (exit_flag.ne.NoError) THEN
            IF (Master) THEN
	      WRITE (stdout,60) Rerror(exit_flag), exit_flag
            END IF
            RETURN
          END IF

        END DO TL_LOOP

      END DO NEST_LOOP
!
 10   FORMAT (/,a,i2,/)
 20   FORMAT (/,'AD ROMS/TOMS: started time-stepping:',                 &
     &        ' (TimeSteps: ',i8.8,' - ',i8.8,')',/,14x,                &
     &        'adjoint forcing time range: ',f12.4,' - ',f12.4 ,/)
 30   FORMAT (/,' Out of range adjoint forcing time, ',a,f12.4,/,       &
     &        ' It must be between ',f12.4,' and ',f12.4)
 40   FORMAT(/,'AD ROMS/TOMS: convolving final Adjoint solution:',/)
 50   FORMAT (/,'TL ROMS/TOMS: started time-stepping:',                 &
     &            '( TimeSteps: ',i8.8,' - ',i8.8,')',/)
 60   FORMAT (/,a,i3,/)

      RETURN
      END SUBROUTINE run

      SUBROUTINE finalize
!
!=======================================================================
!                                                                      !
!  This routine terminates ROMS/TOMS nonlinear and adjoint models      !
!  execution.                                                          !
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
!  If cycling restart records, write solution into the next record.
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
 20     FORMAT (/,' Elapsed CPU time (seconds):',/)
      END IF

      DO ng=1,Ngrids
!$OMP PARALLEL DO PRIVATE(thread) SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          CALL wclock_off (ng, iADM, 0)
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
