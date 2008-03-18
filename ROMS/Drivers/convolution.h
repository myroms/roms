      MODULE ocean_control_mod
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  ROMS/TOMS Convolution Driver:                                       !
!                                                                      !
!  This driver executes ROMS/TOMS standard nonlinear model.  It        !
!  controls the initialization, time-stepping, and finalization        !
!  of the nonlinear model execution following ESMF conventions:        !
!                                                                      !
!     initialize                                                       !
!     run                                                              !
!     finalize                                                         !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: ROMS_initialize
      PUBLIC  :: ROMS_run
      PUBLIC  :: ROMS_finalize

      CONTAINS

      SUBROUTINE ROMS_initialize (first, MyCOMM)
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

#ifdef DISTRIBUTE
!
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
      END SUBROUTINE ROMS_initialize

      SUBROUTINE ROMS_run (Tstr, Tend)
!
!=======================================================================
!                                                                      !
!  This routine time-steps ROMS/TOMS correlation diffusion operator.   !
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
!!    USE ad_balance_mod, ONLY: ad_balance
#endif
      USE ad_convolution_mod, ONLY : ad_convolution
      USE ad_variability_mod, ONLY : ad_variability
#ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti
#endif
      USE impulse_mod, ONLY : impulse
      USE ini_adjust_mod, ONLY : load_ADtoTL
      USE ini_adjust_mod, ONLY : load_TLtoAD
      USE normalization_mod, ONLY : normalization
#ifdef BALANCE_OPERATOR
!!    USE tl_balance_mod, ONLY: tl_balance
#endif
      USE tl_convolution_mod, ONLY : tl_convolution
      USE tl_variability_mod, ONLY : tl_variability
!
!  Imported variable declarations
!
      integer, dimension(Ngrids) :: Tstr
      integer, dimension(Ngrids) :: Tend
!
!  Local variable declarations.
!
      logical :: add, Lweak, outer_impulse

      integer :: IADrec, Nrec, i, ng, nvd, subs, tile, thread

      integer, dimension(4) :: Vsize

      character (len=80) :: fname
!
!-----------------------------------------------------------------------
!  Run model for all nested grids, if any.
!-----------------------------------------------------------------------
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
!  Compute or read in background-error correlations normalization
!  factors.
!-----------------------------------------------------------------------
!
!  If computing, write out factors to NetCDF. This is an expensive
!  computation and needs to be computed once for a particular
!  application grid.
!
        IF (LwrtNRM(ng)) THEN
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
!-----------------------------------------------------------------------
!  Read in adjoint solution (IADname), apply spatial convolution,
!  and then write output NetCDF file (ADJname).
!-----------------------------------------------------------------------
!
!  Inquire about the number of records in input NetCDF.
!
        IF (InpThread) THEN
          CALL opencdf (ng, 1, IADname(ng), fname, N(ng), 0, Nrec, nvd, &
     &                  Vsize)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Create convoluted adjoint solution NetCDF file.
!
        LdefADJ(ng)=.TRUE.
        CALL ad_def_his (ng, LdefADJ(ng))
        LdefADJ(ng)=.FALSE.
        LwrtADJ(ng)=.TRUE.
!
!  Process all available adjoint records.
!
        Lnew(ng)=1
        tADJindx(ng)=0
        NrecADJ(ng)=0
        LwrtState2d(ng)=.TRUE.
!
!  Proccess each time record of current adjoint solution in ADJname.
!
        DO i=1,Nrec
!
!  Set switch to scale model error covariace with background error
!  covariance factor Cfscale(:).
!
          IF (i.eq.Nrec) THEN
            Lweak=.FALSE.
          ELSE
            Lweak=.TRUE.
          END IF
!
!  Read adjoint solution. Since routine "get_state" loads data into the
!  ghost points, the adjoint solution is read in the tangent linear
!  state arrays by using iTLM instead of iADM in the calling arguments.
!
          IADrec=i
          CALL get_state (ng, iTLM, 4, IADname(ng), IADrec, Lnew(ng))
!
!  Load interior solution, read above, into adjoint state arrays. 
!  Then, multiply adjoint solution by the background-error standard
!  deviations. Next, convolve resulting adjoint solution with the 
!  squared-root adjoint diffusion operator which impose the model-error
!  spatial correlations. Notice that the spatial convolution is only
!  done for half of the diffusion steps (squared-root filter). Clear
!  tangent linear state arrays when done.
!
          add=.FALSE.
!$OMP PARALLEL DO PRIVATE(ng,thread,subs,tile)                          &
!$OMP&            SHARED(inner,add,numthreads)
          DO thread=0,numthreads-1
            subs=NtileX(ng)*NtileE(ng)/numthreads
            DO tile=subs*thread,subs*(thread+1)-1
              CALL load_TLtoAD (ng, TILE, Lnew(ng), Lnew(ng), add)
#ifdef BALANCE_OPERATOR
              CALL ad_balance (ng, TILE, Lnew(ng))
#endif
              CALL ad_variability (ng, TILE, Lnew(ng), Lweak)
              CALL ad_convolution (ng, TILE, Lnew(ng), 2)
              CALL initialize_ocean (ng, TILE, iTLM)
            END DO
          END DO
!$OMP END PARALLEL DO
!
!  To insure symmetry, convolve resulting filtered adjoint solution
!  from above with the squared-root (half of steps) tangent linear
!  diffusion operator. Then, multiply result with its corresponding
!  background-error standard deviations.  Since the convolved solution
!  is in the adjoint state arrays, first copy to tangent linear state
!  arrays including the ghosts points. Copy back to adjoint state
!  arrays when done with the convolution for output purposes.
!
          add=.FALSE.
!$OMP PARALLEL DO PRIVATE(ng,thread,subs,tile)                          &
!$OMP&            SHARED(inner,add,numthreads)
          DO thread=0,numthreads-1
            subs=NtileX(ng)*NtileE(ng)/numthreads
            DO tile=subs*thread,subs*(thread+1)-1,+1
              CALL load_ADtoTL (ng, TILE, Lnew(ng), Lnew(ng), add)
              CALL tl_convolution (ng, TILE, Lnew(ng), 2)
              CALL tl_variability (ng, TILE, Lnew(ng), Lweak)
#ifdef BALANCE_OPERATOR
              CALL tl_balance (ng, TILE, Lnew(ng))
#endif
              CALL load_TLtoAD (ng, TILE, Lnew(ng), Lnew(ng), add)
            END DO
          END DO
!$OMP END PARALLEL DO
!
!  Write out convolved solution to adjoint history NetCDF file.
!
          kstp(ng)=Lnew(ng)
#ifdef SOLVE3D
          nstp(ng)=Lnew(ng)
#endif
          CALL ad_wrt_his (ng)
          IF (exit_flag.ne.NoError) RETURN

        END DO
        LwrtState2d(ng)=.FALSE.
!
!  Create TLM/RPM impulse NetCDF file. Convert convolved adjoint
!  solution to impulse forcing. Write out impulse forcing into TLFname
!  NetCDF file. To facilitate the forcing by the TLM and RPM, the
!  forcing is process and written in increasing time coordinates.
!
        LdefTLF(ng)=.TRUE.
        tTLFindx(ng)=0
        CALL def_impulse (ng)
        outer_impulse=.FALSE.
#ifdef DISTRIBUTE
        tile=MyRank
#else
        tile=-1
#endif
        CALL impulse (ng, tile, iADM, outer_impulse, ADJname(ng))

      END DO NEST_LOOP

      RETURN
      END SUBROUTINE ROMS_run

      SUBROUTINE ROMS_finalize
!
!=======================================================================
!                                                                      !
!  This routine terminates ROMS/TOMS nonlinear model execution.        !
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
          CALL wclock_off (ng, iNLM, 0)
        END DO
!$OMP END PARALLEL DO
      END DO
!
!  Close IO files.
!
      CALL close_io

      RETURN
      END SUBROUTINE ROMS_finalize

      END MODULE ocean_control_mod
