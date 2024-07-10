      MODULE roms_kernel_mod
!
!git $Id$
!svn $Id: tl_roms.h 1210 2024-01-03 22:03:03Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2024 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  ROMS/TOMS Tangent Linear Model Driver:                              !
!                                                                      !
!  This driver executes ROMS/TOMS tangent linear model.  It  controls  !
!  the initialization, time-stepping, and finalization of the tangent  !
!  model execution following ESMF conventions:                         !
!                                                                      !
!     ROMS_initialize                                                  !
!     ROMS_run                                                         !
!     ROMS_finalize                                                    !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_arrays
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
!
      USE close_io_mod,      ONLY : close_inp, close_out
      USE inp_par_mod,       ONLY : inp_par
#ifdef MCT_LIB
# ifdef ATM_COUPLING
      USE mct_coupler_mod,   ONLY : initialize_ocn2atm_coupling
# endif
# ifdef WAV_COUPLING
      USE mct_coupler_mod,   ONLY : initialize_ocn2wav_coupling
# endif
#endif
      USE stdout_mod,        ONLY : Set_StdOutUnit, stdout_unit
      USE strings_mod,       ONLY : FoundError
      USE wrt_rst_mod,       ONLY : wrt_rst
!
      implicit none
!
      PUBLIC  :: ROMS_initialize
      PUBLIC  :: ROMS_run
      PUBLIC  :: ROMS_finalize
!
      CONTAINS
!
      SUBROUTINE ROMS_initialize (first, mpiCOMM)
!
!=======================================================================
!                                                                      !
!  This routine allocates and initializes ROMS/TOMS state variables    !
!  and internal and external parameters.                               !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      logical, intent(inout) :: first
!
      integer, intent(in), optional :: mpiCOMM
!
!  Local variable declarations.
!
      logical :: allocate_vars = .TRUE.

#ifdef DISTRIBUTE
      integer :: MyError, MySize
#endif
      integer :: chunk_size, ng, thread
#ifdef _OPENMP
      integer :: my_threadnum
#endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", ROMS_initialize"

#ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Set distribute-memory (mpi) world communictor.
!-----------------------------------------------------------------------
!
      IF (PRESENT(mpiCOMM)) THEN
        OCN_COMM_WORLD=mpiCOMM
      ELSE
        OCN_COMM_WORLD=MPI_COMM_WORLD
      END IF
      CALL mpi_comm_rank (OCN_COMM_WORLD, MyRank, MyError)
      CALL mpi_comm_size (OCN_COMM_WORLD, MySize, MyError)
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
!  Initialize parallel control switches. These scalars switches are
!  independent from standard input parameters.
!
        CALL initialize_parallel
!
!  Set the ROMS standard output unit to write verbose execution info.
!  Notice that the default standard out unit in Fortran is 6.
!
!  In some applications like coupling or disjointed mpi-communications,
!  it is advantageous to write standard output to a specific filename
!  instead of the default Fortran standard output unit 6. If that is
!  the case, it opens such formatted file for writing.
!
        IF (Set_StdOutUnit) THEN
          stdout=stdout_unit(Master)
          Set_StdOutUnit=.FALSE.
        END IF
!
!  Read in model tunable parameters from standard input. Allocate and
!  initialize variables in several modules after the number of nested
!  grids and dimension parameters are known.
!
        CALL inp_par (iTLM)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Set domain decomposition tile partition range.  This range is
!  computed only once since the "first_tile" and "last_tile" values
!  are private for each parallel thread/node.
!
!$OMP PARALLEL
#if defined _OPENMP
      MyThread=my_threadnum()
#elif defined DISTRIBUTE
      MyThread=MyRank
#else
      MyThread=0
#endif
      DO ng=1,Ngrids
        chunk_size=(NtileX(ng)*NtileE(ng)+numthreads-1)/numthreads
        first_tile(ng)=MyThread*chunk_size
        last_tile (ng)=first_tile(ng)+chunk_size-1
      END DO
!$OMP END PARALLEL
!
!  Initialize internal wall clocks. Notice that the timings does not
!  includes processing standard input because several parameters are
!  needed to allocate clock variables.
!
        IF (Master) THEN
          WRITE (stdout,10)
 10       FORMAT (/,' Process Information:',/)
        END IF
!
        DO ng=1,Ngrids
!$OMP PARALLEL
          DO thread=THREAD_RANGE
            CALL wclock_on (ng, iTLM, 0, __LINE__, MyFile)
          END DO
!$OMP END PARALLEL
        END DO
!
!  Allocate and initialize modules variables.
!
!$OMP PARALLEL
        CALL ROMS_allocate_arrays (allocate_vars)
        CALL ROMS_initialize_arrays
!$OMP END PARALLEL
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      END IF

#if defined MCT_LIB && (defined ATM_COUPLING || defined WAV_COUPLING)
!
!-----------------------------------------------------------------------
!  Initialize coupling streams between model(s).
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
# ifdef ATM_COUPLING
        CALL initialize_ocn2atm_coupling (ng, MyRank)
# endif
# ifdef WAV_COUPLING
        CALL initialize_ocn2wav_coupling (ng, MyRank)
# endif
      END DO
#endif
!
!-----------------------------------------------------------------------
!  Initialize tangent linear model state variables over all nested
!  grids, if applicable.
!-----------------------------------------------------------------------
!
#ifdef FORWARD_FLUXES
!
!  Set the BLK structure to contain the nonlinear model surface fluxes
!  needed by the tangent linear and adjoint models. Also, set switches
!  to process that structure in routine "check_multifile". Notice that
!  it is possible to split the solution into multiple NetCDF files to
!  reduce their size.
!
!  The switch LreadFRC is deactivated because all the atmospheric
!  forcing, including shortwave radiation, is read from the NLM
!  surface fluxes or is assigned during ESM coupling.  Such fluxes
!  are available from the QCK structure. There is no need for reading
!  and processing from the FRC structure input forcing-files.
!
      CALL edit_multifile ('QCK2BLK')
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      DO ng=1,Ngrids
        LreadBLK(ng)=.TRUE.
        LreadFRC(ng)=.FALSE.
      END DO
#endif
!
!  Initialize perturbation tangent linear model.
!
      Lstiffness=.FALSE.
      DO ng=1,Ngrids
        LreadFWD(ng)=.TRUE.
!$OMP PARALLEL
        CALL tl_initial (ng)
!$OMP END PARALLEL
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END DO
!
!  Initialize run or ensemble counter.
!
      Nrun=1
!
!  Activate tangent linear output.
!
      DO ng=1,Ngrids
        LdefTLM(ng)=.TRUE.
        LwrtTLM(ng)=.TRUE.
        LcycleTLM(ng)=.FALSE.
      END DO
!
      RETURN
      END SUBROUTINE ROMS_initialize
!
      SUBROUTINE ROMS_run (RunInterval)
!
!=======================================================================
!                                                                      !
!  This routine runs ROMS/TOMS tangent linear model for the specified  !
!  time interval (seconds), RunInterval.                               !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      real(dp), intent(in) :: RunInterval            ! seconds
!
!  Local variable declarations.
!
      integer :: ng
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", ROMS_run"
!
!-----------------------------------------------------------------------
!  Time-step tangent linear model over all nested grids, if applicable.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        IF (Master) THEN
          WRITE (stdout,10) 'TL', ng, ntstart(ng), ntend(ng)
        END IF
      END DO

!$OMP PARALLEL
#ifdef SOLVE3D
      CALL tl_main3d (RunInterval)
#else
      CALL tl_main2d (RunInterval)
#endif
!$OMP END PARALLEL

      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
 10   FORMAT (/,1x,a,1x,'ROMS/TOMS: started time-stepping:',            &
     &        ' (Grid: ',i2.2,' TimeSteps: ',i8.8,' - ',i8.8,')',/)
!
      RETURN
      END SUBROUTINE ROMS_run
!
      SUBROUTINE ROMS_finalize
!
!=======================================================================
!                                                                      !
!  This routine terminates ROMS/TOMS tangent linear model execution.   !
!                                                                      !
!=======================================================================
!
!  Local variable declarations.
!
      integer :: Fcount, ng, thread
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", ROMS_finalize"
!
!-----------------------------------------------------------------------
!  If blowing-up, save latest model state into RESTART NetCDF file.
!-----------------------------------------------------------------------
!
!  If cycling restart records, write solution into the next record.
!
      IF (exit_flag.eq.1) THEN
        DO ng=1,Ngrids
          IF (LwrtRST(ng)) THEN
            IF (Master) WRITE (stdout,10)
 10         FORMAT (/,' Blowing-up: Saving latest model state into ',   &
     &                ' RESTART file',/)
            IF (LcycleRST(ng).and.(RST(ng)%Nrec(Fcount).ge.2)) THEN
              RST(ng)%Rindex=2
              LcycleRST(ng)=.FALSE.
            END IF
            blowup=exit_flag
            exit_flag=NoError
#ifdef DISTRIBUTE
            CALL wrt_rst (ng, MyRank)
#else
            CALL wrt_rst (ng, -1)
#endif
          END IF
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Stop model and time profiling clocks, report memory requirements, and
!  close output NetCDF files.
!-----------------------------------------------------------------------
!
!  Stop time clocks.
!
      IF (Master) THEN
        WRITE (stdout,20)
 20     FORMAT (/,'Elapsed wall CPU time for each process (seconds):',/)
      END IF

      DO ng=1,Ngrids
!$OMP PARALLEL
        DO thread=THREAD_RANGE
          CALL wclock_off (ng, iTLM, 0, __LINE__, MyFile)
        END DO
!$OMP END PARALLEL
      END DO
!
!  Report dynamic memory and automatic memory requirements.
!
!$OMP PARALLEL
      CALL memory
!$OMP END PARALLEL
!
!  Close IO files.
!
      DO ng=1,Ngrids
        CALL close_inp (ng, iTLM)
      END DO
      CALL close_out
!
      RETURN
      END SUBROUTINE ROMS_finalize

      END MODULE roms_kernel_mod
