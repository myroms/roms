      PROGRAM mct_driver
!
!git $Id$
!=======================================================================
!  Copyright (c) 2002-2024 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                            Hernan G. Arango   !
!==================================================== John C. Warner ===
!                                                                      !
!  Master program to couple ROMS/TOMS to other models using the Model  !
!  Coupling Toolkit (MCT) library.                                     !
!                                                                      !
!  The following models are coupled to ROMS/TOMS:                      !
!                                                                      !
#ifdef WRF_COUPLING
!  WRF, Weather Research and Forecasting model:                        !
!       http://www.wrf-model.org                                       !
!                                                                      !
#endif
#ifdef SWAN_COUPLING
!  SWAN, Simulating WAves Nearshore model:                             !
!        http://vlm089.citg.tudelft.nl/swan/index.htm                  !
!                                                                      !
#endif
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_coupler
      USE mod_iounits
      USE mod_scalars
!
      USE m_MCTWorld,        ONLY : MCTWorld_clean => clean

#ifdef WRF_COUPLING
      USE mct_coupler_mod,   ONLY : finalize_ocn2atm_coupling
#endif
#if defined SWAN_COUPLING || defined REFDIF_COUPLING
      USE mct_coupler_mod,   ONLY : finalize_ocn2wav_coupling
#endif
      USE roms_kernel_mod,   ONLY : ROMS_initialize
      USE roms_kernel_mod,   ONLY : ROMS_run
      USE roms_kernel_mod,   ONLY : ROMS_finalize
#if defined PIO_LIB && defined DISTRIBUTE
      USE set_pio_mod,       ONLY : finalize_pio
#endif
!
      implicit none
!
!  Local variable declarations.
!
      logical, save :: first

      integer :: MyColor, MyCOMM, MyError, MyKey, Nnodes
      integer :: ng

      real(r4) :: CouplingTime             ! single precision
!
!-----------------------------------------------------------------------
!  Initialize distributed-memory (MPI) configuration
!-----------------------------------------------------------------------
!
!  Initialize MPI execution environment.
!
      CALL mpi_init (MyError)
!
!  Get rank of the local process in the group associated with the
!  comminicator.
!
      CALL mpi_comm_size (MPI_COMM_WORLD, Nnodes, MyError)
      CALL mpi_comm_rank (MPI_COMM_WORLD, MyRank, MyError)
!
!  Set temporarily the ocean communicator to current handle before
!  splitting so the input coupling script name can be broadcasted to
!  all the nodes.
!
      OCN_COMM_WORLD=MPI_COMM_WORLD
!
!  Read in coupled model parameters from standard input.
!
      CALL read_CouplePar (iNLM)
!
!  Allocate several coupling variables.
!
      CALL allocate_coupler (Nnodes)
!
!  Split the communicator into coupled models sub-groups based
!  on color and key.
!
      MyKey=0
      IF ((pets(Iocean)%val(1).le.MyRank).and.                          &
     &    (MyRank.le.pets(Iocean)%val(Nthreads(Iocean)))) THEN
        MyColor=OCNid
      END IF
#ifdef ATM_COUPLING
      IF ((pets(Iatmos)%val(1).le.MyRank).and.                          &
     &    (MyRank.le.pets(Iatmos)%val(Nthreads(Iatmos)))) THEN
        MyColor=ATMid
      END IF
#endif
#ifdef WAV_COUPLING
      IF ((pets(Iwaves)%val(1).le.MyRank).and.                          &
     &    (MyRank.le.pets(Iwaves)%val(Nthreads(Iwaves)))) THEN
        MyColor=WAVid
      END IF
#endif
      CALL mpi_comm_split (MPI_COMM_WORLD, MyColor, MyKey, MyCOMM,      &
     &                     MyError)
!
!-----------------------------------------------------------------------
!  Run coupled models according to the processor rank.
!-----------------------------------------------------------------------
!
#if defined SWAN_COUPLING
      IF (MyColor.eq.WAVid) THEN
        CouplingTime=REAL(TimeInterval(Iocean,Iwaves))
        CALL SWAN_INITIALIZE (MyCOMM, INPname(Iwaves))
        CALL SWAN_RUN (CouplingTime)
        CALL SWAN_FINALIZE
      END IF
#elif defined REFDIF_COUPLING
      IF (MyColor.eq.WAVid) THEN
        CouplingTime=REAL(TimeInterval(Iocean,Iwaves))
        CALL refdif_initialize (MyCOMM)
        CALL refdif_run (CouplingTime, INPname(Iwaves))
        CALL refdif_finalize
      END IF
#endif
#ifdef WRF_COUPLING
      IF (MyColor.eq.ATMid) THEN
        CouplingTime=REAL(TimeInterval(Iocean,Iatmos))
!!      CALL module_wrf_top_mp_wrf_init (MyCOMM)
!!      CALL module_wrf_top_mp_wrf_run (TimeInterval(Iocean,Iwaves))
!!      CALL module_wrf_top_mp_wrf_finalize
        CALL module_wrf_top_wrf_init (MyCOMM)
        CALL module_wrf_top_wrf_run (CouplingTime)
        CALL module_wrf_top_wrf_finalize
      END IF
#endif
      IF (MyColor.eq.OCNid) THEN
        first=.TRUE.
        IF (exit_flag.eq.NoError) THEN
          CALL ROMS_initialize (first, mpiCOMM=MyCOMM)
          run_time=0.0_r8
          DO ng=1,Ngrids
            run_time=MAX(run_time, dt(ng)*ntimes(ng))
          END DO
        END IF
        IF (exit_flag.eq.NoError) THEN
          CALL ROMS_run (run_time)
        END IF
        CALL ROMS_finalize
#if defined PIO_LIB && defined DISTRIBUTE
        CALL finalize_pio
#endif
#if defined SWAN_COUPLING || defined REFDIF_COUPLING
        CALL finalize_ocn2wav_coupling
#endif
#ifdef WRF_COUPLING
        CALL finalize_ocn2atm_coupling
#endif
      END IF
!
!-----------------------------------------------------------------------
!  Terminates all the mpi-processing and coupling.
!-----------------------------------------------------------------------
!
      CALL mpi_barrier (MPI_COMM_WORLD, MyError)
      CALL MCTWorld_clean ()
      CALL mpi_finalize (MyError)

      STOP

      END PROGRAM mct_driver
