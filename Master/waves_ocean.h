      PROGRAM waves_ocean
!
!==================================================== John C. Warner ===
!  Copyright (c) 2006 ROMS/TOMS                                        !
!================================================== Hernan G. Arango ===
!                                                                      !
!  Master program to couple Atmosphere, waves, ocean models.           !
!                                                                      !
!  Atmophere Model:  to be determined.                                 !
!                                                                      !
!  Waves Model: SWAN (Simulating WAves Nearshore), Version 40.41AB     !
!                                                                      !
!     http://vlm089.citg.tudelft.nl/swan/index.htm                     !
!                                                                      ! 
!  Ocean Model: ROMS (Regional Ocean Model System), Version 2.3        !
!                                                                      !
!     http://marine.rutgers.edu/roms                                   !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      USE ocean_control_mod, ONLY : roms_initialize => initialize
      USE ocean_control_mod, ONLY : roms_run => run
      USE ocean_control_mod, ONLY : roms_finalize => finalize
!
      implicit none
!
!  Local variable declarations.
!
      logical, save :: first

      integer :: MyColor, MyCOMM, MyError, MyKey, MySize, MyValue
!      integer :: peOCN_frst, peOCN_last, peWAV_frst, peWAV_last
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
      CALL mpi_comm_size (MPI_COMM_WORLD, MySize, MyError)
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
!  Assign processors to ocean and wave models.
!
      peOCN_frst=0
      peOCN_last=peOCN_frst+NnodesOCN-1
      peWAV_frst=peOCN_last+1
      peWAV_last=peWAV_frst+NnodesWAV-1
      IF (peWAV_last.ne.MySize-1) THEN
        IF (Master) THEN
          WRITE (stdout,10) peWAV_last, MySize
 10       FORMAT (/,' WAVES_OCEAN - Number assigned processors: ',i3.3, &
     &            /,15x,'not equal to spawned MPI nodes: ',i3.3)
        END IF
        STOP
      ELSE
         WRITE (stdout,20) peOCN_frst, peOCN_last,                      &
     &                     peWAV_frst, peWAV_last
 20      FORMAT (/,' Waves-Ocean Models Coupling: ',/,                  &
     &           /,7x,'Ocean Model MPI nodes: ',i3.3,' - ', i3.3,/,     &
     &           /,7x,'Waves Model MPI nodes: ',i3.3,' - ', i3.3)
      END IF
!
!  Split the communicator into SWAN and ROMS subgroups based on color
!  and key.
!
      MyKey=0
      IF ((peOCN_frst.le.MyRank).and.(MyRank.le.peOCN_last)) THEN
        MyColor=OCNid
      END IF
      IF ((peWAV_frst.le.MyRank).and.(MyRank.le.peWAV_last)) THEN
        MyColor=WAVid
      END IF
      CALL mpi_comm_split (MPI_COMM_WORLD, MyColor, MyKey, MyCOMM,      &
     &                     MyError)
!
!-----------------------------------------------------------------------
!  Run either SWAN or ROMS according to the processor rank.  Notice that
!  in ensemble forecasting, ROMS is run over an ensemble loop. Also, in
!  variational data assimilation ROMS is run over outer and inner loops.
!  This requires a different code structure than the simple one below.
!  For now, the outside loop is deactivated and "Nrun" is set to unity.
!
!  In ensemble forecasting, a full waves-ocean coupling is possible
!  but each member of the ensemble needs to be run on different parallel
!  nodes. Variational data assimilation (4DVAR) is more complicated and
!  requires more thinking.
!-----------------------------------------------------------------------
!
      IF (MyColor.eq.WAVid) THEN
        CALL SWINITMPI (MyCOMM)
        CALL SWMAIN (REAL(TI_WAV_OCN))
        CALL SWEXITMPI
      END IF
      IF (MyColor.eq.OCNid) THEN
        first=.TRUE.
        Nrun=1
        IF (exit_flag.eq.NoError) THEN
          CALL roms_initialize (first, MyCOMM)
        END IF
        IF (exit_flag.eq.NoError) THEN
          CALL roms_run
        END IF
        CALL roms_finalize
      END IF
!
!-----------------------------------------------------------------------
!  Terminates all the MPI processing.
!-----------------------------------------------------------------------
!
      CALL mpi_finalize (MyError)

      END PROGRAM waves_ocean
