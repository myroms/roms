#include "cppdefs.h"
      PROGRAM air_ocean
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2007 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!======================================================================= 
!                                                                      !
!  Wheather Research and Forcasting (WRF) model, Version 1.3           !
!                                                                      !
!     http://www.wrf-model.org                                         !
!                                                                      ! 
!  Regional Ocean Model System (ROMS), Version 3.0                     !
!  Terrain-following Ocean Model System (TOMS), Version 3.0            !
!                                                                      !
!     http://marine.rutgers.edu/po/index.php?model=roms                !
!     http://marine.rutgers.edu/po/index.php?model=toms                !
!                                                                      !
!  Master program to execute WRF and ROMS/TOMS in sequential or        !
!  concurrent modes.                                                   !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_scalars
!
      USE ocean_control_mod, ONLY : roms_init => initialize
      USE ocean_control_mod, ONLY : roms_run => run
      USE ocean_control_mod, ONLY : roms_finalize => finalize
!
      implicit none

#ifdef MPI
      include 'mpif.h'
#endif
!
!  Local variable declarations.
!
      logical, save :: first

      integer :: MyColor, MyCOMM, MyError, MyKey, MyRank, MyValue

      character (len=132) :: MyString
!
!-----------------------------------------------------------------------
!  Initialize distributed-memory (MPI) configuration
!-----------------------------------------------------------------------
!
!  Initialize MPI execution environment.
! 
      CALL mpi_init (MyError)
      CALL wrf_termio_dup
!
!  Get rank of the local process in the group associated with the
!  comminicator.
!
      CALL mpi_comm_rank (MPI_COMM_WORLD, MyRank, MyError)
!
!  Read in atmosphere-ocean coupling parameters from standard input.
!  Set temporarily the ocean communicator to current handle before
!  splitting so the input coupling script name can be broadcasted to
!  all the nodes.
!
      OCN_COMM_WORLD=MPI_COMM_WORLD
!
      CALL read_CouplePar (iNLM)
!
!  Split the communicator into WRF and ROMS subgroups based on color
!  and key.
!
      MyKey=0
      IF (MyRank.le.peATM_last) THEN
        MyColor=1
        MyString="COMPONENT_ID=1,COMPONENT_NAME=wrf"
      ELSE
        MyColor=2
        MyString="COMPONENT_ID=2,COMPONENT_NAME=roms"
      END IF
      CALL ext_mct_ioinit (MyString, MyError)
      CALL mpi_comm_split (MPI_COMM_WORLD, MyColor, MyKey, MyCOMM,      &
     &                     MyError)
!
!-----------------------------------------------------------------------
!  Run either WRF or ROMS according to the processor rank.  Notice that
!  in ensemble forecasting, ROMS is run over an ensemble loop. Also, in
!  variational data assimilation ROMS is run over outer and inner loops.
!  This requires a different code structure than the simple one below.
!  For now, the outside loop is deactivated and "Nrun" is set to unity.
!
!  In ensemble forecasting, a full atmosphere-ocean coupling is possible
!  but each member of the ensemble needs to be run on different parallel
!  nodes. Variational data assimilation (4DVAR) is more complicated and
!  requires more thinking.
!-----------------------------------------------------------------------
!
      IF (MyRank.le.peATM_last) THEN
        CALL wrf_init (MyCOMM)
        CALL wrf_run (nATM_steps, MyValue)
        CALL wrf_finalize
      ELSE
        first=.TRUE.
        Nrun=1
        IF (exit_flag.eq.NoError) THEN
          CALL roms_init (first, MyCOMM)
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
      CALL ext_mct_ioexit (MyError)
      CALL mpi_finalize (MyError)

      END PROGRAM air_ocean
