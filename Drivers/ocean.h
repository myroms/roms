#include "cppdefs.h"
      PROGRAM ocean
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  Regional Ocean Model System (ROMS), Version 2.1                     !
!  Terrain-following Ocean Model System (TOMS), Version 2.1            !
!                                                                      !
!  Master program to execute ROMS/TOMS in ocean mode only without      !
!  coupling (sequential or concurrent) to any atmospheric model.       !
!                                                                      !
!  This ocean model solves the free surface, hydrostatic, primitive    !
!  equations  over  variable  topography  using  stretched terrain-    !
!  following coordinates in the vertical and orthogonal curvilinear    !
!  coordinates in the horizontal.                                      !
!                                                                      !
!  Developers:                                                         !
!                                                                      !
!  Dr. Hernan G. Arango                                                !
!    Institute of Marine and Coastal Sciences                          !
!    Rutgers University, New Brunswick, NJ, USA                        !
!    (arango@imcs.rutgers.edu)                                         !
!                                                                      !
!  Dr. Alexander F. Shchepetkin                                        !
!    Institute of Geophysics and Planetary Physics                     !
!    UCLA, Los Angeles, CA, USA                                        !
!    (alex@atmos.ucla.edu)                                             !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      USE ocean_control_mod, ONLY : initialize
      USE ocean_control_mod, ONLY : run
      USE ocean_control_mod, ONLY : finalize
!
      implicit none
!
!  Local variable declarations.
!
      logical, save :: first

      integer :: MyError

#ifdef DISTRIBUTE
# ifdef MPI
!
!-----------------------------------------------------------------------
!  Initialize distributed-memory MPI configuration.
!-----------------------------------------------------------------------
!
      CALL mpi_init (MyError)
      IF (MyError.ne.0) THEN
        WRITE (stdout,10)
  10    FORMAT (/,' ROMS/TOMS - Unable to initialize MPI.')
        exit_flag=6
        STOP
      END IF
!
!  Get rank of the local process in the group associated with the
!  comunicator.
!
      CALL mpi_comm_rank (MPI_COMM_WORLD, MyRank, MyError)
      IF (MyError.ne.0) THEN
        WRITE (stdout,20)
  20    FORMAT (/,' ROMS/TOMS - Unable to inquire rank of local',       &
     &              ' processor.')
        exit_flag=6
        STOP
      END IF
# endif
#endif
!
!-----------------------------------------------------------------------
!  Initialize internal and external parameters and state variables.
!-----------------------------------------------------------------------
!
      first=.TRUE.
      CALL initialize (first)
!
!-----------------------------------------------------------------------
!  Time-step ocean model once or over an ensemble loop, if any.
!-----------------------------------------------------------------------
!
      DO Nrun=ERstr,ERend
        CALL run
      END DO
!
!-----------------------------------------------------------------------
!  Terminate model execution: flush and close all IO files.
!-----------------------------------------------------------------------
!
      CALL finalize
#if defined DISTRIBUTE && defined MPI
      CALL mpi_finalize (MyError)
#endif

      END PROGRAM ocean
