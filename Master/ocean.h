      PROGRAM ocean
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2007 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!======================================================================= 
!                                                                      !
!  Regional Ocean Model System (ROMS)                                  !
!  Terrain-following Ocean Model System (TOMS)                         !
!                                                                      !
!  Master program to execute  ROMS/TOMS  drivers in ocean mode only    !
!  without coupling (sequential or concurrent) to  any  atmospheric    !
!  model.                                                              !
!                                                                      !
!  This ocean model solves the free surface, hydrostatic, primitive    !
!  equations  over  variable  topography  using  stretched terrain-    !
!  following coordinates in the vertical and orthogonal curvilinear    !
!  coordinates in the horizontal.                                      !
!                                                                      !
!  Nonlinear Model Developers:                                         !
!                                                                      !
!  Dr. Hernan G. Arango                                                !
!    Institute of Marine and Coastal Sciences                          !
!    Rutgers University, New Brunswick, NJ, USA                        !
!    (arango@marine.rutgers.edu)                                       !
!                                                                      !
!  Dr. Alexander F. Shchepetkin                                        !
!    Institute of Geophysics and Planetary Physics                     !
!    UCLA, Los Angeles, CA, USA                                        !
!    (alex@atmos.ucla.edu)                                             !
!                                                                      !
!  Tangent linear and Adjoint Models Developers:                       !
!                                                                      !
!    Dr. Hernan G. Arango    (arango@marine.rutgers.edu)               !
!    Dr. Bruce Cornuelle     (bcornuelle@ucsd.edu)                     !
!    Dr. Emanuele Di Lorenzo (edl@eas.gatech.edu)                      !
!    Dr. Arthur J. Miller    (ajmiller@ucsd.edu)                       !
!    Dr. Andrew M. Moore     (ammoore@ucsc.edu)                        !
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
      END IF
# endif
#endif
!
!-----------------------------------------------------------------------
!  Initialize ocean internal and external parameters and state
!  variables.
!-----------------------------------------------------------------------
!
      IF (exit_flag.eq.NoError) THEN
        first=.TRUE.
        CALL initialize (first)
      END IF
!
!-----------------------------------------------------------------------
!  Time-step ocean model.
!-----------------------------------------------------------------------
!
      IF (exit_flag.eq.NoError) THEN
        CALL run
      END IF
!
!-----------------------------------------------------------------------
!  Terminate ocean model execution: flush and close all IO files.
!-----------------------------------------------------------------------
!
      CALL finalize
#if defined DISTRIBUTE && defined MPI
      CALL mpi_finalize (MyError)
#endif

      END PROGRAM ocean
