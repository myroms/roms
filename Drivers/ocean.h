#include "cppdefs.h"
      PROGRAM ocean
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  Regional Ocean Model System (ROMS), Version 2.0                     !
!  Terrain-following Ocean Model System (TOMS), Version 2.0            !
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
      USE mod_ncparam
      USE mod_scalars
!
      implicit none

      logical :: allocate_vars = .true.

      integer :: my_iic, ng
!
!-----------------------------------------------------------------------
!  Initialize model internal parameters.
!-----------------------------------------------------------------------
!
      CALL initialize_param
      CALL initialize_parallel
      CALL initialize_scalars
!
!-----------------------------------------------------------------------
!  Read in model tunable parameters from standard input.
!-----------------------------------------------------------------------
!
      CALL inp_par
!
!-----------------------------------------------------------------------
!  Allocate and initialize model variables for each nested grid..
!-----------------------------------------------------------------------
!
      CALL mod_arrays (allocate_vars)
!
!-----------------------------------------------------------------------
!  Run model for all nested grids, if any.  The model also can be run
!  once or over an ensemble loop.
!-----------------------------------------------------------------------
!
      ENSEMBLE : DO Nrun=ERstr,ERend
!
!  Initialize model state variables.
!
        DO ng=1,Ngrids
          CALL initial (ng)
          IF (exit_flag.ne.0) THEN
            write(stdout,10) Rerror(exit_flag), exit_flag
 10         format(/,a,i3,/)
            EXIT ENSEMBLE
          END IF
        END DO
!
!  Time step ocean model.
!
        NEST : DO ng=1,Ngrids
          IF (Master) WRITE (stdout,20) ntstart, ntend
 20       FORMAT (/,' ROMS/TOMS: started time-stepping:',               &
     &              '( TimeSteps: ',i8.8,' - ',i8.8,')',/)
!
          time(ng)=time(ng)-dt(ng)
          STEP : DO my_iic=ntstart,ntend+1
            iic(ng)=my_iic
#ifdef SOLVE3D
            CALL main3d (ng)
#else
            CALL main2d (ng)
#endif
            IF (exit_flag.ne.0) THEN
              write(stdout,10) Rerror(exit_flag), exit_flag
              exit NEST
            END IF
          END DO STEP
        END DO NEST
      END DO ENSEMBLE
!
!-----------------------------------------------------------------------
!  Close IO files.
!-----------------------------------------------------------------------
!
      CALL close_io
      STOP
      END PROGRAM ocean
