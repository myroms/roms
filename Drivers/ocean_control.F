#include "cppdefs.h"
      MODULE ocean_control_mod
!
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002 ROMS/TOMS Group              Christopher Moore   !
!=======================================================================
!                                                                      !
!  This subroutines control the initialization, time-stepping, and     !
!  finalization of ROMS/TOMS execution following ESMF conventions:     !
!                                                                      !
!     initialize                                                       !
!     run                                                              !
!     finalize                                                         !
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
!  This routine initialize ROMS/TOMS internal and external parameters  !
!  and state variables.                                                !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
#ifdef AIR_OCEAN 
      USE atm_coupler_mod, ONLY : initialize_coupling
#endif
!
!  Imported variable declarations.
!
      logical, intent(inout) :: first

      integer, intent(in), optional :: MyCOMM
!
!  Local variable declarations.
!
      logical :: allocate_vars = .true.

      integer :: ng

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
!  Initialize model internal parameters.
!
        CALL initialize_param
        CALL initialize_parallel
        CALL initialize_scalars

#ifdef AIR_OCEAN 
!
!  Initialize coupling streams between atmosphere and ocean using the
!  Model Coupling Toolkit (MCT).
!
        CALL initialize_coupling (MyRank)
#endif
!
!  Read in model tunable parameters from standard input.
!
        CALL inp_par
!
!  Allocate and initialize model variables.
!
        CALL init_arrays (allocate_vars)
      END IF
!
!-----------------------------------------------------------------------
!  Initialize model state variables for all nested/composed grids.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        CALL initial (ng)
        IF (exit_flag.ne.0) THEN
          WRITE (stdout,10) Rerror(exit_flag), exit_flag
 10       FORMAT (/,a,i3,/)
          RETURN
        END IF
      END DO

      RETURN
      END SUBROUTINE initialize

      SUBROUTINE run
!
!=======================================================================
!                                                                      !
!  This routine time-steps ocean model ROMS/TOMS.                      !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
!  Local variable declarations.
!
      integer :: ng, my_iic
!
!-----------------------------------------------------------------------
!  Run model for all nested grids, if any.
!-----------------------------------------------------------------------
!
      NEST : DO ng=1,Ngrids

        IF (Master) WRITE (stdout,20) ntstart, ntend
 20     FORMAT (/,' ROMS/TOMS: started time-stepping:',               &
     &            '( TimeSteps: ',i8.8,' - ',i8.8,')',/)
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
            WRITE (stdout,30) Rerror(exit_flag), exit_flag
 30         FORMAT (/,a,i3,/)
            RETURN
          END IF

        END DO STEP

      END DO NEST

      RETURN
      END SUBROUTINE run

      SUBROUTINE finalize
!
!=======================================================================
!                                                                      !
!  This routine terminates ocean model ROMS/TOMS execution.            !
!                                                                      !
!=======================================================================
!
!  Close IO files.
!
      CALL close_io

      RETURN
      END SUBROUTINE finalize

      END MODULE ocean_control_mod
