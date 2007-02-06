#include "cppdefs.h"

      MODULE ocean_control_mod
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2007 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  ROMS/TOMS Adjoint Sensitivity Analysis Driver:                      !
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
      USE mod_iounits
      USE mod_scalars

#ifdef AIR_OCEAN 
!
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
      logical :: allocate_vars = .TRUE.

      integer :: ng, thread

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
!  Read in model tunable parameters from standard input.
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
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
!
!  Local variable declarations.
!
      integer :: i, my_iic, ng

      real (r8) :: str_day, end_day
!
!=======================================================================
!  Run model for all nested grids, if any.
!=======================================================================
!
      NEST_LOOP : DO ng=1,Ngrids
!
!  Initialize adjoint model and define sensitivity functional.
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
!  Time-step adjoint model: Compute the gradient or index, dJ/dS, of
!  the sensitivity functional.
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

      END DO NEST_LOOP
!
 10   FORMAT (/,a,i2,/)
 20   FORMAT (/,'AD ROMS/TOMS: started time-stepping:',                 &
     &        ' (TimeSteps: ',i8.8,' - ',i8.8,')',/,14x,                &
     &        'adjoint forcing time range: ',f12.4,' - ',f12.4 ,/)
 30   FORMAT (/,' Out of range adjoint forcing time, ',a,f12.4,/,       &
     &        ' It must be between ',f12.4,' and ',f12.4)

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
