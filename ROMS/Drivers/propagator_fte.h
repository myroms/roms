      SUBROUTINE propagator (RunInterval, state, tl_state)
!
!svn $Id$
!************************************************** Hernan G. Arango ***
!  Copyright (c) 2002-2012 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!***********************************************************************
!                                                                      !
!  Finite Time Eigenvalues Propagator:                                 !
!                                                                      !
!  This routine it is used to compute the  Finite  Time  Eigenmodes    !
!  (FTEs) of the propagator R(0,t) linearized about a time evolving    !
!  circulation.  It involves a  single integration  of an arbitrary    !
!  perturbation state vector forward in time over the interval [0,t]   !
!  by the tangent linear model.                                        !
!                                                                      !
!   Reference:                                                         !
!                                                                      !
!     Moore, A.M. et al., 2004: A comprehensive ocean prediction and   !
!       analysis system based on the tangent linear and adjoint of a   !
!       regional ocean model, Ocean Modelling, 7, 227-258.             !
!                                                                      !
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
#ifdef SOLVE3D
      USE mod_coupling
#endif
      USE mod_iounits
      USE mod_ocean
      USE mod_scalars
      USE mod_stepping
!
      USE dotproduct_mod, ONLY : tl_statenorm
      USE packing_mod, ONLY : tl_unpack, tl_pack
#ifdef SOLVE3D
      USE set_depth_mod, ONLY: set_depth
#endif
!
!  Imported variable declarations.
!
      real(r8), intent(in) :: RunInterval

      TYPE (T_GST), intent(in) :: state(Ngrids)
      TYPE (T_GST), intent(inout) :: tl_state(Ngrids)
!
!  Local variable declarations.
!
#ifdef SOLVE3D
      logical :: FirstPass = .TRUE.
#endif
      integer :: ng, subs, tile, thread

      real(r8) :: StateNorm(Ngrids)
!
!=======================================================================
!  Forward integration of the tangent linear model.
!=======================================================================
!
      Nrun=Nrun+1
      IF (Master) THEN
        DO ng=1,Ngrids
          WRITE (stdout,10) ' PROPAGATOR - Grid: ', ng,                 &
     &                      ',  Iteration: ', Nrun,                     &
     &                      ',  number converged RITZ values: ',        &
     &                      Nconv(ng)
        END DO
      END IF
!
!  Initialize time stepping indices and counters.
!
      DO ng=1,Ngrids
        iif(ng)=1
        indx1(ng)=1
        kstp(ng)=1
        krhs(ng)=1
        knew(ng)=1
        PREDICTOR_2D_STEP(ng)=.FALSE.
        synchro_flag(ng)=.TRUE.
!
        iic(ng)=0
        nstp(ng)=1
        nrhs(ng)=1
        nnew(ng)=1
!
        tdays(ng)=dstart
        time(ng)=tdays(ng)*day2sec
        ntstart(ng)=INT((time(ng)-dstart*day2sec)/dt(ng))+1
        ntend(ng)=ntimes(ng)
        ntfirst(ng)=ntstart(ng)
      END DO
!
!-----------------------------------------------------------------------
!  Clear tangent linear state variables. There is not need to clean
!  the basic state arrays since they were zeroth out at initialization
!  and bottom of previous iteration.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1,+1
            CALL initialize_ocean (ng, TILE, iTLM)
          END DO
        END DO
!$OMP END PARALLEL DO
      END DO

#ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Compute basic state initial level thicknesses used for state norm
!  scaling. It uses zero time averaged free-surface (rest state).
!  Therefore, the norm scaling is time invariant.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*(thread+1)-1,subs*thread,-1
            CALL set_depth (ng, TILE)
          END DO
        END DO
!$OMP END PARALLEL DO
      END DO
#endif
!
!-----------------------------------------------------------------------
!  Unpack tangent linear initial conditions from state vector.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
!$OMP PARALLEL DO PRIVATE(thread,subs,tile)                             &
!$OMP&            SHARED(numthreads,Nstr,Nend,state)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1,+1
            CALL tl_unpack (ng, TILE, Nstr(ng), Nend(ng),               &
     &                      state(ng)%vector)
          END DO
        END DO
!$OMP END PARALLEL DO
      END DO
!
!-----------------------------------------------------------------------
!  Compute initial tangent linear state dot product norm.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
!$OMP PARALLEL DO PRIVATE(thread,subs,tile)                             &
!$OMP&            SHARED(numthreads,krhs,nstp,StateNorm)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*(thread+1)-1,subs*thread,-1
            CALL tl_statenorm (ng, TILE, kstp(ng), nstp(ng),            &
     &                         StateNorm(ng))
          END DO
        END DO
!$OMP END PARALLEL DO

        IF (Master) THEN
          WRITE (stdout,20) ' PROPAGATOR - Grid: ', ng,                 &
     &                      ',  Tangent Initial Norm: ', StateNorm(ng)
        END IF
      END DO
!
!-----------------------------------------------------------------------
!  Read in initial forcing, climatology and assimilation data from
!  input NetCDF files.  It loads the first relevant data record for
!  the time-interpolation between snapshots.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        CALL close_inp (ng, iTLM)
        IF (exit_flag.ne.NoError) RETURN
#ifdef TIMELESS_DATA
        CALL tl_get_idata (ng)
        IF (exit_flag.ne.NoError) RETURN
#endif
        CALL tl_get_data (ng)
        IF (exit_flag.ne.NoError) RETURN
      END DO
!
!-----------------------------------------------------------------------
!  Time-step the tangent linear model.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        IF (Master) THEN
          WRITE (stdout,30) 'TL', ng, ntstart(ng), ntend(ng)
        END IF
        iic(ng)=ntstart(ng)-1
        time(ng)=time(ng)-dt(ng)
      END DO

#ifdef SOLVE3D
      CALL tl_main3d (RunInterval)
#else
      CALL tl_main2d (RunInterval)
#endif
      IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Clear nonlinear state (basic state) variables and insure that the
!  time averaged free-surface is zero for scaling below and next
!  iteration.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1,+1
            CALL initialize_ocean (ng, TILE, iNLM)
#ifdef SOLVE3D
            CALL initialize_coupling (ng, TILE, 0)
#endif
          END DO
        END DO
!$OMP END PARALLEL DO
      END DO

#ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Compute basic state final level thicknesses used for state norm
!  scaling. It uses zero time averaged free-surface (rest state).
!  Therefore, the norm scaling is time invariant.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*(thread+1)-1,subs*thread,-1
            CALL set_depth (ng, TILE)
          END DO
        END DO
!$OMP END PARALLEL DO
      END DO
#endif
!
!-----------------------------------------------------------------------
!  Compute final tangent linear state dot product norm.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
!$OMP PARALLEL DO PRIVATE(thread,subs,tile)                             &
!$OMP&            SHARED(numthreads,krhs,nstp,StateNorm)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1,+1
            CALL tl_statenorm (ng, TILE, knew(ng), nstp(ng),            &
     &                         StateNorm(ng))
          END DO
        END DO
!$OMP END PARALLEL DO

        IF (Master) THEN
          WRITE (stdout,20) ' PROPAGATOR - Grid: ', ng,                 &
     &                      ',  Tangent   Final Norm: ', StateNorm(ng)
        END IF
      END DO
!
!-----------------------------------------------------------------------
!  Pack final tangent linear solution into tangent state vector.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
!$OMP PARALLEL DO PRIVATE(thread,subs,tile)                             &
!$OMP&            SHARED(numthreads,Nstr,Nend,tl_state)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*(thread+1)-1,subs*thread,-1
            CALL tl_pack (ng, TILE, Nstr(ng), Nend(ng),                 &
     &                    tl_state(ng)%vector)
          END DO
        END DO
!$OMP END PARALLEL DO
      END DO
!
 10   FORMAT (/,a,i2.2,a,i3.3,a,i3.3/)
 20   FORMAT (/,a,i2.2,a,1p,e15.6,/)
 30   FORMAT (/,1x,a,1x,'ROMS/TOMS: started time-stepping:',            &
     &        ' (Grid: ',i2.2,' TimeSteps: ',i8.8,' - ',i8.8,')')

      RETURN
      END SUBROUTINE propagator
