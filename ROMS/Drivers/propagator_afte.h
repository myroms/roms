      SUBROUTINE propagator (RunInterval, state, ad_state)
!
!svn $Id$
!************************************************** Hernan G. Arango ***
!  Copyright (c) 2002-2011 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!***********************************************************************
!                                                                      !
!  Adjoint Finite Time Eigenvalues Propagator:                         !
!                                                                      !
!  This routine is used during the computation of the eigenvectors of  !
!  the adjoint propagator, transpose[R(t,0)]. They are computed in an  !
!  analogous way to those of  R(t,0).  A  single  integration  of  an  !
!  arbitrary perturbation state vector  "u" backward in time over the  !
!  interval [t,0] by the adjoint model: transpose[R(t,0)]*u.           !
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
      USE dotproduct_mod, ONLY : ad_statenorm
      USE packing_mod, ONLY : ad_unpack, ad_pack
#ifdef SOLVE3D
      USE set_depth_mod, ONLY: set_depth
#endif
!
!  Imported variable declarations.
!
      real(r8), intent(in) :: RunInterval

      TYPE (T_GST), intent(in) :: state(Ngrids)
      TYPE (T_GST), intent(inout) :: ad_state(Ngrids)
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
        krhs(ng)=3
        knew(ng)=2
        PREDICTOR_2D_STEP(ng)=.FALSE.
        synchro_flag(ng)=.TRUE.
!
        iic(ng)=0
        nstp(ng)=1
        nrhs(ng)=1
        nnew(ng)=2
!
        tdays(ng)=dstart+dt(ng)*FLOAT(ntimes(ng))*sec2day
        time(ng)=tdays(ng)*day2sec
        ntstart(ng)=ntimes(ng)+1
        ntend(ng)=1
        ntfirst(ng)=ntend(ng)
      END DO
!
!-----------------------------------------------------------------------
!  Clear adjoint state variables.  There is not need to clean the basic
!  state arrays since they were zeroth out at initialization and bottom
!  of previous iteration.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1,+1
            CALL initialize_ocean (ng, TILE, iADM)
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
!  Unpack adjoint initial conditions from state vector.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
!$OMP PARALLEL DO PRIVATE(thread,subs,tile)                             &
!$OMP&            SHARED(numthreads,Nstr,Nend,state)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1,+1
            CALL ad_unpack (ng, TILE, Nstr(ng), Nend(ng),               &
     &                      state(ng)%vector)
          END DO
        END DO
!$OMP END PARALLEL DO
      END DO
!
!-----------------------------------------------------------------------
!  Compute initial adjoint state dot product norm.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
!$OMP PARALLEL DO PRIVATE(thread,subs,tile)                             &
!$OMP&            SHARED(numthreads,krhs,nstp,StateNorm)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*(thread+1)-1,subs*thread,-1
            CALL ad_statenorm (ng, TILE, knew(ng), nstp(ng),            &
     &                         StateNorm(ng))
          END DO
        END DO
!$OMP END PARALLEL DO

        IF (Master) THEN
          WRITE (stdout,20) ' PROPAGATOR - Grid: ', ng,                 &
     &                      ',  Adjoint Initial Norm: ', StateNorm(ng)
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
        CALL close_inp (ng, iADM)
        IF (exit_flag.ne.NoError) RETURN
#ifdef TIMELESS_DATA
        CALL ad_get_idata (ng)
        IF (exit_flag.ne.NoError) RETURN
#endif
        CALL ad_get_data (ng)
        IF (exit_flag.ne.NoError) RETURN
      END DO
!
!-----------------------------------------------------------------------
!  Time-step the adjoint model backwards.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        IF (Master) THEN
          WRITE (stdout,30) 'AD', ng, ntstart(ng), ntend(ng)
        END IF
        iic(ng)=ntstart(ng)+1
        time(ng)=time(ng)+dt(ng)
      END DO

#ifdef SOLVE3D
      CALL ad_main3d (RunInterval)
#else
      CALL ad_main2d (RunInterval)
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
!  Compute final adjoint state dot product norm.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
!$OMP PARALLEL DO PRIVATE(thread,subs,tile)                             &
!$OMP&            SHARED(numthreads,krhs,nstp,StateNorm)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1,+1
            CALL ad_statenorm (ng, TILE, knew(ng), nstp(ng),            &
     &                         StateNorm(ng))
          END DO
        END DO
!$OMP END PARALLEL DO

        IF (Master) THEN
          WRITE (stdout,20) ' PROPAGATOR - Grid: ', ng,                 &
     &                      ',  Adjoint   Final Norm: ', StateNorm(ng)
        END IF
      END DO
!
!-----------------------------------------------------------------------
!  Pack final adjoint solution into adjoint state vector.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
!$OMP PARALLEL DO PRIVATE(thread,subs,tile)                             &
!$OMP&            SHARED(numthreads,Nstr,Nend,ad_state)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*(thread+1)-1,subs*thread,-1
            CALL ad_pack (ng, TILE, Nstr(ng), Nend(ng),                 &
     &                    ad_state(ng)%vector)
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
