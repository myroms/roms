      MODULE ocean_control_mod
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  ROMS/TOMS Tangent Linear Model Linearization Test:                  !
!                                                                      !
!  This driver is used to check the "linearization" of the tangent     !
!  linear model using a structure similar to the gradient test.        !
!                                                                      !
!  The subroutines in this driver control the initialization, time-    !
!  stepping,  and  finalization of ROMS/TOMS  model following ESMF     !
!  conventions:                                                        !
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
!  This routine allocates and initializes ROMS/TOMS state variables    !
!  and internal and external parameters.                               !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_iounits
      USE mod_scalars
!
#ifdef AIR_OCEAN 
      USE ocean_coupler_mod, ONLY : initialize_atmos_coupling
#endif
#ifdef WAVES_OCEAN
      USE ocean_coupler_mod, ONLY : initialize_waves_coupling
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
!  Initialize parallel parameters.
!
        CALL initialize_parallel
!
!  Initialize wall clocks.
!
        IF (Master) THEN
          WRITE (stdout,10)
 10       FORMAT (' Process Information:',/)
        END IF
        DO ng=1,NGRIDS
!$OMP PARALLEL DO PRIVATE(thread) SHARED(ng,numthreads)
          DO thread=0,numthreads-1
            CALL wclock_on (ng, iNLM, 0)
          END DO
!$OMP END PARALLEL DO
        END DO

#if defined AIR_OCEAN || defined WAVES_OCEAN
!
!  Initialize coupling streams between model(s).
!
        DO ng=1,Ngrids
# ifdef AIR_OCEAN
          CALL initialize_atmos_coupling (ng, MyRank)
# endif
# ifdef WAVES_OCEAN
          CALL initialize_waves_coupling (ng, MyRank)
# endif
        END DO
#endif
!
!  Read in model tunable parameters from standard input. Initialize
!  "mod_param", "mod_ncparam" and "mod_scalar" modules.
!
        CALL inp_par (iNLM)
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
!
!  Allocate and initialize observation arrays.
!
        CALL initialize_fourdvar
      END IF

      RETURN
      END SUBROUTINE initialize

      SUBROUTINE run
!
!=======================================================================
!                                                                      !
!  This routine time-steps ROMS/TOMS nonlinear, tangent linear and     !
!  adjoint models.                                                     !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
      USE mod_stepping
!
      USE dotproduct_mod, ONLY : ad_dotproduct
!
!  Local variable declarations.
!
      integer :: IniRec, i, lstr, my_iic, ng, status
      integer :: subs, tile, thread

      real(r8) :: gp, hp, p
!
!=======================================================================
!  Run model for all nested grids, if any.
!=======================================================================
!
      NEST_LOOP : DO ng=1,Ngrids
!
!-----------------------------------------------------------------------
!  Run nonlinear and adjoint models.
!-----------------------------------------------------------------------
!
!  Initialize relevant parameters.
!
        Lold(ng)=1
        Lnew(ng)=1
        Nrun=1
        Ipass=1
        ERstr=1
        ERend=NstateVar(ng)+1
        IniRec=1
        ig1count=0
        ig2count=0
        nTLM(ng)=nHIS(ng)                      ! to allow IO comparison
        LcycleTLM=.FALSE.
!
!  Initialize nonlinear model with first guess initial conditions.
!
        CALL initial (ng)
        IF (exit_flag.ne.NoError) THEN
          IF (Master) THEN
            WRITE (stdout,10) Rerror(exit_flag), exit_flag
          END IF
          RETURN
        END IF
!
!  Run nonlinear model. Extract and store nonlinear model values at
!  observation locations.
!
        IF (Master) THEN
          WRITE (stdout,20) 'NL', ntstart(ng), ntend(ng)
        END IF

        wrtNLmod(ng)=.TRUE.
        wrtTLmod(ng)=.FALSE.
        time(ng)=time(ng)-dt(ng)

        NL_LOOP1 : DO my_iic=ntstart(ng),ntend(ng)+1

          iic(ng)=my_iic
#ifdef SOLVE3D
          CALL main3d (ng)
#else
          CALL main2d (ng)
#endif
          IF (exit_flag.ne.NoError) THEN
            IF (Master) THEN
              WRITE (stdout,10) Rerror(exit_flag), exit_flag
            END IF
            RETURN
          END IF

        END DO NL_LOOP1
!
!  Close current nonlinear model history file.
!
        status=nf90_close(ncHISid(ng))
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,30) TRIM(HISname(ng))
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        wrtNLmod(ng)=.FALSE.
        wrtTLmod(ng)=.TRUE.
!
!  Save and Report cost function between nonlinear model and
!  observations.
!
        DO i=0,NstateVar(ng)
          FOURDVAR(ng)%CostFunOld(i)=FOURDVAR(ng)%CostFun(i)
        END DO
        IF (Master) THEN        
          WRITE (stdout,40) FOURDVAR(ng)%CostFunOld(0)
          DO i=1,NstateVar(ng)
            IF (FOURDVAR(ng)%CostFunOld(i).gt.0.0_r8) THEN
              IF (i.eq.1) THEN
                WRITE (stdout,50) FOURDVAR(ng)%CostFunOld(i),           &
     &                            TRIM(Vname(1,idSvar(i)))
              ELSE
                WRITE (stdout,60) FOURDVAR(ng)%CostFunOld(i),           &
     &                            TRIM(Vname(1,idSvar(i)))
              END IF
            END IF
          END DO
        END IF
!
!  Initialize the adjoint model from rest.
!
        CALL ad_initial (ng)
        IF (exit_flag.ne.NoError) THEN
          IF (Master) THEN
            WRITE (stdout,10) Rerror(exit_flag), exit_flag
          END IF
          RETURN
        END IF
!
!  Time-step adjoint model: Compute model state gradient, GRAD(J).
!  Force the adjoint model with the adjoint misfit between nonlinear
!  model and observations.
!
        IF (Master) THEN
          WRITE (stdout,20) 'AD', ntstart(ng), ntend(ng)
        END IF

        time(ng)=time(ng)+dt(ng)

        AD_LOOP : DO my_iic=ntstart(ng),ntend(ng),-1

          iic(ng)=my_iic
#ifdef SOLVE3D
          CALL ad_main3d (ng)
#else
          CALL ad_main2d (ng)
#endif
          IF (exit_flag.ne.NoError) THEN
            IF (Master) THEN
              WRITE (stdout,10) Rerror(exit_flag), exit_flag
            END IF
            RETURN
          END IF

        END DO AD_LOOP
!
!-----------------------------------------------------------------------
!  Perturb each tangent linear state variable using the steepest decent
!  direction by grad(J). If ndefTLM < 0, suppress IO for both nonlinear
!  and tangent linear models in the outer and inner loops.
!-----------------------------------------------------------------------
!
!  Load adjoint solution.
!
        CALL get_state (ng, iADM, 3, ADJname(ng), tADJindx(ng),         &
     &                  Lnew(ng))
!
!  Compute adjoint solution dot product for scaling purposes.
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile)                             &
!$OMP&            SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1
            CALL ad_dotproduct (ng, TILE, Lnew(ng))
          END DO
        END DO
!$OMP END PARALLEL DO

#ifdef SOLVE3D
!
!  OUTER LOOP: First, Perturb all state variables (zeta, u, v, t) at
!  ==========  once. Then, perturb one state variable at the time.
!              Notice, that ubar and vbar are not perturbed. They are
!              computed by vertically integarting u and v.
!
#else
!
!  OUTER LOOP: First, perturb all state variables (zeta, ubar, vbar) at
!  ==========  once. Then, perturb one state variable at the time.
!
#endif
        OUTER_LOOP : DO outer=ERstr,ERend

          CALL get_state (ng, iNLM, 1, FWDname(ng), IniRec, Lnew(ng))
!
!  INNER LOOP: scale perturbation amplitude by selecting "p" scalar,
!  ==========  such that:
!                              p = 10 ** FLOAT(-inner) 
!
          INNER_LOOP : DO inner=1,Ninner
!
!  Add a perturbation to the nonlinear state initial conditions
!  according with the outer and inner loop iterations. The added
!  term is a function of the steepest descent direction defined
!  by grad(J) times the perturbation amplitude "p".
!
            CALL initial (ng)
            IF (exit_flag.ne.NoError) THEN
              IF (Master) THEN
                WRITE (stdout,10) Rerror(exit_flag), exit_flag
              END IF
              RETURN
            END IF
            IF (Master) THEN
              lstr=LEN_TRIM(HISbase(ng))
              WRITE (HISname,70) HISbase(ng)(1:lstr-3), Nrun
            END IF
            IF (ndefTLM(ng).lt.0) THEN
              LdefHIS(ng)=.FALSE.              ! suppress IO
              LwrtHIS(ng)=.FALSE.
            ELSE
              LdefHIS(ng)=.TRUE.
            END IF
            wrtNLmod(ng)=.TRUE.
            wrtTLmod(ng)=.FALSE.
!
!  Time-step nonlinear model: compute perturbed nonlinear state.
!
            IF (Master) THEN
              WRITE (stdout,20) 'NL', ntstart(ng), ntend(ng)
            END IF

            time(ng)=time(ng)-dt(ng)

            NL_LOOP2 : DO my_iic=ntstart(ng),ntend(ng)+1

              iic(ng)=my_iic
#ifdef SOLVE3D
              CALL main3d (ng)
#else
              CALL main2d (ng)
#endif
              IF (exit_flag.ne.NoError) THEN
                IF (Master) THEN
                  WRITE (stdout,10) Rerror(exit_flag), exit_flag
                END IF
                RETURN
              END IF

            END DO NL_LOOP2
!
!  Get current nonlinear model trajectory.
!
            FWDname(ng)=HISbase(ng)
            CALL get_state (ng, iNLM, 1, FWDname(ng), IniRec, Lnew(ng))
!
!  Initialize tangent linear with the steepest decent direction
!  (adjoint state, GRAD(J)) times the perturbation amplitude "p".
! 
            CALL tl_initial (ng)
            IF (exit_flag.ne.NoError) THEN
              IF (Master) THEN
                WRITE (stdout,10) Rerror(exit_flag), exit_flag
              END IF
              RETURN
            END IF
            IF (Master) THEN
              lstr=LEN_TRIM(TLMbase(ng))
              WRITE (TLMname(ng),70) TLMbase(ng)(1:lstr-3), Nrun
            END IF
            IF (ndefTLM(ng).lt.0) THEN
              LdefTLM(ng)=.FALSE.              ! suppress IO
              LwrtTLM(ng)=.FALSE.
            ELSE
              LdefTLM(ng)=.TRUE.
            END IF
!
!  Time-step tangent linear model:  Compute misfit cost function
!  between model (nonlinear + tangent linear) and observations.
!
            IF (Master) THEN
              WRITE (stdout,20) 'TL', ntstart(ng), ntend(ng)
            END IF

            time(ng)=time(ng)-dt(ng)

            TL_LOOP : DO my_iic=ntstart(ng),ntend(ng)+1

              iic(ng)=my_iic
#ifdef SOLVE3D
              CALL tl_main3d (ng)
#else
              CALL tl_main2d (ng)
#endif
              IF (exit_flag.ne.NoError) THEN
                IF (Master) THEN
                  WRITE (stdout,10) Rerror(exit_flag), exit_flag
                END IF
                RETURN
              END IF

            END DO TL_LOOP
!
!  Close current tangent linear model history file.
!
            status=nf90_close(ncTLMid(ng))
            IF (status.ne.nf90_noerr) THEN
              WRITE (stdout,30) TRIM(TLMname(ng))
              exit_flag=3
              ioerror=status
              RETURN
            END IF
!
!  Advance model run counter.
!
            Nrun=Nrun+1

          END DO INNER_LOOP

        END DO OUTER_LOOP
!
!  Report dot products.
!
        IF (Master) THEN
          WRITE (stdout,80)                                             &
     &      'TLM Test - Dot Products Summary: p, g1, g2, (g1-g2)/g1'
          inner=1
          DO i=1,MIN(ig1count,ig2count)
            p=10.0_r8**FLOAT(-inner)
            IF (MOD(i,1+ntimes(ng)/nTLM(ng)).eq.0) inner=inner+1
            WRITE (stdout,90) i, p, g1(i), g2(i), (g1(i)-g2(i))/g1(i)
            IF ((MOD(i,1+ntimes(ng)/nTLM(ng)).eq.0).and.                &
     &          (MOD(i,Ninner*(1+ntimes(ng)/nTLM(ng))).ne.0)) THEN
              WRITE (stdout,100)
            ELSE IF (MOD(i,Ninner*(1+ntimes(ng)/nTLM(ng))).eq.0) THEN
              inner=1
              WRITE (stdout,110)
            END IF
          END DO
        END IF

      END DO NEST_LOOP
!
 10   FORMAT (/,a,i3,/)
 20   FORMAT (/,1x,a,1x,'ROMS/TOMS : started time-stepping:',           &
     &        '( TimeSteps: ',i8.8,' - ',i8.8,')',/)
 30   FORMAT (/,' TLCHECK_OCEAN - unable to close file : ',a)
 40   FORMAT (/,' Nonlinear Model Cost Function = ',1p,e21.14)
 50   FORMAT (' --------------- ','cost function = ',1p,e21.14,2x,a)
 60   FORMAT (17x,'cost function = ',1p,e21.14,2x,a)
 70   FORMAT (a,'_',i3.3,'.nc')
 80   FORMAT (/,a,/)
 90   FORMAT (i4,2x,1pe8.1,3(1x,1p,e20.12,0p))
100   FORMAT (77('.'))
110   FORMAT (77('-'))

      RETURN
      END SUBROUTINE run

      SUBROUTINE finalize
!
!=======================================================================
!                                                                      !
!  This routine terminates ROMS/TOMS nonlinear, tangent linear, and    !
!  adjoint models execution.                                           !
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
          CALL wclock_off (ng, iNLM, 0)
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
