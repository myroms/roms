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
!  ROMS/TOMS Weak Constraint Variational (4DVar) Data Assimilation     !
!      Driver: Physical-space Statistical Analysis System (PSAS)       !
!                                                                      !
!  This driver is used for weak constraint 4DVar where errors are      !
!  considered  in  both model and observations.  It is similar to      !
!  to the indirect representer driver but here the tangent linear      !
!  representer model is replaced with the nonlinear model.             !
!                                                                      !
!  The routines in this driver control the initialization,  time-      !
!  stepping, and finalization of  ROMS/TOMS  model following ESMF      !
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
      USE atm_coupler_mod, ONLY : initialize_coupling
#endif
#ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti
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

      integer :: STDrec, ng, thread

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
            CALL wclock_on (ng, iNLM, 0)
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
!
!  Read in model-error standard deviation factors and spatial
!  convolution diffusion convolution.
!  
        STDrec=1
        DO ng=1,Ngrids
          CALL get_state (ng, 6, 6, STDname(ng), STDrec, 1)
#ifdef DISTRIBUTE
          CALL mp_bcasti (ng, iNLM, exit_flag, 1)
#endif
          IF (exit_flag.ne.NoError) RETURN
        END DO

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
      USE mod_ocean
      USE mod_scalars
      USE mod_stepping
!
#ifdef BALANCE_OPERATOR
      USE ad_balance_mod, ONLY: ad_balance
#endif
      USE ad_convolution_mod, ONLY : ad_convolution
      USE ad_variability_mod, ONLY : ad_variability
#ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti
#endif
      USE impulse_mod, ONLY : impulse
      USE ini_adjust_mod, ONLY : ini_adjust
      USE ini_fields_mod, ONLY : ini_fields
      USE ini_adjust_mod, ONLY : load_ADtoTL
      USE ini_adjust_mod, ONLY : load_TLtoAD
      USE mod_ocean, ONLY : initialize_ocean
      USE normalization_mod, ONLY : normalization
#ifdef BALANCE_OPERATOR
      USE tl_balance_mod, ONLY: tl_balance
#endif
      USE tl_convolution_mod, ONLY : tl_convolution
      USE tl_variability_mod, ONLY : tl_variability
!
!  Local variable declarations.
!
      logical :: add, converged
      integer :: ADrec, Lbck, Lini
      integer :: i, lstr, my_iic, ng, rec, status, subs, tile, thread
      integer, save :: Nrec

      real(r8) :: LB_time, UB_time
!
!=======================================================================
!  Run model for all nested grids, if any.
!=======================================================================
!
      NEST_LOOP : DO ng=1,Ngrids
!
!  Initialize relevant parameters.
!
        Lold(ng)=1          ! old minimization time index
        Lnew(ng)=2          ! new minimization time index
        Lini=1              ! NLM initial conditions record in INIname
        Lbck=2              ! background record in INIname
        Nrun=1
        Ipass=1
        outer=0
        inner=0
        ERstr=1
        ERend=Nouter
!
!-----------------------------------------------------------------------
!  Configure weak constraint 4DVAR algorithm: Indirect Representer
!  Approach.
!-----------------------------------------------------------------------
!
!  Initialize and set nonlinear model initial conditions.
!
        wrtNLmod(ng)=.TRUE.
        wrtRPmod(ng)=.FALSE.
        wrtTLmod(ng)=.FALSE.
        CALL initial (ng)
        IF (exit_flag.ne.NoError) THEN
          IF (Master) THEN
            WRITE (stdout,10) Rerror(exit_flag), exit_flag
          END IF
          RETURN
        END IF
!
!  Save nonlinear initial conditions (currently in time index 1,
!  background) into record "Lbck" of INIname NetCDF file. The record
!  "Lbck" becomes the background state record and the record "Lini"
!  becomes current nonlinear initial conditions.
!
        IF (LcycleINI(ng)) THEN
          tINIindx(ng)=1
          NrecINI(ng)=1
        END IF
        CALL wrt_ini (ng, 1)
#ifdef DISTRIBUTE
        CALL mp_bcasti (ng, iNLM, exit_flag, 1)
#endif
!
!  Set nonlinear output history file as the initial basic state
!  trajectory.
!
        LdefHIS(ng)=.TRUE.
        LwrtHIS(ng)=.TRUE.
        lstr=LEN_TRIM(FWDbase(ng))
        WRITE (HISname(ng),20) FWDbase(ng)(1:lstr-3), 1
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Model-error covariance normalization and stardard deviation factors.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!  Compute or read in the model-error correlation normalization factors.
!  If computing, write out factors to NetCDF. This is an expensive
!  computation that needs to be computed only once for a particular
!  application grid and decorrelation scales.
!  
        IF (LwrtNRM(ng)) THEN
          CALL def_norm (ng)
!$OMP PARALLEL DO PRIVATE(ng,thread,subs,tile) SHARED(numthreads)
          DO thread=0,numthreads-1
            subs=NtileX(ng)*NtileE(ng)/numthreads
            DO tile=subs*thread,subs*(thread+1)-1
              CALL normalization (ng, TILE, 2)
            END DO
          END DO
!$OMP END PARALLEL DO
          LdefNRM(ng)=.FALSE.
          LwrtNRM(ng)=.FALSE.
        ELSE
          tNRMindx(ng)=1
          CALL get_state (ng, 5, 5, NRMname(ng), tNRMindx(ng), 1)
#ifdef DISTRIBUTE
          CALL mp_bcasti (ng, iNLM, exit_flag, 1)
#endif
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define TLM/RPM impulse forcing NetCDF file.
!
        LdefTLF(ng)=.TRUE.
        CALL def_impulse (ng)
!
!  Define output 4DVAR NetCDF file containing all processed data
!  at observation locations.
!
        LdefMOD(ng)=.TRUE.
        CALL def_mod (ng)
!
!  Define TLM initial conditions file.
!
        LdefITL(ng)=.FALSE.
        CALL tl_def_ini (ng)
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Run nonlinear model and compute basic state trajectory, X_n-1(t).
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
        IF (Master) THEN
          WRITE (stdout,30) 'NL', ntstart, ntend
        END IF

        SporadicImpulse=.FALSE.
        FrequentImpulse=.FALSE.
        time(ng)=time(ng)-dt(ng)

        NL_LOOP1 : DO my_iic=ntstart,ntend+1

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
        wrtNLmod(ng)=.FALSE.
!
!  Set model basic state trajectory file to nonlinear file.
!
        FWDname(ng)=HISname(ng)
        ncFWDid(ng)=ncHISid(ng)
!
!-----------------------------------------------------------------------
!  Solve the system:
!
!              (R_n + Cobs) * Beta_n = h_n
!
!              h_n = Xo - H * X_n
!
!  where R_n is the representer matrix, Cobs is the observation-error
!  covariance, Beta_n are the representer coefficients, h_n is the
!  misfit between observations (Xo) and model (H * X_n), and H is
!  the linearized observation operator. Here, _n denotes a sequence
!  of estimates.
!
!  The system does not need to be solved explicitly by inverting the
!  symmetric stabilized representer matrix, P_n:
!
!              P_n = R_n + Cobs
!
!  but by computing the action of P_n on any vector PSI, such that
!
!              P_n * PSI = R_n * PSI + Cobs * PSI
!
!  The representer matrix is not explicitly computed but evaluated by
!  one integration backward of the adjoint model and one integration
!  forward of the tangent linear model for any forcing vector PSI.
!
!  A preconditioned conjugate gradient algorithm is used to compute
!  an approximation PSI for Beta_n.
!
!-----------------------------------------------------------------------
!
        OUTER_LOOP : DO outer=1,Nouter
!
!  Set approximation vector PSI to representer coefficients Beta_n.
!  Here, PSI is set to misfit between observations and model, H_n.
!
          CALL congrad (ng, outer, 0, Ninner, converged)
!
!  Set basic state trajectory.
!
          lstr=LEN_TRIM(FWDbase(ng))
          WRITE (FWDname(ng),20) FWDbase(ng)(1:lstr-3), outer

          INNER_LOOP : DO inner=1,Ninner
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Integrate adjoint model forced with any vector PSI at the observation
!  locations and generate adjoint trajectory, Lambda_n(t).
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
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
            wrtMisfit(ng)=.FALSE.
!
!  Set adjoint history NetCDF parameters.  Define adjoint history
!  file one to avoid opening to many files.
!
            IF (Nrun.gt.1) LdefADJ(ng)=.FALSE.
            NrecADJ(ng)=0
            tADJindx(ng)=0
!
!  Time-step adjoint model backwards forced with current PSI vector.
!
            IF (Master) THEN
              WRITE (stdout,30) 'AD', ntstart, ntend
            END IF

            time(ng)=time(ng)+dt(ng)

            AD_LOOP1 : DO my_iic=ntstart,ntend,-1

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

            END DO AD_LOOP1

#ifdef CONVOLVE
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Convolve adjoint trajectory with model-error covariance and convert
!  to impulse forcing.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
            Nrec=NrecADJ(ng)
            NrecADJ(ng)=0
            tADJindx(ng)=0
            LwrtState2d(ng)=.TRUE.
            IF (Master) THEN
              WRITE (stdout,40) outer, inner
            END IF
!
!  Clear adjoint state arrays.
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile), SHARED(numthreads)
            DO thread=0,numthreads-1
              subs=NtileX(ng)*NtileE(ng)/numthreads
              DO tile=subs*thread,subs*(thread+1)-1
                CALL initialize_ocean (ng, TILE, iADM)
              END DO
            END DO
!$OMF END PARALLEL DO
!
!  Proccess each time record of current adjoint solution in ADJname.
!
            DO rec=1,Nrec
!
!  Read adjoint solution. Since routine "get_state" loads data into the
!  ghost points, the adjoint solution is read in the tangent linear
!  state arrays by using iTLM instead of iADM in the calling arguments.
!
              ADrec=rec
              CALL get_state (ng, iTLM, 4, ADJname(ng), ADrec, Lold(ng))
!
!  Load interior solution, read above, into adjoint state arrays. 
!  Then, multiply adjoint solution by the background-error standard
!  deviations. Next, convolve resulting adjoint solution with the 
!  squared-root adjoint diffusion operator which impose the model-error
!  spatial correlations. Notice that the spatial convolution is only
!  done for half of the diffusion steps (squared-root filter). Clear
!  tangent linear state arrays when done.
!
!  WARNING: in weak constraint we need to use  different statistics
!           and correlation normalization coefficients for initial
!           impulse (last adjoint record, rec=Nrec) versus other
!           impulses (adjoint records 1:Nrec-1).
!
              add=.FALSE.
!$OMP PARALLEL DO PRIVATE(ng,thread,subs,tile)                          &
!$OMP&            SHARED(inner,add,numthreads)
              DO thread=0,numthreads-1
                subs=NtileX(ng)*NtileE(ng)/numthreads
                DO tile=subs*thread,subs*(thread+1)-1,+1
                  CALL load_TLtoAD (ng, TILE, Lold(ng), Lold(ng), add)
#ifdef BALANCE_OPERATOR
                  CALL ad_balance (ng, TILE, Lold(ng))
#endif
                  CALL ad_variability (ng, TILE, Lold(ng))
                  CALL ad_convolution (ng, TILE, Lold(ng), 2)
                  CALL initialize_ocean (ng, TILE, iTLM)
                END DO
              END DO
!$OMP END PARALLEL DO
!
!  To insure symmetry, convolve resulting filtered adjoint solution
!  from above with the squared-root (half of steps) tangent linear
!  diffusion operator. Then, multiply result with its corresponding
!  background-error standard deviations.  Since the convolved solution
!  is in the adjoint state arrays, first copy to tangent linear state
!  arrays including the ghosts points. Copy back to adjoint state
!  arrays when done with the convolution for output purposes.
!
              add=.FALSE.
!$OMP PARALLEL DO PRIVATE(ng,thread,subs,tile)                          &
!$OMP&            SHARED(inner,add,numthreads)
              DO thread=0,numthreads-1
                subs=NtileX(ng)*NtileE(ng)/numthreads
                DO tile=subs*thread,subs*(thread+1)-1
                  CALL load_ADtoTL (ng, TILE, Lold(ng), Lold(ng), add)
                  CALL tl_convolution (ng, TILE, Lold(ng), 2)
                  CALL tl_variability (ng, TILE, Lold(ng))
#ifdef BALANCE_OPERATOR
                  CALL tl_balance (ng, TILE, Lold(ng))
#endif
                  CALL load_TLtoAD (ng, TILE, Lold(ng), Lold(ng), add)
                END DO
              END DO
!$OMP END PARALLEL DO
!
!  Overwrite ADJname history NetCDF file with convolved adjoint
!  solution.
!
              kstp(ng)=Lold(ng)
# ifdef SOLVE3D
              nstp(ng)=Lold(ng)
# endif
              CALL ad_wrt_his (ng)
              IF (exit_flag.ne.NoError) RETURN
            END DO
            LwrtState2d(ng)=.FALSE.
#endif
!
!  Convert the current adjoint solution in ADJname to impulse forcing.
!  Write out impulse forcing into TLFname NetCDF file. To facilitate
!  the forcing to the TLM and RPM, the forcing is processed and written
!  in increasing time coordinates (recall that the adjoint solution
!  in ADJname is backwards in time).
!
            IF (Master) THEN
              WRITE (stdout,50) outer, inner
            END IF
            tTLFindx(ng)=0
            CALL impulse (ng, iADM, ADJname(ng))
            NrecFrc(ng)=Nrec
!
!  Write new tangent linear model initial conditions into ITLname
!  file.  Use initial impulse, FrcRec=1.
!
            FrcRec(ng)=1
            CALL get_state (ng, iTLM, 7, TLFname(ng), FrcRec(ng),       &
     &                      Lold(ng))
            CALL tl_wrt_ini (ng, Lold(ng), 1)
#ifdef DISTRIBUTE
            CALL mp_bcasti (ng, iTLM, exit_flag, 1)
#endif
            IF (exit_flag.ne.NoError) RETURN
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Integrate tangent linear model forced by the convolved adjoint
!  trajectory (impulse forcing) to compute R_n * PSI at observation
!  points.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!  Initialize tangent linear model from initial impulse which is now
!  stored in file ITLname.
! 
            wrtNLmod(ng)=.FALSE.
            wrtTLmod(ng)=.TRUE.
            CALL tl_initial (ng)
            IF (exit_flag.ne.NoError) THEN
              IF (Master) THEN
                WRITE (stdout,10) Rerror(exit_flag), exit_flag
              END IF
              RETURN
            END IF
!
!  Activate switch to write out initial misfit between model and
!  observations.
!
            IF ((outer.eq.1).and.(inner.eq.1)) THEN
              wrtMisfit(ng)=.TRUE.
            END IF
!
!  Set tangent linear history NetCDF parameters.  Define tangent linear
!  history file at the beggining of each inner loop  to avoid opening
!  too many NetCDF files.
!
            IF (inner.gt.1) LdefTLM(ng)=.FALSE.
            NrecTLM(ng)=0
            tTLMindx(ng)=0
!
!  Run tangent linear model forward and force with convolved adjoint
!  trajectory impulses. Compute R_n * PSI at observation points which
!  are used in the conjugate gradient algorithm.
!
            IF (Master) THEN
              WRITE (stdout,30) 'TL', ntstart, ntend
            END IF

            time(ng)=time(ng)-dt(ng)

            TL_LOOP : DO my_iic=ntstart,ntend+1

              iic(ng)=my_iic
!
!  Set impulse forcing switches.
!
              LB_time=time(ng)+0.5_r8*dt(ng)
              UB_time=time(ng)+1.5_r8*dt(ng)
              SporadicImpulse=(LB_time.le.FrcTime(ng)).and.             &
     &                        (FrcTime(ng).lt.UB_time).and.             &
     &                        (FrcRec(ng).gt.1)
              FrequentImpulse=.FALSE.
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
            wrtNLmod(ng)=.FALSE.
            wrtTLmod(ng)=.FALSE.
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Use conjugate gradient algorithm to find a better approximation
!  PSI to representer coefficients Beta_n. Exit inner loop if
!  convergence is achieved.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
            Nrun=Nrun+1
            CALL congrad (ng, outer, inner, Ninner, converged)
            IF (converged) EXIT INNER_LOOP

          END DO INNER_LOOP
!
!  Close tangent linear NetCDF file.
!
          status=nf_close(ncTLMid(ng))
          ncTLMid(ng)=-1
!
!-----------------------------------------------------------------------
!  Once that the representer coefficients, Beta_n, have been
!  approximated with sufficient accuracy, compute estimates of
!  Lambda_n and Xhat_n by carrying out one backward intergration
!  of the adjoint model and one forward itegration of the representer
!  model.
!-----------------------------------------------------------------------
!
!  Initialize the adjoint model always from rest.
!
          CALL ad_initial (ng)
          IF (exit_flag.ne.NoError) THEN
            IF (Master) THEN
              WRITE (stdout,10) Rerror(exit_flag), exit_flag
            END IF
            RETURN
          END IF
!
!  Set adjoint history NetCDF parameters.  Define adjoint history
!  file one to avoid opening to many files.
!
          IF (Nrun.gt.1) LdefADJ(ng)=.FALSE.
          NrecADJ(ng)=0
          tADJindx(ng)=0
!
!  Time-step adjoint model backwards forced with estimated representer
!  coefficients, Beta_n.
!
          IF (Master) THEN
            WRITE (stdout,30) 'AD', ntstart, ntend
          END IF

          time(ng)=time(ng)+dt(ng)

          AD_LOOP2 : DO my_iic=ntstart,ntend,-1

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

          END DO AD_LOOP2

#ifdef CONVOLVE
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Convolve adjoint trajectory with model-error covariance and convert
!  to impulse forcing.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
          Nrec=NrecADJ(ng)
          NrecADJ(ng)=0
          tADJindx(ng)=0
          LwrtState2d(ng)=.TRUE.
          IF (Master) THEN
            WRITE (stdout,40) outer, 0
          END IF
!
!  Clear tangent linear and adjoint state arrays.
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile), SHARED(numthreads)
          DO thread=0,numthreads-1
            subs=NtileX(ng)*NtileE(ng)/numthreads
            DO tile=subs*thread,subs*(thread+1)-1
              CALL initialize_ocean (ng, TILE, iADM)
            END DO
          END DO
!$OMF END PARALLEL DO
!
!  Proccess each time record of current adjoint solution in ADJname.
!
          DO rec=1,Nrec
!
!  Read adjoint solution. Since routine "get_state" loads data into the
!  ghost points, the adjoint solution is read in the tangent linear
!  state arrays by using iTLM instead of iADM in the calling arguments.
!
            ADrec=rec
            CALL get_state (ng, iTLM, 4, ADJname(ng), ADrec, Lold(ng))
!
!  Load interior solution, read above, into adjoint state arrays. 
!  Then, multiply adjoint solution by the background-error standard
!  deviations. Next, convolve resulting adjoint solution with the 
!  squared-root adjoint diffusion operator which impose the model-error
!  spatial correlations. Notice that the spatial convolution is only
!  done for half of the diffusion steps (squared-root filter). Clear
!  tangent linear state arrays when done.
!
!  WARNING: in weak constraint we need to use  different statistics
!           and correlation normalization coefficients for initial
!           impulse (last adjoint record, rec=Nrec) versus other
!           impulses (adjoint records 1:Nrec-1).
!
            add=.FALSE.
!$OMP PARALLEL DO PRIVATE(ng,thread,subs,tile)                          &
!$OMP&            SHARED(inner,add,numthreads)
            DO thread=0,numthreads-1
              subs=NtileX(ng)*NtileE(ng)/numthreads
              DO tile=subs*thread,subs*(thread+1)-1
                CALL load_TLtoAD (ng, TILE, Lold(ng), Lold(ng), add)
#ifdef BALANCE_OPERATOR
                CALL ad_balance (ng, TILE, Lold(ng))
#endif
                CALL ad_variability (ng, TILE, Lold(ng))
                CALL ad_convolution (ng, TILE, Lold(ng), 2)
                CALL initialize_ocean (ng, TILE, iTLM)
              END DO
            END DO
!$OMP END PARALLEL DO
!
!  To insure symmetry, convolve resulting filtered adjoint solution
!  from above with the squared-root (half of steps) tangent linear
!  diffusion operator. Then, multiply result with its corresponding
!  background-error standard deviations.  Since the convolved solution
!  is in the adjoint state arrays, first copy to tangent linear state
!  arrays including the ghosts points. Copy back to adjoint state
!  arrays when done with the convolution for outpur purposes.
!
            add=.FALSE.
!$OMP PARALLEL DO PRIVATE(ng,thread,subs,tile)                          &
!$OMP&            SHARED(inner,add,numthreads)
            DO thread=0,numthreads-1
              subs=NtileX(ng)*NtileE(ng)/numthreads
              DO tile=subs*thread,subs*(thread+1)-1,+1
                CALL load_ADtoTL (ng, TILE, Lold(ng), Lold(ng), add)
                CALL tl_convolution (ng, TILE, Lold(ng), 2)
                CALL tl_variability (ng, TILE, Lold(ng))
#ifdef BALANCE_OPERATOR
                CALL tl_balance (ng, TILE, Lold(ng))
#endif
                CALL load_TLtoAD (ng, TILE, Lold(ng), Lold(ng), add)
              END DO
            END DO
!$OMP END PARALLEL DO
!
!  Overwrite ADJname history NetCDF file with convolved adjoint
!  solution.
!
            kstp(ng)=Lold(ng)
# ifdef SOLVE3D
            nstp(ng)=Lold(ng)
# endif
            CALL ad_wrt_his (ng)
            IF (exit_flag.ne.NoError) RETURN
          END DO
          LwrtState2d(ng)=.FALSE.
#endif
!
!  Convert convolved adjoint solution to impulse forcing. Write out
!  impulse forcing into TLFname NetCDF file. To facilitate the forcing
!  by the TLM and RPM, the forcing is process and written in
!  increasing time coordinates.
!
          IF (Master) THEN
            WRITE (stdout,50) outer, 0
          END IF
          tTLFindx(ng)=0
          CALL impulse (ng, iADM, ADJname(ng))
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Run nonlinear model and compute a "new estimate" of the state
!  trajectory, X_n(t).
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!  Set new basic state trajectory.
!
          LdefHIS(ng)=.TRUE.
          LwrtHIS(ng)=.TRUE.
          lstr=LEN_TRIM(FWDbase(ng))
          WRITE (HISname(ng),20) FWDbase(ng)(1:lstr-3), outer+1
!
!  Initialize nonlinear model with background or reference state. The
!  initial contribution from the adjoint model, Beta_n(0), will be added
!  as impulse forcing in "nl_forcing".
!
          wrtNLmod(ng)=.TRUE.
          wrtTLmod(ng)=.FALSE.
          CALL initial (ng)
          IF (exit_flag.ne.NoError) THEN
            IF (Master) THEN
              WRITE (stdout,10) Rerror(exit_flag), exit_flag
            END IF
            RETURN
          END IF
!
!  Activate switch to write out final misfit between model and
!  observations.
!
          IF (outer.eq.Nouter) THEN
            wrtMisfit(ng)=.TRUE.
          END IF
!
!  Run nonlinear forced by convolved adjoint trajectory impulses and
!  compute new basic state trajectory X_n.
!
          IF (Master) THEN
            WRITE (stdout,30) 'NL', ntstart, ntend
          END IF

          time(ng)=time(ng)-dt(ng)

          NL_LOOP2 : DO my_iic=ntstart,ntend+1

            iic(ng)=my_iic
!
!  Set impulse forcing switches.
!
            SporadicImpulse=.FALSE.
            FrequentImpulse=FrcRec(ng).gt.1

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
          wrtNLmod(ng)=.FALSE.
          wrtTLmod(ng)=.FALSE.
!
!  Set basic state trajectory file for next outer loop iteration.  Close
!  current forward NetCDF file.
!
          FWDname(ng)=HISbase(ng)
          status=nf_close(ncFWDid(ng))
          ncFWDid(ng)=-1

        END DO OUTER_LOOP
!
!-----------------------------------------------------------------------
!  Compute new nonlinear model initial conditions by adding last
!  convolved adjoint solution (Beta(t=0), currently saved in
!  record Nrec of ADJname) to the background state.
!-----------------------------------------------------------------------
!
!  Read in adjoint (t=0) and background state.
!
        ADrec=Nrec
        CALL get_state (ng, iADM, 4, ADJname(ng), ADrec, Lold(ng))
        CALL get_state (ng, iNLM, 1, INIname(ng), Lbck, Lini)
!
!  Notice that "ini_fields" is called below for output purposes only. 
!  It computes the vertically integrated momentum in 3D applications.
!  In order to use the correct fields, the model time indices are set
!  to Lini.
!
          kstp(ng)=Lini
#ifdef SOLVE3D
          nstp(ng)=Lini
#endif
!
!  Compute nonlinear model initial conditions.
!
!$OMP PARALLEL DO PRIVATE(ng,thread,subs,tile)                          &
!$OMP&            SHARED(numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1
            CALL ini_adjust (ng, TILE, Lold(ng), Lini)
            CALL ini_fields (ng, TILE, iNLM)
          END DO
        END DO
!$OMP END PARALLEL DO
!
!  Write out new nonlinear model initial conditions.
!
        IF (LcycleINI(ng)) THEN
          tINIindx(ng)=0
          NrecINI(ng)=1
        END IF
        CALL wrt_ini (ng, Lini)
#ifdef DISTRIBUTE
        CALL mp_bcasti (ng, iNLM, exit_flag, 1)
#endif
        IF (exit_flag.ne.NoError) RETURN
!
!  Compute and report model-observation comparison statistics.
!
        CALL stats_modobs (ng)

      END DO NEST_LOOP
!
 10   FORMAT (/,a,i3,/)
 20   FORMAT (a,'_',i3.3,'.nc')
 30   FORMAT (/,1x,a,1x,'ROMS/TOMS: started time-stepping:',            &
     &        '( TimeSteps: ',i8.8,' - ',i8.8,')',/)
 40   FORMAT (/,' Convolving Adjoint Tracjectory: Outer = ',i3.3,       &
     &          ' Inner = ',i3.3)
 50   FORMAT (/,' Converting Convolved Adjoint Tracjectory to',         &
     &          ' Impulses: Outer = ',i3.3,' Inner = ',i3.3,/)

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
!  If cycling restart records, write solution into record 3.
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
 20     FORMAT (/,'Elapsed CPU time (seconds):',/)
      END IF

      DO ng=1,Ngrids
!$OMP PARALLEL DO PRIVATE(ng,thread) SHARED(numthreads)
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
