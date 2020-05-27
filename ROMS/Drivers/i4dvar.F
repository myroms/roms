#include "cppdefs.h"
      MODULE i4dvar_mod

#ifdef I4DVAR
!
!git $Id$
!svn $Id: i4dvar.F 1027 2020-05-27 17:06:32Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2020 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module splits the I4D-Var data assimilation algorithm into     !
!  its logical components routines:                                    !
!                                                                      !
!    background          Xb, nonlinear trajectory used to linearize    !
!                            the tangent linear and adjoint models     !
!                                                                      !
!    increment           dXa, 4D-Var data assimilation increment:      !
!                             inner loops and minimization solver      !
!                                                                      !
!    analysis            Xa, analyzed state: Xa = Xb + dXa             !
!                                                                      !
!    posterior_analysis  Initialize NLM with estimated 4D-Var state    !
!                        and interpolates solution at observation      !
!                        locations                                     !
!                                                                      !
!    prior_error         Process prior background error covariance     !
!                        and its normalization coefficients            !
!                                                                      !
!  References:                                                         !
!                                                                      !
!    Moore, A.M., H.G. Arango, G. Broquet, B.S. Powell, A.T. Weaver,   !
!      and J. Zavala-Garay, 2011: The Regional Ocean Modeling System   !
!      (ROMS)  4-dimensional variational data assimilations systems,   !
!      Part I - System overview and formulation, Prog. Oceanogr., 91,  !
!      34-49, doi:10.1016/j.pocean.2011.05.004.                        !
!                                                                      !
!    Moore, A.M., H.G. Arango, G. Broquet, C. Edward, M. Veneziani,    !
!      B. Powell, D. Foley, J.D. Doyle, D. Costa, and P. Robinson,     !
!      2011: The Regional Ocean Modeling System (ROMS) 4-dimensional   !
!      variational data assimilations systems, Part II - Performance   !
!      and application to the California Current System, Prog.         !
!      Oceanogr., 91, 50-73, doi:10.1016/j.pocean.2011.05.003.         !
!                                                                      !
!=======================================================================
!
      USE mod_kinds
!
      implicit none
!
      PUBLIC  :: background
      PUBLIC  :: increment
      PUBLIC  :: analysis
      PUBLIC  :: posterior_analysis
      PUBLIC  :: prior_error
!
!  Set module internal parameters.
!
      integer :: LTLM1 = 1      ! trial x-space TLM IC record in ITL
      integer :: LTLM2 = 2      ! previous v-space TLM IC record in ITL
      integer :: LTLM3 = 3      ! trial v-space TLM IC record in ITL
      integer :: LADJ1 = 1      ! initial cost gradient
      integer :: LADJ2 = 2      ! new cost gradient (not normalized)
      integer :: Lini  = 1      ! NLM initial conditions record in INI
      integer :: Lbck  = 2      ! background record in INI
      integer :: Rec1  = 1
      integer :: Rec2  = 2
      integer :: Rec3  = 3
      integer :: Rec4  = 4
!
      CONTAINS
!
      SUBROUTINE background (my_outer, RunInterval)
!
!=======================================================================
!                                                                      !
!  This routine computes the backgound state trajectory, Xb_n-1(t),    !
!  used to linearize the tangent linear and adjoint models in the      !
!  inner loops. It interpolates the background at the observations     !
!  locations, and computes the accept/reject quality control flag,     !
!  ObsScale.                                                           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     my_outer        Outer-loop counter (integer)                     !
!     RunInterval     NLM kernel time stepping window (seconds)        !
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
      USE mod_mixing,       ONLY : initialize_mixing
      USE strings_mod,      ONLY : FoundError
# if defined ADJUST_BOUNDARY || defined ADJUST_STFLUX || \
     defined ADJUST_WSTRESS
      USE wrt_ini_mod,      ONLY : wrt_frc
# endif
      USE wrt_ini_mod,      ONLY : wrt_ini
# if defined BALANCE_OPERATOR && defined ZETA_ELLIPTIC
      USE zeta_balance_mod, ONLY : balance_ref, biconj
# endif
!
!  Imported variable declarations
!
      integer, intent(in)  :: my_outer
!
      real(dp), intent(in) :: RunInterval
!
!  Local variable declarations.
!
      integer :: i, lstr, ng, tile
      integer :: Fcount
# ifdef PROFILE
      integer :: thread
# endif
!
      SourceFile=__FILE__
!
!-----------------------------------------------------------------------
!  Initialize and set nonlinear model initial conditions.
!-----------------------------------------------------------------------

# ifdef PROFILE
!
!  Start profile clock.
!
      DO ng=1,Ngrids
        DO thread=THREAD_RANGE
          CALL wclock_on (ng, iNLM, 86, __LINE__,                       &
     &                    __FILE__)
        END DO
      END DO
# endif
!
!  Set nonlinear output history file name. Create a basic state file
!  for each outher loop.
!
      DO ng=1,Ngrids
        idefHIS(ng)=-1
        LdefHIS(ng)=.TRUE.
        LwrtHIS(ng)=.TRUE.
# if defined BULK_FLUXES && defined NL_BULK_FLUXES
        LreadBLK(ng)=.FALSE.
# endif
        LreadFWD(ng)=.FALSE.
        WRITE (HIS(ng)%name,10) TRIM(FWD(ng)%head), my_outer-1
        lstr=LEN_TRIM(HIS(ng)%name)
        HIS(ng)%base=HIS(ng)%name(1:lstr-3)
      END DO
!
!  Clear nonlinear mixing arrays.
!
      DO ng=1,Ngrids
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL initialize_mixing (ng, tile, iNLM)
        END DO
      END DO
!
!  Initialize nonlinear model. If outer=1, the model is initialized
!  with the background or reference state. Otherwise, the model is
!  initialized with the estimated initial conditions from previous
!  iteration, X(0) = X(0) + deltaX(0).
!
      DO ng=1,Ngrids
        wrtNLmod(ng)=.TRUE.
        wrtTLmod(ng)=.FALSE.
        RST(ng)%Rindex=0
        Fcount=RST(ng)%load
        RST(ng)%Nrec(Fcount)=0
      END DO
!
      CALL initial
      IF (FoundError(exit_flag, NoError, __LINE__,                      &
     &               __FILE__)) RETURN
!
!  If first pass, save nonlinear initial conditions (currently in time
!  index 1, background) into next record (Lbck) of INI(ng)%name NetCDF
!  file. The record "Lbck" becomes the background state record and the
!  record "Lini" becomes current nonlinear initial conditions.  Both
!  records are used in the algorithm below.
!
      IF (my_outer.eq.1) THEN
        DO ng=1,Ngrids
          INI(ng)%Rindex=1
          Fcount=INI(ng)%load
          INI(ng)%Nrec(Fcount)=1
          CALL wrt_ini (ng, 1)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END DO
      END IF
!
!  Create 4D-Var analysis file that used as initial conditions for the
!  next data assimilation cycle now that grid information is known.
!
      IF (my_outer.eq.1) THEN
        DO ng=1,Ngrids
          LdefDAI(ng)=.TRUE.
          CALL def_dai (ng)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END DO
      END IF

# if defined BALANCE_OPERATOR && defined ZETA_ELLIPTIC
!
!  Compute the reference zeta and biconjugate gradient arrays
!  required for the balance of free surface.
!
      IF (balance(isFsur)) THEN
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL balance_ref (ng, tile, Lini)
            CALL biconj (ng, tile, iNLM, Lini)
          END DO
          wrtZetaRef(ng)=.TRUE.
        END DO
      END IF
# endif
!
!  If first pass, define output 4DVAR NetCDF file containing all
!  processed data at observation locations.
!
      IF (my_outer.eq.1) THEN
        DO ng=1,Ngrids
          LdefMOD(ng)=.TRUE.
          CALL def_mod (ng)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END DO
      END IF
!
!  Run nonlinear model. Save nonlinear tracjectory needed by the
!  adjoint and tangent linear models. Interpolate nonlinear model
!  to observation locations (compute and save H x). It processes
!  and writes the observations accept/reject flag (ObsScale) once
!  to allow background quality control, if any.
!
      DO ng=1,Ngrids
# ifdef AVERAGES
        idefAVG(ng)=-1
        LdefAVG(ng)=.TRUE.
        LwrtAVG(ng)=.TRUE.
        WRITE (AVG(ng)%name,10) TRIM(AVG(ng)%head), my_outer
        lstr=LEN_TRIM(AVG(ng)%name)
        AVG(ng)%base=AVG(ng)%name(1:lstr-3)
# endif
# ifdef DIAGNOSTICS
        idefDIA(ng)=-1
        LdefDIA(ng)=.TRUE.
        LwrtDIA(ng)=.TRUE.
        WRITE (DIA(ng)%name,10) TRIM(DIA(ng)%head), my_outer
        lstr=LEN_TRIM(DIA(ng)%name)
        DIA(ng)%base=DIA(ng)%name(1:lstr-3)
# endif
        wrtObsScale(ng)=.TRUE.
        IF (Master) THEN
          WRITE (stdout,20) 'NL', ng, ntstart(ng), ntend(ng)
        END IF
      END DO
!
# ifdef SOLVE3D
      CALL main3d (RunInterval)
# else
      CALL main2d (RunInterval)
# endif
      IF (FoundError(exit_flag, NoError, __LINE__,                      &
     &               __FILE__)) RETURN
!
      DO ng=1,Ngrids
# ifdef AVERAGES
        LdefAVG(ng)=.FALSE.
        LwrtAVG(ng)=.FALSE.
# endif
# ifdef DIAGNOSTICS
        LdefDIA(ng)=.FALSE.
        LwrtDIA(ng)=.FALSE.
# endif
        wrtNLmod(ng)=.FALSE.
        wrtObsScale(ng)=.FALSE.
        wrtTLmod(ng)=.TRUE.
      END DO
!
!  Set structure for the nonlinear forward trajectory to be processed
!  by the tangent linear and adjoint models. Also, set switches to
!  process the FWD structure in routine "check_multifile". Notice that
!  it is possible to split solution into multiple NetCDF files to reduce
!  their size.
!
      CALL edit_multifile ('HIS2FWD')
      IF (FoundError(exit_flag, NoError, __LINE__,                      &
     &               __FILE__)) RETURN
      DO ng=1,Ngrids
        LreadFWD(ng)=.TRUE.
      END DO

# if defined BULK_FLUXES && defined NL_BULK_FLUXES
!
!  Set structure for the nonlinear surface fluxes to be processed by
!  by the tangent linear and adjoint models. Also, set switches to
!  process the BLK structure in routine "check_multifile".  Notice that
!  it is possible to split solution into multiple NetCDF files to reduce
!  their size.
!
      CALL edit_multifile ('HIS2BLK')
      IF (FoundError(exit_flag, NoError, __LINE__,                      &
     &               __FILE__)) RETURN
      DO ng=1,Ngrids
        LreadBLK(ng)=.TRUE.
      END DO
# endif

# if defined ADJUST_BOUNDARY || defined ADJUST_STFLUX || \
     defined ADJUST_WSTRESS
!
!  Write out initial and background surface forcing into initial
!  INI(ng)%name NetCDF file for latter use.
!
      DO ng=1,Ngrids
        CALL wrt_frc (ng, Lfout(ng), Lini)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN

        IF (my_outer.eq.1) THEN
          CALL wrt_frc (ng, Lfout(ng), Lbck)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
      END DO
# endif
!
!  Write out nonlinear model misfit cost function into DAV(ng)%name
!  NetCDF file.
!
      SourceFile=__FILE__ // ", background"
      DO ng=1,Ngrids
        CALL netcdf_put_fvar (ng, iNLM, DAV(ng)%name,                   &
     &                        'NLcost_function',                        &
     &                        FOURDVAR(ng)%NLobsCost(0:),               &
     &                        (/1,my_outer/), (/NobsVar(ng)+1,1/),      &
     &                        ncid = DAV(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN
      END DO

# ifdef PROFILE
!
!  Stop profile clock
!
      DO ng=1,Ngrids
        DO thread=THREAD_RANGE
          CALL wclock_off (ng, iNLM, 86, __LINE__,                      &
     &                     __FILE__)
        END DO
      END DO
# endif
!
 10   FORMAT (a,'_outer',i0,'.nc')
 20   FORMAT (/,1x,a,1x,'ROMS/TOMS: started time-stepping:',            &
     &        ' (Grid: ',i2.2,' TimeSteps: ',i8.8,' - ',i8.8,')',/)
!
      RETURN
      END SUBROUTINE background
!
      SUBROUTINE increment (my_outer, RunInterval)
!
!=======================================================================
!                                                                      !
!  This routine computes the 4D-Var data assimilation state increment, !
!  dXa, by iterating the inner loops and minimizing the cost function. !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     my_outer        Outer-loop counter (integer)                     !
!     RunInterval     TLM/ADM kernels time stepping window (seconds)   !
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
# ifdef BALANCE_OPERATOR
      USE ad_balance_mod,     ONLY : ad_balance
# endif
      USE ad_convolution_mod, ONLY : ad_convolution
      USE ad_variability_mod, ONLY : ad_variability
      USE back_cost_mod,      ONLY : back_cost
      USE cgradient_mod,      ONLY : cgradient
      USE cost_grad_mod,      ONLY : cost_grad
# ifdef ADJUST_BOUNDARY
      USE mod_boundary,       ONLY : initialize_boundary
# endif
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      USE mod_forces,         ONLY : initialize_forces
# endif
      USE mod_ocean,          ONLY : initialize_ocean
      USE strings_mod,        ONLY : FoundError
# ifdef BALANCE_OPERATOR
      USE tl_balance_mod,     ONLY : tl_balance
# endif
      USE tl_convolution_mod, ONLY : tl_convolution
      USE tl_variability_mod, ONLY : tl_variability
!
!  Imported variable declarations
!
      logical :: Lweak = .FALSE.
!
      integer, intent(in) :: my_outer
!
      real(dp), intent(in) :: RunInterval
!
!  Local variable declarations.
!
      integer :: i, lstr, my_inner, ng, tile
      integer :: Fcount, Lcon, Lsav
# ifdef PROFILE
      integer :: thread
# endif
!
      real(r8) :: rate
!
      SourceFile=__FILE__
!
!=======================================================================
!  Compute 4D-Var increment.
!=======================================================================

# ifdef PROFILE
!
!  Start profile clock.
!
      DO ng=1,Ngrids
        DO thread=THREAD_RANGE
          CALL wclock_on (ng, iTLM, 87, __LINE__,                       &
     &                    __FILE__)
        END DO
      END DO
# endif
!
!  The minimization algorithm requires to save all the gradient
!  solutions for each inner loop iteration.  They are used for
!  orthogonalization in the conjugate gradient algorithm.  Thus,
!  we need to reset adjoint file record indices.
!
      DO ng=1,Ngrids
        ADM(ng)%Rindex=0
        Fcount=ADM(ng)%load
        ADM(ng)%Nrec(Fcount)=0
      END DO
!
!  An adjoint NetCDF is created for each outer loop.
!
      DO ng=1,Ngrids
        idefADJ(ng)=-1
        LdefADJ(ng)=.TRUE.
        WRITE (ADM(ng)%name,10) TRIM(ADM(ng)%head), my_outer
        lstr=LEN_TRIM(ADM(ng)%name)
        ADM(ng)%base=ADM(ng)%name(1:lstr-3)
      END DO
!
!  Define output Hessian NetCDF file containing the eigenvectors
!  approximation to the Hessian matrix computed from the Lanczos
!  algorithm. Notice that the file name is a function of the
!  outer loop. That is, a file is created for each outer loop.
!
      DO ng=1,Ngrids
        LdefHSS(ng)=.TRUE.
        WRITE (HSS(ng)%name,10) TRIM(HSS(ng)%head), my_outer
        lstr=LEN_TRIM(HSS(ng)%name)
        HSS(ng)%base=HSS(ng)%name(1:lstr-3)
        CALL def_hessian (ng)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN
      END DO
!
!  Notice that inner loop iteration start from zero. This is needed to
!  compute the minimization initial increment deltaX(0), its associated
!  gradient G(0), and descent direction d(0) used in the conjugate
!  gradient algorithm.
!
      INNER_LOOP : DO my_inner=0,Ninner
        inner=my_inner
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Time-step tangent linear model: compute cost function.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!  If first pass inner=0, initialize tangent linear state (increments,
!  deltaX) from rest. Otherwise, use trial initial conditions estimated
!  by the conjugate gradient algorithm in previous inner loop. The TLM
!  initial conditions are read from ITL(ng)%name, record 1.
!
        DO ng=1,Ngrids
          ITL(ng)%Rindex=1
          CALL tl_initial (ng)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END DO
!
!  On first pass, initialize records 2, 3 and 4 of the ITL file to zero.
!
        IF ((my_inner.eq.0).and.(my_outer.eq.1)) THEN
          DO ng=1,Ngrids
            CALL tl_wrt_ini (ng, LTLM1, Rec2)
            IF (FoundError(exit_flag, NoError, __LINE__,                &
     &                     __FILE__)) RETURN
            CALL tl_wrt_ini (ng, LTLM1, Rec3)
            IF (FoundError(exit_flag, NoError, __LINE__,                &
     &                     __FILE__)) RETURN
            CALL tl_wrt_ini (ng, LTLM1, Rec4)
            IF (FoundError(exit_flag, NoError, __LINE__,                &
     &                     __FILE__)) RETURN
          END DO
        END IF

# ifdef MULTIPLE_TLM
!
!  If multiple TLM history NetCDF files, activate writing and determine
!  output filename. The multiple file option is used to perturb initial
!  state and create ensembles.  The TLM final trajectory is written for
!  each inner loop on separated NetCDF files.
!
        DO ng=1,Ngrids
          idefTLM(ng)=-1
          LdefTLM(ng)=.TRUE.
          LwrtTLM(ng)=.TRUE.
          WRITE (TLM(ng)%name,20) TRIM(TLM(ng)%head), Nrun
          lstr=LEN_TRIM(TLM(ng)%name)
          TLM(ng)%base=TLM(ng)%name(1:lstr-3)
        END DO
# endif
!
!  Activate switch to write out initial and final misfit between
!  model and observations.
!
        DO ng=1,Ngrids
          wrtMisfit(ng)=.FALSE.
          IF (((my_outer.eq.1).and.(my_inner.eq.0)).or.                 &
     &        ((my_outer.eq.Nouter).and.(my_inner.eq.Ninner))) THEN
            wrtMisfit(ng)=.TRUE.
          END IF
        END DO
!
!  Run tangent linear model. Compute misfit observation cost function,
!  Jo.
!
        DO ng=1,Ngrids
          IF (Master) THEN
            WRITE (stdout,30) 'TL', ng, ntstart(ng), ntend(ng)
          END IF
        END DO
!
# ifdef SOLVE3D
        CALL tl_main3d (RunInterval)
# else
        CALL tl_main2d (RunInterval)
# endif
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN

# ifdef EVOLVED_LCZ
!
!  Write evolved tangent Lanczos vector into hessian netcdf file for use
!  later.
!
!  NOTE: When using this option, it is important to set LhessianEV and
!  Lprecond to FALSE in s4dvar.in, otherwise the evolved Lanczos vectors
!  with be overwritten by the Hessian eigenvectors. The fix to this
!  is to define a new NetCDF file that contains the evolved Lanczos
!  vectors.
!
        IF (my_inner.ne.0) THEN
          DO ng=1,Ngrids
            CALL wrt_evolved (ng, kstp(ng), nrhs(ng))
            IF (FoundError(exit_flag, NoError, __LINE__,                &
     &                     __FILE__)) RETURN
          END DO
        END IF
# endif

# ifdef MULTIPLE_TLM
!
!  If multiple TLM history NetCDF files, close current NetCDF file.
!
        SourceFile=__FILE__ // ", increment"
        DO ng=1,Ngrids
          IF (TLM(ng)%ncid.ne.-1) THEN
            CALL netcdf_close (ng, iTLM, TLM(ng)%ncid)
            IF (FoundError(exit_flag, NoError, __LINE__,                &
     &                     __FILE__)) RETURN
          END IF
        END DO
# endif
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Time step adjoint model backwards: compute cost function gradient.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!  Initialize the adjoint model always from rest.
!
        DO ng=1,Ngrids
          CALL ad_initial (ng)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END DO
!
!  Time-step adjoint model backwards. The adjoint model is forced with
!  the adjoint of the observation misfit (Jo) term.
!
        DO ng=1,Ngrids
          IF (Master) THEN
            WRITE (stdout,30) 'AD', ng, ntstart(ng), ntend(ng)
          END IF
        END DO
!
# ifdef SOLVE3D
        CALL ad_main3d (RunInterval)
# else
        CALL ad_main2d (RunInterval)
# endif
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN
!
!  Clear adjoint arrays.
!
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL initialize_ocean (ng, tile, iADM)
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
            CALL initialize_forces (ng, tile, iADM)
# endif
# ifdef ADJUST_BOUNDARY
            CALL initialize_boundary (ng, tile, iADM)
# endif
          END DO
        END DO
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Descent algorithm.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!  Read TLM v-space initial conditions, record 3 in ITL(ng)%name, and
!  load it into time index LTLM1. This is needed to compute background
!  cost function. Also read in new (x-space) gradient vector, GRADx(Jo),
!  from adjoint history file ADM(ng)%name.  Read in the sum of all the
!  previous outer-loop increments which are always in record 4 of
!  the ITL file.
!
        DO ng=1,Ngrids
          IF (my_inner.eq.0) THEN
            CALL get_state (ng, iTLM, 8, ITL(ng)%name, Rec1, LTLM1)
            IF (FoundError(exit_flag, NoError, __LINE__,                &
     &                     __FILE__)) RETURN
          ELSE
            CALL get_state (ng, iTLM, 8, ITL(ng)%name, Rec3, LTLM1)
            IF (FoundError(exit_flag, NoError, __LINE__,                &
     &                     __FILE__)) RETURN
          END IF

          CALL get_state (ng, iTLM, 8, ITL(ng)%name, Rec4, LTLM2)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN

          CALL get_state (ng, iADM, 4, ADM(ng)%name, ADM(ng)%Rindex,    &
     &                    LADJ2)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN

# ifdef BALANCE_OPERATOR
          CALL get_state (ng, iNLM, 2, INI(ng)%name, Lini, Lini)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
          nrhs(ng)=Lini
# endif
        END DO
!
!  Convert observation cost function gradient, GRADx(Jo), from model
!  space (x-space) to minimization space (v-space):
!
!     GRADv(Jo) = B^(T/2) GRADx(Jo),  operator: S G L^(T/2) W^(-1/2)
!
!  First, multiply the adjoint solution, GRADx(Jo), by the background-
!  error standard deviations, S.  Second, convolve result with the
!  adjoint diffusion operator, G L^(T/2) W^(-1/2). Then, backgound
!  cost function contribution (BackCost) and cost function gradient
!  (v-space) by adding background and observation contributions:
!
!     GRADv(J) = GRADv(Jb) + GRADv(Jo) = deltaV + GRADv(Jo)
!
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1
# ifdef BALANCE_OPERATOR
            CALL ad_balance (ng, tile, Lini, LADJ2)
# endif
            CALL ad_variability (ng, tile, LADJ2, Lweak)
            CALL ad_convolution (ng, tile, LADJ2, Lweak, 2)
            CALL cost_grad (ng, tile, LTLM1, LTLM2, LADJ2)
          END DO
        END DO
!
!  Compute current total cost function.
!
        DO ng=1,Ngrids
          IF (Nrun.eq.1) THEN
            DO i=0,NobsVar(ng)
              FOURDVAR(ng)%CostFunOld(i)=FOURDVAR(ng)%CostNorm(i)
              FOURDVAR(ng)%CostFun(i)=FOURDVAR(ng)%CostNorm(i)
            END DO
          ELSE
            DO i=0,NobsVar(ng)
              FOURDVAR(ng)%CostFunOld(i)=FOURDVAR(ng)%CostFun(i)
            END DO
          END IF
        END DO
!
!  Prepare for background cost function (Jb) calculation:
!
!  Read the convolved gradient from inner=0 (which is permanently
!  saved in record 1 of the adjoint file)  ALWAYS into record 1.
!
        IF (my_inner.gt.0) THEN
          DO ng=1,Ngrids
            CALL get_state (ng, iADM, 3, ADM(ng)%name, LADJ1, LADJ1)
            IF (FoundError(exit_flag, NoError, __LINE__,                &
     &                     __FILE__)) RETURN
          END DO
        END IF
!
!  Compute background cost function (Jb) for inner=0:
!
!  If first pass of inner loop, read in the sum of previous v-space
!  gradients from record 4 of ITL file using the TLM model variables
!  as temporary storage. Also add background cost function to Cost0.
!
        IF (my_inner.eq.0) THEN
          DO ng=1,Ngrids
            CALL get_state (ng, iTLM, 2, ITL(ng)%name, Rec4, LTLM2)
            IF (FoundError(exit_flag, NoError, __LINE__,                &
     &                     __FILE__)) RETURN
!
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL back_cost (ng, tile, LTLM2)
            END DO
!
            FOURDVAR(ng)%Cost0(my_outer)=FOURDVAR(ng)%Cost0(my_outer)+  &
     &                                   FOURDVAR(ng)%BackCost(0)
          END DO
        END IF
!
!  Compute current total cost function.
!
        DO ng=1,Ngrids
          IF (Nrun.eq.1) THEN
            DO i=0,NobsVar(ng)
              FOURDVAR(ng)%CostNorm(i)=FOURDVAR(ng)%CostNorm(i)+        &
     &                                 FOURDVAR(ng)%BackCost(i)
              FOURDVAR(ng)%CostFunOld(i)=FOURDVAR(ng)%CostNorm(i)
              FOURDVAR(ng)%CostFun(i)=FOURDVAR(ng)%CostNorm(i)
            END DO
          ELSE
            DO i=0,NobsVar(ng)
              FOURDVAR(ng)%CostFunOld(i)=FOURDVAR(ng)%CostFun(i)
            END DO
          END IF
        END DO
!
!  Determine the descent direction in which the quadractic total cost
!  function decreases. Then, determine the TLM initial conditions,
!  deltaV(LTLM1), and its gradient, GRADv{J(Lnew)} at the new
!  direction.  Also, Compute TLM v-space trial initial conditions for
!  next inner loop, deltaV(LTLM2). The new gradient minimize the
!  quadratic cost function spanned by current and previous inner loop
!  iterations.  This is achieved by orthogonalizing (Gramm-Schmidt
!  algorithm) against all previous inner loop gradients.
!
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL cgradient (ng, tile, iTLM, my_inner, my_outer)
          END DO
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END DO
!
!  Report background (Jb) and observations (Jo) cost function values
!  normalized by their first minimization value. It also reports the
!  percentage change on total cost function value with respect to
!  previous iteration. Compute the optimality of the minimization to
!  check the statistical hypotheses between the background and
!  observations errors: the cost function value at the minimum, Jmin,
!  is idealy equal to half the number of observations assimilated
!  (Optimality=1=2*Jmin/Nobs), for a linear system.
!
        IF (Master) THEN
          DO ng=1,Ngrids
            IF (Nrun.gt.1) THEN
              rate=100.0_r8*ABS(FOURDVAR(ng)%CostFun(0)-                &
     &                          FOURDVAR(ng)%CostFunOld(0))/            &
     &                      FOURDVAR(ng)%CostFunOld(0)
            ELSE
              rate=0.0_r8
            END IF
            Optimality(ng)=2.0_r8*FOURDVAR(ng)%CostFun(0)/              &
     &                     (FOURDVAR(ng)%ObsCount(0)-                   &
     &                      FOURDVAR(ng)%ObsReject(0))
            WRITE (stdout,40) my_outer, my_inner,                       &
     &                        FOURDVAR(ng)%BackCost(0)/                 &
     &                        FOURDVAR(ng)%CostNorm(0),                 &
     &                        FOURDVAR(ng)%ObsCost(0)/                  &
     &                        FOURDVAR(ng)%CostNorm(0),                 &
     &                        rate
            IF (my_inner.eq.0) THEN
              DO i=0,NobsVar(ng)
                IF (FOURDVAR(ng)%NLobsCost(i).ne.0.0_r8) THEN
                  IF (i.eq.0) THEN
                    WRITE (stdout,50) my_outer, my_inner,               &
     &                                FOURDVAR(ng)%NLobsCost(i)/        &
     &                                FOURDVAR(ng)%CostNorm(i)
                  ELSE
                    WRITE (stdout,60) my_outer, my_inner,               &
     &                                FOURDVAR(ng)%NLobsCost(i)/        &
     &                                FOURDVAR(ng)%CostNorm(i),         &
     &                                TRIM(ObsName(i))
                  END IF
                END IF
                FOURDVAR(ng)%NLobsCost(i)=0.0
              END DO
            END IF
            WRITE (stdout,70) my_outer, my_inner, Optimality(ng)
          END DO
        END IF
!
!  Save total v-space cost function gradient, GRADv{J(Lnew)}, into
!  ADM(ng)%name history NetCDF file. Noticed that the lastest adjoint
!  solution record is over-written in the NetCDF file for future use.
!  The switch "LwrtState2d" is activated to write out state arrays
!  instead ad_*_sol arrays.
!
        DO ng=1,Ngrids
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
          Lfout(ng)=LADJ2
# endif
# ifdef ADJUST_BOUNDARY
          Lbout(ng)=LADJ2
# endif
          kstp(ng)=LADJ2
# ifdef SOLVE3D
          nstp(ng)=LADJ2
# endif
          ADM(ng)%Rindex=ADM(ng)%Rindex-1
          LwrtState2d(ng)=.TRUE.
          CALL ad_wrt_his (ng)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
          LwrtState2d(ng)=.FALSE.
        END DO
!
!  Write out trial v-space TLM initial conditions, currently in time
!  index LTM2, into record 3 of ITL(ng)%name NetCDF file.
!
        DO ng=1,Ngrids
          CALL tl_wrt_ini (ng, LTLM2, Rec3)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END DO
!
!  Read current outer loop nonlinear model initial conditions and
!  background state vectors.
!
        DO ng=1,Ngrids
          CALL get_state (ng, iNLM, 2, INI(ng)%name, Lini, Lini)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN

          CALL get_state (ng, iNLM, 9, INI(ng)%name, Lbck, Lbck)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END DO
!
!  Convert increment vector, deltaV, from minimization space (v-space)
!  to model space (x-space):
!
!     deltaX = B^(1/2) deltaV
!  or
!     deltaX = W^(1/2) L^(1/2) G S
!
!  First, convolve estimated increment vector (v-space) by with the
!  tangent linear diffusion operator, W^(1/2) L^(1/2) G.  Second,
!  multiply result by the background-error standard deviation, S.
!
        Lcon=LTLM2
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL tl_convolution (ng, tile, Lcon, Lweak, 2)
            CALL tl_variability (ng, tile, Lcon, Lweak)
# ifdef BALANCE_OPERATOR
            CALL tl_balance (ng, tile, Lini, Lcon)
# endif
          END DO
        END DO
!
!  Write out trial x-space (convolved) TLM initial conditions, currently
!  in time index Lcon, into record 1 of ITL(ng)%name NetCDF file.
!
        DO ng=1,Ngrids
          CALL tl_wrt_ini (ng, Lcon, Rec1)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END DO
!
!-----------------------------------------------------------------------
!  Update counters.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          Lsav=Lnew(ng)
          Lnew(ng)=Lold(ng)
          Lold(ng)=Lsav
          Nrun=Nrun+1
        END DO
!
      END DO INNER_LOOP
!
!  Close adjoint NetCDF file.
!
      SourceFile=__FILE__ // ", increment"
      DO ng=1,Ngrids
        IF (ADM(ng)%ncid.ne.-1) THEN
          CALL netcdf_close (ng, iADM, ADM(ng)%ncid)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
      END DO
!
!  Close Hessian NetCDF file.
!
      DO ng=1,Ngrids
        IF (HSS(ng)%ncid.ne.-1) THEN
          CALL netcdf_close (ng, iADM, HSS(ng)%ncid)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
      END DO

# ifdef PROFILE
!
!  Stop profile clock
!
      DO ng=1,Ngrids
        DO thread=THREAD_RANGE
          CALL wclock_off (ng, iTLM, 87, __LINE__,                      &
     &                     __FILE__)
        END DO
      END DO
# endif
!
 10   FORMAT (a,'_outer',i0,'.nc')
 20   FORMAT (a,'_member',i3.3,'.nc')
 30   FORMAT (/,1x,a,1x,'ROMS/TOMS: started time-stepping:',            &
     &        ' (Grid: ',i2.2,' TimeSteps: ',i8.8,' - ',i8.8,')',/)
 40   FORMAT (/,' (',i3.3,',',i3.3,'): TLM Cost Jb, J  = ',             &
     &        1p,e17.10,0p,1x,1p,e17.10,0p,t68,1p,e11.4,' %')
 50   FORMAT (/,'>(',i3.3,',',i3.3,'): NLM Cost     J  = ',             &
     &        18x,1p,e17.10,0p)
 60   FORMAT (' (',i3.3,',',i3.3,'): NLM Cost     J  = ',               &
     &        18x,1p,e17.10,0p,t69,a)
 70   FORMAT (/,1x,'(',i3.3,',',i3.3,'): Optimality (2*J/Nobs) = ',     &
     &        1p,e17.10,/)
!
      RETURN
      END SUBROUTINE increment
!
      SUBROUTINE analysis (my_outer, RunInterval)
!
!=======================================================================
!                                                                      !
!  This routine computes 4D-Var data assimilation analysis, Xa. The    !
!  nonlinear model initial conditions are computed by adding the       !
!  4D-Var increments to the current background:  Xa = Xb + dXa.        !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     my_outer        Outer-loop counter (integer)                     !
!     RunInterval     NLM kernel time stepping window (seconds)        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_iounits
      USE mod_scalars
      USE mod_stepping
!
      USE ini_adjust_mod,     ONLY : ini_adjust
      USE ini_fields_mod,     ONLY : ini_fields
      USE mod_ocean,          ONLY : initialize_ocean
      USE strings_mod,        ONLY : FoundError
      USE sum_grad_mod,       ONLY : sum_grad
      USE tl_convolution_mod, ONLY : tl_convolution
      USE tl_variability_mod, ONLY : tl_variability
      USE wrt_ini_mod,        ONLY : wrt_ini
!
!  Imported variable declarations
!
      logical :: Lweak = .FALSE.
!
      integer, intent(in) :: my_outer
!
      real(dp), intent(in) :: RunInterval
!
!  Local variable declarations.
!
      integer :: ng, tile
      integer :: Fcount, Lcon
# ifdef PROFILE
      integer :: thread
# endif
!
      SourceFile=__FILE__
!
!-----------------------------------------------------------------------
!  Compute new nonlinear initial conditions by adding minimization
!  increment to previous outer loop initial conditions:
!
!         Xi(outer+1) = Xi(outer) + deltaX(Lcon)
!
!-----------------------------------------------------------------------

# ifdef PROFILE
!
!  Start profile clock.
!
      DO ng=1,Ngrids
        DO thread=THREAD_RANGE
          CALL wclock_on (ng, iNLM, 88, __LINE__,                       &
     &                    __FILE__)
        END DO
      END DO
# endif
!
!  Clear nonlinear state variables.
!
      DO ng=1,Ngrids
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL initialize_ocean (ng, tile, iNLM)
        END DO
      END DO
!
!  Notice that "ini_fields" is called here for output purposes only.
!  It computes the vertically integrated momentum in 3D applications.
!  In order to use the correct fields, the model time indices are set
!  to Lini.
!  The appropriate tl correction for the NL model resides in record 1
!  of the ITL file.
!
      DO ng=1,Ngrids
        kstp(ng)=Lini
# ifdef SOLVE3D
        nstp(ng)=Lini
# endif
        CALL get_state (ng, iNLM, 1, INI(ng)%name, Lini, Lini)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN

        CALL get_state (ng, iTLM, 1, ITL(ng)%name, LTLM1, LTLM1)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN
!
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL ini_adjust (ng, tile, LTLM1, Lini)
          CALL ini_fields (ng, tile, iNLM)
        END DO
      END DO
!
!  Write out new nonlinear model initial conditions into record Lini
!  of INI(ng)%name.
!
      DO ng=1,Ngrids
        INI(ng)%Rindex=0
        Fcount=INI(ng)%load
        INI(ng)%Nrec(Fcount)=1
        CALL wrt_ini (ng, Lini)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN
      END DO
!
! Gather the v-space increments from the final inner-loop and
! save in record 4 of the ITL file. The current v-space increment
! is in record 3 and the sum so far is in record 4.
!
      DO ng=1,Ngrids
        CALL get_state (ng, iTLM, 8, ITL(ng)%name, Rec3, LTLM1)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN

        CALL get_state (ng, iTLM, 8, ITL(ng)%name, Rec4, LTLM2)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN
!
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL sum_grad (ng, tile, LTLM1, LTLM2)
        END DO
      END DO
!
! Write the current sum into record 4 of the ITL file.
!
      DO ng=1,Ngrids
        CALL tl_wrt_ini (ng, LTLM2, Rec4)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN
      END DO

# if defined ADJUST_STFLUX   || defined ADJUST_WSTRESS || \
     defined ADJUST_BOUNDARY
!
!  Set index containing the surface forcing increments used the run
!  the nonlinear model in the outer loop and read the forcing
!  increments. For bulk fluxes, we read Rec1 because the stress
!  fluxes change by virtue of the changing initial conditions.
!  When not using bulk fluxes, we read Rec4 because the background
!  stress and flux is prescribed by input files which are not
!  overwritten so we need to correct the background using the
!  sum of the increments from all previous outer-loops.
!  If using Rec4 we need to convert from v-space to x-space
!  by applying the convolution.
!  Note that Lfinp=Lbinp so the the forcing and boundary
!  adjustments are both processsed correctly.
#  ifdef BALANCE_OPERATOR
!  Currently, We do not need the call to tl_balance below, but we
!  might later if we impose a balance constraint on the wind stress
!  corrections.
#  endif
!
!  AMM: CHECK WHAT HAPPENS WITH SECONDARY PRECONDITIONING.
!
      DO ng=1,Ngrids
        Lfinp(ng)=LTLM1
#  if defined BULK_FLUXES && !defined NL_BULK_FLUXES
        CALL get_state (ng, iTLM, 1, ITL(ng)%name, Rec1, Lfinp(ng))
#  endif
#  if defined NL_BULK_FLUXES || !defined BULK_FLUXES
        CALL get_state (ng, iTLM, 1, ITL(ng)%name, Rec4, Lfinp(ng))
        Lcon=Lfinp(ng)
!
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL tl_convolution (ng, tile, Lcon, Lweak, 2)
          CALL tl_variability (ng, tile, Lcon, Lweak)
#   ifdef BALANCE_OPERATOR
!!        CALL tl_balance (ng, tile, Lini, Lcon)
#   endif
        END DO
#  endif
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN
      END DO
# endif
!
!  Clear tangent linear state variables.
!
      DO ng=1,Ngrids
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL initialize_ocean (ng, tile, iTLM)
        END DO
      END DO
!
!  Close current forward NetCDF file.
!
      SourceFile=__FILE__ // ", analysis"
      DO ng=1,Ngrids
!!      CALL netcdf_close (ng, iNLM, FWD(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN
      END DO

# ifdef PROFILE
!
!  Stop profile clock
!
      DO ng=1,Ngrids
        DO thread=THREAD_RANGE
          CALL wclock_off (ng, iNLM, 88, __LINE__,                      &
     &                     __FILE__)
        END DO
      END DO
# endif
!
      RETURN
      END SUBROUTINE analysis
!
      SUBROUTINE posterior_analysis (RunInterval)
!
!=======================================================================
!                                                                      !
!  This routine initialize the NLM with estimated 4D-Var state and     !
!  interpolates solution at observation locations for posterior        !
!  analysis.                                                           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     RunInterval     NLM kernel time stepping window (seconds)        !
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
!
      USE mod_mixing,  ONLY : initialize_mixing
      USE strings_mod, ONLY : FoundError
!
!  Imported variable declarations
!
      real(dp), intent(in) :: RunInterval
!
!  Local variable declarations.
!
      integer :: i, lstr, ng, tile
      integer :: Fcount
!
      SourceFile=__FILE__
!
!-----------------------------------------------------------------------
!  Run NLM initialized from estimated 4D-Var state and interpolate
!  solution at observation locations.
!-----------------------------------------------------------------------
!
!  Set nonlinear output history file name. Create a basic state file
!  for each outher loop. Notice that the LreadBLK and LreadFWD switches
!  are turned off to suppress processing of the structures when
!  "check_multifile" during nonlinear model execution.
!
      DO ng=1,Ngrids
        idefHIS(ng)=-1
        LdefHIS(ng)=.TRUE.
        LwrtHIS(ng)=.TRUE.
# if defined BULK_FLUXES && defined NL_BULK_FLUXES
        LreadBLK(ng)=.FALSE.
# endif
        LreadFWD(ng)=.FALSE.
        HIS(ng)%Rindex=0
        Fcount=HIS(ng)%load
        HIS(ng)%Nrec(Fcount)=0
        WRITE (HIS(ng)%name,10) TRIM(FWD(ng)%head), Nouter
        lstr=LEN_TRIM(HIS(ng)%name)
        HIS(ng)%base=HIS(ng)%name(1:lstr-3)
      END DO
!
!  Clear nonlinear mixing arrays.
!
      DO ng=1,Ngrids
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL initialize_mixing (ng, tile, iNLM)
        END DO
      END DO
!
!  Initialize nonlinear model with estimated initial conditions.
!
      DO ng=1,Ngrids
        wrtNLmod(ng)=.TRUE.
        wrtTLmod(ng)=.FALSE.
        wrtMisfit(ng)=.FALSE.
        RST(ng)%Rindex=0
        Fcount=RST(ng)%load
        RST(ng)%Nrec(Fcount)=0
      END DO
!
      CALL initial
      IF (FoundError(exit_flag, NoError, __LINE__,                      &
     &               __FILE__)) RETURN
!
!  Clear NLobsCost.
!
      DO ng=1,Ngrids
        DO i=0,NobsVar(ng)
          FOURDVAR(ng)%NLobsCost(i)=0.0_r8
        END DO
      END DO
!
!  Run nonlinear model. Interpolate nonlinear model to observation
!  locations.
!
      DO ng=1,Ngrids
# ifdef AVERAGES
        idefAVG(ng)=-1
        LdefAVG(ng)=.TRUE.
        LwrtAVG(ng)=.TRUE.
        WRITE (AVG(ng)%name,10) TRIM(AVG(ng)%head), outer
        lstr=LEN_TRIM(AVG(ng)%name)
        AVG(ng)%base=AVG(ng)%name(1:lstr-3)
# endif
# ifdef DIAGNOSTICS
        idefDIA(ng)=-1
        LdefDIA(ng)=.TRUE.
        LwrtDIA(ng)=.TRUE.
        WRITE (DIA(ng)%name,10) TRIM(DIA(ng)%head), outer
        lstr=LEN_TRIM(DIA(ng)%name)
        DIA(ng)%base=DIA(ng)%name(1:lstr-3)
# endif
        IF (Master) THEN
          WRITE (stdout,20) 'NL', ng, ntstart(ng), ntend(ng)
        END IF
      END DO
!
# ifdef SOLVE3D
      CALL main3d (RunInterval)
# else
      CALL main2d (RunInterval)
# endif
      IF (FoundError(exit_flag, NoError, __LINE__,                      &
     &               __FILE__)) RETURN
!
!  Write out nonlinear model final misfit cost function into
!  DAV(ng)%name NetCDF file. Notice that it is written in the
!  Nouter+1 record.
!
      SourceFile=__FILE__ // ", posterior_analysis"
      DO ng=1,Ngrids
        CALL netcdf_put_fvar (ng, iNLM, DAV(ng)%name, 'NLcost_function',&
     &                        FOURDVAR(ng)%NLobsCost(0:),               &
     &                        (/1,Nouter+1/), (/NobsVar(ng)+1,1/),      &
     &                        ncid = DAV(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN
      END DO
!
!  Report the final value of the nonlinear model misfit cost function.
!
      IF (Master) THEN
        DO ng=1,Ngrids
          DO i=0,NobsVar(ng)
            IF (FOURDVAR(ng)%NLobsCost(i).ne.0.0_r8) THEN
              IF (i.eq.0) THEN
                WRITE (stdout,30) outer, inner,                         &
     &                            FOURDVAR(ng)%NLobsCost(i)/            &
     &                            FOURDVAR(ng)%CostNorm(i)
              ELSE
                WRITE (stdout,40) outer, inner,                         &
     &                            FOURDVAR(ng)%NLobsCost(i)/            &
     &                            FOURDVAR(ng)%CostNorm(i),             &
     &                            TRIM(ObsName(i))
              END IF
            END IF
          END DO
        END DO
      END IF
!
!  Done.  Set history file ID to closed state since we manipulated
!  its indices with the forward file ID which was closed above.
!
      DO ng=1,Ngrids
        HIS(ng)%ncid=-1
      END DO
!
 10   FORMAT (a,'_outer',i0,'.nc')
 20   FORMAT (/,1x,a,1x,'ROMS/TOMS: started time-stepping:',            &
     &        ' (Grid: ',i2.2,' TimeSteps: ',i8.8,' - ',i8.8,')',/)
 30   FORMAT (/,'>(',i3.3,',',i3.3,'): NLM Cost     J  = ',             &
     &        18x,1p,e17.10,0p)
 40   FORMAT (' (',i3.3,',',i3.3,'): NLM Cost     J  = ',               &
     &        18x,1p,e17.10,0p,t69,a)
!
      RETURN
      END SUBROUTINE posterior_analysis
!
      SUBROUTINE prior_error (ng)
!
!=======================================================================
!                                                                      !
!  This routine processes background prior error covariance standard   !
!  deviations and normalization coefficients.                          !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng              Nested grid number                               !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      USE normalization_mod, ONLY : normalization
      USE strings_mod,       ONLY : FoundError
!
!  Imported variable declarations
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      integer :: tile
      integer :: NRMrec, STDrec, Tindex
!
      SourceFile=__FILE__
!
!-----------------------------------------------------------------------
!  Read in standard deviation factors for error covariance.
!-----------------------------------------------------------------------
!
!  Initial conditions standard deviation. They are loaded in Tindex=1
!  of the e_var(...,Tindex) state variables.
!
      STDrec=1
      Tindex=1
      CALL get_state (ng, 10, 10, STD(1,ng)%name, STDrec, Tindex)
      IF (FoundError(exit_flag, NoError, __LINE__,                      &
     &               __FILE__)) RETURN

# ifdef ADJUST_BOUNDARY
!
!  Open boundary conditions standard deviation.
!
      STDrec=1
      Tindex=1
      CALL get_state (ng, 12, 12, STD(3,ng)%name, STDrec, Tindex)
      IF (FoundError(exit_flag, NoError, __LINE__,                      &
     &                 __FILE__)) RETURN
# endif

# if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
!
!  Surface forcing standard deviation.
!
      STDrec=1
      Tindex=1
      CALL get_state (ng, 13, 13, STD(4,ng)%name, STDrec, Tindex)
      IF (FoundError(exit_flag, NoError, __LINE__,                      &
     &               __FILE__)) RETURN
# endif
!
!-----------------------------------------------------------------------
!  Error covariance normalization coefficients.
!-----------------------------------------------------------------------
!
!  Compute or read in background-error covariance normalization factors.
!  If computing, write out factors to NetCDF. This is an expensive
!  computation that needs to be computed only once for a particular
!  application grid.
!
      IF (ANY(LwrtNRM(:,ng))) THEN
        CALL def_norm (ng, iNLM, 1)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN

# ifdef ADJUST_BOUNDARY
        CALL def_norm (ng, iNLM, 3)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN
# endif

# if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
        CALL def_norm (ng, iNLM, 4)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN
# endif

        DO tile=first_tile(ng),last_tile(ng),+1
          CALL normalization (ng, tile, 2)
        END DO
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN
        LdefNRM(1:4,ng)=.FALSE.
        LwrtNRM(1:4,ng)=.FALSE.

      ELSE

        NRMrec=1
        CALL get_state (ng, 14, 14, NRM(1,ng)%name, NRMrec, 1)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN

# ifdef ADJUST_BOUNDARY
        CALL get_state (ng, 16, 16, NRM(3,ng)%name, NRMrec, 1)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN
# endif

# if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
        CALL get_state (ng, 17, 17, NRM(4,ng)%name, NRMrec, 1)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN
# endif

      END IF
!
      RETURN
      END SUBROUTINE prior_error
#endif
      END MODULE i4dvar_mod
