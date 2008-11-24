#include "cppdefs.h"

#if defined W4DVAR || defined W4DPSAS
      SUBROUTINE congrad (ng, model, outLoop, innLoop, NinnLoop,        &
     &                    converged)
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Weak Constraint 4-Dimensional Variational (4DVar) Pre-conditioned   !
!                     Conjugate Gradient Algorithm                     !
!                                                                      !
!  The indirect representer method solves the system:                  !
!                                                                      !
!              (R_n + Cobs) * Beta_n = h_n                             !
!                                                                      !
!              h_n = Xo - H * X_n                                      !
!                                                                      !
!  where R_n is the representer matrix, Cobs is the observation-error  !
!  covariance,  Beta_n  are the representer coefficients,  h_n is the  !
!  misfit between observations (Xo) and model (H*X_n),  and  H is the  !
!  linearized observation operator. Here, _n denotes iteration.        !
!                                                                      !
!  This system does not need to be solved explicitly by inverting the  !
!  symmetric stabilized representer matrix, P_n:                       !
!                                                                      !
!              P_n = R_n + Cobs                                        !
!                                                                      !
!  but by computing the action of P_n on any vector PSI, such that     !
!                                                                      !
!              P_n * PSI = R_n * PSI + Cobs * PSI                      !
!                                                                      !
!  The representer matrix is not explicitly computed but evaluated by  !
!  one integration backward of the adjoint model  and one integration  !
!  forward of the tangent linear model for any forcing vector PSI.     !
!                                                                      !
!  Notice that "ObsScale" vector is used for screenning observations.  !
!  This scale is one (zero) for good (bad) observations.               !
!                                                                      !
!  Currently, parallelization of this algorithm is not needed because  !
!  each parallel node has a full copy of the assimilation vectors.     !
!                                                                      !
!  This code solves Ax=b by minimizing the cost function               !
!  0.5*x*A*x-x*b, assuming an initial guess of x=0. In this case the   !
!  gradient is Ax-b and the Hessian is A.                              !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    Chua, B. S. and A. F. Bennett,  2001:  An inverse ocean modeling  !
!      sytem, Ocean Modelling, 3, 137-165.                             !
!                                                                      !
!  Lanczos Algorithm Reference:                                        !
!                                                                      !
!    Fisher, M., 1998: Minimization Algorithms for Variational Data    !
!      Assimilation. In Seminar on Recent Developments in Numerical    !
!      Methods for Atmospheric Modelling, 1998.                        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_iounits
      USE mod_scalars

# ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcastf, mp_bcasti, mp_bcastl
# endif
      implicit none
!
!  Imported variable declarations
!
      integer, intent(in) :: ng, model, outLoop, innLoop, NinnLoop
      logical, intent(out) :: converged
!
!  Local variable declarations.
!
      integer :: iobs, ivec
      real(r8), dimension(NinnLoop) :: zu, zgam
      real(r8), dimension(Ndatum(ng)) :: px, pgrad, zw
      real(r8), dimension(NinnLoop,3) :: zwork
      real(r8) :: dla, zbet, eps, preduc
      real(r8) :: Jopt, Jf, Jmod, Jdata
!
!-----------------------------------------------------------------------
!  Weak constraint 4DVar conjugate gradient, Lanczos-based algorithm.
!-----------------------------------------------------------------------
!
!  This conjugate gradient algorithm is not run in parallel since the
!  weak constraint is done in observation space. Mostly all the
!  variables are 1D arrays. Therefore, in parallel applications (only
!  distributed-memory is possible) the master node does the computations
!  and then broadcasts results to remaining nodes in the communicator.
!
      IF (Master) THEN
!
!  Initialize internal parameters.  This is needed here for output
!  reasons.
!
        Jopt=0.0_r8
        Jf=0.0_r8
        Jdata=0.0_r8
        Jmod=0.0_r8
        preduc=0.0_r8

        IF (innLoop.eq.0) THEN
!
!  If this is the first pass of the inner loop, set up the vectors and
!  store the first guess. The initial starting guess is assumed to be
!  zero in which case the gradient is just: -(obs-model).
!
          zgnorm=0.0
          DO iobs=1,Ndatum(ng)
# ifdef W4DPSAS
            pgrad(iobs)=-ObsScale(iobs)*                                &
     &                  (ObsVal(iobs)-NLmodVal(iobs))
# else
            pgrad(iobs)=-ObsScale(iobs)*                                &
     &                  (ObsVal(iobs)-TLmodVal(iobs))
# endif
            zgrad0(iobs)=pgrad(iobs)
            zgnorm=zgnorm+zgrad0(iobs)*zgrad0(iobs)
          END DO
          zgnorm=SQRT(zgnorm)
          DO iobs=1,Ndatum(ng)
            zcglwk(iobs,1)=pgrad(iobs)/zgnorm
            ADmodVal(iobs)=zcglwk(iobs,1)
          END DO
          zqg0(1)=0.0_r8
          DO iobs=1,Ndatum(ng)
            zqg0(1)=zqg0(1)+zcglwk(iobs,1)*zgrad0(iobs)
          END DO

        ELSE
!
!  After the initialization, for every other inner loop, calculate a
!  new Lanczos vector, store it in the matrix, and update search.
!
          zdelta(innLoop)=0.0_r8
          DO iobs=1,Ndatum(ng)
            pgrad(iobs)=ObsScale(iobs)*                                 &
     &                 (TLmodVal(iobs)+                                 &
     &                  ObsErr(iobs)*zcglwk(iobs,innLoop))
            zdelta(innLoop)=zdelta(innLoop)+                            &
     &                      zcglwk(iobs,innLoop)*pgrad(iobs)
          END DO
!
!  Exit, if not a positive definity algorithm.
!
          IF (zdelta(innLoop).le.0.0_r8) THEN
            WRITE (stdout,20) zdelta(innLoop), outLoop, innLoop
            exit_flag=7
            GO TO 10
          END IF
!
!  Compute the new Lanczos vector.
!
          DO iobs=1,Ndatum(ng)
            pgrad(iobs)=pgrad(iobs)-                                    &
     &                  zdelta(innLoop)*zcglwk(iobs,innLoop)
          END DO
          IF (innLoop.gt.1) THEN
            DO iobs=1,Ndatum(ng)
              pgrad(iobs)=pgrad(iobs)-                                  &
     &                    zbeta(innLoop)*zcglwk(iobs,innLoop-1)
            END DO
          END IF
!
!  Orthonormalize against previous Lanczos vectors.
!
          DO ivec=innLoop,1,-1
            dla=0.0_r8
            DO iobs=1,Ndatum(ng)
              dla=dla+pgrad(iobs)*zcglwk(iobs,ivec)
            END DO
            DO iobs=1,Ndatum(ng)
              pgrad(iobs)=pgrad(iobs)-dla*zcglwk(iobs,ivec)
            END DO
          END DO
!
          zbeta(innLoop+1)=0.0_r8
          DO iobs=1,Ndatum(ng)
            zbeta(innLoop+1)=zbeta(innLoop+1)+pgrad(iobs)*pgrad(iobs)
          END DO
          zbeta(innLoop+1)=SQRT(zbeta(innLoop+1))
!
          DO iobs=1,Ndatum(ng)
            zcglwk(iobs,innLoop+1)=pgrad(iobs)/zbeta(innLoop+1)
          END DO
!
          zqg0(innLoop+1)=0.0_r8
          DO iobs=1,Ndatum(ng)
            zqg0(innLoop+1)=zqg0(innLoop+1)+                            &
     &                      zcglwk(iobs,innLoop+1)*zgrad0(iobs)
          END DO
!
!  Calculate the new solution based upon the new, orthonormalized
!  Lanczos vector. First, the tridiagonal system is solved by
!  decomposition and forward/back substitution.
!
          zbet=zdelta(1)
          zu(1)=-zqg0(1)/zbet
          DO ivec=2,innLoop
            zgam(ivec)=zbeta(ivec)/zbet
            zbet=zdelta(ivec)-zbeta(ivec)*zgam(ivec)
            zu(ivec)=(-zqg0(ivec)-zbeta(ivec)*zu(ivec-1))/zbet
          END DO
          zwork(innLoop,3)=zu(innLoop)

          DO ivec=innLoop-1,1,-1
            zu(ivec)=zu(ivec)-zgam(ivec+1)*zu(ivec+1)
            zwork(ivec,3)=zu(ivec)
          END DO

          DO iobs=1,Ndatum(ng)
            zw(iobs)=zgrad0(iobs)+                                      &
     &               zbeta(innLoop+1)*zcglwk(iobs,innLoop+1)*           &
     &               zwork(innLoop,3)
          END DO

          DO iobs=1,Ndatum(ng)
            px(iobs)=0.0
            DO ivec=1,innLoop
              px(iobs)=px(iobs)+zcglwk(iobs,ivec)*zwork(ivec,3)
              zw(iobs)=zw(iobs)-zcglwk(iobs,ivec)*zqg0(ivec)
            END DO
          END DO
!
!  Compute the reduction in the gradient Ax-b.
!
          preduc=0.0_r8
          DO iobs=1,Ndatum(ng)
            preduc=preduc+zw(iobs)*zw(iobs)
          END DO
          preduc=SQRT(preduc)/zgnorm
!
!  Estimate the residual by: ||Ax-b|| = zbeta(k+1)
!
          eps=ABS(zbeta(innLoop+1))
!
!  Estimate the various penalty function components based on
!  Bennett (2002), page 22 and pages 42-44:
!
!     Jopt  = value of the penalty function when true px
!             has been identified. 
!     Jf    = initial data misfit for the first-guess
!     Jdata = data misfit for the state estimate
!     Jmod  = the model penalty function (the sum of the
!             dynamical, initial and boundary penalties)     
!
!  NOTE:  Unlike IS4DVAR where we compute the actual cost function
!         for each iteration that decreases with increasing 
!         number of iterations, the various penalty terms 
!         represent the optimal values (i.e. the values when
!         optimal state estimate has been identified). Therefore
!         the various penalty terms will only be correct when the
!         optimal state has been found. Therefore, as the
!         W4DVAR proceeds, Jopt, Jdata and Jmod will asymptote
!         to their final values and WILL NOT necessary decrease
!         with increasing iteration number.
!
          DO iobs=1,Ndatum(ng)
            Jopt=Jopt-px(iobs)*zgrad0(iobs)
            IF (ObsErr(iobs).ne.0.0_r8) THEN
              Jf=Jf+zgrad0(iobs)*zgrad0(iobs)/ObsErr(iobs)
            END IF
            Jdata=Jdata+px(iobs)*px(iobs)*ObsErr(iobs)
          END DO
          Jmod=Jopt-Jdata
!
!  Check the convergence criteria: The Lanczos algorithm should not
!  be run for more loops than the size of the data. This will cause
!  problems. If we have performed more loops than we have datum,
!  then we set converged to TRUE so that the inner loop stops. 
!
          converged=.FALSE.
          IF ((preduc.le.CGeps).or.(innLoop.eq.Ndatum(ng))) THEN
            converged=.TRUE.
          END IF
!
!  Put the new trial solution into the adjoint vector for the next loop
!  Put the final solution into the adjoint vector when converged
!  of on the final inner-loop.
!
          IF ((innLoop.eq.NinnLoop).or.converged) THEN
            DO iobs=1,Ndatum(ng)
              ADmodVal(iobs)=px(iobs)
            END DO
          ELSE
            DO iobs=1,Ndatum(ng)
              ADmodVal(iobs)=zcglwk(iobs,innLoop+1)
            END DO
          END IF
        END IF
!
!  Report minimization progress.
!
        WRITE (stdout,30) Ndatum(ng), converged,                         &
     &                    outLoop, innLoop, eps,                         &
     &                    outLoop, innLoop, preduc,                      &
     &                    outLoop, innLoop, Jf,                          &
     &                    outLoop, innLoop, Jdata,                       &
     &                    outLoop, innLoop, Jmod,                        &
     &                    outLoop, innLoop, Jopt,                        &
     &                    outLoop, innLoop
        DO ivec=1,innLoop
          WRITE (stdout,40) ivec, zdelta(ivec), zbeta(ivec),             &
     &                      zwork(ivec,3)
        END DO
        WRITE (stdout,50) innLoop+1, zbeta(innLoop+1)

      END IF
!
!  Come here if not a possite definite algorithm.
!
 10   CONTINUE
# ifdef DISTRIBUTE
      CALL mp_bcasti (ng, model, exit_flag, 1)
# endif
      IF (exit_flag.ne.NoError) RETURN
   
# ifdef DISTRIBUTE
!
!  Broadcast new solution to other nodes.
!
      CALL mp_bcastf (ng, iTLM, ADmodVal, Mobs)
      CALL mp_bcastl (ng, iTLM, converged, 1)
# endif
!!
!! Write out conjugate gradient vectors into 4DVAR NetCDF file.
!!
!!    CALL cg_write (ng, converged)
!
 20   FORMAT (/,' CONGRAD - Fatal error, not possitive definite',       &
     &        ' algorithm:',/,                                          &
     &        /,11x,'zdelta = ',1p,e15.8,0p,3x,'(',i3.3,', ',i3.3,')')
 30   FORMAT (/,' CONGRAD - Conjugate Gradient Information:',/,         &
     &        /,11x,'Ndatum     = ',i7.7,                               &
     &        /,11x,'converged  = ',l1,/,/,                             &
     &        1x,'(',i3.3,',',i3.3,'): ',                               &
     &        'Residual estimate,                  eps = ',1p,e14.7,/,  &
     &        1x,'(',i3.3,',',i3.3,'): ',                               &
     &        'Reduction in gradient norm,      Greduc = ',1p,e14.7,/,  &
     &        1x,'(',i3.3,',',i3.3,'): ',                               &
     &        'First guess initial data misfit,     Jf = ',1p,e14.7,/,  &
     &        1x,'(',i3.3,',',i3.3,'): ',                               &
     &        'State estimate data misfit,       Jdata = ',1p,e14.7,/,  &
     &        1x,'(',i3.3,',',i3.3,'): ',                               &
     &        'Model penalty function,            Jmod = ',1p,e14.7,/,  &
     &        1x,'(',i3.3,',',i3.3,'): ',                               &
     &        'Optimal penalty function,          Jopt = ',1p,e14.7,/,  &
     &        /,1x,'(',i3.3,',',i3.3,'): ',                             &
     &        'Lanczos vectors - zdelta, zbeta, zwork:',/)
 40   FORMAT (6x,i3.3,4x,3(1p,e15.8,5x))
 50   FORMAT (6x,i3.3,24x,1p,e15.8)

      RETURN
      END SUBROUTINE congrad

      SUBROUTINE cg_write (ng, converged)
!
!=======================================================================
!                                                                      !
!  This routine writes conjugate gradient vectors into 4DVAR NetCDF    !
!  for restart purposes.                                               !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      Use mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations
!
      logical, intent(in) :: converged

      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      integer :: status, varid
      integer :: start(2), total(2)
!
!-----------------------------------------------------------------------
!  Write out conjugate gradient vectors.
!-----------------------------------------------------------------------
!
      IF (OutThread) THEN
!
!  Write out outer and inner iteration.
!
        status=nf90_inq_varid(ncMODid(ng), 'outer', varid)
        status=nf90_put_var(ncMODid(ng), varid, outer)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'outer', TRIM(MODname(ng))
          exit_flag=3
          ioerror=status
          RETURN
        END IF

        status=nf90_inq_varid(ncMODid(ng), 'inner', varid)
        status=nf90_put_var(ncMODid(ng), varid, inner)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'inner', TRIM(MODname(ng))
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
!  Write out representer coefficients.
!
        start(1)=1
        total(1)=Ndatum(ng)
        status=nf90_put_var(ncMODid(ng), modVid(idRepC,ng),             &
     &                      ADmodVal, start, total)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,10) TRIM(Vname(1,idRepC)), TRIM(MODname(ng))
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
!  Write out conjugate gradient vectors.
!
!  HERNAN -- Here we need to write out zdelta, zbeta, and zcglwk
!
       END IF
!
!-----------------------------------------------------------------------
!  Synchronize observations NetCDF file to disk.
!-----------------------------------------------------------------------
!
      IF (OutThread) THEN
        status=nf90_sync(ncMODid(ng))
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,20)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

  10  FORMAT (/,' CG_WRITE - error while writing variable: ',a,/,       &
     &        12x,'into NetCDF file: ',a)
  20  FORMAT (/,' CG_WRITE - unable to synchronize 4DVAR',              &
     &        1x,'NetCDF file to disk.')

      END SUBROUTINE cg_write
#else
      SUBROUTINE congrad
      RETURN
      END SUBROUTINE congrad
#endif
