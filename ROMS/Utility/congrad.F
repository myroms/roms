#include "cppdefs.h"

#if defined W4DVAR || defined W4DPSAS
      SUBROUTINE congrad (ng, outer, inner, Ninner, converged)
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2007 The ROMS/TOMS Group       Andrew M. Moore   !
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
!  Congugate gradient variables:                                       !
!                                                                      !
!    cg_r      Current PSI approximation for Beta_n.                   !
!                                                                      !
!    cg_s      Right-hand-side terms:  R_n * PSI + Cobs * PSI   where  !
!              R_n * PSI are the TLM or RPM values at the observation  !
!              points.                                                 !
!                                                                      !
!  Notice that "ObsScale" vector is used for screenning observations.  !
!  This scale is one (zero) for good (bad) observations.               !
!                                                                      !
!  Currently, parallelization of this algorithm is not needed because  !
!  each parallel node has a full copy of the assimilation vectors.     !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    Chua, B. S. and A. F. Bennett,  2001:  An inverse ocean modeling  !
!      sytem, Ocean Modelling, 3, 137-165.                             !
!                                                                      !
!  This code was adapted from IOM routine "congrad0".                  !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
!
      implicit none
!
!  Imported variable declarations
!
      integer, intent(in) :: ng, outer, inner, Ninner

      logical, intent(out) :: converged
!
!  Local variable declarations.
!
      integer :: iobs

      real(r8) :: cg_alpha, cg_beta, cg_delta, cg_eps
      real(r8) :: cg_gamma, cg_sigma
!
!-----------------------------------------------------------------------
!  If entering inner loop, initialize vector PSI needed to force the
!  adjoint model to the misfit between observations and representer
!  model, h_n.
!-----------------------------------------------------------------------
!
      IF (inner.eq.0) THEN
        DO iobs=1,Ndatum(ng)
          cg_p(iobs)=0.0_r8
          cg_s(iobs)=0.0_r8
          cg_v(iobs)=0.0_r8
          cg_x(iobs)=0.0_r8
          cg_z(iobs)=0.0_r8
# ifdef W4DPSAS
          cg_r(iobs)=ObsScale(iobs)*                                    &
     &               (ObsVal(iobs)-NLmodVal(iobs))
# else
          cg_r(iobs)=ObsScale(iobs)*                                    &
     &               (ObsVal(iobs)-TLmodVal(iobs))
# endif
        END DO
        cg_eps=0.0_r8
        cg_gamma=0.0_r8
        cg_sigma=0.0_r8
        cg_rnorm=0.0_r8
!
!-----------------------------------------------------------------------
!  First inner loop iteration.
!-----------------------------------------------------------------------
!
      ELSE IF (inner.eq.1) THEN
!
!  Initialize.
!
        DO iobs=1,Ndatum(ng)
          cg_z(iobs)=cg_r(iobs)
          cg_s(iobs)=ObsScale(iobs)*                                    &
     &               (TLmodVal(iobs)+ObsErr(iobs)*cg_r(iobs))
          cg_v(iobs)=cg_s(iobs)
          cg_p(iobs)=cg_z(iobs)                              ! p_1 = z_1
        END DO
!
!  Compute dot products:  assumes that z_1 = M^(-1) * r_1
!                                      v_1 = A * z_1
!
        cg_gamma=0.0_r8
        cg_sigma=0.0_r8
        cg_rnorm=0.0_r8
        DO iobs=1,Ndatum(ng)
          cg_gamma=cg_gamma+cg_z(iobs)*cg_r(iobs)            ! <z_1,r_1>
          cg_sigma=cg_sigma+cg_p(iobs)*cg_v(iobs)            ! <p_1,v_1>
          cg_rnorm=cg_rnorm+cg_r(iobs)*cg_r(iobs)            ! <b  ,b  >
        END DO
        cg_gammam1=cg_gamma
        cg_sigmam1=cg_sigma
        cg_eps=ABS(cg_gamma/cg_rnorm)
        cg_alpha=cg_gamma/cg_sigma
!
!  Update:   x_2 = x_1 + alpha * p_1,      recall x_1 = 0
!            r_2 = r_1 - alpha * v_1
!
        DO iobs=1,Ndatum(ng)
          cg_x(iobs)=cg_alpha*cg_p(iobs)
          cg_r(iobs)=cg_r(iobs)-cg_alpha*cg_v(iobs)
        END DO
!
!-----------------------------------------------------------------------
!  Compute a new approximation of vector PSI for Beta_n (cg_r) based on
!  previous iteration conjugate gradient vectors and new righ-hand-side
!  term: cg_s = R_n*PSI + Cobs*PSI evaluated by integrating the adjoint
!  and tangent linear model in the inner loop.
!-----------------------------------------------------------------------
!
      ELSE IF (inner.gt.1) THEN
!
!  Initialize.
!
        DO iobs=1,Ndatum(ng)
          cg_z(iobs)=cg_r(iobs)
          cg_s(iobs)=ObsScale(iobs)*                                    &
     &               (TLmodVal(iobs)+ObsErr(iobs)*cg_r(iobs))
        END DO
!
!  Compute dot products:  assumes that z_n = M^(-1) * r_n
!                                      s_n = A * z_n
!
        cg_gamma=0.0_r8
        cg_delta=0.0_r8
        DO iobs=1,Ndatum(ng)
          cg_gamma=cg_gamma+cg_z(iobs)*cg_r(iobs)            ! <z_n,r_n>
          cg_delta=cg_delta+cg_z(iobs)*cg_s(iobs)            ! <z_n,s_n>
        END DO
        cg_eps=ABS(cg_gamma/cg_rnorm)
        cg_beta=cg_gamma/cg_gammam1
        cg_sigma=cg_delta-cg_beta*cg_beta*cg_sigmam1
        cg_alpha=cg_gamma/cg_sigma
        cg_gammam1=cg_gamma
        cg_sigmam1=cg_sigma
!
!  Update.
!
        DO iobs=1,Ndatum(ng)
          cg_p(iobs)=cg_z(iobs)+cg_beta*cg_p(iobs)
          cg_v(iobs)=cg_s(iobs)+cg_beta*cg_v(iobs)
          cg_x(iobs)=cg_x(iobs)+cg_alpha*cg_p(iobs)
          cg_r(iobs)=cg_r(iobs)-cg_alpha*cg_v(iobs)
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Check convergence and load representer coefficients (or its
!  approximation into adjoint forcing vector ADmodVal.
!-----------------------------------------------------------------------
!
      IF (((cg_eps.lt.CGeps).and.(inner.gt.0)).or.                      &
     &    (inner.eq.Ninner)) THEN
        converged=.TRUE.
         DO iobs=1,Ndatum(ng)
          ADmodVal(iobs)=cg_x(iobs)
        END DO
      ELSE
        converged=.FALSE.
        DO iobs=1,Ndatum(ng)
          ADmodVal(iobs)=cg_r(iobs)
        END DO
      END IF
!
!  Write out conjugate gradient vectors into 4DVAR NetCDF file.
!
      CALL cg_write (ng, converged)
!
!  Report minimization parameters.
!
      IF (Master) THEN
        PRINT 10, outer, inner, Ndatum(ng), converged,                  &
     &            cg_eps,   outer, inner,                               &
     &            CGeps,    outer, inner,                               &
     &            cg_gamma, outer, inner,                               &
     &            cg_sigma, outer, inner,                               &
     &            cg_rnorm, outer, inner,                               &
     &            MINVAL(cg_x), MAXVAL(cg_x),                           &
     &            MINVAL(cg_z), MAXVAL(cg_z),                           &
     &            MINVAL(cg_p), MAXVAL(cg_p),                           &
     &            MINVAL(cg_v), MAXVAL(cg_v),                           &
     &            MINVAL(cg_s), MAXVAL(cg_s),                           &
     &            MINVAL(cg_r), MAXVAL(cg_r)
 10     FORMAT (/,' CONGRAD - Conjugate Gradient Information: ',/,      &
     &          /,11x,'Outer Loop = ',i3.3,                             &
     &          /,11x,'Inner Loop = ',i3.3,                             &
     &          /,11x,'Ndatum     = ',i7.7,                             &
     &          /,11x,'converged  = ',l1,/,                             &
     &          /,11x,'cg_eps   = ',1p,e15.8,0p,                        &
     &            3x,'(',i3.3,', ',i3.3,')',                            &
     &          /,11x,'cg_tol   = ',1p,e15.8,0p,                        &
     &            3x,'(',i3.3,', ',i3.3,')',                            &
     &          /,11x,'cg_gamma = ',1p,e15.8,0p,                        &
     &            3x,'(',i3.3,', ',i3.3,')',                            &
     &          /,11x,'cg_sigma = ',1p,e15.8,0p,                        &
     &            3x,'(',i3.3,', ',i3.3,')',                            &
     &          /,11x,'cg_rnorm = ',1p,e15.8,0p,                        &
     &            3x,'(',i3.3,', ',i3.3,')',/,                          &
     &          /,11x,'Min cg_x = ',1p,e15.8,0p,2x,                     &
     &                'Max cg_x = ',1p,e15.8,0p,                        &
     &          /,11x,'Min cg_z = ',1p,e15.8,0p,2x,                     &
     &                'Max cg_z = ',1p,e15.8,0p,                        &
     &          /,11x,'Min cg_p = ',1p,e15.8,0p,2x,                     &
     &                'Max cg_p = ',1p,e15.8,0p,                        &
     &          /,11x,'Min cg_v = ',1p,e15.8,0p,2x,                     &
     &                'Max cg_v = ',1p,e15.8,0p,                        &
     &          /,11x,'Min cg_s = ',1p,e15.8,0p,2x,                     &
     &                'Max cg_s = ',1p,e15.8,0p,                        &
     &          /,11x,'Min cg_r = ',1p,e15.8,0p,2x,                     &
     &                'Max cg_r = ',1p,e15.8,0p)
      END IF

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
        status=nf_inq_varid(ncMODid(ng),'outer',varid)
        status=nf_put_var1_int(ncMODid(ng),varid,1,outer)
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,10) 'outer', TRIM(MODname(ng))
          exit_flag=3
          ioerror=status
          RETURN
        END IF

        status=nf_inq_varid(ncMODid(ng),'inner',varid)
        status=nf_put_var1_int(ncMODid(ng),varid,1,inner)
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,10) 'inner', TRIM(MODname(ng))
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
!  Write out norms.
!
        status=nf_inq_varid(ncMODid(ng),'cg_gamma',varid)
        status=nf_put_var1_TYPE(ncMODid(ng),varid,1,cg_gammam1)
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,10) 'cg_gamma', TRIM(MODname(ng))
          exit_flag=3
          ioerror=status
          RETURN
        END IF

        status=nf_inq_varid(ncMODid(ng),'cg_sigma',varid)
        status=nf_put_var1_TYPE(ncMODid(ng),varid,1,cg_sigmam1)
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,10) 'cg_sigma', TRIM(MODname(ng))
          exit_flag=3
          ioerror=status
          RETURN
        END IF

        status=nf_inq_varid(ncMODid(ng),'cg_rnorm',varid)
        status=nf_put_var1_TYPE(ncMODid(ng),varid,1,cg_rnorm)
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,10) 'cg_rnorm', TRIM(MODname(ng))
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
!  Write out representer coefficients.
!
        start(1)=1
        total(1)=Ndatum(ng)
        status=nf_put_vara_TYPE(ncMODid(ng), modVid(idRepC,ng),         &
     &                          start, total, ADmodVal)
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,10) TRIM(Vname(1,idRepC)), TRIM(MODname(ng))
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
!  Write out conjugate gradient vectors.
!
        start(1)=1
        total(1)=Ndatum(ng)
        status=nf_put_vara_TYPE(ncMODid(ng), modVid(idCG_P,ng),         &
     &                          start, total, cg_p)
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,10) TRIM(Vname(1,idCG_P)), TRIM(MODname(ng))
          exit_flag=3
          ioerror=status
          RETURN
        END IF

        status=nf_put_vara_TYPE(ncMODid(ng), modVid(idCG_R,ng),         &
     &                        start, total, cg_r)
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,10) TRIM(Vname(1,idCG_R)), TRIM(MODname(ng))
          exit_flag=3
          ioerror=status
          RETURN
        END IF

        status=nf_put_vara_TYPE(ncMODid(ng), modVid(idCG_S,ng),         &
     &                          start, total, cg_s)
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,10) TRIM(Vname(1,idCG_S)), TRIM(MODname(ng))
          exit_flag=3
          ioerror=status
          RETURN
        END IF

        status=nf_put_vara_TYPE(ncMODid(ng), modVid(idCG_V,ng),         &
     &                          start, total, cg_v)
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,10) TRIM(Vname(1,idCG_V)), TRIM(MODname(ng))
          exit_flag=3
          ioerror=status
          RETURN
        END IF

        status=nf_put_vara_TYPE(ncMODid(ng), modVid(idCG_X,ng),         &
     &                          start, total, cg_x)
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,10) TRIM(Vname(1,idCG_X)), TRIM(MODname(ng))
          exit_flag=3
          ioerror=status
          RETURN
        END IF

        status=nf_put_vara_TYPE(ncMODid(ng), modVid(idCG_Z,ng),         &
     &                          start, total, cg_z)
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,10) TRIM(Vname(1,idCG_Z)), TRIM(MODname(ng))
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Synchronize observations NetCDF file to disk.
!-----------------------------------------------------------------------
!
      IF (OutThread) THEN
        status=nf_sync(ncMODid(ng))
        IF (status.ne.nf_noerr) THEN
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
