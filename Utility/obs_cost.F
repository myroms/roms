#include "cppdefs.h"
#ifdef FOUR_DVAR
      SUBROUTINE obs_cost (ng)
!
!=================================================== Andrew M. Moore ===
!  Copyright (c) 2005 ROMS/TOMS Adjoint Group                          !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine computes the observation cost function (Jo) as the     !
!  misfit (squared difference) between model and observations:         !
!                                                                      !
!         Jo = transpose(X - Xo) * O^(-1) * (X -Xo)                    !
!                                                                      !
!  where                                                               !
!                                                                      !
!         Xo : observations                                            !
!         X  : model at observation points                             !
!         O  : observations error covariance (assigned weight)         !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      integer :: iobs, ivar

      real(r8) ::  cff
!
!-----------------------------------------------------------------------
!  Compute observation norm scale (ObsNorm).
!-----------------------------------------------------------------------
!
      CALL obs_norm (ng)
!
!-----------------------------------------------------------------------
!  Compute observation misfit cost function (CostFun).
!-----------------------------------------------------------------------
!
      IF (Ipass.eq.1) THEN
        DO iobs=1,Nobs(ng)
          ivar=ObsType(iobs)
# if defined S4DVAR
          cff=ObsNorm(iobs)*ObsErr(iobs)*                               &
     &        (NLmodVal(iobs,Ipass)-ObsVal(iobs))**2
# elif defined TLM_OBS
          cff=ObsNorm(iobs)*ObsErr(iobs)*                               &
     &        (NLmodVal(iobs,1)+TLmodVal(iobs,Ipass)-ObsVal(iobs))**2
# endif
          FOURDVAR(ng)%CostFun(0)=FOURDVAR(ng)%CostFun(0)+cff
          FOURDVAR(ng)%CostFun(ivar)=FOURDVAR(ng)%CostFun(ivar)+cff
        END DO
      END IF
# if defined S4DVAR || defined IS4DVAR
!
!-----------------------------------------------------------------------
!  Compute Conjugate gradient algorithm step size dot products,
!  StepSizeT, StepSizeB:
!
!               CGstepR=CGstepF*(StepSizeT/StepSizeB)
!-----------------------------------------------------------------------
!
      IF (Ipass.eq.2) THEN
        StepSizeT=0.0_r8
        StepSizeB=0.0_r8
        DO iobs=1,Nobs(ng)
          cff=ObsNorm(iobs)*ObsErr(iobs)
#  if defined S4DVAR
          StepSizeT=StepSizeT-                                          &
     &              cff*                                                &
     &              (NLmodVal(iobs,Ipass-1)-ObsVal(iobs))*              &
     &              (NLmodVal(iobs,Ipass)-NLmodVal(iobs,Ipass-1))
          StepSizeB=StepSizeB+                                          &
     &              cff*                                                &
     &              (NLmodVal(iobs,Ipass)-NLmodVal(iobs,Ipass-1))**2
#  elif defined IS4DVAR
          StepSizeT=StepSizeT-                                          &
     &              cff*                                                &
     &              (NLmodVal(iobs,1)+TLmodVal(iobs,Ipass-1)-           &
     &               ObsVal(iobs))*                                     &
     &              (TLmodVal(iobs,Ipass)-TLmodVal(iobs,Ipass-1))
          StepSizeB=StepSizeB+                                          &
     &              cff*                                                &
     &              (TLmodVal(iobs,Ipass)-TLmodVal(iobs,Ipass-1))**2
#  endif 
        END DO
      END IF
# endif
      RETURN
      END SUBROUTINE obs_cost
#else
      SUBROUTINE obs_cost
      RETURN
      END SUBROUTINE obs_cost
#endif
