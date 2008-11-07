#include "cppdefs.h"
#if defined FOUR_DVAR && defined OBSERVATIONS && \
    !(defined WEAK_CONSTRAINT || defined IOM)
      SUBROUTINE obs_cost (ng, model)
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  If minimization first pass, this routine computes the observation   !
!  cost function (Jo) as the misfit (squared difference) between the   !
!  model and observations.                                             !
!                                                                      !
!  If conventional strong contraint 4DVAR:                             !
!                                                                      !
!         Jo = 1/2 transpose(H X - Xo) * O^(-1) * (H X - Xo)           !
!                                                                      !
!  or if incremental strong contraint 4DVAR:                           !
!                                                                      !
!         Jo = 1/2 transpose(H deltaX - d) * O^(-1) * (H deltaX - d)   !
!                                                                      !
!  where                                                               !
!                                                                      !
!          d = Xo - H Xb                                               !
!                                                                      !
!         d  : innovation vector                                       !
!         H  : observation operator (linearized if incremental)        !
!       H Xb : background at observation points previous forecast)     !
!         Xo : observations vector                                     !
!       H X  : nonlinear model at observation points                   !
!  H deltaX  : increment at observation point                          !
!         O  : observations error covariance                           !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
!
!  Local variable declarations.
!
      integer :: NSUB, iobs, ivar

      real(r8) ::  cff

      real(r8), dimension(0:NstateVar(ng)) :: my_ObsCost
!
!-----------------------------------------------------------------------
!  Compute observation misfit cost function (ObsCost).
!-----------------------------------------------------------------------
!
!  Compute tangent linear model cost function.
!
      IF (model.eq.iTLM) THEN
        DO ivar=0,NstateVar(ng)
          my_ObsCost(ivar)=0.0_r8
        END DO
        DO iobs=1,Nobs(ng)
          ivar=ObsType(iobs)
          cff=0.5_r8*ObsScale(iobs)*ObsErr(iobs)*                       &
     &        (NLmodVal(iobs)+TLmodVal(iobs)-ObsVal(iobs))**2
          my_ObsCost(0)=my_ObsCost(0)+cff
          my_ObsCost(ivar)=my_ObsCost(ivar)+cff
        END DO
!
!  Compute nonlinear model cost function.
!
      ELSE IF (model.eq.iNLM) THEN
        DO ivar=0,NstateVar(ng)
          my_ObsCost(ivar)=0.0_r8
        END DO
        DO iobs=1,Nobs(ng)
          ivar=ObsType(iobs)
          cff=0.5_r8*ObsScale(iobs)*ObsErr(iobs)*                       &
     &        (NLmodVal(iobs)-ObsVal(iobs))**2
          my_ObsCost(0)=my_ObsCost(0)+cff
          my_ObsCost(ivar)=my_ObsCost(ivar)+cff
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Load global values.  Notice that there is not need for a global
!  reduction here since all the threads have the same copy of all
!  the vectors used.
!-----------------------------------------------------------------------
!
      NSUB=NtileX(ng)*NtileE(ng)
!$OMP CRITICAL (COST_FUN)
      tile_count=tile_count+1
      IF (tile_count.eq.NSUB) THEN
        tile_count=0
        IF (model.eq.iTLM) THEN
          DO ivar=0,NstateVar(ng)
            FOURDVAR(ng)%ObsCost(ivar)=FOURDVAR(ng)%ObsCost(ivar)+      &
     &                                 my_ObsCost(ivar)
          END DO
        ELSE IF (model.eq.iNLM) THEN
          DO ivar=0,NstateVar(ng)
            FOURDVAR(ng)%NLobsCost(ivar)=FOURDVAR(ng)%NLobsCost(ivar)+  &
     &                                   my_ObsCost(ivar)
          END DO
        END IF
      END IF
!$OMP END CRITICAL (COST_FUN)
!
!  If start of minimization, set cost function scales used to report
!  normalized values.
!
      IF ((Nrun.eq.1).and.(model.eq.iTLM)) THEN
        DO ivar=0,NstateVar(ng)
          FOURDVAR(ng)%CostNorm(ivar)=FOURDVAR(ng)%ObsCost(ivar)
        END DO
      END IF
!
!  Save initial inner loop cost function.
!
      IF ((inner.eq.0).and.(model.eq.iTLM)) THEN
        FOURDVAR(ng)%Cost0(outer)=FOURDVAR(ng)%ObsCost(0)
      END IF

      RETURN
      END SUBROUTINE obs_cost
#else
      SUBROUTINE obs_cost
      RETURN
      END SUBROUTINE obs_cost
#endif
