#include "cppdefs.h"
#ifdef SOLVE3D
      SUBROUTINE set_scoord (ng)
!
!========================================== Alexander F. Shchepetkin ===
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
! This routine sets and initializes relevant variables associated      !
! with the vertical transformation to nondimensional S-coordinate.     !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_grid
      USE mod_iounits
      USE mod_scalars

      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      integer :: k
      real(r8) :: cff, cff1, cff2, cff3
!
!-----------------------------------------------------------------------
!  Define S-Curves in domain [-1 < sc < 0] at vertical W- and
!  RHO-points.
!-----------------------------------------------------------------------
!
      hc(ng)=MIN(hmin(ng),Tcline(ng))
      IF (theta_s(ng).ne.0.0_r8) THEN
        cff1=1.0_r8/SINH(theta_s(ng))
        cff2=0.5_r8/TANH(0.5_r8*theta_s(ng))
      END IF
      SCALARS(ng)%sc_w(0)=-1.0_r8
      SCALARS(ng)%Cs_w(0)=-1.0_r8
      cff=1.0_r8/REAL(N(ng),r8)
      DO k=1,N(ng)
        SCALARS(ng)%sc_w(k)=cff*REAL(k-N(ng),r8)
        SCALARS(ng)%sc_r(k)=cff*(REAL(k-N(ng),r8)-0.5_r8)
        IF (theta_s(ng).ne.0.0_r8) THEN
          SCALARS(ng)%Cs_w(k)=(1.0_r8-theta_b(ng))*                     &
     &                        cff1*SINH(theta_s(ng)*                    &
     &                                  SCALARS(ng)%sc_w(k))+           &
     &                        theta_b(ng)*                              &
     &                        (cff2*TANH(theta_s(ng)*                   &
     &                                   (SCALARS(ng)%sc_w(k)+0.5_r8))- &
     &                         0.5_r8)
          SCALARS(ng)%Cs_r(k)=(1.0_r8-theta_b(ng))*                     &
     &                        cff1*SINH(theta_s(ng)*                    &
     &                                  SCALARS(ng)%sc_r(k))+           &
     &                        theta_b(ng)*                              &
     &                        (cff2*TANH(theta_s(ng)*                   &
     &                                   (SCALARS(ng)%sc_r(k)+0.5_r8))- &
     &                         0.5_r8)
        ELSE
          SCALARS(ng)%Cs_w(k)=SCALARS(ng)%sc_w(k)
          SCALARS(ng)%Cs_r(k)=SCALARS(ng)%sc_r(k)
        END IF
      END DO
!
!-----------------------------------------------------------------------
!  Report information about vertical S-levels.
!-----------------------------------------------------------------------
!
      IF (Master) THEN
        WRITE (stdout,10)
  10    FORMAT (/,' Vertical S-coordinate System: ',/,/,                &
     &            ' level   S-coord     Cs-curve',10x,                  &
     &            'at_hmin  over_slope     at_hmax',/)
        DO k=N(ng),0,-1
          cff1=SCALARS(ng)%sc_w(k)*hc(ng)+                              &
     &         (hmin(ng)-hc(ng))*SCALARS(ng)%Cs_w(k)
          cff2=SCALARS(ng)%sc_w(k)*hc(ng)+                              &
     &         (0.5*(hmin(ng)+hmax(ng))-hc(ng))*SCALARS(ng)%Cs_w(k)
          cff3=SCALARS(ng)%sc_w(k)*hc(ng)+                              &
     &         (hmax(ng)-hc(ng))*SCALARS(ng)%Cs_w(k)
          WRITE (stdout,20) k, SCALARS(ng)%sc_w(k),                     &
     &                         SCALARS(ng)%Cs_w(k), cff1, cff2, cff3
  20      FORMAT (i6,2f12.7,4x,3f12.3)
        END DO
      END IF
#else
      SUBROUTINE set_scoord
#endif
      RETURN
      END SUBROUTINE set_scoord
