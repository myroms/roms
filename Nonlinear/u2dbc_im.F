#include "cppdefs.h"
      MODULE u2dbc_mod
!
!========================================== Alexander F. Shchepetkin ===
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This subroutine sets lateral boundary conditions for vertically     !
!  integrated U-velocity.                                              !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC u2dbc_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE u2dbc (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_ocean
      USE mod_stepping
!
      integer, intent(in) :: ng, tile

#include "tile.h"
!
      CALL u2dbc_tile (ng, Istr, Iend, Jstr, Jend,                      &
     &                 LBi, UBi, LBj, UBj,                              &
     &                 krhs(ng), kstp(ng), knew(ng),                    &
     &                 OCEAN(ng) % zeta,                                &
     &                 OCEAN(ng) % ubar)
      RETURN
      END SUBROUTINE u2dbc
!
!***********************************************************************
      SUBROUTINE u2dbc_tile (ng, Istr, Iend, Jstr, Jend,                &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       krhs, kstp, knew,                          &
     &                       zeta, ubar)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
      USE mod_grid
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: krhs, kstp, knew

#ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: zeta(LBi:,LBj:,:)
      real(r8), intent(inout) :: ubar(LBi:,LBj:,:)
#else
      real(r8), intent(in) :: zeta(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ubar(LBi:UBi,LBj:UBj,3)
#endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j, know

      real(r8), parameter :: eps = 1.0E-20_r8

      real(r8) :: Ce, Cx, cff, dt2d, dUde, dUdt, dUdx, tau

      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: grad

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set time-indices
!-----------------------------------------------------------------------
!
      IF (FIRST_2D_STEP) THEN
        know=krhs
        dt2d=dtfast(ng)
      ELSE IF (PREDICTOR_2D_STEP(ng)) THEN
        know=krhs
        dt2d=2.0_r8*dtfast(ng)
      ELSE
        know=kstp
        dt2d=dtfast(ng)
      END IF

#ifndef EW_PERIODIC
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the western edge.
!-----------------------------------------------------------------------
!
      IF (WESTERN_EDGE) THEN

# if defined WEST_M2RADIATION
!
!  Western edge, implicit upstream radiation condition.
!
        DO j=Jstr,Jend+1
          grad(1,j)=ubar(1,j  ,know)-                                   &
     &              ubar(1,j-1,know)
          grad(2,j)=ubar(2,j  ,know)-                                   &
     &              ubar(2,j-1,know)
        END DO
        DO j=Jstr,Jend
          dUdt=ubar(2,j,know)-ubar(2,j,knew)
          dUdx=ubar(2,j,knew)-ubar(3,j,knew)
#  ifdef WEST_M2NUDGING
          tau=M2obc_out(ng,iwest)
          IF ((dUdt*dUdx).lt.0.0_r8) tau=M2obc_in(ng,iwest)
          tau=tau*dt2d
#  endif
          IF ((dUdt*dUdx).lt.0.0_r8) dUdt=0.0_r8
          IF ((dUdt*(grad(2,j)+grad(2,j+1))).gt.0.0_r8) THEN
            dUde=grad(2,j  )
          ELSE
            dUde=grad(2,j+1)
          END IF
          cff=MAX(dUdx*dUdx+dUde*dUde,eps)
          Cx=dUdt*dUdx
#  ifdef RADIATION_2D
          Ce=MIN(cff,MAX(dUdt*dUde,-cff))
#  else
          Ce=0.0_r8
#  endif
          ubar(1,j,knew)=(cff*ubar(1,j,know)+                           &
     &                    Cx*ubar(2,j,knew)-                            &
     &                    MAX(Ce,0.0_r8)*grad(1,j  )-                   &
     &                    MIN(Ce,0.0_r8)*grad(1,j+1))/                  &
     &                   (cff+Cx)
#  ifdef WEST_M2NUDGING
         ubar(1,j,knew)=ubar(1,j,knew)+                                 &
     &                  tau*(BOUNDARY(ng)%ubar_west(j)-                 &
     &                       ubar(1,j,know))
#  endif
#  ifdef MASKING
          ubar(1,j,knew)=ubar(1,j,knew)*                                &
     &                   GRID(ng)%umask(1,j)
#  endif
        END DO
# elif defined WEST_M2FLATHER
!
!  Western edge, Flather boundary condition.
!
        DO j=Jstr,Jend
          ubar(1,j,knew)=BOUNDARY(ng)%ubar_west(j)-                     &
     &                   SQRT(g/(0.5_r8*(GRID(ng)%h(0,j)+               &
     &                                   GRID(ng)%h(1,j))))*            &
     &                   (0.5_r8*(zeta(0,j,know)+                       &
     &                            zeta(1,j,know))-                      &
     &                    BOUNDARY(ng)%zeta_west(j))
#  ifdef MASKING
          ubar(1,j,knew)=ubar(1,j,knew)*                                &
     &                   GRID(ng)%umask(1,j)
#  endif
        END DO
# elif defined WEST_M2CLAMPED
!
!  Western edge, clamped boundary condition.
!
        DO j=Jstr,Jend
          ubar(1,j,knew)=BOUNDARY(ng)%ubar_west(j)
#  ifdef MASKING
          ubar(1,j,knew)=ubar(1,j,knew)*                                &
     &                   GRID(ng)%umask(1,j)
#  endif
        END DO
# elif defined WEST_M2GRADIENT
!
!  Western edge, gradient boundary condition.
!
        DO j=Jstr,Jend
          ubar(1,j,knew)=ubar(2,j,knew)
#  ifdef MASKING
          ubar(1,j,knew)=ubar(1,j,knew)*                                &
     &                   GRID(ng)%umask(1,j)
#  endif
        END DO
# else
!
!  Western edge, closed boundary condition.
!
        DO j=Jstr,Jend
          ubar(1,j,knew)=0.0_r8
        END DO
# endif
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the eastern edge.
!-----------------------------------------------------------------------
!
      IF (EASTERN_EDGE) THEN

# if defined EAST_M2RADIATION
!
!  Eastern edge, implicit upstream radiation condition.
!
        DO j=Jstr,Jend+1
          grad(Lm(ng),j)=ubar(Lm(ng),j  ,know)-                         &
     &                   ubar(Lm(ng),j-1,know)
          grad(Lm(ng)+1,j)=ubar(Lm(ng)+1,j  ,know)-                     &
     &                     ubar(Lm(ng)+1,j-1,know)
        END DO
        DO j=Jstr,Jend
          dUdt=ubar(Lm(ng),j,know)-ubar(Lm(ng),j,knew)
          dUdx=ubar(Lm(ng),j,knew)-ubar(Lm(ng)-1,j,knew)
#  ifdef EAST_M2NUDGING
          tau=M2obc_out(ng,ieast)
          IF ((dUdt*dUdx).lt.0.0_r8) tau=M2obc_in(ng,ieast)
          tau=tau*dt2d
#  endif
          IF ((dUdt*dUdx).lt.0.0_r8) dUdt=0.0_r8
          IF ((dUdt*(grad(Lm(ng),j)+grad(Lm(ng),j+1))).gt.0.0_r8) THEN
            dUde=grad(Lm(ng),j)
          ELSE
            dUde=grad(Lm(ng),j+1)
          END IF
          cff=MAX(dUdx*dUdx+dUde*dUde,eps)
          Cx=dUdt*dUdx
#  ifdef RADIATION_2D
          Ce=MIN(cff,MAX(dUdt*dUde,-cff))
#  else
          Ce=0.0_r8
#  endif
          ubar(Lm(ng)+1,j,knew)=(cff*ubar(Lm(ng)+1,j,know)+             &
     &                           Cx*ubar(Lm(ng),j,knew)-                &
     &                           MAX(Ce,0.0_r8)*grad(Lm(ng)+1,j  )-     &
     &                           MIN(Ce,0.0_r8)*grad(Lm(ng)+1,j+1))/    &
     &                          (cff+Cx)
#  ifdef EAST_M2NUDGING
          ubar(Lm(ng)+1,j,knew)=ubar(Lm(ng)+1,j,knew)+                  &
     &                          tau*(BOUNDARY(ng)%ubar_east(j)-         &
     &                               ubar(Lm(ng)+1,j,know))
#  endif
#  ifdef MASKING
          ubar(Lm(ng)+1,j,knew)=ubar(Lm(ng)+1,j,knew)*                  &
     &                          GRID(ng)%umask(Lm(ng)+1,j)
#  endif
        END DO
# elif defined EAST_M2FLATHER
!
!  Eastern edge, Flather boundary condition.
!
        DO j=Jstr,Jend
          ubar(Lm(ng)+1,j,knew)=(BOUNDARY(ng)%ubar_east(j)+             &
     &                           SQRT(g/(0.5_r8*                        &
     &                                   (GRID(ng)%h(Lm(ng)  ,j)+       &
     &                                    GRID(ng)%h(Lm(ng)+1,j))))*    &
     &                           (0.5_r8*(zeta(Lm(ng)  ,j,know)+        &
     &                                    zeta(Lm(ng)+1,j,know))-       &
     &                            BOUNDARY(ng)%zeta_east(j)))
#  ifdef MASKING
          ubar(Lm(ng)+1,j,knew)=ubar(Lm(ng)+1,j,knew)*                  &
     &                          GRID(ng)%umask(Lm(ng)+1,j)
#  endif
        END DO
# elif defined EAST_M2CLAMPED
!
!  Eastern edge, clamped boundary condition.
!
        DO j=Jstr,Jend
          ubar(Lm(ng)+1,j,knew)=BOUNDARY(ng)%ubar_east(j)
#  ifdef MASKING
          ubar(Lm(ng)+1,j,knew)=ubar(Lm(ng)+1,j,knew)*                  &
     &                          GRID(ng)%umask(Lm(ng)+1,j)
#  endif
        END DO
# elif defined EAST_M2GRADIENT
!
!  Eastern edge, gradient boundary condition.
!
        DO j=Jstr,Jend
          ubar(Lm(ng)+1,j,knew)=ubar(Lm(ng),j,knew)
#  ifdef MASKING
          ubar(Lm(ng)+1,j,knew)=ubar(Lm(ng)+1,j,knew)*                  &
     &                          GRID(ng)%umask(Lm(ng)+1,j)
#  endif
        END DO
# else
!
!  Eastern edge, closed boundary condition.
!
        DO j=Jstr,Jend
          ubar(Lm(ng)+1,j,knew)=0.0_r8
        END DO
# endif
      END IF
#endif
#ifndef NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the southern edge.
!-----------------------------------------------------------------------
!
      IF (SOUTHERN_EDGE) THEN

# if defined SOUTH_M2RADIATION
!
!  Southern edge, implicit upstream radiation condition.
!
        DO i=IstrU-1,Iend
          grad(i,0)=ubar(i+1,0,know)-                                   &
     &              ubar(i  ,0,know)
          grad(i,1)=ubar(i+1,1,know)-                                   &
     &              ubar(i  ,1,know)
        END DO
        DO i=IstrU,Iend
          dUdt=ubar(i,1,know)-ubar(i,1,knew)
          dUde=ubar(i,1,knew)-ubar(i,2,knew)
#  ifdef SOUTH_M2NUDGING
          tau=M2obc_out(ng,isouth)
          IF ((dUdt*dUde).lt.0.0_r8) tau=M2obc_in(ng,isouth)
          tau=tau*dt2d
#  endif
          IF ((dUdt*dUde).lt.0.0_r8) dUdt=0.0_r8
          IF ((dUdt*(grad(i-1,1)+grad(i,1))).gt.0.0_r8) THEN
            dUdx=grad(i-1,1)
          ELSE
            dUdx=grad(i  ,1)
          END IF
          cff=MAX(dUdx*dUdx+dUde*dUde,eps)
#  ifdef RADIATION_2D
          Cx=MIN(cff,MAX(dUdt*dUdx,-cff))
#  else
          Cx=0.0_r8
#  endif
          Ce=dUdt*dUde
          ubar(i,0,knew)=(cff*ubar(i,0,know)+                           &
     &                    Ce*ubar(i,1,knew)-                            &
     &                    MAX(Cx,0.0_r8)*grad(i-1,0)-                   &
     &                    MIN(Cx,0.0_r8)*grad(i  ,0))/                  &
     &                   (cff+Ce)
#  ifdef SOUTH_M2NUDGING
          ubar(i,0,knew)=ubar(i,0,knew)+                                &
     &                   tau*(BOUNDARY(ng)%ubar_south(i)-               &
     &                        ubar(i,0,know))
#  endif
#  ifdef MASKING
          ubar(i,0,knew)=ubar(i,0,knew)*                                &
     &                   GRID(ng)%umask(i,0)
#  endif
        END DO
# elif defined SOUTH_M2FLATHER
!
!  Southern edge, Chapman boundary condition.
!
        DO i=IstrU,Iend
          Ce=dt2d*0.5_r8*(GRID(ng)%pn(i-1,1)+                           &
     &                    GRID(ng)%pn(i  ,1))*                          &
     &       SQRT(g*0.5_r8*(GRID(ng)%h(i-1,1)+                          &
     &                      GRID(ng)%h(i  ,1)))
          ubar(i,0,knew)=(ubar(i,0,know)+                               &
     &                    Ce*ubar(i,1,knew))/                           &
     &                   (1.0_r8+Ce)
#  ifdef MASKING
          ubar(i,0,knew)=ubar(i,0,knew)*                                &
     &                   GRID(ng)%umask(i,0)
#  endif
        END DO
# elif defined SOUTH_M2CLAMPED
!
!  Southern edge, clamped boundary condition.
!
        DO i=IstrU,Iend
          ubar(i,0,knew)=BOUNDARY(ng)%ubar_south(i)
#  ifdef MASKING
          ubar(i,0,knew)=ubar(i,0,knew)*                                &
     &                   GRID(ng)%umask(i,0)
#  endif
        END DO
# elif defined SOUTH_M2GRADIENT
!
!  Southern edge, gradient boundary condition.
!
        DO i=IstrU,Iend
          ubar(i,0,knew)=ubar(i,1,knew)
#  ifdef MASKING
          ubar(i,0,knew)=ubar(i,0,knew)*                                &
     &                   GRID(ng)%umask(i,0)
#  endif
        END DO
# else
!
!  Southern edge, closed boundary condition: free slip (gamma2=1)  or
!                                            no   slip (gamma2=-1).
!
#  ifdef EW_PERIODIC
#   define I_RANGE IstrU,Iend
#  else
#   define I_RANGE Istr,IendR
#  endif
        DO i=I_RANGE
          ubar(i,0,knew)=gamma2*ubar(i,1,knew)
#  ifdef MASKING
          ubar(i,0,knew)=ubar(i,0,knew)*                                &
     &                   GRID(ng)%umask(i,0)
#  endif
        END DO
#  undef I_RANGE
# endif
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the northern edge.
!-----------------------------------------------------------------------
!
      IF (NORTHERN_EDGE) THEN

# if defined NORTH_M2RADIATION
!
!  Northern edge, implicit upstream radiation condition.
!
        DO i=IstrU-1,Iend
          grad(i,Mm(ng))=ubar(i+1,Mm(ng),know)-                         &
     &                   ubar(i  ,Mm(ng),know)
          grad(i,Mm(ng)+1)=ubar(i+1,Mm(ng)+1,know)-                     &
     &                     ubar(i  ,Mm(ng)+1,know)
        END DO
        DO i=IstrU,Iend
          dUdt=ubar(i,Mm(ng),know)-ubar(i,Mm(ng),knew)
          dUde=ubar(i,Mm(ng),knew)-ubar(i,Mm(ng)-1,knew)
#  ifdef NORTH_M2NUDGING
          tau=M2obc_out(ng,inorth)
          IF ((dUdt*dUde).lt.0.0_r8) tau=M2obc_in(ng,inorth)
          tau=tau*dt2d
#  endif
          IF ((dUdt*dUde).lt.0.0_r8) dUdt=0.0_r8
          IF ((dUdt*(grad(i-1,Mm(ng))+grad(i,Mm(ng)))).gt.0.0_r8) THEN
            dUdx=grad(i-1,Mm(ng))
          ELSE
            dUdx=grad(i  ,Mm(ng))
          END IF
          cff=MAX(dUdx*dUdx+dUde*dUde,eps)
#  ifdef RADIATION_2D
          Cx=MIN(cff,MAX(dUdt*dUdx,-cff))
#  else
          Cx=0.0_r8
#  endif
          Ce=dUdt*dUde
          ubar(i,Mm(ng)+1,knew)=(cff*ubar(i,Mm(ng)+1,know)+             &
     &                           Ce*ubar(i,Mm(ng),knew)-                &
     &                           MAX(Cx,0.0_r8)*grad(i-1,Mm(ng)+1)-     &
     &                           MIN(Cx,0.0_r8)*grad(i  ,Mm(ng)+1))/    &
     &                          (cff+Ce)
#  ifdef NORTH_M2NUDGING
          ubar(i,Mm(ng)+1,knew)=ubar(i,Mm(ng)+1,knew)+                  &
     &                          tau*(BOUNDARY(ng)%ubar_north(i)-        &
     &                               ubar(i,Mm(ng)+1,know))
#  endif
#  ifdef MASKING
          ubar(i,Mm(ng)+1,knew)=ubar(i,Mm(ng)+1,knew)*                  &
     &                          GRID(ng)%umask(i,Mm(ng)+1)
#  endif
        END DO
# elif defined NORTH_M2FLATHER
!
!  Northern edge, Chapman boundary condition.
!
        DO i=IstrU,Iend
          Ce=dt2d*0.5_r8*(GRID(ng)%pn(i-1,Mm(ng))+                      &
     &                    GRID(ng)%pn(i  ,Mm(ng)))*                     &
     &       SQRT(g*0.5_r8*(GRID(ng)%h(i-1,Mm(ng))+                     &
     &                      GRID(ng)%h(i  ,Mm(ng))))
          ubar(i,Mm(ng)+1,knew)=(ubar(i,Mm(ng)+1,know)+                 &
     &                           Ce*ubar(i,Mm(ng),knew))/               &
     &                          (1.0_r8+Ce)
#  ifdef MASKING
          ubar(i,Mm(ng)+1,knew)=ubar(i,Mm(ng)+1,knew)*                  &
     &                          GRID(ng)%umask(i,Mm(ng)+1)
#  endif
        END DO
# elif defined NORTH_M2CLAMPED
!
!  Northern edge, clamped boundary condition.
!
        DO i=IstrU,Iend
          ubar(i,Mm(ng)+1,knew)=BOUNDARY(ng)%ubar_north(i)
#  ifdef MASKING
          ubar(i,Mm(ng)+1,knew)=ubar(i,Mm(ng)+1,knew)*                  &
     &                          GRID(ng)%umask(i,Mm(ng)+1)
#  endif
        END DO
# elif defined NORTH_M2GRADIENT
!
!  Northern edge, gradient boundary condition.
!
        DO i=IstrU,Iend
          ubar(i,Mm(ng)+1,knew)=ubar(i,Mm(ng),knew)
#  ifdef MASKING
          ubar(i,Mm(ng)+1,knew)=ubar(i,Mm(ng)+1,knew)*                  &
     &                          GRID(ng)%umask(i,Mm(ng)+1)
#  endif
        END DO
# else
!
!  Northern edge, closed boundary condition: free slip (gamma2=1)  or
!                                            no   slip (gamma2=-1).
!
#  ifdef EW_PERIODIC
#   define I_RANGE IstrU,Iend
#  else
#   define I_RANGE Istr,IendR
#  endif
        DO i=I_RANGE
          ubar(i,Mm(ng)+1,knew)=gamma2*ubar(i,Mm(ng),knew)
#  ifdef MASKING
          ubar(i,Mm(ng)+1,knew)=ubar(i,Mm(ng)+1,knew)*                  &
     &                          GRID(ng)%umask(i,Mm(ng)+1)
#  endif
        END DO
#  undef I_RANGE
# endif
      END IF
#endif
#if !defined EW_PERIODIC && !defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (SOUTHERN_EDGE.and.WESTERN_EDGE) THEN
        ubar(1,0,knew)=0.5_r8*(ubar(2,0,knew)+                          &
     &                         ubar(1,1,knew))
      END IF
      IF (SOUTHERN_EDGE.and.EASTERN_EDGE) THEN
        ubar(Lm(ng)+1,0,knew)=0.5_r8*(ubar(Lm(ng)  ,0,knew)+            &
     &                                ubar(Lm(ng)+1,1,knew))
      END IF
      IF (NORTHERN_EDGE.and.WESTERN_EDGE) THEN
        ubar(1,Mm(ng)+1,knew)=0.5_r8*(ubar(2,Mm(ng)+1,knew)+            &
     &                                ubar(1,Mm(ng)  ,knew))
      END IF
      IF (NORTHERN_EDGE.and.EASTERN_EDGE) THEN
        ubar(Lm(ng)+1,Mm(ng)+1,knew)=0.5_r8*                            &
     &                               (ubar(Lm(ng)  ,Mm(ng)+1,knew)+     &
     &                                ubar(Lm(ng)+1,Mm(ng)  ,knew))
      END IF
#endif
      RETURN
      END SUBROUTINE u2dbc_tile
      END MODULE u2dbc_mod
