#include "cppdefs.h"
      MODULE v2dbc_mod
!
!========================================== Alexander F. Shchepetkin ===
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This subroutine sets lateral boundary conditions for vertically     !
!  integrated V-velocity.                                              !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC v2dbc_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE v2dbc (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_ocean
      USE mod_stepping
!
      integer, intent(in) :: ng, tile

#include "tile.h"
!
      CALL v2dbc_tile (ng, Istr, Iend, Jstr, Jend,                      &
     &                 LBi, UBi, LBj, UBj,                              &
     &                 krhs(ng), kstp(ng), knew(ng),                    &
     &                 OCEAN(ng) % zeta,                                &
     &                 OCEAN(ng) % vbar)
      RETURN
      END SUBROUTINE v2dbc
!
!***********************************************************************
      SUBROUTINE v2dbc_tile (ng, Istr, Iend, Jstr, Jend,                &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       krhs, kstp, knew,                          &
     &                       zeta, vbar)
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
      real(r8), intent(inout) :: vbar(LBi:,LBj:,:)
#else
      real(r8), intent(in) :: zeta(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: vbar(LBi:UBi,LBj:UBj,3)
#endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j, know

      real(r8), parameter :: eps = 1.0E-20_r8

      real(r8) :: Ce, Cx, cff, dt2d, dVde, dVdt, dVdx, tau

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
        DO j=JstrV-1,Jend
          grad(0,j)=vbar(0,j+1,know)-                                   &
     &              vbar(0,j  ,know)
          grad(1,j)=vbar(1,j+1,know)-                                   &
     &              vbar(1,j  ,know)
        END DO
        DO j=JstrV,Jend
          dVdt=vbar(1,j,know)-vbar(1,j,knew)
          dVdx=vbar(1,j,knew)-vbar(2,j,knew)
#  ifdef WEST_M2NUDGING
          tau=M2obc_out(ng,iwest)
          IF ((dVdt*dVdx).lt.0.0_r8) tau=M2obc_in(ng,iwest)
          tau=tau*dt2d
#  endif
          IF ((dVdt*dVdx).lt.0.0_r8) dVdt=0.0_r8
          IF ((dVdt*(grad(1,j-1)+grad(1,j))).gt.0.0_r8) THEN
            dVde=grad(1,j-1)
          ELSE
            dVde=grad(1,j  )
          END IF
          cff=MAX(dVdx*dVdx+dVde*dVde,eps)
          Cx=dVdt*dVdx
#  ifdef RADIATION_2D
          Ce=MIN(cff,MAX(dVdt*dVde,-cff))
#  else
          Ce=0.0_r8
#  endif
          vbar(0,j,knew)=(cff*vbar(0,j,know)+                           &
     &                    Cx*vbar(1,j,knew)-                            &
     &                    MAX(Ce,0.0_r8)*grad(0,j-1)-                   &
     &                    MIN(Ce,0.0_r8)*grad(0,j  ))/                  &
     &                   (cff+Cx)
#  ifdef WEST_M2NUDGING
          vbar(0,j,knew)=vbar(0,j,knew)+                                &
     &                   tau*(BOUNDARY(ng)%vbar_west(j)-                &
     &                        vbar(0,j,know))
#   endif
#  ifdef MASKING
          vbar(0,j,knew)=vbar(0,j,knew)*                                &
     &                   GRID(ng)%vmask(0,j)
#  endif
        END DO
# elif defined WEST_M2FLATHER
!
!  Western edge, Chapman boundary condition.
!
        DO j=JstrV,Jend
          Cx=dt2d*0.5_r8*(GRID(ng)%pm(1,j-1)+                           &
     &                    GRID(ng)%pm(1,j  ))*                          &
     &       SQRT(g*0.5_r8*(GRID(ng)%h(1,j-1)+                          &
     &                      GRID(ng)%h(1,j  )))
          vbar(0,j,knew)=(vbar(0,j,know)+                               &
     &                    Cx*vbar(1,j,knew))/                           &
     &                   (1.0_r8+Cx)
#  ifdef MASKING
          vbar(0,j,knew)=vbar(0,j,knew)*                                &
     &                   GRID(ng)%vmask(0,j)
#  endif
        END DO
# elif defined WEST_M2CLAMPED
!
!  Western edge, clamped boundary condition.
!
        DO j=JstrV,Jend
          vbar(0,j,knew)=BOUNDARY(ng)%vbar_west(j)
#  ifdef MASKING
          vbar(0,j,knew)=vbar(0,j,knew)*                                &
     &                   GRID(ng)%vmask(0,j)
#  endif
        END DO
# elif defined WEST_M2GRADIENT
!
!  Western edge, gradient boundary condition.
!
        DO j=JstrV,Jend
          vbar(0,j,knew)=vbar(1,j,knew)
#  ifdef MASKING
          vbar(0,j,knew)=vbar(0,j,knew)*                                &
     &                   GRID(ng)%vmask(0,j)
#  endif
        END DO
# else
!
!  Western edge, closed boundary condition: free slip (gamma2=1)  or
!                                           no   slip (gamma2=-1).
!
#  ifdef NS_PERIODIC
#   define J_RANGE JstrV,Jend
#  else
#   define J_RANGE Jstr,JendR
#  endif
        DO j=J_RANGE
          vbar(0,j,knew)=gamma2*vbar(1,j,knew)
#  ifdef MASKING
          vbar(0,j,knew)=vbar(0,j,knew)*                                &
     &                   GRID(ng)%vmask(0,j)
#  endif
        END DO
#  undef J_RANGE
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
        DO j=JstrV-1,Jend
          grad(Lm(ng),j)=vbar(Lm(ng),j+1,know)-                         &
     &                   vbar(Lm(ng),j  ,know)
          grad(Lm(ng)+1,j)=vbar(Lm(ng)+1,j+1,know)-                     &
     &                     vbar(Lm(ng)+1,j  ,know)
        END DO
        DO j=JstrV,Jend
          dVdt=vbar(Lm(ng),j,know)-vbar(Lm(ng),j,knew)
          dVdx=vbar(Lm(ng),j,knew)-vbar(Lm(ng)-1,j,knew)
#  ifdef EAST_M2NUDGING
          tau=M2obc_out(ng,ieast)
          IF ((dVdt*dVdx).lt.0.0_r8) tau=M2obc_in(ng,ieast)
          tau=tau*dt2d
#  endif
          IF ((dVdt*dVdx).lt.0.0_r8) dVdt=0.0_r8
          IF ((dVdt*(grad(Lm(ng),j-1)+grad(Lm(ng),j))).gt.0.0_r8) THEN
            dVde=grad(Lm(ng),j-1)
          ELSE
            dVde=grad(Lm(ng),j  )
          END IF
          cff=MAX(dVdx*dVdx+dVde*dVde,eps)
          Cx=dVdt*dVdx
#  ifdef RADIATION_2D
          Ce=MIN(cff,MAX(dVdt*dVde,-cff))
#  else
          Ce=0.0_r8
#  endif
          vbar(Lm(ng)+1,j,knew)=(cff*vbar(Lm(ng)+1,j,know)+             &
     &                           Cx*vbar(Lm(ng),j,knew)-                &
     &                           MAX(Ce,0.0_r8)*grad(Lm(ng)+1,j-1)-     &
     &                           MIN(Ce,0.0_r8)*grad(Lm(ng)+1,j  ))/    &
     &                          (cff+Cx)
#  ifdef EAST_M2NUDGING
          vbar(Lm(ng)+1,j,knew)=vbar(Lm(ng)+1,j,knew)+                  &
     &                          tau*(BOUNDARY(ng)%vbar_east(j)-         &
     &                               vbar(Lm(ng)+1,j,know))
#  endif
#  ifdef MASKING
          vbar(Lm(ng)+1,j,knew)=vbar(Lm(ng)+1,j,knew)*                  &
     &                          GRID(ng)%vmask(Lm(ng)+1,j)
#  endif
        END DO
# elif defined EAST_M2FLATHER
!
!  Eastern edge, Chapman boundary condition.
!
        DO j=JstrV,Jend
          Cx=dt2d*0.5_r8*(GRID(ng)%pm(Lm(ng),j-1)+                      &
     &                    GRID(ng)%pm(Lm(ng),j  ))*                     &
     &       SQRT(g*0.5_r8*(GRID(ng)%h(Lm(ng),j-1)+                     &
     &                      GRID(ng)%h(Lm(ng),j  )))
          vbar(Lm(ng)+1,j,knew)=(vbar(Lm(ng)+1,j,know)+                 &
     &                           Cx*vbar(Lm(ng),j,knew))/               &
     &                          (1.0_r8+Cx)
#  ifdef MASKING
          vbar(Lm(ng)+1,j,knew)=vbar(Lm(ng)+1,j,knew)*                  &
     &                          GRID(ng)%vmask(Lm(ng)+1,j)
#  endif
        END DO
# elif defined EAST_M2CLAMPED
!
!  Eastern edge, clamped boundary condition.
!
        DO j=JstrV,Jend
          vbar(Lm(ng)+1,j,knew)=BOUNDARY(ng)%vbar_east(j)
#  ifdef MASKING
          vbar(Lm(ng)+1,j,knew)=vbar(Lm(ng)+1,j,knew)*                  &
     &                          GRID(ng)%vmask(Lm(ng)+1,j)
#  endif
        END DO
# elif defined EAST_M2GRADIENT
!
!  Eastern edge, gradient boundary condition.
!
        DO j=JstrV,Jend
          vbar(Lm(ng)+1,j,knew)=vbar(Lm(ng),j,knew)
#  ifdef MASKING
          vbar(Lm(ng)+1,j,knew)=vbar(Lm(ng)+1,j,knew)*                        &
     &                       GRID(ng)%vmask(Lm(ng)+1,j)
#  endif
        END DO
# else
!
!  Eastern edge, closed boundary condition: free slip (gamma2=1)  or
!                                           no   slip (gamma2=-1).
!
#  ifdef NS_PERIODIC
#   define J_RANGE JstrV,Jend
#  else
#   define J_RANGE Jstr,JendR
#  endif
        DO j=J_RANGE
          vbar(Lm(ng)+1,j,knew)=gamma2*vbar(Lm(ng),j,knew)
#  ifdef MASKING
          vbar(Lm(ng)+1,j,knew)=vbar(Lm(ng)+1,j,knew)*                  &
     &                          GRID(ng)%vmask(Lm(ng)+1,j)
#  endif
        END DO
#  undef J_RANGE
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
        DO i=Istr,Iend+1
          grad(i,1)=vbar(i  ,1,know)-                                   &
     &              vbar(i-1,1,know)
          grad(i,2)=vbar(i  ,2,know)-                                   &
     &              vbar(i-1,2,know)
        END DO
        DO i=Istr,Iend
          dVdt=vbar(i,2,know)-vbar(i,2,knew)
          dVde=vbar(i,2,knew)-vbar(i,3,knew)
#  ifdef SOUTH_M2NUDGING
          tau=M2obc_out(ng,isouth)
          IF ((dVdt*dVde).lt.0.0_r8) tau=M2obc_in(ng,isouth)
          tau=tau*dt2d
#  endif
          IF ((dVdt*dVde).lt.0.0_r8) dVdt=0.0_r8
          IF ((dVdt*(grad(i,2)+grad(i+1,2))).gt.0.0_r8) THEN
            dVdx=grad(i  ,2)
          ELSE
            dVdx=grad(i+1,2)
          END IF
          cff=MAX(dVdx*dVdx+dVde*dVde,eps)
#  ifdef RADIATION_2D
          Cx=MIN(cff,MAX(dVdt*dVdx,-cff))
#  else
          Cx=0.0_r8
#  endif
          Ce=dVdt*dVde
          vbar(i,1,knew)=(cff*vbar(i,1,know)+                           &
     &                    Ce*vbar(i,2,knew)-                            &
     &                    MAX(Cx,0.0_r8)*grad(i  ,1)-                   &
     &                    MIN(Cx,0.0_r8)*grad(i+1,1))/                  &
     &                   (cff+Ce)
#  ifdef SOUTH_M2NUDGING
          vbar(i,1,knew)=vbar(i,1,knew)+                                &
     &                   tau*(BOUNDARY(ng)%vbar_south(i)-               &
     &                        vbar(i,1,know))
#  endif
#  ifdef MASKING
          vbar(i,1,knew)=vbar(i,1,knew)*                                &
     &                   GRID(ng)%vmask(i,1)
#  endif
        END DO
# elif defined SOUTH_M2FLATHER
!
!  Southern edge, Flather boundary condition.
!
        DO i=Istr,Iend
          vbar(i,1,knew)=BOUNDARY(ng)%vbar_south(i)-                    &
     &                   SQRT(g/(0.5_r8*(GRID(ng)%h(i,0)+               &
     &                                   GRID(ng)%h(i,1))))*            &
     &                   (0.5_r8*(zeta(i,0,know)+zeta(i,1,know))-       &
     &                    BOUNDARY(ng)%zeta_south(i))
#  ifdef MASKING
          vbar(i,1,knew)=vbar(i,1,knew)*                                &
     &                   GRID(ng)%vmask(i,1)
#  endif
        END DO
# elif defined SOUTH_M2CLAMPED
!
!  Southern edge, clamped boundary condition.
!
        DO i=Istr,Iend
          vbar(i,1,knew)=BOUNDARY(ng)%vbar_south(i)
#  ifdef MASKING
          vbar(i,1,knew)=vbar(i,1,knew)*                                &
     &                   GRID(ng)%vmask(i,1)
#  endif
        END DO
# elif defined SOUTH_M2GRADIENT
!
!  Southern edge, gradient boundary condition.
!
        DO i=Istr,Iend
          vbar(i,1,knew)=vbar(i,2,knew)
#  ifdef MASKING
          vbar(i,1,knew)=vbar(i,1,knew)*                                &
     &                   GRID(ng)%vmask(i,1)
#  endif
        END DO
# else
!
!  Southern edge, closed boundary condition.
!
        DO i=Istr,Iend
          vbar(i,1,knew)=0.0_r8
        END DO
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
        DO i=Istr,Iend+1
          grad(i,Mm(ng))=vbar(i  ,Mm(ng),know)-                         &
     &                   vbar(i-1,Mm(ng),know)
          grad(i,Mm(ng)+1)=vbar(i  ,Mm(ng)+1,know)-                     &
     &                     vbar(i-1,Mm(ng)+1,know)
        END DO
        DO i=Istr,Iend
          dVdt=vbar(i,Mm(ng),know)-vbar(i,Mm(ng),knew)
          dVde=vbar(i,Mm(ng),knew)-vbar(i,Mm(ng)-1,knew)
#  ifdef NORTH_M2NUDGING
          tau=M2obc_out(ng,inorth)
          IF ((dVdt*dVde).lt.0.0_r8) tau=M2obc_in(ng,inorth)
          tau=tau*dt2d
#  endif
          IF ((dVdt*dVde).lt.0.0_r8) dVdt=0.0_r8
          IF ((dVdt*(grad(i,Mm(ng))+grad(i,Mm(ng)+1))).gt.0.0_r8) THEN
            dVdx=grad(i  ,Mm(ng))
          ELSE
            dVdx=grad(i+1,Mm(ng))
          END IF
          cff=MAX(dVdx*dVdx+dVde*dVde,eps)
#  ifdef RADIATION_2D
          Cx=MIN(cff,MAX(dVdt*dVdx,-cff))
#  else
          Cx=0.0_r8
#  endif
          Ce=dVdt*dVde
          vbar(i,Mm(ng)+1,knew)=(cff*vbar(i,Mm(ng)+1,know)+             &
     &                           Ce*vbar(i,Mm(ng),knew)-                &
     &                           MAX(Cx,0.0_r8)*grad(i  ,Mm(ng)+1)-     &
     &                           MIN(Cx,0.0_r8)*grad(i+1,Mm(ng)+1))/    &
     &                          (cff+Ce)
#  ifdef NORTH_M2NUDGING
          vbar(i,Mm(ng)+1,knew)=vbar(i,Mm(ng)+1,knew)+                  &
     &                          tau*(BOUNDARY(ng)%vbar_north(i)-        &
     &                               vbar(i,Mm(ng)+1,know))
#  endif
#  ifdef MASKING
          vbar(i,Mm(ng)+1,knew)=vbar(i,Mm(ng)+1,knew)*                  &
     &                          GRID(ng)%vmask(i,Mm(ng)+1)
#  endif
        END DO
# elif defined NORTH_M2FLATHER
!
!  Northern edge, Flather boundary condition.
!
        DO i=Istr,Iend
          vbar(i,Mm(ng)+1,knew)=BOUNDARY(ng)%vbar_north(i)+             &
     &                          SQRT(g/(0.5_r8*                         &
     &                                  (GRID(ng)%h(i,Mm(ng))+          &
     &                                   GRID(ng)%h(i,Mm(ng)+1))))*     &
     &                          (0.5_r8*(zeta(i,Mm(ng),know)+           &
     &                                   zeta(i,Mm(ng)+1 ,know))-       &
     &                           BOUNDARY(ng)%zeta_north(i))
#  ifdef MASKING
          vbar(i,Mm(ng)+1,knew)=vbar(i,Mm(ng)+1,knew)*                  &
     &                          GRID(ng)%vmask(i,Mm(ng)+1)
#  endif
        END DO
# elif defined NORTH_M2CLAMPED
!
!  Northern edge, clamped boundary condition.
!
        DO i=Istr,Iend
          vbar(i,Mm(ng)+1,knew)=BOUNDARY(ng)%vbar_north(i)
#  ifdef MASKING
          vbar(i,Mm(ng)+1,knew)=vbar(i,Mm(ng)+1,knew)*                  &
     &                          GRID(ng)%vmask(i,Mm(ng)+1)
#  endif
        END DO
# elif defined NORTH_M2GRADIENT
!
!  Northern edge, gradient boundary condition.
!
        DO i=Istr,Iend
          vbar(i,Mm(ng)+1,knew)=vbar(i,Mm(ng),knew)
#  ifdef MASKING
          vbar(i,Mm(ng)+1,knew)=vbar(i,Mm(ng)+1,knew)*                  &
     &                          GRID(ng)%vmask(i,Mm(ng)+1)
#  endif
        END DO
# else
!
!  Northern edge, closed boundary condition.
!
        DO i=Istr,Iend
          vbar(i,Mm(ng)+1,knew)=0.0_r8
        END DO
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
        vbar(0,1,knew)=0.5_r8*(vbar(0,2,knew)+                          &
     &                         vbar(1,1,knew))
      END IF
      IF (SOUTHERN_EDGE.and.EASTERN_EDGE) THEN
        vbar(Lm(ng)+1,1,knew)=0.5_r8*(vbar(Lm(ng)  ,1,knew)+            &
     &                                vbar(Lm(ng)+1,2,knew))
      END IF
      IF (NORTHERN_EDGE.and.WESTERN_EDGE) THEN
        vbar(0,Mm(ng)+1,knew)=0.5_r8*(vbar(0,Mm(ng)  ,knew)+            &
     &                                vbar(1,Mm(ng)+1,knew))
      END IF
      IF (NORTHERN_EDGE.and.EASTERN_EDGE) THEN
        vbar(Lm(ng)+1,Mm(ng)+1,knew)=0.5_r8*                            &
     &                               (vbar(Lm(ng)+1,Mm(ng)  ,knew)+     &
     &                                vbar(Lm(ng)  ,Mm(ng)+1,knew))
      END IF
#endif
      RETURN
      END SUBROUTINE v2dbc_tile
      END MODULE v2dbc_mod
