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
      PUBLIC  :: v2dbc_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE v2dbc (ng, kout, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_ocean
      USE mod_stepping
!
      integer, intent(in) :: ng, kout, tile

#include "tile.h"
!
      CALL v2dbc_tile (ng, Istr, Iend, Jstr, Jend,                      &
     &                 LBi, UBi, LBj, UBj,                              &
     &                 krhs(ng), kstp(ng), kout,                        &
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
          grad(Istr-1,j)=vbar(Istr-1,j+1,know)-                         &
     &                   vbar(Istr-1,j  ,know)
          grad(Istr  ,j)=vbar(Istr  ,j+1,know)-                         &
     &                   vbar(Istr  ,j  ,know)
        END DO
        DO j=JstrV,Jend
          dVdt=vbar(Istr,j,know)-vbar(Istr  ,j,knew)
          dVdx=vbar(Istr,j,knew)-vbar(Istr+1,j,knew)
#  ifdef WEST_M2NUDGING
          tau=M2obc_out(ng,iwest)
          IF ((dVdt*dVdx).lt.0.0_r8) tau=M2obc_in(ng,iwest)
          tau=tau*dt2d
#  endif
          IF ((dVdt*dVdx).lt.0.0_r8) dVdt=0.0_r8
          IF ((dVdt*(grad(Istr,j-1)+grad(Istr,j))).gt.0.0_r8) THEN
            dVde=grad(Istr,j-1)
          ELSE
            dVde=grad(Istr,j  )
          END IF
          cff=MAX(dVdx*dVdx+dVde*dVde,eps)
          Cx=dVdt*dVdx
#  ifdef RADIATION_2D
          Ce=MIN(cff,MAX(dVdt*dVde,-cff))
#  else
          Ce=0.0_r8
#  endif
          vbar(Istr-1,j,knew)=(cff*vbar(Istr-1,j,know)+                 &
     &                         Cx*vbar(Istr,j,knew)-                    &
     &                         MAX(Ce,0.0_r8)*grad(Istr-1,j-1)-         &
     &                         MIN(Ce,0.0_r8)*grad(Istr-1,j  ))/        &
     &                        (cff+Cx)
#  ifdef WEST_M2NUDGING
          vbar(Istr-1,j,knew)=vbar(Istr-1,j,knew)+                      &
     &                        tau*(BOUNDARY(ng)%vbar_west(j)-           &
     &                        vbar(Istr-1,j,know))
#   endif
#  ifdef MASKING
          vbar(Istr-1,j,knew)=vbar(Istr-1,j,knew)*                      &
     &                        GRID(ng)%vmask(Istr-1,j)
#  endif
        END DO

# elif defined WEST_M2FLATHER
!
!  Western edge, Chapman boundary condition.
!
        DO j=JstrV,Jend
          Cx=dt2d*0.5_r8*(GRID(ng)%pm(Istr,j-1)+                        &
     &                    GRID(ng)%pm(Istr,j  ))*                       &
     &       SQRT(g*0.5_r8*(GRID(ng)%h(Istr,j-1)+                       &
     &                      GRID(ng)%h(Istr,j  )))
          vbar(Istr-1,j,knew)=(vbar(Istr-1,j,know)+                     &
     &                         Cx*vbar(Istr,j,knew))/                   &
     &                        (1.0_r8+Cx)
#  ifdef MASKING
          vbar(Istr-1,j,knew)=vbar(Istr-1,j,knew)*                      &
     &                        GRID(ng)%vmask(Istr-1,j)
#  endif
        END DO

# elif defined WEST_M2CLAMPED
!
!  Western edge, clamped boundary condition.
!
        DO j=JstrV,Jend
          vbar(Istr-1,j,knew)=BOUNDARY(ng)%vbar_west(j)
#  ifdef MASKING
          vbar(Istr-1,j,knew)=vbar(Istr-1,j,knew)*                      &
     &                        GRID(ng)%vmask(Istr-1,j)
#  endif
        END DO

# elif defined WEST_M2GRADIENT
!
!  Western edge, gradient boundary condition.
!
        DO j=JstrV,Jend
          vbar(Istr-1,j,knew)=vbar(Istr,j,knew)
#  ifdef MASKING
          vbar(Istr-1,j,knew)=vbar(Istr-1,j,knew)*                      &
     &                        GRID(ng)%vmask(Istr-1,j)
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
          vbar(Istr-1,j,knew)=gamma2*vbar(Istr,j,knew)
#  ifdef MASKING
          vbar(Istr-1,j,knew)=vbar(Istr-1,j,knew)*                      &
     &                        GRID(ng)%vmask(Istr-1,j)
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
          grad(Iend  ,j)=vbar(Iend  ,j+1,know)-                         &
     &                   vbar(Iend  ,j  ,know)
          grad(Iend+1,j)=vbar(Iend+1,j+1,know)-                         &
     &                   vbar(Iend+1,j  ,know)
        END DO
        DO j=JstrV,Jend
          dVdt=vbar(Iend,j,know)-vbar(Iend  ,j,knew)
          dVdx=vbar(Iend,j,knew)-vbar(Iend-1,j,knew)
#  ifdef EAST_M2NUDGING
          tau=M2obc_out(ng,ieast)
          IF ((dVdt*dVdx).lt.0.0_r8) tau=M2obc_in(ng,ieast)
          tau=tau*dt2d
#  endif
          IF ((dVdt*dVdx).lt.0.0_r8) dVdt=0.0_r8
          IF ((dVdt*(grad(Iend,j-1)+grad(Iend,j))).gt.0.0_r8) THEN
            dVde=grad(Iend,j-1)
          ELSE
            dVde=grad(Iend,j  )
          END IF
          cff=MAX(dVdx*dVdx+dVde*dVde,eps)
          Cx=dVdt*dVdx
#  ifdef RADIATION_2D
          Ce=MIN(cff,MAX(dVdt*dVde,-cff))
#  else
          Ce=0.0_r8
#  endif
          vbar(Iend+1,j,knew)=(cff*vbar(Iend+1,j,know)+                 &
     &                         Cx*vbar(Iend,j,knew)-                    &
     &                         MAX(Ce,0.0_r8)*grad(Iend+1,j-1)-         &
     &                         MIN(Ce,0.0_r8)*grad(Iend+1,j  ))/        &
     &                        (cff+Cx)
#  ifdef EAST_M2NUDGING
          vbar(Iend+1,j,knew)=vbar(Iend+1,j,knew)+                      &
     &                        tau*(BOUNDARY(ng)%vbar_east(j)-           &
     &                             vbar(Iend+1,j,know))
#  endif
#  ifdef MASKING
          vbar(Iend+1,j,knew)=vbar(Iend+1,j,knew)*                      &
     &                        GRID(ng)%vmask(Iend+1,j)
#  endif
        END DO

# elif defined EAST_M2FLATHER
!
!  Eastern edge, Chapman boundary condition.
!
        DO j=JstrV,Jend
          Cx=dt2d*0.5_r8*(GRID(ng)%pm(Iend,j-1)+                        &
     &                    GRID(ng)%pm(Iend,j  ))*                       &
     &       SQRT(g*0.5_r8*(GRID(ng)%h(Iend,j-1)+                       &
     &                      GRID(ng)%h(Iend,j  )))
          vbar(Iend+1,j,knew)=(vbar(Iend+1,j,know)+                     &
     &                         Cx*vbar(Iend,j,knew))/                   &
     &                        (1.0_r8+Cx)
#  ifdef MASKING
          vbar(Iend+1,j,knew)=vbar(Iend+1,j,knew)*                      &
     &                        GRID(ng)%vmask(Iend+1,j)
#  endif
        END DO

# elif defined EAST_M2CLAMPED
!
!  Eastern edge, clamped boundary condition.
!
        DO j=JstrV,Jend
          vbar(Iend+1,j,knew)=BOUNDARY(ng)%vbar_east(j)
#  ifdef MASKING
          vbar(Iend+1,j,knew)=vbar(Iend+1,j,knew)*                      &
     &                        GRID(ng)%vmask(Iend+1,j)
#  endif
        END DO

# elif defined EAST_M2GRADIENT
!
!  Eastern edge, gradient boundary condition.
!
        DO j=JstrV,Jend
          vbar(Iend+1,j,knew)=vbar(Iend,j,knew)
#  ifdef MASKING
          vbar(Iend+1,j,knew)=vbar(Iend+1,j,knew)*                      &
     &                        GRID(ng)%vmask(Iend+1,j)
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
          vbar(Iend+1,j,knew)=gamma2*vbar(Iend,j,knew)
#  ifdef MASKING
          vbar(Iend+1,j,knew)=vbar(Iend+1,j,knew)*                      &
     &                        GRID(ng)%vmask(Iend+1,j)
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
          grad(i,Jstr  )=vbar(i  ,Jstr  ,know)-                          &
     &                   vbar(i-1,Jstr  ,know)
          grad(i,Jstr+1)=vbar(i  ,Jstr+1,know)-                          &
     &                   vbar(i-1,Jstr+1,know)
        END DO
        DO i=Istr,Iend
          dVdt=vbar(i,Jstr+1,know)-vbar(i,Jstr+1,knew)
          dVde=vbar(i,Jstr+1,knew)-vbar(i,Jstr+2,knew)
#  ifdef SOUTH_M2NUDGING
          tau=M2obc_out(ng,isouth)
          IF ((dVdt*dVde).lt.0.0_r8) tau=M2obc_in(ng,isouth)
          tau=tau*dt2d
#  endif
          IF ((dVdt*dVde).lt.0.0_r8) dVdt=0.0_r8
          IF ((dVdt*(grad(i,Jstr+1)+grad(i+1,Jstr+1))).gt.0.0_r8) THEN
            dVdx=grad(i  ,Jstr+1)
          ELSE
            dVdx=grad(i+1,Jstr+1)
          END IF
          cff=MAX(dVdx*dVdx+dVde*dVde,eps)
#  ifdef RADIATION_2D
          Cx=MIN(cff,MAX(dVdt*dVdx,-cff))
#  else
          Cx=0.0_r8
#  endif
          Ce=dVdt*dVde
          vbar(i,Jstr,knew)=(cff*vbar(i,Jstr  ,know)+                   &
     &                       Ce *vbar(i,Jstr+1,knew)-                   &
     &                       MAX(Cx,0.0_r8)*grad(i  ,Jstr)-             &
     &                       MIN(Cx,0.0_r8)*grad(i+1,Jstr))/            &
     &                      (cff+Ce)
#  ifdef SOUTH_M2NUDGING
          vbar(i,Jstr,knew)=vbar(i,Jstr,knew)+                          &
     &                      tau*(BOUNDARY(ng)%vbar_south(i)-            &
     &                      vbar(i,Jstr,know))
#  endif
#  ifdef MASKING
          vbar(i,Jstr,knew)=vbar(i,Jstr,knew)*                          &
     &                      GRID(ng)%vmask(i,Jstr)
#  endif
        END DO

# elif defined SOUTH_M2FLATHER
!
!  Southern edge, Flather boundary condition.
!
        DO i=Istr,Iend
          vbar(i,Jstr,knew)=BOUNDARY(ng)%vbar_south(i)-                 &
     &                      SQRT(g/(0.5_r8*(GRID(ng)%h(i,Jstr-1)+       &
     &                                      GRID(ng)%h(i,Jstr  ))))*    &
     &                      (0.5_r8*(zeta(i,Jstr-1,know)+               &
     &                               zeta(i,Jstr  ,know))-              &
     &                       BOUNDARY(ng)%zeta_south(i))
#  ifdef MASKING
          vbar(i,Jstr,knew)=vbar(i,Jstr,knew)*                          &
     &                      GRID(ng)%vmask(i,Jstr)
#  endif
        END DO

# elif defined SOUTH_M2CLAMPED
!
!  Southern edge, clamped boundary condition.
!
        DO i=Istr,Iend
          vbar(i,Jstr,knew)=BOUNDARY(ng)%vbar_south(i)
#  ifdef MASKING
          vbar(i,Jstr,knew)=vbar(i,Jstr,knew)*                          &
     &                      GRID(ng)%vmask(i,Jstr)
#  endif
        END DO

# elif defined SOUTH_M2GRADIENT
!
!  Southern edge, gradient boundary condition.
!
        DO i=Istr,Iend
          vbar(i,Jstr,knew)=vbar(i,Jstr+1,knew)
#  ifdef MASKING
          vbar(i,Jstr,knew)=vbar(i,Jstr,knew)*                          &
     &                      GRID(ng)%vmask(i,Jstr)
#  endif
        END DO

# else
!
!  Southern edge, closed boundary condition.
!
        DO i=Istr,Iend
          vbar(i,Jstr,knew)=0.0_r8
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
          grad(i,Jend  )=vbar(i  ,Jend  ,know)-                         &
     &                   vbar(i-1,Jend  ,know)
          grad(i,Jend+1)=vbar(i  ,Jend+1,know)-                         &
     &                   vbar(i-1,Jend+1,know)
        END DO
        DO i=Istr,Iend
          dVdt=vbar(i,Jend,know)-vbar(i,Jend  ,knew)
          dVde=vbar(i,Jend,knew)-vbar(i,Jend-1,knew)
#  ifdef NORTH_M2NUDGING
          tau=M2obc_out(ng,inorth)
          IF ((dVdt*dVde).lt.0.0_r8) tau=M2obc_in(ng,inorth)
          tau=tau*dt2d
#  endif
          IF ((dVdt*dVde).lt.0.0_r8) dVdt=0.0_r8
          IF ((dVdt*(grad(i,Jend)+grad(i,Jend+1))).gt.0.0_r8) THEN
            dVdx=grad(i  ,Jend)
          ELSE
            dVdx=grad(i+1,Jend)
          END IF
          cff=MAX(dVdx*dVdx+dVde*dVde,eps)
#  ifdef RADIATION_2D
          Cx=MIN(cff,MAX(dVdt*dVdx,-cff))
#  else
          Cx=0.0_r8
#  endif
          Ce=dVdt*dVde
          vbar(i,Jend+1,knew)=(cff*vbar(i,Jend+1,know)+                 &
     &                         Ce *vbar(i,Jend  ,knew)-                 &
     &                         MAX(Cx,0.0_r8)*grad(i  ,Jend+1)-         &
     &                         MIN(Cx,0.0_r8)*grad(i+1,Jend+1))/        &
     &                        (cff+Ce)
#  ifdef NORTH_M2NUDGING
          vbar(i,Jend+1,knew)=vbar(i,Jend+1,knew)+                      &
     &                         tau*(BOUNDARY(ng)%vbar_north(i)-         &
     &                         vbar(i,Jend+1,know))
#  endif
#  ifdef MASKING
          vbar(i,Jend+1,knew)=vbar(i,Jend+1,knew)*                      &
     &                        GRID(ng)%vmask(i,Jend+1)
#  endif
        END DO

# elif defined NORTH_M2FLATHER
!
!  Northern edge, Flather boundary condition.
!
        DO i=Istr,Iend
          vbar(i,Jend+1,knew)=BOUNDARY(ng)%vbar_north(i)+               &
     &                        SQRT(g/(0.5_r8*(GRID(ng)%h(i,Jend  )+     &
     &                                        GRID(ng)%h(i,Jend+1))))*  &
     &                        (0.5_r8*(zeta(i,Jend  ,know)+             &
     &                                 zeta(i,Jend+1,know))-            &
     &                         BOUNDARY(ng)%zeta_north(i))
#  ifdef MASKING
          vbar(i,Jend+1,knew)=vbar(i,Jend+1,knew)*                      &
     &                        GRID(ng)%vmask(i,Jend+1)
#  endif
        END DO

# elif defined NORTH_M2CLAMPED
!
!  Northern edge, clamped boundary condition.
!
        DO i=Istr,Iend
          vbar(i,Jend+1,knew)=BOUNDARY(ng)%vbar_north(i)
#  ifdef MASKING
          vbar(i,Jend+1,knew)=vbar(i,Jend+1,knew)*                      &
     &                        GRID(ng)%vmask(i,Jend+1)
#  endif
        END DO

# elif defined NORTH_M2GRADIENT
!
!  Northern edge, gradient boundary condition.
!
        DO i=Istr,Iend
          vbar(i,Jend+1,knew)=vbar(i,Jend,knew)
#  ifdef MASKING
          vbar(i,Jend+1,knew)=vbar(i,Jend+1,knew)*                      &
     &                        GRID(ng)%vmask(i,Jend+1)
#  endif
        END DO
# else
!
!  Northern edge, closed boundary condition.
!
        DO i=Istr,Iend
          vbar(i,Jend+1,knew)=0.0_r8
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
        vbar(Istr-1,Jstr,knew)=0.5_r8*(vbar(Istr-1,Jstr+1,knew)+        &
     &                                 vbar(Istr  ,Jstr  ,knew))
      END IF
      IF (SOUTHERN_EDGE.and.EASTERN_EDGE) THEN
        vbar(Iend+1,Jstr,knew)=0.5_r8*(vbar(Iend  ,Jstr  ,knew)+        &
     &                                 vbar(Iend+1,Jstr+1,knew))
      END IF
      IF (NORTHERN_EDGE.and.WESTERN_EDGE) THEN
        vbar(Istr-1,Jend+1,knew)=0.5_r8*(vbar(Istr-1,Jend  ,knew)+      &
     &                                   vbar(Istr  ,Jend+1,knew))
      END IF
      IF (NORTHERN_EDGE.and.EASTERN_EDGE) THEN
        vbar(Iend+1,Jend+1,knew)=0.5_r8*(vbar(Iend  ,Jend+1,knew)+      &
     &                                   vbar(Iend+1,Jend  ,knew))
      END IF
#endif
      RETURN

      END SUBROUTINE v2dbc_tile
      END MODULE v2dbc_mod
