#include "cppdefs.h"
      MODULE u3dbc_mod
#ifdef SOLVE3D
!
!========================================== Alexander F. Shchepetkin ===
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This subroutine sets lateral boundary conditions for total 3D       !
!  U-velocity.                                                         !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC u3dbc_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE u3dbc (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_ocean
      USE mod_stepping
!
      integer, intent(in) :: ng, tile

#include "tile.h"
!
      CALL u3dbc_tile (ng, Istr, Iend, Jstr, Jend,                      &
     &                 LBi, UBi, LBj, UBj, N(ng),                       &
     &                 nstp(ng), nnew(ng),                              &
     &                 OCEAN(ng) % u)
      RETURN
      END SUBROUTINE u3dbc
!
!***********************************************************************
      SUBROUTINE u3dbc_tile (ng, Istr, Iend, Jstr, Jend,                &
     &                       LBi, UBi, LBj, UBj, UBk,                   &
     &                       nstp, nnew,                                &
     &                       u)
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
      integer, intent(in) :: LBi, UBi, LBj, UBj, UBk
      integer, intent(in) :: nstp, nnew

# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: u(LBi:,LBj:,:,:)
# else
      real(r8), intent(inout) :: u(LBi:UBi,LBj:UBj,UBk,2)
# endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j, k

      real(r8), parameter :: eps = 1.0E-20_r8

      real(r8) :: Ce, Cx, cff, dUde, dUdt, dUdx, tau

      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: grad

# include "set_bounds.h"

# ifndef EW_PERIODIC
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the western edge.
!-----------------------------------------------------------------------
!
      IF (WESTERN_EDGE) THEN

#  if defined WEST_M3RADIATION
!
!  Western edge, implicit upstream radiation condition.
!
        DO k=1,N(ng)
          DO j=Jstr,Jend+1
            grad(1,j)=u(1,j  ,k,nstp)-                                  &
     &                u(1,j-1,k,nstp)
            grad(2,j)=u(2,j  ,k,nstp)-                                  &
     &                u(2,j-1,k,nstp)
          END DO
          DO j=Jstr,Jend
            dUdt=u(2,j,k,nstp)-u(2,j,k,nnew)
            dUdx=u(2,j,k,nnew)-u(3,j,k,nnew)
#   ifdef WEST_M3NUDGING
            tau=M3obc_out(ng,iwest)
            IF ((dUdt*dUdx).lt.0.0_r8) tau=M3obc_in(ng,iwest)
            tau=tau*dt(ng)
#   endif
            IF ((dUdt*dUdx).lt.0.0_r8) dUdt=0.0_r8
            IF ((dUdt*(grad(2,j)+grad(2,j+1))).gt.0.0_r8) THEN
              dUde=grad(2,j)
            ELSE
              dUde=grad(2,j+1)
            END IF
            cff=MAX(dUdx*dUdx+dUde*dUde,eps)
            Cx=dUdt*dUdx
#   ifdef RADIATION_2D
            Ce=MIN(cff,MAX(dUdt*dUde,-cff))
#   else
            Ce=0.0_r8
#   endif
            u(1,j,k,nnew)=(cff*u(1,j,k,nstp)+                           &
     &                     Cx*u(2,j,k,nnew)-                            &
     &                     MAX(Ce,0.0_r8)*grad(1,j  )-                  &
     &                     MIN(Ce,0.0_r8)*grad(1,j+1))/                 &
     &                    (cff+Cx)
#   ifdef WEST_M3NUDGING
            u(1,j,k,nnew)=u(1,j,k,nnew)+                                &
     &                    tau*(BOUNDARY(ng)%u_west(j,k)-                &
     &                         u(1,j,k,nstp))
#   endif
#   ifdef MASKING
            u(1,j,k,nnew)=u(1,j,k,nnew)*                                &
     &                    GRID(ng)%umask(1,j)
#   endif
          END DO
        END DO
#  elif defined WEST_M3CLAMPED
!
!  Western edge, clamped boundary condition.
!
        DO k=1,N(ng)
          DO j=Jstr,Jend
            u(1,j,k,nnew)=BOUNDARY(ng)%u_west(j,k)
#   ifdef MASKING
            u(1,j,k,nnew)=u(1,j,k,nnew)*                                &
     &                    GRID(ng)%umask(1,j)
#   endif
          END DO
        END DO
#  elif defined WEST_M3GRADIENT
!
!  Western edge, gradient boundary condition.
!
        DO k=1,N(ng)
          DO j=Jstr,Jend
            u(1,j,k,nnew)=u(2,j,k,nnew)
#   ifdef MASKING
            u(1,j,k,nnew)=u(1,j,k,nnew)*                                &
     &                    GRID(ng)%umask(1,j)
#   endif
          END DO
        END DO
#  else
!
!  Western edge, closed boundary condition.
!
        DO k=1,N(ng)
          DO j=Jstr,Jend
            u(1,j,k,nnew)=0.0_r8
          END DO
        END DO
#  endif
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the eastern edge.
!-----------------------------------------------------------------------
!
      IF (EASTERN_EDGE) THEN

#  if defined EAST_M3RADIATION
!
!  Eastern edge, implicit upstream radiation condition.
!
        DO k=1,N(ng)
          DO j=Jstr,Jend+1
            grad(Lm(ng),j)=u(Lm(ng),j  ,k,nstp)-                        &
     &                     u(Lm(ng),j-1,k,nstp)
            grad(Lm(ng)+1,j)=u(Lm(ng)+1,j  ,k,nstp)-                    &
     &                       u(Lm(ng)+1,j-1,k,nstp)
          END DO
          DO j=Jstr,Jend
            dUdt=u(Lm(ng),j,k,nstp)-u(Lm(ng),j,k,nnew)
            dUdx=u(Lm(ng),j,k,nnew)-u(Lm(ng)-1,j,k,nnew)
#   ifdef EAST_M3NUDGING
            tau=M3obc_out(ng,ieast)
            IF ((dUdt*dUdx).lt.0.0_r8) tau=M3obc_in(ng,ieast)
            tau=tau*dt(ng)
#   endif
            IF ((dUdt*dUdx).lt.0.0_r8) dUdt=0.0_r8
            IF ((dUdt*(grad(Lm(ng),j)+grad(Lm(ng),j+1))).gt.0.0_r8) THEN
              dUde=grad(Lm(ng),j  )
            ELSE
              dUde=grad(Lm(ng),j+1)
            END IF
            cff=MAX(dUdx*dUdx+dUde*dUde,eps)
            Cx=dUdt*dUdx
#   ifdef RADIATION_2D
            Ce=MIN(cff,MAX(dUdt*dUde,-cff))
#   else
            Ce=0.0_r8
#   endif
            u(Lm(ng)+1,j,k,nnew)=(cff*u(Lm(ng)+1,j,k,nstp)+             &
     &                            Cx*u(Lm(ng),j,k,nnew)-                &
     &                            MAX(Ce,0.0_r8)*grad(Lm(ng)+1,j  )-    &
     &                            MIN(Ce,0.0_r8)*grad(Lm(ng)+1,j+1))/   &
     &                           (cff+Cx)
#   ifdef EAST_M3NUDGING
            u(Lm(ng)+1,j,k,nnew)=u(Lm(ng)+1,j,k,nnew)+                  &
     &                        tau*(BOUNDARY(ng)%u_east(j,k)-            &
     &                             u(Lm(ng)+1,j,k,nstp))
#   endif
#   ifdef MASKING
            u(Lm(ng)+1,j,k,nnew)=u(Lm(ng)+1,j,k,nnew)*                  &
     &                           GRID(ng)%umask(Lm(ng)+1,j)
#   endif
          END DO
        END DO
#  elif defined EAST_M3CLAMPED
!
!  Eastern edge, clamped boundary condition.
!
        DO k=1,N(ng)
          DO j=Jstr,Jend
            u(Lm(ng)+1,j,k,nnew)=BOUNDARY(ng)%u_east(j,k)
#   ifdef MASKING
            u(Lm(ng)+1,j,k,nnew)=u(Lm(ng)+1,j,k,nnew)*                  &
     &                           GRID(ng)%umask(Lm(ng)+1,j)
#   endif
          END DO
        END DO
#  elif defined EAST_M3GRADIENT
!
!  Eastern edge, gradient boundary condition.
!
      DO k=1,N(ng)
        DO j=Jstr,Jend
          u(Lm(ng)+1,j,k,nnew)=u(Lm(ng),j,k,nnew)
#   ifdef MASKING
          u(Lm(ng)+1,j,k,nnew)=u(Lm(ng)+1,j,k,nnew)*                    &
     &                         GRID(ng)%umask(Lm(ng)+1,j)
#   endif
        END DO
      END DO
#  else
!
!  Eastern edge, closed boundary condition.
!
        DO k=1,N(ng)
          DO j=Jstr,Jend
            u(Lm(ng)+1,j,k,nnew)=0.0_r8
          END DO
        END DO
#  endif
      END IF
# endif
# ifndef NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the southern edge.
!-----------------------------------------------------------------------
!
      IF (SOUTHERN_EDGE) THEN

#  if defined SOUTH_M3RADIATION
!
!  Southern edge, implicit upstream radiation condition.
!
        DO k=1,N(ng)
          DO i=IstrU-1,Iend
            grad(i,0)=u(i+1,0,k,nstp)-                                  &
     &                u(i  ,0,k,nstp)
            grad(i,1)=u(i+1,1,k,nstp)-                                  &
     &                u(i  ,1,k,nstp)
          END DO
          DO i=IstrU,Iend
            dUdt=u(i,1,k,nstp)-u(i,1,k,nnew)
            dUde=u(i,1,k,nnew)-u(i,2,k,nnew)
#   ifdef SOUTH_M3NUDGING
            tau=M3obc_out(ng,isouth)
            IF ((dUdt*dUde).lt.0.0_r8) tau=M3obc_in(ng,isouth)
            tau=tau*dt(ng)
#   endif
            IF ((dUdt*dUde).lt.0.0_r8) dUdt=0.0_r8
            IF ((dUdt*(grad(i-1,1)+grad(i,1))).gt.0.0_r8) THEN
              dUdx=grad(i-1,1)
            ELSE
              dUdx=grad(i  ,1)
            END IF
            cff=MAX(dUdx*dUdx+dUde*dUde,eps)
#   ifdef RADIATION_2D
            Cx=MIN(cff,MAX(dUdt*dUdx,-cff))
#   else
            Cx=0.0_r8
#   endif
            Ce=dUdt*dUde
            u(i,0,k,nnew)=(cff*u(i,0,k,nstp)+                           &
     &                     Ce*u(i,1,k,nnew)-                            &
     &                     MAX(Cx,0.0_r8)*grad(i-1,0)-                  &
     &                     MIN(Cx,0.0_r8)*grad(i  ,0))/                 &
     &                    (cff+Ce)
#   ifdef SOUTH_M3NUDGING
            u(i,0,k,nnew)=u(i,0,k,nnew)+                                &
     &                    tau*(BOUNDARY(ng)%u_south(i,k)-               &
     &                         u(i,0,k,nstp))
#   endif
#   ifdef MASKING
            u(i,0,k,nnew)=u(i,0,k,nnew)*                                &
     &                    GRID(ng)%umask(i,0)
#   endif
          END DO
        END DO
#  elif defined SOUTH_M3CLAMPED
!
!  Southern edge, clamped boundary condition.
!
        DO k=1,N(ng)
          DO i=IstrU,Iend
            u(i,0,k,nnew)=BOUNDARY(ng)%u_south(i,k)
#   ifdef MASKING
            u(i,0,k,nnew)=u(i,0,k,nnew)*                                &
     &                    GRID(ng)%umask(i,0)
#   endif
          END DO
        END DO
#  elif defined SOUTH_M3GRADIENT
!
!  Southern edge, gradient boundary condition.
!
        DO k=1,N(ng)
          DO i=IstrU,Iend
            u(i,0,k,nnew)=u(i,1,k,nnew)
#   ifdef MASKING
            u(i,0,k,nnew)=u(i,0,k,nnew)*                                &
     &                    GRID(ng)%umask(i,0)
#   endif
          END DO
        END DO
#  else
!
!  Southern edge, closed boundary condition: free slip (gamma2=1)  or
!                                            no   slip (gamma2=-1).
!
#   ifdef EW_PERIODIC
#    define I_RANGE IstrU,Iend
#   else
#    define I_RANGE Istr,IendR
#   endif
        DO k=1,N(ng)
          DO i=I_RANGE
            u(i,0,k,nnew)=gamma2*u(i,1,k,nnew)
#   ifdef MASKING
            u(i,0,k,nnew)=u(i,0,k,nnew)*                                &
     &                    GRID(ng)%umask(i,0)
#   endif
          END DO
        END DO
#   undef I_RANGE
#  endif
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the northern edge.
!-----------------------------------------------------------------------
!
      IF (NORTHERN_EDGE) THEN

#  if defined NORTH_M3RADIATION
!
!  Northern edge, implicit upstream radiation condition.
!
        DO k=1,N(ng)
          DO i=IstrU-1,Iend
            grad(i,Mm(ng))=u(i+1,Mm(ng),k,nstp)-                        &
     &                     u(i  ,Mm(ng),k,nstp)
            grad(i,Mm(ng)+1)=u(i+1,Mm(ng)+1,k,nstp)-                    &
     &                       u(i  ,Mm(ng)+1,k,nstp)
          END DO
          DO i=IstrU,Iend
            dUdt=u(i,Mm(ng),k,nstp)-u(i,Mm(ng),k,nnew)
            dUde=u(i,Mm(ng),k,nnew)-u(i,Mm(ng)-1,k,nnew)
#   ifdef NORTH_M3NUDGING
            tau=M3obc_out(ng,inorth)
            IF ((dUdt*dUde).lt.0.0_r8) tau=M3obc_in(ng,inorth)
            tau=tau*dt(ng)
#   endif
            IF ((dUdt*dUde).lt.0.0_r8) dUdt=0.0_r8
            IF ((dUdt*(grad(i-1,Mm(ng))+grad(i,Mm(ng)))).gt.0.0_r8) THEN
              dUdx=grad(i-1,Mm(ng))
            ELSE
              dUdx=grad(i  ,Mm(ng))
            END IF
            cff=MAX(dUdx*dUdx+dUde*dUde,eps)
#   ifdef RADIATION_2D
            Cx=MIN(cff,MAX(dUdt*dUdx,-cff))
#   else
            Cx=0.0_r8
#   endif
            Ce=dUdt*dUde
            u(i,Mm(ng)+1,k,nnew)=(cff*u(i,Mm(ng)+1,k,nstp)+             &
     &                            Ce*u(i,Mm(ng),k,nnew)-                &
     &                            MAX(Cx,0.0_r8)*grad(i-1,Mm(ng)+1)-    &
     &                            MIN(Cx,0.0_r8)*grad(i  ,Mm(ng)+1))/   &
     &                           (cff+Ce)
#   ifdef NORTH_M3NUDGING
            u(i,Mm(ng)+1,k,nnew)=u(i,Mm(ng)+1,k,nnew)+                  &
     &                           tau*(BOUNDARY(ng)%u_north(i,k)-        &
     &                                u(i,Mm(ng)+1,k,nstp))
#    endif
#   ifdef MASKING
            u(i,Mm(ng)+1,k,nnew)=u(i,Mm(ng)+1,k,nnew)*                  &
     &                           GRID(ng)%umask(i,Mm(ng)+1)
#   endif
          END DO
        END DO
#  elif defined NORTH_M3CLAMPED
!
!  Northern edge, clamped boundary condition.
!
        DO k=1,N(ng)
          DO i=IstrU,Iend
            u(i,Mm(ng)+1,k,nnew)=BOUNDARY(ng)%u_north(i,k)
#   ifdef MASKING
            u(i,Mm(ng)+1,k,nnew)=u(i,Mm(ng)+1,k,nnew)*                  &
     &                           GRID(ng)%umask(i,Mm(ng)+1)
#   endif
          END DO
        END DO
#  elif defined NORTH_M3GRADIENT
!
!  Northern edge, gradient boundary condition.
!
        DO k=1,N(ng)
          DO i=IstrU,Iend
            u(i,Mm(ng)+1,k,nnew)=u(i,Mm(ng),k,nnew)
#   ifdef MASKING
            u(i,Mm(ng)+1,k,nnew)=u(i,Mm(ng)+1,k,nnew)*                  &
     &                           GRID(ng)%umask(i,Mm(ng)+1)
#   endif
          END DO
        END DO
#  else
!
!  Northern edge, closed boundary condition: free slip (gamma2=1)  or
!                                            no   slip (gamma2=-1).
!
#   ifdef EW_PERIODIC
#    define I_RANGE IstrU,Iend
#   else
#    define I_RANGE Istr,IendR
#   endif
        DO k=1,N(ng)
          DO i=I_RANGE
            u(i,Mm(ng)+1,k,nnew)=gamma2*u(i,Mm(ng),k,nnew)
#   ifdef MASKING
            u(i,Mm(ng)+1,k,nnew)=u(i,Mm(ng)+1,k,nnew)*                  &
     &                           GRID(ng)%umask(i,Mm(ng)+1)
#   endif
          END DO
        END DO
#   undef I_RANGE
#  endif
      END IF
# endif
# if !defined EW_PERIODIC && !defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (SOUTHERN_EDGE.and.WESTERN_EDGE) THEN
        DO k=1,N(ng)
          u(1,0,k,nnew)=0.5_r8*(u(2,0,k,nnew)+                          &
     &                          u(1,1,k,nnew))
        END DO
      END IF
      IF (SOUTHERN_EDGE.and.EASTERN_EDGE) THEN
        DO k=1,N(ng)
          u(Lm(ng)+1,0,k,nnew)=0.5_r8*(u(Lm(ng)  ,0,k,nnew)+            &
     &                                 u(Lm(ng)+1,1,k,nnew))
        END DO
      END IF
      IF (NORTHERN_EDGE.and.WESTERN_EDGE) THEN
        DO k=1,N(ng)
          u(1,Mm(ng)+1,k,nnew)=0.5_r8*(u(2,Mm(ng)+1,k,nnew)+            &
     &                                 u(1,Mm(ng),k,nnew))
        END DO
      END IF
      IF (NORTHERN_EDGE.and.EASTERN_EDGE) THEN
        DO k=1,N(ng)
          u(Lm(ng)+1,Mm(ng)+1,k,nnew)=0.5_r8*                           &
     &                                (u(Lm(ng)  ,Mm(ng)+1,k,nnew)+     &
     &                                 u(Lm(ng)+1,Mm(ng)  ,k,nnew))
        END DO
      END IF
# endif
      RETURN
      END SUBROUTINE u3dbc_tile
#endif
      END MODULE u3dbc_mod
