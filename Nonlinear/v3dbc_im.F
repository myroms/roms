#include "cppdefs.h"
      MODULE v3dbc_mod
#ifdef SOLVE3D
!
!========================================== Alexander F. Shchepetkin ===
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This subroutine sets lateral boundary conditions for total 3D       !
!  V-velocity.                                                         !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC v3dbc_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE v3dbc (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_ocean
      USE mod_stepping
!
      integer, intent(in) :: ng, tile

#include "tile.h"
!
      CALL v3dbc_tile (ng, Istr, Iend, Jstr, Jend,                      &
     &                 LBi, UBi, LBj, UBj, N(ng),                       &
     &                 nstp(ng), nnew(ng),                              &
     &                 OCEAN(ng) % v)
      RETURN
      END SUBROUTINE v3dbc
!
!***********************************************************************
      SUBROUTINE v3dbc_tile (ng, Istr, Iend, Jstr, Jend,                &
     &                       LBi, UBi, LBj, UBj, UBk,                   &
     &                       nstp, nnew,                                &
     &                       v)
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
      real(r8), intent(inout) :: v(LBi:,LBj:,:,:)
# else
      real(r8), intent(inout) :: v(LBi:UBi,LBj:UBj,UBk,2)
# endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j, k

      real(r8), parameter :: eps = 1.0E-20_r8

      real(r8) :: Ce, Cx, cff, dVde, dVdt, dVdx, tau

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
          DO j=JstrV-1,Jend
            grad(0,j)=v(0,j+1,k,nstp)-                                  &
     &                v(0,j  ,k,nstp)
            grad(1,j)=v(1,j+1,k,nstp)-                                  &
     &                v(1,j  ,k,nstp)
          END DO
          DO j=JstrV,Jend
            dVdt=v(1,j,k,nstp)-v(1,j,k,nnew)
            dVdx=v(1,j,k,nnew)-v(2,j,k,nnew)
#   ifdef WEST_M3NUDGING
            tau=M3obc_out(ng,iwest)
            IF ((dVdt*dVdx).lt.0.0_r8) tau=M3obc_in(ng,iwest)
            tau=tau*dt(ng)
#   endif
            IF ((dVdt*dVdx).lt.0.0_r8) dVdt=0.0_r8
            IF ((dVdt*(grad(1,j-1)+grad(1,j))).gt.0.0_r8) THEN
              dVde=grad(1,j-1)
            ELSE
              dVde=grad(1,j  )
            END IF
            cff=MAX(dVdx*dVdx+dVde*dVde,eps)
            Cx=dVdt*dVdx
#   ifdef RADIATION_2D
            Ce=MIN(cff,MAX(dVdt*dVde,-cff))
#   else
            Ce=0.0_r8
#   endif
            v(0,j,k,nnew)=(cff*v(0,j,k,nstp)+                           &
     &                     Cx*v(1,j,k,nnew)-                            &
     &                     MAX(Ce,0.0_r8)*grad(0,j-1)-                  &
     &                     MIN(Ce,0.0_r8)*grad(0,j  ))/                 &
     &                    (cff+Cx)
#   ifdef WEST_M3NUDGING
            v(0,j,k,nnew)=v(0,j,k,nnew)+                                &
     &                    tau*(BOUNDARY(ng)%v_west(j,k)-                &
     &                         v(0,j,k,nstp))
#    endif
#   ifdef MASKING
            v(0,j,k,nnew)=v(0,j,k,nnew)*                                &
     &                    GRID(ng)%vmask(0,j)
#   endif
          END DO
        END DO
#  elif defined WEST_M3CLAMPED
!
!  Western edge, clamped boundary condition.
!
        DO k=1,N(ng)
          DO j=JstrV,Jend
            v(0,j,k,nnew)=BOUNDARY(ng)%v_west(j,k)
#   ifdef MASKING
            v(0,j,k,nnew)=v(0,j,k,nnew)*                                &
     &                    GRID(ng)%vmask(0,j)
#   endif
          END DO
        END DO
#  elif defined WEST_M3GRADIENT
!
!  Western edge, gradient boundary condition.
!
        DO k=1,N(ng)
          DO j=JstrV,Jend
            v(0,j,k,nnew)=v(1,j,k,nnew)
#   ifdef MASKING
            v(0,j,k,nnew)=v(0,j,k,nnew)*                                &
     &                    GRID(ng)%vmask(0,j)
#   endif
          END DO
        END DO
#  else
!
!  Western edge, closed boundary condition: free slip (gamma2=1)  or
!                                           no   slip (gamma2=-1).
!
#   ifdef NS_PERIODIC
#    define J_RANGE JstrV,Jend
#   else
#    define J_RANGE Jstr,JendR
#   endif
        DO k=1,N(ng)
          DO j=J_RANGE
            v(0,j,k,nnew)=gamma2*v(1,j,k,nnew)
#   ifdef MASKING
            v(0,j,k,nnew)=v(0,j,k,nnew)*                                &
     &                    GRID(ng)%vmask(0,j)
#   endif
          END DO
        END DO
#   undef J_RANGE
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
          DO j=JstrV-1,Jend
            grad(Lm(ng),j)=v(Lm(ng),j+1,k,nstp)-                        &
     &                     v(Lm(ng),j  ,k,nstp)
            grad(Lm(ng)+1,j)=v(Lm(ng)+1,j+1,k,nstp)-                    &
     &                       v(Lm(ng)+1,j  ,k,nstp)
          END DO
          DO j=JstrV,Jend
            dVdt=v(Lm(ng),j,k,nstp)-v(Lm(ng),j,k,nnew)
            dVdx=v(Lm(ng),j,k,nnew)-v(Lm(ng)-1,j,k,nnew)
#   ifdef EAST_M3NUDGING
            tau=M3obc_out(ng,ieast)
            IF ((dVdt*dVdx).lt.0.0_r8) tau=M3obc_in(ng,ieast)
            tau=tau*dt(ng)
#   endif
            IF ((dVdt*dVdx).lt.0.0_r8) dVdt=0.0_r8
            IF ((dVdt*(grad(Lm(ng),j-1)+grad(Lm(ng),j))).gt.0.0_r8) THEN
              dVde=grad(Lm(ng),j-1)
            ELSE
              dVde=grad(Lm(ng),j  )
            END IF
            cff=MAX(dVdx*dVdx+dVde*dVde,eps)
            Cx=dVdt*dVdx
#   ifdef RADIATION_2D
            Ce=MIN(cff,MAX(dVdt*dVde,-cff))
#   else
            Ce=0.0_r8
#   endif
            v(Lm(ng)+1,j,k,nnew)=(cff*v(Lm(ng)+1,j,k,nstp)+             &
     &                            Cx*v(Lm(ng),j,k,nnew)-                &
     &                            MAX(Ce,0.0_r8)*grad(Lm(ng)+1,j-1)-    &
     &                            MIN(Ce,0.0_r8)*grad(Lm(ng)+1,j  ))/   &
     &                           (cff+Cx)
#   ifdef EAST_M3NUDGING
            v(Lm(ng)+1,j,k,nnew)=v(Lm(ng)+1,j,k,nnew)+                  &
     &                           tau*(BOUNDARY(ng)%v_east(j,k)-         &
     &                                v(Lm(ng)+1,j,k,nstp))
#   endif
#   ifdef MASKING
            v(Lm(ng)+1,j,k,nnew)=v(Lm(ng)+1,j,k,nnew)*                  &
     &                           GRID(ng)%vmask(Lm(ng)+1,j)
#   endif
          END DO
        END DO
#  elif defined EAST_M3CLAMPED
!
!  Eastern edge, clamped boundary condition.
!
        DO k=1,N(ng)
          DO j=JstrV,Jend
            v(Lm(ng)+1,j,k,nnew)=BOUNDARY(ng)%v_east(j,k)
#   ifdef MASKING
            v(Lm(ng)+1,j,k,nnew)=v(Lm(ng)+1,j,k,nnew)*                  &
     &                           GRID(ng)%vmask(Lm(ng)+1,j)
#   endif
          END DO
        END DO
#  elif defined EAST_M3GRADIENT
!
!  Eastern edge, gradient boundary condition.
!
        DO k=1,N(ng)
          DO j=JstrV,Jend
            v(Lm(ng)+1,j,k,nnew)=v(Lm(ng),j,k,nnew)
#   ifdef MASKING
            v(Lm(ng)+1,j,k,nnew)=v(Lm(ng)+1,j,k,nnew)*                  &
     &                           GRID(ng)%vmask(Lm(ng)+1,j)
#   endif
          END DO
        END DO
#  else
!
!  Eastern edge, closed boundary condition: free slip (gamma2=1)  or
!                                           no   slip (gamma2=-1).
!
#   ifdef NS_PERIODIC
#    define J_RANGE JstrV,Jend
#   else
#    define J_RANGE Jstr,JendR
#   endif
        DO k=1,N(ng)
          DO j=J_RANGE
            v(Lm(ng)+1,j,k,nnew)=gamma2*v(Lm(ng),j,k,nnew)
#   ifdef MASKING
            v(Lm(ng)+1,j,k,nnew)=v(Lm(ng)+1,j,k,nnew)*                  &
     &                           GRID(ng)%vmask(Lm(ng)+1,j)
#   endif
          END DO
        END DO
#   undef J_RANGE
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
          DO i=Istr,Iend+1
            grad(i,1)=v(i  ,1,k,nstp)-                                  &
     &                v(i-1,1,k,nstp)
            grad(i,2)=v(i  ,2,k,nstp)-                                  &
     &                v(i-1,2,k,nstp)
          END DO
          DO i=Istr,Iend
            dVdt=v(i,2,k,nstp)-v(i,2,k,nnew)
            dVde=v(i,2,k,nnew)-v(i,3,k,nnew)
#   ifdef SOUTH_M3NUDGING
            tau=M3obc_out(ng,isouth)
            IF ((dVdt*dVde).lt.0.0_r8) tau=M3obc_in(ng,isouth)
            tau=tau*dt(ng)
#   endif
            IF ((dVdt*dVde).lt.0.0_r8) dVdt=0.0_r8
            IF ((dVdt*(grad(i,2)+grad(i+1,2))).gt.0.0_r8) THEN
              dVdx=grad(i  ,2)
            ELSE
              dVdx=grad(i+1,2)
            END IF
            cff=MAX(dVdx*dVdx+dVde*dVde,eps)
#   ifdef RADIATION_2D
            Cx=MIN(cff,MAX(dVdt*dVdx,-cff))
#   else
            Cx=0.0_r8
#   endif
            Ce=dVdt*dVde
            v(i,1,k,nnew)=(cff*v(i,1,k,nstp)+                           &
     &                     Ce*v(i,2,k,nnew)-                            &
     &                     MAX(Cx,0.0_r8)*grad(i  ,1)-                  &
     &                     MIN(Cx,0.0_r8)*grad(i+1,1))/                 &
     &                    (cff+Ce)
#   ifdef SOUTH_M3NUDGING
            v(i,1,k,nnew)=v(i,1,k,nnew)+                                &
     &                    tau*(BOUNDARY(ng)%v_south(i,k)-               &
     &                         v(i,1,k,nstp))
#   endif
#   ifdef MASKING
            v(i,1,k,nnew)=v(i,1,k,nnew)*                                &
     &                    GRID(ng)%vmask(i,1)
#   endif
          END DO
        END DO

#  elif defined SOUTH_M3CLAMPED
!
!  Southern edge, clamped boundary condition.
!
        DO k=1,N(ng)
          DO i=Istr,Iend
            v(i,1,k,nnew)=BOUNDARY(ng)%v_south(i,k)
#   ifdef MASKING
            v(i,1,k,nnew)=v(i,1,k,nnew)*                                &
     &                    GRID(ng)%vmask(i,1)
#   endif
          END DO
        END DO
#  elif defined SOUTH_M3GRADIENT
!
!  Southern edge, gradient boundary condition.
!
        DO k=1,N(ng)
          DO i=Istr,Iend
            v(i,1,k,nnew)=v(i,2,k,nnew)
#   ifdef MASKING
            v(i,1,k,nnew)=v(i,1,k,nnew)*                                &
     &                    GRID(ng)%vmask(i,1)
#   endif
          END DO
        END DO
#  else
!
!  Southern edge, closed boundary condition.
!
        DO k=1,N(ng)
          DO i=Istr,Iend
            v(i,1,k,nnew)=0.0_r8
          END DO
        END DO
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
          DO i=Istr,Iend+1
            grad(i,Mm(ng))=v(i  ,Mm(ng),k,nstp)-                        &
     &                     v(i-1,Mm(ng),k,nstp)
            grad(i,Mm(ng)+1)=v(i  ,Mm(ng)+1,k,nstp)-                    &
     &                       v(i-1,Mm(ng)+1,k,nstp)
          END DO
          DO i=Istr,Iend
            dVdt=v(i,Mm(ng),k,nstp)-v(i,Mm(ng),k,nnew)
            dVde=v(i,Mm(ng),k,nnew)-v(i,Mm(ng)-1,k,nnew)
#   ifdef NORTH_M3NUDGING
            tau=M3obc_out(ng,inorth)
            IF ((dVdt*dVde).lt.0.0_r8) tau=M3obc_in(ng,inorth)
            tau=tau*dt(ng)
#   endif
            IF ((dVdt*dVde).lt.0.0_r8) dVdt=0.0_r8
            IF ((dVdt*(grad(i,Mm(ng))+grad(i+1,Mm(ng)))).gt.0.0_r8) THEN
              dVdx=grad(i  ,Mm(ng))
            ELSE
              dVdx=grad(i+1,Mm(ng))
            END IF
            cff=MAX(dVdx*dVdx+dVde*dVde,eps)
#   ifdef RADIATION_2D
            Cx=MIN(cff,MAX(dVdt*dVdx,-cff))
#   else
            Cx=0.0_r8
#   endif
            Ce=dVdt*dVde
            v(i,Mm(ng)+1,k,nnew)=(cff*v(i,Mm(ng)+1,k,nstp)+             &
     &                            Ce*v(i,Mm(ng),k,nnew)-                &
     &                            MAX(Cx,0.0_r8)*grad(i  ,Mm(ng)+1)-    &
     &                            MIN(Cx,0.0_r8)*grad(i+1,Mm(ng)+1))/   &
     &                           (cff+Ce)
#   ifdef NORTH_M3NUDGING
            v(i,Mm(ng)+1,k,nnew)=v(i,Mm(ng)+1,k,nnew)+                  &
     &                           tau*(BOUNDARY(ng)%v_north(i,k)-        &
     &                                v(i,Mm(ng)+1,k,nstp))
#   endif
#   ifdef MASKING
            v(i,Mm(ng)+1,k,nnew)=v(i,Mm(ng)+1,k,nnew)*                  &
     &                           GRID(ng)%vmask(i,Mm(ng)+1)
#   endif
          END DO
        END DO
#  elif defined NORTH_M3CLAMPED
!
!  Northern edge, clamped boundary condition.
!
        DO k=1,N(ng)
          DO i=Istr,Iend
            v(i,Mm(ng)+1,k,nnew)=BOUNDARY(ng)%v_north(i,k)
#   ifdef MASKING
            v(i,Mm(ng)+1,k,nnew)=v(i,Mm(ng)+1,k,nnew)*                  &
     &                           GRID(ng)%vmask(i,Mm(ng)+1)
#   endif
          END DO
        END DO
#  elif defined NORTH_M3GRADIENT
!
!  Northern edge, gradient boundary condition.
!
        DO k=1,N(ng)
          DO i=Istr,Iend
            v(i,Mm(ng)+1,k,nnew)=v(i,Mm(ng),k,nnew)
#   ifdef MASKING
            v(i,Mm(ng)+1,k,nnew)=v(i,Mm(ng)+1,k,nnew)*                  &
     &                           GRID(ng)%vmask(i,Mm(ng)+1)
#   endif
          END DO
        END DO
#  else
!
!  Northern edge, closed boundary condition.
!
        DO k=1,N(ng)
          DO i=Istr,Iend
            v(i,Mm(ng)+1,k,nnew)=0.0_r8
          END DO
        END DO
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
          v(0,1,k,nnew)=0.5_r8*(v(0,2,k,nnew)+                          &
     &                          v(1,1,k,nnew))
        END DO
      END IF
      IF (SOUTHERN_EDGE.and.EASTERN_EDGE) THEN
        DO k=1,N(ng)
          v(Lm(ng)+1,1,k,nnew)=0.5_r8*(v(Lm(ng)  ,1,k,nnew)+            &
     &                                 v(Lm(ng)+1,2,k,nnew))
        END DO
      END IF
      IF (NORTHERN_EDGE.and.WESTERN_EDGE) THEN
        DO k=1,N(ng)
          v(0,Mm(ng)+1,k,nnew)=0.5_r8*(v(0,Mm(ng)  ,k,nnew)+            &
     &                                 v(1,Mm(ng)+1,k,nnew))
        END DO
      END IF
      IF (NORTHERN_EDGE.and.EASTERN_EDGE) THEN
        DO k=1,N(ng)
          v(Lm(ng)+1,Mm(ng)+1,k,nnew)=0.5_r8*                           &
     &                                (v(Lm(ng)+1,Mm(ng)  ,k,nnew)+     &
     &                                 v(Lm(ng)  ,Mm(ng)+1,k,nnew))
        END DO
      END IF
# endif
      RETURN
      END SUBROUTINE v3dbc_tile
#endif
      END MODULE v3dbc_mod
