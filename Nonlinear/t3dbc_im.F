#include "cppdefs.h"
      MODULE t3dbc_mod
#ifdef SOLVE3D
!
!========================================== Alexander F. Shchepetkin ===
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This subroutine sets lateral boundary conditions for the ITRC-th    !
!  tracer field.                                                       !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC t3dbc_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE t3dbc (ng, tile, indx, itrc)
!***********************************************************************
!
      USE mod_param
      USE mod_ocean
      USE mod_stepping
!
      integer, intent(in) :: ng, tile, indx, itrc

# include "tile.h"
!
      CALL t3dbc_tile (ng, Istr, Iend, Jstr, Jend, itrc,                &
     &                 LBi, UBi, LBj, UBj, N(ng), NT(ng),               &
     &                 indx, nstp(ng),                                  &
     &                 OCEAN(ng)% t)
      RETURN
      END SUBROUTINE t3dbc
!
!***********************************************************************
      SUBROUTINE t3dbc_tile (ng, Istr, Iend, Jstr, Jend, itrc,          &
     &                       LBi, UBi, LBj, UBj, UBk, UBt,              &
     &                       indx, nstp,                                &
     &                       t)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
      USE mod_grid
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr, itrc
      integer, intent(in) :: LBi, UBi, LBj, UBj, UBk, UBt
      integer, intent(in) :: indx, nstp

# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
# else
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,UBk,3,UBt)
# endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j, k

      real(r8), parameter :: eps =1.0E-20_r8

      real(r8) :: Ce, Cx, cff, dTde, dTdt, dTdx, tau

      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: grad

# include "set_bounds.h"

# ifndef EW_PERIODIC
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the western edge.
!-----------------------------------------------------------------------
!
      IF (WESTERN_EDGE) THEN

#  if defined WEST_TRADIATION
!
!  Western edge, implicit upstream radiation condition.
!
        DO k=1,N(ng)
          DO j=Jstr,Jend+1
            grad(0,j)=t(0,j  ,k,nstp,itrc)-                             &
     &                t(0,j-1,k,nstp,itrc)
#   ifdef MASKING
            grad(0,j)=grad(0,j)*GRID(ng)%vmask(0,j)
#   endif
            grad(1,j)=t(1,j  ,k,nstp,itrc)-                             &
     &                t(1,j-1,k,nstp,itrc)
#   ifdef MASKING
            grad(1,j)=grad(1,j)*GRID(ng)%vmask(1,j)
#   endif
          END DO
          DO j=Jstr,Jend
            dTdt=t(1,j,k,nstp,itrc)-t(1,j,k,indx,itrc)
            dTdx=t(1,j,k,indx,itrc)-t(2,j,k,indx,itrc)
#   ifdef WEST_TNUDGING
            tau=Tobc_out(itrc,ng,iwest)
            IF ((dTdt*dTdx).lt.0.0_r8) tau=Tobc_in(itrc,ng,iwest)
            tau=tau*dt(ng)
#   endif
            IF ((dTdt*dTdx).lt.0.0_r8) dTdt=0.0_r8
            IF ((dTdt*(grad(1,j)+grad(1,j+1))).gt.0.0_r8) THEN
              dTde=grad(1,j  )
            ELSE
              dTde=grad(1,j+1)
            END IF
            cff=MAX(dTdx*dTdx+dTde*dTde,eps)
            Cx=dTdt*dTdx
#   ifdef RADIATION_2D
            Ce=MIN(cff,MAX(dTdt*dTde,-cff))
#   else
            Ce=0.0_r8
#   endif
            t(0,j,k,indx,itrc)=(cff*t(0,j,k,nstp,itrc)+                 &
     &                          Cx*t(1,j,k,indx,itrc )-                 &
     &                          MAX(Ce,0.0_r8)*grad(0,j  )-             &
     &                          MIN(Ce,0.0_r8)*grad(0,j+1))/            &
     &                         (cff+Cx)
#   ifdef WEST_TNUDGING
            t(0,j,k,indx,itrc)=t(0,j,k,indx,itrc)+                      &
     &                         tau*(BOUNDARY(ng)%t_west(j,k,itrc)-      &
     &                              t(0,j,k,nstp,itrc))
#   endif
#   ifdef MASKING
            t(0,j,k,indx,itrc)=t(0,j,k,indx,itrc)*                      &
     &                         GRID(ng)%rmask(0,j)
#   endif
          END DO
        END DO
#  elif defined WEST_TCLAMPED
!
!  Western edge, clamped boundary condition.
!
        DO k=1,N(ng)
          DO j=Jstr,Jend
            t(0,j,k,indx,itrc)=BOUNDARY(ng)%t_west(j,k,itrc)
#   ifdef MASKING
            t(0,j,k,indx,itrc)=t(0,j,k,indx,itrc)*                      &
     &                         GRID(ng)%rmask(0,j)
#   endif
          END DO
        END DO
#  elif defined WEST_TGRADIENT
!
!  Western edge, gradient boundary condition.
!
        DO k=1,N(ng)
          DO j=Jstr,Jend
            t(0,j,k,indx,itrc)=t(1,j,k,indx,itrc)
#   ifdef MASKING
            t(0,j,k,indx,itrc)=t(0,j,k,indx,itrc)*                      &
     &                         GRID(ng)%rmask(0,j)
#   endif
          END DO
        END DO
#  else
!
!  Western edge, closed boundary condition.
!
        DO k=1,N(ng)
          DO j=Jstr,Jend
            t(0,j,k,indx,itrc)=t(1,j,k,indx,itrc)
#   ifdef MASKING
            t(0,j,k,indx,itrc)=t(0,j,k,indx,itrc)*                      &
     &                         GRID(ng)%rmask(0,j)
#   endif
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

#  if defined EAST_TRADIATION
!
!  Eastern edge, implicit upstream radiation condition.
!
        DO k=1,N(ng)
          DO j=Jstr,Jend+1
           grad(Lm(ng),j)=t(Lm(ng),j  ,k,nstp,itrc)-                    &
     &                    t(Lm(ng),j-1,k,nstp,itrc)
#   ifdef MASKING
           grad(Lm(ng),j)=grad(Lm(ng),j)*GRID(ng)%vmask(Lm(ng),j)
#   endif
           grad(Lm(ng)+1,j)=t(Lm(ng)+1,j  ,k,nstp,itrc)-                &
     &                      t(Lm(ng)+1,j-1,k,nstp,itrc)
#   ifdef MASKING
           grad(Lm(ng)+1,j)=grad(Lm(ng)+1,j)*GRID(ng)%vmask(Lm(ng)+1,j)
#   endif
          END DO
          DO j=Jstr,Jend
            dTdt=t(Lm(ng),j,k,nstp,itrc)-t(Lm(ng),j,k,indx,itrc)
            dTdx=t(Lm(ng),j,k,indx,itrc)-t(Lm(ng)-1,j,k,indx,itrc)
#   ifdef EAST_TNUDGING
            tau=Tobc_out(itrc,ng,ieast)
            IF ((dTdt*dTdx).lt.0.0_r8) tau=Tobc_in(itrc,ng,ieast)
            tau=tau*dt(ng)
#   endif
            IF ((dTdt*dTdx).lt.0.0_r8) dTdt=0.0_r8
            IF ((dTdt*(grad(Lm(ng),j)+grad(Lm(ng),j+1))).gt.0.0_r8) THEN
              dTde=grad(Lm(ng),j  )
            ELSE
              dTde=grad(Lm(ng),j+1)
            END IF
            cff=MAX(dTdx*dTdx+dTde*dTde,eps)
            Cx=dTdt*dTdx
#   ifdef RADIATION_2D
            Ce=MIN(cff,MAX(dTdt*dTde,-cff))
#   else
            Ce=0.0_r8
#   endif
            t(Lm(ng)+1,j,k,indx,itrc)=(cff*t(Lm(ng)+1,j,k,nstp,itrc)+   &
     &                                 Cx*t(Lm(ng),j,k,indx,itrc)-      &
     &                                 MAX(Ce,0.0_r8)*                  &
     &                                    grad(Lm(ng)+1,j  )-           &
     &                                 MIN(Ce,0.0_r8)*                  &
     &                                    grad(Lm(ng)+1,j+1))/          &
     &                                (cff+Cx)
#   ifdef EAST_TNUDGING
            t(Lm(ng)+1,j,k,indx,itrc)=t(Lm(ng)+1,j,k,indx,itrc)+        &
     &                                tau*                              &
     &                                (BOUNDARY(ng)%t_east(j,k,itrc)-   &
     &                                 t(Lm(ng)+1,j,k,nstp,itrc))
#   endif
#   ifdef MASKING
            t(Lm(ng)+1,j,k,indx,itrc)=t(Lm(ng)+1,j,k,indx,itrc)*        &
     &                                GRID(ng)%rmask(Lm(ng)+1,j)
#   endif
          END DO
        END DO
#  elif defined EAST_TCLAMPED
!
!  Eastern edge, clamped boundary condition.
!
        DO k=1,N(ng)
          DO j=Jstr,Jend
            t(Lm(ng)+1,j,k,indx,itrc)=BOUNDARY(ng)%t_east(j,k,itrc)
#   ifdef MASKING
            t(Lm(ng)+1,j,k,indx,itrc)=t(Lm(ng)+1,j,k,indx,itrc)*        &
     &                                GRID(ng)%rmask(Lm(ng)+1,j)
#   endif
          END DO
        END DO
#  elif defined EAST_TGRADIENT
!
!  Eastern edge, gradient boundary condition.
!
        DO k=1,N(ng)
          DO j=Jstr,Jend
            t(Lm(ng)+1,j,k,indx,itrc)=t(Lm(ng),j,k,indx,itrc)
#   ifdef MASKING
            t(Lm(ng)+1,j,k,indx,itrc)=t(Lm(ng)+1,j,k,indx,itrc)*        &
     &                                GRID(ng)%rmask(Lm(ng)+1,j)
#   endif
          END DO
        END DO
#  else
!
!  Eastern edge, closed boundary condition.
!
        DO k=1,N(ng)
          DO j=Jstr,Jend
            t(Lm(ng)+1,j,k,indx,itrc)=t(Lm(ng),j,k,indx,itrc)
#   ifdef MASKING
            t(Lm(ng)+1,j,k,indx,itrc)=t(Lm(ng)+1,j,k,indx,itrc)*        &
     &                                GRID(ng)%rmask(Lm(ng)+1,j)
#   endif
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

#  if defined SOUTH_TRADIATION
!
!  Southern edge, implicit upstream radiation condition.
!
        DO k=1,N(ng)
          DO i=Istr,Iend+1
            grad(i,1)=t(i  ,1,k,nstp,itrc)-                             &
     &                t(i-1,1,k,nstp,itrc)
#   ifdef MASKING
            grad(i,1)=grad(i,1)*GRID(ng)%umask(i,1)
#   endif
            grad(i,0)=t(i  ,0,k,nstp,itrc)-                             &
     &                t(i-1,0,k,nstp,itrc)
#   ifdef MASKING
            grad(i,0)=grad(i,0)*GRID(ng)%umask(i,0)
#   endif
          END DO
          DO i=Istr,Iend
            dTdt=t(i,1,k,nstp,itrc)-t(i,1,k,indx,itrc)
            dTde=t(i,1,k,indx,itrc)-t(i,2,k,indx,itrc)
#   ifdef SOUTH_TNUDGING
            tau=Tobc_out(itrc,ng,isouth)
            IF ((dTdt*dTde).lt.0.0_r8) tau=Tobc_in(itrc,ng,isouth)
            tau=tau*dt(ng)
#   endif
            IF ((dTdt*dTde).lt.0.0_r8) dTdt=0.0_r8
            IF ((dTdt*(grad(i,1)+grad(i+1,1))).gt.0.0_r8) THEN
              dTdx=grad(i  ,1)
            ELSE
              dTdx=grad(i+1,1)
            END IF
            cff=MAX(dTdx*dTdx+dTde*dTde,eps)
#   ifdef RADIATION_2D
            Cx=MIN(cff,MAX(dTdt*dTdx,-cff))
#   else
            Cx=0.0_r8
#   endif
            Ce=dTdt*dTde
            t(i,0,k,indx,itrc)=(cff*t(i,0,k,nstp,itrc)+                 &
     &                          Ce*t(i,1,k,indx,itrc )-                 &
     &                          MAX(Cx,0.0_r8)*grad(i  ,0)-             &
     &                          MIN(Cx,0.0_r8)*grad(i+1,0))/            &
     &                         (cff+Ce)
#   ifdef SOUTH_TNUDGING
            t(i,0,k,indx,itrc)=t(i,0,k,indx,itrc)+                      &
     &                         tau*(BOUNDARY(ng)%t_south(i,k,itrc)-     &
     &                              t(i,0,k,nstp,itrc))
#   endif
#   ifdef MASKING
            t(i,0,k,indx,itrc)=t(i,0,k,indx,itrc)*                      &
     &                         GRID(ng)%rmask(i,0)
#   endif
          END DO
        END DO
#  elif defined SOUTH_TCLAMPED
!
!  Southern edge, clamped boundary condition.
!
        DO k=1,N(ng)
          DO i=Istr,Iend
            t(i,0,k,indx,itrc)=BOUNDARY(ng)%t_south(i,k,itrc)
#   ifdef MASKING
            t(i,0,k,indx,itrc)=t(i,0,k,indx,itrc)*                      &
     &                         GRID(ng)%rmask(i,0)
#   endif
          END DO
        END DO
#  elif defined SOUTH_TGRADIENT
!
!  Southern edge, gradient boundary condition.
!
        DO k=1,N(ng)
          DO i=Istr,Iend
            t(i,0,k,indx,itrc)=t(i,1,k,indx,itrc)
#   ifdef MASKING
            t(i,0,k,indx,itrc)=t(i,0,k,indx,itrc)*                      &
     &                         GRID(ng)%rmask(i,0)
#   endif
          END DO
        END DO
#  else
!
!  Southern edge, closed boundary condition.
!
        DO k=1,N(ng)
          DO i=Istr,Iend
            t(i,0,k,indx,itrc)=t(i,1,k,indx,itrc)
#   ifdef MASKING
            t(i,0,k,indx,itrc)=t(i,0,k,indx,itrc)*                      &
     &                         GRID(ng)%rmask(i,0)
#   endif
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

#  if defined NORTH_TRADIATION
!
!  Northern edge, implicit upstream radiation condition.
!
        DO k=1,N(ng)
          DO i=Istr,Iend+1
            grad(i,Mm(ng))=t(i  ,Mm(ng),k,nstp,itrc)-                   &
     &                     t(i-1,Mm(ng),k,nstp,itrc)
#   ifdef MASKING
            grad(i,Mm(ng))=grad(i,Mm(ng))*GRID(ng)%umask(i,Mm(ng))
#   endif
            grad(i,Mm(ng)+1)=t(i  ,Mm(ng)+1,k,nstp,itrc)-               &
     &                       t(i-1,Mm(ng)+1,k,nstp,itrc)
#   ifdef MASKING
            grad(i,Mm(ng)+1)=grad(i,Mm(ng)+1)*GRID(ng)%umask(i,Mm(ng)+1)
#   endif
          END DO
          DO i=Istr,Iend
            dTdt=t(i,Mm(ng),k,nstp,itrc)-t(i,Mm(ng),k,indx,itrc)
            dTde=t(i,Mm(ng),k,indx,itrc)-t(i,Mm(ng)-1,k,indx,itrc)
#   ifdef NORTH_TNUDGING
            tau=Tobc_out(itrc,ng,inorth)
            IF ((dTdt*dTde).lt.0.0_r8) tau=Tobc_in(itrc,ng,inorth)
            tau=tau*dt(ng)
#   endif
            IF ((dTdt*dTde).lt.0.0_r8) dTdt=0.0_r8
            IF ((dTdt*(grad(i,Mm(ng))+grad(i+1,Mm(ng)))).gt.0.0_r8) THEN
              dTdx=grad(i  ,Mm(ng))
            ELSE
              dTdx=grad(i+1,Mm(ng))
            END IF
            cff=MAX(dTdx*dTdx+dTde*dTde,eps)
#   ifdef RADIATION_2D
            Cx=MIN(cff,MAX(dTdt*dTdx,-cff))
#   else
            Cx=0.0_r8
#   endif
            Ce=dTdt*dTde
            t(i,Mm(ng)+1,k,indx,itrc)=(cff*t(i,Mm(ng)+1,k,nstp,itrc)+   &
     &                                 Ce*t(i,Mm(ng),k,indx,itrc)-      &
     &                                 MAX(Cx,0.0_r8)*                  &
     &                                    grad(i  ,Mm(ng)+1)-           &
     &                                 MIN(Cx,0.0_r8)*                  &
     &                                    grad(i+1,Mm(ng)+1))/          &
     &                                (cff+Ce)
#   ifdef NORTH_TNUDGING
            t(i,Mm(ng)+1,k,indx,itrc)=t(i,Mm(ng)+1,k,indx,itrc)+        &
     &                                tau*                              &
     &                                (BOUNDARY(ng)%t_north(i,k,itrc)-  &
     &                                 t(i,Mm(ng)+1,k,nstp,itrc))
#   endif
#   ifdef MASKING
            t(i,Mm(ng)+1,k,indx,itrc)=t(i,Mm(ng)+1,k,indx,itrc)*        &
     &                                GRID(ng)%rmask(i,Mm(ng)+1)
#   endif
          END DO
        END DO
#  elif defined NORTH_TCLAMPED
!
!  Northern edge, clamped boundary condition.
!
        DO k=1,N(ng)
          DO i=Istr,Iend
            t(i,Mm(ng)+1,k,indx,itrc)=BOUNDARY(ng)%t_north(i,k,itrc)
#   ifdef MASKING
            t(i,Mm(ng)+1,k,indx,itrc)=t(i,Mm(ng)+1,k,indx,itrc)*        &
     &                                GRID(ng)%rmask(i,Mm(ng)+1)
#   endif
          END DO
        END DO
#  elif defined NORTH_TGRADIENT
!
!  Northern edge, gradient boundary condition.
!
        DO k=1,N(ng)
          DO i=Istr,Iend
            t(i,Mm(ng)+1,k,indx,itrc)=t(i,Mm(ng),k,indx,itrc)
#   ifdef MASKING
            t(i,Mm(ng)+1,k,indx,itrc)=t(i,Mm(ng)+1,k,indx,itrc)*        &
     &                                GRID(ng)%rmask(i,Mm(ng)+1)
#   endif
          END DO
        END DO
#  else
!
!  Northern edge, closed boundary condition.
!
        DO k=1,N(ng)
          DO i=Istr,Iend
            t(i,Mm(ng)+1,k,indx,itrc)=t(i,Mm(ng),k,indx,itrc)
#   ifdef MASKING
            t(i,Mm(ng)+1,k,indx,itrc)=t(i,Mm(ng)+1,k,indx,itrc)*        &
     &                                GRID(ng)%rmask(i,Mm(ng)+1)
#   endif
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
          t(0,0,k,indx,itrc)=0.5_r8*(t(1,0,k,indx,itrc)+                &
     &                               t(0,1,k,indx,itrc))
        END DO
      END IF
      IF (SOUTHERN_EDGE.and.EASTERN_EDGE) THEN
        DO k=1,N(ng)
          t(Lm(ng)+1,0,k,indx,itrc)=0.5_r8*(t(Lm(ng)  ,0,k,indx,itrc)+  &
     &                                      t(Lm(ng)+1,1,k,indx,itrc))
        END DO
      END IF
      IF (NORTHERN_EDGE.and.WESTERN_EDGE) THEN
        DO k=1,N(ng)
          t(0,Mm(ng)+1,k,indx,itrc)=0.5_r8*(t(1,Mm(ng)+1,k,indx,itrc)+  &
     &                                      t(0,Mm(ng)  ,k,indx,itrc))
        END DO
      END IF
      IF (NORTHERN_EDGE.and.EASTERN_EDGE) THEN
        DO k=1,N(ng)
           t(Lm(ng)+1,Mm(ng)+1,k,indx,itrc)=0.5_r8*                     &
     &                              (t(Lm(ng)  ,Mm(ng)+1,k,indx,itrc)+  &
     &                               t(Lm(ng)+1,Mm(ng)  ,k,indx,itrc))
        END DO
      END IF
# endif
      RETURN
      END SUBROUTINE t3dbc_tile
#endif
      END MODULE t3dbc_mod
