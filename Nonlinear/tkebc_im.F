#include "cppdefs.h"
      MODULE tkebc_mod
#if defined SOLVE3D && (defined MY25_MIXING || defined GLS_MIXING)
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This subroutine sets lateral boundary conditions for turbulent      !
!  kinetic energy and turbulent length scale variables associated      !
!  with the Mellor and Yamada or GOTM closures.                        !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC tkebc_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE tkebc (ng, tile, indx)
!***********************************************************************
!
      USE mod_param
      USE mod_mixing
      USE mod_stepping
!
      integer, intent(in) :: ng, tile, indx

# include "tile.h"
!
      CALL tkebc_tile (ng, Istr, Iend, Jstr, Jend,                      &
     &                 LBi, UBi, LBj, UBj, N(ng),                       &
     &                 indx, nstp(ng),                                  &
     &                 MIXING(ng)% gls,                                 &
     &                 MIXING(ng)% tke)
      RETURN
      END SUBROUTINE tkebc
!
!***********************************************************************
      SUBROUTINE tkebc_tile (ng, Istr, Iend, Jstr, Jend,                &
     &                       LBi, UBi, LBj, UBj, UBk,                   &
     &                       indx, nstp,                                &
     &                       gls, tke)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj, UBk
      integer, intent(in) :: indx, nstp

# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: gls(LBi:,LBj:,0:,:)
      real(r8), intent(inout) :: tke(LBi:,LBj:,0:,:)
# else
      real(r8), intent(inout) :: gls(LBi:UBi,LBj:UBj,0:UBk,3)
      real(r8), intent(inout) :: tke(LBi:UBi,LBj:UBj,0:UBk,3)
# endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j, k

      real(r8), parameter :: eps = 1.0e-20_r8

      real(r8) :: Ce, Cx, cff, dKde, dKdt, dKdx

      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: grad
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: gradL

# include "set_bounds.h"

# ifndef EW_PERIODIC
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the western edge.
!-----------------------------------------------------------------------
!
      IF (WESTERN_EDGE) THEN
!
#  if defined WEST_KRADIATION
!
!  Western edge, implicit upstream radiation condition.
!
        DO k=0,N(ng)
          DO j=Jstr,Jend+1
            grad(0,j)=tke(0,j  ,k,nstp)-                                &
     &                tke(0,j-1,k,nstp)
#   ifdef MASKING
            grad(0,j)=grad(0,j)*GRID(ng)%vmask(0,j)
#   endif
            grad(1,j)=tke(1,j  ,k,nstp)-                                &
     &                tke(1,j-1,k,nstp)
#   ifdef MASKING
            grad(1,j)=grad(1,j)*GRID(ng)%vmask(1,j)
#   endif
            gradL(0,j)=gls(0,j  ,k,nstp)-                               &
     &                 gls(0,j-1,k,nstp)
#   ifdef MASKING
            gradL(0,j)=gradL(0,j)*GRID(ng)%vmask(0,j)
#   endif
            gradL(1,j)=gls(1,j  ,k,nstp)-                               &
     &                 gls(1,j-1,k,nstp)
#   ifdef MASKING
            gradL(1,j)=gradL(1,j)*GRID(ng)%vmask(1,j)
#   endif
          END DO
          DO j=Jstr,Jend
            dKdt=tke(1,j,k,nstp)-tke(1,j,k,indx)
            dKdx=tke(1,j,k,indx)-tke(2,j,k,indx)
            IF ((dKdt*dKdx).lt.0.0_r8) dKdt=0.0_r8
            IF ((dKdt*(grad(1,j)+grad(1,j+1))).gt.0.0_r8) THEN
              dKde=grad(1,j  )
            ELSE
              dKde=grad(1,j+1)
            END IF
            cff=MAX(dKdx*dKdx+dKde*dKde,eps)
            Cx=dKdt*dKdx
#   ifdef RADIATION_2D
            Ce=MIN(cff,MAX(dKdt*dKde,-cff))
#   else
            Ce=0.0_r8
#   endif
            tke(0,j,k,indx)=(cff*tke(0,j,k,nstp)+                       &
     &                       Cx*tke(1,j,k,indx)-                        &
     &                       MAX(Ce,0.0_r8)*                            &
     &                          grad(0,j  )-                            &
     &                       MIN(Ce,0.0_r8)*                            &
     &                          grad(0,j+1))/                           &
     &                      (cff+Cx)
#   ifdef MASKING
            tke(0,j,k,indx)=tke(0,j,k,indx)*                            &
     &                      GRID(ng)%rmask(0,j)
#   endif
            dKdt=gls(1,j,k,nstp)-gls(1,j,k,indx)
            dKdx=gls(1,j,k,indx)-gls(2,j,k,indx)
            IF ((dKdt*dKdx).lt.0.0_r8) dKdt=0.0_r8
            IF ((dKdt*(gradL(1,j)+gradL(1,j+1))).gt.0.0_r8) THEN
              dKde=gradL(1,j  )
            ELSE
              dKde=gradL(1,j+1)
            END IF
            cff=MAX(dKdx*dKdx+dKde*dKde,eps)
            Cx=dKdt*dKdx
#   ifdef RADIATION_2D
            Ce=MIN(cff,MAX(dKdt*dKde,-cff))
#   else
            Ce=0.0_r8
#   endif
            gls(0,j,k,indx)=(cff*gls(0,j,k,nstp)+                       &
     &                       Cx*gls(1,j,k,indx)-                        &
     &                       MAX(Ce,0.0_r8)*                            &
     &                          gradL(0,j  )-                           &
     &                       MIN(Ce,0.0_r8)*                            &
     &                          gradL(0,j+1))/                          &
     &                      (cff+Cx)
#   ifdef MASKING
            gls(0,j,k,indx)=gls(0,j,k,indx)*                            &
     &                      GRID(ng)%rmask(0,j)
#   endif
          END DO
        END DO
#  elif defined WEST_KGRADIENT
!
!  Western edge, gradient boundary condition.
!
        DO k=0,N(ng)
          DO j=Jstr,Jend
            tke(0,j,k,indx)=tke(1,j,k,indx)
#   ifdef MASKING
            tke(0,j,k,indx)=tke(0,j,k,indx)*                            &
     &                      GRID(ng)%rmask(0,j)
#   endif
            gls(0,j,k,indx)=gls(1,j,k,indx)
#   ifdef MASKING
            gls(0,j,k,indx)=gls(0,j,k,indx)*                            &
     &                      GRID(ng)%rmask(0,j)
#   endif
          END DO
        END DO
#  else
!
!  Western edge, closed boundary condition.
!
        DO k=0,N(ng)
          DO j=Jstr,Jend
            tke(0,j,k,indx)=tke(1,j,k,indx)
#   ifdef MASKING
            tke(0,j,k,indx)=tke(0,j,k,indx)*                            &
     &                      GRID(ng)%rmask(0,j)
#   endif
            gls(0,j,k,indx)=gls(1,j,k,indx)
#   ifdef MASKING
            gls(0,j,k,indx)=gls(0,j,k,indx)*                            &
     &                      GRID(ng)%rmask(0,j)
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
!
#  if defined EAST_KRADIATION
!
!  Eastern edge, implicit upstream radiation condition.
!
        DO k=0,N(ng)
          DO j=Jstr,Jend+1
           grad(Lm(ng),j)=tke(Lm(ng),j  ,k,nstp)-                       &
     &                    tke(Lm(ng),j-1,k,nstp)
#   ifdef MASKING
           grad(Lm(ng),j)=grad(Lm(ng),j)*                               &
     &                    GRID(ng)%vmask(Lm(ng),j)
#   endif
           grad(Lm(ng)+1,j)=tke(Lm(ng)+1,j  ,k,nstp)-                   &
     &                      tke(Lm(ng)+1,j-1,k,nstp)
#   ifdef MASKING
           grad(Lm(ng)+1,j)=grad(Lm(ng)+1,j)*                           &
     &                      GRID(ng)%vmask(Lm(ng)+1,j)
#   endif
           gradL(Lm(ng),j)=gls(Lm(ng),j  ,k,nstp)-                      &
     &                     gls(Lm(ng),j-1,k,nstp)
#   ifdef MASKING
           gradL(Lm(ng),j)=gradL(Lm(ng),j)*                             &
     &                     GRID(ng)%vmask(Lm(ng),j)
#   endif
           gradL(Lm(ng)+1,j)=gls(Lm(ng)+1,j  ,k,nstp)-                  &
     &                       gls(Lm(ng)+1,j-1,k,nstp)
#   ifdef MASKING
           gradL(Lm(ng)+1,j)=gradL(Lm(ng)+1,j)*                         &
     &                       GRID(ng)%vmask(Lm(ng)+1,j)
#   endif
          END DO
          DO j=Jstr,Jend
            dKdt=tke(Lm(ng),j,k,nstp)-tke(Lm(ng),j,k,indx)
            dKdx=tke(Lm(ng),j,k,indx)-tke(Lm(ng)-1,j,k,indx)
            IF ((dKdt*dKdx).lt.0.0_r8) dKdt=0.0_r8
            IF ((dKdt*(grad(Lm(ng),j)+grad(Lm(ng),j+1))).gt.0.0_r8) THEN
              dKde=grad(Lm(ng),j  )
            ELSE
              dKde=grad(Lm(ng),j+1)
            END IF
            cff=MAX(dKdx*dKdx+dKde*dKde,eps)
            Cx=dKdt*dKdx
#   ifdef RADIATION_2D
            Ce=MIN(cff,MAX(dKdt*dKde,-cff))
#   else
            Ce=0.0_r8
#   endif
            tke(Lm(ng)+1,j,k,indx)=(cff*tke(Lm(ng)+1,j,k,nstp)+         &
     &                              Cx*tke(Lm(ng),j,k,indx)-            &
     &                              MAX(Ce,0.0_r8)*                     &
     &                                 grad(Lm(ng)+1,j  )-              &
     &                              MIN(Ce,0.0_r8)*                     &
     &                                 grad(Lm(ng)+1,j+1))/             &
     &                             (cff+Cx)
#   ifdef MASKING
            tke(Lm(ng)+1,j,k,indx)=tke(Lm(ng)+1,j,k,indx)*              &
     &                             GRID(ng)%rmask(Lm(ng)+1,j)
#   endif
            dKdt=gls(Lm(ng),j,k,nstp)-gls(Lm(ng),j,k,indx)
            dKdx=gls(Lm(ng),j,k,indx)-gls(Lm(ng)-1,j,k,indx)
            IF ((dKdt*dKdx).lt.0.0_r8) dKdt=0.0_r8
            IF ((dKdt*(gradL(Lm(ng),j)+                                 &
     &                 gradL(Lm(ng),j+1))).gt.0.0_r8) THEN
              dKde=gradL(Lm(ng),j  )
            ELSE
              dKde=gradL(Lm(ng),j+1)
            END IF
            cff=MAX(dKdx*dKdx+dKde*dKde,eps)
            Cx=dKdt*dKdx
#   ifdef RADIATION_2D
            Ce=MIN(cff,MAX(dKdt*dKde,-cff))
#   else
            Ce=0.0_r8
#   endif
            gls(Lm(ng)+1,j,k,indx)=(cff*gls(Lm(ng)+1,j,k,nstp)+         &
     &                              Cx*gls(Lm(ng),j,k,indx)-            &
     &                              MAX(Ce,0.0_r8)*                     &
     &                                 gradL(Lm(ng)+1,j  )-             &
     &                              MIN(Ce,0.0_r8)*                     &
     &                                 gradL(Lm(ng)+1,j+1))/            &
     &                             (cff+Cx)
#   ifdef MASKING
            gls(Lm(ng)+1,j,k,indx)=gls(Lm(ng)+1,j,k,indx)*              &
     &                             GRID(ng)%rmask(Lm(ng)+1,j)
#   endif
          END DO
        END DO
#  elif defined EAST_KGRADIENT
!
!  Eastern edge, gradient boundary condition.
!
        DO k=0,N(ng)
          DO j=Jstr,Jend
            tke(Lm(ng)+1,j,k,indx)=tke(Lm(ng),j,k,indx)
#   ifdef MASKING
            tke(Lm(ng)+1,j,k,indx)=tke(Lm(ng)+1,j,k,indx)*              &
     &                             GRID(ng)%rmask(Lm(ng)+1,j)
#   endif
            gls(Lm(ng)+1,j,k,indx)=gls(Lm(ng),j,k,indx)
#   ifdef MASKING
            gls(Lm(ng)+1,j,k,indx)=gls(Lm(ng)+1,j,k,indx)*              &
     &                             GRID(ng)%rmask(Lm(ng)+1,j)
#   endif
          END DO
        END DO
#  else
!
!  Eastern edge, closed boundary condition.
!
        DO k=0,N(ng)
          DO j=Jstr,Jend
            tke(Lm(ng)+1,j,k,indx)=tke(Lm(ng),j,k,indx)
#   ifdef MASKING
            tke(Lm(ng)+1,j,k,indx)=tke(Lm(ng)+1,j,k,indx)*              &
     &                             GRID(ng)%rmask(Lm(ng)+1,j)
#   endif
            gls(Lm(ng)+1,j,k,indx)=gls(Lm(ng),j,k,indx)
#   ifdef MASKING
            gls(Lm(ng)+1,j,k,indx)=gls(Lm(ng)+1,j,k,indx)*              &
     &                             GRID(ng)%rmask(Lm(ng)+1,j)
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
!
#  if defined SOUTH_KRADIATION
!
!  Southern edge, implicit upstream radiation condition.
!
        DO k=0,N(ng)
          DO i=Istr,Iend+1
            grad(i,1)=tke(i  ,1,k,nstp)-                                &
     &                tke(i-1,1,k,nstp)
#   ifdef MASKING
            grad(i,1)=grad(i,1)*GRID(ng)%umask(i,1)
#   endif
            grad(i,0)=tke(i  ,0,k,nstp)-                                &
     &                tke(i-1,0,k,nstp)
#   ifdef MASKING
            grad(i,0)=grad(i,0)*GRID(ng)%umask(i,0)
#   endif
            gradL(i,1)=gls(i  ,1,k,nstp)-                               &
     &                 gls(i-1,1,k,nstp)
#   ifdef MASKING
            gradL(i,1)=gradL(i,1)*GRID(ng)%umask(i,1)
#   endif
            gradL(i,0)=gls(i  ,0,k,nstp)-                               &
     &                 gls(i-1,0,k,nstp)
#   ifdef MASKING
            gradL(i,0)=gradL(i,0)*GRID(ng)%umask(i,0)
#   endif
          END DO
          DO i=Istr,Iend
            dKdt=tke(i,1,k,nstp)-tke(i,1,k,indx)
            dKde=tke(i,1,k,indx)-tke(i,2,k,indx)
            IF ((dKdt*dKde).lt.0.0_r8) dKdt=0.0_r8
            IF ((dKdt*(grad(i,1)+grad(i+1,1))).gt.0.0_r8) THEN
              dKdx=grad(i  ,1)
            ELSE
              dKdx=grad(i+1,1)
            END IF
            cff=MAX(dKdx*dKdx+dKde*dKde, eps)
#   ifdef RADIATION_2D
            Cx=MIN(cff,MAX(dKdt*dKdx,-cff))
#   else
            Cx=0.0_r8
#   endif
            Ce=dKdt*dKde
            tke(i,0,k,indx)=(cff*tke(i,0,k,nstp)+                       &
     &                       Ce*tke(i,1,k,indx)-                        &
     &                       MAX(Cx,0.0_r8)*                            &
     &                          grad(i  ,0)-                            &
     &                       MIN(Cx,0.0_r8)*                            &
     &                          grad(i+1,0))/                           &
     &                      (cff+Ce)
#   ifdef MASKING
            tke(i,0,k,indx)=tke(i,0,k,indx)*                            &
     &                      GRID(ng)%rmask(i,0)
#   endif
            dKdt=gls(i,1,k,nstp)-gls(i,1,k,indx)
            dKde=gls(i,1,k,indx)-gls(i,2,k,indx)
            IF ((dKdt*dKde).lt.0.0_r8) dKdt=0.0_r8
            IF ((dKdt*(gradL(i,1)+gradL(i+1,1))).gt.0.0_r8) THEN
              dKdx=gradL(i  ,1)
            ELSE
              dKdx=gradL(i+1,1)
            END IF
            cff=MAX(dKdx*dKdx+dKde*dKde,eps)
#   ifdef RADIATION_2D
            Cx=MIN(cff,MAX(dKdt*dKdx,-cff))
#   else
            Cx=0.0_r8
#   endif
            Ce=dKdt*dKde
            gls(i,0,k,indx)=(cff*gls(i,0,k,nstp)+                       &
     &                       Ce*gls(i,1,k,indx)-                        &
     &                       MAX(Cx,0.0_r8)*                            &
     &                          gradL(i  ,0)-                           &
     &                       MIN(Cx,0.0_r8)*                            &
     &                          gradL(i+1,0))/                          &
     &                      (cff+Ce)
#   ifdef MASKING
            gls(i,0,k,indx)=gls(i,0,k,indx)*                            &
     &                      GRID(ng)%rmask(i,0)
#   endif
          END DO
        END DO
#  elif defined SOUTH_KGRADIENT
!
!  Southern edge, gradient boundary condition.
!
        DO k=0,N(ng)
          DO i=Istr,Iend
            tke(i,0,k,indx)=tke(i,1,k,indx)
#   ifdef MASKING
            tke(i,0,k,indx)=tke(i,0,k,indx)*                            &
     &                      GRID(ng)%rmask(i,0)
#   endif
            gls(i,0,k,indx)=gls(i,1,k,indx)
#   ifdef MASKING
            gls(i,0,k,indx)=gls(i,0,k,indx)*                            &
     &                      GRID(ng)%rmask(i,0)
#   endif
          END DO
        END DO
#  else
!
!  Southern edge, closed boundary condition.
!
        DO k=0,N(ng)
          DO i=Istr,Iend
            tke(i,0,k,indx)=tke(i,1,k,indx)
#   ifdef MASKING
            tke(i,0,k,indx)=tke(i,0,k,indx)*                            &
     &                      GRID(ng)%rmask(i,0)
#   endif
            gls(i,0,k,indx)=gls(i,1,k,indx)
#   ifdef MASKING
            gls(i,0,k,indx)=gls(i,0,k,indx)*                            &
     &                      GRID(ng)%rmask(i,0)
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
!
#  if defined NORTH_KRADIATION
!
!  Northern edge, implicit upstream radiation condition.
!
        DO k=0,N(ng)
          DO i=Istr,Iend+1
            grad(i,Mm(ng))=tke(i  ,Mm(ng),k,nstp)-                      &
     &                     tke(i-1,Mm(ng),k,nstp)
#   ifdef MASKING
            grad(i,Mm(ng))=grad(i,Mm(ng))*                              &
     &                     GRID(ng)%umask(i,Mm(ng))
#   endif
            grad(i,Mm(ng)+1)=tke(i  ,Mm(ng)+1,k,nstp)-                  &
     &                       tke(i-1,Mm(ng)+1,k,nstp)
#   ifdef MASKING
            grad(i,Mm(ng)+1)=grad(i,Mm(ng)+1)*                          &
     &                       GRID(ng)%umask(i,Mm(ng)+1)
#   endif
            gradL(i,Mm(ng))=gls(i  ,Mm(ng),k,nstp)-                     &
     &                      gls(i-1,Mm(ng),k,nstp)
#   ifdef MASKING
            gradL(i,Mm(ng))=gradL(i,Mm(ng))*                            &
     &                      GRID(ng)%umask(i,Mm(ng))
#   endif
            gradL(i,Mm(ng)+1)=gls(i  ,Mm(ng)+1,k,nstp)-                 &
     &                        gls(i-1,Mm(ng)+1,k,nstp)
#   ifdef MASKING
            gradL(i,Mm(ng)+1)=gradL(i,Mm(ng)+1)*                        &
     &                        GRID(ng)%umask(i,Mm(ng)+1)
#   endif
          END DO
          DO i=Istr,Iend
            dKdt=tke(i,Mm(ng),k,nstp)-tke(i,Mm(ng),k,indx)
            dKde=tke(i,Mm(ng),k,indx)-tke(i,Mm(ng)-1,k,indx)
            IF ((dKdt*dKde).lt.0.0_r8) dKdt=0.0_r8
            IF ((dKdt*(grad(i,Mm(ng))+grad(i+1,Mm(ng)))).gt.0.0_r8) THEN
              dKdx=grad(i  ,Mm(ng))
            ELSE
              dKdx=grad(i+1,Mm(ng))
            END IF
            cff=MAX(dKdx*dKdx+dKde*dKde,eps)
#   ifdef RADIATION_2D
            Cx=MIN(cff,MAX(dKdt*dKdx,-cff))
#   else
            Cx=0.0_r8
#   endif
            Ce=dKdt*dKde
            tke(i,Mm(ng)+1,k,indx)=(cff*tke(i,Mm(ng)+1,k,nstp)+         &
     &                              Ce*tke(i,Mm(ng),k,indx)-            &
     &                              MAX(Cx,0.0_r8)*                     &
     &                                 grad(i  ,Mm(ng)+1)-              &
     &                              MIN(Cx,0.0_r8)*                     &
     &                                 grad(i+1,Mm(ng)+1))/             &
     &                             (cff+Ce)
#   ifdef MASKING
            tke(i,Mm(ng)+1,k,indx)=tke(i,Mm(ng)+1,k,indx)*              &
     &                             GRID(ng)%rmask(i,Mm(ng)+1)
#   endif
            dKdt=gls(i,Mm(ng),k,nstp)-gls(i,Mm(ng),k,indx)
            dKde=gls(i,Mm(ng),k,indx)-gls(i,Mm(ng)-1,k,indx)
            IF ((dKdt*dKde).lt.0.0_r8) dKdt=0.0_r8
            IF ((dKdt*(gradL(i  ,Mm(ng))+                               &
     &                 gradL(i+1,Mm(ng)))).gt.0.0_r8) THEN
              dKdx=gradL(i  ,Mm(ng))
            ELSE
              dKdx=gradL(i+1,Mm(ng))
            END IF
            cff=MAX(dKdx*dKdx+dKde*dKde,eps)
#   ifdef RADIATION_2D
            Cx=MIN(cff,MAX(dKdt*dKdx,-cff))
#   else
            Cx=0.0_r8
#   endif
            Ce=dKdt*dKde
            gls(i,Mm(ng)+1,k,indx)=(cff*gls(i,Mm(ng)+1,k,nstp)+         &
     &                              Ce*gls(i,Mm(ng),k,indx)-            &
     &                              MAX(Cx,0.0_r8)*                     &
     &                                 gradL(i  ,Mm(ng)+1)-             &
     &                              MIN(Cx,0.0_r8)*                     &
     &                                 gradL(i+1,Mm(ng)+1))/            &
     &                             (cff+Ce)
#   ifdef MASKING
            gls(i,Mm(ng)+1,k,indx)=gls(i,Mm(ng)+1,k,indx)*              &
     &                             GRID(ng)%rmask(i,Mm(ng)+1)
#   endif
          END DO
        END DO
#  elif defined NORTH_KGRADIENT
!
!  Northern edge, gradient boundary condition.
!
        DO k=0,N(ng)
          DO i=Istr,Iend
            tke(i,Mm(ng)+1,k,indx)=tke(i,Mm(ng),k,indx)
#   ifdef MASKING
            tke(i,Mm(ng)+1,k,indx)=tke(i,Mm(ng)+1,k,indx)*              &
     &                             GRID(ng)%rmask(i,Mm(ng)+1)
#   endif
            gls(i,Mm(ng)+1,k,indx)=gls(i,Mm(ng),k,indx)
#   ifdef MASKING
            gls(i,Mm(ng)+1,k,indx)=gls(i,Mm(ng)+1,k,indx)*              &
     &                             GRID(ng)%rmask(i,Mm(ng)+1)
#   endif
          END DO
        END DO
#  else
!
!  Northern edge, closed boundary condition.
!
        DO k=0,N(ng)
          DO i=Istr,Iend
            tke(i,Mm(ng)+1,k,indx)=tke(i,Mm(ng),k,indx)
#   ifdef MASKING
            tke(i,Mm(ng)+1,k,indx)=tke(i,Mm(ng)+1,k,indx)*              &
     &                             GRID(ng)%rmask(i,Mm(ng)+1)
#   endif
            gls(i,Mm(ng)+1,k,indx)=gls(i,Mm(ng),k,indx)
#   ifdef MASKING
            gls(i,Mm(ng)+1,k,indx)=gls(i,Mm(ng)+1,k,indx)*              &
     &                             GRID(ng)%rmask(i,Mm(ng)+1)
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
        DO k=0,N(ng)
          tke(0,0,k,indx)=0.5_r8*(tke(1,0,k,indx)+                      &
     &                            tke(0,1,k,indx))
          gls(0,0,k,indx)=0.5_r8*(gls(1,0,k,indx)+                      &
     &                            gls(0,1,k,indx))
        END DO
      END IF
      IF (SOUTHERN_EDGE.and.EASTERN_EDGE) THEN
        DO k=0,N(ng)
          tke(Lm(ng)+1,0,k,indx)=0.5_r8*(tke(Lm(ng)  ,0,k,indx)+        &
     &                                   tke(Lm(ng)+1,1,k,indx))
          gls(Lm(ng)+1,0,k,indx)=0.5_r8*(gls(Lm(ng)  ,0,k,indx)+        &
     &                                   gls(Lm(ng)+1,1,k,indx))
        END DO
      END IF
      IF (NORTHERN_EDGE.and.WESTERN_EDGE) THEN
        DO k=0,N(ng)
          tke(0,Mm(ng)+1,k,indx)=0.5_r8*(tke(1,Mm(ng)+1,k,indx)+        &
     &                                   tke(0,Mm(ng)  ,k,indx))
          gls(0,Mm(ng)+1,k,indx)=0.5_r8*(gls(1,Mm(ng)+1,k,indx)+        &
     &                                   gls(0,Mm(ng)  ,k,indx))
        END DO
      END IF
      IF (NORTHERN_EDGE.and.EASTERN_EDGE) THEN
        DO k=0,N(ng)
          tke(Lm(ng)+1,Mm(ng)+1,k,indx)=0.5_r8*                         &
     &                                  (tke(Lm(ng)  ,Mm(ng)+1,k,indx)+ &
     &                                   tke(Lm(ng)+1,Mm(ng)  ,k,indx))
          gls(Lm(ng)+1,Mm(ng)+1,k,indx)=0.5_r8*                         &
     &                                  (gls(Lm(ng)  ,Mm(ng)+1,k,indx)+ &
     &                                   gls(Lm(ng)+1,Mm(ng)  ,k,indx))
        END DO
      END IF
# endif
      RETURN
      END SUBROUTINE tkebc_tile
#endif
      END MODULE tkebc_mod
