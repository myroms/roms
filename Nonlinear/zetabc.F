#include "cppdefs.h"
      MODULE zetabc_mod
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine sets lateral boundary conditions for free-surface.     !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC zetabc_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE zetabc (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_ocean
      USE mod_stepping
!
      integer, intent(in) :: ng, tile

# include "tile.h"
!
      CALL zetabc_tile (ng, Istr, Iend, Jstr, Jend,                     &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  krhs(ng), kstp(ng), knew(ng),                   &
     &                  OCEAN(ng) % zeta)
      RETURN
      END SUBROUTINE zetabc
!
!***********************************************************************
      SUBROUTINE zetabc_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        krhs, kstp, knew,                         &
     &                        zeta)
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
      real(r8), intent(inout) :: zeta(LBi:,LBj:,:)
#else
      real(r8), intent(inout) :: zeta(LBi:UBi,LBj:UBj,3)
#endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j, know

      real(r8), parameter :: eps =1.0E-20_r8

      real(r8) :: Ce, Cx, cff, dt2d, dZde, dZdt, dZdx, tau

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

# if defined WEST_FSRADIATION
!
!  Western edge, implicit upstream radiation condition.
!
        DO j=Jstr,Jend+1
          grad(0,j)=zeta(0,j  ,know)-                                   &
     &              zeta(0,j-1,know)
#  ifdef MASKING
          grad(0,j)=grad(0,j)*GRID(ng)%vmask(0,j)
#  endif
          grad(1,j)=zeta(1,j  ,know)-                                   &
     &              zeta(1,j-1,know)
#  ifdef MASKING
          grad(1,j)=grad(1,j)*GRID(ng)%vmask(1,j)
#  endif
        END DO
        DO j=Jstr,Jend
          dZdt=zeta(1,j,know)-zeta(1,j,knew)
          dZdx=zeta(1,j,knew)-zeta(2,j,knew)
#  ifdef WEST_FSNUDGING
          tau=FSobc_out(ng,iwest)
          IF ((dZdt*dZdx).lt.0.0_r8) tau=FSobc_in(ng,iwest)
          tau=tau*dt2d
#  endif
          IF ((dZdt*dZdx).lt.0.0_r8) dZdt=0.0_r8
          IF ((dZdt*(grad(1,j)+grad(1,j+1))).gt.0.0_r8) THEN
            dZde=grad(1,j  )
          ELSE
            dZde=grad(1,j+1)
          END IF
          cff=MAX(dZdx*dZdx+dZde*dZde,eps)
          Cx=dZdt*dZdx
#  ifdef RADIATION_2D
          Ce=MIN(cff,MAX(dZdt*dZde,-cff))
#  else
          Ce=0.0_r8
#  endif
          zeta(0,j,knew)=(cff*zeta(0,j,know)+                           &
     &                    Cx*zeta(1,j,knew)-                            &
     &                    MAX(Ce,0.0_r8)*grad(0,j  )-                   &
     &                    MIN(Ce,0.0_r8)*grad(0,j+1))/                  &
     &                   (cff+Cx)
#  ifdef WEST_FSNUDGING
          zeta(0,j,knew)=zeta(0,j,knew)+                                &
     &                   tau*(BOUNDARY(ng)%zeta_west(j)-                &
     &                        zeta(0,j,know))
#  endif
#  ifdef MASKING
          zeta(0,j,knew)=zeta(0,j,knew)*                                &
     &                   GRID(ng)%rmask(0,j)
#  endif
        END DO
# elif defined WEST_FSCHAPMAN
!
!  Western edge, Chapman boundary condition.
!
        DO j=Jstr,Jend
          Cx=dt2d*GRID(ng)%pm(1,j)*SQRT(g*GRID(ng)%h(1,j))
          zeta(0,j,knew)=(zeta(0,j,know)+                               &
     &                    Cx*zeta(1,j,knew))/                           &
     &                   (1.0_r8+Cx)
#  ifdef MASKING
          zeta(0,j,knew)=zeta(0,j,knew)*                                &
     &                   GRID(ng)%rmask(0,j)
#  endif
        END DO
# elif defined WEST_FSCLAMPED
!
!  Western edge, clamped boundary condition.
!
        DO j=Jstr,Jend
          zeta(0,j,knew)=BOUNDARY(ng)%zeta_west(j)
#  ifdef MASKING
          zeta(0,j,knew)=zeta(0,j,knew)*                                &
     &                   GRID(ng)%rmask(0,j)
#  endif
        END DO
# elif defined WEST_FSGRADIENT
!
!  Western edge, gradient boundary condition.
!
        DO j=Jstr,Jend
          zeta(0,j,knew)=zeta(1,j,knew)
#  ifdef MASKING
          zeta(0,j,knew)=zeta(0,j,knew)*                                &
     &                   GRID(ng)%rmask(0,j)
#  endif
        END DO
# else
!
!  Western edge, closed boundary condition.
!
        DO j=Jstr,Jend
          zeta(0,j,knew)=zeta(1,j,knew)
#  ifdef MASKING
          zeta(0,j,knew)=zeta(0,j,knew)*                                &
     &                   GRID(ng)%rmask(0,j)
#  endif
        END DO
# endif
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the eastern edge.
!-----------------------------------------------------------------------
!
      IF (EASTERN_EDGE) THEN

# if defined EAST_FSRADIATION
!
!  Eastern edge, implicit upstream radiation condition.
!
        DO j=Jstr,Jend+1
          grad(Lm(ng),j)=zeta(Lm(ng),j  ,know)-                         &
     &                   zeta(Lm(ng),j-1,know)
#  ifdef MASKING
          grad(Lm(ng),j)=grad(Lm(ng),j)*GRID(ng)%vmask(Lm(ng),j)
#  endif
          grad(Lm(ng)+1,j)=zeta(Lm(ng)+1,j  ,know)-                     &
     &                     zeta(Lm(ng)+1,j-1,know)
#  ifdef MASKING
          grad(Lm(ng)+1,j)=grad(Lm(ng)+1,j)*GRID(ng)%vmask(Lm(ng)+1,j)
#  endif
        END DO
        DO j=Jstr,Jend
          dZdt=zeta(Lm(ng),j,know)-zeta(Lm(ng),j,knew)
          dZdx=zeta(Lm(ng),j,knew)-zeta(Lm(ng)-1,j,knew)
#  ifdef EAST_FSNUDGING
          tau=FSobc_out(ng,ieast)
          IF ((dZdt*dZdx).lt.0.0_r8) tau=FSobc_in(ng,ieast)
          tau=tau*dt2d
#  endif
          IF ((dZdt*dZdx).lt.0.0_r8) dZdt=0.0_r8
          IF ((dZdt*(grad(Lm(ng),j)+grad(Lm(ng),j+1))).gt.0.0_r8) THEN
            dZde=grad(Lm(ng),j  )
          ELSE
            dZde=grad(Lm(ng),j+1)
          END IF
          cff=MAX(dZdx*dZdx+dZde*dZde,eps)
          Cx=dZdt*dZdx
#  ifdef RADIATION_2D
          Ce=MIN(cff,MAX(dZdt*dZde,-cff))
#  else
          Ce=0.0_r8
#  endif
          zeta(Lm(ng)+1,j,knew)=(cff*zeta(Lm(ng)+1,j,know)+             &
     &                           Cx*zeta(Lm(ng),j,knew)-                &
     &                           MAX(Ce,0.0_r8)*grad(Lm(ng)+1,j  )-     &
     &                           MIN(Ce,0.0_r8)*grad(Lm(ng)+1,j+1))/    &
     &                          (cff+Cx)
#  ifdef EAST_FSNUDGING
          zeta(Lm(ng)+1,j,knew)=zeta(Lm(ng)+1,j,knew)+                  &
     &                          tau*(BOUNDARY(ng)%zeta_east(j)-         &
     &                               zeta(Lm(ng)+1,j,know))
#  endif
#  ifdef MASKING
          zeta(Lm(ng)+1,j,knew)=zeta(Lm(ng)+1,j,knew)*                  &
     &                          GRID(ng)%rmask(Lm(ng)+1,j)
#  endif
        END DO
# elif defined EAST_FSCHAPMAN
!
!  Eastern edge, Chapman boundary condition.
!
        DO j=Jstr,Jend
          Cx=dt2d*GRID(ng)%pm(Lm(ng),j)*SQRT(g*GRID(ng)%h(Lm(ng),j))
          zeta(Lm(ng)+1,j,knew)=(zeta(Lm(ng)+1,j,know)+                 &
     &                           Cx*zeta(Lm(ng),j,knew))/               &
     &                          (1.0_r8+Cx)
#  ifdef MASKING
          zeta(Lm(ng)+1,j,knew)=zeta(Lm(ng)+1,j,knew)*                  &
     &                          GRID(ng)%rmask(Lm(ng)+1,j)
#  endif
        END DO
# elif defined EAST_FSCLAMPED
!
!  Eastern edge, clamped boundary condition.
!
        DO j=Jstr,Jend
          zeta(Lm(ng)+1,j,knew)=BOUNDARY(ng)%zeta_east(j)
#  ifdef MASKING
          zeta(Lm(ng)+1,j,knew)=zeta(Lm(ng)+1,j,knew)*                  &
     &                          GRID(ng)%rmask(Lm(ng)+1,j)
#  endif
        END DO
# elif defined EAST_FSGRADIENT
!
!  Eastern edge, gradient boundary condition.
!
        DO j=Jstr,Jend
          zeta(Lm(ng)+1,j,knew)=zeta(Lm(ng),j,knew)
#  ifdef MASKING
          zeta(Lm(ng)+1,j,knew)=zeta(Lm(ng)+1,j,knew)*                  &
     &                          GRID(ng)%rmask(Lm(ng)+1,j)
#  endif
        END DO
# else
!
!  Eastern edge, closed boundary condition.
!
        DO j=Jstr,Jend
          zeta(Lm(ng)+1,j,knew)=zeta(Lm(ng),j,knew)
#  ifdef MASKING
          zeta(Lm(ng)+1,j,knew)=zeta(Lm(ng)+1,j,knew)*                  &
     &                          GRID(ng)%rmask(Lm(ng)+1,j)
#  endif
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

# if defined SOUTH_FSRADIATION
!
!  Southern edge, implicit upstream radiation condition.
!
        DO i=Istr,Iend+1
          grad(i,1)=zeta(i  ,1,know)-                                   &
     &              zeta(i-1,1,know)
#  ifdef MASKING
          grad(i,1)=grad(i,1)*GRID(ng)%umask(i,1)
#  endif
          grad(i,0)=zeta(i  ,0,know)-                                   &
     &              zeta(i-1,0,know)
#  ifdef MASKING
          grad(i,0)=grad(i,0)*GRID(ng)%umask(i,0)
#  endif
        END DO
        DO i=Istr,Iend
          dZdt=zeta(i,1,know)-zeta(i,1,knew)
          dZde=zeta(i,1,knew)-zeta(i,2,knew)
#  ifdef SOUTH_FSNUDGING
          tau=FSobc_out(ng,isouth)
          IF ((dZdt*dZde).lt.0.0_r8) tau=FSobc_in(ng,isouth)
          tau=tau*dt2d
#  endif
          IF ((dZdt*dZde).lt.0.0_r8) dZdt=0.0_r8
          IF ((dZdt*(grad(i,1)+grad(i+1,1))).gt.0.0_r8) THEN
            dZdx=grad(i  ,1)
          ELSE
            dZdx=grad(i+1,1)
          END IF
          cff=MAX(dZdx*dZdx+dZde*dZde,eps)
#  ifdef RADIATION_2D
          Cx=MIN(cff,MAX(dZdt*dZdx,-cff))
#  else
          Cx=0.0_r8
#  endif
          Ce=dZdt*dZde
          zeta(i,0,knew)=(cff*zeta(i,0,know)+                           &
     &                    Ce*zeta(i,1,knew)-                            &
     &                    MAX(Cx,0.0_r8)*grad(i  ,0)-                   &
     &                    MIN(Cx,0.0_r8)*grad(i+1,0))/                  &
     &                   (cff+Ce)
#  ifdef SOUTH_FSNUDGING
          zeta(i,0,knew)=zeta(i,0,knew)+                                &
     &                   tau*(BOUNDARY(ng)%zeta_south(i)-               &
     &                        zeta(i,0,know))
#  endif
#  ifdef MASKING
          zeta(i,0,knew)=zeta(i,0,knew)*                                &
     &                   GRID(ng)%rmask(i,0)
#  endif
        END DO
# elif defined SOUTH_FSCHAPMAN
!
!  Southern edge, Chapman boundary condition.
!
        DO i=Istr,Iend
          Ce=dt2d*GRID(ng)%pn(i,1)*SQRT(g*GRID(ng)%h(i,1))
          zeta(i,0,knew)=(zeta(i,0,know)+                               &
     &                    Ce*zeta(i,1,knew))/                           &
     &                   (1.0_r8+Ce)
#  ifdef MASKING
          zeta(i,0,knew)=zeta(i,0,knew)*                                &
     &                   GRID(ng)%rmask(i,0)
#  endif
        END DO
# elif defined SOUTH_FSCLAMPED
!
!  Southern edge, clamped boundary condition.
!
        DO i=Istr,Iend
          zeta(i,0,knew)=BOUNDARY(ng)%zeta_south(i)
#  ifdef MASKING
          zeta(i,0,knew)=zeta(i,0,knew)*                                &
     &                   GRID(ng)%rmask(i,0)
#  endif
        END DO
# elif defined SOUTH_FSGRADIENT
!
!  Southern edge, gradient boundary condition.
!
        DO i=Istr,Iend
          zeta(i,0,knew)=zeta(i,1,knew)
#  ifdef MASKING
          zeta(i,0,knew)=zeta(i,0,knew)*                                &
     &                   GRID(ng)%rmask(i,0)
#  endif
        END DO
# else
!
!  Southern edge, closed boundary condition.
!
        DO i=Istr,Iend
          zeta(i,0,knew)=zeta(i,1,knew)
#  ifdef MASKING
          zeta(i,0,knew)=zeta(i,0,knew)*                                &
     &                   GRID(ng)%rmask(i,0)
#  endif
        END DO
# endif
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the northern edge.
!-----------------------------------------------------------------------
!
      IF (NORTHERN_EDGE) THEN

# if defined NORTH_FSRADIATION
!
!  Northern edge, implicit upstream radiation condition.
!
        DO i=Istr,Iend+1
          grad(i,Mm(ng))=zeta(i  ,Mm(ng),know)-                         &
     &                   zeta(i-1,Mm(ng),know)
#  ifdef MASKING
          grad(i,Mm(ng))=grad(i,Mm(ng))*GRID(ng)%umask(i,Mm(ng))
#  endif
          grad(i,Mm(ng)+1)=zeta(i  ,Mm(ng)+1,know)-                     &
     &                     zeta(i-1,Mm(ng)+1,know)
#  ifdef MASKING
          grad(i,Mm(ng)+1)=grad(i,Mm(ng)+1)*GRID(ng)%umask(i,Mm(ng)+1)
#  endif
        END DO
        DO i=Istr,Iend
          dZdt=zeta(i,Mm(ng),know)-zeta(i,Mm(ng),knew)
          dZde=zeta(i,Mm(ng),knew)-zeta(i,Mm(ng)-1,knew)
#  ifdef NORTH_FSNUDGING
          tau=FSobc_out(ng,inorth)
          IF ((dZdt*dZde).lt.0.0_r8) tau=FSobc_in(ng,inorth)
          tau=tau*dt2d
#  endif
          IF ((dZdt*dZde).lt.0.0_r8) dZdt=0.0_r8
          IF ((dZdt*(grad(i,Mm(ng))+grad(i+1,Mm(ng)))).gt.0.0_r8) THEN
            dZdx=grad(i  ,Mm(ng))
          ELSE
            dZdx=grad(i+1,Mm(ng))
          END IF
          cff=MAX(dZdx*dZdx+dZde*dZde,eps)
#  ifdef RADIATION_2D
          Cx=MIN(cff,MAX(dZdt*dZdx,-cff))
#  else
          Cx=0.0_r8
#  endif
          Ce=dZdt*dZde
          zeta(i,Mm(ng)+1,knew)=(cff*zeta(i,Mm(ng)+1,know)+             &
     &                           Ce*zeta(i,Mm(ng),knew)-                &
     &                           MAX(Cx,0.0_r8)*grad(i  ,Mm(ng)+1)-     &
     &                           MIN(Cx,0.0_r8)*grad(i+1,Mm(ng)+1))/    &
     &                          (cff+Ce)
#  ifdef NORTH_FSNUDGING
          zeta(i,Mm(ng)+1,knew)=zeta(i,Mm(ng)+1,knew)+                  &
     &                          tau*(BOUNDARY(ng)%zeta_north(i)-        &
     &                               zeta(i,Mm(ng)+1,know))
#  endif
#  ifdef MASKING
          zeta(i,Mm(ng)+1,knew)=zeta(i,Mm(ng)+1,knew)*                  &
     &                          GRID(ng)%rmask(i,Mm(ng)+1)
#  endif
        END DO
# elif defined NORTH_FSCHAPMAN
!
!  Northern edge, Chapman boundary condition.
!
        DO i=Istr,Iend
          Ce=dt2d*GRID(ng)%pn(i,Mm(ng))*SQRT(g*GRID(ng)%h(i,Mm(ng)))
          zeta(i,Mm(ng)+1,knew)=(zeta(i,Mm(ng)+1,know)+                 &
     &                           Ce*zeta(i,Mm(ng),knew))/               &
     &                          (1.0_r8+Ce)
#  ifdef MASKING
          zeta(i,Mm(ng)+1,knew)=zeta(i,Mm(ng)+1,knew)*                  &
     &                          GRID(ng)%rmask(i,Mm(ng)+1)
#  endif
        END DO
# elif defined NORTH_FSCLAMPED
!
!  Northern edge, clamped boundary condition.
!
        DO i=Istr,Iend
          zeta(i,Mm(ng)+1,knew)=BOUNDARY(ng)%zeta_north(i)
#  ifdef MASKING
          zeta(i,Mm(ng)+1,knew)=zeta(i,Mm(ng)+1,knew)*                  &
     &                          GRID(ng)%rmask(i,Mm(ng)+1)
#  endif
        END DO
# elif defined NORTH_FSGRADIENT
!
!  Northern edge, gradient boundary condition.
!
        DO i=Istr,Iend
          zeta(i,Mm(ng)+1,knew)=zeta(i,Mm(ng),knew)
#  ifdef MASKING
          zeta(i,Mm(ng)+1,knew)=zeta(i,Mm(ng)+1,knew)*                  &
     &                          GRID(ng)%rmask(i,Mm(ng)+1)
#  endif
        END DO
# else
!
!  Northern edge, closed boundary condition.
!
        DO i=Istr,Iend
          zeta(i,Mm(ng)+1,knew)=zeta(i,Mm(ng),knew)
#  ifdef MASKING
          zeta(i,Mm(ng)+1,knew)=zeta(i,Mm(ng)+1,knew)*                  &
     &                          GRID(ng)%rmask(i,Mm(ng)+1)
#  endif
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
        zeta(0,0,knew)=0.5_r8*(zeta(1,0,knew)+                          &
     &                         zeta(0,1,knew))
      END IF
      IF (SOUTHERN_EDGE.and.EASTERN_EDGE) THEN
        zeta(Lm(ng)+1,0,knew)=0.5_r8*(zeta(Lm(ng)+1,1,knew)+            &
     &                                zeta(Lm(ng)  ,0,knew))
      END IF
      IF (NORTHERN_EDGE.and.WESTERN_EDGE) THEN
        zeta(0,Mm(ng)+1,knew)=0.5_r8*(zeta(0,Mm(ng)  ,knew)+            &
     &                                zeta(1,Mm(ng)+1,knew))
      END IF
      IF (NORTHERN_EDGE.and.EASTERN_EDGE) THEN
        zeta(Lm(ng)+1,Mm(ng)+1,knew)=0.5_r8*                            &
     &                               (zeta(Lm(ng)+1,Mm(ng)  ,knew)+     &
     &                                zeta(Lm(ng)  ,Mm(ng)+1,knew))
      END IF
#endif
      RETURN
      END SUBROUTINE zetabc_tile
      END MODULE zetabc_mod
