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
      PUBLIC  :: zetabc_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE zetabc (ng, kout, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_ocean
      USE mod_stepping
!
      integer, intent(in) :: ng, kout, tile

# include "tile.h"
!
      CALL zetabc_tile (ng, Istr, Iend, Jstr, Jend,                     &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  krhs(ng), kstp(ng), kout,                       &
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
          grad(Istr-1,j)=zeta(Istr-1,j  ,kstp)-                         &
     &                   zeta(Istr-1,j-1,kstp)
#  ifdef MASKING
          grad(Istr-1,j)=grad(Istr-1,j)*GRID(ng)%vmask(Istr-1,j)
#  endif
          grad(Istr,j)=zeta(Istr,j  ,kstp)-                             &
     &                 zeta(Istr,j-1,kstp)
#  ifdef MASKING
          grad(Istr,j)=grad(Istr,j)*GRID(ng)%vmask(Istr,j)
#  endif
        END DO
        DO j=Jstr,Jend
          dZdt=zeta(Istr,j,know)-zeta(Istr  ,j,knew)
          dZdx=zeta(Istr,j,knew)-zeta(Istr+1,j,knew)
#  ifdef WEST_FSNUDGING
          tau=FSobc_out(ng,iwest)
          IF ((dZdt*dZdx).lt.0.0_r8) tau=FSobc_in(ng,iwest)
          tau=tau*dt2d
#  endif
          IF ((dZdt*dZdx).lt.0.0_r8) dZdt=0.0_r8
          IF ((dZdt*(grad(Istr,j)+grad(Istr,j+1))).gt.0.0_r8) THEN
            dZde=grad(Istr,j  )
          ELSE
            dZde=grad(Istr,j+1)
          END IF
          cff=MAX(dZdx*dZdx+dZde*dZde,eps)
          Cx=dZdt*dZdx
#  ifdef RADIATION_2D
          Ce=MIN(cff,MAX(dZdt*dZde,-cff))
#  else
          Ce=0.0_r8
#  endif
          zeta(Istr-1,j,knew)=(cff*zeta(Istr-1,j,know)+                 &
     &                         Cx *zeta(Istr  ,j,knew)-                 &
     &                         MAX(Ce,0.0_r8)*grad(Istr-1,j  )-         &
     &                         MIN(Ce,0.0_r8)*grad(Istr-1,j+1))/        &
     &                        (cff+Cx)
#  ifdef WEST_FSNUDGING
          zeta(Istr-1,j,knew)=zeta(Istr-1,j,knew)+                      &
     &                        tau*(BOUNDARY(ng)%zeta_west(j)-           &
     &                             zeta(Istr-1,j,know))
#  endif
#  ifdef MASKING
          zeta(Istr-1,j,knew)=zeta(Istr-1,j,knew)*                      &
     &                        GRID(ng)%rmask(Istr-1,j)
#  endif
        END DO

# elif defined WEST_FSCHAPMAN
!
!  Western edge, Chapman boundary condition.
!
        DO j=Jstr,Jend
          Cx=dt2d*GRID(ng)%pm(Istr,j)*SQRT(g*GRID(ng)%h(Istr,j))
          zeta(Istr-1,j,knew)=(zeta(Istr-1,j,know)+                     &
     &                         Cx*zeta(Istr,j,knew))/                   &
     &                        (1.0_r8+Cx)
#  ifdef MASKING
          zeta(Istr-1,j,knew)=zeta(Istr-1,j,knew)*                      &
     &                        GRID(ng)%rmask(Istr-1,j)
#  endif
        END DO

# elif defined WEST_FSCLAMPED
!
!  Western edge, clamped boundary condition.
!
        DO j=Jstr,Jend
          zeta(Istr-1,j,knew)=BOUNDARY(ng)%zeta_west(j)
#  ifdef MASKING
          zeta(Istr-1,j,knew)=zeta(Istr-1,j,knew)*                      &
     &                        GRID(ng)%rmask(Istr-1,j)
#  endif
        END DO

# elif defined WEST_FSGRADIENT
!
!  Western edge, gradient boundary condition.
!
        DO j=Jstr,Jend
          zeta(Istr-1,j,knew)=zeta(Istr,j,knew)
#  ifdef MASKING
          zeta(Istr-1,j,knew)=zeta(Istr-1,j,knew)*                      &
     &                   GRID(ng)%rmask(Istr-1,j)
#  endif
        END DO

# else
!
!  Western edge, closed boundary condition.
!
        DO j=Jstr,Jend
          zeta(Istr-1,j,knew)=zeta(Istr,j,knew)
#  ifdef MASKING
          zeta(Istr-1,j,knew)=zeta(Istr-1,j,knew)*                      &
     &                        GRID(ng)%rmask(Istr-1,j)
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
          grad(Iend  ,j)=zeta(Iend  ,j  ,know)-                         &
     &                   zeta(Iend  ,j-1,know)
#  ifdef MASKING
          grad(Iend  ,j)=grad(Iend  ,j)*GRID(ng)%vmask(Iend  ,j)
#  endif
          grad(Iend+1,j)=zeta(Iend+1,j  ,know)-                         &
     &                   zeta(Iend+1,j-1,know)
#  ifdef MASKING
          grad(Iend+1,j)=grad(Iend+1,j)*GRID(ng)%vmask(Iend+1,j)
#  endif
        END DO
        DO j=Jstr,Jend
          dZdt=zeta(Iend,j,know)-zeta(Iend  ,j,knew)
          dZdx=zeta(Iend,j,knew)-zeta(Iend-1,j,knew)
#  ifdef EAST_FSNUDGING
          tau=FSobc_out(ng,ieast)
          IF ((dZdt*dZdx).lt.0.0_r8) tau=FSobc_in(ng,ieast)
          tau=tau*dt2d
#  endif
          IF ((dZdt*dZdx).lt.0.0_r8) dZdt=0.0_r8
          IF ((dZdt*(grad(Iend,j)+grad(Iend,j+1))).gt.0.0_r8) THEN
            dZde=grad(Iend,j  )
          ELSE
            dZde=grad(Iend,j+1)
          END IF
          cff=MAX(dZdx*dZdx+dZde*dZde,eps)
          Cx=dZdt*dZdx
#  ifdef RADIATION_2D
          Ce=MIN(cff,MAX(dZdt*dZde,-cff))
#  else
          Ce=0.0_r8
#  endif
          zeta(Iend+1,j,knew)=(cff*zeta(Iend+1,j,know)+                 &
     &                         Cx *zeta(Iend  ,j,knew)-                 &
     &                         MAX(Ce,0.0_r8)*grad(Iend+1,j  )-         &
     &                         MIN(Ce,0.0_r8)*grad(Iend+1,j+1))/        &
     &                        (cff+Cx)
#  ifdef EAST_FSNUDGING
          zeta(Iend+1,j,knew)=zeta(Iend+1,j,knew)+                      &
     &                        tau*(BOUNDARY(ng)%zeta_east(j)-           &
     &                             zeta(Iend+1,j,know))
#  endif
#  ifdef MASKING
          zeta(Iend+1,j,knew)=zeta(Iend+1,j,knew)*                      &
     &                        GRID(ng)%rmask(Iend+1,j)
#  endif
        END DO

# elif defined EAST_FSCHAPMAN
!
!  Eastern edge, Chapman boundary condition.
!
        DO j=Jstr,Jend
          Cx=dt2d*GRID(ng)%pm(Iend,j)*SQRT(g*GRID(ng)%h(Iend,j))
          zeta(Iend+1,j,knew)=(zeta(Iend+1,j,know)+                     &
     &                         Cx*zeta(Iend,j,knew))/                   &
     &                        (1.0_r8+Cx)
#  ifdef MASKING
          zeta(Iend+1,j,knew)=zeta(Iend+1,j,knew)*                      &
     &                        GRID(ng)%rmask(Iend+1,j)
#  endif
        END DO

# elif defined EAST_FSCLAMPED
!
!  Eastern edge, clamped boundary condition.
!
        DO j=Jstr,Jend
          zeta(Iend+1,j,knew)=BOUNDARY(ng)%zeta_east(j)
#  ifdef MASKING
          zeta(Iend+1,j,knew)=zeta(Iend+1,j,knew)*                      &
     &                        GRID(ng)%rmask(Iend+1,j)
#  endif
        END DO

# elif defined EAST_FSGRADIENT
!
!  Eastern edge, gradient boundary condition.
!
        DO j=Jstr,Jend
          zeta(Iend+1,j,knew)=zeta(Iend,j,knew)
#  ifdef MASKING
          zeta(Iend+1,j,knew)=zeta(Iend+1,j,knew)*                      &
     &                        GRID(ng)%rmask(Iend+1,j)
#  endif
        END DO

# else
!
!  Eastern edge, closed boundary condition.
!
        DO j=Jstr,Jend
          zeta(Iend+1,j,knew)=zeta(Iend,j,knew)
#  ifdef MASKING
          zeta(Iend+1,j,knew)=zeta(Iend+1,j,knew)*                      &
     &                        GRID(ng)%rmask(Iend+1,j)
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
          grad(i,Jstr  )=zeta(i  ,Jstr,know)-                           &
     &                   zeta(i-1,Jstr,know)
#  ifdef MASKING
          grad(i,Jstr  )=grad(i,Jstr  )*GRID(ng)%umask(i,Jstr  )
#  endif
          grad(i,Jstr-1)=zeta(i  ,Jstr-1,know)-                         &
     &                   zeta(i-1,Jstr-1,know)
#  ifdef MASKING
          grad(i,Jstr-1)=grad(i,Jstr-1)*GRID(ng)%umask(i,Jstr-1)
#  endif
        END DO
        DO i=Istr,Iend
          dZdt=zeta(i,Jstr,know)-zeta(i,Jstr  ,knew)
          dZde=zeta(i,Jstr,knew)-zeta(i,Jstr-1,knew)
#  ifdef SOUTH_FSNUDGING
          tau=FSobc_out(ng,isouth)
          IF ((dZdt*dZde).lt.0.0_r8) tau=FSobc_in(ng,isouth)
          tau=tau*dt2d
#  endif
          IF ((dZdt*dZde).lt.0.0_r8) dZdt=0.0_r8
          IF ((dZdt*(grad(i,Jstr)+grad(i+1,Jstr))).gt.0.0_r8) THEN
            dZdx=grad(i  ,Jstr)
          ELSE
            dZdx=grad(i+1,Jstr)
          END IF
          cff=MAX(dZdx*dZdx+dZde*dZde,eps)
#  ifdef RADIATION_2D
          Cx=MIN(cff,MAX(dZdt*dZdx,-cff))
#  else
          Cx=0.0_r8
#  endif
          Ce=dZdt*dZde
          zeta(i,Jstr-1,knew)=(cff*zeta(i,Jstr-1,know)+                 &
     &                         Ce *zeta(i,Jstr  ,knew)-                 &
     &                         MAX(Cx,0.0_r8)*grad(i  ,Jstr)-           &
     &                         MIN(Cx,0.0_r8)*grad(i+1,Jstr))/          &
     &                        (cff+Ce)
#  ifdef SOUTH_FSNUDGING
          zeta(i,Jstr-1,knew)=zeta(i,Jstr-1,knew)+                      &
     &                        tau*(BOUNDARY(ng)%zeta_south(i)-          &
     &                             zeta(i,Jstr-1,know))
#  endif
#  ifdef MASKING
          zeta(i,Jstr-1,knew)=zeta(i,Jstr-1,knew)*                      &
     &                        GRID(ng)%rmask(i,Jstr-1)
#  endif
        END DO

# elif defined SOUTH_FSCHAPMAN
!
!  Southern edge, Chapman boundary condition.
!
        DO i=Istr,Iend
          Ce=dt2d*GRID(ng)%pn(i,1)*SQRT(g*GRID(ng)%h(i,1))
          zeta(i,Jstr-1,knew)=(zeta(i,Jstr-1,know)+                     &
     &                         Ce*zeta(i,Jstr,knew))/                   &
     &                        (1.0_r8+Ce)
#  ifdef MASKING
          zeta(i,Jstr-1,knew)=zeta(i,Jstr-1,knew)*                      &
     &                        GRID(ng)%rmask(i,Jstr-1)
#  endif
        END DO

# elif defined SOUTH_FSCLAMPED
!
!  Southern edge, clamped boundary condition.
!
        DO i=Istr,Iend
          zeta(i,Jstr-1,knew)=BOUNDARY(ng)%zeta_south(i)
#  ifdef MASKING
          zeta(i,Jstr-1,knew)=zeta(i,Jstr-1,knew)*                      &
     &                        GRID(ng)%rmask(i,Jstr-1)
#  endif
        END DO

# elif defined SOUTH_FSGRADIENT
!
!  Southern edge, gradient boundary condition.
!
        DO i=Istr,Iend
          zeta(i,Jstr-1,knew)=zeta(i,Jstr,knew)
#  ifdef MASKING
          zeta(i,Jstr-1,knew)=zeta(i,Jstr-1,knew)*                      &
     &                        GRID(ng)%rmask(i,Jstr-1)
#  endif
        END DO

# else
!
!  Southern edge, closed boundary condition.
!
        DO i=Istr,Iend
          zeta(i,Jstr-1,knew)=zeta(i,Jstr,knew)
#  ifdef MASKING
          zeta(i,Jstr-1,knew)=zeta(i,Jstr-1,knew)*                      &
     &                        GRID(ng)%rmask(i,Jstr-1)
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
          grad(i,Jend  )=zeta(i  ,Jend  ,know)-                         &
     &                   zeta(i-1,Jend  ,know)
#  ifdef MASKING
          grad(i,Jend  )=grad(i,Jend  )*GRID(ng)%umask(i,Jend  )
#  endif
          grad(i,Jend+1)=zeta(i  ,Jend+1,know)-                         &
     &                   zeta(i-1,Jend+1,know)
#  ifdef MASKING
          grad(i,Jend+1)=grad(i,Jend+1)*GRID(ng)%umask(i,Jend+1)
#  endif
        END DO
        DO i=Istr,Iend
          dZdt=zeta(i,Jend,know)-zeta(i,Jend  ,knew)
          dZde=zeta(i,Jend,knew)-zeta(i,Jend-1,knew)
#  ifdef NORTH_FSNUDGING
          tau=FSobc_out(ng,inorth)
          IF ((dZdt*dZde).lt.0.0_r8) tau=FSobc_in(ng,inorth)
          tau=tau*dt2d
#  endif
          IF ((dZdt*dZde).lt.0.0_r8) dZdt=0.0_r8
          IF ((dZdt*(grad(i,Jend)+grad(i+1,Jend))).gt.0.0_r8) THEN
            dZdx=grad(i  ,Jend)
          ELSE
            dZdx=grad(i+1,Jend)
          END IF
          cff=MAX(dZdx*dZdx+dZde*dZde,eps)
#  ifdef RADIATION_2D
          Cx=MIN(cff,MAX(dZdt*dZdx,-cff))
#  else
          Cx=0.0_r8
#  endif
          Ce=dZdt*dZde
          zeta(i,Jend+1,knew)=(cff*zeta(i,Jend+1,know)+                 &
     &                         Ce *zeta(i,Jend  ,knew)-                 &
     &                         MAX(Cx,0.0_r8)*grad(i  ,Jend+1)-         &
     &                         MIN(Cx,0.0_r8)*grad(i+1,Jend+1))/        &
     &                        (cff+Ce)
#  ifdef NORTH_FSNUDGING
          zeta(i,Jend+1,knew)=zeta(i,Jend+1,knew)+                      &
     &                        tau*(BOUNDARY(ng)%zeta_north(i)-          &
     &                             zeta(i,Jend+1,know))
#  endif
#  ifdef MASKING
          zeta(i,Jend+1,knew)=zeta(i,Jend+1,knew)*                      &
     &                        GRID(ng)%rmask(i,Jend+1)
#  endif
        END DO

# elif defined NORTH_FSCHAPMAN
!
!  Northern edge, Chapman boundary condition.
!
        DO i=Istr,Iend
          Ce=dt2d*GRID(ng)%pn(i,Jend)*SQRT(g*GRID(ng)%h(i,Jend))
          zeta(i,Jend+1,knew)=(zeta(i,Jend+1,know)+                     &
     &                         Ce*zeta(i,Jend,knew))/                   &
     &                        (1.0_r8+Ce)
#  ifdef MASKING
          zeta(i,Jend+1,knew)=zeta(i,Jend+1,knew)*                      &
     &                        GRID(ng)%rmask(i,Jend+1)
#  endif
        END DO

# elif defined NORTH_FSCLAMPED
!
!  Northern edge, clamped boundary condition.
!
        DO i=Istr,Iend
          zeta(i,Jend+1,knew)=BOUNDARY(ng)%zeta_north(i)
#  ifdef MASKING
          zeta(i,Jend+1,knew)=zeta(i,Jend+1,knew)*                      &
     &                        GRID(ng)%rmask(i,Jend+1)
#  endif
        END DO

# elif defined NORTH_FSGRADIENT
!
!  Northern edge, gradient boundary condition.
!
        DO i=Istr,Iend
          zeta(i,Jend+1,knew)=zeta(i,Jend,knew)
#  ifdef MASKING
          zeta(i,Jend+1,knew)=zeta(i,Jend+1,knew)*                      &
     &                        GRID(ng)%rmask(i,Jend+1)
#  endif
        END DO

# else
!
!  Northern edge, closed boundary condition.
!
        DO i=Istr,Iend
          zeta(i,Jend+1,knew)=zeta(i,Jend,knew)
#  ifdef MASKING
          zeta(i,Jend+1,knew)=zeta(i,Jend+1,knew)*                      &
     &                        GRID(ng)%rmask(i,Jend+1)
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
        zeta(Istr-1,Jstr-1,knew)=0.5_r8*(zeta(Istr  ,Jstr-1,knew)+      &
     &                                   zeta(Istr-1,Jstr  ,knew))
      END IF
      IF (SOUTHERN_EDGE.and.EASTERN_EDGE) THEN
        zeta(Iend+1,Jstr-1,knew)=0.5_r8*(zeta(Iend+1,Jstr  ,knew)+      &
     &                                   zeta(Iend  ,Jstr-1,knew))
      END IF
      IF (NORTHERN_EDGE.and.WESTERN_EDGE) THEN
        zeta(Istr-1,Jend+1,knew)=0.5_r8*(zeta(Istr-1,Jend  ,knew)+      &
     &                                   zeta(Istr,Jend+1,knew))
      END IF
      IF (NORTHERN_EDGE.and.EASTERN_EDGE) THEN
        zeta(Iend+1,Jend+1,knew)=0.5_r8*(zeta(Iend+1,Jend  ,knew)+      &
     &                                   zeta(Iend  ,Jend+1,knew))
      END IF
#endif
      RETURN

      END SUBROUTINE zetabc_tile
      END MODULE zetabc_mod
