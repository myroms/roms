#include "cppdefs.h"
      MODULE bc_2d_mod
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This package applies gradient or periodic boundary conditions for   !
!  generic 2D fields.                                                  !
!                                                                      !
!  Routines:                                                           !
!                                                                      !
!    bc_r2d_tile        Boundary conditions for field at RHO-points    !
!    bc_u2d_tile        Boundary conditions for field at U-points      !
!    bc_v2d_tile        Boundary conditions for field at V-points      !
!                                                                      !
!=======================================================================
!
      implicit none

      CONTAINS
!
      SUBROUTINE bc_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, A)
!
!=======================================================================
!                                                                      !
!  This routine sets gradient or periodic boundary conditions for      !
!  2D variables at RHO-points.                                         !
!                                                                      !
!=======================================================================
!
      USE mod_param

#if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
!
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj

#ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: A(LBi:,LBj:)
#else
      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj)
#endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j

#include "set_bounds.h"

#ifndef EW_PERIODIC
!
!-----------------------------------------------------------------------
!  East-West gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (EASTERN_EDGE) THEN
        DO j=Jstr,Jend
          A(Lm(ng)+1,j)=A(Lm(ng),j)
        END DO
      END IF
      IF (WESTERN_EDGE) THEN
        DO j=Jstr,Jend
          A(0,j)=A(1,j)
        END DO
      END IF
#endif
#ifndef NS_PERIODIC
!
!-----------------------------------------------------------------------
!  North-South gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (NORTHERN_EDGE) THEN
        DO i=Istr,Iend
          A(i,Mm(ng)+1)=A(i,Mm(ng))
        END DO
      END IF
      IF (SOUTHERN_EDGE) THEN
        DO i=Istr,Iend
          A(i,0)=A(i,1)
        END DO
      END IF
#endif
#if !defined EW_PERIODIC && !defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (SOUTHERN_EDGE.and.WESTERN_EDGE) THEN
        A(0,0)=0.5_r8*(A(1,0)+A(0,1))
      END IF
      IF (SOUTHERN_EDGE.and.EASTERN_EDGE) THEN
        A(Lm(ng)+1,0)=0.5_r8*(A(Lm(ng)+1,1)+A(Lm(ng),0))
      END IF
      IF (NORTHERN_EDGE.and.WESTERN_EDGE) THEN
        A(0,Mm(ng)+1)=0.5_r8*(A(0,Mm(ng))+A(1,Mm(ng)+1))
      END IF
      IF (NORTHERN_EDGE.and.EASTERN_EDGE) THEN
        A(Lm(ng)+1,Mm(ng)+1)=0.5_r8*(A(Lm(ng)+1,Mm(ng)  )+              &
     &                               A(Lm(ng)  ,Mm(ng)+1))
      END IF
#endif
#if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Exchange boundary data.
!-----------------------------------------------------------------------
!
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, A)
#endif
      RETURN
      END SUBROUTINE bc_r2d_tile

      SUBROUTINE bc_u2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, A)
!
!=======================================================================
!                                                                      !
!  This routine sets gradient or periodic boundary conditions for      !
!  2D variables at U-points.                                           !
!                                                                      !
!=======================================================================
!
      USE mod_param

#if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
!
      USE exchange_2d_mod, ONLY : exchange_u2d_tile
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj

#ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: A(LBi:,LBj:)
#else
      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj)
#endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j

#include "set_bounds.h"

#ifndef EW_PERIODIC
!
!-----------------------------------------------------------------------
!  East-West gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (EASTERN_EDGE) THEN
        DO j=Jstr,Jend
          A(Lm(ng)+1,j)=A(Lm(ng),j)
        END DO
      END IF
      IF (WESTERN_EDGE) THEN
        DO j=Jstr,Jend
          A(1,j)=A(2,j)
        END DO
      END IF
#endif
#ifndef NS_PERIODIC
!
!-----------------------------------------------------------------------
!  North-South gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (NORTHERN_EDGE) THEN
        DO i=IstrU,Iend
          A(i,Mm(ng)+1)=A(i,Mm(ng))
        END DO
      END IF
      IF (SOUTHERN_EDGE) THEN
        DO i=IstrU,Iend
          A(i,0)=A(i,1)
        END DO
      END IF
#endif
#if !defined EW_PERIODIC && !defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (SOUTHERN_EDGE.and.WESTERN_EDGE) THEN
        A(1,0)=0.5_r8*(A(2,0)+A(1,1))
      END IF
      IF (SOUTHERN_EDGE.and.EASTERN_EDGE) THEN
        A(Lm(ng)+1,0)=0.5_r8*(A(Lm(ng),0)+A(Lm(ng)+1,1))
      END IF
      IF (NORTHERN_EDGE.and.WESTERN_EDGE) THEN
        A(1,Mm(ng)+1)=0.5_r8*(A(2,Mm(ng)+1)+A(1,Mm(ng)))
      END IF
      IF (NORTHERN_EDGE.and.EASTERN_EDGE) THEN
        A(Lm(ng)+1,Mm(ng)+1)=0.5_r8*(A(Lm(ng)  ,Mm(ng)+1)+              &
     &                               A(Lm(ng)+1,Mm(ng)  ))
      END IF
#endif
#if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Exchange boundary data.
!-----------------------------------------------------------------------
!
      CALL exchange_u2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, A)
#endif
      RETURN
      END SUBROUTINE bc_u2d_tile

      SUBROUTINE bc_v2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, A)
!
!=======================================================================
!                                                                      !
!  This routine sets gradient or periodic boundary conditions for      !
!  2D variables at V-points.                                           !
!                                                                      !
!=======================================================================
!
      USE mod_param

#if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
!
      USE exchange_2d_mod, ONLY : exchange_v2d_tile
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj

#ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: A(LBi:,LBj:)
#else
      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj)
#endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j

#include "set_bounds.h"

#ifndef EW_PERIODIC
!
!-----------------------------------------------------------------------
!  East-West gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (EASTERN_EDGE) THEN
        DO j=JstrV,Jend
          A(Lm(ng)+1,j)=A(Lm(ng),j)
        END DO
      END IF
      IF (WESTERN_EDGE) THEN
        DO j=JstrV,Jend
          A(0,j)=A(1,j)
        END DO
      END IF
#endif
#ifndef NS_PERIODIC
!
!-----------------------------------------------------------------------
!  North-South periodic boundary conditions.
!-----------------------------------------------------------------------
!
      IF (NORTHERN_EDGE) THEN
        DO i=Istr,Iend
          A(i,Mm(ng)+1)=A(i,Mm(ng))
        END DO
      END IF
      IF (SOUTHERN_EDGE) THEN
        DO i=Istr,Iend
          A(i,1)=A(i,2)
        END DO
      END IF
#endif
#if !defined EW_PERIODIC && !defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (SOUTHERN_EDGE.and.WESTERN_EDGE) THEN
        A(0,1)=0.5_r8*(A(0,2)+A(1,1))
      END IF
      IF (SOUTHERN_EDGE.and.EASTERN_EDGE) THEN
        A(Lm(ng)+1,0)=0.5_r8*(A(Lm(ng)+1,1)+A(Lm(ng),0))
      END IF
      IF (NORTHERN_EDGE.and.WESTERN_EDGE) THEN
        A(0,Mm(ng)+1)=0.5_r8*(A(0,Mm(ng))+A(1,Mm(ng)+1))
      END IF
      IF (NORTHERN_EDGE.and.EASTERN_EDGE) THEN
        A(Lm(ng)+1,Mm(ng)+1)=0.5_r8*(A(Lm(ng)+1,Mm(ng)  )+              &
     &                               A(Lm(ng)  ,Mm(ng)+1))
      END IF
#endif
#if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Exchange boundary data.
!-----------------------------------------------------------------------
!
      CALL exchange_v2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, A)
#endif
      RETURN
      END SUBROUTINE bc_v2d_tile
      END MODULE bc_2d_mod
