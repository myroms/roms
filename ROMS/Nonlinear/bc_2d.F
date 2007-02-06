#include "cppdefs.h"
      MODULE bc_2d_mod
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2007 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
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
!***********************************************************************
      SUBROUTINE bc_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        A)
!***********************************************************************
!
      USE mod_param

#if defined EW_PERIODIC || defined NS_PERIODIC
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
          A(Iend+1,j)=A(Iend,j)
        END DO
      END IF
      IF (WESTERN_EDGE) THEN
        DO j=Jstr,Jend
          A(Istr-1,j)=A(Istr,j)
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
          A(i,Jend+1)=A(i,Jend)
        END DO
      END IF
      IF (SOUTHERN_EDGE) THEN
        DO i=Istr,Iend
          A(i,Jstr-1)=A(i,Jstr)
        END DO
      END IF
#endif

#if !defined EW_PERIODIC && !defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF ((SOUTHERN_EDGE).and.(WESTERN_EDGE)) THEN
        A(Istr-1,Jstr-1)=0.5_r8*(A(Istr  ,Jstr-1)+                      &
     &                           A(Istr-1,Jstr  ))
      END IF
      IF ((SOUTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        A(Iend+1,Jstr-1)=0.5_r8*(A(Iend  ,Jstr-1)+                      &
     &                           A(Iend+1,Jstr  ))
      END IF
      IF ((NORTHERN_EDGE).and.(WESTERN_EDGE)) THEN
        A(Istr-1,Jend+1)=0.5_r8*(A(Istr-1,Jend  )+                      &
     &                           A(Istr  ,Jend+1))
      END IF
      IF ((NORTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        A(Iend+1,Jend+1)=0.5_r8*(A(Iend+1,Jend  )+                      &
     &                           A(Iend  ,Jend+1))
      END IF
#endif

#if defined EW_PERIODIC || defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Exchange boundary data.
!-----------------------------------------------------------------------
!
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        A)
#endif

      RETURN
      END SUBROUTINE bc_r2d_tile

!
!***********************************************************************
      SUBROUTINE bc_u2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        A)
!***********************************************************************
!
      USE mod_param

#if defined EW_PERIODIC || defined NS_PERIODIC
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
          A(Iend+1,j)=A(Iend,j)
        END DO
      END IF
      IF (WESTERN_EDGE) THEN
        DO j=Jstr,Jend
          A(Istr,j)=A(Istr+1,j)
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
          A(i,Jend+1)=A(i,Jend)
        END DO
      END IF
      IF (SOUTHERN_EDGE) THEN
        DO i=IstrU,Iend
          A(i,Jstr-1)=A(i,Jstr)
        END DO
      END IF
#endif

#if !defined EW_PERIODIC && !defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF ((SOUTHERN_EDGE).and.(WESTERN_EDGE)) THEN
        A(Istr  ,Jstr-1)=0.5_r8*(A(Istr+1,Jstr-1)+                      &
     &                           A(Istr  ,Jstr  ))
      END IF
      IF ((SOUTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        A(Iend+1,Jstr-1)=0.5_r8*(A(Iend  ,Jstr-1)+                      &
     &                           A(Iend+1,Jstr  ))
      END IF
      IF ((NORTHERN_EDGE).and.(WESTERN_EDGE)) THEN
        A(Istr  ,Jend+1)=0.5_r8*(A(Istr  ,Jend  )+                      &
     &                           A(Istr+1,Jend+1)) 
      END IF
      IF ((NORTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        A(Iend+1,Jend+1)=0.5_r8*(A(Iend+1,Jend  )+                      &
     &                           A(Iend  ,Jend+1))
      END IF
#endif

#if defined EW_PERIODIC || defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Exchange boundary data.
!-----------------------------------------------------------------------
!
      CALL exchange_u2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        A)
#endif

      RETURN
      END SUBROUTINE bc_u2d_tile

!
!***********************************************************************
      SUBROUTINE bc_v2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        A)
!***********************************************************************
!
      USE mod_param

#if defined EW_PERIODIC || defined NS_PERIODIC
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
          A(Iend+1,j)=A(Iend,j)
        END DO
      END IF
      IF (WESTERN_EDGE) THEN
        DO j=JstrV,Jend
          A(Istr-1,j)=A(Istr,j)
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
          A(i,Jend+1)=A(i,Jend)
        END DO
      END IF
      IF (SOUTHERN_EDGE) THEN
        DO i=Istr,Iend
          A(i,Jstr)=A(i,Jstr+1)
        END DO
      END IF
#endif

#if !defined EW_PERIODIC && !defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF ((SOUTHERN_EDGE).and.(WESTERN_EDGE)) THEN
        A(Istr-1,Jstr  )=0.5_r8*(A(Istr  ,Jstr  )+                      &
     &                           A(Istr-1,Jstr+1))
      END IF
      IF ((SOUTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        A(Iend+1,Jstr  )=0.5_r8*(A(Iend  ,Jstr  )+                      &
     &                           A(Iend+1,Jstr+1))
      END IF
      IF ((NORTHERN_EDGE).and.(WESTERN_EDGE)) THEN
        A(Istr-1,Jend+1)=0.5_r8*(A(Istr-1,Jend  )+                      &
     &                           A(Istr  ,Jend+1))
      END IF
      IF ((NORTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        A(Iend+1,Jend+1)=0.5_r8*(A(Iend+1,Jend  )+                      &
     &                           A(Iend  ,Jend+1))
      END IF
#endif

#if defined EW_PERIODIC || defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Exchange boundary data.
!-----------------------------------------------------------------------
!
      CALL exchange_v2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        A)
#endif

      RETURN
      END SUBROUTINE bc_v2d_tile

      END MODULE bc_2d_mod
