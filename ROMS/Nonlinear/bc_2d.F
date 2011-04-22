#include "cppdefs.h"
      MODULE bc_2d_mod
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2011 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  These routines apply close, gradient or periodic boundary           !
!  conditions to generic 2D fields.                                    !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng                Nested grid number.                            !
!     tile              Domain partition.                              !
!     LBi               I-dimension Lower bound.                       !
!     UBi               I-dimension Upper bound.                       !
!     LBj               J-dimension Lower bound.                       !
!     UBj               J-dimension Upper bound.                       !
!     A                 2D field.                                      !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A                 Processed 2D field.                            !
!                                                                      !
!  Routines:                                                           !
!                                                                      !
!     bc_r2d_tile       Boundary conditions for field at RHO-points    !
!     bc_u2d_tile       Boundary conditions for field at U-points      !
!     bc_v2d_tile       Boundary conditions for field at V-points      !
!                                                                      !
!=======================================================================
!
      implicit none

      CONTAINS
!
!***********************************************************************
      SUBROUTINE bc_r2d_tile (ng, tile,                                 &
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
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj

#ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: A(LBi:,LBj:)
#else
      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj)
#endif
!
!  Local variable declarations.
!
      integer :: i, j

#include "set_bounds.h"

#ifndef EW_PERIODIC
!
!-----------------------------------------------------------------------
!  East-West gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
        DO j=Jstr,Jend
          A(Iend+1,j)=A(Iend,j)
        END DO
      END IF
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
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
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
        DO i=Istr,Iend
          A(i,Jend+1)=A(i,Jend)
        END DO
      END IF
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
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
      IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
        A(Istr-1,Jstr-1)=0.5_r8*(A(Istr  ,Jstr-1)+                      &
     &                           A(Istr-1,Jstr  ))
      END IF
      IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
        A(Iend+1,Jstr-1)=0.5_r8*(A(Iend  ,Jstr-1)+                      &
     &                           A(Iend+1,Jstr  ))
      END IF
      IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
        A(Istr-1,Jend+1)=0.5_r8*(A(Istr-1,Jend  )+                      &
     &                           A(Istr  ,Jend+1))
      END IF
      IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
        A(Iend+1,Jend+1)=0.5_r8*(A(Iend+1,Jend  )+                      &
     &                           A(Iend  ,Jend+1))
      END IF
#endif

#if defined EW_PERIODIC || defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Apply periodic boundary conditions.
!-----------------------------------------------------------------------
!
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        A)
#endif

      RETURN
      END SUBROUTINE bc_r2d_tile

!
!***********************************************************************
      SUBROUTINE bc_u2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        A)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_scalars

#if defined EW_PERIODIC || defined NS_PERIODIC
!
      USE exchange_2d_mod, ONLY : exchange_u2d_tile
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj

#ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: A(LBi:,LBj:)
#else
      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj)
#endif
!
!  Local variable declarations.
!
      integer :: i, j

#include "set_bounds.h"

#ifndef EW_PERIODIC
!
!-----------------------------------------------------------------------
!  East-West boundary conditions: Closed or gradient
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
        DO j=Jstr,Jend
# ifdef EASTERN_WALL
          A(Iend+1,j)=0.0_r8
# else
          A(Iend+1,j)=A(Iend,j)
# endif
        END DO
      END IF
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
        DO j=Jstr,Jend
# ifdef WESTERN_WALL
          A(Istr,j)=0.0_r8
# else
          A(Istr,j)=A(Istr+1,j)
# endif
        END DO
      END IF
#endif

#ifndef NS_PERIODIC
!
!-----------------------------------------------------------------------
!  North-South boundary conditions: Closed (free-slip/no-slip) or
!  gradient.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
# ifdef NORTHERN_WALL
#  ifdef EW_PERIODIC
#   define I_RANGE IstrU,Iend
#  else
#   define I_RANGE Istr,IendR
#  endif
        DO i=I_RANGE
          A(i,Jend+1)=gamma2(ng)*A(i,Jend)
#  ifdef MASKING
          A(i,Jend+1)=A(i,Jend+1)*GRID(ng)%umask(i,Jend+1)
#  endif
        END DO
#  undef I_RANGE
# else
        DO i=IstrU,Iend
          A(i,Jend+1)=A(i,Jend)
        END DO
# endif
      END IF

      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
# ifdef SOUTHERN_WALL
#  ifdef EW_PERIODIC
#   define I_RANGE IstrU,Iend
#  else
#   define I_RANGE Istr,IendR
#  endif
        DO i=I_RANGE
          A(i,Jstr-1)=gamma2(ng)*A(i,Jstr)
#  ifdef MASKING
          A(i,Jstr-1)=A(i,Jstr-1)*GRID(ng)%umask(i,Jstr-1)
#  endif
        END DO
#  undef I_RANGE
# else
        DO i=IstrU,Iend
          A(i,Jstr-1)=A(i,Jstr)
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
      IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
        A(Istr  ,Jstr-1)=0.5_r8*(A(Istr+1,Jstr-1)+                      &
     &                           A(Istr  ,Jstr  ))
      END IF
      IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
        A(Iend+1,Jstr-1)=0.5_r8*(A(Iend  ,Jstr-1)+                      &
     &                           A(Iend+1,Jstr  ))
      END IF
      IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
        A(Istr  ,Jend+1)=0.5_r8*(A(Istr  ,Jend  )+                      &
     &                           A(Istr+1,Jend+1))
      END IF
      IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
        A(Iend+1,Jend+1)=0.5_r8*(A(Iend+1,Jend  )+                      &
     &                           A(Iend  ,Jend+1))
      END IF
#endif

#if defined EW_PERIODIC || defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Apply periodic boundary conditions.
!-----------------------------------------------------------------------
!
      CALL exchange_u2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        A)
#endif

      RETURN
      END SUBROUTINE bc_u2d_tile

!
!***********************************************************************
      SUBROUTINE bc_v2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        A)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_scalars

#if defined EW_PERIODIC || defined NS_PERIODIC
!
      USE exchange_2d_mod, ONLY : exchange_v2d_tile
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj

#ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: A(LBi:,LBj:)
#else
      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj)
#endif
!
!  Local variable declarations.
!
      integer :: i, j

#include "set_bounds.h"

#ifndef EW_PERIODIC
!
!-----------------------------------------------------------------------
!  East-West boundary conditions: Closed (free-slip/no-slip) or
!  gradient.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
# ifdef EASTERN_WALL
#  ifdef NS_PERIODIC
#   define J_RANGE JstrV,Jend
#  else
#   define J_RANGE Jstr,JendR
#  endif
        DO j=J_RANGE
          A(Iend+1,j)=gamma2(ng)*A(Iend,j)
#  ifdef MASKING
          A(Iend+1,j)=A(Iend+1,j)*GRID(ng)%vmask(Iend+1,j)
#  endif
        END DO
#  undef J_RANGE
# else
        DO j=JstrV,Jend
          A(Iend+1,j)=A(Iend,j)
        END DO
# endif
      END IF

      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
# ifdef WESTERN_WALL
#  ifdef NS_PERIODIC
#   define J_RANGE JstrV,Jend
#  else
#   define J_RANGE Jstr,JendR
#  endif
        DO j=J_RANGE
          A(Istr-1,j)=gamma2(ng)*A(Istr,j)
#  ifdef MASKING
          A(Istr-1,j)=A(Istr-1,j)*GRID(ng)%vmask(Istr-1,j)
#  endif
        END DO
#  undef J_RANGE
# else
        DO j=JstrV,Jend
          A(Istr-1,j)=A(Istr,j)
        END DO
# endif
      END IF
#endif

#ifndef NS_PERIODIC
!
!-----------------------------------------------------------------------
!  North-South boundary conditions: Closed or Gradient.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
        DO i=Istr,Iend
# ifdef NORTHERN_WALL
          A(i,Jend+1)=0.0_r8
# else
          A(i,Jend+1)=A(i,Jend)
# endif
        END DO
      END IF
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
        DO i=Istr,Iend
# ifdef SOUTHERN_WALL
          A(i,Jstr)=0.0_r8
# else
          A(i,Jstr)=A(i,Jstr+1)
# endif
        END DO
      END IF
#endif

#if !defined EW_PERIODIC && !defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
        A(Istr-1,Jstr  )=0.5_r8*(A(Istr  ,Jstr  )+                      &
     &                           A(Istr-1,Jstr+1))
      END IF
      IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
        A(Iend+1,Jstr  )=0.5_r8*(A(Iend  ,Jstr  )+                      &
     &                           A(Iend+1,Jstr+1))
      END IF
      IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
        A(Istr-1,Jend+1)=0.5_r8*(A(Istr-1,Jend  )+                      &
     &                           A(Istr  ,Jend+1))
      END IF
      IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
        A(Iend+1,Jend+1)=0.5_r8*(A(Iend+1,Jend  )+                      &
     &                           A(Iend  ,Jend+1))
      END IF
#endif

#if defined EW_PERIODIC || defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Apply periodic boundary conditions.
!-----------------------------------------------------------------------
!
      CALL exchange_v2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        A)
#endif

      RETURN
      END SUBROUTINE bc_v2d_tile

      END MODULE bc_2d_mod
