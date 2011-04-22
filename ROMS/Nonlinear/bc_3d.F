#include "cppdefs.h"
      MODULE bc_3d_mod
#ifdef SOLVE3D
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2011 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  These routines apply close, gradient or periodic boundary           !
!  conditions to generic 3D fields.                                    !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng                Nested grid number.                            !
!     tile              Domain partition.                              !
!     LBi               I-dimension Lower bound.                       !
!     UBi               I-dimension Upper bound.                       !
!     LBj               J-dimension Lower bound.                       !
!     UBj               J-dimension Upper bound.                       !
!     LBk               K-dimension Lower bound.                       !
!     UBk               K-dimension Upper bound.                       !
!     A                 3D field.                                      !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A                 Processed 3D field.                            !
!                                                                      !
!  Routines:                                                           !
!                                                                      !
!     bc_r3d_tile       Boundary conditions for field at RHO-points    !
!     bc_u3d_tile       Boundary conditions for field at U-points      !
!     bc_v3d_tile       Boundary conditions for field at V-points      !
!     bc_w3d_tile       Boundary conditions for field at W-points      !
!                                                                      !
!=======================================================================
!
      implicit none

      CONTAINS
!
!***********************************************************************
      SUBROUTINE bc_r3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, LBk, UBk,             &
     &                        A)
!***********************************************************************
!
      USE mod_param

# if defined EW_PERIODIC || defined NS_PERIODIC
!
      USE exchange_3d_mod, ONLY : exchange_r3d_tile
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: A(LBi:,LBj:,LBk:)
# else
      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj,LBk:UBk)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k

# include "set_bounds.h"

# ifndef EW_PERIODIC
!
!-----------------------------------------------------------------------
!  East-West gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
        DO k=LBk,UBk
          DO j=Jstr,Jend
            A(Iend+1,j,k)=A(Iend,j,k)
          END DO
        END DO
      END IF
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
        DO k=LBk,UBk
          DO j=Jstr,Jend
            A(Istr-1,j,k)=A(Istr,j,k)
          END DO
        END DO
      END IF
# endif

# ifndef NS_PERIODIC
!
!-----------------------------------------------------------------------
!  North-South gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
        DO k=LBk,UBk
          DO i=Istr,Iend
            A(i,Jend+1,k)=A(i,Jend,k)
          END DO
        END DO
      END IF
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
        DO k=LBk,UBk
          DO i=Istr,Iend
            A(i,Jstr-1,k)=A(i,Jstr,k)
          END DO
        END DO
      END IF
# endif

# if !defined EW_PERIODIC && !defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
        DO k=LBk,UBk
          A(Istr-1,Jstr-1,k)=0.5_r8*(A(Istr  ,Jstr-1,k)+                &
     &                               A(Istr-1,Jstr  ,k))
        END DO
      END IF
      IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
        DO k=LBk,UBk
          A(Iend+1,Jstr-1,k)=0.5_r8*(A(Iend  ,Jstr-1,k)+                &
     &                               A(Iend+1,Jstr  ,k))
        END DO
      END IF
      IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
        DO k=LBk,UBk
          A(Istr-1,Jend+1,k)=0.5_r8*(A(Istr-1,Jend  ,k)+                &
     &                               A(Istr  ,Jend+1,k))
        END DO
      END IF
      IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
        DO k=LBk,UBk
          A(Iend+1,Jend+1,k)=0.5_r8*(A(Iend+1,Jend  ,k)+                &
     &                               A(Iend  ,Jend+1,k))
        END DO
      END IF
# endif

# if defined EW_PERIODIC || defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Apply periodic boundary conditions.
!-----------------------------------------------------------------------
!
      CALL exchange_r3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, LBk, UBk,             &
     &                        A)
# endif

      RETURN
      END SUBROUTINE bc_r3d_tile

!
!***********************************************************************
      SUBROUTINE bc_u3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, LBk, UBk,             &
     &                        A)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_scalars

# if defined EW_PERIODIC || defined NS_PERIODIC
!
      USE exchange_3d_mod, ONLY : exchange_u3d_tile
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: A(LBi:,LBj:,LBk:)
# else
      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj,LBk:UBk)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k

# include "set_bounds.h"

# ifndef EW_PERIODIC
!
!-----------------------------------------------------------------------
!  East-West boundary conditions: Closed or gradient.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
        DO k=LBk,UBk
          DO j=Jstr,Jend
#  ifdef EASTERN_WALL
            A(Iend+1,j,k)=0.0_r8
#  else
            A(Iend+1,j,k)=A(Iend,j,k)
#  endif
          END DO
        END DO
      END IF
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
        DO k=LBk,UBk
          DO j=Jstr,Jend
#  ifdef WESTERN_WALL
            A(Istr,j,k)=0.0_r8
#  else
            A(Istr,j,k)=A(Istr+1,j,k)
#  endif
          END DO
        END DO
      END IF
# endif

# ifndef NS_PERIODIC
!
!-----------------------------------------------------------------------
!  North-South boundary conditions: Closed (free-slip/no-slip) or
!  gradient.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
#  ifdef NORTHERN_WALL
#   ifdef EW_PERIODIC
#    define I_RANGE IstrU,Iend
#   else
#    define I_RANGE Istr,IendR
#   endif
        DO k=LBk,UBk
          DO i=I_RANGE
            A(i,Jend+1,k)=gamma2(ng)*A(i,Jend,k)
#   ifdef MASKING
            A(i,Jend+1,k)=A(i,Jend+1,k)*GRID(ng)%umask(i,Jend+1)
#   endif
          END DO
        END DO
#   undef I_RANGE
#  else
        DO k=LBk,UBk
          DO i=IstrU,Iend
            A(i,Jend+1,k)=A(i,Jend,k)
          END DO
        END DO
#  endif
      END IF

      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
#  ifdef SOUTHERN_WALL
#   ifdef EW_PERIODIC
#    define I_RANGE IstrU,Iend
#   else
#    define I_RANGE Istr,IendR
#   endif
        DO k=LBk,UBk
          DO i=I_RANGE
            A(i,Jstr-1,k)=gamma2(ng)*A(i,Jstr,k)
#   ifdef MASKING
            A(i,Jstr-1,k)=A(i,Jstr-1,k)*GRID(ng)%umask(i,Jstr-1)
#   endif
          END DO
        END DO
#   undef I_RANGE
#  else
        DO k=LBk,UBk
          DO i=IstrU,Iend
            A(i,Jstr-1,k)=A(i,Jstr,k)
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
      IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
        DO k=LBk,UBk
          A(Istr  ,Jstr-1,k)=0.5_r8*(A(Istr+1,Jstr-1,k)+                &
     &                               A(Istr  ,Jstr  ,k))
        END DO
      END IF
      IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
        DO k=LBk,UBk
          A(Iend+1,Jstr-1,k)=0.5_r8*(A(Iend  ,Jstr-1,k)+                &
     &                               A(Iend+1,Jstr  ,k))
        END DO
      END IF
      IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
        DO k=LBk,UBk
          A(Istr  ,Jend+1,k)=0.5_r8*(A(Istr  ,Jend  ,k)+                &
     &                               A(Istr+1,Jend+1,k))
        END DO
      END IF
      IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
        DO k=LBk,UBk
          A(Iend+1,Jend+1,k)=0.5_r8*(A(Iend+1,Jend  ,k)+                &
     &                               A(Iend  ,Jend+1,k))
        END DO
      END IF
# endif

# if defined EW_PERIODIC || defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Apply periodic boundary conditions.
!-----------------------------------------------------------------------
!
      CALL exchange_u3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, LBk, UBk,             &
     &                        A)
# endif

      RETURN
      END SUBROUTINE bc_u3d_tile

!
!***********************************************************************
      SUBROUTINE bc_v3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, LBk, UBk,             &
     &                        A)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_scalars

# if defined EW_PERIODIC || defined NS_PERIODIC
!
      USE exchange_3d_mod, ONLY : exchange_v3d_tile
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: A(LBi:,LBj:,:)
# else
      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj,LBk:UBk)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k

# include "set_bounds.h"

# ifndef EW_PERIODIC
!
!-----------------------------------------------------------------------
!  East-West boundary conditions: Closed (free-slip/no-slip) or
!  gradient.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
#  ifdef EASTERN_WALL
#   ifdef NS_PERIODIC
#    define J_RANGE JstrV,Jend
#   else
#    define J_RANGE Jstr,JendR
#   endif
        DO k=LBk,UBk
          DO j=J_RANGE
            A(Iend+1,j,k)=gamma2(ng)*A(Iend,j,k)
#   ifdef MASKING
            A(Iend+1,j,k)=A(Iend+1,j,k)*GRID(ng)%vmask(Iend+1,j)
#   endif
          END DO
        END DO
#   undef J_RANGE
#  else
        DO k=LBk,UBk
          DO j=JstrV,Jend
            A(Iend+1,j,k)=A(Iend,j,k)
          END DO
        END DO
#  endif
      END IF

      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
#  ifdef WESTERN_WALL
#   ifdef NS_PERIODIC
#    define J_RANGE JstrV,Jend
#   else
#    define J_RANGE Jstr,JendR
#   endif
        DO k=LBk,UBk
          DO j=J_RANGE
            A(Istr-1,j,k)=gamma2(ng)*A(Istr,j,k)
#   ifdef MASKING
            A(Istr-1,j,k)=A(Istr-1,j,k)*GRID(ng)%vmask(Istr-1,j)
#   endif
          END DO
        END DO
#   undef J_RANGE
#  else
        DO k=LBk,UBk
          DO j=JstrV,Jend
            A(Istr-1,j,k)=A(Istr,j,k)
          END DO
        END DO
#  endif
      END IF
# endif

# ifndef NS_PERIODIC
!
!-----------------------------------------------------------------------
!  North-South boundary conditions: Closed or gradient.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
        DO k=LBk,UBk
          DO i=Istr,Iend
#  ifdef NORTHERN_WALL
            A(i,Jend+1,k)=0.0_r8
#  else
            A(i,Jend+1,k)=A(i,Jend,k)
#  endif
          END DO
        END DO
      END IF
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
        DO k=LBk,UBk
          DO i=Istr,Iend
#  ifdef SOUTHERN_WALL
            A(i,Jstr,k)=0.0_r8
#  else
            A(i,Jstr,k)=A(i,Jstr+1,k)
#  endif
          END DO
        END DO
      END IF
# endif

# if !defined EW_PERIODIC && !defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
        DO k=LBk,UBk
          A(Istr-1,Jstr  ,k)=0.5_r8*(A(Istr  ,Jstr  ,k)+                &
     &                               A(Istr-1,Jstr+1,k))
        END DO
      END IF
      IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
        DO k=LBk,UBk
          A(Iend+1,Jstr  ,k)=0.5_r8*(A(Iend  ,Jstr  ,k)+                &
     &                               A(Iend+1,Jstr+1,k))
        END DO
      END IF
      IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
        DO k=LBk,UBk
          A(Istr-1,Jend+1,k)=0.5_r8*(A(Istr-1,Jend  ,k)+                &
     &                               A(Istr  ,Jend+1,k))
        END DO
      END IF
      IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
        DO k=LBk,UBk
          A(Iend+1,Jend+1,k)=0.5_r8*(A(Iend+1,Jend  ,k)+                &
     &                               A(Iend  ,Jend+1,k))
        END DO
      END IF
# endif

# if defined EW_PERIODIC || defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Apply periodic boundary conditions.
!-----------------------------------------------------------------------
!
      CALL exchange_v3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, LBk, UBk,             &
     &                        A)
# endif

      RETURN
      END SUBROUTINE bc_v3d_tile

!
!***********************************************************************
      SUBROUTINE bc_w3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, LBk, UBk,             &
     &                        A)
!***********************************************************************
!
      USE mod_param

# if defined EW_PERIODIC || defined NS_PERIODIC
!
      USE exchange_3d_mod, ONLY : exchange_w3d_tile
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: A(LBi:,LBj:,LBk:)
# else
      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj,LBk:UBk)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k

# include "set_bounds.h"

# ifndef EW_PERIODIC
!
!-----------------------------------------------------------------------
!  East-West gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
        DO k=LBk,UBk
          DO j=Jstr,Jend
            A(Iend+1,j,k)=A(Iend,j,k)
          END DO
        END DO
      END IF
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
        DO k=LBk,UBk
          DO j=Jstr,Jend
            A(Istr-1,j,k)=A(Istr,j,k)
          END DO
        END DO
      END IF
# endif

# ifndef NS_PERIODIC
!
!-----------------------------------------------------------------------
!  North-South gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
        DO k=LBk,UBk
          DO i=Istr,Iend
            A(i,Jend+1,k)=A(i,Jend,k)
          END DO
        END DO
      END IF
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
        DO k=LBk,UBk
          DO i=Istr,Iend
            A(i,Jstr-1,k)=A(i,Jstr,k)
          END DO
        END DO
      END IF
# endif

# if !defined EW_PERIODIC && !defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
        DO k=LBk,UBk
          A(Istr-1,Jstr-1,k)=0.5_r8*(A(Istr  ,Jstr-1,k)+                &
     &                               A(Istr-1,Jstr  ,k))
        END DO
      END IF
      IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
        DO k=LBk,UBk
          A(Iend+1,Jstr-1,k)=0.5_r8*(A(Iend  ,Jstr-1,k)+                &
     &                               A(Iend+1,Jstr  ,k))
        END DO
      END IF
      IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
        DO k=LBk,UBk
          A(Istr-1,Jend+1,k)=0.5_r8*(A(Istr-1,Jend  ,k)+                &
     &                               A(Istr  ,Jend+1,k))
        END DO
      END IF
      IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
        DO k=LBk,UBk
          A(Iend+1,Jend+1,k)=0.5_r8*(A(Iend+1,Jend  ,k)+                &
     &                               A(Iend  ,Jend+1,k))
        END DO
      END IF
# endif

# if defined EW_PERIODIC || defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Apply periodic boundary conditions.
!-----------------------------------------------------------------------
!
      CALL exchange_w3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, LBk, UBk,             &
     &                        A)
# endif

      RETURN
      END SUBROUTINE bc_w3d_tile

#endif
      END MODULE bc_3d_mod
