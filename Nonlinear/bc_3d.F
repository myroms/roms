#include "cppdefs.h"
      MODULE bc_3d_mod
#ifdef SOLVE3D
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This package applies gradient or periodic boundary conditions for   !
!  generic 3D fields.                                                  !
!                                                                      !
!  Routines:                                                           !
!                                                                      !
!    bc_r3d_tile        Boundary conditions for field at RHO-points    !
!    bc_u3d_tile        Boundary conditions for field at U-points      !
!    bc_v3d_tile        Boundary conditions for field at V-points      !
!    bc_w3d_tile        Boundary conditions for field at W-points      !
!                                                                      !
!=======================================================================
!
      implicit none

      CONTAINS
!
      SUBROUTINE bc_r3d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, A)
!
!=======================================================================
!                                                                      !
!  This routine sets gradient or periodic boundary conditions for      !
!  3D variables at RHO-points.                                         !
!                                                                      !
!=======================================================================
!
      USE mod_param

# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
!
      USE exchange_3d_mod, ONLY : exchange_r3d_tile
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj

# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: A(LBi:,LBj:,:)
# else
      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj,N(ng))
# endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j, k

# include "set_bounds.h"

# ifndef EW_PERIODIC
!
!-----------------------------------------------------------------------
!  East-West gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (EASTERN_EDGE) THEN
        DO k=1,N(ng)
          DO j=Jstr,Jend
            A(Lm(ng)+1,j,k)=A(Lm(ng),j,k)
          END DO
        END DO
      END IF
      IF (WESTERN_EDGE) THEN
        DO k=1,N(ng)
          DO j=Jstr,Jend
            A(0,j,k)=A(1,j,k)
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
      IF (NORTHERN_EDGE) THEN
        DO k=1,N(ng)
          DO i=Istr,Iend
            A(i,Mm(ng)+1,k)=A(i,Mm(ng),k)
          END DO
        END DO
      END IF
      IF (SOUTHERN_EDGE) THEN
        DO k=1,N(ng)
          DO i=Istr,Iend
            A(i,0,k)=A(i,1,k)
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
      IF (SOUTHERN_EDGE.and.WESTERN_EDGE) THEN
        DO k=1,N(ng)
          A(0,0,k)=0.5_r8*(A(1,0,k)+A(0,1,k))
        END DO
      END IF
      IF (SOUTHERN_EDGE.and.EASTERN_EDGE) THEN
        DO k=1,N(ng)
          A(Lm(ng)+1,0,k)=0.5_r8*(A(Lm(ng)+1,1,k)+A(Lm(ng),0,k))
        END DO
      END IF
      IF (NORTHERN_EDGE.and.WESTERN_EDGE) THEN
        DO k=1,N(ng)
          A(0,Mm(ng)+1,k)=0.5_r8*(A(0,Mm(ng),k)+A(1,Mm(ng)+1,k))
        END DO
      END IF
      IF (NORTHERN_EDGE.and.EASTERN_EDGE) THEN
        DO k=1,N(ng)
          A(Lm(ng)+1,Mm(ng)+1,k)=0.5_r8*(A(Lm(ng)+1,Mm(ng)  ,k)+        &
     &                                   A(Lm(ng)  ,Mm(ng)+1,k))
        END DO
      END IF
# endif
# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Exchange boundary data.
!-----------------------------------------------------------------------
!
      CALL exchange_r3d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, A)
# endif
      RETURN
      END SUBROUTINE bc_r3d_tile

      SUBROUTINE bc_u3d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, A)
!
!=======================================================================
!                                                                      !
!  This routine sets gradient or periodic boundary conditions for      !
!  3D variables at U-points.                                           !
!                                                                      !
!=======================================================================
!
      USE mod_param

# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
!
      USE exchange_3d_mod, ONLY : exchange_u3d_tile
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj

# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: A(LBi:,LBj:,:)
# else
      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj,N(ng))
# endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j, k

# include "set_bounds.h"

# ifndef EW_PERIODIC
!
!-----------------------------------------------------------------------
!  East-West gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (EASTERN_EDGE) THEN
        DO k=1,N(ng)
          DO j=Jstr,Jend
            A(Lm(ng)+1,j,k)=A(Lm(ng),j,k)
          END DO
        END DO
      END IF
      IF (WESTERN_EDGE) THEN
        DO k=1,N(ng)
          DO j=Jstr,Jend
            A(1,j,k)=A(2,j,k)
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
      IF (NORTHERN_EDGE) THEN
        DO k=1,N(ng)
          DO i=IstrU,Iend
            A(i,Mm(ng)+1,k)=A(i,Mm(ng),k)
          END DO
        END DO
      END IF
      IF (SOUTHERN_EDGE) THEN
        DO k=1,N(ng)
          DO i=IstrU,Iend
            A(i,0,k)=A(i,1,k)
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
      IF (SOUTHERN_EDGE.and.WESTERN_EDGE) THEN
        DO k=1,N(ng)
          A(1,0,k)=0.5_r8*(A(2,0,k)+A(1,1,k))
        END DO
      END IF
      IF (SOUTHERN_EDGE.and.EASTERN_EDGE) THEN
        DO k=1,N(ng)
          A(Lm(ng)+1,0,k)=0.5_r8*(A(Lm(ng),0,k)+A(Lm(ng)+1,1,k))
        END DO
      END IF
      IF (NORTHERN_EDGE.and.WESTERN_EDGE) THEN
        DO k=1,N(ng)
          A(1,Mm(ng)+1,k)=0.5_r8*(A(2,Mm(ng)+1,k)+A(1,Mm(ng),k))
        END DO
      END IF
      IF (NORTHERN_EDGE.and.EASTERN_EDGE) THEN
        DO k=1,N(ng)
          A(Lm(ng)+1,Mm(ng)+1,k)=0.5_r8*(A(Lm(ng)  ,Mm(ng)+1,k)+        &
     &                                   A(Lm(ng)+1,Mm(ng)  ,k))
        END DO
      END IF
# endif
# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Exchange boundary data.
!-----------------------------------------------------------------------
!
      CALL exchange_u3d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, A)
# endif
      RETURN
      END SUBROUTINE bc_u3d_tile

      SUBROUTINE bc_v3d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, A)
!
!=======================================================================
!                                                                      !
!  This routine sets gradient or periodic boundary conditions for      !
!  3D variables at V-points.                                           !
!                                                                      !
!=======================================================================
!
      USE mod_param

# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
!
      USE exchange_3d_mod, ONLY : exchange_v3d_tile
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj

# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: A(LBi:,LBj:,:)
# else
      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj,N(ng))
# endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j, k

# include "set_bounds.h"

# ifndef EW_PERIODIC
!
!-----------------------------------------------------------------------
!  East-West gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (EASTERN_EDGE) THEN
        DO k=1,N(ng)
          DO j=JstrV,Jend
            A(Lm(ng)+1,j,k)=A(Lm(ng),j,k)
          END DO
        END DO
      END IF
      IF (WESTERN_EDGE) THEN
        DO k=1,N(ng)
          DO j=JstrV,Jend
            A(0,j,k)=A(1,j,k)
          END DO
        END DO
      END IF
# endif
# ifndef NS_PERIODIC
!
!-----------------------------------------------------------------------
!  North-South periodic boundary conditions.
!-----------------------------------------------------------------------
!
      IF (NORTHERN_EDGE) THEN
        DO k=1,N(ng)
          DO i=Istr,Iend
            A(i,Mm(ng)+1,k)=A(i,Mm(ng),k)
          END DO
        END DO
      END IF
      IF (SOUTHERN_EDGE) THEN
        DO k=1,N(ng)
          DO i=Istr,Iend
            A(i,1,k)=A(i,2,k)
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
      IF (SOUTHERN_EDGE.and.WESTERN_EDGE) THEN
        DO k=1,N(ng)
          A(0,1,k)=0.5_r8*(A(0,2,k)+A(1,1,k))
        END DO
      END IF
      IF (SOUTHERN_EDGE.and.EASTERN_EDGE) THEN
        DO k=1,N(ng)
          A(Lm(ng)+1,0,k)=0.5_r8*(A(Lm(ng)+1,1,k)+A(Lm(ng),0,k))
        END DO
      END IF
      IF (NORTHERN_EDGE.and.WESTERN_EDGE) THEN
        DO k=1,N(ng)
          A(0,Mm(ng)+1,k)=0.5_r8*(A(0,Mm(ng),k)+A(1,Mm(ng)+1,k))
        END DO
      END IF
      IF (NORTHERN_EDGE.and.EASTERN_EDGE) THEN
        DO k=1,N(ng)
          A(Lm(ng)+1,Mm(ng)+1,k)=0.5_r8*(A(Lm(ng)+1,Mm(ng)  ,k)+        &
     &                                   A(Lm(ng)  ,Mm(ng)+1,k))
        END DO
      END IF
# endif
# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Exchange boundary data.
!-----------------------------------------------------------------------
!
      CALL exchange_v3d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, A)
# endif
      RETURN
      END SUBROUTINE bc_v3d_tile

      SUBROUTINE bc_w3d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, A)
!
!=======================================================================
!                                                                      !
!  This routine sets gradient or periodic boundary conditions for      !
!  3D variables at W-points.                                           !
!                                                                      !
!=======================================================================
!
      USE mod_param

# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
!
      USE exchange_3d_mod, ONLY : exchange_w3d_tile
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj

# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: A(LBi:,LBj:,0:)
# else
      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj,0:N(ng))
# endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j, k

# include "set_bounds.h"

# ifndef EW_PERIODIC
!
!-----------------------------------------------------------------------
!  East-West gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (EASTERN_EDGE) THEN
        DO k=0,N(ng)
          DO j=Jstr,Jend
            A(Lm(ng)+1,j,k)=A(Lm(ng),j,k)
          END DO
        END DO
      END IF
      IF (WESTERN_EDGE) THEN
        DO k=0,N(ng)
          DO j=Jstr,Jend
            A(0,j,k)=A(1,j,k)
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
      IF (NORTHERN_EDGE) THEN
        DO k=0,N(ng)
          DO i=Istr,Iend
            A(i,Mm(ng)+1,k)=A(i,Mm(ng),k)
          END DO
        END DO
      END IF
      IF (SOUTHERN_EDGE) THEN
        DO k=0,N(ng)
          DO i=Istr,Iend
            A(i,0,k)=A(i,1,k)
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
      IF (SOUTHERN_EDGE.and.WESTERN_EDGE) THEN
        DO k=0,N(ng)
          A(0,0,k)=0.5_r8*(A(1,0,k)+A(0,1,k))
        END DO
      END IF
      IF (SOUTHERN_EDGE.and.EASTERN_EDGE) THEN
        DO k=0,N(ng)
          A(Lm(ng)+1,0,k)=0.5_r8*(A(Lm(ng)+1,1,k)+A(Lm(ng),0,k))
        END DO
      END IF
      IF (NORTHERN_EDGE.and.WESTERN_EDGE) THEN
        DO k=0,N(ng)
          A(0,Mm(ng)+1,k)=0.5_r8*(A(0,Mm(ng),k)+A(1,Mm(ng)+1,k))
        END DO
      END IF
      IF (NORTHERN_EDGE.and.EASTERN_EDGE) THEN
        DO k=0,N(ng)
          A(Lm(ng)+1,Mm(ng)+1,k)=0.5_r8*(A(Lm(ng)+1,Mm(ng)  ,k)+        &
     &                                   A(Lm(ng)  ,Mm(ng)+1,k))
        END DO
      END IF
# endif
# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Exchange boundary data.
!-----------------------------------------------------------------------
!
      CALL exchange_w3d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, A)
# endif
      RETURN
      END SUBROUTINE bc_w3d_tile
#endif
      END MODULE bc_3d_mod
