#include "cppdefs.h"
      MODULE exchange_3d_mod
#ifdef SOLVE3D
!
!========================================== Alexander F. Shchepetkin ===
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This package contains periodic boundary conditions and parallel     !
!  exchage (distributed-memory only) routines for 3D variables.        !
!                                                                      !
!  Routines:                                                           !
!                                                                      !
!    exchange_p3d_tile    periodic conditions/exchange at PSI-points   !
!    exchange_r3d_tile    periodic conditions/exchange at RHO-points   !
!    exchange_u3d_tile    periodic conditions/exchange at U-points     !
!    exchange_v3d_tile    periodic conditions/exchange at V-points     !
!    exchange_w3d_tile    periodic conditions/exchange at W-points     !
!                                                                      !
!=======================================================================
!
      implicit none

      CONTAINS
!
!***********************************************************************
      SUBROUTINE exchange_p3d_tile (ng, Istr, Iend, Jstr, Jend,         &
     &                              LBi, UBi, LBj, UBj, LBk, UBk, A)
!***********************************************************************
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk

# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: A(LBi:,LBj:,LBk:)
# else
      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj,LBk:UBk)
# endif
!
!  Local variable declarations.
!
# ifdef DISTRIBUTE
      logical :: EWperiodic, NSperiodic
# endif
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j, k

# include "set_bounds.h"
!
# ifdef DISTRIBUTE
      EWperiodic=.false.
      NSperiodic=.false.
# endif

# ifdef EW_PERIODIC
#  ifdef NS_PERIODIC
#   define J_RANGE Jstr,Jend
#  else
#   define J_RANGE Jstr,JendR
#  endif
!
!-----------------------------------------------------------------------
!  East-West periodic boundary conditions.
!-----------------------------------------------------------------------
!
#  ifdef DISTRIBUTE
      EWperiodic=.true.
      IF (NtileI(ng).eq.1) THEN
#  endif
        IF (WESTERN_EDGE) THEN
          DO k=LBk,UBk
            DO j=J_RANGE
              A(Lm(ng)+1,j,k)=A(1,j,k)
              A(Lm(ng)+2,j,k)=A(2,j,k)
            END DO
          END DO
        END IF
        IF (EASTERN_EDGE) THEN
          DO k=LBk,UBk
            DO j=J_RANGE
              A(-2,j,k)=A(Lm(ng)-2,j,k)
              A(-1,j,k)=A(Lm(ng)-1,j,k)
              A( 0,j,k)=A(Lm(ng)  ,j,k)
            END DO
          END DO
        END IF
#  ifdef DISTRIBUTE
      END IF      
#  endif
#  undef J_RANGE
# endif

# ifdef NS_PERIODIC
#  ifdef EW_PERIODIC
#   define I_RANGE Istr,Iend
#  else
#   define I_RANGE Istr,IendR
#  endif
!
!-----------------------------------------------------------------------
!  North-South periodic boundary conditions.
!-----------------------------------------------------------------------
!
#  ifdef DISTRIBUTE
      NSperiodic=.true.
      IF (NtileJ(ng).eq.1) THEN
#  endif
        IF (SOUTHERN_EDGE) THEN
          DO k=LBk,UBk
            DO i=I_RANGE
              A(i,Mm(ng)+1,k)=A(i,1,k)
              A(i,Mm(ng)+2,k)=A(i,2,k)
            END DO
          END DO
        END IF
        IF (NORTHERN_EDGE) THEN
          DO k=LBk,UBk
            DO i=I_RANGE
              A(i,-2,k)=A(i,Mm(ng)-2,k)
              A(i,-1,k)=A(i,Mm(ng)-1,k)
              A(i, 0,k)=A(i,Mm(ng)  ,k)
            END DO
          END DO
        END IF
#  ifdef DISTRIBUTE
      END IF      
#  endif
#  undef I_RANGE
# endif

# if defined EW_PERIODIC && defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
#  ifdef DISTRIBUTE
      IF ((NtileI(ng).eq.1).and.(NtileJ(ng).eq.1)) THEN
#  endif
        IF (WESTERN_EDGE.and.SOUTHERN_EDGE) THEN
          DO k=LBk,UBk
            A(Lm(ng)+1,Mm(ng)+1,k)=A(1,1,k)
            A(Lm(ng)+1,Mm(ng)+2,k)=A(1,2,k)
            A(Lm(ng)+2,Mm(ng)+1,k)=A(2,1,k)
            A(Lm(ng)+2,Mm(ng)+2,k)=A(2,2,k)
          END DO
        END IF
        IF (EASTERN_EDGE.and.SOUTHERN_EDGE) THEN
          DO k=LBk,UBk
            A(-2,Mm(ng)+1,k)=A(Lm(ng)-2,1,k)
            A(-1,Mm(ng)+1,k)=A(Lm(ng)-1,1,k)
            A( 0,Mm(ng)+1,k)=A(Lm(ng)  ,1,k)
            A(-2,Mm(ng)+2,k)=A(Lm(ng)-2,2,k)
            A(-1,Mm(ng)+2,k)=A(Lm(ng)-1,2,k)
            A( 0,Mm(ng)+2,k)=A(Lm(ng)  ,2,k)
          END DO
        END IF
        IF (WESTERN_EDGE.and.NORTHERN_EDGE) THEN
          DO k=LBk,UBk
            A(Lm(ng)+1,-2,k)=A(1,Mm(ng)-2,k)
            A(Lm(ng)+1,-1,k)=A(1,Mm(ng)-1,k)
            A(Lm(ng)+1, 0,k)=A(1,Mm(ng)  ,k)
            A(Lm(ng)+2,-2,k)=A(2,Mm(ng)-2,k)
            A(Lm(ng)+2,-1,k)=A(2,Mm(ng)-1,k)
            A(Lm(ng)+2, 0,k)=A(2,Mm(ng)  ,k)
          END DO
        END IF
        IF (EASTERN_EDGE.and.NORTHERN_EDGE) THEN
          DO k=LBk,UBk
            A(-2,-2,k)=A(Lm(ng)-2,Mm(ng)-2,k)
            A(-2,-1,k)=A(Lm(ng)-2,Mm(ng)-1,k)
            A(-2, 0,k)=A(Lm(ng)-2,Mm(ng)  ,k)
            A(-1,-2,k)=A(Lm(ng)-1,Mm(ng)-2,k)
            A(-1,-1,k)=A(Lm(ng)-1,Mm(ng)-1,k)
            A(-1, 0,k)=A(Lm(ng)-1,Mm(ng)  ,k)
            A( 0,-2,k)=A(Lm(ng)  ,Mm(ng)-2,k)
            A( 0,-1,k)=A(Lm(ng)  ,Mm(ng)-1,k)
            A( 0, 0,k)=A(Lm(ng)  ,Mm(ng)  ,k)
          END DO
        END IF
#  ifdef DISTRIBUTE
      END IF      
#  endif
# endif
# ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Send/recieve tile boundary data.
!-----------------------------------------------------------------------
!
      CALL mp_exchange (ng, Istr, Iend, Jstr, Jend,                     &
     &                  LBi, UBi, LBj, UBj, LBk, UBk,                   &
     &                  GHOST_POINTS, EWperiodic, NSperiodic, A)
# endif
      RETURN
      END SUBROUTINE exchange_p3d_tile

!
!***********************************************************************
      SUBROUTINE exchange_r3d_tile (ng, Istr, Iend, Jstr, Jend,         &
     &                              LBi, UBi, LBj, UBj, LBk, UBk, A)
!***********************************************************************
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk

# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: A(LBi:,LBj:,LBk:)
# else
      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj,LBk:UBk)
# endif
!
!  Local variable declarations.
!
# ifdef DISTRIBUTE
      logical :: EWperiodic, NSperiodic
# endif
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j, k

# include "set_bounds.h"
!
# ifdef DISTRIBUTE
      EWperiodic=.false.
      NSperiodic=.false.
# endif

# ifdef EW_PERIODIC
#  ifdef NS_PERIODIC
#   define J_RANGE Jstr,Jend
#  else
#   define J_RANGE JstrR,JendR
#  endif
!
!-----------------------------------------------------------------------
!  East-West periodic boundary conditions.
!-----------------------------------------------------------------------
!
#  ifdef DISTRIBUTE
      EWperiodic=.true.
      IF (NtileI(ng).eq.1) THEN
#  endif
        IF (WESTERN_EDGE) THEN
          DO k=LBk,UBk
            DO j=J_RANGE
              A(Lm(ng)+1,j,k)=A(1,j,k)
              A(Lm(ng)+2,j,k)=A(2,j,k)
            END DO
          END DO
        END IF
        IF (EASTERN_EDGE) THEN
          DO k=LBk,UBk
            DO j=J_RANGE
              A(-2,j,k)=A(Lm(ng)-2,j,k)
              A(-1,j,k)=A(Lm(ng)-1,j,k)
              A( 0,j,k)=A(Lm(ng)  ,j,k)
            END DO
          END DO
        END IF
#  ifdef DISTRIBUTE
      END IF      
#  endif
#  undef J_RANGE
# endif

# ifdef NS_PERIODIC
#  ifdef EW_PERIODIC
#   define I_RANGE Istr,Iend
#  else
#   define I_RANGE IstrR,IendR
#  endif
!
!-----------------------------------------------------------------------
!  North-South periodic boundary conditions.
!-----------------------------------------------------------------------
!
#  ifdef DISTRIBUTE
      NSperiodic=.true.
      IF (NtileJ(ng).eq.1) THEN
#  endif
        IF (SOUTHERN_EDGE) THEN
          DO k=LBk,UBk
            DO i=I_RANGE
              A(i,Mm(ng)+1,k)=A(i,1,k)
              A(i,Mm(ng)+2,k)=A(i,2,k)
            END DO
          END DO
        END IF
        IF (NORTHERN_EDGE) THEN
          DO k=LBk,UBk
            DO i=I_RANGE
              A(i,-2,k)=A(i,Mm(ng)-2,k)
              A(i,-1,k)=A(i,Mm(ng)-1,k)
              A(i, 0,k)=A(i,Mm(ng)  ,k)
            END DO
          END DO
        END IF
#  ifdef DISTRIBUTE
      END IF      
#  endif
#  undef I_RANGE
# endif
# if defined EW_PERIODIC && defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
#  ifdef DISTRIBUTE
      IF ((NtileI(ng).eq.1).and.(NtileJ(ng).eq.1)) THEN
#  endif
        IF (WESTERN_EDGE.and.SOUTHERN_EDGE) THEN
          DO k=LBk,UBk
            A(Lm(ng)+1,Mm(ng)+1,k)=A(1,1,k)
            A(Lm(ng)+1,Mm(ng)+2,k)=A(1,2,k)
            A(Lm(ng)+2,Mm(ng)+1,k)=A(2,1,k)
            A(Lm(ng)+2,Mm(ng)+2,k)=A(2,2,k)
          END DO
        END IF
        IF (EASTERN_EDGE.and.SOUTHERN_EDGE) THEN
          DO k=LBk,UBk
            A(-2,Mm(ng)+1,k)=A(Lm(ng)-2,1,k)
            A(-1,Mm(ng)+1,k)=A(Lm(ng)-1,1,k)
            A( 0,Mm(ng)+1,k)=A(Lm(ng)  ,1,k)
            A(-2,Mm(ng)+2,k)=A(Lm(ng)-2,2,k)
            A(-1,Mm(ng)+2,k)=A(Lm(ng)-1,2,k)
            A( 0,Mm(ng)+2,k)=A(Lm(ng)  ,2,k)
          END DO
        END IF
        IF (WESTERN_EDGE.and.NORTHERN_EDGE) THEN
          DO k=LBk,UBk
            A(Lm(ng)+1,-2,k)=A(1,Mm(ng)-2,k)
            A(Lm(ng)+1,-1,k)=A(1,Mm(ng)-1,k)
            A(Lm(ng)+1, 0,k)=A(1,Mm(ng)  ,k)
            A(Lm(ng)+2,-2,k)=A(2,Mm(ng)-2,k)
            A(Lm(ng)+2,-1,k)=A(2,Mm(ng)-1,k)
            A(Lm(ng)+2, 0,k)=A(2,Mm(ng)  ,k)
          END DO
        END IF
        IF (EASTERN_EDGE.and.NORTHERN_EDGE) THEN
          DO k=LBk,UBk
            A(-2,-2,k)=A(Lm(ng)-2,Mm(ng)-2,k)
            A(-2,-1,k)=A(Lm(ng)-2,Mm(ng)-1,k)
            A(-2, 0,k)=A(Lm(ng)-2,Mm(ng)  ,k)
            A(-1,-2,k)=A(Lm(ng)-1,Mm(ng)-2,k)
            A(-1,-1,k)=A(Lm(ng)-1,Mm(ng)-1,k)
            A(-1, 0,k)=A(Lm(ng)-1,Mm(ng)  ,k)
            A( 0,-2,k)=A(Lm(ng)  ,Mm(ng)-2,k)
            A( 0,-1,k)=A(Lm(ng)  ,Mm(ng)-1,k)
            A( 0, 0,k)=A(Lm(ng)  ,Mm(ng)  ,k)
          END DO
        END IF
#  ifdef DISTRIBUTE
      END IF      
#  endif
# endif
# ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Send/recieve tile boundary data.
!-----------------------------------------------------------------------
!
      CALL mp_exchange (ng, Istr, Iend, Jstr, Jend,                     &
     &                  LBi, UBi, LBj, UBj, LBk, UBk,                   &
     &                  GHOST_POINTS, EWperiodic, NSperiodic, A)
# endif
      RETURN
      END SUBROUTINE exchange_r3d_tile

!
!***********************************************************************
      SUBROUTINE exchange_u3d_tile (ng, Istr, Iend, Jstr, Jend,         &
     &                              LBi, UBi, LBj, UBj, LBk, UBk, A)
!***********************************************************************
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk

# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: A(LBi:,LBj:,LBk:)
# else
      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj,LBk:UBk)
# endif
!
!  Local variable declarations.
!
# ifdef DISTRIBUTE
      logical :: EWperiodic, NSperiodic
# endif
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j, k

# include "set_bounds.h"
!
# ifdef DISTRIBUTE
      EWperiodic=.false.
      NSperiodic=.false.
# endif

# ifdef EW_PERIODIC
#  ifdef NS_PERIODIC
#   define J_RANGE Jstr,Jend
#  else
#   define J_RANGE JstrR,JendR
#  endif
!
!-----------------------------------------------------------------------
!  East-West periodic boundary conditions.
!-----------------------------------------------------------------------
!
#  ifdef DISTRIBUTE
      EWperiodic=.true.
      IF (NtileI(ng).eq.1) THEN
#  endif
        IF (WESTERN_EDGE) THEN
          DO k=LBk,UBk
            DO j=J_RANGE
              A(Lm(ng)+1,j,k)=A(1,j,k)
              A(Lm(ng)+2,j,k)=A(2,j,k)
            END DO
          END DO
        END IF
        IF (EASTERN_EDGE) THEN
          DO k=LBk,UBk
            DO j=J_RANGE
              A(-2,j,k)=A(Lm(ng)-2,j,k)
              A(-1,j,k)=A(Lm(ng)-1,j,k)
              A( 0,j,k)=A(Lm(ng)  ,j,k)
            END DO
          END DO
        END IF
#  ifdef DISTRIBUTE
      END IF      
#  endif
#  undef J_RANGE
# endif

# ifdef NS_PERIODIC
#  ifdef EW_PERIODIC
#   define I_RANGE Istr,Iend
#  else
#   define I_RANGE Istr,IendR
#  endif
!
!-----------------------------------------------------------------------
!  North-South periodic boundary conditions.
!-----------------------------------------------------------------------
!
#  ifdef DISTRIBUTE
      NSperiodic=.true.
      IF (NtileJ(ng).eq.1) THEN
#  endif
        IF (SOUTHERN_EDGE) THEN
          DO k=LBk,UBk
            DO i=I_RANGE
              A(i,Mm(ng)+1,k)=A(i,1,k)
              A(i,Mm(ng)+2,k)=A(i,2,k)
            END DO
          END DO
        END IF
        IF (NORTHERN_EDGE) THEN
          DO k=LBk,UBk
            DO i=I_RANGE
              A(i,-2,k)=A(i,Mm(ng)-2,k)
              A(i,-1,k)=A(i,Mm(ng)-1,k)
              A(i, 0,k)=A(i,Mm(ng)  ,k)
            END DO
          END DO
        END IF
#  ifdef DISTRIBUTE
      END IF      
#  endif
#  undef I_RANGE
# endif

# if defined EW_PERIODIC && defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
#  ifdef DISTRIBUTE
      IF ((NtileI(ng).eq.1).and.(NtileJ(ng).eq.1)) THEN
#  endif
        IF (WESTERN_EDGE.and.SOUTHERN_EDGE) THEN
          DO k=LBk,UBk
            A(Lm(ng)+1,Mm(ng)+1,k)=A(1,1,k)
            A(Lm(ng)+1,Mm(ng)+2,k)=A(1,2,k)
            A(Lm(ng)+2,Mm(ng)+1,k)=A(2,1,k)
            A(Lm(ng)+2,Mm(ng)+2,k)=A(2,2,k)
          END DO
        END IF
        IF (EASTERN_EDGE.and.SOUTHERN_EDGE) THEN
          DO k=LBk,UBk
            A(-2,Mm(ng)+1,k)=A(Lm(ng)-2,1,k)
            A(-1,Mm(ng)+1,k)=A(Lm(ng)-1,1,k)
            A( 0,Mm(ng)+1,k)=A(Lm(ng)  ,1,k)
            A(-2,Mm(ng)+2,k)=A(Lm(ng)-2,2,k)
            A(-1,Mm(ng)+2,k)=A(Lm(ng)-1,2,k)
            A( 0,Mm(ng)+2,k)=A(Lm(ng)  ,2,k)
          END DO
        END IF
        IF (WESTERN_EDGE.and.NORTHERN_EDGE) THEN
          DO k=LBk,UBk
            A(Lm(ng)+1,-2,k)=A(1,Mm(ng)-2,k)
            A(Lm(ng)+1,-1,k)=A(1,Mm(ng)-1,k)
            A(Lm(ng)+1, 0,k)=A(1,Mm(ng)  ,k)
            A(Lm(ng)+2,-2,k)=A(2,Mm(ng)-2,k)
            A(Lm(ng)+2,-1,k)=A(2,Mm(ng)-1,k)
            A(Lm(ng)+2, 0,k)=A(2,Mm(ng)  ,k)
          END DO
        END IF
        IF (EASTERN_EDGE.and.NORTHERN_EDGE) THEN
          DO k=LBk,UBk
            A(-2,-2,k)=A(Lm(ng)-2,Mm(ng)-2,k)
            A(-2,-1,k)=A(Lm(ng)-2,Mm(ng)-1,k)
            A(-2, 0,k)=A(Lm(ng)-2,Mm(ng)  ,k)
            A(-1,-2,k)=A(Lm(ng)-1,Mm(ng)-2,k)
            A(-1,-1,k)=A(Lm(ng)-1,Mm(ng)-1,k)
            A(-1, 0,k)=A(Lm(ng)-1,Mm(ng)  ,k)
            A( 0,-2,k)=A(Lm(ng)  ,Mm(ng)-2,k)
            A( 0,-1,k)=A(Lm(ng)  ,Mm(ng)-1,k)
            A( 0, 0,k)=A(Lm(ng)  ,Mm(ng)  ,k)
          END DO
        END IF
#  ifdef DISTRIBUTE
      END IF      
#  endif
# endif
# ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Send/recieve tile boundary data.
!-----------------------------------------------------------------------
!
      CALL mp_exchange (ng, Istr, Iend, Jstr, Jend,                     &
     &                  LBi, UBi, LBj, UBj, LBk, UBk,                   &
     &                  GHOST_POINTS, EWperiodic, NSperiodic, A)
# endif
      RETURN
      END SUBROUTINE exchange_u3d_tile

!
!***********************************************************************
      SUBROUTINE exchange_v3d_tile (ng, Istr, Iend, Jstr, Jend,         &
     &                              LBi, UBi, LBj, UBj, LBk, UBk, A)
!***********************************************************************
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk

# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: A(LBi:,LBj:,LBk:)
# else
      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj,LBk:UBk)
# endif
!
!  Local variable declarations.
!
# ifdef DISTRIBUTE
      logical :: EWperiodic, NSperiodic
# endif
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j, k

# include "set_bounds.h"
!
# ifdef DISTRIBUTE
      EWperiodic=.false.
      NSperiodic=.false.
# endif

# ifdef EW_PERIODIC
#  ifdef NS_PERIODIC
#   define J_RANGE Jstr,Jend
#  else
#   define J_RANGE Jstr,JendR
#  endif
!
!-----------------------------------------------------------------------
!  East-West periodic boundary conditions.
!-----------------------------------------------------------------------
!
#  ifdef DISTRIBUTE
      EWperiodic=.true.
      IF (NtileI(ng).eq.1) THEN
#  endif
        IF (WESTERN_EDGE) THEN
          DO k=LBk,UBk
            DO j=J_RANGE
              A(Lm(ng)+1,j,k)=A(1,j,k)
              A(Lm(ng)+2,j,k)=A(2,j,k)
            END DO
          END DO
        END IF
        IF (EASTERN_EDGE) THEN
          DO k=LBk,UBk
            DO j=J_RANGE
              A(-2,j,k)=A(Lm(ng)-2,j,k)
              A(-1,j,k)=A(Lm(ng)-1,j,k)
              A( 0,j,k)=A(Lm(ng)  ,j,k)
            END DO
          END DO
        END IF
#  ifdef DISTRIBUTE
      END IF      
#  endif
#  undef J_RANGE
# endif

# ifdef NS_PERIODIC
#  ifdef EW_PERIODIC
#   define I_RANGE Istr,Iend
#  else
#   define I_RANGE IstrR,IendR
#  endif
!
!-----------------------------------------------------------------------
!  North-South periodic boundary conditions.
!-----------------------------------------------------------------------
!
#  ifdef DISTRIBUTE
      NSperiodic=.true.
      IF (NtileJ(ng).eq.1) THEN
#  endif
        IF (SOUTHERN_EDGE) THEN
          DO k=LBk,UBk
            DO i=I_RANGE
              A(i,Mm(ng)+1,k)=A(i,1,k)
              A(i,Mm(ng)+2,k)=A(i,2,k)
            END DO
          END DO
        END IF
        IF (NORTHERN_EDGE) THEN
          DO k=LBk,UBk
            DO i=I_RANGE
              A(i,-2,k)=A(i,Mm(ng)-2,k)
              A(i,-1,k)=A(i,Mm(ng)-1,k)
              A(i, 0,k)=A(i,Mm(ng)  ,k)
            END DO
          END DO
        END IF
#  ifdef DISTRIBUTE
      END IF      
#  endif
#  undef I_RANGE
# endif

# if defined EW_PERIODIC && defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
#  ifdef DISTRIBUTE
      IF ((NtileI(ng).eq.1).and.(NtileJ(ng).eq.1)) THEN
#  endif
        IF (WESTERN_EDGE.and.SOUTHERN_EDGE) THEN
          DO k=LBk,UBk
            A(Lm(ng)+1,Mm(ng)+1,k)=A(1,1,k)
            A(Lm(ng)+1,Mm(ng)+2,k)=A(1,2,k)
            A(Lm(ng)+2,Mm(ng)+1,k)=A(2,1,k)
            A(Lm(ng)+2,Mm(ng)+2,k)=A(2,2,k)
          END DO
        END IF
        IF (EASTERN_EDGE.and.SOUTHERN_EDGE) THEN
          DO k=LBk,UBk
            A(-2,Mm(ng)+1,k)=A(Lm(ng)-2,1,k)
            A(-1,Mm(ng)+1,k)=A(Lm(ng)-1,1,k)
            A( 0,Mm(ng)+1,k)=A(Lm(ng)  ,1,k)
            A(-2,Mm(ng)+2,k)=A(Lm(ng)-2,2,k)
            A(-1,Mm(ng)+2,k)=A(Lm(ng)-1,2,k)
            A( 0,Mm(ng)+2,k)=A(Lm(ng)  ,2,k)
          END DO
        END IF
        IF (WESTERN_EDGE.and.NORTHERN_EDGE) THEN
          DO k=LBk,UBk
            A(Lm(ng)+1,-2,k)=A(1,Mm(ng)-2,k)
            A(Lm(ng)+1,-1,k)=A(1,Mm(ng)-1,k)
            A(Lm(ng)+1, 0,k)=A(1,Mm(ng)  ,k)
            A(Lm(ng)+2,-2,k)=A(2,Mm(ng)-2,k)
            A(Lm(ng)+2,-1,k)=A(2,Mm(ng)-1,k)
            A(Lm(ng)+2, 0,k)=A(2,Mm(ng)  ,k)
          END DO
        END IF
        IF (EASTERN_EDGE.and.NORTHERN_EDGE) THEN
          DO k=LBk,UBk
            A(-2,-2,k)=A(Lm(ng)-2,Mm(ng)-2,k)
            A(-2,-1,k)=A(Lm(ng)-2,Mm(ng)-1,k)
            A(-2, 0,k)=A(Lm(ng)-2,Mm(ng)  ,k)
            A(-1,-2,k)=A(Lm(ng)-1,Mm(ng)-2,k)
            A(-1,-1,k)=A(Lm(ng)-1,Mm(ng)-1,k)
            A(-1, 0,k)=A(Lm(ng)-1,Mm(ng)  ,k)
            A( 0,-2,k)=A(Lm(ng)  ,Mm(ng)-2,k)
            A( 0,-1,k)=A(Lm(ng)  ,Mm(ng)-1,k)
            A( 0, 0,k)=A(Lm(ng)  ,Mm(ng)  ,k)
          END DO
        END IF
#  ifdef DISTRIBUTE
      END IF      
#  endif
# endif
# ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Send/recieve tile boundary data.
!-----------------------------------------------------------------------
!
      CALL mp_exchange (ng, Istr, Iend, Jstr, Jend,                     &
     &                  LBi, UBi, LBj, UBj, LBk, UBk,                   &
     &                  GHOST_POINTS, EWperiodic, NSperiodic, A)
# endif
      RETURN
      END SUBROUTINE exchange_v3d_tile

!
!***********************************************************************
      SUBROUTINE exchange_w3d_tile (ng, Istr, Iend, Jstr, Jend,         &
     &                              LBi, UBi, LBj, UBj, LBk, UBk, A)
!***********************************************************************
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk

# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: A(LBi:,LBj:,LBk:)
# else
      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj,LBk:UBk)
# endif
!
!  Local variable declarations.
!
# ifdef DISTRIBUTE
      logical :: EWperiodic, NSperiodic
# endif
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j, k

# include "set_bounds.h"
!
# ifdef DISTRIBUTE
      EWperiodic=.false.
      NSperiodic=.false.
# endif

# ifdef EW_PERIODIC
#  ifdef NS_PERIODIC
#   define J_RANGE Jstr,Jend
#  else
#   define J_RANGE JstrR,JendR
#  endif
!
!-----------------------------------------------------------------------
!  East-West periodic boundary conditions.
!-----------------------------------------------------------------------
!
#  ifdef DISTRIBUTE
      EWperiodic=.true.
      IF (NtileI(ng).eq.1) THEN
#  endif
        IF (WESTERN_EDGE) THEN
          DO k=LBk,UBk
            DO j=J_RANGE
              A(Lm(ng)+1,j,k)=A(1,j,k)
              A(Lm(ng)+2,j,k)=A(2,j,k)
            END DO
          END DO
        END IF
        IF (EASTERN_EDGE) THEN
          DO k=LBk,UBk
            DO j=J_RANGE
              A(-2,j,k)=A(Lm(ng)-2,j,k)
              A(-1,j,k)=A(Lm(ng)-1,j,k)
              A( 0,j,k)=A(Lm(ng)  ,j,k)
            END DO
          END DO
        END IF
#  ifdef DISTRIBUTE
      END IF      
#  endif
#  undef J_RANGE
# endif

# ifdef NS_PERIODIC
#  ifdef EW_PERIODIC
#   define I_RANGE Istr,Iend
#  else
#   define I_RANGE IstrR,IendR
#  endif
!
!-----------------------------------------------------------------------
!  North-South periodic boundary conditions.
!-----------------------------------------------------------------------
!
#  ifdef DISTRIBUTE
      NSperiodic=.true.
      IF (NtileJ(ng).eq.1) THEN
#  endif
        IF (SOUTHERN_EDGE) THEN
          DO k=LBk,UBk
            DO i=I_RANGE
              A(i,Mm(ng)+1,k)=A(i,1,k)
              A(i,Mm(ng)+2,k)=A(i,2,k)
            END DO
          END DO
        END IF
        IF (NORTHERN_EDGE) THEN
          DO k=LBk,UBk
            DO i=I_RANGE
              A(i,-2,k)=A(i,Mm(ng)-2,k)
              A(i,-1,k)=A(i,Mm(ng)-1,k)
              A(i, 0,k)=A(i,Mm(ng)  ,k)
            END DO
          END DO
        END IF
#  ifdef DISTRIBUTE
      END IF      
#  endif
#  undef I_RANGE
# endif
# if defined EW_PERIODIC && defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
#  ifdef DISTRIBUTE
      IF ((NtileI(ng).eq.1).and.(NtileJ(ng).eq.1)) THEN
#  endif
        IF (WESTERN_EDGE.and.SOUTHERN_EDGE) THEN
          DO k=LBk,UBk
            A(Lm(ng)+1,Mm(ng)+1,k)=A(1,1,k)
            A(Lm(ng)+1,Mm(ng)+2,k)=A(1,2,k)
            A(Lm(ng)+2,Mm(ng)+1,k)=A(2,1,k)
            A(Lm(ng)+2,Mm(ng)+2,k)=A(2,2,k)
          END DO
        END IF
        IF (EASTERN_EDGE.and.SOUTHERN_EDGE) THEN
          DO k=LBk,UBk
            A(-2,Mm(ng)+1,k)=A(Lm(ng)-2,1,k)
            A(-1,Mm(ng)+1,k)=A(Lm(ng)-1,1,k)
            A( 0,Mm(ng)+1,k)=A(Lm(ng)  ,1,k)
            A(-2,Mm(ng)+2,k)=A(Lm(ng)-2,2,k)
            A(-1,Mm(ng)+2,k)=A(Lm(ng)-1,2,k)
            A( 0,Mm(ng)+2,k)=A(Lm(ng)  ,2,k)
          END DO
        END IF
        IF (WESTERN_EDGE.and.NORTHERN_EDGE) THEN
          DO k=LBk,UBk
            A(Lm(ng)+1,-2,k)=A(1,Mm(ng)-2,k)
            A(Lm(ng)+1,-1,k)=A(1,Mm(ng)-1,k)
            A(Lm(ng)+1, 0,k)=A(1,Mm(ng)  ,k)
            A(Lm(ng)+2,-2,k)=A(2,Mm(ng)-2,k)
            A(Lm(ng)+2,-1,k)=A(2,Mm(ng)-1,k)
            A(Lm(ng)+2, 0,k)=A(2,Mm(ng)  ,k)
          END DO
        END IF
        IF (EASTERN_EDGE.and.NORTHERN_EDGE) THEN
          DO k=LBk,UBk
            A(-2,-2,k)=A(Lm(ng)-2,Mm(ng)-2,k)
            A(-2,-1,k)=A(Lm(ng)-2,Mm(ng)-1,k)
            A(-2, 0,k)=A(Lm(ng)-2,Mm(ng)  ,k)
            A(-1,-2,k)=A(Lm(ng)-1,Mm(ng)-2,k)
            A(-1,-1,k)=A(Lm(ng)-1,Mm(ng)-1,k)
            A(-1, 0,k)=A(Lm(ng)-1,Mm(ng)  ,k)
            A( 0,-2,k)=A(Lm(ng)  ,Mm(ng)-2,k)
            A( 0,-1,k)=A(Lm(ng)  ,Mm(ng)-1,k)
            A( 0, 0,k)=A(Lm(ng)  ,Mm(ng)  ,k)
          END DO
        END IF
#  ifdef DISTRIBUTE
      END IF      
#  endif
# endif
# ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Send/recieve tile boundary data.
!-----------------------------------------------------------------------
!
      CALL mp_exchange (ng, Istr, Iend, Jstr, Jend,                     &
     &                  LBi, UBi, LBj, UBj, LBk, UBk,                   &
     &                  GHOST_POINTS, EWperiodic, NSperiodic, A)
# endif
      RETURN
      END SUBROUTINE exchange_w3d_tile
#endif

      END MODULE exchange_3d_mod
