#include "cppdefs.h"
      MODULE exchange_2d_mod
#if defined EW_PERIODIC || defined NS_PERIODIC
!
!svn $Id$
!=======================================================================
!  Copyright (c) 2002-2011 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  These routines apply periodic boundary conditions to generic        !
!  2D fields.                                                          !
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
!     A                 Processed 2D field                             !
!                                                                      !
!  Routines:                                                           !
!                                                                      !
!    exchange_p2d_tile  periodic conditions/exchange at PSI-points     !
!    exchange_r2d_tile  periodic conditions/exchange at RHO-points     !
!    exchange_u2d_tile  periodic conditions/exchange at U-points       !
!    exchange_v2d_tile  periodic conditions/exchange at V-points       !
!                                                                      !
!=======================================================================
!
      implicit none

      CONTAINS
!
!***********************************************************************
      SUBROUTINE exchange_p2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              A)
!***********************************************************************
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: A(LBi:,LBj:)
# else
      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
# ifdef EW_PERIODIC
      logical :: EW_exchange
# endif
# ifdef NS_PERIODIC
      logical :: NS_exchange
# endif
      integer :: i, j

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Determine processing switches.
!-----------------------------------------------------------------------
!
# ifdef EW_PERIODIC
#  ifdef DISTRIBUTE
      EW_exchange=NtileI(ng).eq.1
#  else
      EW_exchange=.TRUE.
#  endif
# endif
# ifdef NS_PERIODIC
#  ifdef DISTRIBUTE
      NS_exchange=NtileJ(ng).eq.1
#  else
      NS_exchange=.TRUE.
#  endif
# endif

# ifdef EW_PERIODIC
!
!-----------------------------------------------------------------------
!  East-West periodic boundary conditions.
!-----------------------------------------------------------------------
!
#  ifdef NS_PERIODIC
#   define J_RANGE Jstr,Jend
#  else
#   define J_RANGE Jstr,JendR
#  endif

      IF (EW_exchange) THEN
        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          DO j=J_RANGE
            A(Lm(ng)+1,j)=A(1,j)
            A(Lm(ng)+2,j)=A(2,j)
#  ifdef THREE_GHOST
            A(Lm(ng)+3,j)=A(3,j)
#  endif
          END DO
        END IF
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO j=J_RANGE
            A(-2,j)=A(Lm(ng)-2,j)
            A(-1,j)=A(Lm(ng)-1,j)
            A( 0,j)=A(Lm(ng)  ,j)
          END DO
        END IF
      END IF
#  undef J_RANGE
# endif

# ifdef NS_PERIODIC
!
!-----------------------------------------------------------------------
!  North-South periodic boundary conditions.
!-----------------------------------------------------------------------
!
#  ifdef EW_PERIODIC
#   define I_RANGE Istr,Iend
#  else
#   define I_RANGE Istr,IendR
#  endif

      IF (NS_exchange) THEN
        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO i=I_RANGE
            A(i,Mm(ng)+1)=A(i,1)
            A(i,Mm(ng)+2)=A(i,2)
#  ifdef THREE_GHOST
            A(i,Mm(ng)+3)=A(i,3)
#  endif
          END DO
        END IF
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          DO i=I_RANGE
            A(i,-2)=A(i,Mm(ng)-2)
            A(i,-1)=A(i,Mm(ng)-1)
            A(i, 0)=A(i,Mm(ng)  )
          END DO
        END IF
      END IF
#  undef I_RANGE
# endif

# if defined EW_PERIODIC && defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (EW_exchange.and.NS_exchange) THEN
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          A(Lm(ng)+1,Mm(ng)+1)=A(1,1)
          A(Lm(ng)+1,Mm(ng)+2)=A(1,2)
#  ifdef THREE_GHOST
          A(Lm(ng)+1,Mm(ng)+3)=A(1,3)
#  endif
          A(Lm(ng)+2,Mm(ng)+1)=A(2,1)
          A(Lm(ng)+2,Mm(ng)+2)=A(2,2)
#  ifdef THREE_GHOST
          A(Lm(ng)+2,Mm(ng)+3)=A(2,3)
          A(Lm(ng)+3,Mm(ng)+1)=A(3,1)
          A(Lm(ng)+3,Mm(ng)+2)=A(3,2)
          A(Lm(ng)+3,Mm(ng)+3)=A(3,3)
#  endif
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          A(-2,Mm(ng)+1)=A(Lm(ng)-2,1)
          A(-1,Mm(ng)+1)=A(Lm(ng)-1,1)
          A( 0,Mm(ng)+1)=A(Lm(ng)  ,1)
          A(-2,Mm(ng)+2)=A(Lm(ng)-2,2)
          A(-1,Mm(ng)+2)=A(Lm(ng)-1,2)
          A( 0,Mm(ng)+2)=A(Lm(ng)  ,2)
#  ifdef THREE_GHOST
          A(-2,Mm(ng)+3)=A(Lm(ng)-2,3)
          A(-1,Mm(ng)+3)=A(Lm(ng)-1,3)
          A( 0,Mm(ng)+3)=A(Lm(ng)  ,3)
#  endif
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          A(Lm(ng)+1,-2)=A(1,Mm(ng)-2)
          A(Lm(ng)+1,-1)=A(1,Mm(ng)-1)
          A(Lm(ng)+1, 0)=A(1,Mm(ng)  )
          A(Lm(ng)+2,-2)=A(2,Mm(ng)-2)
          A(Lm(ng)+2,-1)=A(2,Mm(ng)-1)
          A(Lm(ng)+2, 0)=A(2,Mm(ng)  )
#  ifdef THREE_GHOST
          A(Lm(ng)+3,-2)=A(3,Mm(ng)-2)
          A(Lm(ng)+3,-1)=A(3,Mm(ng)-1)
          A(Lm(ng)+3, 0)=A(3,Mm(ng)  )
#  endif
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          A(-2,-2)=A(Lm(ng)-2,Mm(ng)-2)
          A(-2,-1)=A(Lm(ng)-2,Mm(ng)-1)
          A(-2, 0)=A(Lm(ng)-2,Mm(ng)  )
          A(-1,-2)=A(Lm(ng)-1,Mm(ng)-2)
          A(-1,-1)=A(Lm(ng)-1,Mm(ng)-1)
          A(-1, 0)=A(Lm(ng)-1,Mm(ng)  )
          A( 0,-2)=A(Lm(ng)  ,Mm(ng)-2)
          A( 0,-1)=A(Lm(ng)  ,Mm(ng)-1)
          A( 0, 0)=A(Lm(ng)  ,Mm(ng)  )
        END IF
      END IF
# endif
      RETURN
      END SUBROUTINE exchange_p2d_tile

!
!***********************************************************************
      SUBROUTINE exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              A)
!***********************************************************************
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: A(LBi:,LBj:)
# else
      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
# ifdef EW_PERIODIC
      logical :: EW_exchange
# endif
# ifdef NS_PERIODIC
      logical :: NS_exchange
# endif
      integer :: i, j

#  include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Determine processing switches.
!-----------------------------------------------------------------------
!
# ifdef EW_PERIODIC
#  ifdef DISTRIBUTE
      EW_exchange=NtileI(ng).eq.1
#  else
      EW_exchange=.TRUE.
#  endif
# endif
# ifdef NS_PERIODIC
#  ifdef DISTRIBUTE
      NS_exchange=NtileJ(ng).eq.1
#  else
      NS_exchange=.TRUE.
#  endif
# endif

# ifdef EW_PERIODIC
!
!-----------------------------------------------------------------------
!  East-West periodic boundary conditions.
!-----------------------------------------------------------------------
!
#  ifdef NS_PERIODIC
#   define J_RANGE Jstr,Jend
#  else
#   define J_RANGE JstrR,JendR
#  endif

      IF (EW_exchange) THEN
        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          DO j=J_RANGE
            A(Lm(ng)+1,j)=A(1,j)
            A(Lm(ng)+2,j)=A(2,j)
#  ifdef THREE_GHOST
            A(Lm(ng)+3,j)=A(3,j)
#  endif
          END DO
        END IF
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO j=J_RANGE
            A(-2,j)=A(Lm(ng)-2,j)
            A(-1,j)=A(Lm(ng)-1,j)
            A( 0,j)=A(Lm(ng)  ,j)
          END DO
        END IF
      END IF
#  undef J_RANGE
# endif

# ifdef NS_PERIODIC
!
!-----------------------------------------------------------------------
!  North-South periodic boundary conditions.
!-----------------------------------------------------------------------
!
#  ifdef EW_PERIODIC
#   define I_RANGE Istr,Iend
#  else
#   define I_RANGE IstrR,IendR
#  endif

      IF (NS_exchange) THEN
        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO i=I_RANGE
            A(i,Mm(ng)+1)=A(i,1)
            A(i,Mm(ng)+2)=A(i,2)
#  ifdef THREE_GHOST
            A(i,Mm(ng)+3)=A(i,3)
#  endif
          END DO
        END IF
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          DO i=I_RANGE
            A(i,-2)=A(i,Mm(ng)-2)
            A(i,-1)=A(i,Mm(ng)-1)
            A(i, 0)=A(i,Mm(ng)  )
          END DO
        END IF
      END IF
#  undef I_RANGE
# endif

# if defined EW_PERIODIC && defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (EW_exchange.and.NS_exchange) THEN
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          A(Lm(ng)+1,Mm(ng)+1)=A(1,1)
          A(Lm(ng)+1,Mm(ng)+2)=A(1,2)
#  ifdef THREE_GHOST
          A(Lm(ng)+1,Mm(ng)+3)=A(1,3)
#  endif
          A(Lm(ng)+2,Mm(ng)+1)=A(2,1)
          A(Lm(ng)+2,Mm(ng)+2)=A(2,2)
#  ifdef THREE_GHOST
          A(Lm(ng)+2,Mm(ng)+3)=A(2,3)
          A(Lm(ng)+3,Mm(ng)+1)=A(3,1)
          A(Lm(ng)+3,Mm(ng)+2)=A(3,2)
          A(Lm(ng)+3,Mm(ng)+3)=A(3,3)
#  endif
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          A(-2,Mm(ng)+1)=A(Lm(ng)-2,1)
          A(-1,Mm(ng)+1)=A(Lm(ng)-1,1)
          A( 0,Mm(ng)+1)=A(Lm(ng)  ,1)
          A(-2,Mm(ng)+2)=A(Lm(ng)-2,2)
          A(-1,Mm(ng)+2)=A(Lm(ng)-1,2)
          A( 0,Mm(ng)+2)=A(Lm(ng)  ,2)
#  ifdef THREE_GHOST
          A(-2,Mm(ng)+3)=A(Lm(ng)-2,3)
          A(-1,Mm(ng)+3)=A(Lm(ng)-1,3)
          A( 0,Mm(ng)+3)=A(Lm(ng)  ,3)
#  endif
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          A(Lm(ng)+1,-2)=A(1,Mm(ng)-2)
          A(Lm(ng)+1,-1)=A(1,Mm(ng)-1)
          A(Lm(ng)+1, 0)=A(1,Mm(ng)  )
          A(Lm(ng)+2,-2)=A(2,Mm(ng)-2)
          A(Lm(ng)+2,-1)=A(2,Mm(ng)-1)
          A(Lm(ng)+2, 0)=A(2,Mm(ng)  )
#  ifdef THREE_GHOST
          A(Lm(ng)+3,-2)=A(3,Mm(ng)-2)
          A(Lm(ng)+3,-1)=A(3,Mm(ng)-1)
          A(Lm(ng)+3, 0)=A(3,Mm(ng)  )
#  endif
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          A(-2,-2)=A(Lm(ng)-2,Mm(ng)-2)
          A(-2,-1)=A(Lm(ng)-2,Mm(ng)-1)
          A(-2, 0)=A(Lm(ng)-2,Mm(ng)  )
          A(-1,-2)=A(Lm(ng)-1,Mm(ng)-2)
          A(-1,-1)=A(Lm(ng)-1,Mm(ng)-1)
          A(-1, 0)=A(Lm(ng)-1,Mm(ng)  )
          A( 0,-2)=A(Lm(ng)  ,Mm(ng)-2)
          A( 0,-1)=A(Lm(ng)  ,Mm(ng)-1)
          A( 0, 0)=A(Lm(ng)  ,Mm(ng)  )
        END IF
      END IF
# endif
      RETURN
      END SUBROUTINE exchange_r2d_tile

!
!***********************************************************************
      SUBROUTINE exchange_u2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              A)
!***********************************************************************
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: A(LBi:,LBj:)
# else
      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
# ifdef EW_PERIODIC
      logical :: EW_exchange
# endif
# ifdef NS_PERIODIC
      logical :: NS_exchange
# endif
      integer :: i, j

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Determine processing switches.
!-----------------------------------------------------------------------
!
# ifdef EW_PERIODIC
#  ifdef DISTRIBUTE
      EW_exchange=NtileI(ng).eq.1
#  else
      EW_exchange=.TRUE.
#  endif
# endif
# ifdef NS_PERIODIC
#  ifdef DISTRIBUTE
      NS_exchange=NtileJ(ng).eq.1
#  else
      NS_exchange=.TRUE.
#  endif
# endif

# ifdef EW_PERIODIC
!
!-----------------------------------------------------------------------
!  East-West periodic boundary conditions.
!-----------------------------------------------------------------------
!
#  ifdef NS_PERIODIC
#   define J_RANGE Jstr,Jend
#  else
#   define J_RANGE JstrR,JendR
#  endif

      IF (EW_exchange) THEN
        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          DO j=J_RANGE
            A(Lm(ng)+1,j)=A(1,j)
            A(Lm(ng)+2,j)=A(2,j)
#  ifdef THREE_GHOST
            A(Lm(ng)+3,j)=A(3,j)
#  endif
          END DO
        END IF
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO j=J_RANGE
            A(-2,j)=A(Lm(ng)-2,j)
            A(-1,j)=A(Lm(ng)-1,j)
            A( 0,j)=A(Lm(ng)  ,j)
          END DO
        END IF
      END IF
#  undef J_RANGE
# endif

# ifdef NS_PERIODIC
!
!-----------------------------------------------------------------------
!  North-South periodic boundary conditions.
!-----------------------------------------------------------------------
!
#  ifdef EW_PERIODIC
#   define I_RANGE Istr,Iend
#  else
#   define I_RANGE Istr,IendR
#  endif

      IF (NS_exchange) THEN
        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO i=I_RANGE
            A(i,Mm(ng)+1)=A(i,1)
            A(i,Mm(ng)+2)=A(i,2)
#  ifdef THREE_GHOST
            A(i,Mm(ng)+3)=A(i,3)
#  endif
          END DO
        END IF
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          DO i=I_RANGE
            A(i,-2)=A(i,Mm(ng)-2)
            A(i,-1)=A(i,Mm(ng)-1)
            A(i, 0)=A(i,Mm(ng)  )
          END DO
        END IF
      END IF
#  undef I_RANGE
# endif

# if defined EW_PERIODIC && defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (EW_exchange.and.NS_exchange) THEN
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          A(Lm(ng)+1,Mm(ng)+1)=A(1,1)
          A(Lm(ng)+1,Mm(ng)+2)=A(1,2)
#  ifdef THREE_GHOST
          A(Lm(ng)+1,Mm(ng)+3)=A(1,3)
#  endif
          A(Lm(ng)+2,Mm(ng)+1)=A(2,1)
          A(Lm(ng)+2,Mm(ng)+2)=A(2,2)
#  ifdef THREE_GHOST
          A(Lm(ng)+2,Mm(ng)+3)=A(2,3)
          A(Lm(ng)+3,Mm(ng)+1)=A(3,1)
          A(Lm(ng)+3,Mm(ng)+2)=A(3,2)
          A(Lm(ng)+3,Mm(ng)+3)=A(3,3)
#  endif
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          A(-2,Mm(ng)+1)=A(Lm(ng)-2,1)
          A(-1,Mm(ng)+1)=A(Lm(ng)-1,1)
          A( 0,Mm(ng)+1)=A(Lm(ng)  ,1)
          A(-2,Mm(ng)+2)=A(Lm(ng)-2,2)
          A(-1,Mm(ng)+2)=A(Lm(ng)-1,2)
          A( 0,Mm(ng)+2)=A(Lm(ng)  ,2)
#  ifdef THREE_GHOST
          A(-2,Mm(ng)+3)=A(Lm(ng)-2,3)
          A(-1,Mm(ng)+3)=A(Lm(ng)-1,3)
          A( 0,Mm(ng)+3)=A(Lm(ng)  ,3)
#  endif
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          A(Lm(ng)+1,-2)=A(1,Mm(ng)-2)
          A(Lm(ng)+1,-1)=A(1,Mm(ng)-1)
          A(Lm(ng)+1, 0)=A(1,Mm(ng)  )
          A(Lm(ng)+2,-2)=A(2,Mm(ng)-2)
          A(Lm(ng)+2,-1)=A(2,Mm(ng)-1)
          A(Lm(ng)+2, 0)=A(2,Mm(ng)  )
#  ifdef THREE_GHOST
          A(Lm(ng)+3,-2)=A(3,Mm(ng)-2)
          A(Lm(ng)+3,-1)=A(3,Mm(ng)-1)
          A(Lm(ng)+3, 0)=A(3,Mm(ng)  )
#  endif
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          A(-2,-2)=A(Lm(ng)-2,Mm(ng)-2)
          A(-2,-1)=A(Lm(ng)-2,Mm(ng)-1)
          A(-2, 0)=A(Lm(ng)-2,Mm(ng)  )
          A(-1,-2)=A(Lm(ng)-1,Mm(ng)-2)
          A(-1,-1)=A(Lm(ng)-1,Mm(ng)-1)
          A(-1, 0)=A(Lm(ng)-1,Mm(ng)  )
          A( 0,-2)=A(Lm(ng)  ,Mm(ng)-2)
          A( 0,-1)=A(Lm(ng)  ,Mm(ng)-1)
          A( 0, 0)=A(Lm(ng)  ,Mm(ng)  )
        END IF
      END IF
# endif
      RETURN
      END SUBROUTINE exchange_u2d_tile

!
!***********************************************************************
      SUBROUTINE exchange_v2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              A)
!***********************************************************************
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: A(LBi:,LBj:)
# else
      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
# ifdef EW_PERIODIC
      logical :: EW_exchange
# endif
# ifdef NS_PERIODIC
      logical :: NS_exchange
# endif
      integer :: i, j

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Determine processing switches.
!-----------------------------------------------------------------------
!
# ifdef EW_PERIODIC
#  ifdef DISTRIBUTE
      EW_exchange=NtileI(ng).eq.1
#  else
      EW_exchange=.TRUE.
#  endif
# endif
# ifdef NS_PERIODIC
#  ifdef DISTRIBUTE
      NS_exchange=NtileJ(ng).eq.1
#  else
      NS_exchange=.TRUE.
#  endif
# endif

# ifdef EW_PERIODIC
!
!-----------------------------------------------------------------------
!  East-West periodic boundary conditions.
!-----------------------------------------------------------------------
!
#  ifdef NS_PERIODIC
#   define J_RANGE Jstr,Jend
#  else
#   define J_RANGE Jstr,JendR
#  endif

      IF (EW_exchange) THEN
        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          DO j=J_RANGE
            A(Lm(ng)+1,j)=A(1,j)
            A(Lm(ng)+2,j)=A(2,j)
#  ifdef THREE_GHOST
            A(Lm(ng)+3,j)=A(3,j)
#  endif
          END DO
        END IF
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO j=J_RANGE
            A(-2,j)=A(Lm(ng)-2,j)
            A(-1,j)=A(Lm(ng)-1,j)
            A( 0,j)=A(Lm(ng)  ,j)
          END DO
        END IF
      END IF
#  undef J_RANGE
# endif

# ifdef NS_PERIODIC
!
!-----------------------------------------------------------------------
!  North-South periodic boundary conditions.
!-----------------------------------------------------------------------
!
#  ifdef EW_PERIODIC
#   define I_RANGE Istr,Iend
#  else
#   define I_RANGE IstrR,IendR
#  endif

      IF (NS_exchange) THEN
        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO i=I_RANGE
            A(i,Mm(ng)+1)=A(i,1)
            A(i,Mm(ng)+2)=A(i,2)
#  ifdef THREE_GHOST
            A(i,Mm(ng)+3)=A(i,3)
#  endif
          END DO
        END IF
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          DO i=I_RANGE
            A(i,-2)=A(i,Mm(ng)-2)
            A(i,-1)=A(i,Mm(ng)-1)
            A(i, 0)=A(i,Mm(ng)  )
          END DO
        END IF
      END IF
#  undef I_RANGE
# endif

# if defined EW_PERIODIC && defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (EW_exchange.and.NS_exchange) THEN
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          A(Lm(ng)+1,Mm(ng)+1)=A(1,1)
          A(Lm(ng)+1,Mm(ng)+2)=A(1,2)
#  ifdef THREE_GHOST
          A(Lm(ng)+1,Mm(ng)+3)=A(1,3)
#  endif
          A(Lm(ng)+2,Mm(ng)+1)=A(2,1)
          A(Lm(ng)+2,Mm(ng)+2)=A(2,2)
#  ifdef THREE_GHOST
          A(Lm(ng)+2,Mm(ng)+3)=A(2,3)
          A(Lm(ng)+3,Mm(ng)+1)=A(3,1)
          A(Lm(ng)+3,Mm(ng)+2)=A(3,2)
          A(Lm(ng)+3,Mm(ng)+3)=A(3,3)
#  endif
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          A(-2,Mm(ng)+1)=A(Lm(ng)-2,1)
          A(-1,Mm(ng)+1)=A(Lm(ng)-1,1)
          A( 0,Mm(ng)+1)=A(Lm(ng)  ,1)
          A(-2,Mm(ng)+2)=A(Lm(ng)-2,2)
          A(-1,Mm(ng)+2)=A(Lm(ng)-1,2)
          A( 0,Mm(ng)+2)=A(Lm(ng)  ,2)
#  ifdef THREE_GHOST
          A(-2,Mm(ng)+3)=A(Lm(ng)-2,3)
          A(-1,Mm(ng)+3)=A(Lm(ng)-1,3)
          A( 0,Mm(ng)+3)=A(Lm(ng)  ,3)
#  endif
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          A(Lm(ng)+1,-2)=A(1,Mm(ng)-2)
          A(Lm(ng)+1,-1)=A(1,Mm(ng)-1)
          A(Lm(ng)+1, 0)=A(1,Mm(ng)  )
          A(Lm(ng)+2,-2)=A(2,Mm(ng)-2)
          A(Lm(ng)+2,-1)=A(2,Mm(ng)-1)
          A(Lm(ng)+2, 0)=A(2,Mm(ng)  )
#  ifdef THREE_GHOST
          A(Lm(ng)+3,-2)=A(3,Mm(ng)-2)
          A(Lm(ng)+3,-1)=A(3,Mm(ng)-1)
          A(Lm(ng)+3, 0)=A(3,Mm(ng)  )
#  endif
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          A(-2,-2)=A(Lm(ng)-2,Mm(ng)-2)
          A(-2,-1)=A(Lm(ng)-2,Mm(ng)-1)
          A(-2, 0)=A(Lm(ng)-2,Mm(ng)  )
          A(-1,-2)=A(Lm(ng)-1,Mm(ng)-2)
          A(-1,-1)=A(Lm(ng)-1,Mm(ng)-1)
          A(-1, 0)=A(Lm(ng)-1,Mm(ng)  )
          A( 0,-2)=A(Lm(ng)  ,Mm(ng)-2)
          A( 0,-1)=A(Lm(ng)  ,Mm(ng)-1)
          A( 0, 0)=A(Lm(ng)  ,Mm(ng)  )
        END IF
      END IF
# endif
      RETURN
      END SUBROUTINE exchange_v2d_tile
#endif
      END MODULE exchange_2d_mod
