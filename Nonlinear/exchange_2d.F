#include "cppdefs.h"
      MODULE exchange_2d_mod
#if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
!
!========================================== Alexander F. Shchepetkin ===
!  Copyright (c) 2005 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This package contains periodic boundary conditions and parallel     !
!  exchange (distributed-memory only) routines for 2D variables.       !
!                                                                      !
!  Routines:                                                           !
!                                                                      !
!    exchange_p2d_tile    periodic conditions/exchange at PSI-points   !
!    exchange_r2d_tile    periodic conditions/exchange at RHO-points   !
!    exchange_u2d_tile    periodic conditions/exchange at U-points     !
!    exchange_v2d_tile    periodic conditions/exchange at V-points     !
!                                                                      !
!=======================================================================
!
      implicit none

      CONTAINS
!
!***********************************************************************
      SUBROUTINE exchange_p2d_tile (ng, model, Istr, Iend, Jstr, Jend,  &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              Nghost,                             &
     &                              A)
!***********************************************************************
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Nghost

# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: A(LBi:,LBj:)
# else
      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
#  ifdef DISTRIBUTE
#   ifdef EW_PERIODIC
      logical :: EWperiodic=.TRUE.
#   else
      logical :: EWperiodic=.FALSE.
#   endif
#   ifdef NS_PERIODIC
      logical :: NSperiodic=.TRUE.
#   else
      logical :: NSperiodic=.FALSE.
#   endif
#  endif
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j

# include "set_bounds.h"

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
      IF (NtileI(ng).eq.1) THEN
#  endif
        IF (WESTERN_EDGE) THEN
          DO j=J_RANGE
            A(Lm(ng)+1,j)=A(1,j)
            A(Lm(ng)+2,j)=A(2,j)
#  ifdef THREE_GHOST
            A(Lm(ng)+3,j)=A(3,j)
#  endif
          END DO
        END IF
        IF (EASTERN_EDGE) THEN
          DO j=J_RANGE
            A(-2,j)=A(Lm(ng)-2,j)
            A(-1,j)=A(Lm(ng)-1,j)
            A( 0,j)=A(Lm(ng)  ,j)
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
      IF (NtileJ(ng).eq.1) THEN
#  endif
        IF (SOUTHERN_EDGE) THEN
          DO i=I_RANGE
            A(i,Mm(ng)+1)=A(i,1)
            A(i,Mm(ng)+2)=A(i,2)
#  ifdef THREE_GHOST
            A(i,Mm(ng)+3)=A(i,3)
#  endif
          END DO
        END IF
        IF (NORTHERN_EDGE) THEN
          DO i=I_RANGE
            A(i,-2)=A(i,Mm(ng)-2)
            A(i,-1)=A(i,Mm(ng)-1)
            A(i, 0)=A(i,Mm(ng)  )
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
        IF ((WESTERN_EDGE).and.(SOUTHERN_EDGE)) THEN
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
        IF ((EASTERN_EDGE).and.(SOUTHERN_EDGE)) THEN
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
        IF ((WESTERN_EDGE).and.(NORTHERN_EDGE)) THEN
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
        IF ((EASTERN_EDGE).and.(NORTHERN_EDGE)) THEN
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
      CALL mp_exchange (ng, model, Istr, Iend, Jstr, Jend,              &
     &                  LBi, UBi, LBj, UBj, 1, 1,                       &
     &                  Nghost, EWperiodic, NSperiodic,                 &
     &                  A)
# endif
      RETURN
      END SUBROUTINE exchange_p2d_tile

!
!***********************************************************************
      SUBROUTINE exchange_r2d_tile (ng, model, Istr, Iend, Jstr, Jend,  &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              Nghost,                             &
     &                              A)
!***********************************************************************
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Nghost

# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: A(LBi:,LBj:)
# else
      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
#  ifdef DISTRIBUTE
#   ifdef EW_PERIODIC
      logical :: EWperiodic=.TRUE.
#   else
      logical :: EWperiodic=.FALSE.
#   endif
#   ifdef NS_PERIODIC
      logical :: NSperiodic=.TRUE.
#   else
      logical :: NSperiodic=.FALSE.
#   endif
#  endif
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j

# include "set_bounds.h"

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
      IF (NtileI(ng).eq.1) THEN
#  endif
        IF (WESTERN_EDGE) THEN
          DO j=J_RANGE
            A(Lm(ng)+1,j)=A(1,j)
            A(Lm(ng)+2,j)=A(2,j)
#  ifdef THREE_GHOST
            A(Lm(ng)+3,j)=A(3,j)
#  endif
          END DO
        END IF
        IF (EASTERN_EDGE) THEN
          DO j=J_RANGE
            A(-2,j)=A(Lm(ng)-2,j)
            A(-1,j)=A(Lm(ng)-1,j)
            A( 0,j)=A(Lm(ng)  ,j)
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
      IF (NtileJ(ng).eq.1) THEN
#  endif
        IF (SOUTHERN_EDGE) THEN
          DO i=I_RANGE
            A(i,Mm(ng)+1)=A(i,1)
            A(i,Mm(ng)+2)=A(i,2)
#  ifdef THREE_GHOST
            A(i,Mm(ng)+3)=A(i,3)
#  endif
          END DO
        END IF
        IF (NORTHERN_EDGE) THEN
          DO i=I_RANGE
            A(i,-2)=A(i,Mm(ng)-2)
            A(i,-1)=A(i,Mm(ng)-1)
            A(i, 0)=A(i,Mm(ng)  )
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
        IF ((WESTERN_EDGE).and.(SOUTHERN_EDGE)) THEN
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
        IF ((EASTERN_EDGE).and.(SOUTHERN_EDGE)) THEN
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
        IF ((WESTERN_EDGE).and.(NORTHERN_EDGE)) THEN
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
        IF ((EASTERN_EDGE).and.(NORTHERN_EDGE)) THEN
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
      CALL mp_exchange (ng, model, Istr, Iend, Jstr, Jend,              &
     &                  LBi, UBi, LBj, UBj, 1, 1,                       &
     &                  Nghost, EWperiodic, NSperiodic,                 &
     &                  A)
# endif
      RETURN
      END SUBROUTINE exchange_r2d_tile

!
!***********************************************************************
      SUBROUTINE exchange_u2d_tile (ng, model, Istr, Iend, Jstr, Jend,  &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              Nghost,                             &
     &                              A)
!***********************************************************************
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Nghost

# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: A(LBi:,LBj:)
# else
      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
#  ifdef DISTRIBUTE
#   ifdef EW_PERIODIC
      logical :: EWperiodic=.TRUE.
#   else
      logical :: EWperiodic=.FALSE.
#   endif
#   ifdef NS_PERIODIC
      logical :: NSperiodic=.TRUE.
#   else
      logical :: NSperiodic=.FALSE.
#   endif
#  endif
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j

# include "set_bounds.h"

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
      IF (NtileI(ng).eq.1) THEN
#  endif
        IF (WESTERN_EDGE) THEN
          DO j=J_RANGE
            A(Lm(ng)+1,j)=A(1,j)
            A(Lm(ng)+2,j)=A(2,j)
#  ifdef THREE_GHOST
            A(Lm(ng)+3,j)=A(3,j)
#  endif
          END DO
        END IF
        IF (EASTERN_EDGE) THEN
          DO j=J_RANGE
            A(-2,j)=A(Lm(ng)-2,j)
            A(-1,j)=A(Lm(ng)-1,j)
            A( 0,j)=A(Lm(ng)  ,j)
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
      IF (NtileJ(ng).eq.1) THEN
#  endif
        IF (SOUTHERN_EDGE) THEN
          DO i=I_RANGE
            A(i,Mm(ng)+1)=A(i,1)
            A(i,Mm(ng)+2)=A(i,2)
#  ifdef THREE_GHOST
            A(i,Mm(ng)+3)=A(i,3)
#  endif
          END DO
        END IF
        IF (NORTHERN_EDGE) THEN
          DO i=I_RANGE
            A(i,-2)=A(i,Mm(ng)-2)
            A(i,-1)=A(i,Mm(ng)-1)
            A(i, 0)=A(i,Mm(ng)  )
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
        IF ((WESTERN_EDGE).and.(SOUTHERN_EDGE)) THEN
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
        IF ((EASTERN_EDGE).and.(SOUTHERN_EDGE)) THEN
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
        IF ((WESTERN_EDGE).and.(NORTHERN_EDGE)) THEN
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
        IF ((EASTERN_EDGE).and.(NORTHERN_EDGE)) THEN
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
      CALL mp_exchange (ng, model, Istr, Iend, Jstr, Jend,              &
     &                  LBi, UBi, LBj, UBj, 1, 1,                       &
     &                  Nghost, EWperiodic, NSperiodic,                 &
     &                  A)
# endif
      RETURN
      END SUBROUTINE exchange_u2d_tile

!
!***********************************************************************
      SUBROUTINE exchange_v2d_tile (ng, model, Istr, Iend, Jstr, Jend,  &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              Nghost,                             &
     &                              A)
!***********************************************************************
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Nghost

# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: A(LBi:,LBj:)
# else
      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
#  ifdef DISTRIBUTE
#   ifdef EW_PERIODIC
      logical :: EWperiodic=.TRUE.
#   else
      logical :: EWperiodic=.FALSE.
#   endif
#   ifdef NS_PERIODIC
      logical :: NSperiodic=.TRUE.
#   else
      logical :: NSperiodic=.FALSE.
#   endif
#  endif
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j

# include "set_bounds.h"

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
      IF (NtileI(ng).eq.1) THEN
#  endif
        IF (WESTERN_EDGE) THEN
          DO j=J_RANGE
            A(Lm(ng)+1,j)=A(1,j)
            A(Lm(ng)+2,j)=A(2,j)
#  ifdef THREE_GHOST
            A(Lm(ng)+3,j)=A(3,j)
#  endif
          END DO
        END IF
        IF (EASTERN_EDGE) THEN
          DO j=J_RANGE
            A(-2,j)=A(Lm(ng)-2,j)
            A(-1,j)=A(Lm(ng)-1,j)
            A( 0,j)=A(Lm(ng)  ,j)
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
      IF (NtileJ(ng).eq.1) THEN
#  endif
        IF (SOUTHERN_EDGE) THEN
          DO i=I_RANGE
            A(i,Mm(ng)+1)=A(i,1)
            A(i,Mm(ng)+2)=A(i,2)
#  ifdef THREE_GHOST
            A(i,Mm(ng)+3)=A(i,3)
#  endif
          END DO
        END IF
        IF (NORTHERN_EDGE) THEN
          DO i=I_RANGE
            A(i,-2)=A(i,Mm(ng)-2)
            A(i,-1)=A(i,Mm(ng)-1)
            A(i, 0)=A(i,Mm(ng)  )
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
        IF ((WESTERN_EDGE).and.(SOUTHERN_EDGE)) THEN
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
        IF ((EASTERN_EDGE).and.(SOUTHERN_EDGE)) THEN
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
        IF ((WESTERN_EDGE).and.(NORTHERN_EDGE)) THEN
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
        IF ((EASTERN_EDGE).and.(NORTHERN_EDGE)) THEN
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
      CALL mp_exchange (ng, model, Istr, Iend, Jstr, Jend,              &
     &                  LBi, UBi, LBj, UBj, 1, 1,                       &
     &                  Nghost, EWperiodic, NSperiodic,                 &
     &                  A)
# endif
      RETURN
      END SUBROUTINE exchange_v2d_tile
#endif
      END MODULE exchange_2d_mod
