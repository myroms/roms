#include "cppdefs.h"
      MODULE exchange_2d_mod
!
!========================================== Alexander F. Shchepetkin ===
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This package contains periodic boundary conditions routines for     !
!  2D primitive variables.                                             !
!                                                                      !
!  Routines:                                                           !
!                                                                      !
!    exchange_p2d_tile    periodic conditions at PSI-points            !
!    exchange_r2d_tile    periodic conditions at RHO-points            !
!    exchange_u2d_tile    periodic conditions at U-points              !
!    exchange_v2d_tile    periodic conditions at V-points              !
!                                                                      !
!=======================================================================
!
      implicit none

      CONTAINS
!
!***********************************************************************
      SUBROUTINE exchange_p2d_tile (ng, Istr, Iend, Jstr, Jend,         &
     &                              LBi, UBi, LBj, UBj, A)
!***********************************************************************
!
      USE mod_param
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
#ifdef DISTRIBUTE
      logical :: EWperiodic, NSperiodic
#endif
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j

#include "set_bounds.h"
!
#ifdef DISTRIBUTE
      EWperiodic=.false.
      NSperiodic=.false.
#endif

#ifdef EW_PERIODIC
# ifdef NS_PERIODIC
#  define J_RANGE Jstr,Jend
# else
#  define J_RANGE Jstr,JendR
# endif
!
!-----------------------------------------------------------------------
!  East-West periodic boundary conditions.
!-----------------------------------------------------------------------
!
# ifdef DISTRIBUTE
      EWperiodic=.true.
      IF (NtileI(ng).eq.1) THEN
# endif
        IF (WESTERN_EDGE) THEN
          DO j=J_RANGE
            A(Lm(ng)+1,j)=A(1,j)
            A(Lm(ng)+2,j)=A(2,j)
          END DO
        END IF
        IF (EASTERN_EDGE) THEN
          DO j=J_RANGE
            A(-2,j)=A(Lm(ng)-2,j)
            A(-1,j)=A(Lm(ng)-1,j)
            A( 0,j)=A(Lm(ng)  ,j)
          END DO
        END IF
# ifdef DISTRIBUTE
      END IF      
# endif
# undef J_RANGE
#endif

#ifdef NS_PERIODIC
# ifdef EW_PERIODIC
#  define I_RANGE Istr,Iend
# else
#  define I_RANGE Istr,IendR
# endif
!
!-----------------------------------------------------------------------
!  North-South periodic boundary conditions.
!-----------------------------------------------------------------------
!
# ifdef DISTRIBUTE
      NSperiodic=.true.
      IF (NtileJ(ng).eq.1) THEN
# endif
        IF (SOUTHERN_EDGE) THEN
          DO i=I_RANGE
            A(i,Mm(ng)+1)=A(i,1)
            A(i,Mm(ng)+2)=A(i,2)
          END DO
        END IF
        IF (NORTHERN_EDGE) THEN
          DO i=I_RANGE
            A(i,-2)=A(i,Mm(ng)-2)
            A(i,-1)=A(i,Mm(ng)-1)
            A(i, 0)=A(i,Mm(ng)  )
          END DO
        END IF
# ifdef DISTRIBUTE
      END IF      
# endif
# undef I_RANGE
#endif

#if defined EW_PERIODIC && defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
# ifdef DISTRIBUTE
      IF ((NtileI(ng).eq.1).and.(NtileJ(ng).eq.1)) THEN
# endif
        IF (WESTERN_EDGE.and.SOUTHERN_EDGE) THEN
          A(Lm(ng)+1,Mm(ng)+1)=A(1,1)
          A(Lm(ng)+1,Mm(ng)+2)=A(1,2)
          A(Lm(ng)+2,Mm(ng)+1)=A(2,1)
          A(Lm(ng)+2,Mm(ng)+2)=A(2,2)
        END IF
        IF (EASTERN_EDGE.and.SOUTHERN_EDGE) THEN
          A(-2,Mm(ng)+1)=A(Lm(ng)-2,1)
          A(-1,Mm(ng)+1)=A(Lm(ng)-1,1)
          A( 0,Mm(ng)+1)=A(Lm(ng)  ,1)
          A(-2,Mm(ng)+2)=A(Lm(ng)-2,2)
          A(-1,Mm(ng)+2)=A(Lm(ng)-1,2)
          A( 0,Mm(ng)+2)=A(Lm(ng)  ,2)
        END IF
        IF (WESTERN_EDGE.and.NORTHERN_EDGE) THEN
          A(Lm(ng)+1,-2)=A(1,Mm(ng)-2)
          A(Lm(ng)+1,-1)=A(1,Mm(ng)-1)
          A(Lm(ng)+1, 0)=A(1,Mm(ng)  )
          A(Lm(ng)+2,-2)=A(2,Mm(ng)-2)
          A(Lm(ng)+2,-1)=A(2,Mm(ng)-1)
          A(Lm(ng)+2, 0)=A(2,Mm(ng)  )
        END IF
        IF (EASTERN_EDGE.and.NORTHERN_EDGE) THEN
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
# ifdef DISTRIBUTE
      END IF
# endif
#endif
#ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Send/recieve tile boundary data.
!-----------------------------------------------------------------------
!
      CALL mp_exchange (ng, Istr, Iend, Jstr, Jend,                     &
     &                  LBi, UBi, LBj, UBj, 1, 1,                       &
     &                  GHOST_POINTS, EWperiodic, NSperiodic, A)
#endif
      RETURN
      END SUBROUTINE exchange_p2d_tile

!
!***********************************************************************
      SUBROUTINE exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,         &
     &                              LBi, UBi, LBj, UBj, A)
!***********************************************************************
!
      USE mod_param
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
#ifdef DISTRIBUTE
      logical :: EWperiodic, NSperiodic
#endif
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j

#include "set_bounds.h"
!
#ifdef DISTRIBUTE
      EWperiodic=.false.
      NSperiodic=.false.
#endif

#ifdef EW_PERIODIC
# ifdef NS_PERIODIC
#  define J_RANGE Jstr,Jend
# else
#  define J_RANGE JstrR,JendR
# endif
!
!-----------------------------------------------------------------------
!  East-West periodic boundary conditions.
!-----------------------------------------------------------------------
!
# ifdef DISTRIBUTE
      EWperiodic=.true.
      IF (NtileI(ng).eq.1) THEN
# endif
        IF (WESTERN_EDGE) THEN
          DO j=J_RANGE
            A(Lm(ng)+1,j)=A(1,j)
            A(Lm(ng)+2,j)=A(2,j)
          END DO
        END IF
        IF (EASTERN_EDGE) THEN
          DO j=J_RANGE
            A(-2,j)=A(Lm(ng)-2,j)
            A(-1,j)=A(Lm(ng)-1,j)
            A( 0,j)=A(Lm(ng)  ,j)
          END DO
        END IF
# ifdef DISTRIBUTE
      END IF
# endif
# undef J_RANGE
#endif

#ifdef NS_PERIODIC
# ifdef EW_PERIODIC
#  define I_RANGE Istr,Iend
# else
#  define I_RANGE IstrR,IendR
# endif
!
!-----------------------------------------------------------------------
!  North-South periodic boundary conditions.
!-----------------------------------------------------------------------
!
# ifdef DISTRIBUTE
      NSperiodic=.true.
      IF (NtileJ(ng).eq.1) THEN
# endif
        IF (SOUTHERN_EDGE) THEN
          DO i=I_RANGE
            A(i,Mm(ng)+1)=A(i,1)
            A(i,Mm(ng)+2)=A(i,2)
          END DO
        END IF
        IF (NORTHERN_EDGE) THEN
          DO i=I_RANGE
            A(i,-2)=A(i,Mm(ng)-2)
            A(i,-1)=A(i,Mm(ng)-1)
            A(i, 0)=A(i,Mm(ng)  )
          END DO
        END IF
# ifdef DISTRIBUTE
      END IF      
# endif
# undef I_RANGE
#endif

#if defined EW_PERIODIC && defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
# ifdef DISTRIBUTE
      IF ((NtileI(ng).eq.1).and.(NtileJ(ng).eq.1)) THEN
# endif
        IF (WESTERN_EDGE.and.SOUTHERN_EDGE) THEN
          A(Lm(ng)+1,Mm(ng)+1)=A(1,1)
          A(Lm(ng)+1,Mm(ng)+2)=A(1,2)
          A(Lm(ng)+2,Mm(ng)+1)=A(2,1)
          A(Lm(ng)+2,Mm(ng)+2)=A(2,2)
        END IF
        IF (EASTERN_EDGE.and.SOUTHERN_EDGE) THEN
          A(-2,Mm(ng)+1)=A(Lm(ng)-2,1)
          A(-1,Mm(ng)+1)=A(Lm(ng)-1,1)
          A( 0,Mm(ng)+1)=A(Lm(ng)  ,1)
          A(-2,Mm(ng)+2)=A(Lm(ng)-2,2)
          A(-1,Mm(ng)+2)=A(Lm(ng)-1,2)
          A( 0,Mm(ng)+2)=A(Lm(ng)  ,2)
        END IF
        IF (WESTERN_EDGE.and.NORTHERN_EDGE) THEN
          A(Lm(ng)+1,-2)=A(1,Mm(ng)-2)
          A(Lm(ng)+1,-1)=A(1,Mm(ng)-1)
          A(Lm(ng)+1, 0)=A(1,Mm(ng)  )
          A(Lm(ng)+2,-2)=A(2,Mm(ng)-2)
          A(Lm(ng)+2,-1)=A(2,Mm(ng)-1)
          A(Lm(ng)+2, 0)=A(2,Mm(ng)  )
        END IF
        IF (EASTERN_EDGE.and.NORTHERN_EDGE) THEN
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
# ifdef DISTRIBUTE
      END IF      
# endif
#endif
#ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Send/recieve tile boundary data.
!-----------------------------------------------------------------------
!
      CALL mp_exchange (ng, Istr, Iend, Jstr, Jend,                     &
     &                  LBi, UBi, LBj, UBj, 1, 1,                       &
     &                  GHOST_POINTS, EWperiodic, NSperiodic, A)
#endif
      RETURN
      END SUBROUTINE exchange_r2d_tile

!
!***********************************************************************
      SUBROUTINE exchange_u2d_tile (ng, Istr, Iend, Jstr, Jend,         &
     &                              LBi, UBi, LBj, UBj, A)
!***********************************************************************
!
      USE mod_param
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
#ifdef DISTRIBUTE
      logical :: EWperiodic, NSperiodic
#endif
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j

#include "set_bounds.h"
!
#ifdef DISTRIBUTE
      EWperiodic=.false.
      NSperiodic=.false.
#endif

#ifdef EW_PERIODIC
# ifdef NS_PERIODIC
#  define J_RANGE Jstr,Jend
# else
#  define J_RANGE JstrR,JendR
# endif
!
!-----------------------------------------------------------------------
!  East-West periodic boundary conditions.
!-----------------------------------------------------------------------
!
# ifdef DISTRIBUTE
      EWperiodic=.true.
      IF (NtileI(ng).eq.1) THEN
# endif
        IF (WESTERN_EDGE) THEN
          DO j=J_RANGE
            A(Lm(ng)+1,j)=A(1,j)
            A(Lm(ng)+2,j)=A(2,j)
          END DO
        END IF
        IF (EASTERN_EDGE) THEN
          DO j=J_RANGE
            A(-2,j)=A(Lm(ng)-2,j)
            A(-1,j)=A(Lm(ng)-1,j)
            A( 0,j)=A(Lm(ng)  ,j)
          END DO
        END IF
# ifdef DISTRIBUTE
      END IF      
# endif
# undef J_RANGE
#endif

#ifdef NS_PERIODIC
# ifdef EW_PERIODIC
#  define I_RANGE Istr,Iend
# else
#  define I_RANGE Istr,IendR
# endif
!
!-----------------------------------------------------------------------
!  North-South periodic boundary conditions.
!-----------------------------------------------------------------------
!
# ifdef DISTRIBUTE
      NSperiodic=.true.
      IF (NtileJ(ng).eq.1) THEN
# endif
        IF (SOUTHERN_EDGE) THEN
          DO i=I_RANGE
            A(i,Mm(ng)+1)=A(i,1)
            A(i,Mm(ng)+2)=A(i,2)
          END DO
        END IF
        IF (NORTHERN_EDGE) THEN
          DO i=I_RANGE
            A(i,-2)=A(i,Mm(ng)-2)
            A(i,-1)=A(i,Mm(ng)-1)
            A(i, 0)=A(i,Mm(ng)  )
          END DO
        END IF
# ifdef DISTRIBUTE
      END IF      
# endif
# undef I_RANGE
#endif

#if defined EW_PERIODIC && defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
# ifdef DISTRIBUTE
      IF ((NtileI(ng).eq.1).and.(NtileJ(ng).eq.1)) THEN
# endif
        IF (WESTERN_EDGE.and.SOUTHERN_EDGE) THEN
          A(Lm(ng)+1,Mm(ng)+1)=A(1,1)
          A(Lm(ng)+1,Mm(ng)+2)=A(1,2)
          A(Lm(ng)+2,Mm(ng)+1)=A(2,1)
          A(Lm(ng)+2,Mm(ng)+2)=A(2,2)
        END IF
        IF (EASTERN_EDGE.and.SOUTHERN_EDGE) THEN
          A(-2,Mm(ng)+1)=A(Lm(ng)-2,1)
          A(-1,Mm(ng)+1)=A(Lm(ng)-1,1)
          A( 0,Mm(ng)+1)=A(Lm(ng)  ,1)
          A(-2,Mm(ng)+2)=A(Lm(ng)-2,2)
          A(-1,Mm(ng)+2)=A(Lm(ng)-1,2)
          A( 0,Mm(ng)+2)=A(Lm(ng)  ,2)
        END IF
        IF (WESTERN_EDGE.and.NORTHERN_EDGE) THEN
          A(Lm(ng)+1,-2)=A(1,Mm(ng)-2)
          A(Lm(ng)+1,-1)=A(1,Mm(ng)-1)
          A(Lm(ng)+1, 0)=A(1,Mm(ng)  )
          A(Lm(ng)+2,-2)=A(2,Mm(ng)-2)
          A(Lm(ng)+2,-1)=A(2,Mm(ng)-1)
          A(Lm(ng)+2, 0)=A(2,Mm(ng)  )
        END IF
        IF (EASTERN_EDGE.and.NORTHERN_EDGE) THEN
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
# ifdef DISTRIBUTE
      END IF      
# endif
#endif
#ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Send/recieve tile boundary data.
!-----------------------------------------------------------------------
!
      CALL mp_exchange (ng, Istr, Iend, Jstr, Jend,                     &
     &                  LBi, UBi, LBj, UBj, 1, 1,                       &
     &                  GHOST_POINTS, EWperiodic, NSperiodic, A)
#endif
      RETURN
      END SUBROUTINE exchange_u2d_tile

!
!***********************************************************************
      SUBROUTINE exchange_v2d_tile (ng, Istr, Iend, Jstr, Jend,         &
     &                              LBi, UBi, LBj, UBj, A)
!***********************************************************************
!
      USE mod_param
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
#ifdef DISTRIBUTE
      logical :: EWperiodic, NSperiodic
#endif
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j

#include "set_bounds.h"
!
#ifdef DISTRIBUTE
      EWperiodic=.false.
      NSperiodic=.false.
#endif

#ifdef EW_PERIODIC
# ifdef NS_PERIODIC
#  define J_RANGE Jstr,Jend
# else
#  define J_RANGE Jstr,JendR
# endif
!
!-----------------------------------------------------------------------
!  East-West periodic boundary conditions.
!-----------------------------------------------------------------------
!
# ifdef DISTRIBUTE
      EWperiodic=.true.
      IF (NtileI(ng).eq.1) THEN
# endif
        IF (WESTERN_EDGE) THEN
          DO j=J_RANGE
            A(Lm(ng)+1,j)=A(1,j)
            A(Lm(ng)+2,j)=A(2,j)
          END DO
        END IF
        IF (EASTERN_EDGE) THEN
          DO j=J_RANGE
            A(-2,j)=A(Lm(ng)-2,j)
            A(-1,j)=A(Lm(ng)-1,j)
            A( 0,j)=A(Lm(ng)  ,j)
          END DO
        END IF
# ifdef DISTRIBUTE
      END IF      
# endif
# undef J_RANGE
#endif

#ifdef NS_PERIODIC
# ifdef EW_PERIODIC
#  define I_RANGE Istr,Iend
# else
#  define I_RANGE IstrR,IendR
# endif
!
!-----------------------------------------------------------------------
!  North-South periodic boundary conditions.
!-----------------------------------------------------------------------
!
# ifdef DISTRIBUTE
      NSperiodic=.true.
      IF (NtileJ(ng).eq.1) THEN
# endif
        IF (SOUTHERN_EDGE) THEN
          DO i=I_RANGE
            A(i,Mm(ng)+1)=A(i,1)
            A(i,Mm(ng)+2)=A(i,2)
          END DO
        END IF
        IF (NORTHERN_EDGE) THEN
          DO i=I_RANGE
            A(i,-2)=A(i,Mm(ng)-2)
            A(i,-1)=A(i,Mm(ng)-1)
            A(i, 0)=A(i,Mm(ng)  )
          END DO
        END IF
# ifdef DISTRIBUTE
      END IF      
# endif
# undef I_RANGE
#endif

#if defined EW_PERIODIC && defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
# ifdef DISTRIBUTE
      IF ((NtileI(ng).eq.1).and.(NtileJ(ng).eq.1)) THEN
# endif
        IF (WESTERN_EDGE.and.SOUTHERN_EDGE) THEN
          A(Lm(ng)+1,Mm(ng)+1)=A(1,1)
          A(Lm(ng)+1,Mm(ng)+2)=A(1,2)
          A(Lm(ng)+2,Mm(ng)+1)=A(2,1)
          A(Lm(ng)+2,Mm(ng)+2)=A(2,2)
        END IF
        IF (EASTERN_EDGE.and.SOUTHERN_EDGE) THEN
          A(-2,Mm(ng)+1)=A(Lm(ng)-2,1)
          A(-1,Mm(ng)+1)=A(Lm(ng)-1,1)
          A( 0,Mm(ng)+1)=A(Lm(ng)  ,1)
          A(-2,Mm(ng)+2)=A(Lm(ng)-2,2)
          A(-1,Mm(ng)+2)=A(Lm(ng)-1,2)
          A( 0,Mm(ng)+2)=A(Lm(ng)  ,2)
        END IF
        IF (WESTERN_EDGE.and.NORTHERN_EDGE) THEN
          A(Lm(ng)+1,-2)=A(1,Mm(ng)-2)
          A(Lm(ng)+1,-1)=A(1,Mm(ng)-1)
          A(Lm(ng)+1, 0)=A(1,Mm(ng)  )
          A(Lm(ng)+2,-2)=A(2,Mm(ng)-2)
          A(Lm(ng)+2,-1)=A(2,Mm(ng)-1)
          A(Lm(ng)+2, 0)=A(2,Mm(ng)  )
        END IF
        IF (EASTERN_EDGE.and.NORTHERN_EDGE) THEN
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
# ifdef DISTRIBUTE
      END IF      
# endif
#endif
#ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Send/recieve tile boundary data.
!-----------------------------------------------------------------------
!
      CALL mp_exchange (ng, Istr, Iend, Jstr, Jend,                     &
     &                  LBi, UBi, LBj, UBj, 1, 1,                       &
     &                  GHOST_POINTS, EWperiodic, NSperiodic, A)
#endif
      RETURN
      END SUBROUTINE exchange_v2d_tile

      END MODULE exchange_2d_mod
