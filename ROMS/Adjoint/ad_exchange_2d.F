#include "cppdefs.h"
      MODULE ad_exchange_2d_mod
#if defined ADJOINT && (defined EW_PERIODIC || defined NS_PERIODIC)
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2011 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  These routines apply periodic boundary conditions to generic        !
!  2D adjoint fields.                                                  !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng                      Nested grid number.                      !
!     tile                    Domain partition.                        !
!     LBi                     I-dimension Lower bound.                 !
!     UBi                     I-dimension Upper bound.                 !
!     LBj                     J-dimension Lower bound.                 !
!     UBj                     J-dimension Upper bound.                 !
!     ad_A                    2D adjoint field.                        !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     ad_A                    Processed 2D adjoint field               !
!                                                                      !
!  Routines:                                                           !
!                                                                      !
!     ad_exchange_p2d_tile    periodic conditions at PSI-points        !
!     ad_exchange_r2d_tile    periodic conditions at RHO-points        !
!     ad_exchange_u2d_tile    periodic conditions at U-points          !
!     ad_exchange_v2d_tile    periodic conditions at V-points          !
!                                                                      !
!=======================================================================
!
      implicit none

      CONTAINS
!
!***********************************************************************
      SUBROUTINE ad_exchange_p2d_tile (ng, tile,                        &
     &                                 LBi, UBi, LBj, UBj,              &
     &                                 ad_A)
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
      real(r8), intent(inout) :: ad_A(LBi:,LBj:)
# else
      real(r8), intent(inout) :: ad_A(LBi:UBi,LBj:UBj)
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

# if defined EW_PERIODIC && defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (EW_exchange.and.NS_exchange) THEN
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
!>        tl_A(-2,-2)=tl_A(Lm(ng)-2,Mm(ng)-2)
!>
          ad_A(Lm(ng)-2,Mm(ng)-2)=ad_A(Lm(ng)-2,Mm(ng)-2)+              &
     &                            ad_A(-2,-2)
          ad_A(-2,-2)=0.0_r8
!>        tl_A(-2,-1)=tl_A(Lm(ng)-2,Mm(ng)-1)
!>
          ad_A(Lm(ng)-2,Mm(ng)-1)=ad_A(Lm(ng)-2,Mm(ng)-1)+              &
     &                            ad_A(-2,-1)
          ad_A(-2,-1)=0.0_r8
!>        tl_A(-2, 0)=tl_A(Lm(ng)-2,Mm(ng)  )
!>
          ad_A(Lm(ng)-2,Mm(ng)  )=ad_A(Lm(ng)-2,Mm(ng)  )+              &
     &                            ad_A(-2, 0)
          ad_A(-2, 0)=0.0_r8
!>        tl_A(-1,-2)=tl_A(Lm(ng)-1,Mm(ng)-2)
!>
          ad_A(Lm(ng)-1,Mm(ng)-2)=ad_A(Lm(ng)-1,Mm(ng)-2)+              &
     &                            ad_A(-1,-2)
          ad_A(-1,-2)=0.0_r8
!>        tl_A(-1,-1)=tl_A(Lm(ng)-1,Mm(ng)-1)
!>
          ad_A(Lm(ng)-1,Mm(ng)-1)=ad_A(Lm(ng)-1,Mm(ng)-1)+              &
     &                            ad_A(-1,-1)
          ad_A(-1,-1)=0.0_r8
!>        tl_A(-1, 0)=tl_A(Lm(ng)-1,Mm(ng)  )
!>
          ad_A(Lm(ng)-1,Mm(ng)  )=ad_A(Lm(ng)-1,Mm(ng)  )+              &
     &                            ad_A(-1, 0)
          ad_A(-1, 0)=0.0_r8
!>        tl_A( 0,-2)=tl_A(Lm(ng)  ,Mm(ng)-2)
!>
          ad_A(Lm(ng)  ,Mm(ng)-2)=ad_A(Lm(ng)  ,Mm(ng)-2)+              &
     &                            ad_A( 0,-2)
          ad_A( 0,-2)=0.0_r8
!>        tl_A( 0,-1)=tl_A(Lm(ng)  ,Mm(ng)-1)
!>
          ad_A(Lm(ng)  ,Mm(ng)-1)=ad_A(Lm(ng)  ,Mm(ng)-1)+              &
     &                            ad_A( 0,-1)
          ad_A( 0,-1)=0.0_r8
!>        tl_A( 0, 0)=tl_A(Lm(ng)  ,Mm(ng)  )
!>
          ad_A(Lm(ng)  ,Mm(ng)  )=ad_A(Lm(ng)  ,Mm(ng)  )+              &
     &                            ad_A( 0, 0)
          ad_A( 0, 0)=0.0_r8
        END IF

        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
!>        tl_A(Lm(ng)+1,-2)=tl_A(1 ,Mm(ng)-2)
!>
          ad_A(1 ,Mm(ng)-2)=ad_A(1 ,Mm(ng)-2)+                          &
     &                      ad_A(Lm(ng)+1,-2)
          ad_A(Lm(ng)+1,-2)=0.0_r8
!>        tl_A(Lm(ng)+1,-1)=tl_A(1 ,Mm(ng)-1)
!>
          ad_A(1 ,Mm(ng)-1)=ad_A(1 ,Mm(ng)-1)+                          &
     &                      ad_A(Lm(ng)+1,-1)
          ad_A(Lm(ng)+1,-1)=0.0_r8
!>        tl_A(Lm(ng)+1, 0)=tl_A(1 ,Mm(ng)  )
!>
          ad_A(1 ,Mm(ng)  )=ad_A(1 ,Mm(ng)  )+                          &
     &                      ad_A(Lm(ng)+1, 0)
          ad_A(Lm(ng)+1, 0)=0.0_r8
!>        tl_A(Lm(ng)+2,-2)=tl_A(2 ,Mm(ng)-2)
!>
          ad_A(2 ,Mm(ng)-2)=ad_A(2 ,Mm(ng)-2)+                          &
     &                      ad_A(Lm(ng)+2,-2)
          ad_A(Lm(ng)+2,-2)=0.0_r8
!>        tl_A(Lm(ng)+2,-1)=tl_A(2 ,Mm(ng)-1)
!>
          ad_A(2 ,Mm(ng)-1)=ad_A(2 ,Mm(ng)-1)+                          &
     &                      ad_A(Lm(ng)+2,-1)
          ad_A(Lm(ng)+2,-1)=0.0_r8
!>        tl_A(Lm(ng)+2, 0)=tl_A(2 ,Mm(ng)  )
!>
          ad_A(2 ,Mm(ng)  )=ad_A(2 ,Mm(ng)  )+                          &
     &                      ad_A(Lm(ng)+2, 0)
          ad_A(Lm(ng)+2, 0)=0.0_r8
#  ifdef THREE_GHOST
!>        tl_A(Lm(ng)+3,-2)=tl_A(3 ,Mm(ng)-2)
!>
          ad_A(3 ,Mm(ng)-2)=ad_A(3 ,Mm(ng)-2)+                          &
     &                      ad_A(Lm(ng)+3,-2)
          ad_A(Lm(ng)+3,-2)=0.0_r8
!>        tl_A(Lm(ng)+3,-1)=tl_A(3 ,Mm(ng)-1)
!>
          ad_A(3 ,Mm(ng)-1)=ad_A(3 ,Mm(ng)-1)+                          &
     &                      ad_A(Lm(ng)+3,-1)
          ad_A(Lm(ng)+3,-1)=0.0_r8
!>        tl_A(Lm(ng)+3, 0)=tl_A(3 ,Mm(ng)  )
!>
          ad_A(3 ,Mm(ng)  )=ad_A(3 ,Mm(ng)  )+                          &
     &                      ad_A(Lm(ng)+3, 0)
          ad_A(Lm(ng)+3, 0)=0.0_r8
#  endif
        END IF

        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
!>        tl_A(-2,Mm(ng)+1)=tl_A(Lm(ng)-2, 1)
!>
          ad_A(Lm(ng)-2, 1)=ad_A(Lm(ng)-2, 1)+                          &
     &                      ad_A(-2,Mm(ng)+1)
          ad_A(-2,Mm(ng)+1)=0.0_r8
!>        tl_A(-1,Mm(ng)+1)=tl_A(Lm(ng)-1, 1)
!>
          ad_A(Lm(ng)-1, 1)=ad_A(Lm(ng)-1, 1)+                          &
     &                      ad_A(-1,Mm(ng)+1)
          ad_A(-1,Mm(ng)+1)=0.0_r8
!>        tl_A( 0,Mm(ng)+1)=tl_A(Lm(ng)  , 1)
!>
          ad_A(Lm(ng)  , 1)=ad_A(Lm(ng)  , 1)+                          &
     &                      ad_A( 0,Mm(ng)+1)
          ad_A( 0,Mm(ng)+1)=0.0_r8
!>        tl_A(-2,Mm(ng)+2)=tl_A(Lm(ng)-2, 2)
!>
          ad_A(Lm(ng)-2, 2)=ad_A(Lm(ng)-2, 2)+                          &
     &                      ad_A(-2,Mm(ng)+2)
          ad_A(-2,Mm(ng)+2)=0.0_r8
!>        tl_A(-1,Mm(ng)+2)=tl_A(Lm(ng)-1, 2)
!>
          ad_A(Lm(ng)-1, 2)=ad_A(Lm(ng)-1, 2)+                          &
     &                      ad_A(-1,Mm(ng)+2)
          ad_A(-1,Mm(ng)+2)=0.0_r8
!>        tl_A( 0,Mm(ng)+2)=tl_A(Lm(ng)  , 2)
!>
          ad_A(Lm(ng)  , 2)=ad_A(Lm(ng)  , 2)+                          &
     &                      ad_A( 0,Mm(ng)+2)
          ad_A( 0,Mm(ng)+2)=0.0_r8
#  ifdef THREE_GHOST
!>        tl_A(-2,Mm(ng)+3)=tl_A(Lm(ng)-2, 3)
!>
          ad_A(Lm(ng)-2, 3)=ad_A(Lm(ng)-2, 3)+                          &
     &                      ad_A(-2,Mm(ng)+3)
          ad_A(-2,Mm(ng)+3)=0.0_r8
!>        tl_A(-1,Mm(ng)+3)=tl_A(Lm(ng)-1, 3)
!>
          ad_A(Lm(ng)-1, 3)=ad_A(Lm(ng)-1, 3)+                          &
     &                      ad_A(-1,Mm(ng)+3)
          ad_A(-1,Mm(ng)+3)=0.0_r8
!>        tl_A( 0,Mm(ng)+3)=tl_A(Lm(ng)  , 3)
!>
          ad_A(Lm(ng)  , 3)=ad_A(Lm(ng)  , 3)+                          &
     &                      ad_A( 0,Mm(ng)+3)
          ad_A( 0,Mm(ng)+3)=0.0_r8
#  endif
        END IF

        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
!>        tl_A(Lm(ng)+1,Mm(ng)+1)=tl_A( 1, 1)
!>
          ad_A( 1, 1)=ad_A( 1, 1)+                                      &
     &                ad_A(Lm(ng)+1,Mm(ng)+1)
          ad_A(Lm(ng)+1,Mm(ng)+1)=0.0_r8
!>        tl_A(Lm(ng)+1,Mm(ng)+2)=tl_A( 1, 2)
!>
          ad_A( 1, 2)=ad_A( 1, 2)+                                      &
     &                ad_A(Lm(ng)+1,Mm(ng)+2)
          ad_A(Lm(ng)+1,Mm(ng)+2)=0.0_r8
#  ifdef THREE_GHOST
!>        tl_A(Lm(ng)+1,Mm(ng)+3)=tl_A( 1, 3)
!>
          ad_A( 1, 3)=ad_A( 1, 3)+                                      &
     &                ad_A(Lm(ng)+1,Mm(ng)+3)
          ad_A(Lm(ng)+1,Mm(ng)+3)=0.0_r8
#  endif
!>        tl_A(Lm(ng)+2,Mm(ng)+1)=tl_A( 2, 1)
!>
          ad_A( 2, 1)=ad_A( 2, 1)+                                      &
     &                ad_A(Lm(ng)+2,Mm(ng)+1)
          ad_A(Lm(ng)+2,Mm(ng)+1)=0.0_r8
!>        tl_A(Lm(ng)+2,Mm(ng)+2)=tl_A( 2, 2)
!>
          ad_A( 2, 2)=ad_A( 2, 2)+                                      &
     &                ad_A(Lm(ng)+2,Mm(ng)+2)
          ad_A(Lm(ng)+2,Mm(ng)+2)=0.0_r8
#  ifdef THREE_GHOST
!>        tl_A(Lm(ng)+2,Mm(ng)+3)=tl_A( 2, 3)
!>
          ad_A( 2, 3)=ad_A( 2, 3)+                                      &
     &                ad_A(Lm(ng)+2,Mm(ng)+3)
          ad_A(Lm(ng)+2,Mm(ng)+3)=0.0_r8
!>        tl_A(Lm(ng)+3,Mm(ng)+1)=tl_A( 3, 1)
!>
          ad_A( 3, 1)=ad_A( 3, 1)+                                      &
     &                ad_A(Lm(ng)+3,Mm(ng)+1)
          ad_A(Lm(ng)+3,Mm(ng)+1)=0.0_r8
!>        tl_A(Lm(ng)+3,Mm(ng)+2)=tl_A( 3, 2)
!>
          ad_A( 3, 2)=ad_A( 3, 2)+                                      &
     &                ad_A(Lm(ng)+3,Mm(ng)+2)
          ad_A(Lm(ng)+3,Mm(ng)+2)=0.0_r8
!>        tl_A(Lm(ng)+3,Mm(ng)+3)=tl_A( 3, 3)
!>
          ad_A( 3, 3)=ad_A( 3, 3)+                                      &
     &                ad_A(Lm(ng)+3,Mm(ng)+3)
          ad_A(Lm(ng)+3,Mm(ng)+3)=0.0_r8
#  endif
        END IF
      END IF
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
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          DO i=I_RANGE
!>          tl_A(i,-2)=tl_A(i,Mm(ng)-2)
!>
            ad_A(i,Mm(ng)-2)=ad_A(i,Mm(ng)-2)+                          &
     &                       ad_A(i,-2)
            ad_A(i,-2)=0.0_r8
!>          tl_A(i,-1)=tl_A(i,Mm(ng)-1)
!>
            ad_A(i,Mm(ng)-1)=ad_A(i,Mm(ng)-1)+                          &
     &                       ad_A(i,-1)
            ad_A(i,-1)=0.0_r8
!>          tl_A(i, 0)=tl_A(i,Mm(ng)  )
!>
            ad_A(i,Mm(ng)  )=ad_A(i,Mm(ng)  )+                          &
     &                       ad_A(i, 0)
            ad_A(i, 0)=0.0_r8
          END DO
        END IF
        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO i=I_RANGE
!>          tl_A(i,Mm(ng)+1)=tl_A(i,1)
!>
            ad_A(i,1)=ad_A(i,1)+                                        &
     &                ad_A(i,Mm(ng)+1)
            ad_A(i,Mm(ng)+1)=0.0_r8
!>          tl_A(i,Mm(ng)+2)=tl_A(i,2)
!>
            ad_A(i,2)=ad_A(i,2)+                                        &
     &                ad_A(i,Mm(ng)+2)
            ad_A(i,Mm(ng)+2)=0.0_r8
#  ifdef THREE_GHOST
!>          tl_A(i,Mm(ng)+3)=tl_A(i,3)
!>
            ad_A(i,3)=ad_A(i,3)+                                        &
     &                ad_A(i,Mm(ng)+3)
            ad_A(i,Mm(ng)+3)=0.0_r8
#  endif
          END DO
        END IF
      END IF
#  undef I_RANGE
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
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO j=J_RANGE
!>          tl_A(-2,j)=tl_A(Lm(ng)-2,j)
!>
            ad_A(Lm(ng)-2,j)=ad_A(Lm(ng)-2,j)+                          &
     &                       ad_A(-2,j)
            ad_A(-2,j)=0.0_r8
!>          tl_A(-1,j)=tl_A(Lm(ng)-1,j)
!>
            ad_A(Lm(ng)-1,j)=ad_A(Lm(ng)-1,j)+                          &
     &                       ad_A(-1,j)
            ad_A(-1,j)=0.0_r8
!>          tl_A( 0,j)=tl_A(Lm(ng)  ,j)
!>
            ad_A(Lm(ng)  ,j)=ad_A(Lm(ng)  ,j)+                          &
     &                       ad_A( 0,j)
            ad_A( 0,j)=0.0_r8
          END DO
        END IF

        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          DO j=J_RANGE
!>          tl_A(Lm(ng)+1,j)=tl_A(1,j)
!>
            ad_A(1,j)=ad_A(1,j)+                                        &
     &                ad_A(Lm(ng)+1,j)
            ad_A(Lm(ng)+1,j)=0.0_r8
!>          tl_A(Lm(ng)+2,j)=tl_A(2,j)
!>
            ad_A(2,j)=ad_A(2,j)+                                        &
     &                ad_A(Lm(ng)+2,j)
            ad_A(Lm(ng)+2,j)=0.0_r8
#  ifdef THREE_GHOST
!>          tl_A(Lm(ng)+3,j)=tl_A(3,j)
!>
            ad_A(3,j)=ad_A(3,j)+                                        &
     &                ad_A(Lm(ng)+3,j)
            ad_A(Lm(ng)+3,j)=0.0_r8
#  endif
          END DO
        ENDIF
      END IF
#  undef J_RANGE
# endif

      RETURN
      END SUBROUTINE ad_exchange_p2d_tile

!
!***********************************************************************
      SUBROUTINE ad_exchange_r2d_tile (ng, tile,                        &
     &                                 LBi, UBi, LBj, UBj,              &
     &                                 ad_A)
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
      real(r8), intent(inout) :: ad_A(LBi:,LBj:)
# else
      real(r8), intent(inout) :: ad_A(LBi:UBi,LBj:UBj)
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

# if defined EW_PERIODIC && defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (EW_exchange.and.NS_exchange) THEN
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
!>        tl_A(-2,-2)=tl_A(Lm(ng)-2,Mm(ng)-2)
!>
          ad_A(Lm(ng)-2,Mm(ng)-2)=ad_A(Lm(ng)-2,Mm(ng)-2)+              &
     &                            ad_A(-2,-2)
          ad_A(-2,-2)=0.0_r8
!>        tl_A(-2,-1)=tl_A(Lm(ng)-2,Mm(ng)-1)
!>
          ad_A(Lm(ng)-2,Mm(ng)-1)=ad_A(Lm(ng)-2,Mm(ng)-1)+              &
     &                            ad_A(-2,-1)
          ad_A(-2,-1)=0.0_r8
!>        tl_A(-2, 0)=tl_A(Lm(ng)-2,Mm(ng)  )
!>
          ad_A(Lm(ng)-2,Mm(ng)  )=ad_A(Lm(ng)-2,Mm(ng)  )+              &
     &                            ad_A(-2, 0)
          ad_A(-2, 0)=0.0_r8
!>        tl_A(-1,-2)=tl_A(Lm(ng)-1,Mm(ng)-2)
!>
          ad_A(Lm(ng)-1,Mm(ng)-2)=ad_A(Lm(ng)-1,Mm(ng)-2)+              &
     &                            ad_A(-1,-2)
          ad_A(-1,-2)=0.0_r8
!>        tl_A(-1,-1)=tl_A(Lm(ng)-1,Mm(ng)-1)
!>
          ad_A(Lm(ng)-1,Mm(ng)-1)=ad_A(Lm(ng)-1,Mm(ng)-1)+              &
     &                            ad_A(-1,-1)
          ad_A(-1,-1)=0.0_r8
!>        tl_A(-1, 0)=tl_A(Lm(ng)-1,Mm(ng)  )
!>
          ad_A(Lm(ng)-1,Mm(ng)  )=ad_A(Lm(ng)-1,Mm(ng)  )+              &
     &                            ad_A(-1, 0)
          ad_A(-1, 0)=0.0_r8
!>        tl_A( 0,-2)=tl_A(Lm(ng)  ,Mm(ng)-2)
!>
          ad_A(Lm(ng)  ,Mm(ng)-2)=ad_A(Lm(ng)  ,Mm(ng)-2)+              &
     &                            ad_A( 0,-2)
          ad_A( 0,-2)=0.0_r8
!>        tl_A( 0,-1)=tl_A(Lm(ng)  ,Mm(ng)-1)
!>
          ad_A(Lm(ng)  ,Mm(ng)-1)=ad_A(Lm(ng)  ,Mm(ng)-1)+              &
     &                            ad_A( 0,-1)
          ad_A( 0,-1)=0.0_r8
!>        tl_A( 0, 0)=tl_A(Lm(ng)  ,Mm(ng)  )
!>
          ad_A(Lm(ng)  ,Mm(ng)  )=ad_A(Lm(ng)  ,Mm(ng)  )+              &
     &                            ad_A( 0, 0)
          ad_A( 0, 0)=0.0_r8
        END IF

        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
!>        tl_A(Lm(ng)+1,-2)=tl_A( 1,Mm(ng)-2)
!>
          ad_A( 1,Mm(ng)-2)=ad_A( 1,Mm(ng)-2)+                          &
     &                      ad_A(Lm(ng)+1,-2)
          ad_A(Lm(ng)+1,-2)=0.0_r8
!>        tl_A(Lm(ng)+1,-1)=tl_A( 1,Mm(ng)-1)
!>
          ad_A( 1,Mm(ng)-1)=ad_A( 1,Mm(ng)-1)+                          &
     &                      ad_A(Lm(ng)+1,-1)
          ad_A(Lm(ng)+1,-1)=0.0_r8
!>        tl_A(Lm(ng)+1, 0)=tl_A( 1,Mm(ng)  )
!>
          ad_A( 1,Mm(ng)  )=ad_A( 1,Mm(ng)  )+                          &
     &                      ad_A(Lm(ng)+1, 0)
          ad_A(Lm(ng)+1, 0)=0.0_r8
!>        tl_A(Lm(ng)+2,-2)=tl_A( 2,Mm(ng)-2)
!>
          ad_A( 2,Mm(ng)-2)=ad_A( 2,Mm(ng)-2)+                          &
     &                      ad_A(Lm(ng)+2,-2)
          ad_A(Lm(ng)+2,-2)=0.0_r8
!>        tl_A(Lm(ng)+2,-1)=tl_A( 2,Mm(ng)-1)
!>
          ad_A( 2,Mm(ng)-1)=ad_A( 2,Mm(ng)-1)+                          &
     &                      ad_A(Lm(ng)+2,-1)
          ad_A(Lm(ng)+2,-1)=0.0_r8
!>        tl_A(Lm(ng)+2, 0)=tl_A( 2,Mm(ng)  )
!>
          ad_A( 2,Mm(ng)  )=ad_A( 2,Mm(ng)  )+                          &
     &                      ad_A(Lm(ng)+2, 0)
          ad_A(Lm(ng)+2, 0)=0.0_r8
#  ifdef THREE_GHOST
!>        tl_A(Lm(ng)+3,-2)=tl_A(3 ,Mm(ng)-2)
!>
          ad_A(3 ,Mm(ng)-2)=ad_A(3 ,Mm(ng)-2)+                          &
     &                      ad_A(Lm(ng)+3,-2)
          ad_A(Lm(ng)+3,-2)=0.0_r8
!>        tl_A(Lm(ng)+3,-1)=tl_A(3 ,Mm(ng)-1)
!>
          ad_A(3 ,Mm(ng)-1)=ad_A(3 ,Mm(ng)-1)+                          &
     &                      ad_A(Lm(ng)+3,-1)
          ad_A(Lm(ng)+3,-1)=0.0_r8
!>        tl_A(Lm(ng)+3, 0)=tl_A(3 ,Mm(ng)  )
!>
          ad_A(3 ,Mm(ng)  )=ad_A(3 ,Mm(ng)  )+                          &
     &                      ad_A(Lm(ng)+3, 0)
          ad_A(Lm(ng)+3, 0)=0.0_r8
#  endif
        END IF

        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
!>        tl_A(-2,Mm(ng)+1)=tl_A(Lm(ng)-2, 1)
!>
          ad_A(Lm(ng)-2, 1)=ad_A(Lm(ng)-2, 1)+                          &
     &                      ad_A(-2,Mm(ng)+1)
          ad_A(-2,Mm(ng)+1)=0.0_r8
!>        tl_A(-1,Mm(ng)+1)=tl_A(Lm(ng)-1, 1)
!>
          ad_A(Lm(ng)-1, 1)=ad_A(Lm(ng)-1, 1)+                          &
     &                      ad_A(-1,Mm(ng)+1)
          ad_A(-1,Mm(ng)+1)=0.0_r8
!>        tl_A( 0,Mm(ng)+1)=tl_A(Lm(ng)  , 1)
!>
          ad_A(Lm(ng)  , 1)=ad_A(Lm(ng)  , 1)+                          &
     &                      ad_A( 0,Mm(ng)+1)
          ad_A( 0,Mm(ng)+1)=0.0_r8
!>        tl_A(-2,Mm(ng)+2)=tl_A(Lm(ng)-2, 2)
!>
          ad_A(Lm(ng)-2, 2)=ad_A(Lm(ng)-2, 2)+                          &
     &                      ad_A(-2,Mm(ng)+2)
          ad_A(-2,Mm(ng)+2)=0.0_r8
!>        tl_A(-1,Mm(ng)+2)=tl_A(Lm(ng)-1, 2)
!>
          ad_A(Lm(ng)-1, 2)=ad_A(Lm(ng)-1, 2)+                          &
     &                      ad_A(-1,Mm(ng)+2)
          ad_A(-1,Mm(ng)+2)=0.0_r8
!>        tl_A( 0,Mm(ng)+2)=tl_A(Lm(ng)  , 2)
!>
          ad_A(Lm(ng)  , 2)=ad_A(Lm(ng)  , 2)+                          &
     &                      ad_A( 0,Mm(ng)+2)
          ad_A( 0,Mm(ng)+2)=0.0_r8
#  ifdef THREE_GHOST
!>        tl_A(-2,Mm(ng)+3)=tl_A(Lm(ng)-2, 3)
!>
          ad_A(Lm(ng)-2, 3)=ad_A(Lm(ng)-2, 3)+                          &
     &                      ad_A(-2,Mm(ng)+3)
          ad_A(-2,Mm(ng)+3)=0.0_r8
!>        tl_A(-1,Mm(ng)+3)=tl_A(Lm(ng)-1, 3)
!>
          ad_A(Lm(ng)-1, 3)=ad_A(Lm(ng)-1, 3)+                          &
     &                      ad_A(-1,Mm(ng)+3)
          ad_A(-1,Mm(ng)+3)=0.0_r8
!>        tl_A( 0,Mm(ng)+3)=tl_A(Lm(ng)  , 3)
!>
          ad_A(Lm(ng)  , 3)=ad_A(Lm(ng)  , 3)+                          &
     &                      ad_A( 0,Mm(ng)+3)
          ad_A( 0,Mm(ng)+3)=0.0_r8
#  endif
        END IF

        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
!>        tl_A(Lm(ng)+1,Mm(ng)+1)=tl_A(1  ,1  )
!>
          ad_A(1  ,1  )=ad_A(1  ,1  )+                                  &
     &                  ad_A(Lm(ng)+1,Mm(ng)+1)
          ad_A(Lm(ng)+1,Mm(ng)+1)=0.0_r8
!>        tl_A(Lm(ng)+1,Mm(ng)+2)=tl_A(1  ,2)
!>
          ad_A(1  ,2)=ad_A(1  ,2)+                                      &
     &                ad_A(Lm(ng)+1,Mm(ng)+2)
          ad_A(Lm(ng)+1,Mm(ng)+2)=0.0_r8
#  ifdef THREE_GHOST
!>        tl_A(Lm(ng)+1,Mm(ng)+3)=tl_A( 1, 3)
!>
          ad_A( 1, 3)=ad_A( 1, 3)+                                      &
     &                ad_A(Lm(ng)+1,Mm(ng)+3)
          ad_A(Lm(ng)+1,Mm(ng)+3)=0.0_r8
#  endif
!>        tl_A(Lm(ng)+2,Mm(ng)+1)=tl_A(2,1  )
!>
          ad_A(2,1  )=ad_A(2,1  )+                                      &
     &                ad_A(Lm(ng)+2,Mm(ng)+1)
          ad_A(Lm(ng)+2,Mm(ng)+1)=0.0_r8
!>        A(Lm(ng)+2,Mm(ng)+2)=A(2,2)
!>
          ad_A(2,2)=ad_A(2,2)+                                          &
     &              ad_A(Lm(ng)+2,Mm(ng)+2)
          ad_A(Lm(ng)+2,Mm(ng)+2)=0.0_r8
#  ifdef THREE_GHOST
!>        tl_A(Lm(ng)+2,Mm(ng)+3)=tl_A( 2, 3)
!>
          ad_A( 2, 3)=ad_A( 2, 3)+                                      &
     &                ad_A(Lm(ng)+2,Mm(ng)+3)
          ad_A(Lm(ng)+2,Mm(ng)+3)=0.0_r8
!>        tl_A(Lm(ng)+3,Mm(ng)+1)=tl_A( 3, 1)
!>
          ad_A( 3, 1)=ad_A( 3, 1)+                                      &
     &                ad_A(Lm(ng)+3,Mm(ng)+1)
          ad_A(Lm(ng)+3,Mm(ng)+1)=0.0_r8
!>        tl_A(Lm(ng)+3,Mm(ng)+2)=tl_A( 3, 2)
!>
          ad_A( 3, 2)=ad_A( 3, 2)+                                      &
     &                ad_A(Lm(ng)+3,Mm(ng)+2)
          ad_A(Lm(ng)+3,Mm(ng)+2)=0.0_r8
!>        tl_A(Lm(ng)+3,Mm(ng)+3)=tl_A( 3, 3)
!>
          ad_A( 3, 3)=ad_A( 3, 3)+                                      &
     &                ad_A(Lm(ng)+3,Mm(ng)+3)
          ad_A(Lm(ng)+3,Mm(ng)+3)=0.0_r8
#  endif
        END IF
      END IF
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
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          DO i=I_RANGE
!>          tl_A(i,-2)=tl_A(i,Mm(ng)-2)
!>
            ad_A(i,Mm(ng)-2)=ad_A(i,Mm(ng)-2)+                          &
     &                       ad_A(i,-2)
            ad_A(i,-2)=0.0_r8
!>          tl_A(i,-1)=tl_A(i,Mm(ng)-1)
!>
            ad_A(i,Mm(ng)-1)=ad_A(i,Mm(ng)-1)+                          &
     &                       ad_A(i,-1)
            ad_A(i,-1)=0.0_r8
!>          tl_A(i, 0)=tl_A(i,Mm(ng)  )
!>
            ad_A(i,Mm(ng)  )=ad_A(i,Mm(ng)  )+                          &
     &                       ad_A(i, 0)
            ad_A(i, 0)=0.0_r8
          END DO
        END IF

        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO i=I_RANGE
!>          tl_A(i,Mm(ng)+1)=tl_A(i,1)
!>
            ad_A(i,1)=ad_A(i,1)+                                        &
     &                ad_A(i,Mm(ng)+1)
            ad_A(i,Mm(ng)+1)=0.0_r8
!>          tl_A(i,Mm(ng)+2)=tl_A(i,2)
!>
            ad_A(i,2)=ad_A(i,2)+                                        &
     &                ad_A(i,Mm(ng)+2)
            ad_A(i,Mm(ng)+2)=0.0_r8
#  ifdef THREE_GHOST
!>          tl_A(i,Mm(ng)+3)=tl_A(i,3)
!>
            ad_A(i,3)=ad_A(i,3)+                                        &
     &                ad_A(i,Mm(ng)+3)
            ad_A(i,Mm(ng)+3)=0.0_r8
#  endif
          END DO
        END IF
      END IF
#  undef I_RANGE
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
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO j=J_RANGE
!>          tl_A(-2,j)=tl_A(Lm(ng)-2,j)
!>
            ad_A(Lm(ng)-2,j)=ad_A(Lm(ng)-2,j)+                          &
     &                       ad_A(-2,j)
            ad_A(-2,j)=0.0_r8
!>          tl_A(-1,j)=tl_A(Lm(ng)-1,j)
!>
            ad_A(Lm(ng)-1,j)=ad_A(Lm(ng)-1,j)+                          &
     &                       ad_A(-1,j)
            ad_A(-1,j)=0.0_r8
!>          tl_A( 0,j)=tl_A(Lm(ng)  ,j)
!>
            ad_A(Lm(ng)  ,j)=ad_A(Lm(ng)  ,j)+                          &
     &                       ad_A( 0,j)
            ad_A( 0,j)=0.0_r8
          END DO
        END IF

        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          DO j=J_RANGE
!>          tl_A(Lm(ng)+1,j)=tl_A(1,j)
!>
            ad_A(1,j)=ad_A(1,j)+                                        &
     &                ad_A(Lm(ng)+1,j)
            ad_A(Lm(ng)+1,j)=0.0_r8
!>          tl_A(Lm(ng)+2,j)=tl_A(2,j)
!>
            ad_A(2,j)=ad_A(2,j)+                                        &
     &                ad_A(Lm(ng)+2,j)
            ad_A(Lm(ng)+2,j)=0.0_r8
#  ifdef THREE_GHOST
!>          tl_A(Lm(ng)+3,j)=tl_A(3,j)
!>
            ad_A(3,j)=ad_A(3,j)+                                        &
     &                ad_A(Lm(ng)+3,j)
            ad_A(Lm(ng)+3,j)=0.0_r8
#  endif
          END DO
        END IF
      END IF
#  undef J_RANGE
# endif

      RETURN
      END SUBROUTINE ad_exchange_r2d_tile

!
!***********************************************************************
      SUBROUTINE ad_exchange_u2d_tile (ng, tile,                        &
     &                                 LBi, UBi, LBj, UBj,              &
     &                                 ad_A)
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
      real(r8), intent(inout) :: ad_A(LBi:,LBj:)
# else
      real(r8), intent(inout) :: ad_A(LBi:UBi,LBj:UBj)
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

# if defined EW_PERIODIC && defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (EW_exchange.and.NS_exchange) THEN
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
!>        tl_A(-2,-2)=tl_A(Lm(ng)-2,Mm(ng)-2)
!>
          ad_A(Lm(ng)-2,Mm(ng)-2)=ad_A(Lm(ng)-2,Mm(ng)-2)+              &
     &                            ad_A(-2,-2)
          ad_A(-2,-2)=0.0_r8
!>        tl_A(-2,-1)=tl_A(Lm(ng)-2,Mm(ng)-1)
!>
          ad_A(Lm(ng)-2,Mm(ng)-1)=ad_A(Lm(ng)-2,Mm(ng)-1)+              &
     &                            ad_A(-2,-1)
          ad_A(-2,-1)=0.0_r8
!>        tl_A(-2, 0)=tl_A(Lm(ng)-2,Mm(ng)  )
!>
          ad_A(Lm(ng)-2,Mm(ng)  )=ad_A(Lm(ng)-2,Mm(ng)  )+              &
     &                            ad_A(-2, 0)
          ad_A(-2, 0)=0.0_r8
!>        tl_A(-1,-2)=tl_A(Lm(ng)-1,Mm(ng)-2)
!>
          ad_A(Lm(ng)-1,Mm(ng)-2)=ad_A(Lm(ng)-1,Mm(ng)-2)+              &
     &                            ad_A(-1,-2)
          ad_A(-1,-2)=0.0_r8
!>        tl_A(-1,-1)=tl_A(Lm(ng)-1,Mm(ng)-1)
!>
          ad_A(Lm(ng)-1,Mm(ng)-1)=ad_A(Lm(ng)-1,Mm(ng)-1)+              &
     &                            ad_A(-1,-1)
          ad_A(-1,-1)=0.0_r8
!>        tl_A(-1, 0)=tl_A(Lm(ng)-1,Mm(ng)  )
!>
          ad_A(Lm(ng)-1,Mm(ng)  )=ad_A(Lm(ng)-1,Mm(ng)  )+              &
     &                            ad_A(-1, 0)
          ad_A(-1, 0)=0.0_r8
!>        tl_A( 0,-2)=tl_A(Lm(ng)  ,Mm(ng)-2)
!>
          ad_A(Lm(ng)  ,Mm(ng)-2)=ad_A(Lm(ng)  ,Mm(ng)-2)+              &
     &                            ad_A( 0,-2)
          ad_A( 0,-2)=0.0_r8
!>        tl_A( 0,-1)=tl_A(Lm(ng)  ,Mm(ng)-1)
!>
          ad_A(Lm(ng)  ,Mm(ng)-1)=ad_A(Lm(ng)  ,Mm(ng)-1)+              &
     &                            ad_A( 0,-1)
          ad_A( 0,-1)=0.0_r8
!>        tl_A( 0, 0)=tl_A(Lm(ng)  ,Mm(ng)  )
!>
          ad_A(Lm(ng)  ,Mm(ng)  )=ad_A(Lm(ng)  ,Mm(ng)  )+              &
     &                            ad_A( 0, 0)
          ad_A( 0, 0)=0.0_r8
        END IF

        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
!>        tl_A(Lm(ng)+1,-2)=tl_A( 1,Mm(ng)-2)
!>
          ad_A( 1,Mm(ng)-2)=ad_A( 1,Mm(ng)-2)+                          &
     &                      ad_A(Lm(ng)+1,-2)
          ad_A(Lm(ng)+1,-2)=0.0_r8
!>        tl_A(Lm(ng)+1,-1)=tl_A( 1,Mm(ng)-1)
!>
          ad_A( 1,Mm(ng)-1)=ad_A( 1,Mm(ng)-1)+                          &
     &                      ad_A(Lm(ng)+1,-1)
          ad_A(Lm(ng)+1,-1)=0.0_r8
!>        tl_A(Lm(ng)+1, 0)=tl_A( 1,Mm(ng)  )
!>
          ad_A( 1,Mm(ng)  )=ad_A( 1,Mm(ng)  )+                          &
     &                      ad_A(Lm(ng)+1, 0)
          ad_A(Lm(ng)+1, 0)=0.0_r8
!>        tl_A(Lm(ng)+2,-2)=tl_A( 2,Mm(ng)-2)
!>
          ad_A( 2,Mm(ng)-2)=ad_A( 2,Mm(ng)-2)+                          &
     &                      ad_A(Lm(ng)+2,-2)
          ad_A(Lm(ng)+2,-2)=0.0_r8
!>        tl_A(Lm(ng)+2,-1)=tl_A( 2,Mm(ng)-1)
!>
          ad_A( 2,Mm(ng)-1)=ad_A( 2,Mm(ng)-1)+                          &
     &                      ad_A(Lm(ng)+2,-1)
          ad_A(Lm(ng)+2,-1)=0.0_r8
!>        tl_A(Lm(ng)+2, 0)=tl_A( 2,Mm(ng)  )
!>
          ad_A( 2,Mm(ng)  )=ad_A( 2,Mm(ng)  )+                          &
     &                      ad_A(Lm(ng)+2, 0)
          ad_A(Lm(ng)+2, 0)=0.0_r8
#  ifdef THREE_GHOST
!>        tl_A(Lm(ng)+3,-2)=tl_A(3 ,Mm(ng)-2)
!>
          ad_A(3 ,Mm(ng)-2)=ad_A(3 ,Mm(ng)-2)+                          &
     &                      ad_A(Lm(ng)+3,-2)
          ad_A(Lm(ng)+3,-2)=0.0_r8
!>        tl_A(Lm(ng)+3,-1)=tl_A(3 ,Mm(ng)-1)
!>
          ad_A(3 ,Mm(ng)-1)=ad_A(3 ,Mm(ng)-1)+                          &
     &                      ad_A(Lm(ng)+3,-1)
          ad_A(Lm(ng)+3,-1)=0.0_r8
!>        tl_A(Lm(ng)+3, 0)=tl_A(3 ,Mm(ng)  )
!>
          ad_A(3 ,Mm(ng)  )=ad_A(3 ,Mm(ng)  )+                          &
     &                      ad_A(Lm(ng)+3, 0)
          ad_A(Lm(ng)+3, 0)=0.0_r8
#  endif
        END IF

        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
!>        tl_A(-2,Mm(ng)+1)=tl_A(Lm(ng)-2, 1)
!>
          ad_A(Lm(ng)-2, 1)=ad_A(Lm(ng)-2, 1)+                          &
     &                      ad_A(-2,Mm(ng)+1)
          ad_A(-2,Mm(ng)+1)=0.0_r8
!>        tl_A(-1,Mm(ng)+1)=tl_A(Lm(ng)-1, 1)
!>
          ad_A(Lm(ng)-1, 1)=ad_A(Lm(ng)-1, 1)+                          &
     &                      ad_A(-1,Mm(ng)+1)
          ad_A(-1,Mm(ng)+1)=0.0_r8
!>        tl_A( 0,Mm(ng)+1)=tl_A(Lm(ng)  , 1)
!>
          ad_A(Lm(ng)  , 1)=ad_A(Lm(ng)  , 1)+                          &
     &                      ad_A( 0,Mm(ng)+1)
          ad_A( 0,Mm(ng)+1)=0.0_r8
!>        tl_A(-2,Mm(ng)+2)=tl_A(Lm(ng)-2, 2)
!>
          ad_A(Lm(ng)-2, 2)=ad_A(Lm(ng)-2, 2)+                          &
     &                      ad_A(-2,Mm(ng)+2)
          ad_A(-2,Mm(ng)+2)=0.0_r8
!>        tl_A(-1,Mm(ng)+2)=tl_A(Lm(ng)-1, 2)
!>
          ad_A(Lm(ng)-1, 2)=ad_A(Lm(ng)-1, 2)+                          &
     &                      ad_A(-1,Mm(ng)+2)
          ad_A(-1,Mm(ng)+2)=0.0_r8
!>        tl_A( 0,Mm(ng)+2)=tl_A(Lm(ng)  , 2)
!>
          ad_A(Lm(ng)  , 2)=ad_A(Lm(ng)  , 2)+                          &
     &                      ad_A( 0,Mm(ng)+2)
          ad_A( 0,Mm(ng)+2)=0.0_r8
#  ifdef THREE_GHOST
!>        tl_A(-2,Mm(ng)+3)=tl_A(Lm(ng)-2, 3)
!>
          ad_A(Lm(ng)-2, 3)=ad_A(Lm(ng)-2, 3)+                          &
     &                      ad_A(-2,Mm(ng)+3)
          ad_A(-2,Mm(ng)+3)=0.0_r8
!>        tl_A(-1,Mm(ng)+3)=tl_A(Lm(ng)-1, 3)
!>
          ad_A(Lm(ng)-1, 3)=ad_A(Lm(ng)-1, 3)+                          &
     &                      ad_A(-1,Mm(ng)+3)
          ad_A(-1,Mm(ng)+3)=0.0_r8
!>        tl_A( 0,Mm(ng)+3)=tl_A(Lm(ng)  , 3)
!>
          ad_A(Lm(ng)  , 3)=ad_A(Lm(ng)  , 3)+                          &
     &                      ad_A( 0,Mm(ng)+3)
          ad_A( 0,Mm(ng)+3)=0.0_r8
#  endif
        END IF

        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
!>        tl_A(Lm(ng)+1,Mm(ng)+1)=tl_A( 1, 1)
!>
          ad_A( 1, 1)=ad_A( 1, 1)+                                      &
     &                ad_A(Lm(ng)+1,Mm(ng)+1)
          ad_A(Lm(ng)+1,Mm(ng)+1)=0.0_r8
!>        tl_A(Lm(ng)+1,Mm(ng)+2)=tl_A( 1, 2)
!>
          ad_A( 1, 2)=ad_A( 1, 2)+                                      &
     &                ad_A(Lm(ng)+1,Mm(ng)+2)
          ad_A(Lm(ng)+1,Mm(ng)+2)=0.0_r8
#  ifdef THREE_GHOST
!>        tl_A(Lm(ng)+1,Mm(ng)+3)=tl_A( 1, 3)
!>
          ad_A( 1, 3)=ad_A( 1, 3)+                                      &
     &                ad_A(Lm(ng)+1,Mm(ng)+3)
          ad_A(Lm(ng)+1,Mm(ng)+3)=0.0_r8
#  endif
!>        tl_A(Lm(ng)+2,Mm(ng)+1)=tl_A( 2, 1)
!>
          ad_A( 2, 1)=ad_A( 2, 1)+                                      &
     &                ad_A(Lm(ng)+2,Mm(ng)+1)
          ad_A(Lm(ng)+2,Mm(ng)+1)=0.0_r8
!>        tl_A(Lm(ng)+2,Mm(ng)+2)=tl_A( 2, 2)
!>
          ad_A( 2, 2)=ad_A( 2, 2)+                                      &
     &                ad_A(Lm(ng)+2,Mm(ng)+2)
          ad_A(Lm(ng)+2,Mm(ng)+2)=0.0_r8
#  ifdef THREE_GHOST
!>        tl_A(Lm(ng)+2,Mm(ng)+3)=tl_A( 2, 3)
!>
          ad_A( 2, 3)=ad_A( 2, 3)+                                      &
     &                ad_A(Lm(ng)+2,Mm(ng)+3)
          ad_A(Lm(ng)+2,Mm(ng)+3)=0.0_r8
!>        tl_A(Lm(ng)+3,Mm(ng)+1)=tl_A( 3, 1)
!>
          ad_A( 3, 1)=ad_A( 3, 1)+                                      &
     &                ad_A(Lm(ng)+3,Mm(ng)+1)
          ad_A(Lm(ng)+3,Mm(ng)+1)=0.0_r8
!>        tl_A(Lm(ng)+3,Mm(ng)+2)=tl_A( 3, 2)
!>
          ad_A( 3, 2)=ad_A( 3, 2)+                                      &
     &                ad_A(Lm(ng)+3,Mm(ng)+2)
          ad_A(Lm(ng)+3,Mm(ng)+2)=0.0_r8
!>        tl_A(Lm(ng)+3,Mm(ng)+3)=tl_A( 3, 3)
!>
          ad_A( 3, 3)=ad_A( 3, 3)+                                      &
     &                ad_A(Lm(ng)+3,Mm(ng)+3)
          ad_A(Lm(ng)+3,Mm(ng)+3)=0.0_r8
#  endif
        END IF
      END IF
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
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          DO i=I_RANGE
!>          tl_A(i,-2)=tl_A(i,Mm(ng)-2)
!>
            ad_A(i,Mm(ng)-2)=ad_A(i,Mm(ng)-2)+                          &
     &                       ad_A(i,-2)
            ad_A(i,-2)=0.0_r8
!>          tl_A(i,-1)=tl_A(i,Mm(ng)-1)
!>
            ad_A(i,Mm(ng)-1)=ad_A(i,Mm(ng)-1)+                          &
     &                       ad_A(i,-1)
            ad_A(i,-1)=0.0_r8
!>          tl_A(i, 0)=tl_A(i,Mm(ng)  )
!>
            ad_A(i,Mm(ng)  )=ad_A(i,Mm(ng)  )+                          &
     &                       ad_A(i, 0)
            ad_A(i, 0)=0.0_r8
          END DO
        END IF

        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO i=I_RANGE
!>          tl_A(i,Mm(ng)+1)=tl_A(i,1)
!>
            ad_A(i,1)=ad_A(i,1)+                                        &
     &                ad_A(i,Mm(ng)+1)
            ad_A(i,Mm(ng)+1)=0.0_r8
!>          tl_A(i,Mm(ng)+2)=tl_A(i,2)
!>
            ad_A(i,2)=ad_A(i,2)+                                        &
     &                ad_A(i,Mm(ng)+2)
            ad_A(i,Mm(ng)+2)=0.0_r8
#  ifdef THREE_GHOST
!>          tl_A(i,Mm(ng)+3)=tl_A(i,3)
!>
            ad_A(i,3)=ad_A(i,3)+                                        &
     &                ad_A(i,Mm(ng)+3)
            ad_A(i,Mm(ng)+3)=0.0_r8
#  endif
          END DO
        END IF
      END IF
#  undef I_RANGE
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
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO j=J_RANGE
!>          tl_A(-2,j)=tl_A(Lm(ng)-2,j)
!>
            ad_A(Lm(ng)-2,j)=ad_A(Lm(ng)-2,j)+                          &
     &                       ad_A(-2,j)
            ad_A(-2,j)=0.0_r8
!>          tl_A(-1,j)=tl_A(Lm(ng)-1,j)
!>
            ad_A(Lm(ng)-1,j)=ad_A(Lm(ng)-1,j)+                          &
     &                       ad_A(-1,j)
            ad_A(-1,j)=0.0_r8
!>          tl_A( 0,j)=tl_A(Lm(ng)  ,j)
!>
            ad_A(Lm(ng)  ,j)=ad_A(Lm(ng)  ,j)+                          &
     &                       ad_A( 0,j)
            ad_A( 0,j)=0.0_r8
          END DO
        END IF

        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          DO j=J_RANGE
!>          tl_A(Lm(ng)+1,j)=tl_A(1,j)
!>
            ad_A(1,j)=ad_A(1,j)+                                        &
     &                ad_A(Lm(ng)+1,j)
            ad_A(Lm(ng)+1,j)=0.0_r8
!>          tl_A(Lm(ng)+2,j)=tl_A(2,j)
!>
            ad_A(2,j)=ad_A(2,j)+                                        &
     &                ad_A(Lm(ng)+2,j)
            ad_A(Lm(ng)+2,j)=0.0_r8
#  ifdef THREE_GHOST
!>          tl_A(Lm(ng)+3,j)=tl_A(3,j)
!>
            ad_A(3,j)=ad_A(3,j)+                                        &
     &                ad_A(Lm(ng)+3,j)
            ad_A(Lm(ng)+3,j)=0.0_r8
#  endif
          END DO
        END IF
      END IF
#  undef J_RANGE
# endif

      RETURN
      END SUBROUTINE ad_exchange_u2d_tile

!
!***********************************************************************
      SUBROUTINE ad_exchange_v2d_tile (ng, tile,                        &
     &                                 LBi, UBi, LBj, UBj,              &
     &                                 ad_A)
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
      real(r8), intent(inout) :: ad_A(LBi:,LBj:)
# else
      real(r8), intent(inout) :: ad_A(LBi:UBi,LBj:UBj)
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

# if defined EW_PERIODIC && defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (EW_exchange.and.NS_exchange) THEN
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
!>        tl_A(-2,-2)=tl_A(Lm(ng)-2,Mm(ng)-2)
!>
          ad_A(Lm(ng)-2,Mm(ng)-2)=ad_A(Lm(ng)-2,Mm(ng)-2)+              &
     &                            ad_A(-2,-2)
          ad_A(-2,-2)=0.0_r8
!>        tl_A(-2,-1)=tl_A(Lm(ng)-2,Mm(ng)-1)
!>
          ad_A(Lm(ng)-2,Mm(ng)-1)=ad_A(Lm(ng)-2,Mm(ng)-1)+              &
     &                            ad_A(-2,-1)
          ad_A(-2,-1)=0.0_r8
!>        tl_A(-2, 0)=tl_A(Lm(ng)-2,Mm(ng)  )
!>
          ad_A(Lm(ng)-2,Mm(ng)  )=ad_A(Lm(ng)-2,Mm(ng)  )+              &
     &                            ad_A(-2, 0)
          ad_A(-2, 0)=0.0_r8
!>        tl_A(-1,-2)=tl_A(Lm(ng)-1,Mm(ng)-2)
!>
          ad_A(Lm(ng)-1,Mm(ng)-2)=ad_A(Lm(ng)-1,Mm(ng)-2)+              &
     &                            ad_A(-1,-2)
          ad_A(-1,-2)=0.0_r8
!>        tl_A(-1,-1)=tl_A(Lm(ng)-1,Mm(ng)-1)
!>
          ad_A(Lm(ng)-1,Mm(ng)-1)=ad_A(Lm(ng)-1,Mm(ng)-1)+              &
     &                            ad_A(-1,-1)
          ad_A(-1,-1)=0.0_r8
!>        tl_A(-1, 0)=tl_A(Lm(ng)-1,Mm(ng)  )
!>
          ad_A(Lm(ng)-1,Mm(ng)  )=ad_A(Lm(ng)-1,Mm(ng)  )+              &
     &                            ad_A(-1, 0)
          ad_A(-1, 0)=0.0_r8
!>        tl_A( 0,-2)=tl_A(Lm(ng)  ,Mm(ng)-2)
!>
          ad_A(Lm(ng)  ,Mm(ng)-2)=ad_A(Lm(ng)  ,Mm(ng)-2)+              &
     &                            ad_A( 0,-2)
          ad_A( 0,-2)=0.0_r8
!>        tl_A( 0,-1)=tl_A(Lm(ng)  ,Mm(ng)-1)
!>
          ad_A(Lm(ng)  ,Mm(ng)-1)=ad_A(Lm(ng)  ,Mm(ng)-1)+              &
     &                            ad_A( 0,-1)
          ad_A( 0,-1)=0.0_r8
!>        tl_A( 0, 0)=tl_A(Lm(ng)  ,Mm(ng)  )
!>
          ad_A(Lm(ng)  ,Mm(ng)  )=ad_A(Lm(ng)  ,Mm(ng)  )+              &
     &                            ad_A( 0, 0)
          ad_A( 0, 0)=0.0_r8
        END IF

        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
!>        tl_A(Lm(ng)+1,-2)=tl_A( 1,Mm(ng)-2)
!>
          ad_A( 1,Mm(ng)-2)=ad_A( 1,Mm(ng)-2)+                          &
     &                      ad_A(Lm(ng)+1,-2)
          ad_A(Lm(ng)+1,-2)=0.0_r8
!>        tl_A(Lm(ng)+1,-1)=tl_A( 1,Mm(ng)-1)
!>
          ad_A( 1,Mm(ng)-1)=ad_A( 1,Mm(ng)-1)+                          &
     &                      ad_A(Lm(ng)+1,-1)
          ad_A(Lm(ng)+1,-1)=0.0_r8
!>        tl_A(Lm(ng)+1, 0)=tl_A( 1,Mm(ng)  )
!>
          ad_A( 1,Mm(ng)  )=ad_A( 1,Mm(ng)  )+                          &
     &                      ad_A(Lm(ng)+1, 0)
          ad_A(Lm(ng)+1, 0)=0.0_r8
!>        tl_A(Lm(ng)+2,-2)=tl_A( 2,Mm(ng)-2)
!>
          ad_A( 2,Mm(ng)-2)=ad_A( 2,Mm(ng)-2)+                          &
     &                      ad_A(Lm(ng)+2,-2)
          ad_A(Lm(ng)+2,-2)=0.0_r8
!>        tl_A(Lm(ng)+2,-1)=tl_A( 2,Mm(ng)-1)
!>
          ad_A( 2,Mm(ng)-1)=ad_A( 2,Mm(ng)-1)+                          &
     &                      ad_A(Lm(ng)+2,-1)
          ad_A(Lm(ng)+2,-1)=0.0_r8
!>        tl_A(Lm(ng)+2, 0)=tl_A( 2,Mm(ng)  )
!>
          ad_A( 2,Mm(ng)  )=ad_A( 2,Mm(ng)  )+                          &
     &                      ad_A(Lm(ng)+2, 0)
          ad_A(Lm(ng)+2, 0)=0.0_r8
#  ifdef THREE_GHOST
!>        tl_A(Lm(ng)+3,-2)=tl_A(3 ,Mm(ng)-2)
!>
          ad_A(3 ,Mm(ng)-2)=ad_A(3 ,Mm(ng)-2)+                          &
     &                      ad_A(Lm(ng)+3,-2)
          ad_A(Lm(ng)+3,-2)=0.0_r8
!>        tl_A(Lm(ng)+3,-1)=tl_A(3 ,Mm(ng)-1)
!>
          ad_A(3 ,Mm(ng)-1)=ad_A(3 ,Mm(ng)-1)+                          &
     &                      ad_A(Lm(ng)+3,-1)
          ad_A(Lm(ng)+3,-1)=0.0_r8
!>        tl_A(Lm(ng)+3, 0)=tl_A(3 ,Mm(ng)  )
!>
          ad_A(3 ,Mm(ng)  )=ad_A(3 ,Mm(ng)  )+                          &
     &                      ad_A(Lm(ng)+3, 0)
          ad_A(Lm(ng)+3, 0)=0.0_r8
#  endif
        END IF

        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
!>        tl_A(-2,Mm(ng)+1)=tl_A(Lm(ng)-2, 1)
!>
          ad_A(Lm(ng)-2, 1)=ad_A(Lm(ng)-2, 1)+                          &
     &                      ad_A(-2,Mm(ng)+1)
          ad_A(-2,Mm(ng)+1)=0.0_r8
!>        tl_A(-1,Mm(ng)+1)=tl_A(Lm(ng)-1, 1)
!>
          ad_A(Lm(ng)-1, 1)=ad_A(Lm(ng)-1, 1)+                          &
     &                      ad_A(-1,Mm(ng)+1)
          ad_A(-1,Mm(ng)+1)=0.0_r8
!>        tl_A( 0,Mm(ng)+1)=tl_A(Lm(ng)  , 1)
!>
          ad_A(Lm(ng)  , 1)=ad_A(Lm(ng)  , 1)+                          &
     &                      ad_A( 0,Mm(ng)+1)
          ad_A( 0,Mm(ng)+1)=0.0_r8
!>        tl_A(-2,Mm(ng)+2)=tl_A(Lm(ng)-2, 2)
!>
          ad_A(Lm(ng)-2, 2)=ad_A(Lm(ng)-2, 2)+                          &
     &                      ad_A(-2,Mm(ng)+2)
          ad_A(-2,Mm(ng)+2)=0.0_r8
!>        tl_A(-1,Mm(ng)+2)=tl_A(Lm(ng)-1, 2)
!>
          ad_A(Lm(ng)-1, 2)=ad_A(Lm(ng)-1, 2)+                          &
     &                      ad_A(-1,Mm(ng)+2)
          ad_A(-1,Mm(ng)+2)=0.0_r8
!>        tl_A( 0,Mm(ng)+2)=tl_A(Lm(ng)  , 2)
!>
          ad_A(Lm(ng)  , 2)=ad_A(Lm(ng)  , 2)+                          &
     &                      ad_A( 0,Mm(ng)+2)
          ad_A( 0,Mm(ng)+2)=0.0_r8
#  ifdef THREE_GHOST
!>        tl_A(-2,Mm(ng)+3)=tl_A(Lm(ng)-2, 3)
!>
          ad_A(Lm(ng)-2, 3)=ad_A(Lm(ng)-2, 3)+                          &
     &                      ad_A(-2,Mm(ng)+3)
          ad_A(-2,Mm(ng)+3)=0.0_r8
!>        tl_A(-1,Mm(ng)+3)=tl_A(Lm(ng)-1, 3)
!>
          ad_A(Lm(ng)-1, 3)=ad_A(Lm(ng)-1, 3)+                          &
     &                      ad_A(-1,Mm(ng)+3)
          ad_A(-1,Mm(ng)+3)=0.0_r8
!>        tl_A( 0,Mm(ng)+3)=tl_A(Lm(ng)  , 3)
!>
          ad_A(Lm(ng)  , 3)=ad_A(Lm(ng)  , 3)+                          &
     &                      ad_A( 0,Mm(ng)+3)
          ad_A( 0,Mm(ng)+3)=0.0_r8
#  endif
        END IF

        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
!>        tl_A(Lm(ng)+1,Mm(ng)+1)=tl_A( 1, 1)
!>
          ad_A( 1, 1)=ad_A( 1, 1)+                                      &
     &                ad_A(Lm(ng)+1,Mm(ng)+1)
          ad_A(Lm(ng)+1,Mm(ng)+1)=0.0_r8
!>        tl_A(Lm(ng)+1,Mm(ng)+2)=tl_A( 1 , 2)
!>
          ad_A( 1, 2)=ad_A( 1, 2)+                                      &
     &                ad_A(Lm(ng)+1,Mm(ng)+2)
          ad_A(Lm(ng)+1,Mm(ng)+2)=0.0_r8
#  ifdef THREE_GHOST
!>        tl_A(Lm(ng)+1,Mm(ng)+3)=tl_A( 1, 3)
!>
          ad_A( 1, 3)=ad_A( 1, 3)+                                      &
     &                ad_A(Lm(ng)+1,Mm(ng)+3)
          ad_A(Lm(ng)+1,Mm(ng)+3)=0.0_r8
#  endif
!>        tl_A(Lm(ng)+2,Mm(ng)+1)=tl_A(2,1  )
!>
          ad_A( 2, 1)=ad_A( 2, 1 )+                                     &
     &                ad_A(Lm(ng)+2,Mm(ng)+1)
          ad_A(Lm(ng)+2,Mm(ng)+1)=0.0_r8
!>        tl_A(Lm(ng)+2,Mm(ng)+2)=tl_A( 2, 2)
!>
          ad_A( 2, 2)=ad_A( 2, 2)+                                      &
     &                ad_A(Lm(ng)+2,Mm(ng)+2)
          ad_A(Lm(ng)+2,Mm(ng)+2)=0.0_r8
#  ifdef THREE_GHOST
!>        tl_A(Lm(ng)+2,Mm(ng)+3)=tl_A( 2, 3)
!>
          ad_A( 2, 3)=ad_A( 2, 3)+                                      &
     &                ad_A(Lm(ng)+2,Mm(ng)+3)
          ad_A(Lm(ng)+2,Mm(ng)+3)=0.0_r8
!>        tl_A(Lm(ng)+3,Mm(ng)+1)=tl_A( 3, 1)
!>
          ad_A( 3, 1)=ad_A( 3, 1)+                                      &
     &                ad_A(Lm(ng)+3,Mm(ng)+1)
          ad_A(Lm(ng)+3,Mm(ng)+1)=0.0_r8
!>        tl_A(Lm(ng)+3,Mm(ng)+2)=tl_A( 3, 2)
!>
          ad_A( 3, 2)=ad_A( 3, 2)+                                      &
     &                ad_A(Lm(ng)+3,Mm(ng)+2)
          ad_A(Lm(ng)+3,Mm(ng)+2)=0.0_r8
!>        tl_A(Lm(ng)+3,Mm(ng)+3)=tl_A( 3, 3)
!>
          ad_A( 3, 3)=ad_A( 3, 3)+                                      &
     &                ad_A(Lm(ng)+3,Mm(ng)+3)
          ad_A(Lm(ng)+3,Mm(ng)+3)=0.0_r8
#  endif
        END IF
      END IF
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
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          DO i=I_RANGE
!>          tl_A(i,-2)=tl_A(i,Mm(ng)-2)
!>
            ad_A(i,Mm(ng)-2)=ad_A(i,Mm(ng)-2)+                          &
     &                       ad_A(i,-2)
            ad_A(i,-2)=0.0_r8
!>          tl_A(i,-1)=tl_A(i,Mm(ng)-1)
!>
            ad_A(i,Mm(ng)-1)=ad_A(i,Mm(ng)-1)+                          &
     &                       ad_A(i,-1)
            ad_A(i,-1)=0.0_r8
!>          tl_A(i, 0)=tl_A(i,Mm(ng)  )
!>
            ad_A(i,Mm(ng)  )=ad_A(i,Mm(ng)  )+                          &
     &                       ad_A(i, 0)
            ad_A(i, 0)=0.0_r8
          END DO
        END IF

        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO i=I_RANGE
!>          tl_A(i,Mm(ng)+1)=tl_A(i,1)
!>
            ad_A(i,1)=ad_A(i,1)+                                        &
     &                ad_A(i,Mm(ng)+1)
            ad_A(i,Mm(ng)+1)=0.0_r8
!>          tl_A(i,Mm(ng)+2)=tl_A(i,2)
!>
            ad_A(i,2)=ad_A(i,2)+                                        &
     &                ad_A(i,Mm(ng)+2)
            ad_A(i,Mm(ng)+2)=0.0_r8
#  ifdef THREE_GHOST
!>          tl_A(i,Mm(ng)+3)=tl_A(i,3)
!>
            ad_A(i,3)=ad_A(i,3)+                                        &
     &                ad_A(i,Mm(ng)+3)
            ad_A(i,Mm(ng)+3)=0.0_r8
#  endif
          END DO
        END IF
      END IF
#  undef I_RANGE
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
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO j=J_RANGE
!>          tl_A(-2,j)=tl_A(Lm(ng)-2,j)
!>
            ad_A(Lm(ng)-2,j)=ad_A(Lm(ng)-2,j)+                          &
     &                       ad_A(-2,j)
            ad_A(-2,j)=0.0_r8
!>          tl_A(-1,j)=tl_A(Lm(ng)-1,j)
!>
            ad_A(Lm(ng)-1,j)=ad_A(Lm(ng)-1,j)+                          &
     &                       ad_A(-1,j)
            ad_A(-1,j)=0.0_r8
!>          tl_A( 0,j)=tl_A(Lm(ng)  ,j)
!>
            ad_A(Lm(ng)  ,j)=ad_A(Lm(ng)  ,j)+                          &
     &                       ad_A( 0,j)
            ad_A( 0,j)=0.0_r8
          END DO
        END IF

        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          DO j=J_RANGE
!>          tl_A(Lm(ng)+1,j)=tl_A(1  ,j)
!>
            ad_A(1  ,j)=ad_A(1  ,j)+                                    &
     &                  ad_A(Lm(ng)+1,j)
            ad_A(Lm(ng)+1,j)=0.0_r8
!>          tl_A(Lm(ng)+2,j)=tl_A(2,j)
!>
            ad_A(2,j)=ad_A(2,j)+                                        &
     &                ad_A(Lm(ng)+2,j)
            ad_A(Lm(ng)+2,j)=0.0_r8
#  ifdef THREE_GHOST
!>          tl_A(Lm(ng)+3,j)=tl_A(3,j)
!>
            ad_A(3,j)=ad_A(3,j)+                                        &
     &                ad_A(Lm(ng)+3,j)
            ad_A(Lm(ng)+3,j)=0.0_r8
#  endif
          END DO
        END IF
      END IF
#  undef J_RANGE
# endif

      RETURN
      END SUBROUTINE ad_exchange_v2d_tile
#endif
      END MODULE ad_exchange_2d_mod
