#include "cppdefs.h"
      MODULE ad_bc_2d_mod
#ifdef ADJOINT
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2011 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  These routines apply close, gradient or periodic boundary           !
!  conditions to generic 2D adjoint fields.                            !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng                Nested grid number.                            !
!     tile              Domain partition.                              !
!     LBi               I-dimension Lower bound.                       !
!     UBi               I-dimension Upper bound.                       !
!     LBj               J-dimension Lower bound.                       !
!     UBj               J-dimension Upper bound.                       !
!     ad_A              2D adjoint field.                              !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     ad_A              Processed 2D adjoint field.                    !
!                                                                      !
!  Routines:                                                           !
!                                                                      !
!     ad_bc_r2d_tile    Boundary conditions for field at RHO-points    !
!     ad_bc_u2d_tile    Boundary conditions for field at U-points      !
!     ad_bc_v2d_tile    Boundary conditions for field at V-points      !
!                                                                      !
!=======================================================================
!
      implicit none

      CONTAINS

!
!***********************************************************************
      SUBROUTINE ad_bc_r2d_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           ad_A)
!***********************************************************************
!
      USE mod_param

# if defined EW_PERIODIC || defined NS_PERIODIC
!
      USE ad_exchange_2d_mod, ONLY : ad_exchange_r2d_tile
# endif
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
      integer :: i, j

      real(r8) :: adfac

# include "set_bounds.h"

# if defined EW_PERIODIC || defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Apply adjoint periodic boundary conditons.
!-----------------------------------------------------------------------
!
      CALL ad_exchange_r2d_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           ad_A)
# endif

# if !defined EW_PERIODIC && !defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
!>      tl_A(Iend+1,Jend+1)=0.5_r8*(tl_A(Iend+1,Jend  )+                &
!>   &                              tl_A(Iend  ,Jend+1))
!>
        adfac=0.5_r8*ad_A(Iend+1,Jend+1)
        ad_A(Iend+1,Jend  )=ad_A(Iend+1,Jend  )+adfac
        ad_A(Iend  ,Jend+1)=ad_A(Iend  ,Jend+1)+adfac
        ad_A(Iend+1,Jend+1)=0.0_r8
      END IF
      IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
!>      tl_A(Istr-1,Jend+1)=0.5_r8*(tl_A(Istr-1,Jend  )+                &
!>   &                              tl_A(Istr  ,Jend+1))
!>
        adfac=0.5_r8*ad_A(Istr-1,Jend+1)
        ad_A(Istr-1,Jend  )=ad_A(Istr-1,Jend  )+adfac
        ad_A(Istr  ,Jend+1)=ad_A(Istr  ,Jend+1)+adfac
        ad_A(Istr-1,Jend+1)=0.0_r8
      END IF
      IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
!>      tl_A(Iend+1,Jstr-1)=0.5_r8*(tl_A(Iend  ,Jstr-1)+                &
!>   &                              tl_A(Iend+1,Jstr  ))
!>
        adfac=0.5_r8*ad_A(Iend+1,Jstr-1)
        ad_A(Iend  ,Jstr-1)=ad_A(Iend  ,Jstr-1)+adfac
        ad_A(Iend+1,Jstr  )=ad_A(Iend+1,Jstr  )+adfac
        ad_A(Iend+1,Jstr-1)=0.0_r8
      END IF
      IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
!>      tl_A(Istr-1,Jstr-1)=0.5_r8*(tl_A(Istr  ,Jstr-1)+                &
!>   &                              tl_A(Istr-1,Jstr  ))
!>
        adfac=0.5_r8*ad_A(Istr-1,Jstr-1)
        ad_A(Istr  ,Jstr-1)=ad_A(Istr  ,Jstr-1)+adfac
        ad_A(Istr-1,Jstr  )=ad_A(Istr-1,Jstr  )+adfac
        ad_A(Istr-1,Jstr-1)=0.0_r8
      END IF
# endif

# ifndef NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Adjoint North-South gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
        DO i=Istr,Iend
!>        tl_A(i,Jstr-1)=tl_A(i,Jstr)
!>
          ad_A(i,Jstr  )=ad_A(i,Jstr)+ad_A(i,Jstr-1)
          ad_A(i,Jstr-1)=0.0_r8
        END DO
      END IF
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
        DO i=Istr,Iend
!>        tl_A(i,Jend+1)=tl_A(i,Jend)
!>
          ad_A(i,Jend  )=ad_A(i,Jend)+ad_A(i,Jend+1)
          ad_A(i,Jend+1)=0.0_r8
        END DO
      END IF
# endif

# ifndef EW_PERIODIC
!
!-----------------------------------------------------------------------
!  Adjoint East-West gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
        DO j=Jstr,Jend
!>        tl_A(Istr-1,j)=tl_A(Istr,j)
!>
          ad_A(Istr  ,j)=ad_A(Istr,j)+ad_A(Istr-1,j)
          ad_A(Istr-1,j)=0.0_r8
        END DO
      END IF
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
        DO j=Jstr,Jend
!>        tl_A(Iend+1,j)=tl_A(Iend,j)
!>
          ad_A(Iend  ,j)=ad_A(Iend,j)+ad_A(Iend+1,j)
          ad_A(Iend+1,j)=0.0_r8
        END DO
      END IF
# endif

      RETURN
      END SUBROUTINE ad_bc_r2d_tile

!
!***********************************************************************
      SUBROUTINE ad_bc_u2d_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           ad_A)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_scalars

# if defined EW_PERIODIC || defined NS_PERIODIC
!
      USE ad_exchange_2d_mod, ONLY : ad_exchange_u2d_tile
# endif
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
      integer :: i, j

      real(r8) :: adfac

# include "set_bounds.h"

# if defined EW_PERIODIC || defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Apply adjoint periodic boundary conditons.
!-----------------------------------------------------------------------
!
      CALL ad_exchange_u2d_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           ad_A)
# endif

# if !defined EW_PERIODIC && !defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
!>      tl_A(Iend+1,Jend+1)=0.5_r8*(tl_A(Iend+1,Jend  )+                &
!>   &                              tl_A(Iend  ,Jend+1))
!>
        adfac=0.5_r8*ad_A(Iend+1,Jend+1)
        ad_A(Iend+1,Jend  )=ad_A(Iend+1,Jend  )+adfac
        ad_A(Iend  ,Jend+1)=ad_A(Iend  ,Jend+1)+adfac
        ad_A(Iend+1,Jend+1)=0.0_r8
      END IF
      IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
!>      tl_A(Istr  ,Jend+1)=0.5_r8*(tl_A(Istr  ,Jend  )+                 &
!>   &                              tl_A(Istr+1,Jend+1))
!>
        adfac=0.5_r8*ad_A(Istr,Jend+1)
        ad_A(Istr  ,Jend  )=ad_A(Istr  ,Jend  )+adfac
        ad_A(Istr+1,Jend+1)=ad_A(Istr+1,Jend+1)+adfac
        ad_A(Istr  ,Jend+1)=0.0_r8
      END IF
      IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
!>      tl_A(Iend+1,Jstr-1)=0.5_r8*(tl_A(Iend  ,Jstr-1)+                &
!>   &                              tl_A(Iend+1,Jstr  ))
!>
        adfac=0.5_r8*ad_A(Iend+1,Jstr-1)
        ad_A(Iend  ,Jstr-1)=ad_A(Iend  ,Jstr-1)+adfac
        ad_A(Iend+1,Jstr  )=ad_A(Iend+1,Jstr  )+adfac
        ad_A(Iend+1,Jstr-1)=0.0_r8
      END IF
      IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
!>      tl_A(Istr  ,Jstr-1)=0.5_r8*(tl_A(Istr+1,Jstr-1)+                &
!>   &                              tl_A(Istr  ,Jstr  ))
!>
        adfac=0.5_r8*ad_A(Istr,Jstr-1)
        ad_A(Istr+1,Jstr-1)=ad_A(Istr+1,Jstr-1)+adfac
        ad_A(Istr  ,Jstr  )=ad_A(Istr  ,Jstr  )+adfac
        ad_A(Istr  ,Jstr-1)=0.0_r8
      END IF
# endif

# ifndef NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Adjoint North-South boundary conditions: Closed (free-slip/no-slip)
!  or gradient.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
#  ifdef SOUTHERN_WALL
#   ifdef EW_PERIODIC
#    define I_RANGE IstrU,Iend
#   else
#    define I_RANGE Istr,IendR
#   endif
        DO i=I_RANGE
#   ifdef MASKING
!>        tl_A(i,Jstr-1)=tl_A(i,Jstr-1)*GRID(ng)%umask(i,Jstr-1)
!>
          ad_A(i,Jstr-1)=ad_A(i,Jstr-1)*GRID(ng)%umask(i,Jstr-1)
#   endif
!>        tl_A(i,Jstr-1)=gamma2(ng)*tl_A(i,Jstr)
!>
          ad_A(i,Jstr)=ad_A(i,Jstr)+gamma2(ng)*ad_A(i,Jstr-1)
          ad_A(i,Jstr-1)=0.0_r8
        END DO
#   undef I_RANGE
#  else
        DO i=IstrU,Iend
!>        tl_A(i,Jstr-1)=tl_A(i,Jstr)
!>
          ad_A(i,Jstr  )=ad_A(i,Jstr)+ad_A(i,Jstr-1)
          ad_A(i,Jstr-1)=0.0_r8
        END DO
#  endif
      END IF

      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
#  ifdef NORTHERN_WALL
#   ifdef EW_PERIODIC
#    define I_RANGE IstrU,Iend
#   else
#    define I_RANGE Istr,IendR
#   endif
        DO i=I_RANGE
#   ifdef MASKING
!>        tl_A(i,Jend+1)=tl_A(i,Jend+1)*GRID(ng)%umask(i,Jend+1)
!>
          ad_A(i,Jend+1)=ad_A(i,Jend+1)*GRID(ng)%umask(i,Jend+1)
#   endif
!>        tl_A(i,Jend+1)=gamma2(ng)*tl_A(i,Jend)
!>
          ad_A(i,Jend)=ad_A(i,Jend)+gamma2(ng)*ad_A(i,Jend+1)
          ad_A(i,Jend+1)=0.0_r8
        END DO
#   undef I_RANGE
#  else
        DO i=IstrU,Iend
!>        tl_A(i,Jend+1)=tl_A(i,Jend)
!>
          ad_A(i,Jend  )=ad_A(i,Jend)+ad_A(i,Jend+1)
          ad_A(i,Jend+1)=0.0_r8
        END DO
#  endif
      END IF
# endif

# ifndef EW_PERIODIC
!
!-----------------------------------------------------------------------
!  Adjoint East-West boundary conditions: Closed or gradient.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
        DO j=Jstr,Jend
#  ifdef WESTERN_WALL
!>        tl_A(Istr,j)=0.0_r8
!>
          ad_A(Istr,j)=0.0_r8
#  else
!>        tl_A(Istr,j)=tl_A(Istr+1,j)
!>
          ad_A(Istr+1,j)=ad_A(Istr+1,j)+ad_A(Istr,j)
          ad_A(Istr  ,j)=0.0_r8
#  endif
        END DO
      END IF
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
        DO j=Jstr,Jend
#  ifdef EASTERN_WALL
!>        A(Iend+1,j)=0.0_r8
!>
          ad_A(Iend+1,j)=0.0_r8
#  else
!>        tl_A(Iend+1,j)=tl_A(Iend,j)
!>
          ad_A(Iend  ,j)=ad_A(Iend,j)+ad_A(Iend+1,j)
          ad_A(Iend+1,j)=0.0_r8
#  endif
        END DO
      END IF
# endif

      RETURN
      END SUBROUTINE ad_bc_u2d_tile

!
!***********************************************************************
      SUBROUTINE ad_bc_v2d_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           ad_A)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_scalars

# if defined EW_PERIODIC || defined NS_PERIODIC
!
      USE ad_exchange_2d_mod, ONLY : ad_exchange_v2d_tile
# endif
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
      integer :: i, j

      real(r8) :: adfac

# include "set_bounds.h"

# if defined EW_PERIODIC || defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Apply adjoint periodic boundary conditons.
!-----------------------------------------------------------------------
!
      CALL ad_exchange_v2d_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           ad_A)
# endif

# if !defined EW_PERIODIC && !defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
!>      tl_A(Iend+1,Jend+1)=0.5_r8*(tl_A(Iend+1,Jend  )+                &
!>   &                              tl_A(Iend  ,Jend+1))
!>
        adfac=0.5_r8*ad_A(Iend+1,Jend+1)
        ad_A(Iend+1,Jend  )=ad_A(Iend+1,Jend  )+adfac
        ad_A(Iend  ,Jend+1)=ad_A(Iend  ,Jend+1)+adfac
        ad_A(Iend+1,Jend+1)=0.0_r8
      END IF
      IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
!>      tl_A(Istr-1,Jend+1)=0.5_r8*(tl_A(Istr-1,Jend  )+                &
!>   &                              tl_A(Istr  ,Jend+1))
!>
        adfac=0.5_r8*ad_A(Istr-1,Jend+1)
        ad_A(Istr-1,Jend  )=ad_A(Istr-1,Jend  )+adfac
        ad_A(Istr  ,Jend+1)=ad_A(Istr  ,Jend+1)+adfac
        ad_A(Istr-1,Jend+1)=0.0_r8
      END IF
      IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
!>      tl_A(Iend+1,Jstr  )=0.5_r8*(tl_A(Iend  ,Jstr  )+                &
!>   &                              tl_A(Iend+1,Jstr+1))
!>
        adfac=0.5_r8*ad_A(Iend+1,Jstr  )
        ad_A(Iend  ,Jstr  )=ad_A(Iend  ,Jstr  )+adfac
        ad_A(Iend+1,Jstr+1)=ad_A(Iend+1,Jstr+1)+adfac
        ad_A(Iend+1,Jstr  )=0.0_r8
      END IF
      IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
!>      tl_A(Istr-1,Jstr  )=0.5_r8*(tl_A(Istr  ,Jstr  )+                &
!>   &                              tl_A(Istr-1,Jstr+1))
!>
        adfac=0.5_r8*ad_A(Istr-1,Jstr  )
        ad_A(Istr  ,Jstr  )=ad_A(Istr  ,Jstr  )+adfac
        ad_A(Istr-1,Jstr+1)=ad_A(Istr-1,Jstr+1)+adfac
        ad_A(Istr-1,Jstr  )=0.0_r8
      END IF
# endif

# ifndef NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Adjoint North-South boundary conditions: Closed or Gradient.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
        DO i=Istr,Iend
#  ifdef SOUTHERN_WALL
!>        tl_A(i,Jstr)=0.0_r8
!>
          ad_A(i,Jstr)=0.0_r8
#  else
!>        tl_A(i,Jstr)=tl_A(i,Jstr+1)
!>
          ad_A(i,Jstr+1)=ad_A(i,Jstr+1)+ad_A(i,Jstr)
          ad_A(i,Jstr)=0.0_r8
#  endif
        END DO
      END IF
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
        DO i=Istr,Iend
#  ifdef NORTHERN_WALL
!>        tl_A(i,Jend+1)=0.0_r8
!>
          ad_A(i,Jend+1)=0.0_r8
#  else
!>        tl_A(i,Jend+1)=tl_A(i,Jend)
!>
          ad_A(i,Jend)=ad_A(i,Jend)+ad_A(i,Jend+1)
          ad_A(i,Jend+1)=0.0_r8
#  endif
        END DO
      END IF
# endif

# ifndef EW_PERIODIC
!
!-----------------------------------------------------------------------
!  Adjoint East-West boundary conditions.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
#  ifdef WESTERN_WALL
#   ifdef NS_PERIODIC
#    define J_RANGE JstrV,Jend
#   else
#    define J_RANGE Jstr,JendR
#   endif
        DO j=J_RANGE
#   ifdef MASKING
!>        tl_A(Istr-1,j)=tl_A(Istr-1,j)*GRID(ng)%vmask(Istr-1,j)
!>
          ad_A(Istr-1,j)=ad_A(Istr-1,j)*GRID(ng)%vmask(Istr-1,j)
#   endif
!>        tl_A(Istr-1,j)=gamma2(ng)*tl_A(Istr,j)
!>
          ad_A(Istr  ,j)=ad_A(Istr,j)+gamma2(ng)*ad_A(Istr-1,j)
          ad_A(Istr-1,j)=0.0_r8
        END DO
#   undef J_RANGE
#  else
        DO j=JstrV,Jend
!>        tl_A(Istr-1,j)=tl_A(Istr,j)
!>
          ad_A(Istr  ,j)=ad_A(Istr,j)+ad_A(Istr-1,j)
          ad_A(Istr-1,j)=0.0_r8
        END DO
#  endif
      END IF

      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
#  ifdef EASTERN_WALL
#   ifdef NS_PERIODIC
#    define J_RANGE JstrV,Jend
#   else
#    define J_RANGE Jstr,JendR
#   endif
        DO j=J_RANGE
#   ifdef MASKING
!>        tl_A(Iend+1,j)=tl_A(Iend+1,j)*GRID(ng)%vmask(Iend+1,j)
!>
          ad_A(Iend+1,j)=ad_A(Iend+1,j)*GRID(ng)%vmask(Iend+1,j)
#   endif
!>        tl_A(Iend+1,j)=gamma2(ng)*tl_A(Iend,j)
!>
          ad_A(Iend  ,j)=ad_A(Iend,j)+ad_A(Iend+1,j)
          ad_A(Iend+1,j)=0.0_r8
        END DO
#   undef J_RANGE
#  else
        DO j=JstrV,Jend
!>        tl_A(Iend+1,j)=tl_A(Iend,j)
!>
          ad_A(Iend  ,j)=ad_A(Iend,j)+ad_A(Iend+1,j)
          ad_A(Iend+1,j)=0.0_r8
        END DO
#  endif
      END IF
# endif

      RETURN
      END SUBROUTINE ad_bc_v2d_tile

#endif
      END MODULE ad_bc_2d_mod
