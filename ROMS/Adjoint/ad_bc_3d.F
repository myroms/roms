#include "cppdefs.h"
      MODULE ad_bc_3d_mod
#if defined ADJOINT && defined SOLVE3D
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2007 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This package applies gradient or periodic boundary conditions for   !
!  generic 3D adjoint fields.                                          !
!                                                                      !
!  Routines:                                                           !
!                                                                      !
!    ad_bc_r3d_tile     Boundary conditions for field at RHO-points    !
!    ad_bc_u3d_tile     Boundary conditions for field at U-points      !
!    ad_bc_v3d_tile     Boundary conditions for field at V-points      !
!    ad_bc_w3d_tile     Boundary conditions for field at W-points      !
!                                                                      !
!=======================================================================
!
      implicit none

      CONTAINS
!
!***********************************************************************
      SUBROUTINE ad_bc_r3d_tile (ng, Istr, Iend, Jstr, Jend,            &
     &                           LBi, UBi, LBj, UBj, LBk, UBk,          &
     &                           ad_A)
!***********************************************************************
!
      USE mod_param

# if defined EW_PERIODIC || defined NS_PERIODIC
!
      USE ad_exchange_3d_mod, ONLY : ad_exchange_r3d_tile
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: ad_A(LBi:,LBj:,LBk:)
# else
      real(r8), intent(inout) :: ad_A(LBi:UBi,LBj:UBj,LBk:UBk)
# endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j, k

      real(r8) :: adfac

# include "set_bounds.h"

# if defined EW_PERIODIC || defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Set periodic boundary conditons.
!-----------------------------------------------------------------------
!
      CALL ad_exchange_r3d_tile (ng, Istr, Iend, Jstr, Jend,            &
     &                           LBi, UBi, LBj, UBj, LBk, UBk,          &
     &                           ad_A)
# endif

# if !defined EW_PERIODIC && !defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF ((NORTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        DO k=1,N(ng)
!>        tl_A(Iend+1,Jend+1,k)=0.5_r8*(tl_A(Iend+1,Jend  ,k)+          &
!>   &                                  tl_A(Iend  ,Jend+1,k))
!>
          adfac=0.5_r8*ad_A(Iend+1,Jend+1,k)
          ad_A(Iend+1,Jend  ,k)=ad_A(Iend+1,Jend  ,k)+adfac
          ad_A(Iend  ,Jend+1,k)=ad_A(Iend  ,Jend+1,k)+adfac
          ad_A(Iend+1,Jend+1,k)=0.0_r8
        END DO
      END IF
      IF ((NORTHERN_EDGE).and.(WESTERN_EDGE)) THEN
        DO k=1,N(ng)
!>        tl_A(Istr-1,Jend+1,k)=0.5_r8*(tl_A(Istr-1,Jend  ,k)+          &
!>   &                                  tl_A(Istr  ,Jend+1,k))
!>
          adfac=0.5_r8*ad_A(Istr-1,Jend+1,k)
          ad_A(Istr-1,Jend  ,k)=ad_A(Istr-1,Jend  ,k)+adfac
          ad_A(Istr  ,Jend+1,k)=ad_A(Istr  ,Jend+1,k)+adfac
          ad_A(Istr-1,Jend+1,k)=0.0_r8
        END DO
      END IF
      IF ((SOUTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        DO k=1,N(ng)
!>        tl_A(Iend+1,Jstr-1,k)=0.5_r8*(tl_A(Iend  ,Jstr-1,k)+          &
!>   &                                  tl_A(Iend+1,Jstr  ,k))
!>
          adfac=0.5_r8*ad_A(Iend+1,Jstr-1,k)
          ad_A(Iend  ,Jstr-1,k)=ad_A(Iend  ,Jstr-1,k)+adfac
          ad_A(Iend+1,Jstr  ,k)=ad_A(Iend+1,Jstr  ,k)+adfac
          ad_A(Iend+1,Jstr-1,k)=0.0_r8
        END DO
      END IF
      IF ((SOUTHERN_EDGE).and.(WESTERN_EDGE)) THEN
        DO k=1,N(ng)
!>        tl_A(Istr-1,Jstr-1,k)=0.5_r8*(tl_A(Istr  ,Jstr-1,k)+          &
!>   &                                  tl_A(Istr-1,Jstr  ,k))
!>
          adfac=0.5_r8*ad_A(Istr-1,Jstr-1,k)
          ad_A(Istr  ,Jstr-1,k)=ad_A(Istr  ,Jstr-1,k)+adfac
          ad_A(Istr-1,Jstr  ,k)=ad_A(Istr-1,Jstr  ,k)+adfac
          ad_A(Istr-1,Jstr-1,k)=0.0_r8
        END DO
      END IF
# endif

# ifndef NS_PERIODIC
!
!-----------------------------------------------------------------------
!  North-South gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (SOUTHERN_EDGE) THEN
        DO k=1,N(ng)
          DO i=Istr,Iend
!>          tl_A(i,Jstr-1,k)=tl_A(i,Jstr,k)
!>
            ad_A(i,Jstr  ,k)=ad_A(i,Jstr,k)+ad_A(i,Jstr-1,k)
            ad_A(i,Jstr-1,k)=0.0_r8
          END DO
        END DO
      END IF
      IF (NORTHERN_EDGE) THEN
        DO k=1,N(ng)
          DO i=Istr,Iend
!>          tl_A(i,Jend+1,k)=tl_A(i,Jend,k)
!>
            ad_A(i,Jend  ,k)=ad_A(i,Jend,k)+ad_A(i,Jend+1,k)
            ad_A(i,Jend+1,k)=0.0_r8
          END DO
        END DO
      END IF
# endif

# ifndef EW_PERIODIC
!
!-----------------------------------------------------------------------
!  East-West gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (WESTERN_EDGE) THEN
        DO k=1,N(ng)
          DO j=Jstr,Jend
!>          tl_A(Istr-1,j,k)=tl_A(Istr,j,k)
!>
            ad_A(Istr  ,j,k)=ad_A(Istr,j,k)+ad_A(Istr-1,j,k)
            ad_A(Istr-1,j,k)=0.0_r8
          END DO
        END DO
      END IF
      IF (EASTERN_EDGE) THEN
        DO k=1,N(ng)
          DO j=Jstr,Jend
!>          tl_A(Iend+1,j,k)=tl_A(Iend,j,k)
!>
            ad_A(Iend  ,j,k)=ad_A(Iend,j,k)+ad_A(Iend+1,j,k)
            ad_A(Iend+1,j,k)=0.0_r8
          END DO
        END DO
      END IF
# endif

      RETURN
      END SUBROUTINE ad_bc_r3d_tile

!
!***********************************************************************
      SUBROUTINE ad_bc_u3d_tile (ng, Istr, Iend, Jstr, Jend,            &
     &                           LBi, UBi, LBj, UBj, LBk, UBk,          &
     &                           ad_A)
!***********************************************************************
!
      USE mod_param

# if defined EW_PERIODIC || defined NS_PERIODIC
!
      USE ad_exchange_3d_mod, ONLY : ad_exchange_u3d_tile
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: ad_A(LBi:,LBj:,LBk:)
# else
      real(r8), intent(inout) :: ad_A(LBi:UBi,LBj:UBj,LBk:UBk)
# endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j, k

      real(r8) :: adfac

# include "set_bounds.h"

# if defined EW_PERIODIC || defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Set periodic boundary conditons.
!-----------------------------------------------------------------------
!
      CALL ad_exchange_u3d_tile (ng, Istr, Iend, Jstr, Jend,            &
     &                           LBi, UBi, LBj, UBj, LBk, UBk,          &
     &                           ad_A)
# endif

# if !defined EW_PERIODIC && !defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF ((NORTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        DO k=1,N(ng)
!>        tl_A(Iend+1,Jend+1,k)=0.5_r8*(tl_A(Iend+1,Jend  ,k)+          &
!>   &                                  tl_A(Iend  ,Jend+1,k))
!>
          adfac=0.5_r8*ad_A(Iend+1,Jend+1,k)
          ad_A(Iend+1,Jend  ,k)=ad_A(Iend+1,Jend  ,k)+adfac
          ad_A(Iend  ,Jend+1,k)=ad_A(Iend  ,Jend+1,k)+adfac
          ad_A(Iend+1,Jend+1,k)=0.0_r8
        END DO
      END IF
      IF ((NORTHERN_EDGE).and.(WESTERN_EDGE)) THEN
        DO k=1,N(ng)
!>        tl_A(Istr  ,Jend+1,k)=0.5_r8*(tl_A(Istr  ,Jend  ,k)+          &
!>   &                                  tl_A(Istr+1,Jend+1,k))
!>
          adfac=0.5_r8*ad_A(Istr  ,Jend+1,k)
          ad_A(Istr  ,Jend  ,k)=ad_A(Istr  ,Jend  ,k)+adfac
          ad_A(Istr+1,Jend+1,k)=ad_A(Istr+1,Jend+1,k)+adfac
          ad_A(Istr  ,Jend+1,k)=0.0_r8
        END DO
      END IF
      IF ((SOUTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        DO k=1,N(ng)
!>        tl_A(Iend+1,Jstr-1,k)=0.5_r8*(tl_A(Iend  ,Jstr-1,k)+          &
!>   &                                  tl_A(Iend+1,Jstr  ,k))
!>
          adfac=0.5_r8*ad_A(Iend+1,Jstr-1,k)
          ad_A(Iend  ,Jstr-1,k)=ad_A(Iend  ,Jstr-1,k)+adfac
          ad_A(Iend+1,Jstr  ,k)=ad_A(Iend+1,Jstr  ,k)+adfac
          ad_A(Iend+1,Jstr-1,k)=0.0_r8
        END DO
      END IF
      IF ((SOUTHERN_EDGE).and.(WESTERN_EDGE)) THEN
        DO k=1,N(ng)
!>        tl_A(Istr  ,Jstr-1,k)=0.5_r8*(tl_A(Istr+1,Jstr-1,k)+          &
!>   &                                  tl_A(Istr  ,Jstr  ,k))
!>
          adfac=0.5_r8*ad_A(Istr  ,Jstr-1,k)
          ad_A(Istr+1,Jstr-1,k)=ad_A(Istr+1,Jstr-1,k)+adfac
          ad_A(Istr  ,Jstr  ,k)=ad_A(Istr  ,Jstr  ,k)+adfac
          ad_A(Istr  ,Jstr-1,k)=0.0_r8
        END DO
      END IF
# endif

# ifndef NS_PERIODIC
!
!-----------------------------------------------------------------------
!  North-South gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (SOUTHERN_EDGE) THEN
        DO k=1,N(ng)
          DO i=IstrU,Iend
!>          tl_A(i,Jstr-1,k)=tl_A(i,Jstr,k)
!>
            ad_A(i,Jstr  ,k)=ad_A(i,Jstr,k)+ad_A(i,Jstr-1,k)
            ad_A(i,Jstr-1,k)=0.0_r8
          END DO
        END DO
      END IF
      IF (NORTHERN_EDGE) THEN
        DO k=1,N(ng)
          DO i=IstrU,Iend
!>          tl_A(i,Jend+1,k)=tl_A(i,Jend,k)
!>
            ad_A(i,Jend  ,k)=ad_A(i,Jend,k)+ad_A(i,Jend+1,k)
            ad_A(i,Jend+1,k)=0.0_r8
          END DO
        END DO
      END IF
# endif

# ifndef EW_PERIODIC
!
!-----------------------------------------------------------------------
!  East-West gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (WESTERN_EDGE) THEN
        DO k=1,N(ng)
          DO j=Jstr,Jend
!>          tl_A(Istr,j,k)=tl_A(Istr+1,j,k)
!>
            ad_A(Istr+1,j,k)=ad_A(Istr+1,j,k)+ad_A(Istr,j,k)
            ad_A(Istr  ,j,k)=0.0_r8
          END DO
        END DO
      END IF
      IF (EASTERN_EDGE) THEN
        DO k=1,N(ng)
          DO j=Jstr,Jend
!>          tl_A(Iend+1,j,k)=tl_A(Iend,j,k)
!>
            ad_A(Iend  ,j,k)=ad_A(Iend,j,k)+ad_A(Iend+1,j,k)
            ad_A(Iend+1,j,k)=0.0_r8
          END DO
        END DO
      END IF
# endif

      RETURN
      END SUBROUTINE ad_bc_u3d_tile

!
!***********************************************************************
      SUBROUTINE ad_bc_v3d_tile (ng, Istr, Iend, Jstr, Jend,            &
     &                           LBi, UBi, LBj, UBj, LBk, UBk,          &
     &                           ad_A)
!***********************************************************************
!
      USE mod_param

# if defined EW_PERIODIC || defined NS_PERIODIC
!
      USE ad_exchange_3d_mod, ONLY : ad_exchange_v3d_tile
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: ad_A(LBi:,LBj:,LBk:)
# else
      real(r8), intent(inout) :: ad_A(LBi:UBi,LBj:UBj,LBk:UBk)
# endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j, k

      real(r8) :: adfac

# include "set_bounds.h"

# if defined EW_PERIODIC || defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Set periodic boundary conditons.
!-----------------------------------------------------------------------
!
      CALL ad_exchange_v3d_tile (ng, Istr, Iend, Jstr, Jend,            &
     &                           LBi, UBi, LBj, UBj, LBk, UBk,          &
     &                           ad_A)
# endif

# if !defined EW_PERIODIC && !defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF ((NORTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        DO k=1,N(ng)
!>        tl_A(Iend+1,Jend+1,k)=0.5_r8*(tl_A(Iend+1,Jend  ,k)+          &
!>   &                                  tl_A(Iend  ,Jend+1,k))
!>
          adfac=0.5_r8*ad_A(Iend+1,Jend+1,k)
          ad_A(Iend+1,Jend  ,k)=ad_A(Iend+1,Jend  ,k)+adfac
          ad_A(Iend  ,Jend+1,k)=ad_A(Iend  ,Jend+1,k)+adfac
          ad_A(Iend+1,Jend+1,k)=0.0_r8
        END DO
      END IF
      IF ((NORTHERN_EDGE).and.(WESTERN_EDGE)) THEN
        DO k=1,N(ng)
!>        tl_A(Istr-1,Jend+1,k)=0.5_r8*(tl_A(Istr-1,Jend  ,k)+          &
!>   &                                  tl_A(Istr  ,Jend+1,k))
!>
          adfac=0.5_r8*ad_A(Istr-1,Jend+1,k)
          ad_A(Istr-1,Jend  ,k)=ad_A(Istr-1,Jend  ,k)+adfac
          ad_A(Istr  ,Jend+1,k)=ad_A(Istr  ,Jend+1,k)+adfac
          ad_A(Istr-1,Jend+1,k)=0.0_r8
        END DO
      END IF
      IF ((SOUTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        DO k=1,N(ng)
!>        tl_A(Iend+1,Jstr  ,k)=0.5_r8*(tl_A(Iend  ,Jstr  ,k)+          &
!>   &                                  tl_A(Iend+1,Jstr+1,k))
!>
          adfac=0.5_r8*ad_A(Iend+1,Jstr  ,k)
          ad_A(Iend  ,Jstr  ,k)=ad_A(Iend  ,Jstr  ,k)+adfac
          ad_A(Iend+1,Jstr+1,k)=ad_A(Iend+1,Jstr+1,k)+adfac
          ad_A(Iend+1,Jstr  ,k)=0.0_r8         
        END DO
      END IF
      IF ((SOUTHERN_EDGE).and.(WESTERN_EDGE)) THEN
        DO k=1,N(ng)
!>        tl_A(Istr-1,Jstr  ,k)=0.5_r8*(tl_A(Istr  ,Jstr  ,k)+          &
!>   &                                  tl_A(Istr-1,Jstr+1,k))
!>
          adfac=0.5_r8*ad_A(Istr-1,Jstr  ,k)
          ad_A(Istr  ,Jstr  ,k)=ad_A(Istr  ,Jstr  ,k)+adfac
          ad_A(Istr-1,Jstr+1,k)=ad_A(Istr-1,Jstr+1,k)+adfac
          ad_A(Istr-1,Jstr  ,k)=0.0_r8
        END DO
      END IF
# endif

# ifndef NS_PERIODIC
!
!-----------------------------------------------------------------------
!  North-South periodic boundary conditions.
!-----------------------------------------------------------------------
!
      IF (SOUTHERN_EDGE) THEN
        DO k=1,N(ng)
          DO i=Istr,Iend
!>          tl_A(i,Jstr,k)=tl_A(i,Jstr+1,k)
!>
            ad_A(i,Jstr+1,k)=ad_A(i,Jstr+1,k)+ad_A(i,Jstr,k)
            ad_A(i,Jstr  ,k)=0.0_r8
          END DO
        END DO
      END IF
      IF (NORTHERN_EDGE) THEN
        DO k=1,N(ng)
          DO i=Istr,Iend
!>          tl_A(i,Jend+1,k)=tl_A(i,Jend,k)
!>
            ad_A(i,Jend  ,k)=ad_A(i,Jend,k)+ad_A(i,Jend+1,k)
            ad_A(i,Jend+1,k)=0.0_r8
          END DO
        END DO
      END IF
# endif

# ifndef EW_PERIODIC
!
!-----------------------------------------------------------------------
!  East-West gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (WESTERN_EDGE) THEN
        DO k=1,N(ng)
          DO j=JstrV,Jend
!>          tl_A(Istr-1,j,k)=tl_A(Istr,j,k)
!>
            ad_A(Istr  ,j,k)=ad_A(Istr,j,k)+ad_A(Istr-1,j,k)
            ad_A(Istr-1,j,k)=0.0_r8
          END DO
        END DO
      END IF
      IF (EASTERN_EDGE) THEN
        DO k=1,N(ng)
          DO j=JstrV,Jend
!>          tl_A(Iend+1,j,k)=tl_A(Iend,j,k)
!>
            ad_A(Iend  ,j,k)=ad_A(Iend,j,k)+ad_A(Iend+1,j,k)
            ad_A(Iend+1,j,k)=0.0_r8
          END DO
        END DO
      END IF
# endif

      RETURN
      END SUBROUTINE ad_bc_v3d_tile

!
!***********************************************************************
      SUBROUTINE ad_bc_w3d_tile (ng, Istr, Iend, Jstr, Jend,            &
     &                           LBi, UBi, LBj, UBj, LBk, UBk,          &
     &                           ad_A)
!***********************************************************************
!
      USE mod_param

# if defined EW_PERIODIC || defined NS_PERIODIC
!
      USE ad_exchange_3d_mod, ONLY : ad_exchange_w3d_tile
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: ad_A(LBi:,LBj:,LBk:)
# else
      real(r8), intent(inout) :: ad_A(LBi:UBi,LBj:UBj,LBk:UBk)
# endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j, k

      real(r8) :: adfac

# include "set_bounds.h"

# if defined EW_PERIODIC || defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Set periodic boundary conditons.
!-----------------------------------------------------------------------
!
      CALL ad_exchange_w3d_tile (ng, Istr, Iend, Jstr, Jend,            &
     &                           LBi, UBi, LBj, UBj, LBk, UBk,          &
     &                           ad_A)
# endif

# if !defined EW_PERIODIC && !defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF ((NORTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        DO k=0,N(ng)
!>        tl_A(Iend+1,Jend+1,k)=0.5_r8*(tl_A(Iend+1,Jend  ,k)+          &
!>   &                                  tl_A(Iend  ,Jend+1,k))
!>
          adfac=0.5_r8*ad_A(Iend+1,Jend+1,k)
          ad_A(Iend+1,Jend  ,k)=ad_A(Iend+1,Jend  ,k)+adfac
          ad_A(Iend  ,Jend+1,k)=ad_A(Iend  ,Jend+1,k)+adfac
          ad_A(Iend+1,Jend+1,k)=0.0_r8
        END DO
      END IF
      IF ((NORTHERN_EDGE).and.(WESTERN_EDGE)) THEN
        DO k=0,N(ng)
!>        tl_A(Istr-1,Jend+1,k)=0.5_r8*(tl_A(Istr-1,Jend  ,k)+          &
!>   &                                  tl_A(Istr  ,Jend+1,k))
!>
          adfac=0.5_r8*ad_A(Istr-1,Jend+1,k)
          ad_A(Istr-1,Jend  ,k)=ad_A(Istr-1,Jend  ,k)+adfac
          ad_A(Istr  ,Jend+1,k)=ad_A(Istr  ,Jend+1,k)+adfac
          ad_A(Istr-1,Jend+1,k)=0.0_r8
        END DO
      END IF
      IF ((SOUTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        DO k=0,N(ng)
!>        tl_A(Iend+1,Jstr-1,k)=0.5_r8*(tl_A(Iend+1,Jstr  ,k)+          &
!>                                      tl_A(Iend  ,Jstr-1,k))
!>
          adfac=0.5_r8*ad_A(Iend+1,Jstr-1,k)
          ad_A(Iend+1,Jstr  ,k)=ad_A(Iend+1,Jstr  ,k)+adfac
          ad_A(Iend  ,Jstr-1,k)=ad_A(Iend  ,Jstr-1,k)+adfac
          ad_A(Iend+1,Jstr-1,k)=0.0_r8
        END DO
      END IF
      IF ((SOUTHERN_EDGE).and.(WESTERN_EDGE)) THEN
        DO k=0,N(ng)
!>        tl_A(Istr-1,Jstr-1,k)=0.5_r8*(tl_A(Istr  ,Jstr-1,k)+          &
!>   &                                  tl_A(Istr-1,Jstr  ,k))
!>
          adfac=0.5_r8*ad_A(Istr-1,Jstr-1,k)
          ad_A(Istr  ,Jstr-1,k)=ad_A(Istr  ,Jstr-1,k)+adfac
          ad_A(Istr-1,Jstr  ,k)=ad_A(Istr-1,Jstr  ,k)+adfac
          ad_A(Istr-1,Jstr-1,k)=0.0_r8
        END DO
      END IF
# endif

# ifndef NS_PERIODIC
!
!-----------------------------------------------------------------------
!  North-South gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (SOUTHERN_EDGE) THEN
        DO k=0,N(ng)
          DO i=Istr,Iend
!>          tl_A(i,Jstr-1,k)=tl_A(i,Jstr,k)
!>
            ad_A(i,Jstr  ,k)=ad_A(i,Jstr,k)+ad_A(i,Jstr-1,k)
            ad_A(i,Jstr-1,k)=0.0_r8
          END DO
        END DO
      END IF
      IF (NORTHERN_EDGE) THEN
        DO k=0,N(ng)
          DO i=Istr,Iend
!>          tl_A(i,Jend+1,k)=tl_A(i,Jend,k)
!>
            ad_A(i,Jend  ,k)=ad_A(i,Jend,k)+ad_A(i,Jend+1,k)
            ad_A(i,Jend+1,k)=0.0_r8
          END DO
        END DO
      END IF
# endif

# ifndef EW_PERIODIC
!
!-----------------------------------------------------------------------
!  East-West gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (WESTERN_EDGE) THEN
        DO k=0,N(ng)
          DO j=Jstr,Jend
!>          tl_A(Istr-1,j,k)=tl_A(Istr,j,k)
!>
            ad_A(Istr  ,j,k)=ad_A(Istr,j,k)+ad_A(Istr-1,j,k)
            ad_A(Istr-1,j,k)=0.0_r8
          END DO
        END DO
      END IF
      IF (EASTERN_EDGE) THEN
        DO k=0,N(ng)
          DO j=Jstr,Jend
!>          tl_A(Iend+1,j,k)=tl_A(Iend,j,k)
!>
            ad_A(Iend  ,j,k)=ad_A(Iend,j,k)+ad_A(Iend+1,j,k)
            ad_A(Iend+1,j,k)=0.0_r8
          END DO
        END DO
      END IF
# endif
      RETURN
      END SUBROUTINE ad_bc_w3d_tile

#endif
      END MODULE ad_bc_3d_mod
