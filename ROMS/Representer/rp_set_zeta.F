#include "cppdefs.h"
      MODULE rp_set_zeta_mod

#if defined TL_IOMS && defined SOLVE3D
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2007 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine sets tangent linear free-surface to its fast-time      !
!  averaged value.                                                     !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: rp_set_zeta

      CONTAINS
!
!***********************************************************************
      SUBROUTINE rp_set_zeta (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_coupling
      USE mod_ocean
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iRPM, 12)
# endif
      CALL rp_set_zeta_tile (ng, Istr, Iend, Jstr, Jend,                &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       COUPLING(ng) % tl_Zt_avg1,                 &
     &                       OCEAN(ng) % tl_zeta)
# ifdef PROFILE
      CALL wclock_off (ng, iRPM, 12)
# endif
      RETURN
      END SUBROUTINE rp_set_zeta
!
!***********************************************************************
      SUBROUTINE rp_set_zeta_tile (ng, Istr, Iend, Jstr, Jend,          &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             tl_Zt_avg1, tl_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars

# if defined EW_PERIODIC || defined NS_PERIODIC
!
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
# endif
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: tl_Zt_avg1(LBi:,LBj:)

      real(r8), intent(out) :: tl_zeta(LBi:,LBj:,:)
# else
      real(r8), intent(in) :: tl_Zt_avg1(LBi:UBi,LBj:UBj)

      real(r8), intent(out) :: tl_zeta(LBi:UBi,LBj:UBj,3)
# endif
!
!  Local variable declarations.
!
# ifdef DISTRIBUTE
#  ifdef EW_PERIODIC
      logical :: EWperiodic=.TRUE.
#  else
      logical :: EWperiodic=.FALSE.
#  endif
#  ifdef NS_PERIODIC
      logical :: NSperiodic=.TRUE.
#  else
      logical :: NSperiodic=.FALSE.
#  endif
# endif
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Prepare to time-step 2D equations:  set initial free-surface
!  to its fast-time averaged values (which corresponds to the time
!  step "n").
!-----------------------------------------------------------------------
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
!>        zeta(i,j,1)=Zt_avg1(i,j)
!>        zeta(i,j,2)=Zt_avg1(i,j)
!>
          tl_zeta(i,j,1)=tl_Zt_avg1(i,j)
          tl_zeta(i,j,2)=tl_Zt_avg1(i,j)
        END DO
      END DO
# if defined EW_PERIODIC || defined NS_PERIODIC
!>    CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
!>   &                        LBi, UBi, LBj, UBj,                       &
!>   &                        zeta(:,:,1))
!>
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        tl_zeta(:,:,1))
!>    CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
!>   &                        LBi, UBi, LBj, UBj,                       &
!>   &                        zeta(:,:,2))
!>
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        tl_zeta(:,:,2))
# endif
# ifdef DISTRIBUTE
!>    CALL mp_exchange2d (ng, iNLM, 2, Istr, Iend, Jstr, Jend,          &
!>   &                    LBi, UBi, LBj, UBj,                           &
!>   &                    NghostPoints, EWperiodic, NSperiodic,         &
!>   &                    zeta(:,:,1),                                  &
!>   &                    zeta(:,:,2))
!>
      CALL mp_exchange2d (ng, iNLM, 2, Istr, Iend, Jstr, Jend,          &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    tl_zeta(:,:,1),                               &
     &                    tl_zeta(:,:,2))
# endif

      RETURN
      END SUBROUTINE rp_set_zeta_tile
#endif
      END MODULE rp_set_zeta_mod