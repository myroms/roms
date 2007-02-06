#include "cppdefs.h"
      MODULE rp_set_massflux_mod

#if defined TL_IOMS && defined SOLVE3D

# ifdef EW_PERIODIC
#  define IU_RANGE Istr,Iend
#  define IV_RANGE Istr,Iend
# else
#  define IU_RANGE Istr,IendR
#  define IV_RANGE IstrR,IendR
# endif
# ifdef NS_PERIODIC
#  define JU_RANGE Jstr,Jend
#  define JV_RANGE Jstr,Jend
# else
#  define JU_RANGE JstrR,JendR
#  define JV_RANGE Jstr,JendR
# endif
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2007 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes representers tangent linear horizontal mass   !
!  fluxes, Hz*u/n and Hz*v/m.                                          !
!                                                                      !
!   BASIC STATE variables required:  Hz, u, v                          !
!   Dependend variables:  tl_Huon, tl_Hvom                             !
!   Independend variables:  tl_Hz, tl_u, tl_v                          !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: rp_set_massflux

      CONTAINS
!
!***********************************************************************
      SUBROUTINE rp_set_massflux (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_ocean
      USE mod_stepping
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
      CALL rp_set_massflux_tile (ng, Istr, Iend, Jstr, Jend,            &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           nrhs(ng),                              &
     &                           OCEAN(ng) % u,                         &
     &                           OCEAN(ng) % v,                         &
     &                           OCEAN(ng) % tl_u,                      &
     &                           OCEAN(ng) % tl_v,                      &
     &                           GRID(ng) % Hz,                         &
     &                           GRID(ng) % tl_Hz,                      &
     &                           GRID(ng) % om_v,                       &
     &                           GRID(ng) % on_u,                       &
     &                           GRID(ng) % tl_Huon,                    &
     &                           GRID(ng) % tl_Hvom)
# ifdef PROFILE
      CALL wclock_off (ng, iRPM, 12)
# endif
      RETURN
      END SUBROUTINE rp_set_massflux
!
!***********************************************************************
      SUBROUTINE rp_set_massflux_tile (ng, Istr, Iend, Jstr, Jend,      &
     &                                 LBi, UBi, LBj, UBj,              &
     &                                 nrhs,                            &
     &                                 u, v,                            &
     &                                 tl_u, tl_v,                      &
     &                                 Hz, tl_Hz,                       &
     &                                 om_v, on_u,                      &
     &                                 tl_Huon, tl_Hvom)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars

# if defined EW_PERIODIC || defined NS_PERIODIC
!
      USE exchange_3d_mod
# endif
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange3d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: nrhs
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
      real(r8), intent(in) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: tl_v(LBi:,LBj:,:,:)
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: tl_Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)

      real(r8), intent(out) :: tl_Huon(LBi:,LBj:,:)
      real(r8), intent(out) :: tl_Hvom(LBi:,LBj:,:)
# else
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: tl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: tl_v(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: tl_Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: om_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)

      real(r8), intent(out) :: tl_Huon(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(out) :: tl_Hvom(LBi:UBi,LBj:UBj,N(ng))
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
      integer :: i, j, k

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute horizontal mass fluxes, Hz*u/n and Hz*v/m.
!-----------------------------------------------------------------------
!
      DO k=1,N(ng)
        DO j=JU_RANGE
          DO i=IU_RANGE
!>          Huon(i,j,k)=0.5_r8*(Hz(i,j,k)+Hz(i-1,j,k))*u(i,j,k,nrhs)*   &
!>   &                  on_u(i,j)
!>
            tl_Huon(i,j,k)=0.5_r8*on_u(i,j)*                            &
     &                     ((Hz(i,j,k)+Hz(i-1,j,k))*                    &
     &                      tl_u(i,j,k,nrhs)+                           &
     &                      (tl_Hz(i,j,k)+tl_Hz(i-1,j,k))*              &
     &                      u(i,j,k,nrhs))-                             &
# ifdef TL_IOMS
     &                     0.5_r8*on_u(i,j)*                            &
     &                     (Hz(i,j,k)+Hz(i-1,j,k))*u(i,j,k,nrhs)
# endif
          END DO
        END DO
        DO j=JV_RANGE
          DO i=IV_RANGE
!>          Hvom(i,j,k)=0.5_r8*(Hz(i,j,k)+Hz(i,j-1,k))*v(i,j,k,nrhs)*   &
!>   &                  om_v(i,j)
!>
            tl_Hvom(i,j,k)=0.5_r8*om_v(i,j)*                            &
     &                     ((Hz(i,j,k)+Hz(i,j-1,k))*                    &
     &                      tl_v(i,j,k,nrhs)+                           &
     &                      (tl_Hz(i,j,k)+tl_Hz(i,j-1,k))*              &
     &                      v(i,j,k,nrhs))-                             &
#  ifdef TL_IOMS
     &                     0.5_r8*om_v(i,j)*                            &
     &                     (Hz(i,j,k)+Hz(i,j-1,k))*v(i,j,k,nrhs)
#  endif
          END DO
        END DO
      END DO
# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
!
!  Exchange boundary information.
!
#  if defined EW_PERIODIC || defined NS_PERIODIC
!>    CALL exchange_u3d_tile (ng, Istr, Iend, Jstr, Jend,               &
!>   &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
!>   &                        Huon)
!>
      CALL exchange_u3d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        tl_Huon)
!>    CALL exchange_v3d_tile (ng, Istr, Iend, Jstr, Jend,               &
!>   &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
!>   &                        Hvom)
!>
      CALL exchange_v3d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        tl_Hvom)
#  endif
#  ifdef DISTRIBUTE
!>    CALL mp_exchange3d (ng, iNLM, 2, Istr, Iend, Jstr, Jend,          &
!>   &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
!>   &                    NghostPoints, EWperiodic, NSperiodic,         &
!>   &                    Huon, Hvom)
!>
      CALL mp_exchange3d (ng, iRPM, 2, Istr, Iend, Jstr, Jend,          &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    tl_Huon, tl_Hvom)
#  endif
# endif
      RETURN
      END SUBROUTINE rp_set_massflux_tile
#endif
      END MODULE rp_set_massflux_mod

