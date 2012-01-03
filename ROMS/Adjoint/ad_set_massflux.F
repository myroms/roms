#include "cppdefs.h"
      MODULE ad_set_massflux_mod

#if defined ADJOINT && defined SOLVE3D
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2012 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes adjoint horizontal mass fluxes, Hz*u/n and    !
!  Hz*v/m.                                                             !
!                                                                      !
!   BASIC STATE variables required:  Hz, u, v                          !
!   Dependend variables:  ad_Huon, ad_Hvom                             !
!   Independend variables:  ad_Hz, ad_u, ad_v                          !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC  :: ad_set_massflux
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE ad_set_massflux (ng, tile)
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
      CALL wclock_on (ng, iADM, 12)
# endif
      CALL ad_set_massflux_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           nrhs(ng),                              &
     &                           OCEAN(ng) % u,                         &
     &                           OCEAN(ng) % v,                         &
     &                           OCEAN(ng) % ad_u,                      &
     &                           OCEAN(ng) % ad_v,                      &
# ifdef NEARSHORE_MELLOR
     &                           OCEAN(ng) % u_stokes,                  &
     &                           OCEAN(ng) % v_stokes,                  &
     &                           OCEAN(ng) % ad_u_stokes,               &
     &                           OCEAN(ng) % ad_v_stokes,               &
# endif
     &                           GRID(ng) % Hz,                         &
     &                           GRID(ng) % ad_Hz,                      &
     &                           GRID(ng) % om_v,                       &
     &                           GRID(ng) % on_u,                       &
     &                           GRID(ng) % ad_Huon,                    &
     &                           GRID(ng) % ad_Hvom)
# ifdef PROFILE
      CALL wclock_off (ng, iADM, 12)
# endif
      RETURN
      END SUBROUTINE ad_set_massflux

!
!***********************************************************************
      SUBROUTINE ad_set_massflux_tile (ng, tile,                        &
     &                                 LBi, UBi, LBj, UBj,              &
     &                                 IminS, ImaxS, JminS, JmaxS,      &
     &                                 nrhs,                            &
     &                                 u, v,                            &
     &                                 ad_u, ad_v,                      &
# ifdef NEARSHORE_MELLOR
     &                                 u_stokes, v_stokes,              &
     &                                 ad_u_stokes, ad_v_stokes,        &
# endif
     &                                 Hz, ad_Hz,                       &
     &                                 om_v, on_u,                      &
     &                                 ad_Huon, ad_Hvom)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE ad_exchange_3d_mod
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : ad_mp_exchange3d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
#  ifdef NEARSHORE_MELLOR
      real(r8), intent(in) :: u_stokes(LBi:,LBj:,:)
      real(r8), intent(in) :: v_stokes(LBi:,LBj:,:)
#  endif
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)

      real(r8), intent(inout) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_v(LBi:,LBj:,:,:)
#  ifdef NEARSHORE_MELLOR
      real(r8), intent(inout) :: ad_u_stokes(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_v_stokes(LBi:,LBj:,:)
#  endif
      real(r8), intent(inout) :: ad_Hz(LBi:,LBj:,:)

      real(r8), intent(inout) :: ad_Huon(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_Hvom(LBi:,LBj:,:)
# else
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),2)
#  ifdef NEARSHORE_MELLOR
      real(r8), intent(in) :: u_stokes(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: v_stokes(LBi:UBi,LBj:UBj,N(ng))
#  endif
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: om_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)

      real(r8), intent(inout) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
#  ifdef NEARSHORE_MELLOR
      real(r8), intent(inout) :: ad_u_stokes(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: ad_v_stokes(LBi:UBi,LBj:UBj,N(ng))
#  endif
      real(r8), intent(inout) :: ad_Hz(LBi:UBi,LBj:UBj,N(ng))

      real(r8), intent(inout) :: ad_Huon(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: ad_Hvom(LBi:UBi,LBj:UBj,N(ng))
# endif
!
!  Local variable declarations.
!
      integer :: IminU, IminV, ImaxU, ImaxV
      integer :: JminU, JminV, JmaxU, JmaxV
      integer :: i, j, k

      real(r8) :: adfac, adfac1

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute horizontal mass fluxes, Hz*u/n and Hz*v/m.
!-----------------------------------------------------------------------
!
!  Set local I- and J-ranges.
!
      IF (EWperiodic(ng)) THEN
        IminU=Istr
        ImaxU=Iend
        IminV=Istr
        ImaxV=Iend
      ELSE
        IF (ComposedGrid(ng).or.RefinedGrid(ng)) THEN
          IminU=IstrP+1
          ImaxU=IendT
          IminV=IstrT
          ImaxV=IendT
        ELSE
          IminU=Istr
          ImaxU=IendR
          IminV=IstrR
          ImaxV=IendR
        END IF
      END IF

      IF (NSperiodic(ng)) THEN
        JminU=Jstr
        JmaxU=Jend
        JminV=Jstr
        JmaxV=Jend
      ELSE
        IF (ComposedGrid(ng).or.RefinedGrid(ng)) THEN
          JminU=JstrT
          JmaxU=JendT
          JminV=JstrP+1
          JmaxV=JendT
        ELSE
          JminU=JstrR
          JmaxU=JendR
          JminV=Jstr
          JmaxV=JendR
        END IF
      END IF
!
!  Exchange boundary information.
!
# ifdef DISTRIBUTE
!>    CALL mp_exchange3d (ng, tile, iTLM, 2,                            &
!>   &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
!>   &                    NghostPoints,                                 &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    tl_Huon, tl_Hvom)
!>
      CALL ad_mp_exchange3d (ng, tile, iADM, 2,                         &
     &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_Huon, ad_Hvom)
!
# endif

      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
!>      CALL exchange_v3d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
!>   &                          tl_Hvom)
!>
        CALL ad_exchange_v3d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj, 1, N(ng),        &
     &                             ad_Hvom)
!>      CALL exchange_u3d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
!>   &                          tl_Huon)
!>
        CALL ad_exchange_u3d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj, 1, N(ng),        &
     &                             ad_Huon)
      END IF
!
!  Compute adjoint horizontal mass fluxes.
!
      DO k=1,N(ng)
        DO j=JminV,JmaxV
          DO i=IminV,ImaxV
# ifdef NEARSHORE_MELLOR
!>          tl_Hvom(i,j,k)=tl_Hvom(i,j,k)+                              &
!>   &                     0.5_r8*om_v(i,j)*                            &
!>   &                     ((Hz(i,j,k)+Hz(i,j-1,k))*                    &
!>   &                      tl_v_stokes(i,j,k)+                         &
!>   &                      (tl_Hz(i,j,k)+tl_Hz(i,j-1,k))*              &
!>   &                      v_stokes(i,j,k))
!>
            adfac=0.5_r8*om_v(i,j)*tl_Hvom(i,j,k)
            adfac1=adfac*v_stokes(i,j,k)
            tl_v_stokes(i,j,k)=tl_v_stokes(i,j,k)+                      &
     &                         adfac*(Hz(i,j,k)+Hz(i,j-1,k))
            ad_Hz(i,j-1,k)=ad_Hz(i,j-1,k)+adfac1
            ad_Hz(i,j  ,k)=ad_Hz(i,j  ,k)+adfac1
# endif
!>          tl_Hvom(i,j,k)=0.5_r8*om_v(i,j)*                            &
!>   &                     ((Hz(i,j,k)+Hz(i,j-1,k))*                    &
!>   &                      tl_v(i,j,k,nrhs)+                           &
!>   &                      (tl_Hz(i,j,k)+tl_Hz(i,j-1,k))*              &
!>   &                      v(i,j,k,nrhs))
!>
            adfac=0.5_r8*om_v(i,j)*ad_Hvom(i,j,k)
            adfac1=adfac*v(i,j,k,nrhs)
            ad_v(i,j,k,nrhs)=ad_v(i,j,k,nrhs)+                          &
     &                       adfac*(Hz(i,j,k)+Hz(i,j-1,k))
            ad_Hz(i,j-1,k)=ad_Hz(i,j-1,k)+adfac1
            ad_Hz(i,j  ,k)=ad_Hz(i,j  ,k)+adfac1
            ad_Hvom(i,j,k)=0.0_r8
          END DO
        END DO
        DO j=JminU,JmaxU
          DO i=IminU,ImaxU
# ifdef NEARSHORE_MELLOR
!>          tl_Huon(i,j,k)=tl_Huon(i,j,k)+                              &
!>   &                     0.5_r8*on_u(i,j)*                            &
!>   &                     ((Hz(i,j,k)+Hz(i-1,j,k))*                    &
!>   &                      tl_u_stokes(i,j,k)+                         &
!>   &                      (tl_Hz(i,j,k)+tl_Hz(i-1,j,k))*              &
!>   &                      u_stokes(i,j,k))
!>
            adfac=0.5_r8*on_u(i,j)*ad_Huon(i,j,k)
            adfac1=adfac*u_stokes(i,j,k)
            ad_u_stokes(i,j,k)=ad_u_stokes(i,j,k)+                      &
     &                         adfac*(Hz(i,j,k)+Hz(i-1,j,k))
            ad_Hz(i-1,j,k)=ad_Hz(i-1,j,k)+adfac1
            ad_Hz(i  ,j,k)=ad_Hz(i  ,j,k)+adfac1
# endif
!>          tl_Huon(i,j,k)=0.5_r8*on_u(i,j)*                            &
!>   &                     ((Hz(i,j,k)+Hz(i-1,j,k))*                    &
!>   &                      tl_u(i,j,k,nrhs)+                           &
!>   &                      (tl_Hz(i,j,k)+tl_Hz(i-1,j,k))*              &
!>   &                      u(i,j,k,nrhs))
!>
            adfac=0.5_r8*on_u(i,j)*ad_Huon(i,j,k)
            adfac1=adfac*u(i,j,k,nrhs)
            ad_u(i,j,k,nrhs)=ad_u(i,j,k,nrhs)+                          &
     &                       adfac*(Hz(i,j,k)+Hz(i-1,j,k))
            ad_Hz(i-1,j,k)=ad_Hz(i-1,j,k)+adfac1
            ad_Hz(i  ,j,k)=ad_Hz(i  ,j,k)+adfac1
            ad_Huon(i,j,k)=0.0_r8
          END DO
        END DO
      END DO

      RETURN
      END SUBROUTINE ad_set_massflux_tile
#endif
      END MODULE ad_set_massflux_mod
