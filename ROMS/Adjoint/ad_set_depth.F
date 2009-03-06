#include "cppdefs.h"
      MODULE ad_set_depth_mod
#if defined ADJOINT && defined SOLVE3D
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2009 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes the time evolving depths of the model grid    !
!  and its associated vertical transformation metric (thickness).      !
!                                                                      !
!  Currently, two vertical coordinate transformations are available    !
!  with various possible vertical stretching, C(s), functions, (see    !
!  routine "set_scoord.F" for details).                                !
!                                                                      !
!  BASIC STATE variables needed: NONE                                  !
!  Independent Variables: ad_Hz, ad_z_r, ad_z_w                        !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: ad_set_depth, ad_set_depth_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE ad_set_depth (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_coupling
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
      CALL ad_set_depth_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        nstp(ng), nnew(ng),                       &
     &                        GRID(ng) % h,                             &
     &                        GRID(ng) % ad_h,                          &
# ifdef ICESHELF
     &                        GRID(ng) % zice,                          &
# endif
# if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
     &                        OCEAN(ng) % ad_bed,                       &
     &                        GRID(ng) % ad_bed_thick0,                 &
# endif
     &                        COUPLING(ng) % Zt_avg1,                   &
     &                        COUPLING(ng) % ad_Zt_avg1,                &
     &                        GRID(ng) % ad_Hz,                         &
     &                        GRID(ng) % ad_z_r,                        &
     &                        GRID(ng) % ad_z_w)
      RETURN
      END SUBROUTINE ad_set_depth

!
!***********************************************************************
      SUBROUTINE ad_set_depth_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              IminS, ImaxS, JminS, JmaxS,         &
     &                              nstp, nnew,                         &
     &                              h, ad_h,                            &
# ifdef ICESHELF
     &                              zice,                               &
# endif
# if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
     &                              ad_bed, bed_thick0,                 &
# endif
     &                              Zt_avg1, ad_Zt_avg1,                &
     &                              ad_Hz, ad_z_r, ad_z_w)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
# if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
      USE mod_sediment
# endif
!
# if defined EW_PERIODIC || defined NS_PERIODIC
      USE ad_exchange_2d_mod
      USE ad_exchange_3d_mod
# endif
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : ad_mp_exchange2d, ad_mp_exchange3d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp, nnew
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: Zt_avg1(LBi:,LBj:)
#  ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:,LBj:)
#  endif
#  if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
      real(r8), intent(in) :: ad_bed(LBi:,LBj:,:,:)
      real(r8), intent(inout):: ad_bed_thick0(LBi:,LBj:)
#  endif
      real(r8), intent(inout) :: ad_h(LBi:,LBj:)
      real(r8), intent(inout) :: ad_Zt_avg1(LBi:,LBj:)
      real(r8), intent(inout) :: ad_Hz(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_z_r(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_z_w(LBi:,LBj:,0:)
# else
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Zt_avg1(LBi:UBi,LBj:UBj)
#  ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:UBi,LBj:UBj)
#  endif
#  if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
      real(r8), intent(inout) :: ad_bed(LBi:UBi,LBj:UBj,Nbed,MBEDP)
      real(r8), intent(inout) :: ad_bed_thick0(LBi:UBi,LBj:UBi)
#  endif
      real(r8), intent(inout) :: ad_h(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_Zt_avg1(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: ad_z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: ad_z_w(LBi:UBi,LBj:UBj,0:N(ng))
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
      integer :: i, j, k

      real(r8) :: cff, cff_r, cff1_r, cff2_r, cff_w, cff1_w, cff2_w
      real(r8) :: hinv, hwater, z_r0, z_w0
      real(r8) :: adfac, ad_hinv, ad_hwater, ad_z_r0, ad_z_w0
      real(r8) :: ad_cff2_r, ad_cff2_w

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Initialize adjoint private variables.
!-----------------------------------------------------------------------
!
      ad_cff2_r=0.0_r8
      ad_cff2_w=0.0_r8
      ad_z_r0=0.0_r8
      ad_z_w0=0.0_r8
      ad_hinv=0.0_r8
      ad_hwater=0.0_r8
!
!-----------------------------------------------------------------------
!  Compute time evolving adjoint depths and vertical thicknesses.
!-----------------------------------------------------------------------
# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
!
!  Exchange boundary information.
!
#  ifdef DISTRIBUTE
!>    CALL mp_exchange3d (ng, tile, iTLM, 2,                            &
!>   &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
!>   &                    NghostPoints, EWperiodic, NSperiodic,         &
!>   &                    tl_z_r, tl_Hz)
!>
      CALL ad_mp_exchange3d (ng, tile, iADM, 2,                         &
     &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                       NghostPoints, EWperiodic, NSperiodic,      &
     &                       ad_z_r, ad_Hz)
!>    CALL mp_exchange3d (ng, tile, iTLM, 1,                            &
!>   &                    LBi, UBi, LBj, UBj, 0, N(ng),                 &
!>   &                    NghostPoints, EWperiodic, NSperiodic,         &
!>   &                    tl_z_w)
!>
      CALL ad_mp_exchange3d (ng, tile, iADM, 1,                         &
     &                       LBi, UBi, LBj, UBj, 0, N(ng),              &
     &                       NghostPoints, EWperiodic, NSperiodic,      &
     &                       ad_z_w)
!>    CALL mp_exchange2d (ng, tile, iTLM, 1,                            &
!>   &                    LBi, UBi, LBj, UBj,                           &
!>   &                    NghostPoints, EWperiodic, NSperiodic,         &
!>   &                    tl_h)
!>
      CALL ad_mp_exchange2d (ng, tile, iADM, 1,                         &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       NghostPoints, EWperiodic, NSperiodic,      &
     &                       ad_h)
#  endif
#  if defined EW_PERIODIC || defined NS_PERIODIC

!>    CALL exchange_r3d_tile (ng, tile,                                 &
!>   &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
!>   &                        tl_Hz)
!>
      CALL ad_exchange_r3d_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj, 1, N(ng),          &
     &                           ad_Hz)
!>    CALL exchange_r3d_tile (ng, tile,                                 &
!>   &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
!>   &                        tl_z_r)
!>
      CALL ad_exchange_r3d_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj, 1, N(ng),          &
     &                           ad_z_r)
!>    CALL exchange_w3d_tile (ng, tile,                                 &
!>   &                        LBi, UBi, LBj, UBj, 0, N(ng),             &
!>   &                        tl_z_w)
!>
      CALL ad_exchange_w3d_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj, 0, N(ng),          &
     &                           ad_z_w)
!>    CALL exchange_r2d_tile (ng, tile,                                 &
!>   &                        LBi, UBi, LBj, UBj,                       &
!>   &                        tl_h)
!>
      CALL ad_exchange_r2d_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           ad_h)
#  endif
# endif
!
!-----------------------------------------------------------------------
!  Original formulation: Compute vertical depths (meters, negative) at
!                        RHO- and W-points, and vertical grid
!  thicknesses. Various stretching functions are possible.
!
!         z_w(x,y,s,t) = Zo_w + zeta(x,y,t) * [1.0 + Zo_w / h(x,y)]
!
!                 Zo_w = hc * [s(k) - C(k)] + C(k) * h(x,y)
!
!-----------------------------------------------------------------------
!
      IF (Vtransform(ng).eq.1) THEN
        DO j=JstrR,JendR
          DO k=N(ng),1,-1
            cff_w=hc(ng)*(SCALARS(ng)%sc_w(k)-SCALARS(ng)%Cs_w(k))
            cff_r=hc(ng)*(SCALARS(ng)%sc_r(k)-SCALARS(ng)%Cs_r(k))
            cff1_r=SCALARS(ng)%Cs_r(k)
            cff1_w=SCALARS(ng)%Cs_w(k)
            DO i=IstrR,IendR
              hwater=h(i,j)
# ifdef ICESHELF
              hwater=hwater-ABS(zice(i,j))
# endif
              hinv=1.0_r8/hwater
              z_w0=cff_w+cff1_w*hwater
              z_r0=cff_r+cff1_r*hwater

!>            tl_Hz(i,j,k)=tl_z_w(i,j,k)-tl_z_w(i,j,k-1)
!>
              ad_z_w(i,j,k  )=ad_z_w(i,j,k  )+ad_Hz(i,j,k)
              ad_z_w(i,j,k-1)=ad_z_w(i,j,k-1)-ad_Hz(i,j,k)
              ad_Hz(i,j,k)=0.0_r8
!>            tl_z_r(i,j,k)=tl_z_r0+                                    &
!>   &                      tl_Zt_avg1(i,j)*(1.0_r8+z_r0*hinv)+         &
!>   &                      Zt_avg1(i,j)*(tl_z_r0*hinv+z_r0*tl_hinv)
!>
              adfac=Zt_avg1(i,j)*ad_z_r(i,j,k)
              ad_z_r0=ad_z_r0+hinv*adfac+ad_z_r(i,j,k)
              ad_hinv=ad_hinv+z_r0*adfac
              ad_Zt_avg1(i,j)=ad_Zt_avg1(i,j)+                          &
     &                        (1.0_r8+z_r0*hinv)*ad_z_r(i,j,k)
              ad_z_r(i,j,k)=0.0_r8
!>            tl_z_r0=cff1_r*tl_hwater
!>
              ad_hwater=ad_hwater+cff1_r*ad_z_r0
              ad_z_r0=0.0_r8
!>            tl_z_w(i,j,k)=tl_z_w0+                                    &
!>   &                      tl_Zt_avg1(i,j)*(1.0_r8+z_w0*hinv)+         &
!>   &                      Zt_avg1(i,j)*(tl_z_w0*hinv+z_w0*tl_hinv)
!>
              adfac=Zt_avg1(i,j)*ad_z_w(i,j,k)
              ad_z_w0=ad_z_w0+hinv*adfac+ad_z_w(i,j,k)
              ad_hinv=ad_hinv+z_w0*adfac
              ad_Zt_avg1(i,j)=ad_Zt_avg1(i,j)+                          &
     &                        (1.0_r8+z_w0*hinv)*ad_z_w(i,j,k)
              ad_z_w(i,j,k)=0.0_r8
!>            tl_z_w0=cff1_w*tl_hwater
!>
              ad_hwater=ad_hwater+cff1_w*ad_z_w0
              ad_z_w0=0.0_r8
!>            tl_hinv=-hinv*hinv*tl_hwater
!>
              ad_hwater=ad_hwater-hinv*hinv*ad_hinv
              ad_hinv=0.0_r8
!>            tl_hwater=tl_h(i,j)
!>
              ad_h(i,j)=ad_h(i,j)+ad_hwater
              ad_hwater=0.0_r8
            END DO
          END DO

          DO i=IstrR,IendR
!>          tl_z_w(i,j,0)=-tl_h(i,j)
!>
            ad_h(i,j)=ad_h(i,j)-ad_z_w(i,j,0)
            ad_z_w(i,j,0)=0.0_r8
# if defined WET_DRY
            IF (h(i,j).eq.0.0_r8) THEN
!>            tl_h(i,j)=0.0_r8
!>
              ad_h(i,j)=0.0_r8
            END IF
# endif
# if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
!>          tl_h(i,j)=tl_h(i,j)-                                        &
!>   &                tl_bed_thick(i,j,nstp)+tl_bed_thick(i,j,nnew)
!>
            ad_bed_thick(i,j,nnew)=ad_bed_thick(i,j,nnew)+              &
     &                             ad_h(i,j)

            ad_bed_thick(i,j,nstp)=ad_bed_thick(i,j,nstp)-              &
     &                             ad_h(i,j)
# endif
          END DO
        END DO
!
!-----------------------------------------------------------------------
!  New formulation: Compute vertical depths (meters, negative) at
!                   RHO- and W-points, and vertical grid thicknesses.
!  Various stretching functions are possible.
!
!         z_w(x,y,s,t) = zeta(x,y,t) + [zeta(x,y,t)+ h(x,y)] * Zo_w
!
!                 Zo_w = [hc * s(k) + C(k) * h(x,y)] / [hc + h(x,y)]
!
!-----------------------------------------------------------------------
!
      ELSE IF (Vtransform(ng).eq.2) THEN
        DO j=JstrR,JendR
          DO k=N(ng),1,-1
            cff_r=hc(ng)*(SCALARS(ng)%sc_r(k)-SCALARS(ng)%Cs_r(k))
            cff_w=hc(ng)*(SCALARS(ng)%sc_w(k)-SCALARS(ng)%Cs_w(k))
            cff1_r=SCALARS(ng)%Cs_r(k)
            cff1_w=SCALARS(ng)%Cs_w(k)
            DO i=IstrR,IendR
              hwater=h(i,j)
# ifdef ICESHELF
              hwater=hwater-ABS(zice(i,j))
# endif
              hinv=1.0_r8/hwater
              cff2_r=(cff_r+cff1_r*hwater)*hinv
              cff2_w=(cff_w+cff1_w*hwater)*hinv
!>            tl_Hz(i,j,k)=tl_z_w(i,j,k)-tl_z_w(i,j,k-1)
!>
              ad_z_w(i,j,k  )=ad_z_w(i,j,k  )+ad_Hz(i,j,k)
              ad_z_w(i,j,k-1)=ad_z_w(i,j,k-1)-ad_Hz(i,j,k)
              ad_Hz(i,j,k)=0.0_r8
!>            tl_z_r(i,j,k)=tl_Zt_avg1(i,j)+                            &
!>   &                      (tl_Zt_avg1(i,j)+tl_hwater)*cff2_r+         &
!>   &                      (Zt_avg1(i,j)+hwater)*tl_cff2_r
!>
              adfac=cff2_r*ad_z_r(i,j,k)
              ad_cff2_r=ad_cff2_r+                                      &
     &                  (Zt_avg1(i,j)+hwater)*ad_z_r(i,j,k)
              ad_hwater=ad_hwater+adfac
              ad_Zt_avg1(i,j)=ad_Zt_avg1(i,j)+ad_z_r(i,j,k)+adfac
              ad_z_r(i,j,k)=0.0_r8
!>            tl_cff2_r=cff1_r*tl_hwater*hinv+                          &
!>   &                  (cff_r+cff1_r*hwater)*tl_hinv
!>
              ad_hinv=ad_hinv+                                          &
     &                (cff_r+cff1_r*hwater)*ad_cff2_r
              ad_hwater=ad_hwater+                                      &
     &                  cff1_r*hinv*ad_cff2_r
              ad_cff2_r=0.0_r8
!>            tl_z_w(i,j,k)=tl_Zt_avg1(i,j)+                            &
!>   &                      (tl_Zt_avg1(i,j)+tl_hwater)*cff2_w+         &
!>   &                      (Zt_avg1(i,j)+hwater)*tl_cff2_w
!>
              adfac=cff2_w*ad_z_w(i,j,k)
              ad_cff2_w=ad_cff2_w+                                      &
     &                  (Zt_avg1(i,j)+hwater)*ad_z_w(i,j,k)
              ad_hwater=ad_hwater+adfac
              ad_Zt_avg1(i,j)=ad_Zt_avg1(i,j)+ad_z_w(i,j,k)+adfac
              ad_z_w(i,j,k)=0.0_r8
!>            tl_cff2_w=cff1_w*tl_hwater*hinv+                          &
!>   &                  (cff_w+cff1_w*hwater)*tl_hinv
!>
              ad_hinv=ad_hinv+                                          &
     &                (cff_w+cff1_w*hwater)*ad_cff2_w
              ad_hwater=ad_hwater+                                      &
     &                  cff1_w*hinv*ad_cff2_w
              ad_cff2_w=0.0_r8
!>            tl_hinv=-hinv*hinv*tl_hwater
!>
              ad_hwater=ad_hwater-hinv*hinv*ad_hinv
              ad_hinv=0.0_r8
!>            tl_hwater=tl_h(i,j)
!>
              ad_h(i,j)=ad_h(i,j)+ad_hwater
              ad_hwater=0.0_r8
            END DO
          END DO

          DO i=IstrR,IendR
!>          tl_z_w(i,j,0)=-tl_h(i,j)
!>
            ad_h(i,j)=ad_h(i,j)-ad_z_w(i,j,0)
            ad_z_w(i,j,0)=0.0_r8
# if defined WET_DRY
            IF (h(i,j).eq.0.0_r8) THEN
!>            tl_h(i,j)=0.0_r8
!>
              ad_h(i,j)=0.0_r8
            END IF
# endif
# if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
!>          tl_h(i,j)=tl_h(i,j)-                                        &
!>   &                tl_bed_thick(i,j,nstp)+tl_bed_thick(i,j,nnew)
!>
            ad_bed_thick(i,j,nnew)=ad_bed_thick(i,j,nnew)+              &
     &                             ad_h(i,j)

            ad_bed_thick(i,j,nstp)=ad_bed_thick(i,j,nstp)-              &
     &                             ad_h(i,j)
# endif
          END DO
        END DO
      END IF

      RETURN
      END SUBROUTINE ad_set_depth_tile
#endif
      END MODULE ad_set_depth_mod
