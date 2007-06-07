#include "cppdefs.h"
      MODULE metrics_mod
!
!svn $Id$
!=======================================================================
!  Copyright (c) 2002-2007 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This routine computes various horizontal metric terms.              !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: metrics

      CONTAINS
!
!***********************************************************************
      SUBROUTINE metrics (ng, tile, model)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_ocean
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
#include "tile.h"
!
      CALL metrics_tile (ng, model, Istr, Iend, Jstr, Jend,             &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   nstp(ng), nnew(ng),                            &
     &                   GRID(ng) % f,                                  &
     &                   GRID(ng) % h,                                  &
     &                   GRID(ng) % pm,                                 &
     &                   GRID(ng) % pn,                                 &
#ifdef MASKING
     &                   GRID(ng) % pmask,                              &
     &                   GRID(ng) % rmask,                              &
#endif
#ifdef CURVGRID
     &                   GRID(ng) % angler,                             &
     &                   GRID(ng) % CosAngler,                          &
     &                   GRID(ng) % SinAngler,                          &
#endif
#ifdef SOLVE3D
# ifdef ICESHELF
     &                   GRID(ng) % zice,                               &
# endif
# if defined SEDIMENT && defined SED_MORPH
     &                   GRID(ng) % bed_thick,                          &
# endif
     &                   GRID(ng) % Hz,                                 &
     &                   GRID(ng) % z_r,                                &
     &                   GRID(ng) % z_w,                                &
#endif
#if defined VISC_GRID || defined DIFF_GRID
     &                   GRID(ng) % grdscl,                             &
#endif
     &                   GRID(ng) % om_p,                               &
     &                   GRID(ng) % om_r,                               &
     &                   GRID(ng) % om_u,                               &
     &                   GRID(ng) % om_v,                               &
     &                   GRID(ng) % on_p,                               &
     &                   GRID(ng) % on_r,                               &
     &                   GRID(ng) % on_u,                               &
     &                   GRID(ng) % on_v,                               &
     &                   GRID(ng) % fomn,                               &
     &                   GRID(ng) % omn,                                &
     &                   GRID(ng) % pnom_p,                             &
     &                   GRID(ng) % pnom_r,                             &
     &                   GRID(ng) % pnom_u,                             &
     &                   GRID(ng) % pnom_v,                             &
     &                   GRID(ng) % pmon_p,                             &
     &                   GRID(ng) % pmon_r,                             &
     &                   GRID(ng) % pmon_u,                             &
     &                   GRID(ng) % pmon_v)
      RETURN
      END SUBROUTINE metrics
!
!***********************************************************************
      SUBROUTINE metrics_tile (ng, model, Istr, Iend, Jstr, Jend,       &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         nstp, nnew,                              &
     &                         f, h, pm, pn,                            &
#ifdef MASKING
     &                         pmask, rmask,                            &
#endif
#ifdef CURVGRID
     &                         angler, CosAngler, SinAngler,            &
#endif
#ifdef SOLVE3D
# ifdef ICESHELF
     &                         zice,                                    &
# endif
# if defined SEDIMENT && defined SED_MORPH
     &                         bed_thick,                               &
# endif
     &                         Hz, z_r, z_w,                            &
#endif
#if defined VISC_GRID || defined DIFF_GRID
     &                         grdscl,                                  &
#endif
     &                         om_p, om_r, om_u, om_v,                  &
     &                         on_p, on_r, on_u, on_v,                  &
     &                         fomn, omn,                               &
     &                         pnom_p, pnom_r, pnom_u, pnom_v,          &
     &                         pmon_p, pmon_r, pmon_u, pmon_v)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
#ifdef FOUR_DVAR
      USE mod_fourdvar
#endif
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
      USE mod_iounits
!
#ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_reduce
#endif
#if defined EW_PERIODIC || defined NS_PERIODIC
      USE exchange_2d_mod
#endif
#ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
#endif
#ifdef SOLVE3D
      USE set_depth_mod, ONLY : set_depth_tile
#endif
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: nstp, nnew

#ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: f(LBi:,LBj:)
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)

      real(r8), intent(inout) :: h(LBi:,LBj:)
# ifdef MASKING
      real(r8), intent(inout) :: pmask(LBi:,LBj:)
      real(r8), intent(inout) :: rmask(LBi:,LBj:)
# endif
# ifdef CURVGRID
      real(r8), intent(inout) :: angler(LBi:,LBj:)
# endif
# ifdef SOLVE3D
#  ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:,LBj:)
#  endif
#  if defined SEDIMENT && defined SED_MORPH
      real(r8), intent(in):: bed_thick(LBi:,LBj:,:)
#  endif
# endif
# if defined VISC_GRID || defined DIFF_GRID
      real(r8), intent(inout) :: grdscl(LBi:,LBj:)
# endif
      real(r8), intent(out) :: om_p(LBi:,LBj:)
      real(r8), intent(out) :: om_r(LBi:,LBj:)
      real(r8), intent(out) :: om_u(LBi:,LBj:)
      real(r8), intent(out) :: om_v(LBi:,LBj:)
      real(r8), intent(out) :: on_p(LBi:,LBj:)
      real(r8), intent(out) :: on_r(LBi:,LBj:)
      real(r8), intent(out) :: on_u(LBi:,LBj:)
      real(r8), intent(out) :: on_v(LBi:,LBj:)
      real(r8), intent(out) :: fomn(LBi:,LBj:)
      real(r8), intent(out) :: omn(LBi:,LBj:)
      real(r8), intent(out) :: pnom_p(LBi:,LBj:)
      real(r8), intent(out) :: pnom_r(LBi:,LBj:)
      real(r8), intent(out) :: pnom_u(LBi:,LBj:)
      real(r8), intent(out) :: pnom_v(LBi:,LBj:)
      real(r8), intent(out) :: pmon_p(LBi:,LBj:)
      real(r8), intent(out) :: pmon_r(LBi:,LBj:)
      real(r8), intent(out) :: pmon_u(LBi:,LBj:)
      real(r8), intent(out) :: pmon_v(LBi:,LBj:)
# ifdef CURVGRID
      real(r8), intent(out) :: CosAngler(LBi:,LBj:)
      real(r8), intent(out) :: SinAngler(LBi:,LBj:)
# endif
# ifdef SOLVE3D
      real(r8), intent(out) :: Hz(LBi:,LBj:,:)
      real(r8), intent(out) :: z_r(LBi:,LBj:,:)
      real(r8), intent(out) :: z_w(LBi:,LBj:,0:)
# endif
#else
      real(r8), intent(in) :: f(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)

      real(r8), intent(inout) :: h(LBi:UBi,LBj:UBj)
# ifdef MASKING
      real(r8), intent(inout) :: pmask(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: rmask(LBi:UBi,LBj:UBj)
# endif
# ifdef CURVGRID
      real(r8), intent(inout) :: angler(LBi:UBi,LBj:UBj)
# endif
# ifdef SOLVE3D
#  ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:UBi,LBj:UBj)
#  endif
#  if defined SEDIMENT && defined SED_MORPH
      real(r8), intent(in):: bed_thick(LBi:UBi,LBj:UBj,2)
#  endif
# endif
# if defined VISC_GRID || defined DIFF_GRID
      real(r8), intent(inout) :: grdscl(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(out) :: om_p(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: om_r(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: om_u(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: om_v(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: on_p(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: on_r(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: on_u(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: on_v(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: fomn(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: omn(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: pnom_p(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: pnom_r(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: pnom_u(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: pnom_v(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: pmon_p(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: pmon_r(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: pmon_u(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: pmon_v(LBi:UBi,LBj:UBj)
# ifdef CURVGRID
      real(r8), intent(out) :: CosAngler(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: SinAngler(LBi:UBi,LBj:UBj)
# endif
# ifdef SOLVE3D
      real(r8), intent(out) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(out) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(out) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
# endif
#endif
!
!  Local variable declarations.
!
#ifdef DISTRIBUTE
# ifdef EW_PERIODIC
      logical :: EWperiodic=.TRUE.
# else
      logical :: EWperiodic=.FALSE.
# endif
# ifdef NS_PERIODIC
      logical :: NSperiodic=.TRUE.
# else
      logical :: NSperiodic=.FALSE.
# endif
#endif
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: NSUB, i, is, j, k

      real(r8), parameter :: Large = 1.0E+20_r8

      real(r8) :: cff, cff1, cff2
      real(r8) :: my_DXmax, my_DXmin, my_DYmax, my_DYmin
#ifdef SOLVE3D
      real(r8) :: my_DZmax, my_DZmin
#endif
      real(r8) :: my_Cu_Cor, my_Cu_max, my_Cu_min, my_grdmax

#ifdef DISTRIBUTE
      integer :: size
      real(r8), dimension(10) :: buffer
      character (len=3), dimension(10) :: op_handle
#endif

#ifdef SOLVE3D
      real(r8), dimension(LBi:UBi,LBj:UBj) :: A2d
#endif

#include "set_bounds.h"


#ifdef EW_PERIODIC
# define IU_RANGE Istr,Iend
# define IV_RANGE Istr,Iend
#else
# define IU_RANGE Istr,IendR
# define IV_RANGE IstrR,IendR
#endif
#ifdef NS_PERIODIC
# define JU_RANGE Jstr,Jend
# define JV_RANGE Jstr,Jend
#else
# define JU_RANGE JstrR,JendR
# define JV_RANGE Jstr,JendR
#endif
!
!-----------------------------------------------------------------------
!  Compute 1/m, 1/n, 1/mn, and f/mn at horizontal RHO-points.
!-----------------------------------------------------------------------
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          om_r(i,j)=1.0_r8/pm(i,j)
          on_r(i,j)=1.0_r8/pn(i,j)
          omn(i,j)=1.0_r8/(pm(i,j)*pn(i,j))
          fomn(i,j)=f(i,j)*omn(i,j)
        END DO
      END DO
#if defined EW_PERIODIC || defined NS_PERIODIC
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        om_r)
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        on_r)
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        omn)
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        fomn)
#endif
#ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, model, 4, Istr, Iend, Jstr, Jend,         &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    om_r, on_r, omn, fomn)
#endif
!
!-----------------------------------------------------------------------
!  Compute n/m, and m/n at horizontal RHO-points.
!-----------------------------------------------------------------------
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          pnom_r(i,j)=pn(i,j)/pm(i,j)
          pmon_r(i,j)=pm(i,j)/pn(i,j)
        END DO
      END DO
#if defined EW_PERIODIC || defined NS_PERIODIC
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        pnom_r)
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        pmon_r)
#endif
#ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, model, 2, Istr, Iend, Jstr, Jend,         &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    pnom_r, pmon_r)
#endif
!
!-----------------------------------------------------------------------
!  Compute m/n, 1/m, and 1/n at horizontal U-points.
!-----------------------------------------------------------------------
!
      DO j=JU_RANGE
        DO i=IU_RANGE
          pmon_u(i,j)=(pm(i-1,j)+pm(i,j))/(pn(i-1,j)+pn(i,j))
          pnom_u(i,j)=(pn(i-1,j)+pn(i,j))/(pm(i-1,j)+pm(i,j))
          om_u(i,j)=2.0_r8/(pm(i-1,j)+pm(i,j))
          on_u(i,j)=2.0_r8/(pn(i-1,j)+pn(i,j))
        END DO
      END DO
#if defined EW_PERIODIC || defined NS_PERIODIC
      CALL exchange_u2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        pmon_u)
      CALL exchange_u2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        pnom_u)
      CALL exchange_u2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        om_u)
      CALL exchange_u2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        on_u)
#endif
#ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, model, 4, Istr, Iend, Jstr, Jend,         &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    pmon_u, pnom_u, om_u, on_u)
#endif
!
!-----------------------------------------------------------------------
!  Compute n/m, 1/m, and 1/m at horizontal V-points.
!-----------------------------------------------------------------------
!
      DO j=JV_RANGE
        DO i=IV_RANGE
          pmon_v(i,j)=(pm(i,j-1)+pm(i,j))/(pn(i,j-1)+pn(i,j))
          pnom_v(i,j)=(pn(i,j-1)+pn(i,j))/(pm(i,j-1)+pm(i,j))
          om_v(i,j)=2.0_r8/(pm(i,j-1)+pm(i,j))
          on_v(i,j)=2.0_r8/(pn(i,j-1)+pn(i,j))
        END DO
      END DO
#if defined EW_PERIODIC || defined NS_PERIODIC
      CALL exchange_v2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        pmon_v)
      CALL exchange_v2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        pnom_v)
      CALL exchange_v2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        om_v)
      CALL exchange_v2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        on_v)
#endif
#ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, model, 4, Istr, Iend, Jstr, Jend,         &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    pmon_v, pnom_v, om_v, on_v)
#endif
!
!-----------------------------------------------------------------------
!  Compute n/m and m/n at horizontal PSI-points.
!-----------------------------------------------------------------------
!
      DO j=JV_RANGE
        DO i=IU_RANGE
          pnom_p(i,j)=(pn(i-1,j-1)+pn(i-1,j)+pn(i,j-1)+pn(i,j))/        &
     &                (pm(i-1,j-1)+pm(i-1,j)+pm(i,j-1)+pm(i,j))
          pmon_p(i,j)=(pm(i-1,j-1)+pm(i-1,j)+pm(i,j-1)+pm(i,j))/        &
     &                (pn(i-1,j-1)+pn(i-1,j)+pn(i,j-1)+pn(i,j))
          om_p(i,j)=4.0_r8/(pm(i-1,j-1)+pm(i-1,j)+pm(i,j-1)+pm(i,j))
          on_p(i,j)=4.0_r8/(pn(i-1,j-1)+pn(i-1,j)+pn(i,j-1)+pn(i,j))
        END DO
      END DO
#if defined EW_PERIODIC || defined NS_PERIODIC
      CALL exchange_p2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        pnom_p)
      CALL exchange_p2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        pmon_p)
      CALL exchange_p2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        om_p)
      CALL exchange_p2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        on_p)
#endif
#ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, model, 4, Istr, Iend, Jstr, Jend,         &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    pnom_p, pmon_p, om_p, on_p)
#endif
#ifdef MASKING
!
!-----------------------------------------------------------------------
!  Set slipperiness (no-slip) mask at PSI-points.
!-----------------------------------------------------------------------
!
! Set no-slip boundary conditions on land-mask boundaries regardless of
! supplied value of gamma2. 
! 
      cff1=1.0_r8       ! computation of off-diagonal nonlinear terms 
      cff2=2.0_r8
      DO j=Jstr,Jend
        DO i=Istr,Iend
          IF ((rmask(i-1,j  ).gt.0.5_r8).and.                           &
     &        (rmask(i  ,j  ).gt.0.5_r8).and.                           &
     &        (rmask(i-1,j-1).gt.0.5_r8).and.                           &
     &        (rmask(i  ,j-1).gt.0.5_r8)) THEN 
            pmask(i,j)=1.0_r8
          ELSE IF ((rmask(i-1,j  ).lt.0.5_r8).and.                      &
     &             (rmask(i  ,j  ).gt.0.5_r8).and.                      &
     &             (rmask(i-1,j-1).gt.0.5_r8).and.                      &
     &             (rmask(i  ,j-1).gt.0.5_r8)) THEN 
            pmask(i,j)=cff1
          ELSE IF ((rmask(i-1,j  ).gt.0.5_r8).and.                      &
     &             (rmask(i  ,j  ).lt.0.5_r8).and.                      &
     &             (rmask(i-1,j-1).gt.0.5_r8).and.                      &
     &             (rmask(i  ,j-1).gt.0.5_r8)) THEN 
            pmask(i,j)=cff1
          ELSE IF ((rmask(i-1,j  ).gt.0.5_r8).and.                      &
     &             (rmask(i  ,j  ).gt.0.5_r8).and.                      &
     &             (rmask(i-1,j-1).lt.0.5_r8).and.                      &
     &             (rmask(i  ,j-1).gt.0.5_r8)) THEN 
            pmask(i,j)=cff1
          ELSE IF ((rmask(i-1,j  ).gt.0.5_r8).and.                      &
     &             (rmask(i  ,j  ).gt.0.5_r8).and.                      &
     &             (rmask(i-1,j-1).gt.0.5_r8).and.                      &
     &             (rmask(i  ,j-1).lt.0.5_r8)) THEN 
            pmask(i,j)=cff1
          ELSE IF ((rmask(i-1,j  ).gt.0.5_r8).and.                      &
     &             (rmask(i  ,j  ).lt.0.5_r8).and.                      &
     &             (rmask(i-1,j-1).gt.0.5_r8).and.                      &
     &             (rmask(i  ,j-1).lt.0.5_r8)) THEN 
            pmask(i,j)=cff2
          ELSE IF ((rmask(i-1,j  ).lt.0.5_r8).and.                      &
     &             (rmask(i  ,j  ).gt.0.5_r8).and.                      &
     &             (rmask(i-1,j-1).lt.0.5_r8).and.                      &
     &             (rmask(i  ,j-1).gt.0.5_r8)) THEN 
            pmask(i,j)=cff2
          ELSE IF ((rmask(i-1,j  ).gt.0.5_r8).and.                      &
     &             (rmask(i  ,j  ).gt.0.5_r8).and.                      &
     &             (rmask(i-1,j-1).lt.0.5_r8).and.                      &
     &             (rmask(i  ,j-1).lt.0.5_r8)) THEN 
            pmask(i,j)=cff2
          ELSE IF ((rmask(i-1,j  ).lt.0.5_r8).and.                      &
     &             (rmask(i  ,j  ).lt.0.5_r8).and.                      &
     &             (rmask(i-1,j-1).gt.0.5_r8).and.                      &
     &             (rmask(i  ,j-1).gt.0.5_r8)) THEN
            pmask(i,j)=cff2
          ELSE 
            pmask(i,j)=0.0_r8
          END IF
        END DO
      END DO
# if defined EW_PERIODIC || defined NS_PERIODIC
      CALL exchange_p2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        pmask)
# endif
# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, model, 1, Istr, Iend, Jstr, Jend,         &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    pmask)
# endif
#endif
#ifdef CURVGRID
!
!-----------------------------------------------------------------------
! Compute cosine and sine of curvilinear rotation angle.
!-----------------------------------------------------------------------
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          CosAngler(i,j)=COS(angler(i,j))
          SinAngler(i,j)=SIN(angler(i,j))
        END DO
      END DO
# if defined EW_PERIODIC || defined NS_PERIODIC
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        CosAngler)
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        SinAngler)
# endif
# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, model, 2, Istr, Iend, Jstr, Jend,         &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    CosAngler, SinAngler)
# endif
#endif
#if defined VISC_GRID || defined DIFF_GRID
!
!-----------------------------------------------------------------------
! Determine maximum area of the grid box over the whole curvilinear
! grid of the model. Used for scaling horizontal mixing by the grid
! size.
!-----------------------------------------------------------------------
!
      cff=0.0_r8
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          grdscl(i,j)=SQRT(om_r(i,j)*om_r(i,j)+on_r(i,j)*on_r(i,j))
        END DO
      END DO
# if defined EW_PERIODIC || defined NS_PERIODIC
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        grdscl)
# endif
# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, model, 1, Istr, Iend, Jstr, Jend,         &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    grdscl)
# endif
#endif
!
!-----------------------------------------------------------------------
!  Compute minimum and maximum grid spacing.
!-----------------------------------------------------------------------
#ifdef SOLVE3D
!
!  Compute time invariant depths (use zero free-surface).
!
      DO i=LBi,UBi
        DO j=LBj,UBj
          A2d(i,j)=0.0_r8
        END DO
      END DO

      CALL set_depth_tile (ng, Istr, Iend, Jstr, Jend,                  &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     nstp, nnew,                                  &
     &                     h,                                           &
# ifdef ICESHELF
     &                     zice,                                        &
# endif
# if defined SEDIMENT && defined SED_MORPH
     &                     bed_thick,                                   &
# endif
     &                     A2d,                                         &
     &                     Hz, z_r, z_w)
#endif
!
!  Compute grid spacing range.
!
      my_DXmin= Large
      my_DXmax=-Large
      my_DYmin= Large
      my_DYmax=-Large
#ifdef SOLVE3D
      my_DZmin= Large
      my_DZmax=-Large
#endif
      my_grdmax=-Large
      DO j=JstrR,JendR
        DO i=IstrR,IendR
#if defined VISC_GRID || defined DIFF_GRID
          cff=grdscl(i,j)
#else
          cff=SQRT(om_r(i,j)*om_r(i,j)+on_r(i,j)*on_r(i,j))
#endif
          my_DXmin=MIN(my_DXmin,om_r(i,j))
          my_DXmax=MAX(my_DXmax,om_r(i,j))
          my_DYmin=MIN(my_DYmin,on_r(i,j))
          my_DYmax=MAX(my_DYmax,on_r(i,j))
          my_grdmax=MAX(my_grdmax,cff)
#ifdef SOLVE3D
          DO k=1,N(ng)
            my_DZmin=MIN(my_DZmin,Hz(i,j,k))
            my_DZmax=MAX(my_DZmax,Hz(i,j,k))
          END DO
#endif
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Compute Courant number.
!-----------------------------------------------------------------------
!
!  The Courant number is defined as:
!
!     Cu = c * dt * SQRT (1/dx^2 + 1/dy^2)
!
!  where c=SQRT(g*h) is phase speed for barotropic mode, and dx, dy
!  are grid spacing in each direction.
!
      my_Cu_min= Large
      my_Cu_max=-Large
      my_Cu_Cor=-Large

      DO j=JstrR,JendR
        DO i=IstrR,IendR
#ifdef MASKING
          IF (rmask(i,j).gt.0.0_r8) THEN
#endif
            cff=dtfast(ng)*                                             &
     &          SQRT(g*ABS(h(i,j))*(pm(i,j)*pm(i,j)+pn(i,j)*pn(i,j)))
            my_Cu_min=MIN(my_Cu_min,cff)
            my_Cu_max=MAX(my_Cu_max,cff)

            cff=dt(ng)*ABS(f(i,j))
            my_Cu_Cor=MAX(my_Cu_Cor,cff)

#ifdef MASKING
          END IF
#endif
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Perform global reductions.
!-----------------------------------------------------------------------
!
      IF (SOUTH_WEST_CORNER.and.                                        &
     &    NORTH_EAST_CORNER) THEN
        NSUB=1                           ! non-tiled application
      ELSE
        NSUB=NtileX(ng)*NtileE(ng)       ! tiled application
      END IF
!$OMP CRITICAL (REDUCTIONS)
      IF (tile_count.eq.0) THEN
        Cu_min=my_Cu_min
        Cu_max=my_Cu_max
        Cu_Cor=my_Cu_Cor
        grdmax(ng)=my_grdmax
        DXmin(ng)=my_DXmin
        DXmax(ng)=my_DXmax
        DYmin(ng)=my_DYmin
        DYmax(ng)=my_DYmax
#ifdef SOLVE3D
        DZmin(ng)=my_DZmin
        DZmax(ng)=my_DZmax
#endif
      ELSE
        Cu_min=MIN(Cu_min,my_Cu_min)
        Cu_max=MAX(Cu_max,my_Cu_max)
        Cu_Cor=MAX(Cu_Cor,my_Cu_Cor)
        grdmax(ng)=MAX(grdmax(ng),my_grdmax)
        DXmin(ng)=MIN(DXmin(ng),my_DXmin)
        DXmax(ng)=MAX(DXmax(ng),my_DXmax)
        DYmin(ng)=MIN(DYmin(ng),my_DYmin)
        DYmax(ng)=MAX(DYmax(ng),my_DYmax)
#ifdef SOLVE3D
        DZmin(ng)=MIN(DZmin(ng),my_DZmin)
        DZmax(ng)=MAX(DZmax(ng),my_DZmax)
#endif
      END IF
      tile_count=tile_count+1
      IF (tile_count.eq.NSUB) THEN
        tile_count=0
#ifdef DISTRIBUTE
        buffer(1)=Cu_min
        op_handle(1)='MIN'
        buffer(2)=Cu_max
        op_handle(2)='MAX'
        buffer(3)=Cu_Cor
        op_handle(3)='MAX'
        buffer(4)=grdmax(ng)
        op_handle(4)='MAX'
        buffer(5)=DXmin(ng)
        op_handle(5)='MIN'
        buffer(6)=DXmax(ng)
        op_handle(6)='MAX'
        buffer(7)=DYmin(ng)
        op_handle(7)='MIN'
        buffer(8)=DYmax(ng)
        op_handle(8)='MAX'
# ifdef SOLVE3D
        buffer(9)=DZmin(ng)
        op_handle(9)='MIN'
        buffer(10)=DZmax(ng)
        op_handle(10)='MAX'
        size=10
# else
        size=8
# endif
        CALL mp_reduce (ng, model, size, buffer, op_handle)
        Cu_min=buffer(1)
        Cu_max=buffer(2)
        Cu_Cor=buffer(3)
        grdmax(ng)=buffer(4)
        DXmin(ng)=buffer(5)
        DXmax(ng)=buffer(6)
        DYmin(ng)=buffer(7)
        DYmax(ng)=buffer(8)
# ifdef SOLVE3D
        DZmin(ng)=buffer(9)
        DZmax(ng)=buffer(10)
# endif
#endif
        IF (Master) THEN
          WRITE(stdout,10) DXmin(ng)/1000.0_r8, DXmax(ng)/1000.0_r8,    &
     &                     DYmin(ng)/1000.0_r8, DYmax(ng)/1000.0_r8
  10      FORMAT (/,' Minimum X-grid spacing, DXmin = ',1pe15.8,' km',  &
     &            /,' Maximum X-grid spacing, DXmax = ',1pe15.8,' km',  &
     &            /,' Minimum Y-grid spacing, DYmin = ',1pe15.8,' km',  &
     &            /,' Maximum Y-grid spacing, DYmax = ',1pe15.8,' km')
#ifdef SOLVE3D
          WRITE(stdout,20) DZmin(ng), DZmax(ng)
  20      FORMAT (' Minimum Z-grid spacing, DZmin = ',1pe15.8,' m',/,   &
     &            ' Maximum Z-grid spacing, DZmax = ',1pe15.8,' m')
#endif
          WRITE (stdout,30) Cu_min, Cu_max, Cu_Cor
  30      FORMAT (/,' Minimum barotropic Courant Number = ', 1pe15.8,/, &
     &              ' Maximum barotropic Courant Number = ', 1pe15.8,/, &
     &              ' Maximum Coriolis   Courant Number = ', 1pe15.8,/)
# if defined VISC_GRID || defined DIFF_GRID
          WRITE (stdout,40) grdmax(ng)/1000.0_r8
  40      FORMAT (' Horizontal mixing scaled by grid size,',            &
     &            ' GRDMAX = ',1pe15.8,' km',/)
# endif
        END IF
      END IF
!$OMP END CRITICAL (REDUCTIONS)

#ifdef FOUR_DVAR
!
!-----------------------------------------------------------------------
!  Spatial convolution parameters.
!-----------------------------------------------------------------------
!
!  Determine diffusion operator time-step size using FTCS stability
!  criteria:
!
!       Kh DTsizeH / MIN(DXmin,DYmin)^2  < Hgamma / 4
!
!       Kv DTsizeH / DZmin^2  < Vgamma / 2
!
!  where a Hgamma and Vgamma are used to scale the time-step below
!  its theoretical limit for stability and accurary.
!
      cff=MIN(DXmin(ng),DYmin(ng))
      DO is=1,NstateVar(ng)
        DTsizeH(is)=Hgamma*cff*cff/(4.0_r8*KhMax(ng))
# ifdef SOLVE3D
#  ifdef IMPLICIT_VCONV
        DTsizeV(is)=Vgamma*DZmax(ng)*DZmax(ng)/(2.0_r8*KvMax(ng))
#  else
        DTsizeV(is)=Vgamma*DZmin(ng)*DZmin(ng)/(2.0_r8*KvMax(ng))
#  endif
# endif
      END DO
!
!  Determine FULL number of integeration steps as function of the
!  spatial decorrelation scale (Hdecay, Vdecay). Notice that in
!  the squared-root filter only half of number of step is used.
!  Thefore, the number of steps are forced to be even.
!
      DO is=1,NstateVar(ng)
        NHsteps(is)=NINT(Hdecay(is,ng)*Hdecay(is,ng)/                   &
     &                   (2.0_r8*KhMax(ng)*DTsizeH(is)))
        IF (MOD(NHsteps(is),2).ne.0) THEN
          NHsteps(is)=NHsteps(is)+1
        END IF
# ifdef SOLVE3D
        NVsteps(is)=NINT(Vdecay(is,ng)*Vdecay(is,ng)/                   &
     &                   (2.0_r8*KvMax(ng)*DTsizeV(is)))
#  ifdef IMPLICIT_VCONV
        NVsteps(is)=MAX(1,NVsteps(is))
#  endif
        IF (MOD(NVsteps(is),2).ne.0) THEN
          NVsteps(is)=NVsteps(is)+1
        END IF
# endif
      END DO
!
!  Report convolution parameters.
!
      IF (Master) THEN
        DO is=1,NstateVar(ng)
          WRITE (stdout,50) NHsteps(is), DTsizeH(ng),                   &
     &                      TRIM(Vname(1,idSvar(is)))
  50      FORMAT (' Horizontal convolution, NHsteps, DTsizeH = ',       &
     &            i5,1x,1pe15.8,' s',2x,a)
        END DO
# if defined SOLVE3D && defined VCONVOLUTION
        WRITE (stdout,'(1x)')
        DO is=4,NstateVar(ng)
          WRITE (stdout,60) NVsteps(is), DTsizeV(is),                   &
     &                      TRIM(Vname(1,idSvar(is)))
  60      FORMAT (' Vertical   convolution, NVsteps, DTsizeV = ',       &
     &            i5,1x,1pe15.8,' s',2x,a)
        END DO
# endif
        WRITE (stdout,'(/)')
      END IF
#endif

      RETURN
      END SUBROUTINE metrics_tile

      END MODULE metrics_mod
