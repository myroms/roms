#include "cppdefs.h"
      MODULE metrics_mod
!
!========================================== Alexander F. Shchepetkin ===
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
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
      SUBROUTINE metrics (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
!
      integer, intent(in) :: ng, tile

#include "tile.h"
!
      CALL metrics_tile (ng, Istr, Iend, Jstr, Jend,                    &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   GRID(ng) % f,                                  &
     &                   GRID(ng) % h,                                  &
     &                   GRID(ng) % pm,                                 &
     &                   GRID(ng) % pn,                                 &
#ifdef MASKING
     &                   GRID(ng) % pmask,                              &
     &                   GRID(ng) % rmask,                              &
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
      SUBROUTINE metrics_tile (ng, Istr, Iend, Jstr, Jend,              &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         f, h, pm, pn,                            &
#ifdef MASKING
     &                         pmask, rmask,                            &
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
      USE mod_iounits
      USE mod_scalars
      USE mod_iounits
!
#if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      USE exchange_2d_mod
!
#endif
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj

#ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: f(LBi:,LBj:)
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)

# ifdef MASKING
      real(r8), intent(inout) :: pmask(LBi:,LBj:)
      real(r8), intent(inout) :: rmask(LBi:,LBj:)
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
#else
      real(r8), intent(in) :: f(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)

# ifdef MASKING
      real(r8), intent(inout) :: pmask(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: rmask(LBi:UBi,LBj:UBj)
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
#endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: NSUB, i, j

      real(r8) :: cff
      real(r8) :: my_Cu_Cor, my_Cu_max, my_Cu_min, my_grdmax

#ifdef DISTRIBUTE
      integer :: size
      real(r8), dimension(4) :: buffer
      character (len=3), dimension(4) :: op_handle
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
#if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, om_r)
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, on_r)
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, omn)
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, fomn)
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
#if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, pnom_r)
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, pmon_r)
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
#if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      CALL exchange_u2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, pmon_u)
      CALL exchange_u2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, pnom_u)
      CALL exchange_u2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, om_u)
      CALL exchange_u2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, on_u)
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
#if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      CALL exchange_v2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, pmon_v)
      CALL exchange_v2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, pnom_v)
      CALL exchange_v2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, om_v)
      CALL exchange_v2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, on_v)
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
#if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      CALL exchange_p2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, pnom_p)
      CALL exchange_p2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, pmon_p)
      CALL exchange_p2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, om_p)
      CALL exchange_p2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, on_p)
#endif
#ifdef MASKING
!
!-----------------------------------------------------------------------
!  Set slipperiness (no-slip) mask at PSI-points.
!-----------------------------------------------------------------------
!
!  The slipperiness mask is Sea=1, Land=0, and boundary=1-gamma2.
!
      IF (gamma2.lt.0.0_r8) THEN
        DO j=Jstr,Jend
          DO i=Istr,Iend
            pmask(i,j)=2.0_r8-pmask(i,j)
          END DO
        END DO
# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
        CALL exchange_p2d_tile (ng, Istr, Iend, Jstr, Jend,             &
     &                          LBi, UBi, LBj, UBj, pmask)
# endif
      END IF
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
# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, grdscl)
# endif
#endif
!
!-----------------------------------------------------------------------
!  Compute Courant number.  Perform global reduction operations.
!-----------------------------------------------------------------------
!
!  The Courant number is defined as:
!
!     Cu = c * dt * SQRT (1/dx^2 + 1/dy^2)
!
!  where c=SQRT(g*h) is phase speed for barotropic mode, and dx, dy
!  are grid spacing in each direction.
!
      my_Cu_min= 1.0E+20_r8
      my_Cu_max=-1.0E+20_r8
      my_Cu_Cor=-1.0E+20_r8
      my_grdmax=-1.0E+20_r8

      DO j=JstrR,JendR
        DO i=IstrR,IendR
#ifdef MASKING
          IF (rmask(i,j).gt.0.0_r8) THEN
#endif
            cff=dtfast(ng)*SQRT(g*h(i,j)*                               &
     &          (pm(i,j)*pm(i,j)+pn(i,j)*pn(i,j)))
            my_Cu_min=MIN(my_Cu_min,cff)
            my_Cu_max=MAX(my_Cu_max,cff)

            cff=dt(ng)*ABS(f(i,j))
            my_Cu_Cor=MAX(my_Cu_Cor,cff)

#if defined VISC_GRID || defined DIFF_GRID
            my_grdmax=MAX(my_grdmax,grdscl(i,j))
#endif
#ifdef MASKING
          END IF
#endif
        END DO
      END DO
!
!  Perform global reductions.
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
#if defined VISC_GRID || defined DIFF_GRID
        grdmax(ng)=my_grdmax
#endif
      ELSE
        Cu_min=MIN(Cu_min,my_Cu_min)
        Cu_max=MAX(Cu_max,my_Cu_max)
        Cu_Cor=MAX(Cu_Cor,my_Cu_Cor)
#if defined VISC_GRID || defined DIFF_GRID
        grdmax(ng)=MAX(grdmax(ng),my_grdmax)
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
# if defined VISC_GRID || defined DIFF_GRID
        buffer(4)=grdmax(ng)
        op_handle(4)='MAX'
        size=4
# else
        size=3
# endif
        CALL mp_reduce (ng, size, buffer, op_handle)
        Cu_min=buffer(1)
        Cu_max=buffer(2)
        Cu_Cor=buffer(3)
# if defined VISC_GRID || defined DIFF_GRID
        grdmax(ng)=buffer(4)
# endif
#endif
        IF (Master) WRITE (stdout,10) Cu_min, Cu_max, Cu_Cor
  10    FORMAT (/,' Minimum barotropic Courant Number = ', 1pe15.8,/,   &
     &            ' Maximum barotropic Courant Number = ', 1pe15.8,/,   &
     &            ' Maximum Coriolis   Courant Number = ', 1pe15.8)

# if defined VISC_GRID || defined DIFF_GRID
        IF (Master) WRITE (stdout,20) grdmax(ng)/1000.0_r8
  20    FORMAT (/,' Horizontal mixing scaled by grid size,',            &
     &            ' GRDMAX = ',1pe15.8,' km',/)
# endif
      END IF
!$OMP END CRITICAL (REDUCTIONS)

      RETURN
      END SUBROUTINE metrics_tile

      END MODULE metrics_mod
