#include "cppdefs.h"
      MODULE set_vbc_mod
#ifdef SOLVE3D
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This module sets vertical boundary conditons for momentum and       !
!  tracers.                                                            !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC set_vbc

      CONTAINS
!
!***********************************************************************
      SUBROUTINE set_vbc (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_forces
      USE mod_ocean
      USE mod_stepping
!
      integer, intent(in) :: ng, tile

# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, 6)
# endif
      CALL set_vbc_tile (ng, Istr, Iend, Jstr, Jend,                    &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   nrhs(ng),                                      &
     &                   GRID(ng) % Hz,                                 &
     &                   GRID(ng) % z_r,                                &
     &                   GRID(ng) % z_w,                                &
# if defined ICESHELF
     &                   GRID(ng) % zice,                               &
# endif
     &                   OCEAN(ng) % t,                                 &
     &                   OCEAN(ng) % u,                                 &
     &                   OCEAN(ng) % v,                                 &
# ifdef QCORRECTION
     &                   FORCES(ng) % dqdt,                             &
     &                   FORCES(ng) % sst,                              &
# endif
# if defined SCORRECTION || defined SRELAXATION
     &                   FORCES(ng) % sss,                              &
# endif
# if defined ICESHELF
#  ifdef SHORTWAVE
     &                   FORCES(ng) % srflx,                            &
#  endif
     &                   FORCES(ng) % sustr,                            &
     &                   FORCES(ng) % svstr,                            &
# endif
# ifndef BBL
     &                   FORCES(ng) % bustr,                            &
     &                   FORCES(ng) % bvstr,                            &
# endif
     &                   FORCES(ng) % stflx,                            &
     &                   FORCES(ng) % btflx)
# ifdef PROFILE
      CALL wclock_off (ng, 6)
# endif
      RETURN
      END SUBROUTINE set_vbc
!
!***********************************************************************
      SUBROUTINE set_vbc_tile (ng, Istr, Iend, Jstr, Jend,              &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         nrhs,                                    &
     &                         Hz, z_r, z_w,                            &
# if defined ICESHELF
     &                         zice,                                    &
# endif
     &                         t, u, v,                                 &
# ifdef QCORRECTION
     &                         dqdt, sst,                               &
# endif
# if defined SCORRECTION || defined SRELAXATION
     &                         sss,                                     &
# endif
# if defined ICESHELF
#  ifdef SHORTWAVE
     &                         srflx,                                   &
#  endif
     &                         sustr, svstr,                            &
# endif
# ifndef BBL
     &                         bustr, bvstr,                            &
# endif
     &                         stflx, btflx)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE bc_2d_mod
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: nrhs

# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
#  if defined ICESHELF
      real(r8), intent(in) :: zice(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: t(LBi:,LBj:,:,:,:)
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
#  ifdef QCORRECTION
      real(r8), intent(in) :: dqdt(LBi:,LBj:)
      real(r8), intent(in) :: sst(LBi:,LBj:)
#  endif
#  if defined SCORRECTION || defined SRELAXATION
      real(r8), intent(in) :: sss(LBi:,LBj:)
#  endif
#  if defined ICESHELF
#   ifdef SHORTWAVE
      real(r8), intent(inout) :: srflx(LBi:,LBj:)
#   endif
      real(r8), intent(inout) :: sustr(LBi:,LBj:)
      real(r8), intent(inout) :: svstr(LBi:,LBj:)
#  endif
#  ifndef BBL
      real(r8), intent(inout) :: bustr(LBi:,LBj:)
      real(r8), intent(inout) :: bvstr(LBi:,LBj:)
#  endif
      real(r8), intent(inout) :: stflx(LBi:,LBj:,:)
      real(r8), intent(inout) :: btflx(LBi:,LBj:,:)
# else
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
#  if defined ICESHELF
      real(r8), intent(in) :: zice(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),2)
#  ifdef QCORRECTION
      real(r8), intent(in) :: dqdt(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: sst(LBi:UBi,LBj:UBj)
#  endif
#  if defined SCORRECTION || defined SRELAXATION
      real(r8), intent(in) :: sss(LBi:UBi,LBj:UBj)
#  endif
#  if defined ICESHELF
#   ifdef SHORTWAVE
      real(r8), intent(inout) :: srflx(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(inout) :: sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: svstr(LBi:UBi,LBj:UBj)
#  endif
#  ifndef BBL
      real(r8), intent(inout) :: bustr(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: bvstr(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(inout) :: stflx(LBi:UBi,LBj:UBj,NT(ng))
      real(r8), intent(inout) :: btflx(LBi:UBi,LBj:UBj,NT(ng))
# endif

!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j, itrc

      real(r8) :: cff

      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: wrk

# include "set_bounds.h"

# ifdef QCORRECTION
!
!-----------------------------------------------------------------------
!  Add in flux correction to surface net heat flux (degC m/s).
!-----------------------------------------------------------------------
!
! Add in net heat flux correction.
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          stflx(i,j,itemp)=stflx(i,j,itemp)+                            &
     &                     dqdt(i,j)*(t(i,j,N(ng),nrhs,itemp)-sst(i,j))
        END DO
      END DO
# endif
# ifdef SALINITY
!
!-----------------------------------------------------------------------
!  Multiply fresh water flux with surface salinity. If appropriate,
!  apply correction.
!-----------------------------------------------------------------------
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
#  if defined SCORRECTION
          stflx(i,j,isalt)=stflx(i,j,isalt)*t(i,j,N(ng),nrhs,isalt)-    &
     &                     Tnudg(isalt,ng)*Hz(i,j,N(ng))*               &
     &                     (t(i,j,N(ng),nrhs,isalt)-sss(i,j))
#  elif defined SRELAXATION
          stflx(i,j,isalt)=-Tnudg(isalt,ng)*Hz(i,j,N(ng))*              &
     &                     (t(i,j,N(ng),nrhs,isalt)-sss(i,j))
#  else
          stflx(i,j,isalt)=stflx(i,j,isalt)*t(i,j,N(ng),nrhs,isalt)
#  endif
          btflx(i,j,isalt)=btflx(i,j,isalt)*t(i,j,1,nrhs,isalt)
        END DO
      END DO
# endif
# ifdef ICESHELF
!
!-----------------------------------------------------------------------
!  If ice shelf cavities, zero out for now the surface tracer flux
!  over the ice.
!-----------------------------------------------------------------------
!
      DO itrc=1,NT(ng)
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            IF (zice(i,j).ne.0.0_r8) THEN
              stflx(i,j,itrc)=0.0_r8
            END IF
          END DO
        END DO
      END DO
#  ifdef SHORTWAVE
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          IF (zice(i,j).ne.0.0_r8) THEN
            srflx(i,j)=0.0_r8
          END IF
        END DO
      END DO
#  endif
!
!-----------------------------------------------------------------------
!  If ice shelf cavities, replace surface wind stress with ice shelf
!  cavity stress (m2/s2).
!-----------------------------------------------------------------------
!
      IF (Zob(ng).ne.0.0_r8) THEN
        DO j=JstrV-1,Jend
          DO i=IstrU-1,Iend
            cff=vonKar/LOG((z_w(i,j,N(ng))-z_r(i,j,N(ng)))/Zob(ng))
            wrk(i,j)=MIN(Cdb_max,MAX(Cdb_min,cff*cff))
          END DO
        END DO
        DO j=Jstr,Jend
          DO i=IstrU,Iend
            IF (zice(i,j)*zice(i-1,j).ne.0.0_r8) THEN
              cff=0.25_r8*(v(i  ,j,N(ng),nrhs)+v(i  ,j+1,N(ng),nrhs)+   &
     &                     v(i-1,j,N(ng),nrhs)+v(i-1,j+1,N(ng),nrhs))
              sustr(i,j)=-0.5_r8*(wrk(i-1,j)+wrk(i,j))*                 &
     &                   u(i,j,N(ng),nrhs)*                             &
     &                   SQRT(u(i,j,N(ng),nrhs)*u(i,j,N(ng),nrhs)+      &
     &                        cff*cff)
            END IF
          END DO
        END DO
        DO j=JstrV,Jend
          DO i=Istr,Iend
            IF (zice(i,j)*zice(i,j-1).ne.0.0_r8) THEN
              cff=0.25_r8*(u(i,j  ,N(ng),nrhs)+u(i+1,j  ,N(ng),nrhs)+   &
     &                     u(i,j-1,N(ng),nrhs)+u(i+1,j-1,N(ng),nrhs))
              svstr(i,j)=-0.5_r8*(wrk(i,j-1)+wrk(i,j))*                 &
     &                   v(i,j,N(ng),nrhs)*                             &
     &                   SQRT(cff*cff+                                  &
     &                        v(i,j,N(ng),nrhs)*v(i,j,N(ng),nrhs))
            END IF
          END DO
        END DO
      ELSE IF (rdrg2(ng).ne.0.0_r8) THEN
        DO j=Jstr,Jend
          DO i=IstrU,Iend
            IF (zice(i,j)*zice(i-1,j).ne.0.0_r8) THEN
              cff=0.25_r8*(v(i  ,j,N(ng),nrhs)+v(i  ,j+1,N(ng),nrhs)+   &
     &                     v(i-1,j,N(ng),nrhs)+v(i-1,j+1,N(ng),nrhs))
              sustr(i,j)=-rdrg2(ng)*u(i,j,N(ng),nrhs)*                  &
     &                   SQRT(u(i,j,N(ng),nrhs)*u(i,j,N(ng),nrhs)+      &
     &                        cff*cff)
            END IF
          END DO
        END DO
        DO j=JstrV,Jend
          DO i=Istr,Iend
            IF (zice(i,j)*zice(i,j-1).ne.0.0_r8) THEN
              cff=0.25_r8*(u(i,j  ,N(ng),nrhs)+u(i+1,j  ,N(ng),nrhs)+   &
     &                     u(i,j-1,N(ng),nrhs)+u(i+1,j-1,N(ng),nrhs))
              svstr(i,j)=-rdrg2(ng)*v(i,j,N(ng),nrhs)*                  &
     &                   SQRT(cff*cff+                                  &
     &                        v(i,j,N(ng),nrhs)*v(i,j,N(ng),nrhs))
            END IF
          END DO
        END DO
      ELSE IF (rdrg(ng).ne.0.0_r8) THEN
        DO j=Jstr,Jend
          DO i=IstrU,Iend
            IF (zice(i,j)*zice(i-1,j).ne.0.0_r8) THEN
              sustr(i,j)=-rdrg(ng)*u(i,j,N(ng),nrhs)
            END IF
          END DO
        END DO
        DO j=JstrV,Jend
          DO i=Istr,Iend
            IF (zice(i,j)*zice(i,j-1).ne.0.0_r8) THEN
              svstr(i,j)=-rdrg(ng)*v(i,j,N(ng),nrhs)
            END IF
          END DO
        END DO
      ELSE
        DO j=Jstr,Jend
          DO i=IstrU,Iend
            IF (zice(i,j)*zice(i-1,j).ne.0.0_r8) THEN
              sustr(i,j)=0.0_r8
            END IF
          END DO
        END DO
        DO j=JstrV,Jend
          DO i=Istr,Iend
            IF (zice(i,j)*zice(i,j-1).ne.0.0_r8) THEN
              svstr(i,j)=0.0_r8
            END IF
          END DO
        END DO
      END IF
!
!  Apply periodic or gradient boundary conditions for output
!  purposes only.
!
      CALL bc_u2d_tile (ng, Istr, Iend, Jstr, Jend,                     &
     &                  LBi, UBi, LBj, UBj, sustr)
      CALL bc_v2d_tile (ng, Istr, Iend, Jstr, Jend,                     &
     &                  LBi, UBi, LBj, UBj, svstr)
# endif
# ifndef BBL
!
!-----------------------------------------------------------------------
!  Set kinematic bottom momentum flux (m2/s2).
!-----------------------------------------------------------------------
!
!  Set bottom stress using linear and/or quadratic formulation.
!
      IF (Zob(ng).ne.0.0_r8) THEN
        DO j=JstrV-1,Jend
          DO i=IstrU-1,Iend
            cff=vonKar/LOG((z_r(i,j,1)-z_w(i,j,0))/Zob(ng))
            wrk(i,j)=MIN(Cdb_max,MAX(Cdb_min,cff*cff))
          END DO
        END DO
        DO j=Jstr,Jend
          DO i=IstrU,Iend
            cff=0.25_r8*(v(i  ,j,1,nrhs)+v(i  ,j+1,1,nrhs)+             &
     &                   v(i-1,j,1,nrhs)+v(i-1,j+1,1,nrhs))
            bustr(i,j)=0.5_r8*(wrk(i-1,j)+wrk(i,j))*u(i,j,1,nrhs)*      &
     &                 SQRT(u(i,j,1,nrhs)*u(i,j,1,nrhs)+cff*cff)
          END DO
        END DO
        DO j=JstrV,Jend
          DO i=Istr,Iend
            cff=0.25_r8*(u(i,j  ,1,nrhs)+u(i+1,j,1,nrhs)+               &
     &                   u(i,j-1,1,nrhs)+u(i+1,j-1,1,nrhs))
            bvstr(i,j)=0.5_r8*(wrk(i,j-1)+wrk(i,j))*v(i,j,1,nrhs)*      &
     &                 SQRT(cff*cff+v(i,j,1,nrhs)*v(i,j,1,nrhs))
          END DO
        END DO
      ELSE IF (rdrg2(ng).ne.0.0_r8) THEN
        DO j=Jstr,Jend
          DO i=IstrU,Iend
            cff=0.25_r8*(v(i  ,j,1,nrhs)+v(i  ,j+1,1,nrhs)+             &
     &                   v(i-1,j,1,nrhs)+v(i-1,j+1,1,nrhs))
            bustr(i,j)=rdrg2(ng)*u(i,j,1,nrhs)*                         &
     &                 SQRT(u(i,j,1,nrhs)*u(i,j,1,nrhs)+cff*cff)
          END DO
        END DO
        DO j=JstrV,Jend
          DO i=Istr,Iend
            cff=0.25_r8*(u(i,j  ,1,nrhs)+u(i+1,j,1,nrhs)+               &
     &                   u(i,j-1,1,nrhs)+u(i+1,j-1,1,nrhs))
            bvstr(i,j)=rdrg2(ng)*v(i,j,1,nrhs)*                         &
     &                 SQRT(cff*cff+v(i,j,1,nrhs)*v(i,j,1,nrhs))
          END DO
        END DO
      ELSE IF (rdrg(ng).ne.0.0_r8) THEN
        DO j=Jstr,Jend
          DO i=IstrU,Iend
            bustr(i,j)=rdrg(ng)*u(i,j,1,nrhs)
          END DO
        END DO
        DO j=JstrV,Jend
          DO i=Istr,Iend
            bvstr(i,j)=rdrg(ng)*v(i,j,1,nrhs)
          END DO
        END DO
      END IF
!
!  Apply boundary conditions.
!
      CALL bc_u2d_tile (ng, Istr, Iend, Jstr, Jend,                     &
     &                  LBi, UBi, LBj, UBj, bustr)
      CALL bc_v2d_tile (ng, Istr, Iend, Jstr, Jend,                     &
     &                  LBi, UBi, LBj, UBj, bvstr)
# endif
      RETURN
      END SUBROUTINE set_vbc_tile
#endif
      END MODULE set_vbc_mod
