#include "cppdefs.h"
#ifdef AVERAGES
      SUBROUTINE set_avg (ng, tile)
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This subroutine accumulates and computes output time-averaged       !
!  fields.  Due to synchronization, the time-averaged fields are       !
!  computed in delayed mode. All averages are accumulated at the       !
!  beggining of the next time-step.                                    !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_scalars
      USE mod_stepping
!
      implicit none

      integer, intent(in) :: ng, tile

# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, 5)
# endif
      CALL set_avg_tile (ng, Istr, Iend, Jstr, Jend,                    &
     &                   KOUT, NOUT)
# ifdef PROFILE
      CALL wclock_off (ng, 5)
# endif
      RETURN
      END SUBROUTINE set_avg
!
!***********************************************************************
      SUBROUTINE set_avg_tile (ng, Istr, Iend, Jstr, Jend,              &
     &                         Kout, Nout)
!***********************************************************************
!
      USE mod_param
      USE mod_average
      USE mod_grid
      USE mod_mixing
      USE mod_ocean
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: Kout, Nout
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, it, j, k

      real(r8) :: fac

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
! Initialize time-averaged arrays when appropriate.  Notice that
! fields are initilized twice during re-start.  However, the time-
! averaged fields are computed correctly.
!-----------------------------------------------------------------------
!
      IF (((iic(ng).gt.ntsavg(ng)).and.                                 &
     &     (MOD(iic(ng)-1,navg(ng)).eq.1)).or.                          &
     &    ((nrrec.gt.0).and.(iic(ng).eq.ntstart))) THEN
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            AVERAGE(ng)%avgzeta(i,j)=OCEAN(ng)%zeta(i,j,Kout)
            AVERAGE(ng)%avgu2d (i,j)=OCEAN(ng)%ubar(i,j,Kout)
            AVERAGE(ng)%avgv2d (i,j)=OCEAN(ng)%vbar(i,j,Kout)
          END DO
        END DO
# ifdef SOLVE3D
#  ifdef LMD_SKPP
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            AVERAGE(ng)%avghsbl(i,j)=MIXING(ng)%hsbl(i,j)
          END DO
        END DO
#  endif
#  ifdef LMD_BKPP
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            AVERAGE(ng)%avghbbl(i,j)=MIXING(ng)%hbbl(i,j)
          END DO
        END DO
#  endif
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgu3d(i,j,k)=OCEAN(ng)%u(i,j,k,Nout)
              AVERAGE(ng)%avgv3d(i,j,k)=OCEAN(ng)%v(i,j,k,Nout)
              AVERAGE(ng)%avgrho(i,j,k)=OCEAN(ng)%rho(i,j,k)
            END DO
          END DO
          DO it=1,NT(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgt(i,j,k,it)=OCEAN(ng)%t(i,j,k,Nout,it)
              END DO
            END DO
          END DO
        END DO
        DO k=0,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgw3d(i,j,k)=OCEAN(ng)%W(i,j,k)*             &
     &                                  GRID(ng)%pm(i,j)*               &
     &                                  GRID(ng)%pn(i,j)
            END DO
          END DO
        END DO
#  ifdef AVERAGES_AKV
        DO k=0,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgAKv(i,j,k)=MIXING(ng)%Akv(i,j,k)
            END DO
          END DO
        END DO
#  endif
#  ifdef AVERAGES_AKT
        DO k=0,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgAKt(i,j,k)=MIXING(ng)%Akt(i,j,k,itemp)
            END DO
          END DO
        END DO
#  endif
#  ifdef AVERAGES_AKS
        DO k=0,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgAKs(i,j,k)=MIXING(ng)%Akt(i,j,k,isalt)
            END DO
          END DO
        END DO
#  endif
# endif
!
!-----------------------------------------------------------------------
!  Accumulate time-averaged fields.
!-----------------------------------------------------------------------
!
      ELSE IF (iic(ng).gt.ntsavg(ng)) THEN
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            AVERAGE(ng)%avgzeta(i,j)=AVERAGE(ng)%avgzeta(i,j)+          &
     &                               OCEAN(ng)%zeta(i,j,Kout)
            AVERAGE(ng)%avgu2d (i,j)=AVERAGE(ng)%avgu2d (i,j)+          &
     &                               OCEAN(ng)%ubar(i,j,Kout)
            AVERAGE(ng)%avgv2d (i,j)=AVERAGE(ng)%avgv2d (i,j)+          &
     &                               OCEAN(ng)%vbar(i,j,Kout)
          END DO
        END DO
# ifdef SOLVE3D
#  ifdef LMD_SKPP
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            AVERAGE(ng)%avghsbl(i,j)=AVERAGE(ng)%avghsbl(i,j)+          &
     &                               MIXING(ng)%hsbl(i,j)
          END DO
        END DO
#  endif
#  ifdef LMD_BKPP
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            AVERAGE(ng)%avghbbl(i,j)=AVERAGE(ng)%avghbbl(i,j)+          &
     &                               MIXING(ng)%hbbl(i,j)
          END DO
        END DO
#  endif
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgu3d(i,j,k)=AVERAGE(ng)%avgu3d(i,j,k)+      &
     &                                  OCEAN(ng)%u(i,j,k,Nout)
              AVERAGE(ng)%avgv3d(i,j,k)=AVERAGE(ng)%avgv3d(i,j,k)+      &
     &                                  OCEAN(ng)%v(i,j,k,Nout)
              AVERAGE(ng)%avgrho(i,j,k)=AVERAGE(ng)%avgrho(i,j,k)+      &
     &                                  OCEAN(ng)%rho(i,j,k)
            END DO
          END DO
          DO it=1,NT(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgt(i,j,k,it)=AVERAGE(ng)%avgt(i,j,k,it)+  &
     &                                     OCEAN(ng)%t(i,j,k,Nout,it)
              END DO
            END DO
          END DO
        END DO
        DO k=0,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgw3d(i,j,k)=AVERAGE(ng)%avgw3d(i,j,k)+      &
     &                                  OCEAN(ng)%W(i,j,k)*             &
     &                                  GRID(ng)%pm(i,j)*               &
     &                                  GRID(ng)%pn(i,j)
            END DO
          END DO
        END DO
#  ifdef AVERAGES_AKV
        DO k=0,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgAKv(i,j,k)=AVERAGE(ng)%avgAKv(i,j,k)+      &
     &                                  MIXING(ng)%Akv(i,j,k)
            END DO
          END DO
        END DO
#  endif
#  ifdef AVERAGES_AKT
        DO k=0,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgAKt(i,j,k)=AVERAGE(ng)%avgAKt(i,j,k)+      &
     &                                  MIXING(ng)%Akt(i,j,k,itemp)
            END DO
          END DO
        END DO
#  endif
#  ifdef AVERAGES_AKS
        DO k=0,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgAKs(i,j,k)=AVERAGE(ng)%avgAKs(i,j,k)+      &
     &                                  MIXING(ng)%Akt(i,j,k,isalt)
            END DO
          END DO
        END DO
#  endif
# endif
      END IF
!
!-----------------------------------------------------------------------
!  Convert accumulated sums into time-averages, if appropriate.
!-----------------------------------------------------------------------
!
      IF ((iic(ng).gt.ntsavg(ng)).and.                                  &
     &    (MOD(iic(ng)-1,navg(ng)).eq.0).and.                           &
     &    ((iic(ng).ne.ntstart).or.(nrrec.eq.0))) THEN
        fac=1.0_r8/REAL(navg(ng),r8)
        IF (SOUTH_WEST_TEST) THEN
          avgtime(ng)=avgtime(ng)+REAL(navg(ng),r8)*dt(ng)
        END IF
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            AVERAGE(ng)%avgzeta(i,j)=fac*AVERAGE(ng)%avgzeta(i,j)
            AVERAGE(ng)%avgu2d (i,j)=fac*AVERAGE(ng)%avgu2d (i,j)
            AVERAGE(ng)%avgv2d (i,j)=fac*AVERAGE(ng)%avgv2d (i,j)
          END DO
        END DO
# ifdef SOLVE3D
#  ifdef LMD_SKPP
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            AVERAGE(ng)%avghsbl(i,j)=fac*AVERAGE(ng)%avghsbl(i,j)
          END DO
        END DO
#  endif
#  ifdef LMD_BKPP
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            AVERAGE(ng)%avghbbl(i,j)=fac*AVERAGE(ng)%avghbbl(i,j)
          END DO
        END DO
#  endif
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgu3d(i,j,k)=fac*AVERAGE(ng)%avgu3d(i,j,k)
              AVERAGE(ng)%avgv3d(i,j,k)=fac*AVERAGE(ng)%avgv3d(i,j,k)
              AVERAGE(ng)%avgrho(i,j,k)=fac*AVERAGE(ng)%avgrho(i,j,k)
            END DO
          END DO
          DO it=1,NT(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgt(i,j,k,it)=fac*                         &
     &                                     AVERAGE(ng)%avgt(i,j,k,it)
              END DO
            END DO
          END DO
        END DO
        DO k=0,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgw3d(i,j,k)=fac*AVERAGE(ng)%avgw3d(i,j,k)
            END DO
          END DO
        END DO
#  ifdef AVERAGES_AKV
        DO k=0,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgAKv(i,j,k)=fac*AVERAGE(ng)%avgAKv(i,j,k)
            END DO
          END DO
        END DO
#  endif
#  ifdef AVERAGES_AKT
        DO k=0,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgAKt(i,j,k)=fac*AVERAGE(ng)%avgAKt(i,j,k)
            END DO
          END DO
        END DO
#  endif
#  ifdef AVERAGES_AKS
        DO k=0,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgAKs(i,j,k)=fac*AVERAGE(ng)%avgAKs(i,j,k)
            END DO
          END DO
        END DO
#  endif
# endif
      END IF
      RETURN
      END SUBROUTINE set_avg_tile
#else
      SUBROUTINE set_avg
      RETURN
      END SUBROUTINE set_avg
#endif
