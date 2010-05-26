#include "cppdefs.h"
      MODULE ad_set_avg_mod
#if defined AVERAGES && defined ADJOINT
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine accumulates and computes output time-averaged       !
!  adjoint fields.                                                     !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC :: ad_set_avg

      CONTAINS
!
!***********************************************************************
      SUBROUTINE ad_set_avg (ng, tile)
!***********************************************************************
!
      USE mod_param
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
      CALL wclock_on (ng, iADM, 5)
# endif
      CALL ad_set_avg_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      nstp(ng))
# ifdef PROFILE
      CALL wclock_off (ng, iADM, 5)
# endif
      RETURN
      END SUBROUTINE ad_set_avg
!
!***********************************************************************
      SUBROUTINE ad_set_avg_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            nstp)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_average
      USE mod_forces
# ifdef SOLVE3D
      USE mod_grid
      USE mod_mixing
# endif
      USE mod_ocean
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp
!
!
!  Local variable declarations.
!
      integer :: i, it, j, k

      real(r8) :: fac

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Return if time-averaging window is zero.
!-----------------------------------------------------------------------
!
      IF (nAVG(ng).eq.0) RETURN
!
!-----------------------------------------------------------------------
!  Initialize time-averaged arrays when appropriate.  Notice that
!  fields are initilized twice during re-start.  However, the time-
!  averaged fields are computed correctly.
!-----------------------------------------------------------------------
!
      IF (((iic(ng).lt.ntsAVG(ng)).and.                                 &
     &     (MOD(iic(ng),nAVG(ng)).eq.0)).or.                            &
     &    ((nrrec(ng).gt.0).and.(iic(ng).eq.ntstart(ng)))) THEN
!
!  Initialize adjoint state variables.
!
        IF (Aout(idFsur,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgzeta(i,j)=OCEAN(ng)%ad_zeta_sol(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idUbar,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgu2d(i,j)=OCEAN(ng)%ad_ubar_sol(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idVbar,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgv2d(i,j)=OCEAN(ng)%ad_vbar_sol(i,j)
            END DO
          END DO
        END IF

# ifdef SOLVE3D
        IF (Aout(idUvel,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgu3d(i,j,k)=OCEAN(ng)%ad_u(i,j,k,nstp)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idVvel,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgv3d(i,j,k)=OCEAN(ng)%ad_v(i,j,k,nstp)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idOvel,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgw3d(i,j,k)=OCEAN(ng)%ad_W(i,j,k)*        &
     &                                    GRID(ng)%pm(i,j)*             &
     &                                    GRID(ng)%pn(i,j)
              END DO
            END DO
          END DO
        END IF

        IF (Aout(idDano,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgrho(i,j,k)=OCEAN(ng)%ad_rho(i,j,k)
              END DO
            END DO
          END DO
        END IF
        DO it=1,NT(ng)
          IF (Aout(idTvar(it),ng)) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgt(i,j,k,it)=OCEAN(ng)%ad_t(i,j,k,      &
     &                                                      nstp,it)
                END DO
              END DO
            END DO
          END IF
        END DO

#  if defined LMD_MIXING || defined MY25_MIXING || defined GLS_MIXING
        IF (Aout(idVvis,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgAKv(i,j,k)=MIXING(ng)%ad_Akv(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idTdif,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgAKt(i,j,k)=MIXING(ng)%ad_Akt(i,j,k,itemp)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idSdif,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgAKs(i,j,k)=MIXING(ng)%ad_Akt(i,j,k,isalt)
              END DO
            END DO
          END DO
        END IF
#  endif
# endif
!
!  Initialize adjoint surface and bottom fluxes.
!
        IF (Aout(idUsms,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgsus(i,j)=FORCES(ng)%ad_sustr(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idVsms,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgsvs(i,j)=FORCES(ng)%ad_svstr(i,j)
            END DO
          END DO
        END IF

        IF (Aout(idUbms,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgbus(i,j)=FORCES(ng)%ad_bustr(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idVbms,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgbvs(i,j)=FORCES(ng)%ad_bvstr(i,j)
            END DO
          END DO
        END IF

# ifdef SOLVE3D
        IF (Aout(idTsur(itemp),ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgstf(i,j)=FORCES(ng)%ad_stflx(i,j,itemp)
            END DO
          END DO
        END IF
        IF (Aout(idTsur(isalt),ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgswf(i,j)=FORCES(ng)%ad_stflx(i,j,isalt)
            END DO
          END DO
        END IF
#  ifdef SHORTWAVE
        IF (Aout(idSrad,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgsrf(i,j)=FORCES(ng)%ad_srflx(i,j)
            END DO
          END DO
        END IF
#  endif
#  ifdef BULK_FLUXES
        IF (Aout(idLhea,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avglhf(i,j)=FORCES(ng)%ad_lhflx(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idLrad,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avglrf(i,j)=FORCES(ng)%ad_lrflx(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idShea,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgshf(i,j)=FORCES(ng)%ad_shflx(i,j)
            END DO
          END DO
        END IF
#   ifdef EMINUSP
        IF (Aout(idevap,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgevap(i,j)=FORCES(ng)%ad_evap(i,j)
            END DO
          END DO
        END IF
#   endif
#  endif
# endif
!
!  Initialize adjoint of quadratic fields.
!
        IF (Aout(idZZav,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgZZ(i,j)=OCEAN(ng)%ad_zeta_sol(i,j)*        &
     &                               OCEAN(ng)%ad_zeta_sol(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idU2av,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgU2(i,j)=OCEAN(ng)%ad_ubar_sol(i,j)*        &
     &                               OCEAN(ng)%ad_ubar_sol(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idV2av,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgV2(i,j)=OCEAN(ng)%ad_vbar_sol(i,j)*        &
     &                               OCEAN(ng)%ad_vbar_sol(i,j)
            END DO
          END DO
        END IF

# ifdef SOLVE3D
        IF (Aout(idUUav,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgUU(i,j,k)=OCEAN(ng)%ad_u(i,j,k,nstp)*    &
     &                                   OCEAN(ng)%ad_u(i,j,k,nstp)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idVVav,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgVV(i,j,k)=OCEAN(ng)%ad_v(i,j,k,nstp)*     &
     &                                   OCEAN(ng)%ad_v(i,j,k,nstp)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idUVav,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              DO i=Istr,Iend
                AVERAGE(ng)%avgUV(i,j,k)=0.25_r8*                       &
     &                                   (OCEAN(ng)%ad_u(i  ,j,k,nstp)+ &
     &                                    OCEAN(ng)%ad_u(i+1,j,k,nstp))*&
     &                                   (OCEAN(ng)%ad_v(i,j  ,k,nstp)+ &
     &                                    OCEAN(ng)%ad_v(i,j+1,k,nstp))
              END DO
            END DO
          END DO
        END IF

        DO it=1,NAT
          IF (Aout(idTTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgTT(i,j,k,it)=OCEAN(ng)%ad_t(i,j,k,     &
     &                                                       nstp,it)*  &
     &                                        OCEAN(ng)%ad_t(i,j,k,     &
     &                                                       nstp,it)
                END DO
              END DO
            END DO
          END IF
          IF (Aout(idUTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=Istr,Iend
                  AVERAGE(ng)%avgUT(i,j,k,it)=0.5_r8*                   &
     &                                        OCEAN(ng)%ad_u(i,j,k,     &
     &                                                       nstp)*     &
     &                                        (OCEAN(ng)%ad_t(i-1,j,k,  &
     &                                                        nstp,it)+ &
     &                                         OCEAN(ng)%ad_t(i  ,j,k,  &
     &                                                        nstp,it))
                END DO
              END DO
            END DO
          END IF
          IF (Aout(idVTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=Jstr,Jend
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgVT(i,j,k,it)=0.5_r8*                   &
     &                                        OCEAN(ng)%ad_v(i,j,k,     &
     &                                                       nstp)*     &
     &                                        (OCEAN(ng)%ad_t(i,j-1,k,  &
     &                                                        nstp,it)+ &
     &                                         OCEAN(ng)%ad_t(i,j  ,k,  &
     &                                                        nstp,it))
                END DO
              END DO
            END DO
          END IF
        END DO
# endif
!
!-----------------------------------------------------------------------
!  Accumulate time-averaged fields.
!-----------------------------------------------------------------------
!
      ELSE IF ((iic(ng)-1).le.ntsAVG(ng)) THEN
!
!  Accumulate adjoint state variables.
!
        IF (Aout(idFsur,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgzeta(i,j)=AVERAGE(ng)%avgzeta(i,j)+        &
     &                                 OCEAN(ng)%ad_zeta_sol(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idUbar,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgu2d(i,j)=AVERAGE(ng)%avgu2d(i,j)+          &
     &                                OCEAN(ng)%ad_ubar_sol(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idVbar,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgv2d(i,j)=AVERAGE(ng)%avgv2d(i,j)+          &
     &                                OCEAN(ng)%ad_vbar_sol(i,j)
            END DO
          END DO
        END IF

# ifdef SOLVE3D
        IF (Aout(idUvel,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgu3d(i,j,k)=AVERAGE(ng)%avgu3d(i,j,k)+    &
     &                                    OCEAN(ng)%ad_u(i,j,k,nstp)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idVvel,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgv3d(i,j,k)=AVERAGE(ng)%avgv3d(i,j,k)+    &
     &                                    OCEAN(ng)%ad_v(i,j,k,nstp)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idOvel,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgw3d(i,j,k)=AVERAGE(ng)%avgw3d(i,j,k)+    &
     &                                    OCEAN(ng)%ad_W(i,j,k)*        &
     &                                    GRID(ng)%pm(i,j)*             &
     &                                    GRID(ng)%pn(i,j)
              END DO
            END DO
          END DO
        END IF

        IF (Aout(idDano,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgrho(i,j,k)=AVERAGE(ng)%avgrho(i,j,k)+    &
     &                                    OCEAN(ng)%ad_rho(i,j,k)
              END DO
            END DO
          END DO
        END IF
        DO it=1,NT(ng)
          IF (Aout(idTvar(it),ng)) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgt(i,j,k,it)=AVERAGE(ng)%avgt(i,j,k,it)+&
     &                                       OCEAN(ng)%ad_t(i,j,k,      &
     &                                                      nstp,it)
                END DO
              END DO
            END DO
          END IF
        END DO

#  if defined LMD_MIXING || defined MY25_MIXING || defined GLS_MIXING
        IF (Aout(idVvis,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgAKv(i,j,k)=AVERAGE(ng)%avgAKv(i,j,k)+    &
     &                                    MIXING(ng)%ad_Akv(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idTdif,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgAKt(i,j,k)=AVERAGE(ng)%avgAKt(i,j,k)+    &
     &                                    MIXING(ng)%ad_Akt(i,j,k,itemp)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idSdif,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgAKs(i,j,k)=AVERAGE(ng)%avgAKs(i,j,k)+    &
     &                                    MIXING(ng)%ad_Akt(i,j,k,isalt)
              END DO
            END DO
          END DO
        END IF
#  endif
# endif
!
!  Accumulate adjoint surface and bottom fluxes.
!
        IF (Aout(idUsms,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgsus(i,j)=AVERAGE(ng)%avgsus(i,j)+          &
     &                                FORCES(ng)%ad_sustr(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idVsms,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgsvs(i,j)=AVERAGE(ng)%avgsvs(i,j)+          &
     &                                FORCES(ng)%ad_svstr(i,j)
            END DO
          END DO
        END IF

        IF (Aout(idUbms,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgbus(i,j)=AVERAGE(ng)%avgbus(i,j)+          &
     &                                FORCES(ng)%ad_bustr(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idVbms,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgbvs(i,j)=AVERAGE(ng)%avgbvs(i,j)+          &
     &                                FORCES(ng)%ad_bvstr(i,j)
            END DO
          END DO
        END IF

# ifdef SOLVE3D
        IF (Aout(idTsur(itemp),ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgstf(i,j)=AVERAGE(ng)%avgstf(i,j)+          &
     &                                FORCES(ng)%ad_stflx(i,j,itemp)
            END DO
          END DO
        END IF
        IF (Aout(idTsur(isalt),ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgswf(i,j)=AVERAGE(ng)%avgswf(i,j)+          &
     &                                FORCES(ng)%ad_stflx(i,j,isalt)
            END DO
          END DO
        END IF
#  ifdef SHORTWAVE
        IF (Aout(idSrad,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgsrf(i,j)=AVERAGE(ng)%avgsrf(i,j)+          &
     &                                FORCES(ng)%ad_srflx(i,j)
            END DO
          END DO
        END IF
#  endif
#  ifdef BULK_FLUXES
        IF (Aout(idLhea,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avglhf(i,j)=AVERAGE(ng)%avglhf(i,j)+          &
     &                                FORCES(ng)%ad_lhflx(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idLrad,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avglrf(i,j)=AVERAGE(ng)%avglrf(i,j)+          &
     &                                FORCES(ng)%ad_lrflx(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idShea,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgshf(i,j)=AVERAGE(ng)%avgshf(i,j)+          &
     &                                FORCES(ng)%ad_shflx(i,j)
            END DO
          END DO
        END IF
#   ifdef EMINUSP
        IF (Aout(idevap,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgevap(i,j)=AVERAGE(ng)%avgevap(i,j)+        &
     &                                 FORCES(ng)%ad_evap(i,j)
            END DO
          END DO
        END IF
#   endif
#  endif
# endif
!
!  Accumulate adjoint quadratic fields.
!
        IF (Aout(idZZav,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgZZ(i,j)=AVERAGE(ng)%avgZZ(i,j)+            &
     &                               OCEAN(ng)%ad_zeta_sol(i,j)*        &
     &                               OCEAN(ng)%ad_zeta_sol(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idU2av,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgU2(i,j)=AVERAGE(ng)%avgU2(i,j)+            &
     &                               OCEAN(ng)%ad_ubar_sol(i,j)*        &
     &                               OCEAN(ng)%ad_ubar_sol(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idV2av,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgV2(i,j)=AVERAGE(ng)%avgV2(i,j)+            &
     &                               OCEAN(ng)%ad_vbar_sol(i,j)*        &
     &                               OCEAN(ng)%ad_vbar_sol(i,j)
            END DO
          END DO
        END IF

# ifdef SOLVE3D
        IF (Aout(idUUav,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgUU(i,j,k)=AVERAGE(ng)%avgUU(i,j,k)+      &
     &                                   OCEAN(ng)%ad_u(i,j,k,nstp)*    &
     &                                   OCEAN(ng)%ad_u(i,j,k,nstp)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idVVav,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgVV(i,j,k)=AVERAGE(ng)%avgVV(i,j,k)+       &
     &                                   OCEAN(ng)%ad_v(i,j,k,nstp)*     &
     &                                   OCEAN(ng)%ad_v(i,j,k,nstp)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idUVav,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              DO i=Istr,Iend
                AVERAGE(ng)%avgUV(i,j,k)=AVERAGE(ng)%avgUV(i,j,k)+      &
     &                                   0.25_r8*                       &
     &                                   (OCEAN(ng)%ad_u(i  ,j,k,nstp)+ &
     &                                    OCEAN(ng)%ad_u(i+1,j,k,nstp))*&
     &                                   (OCEAN(ng)%ad_v(i,j  ,k,nstp)+ &
     &                                    OCEAN(ng)%ad_v(i,j+1,k,nstp))
              END DO
            END DO
          END DO
        END IF

        DO it=1,NAT
          IF (Aout(idTTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgTT(i,j,k,it)=AVERAGE(ng)%avgTT(i,j,k,  &
     &                                                          it)+    &
     &                                        OCEAN(ng)%ad_t(i,j,k,     &
     &                                                       nstp,it)*  &
     &                                        OCEAN(ng)%ad_t(i,j,k,     &
     &                                                       nstp,it)
                END DO
              END DO
            END DO
          END IF
          IF (Aout(idUTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=Istr,Iend
                  AVERAGE(ng)%avgUT(i,j,k,it)=AVERAGE(ng)%avgUT(i,j,k,  &
     &                                                          it)+    &
     &                                        0.5_r8*                   &
     &                                        OCEAN(ng)%ad_u(i,j,k,     &
     &                                                       nstp)*     &
     &                                        (OCEAN(ng)%ad_t(i-1,j,k,  &
     &                                                        nstp,it)+ &
     &                                         OCEAN(ng)%ad_t(i  ,j,k,  &
     &                                                        nstp,it))
                END DO
              END DO
            END DO
          END IF
          IF (Aout(idVTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=Jstr,Jend
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgVT(i,j,k,it)=AVERAGE(ng)%avgVT(i,j,k,  &
     &                                                          it)+    &
     &                                        0.5_r8*                   &
     &                                        OCEAN(ng)%ad_v(i,j,k,     &
     &                                                       nstp)*     &
     &                                        (OCEAN(ng)%ad_t(i,j-1,k,  &
     &                                                        nstp,it)+ &
     &                                         OCEAN(ng)%ad_t(i,j  ,k,  &
     &                                                        nstp,it))
                END DO
              END DO
            END DO
          END IF
        END DO
# endif
      END IF
!
!-----------------------------------------------------------------------
!  Convert accumulated sums into time-averages, if appropriate.
!-----------------------------------------------------------------------
!
      IF ((iic(ng).lt.ntsAVG(ng)).and.                                  &
     &    (MOD(iic(ng)-1,nAVG(ng)).eq.0).and.                           &
     &    ((iic(ng).ne.ntstart(ng)).or.(nrrec(ng).eq.0))) THEN
# if defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
     defined OPT_OBSERVATIONS || defined SENSITIVITY_4DVAR   || \
     defined SO_SEMI
        fac=1.0_r8
# else
        fac=1.0_r8/REAL(nAVG(ng),r8)
# endif
        IF (SOUTH_WEST_TEST) THEN
          AVGtime(ng)=AVGtime(ng)-REAL(nAVG(ng),r8)*dt(ng)
        END IF
!
!  Process adjoint state variables.
!
        IF (Aout(idFsur,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgzeta(i,j)=fac*AVERAGE(ng)%avgzeta(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idUbar,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgu2d(i,j)=fac*AVERAGE(ng)%avgu2d(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idVbar,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgv2d(i,j)=fac*AVERAGE(ng)%avgv2d(i,j)
            END DO
          END DO
        END IF

# ifdef SOLVE3D
        IF (Aout(idUvel,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgu3d(i,j,k)=fac*AVERAGE(ng)%avgu3d(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idVvel,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgv3d(i,j,k)=fac*AVERAGE(ng)%avgv3d(i,j,k)
              END DO
            END DO
          END DO
        END IF

        IF (Aout(idOvel,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgw3d(i,j,k)=fac*AVERAGE(ng)%avgw3d(i,j,k)
              END DO
            END DO
          END DO
        END IF

        IF (Aout(idDano,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgrho(i,j,k)=fac*AVERAGE(ng)%avgrho(i,j,k)
              END DO
            END DO
          END DO
        END IF
        DO it=1,NT(ng)
          IF (Aout(idTvar(it),ng)) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgt(i,j,k,it)=fac*                       &
     &                                       AVERAGE(ng)%avgt(i,j,k,it)
                END DO
              END DO
            END DO
          END IF
        END DO

#  if defined LMD_MIXING || defined MY25_MIXING || defined GLS_MIXING
        IF (Aout(idVvis,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgAKv(i,j,k)=fac*AVERAGE(ng)%avgAKv(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idTdif,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgAKt(i,j,k)=fac*AVERAGE(ng)%avgAKt(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idSdif,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgAKs(i,j,k)=fac*AVERAGE(ng)%avgAKs(i,j,k)
              END DO
            END DO
          END DO
        END IF
#  endif
# endif
!
!  Process adjoint surface and bottom fluxes.
!
        IF (Aout(idUsms,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgsus(i,j)=fac*AVERAGE(ng)%avgsus(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idVsms,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgsvs(i,j)=fac*AVERAGE(ng)%avgsvs(i,j)
            END DO
          END DO
        END IF

        IF (Aout(idUbms,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgbus(i,j)=fac*AVERAGE(ng)%avgbus(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idVbms,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgbvs(i,j)=fac*AVERAGE(ng)%avgbvs(i,j)
            END DO
          END DO
        END IF

# ifdef SOLVE3D
        IF (Aout(idTsur(itemp),ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgstf(i,j)=fac*AVERAGE(ng)%avgstf(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idTsur(isalt),ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgswf(i,j)=fac*AVERAGE(ng)%avgswf(i,j)
            END DO
          END DO
        END IF
#  ifdef SHORTWAVE
        IF (Aout(idSrad,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgsrf(i,j)=fac*AVERAGE(ng)%avgsrf(i,j)
            END DO
          END DO
        END IF
#  endif
#  ifdef BULK_FLUXES
        IF (Aout(idLhea,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avglhf(i,j)=fac*AVERAGE(ng)%avglhf(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idShea,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgshf(i,j)=fac*AVERAGE(ng)%avgshf(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idLrad,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avglrf(i,j)=fac*AVERAGE(ng)%avglrf(i,j)
            END DO
          END DO
        END IF
#   ifdef EMINUSP
        IF (Aout(idevap,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgevap(i,j)=fac*AVERAGE(ng)%avgevap(i,j)
            END DO
          END DO
        END IF
#   endif
#  endif
# endif
!
!  Process adjoint quadratic fields.
!
        IF (Aout(idZZav,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgZZ(i,j)=fac*AVERAGE(ng)%avgZZ(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idU2av,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgU2(i,j)=fac*AVERAGE(ng)%avgU2(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idV2av,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgV2(i,j)=fac*AVERAGE(ng)%avgV2(i,j)
            END DO
          END DO
        END IF

# ifdef SOLVE3D
        IF (Aout(idUUav,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgUU(i,j,k)=fac*AVERAGE(ng)%avgUU(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idVVav,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgVV(i,j,k)=fac*AVERAGE(ng)%avgVV(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idUVav,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              DO i=Istr,Iend
                AVERAGE(ng)%avgUV(i,j,k)=fac*AVERAGE(ng)%avgUV(i,j,k)
              END DO
            END DO
          END DO
        END IF

        DO it=1,NAT
          IF (Aout(idTTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgTT(i,j,k,it)=fac*                      &
     &                                       AVERAGE(ng)%avgTT(i,j,k,it)
                END DO
              END DO
            END DO
          END IF
          IF (Aout(idUTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=Istr,Iend
                  AVERAGE(ng)%avgUT(i,j,k,it)=fac*                      &
     &                                       AVERAGE(ng)%avgUT(i,j,k,it)
                END DO
              END DO
            END DO
          END IF
          IF (Aout(idVTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=Jstr,Jend
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgVT(i,j,k,it)=fac*                      &
     &                                       AVERAGE(ng)%avgVT(i,j,k,it)
                END DO
              END DO
            END DO
          END IF
        END DO
# endif
      END IF

      RETURN
      END SUBROUTINE ad_set_avg_tile
#endif
      END MODULE ad_set_avg_mod
