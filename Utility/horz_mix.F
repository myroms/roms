#include "cppdefs.h"
      MODULE horz_mix_mod
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine rescales horizontal mixing coefficients according      !
!  to the grid size.  Also,  if applicable,  increases horizontal      !
!  in sponge areas.                                                    !
!                                                                      !
!  WARNING:   All biharmonic coefficients are assumed to have the      !
!             square root taken and have  m^2 s^-1/2 units.  This      !
!             will allow multiplying the  biharmonic  coefficient      !
!             to harmonic operator.                                    !
!                                                                      !
!=======================================================================
!
      implicit none
      
      CONTAINS
!
!***********************************************************************
      SUBROUTINE horz_mix (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_mixing
      USE mod_scalars
!
      integer, intent(in) :: ng, tile

#include "tile.h"

      CALL horz_mix_tile (ng, Istr, Iend, Jstr, Jend,                   &
     &                    LBi, UBi, LBj, UBj,                           &
#ifdef TS_DIF2
     &                    MIXING(ng) % diff2,                           &
#endif
#ifdef TS_DIF4
     &                    MIXING(ng) % diff4,                           &
#endif
#ifdef UV_VIS2
     &                    MIXING(ng) % visc2_p,                         &
     &                    MIXING(ng) % visc2_r,                         &
#endif
#ifdef UV_VIS4
     &                    MIXING(ng) % visc4_p,                         &
     &                    MIXING(ng) % visc4_r,                         &
#endif
     &                    GRID(ng) % grdscl,                            &
     &                    GRID(ng) % xr,                                &
     &                    GRID(ng) % yr)
      RETURN
      END SUBROUTINE horz_mix
!
!***********************************************************************
      SUBROUTINE horz_mix_tile (ng, Istr, Iend, Jstr, Jend,             &
     &                          LBi, UBi, LBj, UBj,                     &
#ifdef TS_DIF2
     &                          diff2,                                  &
#endif
#ifdef TS_DIF4
     &                          diff4,                                  &
#endif
#ifdef UV_VIS2
     &                          visc2_p,                                &
     &                          visc2_r,                                &
#endif
#ifdef UV_VIS4
     &                          visc4_p,                                &
     &                          visc4_r,                                &
#endif
     &                          grdscl, xr, yr)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
#if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      USE exchange_2d_mod
!
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj

#ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: grdscl(LBi:,LBj:)
      real(r8), intent(in) :: xr(LBi:,LBj:)
      real(r8), intent(in) :: yr(LBi:,LBj:)
# ifdef TS_DIF2
      real(r8), intent(inout) :: diff2(LBi:,LBj:,:)
# endif
# ifdef TS_DIF4
      real(r8), intent(inout) :: diff4(LBi:,LBj:,:)
# endif
# ifdef UV_VIS2
      real(r8), intent(inout) :: visc2_p(LBi:,LBj:)
      real(r8), intent(inout) :: visc2_r(LBi:,LBj:)
# endif
# ifdef UV_VIS4
      real(r8), intent(inout) :: visc4_p(LBi:,LBj:)
      real(r8), intent(inout) :: visc4_r(LBi:,LBj:)
# endif
#else
      real(r8), intent(in) :: grdscl(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: xr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: yr(LBi:UBi,LBj:UBj)
# ifdef TS_DIF2
      real(r8), intent(inout) :: diff2(LBi:UBi,LBj:UBj,NT(ng))
# endif
# ifdef TS_DIF4
      real(r8), intent(inout) :: diff4(LBi:UBi,LBj:UBj,NT(ng))
# endif
# ifdef UV_VIS2
      real(r8), intent(inout) :: visc2_p(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: visc2_r(LBi:UBi,LBj:UBj)
# endif
# ifdef UV_VIS4
      real(r8), intent(inout) :: visc4_p(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: visc4_r(LBi:UBi,LBj:UBj)
# endif
#endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: Iwrk, i, j, itrc
      real(r8) :: cff, cff1, cff2, fac

#include "set_bounds.h"

#ifdef VISC_GRID
!
!-----------------------------------------------------------------------
!  Scale horizontal viscosity according to the grid size.
!-----------------------------------------------------------------------
!
# ifdef UV_VIS2
      cff=visc2(ng)/grdmax(ng)
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          visc2_r(i,j)=cff*grdscl(i,j)
        END DO
      END DO
      cff=0.25_r8*cff
      DO j=Jstr,JendR
        DO i=Istr,IendR
          visc2_p(i,j)=cff*(grdscl(i,j  )+grdscl(i-1,j  )+              &
     &                      grdscl(i,j-1)+grdscl(i-1,j-1))
        END DO
      END DO
# endif
# ifdef UV_VIS4
      cff=visc4(ng)/(grdmax(ng)**3)
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          visc4_r(i,j)=cff*grdscl(i,j)**3
        END DO
      END DO
      cff=0.25_r8*cff
      DO j=Jstr,JendR
        DO i=Istr,IendR
          visc4_p(i,j)=cff*(grdscl(i,j  )**3+grdscl(i-1,j  )**3+        &
     &                      grdscl(i,j-1)**3+grdscl(i-1,j-1)**3)
        END DO
      END DO
# endif
#endif
#ifdef DIFF_GRID
!
!-----------------------------------------------------------------------
!  Scale horizontal diffusion according to the grid size.
!-----------------------------------------------------------------------
!
# ifdef TS_DIF2
      DO itrc=1,NT(ng)
        cff=tnu2(itrc,ng)/grdmax(ng)
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            diff2(i,j,itrc)=cff*grdscl(i,j)
          END DO
        END DO
      END DO
# endif
# ifdef TS_DIF4
      DO itrc=1,NT(ng)
        cff=tnu4(itrc,ng)/(grdmax(ng)**3)
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            diff4(i,j,itrc)=cff*grdscl(i,j)**3
          END DO
        END DO
      END DO
# endif
#endif
#ifdef SPONGE
!
!-----------------------------------------------------------------------
!  Increase horizontal mixing in the sponge areas.
!-----------------------------------------------------------------------
!
# if defined ADRIATIC1
!
!  Adriatic Sea southern sponge areas.
!
      fac=4.0_r8
#  if defined UV_VIS2
      DO i=IstrR,IendR
        DO j=JstrR,MIN(6,JendR)
          cff=visc2(ng)+REAL(6-j,r8)*(fac*visc2(ng)-visc2(ng))/6.0_r8
          visc2_r(i,j)=cff
          visc2_p(i,j)=cff
        END DO
        DO j=MAX(JstrR,7),JendR
          visc2_r(i,j)=0.0_r8
          visc2_p(i,j)=0.0_r8
        END DO
      END DO
#  endif
#  if defined TS_DIF2
      DO i=IstrR,IendR
        DO j=JstrR,MIN(6,JendR)
          cff1=tnu2(itemp,ng)+                                          &
     &         REAL(6-j,r8)*(fac*tnu2(itemp,ng)-tnu2(itemp,ng))/6.0_r8
          cff2=tnu2(isalt,ng)+                                          &
     &         REAL(6-j,r8)*(fac*tnu2(isalt,ng)-tnu2(isalt,ng))/6.0_r8
          diff2(i,j,itemp)=cff1
          diff2(i,j,isalt)=cff2
        END DO
        DO j=MAX(JstrR,7),JendR
          diff2(i,j,itemp)=0.0_r8
          diff2(i,j,isalt)=0.0_r8
        END DO
      END DO
#  endif
# elif defined CALCOFI
!
!  Southern California Bight sponge areas.
!
      fac=4.0_r8
#  if defined UV_VIS2
      DO j=JstrR,MIN(6,JendR)
        cff=visc2(ng)+REAL(6-j,r8)*(fac*visc2(ng)-visc2(ng))/6.0_r8
        DO i=IstrR,IendR
          visc2_r(i,j)=cff
          visc2_p(i,j)=cff
        END DO
      END DO
      DO j=MAX(JstrR,Mm(ng)+1-6),JendR
        cff=fac*visc2(ng)+                                              &
     &      REAL(Mm(ng)+1-j,r8)*(visc2(ng)-fac*visc2(ng))/6.0_r8
        DO i=IstrR,IendR
          visc2_r(i,j)=cff
          visc2_p(i,j)=cff
        END DO
      END DO
      DO i=IstrR,MIN(6,IendR)
        DO j=MAX(JstrR,i),MIN(Mm(ng)+1-i,JendR)
          cff=visc2(ng)+REAL(6-i,r8)*(fac*visc2(ng)-visc2(ng))/6.0_r8
          visc2_r(i,j)=cff
          visc2_p(i,j)=cff
        END DO
      END DO
#  endif
#  if defined TS_DIF2
      DO j=JstrR,MIN(6,JendR)
        cff1=tnu2(itemp,ng)+                                            &
     &       REAL(6-j,r8)*(fac*tnu2(itemp,ng)-tnu2(itemp,ng))/6.0_r8
        cff2=tnu2(isalt,ng)+                                            &
     &       REAL(6-j,r8)*(fac*tnu2(isalt,ng)-tnu2(isalt,ng))/6.0_r8
        DO i=IstrR,IendR
          diff2(i,j,itemp)=cff1
          diff2(i,j,isalt)=cff2
        END DO
      END DO
      DO j=MAX(JstrR,Mm(ng)+1-6),JendR
        cff1=fac*tnu2(itemp,ng)+                                        &
     &       REAL(Mm(ng)+1-j,r8)*(tnu2(itemp,ng)-                       &
     &                            fac*tnu2(itemp,ng))/6.0_r8
        cff2=fac*tnu2(isalt,ng)+                                        &
     &       REAL(Mm(ng)+1-j,r8)*(tnu2(isalt,ng)-                       &
     &                            fac*tnu2(isalt,ng))/6.0_r8
        DO i=IstrR,IendR
          diff2(i,j,itemp)=cff1
          diff2(i,j,isalt)=cff2
        END DO
      END DO
      DO i=IstrR,MIN(6,IendR)
        DO j=MAX(JstrR,i),MIN(Mm(ng)+1-i,JendR)
          cff1=tnu2(itemp,ng)+                                          &
     &         REAL(6-i,r8)*(fac*tnu2(itemp,ng)-tnu2(itemp,ng))/6.0_r8
          cff2=tnu2(isalt,ng)+                                          &
     &         REAL(6-i,r8)*(fac*tnu2(isalt,ng)-tnu2(isalt,ng))/6.0_r8
          diff2(i,j,itemp)=cff1
          diff2(i,j,isalt)=cff2
        END DO
      END DO
#  endif
# elif defined BISCAY
!
!  Bay of Biscay sponge areas.
!
      fac=8.0_r8
#  if defined UV_VIS2
      DO j=JstrR,MIN(6,JendR)
        cff=visc2(ng)+REAL(6-j,r8)*(fac*visc2(ng)-visc2(ng))/6.0_r8
        DO i=IstrR,IendR
          visc2_r(i,j)=cff
          visc2_p(i,j)=cff
        END DO
      END DO
      DO j=MAX(JstrR,Mm(ng)+1-6),JendR
        cff=fac*visc2(ng)+                                              &
     &      REAL(Mm(ng)+1-j,r8)*(visc2(ng)-fac*visc2(ng))/6.0_r8
        DO i=IstrR,IendR
          visc2_r(i,j)=cff
          visc2_p(i,j)=cff
        END DO
      END DO
      DO i=IstrR,MIN(6,IendR)
        DO j=MAX(JstrR,i),MIN(Mm(ng)+1-i,JendR)
          cff=visc2(ng)+REAL(6-i,r8)*(fac*visc2(ng)-visc2(ng))/6.0_r8
          visc2_r(i,j)=cff
          visc2_p(i,j)=cff
        END DO
      END DO
#  endif
#  if defined TS_DIF2
      DO j=JstrR,MIN(6,JendR)
        cff1=tnu2(itemp,ng)+                                            &
     &       REAL(6-j,r8)*(fac*tnu2(itemp,ng)-tnu2(itemp,ng))/6.0_r8
        cff2=tnu2(isalt,ng)+                                            &
     &       REAL(6-j,r8)*(fac*tnu2(isalt,ng)-tnu2(isalt,ng))/6.0_r8
        DO i=IstrR,IendR
          diff2(i,j,itemp)=cff1
          diff2(i,j,isalt)=cff2
        END DO
      END DO
      DO j=MAX(JstrR,Mm(ng)+1-6),JendR
        cff1=8.0_r8*tnu2(itemp,ng)+                                     &
     &       REAL(Mm(ng)+1-j,r8)*(tnu2(itemp,ng)-                       &
     &                            fac*tnu2(itemp,ng))/6.0_r8
        cff2=8.0_r8*tnu2(isalt,ng)+                                     &
     &       REAL(Mm(ng)+1-j,r8)*(tnu2(isalt,ng)-                       &
     &                            fac*tnu2(isalt,ng))/6.0_r8
        DO i=IstrR,IendR
          diff2(i,j,itemp)=cff1
          diff2(i,j,isalt)=cff2
        END DO
      END DO
      DO i=IstrR,MIN(6,IendR)
        DO j=MAX(JstrR,i),MIN(Mm(ng)+1-i,JendR)
          cff1=tnu2(itemp,ng)+                                          &
     &         REAL(6-j,r8)*(fac*tnu2(itemp,ng)-tnu2(itemp,ng))/6.0_r8
          cff2=tnu2(isalt,ng)+                                          &
     &         REAL(6-j,r8)*(fac*tnu2(isalt,ng)-tnu2(isalt,ng))/6.0_r8
          diff2(i,j,itemp)=cff1
          diff2(i,j,isalt)=cff2
        END DO
      END DO
#  endif
# elif defined USWEST
!
!  US West Coast sponge areas.
!
      cff1=150000.0_r8                    ! width (m) of sponge layer
      Iwrk=INT(cff1*REAL(Lm(ng),r8)/                                    &
     &         (xr(Lm(ng)+1,INT(0.5_r8*REAL(Mm(ng),r8)))-               &
     &          xr(1,INT(0.5_r8*REAL(Mm(ng),r8)))))
#  if defined UV_VIS2
      DO j=JstrR,JendR
        DO i=IstrR,MIN(Iwrk,IendR)
          cff=visc2(ng)*0.5_r8*(1.0_r8+COS(pi*REAL(i,r8)/REAL(Iwrk,r8)))
          visc2_r(i,j)=cff
          visc2_p(i,j)=cff
        END DO
        DO i=MAX(IstrR,Lm(ng)+1-Iwrk),IendR
          cff=visc2(ng)*0.5_r8*(1.0_r8+COS(pi*REAL(Lm(ng)+1-i,r8)/          &
     &                                     REAL(Iwrk,r8)))
          visc2_r(i,j)=cff
          visc2_p(i,j)=cff
        END DO
      END DO
      DO j=MAX(JstrR,Mm(ng)+1-Iwrk),MIN(Mm(ng)+1,JendR)
        DO i=MAX(IstrR,Mm(ng)+1-j),MIN(Lm(ng)+1-Mm(ng)+1+j,IendR)
          cff=visc2(ng)*0.5_r8*(1.0_r8+COS(pi*REAL(Mm(ng)+1-j,r8)/      &
     &                                     REAL(Iwrk,r8)))
          visc2_r(i,j)=cff
          visc2_p(i,j)=cff
        END DO
      END DO
!!    DO j=JstrR,MIN(Iwrk,JendR)
!!      DO i=MAX(IstrR,j),MIN(Lm(ng)+1-j,IendR)
!!        cff=visc2(ng)*0.5_r8*(1.0_r8+COS(pi*REAL(j,r8)/REAL(Iwrk,r8)))
!!        visc2_r(i,j)=cff
!!        visc2_p(i,j)=cff
!!      END DO
!!    END DO
#  endif
#  if defined TS_DIF2
      DO itrc=1,NT(ng)
        cff=tnu2(itrc,ng)*0.5_r8
        DO j=JstrR,JendR
          DO i=IstrR,MIN(Iwrk,IendR)
            diff2(i,j,itrc)=cff*                                        &
     &                      (1.0_r8+COS(pi*REAL(i,r8)/REAL(Iwrk,r8)))
          END DO
          DO i=MAX(IstrR,Lm(ng)+1-Iwrk),IendR
            diff2(i,j,itrc)=cff*                                        &
     &                      (1.0_r8+COS(pi*REAL(Lm(ng)+1-i,r8)/         &
     &                                  REAL(Iwrk,r8)))
          END DO
        END DO
        DO j=MAX(JstrR,Mm(ng)+1-Iwrk),MIN(Mm(ng)+1,JendR)
          DO i=MAX(IstrR,Mm(ng)+1-j),MIN(Lm(ng)+1-Mm(ng)+1+j,IendR)
            diff2(i,j,itrc)=cff*                                        &
     &                      (1.0_r8+COS(pi*REAL(Mm(ng)+1-j,r8)/         &
     &                                  REAL(Iwrk,r8)))
          END DO
        END DO
!!      DO j=JstrR,MIN(Iwrk,JendR)
!!        DO i=MAX(IstrR,j),MIN(Lm(ng)+1-j,IendR)
!!          diff2(i,j,itrc)=cff*
!!   &                      (1.0_r8+COS(pi*REAL(j,r8)/REAL(Iwrk,r8)))
!!        END DO
!!      END DO
      END DO
#  endif
# endif
#endif
#if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Exchange boundary data.
!-----------------------------------------------------------------------
!
# ifdef UV_VIS2
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, visc2_r)
      CALL exchange_p2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, visc2_p)
# endif
# ifdef UV_VIS4
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, visc4_r)
      CALL exchange_p2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, visc4_p)
# endif
# ifdef TS_DIF2
      DO itrc=1,NT(ng)
        CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,             &
     &                          LBi, UBi, LBj, UBj, diff2(:,:,itrc))
      END DO
# endif
# ifdef TS_DIF4
      DO itrc=1,NT(ng)
        CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,             &
     &                          LBi, UBi, LBj, UBj, diff4(:,:,itrc))
      END DO
# endif
#endif
      RETURN
      END SUBROUTINE horz_mix_tile
      END MODULE horz_mix_mod
