      SUBROUTINE t3dmix2 (ng, tile)
!
!svn $Id$
!***********************************************************************
!  Copyright (c) 2002-2007 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!****************************************** Alexander F. Shchepetkin ***
!                                                                      !
!  This subroutine computes horizontal harmonic mixing of tracers      !
!  along S-coordinate levels surfaces.                                 !
!                                                                      !
!***********************************************************************
!
      USE mod_param
#ifdef DIAGNOSTICS_TS
      USE mod_diags
#endif
      USE mod_grid
      USE mod_mixing
      USE mod_ocean
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
#include "tile.h"
!
#ifdef PROFILE
      CALL wclock_on (ng, iNLM, 24)
#endif
      CALL t3dmix2_tile (ng, Istr, Iend, Jstr, Jend,                    &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   nrhs(ng), nnew(ng),                            &
#ifdef MASKING
     &                   GRID(ng) % umask,                              &
     &                   GRID(ng) % vmask,                              &
#endif
     &                   GRID(ng) % Hz,                                 &
     &                   GRID(ng) % pmon_u,                             &
     &                   GRID(ng) % pnom_v,                             &
     &                   GRID(ng) % pm,                                 &
     &                   GRID(ng) % pn,                                 &
     &                   MIXING(ng) % diff2,                            &
#ifdef DIAGNOSTICS_TS
     &                   DIAGS(ng) % DiaTwrk,                           &
#endif
     &                   OCEAN(ng) % t)
#ifdef PROFILE
      CALL wclock_off (ng, iNLM, 24)
#endif
      RETURN
      END SUBROUTINE t3dmix2
!
!***********************************************************************
      SUBROUTINE t3dmix2_tile (ng, Istr, Iend, Jstr, Jend,              &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         nrhs, nnew,                              &
#ifdef MASKING
     &                         umask, vmask,                            &
#endif
     &                         Hz, pmon_u, pnom_v, pm, pn,              &
     &                         diff2,                                   &
#ifdef DIAGNOSTICS_TS
     &                         DiaTwrk,                                 &
#endif
     &                         t)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: nrhs, nnew

#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
# endif
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: pmon_u(LBi:,LBj:)
      real(r8), intent(in) :: pnom_v(LBi:,LBj:)
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: diff2(LBi:,LBj:,:)
# ifdef DIAGNOSTICS_TS
      real(r8), intent(inout) :: DiaTwrk(LBi:,LBj:,:,:,:)
# endif
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
#else
# ifdef MASKING
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: pmon_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pnom_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: diff2(LBi:UBi,LBj:UBj,NT(ng))
# ifdef DIAGNOSTICS_TS
      real(r8), intent(inout) :: DiaTwrk(LBi:UBi,LBj:UBj,N(ng),NT(ng),  &
     &                                   NDT)
# endif
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
#endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, itrc, j, k

      real(r8) :: cff, cff1

      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: FE
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: FX

#include "set_bounds.h"
!
!----------------------------------------------------------------------
!  Compute horizontal harmonic diffusion along constant S-surfaces.
!----------------------------------------------------------------------
!
      DO itrc=1,NT(ng)
        DO k=1,N(ng)
!
!  Compute XI- and ETA-components of diffusive tracer flux (T m3/s).
!
          DO j=Jstr,Jend
            DO i=Istr,Iend+1
              cff=0.25_r8*(diff2(i,j,itrc)+diff2(i-1,j,itrc))*          &
     &            pmon_u(i,j)
              FX(i,j)=cff*                                              &
     &                (Hz(i,j,k)+Hz(i-1,j,k))*                          &
     &                (t(i,j,k,nrhs,itrc)-t(i-1,j,k,nrhs,itrc))
#ifdef MASKING
              FX(i,j)=FX(i,j)*umask(i,j)
#endif
            END DO
          END DO
          DO j=Jstr,Jend+1
            DO i=Istr,Iend
              cff=0.25_r8*(diff2(i,j,itrc)+diff2(i,j-1,itrc))*          &
     &            pnom_v(i,j)
              FE(i,j)=cff*                                              &
     &                (Hz(i,j,k)+Hz(i,j-1,k))*                          &
     &                (t(i,j,k,nrhs,itrc)-t(i,j-1,k,nrhs,itrc))
#ifdef MASKING
              FE(i,j)=FE(i,j)*vmask(i,j)
#endif
            END DO
          END DO
!
! Time-step harmonic, S-surfaces diffusion term (m Tunits).
!
          DO j=Jstr,Jend
            DO i=Istr,Iend
              cff=dt(ng)*pm(i,j)*pn(i,j)*                               &
     &                   (FX(i+1,j)-FX(i,j)+                            &
     &                    FE(i,j+1)-FE(i,j))
              t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)+cff
#ifdef TS_MPDATA
              cff1=1.0_r8/Hz(i,j,k)
              t(i,j,k,3,itrc)=cff1*t(i,j,k,nnew,itrc)
#endif
#ifdef DIAGNOSTICS_TS
              DiaTwrk(i,j,k,itrc,iThdif)=cff
#endif
            END DO
          END DO
        END DO
      END DO
      RETURN
      END SUBROUTINE t3dmix2_tile
