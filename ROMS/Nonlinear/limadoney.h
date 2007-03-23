      SUBROUTINE biology (ng,tile)
!
!svn $Id$
!************************************************** Hernan G. Arango ***
!  Copyright (c) 2002-2007 The ROMS/TOMS Group          Katja Fennel   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!***********************************************************************
!                                                                      !
!  This routine computes the biological sources and sinks and adds     !
!  them to the global biological fields.                               !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!      Lima I.D., S.C. Doney (2004) A three-dimensional, multinutrient !
!      and size-structured ecosystem model for the North Atlantic,     !
!      Global Biogeochemical Cycles 18, GB3019                         !
!      doi:10.1029/2003GB002146.                                       !
!                                                                      !
!  Adapted from Ivan Limas code by Katja Fennel.                       !
!                                                                      !
!***********************************************************************
!
      USE mod_param
#ifdef DIAGNOSTICS_BIO
      USE mod_diags
#endif
      USE mod_forces
      USE mod_grid
      USE mod_ocean
      USE mod_stepping
!
      integer, intent(in) :: ng, tile

#include "tile.h"
!
#ifdef PROFILE
      CALL wclock_on (ng, iNLM, 15)
#endif
      CALL biology_tile (ng, Istr, Iend, Jstr, Jend,                    &
     &                   LBi, UBi, LBj, UBj, N(ng), NT(ng),             &
     &                   nnew(ng),                                      &
#ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#endif
     &                   GRID(ng) % Hz,                                 &
     &                   GRID(ng) % z_r,                                &
     &                   GRID(ng) % z_w,                                &
     &                   FORCES(ng) % srflx,                            &
#ifdef DIAGNOSTICS_BIO
     &                   DIAGS(ng) % DiaBio2d,                          &
     &                   DIAGS(ng) % DiaBio3d,                          &
#endif
     &                   OCEAN(ng) % t)

#ifdef PROFILE
      CALL wclock_off (ng, iNLM, 15)
#endif
      RETURN
      END SUBROUTINE biology
!
!-----------------------------------------------------------------------
      SUBROUTINE biology_tile (ng, Istr, Iend, Jstr, Jend,              &
     &                         LBi, UBi, LBj, UBj, UBk, UBt,            &
     &                         nnew,                                    &
#ifdef MASKING
     &                         rmask,                                   &
#endif
     &                         Hz, z_r, z_w, srflx,                     &
#ifdef DIAGNOSTICS_BIO
     &                         DiaBio2d, DiaBio3d,                      &
#endif
     &                         t)
!-----------------------------------------------------------------------
!
      USE mod_param
      USE mod_biology
      USE mod_ncparam
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj, UBk, UBt
      integer, intent(in) :: nnew

#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
# endif
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: srflx(LBi:,LBj:)
# ifdef DIAGNOSTICS_BIO
      real(r8), intent(inout) :: DiaBio2d(LBi:,LBj:,:)
      real(r8), intent(inout) :: DiaBio3d(LBi:,LBj:,:,:)
# endif
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
#else
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:UBk)
      real(r8), intent(in) :: srflx(LBi:UBi,LBj:UBj)
# ifdef DIAGNOSTICS_BIO
      real(r8), intent(inout) :: DiaBio2d(LBi:UBi,LBj:UBj,NDbio2d)
      real(r8), intent(inout) :: DiaBio3d(LBi:UBi,LBj:UBj,UBk,NDbio3d)
# endif
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,UBk,3,UBt)
#endif
!
!  Local variable declarations.
!
      integer, parameter :: Nsink = 2

      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: Iter, i, indx, isink, ibio, ivar, j, k, ks

      integer, dimension(Nsink) :: idsink

      real(r8), parameter :: eps = 1.0e-20_r8

      real(r8), parameter :: Kelvin = 273.15_r8

      real(r8) :: Att, dtdays, PAR

      real(r8) :: FTEMP, QNO3_1, QNO3_2, QNH4_1, QNH4_2
      real(r8) :: VNCMAX1, VNCMAX2, F_NIT1, F_NIT2, F_NUT2
      real(r8) :: Q1, Q2, QZOO, QDET1, QDET2, UQ1, UQ2
      real(r8) :: THETA_N1, THETA_N2, THETA_C1, THETA_C2
      real(r8) :: FRATIO, FALLOC, ZD_ALLOC
      real(r8) :: growth1_NO3, growth1_NH4, growth2_NO3, growth2_NH4 
      real(r8) :: PHOTO1, PHOTO2, accl1, accl2, graze1, graze2
      real(r8) :: cff, cff1, cff2, cff3, cff4, cff5
      real(r8) :: fac1, fac2, fac3
      real(r8) :: cffL, cffR, cu, dltL, dltR

      real(r8), dimension(Nsink) :: Wbio

      integer, dimension(PRIVATE_1D_SCRATCH_ARRAY,N(ng)) :: ksource

      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY) :: PARsur

      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,N(ng),NT(ng)) :: Bio

      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N(ng)) :: FC

      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,N(ng)) :: Hz_inv
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,N(ng)) :: Hz_inv2
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,N(ng)) :: Hz_inv3
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,N(ng)) :: WL
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,N(ng)) :: WR
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,N(ng)) :: bL
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,N(ng)) :: bR
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,N(ng)) :: qc

#include "set_bounds.h"

#ifdef DIAGNOSTICS_BIO
!
!-----------------------------------------------------------------------
! If appropriate, initialize time-averaged diagnostic arrays.
!-----------------------------------------------------------------------
!
      IF (((iic(ng).gt.ntsDIA(ng)).and.                                 &
     &     (MOD(iic(ng),nDIA(ng)).eq.1)).or.                            &
     &    ((nrrec.gt.0).and.(iic(ng).eq.ntstart))) THEN
        DO ivar=1,NDbio2d
          DO j=Jstr,Jend
            DO i=Istr,Iend
              DiaBio2d(i,j,ivar)=0.0_r8
            END DO
          END DO
        END DO
        DO ivar=1,NDbio3d
          DO k=1,N(ng)
            DO j=Jstr,Jend
              DO i=Istr,Iend
                DiaBio3d(i,j,k,ivar)=0.0_r8
              END DO
            END DO
          END DO
        END DO
      END IF
#endif
!
!-----------------------------------------------------------------------
!  Add biological Source/Sink terms.
!-----------------------------------------------------------------------
!
!  Set time-stepping according to the number of iterations.
!
      dtdays=dt(ng)*sec2day/REAL(BioIter(ng),r8)
!
!  Set vertical sinking indentification vector.
!
      idsink(1)=iDe2N
      idsink(2)=iDe2C
!
!  Set vertical sinking velocity vector in the same order as the
!  identification vector, IDSINK.
!
      Wbio(1)=wDet2(ng)               ! large/sinking detritus
      Wbio(2)=wDet2(ng)
!
!  Compute inverse thickness to avoid repeated divisions.
!
      J_LOOP : DO j=Jstr,Jend
        DO k=1,N(ng)
          DO i=Istr,Iend
            Hz_inv(i,k)=1.0_r8/Hz(i,j,k)
          END DO
        END DO
        DO k=1,N(ng)-1
          DO i=Istr,Iend
            Hz_inv2(i,k)=1.0_r8/(Hz(i,j,k)+Hz(i,j,k+1))
          END DO
        END DO
        DO k=2,N(ng)-1
          DO i=Istr,Iend
            Hz_inv3(i,k)=1.0_r8/(Hz(i,j,k-1)+Hz(i,j,k)+Hz(i,j,k+1))
          END DO
        END DO
!
!  Extract biological variables from tracer arrays, place them into
!  scratch arrays, and restrict their values to be non-negative.
!
        DO ibio=1,NBT
          indx=idbio(ibio)
          DO k=1,N(ng)
            DO i=Istr,Iend
              Bio(i,k,indx)=MAX(t(i,j,k,nnew,indx)*Hz_inv(i,k),0.0_r8)
            END DO
          END DO
        END DO
!
!  Extract potential temperature.
!
        DO k=1,N(ng)
          DO i=Istr,Iend
            Bio(i,k,itemp)=t(i,j,k,nnew,itemp)*Hz_inv(i,k)
          END DO
        END DO
!
!  Calculate surface Photosynthetically Available Radiation (PAR).  The
!  net shortwave radiation is scaled back to Watts/m2 and multiplied by
!  the fraction that is photosynthetically available, PARfrac.
!
        DO i=Istr,Iend 
          PARsur(i)=MAX(PARfrac(ng)*srflx(i,j)*rho0*Cp,0.0_r8)
        END DO
!
!=======================================================================
!  Start internal iterations to achieve convergence of the nonlinear
!  backward-implicit solution.
!=======================================================================
!
!  During the iterative procedure a series of fractional time steps are
!  performed in a chained mode (splitting by different biological
!  conversion processes) in sequence of the main food chain.  In all
!  stages the concentration of the component being consumed is treated
!  in fully implicit manner, so the algorithm guarantees non-negative
!  values, no matter how strong s the concentration of active consuming
!  component (Phytoplankton or Zooplankton).  The overall algorithm,
!  as well as any stage of it, is formulated in conservative form
!  (except explicit sinking) in sense that the sum of concentration of
!  all components is conserved.
!
        ITER_LOOP: DO Iter=1,BioIter(ng)
!
          DO i=Istr,Iend 
            PAR=PARsur(i)
            DO k=N(ng),1,-1
!
!  Attenuate the light to the center of the grid cell.
!
              Att=EXP(-0.5_r8*(AttSW(ng)+                               &
     &                         AttChl(ng)*(Bio(i,k,iChl1)+              &
     &                                     Bio(i,k,iChl2)))*            &
     &                         (z_w(i,j,k)-z_w(i,j,k-1)))
              PAR=PAR*Att
!
              FTEMP=MIN(1.0_r8,                                         &
     &                  EXP(-e_a(ng)*(1.0_r8/(Bio(i,k,itemp)+Kelvin)-   &
     &                                1.0_r8/(Tref(ng)+Kelvin))))
!
!  Intracellular ratios
!
              THETA_N1=Bio(i,k,iChl1)/(Bio(i,k,iPh1N)+eps)
              THETA_C1=Bio(i,k,iChl1)/(Bio(i,k,iPh1C)+eps)
              Q1=Bio(i,k,iPh1N)/(Bio(i,k,iPh1C)+eps)
!
              THETA_N2=Bio(i,k,iChl2)/(Bio(i,k,iPh2N)+eps)
              THETA_C2=Bio(i,k,iChl2)/(Bio(i,k,iPh2C)+eps)
              Q2=Bio(i,k,iPh2N)/(Bio(i,k,iPh2C)+eps)
!
              QZOO =Bio(i,k,iZooN)/(Bio(i,k,iZooC)+eps)
              QDET1=Bio(i,k,iDe1N)/(Bio(i,k,iDe1C)+eps)
              QDET2=Bio(i,k,iDe2N)/(Bio(i,k,iDe2C)+eps)
!
              UQ1=MIN(Q1 ,q_max(ng))
              UQ1=MAX(UQ1,q_min(ng))
              UQ2=MIN(Q2 ,q_max(ng))
              UQ2=MAX(UQ2,q_min(ng))
!
              F_NIT1 = (UQ1-q_min(ng))/(q_max(ng)-q_min(ng))
              VNCMAX1 = vncref(ng)*(1.0_r8-F_NIT1)/                     &
     &                  (1.0_r8-F_NIT1+0.015_r8)*FTEMP
!
              F_NIT2 = (UQ2-q_min(ng))/(q_max(ng)-q_min(ng))
              VNCMAX2 = vncref(ng)*(1.0_r8-F_NIT2)/                     &
     &                  (1.0_r8-F_NIT2+0.015_r8)*FTEMP
!
!  Determine diatom limitation factor for carbon fixation.
!
              F_NUT2 = F_NIT2 ! no Si yet
!
!  Compute nitrate and ammonia uptake/inhibition rates
!
              cff1 = Bio(i,k,iNO3_)/p1_kno3(ng)
              cff2 = Bio(i,k,iNH4_)/p1_knh4(ng)
              QNO3_1 = cff1/(1.0_r8+cff1+cff2)
              QNH4_1 = cff2/(1.0_r8+cff1+cff2)
!
              cff1 = Bio(i,k,iNO3_)/p2_kno3(ng)
              cff2 = Bio(i,k,iNH4_)/p2_knh4(ng)
              QNO3_2 = cff1/(1.0_r8+cff1+cff2)
              QNH4_2 = cff2/(1.0_r8+cff1+cff2)
!
!  Phytoplankton growth rates.
!
              growth1_NO3 = Bio(i,k,iPh1N)*( (VNCMAX1*QNO3_1)/          &
     &                           (Q1+eps)-r_ref(ng)*FTEMP )
              growth1_NH4 = Bio(i,k,iPh1N)*( (VNCMAX1*QNH4_1)/          &
     &                           (Q1+eps)-r_ref(ng)*FTEMP )
              growth2_NO3 = Bio(i,k,iPh2N)*( (VNCMAX2*QNO3_2)/          &
     &                           (Q2+eps)-r_ref(ng)*FTEMP )
              growth2_NH4 = Bio(i,k,iPh2N)*( (VNCMAX2*QNH4_2)/          &
     &                           (Q2+eps)-r_ref(ng)*FTEMP )
!
              cff1 = growth1_NO3+growth2_NO3
              cff2 = growth1_NH4+growth2_NH4
              FRATIO = MIN(cff1/(cff1+cff2+eps),1.0_r8)
!
              Bio(i,k,iPh1N)=Bio(i,k,iPh1N)+dtdays*                     &
     &                            (growth1_NO3+growth1_NH4)
              Bio(i,k,iPh2N)=Bio(i,k,iPh2N)+dtdays*                     &
     &                            (growth2_NO3+growth2_NH4)
              Bio(i,k,iNO3_)=Bio(i,k,iNO3_)-dtdays*                     &
     &                            (growth1_NO3+growth2_NO3)
              Bio(i,k,iNH4_)=Bio(i,k,iNH4_)-dtdays*                     &
     &                            (growth1_NH4+growth2_NH4)
!
!  Compute photosynthetic rate (C fixation rate) and photoacclimation.
!
              fac1 = (QNO3_1 + QNH4_1)*VNCMAX1
              fac2 = (QNO3_2 + QNH4_2)*VNCMAX2
!
              cff = p_cref(ng)*F_NIT1*FTEMP
              cff1 = cff*( 1.0_r8-EXP((-alphaPI(ng)*THETA_C1*PAR)/      &
     &              (cff+eps)) )
              PHOTO1 =(cff1-biolambda(ng)*fac1-r_ref(ng))*Bio(i,k,iPh1C)
!
              cff = p_cref(ng)*F_NUT2*FTEMP
              cff2 = cff*( 1.0_r8-EXP((-alphaPI(ng)*THETA_C2*PAR)/      & 
     &              (cff+eps)) )
              PHOTO2 =(cff2-biolambda(ng)*fac2-r_ref(ng))*Bio(i,k,iPh2C)
!
              cff = thetaN0(ng)*cff1/(alphaPI(ng)*PAR*THETA_C1+eps)
              accl1 = ((cff*fac1)/(THETA_C1+eps)-r_ref(ng))             &
     &               *Bio(i,k,iChl1)
!
              cff = thetaN0(ng)*cff2/(alphaPI(ng)*PAR*THETA_C2+eps)
              accl2 = ((cff*fac2)/(THETA_C2+eps)-r_ref(ng))             &
     &               *Bio(i,k,iChl2)
!
              Bio(i,k,iPh1C)=Bio(i,k,iPh1C)+dtdays*PHOTO1
              Bio(i,k,iPh2C)=Bio(i,k,iPh2C)+dtdays*PHOTO2
              Bio(i,k,iChl1)=Bio(i,k,iChl1)+dtdays*accl1
              Bio(i,k,iChl2)=Bio(i,k,iChl2)+dtdays*accl2
!
!  Nitrification only below 1% of surface light or in the dark
! 
              if (PAR.le.0.01_r8*PARsur(i)) then
                 cff = p_nitr(ng)*Bio(i,k,iNH4_)
                 Bio(i,k,iNO3_)=Bio(i,k,iNO3_)+dtdays*cff
                 Bio(i,k,iNH4_)=Bio(i,k,iNH4_)-dtdays*cff
              endif
!              
!  Compute zooplankton grazing and losses.
!
              graze1 = Bio(i,k,iZooN)*z_umax1(ng)*FTEMP*Bio(i,k,iPh1N)* &
     &                 Bio(i,k,iPh1N)/(Bio(i,k,iPh1N)*Bio(i,k,iPh1N)+   &
     &                 z_grz(ng)*z_grz(ng))
!
              graze2 = Bio(i,k,iZooN)*z_umax2(ng)*FTEMP*Bio(i,k,iPh2N)* &
     &                 Bio(i,k,iPh2N)/(Bio(i,k,iPh2N)*Bio(i,k,iPh2N)+   &
     &                 z_grz(ng)*z_grz(ng))
!
              FALLOC = falloc0(ng)+(1.0_r8-falloc0(ng))*(1.0_r8-FRATIO)
              ZD_ALLOC = (0.7_r8*graze2+0.3_r8*graze1)/                 &
     &                   (graze1+graze2+eps)
!
              cff1 = p1_excr(ng)*Bio(i,k,iPh1N)
              cff2 = p2_excr(ng)*Bio(i,k,iPh2N)
              cff3 = p1_mort2(ng)*Bio(i,k,iPh1N)*Bio(i,k,iPh1N)
              cff4 = p2_mort2(ng)*Bio(i,k,iPh2N)*Bio(i,k,iPh2N)
!
              Bio(i,k,iPh1N)=Bio(i,k,iPh1N)-dtdays*(graze1+cff1+cff3) 
              Bio(i,k,iPh2N)=Bio(i,k,iPh2N)-dtdays*(graze2+cff2+cff4)
              Bio(i,k,iPh1C)=Bio(i,k,iPh1C)-dtdays*(graze1+cff1+cff3)/  &
     &                       (Q1+eps) 
              Bio(i,k,iPh2C)=Bio(i,k,iPh2C)-dtdays*(graze2+cff2+cff4)/  &
     &                       (Q2+eps)
              Bio(i,k,iChl1)=Bio(i,k,iChl1)-dtdays*(graze1+cff1+cff3)*  &
     &                       THETA_N1 
              Bio(i,k,iChl2)=Bio(i,k,iChl2)-dtdays*(graze2+cff2+cff4)*  &
     &                       THETA_N2
!
              Bio(i,k,iDe1N)=Bio(i,k,iDe1N)+dtdays*(cff1+cff3)
              Bio(i,k,iDe2N)=Bio(i,k,iDe2N)+dtdays*(cff2+cff4)
              Bio(i,k,iDe1C)=Bio(i,k,iDe1C)+dtdays*(cff1+cff3)/(Q1+eps)
              Bio(i,k,iDe2C)=Bio(i,k,iDe2C)+dtdays*(cff2+cff4)/(Q2+eps)
!
!   Zooplankton losses:
!
              cff = z_mort(ng)*FTEMP*Bio(i,k,iZooN)+                    &
     &                  z_mort2(ng)*Bio(i,k,iZooN)*Bio(i,k,iZooN)
!
              Bio(i,k,iZooN)=Bio(i,k,iZooN)+dtdays*                     &
     &                  ( (1.0_r8-zegest0(ng))*(graze1+graze2) - cff)
              Bio(i,k,iZooC)=Bio(i,k,iZooC)+dtdays*                     &
     &           ( (1.0_r8-zegest0(ng))*(graze1/(Q1+eps)                &
     &           +graze2/(Q2+eps)) - cff/(QZOO+eps) )
! 
              Bio(i,k,iNH4_)=Bio(i,k,iNH4_)+dtdays*(                    &
     &             FALLOC*zegest0(ng)*(graze1+graze2)                   &
     &            +d1_rem(ng)*FTEMP*Bio(i,k,iDe1N)                      &
     &            +d2_rem(ng)*FTEMP*Bio(i,k,iDe2N)  )
!
              Bio(i,k,iDe1N)=Bio(i,k,iDe1N)+dtdays*(                    &
     &             (1.0_r8-ZD_ALLOC)*cff                                &
     &            +(1.0_r8-FALLOC)*zegest0(ng)*(graze1+graze2)          &
     &            -d1_rem(ng)*FTEMP*Bio(i,k,iDe1N) )
              Bio(i,k,iDe2N)=Bio(i,k,iDe2N)+dtdays*(ZD_ALLOC*cff        &
     &            -d2_rem(ng)*FTEMP*Bio(i,k,iDe2N) )         
              Bio(i,k,iDe1C)=Bio(i,k,iDe1C)+dtdays*(                    &
     &             (1.0_r8-ZD_ALLOC)*cff/(QZOO+eps)                     &
     &            +(1.0_r8-FALLOC)*zegest0(ng)*(graze1/(Q1+eps)         &
     &            +graze2/(Q2+eps))                                     &
     &            -d1_rem(ng)*FTEMP*Bio(i,k,iDe1C))
              Bio(i,k,iDe2C)=Bio(i,k,iDe2C)+dtdays*(                    &
     &             ZD_ALLOC*cff/(QZOO+eps)                              &
     &            -d2_rem(ng)*FTEMP*Bio(i,k,iDe2C))
!!            Bio(i,k,iDe1C)=Bio(i,k,iDe1C)+dtdays*(                    &
!!   &             (1.0_r8-ZD_ALLOC)*cff/(QZOO+eps)                     &
!!   &            +(1.0_r8-FALLOC)*zegest0(ng)*(graze1/(Q1+eps)         &
!!   &            +graze2/(Q2+eps))                                     &
!!   &            -d1_rem(ng)*FTEMP*Bio(i,k,iDe1N)/(QDET1+eps))
!!            Bio(i,k,iDe2C)=Bio(i,k,iDe2C)+dtdays*(                    &
!!   &             ZD_ALLOC*cff/(QZOO+eps)                              &
!!   &            -d2_rem(ng)*FTEMP*Bio(i,k,iDe2N)/(QDET2+eps))
!
!  Attenuate the light to the bottom of the grid cell.
!
              PAR=PAR*Att
!
            END DO
          END DO
!
!-----------------------------------------------------------------------
!  Vertical sinking terms.
!-----------------------------------------------------------------------
!
!  Reconstruct vertical profile of selected biological constituents
!  "Bio(:,:,isink)" in terms of a set of parabolic segments within each
!  grid box. Then, compute semi-Lagrangian flux due to sinking.
!
          SINK_LOOP: DO isink=1,Nsink
            indx=idsink(isink)
!
!  Copy concentration of biological particulates into scratch array
!  "qc" (q-central, restrict it to be positive) which is hereafter
!  interpreted as a set of grid-box averaged values for biogeochemical
!  constituent concentration.
!
            DO k=1,N(ng)
              DO i=Istr,Iend
                qc(i,k)=Bio(i,k,indx)
              END DO
            END DO
!
            DO k=N(ng)-1,1,-1
              DO i=Istr,Iend
                FC(i,k)=(qc(i,k+1)-qc(i,k))*Hz_inv2(i,k)
              END DO
            END DO
            DO k=2,N(ng)-1
              DO i=Istr,Iend
                dltR=Hz(i,j,k)*FC(i,k)
                dltL=Hz(i,j,k)*FC(i,k-1)
                cff=Hz(i,j,k-1)+2.0_r8*Hz(i,j,k)+Hz(i,j,k+1)
                cffR=cff*FC(i,k)
                cffL=cff*FC(i,k-1)
!
!  Apply PPM monotonicity constraint to prevent oscillations within the
!  grid box.
!
                IF ((dltR*dltL).le.0.0_r8) THEN
                  dltR=0.0_r8
                  dltL=0.0_r8
                ELSE IF (ABS(dltR).gt.ABS(cffL)) THEN
                  dltR=cffL
                ELSE IF (ABS(dltL).gt.ABS(cffR)) THEN
                  dltL=cffR
                END IF
!
!  Compute right and left side values (bR,bL) of parabolic segments
!  within grid box Hz(k); (WR,WL) are measures of quadratic variations. 
!
!  NOTE: Although each parabolic segment is monotonic within its grid
!        box, monotonicity of the whole profile is not guaranteed,
!        because bL(k+1)-bR(k) may still have different sign than
!        qc(i,k+1)-qc(i,k).  This possibility is excluded,
!        after bL and bR are reconciled using WENO procedure.
!
                cff=(dltR-dltL)*Hz_inv3(i,k)
                dltR=dltR-cff*Hz(i,j,k+1)
                dltL=dltL+cff*Hz(i,j,k-1)
                bR(i,k)=qc(i,k)+dltR
                bL(i,k)=qc(i,k)-dltL
                WR(i,k)=(2.0_r8*dltR-dltL)**2
                WL(i,k)=(dltR-2.0_r8*dltL)**2
              END DO
            END DO
            cff=1.0E-14_r8
            DO k=2,N(ng)-2
              DO i=Istr,Iend
                dltL=MAX(cff,WL(i,k  ))
                dltR=MAX(cff,WR(i,k+1))
                bR(i,k)=(dltR*bR(i,k)+dltL*bL(i,k+1))/(dltR+dltL)
                bL(i,k+1)=bR(i,k)
              END DO
            END DO
            DO i=Istr,Iend
              FC(i,N(ng))=0.0_r8            ! NO-flux boundary condition
#if defined LINEAR_CONTINUATION
              bL(i,N(ng))=bR(i,N(ng)-1)
              bR(i,N(ng))=2.0_r8*qc(i,N(ng))-bL(i,N(ng))
#elif defined NEUMANN
              bL(i,N(ng))=bR(i,N(ng)-1)
              bR(i,N(ng))=1.5*qc(i,N(ng))-0.5_r8*bL(i,N(ng))
#else
              bR(i,N(ng))=qc(i,N(ng))       ! default strictly monotonic
              bL(i,N(ng))=qc(i,N(ng))       ! conditions
              bR(i,N(ng)-1)=qc(i,N(ng))
#endif
#if defined LINEAR_CONTINUATION 
              bR(i,1)=bL(i,2)
              bL(i,1)=2.0_r8*qc(i,1)-bR(i,1)
#elif defined NEUMANN
              bR(i,1)=bL(i,2)
              bL(i,1)=1.5_r8*qc(i,1)-0.5_r8*bR(i,1)
#else  
              bL(i,2)=qc(i,1)               ! bottom grid boxes are
              bR(i,1)=qc(i,1)               ! re-assumed to be
              bL(i,1)=qc(i,1)               ! piecewise constant.
#endif
            END DO
!
!  Apply monotonicity constraint again, since the reconciled interfacial
!  values may cause a non-monotonic behavior of the parabolic segments
!  inside the grid box.
!
            DO k=1,N(ng)
              DO i=Istr,Iend
                dltR=bR(i,k)-qc(i,k)
                dltL=qc(i,k)-bL(i,k)
                cffR=2.0_r8*dltR
                cffL=2.0_r8*dltL
                IF ((dltR*dltL).lt.0.0_r8) THEN
                  dltR=0.0_r8
                  dltL=0.0_r8
                ELSE IF (ABS(dltR).gt.ABS(cffL)) THEN
                  dltR=cffL
                ELSE IF (ABS(dltL).gt.ABS(cffR)) THEN
                  dltL=cffR
                END IF
                bR(i,k)=qc(i,k)+dltR
                bL(i,k)=qc(i,k)-dltL
              END DO
            END DO
!
!  After this moment reconstruction is considered complete. The next
!  stage is to compute vertical advective fluxes, FC. It is expected
!  that sinking may occurs relatively fast, the algorithm is designed
!  to be free of CFL criterion, which is achieved by allowing
!  integration bounds for semi-Lagrangian advective flux to use as
!  many grid boxes in upstream direction as necessary.
!
!  In the two code segments below, WL is the z-coordinate of the
!  departure point for grid box interface z_w with the same indices;
!  FC is the finite volume flux; ksource(:,k) is index of vertical
!  grid box which contains the departure point (restricted by N(ng)). 
!  During the search: also add in content of whole grid boxes
!  participating in FC.
!
            cff=dtdays*ABS(Wbio(isink))
            DO k=1,N(ng)
              DO i=Istr,Iend
                FC(i,k-1)=0.0_r8
                WL(i,k)=z_w(i,j,k-1)+cff
                WR(i,k)=Hz(i,j,k)*qc(i,k)
                ksource(i,k)=k
              END DO
            END DO
            DO k=1,N(ng)
              DO ks=k,N(ng)-1
                DO i=Istr,Iend
                  IF (WL(i,k).gt.z_w(i,j,ks)) THEN
                    ksource(i,k)=ks+1
                    FC(i,k-1)=FC(i,k-1)+WR(i,ks)
                  END IF
                END DO
              END DO
            END DO
!
!  Finalize computation of flux: add fractional part.
!
            DO k=1,N(ng)
              DO i=Istr,Iend
                ks=ksource(i,k)
                cu=MIN(1.0_r8,(WL(i,k)-z_w(i,j,ks-1))*Hz_inv(i,ks))
                FC(i,k-1)=FC(i,k-1)+                                    &
     &                    Hz(i,j,ks)*cu*                                &
     &                    (bL(i,ks)+                                    &
     &                     cu*(0.5_r8*(bR(i,ks)-bL(i,ks))-              &
     &                         (1.5_r8-cu)*                             &
     &                         (bR(i,ks)+bL(i,ks)-                      &
     &                          2.0_r8*qc(i,ks))))
              END DO
            END DO
            DO k=1,N(ng)
              DO i=Istr,Iend
                Bio(i,k,indx)=qc(i,k)+(FC(i,k)-FC(i,k-1))*Hz_inv(i,k)
              END DO
            END DO

#ifdef BIO_SEDIMENT
!
!  Particulate flux reaching the seafloor is remineralized and returned
!  to the dissolved nitrate pool. Without this conversion, particulate
!  material falls out of the system. This is a temporary fix to restore
!  total nitrogen conservation. It will be replaced later by a
!  parameterization that includes the time delay of remineralization
!  and dissolved oxygen.
!  
            cff2=4.0_r8/16.0_r8
            IF (indx.eq.iDe2N) THEN
              DO i=Istr,Iend
                cff1=FC(i,0)*Hz_inv(i,1)
# ifdef DENITRIFICATION
                Bio(i,1,iNH4_)=Bio(i,1,iNH4_)+cff1*cff2
# else
                Bio(i,1,iNH4_)=Bio(i,1,iNH4_)+cff1
# endif
# ifdef DIAGNOSTICS_BIO
                DiaBio2d(i,j,iDNIT)=DiaBio2d(i,j,iDNIT)+                &
     &                              (1.0_r8-cff2)*cff1
# endif
              END DO
            END IF
#endif
          END DO SINK_LOOP
        END DO ITER_LOOP
!
!-----------------------------------------------------------------------
!  Update global tracer variables.
!-----------------------------------------------------------------------
!
        DO ibio=1,NBT
          indx=idbio(ibio)
          DO k=1,N(ng)
            DO i=Istr,Iend
              t(i,j,k,nnew,indx)=MIN(t(i,j,k,nnew,indx),0.0_r8)+        &
     &                           Hz(i,j,k)*Bio(i,k,indx)
!              t(i,j,k,nnew,indx)= Hz(i,j,k)*Bio(i,k,indx)
#ifdef TS_MPDATA
               t(i,j,k,3,indx)=t(i,j,k,nnew,indx)
#endif
            END DO
          END DO
        END DO
      END DO J_LOOP

      RETURN
      END SUBROUTINE biology_tile
