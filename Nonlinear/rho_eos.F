#include "cppdefs.h"
      MODULE rho_eos_mod
#ifdef SOLVE3D
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine computes  "in situ" density and other associated       !
!  quantitites as a function of potential temperature,  salinity,      !
!  and pressure from a polynomial expression (Jackett & McDougall,     !
!  1992). The polynomial expression was found from fitting to 248      !
!  values  in the  oceanographic  ranges of  salinity,  potential      !
!  temperature,  and pressure.  It  assumes no pressure variation      !
!  along geopotential surfaces, that is, depth (meters; negative)      !
!  and pressure (dbar; assumed negative here) are interchangeable.     !
!                                                                      !
!  Check Values: (T=3 C, S=35.5 PSU, Z=-5000 m)                        !
!                                                                      !
!     alpha = 2.1014611551470d-04 (1/Celsius)                          !
!     beta  = 7.2575037309946d-04 (1/PSU)                              !
!     gamma = 3.9684764511766d-06 (1/Pa)                               !
!     den   = 1050.3639165364     (kg/m3)                              !
!     den1  = 1028.2845117925     (kg/m3)                              !
!     sound = 1548.8815240223     (m/s)                                !
!     bulk  = 23786.056026320     (Pa)                                 !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!  Jackett, D. R. and T. J. McDougall, 1995, Minimal Adjustment of     !
!    Hydrostatic Profiles to Achieve Static Stability, J. of Atmos.    !
!    and Oceanic Techn., vol. 12, pp. 381-389.                         !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC rho_eos

      CONTAINS
!
!***********************************************************************
      SUBROUTINE rho_eos (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_coupling
      USE mod_grid
      USE mod_mixing
      USE mod_ocean
      USE mod_stepping
!
      integer, intent(in) :: ng, tile

# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, 14)
# endif
      CALL rho_eos_tile (ng, Istr, Iend, Jstr, Jend,                    &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   nrhs(ng),                                      &
# ifdef MASKING
     &                   GRID(ng) % rmask,                              &
# endif
     &                   GRID(ng) % Hz,                                 &
     &                   GRID(ng) % z_r,                                &
     &                   GRID(ng) % z_w,                                &
     &                   OCEAN(ng) % t,                                 &
# ifdef VAR_RHO_2D
     &                   COUPLING(ng) % rhoA,                           &
     &                   COUPLING(ng) % rhoS,                           &
# endif
# ifdef BV_FREQUENCY
     &                   MIXING(ng) % bvf,                              &
# endif
# if defined LMD_SKPP || defined LMD_BKPP || defined BULK_FLUXES
     &                   MIXING(ng) % alpha,                            &
     &                   MIXING(ng) % beta,                             &
#  ifdef LMD_DDMIX
     &                   MIXING(ng) % alfaobeta,                        &
#  endif
# endif
# if defined MIX_ISO_TS || defined GENT_McWILLIAMS
     &                   MIXING(ng) % neutral,                          &
# endif
# if defined LMD_SKPP || defined LMD_BKPP
     &                   OCEAN(ng) % pden,                              &
# endif
     &                   OCEAN(ng) % rho)
# ifdef PROFILE
      CALL wclock_off (ng, 14)
# endif
      RETURN
      END SUBROUTINE rho_eos
# ifdef NONLIN_EOS
!
!***********************************************************************
      SUBROUTINE rho_eos_tile (ng, Istr, Iend, Jstr, Jend,              &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         nrhs,                                    &
# ifdef MASKING
     &                         rmask,                                   &
# endif
     &                         Hz, z_r, z_w, t,                         &
# ifdef VAR_RHO_2D
     &                         rhoA, rhoS,                              &
# endif
# ifdef BV_FREQUENCY
     &                         bvf,                                     &
# endif
# if defined LMD_SKPP || defined LMD_BKPP || defined BULK_FLUXES
     &                         alpha, beta,                             &
#  ifdef LMD_DDMIX
     &                         alfaobeta,                               &
#  endif
# endif
# if defined MIX_ISO_TS || defined GENT_McWILLIAMS
     &                         neutral,                                 &
# endif
# if defined LMD_SKPP || defined LMD_BKPP
     &                         pden,                                    &
# endif
     &                         rho)
!***********************************************************************
!
      USE mod_param
      USE mod_eoscoef
      USE mod_scalars
!
# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      USE exchange_2d_mod
      USE exchange_3d_mod
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: nrhs

# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: t(LBi:,LBj:,:,:,:)
#  ifdef VAR_RHO_2D
      real(r8), intent(out) :: rhoA(LBi:,LBj:)
      real(r8), intent(out) :: rhoS(LBi:,LBj:)
#  endif
#  ifdef BV_FREQUENCY
      real(r8), intent(out) :: bvf(LBi:,LBj:,0:)
#  endif
#  if defined LMD_SKPP || defined LMD_BKPP || defined BULK_FLUXES
      real(r8), intent(out) :: alpha(LBi:,LBj:)
      real(r8), intent(out) :: beta(LBi:,LBj:)
#   ifdef LMD_DDMIX
      real(r8), intent(out) :: alfaobeta(LBi:,LBj:,0:)
#   endif
#  endif
#  if defined MIX_ISO_TS || defined GENT_McWILLIAMS
      real(r8), intent(out) :: neutral(LBi:,LBj:,0:)
#  endif
#  if defined LMD_SKPP || defined LMD_BKPP
      real(r8), intent(out) :: pden(LBi:,LBj:,:)
#  endif
      real(r8), intent(out) :: rho(LBi:,LBj:,:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(in) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
#  ifdef VAR_RHO_2D
      real(r8), intent(out) :: rhoA(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: rhoS(LBi:UBi,LBj:UBj)
#  endif
#  ifdef BV_FREQUENCY
      real(r8), intent(out) :: bvf(LBi:UBi,LBj:UBj,0:N(ng))
#  endif
#  if defined LMD_SKPP || defined LMD_BKPP || defined BULK_FLUXES
      real(r8), intent(out) :: alpha(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: beta(LBi:UBi,LBj:UBj)
#   ifdef LMD_DDMIX
      real(r8), intent(out) :: alfaobeta(LBi:UBi,LBj:UBj,0:N(ng))
#   endif
#  endif
#  if defined MIX_ISO_TS || defined GENT_McWILLIAMS
      real(r8), intent(out) :: neutral(LBi:UBi,LBj:UBj,0:N(ng))
#  endif
#  if defined LMD_SKPP || defined LMD_BKPP
      real(r8), intent(out) :: pden(LBi:UBi,LBj:UBj,N(ng))
#  endif
      real(r8), intent(out) :: rho(LBi:UBi,LBj:UBj,N(ng))
# endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j, k

      real(r8) :: DbulkDP, DbulkDS, DbulkDT, Dden1DS, Dden1DT
      real(r8) :: Tp, Tpr10, Ts, Tt
      real(r8) :: bulk_dn, bulk_up, cff1, cff2, den_dn, den_up
      real(r8) :: rhogamma, sound2, sqrtTs

      real(r8), dimension(23) :: cff

      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N(ng)) :: Scof
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N(ng)) :: Tcof
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N(ng)) :: bulk
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N(ng)) :: bulk0
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N(ng)) :: bulk1
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N(ng)) :: bulk2
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N(ng)) :: den
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N(ng)) :: den1
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N(ng)) :: wrk

#  include "set_bounds.h"
!
!=======================================================================
!  Nonlinear equation of state.  Notice that this equation of state
!  is only valid for potential temperature range of -2C to 40C and
!  a salinity range of 0 PSU to 42 PSU.
!=======================================================================
!
      DO j=JstrR,JendR
        DO k=1,N(ng)
          DO i=IstrR,IendR
            Tt=MAX(-2.0_r8,t(i,j,k,nrhs,itemp))
#  ifdef SALINITY
            Ts=MAX(0.0_r8,t(i,j,k,nrhs,isalt))
            sqrtTs=SQRT(Ts)
#  else
            Ts=0.0_r8
            sqrtTs=0.0_r8
#  endif
            Tp=z_r(i,j,k)
            Tpr10=0.1_r8*Tp
!
!-----------------------------------------------------------------------
!  Compute density (kg/m3) at standard one atmosphere pressure.
!-----------------------------------------------------------------------
!
            cff(1)=Q00+Tt*(Q01+Tt*(Q02+Tt*(Q03+Tt*(Q04+Tt*Q05))))
            cff(2)=U00+Tt*(U01+Tt*(U02+Tt*(U03+Tt*U04)))
            cff(3)=V00+Tt*(V01+Tt*V02)
!
            den1(i,k)=cff(1)+Ts*(cff(2)+sqrtTs*cff(3)+Ts*W00)
!
!-----------------------------------------------------------------------
!  Compute secant bulk modulus (bulk = bulk0 - bulk1*z + bulk2*z*z).
!-----------------------------------------------------------------------
!
            cff(4)=A00+Tt*(A01+Tt*(A02+Tt*(A03+Tt*A04)))
            cff(5)=B00+Tt*(B01+Tt*(B02+Tt*B03))
            cff(6)=D00+Tt*(D01+Tt*D02)
            cff(7)=E00+Tt*(E01+Tt*(E02+Tt*E03))
            cff(8)=F00+Tt*(F01+Tt*F02)
            cff(9)=G01+Tt*(G02+Tt*G03)
            cff(10)=H00+Tt*(H01+Tt*H02)
!
            bulk0(i,k)=cff(4)+Ts*(cff(5)+sqrtTs*cff(6))
            bulk1(i,k)=cff(7)+Ts*(cff(8)+sqrtTs*G00)
            bulk2(i,k)=cff(9)+Ts*cff(10)
!
            bulk(i,k)=bulk0(i,k)-Tp*(bulk1(i,k)-Tp*bulk2(i,k))
!
!-----------------------------------------------------------------------
!  Compute "in situ" density anomaly (kg/m3).
!-----------------------------------------------------------------------
!
            den(i,k)=(den1(i,k)*bulk(i,k))/(bulk(i,k)+Tpr10)
            rho(i,j,k)=(den(i,k)-1000.0_r8)
!!          rho(i,j,k)=(den1(i,k)-1000.0_r8)
#  ifdef MASKING
            rho(i,j,k)=rho(i,j,k)*rmask(i,j)
#  endif
#  if defined LMD_SKPP || defined LMD_BKPP
!
!-----------------------------------------------------------------------
!  Compute "in situ" density anomaly (kg/m3) referenced to the
!  surface.
!-----------------------------------------------------------------------
!
            pden(i,j,k)=(den1(i,k)-1000.0_r8)
#   ifdef MASKING
            pden(i,j,k)=pden(i,j,k)*rmask(i,j)
#   endif
#  endif
          END DO
        END DO
#  ifdef VAR_RHO_2D
!
!-----------------------------------------------------------------------
!  Compute vertical averaged density and density perturbation used
!  in barotropic pressure gradient.
!-----------------------------------------------------------------------
!
        DO i=IstrR,IendR
          cff1=rho(i,j,N(ng))*Hz(i,j,N(ng))
          rhoS(i,j)=0.5_r8*cff1*Hz(i,j,N(ng))
          rhoA(i,j)=cff1
        END DO
        DO k=N(ng)-1,1,-1
          DO i=IstrR,IendR
            cff1=rho(i,j,k)*Hz(i,j,k)
            rhoS(i,j)=rhoS(i,j)+Hz(i,j,k)*(rhoA(i,j)+0.5_r8*cff1)
            rhoA(i,j)=rhoA(i,j)+cff1
          END DO
        END DO
        cff2=1.0_r8/rho0
        DO i=IstrR,IendR
          cff1=1.0_r8/(z_w(i,j,N(ng))-z_w(i,j,0))
          rhoA(i,j)=cff2*cff1*rhoA(i,j)
          rhoS(i,j)=2.0_r8*cff1*cff1*cff2*rhoS(i,j)
        END DO
#  endif
#  ifdef BV_FREQUENCY
!
!-----------------------------------------------------------------------
!  Compute Brunt-Vaisala frequency (1/s2) at horizontal RHO-points
!  and vertical W-points:
!
!                  bvf = - g/rho d(rho)/d(z).
!
!  The density anomaly difference is computed by lowering/rising the
!  water parcel above/below adiabatically at W-point depth "z_w".
!-----------------------------------------------------------------------
!
        DO k=1,N(ng)-1
          DO i=IstrR,IendR
            bulk_up=bulk0(i,k+1)-                                       &
     &              z_w(i,j,k)*(bulk1(i,k+1)-                           &
     &                          bulk2(i,k+1)*z_w(i,j,k))
            bulk_dn=bulk0(i,k  )-                                       &
     &              z_w(i,j,k)*(bulk1(i,k  )-                           &
     &                          bulk2(i,k  )*z_w(i,j,k))
            den_up=(den1(i,k+1)*bulk_up)/(bulk_up+0.1_r8*z_w(i,j,k))
            den_dn=(den1(i,k  )*bulk_dn)/(bulk_dn+0.1_r8*z_w(i,j,k))
            bvf(i,j,k)=-g*(den_up-den_dn)/                              &
     &                 (0.5_r8*(den_up+den_dn)*                         &
     &                  (z_r(i,j,k+1)-z_r(i,j,k)))
          END DO
        END DO
        DO i=IstrR,IendR
!!        bvf(i,j,0)=bvf(i,j,1)
!!        bvf(i,j,N(ng))=bvf(i,j,N(ng)-1)
          bvf(i,j,0)=0.0_r8
          bvf(i,j,N(ng))=0.0_r8
        END DO
#  endif
#  if defined LMD_SKPP || defined LMD_BKPP || defined BULK_FLUXES
!
!-----------------------------------------------------------------------
!  Compute thermal expansion (1/Celsius) and saline contraction
!  (1/PSU) coefficients.
!-----------------------------------------------------------------------
!
#   ifdef LMD_DDMIX
        DO k=1,N(ng)
#   else
        k=N(ng)
#   endif /* LMD_DDMIX */
          DO i=IstrR,IendR
            Tt=MAX(-2.0_r8,t(i,j,k,nrhs,itemp))
#  ifdef SALINITY
            Ts=MAX(0.0_r8,t(i,j,k,nrhs,isalt))
            sqrtTs=SQRT(Ts)
#   else
            Ts=0.0_r8
            sqrtTs=0.0_r8
#   endif
            Tp=z_r(i,j,k)
            Tpr10=0.1_r8*Tp
!
!  Compute d(den1)/d(S) and d(den1)/d(T) derivatives.
!
            cff(2)=U00+Tt*(U01+Tt*(U02+Tt*(U03+Tt*U04)))
            cff(3)=V00+Tt*(V01+Tt*V02)
            cff(11)=Q01+Tt*(2.0_r8*Q02+                                 &
     &                      Tt*(3.0_r8*Q03+                             &
     &                          Tt*(4.0_r8*Q04+Tt*5.0_r8*Q05)))
            cff(12)=U01+Tt*(2.0_r8*U02+                                 &
     &                      Tt*(3.0_r8*U03+Tt*4.0_r8*U04))
            cff(13)=V01+Tt*2.0_r8*V02
!
            Dden1DS=cff(2)+1.5_r8*cff(3)*sqrtTs+2.0_r8*W00*Ts
            Dden1DT=cff(11)+Ts*(cff(12)+sqrtTs*cff(13))
!
!  Compute d(bulk)/d(S) and d(bulk)/d(T) derivatives.
!
            cff(5)=B00+Tt*(B01+Tt*(B02+Tt*B03))
            cff(6)=D00+Tt*(D01+Tt*D02)
            cff(8)=F00+Tt*(F01+Tt*F02)
            cff(10)=H00+Tt*(H01+Tt*H02)
            cff(14)=A01+Tt*(2.0_r8*A02+Tt*(3.0_r8*A03+Tt*4.0_r8*A04))
            cff(15)=B01+Tt*(2.0_r8*B02+Tt*3.0_r8*B03)
            cff(16)=D01+Tt*2.0_r8*D02
            cff(17)=E01+Tt*(2.0_r8*E02+Tt*3.0_r8*E03)
            cff(18)=F01+Tt*2.0_r8*F02
            cff(19)=G02+Tt*2.0_r8*G03
            cff(20)=H01+Tt*2.0_r8*H02
!
            DbulkDS=cff(5)+sqrtTs*1.5_r8*cff(6)-                        &
     &              Tp*(cff(8)+sqrtTs*1.5_r8*G00-Tp*cff(10))
            DbulkDT=cff(14)+Ts*(cff(15)+sqrtTs*cff(16))-                &
     &              Tp*(cff(17)+Ts*cff(18)-Tp*(cff(19)+Ts*cff(20)))
!
!  Compute thermal expansion and saline contraction coefficients.
!
            cff(21)=Tpr10*den1(i,k)
            cff(22)=bulk(i,k)*(bulk(i,k)+Tpr10)
            wrk(i,k)=den(i,k)*(bulk(i,k)+Tpr10)*(bulk(i,k)+Tpr10)
            Tcof(i,k)=-(DbulkDT*cff(21)+Dden1DT*cff(22))
            Scof(i,k)= (DbulkDS*cff(21)+Dden1DS*cff(22))
#   ifdef LMD_DDMIX
            alfaobeta(i,j,k)=Tcof(i,k)/Scof(i,k)
#   endif
          END DO
          IF (k.eq.N(ng)) THEN
            DO i=IstrR,IendR
              alpha(i,j)=Tcof(i,N(ng))/wrk(i,N(ng))
              beta (i,j)=Scof(i,N(ng))/wrk(i,N(ng))
            END DO
          END IF
#   ifdef LMD_DDMIX
        END DO
#   endif
#  endif
#  if defined MIX_ISO_TS || defined GENT_McWILLIAMS
!
!-----------------------------------------------------------------------
!  Compute factor to convert isopycnal, "in situ" density slopes to
!  to isopycnal, neutral slopes:
!
!    neutral = [1 + g^2 / (c^2 * bvf)];    c^2 = 10000/(rho*gamma)
!
!  were "c" is the sound speed (m/s) and "gamma" is the adiabatic and
!  isentropic compressibility coefficient (1/Pa).
!-----------------------------------------------------------------------
!
        DO k=1,N(ng)
          DO i=IstrR,IendR
            Tp=z_r(i,j,k)
            Tpr10=0.1_r8*Tp
            cff(23)=(bulk(i,k)+Tpr10)*(bulk(i,k)+Tpr10)
            DbulkDP=-bulk1(i,k)+Tp*2.0_r8*bulk2(i,k)
            rhogamma=ABS(0.1_r8*den1(i,k)*                              &
     &                   ((Tp*DbulkDP-bulk(i,k))/cff(23)))
            sound2=10000.0_r8/rhogamma
            neutral(i,j,k)=1.0_r8+g*g/(sound2*bvf(i,j,k))
          END DO
        END DO
#  endif
      END DO
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
!
!------------------------------------------------------------------------
!  Exchange boundary data.
!------------------------------------------------------------------------
!
      CALL exchange_r3d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, rho)
#   if defined LMD_SKPP || defined LMD_BKPP
      CALL exchange_r3d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, pden)
#   endif
#   if defined LMD_SKPP || defined LMD_BKPP || defined BULK_FLUXES
#    ifdef LMD_DDMIX
      CALL exchange_w3d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, alfaobeta)
#    endif
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, alpha)
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, beta)
#   endif    
#   if defined MIX_ISO_TS || defined GENT_McWILLIAMS
      CALL exchange_r3d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, neutral)
#   endif
#   ifdef VAR_RHO_2D
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, rhoA)
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, rhoS)
#   endif
#   ifdef BV_FREQUENCY
      CALL exchange_w3d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, bvf)
#   endif
#  endif
      RETURN
      END SUBROUTINE rho_eos_tile
# endif
# ifndef NONLIN_EOS
!
!***********************************************************************
      SUBROUTINE rho_eos_tile (ng, Istr, Iend, Jstr, Jend,              &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         nrhs,                                    &
# ifdef MASKING
     &                         rmask,                                   &
# endif
     &                         Hz, z_r, z_w, t,                         &
# ifdef VAR_RHO_2D
     &                         rhoA, rhoS,                              &
# endif
# ifdef BV_FREQUENCY
     &                         bvf,                                     &
# endif
# if defined LMD_SKPP || defined LMD_BKPP || defined BULK_FLUXES
     &                         alpha, beta,                             &
#  ifdef LMD_DDMIX
     &                         alfaobeta,                               &
#  endif
# endif
# if defined MIX_ISO_TS || defined GENT_McWILLIAMS
     &                         neutral,                                 &
# endif
# if defined LMD_SKPP || defined LMD_BKPP
     &                         pden,                                    &
# endif
     &                         rho)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      USE exchange_2d_mod
      USE exchange_3d_mod
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: nrhs

# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: t(LBi:,LBj:,:,:,:)

#  ifdef VAR_RHO_2D
      real(r8), intent(out) :: rhoA(LBi:,LBj:)
      real(r8), intent(out) :: rhoS(LBi:,LBj:)
#  endif
#  ifdef BV_FREQUENCY
      real(r8), intent(out) :: bvf(LBi:,LBj:,0:)
#  endif
#  if defined LMD_SKPP || defined LMD_BKPP || defined BULK_FLUXES
      real(r8), intent(out) :: alpha(LBi:,LBj:)
      real(r8), intent(out) :: beta(LBi:,LBj:)
#   ifdef LMD_DDMIX
      real(r8), intent(out) :: alfaobeta(LBi:,LBj:,0:)
#   endif
#  endif
#  if defined MIX_ISO_TS || defined GENT_McWILLIAMS
      real(r8), intent(out) :: neutral(LBi:,LBj:,0:)
#  endif
#  if defined LMD_SKPP || defined LMD_BKPP
      real(r8), intent(out) :: pden(LBi:,LBj:,:)
#  endif
      real(r8), intent(out) :: rho(LBi:,LBj:,:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(in) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
#  ifdef VAR_RHO_2D
      real(r8), intent(out) :: rhoA(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: rhoS(LBi:UBi,LBj:UBj)
#  endif
#  ifdef BV_FREQUENCY
      real(r8), intent(out) :: bvf(LBi:UBi,LBj:UBj,0:N(ng))
#  endif
#  if defined LMD_SKPP || defined LMD_BKPP || defined BULK_FLUXES
      real(r8), intent(out) :: alpha(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: beta(LBi:UBi,LBj:UBj)
#   ifdef LMD_DDMIX
      real(r8), intent(out) :: alfaobeta(LBi:UBi,LBj:UBj,0:N(ng))
#   endif
#  endif
#  if defined MIX_ISO_TS || defined GENT_McWILLIAMS
      real(r8), intent(out) :: neutral(LBi:UBi,LBj:UBj,0:N(ng))
#  endif
#  if defined LMD_SKPP || defined LMD_BKPP
      real(r8), intent(out) :: pden(LBi:UBi,LBj:UBj,N(ng))
#  endif
      real(r8), intent(out) :: rho(LBi:UBi,LBj:UBj,N(ng))
# endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j, k
      real(r8) :: cff, cff1, cff2

#  include "set_bounds.h"
!
!=======================================================================
!  Linear equation of state.
!=======================================================================
!
!-----------------------------------------------------------------------
!  Compute "in situ" density anomaly (kg/m3).
!-----------------------------------------------------------------------
!
      DO j=JstrR,JendR
        DO k=1,N(ng)
          DO i=IstrR,IendR
            rho(i,j,k)=R0(ng)+Tcoef(ng)*t(i,j,k,nrhs,itemp)
#  ifdef SALINITY
            rho(i,j,k)=rho(i,j,k)+Scoef(ng)*t(i,j,k,nrhs,isalt)
#  endif
#  ifdef MASKING
            rho(i,j,k)=rho(i,j,k)*rmask(i,j)
#  endif
#  if defined LMD_SKPP || defined LMD_BKPP
            pden(i,j,k)=rho(i,j,k)
#  endif
          END DO
        END DO
#  ifdef VAR_RHO_2D
!
!-----------------------------------------------------------------------
!  Compute vertical averaged density and density perturbation used
!  in barotropic pressure gradient.
!-----------------------------------------------------------------------
!
        DO i=IstrR,IendR
          cff1=rho(i,j,N(ng))*Hz(i,j,N(ng))
          rhoS(i,j)=0.5_r8*cff1*Hz(i,j,N(ng))
          rhoA(i,j)=cff1
        END DO
        DO k=N(ng)-1,1,-1
          DO i=IstrR,IendR
            cff1=rho(i,j,k)*Hz(i,j,k)
            rhoS(i,j)=rhoS(i,j)+Hz(i,j,k)*(rhoA(i,j)+0.5_r8*cff1)
            rhoA(i,j)=rhoA(i,j)+cff1
          END DO
        END DO
        cff2=1.0_r8/rho0
        DO i=IstrR,IendR
          cff1=1.0_r8/(z_w(i,j,N(ng))-z_w(i,j,0))
          rhoA(i,j)=cff2*cff1*rhoA(i,j)
          rhoS(i,j)=2.0_r8*cff1*cff1*cff2*rhoS(i,j)
        END DO
#  endif
#  ifdef BV_FREQUENCY
!
!-----------------------------------------------------------------------
!  Compute Brunt-Vaisala frequency (1/s2) at horizontal RHO-points
!  and vertical W-points.
!-----------------------------------------------------------------------
!
        DO k=1,N(ng)-1
          DO i=IstrR,IendR
            bvf(i,j,k)=-gorho0*(rho(i,j,k+1)-rho(i,j,k))/               &
     &                         (z_r(i,j,k+1)-z_r(i,j,k))
          END DO
        END DO
#  endif
#  if defined LMD_SKPP || defined LMD_BKPP || defined BULK_FLUXES
!
!-----------------------------------------------------------------------
!  Compute thermal expansion (1/Celsius) and saline contraction
!  (1/PSU) coefficients.
!-----------------------------------------------------------------------
!
        DO i=IstrR,IendR
          alpha(i,j)=ABS(Tcoef(ng))
#   ifdef SALINITY
          beta(i,j)=ABS(Scoef(ng))
#   else
          beta(i,j)=0.0_r8
#   endif
        END DO
#   ifdef LMD_DDMIX
!
!  Compute ratio of thermal expansion and saline contraction
!  coefficients.
!
        IF (Scoef(ng).eq.0.0_r8) THEN
          cff=1.0_r8
        ELSE
          cff=1.0_r8/Scoef(ng)
        END IF
        DO k=1,N(ng)
          DO i=IstrR,IendR
            alfaobeta(i,j,k)=cff*Tcoef(ng)
          END DO
        END DO
#   endif
#  endif
      END DO
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
!
!------------------------------------------------------------------------
!  Exchange boundary data.
!------------------------------------------------------------------------
!
      CALL exchange_r3d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, rho)
#   if defined LMD_SKPP || defined LMD_BKPP
      CALL exchange_r3d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, pden)
#   endif
#   if defined LMD_SKPP || defined LMD_BKPP || defined BULK_FLUXES
#    ifdef LMD_DDMIX
      CALL exchange_w3d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, alfaobeta)
#    endif
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, alpha)
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, beta)
#   endif    
#   if defined MIX_ISO_TS || defined GENT_McWILLIAMS
      CALL exchange_r3d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, neutral)
#   endif
#   ifdef VAR_RHO_2D
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, rhoA)
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, rhoS)
#   endif
#   ifdef BV_FREQUENCY
      CALL exchange_w3d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, bvf)
#   endif
#  endif
      RETURN
      END SUBROUTINE rho_eos_tile
# endif
#endif
      END MODULE rho_eos_mod
