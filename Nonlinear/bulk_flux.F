#include "cppdefs.h"
      MODULE bulk_flux_mod
#ifdef BULK_FLUXES
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine computes the bulk parameterization of surface wind     !
!  stress and surface net heat fluxes.                                 !
!                                                                      !
!  References:                                                         !
!                                                                      !
!    Fairall, C.W., E.F. Bradley, D.P. Rogers, J.B. Edson and G.S.     !
!      Young, 1996:  Bulk parameterization of air-sea fluxes for       !
!      tropical ocean-global atmosphere Coupled-Ocean Atmosphere       !
!      Response Experiment, JGR, 101, 3747-3764.                       !
!                                                                      !
!    Fairall, C.W., E.F. Bradley, J.S. Godfrey, G.A. Wick, J.B.        !
!      Edson, and G.S. Young, 1996:  Cool-skin and warm-layer          !
!      effects on sea surface temperature, JGR, 101, 1295-1308.        !
!                                                                      !
!    Liu, W.T., K.B. Katsaros, and J.A. Businger, 1979:  Bulk          !
!        parameterization of the air-sea exchange of heat and          !
!        water vapor including the molecular constraints at            !
!        the interface, J. Atmos. Sci, 36, 1722-1735.                  !
!                                                                      !
!  Adapted from COARE code written originally by David Rutgers and     !
!  Frank Bradley.                                                      !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC bulk_flux

      CONTAINS
!
!***********************************************************************
      SUBROUTINE bulk_flux (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_forces
      USE mod_grid
      USE mod_mixing
      USE mod_ocean
      USE mod_stepping
!
      integer, intent(in) :: ng, tile

# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, 17)
# endif
      CALL bulk_flux_tile (ng, Istr, Iend, Jstr, Jend,                  &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     nrhs(ng),                                    &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
     &                     GRID(ng) % umask,                            &
     &                     GRID(ng) % vmask,                            &
# endif
     &                     MIXING(ng) % alpha,                          &
     &                     MIXING(ng) % beta,                           &
     &                     OCEAN(ng) % rho,                             &
     &                     OCEAN(ng) % t,                               &
     &                     FORCES(ng) % Hair,                           &
     &                     FORCES(ng) % Pair,                           &
     &                     FORCES(ng) % Tair,                           &
     &                     FORCES(ng) % Uwind,                          &
     &                     FORCES(ng) % Vwind,                          &
     &                     FORCES(ng) % cloud,                          &
     &                     FORCES(ng) % rain,                           &
     &                     FORCES(ng) % lhflx,                          &
     &                     FORCES(ng) % lrflx,                          &
     &                     FORCES(ng) % shflx,                          &
     &                     FORCES(ng) % srflx,                          &
     &                     FORCES(ng) % stflx,                          &
     &                     FORCES(ng) % sustr,                          &
     &                     FORCES(ng) % svstr)
# ifdef PROFILE
      CALL wclock_off (ng, 17)
# endif
      RETURN
      END SUBROUTINE bulk_flux
!
!***********************************************************************
      SUBROUTINE bulk_flux_tile (ng, Istr, Iend, Jstr, Jend,            &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           nrhs,                                  &
# ifdef MASKING
     &                           rmask, umask, vmask,                   &
# endif
     &                           alpha, beta, rho, t,                   &
     &                           Hair, Pair, Tair, Uwind, Vwind,        &
     &                           cloud, rain, lhflx, lrflx, shflx,      &
     &                           srflx, stflx, sustr, svstr)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars

# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
!
      USE exchange_2d_mod
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
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: alpha(LBi:,LBj:)
      real(r8), intent(in) :: beta(LBi:,LBj:)
      real(r8), intent(in) :: rho(LBi:,LBj:,:)
      real(r8), intent(in) :: t(LBi:,LBj:,:,:,:)
      real(r8), intent(in) :: Hair(LBi:,LBj:)
      real(r8), intent(in) :: Pair(LBi:,LBj:)
      real(r8), intent(in) :: Tair(LBi:,LBj:)
      real(r8), intent(in) :: Uwind(LBi:,LBj:)
      real(r8), intent(in) :: Vwind(LBi:,LBj:)
      real(r8), intent(in) :: cloud(LBi:,LBj:)
      real(r8), intent(in) :: rain(LBi:,LBj:)

      real(r8), intent(inout) :: lhflx(LBi:,LBj:)
      real(r8), intent(inout) :: lrflx(LBi:,LBj:)
      real(r8), intent(inout) :: shflx(LBi:,LBj:)
      real(r8), intent(inout) :: srflx(LBi:,LBj:)
      real(r8), intent(inout) :: stflx(LBi:,LBj:,:)

      real(r8), intent(out) :: sustr(LBi:,LBj:)
      real(r8), intent(out) :: svstr(LBi:,LBj:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: alpha(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: beta(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: rho(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(in) :: Hair(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Pair(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Tair(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Uwind(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Vwind(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: cloud(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: rain(LBi:UBi,LBj:UBj)

      real(r8), intent(inout) :: lhflx(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: lrflx(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: shflx(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: srflx(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: stflx(LBi:UBi,LBj:UBj,NT(ng))

      real(r8), intent(out) :: sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: svstr(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: Iter, i, j, k

      integer, parameter :: IterMax = 20

      real(r8), parameter :: r3 = 1.0_r8/3.0_r8

      real(r8) :: Bf, Cd, Ce, Ch, Hl, Hlv, Hlw, Hscale, Hs, Hsr, IER
      real(r8) :: PairM, Q, Qair, Qpsi, Qsea, Qstar, RH, Ri, Rq, Rr
      real(r8) :: Rt, Scff, TairC, TairK, Taur, Tcff, Tpsi, TseaC
      real(r8) :: TseaK, Tstar, TVstar, VisAir, Wgus, Wmag, Wpsi
      real(r8) :: Wspeed, Wstar, ZQ, ZQoL, ZT, ZToL, ZW, ZWoL, ZoQ
      real(r8) :: ZoT, ZoW

      real(r8) :: cff, cff1, cff2, delQ, delQc, delT, delTc, delW
      real(r8) :: diffh, diffw, oL, rainfall, rhoAir, rhoSea, upvel
      real(r8) :: wet_bulb
# ifdef LONGWAVE
      real(r8) :: e_sat, vap_p
# endif
# ifdef COOL_SKIN
      real(r8) :: Clam, Cwet, Fc, Hcool, Hsb, Hlb, Qbouy, Qcool, lambd
# endif

      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: LHeat
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: LRad
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: SHeat
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: SRad
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: Taux
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: Tauy

# include "set_bounds.h"
!
!=======================================================================
!  Atmosphere-Ocean bulk fluxes parameterization.
!=======================================================================
!
      Hscale=rho0*Cp
      DO j=MIN(Jstr-1,JstrR),JendR
        DO i=MIN(Istr-1,IstrR),IendR
!
!  Input bulk parameterization fields.
!
          Wmag=SQRT(Uwind(i,j)*Uwind(i,j)+Vwind(i,j)*Vwind(i,j))
          PairM=Pair(i,j)
          TairC=Tair(i,j)
          TairK=TairC+273.16_r8
          TseaC=t(i,j,N(ng),nrhs,itemp)
          TseaK=TseaC+273.16_r8
          rhoSea=rho(i,j,N(ng))+1000.0_r8
          RH=Hair(i,j)
          SRad(i,j)=srflx(i,j)*Hscale
          rainfall=rain(i,j)
          Tcff=alpha(i,j)
          Scff=beta(i,j)
          ZT=blk_ZT
          ZQ=blk_ZQ
          ZW=blk_ZW
!
!  Initialize.
!
          IER=0
          delTc=0.0_r8
          delQc=0.0_r8
          LHeat(i,j)=lhflx(i,j)*Hscale
          SHeat(i,j)=shflx(i,j)*Hscale
          Taur=0.0_r8
          Taux(i,j)=0.0_r8
          Tauy(i,j)=0.0_r8
!
!-----------------------------------------------------------------------
!  Compute net longwave radiation (W/m2), LRad.
!-----------------------------------------------------------------------
!
# ifdef LONGWAVE
!  Use Berliand (1952) formula.
!
          cff=(0.7859_r8+0.03477_r8*TairC)/(1.0_r8+0.00412_r8*TairC)
          e_sat=10.0_r8**cff
          vap_p=e_sat*RH
          cff2=TairK*TairK*TairK
          cff1=cff2*TairK
          LRad(i,j)=-emmiss*StefBo*                                     &
     &              (cff1*(0.39_r8-0.5_r8*SQRT(0.01_r8*vap_p))*         &
     &                    (1.0_r8-0.6823_r8*cloud(i,j)*cloud(i,j))+     &
     &               cff2*4.0_r8*(TseaK-TairK))
# else
          LRad(i,j)=lrflx(i,j)*Hscale
# endif
# ifdef MASKING
          LRad(i,j)=LRad(i,j)*rmask(i,j)
# endif
!
!-----------------------------------------------------------------------
!  Compute specific humidities (kg/kg).
!-----------------------------------------------------------------------
!
!  Compute air saturation vapor pressure (mb), using Teten formula.
!
          Qair=(1.0007_r8+3.46E-6_r8*PairM)*6.1121_r8*                  &
     &         EXP(17.502_r8*TairC/(240.97_r8+TairC))
!
!  Compute specific humdity, Q (kg/kg).  Check whether input relative
!  humidity is in kg/kg or g/kg.
!
          IF (RH.lt.2.0_r8) THEN
            cff=Qair*RH
            Q=0.62197_r8*(cff/(PairM-0.378_r8*cff))
          ELSE
            Q=RH/1000.0_r8
          END IF
          Qair=0.62197_r8*(Qair/(PairM-0.378_r8*Qair))
!
!  Compute water saturation vapor pressure (mb), using Teten formula.
!
          Qsea=(1.0007_r8+3.46E-6_r8*PairM)*6.1121_r8*                  &
     &         EXP(17.502_r8*TseaC/(240.97_r8+TseaC))
!
!  Qsea reduced for salinity (Kraus, 1972, pp 46).
!
          Qsea=Qsea*0.98_r8
!
!  Convert Qsea from mb to specific humidity (kg/kg).
!
          Qsea=0.62197_r8*(Qsea/(PairM-0.378_r8*Qsea))
!
!-----------------------------------------------------------------------
!  Compute Monin-Obukhov similarity parameters for wind (Wstar),
!  heat (Tstar), and moisture (Qstar), Liu et al. (1979).
!-----------------------------------------------------------------------
!
!  Moist air density (kg/m3).
!
          rhoAir=PairM*100.0_r8/(blk_Rgas*TairK*(1.0_r8+0.61_r8*Q))
!
!  Kinematic viscosity of dry air (m2/s), Andreas (1989).
!
          VisAir=1.326E-5_r8*(1.0_r8+TairC*(6.542E-3_r8+TairC*          &
     &                        (8.301E-6_r8-4.84E-9_r8*TairC)))
!
!  Compute latent heat of vaporization (J/kg) at sea surface, Hlv.
!
          Hlv=(2.501_r8-0.00237_r8*TseaC)*1.0E+6_r8
# ifdef COOL_SKIN
!
!  Cool skin correction constants. Clam: part of Saunders constant
!  lambda; Cwet: slope of saturation vapor.
!
          Clam=16.0_r8*g*blk_Cpw*(rhoSea*blk_visw)**3/                  &
     &         (blk_tcw*blk_tcw*rhoAir*rhoAir)
          Cwet=0.622_r8*Hlv*Qsea/(blk_Rgas*TseaK*TseaK)
!
!  Set initial guesses for cool-skin layer thickness (Hcool).
!
          Hcool=0.001_r8
# endif
!
!  Assume that wind is measured relative to sea surface and include
!  gustiness.
!
          Wgus=0.5_r8
          delW=SQRT(Wmag*Wmag+Wgus*Wgus)
          delQ=Q-Qsea
          delT=TairC-TseaC
!
!  Initial guesses for Monon-Obukhov similarity scales.
!
          ZWoL=0.0_r8
          ZoW=0.0005_r8
          Wstar=0.04_r8*delW
          Tstar=0.04_r8*delT
          Qstar=0.04_r8*delQ
          TVstar=Tstar*(1.0_r8+0.61_r8*Q)+0.61_r8*TairK*Qstar
!
!  Compute Richardson number.
!
          IF (delW.ne.0.0_r8) THEN
            Ri=g*ZW*(delT+0.61_r8*TairK*delQ)/(TairK*delW*delW)
!           if (Ri.gt.0.25_r8) IER=-1
          ELSE
            IER=-1
          END IF
!
!  Iterate until convergence when IER is negative.  It usually
!  converges within four iterations.
!
          DO Iter=1,IterMax
            IF (IER.ge.0) THEN
!
!  Compute Monin-Obukhov stability parameter, Z/L.
!
              oL=g*vonKar*TVstar/                                       &
     &           (TairK*(1.0_r8+0.61_r8*Q)*Wstar*Wstar)
              ZWoL=ZW*oL
              ZToL=ZT*oL
              ZQoL=ZQ*oL
!
!  Evaluate stability functions at Z/L.
!
              Wpsi=bulk_psi(1,ZWoL)
              Tpsi=bulk_psi(2,ZToL)
              Qpsi=bulk_psi(2,ZQoL)
!
!  Compute wind scaling parameters, Wstar.
!
              ZoW=0.011_r8*Wstar*Wstar/g+0.11_r8*VisAir/Wstar
              Wstar=delW*vonKar/(LOG(ZW/ZoW)-Wpsi)
!
!  Computes roughness Reynolds number for wind (Rr), heat (Rt),
!  and moisture (Rq). Use Liu et al. (1976) look-up table to
!  compute "Rt" and "Rq" as function of "Rr".
!
              Rr=ZoW*Wstar/VisAir
              IF ((Rr.ge.0.0_r8).and.(Rr.lt.1000.0_r8)) THEN
                DO k=1,8
                  IF ((Liu_Rr(k).le.Rr).and.(Rr.lt.Liu_Rr(k+1))) THEN
                    Rt=Liu_a(k,1)*Rr**Liu_b(k,1)
                    Rq=Liu_a(k,2)*Rr**Liu_b(k,2)
                  END IF
                END DO
!
!  Compute heat and moisture scaling parameters, Tstar and Qstar.
!
                cff=VisAir/Wstar
                ZoT=Rt*cff
                ZoQ=Rq*cff
                cff=vonKar*blk_fdg
                Tstar=(delT+delTc)*cff/(LOG(ZT/ZoT)-Tpsi)
                Qstar=(delQ+delQc)*cff/(LOG(ZQ/ZoQ)-Qpsi)
!
!  Compute gustiness in wind speed.
!
                TVstar=Tstar*(1.0_r8+0.61_r8*Q)+0.61_r8*TairK*Qstar
                Bf=-g/TairK*Wstar*TVstar
                IF (Bf.gt.0) THEN
                  Wgus=blk_beta*(Bf*blk_Zabl)**r3
                ELSE
                  Wgus=0.0_r8
                END IF
                delW=SQRT(Wmag*Wmag+Wgus*Wgus)
# ifdef COOL_SKIN
!
!-----------------------------------------------------------------------
!  Cool Skin correction.
!-----------------------------------------------------------------------
!
!  Backgound sensible and latent heat.
!
                Hsb=-rhoAir*blk_Cpa*Wstar*Tstar
                Hlb=-rhoAir*Hlv*Wstar*Qstar
!
!  Mean absoption in cool-skin layer.
!
                Fc=0.137_r8+11.0_r8*Hcool-                              &
     &             (1.0_r8-EXP(-Hcool/8.0e-4))*6.6E-5_r8/Hcool
!
!  Total cooling at the interface.
!
                Qcool=LRad(i,j)+Hsb+Hlb-SRad(i,j)*Fc
                Qbouy=Tcff*Qcool+Scff*Hlb*blk_Cpw/Hlv
!
!  Compute temperature and moisture change.
!
                IF ((Qcool.gt.0.0_r8).and.(Qbouy.gt.0.0_r8)) THEN
                  lambd=6.0_r8/                                         &
     &                  (1.0_r8+(Clam*Qbouy/Wstar**4)**0.75_r8)**r3
                  Hcool=lambd*blk_visw/(SQRT(rhoAir/rhoSea)*Wstar)
                  delTc=Qcool*Hcool/blk_tcw
                ELSE
                  delTc=0.0_r8
                END IF
                delQc=Cwet*delTc
# endif
              ELSE
                IER=-2
              END IF
            END IF
          END DO
!
!---------------------------------------------------------------------
!  Compute Atmosphere/Ocean fluxes.
!---------------------------------------------------------------------
!
          IF (IER.ge.0) THEN
!
!  Compute transfer coefficients for momentun (Cd), heat (Ch),
!  and moisture (Ce).
!
            Wspeed=SQRT(Wmag*Wmag+Wgus*Wgus)
            Cd=Wstar*Wstar/(Wspeed*Wspeed)
            Ch=Wstar*Tstar/(Wspeed*(TairC-TseaC+delTc))
            Ce=Wstar*Qstar/(Wspeed*(Q-Qsea+delQc))
!
!  Compute turbulent sensible heat flux (W/m2), Hs.
!
            Hs=-blk_Cpa*rhoAir*Wstar*Tstar
!
!  Compute sensible heat flux (W/m2) due to rainfall (kg/m2/s), Hsr.
!
            diffw=2.11E-5_r8*(TairK/273.16_r8)**1.94_r8
            diffh=0.02411_r8*                                           &
     &            (1.0_r8+TairC*(3.309E-3_r8-1.44E-6_r8*TairC))/        &
     &            (rhoAir*blk_Cpa)
            cff=Qair*Hlv/(blk_Rgas*TairK*TairK)
            wet_bulb=1.0_r8/(1.0_r8+0.622_r8*(cff*Hlv*diffw)/           &
     &                                       (blk_Cpa*diffh))
            Hsr=rainfall*wet_bulb*blk_Cpw*                              &
     &          ((TseaC-TairC)+(Qsea-Q)*Hlv/blk_Cpa)
            SHeat(i,j)=(Hs+Hsr)
# ifdef MASKING
            SHeat(i,j)=SHeat(i,j)*rmask(i,j)
# endif
!
!  Compute turbulent latent heat flux (W/m2), Hl.
!
            Hl=-Hlv*rhoAir*Wstar*Qstar
!
!  Compute Webb correction (Webb effect) to latent heat flux, Hlw.
!
            upvel=-1.61_r8*Wstar*Qstar-                                 &
     &            (1.0_r8+1.61_r8*Q)*Wstar*Tstar/TairK
            Hlw=rhoAir*Hlv*upvel*Q
            LHeat(i,j)=(Hl+Hlw)
# ifdef MASKING
            LHeat(i,j)=LHeat(i,j)*rmask(i,j)
# endif
!
!  Compute momentum flux (N/m2) due to rainfall (kg/m2/s).
!
            Taur=0.85_r8*rainfall*Wmag
!
!  Compute wind stress components (N/m2), Tau.
!
            cff=rhoAir*Cd*Wspeed
            Taux(i,j)=(cff*Uwind(i,j)+Taur*SIGN(1.0_r8,Uwind(i,j)))
# ifdef MASKING
            Taux(i,j)=Taux(i,j)*rmask(i,j)
# endif
            Tauy(i,j)=(cff*Vwind(i,j)+Taur*SIGN(1.0_r8,Vwind(i,j)))
# ifdef MASKING
            Tauy(i,j)=Tauy(i,j)*rmask(i,j)
# endif
          END IF
        END DO
      END DO
!
!=======================================================================
!  Compute surface net heat flux and surface wind stress.
!=======================================================================
!
!  Compute kinematic, surface, net heat flux (degC m/s).  Notice that
!  the signs of latent and sensible fluxes are reversed because fluxes
!  from the bulk formulation are positive out of the ocean.
!
      Hscale=1.0_r8/(rho0*Cp)
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          lrflx(i,j)=LRad(i,j)*Hscale
          lhflx(i,j)=-LHeat(i,j)*Hscale
          shflx(i,j)=-SHeat(i,j)*Hscale
          stflx(i,j,itemp)=(srflx(i,j)+lrflx(i,j)+                      &
     &                      lhflx(i,j)+shflx(i,j))
# ifdef MASKING
          stflx(i,j,itemp)=stflx(i,j,itemp)*rmask(i,j)
# endif
        END DO
      END DO
!
!  Compute kinematic, surface wind stress (m2/s2).
!
      cff=0.5_r8/rho0
      DO j=JstrR,JendR
        DO i=Istr,IendR
          sustr(i,j)=cff*(Taux(i-1,j)+Taux(i,j))
# ifdef MASKING
          sustr(i,j)=sustr(i,j)*umask(i,j)
# endif
        END DO
      END DO
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          svstr(i,j)=cff*(Tauy(i,j-1)+Tauy(i,j))
# ifdef MASKING
          svstr(i,j)=svstr(i,j)*vmask(i,j)
# endif
        END DO
      END DO
# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Apply periodic boundary conditions.
!-----------------------------------------------------------------------
!
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, lrflx)
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, lhflx)
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, shflx)
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, stflx(:,:,itemp))
      CALL exchange_u2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, sustr)
      CALL exchange_v2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, svstr)
# endif
      RETURN
      END SUBROUTINE bulk_flux_tile
      FUNCTION bulk_psi (iflag, ZoL)
!
!=======================================================================
!                                                                      !
!  This function evaluates the stability function, PSI, for wind       !
!  speed (iflag=1) or for air temperature and moisture (iflag=2)       !
!  profiles as function of the stability parameter, ZoL (z/L).         !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    Liu, W.T., K.B. Katsaros, and J.A. Businger, 1979:  Bulk          !
!        parameterization of the air-sea exchange of heat and          !
!        water vapor including the molecular constraints at            !
!        the interface, J. Atmos. Sci, 36, 1722-1735.                  !
!                                                                      !
!=======================================================================
!
      USE mod_kinds
!
!  Function result
!
      real(r8) :: bulk_psi
!
!  Imported variable declarations.
!
      integer, intent(in) :: iflag
      real(r8), intent(in) :: ZoL
!
!  Local variable declarations.
!
      real(r8), parameter :: r3 = 1.0_r8/3.0_r8

      real(r8) :: Fw, chic, chik, cff, pi, psic, psik 
!
!-----------------------------------------------------------------------
!  Compute stability function, PSI.
!-----------------------------------------------------------------------
!
!  Initialize for the zero "ZoL" case.
!
      bulk_psi=0.0_r8
!
!  Unstable conditions.
!
      IF (ZoL.lt.0.0_r8) THEN
        pi=4.0_r8*ATAN(1.0_r8)
        chik=(1.0_r8-16.0_r8*ZoL)**0.25_r8
        IF (iflag.eq.1) THEN
          psik=2.0_r8*LOG(0.5_r8*(1.0_r8+chik))+                        &
     &                LOG(0.5_r8*(1.0_r8+chik*chik))-                   &
     &         2.0_r8*ATAN(chik)+                                       &
     &         0.5_r8*pi
        ELSE IF (iflag.eq.2) THEN
          psik=2.0_r8*LOG(0.5_r8*(1.0_r8+chik*chik))
        END IF
!
!  For very unstable conditions, use free-convection (Fairall).
!
        cff=SQRT(3.0_r8)
        chic=(1.0_r8-12.87_r8*ZoL)**r3
        psic=1.5_r8*LOG(r3*(1.0_r8+chic+chic*chic))-                    &
     &       cff*ATAN((1.0_r8+2.0_r8*chic)/cff)+                        &
     &       pi/cff
!
!  Match Kansas and free-convection forms with weighting Fw.
!
        Fw=1.0_r8/(1.0_r8+ZoL*ZoL)
        bulk_psi=Fw*psik+(1.0_r8-Fw)*psic
!
!  Stable conditions.
!
      ELSE IF (ZoL.gt.0.0_r8) THEN
        bulk_psi=-4.7_r8*ZoL
      END IF
      RETURN
      END FUNCTION bulk_psi
#endif
      END MODULE bulk_flux_mod
