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
      PUBLIC  :: bulk_flux

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
# ifdef CLOUDS
     &                     FORCES(ng) % cloud,                          &
# endif
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
# ifdef CLOUDS
     &                           cloud,                                 &
# endif
     &                           rain, lhflx, lrflx, shflx,             &
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
#  ifdef CLOUDS
      real(r8), intent(in) :: cloud(LBi:,LBj:)
#  endif
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
#  ifdef CLOUDS
      real(r8), intent(in) :: cloud(LBi:UBi,LBj:UBj)
#  endif
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

      real(r8), parameter :: eps = 1.0E-20_r8
      real(r8), parameter :: r3 = 1.0_r8/3.0_r8

      real(r8) :: Bf, Cd, Hl, Hlw, Hscale, Hs, Hsr, IER
      real(r8) :: PairM,  RH, Ri, Scff, Taur, Tcff
      real(r8) :: Wspeed, ZQoL, ZToL, ZoQ, ZoT

      real(r8) :: cff, cff1, cff2, diffh, diffw, oL, upvel
      real(r8) :: wet_bulb
# ifdef LONGWAVE
      real(r8) :: e_sat, vap_p
# endif
# ifdef COOL_SKIN
      real(r8) :: Clam, Cwet, Fc, Hcool, Hsb, Hlb, Qbouy, Qcool, lambd
# endif

      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY) :: delQ
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY) :: delQc
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY) :: delT
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY) :: delTc
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY) :: delW
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY) :: Hlv
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY) :: Q
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY) :: Qair
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY) :: Qpsi
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY) :: Qsea
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY) :: Qstar
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY) :: rhoAir
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY) :: rhoSea
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY) :: Rq
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY) :: Rr
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY) :: Rt
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY) :: TairC
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY) :: TairK
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY) :: Tpsi
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY) :: TseaC
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY) :: TseaK
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY) :: Tstar
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY) :: TVstar
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY) :: VisAir
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY) :: Wgus
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY) :: Wmag
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY) :: Wpsi
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY) :: Wstar
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY) :: ZoW
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY) :: ZQ
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY) :: ZT
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY) :: ZW
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY) :: ZWoL

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
          Wmag(i)=SQRT(Uwind(i,j)*Uwind(i,j)+Vwind(i,j)*Vwind(i,j))
          PairM=Pair(i,j)
          TairC(i)=Tair(i,j)
          TairK(i)=TairC(i)+273.16_r8
          TseaC(i)=t(i,j,N(ng),nrhs,itemp)
          TseaK(i)=TseaC(i)+273.16_r8
          rhoSea(i)=rho(i,j,N(ng))+1000.0_r8
          RH=Hair(i,j)
          SRad(i,j)=srflx(i,j)*Hscale
          Tcff=alpha(i,j)
          Scff=beta(i,j)
          ZT(i)=blk_ZT
          ZQ(i)=blk_ZQ
          ZW(i)=blk_ZW
!
!  Initialize.
!
          delTc(i)=0.0_r8
          delQc(i)=0.0_r8
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
# if defined LONGWAVE

!  Use Berliand (1952) formula.
!
          cff=(0.7859_r8+0.03477_r8*TairC(i))/                          &
     &        (1.0_r8+0.00412_r8*TairC(i))
          e_sat=10.0_r8**cff
          vap_p=e_sat*RH
          cff2=TairK(i)*TairK(i)*TairK(i)
          cff1=cff2*TairK(i)
          LRad(i,j)=-emmiss*StefBo*                                     &
     &              (cff1*(0.39_r8-0.5_r8*SQRT(0.01_r8*vap_p))*         &
     &                    (1.0_r8-0.6823_r8*cloud(i,j)*cloud(i,j))+     &
     &               cff2*4.0_r8*(TseaK(i)-TairK(i)))

# elif defined LONGWAVE_OUT

!  Treat input longwave data as downwelling radiation only and add
!  outgoing IR from model sea surface temperature.
!
          LRad(i,j)=lrflx(i,j)*Hscale-                                  &
     &              emmiss*StefBo*TseaK(i)*TseaK(i)*TseaK(i)*TseaK(i)
 
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
          Qair(i)=(1.0007_r8+3.46E-6_r8*PairM)*6.1121_r8*               &
     &            EXP(17.502_r8*TairC(i)/(240.97_r8+TairC(i)))
!
!  Compute specific humdity, Q (kg/kg).  Check whether input relative
!  humidity is in kg/kg or g/kg.
!
          IF (RH.lt.2.0_r8) THEN
            cff=Qair(i)*RH
            Q(i)=0.62197_r8*(cff/(PairM-0.378_r8*cff))
          ELSE
            Q(i)=RH/1000.0_r8
          END IF
          Qair(i)=0.62197_r8*(Qair(i)/(PairM-0.378_r8*Qair(i)))
!
!  Compute water saturation vapor pressure (mb), using Teten formula.
!
          Qsea(i)=(1.0007_r8+3.46E-6_r8*PairM)*6.1121_r8*               &
     &            EXP(17.502_r8*TseaC(i)/(240.97_r8+TseaC(i)))
!
!  Qsea reduced for salinity (Kraus, 1972, pp 46).
!
          Qsea(i)=Qsea(i)*0.98_r8
!
!  Convert Qsea from mb to specific humidity (kg/kg).
!
          Qsea(i)=0.62197_r8*(Qsea(i)/(PairM-0.378_r8*Qsea(i)))
!
!-----------------------------------------------------------------------
!  Compute Monin-Obukhov similarity parameters for wind (Wstar),
!  heat (Tstar), and moisture (Qstar), Liu et al. (1979).
!-----------------------------------------------------------------------
!
!  Moist air density (kg/m3).
!
          rhoAir(i)=PairM*100.0_r8/(blk_Rgas*TairK(i)*                  &
     &                              (1.0_r8+0.61_r8*Q(i)))
!
!  Kinematic viscosity of dry air (m2/s), Andreas (1989).
!
          VisAir(i)=1.326E-5_r8*                                        &
     &              (1.0_r8+TairC(i)*(6.542E-3_r8+TairC(i)*             &
     &               (8.301E-6_r8-4.84E-9_r8*TairC(i))))
!
!  Compute latent heat of vaporization (J/kg) at sea surface, Hlv.
!
          Hlv(i)=(2.501_r8-0.00237_r8*TseaC(i))*1.0E+6_r8
!
!  Assume that wind is measured relative to sea surface and include
!  gustiness.
!
          Wgus(i)=0.5_r8
          delW(i)=SQRT(Wmag(i)*Wmag(i)+Wgus(i)*Wgus(i))
          delQ(i)=Q(i)-Qsea(i)
          delT(i)=TairC(i)-TseaC(i)
!
!  Initial guesses for Monon-Obukhov similarity scales.
!
          ZWoL(i)=0.0_r8
          ZoW(i)=0.0005_r8
          Wstar(i)=0.04_r8*delW(i)
          Tstar(i)=0.04_r8*delT(i)
          Qstar(i)=0.04_r8*delQ(i)
          TVstar(i)=Tstar(i)*(1.0_r8+0.61_r8*Q(i))+                     &
     &              0.61_r8*TairK(i)*Qstar(i)
!
!  Compute Richardson number.
!
          Ri=g*ZW(i)*(delT(i)+0.61_r8*TairK(i)*delQ(i))/                &
     &       (TairK(i)*delW(i)*delW(i))
          Rt(i)=0.0_r8
          Rq(i)=0.0_r8
        END DO
!
!  Iterate until convergence. It usually converges within four
!  iterations.
!
        DO Iter=1,IterMax
          DO i=MIN(Istr-1,IstrR),IendR
!
!  Compute Monin-Obukhov stability parameter, Z/L.
!
            oL=g*vonKar*TVstar(i)/                                      &
     &         (TairK(i)*(1.0_r8+0.61_r8*Q(i))*Wstar(i)*Wstar(i)+eps)
            ZWoL(i)=ZW(i)*oL
            ZToL=ZT(i)*oL
            ZQoL=ZQ(i)*oL
!
!  Evaluate stability functions at Z/L.
!
            Wpsi(i)=bulk_psi(1,ZWoL(i))
            Tpsi(i)=bulk_psi(2,ZToL)
            Qpsi(i)=bulk_psi(2,ZQoL)
!
!  Compute wind scaling parameters, Wstar.
!
            ZoW(i)=0.011_r8*Wstar(i)*Wstar(i)/g+                        &
     &             0.11_r8*VisAir(i)/(Wstar(i)+eps)
            Wstar(i)=MAX(eps,delW(i)*vonKar/(LOG(ZW(i)/ZoW(i))-Wpsi(i)))
!
!  Computes roughness Reynolds number for wind (Rr), heat (Rt),
!  and moisture (Rq). Use Liu et al. (1976) look-up table to
!  compute "Rt" and "Rq" as function of "Rr".
!
            Rr(i)=ZoW(i)*Wstar(i)/VisAir(i)
            Rr(i)=MAX(MIN(1000.0_r8,Rr(i)),0.0_r8)
          END DO
          DO k=1,8
            DO i=MIN(Istr-1,IstrR),IendR
              IF ((Liu_Rr(k).le.Rr(i)).and.(Rr(i).le.Liu_Rr(k+1))) THEN
                Rt(i)=Liu_a(k,1)*Rr(i)**Liu_b(k,1)
                Rq(i)=Liu_a(k,2)*Rr(i)**Liu_b(k,2)
              END IF
            END DO
          END DO
!
!  Compute heat and moisture scaling parameters, Tstar and Qstar.
!
          DO i=MIN(Istr-1,IstrR),IendR
            cff=VisAir(i)/(Wstar(i)+eps)
            ZoT=Rt(i)*cff
            ZoQ=Rq(i)*cff
            cff=vonKar*blk_fdg
            Tstar(i)=(delT(i)+delTc(i))*                                &
     &               cff/(LOG(ZT(i)/ZoT)-Tpsi(i))
            Qstar(i)=(delQ(i)+delQc(i))*                                &
     &               cff/(LOG(ZQ(i)/ZoQ)-Qpsi(i))
!
!  Compute gustiness in wind speed.
!
            TVstar(i)=Tstar(i)*(1.0_r8+0.61_r8*Q(i))+                   &
     &                0.61_r8*TairK(i)*Qstar(i)
            Bf=-g/TairK(i)*Wstar(i)*TVstar(i)
            IF (Bf.gt.0) THEN
              Wgus(i)=blk_beta*(Bf*blk_Zabl)**r3
            ELSE
              Wgus(i)=0.0_r8
            END IF
            delW(i)=SQRT(Wmag(i)*Wmag(i)+Wgus(i)*Wgus(i))
# ifdef COOL_SKIN
!
!-----------------------------------------------------------------------
!  Cool Skin correction.
!-----------------------------------------------------------------------
!
!  Cool skin correction constants. Clam: part of Saunders constant
!  lambda; Cwet: slope of saturation vapor.
!
            Clam=16.0_r8*g*blk_Cpw*(rhoSea(i)*blk_visw)**3/             &
     &           (blk_tcw*blk_tcw*rhoAir(i)*rhoAir(i))
            Cwet=0.622_r8*Hlv(i)*Qsea(i)/                               &
     &           (blk_Rgas*TseaK(i)*TseaK(i))
!
!  Set initial guesses for cool-skin layer thickness (Hcool).
!
            Hcool=0.001_r8
!
!  Backgound sensible and latent heat.
!
            Hsb=-rhoAir(i)*blk_Cpa*Wstar(i)*Tstar(i)
            Hlb=-rhoAir(i)*Hlv(i)*Wstar(i)*Qstar(i)
!
!  Mean absoption in cool-skin layer.
!
            Fc=0.137_r8+11.0_r8*Hcool-                                  &
     &         (1.0_r8-EXP(-Hcool/8.0e-4))*6.6E-5_r8/Hcool
!
!  Total cooling at the interface.
!
            Qcool=LRad(i,j)+Hsb+Hlb-SRad(i,j)*Fc
            Qbouy=Tcff*Qcool+Scff*Hlb*blk_Cpw/Hlv(i)
!
!  Compute temperature and moisture change.
!
            IF ((Qcool.gt.0.0_r8).and.(Qbouy.gt.0.0_r8)) THEN
              lambd=6.0_r8/                                             &
     &              (1.0_r8+(Clam*Qbouy/(Wstar(i)+eps)**4)**0.75_r8)**r3
              Hcool=lambd*blk_visw/(SQRT(rhoAir(i)/rhoSea(i))*          &
     &                              Wstar(i)+eps)
              delTc(i)=Qcool*Hcool/blk_tcw
            ELSE
              delTc(i)=0.0_r8
            END IF
            delQc(i)=Cwet*delTc(i)
# endif
          END DO
        END DO
!
!-----------------------------------------------------------------------
!  Compute Atmosphere/Ocean fluxes.
!-----------------------------------------------------------------------
!
        DO i=MIN(Istr-1,IstrR),IendR
!
!  Compute transfer coefficients for momentun (Cd).
!
          Wspeed=SQRT(Wmag(i)*Wmag(i)+Wgus(i)*Wgus(i))
          Cd=Wstar(i)*Wstar(i)/(Wspeed*Wspeed+eps)
!
!  Compute turbulent sensible heat flux (W/m2), Hs.
!
            Hs=-blk_Cpa*rhoAir(i)*Wstar(i)*Tstar(i)
!
!  Compute sensible heat flux (W/m2) due to rainfall (kg/m2/s), Hsr.
!
            diffw=2.11E-5_r8*(TairK(i)/273.16_r8)**1.94_r8
            diffh=0.02411_r8*(1.0_r8+TairC(i)*                          &
     &                        (3.309E-3_r8-1.44E-6_r8*TairC(i)))/       &
     &            (rhoAir(i)*blk_Cpa)
            cff=Qair(i)*Hlv(i)/(blk_Rgas*TairK(i)*TairK(i))
            wet_bulb=1.0_r8/(1.0_r8+0.622_r8*(cff*Hlv(i)*diffw)/        &
     &                                       (blk_Cpa*diffh))
            Hsr=rain(i,j)*wet_bulb*blk_Cpw*                             &
     &          ((TseaC(i)-TairC(i))+(Qsea(i)-Q(i))*Hlv(i)/blk_Cpa)
            SHeat(i,j)=(Hs+Hsr)
# ifdef MASKING
            SHeat(i,j)=SHeat(i,j)*rmask(i,j)
# endif
!
!  Compute turbulent latent heat flux (W/m2), Hl.
!
            Hl=-Hlv(i)*rhoAir(i)*Wstar(i)*Qstar(i)
!
!  Compute Webb correction (Webb effect) to latent heat flux, Hlw.
!
            upvel=-1.61_r8*Wstar(i)*Qstar(i)-                         &
     &            (1.0_r8+1.61_r8*Q(i))*Wstar(i)*Tstar(i)/TairK(i)
            Hlw=rhoAir(i)*Hlv(i)*upvel*Q(i)
            LHeat(i,j)=(Hl+Hlw)
# ifdef MASKING
            LHeat(i,j)=LHeat(i,j)*rmask(i,j)
# endif
!
!  Compute momentum flux (N/m2) due to rainfall (kg/m2/s).
!
            Taur=0.85_r8*rain(i,j)*Wmag(i)
!
!  Compute wind stress components (N/m2), Tau.
!
            cff=rhoAir(i)*Cd*Wspeed
            Taux(i,j)=(cff*Uwind(i,j)+Taur*SIGN(1.0_r8,Uwind(i,j)))
# ifdef MASKING
            Taux(i,j)=Taux(i,j)*rmask(i,j)
# endif
            Tauy(i,j)=(cff*Vwind(i,j)+Taur*SIGN(1.0_r8,Vwind(i,j)))
# ifdef MASKING
            Tauy(i,j)=Tauy(i,j)*rmask(i,j)
# endif
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
