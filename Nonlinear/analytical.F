#include "cppdefs.h"
      MODULE analytical_mod

      implicit none

      CONTAINS
!!
!!======================================================================
!! Copyright (c) 2002 ROMS/TOMS Group                                  !
!!================================================= Hernan G. Arango ===
!!                                                                     !
!! ANALYTICAL PACKAGE:                                                 !
!!                                                                     !
!! This package is used to provide various analytical fields to the    !
!! model when appropriate.                                             !
!!                                                                     !
!! Routines:                                                           !
!!                                                                     !
!! ana_biology         Analytical initial conditions for biological    !
!!                       tracers.                                      !
!! ana_btflux          Analytical kinematic bottom flux of tracer      !
!!                       type variables.                               !
!! ana_cloud           Analytical cloud fraction.                      !
!! ana_diag            Customized diagnostics.                         !
!! ana_fsobc           Analytical free-surface boundary conditions.    !
!! ana_grid            Analytical model grid set-up.                   !
!! ana_humid           Analytical surface air humidity.                !
!! ana_initial         Analytical initial conditions for momentum,     !
!!                       free surface and tracers.                     !
!! ana_m2clima         Analytical 2D momentum climatology.             !
!! ana_m2obc           Analytical 2D momentum boundary conditions.     !
!! ana_m3clima         Analytical 3D momentum climatology.             !
!! ana_m3obc           Analytical 3D momentum boundary conditions.     !
!! ana_mask            Analytical Land/Sea masking.                    !
!! ana_pair            Analytical surface air pressure.                !
!! ana_psource         Analytical mass/tracer point sources/sinks.     !
!! ana_rain            Analytical rain fall rate.                      !
!! ana_sediment        Analytical initial conditions for sediment      !
!!                       tracers and bottom/bed properties.            !
!! ana_smflux          Analytical kinematic surface momentum flux      !
!!                       (wind stress).                                !
!! ana_spinning        Analytical time variable rotation forces.       !
!! ana_srflux          Analytical kinematic surface shortwave          !
!!                       radiation flux.                               !
!! ana_specir          Analytical calculation of spectral downwelling  !
!!                       irradiance.                                   !
!! ana_ssh             Analytical sea surface height climatology.      !
!! ana_sst             Analytical SST and dQdSST which are used        !
!!                       for heat flux correction.                     !
!! ana_sss             Analytical sea surface salinity.                !
!! ana_stflux          Analytical kinematic surface flux of tracer     !
!!                       type variables.                               !
!! ana_tair            Analytical surface air temperature.             !
!! ana_tclima          Analytical tracer climatology fields.           !
!! ana_tobc            Analytical tracer boundary conditions.          !
!! ana_vmix            Analytical vertical mixing coefficients for     !
!!                       momentum and tracers.                         !
!! ana_winds           Analytical surface winds.                       !
!! ana_wwave           Analytical wind induced wave amplitude,         !
!!                       direction and period.                         !
!!                                                                     !
!!======================================================================
!!
#if defined ANA_BIOLOGY && defined BIOLOGY
      SUBROUTINE ana_biology (ng, tile)
!
!=======================================================================
!                                                                      !
!  This routine sets initial conditions for biological tracer fields   !
!  using analytical expressions.                                       !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_ocean
!    
      integer, intent(in) :: ng, tile

# include "tile.h"
!
      CALL ana_biology_tile (ng, Istr, Iend, Jstr, Jend,                &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       OCEAN(ng) % t)
      RETURN
      END SUBROUTINE ana_biology
!
!***********************************************************************
      SUBROUTINE ana_biology_tile (ng, Istr, Iend, Jstr, Jend,          &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             t)
!***********************************************************************
!
      USE mod_param
      USE mod_biology
      USE mod_scalars
!
# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      USE exchange_3d_mod, ONLY : exchange_r3d_tile
!
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj

# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
# else
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
# endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, is, itrc, j, k

#ifdef BIO_FASHAM
      real(r8) :: SiO4, cff1, cff2, temp
#elif defined ECOSIM
      real(r8) :: cff1, cff2, cff3, cff4, cff5, cff6, cff7, cff8, cff9
      real(r8) :: cff10, cff11, cff12, cff13, cff14, cff15
      real(r8) :: salt, sftm, temp
#endif

# include "set_bounds.h"

# ifdef BIO_FASHAM
!-----------------------------------------------------------------------
!  Fasham type, nitrogen-based biology model.
!-----------------------------------------------------------------------
!
      cff1=20.0_r8/3.0_r8
      cff2= 2.0_r8/3.0_r8
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            temp=t(i,j,k,1,itemp)
            IF (temp.lt.8.0_r8) THEN
              SiO4=30.0_r8
            ELSE IF ((temp.ge.8.0_r8).and.(temp.le.11.0_r8)) THEN
              SiO4=30.0_r8-((temp-8.0_r8)*cff1)
            ELSE IF ((temp.gt.11.0_r8).and.(temp.le.13.0_r8)) THEN
              SiO4=10.0_r8-((temp-11.0_r8)*4.0_r8)
            ELSE IF ((temp.gt.13.0_r8).and.(temp.le.16.0_r8)) THEN
              SiO4=2.0_r8-((temp-13.0_r8)*cff2)
            ELSE IF (temp.gt.16.0_r8) THEN
              SiO4=0.0_r8
            END IF
            t(i,j,k,1,iNO3_)=1.67_r8+0.5873_r8*SiO4+                    &
     &                               0.0144_r8*SiO4**2+                 &
     &                               0.0003099_r8*SiO4**3
            t(i,j,k,1,iPhyt)=0.08_r8
            t(i,j,k,1,iZoop)=0.06_r8
            t(i,j,k,1,iNH4_)=0.1_r8
            t(i,j,k,1,iLDeN)=0.02_r8
            t(i,j,k,1,iSDeN)=0.04_r8
            t(i,j,k,1,iChlo)=0.02_r8
# ifdef CARBON
            t(i,j,k,1,iTIC_)=2100.0_r8
            t(i,j,k,1,iTAlk)=2350.0_r8
            t(i,j,k,1,iLDeC)=0.002_r8
            t(i,j,k,1,iSDeC)=0.06_r8
# endif
            t(i,j,k,2,iNO3_)=t(i,j,k,1,iNO3_)
            t(i,j,k,2,iPhyt)=t(i,j,k,1,iPhyt)
            t(i,j,k,2,iNH4_)=t(i,j,k,1,iNH4_)
            t(i,j,k,2,iLDeN)=t(i,j,k,1,iLDeN)
            t(i,j,k,2,iSDeN)=t(i,j,k,1,iSDeN)
            t(i,j,k,2,iChlo)=t(i,j,k,1,iChlo)
# ifdef CARBON
            t(i,j,k,2,iTIC_)=t(i,j,k,1,iTIC_)
            t(i,j,k,2,iTAlk)=t(i,j,k,1,iTAlk)
            t(i,j,k,2,iLDeC)=t(i,j,k,1,iLDeC)
            t(i,j,k,2,iSDeC)=t(i,j,k,1,iSDeC)
# endif
          END DO
        END DO
      END DO
# elif defined ECOSIM
!
!---------------------------------------------------------------------
!  EcoSim initial fields.
!---------------------------------------------------------------------
!
! Assumed maximum temperature gradient.
!
      cff3=1.0_r8/14.0_r8
      cff4=1.0_r8/16.0_r8
      cff5=32.0_r8
      cff7=1.0_r8/0.0157_r8
      cff8=1.0_r8/6.625_r8
      cff9=1.0_r8/16.0_r8
      cff10=1.0_r8/15.0_r8
      cff11=1.0_r8/8.0_r8
      cff12=1.0_r8/128.0_r8
      cff13=1.0_r8/1000.0_r8
      cff14=1.0_r8/12.0_r8
      cff15=cff5*cff8*cff14                  ! mole N : gram Chl

      DO k=N(ng),1,-1
        DO j=JstrR,JendR
          DO i=IstrR,IendR
!
! Initialization of surface chlorophyll.
!
            sftm=t(i,j,N(ng),1,itemp)
            temp=t(i,j,k,1,itemp)
            salt=t(i,j,k,1,isalt)
            cff1=-0.0827_r8*sftm+2.6386_r8
            cff2=MAX(0.00001_r8,cff1*(1.0_r8-(sftm-temp)*cff3))
!
! Initialization of nutrients.
!
            t(i,j,k,1,iNH4_)=0.053_r8*temp+0.7990_r8
            t(i,j,k,1,iNO3_)=8.5_r8-cff2*cff15-t(i,j,k,1,iNH4_)
            t(i,j,k,1,iPO4_)=(t(i,j,k,1,iNH4_)+t(i,j,k,1,iNO3_))*cff4
            t(i,j,k,1,iFeO_)=1.0_r8
!
! Assuming diatoms are 75% of initialized chlorophyll.
!
            t(i,j,k,1,iSiO_)=5.5_r8-(cff2*0.75_r8)*cff15*1.20_r8
            t(i,j,k,1,iDIC_)=2000.0_r8
!
! Bacteria Initialization.
!
            DO is=1,Nbac
              t(i,j,k,1,iBacC(is))=0.85_r8
              t(i,j,k,1,iBacN(is))=t(i,j,k,1,iBacC(is))*N2cBAC(ng)
              t(i,j,k,1,iBacP(is))=t(i,j,k,1,iBacC(is))*P2cBAC(ng)
              t(i,j,k,1,iBacF(is))=t(i,j,k,1,iBacC(is))*Fe2cBAC(ng)
            END DO
!
! Initialize phytoplankton populations.
!
            t(i,j,k,1,iPhyC(1))=MAX(0.02_r8,                            &
     &                              0.75_r8*0.75_r8*cff5*cff2*cff14)
            t(i,j,k,1,iPhyC(2))=MAX(0.02_r8,                            &
     &                              0.75_r8*0.25_r8*cff5*cff2*cff14)
            t(i,j,k,1,iPhyC(3))=MAX(0.02_r8,                            &
     &                              0.125_r8*cff5*cff2*cff14)
            t(i,j,k,1,iPhyC(4))=t(i,j,k,1,iPhyC(3))
            DO is=1,Nphy
              t(i,j,k,1,iPhyN(is))=t(i,j,k,1,iPhyC(is))*cff8
              t(i,j,k,1,iPhyP(is))=t(i,j,k,1,iPhyN(is))*cff4
              t(i,j,k,1,iPhyF(is))=t(i,j,k,1,iPhyC(is))*cff13
              IF (iPhyS(is).gt.0) THEN
                t(i,j,k,1,iPhyS(is))=t(i,j,k,1,iPhyN(is))*1.20_r8
              END IF
!
!  Initialize Pigments in ugrams/liter (not umole/liter).
!  Chlorophyll-a
!
              cff6=12.0_r8/cff5
              t(i,j,k,1,iPigs(is,1))=cff6*t(i,j,k,1,iPhyC(is))
!
!  Chlorophyll-b.
!
              cff6=cff5-b_C2Cl(is,ng)
              IF (iPigs(is,2).gt.0) THEN
                 t(i,j,k,1,iPigs(is,2))=t(i,j,k,1,iPigs(is,1))*         &
     &                                  (b_ChlB(is,ng)+                 &
     &                                   mxChlB(is,ng)*cff6)
              END IF
!
!  Chlorophyll-c.
!
              IF (iPigs(is,3).gt.0) THEN
                 t(i,j,k,1,iPigs(is,3))=t(i,j,k,1,iPigs(is,1))*         &
     &                                  (b_ChlC(is,ng)+                 &
     &                                   mxChlC(is,ng)*cff6) 
              END IF
!
!  Photosynthetic Carotenoids.
!
              IF (iPigs(is,4).gt.0) THEN
                 t(i,j,k,1,iPigs(is,4))=t(i,j,k,1,iPigs(is,1))*         &
     &                                  (b_PSC(is,ng)+                  &
     &                                   mxPSC(is,ng)*cff6)
              END IF
!
!  Photoprotective Carotenoids.
!
              IF (iPigs(is,5).gt.0) THEN
                 t(i,j,k,1,iPigs(is,5))=t(i,j,k,1,iPigs(is,1))*         &
     &                                  (b_PPC(is,ng)+                  &
     &                                   mxPPC(is,ng)*cff6) 
              END IF
!
!  Low Urobilin Phycoeurythin Carotenoids.
!
              IF (iPigs(is,6).gt.0) THEN
                 t(i,j,k,1,iPigs(is,6))=t(i,j,k,1,iPigs(is,1))*         &
     &                                  (b_LPUb(is,ng)+                 &
     &                                   mxLPUb(is,ng)*cff6) 
              END IF
!
!  High Urobilin Phycoeurythin Carotenoids.
!
              IF (iPigs(is,7).gt.0) THEN
                 t(i,j,k,1,iPigs(is,7))=t(i,j,k,1,iPigs(is,1))*         &
     &                                  (b_HPUb(is,ng)+                 &
     &                                   mxHPUb(is,ng)*cff6) 
              END IF
            END DO
!
! DOC initialization.
!
            cff6=MAX(0.001_r8,-0.9833_r8*salt+33.411_r8)
            t(i,j,k,1,iDOMC(1))=0.1_r8
            t(i,j,k,1,iDOMN(1))=t(i,j,k,1,iDOMC(1))*cff8
            t(i,j,k,1,iDOMP(1))=t(i,j,k,1,iDOMN(1))*cff9
            t(i,j,k,1,iCDMC(1))=t(i,j,k,1,iDOMC(1))*cDOCfrac_c(1,ng)
            t(i,j,k,1,iDOMC(2))=15.254_r8*cff6+70.0_r8
            t(i,j,k,1,iDOMN(2))=t(i,j,k,1,iDOMC(2))*cff10
            t(i,j,k,1,iDOMP(2))=0.0_r8
            t(i,j,k,1,iCDMC(2))=(0.243_r8*cff6+0.055_r8)*cff7
!
! Fecal Initialization.
!
            DO is=1,Nfec
              t(i,j,k,1,iFecC(is))=0.002_r8
              t(i,j,k,1,iFecN(is))=t(i,j,k,1,iFecC(is))*cff11
              t(i,j,k,1,iFecP(is))=t(i,j,k,1,iFecC(is))*cff12
              t(i,j,k,1,iFecF(is))=t(i,j,k,1,iFecC(is))*cff13
              t(i,j,k,1,iFecS(is))=t(i,j,k,1,iFecC(is))*cff11
            END DO
          END DO
        END DO
      END DO
# endif
# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      DO itrc=NAT+1,NAT+NBT
        CALL exchange_r3d_tile (ng, Istr, Iend, Jstr, Jend,             &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          t(:,:,:,1,itrc))
        CALL exchange_r3d_tile (ng, Istr, Iend, Jstr, Jend,             &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          t(:,:,:,2,itrc))
      END DO
# endif
      RETURN
      END SUBROUTINE ana_biology_tile
#endif

#ifdef SOLVE3D
# if defined ANA_BTFLUX || defined ANA_BSFLUX || defined ANA_BPFLUX
      SUBROUTINE ana_btflux (ng, tile, itrc)
!
!=======================================================================
!                                                                      !
!  This routine sets kinematic bottom flux of tracer type variables    !
!  (tracer units m/s).                                                 !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_forces
!
      integer, intent(in) :: ng, tile, itrc

# include "tile.h"
!
      CALL ana_btflux_tile (ng, Istr, Iend, Jstr, Jend, itrc,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      FORCES(ng) % btflx)
      RETURN
      END SUBROUTINE ana_btflux
!
!***********************************************************************
      SUBROUTINE ana_btflux_tile (ng, Istr, Iend, Jstr, Jend, itrc,     &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            btflx)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer,  intent(in) :: ng, Iend, Istr, Jend, Jstr, itrc
      integer, intent(in) :: LBi, UBi, LBj, UBj

# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: btflx(LBi:,LBj:,:)
# else
      real(r8), intent(inout) :: btflx(LBi:UBi,LBj:UBj,NT(ng))
# endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set kinematic bottom heat flux (degC m/s) at horizontal RHO-points.
!-----------------------------------------------------------------------
!
      IF (itrc.eq.itemp) THEN
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            btflx(i,j,itrc)=0.0_r8
          END DO
        END DO
!
!-----------------------------------------------------------------------
!  Set kinematic bottom salt flux (m/s) at horizontal RHO-points,
!  scaling by bottom salinity is done elsewhere.
!-----------------------------------------------------------------------
!
      ELSE IF (itrc.eq.isalt) THEN
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            btflx(i,j,itrc)=0.0_r8
          END DO
        END DO
!
!-----------------------------------------------------------------------
!  Set kinematic bottom flux (T m/s) of passive tracers, if any.
!-----------------------------------------------------------------------
!
      ELSE
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            btflx(i,j,itrc)=0.0_r8
          END DO
        END DO
      END IF
      RETURN
      END SUBROUTINE ana_btflux_tile
# endif
#endif

#if defined ANA_CLOUD && defined CLOUDS
      SUBROUTINE ana_cloud (ng, tile)
!
!=======================================================================
!                                                                      !
!  This routine sets cloud fraction using an analytical expression.    !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_forces
!
      integer, intent(in) :: ng, tile

# include "tile.h"
!
      CALL ana_cloud_tile (ng, Istr, Iend, Jstr, Jend,                  &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng) % cloud)
      RETURN
      END SUBROUTINE ana_cloud
!
!***********************************************************************
      SUBROUTINE ana_cloud_tile (ng, Istr, Iend, Jstr, Jend,            &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           cloud)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars

# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
!
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj

# ifdef ASSUMED_SHAPE
      real(r8), intent(out) :: cloud(LBi:,LBj:)
# else
      real(r8), intent(out) :: cloud(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: iday, i, j, month, year
      real(r8) :: Cval, hour, yday

# ifdef PAPA_CLM
      real(r8), dimension(14) :: Coktas =                               &
     &         (/ 6.29_r8, 6.26_r8, 6.31_r8, 6.31_r8, 6.32_r8,          &
     &            6.70_r8, 7.12_r8, 7.26_r8, 6.93_r8, 6.25_r8,          &
     &            6.19_r8, 6.23_r8, 6.31_r8, 6.29_r8          /)

      real(r8), dimension(14) :: Cyday =                                &
     &          (/  0.0_r8,  16.0_r8,  46.0_r8,  75.0_r8, 105.0_r8,     &
     &            136.0_r8, 166.0_r8, 197.0_r8, 228.0_r8, 258.0_r8,     &
     &            289.0_r8, 319.0_r8, 350.0_r8, 365.0_r8           /)
# endif

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set analytical cloud fraction (%/100): 0=clear sky, 1:overcast sky.
!-----------------------------------------------------------------------

# if defined PAPA_CLM
!
!  OWS Papa cloud climatology.
!
      CALL caldate (r_date, tdays(ng), year, yday, month, iday, hour)
      DO i=1,13
        IF ((yday.ge.Cyday(i)).and.(yday.le.Cyday(i+1))) THEN
          Cval=0.125_r8*(Coktas(i  )*(Cyday(i+1)-yday)+                 &
     &                   Coktas(i+1)*(yday-Cyday(i)))/                  &
     &                  (Cyday(i+1)-Cyday(i))
        END IF
      END DO
# elif defined ADRIATIC1 || defined ADRIATIC2
      Cval=0.4_r8
# elif defined BENCHMARK1 || defined BENCHMARK2 || defined BENCHMARK3
      Cval=0.6_r8
# elif defined NJ_BIGHT
      Cval=0.3_r8
# else
      Cval=0.0_r8
# endif
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          cloud(i,j)=Cval
        END DO
      END DO
# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, cloud)
# endif
      RETURN
      END SUBROUTINE ana_cloud_tile
#endif

#ifdef ANA_DIAG
      SUBROUTINE ana_diag (ng, tile)
!
!=======================================================================
!                                                                      !
!  This routine is provided so the USER can compute any specialized    !
!  diagnostics.  If activated, this routine is call at end of every    !
!  3D-equations timestep.                                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_ocean
!
      integer, intent(in) :: ng, tile

# include "tile.h"
!
      CALL ana_diag_tile (ng, Istr, Iend, Jstr, Jend,                   &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    OCEAN(ng) % ubar,                             &
     &                    OCEAN(ng) % vbar,                             &
     &                    OCEAN(ng) % u,                                &
     &                    OCEAN(ng) % v)
      RETURN
      END SUBROUTINE ana_diag
!
!***********************************************************************
      SUBROUTINE ana_diag_tile (ng, Istr, Iend, Jstr, Jend,             &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          ubar, vbar, u, v)
!***********************************************************************
!
      USE mod_param
      USE mod_iounits
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj

# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: vbar(LBi:,LBj:,:)
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
# else
      real(r8), intent(in) :: ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),2)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k
      real(r8) :: umax, ubarmax, vmax, vbarmax

# ifdef SEAMOUNT
!
!  Open USER file.
!
      IF (iic.eq.ntstart) THEN
        OPEN (usrout,file=usrname,form='formatted',status='unknown',    &
     &        err=40)
        GO TO 60
  40    WRITE (stdout,50) usrname
  50    FORMAT (' ANA_DIAG - unable to open output file: ',a)
        exit_flag=1
  60    CONTINUE
      END IF
!
!  Write out maximum values of velocity.
!
      umax=0.0_r8
      vmax=0.0_r8
      ubarmax=0.0_r8
      vbarmax=0.0_r8
      DO k=1,N(ng)
        DO j=0,Mm(ng)+1
          DO i=1,Lm(ng)+1
            umax=MAX(umax,u(i,j,k,nnew))
          END DO
        END DO
        DO j=1,Mm(ng)+1
          DO i=0,Lm(ng)+1
            vmax=MAX(vmax,v(i,j,k,nnew))
          END DO
        END DO
      END DO
      DO j=0,Mm(ng)+1
        DO i=1,Lm(ng)+1
          ubarmax=MAX(ubarmax,ubar(i,j,knew))
        END DO
      END DO
      DO j=1,Mm(ng)+1
        DO i=0,Lm(ng)+1
          vbarmax=MAX(vbarmax,vbar(i,j,knew))
        END DO
      END DO
!
!  Write out maximum values on velocity.
!
      WRITE (usrout,70) tdays(ng), ubarmax, vbarmax, umax, vmax
  70  FORMAT (2x,f13.6,2x,1pe13.6,2x,1pe13.6,2x,1pe13.6,2x,1pe13.6)
# endif
      RETURN
      END SUBROUTINE ana_diag_tile
#endif

#ifdef ANA_FSOBC
      SUBROUTINE ana_fsobc (ng, tile)
!
!=======================================================================
!                                                                      !
!  This routine sets free-surface open boundary conditions using       !
!  analytical expressions.                                             !
!                                                                      !
!=======================================================================
!
      USE mod_param
!
      integer, intent(in) :: ng, tile

# include "tile.h"
!
      CALL ana_fsobc_tile (ng, Istr, Iend, Jstr, Jend,                  &
     &                     LBi, UBi, LBj, UBj)
      RETURN
      END SUBROUTINE ana_fsobc
!
!***********************************************************************
      SUBROUTINE ana_fsobc_tile (ng, Istr, Iend, Jstr, Jend,            &
     &                           LBi, UBi, LBj, UBj)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
      USE mod_grid
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j
      real(r8) :: cff, fac, omega, phase, val

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Free-surface open boundary conditions.
!-----------------------------------------------------------------------
!
# if defined KELVIN
      fac=1.0_r8                                ! zeta0
      omega=2.0_r8*pi/(12.42_r8*3600.0_r8)      ! M2 Tide period
      IF (WESTERN_EDGE) THEN
        DO j=JstrR,JendR
          val=fac*EXP(-GRID(ng)%f(Istr-1,j)*GRID(ng)%yp(Istr-1,j)/      &
     &                SQRT(g*GRID(ng)%h(Istr-1,j)))
          BOUNDARY(ng)%zeta_west(j)=val*COS(omega*time(ng))
        END DO
      END IF
      IF (EASTERN_EDGE) THEN
        DO j=JstrR,JendR
          cff=1.0_r8/SQRT(g*GRID(ng)%h(Istr-1,j))
     &    val=fac*EXP(-GRID(ng)%f(Istr-1,j)*GRID(ng)%yp(Iend,j)*cff
          BOUNDARY(ng)%zeta_east(j)=val*COS(omega*GRID(ng)%xp(Iend,j)*  &
     &                                      cff-omega*time(ng))
        END DO
      END IF
# elif defined SED_TEST1
      IF (WESTERN_EDGE) THEN
        fac=100.0_r8
        DO j=JstrR,JendR
          BOUNDARY(ng)%zeta_west(j)=9.0E-06_r8*fac
        END DO
      END IF
      IF (EASTERN_EDGE) THEN
        fac=100.0_r8
        DO j=JstrR,JendR
          BOUNDARY(ng)%zeta_east(j)=9.0E-06_r8*REAL(Iend+1,r8)*fac
        END DO
      END IF
# elif defined TEST_CHAN
      IF (WESTERN_EDGE) THEN
        DO j=JstrR,JendR
          BOUNDARY(ng)%zeta_west(j)=-0.396_r8
        END DO
      END IF
# elif defined WEDDELL
      IF (WESTERN_EDGE) THEN
        fac=TANH((tdays(ng)-dstart)/1.0_r8)
        omega=2.0_r8*pi*time(ng)/(12.42_r8*3600.0_r8)  !  M2 Tide period
        val=0.53_r8+(0.53_r8-0.48_r8)/REAL(Iend+1,r8)
        phase=(277.0_r8+(277.0_r8-240.0_r8)/REAL(Iend+1),r8)*deg2rad
        DO j=JstrR,JendR
          BOUNDARY(ng)%zeta_west(j)=fac*val*COS(omega-phase)
        END DO
      END IF
      IF (EASTERN_EDGE) THEN
        fac=TANH((tdays(ng)-dstart)/1.0_r8)
        omega=2.0_r8*pi*time(ng)/(12.42_r8*3600.0_r8)  !  M2 Tide period
        val=0.53_r8+(0.53_r8-0.48_r8)
        phase=(277.0_r8+(277.0_r8-240.0_r8))*deg2rad
        DO j=JstrR,JendR
          BOUNDARY(ng)%zeta_east(j)=fac*val*COS(omega-phase)
        END DO
      END IF
# else
      IF (EASTERN_EDGE) THEN
        DO j=JstrR,JendR
          BOUNDARY(ng)%zeta_east(j)=0.0_r8
        END DO
      END IF
      IF (WESTERN_EDGE) THEN
        DO j=JstrR,JendR
          BOUNDARY(ng)%zeta_west(j)=0.0_r8
        END DO
      END IF
      IF (SOUTHERN_EDGE) THEN
        DO i=IstrR,IendR
          BOUNDARY(ng)%zeta_south(i)=0.0_r8
        END DO
      END IF
      IF (NORTHERN_EDGE) THEN
        DO i=IstrR,IendR
          BOUNDARY(ng)%zeta_north(i)=0.0_r8
        END DO
      END IF
# endif
      RETURN
      END SUBROUTINE ana_fsobc_tile
#endif

#ifdef ANA_GRID
      SUBROUTINE ana_grid (ng, tile)
!
!=======================================================================
!                                                                      !
!  This routine sets model grid using an analytical expressions.       !
!                                                                      !
!  On Output:  stored in common blocks:                                !
!                                                                      !
!                           "grid"    (file grid.h)                    !
!                           "scalars" (file scalar.h)                  !
!                                                                      !
!     el       Length (m) of domain box in the ETA-direction.          !
!     f        Coriolis parameter (1/seconds) at RHO-points.           !
!     h        Bathymetry (meters; positive) at RHO-points.            !
!     hmin     Minimum depth of bathymetry (m).                        !
!     hmax     Maximum depth of bathymetry (m).                        !
!     pm       Coordinate transformation metric "m" (1/meters)         !
!              associated with the differential distances in XI        !
!              at RHO-points.                                          !
!     pn       Coordinate transformation metric "n" (1/meters)         !
!              associated with the differential distances in ETA.      !
!              at RHO-points.                                          !
!     xl       Length (m) of domain box in the XI-direction.           !
!     xp       XI-coordinates (m) at PSI-points.                       !
!     xr       XI-coordinates (m) at RHO-points.                       !
!     yp       ETA-coordinates (m) at PSI-points.                      !
!     yr       ETA-coordinates (m) at RHO-points.                      !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
!
      integer, intent(in) :: ng, tile

# include "tile.h"
!
      CALL ana_grid_tile (ng, Istr, Iend, Jstr, Jend,                   &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    GRID(ng) % angler,                            &
# if defined CURVGRID && defined UV_ADV
     &                    GRID(ng) % dmde,                              &
     &                    GRID(ng) % dndx,                              &
# endif
# ifdef ICESHELF
     &                    GRID(ng) % zice,                              &
# endif
# ifdef SPHERICAL
     &                    GRID(ng) % lonr,                              &
     &                    GRID(ng) % latr,                              &
# else
     &                    GRID(ng) % xp,                                &
     &                    GRID(ng) % xr,                                &
     &                    GRID(ng) % yp,                                &
     &                    GRID(ng) % yr,                                &
# endif
     &                    GRID(ng) % pn,                                &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % f,                                 &
     &                    GRID(ng) % h)
      RETURN
      END SUBROUTINE ana_grid
!
!***********************************************************************
      SUBROUTINE ana_grid_tile (ng, Istr, Iend, Jstr, Jend,             &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          angler,                                 &
# if defined CURVGRID && defined UV_ADV
     &                          dmde, dndx,                             &
# endif
# ifdef ICESHELF
     &                          zice,                                   &
# endif
# ifdef SPHERICAL
     &                          lonr, latr,                             &
# else
     &                          xp, xr, yp, yr,                         &
# endif
     &                          pn, pm, f, h)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_scalars

# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
!
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj

# ifdef ASSUMED_SHAPE
      real(r8), intent(out) :: angler(LBi:,LBj:)
#  if defined CURVGRID && defined UV_ADV
      real(r8), intent(out) :: dmde(LBi:,LBj:)
      real(r8), intent(out) :: dndx(LBi:,LBj:)
#  endif
#  ifdef ICESHELF
      real(r8), intent(out) :: zice(LBi:,LBj:)
#  endif
#  ifdef SPHERICAL
      real(r8), intent(out) :: lonr(LBi:,LBj:)
      real(r8), intent(out) :: latr(LBi:,LBj:)
#  else
      real(r8), intent(out) :: xp(LBi:,LBj:)
      real(r8), intent(out) :: xr(LBi:,LBj:)
      real(r8), intent(out) :: yp(LBi:,LBj:)
      real(r8), intent(out) :: yr(LBi:,LBj:)
#  endif
      real(r8), intent(out) :: pn(LBi:,LBj:)
      real(r8), intent(out) :: pm(LBi:,LBj:)
      real(r8), intent(out) :: f(LBi:,LBj:)
      real(r8), intent(out) :: h(LBi:,LBj:)
# else
      real(r8), intent(out) :: angler(LBi:UBi,LBj:UBj)
#  if defined CURVGRID && defined UV_ADV
      real(r8), intent(out) :: dmde(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: dndx(LBi:UBi,LBj:UBj)
#  endif
#  ifdef ICESHELF
      real(r8), intent(out) :: zice(LBi:UBi,LBj:UBj)
#  endif
#  ifdef SPHERICAL
      real(r8), intent(out) :: lonr(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: latr(LBi:UBi,LBj:UBj)
#  else
      real(r8), intent(out) :: xp(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: xr(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: yp(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: yr(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(out) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: f(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: h(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: Imin, Imax, Jmin, Jmax
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: NSUB, i, j, k

      real(r8), parameter :: twopi = 2.0_r8*pi

      real(r8) :: Esize, Xsize, beta, cff, depth, dth
      real(r8) :: dx, dy, f0, my_min, my_max, r, theta, val1, val2

# ifdef DISTRIBUTE
      real(r8), dimension(2) :: buffer
      character (len=3), dimension(2) :: op_handle
# endif
# ifdef WEDDELL
      real(r8) :: hwrk(-1:235), xwrk(-1:235), zwrk
# endif
      real(r8) :: wrkX(PRIVATE_2D_SCRATCH_ARRAY)
      real(r8) :: wrkY(PRIVATE_2D_SCRATCH_ARRAY)
!
# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set grid parameters:
!
!     Xsize    Length (m) of domain box in the XI-direction.
!     Esize    Length (m) of domain box in the ETA-direction.
!     depth    Maximum depth of bathymetry (m).
!     f0       Coriolis parameter, f-plane constant (1/s).
!     beta     Coriolis parameter, beta-plane constant (1/s/m).
!-----------------------------------------------------------------------
!
# if defined BASIN
      Xsize=3600.0E+03_r8
      Esize=2800.0E+03_r8
      depth=5000.0_r8
      f0=1.0E-04_r8
      beta=2.0E-11_r8
# elif defined BENCHMARK1 || defined BENCHMARK2 || defined BENCHMARK3
      Xsize=360.0_r8              ! degrees of longitude
      Esize=20.0_r8               ! degrees of latitude
      depth=4000.0_r8
      f0=-1.0E-04_r8
      beta=2.0E-11_r8
# elif defined BL_TEST
      Xsize=100.0E+03_r8
      Esize=5.0E+03_r8
      depth=47.5_r8
      f0=9.25E-04_r8
      beta=0.0_r8
# elif defined CANYON_A || defined CANYON_B
      Xsize=128.0E+03_r8
      Esize=96.0E+03_r8
      depth=4000.0_r8
      f0=1.0E-04_r8
      beta=0.0_r8
# elif defined COUPLING_TEST
      Xsize=6000.0_r8*REAL(Lm(ng),r8)
      Esize=6000.0_r8*REAL(Mm(ng),r8)
      depth=1500.0_r8
      f0=5.0E-05_r8
      beta=0.0_r8
# elif defined DOUBLE_GYRE
      Xsize=1000.0E+03_r8
      Esize=2000.0E+03_r8
      depth=500.0_r8
      f0=7.3E-05_r8
      beta=2.0E-11_r8
# elif defined ESTUARY_TEST
      Xsize=100000.0_r8
      Esize=300.0_r8
      depth=10.0_r8
      f0=0.0_r8
      beta=0.0_r8
# elif defined KELVIN
      Xsize=20000.0_r8*REAL(Lm(ng),r8)
      Esize=20000.0_r8*REAL(Mm(ng),r8)
      depth=100.0_r8
      f0=1.0E-04_r8
      beta=0.0_r8
# elif defined FLT_TEST
      Xsize=1.0E+03_r8*REAL(Lm(ng),r8)
      Esize=1.0E+03_r8*REAL(Mm(ng),r8)
      depth=10.0_r8
      f0=0.0_r8
      beta=0.0_r8
# elif defined GRAV_ADJ
      Xsize=64.0E+03_r8
      Esize=2.0E+03_r8
      depth=20.0_r8
      f0=0.0_r8
      beta=0.0_r8
# elif defined LAB_CANYON
      Xsize=0.55_r8                  ! width of annulus
      Esize=2.0_r8*pi                ! azimuthal length (radians)
      f0=4.0_r8*pi/25.0_r8
      beta=0.0_r8
# elif defined LAKE_SIGNELL
      Xsize=50.0e3_r8
      Esize=10.0e3_r8
      depth=18.0_r8
      f0=0.0E-04_r8
      beta=0.0_r8
# elif defined LMD_TEST
      Xsize=100.0E+03_r8
      Esize=100.0E+03_r8
      depth=50.0_r8
      f0=1.09E-04_r8
      beta=0.0_r8
# elif defined OVERFLOW
      Xsize=4.0E+03_r8
      Esize=200.0E+03_r8
      depth=4000.0_r8
      f0=0.0_r8
      beta=0.0_r8
# elif defined RIVERPLUME
      Xsize=58.5E+03_r8
      Esize=201.0E+03_r8
      depth=150.0_r8
      f0=1.0E-04_r8
      beta=0.0_r8
# elif defined SEAMOUNT
      Xsize=320.0E+03_r8
      Esize=320.0E+03_r8
      depth=5000.0_r8
      f0=1.0E-04_r8
      beta=0.0_r8
# elif defined SHELFRONT
      Xsize=20.0E+03_r8
      Esize=200.0E+03_r8
      depth=1660.0_r8
      f0=1.0E-04_r8
      beta=0.0_r8
# elif defined SOLITON
!!    Xsize=0.5_r8*REAL(Lm(ng),r8)
!!    Esize=0.5_r8*REAL(Mm(ng),r8)
      Xsize=48.0_r8
      Esize=16.0_r8
      depth=1.0_r8
      f0=0.0_r8
      beta=1.0_r8
      g=1.0_r8
# elif defined SED_TEST1
      Xsize=10000.0_r8
      Esize=500.0_r8
      depth=10.0_r8
      f0=0.0_r8
      beta=0.0_r8
# elif defined TEST_CHAN
      Xsize=10000.0_r8
      Esize=1000.0_r8
      depth=10.0_r8
      f0=0.0_r8
      beta=0.0_r8
# elif defined UPWELLING
      Xsize=1000.0_r8*REAL(Lm(ng),r8)
      Esize=1000.0_r8*REAL(Mm(ng),r8)
      depth=150.0_r8
      f0=-8.26E-05_r8
      beta=0.0_r8
# elif defined WEDDELL
      Xsize=4000.0_r8*REAL(Lm(ng),r8)
      Esize=4000.0_r8*REAL(Mm(ng),r8)
      depth=4500.0_r8
      f0=0.0_r8
      beta=0.0_r8
# elif defined WINDBASIN
      Xsize=2000.0_r8*REAL(Lm(ng),r8)
      Esize=1000.0_r8*REAL(Mm(ng),r8)
      depth=50.0_r8
      f0=1.0E-04_r8
      beta=0.0_r8
# else
      ANA_GRID: no values provided for Xsize, Esize, depth, f0, and beta.
# endif
!
!  Load grid parameters to global storage.
!
      IF (SOUTH_WEST_TEST) THEN
        xl(ng)=Xsize
        el(ng)=Esize
      END IF
!
!-----------------------------------------------------------------------
!  Compute the (XI,ETA) coordinates at PSI- and RHO-points.
!  Set grid spacing (m).
!-----------------------------------------------------------------------
!
!  Determine I- and J-ranges for computing grid data.  This ranges
!  are special in periodic boundary conditons since periodicity cannot
!  be imposed in the grid coordinates.
! 
      IF (WESTERN_EDGE) THEN
        Imin=Istr-1
      ELSE
        Imin=Istr
      END IF
      IF (EASTERN_EDGE) THEN
        Imax=Iend+1
      ELSE
        Imax=Iend
      END IF
      IF (SOUTHERN_EDGE) THEN
        Jmin=Jstr-1
      ELSE
        Jmin=Jstr
      END IF
      IF (NORTHERN_EDGE) THEN
        Jmax=Jend+1
      ELSE
        Jmax=Jend
      END IF

# if defined BENCHMARK1 || defined BENCHMARK2 || defined BENCHMARK3

!  Spherical coordinates set-up.
!
      dx=Xsize/REAL(Lm(ng),r8)
      dy=Esize/REAL(Mm(ng),r8)
      spherical=.true.
      DO j=Jmin,Jmax
        cff=-70.0_r8+dy*(REAL(j,r8)-0.5_r8)
        DO i=Imin,Imax
          lonr(i,j)=dx*(REAL(i,r8)-0.5_r8)
          latr(i,j)=cff
        END DO
      END DO
# elif defined LAB_CANYON

!  Polar coordinates set-up.
!
      dx=Xsize/REAL(Lm(ng),r8)
      dy=Esize/REAL(Mm(ng),r8)
!!    dth=twopi/REAL(Mm(ng),r8)               ! equal azimultal spacing
      dth=0.01_r8                             ! azimultal spacing
      cff=(4.0_r8*pi/(dth*REAL(Mm(ng),r8)))-1.0_r8   ! F
      DO j=Jmin,Jmax
        DO i=Imin,Imax
          r=0.35_r8+dx*REAL(i-1,r8)
          theta=-pi+                                                    &
     &          0.5_r8*dth*((cff+1.0_r8)*REAL(j-1,r8)+                  &
     &                      (cff-1.0_r8)*(REAL(Mm(ng),r8)/twopi)*       &
     &                      SIN(twopi*REAL(j-1,r8)/REAL(Mm(ng),r8)))
          xp(i,j)=r*COS(theta)
          yp(i,j)=r*SIN(theta)
          r=0.35_r8+dx*(REAL(i-1,r8)+0.5_r8)
          theta=-pi+                                                    &
     &          0.5_r8*dth*((cff+1.0_r8)*(REAL(j-1,r8)+0.5_r8)+         &
     &                      (cff-1.0_r8)*(REAL(Mm(ng),r8)/twopi)*       &
     &                      SIN(twopi*(REAL(j-1,r8)+0.5_r8)/            &
     &                          REAL(Mm(ng),r8)))
          xr(i,j)=r*COS(theta)
          yr(i,j)=r*SIN(theta)
        END DO
      END DO
# else
      dx=Xsize/REAL(Lm(ng),r8)
      dy=Esize/REAL(Mm(ng),r8)
      DO j=Jmin,Jmax
        DO i=Imin,Imax
#  ifdef BL_TEST
          dx=0.5_r8*(4000.0_r8/REAL(Lm(ng)+1,r8))*REAL(i,r8)+675.0_r8
#  endif
          xp(i,j)=dx*REAL(i-1,r8)
          xr(i,j)=dx*(REAL(i-1,r8)+0.5_r8)
          yp(i,j)=dy*REAL(j-1,r8)
          yr(i,j)=dy*(REAL(j-1,r8)+0.5_r8)
        END DO
      END DO
# endif
# ifdef DISTRIBUTE
#  ifdef SPHERICAL
      CALL mp_exchange (ng, Istr, Iend, Jstr, Jend,                     &
     &                  LBi, UBi, LBj, UBj, 1, 1,                       &
     &                  GHOST_POINTS, .false., .false., lonr(LBi,LBj)) 
      CALL mp_exchange (ng, Istr, Iend, Jstr, Jend,                     &
     &                  LBi, UBi, LBj, UBj, 1, 1,                       &
     &                  GHOST_POINTS, .false., .false., latr(LBi,LBj))
#  else
      CALL mp_exchange (ng, Istr, Iend, Jstr, Jend,                     &
     &                  LBi, UBi, LBj, UBj, 1, 1,                       &
     &                  GHOST_POINTS, .false., .false., xp(LBi,LBj))
      CALL mp_exchange (ng, Istr, Iend, Jstr, Jend,                     &
     &                  LBi, UBi, LBj, UBj, 1, 1,                       &
     &                  GHOST_POINTS, .false., .false., xr(LBi,LBj))
      CALL mp_exchange (ng, Istr, Iend, Jstr, Jend,                     &
     &                  LBi, UBi, LBj, UBj, 1, 1,                       &
     &                  GHOST_POINTS, .false., .false., yp(LBi,LBj))
      CALL mp_exchange (ng, Istr, Iend, Jstr, Jend,                     &
     &                  LBi, UBi, LBj, UBj, 1, 1,                       &
     &                  GHOST_POINTS, .false., .false., yr(LBi,LBj))
#  endif   
# endif
!
!-----------------------------------------------------------------------
! Compute coordinate transformation metrics at RHO-points "pm" and
! "pn"  (1/m) associated with the differential distances in XI and
! ETA, respectively.
!-----------------------------------------------------------------------
!
# define J_RANGE MIN(JstrR,Jstr-1),MAX(Jend+1,JendR)
# define I_RANGE MIN(IstrR,Istr-1),MAX(Iend+1,IendR)

# if defined BENCHMARK1 || defined BENCHMARK2 || defined BENCHMARK3
!
!  Spherical coordinates set-up.
!
      val1=REAL(Lm(ng),r8)/(2.0_r8*pi*Eradius)
      val2=REAL(Mm(ng),r8)*360.0_r8/(2.0_r8*pi*Eradius*Esize)
      DO j=J_RANGE
         cff=1.0_r8/COS((-70.0_r8+dy*(REAL(j,r8)-0.5_r8))*deg2rad)
        DO i=I_RANGE
          wrkX(i,j)=val1*cff
          wrkY(i,j)=val2
        END DO
      END DO
# elif defined LAB_CANYON
!
!  Polar coordinates set-up.
!
      DO j=J_RANGE
        DO i=I_RANGE
          r=0.35_r8+dx*(REAL(i-1,r8)+0.5_r8)
          theta=0.5_r8*dth*((cff+1.0_r8)+                               &
     &                      (cff-1.0_r8)*                               &
     &                      COS(twopi*REAL(j-1,r8)/REAL(Mm(ng),r8)))
          wrkX(i,j)=1.0_r8/dx
          wrkY(i,j)=1.0_r8/(r*theta)
        END DO
      END DO
# else
      DO j=J_RANGE
        DO i=I_RANGE
#  ifdef BL_TEST
          dx=0.5_r8*(4000.0_r8/REAL(Lm(ng)+1,r8))*REAL(i,r8)+675.0_r8
#  endif
          wrkX(i,j)=1.0_r8/dx
          wrkY(i,j)=1.0_r8/dy
        END DO
      END DO
# endif
# undef J_RANGE
# undef I_RANGE
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          pm(i,j)=wrkX(i,j)
          pn(i,j)=wrkY(i,j)
        END DO
      END DO
# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, pm)
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, pn)
# endif
# if (defined CURVGRID && defined UV_ADV)
!
!-----------------------------------------------------------------------
!  Compute d(1/n)/d(xi) and d(1/m)/d(eta) at RHO-points.
!-----------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          dndx(i,j)=0.5_r8*((1.0_r8/wrkY(i+1,j  ))-                     &
     &                      (1.0_r8/wrkY(i-1,j  )))
          dmde(i,j)=0.5_r8*((1.0_r8/wrkX(i  ,j+1))-                     &
     &                      (1.0_r8/wrkX(i  ,j-1)))
        END DO
      END DO
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, dndx)
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, dmde)
#  endif
# endif
!
!-----------------------------------------------------------------------
! Angle (radians) between XI-axis and true EAST at RHO-points.
!-----------------------------------------------------------------------
!
# if defined LAB_CANYON
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          theta=-pi+                                                    &
     &          0.5_r8*dth*((cff+1.0_r8)*(REAL(j-1,r8)+0.5_r8)+         &
     &                      (cff-1.0_r8)*(REAL(Mm(ng),r8)/twopi)*       &
     &                      SIN(twopi*(REAL(j-1,r8)+0.5_r8)/            &
     &                          REAL(Mm(ng),r8)))
          angler(i,j)=theta
        END DO
      END DO
# elif defined WEDDELL
      val1=90.0_r8*deg2rad
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          angler(i,j)=val1
        END DO
      END DO
# else
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          angler(i,j)=0.0_r8
        END DO
      END DO
# endif
# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, angler)
# endif
!
!-----------------------------------------------------------------------
!  Compute Coriolis parameter (1/s) at RHO-points.
!-----------------------------------------------------------------------
!
# if defined BENCHMARK1 || defined BENCHMARK2 || defined BENCHMARK3
      val1=2.0_r8*(2.0_r8*pi*366.25_r8/365.25_r8)/86400.0_r8
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          f(i,j)=val1*SIN(latr(i,j)*deg2rad)
        END DO
      END DO
# elif defined WEDDELL
      val1=10.4_r8/REAL(Lm(ng),r8)
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          f(i,j)=2.0_r8*7.2E-05_r8*                                     &
     &           SIN((-79.0_r8+REAL(i-1,r8)*val1)*deg2rad)
        END DO
      END DO
# else
      val1=0.5_r8*Esize
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          f(i,j)=f0+beta*(yr(i,j)-val1)
        END DO
      END DO
# endif
# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, f)
# endif
!
!-----------------------------------------------------------------------
!  Set bathymetry (meters; positive) at RHO-points.
!-----------------------------------------------------------------------
!
# if defined BENCHMARK1 || defined BENCHMARK2 || defined BENCHMARK3
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          h(i,j)=500.0_r8+1750.0_r8*(1.0+TANH((68.0_r8+latr(i,j))/dy))
        END DO
      END DO
# elif defined BL_TEST
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          val1=(xr(i,j)+500.0_r8)/15000.0_r8
          h(i,j)=14.0_r8+                                               &
     &           25.0_r8*(1.0_r8-EXP(-pi*xr(i,j)*1.0E-05_r8))-          &
     &           8.0_r8*exp(-val1*val1)
        END DO
      END DO
# elif defined CANYON_A || defined CANYON_B
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          val1=32000.0_r8-16000.0_r8*(SIN(pi*xr(i,j)/Xsize))**24
          h(i,j)=20.0_r8+0.5_r8*(depth-20.0_r8)*                        &
     &           (1.0_r8+TANH((yr(i,j)-val1)/10000.0_r8))
        END DO
      END DO
# elif defined ESTUARY_TEST
      DO j=JstrR,JendR
        DO i=IstrR,IendR
!!        h(i,j)=10.0_r8
          h(i,j)=5.0_r8+(Xsize-xr(i,j))/Xsize*5.0_r8
        END DO
      END DO 
# elif defined LAB_CANYON
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          r=0.35_r8+dx*(REAL(i-1,r8)+0.5_r8)
          theta=-pi+                                                    &
     &           0.5_r8*dth*((cff+1.0_r8)*(REAL(j-1,r8)+0.5_r8)+        &
     &                       (cff-1.0_r8)*(REAL(Mm(ng),r8)/twopi)*      &
     &                       SIN(dth*(REAL(j-1,r8)+0.5_r8)/             &
     &                           REAL(Mm(ng),r8)))
          val1=0.55_r8-0.15_r8*(COS(pi*theta*0.55_r8/0.2_r8)**2) !r_small
          val2=0.15_r8+0.15_r8*(COS(pi*theta*0.55_r8/0.2_r8)**2) !lambda
          IF (ABS(theta).ge.0.181818181818_r8) THEN
            IF (r.le.0.55_r8) THEN
              h(i,j)=0.025_r8                      ! shelf
            ELSE IF (r.ge.0.7_r8) THEN
              h(i,j)=0.125_r8                      ! deep
            ELSE
              h(i,j)=0.125_r8-0.1_r8*                                   &
     &               (COS(0.5_r8*pi*(r-0.55_r8)/0.15_r8)**2)
            END IF
          ELSE
            IF (r.le.val1) THEN
              h(i,j)=0.025_r8                      ! shelf
            ELSE IF (r.ge.0.7_r8) THEN
              h(i,j)=0.125_r8                      ! deep
            ELSE
              h(i,j)=0.125_r8-0.1_r8*                                   &
     &               (COS(0.5_r8*pi*(r-val1)/val2)**2)
            END IF
          END IF
        END DO
      END DO
# elif defined LAKE_SIGNELL
      DO j=JstrR,JendR
        DO i=IstrR,IendR
!!        h(i,j)=18.0_r8
          h(i,j)=18.0_r8-16.0_r8*FLOAT(Mm(ng)-j)/FLOAT(Mm(ng)-1)
        END DO
      END DO
# elif defined OVERFLOW
      val1=200.0_r8
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          h(i,j)=val1+0.5_r8*(depth-val1)*                              &
     &           (1.0_r8+TANH((yr(i,j)-100000.0_r8)/20000.0_r8))
        END DO
      END DO
# elif defined RIVERPLUME
      DO j=JstrR,JendR
        DO i=IstrR,MIN(5,IendR)
          h(i,j)=15.0_r8
        END DO
        DO i=MAX(6,IstrR),IendR
          h(i,j)=depth+REAL(Lm(ng)-i,r8)*(15.0_r8-depth)/               &
     &                 REAL(Lm(ng)-6,r8)
        END DO
      END DO
# elif defined SEAMOUNT
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          val1=(xr(i,j)-xr(Lm(ng)/2,Mm(ng)/2))/40000.0_r8
          val2=(yr(i,j)-yr(Lm(ng)/2,Mm(ng)/2))/40000.0_r8
          h(i,j)=depth-4500.0_r8*EXP(-(val1*val1+val2*val2))
        END DO
      END DO
# elif defined SHELFRONT
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          val1=yr(i,j)/1000.0_r8
          IF (val1.lt.50.0_r8) THEN
            h(i,j)=50.0_r8+2.0_r8*val1
          ELSE IF (val1.lt.60.0_r8) THEN
            h(i,j)=160.0_r8+1.5_r8*(val1-50.0_r8)**2                    &
     &                      0.1_r8*(val1-60.0_r8)**2
          ELSE IF (val1.lt.100.0_r8) THEN
            h(i,j)=310.0_r8+30.0_r8*(val1-60.0_r8)
          ELSE IF (val1.lt.110.0_r8) THEN
            h(i,j)=1660.0_r8-1.5_r8*(val1-110.0_r8)**2
          ELSE
            h(i,j)=1660.0_r8
          END IF
        END DO
      END DO
# elif defined TEST_CHAN
      DO j=JstrR,JendR
        DO i=IstrR,IendR
!!        h(i,j)=10.0_r8
          h(i,j)=10.396_r8-0.396_r8*REAL(i-1,r8)/REAL(Lm(ng)-1,r8)
        END DO
      END DO 
# elif defined UPWELLING
      DO j=JstrR,JendR
        IF (j.le.Mm(ng)/2) THEN
          val1=REAL(j,r8)
        ELSE
          val1=REAL(Mm(ng)+1-j,r8)
        END IF
        val2=MIN(depth,84.5_r8+66.526_r8*TANH((val1-10.0_r8)/7.0_r8))
        DO i=IstrR,IendR
          h(i,j)=val2
        END DO
      END DO
# elif defined WEDDELL
      val1=98.80_r8
      val2=0.8270_r8
      DO k=-1,26
        xwrk(k)=REAL(k-1,r8)*15.0_r8*1000.0_r8
        hwrk(k)=375.0_r8
      END DO
      DO k=27,232
        zwrk=-2.0_r8+REAL(k-1,r8)*0.020_r8
        xwrk(k)=(520.0_r8+val1+zwrk*val1+                               &
     &           val1*val2*LOG(COSH(zwrk)))*1000.0_r8
        hwrk(k)=-75.0_r8+2198.0_r8*(1.0_r8+val2*TANH(zwrk))
      END DO
      DO k=233,235
        xwrk(k)=(850.0_r8+REAL(k-228,r8)*50.0_r8)*1000.0_r8
        hwrk(k)=4000.0_r8
      END DO
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          DO k=1,234
            IF ((xwrk(k).le.xr(i,1)).and.(xr(i,1).lt.xwrk(k+1))) THEN
               cff=1.0_r8/(xwrk(k+1)-xwrk(k))
               h(i,j)=cff*(xwrk(k+1)-xr(i,j))*hwrk(k  )+                &
     &                cff*(xr(i,j)-xwrk(k  ))*hwrk(k+1)
            END IF
          END DO
        END DO
      END DO
# elif defined WINDBASIN
      DO i=IstrR,IendR
        val1=1;
        IF ((i-IstrR).lt.(IFIX(0.03_r8*(IendR-IstrR)))) THEN
          val1=1.0_r8-(REAL((i-IstrR+1)-                                &
     &                      IFIX(0.03_r8*(IendR-IstrR)),r8)/            &
     &                 REAL(0.03_r8*(IendR-IstrR),r8))**2
        END IF
        IF ((IendR-i).lt.(IFIX(0.03_r8*(IendR-IstrR)))) THEN
          val1=1.0_r8-(REAL((IendR-i+1)-                                &
     &                      IFIX(0.03_r8*(IendR-IstrR)),r8)/            &
     &                 REAL(0.03_r8*(IendR-IstrR),r8))**2
        END IF
        DO j=JstrR,JendR
         val2=2.0_r8*REAL(j-(Mm(ng)+1)/2,r8)/REAL(Mm(ng)+1,r8)
         h(i,j)=depth*(0.08_r8+0.92_r8*val1*(1.0_r8-val2*val2))
        END DO
      END DO
# else
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          h(i,j)=depth
        END DO
      END DO
# endif
# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, h)
# endif
!
! Determine minimum depth: first, determine minimum values of depth
! within each subdomain (stored as private variable cff), then
! determine global minimum by comparing these  subdomain minima.
!
      my_min=h(IstrR,JstrR)
      my_max=h(IstrR,JstrR)
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          my_min=MIN(my_min,h(i,j))
          my_max=MAX(my_max,h(i,j))
        END DO
      END DO
      IF (SOUTH_WEST_CORNER.and.                                        &
     &    NORTH_EAST_CORNER) THEN
        NSUB=1                           ! non-tiled application
      ELSE
        NSUB=NtileX(ng)*NtileE(ng)       ! tiled application
      END IF
!$OMP CRITICAL (H_RANGE)
      IF (tile_count.eq.0) THEN
        hmin(ng)=my_min
        hmax(ng)=my_max
      ELSE
        hmin(ng)=MIN(hmin(ng),my_min)
        hmax(ng)=MIN(hmax(ng),my_max)
      END IF
      tile_count=tile_count+1
      IF (tile_count.eq.NSUB) THEN
        tile_count=0
# ifdef DISTRIBUTE
        buffer(1)=hmin(ng)
        buffer(2)=hmax(ng)
        op_handle(1)='MIN'
        op_handle(2)='MAX'
        CALL mp_reduce (ng, 2, buffer, op_handle)
        hmin(ng)=buffer(1)
        hmax(ng)=buffer(2)
# endif
      END IF
!$OMP END CRITICAL (H_RANGE)
# ifdef ICESHELF
!
!-----------------------------------------------------------------------
!  Set depth of ice shelf (meters; negative) at RHO-points.
!-----------------------------------------------------------------------
!
#  ifdef WEDDELL
      val1=340.0_r8/16.0_r8
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          IF (i.gt.20) THEN
            zice(i,j)=0.0_r8
          ELSE IF (i.gt.4) THEN
            zice(i,j)=-340.0_r8+REAL(i-1,r8)*val1
          ELSE
            zice(i,j)=-340.0_r8
          END IF
        END DO
      END DO
#  else
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          zice(i,j)=0.0_r8
        END DO
      END DO
#  endif
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, zice)
#  endif
# endif
      RETURN
      END SUBROUTINE ana_grid_tile
#endif

#if defined ANA_HUMIDITY && \
  ( defined BULK_FLUXES || defined ECOSIM || \
   (defined ANA_SRFLUX  && defined ALBEDO) )
      SUBROUTINE ana_humid (ng, tile)
!
!=======================================================================
!                                                                      !
!  This routine sets surface air humidity (moisture) using an          !
!  analytical expression.  There three types of humidity:              !
!                                                                      !
!     1) Absolute humidity: density of water vapor.                    !
!     2) Specific humidity: ratio of the mass of water vapor to        !
!        the mass of moist air cointaining the vapor (g/kg)            !
!     3) Relative humidity: ratio of the actual mixing ratio to        !
!        saturation mixing ratio of the air at given temperature       !
!        and pressure (percentage).                                    !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_forces
!
      integer, intent(in) :: ng, tile

# include "tile.h"
!
      CALL ana_humid_tile (ng, Istr, Iend, Jstr, Jend,                  &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng) % Hair)
      RETURN
      END SUBROUTINE ana_humid
!
!***********************************************************************
      SUBROUTINE ana_humid_tile (ng, Istr, Iend, Jstr, Jend,            &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           Hair)
!***********************************************************************
!
      USE mod_param

# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
!
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj

# ifdef ASSUMED_SHAPE
      real(r8), intent(out) :: Hair(LBi:,LBj:)
# else
      real(r8), intent(out) :: Hair(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set analytical surface air humidity.
!-----------------------------------------------------------------------
!
# if defined ADRIATIC1 || defined ADRIATIC2
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          Hair(i,j)=0.7_r8
        END DO
      END DO
# elif defined BENCHMARK1 || defined BENCHMARK2 || defined BENCHMARK3
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          Hair(i,j)=0.8_r8
        END DO
      END DO
# elif defined BL_TEST
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          Hair(i,j)=0.776_r8
        END DO
      END DO
# else
      ANA_HUMIDITY: no values provided for HAIR.
# endif
# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, Hair)
# endif
      RETURN
      END SUBROUTINE ana_humid_tile
#endif

#ifdef ANA_INITIAL
      SUBROUTINE ana_initial (ng, tile)
!
!=======================================================================
!                                                                      !
!  This subroutine sets initial conditions for momentum and tracer     !
!  type variables using analytical expressions.                        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_ocean
!
      integer, intent(in) :: ng, tile

# include "tile.h"
!
      CALL ana_initial_tile (ng, Istr, Iend, Jstr, Jend,                &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       OCEAN(ng) % ubar,                          &
     &                       OCEAN(ng) % vbar,                          &
     &                       OCEAN(ng) % zeta,                          &
# ifdef SOLVE3D
     &                       OCEAN(ng) % u,                             &
     &                       OCEAN(ng) % v,                             &
     &                       OCEAN(ng) % t,                             &
     &                       GRID(ng) % z_r,                            &
# endif
# ifdef SPHERICAL
     &                       GRID(ng) % lonr,                           &
     &                       GRID(ng) % latr,                           &
# else
     &                       GRID(ng) % xr,                             &
     &                       GRID(ng) % yr,                             &
# endif
     &                       GRID(ng) % h)
      RETURN
      END SUBROUTINE ana_initial
!
!***********************************************************************
      SUBROUTINE ana_initial_tile (ng, Istr, Iend, Jstr, Jend,          &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             ubar, vbar, zeta,                    &
# ifdef SOLVE3D
     &                             u, v, t, z_r,                        &
# endif
# ifdef SPHERICAL
     &                             lonr, latr,                          &
# else
     &                             xr, yr,                              &
# endif
     &                             h)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars

# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
!
      USE exchange_2d_mod
#  ifdef SOLVE3D
      USE exchange_3d_mod
#  endif
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj

# ifdef ASSUMED_SHAPE
      real(r8), intent(out) :: ubar(LBi:,LBj:,:)
      real(r8), intent(out) :: vbar(LBi:,LBj:,:)
      real(r8), intent(out) :: zeta(LBi:,LBj:,:)
#  ifdef SOLVE3D
      real(r8), intent(out) :: u(LBi:,LBj:,:,:)
      real(r8), intent(out) :: v(LBi:,LBj:,:,:)
      real(r8), intent(out) :: t(LBi:,LBj:,:,:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
#  endif
#  ifdef SPHERICAL
      real(r8), intent(in) :: lonr(LBi:,LBj:)
      real(r8), intent(in) :: latr(LBi:,LBj:)
#  else
      real(r8), intent(in) :: xr(LBi:,LBj:)
      real(r8), intent(in) :: yr(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: h(LBi:,LBj:)
# else
      real(r8), intent(out) :: ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(out) :: vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(out) :: zeta(LBi:UBi,LBj:UBj,3)
#  ifdef SOLVE3D
      real(r8), intent(out) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(out) :: v(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(out) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
#  endif
#  ifdef SPHERICAL
      real(r8), intent(in) :: lonr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: latr(LBi:UBi,LBj:UBj)
#  else
      real(r8), intent(in) :: xr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: yr(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: Iless, Iplus, i, itrc, j, k
      real(r8) :: depth, dx, val1, val2, val3, val4, x, x0, y, y0

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Initial conditions for 2D momentum (m/s) components.
!-----------------------------------------------------------------------
!
# if defined SOLITON
      x0=2.0_r8*xl(ng)/3.0_r8
      y0=0.5_r8*el(ng)
      val1=0.395_r8
      val2=0.771_r8*(val1*val1)
      DO j=JstrR,JendR
        DO i=Istr,IendR
          x=0.5_r8*(xr(i-1,j)+xr(i,j))-x0
          y=0.5_r8*(yr(i-1,j)+yr(i,j))-y0
          val3=EXP(-val1*x)
          val4=val2*((2.0_r8*val3/(1.0_r8+(val3*val3)))**2)
          ubar(i,j,1)=0.25_r8*val4*(6.0_r8*y*y-9.0_r8)*                 &
     &                EXP(-0.5_r8*y*y)
          ubar(i,j,2)=ubar(i,j,1)
        END DO
      END DO
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          x=0.5_r8*(xr(i,j-1)+xr(i,j))-x0
          y=0.5_r8*(yr(i,j-1)+yr(i,j))-y0
          val3=EXP(-val1*x)
          val4=val2*((2.0_r8*val3/(1.0_r8+(val3*val3)))**2)
          vbar(i,j,1)=2.0_r8*val4*y*(-2.0_r8*val1*TANH(val1*x))*        &
     &                EXP(-0.5_r8*y*y)
          vbar(i,j,2)=vbar(i,j,1)
        END DO
      END DO
# elif defined LAKE_SIGNELL
      DO j=JstrR,JendR
        DO i=Istr,IendR
          ubar(i,j,1)=0.001_r8
          ubar(i,j,2)=ubar(i,j,1)
        END DO
      END DO
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          vbar(i,j,1)=0.001_r8
          vbar(i,j,2)=vbar(i,j,1)
        END DO
      END DO
# elif defined SED_TEST1
      val1=100.0_r8
      DO j=JstrR,JendR
        DO i=Istr,IendR
          ubar(i,j,1)=-10.0_r8/(10.0_r8+9.0E-06_r8*REAL(i,r8)*val1)
          ubar(i,j,2)=ubar(i,j,1)
        END DO
      END DO
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          vbar(i,j,1)=0.0_r8
          vbar(i,j,2)=vbar(i,j,1)
        END DO
      END DO
# else
      DO j=JstrR,JendR
        DO i=Istr,IendR
          ubar(i,j,1)=0.0_r8
          ubar(i,j,2)=ubar(i,j,1)
        END DO
      END DO
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          vbar(i,j,1)=0.0_r8
          vbar(i,j,2)=vbar(i,j,1)
        END DO
      END DO
# endif
# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      DO i=1,2
        CALL exchange_u2d_tile (ng, Istr, Iend, Jstr, Jend,             &
     &                          LBi, UBi, LBj, UBj, ubar(:,:,i))
        CALL exchange_v2d_tile (ng, Istr, Iend, Jstr, Jend,             &
     &                          LBi, UBi, LBj, UBj, vbar(:,:,i))
      END DO
# endif
!
!-----------------------------------------------------------------------
!  Initial conditions for free-surface (m).
!-----------------------------------------------------------------------
!
# if defined KELVIN
!!    val1=1.0_r8                               ! zeta0
!!    val2=2.0_r8*pi/(12.42_r8*3600.0_r8)       ! M2 Tide period
      DO j=JstrR,JendR
        DO i=IstrR,IendR
!!        zeta(i,j,1)=val1*                                             &
!!   &                EXP(-GRID(ng)%f(i,j)*GRID(ng)%yp(i,j)/            &
!!   &                    SQRT(g*GRID(ng)%h(i,j)))*                     &
!!   &                COS(val2*GRID(ng)%xp(i,j)/                        &
!!   &                    SQRT(g*GRID(ng)%h(i,j)))
          zeta(i,j,1)=0.0_r8
          zeta(i,j,2)=zeta(i,j,1)
        END DO
      END DO
# elif defined SOLITON
      x0=2.0_r8*xl(ng)/3.0_r8
      y0=0.5_r8*el(ng)
      val1=0.395_r8
      val2=0.771_r8*(val1*val1)
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          x=xr(i,j)-x0
          y=yr(i,j)-y0
          val3=EXP(-val1*x)
          val4=val2*((2.0_r8*val3/(1.0_r8+(val3*val3)))**2)
          zeta(i,j,1)=0.25_r8*val4*(6.0_r8*y*y+3.0_r8)*                 &
     &                EXP(-0.5_r8*y*y)
          zeta(i,j,2)=zeta(i,j,1)
        END DO
      END DO
# elif defined SED_TEST1
      val1=100.0_r8
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          zeta(i,j,1)=9.0E-06_r8*REAL(i,r8)*val1
          zeta(i,j,2)=zeta(i,j,1)
        END DO
      END DO
# elif defined TEST_CHAN
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          zeta(i,j,1)=-0.396_r8+0.396_r8*REAL(i-1,r8)/REAL(Lm(ng)-1,r8)
          zeta(i,j,2)=zeta(i,j,1)
        END DO
      END DO
# else
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          zeta(i,j,1)=0.0_r8
          zeta(i,j,2)=zeta(i,j,1)
        END DO
      END DO
# endif
# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      DO i=1,2
        CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,             &
     &                          LBi, UBi, LBj, UBj, zeta(:,:,i))
      END DO
# endif
# ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Initial conditions for 3D momentum components (m/s).
!-----------------------------------------------------------------------
!
#  if defined LAKE_SIGNELL
      DO k=1,N(ng)
       DO j=JstrR,JendR
         DO i=Istr,IendR
            u(i,j,k,1)=0.001_r8
            u(i,j,k,2)=u(i,j,k,1)
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            v(i,j,k,1)=0.001_r8
            v(i,j,k,2)=v(i,j,k,1)
          END DO
        END DO
      END DO
#  elif defined SED_TEST1
      DO k=1,N(ng)
       DO j=JstrR,JendR
         DO i=Istr,IendR
            u(i,j,k,1)=-1.0_r8*LOG((h(i,j)+z_r(i,j,k))/Zob(ng))/        &
     &                 (LOG(h(i,j)/Zob(ng))-1.0_r8+Zob(ng)/h(i,j))
            u(i,j,k,2)=u(i,j,k,1)
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            v(i,j,k,1)=0.0_r8
            v(i,j,k,2)=v(i,j,k,1)
          END DO
        END DO
      END DO
#  else
      DO k=1,N(ng)
       DO j=JstrR,JendR
         DO i=Istr,IendR
            u(i,j,k,1)=0.0_r8
            u(i,j,k,2)=u(i,j,k,1)
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            v(i,j,k,1)=0.0_r8
            v(i,j,k,2)=v(i,j,k,1)
          END DO
        END DO
      END DO
#  endif
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      DO i=1,2
        CALL exchange_u3d_tile (ng, Istr, Iend, Jstr, Jend,             &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          u(:,:,:,i))
        CALL exchange_v3d_tile (ng, Istr, Iend, Jstr, Jend,             &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          v(:,:,:,i))
      END DO
#  endif
!
!-----------------------------------------------------------------------
!  Initial conditions for tracer type variables.
!-----------------------------------------------------------------------
!
!  Set initial conditions for potential temperature (Celsius) and
!  salinity (PSU).
!
#  if defined BENCHMARK1 || defined BENCHMARK2 || defined BENCHMARK3
      val1=(44.69_r8/39.382_r8)**2
      val2=val1*(rho0*800.0_r8/g)*(5.0E-05_r8/((42.689_r8/44.69_r8)**2))
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            t(i,j,k,1,itemp)=val2*EXP(z_r(i,j,k)/800.0_r8)*             &
     &                       (0.6_r8-0.4_r8*TANH(z_r(i,j,k)/800.0_r8))
            t(i,j,k,1,isalt)=35.0_r8
            t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
            t(i,j,k,2,isalt)=t(i,j,k,1,isalt)
          END DO
        END DO
      END DO
#  elif defined BASIN
      val1=(44.69_r8/39.382_r8)**2
      val2=val1*(rho0*800.0_r8/g)*(5.0E-05_r8/((42.689_r8/44.69_r8)**2))
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            t(i,j,k,1,itemp)=val2*EXP(z_r(i,j,k)/800.0_r8)*             &
     &                       (0.6_r8-0.4_r8*TANH(z_r(i,j,k)/800.0_r8))
            t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
          END DO
        END DO
      END DO
#  elif defined BL_TEST
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            val1=TANH(1.1_r8*z_r(i,j,k)+11.0_r8)
            t(i,j,k,1,itemp)=17.25_r8+6.25_r8*val1
            t(i,j,k,1,isalt)=31.0_r8-0.75_r8*val1
            t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
            t(i,j,k,2,isalt)=t(i,j,k,1,isalt)
          END DO
        END DO
      END DO
#  elif defined CANYON_A
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            t(i,j,k,1,itemp)=T0(ng)
            t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
          END DO
        END DO
      END DO
#  elif defined CANYON_B
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            t(i,j,k,1,itemp)=3.488_r8*EXP(z_r(i,j,k)/800.0_r8)*         &
     &                       (1.0_r8-(2.0_r8/3.0_r8)*                   &
     &                               TANH(z_r(i,j,k)/800.0_r8))
            t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
          END DO
        END DO
      END DO
#  elif defined COUPLING_TEST
      val1=40.0_r8
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            t(i,j,k,1,itemp)=val1*EXP(z_r(i,j,k)/800.0_r8)*             &
     &                       (0.6_r8-0.4_r8*TANH(z_r(i,j,k)/800.0_r8))+ &
     &                       1.5_r8
            t(i,j,k,1,isalt)=35.0_r8
            t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
            t(i,j,k,2,isalt)=t(i,j,k,1,isalt)
          END DO
        END DO
      END DO
#  elif defined DOUBLE_GYRE
      val1=(44.69_r8/39.382_r8)**2
      val2=val1*(rho0*100.0_r8/g)*(5.0E-5_r8/((42.689_r8/44.69_r8)**2))
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            val3=T0(ng)+val2*EXP(z_r(i,j,k)/100.0_r8)*                  &
     &           (10.0_r8-0.4_r8*TANH(z_r(i,j,k)/100.0_r8))
            t(i,j,k,1,itemp)=val3*(1.0_r8-3.0_r8*yr(i,j)/(val3*el(ng)))
            t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
          END DO
        END DO
      END DO
#  elif defined ESTUARY_TEST
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            t(i,j,k,1,itemp)=10.0_r8
            IF (xr(i,j).le.30000.0_r8) then
              t(i,j,k,1,isalt)=30.0_r8
            ELSEIF (xr(i,j).le.80000.0_r8) then
              t(i,j,k,1,isalt)=(80000.0_r8-xr(i,j))/50000.0_r8*30.0_r8
            ELSE
              t(i,j,k,1,isalt)=0.0_r8
            END IF
            t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
            t(i,j,k,2,isalt)=t(i,j,k,1,isalt)
          END DO
        END DO
      END DO
#  elif defined FLT_TEST
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            t(i,j,k,1,itemp)=T0(ng)
            t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
          END DO
        END DO
      END DO
#  elif defined GRAV_ADJ
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=IstrR,MIN((Lm(ng)+1)/2,IendR)
            t(i,j,k,1,itemp)=T0(ng)+5.0_r8
            t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
          END DO
          DO i=MAX((Lm(ng)+1)/2+1,IstrR),IendR
            t(i,j,k,1,itemp)=T0(ng)
            t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
          END DO
!!        do i=IstrR,IendR
!!         if (i.lt.Lm(ng)/2) then
!!            t(i,j,k,1,itemp)=T0(ng)+5.0_r8
!!            t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
!!         elseif (i.eq.Lm(ng)/2) then
!!            t(i,j,k,1,itemp)=T0(ng)+4.0_r8
!!            t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
!!          elseif (i.eq.Lm(ng)/2+1) then
!!            t(i,j,k,1,itemp)=T0(ng)+1.0_r8
!!            t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
!!          else
!!            t(i,j,k,1,itemp)=T0(ng)
!!            t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
!!          endif
!!        enddo
        END DO
      END DO
#  elif defined LAB_CANYON
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            t(i,j,k,1,itemp)=-659.34183_r8*z_r(i,j,k)
          END DO
        END DO
      END DO
#  elif defined LAKE_SIGNELL
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            t(i,j,k,1,itemp)=10.0_r8
            t(i,j,k,1,isalt)=0.0_r8
            t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
            t(i,j,k,2,isalt)=t(i,j,k,1,isalt)
          END DO
        END DO
      END DO
#  elif defined LMD_TEST
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            t(i,j,k,1,itemp)=MIN(13.0_r8,                               &
     &                           7.0_r8+0.2_r8*(z_r(i,j,k)+50.0_r8))
            t(i,j,k,1,isalt)=35.0_r8
            t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
            t(i,j,k,2,isalt)=t(i,j,k,1,isalt)
          END DO
        END DO
      END DO
#  elif defined NJ_BIGHT
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            depth=z_r(i,j,k)
            IF (depth.ge.-15.0_r8) THEN
              t(i,j,k,1,itemp)= 2.049264257728403E+01_r8-depth*         &
     &                         (2.640850848793918E-01_r8+depth*         &
     &                         (2.751125328535212E-01_r8+depth*         &
     &                         (9.207489761648872E-02_r8+depth*         &
     &                         (1.449075725742839E-02_r8+depth*         &
     &                         (1.078215685912076E-03_r8+depth*         &
     &                         (3.240318053903974E-05_r8+               &
     &                          1.262826857690271E-07_r8*depth))))))
              t(i,j,k,1,isalt)= 3.066489149193135E+01_r8-depth*         &
     &                         (1.476725262946735E-01_r8+depth*         &
     &                         (1.126455760313399E-01_r8+depth*         &
     &                         (3.900923281871022E-02_r8+depth*         &
     &                         (6.939014937447098E-03_r8+depth*         &
     &                         (6.604436696792939E-04_r8+depth*         &
     &                         (3.191792361954220E-05_r8+               &
     &                          6.177352634409320E-07_r8*depth))))))
            ELSE
               t(i,j,k,1,itemp)=14.6_r8+                                &
     &                          6.70_r8*TANH(1.1_r8*depth+15.9_r8)
               t(i,j,k,1,isalt)=31.3_r8-                                &
     &                          0.55_r8*TANH(1.1_r8*depth+15.9_r8)
            END IF
            t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
            t(i,j,k,2,isalt)=t(i,j,k,1,isalt)
          END DO
        END DO
      END DO
#  elif defined OVERFLOW
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            t(i,j,k,1,itemp)=T0(ng)-0.5_r8*T0(ng)*(1.0_r8+              &
     &                       TANH((yr(i,j)-60000.0_r8)/2000.0_r8))
            t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
          END DO
        END DO
      END DO
#  elif defined RIVERPLUME
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            t(i,j,k,1,itemp)=T0(ng)+0.01_r8*REAL(k,r8)
            t(i,j,k,1,isalt)=S0(ng)
            t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
            t(i,j,k,2,isalt)=t(i,j,k,1,isalt)
          END DO
        END DO
      END DO
#  elif defined SEAMOUNT
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            t(i,j,k,1,itemp)=T0(ng)+7.5_r8*EXP(z_r(i,j,k)/1000.0_r8)
            t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
          END DO
        END DO
      END DO
#  elif defined SHELFRONT
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            t(i,j,k,1,itemp)=T0(ng)+2.5_r8*                             &
     &                       TANH((yr(i,j)-50000.0_r8)/20000.0_r8)
            t(i,j,k,1,isalt)=S0(ng)+                                    &
     &                       TANH((yr(i,j)-50000.0_r8)/20000.0_r8)
            t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
            t(i,j,k,2,isalt)=t(i,j,k,1,isalt)
          END DO
        END DO
      END DO
#  elif defined SED_TEST1
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            t(i,j,k,1,itemp)=20.0_r8
            t(i,j,k,1,isalt)=0.0_r8
            t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
            t(i,j,k,2,isalt)=t(i,j,k,1,isalt)
          END DO
        END DO
      END DO
#  elif defined UPWELLING
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            t(i,j,k,1,itemp)=14.0_r8+8.0_r8*EXP(z_r(i,j,k)/50.0_r8)-    &
     &                       T0(ng)
!!          t(i,j,k,1,itemp)=14.0_r8+(z_r(i,j,k)+75.0_r8)/150.0_r8+
!!   &                       4.0_r8*(1.0_r8+TANH((z_r(i,j,k)+35.0_r8)/
!!   &                                           6.5_r8))
            t(i,j,k,1,isalt)=1.0E-04_r8*yr(i,j)-S0(ng)
!!          t(i,j,k,1,isalt)=S0(ng)
!!          if (j.lt.Mm(ng)/2) then
!!            t(i,j,k,1,isalt)=0.0_r8
!!          elseif (j.eq.Mm(ng)/2) then
!!            t(i,j,k,1,isalt)=0.5_r8
!!          elseif (j.gt.Mm(ng)/2) then
!!            t(i,j,k,1,isalt)=1.0_r8
!!          endif
            t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
            t(i,j,k,2,isalt)=t(i,j,k,1,isalt)
          END DO
        END DO
      END DO
#  elif defined WINDBASIN
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            t(i,j,k,1,itemp)=20.0_r8                ! homogeneous
!!          t(i,j,k,1,itemp)=14.0_r8+8.0_r8*EXP(z_r(i,j,k)/50.0_r8)-    &
!!   &                       T0(ng)                 ! stratified
            t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
          END DO
        END DO
      END DO
#  else
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            t(i,j,k,1,itemp)=T0(ng)
            t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
#   ifdef SALINITY
            t(i,j,k,1,isalt)=S0(ng)
            t(i,j,k,2,isalt)=t(i,j,k,1,isalt)
#   endif
          END DO
        END DO
      END DO
#  endif
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      DO itrc=1,NAT
        DO i=1,2
          CALL exchange_r3d_tile (ng, Istr, Iend, Jstr, Jend,           &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            t(:,:,:,i,itrc))
        END DO
      END DO
#  endif
# endif
      RETURN
      END SUBROUTINE ana_initial_tile
#endif

#if defined ANA_M2CLIMA && defined M2CLIMATOLOGY
      SUBROUTINE ana_m2clima (ng, tile)
!
!=======================================================================
!                                                                      !
!  This routine sets analytical 2D momentum climatology fields.        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_clima
!
      integer, intent(in) :: ng, tile

# include "tile.h"
!
      CALL ana_m2clima_tile (ng, Istr, Iend, Jstr, Jend,                &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       CLIMA(ng) % ubarclm,                       &
     &                       CLIMA(ng) % vbarclm)
      RETURN
      END SUBROUTINE ana_m2clima
!
!***********************************************************************
      SUBROUTINE ana_m2clima_tile (ng, Istr, Iend, Jstr, Jend,          &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             ubarclm, vbarclm)
!***********************************************************************
!
      USE mod_param

# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
!
      USE exchange_2d_mod
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj

# ifdef ASSUMED_SHAPE
      real(r8), intent(out) :: ubarclm(LBi:,LBj:)
      real(r8), intent(out) :: vbarclm(LBi:,LBj:)
# else
      real(r8), intent(out) :: ubarclm(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: vbarclm(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set 2D momentum climatology.
!-----------------------------------------------------------------------
!
      DO j=JstrR,JendR
        DO i=Istr,IendR
          ubarclm(i,j)=???
        END DO
      END DO
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          vbarclm(i,j)=???
        END DO
      END DO
# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      CALL exchange_u2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, ubarclm)
      CALL exchange_v2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, vbarclm)
# endif
      RETURN
      END SUBROUTINE ana_m2clima_tile
#endif

#ifdef ANA_M2OBC
      SUBROUTINE ana_m2obc (ng, tile)
!
!=======================================================================
!                                                                      !
!  This routine sets 2D momentum open boundary conditions using        !
!  analytical expressions.                                             !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_ocean
      USE mod_stepping
!
      integer, intent(in) :: ng, tile

# include "tile.h"
!    
      CALL ana_m2obc_tile (ng, Istr, Iend, Jstr, Jend,                  &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     knew(ng),                                    &
     &                     GRID(ng) % angler,                           &
     &                     GRID(ng) % h,                                &
     &                     GRID(ng) % pm,                               &
     &                     GRID(ng) % pn,                               &
     &                     OCEAN(ng) % zeta)
      RETURN
      END SUBROUTINE ana_m2obc
!
!***********************************************************************
      SUBROUTINE ana_m2obc_tile (ng, Istr, Iend, Jstr, Jend,            &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           knew,                                  &
     &                           angler, h, pm, pn, zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
      USE mod_grid
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: knew

# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: angler(LBi:,LBj:)   
      real(r8), intent(in) :: h(LBi:,LBj:)   
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: zeta(LBi:,LBj:,:)
# else
      real(r8), intent(in) :: angler(LBi:UBi,LBj:UBj)   
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)   
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: zeta(LBi:UBi,LBj:UBj,3)
# endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j
      real(r8) :: angle, fac, major, minor, omega, phase, val
#if defined ESTUARY_TEST
      real(r8) :: my_area, my_flux, tid_flow, riv_flow, cff1, cff2,     &
     &            model_flux, cff      
#endif

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  2D momentum open boundary conditions.
!-----------------------------------------------------------------------
!
#if defined ESTUARY_TEST
        cff1=0.40_r8
        cff2=0.08_r8
        tid_flow=cff1*300.0_r8*10.0_r8
        riv_flow=cff2*300.0_r8*5.0_r8
        IF (WESTERN_EDGE) THEN
          my_area=0._r8
          my_flux=0._r8
          DO j=Jstr,Jend
            cff=0.5_r8*(zeta(Istr-1,j,knew)+h(Istr-1,j)+                &
     &                  zeta(Istr  ,j,knew)+h(Istr  ,j))/pn(Istr,j)
            my_area=my_area+cff
          END DO
          my_flux=-tid_flow*SIN(2.0_r8*pi*time(ng)/                     &
     &                          (12.0_r8*3600.0_r8))-riv_flow
          DO j=Jstr,Jend
            BOUNDARY(ng)%ubar_west(j)=my_flux/my_area
            BOUNDARY(ng)%vbar_west(j)=0.0_r8
          END DO
        END IF
        IF (EASTERN_EDGE) THEN
          my_area=0._r8
          my_flux=0._r8
          DO j=Jstr,Jend
            cff=0.5_r8*(zeta(Iend  ,j,knew)+h(Iend  ,j)+                &  
     &                  zeta(Iend+1,j,knew)+h(Iend+1,j))/pn(Iend,j)
            my_area=my_area+cff
          END DO
          my_flux=-riv_flow
          DO j=Jstr,Jend
            BOUNDARY(ng)%ubar_east(j)=my_flux/my_area
            BOUNDARY(ng)%vbar_east(j)=0.0_r8
          END DO
        END IF
# elif defined KELVIN
      fac=1.0_r8                                ! zeta0
      omega=2.0_r8*pi/(12.42_r8*3600.0_r8)      ! M2 Tide period
      val=fac*SIN(omega*time(ng))
      IF (WESTERN_EDGE) THEN
        DO j=JstrR,JendR
          cff=SQRT(g/GRID(ng)%h(Istr-1,j))
          BOUNDARY(ng)%ubar_west(j)=val*cff*EXP(-GRID(ng)%f(Istr-1,j)*  &
     &                                           GRID(ng)%yp(Istr-1,j)/ &
     &                                          cff)
        END DO
        DO j=Jstr,JendR
          vbar_west(j)=0.0_r8
        END DO
      END IF
      IF (EASTERN_EDGE) THEN
        DO j=JstrR,JendR
          cff=1.0_r8/SQRT(g/GRID(ng)%h(Iend,j))
          val=fac*EXP(-GRID(ng)%f(Iend,j)*GRID(ng)%yp(Istr-1,j)*cff)
          BOUNDARY(ng)%ubar_east(j)=val*SIN(omega*GRID(ng)%xp(Iend,j)*  &
     &                                      cff-omega*time(ng))
        END DO
        DO j=Jstr,JendR
          vbar_east(j)=0.0_r8
        END DO
      END IF
# elif defined SED_TEST1
      IF (WESTERN_EDGE) THEN
        DO j=JstrR,JendR
          val=0.5_r8*(zeta(Istr-1,j,knew)+h(Istr-1,j)+                  &
     &                zeta(Istr  ,j,knew)+h(Istr  ,j))
          BOUNDARY(ng)%ubar_west(j)=-10.0_r8/val
        END DO
        DO j=Jstr,JendR
          vbar_west(j)=0.0_r8
        END DO
      END IF
      IF (EASTERN_EDGE) THEN
        DO j=JstrR,JendR
          val=0.5_r8*(zeta(Iend  ,j,knew)+h(Iend  ,j)+                  &
     &                zeta(Iend+1,j,knew)+h(Iend+1,j))
          BOUNDARY(ng)%ubar_east(j)=-10.0_r8/val
        END DO
        DO j=Jstr,JendR
          BOUNDARY(ng)%vbar_east(j)=0.0_r8
        END DO
      END IF
# elif defined TEST_CHAN
      IF (EASTERN_EDGE) THEN
        DO j=JstrR,JendR
          val=0.5_r8*(zeta(Iend  ,j,knew)+h(Iend  ,j)+                  &
     &                zeta(Iend+1,j,knew)+h(Iend+1,j))
          BOUNDARY(ng)%ubar_east(j)=-10.0_r8/val
        END DO
        DO j=Jstr,JendR
          BOUNDARY(ng)%vbar_east(j)=0.0_r8
        END DO
      END IF
# elif defined WEDDELL
      IF (WESTERN_EDGE) THEN
        fac=TANH((tdays(ng)-dstart)/1.0_r8)
        omega=2.0_r8*pi*time(ng)/(12.42_r8*3600.0_r8)  !  M2 Tide period
        minor=0.0143_r8+(0.0143_r8+0.010_r8)/REAL(Iend+1,r8)
        major=0.1144_r8+(0.1144_r8-0.013_r8)/REAL(Iend+1,r8)
        phase=(318.0_r8+(318.0_r8-355.0_r8)/REAL(Iend+1,r8))*deg2rad
        angle=(125.0_r8+(125.0_r8- 25.0_r8)/REAL(Iend+1,r8))*deg2rad
        DO j=JstrR,JendR
          val=0.5_r8*(angler(Istr-1,j)+angler(Istr,j))
          BOUNDARY(ng)%ubar_west(j)=fac*(major*COS(angle-val)*          &
     &                                         COS(omega-phase)-        &
     &                                   minor*SIN(angle-val)*          &
     &                                         SIN(omega-phase))
        END DO
        DO j=Jstr,JendR
          val=0.5_r8*(angler(Istr-1,j-1)+angler(Istr-1,j))
          BOUNDARY(ng)%vbar_west(j)=fac*(major*SIN(angle-val)*          &
     &                                         COS(omega-phase)-        &
     &                                   minor*SIN(angle-val)*          &
     &                                         COS(omega-phase))
        END DO
      END IF
      IF (EASTERN_EDGE) THEN
        fac=TANH((tdays(ng)-dstart)/1.0_r8)
        omega=2.0_r8*pi*time(ng)/(12.42_r8*3600.0_r8)  !  M2 Tide period
        minor=0.0143_r8+(0.0143_r8+0.010_r8)
        major=0.1144_r8+(0.1144_r8-0.013_r8)
        phase=(318.0_r8+(318.0_r8-355.0_r8))*deg2rad
        angle=(125.0_r8+(125.0_r8- 25.0_r8))*deg2rad
        DO j=JstrR,JendR
          val=0.5_r8*(angler(Iend,j)+angler(Iend+1,j))
          BOUNDARY(ng)%ubar_east(j)=fac*(major*COS(angle-val)*          &
     &                                         COS(omega-phase)-        &
     &                                   minor*SIN(angle-val)*          &
     &                                         SIN(omega-phase))
        END DO
        DO j=Jstr,JendR
          val=0.5_r8*(angler(Iend+1,j-1)+angler(Iend+1,j))
          BOUNDARY(ng)%vbar_east(j)=fac*(major*SIN(angle-val)*          &
     &                                         COS(omega-phase)-        &
     &                                   minor*SIN(angle-val)*          &
     &                                         COS(omega-phase))
        END DO
      END IF
# else
      IF (EASTERN_EDGE) THEN
        DO j=JstrR,JendR
          BOUNDARY(ng)%ubar_east(j)=0.0_r8
        END DO
        DO j=Jstr,JendR
          BOUNDARY(ng)%vbar_east(j)=0.0_r8
        END DO
      END IF
      IF (WESTERN_EDGE) THEN
        DO j=JstrR,JendR
          BOUNDARY(ng)%ubar_west(j)=0.0_r8
        END DO
        DO j=Jstr,JendR
          BOUNDARY(ng)%vbar_west(j)=0.0_r8
        END DO
      END IF
      IF (SOUTHERN_EDGE) THEN
        DO i=Istr,IendR
          BOUNDARY(ng)%ubar_south(i)=0.0_r8
        END DO
        DO i=IstrR,IendR
          BOUNDARY(ng)%vbar_south(i)=0.0_r8
        END DO
      END IF
      IF (NORTHERN_EDGE) THEN
        DO i=Istr,IendR
          BOUNDARY(ng)%ubar_north(i)=0.0_r8
        END DO
        DO i=IstrR,IendR
          BOUNDARY(ng)%vbar_north(i)=0.0_r8
        END DO
      END IF
# endif
      RETURN
      END SUBROUTINE ana_m2obc_tile
#endif

#if defined ANA_M3CLIMA && defined M3CLIMATOLOGY
      SUBROUTINE ana_m3clima (ng, tile)
!
!=======================================================================
!                                                                      !
!  This routine sets analytical 3D momentum climatology fields.        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_clima
!
      integer, intent(in) :: ng, tile

# include "tile.h"
!
      CALL ana_m3clima_tile (ng, Istr, Iend, Jstr, Jend,                &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       CLIMA(ng) % uclm,                          &
     &                       CLIMA(ng) % vclm)
      RETURN
      END SUBROUTINE ana_m3clima
!
!***********************************************************************
      SUBROUTINE ana_m3clima_tile (ng, Istr, Iend, Jstr, Jend,          &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             uclm, vclm)
!***********************************************************************
!
      USE mod_param

# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
!
      USE exchange_3d_mod
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj

# ifdef ASSUMED_SHAPE
      real(r8), intent(out) :: uclm(LBi:,LBj:,:)
      real(r8), intent(out) :: vclm(LBi:,LBj:,:)
# else
      real(r8), intent(out) :: uclm(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(out) :: vclm(LBi:UBi,LBj:UBj,N(ng))
# endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j, k

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set 3D momentum climatology.
!-----------------------------------------------------------------------
!
# ifdef CBLAST
      DO k=1,N
        DO j=JstrR,JendR
          DO i=Istr,IendR
            uclm(i,j,k)=0.0_r8
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            vclm(i,j,k)=0.0_r8
          END DO
        END DO
      END DO
# else
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=Istr,IendR
            uclm(i,j,k)=???
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            vclm(i,j,k)=???
          END DO
        END DO
      END DO
# endif
# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      CALL exchange_u3d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, 1, N(ng), uclm)
      CALL exchange_v3d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, 1, N(ng), vclm)
# endif
      RETURN
      END SUBROUTINE ana_m3clima_tile
#endif

#ifdef ANA_M3OBC
      SUBROUTINE ana_m3obc (ng, tile)
!
!=======================================================================
!                                                                      !
!  This routine sets 3D momentum open boundary conditions using        !
!  analytical expressions.                                             !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_boundary
!
      integer, intent(in) :: ng, tile

# include "tile.h"
!
      CALL ana_m3obc_tile (ng, Istr, Iend, Jstr, Jend,                  &
     &                     LBi, UBi, LBj, UBj)
      RETURN
      END SUBROUTINE ana_m3obc
!
!***********************************************************************
      SUBROUTINE ana_m3obc_tile (ng, Istr, Iend, Jstr, Jend,            &
     &                           LBi, UBi, LBj, UBj)
!***********************************************************************
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j, k
      real(r8) :: fac, val

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  3D momentum open boundary conditions.
!-----------------------------------------------------------------------
!
# if defined SED_TEST1
      IF (WESTERN_EDGE) THEN
        fac=5.0E-06_r8
        DO k=1,N(ng)
          DO j=JstrR,JendR
            val=0.5_r8*(zeta(0 ,j,knew)+h(0 ,j)+                        &
     &                  zeta(1 ,j,knew)+h(1 ,j))
            BOUNDARY(ng)%u_west(j,k)=-LOG((val+0.5*(z_r(Istr-1,j,k)+    &
     &                                              z_r(Istr  ,j,k)))/  &
     &                                    fac)/                         &
     &                               (LOG(val/fac)-1.0_r8+fac/val)
          END DO
          DO j=Jstr,JendR
            BOUNDARY(ng)%v_west(j,k)=0.0_r8
          END DO
        END DO
      END IF
      IF (EASTERN_EDGE) THEN
        fac=5.0E-06_r8
        DO k=1,N(ng)
          DO j=JstrR,JendR
            val=0.5_r8*(zeta(Iend  ,j,knew)+h(Iend  ,j)+                &
     &                  zeta(Iend+1,j,knew)+h(Iend+1,j))
            BOUNDARY(ng)%u_east(j,k)=-LOG((val+0.5*(z_r(Iend  ,j,k)+    &
     &                                              z_r(Iend+1,j,k)))/  &
     &                                    fac)/                         &
     &                               (LOG(val/fac)-1.0_r8+fac/val)
          END DO
          DO j=Jstr,JendR
            BOUNDARY(ng)%v_east(j,k)=0.0_r8
          END DO
        END DO
      END IF
# else
      IF (EASTERN_EDGE) THEN
        DO k=1,N(ng)
          DO j=JstrR,JendR
            BOUNDARY(ng)%u_east(j,k)=0.0_r8
          END DO
          DO j=Jstr,JendR
            BOUNDARY(ng)%v_east(j,k)=0.0_r8
          END DO
        END DO
      END IF
      IF (WESTERN_EDGE) THEN
        DO k=1,N(ng)
          DO j=JstrR,JendR
            BOUNDARY(ng)%u_west(j,k)=0.0_r8
          END DO
          DO j=Jstr,JendR
            BOUNDARY(ng)%v_west(j,k)=0.0_r8
          END DO
        END DO
      END IF
      IF (SOUTHERN_EDGE) THEN
        DO k=1,N(ng)
          DO i=Istr,IendR
            BOUNDARY(ng)%u_south(i,k)=0.0_r8
          END DO
          DO i=IstrR,IendR
            BOUNDARY(ng)%v_south(i,k)=0.0_r8
          END DO
        END DO
      END IF
      IF (NORTHERN_EDGE) THEN
        DO k=1,N(ng)
          DO i=Istr,IendR
            BOUNDARY(ng)%u_north(i,k)=0.0_r8
          END DO
          DO i=IstrR,IendR
            BOUNDARY(ng)%v_north(i,k)=0.0_r8
          END DO
        END DO
      END IF
# endif
      RETURN
      END SUBROUTINE ana_m3obc_tile
#endif

#if defined ANA_GRID && defined MASKING
      SUBROUTINE ana_mask (ng, tile)
!
!=======================================================================
!                                                                      !
!  This subroutine sets analytical Land/Sea masking.                   !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
!
      integer, intent(in) :: ng, tile

# include "tile.h"
!
      CALL ana_mask_tile (ng, Istr, Iend, Jstr, Jend,                   &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    GRID(ng) % pmask,                             &
     &                    GRID(ng) % rmask,                             &
     &                    GRID(ng) % umask,                             &
     &                    GRID(ng) % vmask)
      RETURN
      END SUBROUTINE ana_mask
!
!***********************************************************************
      SUBROUTINE ana_mask_tile (ng, Istr, Iend, Jstr, Jend,             &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          pmask, rmask, umask, vmask)
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

# ifdef ASSUMED_SHAPE
      real(r8), intent(out) :: pmask(LBi:,LBj:)
      real(r8), intent(out) :: rmask(LBi:,LBj:)
      real(r8), intent(out) :: umask(LBi:,LBj:)
      real(r8), intent(out) :: vmask(LBi:,LBj:)
# else
      real(r8), intent(out) :: pmask(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: vmask(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j
      real(r8) :: mask(PRIVATE_2D_SCRATCH_ARRAY)

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set Land/Sea mask of RHO-points: Land=0, Sea=1.
!-----------------------------------------------------------------------
!
!  Notice that private scratch array "mask" is used to allow
!  computation within a parallel loop.
!
# if defined FLT_TEST
      DO j=Jstr-2,Jend+2
        DO i=Istr-2,Iend+2
          mask(i,j)=1.0_r8
          IF (j.eq.1 ) mask(i,j)=0.0_r8
          IF (j.eq.Mm(ng)) mask(i,j)=0.0_r8
          IF ((i.ge.((Lm(ng)+1)/2)).and.                                &
     &        (i.le.((Lm(ng)+1)/2+1)).and.                              &
     &        (j.ge.((Mm(ng)+1)/2)).and.                                &
     &        (j.le.((Mm(ng)+1)/2+1))) mask(i,j)=0.0_r8
        END DO
      END DO
# elif defined RIVERPLUME
      DO j=Jstr-2,Jend+2
        DO i=Istr-2,Iend+2
          mask(i,j)=1.0_r8
        END DO
      END DO
      DO i=Istr-2,MIN(5,Iend+2)
        DO j=Jstr-2,MIN(Mm(ng)-18,Jend+2)
          mask(i,j)=0.0_r8
        END DO
        DO j=MAX(Jstr-2,Mm(ng)-16),Jend+2
          mask(i,j)=0.0_r8
        END DO
      END DO
# else
      ANA_MASK: no values provided for RMASK.
# endif
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          rmask(i,j)=mask(i,j)
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Compute Land/Sea mask of U- and V-points.
!-----------------------------------------------------------------------
!
      DO j=JstrR,JendR
        DO i=Istr,IendR
          umask(i,j)=mask(i-1,j)*mask(i,j)
        END DO
      END DO
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          vmask(i,j)=mask(i,j-1)*mask(i,j)
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Compute Land/Sea mask of PSI-points.  If applicable,  set
!  slipperiness mask (Sea=1, Land=0, boundary=1-gamma2).
!-----------------------------------------------------------------------
!
      DO j=Jstr,JendR
        DO i=Istr,IendR
          pmask(i,j)=mask(i-1,j-1)*mask(i,j-1)*                         &
     &               mask(i-1,j  )*mask(i,j  )
        END DO
      END DO
      IF (gamma2.lt.0.0_r8) THEN
        DO j=Jstr,JendR
          DO i=Istr,IendR
            pmask(i,j)=2.0_r8-pmask(i,j)
          END DO
        END DO
      END IF
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, rmask)
      CALL exchange_p2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, pmask)
      CALL exchange_u2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, umask)
      CALL exchange_v2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, vmask)
#  endif
      RETURN
      END SUBROUTINE ana_mask_tile
#endif

#if defined ANA_PAIR && (defined BULK_FLUXES || defined ECOSIM)
      SUBROUTINE ana_pair (ng, tile)
!
!=======================================================================
!                                                                      !
!  This routine sets surface air pressure (mb) using an analytical     !
!  expression.                                                         !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_forces
!
      integer, intent(in) :: ng, tile

# include "tile.h"
!
      CALL ana_pair_tile (ng, Istr, Iend, Jstr, Jend,                   &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    FORCES(ng) % Pair)
      RETURN
      END SUBROUTINE ana_pair
!
!***********************************************************************
      SUBROUTINE ana_pair_tile (ng, Istr, Iend, Jstr, Jend,             &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Pair)
!***********************************************************************
!
      USE mod_param

# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
!
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj

# ifdef ASSUMED_SHAPE
      real(r8), intent(out) :: Pair(LBi:,LBj:)
# else
      real(r8), intent(out) :: Pair(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set analytical surface air pressure (mb).
!  (1 mb = 100 Pa = 1 hPa,  1 bar = 1.0e+5 N/m2 = 1.0e+5 dynes/cm2).
!-----------------------------------------------------------------------
!
# if defined BENCHMARK1 || defined BENCHMARK2 || defined BENCHMARK3
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          Pair(i,j)=1025.0_r8
        END DO
      END DO
# elif defined BL_TEST
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          Pair(i,j)=1013.48_r8
        END DO
      END DO
# else
      ANA_PAIR: no values provided for PAIR.
# endif
# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, Pair)
# endif
      RETURN
      END SUBROUTINE ana_pair_tile
#endif

#if (defined TS_PSOURCE || defined UV_PSOURCE) && defined ANA_PSOURCE
      SUBROUTINE ana_psource (ng)
!
!=======================================================================
!                                                                      !
!  This subroutine sets analytical tracer and mass point Sources       !
!  and/or Sinks.  River runoff can be consider as a point source.      !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_ocean
      USE mod_sources
      USE mod_stepping
!
      integer, intent(in) :: ng

      integer :: LBi, UBi, LBj, UBj
!
      LBi=lbound(GRID(ng)%h,DIM=1)
      UBi=ubound(GRID(ng)%h,DIM=1)
      LBj=lbound(GRID(ng)%h,DIM=2)
      UBj=ubound(GRID(ng)%h,DIM=2)
!
      CALL ana_psource_grid (ng, LBi, UBi, LBj, UBj,                    &
     &                       nnew(ng), knew(ng), Nsrc(ng),              &
     &                       OCEAN(ng) % zeta,                          &
     &                       OCEAN(ng) % ubar,                          &
     &                       OCEAN(ng) % vbar,                          &
# ifdef SOLVE3D
     &                       OCEAN(ng) % u,                             &
     &                       OCEAN(ng) % v,                             &
     &                       GRID(ng) % z_w,                            &
# endif
     &                       GRID(ng) % h,                              &
     &                       GRID(ng) % on_u,                           &
     &                       SOURCES(ng) % Isrc,                        &
     &                       SOURCES(ng) % Jsrc,                        &
     &                       SOURCES(ng) % Lsrc,                        &
     &                       SOURCES(ng) % Dsrc,                        &
# ifdef SOLVE3D
#  ifdef UV_PSOURCE
     &                       SOURCES(ng) % Qshape,                      &
     &                       SOURCES(ng) % Qsrc,                        &
#  endif
#  ifdef TS_PSOURCE
     &                       SOURCES(ng) % Tsrc,                        &
#  endif
# endif
     &                       SOURCES(ng) % Qbar)
      RETURN
      END SUBROUTINE ana_psource
!
!***********************************************************************
      SUBROUTINE ana_psource_grid (ng, LBi, UBi, LBj, UBj,              &
     &                             nnew, knew, Nsrc,                    &
     &                             zeta, ubar, vbar,                    &
# ifdef SOLVE3D
     &                             u, v, z_w,                           &
# endif
     &                             h, on_u,                             &
     &                             Isrc, Jsrc, Lsrc, Dsrc,              &
# ifdef SOLVE3D
#  ifdef UV_PSOURCE
     &                             Qshape, Qsrc,                        &
#  endif
#  ifdef TS_PSOURCE
     &                             Tsrc,                                &
#  endif
# endif
     &                             Qbar)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, LBi, UBi, LBj, UBj
      integer, intent(in) :: nnew, knew

      integer, intent(out) :: Nsrc

# ifdef ASSUMED_SHAPE
      logical, intent(out) :: Lsrc(:,:)

      integer, intent(out) :: Isrc(:)
      integer, intent(out) :: Jsrc(:)

      real(r8), intent(in) :: zeta(LBi:,LBj:,:)
      real(r8), intent(in) :: ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: vbar(LBi:,LBj:,:)
#  ifdef SOLVE3D
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
#  endif
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)

      real(r8), intent(out) :: Dsrc(:)
      real(r8), intent(out) :: Qbar(:)
#  ifdef SOLVE3D
#   ifdef UV_PSOURCE
      real(r8), intent(out) :: Qshape(:,:)
      real(r8), intent(out) :: Qsrc(:,:)
#   endif
#   ifdef TS_PSOURCE
      real(r8), intent(out) :: Tsrc(:,:,:)
#   endif
#  endif
# else
      logical, intent(out) :: Lsrc(Msrc,NT(ng))

      integer, intent(out) :: Isrc(Msrc)
      integer, intent(out) :: Jsrc(Msrc)

      real(r8), intent(in) :: zeta(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: vbar(LBi:UBi,LBj:UBj,3)
#  ifdef SOLVE3D
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
#  endif
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)

      real(r8), intent(out) :: Dsrc(Msrc)
      real(r8), intent(out) :: Qbar(Msrc)
#  ifdef SOLVE3D
#   ifdef UV_PSOURCE
      real(r8), intent(out) :: Qshape(Msrc,N(ng))
      real(r8), intent(out) :: Qsrc(Msrc,N(ng))
#   endif
#   ifdef TS_PSOURCE
      real(r8), intent(out) :: Tsrc(Msrc,N(ng),NT(ng))
#   endif
#  endif
# endif
!
!  Local variable declarations.
!
      integer :: is, i, j, k
      real(r8) :: fac, my_area
!
!-----------------------------------------------------------------------
!  Set tracer and/or mass point sources and/or sink.
!-----------------------------------------------------------------------
!
      IF (iic(ng).eq.ntstart) THEN
!
!  Set-up point Sources/Sink number (Nsrc), direction (Dsrc), I- and
!  J-grid locations (Isrc,Jsrc), and logical switch for type of tracer
!  to apply (Lsrc). Currently, the direction can be along XI-direction
!  (Dsrc = 0) or along ETA-direction (Dsrc > 0).  The mass sources are
!  located at U- or V-points so the grid locations should range from
!  1 =< Isrc =< L  and  1 =< Jsrc =< M.
!
        Lsrc(:,:)=.false.
# ifdef RIVERPLUME
        Nsrc=1
        Dsrc(Nsrc)=0.0_r8
        Isrc(Nsrc)=1
        Jsrc(Nsrc)=50
        Lsrc(Nsrc,itemp)=.true.
        Lsrc(Nsrc,isalt)=.true.
# elif defined SED_TEST1
        Nsrc=Mm(ng)*2
        DO is=1,Nsrc/2
          Dsrc(is)=0.0_r8
          Isrc(is)=1
          Jsrc(is)=is
        END DO
        DO is=Nsrc/2+1,Nsrc
          Dsrc(is)=0.0_r8
          Isrc(is)=Lm(ng)+1
          Jsrc(is)=is-Mm(ng)
        END DO
# else
        ANA_PSOURCE: no values provided for Dsrc, Isrc, Jsrc, Lsrc.
# endif
      END IF
# ifdef SOLVE3D
!
!  If appropriate, set-up nondimensional shape function to distribute
!  mass point sources/sinks vertically.  It most add to unity!!.
!
#  if defined SED_TEST1
        DO k=1,N(ng)
          DO is=1,Nsrc
            i=Isrc(is)
            j=Jsrc(is)
            Qshape(is,k)=ABS(u(i,j,k,nnew)/ubar(i,j,knew))*             &
     &                   (z_w(i-1,Mm(ng)/2,k)-z_w(i-1,Mm(ng)/2,k-1)+    &
     &                    z_w(i  ,Mm(ng)/2,k)-z_w(i  ,Mm(ng)/2,k-1))/   &
     &                   (z_w(i-1,Mm(ng)/2,N(ng))-z_w(i-1,Mm(ng)/2,0)+  &
     &                    z_w(i  ,Mm(ng)/2,N(ng))-z_w(i  ,Mm(ng)/2,0))
          END DO
        END DO
#  else
        DO k=1,N(ng)
          DO is=1,Nsrc
            Qshape(is,k)=1.0_r8/REAL(N(ng),r8)
          END DO
        END DO
#  endif
# endif
# ifdef UV_PSOURCE
!
!  Set-up vertically integrated mass transport (m3/s) of point
!  Sources/Sinks (positive in the positive U- or V-direction and
!  viceversa).
!
#  if defined RIVERPLUME
      IF ((tdays(ng)-dstart).lt.0.5_r8) THEN
        fac=1.0_r8+TANH((time(ng)-43200.0_r8)/43200.0_r8)
      ELSE
        fac=1.0_r8
      END IF
      DO is=1,Nsrc
        Qbar(is)=fac*1500.0_r8
      END DO
#  elif defined SED_TEST1
      my_area=0.0_r8         ! West end
      DO is=1,Nsrc/2
        i=Isrc(is)
        j=Jsrc(is)
        my_area=my_area+0.5_r8*(zeta(i-1,j,knew)+h(i-1,j)+              &
     &                          zeta(i  ,j,knew)+h(i  ,j))*on_u(i,j)
      END DO
      fac=-1000.0_r8*10.0_r8*1.0_r8
      DO is=1,Nsrc/2
        i=Isrc(is)
        j=Jsrc(is)
        Qbar(is)=fac*(0.5_r8*(zeta(i-1,j,knew)+h(i-1,j)+                &
     &                        zeta(i  ,j,knew)+h(i  ,j)))*              &
     &           on_u(i,j)/my_area
      END DO
      my_area=0.0_r8         ! East end
      DO is=Nsrc/2+1,Nsrc
        i=Isrc(is)
        j=Jsrc(is)
        my_area=my_area+0.5_r8*(zeta(i-1,j,knew)+h(i-1,j)+              &
     &                          zeta(i  ,j,knew)+h(i  ,j))*on_u(i,j)
      END DO
      fac=-1000.0_r8*10.0_r8*1.0_r8
      DO is=Nsrc/2+1,Nsrc
        i=Isrc(is)
        j=Jsrc(is)
        Qbar(is)=fac*(0.5_r8*(zeta(i-1,j,knew)+h(i-1,j)+                &
     &                        zeta(i  ,j,knew)+h(i  ,j)))*              &
     &           on_u(i,j)/my_area
      END DO
#  else
      ANA_PSOURCE: no values provided for Qbar.
#  endif
#  ifdef SOLVE3D
!
!  Set-up mass transport profile (m3/s) of point Sources/Sinks.
!
      DO k=1,N(ng)
        DO is=1,Nsrc
          Qsrc(is,k)=Qbar(is)*Qshape(is,k)
        END DO
      END DO
#  endif
# endif
# if defined TS_PSOURCE && defined SOLVE3D
!
!  Set-up tracer (tracer units) point Sources/Sinks.
!
#  ifdef RIVERPLUME
      DO k=1,N(ng)
        DO is=1,Nsrc
          Tsrc(is,k,itemp)=T0(ng)
          Tsrc(is,k,isalt)=0.0_r8
        END DO
      END DO
#  else
      ANA_PSOURCE: no values provided for Tsrc.
#  endif
# endif
      RETURN
      END SUBROUTINE ana_psource_grid
#endif

#if defined ANA_RAIN && defined BULK_FLUXES
      SUBROUTINE ana_rain (ng, tile)
!
!=======================================================================
!                                                                      !
!  This routine sets precipitation rate (kg/m2/s) using an             !
!  analytical expression.                                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_forces
!
      integer, intent(in) :: ng, tile

# include "tile.h"
!
      CALL ana_rain_tile (ng, Istr, Iend, Jstr, Jend,                   &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    FORCES(ng) % rain)
      RETURN
      END SUBROUTINE ana_rain
!
!***********************************************************************
      SUBROUTINE ana_rain_tile (ng, Istr, Iend, Jstr, Jend,             &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          rain)
!***********************************************************************
!
      USE mod_param

# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
!
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj

# ifdef ASSUMED_SHAPE
      real(r8), intent(out) :: rain(LBi:,LBj:)
# else
      real(r8), intent(out) :: rain(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set analytical precipitation rate (kg/m2/s).
!-----------------------------------------------------------------------
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          rain(i,j)=0.0_r8
        END DO
      END DO
# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, rain)
# endif
      RETURN
      END SUBROUTINE ana_rain_tile
#endif

#if defined ANA_SEDIMENT && (defined SEDIMENT || defined BBL_MODEL)
      SUBROUTINE ana_sediment (ng, tile)
!
!=======================================================================
!                                                                      !
!  This routine sets initial conditions for  sedimen t tracer fields   !
!  concentrations  (kg/m3) using analytical expressions for sediment   !
!  and/or bottom boundary layer configurations. It also sets initial   !
!  bed conditions in each sediment layer.                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_ocean
!
      integer, intent(in) :: ng, tile

# include "tile.h"
!
      CALL ana_sediment_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        GRID(ng) % pm,                            &
     &                        GRID(ng) % pn,                            &
# if defined BBL_MODEL && defined MB_BBL
     &                        OCEAN(ng) % rho,                          &
# endif
# ifdef SEDIMENT
     &                        OCEAN(ng) % t,                            &
     &                        OCEAN(ng) % bed,                          &
     &                        OCEAN(ng) % bed_frac,                     &
     &                        OCEAN(ng) % bed_mass,                     &
# endif
     &                        OCEAN(ng) % bottom)
      RETURN
      END SUBROUTINE ana_sediment
!
!***********************************************************************
      SUBROUTINE ana_sediment_tile (ng, Istr, Iend, Jstr, Jend,         &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              pm, pn,                             &
# if defined BBL_MODEL && defined MB_BBL
     &                              rho,                                &
# endif
# ifdef SEDIMENT
     &                              t,                                  &
     &                              bed, bed_frac, bed_mass,            &
# endif
     &                              bottom)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_sediment
!
# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      USE exchange_3d_mod, ONLY : exchange_r3d_tile
!
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj

# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
#  if defined BBL_MODEL && defined MB_BBL
      real(r8), intent(in) :: rho(LBi:,LBj:,:)
#  endif
#  ifdef SEDIMENT
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
      real(r8), intent(out) :: bed(LBi:,LBj:,:,:)
      real(r8), intent(out) :: bed_frac(LBi:,LBj:,:,:)
      real(r8), intent(out) :: bed_mass(LBi:,LBj:,:,:)
#  endif
      real(r8), intent(inout) :: bottom(LBi:,LBj:,:)
# else
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
#  if defined BBL_MODEL && defined MB_BBL
      real(r8), intent(in) :: rho(LBi:,LBj:,:)
#  endif
#  ifdef SEDIMENT
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(out) :: bed(LBi:UBi,LBj:UBj,Nbed,MBEDP)
      real(r8), intent(out) :: bed_frac(LBi:UBi,LBj:UBj,Nbed,NST)
      real(r8), intent(out) :: bed_mass(LBi:UBi,LBj:UBj,Nbed,NST)
#  endif
      real(r8), intent(inout) :: bottom(LBi:UBi,LBj:UBj,MBOTP)
# endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, ised, j, k
      real(r8) :: cff1, cff2, cff3, cff4

# include "set_bounds.h"

# if defined BBL_MODEL && !defined SEDIMENT
!
!-----------------------------------------------------------------------
!  If only bottom boundary layer and not sediment model, set bottom
!  sediment grain diameter (m) and density (kg/m3).
!-----------------------------------------------------------------------
!
#  if defined BL_TEST || defined NJ_BIGHT || defined CBLAST
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          bottom(i,j,isd50)=0.0005_r8
          bottom(i,j,idens)=2650.0_r8
        END DO
      END DO
#  elif defined LAKE_SIGNELL
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          bottom(i,j,isd50)=0.0005_r8
          bottom(i,j,idens)=2650.0_r8
        END DO
      END DO
#  else
      ANA_SEDIMENT: no values provided for bottom sediment properties.
#  endif
#  ifdef MB_BBL
#   undef YALIN
!
!-----------------------------------------------------------------------
!  If only Blass bottom boundary layer and not sediment model, set
!  set critical (threshold) bedload stress (m2/s2).
!-----------------------------------------------------------------------
!
#   ifdef YALIN

!  For more accurate estime of critical bedload stress, consider the
!  Yalin method (Miller et. al, 1977).
!
      Kvisc=0.0013_r8/rho0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          rhoWater=rho(i,j,1)+1000.0_r8
          cff=SQRT((bottom(i,j,idens)-rhoWater)*                        &
     &             g*bottom(i,j,isd50)*bottom(i,j,isd50)*               &
     &               bottom(i,j,isd50)/rhoWater)/Kvisc
!!        D=bottom(i,j,isd50)*g*                                        &
!!   &      ((bottom(i,j,idens)/rho0-1.0_r8)/Kvisk)**(1.0_r8/3.0_r8)
!!        theta_cr=0.3_r8./(1.0_r8+1.2_r8*D)+                           &
!!   &             0.055_r8*(1.0_r8-EXP(-0.02_r8*D))          
          IF (cff.lt.100.0_r8) THEN
            theta_cb=0.041_r8*(LOG(cff)**2)-0.356_r8*LOG(cff)-0.977_r8
!!          theta_cb=10**theta_cr
          ELSE IF (cff.gt.3000.0_r8) THEN
            theta_cb=0.045_r8
          ELSE
            theta_cb=0.132_r8*LOG(cff)-1.804_r8
!!          theta_cb=10.0_r8**theta_cr
          ENDIF
          bottom(i,j,itauc)=(bottom(i,j,idens)-rho0)*g*                 &
     &                       bottom(i,j,isd50)*theta_cb/rho0
        END DO
      END DO
#   else
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          bottom(i,j,itauc)=0.15_r8/rho0
        END DO
      END DO
#   endif
#  endif
!
!-----------------------------------------------------------------------
!  If only Blass bottom boundary layer and not sediment model, set
!  sediiment settling velocity (m/s).
!-----------------------------------------------------------------------
!
      Kvisc=0.0013_r8/rho0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          bottom(i,j,iwsed)=0.02_r8
!!
!! Consider Souslby (1997) estimate of settling velocity.
!!
!!        D=bottom(i,j,isd50)*g*                                        &
!!   &      ((bottom(i,j,idens)/rho0-1.0)/Kvisk)**(1.0_r8/3.0_r8)
!!        bottom(i,j,iwsed)=Kvisc*(SQRT(10.36_r8*10.36_r8+
!!   &                      1.049_r8*D*D*D)-10.36_r8)/bottom(i,j,isd50)
        END DO
      END DO
# endif

# ifdef SEDIMENT
!
!-----------------------------------------------------------------------
!  Initial sediment concentrations in the water column.
!-----------------------------------------------------------------------
!
      DO ised=1,NST
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              t(i,j,k,1,idsed(ised))=Csed(ised,ng)
              t(i,j,k,2,idsed(ised))=t(i,j,k,1,idsed(ised))
            END DO
          END DO
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Initial sediment bed layer properties.
!-----------------------------------------------------------------------
!
#  if defined LAKE_SIGNELL || defined ESTUARY_TEST
      DO k=1,Nbed
        DO j=JstrR,JendR
          DO i=IstrR,IendR
             bed(i,j,k,iaged)=time(ng)
             bed(i,j,k,ithck)=1.0_r8
             bed(i,j,k,iporo)=0.00_r8
             DO ised=1,NST
               bed_frac(i,j,k,ised)=1.0_r8/FLOAT(NST)
             ENDDO
!
!  Calculate mass so it is consistent with density, thickness, and
!  porosity.
!
             DO ised=1,NST
               bed_mass(i,j,k,ised)=bed(i,j,k,ithck)*                   &
    &                   Srho(ised,ng)*(1.0_r8-bed(i,j,k,iporo))*        &
    &                   bed_frac(i,j,k,ised)
             END DO
          END DO
        END DO
      END DO
!
!  Set exposed sediment layer properties.
!
      DO j=JstrR,JendR
        DO i=Istr,Iend
          cff1=1.0_r8
          cff2=1.0_r8
          cff3=1.0_r8
          cff4=1.0_r8
          DO ised=1,NST
            cff1=cff1*Sd50(ised,ng)**bed_frac(i,j,1,ised)
            cff2=cff2*Srho(ised,ng)**bed_frac(i,j,1,ised)
            cff3=cff3*wsed(ised,ng)**bed_frac(i,j,1,ised)
            cff3=cff4*tau_ce(ised,ng)**bed_frac(i,j,1,ised)
          END DO
          bottom(i,j,isd50)=cff1
          bottom(i,j,idens)=cff2
          bottom(i,j,iwsed)=cff3
          bottom(i,j,itauc)=cff4
          bottom(i,j,izdef)=Zob(ng)
        END DO
      END DO
#  else
      ANA_SEDIMENT: no values provided for sediment fields.
#  endif
# endif
# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Exchange boundary data.
!-----------------------------------------------------------------------
!
#  ifdef SEDIMENT
      DO ised=1,NST
        CALL exchange_r3d_tile (ng, Istr, Iend, Jstr, Jend,             &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          t(:,:,:,1,idsed(ised)))
        CALL exchange_r3d_tile (ng, Istr, Iend, Jstr, Jend,             &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          t(:,:,:,2,idsed(ised)))
        CALL exchange_r3d_tile (ng, Istr, Iend, Jstr, Jend,             &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          bed_frac(:,:,:,ised))
        CALL exchange_r3d_tile (ng, Istr, Iend, Jstr, Jend,             &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          bed_mass(:,:,:,ised))
      END DO
      CALL exchange_r3d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, 1, MBEDP,             &
     &                        bed)
#  endif
      CALL exchange_r3d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, 1, MBOTP,             &
     &                        bottom)
# endif

      RETURN
      END SUBROUTINE ana_sediment_tile
#endif

#ifdef ANA_SMFLUX
      SUBROUTINE ana_smflux (ng, tile)
!
!=======================================================================
!                                                                      !
!  This routine sets kinematic surface momentum flux (wind stress)     !
!  "sustr" and "svstr" (m2/s2) using an analytical expression.         !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_forces
      USE mod_grid
!
      integer, intent(in) :: ng, tile

# include "tile.h"
!
      CALL ana_smflux_tile (ng, Istr, Iend, Jstr, Jend,                 &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      GRID(ng) % angler,                          &
# ifdef SPHERICAL
     &                      GRID(ng) % lonr,                            &
     &                      GRID(ng) % latr,                            &
# else
     &                      GRID(ng) % xr,                              &
     &                      GRID(ng) % yr,                              &
# endif
     &                      FORCES(ng) % sustr,                         &
     &                      FORCES(ng) % svstr)
      RETURN
      END SUBROUTINE ana_smflux
!
!***********************************************************************
      SUBROUTINE ana_smflux_tile (ng, Istr, Iend, Jstr, Jend,           &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            angler,                               &
# ifdef SPHERICAL
     &                            lonr, latr,                           &
# else
     &                            xr, yr,                               &
# endif
     &                            sustr, svstr)
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

# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: angler(LBi:,LBj:)
#  ifdef SPHERICAL
      real(r8), intent(in) :: lonr(LBi:,LBj:)
      real(r8), intent(in) :: latr(LBi:,LBj:)
#  else
      real(r8), intent(in) :: xr(LBi:,LBj:)
      real(r8), intent(in) :: yr(LBi:,LBj:)
#  endif
      real(r8), intent(out) :: sustr(LBi:,LBj:)
      real(r8), intent(out) :: svstr(LBi:,LBj:)
# else
      real(r8), intent(in) :: angler(LBi:UBi,LBj:UBj)
#  ifdef SPHERICAL
      real(r8), intent(in) :: lonr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: latr(LBi:UBi,LBj:UBj)
#  else
      real(r8), intent(in) :: xr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: yr(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(out) :: sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: svstr(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j
      real(r8) :: Ewind, Nwind, val1, val2, windamp, winddir
# ifdef LAKE_SIGNELL
      real(r8) :: mxst, mxtm, rmptm
# endif

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set kinematic surface momentum flux (wind stress) component in the
!  XI-direction (m2/s2) at horizontal U-points.
!-----------------------------------------------------------------------
!
# ifdef BASIN
      val1=5.0E-05_r8*(1.0_r8+TANH((time(ng)-6.0_r8*86400.0_r8)/        &
     &                 (3.0_r8*86400.0_r8)))
      val2=2.0_r8*pi/el(ng)
      DO j=JstrR,JendR
        DO i=Istr,IendR
          sustr(i,j)=-val1*COS(val2*yr(i,j))
        END DO
      END DO
# elif defined BL_TEST
      Ewind=0.0_r8/rho0
      Nwind=0.3_r8/rho0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          sustr(i,j)=Ewind
        END DO
      END DO
# elif defined CANYON_A || defined CANYON_B
      DO j=JstrR,JendR
        DO i=Istr,IendR
          sustr(i,j)=5.0E-05_r8*SIN(2.0_r8*pi*tdays(ng)/10.0_r8)*       &
     &               (1.0_r8-TANH((yr(i,j)-0.5_r8*el(ng))/10000.0_r8))
        END DO
      END DO
# elif defined DOUBLE_GYRE
      IF ((tdays(ng)-dstart).le.2.0_r8) THEN
!!      windamp=-0.05_r8*SIN(pi*tdays(ng)/4.0_r8)/rho0
        windamp=user(1)*SIN(pi*tdays(ng)/4.0_r8)/rho0
      ELSE
!!      windamp=-0.05_r8/rho0
        windamp=user(1)/rho0
      END IF
      val1=2.0_r8*pi/el(ng)
      DO j=JstrR,JendR
        DO i=Istr,IendR
          sustr(i,j)=windamp*COS(val1*yr(i,j))
        END DO
      END DO
# elif defined FLT_TEST
      DO j=JstrR,JendR
        DO i=Istr,IendR
          sustr(i,j)=1.0E-03_r8
        END DO
      END DO
# elif defined LAKE_SIGNELL
      mxst=0.2515_r8          ! N/m2
      mxtm=48.0_r8            ! hours
      rmptm=3.0_r8            ! hours
      DO j=JstrR,JendR
         DO i=Istr,IendR
!!         sustr(i,j)=0.0_r8
           sustr(i,j)=MIN(mxst/rho0*                                    &
     &                    (TANH((time(ng)/60.0_r8-rmptm*150.0_r8)/      &
     &                           150.0_r8)+1.0_r8)*0.5_r8,              &
     &                    mxst/rho0*                                    &
     &                    (TANH(((mxtm*3600.0_r8-time(ng))/             &
     &                           60.0_r8-rmptm*150.0_r8)/150.0_r8)+     &
     &                     1.0_r8)*0.5_r8)
         END DO
      END DO
# elif defined LMD_TEST
      IF (time(ng).le.57600.0_r8) THEN
        windamp=-0.6_r8*SIN(pi*time(ng)/57600.0_r8)*                    &
     &                  SIN(2.0_r8*pi/57600.0_r8)/rho0
      ELSE
        windamp=0.0_r8
      END IF
      DO j=JstrR,JendR
        DO i=Istr,IendR
          sustr(i,j)=windamp
        END DO
      END DO
# elif defined NJ_BIGHT
!!    windamp=0.086824313_r8
!!    winddir=0.5714286_r8
!!    if ((tdays(ng)-dstart).le.0.5_r8) then
!!      Ewind=windamp*winddir*SIN(pi*(tdays(ng)-dstart))/rho0
!!      Nwind=windamp*SIN(pi*(tdays(ng)-dstart))/rho0
!!    else
!!      Ewind=windamp*winddir/rho0
!!      Nwind=windamp/rho0
!!    endif
      IF ((tdays(ng)-dstart).le.3.0_r8) THEN
         winddir=60.0_r8
         windamp=0.1_r8
      ELSE IF (((tdays(ng)-dstart).gt.3.0_r8).and.                      &
     &        ((tdays(ng)-dstart).le.4.0_r8)) THEN
         winddir= 60.0_r8*((tdays(ng)-dstart)-2.0_r8)-                  &
     &           120.0_r8*((tdays(ng)-dstart)-2.0_r8)
         windamp=0.0_r8
      ELSE
         winddir=-120.0_r8
         windamp=0.0_r8
      END IF
      Ewind=windamp*COS(pi*winddir/180.0_r8)/rho0
      Nwind=windamp*SIN(pi*winddir/180.0_r8)/rho0
      DO j=JstrR,JendR
        DO i=Istr,IendR
          val1=0.5_r8*(angler(i-1,j)+angler(i,j))
          sustr(i,j)=Ewind*COS(val1)+Nwind*SIN(val1)
        END DO
      END DO
# elif defined UPWELLING
      IF ((tdays(ng)-dstart).le.2.0_r8) THEN
        windamp=-0.1_r8*SIN(pi*(tdays(ng)-dstart)/4.0_r8)/rho0
      ELSE
        windamp=-0.1_r8/rho0
      END IF
      DO j=JstrR,JendR
        DO i=Istr,IendR
          sustr(i,j)=windamp
        END DO
      END DO
# elif defined USWEST
      DO j=JstrR,JendR
        DO i=Istr,IendR
!!        val1=(latr(i,j)-latr(Lm(ng)/2,Mm(ng)/2))/20.0_r8
!!        sustr(i,j)=1.0E-04_r8*val1
!!        sustr(i,j)=-1.0E-04_r8
          sustr(i,j)=0.0_r8
        END DO
      END DO
# elif defined WINDBASIN
      IF ((tdays(ng)-dstart).le.2.0_r8) THEN
        windamp=-0.1_r8*SIN(pi*(tdays(ng)-dstart)/4.0_r8)/rho0
      ELSE
        windamp=-0.1_r8/rho0
      END IF
      DO j=JstrR,JendR
        DO i=Istr,IendR
          sustr(i,j)=windamp
        END DO
      END DO
# else
      DO j=JstrR,JendR
        DO i=Istr,IendR
          sustr(i,j)=0.0_r8
        END DO
      END DO
# endif
!
!-----------------------------------------------------------------------
!  Set kinematic surface momentum flux (wind stress) component in the
!  ETA-direction (m2/s2) at horizontal V-points.
!-----------------------------------------------------------------------
!
# if defined BL_TEST
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          svstr(i,j)=Nwind
        END DO
      END DO
# elif defined LMD_TEST
      IF (time(ng).le.57600.0_r8) THEN
        windamp=-0.6_r8*SIN(pi*time(ng)/57600.0_r8)*                    &
     &                  COS(2.0_r8*pi/57600.0_r8)/rho0
      ELSE
        windamp=0.0_r8
      END IF
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          svstr(i,j)=windamp
        END DO
      END DO
# elif defined NJ_BIGHT
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          val1=0.5_r8*(angler(i,j)+angler(i,j-1))
          svstr(i,j)=-Ewind*SIN(val1)+Nwind*COS(val1)
        END DO
      END DO
# elif defined USWEST
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          svstr(i,j)=-1.0E-04_r8
        END DO
      END DO
# else
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          svstr(i,j)=0.0_r8
        END DO
      END DO
# endif
# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      CALL exchange_u2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, sustr)
      CALL exchange_v2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, svstr)
# endif
      RETURN
      END SUBROUTINE ana_smflux_tile
#endif

#ifdef ANA_SPINNING
      SUBROUTINE ana_spinning (ng, tile)
!
!=======================================================================
!                                                                      !
!  This subroutine sets time-variable rotation force as the sum of     !
!  Coriolis and Centripetal accelerations.  This is used in polar      !
!  coordinate applications (annulus grid).                             !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
!
      integer, intent(in) :: ng, tile

# include "tile.h"
!
      CALL ana_spinning_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
# ifdef SPHERICAL
     &                        GRID(ng) % lonr,                          &
     &                        GRID(ng) % latr,                          &
# else
     &                        GRID(ng) % xr,                            &
     &                        GRID(ng) % yr,                            &
# endif
     &                        GRID(ng) % f,                             &
     &                        GRID(ng) % omn,                           &
     &                        GRID(ng) % fomn)
      RETURN
      END SUBROUTINE ana_spinning
!
!***********************************************************************
      SUBROUTINE ana_spinning_tile (ng, Istr, Iend, Jstr, Jend,         &
     &                              LBi, UBi, LBj, UBj,                 &
# ifdef SPHERICAL
     &                              lonr, latr                          &
# else
     &                              xr, yr,                             &
# endif
     &                              f, omn, fomn)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj

# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: f(LBi:,LBj:)
      real(r8), intent(in) :: omn(LBi:,LBj:)
#  ifdef SPHERICAL
      real(r8), intent(in) :: lonr(LBi:,LBj:)
      real(r8), intent(in) :: latr(LBi:,LBj:)
#  else
      real(r8), intent(in) :: xr(LBi:,LBj:)
      real(r8), intent(in) :: yr(LBi:,LBj:)
#  endif
      real(r8), intent(out) :: fomn(LBi:,LBj:)
# else
      real(r8), intent(in) :: f(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: omn(LBi:UBi,LBj:UBj)
#  ifdef SPHERICAL
      real(r8), intent(in) :: lonr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: latr(LBi:UBi,LBj:UBj)
#  else
      real(r8), intent(in) :: xr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: yr(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(out) :: fomn(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV

# ifdef LAB_CANYON
      real(r8), parameter :: Omega0 = 2.0_r8*pi/25.0_r8
      real(r8), parameter :: Width = 0.20_r8
      real(r8), parameter :: Ro = 0.10_r8
      real(r8), parameter :: Rs = 0.55_r8
      real(r8), parameter :: little_omega = 2.0_r8*pi/24.0_r8
      real(r8), parameter :: Bu = 10.0_r8
      real(r8), parameter :: hd = 0.125_r8

      real(r8) :: Omega1, Omega1_of_t, Ro_t
      real(r8) :: fcor, d_rho_dz, d_Omega1_dt, time_fac
# endif

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute time-varying rotation force: Coriolis plus Centripetal
!  accelerations.
!-----------------------------------------------------------------------
!
# ifdef LAB_CANYON
      fcor=2.0_r8*Omega0
      Omega1=fcor*Width*Ro/Rs
      Ro_t=little_omega/fcor
      d_rho_dz=(1000.0_r8*Bu/g)*(fcor*Width/hd)**2
      time_fac=1.0_r8+(Omega1/Omega0)*SIN(little_omega*time(ng))
      Omega1_of_t=Omega1*SIN(little_omega*time(ng))
      d_Omega1_dt=Omega1*little_omega*COS(little_omega*time(ng))
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          fomn(i,j)=(f(i,j)*time_fac+                                   &
     &               SQRT(xr(i,j)*xr(i,j)+yr(i,j)*yr(i,j))*             &
     &               ((2.0_r8*Omega0+Omega1_of_t)*Omega1_of_t))*        &
     &              omn(i,j)
        END DO
      END DO
# endif
      RETURN
      END SUBROUTINE ana_spinning_tile
#endif

#if (defined ANA_SRFLUX && defined SOLVE3D) || defined DIURNAL_SRFLUX 
      SUBROUTINE ana_srflux (ng, tile)
!
!=======================================================================
!                                                                      !
!  This subroutine sets kinematic surface solar shortwave radiation    !
!  flux "srflx" (degC m/s) using an analytical expression.             !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_forces
      USE mod_grid
!
      integer, intent(in) :: ng, tile

# include "tile.h"
!
      CALL ana_srflux_tile (ng, Istr, Iend, Jstr, Jend,                 &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      GRID(ng) % lonr,                            &
     &                      GRID(ng) % latr,                            &
# ifdef ALBEDO
     &                      FORCES(ng) % cloud,                         &
     &                      FORCES(ng) % Hair,                          &
     &                      FORCES(ng) % Tair,                          &
# endif
     &                      FORCES(ng) % srflx)
      RETURN
      END SUBROUTINE ana_srflux
!
!***********************************************************************
      SUBROUTINE ana_srflux_tile (ng, Istr, Iend, Jstr, Jend,           &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            lonr, latr,                           &
# ifdef ALBEDO 
     &                            cloud, Hair, Tair,                    &
# endif
     &                            srflx)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars

# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
!
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj

# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: lonr(LBi:,LBj:)
      real(r8), intent(in) :: latr(LBi:,LBj:)
#  ifdef ALBEDO
      real(r8), intent(in) :: cloud(LBi:,LBj:)
      real(r8), intent(in) :: Hair(LBi:,LBj:)
      real(r8), intent(in) :: Tair(LBi:,LBj:)
#  endif
      real(r8), intent(out) :: srflx(LBi:,LBj:)
# else
      real(r8), intent(in) :: lonr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: latr(LBi:UBi,LBj:UBj)
#  ifdef ALBEDO
      real(r8), intent(in) :: cloud(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Hair(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Tair(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(out) :: srflx(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j
# if defined ALBEDO || defined DIURNAL_SRFLUX
      integer :: iday, month, year
      real(r8) :: Dangle, Hangle, LatRad
      real(r8) :: cff, cff1, cff2, hour, yday
#  ifdef ALBEDO
      real(r8) :: Rsolar, e_sat, vap_p, zenith
#  endif
# endif

# include "set_bounds.h"

# if defined ALBEDO || defined DIURNAL_SRFLUX
!
!-----------------------------------------------------------------------
!  Compute shortwave radiation (degC m/s):
!
!  ALBEDO option: Compute shortwave radiation flux using the Laevastu
!                 cloud correction to the Zillman equation for cloudless
!  radiation (Parkinson and Washington 1979, JGR, 84, 311-337).  Notice
!  that flux is scaled from W/m2 to degC m/s by dividing by (rho0*Cp).
!
!  DIURNAL_SRFLUX option: Modulate shortwave radiation SRFLX (which
!                         read and interpolated elsewhere) by the local
!  diurnal cycle (a function of longitude, latitude and day-of-year).
!  This option is provided for cases where SRFLX computed by SET_DATA is
!  an average over >= 24 hours.  For DIURNAL_SRFLUX to work ANA_SRFLUX
!  must be undefined. If you want a strictly analytical diurnal cycle
!  enter it explicitly at the end of this subroutine or use the ALBEDO
!  option.
!-----------------------------------------------------------------------
!
!  Assume time is in modified Julian day.  Get hour and year day.
!
      CALL caldate (r_date, tdays(ng), year, yday, month, iday, hour)
!
!  Estimate solar declination angle (radians).
!
      Dangle=23.44_r8*COS((172.0_r8-yday)*2.0_r8*pi/365.25_r8)
      Dangle=Dangle*deg2rad
!
!  Compute hour angle (radians).
!
      Hangle=(12.0_r8-hour)*pi/12.0_r8
!
#  ifdef ALBEDO
      Rsolar=Csolar/(rho0*Cp)
#  endif
      DO j=JstrR,JendR
        DO i=IstrR,IendR
!
!  Local daylight is a function of the declination (Dangle) and hour 
!  angle adjusted for the local meridian (Hangle-lonr(i,j)*deg2rad). 
!
          LatRad=latr(i,j)*deg2rad
          cff1=SIN(LatRad)*SIN(Dangle)
          cff2=COS(LatRad)*COS(Dangle)
# if defined ALBEDO
!
!  Estimate variation in optical thickness of the atmosphere over
!  the course of a day under cloudless skies. To obtain net incoming 
!  shortwave radiation multiply by (1.0-0.6*c**3), where c is the 
!  fractional cloud cover.
!
          srflx(i,j)=0.0_r8
          zenith=cff1+cff2*COS(Hangle-lonr(i,j)*deg2rad)
          IF (zenith.gt.0.0_r8) THEN
            cff=(0.7859_r8+0.03477_r8*Tair(i,j))/                       &
     &          (1.0_r8+0.00412_r8*Tair(i,j))
            e_sat=10.0_r8**cff
            vap_p=e_sat*Hair(i,j)
            srflx(i,j)=Rsolar*zenith*zenith*                            &
     &                 (1.0_r8-0.6_r8*cloud(i,j)**3)/                   &
     &                 ((zenith+2.7_r8)*vap_p*1.0E-5_r8+                &
     &                  1.085_r8*zenith+0.1_r8)
          END IF
# elif defined DIURNAL_SRFLUX
!
!  SRFLX is reset on each time step in subroutine SET_DATA which 
!  interpolates values in the forcing file to the current date.
!  This DIURNAL_SRFLUX option is provided so that SRFLX values
!  corresponding to a greater or equal daily average can be modulated
!  by the local length of day to produce a diurnal cycle with the 
!  same daily average as the original data.  This approach assumes 
!  the net effect of clouds is incorporated into the SRFLX data. 
!
!  Normalization factor = INTEGRAL{ABS(a+b*COS(t)) dt} from 0 to 2*pi 
!                       = (a*ARCCOS(-a/b)+SQRT(b**2-a**2))/pi
!  
          cff=(cff1*ACOS(-cff1/cff2)+SQRT(cff2*cff2-cff1*cff1))/pi
          srflx(i,j)=MAX(0.0_r8,                                        &
     &                   srflx(i,j)/cff*                                &
     &                   (cff1+cff2*COS(Hangle-lonr(i,j)*deg2rad)))
# endif
        END DO
      END DO
# else
!
!-----------------------------------------------------------------------
!  Set incoming solar shortwave radiation (W/m2).
!-----------------------------------------------------------------------
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          srflx(i,j)=0.0_r8
        END DO
      END DO
# endif
# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, srflx)
# endif
      RETURN
      END SUBROUTINE ana_srflux_tile
#endif

#if defined ECOSIM && defined SOLVE3D
      SUBROUTINE ana_specir (ng, tile)
!
!=======================================================================
!                                                                      !
!  This subroutine sets surface solar downwelling spectral irradiance  !
!  at just beneath the sea surface, Ed(lambda,0-), in micromol quanta  !
!  per meter squared per second.                                       !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    Gregg, W.W. and K.L. Carder, 1990:  A simple spectral solar       !
!           irradiance model for cloudless maritime atmospheres,       !
!           Limmol. Oceanogr., 35(8), 1657-1675.                       !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_forces
      USE mod_grid
!
      implicit none

      integer, intent(in) :: ng, tile

# include "tile.h"
!
      CALL ana_specir_tile (ng, Istr, Iend, Jstr, Jend,                 &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      GRID(ng) % lonr,                            &
     &                      GRID(ng) % latr,                            &
     &                      FORCES(ng) % cloud,                         &
     &                      FORCES(ng) % Hair,                          &
     &                      FORCES(ng) % Tair,                          &
     &                      FORCES(ng) % Pair,                          &
     &                      FORCES(ng) % Uwind,                         &
     &                      FORCES(ng) % Vwind,                         &
     &                      FORCES(ng) % SpecIr,                        &
     &                      FORCES(ng) % avcos)
      RETURN
      END SUBROUTINE ana_specir
!
!***********************************************************************
      SUBROUTINE ana_specir_tile (ng, Istr, Iend, Jstr, Jend,           &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            lonr, latr,                           &
     &                            cloud, Hair, Tair, Pair,              &
     &                            Uwind, Vwind,                         &
     &                            SpecIr, avcos)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_eclight

# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
!
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
# endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: lonr(LBi:,LBj:)
      real(r8), intent(in) :: latr(LBi:,LBj:)
      real(r8), intent(in) :: cloud(LBi:,LBj:)
      real(r8), intent(in) :: Hair(LBi:,LBj:)
      real(r8), intent(in) :: Tair(LBi:,LBj:)
      real(r8), intent(in) :: Pair(LBi:,LBj:)
      real(r8), intent(in) :: Uwind(LBi:,LBj:)
      real(r8), intent(in) :: Vwind(LBi:,LBj:)
      real(r8), intent(out) :: SpecIr(LBi:,LBj:,:)
      real(r8), intent(out) :: avcos(LBi:,LBj:,:)
# else
      real(r8), intent(in) :: lonr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: latr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: cloud(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Hair(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Tair(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Pair(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: UWind(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Vwind(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: SpecIr(LBi:UBi,LBj:UBj,NBands)
      real(r8), intent(out) :: avcos(LBi:UBi,LBj:UBj,NBands)
# endif
!
!  Local constant declarations.
!
      real(r8) :: am = 1.0_r8        ! Aerosol type 1-10: ocean to land
      real(r8) :: betalam = 0.55_r8
      real(r8) :: p0 = 29.92_r8      ! Standard pressure (inches of Hg)
      real(r8) :: rex = -1.6364_r8
      real(r8) :: roair = 1200.0_r8  ! Density of air (g/m3)
      real(r8) :: rn = 1.341_r8      ! Index of refraction of pure seawater
      real(r8) :: vis = 15.0_r8      ! Visibility (km)
      real(r8) :: wv = 1.5_r8        ! Precipitable water (cm): water vapor
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, iband, ic, j, nc
      integer :: iday, month, year

      real(r8) :: Dangle, Hangle, LatRad, LonRad
      real(r8) :: cff, cff1, cff2, hour, yday
      real(r8) :: alpha, beta, gamma, theta, rtheta, rthetar
      real(r8) :: atra, gtra, otra, rtra, wtra
      real(r8) :: alg, arg, asymp, cosunz, Fa
      real(r8) :: frh, rh, rlam, rlogc
      real(r8) :: rm, rmin, rmo, rmp, rod, rof
      real(r8) :: ros, rospd, rosps, rpls 
      real(r8) :: sumx, sumx2, sumxy, sumy
      real(r8) :: taa, tas, to3, wa, wspeed, zenith

      real(r8), dimension(NBands) :: Fo, Edir, Edif, Ed, qlam

      real(r8), dimension(3) :: a_arr, dndr
      real(r8), dimension(3) :: ro     = (/ 0.03_r8, 0.24_r8, 2.00_r8 /)
      real(r8), dimension(3) :: r_arr  = (/ 0.10_r8, 1.00_r8, 10.0_r8 /)
!
# include "set_bounds.h"

!
!-----------------------------------------------------------------------
!  Compute spectral irradiance: Using RADTRAN formulations.
!-----------------------------------------------------------------------
!
!  Assume time is in modified Julian day.  Get hour and year day.
!
      CALL caldate (r_date, tdays(ng), year, yday, month, iday, hour)
!
!  Estimate solar declination angle (radians).
!
      Dangle=23.44_r8*COS((172.0_r8-yday)*2.0_r8*pi/365.25_r8)
      Dangle=Dangle*deg2rad
!
!  Compute hour angle (radians).
!
      Hangle=(12.0_r8-hour)*pi/12.0_r8
!
!  Conversion constant from E to micromol quanta.
!   1/(Plank*c*(avos#*1.0e6).
!
      cff=1.0E-9_r8/(6.6256E-34_r8*2.998E8_r8*6.023E17_r8)
      DO iband=1,NBands
        qlam(iband)=ec_wave_ab(iband)*cff
      END DO
!
!  Correct solar constant for Earth-Sun distance
!
      cff=(1.0_r8+0.0167_r8*COS(2.0_r8*pi*(yday-3.0_r8)/365.0_r8))**2
      DO iband=1,NBands
        Fo(iband)=ec_Fobar(iband)*cff
      END DO
!
!  Compute spectral irradiance.
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR

          LatRad=latr(i,j)*deg2rad
          LonRad=lonr(i,j)*deg2rad
!
!  Compute Climatological Ozone.
!
          to3=235.0_r8+(150.0_r8+40.0_r8*                               &
     &        SIN(0.9865_r8*(yday-30.0_r8)*deg2rad)+                    &
     &        20.0_r8*SIN(3.0_r8*LonRad))*                              &
     &        SIN(1.28_r8*LatRad)*SIN(1.28_r8*LatRad)*                  &
     &        0.001_r8                                 ! sco3 conversion
!
!  Local daylight is a function of the declination (Dangle) and hour 
!  angle adjusted for the local meridian (Hangle-lonr(i,j)*deg2rad). 
!
          zenith=SIN(LatRad)*SIN(Dangle)+                               &
     &           COS(LatRad)*COS(Dangle)*COS(Hangle-lonr(i,j)*deg2rad)
          theta=zenith*rad2deg
          cosunz=COS(zenith)
!
!  Model for atmospheric transmittance of solar irradiance through
!  a maritime atmosphere.  Computes direct and diffuse separately.
!  Includes water vapor and oxygen absorption.
!
!  Compute atmospheric path lengths (air mass); pressure-corrected
!
          IF ((theta.ge.0.0_r8).and.(theta.le.90.0_r8)) THEN
!
!  Modified March, 1994 according to Kasten and Young 1989.
!
            rm=1.0_r8/(cosunz+0.50572_r8*(96.07995_r8-theta)**rex)
            rmp=rm*(Pair(i,j)*0.02952756_r8)/p0
            rmo=(1.0_r8+22.0_r8/6370.0_r8)/                             &
     &          SQRT(cosunz*cosunz+44.0_r8/6370.0_r8)
!
!  Computes aerosol parameters according to a simplified version
!  of the Navy marine aerosol model.
!
!  Compute wind speed (24 hour mean is equal to current wind).
!
            wspeed=SQRT(Uwind(i,j)*Uwind(i,j)+Vwind(i,j)*Vwind(i,j))
!
!  Relative humidity factor, frh.
!
            rh=Hair(i,j)
            IF (rh.ge.100.0_r8) rh=99.9_r8
            frh=((2.0_r8-rh*0.01_r8)/                                   &
                 (6.0_r8*(1.0_r8-rh*0.01_r8)))**0.333_r8
!
!  Size distribution amplitude components.
!
            a_arr(1)=2000.0_r8*am*am
            a_arr(2)=5.866_r8*(wspeed-2.2_r8)
            a_arr(3)=0.01527_r8*(wspeed-2.2_r8)*0.05_r8   !from Hughes 1987
            IF (a_arr(2).lt.0.5_r8) a_arr(2)=0.5_r8
            IF (a_arr(3).lt.0.000014_r8) a_arr(3)=0.000014_r8
!
!  Compute size distribution at three selected radii according to
!  Navy method.
!
            cff=1.0_r8/frh
            DO nc=1,3
              dndr(nc)=0.0_r8
              DO ic=1,3
                arg=LOG(r_arr(nc)/(frh*ro(ic)))
                dndr(nc)=dndr(nc)+a_arr(ic)*EXP(-arg*arg)*cff
              END DO
            END DO
!
!  Least squares approximation
!
            sumx=0.0_r8
            sumy=0.0_r8
            sumxy=0.0_r8
            sumx2=0.0_r8
            DO ic=1,3
              cff1=LOG10(r_arr(ic))
              cff2=LOG10(dndr(ic))
              sumx=sumx+cff1 
              sumy=sumy+cff2 
              sumxy=sumxy+cff1*cff2
              sumx2=sumx2+cff1*cff1
            END DO
            gamma=sumxy/sumx2
            rlogc=sumy/3.0_r8-gamma*sumx/3.0_r8          ! no used
            alpha=-(gamma+3.0_r8)
!
!  Compute aerosol turbity coefficient, beta.
!
            beta=(3.91_r8/vis)*betalam**alpha
!
!  Compute asymmetry parameter -- a function of alpha.
!
            IF (alpha.gt.1.2_r8) THEN
              asymp=0.65_r8
            ELSE IF (alpha .lt. 0.0_r8) THEN
              asymp=0.82_r8
            ELSE
              asymp=-0.14167_r8*alpha+0.82_r8
            END IF
!
!  Single scattering albedo at 550; function of RH.
!
            wa=(-0.0032_r8*am+0.972_r8)*EXP(0.000306_r8*rh)
!
!  Forward scattering probability.
!
            alg=LOG(1.0_r8-asymp)
            Fa=1.0_r8-0.5_r8*                                           &
     &         EXP((alg*(1.459_r8+alg*(0.1595_r8+alg*0.4129_r8))+       &
     &              alg*(0.0783_r8+alg*(-0.3824_r8-alg*0.5874_r8))*     &
     &             cosunz)*cosunz)
!
!  Compute surface reflectance for direct (rod) and diffuse (ros)
!  components separately, as a function of theta, wind speed or
!  stress.
!
!  Foam and diffuse reflectance
!
            IF (wspeed.gt.4.0_r8) THEN
              IF (wspeed.le.7.0_r8) THEN
                rof=roair*(0.00062_r8+0.00156_r8/wspeed)*               &
     &              0.000022_r8*wspeed*wspeed-0.00040_r8
              ELSE
                rof=(roair*(0.00049_r8+0.000065_r8*wspeed)*             &
     &               0.000045_r8-0.000040_r8)*wspeed*wspeed
              END IF
              rosps=0.057_r8
            ELSE
              rof=0.0_r8
              rosps=0.066_r8
            END IF
!
!  Direct Fresnel reflectance for theta < 40, wspeed < 2 m/s.
!
            IF ((theta.lt.40.0_r8).or.(wspeed.lt.2.0_r8)) THEN
              IF (theta.eq.0.0_r8) THEN
                rospd=0.0211_r8
              ELSE
                rtheta=zenith
                rthetar=ASIN(SIN(rtheta)/rn)
                rmin=rtheta-rthetar
                rpls=rtheta+rthetar
                rospd=0.5_r8*((SIN(rmin)*SIN(rmin))/                    &
     &                        (SIN(rpls)*SIN(rpls))+                    &
     &                        (TAN(rmin)*TAN(rmin))/                    &
     &                        (TAN(rpls)*TAN(rpls)))
              END IF
!
!  Empirical fit otherwise.
!
            ELSE
              rospd=0.0253_r8*EXP((-0.000714_r8*wspeed+0.0618_r8)*      &
     &                            (theta-40.0_r8))
            END IF
!
!  Reflectance totals.
!
            rod=rospd+rof
            ros=rosps+rof
!
!  Compute spectral irradiance for selected spectral bands.
!
            DO iband=1,NBands
              rlam=ec_wave_ab(iband)*0.001_r8
!
!  Transmittance, Rayleigh, by Bird method.
!
              rtra=EXP(-rmp/(115.6406_r8*rlam**4-1.335_r8*rlam**2))
!
!  Ozone.
!
              otra=EXP(-ec_aoz(iband)*to3*rmo)
!
!  Aerosols.
!
              arg=beta*rm*rlam**(-alpha)
              atra=EXP(-arg)
              taa=EXP(-(1.0_r8-wa)*arg)
              tas=EXP(-wa*arg)
!
!  Oxygen/gases.
!
              gtra=EXP((-1.41_r8*ec_ag(iband)*rmp)/                     &
     &                 ((1.0_r8+118.3_r8*ec_ag(iband)*rmp)**0.45_r8))
!
!  Water Vapor.
!
              wtra=EXP((-0.2385_r8*ec_aw(iband)*wv*rm)/                 &
     &                 ((1.0_r8+20.07_r8*ec_aw(iband)*wv*rm)**0.45_r8))
!
!  Direct irradiance.
!
              Edir(iband)=Fo(iband)*cosunz*rtra*otra*atra*gtra*         &
     &                    wtra*(1.0_r8-rod)
!
!  Total diffuse irradiance.
!
              Edif(iband)=(1.0_r8-ros)*                                 &
     &                    Fo(iband)*cosunz*gtra*wtra*otra*              &
     &                    (taa*0.5_r8*(1.0_r8-rtra**0.95_r8)+           &
     &                     taa*Fa*(1.0_r8-tas)*rtra**1.5_r8)
!
!  Cloud effects approximations, Kasten and Czeplak (1980).
!  (See Hydrolight Technical Notes).
!
              IF (cloud(i,j).gt.0.25_r8) THEN
                Ed(iband)=(Edir(iband)+Edif(iband))*                    &
     &                    (1.0_r8-0.75_r8*cloud(i,j)**3.4_r8)
                Edif(iband)=Ed(iband)*                                  &
     &                      (0.3_r8+0.7_r8*cloud(i,j)**2.0_r8)
              ELSE
                Ed(iband)=Edir(iband)+Edif(iband)
              END IF
!
!  Convert to micromole quanta.
!
              SpecIr(i,j,iband)=Ed(iband)*qlam(iband)
!
!  Correction of zenith angle after crossing air/sea interface.
!
              cff1=COS(ASIN((SIN(zenith))/rn))
              cff2=Edif(iband)/Ed(iband)
              avcos(i,j,iband)=cff1*(1.0_r8-cff2)+0.86_r8*cff2
            END DO
          ELSE
            DO iband=1,NBands
              SpecIr(i,j,iband)=0.0_r8
              avcos(i,j,iband)=0.66564_r8
            END DO
          END IF
        END DO
      END DO
# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      DO iband=1,NBands
        CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,             &
     &                          LBi, UBi, LBj, UBj, SpecIr(:,:,iband))
        CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,             &
     &                          LBi, UBi, LBj, UBj, avcos(:,:,iband))
      END DO
# endif
      RETURN
      END SUBROUTINE ana_specir_tile
#endif

#if defined ANA_SSH && defined ZCLIMATOLOGY
      SUBROUTINE ana_ssh (ng, tile)
!
!=======================================================================
!                                                                      !
!  This routine sets analytical sea surface height climatology.        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_clima
!
      integer, intent(in) :: ng, tile

# include "tile.h"
!
      CALL ana_ssh_tile (ng, Istr, Iend, Jstr, Jend,                    &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   CLIMA(ng) % ssh)
      RETURN
      END SUBROUTINE ana_ssh
!
!***********************************************************************
      SUBROUTINE ana_ssh_tile (ng, Istr, Iend, Jstr, Jend,              &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         ssh)
!***********************************************************************
!
      USE mod_param

# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
!
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj

# ifdef ASSUMED_SHAPE
      real(r8), intent(out) :: ssh(LBi:,LBj:)
# else
      real(r8), intent(out) :: ssh(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set sea surface height (meters).
!-----------------------------------------------------------------------
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ssh(i,j)=???
        END DO
      END DO
# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, ssh)
# endif
      RETURN
      END SUBROUTINE ana_ssh_tile
#endif

#if defined SALINITY    && defined ANA_SSS && \
   (defined SCORRECTION || defined SRELAXATION)
      SUBROUTINE ana_sss (ng, tile)
!
!=======================================================================
!                                                                      !
!  This subroutine sets sea surface salinity SST (PSU) which is        !
!  used for surface water flux correction.                             !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_forces
!
      integer, intent(in) :: ng, tile

# include "tile.h"
!
      CALL ana_sss_tile (ng, Istr, Iend, Jstr, Jend,                    &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   FORCES(ng) % sss)
      RETURN
      END SUBROUTINE ana_sss
!
!***********************************************************************
      SUBROUTINE ana_sss_tile (ng, Istr, Iend, Jstr, Jend,              &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         sss)
!***********************************************************************
!
      USE mod_param

# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
!
      USE exchange_2d_mod
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj

# ifdef ASSUMED_SHAPE
      real(r8), intent(out) :: sss(LBi:,LBj:)
# else
      real(r8), intent(out) :: sss(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set sea surface salinity (PSU).
!-----------------------------------------------------------------------
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          sss(i,j)=???
        END DO
      END DO
      ANA_SST: no values provided for SST.
# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, sss)
# endif
      RETURN
      END SUBROUTINE ana_sss_tile
#endif

#if defined ANA_SST && defined QCORRECTION
      SUBROUTINE ana_sst (ng, tile)
!
!=======================================================================
!                                                                      !
!  This subroutine sets sea surface temperature SST  (Celsius)  and    !
!  surface net heat flux sensitivity dQdSTT to SST using analytical    !
!  expressions.  The forcing dQdSTT is usually computed in units of    !
!  (Watts/m2/degC).  It needs to be scaled to (m/s) by dividing by     !
!  rho0*Cp.  These forcing fields are used  when flux correction is    !
!  activated:                                                          !
!                                                                      !
!       Q_model ~ Q + dQdSST * (T_model - SST)                         !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_forces
!
      integer, intent(in) :: ng, tile

# include "tile.h"
!
      CALL ana_sst_tile (ng, Istr, Iend, Jstr, Jend,                    &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   FORCES(ng) % sst,                              &
     &                   FORCES(ng) % dqdt)
      RETURN
      END SUBROUTINE ana_sst
!
!***********************************************************************
      SUBROUTINE ana_sst_tile (ng, Istr, Iend, Jstr, Jend,              &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         sst, dqdt)
!***********************************************************************
!
      USE mod_param

# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
!
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj

# ifdef ASSUMED_SHAPE
      real(r8), intent(out) :: sst(LBi:,LBj:)
      real(r8), intent(out) :: dqdt(LBi:,LBj:)
# else
      real(r8), intent(out) :: sst(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: dqdt(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set sea surface temperature (Celsius) and heat flux sensitivity to
!  SST (Watts/m2).
!-----------------------------------------------------------------------
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          sst(i,j)=???
          dqdt(i,j)=???
        END DO
      END DO
      ANA_SST: no values provided for SST and DQDT.
# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, sst)
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, dqdt)
# endif
      RETURN
      END SUBROUTINE ana_sst_tile
#endif

#if defined ANA_STFLUX || defined ANA_SSFLUX || defined ANA_SPFLUX
      SUBROUTINE ana_stflux (ng, tile, itrc)
!
!=======================================================================
!                                                                      !
!  This routine sets kinematic surface flux of tracer type variables   !
!  "stflx" (tracer units m/s) using analytical expressions.            !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_forces
!
      integer, intent(in) :: ng, tile, itrc
!
# include "tile.h"
!
      CALL ana_stflux_tile (ng, Istr, Iend, Jstr, Jend, itrc,           &
     &                      LBi, UBi, LBj, UBj,                         &
# ifdef SHORTWAVE
     &                      FORCES(ng) % srflx,                         &
# endif
     &                      FORCES(ng) % stflx)
      RETURN
      END SUBROUTINE ana_stflux
!
!***********************************************************************
      SUBROUTINE ana_stflux_tile (ng, Istr, Iend, Jstr, Jend, itrc,     &
     &                            LBi, UBi, LBj, UBj,                   &
# ifdef SHORTWAVE
     &                            srflx,                                &
# endif
     &                            stflx)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars

# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
!
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj

# ifdef ASSUMED_SHAPE
#  ifdef SHORTWAVE
      real(r8), intent(in) :: srflx(LBi:,LBj:)
#  endif
      real(r8), intent(inout) :: stflx(LBi:,LBj:,:)
# else
#  ifdef SHORTWAVE
      real(r8), intent(in) :: srflx(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(inout) :: stflx(LBi:UBi,LBj:UBj,NT(ng))
# endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, itrc, j

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set kinematic surface heat flux (degC m/s) at horizontal
!  RHO-points.
!-----------------------------------------------------------------------
!
      IF (itrc.eq.itemp) THEN
        DO j=JstrR,JendR
          DO i=IstrR,IendR
# ifdef BL_TEST
            stflx(i,j,itrc)=srflx(i,j)
# else
            stflx(i,j,itrc)=0.0_r8
# endif
          END DO
        END DO
!
!-----------------------------------------------------------------------
!  Set kinematic surface freshwater flux (m/s) at horizontal
!  RHO-points, scaling by surface salinity is done in STEP3D.
!-----------------------------------------------------------------------
!
      ELSE IF (itrc.eq.isalt) THEN
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            stflx(i,j,itrc)=0.0_r8
          END DO
        END DO
!
!-----------------------------------------------------------------------
!  Set kinematic surface flux (T m/s) of passive tracers, if any.
!-----------------------------------------------------------------------
!
      ELSE
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            stflx(i,j,itrc)=0.0_r8
          END DO
        END DO
      END IF
# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, stflx(:,:,itrc))
# endif
      RETURN
      END SUBROUTINE ana_stflux_tile
#endif

#if defined ANA_TAIR && \
  ( defined BULK_FLUXES || defined ECOSIM || \
   (defined ANA_SRFLUX && defined ALBEDO) )
      SUBROUTINE ana_tair (ng, tile)
!
!=======================================================================
!                                                                      !
!  This routine sets surface air temperature (degC) using an           !
!  analytical expression.                                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_forces
!
      integer, intent(in) :: ng, tile

# include "tile.h"
!
      CALL ana_tair_tile (ng, Istr, Iend, Jstr, Jend,                   &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    FORCES(ng) % Tair)
      RETURN
      END SUBROUTINE ana_tair
!
!***********************************************************************
      SUBROUTINE ana_tair_tile (ng, Istr, Iend, Jstr, Jend,             &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Tair)
!***********************************************************************
!
      USE mod_param

# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
!
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj

# ifdef ASSUMED_SHAPE
      real(r8), intent(out) :: Tair(LBi:,LBj:)
# else
      real(r8), intent(out) :: Tair(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set analytical surface air temperature (degC).
!-----------------------------------------------------------------------
!
# if defined ADRIATIC1 || defined ADRIATIC2
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          Tair(i,j)=12.0_r8
        END DO
      END DO
# elif defined BENCHMARK1 || defined BENCHMARK2 || defined BENCHMARK3
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          Tair(i,j)=4.0_r8
        END DO
      END DO
# elif defined BL_TEST
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          Tair(i,j)=23.567_r8
        END DO
      END DO
# else
      ANA_TAIR: no values provided for TAIR.
# endif
# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, Tair)
# endif
      RETURN
      END SUBROUTINE ana_tair_tile
#endif

#if defined ANA_TCLIMA && defined TCLIMATOLOGY
      SUBROUTINE ana_tclima (ng, tile)
!
!=======================================================================
!                                                                      !
!  This routine sets analytical tracer climatology fields.             !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_clima
!
      integer, intent(in) :: ng, tile

# include  "tile.h"
!
      CALL ana_tclima_tile (ng, Istr, Iend, Jstr, Jend,                 &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      CLIMA(ng) % tclm)
      RETURN
      END SUBROUTINE ana_tclima
!
!***********************************************************************
      SUBROUTINE ana_tclima_tile (ng, Istr, Iend, Jstr, Jend,           &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            tclm)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_scalars

# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
!
      USE exchange_3d_mod, ONLY : exchange_r3d_tile
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj

# ifdef ASSUMED_SHAPE
      real(r8), intent(out) :: tclm(LBi:,LBj:,:,:)
# else
      real(r8), intent(out) :: tclm(LBi:UBi,LBj:UBj,N(ng),NT(ng))
# endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j, k
      real(r8) :: val1, val2, val3

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set tracer climatology.
!-----------------------------------------------------------------------
!
# if defined DOUBLE_GYRE
      val1=(44.69_r8/39.382_r8)**2
      val2=val1*(rho0*100.0_r8/g)*(5.0E-5_r8/((42.689_r8/44.69_r8)**2))
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            val3=T0(ng)+val2*EXP(GRID(ng)%z_r(i,j,k)/100.0_r8)*         &
     &           (10.0_r8-0.4_r8*TANH(GRID(ng)%z_r(i,j,k)/100.0_r8))
            tclm(i,j,k,itemp)=val3*(1.0_r8-3.0_r8*GRID(ng)%yr(i,j)/     &
     &                             (val3*el(ng)))
          END DO
        END DO
      END DO
# else
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            tclm(i,j,k,itemp)=???
            tclm(i,j,k,isalt)=???
          END DO
        END DO
      END DO
# endif
# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      DO itrc=1,NAT
        CALL exchange_r3d_tile (ng, Istr, Iend, Jstr, Jend,             &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          tclm(:,:,:,itrc))
      END DO
# endif
      RETURN
      END SUBROUTINE ana_tclima_tile
#endif

#ifdef ANA_TOBC
      SUBROUTINE ana_tobc (ng, tile)
!
!=======================================================================
!                                                                      !
!  This routine sets tracer-type variables open boundary conditions    !
!  using analytical expressions.                                       !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_boundary
      USE mod_grid
!
      integer, intent(in) :: ng, tile

# include "tile.h"
!
      CALL ana_tobc_tile (ng, Istr, Iend, Jstr, Jend,                   &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    GRID(ng) % z_r)
      RETURN
      END SUBROUTINE ana_tobc
!
!***********************************************************************
      SUBROUTINE ana_tobc_tile (ng, Istr, Iend, Jstr, Jend,             &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          z_r)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_boundary
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj

# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
# else
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
# endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, ised, itrc, j, k
      real(r8) :: cff

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Tracers open boundary conditions.
!-----------------------------------------------------------------------
!
# ifdef ESTUARY_TEST
      IF (EASTERN_EDGE) THEN
        DO k=1,N(ng)
          DO j=JstrR,JendR
            BOUNDARY(ng)%t_east(j,k,itemp)=T0(ng)
            BOUNDARY(ng)%t_east(j,k,isalt)=0.0_r8
#  if defined SEDIMENT
            DO ised=1,NST
              BOUNDARY(ng)%t_east(j,k,idsed(ised))=0.0_r8
            END DO
#  endif
          END DO
        END DO
      END IF
      IF (WESTERN_EDGE) THEN
        DO k=1,N(ng)
          DO j=JstrR,JendR
            BOUNDARY(ng)%t_west(j,k,itemp)=T0(ng)
            BOUNDARY(ng)%t_west(j,k,isalt)=30.0_r8
#  if defined SEDIMENT
            DO ised=1,NST
              BOUNDARY(ng)%t_west(j,k,idsed(ised))=0.0_r8
            END DO
#  endif
          END DO
        END DO
      END IF
# elif defined NJ_BIGHT
      IF (EASTERN_EDGE) THEN
        DO k=1,N(ng)
          DO j=JstrR,JendR
            IF (z_r(Iend+1,j,k).ge.-15.0_r8) THEN
              BOUNDARY(ng)%t_east(j,k,itemp)=2.04926425772840E+01_r8-   &
     &                                       z_r(Iend+1,j,k)*           &
     &                                       (2.64085084879392E-01_r8+  &
     &                                        z_r(Iend+1,j,k)*          &
     &                                        (2.75112532853521E-01_r8+ &
     &                                         z_r(Iend+1,j,k)*         &
     &                                        (9.20748976164887E-02_r8+ &
     &                                         z_r(Iend+1,j,k)*         &
     &                                        (1.44907572574284E-02_r8+ &
     &                                         z_r(Iend+1,j,k)*         &
     &                                        (1.07821568591208E-03_r8+ &
     &                                         z_r(Iend+1,j,k)*         &
     &                                        (3.24031805390397E-05_r8+ &
     &                                         1.26282685769027E-07_r8*
     &                                         z_r(Iend+1,j,k)))))))
              BOUNDARY(ng)%t_east(j,k,isalt)=3.06648914919313E+01_r8-   &
     &                                       z_r(Iend+1,j,k)*           &
     &                                       (1.47672526294673E-01_r8+  &
     &                                        z_r(Iend+1,j,k)*          &
     &                                        (1.12645576031340E-01_r8+ &
     &                                         z_r(Iend+1,j,k)*         &
     &                                        (3.90092328187102E-02_r8+ &
     &                                         z_r(Iend+1,j,k)*         &
     &                                        (6.93901493744710E-03_r8+ &
     &                                         z_r(Iend+1,j,k)*         &
     &                                        (6.60443669679294E-04_r8+ &
     &                                         z_r(Iend+1,j,k)*         &
     &                                        (3.19179236195422E-05_r8+ &
     &                                         6.17735263440932E-07_r8*
     &                                         z_r(Iend+1,j,k)))))))
            ELSE
              cff=TANH(1.1_r8*z_r(Iend+1,j,k)+15.9_r8)
              t_east(j,k,itemp)=14.6_r8+6.70_r8*cff
              t_east(j,k,isalt)=31.3_r8-0.55_r8*cff
            END IF
          END DO
        END DO
      END IF
      IF (SOUTHERN_EDGE) THEN
        DO k=1,N(ng)
          DO i=IstrR,IendR
            IF (z_r(i,Jstr-1,k).ge.-15.0_r8) THEN
              BOUNDARY(ng)%t_south(i,k,itemp)=2.04926425772840E+01_r8-  &
     &                                        z_r(i,Jstr-1,k)*          &
     &                                        (2.64085084879392E-01_r8+ &
     &                                         z_r(i,Jstr-1,k)*         &
     &                                         (2.75112532853521E-01_r8+&
     &                                          z_r(i,Jstr-1,k)*        &
     &                                         (9.20748976164887E-02_r8+&
     &                                          z_r(i,Jstr-1,k)*        &
     &                                         (1.44907572574284E-02_r8+&
     &                                          z_r(i,Jstr-1,k)*        &
     &                                         (1.07821568591208E-03_r8+&
     &                                          z_r(i,Jstr-1,k)*        &
     &                                         (3.24031805390397E-05_r8+&
     &                                          1.26282685769027E-07_r8*
     &                                          z_r(i,Jstr-1,k)))))))
              BOUNDARY(ng)%t_south(i,k,isalt)=3.06648914919313E+01_r8-  &
     &                                        z_r(i,Jstr-1,k)*          &
     &                                        (1.47672526294673E-01_r8+ &
     &                                         z_r(i,Jstr-1,k)*         &
     &                                         (1.12645576031340E-01_r8+&
     &                                          z_r(i,Jstr-1,k)*        &
     &                                         (3.90092328187102E-02_r8+&
     &                                          z_r(i,Jstr-1,k)*        &
     &                                         (6.93901493744710E-03_r8+&
     &                                          z_r(i,Jstr-1,k)*        &
     &                                         (6.60443669679294E-04_r8+&
     &                                          z_r(i,Jstr-1,k)*        &
     &                                         (3.19179236195422E-05_r8+&
     &                                          6.17735263440932E-07_r8*
     &                                          z_r(i,Jstr-1,k)))))))
            ELSE
              cff=TANH(1.1_r8*depth+15.9_r8)
              BOUNDARY(ng)%t_south(i,k,itemp)=14.6_r8+6.70_r8*cff
              BOUNDARY(ng)%t_south(i,k,isalt)=31.3_r8-0.55_r8*cff
            END IF
          END DO
        END DO
      END IF
# elif defined SED_TEST1
      IF (EASTERN_EDGE) THEN
        DO k=1,N(ng)
          DO j=JstrR,JendR
            BOUNDARY(ng)%t_east(j,k,itemp)=20.0_r8
            BOUNDARY(ng)%t_east(j,k,isalt)=0.0_r8
          END DO
        END DO
      END IF
# else
      IF (EASTERN_EDGE) THEN
        DO itrc=1,NT(ng)
          DO k=1,N(ng)
            DO j=JstrR,JendR
              BOUNDARY(ng)%t_east(j,k,itrc)=0.0_r8
            END DO
          END DO
        END DO
      END IF
      IF (WESTERN_EDGE) THEN
        DO itrc=1,NT(ng)
          DO k=1,N(ng)
            DO j=JstrR,JendR
              BOUNDARY(ng)%t_west(j,k,itrc)=0.0_r8
            END DO
          END DO
        END DO
      END IF
      IF (SOUTHERN_EDGE) THEN
        DO itrc=1,NT(ng)
          DO k=1,N(ng)
            DO i=IstrR,IendR
              BOUNDARY(ng)%t_south(i,k,itrc)=0.0_r8
            END DO
          END DO
        END DO
      END IF
      IF (NORTHERN_EDGE) THEN
        DO itrc=1,NT(ng)
          DO k=1,N(ng)
            DO i=IstrR,IendR
              BOUNDARY(ng)%t_north(i,k,itrc)=0.0_r8
            END DO
          END DO
        END DO
      END IF
# endif
      RETURN
      END SUBROUTINE ana_tobc_tile
#endif

#ifdef ANA_VMIX
      SUBROUTINE ana_vmix (ng, tile)
!
!=======================================================================
!                                                                      !
!  This routine sets vertical mixing coefficients for momentum "Akv"   !
!  and tracers "Akt" (m2/s) using analytical expressions.              !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_mixing
      USE mod_ocean
!
      integer, intent(in) :: ng, tile

# include "tile.h"
!
      CALL ana_vmix_tile (ng, Istr, Iend, Jstr, Jend,                   &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    GRID(ng) % h,                                 &
     &                    GRID(ng) % z_r,                               &
     &                    GRID(ng) % z_w,                               &
     &                    OCEAN(ng) % zeta,                             &
     &                    MIXING(ng) % Akv,                             &
     &                    MIXING(ng) % Akt)
      RETURN
      END SUBROUTINE ana_vmix
!
!***********************************************************************
      SUBROUTINE ana_vmix_tile (ng, Istr, Iend, Jstr, Jend,             &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          h, z_r, z_w, zeta, Akv, Akt)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars

# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
!
      USE exchange_3d_mod, ONLY : exchange_w3d_tile
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj

# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: zeta(LBi:,LBj:,:)
      real(r8), intent(out) :: Akv(LBi:,LBj:,0:)
      real(r8), intent(out) :: Akt(LBi:,LBj:,0:,:)
# else
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(in) :: zeta(LBi:UBi,LBj:UBj,3)
      real(r8), intent(out) :: Akv(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(out) :: Akt(LBi:UBi,LBj:UBj,0:N(ng),NAT)
# endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, itrc, j, k

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set vertical viscosity coefficient (m2/s).
!-----------------------------------------------------------------------
!
# if defined CANYON_B
      DO k=1,N(ng)-1
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            Akv(i,j,k)=1.0E-03_r8+95.0E-04_r8*EXP(z_w(i,j,k)/50.0_r8)+  &
     &                 95.0E-04_r8*EXP(-(z_w(i,j,k)+h(i,j))/50.0_r8)
          END DO
        END DO
      END DO
# elif defined COUPLING_TEST
      DO k=1,N(ng)-1
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            Akv(i,j,k)=2.0E-03_r8+8.0E-03_r8*EXP(z_w(i,j,k)/1500.0_r8)
          END DO
        END DO
      END DO
# elif defined ESTUARY_TEST
      DO k=1,N(ng)-1
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            Akv(i,j,k)=0.002_r8
          END DO
        END DO
      END DO
# elif defined LAKE_SIGNELL
      do k=1,N(ng)-1
        do j=JstrR,JendR
          do i=IstrR,IendR
            Akv(i,j,k)=0.0005_r8
          enddo
        enddo
      enddo
# elif defined NJ_BIGHT
      DO k=1,N(ng)-1
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            Akv(i,j,k)=1.0E-03_r8+2.0E-04_r8*EXP(z_r(i,j,k)/10.0_r8)
          END DO
        END DO
      END DO
# elif defined SED_TEST1
      DO k=1,N(ng)-1                         !  vonkar*ustar*z*(1-z/D)
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            Akv(i,j,k)=0.025_r8*(h(i,j)+z_w(i,j,k))*                    &
     &                 (1.0_r8-(h(i,j)+z_w(i,j,k))/                     &
     &                  (h(i,j)+zeta(i,j,knew)))
            Akt(i,j,k,itemp)=Akv(i,j,k)*0.49_r8/0.39_r8
            Akt(i,j,k,isalt)=Akt(i,j,k,itemp)
          END DO
        END DO
      END DO
# elif defined UPWELLING
      DO k=1,N(ng)-1
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            Akv(i,j,k)=2.0E-03_r8+8.0E-03_r8*EXP(z_w(i,j,k)/150.0_r8)
          END DO
        END DO
      END DO
# else
      ANA_VMIX: no values provided for AKV.
# endif
# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      CALL exchange_w3d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, 0, N(ng), Akv)
# endif
!
!-----------------------------------------------------------------------
!  Set vertical diffusion coefficient (m2/s).
!-----------------------------------------------------------------------
!
# if defined CANYON_B
      DO k=1,N(ng)-1
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            Akt(i,j,k,itemp)=Akt_bak(itemp,ng)
          END DO
        END DO
      END DO
# elif defined COUPLING_TEST
      DO k=1,N(ng)-1
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            Akt(i,j,k,itemp)=Akt_bak(itemp,ng)
            Akt(i,j,k,isalt)=Akt_bak(isalt,ng)
          END DO
        END DO
      END DO
# elif defined ESTUARY_TEST
      DO k=1,N(ng)-1
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            Akt(i,j,k,itemp)=Akv(i,j,k)
            Akt(i,j,k,isalt)=Akv(i,j,k)
          END DO
        END DO
      END DO
# elif defined LAKE_SIGNELL
      DO k=1,N(ng)-1
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            Akt(i,j,k,itemp)=Akt_bak(itemp,ng)
            Akt(i,j,k,isalt)=Akt_bak(isalt,ng)
          END DO
        END DO
      END DO
# elif defined NJ_BIGHT
      DO k=1,N(ng)-1
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            Akt(i,j,k,itemp)=1.0E-05_r8+                                &
     &                       2.0E-06_r8*EXP(z_r(i,j,k)/10.0_r8)
            Akt(i,j,k,isalt)=Akt(i,j,k,itemp)
          END DO
        END DO
      END DO
# elif defined UPWELLING
      DO k=1,N(ng)-1
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            Akt(i,j,k,itemp)=Akt_bak(itemp,ng)
            Akt(i,j,k,isalt)=Akt_bak(isalt,ng)
          END DO
        END DO
      END DO
# else
      ANA_VMIX: no values provided for AKT.
# endif
# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      DO itrc=1,NAT
        CALL exchange_w3d_tile (ng, Istr, Iend, Jstr, Jend,             &
     &                          LBi, UBi, LBj, UBj, 0, N(ng),           &
     &                          Akt(:,:,:,itrc))
      END DO
# endif
      RETURN
      END SUBROUTINE ana_vmix_tile
#endif

#if defined ANA_WINDS && (defined BULK_FLUXES || defined ECOSIM)
      SUBROUTINE ana_winds (ng, tile)
!
!=======================================================================
!                                                                      !
!  This routine sets surface wind components using an analytical       !
!  expression.                                                         !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_forces
      USE mod_grid
!
      integer, intent(in) :: ng, tile

# include "tile.h"
!
      CALL ana_winds_tile (ng, Istr, Iend, Jstr, Jend,                  &
     &                     LBi, UBi, LBj, UBj,                          &
# ifdef SPHERICAL
     &                     GRID(ng) % lonr,                             &
     &                     GRID(ng) % latr,                             &
# else
     &                     GRID(ng) % xr,                               &
     &                     GRID(ng) % yr,                               &
# endif
     &                     FORCES(ng) % Uwind,                          &
     &                     FORCES(ng) % Vwind)
      RETURN
      END SUBROUTINE ana_winds
!
!***********************************************************************
      SUBROUTINE ana_winds_tile (ng, Istr, Iend, Jstr, Jend,            &
     &                           LBi, UBi, LBj, UBj,                    &
# ifdef SPHERICAL
     &                           lonr, latr,                            &
# else
     &                           xr, yr,                                &
# endif
     &                           Uwind, Vwind)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars

# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
!
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj

# ifdef ASSUMED_SHAPE
#  ifdef SPHERICAL
      real(r8), intent(in) :: lonr(LBi:,LBj:)
      real(r8), intent(in) :: latr(LBi:,LBj:)
#  else
      real(r8), intent(in) :: xr(LBi:,LBj:)
      real(r8), intent(in) :: yr(LBi:,LBj:)
#  endif
      real(r8), intent(out) :: Uwind(LBi:,LBj:)
      real(r8), intent(out) :: Vwind(LBi:,LBj:)
# else
#  ifdef SPHERICAL
      real(r8), intent(in) :: lonr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: latr(LBi:UBi,LBj:UBj)
#  else
      real(r8), intent(in) :: xr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: yr(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(out) :: Uwind(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: Vwind(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j
      real(r8) :: Wdir, Wmag, cff, u_wind, v_wind

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set surface wind components (m/s) at RHO-points.
!-----------------------------------------------------------------------
!
# if defined BENCHMARK1 || defined BENCHMARK2 || defined BENCHMARK3
      Wmag=15.0_r8
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          cff=0.2_r8*(60.0_r8+latr(i,j))
          Uwind(i,j)=Wmag*EXP(-cff*cff)
          Vwind(i,j)=0.0_r8
        END DO
      END DO
# elif defined BL_TEST
      IF ((tdays(ng)-dstart).le.6.0_r8) THEN
        u_wind=0.0_r8
!!      v_wind=4.7936_r8
        v_wind=10.0_r8
      END IF
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          Uwind(i,j)=u_wind
          Vwind(i,j)=v_wind
        END DO
      END DO
# else
      ANA_WINDS: no values provided for UWIND and VWIND.
# endif
# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, Uwind)
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, Vwind)
# endif
      RETURN
      END SUBROUTINE ana_winds_tile
#endif

#if defined ANA_WWAVE && defined BBL_MODEL
      SUBROUTINE ana_wwave (ng, tile)
!
!=======================================================================
!                                                                      !
!  This subroutine sets wind induced wave amplitude, direction and     !
!  period to be used in the bottom boundary layer formulation.         !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_forces
      USE mod_grid
!
      integer, intent(in) :: ng, tile

# include "tile.h"
!
      CALL ana_wwave_tile (ng, Istr, Iend, Jstr, Jend,                  &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     GRID(ng) % angler,                           &
     &                     FORCES(ng) % Awave,                          &
     &                     FORCES(ng) % Dwave,                          &
     &                     FORCES(ng) % Pwave)
      RETURN
      END SUBROUTINE ana_wwave
!
!***********************************************************************
      SUBROUTINE ana_wwave_tile (ng, Istr, Iend, Jstr, Jend,            &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           angler, Awave, Dwave, Pwave)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars

# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
!
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj

# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: angler(LBi:,LBj:)
      real(r8), intent(out) :: Awave(LBi:,LBj:)
      real(r8), intent(out) :: Dwave(LBi:,LBj:)
      real(r8), intent(out) :: Pwave(LBi:,LBj:)
# else
      real(r8), intent(in) :: angler(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: Awave(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: Dwave(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: Pwave(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j
      real(r8) :: cff, wdir
#if defined LAKE_SIGNELL
      real(r8) :: mxtm, rmptm, mxst
#endif

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set wind induced wave amplitude (m), direction (radians) and
!  period (s) at RHO-points.
!-----------------------------------------------------------------------
!
# if defined BL_TEST
      wdir=210.0_r8*deg2rad
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          Awave(i,j)=0.5_r8
          Dwave(i,j)=wdir
          Pwave(i,j)=8.0_r8
        END DO
      END DO
# elif defined LAKE_SIGNELL
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          Dwave(i,j)=0.0_r8*deg2rad
!
          mxtm=48.0_r8                  ! hours
          rmptm=3.0_r8                  ! hours
          mxst=0.5_r8                   ! m
          cff=MIN(mxst*(TANH((time(ng)/60.0_r8-rmptm*150.0_r8)/         &
     &                       150.0_r8)+1.0_r8)*0.5_r8,                  &
     &            mxst*(TANH(((mxtm*3600.0_r8-time(ng))/                &
     &                        60.0_r8-rmptm*150.0_r8)/150.0_r8)+        &
     &                  1.0_r8)*0.5_r8)
          Awave(i,j)=MAX(cff,0.01_r8)
!
          mxst=5.0_r8                   ! sec
          cff=MIN(mxst*(TANH((time(ng)/60.0_r8-rmptm*150.0_r8)/         &
     &                       150.0_r8)+1.0_r8)*0.5_r8,                  &
     &            mxst*(TANH(((mxtm*3600.0_r8-time(ng))/                &
     &                        60.0_r8-rmptm*150.0_r8)/150.0_r8)+        &
     &                  1.0_r8)/2.0_r8)
          Pwave(i,j)=MAX(cff,2.0_r8)
        END DO
      END DO
# elif defined NJ_BIGHT
!!    wdir=210.0_r8*deg2rad
      wdir=150.0_r8*deg2rad
      IF ((tdays(ng)-dstart).lt.1.5_r8) THEN
        cff=TANH(0.5_r8*(tdays(ng)-dstart))
        cff=1.0_r8
      ELSE
        cff=1.0_r8
      END IF
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          Awave(i,j)=0.5_r8
          Dwave(i,j)=wdir-angler(i,j)
          Pwave(i,j)=10.0_r8
        END DO
      END DO
# else
      ANA_WWAVE: no values provided for AWAVE, DWAVE, and PWAVE.
# endif
# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, Awave)
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, Dwave)
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, Pwave)
# endif
      RETURN
      END SUBROUTINE ana_wwave_tile
#endif
      END MODULE analytical_mod
