      SUBROUTINE ad_biology (ng,tile)
!
!git $Id$
!************************************************** Hernan G. Arango ***
!  Copyright (c) 2002-2024 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!***********************************************************************
!                                                                      !
!  IOP-based, CDOM (Colored Dissolved Organic Matter) Model            !
!                                                                      !
!  This routine computes the biological sources and sinks and adds     !
!  then the global biological fields.                                  !
!                                                                      !
!  It computes the absorption, scattering, and backscattering from     !
!  the ecosystem variables.  This  facilitates  direct comparisons     !
!  between measuared and models IOPs (Inherent Optical Properties).    !
!  It also possible to calculate remote sensing  reflectances from     !
!  model output  for direct comparisons with remotely sensed ocean     !
!  color (SeaWIFS, MODIS, etc) and data assimilation.                  !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    Fennel, K., ...                                                   !
!                                                                      !
!***********************************************************************
!
      USE mod_param
      USE mod_forces
      USE mod_grid
      USE mod_ncparam
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
!  Set header file name.
!
#ifdef DISTRIBUTE
      IF (Lbiofile(iADM)) THEN
#else
      IF (Lbiofile(iADM).and.(tile.eq.0)) THEN
#endif
        Lbiofile(iADM)=.FALSE.
        BIONAME(iADM)=__FILE__
      END IF
!
#ifdef PROFILE
      CALL wclock_on (ng, iADM, 15)
#endif
      CALL ad_biology_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj, N(ng), NT(ng),          &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      nstp(ng), nnew(ng),                         &
#ifdef MASKING
     &                      GRID(ng) % rmask,                           &
#endif
     &                      GRID(ng) % Hz,                              &
     &                      GRID(ng) % ad_Hz,                           &
     &                      GRID(ng) % z_w,                             &
     &                      GRID(ng) % ad_z_w,                          &
     &                      FORCES(ng) % srflx,                         &
     &                      FORCES(ng) % ad_srflx,                      &
     &                      OCEAN(ng) % t,                              &
     &                      OCEAN(ng) % ad_t)

#ifdef PROFILE
      CALL wclock_off (ng, iADM, 15)
#endif

      RETURN
      END SUBROUTINE ad_biology
!
!-----------------------------------------------------------------------
      SUBROUTINE ad_biology_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, UBk, UBt,         &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            nstp, nnew,                           &
#ifdef MASKING
     &                            rmask,                                &
#endif
     &                            Hz, ad_Hz,                            &
     &                            z_w, ad_z_w,                          &
     &                            srflx, ad_srflx,                      &
     &                            t, ad_t)
!-----------------------------------------------------------------------
!
      USE mod_param
      USE mod_biology
      USE mod_ncparam
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, UBk, UBt
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp, nnew

#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
# endif
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: srflx(LBi:,LBj:)
      real(r8), intent(in) :: t(LBi:,LBj:,:,:,:)

      real(r8), intent(inout) :: ad_Hz(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_z_w(LBi:,LBj:,0:)
      real(r8), intent(inout) :: ad_srflx(LBi:,LBj:)
      real(r8), intent(inout) :: ad_t(LBi:,LBj:,:,:,:)
#else
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:UBk)
      real(r8), intent(in) :: srflx(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: t(LBi:UBi,LBj:UBj,UBk,3,UBt)

      real(r8), intent(inout) :: ad_Hz(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(inout) :: ad_z_w(LBi:UBi,LBj:UBj,0:UBk)
      real(r8), intent(inout) :: ad_srflx(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_t(LBi:UBi,LBj:UBj,UBk,3,UBt)
#endif
!
!  Local variable declarations.
!
      integer :: Iter, i, ibio, itrc, j, k, kk, nb
      integer :: Iteradj

      real(r8) :: Att, ExpAtt, Itop, PAR, PAR1
      real(r8) :: ad_Att, ad_ExpAtt, ad_Itop, ad_PAR
      real(r8) :: cff, cff1, dtdays
      real(r8) :: ad_cff
      real(r8) :: adfac, adfac1

      real(r8), dimension(IminS:ImaxS) :: PARsur
      real(r8), dimension(IminS:ImaxS) :: ad_PARsur

      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: Bio
      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: Bio1
      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: Bio_old

      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: ad_Bio
      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: ad_Bio_old

      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Light

      real(r8), dimension(IminS:ImaxS,N(ng)) :: ad_Hz_inv
      real(r8), dimension(IminS:ImaxS,N(ng)) :: ad_Light

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Add biological Source/Sink terms.
!-----------------------------------------------------------------------
!
!  Set time-stepping size (days) according to the number of iterations.
!
      dtdays=dt(ng)*sec2day/REAL(BioIter(ng),r8)
!
!  Initialize private scratch variables and array.
!
      ad_PAR=0.0_r8
      ad_Att=0.0_r8
      ad_ExpAtt=0.0_r8
      ad_Itop=0.0_r8
      ad_cff=0.0_r8
!
      DO k=1,N(ng)
        DO i=IminS,ImaxS
          ad_Hz_inv(i,k)=0.0_r8
          ad_Light(i,k)=0.0_r8
        END DO
      END DO
      DO i=IminS,ImaxS
        ad_PARsur(i)=0.0_r8
      END DO
      DO itrc=1,NBT
        ibio=idbio(itrc)
        DO k=1,N(ng)
          DO i=Istr,Iend
            Bio(i,k,ibio)=0.0_r8
            Bio1(i,k,ibio)=0.0_r8
            ad_Bio(i,k,ibio)=0.0_r8
            ad_Bio_old(i,k,ibio)=0.0_r8
          END DO
        END DO
      END DO
!
!  Start pipelined J-loop.
!
      J_LOOP : DO j=Jstr,Jend
!
!  Compute inverse thickness to avoid repeated divisions.
!
        DO k=1,N(ng)
          DO i=Istr,Iend
            Hz_inv(i,k)=1.0_r8/Hz(i,j,k)
          END DO
        END DO
!
!  Compute the required basic state arrays.
!
!  Extract biological variables from tracer arrays, place them into
!  scratch arrays, and restrict their values to be positive definite.
!  At input, all tracers (index nnew) from predictor step have
!  transport units (m Tunits) since we do not have yet the new
!  values for zeta and Hz. These are known after the 2D barotropic
!  time-stepping.
!
        DO itrc=1,NBT
          ibio=idbio(itrc)
          DO k=1,N(ng)
            DO i=Istr,Iend
              Bio(i,k,ibio)=MAX(0.0_r8,t(i,j,k,nstp,ibio))
              Bio_old(i,k,ibio)=Bio(i,k,ibio)
            END DO
          END DO
        END DO
!
!  Calculate surface Photosynthetically Available Radiation (PAR).  The
!  net shortwave radiation is scaled back to Watts/m2 and multiplied by
!  the fraction that is photosynthetically available, PARfrac.
!
        DO i=Istr,Iend
          PARsur(i)=PARfrac(ng)*srflx(i,j)*rho0*Cp
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
!  in a fully implicit manner, so the algorithm guarantees non-negative
!  values, no matter how strong the concentration of active consuming
!  component (Phytoplankton or Zooplankton).  The overall algorithm,
!  as well as any stage of it, is formulated in conservative form
!  (except explicit sinking) in sense that the sum of concentration of
!  all components is conserved.
!
!  In the implicit algorithm, we have for example (N: nutrient,
!                                                  P: phytoplankton),
!
!     N(new) = N(old) - uptake * P(old)     uptake = mu * N / (Kn + N)
!                                                    {Michaelis-Menten}
!  below, we set
!                                           The N in the numerator of
!     cff = mu * P(old) / (Kn + N(old))     uptake is treated implicitly
!                                           as N(new)
!
!  so the time-stepping of the equations becomes:
!
!     N(new) = N(old) / (1 + cff)     (1) when substracting a sink term,
!                                         consuming, divide by (1 + cff)
!  and
!
!     P(new) = P(old) + cff * N(new)  (2) when adding a source term,
!                                         growing, add (cff * source)
!
!  Notice that if you substitute (1) in (2), you will get:
!
!     P(new) = P(old) + cff * N(old) / (1 + cff)    (3)
!
!  If you add (1) and (3), you get
!
!     N(new) + P(new) = N(old) + P(old)
!
!  implying conservation regardless how "cff" is computed. Therefore,
!  this scheme is unconditionally stable regardless of the conversion
!  rate. It does not generate negative values since the constituent
!  to be consumed is always treated implicitly. It is also biased
!  toward damping oscillations.
!
!  The iterative loop below is to iterate toward an universal Backward-
!  Euler treatment of all terms. So if there are oscillations in the
!  system, they are only physical oscillations. These iterations,
!  however, do not improve the accuaracy of the solution.
!
        ITER_LOOP: DO Iter=1,BioIter(ng)
!
!  Compute light attenuation as function of depth.
!
          DO i=Istr,Iend
            PAR=PARsur(i)
            IF (PARsur(i).gt.0.0_r8) THEN              ! day time
              DO k=N(ng),1,-1
!
!  Attenuate the light to the center of the grid cell. Here, AttSW is
!  the light attenuation due to seawater.
!
                Att=(AttSW(ng)+                                         &
     &               CDOM_LightAtt(ng)*Bio(i,k,aCDOM(i440n)))*          &
     &              (z_w(i,j,k)-z_w(i,j,k-1))
                ExpAtt=EXP(-Att)
                Itop=PAR
                PAR=Itop*(1.0_r8-ExpAtt)/Att    ! average at cell center
                Light(i,k)=PAR
!
!  Light attenuation at the bottom of the grid cell. It is the starting
!  PAR value for the next (deeper) vertical grid cell.
!
                PAR=Itop*ExpAtt
              END DO
            ELSE                                       ! night time
              DO k=1,N(ng)
                Light(i,k)=0.0_r8
              END DO
            END IF
          END DO
!
!  Degradation of CDOM absorption, aCDOM.
!
          cff1=dtdays*CDOM_sigma(ng)
          DO k=1,N(ng)
            DO i=Istr,Iend
              cff=cff1*Light(i,k)
              DO nb=1,NBands
                Bio(i,k,aCDOM(nb))=Bio(i,k,aCDOM(nb))/(1.0_r8+cff)
              END DO
            END DO
          END DO

        END DO ITER_LOOP
!
!-----------------------------------------------------------------------
!  Adjoint of update global tracer variables: Add increment due to BGC
!  processes to tracer array in time index "nnew". Index "nnew" is
!  solution after advection and mixing and has transport units
!  (m Tunits) hence the increment is multiplied by Hz.  Notice that
!  we need to subtract original values "Bio_old" at the top of the
!  routine to just account for the concentractions affected by BGC
!  processes. This also takes into account any constraints
!  (non-negative concentrations, carbon concentration range) specified
!  before entering BGC kernel. If "Bio" were unchanged by BGC
!  processes, the increment would be exactly zero. Notice that final
!  tracer values, t(:,:,:,nnew,:) are not bounded >=0 so that we can
!  preserve total inventory of nutrients when advection causes tracer
!  concentration to go negative.
!-----------------------------------------------------------------------
!
        DO itrc=1,NBT
          ibio=idbio(itrc)
          DO k=1,N(ng)
            DO i=Istr,Iend
              cff=Bio(i,k,ibio)-Bio_old(i,k,ibio)
!^            tl_t(i,j,k,nnew,ibio)=tl_t(i,j,k,nnew,ibio)+              &
!^   &                              tl_cff*Hz(i,j,k)+cff*tl_Hz(i,j,k)
!^
              ad_Hz(i,j,k)=ad_Hz(i,j,k)+cff*ad_t(i,j,k,nnew,ibio)
              ad_cff=add_cff+Hz(i,j,k)*ad_t(i,j,k,nnew,ibio)
!^            tl_cff=tl_Bio(i,k,ibio)-tl_Bio_old(i,k,ibio)
!^
              ad_Bio_old(i,k,ibio)=ad_Bio_old(i,k,ibio)-ad_cff
              ad_Bio(i,k,ibio)=ad_Bio(i,k,ibio)+ad_cff
              ad_cff=0.0_r8
            END DO
          END DO
        END DO
!
!=======================================================================
!  Start internal iterations to achieve convergence of the nonlinear
!  backward-implicit solution.
!=======================================================================
!
        ITER_LOOP1: DO Iter=BioIter(ng),1,-1
!
!  Compute appropriate basic state arrays.
!
!
!  Extract biological variables from tracer arrays, place them into
!  scratch arrays, and restrict their values to be positive definite.
!  At input, all tracers (index nnew) from predictor step have
!  transport units (m Tunits) since we do not have yet the new
!  values for zeta and Hz. These are known after the 2D barotropic
!  time-stepping.
!
          DO itrc=1,NBT
            ibio=idbio(itrc)
            DO k=1,N(ng)
              DO i=Istr,Iend
                Bio(i,k,ibio)=MAX(0.0_r8,t(i,j,k,nstp,ibio))
                Bio_old(i,k,ibio)=Bio(i,k,ibio)
              END DO
            END DO
          END DO
!
!  Calculate surface Photosynthetically Available Radiation (PAR).  The
!  net shortwave radiation is scaled back to Watts/m2 and multiplied by
!  the fraction that is photosynthetically available, PARfrac.
!
          DO i=Istr,Iend
            PARsur(i)=PARfrac(ng)*srflx(i,j)*rho0*Cp
          END DO
!
!=======================================================================
!  Start internal iterations to achieve convergence of the nonlinear
!  backward-implicit solution.
!=======================================================================
!
          DO Iteradj=1,Iter
!
!  Compute light attenuation as function of depth.
!
            DO i=Istr,Iend
              PAR=PARsur(i)
              IF (PARsur(i).gt.0.0_r8) THEN            ! day time
                DO k=N(ng),1,-1
!
!  Attenuate the light to the center of the grid cell. Here, AttSW is
!  the light attenuation due to seawater.
!
                  Att=(AttSW(ng)+                                       &
     &                 CDOM_LightAtt(ng)*Bio(i,k,aCDOM(i440n)))*        &
     &                (z_w(i,j,k)-z_w(i,j,k-1))
                  ExpAtt=EXP(-Att)
                  Itop=PAR
                  PAR=Itop*(1.0_r8-ExpAtt)/Att  ! average at cell center
                  Light(i,k)=PAR
!
!  Light attenuation at the bottom of the grid cell. It is the starting
!  PAR value for the next (deeper) vertical grid cell.
!
                  PAR=Itop*ExpAtt
                END DO
              ELSE                                     ! night time
                DO k=1,N(ng)
                  Light(i,k)=0.0_r8
                END DO
              END IF
            END DO
!
!  Degradation of CDOM absorption, aCDOM.
!
            cff1=dtdays*CDOM_sigma(ng)
            DO k=1,N(ng)
              DO i=Istr,Iend
                cff=cff1*Light(i,k)
                DO nb=1,NBands
                  Bio1(i,k,aCDOM(nb))=Bio(i,k,aCDOM(nb))
                  Bio(i,k,aCDOM(nb))=Bio(i,k,aCDOM(nb))/(1.0_r8+cff)
                END DO
              END DO
            END DO
          END DO
!
!  End of compute basic state arrays.
!
!  Adjoint of degradation of CDOM absorption, aCDOM.
!
          cff1=dtdays*CDOM_sigma(ng)
          DO k=1,N(ng)
            DO i=Istr,Iend
              cff=cff1*Light(i,k)
              DO nb=1,NBands
!^              tl_Bio(i,k,aCDOM(nb))=(tl_Bio(i,k,aCDOM(nb))-           &
!^   &                                 tl_cff*Bio(i,k,aCDOM(nb)))/      &
!^   &                                (1.0_r8+cff)
!^
                adfac=ad_Bio(i,k,aCDOM(nb))/(1.0_r8+cff)
                ad_cff=ad_cff-Bio(i,k,aCDOM(nb))*adfac
                ad_Bio(i,k,aCDOM(nb))=adfac
              END DO
!^            tl_cff=cff1*tl_Light(i,k)
!^
               ad_Light(i,k)=ad_Light(i,k)+cff1*ad_cff
               ad_cff=0.0_r8
            END DO
          END DO
!
!  Compute adjoint light attenuation as function of depth.
!
          DO i=Istr,Iend
            PAR=PARsur(i)
            IF (PARsur(i).gt.0.0_r8) THEN              ! day time
              DO k=1,N(ng)
!
! Compute the basic state PAR appropriate for each level.
!
                PAR=PARsur(i)
                DO kk=N(ng),k,-1
!
!  Compute average light attenuation for each grid cell. Here, AttSW is
!  the light attenuation due to seawater and AttPhy is the attenuation
!  due to phytoplankton (self-shading coefficient).
!
                  Att=(AttSW(ng)+                                       &
     &                 CDOM_LightAtt(ng)*Bio1(i,kk,aCDOM(i440n)))*      &
     &                (z_w(i,j,kk)-z_w(i,j,kk-1))
                  ExpAtt=EXP(-Att)
                  Itop=PAR
                  PAR=Itop*(1.0_r8-ExpAtt)/Att      ! average at center
                  PAR1=PAR
!
!  Light attenuation at the bottom of the grid cell. It is the starting
!  PAR value for the next (deeper) vertical grid cell.
!
                  PAR=Itop*ExpAtt
                END DO
!
!  Adjoint of light attenuation at the bottom of the grid cell. It is
!  the starting PAR value for the next (deeper) vertical grid cell.
!
!^              tl_PAR=tl_Itop*ExpAtt+Itop*tl_ExpAtt
!^
                ad_ExpAtt=ad_ExpAtt+Itop*ad_PAR
                ad_Itop=ad_Itop+ExpAtt*ad_PAR
                ad_PAR=0.0_r8
!
!  Adjoint of compute average light attenuation for each grid cell.
!  Here, AttSW is the light attenuation due to seawater and AttPhy is
!  the attenuation due to phytoplankton (self-shading coefficient).
!
!^              tl_Light(i,k)=tl_PAR
!^
                ad_PAR=ad_PAR+ad_Light(i,k)
                ad_Light(i,k)=0.0_r8
!^              tl_PAR=(-tl_Att*PAR1+tl_Itop*(1.0_r8-ExpAtt)-           &
!^   &                  Itop*tl_ExpAtt)/Att
!^
                adfac=ad_PAR/Att
                ad_Att=ad_Att-PAR1*adfac
                ad_ExpAtt=ad_ExpAtt-Itop*adfac
                ad_Itop=ad_Itop+(1.0_r8-ExpAtt)*adfac
                ad_PAR=0.0_r8
!^              tl_Itop=tl_PAR
!^
                ad_PAR=ad_PAR+ad_Itop
                ad_Itop=0.0_r8
!^              tl_ExpAtt=-ExpAtt*tl_Att
!^
                ad_Att=ad_Att-ExpAtt*ad_ExpAtt
                ad_ExpAtt=0.0_r8
!^              tl_Att=CDOM_LightAtt(ng)*tl_Bio(i,k,aCDOM(i440n))*      &
!^   &                 (z_w(i,j,k)-z_w(i,j,k-1))+                       &
!^   &                 (AttSW(ng)+                                      &
!^   &                  CDOM_LightAtt(ng)*Bio1(i,k,aCDOM(i440n)))*      &
!^   &                 (tl_z_w(i,j,k)-tl_z_w(i,j,k-1))
!^
                adfac=(AttSW(ng)+                                       &
     &                 CDOM_LightAtt(ng)*Bio1(i,k,aCDOM(i440n)))*       &
     &                ad_Att
                ad_Bio(i,k,aCDOM(i440n))=ad_Bio(i,k,aCDOM(i440n))+      &
     &                                   CDOM_LightAtt(ng)*             &
     &                                   (z_w(i,j,k)-z_w(i,j,k-1))*     &
     &                                   ad_Att
                ad_z_w(i,j,k-1)=ad_z_w(i,j,k-1)-adfac
                ad_z_w(i,j,k  )=ad_z_w(i,j,k  )+adfac
                ad_Att=0.0_r8
              END DO
            ELSE                                       ! night time
              DO k=1,N(ng)
!^              tl_Light(i,k)=0.0_r8
!^
                ad_Light(i,k)=0.0_r8
              END DO
            END IF
!^          tl_PAR=tl_PARsur(i)
!^
            ad_PARsur(i)=ad_PARsur(i)+ad_PAR
            ad_PAR=0.0_r8
          END DO

        END DO ITER_LOOP1
!
!  Calculate adjoint surface Photosynthetically Available Radiation
!  (PAR).  The net shortwave radiation is scaled back to Watts/m2
!  and multiplied by the fraction that is photosynthetically
!  available, PARfrac.
!
        DO i=Istr,Iend
!^        tl_PARsur(i)=(tl_PARfrac(ng)*srflx(i,j)+                      &
!^   &                  PARfrac(ng)*tl_srflx(i,j))*rho0*Cp
!^
          adfac=rho0*Cp*ad_PARsur(i)
          ad_srflx(i,j)=ad_srflx(i,j)+PARfrac(ng)*adfac
          ad_PARfrac(ng)=ad_PARfrac(ng)+srflx(i,j)*adfac
          ad_PARsur(i)=0.0_r8
        END DO
!
!  Extract biological variables from tracer arrays, place them into
!  scratch arrays, and restrict their values to be positive definite.
!
        DO itrc=1,NBT
          ibio=idbio(itrc)
          DO k=1,N(ng)
            DO i=Istr,Iend
!^            tl_Bio_old(i,k,ibio)=tl_Bio(i,k,ibio)
!^
              ad_Bio(i,k,ibio)=ad_Bio(i,k,ibio)+                        &
     &                         ad_Bio_old(i,k,ibio)
              ad_Bio_old(i,k,ibio)=0.0_r8
!^            tl_Bio(i,k,ibio)=(0.5_r8-                                 &
!^   &                          SIGN(0.5_r8,-t(i,j,k,nstp,ibio)))*      &
!^   &                         tl_t(i,j,k,nstp,ibio)
!^
              ad_t(i,j,k,nstp,ibio)=ad_t(i,j,k,nstp,ibio)+              &
     &                              (0.5_r8-                            &
     &                               SIGN(0.5_r8,-t(i,j,k,nstp,ibio)))* &
     &                              ad_Bio(i,k,ibio)
              ad_Bio(i,k,ibio)=0.0_r8
            END DO
          END DO
        END DO
!
!  Adjoint inverse thickness to avoid repeated divisions.
!
        DO k=1,N(ng)
          DO i=Istr,Iend
!^          tl_Hz_inv(i,k)=-Hz_inv(i,k)*Hz_inv(i,k)*tl_Hz(i,j,k)
!^
            ad_Hz(i,j,k)=ad_Hz(i,j,k)-                                  &
     &                   Hz_inv(i,k)*Hz_inv(i,k)*ad_Hz_inv(i,k)
            ad_Hz_inv(i,k)=0.0_r8
          END DO
        END DO

      END DO J_LOOP

      RETURN
      END SUBROUTINE ad_biology_tile
