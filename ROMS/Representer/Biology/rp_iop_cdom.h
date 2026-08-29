      SUBROUTINE rp_biology (ng,tile)
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
      IF (Lbiofile(iRPM)) THEN
#else
      IF (Lbiofile(iRPM).and.(tile.eq.0)) THEN
#endif
        Lbiofile(iRPM)=.FALSE.
        BIONAME(iRPM)=__FILE__
      END IF
!
#ifdef PROFILE
      CALL wclock_on (ng, iRPM, 15)
#endif
      CALL rp_biology_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj, N(ng), NT(ng),          &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      nstp(ng), nnew(ng),                         &
#ifdef MASKING
     &                      GRID(ng) % rmask,                           &
#endif
     &                      GRID(ng) % Hz,                              &
     &                      GRID(ng) % tl_Hz,                           &
     &                      GRID(ng) % z_w,                             &
     &                      GRID(ng) % tl_z_w,                          &
     &                      FORCES(ng) % srflx,                         &
     &                      FORCES(ng) % tl_srflx,                      &
     &                      OCEAN(ng) % t,                              &
     &                      OCEAN(ng) % tl_t)

#ifdef PROFILE
      CALL wclock_off (ng, iRPM, 15)
#endif

      RETURN
      END SUBROUTINE rp_biology
!
!-----------------------------------------------------------------------
      SUBROUTINE rp_biology_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, UBk, UBt,         &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            nstp, nnew,                           &
#ifdef MASKING
     &                            rmask,                                &
#endif
     &                            Hz, tl_Hz,                            &
     &                            z_w, tl_z_w,                          &
     &                            srflx, tl_srflx,                      &
     &                            t, tl_t)
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

      real(r8), intent(in) :: tl_Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: tl_z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: tl_srflx(LBi:,LBj:)
      real(r8), intent(inout) :: tl_t(LBi:,LBj:,:,:,:)
#else
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:UBk)
      real(r8), intent(in) :: srflx(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: t(LBi:UBi,LBj:UBj,UBk,3,UBt)

      real(r8), intent(in) :: tl_Hz(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(in) :: tl_z_w(LBi:UBi,LBj:UBj,0:UBk)
      real(r8), intent(in) :: tl_srflx(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: tl_t(LBi:UBi,LBj:UBj,UBk,3,UBt)
#endif
!
!  Local variable declarations.
!
      integer :: Iter, i, ibio, itrc, j, k, nb
      integer :: Iteradj

      real(r8) :: Att, ExpAtt, Itop, PAR
      real(r8) :: tl_Att, tl_ExpAtt, tl_Itop, tl_PAR
      real(r8) :: cff, cff1, dtdays
      real(r8) :: tl_cff

      real(r8), dimension(IminS:ImaxS) :: PARsur
      real(r8), dimension(IminS:ImaxS) :: tl_PARsur

      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: Bio
      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: Bio_old

      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: tl_Bio
      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: tl_Bio_old

      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Light

      real(r8), dimension(IminS:ImaxS,N(ng)) :: tl_Hz_inv
      real(r8), dimension(IminS:ImaxS,N(ng)) :: tl_Light

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
!  Initialize local Bio scratch arrays.
!
      DO itrc=1,NBT
        ibio=idbio(itrc)
        DO k=1,N(ng)
          DO i=Istr,Iend
            Bio(i,k,ibio)=0.0_r8
            tl_Bio(i,k,ibio)=0.0_r8
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
            tl_Hz_inv(i,k)=-Hz_inv(i,k)*Hz_inv(i,k)*tl_Hz(i,j,k)+       &
#ifdef TL_IOMS
     &                     2.0_r8*Hz_inv(i,k)
#endif
          END DO
        END DO
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
              tl_Bio(i,k,ibio)=(0.5_r8-                                 &
     &                          SIGN(0.5_r8,-t(i,j,k,nstp,ibio)))*      &
     &                         tl_t(i,j,k,nstp,ibio)
              Bio_old(i,k,ibio)=Bio(i,k,ibio)
              tl_Bio_old(i,k,ibio)=tl_Bio(i,k,ibio)
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
          tl_PARsur(i)=(tl_PARfrac(ng)*srflx(i,j)+                      &
     &                  PARfrac(ng)*tl_srflx(i,j))*rho0*Cp-             &
#ifdef TL_IOMS
     &                 PARsur(i)
#endif
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
#ifdef TL_IOMS
            tl_PAR=PARsur(i)
#else
            tl_PAR=tl_PARsur(i)
#endif
            IF (PARsur(i).gt.0.0_r8) THEN              ! day time
              DO k=N(ng),1,-1
!
!  Attenuate the light to the center of the grid cell. Here, AttSW is
!  the light attenuation due to seawater.
!
                Att=(AttSW(ng)+                                         &
     &               CDOM_LightAtt(ng)*Bio(i,k,aCDOM(i440n)))*          &
     &              (z_w(i,j,k)-z_w(i,j,k-1))
                tl_Att=CDOM_LightAtt(ng)*tl_Bio(i,k,aCDOM(i440n))*      &
     &                 (z_w(i,j,k)-z_w(i,j,k-1))+                       &
     &                 (AttSW(ng)+                                      &
     &                  CDOM_LightAtt(ng)*Bio(i,k,aCDOM(i440n)))*       &
     &                 (tl_z_w(i,j,k)-tl_z_w(i,j,k-1))-                 &
#ifdef TL_IOMS
     &                 CDOM_LightAtt(ng)*Bio1(i,k,aCDOM(i440n))*        &
     &                 (z_w(i,j,k)-z_w(i,j,k-1))
#endif
                ExpAtt=EXP(-Att)
                tl_ExpAtt=-ExpAtt*tl_Att+                               &
#ifdef TL_IOMS
     &                    (1.0_r8+Att)*ExpAtt
#endif
                Itop=PAR
                tl_Itop=tl_PAR
                PAR=Itop*(1.0_r8-ExpAtt)/Att    ! average at center
                tl_PAR=(-tl_Att*PAR+tl_Itop*(1.0_r8-ExpAtt)-            &
     &                  Itop*tl_ExpAtt)/Att+                            &
#ifdef TL_IOMS
     &                 Itop/Att
#endif
                Light(i,k)=PAR
                tl_Light(i,k)=tl_PAR
!
!  Light attenuation at the bottom of the grid cell. It is the starting
!  PAR value for the next (deeper) vertical grid cell.
!
                PAR=Itop*ExpAtt
                tl_PAR=tl_Itop*ExpAtt+Itop*tl_ExpAtt-                   &
#ifdef TL_IOMS
     &                 PAR
#endif
              END DO
            ELSE                                       ! night time
              DO k=1,N(ng)
                Light(i,k)=0.0_r8
                tl_Light(i,k)=0.0_r8
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
              tl_cff=cff1*tl_Light(i,k)
              DO nb=1,NBands
                Bio(i,k,aCDOM(nb))=Bio(i,k,aCDOM(nb))/(1.0_r8+cff)
                tl_Bio(i,k,aCDOM(nb))=(tl_Bio(i,k,aCDOM(nb))-           &
     &                                 tl_cff*Bio(i,k,aCDOM(nb)))/      &
     &                                (1.0_r8+cff)+                     &
#ifdef TL_IOMS
     &                                cff*Bio(i,k,aCDOM(nb))/           &
     &                                (1.0_r8+cff)
#endif
              END DO
            END DO
          END DO

        END DO ITER_LOOP
!
!-----------------------------------------------------------------------
!  Update global tracer variables: Add increment due to BGC processes
!  to tracer array in time index "nnew". Index "nnew" is solution after
!  advection and mixing and has transport units (m Tunits) hence the
!  increment is multiplied by Hz.  Notice that we need to subtract
!  original values "Bio_old" at the top of the routine to just account
!  for the concentractions affected by BGC processes. This also takes
!  into account any constraints (non-negative concentrations, carbon
!  concentration range) specified before entering BGC kernel. If "Bio"
!  were unchanged by BGC processes, the increment would be exactly
!  zero. Notice that final tracer values, t(:,:,:,nnew,:) are not
!  bounded >=0 so that we can preserve total inventory of nutrients
!  when advection causes tracer concentration to go negative.
!-----------------------------------------------------------------------
!
        DO itrc=1,NBT
          ibio=idbio(itrc)
          DO k=1,N(ng)
            DO i=Istr,Iend
              cff=Bio(i,k,ibio)-Bio_old(i,k,ibio)
              tl_cff=tl_Bio(i,k,ibio)-tl_Bio_old(i,k,ibio)
!^            t(i,j,k,nnew,ibio)=t(i,j,k,nnew,ibio)+cff*Hz(i,j,k)
!^
              tl_t(i,j,k,nnew,ibio)=tl_t(i,j,k,nnew,ibio)+              &
     &                              tl_cff*Hz(i,j,k)+cff*tl_Hz(i,j,k)-  &
#ifdef TL_IOMS
     &                              cff*Hz(i,j,k)
#endif
            END DO
          END DO
        END DO

      END DO J_LOOP

      RETURN
      END SUBROUTINE rp_biology_tile
