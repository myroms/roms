#include "cppdefs.h"
#undef NEUMANN
#undef LINEAR_CONTINUATION
#undef REMIX_BED

      MODULE sediment_mod
#if defined NONLINEAR && defined SEDIMENT
!
!svn $Id$
!=======================================================================
!  Copyright (c) 2002-2007 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license           Hernan G. Arango   !
!    See License_ROMS.txt                   Alexander F. Shchepetkin   !
!==================================================== John C. Warner ===
!                                                                      !
!  This  routine computes the sediment sources and sinks and adds      !
!  then the global sediment tracer fields. Currently, it includes      !
!  the following:                                                      !
!                                                                      !
!  * Vertical settling of sediment in the water column.                !
!  * Erosive and depositional flux interactions of sediment            !
!    between water column and the bed.                                 !
!  * Transport of multiple grain sizes.                                !
!  * Bed layer stratigraphy.                                           !
!  * Seawater/sediment vertical level distribution:                    !
!                                                                      !
!         W-level  RHO-level                                           !
!                                                                      !
!            N     _________                                           !
!                 |         |                                          !
!                 |    N    |                                          !
!          N-1    |_________|  S                                       !
!                 |         |  E                                       !
!                 |   N-1   |  A                                       !
!            2    |_________|  W                                       !
!                 |         |  A                                       !
!                 |    2    |  T                                       !
!            1    |_________|  E                                       !
!                 |         |  R                                       !
!                 |    1    |                                          !
!            0    |_________|_____ bathymetry                          !
!                 |/////////|                                          !
!                 |    1    |                                          !
!            1    |_________|  S                                       !
!                 |         |  E                                       !
!                 |    2    |  D                                       !
!            2    |_________|  I                                       !
!                 |         |  M                                       !
!                 |  Nbed-1 |  E                                       !
!        Nbed-1   |_________|  N                                       !
!                 |         |  T                                       !
!                 |  Nbed   |                                          !
!         Nbed    |_________|                                          !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: sediment

      CONTAINS
!
!***********************************************************************
      SUBROUTINE sediment (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_forces
      USE mod_grid
      USE mod_ocean
      USE mod_stepping
# ifdef BBL_MODEL
      USE mod_bbl
# endif
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
      CALL wclock_on (ng, iNLM, 16)
# endif
      CALL sediment_tile (ng, Istr, Iend, Jstr, Jend,                   &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    nstp(ng), nnew(ng),                           &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
# ifdef MASKING
     &                    GRID(ng) % rmask,                             &
     &                    GRID(ng) % umask,                             &
     &                    GRID(ng) % vmask,                             &
# endif
     &                    GRID(ng) % Hz,                                &
     &                    GRID(ng) % z_w,                               &
# ifdef BBL_MODEL
     &                    BBL(ng) % bustrc,                             &
     &                    BBL(ng) % bvstrc,                             &
     &                    BBL(ng) % bustrcwmax,                         &
     &                    BBL(ng) % bvstrcwmax,                         &
# endif
     &                    FORCES(ng) % bustr,                           &
     &                    FORCES(ng) % bvstr,                           &
     &                    OCEAN(ng) % t,                                &
# ifdef BEDLOAD
     &                    GRID(ng) % h,                                 &
     &                    GRID(ng) % om_r,                              &
     &                    GRID(ng) % on_r,                              &
     &                    OCEAN(ng) % bedldu,                           &
     &                    OCEAN(ng) % bedldv,                           &
# endif
     &                    OCEAN(ng) % bed,                              &
     &                    OCEAN(ng) % bed_frac,                         &
     &                    OCEAN(ng) % bed_mass,                         &
     &                    OCEAN(ng) % bottom)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 16)
# endif
      RETURN
      END SUBROUTINE sediment
!
!***********************************************************************
      SUBROUTINE sediment_tile (ng, Istr, Iend, Jstr, Jend,             &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          nstp, nnew,                             &
     &                          pm, pn,                                 &
# ifdef MASKING
     &                          rmask, umask, vmask,                    &
# endif
     &                          Hz, z_w,                                &
# ifdef BBL_MODEL
     &                          bustrc, bvstrc,                         &
     &                          bustrcwmax, bvstrcwmax,                 &
# endif
     &                          bustr, bvstr,                           &
     &                          t,                                      &
# ifdef BEDLOAD
     &                          h, om_r, on_r,                          &
     &                          bedldu, bedldv,                         &
# endif
     &                          bed, bed_frac, bed_mass,                &
     &                          bottom)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_sediment
!
      USE bc_3d_mod, ONLY : bc_r3d_tile
# ifdef BEDLOAD
#  if defined EW_PERIODIC || defined NS_PERIODIC
      USE exchange_2d_mod, ONLY : exchange_u2d_tile, exchange_v2d_tile
#  endif
# endif
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange3d, mp_exchange4d
# endif
# ifdef TS_MPDATA
#  if defined EW_PERIODIC || defined NS_PERIODIC
      USE exchange_3d_mod, ONLY : exchange_r3d_tile
#  endif
      USE t3dbc_mod, ONLY : t3dbc_tile
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: nstp, nnew
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
#  ifdef BBL_MODEL
      real(r8), intent(in) :: bustrc(LBi:,LBj:)
      real(r8), intent(in) :: bvstrc(LBi:,LBj:)
      real(r8), intent(in) :: bustrcwmax(LBi:,LBj:)
      real(r8), intent(in) :: bvstrcwmax(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: bustr(LBi:,LBj:)
      real(r8), intent(in) :: bvstr(LBi:,LBj:)
#  ifdef BEDLOAD
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: om_r(LBi:,LBj:)
      real(r8), intent(in) :: on_r(LBi:,LBj:)
      real(r8), intent(inout) :: bedldu(LBi:,LBj:,:)
      real(r8), intent(inout) :: bedldv(LBi:,LBj:,:)
#  endif
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:) 
      real(r8), intent(inout) :: bed(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: bed_frac(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: bed_mass(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: bottom(LBi:,LBj:,:)
# else
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
#  ifdef BBL_MODEL
      real(r8), intent(in) :: bustrc(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: bvstrc(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: bustrcwmax(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: bvstrcwmax(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: bustr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: bvstr(LBi:UBi,LBj:UBj)
#  ifdef BEDLOAD
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_r(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_r(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: bedldu(LBi:UBi,LBj:UBj,NST)
      real(r8), intent(inout) :: bedldv(LBi:UBi,LBj:UBj,NST)
#  endif
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: bed(LBi:UBi,LBj:UBj,Nbed,MBEDP)
      real(r8), intent(inout) :: bed_frac(LBi:UBi,LBj:UBj,Nbed,NST)
      real(r8), intent(inout) :: bed_mass(LBi:UBi,LBj:UBj,Nbed,1:2,NST)
      real(r8), intent(inout) :: bottom(LBi:UBi,LBj:UBj,MBOTP)
# endif
!
!  Local variable declarations.
!
# ifdef DISTRIBUTE
#  ifdef EW_PERIODIC
      logical :: EWperiodic=.TRUE.
#  else
      logical :: EWperiodic=.FALSE.
#  endif
#  ifdef NS_PERIODIC
      logical :: NSperiodic=.TRUE.
#  else
      logical :: NSperiodic=.FALSE.
#  endif
# endif
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: Ksed, i, indx, ised, j, k, ks
      integer :: bnew

      real(r8), parameter :: eps = 1.0E-14_r8

      real(r8) :: cff, cff1, cff2, cff3, cffL, cffR, dltL, dltR
      real(r8) :: cu, ero_flux, cff4
      real(r8) :: thck_avail, thck_to_add

      integer, dimension(PRIVATE_1D_SCRATCH_ARRAY,N(ng)) :: ksource

      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N(ng)) :: FC

      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,N(ng)) :: Hz_inv
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,N(ng)) :: Hz_inv2
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,N(ng)) :: Hz_inv3
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,N(ng)) :: qc
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,N(ng)) :: qR
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,N(ng)) :: qL
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,N(ng)) :: WR
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,N(ng)) :: WL

      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,NST) :: dep_mass
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: tau_w

# ifdef BEDLOAD
      real(r8) :: a_slopex, a_slopey, angleu, anglev
      real(r8) :: bedld, bedld_mass
      real(r8) :: dzdx, dzdy, orsmgd, slope_coeff
      real(r8) :: smgd, smgdr, Umag

      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: FX
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: FE
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: FX_r
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: FE_r
# endif

# include "set_bounds.h"

# ifdef BEDLOAD
      bnew=nnew
# else
      bnew=nstp
# endif
!
!-----------------------------------------------------------------------
! Compute bottom stresses.
!-----------------------------------------------------------------------
!
# ifdef BBL_MODEL
      DO j=Jstr-1,Jend+1
        DO i=Istr-1,Iend+1
          tau_w(i,j)=SQRT(bustrcwmax(i,j)*bustrcwmax(i,j)+              &
     &                    bvstrcwmax(i,j)*bvstrcwmax(i,j))
        END DO
      END DO
# else
#  ifdef EW_PERIODIC
#   define I_RANGE Istr-1,Iend+1
#  else
#   define I_RANGE MAX(Istr-1,1),MIN(Iend+1,Lm(ng))
#  endif
#  ifdef NS_PERIODIC
#   define J_RANGE Jstr-1,Jend+1
#  else
#   define J_RANGE MAX(Jstr-1,1),MIN(Jend+1,Mm(ng))
#  endif
      DO i=I_RANGE
        DO j=J_RANGE
          tau_w(i,j)=0.5_r8*SQRT((bustr(i,j)+bustr(i+1,j))*             &
     &                           (bustr(i,j)+bustr(i+1,j))+             &
     &                           (bvstr(i,j)+bvstr(i,j+1))*             &
     &                           (bvstr(i,j)+bvstr(i,j+1)))
        END DO
      END DO
#  undef I_RANGE
#  undef J_RANGE
# endif

!
# ifdef BEDLOAD
!
!-----------------------------------------------------------------------
!  Compute bedload sediment transport.
!-----------------------------------------------------------------------
!
#  ifdef EW_PERIODIC
#   define I_RANGE Istr-1,Iend+1
#  else
#   define I_RANGE MAX(Istr-1,1),MIN(Iend+1,Lm(ng))
#  endif
#  ifdef NS_PERIODIC
#   define J_RANGE Jstr-1,Jend+1
#  else
#   define J_RANGE MAX(Jstr-1,1),MIN(Jend+1,Mm(ng))
#  endif

      DO ised=1,NST
        smgd=(Srho(ised,ng)/rho0-1.0_r8)*g*Sd50(ised,ng)
        smgdr=SQRT(smgd)*Sd50(ised,ng)*Srho(ised,ng)
        orsmgd=1.0_r8/(rho0*smgd)
!
        DO j=J_RANGE
          DO i=I_RANGE
!
! Compute stress components at rho points.
!
#  ifdef BBL_MODEL
            cff1=bustrc(i,j)
            cff2=bvstrc(i,j)
            Umag=SQRT(cff1*cff1+cff2*cff2)+eps
            angleu=cff1/Umag
            anglev=cff2/Umag
#  else
            cff1=0.5_r8*(bustr(i,j)+bustr(i+1,j))
            cff2=0.5_r8*(bvstr(i,j)+bvstr(i,j+1))
            Umag=SQRT(cff1*cff1+cff2*cff2)+eps
            angleu=cff1/Umag
            anglev=cff2/Umag
#  endif
!
! Magnitude of bed load at rho points. Meyer-Peter Muller formulation.
! bedld has dimensions of kg m-1 s-1
!
            bedld=8.0_r8*(MAX((tau_w(i,j)*rho0*orsmgd-0.047_r8),        &
     &                        0.0_r8)**1.5_r8)*smgdr
!
! Partition bedld into xi and eta directions, still at rho points.
! FX_r and FE_r have dimensions of kg.
!
            FX_r(i,j)=angleu*bedld*on_r(i,j)*dt(ng)
            FE_r(i,j)=anglev*bedld*om_r(i,j)*dt(ng)
!
! Correct for along-direction slope.
!
            cff1=0.5_r8*(1.0_r8+SIGN(1.0_r8,FX_r(i,j)))
            cff2=0.5_r8*(1.0_r8-SIGN(1.0_r8,FX_r(i,j)))
            dzdx=cff1*(h(i+1,j)-h(i,j))*0.5_r8*(pm(i+1,j)+pm(i,j))+     &
     &           cff2*(h(i,j)-h(i-1,j))*0.5_r8*(pm(i-1,j)+pm(i,j))
            cff1=0.5_r8*(1.0_r8+SIGN(1.0_r8,FE_r(i,j)))
            cff2=0.5_r8*(1.0_r8-SIGN(1.0_r8,FE_r(i,j)))
            dzdy=cff1*(h(i,j+1)-h(i,j))*0.5_r8*(pn(i,j+1)+pn(i,j))+     &
     &           cff2*(h(i,j)-h(i,j-1))*0.5_r8*(pn(i,j-1)+pn(i,j))
            slope_coeff=1.0_r8
            a_slopex=1.0_r8+slope_coeff*                                &
     &               (TAN(30.0_r8*pi/180.0_r8)/                         &
     &               (COS(ATAN(dzdx))*(tan(30.0_r8*pi/180.0_r8)-        &
     &                dzdx))-1.0_r8)
            a_slopey=1.0_r8+slope_coeff*                                &
     &               (TAN(30.0_r8*pi/180.0_r8)/                         &
     &               (COS(ATAN(dzdy))*(tan(30.0_r8*pi/180.0_r8)-        &
     &                dzdy))-1.0_r8)
            FX_r(i,j)=FX_r(i,j)*a_slopex
            FE_r(i,j)=FE_r(i,j)*a_slopey
!
            bedld_mass=ABS(FX_r(i,j))+ABS(FE_r(i,j))+eps
!
! Limit bed load to available bed mass.
!
            FX_r(i,j)=MIN(ABS(FX_r(i,j)),                               &
     &                    bed_mass(i,j,1,nstp,ised)*                    &
     &                    om_r(i,j)*on_r(i,j)*ABS(FX_r(i,j))/           &
     &                    bedld_mass)*                                  &
     &                SIGN(1.0_r8,FX_r(i,j))
            FE_r(i,j)=MIN(ABS(FE_r(i,j)),                               &
     &                    bed_mass(i,j,1,nstp,ised)*                    &
     &                    om_r(i,j)*on_r(i,j)*ABS(FE_r(i,j))/           &
     &                    bedld_mass)*                                  &
     &                SIGN(1.0_r8,FE_r(i,j))
          END DO
        END DO
#  ifndef EW_PERIODIC
        IF (WESTERN_EDGE) THEN
          DO j=J_RANGE
            FX_r(Istr-1,j)=FX_r(Istr,j)
            FE_r(Istr-1,j)=FE_r(Istr,j)
          END DO
        END IF
        IF (EASTERN_EDGE) THEN
          DO j=J_RANGE
            FX_r(Iend+1,j)=FX_r(Iend,j)
            FE_r(Iend+1,j)=FE_r(Iend,j)
          END DO
        END IF
#  endif
#  ifndef NS_PERIODIC
        IF (SOUTHERN_EDGE) THEN
          DO i=I_RANGE
            FX_r(i,Jstr-1)=FX_r(i,Jstr)
            FE_r(i,Jstr-1)=FE_r(i,Jstr)
          END DO
        END IF
        IF (NORTHERN_EDGE) THEN
          DO i=I_RANGE
            FX_r(i,Jend+1)=FX_r(i,Jend)
            FE_r(i,Jend+1)=FE_r(i,Jend)
          END DO
        END IF
#  endif
#  undef I_RANGE
#  undef J_RANGE
#  if !defined EW_PERIODIC && !defined NS_PERIODIC
        IF ((SOUTHERN_EDGE).and.(WESTERN_EDGE)) THEN
          FX_r(Istr-1,Jstr-1)=0.5_r8*(FX_r(Istr  ,Jstr-1)+              &
     &                                FX_r(Istr-1,Jstr  ))
          FE_r(Istr-1,Jstr-1)=0.5_r8*(FE_r(Istr  ,Jstr-1)+              &
     &                                FE_r(Istr-1,Jstr  ))
        END IF
        IF ((SOUTHERN_EDGE).and.(EASTERN_EDGE)) THEN
          FX_r(Iend+1,Jstr-1)=0.5_r8*(FX_r(Iend  ,Jstr-1)+              &
     &                                FX_r(Iend+1,Jstr  ))
          FE_r(Iend+1,Jstr-1)=0.5_r8*(FE_r(Iend  ,Jstr-1)+              &
     &                                FE_r(Iend+1,Jstr  ))
        END IF
        IF ((NORTHERN_EDGE).and.(WESTERN_EDGE)) THEN
          FX_r(Istr-1,Jend+1)=0.5_r8*(FX_r(Istr-1,Jend  )+              &
     &                                FX_r(Istr  ,Jend+1))
          FE_r(Istr-1,Jend+1)=0.5_r8*(FE_r(Istr-1,Jend  )+              &
     &                                FE_r(Istr  ,Jend+1))
        END IF
        IF ((NORTHERN_EDGE).and.(EASTERN_EDGE)) THEN
          FX_r(Iend+1,Jend+1)=0.5_r8*(FX_r(Iend+1,Jend  )+              &
     &                                FX_r(Iend  ,Jend+1))
          FE_r(Iend+1,Jend+1)=0.5_r8*(FE_r(Iend+1,Jend  )+              &
     &                                FE_r(Iend  ,Jend+1))
        END IF
#  endif
!
! Upwind shift FX_r and FE_r to u and v points.
!
        DO j=Jstr-1,Jend+1
          DO i=Istr,Iend+1
            cff1=0.5_r8*(1.0_r8+SIGN(1.0_r8,FX_r(i,j)))
            cff2=0.5_r8*(1.0_r8-SIGN(1.0_r8,FX_r(i,j)))
            FX(i,j)=0.5_r8*(1.0_r8+SIGN(1.0_r8,FX_r(i-1,j)))*           &
     &              (cff1*FX_r(i-1,j)+                                  &
     &               cff2*0.5_r8*(FX_r(i-1,j)+FX_r(i,j)))+              &
     &              0.5_r8*(1.0_r8-SIGN(1.0_r8,FX_r(i-1,j)))*           &
     &              (cff2*FX_r(i  ,j)+                                  &
     &               cff1*0.5_r8*(FX_r(i-1,j)+FX_r(i,j)))
#  ifdef MASKING
            FX(i,j)=FX(i,j)*umask(i,j)
#  endif
          END DO
        END DO
        DO j=Jstr,Jend+1
          DO i=Istr-1,Iend+1
            cff1=0.5_r8*(1.0_r8+SIGN(1.0_r8,FE_r(i,j)))
            cff2=0.5_r8*(1.0_r8-SIGN(1.0_r8,FE_r(i,j)))
            FE(i,j)=0.5_r8*(1.0_r8+SIGN(1.0_r8,FE_r(i,j-1)))*           &
     &              (cff1*FE_r(i,j-1)+                                  &
     &               cff2*0.5_r8*(FE_r(i,j-1)+FE_r(i,j)))+              &
     &              0.5_r8*(1.0_r8-SIGN(1.0_r8,FE_r(i,j-1)))*           &
     &              (cff2*FE_r(i  ,j)+                                  &
     &               cff1*0.5_r8*(FE_r(i,j-1)+FE_r(i,j)))
#  ifdef MASKING
            FE(i,j)=FE(i,j)*vmask(i,j)
#  endif
          END DO
        END DO
!
# ifndef EW_PERIODIC
        IF (WESTERN_EDGE) THEN
          DO j=Jstr-1,Jend+1
#  ifdef WESTERN_WALL
            FX(Istr,j)=0.0_r8
#  endif
          END DO
        END IF
        IF (EASTERN_EDGE) THEN
          DO j=Jstr-1,Jend+1
#  ifdef EASTERN_WALL
            FX(Iend+1,j)=0.0_r8
#  endif
          END DO
        END IF
#  endif
#  ifndef NS_PERIODIC
        IF (SOUTHERN_EDGE) THEN
          DO i=Istr-1,Iend+1
#  ifdef SOUTHERN_WALL
            FE(i,Jstr)=0.0_r8
#  endif
          END DO
        END IF
        IF (NORTHERN_EDGE) THEN
          DO i=Istr-1,Iend+1
#  ifdef NORTHERN_WALL
            FE(i,Jend+1)=0.0_r8
#  endif
          END DO
        END IF
#  endif
!
!  Determine flux divergence and evaluate change in bed properties.
!
        DO j=Jstr,Jend
          DO i=Istr,Iend
            cff=(FX(i+1,j)-FX(i,j)+                                     &
     &           FE(i,j+1)-FE(i,j))*pm(i,j)*pn(i,j)
            bed_mass(i,j,1,nnew,ised)=MAX(bed_mass(i,j,1,nstp,ised)-    &
     &                                    cff,0.0_r8)
#if !defined SUSPLOAD
            DO k=2,Nbed
              bed_mass(i,j,k,nnew,ised)=bed_mass(i,j,k,nstp,ised)
            END DO
#endif
            bed(i,j,1,ithck)=MAX(bed(i,j,1,ithck)-                      &
     &                           cff/(Srho(ised,ng)*                    &
     &                                (1.0_r8-bed(i,j,1,iporo))),       &
     &                           0.0_r8)
#  ifdef MASKING
            bed(i,j,1,ithck)=bed(i,j,1,ithck)*rmask(i,j)
#  endif
          END DO
        END DO
!
!-----------------------------------------------------------------------
!  Output bedload fluxes.
!-----------------------------------------------------------------------
!
        cff=0.5_r8/dt(ng)
        DO j=JstrR,JendR
          DO i=Istr,IendR
            bedldu(i,j,ised)=FX(i,j)*(pn(i-1,j)+pn(i,j))*cff
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            bedldv(i,j,ised)=FE(i,j)*(pm(i,j-1)+pm(i,j))*cff
          END DO
        END DO
      END DO
!
!  Update mean surface properties.
!  Sd50 must be positive definite, due to BBL routines.
!  Srho must be >1000, due to (s-1) in BBL routines.
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          cff3=0.0_r8
          DO ised=1,NST
            cff3=cff3+bed_mass(i,j,1,nnew,ised)
          END DO
          DO ised=1,NST
            bed_frac(i,j,1,ised)=bed_mass(i,j,1,nnew,ised)/             &
     &                           MAX(cff3,eps)
          END DO
!
          cff1=1.0_r8
          cff2=1.0_r8
          cff3=1.0_r8
          cff4=1.0_r8
          DO ised=1,NST
            cff1=cff1*tau_ce(ised,ng)**bed_frac(i,j,1,ised)
            cff2=cff2*Sd50(ised,ng)**bed_frac(i,j,1,ised)
            cff3=cff3*(wsed(ised,ng)+eps)**bed_frac(i,j,1,ised)
            cff4=cff4*Srho(ised,ng)**bed_frac(i,j,1,ised)
          END DO
          bottom(i,j,itauc)=cff1
          bottom(i,j,isd50)=cff2
          bottom(i,j,iwsed)=cff3
          bottom(i,j,idens)=MAX(cff4,1050.0_r8)
        END DO
      END DO
# endif
!
# ifdef SUSPLOAD
!
!-----------------------------------------------------------------------
!  Update bed layer stratigraphy.
!-----------------------------------------------------------------------
!
#  ifdef REMIX_BED
      DO j=Jstr,Jend
        DO i=Istr,Iend 
!
!  Calculate active layer thickness, bottom(i,j,iactv).
!
          bottom(i,j,iactv)=MAX(0.0_r8,                                 &
     &                          0.007_r8*                               &
     &                          (tau_w(i,j)-bottom(i,j,itauc))*rho0)+   &
     &                      6.0_r8*bottom(i,j,isd50)
!
!  Ensure top bed layer thickness is greater or equal than the active
!  layer thickness. If need to add sediment to top layer, then entrain
!  from lower levels. Create new layers at bottom to maintain Nbed.
!
          IF (bottom(i,j,iactv).gt.bed(i,j,1,ithck)) THEN
            IF (Nbed.eq.1) THEN
              bottom(i,j,iactv)=bed(i,j,1,ithck)
            ELSE
              thck_to_add=bottom(i,j,iactv)-bed(i,j,1,ithck)
              thck_avail=0.0_r8
              Ksed=1                                         ! initialize
              DO k=2,Nbed
                IF (thck_avail.lt.thck_to_add) THEN
                  thck_avail=thck_avail+bed(i,j,k,ithck)
                  Ksed=k
                END IF
              END DO
!
!  Catch here if there was not enough bed material.
!
              IF (thck_avail.lt.thck_to_add) THEN
                bottom(i,j,iactv)=bed(i,j,1,ithck)+thck_avail
                thck_to_add=thck_avail
              END IF
!
!  Upate bed mass of top layer and fractional layer.
!
              cff2=MAX(thck_avail-thck_to_add,0.0_r8)/                  &
     &             MAX(bed(i,j,Ksed,ithck),eps)
              DO ised=1,NST
                cff1=0.0_r8 
                DO k=1,Ksed
                  cff1=cff1+bed_mass(i,j,k,ised)
                END DO
                bed_mass(i,j,1   ,ised)=cff1-                           &
     &                                  bed_mass(i,j,Ksed,ised)*cff2
                bed_mass(i,j,Ksed,ised)=bed_mass(i,j,Ksed,ised)*cff2
              END DO
!
!  Update thickness of fractional layer Ksed.
!
              bed(i,j,Ksed,ithck)=MAX(thck_avail-thck_to_add,0.0_r8)
!
!  Upate bed fraction of top layer.
!
              cff3=0.0_r8
              DO ised=1,NST
                cff3=cff3+bed_mass(i,j,1,ised)
              END DO
              DO ised=1,NST
                bed_frac(i,j,1,ised)=bed_mass(i,j,1,ised)/              &
     &                               MAX(cff3,eps)
              END DO
!
!  Update bed thickness of top layer.
!
              bed(i,j,1,ithck)=bottom(i,j,iactv)
!
!  Pull all layers closer to the surface.
!
              DO k=Ksed,Nbed
                ks=Ksed-2
                bed(i,j,k-ks,ithck)=bed(i,j,k,ithck)
                bed(i,j,k-ks,iporo)=bed(i,j,k,iporo)
                bed(i,j,k-ks,iaged)=bed(i,j,k,iaged)
                DO ised=1,NST
                  bed_frac(i,j,k-ks,ised)=bed_frac(i,j,k,ised)
                  bed_mass(i,j,k-ks,ised)=bed_mass(i,j,k,ised)
                END DO
              END DO
!
!  Add new layers onto the bottom. Split what was in the bottom layer to
!  fill these new empty cells. ("ks" is the number of new layers).
!
              ks=Ksed-2
              cff=1.0_r8/REAL(ks+1,r8)
              DO k=Nbed,Nbed-ks,-1
                bed(i,j,k,ithck)=bed(i,j,Nbed-ks,ithck)*cff
!!              bed(i,j,k,iporo)=bed(i,j,Nbed-ks,iporo)
                bed(i,j,k,iaged)=bed(i,j,Nbed-ks,iaged)
                DO ised=1,NST
                  bed_frac(i,j,k,ised)=bed_frac(i,j,Nbed-ks,ised)
                  bed_mass(i,j,k,ised)=bed_mass(i,j,Nbed-ks,ised)*cff
                END DO
              END DO
            END IF                             ! Nbed > 1
          END IF                               ! increase top bed layer
        END DO
      END DO
#  endif
!
!-----------------------------------------------------------------------
!  Add sediment Source/Sink terms.
!-----------------------------------------------------------------------
!
!  Compute inverse thicknessed to avoid repeated divisions.
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
!  Copy concentration of suspended sediment into scratch array "qc"
!  (q-central, restrict it to be positive) which is hereafter
!  interpreted as a set of grid-box averaged values for sediment
!  concentration.
!
        SED_LOOP: DO ised=1,NST
          indx=idsed(ised)
          DO k=1,N(ng)
            DO i=Istr,Iend
              qc(i,k)=t(i,j,k,nnew,indx)*Hz_inv(i,k)
            END DO
          END DO
!
!-----------------------------------------------------------------------
!  Vertical sinking of suspended sediment.
!-----------------------------------------------------------------------
!
!  Reconstruct vertical profile of suspended sediment "qc" in terms
!  of a set of parabolic segments within each grid box. Then, compute
!  semi-Lagrangian flux due to sinking.
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
!  Compute right and left side values (qR,qL) of parabolic segments
!  within grid box Hz(k); (WR,WL) are measures of quadratic variations. 
!
!  NOTE: Although each parabolic segment is monotonic within its grid
!        box, monotonicity of the whole profile is not guaranteed,
!        because qL(k+1)-qR(k) may still have different sign than
!        qc(k+1)-qc(k).  This possibility is excluded, after qL and qR
!        are reconciled using WENO procedure.
!
              cff=(dltR-dltL)*Hz_inv3(i,k)
              dltR=dltR-cff*Hz(i,j,k+1)
              dltL=dltL+cff*Hz(i,j,k-1)
              qR(i,k)=qc(i,k)+dltR
              qL(i,k)=qc(i,k)-dltL
              WR(i,k)=(2.0_r8*dltR-dltL)**2
              WL(i,k)=(dltR-2.0_r8*dltL)**2
            END DO
          END DO
          cff=1.0E-14_r8
          DO k=2,N(ng)-2
            DO i=Istr,Iend
              dltL=MAX(cff,WL(i,k  ))
              dltR=MAX(cff,WR(i,k+1))
              qR(i,k)=(dltR*qR(i,k)+dltL*qL(i,k+1))/(dltR+dltL)
              qL(i,k+1)=qR(i,k)
            END DO
          END DO
          DO i=Istr,Iend
            FC(i,N(ng))=0.0_r8              ! no-flux boundary condition
# if defined LINEAR_CONTINUATION
            qL(i,N(ng))=qR(i,N(ng)-1)
            qR(i,N(ng))=2.0_r8*qc(i,N(ng))-qL(i,N(ng))
# elif defined NEUMANN
            qL(i,N(ng))=qR(i,N(ng)-1)
            qR(i,N(ng))=1.5*qc(i,N(ng))-0.5_r8*qL(i,N(ng))
# else
            qR(i,N(ng))=qc(i,N(ng))         ! default strictly monotonic
            qL(i,N(ng))=qc(i,N(ng))         ! conditions
            qR(i,N(ng)-1)=qc(i,N(ng))
# endif
# if defined LINEAR_CONTINUATION 
            qR(i,1)=qL(i,2)
            qL(i,1)=2.0_r8*qc(i,1)-qR(i,1)
# elif defined NEUMANN
            qR(i,1)=qL(i,2)
            qL(i,1)=1.5_r8*qc(i,1)-0.5_r8*qR(i,1)
# else  
            qL(i,2)=qc(i,1)                 ! bottom grid boxes are
            qR(i,1)=qc(i,1)                 ! re-assumed to be
            qL(i,1)=qc(i,1)                 ! piecewise constant.
# endif
          END DO
!
!  Apply monotonicity constraint again, since the reconciled interfacial
!  values may cause a non-monotonic behavior of the parabolic segments
!  inside the grid box.
!
          DO k=1,N(ng)
            DO i=Istr,Iend
              dltR=qR(i,k)-qc(i,k)
              dltL=qc(i,k)-qL(i,k)
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
              qR(i,k)=qc(i,k)+dltR
              qL(i,k)=qc(i,k)-dltL
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
          cff=dt(ng)*ABS(Wsed(ised,ng))
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
              FC(i,k-1)=FC(i,k-1)+                                      &
     &                  Hz(i,j,ks)*cu*                                  &
     &                  (qL(i,ks)+                                      &
     &                   cu*(0.5_r8*(qR(i,ks)-qL(i,ks))-                &
     &                       (1.5_r8-cu)*                               &
     &                       (qR(i,ks)+qL(i,ks)-2.0_r8*qc(i,ks))))
            END DO
          END DO
          DO k=1,N(ng)
            DO i=Istr,Iend
              qc(i,k)=qc(i,k)+(FC(i,k)-FC(i,k-1))*Hz_inv(i,k)
            END DO
          END DO
!
!-----------------------------------------------------------------------
!  Sediment deposition and resuspension near the bottom.
!-----------------------------------------------------------------------
!
!  The deposition and resuspension of sediment on the bottom "bed"
!  is due to precepitation flux FC(:,0), already computed, and the
!  resuspension (erosion, hence called ero_flux). The resuspension is
!  applied to the bottom-most grid box value qc(:,1) so the total mass
!  is conserved. Restrict "ero_flux" so that "bed" cannot go negative
!  after both fluxes are applied.
!
          cff=1.0_r8/tau_ce(ised,ng)
          DO i=Istr,Iend
            dep_mass(i,ised)=0.0_r8
!
!  Compute erosion, ero_flux (kg/m2).
!
            cff1=(1.0_r8-bed(i,j,1,iporo))*bed_frac(i,j,1,ised)
            cff2=dt(ng)*Erate(ised,ng)*cff1
            cff3=Srho(ised,ng)*cff1
            ero_flux=MIN(MAX(0.0_r8,cff2*(cff*tau_w(i,j)-1.0_r8)),      &
     &                   MIN(cff3*bottom(i,j,iactv),                    &
     &                   bed_mass(i,j,1,bnew,ised))+FC(i,0))
!
!  Check if depositional.
!
            IF ((ero_flux-FC(i,0)).lt.0.0_r8) THEN
!
!  If first time step of deposit, then store deposit material in
!  temporary array, dep_mass.
!
              IF ((time(ng).gt.(bed(i,j,1,iaged)+1.1_r8*dt(ng))).and.   &
                  (bed(i,j,1,ithck).gt.0.005_r8))THEN
                dep_mass(i,ised)=-(ero_flux-FC(i,0))+eps
              ELSE
!
!  If it is not first time step of deposit, update bed thickness of
!  top layer.
!  
                bed(i,j,1,ithck)=MAX(bed(i,j,1,ithck)-                  &
     &                               (ero_flux-FC(i,0))/                &
     &                               (Srho(ised,ng)*                    &
     &                               (1.0_r8-bed(i,j,1,iporo))),0.0_r8)
# ifdef MASKING
                bed(i,j,1,ithck)=bed(i,j,1,ithck)*rmask(i,j)
# endif 
!
!  Upate bed mass of top layer.
!
                bed_mass(i,j,1,nnew,ised)=MAX(bed_mass(i,j,1,bnew,ised)-&
          &                                   (ero_flux-FC(i,0)),0.0_r8)
                bed(i,j,1,iaged)=time(ng)
              END IF
!
!  Else, if erosional.
!
            ELSE
!
!  Update bed thickness of top layer.
!  
              bed(i,j,1,ithck)=MAX(bed(i,j,1,ithck)-                    &
     &                             (ero_flux-FC(i,0))/                  &
     &                             (Srho(ised,ng)*                      &
     &                             (1.0_r8-bed(i,j,1,iporo))),0.0_r8)
# ifdef MASKING
              bed(i,j,1,ithck)=bed(i,j,1,ithck)*rmask(i,j)
# endif 
!
!  Upate bed mass of top layer.
!
              bed_mass(i,j,1,nnew,ised)=MAX(bed_mass(i,j,1,bnew,ised)-  &
     &                                      (ero_flux-FC(i,0)),0.0_r8)
            END IF
!
!  Compute new sedidiment concentrations.
!
            qc(i,1)=qc(i,1)+ero_flux*Hz_inv(i,1)
!
!  Update bed mass arrays.
!
            DO k=2,Nbed
              bed_mass(i,j,k,nnew,ised)=bed_mass(i,j,k,nstp,ised)
            END DO
          END DO
!
!-----------------------------------------------------------------------
!  Update global tracer variables (m Tunits).
!-----------------------------------------------------------------------
!
          DO k=1,N(ng)
            DO i=Istr,Iend
              t(i,j,k,nnew,indx)=Hz(i,j,k)*qc(i,k)
# ifdef TS_MPDATA
              t(i,j,k,3,indx)=Hz(i,j,k)*qc(i,k)
# endif
            END DO
          END DO
        END DO SED_LOOP
!
!  Upate bed fraction of top layer.
!
        DO i=Istr,Iend
          cff3=0.0_r8
          DO ised=1,NST
             cff3=cff3+bed_mass(i,j,1,nnew,ised)
          END DO
          DO ised=1,NST
            bed_frac(i,j,1,ised)=bed_mass(i,j,1,nnew,ised)/MAX(cff3,eps)
          END DO
        END DO
!
!  If first time step of deposit, create new layer and combine bottom
!  two bed layers.
!
          DO i=Istr,Iend
            cff=0.0_r8
!
!  Determine if deposition ocurred here.
!
            DO ised=1,NST
               cff=cff+dep_mass(i,ised)
            END DO
            IF (cff.gt.0.0_r8) THEN
              IF (NBED.gt.1) THEN
!
!  Combine bottom layers.
!
                bed(i,j,Nbed,ithck)=bed(i,j,Nbed-1,ithck)+              &
                                    bed(i,j,Nbed  ,ithck)
                bed(i,j,Nbed,iporo)=0.5_r8*(bed(i,j,Nbed-1,iporo)+      &
     &                                      bed(i,j,Nbed,iporo))
                bed(i,j,Nbed,iaged)=0.5_r8*(bed(i,j,Nbed-1,iaged)+      &
     &                                      bed(i,j,Nbed,iaged))
                cff3=0.0_r8
                DO ised=1,NST
                  bed_mass(i,j,Nbed,nnew,ised)=                         &
     &                               bed_mass(i,j,Nbed-1,nnew,ised)+    &
     &                               bed_mass(i,j,Nbed  ,nnew,ised)
                  cff3=cff3+bed_mass(i,j,Nbed,nnew,ised)
                END DO
                cff3=1.0_r8/MAX(cff3,eps)
                DO ised=1,NST
                  bed_frac(i,j,Nbed,ised)=bed_mass(i,j,Nbed,nnew,ised)* &
     &                                    cff3
                END DO
!
!  Push layers down.
!
                DO k=Nbed-1,2,-1
                  bed(i,j,k,ithck)=bed(i,j,k-1,ithck)
                  bed(i,j,k,iporo)=bed(i,j,k-1,iporo)
                  bed(i,j,k,iaged)=bed(i,j,k-1,iaged)
                  DO ised =1,NST
                    bed_frac(i,j,k,ised)=bed_frac(i,j,k-1,ised)
                    bed_mass(i,j,k,nnew,ised)=bed_mass(i,j,k-1,nnew,ised)
                  END DO
                END DO
!
!  Set new top layer parameters.
!
                bed(i,j,1,ithck)=0.0_r8
                DO ised=1,NST
                  bed_mass(i,j,1,nnew,ised)=0.0_r8
                END DO
              END IF !NBED=1
              cff3=0.0_r8
              DO ised=1,NST
                bed_mass(i,j,1,nnew,ised)=bed_mass(i,j,1,nnew,ised)+    &
     &                                    dep_mass(i,ised)
                cff3=cff3+bed_mass(i,j,1,nnew,ised)
                bed(i,j,1,ithck)=bed(i,j,1,ithck)+                      &
     &                           dep_mass(i,ised)/                      &
     &                           (Srho(ised,ng)*                        &
     &                            (1.0_r8-bed(i,j,1,iporo)))
              END DO
              bed(i,j,1,iaged)=time(ng)
              DO ised=1,NST
                bed_frac(i,j,1,ised)=bed_mass(i,j,1,nnew,ised)/         &
     &                               MAX(cff3,eps)
!
!  Reset deposition mass.
!
                dep_mass(i,ised)=0.0_r8
              END DO
            END IF
!
          END DO
      END DO J_LOOP
# endif
!
!  End of Suspended Sediment only section.
!
!  Determine top layer to be at least active layer thickness.
!
!  Ensure top bed layer thickness is greater or equal than active layer
!  thickness. If need to add sed to top layer, then entrain from lower
!  levels. Create new layers at bottom to maintain Nbed.
!
      J_LOOP2 : DO j=Jstr,Jend
        DO i=Istr,Iend
!
!  Calculate active layer thickness, bottom(i,j,iactv).
!
          bottom(i,j,iactv)=MAX(0.0_r8,                                 &
     &                          0.007_r8*                               &
     &                          (tau_w(i,j)-bottom(i,j,itauc))*rho0)+   &
     &                      6.0_r8*bottom(i,j,isd50)
!
          IF (bottom(i,j,iactv).gt.bed(i,j,1,ithck)) THEN
            IF (Nbed.eq.1) THEN
              bottom(i,j,iactv)=bed(i,j,1,ithck)
            ELSE
              thck_to_add=bottom(i,j,iactv)-bed(i,j,1,ithck)
              thck_avail=0.0_r8
              Ksed=1                                        ! initialize
              DO k=2,Nbed
                IF (thck_avail.lt.thck_to_add) THEN
                  thck_avail=thck_avail+bed(i,j,k,ithck)
                  Ksed=k
                END IF
              END DO
!
!  Catch here if there was not enough bed material.
!
              IF (thck_avail.lt.thck_to_add) THEN
                bottom(i,j,iactv)=bed(i,j,1,ithck)+thck_avail
                thck_to_add=thck_avail
              END IF
!
!  Update bed mass of top layer and fractional layer.
!
              cff2=MAX(thck_avail-thck_to_add,0.0_r8)/                  &
     &             MAX(bed(i,j,Ksed,ithck),eps)
              DO ised=1,NST
                cff1=0.0_r8
                DO k=1,Ksed
                  cff1=cff1+bed_mass(i,j,k,nnew,ised)
                END DO
                bed_mass(i,j,1   ,nnew,ised)=cff1-                      &
     &                                bed_mass(i,j,Ksed,nnew,ised)*cff2
                bed_mass(i,j,Ksed,nnew,ised)=                           &
     &                                bed_mass(i,j,Ksed,nnew,ised)*cff2
              END DO
!
!  Update thickness of fractional layer ksource_sed.
!
              bed(i,j,Ksed,ithck)=MAX(thck_avail-thck_to_add,0.0_r8)
!
!  Upate bed fraction of top layer.
!
              cff3=0.0_r8
              DO ised=1,NST
                cff3=cff3+bed_mass(i,j,1,nnew,ised)
              END DO
              DO ised=1,NST
                bed_frac(i,j,1,ised)=bed_mass(i,j,1,nnew,ised)/         &
     &                               MAX(cff3,eps)
              END DO
!
!  Upate bed thickness of top layer.
!
              bed(i,j,1,ithck)=bottom(i,j,iactv)
!
!  Pull all layers closer to the surface.
!
              DO k=Ksed,Nbed
                ks=Ksed-2
                bed(i,j,k-ks,ithck)=bed(i,j,k,ithck)
                bed(i,j,k-ks,iporo)=bed(i,j,k,iporo)
                bed(i,j,k-ks,iaged)=bed(i,j,k,iaged)
                DO ised=1,NST
                  bed_frac(i,j,k-ks,ised)=bed_frac(i,j,k,ised)
                  bed_mass(i,j,k-ks,nnew,ised)=bed_mass(i,j,k,nnew,ised)
                END DO
              END DO
!
!  Add new layers onto the bottom. Split what was in the bottom layer to
!  fill these new empty cells. ("ks" is the number of new layers).
!
              ks=Ksed-2
              cff=1.0_r8/REAL(ks+1,r8)
              DO k=Nbed,Nbed-ks,-1
                bed(i,j,k,ithck)=bed(i,j,Nbed-ks,ithck)*cff
                bed(i,j,k,iaged)=bed(i,j,Nbed-ks,iaged)
                DO ised=1,NST
                  bed_frac(i,j,k,ised)=bed_frac(i,j,Nbed-ks,ised)
                  bed_mass(i,j,k,nnew,ised)=                            &
     &                             bed_mass(i,j,Nbed-ks,nnew,ised)*cff
                END DO
              END DO
            END IF  ! Nbed > 1
          END IF  ! increase top bed layer
!
!  Update mean surface properties.
!  Sd50 must be positive definite, due to BBL routines.
!  Srho must be >1000, due to (s-1) in BBL routines
!
          cff1=1.0_r8
          cff2=1.0_r8
          cff3=1.0_r8
          cff4=1.0_r8
          DO ised=1,NST
            cff1=cff1*tau_ce(ised,ng)**bed_frac(i,j,1,ised)
            cff2=cff2*Sd50(ised,ng)**bed_frac(i,j,1,ised)
            cff3=cff3*(wsed(ised,ng)+eps)**bed_frac(i,j,1,ised)
            cff4=cff4*Srho(ised,ng)**bed_frac(i,j,1,ised)
          END DO
          bottom(i,j,itauc)=cff1
          bottom(i,j,isd50)=cff2
          bottom(i,j,iwsed)=cff3
          bottom(i,j,idens)=MAX(cff4,1050.0_r8)
        END DO
      END DO J_LOOP2
!
!-----------------------------------------------------------------------
!  Apply periodic or gradient boundary conditions to property arrays.
!-----------------------------------------------------------------------
!
      DO ised=1,NST
        CALL bc_r3d_tile (ng, Istr, Iend, Jstr, Jend,                   &
     &                    LBi, UBi, LBj, UBj, 1, Nbed,                  &
     &                    bed_frac(:,:,:,ised))
        CALL bc_r3d_tile (ng, Istr, Iend, Jstr, Jend,                   &
     &                    LBi, UBi, LBj, UBj, 1, Nbed,                  &
     &                    bed_mass(:,:,:,nnew,ised))
# if defined BEDLOAD && (defined EW_PERIODIC || defined NS_PERIODIC)
        CALL exchange_u2d_tile (ng, Istr, Iend, Jstr, Jend,             &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          bedldu(:,:,ised))
        CALL exchange_v2d_tile (ng, Istr, Iend, Jstr, Jend,             &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          bedldv(:,:,ised))
# endif
      END DO
# ifdef DISTRIBUTE
      CALL mp_exchange4d (ng, iNLM, 2, Istr, Iend, Jstr, Jend,          &
     &                    LBi, UBi, LBj, UBj, 1, Nbed, 1, NST,          &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    bed_frac,                                     &
     &                    bed_mass(:,:,:,nnew,:)
#  if defined BEDLOAD && (defined EW_PERIODIC || defined NS_PERIODIC)
      CALL mp_exchange3d (ng, iNLM, 2, Istr, Iend, Jstr, Jend,          &
     &                    LBi, UBi, LBj, UBj, 1, NST,                   &
     &                    bedldu, bedldv)
#  endif
# endif

      DO i=1,MBEDP
        CALL bc_r3d_tile (ng, Istr, Iend, Jstr, Jend,                   &
     &                    LBi, UBi, LBj, UBj, 1, Nbed,                  &
     &                    NghostPoints,                                 &
     &                    bed(:,:,:,i))
      END DO
# ifdef DISTRIBUTE
      CALL mp_exchange4d (ng, iNLM, 1, Istr, Iend, Jstr, Jend,          &
     &                    LBi, UBi, LBj, UBj, 1, Nbed, 1, MBEDP,        &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    bed)
# endif

      CALL bc_r3d_tile (ng, iNLM, Istr, Iend, Jstr, Jend,               &
     &                  LBi, UBi, LBj, UBj, 1, MBOTP,                   &
     &                  NghostPoints,                                   &
     &                  bottom)
# ifdef DISTRIBUTE
      CALL mp_exchange3d (ng, iNLM, 1, Istr, Iend, Jstr, Jend,          &
     &                    LBi, UBi, LBj, UBj, 1, MBOTP,                 &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    bottom)
# endif

      RETURN
      END SUBROUTINE sediment_tile
#endif
      END MODULE sediment_mod
