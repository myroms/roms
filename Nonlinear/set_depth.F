#include "cppdefs.h"
      MODULE set_depth_mod
#ifdef SOLVE3D
!
!========================================== Alexander F. Shchepetkin ===
!  Copyright (c) 2005 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine computes the time evolving depths of the model grid    !
!  and its associated vertical transformation metric (thickness).      !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: set_depth, set_depth_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE set_depth (ng, tile, kini, FirstPass)
!***********************************************************************
!
      USE mod_param
      USE mod_coupling
      USE mod_grid
      USE mod_ocean
!
      integer, intent(in) :: ng, tile, kini
      logical, intent(in) :: FirstPass

# include "tile.h"
!
      CALL set_depth_tile (ng, Istr, Iend, Jstr, Jend,                  &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     kini, FirstPass,                             &
     &                     OCEAN(ng) % zeta,                            &
     &                     GRID(ng) % h,                                &
# ifdef ICESHELF
     &                     GRID(ng) % zice,                             &
# endif
# if defined SEDIMENT && defined SED_MORPH
     &                     OCEAN(ng) % bed,                             &
     &                     GRID(ng) % bed_thick0,                       &
# endif
     &                     COUPLING(ng) % Zt_avg1,                      &
     &                     GRID(ng) % Hz,                               &
     &                     GRID(ng) % z_r,                              &
     &                     GRID(ng) % z_w)
      RETURN
      END SUBROUTINE set_depth
!
!***********************************************************************
      SUBROUTINE set_depth_tile (ng, Istr, Iend, Jstr, Jend,            &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           kini, FirstPass,                       &
     &                           zeta, h,                               &
# ifdef ICESHELF
     &                           zice,                                  &
# endif
# if defined SEDIMENT && defined SED_MORPH
     &                           bed, bed_thick0,                       &
# endif
     &                           Zt_avg1,                               &
     &                           Hz, z_r, z_w)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
# if defined SEDIMENT && defined SED_MORPH
      USE mod_sediment
# endif
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
      integer, intent(in) :: kini
      logical, intent(in) :: FirstPass

# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: zeta(LBi:,LBj:,:)
      real(r8), intent(in) :: h(LBi:,LBj:)
#  ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:,LBj:)
#  endif
#  if defined SEDIMENT && defined SED_MORPH
      real(r8), intent(in) :: bed(LBi:,LBj:,:,:)
      real(r8), intent(inout):: bed_thick0(LBi:,LBj:)
#  endif
      real(r8), intent(inout) :: Zt_avg1(LBi:,LBj:)
      real(r8), intent(out) :: Hz(LBi:,LBj:,:)
      real(r8), intent(out) :: z_r(LBi:,LBj:,:)
      real(r8), intent(out) :: z_w(LBi:,LBj:,0:)
# else
      real(r8), intent(in) :: zeta(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
#  ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:UBi,LBj:UBj)
#  endif
#  if defined SEDIMENT && defined SED_MORPH

      real(r8), intent(in) :: bed(LBi:UBi,LBj:UBj,Nbed,MBEDP)
      real(r8), intent(inout):: bed_thick0(LBi:UBi,LBj:UBi)
#  endif
      real(r8), intent(inout) :: Zt_avg1(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(out) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(out) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
# endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j, k, kbed

      real(r8) :: cff_r, cff1_r, cff2_r, cff_w, cff1_w, cff2_w
      real(r8) :: hinv, hwater, z_r0, z_w0

# if defined SEDIMENT && defined SED_MORPH
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: BedThick
# endif

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute time evolving depths and vertical thicknesses.
!-----------------------------------------------------------------------
!
!  On first pass during initialization, load initial free-surface field
!  (kini index) into fast-time averaged array. If appropriate, also
!  compute initial sum of all bed layers thicknesses.
!
      IF (FirstPass) THEN
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            Zt_avg1(i,j)=zeta(i,j,kini)
# if defined SEDIMENT && defined SED_MORPH
            bed_thick0(i,j)=0.0_r8
            DO kbed=1,Nbed
              bed_thick0(i,j)=bed_thick0(i,j)+bed(i,j,kbed,ithck)
            END DO
# endif
          END DO
        END DO
# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
#  if defined SEDIMENT && defined SED_MORPH
        CALL exchange_r2d_tile (ng, iNLM, Istr, Iend, Jstr, Jend,       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          bed_thick0)
#  endif
# endif
      END IF
!
!  Compute vertical depths (meters, negative) at RHO- and W-points,
!  and vertical grid thicknesses.
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          z_w(i,j,0)=-h(i,j)
# if defined SEDIMENT && defined SED_MORPH
          BedThick(i,j)=0.0_r8
          DO kbed=1,Nbed
            BedThick(i,j)=BedThick(i,j)+bed(i,j,kbed,ithck)
          END DO
          BedThick(i,j)=BedThick(i,j)-bed_thick0(i,j)
          z_w(i,j,0)=z_w(i,j,0)+BedThick(i,j)
# endif
        END DO
        DO k=1,N(ng)
          cff_w=hc(ng)*(SCALARS(ng)%sc_w(k)-SCALARS(ng)%Cs_w(k))
          cff1_w=SCALARS(ng)%Cs_w(k)
          cff2_w=SCALARS(ng)%sc_w(k)+1.0_r8
          cff_r=hc(ng)*(SCALARS(ng)%sc_r(k)-SCALARS(ng)%Cs_r(k))
          cff1_r=SCALARS(ng)%Cs_r(k)
          cff2_r=SCALARS(ng)%sc_r(k)+1.0_r8
          DO i=IstrR,IendR
            hwater=h(i,j)
# ifdef ICESHELF
            hwater=hwater-ABS(zice(i,j))
# endif
# if defined SEDIMENT && defined SED_MORPH
            hwater=hwater-BedThick(i,j)
# endif
            hinv=1.0_r8/hwater
            z_w0=cff_w+cff1_w*hwater
            z_w(i,j,k)=z_w0+Zt_avg1(i,j)*(1.0_r8+z_w0*hinv)
            z_r0=cff_r+cff1_r*hwater
            z_r(i,j,k)=z_r0+Zt_avg1(i,j)*(1.0_r8+z_r0*hinv)
# ifdef ICESHELF
            z_w(i,j,k)=z_w(i,j,k)-ABS(zice(i,j))
            z_r(i,j,k)=z_r(i,j,k)-ABS(zice(i,j))
# endif
            Hz(i,j,k)=z_w(i,j,k)-z_w(i,j,k-1)
          END DO
        END DO
      END DO
# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
!
!  Exchange boundary information.
!
      CALL exchange_w3d_tile (ng, iNLM, Istr, Iend, Jstr, Jend,         &
     &                        LBi, UBi, LBj, UBj, 0, N(ng),             &
     &                        NghostPoints,                             &
     &                        z_w)
      CALL exchange_r3d_tile (ng, iNLM, Istr, Iend, Jstr, Jend,         &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        NghostPoints,                             &
     &                        z_r)
      CALL exchange_r3d_tile (ng, iNLM, Istr, Iend, Jstr, Jend,         &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        NghostPoints,                             &
     &                        Hz)
# endif
      RETURN
      END SUBROUTINE set_depth_tile
#endif
      END MODULE set_depth_mod
