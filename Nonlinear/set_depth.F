#include "cppdefs.h"
      MODULE set_depth_mod
#ifdef SOLVE3D
!
!========================================== Alexander F. Shchepetkin ===
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
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
      SUBROUTINE set_depth (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_coupling
      USE mod_grid
      USE mod_ocean
!
      integer, intent(in) :: ng, tile

# include "tile.h"
!
      CALL set_depth_tile (ng, Istr, Iend, Jstr, Jend,                  &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     OCEAN(ng) % zeta,                            &
     &                     GRID(ng) % h,                                &
# ifdef ICESHELF
     &                     GRID(ng) % zice,                             &
# endif
     &                     GRID(ng) % hinv,                             &
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
     &                           zeta, h,                               &
# ifdef ICESHELF
     &                           zice,                                  &
# endif
     &                           hinv, Zt_avg1,                         &
     &                           Hz, z_r, z_w)
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

# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: zeta(LBi:,LBj:,:)
      real(r8), intent(in) :: h(LBi:,LBj:)
#  ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:,LBj:)
#  endif
      real(r8), intent(inout) :: hinv(LBi:,LBj:)
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
      real(r8), intent(inout) :: hinv(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: Zt_avg1(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(out) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(out) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
# endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j, k

      real(r8) :: cff_r, cff1_r, cff2_r, cff_w, cff1_w, cff2_w
      real(r8) :: z_r0, z_w0

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute time evolving depths and vertical thicknesses.
!-----------------------------------------------------------------------
!
!  At initialization, load initial free-surface field into fast-time
!  averaged array.
!
      IF (iic(ng).eq.0) THEN
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            Zt_avg1(i,j)=zeta(i,j,1)
            hinv(i,j)=1.0_r8/h(i,j)
          END DO
        END DO
# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
        CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,             &
     &                          LBi, UBi, LBj, UBj, hinv)
# endif
      END IF
!
!  Compute vertical depths (meters, negative) at RHO- and W-points,
!  and vertical grid thicknesses.
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          z_w(i,j,0)=-h(i,j)
# ifdef ICESHELF
          z_w(i,j,0)=z_w(i,j,0)+zice(i,j)
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
# ifdef ICESHELF
            z_w0=cff_w+cff1_w*h(i,j)+zice(i,j)
            z_w(i,j,k)=z_w0+                                            &
     &                 Zt_avg1(i,j)*(1.0_r8+(z_w0-zice(i,j))*hinv(i,j))
            z_r0=cff_r+cff1_r*h(i,j)+zice(i,j)
            z_r(i,j,k)=z_r0+                                            &
     &                 Zt_avg1(i,j)*(1.0_r8+(z_r0-zice(i,j))*hinv(i,j))
# else
            z_w0=cff_w+cff1_w*h(i,j)
            z_w(i,j,k)=z_w0+Zt_avg1(i,j)*(1.0_r8+z_w0*hinv(i,j))
            z_r0=cff_r+cff1_r*h(i,j)
            z_r(i,j,k)=z_r0+Zt_avg1(i,j)*(1.0_r8+z_r0*hinv(i,j))
# endif
            Hz(i,j,k)=z_w(i,j,k)-z_w(i,j,k-1)
          END DO
        END DO
      END DO
# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      CALL exchange_w3d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, 0, N(ng), z_w)
      CALL exchange_r3d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, 1, N(ng), z_r)
      CALL exchange_r3d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, 1, N(ng), Hz)
# endif
      RETURN
      END SUBROUTINE set_depth_tile
#endif
      END MODULE set_depth_mod
