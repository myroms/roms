#include "cppdefs.h"
      MODULE bvf_mix_mod
#if defined NONLINEAR && defined BVF_MIXING
!
!=======================================================================
!  Copyright (c) 2005 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine computes tracer vertical mixing as a function of the   !
!  Brunt-Vaisala frequency.  The vertical mixing of momentum  is set   !
!  to its background value.  If static unstable regime, the vertical   !
!  mixing is set to "bvf_nu0c".                                        !
!                                                                      !
!======================================================================!
!
      implicit none

      PRIVATE
      PUBLIC  :: bvf_mix

      CONTAINS
!
!***********************************************************************
      SUBROUTINE bvf_mix (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_mixing

      integer, intent(in) :: ng, tile

# include "tile.h"

      CALL bvf_mix_tile (ng, Istr, Iend, Jstr, Jend,                    &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   MIXING(ng) % bvf,                              &
     &                   MIXING(ng) % AKt,                              &
     &                   MIXING(ng) % AKv)
      RETURN
      END SUBROUTINE bvf_mix
!
!***********************************************************************
      SUBROUTINE bvf_mix_tile (ng, Istr, Iend, Jstr, Jend,              &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         bvf, Akt, Akv)
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
      real(r8), intent(in) :: bvf(LBi:,LBj:,0:)
      real(r8), intent(out) :: Akt(LBi:,LBj:,0:,:)
      real(r8), intent(out) :: Akv(LBi:,LBj:,0:)
# else
      real(r8), intent(in) :: bvf(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(out) :: Akt(LBi:UBi,LBj:UBj,0:N(ng),NAT)
      real(r8), intent(out) :: Akv(LBi:UBi,LBj:UBj,0:N(ng))
# endif
!
!  Local variable declarations.
!
      integer :: i, itrc, j, k

      real(r8) :: cff
!
!-----------------------------------------------------------------------
!  Set tracer diffusivity as function of the Brunt-vaisala frequency.
!  Set vertical viscosity to its background value.
!-----------------------------------------------------------------------
!
      DO k=1,N(ng)-1
        DO j=Jstr,Jend
          DO i=Istr,Iend
            Akv(i,j,k)=Akv_bak(ng)
            IF (bvf(i,j,k).lt.0.0_r8) THEN
              Akv(i,j,k)=bvf_nu0c
              Akt(i,j,k,itemp)=bvf_nu0c
# ifdef SALINITY
              Akt(i,j,k,isalt)=bvf_nu0c
# endif
            ELSE IF (bvf(i,j,k).eq.0.0_r8) THEN
              Akv(i,j,k)=Akv_bak(ng)
              Akt(i,j,k,itemp)=Akt_bak(itemp,ng)
# ifdef SALINITY
              Akt(i,j,k,isalt)=Akt_bak(isalt,ng)
# endif
            ELSE
              cff=bvf_nu0/SQRT(bvf(i,j,k))
              Akt(i,j,k,itemp)=MIN(bvf_numax,MAX(bvf_numin,cff))
              Akv(i,j,k)=Akt(i,j,k,itemp)
# ifdef SALINITY
              Akt(i,j,k,isalt)=Akt(i,j,k,itemp)
# endif
            END IF
          END DO
        END DO
      END DO
# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
!
! Set periodic boundary conditions.
!
      CALL exchange_w3d_tile (ng, iNLM, Istr, Iend, Jstr, Jend,         &
     &                        LBi, UBi, LBj, UBj, 0, N(ng),             &
     &                        NghostPoints,                             &
     &                        Akv)
      DO itrc=1,NAT
        CALL exchange_w3d_tile (ng, iNLM, Istr, Iend, Jstr, Jend,       &
     &                          LBi, UBi, LBj, UBj, 0, N(ng),           &
     &                          NghostPoints,                           &
     &                          Akt(:,:,:,itemp))
      END DO
# endif
      RETURN
      END SUBROUTINE bvf_mix_tile
#endif
      END MODULE bvf_mix_mod
