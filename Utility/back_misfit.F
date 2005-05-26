#include "cppdefs.h"
      MODULE back_misfit_mod
#ifdef BACKGROUND
!
!=================================================== Andrew M. Moore ===
!  Copyright (c) ROMS/TOMS Adjoint Group                               !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine computes the model minus background adjoint misfit     !
!  forcing for each state variable.                                    !
!                                                                      !
!======================================================================!
!
      implicit none

      CONTAINS
!
!***********************************************************************
      SUBROUTINE back_misfit (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_ocean
      USE mod_stepping
!
      integer, intent(in) :: ng, tile

# include "tile.h"
!
      CALL back_misfit_tile (ng, Istr, Iend, Jstr, Jend,                &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       knew(ng),                                  &
# ifdef SOLVE3D
     &                       nstp(ng)                                   &
# endif
     &                       Lnew(ng), Lbck(ng),                        &
# ifdef SOLVE3D
     &                       OCEAN(ng) % u,                             &
     &                       OCEAN(ng) % v,                             &
     &                       OCEAN(ng) % t,                             &
     &                       OCEAN(ng) % b_u,                           &
     &                       OCEAN(ng) % b_v,                           &
     &                       OCEAN(ng) % b_t,                           &
     &                       OCEAN(ng) % ad_u,                          &
     &                       OCEAN(ng) % ad_v,                          &
     &                       OCEAN(ng) % ad_t,                          &
# endif
     &                       OCEAN(ng) % ubar,                          &
     &                       OCEAN(ng) % vbar,                          &
     &                       OCEAN(ng) % zeta,                          &
     &                       OCEAN(ng) % b_ubar,                        &
     &                       OCEAN(ng) % b_vbar,                        &
     &                       OCEAN(ng) % b_zeta,                        &
     &                       OCEAN(ng) % ad_ubar,                       &
     &                       OCEAN(ng) % ad_vbar,                       &
     &                       OCEAN(ng) % ad_zeta)
      RETURN
      END SUBROUTINE back_misfit
!
!***********************************************************************
      SUBROUTINE back_misfit_tile (ng, Istr, Iend, Jstr, Jend,          &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             knew,                                &
# ifdef SOLVE3D
     &                             nstp,                                &
# endif
     &                             Lnew, Lbck,                          &
# ifdef SOLVE3D
     &                             u, v, t,                             &
     &                             b_u, b_v, b_t,                       &
     &                             ad_u, ad_v, ad_t,                    &
# endif
     &                             ubar, vbar, zeta,                    &
     &                             b_ubar, b_vbar, b_zeta,              &
     &                             ad_ubar, ad_vbar, ad_zeta)
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
      integer, intent(in) :: knew
# ifdef SOLVE3D
      integer, intent(in) :: nstp
# endif
      integer, intent(in) :: Lnew, Lbck

# ifdef ASSUMED_SHAPE
#  ifdef SOLVE3D
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
      real(r8), intent(in) :: t(LBi:,LBj:,:,:,:)
      real(r8), intent(in) :: b_u(LBi:,LBj:,:)
      real(r8), intent(in) :: b_v(LBi:,LBj:,:)
      real(r8), intent(in) :: b_t(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_v(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_t(LBi:,LBj:,:,:,:)
#  endif
      real(r8), intent(in) :: ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: vbar(LBi:,LBj:,:)
      real(r8), intent(in) :: zeta(LBi:,LBj:,:)
      real(r8), intent(in) :: b_ubar(LBi:,LBj:)
      real(r8), intent(in) :: b_vbar(LBi:,LBj:)
      real(r8), intent(in) :: b_zeta(LBi:,LBj:)
      real(r8), intent(inout) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_vbar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_zeta(LBi:,LBj:,:)
# else
#  ifdef SOLVE3D
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(in) :: b_u(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: b_v(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: b_t(LBi:UBi,LBj:UBj,N(ng),NT(ng))
      real(r8), intent(inout) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
#  endif
      real(r8), intent(in) :: ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: zeta(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: b_ubar(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: b_vbar(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: b_zeta(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_zeta(LBi:UBi,LBj:UBj,3)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, kfrc
# ifdef SOLVE3D
      integer :: itrc, k, nfrc
# endif
!
# include "set_bounds.h"
!
!----------------------------------------------------------------------
!  Compute model minus background adjoint misfit forcing.
!----------------------------------------------------------------------
!
!  Set adjoint time indices to update.
!
      IF (iic(ng).eq.1) THEN
        kfrc=knew
      ELSE
        kfrc=1
      END IF
# ifdef SOLVE3D
      nfrc=nstp
# endif
!
!  Free-surface.
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ad_zeta(i,j,kfrc)=ad_zeta(i,j,kfrc)+                          &
     &                      2.0_r8*(zeta(i,j,Lnew)-zeta(i,j,Lbck))/     &
     &                      b_zeta(i,j)
        END DO
      END DO
!
!  2D momentum contribution.
!
      DO j=JstrR,JendR
        DO i=Istr,IendR
          ad_ubar(i,j,kfrc)=ad_ubar(i,j,kfrc)+                          &
     &                      2.0_r8*(ubar(i,j,Lnew)-ubar(i,j,Lbck))/     &
     &                      b_ubar(i,j)
        END DO
      END DO
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          ad_vbar(i,j,kfrc)=ad_vbar(i,j,kfrc)+                          &
     &                      2.0_r8*(vbar(i,j,Lnew)-vbar(i,j,Lbck))/     &
     &                      b_vbar(i,j)
        END DO
      END DO
# ifdef SOLVE3D
!
!  3D momentum contribution.
!
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=Istr,IendR
            ad_u(i,j,k,nfrc)=ad_u(i,j,k,nfrc)+                          &
     &                       2.0_r8*(u(i,j,k,Lnew)-u(i,j,k,Lbck))/      &
     &                       b_u(i,j,k)
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            ad_v(i,j,k,nfrc)=ad_v(i,j,k,nfrc)+                          &
     &                       2.0_r8*(v(i,j,k,Lnew)-v(i,j,k,Lbck))/      &
     &                       b_v(i,j,k)
          END DO
        END DO
      END DO
!
!  Tracers contribution.
!
      DO itrc=1,NT(ng)
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              ad_t(i,j,k,nfrc,itrc)=ad_t(i,j,k,nfrc,itrc)+              &
     &                              2.0_r8*(t(i,j,k,Lnew,itrc)-         &
     &                                      t(i,j,k,Lbck,itrc))/        &
     &                              b_t(i,j,k,itrc)
            END DO
          END DO
        END DO
      END DO
# endif
      IF (SOUTH_WEST_CORNER) THEN
        IF (Master) THEN
          WRITE (stdout,10) iic(ng)
  10      FORMAT (6x,'BACK_MISFIT  - Added background adjoint misfit ', &
     &            'forcing at iic = ',i8)
        END IF
      END IF

      RETURN
      END SUBROUTINE back_misfit_tile
#endif
      END MODULE back_misfit_mod

