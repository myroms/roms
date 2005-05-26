#include "cppdefs.h"
      MODULE back_cost_mod
#ifdef BACKGROUND
!
!=================================================== Andrew M. Moore ===
!  Copyright (c) ROMS/TOMS Adjoint Group                               !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine computes the background cost function (Jb) as the      !
!  misfit (squared difference) between model and background state:     !
!                                                                      !
!         Jb = transpose(X - Xb) * B^(-1) * (X -Xb)                    !
!                                                                      !
!  where                                                               !
!                                                                      !
!         Xb : background state (first guess)                          !
!         X  : model state                                             !
!         B  : background error covariance                             !
!                                                                      !
!  The  background  state  is usually a model forecast initialized     !
!  from a previous analysis.                                           !
!                                                                      !
!======================================================================!
!
      implicit none

      CONTAINS
!
!***********************************************************************
      SUBROUTINE back_cost (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_ocean
      USE mod_stepping
!
      integer, intent(in) :: ng, tile

# include "tile.h"
!
      CALL back_cost_tile (ng, Istr, Iend, Jstr, Jend,                  &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     Lnew(ng), Lbck(ng),                          &
# ifdef SOLVE3D
     &                     OCEAN(ng) % u,                               &
     &                     OCEAN(ng) % v,                               &
     &                     OCEAN(ng) % t,                               &
     &                     OCEAN(ng) % b_u,                             &
     &                     OCEAN(ng) % b_v,                             &
     &                     OCEAN(ng) % b_t,                             &
# endif
     &                     OCEAN(ng) % ubar,                            &
     &                     OCEAN(ng) % vbar,                            &
     &                     OCEAN(ng) % zeta,                            &
     &                     OCEAN(ng) % b_ubar,                          &
     &                     OCEAN(ng) % b_vbar,                          &
     &                     OCEAN(ng) % b_zeta)
      RETURN
      END SUBROUTINE back_cost
!
!***********************************************************************
      SUBROUTINE back_cost_tile (ng, Istr, Iend, Jstr, Jend,            &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           Lnew, Lbak,                            &
# ifdef SOLVE3D
     &                           u, v, t,                               &
     &                           b_u, b_v, b_t,                         &
# endif
     &                           ubar, vbar, zeta,                      &
     &                           b_ubar, b_vbar, b_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Lnew, Lbck

# ifdef ASSUMED_SHAPE
#  ifdef SOLVE3D
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
      real(r8), intent(in) :: t(LBi:,LBj:,:,:,:)
      real(r8), intent(in) :: b_u(LBi:,LBj:,:)
      real(r8), intent(in) :: b_v(LBi:,LBj:,:)
      real(r8), intent(in) :: b_t(LBi:,LBj:,:,:)
#  endif
      real(r8), intent(in) :: ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: vbar(LBi:,LBj:,:)
      real(r8), intent(in) :: zeta(LBi:,LBj:,:)
      real(r8), intent(in) :: b_ubar(LBi:,LBj:)
      real(r8), intent(in) :: b_vbar(LBi:,LBj:)
      real(r8), intent(in) :: b_zeta(LBi:,LBj:)
# else
#  ifdef SOLVE3D
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(in) :: b_u(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: b_v(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: b_t(LBi:UBi,LBj:UBj,N(ng),NT(ng))
#  endif
      real(r8), intent(in) :: ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: zeta(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: b_ubar(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: b_vbar(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: b_zeta(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: NSUB, i, j
# ifdef SOLVE3D
      integer :: itrc, k
# endif

      real(r8), dimension(0:NstateVar(ng)) :: my_BackFun

      real(r8) :: cff1, cff2

# ifdef DISTRIBUTE
      character (len=3), dimension(0:NstateVar(ng)) :: op_handle
# endif

# include "set_bounds.h"
!
!----------------------------------------------------------------------
!  Compute background cost function (Jb) as misfit between model and
!  background states at initial time of the assimilation window.
!  Initially, the misfit innovation matrix (X-Xb) is zero. As the
!  assimilation algorithm iterates, Jb becomes greater than zero.
!----------------------------------------------------------------------
!
!
      DO i=0,NstateVar(ng)
        my_BackFun(i)=0.0_r8
      END DO
!
!  Free-surface contribution.
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          cff1=zeta(i,j,Lnew)-zeta(i,j,Lbck)
          cff2=cff1*cff1/b_zeta(i,j)
          my_BackFun(0)=my_BackFun(0)+cff2
          my_BackFun(isFsur)=my_BackFun(isFsur)+cff2
        END DO
      END DO
!
!  2D momentum contribution.
!
      DO j=JstrR,JendR
        DO i=Istr,IendR
          cff1=ubar(i,j,Lnew)-ubar(i,j,Lbck)
          cff2=cff1*cff1/b_ubar(i,j)
          my_BackFun(0)=my_BackFun(0)+cff2
          my_BackFun(isUbar)=my_BackFun(isUbar)+cff2
        END DO
      END DO
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          cff1=vbar(i,j,Lnew)-vbar(i,j,Lbck)
          cff2=cff1*cff1/b_vbar(i,j)
          my_BackFun(0)=my_BackFun(0)+cff2
          my_BackFun(isVbar)=my_BackFun(isVbar)+cff2
        END DO
      END DO
# ifdef SOLVE3D
!
!  3D momentum contribution.
!
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=Istr,IendR
            cff1=u(i,j,k,Lnew)-u(i,j,k,Lbck)
            cff2=cff1*cff1/b_u(i,j,k)
            my_BackFun(0)=my_BackFun(0)+cff2
            my_BackFun(isUvel)=my_BackFun(isUvel)+cff2
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            cff1=v(i,j,k,Lnew)-v(i,j,k,Lbck)
            cff2=cff1*cff1/b_v(i,j,k)
            my_BackFun(0)=my_BackFun(0)+cff2
            my_BackFun(isVvel)=my_BackFun(isVvel)+cff2
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
              cff1=t(i,j,k,Lnew,itrc)-t(i,j,k,Lbck,itrc)
              cff2=cff1*cff1/b_t(i,j,k,itrc)
              my_BackFun(0)=my_BackFun(0)+cff2
              my_BackFun(isTvar(itrc))=my_BackFun(isTvar(itrc))+cff2
            END DO
          END DO
        END DO
      END DO
# endif
!
! Compute global background cost function.
!
      IF (SOUTH_WEST_CORNER.and.                                        &
     &    NORTH_EAST_CORNER) THEN
        NSUB=1                           ! non-tiled application
      ELSE
        NSUB=NtileX(ng)*NtileE(ng)       ! tiled application
      END IF
!$OMP CRITICAL (BACKCOST)
      IF (tile_count.eq.0) THEN
        DO i=0,NstateVar(ng)
          FOURDVAR(ng)%BackFun(i)=0.0_r8
        END DO  
      END IF
      DO i=0,NstateVar(ng)
        FOURDVAR(ng)%BackFun(i)=FOURDVAR(ng)%BackFun(i)+my_BackFun(i)
      END DO
      tile_count=tile_count+1
      IF (tile_count.eq.NSUB) THEN
        tile_count=0
# ifdef DISTRIBUTE
        DO i=0,NstateVar(ng)
          op_handle(i)='SUM'
        END DO 
        CALL mp_reduce (ng, iNLM, NstateVar(ng)+1,                      &
     &                  FOURDVAR(ng)%BackFun(0:),  op_handle(0:))
# endif
      END IF
!$OMP END CRITICAL (BACKCOST)

      RETURN
      END SUBROUTINE back_cost_tile
#endif
      END MODULE back_cost_mod
