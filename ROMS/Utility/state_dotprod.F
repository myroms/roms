#include "cppdefs.h"

      SUBROUTINE state_dotprod (ng, model, Istr, Iend, Jstr, Jend,      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NstateVars, DotProd,                    &
#ifdef MASKING
     &                          rmask, umask, vmask,                    &
#endif
#ifdef SOLVE3D
     &                          s1_t, s2_t,                             &
     &                          s1_u, s2_u,                             &
     &                          s1_v, s2_v,                             &
#endif
     &                          s1_ubar, s2_ubar,                       &
     &                          s1_vbar, s2_vbar,                       &
     &                          s1_zeta, s2_zeta)
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2007 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes the dot product between two model states:     !
!                                                                      !
!      DotProd(0:NstateVars) = < s1, s2 >                              !
!                                                                      !
!  where                                                               !
!                                                                      !
!      DotProd(0)           All state variable dot product             !
!      DotProd(isFsur)      Free-surface contribution                  !
!      DotProd(isUbar)      2D U-momentum contribution                 !
!      DotProd(isVbar)      2D V-momentum contribution                 !
!      DotProd(isUvel)      3D U-momentum contribution                 !
!      DotProd(isVvel)      3D V-momentum contribution                 !
!      DotProd(isTvar(:))   Tracer-type variables contribution         !   
!                                                                      !
!-----------------------------------------------------------------------
!
      USE mod_param
      USE mod_parallel
      USE mod_ncparam
#ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_reduce
#endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: NstateVars
!
#ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#endif
#ifdef SOLVE3D
      real(r8), intent(in) :: s1_t(LBi:UBi,LBj:UBj,N(ng),NT(ng))
      real(r8), intent(in) :: s2_t(LBi:UBi,LBj:UBj,N(ng),NT(ng))
      real(r8), intent(in) :: s1_u(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: s2_u(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: s1_v(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: s2_v(LBi:UBi,LBj:UBj,N(ng))
#endif
      real(r8), intent(in) :: s1_ubar(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: s2_ubar(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: s1_vbar(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: s2_vbar(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: s1_zeta(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: s2_zeta(LBi:UBi,LBj:UBj)
!
      real(r8), intent(out), dimension(0:NstateVars) :: DotProd
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: NSUB, i, j
#ifdef SOLVE3D
      integer :: itrc, k
#endif
      real(r8) :: cff
      real(r8), dimension(0:NstateVars) :: my_DotProd
#ifdef DISTRIBUTE
      character (len=3), dimension(0:NstateVars) :: op_handle
#endif

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute dot product between S1 and S2 model state trajectories.
!-----------------------------------------------------------------------
!
      DO i=0,NstateVars
        my_DotProd(i)=0.0_r8
      END DO
!
!  2D state variables.
!
#ifndef SOLVE3D
      DO j=JstrR,JendR
        DO i=Istr,IendR
          cff=s1_ubar(i,j)*s2_ubar(i,j)
# ifdef MASKING
          cff=cff*umask(i,j)
# endif
          my_DotProd(0)=my_DotProd(0)+cff
          my_DotProd(isUbar)=my_DotProd(isUbar)+cff
        END DO
      END DO
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          cff=s1_vbar(i,j)*s2_vbar(i,j)
# ifdef MASKING
          cff=cff*vmask(i,j)
# endif
          my_DotProd(0)=my_DotProd(0)+cff
          my_DotProd(isVbar)=my_DotProd(isVbar)+cff
        END DO
      END DO
#endif
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          cff=s1_zeta(i,j)*s2_zeta(i,j)
#ifdef MASKING
          cff=cff*rmask(i,j)
#endif
          my_DotProd(0)=my_DotProd(0)+cff
          my_DotProd(isFsur)=my_DotProd(isFsur)+cff
        END DO
      END DO
#ifdef SOLVE3D
!
!  3D state variables.
!
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=Istr,IendR
            cff=s1_u(i,j,k)*s2_u(i,j,k)
# ifdef MASKING
            cff=cff*umask(i,j)
# endif
            my_DotProd(0)=my_DotProd(0)+cff
            my_DotProd(isUvel)=my_DotProd(isUvel)+cff
          END DO
        END DO
      END DO
      DO k=1,N(ng)
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            cff=s1_v(i,j,k)*s2_v(i,j,k)
# ifdef MASKING
            cff=cff*vmask(i,j)
# endif
            my_DotProd(0)=my_DotProd(0)+cff
            my_DotProd(isVvel)=my_DotProd(isVvel)+cff
          END DO
        END DO
      END DO
      DO itrc=1,NT(ng)
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              cff=s1_t(i,j,k,itrc)*s2_t(i,j,k,itrc)
# ifdef MASKING
              cff=cff*rmask(i,j)
# endif
              my_DotProd(0)=my_DotProd(0)+cff
              my_DotProd(isTvar(itrc))=my_DotProd(isTvar(itrc))+cff
            END DO
          END DO
        END DO
      END DO
#endif
!
!-----------------------------------------------------------------------
!  Perform parallel global reduction operations.
!-----------------------------------------------------------------------
!
      IF (SOUTH_WEST_CORNER.and.                                        &
     &    NORTH_EAST_CORNER) THEN
        NSUB=1                           ! non-tiled application
      ELSE
        NSUB=NtileX(ng)*NtileE(ng)       ! tiled application
      END IF
!$OMP CRITICAL (DOT_PROD)
      IF (tile_count.eq.0) THEN
        DO i=0,NstateVars
          DotProd(i)=0.0_r8
        END DO
      END IF
      DO i=0,NstateVars
        DotProd(i)=DotProd(i)+my_DotProd(i)
      END DO
      tile_count=tile_count+1
      IF (tile_count.eq.NSUB) THEN
        tile_count=0
#ifdef DISTRIBUTE
        DO i=0,NstateVars
          op_handle(i)='SUM'
        END DO
        CALL mp_reduce (ng, model, NstateVars+1, DotProd(0:),           &
     &                  op_handle(0:))
#endif
      END IF
!$OMP END CRITICAL (DOT_PROD)

      RETURN
      END SUBROUTINE state_dotprod
