#include "cppdefs.h"
      MODULE forcing_mod

#if defined W4DPSAS || defined NLM_OUTER
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2007 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine is used to force the nonlinear state equations in      !
!  weak constraint variational data assimilation.                      !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: forcing

      CONTAINS
!
!***********************************************************************
      SUBROUTINE forcing (ng, tile, Kfrc, Nfrc)
!***********************************************************************
!
      USE mod_param
      USE mod_ocean
# ifdef SOLVE3D
      USE mod_coupling
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, Kfrc, Nfrc
!
!  Local variable declarations.
!
# include "tile.h"
!
      CALL forcing_tile (ng, Istr, Iend, Jstr, Jend,                    &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   Kfrc, Nfrc,                                    &
# ifdef SOLVE3D
     &                   OCEAN(ng) % f_t,                               &
     &                   OCEAN(ng) % f_u,                               &
     &                   OCEAN(ng) % f_v,                               &
# else
     &                   OCEAN(ng) % f_ubar,                            &
     &                   OCEAN(ng) % f_vbar,                            &
# endif
     &                   OCEAN(ng) % f_zeta,                            &
# ifdef SOLVE3D
     &                   OCEAN(ng) % t,                                 &
     &                   OCEAN(ng) % u,                                 &
     &                   OCEAN(ng) % v,                                 &
# else
     &                   OCEAN(ng) % ubar,                              &
     &                   OCEAN(ng) % vbar,                              &
# endif
# ifdef SOLVE3D
     &                   COUPLING(ng) % Zt_avg1)
# else
     &                   OCEAN(ng) % zeta)
# endif

      RETURN
      END SUBROUTINE forcing
!
!***********************************************************************
      SUBROUTINE forcing_tile (ng, Istr, Iend, Jstr, Jend,              &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         Kfrc, Nfrc,                              &
# ifdef SOLVE3D
     &                         f_t, f_u, f_v,                           &
# else
     &                         f_ubar, f_vbar,                          &
# endif
     &                         f_zeta,                                  &
# ifdef SOLVE3D
     &                         t, u, v,                                 &
# else
     &                         ubar, vbar,                              &
# endif
# ifdef SOLVE3D
     &                         Zt_avg1)
# else
     &                         zeta)
# endif
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Kfrc
# ifdef SOLVE3D
      integer, intent(in) :: Nfrc
# endif
!
# ifdef ASSUMED_SHAPE
#  ifdef SOLVE3D
      real(r8), intent(in) :: f_t(LBi:,LBj:,:,:)
      real(r8), intent(in) :: f_u(LBi:,LBj:,:)
      real(r8), intent(in) :: f_v(LBi:,LBj:,:)
#  else
      real(r8), intent(in) :: f_ubar(LBi:,LBj:)
      real(r8), intent(in) :: f_vbar(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: f_zeta(LBi:,LBj:)
#  ifdef SOLVE3D
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: v(LBi:,LBj:,:,:)
#  else
      real(r8), intent(inout) :: ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: vbar(LBi:,LBj:,:)
#  endif
#  ifdef SOLVE3D
      real(r8), intent(inout) :: Zt_avg1(LBi:,LBj:)
#  else
      real(r8), intent(inout) :: zeta(LBi:,LBj:,:)
#  endif
# else
#  ifdef SOLVE3D
      real(r8), intent(in) :: t(LBi:UBi,LBj:UBj,N(ng),NT(ng))
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng))
#  else
      real(r8), intent(in) :: ubar(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vbar(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: zeta(LBi:UBi,LBj:UBj)
#  ifdef SOLVE3D
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: v(LBi:UBi,LBj:UBj,N(ng),2)
#  else
      real(r8), intent(inout) :: ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: vbar(LBi:UBi,LBj:UBj,3)
#  endif
#  ifdef SOLVE3D
      real(r8), intent(inout) :: Zt_avg1(LBi:UBi,LBj:UBj)
#  else
      real(r8), intent(inout) :: zeta(LBi:UBi,LBj:UBj,3)
#  endif
# endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j
# ifdef SOLVE3D
      integer :: itrc, k
# endif

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Add adjoint impulse forcing to nonlinear linear state.
!-----------------------------------------------------------------------
!
      IF (SOUTH_WEST_CORNER) THEN
        IF (Master) WRITE (stdout,10) tdays(ng)
 10     FORMAT (3x,'NL_FORCING   - added convolved adjoint impulse,',   &
     &          t64,'t = ', f12.4)
      END IF
!
!  Free-surface.
!
# ifdef SOLVE3D
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          Zt_avg1(i,j)=Zt_avg1(i,j)+f_zeta(i,j)
        END DO
      END DO
# else
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          zeta(i,j,Kfrc)=zeta(i,j,Kfrc)+f_zeta(i,j)
        END DO
      END DO
# endif

# ifndef SOLVE3D
!
!  2D momentum.
!
      DO j=JstrR,JendR
        DO i=Istr,IendR
          ubar(i,j,Kfrc)=ubar(i,j,Kfrc)+f_ubar(i,j)
        END DO
      END DO
!
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          vbar(i,j,Kfrc)=vbar(i,j,Kfrc)+f_vbar(i,j)
        END DO
      END DO

# else
!
!  3D momentum.
!
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=Istr,IendR
            u(i,j,k,Nfrc)=u(i,j,k,Nfrc)+f_u(i,j,k)
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            v(i,j,k,Nfrc)=v(i,j,k,Nfrc)+f_v(i,j,k)
          END DO
        END DO
      END DO
!
!  Tracers.
!
      DO itrc=1,NT(ng)
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              t(i,j,k,Nfrc,itrc)=t(i,j,k,Nfrc,itrc)+                    &
     &                           f_t(i,j,k,itrc)
            END DO
          END DO
        END DO
      END DO
# endif

      RETURN
      END SUBROUTINE forcing_tile
#endif
      END MODULE forcing_mod