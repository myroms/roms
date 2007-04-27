#include "cppdefs.h"

      SUBROUTINE state_addition (ng, Istr, Iend, Jstr, Jend,            &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           Linp1, Linp2, Lout,                    &
     &                           fac1, fac2,                            &
#ifdef MASKING
     &                           rmask, umask, vmask,                   &
#endif
#ifdef ADJUST_WSTRESS
     &                           s1_sustr, s2_sustr,                    &
     &                           s1_svstr, s2_svstr,                    &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                           s1_stflx, s2_stflx,                    &
# endif
     &                           s1_t, s2_t,                            &
     &                           s1_u, s2_u,                            &
     &                           s1_v, s2_v,                            &
#else
     &                           s1_ubar, s2_ubar,                      &
     &                           s1_vbar, s2_vbar,                      &
#endif
     &                           s1_zeta, s2_zeta)
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2007 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes the following model state addition:           !
!                                                                      !
!      s1_var(...,Lout) = fac1 * s1_var(...,Linp1) +                   !
!                         fac2 * s2_var(...,Linp2)                     !
!                                                                      !
!  where alpha and beta are scalars.                                   !
!                                                                      !
!-----------------------------------------------------------------------
!
      USE mod_param
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Linp1, Linp2, Lout
!
      real(r8), intent(in) :: fac1, fac2
!
#ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#endif
#ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: s2_sustr(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: s2_svstr(LBi:UBi,LBj:UBj,2)
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: s2_stflx(LBi:UBi,LBj:UBj,2,NT(ng))
# endif
      real(r8), intent(inout) :: s2_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: s2_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: s2_v(LBi:UBi,LBj:UBj,N(ng),2)
#else
      real(r8), intent(inout) :: s2_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: s2_vbar(LBi:UBi,LBj:UBj,3)
#endif
      real(r8), intent(inout) :: s2_zeta(LBi:UBi,LBj:UBj,3)

#ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: s1_sustr(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: s1_svstr(LBi:UBi,LBj:UBj,2)
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: s1_stflx(LBi:UBi,LBj:UBj,2,NT(ng))
# endif
      real(r8), intent(inout) :: s1_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: s1_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: s1_v(LBi:UBi,LBj:UBj,N(ng),2)
#else
      real(r8), intent(inout) :: s1_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: s1_vbar(LBi:UBi,LBj:UBj,3)
#endif
      real(r8), intent(inout) :: s1_zeta(LBi:UBi,LBj:UBj,3)
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j
#ifdef SOLVE3D
      integer :: itrc, k
#endif

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute the following operation between S1 and S2 model state
!  trajectories:
!                 S1(Lout) = fac1 * S1(Linp1) + fac2 * S2(Linp2)
!-----------------------------------------------------------------------
!
!  Free-surface.
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          s1_zeta(i,j,Lout)=fac1*s1_zeta(i,j,Linp1)+                    &
     &                      fac2*s2_zeta(i,j,Linp2)
#ifdef MASKING
          s1_zeta(i,j,Lout)=s1_zeta(i,j,Lout)*rmask(i,j)
#endif
        END DO
      END DO

#ifndef SOLVE3D
!
!  2D momentum.
!
      DO j=JstrR,JendR
        DO i=Istr,IendR
          s1_ubar(i,j,Lout)=fac1*s1_ubar(i,j,Linp1)+                    &
     &                      fac2*s2_ubar(i,j,Linp2)
# ifdef MASKING
          s1_ubar(i,j,Lout)=s1_ubar(i,j,Lout)*umask(i,j)
# endif
        END DO
      END DO
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          s1_vbar(i,j,Lout)=fac1*s1_vbar(i,j,Linp1)+                    &
     &                      fac2*s2_vbar(i,j,Linp2)
# ifdef MASKING
          s1_vbar(i,j,Lout)=s1_vbar(i,j,Lout)*vmask(i,j)
# endif
        END DO
      END DO
#endif

#ifdef ADJUST_STRESS
!
!  Surface momentum stress.
!
      DO j=JstrR,JendR
        DO i=Istr,IendR
          s1_sustr(i,j,Lout)=fac1*s1_sustr(i,j,Linp1)+                  &
     &                       fac2*s2_sustr(i,j,Linp2)
# ifdef MASKING
          s1_sustr(i,j,Lout)=s1_sustr(i,j,Lout)*umask(i,j)
# endif
        END DO
      END DO
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          s1_svstr(i,j,Lout)=fac1*s1_svstr(i,j,Linp1)+                  &
     &                       fac2*s2_svstr(i,j,Linp2)
# ifdef MASKING
          s1_svstr(i,j,Lout)=s1_svstr(i,j,Lout)*vmask(i,j)
# endif
        END DO
      END DO
#endif

#ifdef SOLVE3D
!
!  3D momentum.
!
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=Istr,IendR
            s1_u(i,j,k,Lout)=fac1*s1_u(i,j,k,Linp1)+                    &
     &                       fac2*s2_u(i,j,k,Linp2)
# ifdef MASKING
            s1_u(i,j,k,Lout)=s1_u(i,j,k,Lout)*umask(i,j)
# endif
          END DO
        END DO
      END DO
      DO k=1,N(ng)
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            s1_v(i,j,k,Lout)=fac1*s1_v(i,j,k,Linp1)+                    &
     &                       fac2*s2_v(i,j,k,Linp2)
# ifdef MASKING
            s1_v(i,j,k,Lout)=s1_v(i,j,k,Lout)*vmask(i,j)
# endif
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
              s1_t(i,j,k,Lout,itrc)=fac1*s1_t(i,j,k,Linp1,itrc)+        &
     &                              fac2*s2_t(i,j,k,Linp2,itrc)
# ifdef MASKING
              s1_t(i,j,k,Lout,itrc)=s1_t(i,j,k,Lout,itrc)*rmask(i,j)
# endif
            END DO
          END DO
        END DO
      END DO

# ifdef ADJUST_STFLUX
!
!  Surface tracers flux.
!
      DO itrc=1,NT(ng)
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            s1_stflx(i,j,Lout,itrc)=fac1*s1_stflx(i,j,Linp1,itrc)+      &
     &                              fac2*s2_stflx(i,j,Linp2,itrc)
#  ifdef MASKING
            s1_stflx(i,j,Lout,itrc)=s1_t(i,j,Lout,itrc)*rmask(i,j)
#  endif
          END DO
        END DO
      END DO
# endif
#endif

      RETURN
      END SUBROUTINE state_addition