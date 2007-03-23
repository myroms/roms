!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2007 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module minimizes a quadratic cost function using the conjugate !
!  gradient algorithm proposed by Mike Fisher (ECMWF).                 !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    Fisher, M., 1997: Efficient Minimization of Quadratic Penalty     !
!      funtions, unpublish manuscript, 1-14.                           !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC :: cgradient

      CONTAINS
!
!***********************************************************************
      SUBROUTINE cgradient (ng, tile, model, Iter)
!***********************************************************************
!
      USE mod_param
#ifdef SOLVE3D
      USE mod_coupling
#endif
      USE mod_grid
      USE mod_ocean
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model, Iter
!
!  Local variable declarations.
!
#include "tile.h"
!
#ifdef PROFILE
      CALL wclock_on (ng, model, 36)
#endif
      CALL cgradient_tile (ng, model, Istr, Iend, Jstr, Jend,           &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     Lold(ng), Lnew(ng), Iter,                    &
#ifdef MASKING
     &                     GRID(ng) % rmask,                            &
     &                     GRID(ng) % umask,                            &
     &                     GRID(ng) % vmask,                            &
#endif
#ifdef SOLVE3D
     &                     OCEAN(ng) % tl_t,                            &
     &                     OCEAN(ng) % tl_u,                            &
     &                     OCEAN(ng) % tl_v,                            &
#endif
     &                     OCEAN(ng) % tl_ubar,                         &
     &                     OCEAN(ng) % tl_vbar,                         &
     &                     OCEAN(ng) % tl_zeta,                         &
#ifdef SOLVE3D
     &                     OCEAN(ng) % d_t,                             &
     &                     OCEAN(ng) % d_u,                             &
     &                     OCEAN(ng) % d_v,                             &
#endif
     &                     OCEAN(ng) % d_ubar,                          &
     &                     OCEAN(ng) % d_vbar,                          &
     &                     OCEAN(ng) % d_zeta,                          &
#ifdef SOLVE3D
     &                     OCEAN(ng) % ad_t,                            &
     &                     OCEAN(ng) % ad_u,                            &
     &                     OCEAN(ng) % ad_v,                            &
#endif
     &                     OCEAN(ng) % ad_ubar,                         &
     &                     OCEAN(ng) % ad_vbar,                         &
     &                     OCEAN(ng) % ad_zeta)
#ifdef PROFILE
      CALL wclock_on (ng, model, 36)
#endif
      RETURN
      END SUBROUTINE cgradient
!
!***********************************************************************
      SUBROUTINE cgradient_tile (ng, model, Istr, Iend, Jstr, Jend,     &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           Lold, Lnew, Iter,                      &
#ifdef MASKING
     &                           rmask, umask, vmask,                   &
#endif
#ifdef SOLVE3D
     &                           tl_t, tl_u, tl_v,                      &
#endif
     &                           tl_ubar, tl_vbar, tl_zeta,             &
#ifdef SOLVE3D
     &                           d_t, d_u, d_v,                         &
#endif
     &                           d_ubar, d_vbar, d_zeta,                &
#ifdef SOLVE3D
     &                           ad_t, ad_u, ad_v,                      &
#endif
     &                           ad_ubar, ad_vbar, ad_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_iounits
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Lold, Lnew, Iter
!
#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
# endif
# ifdef SOLVE3D
      real(r8), intent(inout) :: ad_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_v(LBi:,LBj:,:,:)
# endif
      real(r8), intent(inout) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_vbar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_zeta(LBi:,LBj:,:)
# ifdef SOLVE3D
      real(r8), intent(inout) :: d_t(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: d_u(LBi:,LBj:,:)
      real(r8), intent(inout) :: d_v(LBi:,LBj:,:)
# endif
      real(r8), intent(inout) :: d_ubar(LBi:,LBj:)
      real(r8), intent(inout) :: d_vbar(LBi:,LBj:)
      real(r8), intent(inout) :: d_zeta(LBi:,LBj:)
# ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_v(LBi:,LBj:,:,:)
# endif
      real(r8), intent(inout) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_vbar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_zeta(LBi:,LBj:,:)
#else
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
# endif
# ifdef SOLVE3D
      real(r8), intent(inout) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
# endif
      real(r8), intent(inout) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_zeta(LBi:UBi,LBj:UBj,3)
# ifdef SOLVE3D
      real(r8), intent(inout) :: d_t(LBi:UBi,LBj:UBj,N(ng),NT(ng))
      real(r8), intent(inout) :: d_u(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: d_v(LBi:UBi,LBj:UBj,N(ng))
# endif
      real(r8), intent(inout) :: d_ubar(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: d_vbar(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: d_zeta(LBi:UBi,LBj:UBj)
# ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: tl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: tl_v(LBi:UBi,LBj:UBj,N(ng),2)
# endif
      real(r8), intent(inout) :: tl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: tl_vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: tl_zeta(LBi:UBi,LBj:UBj,3)
#endif
!
!  Local variable declarations.
!
      integer :: Linp, Lout, Lwrk, Lwrk1, i

      real(r8) :: norm, zgnorm

      real(r8), save :: alphaK, betaK, tauK, zbet

      real(r8), dimension(0:NstateVar(ng)) :: Adjust
      real(r8), dimension(0:NstateVar(ng)) :: dot_old, dot_new
      real(r8), dimension(0:NstateVar(ng)) :: old_dot, new_dot
      real(r8), dimension(1:Ninner,3) :: zwork
      real(r8), dimension(1:Ninner) :: zdelta
      real(r8), dimension(1:Ninner+1) :: zbeta
      real(r8), dimension(1:Ninner+1) :: zqg
      real(r8), dimension(1:Ninner) :: zu
      real(r8), dimension(1:Ninner) :: zgam
!
!-----------------------------------------------------------------------
!  Initialize trial step size.
!-----------------------------------------------------------------------
!
      tauK=CGstepI            
      alphaK=tauK
      IF (Iter.eq.0) THEN
        DO i=0,NstateVar(ng)
          dot_old(i)=0.0_r8
          dot_new(i)=0.0_r8
          old_dot(i)=0.0_r8
          new_dot(i)=0.0_r8
          FOURDVAR(ng)%CostGradDot(i)=0.0_r8
        END DO
      END IF
      WRITE (stdout,10)
 10   FORMAT (/,' <<<< Descent Algorithm >>>>')
!
!   Estimate the Hessian and save the starting vector in ad_*(Lold).
!
      IF (Iter.gt.0) THEN
       Lwrk=2
       Linp=1
       Lout=2
      CALL hessian (ng, model, Istr, Iend, Jstr, Jend,                  &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   Linp, Lout, Lwrk, Iter, tauK, zdelta,          &
#ifdef MASKING
     &                   rmask, umask, vmask,                           &
#endif
#ifdef SOLVE3D
     &                   ad_t, ad_u, ad_v,                              &
#endif
     &                   ad_ubar, ad_vbar, ad_zeta,                     &
#ifdef SOLVE3D
     &                   tl_t, tl_u, tl_v,                              &
#endif
     &                   tl_ubar, tl_vbar, tl_zeta)
!
     END IF
!
!    Apply the Lanczos recurrence and orthonormalize.
!
      Lwrk=2
      Linp=1
      Lout=2
      CALL lanczos (ng, model, Istr, Iend, Jstr, Jend,                  &
     &                      LBi, UBi, LBj, UBj,                         &
     &              Linp, Lout, Lwrk, Iter, zdelta, zbeta, zqg, zgnorm, &
#  ifdef MASKING
     &                      rmask, umask, vmask,                        &
#  endif
#  ifdef SOLVE3D
     &                      tl_t, tl_u, tl_v,                           &
#  endif
     &                      tl_ubar, tl_vbar, tl_zeta,                  &
#  ifdef SOLVE3D
     &                      ad_t, ad_u, ad_v,                           &
#  endif
     &                      ad_ubar, ad_vbar, ad_zeta)
!
      betaK=0.0_r8
!
!  Compute new direction, d(k+1).
!
      CALL new_direction (ng, model, Istr, Iend, Jstr, Jend,            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    Linp, Lout, betaK,                            &
#ifdef MASKING
     &                    rmask, umask, vmask,                          &
#endif
#ifdef SOLVE3D
     &                    ad_t, ad_u, ad_v,                             &
#endif
     &                    ad_ubar, ad_vbar, ad_zeta,                    &
#ifdef SOLVE3D
     &                    d_t, d_u, d_v,                                &
#endif
     &                    d_ubar, d_vbar, d_zeta)
!
!   Calculate the reduction in the gradient norm.
!
!   Set up the tridiagonal matrix coefficients.
!
      IF (Iter.gt.1) THEN
!
      DO i=1,Iter
        zwork(i,1)=zdelta(i)
        zwork(i,3)=-zqg(i)
      END DO
! AMM: What is the range for zbeta? zwork(2:iter,2)=zbeta(2:iter+1)
      DO i=2,Iter
        zwork(i,2)=zbeta(i)
      END DO
!
!  Solve tridiagonal system.
!
!      CALL DPTSV(Iter,1,zwork(1,1),zwork(2,2),zwork(1,3),Ninner,info)
!
!      IF (info.ne.0) THEN
!        print *,'Error in DPTSV: info=',info
!        stop
!      END IF
!
!  Solve tridiagonal system.
!
!   Decomposition and forward substitution
!
      zbet=zdelta(1)
      zu(1)=-zqg(1)/zbet
      DO i=2,Iter
        zgam(i)=zbeta(i)/zbet
        zbet=zdelta(i)-zbeta(i)*zgam(i)
        zu(i)=(-zqg(i)-zbeta(i)*zu(i-1))/zbet
      END DO
!
      zwork(Iter,3)=zu(Iter)
!
!   Back substitution.
!
      DO i=Iter-1,1,-1
        zu(i)=zu(i)-zgam(i+1)*zu(i+1)
        zwork(i,3)=zu(i)
      END DO
!
!  Compute gradient norm using tl*(:,:,1) and tl_*(:,:,2) as temporary storage.
!
      Lwrk=2
      Lwrk1=1
      Linp=1
      Lout=2
      CALL new_gradient (ng, model, Istr, Iend, Jstr, Jend,             &
     &                      LBi, UBi, LBj, UBj,                         &
     &              Linp, Lout, Lwrk, Lwrk1, Iter,                      &
     &              zgnorm, zbeta, zwork, zqg,                          &
#  ifdef MASKING
     &                      rmask, umask, vmask,                        &
#  endif
#  ifdef SOLVE3D
     &                      tl_t, tl_u, tl_v,                           &
#  endif
     &                      tl_ubar, tl_vbar, tl_zeta,                  &
#  ifdef SOLVE3D
     &                      ad_t, ad_u, ad_v,                           &
#  endif
     &                      ad_ubar, ad_vbar, ad_zeta)
      END IF
!
!
!-----------------------------------------------------------------------
!  Set TLM initial conditions for next inner loop, X(k+1).
!-----------------------------------------------------------------------
!
!    X(k+1) = tau(k+1) * d(k+1)
!
!    For the Lanczos algorithm, X(Linp) is ALWAYS the starting TL initial 
!    condition which for IS4DVAR is zero.
!
      Lout=2
      CALL tl_new_state_lanczos (ng, Istr, Iend, Jstr, Jend,            &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   Lout, tauK,                                    &
#ifdef MASKING
     &                   rmask, umask, vmask,                           &
#endif
#ifdef SOLVE3D
     &                   d_t, d_u, d_v,                                 &
#endif
     &                   d_ubar, d_vbar, d_zeta,                        &
#ifdef SOLVE3D
     &                   tl_t, tl_u, tl_v,                              &
#endif
     &                   tl_ubar, tl_vbar, tl_zeta)
!
!-----------------------------------------------------------------------
!  Report descent algorithm parameters.
!-----------------------------------------------------------------------
!
      IF (Master) THEN
        WRITE (stdout,20) outer,inner,tauK,alphaK,betaK,                &
     &                    outer,MAX(0,inner-1),Adjust(0),               &
     &                    outer,inner,                                  &
     &                    'dot product',inner,inner,dot_old(0),'alpha', &
     &                    'dot product',inner,inner,dot_new(0),'alpha', &
     &                    'dot product',inner,inner,old_dot(0),'beta',  &
     &                    'dot product',inner+1,inner+1,new_dot(0),     &
     &                    'beta'
 20     FORMAT (/,1x,'(',i3.3,',',i3.3,'): ',                           &
     &          'tau = ',1p,e14.7,                                      &
     &          ', alpha = ',1p,e14.7,                                  &
     &          ', Beta = ',1p,e14.7,                                   &
     &          /,1x,'(',i3.3,',',i3.3,'): ',                           &
     &          'Total COST Function Adjustment = ',1p,e19.12,          &
     &          /,1x,'(',i3.3,',',i3.3,'): ',                           &
     &          a,' <d(',i3.3,'),G(',i3.3,')> = ',1p,e19.12,3x,a,/,12x, &
     &          a,' <d(',i3.3,'),g(',i3.3,')> = ',1p,e19.12,3x,a,/,12x, &
     &          a,' <G(',i3.3,'),G(',i3.3,')> = ',1p,e19.12,3x,a,/,12x, &
     &          a,' <G(',i3.3,'),G(',i3.3,')> = ',1p,e19.12,3x,a,/)
      END IF

      RETURN 
      END SUBROUTINE cgradient_tile

!
!***********************************************************************
      SUBROUTINE tl_new_state_lanczos (ng, Istr, Iend, Jstr, Jend,      &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         Lout, alphaK,                            &
#ifdef MASKING
     &                         rmask, umask, vmask,                     &
#endif
#ifdef SOLVE3D
     &                         d_t, d_u, d_v,                           &
#endif
     &                         d_ubar, d_vbar, d_zeta,                  &
#ifdef SOLVE3D
     &                         tl_t, tl_u, tl_v,                        &
#endif
     &                         tl_ubar, tl_vbar, tl_zeta)
!***********************************************************************
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Lout

      real(r8), intent(in) :: alphaK
!
#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
# endif
# ifdef SOLVE3D
      real(r8), intent(inout) :: d_t(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: d_u(LBi:,LBj:,:)
      real(r8), intent(inout) :: d_v(LBi:,LBj:,:)
# endif
      real(r8), intent(inout) :: d_ubar(LBi:,LBj:)
      real(r8), intent(inout) :: d_vbar(LBi:,LBj:)
      real(r8), intent(inout) :: d_zeta(LBi:,LBj:)
# ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_v(LBi:,LBj:,:,:)
# endif
      real(r8), intent(inout) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_vbar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_zeta(LBi:,LBj:,:)
#else
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
# endif
# ifdef SOLVE3D
      real(r8), intent(inout) :: d_t(LBi:UBi,LBj:UBj,N(ng),NT(ng))
      real(r8), intent(inout) :: d_u(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: d_v(LBi:UBi,LBj:UBj,N(ng))
# endif
      real(r8), intent(inout) :: d_ubar(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: d_vbar(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: d_zeta(LBi:UBi,LBj:UBj)
# ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: tl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: tl_v(LBi:UBi,LBj:UBj,N(ng),2)
# endif
      real(r8), intent(inout) :: tl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: tl_vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: tl_zeta(LBi:UBi,LBj:UBj,3)
#endif
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
!  Compute new starting tangent linear state vector, X(k+1).
!-----------------------------------------------------------------------
!
!  2D state variables.
!
#ifndef SOLVE3D
      DO j=JstrR,JendR
        DO i=Istr,IendR
          tl_ubar(i,j,Lout)=alphaK*d_ubar(i,j)
# ifdef MASKING
          tl_ubar(i,j,Lout)=tl_ubar(i,j,Lout)*umask(i,j)
# endif
        END DO
      END DO
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          tl_vbar(i,j,Lout)=alphaK*d_vbar(i,j)
# ifdef MASKING
          tl_vbar(i,j,Lout)=tl_vbar(i,j,Lout)*vmask(i,j)
# endif
        END DO
      END DO
#endif
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          tl_zeta(i,j,Lout)=alphaK*d_zeta(i,j)
#ifdef MASKING
     &                      *rmask(i,j)
#endif
        END DO
      END DO
#ifdef SOLVE3D
!
!  3D state variables.
!
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=Istr,IendR
            tl_u(i,j,k,Lout)=alphaK*d_u(i,j,k)
# ifdef MASKING
     &                       *umask(i,j)
# endif
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            tl_v(i,j,k,Lout)=alphaK*d_v(i,j,k)
# ifdef MASKING
     &                       *vmask(i,j)
# endif
          END DO
        END DO
      END DO
!
      DO itrc=1,NT(ng)
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              tl_t(i,j,k,Lout,itrc)=alphaK*d_t(i,j,k,itrc)
# ifdef MASKING
     &                              *rmask(i,j)
# endif
            END DO
          END DO          
        END DO
      END DO
#endif

      RETURN
      END SUBROUTINE tl_new_state_lanczos
!
!***********************************************************************
      SUBROUTINE ad_new_state (ng, Istr, Iend, Jstr, Jend,              &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         Lold, Lnew, alphaK, tauK,                &
#ifdef MASKING
     &                         rmask, umask, vmask,                     &
#endif
#ifdef SOLVE3D
     &                         ad_t, ad_u, ad_v,                        &
#endif
     &                         ad_ubar, ad_vbar, ad_zeta)
!***********************************************************************
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Lold, Lnew

      real(r8), intent(in) :: alphaK, tauK
!
#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
# endif
# ifdef SOLVE3D
      real(r8), intent(inout) :: ad_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_v(LBi:,LBj:,:,:)
# endif
      real(r8), intent(inout) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_vbar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_zeta(LBi:,LBj:,:)
#else
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
# endif
# ifdef SOLVE3D
      real(r8), intent(inout) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
# endif
      real(r8), intent(inout) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_zeta(LBi:UBi,LBj:UBj,3)
#endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j
#ifdef SOLVE3D
      integer :: itrc, k
#endif
      real(r8) :: fac

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Estimate the gradient for the new state vector, G(k+1). Notice that
!  the Lnew record is overwritten.
!-----------------------------------------------------------------------
!
      fac=alphaK/tauK
!
!  2D state variables.
!
#ifndef SOLVE3D
      DO j=JstrR,JendR
        DO i=Istr,IendR
          ad_ubar(i,j,Lnew)=ad_ubar(i,j,Lold)+                          &
# ifdef MASKING
     &                      umask(i,j)*                                 &
# endif
     &                      fac*(ad_ubar(i,j,Lnew)-                     &
     &                           ad_ubar(i,j,Lold))
          ad_ubar(i,j,Lold)=ad_ubar(i,j,Lnew)
        END DO
      END DO
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          ad_vbar(i,j,Lnew)=ad_vbar(i,j,Lold)+                          &
# ifdef MASKING
     &                      vmask(i,j)*                                 &
# endif
     &                      fac*(ad_vbar(i,j,Lnew)-                     &
     &                           ad_vbar(i,j,Lold))
          ad_vbar(i,j,Lold)=ad_vbar(i,j,Lnew)
        END DO
      END DO
#endif
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ad_zeta(i,j,Lnew)=ad_zeta(i,j,Lold)+                          &
#ifdef MASKING
     &                      rmask(i,j)*                                 &
#endif
     &                      fac*(ad_zeta(i,j,Lnew)-                     &
     &                           ad_zeta(i,j,Lold))
          ad_zeta(i,j,Lold)=ad_zeta(i,j,Lnew)
        END DO
      END DO
#ifdef SOLVE3D
!
!  3D state variables.
!
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=Istr,IendR
            ad_u(i,j,k,Lnew)=ad_u(i,j,k,Lold)+                          &
# ifdef MASKING
     &                       umask(i,j)*                                &
# endif
     &                       fac*(ad_u(i,j,k,Lnew)-                     &
     &                            ad_u(i,j,k,Lold))
            ad_u(i,j,k,Lold)=ad_u(i,j,k,Lnew)
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            ad_v(i,j,k,Lnew)=ad_v(i,j,k,Lold)+                          &
# ifdef MASKING
     &                       vmask(i,j)*                                &
# endif
     &                       fac*(ad_v(i,j,k,Lnew)-                     &
     &                            ad_v(i,j,k,Lold))
            ad_v(i,j,k,Lold)=ad_v(i,j,k,Lnew)
          END DO
        END DO
      END DO
!
      DO itrc=1,NT(ng)
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              ad_t(i,j,k,Lnew,itrc)=ad_t(i,j,k,Lold,itrc)+              &
# ifdef MASKING
     &                              rmask(i,j)*                         &
# endif
     &                              fac*(ad_t(i,j,k,Lnew,itrc)-         &
     &                                   ad_t(i,j,k,Lold,itrc))
              ad_t(i,j,k,Lold,itrc)=ad_t(i,j,k,Lnew,itrc)
            END DO
          END DO
        END DO
      END DO
#endif

      RETURN
      END SUBROUTINE ad_new_state
!
!***********************************************************************
      SUBROUTINE orthogonalize (ng, model, Istr, Iend, Jstr, Jend,      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Lold, Lnew, Lwrk, Iter,                 &
#ifdef MASKING
     &                          rmask, umask, vmask,                    &
#endif
#ifdef SOLVE3D
     &                          tl_t, tl_u, tl_v,                       &
#endif
     &                          tl_ubar, tl_vbar, tl_zeta,              &
#ifdef SOLVE3D
     &                          ad_t, ad_u, ad_v,                       &
#endif
     &                          ad_ubar, ad_vbar, ad_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_iounits
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Lold, Lnew, Lwrk, Iter
!
#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
# endif
# ifdef SOLVE3D
      real(r8), intent(inout) :: ad_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_v(LBi:,LBj:,:,:)
# endif
      real(r8), intent(inout) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_vbar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_zeta(LBi:,LBj:,:)
# ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_v(LBi:,LBj:,:,:)
# endif
      real(r8), intent(inout) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_vbar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_zeta(LBi:,LBj:,:)
#else
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
# endif
# ifdef SOLVE3D
      real(r8), intent(inout) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
# endif
      real(r8), intent(inout) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_zeta(LBi:UBi,LBj:UBj,3)
# ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: tl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: tl_v(LBi:UBi,LBj:UBj,N(ng),2)
# endif
      real(r8), intent(inout) :: tl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: tl_vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: tl_zeta(LBi:UBi,LBj:UBj,3)
#endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j, lstr, rec
#ifdef SOLVE3D
      integer :: itrc, k
#endif
      real(r8) :: fac

      real(r8), dimension(0:NstateVar(ng)) :: dot
      real(r8), dimension(Iter) :: DotProd, dot_new, dot_old

      character (len=80) :: ncname

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Orthogonalize current gradient, G(k+1), against all previous
!  gradients (reverse order) using Gramm-Schmidt procedure.
!-----------------------------------------------------------------------
!
!  We can overwrite adjoint arrays at index Lnew each time around the
!  the following loop because the preceding gradient vectors that we
!  read are orthogonal to each other. The reversed order of the loop
!  is important for the Lanczos vector calculations.
!
      DO rec=Iter,1,-1
!
!  Determine adjoint file to process.
!
        IF (ndefADJ(ng).gt.0) THEN
          lstr=LEN_TRIM(ADJbase(ng))
          WRITE (ncname,10) ADJbase(ng)(1:lstr-3), rec
 10       FORMAT (a,'_',i3.3,'.nc')
        ELSE
          ncname=ADJname(ng)
        END IF
!
!  Read in each previous gradient state solutions, G(0) to G(k), and
!  compute its associated dot angaint curret G(k+1). Each gradient
!  solution is loaded into TANGENT LINEAR STATE ARRAYS at index Lwrk.
!
        CALL get_gradient (ng, model, Istr, Iend, Jstr, Jend,           &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     Lwrk, rec, ncname,                           &
#ifdef MASKING
     &                     rmask, umask, vmask,                         &
#endif
#ifdef SOLVE3D
     &                     tl_t, tl_u, tl_v,                            &
#endif
     &                     tl_ubar, tl_vbar, tl_zeta)
!
!  Compute dot product <G(k+1), G(rec)>.
!
        CALL state_dotprod (ng, model, Istr, Iend, Jstr, Jend,          &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NstateVar(ng), dot(0:),                     &
#ifdef MASKING
     &                      rmask, umask, vmask,                        &
#endif
#ifdef SOLVE3D
     &                      ad_t(:,:,:,Lnew,:), tl_t(:,:,:,Lwrk,:),     &
     &                      ad_u(:,:,:,Lnew), tl_u(:,:,:,Lwrk),         &
     &                      ad_v(:,:,:,Lnew), tl_v(:,:,:,Lwrk),         &
#endif
     &                      ad_ubar(:,:,Lnew), tl_ubar(:,:,Lwrk),       &
     &                      ad_vbar(:,:,Lnew), tl_vbar(:,:,Lwrk),       &
     &                      ad_zeta(:,:,Lnew), tl_zeta(:,:,Lwrk))
        dot_new(rec)=dot(0)
!
!  Compute dot product <G(rec), G(rec)>.
!
        CALL state_dotprod (ng, model, Istr, Iend, Jstr, Jend,          &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NstateVar(ng), dot(0:),                     &
#ifdef MASKING
     &                      rmask, umask, vmask,                        &
#endif
#ifdef SOLVE3D
     &                      tl_t(:,:,:,Lwrk,:), tl_t(:,:,:,Lwrk,:),     &
     &                      tl_u(:,:,:,Lwrk), tl_u(:,:,:,Lwrk),         &
     &                      tl_v(:,:,:,Lwrk), tl_v(:,:,:,Lwrk),         &
#endif
     &                      tl_ubar(:,:,Lwrk), tl_ubar(:,:,Lwrk),       &
     &                      tl_vbar(:,:,Lwrk), tl_vbar(:,:,Lwrk),       &
     &                      tl_zeta(:,:,Lwrk), tl_zeta(:,:,Lwrk))
        dot_old(rec)=dot(0)
!
!  Compute Gramm-Schmidt scaling coefficient.
!
        DotProd(rec)=dot_new(rec)/dot_old(rec)
!
!  Gramm-Schmidt orthonormalization, 2D state gradient.
!
#ifndef SOLVE3D
        DO j=JstrR,JendR
          DO i=Istr,IendR
            ad_ubar(i,j,Lnew)=ad_ubar(i,j,Lnew)-                        &
# ifdef MASKING
     &                        umask(i,j)*                               &
# endif
     &                        DotProd(rec)*tl_ubar(i,j,Lwrk)
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            ad_vbar(i,j,Lnew)=ad_vbar(i,j,Lnew)-                        &
# ifdef MASKING
     &                        vmask(i,j)*                               &
# endif
     &                        DotProd(rec)*tl_vbar(i,j,Lwrk)
          END DO
        END DO
#endif
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            ad_zeta(i,j,Lnew)=ad_zeta(i,j,Lnew)-                        &
# ifdef MASKING
     &                        rmask(i,j)*                               &
# endif
     &                        DotProd(rec)*tl_zeta(i,j,Lwrk)
          END DO
        END DO
#ifdef SOLVE3D
!
!  Gramm-Schmidt orthonormalization, 3D state gradient.
!
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=Istr,IendR
              ad_u(i,j,k,Lnew)=ad_u(i,j,k,Lnew)-                        &
# ifdef MASKING
     &                         umask(i,j)*                              &
# endif
     &                         DotProd(rec)*tl_u(i,j,k,Lwrk)
            END DO
          END DO
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              ad_v(i,j,k,Lnew)=ad_v(i,j,k,Lnew)-                        &
# ifdef MASKING
     &                         vmask(i,j)*                              &
# endif
     &                         DotProd(rec)*tl_v(i,j,k,Lwrk)
            END DO
          END DO
        END DO
!
        DO itrc=1,NT(ng)
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                ad_t(i,j,k,Lnew,itrc)=ad_t(i,j,k,Lnew,itrc)-            &
# ifdef MASKING
     &                                rmask(i,j)*                       &
# endif
     &                                DotProd(rec)*tl_t(i,j,k,Lwrk,itrc)
              END DO
            END DO
          END DO
        END DO
#endif
      END DO
#ifdef NORMALIZATION
!
!-----------------------------------------------------------------------
!  Normalize current orthogonal gradient vector.
!-----------------------------------------------------------------------
!
      CALL state_dotprod (ng, model, Istr, Iend, Jstr, Jend,            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NstateVar(ng), dot(0:),                       &
# ifdef MASKING
     &                    rmask, umask, vmask,                          &
# endif
# ifdef SOLVE3D
     &                    ad_t(:,:,:,Lnew,:), ad_t(:,:,:,Lnew,:),       &
     &                    ad_u(:,:,:,Lnew), ad_u(:,:,:,Lnew),           &
     &                    ad_v(:,:,:,Lnew), ad_v(:,:,:,Lnew),           &
# endif
     &                    ad_ubar(:,:,Lnew), ad_ubar(:,:,Lnew),         &
     &                    ad_vbar(:,:,Lnew), ad_vbar(:,:,Lnew),         &
     &                    ad_zeta(:,:,Lnew), ad_zeta(:,:,Lnew))
!
!  Compute normaliztion factor.
!
      fac=1.0_r8/SQRT(dot(0))
!
!  Normalize current 2D state gradient vector.
!
# ifndef SOLVE3D
      DO j=JstrR,JendR
        DO i=Istr,IendR
          ad_ubar(i,j,Lnew)=fac*ad_ubar(i,j,Lnew)
        END DO
      END DO
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          ad_vbar(i,j,Lnew)=fac*ad_vbar(i,j,Lnew)
        END DO
      END DO
# endif
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ad_zeta(i,j,Lnew)=fac*ad_zeta(i,j,Lnew)
        END DO
      END DO
# ifdef SOLVE3D
!
!  Normalize current 3D state gradient vector.
!
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=Istr,IendR
            ad_u(i,j,k,Lnew)=fac*ad_u(i,j,k,Lnew)
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            ad_v(i,j,k,Lnew)=fac*ad_v(i,j,k,Lnew)
          END DO
        END DO
      END DO
!
      DO itrc=1,NT(ng)
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              ad_t(i,j,k,Lnew,itrc)=fac*ad_t(i,j,k,Lnew,itrc)
            END DO
          END DO
        END DO
      END DO
# endif
#endif
#ifdef TEST_ORTHOGONALIZATION
!
!-----------------------------------------------------------------------
!  Test orthogonal properties of the new gradient.
!-----------------------------------------------------------------------
!
      DO rec=Iter,1,-1
!
!  Determine adjoint file to process.
!
        IF (ndefADJ(ng).gt.0) THEN
          lstr=LEN_TRIM(ADJbase(ng))
          WRITE (ncname,10) ADJbase(ng)(1:lstr-3), rec
        ELSE
          ncname=ADJname(ng)
        END IF
!
!  Read in each previous gradient state solutions, G(0) to G(k), and
!  compute its associated dot angaint orthogonalized G(k+1). Again, 
!  each gradient solution is loaded into TANGENT LINEAR STATE ARRAYS
!  at index Lwrk.
!
        CALL get_gradient (ng, model, Istr, Iend, Jstr, Jend,           &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     Lwrk, rec, ncname,                           &
#ifdef MASKING
     &                     rmask, umask, vmask,                         &
#endif
#ifdef SOLVE3D
     &                     tl_t, tl_u, tl_v,                            &
#endif
     &                     tl_ubar, tl_vbar, tl_zeta)
!
        CALL state_dotprod (ng, model, Istr, Iend, Jstr, Jend,          &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NstateVar(ng), dot(0:),                     &
#ifdef MASKING
     &                      rmask, umask, vmask,                        &
#endif
#ifdef SOLVE3D
     &                      ad_t(:,:,:,Lnew,:), tl_t(:,:,:,Lwrk,:),     &
     &                      ad_u(:,:,:,Lnew), tl_u(:,:,:,Lwrk),         &
     &                      ad_v(:,:,:,Lnew), tl_v(:,:,:,Lwrk),         &
#endif
     &                      ad_ubar(:,:,Lnew), tl_ubar(:,:,Lwrk),       &
     &                      ad_vbar(:,:,Lnew), tl_vbar(:,:,Lwrk),       &
     &                      ad_zeta(:,:,Lnew), tl_zeta(:,:,Lwrk))
        dot_new(rec)=dot(0)
      END DO
!
!  Report dot products. If everything is working correctly, at the
!  end of the orthogonalization dot_new(rec) << dot_old(rec).
!
      IF (Master) THEN
        WRITE (stdout,20) outer, inner
        DO rec=Iter,1,-1
          WRITE (stdout,30) DotProd(rec), rec-1
        END DO
        WRITE (stdout,*) ' '
        DO rec=Iter,1,-1
          WRITE (stdout,40) Iter, rec-1, dot_new(rec),                  &
     &                      rec-1, rec-1, dot_old(rec)
        END DO
 20     FORMAT (/,1x,'(',i3.3,',',i3.3,'): ',                           &
     &          'Gramm-Schmidt Orthogonalization:',/)
 30     FORMAT (12x,'Orthogonalization Factor = ',1p,e19.12,3x,         &
     &          '(Iter=',i3.3,')')
 40     FORMAT (2x,'Ortho Test: ',                                      &
     &          '<G(',i3.3,'),G(',i3.3,')> = ',1p,e15.8,1x,             &
     &          '<G(',i3.3,'),G(',i3.3,')> = ',1p,e15.8)
      END IF
#endif

      RETURN
      END SUBROUTINE orthogonalize
!
!***********************************************************************
      SUBROUTINE get_gradient (ng, model, Istr, Iend, Jstr, Jend,       &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         Lwrk, rec, ncname,                       &
#ifdef MASKING
     &                         rmask, umask, vmask,                     &
#endif
#ifdef SOLVE3D
     &                         tl_t, tl_u, tl_v,                        &
#endif
     &                         tl_ubar, tl_vbar, tl_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Lwrk, rec

      character (len=*), intent(in) :: ncname
!
#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
# endif
# ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_v(LBi:,LBj:,:,:)
# endif
      real(r8), intent(inout) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_vbar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_zeta(LBi:,LBj:,:)
#else
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
# endif
# ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t(LBi:UBi,LBj:UBj,N(ng),NT(ng))
      real(r8), intent(inout) :: tl_u(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: tl_v(LBi:UBi,LBj:UBj,N(ng))
# endif
      real(r8), intent(inout) :: tl_ubar(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: tl_vbar(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: tl_zeta(LBi:UBi,LBj:UBj)
#endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j
#ifdef SOLVE3D
      integer :: itrc, k
#endif
      integer :: ncid, status
      integer, dimension(NV) :: vid
      integer, dimension(4) :: Vsize

      integer :: nf_fread2d
#ifdef SOLVE3D
      integer :: nf_fread3d
#endif

      real(r8) :: Fmin, Fmax, scale

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Read in requested gradient record. Load gradient solution into
!  tangent linear state arrays at index Lwrk.
!-----------------------------------------------------------------------
!
!  Determine file and variables ids.
!
      IF (ndefADJ(ng).gt.0) THEN
        IF (InpThread) THEN
          status=nf_open(TRIM(ncname), nf_nowrite, ncid)
          IF (status.ne.nf_noerr) THEN
            WRITE (stdout,10) TRIM(ncname)
            exit_flag=2
            ioerror=status
            RETURN
          END IF            
#ifndef SOLVE3D
          status=nf_inq_varid(ncid, TRIM(Vname(1,idUbar)), vid(idUbar))
          status=nf_inq_varid(ncid, TRIM(Vname(1,idVbar)), vid(idVbar))
#endif
          status=nf_inq_varid(ncid, TRIM(Vname(1,idFsur)), vid(idFsur))
#ifdef SOLVE3D
          status=nf_inq_varid(ncid, TRIM(Vname(1,idUvel)), vid(idUvel))
          status=nf_inq_varid(ncid, TRIM(Vname(1,idVvel)), vid(idVvel))
          DO itrc=1,NT(ng)
            status=nf_inq_varid(ncid, TRIM(Vname(1,idTvar(itrc))),      &
     &                          vid(idTvar(itrc)))
          END DO
#endif
        END IF
      ELSE
        ncid=ncADJid(ng)
#ifndef SOLVE3D
        vid(idUbar)=adjVid(idUbar,ng)
        vid(idVbar)=adjVid(idVbar,ng)
#endif
        vid(idFsur)=adjVid(idFsur,ng)
#ifdef SOLVE3D
        vid(idUvel)=adjVid(idUvel,ng)
        vid(idVvel)=adjVid(idVvel,ng)
        DO itrc=1,NT(ng)
          vid(idTvar(itrc))=adjTid(itrc,ng)
        END DO
#endif
      END IF
      DO i=1,4
        Vsize(i)=0
      END DO
      scale=1.0_r8
#ifndef SOLVE3D
!
!  Read in 2D adjoint momentum.
!
      status=nf_fread2d(ng, iTLM, ncid, vid(idUbar), rec, u2dvar,       &
     &                  Vsize, LBi, UBi, LBj, UBj,                      &
     &                  scale, Fmin, Fmax,                              &
# ifdef MASKING
     &                  umask(LBi,LBj),                                 &
# endif
     &                  tl_ubar(LBi,LBj,Lwrk))
      IF (status.ne.nf_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idUbar)), rec, TRIM(ncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
      status=nf_fread2d(ng, iTLM, ncid, vid(idVbar), rec, v2dvar,       &
     &                  Vsize, LBi, UBi, LBj, UBj,                      &
     &                  scale, Fmin, Fmax,                              &
# ifdef MASKING
     &                  vmask(LBi,LBj),                                 &
# endif
     &                  tl_vbar(LBi,LBj,Lwrk))
      IF (status.ne.nf_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idVbar)), rec, TRIM(ncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
#endif
!
!  Read in adjoint free-surface
!
      status=nf_fread2d(ng, iTLM, ncid, vid(idFsur), rec, r2dvar,       &
     &                  Vsize, LBi, UBi, LBj, UBj,                      &
     &                  scale, Fmin, Fmax,                              &
#ifdef MASKING
     &                  rmask(LBi,LBj),                                 &
#endif
     &                  tl_zeta(LBi,LBj,Lwrk))
      IF (status.ne.nf_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idFsur)), rec, TRIM(ncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
#ifdef SOLVE3D
!
!  Read in adjoint 3D momentum.
!
      status=nf_fread3d(ng, iTLM, ncid, vid(idUvel), rec, u3dvar,       &
     &                  Vsize, LBi, UBi, LBj, UBj, 1, N(ng),            &
     &                  scale, Fmin, Fmax,                              &
# ifdef MASKING
     &                  umask(LBi,LBj),                                 &
# endif
     &                  tl_u(LBi,LBj,1,Lwrk))
      IF (status.ne.nf_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idUvel)), rec, TRIM(ncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
      status=nf_fread3d(ng, iTLM, ncid, vid(idVvel), rec, v3dvar,       &
     &                  Vsize, LBi, UBi, LBj, UBj, 1, N(ng),            &
     &                  scale, Fmin, Fmax,                              &
# ifdef MASKING
     &                  vmask(LBi,LBj),                                 &
# endif
     &                  tl_v(LBi,LBj,1,Lwrk))
      IF (status.ne.nf_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idVvel)), rec, TRIM(ncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Read in adjoint tracers.
!
      DO itrc=1,NT(ng)
        status=nf_fread3d(ng, iTLM, ncid, vid(idTvar(itrc)), rec,       &
     &                    r3dvar, Vsize, LBi, UBi, LBj, UBj, 1, N(ng),  &
     &                    scale, Fmin, Fmax,                            &
# ifdef MASKING
     &                    rmask(LBi,LBj),                               &
# endif
     &                    tl_t(LBi,LBj,1,Lwrk,itrc))
        IF (status.ne.nf_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idTvar(itrc))), rec,         &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END DO
#endif
!
!  If multiple files, close adjoint history file.
!
      IF (ndefADJ(ng).gt.0) THEN
        status=nf_close(ncid)
      END IF
!
 10   FORMAT (' GET_GRADIENT - unable to open NetCDF file: ',a)
 20   FORMAT (' GET_GRADIENT - error while reading variable: ',a,2x,    &
     &        'at time record = ',i3,/,16x,'in NetCDF file: ',a)

      RETURN
      END SUBROUTINE get_gradient
!
!***********************************************************************
      SUBROUTINE new_direction (ng, model, Istr, Iend, Jstr, Jend,      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Lold, Lnew, betaK,                      &
#ifdef MASKING
     &                          rmask, umask, vmask,                    &
#endif
#ifdef SOLVE3D
     &                          ad_t, ad_u, ad_v,                       &
#endif
     &                          ad_ubar, ad_vbar, ad_zeta,              &
#ifdef SOLVE3D
     &                          d_t, d_u, d_v,                          &
#endif
     &                          d_ubar, d_vbar, d_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Lold, Lnew

      real(r8), intent(in) :: betaK      
!
#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
# endif
# ifdef SOLVE3D
      real(r8), intent(inout) :: ad_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_v(LBi:,LBj:,:,:)
# endif
      real(r8), intent(inout) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_vbar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_zeta(LBi:,LBj:,:)
# ifdef SOLVE3D
      real(r8), intent(inout) :: d_t(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: d_u(LBi:,LBj:,:)
      real(r8), intent(inout) :: d_v(LBi:,LBj:,:)
# endif
      real(r8), intent(inout) :: d_ubar(LBi:,LBj:)
      real(r8), intent(inout) :: d_vbar(LBi:,LBj:)
      real(r8), intent(inout) :: d_zeta(LBi:,LBj:)
#else
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
# endif
# ifdef SOLVE3D
      real(r8), intent(inout) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
# endif
      real(r8), intent(inout) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_zeta(LBi:UBi,LBj:UBj,3)
# ifdef SOLVE3D
      real(r8), intent(inout) :: d_t(LBi:UBi,LBj:UBj,N(ng),NT(ng))
      real(r8), intent(inout) :: d_u(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: d_v(LBi:UBi,LBj:UBj,N(ng))
# endif
      real(r8), intent(inout) :: d_ubar(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: d_vbar(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: d_zeta(LBi:UBi,LBj:UBj)
#endif
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
!  Compute new conjugate descent direction, d(k+1). Notice that the old
!  descent direction is overwritten. Also the initial value is just
!  d(0)=-G(0) since betaK=0 when Iter=0.
!-----------------------------------------------------------------------
!
!  2D state variables.
!
#ifndef SOLVE3D
      DO j=JstrR,JendR
        DO i=Istr,IendR
          d_ubar(i,j)=-ad_ubar(i,j,Lnew)+betaK*d_ubar(i,j)
# ifdef MASKING
          d_ubar(i,j)=d_ubar(i,j)*umask(i,j)
# endif
        END DO
      END DO
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          d_vbar(i,j)=-ad_vbar(i,j,Lnew)+betaK*d_vbar(i,j)
# ifdef MASKING
          d_vbar(i,j)=d_vbar(i,j)*vmask(i,j)
# endif
        END DO
      END DO
#endif
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          d_zeta(i,j)=-ad_zeta(i,j,Lnew)+betaK*d_zeta(i,j)
# ifdef MASKING
          d_zeta(i,j)=d_zeta(i,j)*rmask(i,j)
# endif
        END DO
      END DO
#ifdef SOLVE3D
!
!  3D state variables.
!
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=Istr,IendR
            d_u(i,j,k)=-ad_u(i,j,k,Lnew)+betaK*d_u(i,j,k)
# ifdef MASKING
            d_u(i,j,k)=d_u(i,j,k)*umask(i,j)
# endif
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            d_v(i,j,k)=-ad_v(i,j,k,Lnew)+betaK*d_v(i,j,k)
# ifdef MASKING
            d_v(i,j,k)=d_v(i,j,k)*vmask(i,j)
# endif
          END DO
        END DO
      END DO
!
      DO itrc=1,NT(ng)
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              d_t(i,j,k,itrc)=-ad_t(i,j,k,Lnew,itrc)+                   &
     &                        betaK*d_t(i,j,k,itrc)
# ifdef MASKING
              d_t(i,j,k,itrc)=d_t(i,j,k,itrc)*rmask(i,j)
# endif
            END DO
          END DO
        END DO
      END DO
#endif

      RETURN
      END SUBROUTINE new_direction
!
!***********************************************************************
      SUBROUTINE hessian (ng, model, Istr, Iend, Jstr, Jend,            &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         Lold, Lnew, Lwrk, Iter, tauK, zdelta,    &
#ifdef MASKING
     &                         rmask, umask, vmask,                     &
#endif
#ifdef SOLVE3D
     &                         ad_t, ad_u, ad_v,                        &
#endif
     &                         ad_ubar, ad_vbar, ad_zeta,               &
#ifdef SOLVE3D
     &                         tl_t, tl_u, tl_v,                        &
#endif
     &                         tl_ubar, tl_vbar, tl_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_fourdvar
      USE mod_iounits
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Lold, Lnew, Iter

      real(r8), intent(in) :: tauK
      real(r8), intent(inout) :: zdelta(1:Ninner)
!
#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
# endif
# ifdef SOLVE3D
      real(r8), intent(inout) :: ad_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_v(LBi:,LBj:,:,:)
# endif
      real(r8), intent(inout) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_vbar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_zeta(LBi:,LBj:,:)
# ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_v(LBi:,LBj:,:,:)
# endif
      real(r8), intent(inout) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_vbar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_zeta(LBi:,LBj:,:)
#else
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
# endif
# ifdef SOLVE3D
      real(r8), intent(inout) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
# endif
      real(r8), intent(inout) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_zeta(LBi:UBi,LBj:UBj,3)
# ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: tl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: tl_v(LBi:UBi,LBj:UBj,N(ng),2)
# endif
      real(r8), intent(inout) :: tl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: tl_vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: tl_zeta(LBi:UBi,LBj:UBj,3)
#endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j
      integer :: Lwrk, lstr, rec
#ifdef SOLVE3D
      integer :: itrc, k
#endif
      real(r8), dimension(0:NstateVar(ng)) :: dot

      character (len=80) :: ncname

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Estimate the Hessian.
!-----------------------------------------------------------------------
!
!  2D state variables.
!
      dot=0.0_r8
!
#ifndef SOLVE3D
      DO j=JstrR,JendR
        DO i=Istr,IendR
          ad_ubar(i,j,Lnew)=(ad_ubar(i,j,Lnew)-ad_ubar(i,j,Lold))       &
# ifdef MASKING
     &                      *umask(i,j)                                 &
# endif
     &                                 /tauK
        END DO
      END DO
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          ad_vbar(i,j,Lnew)=(ad_vbar(i,j,Lnew)-ad_vbar(i,j,Lold))       &
# ifdef MASKING
     &                      *vmask(i,j)                                 &
# endif
     &                                 /tauK
        END DO
      END DO
#endif
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ad_zeta(i,j,Lnew)=(ad_zeta(i,j,Lnew)-ad_zeta(i,j,Lold))       &
#ifdef MASKING
     &                      *rmask(i,j)                                 &
#endif
     &                                 /tauK
        END DO
      END DO
#ifdef SOLVE3D
!
!  3D state variables.
!
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=Istr,IendR
            ad_u(i,j,k,Lnew)=(ad_u(i,j,k,Lnew)-ad_u(i,j,k,Lold))        &
# ifdef MASKING
     &                       *umask(i,j)                                &
# endif
     &                                 /tauK
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            ad_v(i,j,k,Lnew)=(ad_v(i,j,k,Lnew)-ad_v(i,j,k,Lold))        &
# ifdef MASKING
     &                       *vmask(i,j)                                &
# endif
     &                                 /tauK
          END DO
        END DO
      END DO
!
      DO itrc=1,NT(ng)
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              ad_t(i,j,k,Lnew,itrc)=(ad_t(i,j,k,Lnew,itrc)-             &
     &                                     ad_t(i,j,k,Lold,itrc))       &
# ifdef MASKING
     &                              *rmask(i,j)                         &
# endif
     &                                 /tauK
            END DO
          END DO
        END DO
      END DO
#endif
!
!  Determine adjoint file to process.
!
        IF (ndefADJ(ng).gt.0) THEN
          lstr=LEN_TRIM(ADJbase(ng))
          WRITE (ncname,10) ADJbase(ng)(1:lstr-3), rec
 10       FORMAT (a,'_',i3.3,'.nc')
        ELSE
          ncname=ADJname(ng)
        END IF
!
        CALL get_gradient (ng, model, Istr, Iend, Jstr, Jend,           &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     Lwrk, Iter, ncname,                          &
#ifdef MASKING
     &                     rmask, umask, vmask,                         &
#endif
#ifdef SOLVE3D
     &                     tl_t, tl_u, tl_v,                            &
#endif
     &                     tl_ubar, tl_vbar, tl_zeta)
!
        CALL state_dotprod (ng, model, Istr, Iend, Jstr, Jend,          &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NstateVar(ng), dot(0:),                     &
#ifdef MASKING
     &                      rmask, umask, vmask,                        &
#endif
#ifdef SOLVE3D
     &                      ad_t(:,:,:,Lnew,:), tl_t(:,:,:,Lwrk,:),     &
     &                      ad_u(:,:,:,Lnew), tl_u(:,:,:,Lwrk),         &
     &                      ad_v(:,:,:,Lnew), tl_v(:,:,:,Lwrk),         &
#endif
     &                      ad_ubar(:,:,Lnew), tl_ubar(:,:,Lwrk),       &
     &                      ad_vbar(:,:,Lnew), tl_vbar(:,:,Lwrk),       &
     &                      ad_zeta(:,:,Lnew), tl_zeta(:,:,Lwrk))
        zdelta(Iter)=dot(0)
!
      IF (zdelta(Iter).le.0.0_r8) THEN
        print *,'ZDELTA not positive'
        print *, 'ZDELTA = ',Iter, zdelta(Iter)
        stop
      END IF
!

      RETURN
      END SUBROUTINE hessian
!
!***********************************************************************
      SUBROUTINE lanczos (ng, model, Istr, Iend, Jstr, Jend,            &
     &                          LBi, UBi, LBj, UBj,                     &
     &                     Lold, Lnew, Lwrk, Iter, zdelta, zbeta, zqg,  &
     &                     zgnorm,                                      &
#ifdef MASKING
     &                          rmask, umask, vmask,                    &
#endif
#ifdef SOLVE3D
     &                          tl_t, tl_u, tl_v,                       &
#endif
     &                          tl_ubar, tl_vbar, tl_zeta,              &
#ifdef SOLVE3D
     &                          ad_t, ad_u, ad_v,                       &
#endif
     &                          ad_ubar, ad_vbar, ad_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_iounits
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Lold, Lnew, Lwrk, Iter
!
      real(r8), intent(in) :: zdelta(1:Ninner)
      real(r8), intent(inout) :: zgnorm
      real(r8), intent(inout) :: zbeta(1:Ninner+1)
      real(r8), intent(inout) :: zqg(1:Ninner+1)
!
#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
# endif
# ifdef SOLVE3D
      real(r8), intent(inout) :: ad_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_v(LBi:,LBj:,:,:)
# endif
      real(r8), intent(inout) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_vbar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_zeta(LBi:,LBj:,:)
# ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_v(LBi:,LBj:,:,:)
# endif
      real(r8), intent(inout) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_vbar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_zeta(LBi:,LBj:,:)
#else
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
# endif
# ifdef SOLVE3D
      real(r8), intent(inout) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
# endif
      real(r8), intent(inout) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_zeta(LBi:UBi,LBj:UBj,3)
# ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: tl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: tl_v(LBi:UBi,LBj:UBj,N(ng),2)
# endif
      real(r8), intent(inout) :: tl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: tl_vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: tl_zeta(LBi:UBi,LBj:UBj,3)
#endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j, lstr, rec
#ifdef SOLVE3D
      integer :: itrc, k
#endif
      real(r8) :: fac

      real(r8), dimension(0:NstateVar(ng)) :: dot
      real(r8), dimension(Iter) :: DotProd, dot_new, dot_old

      character (len=80) :: ncname

#include "set_bounds.h"
!
!   Apply Lanczos algorithm.
!
   IF (Iter.gt.0) THEN
!
!  Determine adjoint file to process.
!
        IF (ndefADJ(ng).gt.0) THEN
          lstr=LEN_TRIM(ADJbase(ng))
          WRITE (ncname,10) ADJbase(ng)(1:lstr-3), rec
 10       FORMAT (a,'_',i3.3,'.nc')
        ELSE
          ncname=ADJname(ng)
        END IF
!
!  At this point the previous orthonormal Lanczos vector is still in
!  tl_*(Lwrk).
!
#ifndef SOLVE3D
        DO j=JstrR,JendR
          DO i=Istr,IendR
            ad_ubar(i,j,Lnew)=ad_ubar(i,j,Lnew)-                        &
# ifdef MASKING
     &                        umask(i,j)*                               &
# endif
     &                        zdelta(Iter)*tl_ubar(i,j,Lwrk)
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            ad_vbar(i,j,Lnew)=ad_vbar(i,j,Lnew)-                        &
# ifdef MASKING
     &                        vmask(i,j)*                               &
# endif
     &                        zdelta(Iter)*tl_vbar(i,j,Lwrk)
          END DO
        END DO
#endif
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            ad_zeta(i,j,Lnew)=ad_zeta(i,j,Lnew)-                        &
# ifdef MASKING
     &                        rmask(i,j)*                               &
# endif
     &                        zdelta(Iter)*tl_zeta(i,j,Lwrk)
          END DO
        END DO
#ifdef SOLVE3D
!
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=Istr,IendR
              ad_u(i,j,k,Lnew)=ad_u(i,j,k,Lnew)-                        &
# ifdef MASKING
     &                         umask(i,j)*                              &
# endif
     &                         zdelta(Iter)*tl_u(i,j,k,Lwrk)
            END DO
          END DO
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              ad_v(i,j,k,Lnew)=ad_v(i,j,k,Lnew)-                        &
# ifdef MASKING
     &                         vmask(i,j)*                              &
# endif
     &                         zdelta(Iter)*tl_v(i,j,k,Lwrk)
            END DO
          END DO
        END DO
!
        DO itrc=1,NT(ng)
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                ad_t(i,j,k,Lnew,itrc)=ad_t(i,j,k,Lnew,itrc)-            &
# ifdef MASKING
     &                                rmask(i,j)*                       &
# endif
     &                                zdelta(Iter)*tl_t(i,j,k,Lwrk,itrc)
              END DO
            END DO
          END DO
        END DO
      END IF
#endif 
!
!  Read in the previous minus 1 orthonormal Lanczos vector.
!
      IF (Iter.gt.1) THEN
        CALL get_gradient (ng, model, Istr, Iend, Jstr, Jend,           &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     Lwrk, Iter-1, ncname,                          &
#ifdef MASKING
     &                     rmask, umask, vmask,                         &
#endif
#ifdef SOLVE3D
     &                     tl_t, tl_u, tl_v,                            &
#endif
     &                     tl_ubar, tl_vbar, tl_zeta)
!
!
#ifndef SOLVE3D
        DO j=JstrR,JendR
          DO i=Istr,IendR
            ad_ubar(i,j,Lnew)=ad_ubar(i,j,Lnew)-                        &
# ifdef MASKING
     &                        umask(i,j)*                               &
# endif
     &                        zbeta(Iter)*tl_ubar(i,j,Lwrk)
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            ad_vbar(i,j,Lnew)=ad_vbar(i,j,Lnew)-                        &
# ifdef MASKING
     &                        vmask(i,j)*                               &
# endif
     &                        zbeta(Iter)*tl_vbar(i,j,Lwrk)
          END DO
        END DO
#endif
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            ad_zeta(i,j,Lnew)=ad_zeta(i,j,Lnew)-                        &
# ifdef MASKING
     &                        rmask(i,j)*                               &
# endif
     &                        zbeta(Iter)*tl_zeta(i,j,Lwrk)
          END DO
        END DO
#ifdef SOLVE3D
!
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=Istr,IendR
              ad_u(i,j,k,Lnew)=ad_u(i,j,k,Lnew)-                        &
# ifdef MASKING
     &                         umask(i,j)*                              &
# endif
     &                         zbeta(Iter)*tl_u(i,j,k,Lwrk)
            END DO
          END DO
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              ad_v(i,j,k,Lnew)=ad_v(i,j,k,Lnew)-                        &
# ifdef MASKING
     &                         vmask(i,j)*                              &
# endif
     &                         zbeta(Iter)*tl_v(i,j,k,Lwrk)
            END DO
          END DO
        END DO
!
        DO itrc=1,NT(ng)
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                ad_t(i,j,k,Lnew,itrc)=ad_t(i,j,k,Lnew,itrc)-            &
# ifdef MASKING
     &                                rmask(i,j)*                       &
# endif
     &                                zbeta(Iter)*tl_t(i,j,k,Lwrk,itrc)
              END DO
            END DO
          END DO
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Orthogonalize current gradient, G(k+1), against all previous
!  gradients (reverse order) using Gramm-Schmidt procedure.
!-----------------------------------------------------------------------
!
!  We can overwrite adjoint arrays at index Lnew each time around the
!  the following loop because the preceding gradient vectors that we
!  read are orthogonal to each other. The reversed order of the loop
!  is important for the Lanczos vector calculations.
!
      DO rec=Iter,1,-1
!
!  Determine adjoint file to process.
!
        IF (ndefADJ(ng).gt.0) THEN
          lstr=LEN_TRIM(ADJbase(ng))
          WRITE (ncname,10) ADJbase(ng)(1:lstr-3), rec
        ELSE
          ncname=ADJname(ng)
        END IF
!
!  Read in each previous gradient state solutions, G(0) to G(k), and
!  compute its associated dot angaint curret G(k+1). Each gradient
!  solution is loaded into TANGENT LINEAR STATE ARRAYS at index Lwrk.
!
        CALL get_gradient (ng, model, Istr, Iend, Jstr, Jend,           &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     Lwrk, rec, ncname,                           &
#ifdef MASKING
     &                     rmask, umask, vmask,                         &
#endif
#ifdef SOLVE3D
     &                     tl_t, tl_u, tl_v,                            &
#endif
     &                     tl_ubar, tl_vbar, tl_zeta)
!
!  Compute dot product <G(k+1), G(rec)>.
!
        CALL state_dotprod (ng, model, Istr, Iend, Jstr, Jend,          &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NstateVar(ng), dot(0:),                     &
#ifdef MASKING
     &                      rmask, umask, vmask,                        &
#endif
#ifdef SOLVE3D
     &                      ad_t(:,:,:,Lnew,:), tl_t(:,:,:,Lwrk,:),     &
     &                      ad_u(:,:,:,Lnew), tl_u(:,:,:,Lwrk),         &
     &                      ad_v(:,:,:,Lnew), tl_v(:,:,:,Lwrk),         &
#endif
     &                      ad_ubar(:,:,Lnew), tl_ubar(:,:,Lwrk),       &
     &                      ad_vbar(:,:,Lnew), tl_vbar(:,:,Lwrk),       &
     &                      ad_zeta(:,:,Lnew), tl_zeta(:,:,Lwrk))
!
!  Compute Gramm-Schmidt scaling coefficient.
!
        DotProd(rec)=dot(0)
!
!  Gramm-Schmidt orthonormalization, 2D state gradient.
!
#ifndef SOLVE3D
        DO j=JstrR,JendR
          DO i=Istr,IendR
            ad_ubar(i,j,Lnew)=ad_ubar(i,j,Lnew)-                        &
# ifdef MASKING
     &                        umask(i,j)*                               &
# endif
     &                        DotProd(rec)*tl_ubar(i,j,Lwrk)
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            ad_vbar(i,j,Lnew)=ad_vbar(i,j,Lnew)-                        &
# ifdef MASKING
     &                        vmask(i,j)*                               &
# endif
     &                        DotProd(rec)*tl_vbar(i,j,Lwrk)
          END DO
        END DO
#endif
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            ad_zeta(i,j,Lnew)=ad_zeta(i,j,Lnew)-                        &
# ifdef MASKING
     &                        rmask(i,j)*                               &
# endif
     &                        DotProd(rec)*tl_zeta(i,j,Lwrk)
          END DO
        END DO
#ifdef SOLVE3D
!
!  Gramm-Schmidt orthonormalization, 3D state gradient.
!
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=Istr,IendR
              ad_u(i,j,k,Lnew)=ad_u(i,j,k,Lnew)-                        &
# ifdef MASKING
     &                         umask(i,j)*                              &
# endif
     &                         DotProd(rec)*tl_u(i,j,k,Lwrk)
            END DO
          END DO
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              ad_v(i,j,k,Lnew)=ad_v(i,j,k,Lnew)-                        &
# ifdef MASKING
     &                         vmask(i,j)*                              &
# endif
     &                         DotProd(rec)*tl_v(i,j,k,Lwrk)
            END DO
          END DO
        END DO
!
        DO itrc=1,NT(ng)
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                ad_t(i,j,k,Lnew,itrc)=ad_t(i,j,k,Lnew,itrc)-            &
# ifdef MASKING
     &                                rmask(i,j)*                       &
# endif
     &                                DotProd(rec)*tl_t(i,j,k,Lwrk,itrc)
              END DO
            END DO
          END DO
        END DO
#endif
      END DO
!
!-----------------------------------------------------------------------
!  Normalize current orthogonal gradient vector.
!-----------------------------------------------------------------------
!
      CALL state_dotprod (ng, model, Istr, Iend, Jstr, Jend,            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NstateVar(ng), dot(0:),                       &
# ifdef MASKING
     &                    rmask, umask, vmask,                          &
# endif
# ifdef SOLVE3D
     &                    ad_t(:,:,:,Lnew,:), ad_t(:,:,:,Lnew,:),       &
     &                    ad_u(:,:,:,Lnew), ad_u(:,:,:,Lnew),           &
     &                    ad_v(:,:,:,Lnew), ad_v(:,:,:,Lnew),           &
# endif
     &                    ad_ubar(:,:,Lnew), ad_ubar(:,:,Lnew),         &
     &                    ad_vbar(:,:,Lnew), ad_vbar(:,:,Lnew),         &
     &                    ad_zeta(:,:,Lnew), ad_zeta(:,:,Lnew))
!
!  Compute normalization factor.
!
      IF (Iter.eq.0) THEN
         zgnorm=SQRT(dot(0))
         PRINT *, 'ZGNORM = ', zgnorm
      ELSE
         zbeta(Iter+1)=SQRT(dot(0))
      END IF
!
      fac=1.0_r8/SQRT(dot(0))
!
!  Normalize current 2D state gradient vector.
!
# ifndef SOLVE3D
      DO j=JstrR,JendR
        DO i=Istr,IendR
          ad_ubar(i,j,Lnew)=fac*ad_ubar(i,j,Lnew)
        END DO
      END DO
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          ad_vbar(i,j,Lnew)=fac*ad_vbar(i,j,Lnew)
        END DO
      END DO
# endif
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ad_zeta(i,j,Lnew)=fac*ad_zeta(i,j,Lnew)
        END DO
      END DO
# ifdef SOLVE3D
!
!  Normalize current 3D state gradient vector.
!
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=Istr,IendR
            ad_u(i,j,k,Lnew)=fac*ad_u(i,j,k,Lnew)
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            ad_v(i,j,k,Lnew)=fac*ad_v(i,j,k,Lnew)
          END DO
        END DO
      END DO
!
      DO itrc=1,NT(ng)
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              ad_t(i,j,k,Lnew,itrc)=fac*ad_t(i,j,k,Lnew,itrc)
            END DO
          END DO
        END DO
      END DO
# endif
!
!-----------------------------------------------------------------------
!  Compute dot product of new Lanczos vector with gradient.
!-----------------------------------------------------------------------
!
      CALL state_dotprod (ng, model, Istr, Iend, Jstr, Jend,            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NstateVar(ng), dot(0:),                       &
# ifdef MASKING
     &                    rmask, umask, vmask,                          &
# endif
# ifdef SOLVE3D
     &                    ad_t(:,:,:,Lold,:), ad_t(:,:,:,Lnew,:),       &
     &                    ad_u(:,:,:,Lold), ad_u(:,:,:,Lnew),           &
     &                    ad_v(:,:,:,Lold), ad_v(:,:,:,Lnew),           &
# endif
     &                    ad_ubar(:,:,Lold), ad_ubar(:,:,Lnew),         &
     &                    ad_vbar(:,:,Lold), ad_vbar(:,:,Lnew),         &
     &                    ad_zeta(:,:,Lold), ad_zeta(:,:,Lnew))
!
      zqg(Iter+1)=dot(0)
#ifdef TEST_ORTHOGONALIZATION
!
!-----------------------------------------------------------------------
!  Test orthogonal properties of the new gradient.
!-----------------------------------------------------------------------
!
      DO rec=Iter,1,-1
!
!  Determine adjoint file to process.
!
        IF (ndefADJ(ng).gt.0) THEN
          lstr=LEN_TRIM(ADJbase(ng))
          WRITE (ncname,10) ADJbase(ng)(1:lstr-3), rec
        ELSE
          ncname=ADJname(ng)
        END IF
!
!  Read in each previous gradient state solutions, G(0) to G(k), and
!  compute its associated dot angaint orthogonalized G(k+1). Again, 
!  each gradient solution is loaded into TANGENT LINEAR STATE ARRAYS
!  at index Lwrk.
!
        CALL get_gradient (ng, model, Istr, Iend, Jstr, Jend,           &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     Lwrk, rec, ncname,                           &
#ifdef MASKING
     &                     rmask, umask, vmask,                         &
#endif
#ifdef SOLVE3D
     &                     tl_t, tl_u, tl_v,                            &
#endif
     &                     tl_ubar, tl_vbar, tl_zeta)
!
        CALL state_dotprod (ng, model, Istr, Iend, Jstr, Jend,          &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NstateVar(ng), dot(0:),                     &
#ifdef MASKING
     &                      rmask, umask, vmask,                        &
#endif
#ifdef SOLVE3D
     &                      ad_t(:,:,:,Lnew,:), tl_t(:,:,:,Lwrk,:),     &
     &                      ad_u(:,:,:,Lnew), tl_u(:,:,:,Lwrk),         &
     &                      ad_v(:,:,:,Lnew), tl_v(:,:,:,Lwrk),         &
#endif
     &                      ad_ubar(:,:,Lnew), tl_ubar(:,:,Lwrk),       &
     &                      ad_vbar(:,:,Lnew), tl_vbar(:,:,Lwrk),       &
     &                      ad_zeta(:,:,Lnew), tl_zeta(:,:,Lwrk))
        dot_new(rec)=dot(0)
      END DO
!
!  Report dot products. If everything is working correctly, at the
!  end of the orthogonalization dot_new(rec) << dot_old(rec).
!
      IF (Master) THEN
        WRITE (stdout,20) outer, inner
        DO rec=Iter,1,-1
          WRITE (stdout,30) DotProd(rec), rec-1
        END DO
        WRITE (stdout,*) ' '
        DO rec=Iter,1,-1
          WRITE (stdout,40) Iter, rec-1, dot_new(rec),                  &
     &                      rec-1, rec-1, dot_old(rec)
        END DO
 20     FORMAT (/,1x,'(',i3.3,',',i3.3,'): ',                           &
     &          'Gramm-Schmidt Orthogonalization:',/)
 30     FORMAT (12x,'Orthogonalization Factor = ',1p,e19.12,3x,         &
     &          '(Iter=',i3.3,')')
 40     FORMAT (2x,'Ortho Test: ',                                      &
     &          '<G(',i3.3,'),G(',i3.3,')> = ',1p,e15.8,1x,             &
     &          '<G(',i3.3,'),G(',i3.3,')> = ',1p,e15.8)
      END IF
#endif

      RETURN
      END SUBROUTINE lanczos
!
!***********************************************************************
      SUBROUTINE new_gradient (ng, model, Istr, Iend, Jstr, Jend,       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                           Lold, Lnew, Lwrk, Lwrk1, Iter,         &
     &                           zgnorm, zbeta, zwork, zqg,             &
#ifdef MASKING
     &                          rmask, umask, vmask,                    &
#endif
#ifdef SOLVE3D
     &                          tl_t, tl_u, tl_v,                       &
#endif
     &                          tl_ubar, tl_vbar, tl_zeta,              &
#ifdef SOLVE3D
     &                          ad_t, ad_u, ad_v,                       &
#endif
     &                          ad_ubar, ad_vbar, ad_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_iounits
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Lold, Lnew, Lwrk, Lwrk1, Iter
      real(r8), intent(in) :: zgnorm
      real(r8), intent(in) :: zbeta(1:Ninner+1)
      real(r8), intent(in) :: zqg(1:Ninner+1)
      real(r8), intent(in) :: zwork(1:Ninner,3)
!
#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
# endif
# ifdef SOLVE3D
      real(r8), intent(inout) :: ad_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_v(LBi:,LBj:,:,:)
# endif
      real(r8), intent(inout) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_vbar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_zeta(LBi:,LBj:,:)
# ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_v(LBi:,LBj:,:,:)
# endif
      real(r8), intent(inout) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_vbar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_zeta(LBi:,LBj:,:)
#else
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
# endif
# ifdef SOLVE3D
      real(r8), intent(inout) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
# endif
      real(r8), intent(inout) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_zeta(LBi:UBi,LBj:UBj,3)
# ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: tl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: tl_v(LBi:UBi,LBj:UBj,N(ng),2)
# endif
      real(r8), intent(inout) :: tl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: tl_vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: tl_zeta(LBi:UBi,LBj:UBj,3)
#endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j, lstr, rec
#ifdef SOLVE3D
      integer :: itrc, k
#endif
      real(r8) :: fac, preduc

      real(r8), dimension(0:NstateVar(ng)) :: dot
      real(r8), dimension(Iter) :: DotProd, dot_new, dot_old

      character (len=80) :: ncname

#include "set_bounds.h"
!
#ifndef SOLVE3D
        DO j=JstrR,JendR
          DO i=Istr,IendR
            tl_ubar(i,j,Lwrk1)=ad_ubar(i,j,Lwrk1)+                      &
# ifdef MASKING
     &                        umask(i,j)*                               &
# endif
     &                   zbeta(Iter+1)*ad_ubar(i,j,Lnew)*zwork(Iter,3)
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            tl_vbar(i,j,Lwrk1)=ad_vbar(i,j,Lwrk1)+                      &
# ifdef MASKING
     &                        vmask(i,j)*                               &
# endif
     &                   zbeta(Iter+1)*ad_vbar(i,j,Lnew)*zwork(Iter,3)
          END DO
        END DO
#endif
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            tl_zeta(i,j,Lwrk1)=ad_zeta(i,j,Lwrk1)+                      &
# ifdef MASKING
     &                        rmask(i,j)*                               &
# endif
     &                   zbeta(Iter+1)*ad_zeta(i,j,Lnew)*zwork(Iter,3)
          END DO
        END DO
#ifdef SOLVE3D
!
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=Istr,IendR
              tl_u(i,j,k,Lwrk1)=ad_u(i,j,k,Lwrk1)+                      &
# ifdef MASKING
     &                         umask(i,j)*                              &
# endif
     &                   zbeta(Iter+1)*ad_u(i,j,k,Lnew)*zwork(Iter,3)
            END DO
          END DO
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              tl_v(i,j,k,Lwrk1)=ad_v(i,j,k,Lwrk1)+                      &
# ifdef MASKING
     &                         vmask(i,j)*                              &
# endif
     &                   zbeta(Iter+1)*ad_v(i,j,k,Lnew)*zwork(Iter,3)
            END DO
          END DO
        END DO
!
        DO itrc=1,NT(ng)
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                tl_t(i,j,k,Lwrk1,itrc)=ad_t(i,j,k,Lwrk1,itrc)+          &
# ifdef MASKING
     &                                rmask(i,j)*                       &
# endif
     &                   zbeta(Iter+1)*ad_t(i,j,k,Lnew,itrc)*zwork(Iter,3)
              END DO
            END DO
          END DO
        END DO
#endif
!
!
      DO rec=1,Iter
!
!  Determine adjoint file to process.
!
        IF (ndefADJ(ng).gt.0) THEN
          lstr=LEN_TRIM(ADJbase(ng))
          WRITE (ncname,10) ADJbase(ng)(1:lstr-3), rec
 10       FORMAT (a,'_',i3.3,'.nc')
        ELSE
          ncname=ADJname(ng)
        END IF
!
!  Read in each previous gradient state solutions, G(0) to G(k), and
!  compute its associated dot angaint curret G(k+1). Each gradient
!  solution is loaded into TANGENT LINEAR STATE ARRAYS at index Lwrk.
!
        CALL get_gradient (ng, model, Istr, Iend, Jstr, Jend,           &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     Lwrk, rec, ncname,                           &
#ifdef MASKING
     &                     rmask, umask, vmask,                         &
#endif
#ifdef SOLVE3D
     &                     tl_t, tl_u, tl_v,                            &
#endif
     &                     tl_ubar, tl_vbar, tl_zeta)
!
!
#ifndef SOLVE3D
        DO j=JstrR,JendR
          DO i=Istr,IendR
            tl_ubar(i,j,Lwrk1)=tl_ubar(i,j,Lwrk1)-                      &
# ifdef MASKING
     &                        umask(i,j)*                               &
# endif
     &                        zqg(rec)*tl_ubar(i,j,Lwrk)
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            tl_vbar(i,j,Lwrk1)=tl_vbar(i,j,Lwrk1)-                      &
# ifdef MASKING
     &                        vmask(i,j)*                               &
# endif
     &                        zqg(rec)*tl_vbar(i,j,Lwrk)
          END DO
        END DO
#endif
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            tl_zeta(i,j,Lwrk1)=tl_zeta(i,j,Lwrk1)-                      &
# ifdef MASKING
     &                        rmask(i,j)*                               &
# endif
     &                        zqg(rec)*tl_zeta(i,j,Lwrk)
          END DO
        END DO
#ifdef SOLVE3D
!
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=Istr,IendR
              tl_u(i,j,k,Lwrk1)=tl_u(i,j,k,Lwrk1)-                      &
# ifdef MASKING
     &                         umask(i,j)*                              &
# endif
     &                        zqg(rec)*tl_u(i,j,k,Lwrk)
            END DO
          END DO
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              tl_v(i,j,k,Lwrk1)=tl_v(i,j,k,Lwrk1)-                      &
# ifdef MASKING
     &                         vmask(i,j)*                              &
# endif
     &                        zqg(rec)*tl_v(i,j,k,Lwrk)
            END DO
          END DO
        END DO
!
        DO itrc=1,NT(ng)
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                tl_t(i,j,k,Lwrk1,itrc)=tl_t(i,j,k,Lwrk1,itrc)-          &
# ifdef MASKING
     &                                rmask(i,j)*                       &
# endif
     &                                zqg(rec)*tl_t(i,j,k,Lwrk,itrc)
              END DO
            END DO
          END DO
        END DO
#endif
      END DO
!
      CALL state_dotprod (ng, model, Istr, Iend, Jstr, Jend,            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NstateVar(ng), dot(0:),                       &
# ifdef MASKING
     &                    rmask, umask, vmask,                          &
# endif
# ifdef SOLVE3D
     &                    tl_t(:,:,:,Lwrk1,:), tl_t(:,:,:,Lwrk1,:),     &
     &                    tl_u(:,:,:,Lwrk1), tl_u(:,:,:,Lwrk1),         &
     &                    tl_v(:,:,:,Lwrk1), tl_v(:,:,:,Lwrk1),         &
# endif
     &                    tl_ubar(:,:,Lwrk1), tl_ubar(:,:,Lwrk1),       &
     &                    tl_vbar(:,:,Lwrk1), tl_vbar(:,:,Lwrk1),       &
     &                    tl_zeta(:,:,Lwrk1), tl_zeta(:,:,Lwrk1))
!
      preduc=SQRT(dot(0))/zgnorm
!
      print *,'preduc=',preduc
!
#endif

      RETURN
      END SUBROUTINE new_gradient
