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
!  Given an initial model state X(0), gradient G(0), descent direction !
!  d(0), and trial step size tau(1), the minimization algorithm at the !
!  k-iteration is :                                                    !
!                                                                      !
!  (1) Run tangent linear model initialized with trial step, Xhat(k):  !
!                                                                      !
!      Xhat(k) = X(k) + tau(k) * d(k)                          (Eq 5a) !
!                                                                      !
!  (2) Run adjoint model to compute gradient at trial point, Ghat(k):  !
!                                                                      !
!      Ghat(k) = GRAD[ f(Xhat(k)) ]                            (Eq 5b) !
!                                                                      !
!  (3) Compute optimum step size, alpha(k):                            !
!                                                                      !
!      alpha(k) = tau(k) * <d(k),G(k)> / (<d(k),G(k)> - <d(k),Ghat(k)>)!
!                                                                      !
!      here <...> denotes dot product                          (Eq 5c) !
!                                                                      !
!  (4) Compute new starting point (TLM increments), X(k+1):            !
!                                                                      !
!      X(k+1) = X(k) + alpha(k) * d(k)                         (Eq 5d) !
!                                                                      !
!  (5) Compute gradient at new point, G(k+1):                          !
!                                                                      !
!      G(k+1) = G(k) + (alpha(k) / tau(k)) * (Ghat(k) - G(k))  (Eq 5e) !
!                                                                      !
!      overwrite G(k+1) in the NetCDF for latter use.                  !
!                                                                      !
!  (6) Orthogonalize new gradient, G(k+1), against all previous        !
!      gradients [G(k), ..., G(0)], in reverse order, using the        !
!      modified Gramm-Schmidt algorithm. Notice that we need to        !
!      all inner loop gradient solutions.                              !
!                                                                      !
!  (7) Compute new descent direction, d(k+1):                          !
!                                                                      !
!      beta(k+1) = <G(k+1),G(k+1)> / <G(k),G(k)>               (Eq 5g) !
!                                                                      !
!      d(k+1) = - G(k+1) + beta(k+1) * d(k)                    (Eq 5f) !
!                                                                      !
!  After the first iteration, the trial step size is:                  !
!                                                                      !
!      tau(k) = alpha(k-1)                                             !
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
#ifdef ADJUST_WSTRESS
     &                     FORCES(ng) % tl_ustr,                        &
     &                     FORCES(ng) % tl_vstr,                        &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                     FORCES(ng) % tl_tflux,                       &
# endif
     &                     OCEAN(ng) % tl_t,                            &
     &                     OCEAN(ng) % tl_u,                            &
     &                     OCEAN(ng) % tl_v,                            &
#else
     &                     OCEAN(ng) % tl_ubar,                         &
     &                     OCEAN(ng) % tl_vbar,                         &
#endif
     &                     OCEAN(ng) % tl_zeta,                         &
#ifdef ADJUST_WSTRESS
     &                     OCEAN(ng) % d_sustr,                         &
     &                     OCEAN(ng) % d_svstr,                         &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                     OCEAN(ng) % d_stflx,                         &
# endif
     &                     OCEAN(ng) % d_t,                             &
     &                     OCEAN(ng) % d_u,                             &
     &                     OCEAN(ng) % d_v,                             &
#else
     &                     OCEAN(ng) % d_ubar,                          &
     &                     OCEAN(ng) % d_vbar,                          &
#endif
     &                     OCEAN(ng) % d_zeta,                          &
#ifdef ADJUST_WSTRESS
     &                     FORCES(ng) % ad_ustr,                        &
     &                     FORCES(ng) % ad_vstr,                        &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                     FORCES(ng) % ad_stflx,                       &
# endif
     &                     OCEAN(ng) % ad_t,                            &
     &                     OCEAN(ng) % ad_u,                            &
     &                     OCEAN(ng) % ad_v,                            &
#else
     &                     OCEAN(ng) % ad_ubar,                         &
     &                     OCEAN(ng) % ad_vbar,                         &
#endif
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
#ifdef ADJUST_WSTRESS
     &                           tl_ustr, tl_vstr,                      &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                           tl_tflux,                              &
# endif
     &                           tl_t, tl_u, tl_v,                      &
#else
     &                           tl_ubar, tl_vbar,                      &
#endif
     &                           tl_zeta,                               &
#ifdef ADJUST_WSTRESS
     &                           d_sustr, d_svstr,                      &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                           d_stflx,                               &
# endif
     &                           d_t, d_u, d_v,                         &
#else
     &                           d_ubar, d_vbar,                        &
#endif
     &                           d_zeta,                                &
#ifdef ADJUST_WSTRESS
     &                           ad_ustr, ad_vstr,                      &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                           ad_stflx,                              &
# endif
     &                           ad_t, ad_u, ad_v,                      &
#else
     &                           ad_ubar, ad_vbar,                      &
#endif
     &                           ad_zeta)
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
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: ad_ustr(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_vstr(LBi:,LBj:,:)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: ad_stflx(LBi:,LBj:,:,:)
#  endif
      real(r8), intent(inout) :: ad_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_v(LBi:,LBj:,:,:)
# else
      real(r8), intent(inout) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_vbar(LBi:,LBj:,:)
# endif
      real(r8), intent(inout) :: ad_zeta(LBi:,LBj:,:)
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: d_sustr(LBi:,LBj:)
      real(r8), intent(inout) :: d_svstr(LBi:,LBj:)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: d_stflx(LBi:,LBj:,:)
#  endif
      real(r8), intent(inout) :: d_t(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: d_u(LBi:,LBj:,:)
      real(r8), intent(inout) :: d_v(LBi:,LBj:,:)
# else
      real(r8), intent(inout) :: d_ubar(LBi:,LBj:)
      real(r8), intent(inout) :: d_vbar(LBi:,LBj:)
# endif
      real(r8), intent(inout) :: d_zeta(LBi:,LBj:)
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: tl_ustr(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_vstr(LBi:,LBj:,:)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: tl_tflux(LBi:,LBj:,:,:)
#  endif
      real(r8), intent(inout) :: tl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_v(LBi:,LBj:,:,:)
# else
      real(r8), intent(inout) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_vbar(LBi:,LBj:,:)
# endif
      real(r8), intent(inout) :: tl_zeta(LBi:,LBj:,:)
#else
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
# endif
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: ad_ustr(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: ad_vstr(LBi:UBi,LBj:UBj,2)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: ad_stflx(LBi:UBi,LBj:UBj,2,NT(ng))
#  endif
      real(r8), intent(inout) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
# else
      real(r8), intent(inout) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_vbar(LBi:UBi,LBj:UBj,3)
# endif
      real(r8), intent(inout) :: ad_zeta(LBi:UBi,LBj:UBj,3)
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: d_sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: d_svstr(LBi:UBi,LBj:UBj)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: d_stflx(LBi:UBi,LBj:UBj,NT(ng))
#  endif
      real(r8), intent(inout) :: d_t(LBi:UBi,LBj:UBj,N(ng),NT(ng))
      real(r8), intent(inout) :: d_u(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: d_v(LBi:UBi,LBj:UBj,N(ng))
# else
      real(r8), intent(inout) :: d_ubar(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: d_vbar(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(inout) :: d_zeta(LBi:UBi,LBj:UBj)
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: tl_ustr(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: tl_vstr(LBi:UBi,LBj:UBj,2)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: tl_tflux(LBi:UBi,LBj:UBj,2,NT(ng))
#  endif
      real(r8), intent(inout) :: tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: tl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: tl_v(LBi:UBi,LBj:UBj,N(ng),2)
# else
      real(r8), intent(inout) :: tl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: tl_vbar(LBi:UBi,LBj:UBj,3)
# endif
      real(r8), intent(inout) :: tl_zeta(LBi:UBi,LBj:UBj,3)
#endif
!
!  Local variable declarations.
!
      integer :: Linp, Lout, Lwrk, i

      real(r8) :: norm

      real(r8), save :: alphaK, betaK, tauK

      real(r8), dimension(0:NstateVar(ng)) :: Adjust
      real(r8), dimension(0:NstateVar(ng)) :: dot_old, dot_new
      real(r8), dimension(0:NstateVar(ng)) :: old_dot, new_dot
!
!-----------------------------------------------------------------------
!  Initialize trial step size.
!-----------------------------------------------------------------------
!
      IF (Iter.eq.0) THEN
        tauK=CGstepI              ! initial value
        alphaK=tauK
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
!-----------------------------------------------------------------------
!  Compute conjugate gradient optimum step size, alpha(k).
!-----------------------------------------------------------------------
!
      IF (Iter.gt.0) THEN
!
!  Compute old dot product, <d(k), G(k)>.
!
        CALL state_dotprod (ng, model, Istr, Iend, Jstr, Jend,          &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NstateVar(ng), dot_old(0:),                 &
#ifdef MASKING
     &                      rmask, umask, vmask,                        &
#endif
#ifdef ADJUST_WSTRESS
     &                      d_sustr, ad_ustr(:,:,Lold),                 &
     &                      d_svstr, ad_vstr(:,:,Lold),                 &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                      d_stflx, ad_stflx(:,:,:,Lold,:),            &
# endif
     &                      d_t, ad_t(:,:,:,Lold,:),                    &
     &                      d_u, ad_u(:,:,:,Lold),                      &
     &                      d_v, ad_v(:,:,:,Lold),                      &
#else
     &                      d_ubar, ad_ubar(:,:,Lold),                  &
     &                      d_vbar, ad_vbar(:,:,Lold),                  &
#endif
     &                      d_zeta, ad_zeta(:,:,Lold))
!
!  Compute new dot product, <d(k), Ghat(k)>.
!
        CALL state_dotprod (ng, model, Istr, Iend, Jstr, Jend,          &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NstateVar(ng), dot_new(0:),                 &
#ifdef MASKING
     &                      rmask, umask, vmask,                        &
#endif
#ifdef ADJUST_WSTRESS
     &                      d_sustr, ad_ustr(:,:,Lnew),                 &
     &                      d_svstr, ad_vstr(:,:,Lnew),                 &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                      d_stflx, ad_stflx(:,:,Lnew,:),              &
# endif
     &                      d_t, ad_t(:,:,:,Lnew,:),                    &
     &                      d_u, ad_u(:,:,:,Lnew),                      &
     &                      d_v, ad_v(:,:,:,Lnew),                      &
#else
     &                      d_ubar, ad_ubar(:,:,Lnew),                  &
     &                      d_vbar, ad_vbar(:,:,Lnew),                  &
#endif
     &                      d_zeta, ad_zeta(:,:,Lnew))
!
!  Compute new optimal step size.
!
        tauK=alphaK
        alphaK=tauK*dot_old(0)/(dot_old(0)-dot_new(0))
      END IF
!
!  Adjust the cost function for the previous inner-loop iteration.
!  This is based on a first-order Taylor expansion of the cost function.
!  Let vhat=v+tauK*d. During each inner-loop the tangent linear
!  model provides J(vhat). What we require is J(v). Using a 1st-order
!  Taylor expansion we have: J(vhat)=J(v)+tauK*<d,grad> where grad is
!  the cost function gradient computed during the last inner-loop
!  immediately prior to the orthogonalization. Rearranging this
!  equation we have: J(v)=J(vhat)-tauK*<d,grad>. In the code
!  J(vhat)=CostFun(:) and <d,grad>=CostFunDot(:). Remember though
!  that J(v) is the cost function associated with v from the previous
!  inner-loop.
!
      DO i=0,NstateVar(ng)
        Adjust(i)=tauK*FOURDVAR(ng)%CostGradDot(i)
        FOURDVAR(ng)%CostFun(i)=FOURDVAR(ng)%CostFun(i)-Adjust(i)
      END DO
!
!-----------------------------------------------------------------------
!  Estimate the gradient for the new state vector, G(k+1).
!-----------------------------------------------------------------------
!
!  Compute old dot product, <G(k), G(k)>, here since ad_*(Lold) will be
!  used a as temporary storage after this.
!
      CALL state_dotprod (ng, model, Istr, Iend, Jstr, Jend,            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NstateVar(ng), old_dot(0:),                   &
#ifdef MASKING
     &                    rmask, umask, vmask,                          &
#endif
#ifdef ADJUST_WSTRESS
     &                    ad_ustr(:,:,Lold), ad_ustr(:,:,Lold),         &
     &                    ad_vstr(:,:,Lold), ad_vstr(:,:,Lold),         &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                    ad_stflx(:,:,Lold,:), ad_stflx(:,:,Lold,:),   &
# endif
     &                    ad_t(:,:,:,Lold,:), ad_t(:,:,:,Lold,:),       &
     &                    ad_u(:,:,:,Lold), ad_u(:,:,:,Lold),           &
     &                    ad_v(:,:,:,Lold), ad_v(:,:,:,Lold),           &
#else
     &                    ad_ubar(:,:,Lold), ad_ubar(:,:,Lold),         &
     &                    ad_vbar(:,:,Lold), ad_vbar(:,:,Lold),         &
#endif
     &                    ad_zeta(:,:,Lold), ad_zeta(:,:,Lold))
!
!
!  Notice that the current gradient Ghat(k) in time index Lnew is
!  overwritten with the new gradient G(k+1).
!
!    G(k+1) = G(k) + (alpha(k) / tau(k)) * (Ghat(k) - G(k))
!    Lnew     Lold                          Lnew      Lold      index
!
!  Also save G(k+1) in time index Lold as a non-orthogonalized new
!  gradient.
!
      CALL ad_new_state (ng, Istr, Iend, Jstr, Jend,                    &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   Lold, Lnew, alphaK, tauK,                      &
#ifdef MASKING
     &                   rmask, umask, vmask,                           &
#endif
#ifdef ADJUST_WSTRESS
     &                   ad_ustr, ad_vstr,                              &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                   ad_stflx,                                      &
# endif
     &                   ad_t, ad_u, ad_v,                              &
#else
     &                   ad_ubar, ad_vbar,                              &
#endif
     &                   ad_zeta)

#ifdef ORTHOGONALIZATION
!
!  Orthogonalize new gradient, G(k+1), against all previous gradients
!  G(0) to G(k). Use TLM state arrays at time index Lwrk=2, to load
!  each of the previous gradients.
!
      IF (Iter.gt.0) THEN
        Lwrk=2
        CALL orthogonalize (ng, model, Istr, Iend, Jstr, Jend,          &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      Lold, Lnew, Lwrk, Iter,                     &
# ifdef MASKING
     &                      rmask, umask, vmask,                        &
# endif
# ifdef ADJUST_WSTRESS
     &                      tl_ustr, tl_vstr,                           &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                      tl_tflux,                                   &
#  endif
     &                      tl_t, tl_u, tl_v,                           &
# else
     &                      tl_ubar, tl_vbar,                           &
# endif
     &                      tl_zeta,                                    &
# ifdef ADJUST_WSTRESS
     &                      ad_ustr, ad_vstr,                           &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                      ad_stflx,                                   &
#  endif
     &                      ad_t, ad_u, ad_v,                           &
# else
     &                      ad_ubar, ad_vbar,                           &
# endif
     &                      ad_zeta)
      END IF
#endif
!
!-----------------------------------------------------------------------
!  Compute new starting tangent linear state vector, X(k+1).
!-----------------------------------------------------------------------
!
!  Here we are doing step (4), equation 5d, the new TLM increment for
!  the initial conditions are always saved at time level Lout=1.
!
!    X(k+1) = X(k) + alpha(k) * d(k)
!    Lout     Linp                      index
!
      IF (Iter.gt.0) THEN
        Linp=1
        Lout=1
        CALL tl_new_state (ng, Istr, Iend, Jstr, Jend,                  &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     Linp, Lout, alphaK,                          &
#ifdef MASKING
     &                     rmask, umask, vmask,                         &
#endif
#ifdef ADJUST_WSTRESS
     &                     d_sustr, d_svstr,                            &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                     d_stflx,                                     &
# endif
     &                     d_t, d_u, d_v,                               &
#else
     &                     d_ubar, d_vbar,                              &
#endif
     &                     d_zeta,                                      &
#ifdef ADJUST_WSTRESS
     &                     tl_ustr, tl_vstr,                            &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                     tl_tflux,                                    &
# endif
     &                     tl_t, tl_u, tl_v,                            &
#else
     &                     tl_ubar, tl_vbar,                            &
#endif
     &                     tl_zeta)
!
!  If last iteration of inner loop, skip remaining computations. The
!  TLM increments computed here are the ones that are needed update
!  the NLM model initial conditions.
!
!!      IF (Iter.eq.Ninner) RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Compute new conjugate descent direction, d(k+1).
!-----------------------------------------------------------------------
!
      IF (Iter.gt.0) THEN
!
!  Compute new dot product, <G(k+1), G(k+1)>.
!
        CALL state_dotprod (ng, model, Istr, Iend, Jstr, Jend,          &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NstateVar(ng), new_dot(0:),                 &
#ifdef MASKING
     &                      rmask, umask, vmask,                        &
#endif
#ifdef ADJUST_WSTRESS
     &                      ad_ustr(:,:,Lnew), ad_ustr(:,:,Lnew),       &
     &                      ad_vstr(:,:,Lnew), ad_vstr(:,:,Lnew),       &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                      ad_stflx(:,:,Lnew,:), ad_stflx(:,:,Lnew,:), &
# endif
     &                      ad_t(:,:,:,Lnew,:), ad_t(:,:,:,Lnew,:),     &
     &                      ad_u(:,:,:,Lnew), ad_u(:,:,:,Lnew),         &
     &                      ad_v(:,:,:,Lnew), ad_v(:,:,:,Lnew),         &
#else
     &                      ad_ubar(:,:,Lnew), ad_ubar(:,:,Lnew),       &
     &                      ad_vbar(:,:,Lnew), ad_vbar(:,:,Lnew),       &
#endif
     &                      ad_zeta(:,:,Lnew), ad_zeta(:,:,Lnew))
!
!  Compute conjugate direction coefficient, beta(k+1).
!
        betaK=new_dot(0)/old_dot(0)
      ELSE
        betaK=0.0_r8
      END IF
!
!  Compute new conjugate direction, d(k+1).
!
      CALL new_direction (ng, model, Istr, Iend, Jstr, Jend,            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    Lold, Lnew, betaK,                            &
#ifdef MASKING
     &                    rmask, umask, vmask,                          &
#endif
#ifdef ADJUST_WSTRESS
     &                    ad_ustr, ad_vstr,                             &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                    ad_stflx,                                     &
# endif
     &                    ad_t, ad_u, ad_v,                             &
#else
     &                    ad_ubar, ad_vbar,                             &
#endif
     &                    ad_zeta,                                      &
#ifdef ADJUST_WSTRESS
     &                    d_sustr, d_svstr,                             &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                    d_stflx,                                      &
# endif
     &                    d_t, d_u, d_v,                                &
#else
     &                    d_ubar, d_vbar,                               &
#endif
     &                    d_zeta)
!
!  Compute next iteration dot product, <d(k), G(k)>, using new d(k+1)
!  and non-orthogonalized G(k+1) used to adjust cost function.
!
      CALL state_dotprod (ng, model, Istr, Iend, Jstr, Jend,            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NstateVar(ng), FOURDVAR(ng)%CostGradDot(0:),  &
#ifdef MASKING
     &                    rmask, umask, vmask,                          &
#endif
#ifdef ADJUST_WSTRESS
     &                    d_sustr, ad_ustr(:,:,Lold),                   &
     &                    d_svstr, ad_vstr(:,:,Lold),                   &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                    d_stflx, ad_stflx(:,:,:,Lold,:),              &
# endif
     &                    d_t, ad_t(:,:,:,Lold,:),                      &
     &                    d_u, ad_u(:,:,:,Lold),                        &
     &                    d_v, ad_v(:,:,:,Lold),                        &
#else
     &                    d_ubar, ad_ubar(:,:,Lold),                    &
     &                    d_vbar, ad_vbar(:,:,Lold),                    &
#endif
     &                    d_zeta, ad_zeta(:,:,Lold))
!
!-----------------------------------------------------------------------
!  Set TLM initial conditions for next inner loop, Xhat(k+1).
!-----------------------------------------------------------------------
!
!  Here we are doing step (1), equation 5a, the new TLM initial
!  conditions for the next inner loop are always saved at Lout=2.
!
!    Xhat(k+1) = X(k+1) + tau(k+1) * d(k+1),  where  tau(k+1)=alpha(k)
!    Lout        Linp                         index
!
      Linp=1
      Lout=2
      CALL tl_new_state (ng, Istr, Iend, Jstr, Jend,                    &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   Linp, Lout, alphaK,                            &
#ifdef MASKING
     &                   rmask, umask, vmask,                           &
#endif
#ifdef ADJUST_WSTRESS
     &                   d_sustr, d_svstr,                              &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                   d_stflx,                                       &
# endif
     &                   d_t, d_u, d_v,                                 &
#else
     &                   d_ubar, d_vbar,                                &
#endif
     &                   d_zeta,                                        &
#ifdef ADJUST_WSTRESS
     &                   tl_ustr, tl_vstr,                              &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                   tl_tflux,                                      &
# endif
     &                   tl_t, tl_u, tl_v,                              &
#else
     &                   tl_ubar, tl_vbar,                              &
#endif
     &                   tl_zeta)
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
      SUBROUTINE tl_new_state (ng, Istr, Iend, Jstr, Jend,              &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         Linp, Lout, alphaK,                      &
#ifdef MASKING
     &                         rmask, umask, vmask,                     &
#endif
#ifdef ADJUST_WSTRESS
     &                         d_sustr, d_svstr,                        &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                         d_stflx,                                 &
# endif
     &                         d_t, d_u, d_v,                           &
#else
     &                         d_ubar, d_vbar,                          &
#endif
     &                         d_zeta,                                  &
#ifdef ADJUST_WSTRESS
     &                         tl_ustr, tl_vstr,                        &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                         tl_tflux,                                &
# endif
     &                         tl_t, tl_u, tl_v,                        &
#else
     &                         tl_ubar, tl_vbar,                        &
#endif
     &                         tl_zeta)
!***********************************************************************
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Linp, Lout

      real(r8), intent(in) :: alphaK
!
#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
# endif
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: d_sustr(LBi:,LBj:)
      real(r8), intent(inout) :: d_svstr(LBi:,LBj:)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: d_stflx(LBi:,LBj:,:)
#  endif
      real(r8), intent(inout) :: d_t(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: d_u(LBi:,LBj:,:)
      real(r8), intent(inout) :: d_v(LBi:,LBj:,:)
# else
      real(r8), intent(inout) :: d_ubar(LBi:,LBj:)
      real(r8), intent(inout) :: d_vbar(LBi:,LBj:)
# endif
      real(r8), intent(inout) :: d_zeta(LBi:,LBj:)
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: tl_ustr(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_vstr(LBi:,LBj:,:)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: tl_tflux(LBi:,LBj:,:,:)
#  endif
      real(r8), intent(inout) :: tl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_v(LBi:,LBj:,:,:)
# else
      real(r8), intent(inout) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_vbar(LBi:,LBj:,:)
# endif
      real(r8), intent(inout) :: tl_zeta(LBi:,LBj:,:)
#else
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
# endif
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: d_sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: d_svstr(LBi:UBi,LBj:UBj)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: d_stflx(LBi:UBi,LBj:UBj,NT(ng))
#  endif
      real(r8), intent(inout) :: d_t(LBi:UBi,LBj:UBj,N(ng),NT(ng))
      real(r8), intent(inout) :: d_u(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: d_v(LBi:UBi,LBj:UBj,N(ng))
# else
      real(r8), intent(inout) :: d_ubar(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: d_vbar(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(inout) :: d_zeta(LBi:UBi,LBj:UBj)
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: tl_ustr(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: tl_vstr(LBi:UBi,LBj:UBj,2)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: tl_tflux(LBi:UBi,LBj:UBj,2,NT(ng))
#  endif
      real(r8), intent(inout) :: tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: tl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: tl_v(LBi:UBi,LBj:UBj,N(ng),2)
# else
      real(r8), intent(inout) :: tl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: tl_vbar(LBi:UBi,LBj:UBj,3)
# endif
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
!  Free-surface.
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          tl_zeta(i,j,Lout)=tl_zeta(i,j,Linp)+                          &
     &                      alphaK*d_zeta(i,j)
#ifdef MASKING
          tl_zeta(i,j,Lout)=tl_zeta(i,j,Lout)*rmask(i,j)
#endif
        END DO
      END DO

#ifndef SOLVE3D
!
!  2D momentum.
!
      DO j=JstrR,JendR
        DO i=Istr,IendR
          tl_ubar(i,j,Lout)=tl_ubar(i,j,Linp)+                          &
     &                      alphaK*d_ubar(i,j)
# ifdef MASKING
          tl_ubar(i,j,Lout)=tl_ubar(i,j,Lout)*umask(i,j)
# endif
        END DO
      END DO
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          tl_vbar(i,j,Lout)=tl_vbar(i,j,Linp)+                          &
     &                      alphaK*d_vbar(i,j)
# ifdef MASKING
          tl_vbar(i,j,Lout)=tl_vbar(i,j,Lout)*vmask(i,j)
# endif
        END DO
      END DO
#endif
#ifdef ADJUST_WSTRESS
!
!  Surface momentum stress.
!
      DO j=JstrR,JendR
        DO i=Istr,IendR
          tl_ustr(i,j,Lout)=tl_ustr(i,j,Linp)+                          &
     &                      alphaK*d_sustr(i,j)
# ifdef MASKING
          tl_ustr(i,j,Lout)=tl_ustr(i,j,Lout)*umask(i,j)
# endif
        END DO
      END DO
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          tl_vstr(i,j,Lout)=tl_vstr(i,j,Linp)+                          &
     &                      alphaK*d_vbar(i,j)
# ifdef MASKING
          tl_vstr(i,j,Lout)=tl_vstr(i,j,Lout)*vmask(i,j)
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
            tl_u(i,j,k,Lout)=tl_u(i,j,k,Linp)+                          &
     &                       alphaK*d_u(i,j,k)
# ifdef MASKING
            tl_u(i,j,k,Lout)=tl_u(i,j,k,Lout)*umask(i,j)
# endif
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            tl_v(i,j,k,Lout)=tl_v(i,j,k,Linp)+                          &
     &                       alphaK*d_v(i,j,k)
# ifdef MASKING
            tl_v(i,j,k,Lout)=tl_v(i,j,k,Lout)*vmask(i,j)
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
              tl_t(i,j,k,Lout,itrc)=tl_t(i,j,k,Linp,itrc)+              &
     &                              alphaK*d_t(i,j,k,itrc)
# ifdef MASKING
              tl_t(i,j,k,Lout,itrc)=tl_t(i,j,k,Lout,itrc)*rmask(i,j)
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
            tl_tflux(i,j,Lout,itrc)=tl_tflux(i,j,Linp,itrc)+            &
     &                              alphaK*d_t(i,j,k,itrc)
#  ifdef MASKING
            tl_tflux(i,j,Lout,itrc)=tl_tflux(i,j,Lout,itrc)*rmask(i,j)
#  endif
          END DO          
        END DO
      END DO
# endif
#endif

      RETURN
      END SUBROUTINE tl_new_state
!
!***********************************************************************
      SUBROUTINE ad_new_state (ng, Istr, Iend, Jstr, Jend,              &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         Lold, Lnew, alphaK, tauK,                &
#ifdef MASKING
     &                         rmask, umask, vmask,                     &
#endif
#ifdef ADJUST_WSTRESS
     &                         ad_ustr, ad_vstr,                        &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                         ad_stflx,                                &
# endif
     &                         ad_t, ad_u, ad_v,                        &
#else
     &                         ad_ubar, ad_vbar,                        &
#endif
     &                         ad_zeta)
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
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: ad_ustr(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_vstr(LBi:,LBj:,:)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: ad_stflx(LBi:,LBj:,:,:)
#  endif
      real(r8), intent(inout) :: ad_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_v(LBi:,LBj:,:,:)
# else
      real(r8), intent(inout) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_vbar(LBi:,LBj:,:)
# endif
      real(r8), intent(inout) :: ad_zeta(LBi:,LBj:,:)
#else
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
# endif
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: ad_ustr(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: ad_vstr(LBi:UBi,LBj:UBj,2)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: ad_stflx(LBi:UBi,LBj:UBj,2,NT(ng))
#  endif
      real(r8), intent(inout) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
# else
      real(r8), intent(inout) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_vbar(LBi:UBi,LBj:UBj,3)
# endif
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
!  Free-surface.
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ad_zeta(i,j,Lnew)=ad_zeta(i,j,Lold)+                          &
     &                      fac*(ad_zeta(i,j,Lnew)-                     &
     &                           ad_zeta(i,j,Lold))
#ifdef MASKING
          ad_zeta(i,j,Lnew)=ad_zeta(i,j,Lnew)*rmask(i,j)
#endif
          ad_zeta(i,j,Lold)=ad_zeta(i,j,Lnew)
        END DO
      END DO

#ifndef SOLVE3D
!
!  2D momentum.
!
      DO j=JstrR,JendR
        DO i=Istr,IendR
          ad_ubar(i,j,Lnew)=ad_ubar(i,j,Lold)+                          &
     &                      fac*(ad_ubar(i,j,Lnew)-                     &
     &                           ad_ubar(i,j,Lold))
# ifdef MASKING
          ad_ubar(i,j,Lnew)=ad_ubar(i,j,Lnew)*umask(i,j)
# endif
          ad_ubar(i,j,Lold)=ad_ubar(i,j,Lnew)
        END DO
      END DO
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          ad_vbar(i,j,Lnew)=ad_vbar(i,j,Lold)+                          &
     &                      fac*(ad_vbar(i,j,Lnew)-                     &
     &                           ad_vbar(i,j,Lold))
# ifdef MASKING
          ad_vbar(i,j,Lnew)=ad_vbar(i,j,Lnew)*vmask(i,j)
# endif
          ad_vbar(i,j,Lold)=ad_vbar(i,j,Lnew)
        END DO
      END DO
#endif

#ifdef ADJUST_WSTRESS
!
!  Surface momentum stress.
!
      DO j=JstrR,JendR
        DO i=Istr,IendR
          ad_ustr(i,j,Lnew)=ad_ustr(i,j,Lold)+                          &
     &                      fac*(ad_ustr(i,j,Lnew)-                     &
     &                           ad_ustr(i,j,Lold))
# ifdef MASKING
          ad_ustr(i,j,Lnew)=ad_ustr(i,j,Lnew)*umask(i,j)
# endif
          ad_ustr(i,j,Lold)=ad_ustr(i,j,Lnew)
        END DO
      END DO
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          ad_vstr(i,j,Lnew)=ad_vstr(i,j,Lold)+                          &
     &                      fac*(ad_vstr(i,j,Lnew)-                     &
     &                          ad_vstr(i,j,Lold))
# ifdef MASKING
          ad_vstr(i,j,Lnew)=ad_vstr(i,j,Lnew)*vmask(i,j)
# endif
          ad_vstr(i,j,Lold)=ad_vstr(i,j,Lnew)
        END DO
      END DO
#endif
#ifdef SOLVE3D
!
!  3D state variables.
!
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=Istr,IendR
            ad_u(i,j,k,Lnew)=ad_u(i,j,k,Lold)+                          &
     &                       fac*(ad_u(i,j,k,Lnew)-                     &
     &                            ad_u(i,j,k,Lold))
# ifdef MASKING
            ad_u(i,j,k,Lnew)=ad_u(i,j,k,Lnew)*umask(i,j)
# endif
            ad_u(i,j,k,Lold)=ad_u(i,j,k,Lnew)
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            ad_v(i,j,k,Lnew)=ad_v(i,j,k,Lold)+                          &
     &                       fac*(ad_v(i,j,k,Lnew)-                     &
     &                            ad_v(i,j,k,Lold))
# ifdef MASKING
            ad_v(i,j,k,Lnew)=ad_v(i,j,k,Lnew)*vmask(i,j)
# endif
            ad_v(i,j,k,Lold)=ad_v(i,j,k,Lnew)
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
              ad_t(i,j,k,Lnew,itrc)=ad_t(i,j,k,Lold,itrc)+              &
     &                              fac*(ad_t(i,j,k,Lnew,itrc)-         &
     &                                   ad_t(i,j,k,Lold,itrc))
# ifdef MASKING
              ad_t(i,j,k,Lnew,itrc)=ad_t(i,j,k,Lnew,itrc)*rmask(i,j)
# endif
              ad_t(i,j,k,Lold,itrc)=ad_t(i,j,k,Lnew,itrc)
            END DO
          END DO
        END DO
      END DO

# ifdef ADJUST_STFLUX
!
!  Tracers.
!
      DO itrc=1,NT(ng)
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            ad_stflx(i,j,Lnew,itrc)=ad_stflx(i,j,Lold,itrc)+            &
     &                              fac*(ad_stflx(i,j,Lnew,itrc)-       &
     &                                   ad_stflx(i,j,Lold,itrc))
#  ifdef MASKING
            ad_stflx(i,j,Lnew,itrc)=ad_stflx(i,j,Lnew,itrc)*rmask(i,j)
#  endif
            ad_stflx(i,j,Lold,itrc)=ad_stflx(i,j,Lnew,itrc)
          END DO
        END DO
      END DO
# endif
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
#ifdef ADJUST_WSTRESS
     &                          tl_ustr, tl_vstr,                       &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                          tl_tflux,                               &
# endif
     &                          tl_t, tl_u, tl_v,                       &
#else
     &                          tl_ubar, tl_vbar,                       &
#endif
     &                          tl_zeta,                                &
#ifdef ADJUST_WSTRESS
     &                          ad_ustr, ad_vstr,                       &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                          ad_stflx,                               &
# endif
     &                          ad_t, ad_u, ad_v,                       &
#else
     &                          ad_ubar, ad_vbar,                       &
#endif
     &                          ad_zeta)
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
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: ad_ustr(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_vstr(LBi:,LBj:,:)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: ad_stflx(LBi:,LBj:,:,:)
#  endif
      real(r8), intent(inout) :: ad_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_v(LBi:,LBj:,:,:)
# else
      real(r8), intent(inout) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_vbar(LBi:,LBj:,:)
# endif
      real(r8), intent(inout) :: ad_zeta(LBi:,LBj:,:)
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: tl_ustr(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_vstr(LBi:,LBj:,:)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: tl_tflux(LBi:,LBj:,:,:)
#  endif
      real(r8), intent(inout) :: tl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_v(LBi:,LBj:,:,:)
# else
      real(r8), intent(inout) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_vbar(LBi:,LBj:,:)
# endif
      real(r8), intent(inout) :: tl_zeta(LBi:,LBj:,:)
#else
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
# endif
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: ad_ustr(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: ad_vstr(LBi:UBi,LBj:UBj,2)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: ad_stflx(LBi:UBi,LBj:UBj,2,NT(ng))
#  endif
      real(r8), intent(inout) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
# else
      real(r8), intent(inout) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_vbar(LBi:UBi,LBj:UBj,3)
# endif
      real(r8), intent(inout) :: ad_zeta(LBi:UBi,LBj:UBj,3)
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: tl_ustr(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: tl_vstr(LBi:UBi,LBj:UBj,2)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: tl_tflux(LBi:UBi,LBj:UBj,2,NT(ng))
#  endif
      real(r8), intent(inout) :: tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: tl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: tl_v(LBi:UBi,LBj:UBj,N(ng),2)
# else
      real(r8), intent(inout) :: tl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: tl_vbar(LBi:UBi,LBj:UBj,3)
# endif
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
      real(r8) :: fac, fac1, fac2

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
#ifdef ADJUST_WSTRESS
     &                     tl_ustr, tl_vstr,                            &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                     tl_stlfx,                                    &
# endif
     &                     tl_t, tl_u, tl_v,                            &
#else
     &                     tl_ubar, tl_vbar,                            &
#endif
     &                     tl_zeta)
!
!  Compute dot product <G(k+1), G(rec)>.
!
        CALL state_dotprod (ng, model, Istr, Iend, Jstr, Jend,          &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NstateVar(ng), dot(0:),                     &
#ifdef MASKING
     &                      rmask, umask, vmask,                        &
#endif
#ifdef ADJUST_WSTRESS
     &                      ad_ustr(:,:,Lnew), tl_ustr(:,:,Lwrk),       &
     &                      ad_vstr(:,:,Lnew), tl_vstr(:,:,Lwrk),       &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                      ad_stflx(:,:,Lnew,:), tl_tflux(:,:,Lwrk,:), &
# endif
     &                      ad_t(:,:,:,Lnew,:), tl_t(:,:,:,Lwrk,:),     &
     &                      ad_u(:,:,:,Lnew), tl_u(:,:,:,Lwrk),         &
     &                      ad_v(:,:,:,Lnew), tl_v(:,:,:,Lwrk),         &
#else
     &                      ad_ubar(:,:,Lnew), tl_ubar(:,:,Lwrk),       &
     &                      ad_vbar(:,:,Lnew), tl_vbar(:,:,Lwrk),       &
#endif
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
#ifdef ADJUST_WSTRESS
     &                      tl_ustr(:,:,Lwrk), tl_ustr(:,:,Lwrk),       &
     &                      tl_vstr(:,:,Lwrk), tl_vstr(:,:,Lwrk),       &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                      tl_tflux(:,:,Lwrk,:), tl_tflux(:,:,Lwrk,:), &
# endif
     &                      tl_t(:,:,:,Lwrk,:), tl_t(:,:,:,Lwrk,:),     &
     &                      tl_u(:,:,:,Lwrk), tl_u(:,:,:,Lwrk),         &
     &                      tl_v(:,:,:,Lwrk), tl_v(:,:,:,Lwrk),         &
#else
     &                      tl_ubar(:,:,Lwrk), tl_ubar(:,:,Lwrk),       &
     &                      tl_vbar(:,:,Lwrk), tl_vbar(:,:,Lwrk),       &
#endif
     &                      tl_zeta(:,:,Lwrk), tl_zeta(:,:,Lwrk))
        dot_old(rec)=dot(0)
!
!  Compute Gramm-Schmidt scaling coefficient.
!
        DotProd(rec)=dot_new(rec)/dot_old(rec)
!
!  Gramm-Schmidt orthonormalization, 2D state gradient.
!
        fac1=1.0_r8
        fac2=-DotProd(rec)

        CALL state_addition (ng, Istr, Iend, Jstr, Jend,                &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       Lnew, Lwrk, Lnew, fac1, fac2,              &
#ifdef MASKING
     &                       rmask, umask, vmask,                       &
#endif
#ifdef ADJUST_WSTRESS
     &                       ad_ustr, tl_vstr,                          &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                       ad_stflx, tl_tflux,                        &
# endif
     &                       ad_t, tl_t,                                &
     &                       ad_u, tl_u,                                &
     &                       ad_v, tl_v,                                &
#else
     &                       ad_ubar, tl_ubar,                          &
     &                       ad_vbar, tl_vbar,                          &
#endif
     &                       ad_zeta, tl_zeta)
      END DO

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
#ifdef ADJUST_WSTRESS
     &                     tl_ustr, tl_vstr,                            &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                     tl_stlfx,                                    &
# endif
     &                     tl_t, tl_u, tl_v,                            &
#else
     &                     tl_ubar, tl_vbar,                            &
#endif
     &                     tl_zeta)
!
        CALL state_dotprod (ng, model, Istr, Iend, Jstr, Jend,          &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NstateVar(ng), dot(0:),                     &
#ifdef MASKING
     &                      rmask, umask, vmask,                        &
#endif
#ifdef ADJUST_WSTRESS
     &                      ad_ustr(:,:,Lnew), tl_ustr(:,:,Lwrk),       &
     &                      ad_vstr(:,:,Lnew), tl_vstr(:,:,Lwrk),       &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                      ad_stflx(:,:,Lnew,:), tl_tflux(:,:,Lwrk,:), &
# endif
     &                      ad_t(:,:,:,Lnew,:), tl_t(:,:,:,Lwrk,:),     &
     &                      ad_u(:,:,:,Lnew), tl_u(:,:,:,Lwrk),         &
     &                      ad_v(:,:,:,Lnew), tl_v(:,:,:,Lwrk),         &
#else
     &                      ad_ubar(:,:,Lnew), tl_ubar(:,:,Lwrk),       &
     &                      ad_vbar(:,:,Lnew), tl_vbar(:,:,Lwrk),       &
#endif
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
#ifdef ADJUST_WSTRESS
     &                         s_ustr, s_vstr,                          &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                         s_tflux,                                 &
# endif
     &                         s_t, s_u, s_v,                           &
#else
     &                         s_ubar, s_vbar,                          &
#endif
     &                         s_zeta)
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
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: s_ustr(LBi:,LBj:,:)
      real(r8), intent(inout) :: s_vstr(LBi:,LBj:,:)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: s_tflux(LBi:,LBj:,:,:)
#  endif
      real(r8), intent(inout) :: s_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: s_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: s_v(LBi:,LBj:,:,:)
# else
      real(r8), intent(inout) :: s_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: s_vbar(LBi:,LBj:,:)
# endif
      real(r8), intent(inout) :: s_zeta(LBi:,LBj:,:)
#else
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
# endif
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: s_ustr(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: s_vstr(LBi:UBi,LBj:UBj,2)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: s_tflux(LBi:UBi,LBj:UBj,2,NT(ng))
#  endif
      real(r8), intent(inout) :: s_t(LBi:UBi,LBj:UBj,N(ng),NT(ng))
      real(r8), intent(inout) :: s_u(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: s_v(LBi:UBi,LBj:UBj,N(ng))
# else
      real(r8), intent(inout) :: s_ubar(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: s_vbar(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(inout) :: s_zeta(LBi:UBi,LBj:UBj)
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
#ifdef ADJUST_WSTRESS
          status=nf_inq_varid(ncid, TRIM(Vname(1,idUsms)), vid(idUsms))
          status=nf_inq_varid(ncid, TRIM(Vname(1,idVsms)), vid(idVsms))
#endif
#ifdef SOLVE3D
          status=nf_inq_varid(ncid, TRIM(Vname(1,idUvel)), vid(idUvel))
          status=nf_inq_varid(ncid, TRIM(Vname(1,idVvel)), vid(idVvel))
          DO itrc=1,NT(ng)
            status=nf_inq_varid(ncid, TRIM(Vname(1,idTvar(itrc))),      &
     &                          vid(idTvar(itrc)))
# ifdef ADJUST_STFLUX
            status=nf_inq_varid(ncid, TRIM(Vname(1,idTsur(itrc))),      &
     &                          vid(idTsur(itrc)))
# endif
          END DO
#else
          status=nf_inq_varid(ncid, TRIM(Vname(1,idUbar)), vid(idUbar))
          status=nf_inq_varid(ncid, TRIM(Vname(1,idVbar)), vid(idVbar))
#endif
          status=nf_inq_varid(ncid, TRIM(Vname(1,idFsur)), vid(idFsur))
        END IF
      ELSE
        ncid=ncADJid(ng)
#ifdef ADJUST_WSTRESS
        vid(idUsms)=adjVid(idUsms,ng)
        vid(idVsms)=adjVid(idVsms,ng)
#endif
#ifdef SOLVE3D
        vid(idUvel)=adjVid(idUvel,ng)
        vid(idVvel)=adjVid(idVvel,ng)
        DO itrc=1,NT(ng)
          vid(idTvar(itrc))=adjTid(itrc,ng)
# ifdef ADJUST_STFLUX
          vid(idTsur(itrc))=adjVid(idTsur(itrc),ng)
# endif
        END DO
#else
        vid(idUbar)=adjVid(idUbar,ng)
        vid(idVbar)=adjVid(idVbar,ng)
#endif
        vid(idFsur)=adjVid(idFsur,ng)
      END IF
      DO i=1,4
        Vsize(i)=0
      END DO
      scale=1.0_r8
!
!  Read in free-surface.
!
      status=nf_fread2d(ng, iTLM, ncid, vid(idFsur), rec, r2dvar,       &
     &                  Vsize, LBi, UBi, LBj, UBj,                      &
     &                  scale, Fmin, Fmax,                              &
#ifdef MASKING
     &                  rmask(LBi,LBj),                                 &
#endif
     &                  s_zeta(LBi,LBj,Lwrk))
      IF (status.ne.nf_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idFsur)), rec, TRIM(ncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF

#ifndef SOLVE3D
!
!  Read in 2D momentum.
!
      status=nf_fread2d(ng, iTLM, ncid, vid(idUbar), rec, u2dvar,       &
     &                  Vsize, LBi, UBi, LBj, UBj,                      &
     &                  scale, Fmin, Fmax,                              &
# ifdef MASKING
     &                  umask(LBi,LBj),                                 &
# endif
     &                  s_ubar(LBi,LBj,Lwrk))
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
     &                  s_vbar(LBi,LBj,Lwrk))
      IF (status.ne.nf_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idVbar)), rec, TRIM(ncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
#endif

#ifdef ADJUST_WSTRESS
!
!  Read in surface momentum stress.
!
      status=nf_fread2d(ng, iTLM, ncid, vid(idUsms), rec, u2dvar,       &
     &                  Vsize, LBi, UBi, LBj, UBj,                      &
     &                  scale, Fmin, Fmax,                              &
# ifdef MASKING
     &                  umask(LBi,LBj),                                 &
# endif
     &                  s_ustr(LBi,LBj,Lwrk))
      IF (status.ne.nf_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idUsms)), rec, TRIM(ncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
      status=nf_fread2d(ng, iTLM, ncid, vid(idVsms), rec, v2dvar,       &
     &                  Vsize, LBi, UBi, LBj, UBj,                      &
     &                  scale, Fmin, Fmax,                              &
# ifdef MASKING
     &                  vmask(LBi,LBj),                                 &
# endif
     &                  s_vstr(LBi,LBj,Lwrk))
      IF (status.ne.nf_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idVsms)), rec, TRIM(ncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
#endif

#ifdef SOLVE3D
!
!  Read in 3D momentum.
!
      status=nf_fread3d(ng, iTLM, ncid, vid(idUvel), rec, u3dvar,       &
     &                  Vsize, LBi, UBi, LBj, UBj, 1, N(ng),            &
     &                  scale, Fmin, Fmax,                              &
# ifdef MASKING
     &                  umask(LBi,LBj),                                 &
# endif
     &                  s_u(LBi,LBj,1,Lwrk))
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
     &                  s_v(LBi,LBj,1,Lwrk))
      IF (status.ne.nf_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idVvel)), rec, TRIM(ncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Read in tracers.
!
      DO itrc=1,NT(ng)
        status=nf_fread3d(ng, iTLM, ncid, vid(idTvar(itrc)), rec,       &
     &                    r3dvar, Vsize, LBi, UBi, LBj, UBj, 1, N(ng),  &
     &                    scale, Fmin, Fmax,                            &
# ifdef MASKING
     &                    rmask(LBi,LBj),                               &
# endif
     &                    s_t(LBi,LBj,1,Lwrk,itrc))
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

# ifdef ADJUST_STFLUX
!
!  Read in surface tracers flux.
!
      DO itrc=1,NT(ng)
        status=nf_fread2d(ng, iTLM, ncid, vid(idTsur(itrc)), rec,       &
     &                    r3dvar, Vsize, LBi, UBi, LBj, UBj,            &
     &                    scale, Fmin, Fmax,                            &
#  ifdef MASKING
     &                    rmask(LBi,LBj),                               &
#  endif
     &                    s_tflux(LBi,LBj,Lwrk,itrc))
        IF (status.ne.nf_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idTsur(itrc))), rec,         &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END DO
# endif
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
#ifdef ADJUST_WSTRESS
     &                          ad_ustr, ad_vstr,                       &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                          ad_stflx,                               &
# endif
     &                          ad_t, ad_u, ad_v,                       &
#else
     &                          ad_ubar, ad_vbar,                       &
#endif
     &                          ad_zeta,                                &
#ifdef ADJUST_WSTRESS
     &                          d_sustr, d_svstr,                       &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                          d_stflx,                                &
# endif
     &                          d_t, d_u, d_v,                          &
#else
     &                          d_ubar, d_vbar,                         &
#endif
     &                          d_zeta)
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
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: ad_ustr(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_vstr(LBi:,LBj:,:)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: ad_stflx(LBi:,LBj:,:,:)
#  endif
      real(r8), intent(inout) :: ad_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_v(LBi:,LBj:,:,:)
# else
      real(r8), intent(inout) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_vbar(LBi:,LBj:,:)
# endif
      real(r8), intent(inout) :: ad_zeta(LBi:,LBj:,:)
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: d_sustr(LBi:,LBj:)
      real(r8), intent(inout) :: d_svstr(LBi:,LBj:)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: d_stflx(LBi:,LBj:,:)
#  endif
      real(r8), intent(inout) :: d_t(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: d_u(LBi:,LBj:,:)
      real(r8), intent(inout) :: d_v(LBi:,LBj:,:)
# else
      real(r8), intent(inout) :: d_ubar(LBi:,LBj:)
      real(r8), intent(inout) :: d_vbar(LBi:,LBj:)
# endif
      real(r8), intent(inout) :: d_zeta(LBi:,LBj:)
#else
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
# endif
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: ad_ustr(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: ad_vstr(LBi:UBi,LBj:UBj,2)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: ad_stflx(LBi:UBi,LBj:UBj,2,NT(ng))
#  endif
      real(r8), intent(inout) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
# else
      real(r8), intent(inout) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_vbar(LBi:UBi,LBj:UBj,3)
# endif
      real(r8), intent(inout) :: ad_zeta(LBi:UBi,LBj:UBj,3)
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: d_sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: d_svstr(LBi:UBi,LBj:UBj)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: d_stflx(LBi:UBi,LBj:UBj,NT(ng))
#  endif
      real(r8), intent(inout) :: d_t(LBi:UBi,LBj:UBj,N(ng),NT(ng))
      real(r8), intent(inout) :: d_u(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: d_v(LBi:UBi,LBj:UBj,N(ng))
# else
      real(r8), intent(inout) :: d_ubar(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: d_vbar(LBi:UBi,LBj:UBj)
# endif
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
!  Free-surface.
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          d_zeta(i,j)=-ad_zeta(i,j,Lnew)+betaK*d_zeta(i,j)
#ifdef MASKING
          d_zeta(i,j)=d_zeta(i,j)*rmask(i,j)
#endif
        END DO
      END DO

#ifndef SOLVE3D
!
!  2D momentum.
!
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

#ifdef ADJUST_WSTRESS
!
!  Surface momentum stress.
!
      DO j=JstrR,JendR
        DO i=Istr,IendR
          d_sustr(i,j)=-ad_ustr(i,j,Lnew)+betaK*d_sustr(i,j)
# ifdef MASKING
          d_sustr(i,j)=d_sustr(i,j)*umask(i,j)
# endif
        END DO
      END DO
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          d_svstr(i,j)=-ad_vstr(i,j,Lnew)+betaK*d_svstr(i,j)
# ifdef MASKING
          d_svstr(i,j)=d_svstr(i,j)*vmask(i,j)
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
!  Tracers.
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

# ifdef ADJUST_STFLUX
!
!  Surface tracers flux.
!
      DO itrc=1,NT(ng)
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            d_stflx(i,j,itrc)=-ad_stflx(i,j,Lnew,itrc)+                 &
     &                        betaK*d_stflx(i,j,itrc)
#  ifdef MASKING
            d_stflx(i,j,itrc)=d_stflx(i,j,itrc)*rmask(i,j)
#  endif
          END DO
        END DO
      END DO
# endif
#endif

      RETURN
      END SUBROUTINE new_direction
