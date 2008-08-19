!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module minimizes  incremental  4Dvar  quadratic cost function  !
!  using a preconditioned version of the conjugate gradient algorithm  !
!  proposed by Mike Fisher (ECMWF) and modified  by  Tshimanga et al.  !
!  (2008).                                                             !
!                                                                      !
!  In the following,  M represents the preconditioner.  Specifically,  !
!                                                                      !
!    M = I + SUM_i [ (mu_i-1) h_i (h_i)'],                             !
!                                                                      !
!  where mu_i can take the following values:                           !
!                                                                      !
!    Lscale=-1:    mu_i = lambda_i                                     !
!    Lscale= 1:    mu_i = 1 / lambda_i                                 !
!    Lscale=-2:    mu_i = SQRT (lambda_i)                              !
!    Lscale= 2:    mu_i = 1 / SQRT(lambda_i)                           !
!                                                                      !
!  where lambda_i are the Hessian eigenvalues and h_i are the Hessian  !
!  eigenvectors.                                                       !
!                                                                      !
!  For Lscale= 1 spectral LMP is used as the preconditioner.           !
!  For Lscale=-1 inverse spectral LMP is used as the preconditioner.   !
!  For Lscale= 2 SQRT spectral LMP is used as the preconditioner.      !
!  For Lscale=-2 inverse SQRT spectral LMP is the preconditioner.      !
!                                                                      !
!  If Lritz=.TRUE. then Ritz LMP is used and the expressions for       !
!   mu_i are more complicated.                                         !
!                                                                      !
!  For some operations the tranpose of the preconditioner is require.  !
!  For spectral LMP the preconditioner and its tranpose are identical. !
!  For Ritz LMP the preconditioner and its tranpose differ.            !
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
!      alpha(k) = tau(k) * <d(k), G(k)> /                              !
!                                                                      !
!                 (<d(k),G(k)> - <d(k), Ghat(k)>)              (Eq 5c) !
!                                                                      !
!      here <...> denotes dot product                                  !
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
!      save all inner loop gradient solutions.                         !
!                                                                      !
!      For the preconditioned case the appropriate inner-product       !
!      for the orthonormalizatio is <G,G>.                             !
!                                                                      !
!  (7) Compute new descent direction, d(k+1):                          !
!                                                                      !
!      beta(k+1) = <G(k+1), G(k+1)> / <G(k), G(k)>             (Eq 5g) !
!                                                                      !
!      d(k+1) = - G(k+1) + beta(k+1) * d(k)                    (Eq 5f) !
!                                                                      !
!  After the first iteration, the trial step size is:                  !
!                                                                      !
!      tau(k) = alpha(k-1)                                             !
!                                                                      !
!  On entering cgradient, all variable are in v-space, the variable    !
!  results from the first level of preconditioning. Within cgradient   !
!  the minimization will be performed in y-space where v=My, and M is  !
!  the square root spectral preconditioner. The congrad algorithm is   !
!  is the same in y-space, but various transformations from v-space    !
!  to y-space and vice-versa are required.                             !
!                                                                      !
!  NOTE: In all of the following computations we are using the NLM     !
!        state variable arrays as temporary arrays.                    !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    Fisher, M., 1997: Efficient Minimization of Quadratic Penalty     !
!      funtions, unpublish manuscript, 1-14.                           !
!                                                                      !
!    Tshimanga, J., S. Gratton, A. Weaver, and A. Sartenaer, 2008:     !
!      Limited-memory preconditioners with application to              !
!      incremental four-dimensional ocean data assimilation,           !
!      Q. J. R. Meteorol. Soc., 134, 753-771.                          !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC :: cgradient

      CONTAINS
!
!***********************************************************************
      SUBROUTINE cgradient (ng, tile, model, innLoop, outLoop)
!***********************************************************************
!
      USE mod_param
#ifdef SOLVE3D
      USE mod_coupling
#endif
#if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      USE mod_forces
#endif
      USE mod_grid
      USE mod_ocean
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model, innLoop, outLoop
!
!  Local variable declarations.
!
#include "tile.h"
!
#ifdef PROFILE
      CALL wclock_on (ng, model, 36)
#endif
      CALL cgradient_tile (ng, tile, model,                             &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     Lold(ng), Lnew(ng),                          &
     &                     innLoop, outLoop,                            &
#ifdef MASKING
     &                     GRID(ng) % rmask,                            &
     &                     GRID(ng) % umask,                            &
     &                     GRID(ng) % vmask,                            &
#endif
#ifdef ADJUST_WSTRESS
     &                     FORCES(ng) % ustr,                           &
     &                     FORCES(ng) % vstr,                           &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                     FORCES(ng) % tflux,                          &
# endif
     &                     OCEAN(ng) % t,                               &
     &                     OCEAN(ng) % u,                               &
     &                     OCEAN(ng) % v,                               &
#else
     &                     OCEAN(ng) % ubar,                            &
     &                     OCEAN(ng) % vbar,                            &
#endif
     &                     OCEAN(ng) % zeta,                            &
#ifdef ADJUST_WSTRESS
     &                     FORCES(ng) % d_sustr,                        &
     &                     FORCES(ng) % d_svstr,                        &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                     FORCES(ng) % d_stflx,                        &
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
     &                     FORCES(ng) % ad_tflux,                       &
# endif
     &                     OCEAN(ng) % ad_t,                            &
     &                     OCEAN(ng) % ad_u,                            &
     &                     OCEAN(ng) % ad_v,                            &
#else
     &                     OCEAN(ng) % ad_ubar,                         &
     &                     OCEAN(ng) % ad_vbar,                         &
#endif
     &                     OCEAN(ng) % ad_zeta,                         &
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
     &                     OCEAN(ng) % tl_zeta)
#ifdef PROFILE
      CALL wclock_on (ng, model, 36)
#endif
      RETURN
      END SUBROUTINE cgradient
!
!***********************************************************************
      SUBROUTINE cgradient_tile (ng, tile, model,                       &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           Lold, Lnew,                            &
     &                           innLoop, outLoop,                      &
#ifdef MASKING
     &                           rmask, umask, vmask,                   &
#endif
#ifdef ADJUST_WSTRESS
     &                           nl_ustr, nl_vstr,                      &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                           nl_tflux,                              &
# endif
     &                           nl_t, nl_u, nl_v,                      &
#else
     &                           nl_ubar, nl_vbar,                      &
#endif
     &                           nl_zeta,                               &
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
     &                           ad_tflux,                              &
# endif
     &                           ad_t, ad_u, ad_v,                      &
#else
     &                           ad_ubar, ad_vbar,                      &
#endif
     &                           ad_zeta,                               &
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
     &                           tl_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
!
      USE state_copy_mod, ONLY : state_copy
      USE state_dotprod_mod, ONLY : state_dotprod
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Lold, Lnew
      integer, intent(in) :: innLoop, outLoop
!
#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
# endif
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: ad_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_vstr(LBi:,LBj:,:,:)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: ad_tflux(LBi:,LBj:,:,:,:)
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
      real(r8), intent(inout) :: d_sustr(LBi:,LBj:,:)
      real(r8), intent(inout) :: d_svstr(LBi:,LBj:,:)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: d_stflx(LBi:,LBj:,:,:)
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
      real(r8), intent(inout) :: nl_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: nl_vstr(LBi:,LBj:,:,:)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: nl_tflux(LBi:,LBj:,:,:,:)
#  endif
      real(r8), intent(inout) :: nl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: nl_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: nl_v(LBi:,LBj:,:,:)
# else
      real(r8), intent(inout) :: nl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: nl_vbar(LBi:,LBj:,:)
# endif
      real(r8), intent(inout) :: nl_zeta(LBi:,LBj:,:)
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: tl_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_vstr(LBi:,LBj:,:,:)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: tl_tflux(LBi:,LBj:,:,:,:)
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
      real(r8), intent(inout) :: ad_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: ad_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: ad_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
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
      real(r8), intent(inout) :: d_sustr(LBi:UBi,LBj:UBj,Nfrec(ng))
      real(r8), intent(inout) :: d_svstr(LBi:UBi,LBj:UBj,Nfrec(ng))
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: d_stflx(LBi:UBi,LBj:UBj,               &
     &                                   Nfrec(ng),NT(ng))
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
      real(r8), intent(inout) :: nl_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: nl_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: nl_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
#  endif
      real(r8), intent(inout) :: nl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: nl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: nl_v(LBi:UBi,LBj:UBj,N(ng),2)
# else
      real(r8), intent(inout) :: nl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: nl_vbar(LBi:UBi,LBj:UBj,3)
# endif
      real(r8), intent(inout) :: nl_zeta(LBi:UBi,LBj:UBj,3)
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: tl_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: tl_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: tl_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
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
      logical :: Lritz, Ltrans

      integer :: Linp, Lout, Lwrk, L1, L2, i, Lscale
      integer :: status, varid
      integer :: start(4), total(4)

      real(r8) :: norm

      real(r8), dimension(0:NstateVar(ng)) :: Adjust
      real(r8), dimension(0:NstateVar(ng)) :: dot_old, dot_new
      real(r8), dimension(0:NstateVar(ng)) :: old_dot, new_dot

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Initialize trial step size.
!-----------------------------------------------------------------------
!
      Lritz=.TRUE.
      Ltrans=.FALSE.
!
      IF (innLoop.eq.0) THEN
        cg_tau(innLoop,outLoop)=CGstepI
        cg_alpha(innLoop,outLoop)=cg_tau(innLoop,outLoop)
        DO i=0,NstateVar(ng)
          dot_old(i)=0.0_r8
          dot_new(i)=0.0_r8
          old_dot(i)=0.0_r8
          new_dot(i)=0.0_r8
        END DO
      END IF
      IF (Master) THEN
        WRITE (stdout,10)
 10     FORMAT (/,' <<<< Descent Algorithm >>>>')
      END IF
!
!  If preconditioning, read in number of converged eigenvectors and their
!  associated eigenvalues.
!
      IF (Lprecond.and.((innLoop.eq.0).and.(outLoop.eq.1))) THEN

        CALL netcdf_get_ivar (ng, model, HSSname(ng), 'nConvRitz',      &
     &                        nConvRitz,                                &
     &                        ncid = ncHSSid(ng))
        IF (exit_flag.ne.NoError) RETURN

        CALL netcdf_get_fvar (ng, model, HSSname(ng), 'Ritz',           &
     &                        Ritz,                                     &
     &                        ncid = ncHSSid(ng),                       &
     &                        start = (/1/), total = (/nConvRitz/))
        IF (exit_flag.ne. NoError) RETURN
!
!  Reset number of eigenpairs to use to specified value.
!
        nConvRitz=NritzEV

      END IF
!
      IF (Lprecond) THEN
!
!  Convert ad_var(Lnew) from v-space to y-space.
!
        Lscale=2
        Lwrk=1
        Ltrans=.TRUE.
!
!  Copy ad_var(Lnew) into nl_var(1)
!
        CALL state_copy (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   Lnew, Lwrk,                                    &
#ifdef ADJUST_WSTRESS
     &                   nl_ustr, ad_ustr,                              &
     &                   nl_vstr, ad_vstr,                              &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                   nl_tflux, ad_tflux,                            &
# endif
     &                   nl_t, ad_t,                                    &
     &                   nl_u, ad_u,                                    &
     &                   nl_v, ad_v,                                    &
#else
     &                   nl_ubar, ad_ubar,                              &
     &                   nl_vbar, ad_vbar,                              &
#endif
     &                   nl_zeta, ad_zeta)
!
        CALL precond (ng, tile, model,                                  &
     &                LBi, UBi, LBj, UBj,                               &
     &                IminS, ImaxS, JminS, JmaxS,                       &
     &                innLoop, outLoop,                                 &
     &                NstateVar(ng), Lscale,                            &
     &                Lritz, Ltrans,                                    &
     &                nConvRitz, Ritz,                                  &
#ifdef MASKING
     &                rmask, umask, vmask,                              &
#endif
#ifdef ADJUST_WSTRESS
     &                nl_ustr, nl_vstr,                                 &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                nl_tflux,                                         &
# endif
     &                nl_t, nl_u, nl_v,                                 &
#else
     &                nl_ubar, nl_vbar,                                 &
#endif
     &                nl_zeta)
        IF (exit_flag.ne.NoError) RETURN
!
!  Copy nl_var(1) into ad_var(Lnew)
!
        Lwrk=1
        CALL state_copy (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   Lwrk, Lnew,                                    &
#ifdef ADJUST_WSTRESS
     &                   ad_ustr, nl_ustr,                              &
     &                   ad_vstr, nl_vstr,                              &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                   ad_tflux, nl_tflux,                            &
# endif
     &                   ad_t, nl_t,                                    &
     &                   ad_u, nl_u,                                    &
     &                   ad_v, nl_v,                                    &
#else
     &                   ad_ubar, nl_ubar,                              &
     &                   ad_vbar, nl_vbar,                              &
#endif
     &                   ad_zeta, nl_zeta)

      END IF
!
!-----------------------------------------------------------------------
!  Compute conjugate gradient optimum step size, alpha(k).
!-----------------------------------------------------------------------
!
      IF (innLoop.gt.0) THEN
        CALL state_dotprod (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NstateVar(ng), dot_old(0:),                 &
#ifdef MASKING
     &                      rmask, umask, vmask,                        &
#endif
#ifdef ADJUST_WSTRESS
     &                      d_sustr, ad_ustr(:,:,:,Lold),               &
     &                      d_svstr, ad_vstr(:,:,:,Lold),               &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                      d_stflx, ad_tflux(:,:,:,Lold,:),            &
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
        CALL state_dotprod (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NstateVar(ng), dot_new(0:),                 &
#ifdef MASKING
     &                      rmask, umask, vmask,                        &
#endif
#ifdef ADJUST_WSTRESS
     &                      d_sustr, ad_ustr(:,:,:,Lnew),               &
     &                      d_svstr, ad_vstr(:,:,:,Lnew),               &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                      d_stflx, ad_tflux(:,:,:,Lnew,:),            &
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
        cg_tau(innLoop,outLoop)=cg_alpha(innLoop-1,outLoop)
        cg_alpha(innLoop,outLoop)=cg_tau(innLoop,outLoop)*              &
     &                            dot_old(0)/(dot_old(0)-dot_new(0))
      END IF
!
!-----------------------------------------------------------------------
!  Estimate the gradient for the new state vector, G(k+1).
!-----------------------------------------------------------------------
!
!  Compute old dot product, <G(k), G(k)>. The ADM arrays, index Lold,
!  will be used a as temporary storage after this.
!
        CALL state_dotprod (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NstateVar(ng), old_dot(0:),                 &
#ifdef MASKING
     &                      rmask, umask, vmask,                        &
#endif
#ifdef ADJUST_WSTRESS
     &                      ad_ustr(:,:,:,Lold), ad_ustr(:,:,:,Lold),   &
     &                      ad_vstr(:,:,:,Lold), ad_vstr(:,:,:,Lold),   &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                      ad_tflux(:,:,:,Lold,:),                     &
     &                      ad_tflux(:,:,:,Lold,:),                     &
# endif
     &                      ad_t(:,:,:,Lold,:), ad_t(:,:,:,Lold,:),     &
     &                      ad_u(:,:,:,Lold), ad_u(:,:,:,Lold),         &
     &                      ad_v(:,:,:,Lold), ad_v(:,:,:,Lold),         &
#else
     &                      ad_ubar(:,:,Lold), ad_ubar(:,:,Lold),       &
     &                      ad_vbar(:,:,Lold), ad_vbar(:,:,Lold),       &
#endif
     &                      ad_zeta(:,:,Lold), ad_zeta(:,:,Lold))
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
      CALL ad_new_state (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   Lold, Lnew,                                    &
     &                   cg_alpha(innLoop,outLoop),                     &
     &                   cg_tau(innLoop,outLoop),                       &
#ifdef MASKING
     &                   rmask, umask, vmask,                           &
#endif
#ifdef ADJUST_WSTRESS
     &                   ad_ustr, ad_vstr,                              &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                   ad_tflux,                                      &
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
!  NOTE: All of the gradient vectors saved in the adjoint netcdf
!        are in y-space.
!
      IF (innLoop.gt.0) THEN
        Lwrk=2
        CALL orthogonalize (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      Lold, Lnew, Lwrk,                           &
     &                      innLoop, outLoop,                           &
# ifdef MASKING
     &                      rmask, umask, vmask,                        &
# endif
# ifdef ADJUST_WSTRESS
     &                      nl_ustr, nl_vstr,                           &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                      nl_tflux,                                   &
#  endif
     &                      nl_t, nl_u, nl_v,                           &
# else
     &                      nl_ubar, nl_vbar,                           &
# endif
     &                      nl_zeta,                                    &
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
     &                      ad_tflux,                                   &
#  endif
     &                      ad_t, ad_u, ad_v,                           &
# else
     &                      ad_ubar, ad_vbar,                           &
# endif
     &                      ad_zeta)
        IF (exit_flag.ne.NoError) RETURN
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
!  If preconditioning, convert tl_var(Linp) from v-space into y-space.
!
      IF (innLoop.gt.0) THEN
        IF (Lprecond) THEN
          Lscale=-2
          Linp=1
          Lwrk=1
          Ltrans=.FALSE.
!
!  Copy tl_var(Linp) into nl_var(1)
!
          CALL state_copy (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     Linp, Lwrk,                                  &
#ifdef ADJUST_WSTRESS
     &                     nl_ustr, tl_ustr,                            &
     &                     nl_vstr, tl_vstr,                            &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                     nl_tflux, tl_tflux,                          &
# endif
     &                     nl_t, tl_t,                                  &
     &                     nl_u, tl_u,                                  &
     &                     nl_v, tl_v,                                  &
#else
     &                     nl_ubar, tl_ubar,                            &
     &                     nl_vbar, tl_vbar,                            &
#endif
     &                     nl_zeta, tl_zeta)
!
          CALL precond (ng, tile, model,                                &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  IminS, ImaxS, JminS, JmaxS,                     &
     &                  innLoop, outLoop,                               &
     &                  NstateVar(ng), Lscale,                          &
     &                  Lritz, Ltrans,                                  &
     &                  nConvRitz, Ritz,                                &
#ifdef MASKING 
     &                  rmask, umask, vmask,                            &
#endif 
#ifdef ADJUST_WSTRESS
     &                  nl_ustr, nl_vstr,                               &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                  nl_tflux,                                       &
# endif
     &                  nl_t, nl_u, nl_v,                               &
#else
     &                  nl_ubar, nl_vbar,                               &
#endif
     &                  nl_zeta)
          IF (exit_flag.ne.NoError) RETURN
!
!  Copy nl_var(1) into tl_var(Linp)
!
          Lwrk=1
          CALL state_copy (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     Lwrk, Linp,                                  &
#ifdef ADJUST_WSTRESS
     &                     tl_ustr, nl_ustr,                            &
     &                     tl_vstr, nl_vstr,                            &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                     tl_tflux, nl_tflux,                          &
# endif
     &                     tl_t, nl_t,                                  &
     &                     tl_u, nl_u,                                  &
     &                     tl_v, nl_v,                                  &
#else
     &                     tl_ubar, nl_ubar,                            &
     &                     tl_vbar, nl_vbar,                            &
#endif
     &                     tl_zeta, nl_zeta)

        END IF
!
        Linp=1
        Lout=1
!
!  tl_var(Lout) is now in y-space if preconditioning.
!
        CALL tl_new_state (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     Linp, Lout,                                  &
     &                     cg_alpha(innLoop,outLoop),                   &
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
!  Compute the new cost function.
!
        CALL new_cost (ng, tile, model,                                 &
     &                 LBi, UBi, LBj, UBj,                              &
     &                 IminS, ImaxS, JminS, JmaxS,                      &
     &                 Lout,                                            &
#ifdef MASKING
     &                 rmask, umask, vmask,                             &
#endif
#ifdef ADJUST_WSTRESS
     &                 nl_ustr, nl_vstr,                                &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                 nl_tflux,                                        &
# endif
     &                 nl_t, nl_u, nl_v,                                &
#else
     &                 nl_ubar, nl_vbar,                                &
#endif
     &                 nl_zeta,                                         &
#ifdef ADJUST_WSTRESS
     &                 tl_ustr, tl_vstr,                                &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                 tl_tflux,                                        &
# endif
     &                 tl_t, tl_u, tl_v,                                &
#else
     &                 tl_ubar, tl_vbar,                                &
#endif
     &                 tl_zeta) 
        IF (exit_flag.ne.NoError) RETURN
!
!  If last iteration of inner loop, skip remaining computations. The
!  TLM increments computed here are the ones that are needed update
!  the NLM model initial conditions.
!
!!      IF (innLoop.eq.Ninner) RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Compute new conjugate descent direction, d(k+1).
!-----------------------------------------------------------------------
!
      IF (innLoop.gt.0) THEN
!
!  Compute new dot product, <G(k+1), G(k+1)>.
!
        CALL state_dotprod (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NstateVar(ng), new_dot(0:),                 &
#ifdef MASKING
     &                      rmask, umask, vmask,                        &
#endif
#ifdef ADJUST_WSTRESS
     &                      ad_ustr(:,:,:,Lnew), ad_ustr(:,:,:,Lnew),   &
     &                      ad_vstr(:,:,:,Lnew), ad_vstr(:,:,:,Lnew),   &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                      ad_tflux(:,:,:,Lnew,:),                     &
     &                      ad_tflux(:,:,:,Lnew,:),                     &
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
        cg_beta(innLoop,outLoop)=new_dot(0)/old_dot(0)
      ELSE
        cg_beta(innLoop,outLoop)=0.0_r8
      END IF
!
!  Compute new conjugate direction, d(k+1).
!
      CALL new_direction (ng, tile, model,                              &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    Lold, Lnew,                                   &
     &                    cg_beta(innLoop,outLoop),                     &
#ifdef MASKING
     &                    rmask, umask, vmask,                          &
#endif
#ifdef ADJUST_WSTRESS
     &                    ad_ustr, ad_vstr,                             &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                    ad_tflux,                                     &
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
!-----------------------------------------------------------------------
!  Set TLM initial conditions for next inner loop, Xhat(k+1).
!-----------------------------------------------------------------------
!
!  Here we are doing step (1), equation 5a, the new TLM initial
!  conditions for the next inner loop are always saved at Lout=2.
!
!  NOTE: If preconditioning, then X(k+1) is already in y-space at this
!        point.
!
!    Xhat(k+1) = X(k+1) + tau(k+1) * d(k+1),  where  tau(k+1)=alpha(k)
!    Lout        Linp                         index
!
      Linp=1
      Lout=2
      CALL tl_new_state (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   Linp, Lout,                                    &
     &                   cg_alpha(innLoop,outLoop),                     &
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
!  Convert the new tl_var(Lout) and tl_var(Linp) back to v-space.
!  Note that we leave ad_var(Lnew) in y-space since we want all of the
!  gradient vectors in the adjoint netcdf file to be saved in y-space
!  to avoid costly applications of the spectral LMP in the
!  routine orthogonalize.
!
      IF (Lprecond) THEN
        Lscale=2
        Lwrk=1
        Ltrans=.FALSE.

!
!  Copy tl_var(Lout) into nl_var(1)
!
        CALL state_copy (ng, tile,                                      &
     &                 LBi, UBi, LBj, UBj,                              &
     &                 Lout, Lwrk,                                      &
#ifdef ADJUST_WSTRESS
     &                 nl_ustr, tl_ustr,                                &
     &                 nl_vstr, tl_vstr,                                &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                 nl_tflux, tl_tflux,                              &
# endif
     &                 nl_t, tl_t,                                      &
     &                 nl_u, tl_u,                                      &
     &                 nl_v, tl_v,                                      &
#else
     &                 nl_ubar, tl_ubar,                                &
     &                 nl_vbar, tl_vbar,                                &
#endif
     &                 nl_zeta, tl_zeta)
!
        CALL precond (ng, tile, model,                                  &
     &                LBi, UBi, LBj, UBj,                               &
     &                IminS, ImaxS, JminS, JmaxS,                       &
     &                innLoop, outLoop,                                 &
     &                NstateVar(ng), Lscale,                            &
     &                Lritz, Ltrans,                                    &
     &                nConvRitz, Ritz,                                  &
#ifdef MASKING
     &                rmask, umask, vmask,                              &
#endif
#ifdef ADJUST_WSTRESS
     &                nl_ustr, nl_vstr,                                 &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                nl_tflux,                                         &
# endif
     &                nl_t, nl_u, nl_v,                                 &
#else
     &                nl_ubar, nl_vbar,                                 &
#endif
     &                nl_zeta)
        IF (exit_flag.ne.NoError) RETURN
!
!  Copy nl_var(1) into tl_var(Lout)
!
        Lwrk=1
        CALL state_copy (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   Lwrk, Lout,                                    &
#ifdef ADJUST_WSTRESS
     &                   tl_ustr, nl_ustr,                              &
     &                   tl_vstr, nl_vstr,                              &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                   tl_tflux, nl_tflux,                            &
# endif
     &                   tl_t, nl_t,                                    &
     &                   tl_u, nl_u,                                    &
     &                   tl_v, nl_v,                                    &
#else
     &                   tl_ubar, nl_ubar,                              &
     &                   tl_vbar, nl_vbar,                              &
#endif
     &                   tl_zeta, nl_zeta)
! 
!  Copy tl_var(Linp) into nl_var(1)
!
        CALL state_copy (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   Linp, Lwrk,                                    &
#ifdef ADJUST_WSTRESS
     &                   nl_ustr, tl_ustr,                              &
     &                   nl_vstr, tl_vstr,                              &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                   nl_tflux, tl_tflux,                            &
# endif
     &                   nl_t, tl_t,                                    &
     &                   nl_u, tl_u,                                    &
     &                   nl_v, tl_v,                                    &
#else
     &                   nl_ubar, tl_ubar,                              &
     &                   nl_vbar, tl_vbar,                              &
#endif
     &                   nl_zeta, tl_zeta)
!
        CALL precond (ng, tile, model,                                  &
     &                LBi, UBi, LBj, UBj,                               &
     &                IminS, ImaxS, JminS, JmaxS,                       &
     &                innLoop, outLoop,                                 &
     &                NstateVar(ng), Lscale,                            &
     &                Lritz, Ltrans,                                    &
     &                nConvRitz, Ritz,                                  &
#ifdef MASKING
     &                rmask, umask, vmask,                              &
#endif
#ifdef ADJUST_WSTRESS
     &                nl_ustr, nl_vstr,                                 &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                nl_tflux,                                         &
# endif
     &                nl_t, nl_u, nl_v,                                 &
#else
     &                nl_ubar, nl_vbar,                                 &
#endif
     &                nl_zeta)
        IF (exit_flag.ne.NoError) RETURN
!
!  Copy nl_var(1) into tl_var(Linp)
!
        Lwrk=1
        CALL state_copy (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   Lwrk, Linp,                                    &
#ifdef ADJUST_WSTRESS
     &                   tl_ustr, nl_ustr,                              &
     &                   tl_vstr, nl_vstr,                              &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                   tl_tflux, nl_tflux,                            &
# endif
     &                   tl_t, nl_t,                                    &
     &                   tl_u, nl_u,                                    &
     &                   tl_v, nl_v,                                    &
#else
     &                   tl_ubar, nl_ubar,                              &
     &                   tl_vbar, nl_vbar,                              &
#endif
     &                   tl_zeta, nl_zeta)
! 
      END IF
!
!-----------------------------------------------------------------------
!  Write out conjugate gradient information into NetCDF file.
!-----------------------------------------------------------------------
!
      CALL cg_write (ng, model, innLoop, outLoop)
      IF (exit_flag.ne.NoError) RETURN
!
!  Report  algorithm parameters.
!
      IF (Master) THEN
        WRITE (stdout,30) outLoop, innLoop,                             &
     &                    cg_tau(innLoop,outLoop),                      &
     &                    cg_alpha(innLoop,outLoop),                    &
     &                    cg_beta(innLoop,outLoop),                     &
     &                    outLoop, MAX(0,innLoop-1), Adjust(0),         &
     &                    outLoop, innLoop,                             &
     &                    'dot product', innLoop, innLoop,              &
     &                    dot_old(0), 'alpha',                          &
     &                    'dot product', innLoop, innLoop,              &
     &                    dot_new(0), 'alpha',                          &
     &                    'dot product', innLoop, innLoop,              &
     &                    old_dot(0), 'beta',                           &
     &                    'dot product', innLoop+1, innLoop+1,          &
     &                     new_dot(0), 'beta'
 30     FORMAT (/,1x,'(',i3.3,',',i3.3,'): ',                           &
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
      SUBROUTINE tl_new_state (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
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
#if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      USE mod_scalars
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
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
      real(r8), intent(inout) :: d_sustr(LBi:,LBj:,:)
      real(r8), intent(inout) :: d_svstr(LBi:,LBj:,:)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: d_stflx(LBi:,LBj:,:,:)
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
      real(r8), intent(inout) :: tl_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_vstr(LBi:,LBj:,:,:)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: tl_tflux(LBi:,LBj:,:,:,:)
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
      real(r8), intent(inout) :: d_sustr(LBi:UBi,LBj:UBj,Nfrec(ng))
      real(r8), intent(inout) :: d_svstr(LBi:UBi,LBj:UBj,Nfrec(ng))
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: d_stflx(LBi:UBi,LBj:UBj,               &
     &                                   Nfrec(ng),NT(ng))
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
      real(r8), intent(inout) :: tl_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: tl_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: tl_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
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
      integer :: i, j, k
#ifdef SOLVE3D
      integer :: itrc
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
      DO k=1,Nfrec(ng)
        DO j=JstrR,JendR
          DO i=Istr,IendR
            tl_ustr(i,j,k,Lout)=tl_ustr(i,j,k,Linp)+                    &
     &                          alphaK*d_sustr(i,j,k)
# ifdef MASKING
            tl_ustr(i,j,k,Lout)=tl_ustr(i,j,k,Lout)*umask(i,j)
# endif
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            tl_vstr(i,j,k,Lout)=tl_vstr(i,j,k,Linp)+                    &
     &                         alphaK*d_svstr(i,j,k)
# ifdef MASKING
            tl_vstr(i,j,k,Lout)=tl_vstr(i,j,k,Lout)*vmask(i,j)
# endif
          END DO
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
        DO k=1,Nfrec(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              tl_tflux(i,j,k,Lout,itrc)=tl_tflux(i,j,k,Linp,itrc)+      &
     &                                  alphaK*d_stflx(i,j,k,itrc)
#  ifdef MASKING
              tl_tflux(i,j,k,Lout,itrc)=tl_tflux(i,j,k,Lout,itrc)*      &
     &                                  rmask(i,j)
#  endif
            END DO
          END DO
        END DO
      END DO
# endif
#endif

      RETURN
      END SUBROUTINE tl_new_state
!
!***********************************************************************
      SUBROUTINE ad_new_state (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         Lold, Lnew, alphaK, tauK,                &
#ifdef MASKING
     &                         rmask, umask, vmask,                     &
#endif
#ifdef ADJUST_WSTRESS
     &                         ad_ustr, ad_vstr,                        &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                         ad_tflux,                                &
# endif
     &                         ad_t, ad_u, ad_v,                        &
#else
     &                         ad_ubar, ad_vbar,                        &
#endif
     &                         ad_zeta)
!***********************************************************************
!
      USE mod_param
#if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      USE mod_scalars
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
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
      real(r8), intent(inout) :: ad_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_vstr(LBi:,LBj:,:,:)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: ad_tflux(LBi:,LBj:,:,:,:)
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
      real(r8), intent(inout) :: ad_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: ad_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: ad_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
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
      integer :: i, j, k
#ifdef SOLVE3D
      integer :: itrc
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
      DO k=1,Nfrec(ng)
        DO j=JstrR,JendR
          DO i=Istr,IendR
            ad_ustr(i,j,k,Lnew)=ad_ustr(i,j,k,Lold)+                    &
     &                          fac*(ad_ustr(i,j,k,Lnew)-               &
     &                               ad_ustr(i,j,k,Lold))
# ifdef MASKING
            ad_ustr(i,j,k,Lnew)=ad_ustr(i,j,k,Lnew)*umask(i,j)
# endif
            ad_ustr(i,j,k,Lold)=ad_ustr(i,j,k,Lnew)
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            ad_vstr(i,j,k,Lnew)=ad_vstr(i,j,k,Lold)+                    &
     &                          fac*(ad_vstr(i,j,k,Lnew)-               &
     &                               ad_vstr(i,j,k,Lold))
# ifdef MASKING
            ad_vstr(i,j,k,Lnew)=ad_vstr(i,j,k,Lnew)*vmask(i,j)
# endif
            ad_vstr(i,j,k,Lold)=ad_vstr(i,j,k,Lnew)
          END DO
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
        DO k=1,Nfrec(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              ad_tflux(i,j,k,Lnew,itrc)=ad_tflux(i,j,k,Lold,itrc)+      &
     &                                  fac*(ad_tflux(i,j,k,Lnew,itrc)- &
     &                                       ad_tflux(i,j,k,Lold,itrc))
#  ifdef MASKING
              ad_tflux(i,j,k,Lnew,itrc)=ad_tflux(i,j,k,Lnew,itrc)*      &
     &                                  rmask(i,j)
#  endif
              ad_tflux(i,j,k,Lold,itrc)=ad_tflux(i,j,k,Lnew,itrc)
            END DO
          END DO
        END DO
      END DO
# endif
#endif

      RETURN
      END SUBROUTINE ad_new_state
!
!***********************************************************************
      SUBROUTINE orthogonalize (ng, tile, model,                        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          Lold, Lnew, Lwrk,                       &
     &                          innLoop, outLoop,                       &
#ifdef MASKING
     &                          rmask, umask, vmask,                    &
#endif
#ifdef ADJUST_WSTRESS
     &                          nl_ustr, nl_vstr,                       &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                          nl_tflux,                               &
# endif
     &                          nl_t, nl_u, nl_v,                       &
#else
     &                          nl_ubar, nl_vbar,                       &
#endif
     &                          nl_zeta,                                &
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
     &                          ad_tflux,                               &
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
      USE mod_ncparam
      USE mod_scalars
!
      USE state_addition_mod, ONLY : state_addition
      USE state_dotprod_mod, ONLY : state_dotprod
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Lold, Lnew, Lwrk
      integer, intent(in) :: innLoop, outLoop
!
#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
# endif
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: ad_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_vstr(LBi:,LBj:,:,:)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: ad_tflux(LBi:,LBj:,:,:,:)
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
      real(r8), intent(inout) :: nl_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: nl_vstr(LBi:,LBj:,:,:)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: nl_tflux(LBi:,LBj:,:,:,:)
#  endif
      real(r8), intent(inout) :: nl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: nl_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: nl_v(LBi:,LBj:,:,:)
# else
      real(r8), intent(inout) :: nl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: nl_vbar(LBi:,LBj:,:)
# endif
      real(r8), intent(inout) :: nl_zeta(LBi:,LBj:,:)
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: tl_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_vstr(LBi:,LBj:,:,:)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: tl_tflux(LBi:,LBj:,:,:,:)
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
      real(r8), intent(inout) :: ad_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: ad_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: ad_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
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
      real(r8), intent(inout) :: nl_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: nl_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: nl_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
#  endif
      real(r8), intent(inout) :: nl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: nl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: nl_v(LBi:UBi,LBj:UBj,N(ng),2)
# else
      real(r8), intent(inout) :: nl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: nl_vbar(LBi:UBi,LBj:UBj,3)
# endif
      real(r8), intent(inout) :: nl_zeta(LBi:UBi,LBj:UBj,3)
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: tl_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: tl_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: tl_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
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
      integer :: i, j, k, lstr, rec, Lscale
#ifdef SOLVE3D
      integer :: itrc
#endif
      integer :: L1 = 1
      integer :: L2 = 2

      real(r8) :: fac, fac1, fac2

      real(r8), dimension(0:NstateVar(ng)) :: dot
      real(r8), dimension(0:Ninner) :: DotProd, dot_new, dot_old

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
      DO rec=innLoop,1,-1
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
!  compute its associated dot against current G(k+1). Each gradient
!  solution is loaded into TLM (index Lwrk, if not preconditioning)
!  state arrays.
!
        CALL read_state (ng, tile, model,                               &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   Lwrk, rec,                                     &
     &                   ndefADJ(ng), ncADJid(ng), ncname,              &
#ifdef MASKING
     &                   rmask, umask, vmask,                           &
#endif
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
        IF (exit_flag.ne.NoError) RETURN
!
!  Compute <G(k+1), G(rec)>. Recall that the TLM (index Lwrk) contains
!  G(rec).
!
        CALL state_dotprod (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NstateVar(ng), dot(0:),                     &
#ifdef MASKING
     &                      rmask, umask, vmask,                        &
#endif
#ifdef ADJUST_WSTRESS
     &                      ad_ustr(:,:,:,Lnew), tl_ustr(:,:,:,Lwrk),   &
     &                      ad_vstr(:,:,:,Lnew), tl_vstr(:,:,:,Lwrk),   &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                      ad_tflux(:,:,:,Lnew,:),                     &
     &                      tl_tflux(:,:,:,Lwrk,:),                     &
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
        CALL state_dotprod (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NstateVar(ng), dot(0:),                     &
#ifdef MASKING
     &                      rmask, umask, vmask,                        &
#endif
#ifdef ADJUST_WSTRESS
     &                      tl_ustr(:,:,:,Lwrk), tl_ustr(:,:,:,Lwrk),   &
     &                      tl_vstr(:,:,:,Lwrk), tl_vstr(:,:,:,Lwrk),   &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                      tl_tflux(:,:,:,Lwrk,:),                     &
     &                      tl_tflux(:,:,:,Lwrk,:),                     &
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

        fac1=1.0_r8
        fac2=-DotProd(rec)
!
!  Perform Gramm-Schmidt orthonormalization as:
!
!    ad_var(Lnew) = fac1 * ad_var(Lnew) + fac2 * tl_var(Lwrk)
!
        CALL state_addition (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       Lnew, Lwrk, Lnew, fac1, fac2,              &
#ifdef MASKING
     &                       rmask, umask, vmask,                       &
#endif
#ifdef ADJUST_WSTRESS
     &                       ad_ustr, tl_ustr,                          &
     &                       ad_vstr, tl_vstr,                          &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                       ad_tflux, tl_tflux,                        &
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
      DO rec=innLoop,1,-1
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
        CALL read_state (ng, tile, model,                               &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   Lwrk, rec,                                     &
     &                   ndefADJ(ng), ncADJid(ng), ncname,              &
#ifdef MASKING
     &                   rmask, umask, vmask,                           &
#endif
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
        IF (exit_flag.ne.NoError) RETURN
!
        CALL state_dotprod (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NstateVar(ng), dot(0:),                     &
#ifdef MASKING
     &                      rmask, umask, vmask,                        &
#endif
#ifdef ADJUST_WSTRESS
     &                      ad_ustr(:,:,:,Lnew), tl_ustr(:,:,:,Lwrk),   &
     &                      ad_vstr(:,:,:,Lnew), tl_vstr(:,:,:,Lwrk),   &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                      ad_tflux(:,:,:,Lnew,:),                     &
     &                      tl_tflux(:,:,:,Lwrk,:),                     &
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
        DO rec=innLoop,1,-1
          WRITE (stdout,30) DotProd(rec), rec-1
        END DO
        WRITE (stdout,*) ' '
        DO rec=innLoop,1,-1
          WRITE (stdout,40) innLoop, rec-1, dot_new(rec),               &
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
      SUBROUTINE new_direction (ng, tile, model,                        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          Lwrk, Lnew, betaK,                      &
#ifdef MASKING
     &                          rmask, umask, vmask,                    &
#endif
#ifdef ADJUST_WSTRESS
     &                          ad_ustr, ad_vstr,                       &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                          ad_tflux,                               &
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
#if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      USE mod_scalars
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Lwrk, Lnew

      real(r8), intent(in) :: betaK
!
#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
# endif
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: ad_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_vstr(LBi:,LBj:,:,:)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: ad_tflux(LBi:,LBj:,:,:,:)
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
      real(r8), intent(inout) :: d_sustr(LBi:,LBj:,:)
      real(r8), intent(inout) :: d_svstr(LBi:,LBj:,:)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: d_stflx(LBi:,LBj:,:,:)
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
      real(r8), intent(inout) :: ad_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: ad_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: ad_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
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
      real(r8), intent(inout) :: d_sustr(LBi:UBi,LBj:UBj,Nfrec(ng))
      real(r8), intent(inout) :: d_svstr(LBi:UBi,LBj:UBj,Nfrec(ng))
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: d_stflx(LBi:UBi,LBj:UBj,               &
     &                                   Nfrec(ng),NT(ng))
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
      integer :: i, j, k
#ifdef SOLVE3D
      integer :: itrc
#endif

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute new conjugate descent direction, d(k+1). Notice that the old
!  descent direction is overwritten. Also the initial value is just
!  d(0)=-G(0) since betaK=0 when inner=0.
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
      DO k=1,Nfrec(ng)
        DO j=JstrR,JendR
          DO i=Istr,IendR
            d_sustr(i,j,k)=-ad_ustr(i,j,k,Lnew)+betaK*d_sustr(i,j,k)
# ifdef MASKING
            d_sustr(i,j,k)=d_sustr(i,j,k)*umask(i,j)
# endif
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            d_svstr(i,j,k)=-ad_vstr(i,j,k,Lnew)+betaK*d_svstr(i,j,k)
# ifdef MASKING
            d_svstr(i,j,k)=d_svstr(i,j,k)*vmask(i,j)
# endif
          END DO
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
        DO k=1,Nfrec(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              d_stflx(i,j,k,itrc)=-ad_tflux(i,j,k,Lnew,itrc)+           &
     &                            betaK*d_stflx(i,j,k,itrc)
#  ifdef MASKING
              d_stflx(i,j,k,itrc)=d_stflx(i,j,k,itrc)*rmask(i,j)
#  endif
            END DO
          END DO
        END DO
      END DO
# endif
#endif

      RETURN
      END SUBROUTINE new_direction

!
!***********************************************************************
      SUBROUTINE precond (ng, tile, model,                              &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    innLoop, outLoop,                             &
     &                    NstateVars, Lscale,                           &
     &                    Lritz, Ltrans,                                &
     &                    nConvRitz, Ritz,                              &
#ifdef MASKING
     &                    rmask, umask, vmask,                          &
#endif
#ifdef ADJUST_WSTRESS
     &                    nl_ustr, nl_vstr,                             &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                    nl_tflux,                                     &
# endif
     &                    nl_t, nl_u, nl_v,                             &
#else
     &                    nl_ubar, nl_vbar,                             &
#endif
     &                    nl_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_ncparam
      USE mod_netcdf
      USE mod_iounits
      USE mod_scalars
!
#ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_reduce
#endif
      USE state_addition_mod, ONLY : state_addition
      USE state_copy_mod, ONLY : state_copy
      USE state_dotprod_mod, ONLY : state_dotprod
!
!  Imported variable declarations.
!
      logical, intent(in) :: Lritz
      logical, intent(in) :: Ltrans

      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: innLoop, outLoop
      integer, intent(in) :: NstateVars, Lscale
      integer, intent(in) :: nConvRitz
!
      real(r8), intent(in) :: Ritz(:)
!
#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
# endif
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: nl_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: nl_vstr(LBi:,LBj:,:,:)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: nl_tflux(LBi:,LBj:,:,:,:)
#  endif
      real(r8), intent(inout) :: nl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: nl_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: nl_v(LBi:,LBj:,:,:)
# else
      real(r8), intent(inout) :: nl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: nl_vbar(LBi:,LBj:,:)
# endif
      real(r8), intent(inout) :: nl_zeta(LBi:,LBj:,:)

#else

# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
# endif
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: nl_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: nl_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: nl_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
#  endif
      real(r8), intent(inout) :: nl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: nl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: nl_v(LBi:UBi,LBj:UBj,N(ng),2)
# else
      real(r8), intent(inout) :: nl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: nl_vbar(LBi:UBi,LBj:UBj,3)
# endif
      real(r8), intent(inout) :: nl_zeta(LBi:UBi,LBj:UBj,3)
#endif
!
!  Local variable declarations.
!
      integer :: NSUB, i, j, k, L1, L2, nvec, ndefLCZ, rec
      integer :: lstr
#ifdef SOLVE3D
      integer :: itrc
#endif
      real(r8) :: cff, fac, fac1, fac2, facritz
      real(r8), dimension(0:NstateVars) :: Dotprod

      character (len=80) :: ncname

#ifdef DISTRIBUTE
      character (len=3) :: op_handle
#endif

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  THIS PRECONDITIONER IS WRITTEN IN PRODUCT FORM AS DESCRIBED BY
!  TSHIMANGA - PhD thesis, page 75, proof of proposition 2.3.1.
!  IT IS THEREFORE IMPORTANT THAT THE EIGENVECTORS/RITZ VECTORS THAT
!  ARE COMPUTED BY IS4DVAR_LANCZOS ARE ORTHONORMALIZED.
!
!  Apply the preconditioner. The approximated Hessian matrix is computed
!  from the eigenvectors computed by the Lanczos algorithm which are
!  stored in HSSname NetCDF file.
!-----------------------------------------------------------------------
!
      L1=1
      L2=2
      fac2=0.0_r8
!
!  If using the Ritz preconditioner, read information from the 
!  Lanczos vector file.
!
      IF (Lritz) THEN

        IF (.not.Ltrans) THEN
!
!  Determine if single or multiple Lanczos vector NetCDF files.
!
          CALL netcdf_get_ivar (ng, iADM, TRIM(LCZname(ng)), 'ndefADJ', &
     &                          ndefLCZ)
          IF (exit_flag.ne.NoError) RETURN
!
!  Determine Lanczos vector file to read.
!
          IF (ndefLCZ.gt.0) THEN
            lstr=LEN_TRIM(LCZname(ng))
            WRITE (ncname,10) LCZname(ng)(1:lstr-8), inner
 10         FORMAT (a,'_',i4.4,'.nc')
          ELSE
            ncname=LCZname(ng)
          END IF
!
!  Read in the Lanczos vector q_k+1 computed from the IS4DVAR algorithm
!  first outer loop, where k=nConvRitz. Load Lanczos vectors into NL
!  state arrays at index L2.
!
          rec=nConvRitz+1
          CALL read_state (ng, tile, model,                             &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     L2, rec,                                     &
     &                     ndefLCZ, ncLCZid(ng), TRIM(ncname),          &
# ifdef MASKING
     &                     rmask, umask, vmask,                         &
# endif
# ifdef ADJUST_WSTRESS
     &                     nl_ustr, nl_vstr,                            &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                     nl_tflux,                                    &
#  endif
     &                     nl_t, nl_u, nl_v,                            &
# else
     &                     nl_ubar, nl_vbar,                            &
# endif
     &                     nl_zeta)
          IF (exit_flag.ne.NoError) RETURN
!
!  Compute the dot-product between the input vector and the nConvRitz+1
!  Lanczos vector.
!
          CALL state_dotprod (ng, tile, model,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NstateVars, Dotprod(0:),                  &
#ifdef MASKING
     &                        rmask, umask, vmask,                      &
#endif
#ifdef ADJUST_WSTRESS
     &                        nl_ustr(:,:,:,L1), nl_ustr(:,:,:,L2),     &
     &                        nl_vstr(:,:,:,L1), nl_vstr(:,:,:,L2),     &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                        nl_tflux(:,:,:,L1,:),                     &
     &                        nl_tflux(:,:,:,L2,:),                     &
# endif
     &                        nl_t(:,:,:,L1,:), nl_t(:,:,:,L2,:),       &
     &                        nl_u(:,:,:,L1), nl_u(:,:,:,L2),           &
     &                        nl_v(:,:,:,L1), nl_v(:,:,:,L2),           &
#else
     &                        nl_ubar(:,:,L1), nl_ubar(:,:,L2),         &
     &                        nl_vbar(:,:,L1), nl_vbar(:,:,L2),         &
#endif
     &                        nl_zeta(:,:,L1), nl_zeta(:,:,L2))

        END IF
!
!  Now read the primitive Ritz vectors cg_v and cg_beta.
!
        CALL netcdf_get_fvar (ng, iADM, LCZname(ng), 'cg_beta',         &
     &                        cg_beta)
        IF (exit_flag.ne. NoError) RETURN
        CALL netcdf_get_fvar (ng, iADM, LCZname(ng), 'cg_zv',           &
     &                        cg_zv)
        IF (exit_flag.ne. NoError) RETURN 
!
        IF (Ltrans) THEN
          facritz=cg_beta(nConvRitz,outLoop)
        ELSE
          facritz=cg_beta(nConvRitz,outLoop)*Dotprod(0)
        END IF
      END IF
!
!  Read the converged Hessian eigenvectors into NLM state array,
!  index L2.
!
      DO nvec=1,nConvRitz
        CALL read_state (ng, tile, model,                               &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   L2, nvec,                                      &
     &                   0, ncHSSid(ng), HSSname(ng),                   &
#ifdef MASKING
     &                   rmask, umask, vmask,                           &
#endif
#ifdef ADJUST_WSTRESS
     &                   nl_ustr, nl_vstr,                              &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                   nl_tflux,                                      &
# endif
     &                   nl_t, nl_u, nl_v,                              &
#else
     &                   nl_ubar, nl_vbar,                              &
#endif
     &                   nl_zeta)
        IF (exit_flag.ne.NoError) RETURN
!
!  Compute dot product between input vector and Hessian eigenvector.
!  The input vector is in nl_var(L1) and the Hessian vector in 
!  nl_var(L2)
!
        CALL state_dotprod (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NstateVars, Dotprod(0:),                    &
#ifdef MASKING
     &                      rmask, umask, vmask,                        &
#endif
#ifdef ADJUST_WSTRESS
     &                      nl_ustr(:,:,:,L1), nl_ustr(:,:,:,L2),       &
     &                      nl_vstr(:,:,:,L1), nl_vstr(:,:,:,L2),       &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                      nl_tflux(:,:,:,L1,:),                       &
     &                      nl_tflux(:,:,:,L2,:),                       &
# endif
     &                      nl_t(:,:,:,L1,:), nl_t(:,:,:,L2,:),         &
     &                      nl_u(:,:,:,L1), nl_u(:,:,:,L2),             &
     &                      nl_v(:,:,:,L1), nl_v(:,:,:,L2),             &
#else
     &                      nl_ubar(:,:,L1), nl_ubar(:,:,L2),           &
     &                      nl_vbar(:,:,L1), nl_vbar(:,:,L2),           &
#endif
     &                      nl_zeta(:,:,L1), nl_zeta(:,:,L2))
!
!    Lscale determines the form of the preconditioner:
!
!       1= Spectral LMP
!      -1= Inverse Spectral LMP
!       2= Square root spectral LMP
!      -2= Inverse square root spectral LMP
!
!    tl_var(Lwrk) = fac1 * tl_var(Lwrk) + fac2 * nl_var(L1)
!
        fac1=1.0_r8

        IF (Lscale.eq.-1) THEN
          fac2=(Ritz(nvec)-1.0_r8)*Dotprod(0)
        ELSE IF (Lscale.eq.1) THEN
          fac2=(1.0_r8/Ritz(nvec)-1.0_r8)*Dotprod(0)
        ELSE IF (Lscale.eq.-2) THEN
          fac2=(SQRT(Ritz(nvec))-1.0_r8)*Dotprod(0)
        ELSE IF (Lscale.eq.2) THEN
          fac2=(1.0_r8/SQRT(Ritz(nvec))-1.0_r8)*Dotprod(0)
        END IF

        IF (Lritz.and.Lscale.eq.-2) THEN
         fac2=-fac2
        END IF

        IF(.not.Ltrans) THEN
          IF (Lritz.and.Lscale.eq.-2) THEN
            fac2=fac2+SQRT(Ritz(nvec))*cg_zv(nConvRitz,nvec)*facritz
          END IF
          IF (Lritz.and.Lscale.eq.2) THEN
            fac2=fac2-cg_zv(nConvRitz,nvec)*facritz/Ritz(nvec)
          END IF
        END IF

        CALL state_addition (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       L1, L2, L1, fac1, fac2,                    &
#ifdef MASKING
     &                       rmask, umask, vmask,                       &
#endif
#ifdef ADJUST_WSTRESS
     &                       nl_ustr, nl_ustr,                          &
     &                       nl_vstr, nl_vstr,                          &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                       nl_tflux, nl_tflux,                        &
# endif
     &                       nl_t, nl_t,                                &
     &                       nl_u, nl_u,                                &
     &                       nl_v, nl_v,                                &
#else
     &                       nl_ubar, nl_ubar,                          &
     &                       nl_vbar, nl_vbar,                          &
#endif
     &                       nl_zeta, nl_zeta)
!
        IF (Lritz.and.Ltrans) THEN
!
!  Determine if single or multiple Lanczos vector NetCDF files.
!
          CALL netcdf_get_ivar (ng, iADM, TRIM(LCZname(ng)), 'ndefADJ', &
     &                          ndefLCZ)
          IF (exit_flag.ne.NoError) RETURN
!
!  Determine Lanczos vector file to read.
!
          IF (ndefLCZ.gt.0) THEN
            lstr=LEN_TRIM(LCZname(ng))
            WRITE (ncname,10) LCZname(ng)(1:lstr-8), inner
          ELSE
            ncname=LCZname(ng)
          END IF
!
!  Read in the Lanczos vector q_k+1 computed from the incremental 4DVar
!  algorithm first outer loop, where k=nConvRitz. Load Lanczos vectors
!  into NLM state arrays at index L2.
!
          rec=nConvRitz+1
          CALL read_state (ng, tile, model,                             &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     L2, rec,                                     &
     &                     ndefLCZ, ncLCZid(ng), TRIM(ncname),          &
# ifdef MASKING
     &                     rmask, umask, vmask,                         &
# endif
# ifdef ADJUST_WSTRESS
     &                     nl_ustr, nl_vstr,                            &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                     nl_tflux,                                    &
#  endif
     &                     nl_t, nl_u, nl_v,                            &
# else
     &                     nl_ubar, nl_vbar,                            &
# endif
     &                     nl_zeta)
          IF (exit_flag.ne.NoError) RETURN
!
          IF (Lscale.eq.2) THEN
            fac2=-cg_zv(nConvRitz,nvec)*facritz*Dotprod(0)/Ritz(nvec)
          END IF
          IF (Lscale.eq.-2) THEN
            fac2=SQRT(Ritz(nvec))*cg_zv(nConvRitz,nvec)*facritz*        &
     &                Dotprod(0)
          END IF
!
          CALL state_addition (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         L1, L2, L1, fac1, fac2,                  &
#ifdef MASKING
     &                         rmask, umask, vmask,                     &
#endif
#ifdef ADJUST_WSTRESS
     &                         nl_ustr, nl_ustr,                        &
     &                         nl_vstr, nl_vstr,                        &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                         nl_tflux, nl_tflux,                      &
# endif
     &                         nl_t, nl_t,                              &
     &                         nl_u, nl_u,                              &
     &                         nl_v, nl_v,                              &
#else
     &                         nl_ubar, nl_ubar,                        &
     &                         nl_vbar, nl_vbar,                        &
#endif
     &                         nl_zeta, nl_zeta)

        END IF
      END DO

      RETURN
      END SUBROUTINE precond

!
!***********************************************************************
      SUBROUTINE new_cost (ng, tile, model,                             &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     Lwrk,                                        &
#ifdef MASKING
     &                     rmask, umask, vmask,                         &
#endif
#ifdef ADJUST_WSTRESS
     &                     nl_ustr, nl_vstr,                            &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                     nl_tflux,                                    &
# endif
     &                     nl_t, nl_u, nl_v,                            &
#else
     &                     nl_ubar, nl_vbar,                            &
#endif
     &                     nl_zeta,                                     &
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
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
!
      USE state_addition_mod, ONLY : state_addition
      USE state_dotprod_mod, ONLY : state_dotprod
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Lwrk
!
#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
# endif
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: nl_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: nl_vstr(LBi:,LBj:,:,:)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: nl_tflux(LBi:,LBj:,:,:,:)
#  endif
      real(r8), intent(inout) :: nl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: nl_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: nl_v(LBi:,LBj:,:,:)
# else
      real(r8), intent(inout) :: nl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: nl_vbar(LBi:,LBj:,:)
# endif
      real(r8), intent(inout) :: nl_zeta(LBi:,LBj:,:)
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: tl_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_vstr(LBi:,LBj:,:,:)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: tl_tflux(LBi:,LBj:,:,:,:)
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
      real(r8), intent(inout) :: nl_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: nl_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: nl_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
#  endif
      real(r8), intent(inout) :: nl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: nl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: nl_v(LBi:UBi,LBj:UBj,N(ng),2)
# else
      real(r8), intent(inout) :: nl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: nl_vbar(LBi:UBi,LBj:UBj,3)
# endif
      real(r8), intent(inout) :: nl_zeta(LBi:UBi,LBj:UBj,3)
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: tl_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: tl_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: tl_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
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
      integer :: i, j, k, lstr, rec, Lscale
#ifdef SOLVE3D
      integer :: itrc
#endif
      integer :: L1 = 1
      integer :: L2 = 2

      real(r8) :: fac, fac1, fac2

      real(r8), dimension(0:NstateVar(ng)) :: dot

      character (len=80) :: ncname

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute the cost function based on the formula of Tshimanga
!  (PhD thesis, p 154, eqn A.15):
!
!    J=J_initial+0.5*r'y where J_initial is the value of the
!    cost function on the first iteration (i.e. CostNorm),
!    r is the initial cost function gradient, and y is the 
!    new increment. Note that even r and x are in y-space,
!    their dot-product is equal to that of the same variables
!    transformed to v-space.
!-----------------------------------------------------------------------
!
!  Determine adjoint file to process.
!
      rec=1
      IF (ndefADJ(ng).gt.0) THEN
        lstr=LEN_TRIM(ADJbase(ng))
        WRITE (ncname,10) ADJbase(ng)(1:lstr-3), rec
 10     FORMAT (a,'_',i3.3,'.nc')
      ELSE
        ncname=ADJname(ng)
      END IF
!
!  Read the initial gradient into NLM index L2.
!
      CALL read_state (ng, tile, model,                                 &
     &                 LBi, UBi, LBj, UBj,                              &
     &                 L2, rec,                                         &
     &                 ndefADJ(ng), ncADJid(ng), ncname,                &
#ifdef MASKING
     &                 rmask, umask, vmask,                             &
#endif
#ifdef ADJUST_WSTRESS
     &                 nl_ustr, nl_vstr,                                &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                 nl_tflux,                                        &
# endif
     &                 nl_t, nl_u, nl_v,                                &
#else
     &                 nl_ubar, nl_vbar,                                &
#endif
     &                 nl_zeta)
      IF (exit_flag.ne.NoError) RETURN
!
!  Compute the dot-product of the initial gradient with the current
!  increment.
!
      CALL state_dotprod (ng, tile, model,                              &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NstateVar(ng), dot(0:),                       &
#ifdef MASKING
     &                    rmask, umask, vmask,                          &
#endif
#ifdef ADJUST_WSTRESS
     &                    nl_ustr(:,:,:,L2), tl_ustr(:,:,:,Lwrk),       &
     &                    nl_vstr(:,:,:,L2), tl_vstr(:,:,:,Lwrk),       &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                    nl_tflux(:,:,:,L2,:),                         &
     &                    tl_tflux(:,:,:,Lwrk,:),                       &
# endif
     &                    nl_t(:,:,:,L2,:), tl_t(:,:,:,Lwrk,:),         &
     &                    nl_u(:,:,:,L2), tl_u(:,:,:,Lwrk),             &
     &                    nl_v(:,:,:,L2), tl_v(:,:,:,Lwrk),             &
#else
     &                    nl_ubar(:,:,L2), tl_ubar(:,:,Lwrk),           &
     &                    nl_vbar(:,:,L2), tl_vbar(:,:,Lwrk),           &
#endif
     &                    nl_zeta(:,:,L2), tl_zeta(:,:,Lwrk))
!
!   Compute the new cost function. Only the total value is meaningful.
!
      FOURDVAR(ng)%CostFun(0)=FOURDVAR(ng)%CostNorm(0)+0.5_r8*dot(0)
      DO i=1,NstateVar(ng)
        FOURDVAR(ng)%CostFun(i)=0.0_r8
      END DO

      RETURN
      END SUBROUTINE new_cost

!
!***********************************************************************
      SUBROUTINE read_state (ng, tile, model,                           &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       Lwrk, rec,                                 &
     &                       ndef, ncfileid, ncname,                    &
#ifdef MASKING
     &                       rmask, umask, vmask,                       &
#endif
#ifdef ADJUST_WSTRESS
     &                       s_ustr, s_vstr,                            &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                       s_tflux,                                   &
# endif
     &                       s_t, s_u, s_v,                             &
#else
     &                       s_ubar, s_vbar,                            &
#endif
     &                       s_zeta)
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
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Lwrk, rec, ndef, ncfileid

      character (len=*), intent(in) :: ncname
!
#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
# endif
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: s_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: s_vstr(LBi:,LBj:,:,:)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: s_tflux(LBi:,LBj:,:,:,:)
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
      real(r8), intent(inout) :: s_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: s_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: s_tflux(LBi:UBi,LBj:UBj,               &
     &                                   Nfrec(ng),2,NT(ng))
#  endif
      real(r8), intent(inout) :: s_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: s_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: s_v(LBi:UBi,LBj:UBj,N(ng),2)
# else
      real(r8), intent(inout) :: s_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: s_vbar(LBi:UBi,LBj:UBj,3)
# endif
      real(r8), intent(inout) :: s_zeta(LBi:UBi,LBj:UBj,3)
#endif
!
!  Local variable declarations.
!
      integer :: i, j
#ifdef SOLVE3D
      integer :: itrc, k
#endif
      integer :: gtype, ncid, status
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
!  Read in requested model state record. Load data into state array
!  index Lwrk.
!-----------------------------------------------------------------------
!
!  Determine file and variables ids.
!
      IF (ndef.gt.0) THEN
        IF (InpThread) THEN
          status=nf90_open(TRIM(ncname), nf90_nowrite, ncid)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) TRIM(ncname)
            exit_flag=2
            ioerror=status
            RETURN
          END IF
        END IF
      ELSE
        ncid=ncfileid
      END IF
      IF (InpThread) THEN
#ifndef SOLVE3D
        status=nf90_inq_varid(ncid, TRIM(Vname(1,idUbar)), vid(idUbar))
        status=nf90_inq_varid(ncid, TRIM(Vname(1,idVbar)), vid(idVbar))
#endif
        status=nf90_inq_varid(ncid, TRIM(Vname(1,idFsur)), vid(idFsur))
#ifdef ADJUST_WSTRESS
        status=nf90_inq_varid(ncid, TRIM(Vname(1,idUsms)), vid(idUsms))
        status=nf90_inq_varid(ncid, TRIM(Vname(1,idVsms)), vid(idVsms))
#endif
#ifdef SOLVE3D
        status=nf90_inq_varid(ncid, TRIM(Vname(1,idUvel)), vid(idUvel))
        status=nf90_inq_varid(ncid, TRIM(Vname(1,idVvel)), vid(idVvel))
        DO itrc=1,NT(ng)
          status=nf90_inq_varid(ncid, TRIM(Vname(1,idTvar(itrc))),      &
     &                          vid(idTvar(itrc)))
# ifdef ADJUST_STFLUX
          status=nf90_inq_varid(ncid, TRIM(Vname(1,idTsur(itrc))),      &
     &                          vid(idTsur(itrc)))
# endif
        END DO
#endif
      END IF
      DO i=1,4
        Vsize(i)=0
      END DO
!
!  Read in free-surface.
!
      gtype=r2dvar
      scale=1.0_r8
      status=nf_fread2d(ng, iTLM, ncid, vid(idFsur), rec, gtype,        &
     &                  Vsize, LBi, UBi, LBj, UBj,                      &
     &                  scale, Fmin, Fmax,                              &
#ifdef MASKING
     &                  rmask(LBi,LBj),                                 &
#endif
     &                  s_zeta(LBi,LBj,Lwrk))
      IF (status.ne.nf90_noerr) THEN
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
      gtype=u2dvar
      scale=1.0_r8
      status=nf_fread2d(ng, iTLM, ncid, vid(idUbar), rec, gtype,        &
     &                  Vsize, LBi, UBi, LBj, UBj,                      &
     &                  scale, Fmin, Fmax,                              &
# ifdef MASKING
     &                  umask(LBi,LBj),                                 &
# endif
     &                  s_ubar(LBi,LBj,Lwrk))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idUbar)), rec, TRIM(ncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF

      gtype=v2dvar
      scale=1.0_r8
      status=nf_fread2d(ng, iTLM, ncid, vid(idVbar), rec, gtype,        &
     &                  Vsize, LBi, UBi, LBj, UBj,                      &
     &                  scale, Fmin, Fmax,                              &
# ifdef MASKING
     &                  vmask(LBi,LBj),                                 &
# endif
     &                  s_vbar(LBi,LBj,Lwrk))
      IF (status.ne.nf90_noerr) THEN
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
!  Read surface momentum stress.
!
      gtype=u3dvar
      scale=1.0_r8/rho0                           ! N/m2 (Pa) to m2/s2
      status=nf_fread3d(ng, iTLM, ncid, vid(idUsms), rec, gtype,        &
     &                  Vsize, LBi, UBi, LBj, UBj, 1, Nfrec(ng),        &
     &                  scale, Fmin, Fmax,                              &
# ifdef MASKING
     &                  umask(LBi,LBj),                                 &
# endif
     &                  s_ustr(LBi,LBj,1,Lwrk))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idUsms)), rec, TRIM(ncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF

      gtype=v3dvar
      scale=1.0_r8/rho0                           ! N/m2 (Pa) to m2/s2
      status=nf_fread3d(ng, iTLM, ncid, vid(idVsms), rec, gtype,        &
     &                  Vsize, LBi, UBi, LBj, UBj, 1, Nfrec(ng),        &
     &                  scale, Fmin, Fmax,                              &
# ifdef MASKING
     &                  vmask(LBi,LBj),                                 &
# endif
     &                  s_vstr(LBi,LBj,1,Lwrk))
      IF (status.ne.nf90_noerr) THEN
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
      gtype=u3dvar
      scale=1.0_r8
      status=nf_fread3d(ng, iTLM, ncid, vid(idUvel), rec, gtype,        &
     &                  Vsize, LBi, UBi, LBj, UBj, 1, N(ng),            &
     &                  scale, Fmin, Fmax,                              &
# ifdef MASKING
     &                  umask(LBi,LBj),                                 &
# endif
     &                  s_u(LBi,LBj,1,Lwrk))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idUvel)), rec, TRIM(ncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF

      gtype=v3dvar
      scale=1.0_r8
      status=nf_fread3d(ng, iTLM, ncid, vid(idVvel), rec, gtype,        &
     &                  Vsize, LBi, UBi, LBj, UBj, 1, N(ng),            &
     &                  scale, Fmin, Fmax,                              &
# ifdef MASKING
     &                  vmask(LBi,LBj),                                 &
# endif
     &                  s_v(LBi,LBj,1,Lwrk))
      IF (status.ne.nf90_noerr) THEN
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
      gtype=r3dvar
      scale=1.0_r8
      DO itrc=1,NT(ng)
        status=nf_fread3d(ng, iTLM, ncid, vid(idTvar(itrc)), rec,       &
     &                    gtype, Vsize, LBi, UBi, LBj, UBj, 1, N(ng),   &
     &                    scale, Fmin, Fmax,                            &
# ifdef MASKING
     &                    rmask(LBi,LBj),                               &
# endif
     &                    s_t(LBi,LBj,1,Lwrk,itrc))
        IF (status.ne.nf90_noerr) THEN
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
      gtype=r3dvar
      DO itrc=1,NT(ng)
        IF (itrc.eq.itemp) THEN
          scale=1.0_r8/(rho0*Cp)                  ! W/m2 to Celsius m/s
        ELSE
          scale=1.0_r8
        END IF
        status=nf_fread3d(ng, iTLM, ncid, vid(idTsur(itrc)), rec,       &
     &                    gtype, Vsize, LBi, UBi, LBj, UBj, 1,Nfrec(ng),&
     &                    scale, Fmin, Fmax,                            &
#  ifdef MASKING
     &                    rmask(LBi,LBj),                               &
#  endif
     &                    s_tflux(LBi,LBj,1,Lwrk,itrc))
        IF (status.ne.nf90_noerr) THEN
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
!  If multiple files, close current file.
!
      IF (ndef.gt.0) THEN
        status=nf90_close(ncid)
      END IF
!
 10   FORMAT (' READ_STATE - unable to open NetCDF file: ',a)
 20   FORMAT (' READ_STATE - error while reading variable: ',a,2x,      &
     &        'at time record = ',i3,/,16x,'in NetCDF file: ',a)

      RETURN
      END SUBROUTINE read_state

      SUBROUTINE cg_write (ng, model, innLoop, outLoop)
!
!=======================================================================
!                                                                      !
!  This routine writes conjugate gradient vectors into 4DVAR NetCDF    !
!  for restart purposes.                                               !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      Use mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
# ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcasti
# endif
!
      implicit none
!
!  Imported variable declarations
!
      integer, intent(in) :: ng, model, innLoop, outLoop
!
!  Local variable declarations.
!
      integer :: status
!
!-----------------------------------------------------------------------
!  Write out conjugate gradient vectors.
!-----------------------------------------------------------------------
!
!  Write out outer and inner iteration.
!
      CALL netcdf_put_ivar (ng, model, MODname(ng), 'outer',            &
     &                      outer, (/0/), (/0/),                        &
     &                      ncid = ncMODid(ng))
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_ivar (ng, model, MODname(ng), 'inner',            &
     &                      inner, (/0/), (/0/),                        &
     &                      ncid = ncMODid(ng))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out number of converged Ritz eigenvalues.
!
      IF ((innLoop.eq.0).and.(outloop.eq.1)) THEN
        CALL netcdf_put_ivar (ng, model, MODname(ng), 'nConvRitz',      &
     &                        nConvRitz, (/0/), (/0/),                  &
     &                        ncid = ncMODid(ng))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out converged Ritz eigenvalues.
!
      IF ((innLoop.eq.0).and.(outloop.eq.1)) THEN
        CALL netcdf_put_fvar (ng, model, MODname(ng), 'Ritz',           &
     &                        Ritz, (/1/), (/nConvRitz/),               &
     &                        ncid = ncMODid(ng))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out conjugate gradient norms.
!
      CALL netcdf_put_fvar (ng, model, MODname(ng), 'cg_alpha',         &
     &                      cg_alpha(0:,:), (/1,1/),                    &
     &                      (/Ninner+1,Nouter/),                        &
     &                      ncid = ncMODid(ng))
      IF (exit_flag.ne.NoError) RETURN
!
      CALL netcdf_put_fvar (ng, model, MODname(ng), 'cg_beta',          &
     &                      cg_beta(0:,:), (/1,1/),                     &
     &                      (/Ninner+1,Nouter/),                        &
     &                      ncid = ncMODid(ng))
      IF (exit_flag.ne.NoError) RETURN
!
      CALL netcdf_put_fvar (ng, model, MODname(ng), 'cg_tau',           &
     &                      cg_tau(0:,:), (/1,1/),                      &
     &                      (/Ninner+1,Nouter/),                        &
     &                      ncid = ncMODid(ng))
      IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Synchronize model/observation NetCDF file to disk.
!-----------------------------------------------------------------------
!
      IF (OutThread) THEN
        status=nf90_sync(ncMODid(ng))
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,10) TRIM(MODname(ng))
          exit_flag=3
          ioerror=status
        END IF
      END IF
#ifdef DISTRIBUTE
      CALL mp_bcasti (ng, model, exit_flag, 1)
#endif

  10  FORMAT (/,' CG_WRITE - unable to synchronize to disk ',           &
     &        ' model/observation file: ',/,12x,a)

      RETURN
      END SUBROUTINE cg_write
