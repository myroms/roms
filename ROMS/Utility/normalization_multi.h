#include "cppdefs.h"
      MODULE normalization_mod
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2026 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  Normalization Coefficients Multiscale Approach.                     !
!                                                                      !
!  This module computes normalization factors to model the spreading   !
!  of the background-error symmetric covariance matrix (B). These      !
!  coefficients ensure that the diagonal elements of B are equal to    !
!  unity. In applications involving land or sea masking, the           !
!  normalization factors can significantly modify the covariance       !
!  structure near boundaries.                                          !
!                                                                      !
!  Normalization coefficients may be computed using either the exact   !
!  method or the approximated randomization method:                    !
!                                                                      !
!  The exact method is computationally intensive. In this approach,    !
!  normalization factors are determined by perturbing each model grid  !
!  cell with a delta function scaled by area (for 2D factors) or       !
!  volume (for 3D factors), then convolving with the square-root       !
!  adjoint and tangent diffusion operators.                            !
!                                                                      !
!  The approximated method is less computationally demanding. In this  !
!  approach, normalization factors are computed using the randomization!
!  technique described by Fisher and Courtier (1995). Factors are      !
!  initialized with random numbers drawn from a normal distribution    !
!  with zero mean and unit variance. These factors are then scaled by  !
!  the inverse-square-root of cell area (for 2D factors) or volume     !
!  (for 3D factors), and convolved with the square-root adjoint and    !
!  tangent diffusion operators over a specified number of iterations,  !
!  denoted as Nrandom.                                                 !
!                                                                      !
!  Multiscale Background-error Covariance (B) Modeling/Spreading:      !
!                                                                      !
!  Spreading is modeled using pseudo-diffusion operators in spatial    !
!  correlation space, where diffusion coefficients (K) are proportional!
!  to the square of the correlation length scale (Daley, 1992).        !
!  Correlation scales may be constant or spatially varying for each    !
!  variable in the control vector. The background-error covariance     !
!  matrix (B) can be represented using either a multiscale aproach.    !
!  approach.                                                           !
!                                                                      !
!  In the multiscale formulation, B is expressed as a linear           !
!  combination of distinct spatial scales, ranging from large to small.!
!  This approach enables the representation of both broad structures   !
!  and fine-scale features, while reducing scale aliasing in the data  !
!  assimilation cost function (Weaver et al., 2016). In contrast, the  !
!  monoscale formulation employs a single scale.                       !
!                                                                      !
!  The multiscale approach uses an implicit horizontal pseudo-         !
!  diffusion operator implemented via Conjugate Gradient (CG) and      !
!  Chebyshev Iterations (CI). Conversely, the monoscale default method !
!  employs an explicit horizontal algorithm. Error correlations are    !
!  considered separable in the horizontal and vertical directions.     !
!  The vertical diffusion operator is also implicit.                   !
!                                                                      !
!  The multiscale concept, expressed as B = SUM(Wi*Bi), represents a   !
!  weighted sum of B values corresponding to various scales. Typically,!
!  two to four different scales are combined in practical applications.!
!  The weight coefficients Wi are constrained to sum to unity. The     !
!  resulting correlation functions belong to the Matýrn-class family,  !
!  which accommodates complex functional shapes.                       !
!                                                                      !
!  Horizontally spatially varying correlation length scales for the    !
!  K-diffusion coefficient are permitted. A horizontal map of isotropic!
!  or anisotropic correlations can be specified and is read from an    !
!  input NetCDF file. For computational efficiency, these correlations !
!  require an implicit diffusion CG/CI solver.                         !
!                                                                      !
!  References:                                                         !
!                                                                      !
!    Weaver, A. and P. Courtier, 2001: Correlation modeling on the     !
!      sphere using a generalized diffusion equation, Q.J.R. Meteorol. !
!      Soc, 127, 1815-1846, doi:10.1002/qj.49712757518.                !
!                                                                      !
!    Weaver, A.T., J. Tshimanga, and A. Piacentini, 2016: Correlation  !
!      operators based on an implicitly formulated diffusion equation  !
!      solved with the Chebyshev iteration, Q.J.R. Meteorol. Soc.,     !
!      142, 455-471, doi:10.1002/qj.2664.                              !
!                                                                      !
!=======================================================================
!
      USE mod_kinds
      USE mod_param
      USE mod_parallel
#ifdef ADJUST_BOUNDARY
      USE mod_boundary
#endif
      USE mod_grid
#if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
      USE mod_forces
#endif
      USE mod_fourdvar
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_mixing
      USE mod_ocean
#if defined PIO_LIB && defined DISTRIBUTE
      USE mod_pio_netcdf
#endif
      USE mod_scalars
#if defined SEDIMENT && defined SED_MORPH && defined SOLVE3D
      USE mod_sedbed
#endif
      USE mod_stepping
!
      USE ad_conv_2d_mod
#ifdef SOLVE3D
      USE ad_conv_3d_mod
#endif
#ifdef ADJUST_BOUNDARY
      USE ad_conv_bry2d_mod
# ifdef SOLVE3D
      USE ad_conv_bry3d_mod
# endif
#endif
      USE bc_2d_mod
      USE ad_bc_2d_mod
#ifdef SOLVE3D
      USE bc_3d_mod
      USE ad_bc_3d_mod
#endif
#ifdef ADJUST_BOUNDARY
      USE bc_bry2d_mod
# ifdef SOLVE3D
      USE bc_bry3d_mod
# endif
#endif
#ifdef DISTRIBUTE
      USE mp_exchange_mod
#endif
      USE set_depth_mod
      USE tl_conv_2d_mod
#ifdef SOLVE3D
      USE tl_conv_3d_mod
#endif
#ifdef ADJUST_BOUNDARY
      USE tl_conv_bry2d_mod
# ifdef SOLVE3D
      USE tl_conv_bry3d_mod
# endif
#endif
      USE white_noise_mod
!
#ifdef DISTRIBUTE
# ifdef ADJUST_BOUNDARY
      USE distribute_mod,      ONLY : mp_collect
# endif
      USE distribute_mod,      ONLY : mp_reduce
#endif
      USE roms_multiscale_mod, ONLY : multiscale         ! CLASS object

      USE nf_fwrite2d_mod,     ONLY : nf_fwrite2d
      USE nf_fwrite3d_mod,     ONLY : nf_fwrite3d
      USE strings_mod,         ONLY : FoundError
!
      implicit none
!
      PUBLIC  :: normalization
      PRIVATE :: normalization_tile
      PRIVATE :: randomization_tile
!
      PRIVATE :: dot_prod2d
#ifdef SOLVE3D
      PRIVATE :: dot_prod3d
#endif
!
      PRIVATE :: wrt_norm2d_nf90
#if defined PIO_LIB && defined DISTRIBUTE
      PRIVATE :: wrt_norm2d_pio
#endif
      PRIVATE :: wrt_norm3d_nf90
#if defined PIO_LIB && defined DISTRIBUTE
      PRIVATE :: wrt_norm3d_pio
#endif
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE normalization (ng, tile, ifac)
!***********************************************************************
!
      USE roms_multiscale_mod, ONLY : B_ms
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, ifac
!
!  Local variable declarations.
!
#include "tile.h"
!
!  Compute background error covariance normalization factors using
!  the very expensive exact method.
!
      IF (Nmethod(ng).eq.0) THEN
        CALL normalization_tile (B_ms(ns), ng, tile,                    &
     &                           LBi, UBi, LBj, UBj, 1, N(ng),          &
     &                           LBij, UBij,                            &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           nstp(ng), nnew(ng), ifac,              &
     &                           GRID(ng) % pm,                         &
     &                           GRID(ng) % om_p,                       &
     &                           GRID(ng) % om_r,                       &
     &                           GRID(ng) % om_u,                       &
     &                           GRID(ng) % om_v,                       &
     &                           GRID(ng) % pn,                         &
     &                           GRID(ng) % on_p,                       &
     &                           GRID(ng) % on_r,                       &
     &                           GRID(ng) % on_u,                       &
     &                           GRID(ng) % on_v,                       &
     &                           GRID(ng) % pmon_p,                     &
     &                           GRID(ng) % pmon_r,                     &
     &                           GRID(ng) % pmon_u,                     &
     &                           GRID(ng) % pnom_p,                     &
     &                           GRID(ng) % pnom_r,                     &
     &                           GRID(ng) % pnom_v,                     &
#ifdef MASKING
     &                           GRID(ng) % pmask,                      &
     &                           GRID(ng) % rmask,                      &
     &                           GRID(ng) % umask,                      &
     &                           GRID(ng) % vmask,                      &
#endif
#ifdef SOLVE3D
     &                           GRID(ng) % h,                          &
# ifdef ICESHELF
     &                           GRID(ng) % zice,                       &
# endif
# if defined SEDIMENT && defined SED_MORPH
     &                           SEDBED(ng) % bed_thick,                &
# endif
     &                           GRID(ng) % Hz,                         &
     &                           GRID(ng) % z_r,                        &
     &                           GRID(ng) % z_w,                        &
#endif
     &                           MIXING(ng) % Kh,                       &
#ifdef SOLVE3D
     &                           MIXING(ng) % Kv,                       &
#endif
#ifdef ADJUST_BOUNDARY
# ifdef SOLVE3D
     &                           BOUNDARY(ng) % b_t_obc,                &
     &                           BOUNDARY(ng) % b_u_obc,                &
     &                           BOUNDARY(ng) % b_v_obc,                &
# endif
     &                           BOUNDARY(ng) % b_ubar_obc,             &
     &                           BOUNDARY(ng) % b_vbar_obc,             &
     &                           BOUNDARY(ng) % b_zeta_obc,             &
#endif
#ifdef ADJUST_WSTRESS
     &                           FORCES(ng) % b_sustr,                  &
     &                           FORCES(ng) % b_svstr,                  &
#endif
#if defined ADJUST_STFLUX && defined SOLVE3D
     &                           FORCES(ng) % b_stflx,                  &
#endif
#ifdef SOLVE3D
     &                           OCEAN(ng) % b_t,                       &
     &                           OCEAN(ng) % b_u,                       &
     &                           OCEAN(ng) % b_v,                       &
#endif
     &                           OCEAN(ng) % b_zeta,                    &
     &                           OCEAN(ng) % b_ubar,                    &
     &                           OCEAN(ng) % b_vbar)
!
!  Compute background error covariance normalization factors using
!  the approximated randomization method.
!
      ELSE IF (Nmethod(ng).eq.1) THEN
        CALL randomization_tile (B_ms(ng), ng, tile,                    &
     &                           LBi, UBi, LBj, UBj, 1, N(ng),          &
     &                           LBij, UBij,                            &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           nstp(ng), nnew(ng), ifac,              &
     &                           GRID(ng) % pm,                         &
     &                           GRID(ng) % om_p,                       &
     &                           GRID(ng) % om_r,                       &
     &                           GRID(ng) % om_u,                       &
     &                           GRID(ng) % om_v,                       &
     &                           GRID(ng) % pn,                         &
     &                           GRID(ng) % on_p,                       &
     &                           GRID(ng) % on_r,                       &
     &                           GRID(ng) % on_u,                       &
     &                           GRID(ng) % on_v,                       &
     &                           GRID(ng) % pmon_p,                     &
     &                           GRID(ng) % pmon_r,                     &
     &                           GRID(ng) % pmon_u,                     &
     &                           GRID(ng) % pnom_p,                     &
     &                           GRID(ng) % pnom_r,                     &
     &                           GRID(ng) % pnom_v,                     &
#ifdef MASKING
     &                           GRID(ng) % pmask,                      &
     &                           GRID(ng) % rmask,                      &
     &                           GRID(ng) % umask,                      &
     &                           GRID(ng) % vmask,                      &
#endif
#ifdef SOLVE3D
     &                           GRID(ng) % h,                          &
# ifdef ICESHELF
     &                           GRID(ng) % zice,                       &
# endif
# if defined SEDIMENT && defined SED_MORPH
     &                           SEDBED(ng) % bed_thick,                &
# endif
     &                           GRID(ng) % Hz,                         &
     &                           GRID(ng) % z_r,                        &
     &                           GRID(ng) % z_w,                        &
#endif
     &                           MIXING(ng) % Kh,                       &
#ifdef SOLVE3D
     &                           MIXING(ng) % Kv,                       &
#endif
#ifdef ADJUST_BOUNDARY
# ifdef SOLVE3D
     &                           BOUNDARY(ng) % b_t_obc,                &
     &                           BOUNDARY(ng) % b_u_obc,                &
     &                           BOUNDARY(ng) % b_v_obc,                &
# endif
     &                           BOUNDARY(ng) % b_ubar_obc,             &
     &                           BOUNDARY(ng) % b_vbar_obc,             &
     &                           BOUNDARY(ng) % b_zeta_obc,             &
#endif
#ifdef ADJUST_WSTRESS
     &                           FORCES(ng) % b_sustr,                  &
     &                           FORCES(ng) % b_svstr,                  &
#endif
#if defined ADJUST_STFLUX && defined SOLVE3D
     &                           FORCES(ng) % b_stflx,                  &
#endif
#ifdef SOLVE3D
     &                           OCEAN(ng) % b_t,                       &
     &                           OCEAN(ng) % b_u,                       &
     &                           OCEAN(ng) % b_v,                       &
#endif
     &                           OCEAN(ng) % b_zeta,                    &
     &                           OCEAN(ng) % b_ubar,                    &
     &                           OCEAN(ng) % b_vbar)
      END IF
!
      RETURN
      END SUBROUTINE normalization
!
!***********************************************************************
      SUBROUTINE normalization_tile (self, ng, tile,                    &
     &                               LBi, UBi, LBj, UBj, LBk, UBk,      &
     &                               LBij, UBij,                        &
     &                               IminS, ImaxS, JminS, JmaxS,        &
     &                               nstp, nnew, ifac,                  &
     &                               pm, om_p, om_r, om_u, om_v,        &
     &                               pn, on_p, on_r, on_u, on_v,        &
     &                               pmon_p, pmon_r, pmon_u,            &
     &                               pnom_p, pnom_r, pnom_v,            &
#ifdef MASKING
     &                               pmask, rmask,                      &
     &                               umask, vmask,                      &
#endif
#ifdef SOLVE3D
     &                               h,                                 &
# ifdef ICESHELF
     &                               zice,                              &
# endif
# if defined SEDIMENT && defined SED_MORPH
     &                               bed_thick,                         &
# endif
     &                               Hz, z_r, z_w,                      &
#endif
     &                               Kh,                                &
#ifdef SOLVE3D
     &                               Kv,                                &
#endif
#ifdef ADJUST_BOUNDARY
# ifdef SOLVE3D
     &                               VnormRobc, VnormUobc, VnormVobc,   &
# endif
     &                               HnormRobc, HnormUobc, HnormVobc,   &
#endif
#ifdef ADJUST_WSTRESS
     &                               HnormSUS, HnormSVS,                &
#endif
#if defined ADJUST_STFLUX && defined SOLVE3D
     &                               HnormSTF,                          &
#endif
#ifdef SOLVE3D
     &                               VnormR, VnormU, VnormV,            &
#endif
     &                               HnormR, HnormU, HnormV)
!***********************************************************************
!
!  Imported variable declarations.
!
      CLASS (multiscale), intent(inout) :: self
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk, LBij, UBij
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp, nnew, ifac
!
#ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: om_p(LBi:,LBj:)
      real(r8), intent(in) :: om_r(LBi:,LBj:)
      real(r8), intent(in) :: om_u(LBi:,LBj:)
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: on_p(LBi:,LBj:)
      real(r8), intent(in) :: on_r(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)
      real(r8), intent(in) :: on_v(LBi:,LBj:)
      real(r8), intent(in) :: pmon_p(LBi:,LBj:)
      real(r8), intent(in) :: pmon_r(LBi:,LBj:)
      real(r8), intent(in) :: pmon_u(LBi:,LBj:)
      real(r8), intent(in) :: pnom_p(LBi:,LBj:)
      real(r8), intent(in) :: pnom_r(LBi:,LBj:)
      real(r8), intent(in) :: pnom_v(LBi:,LBj:)
# ifdef MASKING
      real(r8), intent(in) :: pmask(LBi:,LBj:)
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
# endif
      real(r8), intent(in) :: Kh(LBi:,LBj:)
# ifdef SOLVE3D
      real(r8), intent(in) :: Kv(LBi:,LBj:,0:)
#  ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:,LBj:)
#  endif
#  if defined SEDIMENT && defined SED_MORPH
      real(r8), intent(in):: bed_thick(LBi:,LBj:,:)
#  endif
      real(r8), intent(inout) :: h(LBi:,LBj:)
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
      real(r8), intent(out) :: VnormRobc(LBij:,:,:,:)
      real(r8), intent(out) :: VnormUobc(LBij:,:,:)
      real(r8), intent(out) :: VnormVobc(LBij:,:,:)
#  endif
      real(r8), intent(out) :: HnormRobc(LBij:,:)
      real(r8), intent(out) :: HnormUobc(LBij:,:)
      real(r8), intent(out) :: HnormVobc(LBij:,:)
# endif
# ifdef ADJUST_WSTRESS
      real(r8), intent(out) :: HnormSUS(LBi:,LBj:)
      real(r8), intent(out) :: HnormSVS(LBi:,LBj:)
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(out) :: HnormSTF(LBi:,LBj:,:)
# endif
# ifdef SOLVE3D
      real(r8), intent(out) :: VnormR(LBi:,LBj:,:,:,:)
      real(r8), intent(out) :: VnormU(LBi:,LBj:,:,:)
      real(r8), intent(out) :: VnormV(LBi:,LBj:,:,:)
# endif
      real(r8), intent(out) :: HnormR(LBi:,LBj:,:)
      real(r8), intent(out) :: HnormU(LBi:,LBj:,:)
      real(r8), intent(out) :: HnormV(LBi:,LBj:,:)
# ifdef SOLVE3D
      real(r8), intent(out) :: Hz(LBi:,LBj:,:)
      real(r8), intent(out) :: z_r(LBi:,LBj:,:)
      real(r8), intent(out) :: z_w(LBi:,LBj:,0:)
# endif
#else
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_p(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_r(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_p(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_r(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pmon_p(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pmon_r(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pmon_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pnom_p(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pnom_r(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pnom_v(LBi:UBi,LBj:UBj)
# ifdef MASKING
      real(r8), intent(in) :: pmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in) :: Kh(LBi:UBi,LBj:UBj)
# ifdef SOLVE3D
      real(r8), intent(in) :: Kv(LBi:UBi,LBj:UBj,0:N(ng))
#  ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:UBi,LBj:UBj)
#  endif
#  if defined SEDIMENT && defined SED_MORPH
      real(r8), intent(in):: bed_thick0(LBi:UBi,LBj:UBj,3)
#  endif
      real(r8), intent(inout) :: h(LBi:UBi,LBj:UBj)
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
      real(r8), intent(out) :: VnormRobc(LBij:UBij,N(ng),4,NT(ng))
      real(r8), intent(out) :: VnormUobc(LBij:UBij,N(ng),4)
      real(r8), intent(out) :: VnormVobc(LBij:UBij,N(ng),4)
#  endif
      real(r8), intent(out) :: HnormRobc(LBij:UBij,4)
      real(r8), intent(out) :: HnormUobc(LBij:UBij,4)
      real(r8), intent(out) :: HnormVobc(LBij:UBij,4)
# endif
# ifdef ADJUST_WSTRESS
      real(r8), intent(out) :: HnormSUS(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: HnormSVS(LBi:UBi,LBj:UBj)
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(out) :: HnormSTF(LBi:UBi,LBj:UBj,NT(ng))
# endif
# ifdef SOLVE3D
      real(r8), intent(out) :: VnormR(LBi:UBi,LBj:UBj,N(ng),NSA,NT(ng))
      real(r8), intent(out) :: VnormU(LBi:UBi,LBj:UBj,N(ng),NSA)
      real(r8), intent(out) :: VnormV(LBi:UBi,LBj:UBj,N(ng),NSA)
# endif
      real(r8), intent(out) :: HnormR(LBi:UBi,LBj:UBj,NSA)
      real(r8), intent(out) :: HnormU(LBi:UBi,LBj:UBj,NSA)
      real(r8), intent(out) :: HnormV(LBi:UBi,LBj:UBj,NSA)
# ifdef SOLVE3D
      real(r8), intent(out) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(out) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(out) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
# endif
#endif
!
!  Local variable declarations.
!
#ifdef SOLVE3D
      logical :: Ldiffer, Lsame
#endif
#ifdef ADJUST_BOUNDARY
      logical :: bounded
      logical :: Lconvolve(4)
#endif
#ifdef MULTI_SCALE_B
      logical :: Lweak
#endif
!
      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, ic, ifile, is, j, jc, rec
      integer :: iLap
#ifdef MULTI_SCALE_B
      integer :: it, klev, ms, ns
#endif
      integer :: NSUB
#ifdef SOLVE3D
      integer :: UBt, itrc, k, kc, ntrc
#endif
#ifdef ADJUST_BOUNDARY
      integer :: Bmin, Bmax, IJlen, IJKlen
      integer :: ib, ibry, ifield, kb
!
      real(r8), parameter :: Aspv = 0.0_r8
#endif
      real(dp) :: my_time
      real(r8) :: cff, compute
      real(r8) :: my_dot, Gdotp
!
      real(r8), dimension(LBi:UBi,LBj:UBj) :: A2d
      real(r8), dimension(LBi:UBi,LBj:UBj) :: Hscale
#ifdef SOLVE3D
      real(r8), dimension(LBi:UBi,LBj:UBj,1:N(ng)) :: A3d
      real(r8), dimension(LBi:UBi,LBj:UBj,1:N(ng)) :: Vscale
#endif
#ifdef ADJUST_BOUNDARY
      real(r8), dimension(LBij:UBij) :: B2d
      real(r8), dimension(LBij:UBij) :: HscaleB
# ifdef SOLVE3D
      real(r8), dimension(LBij:UBij,1:N(ng)) :: B3d
      real(r8), dimension(LBij:UBij,1:N(ng)) :: VscaleB
#  ifdef DISTRIBUTE
      real(r8), dimension((UBij-LBij+1)*N(ng)) :: Bwrk
#  endif
# endif
#endif
!
#ifdef DISTRIBUTE
      character (len=3  ) :: op_handle
#endif
      character (len=40 ) :: Text
      character (len=256) :: ncname

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", normalization_tile"

#if defined PIO_LIB && defined DISTRIBUTE
!
      TYPE (IO_Desc_t),  pointer :: ioDesc
#endif

#include "set_bounds.h"
!
      SourceFile=MyFile

      my_time=tdays(ng)*day2sec

#ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Compute time invariant depths (use zero free-surface).
!-----------------------------------------------------------------------
!
      DO i=LBi,UBi
        DO j=LBj,UBj
          A2d(i,j)=0.0_r8
        END DO
      END DO

      CALL set_depth_tile (ng, tile, iNLM,                              &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     nstp, nnew,                                  &
     &                     h,                                           &
# ifdef ICESHELF
     &                     zice,                                        &
# endif
# if defined SEDIMENT && defined SED_MORPH
     &                     bed_thick,                                   &
# endif
     &                     A2d,                                         &
     &                     Hz, z_r, z_w)
#endif
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Compute initial conditions and model error covariance, B,
!  normalization factors using the exact method. It involves
!  computing the filter variance (convolution) at each point
!  independenly.  That is, each point is perturbed with a delta
!  function, scaled by the inverse squared root of the area (2D)
!  or volume (3D), and then convoluted.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      IF (Master) WRITE (stdout,10)

      FILE_LOOP : DO ifile=1,NSA

        IF (LwrtNRM(ifile,ng)) THEN
          Lweak=.FALSE.
          IF (ifile.eq.1) THEN
            Text='initial conditions'
          ELSE IF (ifile.eq.2) THEN
            Text='model'
            Lweak=.TRUE.
          END IF
!
!  Set time record index to write in normalization NetCDF file.
!
          ncname=NRM(ifile,ng)%name
          NRM(ifile,ng)%Rindex=NRM(ifile,ng)%Rindex+1
          NRM(ifile,ng)%Nrec=NRM(ifile,ng)%Nrec+1
!
!  Write out model time (s).
!
          SELECT CASE (NRM(ifile,ng)%IOtype)
            CASE (io_nf90)
              CALL netcdf_put_fvar (ng, iTLM, ncname,                   &
     &                              Vname(1,idtime), my_time,           &
     &                         start = (/NRM(ifile,ng)%Rindex/),        &
     &                         total = (/1/),                           &
     &                         ncid = NRM(ifile,ng)%ncid,               &
     &                         varid = NRM(ifile,ng)%Vid(idtime))

#if defined PIO_LIB && defined DISTRIBUTE
            CASE (io_pio)
              CALL pio_netcdf_put_fvar (ng, iTLM, ncname,               &
     &                                  Vname(1,idtime), my_time,       &
     &                         start = (/NRM(ifile,ng)%Rindex/),        &
     &                         total = (/1/),                           &
     &                         pioFile = NRM(ifile,ng)%pioFile,         &
     &                         pioVar = NRM(ifile,ng)%pioVar(idtime)%vd)
#endif
          END SELECT
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  2D norm at RHO-points.
!-----------------------------------------------------------------------
!
          IF (Cnorm(ifile,isFsur)) THEN
            MS_R2D_LOOP : DO ns=1,Nscale(ng)
              ms=ns
              iLap=Clap(isFsur,ns,ng)
              Imin=1
              Imax=Lm(ng)
              Jmin=1
              Jmax=Mm(ng)
              IF (Master) THEN
                WRITE (stdout,20) TRIM(Text),                           &
     &                          '2D normalization factors at RHO-points'
                FLUSH (stdout)
              END IF
              DO j=JstrT,JendT
                DO i=IstrT,IendT
                  Hscale(i,j)=1.0_r8/SQRT(om_r(i,j)*on_r(i,j))
                END DO
              END DO
              DO jc=Jmin,Jmax
                DO ic=Imin,Imax
#ifdef MASKING
                  compute=0.0_r8
                  IF (((Jstr.le.jc).and.(jc.le.Jend)).and.              &
     &                ((Istr.le.ic).and.(ic.le.Iend))) THEN
                    IF (rmask(ic,jc).gt.0) compute=1.0_r8
                  END IF
# ifdef DISTRIBUTE
                  CALL mp_reduce (ng, iTLM, 1, compute, 'SUM')
# endif
#else
                  compute=1.0_r8
#endif
                  IF (compute.gt.0.0_r8) THEN
                    DO j=LBj,UBj
                      DO i=LBi,UBi
                        A2d(i,j)=0.0_r8
                      END DO
                    END DO
                    IF (((Jstr.le.jc).and.(jc.le.Jend)).and.            &
     &                  ((Istr.le.ic).and.(ic.le.Iend))) THEN
                      A2d(ic,jc)=1.0_r8
                    END IF
!
                    CALL ad_dabc_r2d_tile (ng, tile,                    &
     &                                     LBi, UBi, LBj, UBj,          &
     &                                     A2d)
!
!  Implicit horizontal convolution, CG/CI solver.
!
                    CALL self%ad_CI_2d (ng, tile, iADM, isFsur,         &
     &                                  r2dvar, ms, NiterCI(ns,ng),     &
     &                                  Lweak,                          & 
     &                                  LBi, UBi, LBj, UBj,             &
     &                                  IminS, ImaxS, JminS, JmaxS,     &
     &                                  A2d)
                    DO j=JstrT,JendT
                      DO i=IstrT,IendT
                        A2d(i,j)=A2d(i,j)*Hscale(i,j)
                      END DO
                    END DO
!
                    Gdotp=dot_prod2d (ng, tile, iADM, r2dvar,           &
     &                                LBi, UBi, LBj, UBj,               &
     &                                A2d, A2d)
                    cff=1.0_r8/SQRT(Gdotp)
                  ELSE
                    cff=0.0_r8
                  END IF
                  IF (((Jstr.le.jc).and.(jc.le.Jend)).and.              &
     &                ((Istr.le.ic).and.(ic.le.Iend))) THEN
                    HnormR(ic,jc,ifile)=HnormR(ic,jc,ifile)+            &
     &                                  Bwgt(isFsur,ms,ng)*cff
                  END IF
                END DO
              END DO
            END DO MS_R2D_LOOP
!
            CALL dabc_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          HnormR(:,:,ifile))
#ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iTLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          HnormR(:,:,ifile))
#endif
!
            SELECT CASE (NRM(ifile,ng)%IOtype)
              CASE (io_nf90)
                CALL wrt_norm2d_nf90 (ng, tile, iTLM, ncname,           &
     &                                LBi, UBi, LBj, UBj, idFsur,       &
     &                                NRM(ifile,ng)%ncid,               &
     &                                NRM(ifile,ng)%Vid(idFsur),        &
     &                                NRM(ifile,ng)%Rindex,             &
#ifdef MASKING
     &                                rmask,                            &
#endif
     &                                HnormR(:,:,ifile))

#if defined PIO_LIB && defined DISTRIBUTE
              CASE (io_pio)
                IF (NRM(ifile,ng)%pioVar(idFsur)%dkind.eq.              &
     &              PIO_double) THEN
                  ioDesc => ioDesc_dp_r2dvar(ng)
                ELSE
                  ioDesc => ioDesc_sp_r2dvar(ng)
                END IF
                CALL wrt_norm2d_pio (ng, tile, iTLM, ncname,            &
     &                               LBi, UBi, LBj, UBj, idFsur,        &
     &                               NRM(ifile,ng)%pioFile,             &
     &                               NRM(ifile,ng)%pioVar(idFsur),      &
     &                               NRM(ifile,ng)%Rindex,              &
     &                               ioDesc,                            &
# ifdef MASKING
     &                               rmask,                             &
# endif
     &                               HnormR(:,:,ifile))
#endif
            END SELECT
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
!
!-----------------------------------------------------------------------
!  2D norm at U-points.
!-----------------------------------------------------------------------
!
          IF (Cnorm(ifile,isUbar)) THEN
            MS_U2D_LOOP : DO ns=1,Nscale(ng)
              ms=ns
              iLap=Clap(isUbar,ns,ng)
              IF (EWperiodic(ng)) THEN
                Imin=1
                Imax=Lm(ng)
                Jmin=1
                Jmax=Mm(ng)
              ELSE
                Imin=2
                Imax=Lm(ng)
                Jmin=1
                Jmax=Mm(ng)
              END IF
              IF (Master) THEN
                WRITE (stdout,20) TRIM(Text),                           &
     &                          '2D normalization factors at   U-points'
                FLUSH (stdout)
              END IF
              DO j=JstrT,JendT
                DO i=IstrP,IendT
                  Hscale(i,j)=1.0_r8/SQRT(om_u(i,j)*on_u(i,j))
                END DO
              END DO
              DO jc=Jmin,Jmax
                DO ic=Imin,Imax
#ifdef MASKING
                  compute=0.0_r8
                  IF (((Jstr.le.jc).and.(jc.le.Jend)).and.              &
     &                ((Istr.le.ic).and.(ic.le.Iend))) THEN
                    IF (umask(ic,jc).gt.0) compute=1.0_r8
                  END IF
# ifdef DISTRIBUTE
                  CALL mp_reduce (ng, iTLM, 1, compute, 'SUM')
# endif
#else
                  compute=1.0_r8
#endif
                  IF (compute.gt.0.0_r8) THEN
                    DO j=LBj,UBj
                      DO i=LBi,UBi
                        A2d(i,j)=0.0_r8
                      END DO
                    END DO
                    IF (((Jstr.le.jc).and.(jc.le.Jend)).and.            &
     &                  ((Istr.le.ic).and.(ic.le.Iend))) THEN
                      A2d(ic,jc)=1.0_r8
                    END IF
!
                    CALL ad_dabc_u2d_tile (ng, tile,                    &
     &                                     LBi, UBi, LBj, UBj,          &
     &                                     A2d)
!
!  Implicit horizontal convolution, CG/CI solver.
!
                    CALL self%ad_CI_2d (ng, tile, iADM, isUbar,         &
     &                                  u2dvar, ns, NiterCI(ns,ng),     &
     &                                  Lweak,                          & 
     &                                  LBi, UBi, LBj, UBj,             &
     &                                  IminS, ImaxS, JminS, JmaxS,     &
     &                                  A2d)
                    DO j=JstrT,JendT
                      DO i=IstrP,IendT
                        A2d(i,j)=A2d(i,j)*Hscale(i,j)
                      END DO
                    END DO
                    Gdotp=dot_prod2d (ng, tile, iADM, u2dvar,           &
     &                                LBi, UBi, LBj, UBj,               &
     &                                A2d, A2d)
                    cff=1.0_r8/SQRT(Gdotp)
                  ELSE
                    cff=0.0_r8
                  END IF
                  IF (((Jstr.le.jc).and.(jc.le.Jend)).and.              &
     &                ((Istr.le.ic).and.(ic.le.Iend))) THEN
                    HnormU(ic,jc,ifile)=HnormU(ic,jc,ifile)+            &
     &                                  Bwgt(isUbar,ms,ng)*cff
                  END IF
                END DO
              END DO
            END DO MS_U2D_LOOP
!
            CALL dabc_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          HnormU(:,:,ifile))
#ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iTLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          HnormU(:,:,ifile))
#endif
!
            SELECT CASE (NRM(ifile,ng)%IOtype)
              CASE (io_nf90)
                CALL wrt_norm2d_nf90 (ng, tile, iTLM, ncname,           &
     &                                LBi, UBi, LBj, UBj, idUbar,       &
     &                                NRM(ifile,ng)%ncid,               &
     &                                NRM(ifile,ng)%Vid(idUbar),        &
     &                                NRM(ifile,ng)%Rindex,             &
#ifdef MASKING
     &                                umask,                            &
#endif
     &                                HnormU(:,:,ifile))

#if defined PIO_LIB && defined DISTRIBUTE
              CASE (io_pio)
                IF (NRM(ifile,ng)%pioVar(idUbar)%dkind.eq.              &
     &              PIO_double) THEN
                  ioDesc => ioDesc_dp_u2dvar(ng)
                ELSE
                  ioDesc => ioDesc_sp_u2dvar(ng)
                END IF
                CALL wrt_norm2d_pio (ng, tile, iTLM, ncname,            &
     &                               LBi, UBi, LBj, UBj, idUbar,        &
     &                               NRM(ifile,ng)%pioFile,             &
     &                               NRM(ifile,ng)%pioVar(idUbar),      &
     &                               NRM(ifile,ng)%Rindex,              &
     &                               ioDesc,                            &
# ifdef MASKING
     &                               umask,                             &
# endif
     &                               HnormU(:,:,ifile))
#endif
            END SELECT
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
!
!-----------------------------------------------------------------------
!  2D norm at V-points.
!-----------------------------------------------------------------------
!
          IF (Cnorm(ifile,isVbar)) THEN
            MS_V2D_LOOP : DO ns=1,Nscale(ng)
              ms=ns
              iLap=Clap(isVbar,ns,ng)
              IF (NSperiodic(ng)) THEN
                Imin=1
                Imax=Lm(ng)
                Jmin=1
                Jmax=Mm(ng)
              ELSE
                Imin=1
                Imax=Lm(ng)
                Jmin=2
                Jmax=Mm(ng)
              END IF
              IF (Master) THEN
                WRITE (stdout,20) TRIM(Text),                           &
     &                          '2D normalization factors at   V-points'
                FLUSH (stdout)
              END IF
              DO j=JstrP,JendT
                DO i=IstrT,IendT
                  Hscale(i,j)=1.0_r8/SQRT(om_v(i,j)*on_v(i,j))
                END DO
              END DO
              DO jc=Jmin,Jmax
                DO ic=Imin,Imax
#ifdef MASKING
                  compute=0.0_r8
                  IF (((Jstr.le.jc).and.(jc.le.Jend)).and.              &
     &              ((Istr.le.ic).and.(ic.le.Iend))) THEN
                    IF (vmask(ic,jc).gt.0) compute=1.0_r8
                  END IF
# ifdef DISTRIBUTE
                  CALL mp_reduce (ng, iTLM, 1, compute, 'SUM')
# endif
#else
                  compute=1.0_r8
#endif
                  IF (compute.gt.0.0_r8) THEN
                    DO j=LBj,UBj
                      DO i=LBi,UBi
                        A2d(i,j)=0.0_r8
                      END DO
                    END DO
                    IF (((Jstr.le.jc).and.(jc.le.Jend)).and.            &
     &                  ((Istr.le.ic).and.(ic.le.Iend))) THEN
                      A2d(ic,jc)=1.0_r8
                    END IF
!
                    CALL ad_dabc_v2d_tile (ng, tile,                    &
     &                                     LBi, UBi, LBj, UBj,          &
     &                                     A2d)
!
!  Implicit horizontal convolution, CG/CI solver.
!
                    CALL self%ad_CI_2d (ng, tile, iADM, isVbar,         &
     &                                  v2dvar, ns, NiterCI(ns,ng),     &
     &                                  Lweak,                          & 
     &                                  LBi, UBi, LBj, UBj,             &
     &                                  IminS, ImaxS, JminS, JmaxS,     &
     &                                  A2d)
                    DO j=JstrP,JendT
                      DO i=IstrT,IendT
                        A2d(i,j)=A2d(i,j)*Hscale(i,j)
                      END DO
                    END DO
                    Gdotp=dot_prod2d (ng, tile, iADM, v2dvar,           &
     &                                LBi, UBi, LBj, UBj,               &
     &                                A2d, A2d)
                    cff=1.0_r8/SQRT(Gdotp)
!
                  ELSE
                    cff=0.0_r8
                  END IF
                  IF (((Jstr.le.jc).and.(jc.le.Jend)).and.              &
     &                ((Istr.le.ic).and.(ic.le.Iend))) THEN
                    HnormV(ic,jc,ifile)=HnormV(ic,jc,ifile)+            &
     &                                  Bwgt(isVbar,ms,ng)*cff
                  END IF
                END DO
              END DO
            END DO MS_V2D_LOOP
!
            CALL dabc_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          HnormV(:,:,ifile))
#ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iTLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          HnormV(:,:,ifile))
#endif
!
            SELECT CASE (NRM(ifile,ng)%IOtype)
              CASE (io_nf90)
                CALL wrt_norm2d_nf90 (ng, tile, iTLM, ncname,           &
     &                                LBi, UBi, LBj, UBj, idVbar,       &
     &                                NRM(ifile,ng)%ncid,               &
     &                                NRM(ifile,ng)%Vid(idVbar),        &
     &                                NRM(ifile,ng)%Rindex,             &
#ifdef MASKING
     &                                vmask,                            &
#endif
     &                                HnormV(:,:,ifile))

#if defined PIO_LIB && defined DISTRIBUTE
              CASE (io_pio)
                IF (NRM(ifile,ng)%pioVar(idVbar)%dkind.eq.              &
     &              PIO_double) THEN
                  ioDesc => ioDesc_dp_v2dvar(ng)
                ELSE
                  ioDesc => ioDesc_sp_v2dvar(ng)
                END IF
                CALL wrt_norm2d_pio (ng, tile, iTLM, ncname,            &
     &                               LBi, UBi, LBj, UBj, idVbar,        &
     &                               NRM(ifile,ng)%pioFile,             &
     &                               NRM(ifile,ng)%pioVar(idVbar),      &
     &                               NRM(ifile,ng)%Rindex,              &
     &                               ioDesc,                            &
# ifdef MASKING
     &                               vmask,                             &
# endif
     &                               HnormV(:,:,ifile))
#endif
            END SELECT
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF

#ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  3D norm U-points.
!-----------------------------------------------------------------------
!
          IF (Cnorm(ifile,isUvel)) THEN
            MS_U3D_LOOP : DO ns=1,Nscale(ng)
              ms=ns
              iLap=Clap(isUvel,ns,ng)
              IF (EWperiodic(ng)) THEN
                Imin=1
                Imax=Lm(ng)
                Jmin=1
                Jmax=Mm(ng)
              ELSE
                Imin=2
                Imax=Lm(ng)
                Jmin=1
                Jmax=Mm(ng)
              END IF
              IF (Master) THEN
                WRITE (stdout,20) TRIM(Text),                           &
     &                          '3D normalization factors at   U-points'
                FLUSH (stdout)
              END IF
              DO j=JstrT,JendT
                DO i=IstrP,IendT
                  cff=om_u(i,j)*on_u(i,j)*0.5_r8
                  DO k=1,N(ng)
                    Vscale(i,j,k)=1.0_r8/                               &
     &                            SQRT(cff*(Hz(i-1,j,k)+Hz(i,j,k)))
                  END DO
                END DO
              END DO
              DO kc=1,N(ng)
                DO jc=Jmin,Jmax
                  DO ic=Imin,Imax
# ifdef MASKING
                    compute=0.0_r8
                    IF (((Jstr.le.jc).and.(jc.le.Jend)).and.            &
     &                  ((Istr.le.ic).and.(ic.le.Iend))) THEN
                      IF (umask(ic,jc).gt.0) compute=1.0_r8
                    END IF
#  ifdef DISTRIBUTE
                    CALL mp_reduce (ng, iTLM, 1, compute, 'SUM')
#  endif
# else
                    compute=1.0_r8
# endif
                    IF (compute.gt.0.0_r8) THEN
                      DO k=1,N(ng)
                        DO j=LBj,UBj
                          DO i=LBi,UBi
                            A3d(i,j,k)=0.0_r8
                          END DO
                        END DO
                      END DO
                      IF (((Jstr.le.jc).and.(jc.le.Jend)).and.          &
     &                    ((Istr.le.ic).and.(ic.le.Iend))) THEN
                        A3d(ic,jc,kc)=1.0_r8
                      END IF
!
!  Implicit vertical convolution.
!
                      CALL self%ad_Vdiff (ng, tile, iADM, isUvel,       &
     &                                    u3dvar,                       &
     &                                    NVsteps(ifile,isUvel)/ifac,   &
     &                                    LBi, UBi, LBj, UBj,           &
     &                                    IminS, ImaxS, JminS, JmaxS,   &
     &                                    DTsizeV(ifile,isUvel),        &
     &                                    Kv, A3d)
                      CALL ad_dabc_u3d_tile (ng, tile,                  &
     &                                       LBi, UBi, LBj, UBj,        &
     &                                       A3d)
!
!  Implicit horizontal convolution, CG/CI solver.
!
                      CALL self%ad_CI_3d (ng, tile, iADM, isUvel,       &
     &                                    u3dvar, ns, NiterCI(ns,ng),   &
     &                                    Lweak,                        & 
     &                                    LBi, UBi, LBj, UBj,           &
     &                                    IminS, ImaxS, JminS, JmaxS,   &
     &                                    A3d)
!
                      DO k=1,N(ng)
                        DO j=JstrT,JendT
                          DO i=IstrP,IendT
                            A3d(i,j,k)=A3d(i,j,k)*Vscale(i,j,k)
                          END DO
                        END DO
                      END DO
                      Gdotp=dot_prod3d (ng, tile, iADM, u3dvar,         &
     &                                  LBi, UBi, LBj, UBj, 1, N(ng),   &
     &                                  A3d, A3d)
                      cff=1.0_r8/SQRT(Gdotp)
                    ELSE
                      cff=0.0_r8
                    END IF
                    IF (((Jstr.le.jc).and.(jc.le.Jend)).and.            &
     &                  ((Istr.le.ic).and.(ic.le.Iend))) THEN
                      VnormU(ic,jc,kc,ifile)=VnormU(ic,jc,kc,ifile)+    &
     &                                       Bwgt(isUvel,ms,ng)*cff
                    END IF
                  END DO
                END DO
              END DO
            END DO MS_U3D_LOOP
!
            CALL dabc_u3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          VnormU(:,:,:,ifile))
# ifdef DISTRIBUTE
            CALL mp_exchange3d (ng, tile, iTLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          VnormU(:,:,:,ifile))
# endif
!
            SELECT CASE (NRM(ifile,ng)%IOtype)
              CASE (io_nf90)
                CALL wrt_norm3d_nf90 (ng, tile, iTLM, ncname,           &
     &                                LBi, UBi, LBj, UBj, 1, N(ng),     &
     &                                idUvel, NRM(ifile,ng)%ncid,       &
     &                                NRM(ifile,ng)%Vid(idUvel),        &
     &                                NRM(ifile,ng)%Rindex,             &
# ifdef MASKING
     &                                umask,                            &
# endif
     &                                VnormU(:,:,:,ifile))

# if defined PIO_LIB && defined DISTRIBUTE
              CASE (io_pio)
                IF (NRM(ifile,ng)%pioVar(idUvel)%dkind.eq.              &
     &              PIO_double) THEN
                  ioDesc => ioDesc_dp_u3dvar(ng)
                ELSE
                  ioDesc => ioDesc_sp_u3dvar(ng)
                END IF
                CALL wrt_norm3d_pio (ng, tile, iTLM, ncname,            &
     &                               LBi, UBi, LBj, UBj, 1, N(ng),      &
     &                               idUvel, NRM(ifile,ng)%pioFile,     &
     &                               NRM(ifile,ng)%pioVar(idUvel),      &
     &                               NRM(ifile,ng)%Rindex,              &
     &                               ioDesc,                            &
#  ifdef MASKING
     &                               umask,                             &
#  endif
     &                               VnormU(:,:,:,ifile))
# endif
            END SELECT
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
!
!-----------------------------------------------------------------------
!  3D norm at V-points.
!-----------------------------------------------------------------------
!
          IF (Cnorm(ifile,isVvel)) THEN
            MS_V3D_LOOP : DO ns=1,Nscale(ng)
              ms=ns
              iLap=Clap(isVvel,ns,ng)
              IF (NSperiodic(ng)) THEN
                Imin=1
                Imax=Lm(ng)
                Jmin=1
                Jmax=Mm(ng)
              ELSE
                Imin=1
                Imax=Lm(ng)
                  Jmin=2
                Jmax=Mm(ng)
              END IF
              IF (Master) THEN
                WRITE (stdout,20) TRIM(Text),                           &
     &                          '3D normalization factors at   V-points'
                FLUSH (stdout)
              END IF
              DO j=JstrP,JendT
                DO i=IstrT,IendT
                  cff=om_v(i,j)*on_v(i,j)*0.5_r8
                  DO k=1,N(ng)
                    Vscale(i,j,k)=1.0_r8/                               &
     &                            SQRT(cff*(Hz(i,j-1,k)+Hz(i,j,k)))
                  END DO
                END DO
              END DO
              DO kc=1,N(ng)
                DO jc=Jmin,Jmax
                  DO ic=Imin,Imax
# ifdef MASKING
                    compute=0.0_r8
                    IF (((Jstr.le.jc).and.(jc.le.Jend)).and.            &
     &                  ((Istr.le.ic).and.(ic.le.Iend))) THEN
                      IF (vmask(ic,jc).gt.0) compute=1.0_r8
                    END IF
#  ifdef DISTRIBUTE
                    CALL mp_reduce (ng, iTLM, 1, compute, 'SUM')
#  endif
# else
                    compute=1.0_r8
# endif
                    IF (compute.gt.0.0_r8) THEN
                      DO k=1,N(ng)
                        DO j=LBj,UBj
                          DO i=LBi,UBi
                            A3d(i,j,k)=0.0_r8
                          END DO
                        END DO
                      END DO
                      IF (((Jstr.le.jc).and.(jc.le.Jend)).and.          &
     &                    ((Istr.le.ic).and.(ic.le.Iend))) THEN
                        A3d(ic,jc,kc)=1.0_r8
                      END IF
!
!  Implicit vertical convolution.
!
                      CALL self%ad_Vdiff (ng, tile, iADM, isVvel,       &
     &                                    v3dvar,                       &
     &                                    NVsteps(ifile,isVvel)/ifac,   &
     &                                    LBi, UBi, LBj, UBj,           &
     &                                    IminS, ImaxS, JminS, JmaxS,   &
     &                                    DTsizeV(ifile,isVvel),        &
     &                                    Kv, A3d)
                      CALL ad_dabc_v3d_tile (ng, tile,                  &
     &                                       LBi, UBi, LBj, UBj,        &
     &                                       A3d)
!
!  Implicit horizontal convolution, CG/CI solver.
!
                      CALL self%ad_CI_3d (ng, tile, iADM, isVvel,       &
     &                                    v3dvar, ns, NiterCI(ns,ng),   &
     &                                    Lweak,                        & 
     &                                    LBi, UBi, LBj, UBj,           &
     &                                    IminS, ImaxS, JminS, JmaxS,   &
     &                                    A3d)
!
                      DO k=1,N(ng)
                        DO j=JstrP,JendT
                          DO i=IstrT,IendT
                            A3d(i,j,k)=A3d(i,j,k)*Vscale(i,j,k)
                          END DO
                        END DO
                      END DO
                      Gdotp=dot_prod3d (ng, tile, iADM, v3dvar,         &
     &                                  LBi, UBi, LBj, UBj, 1, N(ng),   &
     &                                  A3d, A3d)
                      cff=1.0_r8/SQRT(Gdotp)
                    ELSE
                      cff=0.0_r8
                    END IF
                    IF (((Jstr.le.jc).and.(jc.le.Jend)).and.            &
     &                  ((Istr.le.ic).and.(ic.le.Iend))) THEN
                      VnormV(ic,jc,kc,ifile)=VnormV(ic,jc,kc,ifile)+    &
     &                                       Bwgt(isVvel,ms,ng)*cff
                    END IF
                  END DO
                END DO
              END DO
            END DO MS_V3D_LOOP
!
            CALL dabc_v3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          VnormV(:,:,:,ifile))
# ifdef DISTRIBUTE
            CALL mp_exchange3d (ng, tile, iTLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          VnormV(:,:,:,ifile))
# endif
!
            SELECT CASE (NRM(ifile,ng)%IOtype)
              CASE (io_nf90)
                CALL wrt_norm3d_nf90 (ng, tile, iTLM, ncname,           &
     &                                LBi, UBi, LBj, UBj, 1, N(ng),     &
     &                                idVvel, NRM(ifile,ng)%ncid,       &
     &                                NRM(ifile,ng)%Vid(idVvel),        &
     &                                NRM(ifile,ng)%Rindex,             &
# ifdef MASKING
     &                                vmask,                            &
# endif
     &                                VnormV(:,:,:,ifile))

# if defined PIO_LIB && defined DISTRIBUTE
              CASE (io_pio)
                IF (NRM(ifile,ng)%pioVar(idVvel)%dkind.eq.              &
     &              PIO_double) THEN
                  ioDesc => ioDesc_dp_v3dvar(ng)
                ELSE
                  ioDesc => ioDesc_sp_v3dvar(ng)
                END IF
                CALL wrt_norm3d_pio (ng, tile, iTLM, ncname,            &
     &                               LBi, UBi, LBj, UBj, 1, N(ng),      &
     &                               idVvel, NRM(ifile,ng)%pioFile,     &
     &                               NRM(ifile,ng)%pioVar(idVvel),      &
     &                               NRM(ifile,ng)%Rindex,              &
     &                               ioDesc,                            &
#  ifdef MASKING
     &                               vmask,                             &
#  endif
     &                               VnormV(:,:,:,ifile))
# endif
            END SELECT
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
!
!-----------------------------------------------------------------------
!  3D norm at RHO-points.
!-----------------------------------------------------------------------
!
          IF (Master) THEN
            Lsame=.FALSE.
            DO itrc=1,NT(ng)
              is=isTvar(itrc)
              IF (Cnorm(ifile,is)) Lsame=.TRUE.
            END DO
            IF (Lsame) THEN
              WRITE (stdout,20) TRIM(Text),                             &
     &                          '3D normalization factors at RHO-points'
              FLUSH (stdout)
            END IF
          END IF
!
!  Check if the decorrelation scales for all the tracers are different.
!  If not, just compute the normalization factors for the first tracer
!  and assign the same value to the rest.  Recall that this computation
!  is very expensive.
!
          Ldiffer=.FALSE.
          DO ns=1,Nscale(ng)
            Imin=1
            Imax=Lm(ng)
            Jmin=1
            Jmax=Mm(ng)
            DO itrc=2,NT(ng)
              IF ((HdecayX(ifile,isTvar(itrc  ),ns,ng).ne.              &
     &             HdecayX(ifile,isTvar(itrc-1),ns,ng)).or.             &
     &            (HdecayY(ifile,isTvar(itrc  ),ns,ng).ne.              &
     &             HdecayY(ifile,isTvar(itrc-1),ns,ng)).or.             &
     &            (Vdecay(ifile,isTvar(itrc  ),ng).ne.                  &
     &             Vdecay(ifile,isTvar(itrc-1),ng))) THEN
                Ldiffer=.TRUE.
              END IF
            END DO
            IF (.not.Ldiffer) THEN
              Lsame=.TRUE.
              UBt=1
            ELSE
            Lsame=.FALSE.
              UBt=NT(ng)
            END IF
          END DO
!
          DO j=JstrT,JendT
            DO i=IstrT,IendT
              cff=om_r(i,j)*on_r(i,j)
              DO k=1,N(ng)
                Vscale(i,j,k)=1.0_r8/SQRT(cff*Hz(i,j,k))
              END DO
            END DO
          END DO
!
          TRACER_LOOP : DO itrc=1,UBt
            is=isTvar(itrc)
            IF (Cnorm(ifile,is)) THEN
              MS_R3D_LOOP : DO ns=1,Nscale(ng)
                ms=ns
                iLap=Clap(is,ns,ng)
                DO kc=1,N(ng)
                  DO jc=Jmin,Jmax
                    DO ic=Imin,Imax
# ifdef MASKING
                      compute=0.0_r8
                      IF (((Jstr.le.jc).and.(jc.le.Jend)).and.          &
     &                  ((Istr.le.ic).and.(ic.le.Iend))) THEN
                        IF (rmask(ic,jc).gt.0) compute=1.0_r8
                      END IF
#  ifdef DISTRIBUTE
                      CALL mp_reduce (ng, iTLM, 1, compute, 'SUM')
#  endif
# else
                      compute=1.0_r8
# endif
                      IF (compute.gt.0.0_r8) THEN
                        DO k=1,N(ng)
                          DO j=LBj,UBj
                            DO i=LBi,UBi
                              A3d(i,j,k)=0.0_r8
                            END DO
                          END DO
                        END DO
                        IF (((Jstr.le.jc).and.(jc.le.Jend)).and.        &
     &                      ((Istr.le.ic).and.(ic.le.Iend))) THEN
                          A3d(ic,jc,kc)=1.0_r8
                        END IF
!
!  Implicit vertical convolution.
!
                        CALL self%ad_Vdiff (ng, tile, iADM, is,         &
     &                                      r3dvar,                     &
     &                                      NVsteps(ifile,is)/ifac,     &
     &                                      LBi, UBi, LBj, UBj,         &
     &                                      IminS, ImaxS, JminS, JmaxS, &
     &                                      DTsizeV(ifile,is),          &
     &                                      Kv, A3d)
                        CALL ad_dabc_r3d_tile (ng, tile,                &
     &                                         LBi, UBi, LBj, UBj,      &
     &                                         A3d)
!
!  Implicit horizontal convolution, CG/CI solver.
!
                        CALL self%ad_CI_3d (ng, tile, iADM, is,         &
     &                                      r3dvar, ns, NiterCI(ns,ng), &
     &                                      Lweak,                      & 
     &                                      LBi, UBi, LBj, UBj,         &
     &                                      IminS, ImaxS, JminS, JmaxS, &
     &                                      A3d)
!
                        DO k=1,N(ng)
                          DO j=JstrT,JendT
                            DO i=IstrT,IendT
                              A3d(i,j,k)=A3d(i,j,k)*Vscale(i,j,k)
                            END DO
                          END DO
                        END DO
                        Gdotp=dot_prod3d (ng, tile, iADM, r3dvar,       &
     &                                    LBi, UBi, LBj, UBj, 1, N(ng), &
     &                                    A3d, A3d)
                        cff=1.0_r8/SQRT(Gdotp)
                      ELSE
                        cff=0.0_r8
                      END IF
!
                      IF (((Jstr.le.jc).and.(jc.le.Jend)).and.          &
     &                    ((Istr.le.ic).and.(ic.le.Iend))) THEN
                        IF (Lsame) THEN
                          DO ntrc=1,NT(ng)
                            VnormR(ic,jc,kc,ifile,ntrc)=                &
     &                                   VnormR(ic,jc,kc,ifile,ntrc)+   &
     &                                   Bwgt(isTvar(ntrc),ms,ng)*cff
                          END DO
                        ELSE
                          VnormR(ic,jc,kc,ifile,itrc)=                  &
     &                                   VnormR(ic,jc,kc,ifile,itrc)+    &
     &                                   Bwgt(is,ms,ng)*cff
                        END IF
                      END IF
                    END DO
                  END DO
                END DO
              END DO MS_R3D_LOOP
            END IF
          END DO TRACER_LOOP
!
          DO itrc=1,NT(ng)
            is=isTvar(itrc)
            IF (Cnorm(ifile,is)) THEN
              CALL dabc_r3d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            VnormR(:,:,:,ifile,itrc))
# ifdef DISTRIBUTE
              CALL mp_exchange3d (ng, tile, iTLM, 1,                    &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            NghostPoints,                         &
     &                            EWperiodic(ng), NSperiodic(ng),       &
     &                            VnormR(:,:,:,ifile,itrc))
# endif
!
              SELECT CASE (NRM(ifile,ng)%IOtype)
                CASE (io_nf90)
                  CALL wrt_norm3d_nf90 (ng, tile, iTLM, ncname,         &
     &                                  LBi, UBi, LBj, UBj, 1, N(ng),   &
     &                                  idTvar(itrc),                   &
     &                              NRM(ifile,ng)%ncid,                 &
     &                              NRM(ifile,ng)%Vid(idTvar(itrc)),    &
     &                              NRM(ifile,ng)%Rindex,               &
# ifdef MASKING
     &                                  rmask,                          &
# endif
     &                                  VnormR(:,:,:,ifile,itrc))

# if defined PIO_LIB && defined DISTRIBUTE
                CASE (io_pio)
                  IF (NRM(ifile,ng)%pioTrc(itrc)%dkind.eq.              &
     &                PIO_double) THEN
                    ioDesc => ioDesc_dp_r3dvar(ng)
                  ELSE
                    ioDesc => ioDesc_sp_r3dvar(ng)
                  END IF
                  CALL wrt_norm3d_pio (ng, tile, iTLM, ncname,          &
     &                                 LBi, UBi, LBj, UBj, 1, N(ng),    &
     &                                 idTvar(itrc),                    &
     &                              NRM(ifile,ng)%pioFile,              &
     &                              NRM(ifile,ng)%pioTrc(itrc),         &
     &                              NRM(ifile,ng)%Rindex,               &
     &                                  ioDesc,                         &
#  ifdef MASKING
     &                                  rmask,                          &
#  endif
     &                                  VnormR(:,:,:,ifile,itrc))
# endif
              END SELECT
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END IF
          END DO
#endif
        END IF
      END DO FILE_LOOP

#ifdef ADJUST_BOUNDARY
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Compute open boundaries error covariance, B, normalization factors
!  using the exact method.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      ifile=3
      IF (LwrtNRM(ifile,ng)) THEN
        Text='boundary conditions'
        IJlen=UBij-LBij+1
# ifdef SOLVE3D
        IJKlen=IJlen*N(ng)
# endif
        Lconvolve(iwest )=DOMAIN(ng)%Western_Edge (tile)
        Lconvolve(ieast )=DOMAIN(ng)%Eastern_Edge (tile)
        Lconvolve(isouth)=DOMAIN(ng)%Southern_Edge(tile)
        Lconvolve(inorth)=DOMAIN(ng)%Northern_Edge(tile)
!
!  Set time record index to write in normalization NetCDF file.
!
        ncname=NRM(ifile,ng)%name
        NRM(ifile,ng)%Rindex=NRM(ifile,ng)%Rindex+1
        NRM(ifile,ng)%Nrec=NRM(ifile,ng)%Nrec+1
!
!  Write out model time (s).
!
        my_time=tdays(ng)*day2sec

        SELECT CASE (NRM(ifile,ng)%IOtype)
          CASE (io_nf90)
            CALL netcdf_put_fvar (ng, iTLM, ncname,                     &
     &                            Vname(1,idtime), my_time,             &
     &                         start = (/NRM(ifile,ng)%Rindex/),        &
     &                         total = (/1/),                           &
     &                         ncid = NRM(ifile,ng)%ncid,               &
     &                         varid = NRM(ifile,ng)%Vid(idtime))

# if defined PIO_LIB && defined DISTRIBUTE
          CASE (io_pio)
            CALL pio_netcdf_put_fvar (ng, iTLM, ncname,                 &
     &                                Vname(1,idtime), my_time,         &
     &                         start = (/NRM(ifile,ng)%Rindex/),        &
     &                         total = (/1/),                           &
     &                         pioFile = NRM(ifile,ng)%pioFile,         &
     &                         pioVar = NRM(ifile,ng)%pioVar(idtime)%vd)
# endif
        END SELECT
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  2D boundary norm at RHO-points.
!-----------------------------------------------------------------------
!
        HnormRobc=Aspv

        IF (Master.and.ANY(CnormB(isFsur,:))) THEN
          WRITE (stdout,20) TRIM(Text),                                 &
     &                      '2D normalization factors at RHO-points'
          FLUSH (stdout)
        END IF

        DO ibry=1,4
          HscaleB=0.0_r8
          IF (CnormB(isFsur,ibry)) THEN
            DO ns=1,Nscale(ng)
              ms=ns
              iLap=Clap(isFsur,ns,ng)
              IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
                i=BOUNDS(ng)%edge(ibry,r2dvar)
                Bmin=1
                Bmax=Mm(ng)
                IF (Lconvolve(ibry)) THEN
                  DO j=JstrT,JendT
                    HscaleB(j)=1.0_r8/SQRT(on_r(i,j))
                  END DO
                END IF
              ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
                j=BOUNDS(ng)%edge(ibry,r2dvar)
                Bmin=1
                Bmax=Lm(ng)
                IF (Lconvolve(ibry)) THEN
                  DO i=IstrT,IendT
                    HscaleB(i)=1.0_r8/SQRT(om_r(i,j))
                  END DO
                END IF
              END IF
              DO ib=Bmin,Bmax
                IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
                  bounded=Lconvolve(ibry).and.                          &
     &                    ((Jstr.le.ib).and.(ib.le.Jend))
                  j=ib
                ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
                  bounded=Lconvolve(ibry).and.                          &
     &                    ((Istr.le.ib).and.(ib.le.Iend))
                  i=ib
                END IF
#  ifdef MASKING
                IF (bounded) THEN
                  compute=rmask(i,j)
                ELSE
                  compute=0.0_r8
                END IF
#   ifdef DISTRIBUTE
                CALL mp_reduce (ng, iTLM, 1, compute, 'SUM')
#   endif
#  else
                compute=1.0_r8
#  endif
                IF (compute.gt.0.0_r8) THEN
                  B2d=0.0_r8
                  IF (bounded) THEN
                    B2d(ib)=1.0_r8
                  END IF
                  CALL ad_CI_bry1d_tile (ng, tile,                      &
     &                                   LBi, UBi, LBj, UBj,            &
     &                                   LBij, UBij,                    &
     &                                   IminS, ImaxS, JminS, JmaxS,    &
     &                                   nstp, nnew, ifac,              &
     &                                   isFsur, ibry,                  &
     &                                   ms, MiterCI(ns,ng), iLap,      &
     &                                   r2dvar, 1, 1,                  &
     &                                   emaxr1d(:,ibry,ns),            &
     &                                   eminr1d(:,ibry,ns),            &
     &                                   ci_rb, ci_pb, ci_qb, ci_xb,    &
# ifdef SOLVE3D
     &                                   Hz,                            &
# endif
     &                                   B2d)
!
!  HscaleB must be applied twice.
!
                  IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
                    DO j=JstrT,JendT
                      B2d(j)=B2d(j)*HscaleB(j)
                    END DO
                  ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
                    DO i=IstrT,IendT
                      B2d(i)=B2d(i)*HscaleB(i)
                    END DO
                  END IF
!
                  IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
                    DO j=JstrT,JendT
                      B2d(j)=B2d(j)*HscaleB(j)
                    END DO
                  ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
                    DO i=IstrT,IendT
                      B2d(i)=B2d(i)*HscaleB(i)
                    END DO
                  END IF
                  CALL tl_CI_bry1d_tile (ng, tile,                      &
     &                                   LBi, UBi, LBj, UBj,            &
     &                                   LBij, UBij,                    &
     &                                   IminS, ImaxS, JminS, JmaxS,    &
     &                                   nstp, nnew, ifac,              &
     &                                   isFsur, ibry,                  &
     &                                   ms, MiterCI(ns,ng), iLap,      &
     &                                   r2dvar, 1, 1,                  &
     &                                   emaxr1d(:,ibry,ns),            &
     &                                   eminr1d(:,ibry,ns),            &
     &                                   ci_rb, ci_pb, ci_qb, ci_xb,    &
# ifdef SOLVE3D
     &                                   Hz,                            &
# endif
     &                                   B2d)
                  IF (bounded) THEN
                    cff=1.0_r8/SQRT(B2d(ib))
                  END IF
                ELSE
                  cff=0.0_r8
                END IF
                IF (bounded) THEN
                  HnormRobc(ib,ibry)=HnormRobc(ib,ibry)+                &
     &                                 Bwgt(isFsur,ms,ng)*cff
                END IF
              END DO
            END DO       ! end of ns-loop.
            CALL bc_r2d_bry_tile (ng, tile, ibry,                       &
     &                            LBij, UBij,                           &
     &                            HnormRobc(:,ibry))
# ifdef DISTRIBUTE
            CALL mp_collect (ng, iTLM, IJlen, Aspv,                     &
     &                       HnormRobc(LBij:,ibry))
# endif
          END IF
        END DO
        IF (ANY(CnormB(isFsur,:))) THEN
          ifield=idSbry(isFsur)

          SELECT CASE (NRM(ifile,ng)%IOtype)
            CASE (io_nf90)
              CALL netcdf_put_fvar (ng, iTLM, ncname,                   &
     &                              Vname(1,ifield),                    &
     &                              HnormRobc(LBij:,:),                 &
     &                         start = (/1,1,NRM(ifile,ng)%Rindex/),    &
     &                         total = (/IJlen,4,1/),                   &
     &                         ncid = NRM(ifile,ng)%ncid,               &
     &                         varid = NRM(ifile,ng)%Vid(ifield))

# if defined PIO_LIB && defined DISTRIBUTE
            CASE (io_pio)
              CALL pio_netcdf_put_fvar (ng, iTLM, ncname,               &
     &                                  Vname(1,ifield),                &
     &                                  HnormRobc(LBij:,:),             &
     &                         start = (/1,1,NRM(ifile,ng)%Rindex/),    &
     &                         total = (/IJlen,4,1/),                   &
     &                         pioFile = NRM(ifile,ng)%pioFile,         &
     &                         pioVar = NRM(ifile,ng)%pioVar(ifield)%vd)

# endif
          END SELECT
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!-----------------------------------------------------------------------
!  2D boundary norm at U-points.
!-----------------------------------------------------------------------
!
        HnormUobc=Aspv

        IF (Master.and.ANY(CnormB(isUbar,:))) THEN
          WRITE (stdout,20) TRIM(Text),                                 &
     &                      '2D normalization factors at   U-points'
          FLUSH (stdout)
        END IF

        DO ibry=1,4
          HscaleB=0.0_r8
          IF (CnormB(isUbar,ibry)) THEN
            DO ns=1,Nscale(ng)
              ms=ns
              iLap=Clap(isUbar,ns,ng)
              IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
                i=BOUNDS(ng)%edge(ibry,u2dvar)
                Bmin=1
                Bmax=Mm(ng)
                IF (Lconvolve(ibry)) THEN
                  DO j=JstrT,JendT
                    HscaleB(j)=1.0_r8/SQRT(on_u(i,j))
                  END DO
                END IF
              ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
                j=BOUNDS(ng)%edge(ibry,u2dvar)
                IF (EWperiodic(ng)) THEN
                  Bmin=1
                  Bmax=Lm(ng)
                ELSE
                  Bmin=2
                  Bmax=Lm(ng)
                END IF
                IF (Lconvolve(ibry)) THEN
                  DO i=IstrP,IendT
                    HscaleB(i)=1.0_r8/SQRT(om_u(i,j))
                  END DO
                END IF
              END IF
              DO ib=Bmin,Bmax
                IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
                  bounded=Lconvolve(ibry).and.                          &
     &                    ((Jstr.le.ib).and.(ib.le.Jend))
                  j=ib
                ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
                  bounded=Lconvolve(ibry).and.                          &
     &                    ((Istr.le.ib).and.(ib.le.Iend))
                  i=ib
                END IF
#  ifdef MASKING
                IF (bounded) THEN
                  compute=umask(i,j)
                ELSE
                  compute=0.0_r8
                END IF
#   ifdef DISTRIBUTE
                CALL mp_reduce (ng, iTLM, 1, compute, 'SUM')
#   endif
#  else
                compute=1.0_r8
#  endif
                IF (compute.gt.0.0_r8) THEN
                  B2d=0.0_r8
                  IF (bounded) THEN
                    B2d(ib)=1.0_r8
                  END IF
                  CALL ad_CI_bry1d_tile (ng, tile,                      &
     &                                   LBi, UBi, LBj, UBj,            &
     &                                   LBij, UBij,                    &
     &                                   IminS, ImaxS, JminS, JmaxS,    &
     &                                   nstp, nnew, ifac,              &
     &                                   isUbar, ibry,                  &
     &                                   ms, MiterCI(ns,ng), iLap,      &
     &                                   u2dvar, 1, 1,                  &
     &                                   emaxu1d(:,ibry,ns),            &
     &                                   eminu1d(:,ibry,ns),            &
     &                                   ci_rb, ci_pb, ci_qb, ci_xb,    &
#  ifdef SOLVE3D
     &                                   Hz,                            &
#  endif
     &                                   B2d)
!
!  HscaleB must be applied twice.
!
                  IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
                    DO j=JstrT,JendT
                      B2d(j)=B2d(j)*HscaleB(j)
                    END DO
                  ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
                    DO i=IstrP,IendT
                      B2d(i)=B2d(i)*HscaleB(i)
                    END DO
                  END IF
!
                  IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
                    DO j=JstrT,JendT
                      B2d(j)=B2d(j)*HscaleB(j)
                    END DO
                  ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
                    DO i=IstrP,IendT
                      B2d(i)=B2d(i)*HscaleB(i)
                    END DO
                  END IF
                  CALL tl_CI_bry1d_tile (ng, tile,                      &
     &                                   LBi, UBi, LBj, UBj,            &
     &                                   LBij, UBij,                    &
     &                                   IminS, ImaxS, JminS, JmaxS,    &
     &                                   nstp, nnew, ifac,              &
     &                                   isUbar, ibry,                  &
     &                                   ms, MiterCI(ns,ng), iLap,      &
     &                                   u2dvar, 1, 1,                  &
     &                                   emaxu1d(:,ibry,ns),            &
     &                                   eminu1d(:,ibry,ns),            &
     &                                   ci_rb, ci_pb, ci_qb, ci_xb,    &
#  ifdef SOLVE3D
     &                                   Hz,                            &
#  endif
     &                                   B2d)
                  IF (bounded) THEN
                    cff=1.0_r8/SQRT(B2d(ib))
                  END IF
                ELSE
                  cff=0.0_r8
                END IF
                IF (bounded) THEN
                  HnormUobc(ib,ibry)=HnormUobc(ib,ibry)+                &
     &                               Bwgt(isUbar,ms,ng)*cff
                END IF
              END DO
            END DO    ! end of ns-loop
            CALL bc_u2d_bry_tile (ng, tile, ibry,                       &
     &                            LBij, UBij,                           &
     &                            HnormUobc(:,ibry))
# ifdef DISTRIBUTE
            CALL mp_collect (ng, iTLM, IJlen, Aspv,                     &
     &                       HnormUobc(LBij:,ibry))
# endif
          END IF
        END DO
        IF (ANY(CnormB(isUbar,:))) THEN
          ifield=idSbry(isUbar)

          SELECT CASE (NRM(ifile,ng)%IOtype)
            CASE (io_nf90)
              CALL netcdf_put_fvar (ng, iTLM, ncname,                   &
     &                              Vname(1,ifield),                    &
     &                              HnormUobc(LBij:,:),                 &
     &                         start = (/1,1,NRM(ifile,ng)%Rindex/),    &
     &                         total = (/IJlen,4,1/),                   &
     &                         ncid = NRM(ifile,ng)%ncid,               &
     &                         varid = NRM(ifile,ng)%Vid(ifield))

# if defined PIO_LIB && defined DISTRIBUTE
            CASE (io_pio)
              CALL pio_netcdf_put_fvar (ng, iTLM, ncname,               &
     &                                  Vname(1,ifield),                &
     &                                  HnormUobc(LBij:,:),             &
     &                         start = (/1,1,NRM(ifile,ng)%Rindex/),    &
     &                         total = (/IJlen,4,1/),                   &
     &                         pioFile = NRM(ifile,ng)%pioFile,         &
     &                         pioVar = NRM(ifile,ng)%pioVar(ifield)%vd)
# endif
          END SELECT
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!-----------------------------------------------------------------------
!  2D boundary norm at V-points.
!-----------------------------------------------------------------------
!
        HnormVobc=Aspv

        IF (Master.and.ANY(CnormB(isVbar,:))) THEN
          WRITE (stdout,20) TRIM(Text),                                 &
     &                      '2D normalization factors at   V-points'
          FLUSH (stdout)
        END IF

        DO ibry=1,4
          HscaleB=0.0_r8
          IF (CnormB(isVbar,ibry)) THEN
            DO ns=1,Nscale(ng)
              ms=ns
              iLap=Clap(isVbar,ns,ng)
              IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
                i=BOUNDS(ng)%edge(ibry,v2dvar)
                IF (NSperiodic(ng)) THEN
                  Bmin=1
                  Bmax=Mm(ng)
                ELSE
                  Bmin=2
                  Bmax=Mm(ng)
                END IF
                IF (Lconvolve(ibry)) THEN
                  DO j=JstrT,JendT
                    HscaleB(j)=1.0_r8/SQRT(on_v(i,j))
                  END DO
                END IF
              ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
                j=BOUNDS(ng)%edge(ibry,v2dvar)
                Bmin=1
                Bmax=Lm(ng)
                IF (Lconvolve(ibry)) THEN
                  DO i=IstrT,IendT
                    HscaleB(i)=1.0_r8/SQRT(om_v(i,j))
                  END DO
                END IF
              END IF
              DO ib=Bmin,Bmax
                IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
                  bounded=Lconvolve(ibry).and.                          &
     &                    ((Jstr.le.ib).and.(ib.le.Jend))
                  j=ib
                ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
                  bounded=Lconvolve(ibry).and.                          &
     &                    ((Istr.le.ib).and.(ib.le.Iend))
                  i=ib
                END IF
#  ifdef MASKING
                IF (bounded) THEN
                  compute=vmask(i,j)
                ELSE
                  compute=0.0_r8
                END IF
#   ifdef DISTRIBUTE
                CALL mp_reduce (ng, iTLM, 1, compute, 'SUM')
#   endif
#  else
                compute=1.0_r8
#  endif
                IF (compute.gt.0.0_r8) THEN
                  B2d=0.0_r8
                  IF (bounded) THEN
                    B2d(ib)=1.0_r8
                  END IF
                  CALL ad_CI_bry1d_tile (ng, tile,                      &
     &                                   LBi, UBi, LBj, UBj,            &
     &                                   LBij, UBij,                    &
     &                                   IminS, ImaxS, JminS, JmaxS,    &
     &                                   nstp, nnew, ifac,              &
     &                                   isVbar, ibry,                  &
     &                                   ms, MiterCI(ns,ng), iLap,      &
     &                                   v2dvar, 1, 1,                  &
     &                                   emaxv1d(:,ibry,ns),            &
     &                                   eminv1d(:,ibry,ns),            &
     &                                   ci_rb, ci_pb, ci_qb, ci_xb,    &
# ifdef SOLVE3D
     &                                   Hz,                            &
# endif
     &                                   B2d)
!
!  HscaleB must be applied twice.
!
                  IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
                    DO j=JstrP,JendT
                      B2d(j)=B2d(j)*HscaleB(j)
                    END DO
                  ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
                    DO i=IstrT,IendT
                      B2d(i)=B2d(i)*HscaleB(i)
                    END DO
                  END IF
!
                  IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
                    DO j=JstrP,JendT
                      B2d(j)=B2d(j)*HscaleB(j)
                    END DO
                  ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
                    DO i=IstrT,IendT
                      B2d(i)=B2d(i)*HscaleB(i)
                    END DO
                  END IF
                  CALL tl_CI_bry1d_tile (ng, tile,                      &
     &                                   LBi, UBi, LBj, UBj,            &
     &                                   LBij, UBij,                    &
     &                                   IminS, ImaxS, JminS, JmaxS,    &
     &                                   nstp, nnew, ifac,              &
     &                                   isVbar, ibry,                  &
     &                                   ms, MiterCI(ns,ng), iLap,      &
     &                                   v2dvar, 1, 1,                  &
     &                                   emaxv1d(:,ibry,ns),            &
     &                                   eminv1d(:,ibry,ns),            &
     &                                   ci_rb, ci_pb, ci_qb, ci_xb,    &
# ifdef SOLVE3D
     &                                   Hz,                            &
# endif
     &                                   B2d)
                  IF (bounded) THEN
                    cff=1.0_r8/SQRT(B2d(ib))
                  END IF
                ELSE
                  cff=0.0_r8
                END IF
                IF (bounded) THEN
                  HnormVobc(ib,ibry)=HnormVobc(ib,ibry)+                &
     &                               Bwgt(isVbar,ms,ng)*cff
                END IF
              END DO
            END DO   ! end of ns-loop
            CALL bc_v2d_bry_tile (ng, tile, ibry,                       &
     &                            LBij, UBij,                           &
     &                            HnormVobc(:,ibry))
# ifdef DISTRIBUTE
            CALL mp_collect (ng, iTLM, IJlen, Aspv,                     &
     &                       HnormVobc(LBij:,ibry))
# endif
          END IF
        END DO
        IF (ANY(CnormB(isVbar,:))) THEN
          ifield=idSbry(isVbar)

          SELECT CASE (NRM(ifile,ng)%IOtype)
            CASE (io_nf90)
              CALL netcdf_put_fvar (ng, iTLM, ncname,                   &
     &                              Vname(1,ifield),                    &
     &                              HnormVobc(LBij:,:),                 &
     &                         start = (/1,1,NRM(ifile,ng)%Rindex/),    &
     &                         total = (/IJlen,4,1/),                   &
     &                         ncid = NRM(ifile,ng)%ncid,               &
     &                         varid = NRM(ifile,ng)%Vid(ifield))

# if defined PIO_LIB && defined DISTRIBUTE
            CASE (io_pio)
              CALL pio_netcdf_put_fvar (ng, iTLM, ncname,               &
     &                                  Vname(1,ifield),                &
     &                                  HnormVobc(LBij:,:),             &
     &                         start = (/1,1,NRM(ifile,ng)%Rindex/),    &
     &                         total = (/IJlen,4,1/),                   &
     &                         pioFile = NRM(ifile,ng)%pioFile,         &
     &                         pioVar = NRM(ifile,ng)%pioVar(ifield)%vd)
# endif
          END SELECT
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

# ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  3D boundary norm at U-points.
!-----------------------------------------------------------------------
!
        VnormUobc=Aspv

        IF (Master.and.ANY(CnormB(isUvel,:))) THEN
          WRITE (stdout,20) TRIM(Text),                                 &
     &                      '3D normalization factors at   U-points'
          FLUSH (stdout)
        END IF

        DO ibry=1,4
          VscaleB=0.0_r8
          IF (CnormB(isUvel,ibry)) THEN
            DO ns=1,Nscale(ng)
              ms=ns
              iLap=Clap(isUvel,ns,ng)
              IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
                i=BOUNDS(ng)%edge(ibry,u2dvar)
                Bmin=1
                Bmax=Mm(ng)
                IF (Lconvolve(ibry)) THEN
                  DO j=JstrT,JendT
                    cff=on_u(i,j)*0.5_r8
                    DO k=1,N(ng)
                      VscaleB(j,k)=1.0_r8/                              &
     &                             SQRT(cff*(Hz(i-1,j,k)+Hz(i,j,k)))
                    END DO
                  END DO
                END IF
              ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
                j=BOUNDS(ng)%edge(ibry,u2dvar)
                IF (EWperiodic(ng)) THEN
                  Bmin=1
                  Bmax=Lm(ng)
                ELSE
                  Bmin=2
                  Bmax=Lm(ng)
                END IF
                IF (Lconvolve(ibry)) THEN
                  DO i=IstrP,IendT
                    cff=om_u(i,j)*0.5_r8
                    DO k=1,N(ng)
                      VscaleB(i,k)=1.0_r8/                              &
     &                             SQRT(cff*(Hz(i-1,j,k)+Hz(i,j,k)))
                    END DO
                  END DO
                END IF
              END IF
              DO kb=1,N(ng)
                DO ib=Bmin,Bmax
                  IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
                    bounded=Lconvolve(ibry).and.                        &
     &                      ((Jstr.le.ib).and.(ib.le.Jend))
                    j=ib
                  ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
                    bounded=Lconvolve(ibry).and.                        &
     &                      ((Istr.le.ib).and.(ib.le.Iend))
                    i=ib
                  END IF
#   ifdef MASKING
                  IF (bounded) THEN
                    compute=umask(i,j)
                  ELSE
                    compute=0.0_r8
                  END IF
#    ifdef DISTRIBUTE
                  CALL mp_reduce (ng, iTLM, 1, compute, 'SUM')
#    endif
#   else
                  compute=1.0_r8
#   endif
                  IF (compute.gt.0.0_r8) THEN
                    B3d=0.0_r8
                    IF (bounded) THEN
                      B3d(ib,kb)=1.0_r8
                    END IF
!
!  First the vertical convolution.
!
                    CALL ad_conv_u3d_bry_tile (ng, tile, iADM, ibry,    &
     &                                       BOUNDS(ng)%edge(:,u2dvar), &
     &                                       LBij, UBij,                &
     &                                       LBi, UBi, LBj, UBj,        &
     &                                       1, N(ng),                  &
     &                                       IminS, ImaxS, JminS, JmaxS,&
     &                                       NghostPoints,              &
     &                                       NHstepsB(ibry,isUvel)/ifac,&
     &                                       NVstepsB(ibry,isUvel)/ifac,&
     &                                       DTsizeHB(ibry,isUvel),     &
     &                                       DTsizeVB(ibry,isUvel),     &
     &                                       Kh, Kv,                    &
     &                                       pm, pn, pmon_r, pnom_p,    &
#  ifdef MASKING
     &                                       umask, pmask,              &
#  endif
     &                                       Hz, z_r,                   &
     &                                       B3d)
!
!  Now the horizontal convolution.
!
                    DO k=1,N(ng)
                      klev=k
                      CALL ad_CI_bry1d_tile (ng, tile,                  &
     &                                       LBi, UBi, LBj, UBj,        &
     &                                       LBij, UBij,                &
     &                                       IminS, ImaxS, JminS, JmaxS,&
     &                                       nstp, nnew, ifac,          &
     &                                       isUvel, ibry,              &
     &                                       ms, MiterCI(ns,ng), iLap,  &
     &                                       u3dvar, klev, 1,           &
     &                                       emaxu1dz(:,klev,ibry,ns),  &
     &                                       eminu1dz(:,klev,ibry,ns),  &
     &                                       ci_rb, ci_pb,              &
     &                                       ci_qb, ci_xb,              &
# ifdef SOLVE3D
     &                                       Hz,                        &
# endif
     &                                       B3d(:,klev))
                    END DO
!
!  VscaleB must be applied twice.
!
                    IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
                      DO k=1,N(ng)
                        DO j=JstrT,JendT
                          B3d(j,k)=B3d(j,k)*VscaleB(j,k)
                        END DO
                      END DO
                    ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
                      DO k=1,N(ng)
                        DO i=IstrP,IendT
                          B3d(i,k)=B3d(i,k)*VscaleB(i,k)
                        END DO
                      END DO
                    END IF
!
                    IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
                      DO k=1,N(ng)
                        DO j=JstrT,JendT
                          B3d(j,k)=B3d(j,k)*VscaleB(j,k)
                        END DO
                      END DO
                    ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
                      DO k=1,N(ng)
                        DO i=IstrP,IendT
                          B3d(i,k)=B3d(i,k)*VscaleB(i,k)
                        END DO
                      END DO
                    END IF
!
!  First the horizontal convolution.
!
                    DO k=1,N(ng)
                      klev=k
                      CALL tl_CI_bry1d_tile (ng, tile,                  &
     &                                       LBi, UBi, LBj, UBj,        &
     &                                       LBij, UBij,                &
     &                                       IminS, ImaxS, JminS, JmaxS,&
     &                                       nstp, nnew, ifac,          &
     &                                       isUvel, ibry,              &
     &                                       ms, MiterCI(ns,ng), iLap,  &
     &                                       u3dvar, klev, 1,           &
     &                                       emaxu1dz(:,klev,ibry,ns),  &
     &                                       eminu1dz(:,klev,ibry,ns),  &
     &                                       ci_rb, ci_pb,              &
     &                                       ci_qb, ci_xb,              &
# ifdef SOLVE3D
     &                                       Hz,                        &
# endif
     &                                       B3d(:,klev))
                    END DO
!
!  Now the vertical convolution.
!
                    CALL tl_conv_u3d_bry_tile (ng, tile, iTLM, ibry,    &
     &                                       BOUNDS(ng)%edge(:,u2dvar), &
     &                                       LBij, UBij,                &
     &                                       LBi, UBi, LBj, UBj,        &
     &                                       1, N(ng),                  &
     &                                       IminS, ImaxS, JminS, JmaxS,&
     &                                       NghostPoints,              &
     &                                       NHstepsB(ibry,isUvel)/ifac,&
     &                                       NVstepsB(ibry,isUvel)/ifac,&
     &                                       DTsizeHB(ibry,isUvel),     &
     &                                       DTsizeVB(ibry,isUvel),     &
     &                                       Kh, Kv,                    &
     &                                       pm, pn, pmon_r, pnom_p,    &
#  ifdef MASKING
     &                                       umask, pmask,              &
#  endif
     &                                       Hz, z_r,                   &
     &                                       B3d)
!
                    IF (bounded) THEN
                      cff=1.0_r8/SQRT(B3d(ib,kb))
                    END IF
                  ELSE
                    cff=0.0_r8
                  END IF
                  IF (bounded) THEN
                    VnormUobc(ib,kb,ibry)=VnormUobc(ib,kb,ibry)+        &
     &                                    Bwgt(isUvel,ms,ng)*cff
                  END IF
                END DO
              END DO
            END DO    ! end of ns-loop
            CALL bc_u3d_bry_tile (ng, tile, ibry,                       &
     &                            LBij, UBij, 1, N(ng),                 &
     &                            VnormUobc(:,:,ibry))
#  ifdef DISTRIBUTE
            Bwrk=RESHAPE(VnormUobc(:,:,ibry), (/IJKlen/))
            CALL mp_collect (ng, iTLM, IJKlen, Aspv, Bwrk)
            ic=0
            DO k=1,N(ng)
              DO ib=LBij,UBij
                ic=ic+1
                VnormUobc(ib,k,ibry)=Bwrk(ic)
              END DO
            END DO
#  endif
          END IF
        END DO
        IF (ANY(CnormB(isUvel,:))) THEN
          ifield=idSbry(isUvel)

          SELECT CASE (NRM(ifile,ng)%IOtype)
            CASE (io_nf90)
              CALL netcdf_put_fvar (ng, iTLM, ncname,                   &
     &                              Vname(1,ifield),                    &
     &                              VnormUobc(LBij:,:,:),               &
     &                         start = (/1,1,1,NRM(ifile,ng)%Rindex/),  &
     &                         total = (/IJlen,N(ng),4,1/),             &
     &                         ncid = NRM(ifile,ng)%ncid,               &
     &                         varid = NRM(ifile,ng)%Vid(ifield))

#  if defined PIO_LIB && defined DISTRIBUTE
            CASE (io_pio)
              CALL pio_netcdf_put_fvar (ng, iTLM, ncname,               &
     &                                  Vname(1,ifield),                &
     &                                  VnormUobc(LBij:,:,:),           &
     &                         start = (/1,1,1,NRM(ifile,ng)%Rindex/),  &
     &                         total = (/IJlen,N(ng),4,1/),             &
     &                         pioFile = NRM(ifile,ng)%pioFile,         &
     &                         pioVar = NRM(ifile,ng)%pioVar(ifield)%vd)
#  endif
          END SELECT
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!-----------------------------------------------------------------------
!  3D boundary norm at V-points.
!-----------------------------------------------------------------------
!
        VnormVobc=Aspv

        IF (Master.and.ANY(CnormB(isVvel,:))) THEN
          WRITE (stdout,20) TRIM(Text),                                 &
     &                      '3D normalization factors at   V-points'
          FLUSH (stdout)
        END IF

        DO ibry=1,4
          VscaleB=0.0_r8
          IF (CnormB(isVvel,ibry)) THEN
            DO ns=1,Nscale(ng)
              ms=ns
              iLap=Clap(isVvel,ns,ng)
              IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
                i=BOUNDS(ng)%edge(ibry,v2dvar)
                IF (NSperiodic(ng)) THEN
                  Bmin=1
                  Bmax=Mm(ng)
                ELSE
                  Bmin=2
                  Bmax=Mm(ng)
                END IF
                IF (Lconvolve(ibry)) THEN
                  DO j=JstrP,JendT
                    cff=on_v(i,j)*0.5_r8
                    DO k=1,N(ng)
                      VscaleB(j,k)=1.0_r8/                              &
     &                             SQRT(cff*(Hz(i,j-1,k)+Hz(i,j,k)))
                    END DO
                  END DO
                END IF
              ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
                j=BOUNDS(ng)%edge(ibry,v2dvar)
                Bmin=1
                Bmax=Lm(ng)
                IF (Lconvolve(ibry)) THEN
                  DO i=IstrT,IendT
                    cff=om_v(i,j)*0.5_r8
                    DO k=1,N(ng)
                      VscaleB(i,k)=1.0_r8/                              &
     &                             SQRT(cff*(Hz(i,j-1,k)+Hz(i,j,k)))
                    END DO
                  END DO
                END IF
              END IF
              DO kb=1,N(ng)
                DO ib=Bmin,Bmax
                  IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
                    bounded=Lconvolve(ibry).and.                        &
     &                      ((Jstr.le.ib).and.(ib.le.Jend))
                    j=ib
                  ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
                    bounded=Lconvolve(ibry).and.                        &
     &                      ((Istr.le.ib).and.(ib.le.Iend))
                    i=ib
                  END IF
#   ifdef MASKING
                  IF (bounded) THEN
                    compute=vmask(i,j)
                  ELSE
                    compute=0.0_r8
                  END IF
#    ifdef DISTRIBUTE
                  CALL mp_reduce (ng, iTLM, 1, compute, 'SUM')
#    endif
#   else
                  compute=1.0_r8
#   endif
                  IF (compute.gt.0.0_r8) THEN
                    B3d=0.0_r8
                    IF (bounded) THEN
                      B3d(ib,kb)=1.0_r8
                    END IF
!
!  First the vertical convolution.
!
                    CALL ad_conv_v3d_bry_tile (ng, tile, iADM, ibry,    &
     &                                       BOUNDS(ng)%edge(:,v2dvar), &
     &                                       LBij, UBij,                &
     &                                       LBi, UBi, LBj, UBj,        &
     &                                       1, N(ng),                  &
     &                                       IminS, ImaxS, JminS, JmaxS,&
     &                                       NghostPoints,              &
     &                                       NHstepsB(ibry,isVvel)/ifac,&
     &                                       NVstepsB(ibry,isVvel)/ifac,&
     &                                       DTsizeHB(ibry,isVvel),     &
     &                                       DTsizeVB(ibry,isVvel),     &
     &                                       Kh, Kv,                    &
     &                                       pm, pn, pmon_p, pnom_r,    &
#  ifdef MASKING
     &                                       vmask, pmask,              &
#  endif
     &                                       Hz, z_r,                   &
     &                                       B3d)
!
!  Now the horizontal convolution.
!
                    DO k=1,N(ng)
                      klev=k
                      CALL ad_CI_bry1d_tile (ng, tile,                  &
     &                                       LBi, UBi, LBj, UBj,        &
     &                                       LBij, UBij,                &
     &                                       IminS, ImaxS, JminS, JmaxS,&
     &                                       nstp, nnew, ifac,          &
     &                                       isVvel, ibry,              &
     &                                       ms, MiterCI(ns,ng), iLap,  &
     &                                       v3dvar, klev, 1,           &
     &                                       emaxv1dz(:,klev,ibry,ns),  &
     &                                       eminv1dz(:,klev,ibry,ns),  &
     &                                       ci_rb, ci_pb,              &
     &                                       ci_qb, ci_xb,              &
# ifdef SOLVE3D
     &                                       Hz,                        &
# endif
     &                                       B3d(:,klev))
                    END DO
!
!  VscaleB must be applied twice.
!
                    IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
                      DO k=1,N(ng)
                        DO j=JstrP,JendT
                          B3d(j,k)=B3d(j,k)*VscaleB(j,k)
                        END DO
                      END DO
                    ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
                      DO k=1,N(ng)
                        DO i=IstrT,IendT
                          B3d(i,k)=B3d(i,k)*VscaleB(i,k)
                        END DO
                      END DO
                    END IF
!
                    IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
                      DO k=1,N(ng)
                        DO j=JstrP,JendT
                          B3d(j,k)=B3d(j,k)*VscaleB(j,k)
                        END DO
                      END DO
                    ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
                      DO k=1,N(ng)
                        DO i=IstrT,IendT
                          B3d(i,k)=B3d(i,k)*VscaleB(i,k)
                        END DO
                      END DO
                    END IF
!
!  First the horizontal convolution.
!
                    DO k=1,N(ng)
                      klev=k
                      CALL tl_CI_bry1d_tile (ng, tile,                  &
     &                                       LBi, UBi, LBj, UBj,        &
     &                                       LBij, UBij,                &
     &                                       IminS, ImaxS, JminS, JmaxS,&
     &                                       nstp, nnew, ifac,          &
     &                                       isVvel, ibry,              &
     &                                       ms, MiterCI(ns,ng), iLap,  &
     &                                       v3dvar, klev, 1,           &
     &                                       emaxv1dz(:,klev,ibry,ns),  &
     &                                       eminv1dz(:,klev,ibry,ns),  &
     &                                       ci_rb, ci_pb,              &
     &                                       ci_qb, ci_xb,              &
# ifdef SOLVE3D
     &                                       Hz,                        &
# endif
     &                                       B3d(:,klev))
                    END DO
!
!  Now the vertical convolution.
!
                    CALL tl_conv_v3d_bry_tile (ng, tile, iTLM, ibry,    &
     &                                       BOUNDS(ng)%edge(:,v2dvar), &
     &                                       LBij, UBij,                &
     &                                       LBi, UBi, LBj, UBj,        &
     &                                       1, N(ng),                  &
     &                                       IminS, ImaxS, JminS, JmaxS,&
     &                                       NghostPoints,              &
     &                                       NHstepsB(ibry,isVvel)/ifac,&
     &                                       NVstepsB(ibry,isVvel)/ifac,&
     &                                       DTsizeHB(ibry,isVvel),     &
     &                                       DTsizeVB(ibry,isVvel),     &
     &                                       Kh, Kv,                    &
     &                                       pm, pn, pmon_p, pnom_r,    &
#  ifdef MASKING
     &                                       vmask, pmask,              &
#  endif
     &                                       Hz, z_r,                   &
     &                                       B3d)
                    IF (bounded) THEN
                      cff=1.0_r8/SQRT(B3d(ib,kb))
                    END IF
                  ELSE
                    cff=0.0_r8
                  END IF
                  IF (bounded) THEN
                    VnormVobc(ib,kb,ibry)=VnormVobc(ib,kb,ibry)+        &
     &                                    Bwgt(isVvel,ms,ng)*cff
                  END IF
                END DO
              END DO
            END DO    ! end of ns-loop
            CALL bc_v3d_bry_tile (ng, tile, ibry,                       &
     &                            LBij, UBij, 1, N(ng),                 &
     &                            VnormVobc(:,:,ibry))
#  ifdef DISTRIBUTE
            Bwrk=RESHAPE(VnormVobc(:,:,ibry), (/IJKlen/))
            CALL mp_collect (ng, iTLM, IJKlen, Aspv, Bwrk)
            ic=0
            DO k=1,N(ng)
              DO ib=LBij,UBij
                ic=ic+1
                VnormVobc(ib,k,ibry)=Bwrk(ic)
              END DO
            END DO
#  endif
          END IF
        END DO
        IF (ANY(CnormB(isVvel,:))) THEN
          ifield=idSbry(isVvel)

          SELECT CASE (NRM(ifile,ng)%IOtype)
            CASE (io_nf90)
              CALL netcdf_put_fvar (ng, iTLM, ncname,                   &
     &                              Vname(1,ifield),                    &
     &                              VnormVobc(LBij:,:,:),               &
     &                         start = (/1,1,1,NRM(ifile,ng)%Rindex/),  &
     &                         total = (/IJlen,N(ng),4,1/),             &
     &                         ncid = NRM(ifile,ng)%ncid,               &
     &                         varid = NRM(ifile,ng)%Vid(ifield))

#  if defined PIO_LIB && defined DISTRIBUTE
            CASE (io_pio)
              CALL pio_netcdf_put_fvar (ng, iTLM, ncname,               &
     &                                  Vname(1,ifield),                &
     &                                  VnormVobc(LBij:,:,:),           &
     &                         start = (/1,1,1,NRM(ifile,ng)%Rindex/),  &
     &                         total = (/IJlen,N(ng),4,1/),             &
     &                         pioFile = NRM(ifile,ng)%pioFile,         &
     &                         pioVar = NRM(ifile,ng)%pioVar(ifield)%vd)
#  endif
          END SELECT
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!-----------------------------------------------------------------------
!  3D boundary norm at RHO-points.
!-----------------------------------------------------------------------
!
        IF (Master) THEN
          DO itrc=1,NT(ng)
            is=isTvar(itrc)
            IF (ANY(CnormB(is,:))) THEN
              Lsame=.TRUE.
              EXIT
            END IF
          END DO
          IF (Lsame) THEN
            WRITE (stdout,20) TRIM(Text),                               &
     &                        '3D normalization factors at RHO-points'
            FLUSH (stdout)
          END IF
        END IF

        DO itrc=1,NT(ng)
          VnormRobc=Aspv
          is=isTvar(itrc)
          DO ibry=1,4
            VscaleB=0.0_r8
            IF (CnormB(is,ibry)) THEN
              DO ns=1,Nscale(ng)
                ms=ns
                iLap=Clap(is,ns,ng)
                IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
                  i=BOUNDS(ng)%edge(ibry,r2dvar)
                  Bmin=2
                  Bmax=Mm(ng)
                  IF (Lconvolve(ibry)) THEN
                    DO j=JstrT,JendT
                      cff=on_r(i,j)
                      DO k=1,N(ng)
                        VscaleB(j,k)=1.0_r8/SQRT(cff*Hz(i,j,k))
                      END DO
                    END DO
                  END IF
                ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
                  j=BOUNDS(ng)%edge(ibry,r2dvar)
                  Bmin=1
                  Bmax=Lm(ng)
                  IF (Lconvolve(ibry)) THEN
                    DO i=IstrT,IendT
                      cff=om_r(i,j)
                      DO k=1,N(ng)
                        VscaleB(i,k)=1.0_r8/SQRT(cff*Hz(i,j,k))
                      END DO
                    END DO
                  END IF
                END IF
                DO kb=1,N(ng)
                  DO ib=Bmin,Bmax
                    IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
                      bounded=Lconvolve(ibry).and.                      &
     &                        ((Jstr.le.ib).and.(ib.le.Jend))
                      j=ib
                    ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
                      bounded=Lconvolve(ibry).and.                      &
     &                        ((Istr.le.ib).and.(ib.le.Iend))
                      i=ib
                    END IF
#   ifdef MASKING
                    IF (bounded) THEN
                      compute=rmask(i,j)
                    ELSE
                      compute=0.0_r8
                    END IF
#    ifdef DISTRIBUTE
                    CALL mp_reduce (ng, iTLM, 1, compute, 'SUM')
#    endif
#   else
                    compute=1.0_r8
#   endif
                    IF (compute.gt.0.0_r8) THEN
                      B3d=0.0_r8
                      IF (bounded) THEN
                        B3d(ib,kb)=1.0_r8
                      END IF
!
!  First the vertical convolution.
!
                      CALL ad_conv_r3d_bry_tile (ng, tile, iADM, ibry,  &
     &                                       BOUNDS(ng)%edge(:,r2dvar), &
     &                                           LBij, UBij,            &
     &                                           LBi, UBi, LBj, UBj,    &
     &                                           1, N(ng),              &
     &                                           IminS, ImaxS,          &
     &                                           JminS, JmaxS,          &
     &                                           NghostPoints,          &
     &                                           NHstepsB(ibry,is)/ifac,&
     &                                           NVstepsB(ibry,is)/ifac,&
     &                                           DTsizeHB(ibry,is),     &
     &                                           DTsizeVB(ibry,is),     &
     &                                           Kh, Kv,                &
     &                                           pm, pn,                &
     &                                           pmon_u, pnom_v,        &
#   ifdef MASKING
     &                                           rmask, umask, vmask,   &
#   endif
     &                                           Hz, z_r,               &
     &                                           B3d)
!
!  Now the horizontal convolution.
!
                      DO k=1,N(ng)
                        klev=k
                        CALL ad_CI_bry1d_tile (ng, tile,                &
     &                                         LBi, UBi, LBj, UBj,      &
     &                                         LBij, UBij,              &
     &                                         IminS, ImaxS,            &
     &                                         JminS, JmaxS,            &
     &                                         nstp, nnew, ifac,        &
     &                                         is, ibry,                &
     &                                         ms, MiterCI(ns,ng),iLap, &
     &                                         r3dvar, klev, itrc,      &
     &                                   emaxr1dz(:,klev,itrc,ibry,ns), &
     &                                   eminr1dz(:,klev,itrc,ibry,ns), &
     &                                         ci_rb, ci_pb,            &
     &                                         ci_qb, ci_xb,            &
# ifdef SOLVE3D
     &                                         Hz,                      &
# endif
     &                                         B3d(:,klev))
                      END DO
!
!  VscaleB must be applied twice.
!
                      IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
                        DO k=1,N(ng)
                          DO j=JstrT,JendT
                            B3d(j,k)=B3d(j,k)*VscaleB(j,k)
                          END DO
                        END DO
                      ELSE IF ((ibry.eq.isouth).or.                     &
     &                         (ibry.eq.inorth)) THEN
                        DO k=1,N(ng)
                          DO i=IstrT,IendT
                            B3d(i,k)=B3d(i,k)*VscaleB(i,k)
                          END DO
                        END DO
                      END IF
                      IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
                        DO k=1,N(ng)
                          DO j=JstrT,JendT
                            B3d(j,k)=B3d(j,k)*VscaleB(j,k)
                          END DO
                        END DO
                      ELSE IF ((ibry.eq.isouth).or.                     &
     &                         (ibry.eq.inorth)) THEN
                        DO k=1,N(ng)
                          DO i=IstrT,IendT
                            B3d(i,k)=B3d(i,k)*VscaleB(i,k)
                          END DO
                        END DO
                      END IF
!
!  First the horizontal convolution.
!
                      DO k=1,N(ng)
                        klev=k
                        CALL tl_CI_bry1d_tile (ng, tile,                &
     &                                         LBi, UBi, LBj, UBj,      &
     &                                         LBij, UBij,              &
     &                                         IminS, ImaxS,            &
     &                                         JminS, JmaxS,            &
     &                                         nstp, nnew, ifac,        &
     &                                         is, ibry,                &
     &                                         ms, MiterCI(ns,ng), iLap,&
     &                                         r3dvar, klev, itrc,      &
     &                                   emaxr1dz(:,klev,itrc,ibry,ns), &
     &                                   eminr1dz(:,klev,itrc,ibry,ns), &
     &                                         ci_rb, ci_pb,            &
     &                                         ci_qb, ci_xb,            &
# ifdef SOLVE3D
     &                                         Hz,                      &
# endif
     &                                         B3d(:,klev))
                      END DO
!
!  Now the vertical convolution.
!
                      CALL tl_conv_r3d_bry_tile (ng, tile, iTLM, ibry,  &
     &                                       BOUNDS(ng)%edge(:,r2dvar), &
     &                                           LBij, UBij,            &
     &                                           LBi, UBi, LBj, UBj,    &
     &                                           1, N(ng),              &
     &                                           IminS, ImaxS,          &
     &                                           JminS, JmaxS,          &
     &                                           NghostPoints,          &
     &                                           NHstepsB(ibry,is)/ifac,&
     &                                           NVstepsB(ibry,is)/ifac,&
     &                                           DTsizeHB(ibry,is),     &
     &                                           DTsizeVB(ibry,is),     &
     &                                           Kh, Kv,                &
     &                                           pm, pn,                &
     &                                           pmon_u, pnom_v,        &
#   ifdef MASKING
     &                                           rmask, umask, vmask,   &
#   endif
     &                                           Hz, z_r,               &
     &                                           B3d)
                      IF (bounded) THEN
                        cff=1.0_r8/SQRT(B3d(ib,kb))
                      END IF
                    ELSE
                      cff=0.0_r8
                    END IF
                    IF (bounded) THEN
                      VnormRobc(ib,kb,ibry,itrc)=                       &
     &                                VnormRobc(ib,kb,ibry,itrc)+       &
     &                                Bwgt(isTvar(itrc),ms,ng)*cff
                    END IF
                  END DO
                END DO
              END DO   !  end of ns-loop
              CALL bc_r3d_bry_tile (ng, tile, ibry,                     &
     &                              LBij, UBij, 1, N(ng),               &
     &                              VnormRobc(:,:,ibry,itrc))
#  ifdef DISTRIBUTE
              Bwrk=RESHAPE(VnormRobc(:,:,ibry,itrc), (/IJKlen/))
              CALL mp_collect (ng, iTLM, IJKlen, Aspv, Bwrk)
              ic=0
              DO k=1,N(ng)
                DO ib=LBij,UBij
                  ic=ic+1
                  VnormRobc(ib,k,ibry,itrc)=Bwrk(ic)
                END DO
              END DO
#  endif
            END IF
          END DO
          IF (ANY(CnormB(is,:))) THEN
            ifield=idSbry(isTvar(itrc))

            SELECT CASE (NRM(ifile,ng)%IOtype)
              CASE (io_nf90)
                CALL netcdf_put_fvar (ng, iTLM, ncname,                 &
     &                                Vname(1,ifield),                  &
     &                                VnormRobc(LBij:,:,:,itrc),        &
     &                         start =(/1,1,1,NRM(ifile,ng)%Rindex/),   &
     &                         total = (/IJlen,N(ng),4,1/),             &
     &                         ncid = NRM(ifile,ng)%ncid,               &
     &                         varid = NRM(ifile,ng)%Vid(ifield))

#  if defined PIO_LIB && defined DISTRIBUTE
              CASE (io_pio)
                CALL pio_netcdf_put_fvar (ng, iTLM, ncname,             &
     &                                    Vname(1,ifield),              &
     &                                    VnormRobc(LBij:,:,:,itrc),    &
     &                         start =(/1,1,1,NRM(ifile,ng)%Rindex/),   &
     &                         total = (/IJlen,N(ng),4,1/),             &
     &                         pioFile = NRM(ifile,ng)%pioFile,         &
     &                         pioVar = NRM(ifile,ng)%pioVar(ifield)%vd)
#  endif
            END SELECT
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
# endif
!
!  Synchronize open boundaries normalization NetCDF file to disk to
!  allow other processes to access data immediately after it is
!  written.
!
        SELECT CASE (NRM(ifile,ng)%IOtype)
          CASE (io_nf90)
            CALL netcdf_sync (ng, iTLM, ncname,                         &
     &                        NRM(ifile,ng)%ncid)
# if defined PIO_LIB && defined DISTRIBUTE
          CASE (io_pio)
            CALL pio_netcdf_sync (ng, iTLM, ncname,                     &
     &                            NRM(ifile,ng)%pioFile)
# endif
        END SELECT
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#endif

#if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Compute surface forcing error covariance, B, normalization factors
!  using the exact method.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      ifile=4
      IF (LwrtNRM(ifile,ng)) THEN
        rec=1
        Text='surface forcing'
!
!  Set time record index to write in normalization NetCDF file.
!
        ncname=NRM(ifile,ng)%name
        NRM(ifile,ng)%Rindex=NRM(ifile,ng)%Rindex+1
        NRM(ifile,ng)%Nrec=NRM(ifile,ng)%Nrec+1
!
!  Write out model time (s).
!
        my_time=tdays(ng)*day2sec

        SELECT CASE (NRM(ifile,ng)%IOtype)
          CASE (io_nf90)
            CALL netcdf_put_fvar (ng, iTLM, ncname,                     &
     &                            Vname(1,idtime), my_time,             &
     &                         start = (/NRM(ifile,ng)%Rindex/),        &
     &                         total = (/1/),                           &
     &                         ncid = NRM(ifile,ng)%ncid,               &
     &                         varid = NRM(ifile,ng)%Vid(idtime))

# if defined PIO_LIB && defined DISTRIBUTE
          CASE (io_pio)
            CALL pio_netcdf_put_fvar (ng, iTLM, ncname,                 &
     &                                Vname(1,idtime), my_time,         &
     &                         start = (/NRM(ifile,ng)%Rindex/),        &
     &                         total = (/1/),                           &
     &                         pioFile = NRM(ifile,ng)%pioFile,         &
     &                         pioVar = NRM(ifile,ng)%pioVar(idtime)%vd)
# endif
        END SELECT
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef ADJUST_WSTRESS
!
!-----------------------------------------------------------------------
!  2D norm at U-stress points.
!-----------------------------------------------------------------------
!
        IF (Cnorm(rec,isUstr)) THEN
          DO ns=1,Nscale(ng)
            ms=ns
            iLap=Clap(isUstr,ns,ng)
            IF (EWperiodic(ng)) THEN
              Imin=1
              Imax=Lm(ng)
              Jmin=1
              Jmax=Mm(ng)
            ELSE
              Imin=2
              Imax=Lm(ng)
              Jmin=1
              Jmax=Mm(ng)
            END IF
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Text),                             &
     &                     '2D normalization factors at U-stress points'
              FLUSH (stdout)
            END IF
            DO j=JstrT,JendT
              DO i=IstrP,IendT
                Hscale(i,j)=1.0_r8/SQRT(om_u(i,j)*on_u(i,j))
              END DO
            END DO
            DO jc=Jmin,Jmax
              DO ic=Imin,Imax
#   ifdef MASKING
                compute=0.0_r8
                IF (((Jstr.le.jc).and.(jc.le.Jend)).and.                &
     &              ((Istr.le.ic).and.(ic.le.Iend))) THEN
                  IF (umask(ic,jc).gt.0) compute=1.0_r8
                END IF
#    ifdef DISTRIBUTE
                CALL mp_reduce (ng, iTLM, 1, compute, 'SUM')
#    endif
#   else
                compute=1.0_r8
#   endif
                IF (compute.gt.0.0_r8) THEN
                  DO j=LBj,UBj
                    DO i=LBi,UBi
                      A2d(i,j)=0.0_r8
                    END DO
                  END DO
                  IF (((Jstr.le.jc).and.(jc.le.Jend)).and.              &
     &                ((Istr.le.ic).and.(ic.le.Iend))) THEN
                    A2d(ic,jc)=1.0_r8
                  END IF
                  CALL ad_CI_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj, LBij, UBij,      &
     &                             IminS, ImaxS, JminS, JmaxS,          &
     &                             nstp, nnew, Lweak, ifac, isUstr,     &
     &                             ms, MiterCI(ns,ng), iLap,            &
     &                             u2dvar, 1, 1,                        &
     &                             emaxu2ds(:,ns), eminu2ds(:,ns),      &
     &                             ci_r, ci_p, ci_q, ci_x,              &
#   ifdef READ_SCALES
     &                             bussclx(:,:,ns), busscly(:,:,ns),    &
#   endif
#   ifdef SOLVE3D
     &                             Hz,                                  &
#   endif
     &                             A2d)
                  DO j=JstrT,JendT
                    DO i=IstrP,IendT
                      A2d(i,j)=A2d(i,j)*Hscale(i,j)
                    END DO
                  END DO
!
                  my_dot=0.0_r8
                  DO j=JstrT,JendT
                    DO i=IstrP,IendT
                      my_dot=my_dot+A2d(i,j)*A2d(i,j)
                    END DO
                  END DO
!
!  Perform parallel global reduction operation: dot product.
!
#  ifdef DISTRIBUTE
                NSUB=1                           ! distributed-memory
#  else
                  IF (DOMAIN(ng)%SouthWest_Corner(tile).and.            &
     &                DOMAIN(ng)%NorthEast_Corner(tile)) THEN
                    NSUB=1                         ! non-tiled application
                  ELSE
                    NSUB=NtileX(ng)*NtileE(ng)     ! tiled application
                  END IF
#  endif
!$OMP CRITICAL (USTR_DOT)
                  IF (tile_count.eq.0) THEN
                    Gdotp=my_dot
                  ELSE
                    Gdotp=Gdotp+my_dot
                  END IF
                  tile_count=tile_count+1
                  IF (tile_count.eq.NSUB) THEN
                    tile_count=0
#  ifdef DISTRIBUTE
                    op_handle='SUM'
                    CALL mp_reduce (ng, iTLM, 1, Gdotp, op_handle)
#  endif
                    cff=1.0_r8/SQRT(Gdotp)
                  END IF
!$OMP END CRITICAL (USTR_DOT)
                ELSE
                  cff=0.0_r8
                END IF
                IF (((Jstr.le.jc).and.(jc.le.Jend)).and.                &
     &            ((Istr.le.ic).and.(ic.le.Iend))) THEN
                  HnormSUS(ic,jc)=HnormSUS(ic,jc)+Bwgt(isUstr,ms,ng)*cff
                END IF
              END DO
            END DO
          END DO   ! end of ns-loop
          CALL dabc_u2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        HnormSUS)
#  ifdef DISTRIBUTE
          CALL mp_exchange2d (ng, tile, iTLM, 1,                        &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        HnormSUS)
#  endif
!
         SELECT CASE (NRM(ifile,ng)%IOtype)
            CASE (io_nf90)
              CALL wrt_norm2d_nf90 (ng, tile, iTLM, ncname,             &
     &                              LBi, UBi, LBj, UBj, idUsms,         &
     &                              NRM(ifile,ng)%ncid,                 &
     &                              NRM(ifile,ng)%Vid(idUsms),          &
     &                              NRM(ifile,ng)%Rindex,               &
#  ifdef MASKING
     &                              umask,                              &
#  endif
     &                              HnormSUS)

#  if defined PIO_LIB && defined DISTRIBUTE
            CASE (io_pio)
              IF (NRM(ifile,ng)%pioVar(idUsms)%dkind.eq.                &
     &            PIO_double) THEN
                ioDesc => ioDesc_dp_u2dvar(ng)
              ELSE
                ioDesc => ioDesc_sp_u2dvar(ng)
              END IF
              CALL wrt_norm2d_pio (ng, tile, iTLM, ncname,              &
     &                             LBi, UBi, LBj, UBj, idUsms,          &
     &                             NRM(ifile,ng)%pioFile,               &
     &                             NRM(ifile,ng)%pioVar(idUsms),        &
     &                             NRM(ifile,ng)%Rindex,                &
     &                             ioDesc,                              &
#   ifdef MASKING
     &                             umask,                               &
#   endif
     &                             HnormSUS)
#  endif
          END SELECT
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!-----------------------------------------------------------------------
!  2D norm at V-stress points.
!-----------------------------------------------------------------------
!
        IF (Cnorm(rec,isVstr)) THEN
          DO ns=1,Nscale(ng)
            ms=ns
            iLap=Clap(isVstr,ns,ng)
            IF (NSperiodic(ng)) THEN
              Imin=1
              Imax=Lm(ng)
              Jmin=1
              Jmax=Mm(ng)
            ELSE
              Imin=1
              Imax=Lm(ng)
              Jmin=2
              Jmax=Mm(ng)
            END IF
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Text),                             &
     &                     '2D normalization factors at V-stress points'
              FLUSH (stdout)
            END IF
            DO j=JstrP,JendT
              DO i=IstrT,IendT
                Hscale(i,j)=1.0_r8/SQRT(om_v(i,j)*on_v(i,j))
              END DO
            END DO
            DO jc=Jmin,Jmax
              DO ic=Imin,Imax
#   ifdef MASKING
                compute=0.0_r8
                IF (((Jstr.le.jc).and.(jc.le.Jend)).and.                &
     &              ((Istr.le.ic).and.(ic.le.Iend))) THEN
                  IF (vmask(ic,jc).gt.0) compute=1.0_r8
                END IF
#    ifdef DISTRIBUTE
                CALL mp_reduce (ng, iTLM, 1, compute, 'SUM')
#    endif
#   else
                compute=1.0_r8
#   endif
                IF (compute.gt.0.0_r8) THEN
                  DO j=LBj,UBj
                    DO i=LBi,UBi
                      A2d(i,j)=0.0_r8
                    END DO
                  END DO
                  IF (((Jstr.le.jc).and.(jc.le.Jend)).and.              &
     &                ((Istr.le.ic).and.(ic.le.Iend))) THEN
                    A2d(ic,jc)=1.0_r8
                  END IF
                  CALL ad_CI_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj, LBij, UBij,      &
     &                             IminS, ImaxS, JminS, JmaxS,          &
     &                             nstp, nnew, Lweak, ifac, isVstr,     &
     &                             ms, MiterCI(ns,ng), iLap,            &
     &                             v2dvar, 1, 1,                        &
     &                             emaxv2ds(:,ns), eminv2ds(:,ns),      &
     &                             ci_r, ci_p, ci_q, ci_x,              &
#    ifdef READ_SCALES
     &                             bvssclx(:,:,ns), bvsscly(:,:,ns),    &
#    endif
#    ifdef SOLVE3D
     &                             Hz,                                  &
#    endif
     &                             A2d)
                  DO j=JstrP,JendT
                    DO i=IstrT,IendT
                      A2d(i,j)=A2d(i,j)*Hscale(i,j)
                    END DO
                  END DO
!
                  my_dot=0.0_r8
                  DO j=JstrP,JendT
                    DO i=IstrT,IendT
                      my_dot=my_dot+A2d(i,j)*A2d(i,j)
                    END DO
                  END DO
!
!  Perform parallel global reduction operation: dot product.
!
#  ifdef DISTRIBUTE
                  NSUB=1                           ! distributed-memory
#  else
                  IF (DOMAIN(ng)%SouthWest_Corner(tile).and.            &
     &              DOMAIN(ng)%NorthEast_Corner(tile)) THEN
                    NSUB=1                         ! non-tiled application
                  ELSE
                    NSUB=NtileX(ng)*NtileE(ng)     ! tiled application
                  END IF
#  endif
!$OMP CRITICAL (VSTR_DOT)
                  IF (tile_count.eq.0) THEN
                    Gdotp=my_dot
                  ELSE
                    Gdotp=Gdotp+my_dot
                  END IF
                  tile_count=tile_count+1
                  IF (tile_count.eq.NSUB) THEN
                    tile_count=0
#  ifdef DISTRIBUTE
                    op_handle='SUM'
                    CALL mp_reduce (ng, iTLM, 1, Gdotp, op_handle)
#  endif
                    cff=1.0_r8/SQRT(Gdotp)
                  END IF
!$OMP END CRITICAL (VSTR_DOT)
                ELSE
                  cff=0.0_r8
                END IF
                IF (((Jstr.le.jc).and.(jc.le.Jend)).and.                &
     &              ((Istr.le.ic).and.(ic.le.Iend))) THEN
                  HnormSVS(ic,jc)=HnormSVS(ic,jc)+Bwgt(isVstr,ms,ng)*cff
                END IF
              END DO
            END DO
          END DO  ! end of ns-loop
          CALL dabc_v2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        HnormSVS)
#  ifdef DISTRIBUTE
          CALL mp_exchange2d (ng, tile, iTLM, 1,                        &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        HnormSVS)
#  endif
!
          SELECT CASE (NRM(ifile,ng)%IOtype)
            CASE (io_nf90)
              CALL wrt_norm2d_nf90 (ng, tile, iTLM, ncname,             &
     &                              LBi, UBi, LBj, UBj, idVsms,         &
     &                              NRM(ifile,ng)%ncid,                 &
     &                              NRM(ifile,ng)%Vid(idVsms),          &
     &                              NRM(ifile,ng)%Rindex,               &
#  ifdef MASKING
     &                              vmask,                              &
#  endif
     &                              HnormSVS)

#  if defined PIO_LIB && defined DISTRIBUTE
            CASE (io_pio)
              IF (NRM(ifile,ng)%pioVar(idVsms)%dkind.eq.                &
     &            PIO_double) THEN
                ioDesc => ioDesc_dp_v2dvar(ng)
              ELSE
                ioDesc => ioDesc_sp_v2dvar(ng)
              END IF
              CALL wrt_norm2d_pio (ng, tile, iTLM, ncname,              &
     &                             LBi, UBi, LBj, UBj, idVsms,          &
     &                             NRM(ifile,ng)%pioFile,               &
     &                             NRM(ifile,ng)%pioVar(idVsms),        &
     &                             NRM(ifile,ng)%Rindex,                &
     &                             ioDesc,                              &
#   ifdef MASKING
     &                             vmask,                               &
#   endif
     &                             HnormSVS)
#  endif
          END SELECT
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
# endif

# if defined ADJUST_STFLUX && defined SOLVE3D
!
!-----------------------------------------------------------------------
!  2D norm at surface tracer fluxes points.
!-----------------------------------------------------------------------
!
        IF (Master) THEN
          Lsame=.FALSE.
          DO itrc=1,NT(ng)
            IF (Lstflux(itrc,ng)) THEN
              is=isTsur(itrc)
              IF (Cnorm(rec,is)) Lsame=.TRUE.
            END IF
          END DO
          IF (Lsame) THEN
            WRITE (stdout,20) TRIM(Text),                               &
                              '2D normalization factors at RHO-points'
            FLUSH (stdout)
          END IF
        END IF
!
!  Check if the decorrelation scales for all the surface tracer fluxes
!  are different. If not, just compute the normalization factors for the
!  first tracer and assign the same value to the rest.  Recall that this
!  computation is very expensive.
!
        Ldiffer=.FALSE.
        DO ns=1,Nscale(ng)
          Imin=1
          Imax=Lm(ng)
          Jmin=1
          Jmax=Mm(ng)
          DO itrc=2,NT(ng)
            IF ((HdecayX(rec,isTsur(itrc  ),ns,ng).ne.                  &
     &           HdecayX(rec,isTsur(itrc-1),ns,ng)).or.                 &
     &          (HdecayY(rec,isTsur(itrc  ),ns,ng).ne.                  &
     &           HdecayY(rec,isTsur(itrc-1),ns,ng))) THEN
              Ldiffer=.TRUE.
            END IF
          END DO
          IF (.not.Ldiffer) THEN
            Lsame=.TRUE.
            UBt=1
          ELSE
            Lsame=.FALSE.
            UBt=NT(ng)
          END IF
        END DO
!
        DO j=JstrT,JendT
          DO i=IstrT,IendT
            Hscale(i,j)=1.0_r8/SQRT(om_r(i,j)*on_r(i,j))
          END DO
        END DO
!
        DO itrc=1,UBt
          IF (Lstflux(itrc,ng)) THEN
            ifield=isTsur(itrc)
            IF (Cnorm(rec,ifield)) THEN
              DO ns=1,Nscale(ng)
                ms=ns
                iLap=Clap(ifield,ns,ng)
                DO jc=Jmin,Jmax
                  DO ic=Imin,Imax
#   ifdef MASKING
                    compute=0.0_r8
                    IF (((Jstr.le.jc).and.(jc.le.Jend)).and.            &
     &                  ((Istr.le.ic).and.(ic.le.Iend))) THEN
                      IF (rmask(ic,jc).gt.0) compute=1.0_r8
                    END IF
#    ifdef DISTRIBUTE
                    CALL mp_reduce (ng, iTLM, 1, compute, 'SUM')
#    endif
#   else
                    compute=1.0_r8
#   endif
                    IF (compute.gt.0.0_r8) THEN
                      DO j=LBj,UBj
                        DO i=LBi,UBi
                          A2d(i,j)=0.0_r8
                        END DO
                      END DO
                      IF (((Jstr.le.jc).and.(jc.le.Jend)).and.          &
     &                    ((Istr.le.ic).and.(ic.le.Iend))) THEN
                        A2d(ic,jc)=1.0_r8
                      END IF
                      CALL ad_CI_tile (ng, tile,                        &
     &                                 LBi, UBi, LBj, UBj, LBij, UBij,  &
     &                                 IminS, ImaxS, JminS, JmaxS,      &
     &                                 nstp, nnew, Lweak, ifac, ifield, &
     &                                 ms, MiterCI(ns,ng), iLap,        &
     &                                 r2dvar, 1, 1,                    &
     &                                 emaxr2ds(:,itrc,ns),             &
     &                                 eminr2ds(:,itrc,ns),             &
     &                                 ci_r, ci_p, ci_q, ci_x,          &
#   ifdef READ_SCALES
     &                                 btfsclx(:,:,itrc,ns),            &
     &                                 btfscly(:,:,itrc,ns),            &
#   endif
#   ifdef SOLVE3D
     &                                 Hz,                              &
#   endif
     &                                 A2d)
                      DO j=JstrT,JendT
                        DO i=IstrT,IendT
                          A2d(i,j)=A2d(i,j)*Hscale(i,j)
                        END DO
                      END DO
                      my_dot=0.0_r8
                      DO j=JstrT,JendT
                        DO i=IstrT,IendT
                          my_dot=my_dot+A2d(i,j)*A2d(i,j)
                        END DO
                      END DO
!
!  Perform parallel global reduction operation: dot product.
!
#  ifdef DISTRIBUTE
                      NSUB=1                       ! distributed-memory
#  else
                      IF (DOMAIN(ng)%SouthWest_Corner(tile).and.        &
     &                    DOMAIN(ng)%NorthEast_Corner(tile)) THEN
                        NSUB=1                     ! non-tiled application
                      ELSE
                        NSUB=NtileX(ng)*NtileE(ng)     ! tiled application
                      END IF
#  endif
!$OMP CRITICAL (STFLX_DOT)
                      IF (tile_count.eq.0) THEN
                        Gdotp=my_dot
                      ELSE
                        Gdotp=Gdotp+my_dot
                      END IF
                      tile_count=tile_count+1
                      IF (tile_count.eq.NSUB) THEN
                        tile_count=0
#  ifdef DISTRIBUTE
                        op_handle='SUM'
                        CALL mp_reduce (ng, iTLM, 1, Gdotp, op_handle)
#  endif
                        cff=1.0_r8/SQRT(Gdotp)
                      END IF
!$OMP END CRITICAL (STFLX_DOT)
                    ELSE
                      cff=0.0_r8
                    END IF
                    IF (((Jstr.le.jc).and.(jc.le.Jend)).and.            &
     &                  ((Istr.le.ic).and.(ic.le.Iend))) THEN
                      IF (Lsame) THEN
                        DO ntrc=1,NT(ng)
                          nfield=isTsur(ntrc)
                          IF (Lstflux(ntrc,ng)) THEN
                            HnormSTF(ic,jc,ntrc)=HnormSTF(ic,jc,ntrc)+  &
     &                                           Bwgt(nfield,ms,ng)*cff
                          END IF
                        END DO
                      ELSE
                        HnormSTF(ic,jc,itrc)=HnormSTF(ic,jc,itrc)+      &
     &                                       Bwgt(ifield,ms,ng)*cff
                      END IF
                    END IF
                  END DO
                END DO
              END DO   ! end of ns-loop
            END IF
          END IF
        END DO
!
        DO itrc=1,NT(ng)
          IF (Lstflux(itrc,ng)) THEN
            is=isTsur(itrc)
            IF (Cnorm(rec,is)) THEN
              CALL dabc_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            HnormSTF(:,:,itrc))
#  ifdef DISTRIBUTE
              CALL mp_exchange2d (ng, tile, iTLM, 1,                    &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            NghostPoints,                         &
     &                            EWperiodic(ng), NSperiodic(ng),       &
     &                            HnormSTF(:,:,itrc))
#  endif
!
              SELECT CASE (NRM(ifile,ng)%IOtype)
                CASE (io_nf90)
                  CALL wrt_norm2d_nf90 (ng, tile, iTLM, ncname,         &
     &                                  LBi, UBi, LBj, UBj,             &
     &                                  idTsur(itrc),                   &
     &                              NRM(ifile,ng)%ncid,                 &
     &                              NRM(ifile,ng)%Vid(idTsur(itrc)),    &
     &                              NRM(ifile,ng)%Rindex,               &
#  ifdef MASKING
     &                                  rmask,                          &
#  endif
     &                                  HnormSTF(:,:,itrc))

#  if defined PIO_LIB && defined DISTRIBUTE
                CASE (io_pio)
                  IF (NRM(ifile,ng)%pioVar(idTsur(itrc))%dkind.eq.      &
     &                PIO_double) THEN
                    ioDesc => ioDesc_dp_r2dvar(ng)
                  ELSE
                    ioDesc => ioDesc_sp_r2dvar(ng)
                  END IF
                  CALL wrt_norm2d_pio (ng, tile, iTLM, ncname,          &
     &                                 LBi, UBi, LBj, UBj,              &
     &                                 idTsur(itrc),                    &
     &                              NRM(ifile,ng)%pioFile,              &
     &                              NRM(ifile,ng)%pioVar(idTsur(itrc)), &
     &                              NRM(ifile,ng)%Rindex,               &
     &                                 ioDesc,                          &
#   ifdef MASKING
     &                                 rmask,                           &
#   endif
     &                                 HnormSTF(:,:,itrc))
#  endif
              END SELECT
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END IF
          END IF
        END DO
# endif
      END IF
#endif
!
      IF (Master) THEN
        WRITE (stdout,30)
      END IF

 10   FORMAT (/,' Error Covariance Normalization Factors: ',            &
     &        'Exact Method',/)
 20   FORMAT (4x,'Computing',1x,a,1x,a)
 30   FORMAT (/)
!
      RETURN
      END SUBROUTINE normalization_tile

!
!***********************************************************************
      SUBROUTINE randomization_tile (ng, tile,                          &
     &                               LBi, UBi, LBj, UBj, LBk, UBk,      &
     &                               LBij, UBij,                        &
     &                               IminS, ImaxS, JminS, JmaxS,        &
     &                               nstp, nnew, ifac,                  &
#ifdef MULTI_SCALE_B
     &                               MiterCI, Clap,                     &
     &                               emaxr2d, eminr2d,                  &
     &                               emaxu2d, eminu2d,                  &
     &                               emaxv2d, eminv2d,                  &
     &                               emaxr3d, eminr3d,                  &
     &                               emaxu3d, eminu3d,                  &
     &                               emaxv3d, eminv3d,                  &
# ifdef ADJUST_WSTRESS
     &                               emaxu2ds, eminu2ds,                &
     &                               emaxv2ds, eminv2ds,                &
# endif
# ifdef ADJUST_STFLUX
     &                               emaxr2ds, eminr2ds,                &
# endif
# ifdef ADJUST_BOUNDARY
     &                               emaxr1d, eminr1d,                  &
     &                               emaxu1d, eminu1d,                  &
     &                               emaxv1d, eminv1d,                  &
     &                               emaxr1dz, eminr1dz,                &
     &                               emaxu1dz, eminu1dz,                &
     &                               emaxv1dz, eminv1dz,                &
# endif
     &                               ci_r, ci_p, ci_q, ci_x,            &
     &                               ci_r3d, ci_p3d, ci_q3d, ci_x3d,    &
# ifdef ADJUST_BOUNDARY
     &                               ci_rb, ci_pb, ci_qb, ci_xb,        &
# endif
# ifdef READ_SCALES
     &                               bzsclx, bzscly,                    &
#  ifdef SOLVE3D
     &                               busclx, buscly,                    &
     &                               bvsclx, bvscly,                    &
     &                               btsclx, btscly,                    &
#  endif
     &                               bubsclx, bubscly,                  &
     &                               bvbsclx, bvbscly,                  &
#  ifdef ADJUST_WSTRESS
     &                               bussclx, busscly,                  &
     &                               bvssclx, bvsscly,                  &
#  endif
#  ifdef ADJUST_STFLUX
     &                               btfsclx, btfscly,                  &
#  endif
# endif
#endif
     &                               pm, om_p, om_r, om_u, om_v,        &
     &                               pn, on_p, on_r, on_u, on_v,        &
     &                               pmon_p, pmon_r, pmon_u,            &
     &                               pnom_p, pnom_r, pnom_v,            &
#ifdef MASKING
     &                               pmask, rmask,                      &
     &                               umask, vmask,                      &
#endif
#ifdef SOLVE3D
     &                               h,                                 &
# ifdef ICESHELF
     &                               zice,                              &
# endif
# if defined SEDIMENT && defined SED_MORPH
     &                               bed_thick,                         &
# endif
     &                               Hz, z_r, z_w,                      &
#endif
     &                               Kh,                                &
#ifdef SOLVE3D
     &                               Kv,                                &
#endif
#ifdef ADJUST_BOUNDARY
# ifdef SOLVE3D
     &                               VnormRobc, VnormUobc, VnormVobc,   &
# endif
     &                               HnormRobc, HnormUobc, HnormVobc,   &
#endif
#ifdef ADJUST_WSTRESS
     &                               HnormSUS, HnormSVS,                &
#endif
#if defined ADJUST_STFLUX && defined SOLVE3D
     &                               HnormSTF,                          &
#endif
#ifdef SOLVE3D
     &                               VnormR, VnormU, VnormV,            &
#endif
     &                               HnormR, HnormU, HnormV)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk, LBij, UBij
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp, nnew, ifac
#ifdef MULTI_SCALE_B
      integer, intent(in) :: MiterCI(MAXVAL(Nscale),Ngrids)
      integer, intent(in) :: Clap(MstateVar,MAXVAL(Nscale),Ngrids)

      real(r8), intent(inout) :: emaxr2d(MAXVAL(Mlap), Nscale(ng))
      real(r8), intent(inout) :: eminr2d(MAXVAL(Mlap), Nscale(ng))
      real(r8), intent(inout) :: emaxu2d(MAXVAL(Mlap), Nscale(ng))
      real(r8), intent(inout) :: eminu2d(MAXVAL(Mlap), Nscale(ng))
      real(r8), intent(inout) :: emaxv2d(MAXVAL(Mlap), Nscale(ng))
      real(r8), intent(inout) :: eminv2d(MAXVAL(Mlap), Nscale(ng))
      real(r8), intent(inout) :: emaxr3d(MAXVAL(Mlap), N(ng), NT(ng),   &
     &                                                     Nscale(ng))
      real(r8), intent(inout) :: eminr3d(MAXVAL(Mlap), N(ng), NT(ng),   &
     &                                                     Nscale(ng))
      real(r8), intent(inout) :: emaxu3d(MAXVAL(Mlap), N(ng),Nscale(ng))
      real(r8), intent(inout) :: eminu3d(MAXVAL(Mlap), N(ng),Nscale(ng))
      real(r8), intent(inout) :: emaxv3d(MAXVAL(Mlap), N(ng),Nscale(ng))
      real(r8), intent(inout) :: eminv3d(MAXVAL(Mlap), N(ng),Nscale(ng))
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: emaxu2ds(MAXVAL(Mlap), Nscale(ng))
      real(r8), intent(inout) :: eminu2ds(MAXVAL(Mlap), Nscale(ng))
      real(r8), intent(inout) :: emaxv2ds(MAXVAL(Mlap), Nscale(ng))
      real(r8), intent(inout) :: eminv2ds(MAXVAL(Mlap), Nscale(ng))
# endif
# ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: emaxr2ds(MAXVAL(Mlap), NT(ng),         &
     &                                                    Nscale(ng))
      real(r8), intent(inout) :: eminr2ds(MAXVAL(Mlap), NT(ng),         &
     &                                                    Nscale(ng))
# endif
# ifdef ADJUST_BOUNDARY
      real(r8), intent(inout) :: emaxr1d(MAXVAL(Mlap), 4, Nscale(ng))
      real(r8), intent(inout) :: eminr1d(MAXVAL(Mlap), 4, Nscale(ng))
      real(r8), intent(inout) :: emaxu1d(MAXVAL(Mlap), 4, Nscale(ng))
      real(r8), intent(inout) :: eminu1d(MAXVAL(Mlap), 4, Nscale(ng))
      real(r8), intent(inout) :: emaxv1d(MAXVAL(Mlap), 4, Nscale(ng))
      real(r8), intent(inout) :: eminv1d(MAXVAL(Mlap), 4, Nscale(ng))
      real(r8), intent(inout) :: emaxr1dz(MAXVAL(Mlap), N(ng), NT(ng),  &
     &                                                   4, Nscale(ng))
      real(r8), intent(inout) :: eminr1dz(MAXVAL(Mlap), N(ng), NT(ng),  &
     &                                                   4, Nscale(ng))
      real(r8), intent(inout) :: emaxu1dz(MAXVAL(Mlap), N(ng),          &
     &                                                   4, Nscale(ng))
      real(r8), intent(inout) :: eminu1dz(MAXVAL(Mlap), N(ng),          &
     &                                                   4, Nscale(ng))
      real(r8), intent(inout) :: emaxv1dz(MAXVAL(Mlap), N(ng),          &
     &                                                   4, Nscale(ng))
      real(r8), intent(inout) :: eminv1dz(MAXVAL(Mlap), N(ng),          &
     &                                                   4, Nscale(ng))
# endif
      real(r8), intent(inout) :: ci_r(LBi:,LBj:)
      real(r8), intent(inout) :: ci_p(LBi:,LBj:)
      real(r8), intent(inout) :: ci_q(LBi:,LBj:)
      real(r8), intent(inout) :: ci_x(LBi:,LBj:)
      real(r8), intent(inout) :: ci_r3d(LBi:,LBj:,LBk:)
      real(r8), intent(inout) :: ci_p3d(LBi:,LBj:,LBk:)
      real(r8), intent(inout) :: ci_q3d(LBi:,LBj:,LBk:)
      real(r8), intent(inout) :: ci_x3d(LBi:,LBj:,LBk:)
# ifdef ADJUST_BOUNDARY
      real(r8), intent(inout) :: ci_rb(LBij:)
      real(r8), intent(inout) :: ci_pb(LBij:)
      real(r8), intent(inout) :: ci_qb(LBij:)
      real(r8), intent(inout) :: ci_xb(LBij:)
# endif
# ifdef READ_SCALES
      real(r8), intent(in) :: bzsclx(LBi:,LBj:,:)
      real(r8), intent(in) :: bzscly(LBi:,LBj:,:)
#  ifdef SOLVE3D
      real(r8), intent(in) :: busclx(LBi:,LBj:,:)
      real(r8), intent(in) :: buscly(LBi:,LBj:,:)
      real(r8), intent(in) :: bvsclx(LBi:,LBj:,:)
      real(r8), intent(in) :: bvscly(LBi:,LBj:,:)
      real(r8), intent(in) :: btsclx(LBi:,LBj:,:,:)
      real(r8), intent(in) :: btscly(LBi:,LBj:,:,:)
#  endif
      real(r8), intent(in) :: bubsclx(LBi:,LBj:,:)
      real(r8), intent(in) :: bubscly(LBi:,LBj:,:)
      real(r8), intent(in) :: bvbsclx(LBi:,LBj:,:)
      real(r8), intent(in) :: bvbscly(LBi:,LBj:,:)
#  ifdef ADJUST_WSTRESS
      real(r8), intent(in) :: bussclx(LBi:,LBj:,:)
      real(r8), intent(in) :: busscly(LBi:,LBj:,:)
      real(r8), intent(in) :: bvssclx(LBi:,LBj:,:)
      real(r8), intent(in) :: bvsscly(LBi:,LBj:,:)
#  endif
#  ifdef ADJUST_STFLUX
      real(r8), intent(in) :: btfsclx(LBi:,LBj:,:,:)
      real(r8), intent(in) :: btfscly(LBi:,LBj:,:,:)
#  endif
# endif
#endif
!
#ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: om_p(LBi:,LBj:)
      real(r8), intent(in) :: om_r(LBi:,LBj:)
      real(r8), intent(in) :: om_u(LBi:,LBj:)
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: on_p(LBi:,LBj:)
      real(r8), intent(in) :: on_r(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)
      real(r8), intent(in) :: on_v(LBi:,LBj:)
      real(r8), intent(in) :: pmon_p(LBi:,LBj:)
      real(r8), intent(in) :: pmon_r(LBi:,LBj:)
      real(r8), intent(in) :: pmon_u(LBi:,LBj:)
      real(r8), intent(in) :: pnom_p(LBi:,LBj:)
      real(r8), intent(in) :: pnom_r(LBi:,LBj:)
      real(r8), intent(in) :: pnom_v(LBi:,LBj:)
# ifdef MASKING
      real(r8), intent(in) :: pmask(LBi:,LBj:)
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
# endif
      real(r8), intent(in) :: Kh(LBi:,LBj:)
# ifdef SOLVE3D
      real(r8), intent(in) :: Kv(LBi:,LBj:,0:)
#  ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:,LBj:)
#  endif
#  if defined SEDIMENT && defined SED_MORPH
      real(r8), intent(in):: bed_thick(LBi:,LBj:,:)
#  endif
      real(r8), intent(inout) :: h(LBi:,LBj:)
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
      real(r8), intent(out) :: VnormRobc(LBij:,:,:,:)
      real(r8), intent(out) :: VnormUobc(LBij:,:,:)
      real(r8), intent(out) :: VnormVobc(LBij:,:,:)
#  endif
      real(r8), intent(out) :: HnormRobc(LBij:,:)
      real(r8), intent(out) :: HnormUobc(LBij:,:)
      real(r8), intent(out) :: HnormVobc(LBij:,:)
# endif
# ifdef ADJUST_WSTRESS
      real(r8), intent(out) :: HnormSUS(LBi:,LBj:)
      real(r8), intent(out) :: HnormSVS(LBi:,LBj:)
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(out) :: HnormSTF(LBi:,LBj:,:)
# endif
# ifdef SOLVE3D
      real(r8), intent(out) :: VnormR(LBi:,LBj:,:,:,:)
      real(r8), intent(out) :: VnormU(LBi:,LBj:,:,:)
      real(r8), intent(out) :: VnormV(LBi:,LBj:,:,:)
# endif
      real(r8), intent(out) :: HnormR(LBi:,LBj:,:)
      real(r8), intent(out) :: HnormU(LBi:,LBj:,:)
      real(r8), intent(out) :: HnormV(LBi:,LBj:,:)
# ifdef SOLVE3D
      real(r8), intent(out) :: Hz(LBi:,LBj:,:)
      real(r8), intent(out) :: z_r(LBi:,LBj:,:)
      real(r8), intent(out) :: z_w(LBi:,LBj:,0:)
# endif
#else
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_p(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_r(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_p(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_r(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pmon_p(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pmon_r(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pmon_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pnom_p(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pnom_r(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pnom_v(LBi:UBi,LBj:UBj)
# ifdef MASKING
      real(r8), intent(in) :: pmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in) :: Kh(LBi:UBi,LBj:UBj)
# ifdef SOLVE3D
      real(r8), intent(in) :: Kv(LBi:UBi,LBj:UBj,0:N(ng))
#  ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:UBi,LBj:UBj)
#  endif
#  if defined SEDIMENT && defined SED_MORPH
      real(r8), intent(in):: bed_thick(LBi:UBi,LBj:UBj,3)
#  endif
      real(r8), intent(inout) :: h(LBi:UBi,LBj:UBj)
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
      real(r8), intent(out) :: VnormRobc(LBij:UBij,N(ng),4,NT(ng))
      real(r8), intent(out) :: VnormUobc(LBij:UBij,N(ng),4)
      real(r8), intent(out) :: VnormVobc(LBij:UBij,N(ng),4)
#  endif
      real(r8), intent(out) :: HnormRobc(LBij:UBij,4)
      real(r8), intent(out) :: HnormUobc(LBij:UBij,4)
      real(r8), intent(out) :: HnormVobc(LBij:UBij,4)
# endif
# ifdef ADJUST_WSTRESS
      real(r8), intent(out) :: HnormSUS(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: HnormSVS(LBi:UBi,LBj:UBj)
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(out) :: HnormSTF(LBi:UBi,LBj:UBj,NT(ng))
# endif
# ifdef SOLVE3D
      real(r8), intent(out) :: VnormR(LBi:UBi,LBj:UBj,N(ng),NSA,NT(ng))
      real(r8), intent(out) :: VnormU(LBi:UBi,LBj:UBj,N(ng),NSA)
      real(r8), intent(out) :: VnormV(LBi:UBi,LBj:UBj,N(ng),NSA)
# endif
      real(r8), intent(out) :: HnormR(LBi:UBi,LBj:UBj,NSA)
      real(r8), intent(out) :: HnormU(LBi:UBi,LBj:UBj,NSA)
      real(r8), intent(out) :: HnormV(LBi:UBi,LBj:UBj,NSA)
# ifdef SOLVE3D
      real(r8), intent(out) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(out) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(out) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
# endif
#endif
!
!  Local variable declarations.
!
#ifdef SOLVE3D
      logical :: Ldiffer, Lsame
#endif
#ifdef ADJUST_BOUNDARY
      logical :: Lconvolve(4)
#endif
#ifdef MULTI_SCALE_B
      logical :: Lweak
#endif
!
      integer :: i, ifile, is, iter, j, rec
      integer :: iLap
#ifdef MULTI_SCALE_B
      integer :: it, klev, ms, ns
#endif
#ifdef SOLVE3D
      integer :: UBt, itrc, k
#endif
#ifdef ADJUST_BOUNDARY
      integer :: IJlen, IJKlen, ib, ibry, ic, ifield
#endif
      integer :: start(4), total(4)
!
      real(dp) :: my_time
      real(r8) :: Aavg, Amax, Amin, Asqr, FacAvg, FacSqr
      real(r8) :: cff, val

#ifdef ADJUST_BOUNDARY
      real(r8) :: Bavg, Bmin, Bmax, Bsqr

      real(r8), parameter :: Aspv = 0.0_r8
#endif
!
      real(r8), dimension(LBi:UBi,LBj:UBj) :: A2d
      real(r8), dimension(LBi:UBi,LBj:UBj) :: A2davg
      real(r8), dimension(LBi:UBi,LBj:UBj) :: A2dsqr
      real(r8), dimension(LBi:UBi,LBj:UBj) :: Hscale
#ifdef ADJUST_BOUNDARY
      real(r8), dimension(LBij:UBij) :: B2d
      real(r8), dimension(LBij:UBij) :: B2davg
      real(r8), dimension(LBij:UBij) :: B2dsqr
      real(r8), dimension(LBij:UBij) :: HscaleB
#endif
#ifdef SOLVE3D
      real(r8), dimension(LBi:UBi,LBj:UBj,1:N(ng)) :: A3d
      real(r8), dimension(LBi:UBi,LBj:UBj,1:N(ng)) :: A3davg
      real(r8), dimension(LBi:UBi,LBj:UBj,1:N(ng)) :: A3dsqr
      real(r8), dimension(LBi:UBi,LBj:UBj,1:N(ng)) :: Vscale
# ifdef ADJUST_BOUNDARY
      real(r8), dimension(LBij:UBij,1:N(ng)) :: B3d
      real(r8), dimension(LBij:UBij,1:N(ng)) :: B3davg
      real(r8), dimension(LBij:UBij,1:N(ng)) :: B3dsqr
      real(r8), dimension(LBij:UBij,1:N(ng)) :: VscaleB
#  ifdef DISTRIBUTE
      real(r8), dimension((UBij-LBij+1)*N(ng)) :: Bwrk
#  endif
# endif
#endif
!
      character (len=40 ) :: Text
      character (len=256) :: ncname

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", randomization_tile"

#if defined PIO_LIB && defined DISTRIBUTE
!
      TYPE (IO_Desc_t),  pointer :: ioDesc
#endif

#include "set_bounds.h"
!
      SourceFile=MyFile

      my_time=tdays(ng)*day2sec

#ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Compute time invariant depths (use zero free-surface).
!-----------------------------------------------------------------------
!
      DO i=LBi,UBi
        DO j=LBj,UBj
          A2d(i,j)=0.0_r8
        END DO
      END DO

      CALL set_depth_tile (ng, tile, iNLM,                              &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     nstp, nnew,                                  &
     &                     h,                                           &
# ifdef ICESHELF
     &                     zice,                                        &
# endif
# if defined SEDIMENT && defined SED_MORPH
     &                     bed_thick,                                   &
# endif
     &                     A2d,                                         &
     &                     Hz, z_r, z_w)
#endif
!
!-----------------------------------------------------------------------
!  Compute initial conditions and model erro covariance, B,
!  normalization factors using the randomization approach of Fisher
!  and Courtier (1995). These factors ensure that the diagonal
!  elements of B are equal to unity. Notice that in applications
!  with land/sea masking, the boundary conditions will produce
!  large changes in the covariance structures near the boundary.
!
!  Initialize factors with randon numbers ("white-noise") having an
!  uniform distribution (zero mean and unity variance). Then, scale
!  by the inverse squared root area (2D) or volume (3D) and "color"
!  with the diffusion operator. Iterate this step over a specified
!  number of ensamble members, Nrandom.
!-----------------------------------------------------------------------
!
      IF (Master) WRITE (stdout,10)

      FILE_LOOP : DO ifile=1,NSA

        IF (LwrtNRM(ifile,ng)) THEN
          Lweak=.FALSE.
          IF (ifile.eq.1) THEN
            Text='initial conditions'
          ELSE IF (ifile.eq.2) THEN
            Text='model'
            Lweak=.TRUE.
          END IF
!
!  Set randomization summation factors.
!
          FacAvg=1.0_r8/REAL(Nrandom,r8)
          FacSqr=SQRT(REAL(Nrandom,r8))
!
!  Set time record index to write in normalization NetCDF file.
!
          ncname=NRM(ifile,ng)%name
          NRM(ifile,ng)%Rindex=NRM(ifile,ng)%Rindex+1
          NRM(ifile,ng)%Nrec=NRM(ifile,ng)%Nrec+1
!
!  Write out model time (s).
!
          SELECT CASE (NRM(ifile,ng)%IOtype)
            CASE (io_nf90)
              CALL netcdf_put_fvar (ng, iTLM, ncname,                   &
     &                              Vname(1,idtime), my_time,           &
     &                         start = (/NRM(ifile,ng)%Rindex/),        &
     &                         total = (/1/),                           &
     &                         ncid = NRM(ifile,ng)%ncid,               &
     &                         varid = NRM(ifile,ng)%Vid(idtime))
#if defined PIO_LIB && defined DISTRIBUTE
            CASE (io_pio)
              CALL pio_netcdf_put_fvar (ng, iTLM, ncname,               &
     &                                  Vname(1,idtime), my_time,       &
     &                         start = (/NRM(ifile,ng)%Rindex/),        &
     &                         total = (/1/),                           &
     &                         pioFile = NRM(ifile,ng)%pioFile,         &
     &                         pioVar = NRM(ifile,ng)%pioVar(idtime)%vd)
#endif
          END SELECT
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  2D norm at RHO-points.
!
          IF (Cnorm(ifile,isFsur)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Text),                             &
     &                          '2D normalization factors at RHO-points'
              FLUSH (stdout)
            END IF
            DO ns=1,Nscale(ng)
              ms=ns
              iLap=Clap(isFsur,ns,ng)
              DO j=JstrT,JendT
                DO i=IstrT,IendT
                  A2davg(i,j)=0.0_r8
                  A2dsqr(i,j)=0.0_r8
                  Hscale(i,j)=1.0_r8/SQRT(om_r(i,j)*on_r(i,j))
                END DO
              END DO
              DO iter=1,Nrandom
                CALL white_noise2d (ng, iTLM, r2dvar, Rscheme(ng),      &
     &                            IstrR, IendR, JstrR, JendR,           &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Amin, Amax, A2d)
                DO j=JstrT,JendT
                  DO i=IstrT,IendT
                    A2d(i,j)=A2d(i,j)*Hscale(i,j)
                  END DO
                END DO
!
                CALL dabc_r2d_tile (ng, tile,                           &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       A2d)
!
                CALL tl_CI_tile (ng, tile,                              &
     &                       LBi, UBi, LBj, UBj, LBij, UBij,            &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       nstp, nnew, Lweak, ifac, isFsur,           &
     &                       ms, MiterCI(ns,ng), iLap, r2dvar, 1, 1,    &
     &                       emaxr2d(:,ns), eminr2d(:,ns),              &
     &                       ci_r, ci_p, ci_q, ci_x,                    &
# ifdef READ_SCALES
     &                       bzsclx(:,:,ns), bzscly(:,:,ns),            &
# endif
# ifdef SOLVE3D
     &                       Hz,                                        &
# endif
     &                       A2d)
                DO j=Jstr,Jend
                  DO i=Istr,Iend
                    A2davg(i,j)=A2davg(i,j)+A2d(i,j)
                    A2dsqr(i,j)=A2dsqr(i,j)+A2d(i,j)*A2d(i,j)
                  END DO
                END DO
              END DO
              DO j=Jstr,Jend
                DO i=Istr,Iend
                  Aavg=FacAvg*A2davg(i,j)
                  Asqr=FacAvg*A2dsqr(i,j)
# ifdef MASKING
                  IF (rmask(i,j).gt.0.0_r8) THEN
                    HnormR(i,j,ifile)=HnormR(i,j,ifile)+                &
     &                                  Bwgt(isFsur,ms,ng)/SQRT(Asqr)
                  ELSE
                    HnormR(i,j,ifile)=0.0_r8
                  END IF
# else
                  HnormR(i,j,ifile)=HnormR(i,j,ifile)+                  &
     &                                  Bwgt(isFsur,ms,ng)/SQRT(Asqr)
# endif
                END DO
              END DO
            END DO  ! End of ns loop.
            CALL dabc_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          HnormR(:,:,ifile))
#ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iTLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          HnormR(:,:,ifile))
#endif
!
            SELECT CASE (NRM(ifile,ng)%IOtype)
              CASE (io_nf90)
                CALL wrt_norm2d_nf90 (ng, tile, iTLM, ncname,           &
     &                                LBi, UBi, LBj, UBj, idFsur,       &
     &                                NRM(ifile,ng)%ncid,               &
     &                                NRM(ifile,ng)%Vid(idFsur),        &
     &                                NRM(ifile,ng)%Rindex,             &
#ifdef MASKING
     &                                rmask,                            &
#endif
     &                                HnormR(:,:,ifile))

#if defined PIO_LIB && defined DISTRIBUTE
              CASE (io_pio)
                IF (NRM(ifile,ng)%pioVar(idFsur)%dkind.eq.              &
     &              PIO_double) THEN
                  ioDesc => ioDesc_dp_r2dvar(ng)
                ELSE
                  ioDesc => ioDesc_sp_r2dvar(ng)
                END IF
                CALL wrt_norm2d_pio (ng, tile, iTLM, ncname,            &
     &                               LBi, UBi, LBj, UBj, idFsur,        &
     &                               NRM(ifile,ng)%pioFile,             &
     &                               NRM(ifile,ng)%pioVar(idFsur),      &
     &                               NRM(ifile,ng)%Rindex,              &
     &                               ioDesc,                            &
# ifdef MASKING
     &                               rmask,                             &
# endif
     &                               HnormR(:,:,ifile))
#endif
            END SELECT
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
!
!  2D norm at U-points.
!
          IF (Cnorm(ifile,isUbar)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Text),                             &
     &                          '2D normalization factors at   U-points'
              FLUSH (stdout)
            END IF
            DO ns=1,Nscale(ng)
              ms=ns
              iLap=Clap(isUbar,ns,ng)
              DO j=JstrT,JendT
                DO i=IstrP,IendT
                  A2davg(i,j)=0.0_r8
                  A2dsqr(i,j)=0.0_r8
                  Hscale(i,j)=1.0_r8/SQRT(om_u(i,j)*on_u(i,j))
                END DO
              END DO
              DO iter=1,Nrandom
                CALL white_noise2d (ng, iTLM, u2dvar, Rscheme(ng),      &
     &                            Istr, IendR, JstrR, JendR,            &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Amin, Amax, A2d)
                DO j=JstrT,JendT
                  DO i=IstrP,IendT
                    A2d(i,j)=A2d(i,j)*Hscale(i,j)
                  END DO
                END DO
!
                CALL dabc_u2d_tile (ng, tile,                           &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       A2d)
!
                CALL tl_CI_tile (ng, tile,                              &
     &                       LBi, UBi, LBj, UBj, LBij, UBij,            &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       nstp, nnew, Lweak, ifac, isUbar,           &
     &                       ms, MiterCI(ns,ng), iLap, u2dvar, 1, 1,    &
     &                       emaxu2d(:,ns), eminu2d(:,ns),              &
     &                       ci_r, ci_p, ci_q, ci_x,                    &
# ifdef READ_SCALES
     &                       bubsclx(:,:,ns), bubscly(:,:,ns),          &
# endif
# ifdef SOLVE3D
     &                       Hz,                                        &
# endif
     &                       A2d)
                DO j=Jstr,Jend
                  DO i=IstrU,Iend
                    A2davg(i,j)=A2davg(i,j)+A2d(i,j)
                    A2dsqr(i,j)=A2dsqr(i,j)+A2d(i,j)*A2d(i,j)
                  END DO
                END DO
              END DO
              DO j=Jstr,Jend
                DO i=IstrU,Iend
                  Aavg=FacAvg*A2davg(i,j)
                  Asqr=FacAvg*A2dsqr(i,j)
# ifdef MASKING
                  IF (umask(i,j).gt.0.0_r8) THEN
                    HnormU(i,j,ifile)=HnormU(i,j,ifile)+                &
     &                                 Bwgt(isUbar,ms,ng)/SQRT(Asqr)
                  ELSE
                    HnormU(i,j,ifile)=0.0_r8
                  END IF
# else
                  HnormU(i,j,ifile)=HnormU(i,j,ifile)+                  &
     &                                 Bwgt(isUbar,ms,ng)/SQRT(Asqr)
# endif
                END DO
              END DO
            END DO   ! End of ns loop
            CALL dabc_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          HnormU(:,:,ifile))
#ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iTLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          HnormU(:,:,ifile))
#endif
!
            SELECT CASE (NRM(ifile,ng)%IOtype)
              CASE (io_nf90)
                CALL wrt_norm2d_nf90 (ng, tile, iTLM, ncname,           &
     &                                LBi, UBi, LBj, UBj, idUbar,       &
     &                                NRM(ifile,ng)%ncid,               &
     &                                NRM(ifile,ng)%Vid(idUbar),        &
     &                                NRM(ifile,ng)%Rindex,             &
#ifdef MASKING
     &                                umask,                            &
#endif
     &                                HnormU(:,:,ifile))

#if defined PIO_LIB && defined DISTRIBUTE
              CASE (io_pio)
                IF (NRM(ifile,ng)%pioVar(idUbar)%dkind.eq.              &
     &              PIO_double) THEN
                  ioDesc => ioDesc_dp_u2dvar(ng)
                ELSE
                  ioDesc => ioDesc_sp_u2dvar(ng)
                END IF
                CALL wrt_norm2d_pio (ng, tile, iTLM, ncname,            &
     &                               LBi, UBi, LBj, UBj, idUbar,        &
     &                               NRM(ifile,ng)%pioFile,             &
     &                               NRM(ifile,ng)%pioVar(idUbar),      &
     &                               NRM(ifile,ng)%Rindex,              &
     &                               ioDesc,                            &
# ifdef MASKING
     &                               umask,                             &
# endif
     &                               HnormU(:,:,ifile))
#endif
            END SELECT
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
!
!  2D norm at V-points.
!
          IF (Cnorm(ifile,isVbar)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Text),                             &
     &                          '2D normalization factors at   V-points'
              FLUSH (stdout)
            END IF
            DO ns=1,Nscale(ng)
              ms=ns
              iLap=Clap(isVbar,ns,ng)
              DO j=JstrP,JendT
                DO i=IstrT,IendT
                  A2davg(i,j)=0.0_r8
                  A2dsqr(i,j)=0.0_r8
                  Hscale(i,j)=1.0_r8/SQRT(om_v(i,j)*on_v(i,j))
                END DO
              END DO
              DO iter=1,Nrandom
                CALL white_noise2d (ng, iTLM, v2dvar, Rscheme(ng),      &
     &                            IstrR, IendR, Jstr, JendR,            &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Amin, Amax, A2d)
                DO j=JstrP,JendT
                  DO i=IstrT,IendT
                    A2d(i,j)=A2d(i,j)*Hscale(i,j)
                  END DO
                END DO
!
                CALL dabc_v2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              A2d)
!
                CALL tl_CI_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj, LBij, UBij,        &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           nstp, nnew, Lweak, ifac, isVbar,       &
     &                           ms, MiterCI(ns,ng), iLap,              &
     &                           v2dvar, 1, 1,                          &
     &                           emaxv2d(:,ns), eminv2d(:,ns),          &
     &                           ci_r, ci_p, ci_q, ci_x,                &
# ifdef READ_SCALES
     &                           bvbsclx(:,:,ns), bvbscly(:,:,ns),      &
# endif
# ifdef SOLVE3D
     &                           Hz,                                    &
# endif
     &                           A2d)
                DO j=JstrV,Jend
                  DO i=Istr,Iend
                    A2davg(i,j)=A2davg(i,j)+A2d(i,j)
                    A2dsqr(i,j)=A2dsqr(i,j)+A2d(i,j)*A2d(i,j)
                  END DO
                END DO
              END DO
              DO j=JstrV,Jend
                DO i=Istr,Iend
                  Aavg=FacAvg*A2davg(i,j)
                  Asqr=FacAvg*A2dsqr(i,j)
# ifdef MASKING
                  IF (vmask(i,j).gt.0.0_r8) THEN
                    HnormV(i,j,ifile)=HnormV(i,j,ifile)+                &
     &                                Bwgt(isVbar,ms,ng)/SQRT(Asqr)
                  ELSE
                    HnormV(i,j,ifile)=0.0_r8
                  END IF
# else
                  HnormV(i,j,ifile)=HnormV(i,j,ifile)+                  &
     &                              Bwgt(isVbar,ms,ng)/SQRT(Asqr)
# endif
                END DO
              END DO
            END DO    ! End of ns loop
            CALL dabc_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          HnormV(:,:,ifile))
#ifdef DISTRIBUTE
            CALL mp_exchange2d (ng, tile, iTLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          HnormV(:,:,ifile))
#endif
!
            SELECT CASE (NRM(ifile,ng)%IOtype)
              CASE (io_nf90)
                CALL wrt_norm2d_nf90 (ng, tile, iTLM, ncname,           &
     &                                LBi, UBi, LBj, UBj, idVbar,       &
     &                                NRM(ifile,ng)%ncid,               &
     &                                NRM(ifile,ng)%Vid(idVbar),        &
     &                                NRM(ifile,ng)%Rindex,             &
#ifdef MASKING
     &                                vmask,                            &
#endif
     &                                HnormV(:,:,ifile))

#if defined PIO_LIB && defined DISTRIBUTE
              CASE (io_pio)
                IF (NRM(ifile,ng)%pioVar(idVbar)%dkind.eq.              &
     &              PIO_double) THEN
                  ioDesc => ioDesc_dp_v2dvar(ng)
                ELSE
                  ioDesc => ioDesc_sp_v2dvar(ng)
                END IF
                CALL wrt_norm2d_pio (ng, tile, iTLM, ncname,            &
     &                               LBi, UBi, LBj, UBj, idVbar,        &
     &                               NRM(ifile,ng)%pioFile,             &
     &                               NRM(ifile,ng)%pioVar(idVbar),      &
     &                               NRM(ifile,ng)%Rindex,              &
     &                               ioDesc,                            &
# ifdef MASKING
     &                               vmask,                             &
# endif
     &                               HnormV(:,:,ifile))
#endif
            END SELECT
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF

#ifdef SOLVE3D
!
!  3D norm U-points.
!
          IF (Cnorm(ifile,isUvel)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Text),                             &
     &                          '3D normalization factors at   U-points'
              FLUSH (stdout)
            END IF
            DO ns=1,Nscale(ng)
              ms=ns
              iLap=Clap(isUvel,ns,ng)
              DO j=JstrT,JendT
                DO i=IstrP,IendT
                  val=om_u(i,j)*on_u(i,j)*0.5_r8
                  DO k=1,N(ng)
                    A3davg(i,j,k)=0.0_r8
                    A3dsqr(i,j,k)=0.0_r8
                    Vscale(i,j,k)=1.0_r8/SQRT(val*(Hz(i-1,j,k)+Hz(i,j,k)))
                  END DO
                END DO
              END DO
              DO iter=1,Nrandom
                CALL white_noise3d (ng, iTLM, u3dvar, Rscheme(ng),      &
     &                              Istr, IendR, JstrR, JendR,          &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              Amin, Amax, A3d)
                DO k=1,N(ng)
                  DO j=JstrT,JendT
                    DO i=IstrP,IendT
                      A3d(i,j,k)=A3d(i,j,k)*Vscale(i,j,k)
                    END DO
                  END DO
                END DO
                DO k=1,N(ng)
                  klev=k
                  CALL dabc_u2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                A3d(:,:,klev))
                END DO
!
                CALL tl_CI_3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              LBij, UBij,                         &
     &                              IminS, ImaxS, JminS, JmaxS,         &
     &                              Lweak, ifac, isUvel,                &
     &                              ms, MiterCI(ns,ng), iLap,           &
     &                              u3dvar, 1,                          &
     &                              emaxu3d(:,:,ns), eminu3d(:,:,ns),   &
     &                              ci_r3d, ci_p3d, ci_q3d, ci_x3d,     &
#  ifdef READ_SCALES
     &                              busclx(:,:,ns), buscly(:,:,ns),     &
#  endif
     &                              Hz, z_r, A3d)
!
                CALL tl_conv_u3d_tile (ng, tile, iTLM,                  &
     &                                 LBi, UBi, LBj, UBj, 1, N(ng),    &
     &                                 IminS, ImaxS, JminS, JmaxS,      &
     &                                 NghostPoints,                    &
     &                                 NHsteps(ifile,isUvel)/ifac,      &
     &                                 NVsteps(ifile,isUvel)/ifac,      &
     &                                 DTsizeH(ifile,isUvel),           &
     &                                 DTsizeV(ifile,isUvel),           &
     &                                 Kh, Kv,                          &
     &                                 pm, pn,                          &
#  ifdef GEOPOTENTIAL_HCONV
     &                                 on_r, om_p,                      &
#  else
     &                                 pmon_r, pnom_p,                  &
#  endif
#  ifdef MASKING
#   ifdef GEOPOTENTIAL_HCONV
     &                                 pmask, rmask, umask, vmask,      &
#   else
     &                                 umask, pmask,                    &
#   endif
#  endif
     &                                 Hz, z_r,                         &
     &                                 A3d)
                DO k=1,N(ng)
                  DO j=Jstr,Jend
                    DO i=IstrU,Iend
                      A3davg(i,j,k)=A3davg(i,j,k)+A3d(i,j,k)
                      A3dsqr(i,j,k)=A3dsqr(i,j,k)+A3d(i,j,k)*A3d(i,j,k)
                    END DO
                  END DO
                END DO
              END DO
              DO k=1,N(ng)
                DO j=Jstr,Jend
                  DO i=IstrU,Iend
                    Aavg=FacAvg*A3davg(i,j,k)
                    Asqr=FacAvg*A3dsqr(i,j,k)
#  ifdef MASKING
                    IF (umask(i,j).gt.0.0_r8) THEN
                      VnormU(i,j,k,ifile)=VnormU(i,j,k,ifile)+          &
     &                                    Bwgt(isUvel,ms,ng)/SQRT(Asqr)
                    ELSE
                      VnormU(i,j,k,ifile)=0.0_r8
                    END IF
#  else
                    VnormU(i,j,k,ifile)=VnormU(i,j,k,ifile)+            &
     &                                  Bwgt(isUvel,ms,ng)/SQRT(Asqr)
#  endif
                  END DO
                END DO
              END DO
            END DO
            CALL dabc_u3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          VnormU(:,:,:,ifile))
# ifdef DISTRIBUTE
            CALL mp_exchange3d (ng, tile, iTLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          VnormU(:,:,:,ifile))
# endif
!
            SELECT CASE (NRM(ifile,ng)%IOtype)
              CASE (io_nf90)
                CALL wrt_norm3d_nf90 (ng, tile, iTLM, ncname,           &
     &                                LBi, UBi, LBj, UBj, 1, N(ng),     &
     &                                idUvel, NRM(ifile,ng)%ncid,       &
     &                                NRM(ifile,ng)%Vid(idUvel),        &
     &                                NRM(ifile,ng)%Rindex,             &
# ifdef MASKING
     &                                umask,                            &
# endif
     &                                VnormU(:,:,:,ifile))

# if defined PIO_LIB && defined DISTRIBUTE
              CASE (io_pio)
                IF (NRM(ifile,ng)%pioVar(idUvel)%dkind.eq.              &
     &              PIO_double) THEN
                  ioDesc => ioDesc_dp_u3dvar(ng)
                ELSE
                  ioDesc => ioDesc_sp_u3dvar(ng)
                END IF
                CALL wrt_norm3d_pio (ng, tile, iTLM, ncname,            &
     &                               LBi, UBi, LBj, UBj, 1, N(ng),      &
     &                               idUvel, NRM(ifile,ng)%pioFile,     &
     &                               NRM(ifile,ng)%pioVar(idUvel),      &
     &                               NRM(ifile,ng)%Rindex,              &
     &                               ioDesc,                            &
#  ifdef MASKING
     &                               umask,                             &
#  endif
     &                               VnormU(:,:,:,ifile))
# endif
            END SELECT
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
!
!  3D norm at V-points.
!
          IF (Cnorm(ifile,isVvel)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Text),                             &
     &                          '3D normalization factors at   V-points'
              FLUSH (stdout)
            END IF
            DO ns=1,Nscale(ng)
              ms=ns
              iLap=Clap(isVvel,ns,ng)
              DO j=JstrP,JendT
                DO i=IstrT,IendT
                  val=om_v(i,j)*on_v(i,j)*0.5_r8
                  DO k=1,N(ng)
                    A3davg(i,j,k)=0.0_r8
                    A3dsqr(i,j,k)=0.0_r8
                    Vscale(i,j,k)=1.0_r8/                               &
     &                            SQRT(val*(Hz(i,j-1,k)+Hz(i,j,k)))
                  END DO
                END DO
              END DO
              DO iter=1,Nrandom
                CALL white_noise3d (ng, iTLM, v3dvar, Rscheme(ng),      &
     &                              IstrR, IendR, Jstr, JendR,          &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              Amin, Amax, A3d)
                DO k=1,N(ng)
                  DO j=JstrP,JendT
                    DO i=IstrT,IendT
                      A3d(i,j,k)=A3d(i,j,k)*Vscale(i,j,k)
                    END DO
                  END DO
                END DO
                DO k=1,N(ng)
                  klev=k
                  CALL dabc_v2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                A3d(:,:,klev))
                END DO
!
                CALL tl_CI_3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              LBij, UBij,                         &
     &                              IminS, ImaxS, JminS, JmaxS,         &
     &                              Lweak, ifac, isVvel,                &
     &                              ms, MiterCI(ns,ng), iLap,           &
     &                              v3dvar, 1,                          &
     &                              emaxv3d(:,:,ns), eminv3d(:,:,ns),   &
     &                              ci_r3d, ci_p3d, ci_q3d, ci_x3d,     &
#ifdef READ_SCALES
     &                              bvsclx(:,:,ns), bvscly(:,:,ns),     &
#endif
     &                              Hz, z_r, A3d)

                CALL tl_conv_v3d_tile (ng, tile, iTLM,                  &
     &                                 LBi, UBi, LBj, UBj, 1, N(ng),    &
     &                                 IminS, ImaxS, JminS, JmaxS,      &
     &                                 NghostPoints,                    &
     &                                 NHsteps(ifile,isVvel)/ifac,      &
     &                                 NVsteps(ifile,isVvel)/ifac,      &
     &                                 DTsizeH(ifile,isVvel),           &
     &                                 DTsizeV(ifile,isVvel),           &
     &                                 Kh, Kv,                          &
     &                                 pm, pn,                          &
#  ifdef GEOPOTENTIAL_HCONV
     &                                 on_p, om_r,                      &
#  else
     &                                 pmon_p, pnom_r,                  &
#  endif
#  ifdef MASKING
#   ifdef GEOPOTENTIAL_HCONV
     &                                 pmask, rmask, umask, vmask,      &
#   else
     &                                 vmask, pmask,                    &
#   endif
#  endif
     &                                 Hz, z_r,                         &
     &                                 A3d)
                DO k=1,N(ng)
                  DO j=JstrV,Jend
                    DO i=Istr,Iend
                      A3davg(i,j,k)=A3davg(i,j,k)+A3d(i,j,k)
                      A3dsqr(i,j,k)=A3dsqr(i,j,k)+A3d(i,j,k)*A3d(i,j,k)
                    END DO
                  END DO
                END DO
              END DO
              DO k=1,N(ng)
                DO j=JstrV,Jend
                  DO i=Istr,Iend
                    Aavg=FacAvg*A3davg(i,j,k)
                    Asqr=FacAvg*A3dsqr(i,j,k)
#  ifdef MASKING
                    IF (vmask(i,j).gt.0.0_r8) THEN
                      VnormV(i,j,k,ifile)=VnormV(i,j,k,ifile)+          &
     &                                    Bwgt(isVvel,ms,ng)/SQRT(Asqr)
                    ELSE
                      VnormV(i,j,k,ifile)=0.0_r8
                    END IF
#  else
                    VnormV(i,j,k,ifile)=VnormV(i,j,k,ifile)+            &
     &                                  Bwgt(isVvel,ms,ng)/SQRT(Asqr)
#  endif
                  END DO
                END DO
              END DO
            END DO
            CALL dabc_v3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          VnormV(:,:,:,ifile))
# ifdef DISTRIBUTE
            CALL mp_exchange3d (ng, tile, iTLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          VnormV(:,:,:,ifile))
# endif
!
            SELECT CASE (NRM(ifile,ng)%IOtype)
              CASE (io_nf90)
                CALL wrt_norm3d_nf90 (ng, tile, iTLM, ncname,           &
     &                                LBi, UBi, LBj, UBj, 1, N(ng),     &
     &                                idVvel, NRM(ifile,ng)%ncid,       &
     &                                NRM(ifile,ng)%Vid(idVvel),        &
     &                                NRM(ifile,ng)%Rindex,             &
# ifdef MASKING
     &                                vmask,                            &
# endif
     &                                VnormV(:,:,:,ifile))

# if defined PIO_LIB && defined DISTRIBUTE
              CASE (io_pio)
                IF (NRM(ifile,ng)%pioVar(idVvel)%dkind.eq.              &
     &              PIO_double) THEN
                  ioDesc => ioDesc_dp_v3dvar(ng)
                ELSE
                  ioDesc => ioDesc_sp_v3dvar(ng)
                END IF
                CALL wrt_norm3d_pio (ng, tile, iTLM, ncname,            &
     &                               LBi, UBi, LBj, UBj, 1, N(ng),      &
     &                               idVvel, NRM(ifile,ng)%pioFile,     &
     &                               NRM(ifile,ng)%pioVar(idVvel),      &
     &                               NRM(ifile,ng)%Rindex,              &
     &                               ioDesc,                            &
#  ifdef MASKING
     &                               vmask,                             &
#  endif
     &                               VnormV(:,:,:,ifile))
# endif
            END SELECT
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
!
!  3D norm at RHO-points.
!
          IF (Master) THEN
            Lsame=.FALSE.
            DO itrc=1,NT(ng)
              is=isTvar(itrc)
              IF (Cnorm(ifile,is)) Lsame=.TRUE.
            END DO
            IF (Lsame) THEN
              WRITE (stdout,20) TRIM(Text),                             &
     &                          '3D normalization factors at RHO-points'
              FLUSH (stdout)
            END IF
          END IF
!
!  Check if the decorrelation scales for all the tracers are different.
!  If not, just compute the normalization factors for the first tracer
!  and assign the same value to the rest.  Recall that this computation
!  is very expensive.
!
          Ldiffer=.FALSE.
          DO ns=1,Nscale(ng)
            DO itrc=2,NT(ng)
              IF ((HdecayX(ifile,isTvar(itrc  ),ns,ng).ne.              &
     &             HdecayX(ifile,isTvar(itrc-1),ns,ng)).or.             &
     &            (HdecayY(ifile,isTvar(itrc  ),ns,ng).ne.              &
     &             HdecayY(ifile,isTvar(itrc-1),ns,ng)).or.             &
     &            (Vdecay(ifile,isTvar(itrc  ),ng).ne.                  &
     &             Vdecay(ifile,isTvar(itrc-1),ng))) THEN
                Ldiffer=.TRUE.
              END IF
            END DO
            IF (.not.Ldiffer) THEN
              Lsame=.TRUE.
              UBt=1
            ELSE
              Lsame=.FALSE.
              UBt=NT(ng)
            END IF
          END DO
!
          DO j=JstrT,JendT
            DO i=IstrT,IendT
              val=om_r(i,j)*on_r(i,j)
              DO k=1,N(ng)
                Vscale(i,j,k)=1.0_r8/SQRT(val*Hz(i,j,k))
              END DO
            END DO
          END DO
          DO itrc=1,UBt
            is=isTvar(itrc)
            IF (Cnorm(ifile,is)) THEN
              DO ns=1,Nscale(ng)
                ms=ns
                iLap=Clap(is,ns,ng)
                DO k=1,N(ng)
                  DO j=JstrT,JendT
                    DO i=IstrT,IendT
                      A3davg(i,j,k)=0.0_r8
                      A3dsqr(i,j,k)=0.0_r8
                    END DO
                  END DO
                END DO
              END DO
              DO iter=1,Nrandom
                CALL white_noise3d (ng, iTLM, r3dvar, Rscheme(ng),      &
     &                              IstrR, IendR, JstrR, JendR,         &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              Amin, Amax, A3d)
                  DO k=1,N(ng)
                    DO j=JstrT,JendT
                      DO i=IstrT,IendT
                        A3d(i,j,k)=A3d(i,j,k)*Vscale(i,j,k)
                      END DO
                    END DO
                  END DO
                  DO k=1,N(ng)
                    klev=k
                    CALL dabc_r2d_tile (ng, tile,                       &
     &                                  LBi, UBi, LBj, UBj,             &
     &                                  A3d(:,:,klev))
                  END DO
!
                  CALL tl_CI_3d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj, 1, N(ng),     &
     &                                LBij, UBij,                       &
     &                                IminS, ImaxS, JminS, JmaxS,       &
     &                                Lweak, ifac, is,                  &
     &                                ms, MiterCI(ns,ng),iLap,          &
     &                                r3dvar, itrc,                     &
     &                                emaxr3d(:,:,itrc,ns),             &
     &                                eminr3d(:,:,itrc,ns),             &
     &                                ci_r3d, ci_p3d, ci_q3d, ci_x3d,   &
#  ifdef READ_SCALES
     &                                btsclx(:,:,itrc,ns),              &
     &                                btscly(:,:,itrc,ns),              &
#  endif
     &                                Hz, z_r, A3d)
!
                  CALL tl_conv_r3d_tile (ng, tile, iTLM,                &
     &                                   LBi, UBi, LBj, UBj, 1, N(ng),  &
     &                                   IminS, ImaxS, JminS, JmaxS,    &
     &                                   NghostPoints,                  &
     &                                   NHsteps(ifile,is)/ifac,        &
     &                                   NVsteps(ifile,is)/ifac,        &
     &                                   DTsizeH(ifile,is),             &
     &                                   DTsizeV(ifile,is),             &
     &                                   Kh, Kv,                        &
     &                                   pm, pn,                        &
#  ifdef GEOPOTENTIAL_HCONV
     &                                   on_u, om_v,                    &
#  else
     &                                   pmon_u, pnom_v,                &
#  endif
#  ifdef MASKING
     &                                   rmask, umask, vmask,           &
#  endif
     &                                   Hz, z_r,                       &
     &                                   A3d)
                  DO k=1,N(ng)
                    DO j=Jstr,Jend
                      DO i=Istr,Iend
                        A3davg(i,j,k)=A3davg(i,j,k)+A3d(i,j,k)
                        A3dsqr(i,j,k)=A3dsqr(i,j,k)+                    &
     &                                A3d(i,j,k)*A3d(i,j,k)
                      END DO
                    END DO
                  END DO
                END DO
                DO k=1,N(ng)
                  DO j=Jstr,Jend
                    DO i=Istr,Iend
                      Aavg=FacAvg*A3davg(i,j,k)
                      Asqr=FacAvg*A3dsqr(i,j,k)
#  ifdef MASKING
                      IF (rmask(i,j).gt.0.0_r8) THEN
                        VnormR(i,j,k,ifile,itrc)=                       &
     &                                    VnormR(i,j,k,ifile,itrc)+     &
     &                                    Bwgt(isTvar(itrc),ms,ng)/     &
     &                                    SQRT(Asqr)
                      ELSE
                        VnormR(i,j,k,ifile,itrc)=0.0_r8
                      END IF
#  else
                      VnormR(i,j,k,ifile,itrc)=                         &
     &                                  VnormR(i,j,k,ifile,itrc)+       &
     &                                  Bwgt(isTvar(itrc),ms,ng)/       &
     &                                  SQRT(Asqr)
#  endif
                    END DO
                  END DO
                END DO
              END DO
            END IF
          END DO
          IF (Lsame) THEN
            DO itrc=2,NT(ng)
              DO k=1,N(ng)
                DO j=Jstr,Jend
                  DO i=Istr,Iend
                    VnormR(i,j,k,ifile,itrc)=VnormR(i,j,k,ifile,1)
                  END DO
                END DO
              END DO
            END DO
          END IF
          DO itrc=1,NT(ng)
            is=isTvar(itrc)
            IF (Cnorm(ifile,is)) THEN
              CALL dabc_r3d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            VnormR(:,:,:,ifile,itrc))
# ifdef DISTRIBUTE
              CALL mp_exchange3d (ng, tile, iTLM, 1,                    &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            NghostPoints,                         &
     &                            EWperiodic(ng), NSperiodic(ng),       &
     &                            VnormR(:,:,:,ifile,itrc))
# endif
!
              SELECT CASE (NRM(ifile,ng)%IOtype)
                CASE (io_nf90)
                  CALL wrt_norm3d_nf90 (ng, tile, iTLM, ncname,         &
     &                                  LBi, UBi, LBj, UBj, 1, N(ng),   &
     &                                  idTvar(itrc),                   &
     &                              NRM(ifile,ng)%ncid,                 &
     &                              NRM(ifile,ng)%Vid(idTvar(itrc)),    &
     &                              NRM(ifile,ng)%Rindex,               &
# ifdef MASKING
     &                                  rmask,                          &
# endif
     &                                  VnormR(:,:,:,ifile,itrc))

# if defined PIO_LIB && defined DISTRIBUTE
                CASE (io_pio)
                  IF (NRM(ifile,ng)%pioTrc(itrc)%dkind.eq.              &
     &                PIO_double) THEN
                    ioDesc => ioDesc_dp_r3dvar(ng)
                  ELSE
                    ioDesc => ioDesc_sp_r3dvar(ng)
                  END IF
                  CALL wrt_norm3d_pio (ng, tile, iTLM, ncname,          &
     &                                 LBi, UBi, LBj, UBj, 1, N(ng),    &
     &                                 idTvar(itrc),                    &
     &                              NRM(ifile,ng)%pioFile,              &
     &                              NRM(ifile,ng)%pioTrc(itrc),         &
     &                              NRM(ifile,ng)%Rindex,               &
     &                                 ioDesc,                          &
#  ifdef MASKING
     &                                 rmask,                           &
#  endif
     &                                 VnormR(:,:,:,ifile,itrc))
# endif
              END SELECT
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END IF
          END DO
#endif
        END IF
      END DO FILE_LOOP

#ifdef ADJUST_BOUNDARY
!
!-----------------------------------------------------------------------
!  Compute open boundaries error covariance, B, normalization factors
!  using the randomization approach of Fisher and Courtier (1995).
!-----------------------------------------------------------------------
!
      ifile=3
      IF (LwrtNRM(ifile,ng)) THEN
        Text='boundary conditions'
        IJlen=UBij-LBij+1
# ifdef SOLVE3D
        IJKlen=IJlen*N(ng)
# endif
        Lconvolve(iwest )=DOMAIN(ng)%Western_Edge (tile)
        Lconvolve(ieast )=DOMAIN(ng)%Eastern_Edge (tile)
        Lconvolve(isouth)=DOMAIN(ng)%Southern_Edge(tile)
        Lconvolve(inorth)=DOMAIN(ng)%Northern_Edge(tile)
!
!  Set randomization summation factors.
!
        FacAvg=1.0_r8/REAL(Nrandom,r8)
        FacSqr=SQRT(REAL(Nrandom,r8))
!
!  Set time record index to write in normalization NetCDF file.
!
        ncname=NRM(ifile,ng)%name
        NRM(ifile,ng)%Rindex=NRM(ifile,ng)%Rindex+1
        NRM(ifile,ng)%Nrec=NRM(ifile,ng)%Nrec+1
!
!  Write out model time (s).
!
        SELECT CASE (NRM(ifile,ng)%IOtype)
          CASE (io_nf90)
            CALL netcdf_put_fvar (ng, iTLM, ncname,                     &
     &                            Vname(1,idtime), my_time,             &
     &                         start = (/NRM(ifile,ng)%Rindex/),        &
     &                         total = (/1/),                           &
     &                         ncid = NRM(ifile,ng)%ncid,               &
     &                         varid = NRM(ifile,ng)%Vid(idtime))

# if defined PIO_LIB && defined DISTRIBUTE
          CASE (io_pio)
            CALL pio_netcdf_put_fvar (ng, iTLM, ncname,                 &
     &                                Vname(1,idtime), my_time,         &
     &                         start = (/NRM(ifile,ng)%Rindex/),        &
     &                         total = (/1/),                           &
     &                         pioFile = NRM(ifile,ng)%pioFile,         &
     &                         pioVar = NRM(ifile,ng)%pioVar(idtime)%vd)
# endif
        END SELECT
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  2D boundary norm at RHO-points.
!
        HnormRobc=Aspv

        IF (Master.and.ANY(CnormB(isFsur,:))) THEN
          WRITE (stdout,20) TRIM(Text),                                 &
     &                      '2D normalization factors at RHO-points'
          FLUSH (stdout)
        END IF

        DO ibry=1,4
          IF (CnormB(isFsur,ibry)) THEN
            DO ns=1,Nscale(ng)
              ms=ns
              iLap=Clap(isFsur,ns,ng)
              HscaleB=0.0_r8
              B2davg=0.0_r8
              B2dsqr=0.0_r8
              B2d=0.0_r8
              IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
                i=BOUNDS(ng)%edge(ibry,r2dvar)
                IF (Lconvolve(ibry)) THEN
                  DO j=JstrT,JendT
                    HscaleB(j)=1.0_r8/SQRT(on_r(i,j))
                  END DO
                END IF
              ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
                j=BOUNDS(ng)%edge(ibry,r2dvar)
                IF (Lconvolve(ibry)) THEN
                  DO i=IstrT,IendT
                    HscaleB(i)=1.0_r8/SQRT(om_r(i,j))
                  END DO
                END IF
              END IF
              DO iter=1,Nrandom
                IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
                  CALL white_noise2d_bry (ng, tile, iTLM, ibry,         &
     &                                  Rscheme(ng),                    &
     &                                  JstrR, JendR,                   &
     &                                  LBij, UBij,                     &
     &                                  Bmin, Bmax, B2d)
                  DO j=JstrT,JendT
                    B2d(j)=B2d(j)*HscaleB(j)
                  END DO
                ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
                  CALL white_noise2d_bry (ng, tile, iTLM, ibry,         &
     &                                    Rscheme(ng),                  &
     &                                    IstrR, IendR,                 &
     &                                    LBij, UBij,                   &
     &                                    Bmin, Bmax, B2d)
                  DO i=IstrT,IendT
                    B2d(i)=B2d(i)*HscaleB(i)
                  END DO
                END IF
                CALL tl_CI_bry1d_tile (ng, tile,                        &
     &                                 LBi, UBi, LBj, UBj, LBij, UBij,  &
     &                                 IminS, ImaxS, JminS, JmaxS,      &
     &                                 nstp, nnew, ifac, isFsur, ibry,  &
     &                                 ms, MiterCI(ns,ng), iLap,        &
     &                                 r2dvar, 1, 1,                    &
     &                                 emaxr1d(:,ibry,ns),              &
     &                                 eminr1d(:,ibry,ns),              &
     &                                 ci_rb, ci_pb, ci_qb, ci_xb,      &
#  ifdef SOLVE3D
     &                                 Hz,                              &
#  endif
     &                                 B2d)
                IF (Lconvolve(ibry)) THEN
                  IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
                    DO j=Jstr,Jend
                      B2davg(j)=B2davg(j)+B2d(j)
                      B2dsqr(j)=B2dsqr(j)+B2d(j)*B2d(j)
                    END DO
                  ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
                    DO i=Istr,Iend
                      B2davg(i)=B2davg(i)+B2d(i)
                      B2dsqr(i)=B2dsqr(i)+B2d(i)*B2d(i)
                    END DO
                  END IF
                END IF
              END DO
              IF (Lconvolve(ibry)) THEN
                IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
                  DO j=Jstr,Jend
                    Bavg=FacAvg*B2davg(j)
                    Bsqr=FacAvg*B2dsqr(j)
#  ifdef MASKING
                    IF (rmask(i,j).gt.0.0_r8) THEN
                    HnormRobc(j,ibry)=HnormRobc(j,ibry)+                &
     &                                Bwgt(isFsur,ms,ng)/SQRT(Bsqr)
                    ELSE
                      HnormRobc(j,ibry)=0.0_r8
                    END IF
#  else
                    HnormRobc(j,ibry)=HnormRobc(j,ibry)+                &
     &                                Bwgt(isFsur,ms,ng)/SQRT(Bsqr)
#  endif
                  END DO
                ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
                  DO i=Istr,Iend
                    Bavg=FacAvg*B2davg(i)
                    Bsqr=FacAvg*B2dsqr(i)
#  ifdef MASKING
                    IF (rmask(i,j).gt.0.0_r8) THEN
                      HnormRobc(i,ibry)=HnormRobc(i,ibry)+              &
     &                                  Bwgt(isFsur,ms,ng)/SQRT(Bsqr)
                    ELSE
                      HnormRobc(i,ibry)=0.0_r8
                    END IF
#  else
                    HnormRobc(i,ibry)=HnormRobc(i,ibry)+                &
     &                                Bwgt(isFsur,ms,ng)/SQRT(Bsqr)
#  endif
                  END DO
                END IF
              END IF
            END DO   ! end of ns-loop
            CALL bc_r2d_bry_tile (ng, tile, ibry,                       &
     &                            LBij, UBij,                           &
     &                            HnormRobc(:,ibry))
# ifdef DISTRIBUTE
            CALL mp_collect (ng, iTLM, IJlen, Aspv,                     &
     &                       HnormRobc(LBij:,ibry))
# endif
          END IF
        END DO
        IF (ANY(CnormB(isFsur,:))) THEN
          ifield=idSbry(isFsur)

          SELECT CASE (NRM(ifile,ng)%IOtype)
            CASE (io_nf90)
              CALL netcdf_put_fvar (ng, iTLM, ncname,                   &
     &                              Vname(1,ifield),                    &
     &                              HnormRobc(LBij:,:),                 &
     &                         start = (/1,1,NRM(ifile,ng)%Rindex/),    &
     &                         total = (/IJlen,4,1/),                   &
     &                         ncid = NRM(ifile,ng)%ncid,               &
     &                         varid = NRM(ifile,ng)%Vid(ifield))

# if defined PIO_LIB && defined DISTRIBUTE
            CASE (io_pio)
              CALL pio_netcdf_put_fvar (ng, iTLM, ncname,               &
     &                                  Vname(1,ifield),                &
     &                                  HnormRobc(LBij:,:),             &
     &                         start = (/1,1,NRM(ifile,ng)%Rindex/),    &
     &                         total = (/IJlen,4,1/),                   &
     &                         pioFile = NRM(ifile,ng)%pioFile,         &
     &                         pioVar = NRM(ifile,ng)%pioVar(ifield)%vd)

# endif
          END SELECT
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  2D boundary norm at U-points.
!
        HnormUobc=Aspv

        IF (Master.and.ANY(CnormB(isUbar,:))) THEN
          WRITE (stdout,20) TRIM(Text),                                 &
     &                      '2D normalization factors at   U-points'
          FLUSH (stdout)
        END IF

        DO ibry=1,4
          IF (CnormB(isUbar,ibry)) THEN
            DO ns=1,Nscale(ng)
              ms=ns
              iLap=Clap(isUbar,ns,ng)
              HscaleB=0.0_r8
              B2davg=0.0_r8
              B2dsqr=0.0_r8
              B2d=0.0_r8
              IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
                i=BOUNDS(ng)%edge(ibry,u2dvar)
                IF (Lconvolve(ibry)) THEN
                  DO j=JstrT,JendT
                    HscaleB(j)=1.0_r8/SQRT(on_u(i,j))
                  END DO
                END IF
              ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
                j=BOUNDS(ng)%edge(ibry,u2dvar)
                IF (Lconvolve(ibry)) THEN
                  DO i=IstrP,IendT
                    HscaleB(i)=1.0_r8/SQRT(om_u(i,j))
                  END DO
                END IF
              END IF
              DO iter=1,Nrandom
                IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
                  CALL white_noise2d_bry (ng, tile, iTLM, ibry,         &
     &                                    Rscheme(ng),                  &
     &                                    JstrR, JendR,                 &
     &                                    LBij, UBij,                   &
     &                                    Bmin, Bmax, B2d)
                  DO j=JstrT,JendT
                    B2d(j)=B2d(j)*HscaleB(j)
                  END DO
                ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
                  CALL white_noise2d_bry (ng, tile, iTLM, ibry,         &
     &                                    Rscheme(ng),                  &
     &                                    Istr, IendR,                  &
     &                                    LBij, UBij,                   &
     &                                    Bmin, Bmax, B2d)
                  DO i=IstrP,IendT
                    B2d(i)=B2d(i)*HscaleB(i)
                  END DO
                END IF
                CALL tl_CI_bry1d_tile (ng, tile,                        &
     &                                 LBi, UBi, LBj, UBj, LBij, UBij,  &
     &                                 IminS, ImaxS, JminS, JmaxS,      &
     &                                 nstp, nnew, ifac, isUbar, ibry,  &
     &                                 ms, MiterCI(ns,ng), iLap,        &
     &                                 u2dvar, 1, 1,                    &
     &                                 emaxu1d(:,ibry,ns),              &
     &                                 eminu1d(:,ibry,ns),              &
     &                                 ci_rb, ci_pb, ci_qb, ci_xb,      &
#   ifdef SOLVE3D
     &                                 Hz,                              &
#   endif
     &                                 B2d)
                IF (Lconvolve(ibry)) THEN
                  IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
                    DO j=Jstr,Jend
                      B2davg(j)=B2davg(j)+B2d(j)
                      B2dsqr(j)=B2dsqr(j)+B2d(j)*B2d(j)
                    END DO
                  ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
                    DO i=IstrU,Iend
                      B2davg(i)=B2davg(i)+B2d(i)
                      B2dsqr(i)=B2dsqr(i)+B2d(i)*B2d(i)
                    END DO
                  END IF
                END IF
              END DO
              IF (Lconvolve(ibry)) THEN
                IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
                  DO j=Jstr,Jend
                    Bavg=FacAvg*B2davg(j)
                    Bsqr=FacAvg*B2dsqr(j)
#  ifdef MASKING
                    IF (umask(i,j).gt.0.0_r8) THEN
                      HnormUobc(j,ibry)=HnormUobc(j,ibry)+              &
     &                                  Bwgt(isUbar,ms,ng)/SQRT(Bsqr)
                    ELSE
                      HnormUobc(j,ibry)=0.0_r8
                    END IF
#  else
                    HnormUobc(j,ibry)=HnormUobc(j,ibry)+                &
     &                                Bwgt(isUbar,ms,ng)/SQRT(Bsqr)
# endif
                  END DO
                ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
                  DO i=IstrU,Iend
                    Bavg=FacAvg*B2davg(i)
                    Bsqr=FacAvg*B2dsqr(i)
#  ifdef MASKING
                    IF (umask(i,j).gt.0.0_r8) THEN
                      HnormUobc(i,ibry)=HnormUobc(i,ibry)+              &
     &                                  Bwgt(isUbar,ms,ng)/SQRT(Bsqr)
                    ELSE
                      HnormUobc(i,ibry)=0.0_r8
                    END IF
#  else
                    HnormUobc(i,ibry)=HnormUobc(i,ibry)+                &
     &                                Bwgt(isUbar,ms,ng)/SQRT(Bsqr)
# endif
                  END DO
                END IF
              END IF
            END DO  ! end of ns-loop.
            CALL bc_u2d_bry_tile (ng, tile, ibry,                       &
     &                            LBij, UBij,                           &
     &                            HnormUobc(:,ibry))
# ifdef DISTRIBUTE
            CALL mp_collect (ng, iTLM, IJlen, Aspv,                     &
     &                       HnormUobc(LBij:,ibry))
# endif
          END IF
        END DO
        IF (ANY(CnormB(isUbar,:))) THEN
          ifield=idSbry(isUbar)

          SELECT CASE (NRM(ifile,ng)%IOtype)
            CASE (io_nf90)
              CALL netcdf_put_fvar (ng, iTLM, ncname,                   &
     &                              Vname(1,ifield),                    &
     &                              HnormUobc(LBij:,:),                 &
     &                         start = (/1,1,NRM(ifile,ng)%Rindex/),    &
     &                         total = (/IJlen,4,1/),                   &
     &                         ncid = NRM(ifile,ng)%ncid,               &
     &                         varid = NRM(ifile,ng)%Vid(ifield))

# if defined PIO_LIB && defined DISTRIBUTE
            CASE (io_pio)
              CALL pio_netcdf_put_fvar (ng, iTLM, ncname,               &
     &                                  Vname(1,ifield),                &
     &                                  HnormUobc(LBij:,:),             &
     &                         start = (/1,1,NRM(ifile,ng)%Rindex/),    &
     &                         total = (/IJlen,4,1/),                   &
     &                         pioFile = NRM(ifile,ng)%pioFile,         &
     &                         pioVar = NRM(ifile,ng)%pioVar(ifield)%vd)
# endif
          END SELECT
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  2D boundary norm at V-points.
!
        HnormVobc=Aspv

        IF (Master.and.ANY(CnormB(isVbar,:))) THEN
          WRITE (stdout,20) TRIM(Text),                                 &
     &                      '2D normalization factors at   V-points'
          FLUSH (stdout)
        END IF

        DO ibry=1,4
          IF (CnormB(isVbar,ibry)) THEN
            DO ns=1,Nscale(ng)
              ms=ns
              iLap=Clap(isVbar,ns,ng)
              HscaleB=0.0_r8
              B2davg=0.0_r8
              B2dsqr=0.0_r8
              B2d=0.0_r8
              IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
                i=BOUNDS(ng)%edge(ibry,v2dvar)
                IF (Lconvolve(ibry)) THEN
                  DO j=JstrP,JendT
                    HscaleB(j)=1.0_r8/SQRT(on_v(i,j))
                  END DO
                END IF
              ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
                j=BOUNDS(ng)%edge(ibry,v2dvar)
                IF (Lconvolve(ibry)) THEN
                  DO i=IstrT,IendT
                    HscaleB(i)=1.0_r8/SQRT(om_v(i,j))
                  END DO
                END IF
              END IF
              DO iter=1,Nrandom
                IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
                  CALL white_noise2d_bry (ng, tile, iTLM, ibry,         &
     &                                    Rscheme(ng),                  &
     &                                    Jstr, JendR,                  &
     &                                    LBij, UBij,                   &
     &                                    Bmin, Bmax, B2d)
                  DO j=JstrP,JendT
                    B2d(j)=B2d(j)*HscaleB(j)
                  END DO
                ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
                  CALL white_noise2d_bry (ng, tile, iTLM, ibry,         &
     &                                    Rscheme(ng),                  &
     &                                    IstrR, IendR,                 &
     &                                    LBij, UBij,                   &
     &                                    Bmin, Bmax, B2d)
                  DO i=IstrT,IendT
                    B2d(i)=B2d(i)*HscaleB(i)
                    END DO
                END IF
                CALL tl_CI_bry1d_tile (ng, tile,                        &
     &                                 LBi, UBi, LBj, UBj, LBij, UBij,  &
     &                                 IminS, ImaxS, JminS, JmaxS,      &
     &                                 nstp, nnew, ifac, isVbar, ibry,  &
     &                                 ms, MiterCI(ns,ng), iLap,        &
     &                                 v2dvar, 1, 1,                    &
     &                                 emaxv1d(:,ibry,ns),              &
     &                                 eminv1d(:,ibry,ns),              &
     &                                 ci_rb, ci_pb, ci_qb, ci_xb,      &
#  ifdef SOLVE3D
     &                                 Hz,                              &
#  endif
     &                                 B2d)

                IF (Lconvolve(ibry)) THEN
                  IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
                    DO j=JstrV,Jend
                      B2davg(j)=B2davg(j)+B2d(j)
                      B2dsqr(j)=B2dsqr(j)+B2d(j)*B2d(j)
                    END DO
                  ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
                    DO i=Istr,Iend
                      B2davg(i)=B2davg(i)+B2d(i)
                      B2dsqr(i)=B2dsqr(i)+B2d(i)*B2d(i)
                    END DO
                  END IF
                END IF
              END DO
              IF (Lconvolve(ibry)) THEN
                IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
                  DO j=JstrV,Jend
                    Bavg=FacAvg*B2davg(j)
                    Bsqr=FacAvg*B2dsqr(j)
#  ifdef MASKING
                    IF (vmask(i,j).gt.0.0_r8) THEN
                      HnormVobc(j,ibry)=HnormVobc(j,ibry)+              &
     &                                  Bwgt(isVbar,ms,ng)/SQRT(Bsqr)
                    ELSE
                      HnormVobc(j,ibry)=0.0_r8
                    END IF
#  else
                    HnormVobc(j,ibry)=HnormVobc(j,ibry)+                &
     &                                Bwgt(isVbar,ms,ng)/SQRT(Bsqr)
#  endif
                  END DO
                ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
                  DO i=Istr,Iend
                    Bavg=FacAvg*B2davg(i)
                    Bsqr=FacAvg*B2dsqr(i)
#  ifdef MASKING
                    IF (vmask(i,j).gt.0.0_r8) THEN
                      HnormVobc(i,ibry)=HnormVobc(i,ibry)+              &
     &                                  Bwgt(isVbar,ms,ng)/SQRT(Bsqr)
                    ELSE
                      HnormVobc(i,ibry)=0.0_r8
                    END IF
#  else
                    HnormVobc(i,ibry)=HnormVobc(i,ibry)+                &
     &                                Bwgt(isVbar,ms,ng)/SQRT(Bsqr)
#  endif
                  END DO
                END IF
              END IF
            END DO  ! end ns-loop.
            CALL bc_v2d_bry_tile (ng, tile, ibry,                       &
     &                            LBij, UBij,                           &
     &                            HnormVobc(:,ibry))
# ifdef DISTRIBUTE
            CALL mp_collect (ng, iTLM, IJlen, Aspv,                     &
     &                       HnormVobc(LBij:,ibry))
# endif
          END IF
        END DO
        IF (ANY(CnormB(isVbar,:))) THEN
          ifield=idSbry(isVbar)

          SELECT CASE (NRM(ifile,ng)%IOtype)
            CASE (io_nf90)
              CALL netcdf_put_fvar (ng, iTLM, ncname,                   &
     &                              Vname(1,ifield),                    &
     &                              HnormVobc(LBij:,:),                 &
     &                         start = (/1,1,NRM(ifile,ng)%Rindex/),    &
     &                         total = (/IJlen,4,1/),                   &
     &                         ncid = NRM(ifile,ng)%ncid,               &
     &                         varid = NRM(ifile,ng)%Vid(ifield))

# if defined PIO_LIB && defined DISTRIBUTE
            CASE (io_pio)
              CALL pio_netcdf_put_fvar (ng, iTLM, ncname,               &
     &                                  Vname(1,ifield),                &
     &                                  HnormVobc(LBij:,:),             &
     &                         start = (/1,1,NRM(ifile,ng)%Rindex/),    &
     &                         total = (/IJlen,4,1/),                   &
     &                         pioFile = NRM(ifile,ng)%pioFile,         &
     &                         pioVar = NRM(ifile,ng)%pioVar(ifield)%vd)
# endif
          END SELECT
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

# ifdef SOLVE3D
!
!  3D boundary norm at U-points.
!
        VnormUobc=Aspv

        IF (Master.and.ANY(CnormB(isUvel,:))) THEN
          WRITE (stdout,20) TRIM(Text),                                 &
     &                      '3D normalization factors at   U-points'
          FLUSH (stdout)
        END IF

        DO ibry=1,4
          IF (CnormB(isUvel,ibry)) THEN
            DO ns=1,Nscale(ng)
              ms=ns
              iLap=Clap(isUvel,ns,ng)
              VscaleB=0.0_r8
              B3davg=0.0_r8
              B3dsqr=0.0_r8
              B3d=0.0_r8
              IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
                i=BOUNDS(ng)%edge(ibry,u2dvar)
                IF (Lconvolve(ibry)) THEN
                  DO j=JstrT,JendT
                    val=on_u(i,j)*0.5_r8
                    DO k=1,N(ng)
                      VscaleB(j,k)=1.0_r8/                              &
     &                             SQRT(val*(Hz(i-1,j,k)+Hz(i,j,k)))
                    END DO
                  END DO
                END IF
              ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
                j=BOUNDS(ng)%edge(ibry,u2dvar)
                IF (Lconvolve(ibry)) THEN
                  DO i=IstrP,IendT
                    val=om_u(i,j)*0.5_r8
                    DO k=1,N(ng)
                      VscaleB(i,k)=1.0_r8/                              &
     &                             SQRT(val*(Hz(i-1,j,k)+Hz(i,j,k)))
                    END DO
                  END DO
                END IF
              END IF
              DO iter=1,Nrandom
                IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
                  CALL white_noise3d_bry (ng, tile, iTLM, ibry,         &
     &                                    Rscheme(ng),                  &
     &                                    JstrR, JendR,                 &
     &                                    LBij, UBij, 1, N(ng),         &
     &                                    Bmin, Bmax, B3d)
                  DO k=1,N(ng)
                    DO j=JstrT,JendT
                      B3d(j,k)=B3d(j,k)*VscaleB(j,k)
                    END DO
                  END DO
                ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
                  CALL white_noise3d_bry (ng, tile, iTLM, ibry,         &
     &                                    Rscheme(ng),                  &
     &                                    Istr, IendR,                  &
     &                                    LBij, UBij, 1, N(ng),         &
     &                                    Bmin, Bmax, B3d)
                  DO k=1,N(ng)
                    DO i=IstrP,IendT
                      B3d(i,k)=B3d(i,k)*VscaleB(i,k)
                    END DO
                  END DO
                END IF
!
!  First the horizontal convolution.
!
                DO k=1,N(ng)
                  klev=k
                  CALL tl_CI_bry1d_tile (ng, tile,                      &
     &                                   LBi, UBi, LBj, UBj,            &
     &                                   LBij, UBij,                    &
     &                                   IminS, ImaxS, JminS, JmaxS,    &
     &                                   nstp, nnew, ifac,              &
     &                                   isUvel, ibry,                  &
     &                                   ms, MiterCI(ns,ng), iLap,      &
     &                                   u3dvar, klev, 1,               &
     &                                   emaxu1dz(:,klev,ibry,ns),      &
     &                                   eminu1dz(:,klev,ibry,ns),      &
     &                                   ci_rb, ci_pb, ci_qb, ci_xb,    &
#  ifdef SOLVE3D
     &                                   Hz,                            &
#  endif
     &                                   B3d(:,klev))
                END DO
!
!  Now the vertical convolution.
!
                CALL tl_conv_u3d_bry_tile (ng, tile, iTLM, ibry,        &
     &                                     BOUNDS(ng)%edge(:,u2dvar),   &
     &                                     LBij, UBij,                  &
     &                                     LBi, UBi, LBj, UBj, 1, N(ng),&
     &                                     IminS, ImaxS, JminS, JmaxS,  &
     &                                     NghostPoints,                &
     &                                     NHstepsB(ibry,isUvel)/ifac,  &
     &                                     NVstepsB(ibry,isUvel)/ifac,  &
     &                                     DTsizeHB(ibry,isUvel),       &
     &                                     DTsizeVB(ibry,isUvel),       &
     &                                     Kh, Kv,                      &
     &                                     pm, pn,                      &
     &                                     pmon_r, pnom_p,              &
#   ifdef MASKING
     &                                     umask, pmask,                &
#   endif
     &                                     Hz, z_r,                     &
     &                                     B3d)
                IF (Lconvolve(ibry)) THEN
                  IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
                    DO k=1,N(ng)
                      DO j=Jstr,Jend
                        B3davg(j,k)=B3davg(j,k)+B3d(j,k)
                        B3dsqr(j,k)=B3dsqr(j,k)+B3d(j,k)*B3d(j,k)
                      END DO
                    END DO
                  ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
                    DO k=1,N(ng)
                      DO i=IstrU,Iend
                        B3davg(i,k)=B3davg(i,k)+B3d(i,k)
                        B3dsqr(i,k)=B3dsqr(i,k)+B3d(i,k)*B3d(i,k)
                      END DO
                    END DO
                  END IF
                END IF
              END DO
              IF (Lconvolve(ibry)) THEN
                IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
                  DO k=1,N(ng)
                    DO j=Jstr,Jend
                      Bavg=FacAvg*B3davg(j,k)
                      Bsqr=FacAvg*B3dsqr(j,k)
#   ifdef MASKING
                      IF (umask(i,j).gt.0.0_r8) THEN
                        VnormUobc(j,k,ibry)=VnormUobc(j,k,ibry)+        &
     &                                      Bwgt(isUvel,ms,ng)/         &
     &                                      SQRT(Bsqr)
                      ELSE
                        VnormUobc(j,k,ibry)=0.0_r8
                      END IF
#   else
                      VnormUobc(j,k,ibry)=VnormUobc(j,k,ibry)+          &
     &                                    Bwgt(isUvel,ms,ng)/           &
     &                                    SQRT(Bsqr)
#   endif
                    END DO
                  END DO
                ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
                  DO k=1,N(ng)
                    DO i=IstrU,Iend
                      Bavg=FacAvg*B3davg(i,k)
                      Bsqr=FacAvg*B3dsqr(i,k)
#   ifdef MASKING
                      IF (umask(i,j).gt.0.0_r8) THEN
                        VnormUobc(i,k,ibry)=VnormUobc(i,k,ibry)+        &
     &                                      Bwgt(isUvel,ms,ng)/         &
     &                                      SQRT(Bsqr)
                      ELSE
                        VnormUobc(i,k,ibry)=0.0_r8
                      END IF
#   else
                      VnormUobc(i,k,ibry)=VnormUobc(i,k,ibry)+          &
     &                                    Bwgt(isUvel,ms,ng)/           &
     &                                    SQRT(Bsqr)
#   endif
                    END DO
                  END DO
                END IF
              END IF
            END DO    !  end of ns-loop
            CALL bc_u3d_bry_tile (ng, tile, ibry,                       &
     &                            LBij, UBij, 1, N(ng),                 &
     &                            VnormUobc(:,:,ibry))
#  ifdef DISTRIBUTE
            Bwrk=RESHAPE(VnormUobc(:,:,ibry), (/IJKlen/))
            CALL mp_collect (ng, iTLM, IJKlen, Aspv, Bwrk)
            ic=0
            DO k=1,N(ng)
              DO ib=LBij,UBij
                ic=ic+1
                VnormUobc(ib,k,ibry)=Bwrk(ic)
              END DO
            END DO
#  endif
          END IF
        END DO
        IF (ANY(CnormB(isUvel,:))) THEN
          ifield=idSbry(isUvel)

          SELECT CASE (NRM(ifile,ng)%IOtype)
            CASE (io_nf90)
              CALL netcdf_put_fvar (ng, iTLM, ncname,                   &
     &                              Vname(1,ifield),                    &
     &                              VnormUobc(LBij:,:,:),               &
     &                         start = (/1,1,1,NRM(ifile,ng)%Rindex/),  &
     &                         total = (/IJlen,N(ng),4,1/),             &
     &                         ncid = NRM(ifile,ng)%ncid,               &
     &                         varid = NRM(ifile,ng)%Vid(ifield))

#  if defined PIO_LIB && defined DISTRIBUTE
            CASE (io_pio)
              CALL pio_netcdf_put_fvar (ng, iTLM, ncname,               &
     &                                  Vname(1,ifield),                &
     &                                  VnormUobc(LBij:,:,:),           &
     &                         start = (/1,1,1,NRM(ifile,ng)%Rindex/),  &
     &                         total = (/IJlen,N(ng),4,1/),             &
     &                         pioFile = NRM(ifile,ng)%pioFile,         &
     &                         pioVar = NRM(ifile,ng)%pioVar(ifield)%vd)
#  endif
          END SELECT
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  3D boundary norm at V-points.
!
        VnormVobc=Aspv

        IF (Master.and.ANY(CnormB(isVvel,:))) THEN
          WRITE (stdout,20) TRIM(Text),                                 &
     &                      '3D normalization factors at   V-points'
          FLUSH (stdout)
        END IF

        DO ibry=1,4
          IF (CnormB(isVvel,ibry)) THEN
            DO ns=1,Nscale(ng)
              ms=ns
              iLap=Clap(isVvel,ns,ng)
              VscaleB=0.0_r8
              B3davg=0.0_r8
              B3dsqr=0.0_r8
              B3d=0.0_r8
              IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
                i=BOUNDS(ng)%edge(ibry,v2dvar)
                IF (Lconvolve(ibry)) THEN
                  DO j=JstrP,JendT
                    val=on_v(i,j)*0.5_r8
                    DO k=1,N(ng)
                      VscaleB(j,k)=1.0_r8/                              &
     &                             SQRT(val*(Hz(i,j-1,k)+Hz(i,j,k)))
                    END DO
                  END DO
                END IF
              ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
                j=BOUNDS(ng)%edge(ibry,v2dvar)
                IF (Lconvolve(ibry)) THEN
                  DO i=IstrT,IendT
                    val=om_v(i,j)*0.5_r8
                    DO k=1,N(ng)
                      VscaleB(i,k)=1.0_r8/                              &
     &                             SQRT(val*(Hz(i,j-1,k)+Hz(i,j,k)))
                    END DO
                  END DO
                END IF
              END IF
              DO iter=1,Nrandom
                IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
                  CALL white_noise3d_bry (ng, tile, iTLM, ibry,         &
     &                                    Rscheme(ng),                  &
     &                                    Jstr, JendR,                  &
     &                                    LBij, UBij, 1, N(ng),         &
     &                                    Bmin, Bmax, B3d)
                  DO k=1,N(ng)
                    DO j=JstrP,JendT
                      B3d(j,k)=B3d(j,k)*VscaleB(j,k)
                    END DO
                  END DO
                ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
                  CALL white_noise3d_bry (ng, tile, iTLM, ibry,         &
     &                                    Rscheme(ng),                  &
     &                                    IstrR, IendR,                 &
     &                                    LBij, UBij, 1, N(ng),         &
     &                                    Bmin, Bmax, B3d)
                  DO k=1,N(ng)
                    DO i=IstrT,IendT
                      B3d(i,k)=B3d(i,k)*VscaleB(i,k)
                    END DO
                  END DO
                END IF
!
!  First the horizontal convolution.
!
                DO k=1,N(ng)
                  klev=k
                  CALL tl_CI_bry1d_tile (ng, tile,                      &
     &                                   LBi, UBi, LBj, UBj,            &
     &                                   LBij, UBij,                    &
     &                                   IminS, ImaxS, JminS, JmaxS,    &
     &                                   nstp, nnew, ifac,              &
     &                                   isVvel, ibry,                  &
     &                                   ms, MiterCI(ns,ng), iLap,      &
     &                                   v3dvar, klev, 1,               &
     &                                   emaxv1dz(:,klev,ibry,ns),      &
     &                                   eminv1dz(:,klev,ibry,ns),      &
     &                                   ci_rb, ci_pb, ci_qb, ci_xb,    &
#  ifdef SOLVE3D
     &                                   Hz,                            &
#  endif
     &                                   B3d(:,klev))
                END DO
!
!  Now the vertical convolution.
!
                CALL tl_conv_v3d_bry_tile (ng, tile, iTLM, ibry,        &
     &                                     BOUNDS(ng)%edge(:,v2dvar),   &
     &                                     LBij, UBij,                  &
     &                                     LBi, UBi, LBj, UBj, 1, N(ng),&
     &                                     IminS, ImaxS, JminS, JmaxS,  &
     &                                     NghostPoints,                &
     &                                     NHstepsB(ibry,isVvel)/ifac,  &
     &                                     NVstepsB(ibry,isVvel)/ifac,  &
     &                                     DTsizeHB(ibry,isVvel),       &
     &                                     DTsizeVB(ibry,isVvel),       &
     &                                     Kh, Kv,                      &
     &                                     pm, pn,                      &
     &                                     pmon_p, pnom_r,              &
#   ifdef MASKING
     &                                     vmask, pmask,                &
#   endif
     &                                     Hz, z_r,                     &
     &                                     B3d)
                IF (Lconvolve(ibry)) THEN
                  IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
                    DO k=1,N(ng)
                      DO j=JstrV,Jend
                        B3davg(j,k)=B3davg(j,k)+B3d(j,k)
                        B3dsqr(j,k)=B3dsqr(j,k)+B3d(j,k)*B3d(j,k)
                      END DO
                    END DO
                  ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
                    DO k=1,N(ng)
                      DO i=Istr,Iend
                        B3davg(i,k)=B3davg(i,k)+B3d(i,k)
                        B3dsqr(i,k)=B3dsqr(i,k)+B3d(i,k)*B3d(i,k)
                      END DO
                    END DO
                  END IF
                END IF
              END DO
              IF (Lconvolve(ibry)) THEN
                IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
                  DO k=1,N(ng)
                    DO j=JstrV,Jend
                      Bavg=FacAvg*B3davg(j,k)
                      Bsqr=FacAvg*B3dsqr(j,k)
#   ifdef MASKING
                      IF (vmask(i,j).gt.0.0_r8) THEN
                        VnormVobc(j,k,ibry)=VnormVobc(j,k,ibry)+        &
     &                                      Bwgt(isVvel,ms,ng)/         &
     &                                      SQRT(Bsqr)
                      ELSE
                        VnormVobc(j,k,ibry)=0.0_r8
                      END IF
#   else
                      VnormVobc(j,k,ibry)=VnormVobc(j,k,ibry)+          &
     &                                    Bwgt(isVvel,ms,ng)/           &
     &                                    SQRT(Bsqr)
#   endif
                    END DO
                  END DO
                ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
                  DO k=1,N(ng)
                    DO i=Istr,Iend
                      Bavg=FacAvg*B3davg(i,k)
                      Bsqr=FacAvg*B3dsqr(i,k)
#   ifdef MASKING
                      IF (vmask(i,j).gt.0.0_r8) THEN
                        VnormVobc(i,k,ibry)=VnormVobc(i,k,ibry)+        &
     &                                      Bwgt(isVvel,ms,ng)/         &
     &                                      SQRT(Bsqr)
                      ELSE
                        VnormVobc(i,k,ibry)=0.0_r8
                      END IF
#   else
                      VnormVobc(i,k,ibry)=VnormVobc(i,k,ibry)+          &
     &                                    Bwgt(isVvel,ms,ng)/           &
     &                                    SQRT(Bsqr)
#   endif
                    END DO
                  END DO
                END IF
              END IF
            END DO    ! end ns-loop
            CALL bc_v3d_bry_tile (ng, tile, ibry,                       &
     &                            LBij, UBij, 1, N(ng),                 &
     &                            VnormVobc(:,:,ibry))
#  ifdef DISTRIBUTE
            Bwrk=RESHAPE(VnormVobc(:,:,ibry), (/IJKlen/))
            CALL mp_collect (ng, iTLM, IJKlen, Aspv, Bwrk)
            ic=0
            DO k=1,N(ng)
              DO ib=LBij,UBij
                ic=ic+1
                VnormVobc(ib,k,ibry)=Bwrk(ic)
              END DO
            END DO
#  endif
          END IF
        END DO
        IF (ANY(CnormB(isVvel,:))) THEN
          ifield=idSbry(isVvel)

          SELECT CASE (NRM(ifile,ng)%IOtype)
            CASE (io_nf90)
              CALL netcdf_put_fvar (ng, iTLM, ncname,                   &
     &                              Vname(1,ifield),                    &
     &                              VnormVobc(LBij:,:,:),               &
     &                         start = (/1,1,1,NRM(ifile,ng)%Rindex/),  &
     &                         total = (/IJlen,N(ng),4,1/),             &
     &                         ncid = NRM(ifile,ng)%ncid,               &
     &                         varid = NRM(ifile,ng)%Vid(ifield))

#  if defined PIO_LIB && defined DISTRIBUTE
            CASE (io_pio)
              CALL pio_netcdf_put_fvar (ng, iTLM, ncname,               &
     &                                  Vname(1,ifield),                &
     &                                  VnormVobc(LBij:,:,:),           &
     &                         start = (/1,1,1,NRM(ifile,ng)%Rindex/),  &
     &                         total = (/IJlen,N(ng),4,1/),             &
     &                         pioFile = NRM(ifile,ng)%pioFile,         &
     &                         pioVar = NRM(ifile,ng)%pioVar(ifield)%vd)
#  endif
          END SELECT
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  3D boundary norm at RHO-points.
!
        IF (Master) THEN
          DO itrc=1,NT(ng)
            is=isTvar(itrc)
            IF (ANY(CnormB(is,:))) THEN
              Lsame=.TRUE.
              EXIT
            END IF
          END DO
          IF (Lsame) THEN
            WRITE (stdout,20) TRIM(Text),                               &
     &                        '3D normalization factors at RHO-points'
            FLUSH (stdout)
          END IF
        END IF

        DO itrc=1,NT(ng)
          VnormRobc=Aspv
          is=isTvar(itrc)
          DO ibry=1,4
            IF (CnormB(is,ibry)) THEN
              DO ns=1,Nscale(ng)
                ms=ns
                iLap=Clap(is,ns,ng)
                VscaleB=0.0_r8
                B3davg=0.0_r8
                B3dsqr=0.0_r8
                B3d=0.0_r8
                IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
                  i=BOUNDS(ng)%edge(ibry,r2dvar)
                  IF (Lconvolve(ibry)) THEN
                    DO j=JstrT,JendT
                      val=on_r(i,j)
                      DO k=1,N(ng)
                        VscaleB(j,k)=1.0_r8/SQRT(val*Hz(i,j,k))
                      END DO
                    END DO
                  END IF
                ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
                  j=BOUNDS(ng)%edge(ibry,r2dvar)
                  IF (Lconvolve(ibry)) THEN
                    DO i=IstrT,IendT
                      val=om_r(i,j)
                      DO k=1,N(ng)
                        VscaleB(i,k)=1.0_r8/SQRT(val*Hz(i,j,k))
                      END DO
                    END DO
                  END IF
                END IF
                DO iter=1,Nrandom
                  IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
                    CALL white_noise3d_bry (ng, tile, iTLM, ibry,       &
     &                                      Rscheme(ng),                &
     &                                      JstrR, JendR,               &
     &                                      LBij, UBij, 1, N(ng),       &
     &                                      Bmin, Bmax, B3d)
                    DO k=1,N(ng)
                      DO j=JstrT,JendT
                        B3d(j,k)=B3d(j,k)*VscaleB(j,k)
                      END DO
                    END DO
                  ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
                    CALL white_noise3d_bry (ng, tile, iTLM, ibry,       &
     &                                      Rscheme(ng),                &
     &                                      IstrR, IendR,               &
     &                                      LBij, UBij, 1, N(ng),       &
     &                                      Bmin, Bmax, B3d)
                    DO k=1,N(ng)
                      DO i=IstrT,IendT
                        B3d(i,k)=B3d(i,k)*VscaleB(i,k)
                      END DO
                    END DO
                  END IF
!
!  First the horizontal convolution.
!
                  DO k=1,N(ng)
                    klev=k
                    CALL tl_CI_bry1d_tile (ng, tile,                    &
     &                                     LBi, UBi, LBj, UBj,          &
     &                                     LBij, UBij,                  &
     &                                     IminS, ImaxS, JminS, JmaxS,  &
     &                                     nstp, nnew, ifac,            &
     &                                     is, ibry,                    &
     &                                     ms, MiterCI(ns,ng), iLap,    &
     &                                     r3dvar, klev, itrc,          &
     &                                   emaxr1dz(:,klev,itrc,ibry,ns), &
     &                                   eminr1dz(:,klev,itrc,ibry,ns), &
     &                                     ci_rb, ci_pb, ci_qb, ci_xb,  &
#  ifdef SOLVE3D
     &                                     Hz,                          &
#  endif
     &                                     B3d(:,klev))
                  END DO
!
!  Now the vertical convolution.
!
                  CALL tl_conv_r3d_bry_tile (ng, tile, iTLM, ibry,      &
     &                                       BOUNDS(ng)%edge(:,r2dvar), &
     &                                       LBij, UBij,                &
     &                                       LBi, UBi, LBj, UBj,        &
     &                                       1, N(ng),                  &
     &                                       IminS, ImaxS, JminS, JmaxS,&
     &                                       NghostPoints,              &
     &                                       NHstepsB(ibry,is)/ifac,    &
     &                                       NVstepsB(ibry,is)/ifac,    &
     &                                       DTsizeHB(ibry,is),         &
     &                                       DTsizeVB(ibry,is),         &
     &                                       Kh, Kv,                    &
     &                                       pm, pn,                    &
     &                                       pmon_u, pnom_v,            &
#   ifdef MASKING
     &                                       rmask, umask, vmask,       &
#   endif
     &                                       Hz, z_r,                   &
     &                                       B3d)
                  IF (Lconvolve(ibry)) THEN
                    IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
                      DO k=1,N(ng)
                        DO j=Jstr,Jend
                          B3davg(j,k)=B3davg(j,k)+B3d(j,k)
                          B3dsqr(j,k)=B3dsqr(j,k)+B3d(j,k)*B3d(j,k)
                        END DO
                      END DO
                    ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
                      DO k=1,N(ng)
                        DO i=Istr,Iend
                          B3davg(i,k)=B3davg(i,k)+B3d(i,k)
                          B3dsqr(i,k)=B3dsqr(i,k)+B3d(i,k)*B3d(i,k)
                        END DO
                      END DO
                    END IF
                  END IF
                END DO
                IF (Lconvolve(ibry)) THEN
                  IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
                    DO k=1,N(ng)
                      DO j=Jstr,Jend
                        Bavg=FacAvg*B3davg(j,k)
                        Bsqr=FacAvg*B3dsqr(j,k)
#   ifdef MASKING
                        IF (rmask(i,j).gt.0.0_r8) THEN
                          VnormRobc(j,k,ibry,itrc)=                     &
     &                                       VnormRobc(j,k,ibry,itrc)+  &
     &                                       Bwgt(isTvar(itrc),ms,ng)/  &
     &                                       SQRT(Bsqr)

                        ELSE
                          VnormRobc(j,k,ibry,itrc)=0.0_r8
                        END IF
#   else
                        VnormRobc(j,k,ibry,itrc)=                       &
     &                                     VnormRobc(j,k,ibry,itrc)+    &
     &                                     Bwgt(isTvar(itrc),ms,ng)/    &
     &                                     SQRT(Bsqr)
#   endif
                      END DO
                    END DO
                  ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
                    DO k=1,N(ng)
                      DO i=Istr,Iend
                        Bavg=FacAvg*B3davg(i,k)
                        Bsqr=FacAvg*B3dsqr(i,k)
#   ifdef MASKING
                        IF (rmask(i,j).gt.0.0_r8) THEN
                          VnormRobc(i,k,ibry,itrc)=                     &
       &                                     VnormRobc(i,k,ibry,itrc)+  &
       &                                     Bwgt(isTvar(itrc),ms,ng)/  &
       &                                     SQRT(Bsqr)
                        ELSE
                          VnormRobc(i,k,ibry,itrc)=0.0_r8
                        END IF
#   else
                        VnormRobc(i,k,ibry,itrc)=                       &
       &                                   VnormRobc(i,k,ibry,itrc)+    &
       &                                   Bwgt(isTvar(itrc),ms,ng)/    &
       &                                   SQRT(Bsqr)
#   endif
                      END DO
                    END DO
                  END IF
                END IF
              END DO    ! end of ns-loop 
              CALL bc_r3d_bry_tile (ng, tile, ibry,                     &
     &                              LBij, UBij, 1, N(ng),               &
     &                              VnormRobc(:,:,ibry,itrc))
#  ifdef DISTRIBUTE
              Bwrk=RESHAPE(VnormRobc(:,:,ibry,itrc), (/IJKlen/))
              CALL mp_collect (ng, iTLM, IJKlen, Aspv, Bwrk)
              ic=0
              DO k=1,N(ng)
                DO ib=LBij,UBij
                  ic=ic+1
                  VnormRobc(ib,k,ibry,itrc)=Bwrk(ic)
                END DO
              END DO
#  endif
            END IF
          END DO
          IF (ANY(CnormB(is,:))) THEN
            ifield=idSbry(isTvar(itrc))

            SELECT CASE (NRM(ifile,ng)%IOtype)
              CASE (io_nf90)
                CALL netcdf_put_fvar (ng, iTLM, ncname,                 &
     &                                Vname(1,ifield),                  &
     &                                VnormRobc(LBij:,:,:,itrc),        &
     &                         start =(/1,1,1,NRM(ifile,ng)%Rindex/),   &
     &                         total = (/IJlen,N(ng),4,1/),             &
     &                         ncid = NRM(ifile,ng)%ncid,               &
     &                         varid = NRM(ifile,ng)%Vid(ifield))

#  if defined PIO_LIB && defined DISTRIBUTE
              CASE (io_pio)
                CALL pio_netcdf_put_fvar (ng, iTLM, ncname,             &
     &                                    Vname(1,ifield),              &
     &                                    VnormRobc(LBij:,:,:,itrc),    &
     &                         start =(/1,1,1,NRM(ifile,ng)%Rindex/),   &
     &                         total = (/IJlen,N(ng),4,1/),             &
     &                         pioFile = NRM(ifile,ng)%pioFile,         &
     &                         pioVar = NRM(ifile,ng)%pioVar(ifield)%vd)
#  endif
            END SELECT
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
# endif
!
!  Synchronize open boundaries normalization NetCDF file to disk to
!  allow other processes to access data immediately after it is
!  written.
!
        SELECT CASE (NRM(ifile,ng)%IOtype)
          CASE (io_nf90)
            CALL netcdf_sync (ng, iTLM, ncname,                         &
     &                        NRM(ifile,ng)%ncid)
# if defined PIO_LIB && defined DISTRIBUTE
          CASE (io_pio)
            CALL pio_netcdf_sync (ng, iTLM, ncname,                     &
     &                            NRM(ifile,ng)%pioFile)
# endif
        END SELECT
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#endif

#if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
!
!-----------------------------------------------------------------------
!  Compute surface forcing error covariance, B, normalization factors
!  using the randomization approach of Fisher and Courtier (1995).
!-----------------------------------------------------------------------
!
      ifile=4
      IF (LwrtNRM(ifile,ng)) THEN
        rec=1
        Text='surface forcing'
!
!  Set randomization summation factors.
!
        FacAvg=1.0_r8/REAL(Nrandom,r8)
        FacSqr=SQRT(REAL(Nrandom,r8))
!
!  Set time record index to write in normalization NetCDF file.
!
        ncname=NRM(ifile,ng)%name
        NRM(ifile,ng)%Rindex=NRM(ifile,ng)%Rindex+1
        NRM(ifile,ng)%Nrec=NRM(ifile,ng)%Nrec+1
!
!  Write out model time (s).
!
        SELECT CASE (NRM(ifile,ng)%IOtype)
          CASE (io_nf90)
            CALL netcdf_put_fvar (ng, iTLM, ncname,                     &
     &                            Vname(1,idtime), my_time,             &
     &                         start = (/NRM(ifile,ng)%Rindex/),        &
     &                         total = (/1/),                           &
     &                         ncid = NRM(ifile,ng)%ncid,               &
     &                         varid = NRM(ifile,ng)%Vid(idtime))

# if defined PIO_LIB && defined DISTRIBUTE
          CASE (io_pio)
            CALL pio_netcdf_put_fvar (ng, iTLM, ncname,                 &
     &                                Vname(1,idtime), my_time,         &
     &                         start = (/NRM(ifile,ng)%Rindex/),        &
     &                         total = (/1/),                           &
     &                         pioFile = NRM(ifile,ng)%pioFile,         &
     &                         pioVar = NRM(ifile,ng)%pioVar(idtime)%vd)
# endif
        END SELECT
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef ADJUST_WSTRESS
!
!  2D norm at U-stress points.
!
        IF (Cnorm(rec,isUstr)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Text),                               &
     &                      '2D normalization factors at U-points'
            FLUSH (stdout)
          END IF
          DO ns=1,Nscale(ng)
            ms=ns
            iLap=Clap(isUstr,ns,ng)
            DO j=JstrT,JendT
              DO i=IstrP,IendT
                A2davg(i,j)=0.0_r8
                A2dsqr(i,j)=0.0_r8
                Hscale(i,j)=1.0_r8/SQRT(om_u(i,j)*on_u(i,j))
              END DO
            END DO
            DO iter=1,Nrandom
              CALL white_noise2d (ng, iTLM, u2dvar, Rscheme(ng),        &
     &                            Istr, IendR, JstrR, JendR,            &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Amin, Amax, A2d)
              DO j=JstrT,JendT
                DO i=IstrP,IendT
                  A2d(i,j)=A2d(i,j)*Hscale(i,j)
                END DO
              END DO
              CALL tl_CI_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj, LBij, UBij,          &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         nstp, nnew, Lweak, ifac, isUstr,         &
     &                         ms, MiterCI(ns,ng), iLap,                &
     &                         u2dvar, 1, 1,                            &
     &                         emaxu2ds(:,ns), eminu2ds(:,ns),          &
     &                         ci_r, ci_p, ci_q, ci_x,                  &
#    ifdef READ_SCALES
     &                         bussclx(:,:,ns), busscly(:,:,ns),        &
#    endif
#    ifdef SOLVE3D
     &                         Hz,                                      &
#    endif
     &                         A2d)
              DO j=Jstr,Jend
                DO i=IstrU,Iend
                  A2davg(i,j)=A2davg(i,j)+A2d(i,j)
                  A2dsqr(i,j)=A2dsqr(i,j)+A2d(i,j)*A2d(i,j)
                END DO
              END DO
            END DO
            DO j=Jstr,Jend
              DO i=IstrU,Iend
                Aavg=FacAvg*A2davg(i,j)
                Asqr=FacAvg*A2dsqr(i,j)
#   ifdef MASKING
                IF (umask(i,j).gt.0.0_r8) THEN
                  HnormSUS(i,j)=HnormSUS(i,j)+                          &
     &                          Bwgt(isUstr,ms,ng)/SQRT(Asqr)
                ELSE
                  HnormSUS(i,j)=0.0_r8
                END IF
#   else
                HnormSUS(i,j)=HnormSUS(i,j)+                            &
     &                        Bwgt(isUstr,ms,ng)/SQRT(Asqr)
#   endif
              END DO
            END DO     ! end ns-loop
          END DO
          CALL dabc_u2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        HnormSUS)
#  ifdef DISTRIBUTE
          CALL mp_exchange2d (ng, tile, iTLM, 1,                        &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        HnormSUS)
#  endif
!
          SELECT CASE (NRM(ifile,ng)%IOtype)
            CASE (io_nf90)
              CALL wrt_norm2d_nf90 (ng, tile, iTLM, ncname,             &
     &                              LBi, UBi, LBj, UBj, idUsms,         &
     &                              NRM(ifile,ng)%ncid,                 &
     &                              NRM(ifile,ng)%Vid(idUsms),          &
     &                              NRM(ifile,ng)%Rindex,               &
#  ifdef MASKING
     &                              umask,                              &
#  endif
     &                              HnormSUS)

#  if defined PIO_LIB && defined DISTRIBUTE
            CASE (io_pio)
              IF (NRM(ifile,ng)%pioVar(idUsms)%dkind.eq.                &
     &            PIO_double) THEN
                ioDesc => ioDesc_dp_u2dvar(ng)
              ELSE
                ioDesc => ioDesc_sp_u2dvar(ng)
              END IF
              CALL wrt_norm2d_pio (ng, tile, iTLM, ncname,              &
     &                             LBi, UBi, LBj, UBj, idUsms,          &
     &                             NRM(ifile,ng)%pioFile,               &
     &                             NRM(ifile,ng)%pioVar(idUsms),        &
     &                             NRM(ifile,ng)%Rindex,                &
     &                             ioDesc,                              &
#   ifdef MASKING
     &                             umask,                               &
#   endif
     &                             HnormSUS)
#  endif
          END SELECT
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  2D norm at V-stress points.
!
        IF (Cnorm(rec,isVstr)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Text),                               &
     &                        '2D normalization factors at V-points'
            FLUSH (stdout)
          END IF
          DO ns=1,Nscale(ng)
            ms=ns
            iLap=Clap(isVstr,ns,ng)
            DO j=JstrP,JendT
              DO i=IstrT,IendT
                A2davg(i,j)=0.0_r8
                A2dsqr(i,j)=0.0_r8
                Hscale(i,j)=1.0_r8/SQRT(om_v(i,j)*on_v(i,j))
              END DO
            END DO
            DO iter=1,Nrandom
              CALL white_noise2d (ng, iTLM, v2dvar, Rscheme(ng),        &
     &                            IstrR, IendR, Jstr, JendR,            &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Amin, Amax, A2d)
              DO j=JstrP,JendT
                DO i=IstrT,IendT
                  A2d(i,j)=A2d(i,j)*Hscale(i,j)
                END DO
              END DO
              CALL tl_CI_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj, LBij, UBij,          &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         nstp, nnew, Lweak, ifac, isVstr,         &
     &                         ms, MiterCI(ns,ng), iLap,                &
     &                         v2dvar, 1, 1,                            &
     &                         emaxv2ds(:,ns), eminv2ds(:,ns),          &
     &                         ci_r, ci_p, ci_q, ci_x,                  &
#    ifdef READ_SCALES
     &                         bvssclx(:,:,ns), bvsscly(:,:,ns),        &
#    endif
#     ifdef SOLVE3D
     &                         Hz,                                      &
#     endif
     &                         A2d)
              DO j=JstrV,Jend
                DO i=Istr,Iend
                  A2davg(i,j)=A2davg(i,j)+A2d(i,j)
                  A2dsqr(i,j)=A2dsqr(i,j)+A2d(i,j)*A2d(i,j)
                END DO
              END DO
            END DO
            DO j=JstrV,Jend
              DO i=Istr,Iend
                Aavg=FacAvg*A2davg(i,j)
                Asqr=FacAvg*A2dsqr(i,j)
#   ifdef MASKING
                IF (vmask(i,j).gt.0.0_r8) THEN
                  HnormSVS(i,j)=HnormSVS(i,j)+                          &
     &                          Bwgt(isVstr,ms,ng)/SQRT(Asqr)
                ELSE
                  HnormSVS(i,j)=0.0_r8
                END IF
#   else
                HnormSVS(i,j)=HnormSVS(i,j)+                            &
     &                        Bwgt(isVstr,ms,ng)/SQRT(Asqr)
#   endif
              END DO
            END DO   ! end ns-loop
          END DO
          CALL dabc_v2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        HnormSVS)
#  ifdef DISTRIBUTE
          CALL mp_exchange2d (ng, tile, iTLM, 1,                        &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        HnormSVS)
#  endif
!
          SELECT CASE (NRM(ifile,ng)%IOtype)
            CASE (io_nf90)
              CALL wrt_norm2d_nf90 (ng, tile, iTLM, ncname,             &
     &                              LBi, UBi, LBj, UBj, idVsms,         &
     &                              NRM(ifile,ng)%ncid,                 &
     &                              NRM(ifile,ng)%Vid(idVsms),          &
     &                              NRM(ifile,ng)%Rindex,               &
#  ifdef MASKING
     &                              vmask,                              &
#  endif
     &                              HnormSVS)

#  if defined PIO_LIB && defined DISTRIBUTE
            CASE (io_pio)
              IF (NRM(ifile,ng)%pioVar(idVsms)%dkind.eq.                &
     &            PIO_double) THEN
                ioDesc => ioDesc_dp_v2dvar(ng)
              ELSE
                ioDesc => ioDesc_sp_v2dvar(ng)
              END IF
              CALL wrt_norm2d_pio (ng, tile, iTLM, ncname,              &
     &                             LBi, UBi, LBj, UBj, idVsms,          &
     &                             NRM(ifile,ng)%pioFile,               &
     &                             NRM(ifile,ng)%pioVar(idVsms),        &
     &                             NRM(ifile,ng)%Rindex,                &
     &                             ioDesc,                              &
#   ifdef MASKING
     &                             vmask,                               &
#   endif
     &                             HnormSVS)
#  endif
          END SELECT
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
!
!  2D norm at surface tracer flux points.
!
        IF (Master) THEN
          Lsame=.FALSE.
          DO itrc=1,NT(ng)
            IF (Lstflux(itrc,ng)) THEN
              is=isTsur(itrc)
              IF (Cnorm(rec,is)) Lsame=.TRUE.
            END IF
          END DO
          IF (Lsame) THEN
            WRITE (stdout,20) TRIM(Text),                               &
     &                        '2D normalization factors at RHO-points'
            FLUSH (stdout)
          END IF
        END IF
!
!  Check if the decorrelation scales for all the surface tracer fluxes
!  are different. If not, just compute the normalization factors for the
!  first tracer and assign the same value to the rest.  Recall that this
!  computation is very expensive.
!
        Ldiffer=.FALSE.
        DO ns=1,Nscale(ng)
          DO itrc=2,NT(ng)
            IF ((HdecayX(rec,isTsur(itrc  ),ns,ng).ne.                  &
     &           HdecayX(rec,isTsur(itrc-1),ns,ng)).or.                 &
     &          (HdecayY(rec,isTsur(itrc  ),ns,ng).ne.                  &
     &           HdecayY(rec,isTsur(itrc-1),ns,ng))) THEN
              Ldiffer=.TRUE.
            END IF
          END DO
          IF (.not.Ldiffer) THEN
            Lsame=.TRUE.
            UBt=1
          ELSE
            Lsame=.FALSE.
            UBt=NT(ng)
          END IF
        END DO
!
        DO j=JstrT,JendT
          DO i=IstrT,IendT
            Hscale(i,j)=1.0_r8/SQRT(om_r(i,j)*on_r(i,j))
          END DO
        END DO
        DO itrc=1,UBt
          IF (Lstflux(itrc,ng)) THEN
            is=isTsur(itrc)
            IF (Cnorm(rec,is)) THEN
              DO ns=1,Nscale(ng)
                ms=ns
                iLap=Clap(is,ns,ng)
                DO j=JstrT,JendT
                  DO i=IstrT,IendT
                    A2davg(i,j)=0.0_r8
                    A2dsqr(i,j)=0.0_r8
                  END DO
                END DO
                DO iter=1,Nrandom
                  CALL white_noise2d (ng, iTLM, r2dvar, Rscheme(ng),    &
     &                                IstrR, IendR, JstrR, JendR,       &
     &                                LBi, UBi, LBj, UBj,               &
     &                                Amin, Amax, A2d)
                  DO j=JstrT,JendT
                    DO i=IstrT,IendT
                      A2d(i,j)=A2d(i,j)*Hscale(i,j)
                    END DO
                  END DO
                  CALL tl_CI_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj, LBij, UBij,      &
     &                             IminS, ImaxS, JminS, JmaxS,          &
     &                             nstp, nnew, Lweak, ifac, is,         &
     &                             ms, MiterCI(ns,ng), iLap,            &
     &                             r2dvar, 1, 1,                        &
     &                             emaxr2ds(:,itrc,ns),                 &
     &                             eminr2ds(:,itrc,ns),                 &
     &                             ci_r, ci_p, ci_q, ci_x,              &
#   ifdef READ_SCALES
     &                             btfsclx(:,:,itrc,ns),                &
     &                             btfscly(:,:,itrc,ns),                &
#   endif
#   ifdef SOLVE3D
     &                             Hz,                                  &
#   endif
     &                             A2d)
                  DO j=Jstr,Jend
                    DO i=Istr,Iend
                      A2davg(i,j)=A2davg(i,j)+A2d(i,j)
                      A2dsqr(i,j)=A2dsqr(i,j)+A2d(i,j)*A2d(i,j)
                    END DO
                  END DO
                END DO
                DO j=Jstr,Jend
                  DO i=Istr,Iend
                    Aavg=FacAvg*A2davg(i,j)
                    Asqr=FacAvg*A2dsqr(i,j)
#   ifdef MASKING
                    IF (rmask(i,j).gt.0.0_r8) THEN
                      HnormSTF(i,j,itrc)=HnormSTF(i,j,itrc)+            &
     &                                   Bwgt(isTsur(itrc),ms,ng)/      &
     &                                   SQRT(Asqr)
                    ELSE
                      HnormSTF(i,j,itrc)=0.0_r8
                    END IF
#   else
                    HnormSTF(i,j,itrc)=HnormSTF(i,j,itrc)+              &
     &                                 Bwgt(isTsur(itrc),ms,ng)/        &
     &                                 SQRT(Asqr)
#   endif
                  END DO
                END DO
              END DO   ! end ns-loop
            END IF
          END IF
        END DO
        IF (Lsame) THEN
          DO itrc=2,NT(ng)
            IF (Lstflux(itrc,ng)) THEN
              DO j=Jstr,Jend
                DO i=Istr,Iend
                  HnormSTF(i,j,itrc)=HnormSTF(i,j,1)
                END DO
              END DO
            END IF
          END DO
        END IF
        DO itrc=1,NT(ng)
          IF (Lstflux(itrc,ng)) THEN
            is=isTsur(itrc)
            IF (Cnorm(rec,is)) THEN
              CALL dabc_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            HnormSTF(:,:,itrc))
#  ifdef DISTRIBUTE
              CALL mp_exchange2d (ng, tile, iTLM, 1,                    &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            NghostPoints,                         &
     &                            EWperiodic(ng), NSperiodic(ng),       &
     &                            HnormSTF(:,:,itrc))
#  endif
!
              SELECT CASE (NRM(ifile,ng)%IOtype)
                CASE (io_nf90)
                  CALL wrt_norm2d_nf90 (ng, tile, iTLM, ncname,         &
     &                                  LBi, UBi, LBj, UBj,             &
     &                                  idTsur(itrc),                   &
     &                              NRM(ifile,ng)%ncid,                 &
     &                              NRM(ifile,ng)%Vid(idTsur(itrc)),    &
     &                              NRM(ifile,ng)%Rindex,               &
#  ifdef MASKING
     &                                  rmask,                          &
#  endif
     &                                  HnormSTF(:,:,itrc))

#  if defined PIO_LIB && defined DISTRIBUTE
                CASE (io_pio)
                  IF (NRM(ifile,ng)%pioVar(idTsur(itrc))%dkind.eq.      &
     &                PIO_double) THEN
                    ioDesc => ioDesc_dp_r2dvar(ng)
                  ELSE
                    ioDesc => ioDesc_sp_r2dvar(ng)
                  END IF
                  CALL wrt_norm2d_pio (ng, tile, iTLM, ncname,          &
     &                                 LBi, UBi, LBj, UBj,              &
     &                                 idTsur(itrc),                    &
     &                              NRM(ifile,ng)%pioFile,              &
     &                              NRM(ifile,ng)%pioVar(idTsur(itrc)), &
     &                              NRM(ifile,ng)%Rindex,               &
     &                                 ioDesc,                          &
#   ifdef MASKING
     &                                 rmask,                           &
#   endif
     &                                 HnormSTF(:,:,itrc))
#  endif
              END SELECT
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END IF
          END IF
        END DO
# endif
      END IF
#endif
!
      IF (Master) THEN
        WRITE (stdout,30)
      END IF

 10   FORMAT (/,' Error Covariance Factors: Randomization Method',/)
 20   FORMAT (4x,'Computing',1x,a,1x,a)
 30   FORMAT (/)
!
      RETURN
      END SUBROUTINE randomization_tile
!
!***********************************************************************
      FUNCTION dot_prod2d (ng, tile, model, ctype,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     A1, A2) RESULT DotProd
!***********************************************************************
!
!  Imported variable declarations.
!
      integer,   intent(in   ) :: ng                  ! nested grid
      integer,   intent(in   ) :: tile                ! domain partition
      integer,   intent(in   ) :: model               ! kernel ID
      integer,   intent(in   ) :: ctype               ! C-grid type
      integer,   intent(in   ) :: LBi, UBi, LBj, UBj
      real (r8), intent(in   ) :: A1(LBi:,LBj:)
      real (r8), intent(in   ) :: A2(LBi:,LBj:)
!
!  Local variable declarations.
!
      integer   :: IstrT, IstrP, IendT, Imin, Imax
      integer   :: JstrT, JstrP, JendT, Jmin, Jmax
      integer   :: NSUB, i, j
      real (r8) :: DotProd
      real (r8) :: cff, my_DotProd
#ifdef DISTRIBUTE
      character (len=3) :: op_handle
#endif
!
!-----------------------------------------------------------------------
!  It computes the dot product between two 2D tiled arrays.
!-----------------------------------------------------------------------
!
!  Initialize.
!
      IstrT=BOUNDS(ng)%IstrT(tile)   ! tile computational range
      IstrP=BOUNDS(ng)%IstrP(tile)
      IendT=BOUNDS(ng)%IendT(tile)
      JstrP=BOUNDS(ng)%JstrP(tile)
      JendT=BOUNDS(ng)%JendT(tile)
!
      Imin=IstrT
      Imax=IendT
      Jmin=JstrT
      Jmax=JendT
      SELECT CASE (ctype)
        CASE (u2dvar, u3dvar)
          Imin=IstrP
        CASE (v2dvar, v3dvar)
          Jmin=JstrP
      END SELECT
!
!  Compute dot product between A1 and A2 arrays.
!
      my_DotProd=0.0_r8
      DO j=Jmin,Jmax
        DO i=Imin,Imax
          cff=A1(i,j)*A2(i,j)
          my_DotProd=my_DotProd+cff
        END DO
      END DO
!
!  Perform parallel global reduction operation.
!
#ifdef DISTRIBUTE
      NSUB=1                             ! distributed-memory
#else
      IF (DOMAIN(ng)%SouthWest_Corner(tile).and.                        &
     &    DOMAIN(ng)%NorthEast_Corner(tile)) THEN
        NSUB=1                           ! non-tiled application
      ELSE
        NSUB=NtileX(ng)*NtileE(ng)       ! tiled application
      END IF
#endif
!$OMP CRITICAL (DOT_PROD)
      IF (tile_count.eq.0) THEN
        DotProd=0.0_r8
      END IF
      DotProd=DotProd+my_DotProd
      tile_count=tile_count+1
      IF (tile_count.eq.NSUB) THEN
        tile_count=0
#ifdef DISTRIBUTE
        op_handle='SUM'
        CALL mp_reduce (ng, model, 1, DotProd, op_handle)
#endif
      END IF
!$OMP END CRITICAL (DOT_PROD)
!
      RETURN
      END SUBROUTINE dot_prod2d

#ifdef SOLVE3D
!
!***********************************************************************
      FUNCTION dot_prod3d (ng, tile, model, ctype,                      &
     &                     LBi, UBi, LBj, UBj, LBk, UBk,                &
     &                     A1, A2) RESULT DotProd
!***********************************************************************
!
!  Imported variable declarations.
!
      integer,   intent(in   ) :: ng                  ! nested grid
      integer,   intent(in   ) :: tile                ! domain partition
      integer,   intent(in   ) :: model               ! kernel ID
      integer,   intent(in   ) :: ctype               ! C-grid type
      integer,   intent(in   ) :: LBi, UBi, LBj, UBj, LBk, UBk
      real (r8), intent(in   ) :: A1(LBi:,LBj:,LBk:)
      real (r8), intent(in   ) :: A2(LBi:,LBj:,LBk:)
!
!  Local variable declarations.
!
      integer   :: IstrT, IstrP, IendT, Imin, Imax
      integer   :: JstrT, JstrP, JendT, Jmin, Jmax
      integer   :: NSUB, i, j, k
      real (r8) :: DotProd
      real (r8) :: cff, my_DotProd
# ifdef DISTRIBUTE
      character (len=3) :: op_handle
# endif
!
!-----------------------------------------------------------------------
!  It computes the dot product between two 2D tiled arrays.
!-----------------------------------------------------------------------
!
!  Initialize.
!
      IstrT=BOUNDS(ng)%IstrT(tile)   ! tile computational range
      IstrP=BOUNDS(ng)%IstrP(tile)
      IendT=BOUNDS(ng)%IendT(tile)
      JstrP=BOUNDS(ng)%JstrP(tile)
      JendT=BOUNDS(ng)%JendT(tile)
!
      Imin=IstrT
      Imax=IendT
      Jmin=JstrT
      Jmax=JendT
      SELECT CASE (ctype)
        CASE (u2dvar, u3dvar)
          Imin=IstrP
        CASE (v2dvar, v3dvar)
          Jmin=JstrP
      END SELECT
!
!  Compute dot product between A1 and A2 arrays.
!
      my_DotProd=0.0_r8
      DO k=LBk,UBk
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            cff=A1(i,j,k)*A2(i,j,k)
            my_DotProd=my_DotProd+cff
          END DO
        END DO
      END DO
!
!  Perform parallel global reduction operation.
!
# ifdef DISTRIBUTE
      NSUB=1                             ! distributed-memory
# else
      IF (DOMAIN(ng)%SouthWest_Corner(tile).and.                        &
     &    DOMAIN(ng)%NorthEast_Corner(tile)) THEN
        NSUB=1                           ! non-tiled application
      ELSE
        NSUB=NtileX(ng)*NtileE(ng)       ! tiled application
      END IF
# endif
!$OMP CRITICAL (DOT_PROD)
      IF (tile_count.eq.0) THEN
        DotProd=0.0_r8
      END IF
      DotProd=DotProd+my_DotProd
      tile_count=tile_count+1
      IF (tile_count.eq.NSUB) THEN
        tile_count=0
# ifdef DISTRIBUTE
        op_handle='SUM'
        CALL mp_reduce (ng, model, 1, DotProd, op_handle)
# endif
      END IF
!$OMP END CRITICAL (DOT_PROD)
!
      RETURN
      END SUBROUTINE dot_prod3d
#endif

!
!***********************************************************************
      SUBROUTINE wrt_norm2d_nf90 (ng, tile, model, ncname,              &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            ifield, ncid, ncvarid, tindex,        &
#ifdef MASKING
     &                            Amask,                                &
#endif
     &                            A)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: ifield, ncid, ncvarid, tindex

      character (len=*), intent(in) :: ncname
!
#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:,LBj:)
# endif
      real(r8), intent(in) :: A(LBi:,LBj:)
#else
# ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in) :: A(LBi:UBi,LBj:UBj)
#endif
!
!  Local variable declarations.
!
      integer :: gfactor, gtype, status
!
      real(dp) :: scale
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", wrt_norm2d_nf90"

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Write out requested 2D field into normalization NetCDF file. Since
!  the computation of normalization coefficients is a very expensive
!  computation, synchronize file to disk.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Set grid type factor to write full (gfactor=1) fields or water
!  points (gfactor=-1) fields only.
!
#if defined WRITE_WATER && defined MASKING
        gfactor=-1
#else
        gfactor=1
#endif
!
!  Write out 2D normalization field.
!
      gtype=gfactor*Iinfo(1,ifield,ng)
      scale=1.0_dp
      status=nf_fwrite2d(ng, model, ncid, ifield,                       &
     &                   ncvarid, tindex, gtype,                        &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#ifdef MASKING
     &                   Amask,                                         &
#endif
     &                   A)
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,ifield)), tindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Synchronize normalization NetCDF file to disk to allow other
!  processes to access data immediately after it is written.
!
      CALL netcdf_sync (ng, model, ncname, ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      IF (Master) WRITE (stdout,20) TRIM(Vname(1,ifield)), tindex
      FLUSH (stdout)
!
  10  FORMAT (/,' WRT_NORM2D_NF90 - error while writing variable: ',a,  &
     &        /,19x 'into normalization NetCDF file for time record: ', &
     &        i0)
  20  FORMAT (7x,'wrote  ',a, t21,'normalization factors into record ', &
     &        i0)
!
      RETURN
      END SUBROUTINE wrt_norm2d_nf90

#if defined PIO_LIB && defined DISTRIBUTE
!
!***********************************************************************
      SUBROUTINE wrt_norm2d_pio (ng, tile, model, ncname,               &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           ifield, pioFile, pioVar, tindex,       &
     &                           pioDesc,                               &
# ifdef MASKING
     &                           Amask,                                 &
# endif
     &                           A)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: ifield, tindex

      character (len=*), intent(in) :: ncname
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: A(LBi:,LBj:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: A(LBi:UBi,LBj:UBj)
# endif
!
      TYPE (File_desc_t), intent(inout) :: pioFile
      TYPE (IO_Desc_t),   intent(inout) :: pioDesc
      TYPE (My_VarDesc),  intent(inout) :: pioVar
!
!  Local variable declarations.
!
      integer :: status
!
      real(dp) :: scale
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", wrt_norm2d_pio"

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Write out requested 2D field into normalization NetCDF file. Since
!  the computation of normalization coefficients is a very expensive
!  computation, synchronize file to disk.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out 2D normalization field.
!
      scale=1.0_dp
      status=nf_fwrite2d(ng, model, pioFile, ifield,                    &
     &                   pioVar, tindex, pioDesc,                       &
     &                   LBi, UBi, LBj, UBj, scale,                     &
# ifdef MASKING
     &                   Amask,                                         &
# endif
     &                   A)
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,ifield)), tindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Synchronize normalization NetCDF file to disk to allow other
!  processes to access data immediately after it is written.
!
      CALL pio_netcdf_sync (ng, model, ncname, pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      IF (Master) WRITE (stdout,20) TRIM(Vname(1,ifield)), tindex
      FLUSH (stdout)
!
  10  FORMAT (/,' WRT_NORM2D_PIO - error while writing variable: ',a,   &
     &        /,18x 'into normalization NetCDF file for time record: ', &
     &        i0)
  20  FORMAT (7x,'wrote  ',a, t21,'normalization factors into record ', &
     &        i0)
!
      RETURN
      END SUBROUTINE wrt_norm2d_pio
#endif

!
!***********************************************************************
      SUBROUTINE wrt_norm3d_nf90 (ng, tile, model, ncname,              &
     &                            LBi, UBi, LBj, UBj, LBk, UBk,         &
     &                            ifield, ncid, ncvarid, tindex,        &
#ifdef MASKING
     &                            Amask,                                &
#endif
     &                            A)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
      integer, intent(in) :: ifield, ncid, ncvarid, tindex

      character (len=*), intent(in) :: ncname
!
#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:,LBj:)
# endif
      real(r8), intent(in) :: A(LBi:,LBj:,LBk:)
#else
# ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in) :: A(LBi:UBi,LBj:UBj,LBk:UBk)
#endif
!
!  Local variable declarations.
!
      integer :: gfactor, gtype, status
!
      real(dp) :: scale
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", wrt_norm3d_nf90"

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Write out requested 3D field into normalization NetCDF file. Since
!  the computation of normalization coefficients is a very expensive
!  computation, synchronize file to disk.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Set grid type factor to write full (gfactor=1) fields or water
!  points (gfactor=-1) fields only.
!
#if defined WRITE_WATER && defined MASKING
        gfactor=-1
#else
        gfactor=1
#endif
!
!  Write out 3D normalization field.
!
      gtype=gfactor*Iinfo(1,ifield,ng)
      scale=1.0_dp
      status=nf_fwrite3d(ng, model, ncid, ifield,                       &
     &                   ncvarid, tindex, gtype,                        &
     &                   LBi, UBi, LBj, UBj, LBk, UBk, scale,           &
#ifdef MASKING
     &                   Amask,                                         &
#endif
     &                   A)
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,ifield)), tindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Synchronize normalization NetCDF file to disk to allow other
!  processes to access data immediately after it is written.
!
      CALL netcdf_sync (ng, model, ncname, ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      IF (Master) WRITE (stdout,20) TRIM(Vname(1,ifield)), tindex
      FLUSH (stdout)
!
  10  FORMAT (/,' WRT_NORM3D_NF90 - error while writing variable: ',a,  &
     &        /,19x,'into normalization NetCDF file for time record: ', &
     &        i0)
  20  FORMAT (7x,'wrote  ',a, t21,'normalization factors into record ', &
     &        i0)
!
      RETURN
      END SUBROUTINE wrt_norm3d_nf90

#if defined PIO_LIB && defined DISTRIBUTE
!
!***********************************************************************
      SUBROUTINE wrt_norm3d_pio (ng, tile, model, ncname,               &
     &                           LBi, UBi, LBj, UBj, LBk, UBk,          &
     &                           ifield, pioFile, pioVar, tindex,       &
     &                           pioDesc,                               &
# ifdef MASKING
     &                           Amask,                                 &
# endif
     &                           A)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
      integer, intent(in) :: ifield, tindex

      character (len=*), intent(in) :: ncname
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: A(LBi:,LBj:,LBk:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: A(LBi:UBi,LBj:UBj,LBk:UBk)
# endif
!
      TYPE (File_desc_t), intent(inout) :: pioFile
      TYPE (IO_Desc_t),   intent(inout) :: pioDesc
      TYPE (My_VarDesc),  intent(inout) :: pioVar
!
!  Local variable declarations.
!
      integer :: status
!
      real(dp) :: scale
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", wrt_norm3d_pio"

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Write out requested 3D field into normalization NetCDF file. Since
!  the computation of normalization coefficients is a very expensive
!  computation, synchronize file to disk.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out 3D normalization field.
!
      scale=1.0_dp
      status=nf_fwrite3d(ng, model, pioFile, ifield,                    &
     &                   pioVar, tindex, pioDesc,                       &
     &                   LBi, UBi, LBj, UBj, LBk, UBk, scale,           &
# ifdef MASKING
     &                   Amask,                                         &
# endif
     &                   A)
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,ifield)), tindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Synchronize normalization NetCDF file to disk to allow other
!  processes to access data immediately after it is written.
!
      CALL pio_netcdf_sync (ng, model, ncname, pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      IF (Master) WRITE (stdout,20) TRIM(Vname(1,ifield)), tindex
      FLUSH (stdout)
!
  10  FORMAT (/,' WRT_NORM3D_PIO - error while writing variable: ',a,   &
     &        /,18x,'into normalization NetCDF file for time record: ', &
     &        i0)
  20  FORMAT (7x,'wrote  ',a, t21,'normalization factors into record ', &
     &        i0)
!
      RETURN
      END SUBROUTINE wrt_norm3d_pio
#endif
      END MODULE normalization_mod
