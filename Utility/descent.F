#include "cppdefs.h"
      MODULE descent_mod

#if defined S4DVAR || defined IS4DVAR

# ifdef EW_PERIODIC
#  define IU_RANGE Istr,Iend
#  define IV_RANGE Istr,Iend
# else
#  define IU_RANGE Istr,IendR
#  define IV_RANGE IstrR,IendR
# endif
# ifdef NS_PERIODIC
#  define JU_RANGE Jstr,Jend
#  define JV_RANGE Jstr,Jend
# else
#  define JU_RANGE JstrR,JendR
#  define JV_RANGE Jstr,JendR
# endif
!
!=================================================== Andrew M. Moore ===
!  Copyright (c) 2005 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine estimates the  "best" state initial conditions using   !
!  a descent algorithm. The scheme used to find the minimum function   !
!  is selected according to the values of parameter ICG:               !
!                                                                      !
!     ICG=0     Fletcher-Reeves scheme                                 !
!     ICG=1     Polak-Ribiere scheme                                   !
!                                                                      !
!  For 2D applications, the control vector is (ubar,vbar,zeta) at      !
!  initial time. For 3D applications, the control vector is            !
!  (u,v,Tracers,zeta) at initial time and (ubar,vbar) are computed     !
!  at initial time from (u,v).                                         !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     tile      Sub-domain partition.                                  !
!     lock      Lock key for parallel reduction operations.            !
!     Iter      Current iteration.                                     !
!     step      Conjugate direction step size (nondimensional).        !
!                                                                      !
!  References:                                                         !
!                                                                      !
!     Fletcher, R. and. C.M. Reeves, 1964: Function minimization       !
!       by conjugate gradients, Comput. J., 7, 149-154.                !
!                                                                      !
!     Polak, E and G. Ribiere, 1969:  Note sur la convergence de       !
!       methodes de directions conjugees, Rev. Fr. Inform. Rech.       !
!       Oper., 16-R1, 35-43.                                           !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: descent

      CONTAINS
!
!***********************************************************************
      SUBROUTINE descent (ng, tile, model, Iter, step)
!***********************************************************************
!
      USE mod_param
# ifdef SOLVE3D
      USE mod_coupling
# endif
      USE mod_grid
      USE mod_ocean
      USE mod_stepping
!
      integer, intent(in) :: ng, tile, model, Iter

      real(r8), intent(in) :: step

# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, model, 36)
# endif
      CALL descent_tile (ng, Istr, Iend, Jstr, Jend,                    &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   Lold(ng), Lnew(ng), Iter, step,                &
# ifdef MASKING
     &                   GRID(ng) % rmask,                              &
     &                   GRID(ng) % umask,                              &
     &                   GRID(ng) % vmask,                              &
# endif
# ifdef SOLVE3D 
#  if defined SEDIMENT && defined SED_MORPH
     &                   OCEAN(ng) %  bed,                              &
     &                   GRID(ng) % bed_thick0,                         &
#  endif
     &                   GRID(ng) % Hz,                                 &
     &                   GRID(ng) % h,                                  &
     &                   GRID(ng) % om_v,                               &
     &                   GRID(ng) % on_u,                               &
#  ifdef ICESHELF
     &                   GRID(ng) % zice,                               &
#  endif
     &                   GRID(ng) % z_r,                                &
     &                   GRID(ng) % z_w,                                &
     &                   COUPLING(ng) % Zt_avg1,                        &
#  ifndef S4DVAR
#   if defined SEDIMENT && defined SED_MORPH
     &                   OCEAN(ng) %  tl_bed,                           &
#   endif
     &                   GRID(ng) % tl_Hz,                              &
     &                   GRID(ng) % tl_h,                               &
     &                   GRID(ng) % tl_z_r,                             &
     &                   GRID(ng) % tl_z_w,                             &
     &                   COUPLING(ng) % tl_Zt_avg1,                     &
#  endif
# endif
     &                   OCEAN(ng) % ubar,                              &
     &                   OCEAN(ng) % vbar,                              &
     &                   OCEAN(ng) % zeta,                              &
# if defined S4DVAR
#  ifdef SOLVE3D
     &                   OCEAN(ng) % t,                                 &
     &                   OCEAN(ng) % u,                                 &
     &                   OCEAN(ng) % v,                                 &
#  endif
# else
#  ifdef SOLVE3D
     &                   OCEAN(ng) % u,                                 &
     &                   OCEAN(ng) % v,                                 &
     &                   OCEAN(ng) % tl_t,                              &
     &                   OCEAN(ng) % tl_u,                              &
     &                   OCEAN(ng) % tl_v,                              &
#  endif
     &                   OCEAN(ng) % tl_ubar,                           &
     &                   OCEAN(ng) % tl_vbar,                           &
     &                   OCEAN(ng) % tl_zeta,                           &
# endif
# ifdef SOLVE3D
     &                   OCEAN(ng) % b_t,                               &
     &                   OCEAN(ng) % b_u,                               &
     &                   OCEAN(ng) % b_v,                               &
# endif
     &                   OCEAN(ng) % b_ubar,                            &
     &                   OCEAN(ng) % b_vbar,                            &
     &                   OCEAN(ng) % b_zeta,                            &
# ifdef SOLVE3D
     &                   OCEAN(ng) % d_t,                               &
     &                   OCEAN(ng) % d_u,                               &
     &                   OCEAN(ng) % d_v,                               &
# endif
     &                   OCEAN(ng) % d_ubar,                            &
     &                   OCEAN(ng) % d_vbar,                            &
     &                   OCEAN(ng) % d_zeta,                            &
# ifdef SOLVE3D
     &                   OCEAN(ng) % ad_t,                              &
     &                   OCEAN(ng) % ad_u,                              &
     &                   OCEAN(ng) % ad_v,                              &
# endif
     &                   OCEAN(ng) % ad_ubar,                           &
     &                   OCEAN(ng) % ad_vbar,                           &
     &                   OCEAN(ng) % ad_zeta)
# ifdef PROFILE
      CALL wclock_on (ng, model, 36)
# endif
      RETURN
      END SUBROUTINE descent
!
!***********************************************************************
      SUBROUTINE descent_tile (ng, Istr, Iend, Jstr, Jend,              &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         Lold, Lnew, Iter, step,                  &
# ifdef MASKING
     &                         rmask, umask, vmask,                     &
# endif
# ifdef SOLVE3D 
#  if defined SEDIMENT && defined SED_MORPH
     &                         bed, bed_thick0,                         &
#  endif
     &                         Hz, h, om_v, on_u,                       &
#  ifdef ICESHELF
     &                         zice,                                    &
#  endif
     &                         z_r, z_w, Zt_avg1,                       &
#  ifndef S4DVAR
#   if defined SEDIMENT && defined SED_MORPH
     &                         tl_bed,                                  &
#   endif
     &                         tl_Hz, tl_h,                             &
     &                         tl_z_r, tl_z_w, tl_Zt_avg1,              &
#  endif
# endif
     &                         ubar, vbar, zeta,                        &
# if defined S4DVAR
#  ifdef SOLVE3D
     &                         t, u, v,                                 &
#  endif
# else
#  ifdef SOLVE3D
     &                         u, v,                                    &
     &                         tl_t, tl_u, tl_v,                        &
#  endif
     &                         tl_ubar, tl_vbar, tl_zeta,               &
# endif
# ifdef SOLVE3D
     &                         b_t, b_u, b_v,                           &
# endif
     &                         b_ubar, b_vbar, b_zeta,                  &
# ifdef SOLVE3D
     &                         d_t, d_u, d_v,                           &
# endif
     &                         d_ubar, d_vbar, d_zeta,                  &
# ifdef SOLVE3D
     &                         ad_t, ad_u, ad_v,                        &
# endif
     &                         ad_ubar, ad_vbar, ad_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_scalars
!
# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      USE exchange_2d_mod
#  ifdef SOLVE3D
      USE exchange_3d_mod
#  endif
# endif
# ifdef SOLVE3D
      USE set_depth_mod, ONLY : set_depth_tile
#  ifndef S4DVAR
      USE tl_set_depth_mod, ONLY : tl_set_depth_tile
#  endif
# endif
# if defined S4DVAR
      USE u2dbc_mod, ONLY : u2dbc_tile
      USE v2dbc_mod, ONLY : v2dbc_tile
      USE zetabc_mod, ONLY : zetabc_tile
#  ifdef SOLVE3D
      USE t3dbc_mod, ONLY : t3dbc_tile
      USE u3dbc_mod, ONLY : u3dbc_tile
      USE v3dbc_mod, ONLY : v3dbc_tile
#  endif
# else
      USE tl_u2dbc_mod, ONLY : tl_u2dbc_tile
      USE tl_v2dbc_mod, ONLY : tl_v2dbc_tile
      USE tl_zetabc_mod, ONLY : tl_zetabc_tile
#  ifdef SOLVE3D
      USE tl_t3dbc_mod, ONLY : tl_t3dbc_tile
      USE tl_u3dbc_mod, ONLY : tl_u3dbc_tile
      USE tl_v3dbc_mod, ONLY : tl_v3dbc_tile
#  endif
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Lold, Lnew, Iter

      real(r8), intent(in) :: step

# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
#  ifdef SOLVE3D 
#   if defined SEDIMENT && defined SED_MORPH
      real(r8), intent(in) :: bed(LBi:,LBj:,:,:)
      real(r8), intent(in) :: bed_thick0(LBi:,LBj:)
#   endif
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)
#   ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:,LBj:)
#   endif
#   if !defined S4DVAR && (defined SEDIMENT && defined SED_MORPH)
      real(r8), intent(in) :: tl_bed(LBi:,LBj:,:,:)
#   endif
#  endif
#  ifdef SOLVE3D
      real(r8), intent(in) :: ad_t(LBi:,LBj:,:,:,:)
      real(r8), intent(in) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: ad_v(LBi:,LBj:,:,:)
#  endif
      real(r8), intent(in) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: ad_vbar(LBi:,LBj:,:)
      real(r8), intent(in) :: ad_zeta(LBi:,LBj:,:)
#  ifdef SOLVE3D
      real(r8), intent(in) :: b_t(LBi:,LBj:,:,:)
      real(r8), intent(in) :: b_u(LBi:,LBj:,:)
      real(r8), intent(in) :: b_v(LBi:,LBj:,:)
#  endif
      real(r8), intent(in) :: b_ubar(LBi:,LBj:)
      real(r8), intent(in) :: b_vbar(LBi:,LBj:)
      real(r8), intent(in) :: b_zeta(LBi:,LBj:)
#  ifndef S4DVAR
      real(r8), intent(in) :: ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: vbar(LBi:,LBj:,:)
      real(r8), intent(in) :: zeta(LBi:,LBj:,:)
#   ifdef SOLVE3D
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
#   endif
#  endif
#  ifdef SOLVE3D
      real(r8), intent(inout) :: Hz(LBi:,LBj:,:)
      real(r8), intent(inout) :: Zt_avg1(LBi:,LBj:)
      real(r8), intent(inout) :: z_r(LBi:,LBj:,:)
      real(r8), intent(inout) :: z_w(LBi:,LBj:,0:)
#   ifndef S4DVAR
      real(r8), intent(inout) :: tl_Hz(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_h(LBi:,LBj:)
      real(r8), intent(inout) :: tl_Zt_avg1(LBi:,LBj:)
      real(r8), intent(inout) :: tl_z_r(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_z_w(LBi:,LBj:,0:)
#   endif
      real(r8), intent(inout) :: d_t(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: d_u(LBi:,LBj:,:)
      real(r8), intent(inout) :: d_v(LBi:,LBj:,:)
#  endif
      real(r8), intent(inout) :: d_ubar(LBi:,LBj:)
      real(r8), intent(inout) :: d_vbar(LBi:,LBj:)
      real(r8), intent(inout) :: d_zeta(LBi:,LBj:)

#  if defined S4DVAR
#   ifdef SOLVE3D
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: v(LBi:,LBj:,:,:)
#   endif
      real(r8), intent(inout) :: ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: vbar(LBi:,LBj:,:)
      real(r8), intent(inout) :: zeta(LBi:,LBj:,:)
#  else
#   ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_v(LBi:,LBj:,:,:)
#   endif
      real(r8), intent(inout) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_vbar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_zeta(LBi:,LBj:,:)
#  endif

# else

#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef SOLVE3D 
#   if defined SEDIMENT && defined SED_MORPH
      real(r8), intent(in) :: bed(LBi:UBi,LBj:UBj,Nbed,MBEDP)
      real(r8), intent(in) :: bed_thick0(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)
#   ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:UBi,LBj:UBj)
#   endif
#   if !defined S4DVAR && (defined SEDIMENT && defined SED_MORPH)
      real(r8), intent(in) :: tl_bed(LBi:UBi,LBj:UBj,Nbed,MBEDP)
#   endif
#  endif
#  ifdef SOLVE3D
      real(r8), intent(in) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(in) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
#  endif
      real(r8), intent(in) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: ad_vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: ad_zeta(LBi:UBi,LBj:UBj,3)
#  ifdef SOLVE3D
      real(r8), intent(in) :: b_t(LBi:UBi,LBj:UBj,N(ng),NT(ng))
      real(r8), intent(in) :: b_u(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: b_v(LBi:UBi,LBj:UBj,N(ng))
#  endif
      real(r8), intent(in) :: b_ubar(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: b_vbar(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: b_zeta(LBi:UBi,LBj:UBj)
#  ifndef S4DVAR
      real(r8), intent(in) :: ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: zeta(LBi:UBi,LBj:UBj,3)
#   ifdef SOLVE3D
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),2)
#   endif
#  endif
#  ifdef SOLVE3D
      real(r8), intent(inout) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: Zt_avg1(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
#   ifndef S4DVAR
      real(r8), intent(inout) :: tl_Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: tl_h(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: tl_Zt_avg1(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: tl_z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: tl_z_w(LBi:UBi,LBj:UBj,0:N(ng))
#   endif
      real(r8), intent(inout) :: d_t(LBi:UBi,LBj:UBj,N(ng),NT(ng))
      real(r8), intent(inout) :: d_u(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: d_v(LBi:UBi,LBj:UBj,N(ng))
#  endif
      real(r8), intent(inout) :: d_ubar(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: d_vbar(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: d_zeta(LBi:UBi,LBj:UBj)

#  if defined S4DVAR
#   ifdef SOLVE3D
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: v(LBi:UBi,LBj:UBj,N(ng),2)
#   endif
      real(r8), intent(inout) :: ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: zeta(LBi:UBi,LBj:UBj,3)
#  else
#   ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: tl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: tl_v(LBi:UBi,LBj:UBj,N(ng),2)
#   endif
      real(r8), intent(inout) :: tl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: tl_vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: tl_zeta(LBi:UBi,LBj:UBj,3)
#  endif
# endif
!
!  Local variable declarations.
!
# ifdef SOLVE3D
      logical :: FirstPass = .TRUE.
# endif
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: NSUB, i, j
# ifdef SOLVE3D
      integer :: itrc, k
# endif
      real(r8) :: BetaK, CGscheme
      real(r8) :: beta1, beta2, my_beta1, my_beta2
      real(r8) :: dot1, dot2,  my_dot1, my_dot2
      real(r8) :: cff, cff1, cff2, cff3
#  ifdef IS4DVAR
      real(r8) :: tl_cff
#  endif

# ifdef SOLVE3D
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N(ng)) :: CF
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N(ng)) :: DC
#  ifdef IS4DVAR
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N(ng)) :: tl_CF
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N(ng)) :: tl_DC
#  endif
# endif
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: work
!
# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  On first pass, compute the dot product between previous and new
!  adjoint solutions. Compute the Polak-Ribiere scaling factor (BetaK)
!  for the conjugate vectors.
!-----------------------------------------------------------------------
!
      IF ((Ipass.eq.1).and.(Iter.gt.1)) THEN
        my_dot1=0.0_r8
        my_dot2=0.0_r8
        my_beta1=0.0_r8
        my_beta2=0.0_r8
        CGscheme=REAL(ICG,r8)
!
!  2D state variables.
!
# ifndef SOLVE3D
        DO j=JstrR,JendR
          DO i=Istr,IendR
            cff1=ad_ubar(i,j,Lnew)*ad_ubar(i,j,Lnew)
            cff2=ad_ubar(i,j,Lold)*ad_ubar(i,j,Lold)
            cff3=1.0_r8/(b_ubar(i,j)*b_ubar(i,j))
#  ifdef MASKING
            cff3=cff3*umask(i,j)
#  endif
            my_dot1=my_dot1+cff1*cff3
            my_dot2=my_dot2+                                            &
     &              ad_ubar(i,j,Lold)*ad_ubar(i,j,Lnew)*cff3
            my_beta1=my_beta1+cff2*cff3
            my_beta2=my_beta2+                                          &
     &               ad_ubar(i,j,Lnew)*                                 &
     &               (-CGscheme*ad_ubar(i,j,Lold)+                      &
     &                ad_ubar(i,j,Lnew))*cff3
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            cff1=ad_vbar(i,j,Lnew)*ad_vbar(i,j,Lnew)
            cff2=ad_vbar(i,j,Lold)*ad_vbar(i,j,Lold)
            cff3=1.0_r8/(b_vbar(i,j)*b_vbar(i,j))
#  ifdef MASKING
            cff3=cff3*vmask(i,j)
#  endif
            my_dot1=my_dot1+cff1*cff3
            my_dot2=my_dot2+                                            &
     &              ad_vbar(i,j,Lold)*ad_vbar(i,j,Lnew)*cff3
            my_beta1=my_beta1+cff2*cff3
            my_beta2=my_beta2+                                          &
     &               ad_vbar(i,j,Lnew)*                                 &
     &               (-CGscheme*ad_vbar(i,j,Lold)+                      &
     &                ad_vbar(i,j,Lnew))*cff3
          END DO
        END DO
# endif
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            cff1=ad_zeta(i,j,Lnew)*ad_zeta(i,j,Lnew)
            cff2=ad_zeta(i,j,Lold)*ad_zeta(i,j,Lold)
            cff3=1.0_r8/(b_zeta(i,j)*b_zeta(i,j))
#  ifdef MASKING
            cff3=cff3*rmask(i,j)
#  endif
            my_dot1=my_dot1+cff1*cff3
            my_dot2=my_dot2+                                            &
     &              ad_zeta(i,j,Lold)*ad_zeta(i,j,Lnew)*cff3
            my_beta1=my_beta1+cff2*cff3
            my_beta2=my_beta2+                                          &
     &               ad_zeta(i,j,Lnew)*                                 &
     &               (-CGscheme*ad_zeta(i,j,Lold)+                      &
     &                ad_zeta(i,j,Lnew))*cff3
          END DO
        END DO
# ifdef SOLVE3D
!
!  3D state variables.
!
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=Istr,IendR
              cff1=ad_u(i,j,k,Lnew)*ad_u(i,j,k,Lnew)
              cff2=ad_u(i,j,k,Lold)*ad_u(i,j,k,Lold)
              cff3=1.0_r8/(b_u(i,j,k)*b_u(i,j,k))
#  ifdef MASKING
              cff3=cff3*umask(i,j)
#  endif
              my_dot1=my_dot1+cff1*cff3
              my_dot2=my_dot2+                                          &
     &                ad_u(i,j,k,Lold)*ad_u(i,j,k,Lnew)*cff3
              my_beta1=my_beta1+cff2*cff3
              my_beta2=my_beta2+                                        &
     &                 ad_u(i,j,k,Lnew)*                                &
     &                 (-CGscheme*ad_u(i,j,k,Lold)+                     &
     &                  ad_u(i,j,k,Lnew))*cff3
            END DO
          END DO
        END DO
        DO k=1,N(ng)
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              cff1=ad_v(i,j,k,Lnew)*ad_v(i,j,k,Lnew)
              cff2=ad_v(i,j,k,Lold)*ad_v(i,j,k,Lold)
              cff3=1.0_r8/(b_v(i,j,k)*b_v(i,j,k))
#  ifdef MASKING
              cff3=cff3*vmask(i,j)
#  endif
              my_dot1=my_dot1+cff1*cff3
              my_dot2=my_dot2+                                          &
     &                ad_v(i,j,k,Lold)*ad_v(i,j,k, Lnew)*cff3
              my_beta1=my_beta1+cff2*cff3
              my_beta2=my_beta2+                                        &
     &                 ad_v(i,j,k,Lnew)*                                &
     &                 (-CGscheme*ad_v(i,j,k,Lold)+                     &
     &                  ad_v(i,j,k,Lnew))*cff3
            END DO
          END DO
        END DO
        DO itrc=1,NT(ng)
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                cff1=ad_t(i,j,k,Lnew,itrc)*                             &
     &               ad_t(i,j,k,Lnew,itrc)
                cff2=ad_t(i,j,k,Lold,itrc)*                             &
     &               ad_t(i,j,k,Lold,itrc)
                cff3=1.0_r8/(b_t(i,j,k,itrc)*b_t(i,j,k,itrc))
#  ifdef MASKING
                cff3=cff3*rmask(i,j)
#  endif
                my_dot1=my_dot1+cff1*cff3
                my_dot2=my_dot2+                                        &
     &                  ad_t(i,j,k,Lold,itrc)*                          &
     &                  ad_t(i,j,k,Lnew,itrc)*cff3
                my_beta1=my_beta1+cff2*cff3
                my_beta2=my_beta2+                                      &
     &                   ad_t(i,j,k,Lnew,itrc)*                         &
     &                   (-CGscheme*                                    &
     &                    ad_t(i,j,k,Lold,itrc)+                        &
     &                    ad_t(i,j,k,Lnew,itrc))*cff3
              END DO
            END DO
          END DO
        END DO
# endif
!
!  Perform parallel global reduction operations.
!
        IF (SOUTH_WEST_CORNER.and.                                      &
     &      NORTH_EAST_CORNER) THEN
          NSUB=1                         ! non-tiled application
        ELSE
          NSUB=NtileX(ng)*NtileE(ng)     ! tiled application
        END IF
!$OMP CRITICAL (TL_DOT)
        IF (tile_count.eq.0) THEN
          dot1=0.0_r8
          dot2=0.0_r8
          beta1=0.0_r8
          beta2=0.0_r8
        END IF
        dot1=dot1+my_dot1
        dot2=dot2+my_dot2
        beta1=beta1+my_beta1
        beta2=beta2+my_beta2
        tile_count=block_count+1
        IF (tile_count.eq.NSUB) THEN
          tile_count=0
#ifdef DISTRIBUTE
          buffer(1)=dot1
          buffer(2)=dot2
          buffer(3)=beta1
          buffer(4)=beta2
          op_handle(1)='SUM'
          op_handle(2)='SUM'
          op_handle(3)='SUM'
          op_handle(4)='SUM'
          CALL mp_reduce (ng, model, 4, buffer, op_handle)
          dot1=buffer(1)
          dot2=buffer(2)
          beta1=buffer(3)
          beta2=buffer(4)
# endif
        END IF
!$OMP END CRITICAL (TL_DOT)
        dot1=CGtol*ABS(dot1)
        dot2=ABS(dot2)
        BetaK=beta2/beta1
!
!  Restart with steepest decent every "NiterSD" iterations.
!
        IF (MOD(Iter-IterSD,NiterSD).eq.0) THEN
          BetaK=0.0
          IterSD=Iter
        END IF
!
!  Perform conjugacy test and perform steepest descent if necessary.
!
        IF (dot2.gt.dot1) THEN
          BetaK=0.0_r8
          IterSD=Iter
        END IF
!
!  If first pass and first iteration, use steepest descent algorithm.
!
      ELSE IF ((Ipass.eq.1).and.(Iter.eq.1)) THEN
        dot1=0.0_r8
        dot2=0.0_r8
        BetaK=0.0_r8
      END IF
!
!-----------------------------------------------------------------------
!  Calculate conjugate vectors (descent directions) and new initial
!  conditions. Notice that the conjugate vectors are only saved on the
!  second pass.
!-----------------------------------------------------------------------
!
      IF (Ipass.eq.1) THEN
        IF (Master) THEN
          PRINT 10, dot1, dot2, BetaK, IterSD
 10       FORMAT (/,' DESCENT - old state dot product, dot1 = ',        &
     &            1p,e15.8,                                             &
     &            /,11x,'new state dot product, dot2 = ',1p,e15.8,      &
     &            /,11x,'Polak-Ribiere factor, BetaK = ',1p,e15.8,      &
     &            /,11x,'Last steepest descent iteration = ',i5.5,/)
        END IF
!
!  First pass, 2D state variables.
!
# ifndef SOLVE3D
        DO j=JstrR,JendR
          DO i=Istr,IendR
            cff=(-ad_ubar(i,j,Lnew)/b_ubar(i,j)+BetaK*d_ubar(i,j))
#  ifdef MASKING
            cff=cff*umask(i,j)
#  endif
#  if defined S4DVAR
            ubar(i,j,Lnew)=ubar(i,j,Lold)+step*cff/b_ubar(i,j)
#  elif defined IS4DVAR
            tl_ubar(i,j,Lnew)=tl_ubar(i,j,Lold)+step*cff/b_ubar(i,j)
#  endif
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            cff=(-ad_vbar(i,j,Lnew)/b_vbar(i,j)+BetaK*d_vbar(i,j))
#  ifdef MASKING
            cff=cff*vmask(i,j)
#  endif
#  if defined S4DVAR
            vbar(i,j,Lnew)=vbar(i,j,Lold)+step*cff/b_vbar(i,j)
#  elif defined IS4DVAR
            tl_vbar(i,j,Lnew)=tl_vbar(i,j,Lold)+step*cff/b_vbar(i,j)
#  endif
          END DO
        END DO
# endif
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            cff=(-ad_zeta(i,j,Lnew)/b_zeta(i,j)+BetaK*d_zeta(i,j))
# ifdef MASKING
            cff=cff*rmask(i,j)
# endif          
# if defined S4DVAR
            zeta(i,j,Lnew)=zeta(i,j,Lold)+step*cff/b_zeta(i,j)
# elif defined IS4DVAR
            tl_zeta(i,j,Lnew)=tl_zeta(i,j,Lold)+step*cff/b_zeta(i,j)
# endif
          END DO
        END DO
# ifdef SOLVE3D
!
!  First pass, 3D state variables.
!
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=Istr,IendR
              cff=(-ad_u(i,j,k,Lnew)/b_u(i,j,k)+BetaK*d_u(i,j,k))
#  ifdef MASKING
              cff=cff*umask(i,j)
#  endif
#  if defined S4DVAR
              u(i,j,k,Lnew)=u(i,j,k,Lold)+step*cff/b_u(i,j,k)
#  elif defined IS4DVAR
              tl_u(i,j,k,Lnew)=tl_u(i,j,k,Lold)+step*cff/b_u(i,j,k)
#  endif
            END DO
          END DO
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              cff=(-ad_v(i,j,k,Lnew)/b_u(i,j,k)+BetaK*d_v(i,j,k))
#  ifdef MASKING
              cff=cff*vmask(i,j)
#  endif
#  if defined S4DVAR 
              v(i,j,k,Lnew)=v(i,j,k,Lold)+step*cff/b_u(i,j,k)
#  elif defined IS4DVAR
              tl_v(i,j,k,Lnew)=tl_v(i,j,k,Lold)+step*cff/b_u(i,j,k)
#  endif
            END DO
          END DO
        END DO
!
        DO itrc=1,NT(ng)
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                cff=(-ad_t(i,j,k,Lnew,itrc)/b_t(i,j,k,itrc)+            &
     &               BetaK*d_t(i,j,k,itrc))
#  ifdef MASKING
                cff=cff*rmask(i,j)
#  endif          
#  if defined S4DVAR
                t(i,j,k,Lnew,itrc)=t(i,j,k,Lold,itrc)+                  &
     &                             step*cff/b_t(i,j,k,itrc)
#  elif defined IS4DVAR
                tl_t(i,j,k,Lnew,itrc)=tl_t(i,j,k,Lold,itrc)+            &
     &                                step*cff/b_t(i,j,k,itrc)
#  endif
              END DO
            END DO          
          END DO
        END DO
# endif
      ELSE IF (Ipass.eq.2 ) THEN
!
!  Second pass, 2D state variables.
!
# ifndef SOLVE3D
        DO j=JstrR,JendR
          DO i=Istr,IendR
            cff=(-ad_ubar(i,j,Lnew)/b_ubar(i,j)+BetaK*d_ubar(i,j))
#  ifdef MASKING
            cff=cff*umask(i,j)
#  endif
            d_ubar(i,j)=cff
#  if defined S4DVAR
            ubar(i,j,Lnew)=ubar(i,j,Lold)+step*cff/b_ubar(i,j)
#  elif defined IS4DVAR
            tl_ubar(i,j,Lnew)=tl_ubar(i,j,Lold)+step*cff/b_ubar(i,j)
#  endif
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            cff=(-ad_vbar(i,j,Lnew)/b_vbar(i,j)+BetaK*d_vbar(i,j))
#  ifdef MASKING
            cff=cff*vmask(i,j)
#  endif
            d_vbar(i,j)=cff
#  if defined S4DVAR
            vbar(i,j,Lnew)=vbar(i,j,Lold)+step*cff/b_vbar(i,j)
#  elif defined IS4DVAR
            tl_vbar(i,j,Lnew)=tl_vbar(i,j,Lold)+step*cff/b_vbar(i,j)
#  endif
          END DO
        END DO
# endif
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            cff=(-ad_zeta(i,j,Lnew)/b_zeta(i,j)+BetaK*d_zeta(i,j))
# ifdef MASKING
            cff=cff*rmask(i,j)
# endif          
            d_zeta(i,j)=cff
# if defined S4DVAR
            zeta(i,j,Lnew)=zeta(i,j,Lold)+step*cff/b_zeta(i,j)
# elif defined IS4DVAR
            tl_zeta(i,j,Lnew)=tl_zeta(i,j,Lold)+step*cff/b_zeta(i,j)
# endif
          END DO
        END DO
# ifdef SOLVE3D
!
!  Second pass, 3D state variables.
!
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=Istr,IendR
              cff=(-ad_u(i,j,k,Lnew)/b_u(i,j,k)+BetaK*d_u(i,j,k))
#  ifdef MASKING
              cff=cff*umask(i,j)
#  endif
              d_u(i,j,k)=cff
#  if defined S4DVAR
              u(i,j,k,Lnew)=u(i,j,k,Lold)+step*cff/b_u(i,j,k)
#  elif defined IS4DVAR
              tl_u(i,j,k,Lnew)=tl_u(i,j,k,Lold)+step*cff/b_u(i,j,k)
#  endif
            END DO
          END DO
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              cff=(-ad_v(i,j,k,Lnew)/b_u(i,j,k)+BetaK*d_v(i,j,k))
#  ifdef MASKING
              cff=cff*vmask(i,j)
#  endif
              d_v(i,j,k)=cff
#  if defined S4DVAR
              v(i,j,k,Lnew)=v(i,j,k,Lold)+step*cff/b_u(i,j,k)
#  elif defined IS4DVAR
              tl_v(i,j,k,Lnew)=tl_v(i,j,k,Lold)+step*cff/b_u(i,j,k)
#  endif
            END DO
          END DO
        END DO
!
        DO itrc=1,NT(ng)
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                cff=(-ad_t(i,j,k,Lnew,itrc)/b_t(i,j,k,itrc)+            &
     &               BetaK*d_t(i,j,k,itrc))
#  ifdef MASKING
                cff=cff*rmask(i,j)
#  endif          
                d_t(i,j,k,itrc)=cff
#  if defined S4DVAR
                t(i,j,k,Lnew,itrc)=t(i,j,k,Lold,itrc)+                  &
     &                             step*cff/b_t(i,j,k,itrc)
#  elif defined IS4DVAR
                tl_t(i,j,k,Lnew,itrc)=tl_t(i,j,k,Lold,itrc)+            &
     &                                step*cff/b_t(i,j,k,itrc)
#  endif
              END DO
            END DO          
          END DO
        END DO
# endif
      END IF     
!
!-----------------------------------------------------------------------
!  Apply lateral boundary conditions to 2D fields.
!-----------------------------------------------------------------------
!
# if defined S4DVAR
      CALL zetabc_tile (ng, Istr, Iend, Jstr, Jend,                     &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Lnew, Lnew, Lnew,                               &
     &                  zeta)
#  ifndef SOLVE3D
      CALL u2dbc_tile (ng, Istr, Iend, Jstr, Jend,                      &
     &                 LBi, UBi, LBj, UBj,                              &
     &                 Lnew, Lnew, Lnew,                                &
     &                 ubar, vbar, zeta)
      CALL v2dbc_tile (ng, Istr, Iend, Jstr, Jend,                      &
     &                 LBi, UBi, LBj, UBj,                              &
     &                 Lnew, Lnew, Lnew,                                &
     &                 ubar, vbar, zeta)
#  endif
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      CALL exchange_r2d_tile (ng, model, Istr, Iend, Jstr, Jend,        &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        zeta(:,:,Lnew))
#   ifndef SOLVE3D
      CALL exchange_u2d_tile (ng, model, Istr, Iend, Jstr, Jend,        &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        ubar(:,:,Lnew))
      CALL exchange_v2d_tile (ng, model, Istr, Iend, Jstr, Jend,        &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        vbar(:,:,Lnew))
#   endif
#  endif
# elif defined IS4DVAR
      CALL tl_zetabc_tile (ng, Istr, Iend, Jstr, Jend,                  &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     Lnew, Lnew, Lnew,                            &
     &                     zeta, tl_zeta)
#  ifndef SOLVE3D
      CALL tl_u2dbc_tile (ng, Istr, Iend, Jstr, Jend,                   &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    Lnew, Lnew, Lnew,                             &
     &                    ubar, vbar, zeta,                             &
     &                    tl_ubar, tl_vbar, tl_zeta)
      CALL tl_v2dbc_tile (ng, Istr, Iend, Jstr, Jend,                   &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    Lnew, Lnew, Lnew,                             &
     &                    ubar, vbar, zeta,                             &
     &                    tl_ubar, tl_vbar, tl_zeta)
#  endif
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      CALL exchange_r2d_tile (ng, model, Istr, Iend, Jstr, Jend,        &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        tl_zeta(:,:,Lnew))
#   ifndef SOLVE3D
      CALL exchange_u2d_tile (ng, model, Istr, Iend, Jstr, Jend,        &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        tl_ubar(:,:,Lnew))
      CALL exchange_v2d_tile (ng, model, Istr, Iend, Jstr, Jend,        &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        tl_vbar(:,:,Lnew))
#   endif
#  endif
# endif
# ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Compute new depths and thicknesses.
!-----------------------------------------------------------------------
!
      CALL set_depth_tile (ng, Istr, Iend, Jstr, Jend,                  &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     Lnew, FirstPass,                             &
     &                     zeta, h,                                     &
#  ifdef ICESHELF
     &                     zice,                                        &
#  endif
#  if defined SEDIMENT && defined SED_MORPH
     &                     bed, bed_thick0,                             &
#  endif
     &                     Zt_avg1,                                     &
     &                     Hz, z_r, z_w)
#  if defined IS4DVAR
      CALL tl_set_depth_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Lnew, FirstPass,                          &
     &                        tl_zeta,                                  &
     &                        h, tl_h,                                  &
#   ifdef ICESHELF
     &                        zice,                                     &
#   endif
#   if defined SEDIMENT && defined SED_MORPH
     &                        tl_bed, bed_thick0,                       &
#   endif
     &                        Zt_avg1, tl_Zt_avg1,                      &
     &                        tl_Hz, tl_z_r, tl_z_w)
#  endif
!
!-----------------------------------------------------------------------
!  Now apply lateral boundary conditions to corrected 3D fields.
!-----------------------------------------------------------------------
!
#  if defined S4DVAR
      CALL u3dbc_tile (ng, Istr, Iend, Jstr, Jend,                      &
     &                 LBi, UBi, LBj, UBj, N(ng),                       &
     &                 Lnew, Lnew, u)
      CALL v3dbc_tile (ng, Istr, Iend, Jstr, Jend,                      &
     &                 LBi, UBi, LBj, UBj, N(ng),                       &
     &                 Lnew, Lnew, v)
      DO itrc=1,NT(ng)
        CALL t3dbc_tile (ng, Istr, Iend, Jstr, Jend, itrc,              &
     &                   LBi, UBi, LBj, UBj, N(ng), NT(ng),             &
     &                   Lnew, Lnew, t)
      END DO
#   if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      CALL exchange_u3d_tile (ng, model, Istr, Iend, Jstr, Jend,        &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        NghostPoints,                             &
     &                        u(:,:,:,Lnew))
      CALL exchange_v3d_tile (ng, model, Istr, Iend, Jstr, Jend,        &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        NghostPoints,                             &
     &                        v(:,:,:,Lnew))
      DO itrc=1,NT(ng)
        CALL exchange_r3d_tile (ng, model, Istr, Iend, Jstr, Jend,      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          t(:,:,:,Lnew,itrc))
      END DO
#   endif
#  elif defined IS4DVAR
      CALL tl_u3dbc_tile (ng, Istr, Iend, Jstr, Jend,                   &
     &                    LBi, UBi, LBj, UBj, N(ng),                    &
     &                    Lnew, Lnew, tl_u)
      CALL tl_v3dbc_tile (ng, Istr, Iend, Jstr, Jend,                   &
     &                    LBi, UBi, LBj, UBj, N(ng),                    &
     &                    Lnew, Lnew, tl_v)
      DO itrc=1,NT(ng)
        CALL tl_t3dbc_tile (ng, Istr, Iend, Jstr, Jend, itrc,           &
     &                      LBi, UBi, LBj, UBj, N(ng), NT(ng),          &
     &                      Lnew, Lnew, tl_t)
      END DO
#   if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      CALL exchange_u3d_tile (ng, model, Istr, Iend, Jstr, Jend,        &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        NghostPoints,                             &
     &                        tl_u(:,:,:,Lnew))
      CALL exchange_v3d_tile (ng, model, Istr, Iend, Jstr, Jend,        &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        NghostPoints,                             &
     &                        tl_v(:,:,:,Lnew))
      DO itrc=1,NT(ng)
        CALL exchange_r3d_tile (ng, Istr, Iend, Jstr, Jend,             &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          tl_t(:,:,:,Lnew,itrc))
      END DO
#   endif
#  endif
!
!-----------------------------------------------------------------------
!  Compute vertically integrated momentum (ubar,vbar) from corrected
!  (u,v).
!-----------------------------------------------------------------------
!
!  Couple velocity component in the XI-direction.
!
      DO j=JU_RANGE
        DO i=IU_RANGE
          DC(i,0)=0.0_r8
          CF(i,0)=0.0_r8
#  ifdef IS4DVAR
          tl_DC(i,0)=0.0_r8
          tl_CF(i,0)=0.0_r8
#  endif
        END DO
!
!  Compute thicknesses of U-boxes DC(i,1:N), total depth of the water
!  column DC(i,0), and incorrect vertical mean CF(i,0).
!
        DO k=1,N(ng)
          DO i=IU_RANGE
            DC(i,k)=0.5_r8*(Hz(i,j,k)+Hz(i-1,j,k))*on_u(i,j)
            DC(i,0)=DC(i,0)+DC(i,k)
            CF(i,0)=CF(i,0)+DC(i,k)*u(i,j,k,Lnew)
#  if defined IS4DVAR
            tl_DC(i,k)=0.5_r8*(tl_Hz(i,j,k)+tl_Hz(i-1,j,k))*on_u(i,j)
            tl_DC(i,0)=tl_DC(i,0)+tl_DC(i,k)
            tl_CF(i,0)=tl_CF(i,0)+tl_DC(i,k)*u(i,j,k,Lnew)+             &
     &                            DC(i,k)*tl_u(i,j,k,Lnew)
#  endif
          END DO
        END DO
        DO i=IU_RANGE
          cff=1.0_r8/DC(i,0)
#  if defined S4DVAR
          ubar(i,j,Lnew)=CF(i,0)*cff
#  elif defined IS4DVAR
          tl_cff=-cff*cff*tl_DC(i,0)
          tl_ubar(i,j,Lnew)=tl_CF(i,0)*cff+CF(i,0)*tl_cff
#  endif
        END DO
!
!  Couple velocity component in the ETA-direction.
!
        IF (j.ge.Jstr) THEN
          DO i=IV_RANGE
            DC(i,0)=0.0_r8
            CF(i,0)=0.0_r8
#  ifdef IS4DVAR
            tl_DC(i,0)=0.0_r8
            tl_CF(i,0)=0.0_r8
#  endif
          END DO
!
!  Compute thicknesses of V-boxes DC(i,1:N), total depth of the water
!  column DC(i,0), and incorrect vertical mean CF(i,0).
!
          DO k=1,N(ng)
            DO i=IV_RANGE
              DC(i,k)=0.5_r8*(Hz(i,j,k)+Hz(i,j-1,k))*om_v(i,j)
              DC(i,0)=DC(i,0)+DC(i,k)
              CF(i,0)=CF(i,0)+DC(i,k)*v(i,j,k,Lnew)
#  if defined IS4DVAR
              tl_DC(i,k)=0.5_r8*(tl_Hz(i,j,k)+tl_Hz(i,j-1,k))*om_v(i,j)
              tl_DC(i,0)=tl_DC(i,0)+tl_DC(i,k)
              tl_CF(i,0)=tl_CF(i,0)+tl_DC(i,k)*v(i,j,k,Lnew)+           &
     &                              DC(i,k)*tl_v(i,j,k,Lnew)
#  endif
            END DO
          END DO
          DO i=IV_RANGE
            cff=1.0_r8/DC(i,0)
#  if defined S4DVAR
            vbar(i,j,Lnew)=CF(i,0)*cff
#  elif defined IS4DVAR
            tl_cff=-cff*cff*tl_DC(i,0)
            tl_vbar(i,j,Lnew)=tl_CF(i,0)*cff+CF(i,0)*tl_cff
#  endif
          END DO
        END IF   
      END DO
!
!  Apply boundary conditions.
!
#  if defined S4DVAR
      CALL u2dbc_tile (ng, Istr, Iend, Jstr, Jend,                      &
     &                 LBi, UBi, LBj, UBj,                              &
     &                 Lnew, Lnew, Lnew,                                &
     &                 ubar, vbar, zeta)
      CALL v2dbc_tile (ng, Istr, Iend, Jstr, Jend,                      &
     &                 LBi, UBi, LBj, UBj,                              &
     &                 Lnew, Lnew, Lnew,                                &
     &                 ubar, vbar, zeta)
#   if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      CALL exchange_u2d_tile (ng, model, Istr, Iend, Jstr, Jend,        &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        ubar(:,:,Lnew))
      CALL exchange_v2d_tile (ng, model, Istr, Iend, Jstr, Jend,        &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        vbar(:,:,Lnew))
#   endif
#  elif defined IS4DVAR
      CALL tl_u2dbc_tile (ng, Istr, Iend, Jstr, Jend,                   &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    Lnew, Lnew, Lnew,                             &
     &                    ubar, vbar, zeta,                             &
     &                    tl_ubar, tl_vbar, tl_zeta)
      CALL tl_v2dbc_tile (ng, Istr, Iend, Jstr, Jend,                   &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    Lnew, Lnew, Lnew,                             &
     &                    ubar, vbar, zeta,                             &
     &                    tl_ubar, tl_vbar, tl_zeta)
#   if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      CALL exchange_u2d_tile (ng, model, Istr, Iend, Jstr, Jend,        &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        tl_ubar(:,:,Lnew))
      CALL exchange_v2d_tile (ng, model, Istr, Iend, Jstr, Jend,        &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        tl_vbar(:,:,Lnew))
#   endif
#  endif
# endif
# undef IU_RANGE
# undef IV_RANGE
# undef JU_RANGE
# undef JV_RANGE 
      RETURN
      END SUBROUTINE descent_tile
#endif
      END MODULE descent_mod
