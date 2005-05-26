#include "cppdefs.h"
      MODULE ini_adjust_mod

#if defined GRADIENT_CHECK || defined IS4DVAR || defined TLM_CHECK || \
    defined PROPAGATOR
     

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
!  Copyright (c) 2005 ROMS/TOMS Adjoint Group                          !
!================================================== Hernan G. Arango ===
!                                                                      !
!  These routines adjust and perturb state initial conditions.         !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
# ifdef IS4DVAR
      PUBLIC :: ini_adjust
# endif
# ifdef TLM_CHECK
      PUBLIC :: ini_perturb
# endif
# if defined GRADIENT_CHECK || defined TLM_CHECK
      PUBLIC :: tl_ini_perturb
# endif
# if defined PROPAGATOR
      PUBLIC :: ad_ini_perturb
# endif

      CONTAINS

# ifdef IS4DVAR
      SUBROUTINE ini_adjust (ng, tile, Linp, Lout)
!
!=================================================== Andrew M. Moore ===
!  Copyright (c) 2005 ROMS/TOMS Adjoint Group                          !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine adds 4DVAR inner loop tangent linear increments to     !
!  nonlinear state initial conditions.                                 !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng        Nested grid number.                                    !
!     Linp      Tangent linear state time index to add.                !
!     Lout      Nonlinear state time index to update.                  !
!     tile      Domain partition.                                      !
!                                                                      !
!=======================================================================
!
!
      USE mod_param
#  ifdef SOLVE3D
      USE mod_coupling
#  endif
      USE mod_grid
      USE mod_ocean
!
      integer, intent(in) :: ng, tile, Linp, Lout

#  include "tile.h"
!
#  ifdef PROFILE
      CALL wclock_on (ng, iNLM, 7)
#  endif
      CALL ini_adjust_tile (ng, Istr, Iend, Jstr, Jend,                 &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      Linp, Lout,                                 &
#  ifdef MASKING
     &                      GRID(ng) % rmask,                           &
     &                      GRID(ng) % umask,                           &
     &                      GRID(ng) % vmask,                           &
#  endif
#  ifdef SOLVE3D
#   if defined SEDIMENT && defined SED_MORPH
     &                      OCEAN(ng) %  bed,                           &
     &                      GRID(ng) % bed_thick0,                      &
#   endif
     &                      GRID(ng) % Hz,                              &
     &                      GRID(ng) % h,                               &
     &                      GRID(ng) % om_v,                            &
     &                      GRID(ng) % on_u,                            &
#   ifdef ICESHELF
     &                      GRID(ng) % zice,                            &
#   endif
     &                      GRID(ng) % z_r,                             &
     &                      GRID(ng) % z_w,                             &
     &                      COUPLING(ng) % Zt_avg1,                     &
     &                      OCEAN(ng) % tl_t,                           &
     &                      OCEAN(ng) % tl_u,                           &
     &                      OCEAN(ng) % tl_v,                           &
     &                      OCEAN(ng) % t,                              &
     &                      OCEAN(ng) % u,                              &
     &                      OCEAN(ng) % v,                              &
#  endif
     &                      OCEAN(ng) % tl_ubar,                        &
     &                      OCEAN(ng) % tl_vbar,                        &
     &                      OCEAN(ng) % tl_zeta,                        &
     &                      OCEAN(ng) % ubar,                           &
     &                      OCEAN(ng) % vbar,                           &
     &                      OCEAN(ng) % zeta)
#  ifdef PROFILE
      CALL wclock_off (ng, iNLM, 7)
#  endif
      RETURN
      END SUBROUTINE ini_adjust

!
!***********************************************************************
      SUBROUTINE ini_adjust_tile (ng, Istr, Iend, Jstr, Jend,           &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Linp, Lout,                           &
#  ifdef MASKING
     &                            rmask, umask, vmask,                  &
#  endif
#  ifdef SOLVE3D
#   if defined SEDIMENT && defined SED_MORPH
     &                            bed, bed_thick0,                      &
#   endif
     &                            Hz, h, om_v, on_u,                    &
#   ifdef ICESHELF
     &                            zice,                                 &
#   endif
     &                            z_r, z_w, Zt_avg1,                    &
     &                            tl_t, tl_u, tl_v,                     &
     &                            t, u, v,                              &
#  endif
     &                            tl_ubar, tl_vbar, tl_zeta,            &
     &                            ubar, vbar, zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      USE exchange_2d_mod
#   ifdef SOLVE3D
      USE exchange_3d_mod
#   endif
#  endif
#  ifdef SOLVE3D
      USE set_depth_mod, ONLY : set_depth_tile
#  endif
      USE u2dbc_mod, ONLY : u2dbc_tile
      USE v2dbc_mod, ONLY : v2dbc_tile
      USE zetabc_mod, ONLY : zetabc_tile
#  ifdef SOLVE3D
      USE t3dbc_mod, ONLY : t3dbc_tile
      USE u3dbc_mod, ONLY : u3dbc_tile
      USE v3dbc_mod, ONLY : v3dbc_tile
#  endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Linp, Lout

#  ifdef ASSUMED_SHAPE
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#   endif
#   ifdef SOLVE3D
#    if defined SEDIMENT && defined SED_MORPH
      real(r8), intent(in) :: bed(LBi:,LBj:,:,:)
#    endif
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)
#    ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:,LBj:)
#    endif
      real(r8), intent(in) :: tl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(in) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: tl_v(LBi:,LBj:,:,:)
#   endif
      real(r8), intent(in) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: tl_vbar(LBi:,LBj:,:)
      real(r8), intent(in) :: tl_zeta(LBi:,LBj:,:)
#   ifdef SOLVE3D
#    if defined SEDIMENT && defined SED_MORPH
      real(r8), intent(inout) :: bed_thick0(LBi:,LBj:)
#    endif
      real(r8), intent(inout) :: Hz(LBi:,LBj:,:)
      real(r8), intent(inout) :: Zt_avg1(LBi:,LBj:)
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: v(LBi:,LBj:,:,:)
#   endif
      real(r8), intent(inout) :: ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: vbar(LBi:,LBj:,:)
      real(r8), intent(inout) :: zeta(LBi:,LBj:,:)
#   ifdef SOLVE3D
      real(r8), intent(out) :: z_r(LBi:,LBj:,:)
      real(r8), intent(out) :: z_w(LBi:,LBj:,0:)
#   endif
#  else
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#   endif
#   ifdef SOLVE3D
#    if defined SEDIMENT && defined SED_MORPH
      real(r8), intent(in) :: bed(LBi:UBi,LBj:UBj,Nbed,MBEDP)
#    endif
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)
#    ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:,LBj:)
#    endif
      real(r8), intent(in) :: tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(in) :: tl_u(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: tl_v(LBi:UBi,LBj:UBj,N(ng))
#   endif
      real(r8), intent(in) :: tl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: tl_vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: tl_zeta(LBi:UBi,LBj:UBj,3)
#   ifdef SOLVE3D
#    if defined SEDIMENT && defined SED_MORPH
      real(r8), intent(inout) :: bed_thick0(LBi:UBi,LBj:UBj)
#    endif
      real(r8), intent(inout) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: Zt_avg1(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: u(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: v(LBi:UBi,LBj:UBj,N(ng))
#   endif
      real(r8), intent(inout) :: ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: zeta(LBi:UBi,LBj:UBj,3)
#   ifdef SOLVE3D
      real(r8), intent(out) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(out) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
#   endif
#  endif
!
!  Local variable declarations.
!
      logical :: FirstPass = .FALSE.

      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: Tindex, i, j
#  ifdef SOLVE3D
      integer :: itrc, k

      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N(ng)) :: CF
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N(ng)) :: DC
#  endif

#  include "set_bounds.h"
!      
!-----------------------------------------------------------------------
!  Adjust initial conditions for 2D state variables.
!-----------------------------------------------------------------------
!
#  ifndef SOLVE3D
      DO j=JstrR,JendR
        DO i=Istr,IendR
          ubar(i,j,Lout)=ubar(i,j,Lout)+tl_ubar(i,j,Linp)
#   ifdef MASKING
          ubar(i,j,Lout)=ubar(i,j,Lout)*umask(i,j)
#   endif
        END DO
      END DO
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          vbar(i,j,Lout)=vbar(i,j,Lout)+tl_vbar(i,j,Linp)
#   ifdef MASKING
          vbar(i,j,Lout)=vbar(i,j,Lout)*vmask(i,j)
#   endif
        END DO
      END DO
#  endif
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          zeta(i,j,Lout)=zeta(i,j,Lout)+tl_zeta(i,j,Linp)
#  ifdef MASKING
          zeta(i,j,Lout)=zeta(i,j,Lout)*zeta(i,j)
#  endif
        END DO
      END DO

#  ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Adjust initial conditions for 2D state variables.
!-----------------------------------------------------------------------
!
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=Istr,IendR
            u(i,j,k,Lout)=u(i,j,k,Lout)+tl_u(i,j,k,Linp)
#   ifdef MASKING
            u(i,j,k,Lout)=u(i,j,k,Lout)*umask(i,j)
#   endif
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            v(i,j,k,Lout)=v(i,j,k,Lout)+tl_v(i,j,k,Linp)
#   ifdef MASKING
            v(i,j,k,Lout)=v(i,j,k,Lout)*vmask(i,j)
#   endif
          END DO
        END DO
      END DO
      DO itrc=1,NT(ng)
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              t(i,j,k,Lout,itrc)=t(i,j,k,Lout,itrc)+                    &
     &                           tl_t(i,j,k,Linp,itrc)
#   ifdef MASKING
              t(i,j,k,Lout,itrc)=t(i,j,k,Lout,itrc)*rmask(i,j)
#   endif
            END DO
          END DO          
        END DO
      END DO
#  endif
!
!-----------------------------------------------------------------------
!  Apply lateral boundary conditions to adjusted 2D fields.
!-----------------------------------------------------------------------
!
      CALL zetabc_tile (ng, Istr, Iend, Jstr, Jend,                     &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Lout, Lout, Lout,                               &
     &                  zeta)
#  ifndef SOLVE3D
      CALL u2dbc_tile (ng, Istr, Iend, Jstr, Jend,                      &
     &                 LBi, UBi, LBj, UBj,                              &
     &                 Lout, Lout, Lout,                                &
     &                 ubar, vbar, zeta)
      CALL v2dbc_tile (ng, Istr, Iend, Jstr, Jend,                      &
     &                 LBi, UBi, LBj, UBj,                              &
     &                 Lout, Lout, Lout,                                &
     &                 ubar, vbar, zeta)
#  endif
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      CALL exchange_r2d_tile (ng, model, Istr, Iend, Jstr, Jend,        &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        zeta(:,:,Lout))
#   ifndef SOLVE3D
      CALL exchange_u2d_tile (ng, model, Istr, Iend, Jstr, Jend,        &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        ubar(:,:,Lout))
      CALL exchange_v2d_tile (ng, model, Istr, Iend, Jstr, Jend,        &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        vbar(:,:,Lout))
#   endif
#  endif
#  ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Compute new depths and thicknesses.
!-----------------------------------------------------------------------
!
      CALL set_depth_tile (ng, Istr, Iend, Jstr, Jend,                  &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     Lout, FirstPass,                             &
     &                     zeta, h,                                     &
#   ifdef ICESHELF
     &                     zice,                                        &
#   endif
#   if defined SEDIMENT && defined SED_MORPH
     &                     bed, bed_thick0,                             &
#   endif
     &                     Zt_avg1,                                     &
     &                     Hz, z_r, z_w)
!
!-----------------------------------------------------------------------
!  Apply lateral boundary conditions to adjusted 3D fields.
!-----------------------------------------------------------------------
!
      CALL u3dbc_tile (ng, Istr, Iend, Jstr, Jend,                      &
     &                 LBi, UBi, LBj, UBj, N(ng),                       &
     &                 Lout, Lout, u)
      CALL v3dbc_tile (ng, Istr, Iend, Jstr, Jend,                      &
     &                 LBi, UBi, LBj, UBj, N(ng),                       &
     &                 Lout, Lout, v)
      DO itrc=1,NT(ng)
        CALL t3dbc_tile (ng, Istr, Iend, Jstr, Jend, itrc,              &
     &                   LBi, UBi, LBj, UBj, N(ng), NT(ng),             &
     &                   Lout, Lout, t)
      END DO
#   if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      CALL exchange_u3d_tile (ng, model, Istr, Iend, Jstr, Jend,        &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        NghostPoints,                             &
     &                        u(:,:,:,Lout))
      CALL exchange_v3d_tile (ng, model, Istr, Iend, Jstr, Jend,        &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        NghostPoints,                             &
     &                        v(:,:,:,Lout))
      DO itrc=1,NT(ng)
        CALL exchange_r3d_tile (ng, Istr, Iend, Jstr, Jend,             &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          t(:,:,:,Lout,itrc))
      END DO
#   endif
!
!-----------------------------------------------------------------------
!  Compute vertically integrated momentum (ubar,vbar) from adjusted
!  (u,v).
!-----------------------------------------------------------------------
!
!  Couple velocity component in the XI-direction.
!
      DO j=JU_RANGE
        DO i=IU_RANGE
          DC(i,0)=0.0_r8
          CF(i,0)=0.0_r8
        END DO
!
!  Compute thicknesses of U-boxes DC(i,1:N), total depth of the water
!  column DC(i,0), and incorrect vertical mean CF(i,0).
!
        DO k=1,N(ng)
          DO i=IU_RANGE
            DC(i,k)=0.5_r8*(Hz(i,j,k)+Hz(i-1,j,k))*on_u(i,j)
            DC(i,0)=DC(i,0)+DC(i,k)
            CF(i,0)=CF(i,0)+DC(i,k)*u(i,j,k,Lout)
          END DO
        END DO
        DO i=IU_RANGE
          ubar(i,j,Lout)=CF(i,0)/DC(i,0)
#   ifdef MASKING
          ubar(i,j,Lout)=ubar(i,j,Lout)*umask(i,j)
#   endif
        END DO
!
!  Couple velocity component in the ETA-direction.
!
        IF (j.ge.Jstr) THEN
          DO i=IV_RANGE
            DC(i,0)=0.0_r8
            CF(i,0)=0.0_r8
          END DO
!
!  Compute thicknesses of V-boxes DC(i,1:N), total depth of the water
!  column DC(i,0), and incorrect vertical mean CF(i,0).
!
          DO k=1,N(ng)
            DO i=IV_RANGE
              DC(i,k)=0.5_r8*(Hz(i,j,k)+Hz(i,j-1,k))*om_v(i,j)
              DC(i,0)=DC(i,0)+DC(i,k)
              CF(i,0)=CF(i,0)+DC(i,k)*v(i,j,k,Lout)
            END DO
          END DO
          DO i=IV_RANGE
            vbar(i,j,Lout)=CF(i,0)/DC(i,0)
#   ifdef MASKING
            vbar(i,j,Lout)=vbar(i,j,Lout)*vmask(i,j)
#   endif
          END DO
        END IF   
      END DO
!
!  Apply boundary conditions.
!
      CALL u2dbc_tile (ng, Istr, Iend, Jstr, Jend,                      &
     &                 LBi, UBi, LBj, UBj,                              &
     &                 Lout, Lout, Lout,                                &
     &                 ubar, vbar, zeta)
      CALL v2dbc_tile (ng, Istr, Iend, Jstr, Jend,                      &
     &                 LBi, UBi, LBj, UBj,                              &
     &                 Lout, Lout, Lout,                                &
     &                 ubar, vbar, zeta)
#   if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      CALL exchange_u2d_tile (ng, model, Istr, Iend, Jstr, Jend,        &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        ubar(:,:,Lout))
      CALL exchange_v2d_tile (ng, model, Istr, Iend, Jstr, Jend,        &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        vbar(:,:,Lout))
#   endif
#  endif
      RETURN
      END SUBROUTINE ini_adjust_tile
# endif       

# if defined TLM_CHECK
      SUBROUTINE ini_perturb (ng, tile, Linp, Lout)
!
!=======================================================================
!                                                                      !
!  Add a perturbation to nonlinear state variables according to the    !
!  outer and inner loop iterations. The added term is a function of    !
!  the steepest descent direction (grad(J)) times the  perturbation    !
!  amplitude "p" (controlled by inner loop).                           !
!                                                                      !
!=======================================================================
!
      USE mod_param
#  ifdef SOLVE3D
      USE mod_coupling
#  endif
      USE mod_grid
      USE mod_ocean
!
      integer, intent(in) :: ng, Linp, Lout, tile

#  include "tile.h"
!
#  ifdef PROFILE
      CALL wclock_on (ng, iNLM, 7)
#  endif
      CALL ini_perturb_tile (ng, Istr, Iend, Jstr, Jend,                &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       Linp, Lout,                                &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
     &                       GRID(ng) % umask,                          &
     &                       GRID(ng) % vmask,                          &
#  endif
#  ifdef SOLVE3D
#   if defined SEDIMENT && defined SED_MORPH
     &                       OCEAN(ng) %  bed,                          &
     &                       GRID(ng) % bed_thick0,                     &
#   endif
     &                       GRID(ng) % Hz,                             &
     &                       GRID(ng) % h,                              &
     &                       GRID(ng) % om_v,                           &
     &                       GRID(ng) % on_u,                           &
#   ifdef ICESHELF
     &                       GRID(ng) % zice,                           &
#   endif
     &                       GRID(ng) % z_r,                            &
     &                       GRID(ng) % z_w,                            &
     &                       COUPLING(ng) % Zt_avg1,                    &
     &                       OCEAN(ng) % ad_t,                          &
     &                       OCEAN(ng) % ad_u,                          &
     &                       OCEAN(ng) % ad_v,                          &
     &                       OCEAN(ng) % t,                             &
     &                       OCEAN(ng) % u,                             &
     &                       OCEAN(ng) % v,                             &
#  endif
     &                       OCEAN(ng) % ad_ubar,                       &
     &                       OCEAN(ng) % ad_vbar,                       &
     &                       OCEAN(ng) % ad_zeta,                       &
     &                       OCEAN(ng) % ubar,                          &
     &                       OCEAN(ng) % vbar,                          &
     &                       OCEAN(ng) % zeta)
#  ifdef PROFILE
      CALL wclock_off (ng, iNLM, 7)
#  endif
      RETURN
      END SUBROUTINE ini_perturb
!
!***********************************************************************
      SUBROUTINE ini_perturb_tile (ng, Istr, Iend, Jstr, Jend,          &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             Linp, Lout,                          &
#  ifdef MASKING
     &                             rmask, umask, vmask,                 &
#  endif
#  ifdef SOLVE3D
#   if defined SEDIMENT && defined SED_MORPH
     &                             bed, bed_thick0,                     &
#   endif
     &                             Hz, h, om_v, on_u,                   &
#   ifdef ICESHELF
     &                             zice,                                &
#   endif
     &                             z_r, z_w, Zt_avg1,                   &
     &                             ad_t, ad_u, ad_v,                    &
     &                             t, u, v,                             &
#  endif
     &                             ad_ubar, ad_vbar, ad_zeta,           &
     &                             ubar, vbar, zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_ncparam
      USE mod_iounits
      USE mod_scalars
!
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      USE exchange_2d_mod
#   ifdef SOLVE3D
      USE exchange_3d_mod
#   endif
#  endif
#  ifdef SOLVE3D
      USE set_depth_mod, ONLY : set_depth_tile
#  endif
      USE u2dbc_mod, ONLY : u2dbc_tile
      USE v2dbc_mod, ONLY : v2dbc_tile
      USE zetabc_mod, ONLY : zetabc_tile
#  ifdef SOLVE3D
      USE t3dbc_mod, ONLY : t3dbc_tile
      USE u3dbc_mod, ONLY : u3dbc_tile
      USE v3dbc_mod, ONLY : v3dbc_tile
#  endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Linp, Lout

#  ifdef ASSUMED_SHAPE
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#   endif
#   ifdef SOLVE3D
#    if defined SEDIMENT && defined SED_MORPH
      real(r8), intent(in) :: bed(LBi:,LBj:,:,:)
#    endif
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)
#    ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:,LBj:)
#    endif
      real(r8), intent(in) :: ad_t(LBi:,LBj:,:,:,:)
      real(r8), intent(in) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: ad_v(LBi:,LBj:,:,:)
#   endif
      real(r8), intent(in) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: ad_vbar(LBi:,LBj:,:)
      real(r8), intent(in) :: ad_zeta(LBi:,LBj:,:)
#   ifdef SOLVE3D
#    if defined SEDIMENT && defined SED_MORPH
      real(r8), intent(inout) :: bed_thick0(LBi:,LBj:)
#    endif
      real(r8), intent(inout) :: Hz(LBi:,LBj:,:)
      real(r8), intent(inout) :: Zt_avg1(LBi:,LBj:)
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: v(LBi:,LBj:,:,:)
#   endif
      real(r8), intent(inout) :: ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: vbar(LBi:,LBj:,:)
      real(r8), intent(inout) :: zeta(LBi:,LBj:,:)
#   ifdef SOLVE3D
      real(r8), intent(out) :: z_r(LBi:,LBj:,:)
      real(r8), intent(out) :: z_w(LBi:,LBj:,0:)
#   endif
#  else
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#   endif
#   ifdef SOLVE3D
#    if defined SEDIMENT && defined SED_MORPH
      real(r8), intent(in) :: bed(LBi:UBi,LBj:UBj,Nbed,MBEDP)
#    endif
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)
#    ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:,LBj:)
#    endif
      real(r8), intent(in) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(in) :: ad_u(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: ad_v(LBi:UBi,LBj:UBj,N(ng))
#   endif
      real(r8), intent(in) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: ad_vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: ad_zeta(LBi:UBi,LBj:UBj,3)
#   ifdef SOLVE3D
#    if defined SEDIMENT && defined SED_MORPH
      real(r8), intent(inout) :: bed_thick0(LBi:UBi,LBj:UBj)
#    endif
      real(r8), intent(inout) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: Zt_avg1(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: u(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: v(LBi:UBi,LBj:UBj,N(ng))
#   endif
      real(r8), intent(inout) :: ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: zeta(LBi:UBi,LBj:UBj,3)
#   ifdef SOLVE3D
      real(r8), intent(out) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(out) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
#   endif
#  endif
!
!  Local variable declarations.
!
      Logical :: FirstPass =.TRUE.

      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j
#  ifdef SOLVE3D
      integer :: itrc, k
#  endif
      integer, dimension(NstateVar(ng)+1) :: StateVar

      real(r8) :: p, scale

#  include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Add a perturbation to nonlinear state variables: steepest descent
!  direction times the perturbation amplitude "p".
!-----------------------------------------------------------------------
!
!  Set state variable to perturb according to outer loop index.
!
#  ifdef SOLVE3D
      StateVar(1)=0
      StateVar(2)=isFsur
      StateVar(3)=isUbar
      StateVar(4)=isVbar
      StateVar(5)=isUvel
      StateVar(6)=isVvel
      DO i=1,NT(ng)
        StateVar(6+i)=isTvar(i)
      END DO
#  else
      StateVar(1)=0
      StateVar(2)=isFsur
      StateVar(3)=isUbar
      StateVar(4)=isVbar
#  endif
!
!  Set perturbation amplitude as a function of the inner loop.
!
      p=10.0_r8**REAL(-inner,r8)
      scale=1.0_r8/SQRT(adDotProduct)
      IF (SOUTH_WEST_CORNER) THEN
        IF (Master) WRITE (stdout,10)
      END IF
!
!  Free-surface.
!
      IF ((StateVar(outer).eq.0).or.(StateVar(outer).eq.isFsur)) THEN
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            zeta(i,j,Lout)=zeta(i,j,Lout)+p*ad_zeta(i,j,Linp)*scale
#  ifdef MASKING
            zeta(i,j,Lout)=zeta(i,j,Lout)*rmask(i,j)
#  endif
          END DO
        END DO
        IF (SOUTH_WEST_TEST) THEN
          IF (Master) THEN
            WRITE (stdout,20) outer, inner,                             &
     &                        TRIM(Vname(1,idSvar(isFsur)))
          END IF
        END IF
      END IF
!
!  2D u-momentum component.
!
      IF ((StateVar(outer).eq.0).or.(StateVar(outer).eq.isUbar)) THEN
        DO j=JstrR,JendR
          DO i=Istr,IendR
            ubar(i,j,Lout)=ubar(i,j,Lout)+p*ad_ubar(i,j,Linp)*scale
#  ifdef MASKING
            ubar(i,j,Lout)=ubar(i,j,Lout)*umask(i,j)
#  endif
          END DO
        END DO
        IF (SOUTH_WEST_TEST) THEN
          IF (Master) THEN
            WRITE (stdout,20) outer, inner,                             &
     &                        TRIM(Vname(1,idSvar(isUbar)))
          END IF
        END IF
      END IF
!
!  2D v-momentum component.
!
      IF ((StateVar(outer).eq.0).or.(StateVar(outer).eq.isVbar)) THEN
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            vbar(i,j,Lout)=vbar(i,j,Lout)+p*ad_vbar(i,j,Linp)*scale
#  ifdef MASKING
            vbar(i,j,Lout)=vbar(i,j,Lout)*vmask(i,j)
#  endif
          END DO
        END DO
        IF (SOUTH_WEST_TEST) THEN
          IF (Master) THEN
            WRITE (stdout,20) outer, inner,                             &
     &                        TRIM(Vname(1,idSvar(isVbar)))
          END IF
        END IF
      END IF
#  ifdef SOLVE3D
!
!  3D u-momentum component.
!
      IF ((StateVar(outer).eq.0).or.(StateVar(outer).eq.isUvel)) THEN
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=Istr,IendR
              u(i,j,k,Lout)=u(i,j,k,Lout)+p*ad_u(i,j,k,Linp)*scale
#   ifdef MASKING
              u(i,j,k,Lout)=u(i,j,k,Lout)*umask(i,j)
#   endif
            END DO
          END DO
        END DO
        IF (SOUTH_WEST_TEST) THEN
          IF (Master) THEN
            WRITE (stdout,20) outer, inner,                             &
     &                        TRIM(Vname(1,idSvar(isUvel)))
          END IF
        END IF
      END IF
!
!  3D v-momentum component.
!
      IF ((StateVar(outer).eq.0).or.(StateVar(outer).eq.isVvel)) THEN
        DO k=1,N(ng)
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              v(i,j,k,Lout)=v(i,j,k,Lout)+p*ad_v(i,j,k,Linp)*scale
#   ifdef MASKING
              v(i,j,k,Lout)=v(i,j,k,Lout)*vmask(i,j)
#   endif
            END DO
          END DO
        END DO
        IF (SOUTH_WEST_TEST) THEN
          IF (Master) THEN
            WRITE (stdout,20) outer, inner,                             &
     &                        TRIM(Vname(1,idSvar(isVvel)))
          END IF
        END IF
      END IF
!
!  Tracers.
!
      DO itrc=1,NT(ng)
        IF ((StateVar(outer).eq.0).or.                                  &
     &      (StateVar(outer).eq.isTvar(itrc))) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                t(i,j,k,Lout,itrc)=t(i,j,k,Lout,itrc)+                  &
     &                             p*ad_t(i,j,k,Linp,itrc)*scale
#   ifdef MASKING
                t(i,j,k,Lout,itrc)=t(i,j,k,Lout,itrc)*rmask(i,j)
#   endif
              END DO
            END DO
          END DO
          IF (SOUTH_WEST_TEST) THEN
            IF (Master) THEN
              WRITE (stdout,20) outer, inner,                           &
     &                          TRIM(Vname(1,idSvar(isTvar(itrc))))
            END IF
          END IF
        END IF
      END DO
#  endif
      IF (Master) WRITE (stdout,'(/)')
!
!-----------------------------------------------------------------------
!  Apply lateral boundary conditions to 2D fields.
!-----------------------------------------------------------------------
!
      CALL zetabc_tile (ng, Istr, Iend, Jstr, Jend,                     &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Lout, Lout, Lout,                               &
     &                  zeta)
#  ifndef SOLVE3D
      CALL u2dbc_tile (ng, Istr, Iend, Jstr, Jend,                      &
     &                 LBi, UBi, LBj, UBj,                              &
     &                 Lout, Lout, Lout,                                &
     &                 ubar, vbar, zeta)
      CALL v2dbc_tile (ng, Istr, Iend, Jstr, Jend,                      &
     &                 LBi, UBi, LBj, UBj,                              &
     &                 Lout, Lout, Lout,                                &
     &                 ubar, vbar, zeta)
#  endif
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      CALL exchange_r2d_tile (ng, model, Istr, Iend, Jstr, Jend,        &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        zeta(:,:,Lout))
#   ifndef SOLVE3D
      CALL exchange_u2d_tile (ng, model, Istr, Iend, Jstr, Jend,        &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        ubar(:,:,Lout))
      CALL exchange_v2d_tile (ng, model, Istr, Iend, Jstr, Jend,        &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        vbar(:,:,Lout))
#   endif
#  endif
#  ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Compute new depths and thicknesses.
!-----------------------------------------------------------------------
!
      CALL set_depth_tile (ng, Istr, Iend, Jstr, Jend,                  &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     Lout, FirstPass,                             &
     &                     zeta, h,                                     &
#   ifdef ICESHELF
     &                     zice,                                        &
#   endif
#   if defined SEDIMENT && defined SED_MORPH
     &                     bed, bed_thick0,                             &
#   endif
     &                     Zt_avg1,                                     &
     &                     Hz, z_r, z_w)
!
!-----------------------------------------------------------------------
!  Apply lateral boundary conditions to perturbed 3D fields.
!-----------------------------------------------------------------------
!
      CALL u3dbc_tile (ng, Istr, Iend, Jstr, Jend,                      &
     &                 LBi, UBi, LBj, UBj, N(ng),                       &
     &                 Lout, Lout, u)
      CALL v3dbc_tile (ng, Istr, Iend, Jstr, Jend,                      &
     &                 LBi, UBi, LBj, UBj, N(ng),                       &
     &                 Lout, Lout, v)
      DO itrc=1,NT(ng)
        CALL t3dbc_tile (ng, Istr, Iend, Jstr, Jend, itrc,              &
     &                   LBi, UBi, LBj, UBj, N(ng), NT(ng),             &
     &                   Lout, Lout, t)
      END DO
#   if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      CALL exchange_u3d_tile (ng, model, Istr, Iend, Jstr, Jend,        &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        NghostPoints,                             &
     &                        u(:,:,:,Lout))
      CALL exchange_v3d_tile (ng, model, Istr, Iend, Jstr, Jend,        &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        NghostPoints,                             &
     &                        v(:,:,:,Lout))
      DO itrc=1,NT(ng)
        CALL exchange_r3d_tile (ng, Istr, Iend, Jstr, Jend,             &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          t(:,:,:,Lout,itrc))
      END DO
#   endif
#  endif
!
 10   FORMAT (/,' Perturbing Nonlinear model Initial Conditions:',/)
 20   FORMAT (' (Outer,Inner) = ','(',i4.4,',',i4.4,')',3x,             &
     &        'Perturbing State Variable: ',a)
      RETURN
      END SUBROUTINE ini_perturb_tile
# endif

# if defined PROPAGATOR && defined ADJOINT

      SUBROUTINE ad_ini_perturb (ng, tile, Kinp, Kout, Ninp, Nout)
!
!=======================================================================
!                                                                      !
!  Initialize adjoint state variables with tangent linear state scaled !
!  by the energy norm, as required by the Generalized stability Theory !
!  propagators.                                                        !
!                                                                      !
!=======================================================================
!
      USE mod_param
#  ifdef SOLVE3D
      USE mod_coupling
#  endif
      USE mod_grid
      USE mod_ocean
!
      integer, intent(in) :: ng, tile, Kinp, Kout, Ninp, Nout

#  include "tile.h"
!
#  ifdef PROFILE
      CALL wclock_on (ng, iADM, 7)
#  endif
      CALL ad_ini_perturb_tile (ng, Istr, Iend, Jstr, Jend,             &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Kinp, Kout, Ninp, Nout,                 &
#  ifdef MASKING
     &                          GRID(ng) % rmask,                       &
     &                          GRID(ng) % umask,                       &
     &                          GRID(ng) % vmask,                       &
#  endif
#  ifdef SOLVE3D
     &                          GRID(ng) % Hz,                          &
#  else
     &                          OCEAN(ng) % tl_ubar,                    &
     &                          OCEAN(ng) % tl_vbar,                    &
#  endif
     &                          OCEAN(ng) % tl_zeta,                    &
#  ifdef SOLVE3D
     &                          OCEAN(ng) % tl_t,                       &
     &                          OCEAN(ng) % tl_u,                       &
     &                          OCEAN(ng) % tl_v,                       &
     &                          OCEAN(ng) % ad_t,                       &
     &                          OCEAN(ng) % ad_u,                       &
     &                          OCEAN(ng) % ad_v,                       &
#  else
     &                          OCEAN(ng) % ad_ubar,                    &
     &                          OCEAN(ng) % ad_vbar,                    &
#  endif
     &                          OCEAN(ng) % ad_zeta)
#  ifdef PROFILE
      CALL wclock_off (ng, iADM, 7)
#  endif
      RETURN
      END SUBROUTINE ad_ini_perturb
!
!***********************************************************************
      SUBROUTINE ad_ini_perturb_tile (ng, Istr, Iend, Jstr, Jend,       &
     &                                LBi, UBi, LBj, UBj,               &
     &                                Kinp, Kout, Ninp, Nout,           &
#  ifdef MASKING
     &                                rmask, umask, vmask,              &
#  endif
#  ifdef SOLVE3D
     &                                Hz,                               &
#  else
     &                                tl_ubar, tl_vbar,                 &
#  endif
     &                                tl_zeta,                          &
#  ifdef SOLVE3D
     &                                tl_t, tl_u, tl_v,                 &
     &                                ad_t, ad_u, ad_v,                 &
#  else
     &                                ad_ubar, ad_vbar,                 &
#  endif
     &                                ad_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Kinp, Kout, Ninp, Nout

#  ifdef ASSUMED_SHAPE
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#   endif
#   ifdef SOLVE3D
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
#   else
      real(r8), intent(in) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: tl_vbar(LBi:,LBj:,:)
#   endif
      real(r8), intent(in) :: tl_zeta(LBi:,LBj:,:)
#   ifdef SOLVE3D
      real(r8), intent(in) :: tl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(in) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: tl_v(LBi:,LBj:,:,:)
#   else
      real(r8), intent(inout) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_vbar(LBi:,LBj:,:)
#   endif
      real(r8), intent(inout) :: ad_zeta(LBi:,LBj:,:)
#   ifdef SOLVE3D
      real(r8), intent(inout) :: ad_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_v(LBi:,LBj:,:,:)
#   endif
#  else
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#   endif
#   ifdef SOLVE3D
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
#   else
      real(r8), intent(in) :: tl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: tl_vbar(LBi:UBi,LBj:UBj,3)
#   endif
      real(r8), intent(in) :: tl_zeta(LBi:UBi,LBj:UBj,3)
#   ifdef SOLVE3D
      real(r8), intent(in) :: tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(in) :: tl_u(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: tl_v(LBi:UBi,LBj:UBj,N(ng))
#   else
      real(r8), intent(inout) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_vbar(LBi:UBi,LBj:UBj,3)
#   endif
      real(r8), intent(inout) :: ad_zeta(LBi:UBi,LBj:UBj,3)
#   ifdef SOLVE3D
      real(r8), intent(inout) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: ad_u(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: ad_v(LBi:UBi,LBj:UBj,N(ng))
#   endif
#  endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j
#  ifdef SOLVE3D
      integer :: itrc, k
#  endif
      real(r8) :: cff, scale
!
#  include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Initialize adjoint state variables with tangent linear state scaled
!  by the energy norm.
!-----------------------------------------------------------------------

#  ifdef FULL_GRID
#   define IR_RANGE IstrR,IendR
#   define IU_RANGE Istr,IendR
#   define JR_RANGE JstrR,JendR
#   define JV_RANGE Jstr,JendR
#  else
#   define IR_RANGE Istr,Iend
#   define IU_RANGE IstrU,Iend
#   define JR_RANGE Jstr,Jend
#   define JV_RANGE JstrV,Jend
#  endif
!
!  Free-surface.
!
      scale=0.5_r8*g*rho0
      DO j=JR_RANGE
        DO i=IR_RANGE
          ad_zeta(i,j,Kout)=scale*tl_zeta(i,j,Kinp)
#  ifdef MASKING
          ad_zeta(i,j,Kout)=ad_zeta(i,j,Kout)*rmask(i,j)
#  endif
        END DO
      END DO

#  ifndef SOLVE3D
!
!  2D u-momentum component.
!
      cff=0.25_r8*rho0
      DO j=JR_RANGE
        DO i=IU_RANGE
          scale=cff*(h(i-1,j)+h(i,j))
          ad_ubar(i,j,Kout)=scale*tl_ubar(i,j,Kinp)
#   ifdef MASKING
          ad_ubar(i,j,Kout)=scale*ad_ubar(i,j,Kout)*umask(i,j)
#   endif
        END DO
      END DO
!
!  2D v-momentum component.
!
      cff=0.25_r8*rho0
      DO j=JV_RANGE
        DO i=IR_RANGE
          scale=cff*(h(i,j-1)+h(i,j))
          ad_vbar(i,j,Kout)=scale*tl_vbar(i,j,Kinp)
#   ifdef MASKING
          ad_vbar(i,j,Kout)=ad_vbar(i,j,Kout)*vmask(i,j)
#   endif
        END DO
      END DO
#  else
!
!  3D u-momentum component.
!
      cff=0.25_r8*rho0
      DO k=1,N(ng)
        DO j=JR_RANGE
          DO i=IU_RANGE
            scale=cff*(Hz(i-1,j,j)+Hz(i,j,k))
            ad_u(i,j,k,Nout)=scale*tl_u(i,j,k,Ninp)
#   ifdef MASKING
            ad_u(i,j,k,Nout)=ad_u(i,j,k,Nout)*umask(i,j)
#   endif
          END DO
        END DO
      END DO
!
!  3D v-momentum component.
!
      cff=0.25_r8*rho0
      DO k=1,N(ng)
        DO j=JV_RANGE
          DO i=IR_RANGE
            scale=cff*(Hz(i,j-1,j)+Hz(i,j,k))
            ad_v(i,j,k,Nout)=scale*tl_v(i,j,k,Ninp)
#   ifdef MASKING
            ad_v(i,j,k,Nout)=ad_v(i,j,k,Nout)*vmask(i,j)
#   endif
          END DO
        END DO
      END DO
!
!  Tracers. For now, use salinity scale for passive tracers.
!
      DO itrc=1,NT(ng)
        IF (itrc.eq.itemp) THEN
          cff=0.5_r8*rho0*Tcoef(ng)*Tcoef(ng)*g*g/bvf_bak
        ELSE IF (itrc.eq.isalt) THEN
          cff=0.5_r8*rho0*Scoef(ng)*Scoef(ng)*g*g/bvf_bak
        ELSE
          cff=0.5_r8*rho0*Scoef(ng)*Scoef(ng)*g*g/bvf_bak
        END IF
        DO k=1,N(ng)
          DO j=JR_RANGE
            DO i=IR_RANGE
              scale=cff*Hz(i,j,k)
              ad_t(i,j,k,Nout,itrc)=scale*tl_t(i,j,k,Ninp,itrc)
#   ifdef MASKING
              ad_t(i,j,k,Nout,itrc)=ad_t(i,j,k,Nout,itrc)*rmask(i,j)
#   endif
            END DO
          END DO
        END DO
      END DO
#  endif

#  undef IR_RANGE
#  undef IU_RANGE
#  undef JR_RANGE
#  undef JV_RANGE

      RETURN
      END SUBROUTINE ad_ini_perturb_tile
# endif

# if defined TLM_CHECK

      SUBROUTINE tl_ini_perturb (ng, tile, Linp, Lout)
!
!=======================================================================
!                                                                      !
!  Initialize tangent linear state variable according to the outer     !
!  and inner loop iterations.  The initial field is  a function of     !
!  the steepest descent direction (grad(J)) times the perturbation     !
!  amplitude "p" (controlled by inner loop).                           !
!                                                                      !
!=======================================================================
!
      USE mod_param
#  ifdef SOLVE3D
      USE mod_coupling
#  endif
      USE mod_grid
      USE mod_ocean
!
      integer, intent(in) :: ng, tile, Linp, Lout

#  include "tile.h"
!
#  ifdef PROFILE
      CALL wclock_on (ng, iTLM, 7)
#  endif
      CALL tl_ini_perturb_tile (ng, Istr, Iend, Jstr, Jend,             &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Linp, Lout,                             &
#  ifdef SOLVE3D
#   if defined SEDIMENT && defined SED_MORPH
     &                          OCEAN(ng) %  tl_bed,                    &
     &                          GRID(ng) % bed_thick0,                  &
#   endif
     &                          GRID(ng) % tl_Hz,                       &
     &                          GRID(ng) % h,                           &
     &                          GRID(ng) % tl_h,                        &
     &                          GRID(ng) % om_v,                        &
     &                          GRID(ng) % on_u,                        &
#   ifdef ICESHELF
     &                          GRID(ng) % zice,                        &
#   endif
     &                          GRID(ng) % tl_z_r,                      &
     &                          GRID(ng) % tl_z_w,                      &
     &                          COUPLING(ng) % Zt_avg1,                 &
     &                          COUPLING(ng) % tl_Zt_avg1,              &
#  endif
     &                          OCEAN(ng) % ubar,                       &
     &                          OCEAN(ng) % vbar,                       &
     &                          OCEAN(ng) % zeta,                       &
#  ifdef SOLVE3D
     &                          OCEAN(ng) % ad_t,                       &
     &                          OCEAN(ng) % ad_u,                       &
     &                          OCEAN(ng) % ad_v,                       &
     &                          OCEAN(ng) % tl_t,                       &
     &                          OCEAN(ng) % tl_u,                       &
     &                          OCEAN(ng) % tl_v,                       &
#  endif
     &                          OCEAN(ng) % ad_ubar,                    &
     &                          OCEAN(ng) % ad_vbar,                    &
     &                          OCEAN(ng) % ad_zeta,                    &
     &                          OCEAN(ng) % tl_ubar,                    &
     &                          OCEAN(ng) % tl_vbar,                    &
     &                          OCEAN(ng) % tl_zeta)
#  ifdef PROFILE
      CALL wclock_off (ng, iTLM, 7)
#  endif
      RETURN
      END SUBROUTINE tl_ini_perturb
!
!***********************************************************************
      SUBROUTINE tl_ini_perturb_tile (ng, Istr, Iend, Jstr, Jend,       &
     &                                LBi, UBi, LBj, UBj,               &
     &                                Linp, Lout,                       &
#  ifdef SOLVE3D
#   if defined SEDIMENT && defined SED_MORPH
     &                                tl_bed, bed_thick0,               &
#   endif
     &                                tl_Hz, h, tl_h,                   &
     &                                om_v, on_u,                       &
#   ifdef ICESHELF
     &                                zice,                             &
#   endif
     &                                tl_z_r, tl_z_w,                   &
     &                                Zt_avg1, tl_Zt_avg1,              &
#  endif
     &                                ubar, vbar, zeta,                 &
#  ifdef SOLVE3D
     &                                ad_t, ad_u, ad_v,                 &
     &                                tl_t, tl_u, tl_v,                 &
#  endif
     &                                ad_ubar, ad_vbar, ad_zeta,        &
     &                                tl_ubar, tl_vbar, tl_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_ncparam
      USE mod_iounits
      USE mod_scalars
!
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      USE exchange_2d_mod
#   ifdef SOLVE3D
      USE exchange_3d_mod
#   endif
#  endif
#  ifdef SOLVE3D
      USE tl_set_depth_mod, ONLY : tl_set_depth_tile
#  endif
      USE tl_u2dbc_mod, ONLY : tl_u2dbc_tile
      USE tl_v2dbc_mod, ONLY : tl_v2dbc_tile
      USE tl_zetabc_mod, ONLY : tl_zetabc_tile
#  ifdef SOLVE3D
      USE tl_t3dbc_mod, ONLY : tl_t3dbc_tile
      USE tl_u3dbc_mod, ONLY : tl_u3dbc_tile
      USE tl_v3dbc_mod, ONLY : tl_v3dbc_tile
#  endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Linp, Lout

#  ifdef ASSUMED_SHAPE
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#   endif
#   ifdef SOLVE3D
#    if defined SEDIMENT && defined SED_MORPH
      real(r8), intent(in) :: tl_bed(LBi:,LBj:,:,:)
#    endif
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)
#    ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:,LBj:)
#    endif
      real(r8), intent(in) :: Zt_avg1(LBi:,LBj:)
#   endif
      real(r8), intent(in) :: ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: vbar(LBi:,LBj:,:)
      real(r8), intent(in) :: zeta(LBi:,LBj:,:)
#   ifdef SOLVE3D
      real(r8), intent(in) :: ad_t(LBi:,LBj:,:,:,:)
      real(r8), intent(in) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: ad_v(LBi:,LBj:,:,:)
#   endif
      real(r8), intent(in) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: ad_vbar(LBi:,LBj:,:)
      real(r8), intent(in) :: ad_zeta(LBi:,LBj:,:)
#   ifdef SOLVE3D
#    if defined SEDIMENT && defined SED_MORPH
      real(r8), intent(inout) :: bed_thick0(LBi:,LBj:)
#    endif
      real(r8), intent(inout) :: tl_Zt_avg1(LBi:,LBj:)
#   endif
#   ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_v(LBi:,LBj:,:,:)
#   endif
      real(r8), intent(inout) :: tl_h(LBi:,LBj:)
      real(r8), intent(inout) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_vbar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_zeta(LBi:,LBj:,:)
#   ifdef SOLVE3D
      real(r8), intent(out) :: tl_Hz(LBi:,LBj:,:)
      real(r8), intent(out) :: tl_z_r(LBi:,LBj:,:)
      real(r8), intent(out) :: tl_z_w(LBi:,LBj:,0:)
#   endif
#  else
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#   endif
#   ifdef SOLVE3D
#    if defined SEDIMENT && defined SED_MORPH
      real(r8), intent(in) :: tl_bed(LBi:UBi,LBj:UBj,Nbed,MBEDP)
#    endif
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)
#    ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:,LBj:)
#    endif
       real(r8), intent(in) :: Zt_avg1(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(in) :: ubar(LBi;UBi,LBj:UBj,3)
      real(r8), intent(in) :: vbar(LBi;UBi,LBj:UBj,3)
      real(r8), intent(in) :: zeta(LBi;UBi,LBj:UBj,3)
#   ifdef SOLVE3D
      real(r8), intent(in) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(in) :: ad_u(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: ad_v(LBi:UBi,LBj:UBj,N(ng))
#   endif
      real(r8), intent(in) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: ad_vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: ad_zeta(LBi:UBi,LBj:UBj,3)
#   ifdef SOLVE3D
#    if defined SEDIMENT && defined SED_MORPH
      real(r8), intent(inout) :: bed_thick0(LBi:UBi,LBj:UBj)
#    endif
      real(r8), intent(inout) :: tl_Zt_avg1(LBi:UBi,LBj:UBj)
#   endif
#   ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: tl_u(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: tl_v(LBi:UBi,LBj:UBj,N(ng))
#   endif
      real(r8), intent(inout) :: tl_h(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: tl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: tl_vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: tl_zeta(LBi:UBi,LBj:UBj,3)
#   ifdef SOLVE3D
      real(r8), intent(out) :: tl_Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(out) :: tl_z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(out) :: tl_z_w(LBi:UBi,LBj:UBj,0:N(ng))
#   endif
#  endif
!
!  Local variable declarations.
!
      Logical :: FirstPass =.TRUE.

      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j
#  ifdef SOLVE3D
      integer :: itrc, k
#  endif
      integer, dimension(NstateVar(ng)+1) :: StateVar

      real(r8) :: p, scale
!
#  include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Initialize tangent linear with the steepest descent direction times
!  the perturbation amplitude "p".
!-----------------------------------------------------------------------
!
!  Set state variable to perturb according to outer loop index.
!
#  ifdef SOLVE3D
      StateVar(1)=0
      StateVar(2)=isFsur
      StateVar(3)=isUbar
      StateVar(4)=isVbar
      StateVar(5)=isUvel
      StateVar(6)=isVvel
      DO i=1,NT(ng)
        StateVar(6+i)=isTvar(i)
      END DO
#  else
      StateVar(1)=0
      StateVar(2)=isFsur
      StateVar(3)=isUbar
      StateVar(4)=isVbar
#  endif
!
!  Set perturbation amplitude as a function of the inner loop.
!
      p=10.0_r8**REAL(-inner,r8)
      scale=1.0_r8/SQRT(adDotProduct)
      IF (SOUTH_WEST_TEST) THEN
        IF (Master) WRITE (stdout,10)
      END IF
!
!  Free-surface.
!
      IF ((StateVar(outer).eq.0).or.(StateVar(outer).eq.isFsur)) THEN
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            tl_zeta(i,j,Lout)=p*ad_zeta(i,j,Linp)*scale
#  ifdef MASKING
            tl_zeta(i,j,Lout)=tl_zeta(i,j,Lout)*rmask(i,j)
#  endif
          END DO
        END DO
        IF (SOUTH_WEST_TEST) THEN
          IF (Master) THEN
            WRITE (stdout,20) outer, inner,                             &
     &                        TRIM(Vname(1,idSvar(isFsur)))
          END IF
        END IF
      END IF
!
!  2D u-momentum component.
!
      IF ((StateVar(outer).eq.0).or.(StateVar(outer).eq.isUbar)) THEN
        DO j=JstrR,JendR
          DO i=Istr,IendR
            tl_ubar(i,j,Lout)=p*ad_ubar(i,j,Linp)*scale
#  ifdef MASKING
            tl_ubar(i,j,Lout)=tl_ubar(i,j,Lout)*umask(i,j)
#  endif
          END DO
        END DO
        IF (SOUTH_WEST_TEST) THEN
          IF (Master) THEN
            WRITE (stdout,20) outer, inner,                             &
     &                        TRIM(Vname(1,idSvar(isUbar)))
          END IF
        END IF
      END IF
!
!  2D v-momentum component.
!
      IF ((StateVar(outer).eq.0).or.(StateVar(outer).eq.isVbar)) THEN
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            tl_vbar(i,j,Lout)=p*ad_vbar(i,j,Linp)*scale
#  ifdef MASKING
            tl_vbar(i,j,Lout)=tl_vbar(i,j,Lout)*vmask(i,j)
#  endif
          END DO
        END DO
        IF (SOUTH_WEST_TEST) THEN
          IF (Master) THEN
            WRITE (stdout,20) outer, inner,                             &
     &                        TRIM(Vname(1,idSvar(isVbar)))
          END IF
        END IF
      END IF
#  ifdef SOLVE3D
!
!  3D u-momentum component.
!
      IF ((StateVar(outer).eq.0).or.(StateVar(outer).eq.isUvel)) THEN
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=Istr,IendR
              tl_u(i,j,k,Lout)=p*ad_u(i,j,k,Linp)*scale
#   ifdef MASKING
              tl_u(i,j,k,Lout)=tl_u(i,j,k,Lout)*umask(i,j)
#   endif
            END DO
          END DO
        END DO
        IF (SOUTH_WEST_TEST) THEN
          IF (Master) THEN
            WRITE (stdout,20) outer, inner,                             &
     &                        TRIM(Vname(1,idSvar(isUvel)))
          END IF
        END IF
      END IF
!
!  3D v-momentum component.
!
      IF ((StateVar(outer).eq.0).or.(StateVar(outer).eq.isVvel)) THEN
        DO k=1,N(ng)
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              tl_v(i,j,k,Lout)=p*ad_v(i,j,k,Linp)*scale
#   ifdef MASKING
              tl_v(i,j,k,Lout)=tl_v(i,j,k,Lout)*vmask(i,j)
#   endif
            END DO
          END DO
        END DO
        IF (SOUTH_WEST_TEST) THEN
          IF (Master) THEN
            WRITE (stdout,20) outer, inner,                             &
     &                        TRIM(Vname(1,idSvar(isVvel)))
          END IF
        END IF
      END IF
!
!  Tracers.
!
      DO itrc=1,NT(ng)
        IF ((StateVar(outer).eq.0).or.                                  &
     &      (StateVar(outer).eq.isTvar(itrc))) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                tl_t(i,j,k,Lout,itrc)=p*ad_t(i,j,k,Linp,itrc)*scale
#   ifdef MASKING
                tl_t(i,j,k,Lout,itrc)=tl_t(i,j,k,Lout,itrc)*rmask(i,j)
#   endif
              END DO
            END DO
          END DO
          IF (SOUTH_WEST_TEST) THEN
            IF (Master) THEN
              WRITE (stdout,20) outer, inner,                           &
     &                          TRIM(Vname(1,idSvar(isTvar(itrc))))
            END IF
          END IF
        END IF
      END DO
#  endif
      IF (Master) WRITE (stdout,'(/)')
!
!-----------------------------------------------------------------------
!  Apply lateral boundary conditions to 2D fields.
!-----------------------------------------------------------------------
!
      CALL tl_zetabc_tile (ng, Istr, Iend, Jstr, Jend,                  &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     Lout, Lout, Lout,                            &
     &                     zeta, tl_zeta)
#  ifndef SOLVE3D
      CALL tl_u2dbc_tile (ng, Istr, Iend, Jstr, Jend,                   &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    Lout, Lout, Lout,                             &
     &                    ubar, vbar, zeta,                             &
     &                    tl_ubar, tl_vbar, tl_zeta)
      CALL tl_v2dbc_tile (ng, Istr, Iend, Jstr, Jend,                   &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    Lout, Lout, Lout,                             &
     &                    ubar, vbar, zeta,                             &
     &                    tl_ubar, tl_vbar, tl_zeta)
#  endif
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      CALL exchange_r2d_tile (ng, model, Istr, Iend, Jstr, Jend,        &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        tl_zeta(:,:,Lout))
#   ifndef SOLVE3D
      CALL exchange_u2d_tile (ng, model, Istr, Iend, Jstr, Jend,        &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        tl_ubar(:,:,Lout))
      CALL exchange_v2d_tile (ng, model, Istr, Iend, Jstr, Jend,        &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        tl_vbar(:,:,Lout))
#   endif
#  endif
#  ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Compute new depths and thicknesses.
!-----------------------------------------------------------------------
!
      CALL tl_set_depth_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Lout, FirstPass,                          &
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
!
!-----------------------------------------------------------------------
!  Apply lateral boundary conditions to perturbed 3D fields.
!-----------------------------------------------------------------------
!
      CALL tl_u3dbc_tile (ng, Istr, Iend, Jstr, Jend,                   &
     &                    LBi, UBi, LBj, UBj, N(ng),                    &
     &                    Lout, Lout, tl_u)
      CALL tl_v3dbc_tile (ng, Istr, Iend, Jstr, Jend,                   &
     &                    LBi, UBi, LBj, UBj, N(ng),                    &
     &                    Lout, Lout, tl_v)
      DO itrc=1,NT(ng)
        CALL tl_t3dbc_tile (ng, Istr, Iend, Jstr, Jend, itrc,           &
     &                      LBi, UBi, LBj, UBj, N(ng), NT(ng),          &
     &                      Lout, Lout, tl_t)
      END DO
#   if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      CALL exchange_u3d_tile (ng, model, Istr, Iend, Jstr, Jend,        &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        NghostPoints,                             &
     &                        tl_u(:,:,:,Lout))
      CALL exchange_v3d_tile (ng, model, Istr, Iend, Jstr, Jend,        &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        NghostPoints,                             &
     &                        tl_v(:,:,:,Lout))
      DO itrc=1,NT(ng)
        CALL exchange_r3d_tile (ng, Istr, Iend, Jstr, Jend,             &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          tl_t(:,:,:,Lout,itrc))
      END DO
#   endif
#  endif

!
 10   FORMAT (/,' Perturbing Tangent Linear Initial Conditions:',/)
 20   FORMAT (' (Outer,Inner) = ','(',i4.4,',',i4.4,')',3x,             &
     &        'Perturbing State Variable: ',a)
      RETURN
      END SUBROUTINE tl_ini_perturb_tile

# elif defined GRADIENT_CHECK

      SUBROUTINE tl_ini_perturb (ng, tile, Linp, Lout)
!
!=======================================================================
!                                                                      !
!  Initialize tangent linear state variable according to the outer     !
!  and inner loop iterations.  The initial field is  a function of     !
!  the steepest descent direction (grad(J)) times the perturbation     !
!  amplitude "p" (controlled by inner loop).                           !
!                                                                      !
!=======================================================================
!
      USE mod_param
#  ifdef SOLVE3D
      USE mod_coupling
#  endif
      USE mod_grid
      USE mod_ocean
!
      integer, intent(in) :: ng, tile, Linp, Lout

#  include "tile.h"
!
#  ifdef PROFILE
      CALL wclock_on (ng, iTLM, 7)
#  endif
      CALL tl_ini_perturb_tile (ng, Istr, Iend, Jstr, Jend,             &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Linp, Lout,                             &
#  ifdef MASKING
     &                          GRID(ng) % rmask,                       &
     &                          GRID(ng) % umask,                       &
     &                          GRID(ng) % vmask,                       &
#  endif
#  ifdef SOLVE3D
#   if defined SEDIMENT && defined SED_MORPH
     &                          OCEAN(ng) %  bed,                       &
     &                          OCEAN(ng) %  tl_bed,                    &
     &                          GRID(ng) % bed_thick0,                  &
#   endif
     &                          GRID(ng) % Hz,                          &
     &                          GRID(ng) % tl_Hz,                       &
     &                          GRID(ng) % h,                           &
     &                          GRID(ng) % tl_h,                        &
     &                          GRID(ng) % om_v,                        &
     &                          GRID(ng) % on_u,                        &
#   ifdef ICESHELF
     &                          GRID(ng) % zice,                        &
#   endif
     &                          GRID(ng) % z_r,                         &
     &                          GRID(ng) % tl_z_r,                      &
     &                          GRID(ng) % z_w,                         &
     &                          GRID(ng) % tl_z_w,                      &
     &                          COUPLING(ng) % Zt_avg1,                 &
     &                          COUPLING(ng) % tl_Zt_avg1,              &
#  endif
     &                          OCEAN(ng) % ubar,                       &
     &                          OCEAN(ng) % vbar,                       &
     &                          OCEAN(ng) % zeta,                       &
#  ifdef SOLVE3D
     &                          OCEAN(ng) % u,                          &
     &                          OCEAN(ng) % v,                          &
     &                          OCEAN(ng) % ad_t,                       &
     &                          OCEAN(ng) % ad_u,                       &
     &                          OCEAN(ng) % ad_v,                       &
     &                          OCEAN(ng) % tl_t,                       &
     &                          OCEAN(ng) % tl_u,                       &
     &                          OCEAN(ng) % tl_v,                       &
#  endif
     &                          OCEAN(ng) % ad_ubar,                    &
     &                          OCEAN(ng) % ad_vbar,                    &
     &                          OCEAN(ng) % ad_zeta,                    &
     &                          OCEAN(ng) % tl_ubar,                    &
     &                          OCEAN(ng) % tl_vbar,                    &
     &                          OCEAN(ng) % tl_zeta)
#  ifdef PROFILE
      CALL wclock_off (ng, iTLM, 7)
#  endif
      RETURN
      END SUBROUTINE tl_ini_perturb
!
!***********************************************************************
      SUBROUTINE tl_ini_perturb_tile (ng, Istr, Iend, Jstr, Jend,       &
     &                                LBi, UBi, LBj, UBj,               &
     &                                Linp, Lout,                       &
#  ifdef MASKING
     &                                rmask, umask, vmask,              &
#  endif
#  ifdef SOLVE3D
#   if defined SEDIMENT && defined SED_MORPH
     &                                bed, tl_bed,                      &
     &                                bed_thick0,                       &
#   endif
     &                                Hz, tl_Hz,                        &
     &                                h, tl_h,                          &
     &                                om_v, on_u,                       &
#   ifdef ICESHELF
     &                                zice,                             &
#   endif
     &                                z_r, tl_z_r,                      &
     &                                z_w, tl_z_w,                      &
     &                                Zt_avg1, tl_Zt_avg1,              &
#  endif
     &                                ubar, vbar, zeta,                 &
#  ifdef SOLVE3D
     &                                u, v,                             &
     &                                ad_t, ad_u, ad_v,                 &
     &                                tl_t, tl_u, tl_v,                 &
#  endif
     &                                ad_ubar, ad_vbar, ad_zeta,        &
     &                                tl_ubar, tl_vbar, tl_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
!
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      USE exchange_2d_mod
#   ifdef SOLVE3D
      USE exchange_3d_mod
#   endif
#  endif
#  ifdef SOLVE3D
      USE set_depth_mod, ONLY : set_depth_tile
      USE tl_set_depth_mod, ONLY : tl_set_depth_tile
#  endif
      USE tl_u2dbc_mod, ONLY : tl_u2dbc_tile
      USE tl_v2dbc_mod, ONLY : tl_v2dbc_tile
      USE tl_zetabc_mod, ONLY : tl_zetabc_tile
#  ifdef SOLVE3D
      USE tl_t3dbc_mod, ONLY : tl_t3dbc_tile
      USE tl_u3dbc_mod, ONLY : tl_u3dbc_tile
      USE tl_v3dbc_mod, ONLY : tl_v3dbc_tile
#  endif

!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Linp, Lout

#  ifdef ASSUMED_SHAPE
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#   endif
#   ifdef SOLVE3D
#    if defined SEDIMENT && defined SED_MORPH
      real(r8), intent(in) :: bed(LBi:,LBj:,:,:)
      real(r8), intent(in) :: tl_bed(LBi:,LBj:,:,:)
#    endif
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)
#    ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:,LBj:)
#    endif
#   endif
      real(r8), intent(in) :: ubar(LBi:,UBi:,:)
      real(r8), intent(in) :: vbar(LBi:,UBi:,:)
      real(r8), intent(in) :: zeta(LBi:,UBi:,:)
#   ifdef SOLVE3D
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
#   endif
      real(r8), intent(in) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: ad_vbar(LBi:,LBj:,:)
      real(r8), intent(in) :: ad_zeta(LBi:,LBj:,:)
#   ifdef SOLVE3D
      real(r8), intent(in) :: ad_t(LBi:,LBj:,:,:,:)
      real(r8), intent(in) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: ad_v(LBi:,LBj:,:,:)
#   endif
      real(r8), intent(inout) :: tl_h(LBi:,LBj:)
      real(r8), intent(inout) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_vbar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_zeta(LBi:,LBj:,:)
#   ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_v(LBi:,LBj:,:,:)
#    if defined SEDIMENT && defined SED_MORPH
      real(r8), intent(inout) :: bed_thick0(LBi:,LBj:)
#    endif
      real(r8), intent(inout) :: Hz(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_Hz(LBi:,LBj:,:)
      real(r8), intent(inout) :: Zt_avg1(LBi:,LBj:)
      real(r8), intent(inout) :: tl_Zt_avg1(LBi:,LBj:)
#   endif
#   ifdef SOLVE3D
      real(r8), intent(out) :: z_r(LBi:,LBj:,:)
      real(r8), intent(out) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(out) :: tl_z_r(LBi:,LBj:,:)
      real(r8), intent(out) :: tl_z_w(LBi:,LBj:,0:)
#   endif
#  else
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#   endif
#   ifdef SOLVE3D
#    if defined SEDIMENT && defined SED_MORPH
      real(r8), intent(in) :: bed(LBi:UBi,LBj:UBj,Nbed,MBEDP)
      real(r8), intent(in) :: tl_bed(LBi:UBi,LBj:UBj,Nbed,MBEDP)
#    endif
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)
#    ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:,LBj:)
#    endif
#   endif
      real(r8), intent(in) :: ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: zeta(LBi:UBi,LBj:UBj,3)
#   ifdef SOLVE3D
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(in) :: ad_u(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: ad_v(LBi:UBi,LBj:UBj,N(ng))
#   endif
      real(r8), intent(in) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: ad_vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: ad_zeta(LBi:UBi,LBj:UBj,3)
#   ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: tl_u(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: tl_v(LBi:UBi,LBj:UBj,N(ng))
#   endif
      real(r8), intent(inout) :: tl_h(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: tl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: tl_vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: tl_zeta(LBi:UBi,LBj:UBj,3)
#   ifdef SOLVE3D
#    if defined SEDIMENT && defined SED_MORPH
      real(r8), intent(inout) :: bed_thick0(LBi:UBi,LBj:UBj)
#    endif
      real(r8), intent(inout) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: tl_Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: Zt_avg1(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: tl_Zt_avg1(LBi:UBi,LBj:UBj)
#   endif
#   ifdef SOLVE3D
      real(r8), intent(out) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(out) :: tl_z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(out) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(out) :: tl_z_w(LBi:UBi,LBj:UBj,0:N(ng))
#   endif
#  endif
!
!  Local variable declarations.
!
      Logical :: FirstPass =.TRUE.

      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: NSUB, i, j
#  ifdef SOLVE3D
      integer :: itrc, k
#  endif
      integer, dimension(NstateVar(ng)+1) :: StateVar

      real(r8) :: my_dot, p, scale
#  ifdef SOLVE3D
      real(r8) :: cff, tl_cff

      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N(ng)) :: CF
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N(ng)) :: DC
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N(ng)) :: tl_CF
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N(ng)) :: tl_DC
#  endif
!
#  include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Initialize tangent linear with the steepest descent direction times
!  the perturbation amplitude "p".
!-----------------------------------------------------------------------
!
!  Set state variable to perturb according to outer loop index.
!
#  ifdef SOLVE3D
      StateVar(1)=0
      StateVar(2)=isFsur
      StateVar(3)=isUvel
      StateVar(4)=isVvel
      DO i=1,NT(ng)
        StateVar(4+i)=isTvar(i)
      END DO
#  else
      StateVar(1)=0
      StateVar(2)=isFsur
      StateVar(3)=isUbar
      StateVar(4)=isVbar
#  endif
!
!  Set perturbation amplitude as a function of the inner loop.
!
      p=10.0_r8**FLOAT(-inner)
      my_dot=0.0_r8
      scale=1.0_r8/SQRT(adDotProduct)
      IF (SOUTH_WEST_CORNER) THEN
        WRITE (stdout,10)
      END IF
!
!  Free-surface.
!
      IF ((StateVar(outer).eq.0).or.(StateVar(outer).eq.isFsur)) THEN
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            tl_zeta(i,j,Lout)=p*ad_zeta(i,j,Linp)*scale
#  ifdef MASKING
            tl_zeta(i,j,Lout)=tl_zeta(i,j,Lout)*rmask(i,j)
#  endif
            my_dot=my_dot+                                              &
     &             tl_zeta(i,j,Lout)*ad_zeta(i,j,Linp)
#  ifdef MASKING
            my_dot=my_dot*rmask(i,j)
#  endif
          END DO
        END DO
        IF (SOUTH_WEST_TEST) THEN
          IF (Master) THEN
            WRITE (stdout,20) outer, inner,                             &
     &                        TRIM(Vname(1,idSvar(isFsur)))
          END IF
        END IF
      END IF
#  ifndef SOLVE3D
!
!  2D u-momentum component.
!
      IF ((StateVar(outer).eq.0).or.(StateVar(outer).eq.isUbar)) THEN
        DO j=JstrR,JendR
          DO i=Istr,IendR
            tl_ubar(i,j,Lout)=p*ad_ubar(i,j,Linp)*scale
#   ifdef MASKING
            tl_ubar(i,j,Lout)=tl_ubar(i,j,Lout)*umask(i,j)
#   endif
            my_dot=my_dot+                                              &
     &             tl_ubar(i,j,Lout)*ad_ubar(i,j,Linp)
#   ifdef MASKING
            my_dot=my_dot*umask(i,j)
#   endif
          END DO
        END DO
        IF (SOUTH_WEST_TEST) THEN
          IF (Master) THEN
            WRITE (stdout,20) outer, inner,                             &
     &                        TRIM(Vname(1,idSvar(isUbar)))
          END IF
        END IF
      END IF
!
!  2D v-momentum component.
!
      IF ((StateVar(outer).eq.0).or.(StateVar(outer).eq.isVbar)) THEN
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            tl_vbar(i,j,Lout)=p*ad_vbar(i,j,Linp)*scale
#   ifdef MASKING
            tl_vbar(i,j,Lout)=tl_vbar(i,j,Lout)*vmask(i,j)
#   endif
            my_dot=my_dot+                                              &
     &             tl_vbar(i,j,Lout)*ad_vbar(i,j,Linp)
#   ifdef MASKING
            my_dot=my_dot*vmask(i,j)
#   endif
          END DO
        END DO
        IF (SOUTH_WEST_TEST) THEN
          IF (Master) THEN
            WRITE (stdout,20) outer, inner,                             &
     &                        TRIM(Vname(1,idSvar(isVbar)))
          END IF
        END IF
      END IF
#  endif
#  ifdef SOLVE3D
!
!  3D u-momentum component.
!
      IF ((StateVar(outer).eq.0).or.(StateVar(outer).eq.isUvel)) THEN
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=Istr,IendR
              tl_u(i,j,k,Lout)=p*ad_u(i,j,k,Linp)*scale
#   ifdef MASKING
              tl_u(i,j,k,Lout)=tl_u(i,j,k,Lout)*umask(i,j)
#   endif
              my_dot=my_dot+                                            &
     &               tl_u(i,j,k,Lout)*ad_u(i,j,k,Linp)
#   ifdef MASKING
              my_dot=my_dot*umask(i,j)
#   endif
            END DO
          END DO
        END DO
        IF (SOUTH_WEST_TEST) THEN
          IF (Master) THEN
            WRITE (stdout,20) outer, inner,                             &
     &                        TRIM(Vname(1,idSvar(isUvel)))
          END IF
        END IF
      END IF
!
!  3D v-momentum component.
!
      IF ((StateVar(outer).eq.0).or.(StateVar(outer).eq.isVvel)) THEN
        DO k=1,N(ng)
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              tl_v(i,j,k,Lout)=p*ad_v(i,j,k,Linp)*scale
#   ifdef MASKING
              tl_v(i,j,k,Lout)=tl_v(i,j,k,Lout)*vmask(i,j)
#   endif
              my_dot=my_dot+                                            &
     &               tl_v(i,j,k,Lout)*ad_v(i,j,k,Linp)
#   ifdef MASKING
              my_dot=my_dot*vmask(i,j)
#   endif
            END DO
          END DO
        END DO
        IF (SOUTH_WEST_TEST) THEN
          IF (Master) THEN
            WRITE (stdout,20) outer, inner,                             &
     &                        TRIM(Vname(1,idSvar(isVvel)))
          END IF
        END IF
      END IF
!
!  Tracers.
!
      DO itrc=1,NT(ng)
        IF ((StateVar(outer).eq.0).or.                                  &
     &      (StateVar(outer).eq.isTvar(itrc))) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                tl_t(i,j,k,Lout,itrc)=p*ad_t(i,j,k,Linp,itrc)*scale
#   ifdef MASKING
                tl_t(i,j,k,Lout,itrc)=tl_t(i,j,k,Lout,itrc)*rmask(i,j)
#   endif
                my_dot=my_dot+                                          &
     &                 tl_t(i,j,k,Lout,itrc)*ad_t(i,j,k,Linp,itrc)
#   ifdef MASKING
                my_dot=my_dot*rmask(i,j)
#   endif
              END DO
            END DO
          END DO
          IF (SOUTH_WEST_TEST) THEN
            IF (Master) THEN
              WRITE (stdout,20) outer, inner,                           &
     &                          TRIM(Vname(1,idSvar(isTvar(itrc))))
            END IF
          END IF
        END IF
      END DO
#  endif
      IF (Master) WRITE (stdout,'(/)')
!
!  Perform parallel global reduction operation: dot product.
!
      IF (SOUTH_WEST_CORNER.and.                                        &
     &    NORTH_EAST_CORNER) THEN
        NSUB=1                           ! non-tiled application
      ELSE
        NSUB=NtileX(ng)*NtileE(ng)       ! tiled application
      END IF
!$OMP CRITICAL (TL_DOT)
      IF (tile_count.eq.0) THEN
        DotProduct=0.0_r8
      END IF
      DotProduct=DotProduct+my_dot
      tile_count=block_count+1
      IF (tile_count.eq.NSUB) THEN
        tile_count=0
#  ifdef DISTRIBUTE
        op_handle='SUM'
        CALL mp_reduce (ng, iTLM, 1, DotProduct, op_handle)
#  endif
      END IF
!
!-----------------------------------------------------------------------
!  Apply lateral boundary conditions to 2D fields.
!-----------------------------------------------------------------------
!
      CALL tl_zetabc_tile (ng, Istr, Iend, Jstr, Jend,                  &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     Lout, Lout, Lout,                            &
     &                     zeta, tl_zeta)
#  ifndef SOLVE3D
      CALL tl_u2dbc_tile (ng, Istr, Iend, Jstr, Jend,                   &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    Lout, Lout, Lout,                             &
     &                    ubar, vbar, zeta,                             &
     &                    tl_ubar, tl_vbar, tl_zeta)
      CALL tl_v2dbc_tile (ng, Istr, Iend, Jstr, Jend,                   &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    Lout, Lout, Lout,                             &
     &                    ubar, vbar, zeta,                             &
     &                    tl_ubar, tl_vbar, tl_zeta)
#  endif
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      CALL exchange_r2d_tile (ng, iTLM, Istr, Iend, Jstr, Jend,         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        tl_zeta(:,:,Lout))
#   ifndef SOLVE3D
      CALL exchange_u2d_tile (ng, iTLM, Istr, Iend, Jstr, Jend,         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        tl_ubar(:,:,Lout))
      CALL exchange_v2d_tile (ng, iTLM, Istr, Iend, Jstr, Jend,         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        tl_vbar(:,:,Lout))
#   endif
#  endif
#  ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Compute new depths and thicknesses.
!-----------------------------------------------------------------------
!
      CALL set_depth_tile (ng, Istr, Iend, Jstr, Jend,                  &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     Lout, FirstPass,                             &
     &                     zeta, h,                                     &
#   ifdef ICESHELF
     &                     zice,                                        &
#   endif
#   if defined SEDIMENT && defined SED_MORPH
     &                     bed, bed_thick0,                             &
#   endif
     &                     Zt_avg1,                                     &
     &                     Hz, z_r, z_w)
      CALL tl_set_depth_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Lout, FirstPass,                          &
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
!
!-----------------------------------------------------------------------
!  Apply lateral boundary conditions to perturbed 3D fields.
!-----------------------------------------------------------------------
!
      CALL tl_u3dbc_tile (ng, Istr, Iend, Jstr, Jend,                   &
     &                    LBi, UBi, LBj, UBj, N(ng),                    &
     &                    Lout, Lout, tl_u)
      CALL tl_v3dbc_tile (ng, Istr, Iend, Jstr, Jend,                   &
     &                    LBi, UBi, LBj, UBj, N(ng),                    &
     &                    Lout, Lout, tl_v)
      DO itrc=1,NT(ng)
        CALL tl_t3dbc_tile (ng, Istr, Iend, Jstr, Jend, itrc,           &
     &                      LBi, UBi, LBj, UBj, N(ng), NT(ng),          &
     &                      Lout, Lout, tl_t)
      END DO
#   if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      CALL exchange_u3d_tile (ng, iTLM, Istr, Iend, Jstr, Jend,         &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        NghostPoints,                             &
     &                        tl_u(:,:,:,Lout))
      CALL exchange_v3d_tile (ng, iTLM, Istr, Iend, Jstr, Jend,         &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        NghostPoints,                             &
     &                        tl_v(:,:,:,Lout))
      DO itrc=1,NT(ng)
        CALL exchange_r3d_tile (ng, iTLM, Istr, Iend, Jstr, Jend,       &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          tl_t(:,:,:,Lout,itrc))
      END DO
#   endif
!
!-----------------------------------------------------------------------
!  Compute vertically integrated momentum (ubar,vbar) from perturbed
!  (u,v).
!-----------------------------------------------------------------------
!
!  Couple velocity component in the XI-direction.
!
      IF ((StateVar(outer).eq.0).or.(StateVar(outer).eq.isUvel)) THEN
        DO j=JU_RANGE
          DO i=IU_RANGE
            DC(i,0)=0.0_r8
            CF(i,0)=0.0_r8
            tl_DC(i,0)=0.0_r8
            tl_CF(i,0)=0.0_r8
          END DO
!
!  Compute thicknesses of U-boxes DC(i,1:N), total depth of the water
!  column DC(i,0), and incorrect vertical mean CF(i,0).
!
          DO k=1,N(ng)
            DO i=IU_RANGE
              DC(i,k)=0.5_r8*(Hz(i,j,k)+Hz(i-1,j,k))*on_u(i,j)
              DC(i,0)=DC(i,0)+DC(i,k)
              CF(i,0)=CF(i,0)+DC(i,k)*u(i,j,k,Linp)
              tl_DC(i,k)=0.5_r8*(tl_Hz(i,j,k)+tl_Hz(i-1,j,k))*on_u(i,j)
              tl_DC(i,0)=tl_DC(i,0)+tl_DC(i,k)
              tl_CF(i,0)=tl_CF(i,0)+tl_DC(i,k)*u(i,j,k,Linp)+           &
     &                              DC(i,k)*tl_u(i,j,k,Lout)
            END DO
          END DO
          DO i=IU_RANGE
            cff=1.0_r8/DC(i,0)
            tl_cff=-cff*cff*tl_DC(i,0)
            tl_ubar(i,j,Lout)=tl_CF(i,0)*cff+CF(i,0)*tl_cff
#   ifdef MASKING
            tl_ubar(i,j,Lout)=tl_ubar(i,j,Lout)*umask(i,j)
#   endif
          END DO
        END DO
        CALL tl_u2dbc_tile (ng, Istr, Iend, Jstr, Jend,                 &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      Lout, Lout, Lout,                           &
     &                      ubar, vbar, zeta,                           &
     &                      tl_ubar, tl_vbar, tl_zeta)

#   if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
        CALL exchange_u2d_tile (ng, model, Istr, Iend, Jstr, Jend,      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          tl_ubar(:,:,Lout))
#   endif
      END IF
!
!  Couple velocity component in the ETA-direction.
!
      IF ((StateVar(outer).eq.0).or.(StateVar(outer).eq.isVvel)) THEN
        DO j=JV_RANGE
          DO i=IV_RANGE
            DC(i,0)=0.0_r8
            CF(i,0)=0.0_r8
            tl_DC(i,0)=0.0_r8
            tl_CF(i,0)=0.0_r8
          END DO
!
!  Compute thicknesses of V-boxes DC(i,1:N), total depth of the water
!  column DC(i,0), and incorrect vertical mean CF(i,0).
!
          DO k=1,N(ng)
            DO i=IV_RANGE
              DC(i,k)=0.5_r8*(Hz(i,j,k)+Hz(i,j-1,k))*om_v(i,j)
              DC(i,0)=DC(i,0)+DC(i,k)
              CF(i,0)=CF(i,0)+DC(i,k)*v(i,j,k,Linp)
              tl_DC(i,k)=0.5_r8*(tl_Hz(i,j,k)+tl_Hz(i,j-1,k))*om_v(i,j)
              tl_DC(i,0)=tl_DC(i,0)+tl_DC(i,k)
              tl_CF(i,0)=tl_CF(i,0)+tl_DC(i,k)*v(i,j,k,Linp)+           &
     &                              DC(i,k)*tl_v(i,j,k,Lout)
            END DO
          END DO
          DO i=IV_RANGE
            cff=1.0_r8/DC(i,0)
            tl_cff=-cff*cff*tl_DC(i,0)
            tl_vbar(i,j,Lout)=tl_CF(i,0)*cff+CF(i,0)*tl_cff
#   ifdef MASKING
            tl_vbar(i,j,Lout)=tl_vbar(i,j,Lout)*vmask(i,j)
#   endif
          END DO
        END DO
        CALL tl_v2dbc_tile (ng, Istr, Iend, Jstr, Jend,                 &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      Lout, Lout, Lout,                           &
     &                      ubar, vbar, zeta,                           &
     &                      tl_ubar, tl_vbar, tl_zeta)
#   if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
        CALL exchange_v2d_tile (ng, iTLM, Istr, Iend, Jstr, Jend,       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          tl_vbar(:,:,Lout))
#   endif
      END IF
#  endif
!
 10   FORMAT (/,' Perturbing Tangent Linear Initial Conditions:',/)
 20   FORMAT (' (Outer,Inner) = ','(',i4.4,',',i4.4,')',3x,             &
     &        'Perturbing State Variable: ',a)

      RETURN
      END SUBROUTINE tl_ini_perturb_tile
# endif
# undef IU_RANGE
# undef IV_RANGE
# undef JU_RANGE
# undef JV_RANGE
#endif
      END MODULE ini_adjust_mod
