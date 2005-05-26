#include "cppdefs.h"
# ifdef FULL_GRID
#  define IR_RANGE IstrR,IendR
#  define IU_RANGE Istr,IendR
#  define JR_RANGE JstrR,JendR
#  define JV_RANGE Jstr,JendR
# else
#  define IR_RANGE Istr,Iend
#  define IU_RANGE IstrU,Iend
#  define JR_RANGE Jstr,Jend
#  define JV_RANGE JstrV,Jend
# endif

      MODULE packing_mod

#ifdef PROPAGATOR
!
!=================================================== Andrew M. Moore ===
!  Copyright (c) 2005 ROMS/TOMS Adjoint Group                          !
!================================================== Hernan G. Arango ===
!                                                                      !
!  These routines pack and unpack model state varaibles into/from a    !
!  single vector to interface with  ARPACKs Arnoldi Method  for the    !
!  computation Ritz eigenfunctions.                                    !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC :: tl_unpack
      PUBLIC :: ad_pack

      CONTAINS

      SUBROUTINE tl_unpack (ng, tile, Nstate, state)
!
!=======================================================================
!                                                                      !
!  This routine unpacks the tangent linear variables from the state    !
!  vector.  If applicable,  the state vector includes only unmasked    !
!  water points.  In 3D applications,  the  2D momentum is computed    !
!  by vertically integrating 3D momentum.                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
# ifdef SOLVE3D
      USE mod_coupling
# endif
      USE mod_grid
      USE mod_ocean
      USE mod_stepping
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: Nstate
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: state(:)
# else
      real(r8), intent(in) :: state(Nstate)
# endif
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iTLM, 1)
# endif
      CALL tl_unpack_tile (ng, Istr, Iend, Jstr, Jend,                  &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     kstp(ng),                                    &
# ifdef SOLVE3D
     &                     nstp(ng),                                    &
# endif
     &                     Nstate, state,                               &
# ifdef MASKING
     &                     GRID(ng) % IJwaterR,                         &
     &                     GRID(ng) % IJwaterU,                         &
     &                     GRID(ng) % IJwaterV,                         &
     &                     GRID(ng) % rmask,                            &
     &                     GRID(ng) % umask,                            &
     &                     GRID(ng) % vmask,                            &
# endif
# ifdef SOLVE3D
#  if defined SEDIMENT && defined SED_MORPH
     &                     OCEAN(ng) % tl_bed,                          &
     &                     GRID(ng) % tl_bed_thick0,                    &
#  endif
     &                     GRID(ng) % h,                                &
     &                     GRID(ng) % tl_h,                             &
     &                     GRID(ng) % Hz,                               &
     &                     GRID(ng) % tl_Hz,                            &
     &                     GRID(ng) % tl_z_r,                           &
     &                     GRID(ng) % tl_z_w,                           &
     &                     COUPLING(ng) % Zt_avg1,                      &
     &                     COUPLING(ng) % tl_Zt_avg1,                   &
#  ifdef ICESHELF
     &                     GRID(ng) % zice,                             &
#  endif
# endif
# ifdef SOLVE3D
     &                     OCEAN(ng) % tl_t,                            &
     &                     OCEAN(ng) % u,                               &
     &                     OCEAN(ng) % tl_u,                            &
     &                     OCEAN(ng) % v,                               &
     &                     OCEAN(ng) % tl_v,                            &
# endif
     &                     OCEAN(ng) % ubar,                            &
     &                     OCEAN(ng) % tl_ubar,                         &
     &                     OCEAN(ng) % vbar,                            &
     &                     OCEAN(ng) % tl_vbar,                         &
     &                     OCEAN(ng) % zeta,                            &
     &                     OCEAN(ng) % tl_zeta)
# ifdef PROFILE
      CALL wclock_off (ng, iTLM, 1)
# endif
      RETURN
      END SUBROUTINE tl_unpack      
!
!***********************************************************************
      SUBROUTINE tl_unpack_tile (ng, Istr, Iend, Jstr, Jend,            &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           kstp,                                  &
# ifdef SOLVE3D
     &                           nstp,                                  &
# endif
     &                           Nstate, state,                         &
# ifdef MASKING
     &                           IJwaterR, IJwaterU, IJwaterV,          &
     &                           rmask, umask, vmask,                   &
# endif
# ifdef SOLVE3D
#  if defined SEDIMENT && defined SED_MORPH
     &                           tl_bed, tl_bed_thick0,                 &
#  endif
     &                           h, tl_h,                               &
     &                           Hz, tl_Hz,                             &
     &                           tl_z_r, tl_z_w,                        &
     &                           Zt_avg1, tl_Zt_avg1,                   &
#  ifdef ICESHELF
     &                           zice,                                  &
#  endif
# endif
# ifdef SOLVE3D
     &                           tl_t,                                  &
     &                           u, tl_u,                               &
     &                           v, tl_v,                               &
# endif
     &                           ubar, tl_ubar,                         &
     &                           vbar, tl_vbar,                         &
     &                           zeta, tl_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_ncparam
      USE mod_scalars
!
# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      USE exchange_2d_mod
#  ifdef SOLVE3D
      USE exchange_3d_mod
#  endif
# endif
# ifdef SOLVE3D
      USE tl_set_depth_mod, ONLY : tl_set_depth_tile
# endif
      USE tl_u2dbc_mod, ONLY : tl_u2dbc_tile
      USE tl_v2dbc_mod, ONLY : tl_v2dbc_tile
      USE tl_zetabc_mod, ONLY : tl_zetabc_tile
# ifdef SOLVE3D
      USE tl_t3dbc_mod, ONLY : tl_t3dbc_tile
      USE tl_u3dbc_mod, ONLY : tl_u3dbc_tile
      USE tl_v3dbc_mod, ONLY : tl_v3dbc_tile
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: kstp
# ifdef SOLVE3D
      integer, intent(in) :: nstp
# endif
      integer, intent(in) :: Nstate
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      integer, intent(in) :: IJwaterR(LBi:,LBj:)
      integer, intent(in) :: IJwaterU(LBi:,LBj:)
      integer, intent(in) :: IJwaterV(LBi:,LBj:)

      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: state(:)
#  ifdef SOLVE3D
#   if defined SEDIMENT && defined SED_MORPH
      real(r8), intent(in) :: tl_bed(LBi:,LBj:,:,:)
      real(r8), intent(in) :: tl_bed_thick0(LBi:,LBj:)
#   endif
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
#   ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:,LBj:)
#   endif
      real(r8), intent(inout) :: Zt_avg1(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: vbar(LBi:,LBj:,:)
      real(r8), intent(in) :: zeta(LBi:,LBj:,:)
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
#  ifdef SOLVE3D
      real(r8), intent(inout) :: tl_Hz(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_h(LBi:,LBj:)
      real(r8), intent(inout) :: tl_Zt_avg1(LBi:,LBj:)
      real(r8), intent(inout) :: tl_z_r(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_z_w(LBi:,LBj:,0:)
#  endif
      real(r8), intent(inout) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_vbar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_zeta(LBi:,LBj:,:)
#  ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_v(LBi:,LBj:,:,:)
#  endif
# else
#  ifdef MASKING
      integer, intent(in) :: IJwaterR(LBi:UBi,LBj:UBj)
      integer, intent(in) :: IJwaterU(LBi:UBi,LBj:UBj)
      integer, intent(in) :: IJwaterV(LBi:UBi,LBj:UBj)

      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: state(Nstate)
#  ifdef SOLVE3D
#   if defined SEDIMENT && defined SED_MORPH
      real(r8), intent(in) :: tl_bed(LBi:UBi,LBj:UBj,Nbed,MBEDP)
      real(r8), intent(in) :: tl_bed_thick0(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
#   ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(inout) :: Zt_avg1(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: zeta(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),2)
#  ifdef SOLVE3D
      real(r8), intent(inout) :: tl_Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: tl_h(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: tl_Zt_avg1(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: tl_z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: tl_z_w(LBi:UBi,LBj:UBj,0:N(ng))
#  endif
      real(r8), intent(inout) :: tl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: tl_vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: tl_zeta(LBi:UBi,LBj:UBj,3)
#  ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: tl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: tl_v(LBi:UBi,LBj:UBj,N(ng),2)
#  endif
# endif
!
!  Local variable declarations.
!
# ifdef SOLVE3D
      logical :: FirstPass = .TRUE.
# endif
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
# ifndef MASKING
      integer :: Imax, Ioff, Jmax, Joff
# endif
      integer :: i, iadd, is, itrc, j, k

      integer, dimension(5+NT(ng)) :: offset

      real(r8) :: cff, tl_cff, scale

# ifdef SOLVE3D
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY) :: CF
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY) :: DC

      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY) :: tl_CF
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY) :: tl_DC
# endif

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Extract state variables from full 1D state vector.
!-----------------------------------------------------------------------
!
!  Determine the index offset for each variable in the state vector.
# ifdef MASKING
!  Notice that in Land/Sea masking application the state vector only
!  contains water points to avoid large null space.
# endif
!
# ifdef SOLVE3D
#  ifdef MASKING
      offset(isFsur)=0
      offset(isUvel)=offset(isFsur)+NwaterR(ng)
      offset(isVvel)=offset(isUvel)+NwaterU(ng)*N(ng)
      iadd=NwaterV(ng)*N(ng)
      DO itrc=1,NT(ng)
        offset(isTvar(itrc))=offset(isTvar(itrc)-1)+iadd
        iadd=NwaterR(ng)*N(ng)
      END DO
#  else
#   ifdef FULL_GRID
      offset(isFsur)=0
      offset(isUvel)=offset(isFsur)+(Lm(ng)+2)*(Mm(ng)+2)
      offset(isVvel)=offset(isUvel)+(Lm(ng)+1)*(Mm(ng)+2)*N(ng)
      iadd=(Lm(ng)+2)*(Mm(ng)+1)*N(ng)
      DO itrc=1,NT(ng)
        offset(isTvar(itrc))=offset(isTvar(itrc)-1)+iadd
        iadd=(Lm(ng)+2)*(Mm(ng)+2)*N(ng)
      END DO
#   else
      offset(isFsur)=0
      offset(isUvel)=offset(isFsur)+Lm(ng)*Mm(ng)
      offset(isVvel)=offset(isUvel)+(Lm(ng)-1)*Mm(ng)*N(ng)
      iadd=Lm(ng)*(Mm(ng)-1)*N(ng)
      DO itrc=1,NT(ng)
        offset(isTvar(itrc))=offset(isTvar(itrc)-1)+iadd
        iadd=Lm(ng)*Mm(ng)*N(ng)
      END DO
#   endif
#  endif
# else
#  ifdef MASKING
      offset(isFsur)=0
      offset(isUbar)=offset(isFsur)+NwaterR(ng)
      offset(isVbar)=offset(isUbar)+NwaterU(ng)
#  else
#   ifdef FULL_GRID
      offset(isFsur)=0
      offset(isUbar)=offset(isFsur)+(Lm(ng)+2)*(Mm(ng)+2)
      offset(isVbar)=offset(isUbar)+(Lm(ng)+1)*(Mm(ng)+2)
#   else
      offset(isFsur)=0
      offset(isUbar)=offset(isFsur)+Lm(ng)*Mm(ng)
      offset(isVbar)=offset(isUbar)+(Lm(ng)-1)*Mm(ng)
#   endif
#  endif
# endif
!
!  Unpack tangent linear free-surface.
!
# ifndef MASKING
#  ifdef FULL_GRID
      Imax=Lm(ng)+2
      Ioff=1
      Joff=0
#  else
      Imax=Lm(ng)
      Ioff=0
      Joff=1
#  endif
# endif
      scale=1.0_r8/SQRT(0.5_r8*g*rho0)
      DO j=JR_RANGE
        DO i=IR_RANGE
# ifdef MASKING
          IF (rmask(i,j).gt.0.0_r8) THEN
            is=IJwaterR(i,j)+offset(isFsur)
            tl_zeta(i,j,kstp)=scale*state(is)
          ELSE
            tl_zeta(i,j,kstp)=0.0_r8
          END IF
# else
          is=(i+Ioff)+(j-Joff)*Imax+offset(isFsur)
          tl_zeta(i,j,kstp)=cff*state(is)
# endif
        END DO
      END DO
# ifndef SOLVE3D
!
!  Unpack tangent linear 2D U-velocity.
!
#  ifndef MASKING
#   ifdef FULL_GRID
      Imax=Lm(ng)+1
      Ioff=0
      Joff=0
#   else
      Imax=Lm(ng)-1
      Ioff=1
      Joff=1
#   endif
#  endif
      cff=0.25_r8*rho0
      DO j=JR_RANGE
        DO i=IU_RANGE
          scale=1.0_r8/SQRT(cff*(h(i-1,j)+h(i,j)))
#  ifdef MASKING
          IF (umask(i,j).gt.0.0_r8) THEN
            is=IJwaterU(i,j)+offset(isUbar)
            tl_ubar(i,j,kstp)=scale*state(is)
          ELSE
            tl_ubar(i,j,kstp)=0.0_r8
          END IF
#  else
          is=(i-Ioff)+(j-Joff)*Imax+offset(isUbar)
          tl_ubar(i,j,kstp)=scale*state(is)
#  endif
        END DO
      END DO
!
!  Unpack tangent linear 2D V-velocity.
!
#  ifndef MASKING
#   ifdef FULL_GRID
      Imax=Lm(ng)+2
      Ioff=1
      Joff=1
#   else
      Imax=Lm(ng)
      Ioff=0
      Joff=2
#   endif
#  endif
      cff=0.25_r8*rho0
      DO j=JV_RANGE
        DO i=IR_RANGE
          scale=1.0_r8/SQRT(cff*(h(i,j-1)+h(i,j)))
#  ifdef MASKING
          IF (vmask(i,j).gt.0.0_r8) THEN
            is=IJwaterV(i,j)+offset(isVbar)
            tl_vbar(i,j,kstp)=scale*state(is)
          ELSE
            tl_vbar(i,j,kstp)=0.0_r8
          END IF
#  else
          is=(i+Ioff)+(j-Joff)*Imax+offset(isVbar)
          tl_vbar(i,j,kstp)=scale*state(is)
#  endif
        END DO
      END DO
# else
!
!  Unpack tangent linear 3D U-velocity.
!
#  ifndef MASKING
#   ifdef FULL_GRID
      Imax=Lm(ng)+1
      Jmax=Mm(ng)+2
      Ioff=0
      Joff=0
#   else
      Imax=Lm(ng)-1
      Ioff=1
      Joff=1
#   endif
#  endif
      cff=0.25_r8*rho0
      DO k=1,N(ng)
#  ifdef MASKING
        iadd=(k-1)*NwaterU(ng)+offset(isUvel)
#  else
        iadd=(k-1)*Imax*Jmax+offset(isUvel)
#  endif
        DO j=JR_RANGE
          DO i=IU_RANGE
#  ifdef MASKING
            IF (umask(i,j).gt.0.0_r8) THEN
              scale=1.0_r8/SQRT(cff*(Hz(i-1,j,k)+Hz(i,j,k)))
              is=IJwaterU(i,j)+iadd
              tl_u(i,j,k,nstp)=scale*state(is)
            ELSE
              tl_u(i,j,k,nstp)=0.0_r8
            END IF
#  else
            scale=1.0_r8/SQRT(cff*(Hz(i-1,j,k)+Hz(i,j,k)))
            is=(i-Ioff)+(j-Joff)*Imax+iadd
            tl_u(i,j,k,nstp)=scale*state(is)
#  endif
          END DO
        END DO
      END DO
!
!  Unpack tangent linear 3D V-velocity.
!
#  ifndef MASKING
#   ifdef FULL_GRID
      Imax=Lm(ng)+2
      Jmax=Mm(ng)+1
      Ioff=1
      Joff=1
#   else
      Imax=Lm(ng)
      Ioff=0
      Joff=2
#   endif
#  endif
      cff=0.25_r8*rho0
      DO k=1,N(ng)
#  ifdef MASKING
        iadd=(k-1)*NwaterV(ng)+offset(isVvel)
#  else
        iadd=(k-1)*Imax*Jmax+offset(isVvel)
#  endif
        DO j=JV_RANGE
          DO i=IR_RANGE
#  ifdef MASKING
            IF (vmask(i,j).gt.0.0_r8) THEN
              scale=1.0_r8/SQRT(cff*(Hz(i,j-1,k)+Hz(i,j,k)))
              is=IJwaterV(i,j)+iadd
              tl_v(i,j,k,nstp)=scale*state(is)
            ELSE
              tl_v(i,j,k,nstp)=0.0_r8
            END IF
#  else
            scale=1.0_r8/SQRT(cff*(Hz(i,j-1,k)+Hz(i,j,k)))
            is=(i+Ioff)+(j-Joff)*Imax+iadd
            tl_v(i,j,k,nstp)=scale*state(is)
#  endif
          END DO
        END DO
      END DO
!
!  Unpack tangent linear tracers variables. For now, use salinity scale
!  for passive tracers.
!
# ifndef MASKING
#  ifdef FULL_GRID
      Imax=Lm(ng)+2
      Jmax=Mm(ng)+2
      Ioff=1
      Joff=0
#  else
      Imax=Lm(ng)
      Ioff=0
      Joff=1
#  endif
# endif
      DO itrc=1,NT(ng)
        IF (itrc.eq.itemp) THEN
          cff=0.5_r8*rho0*Tcoef(ng)*Tcoef(ng)*g*g/bvf_bak
        ELSE IF (itrc.eq.isalt) THEN
          cff=0.5_r8*rho0*Scoef(ng)*Scoef(ng)*g*g/bvf_bak
        ELSE
          cff=0.5_r8*rho0*Scoef(ng)*Scoef(ng)*g*g/bvf_bak
        END IF
        DO k=1,N(ng)
#  ifdef MASKING
          iadd=(k-1)*NwaterR(ng)+offset(isTvar(itrc))
#  else
          iadd=(k-1)*Imax*Jmax+offset(isTvar(itrc))
#  endif
          DO j=JR_RANGE
            DO i=IR_RANGE
#  ifdef MASKING
              IF (rmask(i,j).gt.0.0_r8) THEN
                scale=1.0_r8/SQRT(cff*Hz(i,j,k))
                is=IJwaterR(i,j)+iadd
                tl_t(i,j,k,nstp,itrc)=scale*state(is)
              ELSE
                tl_t(i,j,k,nstp,itrc)=0.0_r8
              END IF
#  else
              scale=1.0_r8/SQRT(cff*Hz(i,j,k))
              is=(i+Ioff)+(j-Joff)*Imax+iadd
              tl_t(i,j,k,nstp,itrc)=scale*state(is)
#  endif
            END DO
          END DO
        END DO
      END DO
# endif
!
!-----------------------------------------------------------------------
!  Apply lateral boundary conditions to 2D fields.
!-----------------------------------------------------------------------
!
      CALL tl_zetabc_tile (ng, Istr, Iend, Jstr, Jend,                  &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     kstp, kstp, kstp,                            &
     &                     zeta, tl_zeta)
# ifndef SOLVE3D
      CALL tl_u2dbc_tile (ng, Istr, Iend, Jstr, Jend,                   &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    kstp, kstp, kstp,                             &
     &                    ubar, vbar, zeta,                             &
     &                    tl_ubar, tl_vbar, tl_zeta)
      CALL tl_v2dbc_tile (ng, Istr, Iend, Jstr, Jend,                   &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    kstp, kstp, kstp,                             &
     &                    ubar, vbar, zeta,                             &
     &                    tl_ubar, tl_vbar, tl_zeta)
# endif
# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      CALL exchange_r2d_tile (ng, iTLM, Istr, Iend, Jstr, Jend,         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        tl_zeta(:,:,kstp))
#  ifndef SOLVE3D
      CALL exchange_u2d_tile (ng, iTLM, Istr, Iend, Jstr, Jend,         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        tl_ubar(:,:,kstp))
      CALL exchange_v2d_tile (ng, iTLM, Istr, Iend, Jstr, Jend,         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        tl_vbar(:,:,kstp))
#  endif
# endif
# ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Compute tangent linear depths and thicknesses.
!-----------------------------------------------------------------------
!
      CALL tl_set_depth_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        kstp, FirstPass,                          &
     &                        tl_zeta,                                  &
     &                        h, tl_h,                                  &
#  ifdef ICESHELF
     &                        zice,                                     &
#  endif
#  if defined SEDIMENT && defined SED_MORPH
     &                        tl_bed, bed_thick0,                       &
#  endif
     &                        Zt_avg1, tl_Zt_avg1,                      &
     &                        tl_Hz, tl_z_r, tl_z_w)
!
!-----------------------------------------------------------------------
!  Apply lateral boundary conditions to 3D fields.
!-----------------------------------------------------------------------
!
      CALL tl_u3dbc_tile (ng, Istr, Iend, Jstr, Jend,                   &
     &                    LBi, UBi, LBj, UBj, N(ng),                    &
     &                    nstp, nstp, tl_u)
      CALL tl_v3dbc_tile (ng, Istr, Iend, Jstr, Jend,                   &
     &                    LBi, UBi, LBj, UBj, N(ng),                    &
     &                    nstp, nstp, tl_v)
      DO itrc=1,NT(ng)
        CALL tl_t3dbc_tile (ng, Istr, Iend, Jstr, Jend, itrc,           &
     &                      LBi, UBi, LBj, UBj, N(ng), NT(ng),          &
     &                      nstp, nstp, tl_t)
      END DO
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      CALL exchange_u3d_tile (ng, iTLM, Istr, Iend, Jstr, Jend,         &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        NghostPoints,                             &
     &                        tl_u(:,:,:,nstp))
      CALL exchange_v3d_tile (ng, iTLM, Istr, Iend, Jstr, Jend,         &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        NghostPoints,                             &
     &                        tl_v(:,:,:,nstp))
      DO itrc=1,NT(ng)
        CALL exchange_r3d_tile (ng, iTLM, Istr, Iend, Jstr, Jend,       &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          tl_t(:,:,:,nstp,itrc))
      END DO
#  endif
!
!-----------------------------------------------------------------------
!  Compute tangent linear 2D momentum by vertically integrating 
!  3D momentum.
!-----------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO i=IstrU,Iend
          CF(i)=0.0_r8
          tl_CF(i)=0.0_r8
          DC(i)=0.0_r8
          tl_DC(i)=0.0_r8
        END DO
        DO k=1,N(ng)
          DO i=IstrU,Iend
            cff=0.5_r8*(Hz(i-1,j,k)+Hz(i,j,k))
            tl_cff=0.5_r8*(tl_Hz(i-1,j,k)+tl_Hz(i,j,k))
            CF(i)=CF(i)+cff
            tl_CF(i)=tl_CF(i)+tl_cff
            DC(i)=DC(i)+u(i,j,k,nstp)*cff
            tl_DC(i)=tl_DC(i)+tl_u(i,j,k,nstp)*cff+                     &
     &                        u(i,j,k,nstp)*tl_cff
          END DO
        END DO
        DO i=IstrU,Iend
          cff=1.0_r8/CF(i)
          tl_cff=-cff*cff*tl_CF(i)
!>        ubar(i,j,kstp)=DC(i)*cff
!>
          tl_ubar(i,j,kstp)=tl_DC(i)*cff+DC(i)*tl_cff
#  ifdef MASKING
!>        ubar(i,j,kstp)=ubar(i,j,kstp)*umask(i,j)
!>
          tl_ubar(i,j,kstp)=tl_ubar(i,j,kstp)*umask(i,j)
#  endif
        END DO
      END DO
!
      DO j=JstrV,Jend
        DO i=Istr,Iend
          CF(i)=0.0_r8
          tl_CF(i)=0.0_r8
          DC(i)=0.0_r8
          tl_DC(i)=0.0_r8
        END DO
        DO k=1,N(ng)
          DO i=Istr,Iend
            cff=0.5_r8*(Hz(i,j-1,k)+Hz(i,j,k))
            tl_cff=0.5_r8*(tl_Hz(i,j-1,k)+tl_Hz(i,j,k))
            CF(i)=CF(i)+cff
            tl_CF(i)=tl_CF(i)+tl_cff
            DC(i)=DC(i)+v(i,j,k,nstp)*cff
            tl_DC(i)=tl_DC(i)+tl_v(i,j,k,nstp)*cff+                     &
     &                        v(i,j,k,nstp)*tl_cff
          END DO
        END DO
        DO i=Istr,Iend
          cff=1.0_r8/CF(i)
          tl_cff=-cff*cff*tl_CF(i)
!>        vbar(i,j,kstp)=DC(i)*cff
!>
          tl_vbar(i,j,kstp)=tl_DC(i)*cff+DC(i)*tl_cff
#  ifdef MASKING
!>        vbar(i,j,kstp)=vbar(i,j,kstp)*vmask(i,j)
!>
          tl_vbar(i,j,kstp)=tl_vbar(i,j,kstp)*vmask(i,j)
#  endif
        END DO
      END DO
!
!  Set lateral boundary conditions
!
      CALL tl_u2dbc_tile (ng, Istr, Iend, Jstr, Jend,                   &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    kstp, kstp, kstp,                             &
     &                    ubar, vbar, zeta,                             &
     &                    tl_ubar, tl_vbar, tl_zeta)
      CALL tl_v2dbc_tile (ng, Istr, Iend, Jstr, Jend,                   &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    kstp, kstp, kstp,                             &
     &                    ubar, vbar, zeta,                             &
     &                    tl_ubar, tl_vbar, tl_zeta)
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      CALL exchange_u2d_tile (ng, iTLM, Istr, Iend, Jstr, Jend,         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        tl_ubar(:,:,kstp))
      CALL exchange_v2d_tile (ng, iTLM, Istr, Iend, Jstr, Jend,         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        tl_vbar(:,:,kstp))
#  endif
# endif

      RETURN
      END SUBROUTINE tl_unpack_tile

      SUBROUTINE ad_pack (ng, tile, Nstate, ad_state)
!
!=======================================================================
!                                                                      !
!  This routine packs the adjoint variables into the state vector.     !
!  The state vector contains only interior water points.  Scale by     !
!  the inverse of the energy norm.                                     !
!                                                                      !
!=======================================================================
!
      USE mod_param
# ifdef SOLVE3D
      USE mod_coupling
# endif
      USE mod_grid
      USE mod_ocean
      USE mod_stepping
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: Nstate
# ifdef ASSUMED_SHAPE
      real(r8), intent(out) :: ad_state(:)
# else
      real(r8), intent(out) :: ad_state(Nstate)
# endif
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iADM, 1)
# endif
      CALL ad_pack_tile (ng, Istr, Iend, Jstr, Jend,                    &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   krhs(ng), kstp(ng), knew(ng),                  &
# ifdef SOLVE3D
     &                   nstp(ng),                                      &
# endif
     &                   Nstate, ad_state,                              &
# ifdef MASKING
     &                   GRID(ng) % IJwaterR,                           &
     &                   GRID(ng) % IJwaterU,                           &
     &                   GRID(ng) % IJwaterV,                           &
     &                   GRID(ng) % rmask,                              &
     &                   GRID(ng) % umask,                              &
     &                   GRID(ng) % vmask,                              &
# endif
# ifdef SOLVE3D
#  if defined SEDIMENT && defined SED_MORPH
     &                   OCEAN(ng) % ad_bed,                            &
     &                   GRID(ng) % ad_bed_thick0,                      &
#  endif
     &                   GRID(ng) % h,                                  &
     &                   GRID(ng) % ad_h,                               &
     &                   GRID(ng) % ad_z_r,                             &
     &                   GRID(ng) % ad_z_w,                             &
     &                   COUPLING(ng) % Zt_avg1,                        &
     &                   COUPLING(ng) % ad_Zt_avg1,                     &
#  ifdef ICESHELF
     &                   GRID(ng) % zice,                               &
#  endif
     &                   GRID(ng) % Hz,                                 &
     &                   GRID(ng) % ad_Hz,                              &
     &                   OCEAN(ng) % ad_t,                              &
     &                   OCEAN(ng) % u,                                 &
     &                   OCEAN(ng) % ad_u,                              &
     &                   OCEAN(ng) % v,                                 &
     &                   OCEAN(ng) % ad_v,                              &
# endif
     &                   OCEAN(ng) % ubar,                              &
     &                   OCEAN(ng) % ad_ubar,                           &
     &                   OCEAN(ng) % vbar,                              &
     &                   OCEAN(ng) % ad_vbar,                           &
     &                   OCEAN(ng) % zeta,                              &
     &                   OCEAN(ng) % ad_zeta)
# ifdef PROFILE
      CALL wclock_off (ng, iADM, 1)
# endif
      RETURN
      END SUBROUTINE ad_pack      
!
!***********************************************************************
      SUBROUTINE ad_pack_tile (ng, Istr, Iend, Jstr, Jend,              &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         krhs, kstp, knew,                        &
# ifdef SOLVE3D
     &                         nstp,                                    &
# endif
     &                         Nstate, ad_state,                        &
# ifdef MASKING
     &                         IJwaterR, IJwaterU, IJwaterV,            &
     &                         rmask, umask, vmask,                     &
# endif
# ifdef SOLVE3D
#  if defined SEDIMENT && defined SED_MORPH
     &                         ad_bed, ad_bed_thick0,                   &
#  endif
     &                         h, ad_h,                                 &
     &                         ad_z_r, ad_z_w,                          &
     &                         Zt_avg1, ad_Zt_avg1,                     &
#  ifdef ICESHELF
     &                         zice,                                    &
#  endif
     &                         Hz, ad_Hz,                               &
     &                         ad_t,                                    &
     &                         u, ad_u,                                 &
     &                         v, ad_v,                                 &
# endif
     &                         ubar, ad_ubar,                           &
     &                         vbar, ad_vbar,                           &
     &                         zeta, ad_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_ncparam
      USE mod_scalars
!
# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      USE exchange_2d_mod
#  ifdef SOLVE3D
      USE exchange_3d_mod
#  endif
# endif
# ifdef SOLVE3D
      USE ad_set_depth_mod, ONLY : ad_set_depth_tile
# endif
      USE ad_u2dbc_mod, ONLY : ad_u2dbc_tile
      USE ad_v2dbc_mod, ONLY : ad_v2dbc_tile
      USE ad_zetabc_mod, ONLY : ad_zetabc_tile
# ifdef SOLVE3D
      USE ad_t3dbc_mod, ONLY : ad_t3dbc_tile
      USE ad_u3dbc_mod, ONLY : ad_u3dbc_tile
      USE ad_v3dbc_mod, ONLY : ad_v3dbc_tile
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: krhs, kstp, knew
# ifdef SOLVE3D
      integer, intent(in) :: nstp
# endif 
      integer, intent(in) :: Nstate
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      integer, intent(in) :: IJwaterR(LBi:,LBj:)
      integer, intent(in) :: IJwaterU(LBi:,LBj:)
      integer, intent(in) :: IJwaterV(LBi:,LBj:)

      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
#  ifdef SOLVE3D
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
#   ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:,LBj:)
#   endif
      real(r8), intent(in) :: Zt_avg1(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: vbar(LBi:,LBj:,:)
      real(r8), intent(in) :: zeta(LBi:,LBj:,:)
#  ifdef SOLVE3D
#   if defined SEDIMENT && defined SED_MORPH
      real(r8), intent(inout) :: ad_bed(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_bed_thick0(LBi:,LBj:)
#   endif
      real(r8), intent(inout) :: ad_h(LBi:,LBj:)
      real(r8), intent(inout) :: ad_Zt_avg1(LBi:,LBj:)
      real(r8), intent(inout) :: ad_Hz(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_z_r(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_z_w(LBi:,LBj:,0:)
      real(r8), intent(inout) :: ad_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_v(LBi:,LBj:,:,:)
#  endif
      real(r8), intent(inout) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_vbar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_zeta(LBi:,LBj:,:)

      real(r8), intent(out) :: ad_state(:)
# else
#  ifdef MASKING
      integer, intent(in) :: IJwaterR(LBi:UBi,LBj:UBj)
      integer, intent(in) :: IJwaterU(LBi:UBi,LBj:UBj)
      integer, intent(in) :: IJwaterV(LBi:UBi,LBj:UBj)

      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef SOLVE3D
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
#   ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(in) :: Zt_avg1(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: zeta(LBi:UBi,LBj:UBj,3)
#  ifdef SOLVE3D
#   if defined SEDIMENT && defined SED_MORPH
      real(r8), intent(inout) :: ad_bed(LBi:UBi,LBj:UBj,Nbed,MBEDP)
      real(r8), intent(inout) :: ad_bed_thick0(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(inout) :: ad_h(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_Zt_avg1(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: ad_z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: ad_z_w(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(inout) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
#  endif
      real(r8), intent(inout) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_zeta(LBi:UBi,LBj:UBj,3)

      real(r8), intent(out) :: ad_state(Nstate)
# endif
!
!  Local variable declarations.
!
# ifdef SOLVE3D
      logical :: FirstPass = .TRUE.
# endif
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
# ifndef MASKING
      integer :: Imax, Ioff, Jmax, Joff
# endif
# ifdef SOLVE3D
      integer :: ILB, IUB
# endif
      integer :: i, iadd, is, itrc, j, k

      integer, dimension(5+NT(ng)) :: offset

      real(r8), parameter :: Aspv = 0.0_r8

      real(r8) :: cff, scale

# ifdef SOLVE3D
      real(r8) :: adfac, ad_cff

      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY) :: CF
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY) :: DC

      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY) :: ad_CF
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY) :: ad_DC
# endif

# include "set_bounds.h"

# ifdef SOLVE3D
      ILB=LBOUND(ad_CF,DIM=1)
      IUB=UBOUND(ad_CF,DIM=1)
!
!-----------------------------------------------------------------------
!  Initialize adjoint private variables.
!-----------------------------------------------------------------------
!
      ad_cff=0.0_r8

      ad_CF(ILB:IUB)=0.0_r8
      ad_DC(ILB:IUB)=0.0_r8
# endif
!
!-----------------------------------------------------------------------
!  Collect adjoint contributions from the initalization procedure.
!-----------------------------------------------------------------------
!
!  Collect contibutions from boundary exchanges. Notice that both times
!  levels "kstp" and "knew" are needed for 2D momentum.
!
# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      CALL ad_exchange_u2d_tile (ng, iADM, Istr, Iend, Jstr, Jend,      &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           NghostPoints,                          &
     &                           ad_ubar(:,:,kstp)
      CALL ad_exchange_v2d_tile (ng, iADM, Istr, Iend, Jstr, Jend,      &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           NghostPoints,                          &
     &                           ad_vbar(:,:,kstp)
      CALL ad_exchange_u2d_tile (ng, iADM, Istr, Iend, Jstr, Jend,      &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           NghostPoints,                          &
     &                           ad_ubar(:,:,knew)
      CALL ad_exchange_v2d_tile (ng, iADM, Istr, Iend, Jstr, Jend,      &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           NghostPoints,                          &
     &                           ad_vbar(:,:,knew)
#  ifdef SOLVE3D
      CALL ad_exchange_u3d_tile (ng, iADM, Istr, Iend, Jstr, Jend,      &
     &                           LBi, UBi, LBj, UBj, 1, N(ng),          &
     &                           NghostPoints,                          &
     &                           ad_u(:,:,:,nstp))
      CALL ad_exchange_v3d_tile (ng, iADM, Istr, Iend, Jstr, Jend,      &
     &                           LBi, UBi, LBj, UBj, 1, N(ng),          &
     &                           NghostPoints,                          &
     &                           ad_v(:,:,:,nstp))
      DO itrc=1,NT(ng)
        CALL ad_exchange_r3d_tile (ng, iADM, Istr, Iend, Jstr, Jend,    &
     &                             LBi, UBi, LBj, UBj, 1, N(ng),        &
     &                             NghostPoints,                        &
     &                             ad_t(:,:,:,nstp,itrc))
      END DO
#  endif
# endif
!
!  Collect contributions from lateral boundary conditons. Notice that
!  both time levels "kstp" and "knew" are needed for 2D momentum.
!
      CALL ad_u2dbc_tile (ng, Istr, Iend, Jstr, Jend,                   &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    krhs, kstp, kstp,                             &
     &                    ubar, vbar, zeta,                             &
     &                    ad_ubar, ad_vbar, ad_zeta)
      CALL ad_v2dbc_tile (ng, Istr, Iend, Jstr, Jend,                   &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    krhs, kstp, kstp,                             &
     &                    ubar, vbar, zeta,                             &
     &                    ad_ubar, ad_vbar, ad_zeta)
      CALL ad_u2dbc_tile (ng, Istr, Iend, Jstr, Jend,                   &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    krhs, kstp, knew,                             &
     &                    ubar, vbar, zeta,                             &
     &                    ad_ubar, ad_vbar, ad_zeta)
      CALL ad_v2dbc_tile (ng, Istr, Iend, Jstr, Jend,                   &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    krhs, kstp, knew,                             &
     &                    ubar, vbar, zeta,                             &
     &                    ad_ubar, ad_vbar, ad_zeta)
# ifdef SOLVE3D
      CALL ad_u3dbc_tile (ng, Istr, Iend, Jstr, Jend,                   &
     &                    LBi, UBi, LBj, UBj, N(ng),                    &
     &                    nstp, nstp,                                   &
     &                    ad_u)
      CALL ad_v3dbc_tile (ng, Istr, Iend, Jstr, Jend,                   &
     &                    LBi, UBi, LBj, UBj, N(ng),                    &
     &                    nstp, nstp,                                   &
     &                    ad_v)
      DO itrc=1,NT(ng)
        CALL ad_t3dbc_tile (ng, Istr, Iend, Jstr, Jend, itrc,           &
     &                      LBi, UBi, LBj, UBj, N(ng), NT(ng),          &
     &                      nstp, nstp,                                 &
     &                      ad_t)
      END DO
# endif
!
!  Add adjoint contribution from the computation of 2D momentun from
!  vertical integration of 3D momentum.  Notice that we need use the
!  sum of levels "kstp" and "knew" when computing the adjoint from
!  "tl_vbar" and "tl_ubar".
!
      DO j=JstrV,Jend
        DO i=Istr,Iend                                   ! BASIC STATE
          CF(i)=0.0_r8
          DC(i)=0.0_r8
        END DO
        DO k=1,N(ng)
          DO i=Istr,Iend
            cff=0.5_r8*(Hz(i,j-1,k)+Hz(i,j,k))
            CF(i)=CF(i)+cff
            DC(i)=DC(i)+v(i,j,k,nstp)*cff
          END DO
        END DO
        DO i=Istr,Iend
          cff=1.0_r8/CF(i)
# ifdef MASKING
!>        tl_vbar(i,j,kstp)=tl_vbar(i,j,kstp)*vmask(i,j)
!>
          ad_vbar(i,j,kstp)=ad_vbar(i,j,kstp)*vmask(i,j)
          ad_vbar(i,j,knew)=ad_vbar(i,j,knew)*vmask(i,j)
# endif
!>        tl_vbar(i,j,kstp)=tl_DC(i)*cff+DC(i)*tl_cff
!>                                                         ! use both
          adfac=ad_vbar(i,j,kstp)+ad_vbar(i,j,knew)        ! time levels
          ad_cff=ad_cff+DC(i)*adfac
          ad_DC(i)=ad_DC(i)+cff*adfac
          ad_vbar(i,j,kstp)=0.0_r8
          ad_vbar(i,j,knew)=0.0_r8
!>        tl_cff=-cff*cff*tl_CF(i)
!>
          ad_CF(i)=ad_CF(i)-cff*cff*ad_cff
          ad_cff=0.0_r8
        END DO
!>      DO k=1,N(ng)
!>
        DO k=N(ng),1,-1
          DO i=Istr,Iend
            cff=0.5_r8*(Hz(i,j-1,k)+Hz(i,j,k))
!>          tl_DC(i)=tl_DC(i)+tl_v(i,j,k,nstp)*cff+                     &
!>   &                        v(i,j,k,nstp)*tl_cff
!>
            ad_cff=ad_cff+v(i,j,k,nstp)*ad_DC(i)
            ad_v(i,j,k,nstp)=ad_v(i,j,k,nstp)+cff*ad_DC(i)
!>          tl_CF(i)=tl_CF(i)+tl_cff
!>
            ad_cff=ad_cff+ad_CF(i)
!>          tl_cff=0.5_r8*(tl_Hz(i,j-1,k)+tl_Hz(i,j,k))
!>
            adfac=0.5_r8*ad_cff
            ad_Hz(i,j-1,k)=ad_Hz(i,j-1,k)+ad_cff
            ad_Hz(i,j  ,k)=ad_Hz(i,j  ,k)+ad_cff
            ad_cff=0.0_r8
          END DO
        END DO
        DO i=Istr,Iend
          ad_DC(i)=0.0_r8
          ad_CF(i)=0.0_r8
        END DO
      END DO
!
      DO j=Jstr,Jend
        DO i=IstrU,Iend                                  ! BASIC STATE
          CF(i)=0.0_r8
          DC(i)=0.0_r8
        END DO
        DO k=1,N(ng)
          DO i=IstrU,Iend
            cff=0.5_r8*(Hz(i-1,j,k)+Hz(i,j,k))
            CF(i)=CF(i)+cff
            DC(i)=DC(i)+u(i,j,k,nstp)*cff
          END DO
        END DO
        DO i=IstrU,Iend
          cff=1.0_r8/CF(i)
# ifdef MASKING
!>        tl_ubar(i,j,kstp)=tl_ubar(i,j,kstp)*umask(i,j)
!>
          ad_ubar(i,j,kstp)=ad_ubar(i,j,kstp)*umask(i,j)
          ad_ubar(i,j,knew)=ad_ubar(i,j,kstp)*umask(i,j)
# endif
!>        tl_ubar(i,j,kstp)=tl_DC(i)*cff+DC(i)*tl_cff
!>                                                         ! use both
          adfac=ad_ubar(i,j,kstp)+ad_ubar(i,j,knew)        ! time levels
          ad_cff=ad_cff+DC(i)*adfac
          ad_DC(i)=ad_DC(i)+cff*adfac
          ad_ubar(i,j,kstp)=0.0_r8
          ad_vbar(i,j,knew)=0.0_r8
!>        tl_cff=-cff*cff*tl_CF(i)
!>
          ad_CF(i)=ad_CF(i)-cff*cff*ad_cff
          ad_cff=0.0_r8
        END DO
!>      DO k=1,N(ng)
!>
        DO k=N(ng),1,-1
          DO i=IstrU,Iend
            cff=0.5_r8*(Hz(i-1,j,k)+Hz(i,j,k))
!>          tl_DC(i)=tl_DC(i)+tl_u(i,j,k,nstp)*cff+                     &
!>   &                        u(i,j,k,nstp)*tl_cff
!>
            ad_cff=ad_cff+u(i,j,k,nstp)*ad_DC(i)
            ad_u(i,j,k,nstp)=ad_u(i,j,k,nstp)+cff*ad_DC(i)            
!>          tl_CF(i)=tl_CF(i)+tl_cff
!>
            ad_cff=ad_cff+ad_CF(i)
!>          tl_cff=0.5_r8*(tl_Hz(i-1,j,k)+tl_Hz(i,j,k))
!>
            adfac=0.5_r8*ad_cff
            ad_Hz(i-1,j,k)=ad_Hz(i-1,j,k)+adfac
            ad_Hz(i  ,j,k)=ad_Hz(i  ,j,k)+adfac
            ad_cff=0.0_r8
          END DO
        END DO
        DO i=IstrU,Iend
          ad_DC(i)=0.0_r8
          ad_CF(i)=0.0_r8
        END DO
      END DO
!
!  Collect contibutions from free-sruface boundary exchanges. Notice
!  that both times levels "kstp" and "knew" are needed.
!
# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      CALL ad_exchange_r2d_tile (ng, iADM, Istr, Iend, Jstr, Jend,      &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           NghostPoints,                          &
     &                           ad_zeta(:,:,kstp)
      CALL ad_exchange_r2d_tile (ng, iADM, Istr, Iend, Jstr, Jend,      &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           NghostPoints,                          &
     &                           ad_zeta(:,:,knew)
# endif
# ifdef SOLVE3D
!
!  Collect contibution from "ad_zeta_avg1".
!  HGA: What is the correct index to collect from?
!
      CALL ad_set_depth_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        kstp, FirstPass,                          &
     &                        ad_zeta,                                  &
     &                        h, ad_h,                                  &
#  ifdef ICESHELF
     &                        zice,                                     &
#  endif
#  if defined SEDIMENT && defined SED_MORPH
     &                        ad_bed, ad_bed_thick0,                    &
#  endif
     &                        Zt_avg1, ad_Zt_avg1,                      &
     &                        ad_Hz, ad_z_r, ad_z_w)
# endif
!
!
!  Collect contibutions from free-sruface lateral boundary conditions.
!  Notice that both times levels "kstp" and "knew" are needed.
!
      CALL ad_zetabc_tile (ng, Istr, Iend, Jstr, Jend,                  &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     krhs, kstp, kstp,                            &
     &                     zeta, ad_zeta)
      CALL ad_zetabc_tile (ng, Istr, Iend, Jstr, Jend,                  &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     krhs, kstp, knew,                            &
     &                     zeta, ad_zeta)

# ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Initialize adjoint state vector with special value to facilitate
!  collection between all nodes.
!-----------------------------------------------------------------------
!
      DO is=1,Nstate
        ad_state(is)=Aspv
      END DO
# endif
!
!-----------------------------------------------------------------------
!  Extract state variables from full 1D state vector.
!-----------------------------------------------------------------------
!
!  Determine the index offset for each variable in the state vector.
# ifdef MASKING
!  Notice that in Land/Sea masking application the state vector only
!  contains water points to avoid large null space.
# endif
!
# ifdef SOLVE3D
#  ifdef MASKING
      offset(isFsur)=0
      offset(isUvel)=offset(isFsur)+NwaterR(ng)
      offset(isVvel)=offset(isUvel)+NwaterU(ng)*N(ng)
      iadd=NwaterV(ng)*N(ng)
      DO itrc=1,NT(ng)
        offset(isTvar(itrc))=offset(isTvar(itrc)-1)+iadd
        iadd=NwaterR(ng)*N(ng)
      END DO
#  else
#   ifdef FULL_GRID
      offset(isFsur)=0
      offset(isUvel)=offset(isFsur)+(Lm(ng)+2)*(Mm(ng)+2)
      offset(isVvel)=offset(isUvel)+(Lm(ng)+1)*(Mm(ng)+2)*N(ng)
      iadd=(Lm(ng)+2)*(Mm(ng)+1)*N(ng)
      DO itrc=1,NT(ng)
        offset(isTvar(itrc))=offset(isTvar(itrc)-1)+iadd
        iadd=(Lm(ng)+2)*(Mm(ng)+2)*N(ng)
      END DO
#   else
      offset(isFsur)=0
      offset(isUvel)=offset(isFsur)+Lm(ng)*Mn(ng)
      offset(isVvel)=offset(isUvel)+(Lm(ng)-1)*Mm(ng)*N(ng)
      iadd=Lm(ng)*(Mm(ng)-1)*N(ng)
      DO itrc=1,NT(ng)
        offset(isTvar(itrc))=offset(isTvar(itrc)-1)+iadd
        iadd=Lm(ng)*Mm(ng)*N(ng)
      END DO
#   endif
#  endif
# else
#  ifdef MASKING
      offset(isFsur)=0
      offset(isUbar)=offset(isFsur)+NwaterR(ng)
      offset(isVbar)=offset(isUbar)+NwaterU(ng)
#  else
#   ifdef FULL_GRID
      offset(isFsur)=0
      offset(isUbar)=offset(isFsur)+(Lm(ng)+2)*(Mm(ng)+2)
      offset(isVbar)=offset(isUbar)+(Lm(ng)+1)*(Mm(ng)+2)
#   else
      offset(isFsur)=0
      offset(isUbar)=offset(isFsur)+Lm(ng)*Mm(ng)
      offset(isVbar)=offset(isUbar)+(Lm(ng)-1)*Mm(ng)
#   endif
#  endif
# endif
!
!  Pack adjoint free-surface.
!
# ifndef MASKING
#  ifdef FULL_GRID
      Imax=Lm(ng)+2
      Ioff=1
      Joff=0
#  else
      Imax=Lm(ng)
      Ioff=0
      Joff=1
#  endif
# endif
      scale=1.0_r8/SQRT(0.5_r8*g*rho0)
      DO j=JR_RANGE
        DO i=IR_RANGE
# ifdef MASKING
          IF (rmask(i,j).gt.0.0_r8) THEN
            is=IJwaterR(i,j)+offset(isFsur)
            ad_state(is)=scale*(ad_zeta(i,j,kstp)+ad_zeta(i,j,knew))
          END IF
# else
          is=(i=Ioff)+(j-Joff)*Imax+offset(isFsur)
          ad_state(is)=scale*(ad_zeta(i,j,kstp)+ad_zeta(i,j,knew))
# endif
        END DO
      END DO
# ifndef SOLVE3D
!
!  Pack adjoint 2D U-velocity.
!
#  ifndef MASKING
#   ifdef FULL_GRID
      Imax=Lm(ng)+1
      Ioff=0
      Joff=0
#   else
      Imax=Lm(ng)-1
      Ioff=1
      Joff=1
#   endif
#  endif
      cff=0.25_r8*rho0
      DO j=JR_RANGE
        DO i=IU_RANGE
          scale=1.0_r8/SQRT(cff*(h(i-1,j)+h(i,j)))
#  ifdef MASKING
          IF (umask(i,j).gt.0.0_r8) THEN
            is=IJwaterU(i,j)+offset(isUbar)
            ad_state(is)=scale*(ad_ubar(i,j,kstp)+ad_ubar(i,j,knew))
          END IF
#  else
          is=(i-Ioff)+(j-Joff)*Imax+offset(isUbar)
          ad_state(is)=scale*(ad_ubar(i,j,kstp)+ad_ubar(i,j,knew))
#  endif
        END DO
      END DO
!
!  Pack adjoint 2D V-velocity.
!
#  ifndef MASKING
#   ifdef FULL_GRID
      Imax=Lm(ng)+2
      Ioff=1
      Joff=1
#   else
      Imax=Lm(ng)
      Ioff=0
      Joff=2
#   endif
#  endif
      cff=0.25_r8*rho0
      DO j=JV_RANGE
        DO i=IR_RANGE
          scale=1.0_r8/SQRT(cff*(h(i,j-1)+h(i,j)))
#  ifdef MASKING
          IF (vmask(i,j).gt.0.0_r8) THEN
            is=IJwaterV(i,j)+offset(isVbar)
            ad_state(is)=scale*(ad_vbar(i,j,kstp)+ad_vbar(i,j,knew))
          END IF
#  else
          is=(i+Ioff)+(j-Joff)*Imax+offset(isVbar)
          ad_state(is)=scale*(ad_vbar(i,j,kstp)+ad_vbar(i,j,knew))
#  endif
        END DO
      END DO
# else
!
!  Pack adjoint 3D U-velocity.
!
#  ifndef MASKING
#   ifdef FULL_GRID
      Imax=Lm(ng)+1
      Jmax=Mm(ng)+2
      Ioff=0
      Joff=0
#   else
      Imax=Lm(ng)-1
      Jmax=Mm(ng)
      Ioff=1
      Joff=1
#   endif
#  endif
      cff=0.25_r8*rho0
      DO k=1,N(ng)
#  ifdef MASKING
        iadd=(k-1)*NwaterU(ng)+offset(isUvel)
#  else
        iadd=(k-1)*Imax*Jmax+offset(isUvel)
#  endif
        DO j=JR_RANGE
          DO i=IU_RANGE
#  ifdef MASKING
            IF (umask(i,j).gt.0.0_r8) THEN
              scale=1.0_r8/SQRT(cff*(Hz(i-1,j,k)+Hz(i,j,k)))
              is=IJwaterU(i,j)+iadd
              ad_state(is)=scale*ad_u(i,j,k,nstp)
            END IF
#  else
            scale=1.0_r8/SQRT(cff*(Hz(i-1,j,k)+Hz(i,j,k)))
            is=(i-Ioff)+(j-Joff)*Imax+iadd
            ad_state(is)=scale*ad_u(i,j,k,nstp)
#  endif
          END DO
        END DO
      END DO
!
!  Pack adjoint 3D V-velocity.
!
#  ifndef MASKING
#   ifdef FULL_GRID
      Imax=Lm(ng)+2
      Jmax=Mm(ng)+1
      Ioff=1
      Joff=1
#   else
      Imax=Lm(ng)
      Jmax=Mm(ng)-1
      Ioff=0
      Joff=2
#   endif
#  endif
      cff=0.25_r8*rho0
      DO k=1,N(ng)
#  ifdef MASKING
        iadd=(k-1)*NwaterV(ng)+offset(isVvel)
#  else
        iadd=(k-1)*Imax*Jmax+offset(isVvel)
#  endif
        DO j=JV_RANGE
          DO i=IR_RANGE
#  ifdef MASKING
            IF (vmask(i,j).gt.0.0_r8) THEN
              scale=1.0_r8/SQRT(cff*(Hz(i,j-1,k)+Hz(i,j,k)))
              is=IJwaterV(i,j)+iadd
              ad_state(is)=scale*ad_v(i,j,k,nstp)
            END IF
#  else
            scale=1.0_r8/SQRT(cff*(Hz(i,j-1,k)+Hz(i,j,k)))
            is=(i+Ioff)+(j-Joff)*Imax+iadd
            ad_state(is)=scale*ad_v(i,j,k,nstp)
#  endif
          END DO
        END DO
      END DO
!
!  Pack adjoint tracers variables. For now, use salinity scale for
!  passive tracers.
!
# ifndef MASKING
#  ifdef FULL_GRID
      Imax=Lm(ng)+2
      Jmax=Mm(ng)+2
      Ioff=1
      Joff=0
#  else
      Imax=Lm(ng)
      Imax=Mm(ng)
      Ioff=0
      Joff=1
#  endif
# endif
      DO itrc=1,NT(ng)
        IF (itrc.eq.itemp) THEN
          cff=0.5_r8*rho0*Tcoef(ng)*Tcoef(ng)*g*g/bvf_bak
        ELSE IF (itrc.eq.isalt) THEN
          cff=0.5_r8*rho0*Scoef(ng)*Scoef(ng)*g*g/bvf_bak
        ELSE
          cff=0.5_r8*rho0*Scoef(ng)*Scoef(ng)*g*g/bvf_bak
        END IF
        DO k=1,N(ng)
#  ifdef MASKING
          iadd=(k-1)*NwaterR(ng)+offset(isTvar(itrc))
#  else
          iadd=(k-1)*Imax*Jmax+offset(isTvar(itrc))
#  endif
          DO j=JR_RANGE
            DO i=IR_RANGE
#  ifdef MASKING
              IF (rmask(i,j).gt.0.0_r8) THEN
                scale=1.0_r8/SQRT(cff*Hz(i,j,k))
                is=IJwaterR(i,j)+iadd
                ad_state(is)=scale*ad_t(i,j,k,nstp,itrc)
              END IF
#  else
              scale=1.0_r8/SQRT(cff*Hz(i,j,k))
              is=(i+Ioff)+(j-Joff)*Imax+iadd
              ad_state(is)=scale*ad_t(i,j,k,nstp,itrc)
#  endif
            END DO
          END DO
        END DO
      END DO
# endif
# ifdef DISTRIBUTE
!
!------------------------------------------------------------------------
!   Collect all packed data.
!------------------------------------------------------------------------
!
      CALL mp_collect (ng, iADM, Nstate, Aspv, ad_state)
# endif
      RETURN
      END SUBROUTINE ad_pack_tile
#endif
#undef IR_RANGE
#undef IU_RANGE
#undef JR_RANGE
#undef JV_RANGE
      END MODULE packing_mod
     
