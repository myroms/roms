#include "cppdefs.h"
      MODULE rhs3d_mod
#ifdef SOLVE3D
!
!========================================== Alexander F. Shchepetkin ===
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This subroutine evaluates the right-hand-side term of the 3D        !
!  momentum and tracers equations.                                     !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: rhs3d

      CONTAINS
!
!***********************************************************************
      SUBROUTINE rhs3d (ng, tile)
!***********************************************************************
!
      USE mod_param
# ifdef CLIMATOLOGY
      USE mod_clima
# endif
      USE mod_coupling
# if defined DIAGNOSTICS_TS || defined DIAGNOSTICS_UV
      USE mod_diags
# endif
      USE mod_forces
      USE mod_grid
# if defined ASSIMILATION || defined NUDGING
      USE mod_obs
# endif
      USE mod_ocean
      USE mod_stepping
!
      USE pre_step3d_mod, ONLY : pre_step3d
      USE prsgrd_mod, ONLY : prsgrd
# ifndef TS_FIXED
#  ifdef TS_DIF2
#   if defined MIX_S_TS
      USE t3dmix2_s_mod, ONLY : t3dmix2_s
#   elif defined MIX_GEO_TS
      USE t3dmix2_geo_mod, ONLY : t3dmix2_geo
#   elif defined MIX_ISO_TS
      USE t3dmix2_iso_mod, ONLY : t3dmix2_iso
#   endif
#  endif
#  ifdef TS_DIF4
#   if defined MIX_S_TS
      USE t3dmix4_s_mod, ONLY : t3dmix4_s
#   elif defined MIX_GEO_TS
      USE t3dmix4_geo_mod, ONLY : t3dmix4_geo
#   elif defined MIX_ISO_TS
      USE t3dmix4_iso_mod, ONLY : t3dmix4_iso
#   endif
#  endif
# endif
# ifdef UV_VIS2
#  if defined MIX_S_UV
      USE uv3dmix2_s_mod, ONLY : uv3dmix2_s
#  elif defined MIX_GEO_UV
      USE uv3dmix2_geo_mod, ONLY: uv3dmix2_geo
#  endif
# endif
# ifdef UV_VIS4
#  if defined MIX_S_UV
      USE uv3dmix4_s_mod, ONLY : uv3dmix4_s
#  elif defined MIX_GEO_UV
      USE uv3dmix4_geo_mod, ONLY : uv3dmix4_geo
#  endif
# endif
!
      integer, intent(in) :: ng, tile

# include "tile.h"
!
!-----------------------------------------------------------------------
!  Initialize computations for new time step of the 3D primitive
!  variables.
!-----------------------------------------------------------------------
!
      CALL pre_step3d (ng, tile)
!
!-----------------------------------------------------------------------
!  Compute baroclinic pressure gradient.
!-----------------------------------------------------------------------
!
      CALL prsgrd (ng, tile)
# ifndef TS_FIXED
#  ifdef TS_DIF2
!
!-----------------------------------------------------------------------
!  Compute horizontal harmonic mixing of tracer type variables.
!-----------------------------------------------------------------------
!
#   if defined MIX_S_TS
      CALL t3dmix2_s (ng, tile)
#   elif defined MIX_GEO_TS
      CALL t3dmix2_geo (ng, tile)
#   elif defined MIX_ISO_TS
      CALL t3dmix2_iso (ng, tile)
#   endif
#  endif
#  ifdef TS_DIF4
!
!-----------------------------------------------------------------------
!  Compute horizontal biharmonic mixing of tracer type variables.
!-----------------------------------------------------------------------
!
#   if defined MIX_S_TS
      CALL t3dmix4_s (ng, tile)
#   elif defined MIX_GEO_TS
      CALL t3dmix4_geo (ng, tile)
#   elif defined MIX_ISO_TS
      CALL t3dmix4_iso (ng, tile)
#   endif
#  endif
# endif
!
!-----------------------------------------------------------------------
!  Compute the right-hand-side term of the 3D momentum equations.
!-----------------------------------------------------------------------
!
# ifdef PROFILE
      CALL wclock_on (ng, 21)
# endif
      CALL rhs3d_tile (ng, Istr, Iend, Jstr, Jend,                      &
     &                 LBi, UBi, LBj, UBj,                              &
     &                 nrhs(ng),                                        &
     &                 GRID(ng) % Hz,                                   &
     &                 GRID(ng) % Huon,                                 &
     &                 GRID(ng) % Hvom,                                 &
# if defined CURVGRID && defined UV_ADV
     &                 GRID(ng) % dmde,                                 &
     &                 GRID(ng) % dndx,                                 &
# endif
     &                 GRID(ng) % fomn,                                 &
     &                 GRID(ng) % om_u,                                 &
     &                 GRID(ng) % om_v,                                 &
     &                 GRID(ng) % on_u,                                 &
     &                 GRID(ng) % on_v,                                 &
     &                 GRID(ng) % pm,                                   &
     &                 GRID(ng) % pn,                                   &
     &                 FORCES(ng) % bustr,                              &
     &                 FORCES(ng) % bvstr,                              &
     &                 FORCES(ng) % sustr,                              &
     &                 FORCES(ng) % svstr,                              &
# ifdef M3CLM_NUDGING
     &                 CLIMA(ng) % M3nudgcof,                           &
     &                 CLIMA(ng) % uclm,                                &
     &                 CLIMA(ng) % vclm,                                &
# endif
# if defined NUDGING_UVsur || defined NUDGING_UV
     &                 OBS(ng) % EobsUV,                                &
     &                 OBS(ng) % Uobs,                                  &
     &                 OBS(ng) % Vobs,                                  &
# endif
     &                 OCEAN(ng) % u,                                   &
     &                 OCEAN(ng) % v,                                   &
     &                 OCEAN(ng) % W,                                   &
     &                 COUPLING(ng) % Zt_avg1,                          &
     &                 COUPLING(ng) % rufrc,                            &
     &                 COUPLING(ng) % rvfrc,                            &
# ifdef DIAGNOSTICS_UV
     &                 DIAGS(ng) % DiaRUfrc,                            &
     &                 DIAGS(ng) % DiaRVfrc,                            &
     &                 DIAGS(ng) % DiaRU,                               &
     &                 DIAGS(ng) % DiaRV,                               &
# endif
     &                 OCEAN(ng) % ru,                                  &
     &                 OCEAN(ng) % rv,                                  &
     &                 OCEAN(ng) % zeta)
# ifdef PROFILE
      CALL wclock_off (ng, 21)
# endif
# ifdef UV_VIS2
!
!-----------------------------------------------------------------------
!  Compute horizontal, harmonic mixing of momentum.
!-----------------------------------------------------------------------
!
#  if defined MIX_S_UV
      CALL uv3dmix2_s (ng, tile)
#  elif defined MIX_GEO_UV
      CALL uv3dmix2_geo (ng, tile)
#  endif
# endif
# ifdef UV_VIS4
!
!-----------------------------------------------------------------------
!  Compute horizontal, biharmonic mixing of momentum.
!-----------------------------------------------------------------------
!
#  if defined MIX_S_UV
      CALL uv3dmix4_s (ng, tile)
#  elif defined MIX_GEO_UV
      CALL uv3dmix4_geo (ng, tile)
#  endif
# endif
      RETURN
      END SUBROUTINE rhs3d
!
!***********************************************************************
      SUBROUTINE rhs3d_tile (ng, Istr, Iend, Jstr, Jend,                &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       nrhs,                                      &
     &                       Hz, Huon, Hvom,                            &
# if defined CURVGRID && defined UV_ADV
     &                       dmde, dndx,                                &
# endif
     &                       fomn, om_u, om_v, on_u, on_v, pm, pn,      &
     &                       bustr, bvstr, sustr, svstr,                &
# ifdef M3CLM_NUDGING
     &                       M3nudgcof, uclm, vclm,                     &
# endif
# if defined NUDGING_UVsur || defined NUDGING_UV
     &                       EobsUV, Uobs, Vobs,                        &
# endif
     &                       u, v, W, Zt_avg1,                          &
     &                       rufrc, rvfrc,                              &
# ifdef DIAGNOSTICS_UV
     &                       DiaRUfrc, DiaRVfrc,                        &
     &                       DiaRU, DiaRV,                              &
# endif
     &                       ru, rv, zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars

# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
!
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: nrhs

# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: Huon(LBi:,LBj:,:)
      real(r8), intent(in) :: Hvom(LBi:,LBj:,:)
#  if defined CURVGRID && defined UV_ADV
      real(r8), intent(in) :: dmde(LBi:,LBj:)
      real(r8), intent(in) :: dndx(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: fomn(LBi:,LBj:)
      real(r8), intent(in) :: om_u(LBi:,LBj:)
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)
      real(r8), intent(in) :: on_v(LBi:,LBj:)
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: bustr(LBi:,LBj:)
      real(r8), intent(in) :: bvstr(LBi:,LBj:)
      real(r8), intent(in) :: sustr(LBi:,LBj:)
      real(r8), intent(in) :: svstr(LBi:,LBj:)
#  ifdef M3CLM_NUDGING
      real(r8), intent(in) :: M3nudgcof(LBi:,LBj:)
      real(r8), intent(in) :: uclm(LBi:,LBj:,:)
      real(r8), intent(in) :: vclm(LBi:,LBj:,:)
#  endif
#  if defined NUDGING_UVsur || defined NUDGING_UV
      real(r8), intent(in) :: EobsUV(LBi:,LBj:,:)
      real(r8), intent(in) :: Uobs(LBi:,LBj:,:)
      real(r8), intent(in) :: Vobs(LBi:,LBj:,:)
#  endif
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
      real(r8), intent(in) :: W(LBi:,LBj:,0:)
      real(r8), intent(in) :: Zt_avg1(LBi:,LBj:)

      real(r8), intent(inout) :: ru(LBi:,LBj:,0:,:)
      real(r8), intent(inout) :: rv(LBi:,LBj:,0:,:)
#  ifdef DIAGNOSTICS_UV
      real(r8), intent(inout) :: DiaRUfrc(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: DiaRVfrc(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: DiaRU(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: DiaRV(LBi:,LBj:,:,:,:)
#  endif

      real(r8), intent(out) :: rufrc(LBi:,LBj:)
      real(r8), intent(out) :: rvfrc(LBi:,LBj:)
      real(r8), intent(out) :: zeta(LBi:,LBj:,:)
# else
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: Huon(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: Hvom(LBi:UBi,LBj:UBj,N(ng))
#  if defined CURVGRID && defined UV_ADV
      real(r8), intent(in) :: dmde(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: dndx(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: fomn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: bustr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: bvstr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: svstr(LBi:UBi,LBj:UBj)
#  ifdef M3CLM_NUDGING
      real(r8), intent(in) :: M3nudgcof(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: uclm(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: vclm(LBi:UBi,LBj:UBj,N(ng))
#  endif
#  if defined NUDGING_UVsur || defined NUDGING_UV
      real(r8), intent(in) :: EobsUV(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: Uobs(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: Vobs(LBi:UBi,LBj:UBj,N(ng))
#  endif
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: W(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(in) :: Zt_avg1(LBi:UBi,LBj:UBj)

      real(r8), intent(inout) :: ru(LBi:UBi,LBj:UBj,0:N(ng),2)
      real(r8), intent(inout) :: rv(LBi:UBi,LBj:UBj,0:N(ng),2)
#  ifdef DIAGNOSTICS_UV
      real(r8), intent(inout) :: DiaRUfrc(LBi:UBi,LBj:UBj,3,NDM2d-1)
      real(r8), intent(inout) :: DiaRVfrc(LBi:UBi,LBj:UBj,3,NDM2d-1)
      real(r8), intent(inout) :: DiaRU(LBi:UBi,LBj:UBj,N(ng),2,NDrhs)
      real(r8), intent(inout) :: DiaRV(LBi:UBi,LBj:UBj,N(ng),2,NDrhs)
#  endif

      real(r8), intent(out) :: rufrc(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: rvfrc(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: zeta(LBi:UBi,LBj:UBj,3)
# endif
!
!  Local variable declarations. Because the adjoint model the
!  equivalence via CPP is removed.  Better usage of private
!  arrays will be done latter.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j, k

      real(r8), parameter :: Gadv = -0.25_r8

      real(r8) :: cff, cff1, cff2, cff3, cff4
      real(r8) :: fac, fac1, fac2

      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N(ng)) :: CF
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N(ng)) :: DC
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N(ng)) :: FC

      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: Huee
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: Huxx
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: Hvee
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: Hvxx
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: UFx
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: UFe
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: Uwrk
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: VFx
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: VFe
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: Vwrk
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: uee
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: uxx
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: vee
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: vxx
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: wrk

# include "set_bounds.h"

# ifdef BODYFORCE
!
!-----------------------------------------------------------------------
!  Apply surface stress as a bodyforce: determine the thickness (m)
!  of the surface layer; then add in surface stress as a bodyfoce.
!-----------------------------------------------------------------------
!
#  ifdef DIAGNOSTICS_UV
      DO k=1,N(ng)
        DO j=JstrV-1,Jend
          DO i=IstrU-1,Iend
            DiaRU(i,j,k,M3vvis)=0.0_r8
            DiaRV(i,j,k,M3vvis)=0.0_r8
          END DO
        END DO
      END DO
#  endif
      DO j=JstrV-1,Jend
        DO i=IstrU-1,Iend
          wrk(i,j)=0.0_r8
        END DO
      END DO
      DO k=N(ng),levsfrc(ng),-1
        DO j=JstrV-1,Jend
          DO i=IstrU-1,Iend
            wrk(i,j)=wrk(i,j)+Hz(i,j,k)
          END DO
        END DO
      END DO
      DO j=Jstr,Jend
        DO i=IstrU,Iend
          Uwrk(i,j)=sustr(i,j)*4.0_r8/((wrk(i-1,j)+wrk(i,j))*           &
     &                                 (pm (i-1,j)+pm (i,j))*           &
     &                                 (pn (i-1,j)+pn (i,j)))
        END DO
      END DO
      DO j=JstrV,Jend
        DO i=Istr,Iend
          Vwrk(i,j)=svstr(i,j)*4.0_r8/((wrk(i,j-1)+wrk(i,j))*           &
     &                                 (pm (i,j-1)+pm (i,j))*           &
     &                                 (pn (i,j-1)+pn (i,j)))
        END DO
      END DO
      DO k=levsfrc(ng),N(ng)
        DO j=Jstr,Jend
          DO i=IstrU,Iend
            fac=Uwrk(i,j)*(Hz(i,j,k)+Hz(i-1,j,k))
            ru(i,j,k,nrhs)=ru(i,j,k,nrhs)+fac
#  ifdef DIAGNOSTICS_UV
            DiaRU(i,j,k,nrhs,M3vvis)=DiaRU(i,j,k,nrhs,M3vvis)+fac
#  endif
          END DO
        END DO
        DO j=JstrV,Jend
          DO i=Istr,Iend
            fac=Vwrk(i,j)*(Hz(i,j,k)+Hz(i,j-1,k))
            rv(i,j,k,nrhs)=rv(i,j,k,nrhs)+fac
#  ifdef DIAGNOSTICS_UV
            DiaRV(i,j,k,nrhs,M3vvis)=DiaRV(i,j,k,nrhs,M3vvis)+fac
#  endif
          END DO
        END DO
      END DO
!
!  Apply bottom stress as a bodyforce: determine the thickness (m)
!  of the bottom layer; then add in bottom stress as a bodyfoce.
!
      DO j=JstrV-1,Jend
        DO i=IstrU-1,Iend
          wrk(i,j)=0.0_r8
        END DO
      END DO
      DO k=1,levbfrc(ng)
        DO j=JstrV-1,Jend
          DO i=IstrU-1,Iend
            wrk(i,j)=wrk(i,j)+Hz(i,j,k)
          END DO
        END DO
      END DO
      DO j=Jstr,Jend
        DO i=IstrU,Iend
          Uwrk(i,j)=bustr(i,j)*4.0_r8/((wrk(i-1,j)+wrk(i,j))*           &
     &                                 (pm (i-1,j)+pm (i,j))*           &
     &                                 (pn (i-1,j)+pn (i,j)))
        END DO
      END DO
      DO j=JstrV,Jend
        DO i=Istr,Iend
          Vwrk(i,j)=bvstr(i,j)*4.0_r8/((wrk(i,j-1)+wrk(i,j))*           &
     &                                 (pm (i,j-1)+pm (i,j))*           &
     &                                 (pn (i,j-1)+pn (i,j)))
        END DO
      END DO
      DO k=1,levbfrc(ng)
        DO j=Jstr,Jend
          DO i=IstrU,Iend
            fac=Uwrk(i,j)*(Hz(i,j,k)+Hz(i-1,j,k))
            ru(i,j,k,nrhs)=ru(i,j,k,nrhs)-fac
#  ifdef DIAGNOSTICS_UV
            DiaRU(i,j,k,nrhs,M3vvis)=DiaRU(i,j,k,nrhs,M3vvis)-fac
#  endif
          END DO
        END DO
        DO j=JstrV,Jend
          DO i=Istr,Iend
            fac=Vwrk(i,j)*(Hz(i,j,k)+Hz(i,j-1,k))
            rv(i,j,k,nrhs)=rv(i,j,k,nrhs)-fac
#  ifdef DIAGNOSTICS_UV
            DiaRV(i,j,k,nrhs,M3vvis)=DiaRV(i,j,k,nrhs,M3vvis)-fac
#  endif
          END DO
        END DO
      END DO
# endif
!
      DO k=1,N(ng)
!
# if defined UV_COR || (defined CURVGRID && defined UV_ADV)
!-----------------------------------------------------------------------
!  Add in Coriolis and curvilinear transformation terms, if any.
!-----------------------------------------------------------------------
!
        DO j=JstrV-1,Jend
          DO i=IstrU-1,Iend
            cff=0.5_r8*Hz(i,j,k)*(                                      &
#  ifdef UV_COR
     &          fomn(i,j)                                               &
#  endif
#  if (defined CURVGRID && defined UV_ADV)
     &         +0.5_r8*((v(i,j,k,nrhs)+v(i,j+1,k,nrhs))*dndx(i,j)-      &
     &                  (u(i,j,k,nrhs)+u(i+1,j,k,nrhs))*dmde(i,j))      &
#  endif
     &          )
            UFx(i,j)=cff*(v(i,j,k,nrhs)+v(i,j+1,k,nrhs))
            VFe(i,j)=cff*(u(i,j,k,nrhs)+u(i+1,j,k,nrhs))
#  ifdef DIAGNOSTICS_UV
#   ifdef UV_COR
            Uwrk(i,j)=0.5_r8*Hz(i,j,k)*fomn(i,j)*                       &
     &                (v(i,j,k,nrhs)+v(i,j+1,k,nrhs))
            Vwrk(i,j)=0.5_r8*Hz(i,j,k)*fomn(i,j)*                       &
     &                (u(i,j,k,nrhs)+u(i+1,j,k,nrhs))
#   else
            Uwrk(i,j)=0.0_r8
            Vwrk(i,j)=0.0_r8
#   endif
#  endif
          END DO
        END DO
        DO j=Jstr,Jend
          DO i=IstrU,Iend
            fac1=0.5_r8*(UFx(i,j)+UFx(i-1,j))
            ru(i,j,k,nrhs)=ru(i,j,k,nrhs)+fac1
#  ifdef DIAGNOSTICS_UV
            fac2=0.5_r8*(Uwrk(i,j)+Uwrk(i-1,j))
#   ifdef UV_COR
            DiaRU(i,j,k,nrhs,M3fcor)=fac2
#   endif
#   if defined UV_ADV && defined CURVGRID
            DiaRU(i,j,k,nrhs,M3hadv)=fac1-fac2
#   endif
#  endif
          END DO
        END DO
        DO j=JstrV,Jend
          DO i=Istr,Iend
            fac1=0.5_r8*(VFe(i,j)+VFe(i,j-1))
            rv(i,j,k,nrhs)=rv(i,j,k,nrhs)-fac1
#  ifdef DIAGNOSTICS_UV
            fac2=0.5_r8*(Vwrk(i,j)+Vwrk(i,j-1))
#   ifdef UV_COR
            DiaRV(i,j,k,nrhs,M3fcor)=-fac2
#   endif
#   if defined UV_ADV && defined CURVGRID
            DiaRV(i,j,k,nrhs,M3hadv)=-fac1+fac2
#   endif
#  endif
          END DO
        END DO
# endif
# ifdef M3CLM_NUDGING
!
!-----------------------------------------------------------------------
!  Add in nudging of 3D momentum climatology.
!-----------------------------------------------------------------------
!
        DO j=Jstr,Jend
          DO i=IstrU,Iend
            ru(i,j,k,nrhs)=ru(i,j,k,nrhs)+                              &
     &                     0.25_r8*(M3nudgcof(i-1,j)+M3nudgcof(i,j))*   &
     &                     (Hz(i-1,j,k)+Hz(i,j,k))*                     &
     &                     (uclm(i,j,k)-u(i,j,k,nrhs))*                 &
     &                     om_u(i,j)*on_u(i,j)
          END DO
        END DO
        DO j=JstrV,Jend
          DO i=Istr,Iend
            rv(i,j,k,nrhs)=rv(i,j,k,nrhs)+                              &
     &                     0.25_r8*(M3nudgcof(i,j-1)+M3nudgcof(i,j))*   &
     &                     (Hz(i,j-1,k)+Hz(i,j,k))*                     &
     &                     (vclm(i,j,k)-v(i,j,k,nrhs))*                 &
     &                     om_v(i,j)*on_v(i,j)
         END DO
       END DO
# endif
# if defined NUDGING_UVsur || defined NUDGING_UV
!
!-----------------------------------------------------------------------
!  Assimilate horizontal currents observations via nudging.
!-----------------------------------------------------------------------
!
!  Notice that the nondimensional error variance "EobsUV" is flipped,
!  so a value of unity represents unavailable data.
!
        IF (update_UV(ng)) THEN
          DO j=Jstr,Jend
            DO i=IstrU,Iend
              cff=MIN(1.0_r8,MAX(0.0_r8,0.5_r8*(EobsUV(i-1,j,k)+        &
     &                                          EobsUV(i  ,j,k))))
              cff=M3nudass(ng)*(1.0_r8-cff)
              ru(i,j,k,nrhs)=ru(i,j,k,nrhs)+                            &
     &                       cff*0.5_r8*(Hz(i,j,k)+Hz(i-1,j,k))*        &
     &                           (Uobs(i,j,k)-u(i,j,k,nrhs))*           &
     &                           om_u(i,j)*on_u(i,j)
            END DO
          END DO
          DO j=JstrV,Jend
            DO i=Istr,Iend
              cff=MIN(1.0_r8,MAX(0.0_r8,0.5_r8*(EobsUV(i,j-1,k)+        &
     &                                          EobsUV(i,j  ,k))))
              cff=M3nudass(ng)*(1.0_r8-cff)
              rv(i,j,k,nrhs)=rv(i,j,k,nrhs)+                            &
     &                       cff*0.5_r8*(Hz(i,j-1,k)+Hz(i,j,k))*        &
     &                           (Vobs(i,j,k)-v(i,j,k,nrhs))*           &
     &                           om_v(i,j)*on_v(i,j)
            END DO
          END DO
        END IF
# endif
# ifdef UV_ADV
!-----------------------------------------------------------------------
!  Add in horizontal advection of momentum.
!-----------------------------------------------------------------------
!
!  Compute diagonal [UFx,VFe] and off-diagonal [UFe,VFx] components
!  of tensor of momentum flux due to horizontal advection.
!
#  ifdef UV_C2ADVECTION
!
!  Second-order, centered differences advection.
!
        DO j=Jstr,Jend
          DO i=IstrU-1,Iend
            UFx(i,j)=0.25_r8*(u(i,j,k,nrhs)+u(i+1,j,k,nrhs))*           &
     &                       (Huon(i,j,k)+Huon(i+1,j,k))
          END DO
        END DO
        DO j=Jstr,Jend+1
          DO i=IstrU,Iend
            UFe(i,j)=0.25_r8*(u(i,j-1,k,nrhs)+u(i,j,k,nrhs))*           &
     &                       (Hvom(i-1,j,k)+Hvom(i,j,k))
          END DO
        END DO
        DO j=JstrV,Jend
          DO i=Istr,Iend+1
            VFx(i,j)=0.25_r8*(v(i-1,j,k,nrhs)+v(i,j,k,nrhs))*           &
     &                       (Huon(i,j-1,k)+Huon(i,j,k))
          END DO
        END DO
        DO j=JstrV-1,Jend
          DO i=Istr,Iend
            VFe(i,j)=0.25_r8*(v(i,j,k,nrhs)+v(i,j+1,k,nrhs))*           &
     &                       (Hvom(i,j,k)+Hvom(i,j+1,k))
          END DO
        END DO
#  else
#   ifdef EW_PERIODIC
#    define IU_RANGE IstrU-1,Iend+1
#   else
#    define IU_RANGE MAX(IstrU-1,2),MIN(Iend+1,Lm(ng))
#   endif
        DO j=Jstr,Jend
          DO i=IU_RANGE
            uxx (i,j)=u(i-1,j,k,nrhs)-2.0_r8*u(i,j,k,nrhs)+             &
     &                u(i+1,j,k,nrhs)
            Huxx(i,j)=Huon(i-1,j,k)-2.0_r8*Huon(i,j,k)+Huon(i+1,j,k)
          END DO
        END DO
#   undef IU_RANGE
#   ifndef EW_PERIODIC
        IF (WESTERN_EDGE) THEN
          DO j=Jstr,Jend
            uxx (Istr,j)=uxx (Istr+1,j)
            Huxx(Istr,j)=Huxx(Istr+1,j)
          END DO
        END IF
        IF (EASTERN_EDGE) THEN
          DO j=Jstr,Jend
            uxx (Iend+1,j)=uxx (Iend,j)
            Huxx(Iend+1,j)=Huxx(Iend,j)
          END DO
        END IF
#   endif
#   ifdef UV_C4ADVECTION
!
!  Fourth-order, centered differences advection.
!
        cff=1.0_r8/6.0_r8
        DO j=Jstr,Jend
          DO i=IstrU-1,Iend
            UFx(i,j)=0.25_r8*(u(i,j,k,nrhs)+u(i+1,j,k,nrhs)-            &
     &                        cff*(uxx (i,j)+uxx (i+1,j)))*             &
     &                       (Huon(i,j,k)+Huon(i+1,j,k)-                &
     &                        cff*(Huxx(i,j)+Huxx(i+1,j)))
          END DO
        END DO
#   else
!
!  Third-order, upstream bias advection with velocity dependent
!  hyperdiffusion.
!
        DO j=Jstr,Jend
          DO i=IstrU-1,Iend
            cff1=u(i,j,k,nrhs)+u(i+1,j,k,nrhs)
            IF (cff1.gt.0.0_r8) THEN
              cff=uxx(i,j)
            ELSE
              cff=uxx(i+1,j)
            END IF
            UFx(i,j)=0.25_r8*(cff1+Gadv*cff)*                           &
     &               (Huon(i,j,k)+Huon(i+1,j,k)+                        &
     &                Gadv*0.5_r8*(Huxx(i,j)+Huxx(i+1,j)))
          END DO
        END DO
#   endif
#   ifdef NS_PERIODIC
#    define JU_RANGE Jstr-1,Jend+1
#   else
#    define JU_RANGE MAX(Jstr-1,1),MIN(Jend+1,Mm(ng))
#   endif
        DO j=JU_RANGE
          DO i=IstrU,Iend
            uee(i,j)=u(i,j-1,k,nrhs)-2.0_r8*u(i,j,k,nrhs)+              &
     &               u(i,j+1,k,nrhs)
          END DO
        END DO
#   undef JU_RANGE
#   ifndef NS_PERIODIC
        IF (SOUTHERN_EDGE) THEN
          DO i=IstrU,Iend
            uee(i,Jstr-1)=uee(i,Jstr)
          END DO
        END IF
        IF (NORTHERN_EDGE) THEN
          DO i=IstrU,Iend
            uee(i,Jend+1)=uee(i,Jend)
          END DO
        END IF
#   endif
        DO j=Jstr,Jend+1
          DO i=IstrU-1,Iend
           Hvxx(i,j)=Hvom(i-1,j,k)-2.0_r8*Hvom(i,j,k)+Hvom(i+1,j,k)
          END DO
        END DO
#   ifdef UV_C4ADVECTION
        cff=1.0_r8/6.0_r8
        DO j=Jstr,Jend+1
          DO i=IstrU,Iend
            UFe(i,j)=0.25_r8*(u(i,j,k,nrhs)+u(i,j-1,k,nrhs)-            &
     &                        cff*(uee (i,j)+uee (i,j-1)))*             &
     &                       (Hvom(i,j,k)+Hvom(i-1,j,k)-                &
     &                        cff*(Hvxx(i,j)+Hvxx(i-1,j)))
          END DO
        END DO
#   else
        DO j=Jstr,Jend+1
          DO i=IstrU,Iend
            cff1=u(i,j,k,nrhs)+u(i,j-1,k,nrhs)
            cff2=Hvom(i,j,k)+Hvom(i-1,j,k)
            IF (cff2.gt.0.0_r8) THEN
              cff=uee(i,j-1)
            ELSE
              cff=uee(i,j)
            END IF
            UFe(i,j)=0.25_r8*(cff1+Gadv*cff)*                           &
     &               (cff2+Gadv*0.5_r8*(Hvxx(i,j)+Hvxx(i-1,j)))
          END DO
        END DO
#   endif
#   ifdef EW_PERIODIC
#    define IV_RANGE Istr-1,Iend+1
#   else
#    define IV_RANGE MAX(Istr-1,1),MIN(Iend+1,Lm(ng))
#   endif
        DO j=JstrV,Jend
          DO i=IV_RANGE
            vxx(i,j)=v(i-1,j,k,nrhs)-2.0_r8*v(i,j,k,nrhs)+              &
     &               v(i+1,j,k,nrhs)
          END DO
        END DO
#   undef IV_RANGE
#   ifndef EW_PERIODIC
        IF (WESTERN_EDGE) THEN
          DO j=JstrV,Jend
            vxx(Istr-1,j)=vxx(Istr,j)
          END DO
        END IF
        IF (EASTERN_EDGE) THEN
          DO j=JstrV,Jend
            vxx(Iend+1,j)=vxx(Iend,j)
          END DO
        END IF
#   endif
        DO j=JstrV-1,Jend
          DO i=Istr,Iend+1
           Huee(i,j)=Huon(i,j-1,k)-2.0_r8*Huon(i,j,k)+Huon(i,j+1,k)
          END DO
        END DO
#   ifdef UV_C4ADVECTION
        cff=1.0_r8/6.0_r8
        DO j=JstrV,Jend
          DO i=Istr,Iend+1
            VFx(i,j)=0.25_r8*(v(i,j,k,nrhs)+v(i-1,j,k,nrhs)-            &
     &                        cff*(vxx (i,j)+vxx (i-1,j)))*             &
     &                       (Huon(i,j,k)+Huon(i,j-1,k)-                &
     &                        cff*(Huee(i,j)+Huee(i,j-1)))
          END DO
        END DO
#   else
        DO j=JstrV,Jend
          DO i=Istr,Iend+1
            cff1=v(i,j,k,nrhs)+v(i-1,j,k,nrhs)
            cff2=Huon(i,j,k)+Huon(i,j-1,k)
            IF (cff2.gt.0.0_r8) THEN
              cff=vxx(i-1,j)
            ELSE
              cff=vxx(i,j)
            END IF
            VFx(i,j)=0.25_r8*(cff1+Gadv*cff)*                           &
     &               (cff2+Gadv*0.5_r8*(Huee(i,j)+Huee(i,j-1)))
          END DO
        END DO
#   endif
#   ifdef NS_PERIODIC
#    define JV_RANGE JstrV-1,Jend+1
#   else
#    define JV_RANGE MAX(JstrV-1,2),MIN(Jend+1,Mm(ng))
#   endif
        DO j=JV_RANGE
          DO i=Istr,Iend
            vee(i,j)=v(i,j-1,k,nrhs)-2.0_r8*v(i,j,k,nrhs)+              &
     &               v(i,j+1,k,nrhs)
            Hvee(i,j)=Hvom(i,j-1,k)-2.0_r8*Hvom(i,j,k)+Hvom(i,j+1,k)
          END DO
        END DO
#   undef JV_RANGE
#   ifndef NS_PERIODIC
        IF (SOUTHERN_EDGE) THEN
          DO i=Istr,Iend
            vee (i,Jstr)=vee (i,Jstr+1)
            Hvee(i,Jstr)=Hvee(i,Jstr+1)
          END DO
        END IF
        IF (NORTHERN_EDGE) THEN
          DO i=Istr,Iend
            vee (i,Jend+1)=vee (i,Jend)
            Hvee(i,Jend+1)=Hvee(i,Jend)
          END DO
        END IF
#   endif
#   ifdef UV_C4ADVECTION
        cff=1.0_r8/6.0_r8
        DO j=JstrV-1,Jend
          DO i=Istr,Iend
            VFe(i,j)=0.25_r8*(v(i,j,k,nrhs)+v(i,j+1,k,nrhs)-            &
     &                        cff*(vee (i,j)+vee (i,j+1)))*             &
     &                       (Hvom(i,j,k)+Hvom(i,j+1,k)-                &
     &                        cff*(Hvee(i,j)+Hvee(i,j+1)))
          END DO
        END DO
#   else
        DO j=JstrV-1,Jend
          DO i=Istr,Iend
            cff1=v(i,j,k,nrhs)+v(i,j+1,k,nrhs)
            IF (cff1.gt.0.0_r8) THEN
              cff=vee(i,j)
            ELSE
              cff=vee(i,j+1)
            END IF
            VFe(i,j)=0.25_r8*(cff1+Gadv*cff)*                           &
     &               (Hvom(i,j,k)+Hvom(i,j+1,k)+                        &
     &                Gadv*0.5_r8*(Hvee(i,j)+Hvee(i,j+1)))
          END DO
        END DO
#   endif
#  endif
!
!  Add in horizontal advection.
!
        DO j=Jstr,Jend
          DO i=IstrU,Iend
            fac=(UFx(i,j)-UFx(i-1,j))+(UFe(i,j+1)-UFe(i,j))
            ru(i,j,k,nrhs)=ru(i,j,k,nrhs)-fac
#  ifdef DIAGNOSTICS_UV
#   ifdef CURVGRID
            DiaRU(i,j,k,nrhs,M3hadv)=DiaRU(i,j,k,nrhs,M3hadv)-fac
#   else
            DiaRU(i,j,k,nrhs,M3hadv)=-fac
#   endif
#  endif
          END DO
        END DO
        DO j=JstrV,Jend
          DO i=Istr,Iend
            fac=(VFx(i+1,j)-VFx(i,j))+(VFe(i,j)-VFe(i,j-1))
            rv(i,j,k,nrhs)=rv(i,j,k,nrhs)-fac
#  ifdef DIAGNOSTICS_UV
#   ifdef CURVGRID
            DiaRV(i,j,k,nrhs,M3hadv)=DiaRV(i,j,k,nrhs,M3hadv)-fac
#   else
            DiaRV(i,j,k,nrhs,M3hadv)=-fac
#   endif
#  endif
          END DO
        END DO
# endif
      END DO
      DO j=Jstr,Jend
# ifdef UV_ADV
!
!-----------------------------------------------------------------------
!  Add in vertical advection.
!-----------------------------------------------------------------------
!
#  ifdef UV_SADVECTION
!
!  Construct conservative parabolic splines for the vertical
!  derivatives "CF" of u-momentum.
!
        cff1=9.0_r8/16.0_r8
        cff2=1.0_r8/16.0_r8
        DO k=1,N(ng)
          DO i=IstrU,Iend
            DC(i,k)=(cff1*(Hz(i  ,j,k)+Hz(i-1,j,k))-                    &
     &               cff2*(Hz(i+1,j,k)+Hz(i-2,j,k)))
          END DO
        END DO
        DO i=IstrU,Iend
          FC(i,0)=0.0_r8
          CF(i,0)=0.0_r8
        END DO
        DO k=1,N(ng)-1
          DO i=IstrU,Iend
            cff=1.0_r8/(2.0_r8*DC(i,k+1)+DC(i,k)*(2.0_r8-FC(i,k-1)))
            FC(i,k)=cff*DC(i,k+1)
            CF(i,k)=cff*(6.0_r8*(u(i,j,k+1,nrhs)-u(i,j,k,nrhs))-        &
     &                   DC(i,k)*CF(i,k-1))
          END DO
        END DO
        DO i=IstrU,Iend
          CF(i,N(ng))=0.0_r8
        END DO
        DO k=N(ng)-1,1,-1
          DO i=IstrU,Iend
            CF(i,k)=CF(i,k)-FC(i,k)*CF(i,k+1)
          END DO
        END DO
!
! Compute spline-interpolated, vertical advective u-momentum flux.
!
        cff3=1.0_r8/3.0_r8
        cff4=1.0_r8/6.0_r8
        DO k=1,N(ng)-1
          DO i=IstrU,Iend
            FC(i,k)=(cff1*(W(i  ,j,k)+W(i-1,j,k))-                      &
     &               cff2*(W(i+1,j,k)+W(i-2,j,k)))*                     &
     &              (u(i,j,k,nrhs)+                                     &
     &               DC(i,k)*(cff3*CF(i,k)+cff4*CF(i,k-1)))
          END DO
        END DO
        DO i=IstrU,Iend
          FC(i,N(ng))=0.0_r8
          FC(i,0)=0.0_r8
        END DO
#  elif defined UV_C2ADVECTION
        DO k=1,N(ng)-1
          DO i=IstrU,Iend
            FC(i,k)=0.25_r8*(u(i,j,k,nrhs)+u(i,j,k+1,nrhs))*            &
     &              (W(i,j,k)+W(i-1,j,k))
          END DO
        END DO
        DO i=IstrU,Iend
          FC(i,0)=0.0_r8
          FC(i,N(ng))=0.0_r8
        END DO
#  elif defined UV_C4ADVECTION
        cff1=9.0_r8/32.0_r8
        cff2=1.0_r8/32.0_r8
        DO k=2,N(ng)-2
          DO i=IstrU,Iend
            FC(i,k)=(cff1*(u(i,j,k  ,nrhs)+u(i,j,k+1,nrhs))-            &
     &               cff2*(u(i,j,k-1,nrhs)+u(i,j,k+2,nrhs)))*           &
     &              (W(i,j,k)+W(i-1,j,k))
          END DO
        END DO
        DO i=IstrU,Iend
          FC(i,N(ng))=0.0_r8
          FC(i,N(ng)-1)=(cff1*(u(i,j,N(ng)-1,nrhs)+                     &
     &                         u(i,j,N(ng)  ,nrhs))-                    &
     &                   cff2*(u(i,j,N(ng)-2,nrhs)+                     &
     &                         u(i,j,N(ng)  ,nrhs)))*                   &
     &                  (W(i,j,N(ng)-1)+W(i-1,j,N(ng)-1))
          FC(i,1)=(cff1*(u(i,j,1,nrhs)+u(i,j,2,nrhs))-                  &
     &             cff2*(u(i,j,1,nrhs)+u(i,j,3,nrhs)))*                 &
     &            (W(i,j,1)+W(i-1,j,1))
          FC(i,0)=0.0_r8
        END DO
#  else
        cff1=9.0_r8/16.0_r8
        cff2=1.0_r8/16.0_r8
        DO k=2,N(ng)-2
          DO i=IstrU,Iend
            FC(i,k)=(cff1*(u(i,j,k  ,nrhs)+u(i,j,k+1,nrhs))-            &
     &               cff2*(u(i,j,k-1,nrhs)+u(i,j,k+2,nrhs)))*           &
     &              (cff1*(W(i  ,j,k)+W(i-1,j,k))-                      &
     &               cff2*(W(i+1,j,k)+W(i-2,j,k)))
          END DO
        END DO
        DO i=IstrU,Iend
          FC(i,N(ng))=0.0_r8
          FC(i,N(ng)-1)=(cff1*(u(i,j,N(ng)-1,nrhs)+                     &
     &                         u(i,j,N(ng)  ,nrhs))-                    &
     &                   cff2*(u(i,j,N(ng)-2,nrhs)+                     &
     &                         u(i,j,N(ng)  ,nrhs)))*                   &
     &                  (cff1*(W(i  ,j,N(ng)-1)+                        &
     &                         W(i-1,j,N(ng)-1))-                       &
     &                   cff2*(W(i+1,j,N(ng)-1)+                        &
     &                         W(i-2,j,N(ng)-1)))
          FC(i,1)=(cff1*(u(i,j,1,nrhs)+u(i,j,2,nrhs))-                  &
     &             cff2*(u(i,j,1,nrhs)+u(i,j,3,nrhs)))*                 &
     &            (cff1*(W(i  ,j,1)+W(i-1,j,1))-                        &
     &             cff2*(W(i+1,j,1)+W(i-2,j,1)))
          FC(i,0)=0.0_r8
        END DO
#  endif
        DO k=1,N(ng)
          DO i=IstrU,Iend
            fac=(FC(i,k)-FC(i,k-1))
            ru(i,j,k,nrhs)=ru(i,j,k,nrhs)-fac
#  ifdef DIAGNOSTICS_UV
            DiaRU(i,j,k,nrhs,M3vadv)=-fac
#  endif
          END DO
        END DO
        IF (j.ge.JstrV) THEN
#  ifdef UV_SADVECTION
!
!  Construct conservative parabolic splines for the vertical
!  derivatives "CF" of v-momentum.
!
          cff1=9.0_r8/16.0_r8
          cff2=1.0_r8/16.0_r8
          DO k=1,N(ng)
            DO i=Istr,Iend
              DC(i,k)=(cff1*(Hz(i,j  ,k)+Hz(i,j-1,k))-                  &
     &                 cff2*(Hz(i,j+1,k)+Hz(i,j-2,k)))
            END DO
          END DO
          DO i=Istr,Iend
            FC(i,0)=0.0_r8
            CF(i,0)=0.0_r8
          END DO
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              cff=1.0_r8/(2.0_r8*DC(i,k+1)+DC(i,k)*(2.0_r8-FC(i,k-1)))
              FC(i,k)=cff*DC(i,k+1)
              CF(i,k)=cff*(6.0_r8*(v(i,j,k+1,nrhs)-v(i,j,k,nrhs))-      &
     &                     DC(i,k)*CF(i,k-1))
            END DO
          END DO
          DO i=Istr,Iend
            CF(i,N(ng))=0.0_r8
          END DO
          DO k=N(ng)-1,1,-1
            DO i=Istr,Iend
              CF(i,k)=CF(i,k)-FC(i,k)*CF(i,k+1)
            END DO
          END DO
!
! Compute spline-interpolated, vertical advective v-momentum flux.
!
          cff3=1.0_r8/3.0_r8
          cff4=1.0_r8/6.0_r8
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              FC(i,k)=(cff1*(W(i,j  ,k)+W(i,j-1,k))-                    &
     &                 cff2*(W(i,j+1,k)+W(i,j-2,k)))*                   &
     &                (v(i,j,k,nrhs)+                                   &
     &                 DC(i,k)*(cff3*CF(i,k)+cff4*CF(i,k-1)))
            END DO
          END DO
          DO i=Istr,Iend
            FC(i,N(ng))=0.0_r8
            FC(i,0)=0.0_r8
          END DO
#  elif defined UV_C2ADVECTION
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              FC(i,k)=0.25_r8*(v(i,j,k,nrhs)+v(i,j,k+1,nrhs))*          &
     &                (W(i,j,k)+W(i,j-1,k))
            END DO
          END DO
          DO i=Istr,Iend
            FC(i,0)=0.0_r8
            FC(i,N(ng))=0.0_r8
          END DO
#  elif defined UV_C4ADVECTION
          cff1=9.0_r8/32.0_r8
          cff2=1.0_r8/32.0_r8
          DO k=2,N(ng)-2
            DO i=Istr,Iend
              FC(i,k)=(cff1*(v(i,j,k  ,nrhs)+v(i,j,k+1,nrhs))-          &
     &                 cff2*(v(i,j,k-1,nrhs)+v(i,j,k+2,nrhs)))*         &
     &                (W(i,j,k)+W(i,j-1,k))
            END DO
          END DO
          DO i=Istr,Iend
            FC(i,N(ng))=0.0_r8
            FC(i,N(ng)-1)=(cff1*(v(i,j,N(ng)-1,nrhs)+                   &
     &                           v(i,j,N(ng)  ,nrhs))-                  &
     &                     cff2*(v(i,j,N(ng)-2,nrhs)+                   &
     &                           v(i,j,N(ng)  ,nrhs)))*                 &
     &                    (W(i,j,N(ng)-1)+W(i,j-1,N(ng)-1))
            FC(i,1)=(cff1*(v(i,j,1,nrhs)+v(i,j,2,nrhs))-                &
     &               cff2*(v(i,j,1,nrhs)+v(i,j,3,nrhs)))*               &
     &              (W(i,j,1)+W(i,j-1,1))
            FC(i,0)=0.0_r8
          END DO
#  else
          cff1=9.0_r8/16.0_r8
          cff2=1.0_r8/16.0_r8
          DO k=2,N(ng)-2
            DO i=Istr,Iend
              FC(i,k)=(cff1*(v(i,j,k,nrhs)+v(i,j,k+1,nrhs))-            &
     &                 cff2*(v(i,j,k-1,nrhs)+v(i,j,k+2,nrhs)))*         &
     &                (cff1*(W(i,j  ,k)+W(i,j-1,k))-                    &
     &                 cff2*(W(i,j+1,k)+W(i,j-2,k)))
            END DO
          END DO
          DO i=Istr,Iend
            FC(i,N(ng))=0.0_r8
            FC(i,N(ng)-1)=(cff1*(v(i,j,N(ng)-1,nrhs)+                   &
     &                           v(i,j,N(ng)  ,nrhs))-                  &
     &                     cff2*(v(i,j,N(ng)-2,nrhs)+                   &
     &                           v(i,j,N(ng)  ,nrhs)))*                 &
     &                    (cff1*(W(i,j  ,N(ng)-1)+                      &
     &                           W(i,j-1,N(ng)-1))-                     &
     &                     cff2*(W(i,j+1,N(ng)-1)+                      &
     &                           W(i,j-2,N(ng)-1)))
            FC(i,1)=(cff1*(v(i,j,1,nrhs)+v(i,j,2,nrhs))-                &
     &               cff2*(v(i,j,1,nrhs)+v(i,j,3,nrhs)))*               &
     &              (cff1*(W(i,j  ,1)+W(i,j-1,1))-                      &
     &               cff2*(W(i,j+1,1)+W(i,j-2,1)))
            FC(i,0)=0.0_r8
          END DO
#  endif
          DO k=1,N(ng)
            DO i=Istr,Iend
              fac=(FC(i,k)-FC(i,k-1))
              rv(i,j,k,nrhs)=rv(i,j,k,nrhs)-fac
#  ifdef DIAGNOSTICS_UV
              DiaRV(i,j,k,nrhs,M3vadv)=-fac
#  endif
            END DO
          END DO
        END IF
# endif
!
!-----------------------------------------------------------------------
!  Compute forcing term for the 2D momentum equations.
!-----------------------------------------------------------------------
!
!  Vertically integrate baroclinic right-hand-side terms. If not
!  body force stresses, add in the difference between surface and
!  bottom stresses.
!
        DO i=IstrU,Iend
          rufrc(i,j)=ru(i,j,1,nrhs)
# ifdef DIAGNOSTICS_UV
          DiaRUfrc(i,j,3,M2pgrd)=DiaRU(i,j,1,nrhs,M3pgrd)
#  ifdef UV_COR
          DiaRUfrc(i,j,3,M2fcor)=DiaRU(i,j,1,nrhs,M3fcor)
#  endif
#  ifdef UV_ADV
          DiaRUfrc(i,j,3,M2hadv)=DiaRU(i,j,1,nrhs,M3hadv)
#  endif
#  if defined UV_VIS2 || defined UV_VIS4
          DiaRUfrc(i,j,3,M2hvis)=0.0_r8
#  endif
#  ifdef BODYFORCE
          DiaRUfrc(i,j,3,M2strs)=DiaRU(i,j,1,nrhs,M3vvis)
#  endif
# endif
        END DO
        DO k=2,N(ng)
          DO i=IstrU,Iend
            rufrc(i,j)=rufrc(i,j)+ru(i,j,k,nrhs)
# ifdef DIAGNOSTICS_UV
            DiaRUfrc(i,j,3,M2pgrd)=DiaRUfrc(i,j,3,M2pgrd)+              &
     &                             DiaRU(i,j,k,nrhs,M3pgrd)
#  ifdef UV_COR
            DiaRUfrc(i,j,3,M2fcor)=DiaRUfrc(i,j,3,M2fcor)+              &
     &                             DiaRU(i,j,k,nrhs,M3fcor)
#  endif
#  ifdef UV_ADV
            DiaRUfrc(i,j,3,M2hadv)=DiaRUfrc(i,j,3,M2hadv)+              &
     &                             DiaRU(i,j,k,nrhs,M3hadv)
#  endif
#  ifdef BODYFORCE
            DiaRUfrc(i,j,3,M2strs)=DiaRUfrc(i,j,3,M2strs)+              &
     &                             DiaRU(i,j,k,nrhs,M3vvis)
#  endif
# endif
          END DO
        END DO
# ifndef BODYFORCE
        DO i=IstrU,Iend
          fac=(sustr(i,j)-bustr(i,j))*om_u(i,j)*on_u(i,j)
          rufrc(i,j)=rufrc(i,j)+fac
#  ifdef DIAGNOSTICS_UV
          DiaRUfrc(i,j,3,M2strs)=fac
#  endif
        END DO
# endif
        IF (j.ge.JstrV) THEN
          DO i=Istr,Iend
            rvfrc(i,j)=rv(i,j,1,nrhs)
# ifdef DIAGNOSTICS_UV
            DiaRVfrc(i,j,3,M2pgrd)=DiaRV(i,j,1,nrhs,M3pgrd)
#  ifdef UV_COR
            DiaRVfrc(i,j,3,M2fcor)=DiaRV(i,j,1,nrhs,M3fcor)
#  endif
#  ifdef UV_ADV
            DiaRVfrc(i,j,3,M2hadv)=DiaRV(i,j,1,nrhs,M3hadv)
#  endif
#  if defined UV_VIS2 || defined UV_VIS4
            DiaRVfrc(i,j,3,M2hvis)=0.0_r8
#  endif
#  ifdef BODYFORCE
            DiaRVfrc(i,j,3,M2strs)=DiaRV(i,j,1,nrhs,M3vvis)
#  endif
# endif
          END DO
          DO k=2,N(ng)
            DO i=Istr,Iend
              rvfrc(i,j)=rvfrc(i,j)+rv(i,j,k,nrhs)
# ifdef DIAGNOSTICS_UV
              DiaRVfrc(i,j,3,M2pgrd)=DiaRVfrc(i,j,3,M2pgrd)+            &
     &                               DiaRV(i,j,k,nrhs,M3pgrd)
#  ifdef UV_COR
              DiaRVfrc(i,j,3,M2fcor)=DiaRVfrc(i,j,3,M2fcor)+            &
     &                               DiaRV(i,j,k,nrhs,M3fcor)
#  endif
#  ifdef UV_ADV
              DiaRVfrc(i,j,3,M2hadv)=DiaRVfrc(i,j,3,M2hadv)+            &
     &                               DiaRV(i,j,k,nrhs,M3hadv)
#  endif
#  ifdef BODYFORCE
              DiaRVfrc(i,j,3,M2strs)=DiaRVfrc(i,j,3,M2strs)+            &
     &                               DiaRV(i,j,k,nrhs,M3vvis)
#  endif
# endif
            END DO
          END DO
# ifndef BODYFORCE
          DO i=Istr,Iend
            fac=(svstr(i,j)-bvstr(i,j))*om_v(i,j)*on_v(i,j)
            rvfrc(i,j)=rvfrc(i,j)+fac
#  ifdef DIAGNOSTICS_UV
            DiaRVfrc(i,j,3,M2strs)=fac
#  endif
          END DO
# endif
        END IF
      END DO
!
!  Prepare to time-step 2D equations:  set initial free-surface
!  to its fast-time averaged values (which corresponds to the time
!  step "n").
!
      DO k=1,2
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            zeta(i,j,k)=Zt_avg1(i,j)
          END DO
        END DO
# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
        CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,             &
     &                          LBi, UBi, LBj, UBj, zeta(:,:,k))
# endif
      END DO
      RETURN
      END SUBROUTINE rhs3d_tile
#endif
      END MODULE rhs3d_mod
