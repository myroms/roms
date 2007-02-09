#include "cppdefs.h"
      MODULE ini_fields_mod
#ifdef NONLINEAR
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2007 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine initializes other time levels for 2D fields. It also   !
!  couples 3D and 2D momentum equations:  it initializes 2D momentum   !
!  (ubar,vbar) to the vertical integral of initial 3D momentum (u,v).  !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC :: ini_fields

      CONTAINS
!
!***********************************************************************
      SUBROUTINE ini_fields (ng, tile, model)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
# ifdef SOLVE3D
      USE mod_coupling
# endif
      USE mod_ocean
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
# include "tile.h"
!
      CALL ini_fields_tile (ng, model, Istr, Iend, Jstr, Jend,          &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      kstp(ng), krhs(ng), knew(ng),               &
# ifdef SOLVE3D
     &                      nstp(ng), nnew(ng),                         &
# endif
# ifdef MASKING
     &                      GRID(ng) % rmask,                           &
     &                      GRID(ng) % umask,                           &
     &                      GRID(ng) % vmask,                           &
# endif
# if defined PERFECT_RESTART && \
    (defined EW_PERIODIC || defined NS_PERIODIC)
#  ifdef SOLVE3D
     &                      OCEAN(ng) % ru,                             &
     &                      OCEAN(ng) % rv,                             &
#  endif
     &                      OCEAN(ng) % rubar,                          &
     &                      OCEAN(ng) % rvbar,                          &
     &                      OCEAN(ng) % rzeta,                          &
# endif
# ifdef SOLVE3D
#  if defined SEDIMENT && defined SED_MORPH
     &                      OCEAN(ng) %  bed,                           &
     &                      GRID(ng) % bed_thick0,                      &
#  endif
     &                      GRID(ng) % Hz,                              &
     &                      GRID(ng) % h,                               &
#  ifdef ICESHELF
     &                      GRID(ng) % zice,                            &
#  endif
     &                      GRID(ng) % z_r,                             &
     &                      GRID(ng) % z_w,                             &
     &                      COUPLING(ng) % Zt_avg1,                     &
     &                      OCEAN(ng) % t,                              &
     &                      OCEAN(ng) % u,                              &
     &                      OCEAN(ng) % v,                              &
# endif 
     &                      OCEAN(ng) % ubar,                           &
     &                      OCEAN(ng) % vbar,                           &
     &                      OCEAN(ng) % zeta)

      RETURN
      END SUBROUTINE ini_fields
!
!***********************************************************************
      SUBROUTINE ini_fields_tile (ng, model, Istr, Iend, Jstr, Jend,    &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            kstp, krhs, knew,                     &
# ifdef SOLVE3D
     &                            nstp, nnew,                           &
# endif
# ifdef MASKING
     &                            rmask, umask, vmask,                  &
# endif
# if defined PERFECT_RESTART && \
    (defined EW_PERIODIC || defined NS_PERIODIC)
#  ifdef SOLVE3D
     &                            ru, rv,                               &
#  endif
     &                            rubar, rvbar, rzeta,                  &
# endif
# ifdef SOLVE3D
#  if defined SEDIMENT && defined SED_MORPH
     &                            bed, bed_thick0,                      &
#  endif
     &                            Hz, h,                                &
#  ifdef ICESHELF
     &                            zice,                                 &
#  endif
     &                            z_r, z_w,                             &
     &                            Zt_avg1,                              &
     &                            t, u, v,                              &
# endif
     &                            ubar, vbar, zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
# if defined SEDIMENT && defined SED_MORPH
      USE mod_sediment
# endif
!
# if defined EW_PERIODIC || defined NS_PERIODIC
      USE exchange_2d_mod
#  ifdef SOLVE3D
      USE exchange_3d_mod
#  endif
# endif
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
#  ifdef SOLVE3D
      USE mp_exchange_mod, ONLY : mp_exchange3d, mp_exchange4d
#  else
#   if defined PERFECT_RESTART && \
      (defined EW_PERIODIC || defined NS_PERIODIC)
      USE mp_exchange_mod, ONLY : mp_exchange3d
#   endif
#  endif
# endif
# ifdef SOLVE3D
      USE set_depth_mod, ONLY : set_depth_tile
      USE t3dbc_mod, ONLY : t3dbc_tile
      USE u3dbc_mod, ONLY : u3dbc_tile
      USE v3dbc_mod, ONLY : v3dbc_tile
# endif
      USE u2dbc_mod, ONLY : u2dbc_tile
      USE v2dbc_mod, ONLY : v2dbc_tile
      USE zetabc_mod, ONLY : zetabc_tile
!
!  Imported variable declarations.
!
# ifdef DISTRIBUTE
#  ifdef EW_PERIODIC
      logical :: EWperiodic=.TRUE.
#  else
      logical :: EWperiodic=.FALSE.
#  endif
#  ifdef NS_PERIODIC
      logical :: NSperiodic=.TRUE.
#  else
      logical :: NSperiodic=.FALSE.
#  endif
# endif
      integer, intent(in) :: ng, model, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: kstp, krhs, knew
# ifdef SOLVE3D
      integer, intent(in) :: nstp, nnew
# endif
!
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
#   ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:,LBj:)
#   endif
#  endif
#  if defined PERFECT_RESTART && \
     (defined EW_PERIODIC || defined NS_PERIODIC)
#   ifdef SOLVE3D
      real(r8), intent(inout) :: ru(LBi:,LBj:,0:,:)
      real(r8), intent(inout) :: rv(LBi:,LBj:,0:,:)
#   endif
      real(r8), intent(inout) :: rubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: rvbar(LBi:,LBj:,:)
      real(r8), intent(inout) :: rzeta(LBi:,LBj:,:)
#  endif
#  ifdef SOLVE3D
      real(r8), intent(inout) :: h(LBi:,LBj:)
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: v(LBi:,LBj:,:,:)

      real(r8), intent(inout) :: Hz(LBi:,LBj:,:)
      real(r8), intent(inout) :: Zt_avg1(LBi:,LBj:)
      real(r8), intent(inout) :: z_r(LBi:,LBj:,:)
      real(r8), intent(inout) :: z_w(LBi:,LBj:,0:)
#  endif
      real(r8), intent(inout) :: ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: vbar(LBi:,LBj:,:)
      real(r8), intent(inout) :: zeta(LBi:,LBj:,:)
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
#   ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:UBi,LBj:UBj)
#   endif
#  endif
#  if defined PERFECT_RESTART && \
     (defined EW_PERIODIC || defined NS_PERIODIC)
#   ifdef SOLVE3D
      real(r8), intent(inout) :: ru(LBi:UBi,LBj:UBj,0:N(ng),2)
      real(r8), intent(inout) :: rv(LBi:UBi,LBj:UBj,0:N(ng),2)
#   endif
      real(r8), intent(inout) :: rubar(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: rvbar(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: rzeta(LBi:UBi,LBj:UBj,2)
#  endif
#  ifdef SOLVE3D
      real(r8), intent(inout) :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: v(LBi:UBi,LBj:UBj,N(ng),2)

      real(r8), intent(inout) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: Zt_avg1(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
#  endif
      real(r8), intent(inout) :: ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: zeta(LBi:UBi,LBj:UBj,3)
# endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, itrc, j, k, kbed

      real(r8) :: cff1, cff2
# ifdef SOLVE3D
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N(ng)) :: CF
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N(ng)) :: DC
# endif

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Initialize other free-surface time-levels.
!-----------------------------------------------------------------------
!
      IF (.not.PerfectRST(ng)) THEN
        DO j=Jstr,Jend
          DO i=Istr,Iend
            cff1=zeta(i,j,kstp)
# ifdef MASKING
            cff1=cff1*rmask(i,j)
# endif
            zeta(i,j,kstp)=cff1
# ifdef SOLVE3D
            zeta(i,j,knew)=cff1
# else
            zeta(i,j,krhs)=cff1
# endif
          END DO
        END DO
!
!  Apply boundary conditions.
!
        CALL zetabc_tile (ng, Istr, Iend, Jstr, Jend,                   &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    krhs, kstp, kstp,                             &
     &                    zeta)
# ifdef SOLVE3D
        CALL zetabc_tile (ng, Istr, Iend, Jstr, Jend,                   &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    krhs, kstp, knew,                             &
     &                    zeta)
# else
        CALL zetabc_tile (ng, Istr, Iend, Jstr, Jend,                   &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    krhs, kstp, krhs,                             &
     &                    zeta)
# endif
      END IF

# if defined EW_PERIODIC || defined NS_PERIODIC
!
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        zeta(:,:,kstp))
#  ifdef SOLVE3D
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        zeta(:,:,knew))
#  else
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        zeta(:,:,krhs))
#  endif
# endif
# ifdef DISTRIBUTE
#  ifdef SOLVE3D
      CALL mp_exchange2d (ng, model, 2, Istr, Iend, Jstr, Jend,         &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    zeta(:,:,kstp),                               &
     &                    zeta(:,:,knew))
#  else
      CALL mp_exchange2d (ng, model, 2, Istr, Iend, Jstr, Jend,         &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    zeta(:,:,kstp),                               &
     &                    zeta(:,:,krhs))
#  endif
# endif

# ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Initialize fast-time averaged free-surface (Zt_avg1) with the inital
!  free-surface
!-----------------------------------------------------------------------
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          Zt_avg1(i,j)=zeta(i,j,kstp)
        END DO
      END DO

#  if defined EW_PERIODIC || defined NS_PERIODIC
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Zt_avg1)
#  endif
#  ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, model, 1, Istr, Iend, Jstr, Jend,         &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    Zt_avg1)
#  endif

#  if defined SEDIMENT && defined SED_MORPH
!
!-----------------------------------------------------------------------
!  Compute initial total thickness for all sediment bed layers.
!-----------------------------------------------------------------------
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          bed_thick0(i,j)=0.0_r8
          DO kbed=1,Nbed
            bed_thick0(i,j)=bed_thick0(i,j)+bed(i,j,kbed,ithck)
          END DO
        END DO
      END DO

#   if defined EW_PERIODIC || defined NS_PERIODIC
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        bed_thick0)
#   endif
#   ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, model, 1, Istr, Iend, Jstr, Jend,         &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    bed_thick0)
#   endif
#  endif
# endif

# ifndef SOLVE3D
!
!-----------------------------------------------------------------------
!  Initialize other 2D momentum time-levels.
!-----------------------------------------------------------------------
!
      IF (.not.PerfectRST(ng)) THEN
        DO j=Jstr,Jend
          DO i=IstrU,Iend
            cff1=ubar(i,j,kstp)
#  ifdef MASKING
            cff1=cff1*umask(i,j)
#  endif
            ubar(i,j,kstp)=cff1
            ubar(i,j,krhs)=cff1
          END DO
          IF (j.ge.JstrV) THEN
            DO i=Istr,Iend
              cff2=vbar(i,j,kstp)
#  ifdef MASKING
              cff2=cff2*vmask(i,j)
#  endif
              vbar(i,j,kstp)=cff2
              vbar(i,j,krhs)=cff2
            END DO
          END IF
        END DO
!
!  Apply boundary conditions.
!
        CALL u2dbc_tile (ng, Istr, Iend, Jstr, Jend,                    &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   krhs, kstp, kstp,                              &
     &                   ubar, vbar, zeta)
        CALL v2dbc_tile (ng, Istr, Iend, Jstr, Jend,                    &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   krhs, kstp, kstp,                              &
     &                   ubar, vbar, zeta)
        CALL u2dbc_tile (ng, Istr, Iend, Jstr, Jend,                    &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   krhs, kstp, krhs,                              &
     &                   ubar, vbar, zeta)
        CALL v2dbc_tile (ng, Istr, Iend, Jstr, Jend,                    &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   krhs, kstp, krhs,                              &
     &                   ubar, vbar, zeta)
      END IF

#  if defined EW_PERIODIC || defined NS_PERIODIC
!
      CALL exchange_u2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        ubar(:,:,kstp))
      CALL exchange_v2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        vbar(:,:,kstp))
      CALL exchange_u2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        ubar(:,:,krhs))
      CALL exchange_v2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        vbar(:,:,krhs))
#  endif
#  ifdef DISTRIBUTE
!
      CALL mp_exchange2d (ng, model, 4, Istr, Iend, Jstr, Jend,         &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    ubar(:,:,kstp), vbar(:,:,kstp),               &
     &                    ubar(:,:,krhs), vbar(:,:,krhs))
#  endif

# else
!
!-----------------------------------------------------------------------
!  Initialize other 3D momentum time-levels.
!-----------------------------------------------------------------------
!
      IF (.not.PerfectRST(ng)) THEN
        DO j=Jstr,Jend
          DO k=1,N(ng)
            DO i=IstrU,Iend
              cff1=u(i,j,k,nstp)
#  ifdef MASKING
              cff1=cff1*umask(i,j)
#  endif
              u(i,j,k,nstp)=cff1
              u(i,j,k,nnew)=cff1
            END DO
          END DO
          IF (j.ge.JstrV) THEN
            DO k=1,N(ng)
              DO i=Istr,Iend
                cff2=v(i,j,k,nstp)
#  ifdef MASKING
                cff2=cff2*vmask(i,j)
#  endif
                v(i,j,k,nstp)=cff2
                v(i,j,k,nnew)=cff2
              END DO
            END DO
          END IF
        END DO
!
!  Apply boundary conditions.
!
        CALL u3dbc_tile (ng, Istr, Iend, Jstr, Jend,                    &
     &                   LBi, UBi, LBj, UBj, N(ng),                     &
     &                   nstp, nstp,                                    &
     &                   u)
        CALL v3dbc_tile (ng, Istr, Iend, Jstr, Jend,                    &
     &                   LBi, UBi, LBj, UBj, N(ng),                     &
     &                   nstp, nstp,                                    &
     &                   v)

        CALL u3dbc_tile (ng, Istr, Iend, Jstr, Jend,                    &
     &                   LBi, UBi, LBj, UBj, N(ng),                     &
     &                   nstp, nnew,                                    &
     &                   u)
        CALL v3dbc_tile (ng, Istr, Iend, Jstr, Jend,                    &
     &                   LBi, UBi, LBj, UBj, N(ng),                     &
     &                   nstp, nnew,                                    &
     &                   v)
      END IF

#  if defined EW_PERIODIC || defined NS_PERIODIC
!
      CALL exchange_u3d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        u(:,:,:,nstp))
      CALL exchange_v3d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        v(:,:,:,nstp))

      CALL exchange_u3d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        u(:,:,:,nnew))
      CALL exchange_v3d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        v(:,:,:,nnew))
#  endif
#  ifdef DISTRIBUTE
!
      CALL mp_exchange3d (ng, model, 4, Istr, Iend, Jstr, Jend,         &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    u(:,:,:,nstp), v(:,:,:,nstp),                 &
     &                    u(:,:,:,nnew), v(:,:,:,nnew))
#  endif
!
!-----------------------------------------------------------------------
!  Initialize other tracers time-levels.
!-----------------------------------------------------------------------
!
      IF (.not.PerfectRST(ng)) THEN
        DO itrc=1,NT(ng)
          DO k=1,N(ng)
            DO j=Jstr,Jend
              DO i=Istr,Iend
                cff1=t(i,j,k,nstp,itrc)
#  ifdef MASKING
                cff1=cff1*rmask(i,j)
#  endif
                t(i,j,k,nstp,itrc)=cff1
                t(i,j,k,nnew,itrc)=cff1
              END DO
            END DO
          END DO
!
!  Apply boundary conditions.
!
          CALL t3dbc_tile (ng, Istr, Iend, Jstr, Jend, itrc,            &
     &                     LBi, UBi, LBj, UBj, N(ng), NT(ng),           &
     &                     nstp, nstp,                                  &
     &                     t)
          CALL t3dbc_tile (ng, Istr, Iend, Jstr, Jend, itrc,            &
     &                     LBi, UBi, LBj, UBj, N(ng), NT(ng),           &
     &                     nstp, nnew,                                  &
     &                     t)
        END DO
      END IF

#  if defined EW_PERIODIC || defined NS_PERIODIC
!
      DO itrc=1,NT(ng)
        CALL exchange_r3d_tile (ng, Istr, Iend, Jstr, Jend,             &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          t(:,:,:,nstp,itrc))
        CALL exchange_r3d_tile (ng, Istr, Iend, Jstr, Jend,             &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          t(:,:,:,nnew,itrc))
      END DO
#  endif
#  ifdef DISTRIBUTE
      CALL mp_exchange4d (ng, model, 2, Istr, Iend, Jstr, Jend,         &
     &                    LBi, UBi, LBj, UBj, 1, N(ng), 1, NT(ng),      &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    t(:,:,:,nstp,:),                              &
     &                    t(:,:,:,nnew,:))
#  endif
!
!-----------------------------------------------------------------------
!  Compute initial depths and thicknesses.
!-----------------------------------------------------------------------
!
      CALL set_depth_tile (ng, Istr, Iend, Jstr, Jend,                  &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     h,                                           &
#  ifdef ICESHELF
     &                     zice,                                        &
#  endif
#  if defined SEDIMENT && defined SED_MORPH
     &                     bed, bed_thick0,                             &
#  endif
     &                     Zt_avg1,                                     &
     &                     Hz, z_r, z_w)
!
!-----------------------------------------------------------------------
!  Compute vertically integrated momentum (ubar,vbar) from initial 3D
!  momentum (u,v).
!-----------------------------------------------------------------------
!
!  Compute adjoint 2D velocity component in the XI-direction.  Here
!  DC(i,1:N) are the thicknesses of U-boxes, DC(i,0) is total depth of
!  the water column, and CF(i,0) is the vertical integral.
!
      IF (.not.PerfectRST(ng)) THEN
        DO j=Jstr,Jend
          DO i=IstrU,Iend
            DC(i,0)=0.0_r8
            CF(i,0)=0.0_r8
          END DO
          DO k=1,N(ng)
            DO i=IstrU,Iend
              DC(i,k)=0.5_r8*(Hz(i,j,k)+Hz(i-1,j,k))
              DC(i,0)=DC(i,0)+DC(i,k)
              CF(i,0)=CF(i,0)+DC(i,k)*u(i,j,k,nstp)
            END DO
          END DO
          DO i=IstrU,Iend
            cff1=1.0_r8/DC(i,0)
            cff2=CF(i,0)*cff1
#  ifdef MASKING
            cff2=cff2*umask(i,j)
#  endif
            ubar(i,j,kstp)=cff2
            ubar(i,j,knew)=cff2
          END DO
!
!  Compute adjoint 2D velocity component in the ETA-direction.  Here
!  DC(i,1:N) are the thicknesses of V-boxes, DC(i,0) is total depth of
!  the water column, and CF(i,0) is the vertical integral.
!
          IF (j.ge.Jstr) THEN
            DO i=Istr,Iend
              DC(i,0)=0.0_r8
              CF(i,0)=0.0_r8
            END DO
            DO k=1,N(ng)
              DO i=Istr,Iend
                DC(i,k)=0.5_r8*(Hz(i,j,k)+Hz(i,j-1,k))
                DC(i,0)=DC(i,0)+DC(i,k)
                CF(i,0)=CF(i,0)+DC(i,k)*v(i,j,k,nstp)
              END DO
            END DO
            DO i=Istr,Iend
              cff1=1.0_r8/DC(i,0)
              cff2=CF(i,0)*cff1
#  ifdef MASKING
              cff2=cff2*vmask(i,j)
#  endif
              vbar(i,j,kstp)=cff2
              vbar(i,j,knew)=cff2
            END DO
          END IF   
        END DO
!
!  Apply boundary conditions.
!
        CALL u2dbc_tile (ng, Istr, Iend, Jstr, Jend,                    &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   krhs, kstp, kstp,                              &
     &                   ubar, vbar, zeta)
        CALL v2dbc_tile (ng, Istr, Iend, Jstr, Jend,                    &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   krhs, kstp, kstp,                              &
     &                   ubar, vbar, zeta)
        CALL u2dbc_tile (ng, Istr, Iend, Jstr, Jend,                    &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   krhs, kstp, knew,                              &
     &                   ubar, vbar, zeta)
        CALL v2dbc_tile (ng, Istr, Iend, Jstr, Jend,                    &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   krhs, kstp, knew,                              &
     &                   ubar, vbar, zeta)
      END IF

#  if defined EW_PERIODIC || defined NS_PERIODIC
!
      CALL exchange_u2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        ubar(:,:,kstp))
      CALL exchange_v2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        vbar(:,:,kstp))
      CALL exchange_u2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        ubar(:,:,knew))
      CALL exchange_v2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        vbar(:,:,knew))
#  endif
#  ifdef DISTRIBUTE
!
      CALL mp_exchange2d (ng, model, 4, Istr, Iend, Jstr, Jend,         &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    ubar(:,:,kstp), vbar(:,:,kstp),               &
     &                    ubar(:,:,knew), vbar(:,:,knew))
#  endif
# endif
# if defined PERFECT_RESTART && \
    (defined EW_PERIODIC || defined NS_PERIODIC)
!
!-----------------------------------------------------------------------
!  If perfect restart, apply periodic boundary condition to
!  right-hand-side terms.
!-----------------------------------------------------------------------
!
      IF (PerfectRST(ng)) THEN
        DO i=1,2
          CALL exchange_u2d_tile (ng, Istr, Iend, Jstr, Jend,           &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            rubar(:,:,i))
          CALL exchange_v2d_tile (ng, Istr, Iend, Jstr, Jend,           &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            rvbar(:,:,i))
          CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,           &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            rzeta(:,:,i))
#  ifdef SOLVE3D
          CALL exchange_u3d_tile (ng, Istr, Iend, Jstr, Jend,           &
     &                            LBi, UBi, LBj, UBj, 0, N(ng),         &
     &                            ru(:,:,:,i))
          CALL exchange_v3d_tile (ng, Istr, Iend, Jstr, Jend,           &
     &                            LBi, UBi, LBj, UBj, 0, N(ng),         &
     &                            rv(:,:,:,i))
#  endif
        END DO
#  ifdef DISTRIBUTE
        CALL mp_exchange3d (ng, model, 3, Istr, Iend, Jstr, Jend,       &
     &                      LBi, UBi, LBj, UBj, 1, 2,                   &
     &                      NghostPoints, EWperiodic, NSperiodic,       &
     &                      rubar, rvbar, rzeta)
#   ifdef SOLVE3D
        CALL mp_exchange4d (ng, model, 2, Istr, Iend, Jstr, Jend,       &
     &                      LBi, UBi, LBj, UBj, 0, N(ng), 1, 2,         &
     &                      NghostPoints, EWperiodic, NSperiodic,       &
     &                      ru, rv)
#   endif
#  endif
      END IF
# endif

      RETURN
      END SUBROUTINE ini_fields_tile
#endif
      END MODULE ini_fields_mod
