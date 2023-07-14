#include "cppdefs.h"
      MODULE tl_zetabc_mod
#ifdef TANGENT
!
!svn $Id: tl_zetabc.F 275 2023-02-09 03:21:09Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine sets tangent linear lateral boundary conditions for    !
!  free-surface.                                                       !
# if defined STEP2D_FB_AB3_AM4 || defined STEP2D_FB_LF_AM3
!
!  Notice that "zetabc_local" is used for either the Forward-Backward  !
!  AB3-AM4 or Forward-Backward LF-AM3 barotropic kernels where the     !
!  boundary conditions are loaded into private array "zeta_new".       !
# endif
!                                                                      !
!  BASIC STATE variables fields needed: zeta                           !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_boundary
      USE mod_grid
      USE mod_ncparam
      USE mod_ocean
      USE mod_scalars
      USE mod_stepping
!
      implicit none
!
      PRIVATE
      PUBLIC  :: tl_zetabc
# if defined STEP2D_FB_AB3_AM4 || defined STEP2D_FB_LF_AM3
      PUBLIC  :: tl_zetabc_local
# endif
      PUBLIC  :: tl_zetabc_tile
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE tl_zetabc (ng, tile, kout)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, kout
!
!  Local variable declarations.
!
# include "tile.h"
!
      CALL tl_zetabc_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     krhs(ng), kstp(ng), kout,                    &
     &                     OCEAN(ng) % zeta,                            &
     &                     OCEAN(ng) % tl_zeta)

      RETURN
      END SUBROUTINE tl_zetabc
!
!***********************************************************************
      SUBROUTINE tl_zetabc_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           krhs, kstp, kout,                      &
     &                           zeta, tl_zeta)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: krhs, kstp, kout
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in   ) :: zeta(LBi:,LBj:,:)

      real(r8), intent(inout) :: tl_zeta(LBi:,LBj:,:)
# else
      real(r8), intent(in   ) :: zeta(LBi:UBi,LBj:UBj,:)

      real(r8), intent(inout) :: tl_zeta(LBi:UBi,LBj:UBj,:)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, know

      real(r8) :: Ce, Cx
      real(r8) :: cff, cff1, cff2, dt2d, tau

      real(r8) :: tl_Ce, tl_Cx
      real(r8) :: tl_cff1, tl_cff2

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_grad

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set time-indices
!-----------------------------------------------------------------------
!
# if defined STEP2D_FB_AB3_AM4 || defined STEP2D_FB_LF_AM3
        know=kstp
        dt2d=dtfast(ng)
# else
      IF (FIRST_2D_STEP) THEN
        know=krhs
        dt2d=dtfast(ng)
      ELSE IF (PREDICTOR_2D_STEP(ng)) THEN
        know=krhs
        dt2d=2.0_r8*dtfast(ng)
      ELSE
        know=kstp
        dt2d=dtfast(ng)
      END IF
# endif
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the western edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
!
!  Western edge, implicit upstream radiation condition.
!  (Not implemented since the TLM requires nonlinear celerity at every
!   timestep. Using gradient instead).
!
        IF (tl_LBC(iwest,isFsur,ng)%radiation) THEN
          IF (iic(ng).ne.0) THEN
# if defined CELERITY_READ && defined FORWARD_READ
            DO j=Jstr,Jend+1
!^            grad(Istr-1,j)=zeta(Istr-1,j  ,know)-                     &
!^   &                       zeta(Istr-1,j-1,know)
!^
              tl_grad(Istr-1,j)=0.0_r8
            END DO
# endif
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%west(j)) THEN
# if defined CELERITY_READ && defined FORWARD_READ
                IF (tl_LBC(iwest,isFsur,ng)%nudging) THEN
                  IF (BOUNDARY(ng)%zeta_west_Cx(j).eq.0.0_r8) THEN
                    tau=FSobc_in(ng,iwest)
                  ELSE
                    tau=FSobc_out(ng,iwest)
                  END IF
                  tau=tau*dt2d
                END IF
                Cx=BOUNDARY(ng)%zeta_west_Cx(j)
#  ifdef RADIATION_2D
                Ce=BOUNDARY(ng)%zeta_west_Ce(j)
#  else
                Ce=0.0_r8
#  endif
                cff=BOUNDARY(ng)%zeta_west_C2(j)
!^              zeta(Istr-1,j,kout)=(cff*zeta(Istr-1,j,know)+           &
!^   &                               Cx *zeta(Istr  ,j,kout)-           &
!^   &                               MAX(Ce,0.0_r8)*grad(Istr-1,j  )-   &
!^   &                               MIN(Ce,0.0_r8)*grad(Istr-1,j+1))/  &
!^   &                              (cff+Cx)
!^
                tl_zeta(Istr-1,j,kout)=(cff*tl_zeta(Istr-1,j,know)+     &
     &                                  Cx *tl_zeta(Istr  ,j,kout)-     &
     &                                  MAX(Ce,0.0_r8)*                 &
     &                                     tl_grad(Istr-1,j  )-         &
     &                                  MIN(Ce,0.0_r8)*                 &
     &                                     tl_grad(Istr-1,j+1))/        &
     &                                 (cff+Cx)

                IF (tl_LBC(iwest,isFsur,ng)%nudging) THEN
!^                zeta(Istr-1,j,kout)=zeta(Istr-1,j,kout)+              &
!^   &                                tau*(BOUNDARY(ng)%zeta_west(j)-   &
!^   &                                     zeta(Istr-1,j,know))
!^
                  tl_zeta(Istr-1,j,kout)=tl_zeta(Istr-1,j,kout)-        &
     &                                   tau*tl_zeta(Istr-1,j,know)
                END IF
# else
!^              zeta(Istr-1,j,kout)=zeta(Istr,j,kout)         ! gradient
!^
                tl_zeta(Istr-1,j,kout)=tl_zeta(Istr,j,kout)
# endif
# ifdef MASKING
!^              zeta(Istr-1,j,kout)=zeta(Istr-1,j,kout)*                &
!^   &                              GRID(ng)%rmask(Istr-1,j)
!^
                tl_zeta(Istr-1,j,kout)=tl_zeta(Istr-1,j,kout)*          &
     &                                 GRID(ng)%rmask(Istr-1,j)
# endif
              END IF
            END DO
          END IF
!
!  Western edge, explicit Chapman boundary condition.
!
        ELSE IF (tl_LBC(iwest,isFsur,ng)%Chapman_explicit) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%west(j)) THEN
              cff=dt2d*GRID(ng)%pm(Istr,j)
              cff1=SQRT(g*(GRID(ng)%h(Istr,j)+                          &
     &                     zeta(Istr,j,know)))
              tl_cff1=0.5_r8*g*(GRID(ng)%tl_h(Istr,j)+                  &
     &                          tl_zeta(Istr,j,know))/cff1
              Cx=cff*cff1
              tl_Cx=cff*tl_cff1
!^            zeta(Istr-1,j,kout)=(1.0_r8-Cx)*zeta(Istr-1,j,know)+      &
!^   &                            Cx*zeta(Istr,j,know)
!^
              tl_zeta(Istr-1,j,kout)=(1.0_r8-Cx)*tl_zeta(Istr-1,j,know)+&
     &                               tl_Cx*(zeta(Istr-1,j,know)+        &
     &                                      zeta(Istr  ,j,know))+       &
     &                               Cx*tl_zeta(Istr,j,know)
# ifdef MASKING
!^            zeta(Istr-1,j,kout)=zeta(Istr-1,j,kout)*                  &
!^   &                            GRID(ng)%rmask(Istr-1,j)
!^
              tl_zeta(Istr-1,j,kout)=tl_zeta(Istr-1,j,kout)*            &
     &                               GRID(ng)%rmask(Istr-1,j)
# endif
            END IF
          END DO
!
!  Western edge, implicit Chapman boundary condition.
!
        ELSE IF (tl_LBC(iwest,isFsur,ng)%Chapman_implicit) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%west(j)) THEN
              cff=dt2d*GRID(ng)%pm(Istr,j)
              cff1=SQRT(g*(GRID(ng)%h(Istr,j)+                          &
     &                     zeta(Istr,j,know)))
              tl_cff1=0.5_r8*g*(GRID(ng)%tl_h(Istr,j)+                  &
     &                          tl_zeta(Istr,j,know))/cff1
              Cx=cff*cff1
              tl_Cx=cff*tl_cff1
              cff2=1.0_r8/(1.0_r8+Cx)
              tl_cff2=-cff2*cff2*tl_Cx
!^            zeta(Istr-1,j,kout)=cff2*(zeta(Istr-1,j,know)+            &
!^   &                                  Cx*zeta(Istr,j,kout))
!^
              tl_zeta(Istr-1,j,kout)=tl_cff2*(zeta(Istr-1,j,know)+      &
     &                                        Cx*zeta(Istr,j,kout))+    &
     &                               cff2*(tl_zeta(Istr-1,j,know)+      &
     &                                     tl_Cx*zeta(Istr,j,kout)+     &
     &                                     Cx*tl_zeta(Istr,j,kout))
# ifdef MASKING
!^            zeta(Istr-1,j,kout)=zeta(Istr-1,j,kout)*                  &
!^   &                            GRID(ng)%rmask(Istr-1,j)
!^
              tl_zeta(Istr-1,j,kout)=tl_zeta(Istr-1,j,kout)*            &
     &                               GRID(ng)%rmask(Istr-1,j)
# endif
            END IF
          END DO
!
!  Western edge, clamped boundary condition.
!
        ELSE IF (tl_LBC(iwest,isFsur,ng)%clamped) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%west(j)) THEN
!^            zeta(Istr-1,j,kout)=BOUNDARY(ng)%zeta_west(j)
!^
# ifdef ADJUST_BOUNDARY
              IF (Lobc(iwest,isFsur,ng)) THEN
                tl_zeta(Istr-1,j,kout)=BOUNDARY(ng)%tl_zeta_west(j)
              ELSE
                tl_zeta(Istr-1,j,kout)=0.0_r8
              END IF
# else
              tl_zeta(Istr-1,j,kout)=0.0_r8
# endif
# ifdef MASKING
!^            zeta(Istr-1,j,kout)=zeta(Istr-1,j,kout)*                  &
!^   &                            GRID(ng)%rmask(Istr-1,j)
!^
              tl_zeta(Istr-1,j,kout)=tl_zeta(Istr-1,j,kout)*            &
     &                               GRID(ng)%rmask(Istr-1,j)
# endif
            END IF
          END DO
!
!  Western edge, gradient boundary condition.
!
        ELSE IF (tl_LBC(iwest,isFsur,ng)%gradient) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%west(j)) THEN
!^            zeta(Istr-1,j,kout)=zeta(Istr,j,kout)
!^
              tl_zeta(Istr-1,j,kout)=tl_zeta(Istr,j,kout)
# ifdef MASKING
!^            zeta(Istr-1,j,kout)=zeta(Istr-1,j,kout)*                  &
!^   &                            GRID(ng)%rmask(Istr-1,j)
!^
              tl_zeta(Istr-1,j,kout)=tl_zeta(Istr-1,j,kout)*            &
     &                               GRID(ng)%rmask(Istr-1,j)
# endif
            END IF
          END DO
!
!  Western edge, closed boundary condition.
!
        ELSE IF (tl_LBC(iwest,isFsur,ng)%closed) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%west(j)) THEN
!^            zeta(Istr-1,j,kout)=zeta(Istr,j,kout)
!^
              tl_zeta(Istr-1,j,kout)=tl_zeta(Istr,j,kout)
# ifdef MASKING
!^            zeta(Istr-1,j,kout)=zeta(Istr-1,j,kout)*                  &
!^   &                            GRID(ng)%rmask(Istr-1,j)
!^
              tl_zeta(Istr-1,j,kout)=tl_zeta(Istr-1,j,kout)*            &
     &                               GRID(ng)%rmask(Istr-1,j)
# endif
            END IF
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the eastern edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
!
!  Eastern edge, implicit upstream radiation condition.
!  (Not implemented since the TLM requires nonlinear celerity at every
!   timestep. Using gradient instead).
!
        IF (tl_LBC(ieast,isFsur,ng)%radiation) THEN
          IF (iic(ng).ne.0) THEN
# if defined CELERITY_READ && defined FORWARD_READ
            DO j=Jstr,Jend+1
!^            grad(Iend+1,j)=zeta(Iend+1,j  ,know)-                     &
!^   &                       zeta(Iend+1,j-1,know)
!^
              tl_grad(Iend+1,j)=0.0_r8
            END DO
# endif
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%east(j)) THEN
# if defined CELERITY_READ && defined FORWARD_READ
                IF (tl_LBC(ieast,isFsur,ng)%nudging) THEN
                  IF (BOUNDARY(ng)%zeta_east_Cx(j).eq.0.0_r8) THEN
                    tau=FSobc_in(ieast)
                  ELSE
                    tau=FSobc_out(ieast)
                  END IF
                  tau=tau*dt2d
                END IF
                Cx=BOUNDARY(ng)%zeta_east_Cx(j)
#  ifdef RADIATION_2D
                Ce=BOUNDARY(ng)%zeta_east_Ce(j)
#  else
                Ce=0.0_r8
#  endif
                cff=BOUNDARY(ng)%zeta_east_C2(j)

!^              zeta(Iend+1,j,kout)=(cff*zeta(Iend+1,j,know)+           &
!^   &                               Cx *zeta(Iend  ,j,kout)-           &
!^   &                               MAX(Ce,0.0_r8)*grad(Iend+1,j  )-   &
!^   &                               MIN(Ce,0.0_r8)*grad(Iend+1,j+1))/  &
!^   &                              (cff+Cx)
!^
                tl_zeta(Iend+1,j,kout)=(cff*tl_zeta(Iend+1,j,know)+     &
     &                                  Cx *tl_zeta(Iend  ,j,kout)-     &
     &                                  MAX(Ce,0.0_r8)*                 &
     &                                     tl_grad(Iend+1,j  )-         &
     &                                  MIN(Ce,0.0_r8)*                 &
     &                                     tl_grad(Iend+1,j+1))/        &
     &                                 (cff+Cx)

                IF (tl_LBC(ieast,isFsur,ng)%nudging) THEN
!^                zeta(Iend+1,j,kout)=zeta(Iend+1,j,kout)+              &
!^   &                                tau*(BOUNDARY(ng)%zeta_east(j)-   &
!^   &                                     zeta(Iend+1,j,know))
!^
                  tl_zeta(Iend+1,j,kout)=tl_zeta(Iend+1,j,kout)-        &
     &                                   tau*tl_zeta(Iend+1,j,know)
                END IF
# else
!^              zeta(Iend+1,j,kout)=zeta(Iend,j,kout)         ! gradient
!^
                tl_zeta(Iend+1,j,kout)=tl_zeta(Iend,j,kout)
# endif
# ifdef MASKING
!^              zeta(Iend+1,j,kout)=zeta(Iend+1,j,kout)*                &
!^   &                              GRID(ng)%rmask(Iend+1,j)
!^
                tl_zeta(Iend+1,j,kout)=tl_zeta(Iend+1,j,kout)*          &
     &                                 GRID(ng)%rmask(Iend+1,j)
# endif
              END IF
            END DO
          END IF
!
!  Eastern edge, explicit Chapman boundary condition.
!
        ELSE IF (tl_LBC(ieast,isFsur,ng)%Chapman_explicit) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%east(j)) THEN
              cff=dt2d*GRID(ng)%pm(Iend,j)
              cff1=SQRT(g*(GRID(ng)%h(Iend,j)+                          &
     &                     zeta(Iend,j,know)))
              tl_cff1=0.5_r8*g*(GRID(ng)%tl_h(Iend,j)+                  &
     &                          tl_zeta(Iend,j,know))/cff1
              Cx=cff*cff1
              tl_Cx=cff*tl_cff1
!^            zeta(Iend+1,j,kout)=(1.0_r8-Cx)*zeta(Iend+1,j,know)+      &
!^   &                            Cx*zeta(Iend,j,know)
!^
              tl_zeta(Iend+1,j,kout)=(1.0_r8-Cx)*tl_zeta(Iend+1,j,know)+&
     &                               tl_Cx*(zeta(Iend+1,j,know)+        &
     &                                      zeta(Iend  ,j,know))+       &
     &                               Cx*tl_zeta(Iend,j,know)
# ifdef MASKING
!^            zeta(Iend+1,j,kout)=zeta(Iend+1,j,kout)*                  &
!^   &                            GRID(ng)%rmask(Iend+1,j)
!^
              tl_zeta(Iend+1,j,kout)=tl_zeta(Iend+1,j,kout)*            &
     &                               GRID(ng)%rmask(Iend+1,j)
# endif
            END IF
          END DO
!
!  Eastern edge, implicit Chapman boundary condition.
!
        ELSE IF (tl_LBC(ieast,isFsur,ng)%Chapman_implicit) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%east(j)) THEN
              cff=dt2d*GRID(ng)%pm(Iend,j)
              cff1=SQRT(g*(GRID(ng)%h(Iend,j)+                          &
     &                     zeta(Iend,j,know)))
              tl_cff1=0.5_r8*g*(GRID(ng)%tl_h(Iend,j)+                  &
     &                          tl_zeta(Iend,j,know))/cff1
              Cx=cff*cff1
              tl_Cx=cff*tl_cff1
              cff2=1.0_r8/(1.0_r8+Cx)
              tl_cff2=-cff2*cff2*tl_Cx
!^            zeta(Iend+1,j,kout)=cff2*(zeta(Iend+1,j,know)+            &
!^   &                                  Cx*zeta(Iend,j,kout))
!^
              tl_zeta(Iend+1,j,kout)=tl_cff2*(zeta(Iend+1,j,know)+      &
     &                                        Cx*zeta(Iend,j,kout))+    &
     &                               cff2*(tl_zeta(Iend+1,j,know)+      &
     &                                     tl_Cx*zeta(Iend,j,kout)+     &
     &                                     Cx*tl_zeta(Iend,j,kout))
# ifdef MASKING
!^            zeta(Iend+1,j,kout)=zeta(Iend+1,j,kout)*                  &
!^   &                            GRID(ng)%rmask(Iend+1,j)
!^
              tl_zeta(Iend+1,j,kout)=tl_zeta(Iend+1,j,kout)*            &
     &                               GRID(ng)%rmask(Iend+1,j)
# endif
            END IF
          END DO
!
!  Eastern edge, clamped boundary condition.
!
        ELSE IF (tl_LBC(ieast,isFsur,ng)%clamped) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%east(j)) THEN
!^            zeta(Iend+1,j,kout)=BOUNDARY(ng)%zeta_east(j)
!^
# ifdef ADJUST_BOUNDARY
              IF (Lobc(ieast,isFsur,ng)) THEN
                tl_zeta(Iend+1,j,kout)=BOUNDARY(ng)%tl_zeta_east(j)
              ELSE
                tl_zeta(Iend+1,j,kout)=0.0_r8
              END IF
# else
              tl_zeta(Iend+1,j,kout)=0.0_r8
# endif
# ifdef MASKING
!^            zeta(Iend+1,j,kout)=zeta(Iend+1,j,kout)*                  &
!^   &                            GRID(ng)%rmask(Iend+1,j)
!^
              tl_zeta(Iend+1,j,kout)=tl_zeta(Iend+1,j,kout)*            &
     &                               GRID(ng)%rmask(Iend+1,j)
# endif
            END IF
          END DO
!
!  Eastern edge, gradient boundary condition.
!
        ELSE IF (tl_LBC(ieast,isFsur,ng)%gradient) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%east(j)) THEN
!^            zeta(Iend+1,j,kout)=zeta(Iend,j,kout)
!^
              tl_zeta(Iend+1,j,kout)=tl_zeta(Iend,j,kout)
# ifdef MASKING
!^            zeta(Iend+1,j,kout)=zeta(Iend+1,j,kout)*                  &
!^   &                            GRID(ng)%rmask(Iend+1,j)
!^
              tl_zeta(Iend+1,j,kout)=tl_zeta(Iend+1,j,kout)*            &
     &                               GRID(ng)%rmask(Iend+1,j)
# endif
            END IF
          END DO
!
!  Eastern edge, closed boundary condition.
!
        ELSE IF (tl_LBC(ieast,isFsur,ng)%closed) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%east(j)) THEN
!^            zeta(Iend+1,j,kout)=zeta(Iend,j,kout)
!^
              tl_zeta(Iend+1,j,kout)=tl_zeta(Iend,j,kout)
# ifdef MASKING
!^            zeta(Iend+1,j,kout)=zeta(Iend+1,j,kout)*                  &
!^   &                            GRID(ng)%rmask(Iend+1,j)
!^
              tl_zeta(Iend+1,j,kout)=tl_zeta(Iend+1,j,kout)*            &
     &                               GRID(ng)%rmask(Iend+1,j)
# endif
            END IF
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the southern edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
!
!  Southern edge, implicit upstream radiation condition.
!  (Not implemented since the TLM requires nonlinear celerity at every
!   timestep. Using gradient instead).
!
        IF (tl_LBC(isouth,isFsur,ng)%radiation) THEN
          IF (iic(ng).ne.0) THEN
# if defined CELERITY_READ && defined FORWARD_READ
            DO i=Istr,Iend+1
!^            grad(i,Jstr)=zeta(i  ,Jstr,know)-                         &
!^   &                     zeta(i-1,Jstr,know)
!^
              tl_grad(i,Jstr)=0.0_r8
            END DO
# endif
            DO i=Istr,Iend
              IF (LBC_apply(ng)%south(i)) THEN
# if defined CELERITY_READ && defined FORWARD_READ
                IF (tl_LBC(isouth,isFsur,ng)%nudging) THEN
                  IF (BOUNDARY(ng)%zeta_south_Ce(i).eq.0.0_r8) THEN
                    tau=FSobc_in(ng,isouth)
                  ELSE
                    tau=FSobc_out(ng,isouth)
                  END IF
                  tau=tau*dt2d
                END IF
#  ifdef RADIATION_2D
                Cx=BOUNDARY(ng)%zeta_south_Cx(i)
#  else
                Cx=0.0_r8
#  endif
                Ce=BOUNDARY(ng)%zeta_south_Ce(i)
                cff=BOUNDARY(ng)%zeta_south_C2(i)

!^              zeta(i,Jstr-1,kout)=(cff*zeta(i,Jstr-1,know)+           &
!^   &                               Ce *zeta(i,Jstr  ,kout)-           &
!^   &                               MAX(Cx,0.0_r8)*grad(i  ,Jstr)-     &
!^   &                               MIN(Cx,0.0_r8)*grad(i+1,Jstr))/    &
!^   &                              (cff+Ce)
!^
                tl_zeta(i,Jstr-1,kout)=(cff*tl_zeta(i,Jstr-1,know)+     &
     &                                  Ce *tl_zeta(i,Jstr  ,kout)-     &
     &                                  MAX(Cx,0.0_r8)*                 &
     &                                     tl_grad(i  ,Jstr-1)-         &
     &                                  MIN(Cx,0.0_r8)*                 &
     &                                     tl_grad(i+1,Jstr-1))/        &
     &                                 (cff+Ce)

                IF (tl_LBC(isouth,isFsur,ng)%nudging) THEN
!^                zeta(i,Jstr-1,kout)=zeta(i,Jstr-1,kout)+              &
!^   &                                tau*(BOUNDARY(ng)%zeta_south(i)-  &
!^   &                                     zeta(i,Jstr-1,know))
!^
                  tl_zeta(i,Jstr-1,kout)=tl_zeta(i,Jstr-1,kout)-        &
     &                                   tau*tl_zeta(i,Jstr-1,know)
                END IF
# else
!^              zeta(i,Jstr-1,kout)=zeta(i,Jstr,kout)         ! gradient
!^
                tl_zeta(i,Jstr-1,kout)=tl_zeta(i,Jstr,kout)
# endif
# ifdef MASKING
!^              zeta(i,Jstr-1,kout)=zeta(i,Jstr-1,kout)*                &
!^   &                              GRID(ng)%rmask(i,Jstr-1)
!^
                tl_zeta(i,Jstr-1,kout)=tl_zeta(i,Jstr-1,kout)*          &
     &                                 GRID(ng)%rmask(i,Jstr-1)
# endif
              END IF
            END DO
          END IF
!
!  Southern edge, explicit Chapman boundary condition.
!
        ELSE IF (tl_LBC(isouth,isFsur,ng)%Chapman_explicit) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%south(i)) THEN
              cff=dt2d*GRID(ng)%pn(i,Jstr)
              cff1=SQRT(g*(GRID(ng)%h(i,Jstr)+                          &
     &                     zeta(i,Jstr,know)))
              tl_cff1=0.5_r8*g*(GRID(ng)%tl_h(i,Jstr)+                  &
     &                          tl_zeta(i,Jstr,know))/cff1
              Ce=cff*cff1
              tl_Ce=cff*tl_cff1
!^            zeta(i,Jstr-1,kout)=(1.0_r8-Ce)*zeta(i,Jstr-1,know)+      &
!^   &                            Ce*zeta(i,Jstr,know)
!^
              tl_zeta(i,Jstr-1,kout)=(1.0_r8-Ce)*tl_zeta(i,Jstr-1,know)+&
     &                               tl_Ce*(zeta(i,Jstr-1,know)+        &
     &                                      zeta(i,Jstr  ,know))+       &
     &                               Ce*tl_zeta(i,Jstr,know)
# ifdef MASKING
!^            zeta(i,Jstr-1,kout)=zeta(i,Jstr-1,kout)*                  &
!^   &                            GRID(ng)%rmask(i,Jstr-1)
!^
              tl_zeta(i,Jstr-1,kout)=tl_zeta(i,Jstr-1,kout)*            &
     &                               GRID(ng)%rmask(i,Jstr-1)
# endif
            END IF
          END DO
!
!  Southern edge, implicit Chapman boundary condition.
!
        ELSE IF (tl_LBC(isouth,isFsur,ng)%Chapman_implicit) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%south(i)) THEN
              cff=dt2d*GRID(ng)%pn(i,Jstr)
              cff1=SQRT(g*(GRID(ng)%h(i,Jstr)+                          &
     &                     zeta(i,Jstr,know)))
              tl_cff1=0.5_r8*g*(GRID(ng)%tl_h(i,Jstr)+                  &
     &                          tl_zeta(i,Jstr,know))/cff1
              Ce=cff*cff1
              tl_Ce=cff*tl_cff1
              cff2=1.0_r8/(1.0_r8+Ce)
              tl_cff2=-cff2*cff2*tl_Ce
!^            zeta(i,Jstr-1,kout)=cff2*(zeta(i,Jstr-1,know)+            &
!^   &                                  Ce*zeta(i,Jstr,kout))
!^
              tl_zeta(i,Jstr-1,kout)=tl_cff2*(zeta(i,Jstr-1,know)+      &
     &                                        Ce*zeta(i,Jstr,kout))+    &
     &                               cff2*(tl_zeta(i,Jstr-1,know)+      &
     &                                     tl_Ce*zeta(i,Jstr,kout)+     &
     &                                     Ce*tl_zeta(i,Jstr,kout))
# ifdef MASKING
!^            zeta(i,Jstr-1,kout)=zeta(i,Jstr-1,kout)*                  &
!^   &                            GRID(ng)%rmask(i,Jstr-1)
!^
              tl_zeta(i,Jstr-1,kout)=tl_zeta(i,Jstr-1,kout)*            &
     &                               GRID(ng)%rmask(i,Jstr-1)
# endif
            END IF
          END DO
!
!  Southern edge, clamped boundary condition.
!
        ELSE IF (tl_LBC(isouth,isFsur,ng)%clamped) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%south(i)) THEN
!^            zeta(i,Jstr-1,kout)=BOUNDARY(ng)%zeta_south(i)
!^
# ifdef ADJUST_BOUNDARY
              IF (Lobc(isouth,isFsur,ng)) THEN
                tl_zeta(i,Jstr-1,kout)=BOUNDARY(ng)%tl_zeta_south(i)
              ELSE
                tl_zeta(i,Jstr-1,kout)=0.0_r8
              END IF
# else
              tl_zeta(i,Jstr-1,kout)=0.0_r8
# endif
# ifdef MASKING
!^            zeta(i,Jstr-1,kout)=zeta(i,Jstr-1,kout)*                  &
!^   &                            GRID(ng)%rmask(i,Jstr-1)
!^
              tl_zeta(i,Jstr-1,kout)=tl_zeta(i,Jstr-1,kout)*            &
     &                               GRID(ng)%rmask(i,Jstr-1)
# endif
            END IF
          END DO
!
!  Southern edge, gradient boundary condition.
!
        ELSE IF (tl_LBC(isouth,isFsur,ng)%gradient) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%south(i)) THEN
!^            zeta(i,Jstr-1,kout)=zeta(i,Jstr,kout)
!^
              tl_zeta(i,Jstr-1,kout)=tl_zeta(i,Jstr,kout)
# ifdef MASKING
!^            zeta(i,Jstr-1,kout)=zeta(i,Jstr-1,kout)*                  &
!^   &                            GRID(ng)%rmask(i,Jstr-1)
!^
              tl_zeta(i,Jstr-1,kout)=tl_zeta(i,Jstr-1,kout)*            &
     &                               GRID(ng)%rmask(i,Jstr-1)
# endif
            END IF
          END DO
!
!  Southern edge, closed boundary condition.
!
        ELSE IF (tl_LBC(isouth,isFsur,ng)%closed) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%south(i)) THEN
!^            zeta(i,Jstr-1,kout)=zeta(i,Jstr,kout)
!^
              tl_zeta(i,Jstr-1,kout)=tl_zeta(i,Jstr,kout)
# ifdef MASKING
!^            zeta(i,Jstr-1,kout)=zeta(i,Jstr-1,kout)*                  &
!^   &                            GRID(ng)%rmask(i,Jstr-1)
!^
              tl_zeta(i,Jstr-1,kout)=tl_zeta(i,Jstr-1,kout)*            &
     &                               GRID(ng)%rmask(i,Jstr-1)
# endif
            END IF
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the northern edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
!
!  Northern edge, implicit upstream radiation condition.
!  (Not implemented since the TLM requires nonlinear celerity at every
!   timestep. Using gradient instead).
!
        IF (tl_LBC(inorth,isFsur,ng)%radiation) THEN
          IF (iic(ng).ne.0) THEN
# if defined CELERITY_READ && defined FORWARD_READ
            DO i=Istr,Iend+1
!^            grad(i,Jend+1)=zeta(i  ,Jend+1,know)-                     &
!^   &                       zeta(i-1,Jend+1,know)
!^
              tl_grad(i,Jend+1)=0.0_r8
            END DO
# endif
            DO i=Istr,Iend
              IF (LBC_apply(ng)%north(i)) THEN
# if defined CELERITY_READ && defined FORWARD_READ
                IF (tl_LBC(inorth,isFsur,ng)%nudging) THEN
                  IF (BOUNDARY(ng)%zeta_north_Ce(i).eq.0.0_r8) THEN
                    tau=FSobc_in(ng,inorth)
                  ELSE
                    tau=FSobc_out(ng,inorth)
                  END IF
                  tau=tau*dt2d
                END IF
#  ifdef RADIATION_2D
                Cx=BOUNDARY(ng)%zeta_north_Cx(i)
#  else
                Cx=0.0_r8
#  endif
                Ce=BOUNDARY(ng)%zeta_north_Ce(i)
                cff=BOUNDARY(ng)%zeta_north_C2(i)

!^              zeta(i,Jend+1,kout)=(cff*zeta(i,Jend+1,know)+           &
!^   &                               Ce *zeta(i,Jend  ,kout)-           &
!^   &                               MAX(Cx,0.0_r8)*grad(i  ,Jend+1)-   &
!^   &                               MIN(Cx,0.0_r8)*grad(i+1,Jend+1))/  &
!^   &                              (cff+Ce)
!^
                tl_zeta(i,Jend+1,kout)=(cff*tl_zeta(i,Jend+1,know)+     &
     &                                  Ce *tl_zeta(i,Jend  ,kout)-     &
     &                                  MAX(Cx,0.0_r8)*                 &
     &                                     tl_grad(i  ,Jend+1)-         &
     &                                  MIN(Cx,0.0_r8)*                 &
     &                                     tl_grad(i+1,Jend+1))/        &
     &                                 (cff+Ce)

                IF (tl_LBC(inorth,isFsur,ng)%nudging) THEN
!^                zeta(i,Jend+1,kout)=zeta(i,Jend+1,kout)+              &
!^   &                                tau*(BOUNDARY(ng)%zeta_north(i)-  &
!^   &                                     zeta(i,Jend+1,know))
!^
                  tl_zeta(i,Jend+1,kout)=tl_zeta(i,Jend+1,kout)-        &
     &                                   tau*tl_zeta(i,Jend+1,know)
                END IF
# else
!^              zeta(i,Jend+1,kout)=zeta(i,Jend,kout)         ! gradient
!^
                tl_zeta(i,Jend+1,kout)=tl_zeta(i,Jend,kout)
# endif
# ifdef MASKING
!^              zeta(i,Jend+1,kout)=zeta(i,Jend+1,kout)*                &
!^   &                              GRID(ng)%rmask(i,Jend+1)
!^
                tl_zeta(i,Jend+1,kout)=tl_zeta(i,Jend+1,kout)*          &
     &                                 GRID(ng)%rmask(i,Jend+1)
# endif
              END IF
            END DO
          END IF
!
!  Northern edge, explicit Chapman boundary condition.
!
        ELSE IF (tl_LBC(inorth,isFsur,ng)%Chapman_explicit) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%north(i)) THEN
              cff=dt2d*GRID(ng)%pn(i,Jend)
              cff1=SQRT(g*(GRID(ng)%h(i,Jend)+                          &
     &                     zeta(i,Jend,know)))
              tl_cff1=0.5_r8*g*(GRID(ng)%tl_h(i,Jend)+                  &
     &                          tl_zeta(i,Jend,know))/cff1
              Ce=cff*cff1
              tl_Ce=cff*tl_cff1
!^            zeta(i,Jend+1,kout)=(1.0_r8-Ce)*zeta(i,Jend+1,know)+      &
!^   &                            Ce*zeta(i,Jend,know)
!^
              tl_zeta(i,Jend+1,kout)=(1.0_r8-Ce)*tl_zeta(i,Jend+1,know)+&
     &                               tl_Ce*(zeta(i,Jend+1,know)+        &
     &                                      zeta(i,Jend  ,know))+       &
     &                               Ce*tl_zeta(i,Jend,know)
# ifdef MASKING
!^            zeta(i,Jend+1,kout)=zeta(i,Jend+1,kout)*                  &
!^   &                            GRID(ng)%rmask(i,Jend+1)
!^
              tl_zeta(i,Jend+1,kout)=tl_zeta(i,Jend+1,kout)*            &
     &                               GRID(ng)%rmask(i,Jend+1)
# endif
            END IF
          END DO
!
!  Northern edge, implicit Chapman boundary condition.
!
        ELSE IF (tl_LBC(inorth,isFsur,ng)%Chapman_implicit) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%north(i)) THEN
              cff=dt2d*GRID(ng)%pn(i,Jend)
              cff1=SQRT(g*(GRID(ng)%h(i,Jend)+                          &
     &                     zeta(i,Jend,know)))
              tl_cff1=0.5_r8*g*(GRID(ng)%tl_h(i,Jend)+                  &
     &                          tl_zeta(i,Jend,know))/cff1
              Ce=cff*cff1
              tl_Ce=cff*tl_cff1
              cff2=1.0_r8/(1.0_r8+Ce)
              tl_cff2=-cff2*cff2*tl_Ce
!^            zeta(i,Jend+1,kout)=cff2*(zeta(i,Jend+1,know)+            &
!^   &                                  Ce*zeta(i,Jend,kout))
!^
              tl_zeta(i,Jend+1,kout)=tl_cff2*(zeta(i,Jend+1,know)+      &
     &                                        Ce*zeta(i,Jend,kout))+    &
     &                               cff2*(tl_zeta(i,Jend+1,know)+      &
     &                                     tl_Ce*zeta(i,Jend,kout)+     &
     &                                     Ce*tl_zeta(i,Jend,kout))
# ifdef MASKING
!^            zeta(i,Jend+1,kout)=zeta(i,Jend+1,kout)*                  &
!^   &                            GRID(ng)%rmask(i,Jend+1)
!^
              tl_zeta(i,Jend+1,kout)=tl_zeta(i,Jend+1,kout)*            &
     &                               GRID(ng)%rmask(i,Jend+1)
# endif
            END IF
          END DO
!
!  Northern edge, clamped boundary condition.
!
        ELSE IF (tl_LBC(inorth,isFsur,ng)%clamped) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%north(i)) THEN
!^            zeta(i,Jend+1,kout)=BOUNDARY(ng)%zeta_north(i)
!^
# ifdef ADJUST_BOUNDARY
              IF (Lobc(inorth,isFsur,ng)) THEN
                tl_zeta(i,Jend+1,kout)=BOUNDARY(ng)%tl_zeta_north(i)
              ELSE
                tl_zeta(i,Jend+1,kout)=0.0_r8
              END IF
# else
              tl_zeta(i,Jend+1,kout)=0.0_r8
# endif
# ifdef MASKING
!^            zeta(i,Jend+1,kout)=zeta(i,Jend+1,kout)*                  &
!^   &                            GRID(ng)%rmask(i,Jend+1)
!^
              tl_zeta(i,Jend+1,kout)=tl_zeta(i,Jend+1,kout)*            &
     &                               GRID(ng)%rmask(i,Jend+1)
# endif
            END IF
          END DO
!
!  Northern edge, gradient boundary condition.
!
        ELSE IF (tl_LBC(inorth,isFsur,ng)%gradient) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%north(i)) THEN
!^            zeta(i,Jend+1,kout)=zeta(i,Jend,kout)
!^
              tl_zeta(i,Jend+1,kout)=tl_zeta(i,Jend,kout)
# ifdef MASKING
!^            zeta(i,Jend+1,kout)=zeta(i,Jend+1,kout)*                  &
!^   &                            GRID(ng)%rmask(i,Jend+1)
!^
              tl_zeta(i,Jend+1,kout)=tl_zeta(i,Jend+1,kout)*            &
     &                               GRID(ng)%rmask(i,Jend+1)
# endif
            END IF
          END DO
!
!  Northern edge, closed boundary condition.
!
        ELSE IF (tl_LBC(inorth,isFsur,ng)%closed) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%north(i)) THEN
!^            zeta(i,Jend+1,kout)=zeta(i,Jend,kout)
!^
              tl_zeta(i,Jend+1,kout)=tl_zeta(i,Jend,kout)
# ifdef MASKING
!^            zeta(i,Jend+1,kout)=zeta(i,Jend+1,kout)*                  &
!^   &                            GRID(ng)%rmask(i,Jend+1)
!^
              tl_zeta(i,Jend+1,kout)=tl_zeta(i,Jend+1,kout)*            &
     &                               GRID(ng)%rmask(i,Jend+1)
# endif
            END IF
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (.not.(EWperiodic(ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Istr-1).and.                          &
     &        LBC_apply(ng)%west (Jstr-1)) THEN
!^          zeta(Istr-1,Jstr-1,kout)=0.5_r8*(zeta(Istr  ,Jstr-1,kout)+  &
!^   &                                       zeta(Istr-1,Jstr  ,kout))
!^
            tl_zeta(Istr-1,Jstr-1,kout)=0.5_r8*                         &
     &                                  (tl_zeta(Istr  ,Jstr-1,kout)+   &
     &                                   tl_zeta(Istr-1,Jstr  ,kout))
          END IF
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jstr-1)) THEN
!^          zeta(Iend+1,Jstr-1,kout)=0.5_r8*(zeta(Iend  ,Jstr-1,kout)+  &
!^   &                                       zeta(Iend+1,Jstr  ,kout))
!^
            tl_zeta(Iend+1,Jstr-1,kout)=0.5_r8*                         &
     &                                  (tl_zeta(Iend  ,Jstr-1,kout)+   &
     &                                   tl_zeta(Iend+1,Jstr  ,kout))
          END IF
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Istr-1).and.                          &
     &        LBC_apply(ng)%west (Jend+1)) THEN
!^          zeta(Istr-1,Jend+1,kout)=0.5_r8*(zeta(Istr-1,Jend  ,kout)+  &
!^   &                                       zeta(Istr  ,Jend+1,kout))
!^
            tl_zeta(Istr-1,Jend+1,kout)=0.5_r8*                         &
     &                                  (tl_zeta(Istr-1,Jend  ,kout)+   &
     &                                   tl_zeta(Istr  ,Jend+1,kout))
          END IF
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jend+1)) THEN
!^          zeta(Iend+1,Jend+1,kout)=0.5_r8*(zeta(Iend+1,Jend  ,kout)+  &
!^   &                                       zeta(Iend  ,Jend+1,kout))
!^
            tl_zeta(Iend+1,Jend+1,kout)=0.5_r8*                         &
     &                                  (tl_zeta(Iend+1,Jend  ,kout)+   &
     &                                   tl_zeta(Iend  ,Jend+1,kout))
          END IF
        END IF
      END IF

# if defined WET_DRY
!
!-----------------------------------------------------------------------
! Ensure that water level on boundary cells is above bed elevation.
!-----------------------------------------------------------------------
!
      cff=Dcrit(ng)-eps
      IF (.not.EWperiodic(ng)) THEN
        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%west(j)) THEN
              IF (zeta(Istr-1,j,kout).le.                               &
     &            (Dcrit(ng)-GRID(ng)%h(Istr-1,j))) THEN
!^              zeta(Istr-1,j,kout)=cff-GRID(ng)%h(Istr-1,j)
!^
                tl_zeta(Istr-1,j,kout)=-GRID(ng)%tl_h(Istr-1,j)
              END IF
            END IF
          END DO
        END IF
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%east(j)) THEN
              IF (zeta(Iend+1,j,kout).le.                               &
     &            (Dcrit(ng)-GRID(ng)%h(Iend+1,j))) THEN
!^              zeta(Iend+1,j,kout)=cff-GRID(ng)%h(Iend+1,j)
!^
                tl_zeta(Iend+1,j,kout)=-GRID(ng)%tl_h(Iend+1,j)
              END IF
            END IF
          END DO
        END IF
      END IF
!
      IF (.not.NSperiodic(ng)) THEN
        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%south(i)) THEN
              IF (zeta(i,Jstr-1,kout).le.                               &
     &            (Dcrit(ng)-GRID(ng)%h(i,Jstr-1))) THEN
!^              zeta(i,Jstr-1,kout)=cff-GRID(ng)%h(i,Jstr-1)
!^
                tl_zeta(i,Jstr-1,kout)=-GRID(ng)%tl_h(i,Jstr-1)
              END IF
            END IF
          END DO
        END IF
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%north(i)) THEN
              IF (zeta(i,Jend+1,kout).le.                               &
     &            (Dcrit(ng)-GRID(ng)%h(i,Jend+1))) THEN
!^              zeta(i,Jend+1,kout)=cff-GRID(ng)%h(i,Jend+1)
!^
                tl_zeta(i,Jend+1,kout)=-GRID(ng)%tl_h(i,Jend+1)
              END IF
            END IF
          END DO
        END IF
      END IF
!
      IF (.not.(EWperiodic(ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Istr-1).and.                          &
     &        LBC_apply(ng)%west (Jstr-1)) THEN
            IF (zeta(Istr-1,Jstr-1,kout).le.                            &
     &          (Dcrit(ng)-GRID(ng)%h(Istr-1,Jstr-1))) THEN
!^            zeta(Istr-1,Jstr-1,kout)=cff-GRID(ng)%h(Istr-1,Jstr-1)
!^
              tl_zeta(Istr-1,Jstr-1,kout)=-GRID(ng)%tl_h(Istr-1,Jstr-1)
            END IF
          END IF
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jstr-1)) THEN
            IF (zeta(Iend+1,Jstr-1,kout).le.                            &
     &          (Dcrit(ng)-GRID(ng)%h(Iend+1,Jstr-1))) THEN
!^            zeta(Iend+1,Jstr-1,kout)=cff-GRID(ng)%h(Iend+1,Jstr-1)
!^
              tl_zeta(Iend+1,Jstr-1,kout)=-GRID(ng)%tl_h(Iend+1,Jstr-1)
            END IF
          END IF
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Istr-1).and.                          &
     &        LBC_apply(ng)%west (Jend+1)) THEN
            IF (zeta(Istr-1,Jend+1,kout).le.                            &
     &          (Dcrit(ng)-GRID(ng)%h(Istr-1,Jend+1))) THEN
!^            zeta(Istr-1,Jend+1,kout)=cff-GRID(ng)%h(Istr-1,Jend+1)
!^
              tl_zeta(Istr-1,Jend+1,kout)=-GRID(ng)%tl_h(Istr-1,Jend+1)
            END IF
          END IF
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jend+1)) THEN
            IF (zeta(Iend+1,Jend+1,kout).le.                            &
     &          (Dcrit(ng)-GRID(ng)%h(Iend+1,Jend+1))) THEN
!^            zeta(Iend+1,Jend+1,kout)=cff-GRID(ng)%h(Iend+1,Jend+1)
!^
              tl_zeta(Iend+1,Jend+1,kout)=-GRID(ng)%tl_h(Iend+1,Jend+1)
            END IF
          END IF
        END IF
      END IF
# endif
!
      RETURN
      END SUBROUTINE tl_zetabc_tile

# if defined STEP2D_FB_AB3_AM4 || defined STEP2D_FB_LF_AM3
!
!***********************************************************************
      SUBROUTINE tl_zetabc_local (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            kstp,                                 &
     &                            zeta, tl_zeta,                        &
     &                            zeta_new, tl_zeta_new)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: kstp
!
#  ifdef ASSUMED_SHAPE
      real(r8), intent(in   ) :: zeta(LBi:,LBj:,:)
      real(r8), intent(in   ) :: tl_zeta(LBi:,LBj:,:)
      real(r8), intent(in   ) :: zeta_new(IminS:,JminS:)

      real(r8), intent(inout) :: tl_zeta_new(IminS:,JminS:)
#  else
      real(r8), intent(in   ) :: zeta(LBi:UBi,LBj:UBj,:)
      real(r8), intent(in   ) :: tl_zeta(LBi:UBi,LBj:UBj,:)
      real(r8), intent(in   ) :: zeta_new(IminS:ImaxS,JminS:JmaxS)

      real(r8), intent(inout) :: tl_zeta_new(IminS:ImaxS,JminS:JmaxS)
#  endif
!
!  Local variable declarations.
!
      integer :: i, j, know

      real(r8) :: Ce, Cx
      real(r8) :: cff, cff1, cff2, dt2d, tau

      real(r8) :: tl_Ce, tl_Cx
      real(r8) :: tl_cff1, tl_cff2

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_grad

#  include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set time-indices
!-----------------------------------------------------------------------
!
      know=kstp
      dt2d=dtfast(ng)
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the western edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
!
!  Western edge, implicit upstream radiation condition.
!  (Not implemented since the TLM requires nonlinear celerity at every
!   timestep. Using gradient instead).
!
        IF (tl_LBC(iwest,isFsur,ng)%radiation) THEN
          IF (iic(ng).ne.0) THEN
#  if defined CELERITY_READ && defined FORWARD_READ
            DO j=JstrV-1,Jend+1
!^            grad(Istr-1,j)=zeta(Istr-1,j  ,know)-                     &
!^   &                       zeta(Istr-1,j-1,know)
!^
              tl_grad(Istr-1,j)=0.0_r8
            END DO
#  endif
            DO j=JstrV-1,Jend
              IF (LBC_apply(ng)%west(j)) THEN
#  if defined CELERITY_READ && defined FORWARD_READ
                IF (tl_LBC(iwest,isFsur,ng)%nudging) THEN
                  IF (BOUNDARY(ng)%zeta_west_Cx(j).eq.0.0_r8) THEN
                    tau=FSobc_in(ng,iwest)
                  ELSE
                    tau=FSobc_out(ng,iwest)
                  END IF
                  tau=tau*dt2d
                END IF
                Cx=BOUNDARY(ng)%zeta_west_Cx(j)
#   ifdef RADIATION_2D
                Ce=BOUNDARY(ng)%zeta_west_Ce(j)
#   else
                Ce=0.0_r8
#   endif
                cff=BOUNDARY(ng)%zeta_west_C2(j)

!^              zeta_new(Istr-1,j)=(cff*zeta(Istr-1,j,know)+            &
!^   &                              Cx *zeta_new(Istr,j)-               &
!^   &                              MAX(Ce,0.0_r8)*grad(Istr-1,j  )-    &
!^   &                              MIN(Ce,0.0_r8)*grad(Istr-1,j+1))/   &
!^   &                             (cff+Cx)
!^
                tl_zeta_new(Istr-1,j)=(cff*tl_zeta(Istr-1,j,know)+      &
     &                                 Cx *tl_zeta_new(Istr,j)-         &
     &                                 MAX(Ce,0.0_r8)*                  &
     &                                    tl_grad(Istr-1,j  )-          &
     &                                 MIN(Ce,0.0_r8)*                  &
     &                                    tl_grad(Istr-1,j+1))/         &
     &                                (cff+Cx)

                IF (tl_LBC(iwest,isFsur,ng)%nudging) THEN
!^                zeta_new(Istr-1,j)=zeta_new(Istr-1,j)+                &
!^   &                               tau*(BOUNDARY(ng)%zeta_west(j)-    &
!^   &                                    zeta(Istr-1,j,know))
!^
                  tl_zeta_new(Istr-1,j)=tl_zeta_new(Istr-1,j)-          &
     &                                  tau*tl_zeta(Istr-1,j,know)
                END IF
#  else
!^              zeta_new(Istr-1,j)=zeta_new(Istr,j)           ! gradient
!^
                tl_zeta_new(Istr-1,j)=tl_zeta_new(Istr,j)
#  endif
#  ifdef MASKING
!^              zeta_new(Istr-1,j)=zeta_new(Istr-1,j)*                  &
!^   &                             GRID(ng)%rmask(Istr-1,j)
!^
                tl_zeta_new(Istr-1,j)=tl_zeta_new(Istr-1,j)*            &
     &                                GRID(ng)%rmask(Istr-1,j)
#  endif
              END IF
            END DO
          END IF
!
!  Western edge, explicit Chapman boundary condition.
!
        ELSE IF (tl_LBC(iwest,isFsur,ng)%Chapman_explicit) THEN
          DO j=JstrV-1,Jend
            IF (LBC_apply(ng)%west(j)) THEN
              cff=dt2d*GRID(ng)%pm(Istr,j)
              cff1=SQRT(g*(GRID(ng)%h(Istr,j)+                          &
     &                     zeta(Istr,j,know)))
              tl_cff1=0.5_r8*g*(GRID(ng)%tl_h(Istr,j)+                  &
     &                          tl_zeta(Istr,j,know))/cff1
              Cx=cff*cff1
              tl_Cx=cff*tl_cff1
!^            zeta_new(Istr-1,j)=(1.0_r8-Cx)*zeta(Istr-1,j,know)+       &
!^   &                           Cx*zeta(Istr,j,know)
!^
              tl_zeta_new(Istr-1,j)=(1.0_r8-Cx)*tl_zeta(Istr-1,j,know)+ &
     &                              tl_Cx*(zeta(Istr-1,j,know)+         &
     &                                     zeta(Istr  ,j,know))+        &
     &                              Cx*tl_zeta(Istr,j,know)
#  ifdef MASKING
!^            zeta_new(Istr-1,j)=zeta_new(Istr-1,j)*                    &
!^   &                           GRID(ng)%rmask(Istr-1,j)
!^
              tl_zeta_new(Istr-1,j)=tl_zeta_new(Istr-1,j)*              &
     &                              GRID(ng)%rmask(Istr-1,j)
#  endif
            END IF
          END DO
!
!  Western edge, implicit Chapman boundary condition.
!
        ELSE IF (tl_LBC(iwest,isFsur,ng)%Chapman_implicit) THEN
          DO j=JstrV-1,Jend
            IF (LBC_apply(ng)%west(j)) THEN
              cff=dt2d*GRID(ng)%pm(Istr,j)
              cff1=SQRT(g*(GRID(ng)%h(Istr,j)+                          &
     &                     zeta(Istr,j,know)))
              tl_cff1=0.5_r8*g*(GRID(ng)%tl_h(Istr,j)+                  &
     &                          tl_zeta(Istr,j,know))/cff1
              Cx=cff*cff1
              tl_Cx=cff*tl_cff1
              cff2=1.0_r8/(1.0_r8+Cx)
              tl_cff2=-cff2*cff2*tl_Cx
!^            zeta_new(Istr-1,j)=cff2*(zeta(Istr-1,j,know)+             &
!^   &                                 Cx*zeta_new(Istr,j))
!^
              tl_zeta_new(Istr-1,j)=tl_cff2*(zeta(Istr-1,j,know)+       &
     &                                       Cx*zeta_new(Istr,j))+      &
     &                              cff2*(tl_zeta(Istr-1,j,know)+       &
     &                                    tl_Cx*zeta_new(Istr,j)+       &
     &                                    Cx*tl_zeta_new(Istr,j))
#  ifdef MASKING
!^            zeta_new(Istr-1,j)=zeta_new(Istr-1,j)*                    &
!^   &                           GRID(ng)%rmask(Istr-1,j)
!^
              tl_zeta_new(Istr-1,j)=tl_zeta_new(Istr-1,j)*              &
     &                              GRID(ng)%rmask(Istr-1,j)
#  endif
            END IF
          END DO
!
!  Western edge, clamped boundary condition.
!
        ELSE IF (tl_LBC(iwest,isFsur,ng)%clamped) THEN
          DO j=JstrV-1,Jend
            IF (LBC_apply(ng)%west(j)) THEN
!^            zeta_new(Istr-1,j)=BOUNDARY(ng)%zeta_west(j)
!^
#  ifdef ADJUST_BOUNDARY
              IF (Lobc(iwest,isFsur,ng)) THEN
                tl_zeta_new(Istr-1,j)=BOUNDARY(ng)%tl_zeta_west(j)
              ELSE
                tl_zeta_new(Istr-1,j)=0.0_r8
              END IF
#  else
              tl_zeta_new(Istr-1,j)=0.0_r8
#  endif
#  ifdef MASKING
!^            zeta_new(Istr-1,j)=zeta_new(Istr-1,j)*                    &
!^   &                           GRID(ng)%rmask(Istr-1,j)
!^
              tl_zeta_new(Istr-1,j)=tl_zeta_new(Istr-1,j)*              &
     &                              GRID(ng)%rmask(Istr-1,j)
#  endif
            END IF
          END DO
!
!  Western edge, gradient boundary condition.
!
        ELSE IF (tl_LBC(iwest,isFsur,ng)%gradient) THEN
          DO j=JstrV-1,Jend
            IF (LBC_apply(ng)%west(j)) THEN
!^            zeta_new(Istr-1,j)=zeta_new(Istr,j)
!^
              tl_zeta_new(Istr-1,j)=tl_zeta_new(Istr,j)
#  ifdef MASKING
!^            zeta_new(Istr-1,j)=zeta_new(Istr-1,j)*                    &
!^   &                           GRID(ng)%rmask(Istr-1,j)
!^
              tl_zeta_new(Istr-1,j)=tl_zeta_new(Istr-1,j)*              &
     &                              GRID(ng)%rmask(Istr-1,j)
#  endif
            END IF
          END DO
!
!  Western edge, closed boundary condition.
!
        ELSE IF (tl_LBC(iwest,isFsur,ng)%closed) THEN
          DO j=JstrV-1,Jend
            IF (LBC_apply(ng)%west(j)) THEN
!^            zeta_new(Istr-1,j)=zeta_new(Istr,j)
!^
              tl_zeta_new(Istr-1,j)=tl_zeta_new(Istr,j)
#  ifdef MASKING
!^            zeta_new(Istr-1,j)=zeta_new(Istr-1,j)*                    &
!^   &                           GRID(ng)%rmask(Istr-1,j)
!^
              tl_zeta_new(Istr-1,j)=tl_zeta_new(Istr-1,j)*              &
     &                              GRID(ng)%rmask(Istr-1,j)
#  endif
            END IF
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the eastern edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
!
!  Eastern edge, implicit upstream radiation condition.
!  (Not implemented since the TLM requires nonlinear celerity at every
!   timestep. Using gradient instead).
!
        IF (tl_LBC(ieast,isFsur,ng)%radiation) THEN
          IF (iic(ng).ne.0) THEN
#  if defined CELERITY_READ && defined FORWARD_READ
            DO j=JstrV-1,Jend+1
!^            grad(Iend+1,j)=zeta(Iend+1,j  ,know)-                     &
!^   &                       zeta(Iend+1,j-1,know)
!^
              tl_grad(Iend+1,j)=0.0_r8
            END DO
#  endif
            DO j=JstrV-1,Jend
              IF (LBC_apply(ng)%east(j)) THEN
#  if defined CELERITY_READ && defined FORWARD_READ
                IF (tl_LBC(ieast,isFsur,ng)%nudging) THEN
                  IF (BOUNDARY(ng)%zeta_east_Cx(j).eq.0.0_r8) THEN
                    tau=FSobc_in(ieast)
                  ELSE
                    tau=FSobc_out(ieast)
                  END IF
                  tau=tau*dt2d
                END IF
                Cx=BOUNDARY(ng)%zeta_east_Cx(j)
#   ifdef RADIATION_2D
                Ce=BOUNDARY(ng)%zeta_east_Ce(j)
#   else
                Ce=0.0_r8
#   endif
                cff=BOUNDARY(ng)%zeta_east_C2(j)

!^              zeta_new(Iend+1,j)=(cff*zeta(Iend+1,j,know)+            &
!^   &                              Cx *zeta_new(Iend,j)-               &
!^   &                              MAX(Ce,0.0_r8)*grad(Iend+1,j  )-    &
!^   &                              MIN(Ce,0.0_r8)*grad(Iend+1,j+1))/   &
!^   &                             (cff+Cx)
!^
                tl_zeta_new(Iend+1,j)=(cff*tl_zeta(Iend+1,j,know)+      &
     &                                 Cx *tl_zeta_new(Iend,j)-         &
     &                                 MAX(Ce,0.0_r8)*                  &
     &                                    tl_grad(Iend+1,j  )-          &
     &                                 MIN(Ce,0.0_r8)*                  &
     &                                    tl_grad(Iend+1,j+1))/         &
     &                                (cff+Cx)

                IF (tl_LBC(ieast,isFsur,ng)%nudging) THEN
!^                zeta_new(Iend+1,j)=zeta_new(Iend+1,j)+                &
!^   &                               tau*(BOUNDARY(ng)%zeta_east(j)-    &
!^   &                                    zeta(Iend+1,j,know))
!^
                  tl_zeta_new(Iend+1,j)=tl_zeta_new(Iend+1,j)-          &
     &                                  tau*tl_zeta(Iend+1,j,know)
                END IF
#  else
!^              zeta_new(Iend+1,j)=zeta_new(Iend,j)           ! gradient
!^
                tl_zeta_new(Iend+1,j)=tl_zeta_new(Iend,j)
#  endif
#  ifdef MASKING
!^              zeta_new(Iend+1,j)=zeta_new(Iend+1,j)*                  &
!^   &                             GRID(ng)%rmask(Iend+1,j)
!^
                tl_zeta_new(Iend+1,j)=tl_zeta_new(Iend+1,j)*            &
     &                                GRID(ng)%rmask(Iend+1,j)
#  endif
              END IF
            END DO
          END IF
!
!  Eastern edge, explicit Chapman boundary condition.
!
        ELSE IF (tl_LBC(ieast,isFsur,ng)%Chapman_explicit) THEN
          DO j=JstrV-1,Jend
            IF (LBC_apply(ng)%east(j)) THEN
              cff=dt2d*GRID(ng)%pm(Iend,j)
              cff1=SQRT(g*(GRID(ng)%h(Iend,j)+                          &
     &                     zeta(Iend,j,know)))
              tl_cff1=0.5_r8*g*(GRID(ng)%tl_h(Iend,j)+                  &
     &                          tl_zeta(Iend,j,know))/cff1
              Cx=cff*cff1
              tl_Cx=cff*tl_cff1
!^            zeta_new(Iend+1,j)=(1.0_r8-Cx)*zeta(Iend+1,j,know)+       &
!^   &                           Cx*zeta(Iend,j,know)
!^
              tl_zeta_new(Iend+1,j)=(1.0_r8-Cx)*tl_zeta(Iend+1,j,know)+ &
     &                              tl_Cx*(zeta(Iend+1,j,know)+         &
     &                                     zeta(Iend  ,j,know))+        &
     &                              Cx*tl_zeta(Iend,j,know)
#  ifdef MASKING
!^            zeta_new(Iend+1,j)=zeta_new(Iend+1,j)*                    &
!^   &                           GRID(ng)%rmask(Iend+1,j)
!^
              tl_zeta_new(Iend+1,j)=tl_zeta_new(Iend+1,j)*              &
     &                              GRID(ng)%rmask(Iend+1,j)
#  endif
            END IF
          END DO
!
!  Eastern edge, implicit Chapman boundary condition.
!
        ELSE IF (tl_LBC(ieast,isFsur,ng)%Chapman_implicit) THEN
          DO j=JstrV-1,Jend
            IF (LBC_apply(ng)%east(j)) THEN
              cff=dt2d*GRID(ng)%pm(Iend,j)
              cff1=SQRT(g*(GRID(ng)%h(Iend,j)+                          &
     &                     zeta(Iend,j,know)))
              tl_cff1=0.5_r8*g*(GRID(ng)%tl_h(Iend,j)+                  &
     &                          tl_zeta(Iend,j,know))/cff1
              Cx=cff*cff1
              tl_Cx=cff*tl_cff1
              cff2=1.0_r8/(1.0_r8+Cx)
              tl_cff2=-cff2*cff2*tl_Cx
!^            zeta_new(Iend+1,j)=cff2*(zeta(Iend+1,j,know)+             &
!^   &                                 Cx*zeta_new(Iend,j))
!^
              tl_zeta_new(Iend+1,j)=tl_cff2*(zeta(Iend+1,j,know)+       &
     &                                       Cx*zeta_new(Iend,j))+      &
     &                              cff2*(tl_zeta(Iend+1,j,know)+       &
     &                                    tl_Cx*zeta_new(Iend,j)+       &
     &                                    Cx*tl_zeta_new(Iend,j))
#  ifdef MASKING
!^            zeta_new(Iend+1,j)=zeta_new(Iend+1,j)*                    &
!^   &                           GRID(ng)%rmask(Iend+1,j)
!^
              tl_zeta_new(Iend+1,j)=tl_zeta_new(Iend+1,j)*              &
     &                              GRID(ng)%rmask(Iend+1,j)
#  endif
            END IF
          END DO
!
!  Eastern edge, clamped boundary condition.
!
        ELSE IF (tl_LBC(ieast,isFsur,ng)%clamped) THEN
          DO j=JstrV-1,Jend
            IF (LBC_apply(ng)%east(j)) THEN
!^            zeta_new(Iend+1,j)=BOUNDARY(ng)%zeta_east(j)
!^
#  ifdef ADJUST_BOUNDARY
              IF (Lobc(ieast,isFsur,ng)) THEN
                tl_zeta_new(Iend+1,j)=BOUNDARY(ng)%tl_zeta_east(j)
              ELSE
                tl_zeta_new(Iend+1,j)=0.0_r8
              END IF
#  else
              tl_zeta_new(Iend+1,j)=0.0_r8
#  endif
#  ifdef MASKING
!^            zeta_new(Iend+1,j)=zeta_new(Iend+1,j)*                    &
!^   &                           GRID(ng)%rmask(Iend+1,j)
!^
              tl_zeta_new(Iend+1,j)=tl_zeta_new(Iend+1,j)*              &
     &                              GRID(ng)%rmask(Iend+1,j)
#  endif
            END IF
          END DO
!
!  Eastern edge, gradient boundary condition.
!
        ELSE IF (tl_LBC(ieast,isFsur,ng)%gradient) THEN
          DO j=JstrV-1,Jend
            IF (LBC_apply(ng)%east(j)) THEN
!^            zeta_new(Iend+1,j)=zeta_new(Iend,j)
!^
              tl_zeta_new(Iend+1,j)=tl_zeta_new(Iend,j)
#  ifdef MASKING
!^            zeta_new(Iend+1,j)=zeta_new(Iend+1,j)*                    &
!^   &                           GRID(ng)%rmask(Iend+1,j)
!^
              tl_zeta_new(Iend+1,j)=tl_zeta_new(Iend+1,j)*              &
     &                              GRID(ng)%rmask(Iend+1,j)
#  endif
            END IF
          END DO
!
!  Eastern edge, closed boundary condition.
!
        ELSE IF (tl_LBC(ieast,isFsur,ng)%closed) THEN
          DO j=JstrV-1,Jend
            IF (LBC_apply(ng)%east(j)) THEN
!^            zeta_new(Iend+1,j)=zeta_new(Iend,j)
!^
              tl_zeta_new(Iend+1,j)=tl_zeta_new(Iend,j)
#  ifdef MASKING
!^            zeta_new(Iend+1,j)=zeta_new(Iend+1,j)*                    &
!^   &                           GRID(ng)%rmask(Iend+1,j)
!^
              tl_zeta_new(Iend+1,j)=tl_zeta_new(Iend+1,j)*              &
     &                              GRID(ng)%rmask(Iend+1,j)
#  endif
            END IF
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the southern edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
!
!  Southern edge, implicit upstream radiation condition.
!  (Not implemented since the TLM requires nonlinear celerity at every
!   timestep. Using gradient instead).
!
        IF (tl_LBC(isouth,isFsur,ng)%radiation) THEN
          IF (iic(ng).ne.0) THEN
#  if defined CELERITY_READ && defined FORWARD_READ
            DO i=IstrU-1,Iend+1
!^            grad(i,Jstr)=zeta(i  ,Jstr,know)-                         &
!^   &                     zeta(i-1,Jstr,know)
!^
              tl_grad(i,Jstr)=0.0_r8
            END DO
#  endif
            DO i=IstrU-1,Iend
              IF (LBC_apply(ng)%south(i)) THEN
#  if defined CELERITY_READ && defined FORWARD_READ
                IF (tl_LBC(isouth,isFsur,ng)%nudging) THEN
                  IF (BOUNDARY(ng)%zeta_south_Ce(i).eq.0.0_r8) THEN
                    tau=FSobc_in(ng,isouth)
                  ELSE
                    tau=FSobc_out(ng,isouth)
                  END IF
                  tau=tau*dt2d
                END IF
#   ifdef RADIATION_2D
                Cx=BOUNDARY(ng)%zeta_south_Cx(i)
#   else
                Cx=0.0_r8
#   endif
                Ce=BOUNDARY(ng)%zeta_south_Ce(i)
                cff=BOUNDARY(ng)%zeta_south_C2(i)

!^              zeta_new(i,Jstr-1)=(cff*zeta(i,Jstr-1,know)+            &
!^   &                              Ce *zeta_new(i,Jstr)-               &
!^   &                              MAX(Cx,0.0_r8)*grad(i  ,Jstr)-      &
!^   &                              MIN(Cx,0.0_r8)*grad(i+1,Jstr))/     &
!^   &                             (cff+Ce)
!^
                tl_zeta_new(i,Jstr-1)=(cff*tl_zeta(i,Jstr-1,know)+      &
     &                                 Ce *tl_zeta_new(i,Jstr)-         &
     &                                 MAX(Cx,0.0_r8)*                  &
     &                                    tl_grad(i  ,Jstr-1)-          &
     &                                 MIN(Cx,0.0_r8)*                  &
     &                                    tl_grad(i+1,Jstr-1))/         &
     &                                (cff+Ce)

                IF (tl_LBC(isouth,isFsur,ng)%nudging) THEN
!^                zeta_new(i,Jstr-1)=zeta_new(i,Jstr-1)+                &
!^   &                               tau*(BOUNDARY(ng)%zeta_south(i)-   &
!^   &                                    zeta(i,Jstr-1,know))
!^
                  tl_zeta_new(i,Jstr-1)=tl_zeta_new(i,Jstr-1)-          &
     &                                  tau*tl_zeta(i,Jstr-1,know)
                END IF
#  else
!^              zeta_new(i,Jstr-1)=zeta_new(i,Jstr)           ! gradient
!^
                tl_zeta_new(i,Jstr-1)=tl_zeta_new(i,Jstr)
#  endif
#  ifdef MASKING
!^              zeta_new(i,Jstr-1)=zeta_new(i,Jstr-1)*                  &
!^   &                             GRID(ng)%rmask(i,Jstr-1)
!^
                tl_zeta_new(i,Jstr-1)=tl_zeta_new(i,Jstr-1)*            &
     &                                GRID(ng)%rmask(i,Jstr-1)
#  endif
              END IF
            END DO
          END IF
!
!  Southern edge, explicit Chapman boundary condition.
!
        ELSE IF (tl_LBC(isouth,isFsur,ng)%Chapman_explicit) THEN
          DO i=IstrU-1,Iend
            IF (LBC_apply(ng)%south(i)) THEN
              cff=dt2d*GRID(ng)%pn(i,Jstr)
              cff1=SQRT(g*(GRID(ng)%h(i,Jstr)+                          &
     &                     zeta(i,Jstr,know)))
              tl_cff1=0.5_r8*g*(GRID(ng)%tl_h(i,Jstr)+                  &
     &                          tl_zeta(i,Jstr,know))/cff1
              Ce=cff*cff1
              tl_Ce=cff*tl_cff1
!^            zeta_new(i,Jstr-1)=(1.0_r8-Ce)*zeta(i,Jstr-1,know)+       &
!^   &                           Ce*zeta(i,Jstr,know)
!^
              tl_zeta_new(i,Jstr-1)=(1.0_r8-Ce)*tl_zeta(i,Jstr-1,know)+ &
     &                              tl_Ce*(zeta(i,Jstr-1,know)+         &
     &                                     zeta(i,Jstr  ,know))+        &
     &                              Ce*tl_zeta(i,Jstr,know)
#  ifdef MASKING
!^            zeta_new(i,Jstr-1)=zeta_new(i,Jstr-1)*                    &
!^   &                           GRID(ng)%rmask(i,Jstr-1)
!^
              tl_zeta_new(i,Jstr-1)=tl_zeta_new(i,Jstr-1)*              &
     &                              GRID(ng)%rmask(i,Jstr-1)
#  endif
            END IF
          END DO
!
!  Southern edge, implicit Chapman boundary condition.
!
        ELSE IF (tl_LBC(isouth,isFsur,ng)%Chapman_implicit) THEN
          DO i=IstrU-1,Iend
            IF (LBC_apply(ng)%south(i)) THEN
              cff=dt2d*GRID(ng)%pn(i,Jstr)
              cff1=SQRT(g*(GRID(ng)%h(i,Jstr)+                          &
     &                     zeta(i,Jstr,know)))
              tl_cff1=0.5_r8*g*(GRID(ng)%tl_h(i,Jstr)+                  &
     &                          tl_zeta(i,Jstr,know))/cff1
              Ce=cff*cff1
              tl_Ce=cff*tl_cff1
              cff2=1.0_r8/(1.0_r8+Ce)
              tl_cff2=-cff2*cff2*tl_Ce
!^            zeta_new(i,Jstr-1)=cff2*(zeta(i,Jstr-1,know)+             &
!^   &                                 Ce*zeta_new(i,Jstr))
!^
              tl_zeta_new(i,Jstr-1)=tl_cff2*(zeta(i,Jstr-1,know)+       &
     &                                       Ce*zeta_new(i,Jstr))+      &
     &                              cff2*(tl_zeta(i,Jstr-1,know)+       &
     &                                    tl_Ce*zeta_new(i,Jstr)+       &
     &                                    Ce*tl_zeta_new(i,Jstr))
#  ifdef MASKING
!^            zeta_new(i,Jstr-1)=zeta_new(i,Jstr-1)*                    &
!^   &                           GRID(ng)%rmask(i,Jstr-1)
!^
              tl_zeta_new(i,Jstr-1)=tl_zeta_new(i,Jstr-1)*              &
     &                              GRID(ng)%rmask(i,Jstr-1)
#  endif
            END IF
          END DO
!
!  Southern edge, clamped boundary condition.
!
        ELSE IF (tl_LBC(isouth,isFsur,ng)%clamped) THEN
          DO i=IstrU-1,Iend
            IF (LBC_apply(ng)%south(i)) THEN
!^            zeta_new(i,Jstr-1)=BOUNDARY(ng)%zeta_south(i)
!^
#  ifdef ADJUST_BOUNDARY
              IF (Lobc(isouth,isFsur,ng)) THEN
                tl_zeta_new(i,Jstr-1)=BOUNDARY(ng)%tl_zeta_south(i)
              ELSE
                tl_zeta_new(i,Jstr-1)=0.0_r8
              END IF
#  else
              tl_zeta_new(i,Jstr-1)=0.0_r8
#  endif
#  ifdef MASKING
!^            zeta_new(i,Jstr-1)=zeta_new(i,Jstr-1)*                    &
!^   &                           GRID(ng)%rmask(i,Jstr-1)
!^
              tl_zeta_new(i,Jstr-1)=tl_zeta_new(i,Jstr-1)*              &
     &                              GRID(ng)%rmask(i,Jstr-1)
#  endif
            END IF
          END DO
!
!  Southern edge, gradient boundary condition.
!
        ELSE IF (tl_LBC(isouth,isFsur,ng)%gradient) THEN
          DO i=IstrU-1,Iend
            IF (LBC_apply(ng)%south(i)) THEN
!^            zeta_new(i,Jstr-1)=zeta_new(i,Jstr)
!^
              tl_zeta_new(i,Jstr-1)=tl_zeta_new(i,Jstr)
#  ifdef MASKING
!^            zeta_new(i,Jstr-1)=zeta_new(i,Jstr-1)*                    &
!^   &                           GRID(ng)%rmask(i,Jstr-1)
!^
              tl_zeta_new(i,Jstr-1)=tl_zeta_new(i,Jstr-1)*              &
     &                              GRID(ng)%rmask(i,Jstr-1)
#  endif
            END IF
          END DO
!
!  Southern edge, closed boundary condition.
!
        ELSE IF (tl_LBC(isouth,isFsur,ng)%closed) THEN
          DO i=IstrU-1,Iend
            IF (LBC_apply(ng)%south(i)) THEN
!^            zeta_new(i,Jstr-1)=zeta_new(i,Jstr)
!^
              tl_zeta_new(i,Jstr-1)=tl_zeta_new(i,Jstr)
#  ifdef MASKING
!^            zeta_new(i,Jstr-1)=zeta_new(i,Jstr-1)*                    &
!^   &                           GRID(ng)%rmask(i,Jstr-1)
!^
              tl_zeta_new(i,Jstr-1)=tl_zeta_new(i,Jstr-1)*              &
     &                              GRID(ng)%rmask(i,Jstr-1)
#  endif
            END IF
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the northern edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
!
!  Northern edge, implicit upstream radiation condition.
!  (Not implemented since the TLM requires nonlinear celerity at every
!   timestep. Using gradient instead).
!
        IF (tl_LBC(inorth,isFsur,ng)%radiation) THEN
          IF (iic(ng).ne.0) THEN
#  if defined CELERITY_READ && defined FORWARD_READ
            DO i=IstrU-1,Iend+1
!^            grad(i,Jend+1)=zeta(i  ,Jend+1,know)-                     &
!^   &                       zeta(i-1,Jend+1,know)
!^
              tl_grad(i,Jend+1)=0.0_r8
            END DO
#  endif
            DO i=IstrU-1,Iend
              IF (LBC_apply(ng)%north(i)) THEN
#  if defined CELERITY_READ && defined FORWARD_READ
                IF (tl_LBC(inorth,isFsur,ng)%nudging) THEN
                  IF (BOUNDARY(ng)%zeta_north_Ce(i).eq.0.0_r8) THEN
                    tau=FSobc_in(ng,inorth)
                  ELSE
                    tau=FSobc_out(ng,inorth)
                  END IF
                  tau=tau*dt2d
                END IF
#   ifdef RADIATION_2D
                Cx=BOUNDARY(ng)%zeta_north_Cx(i)
#   else
                Cx=0.0_r8
#   endif
                Ce=BOUNDARY(ng)%zeta_north_Ce(i)
                cff=BOUNDARY(ng)%zeta_north_C2(i)

!^              zeta_new(i,Jend+1)=(cff*zeta(i,Jend+1,know)+            &
!^   &                              Ce *zeta_new(i,Jend)-               &
!^   &                              MAX(Cx,0.0_r8)*grad(i  ,Jend+1)-    &
!^   &                              MIN(Cx,0.0_r8)*grad(i+1,Jend+1))/   &
!^   &                             (cff+Ce)
!^
                tl_zeta_new(i,Jend+1)=(cff*tl_zeta(i,Jend+1,know)+      &
     &                                 Ce *tl_zeta_new(i,Jend)-         &
     &                                 MAX(Cx,0.0_r8)*                  &
     &                                    tl_grad(i  ,Jend+1)-          &
     &                                 MIN(Cx,0.0_r8)*                  &
     &                                    tl_grad(i+1,Jend+1))/         &
     &                                (cff+Ce)

                IF (tl_LBC(inorth,isFsur,ng)%nudging) THEN
!^                zeta_new(i,Jend+1)=zeta_new(i,Jend+1)+                &
!^   &                               tau*(BOUNDARY(ng)%zeta_north(i)-   &
!^   &                                    zeta(i,Jend+1,know))
!^
                  tl_zeta_new(i,Jend+1)=tl_zeta_new(i,Jend+1)-          &
     &                                  tau*tl_zeta(i,Jend+1,know)
                END IF
#  else
!^              zeta_new(i,Jend+1)=zeta_new(i,Jend)           ! gradient
!^
                tl_zeta_new(i,Jend+1)=tl_zeta_new(i,Jend)
#  endif
#  ifdef MASKING
!^              zeta_new(i,Jend+1)=zeta_new(i,Jend+1)*                  &
!^   &                             GRID(ng)%rmask(i,Jend+1)
!^
                tl_zeta_new(i,Jend+1)=tl_zeta_new(i,Jend+1)*            &
     &                                GRID(ng)%rmask(i,Jend+1)
#  endif
              END IF
            END DO
          END IF
!
!  Northern edge, explicit Chapman boundary condition.
!
        ELSE IF (tl_LBC(inorth,isFsur,ng)%Chapman_explicit) THEN
          DO i=IstrU-1,Iend
            IF (LBC_apply(ng)%north(i)) THEN
              cff=dt2d*GRID(ng)%pn(i,Jend)
              cff1=SQRT(g*(GRID(ng)%h(i,Jend)+                          &
     &                     zeta(i,Jend,know)))
              tl_cff1=0.5_r8*g*(GRID(ng)%tl_h(i,Jend)+                  &
     &                          tl_zeta(i,Jend,know))/cff1
              Ce=cff*cff1
              tl_Ce=cff*tl_cff1
!^            zeta_new(i,Jend+1)=(1.0_r8-Ce)*zeta(i,Jend+1,know)+       &
!^   &                           Ce*zeta(i,Jend,know)
!^
              tl_zeta_new(i,Jend+1)=(1.0_r8-Ce)*tl_zeta(i,Jend+1,know)+ &
     &                              tl_Ce*(zeta(i,Jend+1,know)+         &
     &                                     zeta(i,Jend  ,know))+        &
     &                              Ce*tl_zeta(i,Jend,know)
#  ifdef MASKING
!^            zeta_new(i,Jend+1)=zeta_new(i,Jend+1)*                    &
!^   &                           GRID(ng)%rmask(i,Jend+1)
!^
              tl_zeta_new(i,Jend+1)=tl_zeta_new(i,Jend+1)*              &
     &                              GRID(ng)%rmask(i,Jend+1)
#  endif
            END IF
          END DO
!
!  Northern edge, implicit Chapman boundary condition.
!
        ELSE IF (tl_LBC(inorth,isFsur,ng)%Chapman_implicit) THEN
          DO i=IstrU-1,Iend
            IF (LBC_apply(ng)%north(i)) THEN
              cff=dt2d*GRID(ng)%pn(i,Jend)
              cff1=SQRT(g*(GRID(ng)%h(i,Jend)+                          &
     &                     zeta(i,Jend,know)))
              tl_cff1=0.5_r8*g*(GRID(ng)%tl_h(i,Jend)+                  &
     &                          tl_zeta(i,Jend,know))/cff1
              Ce=cff*cff1
              tl_Ce=cff*tl_cff1
              cff2=1.0_r8/(1.0_r8+Ce)
              tl_cff2=-cff2*cff2*tl_Ce
!^            zeta_new(i,Jend+1)=cff2*(zeta(i,Jend+1,know)+             &
!^   &                                 Ce*zeta_new(i,Jend))
!^
              tl_zeta_new(i,Jend+1)=tl_cff2*(zeta(i,Jend+1,know)+       &
     &                                       Ce*zeta_new(i,Jend))+      &
     &                              cff2*(tl_zeta(i,Jend+1,know)+       &
     &                                    tl_Ce*zeta_new(i,Jend)+       &
     &                                    Ce*tl_zeta_new(i,Jend))
#  ifdef MASKING
!^            zeta_new(i,Jend+1)=zeta_new(i,Jend+1)*                    &
!^   &                           GRID(ng)%rmask(i,Jend+1)
!^
              tl_zeta_new(i,Jend+1)=tl_zeta_new(i,Jend+1)*              &
     &                              GRID(ng)%rmask(i,Jend+1)
#  endif
            END IF
          END DO
!
!  Northern edge, clamped boundary condition.
!
        ELSE IF (tl_LBC(inorth,isFsur,ng)%clamped) THEN
          DO i=IstrU-1,Iend
            IF (LBC_apply(ng)%north(i)) THEN
!^            zeta_new(i,Jend+1)=BOUNDARY(ng)%zeta_north(i)
!^
#  ifdef ADJUST_BOUNDARY
              IF (Lobc(inorth,isFsur,ng)) THEN
                tl_zeta_new(i,Jend+1)=BOUNDARY(ng)%tl_zeta_north(i)
              ELSE
                tl_zeta_new(i,Jend+1)=0.0_r8
              END IF
#  else
              tl_zeta_new(i,Jend+1)=0.0_r8
#  endif
#  ifdef MASKING
!^            zeta_new(i,Jend+1)=zeta_new(i,Jend+1)*                    &
!^   &                           GRID(ng)%rmask(i,Jend+1)
!^
              tl_zeta_new(i,Jend+1)=tl_zeta_new(i,Jend+1)*              &
     &                              GRID(ng)%rmask(i,Jend+1)
#  endif
            END IF
          END DO
!
!  Northern edge, gradient boundary condition.
!
        ELSE IF (tl_LBC(inorth,isFsur,ng)%gradient) THEN
          DO i=IstrU-1,Iend
            IF (LBC_apply(ng)%north(i)) THEN
!^            zeta_new(i,Jend+1)=zeta_new(i,Jend)
!^
              tl_zeta_new(i,Jend+1)=tl_zeta_new(i,Jend)
#  ifdef MASKING
!^            zeta_new(i,Jend+1)=zeta_new(i,Jend+1)*                    &
!^   &                           GRID(ng)%rmask(i,Jend+1)
!^
              tl_zeta_new(i,Jend+1)=tl_zeta_new(i,Jend+1)*              &
     &                              GRID(ng)%rmask(i,Jend+1)
#  endif
            END IF
          END DO
!
!  Northern edge, closed boundary condition.
!
        ELSE IF (tl_LBC(inorth,isFsur,ng)%closed) THEN
          DO i=IstrU-1,Iend
            IF (LBC_apply(ng)%north(i)) THEN
!^            zeta_new(i,Jend+1)=zeta_new(i,Jend)
!^
              tl_zeta_new(i,Jend+1)=tl_zeta_new(i,Jend)
#  ifdef MASKING
!^            zeta_new(i,Jend+1)=zeta_new(i,Jend+1)*                    &
!^   &                           GRID(ng)%rmask(i,Jend+1)
!^
              tl_zeta_new(i,Jend+1)=tl_zeta_new(i,Jend+1)*              &
     &                              GRID(ng)%rmask(i,Jend+1)
#  endif
            END IF
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (.not.(EWperiodic(ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Istr-1).and.                          &
     &        LBC_apply(ng)%west (Jstr-1)) THEN
!^          zeta_new(Istr-1,Jstr-1)=0.5_r8*(zeta_new(Istr  ,Jstr-1)+    &
!^   &                                      zeta_new(Istr-1,Jstr  ))
!^
            tl_zeta_new(Istr-1,Jstr-1)=0.5_r8*                          &
     &                                 (tl_zeta_new(Istr  ,Jstr-1)+     &
     &                                  tl_zeta_new(Istr-1,Jstr  ))
          END IF
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jstr-1)) THEN
!^          zeta_new(Iend+1,Jstr-1)=0.5_r8*(zeta_new(Iend  ,Jstr-1)+    &
!^   &                                      zeta_new(Iend+1,Jstr  ))
!^
            tl_zeta_new(Iend+1,Jstr-1)=0.5_r8*                          &
     &                                 (tl_zeta_new(Iend  ,Jstr-1)+     &
     &                                  tl_zeta_new(Iend+1,Jstr  ))
          END IF
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Istr-1).and.                          &
     &        LBC_apply(ng)%west (Jend+1)) THEN
!^          zeta_new(Istr-1,Jend+1)=0.5_r8*(zeta_new(Istr-1,Jend  )+    &
!^   &                                      zeta_new(Istr  ,Jend+1))
!^
            tl_zeta_new(Istr-1,Jend+1)=0.5_r8*                          &
     &                                 (tl_zeta_new(Istr-1,Jend  )+     &
     &                                  tl_zeta_new(Istr  ,Jend+1))
          END IF
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jend+1)) THEN
!^          zeta_new(Iend+1,Jend+1)=0.5_r8*(zeta_new(Iend+1,Jend  )+    &
!^   &                                      zeta_new(Iend  ,Jend+1))
!^
            tl_zeta_new(Iend+1,Jend+1)=0.5_r8*                          &
     &                                 (tl_zeta_new(Iend+1,Jend  )+     &
     &                                  tl_zeta_new(Iend  ,Jend+1))
          END IF
        END IF
      END IF

#  if defined WET_DRY
!
!-----------------------------------------------------------------------
! Ensure that water level on boundary cells is above bed elevation.
!-----------------------------------------------------------------------
!
      cff=Dcrit(ng)-eps
      IF (.not.EWperiodic(ng)) THEN
        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          DO j=JstrV-1,Jend
            IF (LBC_apply(ng)%west(j)) THEN
              IF (zeta_new(Istr-1,j).le.                                &
     &            (Dcrit(ng)-GRID(ng)%h(Istr-1,j))) THEN
!^              zeta_new(Istr-1,j)=cff-GRID(ng)%h(Istr-1,j)
!^
                tl_zeta_new(Istr-1,j)=-GRID(ng)%tl_h(Istr-1,j)
              END IF
            END IF
          END DO
        END IF
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO j=JstrV-1,Jend
            IF (LBC_apply(ng)%east(j)) THEN
              IF (zeta_new(Iend+1,j).le.                                &
     &            (Dcrit(ng)-GRID(ng)%h(Iend+1,j))) THEN
!^              zeta_new(Iend+1,j)=cff-GRID(ng)%h(Iend+1,j)
!^
                tl_zeta_new(Iend+1,j)=-GRID(ng)%tl_h(Iend+1,j)
              END IF
            END IF
          END DO
        END IF
      END IF
!
      IF (.not.NSperiodic(ng)) THEN
        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO i=IstrU-1,Iend
            IF (LBC_apply(ng)%south(i)) THEN
              IF (zeta_new(i,Jstr-1).le.                                &
     &            (Dcrit(ng)-GRID(ng)%h(i,Jstr-1))) THEN
!^              zeta_new(i,Jstr-1)=cff-GRID(ng)%h(i,Jstr-1)
!^
                tl_zeta_new(i,Jstr-1)=-GRID(ng)%tl_h(i,Jstr-1)
              END IF
            END IF
          END DO
        END IF
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          DO i=IstrU-1,Iend
            IF (LBC_apply(ng)%north(i)) THEN
              IF (zeta_new(i,Jend+1).le.                                &
     &            (Dcrit(ng)-GRID(ng)%h(i,Jend+1))) THEN
!^              zeta_new(i,Jend+1)=cff-GRID(ng)%h(i,Jend+1)
!^
                tl_zeta_new(i,Jend+1)=-GRID(ng)%tl_h(i,Jend+1)
              END IF
            END IF
          END DO
        END IF
      END IF
!
      IF (.not.(EWperiodic(ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Istr-1).and.                          &
     &        LBC_apply(ng)%west (Jstr-1)) THEN
            IF (zeta_new(Istr-1,Jstr-1).le.                             &
     &          (Dcrit(ng)-GRID(ng)%h(Istr-1,Jstr-1))) THEN
!^            zeta_new(Istr-1,Jstr-1)=cff-GRID(ng)%h(Istr-1,Jstr-1)
!^
              tl_zeta_new(Istr-1,Jstr-1)=-GRID(ng)%tl_h(Istr-1,Jstr-1)
            END IF
          END IF
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jstr-1)) THEN
            IF (zeta_new(Iend+1,Jstr-1).le.                             &
     &          (Dcrit(ng)-GRID(ng)%h(Iend+1,Jstr-1))) THEN
!^            zeta_new(Iend+1,Jstr-1)=cff-GRID(ng)%h(Iend+1,Jstr-1)
!^
              tl_zeta_new(Iend+1,Jstr-1)=-GRID(ng)%tl_h(Iend+1,Jstr-1)
            END IF
          END IF
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Istr-1).and.                          &
     &        LBC_apply(ng)%west (Jend+1)) THEN
            IF (zeta_new(Istr-1,Jend+1).le.                             &
     &          (Dcrit(ng)-GRID(ng)%h(Istr-1,Jend+1))) THEN
!^            zeta_new(Istr-1,Jend+1)=cff-GRID(ng)%h(Istr-1,Jend+1)
!^
              tl_zeta_new(Istr-1,Jend+1)=-GRID(ng)%tl_h(Istr-1,Jend+1)
            END IF
          END IF
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jend+1)) THEN
            IF (zeta_new(Iend+1,Jend+1).le.                             &
     &          (Dcrit(ng)-GRID(ng)%h(Iend+1,Jend+1))) THEN
!^            zeta_new(Iend+1,Jend+1)=cff-GRID(ng)%h(Iend+1,Jend+1)
!^
              tl_zeta_new(Iend+1,Jend+1)=-GRID(ng)%tl_h(Iend+1,Jend+1)
            END IF
          END IF
        END IF
      END IF
#  endif
!
      RETURN
      END SUBROUTINE tl_zetabc_local
# endif
#endif
      END MODULE tl_zetabc_mod
