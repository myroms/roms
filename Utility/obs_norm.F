#include "cppdefs.h"
#ifdef FOUR_DVAR
      SUBROUTINE obs_norm (ng)
!
!=================================================== Andrew M. Moore ===
!  Copyright (c) 2005 ROMS/TOMS Adjoint Group                          !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This subroutine computes the observation norm that it is used to    !
!  scale the cost function and adjoint misfit terms.                   !
!                                                                      !
!  WARNING:                                                            !
!  =======                                                             !
!                                                                      !
!  The  adjoint  model was derived using an L2-norm. If a different    !
!  norm is used in the definition of the  cost  function,  then the    !
!  forcing terms for adjoint model variables need to be transformed    !
!  accordingly.                                                        !
!                                                                      !
!  In 2D configurations,  a shallow water energy norm is used.  The    !
!  energy norm adjoint free surface (he) and  the  L2  norm adjoint    !
!  free surface (h2) are related according to g*he=H*h2  where H is    !
!  the undisturbed depth of the water column.                          !
!                                                                      !
!-----------------------------------------------------------------------
!
      USE mod_param
      USE mod_grid
      USE mod_fourdvar
      USE mod_ncparam
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      logical :: k_bound, r_bound, u_bound, v_bound

      integer :: i, iobs, j
# ifdef SOLVE3D
      integer :: itrc, k, order
      real(r8) :: Hzobs, Tnorm(NT(ng)), TrcCoef(NT(ng))
# endif
      real(r8) :: Dobs, cff, fac
!
!-----------------------------------------------------------------------
!  Compute norm coefficients cost function.
!-----------------------------------------------------------------------
!
      fac=0.5_r8*rho0/TotVolume
# ifdef SOLVE3D
#  ifdef N2NORM_PROFILE
      cff=1.0
#  else
      cff=1.0/bvf_bak
#  endif
      DO itrc=1,NT(ng)
        IF (itrc.eq.itemp) THEN
          TrcCoef(itrc)=Tcoef(ng)*Tcoef(ng)*g*g*cff
          Tnorm(itrc)=TrcCoef(itrc)
        ELSE IF (itrc.eq.isalt) THEN
          TrcCoef(itrc)=Scoef(ng)*Scoef(ng)*g*g*cff
          Tnorm(itrc)=TrcCoef(itrc)
        ELSE
          TrcCoef(itrc)=1.0_r8
          Tnorm(itrc)=TrcCoef(itrc)
        END IF
      END DO
# endif
      DO iobs=1,Nobs(ng)
        ObsNorm(iobs)=0.0_r8
        i=INT(Xobs(iobs))
        j=INT(Yobs(iobs))
# ifdef SOLVE3D
        k=INT(Zobs(iobs))
        k_bound=((1.le.k).and.(k.le.N(ng)))
# endif
        r_bound=((rILB(ng).le.i).and.(i.le.rIUB(ng))).and.              &
     &          ((rJLB(ng).le.j).and.(j.le.rJUB(ng)))
        u_bound=((uILB(ng).le.i).and.(i.le.uIUB(ng))).and.              &
     &          ((uJLB(ng).le.j).and.(j.le.uJUB(ng)))
        v_bound=((vILB(ng).le.i).and.(i.le.vIUB(ng))).and.              &
     &          ((vJLB(ng).le.j).and.(j.le.vJUB(ng)))
        IF (r_bound.and.                                                &
     &      (ObsType(iobs).eq.isFsur)) THEN
# if defined ENERGY1_NORM
          Dobs=Hmat(1,iobs)*GRID(ng)%h(i  ,j  )+                        &
     &         Hmat(2,iobs)*GRID(ng)%h(i+1,j  )+                        &
     &         Hmat(3,iobs)*GRID(ng)%h(i+1,j+1)+                        &
     &         Hmat(4,iobs)*GRID(ng)%h(i  ,j+1)
          IF (Dobs.ne.0.0_r8) THEN
            ObsNorm(iobs)=fac*g*GRID(ng)%omn(i,j)/Dobs
          END IF
# elif defined ENERGY2_NORM
          ObsNorm(iobs)=fac*g*GRID(ng)%omn(i,j)
# elif defined ENERGY3_NORM
          Dobs=Hmat(1,iobs)*GRID(ng)%h(i  ,j  )+                        &
     &         Hmat(2,iobs)*GRID(ng)%h(i+1,j  )+                        &
     &         Hmat(3,iobs)*GRID(ng)%h(i+1,j+1)+                        &
     &         Hmat(4,iobs)*GRID(ng)%h(i  ,j+1)
          IF (Dobs.ne.0.0_r8) THEN
            ObsNorm(iobs)=g/Dobs
          END IF
# else
          ObsNorm(iobs)=1.0_r8
# endif
# ifdef MASKING
          ObsNorm(iobs)=ObsNorm(iobs)*GRID(ng)%rmask(i,j)
# endif
        ELSE IF (u_bound.and.                                           &
     &           (ObsType(iobs).eq.isUbar)) THEN
# if defined ENERGY1_NORM
          ObsNorm(iobs)=fac*GRID(ng)%om_u(i,j)*GRID(ng)%on_u(i,j)
# elif defined ENERGY2_NORM
          Dobs=0.5_r8*(Hmat(1,iobs)*(GRID(ng)%h(i-1,j  )+               &
     &                               GRID(ng)%h(i  ,j  ))+              &
     &                 Hmat(2,iobs)*(GRID(ng)%h(i  ,j  )+               &
     &                               GRID(ng)%h(i+1,j  ))+              &
     &                 Hmat(3,iobs)*(GRID(ng)%h(i  ,j+1)+               &
     &                               GRID(ng)%h(i+1,j+1))+              &
     &                 Hmat(4,iobs)*(GRID(ng)%h(i-1,j+1)+               &
     &                               GRID(ng)%h(i  ,j+1)))
          ObsNorm(iobs)=fac*GRID(ng)%om_u(i,j)*GRID(ng)%on_u(i,j)*Dobs
# elif defined ENERGY3_NORM
          ObsNorm(iobs)=1.0_r8
# else
          ObsNorm(iobs)=1.0_r8
# endif
# ifdef MASKING
          ObsNorm(iobs)=ObsNorm(iobs)*GRID(ng)%umask(i,j)
# endif
        ELSE IF (v_bound.and.                                           &
     &           (ObsType(iobs).eq.isVbar)) THEN
# if defined ENERGY1_NORM
          ObsNorm(iobs)=fac*GRID(ng)%om_v(i,j)*GRID(ng)%on_v(i,j)
# elif defined ENERGY2_NORM
          Dobs=0.5_r8*(Hmat(1,iobs)*(GRID(ng)%h(i  ,j-1)+               &
     &                               GRID(ng)%h(i  ,j  ))+              &
     &                 Hmat(2,iobs)*(GRID(ng)%h(i+1,j-1)+               &
     &                               GRID(ng)%h(i+1,j  ))+              &
     &                 Hmat(3,iobs)*(GRID(ng)%h(i+1,j  )+               &
     &                               GRID(ng)%h(i+1,j+1))+              &
     &                 Hmat(4,iobs)*(GRID(ng)%h(i  ,j  )+               &
     &                               GRID(ng)%h(i  ,j+1)))
          ObsNorm(iobs)=fac*GRID(ng)%om_v(i,j)*GRID(ng)%on_v(i,j)*Dobs
# elif defined ENERGY3_NORM
          ObsNorm(iobs)=1.0_r8
# else
          ObsNorm(iobs)=1.0_r8
# endif
# ifdef MASKING
          ObsNorm(iobs)=ObsNorm(iobs)*GRID(ng)%vmask(i,j)
# endif
# ifdef SOLVE3D
        ELSE IF (u_bound.and.k_bound.and.                               &
     &           (ObsType(iobs).eq.isUvel)) THEN
#  if defined ENERGY1_NORM || defined ENERGY2_NORM
          Hzobs=0.5_r8*(Hmat(1,iobs)*(GRID(ng)%Hz(i-1,j  ,k)+           &
     &                                GRID(ng)%Hz(i  ,j  ,k))+          &
     &                  Hmat(2,iobs)*(GRID(ng)%Hz(i  ,j  ,k)+           &
     &                                GRID(ng)%Hz(i+1,j  ,k))+          &
     &                  Hmat(3,iobs)*(GRID(ng)%Hz(i  ,j+1,k)+           &
     &                                GRID(ng)%Hz(i+1,j+1,k))+          &
     &                  Hmat(4,iobs)*(GRID(ng)%Hz(i-1,j+1,k)+           &
     &                                GRID(ng)%Hz(i  ,j+1,k)))
#  endif
#  if defined ENERGY1_NORM
          Dobs=0.5_r8*(Hmat(1,iobs)*(GRID(ng)%h(i-1,j  )+               &
     &                               GRID(ng)%h(i  ,j  ))+              &
     &                 Hmat(2,iobs)*(GRID(ng)%h(i  ,j  )+               &
     &                               GRID(ng)%h(i+1,j  ))+              &
     &                 Hmat(3,iobs)*(GRID(ng)%h(i  ,j+1)+               &
     &                               GRID(ng)%h(i+1,j+1))+              &
     &                 Hmat(4,iobs)*(GRID(ng)%h(i-1,j+1)+               &
     &                               GRID(ng)%h(i  ,j+1)))
          IF (Dobs.ne.0.0_r8) THEN 
            ObsNorm(iobs)=fac*GRID(ng)%om_u(i,j)*GRID(ng)%on_u(i,j)*    &
     &                    Hzobs/Dobs
          END IF
#  elif defined ENERGY2_NORM
          ObsNorm(iobs)=fac*GRID(ng)%om_u(i,j)*GRID(ng)%on_u(i,j)*Hzobs
#  elif defined ENERGY3_NORM
          ObsNorm(iobs)=1.0_r8
#  else
          ObsNorm(iobs)=1.0_r8
#  endif
#  ifdef MASKING
          ObsNorm(iobs)=ObsNorm(iobs)*GRID(ng)%umask(i,j)
#  endif
        ELSE IF (v_bound.and.k_bound.and.                               &
     &           (ObsType(iobs).eq.isVvel)) THEN
#  if defined ENERGY1_NORM || defined ENERGY2_NORM
          Hzobs=0.5_r8*(Hmat(1,iobs)*(GRID(ng)%Hz(i  ,j-1,k)+           &
     &                                GRID(ng)%Hz(i  ,j  ,k))+          &
     &                  Hmat(2,iobs)*(GRID(ng)%Hz(i+1,j-1,k)+           &
     &                                GRID(ng)%Hz(i+1,j  ,k))+          &
     &                  Hmat(3,iobs)*(GRID(ng)%Hz(i+1,j  ,k)+           &
     &                                GRID(ng)%Hz(i+1,j+1,k))+          &
     &                  Hmat(4,iobs)*(GRID(ng)%Hz(i  ,j  ,k)+           &
     &                                GRID(ng)%Hz(i  ,j+1,k)))
#  endif
#  if defined ENERGY1_NORM
          Dobs=0.5_r8*(Hmat(1,iobs)*(GRID(ng)%h(i  ,j-1)+               &
     &                               GRID(ng)%h(i  ,j  ))+              &
     &                 Hmat(2,iobs)*(GRID(ng)%h(i+1,j-1)+               &
     &                               GRID(ng)%h(i+1,j  ))+              &
     &                 Hmat(3,iobs)*(GRID(ng)%h(i+1,j  )+               &
     &                               GRID(ng)%h(i+1,j+1))+              &
     &                 Hmat(4,iobs)*(GRID(ng)%h(i  ,j  )+               &
     &                               GRID(ng)%h(i  ,j+1)))
          IF (Dobs.ne.0.0_r8) THEN
            ObsNorm(iobs)=fac*GRID(ng)%om_v(i,j)*GRID(ng)%on_v(i,j)*    &
     &                    Hzobs/Dobs
          END IF
#  elif defined ENERGY2_NORM
          ObsNorm(iobs)=fac*GRID(ng)%om_v(i,j)*GRID(ng)%on_v(i,j)*Hzobs
#  elif defined ENERGY3_NORM
          ObsNorm(iobs)=1.0_r8
#  else
          ObsNorm(iobs)=1.0_r8
#  endif
#  ifdef MASKING
          ObsNorm(iobs)=ObsNorm(iobs)*GRID(ng)%vmask(i,j)
#  endif
        ELSE
          DO itrc=1,NT(ng)
            IF (r_bound.and.k_bound.and.                                &
     &          (ObsType(iobs).eq.isTvar(itrc))) THEN
#  ifdef N2NORM_PROFILE
              Dobs=Hmat(1,iobs)*GRID(ng)%z_r(i  ,j  ,k)+                &
     &             Hmat(2,iobs)*GRID(ng)%z_r(i+1,j  ,k)+                &
     &             Hmat(3,iobs)*GRID(ng)%z_r(i+1,j+1,k)+                &
     &             Hmat(4,iobs)*GRID(ng)%z_r(i  ,j+1,k)
              IF (Dobs.ne.0.0_r8) THEN
                cff=pcoef_N2(npN2)
                DO order=npN2-1,0,-1
                  cff=Dobs*cff+pcoef_N2(order)
                END DO
!!              cff=1.0E-7_r8+1.0E-5_r8*(1.0_r8+TANH(Dobs*0.002_r8))
                TrcCoef(itrc)=Tnorm(itrc)/ABS(cff)
              ELSE                               ! Ignored observation
                TrcCoef(itrc)=0.0_r8             ! probably land masked
              END IF
#  endif
#  if defined ENERGY1_NORM || defined ENERGY2_NORM
              Hzobs=Hmat(1,iobs)*GRID(ng)%Hz(i  ,j  ,k)+                &
     &              Hmat(2,iobs)*GRID(ng)%Hz(i+1,j  ,k)+                &
     &              Hmat(3,iobs)*GRID(ng)%Hz(i+1,j+1,k)+                &
     &              Hmat(4,iobs)*GRID(ng)%Hz(i  ,j+1,k)
#  endif
#  ifdef ENERGY1_NORM
              Dobs=Hmat(1,iobs)*GRID(ng)%h(i  ,j  )+                    &
     &             Hmat(2,iobs)*GRID(ng)%h(i+1,j  )+                    &
     &             Hmat(3,iobs)*GRID(ng)%h(i+1,j+1)+                    &
     &             Hmat(4,iobs)*GRID(ng)%h(i  ,j+1)
              IF (Dobs.ne.0.0_r8) THEN
                ObsNorm(iobs)=fac*TrcCoef(itrc)*GRID(ng)%omn(i,j)*      &
     &                        Hzobs/Dobs
              END IF
#  elif defined ENERGY2_NORM
              ObsNorm(iobs)=fac*TrcCoef(itrc)*GRID(ng)%omn(i,j)*Hzobs
#  elif defined ENERGY3_NORM
              ObsNorm(iobs)=TrcCoef(itrc)
#  else
              ObsNorm(iobs)=1.0_r8
#  endif
#  ifdef MASKING
              ObsNorm(iobs)=ObsNorm(iobs)*GRID(ng)%rmask(i,j)
#  endif
            END IF
          END DO
# endif
        END IF
      END DO
      RETURN
      END SUBROUTINE obs_norm
#else
      SUBROUTINE obs_norm
      RETURN
      END SUBROUTINE obs_norm
#endif
 
