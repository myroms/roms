      SUBROUTINE prsgrd (ng, tile)
!
!****************************************** Alexander F. Shchepetkin ***
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!***********************************************************************
!                                                                      !
!  This subroutine evaluates the baroclinic, hydrostatic pressure      !
!  gradient term using a nonconservative Density-Jacobian scheme,      !
!  based on cubic polynomial fits for RHO and Z_R as functions of      !
!  nondimensional coordinates (XI,ETA,s), that is, its respective      !
!  array indices.  The cubic polynomials are monotonized by using      !
!  harmonic mean instead of linear averages to interpolate slopes.     !
!  This scheme retains exact anti-symmetry J(rho,z_r)=-J(z_r,rho).     !
!                                                                      !
!  If parameter OneFifth (see below) is set to zero,  the scheme       !
!  becomes identical to standard Jacobian.                             !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    Shchepetkin A.F and J.C. McWilliams, 2003:  A method for          !
!      computing horizontal pressure gradient force in an ocean        !
!      model with non-aligned vertical coordinate, JGR, 108,           !
!      1-34.                                                           !
!                                                                      !
!***********************************************************************
!
      USE mod_param
# ifdef DIAGNOSTICS
      USE mod_diags
# endif
      USE mod_grid
      USE mod_ocean
      USE mod_stepping
!
      integer, intent(in) :: ng, tile

# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, 23)
# endif
      CALL prsgrd_tile (ng, Istr, Iend, Jstr, Jend,                     &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  nrhs(ng),                                       &
# ifdef MASKING
     &                  GRID(ng) % umask,                               &
     &                  GRID(ng) % vmask,                               &
# endif
     &                  GRID(ng) % Hz,                                  &
     &                  GRID(ng) % om_v,                                &
     &                  GRID(ng) % on_u,                                &
     &                  GRID(ng) % z_r,                                 &
     &                  GRID(ng) % z_w,                                 &
# ifdef ICESHELF
     &                  GRID(ng) % IcePress,                            &
     &                  GRID(ng) % zice,                                &
# endif
     &                  OCEAN(ng) % rho,                                &
# ifdef DIAGNOSTICS_UV
     &                  DIAGS(ng) % DiaRU,                              &
     &                  DIAGS(ng) % DiaRV,                              &
# endif
     &                  OCEAN(ng) % ru,                                 &
     &                  OCEAN(ng) % rv)
# ifdef PROFILE
      CALL wclock_off (ng, 23)
# endif
      RETURN
      END SUBROUTINE prsgrd
!
!***********************************************************************
      SUBROUTINE prsgrd_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        nrhs,                                     &
# ifdef MASKING
     &                        umask, vmask,                             &
# endif
     &                        Hz, om_v, on_u, z_r, z_w,                 &
     &                        rho,                                      &
# ifdef ICESHELF
     &                        IcePress, zice,                           &
# endif
# ifdef DIAGNOSTICS_UV
     &                        DiaRU, DiaRV,                             &
# endif
     &                        ru, rv)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: nrhs

# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
#  ifdef ICESHELF
      real(r8), intent(in) :: IcePress(LBi:,LBj:)
      real(r8), intent(in) :: zice(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: rho(LBi:,LBj:,:)
#  ifdef DIAGNOSTICS_UV
      real(r8), intent(inout) :: DiaRU(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: DiaRV(LBi:,LBj:,:,:,:)
#  endif
      real(r8), intent(inout) :: ru(LBi:,LBj:,0:,:)
      real(r8), intent(inout) :: rv(LBi:,LBj:,0:,:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: om_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
#  ifdef ICESHELF
      real(r8), intent(in) :: IcePress(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: zice(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: rho(LBi:UBi,LBj:UBj,N(ng))
#  ifdef DIAGNOSTICS_UV
      real(r8), intent(inout) :: DiaRU(LBi:UBi,LBj:UBj,N(ng),2,NDrhs)
      real(r8), intent(inout) :: DiaRV(LBi:UBi,LBj:UBj,N(ng),2,NDrhs)
#  endif
      real(r8), intent(inout) :: ru(LBi:UBi,LBj:UBj,0:N(ng),2)
      real(r8), intent(inout) :: rv(LBi:UBi,LBj:UBj,0:N(ng),2)
# endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j, k

      real(r8), parameter :: OneFifth = 0.2_r8
      real(r8), parameter :: OneTwelfth = 1.0_r8/12.0_r8
      real(r8), parameter :: eps = 1.0E-10_r8

      real(r8) :: GRho, GRho0,  HalfGRho, cff, cff1

      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY,N(ng)) :: P

      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N(ng)) :: dR
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N(ng)) :: dZ

      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: FC
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: aux
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: dRx
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: dZx
!
# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Preliminary step (same for XI- and ETA-components:
!-----------------------------------------------------------------------
!
      GRho=g/rho0
      GRho0=1000.0_r8*GRho
      HalfGRho=0.5_r8*GRho
!
      DO j=JstrV-1,Jend
        DO k=1,N(ng)-1
          DO i=IstrU-1,Iend
            dR(i,k)=rho(i,j,k+1)-rho(i,j,k)
            dZ(i,k)=z_r(i,j,k+1)-z_r(i,j,k)
          END DO
        END DO
        DO i=IstrU-1,Iend
          dR(i,N(ng))=dR(i,N(ng)-1)
          dZ(i,N(ng))=dZ(i,N(ng)-1)
          dR(i,0)=dR(i,1)
          dZ(i,0)=dZ(i,1)
        END DO
        DO k=N(ng),1,-1
          DO i=IstrU-1,Iend
            cff=2.0_r8*dR(i,k)*dR(i,k-1)
            IF (cff.gt.eps) THEN
              dR(i,k)=cff/(dR(i,k)+dR(i,k-1))
            ELSE
              dR(i,k)=0.0_r8
            END IF
            dZ(i,k)=2.0_r8*dZ(i,k)*dZ(i,k-1)/(dZ(i,k)+dZ(i,k-1))
          END DO
        END DO
        DO i=IstrU-1,Iend
# ifdef ICESHELF
          cff=z_w(i,j,N(ng))-zice(i,j)
          P(i,j,N(ng))=GRho0*cff+                                       &
     &                 GRho*(rho(i,j,N(ng))*cff+IcePress(i,j))+         &
     &                 GRho*(rho(i,j,N(ng))+                            &
     &                       0.5_r8*(rho(i,j,N(ng))-rho(i,j,N(ng)-1))*  &
     &                       (z_w(i,j,N(ng))-z_r(i,j,N(ng)))/           &
     &                       (z_r(i,j,N(ng))-z_r(i,j,N(ng)-1)))*        &
     &                 (z_w(i,j,N(ng))-z_r(i,j,N(ng)))
# else
          P(i,j,N(ng))=GRho0*z_w(i,j,N(ng))+                            &
     &                 GRho*(rho(i,j,N(ng))+                            &
     &                       0.5_r8*(rho(i,j,N(ng))-rho(i,j,N(ng)-1))*  &
     &                       (z_w(i,j,N(ng))-z_r(i,j,N(ng)))/           &
     &                       (z_r(i,j,N(ng))-z_r(i,j,N(ng)-1)))*        &
     &                 (z_w(i,j,N(ng))-z_r(i,j,N(ng)))
# endif
        END DO
        DO k=N(ng)-1,1,-1
          DO i=IstrU-1,Iend
            P(i,j,k)=P(i,j,k+1)+                                        &
     &               HalfGRho*((rho(i,j,k+1)+rho(i,j,k))*               &
     &                          (z_r(i,j,k+1)-z_r(i,j,k))               &
     &              -OneFifth*((dR(i,k+1)-dR(i,k))*                     &
     &                         (z_r(i,j,k+1)-z_r(i,j,k)-                &
     &                         OneTwelfth*(dZ(i,k+1)+dZ(i,k)))-         &
     &                         (dZ(i,k+1)-dZ(i,k))*                     &
     &                         (rho(i,j,k+1)-rho(i,j,k)-                &
     &                         OneTwelfth*(dR(i,k+1)+dR(i,k)))))
          END DO
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Compute XI-component pressure gradient term.
!-----------------------------------------------------------------------
!
      DO k=N(ng),1,-1
        DO j=Jstr,Jend
          DO i=IstrU-1,Iend+1
            aux(i,j)=(z_r(i,j,k)-z_r(i-1,j,k))
# ifdef MASKING
            aux(i,j)=aux(i,j)*umask(i,j)
# endif
            FC(i,j)=(rho(i,j,k)-rho(i-1,j,k))
# ifdef MASKING
            FC(i,j)=FC(i,j)*umask(i,j)
# endif
          END DO
        END DO
!
        DO j=Jstr,Jend
          DO i=IstrU-1,Iend
            cff=2.0_r8*aux(i,j)*aux(i+1,j)
            IF (cff.gt.eps) THEN
              dZx(i,j)=cff/(aux(i,j)+aux(i+1,j))
            ELSE
              dZx(i,j)=0.0_r8
            END IF
            cff1=2.0_r8*FC(i,j)*FC(i+1,j)
            IF (cff1.gt.eps) THEN
              dRx(i,j)=cff1/(FC(i,j)+FC(i+1,j))
            ELSE
              dRx(i,j)=0.0_r8
            END IF
          END DO
        END DO
!
        DO j=Jstr,Jend
          DO i=IstrU,Iend
            ru(i,j,k,nrhs)=0.5_r8*(Hz(i,j,k)+Hz(i-1,j,k))*on_u(i,j)*    &
     &                     (P(i-1,j,k)-P(i,j,k)-                        &
     &                      HalfGRho*((rho(i,j,k)+rho(i-1,j,k))*        &
     &                                (z_r(i,j,k)-z_r(i-1,j,k))-        &
     &                                OneFifth*                         &
     &                                ((dRx(i,j)-dRx(i-1,j))*           &
     &                                 (z_r(i,j,k)-z_r(i-1,j,k)-        &
     &                                  OneTwelfth*                     &
     &                                  (dZx(i,j)+dZx(i-1,j)))-         &
     &                                 (dZx(i,j)-dZx(i-1,j))*           &
     &                                 (rho(i,j,k)-rho(i-1,j,k)-        &
     &                                  OneTwelfth*                     &
     &                                  (dRx(i,j)+dRx(i-1,j))))))
# ifdef DIAGNOSTICS_UV
            DiaRU(i,j,k,nrhs,M3pgrd)=ru(i,j,k,nrhs)
# endif
          END DO
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  ETA-component pressure gradient term.
!-----------------------------------------------------------------------
!
      DO k=N(ng),1,-1
        DO j=JstrV-1,Jend+1
          DO i=Istr,Iend
            aux(i,j)=(z_r(i,j,k)-z_r(i,j-1,k))
# ifdef MASKING
            aux(i,j)=aux(i,j)*vmask(i,j)
# endif
            FC(i,j)=(rho(i,j,k)-rho(i,j-1,k))
# ifdef MASKING
            FC(i,j)=FC(i,j)*vmask(i,j)
# endif
          END DO
        END DO
!
        DO j=JstrV-1,Jend
          DO i=Istr,Iend
            cff=2.0_r8*aux(i,j)*aux(i,j+1)
            IF (cff.gt.eps) THEN
              dZx(i,j)=cff/(aux(i,j)+aux(i,j+1))
            ELSE
              dZx(i,j)=0.0_r8
            END IF
            cff1=2.0_r8*FC(i,j)*FC(i,j+1)
            IF (cff1.gt.eps) THEN
              dRx(i,j)=cff1/(FC(i,j)+FC(i,j+1))
            ELSE
              dRx(i,j)=0.0_r8
            END IF
          END DO
        END DO
!
        DO j=JstrV,Jend
          DO i=Istr,Iend
            rv(i,j,k,nrhs)=0.5_r8*(Hz(i,j,k)+Hz(i,j-1,k))*om_v(i,j)*    &
     &                     (P(i,j-1,k)-P(i,j,k)-                        &
     &                      HalfGRho*((rho(i,j,k)+rho(i,j-1,k))*        &
     &                                (z_r(i,j,k)-z_r(i,j-1,k))-        &
     &                                OneFifth*                         &
     &                                ((dRx(i,j)-dRx(i,j-1))*           &
     &                                 (z_r(i,j,k)-z_r(i,j-1,k)-        &
     &                                  OneTwelfth*                     &
     &                                  (dZx(i,j)+dZx(i,j-1)))-         &
     &                                 (dZx(i,j)-dZx(i,j-1))*           &
     &                                 (rho(i,j,k)-rho(i,j-1,k)-        &
     &                                  OneTwelfth*                     &
     &                                  (dRx(i,j)+dRx(i,j-1))))))
# ifdef DIAGNOSTICS_UV
            DiaRV(i,j,k,nrhs,M3pgrd)=rv(i,j,k,nrhs)
# endif
          END DO
        END DO
      END DO
      RETURN
      END SUBROUTINE prsgrd_tile
