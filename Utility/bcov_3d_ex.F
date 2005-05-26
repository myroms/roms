#include "cppdefs.h"

      MODULE bcov_3d_ex_mod

#if defined S4DVAR || defined IS4DVAR || defined W4DVAR
!
!=================================================== Andrew M. Moore ===
!  Copyright (c) 2005 ROMS/TOMS Adjoint Group                          !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routines compute the background error covariance for data      !
!  assimilation via a space convolution of the diffusion equation      !
!  for 3D state variables. The diffusion equation is solved using      !
!  an explicit (inefficient) algorithm.                                !
!                                                                      !
!  This routines are used to  compute the isotropic normalization      !
!  factor  Anorm via a randomization algorithm and to compute the      !
!  the background error covariance.  The  Anorm factor is used to      !
!  ensure that the  diagonal elements of the covariance are equal      !
!  to unity.                                                           !
!                                                                      !
!  For Gaussian (bell-shaped) covariances,  the space convolution      !
!  of the diffusion operator is an efficient way  to estimate the      !
!  finite domain error covariances.  The  number  of  convolution      !
!  steps determined from the values of Ldecay,  DTdif,  and Kdif.      !
!  Here,  Kdif  is a diffusion coefficient,  assumed constant and      !
!  with a value of 1 m2/s. The number of convolution steps Nsteps      !
!  is forced to be an even number because  its associated adjoint      !
!  operator needs to be convoluted for half of the steps only.         !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     Istr       Starting tile index in the I-direction.               !
!     Iend       Ending   tile index in the I-direction.               !
!     Jstr       Starting tile index in the J-direction.               !
!     Jend       Ending   tile index in the J-direction.               !
!     LBi        I-dimension lower bound.                              !
!     UBi        I-dimension upper bound.                              !
!     LBj        J-dimension lower bound.                              !
!     UBj        J-dimension upper bound.                              !
!     LBk        K-dimension lower bound.                              !
!     UBk        K-dimension upper bound.                              !
!     Nghost     Number of ghost points.                               !
!     DTdif      Time-step for diffusion operation convolution (s).    !
!     Ldecay     Spatial Gaussian decorrelation scale (m).             !
!     Anorm      Isotropic normalization factor.                       !
!     A          3D state variable to diffuse.                         !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A          Diffused 3D state variable.                           !      
!                                                                      !
!  Routines:                                                           !
!                                                                      !
!    bcov_r3d_tile     Background error covariance at RHO-points       !
!    bcov_u3d_tile     Background error covariance at U-points         !
!    bcov_v3d_tile     Background error covariance at V-points         !
!                                                                      !
!=======================================================================
!
      implicit none
 
      CONTAINS
!
!***********************************************************************
      SUBROUTINE bcov_r3d_tile (ng, model, Istr, Iend, Jstr, Jend,      &
     &                          LBi, UBi, LBj, UBj, LBk, UBk, Nghost,   &
     &                          DTdif, Ldecay, Anorm,                   &
     &                          pm, pn, pmon_u, pnom_v,                 &
# ifdef MASKING
     &                          umask, vmask,                           &
# endif
     &                          z_r, Hz,                                &
     &                          A)
!***********************************************************************
!
      USE mod_param
!
      USE bc_3d_mod, ONLY: bc_r3d_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
      integer, intent(in) :: Nghost

      real(r8), intent(in) :: DTdif, Ldecay, Anorm

# ifdef ASSUMED_SHAPE
      real(r8), intent(in)  :: pm(LBi:,LBj:)
      real(r8), intent(in)  :: pn(LBi:,LBj:)
      real(r8), intent(in)  :: pmon_u(LBi:,LBj:)
      real(r8), intent(in)  :: pnom_v(LBi:,LBj:)
#  ifdef MASKING
      real(r8), intent(in)  :: umask(LBi:,LBj:)
      real(r8), intent(in)  :: vmask(LBi:,LBj:)
#  endif
      real(r8), intent(in)  :: Hz(LBi:,LBj:,:)
      real(r8), intent(in)  :: z_r(LBi:,LBj:,:)

      real(r8), intent(inout) :: A(LBi:,LBj:,LBk:)
# else
      real(r8), intent(in)  :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: pmon_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: pnom_v(LBi:UBi,LBj:UBj)
#  ifdef MASKING
      real(r8), intent(in)  :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: vmask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in)  :: Hz(LBi:,LBj:,N(ng))
      real(r8), intent(in)  :: z_r(LBi:,LBj:,N(ng))

      real(r8), intent(inout) :: Ainp(LBi:UBi,LBj:UBj,LBk:UBk)
# endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: ILB, IUB, JLB, JUB
      integer :: Nnew, Nold, Nsteps, i, j, k, step

      real(r8) :: cff

      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY,N(ng),2) :: Awrk

      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N(ng)) :: FS

      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: DC
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: FE
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: FX

# include "set_bounds.h"
!
      ILB=LBOUND(Awrk,DIM=1)
      IUB=UBOUND(Awrk,DIM=1)
      JLB=LBOUND(Awrk,DIM=2)
      JUB=UBOUND(Awrk,DIM=2)
!
!-----------------------------------------------------------------------
!  Space convolution of the diffusion equation for a 2D state variable
!  at RHO-points
!-----------------------------------------------------------------------
!
!  Set integration indices and initial conditions.
!
      Nold=1
      Nnew=2
      DO k=1,N(ng)
        DO j=Jstr-1,Jend+1
          DO i=Istr-1,Iend+1
            Awrk(i,j,k,Nold)=A(i,j,k)*Anorm
          END DO
        END DO
      END DO
!
!  Determine number of integration steps as a function of the spatial
!  decorrelation scale, Ldecay.  A diffusion coefficient of 1 m2/s is
!  assumed.
!
      Nsteps=NINT(Ldecay*Ldecay/(4.0_r8*DTdif))
      IF (MOD(Nsteps,2).ne.0) Nsteps=Nsteps+1
!
!  Integrate diffusion-type operator.
!
      DO step=1,Nsteps
!
!  Compute vertical diffusive flux.
!
        DO j=Jstr-1,Jend+1
          DO k=1,N(ng)-1
            DO i=Istr-1,Iend+1
              cff=1.0_r8/(z_r(i,j,k+1)-z_r(i,j,k))
              FS(i,k)=cff*(Awrk(i,j,k+1,Nold)-Awrk(i,j,k,Nold))
            END DO
          END DO
          DO i=Istr,Iend
            FS(i,0)=0.0_r8
            FS(i,N(ng))=0.0_r8
          END DO
!
!  Time-step vertical diffusive term.
!
          DO i=Istr-1,Iend+1
            DC(i,j)=DTdif*pm(i,j)*pn(i,j)
          END DO
          DO k=1,N(ng)
            DO i=Istr-1,Iend+1
               Awrk(i,j,k,Nnew)=Awrk(i,j,k,Nold)*Hz(i,j,k)+             &
     &                          DC(i,j)*(FS(i,k)-FS(i,k-1))
            END DO
          END DO
        END DO
!
!  Compute XI- and ETA-components of diffusive flux.
!
        DO k=1,N(ng)
          DO j=Jstr,Jend
            DO i=Istr,Iend+1
              FX(i,j)=pmon_u(i,j)*                                      &
     &                (Awrk(i,j,k,Nold)-Awrk(i-1,j,k,Nold))
# ifdef MASKING
              FX(i,j)=FX(i,j)*umask(i,j)
# endif
            END DO
          END DO
          DO j=Jstr,Jend+1
            DO i=Istr,Iend
              FE(i,j)=pnom_v(i,j)*                                      &
     &                (Awrk(i,j,k,Nold)-Awrk(i,j-1,k,Nold))
# ifdef MASKING
              FE(i,j)=FE(i,j)*vmask(i,j)
# endif
            END DO
          END DO
!
!  Time-step horizontal harmonic diffusion term.
!
          DO j=Jstr,Jend
            DO i=Istr,Iend
              Awrk(i,j,k,Nnew)=Awrk(i,j,k,Nnew)/Hz(i,j,k)+              &
     &                         DC(i,j)*(FX(i+1,j)-FX(i,j)+              &
     &                                  FE(i,j+1)-FE(i,j))
            END DO
          END DO
        END DO
!
!  Apply boundary conditions. If applicable, exchange boundary data.
!
        CALL bc_r3d_tile (ng, model, Istr, Iend, Jstr, Jend,            &
     &                    ILB, IUB, JLB, JUB, LBk, UBk,                 &
     &                    Nghost,                                       &
     &                    Awrk(:,:,:,Nnew))
!
!  Update integration indices.
!
        Nold=Nnew
        Nnew=3-Nold
      END DO
!
!  Load solution.
!
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            A(i,j,k)=Awrk(i,j,k,Nold)
          END DO
        END DO
      END DO

      RETURN
      END SUBROUTINE bcov_r3d_tile
!
!***********************************************************************
      SUBROUTINE bcov_u3d_tile (ng, model, Istr, Iend, Jstr, Jend,      &
     &                          LBi, UBi, LBj, UBj, LBk, UBk, Nghost,   &
     &                          DTdif, Ldecay, Anorm,                   &
     &                          pm, pn, pmon_r, pnom_p,                 &
# ifdef MASKING
     &                          pmask,                                  &
# endif
     &                          z_r, Hz,                                &
     &                          A)
!***********************************************************************
!
      USE mod_param
!
      USE bc_3d_mod, ONLY: bc_u3d_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
      integer, intent(in) :: Nghost

      real(r8), intent(in) :: DTdif, Ldecay, Anorm

# ifdef ASSUMED_SHAPE
      real(r8), intent(in)  :: pm(LBi:,LBj:)
      real(r8), intent(in)  :: pn(LBi:,LBj:)
      real(r8), intent(in)  :: pmon_r(LBi:,LBj:)
      real(r8), intent(in)  :: pnom_p(LBi:,LBj:)
#  ifdef MASKING
      real(r8), intent(in)  :: pmask(LBi:,LBj:)
#  endif
      real(r8), intent(in)  :: Hz(LBi:,LBj:,:)
      real(r8), intent(in)  :: z_r(LBi:,LBj:,:)

      real(r8), intent(inout) :: A(LBi:,LBj:,LBk:)
# else
      real(r8), intent(in)  :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: pmon_r(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: pnom_p(LBi:UBi,LBj:UBj)
#  ifdef MASKING
      real(r8), intent(in)  :: pmask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in)  :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in)  :: z_r(LBi:UBi,LBj:UBj,N(ng))

      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj,LBk:UBk)
# endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: ILB, IUB, JLB, JUB
      integer :: Nnew, Nold, Nsteps, i, j, k, step

      real(r8) :: cff

      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY,N(ng),2) :: Awrk

      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N(ng)) :: FS

      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: DC
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: FE
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: FX

# include "set_bounds.h"
!
      ILB=LBOUND(Awrk,DIM=1)
      IUB=UBOUND(Awrk,DIM=1)
      JLB=LBOUND(Awrk,DIM=2)
      JUB=UBOUND(Awrk,DIM=2)
!
!-----------------------------------------------------------------------
!  Space convolution of the diffusion equation for a 3D state variable
!  at U-points
!-----------------------------------------------------------------------
!
!  Set integration indices and initial conditions.
!
      Nold=1
      Nnew=2
      DO k=1,N(ng)
        DO j=Jstr-1,Jend+1
          DO i=IstrU-1,Iend+1
            Awrk(i,j,k,Nold)=A(i,j,k)*Anorm
          END DO
        END DO
      END DO
!
!  Determine number of integration steps as a function of the spatial
!  decorrelation scale, Ldecay.  A diffusion coefficient of 1 m2/s is
!  assumed.
!
      Nsteps=NINT(Ldecay*Ldecay/(4.0_r8*DTdif))
      IF (MOD(Nsteps,2).ne.0) Nsteps=Nsteps+1
!
!  Integrate diffusion-type operator.
!
      DO step=1,Nsteps
!
!  Compute vertical diffusive flux.
!
        DO j=Jstr-1,Jend+1
          DO k=1,N(ng)-1
            DO i=IstrU-1,Iend+1
              cff=1.0_r8/(0.5_r8*(z_r(i-1,j,k+1)+z_r(i,j,k+1)-          &
     &                            z_r(i-1,j,k  )-z_r(i,j,k  )))
              FS(i,k)=cff*(Awrk(i,j,k+1,Nold)-Awrk(i,j,k,Nold))
            END DO
          END DO
          DO i=IstrU-1,Iend+1
            FS(i,0)=0.0_r8
            FS(i,N(ng))=0.0_r8
          END DO
!
!  Time-step vertical diffusive term.
!
          cff=DTdif*0.25_r8
          DO i=IstrU-1,Iend+1
            DC(i,j)=cff*(pm(i-1,j)+pm(i,j))*(pn(i-1,j)+pn(i,j))
          END DO
          DO k=1,N(ng)
            DO i=IstrU-1,Iend+1
               Awrk(i,j,k,Nnew)=Awrk(i,j,k,Nold)*                       &
     &                          0.5_r8*(Hz(i-1,j,k)+Hz(i,j,k))+         &
     &                          DC(i,j)*(FS(i,k)-FS(i,k-1))
            END DO
          END DO
        END DO
!
!  Compute XI- and ETA-components of diffusive flux.
!
        DO k=1,N(ng)
          DO j=Jstr,Jend
            DO i=IstrU-1,Iend
              FX(i,j)=pmon_r(i,j)*                                      &
     &                (Awrk(i+1,j,k,Nold)-Awrk(i,j,k,Nold))
            END DO
          END DO
          DO j=Jstr,Jend+1
            DO i=IstrU,Iend
              FE(i,j)=pnom_p(i,j)*                                      &
     &                (Awrk(i,j,k,Nold)-Awrk(i,j-1,k,Nold))
# ifdef MASKING
              FE(i,j)=FE(i,j)*pmask(i,j)
# endif
            END DO
          END DO
!
!  Time-step horizontal harmonic diffusion operator.
!
          DO j=Jstr,Jend
            DO i=IstrU,Iend
              Awrk(i,j,k,Nnew)=Awrk(i,j,k,Nnew)/                        &
     &                         (0.5_r8*(Hz(i-1,j,k)+Hz(i,j,k)))+        &
     &                         DC(i,j)*(FX(i,j)-FX(i-1,j)+              &
     &                                  FE(i,j+1)-FE(i,j))
            END DO
          END DO
        END DO
!
!  Apply boundary conditions. If applicable, exchange boundary data.
!
        CALL bc_u3d_tile (ng, model, Istr, Iend, Jstr, Jend,            &
     &                    ILB, IUB, JLB, JUB, LBk, UBk,                 &
     &                    Nghost,                                       &
     &                    Awrk(:,:,:,Nnew))
!
!  Update integration indices.
!
        Nold=Nnew
        Nnew=3-Nold
      END DO
!
!  Load solution.
!
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=IstrU-1,IendR
            A(i,j,k)=Awrk(i,j,k,Nold)
          END DO
        END DO
      END DO

      RETURN
      END SUBROUTINE bcov_u3d_tile
!
!***********************************************************************
      SUBROUTINE bcov_v3d_tile (ng, model, Istr, Iend, Jstr, Jend,      &
     &                          LBi, UBi, LBj, UBj, LBk, UBk, Nghost,   &
     &                          DTdif, Ldecay, Anorm,                   &
     &                          pm, pn, pmon_p, pnom_r,                 &
# ifdef MASKING
     &                          pmask,                                  &
# endif
     &                          z_r, Hz,                                &
     &                          A)
!***********************************************************************
!
      USE mod_param
!
      USE bc_3d_mod, ONLY: bc_v3d_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
      integer, intent(in) :: Nghost

      real(r8), intent(in) :: DTdif, Ldecay, Anorm

# ifdef ASSUMED_SHAPE
      real(r8), intent(in)  :: pm(LBi:,LBj:)
      real(r8), intent(in)  :: pn(LBi:,LBj:)
      real(r8), intent(in)  :: pmon_p(LBi:,LBj:)
      real(r8), intent(in)  :: pnom_r(LBi:,LBj:)
#  ifdef MASKING
      real(r8), intent(in)  :: pmask(LBi:,LBj:)
#  endif
      real(r8), intent(in)  :: Hz(LBi:,LBj:,:)
      real(r8), intent(in)  :: z_r(LBi:,LBj:,:)

      real(r8), intent(inout) :: A(LBi:,LBj:,LBk:)
# else
      real(r8), intent(in)  :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: pmon_p(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: pnom_r(LBi:UBi,LBj:UBj)
#  ifdef MASKING
      real(r8), intent(in)  :: pmask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in)  :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in)  :: z_r(LBi:UBi,LBj:UBj,N(ng))

      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj,LBk:UBk)
# endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: ILB, IUB, JLB, JUB
      integer :: Nnew, Nold, Nsteps, i, j, k, step

      real(r8) :: cff

      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY,N(ng),2) :: Awrk

      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N(ng)) :: FS

      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: DC
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: FE
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: FX

# include "set_bounds.h"
!
      ILB=LBOUND(Awrk,DIM=1)
      IUB=UBOUND(Awrk,DIM=1)
      JLB=LBOUND(Awrk,DIM=2)
      JUB=UBOUND(Awrk,DIM=2)
!
!-----------------------------------------------------------------------
!  Space convolution of the diffusion equation for a 3D state variable
!  at V-points
!-----------------------------------------------------------------------
!
!  Set integration indices and initial conditions.
!
      Nold=1
      Nnew=2
      DO k=1,N(ng)
        DO j=JstrV-1,Jend+1
          DO i=Istr-1,Iend+1
            Awrk(i,j,k,Nold)=A(i,j,k)*Anorm
          END DO
        END DO
      END DO
!
!  Determine number of integration steps as a function of the spatial
!  decorrelation scale, Ldecay.  A diffusion coefficient of 1 m2/s is
!  assumed.
!
      Nsteps=NINT(Ldecay*Ldecay/(4.0_r8*DTdif))
      IF (MOD(Nsteps,2).ne.0) Nsteps=Nsteps+1
!
!  Integrate diffusion-type operator.
!
      DO step=1,Nsteps
!
!  Compute vertical diffusive flux.
!
        DO j=JstrV-1,Jend+1
          DO k=1,N(ng)-1
            DO i=Istr-1,Iend+1
              cff=1.0_r8/(0.5_r8*(z_r(i,j-1,k+1)+z_r(i,j,k+1)-          &
     &                            z_r(i,j-1,k  )-z_r(i,j,k  )))
              FS(i,k)=cff*(Awrk(i,j,k+1,Nold)-Awrk(i,j,k,Nold))
            END DO
          END DO
          DO i=Istr-1,Iend+1
            FS(i,0)=0.0_r8
            FS(i,N(ng))=0.0_r8
          END DO
!
!  Time-step vertical diffusive term.
!
          cff=DTdif*0.25_r8
          DO i=Istr-1,Iend+1
            DC(i,j)=cff*(pm(i,j-1)+pm(i,j))*(pn(i,j-1)+pn(i,j))
          END DO
          DO k=1,N(ng)
            DO i=Istr-1,Iend+1
               Awrk(i,j,k,Nnew)=Awrk(i,j,k,Nold)*                       &
     &                          0.5_r8*(Hz(i,j-1,k)+Hz(i,j,k))+         &
     &                          DC(i,j)*(FS(i,k)-FS(i,k-1))
            END DO
          END DO
        END DO
!
!  Compute XI- and ETA-components of diffusive flux.
!
        DO k=1,N(ng)
          DO j=JstrV,Jend
            DO i=Istr,Iend+1
              FX(i,j)=pmon_p(i,j)*                                      &
     &                (Awrk(i,j,k,Nold)-Awrk(i-1,j,k,Nold))
# ifdef MASKING
              FX(i,j)=FX(i,j)*pmask(i,j)
# endif
            END DO
          END DO
          DO j=JstrV-1,Jend
            DO i=Istr,Iend
              FE(i,j)=pnom_r(i,j)*                                      &
     &                (Awrk(i,j+1,k,Nold)-Awrk(i,j,k,Nold))
            END DO
          END DO
!
!  Time-step horizontal harmonic diffusion term.
!
          DO j=JstrV,Jend
            DO i=Istr,Iend
              Awrk(i,j,k,Nnew)=Awrk(i,j,k,Nnew)/                        &
     &                         (0.5_r8*(Hz(i,j-1,k)+Hz(i,j,k)))+        &
     &                         DC(i,j)*(FX(i+1,j)-FX(i,j)+              &
     &                                  FE(i,j)-FE(i,j-1))
            END DO
          END DO
        END DO
!
!  Apply boundary conditions. If applicable, exchange boundary data.
!
        CALL bc_v3d_tile (ng, model, Istr, Iend, Jstr, Jend,            &
     &                    ILB, IUB, JLB, JUB, LBk, UBk,                 &
     &                    Nghost,                                       &
     &                    Awrk(:,:,:,Nnew))
!
!  Update integration indices.
!
        Nold=Nnew
        Nnew=3-Nold
      END DO
!
!  Load solution.
!
      DO k=1,N(ng)
        DO j=JstrV-1,JendR
          DO i=IstrR,IendR
            A(i,j,k)=Awrk(i,j,k,Nold)
          END DO
        END DO
      END DO

      RETURN
      END SUBROUTINE bcov_v3d_tile
#endif
      END MODULE bcov_3d_ex_mod
