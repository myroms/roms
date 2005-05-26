#include "cppdefs.h"
      MODULE extract_obs_mod
#ifdef FOUR_DVAR
!
!=======================================================================
!  Copyright (c) 2005 ROMS/TOMS Adjoint Group                          !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine extracts model data at the requested observations      !
!  positions (Xobs,Yobs,Zobs).  The extraction is done via linear      !
!  interpolation. The (Xobs,Yobs) positions must be in fractional      !
!  grid coordinates.  Zobs can be in fractional  grid coordinates      !
!  (Zobs >= 0) or actual depths (Zobs < 0), if applicable.             !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     LBi        I-dimension Lower bound.                              !
!     UBi        I-dimension Upper bound.                              !
!     LBj        J-dimension Lower bound.                              !
!     UBj        J-dimension Upper bound.                              !
!     ILB        Global lower bound I-coordinate to consider.          !
!     IUB        Global upper bound I-coordinate to consider.          !
!     JLB        Global lower bound J-coordinate to consider.          !
!     JUB        Global upper bound J-coordinate to consider.          !
!     ifield     State field identification to process.                !
!     Mobs       Observation dimension in the calling program.         !
!     Nobs       Number of observations to process.                    !
!     Nweights   Number of interpolation weights per datum.            !
!     ObsType    Observations type.                                    !
!     Xobs       Observations X-locations (grid coordinates).          !
!     Yobs       Observations Y-locations (grid coordinates).          !
!     Zobs       Observations Z-locations (grid coordinates or meters).!
!     Istr       I-coordinate lower bound threshold for requested      !
!                  state field type.                                   !
!     Iend       I-coordinate upper bound threshold for requested      !
!                  state field type.                                   !
!     Jstr       J-coordinate lower bound threshold for requested      !
!                  state field type.                                   !
!     Jend       J-coordinate upper bound threshold for requested      !
!                  state field type.                                   !
!     LBk        Model K-dimension Lower bound.                        !
!     UBk        Model K-dimension Upper bound.                        !
!     A          Model array (2D or 3D) to process.                    !
!     Adepth     Depths (meter; negative).                             !
!     Amask      Land-sea masking.                                     !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Aobs       Extracted model values at observation positions.      ! 
!     Zobs       Observations Z-locations (grid coordinates).          !
!     Hmat       Interpolation weights matrix, Hmat(1:8,:),            !
!                                                                      !
!                               8____________7                         !
!                               /.          /| (i2,j2,k2)              !
!                              / .         / |                         !
!                            5/___________/6 |                         !
!                             |  .        |  |                         !
!                             |  .        |  |         Grid Cell       !
!                             | 4.........|..|3                        !
!                             | .         |  /                         !
!                             |.          | /                          !
!                  (i1,j1,k1) |___________|/                           !
!                             1           2                            !
!                                                                      !
!=======================================================================
!
      USE mod_kinds

      implicit none

      PUBLIC extract_obs2d
# ifdef SOLVE3D
      PUBLIC extract_obs3d
# endif

      CONTAINS
!
!***********************************************************************
      SUBROUTINE extract_obs2d (ng, Istr, Iend, Jstr, Jend,             &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          ILB, IUB, JLB, JUB,                     &
     &                          ifield, Mobs, Nobs, Nweights,           &
     &                          ObsType, Xobs, Yobs,                    &
     &                          A,                                      &
# ifdef MASKING
     &                          Amask,                                  &
# endif
     &                          Aobs, Hmat)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Istr, Iend, Jstr, Jend
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: ILB, IUB, JLB, JUB
      integer, intent(in) :: ifield, Mobs, Nobs, Nweights

# ifdef ASSUMED_SHAPE
      integer, intent(in) :: ObsType(:)

      real(r8), intent(in) :: Xobs(:)
      real(r8), intent(in) :: Yobs(:)
      real(r8), intent(in) :: A(LBi:,LBj:)
#  ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:,LBj:)
#  endif
      real(r8), intent(inout) :: Aobs(:)
      real(r8), intent(inout) :: Hmat(:,:)
# else
      integer, intent(in) :: ObsType(Mobs)

      real(r8), intent(in) :: Xobs(Mobs)
      real(r8), intent(in) :: Yobs(Mobs)
      real(r8), intent(in) :: A(LBi:UBi,LBj:UBj)
#  ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(inout) :: Aobs(Mobs)
      real(r8), intent(inout) :: Hmat(Nweights,Mobs)
# endif
!
!  Local variable declarations.
!
      integer :: ic, iobs, i1, i2, j1, j2

      real(r8) :: p1, p2, q1, q2, wsum
!
!-----------------------------------------------------------------------
!  Interpolate from requested 2D state field when appropriate.
!-----------------------------------------------------------------------
!
      DO iobs=1,Nobs
        i1=INT(Xobs(iobs))
        i2=i1+1
        j1=INT(Yobs(iobs))
        j2=j1+1
        IF ((ObsType(iobs).eq.ifield).and.                              &
     &      ((Istr.le.i1).and.(i1.le.Iend)).and.                        &
     &      ((ILB .le.i1).and.(i1.le.IUB )).and.                        &
     &      ((Jstr.le.j1).and.(j1.le.Jend)).and.                        &
     &      ((JLB .le.j1).and.(j1.le.JUB ))) THEN
          p2=REAL(i2-i1,r8)*(Xobs(iobs)-REAL(i1,r8))
          q2=REAL(j2-j1,r8)*(Yobs(iobs)-REAL(j1,r8))
          p1=1.0_r8-p2
          q1=1.0_r8-q2
          Hmat(1,iobs)=p1*q1
          Hmat(2,iobs)=p2*q1
          Hmat(3,iobs)=p2*q2
          Hmat(4,iobs)=p1*q2
# ifdef MASKING
          Hmat(1,iobs)=Hmat(1,iobs)*Amask(i1,j1)
          Hmat(2,iobs)=Hmat(2,iobs)*Amask(i2,j1)
          Hmat(3,iobs)=Hmat(3,iobs)*Amask(i2,j2)
          Hmat(4,iobs)=Hmat(4,iobs)*Amask(i1,j2)
          wsum=0.0_r8
          DO ic=1,4
            wsum=wsum+Hmat(ic,iobs)
          END DO
          IF (wsum.gt.0.0_r8) THEN
            wsum=1.0_r8/wsum
            DO ic=1,4 
              Hmat(ic,iobs)=Hmat(ic,iobs)*wsum 
            END DO
          END IF          
# endif
          Aobs(iobs)=Hmat(1,iobs)*A(i1,j1)+                             &
     &               Hmat(2,iobs)*A(i2,j1)+                             &
     &               Hmat(3,iobs)*A(i2,j2)+                             &
     &               Hmat(4,iobs)*A(i1,j2)
        END IF
      END DO

      RETURN
      END SUBROUTINE extract_obs2d

# ifdef SOLVE3D
!
!***********************************************************************
      SUBROUTINE extract_obs3d (ng, Istr, Iend, Jstr, Jend,             &
     &                          LBi, UBi, LBj, UBj, LBk, UBk,           &
     &                          ILB, IUB, JLB, JUB,                     &
     &                          ifield, Mobs, Nobs, Nweights,           &
     &                          ObsType, Xobs, Yobs, Zobs,              &
     &                          A, Adepth,                              &
#  ifdef MASKING
     &                          Amask,                                  &
#  endif
     &                          Aobs, Hmat)
!***********************************************************************
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Istr, Iend, Jstr, Jend
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
      integer, intent(in) :: ILB, IUB, JLB, JUB
      integer, intent(in) :: ifield, Mobs, Nobs, Nweights

#  ifdef ASSUMED_SHAPE
      integer, intent(in) :: ObsType(:)

      real(r8), intent(in) :: Xobs(:)
      real(r8), intent(in) :: Yobs(:)
      real(r8), intent(in) :: A(LBi:,LBj:,LBk:)
      real(r8), intent(in) :: Adepth(LBi:,LBj:,LBk:)
#   ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:,LBj:)
#   endif
      real(r8), intent(inout) :: Zobs(:)
      real(r8), intent(inout) :: Aobs(:)
      real(r8), intent(inout) :: Hmat(:,:)
#  else
      integer, intent(in) :: ObsType(Mobs)

      real(r8), intent(in) :: Xobs(Mobs)
      real(r8), intent(in) :: Yobs(Mobs)
      real(r8), intent(in) :: A(LBi:UBi,LBj:UBj,LBk:UBk)
      real(r8), intent(in) :: Adepth(LBi:UBi,LBj:UBj,LBk:UBk)
#   ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(inout) :: Zobs(Mobs)
      real(r8), intent(inout) :: Aobs(Mobs)
      real(r8), intent(inout) :: Hmat(Nweights,Mobs)
#  endif
!
!  Local variable declarations.
!
      integer :: i, ic, iobs, i1, i2, j1, j2, k, k1, k2

      real(r8) :: Zbot, Ztop, dz, p1, p2, q1, q2, r1, r2
      real(r8) :: w11, w12, w21, w22, wsum
!
!-----------------------------------------------------------------------
!  Interpolate from requested 3D state field.
!-----------------------------------------------------------------------
!
      DO iobs=1,Nobs
        i1=INT(Xobs(iobs))
        i2=i1+1
        j1=INT(Yobs(iobs))
        j2=j1+1
        IF ((ObsType(iobs).eq.ifield).and.                              &
     &      ((Istr.le.i1).and.(i1.le.Iend)).and.                        &
     &      ((ILB .le.i1).and.(i1.le.IUB )).and.                        &
     &      ((Jstr.le.j1).and.(j1.le.Jend)).and.                        &
     &      ((JLB .le.j1).and.(j1.le.JUB ))) THEN
          p2=REAL(i2-i1,r8)*(Xobs(iobs)-REAL(i1,r8))
          q2=REAL(j2-j1,r8)*(Yobs(iobs)-REAL(j1,r8))
          p1=1.0_r8-p2
          q1=1.0_r8-q2
          w11=p1*q1
          w21=p2*q1
          w22=p2*q2
          w12=p1*q2
          IF (Zobs(iobs).ge.0.0_r8) THEN
            k1=INT(Zobs(iobs))
            k2=MIN(k1+1,N(ng))
            r2=REAL(k2-k1,r8)*(Zobs(iobs)-REAL(k1,r8))
            r1=1.0_r8-r2
          ELSE
            Ztop=Adepth(i1,j1,N(ng))
            Zbot=Adepth(i1,j1,1    )              
            IF (Zobs(iobs).ge.Ztop) THEN
              k1=N(ng)
              k2=k1
              r1=1.0_r8
              r2=0.0_r8
            ELSE IF (Zbot.ge.Zobs(iobs)) THEN
              k1=1
              k2=k1
              r1=1.0_r8
              r2=0.0_r8
            ELSE                    
              DO k=N(ng),2,-1         
                Ztop=Adepth(i1,j1,k  )
                Zbot=Adepth(i1,j1,k-1)
                IF ((Ztop.gt.Zobs(iobs)).and.(Zobs(iobs).ge.Zbot)) THEN
                  k1=k-1
                  k2=k
                END IF
              END DO
              dz=Adepth(i1,j1,k2)-Adepth(i1,j1,k1)
              r2=(Zobs(iobs)-Adepth(i1,j1,k1))/dz
              r1=1.0_r8-r2
              Zobs(iobs)=REAL(k1,r8)+r2  ! overwrite to grid coordinates
            END IF
          END IF
          Hmat(1,iobs)=w11*r1
          Hmat(2,iobs)=w21*r1
          Hmat(3,iobs)=w22*r1
          Hmat(4,iobs)=w12*r1
          Hmat(5,iobs)=w11*r2
          Hmat(6,iobs)=w21*r2
          Hmat(7,iobs)=w22*r2
          Hmat(8,iobs)=w12*r2
#  ifdef MASKING
          Hmat(1,iobs)=Hmat(1,iobs)*Amask(i1,j1)
          Hmat(2,iobs)=Hmat(2,iobs)*Amask(i2,j1)
          Hmat(3,iobs)=Hmat(3,iobs)*Amask(i2,j2)
          Hmat(4,iobs)=Hmat(4,iobs)*Amask(i1,j2)
          Hmat(5,iobs)=Hmat(5,iobs)*Amask(i1,j1)
          Hmat(6,iobs)=Hmat(6,iobs)*Amask(i2,j1)
          Hmat(7,iobs)=Hmat(7,iobs)*Amask(i2,j2)
          Hmat(8,iobs)=Hmat(8,iobs)*Amask(i1,j2)
          wsum=0.0_r8
          DO ic=1,8
            wsum=wsum+Hmat(ic,iobs)
          END DO
          IF (wsum.gt.0.0_r8) THEN
            wsum=1.0_r8/wsum
            DO ic=1,8
              Hmat(ic,iobs)=Hmat(ic,iobs)*wsum
            END DO
          END IF          
#  endif
          Aobs(iobs)=Hmat(1,iobs)*A(i1,j1,k1)+                          &
     &               Hmat(2,iobs)*A(i2,j1,k1)+                          &
     &               Hmat(3,iobs)*A(i2,j2,k1)+                          &
     &               Hmat(4,iobs)*A(i1,j2,k1)+                          &
     &               Hmat(5,iobs)*A(i1,j1,k2)+                          &
     &               Hmat(6,iobs)*A(i2,j1,k2)+                          &
     &               Hmat(7,iobs)*A(i2,j2,k2)+                          &
     &               Hmat(8,iobs)*A(i1,j2,k2)
        END IF
      END DO
      RETURN

      END SUBROUTINE extract_obs3d
# endif
#endif
      END MODULE extract_obs_mod
