#include "cppdefs.h"
      MODULE step3d_t_mod
#ifdef SOLVE3D
!
!========================================== Alexander F. Shchepetkin ===
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine time-steps tracer equations. Notice that avective      !
!  and diffusive terms are time-steped differently.                    !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC step3d_t

      CONTAINS
!
!***********************************************************************
      SUBROUTINE step3d_t (ng, tile)
!***********************************************************************
!
      USE mod_param
# ifdef CLIMATOLOGY
      USE mod_clima
# endif
# ifdef BIOLOGY
      USE mod_forces
# endif
      USE mod_grid
      USE mod_mixing
# if defined ASSIMILATION || defined NUDGING
      USE mod_obs
# endif
      USE mod_ocean
# ifdef TS_PSOURCE
      USE mod_sources
# endif
      USE mod_stepping
!
      integer, intent(in) :: ng, tile

# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, 35)
# endif
      CALL step3d_t_tile (ng, Istr, Iend, Jstr, Jend,                   &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    nrhs(ng), nstp(ng), nnew(ng),                 &
# ifdef TS_PSOURCE
     &                    Nsrc(ng),                                     &
     &                    SOURCES(ng) % Lsrc,                           &
     &                    SOURCES(ng) % Isrc,                           &
     &                    SOURCES(ng) % Jsrc,                           &
     &                    SOURCES(ng) % Dsrc,                           &
     &                    SOURCES(ng) % Tsrc,                           &
# endif
# ifdef MASKING
     &                    GRID(ng) % rmask,                             &
     &                    GRID(ng) % umask,                             &
     &                    GRID(ng) % vmask,                             &
# endif
     &                    GRID(ng) % Hz,                                &
     &                    GRID(ng) % Huon,                              &
     &                    GRID(ng) % Hvom,                              &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
     &                    GRID(ng) % z_r,                               &
# ifdef BIOLOGY
     &                    GRID(ng) % z_w,                               &
     &                    FORCES(ng) % srflx,                           &
# endif
# if defined TCLM_NUDGING && defined TCLIMATOLOGY
     &                    CLIMA(ng) % Tnudgcof,                         &
     &                    CLIMA(ng) % tclm,                             &
# endif
# if defined NUDGING_SST || defined NUDGING_T
     &                    OBS(ng) % EobsT,                              &
     &                    OBS(ng) % Tobs,                               &
# endif
     &                    MIXING(ng) % Akt,                             &
     &                    OCEAN(ng) % W,                                &
# if defined FLOATS && defined FLOAT_VWALK
     &                    MIXING(ng) % dAktdz,                          &
# endif
     &                    OCEAN(ng) % t)
# ifdef PROFILE
      CALL wclock_off (ng, 35)
# endif
      RETURN
      END SUBROUTINE step3d_t
!
!***********************************************************************
      SUBROUTINE step3d_t_tile (ng, Istr, Iend, Jstr, Jend,             &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          nrhs, nstp, nnew,                       &
# ifdef TS_PSOURCE
     &                          Nsrc, Lsrc, Isrc, Jsrc, Dsrc, Tsrc,     &
# endif
# ifdef MASKING
     &                          rmask, umask, vmask,                    &
# endif
     &                          Hz, Huon, Hvom, pm, pn, z_r,            &
# ifdef BIOLOGY
     &                          z_w, srflx,                             &
# endif
# if defined TCLM_NUDGING && defined TCLIMATOLOGY
     &                          Tnudgcof, tclm,                         &
# endif
# if defined NUDGING_SST || defined NUDGING_T
     &                          EobsT, Tobs,                            &
# endif
     &                          Akt, W,                                 &
# if defined FLOATS && defined FLOAT_VWALK
     &                          dAktdz,                                 &
# endif
     &                          t)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
!
# ifdef BIOLOGY
      USE biology_mod, ONLY : biology_tile
# endif
# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      USE exchange_3d_mod, ONLY : exchange_r3d_tile
# endif
      USE t3dbc_mod, ONLY : t3dbc_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: nrhs, nstp, nnew

# ifdef ASSUMED_SHAPE
#  ifdef TS_PSOURCE
      integer, intent(in) :: Nsrc

      logical, intent(in) :: Lsrc(:,:)

      integer, intent(in) :: Isrc(:)
      integer, intent(in) :: Jsrc(:)

      real(r8), intent(in) :: Dsrc(:)
      real(r8), intent(in) :: Tsrc(:,:,:)
#  endif
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: Huon(LBi:,LBj:,:)
      real(r8), intent(in) :: Hvom(LBi:,LBj:,:)
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
#  ifdef BIOLOGY
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: srflx(LBi:,LBj:)
#  endif
#  if defined TCLM_NUDGING && defined TCLIMATOLOGY
      real(r8), intent(in) :: Tnudgcof(LBi:,LBj:,:)
      real(r8), intent(in) :: tclm(LBi:,LBj:,:,:)
#  endif
#  if defined NUDGING_SST || defined NUDGING_T
      real(r8), intent(in) :: EobsT(LBi:,LBj:,:,:)
      real(r8), intent(in) :: Tobs(LBi:,LBj:,:,:)
#  endif
      real(r8), intent(in) :: Akt(LBi:,LBj:,0:,:)
!!    real(r8), intent(in) :: Akt(LBi:UBi,LBj:UBj,0:N(ng),NAT)     ! explicit
      real(r8), intent(in) :: W(LBi:,LBj:,0:)

      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
!!    real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng)) ! explicit

#  if defined FLOATS && defined FLOAT_VWALK
      real(r8), intent(out) :: dAktdz(LBi:,LBj:,:)
#  endif

# else

#  ifdef TS_PSOURCE
      integer, intent(in) :: Nsrc

      logical, intent(in) :: Lsrc(Msrc,NT(ng))

      integer, intent(in) :: Isrc(Msrc)
      integer, intent(in) :: Jsrc(Msrc)

      real(r8), intent(in) :: Dsrc(Msrc)
      real(r8), intent(in) :: Tsrc(Msrc,N(ng),NT(ng))
#  endif

#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: Huon(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: Hvom(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
#  ifdef BIOLOGY
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(in) :: srflx(LBi:UBi,LBj:UBj)
#  endif
#  if defined TCLM_NUDGING && defined TCLIMATOLOGY
      real(r8), intent(in) :: Tnudgcof(LBi:UBi,LBj:UBj,NT(ng))
      real(r8), intent(in) :: tclm(LBi:UBi,LBj:UBj,N(ng),NT(ng))
#  endif
#  if defined NUDGING_SST || defined NUDGING_T
      real(r8), intent(in) :: EobsT(LBi:UBi,LBj:UBj,N(ng),NT(ng))
      real(r8), intent(in) :: Tobs(LBi:UBi,LBj:UBj,N(ng),NT(ng))
#  endif
      real(r8), intent(in) :: Akt(LBi:UBi,LBj:UBj,0:N(ng),NAT)
      real(r8), intent(in) :: W(LBi:UBi,LBj:UBj,0:N(ng))

      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))

#  if defined FLOATS && defined FLOAT_VWALK
      real(r8), intent(out) :: dAktdz(LBi:UBi,LBj:UBj,N(ng))
#  endif
# endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, is, itrc, j, k, ltrc

      real(r8), parameter :: eps = 1.0E-16_r8

      real(r8) :: cff, cff1, cff2, cff3

      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N(ng)) :: CF
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N(ng)) :: BC
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N(ng)) :: DC
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N(ng)) :: FC
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N(ng)) :: oHz

      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: FE
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: FX


      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: curv
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: grad

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Time-step horizontal advection term.
!-----------------------------------------------------------------------
!
      DO itrc=1,NT(ng)
        DO k=1,N(ng)
# ifdef TS_C2HADVECTION
          DO j=Jstr,Jend
            DO i=Istr,Iend+1
              FX(i,j)=Huon(i,j,k)*                                      &
     &                0.5_r8*(t(i-1,j,k,3,itrc)+t(i,j,k,3,itrc))
            END DO
          END DO
          DO j=Jstr,Jend+1
            DO i=Istr,Iend
              FE(i,j)=Hvom(i,j,k)*                                      &
     &                0.5_r8*(t(i,j-1,k,3,itrc)+t(i,j,k,3,itrc))
            END DO
          END DO
# else
#  ifdef EW_PERIODIC
#   define I_RANGE Istr-1,Iend+2
#  else
#   define I_RANGE MAX(Istr-1,1),MIN(Iend+2,Lm(ng)+1)
#  endif
          DO j=Jstr,Jend
            DO i=I_RANGE
              FX(i,j)=(t(i,j,k,3,itrc)-t(i-1,j,k,3,itrc))
#  ifdef MASKING
              FX(i,j)=FX(i,j)*umask(i,j)
#  endif
            END DO
          END DO
#  undef I_RANGE
#  ifndef EW_PERIODIC
          IF (WESTERN_EDGE) THEN
            DO j=Jstr,Jend
              FX(Istr-1,j)=FX(Istr,j)
            END DO
          END IF
          IF (EASTERN_EDGE) THEN
            DO j=Jstr,Jend
              FX(Iend+2,j)=FX(Iend+1,j)
            END DO
          END IF
#  endif
!
          DO j=Jstr,Jend
            DO i=Istr-1,Iend+1
#  if defined TS_U3HADVECTION
              curv(i,j)=FX(i+1,j)-FX(i,j)
#  elif defined TS_A4HADVECTION
              cff=2.0_r8*FX(i+1,j)*FX(i,j)
              IF (cff.gt.eps) THEN
                grad(i,j)=cff/(FX(i+1,j)+FX(i,j))
              ELSE
                grad(i,j)=0.0_r8
              END IF
#  else
              grad(i,j)=0.5_r8*(FX(i+1,j)+FX(i,j))
#  endif
            END DO
          END DO
!
          cff1=1.0_r8/6.0_r8
          cff2=1.0_r8/3.0_r8
          DO j=Jstr,Jend
            DO i=Istr,Iend+1
#  ifdef TS_U3HADVECTION
              FX(i,j)=Huon(i,j,k)*0.5_r8*                               &
     &                (t(i-1,j,k,3,itrc)+t(i,j,k,3,itrc))-              &
     &                cff1*(curv(i-1,j)*MAX(Huon(i,j,k),0.0_r8)+        &
     &                      curv(i  ,j)*MIN(Huon(i,j,k),0.0_r8))
#  else
              FX(i,j)=Huon(i,j,k)*0.5_r8*                               &
     &                (t(i-1,j,k,3,itrc)+t(i,j,k,3,itrc)-               &
     &                 cff2*(grad(i,j)-grad(i-1,j)))
#  endif
            END DO
          END DO
!
#  ifdef NS_PERIODIC
#   define J_RANGE Jstr-1,Jend+2
#  else
#   define J_RANGE MAX(Jstr-1,1),MIN(Jend+2,Mm(ng)+1)
#  endif
          DO j=J_RANGE
            DO i=Istr,Iend
              FE(i,j)=(t(i,j,k,3,itrc)-t(i,j-1,k,3,itrc))
#  ifdef MASKING
              FE(i,j)=FE(i,j)*vmask(i,j)
#  endif
            END DO
          END DO
#  undef J_RANGE
#  ifndef NS_PERIODIC
          IF (SOUTHERN_EDGE) THEN
            DO i=Istr,Iend
              FE(i,Jstr-1)=FE(i,Jstr)
            END DO
          END IF
          IF (NORTHERN_EDGE) THEN
            DO i=Istr,Iend
              FE(i,Jend+2)=FE(i,Jend+1)
            END DO
          END IF
#  endif
!
          DO j=Jstr-1,Jend+1
            DO i=Istr,Iend
#  if defined TS_U3HADVECTION
              curv(i,j)=FE(i,j+1)-FE(i,j)
#  elif defined TS_A4HADVECTION
              cff=2.0_r8*FE(i,j+1)*FE(i,j)
              IF (cff.gt.eps) THEN
                grad(i,j)=cff/(FE(i,j+1)+FE(i,j))
              ELSE
                grad(i,j)=0.0_r8
              END IF
#  else
              grad(i,j)=0.5_r8*(FE(i,j+1)+FE(i,j))
#  endif
            END DO
          END DO
!
          cff1=1.0_r8/6.0_r8
          cff2=1.0_r8/3.0_r8
          DO j=Jstr,Jend+1
            DO i=Istr,Iend
#  ifdef TS_U3HADVECTION
              FE(i,j)=Hvom(i,j,k)*0.5_r8*                               &
     &                (t(i,j-1,k,3,itrc)+t(i,j,k,3,itrc))-              &
     &                cff1*(curv(i,j-1)*MAX(Hvom(i,j,k),0.0_r8)+        &
     &                      curv(i,j  )*MIN(Hvom(i,j,k),0.0_r8))
#  else
              FE(i,j)=Hvom(i,j,k)*0.5_r8*                               &
     &                (t(i,j-1,k,3,itrc)+t(i,j,k,3,itrc)-               &
     &                 cff2*(grad(i,j)-grad(i,j-1)))
#  endif
            END DO
          END DO
# endif
# ifdef TS_PSOURCE
!
!  Apply tracers point sources to the horizontal advection terms.
!
          DO is=1,Nsrc
            i=Isrc(is)
            j=Jsrc(is)
            IF (INT(Dsrc(is)).eq.0) THEN
              IF (((Istr.le.i).and.(i.le.Iend+1)).and.                  &
     &            ((Jstr.le.j).and.(j.le.Jend))) THEN
                IF (Lsrc(is,itrc)) THEN
                  FX(i,j)=Huon(i,j,k)*Tsrc(is,k,itrc)
#  ifdef MASKING
                ELSE
                  IF ((rmask(i  ,j).eq.0.0_r8).and.                     &
     &                (rmask(i-1,j).eq.1.0_r8)) THEN
                    FX(i,j)=Huon(i,j,k)*t(i-1,j,k,3,itrc)
                  ELSE IF ((rmask(i  ,j).eq.1.0_r8).and.                &
     &                     (rmask(i-1,j).eq.0.0_r8)) THEN
                    FX(i,j)=Huon(i,j,k)*t(i  ,j,k,3,itrc)
                  END IF
#  endif
                END IF
              END IF
            ELSE IF (INT(Dsrc(is)).eq.1) THEN
              IF (((Istr.le.i).and.(i.le.Iend)).and.                    &
     &            ((Jstr.le.j).and.(j.le.Jend+1))) THEN
                IF (Lsrc(is,itrc)) THEN
                  FE(i,j)=Hvom(i,j,k)*Tsrc(is,k,itrc)
#  ifdef MASKING
                ELSE
                  IF ((rmask(i,j  ).eq.0.0_r8).and.                     &
     &                (rmask(i,j-1).eq.1.0_r8)) THEN
                    FE(i,j)=Hvom(i,j,k)*t(i,j-1,k,3,itrc)
                  ELSE IF ((rmask(i,j  ).eq.1.0_r8).and.                &
     &                     (rmask(i,j-1).eq.0.0_r8)) THEN
                    FE(i,j)=Hvom(i,j,k)*t(i,j  ,k,3,itrc)
                  END IF
#  endif
                END IF
              END IF
            END IF
          END DO
# endif
!
!  Time-step horizontal advection term.
!
          DO j=Jstr,Jend
            DO i=Istr,Iend
              t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)-                    &
     &                           dt(ng)*(FX(i+1,j)-FX(i,j)+             &
     &                                   FE(i,j+1)-FE(i,j))*            &
     &                           pm(i,j)*pn(i,j)
            END DO
          END DO
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Time-step vertical advection term.
!-----------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO itrc=1,NT(ng)
          ltrc=MIN(NAT,itrc)
# ifdef TS_SVADVECTION
!
!  Build conservative parabolic splines for the vertical derivatives
!  "FC" of the tracer.  Then, the interfacial "FC" values are
!  converted to vertical advective fluxes.
!
            DO i=Istr,Iend
#  ifdef NEUMANN
              FC(i,0)=1.5_r8*t(i,j,1,3,itrc)
              CF(i,1)=0.5_r8
#  else
              FC(i,0)=2.0_r8*t(i,j,1,3,itrc)
              CF(i,1)=1.0_r8
#  endif
            END DO
            DO k=1,N(ng)-1
              DO i=Istr,Iend
                cff=1.0_r8/(2.0_r8*Hz(i,j,k)+                           &
     &                      Hz(i,j,k+1)*(2.0_r8-CF(i,k)))
                CF(i,k+1)=cff*Hz(i,j,k)
                FC(i,k)=cff*(3.0_r8*(Hz(i,j,k  )*t(i,j,k+1,3,itrc)+     &
     &                               Hz(i,j,k+1)*t(i,j,k  ,3,itrc))-    &
     &                       Hz(i,j,k+1)*FC(i,k-1))
              END DO
            END DO
            DO i=Istr,Iend
#  ifdef NEUMANN
              FC(i,N(ng))=(3.0_r8*t(i,j,N(ng),3,itrc)-FC(i,N(ng)-1))/   &
     &                    (2.0_r8-CF(i,N(ng)))
#  else
              FC(i,N(ng))=(2.0_r8*t(i,j,N(ng),3,itrc)-FC(i,N(ng)-1))/   &
     &                    (1.0_r8-CF(i,N(ng)))
#  endif
            END DO
            DO k=N(ng)-1,0,-1
              DO i=Istr,Iend
                FC(i,k)=FC(i,k)-CF(i,k+1)*FC(i,k+1)
                FC(i,k+1)=W(i,j,k+1)*FC(i,k+1)
              END DO
            END DO
            DO i=Istr,Iend
              FC(i,N(ng))=0.0_r8
              FC(i,0)=0.0_r8
            END DO
# elif defined TS_A4VADVECTION
            DO k=1,N(ng)-1
              DO i=Istr,Iend
                FC(i,k)=t(i,j,k+1,3,itrc)-t(i,j,k,3,itrc)
              END DO
            END DO
            DO i=Istr,Iend
              FC(i,0)=FC(i,1)
              FC(i,N(ng))=FC(i,N(ng)-1)
            END DO
            DO k=1,N(ng)
              DO i=Istr,Iend
                cff=2.0_r8*FC(i,k)*FC(i,k-1)
                IF (cff.gt.eps) THEN
                  CF(i,k)=cff/(FC(i,k)+FC(i,k-1))
                ELSE
                  CF(i,k)=0.0_r8
                END IF
              END DO
            END DO
            cff1=1.0_r8/3.0_r8
            DO k=1,N(ng)-1
              DO i=Istr,Iend
                FC(i,k)=W(i,j,k)*                                       &
     &                  0.5_r8*(t(i,j,k  ,3,itrc)+                      &
     &                          t(i,j,k+1,3,itrc)-                      &
     &                          cff1*(CF(i,k+1)-CF(i,k)))
              END DO
            END DO
            DO i=Istr,Iend
              FC(i,0)=0.0_r8
              FC(i,N(ng))=0.0_r8
            END DO
# elif defined TS_C2VADVECTION
            DO k=1,N(ng)-1
              DO i=Istr,Iend
                FC(i,k)=W(i,j,k)*                                       &
     &                  0.5_r8*(t(i,j,k  ,3,itrc)+                      &
     &                          t(i,j,k+1,3,itrc))
              END DO
            END DO
            DO i=Istr,Iend
              FC(i,0)=0.0_r8
              FC(i,N(ng))=0.0_r8
            END DO
# else
            cff1=0.5_r8
            cff2=7.0_r8/12.0_r8
            cff3=1.0_r8/12.0_r8
            DO k=2,N(ng)-2
              DO i=Istr,Iend
                FC(i,k)=W(i,j,k)*(cff2*(t(i,j,k  ,3,itrc)+              &
     &                                  t(i,j,k+1,3,itrc))-             &
     &                            cff3*(t(i,j,k-1,3,itrc)+              &
     &                                  t(i,j,k+2,3,itrc)))
              END DO
            END DO
            DO i=Istr,Iend
              FC(i,0)=0.0_r8
              FC(i,1)=W(i,j,1)*(cff1*t(i,j,1,3,itrc)+                   &
     &                          cff2*t(i,j,2,3,itrc)-                   &
     &                          cff3*t(i,j,3,3,itrc))
              FC(i,N(ng)-1)=W(i,j,N(ng)-1)*                             &
     &                      (cff1*t(i,j,N(ng)  ,3,itrc)+                &
     &                       cff2*t(i,j,N(ng)-1,3,itrc)-                &
     &                       cff3*t(i,j,N(ng)-2,3,itrc))
              FC(i,N(ng))=0.0_r8
            END DO
# endif
!
!  Time-step vertical advection term.
!
          DO i=Istr,Iend
            CF(i,0)=dt(ng)*pm(i,j)*pn(i,j)
          END DO
          DO k=1,N(ng)
            DO i=Istr,Iend
              t(i,j,k,nnew,itrc)=(t(i,j,k,nnew,itrc)-                   &
     &                            CF(i,0)*(FC(i,k)-FC(i,k-1)))
# ifdef SPLINES
              t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)/Hz(i,j,k)
# endif
            END DO
          END DO
!
!-----------------------------------------------------------------------
!  Time-step vertical diffusion term.
!-----------------------------------------------------------------------
# ifdef SPLINES
!
!  Use conservative, parabolic spline reconstruction of vertical
!  diffusion derivatives.  Then, time step vertical diffusion term
!  implicitly.
!
          DO k=1,N(ng)
            DO i=Istr,Iend
              oHz(i,k)=1.0_r8/Hz(i,j,k)
            END DO
          END DO
!
          cff1=1.0_r8/6.0_r8
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              FC(i,k)=cff1*Hz(i,j,k  )-                                 &
     &                dt(ng)*Akt(i,j,k-1,ltrc)*oHz(i,k  )
              CF(i,k)=cff1*Hz(i,j,k+1)-                                 &
     &                dt(ng)*Akt(i,j,k+1,ltrc)*oHz(i,k+1)
            END DO
          END DO
          DO i=Istr,Iend
            CF(i,0)=0.0_r8
            DC(i,0)=0.0_r8
          END DO
!
          cff1=1.0_r8/3.0_r8
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              BC(i,k)=cff1*(Hz(i,j,k)+Hz(i,j,k+1))+                     &
     &                dt(ng)*Akt(i,j,k,ltrc)*(oHz(i,k)+oHz(i,k+1))
              cff=1.0_r8/(BC(i,k)-FC(i,k)*CF(i,k-1))
              CF(i,k)=cff*CF(i,k)
              DC(i,k)=cff*(t(i,j,k+1,nnew,itrc)-t(i,j,k,nnew,itrc)-     &
     &                     FC(i,k)*DC(i,k-1))
            END DO
          END DO
!
          DO i=Istr,Iend
            DC(i,N(ng))=0.0_r8
          END DO
          DO k=N(ng)-1,1,-1
            DO i=Istr,Iend
              DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1)
            END DO
          END DO
!
          DO k=1,N(ng)
            DO i=Istr,Iend
              DC(i,k)=DC(i,k)*Akt(i,j,k,ltrc)
              t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)+                    &
     &                           dt(ng)*oHz(i,k)*(DC(i,k)-DC(i,k-1))
            END DO
          END DO
# else
!
!  Compute off-diagonal coefficients FC [lambda*dt*Akt/Hz] for the
!  implicit vertical diffusion terms at future time step, located
!  at horizontal RHO-points and vertical W-points.
!  Also set FC at the top and bottom levels.
!
          cff=-dt(ng)*lambda
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              FC(i,k)=cff*Akt(i,j,k,ltrc)/(z_r(i,j,k+1)-z_r(i,j,k))
            END DO
          END DO
          DO i=Istr,Iend
            FC(i,0)=0.0_r8
            FC(i,N(ng))=0.0_r8
          END DO
!
!  Compute diagonal matrix coefficients BC and load right-hand-side
!  terms for the tracer equation into DC.
!
          DO k=1,N(ng)
            DO i=Istr,Iend
              BC(i,k)=Hz(i,j,k)-FC(i,k)-FC(i,k-1)
              DC(i,k)=t(i,j,k,nnew,itrc)
            END DO
          END DO
!
!  Solve the tridiagonal system.
!
          DO i=Istr,Iend
            cff=1.0_r8/BC(i,1)
            CF(i,1)=cff*FC(i,1)
            DC(i,1)=cff*DC(i,1)
          END DO
          DO k=2,N(ng)-1
            DO i=Istr,Iend
              cff=1.0_r8/(BC(i,k)-FC(i,k-1)*CF(i,k-1))
              CF(i,k)=cff*FC(i,k)
              DC(i,k)=cff*(DC(i,k)-FC(i,k-1)*DC(i,k-1))
            END DO
          END DO
          DO i=Istr,Iend
             DC(i,N(ng))=(DC(i,N(ng))-FC(i,N(ng)-1)*DC(i,N(ng)-1))/     &
     &                   (BC(i,N(ng))-FC(i,N(ng)-1)*CF(i,N(ng)-1))
             t(i,j,N(ng),nnew,itrc)=DC(i,N(ng))
          END DO
          DO k=N(ng)-1,1,-1
            DO i=Istr,Iend
              DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1)
              t(i,j,k,nnew,itrc)=DC(i,k)
            END DO
          END DO
# endif
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Apply lateral boundary conditions and, if appropriate, nudge
!  to tracer data and apply Land/Sea mask.
!-----------------------------------------------------------------------
!
      DO itrc=1,NT(ng)
!
!  Set lateral boundary conditions.
!
        CALL t3dbc_tile (ng, Istr, Iend, Jstr, Jend, itrc,              &
     &                   LBi, UBi, LBj, UBj, N(ng), NT(ng),             &
     &                   nnew, nstp,                                    &
     &                   t)
# ifdef BIOLOGY
      END DO
!
!  Apply interactive part of the biological module.
!
      CALL biology_tile (ng, Istr, Iend, Jstr, Jend,                    &
     &                   LBi, UBi, LBj, UBj, N(ng), NT(ng),             &
     &                   nnew,                                          &
#  ifdef MASKING
     &                   rmask,                                         &
#  endif
     &                   Hz, z_r, z_w, srflx, t)
      DO itrc=1,NT(ng)
# endif
# if defined TCLM_NUDGING && defined TCLIMATOLOGY
!
!  Nudge towards tracer climatology.
!
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)+                    &
     &                           dt(ng)*Tnudgcof(i,j,itrc)*             &
     &                           (tclm(i,j,k,itrc)-t(i,j,k,nnew,itrc))
            END DO
          END DO
        END DO
# endif
# if defined NUDGING_SST || defined NUDGING_T
!
!  Assimilate tracer observations via nudging.
!
        IF (update_T(itrc,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                cff=MIN(1.0_r8,MAX(0.0_r8,EobsT(i,j,k,itrc)))
                cff=dt(ng)*Tnudass(itrc,ng)*(1.0_r8-cff)
                t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)+                  &
     &                             cff*(Tobs(i,j,k,itrc)-               &
     &                                  t(i,j,k,nnew,itrc))
              END DO
            END DO
          END DO
        END IF
# endif
# ifdef MASKING
!
!  Apply Land/Sea mask.
!
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)*rmask(i,j)
            END DO
          END DO
        END DO
# endif
# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
!
!  Apply periodic boundary conditions.
!
        CALL exchange_r3d_tile (ng, Istr, Iend, Jstr, Jend,             &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          t(:,:,:,nnew,itrc))
# endif
      END DO
# if defined FLOATS && defined FLOAT_VWALK
!
!  Compute vertical gradient in vertical T-diffusion coefficient for
!  floats random walk.
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          DO k=1,N(ng)
            dAktdz(i,j,k)=(Akt(i,j,k,1)-Akt(i,j,k-1,1))/Hz(i,j,k)
          END DO
        END DO
      END DO
# endif
      RETURN
      END SUBROUTINE step3d_t_tile
#endif
      END MODULE step3d_t_mod
