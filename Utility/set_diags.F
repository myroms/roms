#include "cppdefs.h"
#if defined DIAGNOSTICS_TS || defined DIAGNOSTICS_UV
      SUBROUTINE set_diags (ng, tile)
!
!=======================================================================
!  Copyright (c) 2005 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This subroutine accumulates and computes output time-averaged       !
!  diagnostic fields.  Due to synchronization, the time-averaged       !
!  diagnostic fields are computed in delayed mode. All averages        !
!  are accumulated at the beginning of the next time-step.             !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_scalars
      USE mod_stepping
!
      implicit none

      integer, intent(in) :: ng, tile

# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 5)
# endif
      CALL set_diags_tile (ng, Istr, Iend, Jstr, Jend,                  &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         KOUT, NOUT)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 5)
# endif
      RETURN
      END SUBROUTINE set_diags
!
!***********************************************************************
      SUBROUTINE set_diags_tile (ng, Istr, Iend, Jstr, Jend,            &
     &                               LBi, UBi, LBj, UBj,                &
     &                               Kout, Nout)
!***********************************************************************
!
      USE mod_param
      USE mod_diags
      USE mod_grid
      USE mod_scalars
!
      USE bc_2d_mod
# ifdef SOLVE3D
      USE bc_3d_mod
# endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Kout, Nout
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, it, j, k
      integer :: idiag

      real(r8) :: fac

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
! Initialize time-averaged diagnostic arrays when appropriate.  Notice
! that fields are initilized twice during re-start.  However, the time-
! averaged fields are computed correctly.
!-----------------------------------------------------------------------
!
      IF (((iic(ng).gt.ntsDIA(ng)).and.                                 &
     &     (MOD(iic(ng)-1,nDIA(ng)).eq.1)).or.                          &
     &    ((nrrec.gt.0).and.(iic(ng).eq.ntstart))) THEN
# ifdef DIAGNOSTICS_TS
        DO idiag=1,NDT
          DO it=1,NT(ng)
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  DIAGS(ng)%DiaTrc(i,j,k,it,idiag)=                     &
     &                      DIAGS(ng)%DiaTwrk(i,j,k,it,idiag)
                END DO
              END DO
            END DO
          END DO
        END DO
# endif
# ifdef DIAGNOSTICS_UV
        DO j=JstrR,JendR
          DO idiag=1,NDM2d
            DO i=IstrR,IendR
              DIAGS(ng)%DiaU2d(i,j,idiag)=DIAGS(ng)%DiaU2wrk(i,j,idiag)
              DIAGS(ng)%DiaV2d(i,j,idiag)=DIAGS(ng)%DiaV2wrk(i,j,idiag)
            END DO
          END DO
#  ifdef SOLVE3D
          DO idiag=1,NDM3d
            DO k=1,N(ng)
              DO i=IstrR,IendR
                DIAGS(ng)%DiaU3d(i,j,k,idiag)=                          &
     &                    DIAGS(ng)%DiaU3wrk(i,j,k,idiag)
                DIAGS(ng)%DiaV3d(i,j,k,idiag)=                          &
     &                    DIAGS(ng)%DiaV3wrk(i,j,k,idiag)
              END DO
            END DO
          END DO
#  endif
        END DO
# endif
!
!-----------------------------------------------------------------------
!  Accumulate time-averaged fields.
!-----------------------------------------------------------------------
!
      ELSE IF (iic(ng).gt.ntsDIA(ng)) THEN
# ifdef DIAGNOSTICS_TS
        DO idiag=1,NDT
          DO it=1,NT(ng)
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  DIAGS(ng)%DiaTrc(i,j,k,it,idiag)=                     &
     &                      DIAGS(ng)%DiaTrc(i,j,k,it,idiag)+           &
     &                      DIAGS(ng)%DiaTwrk(i,j,k,it,idiag)
                END DO
              END DO
            END DO
          END DO
        END DO
# endif
# ifdef DIAGNOSTICS_UV
        DO j=JstrR,JendR
          DO idiag=1,NDM2d
            DO i=IstrR,IendR
              DIAGS(ng)%DiaU2d(i,j,idiag)=DIAGS(ng)%DiaU2d(i,j,idiag)+  &
     &                                    DIAGS(ng)%DiaU2wrk(i,j,idiag)
              DIAGS(ng)%DiaV2d(i,j,idiag)=DIAGS(ng)%DiaV2d(i,j,idiag)+  &
     &                                    DIAGS(ng)%DiaV2wrk(i,j,idiag)
            END DO
          END DO
#  ifdef SOLVE3D
          DO idiag=1,NDM3d
            DO k=1,N(ng)
              DO i=IstrR,IendR
                DIAGS(ng)%DiaU3d(i,j,k,idiag)=                          &
     &                    DIAGS(ng)%DiaU3d(i,j,k,idiag)+                &
     &                    DIAGS(ng)%DiaU3wrk(i,j,k,idiag)
                DIAGS(ng)%DiaV3d(i,j,k,idiag)=                          &
     &                    DIAGS(ng)%DiaV3d(i,j,k,idiag)+                &
     &                    DIAGS(ng)%DiaV3wrk(i,j,k,idiag)
              END DO
            END DO
          END DO
#  endif
        END DO
# endif
      END IF
!
!-----------------------------------------------------------------------
!  Apply periodic or gradient boundary conditions and land-mask
!  for output purposes.
!-----------------------------------------------------------------------
!
      IF ((iic(ng).gt.ntsDIA(ng)).and.                                  &
     &    (MOD(iic(ng)-1,nDIA(ng)).eq.0).and.                           &
     &    ((iic(ng).ne.ntstart).or.(nrrec.eq.0))) THEN
# ifdef DIAGNOSTICS_TS
!
!  Apply periodic or gradient boundary conditions for output purposes.
!
        DO idiag=1,NDT
          DO it=1,NT(ng)
            CALL bc_r3d_tile (ng, iNLM, Istr, Iend, Jstr, Jend,         &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        NghostPoints,                             &
     &                        DIAGS(ng)%DiaTrc(:,:,:,it,idiag))
#  ifdef MASKING
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  DIAGS(ng)%DiaTrc(i,j,k,it,idiag)=                     &
     &                      DIAGS(ng)%DiaTrc(i,j,k,it,idiag)*           &
     &                      GRID(ng)%rmask(i,j)
                END DO
              END DO
            END DO
#  endif
          END DO
        END DO
# endif
# ifdef DIAGNOSTICS_UV
!
!  Apply periodic or gradient boundary conditions for output purposes.
!
        DO idiag=1,NDM2d
          CALL bc_u2d_tile (ng, iNLM, Istr, Iend, Jstr, Jend,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints,                               &
     &                      DIAGS(ng)%DiaU2d(:,:,idiag))
          CALL bc_v2d_tile (ng, iNLM, Istr, Iend, Jstr, Jend,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints,                               &
     &                      DIAGS(ng)%DiaV2d(:,:,idiag))
        END DO
#  ifdef SOLVE3D
        DO idiag=1,NDM3d
          CALL bc_u3d_tile (ng, iNLM, Istr, Iend, Jstr, Jend,           &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      NghostPoints,                               &
     &                      DIAGS(ng)%DiaU3d(:,:,:,idiag))
          CALL bc_v3d_tile (ng, iNLM, Istr, Iend, Jstr, Jend,           &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      NghostPoints,                               &
     &                      DIAGS(ng)%DiaV3d(:,:,:,idiag))
        END DO
#  endif
#  ifdef MASKING
        DO j=JstrR,JendR
          DO idiag=1,NDM2d
            DO i=IstrR,IendR
              DIAGS(ng)%DiaU2d(i,j,idiag)=DIAGS(ng)%DiaU2d(i,j,idiag)*  &
     &                                    GRID(ng)%umask(i,j)
              DIAGS(ng)%DiaV2d(i,j,idiag)=DIAGS(ng)%DiaV2d(i,j,idiag)*  &
     &                                    GRID(ng)%vmask(i,j)
            END DO
          END DO
#   ifdef SOLVE3D
          DO idiag=1,NDM3d
            DO k=1,N(ng)
              DO i=IstrR,IendR
                DIAGS(ng)%DiaU3d(i,j,k,idiag)=                          &
     &                                 DIAGS(ng)%DiaU3d(i,j,k,idiag)*   &
     &                                 GRID(ng)%umask(i,j)
                DIAGS(ng)%DiaV3d(i,j,k,idiag)=                          &
     &                                 DIAGS(ng)%DiaV3d(i,j,k,idiag)*   &
     &                                 GRID(ng)%vmask(i,j)
              END DO
            END DO
          END DO
#   endif
        END DO
#  endif
# endif
      END IF
      RETURN
      END SUBROUTINE set_diags_tile
#else
      SUBROUTINE set_diags
      RETURN
      END SUBROUTINE set_diags
#endif
