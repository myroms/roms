#include "cppdefs.h"
#if defined MASKING && (defined READ_WATER || defined WRITE_WATER)
      SUBROUTINE wpoints (ng, tile)
!
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!=======================================================================
!                                                                      !
!  This routine sets variables associate with reading and writing of   !
!  water points data.                                                  !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid

      implicit none

      integer, intent(in) :: ng, tile

# include "tile.h"
!
      CALL wpoints_tile (ng, Istr, Iend, Jstr, Jend,                    &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   GRID(ng) % pmask,                              &
     &                   GRID(ng) % rmask,                              &
     &                   GRID(ng) % umask,                              &
     &                   GRID(ng) % vmask)
      RETURN
      END SUBROUTINE wpoints
!
!***********************************************************************
      SUBROUTINE wpoints_tile (ng, Istr, Iend, Jstr, Jend,              &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         pmask, rmask, umask, vmask)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_ncparam
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj

# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: pmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
# else
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: NSUB, Npts, i, j, my_Nxyp, my_Nxyr, my_Nxyu, my_Nxyv

# ifdef DISTRIBUTE
      integer :: ic, ij
      real(r8), dimension(4) :: buffer
      character (len=3), dimension(4) :: op_handle
#  ifdef READ_WATER
      real(r8), dimension((Lm(ng)+2)*(Mm(ng)+2)) :: mask
#  endif
# endif

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Determine number of water points. 
!-----------------------------------------------------------------------
!
      my_Nxyr=0
      my_Nxyp=0
      my_Nxyu=0
      my_Nxyv=0
      DO j=MAX(0,JstrR),MIN(JendR,Mm(ng)+1)
        DO i=MAX(0,IstrR),MIN(IendR,Lm(ng)+1)
          IF (rmask(i,j).gt.0.0_r8) THEN
            my_Nxyr=my_Nxyr+1
          END IF
        END DO
        DO i=Istr,MIN(IendR,Lm(ng)+1)
          IF (umask(i,j).gt.0.0_r8) THEN
            my_Nxyu=my_Nxyu+1
          END IF
        END DO
      END DO
      DO j=Jstr,MIN(JendR,Mm(ng)+1) 
        DO i=Istr,MIN(IendR,Lm(ng)+1)
          IF (pmask(i,j).eq.1.0_r8) THEN
            my_Nxyp=my_Nxyp+1
          END IF
        END DO
        DO i=MAX(0,IstrR),MIN(IendR,Lm(ng)+1)
          IF (vmask(i,j).gt.0.0_r8) THEN
            my_Nxyv=my_Nxyv+1
          END IF
        END DO
      END DO
!
! Determine global maximum.
!
      IF (SOUTH_WEST_CORNER.and.                                        &
     &    NORTH_EAST_CORNER) THEN
        NSUB=1                           ! non-tiled application
      ELSE
        NSUB=NtileX(ng)*NtileE(ng)       ! tiled application
      END IF
!$OMP CRITICAL (MAX_WATER)
      IF (tile_count.eq.0) THEN
        Nxyp(ng)=0
        Nxyr(ng)=0
        Nxyu(ng)=0
        Nxyv(ng)=0
      END IF
      Nxyp(ng)=Nxyp(ng)+my_Nxyp
      Nxyr(ng)=Nxyr(ng)+my_Nxyr
      Nxyu(ng)=Nxyu(ng)+my_Nxyu
      Nxyv(ng)=Nxyv(ng)+my_Nxyv
      tile_count=tile_count+1
      IF (tile_count.eq.NSUB) THEN
        tile_count=0
# ifdef DISTRIBUTE
        buffer(1)=REAL(Nxyp(ng),r8)
        buffer(2)=REAL(Nxyr(ng),r8)
        buffer(3)=REAL(Nxyu(ng),r8)
        buffer(4)=REAL(Nxyv(ng),r8)
        op_handle(1)='SUM'
        op_handle(2)='SUM'
        op_handle(3)='SUM'
        op_handle(4)='SUM'
        CALL mp_reduce (ng, 4, buffer, op_handle)
        Nxyp(ng)=INT(buffer(1))
        Nxyr(ng)=INT(buffer(2))
        Nxyu(ng)=INT(buffer(3))
        Nxyv(ng)=INT(buffer(4))
# endif
      END IF
!$OMP END CRITICAL (MAX_WATER)

# if defined DISTRIBUTE && defined READ_WATER
!
!-----------------------------------------------------------------------
!  If distributed-memory, set water indices for IO processing.
!-----------------------------------------------------------------------
!
!  RHO-points indices.
!
      CALL mp_gather (ng, LBi, UBi, LBj, UBj, 1, 1, r2dvar, 1.0_r8,     &
     &                rmask(LBi,LBj), rmask(LBi,LBj), Npts, mask)
      CALL mp_bcastf (ng, mask, Npts)

      ic=0
      DO ij=1,Npts
        IF (mask(ij).gt.0.0_r8) THEN
          ic=ic+1
          SCALARS(ng)%IJwater(ic,r2dvar)=ij
        END IF
      END DO
      my_Nxyr=ic
!
!  PSI-points indices.
!
      CALL mp_gather (ng, LBi, UBi, LBj, UBj, 1, 1, p2dvar, 1.0_r8,     &
     &                pmask(LBi,LBj), pmask(LBi,LBj), Npts, mask)
      CALL mp_bcastf (ng, mask, Npts)

      ic=0
      DO ij=1,Npts
        IF (mask(ij).eq.1.0_r8) THEN
          ic=ic+1
          SCALARS(ng)%IJwater(ic,p2dvar)=ij
        END IF
      END DO
      my_Nxyp=ic
!
!  U-points indices.
!
      CALL mp_gather (ng, LBi, UBi, LBj, UBj, 1, 1, u2dvar, 1.0_r8,     &
     &                umask(LBi,LBj), umask(LBi,LBj), Npts, mask)
      CALL mp_bcastf (ng, mask, Npts)

      ic=0
      DO ij=1,Npts
        IF (mask(ij).gt.0.0_r8) THEN
          ic=ic+1
          SCALARS(ng)%IJwater(ic,u2dvar)=ij
        END IF
      END DO
      my_Nxyu=ic
!
!  V-points indices.
!
      CALL mp_gather (ng, LBi, UBi, LBj, UBj, 1, 1, v2dvar, 1.0_r8,     &
     &                vmask(LBi,LBj), vmask(LBi,LBj), Npts, mask)
      CALL mp_bcastf (ng, mask, Npts)

      ic=0
      DO ij=1,Npts
        IF (mask(ij).gt.0.0_r8) THEN
          ic=ic+1
          SCALARS(ng)%IJwater(ic,v2dvar)=ij
        END IF
      END DO
      my_Nxyv=ic
# endif
      RETURN
      END SUBROUTINE wpoints_tile
#else
      SUBROUTINE wpoints
      RETURN
      END SUBROUTINE wpoints
#endif
