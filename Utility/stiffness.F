#include "cppdefs.h"
      MODULE stiffness_mod
#ifdef SOLVE3D
!
!========================================== Alexander F. Shchepetkin ===
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine surveys the 3D grid in order to determine maximum      !
!  grid stiffness ratio:                                               !
!                                                                      !
!             z(i,j,k)-z(i-1,j,k)+z(i,j,k-1)-z(i-1,j,k-1)              !
!      r_x = ---------------------------------------------             !
!             z(i,j,k)+z(i-1,j,k)-z(i,j,k-1)-z(i-1,j,k-1)              !
!                                                                      !
!  This is done for diagnostic purposes and it does not affect the     !
!  computations.                                                       !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: stiffness

      CONTAINS
!
!***********************************************************************      
      SUBROUTINE stiffness (ng, tile)
!***********************************************************************      
!
      USE mod_param
      USE mod_grid
!
      integer, intent(in) :: ng, tile

# include "tile.h"
!
      CALL stiffness_tile (ng, Istr, Iend, Jstr, Jend,                  &
     &                     LBi, UBi, LBj, UBj,                          &
# ifdef MASKING
     &                     GRID(ng) % umask,                            &
     &                     GRID(ng) % vmask,                            &
# endif
     &                     GRID(ng) % z_w)
      RETURN
      END SUBROUTINE stiffness
!
!***********************************************************************
      SUBROUTINE stiffness_tile (ng, Istr, Iend, Jstr, Jend,            &
     &                           LBi, UBi, LBj, UBj,                    &
# ifdef MASKING
     &                           umask, vmask,                          &
# endif
     &                           z_w)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj

# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
# endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: NSUB, i, j, k
      real(r8) :: my_rx0, my_rx1

# ifdef DISTRIBUTE
      real(r8), dimension(2) :: buffer
      character (len=3), dimension(2) :: op_handle
# endif

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute grid stiffness.
!-----------------------------------------------------------------------
!
      my_rx0=0.0_r8
      my_rx1=0.0_r8
!
      DO j=Jstr,Jend
        DO i=IstrU,Iend
# ifdef MASKING
          IF (umask(i,j).gt.0.0_r8) THEN
# endif
            my_rx0=MAX(my_rx0,ABS((z_w(i,j,0)-z_w(i-1,j,0))/            &
     &                            (z_w(i,j,0)+z_w(i-1,j,0))))
            DO k=1,N(ng)
              my_rx1=MAX(my_rx1,ABS((z_w(i,j,k  )-z_w(i-1,j,k  )+       &
     &                               z_w(i,j,k-1)-z_w(i-1,j,k-1))/      &
     &                              (z_w(i,j,k  )+z_w(i-1,j,k  )-       &
     &                               z_w(i,j,k-1)-z_w(i-1,j,k-1))))
            END DO
# ifdef MASKING
          END IF
# endif
        END DO
      END DO
      DO j=JstrV,Jend
        DO i=Istr,Iend
# ifdef MASKING
          IF (vmask(i,j).gt.0.0_r8) THEN
# endif
            my_rx0=MAX(my_rx0,ABS((z_w(i,j,0)-z_w(i,j-1,0))/            &
     &                            (z_w(i,j,0)+z_w(i,j-1,0))))
            DO k=1,N(ng)
              my_rx1=MAX(my_rx1,ABS((z_w(i,j,k  )-z_w(i,j-1,k  )+       &
     &                               z_w(i,j,k-1)-z_w(i,j-1,k-1))/      &
     &                              (z_w(i,j,k  )+z_w(i,j-1,k  )-       &
     &                               z_w(i,j,k-1)-z_w(i,j-1,k-1))))
            END DO
# ifdef MASKING
          END IF
# endif
        END DO
      END DO
!
!  Compute global maximum stiffness.
!
      IF (SOUTH_WEST_CORNER.and.                                        &
     &    NORTH_EAST_CORNER) THEN
        NSUB=1                           ! non-tiled application
      ELSE
        NSUB=NtileX(ng)*NtileE(ng)       ! tiled application
      END IF
!$OMP CRITICAL (R_FACTOR)
      IF (tile_count.eq.0) THEN
        rx0=my_rx0
        rx1=my_rx1
      ELSE
        rx0=MAX(rx0,my_rx0)
        rx1=MAX(rx1,my_rx1)
      END IF
      tile_count=tile_count+1
      IF (tile_count.eq.NSUB) THEN
        tile_count=0
# ifdef DISTRIBUTE
        buffer(1)=rx0
        buffer(2)=rx1
        op_handle(1)='MAX'
        op_handle(2)='MAX'
        CALL mp_reduce (ng, 2, buffer, op_handle)
        rx0=buffer(1)
        rx1=buffer(2)
# endif
        IF (Master) WRITE (stdout,10) rx0, rx1
  10    FORMAT (/,' Maximum grid stiffness ratios:  rx0 = ',1pe14.6,     &
     &            ' (Beckmann and Haidvogel)',/,t34,'rx1 = ',1pe14.6,    &
     &            ' (Haney)',/)
      END IF
!$OMP END CRITICAL (R_FACTOR)
      RETURN
      END SUBROUTINE stiffness_tile
#endif
      END MODULE stiffness_mod

