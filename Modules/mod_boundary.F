#include "cppdefs.h"
      MODULE mod_boundary
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  Open boundary conditions arrays:                                    !
!                                                                      !
!  t_east       Tracer (T units) eastern boundary conditions.          !
!  tG_east      Latest two-time snapshots of input tracer (T units)    !
!                 eastern boundary data.                               !
!  t_north      Tracer (T units) northern boundary conditions.         !
!  tG_north     Latest two-time snapshots of input tracer (T units)    !
!                 northern boundary data.                              !
!  t_south      Tracer (T units) southern boundary conditions.         !
!  tG_south     Latest two-time snapshots of input tracer (T units)    !
!                 southern boundary data.                              !
!  t_west       Tracer (T units) western boundary conditions.          !
!  tG_west      Latest two-time snapshots of input tracer (T units)    !
!                 western boundary data.                               !
!                                                                      !
!  ubar_east    2D u-momentum (m/s) eastern boundary conditions.       !
!  ubarG_east   Latest two-time snapshots of input 2D u-momentum       !
!                 (m/s) eastern boundary data.                         !
!  ubar_north   2D u-momentum (m/s) northern boundary conditions.      !
!  ubarG_north  Latest two-time snapshots of input 2D u-momentum       !
!                 (m/s) northern boundary data.                        !
!  ubar_south   2D u-momentum (m/s) southern boundary conditions.      !
!  ubarG_south  Latest two-time snapshots of input 2D u-momentum       !
!                 (m/s) southern boundary data.                        !
!  ubar_west    2D u-momentum (m/s) western boundary conditions.       !
!  ubarG_west   Latest two-time snapshots of input 2D u-momentum       !
!                 (m/s) western boundary data.                         !
!  vbar_east    2D v-momentum (m/s) eastern boundary conditions.       !
!  vbarG_east   Latest two-time snapshots of input 2D v-momentum       !
!                 (m/s) eastern boundary data.                         !
!  vbar_north   2D v-momentum (m/s) northern boundary conditions.      !
!  vbarG_north  Latest two-time snapshots of input 2D v-momentum       !
!                 (m/s) northern boundary data.                        !
!  vbar_south   2D v-momentum (m/s) southern boundary conditions.      !
!  vbarG_south  Latest two-time snapshots of input 2D v-momentum       !
!                 (m/s) southern boundary data.                        !
!  vbar_west    2D v-momentum (m/s) western boundary conditions.       !
!  vbarG_west   Latest two-time snapshots of input 2D v-momentum       !
!                 (m/s) western boundary data.                         !
!                                                                      !
!  u_east       3D u-momentum (m/s) eastern boundary conditions.       !
!  uG_east      Latest two-time snapshots of input 3D u-momentum       !
!                 (m/s) eastern boundary data.                         !
!  u_north      3D u-momentum (m/s) northern boundary conditions.      !
!  uG_north     Latest two-time snapshots of input 3D u-momentum       !
!                 (m/s) northern boundary data.                        !
!  u_south      3D u-momentum (m/s) southern boundary conditions.      !
!  uG_south     Latest two-time snapshots of input 3D u-momentum       !
!                 (m/s) southern boundary data.                        !
!  u_west       3D u-momentum (m/s) western boundary conditions.       !
!  uG_west      Latest two-time snapshots of input 3D u-momentum       !
!                 (m/s) western boundary data.                         !
!  v_east       3D v-momentum (m/s) eastern boundary conditions.       !
!  vG_east      Latest two-time snapshots of input 3D v-momentum       !
!                 (m/s) eastern boundary data.                         !
!  v_north      3D v-momentum (m/s) northern boundary conditions.      !
!  vG_north     Latest two-time snapshots of input 3D v-momentum       !
!                 (m/s) northern boundary data.                        !
!  v_south      3D v-momentum (m/s) southern boundary conditions.      !
!  vG_south     Latest two-time snapshots of input 3D v-momentum       !
!                 (m/s) southern boundary data.                        !
!  v_west       3D v-momentum (m/s) western boundary conditions.       !
!  vG_west      Latest two-time snapshots of input 3D v-momentum       !
!                 (m/s) western boundary data.                         !
!                                                                      !
!  zeta_east    Free-surface (m) eastern boundary conditions.          !
!  zetaG_east   Latest two-time snapshots of input free-surface (m)    !
!                 eastern boundary data.                               !
!  zeta_north   Free-surface (m) northern boundary conditions.         !
!  zetaG_north  Latest two-time snapshots of input free-surface (m)    !
!                 northern boundary data.                              !
!  zeta_south   Free-surface (m) southern boundary conditions.         !
!  zetaG_south  Latest two-time snapshots of input free-surface (m)    !
!                 southern boundary data.                              !
!  zeta_west    Free-surface (m) western boundary conditions.          !
!  zetaG_west   Latest two-time snapshots of input free-surface (m)    !
!                 western boundary data.                               !
!                                                                      !
!=======================================================================
!
        USE mod_kinds

        implicit none

        TYPE T_BOUNDARY

          real(r8), pointer :: zeta_west(:)
          real(r8), pointer :: zeta_east(:)
          real(r8), pointer :: zeta_south(:)
          real(r8), pointer :: zeta_north(:)

#ifndef ANA_FSOBC
          real(r8), pointer :: zetaG_west(:,:)
          real(r8), pointer :: zetaG_east(:,:)
          real(r8), pointer :: zetaG_south(:,:)
          real(r8), pointer :: zetaG_north(:,:)
#endif

          real(r8), pointer :: ubar_west(:)
          real(r8), pointer :: vbar_west(:)
          real(r8), pointer :: ubar_east(:)
          real(r8), pointer :: vbar_east(:)
          real(r8), pointer :: ubar_south(:)
          real(r8), pointer :: vbar_south(:)
          real(r8), pointer :: ubar_north(:)
          real(r8), pointer :: vbar_north(:)

#ifndef ANA_M2OBC
          real(r8), pointer :: ubarG_west(:,:)
          real(r8), pointer :: vbarG_west(:,:)
          real(r8), pointer :: ubarG_east(:,:)
          real(r8), pointer :: vbarG_east(:,:)
          real(r8), pointer :: ubarG_south(:,:)
          real(r8), pointer :: vbarG_south(:,:)
          real(r8), pointer :: ubarG_north(:,:)
          real(r8), pointer :: vbarG_north(:,:)
#endif

#ifdef SOLVE3D
          real(r8), pointer :: u_west(:,:)
          real(r8), pointer :: v_west(:,:)
          real(r8), pointer :: u_east(:,:)
          real(r8), pointer :: v_east(:,:)

# ifndef ANA_M3OBC
          real(r8), pointer :: uG_west(:,:,:)
          real(r8), pointer :: vG_west(:,:,:)
          real(r8), pointer :: uG_east(:,:,:)
          real(r8), pointer :: vG_east(:,:,:)
# endif

          real(r8), pointer :: u_south(:,:)
          real(r8), pointer :: v_south(:,:)
          real(r8), pointer :: u_north(:,:)
          real(r8), pointer :: v_north(:,:)

# ifndef ANA_M3OBC
          real(r8), pointer :: uG_south(:,:,:)
          real(r8), pointer :: vG_south(:,:,:)
          real(r8), pointer :: uG_north(:,:,:)
          real(r8), pointer :: vG_north(:,:,:)
# endif

          real(r8), pointer :: t_west(:,:,:)
          real(r8), pointer :: t_east(:,:,:)
          real(r8), pointer :: t_south(:,:,:)
          real(r8), pointer :: t_north(:,:,:)

# ifndef ANA_TOBC
          real(r8), pointer :: tG_west(:,:,:,:)
          real(r8), pointer :: tG_east(:,:,:,:)
          real(r8), pointer :: tG_south(:,:,:,:)
          real(r8), pointer :: tG_north(:,:,:,:)
# endif
#endif

        END TYPE T_BOUNDARY

        TYPE (T_BOUNDARY), allocatable ::BOUNDARY(:)

      CONTAINS

      SUBROUTINE allocate_boundary (ng)
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine initializes all variables in the module for all nested !
!  grids.  Currently, there is not parallel tiling in boundary arrays. !
!                                                                      !
!=======================================================================
!
      USE mod_param
!
!  Local variable declarations.
!
      integer, intent(in) :: ng
!
!-----------------------------------------------------------------------
!  Initialize module variables.
!-----------------------------------------------------------------------
!
      IF (ng.eq.1) allocate ( BOUNDARY(Ngrids) )
!
      allocate ( BOUNDARY(ng) % zeta_west(ETA_DIM) )
      allocate ( BOUNDARY(ng) % zeta_east(ETA_DIM) )
      allocate ( BOUNDARY(ng) % zeta_south(XI_DIM) )
      allocate ( BOUNDARY(ng) % zeta_north(XI_DIM) )

#ifndef ANA_FSOBC
      allocate ( BOUNDARY(ng) % zetaG_west(ETA_DIM,2) )
      allocate ( BOUNDARY(ng) % zetaG_east(ETA_DIM,2) )
      allocate ( BOUNDARY(ng) % zetaG_south(XI_DIM,2) )
      allocate ( BOUNDARY(ng) % zetaG_north(XI_DIM,2) )
#endif

      allocate ( BOUNDARY(ng) % ubar_west(ETA_DIM) )
      allocate ( BOUNDARY(ng) % ubar_east(ETA_DIM) )
      allocate ( BOUNDARY(ng) % ubar_south(XI_DIM) )
      allocate ( BOUNDARY(ng) % ubar_north(XI_DIM) )

      allocate ( BOUNDARY(ng) % vbar_west(ETA_DIM) )
      allocate ( BOUNDARY(ng) % vbar_east(ETA_DIM) )
      allocate ( BOUNDARY(ng) % vbar_south(XI_DIM) )
      allocate ( BOUNDARY(ng) % vbar_north(XI_DIM) )

#ifndef ANA_M2OBC
      allocate ( BOUNDARY(ng) % ubarG_west(ETA_DIM,2) )
      allocate ( BOUNDARY(ng) % ubarG_east(ETA_DIM,2) )
      allocate ( BOUNDARY(ng) % ubarG_south(XI_DIM,2) )
      allocate ( BOUNDARY(ng) % ubarG_north(XI_DIM,2) )

      allocate ( BOUNDARY(ng) % vbarG_west(ETA_DIM,2) )
      allocate ( BOUNDARY(ng) % vbarG_east(ETA_DIM,2) )
      allocate ( BOUNDARY(ng) % vbarG_south(XI_DIM,2) )
      allocate ( BOUNDARY(ng) % vbarG_north(XI_DIM,2) )
#endif

#ifdef SOLVE3D

      allocate ( BOUNDARY(ng) % u_west(ETA_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % u_east(ETA_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % u_south(XI_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % u_north(XI_DIM,N(ng)) )

      allocate ( BOUNDARY(ng) % v_west(ETA_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % v_east(ETA_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % v_south(XI_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % v_north(XI_DIM,N(ng)) )

# ifndef ANA_M3OBC
      allocate ( BOUNDARY(ng) % uG_west(ETA_DIM,N(ng),2) )
      allocate ( BOUNDARY(ng) % uG_east(ETA_DIM,N(ng),2) )
      allocate ( BOUNDARY(ng) % uG_south(XI_DIM,N(ng),2) )
      allocate ( BOUNDARY(ng) % uG_north(XI_DIM,N(ng),2) )

      allocate ( BOUNDARY(ng) % vG_west(ETA_DIM,N(ng),2) )
      allocate ( BOUNDARY(ng) % vG_east(ETA_DIM,N(ng),2) )
      allocate ( BOUNDARY(ng) % vG_south(XI_DIM,N(ng),2) )
      allocate ( BOUNDARY(ng) % vG_north(XI_DIM,N(ng),2) )
# endif

      allocate ( BOUNDARY(ng) % t_west(ETA_DIM,N(ng),NT(ng)) )
      allocate ( BOUNDARY(ng) % t_east(ETA_DIM,N(ng),NT(ng)) )
      allocate ( BOUNDARY(ng) % t_south(XI_DIM,N(ng),NT(ng)) )
      allocate ( BOUNDARY(ng) % t_north(XI_DIM,N(ng),NT(ng)) )

# ifndef ANA_TOBC
      allocate ( BOUNDARY(ng) % tG_west(ETA_DIM,N(ng),2,NT(ng)) )
      allocate ( BOUNDARY(ng) % tG_east(ETA_DIM,N(ng),2,NT(ng)) )
      allocate ( BOUNDARY(ng) % tG_south(XI_DIM,N(ng),2,NT(ng)) )
      allocate ( BOUNDARY(ng) % tG_north(XI_DIM,N(ng),2,NT(ng)) )
# endif
#endif

      RETURN
      END SUBROUTINE allocate_boundary

      SUBROUTINE initialize_boundary (ng, tile)
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine initialize all variables in the module using first     !
!  touch distribution policy. In shared-memory configuration, this     !
!  operation actually performs propagation of the  "shared arrays"     !
!  across the cluster, unless another policy is specified to           !
!  override the default.                                               !
!                                                                      !
!=======================================================================
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
      real(r8), parameter :: IniVal = 0.0_r8

#ifndef DISTRIBUTE
# include "tile.h"
#endif
!
!-----------------------------------------------------------------------
!  Initialize module variables.
!-----------------------------------------------------------------------
!
      IF (NORTH_WEST_TEST) THEN
        BOUNDARY(ng) % zeta_west(ETA_DIM) = IniVal
#ifndef ANA_FSOBC
        BOUNDARY(ng) % zetaG_west(ETA_DIM,1:2) = IniVal
#endif
        BOUNDARY(ng) % ubar_west(ETA_DIM) = IniVal
#ifndef ANA_M2OBC
        BOUNDARY(ng) % ubarG_west(ETA_DIM,1:2) = IniVal
#endif
        BOUNDARY(ng) % vbar_west(ETA_DIM) = IniVal
#ifndef ANA_M2OBC
        BOUNDARY(ng) % vbarG_west(ETA_DIM,1:2) = IniVal
#endif
#ifdef SOLVE3D
        BOUNDARY(ng) % u_west(ETA_DIM,1:N(ng)) = IniVal
# ifndef ANA_M3OBC
        BOUNDARY(ng) % uG_west(ETA_DIM,1:N(ng),1:2) = IniVal
# endif
        BOUNDARY(ng) % v_west(ETA_DIM,1:N(ng)) = IniVal
# ifndef ANA_M3OBC
        BOUNDARY(ng) % vG_west(ETA_DIM,1:N(ng),1:2) = IniVal
# endif
        BOUNDARY(ng) % t_west(ETA_DIM,1:N(ng),1:NT(ng)) = IniVal
# ifndef ANA_TOBC
        BOUNDARY(ng) % tG_west(ETA_DIM,1:N(ng),1:2,1:NT(ng)) =          &
     &                        IniVal
# endif
#endif
      END IF

      IF (SOUTH_EAST_TEST) THEN
        BOUNDARY(ng) % zeta_east(ETA_DIM) = IniVal
#ifndef ANA_FSOBC
        BOUNDARY(ng) % zetaG_east(ETA_DIM,1:2) = IniVal
#endif
        BOUNDARY(ng) % ubar_east(ETA_DIM) = IniVal
#ifndef ANA_M2OBC
        BOUNDARY(ng) % ubarG_east(ETA_DIM,1:2) = IniVal
#endif
        BOUNDARY(ng) % vbar_east(ETA_DIM) = IniVal
#ifndef ANA_M2OBC
        BOUNDARY(ng) % vbarG_east(ETA_DIM,1:2) = IniVal
#endif
#ifdef SOLVE3D
        BOUNDARY(ng) % u_east(ETA_DIM,1:N(ng)) = IniVal
# ifndef ANA_M3OBC
        BOUNDARY(ng) % uG_east(ETA_DIM,1:N(ng),1:2) = IniVal
# endif
        BOUNDARY(ng) % v_east(ETA_DIM,1:N(ng)) = IniVal
# ifndef ANA_M3OBC
        BOUNDARY(ng) % vG_east(ETA_DIM,1:N(ng),1:2) = IniVal
# endif
        BOUNDARY(ng) % t_east(ETA_DIM,1:N(ng),1:NT(ng)) = IniVal
# ifndef ANA_TOBC
        BOUNDARY(ng) % tG_east(ETA_DIM,1:N(ng),1:2,1:NT(ng)) =          &
     &                        IniVal
# endif
#endif
      END IF

      IF (SOUTH_WEST_TEST) THEN
        BOUNDARY(ng) % zeta_south(XI_DIM) = IniVal
#ifndef ANA_FSOBC
        BOUNDARY(ng) % zetaG_south(XI_DIM,1:2) = IniVal
#endif
        BOUNDARY(ng) % ubar_south(XI_DIM) = IniVal
#ifndef ANA_M2OBC
        BOUNDARY(ng) % ubarG_south(XI_DIM,1:2) = IniVal
#endif
        BOUNDARY(ng) % vbar_south(XI_DIM) = IniVal
#ifndef ANA_M2OBC
        BOUNDARY(ng) % vbarG_south(XI_DIM,1:2) = IniVal
#endif
#ifdef SOLVE3D
        BOUNDARY(ng) % u_south(XI_DIM,1:N(ng)) = IniVal
# ifndef ANA_M3OBC
        BOUNDARY(ng) % uG_south(XI_DIM,1:N(ng),1:2) = IniVal
# endif
        BOUNDARY(ng) % v_south(XI_DIM,1:N(ng)) = IniVal
# ifndef ANA_M3OBC
        BOUNDARY(ng) % vG_south(XI_DIM,1:N(ng),1:2) = IniVal
# endif
        BOUNDARY(ng) % t_south(XI_DIM,1:N(ng),1:NT(ng)) = IniVal
# ifndef ANA_TOBC
        BOUNDARY(ng) % tG_south(XI_DIM,1:N(ng),1:2,1:NT(ng)) =          &
     &                         IniVal
# endif
#endif
      END IF

      IF (NORTH_EAST_TEST) THEN
        BOUNDARY(ng) % zeta_north(XI_DIM) = IniVal
#ifndef ANA_FSOBC
        BOUNDARY(ng) % zetaG_north(XI_DIM,1:2) = IniVal
#endif
        BOUNDARY(ng) % ubar_north(XI_DIM) = IniVal
#ifndef ANA_M2OBC
        BOUNDARY(ng) % ubarG_north(XI_DIM,1:2) = IniVal
#endif
        BOUNDARY(ng) % vbar_north(XI_DIM) = IniVal
#ifndef ANA_M2OBC
        BOUNDARY(ng) % vbarG_north(XI_DIM,1:2) = IniVal
#endif
#ifdef SOLVE3D
        BOUNDARY(ng) % u_north(XI_DIM,1:N(ng)) = IniVal
# ifndef ANA_M3OBC
        BOUNDARY(ng) % uG_north(XI_DIM,1:N(ng),1:2) = IniVal
# endif
        BOUNDARY(ng) % v_north(XI_DIM,1:N(ng)) = IniVal
# ifndef ANA_M3OBC
        BOUNDARY(ng) % vG_north(XI_DIM,1:N(ng),1:2) = IniVal
# endif
        BOUNDARY(ng) % t_north(XI_DIM,1:N(ng),1:NT(ng)) = IniVal
# ifndef ANA_TOBC
        BOUNDARY(ng) % tG_north(XI_DIM,1:N(ng),1:2,1:NT(ng)) =          &
     &                         IniVal
# endif
#endif
      END IF
      RETURN
      END SUBROUTINE initialize_boundary

      END MODULE mod_boundary
