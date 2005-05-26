#include "cppdefs.h"
      MODULE mod_boundary
!
!=======================================================================
!  Copyright (c) 2005 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  Open boundary conditions arrays:                                    !
!                                                                      !
!  t_east         Tracer (T units) eastern boundary conditions.        !
!  tG_east        Latest two-time snapshots of input tracer (Tunits)   !
!                   eastern boundary data.                             !
!  t_north        Tracer (T units) northern boundary conditions.       !
!  tG_north       Latest two-time snapshots of input tracer (Tunits)   !
!                   northern boundary data.                            !
!  t_south        Tracer (T units) southern boundary conditions.       !
!  tG_south       Latest two-time snapshots of input tracer (Tunits)   !
!                   southern boundary data.                            !
!  t_west         Tracer (T units) western boundary conditions.        !
!  tG_west        Latest two-time snapshots of input tracer (Tunits)   !
!                   western boundary data.                             !
!                                                                      !
!  ubar_east      2D u-momentum (m/s) eastern boundary conditions.     !
!  ubarG_east     Latest two-time snapshots of input 2D u-momentum     !
!                   (m/s) eastern boundary data.                       !
!  ubar_north     2D u-momentum (m/s) northern boundary conditions.    !
!  ubarG_north    Latest two-time snapshots of input 2D u-momentum     !
!                   (m/s) northern boundary data.                      !
!  ubar_south     2D u-momentum (m/s) southern boundary conditions.    !
!  ubarG_south    Latest two-time snapshots of input 2D u-momentum     !
!                   (m/s) southern boundary data.                      !
!  ubar_west      2D u-momentum (m/s) western boundary conditions.     !
!  ubarG_west     Latest two-time snapshots of input 2D u-momentum     !
!                   (m/s) western boundary data.                       !
!  vbar_east      2D v-momentum (m/s) eastern boundary conditions.     !
!  vbarG_east     Latest two-time snapshots of input 2D v-momentum     !
!                   (m/s) eastern boundary data.                       !
!  vbar_north     2D v-momentum (m/s) northern boundary conditions.    !
!  vbarG_north    Latest two-time snapshots of input 2D v-momentum     !
!                   (m/s) northern boundary data.                      !
!  vbar_south     2D v-momentum (m/s) southern boundary conditions.    !
!  vbarG_south    Latest two-time snapshots of input 2D v-momentum     !
!                   (m/s) southern boundary data.                      !
!  vbar_west      2D v-momentum (m/s) western boundary conditions.     !
!  vbarG_west     Latest two-time snapshots of input 2D v-momentum     !
!                   (m/s) western boundary data.                       !
!                                                                      !
!  u_east         3D u-momentum (m/s) eastern boundary conditions.     !
!  uG_east        Latest two-time snapshots of input 3D u-momentum     !
!                   (m/s) eastern boundary data.                       !
!  u_north        3D u-momentum (m/s) northern boundary conditions.    !
!  uG_north       Latest two-time snapshots of input 3D u-momentum     !
!                   (m/s) northern boundary data.                      !
!  u_south        3D u-momentum (m/s) southern boundary conditions.    !
!  uG_south       Latest two-time snapshots of input 3D u-momentum     !
!                   (m/s) southern boundary data.                      !
!  u_west         3D u-momentum (m/s) western boundary conditions.     !
!  uG_west        Latest two-time snapshots of input 3D u-momentum     !
!                   (m/s) western boundary data.                       !
!  v_east         3D v-momentum (m/s) eastern boundary conditions.     !
!  vG_east        Latest two-time snapshots of input 3D v-momentum     !
!                   (m/s) eastern boundary data.                       !
!  v_north        3D v-momentum (m/s) northern boundary conditions.    !
!  vG_north       Latest two-time snapshots of input 3D v-momentum     !
!                   (m/s) northern boundary data.                      !
!  v_south        3D v-momentum (m/s) southern boundary conditions.    !
!  vG_south       Latest two-time snapshots of input 3D v-momentum     !
!                   (m/s) southern boundary data.                      !
!  v_west         3D v-momentum (m/s) western boundary conditions.     !
!  vG_west        Latest two-time snapshots of input 3D v-momentum     !
!                   (m/s) western boundary data.                       !
!                                                                      !
!  zeta_east      Free-surface (m) eastern boundary conditions.        !
!  zetaG_east     Latest two-time snapshots of input free-surface (m)  !
!                   eastern boundary data.                             !
!  zeta_north     Free-surface (m) northern boundary conditions.       !
!  zetaG_north    Latest two-time snapshots of input free-surface (m)  !
!                   northern boundary data.                            !
!  zeta_south     Free-surface (m) southern boundary conditions.       !
!  zetaG_south    Latest two-time snapshots of input free-surface (m)  !
!                   southern boundary data.                            !
!  zeta_west      Free-surface (m) western boundary conditions.        !
!  zetaG_west     Latest two-time snapshots of input free-surface (m)  !
!                   western boundary data.                             !
!                                                                      !
!  Boundary arrays for radiation celerities needed for tangent linear  !
!  and adjoint models applications:                                    !
!                                                                      !
!  t_east_C2      Tracers eastern boundary, quadratic radition term    !
!                   term dTdx*dTdx+dTde*dTde (Tunits squared).         !
!  t_east_Cx      Tracers eastern boundary, celerity in the            !
!                   XI-direction (Tunits).                             !
!  t_east_Ce      Tracers eastern boundary, celerity in the            !
!                   ETA-direction (Tunits).                            !
!                                                                      !
!  t_north_C2     Tracers northern boundary, quadratic radition term   !
!                   dTdx*dTdx+dTde*dTde (Tunits squared).              !
!  t_north_Cx     Tracers northern boundary, celerity in the           !
!                   XI-direction (Tunits).                             !
!  t_north_Ce     Tracers northern boundary, celerity in the           !
!                   ETA-direction (Tunits).                            !
!                                                                      !
!  t_south_C2     Tracers southern boundary, quadratic radition term   !
!                   term dTdx*dTdx+dTde*dTde (Tunits squared).         !
!  t_south_Cx     Tracers southern boundary, celerity in the           !
!                   XI-direction (Tunits).                             !
!  t_south_Ce     Tracers northern boundary, celerity in the           !
!                   ETA-direction (Tunits).                            !
!                                                                      !
!  t_west_C2      Tracers western boundary, quadratic radition term    !
!                   dTdx*dTdx+dTde*dTde (Tunits squared).              !
!  t_west_Cx      Tracers western boundary, celerity in the            !
!                   XI-direction (Tunits).                             !
!  t_west_Ce      Tracers western boundary, celerity in the            !
!                   ETA-direction (Tunits).                            !
!                                                                      !
!  ubar_east_C2   2D U-momentum eastern boundary, quadratic radition   !
!                   term dUdx*dUdx+dUde*dUde (m2/s2).                  !
!  ubar_east_Cx   2D U-momentum eastern boundary, celerity in the      !
!                   XI-direction (m/s).                                !
!  ubar_east_Ce   2D U-momentum eastern boundary, celerity in the      !
!                   ETA-direction (m/s).                               !
!                                                                      !
!  ubar_north_C2  2D U-momentum northern boundary, quadratic radition  !
!                   term dUdx*dUdx+dUde*dUde (m2/s2).                  !
!  ubar_north_Cx  2D U-momentum northern boundary, celerity in the     !
!                   XI-direction (m/s).                                !
!  ubar_north_Ce  2D U-momentum northern boundary, celerity in the     !
!                   ETA-direction (m/s).                               !
!                                                                      !
!  ubar_south_C2  2D U-momentum southern boundary, quadratic radition  !
!                   term dUdx*dUdx+dUde*dUde (m2/s2).                  !
!  ubar_south_Cx  2D U-momentum southern boundary, celerity in the     !
!                   XI-direction (m/s).                                !
!  ubar_south_Ce  2D U-momentum northern boundary, celerity in the     !
!                   ETA-direction (m/s).                               !
!                                                                      !
!  ubar_west_C2   2D U-momentum western boundary, quadratic radition   !
!                   term dUdx*dUdx+dUde*dUde (m2/s2).                  !
!  ubar_west_Cx   2D U-momentum western boundary, celerity in the      !
!                   XI-direction (m/s).                                !
!  ubar_west_Ce   2D U-momentum western boundary, celerity in the      !
!                   ETA-direction (m/s).                               !
!                                                                      !
!  vbar_east_C2   2D V-momentum eastern boundary, quadratic radition   !
!                   term dUdx*dUdx+dUde*dUde (m2/s2).                  !
!  vbar_east_Cx   2D V-momentum eastern boundary, celerity in the      !
!                   XI-direction (m/s).                                !
!  vbar_east_Ce   2D V-momentum eastern boundary, celerity in the      !
!                   ETA-direction (m/s).                               !
!                                                                      !
!  vbar_north_C2  2D V-momentum northern boundary, quadratic radition  !
!                   term dUdx*dUdx+dUde*dUde (m2/s2).                  !
!  vbar_north_Cx  2D V-momentum northern boundary, celerity in the     !
!                   XI-direction (m/s).                                !
!  vbar_north_Ce  2D V-momentum northern boundary, celerity in the     !
!                   ETA-direction (m/s).                               !
!                                                                      !
!  vbar_south_C2  2D V-momentum southern boundary, quadratic radition  !
!                   term dUdx*dUdx+dUde*dUde (m2/s2).                  !
!  vbar_south_Cx  2D V-momentum southern boundary, celerity in the     !
!                   XI-direction (m/s).                                !
!  vbar_south_Ce  2D V-momentum northern boundary, celerity in the     !
!                   ETA-direction (m/s).                               !
!                                                                      !
!  vbar_west_C2   2D V-momentum western boundary, quadratic radition   !
!                   term dUdx*dUdx+dUde*dUde (m2/s2).                  !
!  vbar_west_Cx   2D V-momentum western boundary, celerity in the      !
!                   XI-direction (m/s).                                !
!  vbar_west_Ce   2D V-momentum western boundary, celerity in the      !
!                   ETA-direction (m/s).                               !
!
!  u_east_C2      3D U-momentum eastern boundary, quadratic radition   !
!                   term dUdx*dUdx+dUde*dUde (m2/s2).                  !
!  u_east_Cx      3D U-momentum eastern boundary, celerity in the      !
!                   XI-direction (m/s).                                !
!  u_east_Ce      3D U-momentum eastern boundary, celerity in the      !
!                   ETA-direction (m/s).                               !
!                                                                      !
!  u_north_C2     3D U-momentum northern boundary, quadratic radition  !
!                   term dUdx*dUdx+dUde*dUde (m2/s2).                  !
!  u_north_Cx     3D U-momentum northern boundary, celerity in the     !
!                   XI-direction (m/s).                                !
!  u_north_Ce     3D U-momentum northern boundary, celerity in the     !
!                   ETA-direction (m/s).                               !
!                                                                      !
!  u_south_C2     3D U-momentum southern boundary, quadratic radition  !
!                   term dUdx*dUdx+dUde*dUde (m2/s2).                  !
!  u_south_Cx     3D U-momentum southern boundary, celerity in the     !
!                   XI-direction (m/s).                                !
!  u_south_Ce     3D U-momentum northern boundary, celerity in the     !
!                   ETA-direction (m/s).                               !
!                                                                      !
!  u_west_C2      3D U-momentum western boundary, quadratic radition   !
!                   term dUdx*dUdx+dUde*dUde (m2/s2).                  !
!  u_west_Cx      3D U-momentum western boundary, celerity in the      !
!                   XI-direction (m/s).                                !
!  u_west_Ce      3D U-momentum western boundary, celerity in the      !
!                   ETA-direction (m/s).                               !
!                                                                      !
!  v_east_C2      3D V-momentum eastern boundary, quadratic radition   !
!                   term dUdx*dUdx+dUde*dUde (m2/s2).                  !
!  v_east_Cx      3D V-momentum eastern boundary, celerity in the      !
!                   XI-direction (m/s).                                !
!  v_east_Ce      3D V-momentum eastern boundary, celerity in the      !
!                   ETA-direction (m/s).                               !
!                                                                      !
!  v_north_C2     3D V-momentum northern boundary, quadratic radition  !
!                   term dUdx*dUdx+dUde*dUde (m2/s2).                  !
!  v_north_Cx     3D V-momentum northern boundary, celerity in the     !
!                   XI-direction (m/s).                                !
!  v_north_Ce     3D V-momentum northern boundary, celerity in the     !
!                   ETA-direction (m/s).                               !
!                                                                      !
!  v_south_C2     3D V-momentum southern boundary, quadratic radition  !
!                   term dUdx*dUdx+dUde*dUde (m2/s2).                  !
!  v_south_Cx     3D V-momentum southern boundary, celerity in the     !
!                   XI-direction (m/s).                                !
!  v_south_Ce     3D V-momentum northern boundary, celerity in the     !
!                   ETA-direction (m/s).                               !
!                                                                      !
!  v_west_C2      3D V-momentum western boundary, quadratic radition   !
!                   term dUdx*dUdx+dUde*dUde (m2/s2).                  !
!  v_west_Cx      3D V-momentum western boundary, celerity in the      !
!                   XI-direction (m/s).                                !
!  v_west_Ce      3D V-momentum western boundary, celerity in the      !
!                   ETA-direction (m/s).                               !
!                                                                      !
!  zeta_east_C2   Free-surface eastern boundary, quadratic radition    !
!                   term dZdx*dZdx+dZde*dZde (m2).                     !
!  zeta_east_Cx   Free-surface eastern boundary, celerity in the       !
!                   XI-direction (m).                                  !
!  zeta_east_Ce   Free-surface eastern boundary, celerity in the       !
!                   ETA-direction (m).                                 !
!                                                                      !
!  zeta_north_C2  Free-surface northern boundary, quadratic radition   !
!                   term dZdx*dZdx+dZde*dZde (m2).                     !
!  zeta_north_Cx  Free-surface northern boundary, celerity in the      !
!                   XI-direction (m).                                  !
!  zeta_north_Ce  Free-surface northern boundary, celerity in the      !
!                   ETA-direction (m).                                 !
!                                                                      !
!  zeta_south_C2  Free-surface southern boundary, quadratic radition   !
!                   term dZdx*dZdx+dZde*dZde (m2).                     !
!  zeta_south_Cx  Free-surface southern boundary, celerity in the      !
!                   XI-direction (m).                                  !
!  zeta_south_Ce  Free-surface northern boundary, celerity in the      !
!                   ETA-direction (m).                                 !
!                                                                      !
!  zeta_west_C2   Free-surface western boundary, quadratic radition    !
!                   term dZdx*dZdx+dZde*dZde (m2).                     !
!  zeta_west_Cx   Free-surface western boundary, celerity in the       !
!                   XI-direction (m).                                  !
!  zeta_west_Ce   Free-surface western boundary, celerity in the       !
!                   ETA-direction (m).                                 !
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

#if defined CELERITY_READ || defined CELERITY_WRITE
          real(r8), pointer :: zeta_west_C2(:)
          real(r8), pointer :: zeta_west_Ce(:)
          real(r8), pointer :: zeta_west_Cx(:)

          real(r8), pointer :: zeta_east_C2(:)
          real(r8), pointer :: zeta_east_Ce(:)
          real(r8), pointer :: zeta_east_Cx(:)

          real(r8), pointer :: zeta_south_C2(:)
          real(r8), pointer :: zeta_south_Ce(:)
          real(r8), pointer :: zeta_south_Cx(:)

          real(r8), pointer :: zeta_north_C2(:)
          real(r8), pointer :: zeta_north_Ce(:)
          real(r8), pointer :: zeta_north_Cx(:)
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

#if defined CELERITY_READ || defined CELERITY_WRITE
          real(r8), pointer :: ubar_west_C2(:)
          real(r8), pointer :: ubar_west_Ce(:)
          real(r8), pointer :: ubar_west_Cx(:)

          real(r8), pointer :: ubar_east_C2(:)
          real(r8), pointer :: ubar_east_Ce(:)
          real(r8), pointer :: ubar_east_Cx(:)

          real(r8), pointer :: ubar_south_C2(:)
          real(r8), pointer :: ubar_south_Ce(:)
          real(r8), pointer :: ubar_south_Cx(:)

          real(r8), pointer :: ubar_north_C2(:)
          real(r8), pointer :: ubar_north_Ce(:)
          real(r8), pointer :: ubar_north_Cx(:)

          real(r8), pointer :: vbar_west_C2(:)
          real(r8), pointer :: vbar_west_Ce(:)
          real(r8), pointer :: vbar_west_Cx(:)

          real(r8), pointer :: vbar_east_C2(:)
          real(r8), pointer :: vbar_east_Ce(:)
          real(r8), pointer :: vbar_east_Cx(:)

          real(r8), pointer :: vbar_south_C2(:)
          real(r8), pointer :: vbar_south_Ce(:)
          real(r8), pointer :: vbar_south_Cx(:)

          real(r8), pointer :: vbar_north_C2(:)
          real(r8), pointer :: vbar_north_Ce(:)
          real(r8), pointer :: vbar_north_Cx(:)
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

# if defined CELERITY_READ || defined CELERITY_WRITE
          real(r8), pointer :: u_west_C2(:,:)
          real(r8), pointer :: u_west_Ce(:,:)
          real(r8), pointer :: u_west_Cx(:,:)

          real(r8), pointer :: u_east_C2(:,:)
          real(r8), pointer :: u_east_Ce(:,:)
          real(r8), pointer :: u_east_Cx(:,:)

          real(r8), pointer :: u_south_C2(:,:)
          real(r8), pointer :: u_south_Ce(:,:)
          real(r8), pointer :: u_south_Cx(:,:)

          real(r8), pointer :: u_north_C2(:,:)
          real(r8), pointer :: u_north_Ce(:,:)
          real(r8), pointer :: u_north_Cx(:,:)

          real(r8), pointer :: v_west_C2(:,:)
          real(r8), pointer :: v_west_Ce(:,:)
          real(r8), pointer :: v_west_Cx(:,:)

          real(r8), pointer :: v_east_C2(:,:)
          real(r8), pointer :: v_east_Ce(:,:)
          real(r8), pointer :: v_east_Cx(:,:)

          real(r8), pointer :: v_south_C2(:,:)
          real(r8), pointer :: v_south_Ce(:,:)
          real(r8), pointer :: v_south_Cx(:,:)

          real(r8), pointer :: v_north_C2(:,:)
          real(r8), pointer :: v_north_Ce(:,:)
          real(r8), pointer :: v_north_Cx(:,:)
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

# if defined CELERITY_READ || defined CELERITY_WRITE
          real(r8), pointer :: t_west_C2(:,:,:)
          real(r8), pointer :: t_west_Ce(:,:,:)
          real(r8), pointer :: t_west_Cx(:,:,:)

          real(r8), pointer :: t_east_C2(:,:,:)
          real(r8), pointer :: t_east_Ce(:,:,:)
          real(r8), pointer :: t_east_Cx(:,:,:)

          real(r8), pointer :: t_south_C2(:,:,:)
          real(r8), pointer :: t_south_Ce(:,:,:)
          real(r8), pointer :: t_south_Cx(:,:,:)

          real(r8), pointer :: t_north_C2(:,:,:)
          real(r8), pointer :: t_north_Ce(:,:,:)
          real(r8), pointer :: t_north_Cx(:,:,:)
# endif

#endif

        END TYPE T_BOUNDARY

        TYPE (T_BOUNDARY), allocatable ::BOUNDARY(:)

      CONTAINS

      SUBROUTINE allocate_boundary (ng)
!
!=======================================================================
!  Copyright (c) 2005 ROMS/TOMS Group                                  !
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

#if defined CELERITY_READ || defined CELERITY_WRITE
      allocate ( BOUNDARY(ng) % zeta_west_C2(ETA_DIM) )
      allocate ( BOUNDARY(ng) % zeta_west_Ce(ETA_DIM) )
      allocate ( BOUNDARY(ng) % zeta_west_Cx(ETA_DIM) )

      allocate ( BOUNDARY(ng) % zeta_east_C2(ETA_DIM) )
      allocate ( BOUNDARY(ng) % zeta_east_Ce(ETA_DIM) )
      allocate ( BOUNDARY(ng) % zeta_east_Cx(ETA_DIM) )

      allocate ( BOUNDARY(ng) % zeta_south_C2(XI_DIM) )
      allocate ( BOUNDARY(ng) % zeta_south_Ce(XI_DIM) )
      allocate ( BOUNDARY(ng) % zeta_south_Cx(XI_DIM) )

      allocate ( BOUNDARY(ng) % zeta_north_C2(XI_DIM) )
      allocate ( BOUNDARY(ng) % zeta_north_Ce(XI_DIM) )
      allocate ( BOUNDARY(ng) % zeta_north_Cx(XI_DIM) )
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

#if defined CELERITY_READ || defined CELERITY_WRITE
      allocate ( BOUNDARY(ng) % ubar_west_C2(ETA_DIM) )
      allocate ( BOUNDARY(ng) % ubar_west_Ce(ETA_DIM) )
      allocate ( BOUNDARY(ng) % ubar_west_Cx(ETA_DIM) )

      allocate ( BOUNDARY(ng) % ubar_east_C2(ETA_DIM) )
      allocate ( BOUNDARY(ng) % ubar_east_Ce(ETA_DIM) )
      allocate ( BOUNDARY(ng) % ubar_east_Cx(ETA_DIM) )

      allocate ( BOUNDARY(ng) % ubar_south_C2(XI_DIM) )
      allocate ( BOUNDARY(ng) % ubar_south_Ce(XI_DIM) )
      allocate ( BOUNDARY(ng) % ubar_south_Cx(XI_DIM) )

      allocate ( BOUNDARY(ng) % ubar_north_C2(XI_DIM) )
      allocate ( BOUNDARY(ng) % ubar_north_Ce(XI_DIM) )
      allocate ( BOUNDARY(ng) % ubar_north_Cx(XI_DIM) )

      allocate ( BOUNDARY(ng) % vbar_west_C2(ETA_DIM) )
      allocate ( BOUNDARY(ng) % vbar_west_Ce(ETA_DIM) )
      allocate ( BOUNDARY(ng) % vbar_west_Cx(ETA_DIM) )

      allocate ( BOUNDARY(ng) % vbar_east_C2(ETA_DIM) )
      allocate ( BOUNDARY(ng) % vbar_east_Ce(ETA_DIM) )
      allocate ( BOUNDARY(ng) % vbar_east_Cx(ETA_DIM) )

      allocate ( BOUNDARY(ng) % vbar_south_C2(XI_DIM) )
      allocate ( BOUNDARY(ng) % vbar_south_Ce(XI_DIM) )
      allocate ( BOUNDARY(ng) % vbar_south_Cx(XI_DIM) )

      allocate ( BOUNDARY(ng) % vbar_north_C2(XI_DIM) )
      allocate ( BOUNDARY(ng) % vbar_north_Ce(XI_DIM) )
      allocate ( BOUNDARY(ng) % vbar_north_Cx(XI_DIM) )
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

# if defined CELERITY_READ || defined CELERITY_WRITE
      allocate ( BOUNDARY(ng) % u_west_C2(ETA_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % u_west_Ce(ETA_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % u_west_Cx(ETA_DIM,N(ng)) )

      allocate ( BOUNDARY(ng) % u_east_C2(ETA_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % u_east_Ce(ETA_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % u_east_Cx(ETA_DIM,N(ng)) )

      allocate ( BOUNDARY(ng) % u_south_C2(XI_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % u_south_Ce(XI_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % u_south_Cx(XI_DIM,N(ng)) )

      allocate ( BOUNDARY(ng) % u_north_C2(XI_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % u_north_Ce(XI_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % u_north_Cx(XI_DIM,N(ng)) )

      allocate ( BOUNDARY(ng) % v_west_C2(ETA_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % v_west_Ce(ETA_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % v_west_Cx(ETA_DIM,N(ng)) )

      allocate ( BOUNDARY(ng) % v_east_C2(ETA_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % v_east_Ce(ETA_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % v_east_Cx(ETA_DIM,N(ng)) )

      allocate ( BOUNDARY(ng) % v_south_C2(XI_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % v_south_Ce(XI_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % v_south_Cx(XI_DIM,N(ng)) )

      allocate ( BOUNDARY(ng) % v_north_C2(XI_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % v_north_Ce(XI_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % v_north_Cx(XI_DIM,N(ng)) )
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

# if defined CELERITY_READ || defined CELERITY_WRITE
      allocate ( BOUNDARY(ng) % t_west_C2(ETA_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % t_west_Ce(ETA_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % t_west_Cx(ETA_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % t_east_C2(ETA_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % t_east_Ce(ETA_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % t_east_Cx(ETA_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % t_south_C2(XI_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % t_south_Ce(XI_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % t_south_Cx(XI_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % t_north_C2(XI_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % t_north_Ce(XI_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % t_north_Cx(XI_DIM,N(ng)) )
# endif
#endif

      RETURN
      END SUBROUTINE allocate_boundary

      SUBROUTINE initialize_boundary (ng, tile)
!
!=======================================================================
!  Copyright (c) 2005 ROMS/TOMS Group                                  !
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
!  Western boundary variables.
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
#if defined CELERITY_READ || defined CELERITY_WRITE
        BOUNDARY(ng) % zeta_west_C2(ETA_DIM) = IniVal
        BOUNDARY(ng) % zeta_west_Ce(ETA_DIM) = IniVal
        BOUNDARY(ng) % zeta_west_Cx(ETA_DIM) = IniVal
        BOUNDARY(ng) % ubar_west_C2(ETA_DIM) = IniVal
        BOUNDARY(ng) % ubar_west_Ce(ETA_DIM) = IniVal
        BOUNDARY(ng) % ubar_west_Cx(ETA_DIM) = IniVal
        BOUNDARY(ng) % vbar_west_C2(ETA_DIM) = IniVal
        BOUNDARY(ng) % vbar_west_Ce(ETA_DIM) = IniVal
        BOUNDARY(ng) % vbar_west_Cx(ETA_DIM) = IniVal
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
# if defined CELERITY_READ || defined CELERITY_WRITE
        BOUNDARY(ng) % u_west_C2(ETA_DIM,1:N(ng)) = IniVal
        BOUNDARY(ng) % u_west_Ce(ETA_DIM,1:N(ng)) = IniVal
        BOUNDARY(ng) % u_west_Cx(ETA_DIM,1:N(ng)) = IniVal
        BOUNDARY(ng) % v_west_C2(ETA_DIM,1:N(ng)) = IniVal
        BOUNDARY(ng) % v_west_Ce(ETA_DIM,1:N(ng)) = IniVal
        BOUNDARY(ng) % v_west_Cx(ETA_DIM,1:N(ng)) = IniVal
        BOUNDARY(ng) % t_west_C2(ETA_DIM,1:N(ng),1:NT(ng)) = IniVal
        BOUNDARY(ng) % t_west_Ce(ETA_DIM,1:N(ng),1:NT(ng)) = IniVal
        BOUNDARY(ng) % t_west_Cx(ETA_DIM,1:N(ng),1:NT(ng)) = IniVal
# endif
#endif
      END IF
!
!  Eastern boundary variables.
!
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
#if defined CELERITY_READ || defined CELERITY_WRITE
        BOUNDARY(ng) % zeta_east_C2(ETA_DIM) = IniVal
        BOUNDARY(ng) % zeta_east_Ce(ETA_DIM) = IniVal
        BOUNDARY(ng) % zeta_east_Cx(ETA_DIM) = IniVal
        BOUNDARY(ng) % ubar_east_C2(ETA_DIM) = IniVal
        BOUNDARY(ng) % ubar_east_Ce(ETA_DIM) = IniVal
        BOUNDARY(ng) % ubar_east_Cx(ETA_DIM) = IniVal
        BOUNDARY(ng) % vbar_east_C2(ETA_DIM) = IniVal
        BOUNDARY(ng) % vbar_east_Ce(ETA_DIM) = IniVal
        BOUNDARY(ng) % vbar_east_Cx(ETA_DIM) = IniVal
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
# if defined CELERITY_READ || defined CELERITY_WRITE
        BOUNDARY(ng) % u_east_C2(ETA_DIM,1:N(ng)) = IniVal
        BOUNDARY(ng) % u_east_Ce(ETA_DIM,1:N(ng)) = IniVal
        BOUNDARY(ng) % u_east_Cx(ETA_DIM,1:N(ng)) = IniVal
        BOUNDARY(ng) % v_east_C2(ETA_DIM,1:N(ng)) = IniVal
        BOUNDARY(ng) % v_east_Ce(ETA_DIM,1:N(ng)) = IniVal
        BOUNDARY(ng) % v_east_Cx(ETA_DIM,1:N(ng)) = IniVal
        BOUNDARY(ng) % t_east_C2(ETA_DIM,1:N(ng),1:NT(ng)) = IniVal
        BOUNDARY(ng) % t_east_Ce(ETA_DIM,1:N(ng),1:NT(ng)) = IniVal
        BOUNDARY(ng) % t_east_Cx(ETA_DIM,1:N(ng),1:NT(ng)) = IniVal
# endif
#endif
      END IF
!
!  Southern boundary variables.
!
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
#if defined CELERITY_READ || defined CELERITY_WRITE
        BOUNDARY(ng) % zeta_south_C2(XI_DIM) = IniVal
        BOUNDARY(ng) % zeta_south_Ce(XI_DIM) = IniVal
        BOUNDARY(ng) % zeta_south_Cx(XI_DIM) = IniVal
        BOUNDARY(ng) % ubar_south_C2(XI_DIM) = IniVal
        BOUNDARY(ng) % ubar_south_Ce(XI_DIM) = IniVal
        BOUNDARY(ng) % ubar_south_Cx(XI_DIM) = IniVal
        BOUNDARY(ng) % vbar_south_C2(XI_DIM) = IniVal
        BOUNDARY(ng) % vbar_south_Ce(XI_DIM) = IniVal
        BOUNDARY(ng) % vbar_south_Cx(XI_DIM) = IniVal
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
# if defined CELERITY_READ || defined CELERITY_WRITE
        BOUNDARY(ng) % u_south_C2(XI_DIM,1:N(ng)) = IniVal
        BOUNDARY(ng) % u_south_Ce(XI_DIM,1:N(ng)) = IniVal
        BOUNDARY(ng) % u_south_Cx(XI_DIM,1:N(ng)) = IniVal
        BOUNDARY(ng) % v_south_C2(XI_DIM,1:N(ng)) = IniVal
        BOUNDARY(ng) % v_south_Ce(XI_DIM,1:N(ng)) = IniVal
        BOUNDARY(ng) % v_south_Cx(XI_DIM,1:N(ng)) = IniVal
        BOUNDARY(ng) % t_south_C2(XI_DIM,1:N(ng),1:NT(ng)) = IniVal
        BOUNDARY(ng) % t_south_Ce(XI_DIM,1:N(ng),1:NT(ng)) = IniVal
        BOUNDARY(ng) % t_south_Cx(XI_DIM,1:N(ng),1:NT(ng)) = IniVal
# endif
#endif
      END IF
!
!  Northern boundary variables.
!
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
#if defined CELERITY_READ || defined CELERITY_WRITE
        BOUNDARY(ng) % zeta_north_C2(XI_DIM) = IniVal
        BOUNDARY(ng) % zeta_north_Ce(XI_DIM) = IniVal
        BOUNDARY(ng) % zeta_north_Cx(XI_DIM) = IniVal
        BOUNDARY(ng) % ubar_north_C2(XI_DIM) = IniVal
        BOUNDARY(ng) % ubar_north_Ce(XI_DIM) = IniVal
        BOUNDARY(ng) % ubar_north_Cx(XI_DIM) = IniVal
        BOUNDARY(ng) % vbar_north_C2(XI_DIM) = IniVal
        BOUNDARY(ng) % vbar_north_Ce(XI_DIM) = IniVal
        BOUNDARY(ng) % vbar_north_Cx(XI_DIM) = IniVal
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
# if defined CELERITY_READ || defined CELERITY_WRITE
        BOUNDARY(ng) % u_north_C2(XI_DIM,1:N(ng)) = IniVal
        BOUNDARY(ng) % u_north_Ce(XI_DIM,1:N(ng)) = IniVal
        BOUNDARY(ng) % u_north_Cx(XI_DIM,1:N(ng)) = IniVal
        BOUNDARY(ng) % v_north_C2(XI_DIM,1:N(ng)) = IniVal
        BOUNDARY(ng) % v_north_Ce(XI_DIM,1:N(ng)) = IniVal
        BOUNDARY(ng) % v_north_Cx(XI_DIM,1:N(ng)) = IniVal
        BOUNDARY(ng) % t_north_C2(XI_DIM,1:N(ng),1:NT(ng)) = IniVal
        BOUNDARY(ng) % t_north_Ce(XI_DIM,1:N(ng),1:NT(ng)) = IniVal
        BOUNDARY(ng) % t_north_Cx(XI_DIM,1:N(ng),1:NT(ng)) = IniVal
# endif
#endif
      END IF
      RETURN
      END SUBROUTINE initialize_boundary

      END MODULE mod_boundary
