#include "cppdefs.h"
      MODULE tl_obc_adjust_mod

#ifdef ADJUST_BOUNDARY
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine time-interpolates 4DVar tangent linear model open      !
!  boundary increments.  The increments can be constant (Nbrec=1)      !
!  or time interpolated between snapshots (Nbrec>1).                   !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng        Nested grid number.                                    !
!     tile      Domain partition.                                      !
!     Linp      4DVar increment time index to process.                 !
!                                                                      !
!=======================================================================
!
      implicit none
 
      PRIVATE
      PUBLIC :: tl_obc_adjust

      CONTAINS
!
!***********************************************************************
      SUBROUTINE tl_obc_adjust (ng, tile, Linp)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, Linp
!
!  Local variable declarations.
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iTLM, 7)
# endif
      CALL tl_obc_adjust_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj, LBij, UBij,          &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         Linp,                                    &
# ifdef WEST_FSOBC
     &                         BOUNDARY(ng) % tl_zeta_west,             &
# endif
# ifdef EAST_FSOBC
     &                         BOUNDARY(ng) % tl_zeta_east,             &
# endif
# ifdef SOUTH_FSOBC
     &                         BOUNDARY(ng) % tl_zeta_south,            &
# endif
# ifdef NORTH_FSOBC
     &                         BOUNDARY(ng) % tl_zeta_north,            &
# endif
# ifdef WEST_M2OBC
     &                         BOUNDARY(ng) % tl_ubar_west,             &
     &                         BOUNDARY(ng) % tl_vbar_west,             &
# endif
# ifdef EAST_M2OBC
     &                         BOUNDARY(ng) % tl_ubar_east,             &
     &                         BOUNDARY(ng) % tl_vbar_east,             &
# endif
# ifdef SOUTH_M2OBC
     &                         BOUNDARY(ng) % tl_ubar_south,            &
     &                         BOUNDARY(ng) % tl_vbar_south,            &
# endif
# ifdef NORTH_M2OBC
     &                         BOUNDARY(ng) % tl_ubar_north,            &
     &                         BOUNDARY(ng) % tl_vbar_north,            &
# endif
# ifdef SOLVE3D
#  ifdef WEST_M3OBC
     &                         BOUNDARY(ng) % tl_u_west,                &
     &                         BOUNDARY(ng) % tl_v_west,                &
#  endif
#  ifdef EAST_M3OBC
     &                         BOUNDARY(ng) % tl_u_east,                &
     &                         BOUNDARY(ng) % tl_v_east,                &
#  endif
#  ifdef SOUTH_M3OBC
     &                         BOUNDARY(ng) % tl_u_south,               &
     &                         BOUNDARY(ng) % tl_v_south,               &
#  endif
#  ifdef NORTH_M3OBC
     &                         BOUNDARY(ng) % tl_u_north,               &
     &                         BOUNDARY(ng) % tl_v_north,               &
#  endif
#  ifdef WEST_TOBC
     &                         BOUNDARY(ng) % tl_t_west,                &
#  endif
#  ifdef EAST_TOBC
     &                         BOUNDARY(ng) % tl_t_east,                &
#  endif
#  ifdef SOUTH_TOBC
     &                         BOUNDARY(ng) % tl_t_south,               &
#  endif
#  ifdef NORTH_TOBC
     &                         BOUNDARY(ng) % tl_t_north,               &
#  endif
# endif
# ifdef SOLVE3D
     &                         BOUNDARY(ng) % tl_t_obc,                 &
     &                         BOUNDARY(ng) % tl_u_obc,                 &
     &                         BOUNDARY(ng) % tl_v_obc,                 &
# endif
     &                         BOUNDARY(ng) % tl_ubar_obc,              &
     &                         BOUNDARY(ng) % tl_vbar_obc,              &
     &                         BOUNDARY(ng) % tl_zeta_obc)
# ifdef PROFILE
      CALL wclock_off (ng, iTLM, 7)
# endif
      RETURN
      END SUBROUTINE tl_obc_adjust
!
!***********************************************************************
      SUBROUTINE tl_obc_adjust_tile (ng, tile,                          &
     &                               LBi, UBi, LBj, UBj, LBij, UBij,    &
     &                               IminS, ImaxS, JminS, JmaxS,        &
     &                               Linp,                              &
# ifdef WEST_FSOBC
     &                               tl_zeta_west,                      &
# endif
# ifdef EAST_FSOBC
     &                               tl_zeta_east,                      &
# endif
# ifdef SOUTH_FSOBC
     &                               tl_zeta_south,                     &
# endif
# ifdef NORTH_FSOBC
     &                               tl_zeta_north,                     &
# endif
# ifdef WEST_M2OBC
     &                               tl_ubar_west, tl_vbar_west,        &
# endif
# ifdef EAST_M2OBC
     &                               tl_ubar_east, tl_vbar_east,        &
# endif
# ifdef SOUTH_M2OBC
     &                               tl_ubar_south, tl_vbar_south,      &
# endif
# ifdef NORTH_M2OBC
     &                               tl_ubar_north, tl_vbar_north,      &
# endif
# ifdef SOLVE3D
#  ifdef WEST_M3OBC
     &                               tl_u_west, tl_v_west,              &
#  endif
#  ifdef EAST_M3OBC
     &                               tl_u_east, tl_v_east,              &
#  endif
#  ifdef SOUTH_M3OBC
     &                               tl_u_south, tl_v_south,            &
#  endif
#  ifdef NORTH_M3OBC
     &                               tl_u_north, tl_v_north,            &
#  endif
#  ifdef WEST_TOBC
     &                               tl_t_west,                         &
#  endif
#  ifdef EAST_TOBC
     &                               tl_t_east,                         &
#  endif
#  ifdef SOUTH_TOBC
     &                               tl_t_south,                        &
#  endif
#  ifdef NORTH_TOBC
     &                               tl_t_north,                        &
#  endif
# endif
# ifdef SOLVE3D
     &                               tl_t_obc, tl_u_obc, tl_v_obc,      &
# endif
     &                               tl_ubar_obc, tl_vbar_obc,          &
     &                               tl_zeta_obc)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBij, UBij
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Linp
!
# ifdef ASSUMED_SHAPE
#  ifdef SOLVE3D
      real(r8), intent(in) :: tl_t_obc(LBij:,:,:,:,:,:)
      real(r8), intent(in) :: tl_u_obc(LBij:,:,:,:,:)
      real(r8), intent(in) :: tl_v_obc(LBij:,:,:,:,:)
#  endif
      real(r8), intent(in) :: tl_ubar_obc(LBij:,:,:,:)
      real(r8), intent(in) :: tl_vbar_obc(LBij:,:,:,:)
      real(r8), intent(in) :: tl_zeta_obc(LBij:,:,:,:)
#  ifdef WEST_FSOBC
      real(r8), intent(inout) :: tl_zeta_west(LBj:)
#  endif
#  ifdef EAST_FSOBC
      real(r8), intent(inout) :: tl_zeta_east(LBj:)
#  endif
#  ifdef SOUTH_FSOBC
      real(r8), intent(inout) :: tl_zeta_south(LBi:)
#  endif
#  ifdef NORTH_FSOBC
      real(r8), intent(inout) :: tl_zeta_north(LBi:)
#  endif
#  ifdef WEST_M2OBC
      real(r8), intent(inout) :: tl_ubar_west(LBj:)
      real(r8), intent(inout) :: tl_vbar_west(LBj:)
#  endif
#  ifdef EAST_M2OBC
      real(r8), intent(inout) :: tl_ubar_east(LBj:)
      real(r8), intent(inout) :: tl_vbar_east(LBj:)
#  endif
#  ifdef SOUTH_M2OBC
      real(r8), intent(inout) :: tl_ubar_south(LBi:)
      real(r8), intent(inout) :: tl_vbar_south(LBi:)
#  endif
#  ifdef NORTH_M2OBC
      real(r8), intent(inout) :: tl_ubar_north(LBi:)
      real(r8), intent(inout) :: tl_vbar_north(LBi:)
#  endif
#  ifdef SOLVE3D
#   ifdef WEST_M3OBC
      real(r8), intent(inout) :: tl_u_west(LBj:,:)
      real(r8), intent(inout) :: tl_v_west(LBj:,:)
#   endif
#   ifdef EAST_M3OBC
      real(r8), intent(inout) :: tl_u_east(LBj:,:)
      real(r8), intent(inout) :: tl_v_east(LBj:,:)
#   endif
#   ifdef SOUTH_M3OBC
      real(r8), intent(inout) :: tl_u_south(LBi:,:)
      real(r8), intent(inout) :: tl_v_south(LBi:,:)
#   endif
#   ifdef NORTH_M3OBC
      real(r8), intent(inout) :: tl_u_north(LBi:,:)
      real(r8), intent(inout) :: tl_v_north(LBi:,:)
#   endif
#   ifdef WEST_TOBC
      real(r8), intent(inout) :: tl_t_west(LBj:,:,:)
#   endif
#   ifdef EAST_TOBC
      real(r8), intent(inout) :: tl_t_east(LBj:,:,:)
#   endif
#   ifdef SOUTH_TOBC
      real(r8), intent(inout) :: tl_t_south(LBi:,:,:)
#   endif
#   ifdef NORTH_TOBC
      real(r8), intent(inout) :: tl_t_north(LBi:,:,:)
#   endif
#  endif

# else

#  ifdef SOLVE3D
      real(r8), intent(in) :: tl_t_obc(LBij:UBij,N(ng),4,               &
     &                                 Nbrec(ng),2,NT(ng))
      real(r8), intent(in) :: tl_u_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
      real(r8), intent(in) :: tl_v_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
#  endif
      real(r8), intent(in) :: tl_ubar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(in) :: tl_vbar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(in) :: tl_zeta_obc(LBij:UBij,4,Nbrec(ng),2)
#  ifdef WEST_FSOBC
      real(r8), intent(inout) :: tl_zeta_west(LBj:UBj)
#  endif
#  ifdef EAST_FSOBC
      real(r8), intent(inout) :: tl_zeta_east(LBj:UBj)
#  endif
#  ifdef SOUTH_FSOBC
      real(r8), intent(inout) :: tl_zeta_south(LBi:UBi)
#  endif
#  ifdef NORTH_FSOBC
      real(r8), intent(inout) :: tl_zeta_north(LBi:UBi)
#  endif
#  ifdef WEST_M2OBC
      real(r8), intent(inout) :: tl_ubar_west(LBj:UBj)
      real(r8), intent(inout) :: tl_vbar_west(LBj:UBj)
#  endif
#  ifdef EAST_M2OBC
      real(r8), intent(inout) :: tl_ubar_east(LBj:UBj)
      real(r8), intent(inout) :: tl_vbar_east(LBj:UBj)
#  endif
#  ifdef SOUTH_M2OBC
      real(r8), intent(inout) :: tl_ubar_south(LBi:UBi)
      real(r8), intent(inout) :: tl_vbar_south(LBi:UBi)
#  endif
#  ifdef NORTH_M2OBC
      real(r8), intent(inout) :: tl_ubar_north(LBi:UBi)
      real(r8), intent(inout) :: tl_vbar_north(LBi:UBi)
#  endif
#  ifdef SOLVE3D
#   ifdef WEST_M3OBC
      real(r8), intent(inout) :: tl_u_west(LBj:UBj,N(ng))
      real(r8), intent(inout) :: tl_v_west(LBj:UBj,N(ng))
#   endif
#   ifdef EAST_M3OBC
      real(r8), intent(inout) :: tl_u_east(LBj:UBj,N(ng))
      real(r8), intent(inout) :: tl_v_east(LBj:UBj,N(ng))
#   endif
#   ifdef SOUTH_M3OBC
      real(r8), intent(inout) :: tl_u_south(LBi:UBi,N(ng))
      real(r8), intent(inout) :: tl_v_south(LBi:UBi,N(ng))
#   endif
#   ifdef NORTH_M3OBC
      real(r8), intent(inout) :: tl_u_north(LBi:UBi,N(ng))
      real(r8), intent(inout) :: tl_v_north(LBi:UBi,N(ng))
#   endif
#   ifdef WEST_TOBC
      real(r8), intent(inout) :: tl_t_west(LBj:UBj,N(ng),NT(ng))
#   endif
#   ifdef EAST_TOBC
      real(r8), intent(inout) :: tl_t_east(LBj:UBj,N(ng),NT(ng))
#   endif
#   ifdef SOUTH_TOBC
      real(r8), intent(inout) :: tl_t_south(LBi:UBi,N(ng),NT(ng))
#   endif
#   ifdef NORTH_TOBC
      real(r8), intent(inout) :: tl_t_north(LBi:UBi,N(ng),NT(ng))
#   endif
#  endif
# endif
!
!  Local variable declarations.
!
      integer :: i, it1, it2, j
# ifdef SOLVE3D
      integer :: it, k
# endif
      real(r8) :: fac, fac1, fac2

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Adjust tangent linear open boundary fields with 4DVar increments.
!-----------------------------------------------------------------------
!
!  Set time records and interpolation factor, if any.
!
      IF (Nbrec(ng).eq.1) THEN
        it1=1
        it2=1
        fac1=1.0_r8
        fac2=0.0_r8
      ELSE
        it1=MAX(0,(iic(ng)-1)/nOBC(ng))+1
        it2=MIN(it1+1,Nbrec(ng))
        fac1=OBC_time(it2,ng)-(time(ng)+dt(ng))
        fac2=(time(ng)+dt(ng))-OBC_time(it1,ng)
        fac=1.0_r8/(fac1+fac2)
        fac1=fac*fac1
        fac2=fac*fac2
      END IF
!
!  Free-surface open boundaries.
!
      IF (ANY(Lobc(:,isFsur,ng))) THEN
# ifdef WEST_FSOBC
        IF ((Lobc(iwest,isFsur,ng)).and.WESTERN_EDGE) THEN
          DO j=Jstr,Jend
            tl_zeta_west(j)=fac1*tl_zeta_obc(j,iwest,it1,Linp)+         &
     &                      fac2*tl_zeta_obc(j,iwest,it2,Linp)
          END DO
        END IF
# endif
# ifdef EAST_FSOBC
        IF ((Lobc(ieast,isFsur,ng)).and.EASTERN_EDGE) THEN
          DO j=Jstr,Jend
            tl_zeta_east(j)=fac1*tl_zeta_obc(j,ieast,it1,Linp)+         &
     &                      fac2*tl_zeta_obc(j,ieast,it2,Linp)
          END DO
        END IF
# endif
# ifdef SOUTH_FSOBC
        IF ((Lobc(isouth,isFsur,ng)).and.SOUTHERN_EDGE) THEN
          DO i=Istr,Iend
            tl_zeta_south(i)=fac1*tl_zeta_obc(i,isouth,it1,Linp)+       &
     &                       fac2*tl_zeta_obc(i,isouth,it2,Linp)
          END DO
        END IF
# endif
# ifdef NORTH_FSOBC
        IF ((Lobc(inorth,isFsur,ng)).and.NORTHERN_EDGE) THEN
          DO i=Istr,Iend
            tl_zeta_north(i)=fac1*tl_zeta_obc(i,inorth,it1,Linp)+       &
     &                       fac2*tl_zeta_obc(i,inorth,it2,Linp)
          END DO
        END IF
# endif
      END IF
!
!  2D U-momentum open boundaries.
!
      IF (ANY(Lobc(:,isUbar,ng))) THEN
# ifdef WEST_M2OBC
        IF ((Lobc(iwest,isUbar,ng)).and.WESTERN_EDGE) THEN
          DO j=Jstr,Jend
            tl_ubar_west(j)=fac1*tl_ubar_obc(j,iwest,it1,Linp)+         &
     &                      fac2*tl_ubar_obc(j,iwest,it2,Linp)
          END DO
        END IF
# endif
# ifdef EAST_M2OBC
        IF ((Lobc(ieast,isUbar,ng)).and.EASTERN_EDGE) THEN
          DO j=Jstr,Jend
            tl_ubar_east(j)=fac1*tl_ubar_obc(j,ieast,it1,Linp)+         &
     &                      fac2*tl_ubar_obc(j,ieast,it2,Linp)
          END DO
        END IF
# endif
# ifdef SOUTH_M2OBC
        IF ((Lobc(isouth,isUbar,ng)).and.SOUTHERN_EDGE) THEN
          DO i=IstrU,Iend
            tl_ubar_south(i)=fac1*tl_ubar_obc(i,isouth,it1,Linp)+       &
     &                       fac2*tl_ubar_obc(i,isouth,it2,Linp)
          END DO
        END IF
# endif
# ifdef NORTH_M2OBC
        IF ((Lobc(inorth,isUbar,ng)).and.NORTHERN_EDGE) THEN
          DO i=IstrU,Iend
            tl_ubar_north(i)=fac1*tl_ubar_obc(i,inorth,it1,Linp)+       &
     &                       fac2*tl_ubar_obc(i,inorth,it2,Linp)
          END DO
        END IF
# endif
      END IF
!
!  2D V-momentum open boundaries.
!
      IF (ANY(Lobc(:,isVbar,ng))) THEN
# ifdef WEST_M2OBC
        IF ((Lobc(iwest,isVbar,ng)).and.WESTERN_EDGE) THEN
          DO j=JstrV,Jend
            tl_vbar_west(j)=fac1*tl_vbar_obc(j,iwest,it1,Linp)+         &
     &                      fac2*tl_vbar_obc(j,iwest,it2,Linp)
          END DO
        END IF
# endif
# ifdef EAST_M2OBC
        IF ((Lobc(ieast,isVbar,ng)).and.EASTERN_EDGE) THEN
          DO j=JstrV,Jend
            tl_vbar_east(j)=fac1*tl_vbar_obc(j,ieast,it1,Linp)+         &
     &                      fac2*tl_vbar_obc(j,ieast,it2,Linp)
          END DO
        END IF
# endif
# ifdef SOUTH_M2OBC
        IF ((Lobc(isouth,isVbar,ng)).and.SOUTHERN_EDGE) THEN
          DO i=Istr,Iend
            tl_vbar_south(i)=fac1*tl_vbar_obc(i,isouth,it1,Linp)+       &
     &                       fac2*tl_vbar_obc(i,isouth,it2,Linp)
          END DO
        END IF
# endif
# ifdef NORTH_M2OBC
        IF ((Lobc(inorth,isVbar,ng)).and.NORTHERN_EDGE) THEN
          DO i=Istr,Iend
            tl_vbar_north(i)=fac1*tl_vbar_obc(i,inorth,it1,Linp)+       &
     &                       fac2*tl_vbar_obc(i,inorth,it2,Linp)
          END DO
        END IF
# endif
      END IF

# ifdef SOLVE3D
!
!  3D U-momentum open boundaries.
!
      IF (ANY(Lobc(:,isUvel,ng))) THEN
#  ifdef WEST_M3OBC
        IF ((Lobc(iwest,isUvel,ng)).and.WESTERN_EDGE) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              tl_u_west(j,k)=fac1*tl_u_obc(j,k,iwest,it1,Linp)+         &
     &                       fac2*tl_u_obc(j,k,iwest,it2,Linp)
            END DO
          END DO
        END IF
#  endif
#  ifdef EAST_M3OBC
        IF ((Lobc(ieast,isUvel,ng)).and.EASTERN_EDGE) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              tl_u_east(j,k)=fac1*tl_u_obc(j,k,ieast,it1,Linp)+         &
     &                       fac2*tl_u_obc(j,k,ieast,it2,Linp)
            END DO
          END DO
        END IF
#  endif
#  ifdef SOUTH_M3OBC
        IF ((Lobc(isouth,isUvel,ng)).and.SOUTHERN_EDGE) THEN
          DO k=1,N(ng)
            DO i=IstrU,Iend
              tl_u_south(i,k)=fac1*tl_u_obc(i,k,isouth,it1,Linp)+       &
     &                        fac2*tl_u_obc(i,k,isouth,it2,Linp)
            END DO
          END DO
        END IF
#  endif
#  ifdef NORTH_M3OBC
        IF ((Lobc(inorth,isUvel,ng)).and.NORTHERN_EDGE) THEN
          DO k=1,N(ng)
            DO i=IstrU,Iend
              tl_u_north(i,k)=fac1*tl_u_obc(i,k,inorth,it1,Linp)+       &
     &                        fac2*tl_u_obc(i,k,inorth,it2,Linp)
            END DO
          END DO
        END IF
#  endif
      END IF
!
!  3D V-momentum open boundaries.
!
      IF (ANY(Lobc(:,isVvel,ng))) THEN
#  ifdef WEST_M3OBC
        IF ((Lobc(iwest,isVvel,ng)).and.WESTERN_EDGE) THEN
          DO k=1,N(ng)
            DO j=JstrV,Jend
              tl_v_west(j,k)=fac1*tl_v_obc(j,k,iwest,it1,Linp)+         &
     &                       fac2*tl_v_obc(j,k,iwest,it2,Linp)
            END DO
          END DO
        END IF
#  endif
#  ifdef EAST_M3OBC
        IF ((Lobc(ieast,isVvel,ng)).and.EASTERN_EDGE) THEN
          DO k=1,N(ng)
            DO j=JstrV,Jend
              tl_v_east(j,k)=fac1*tl_v_obc(j,k,ieast,it1,Linp)+         &
     &                       fac2*tl_v_obc(j,k,ieast,it2,Linp)
            END DO
          END DO
        END IF
#  endif
#  ifdef SOUTH_M3OBC
        IF ((Lobc(isouth,isVvel,ng)).and.SOUTHERN_EDGE) THEN
          DO k=1,N(ng)
            DO i=Istr,Iend
              tl_v_south(i,k)=fac1*tl_v_obc(j,k,isouth,it1,Linp)+       &
     &                        fac2*tl_v_obc(j,k,isouth,it2,Linp)
            END DO
          END DO
        END IF
#  endif
#  ifdef NORTH_M3OBC
        IF ((Lobc(inorth,isVvel,ng)).and.NORTHERN_EDGE) THEN
          DO k=1,N(ng)
            DO i=Istr,Iend
              tl_v_north(i,k)=fac1*tl_v_obc(j,k,inorth,it1,Linp)+       &
     &                        fac2*tl_v_obc(j,k,inorth,it2,Linp)
            END DO
          END DO
        END IF
#  endif
      END IF
!
!  Tracers open boundaries.
!
      DO it=1,NT(ng)
        IF (ANY(Lobc(:,isTvar(it),ng))) THEN
#  ifdef WEST_TOBC
          IF ((Lobc(iwest,isTvar(it),ng)).and.WESTERN_EDGE) THEN
            DO k=1,N(ng)
              DO j=Jstr,Jend
                tl_t_west(j,k,it)=fac1*                                 &
     &                            tl_t_obc(j,k,iwest,it1,Linp,it)+      &
     &                            fac2*                                 &
     &                            tl_t_obc(j,k,iwest,it2,Linp,it)
              END DO
            END DO
          END IF
#  endif
#  ifdef EAST_TOBC
          IF ((Lobc(ieast,isTvar(it),ng)).and.EASTERN_EDGE) THEN
            DO k=1,N(ng)
              DO j=Jstr,Jend
                tl_t_east(j,k,it)=fac1*                                 &
     &                            tl_t_obc(j,k,ieast,it1,Linp,it)+      &
     &                            fac2*                                 &
     &                            tl_t_obc(j,k,ieast,it2,Linp,it)
              END DO
            END DO
          END IF
#  endif
#  ifdef SOUTH_TOBC
          IF ((Lobc(isouth,isTvar(it),ng)).and.SOUTHERN_EDGE) THEN
            DO k=1,N(ng)
              DO i=Istr,Iend
                tl_t_south(i,k,it)=fac1*                                &
     &                             tl_t_obc(i,k,isouth,it1,Linp,it)+    &
     &                             fac2*                                &
     &                             tl_t_obc(i,k,isouth,it2,Linp,it)
              END DO
            END DO
          END IF
#  endif
#  ifdef NORTH_TOBC
          IF ((Lobc(inorth,isTvar(it),ng)).and.NORTHERN_EDGE) THEN
            DO k=1,N(ng)
              DO i=Istr,Iend
                tl_t_north(i,k,it)=fac1*                                &
     &                             tl_t_obc(i,k,inorth,it1,Linp,it)+    &
     &                             fac2*                                &
     &                             tl_t_obc(i,k,inorth,it2,Linp,it)
              END DO
            END DO
          END IF
#  endif
        END IF
      END DO
# endif

      RETURN
      END SUBROUTINE tl_obc_adjust_tile
#endif
      END MODULE tl_obc_adjust_mod
