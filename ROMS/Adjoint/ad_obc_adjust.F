#include "cppdefs.h"
      MODULE ad_obc_adjust_mod

#ifdef ADJUST_BOUNDARY
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine is the adjoint of time-interpolation of 4DVar open     !
!  boundary increments.  The increments can be constant  (Nbrec=1)     !
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
      PUBLIC :: ad_obc_adjust

      CONTAINS
!
!***********************************************************************
      SUBROUTINE ad_obc_adjust (ng, tile, Linp)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
      USE mod_forces
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
      CALL wclock_on (ng, iADM, 7)
# endif
      CALL ad_obc_adjust_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj, LBij, UBij,          &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         Linp,                                    &
# ifdef WEST_FSOBC
     &                         BOUNDARY(ng) % ad_zeta_west,             &
# endif
# ifdef EAST_FSOBC
     &                         BOUNDARY(ng) % ad_zeta_east,             &
# endif
# ifdef SOUTH_FSOBC
     &                         BOUNDARY(ng) % ad_zeta_south,            &
# endif
# ifdef NORTH_FSOBC
     &                         BOUNDARY(ng) % ad_zeta_north,            &
# endif
# ifdef WEST_M2OBC
     &                         BOUNDARY(ng) % ad_ubar_west,             &
     &                         BOUNDARY(ng) % ad_vbar_west,             &
# endif
# ifdef EAST_M2OBC
     &                         BOUNDARY(ng) % ad_ubar_east,             &
     &                         BOUNDARY(ng) % ad_vbar_east,             &
# endif
# ifdef SOUTH_M2OBC
     &                         BOUNDARY(ng) % ad_ubar_south,            &
     &                         BOUNDARY(ng) % ad_vbar_south,            &
# endif
# ifdef NORTH_M2OBC
     &                         BOUNDARY(ng) % ad_ubar_north,            &
     &                         BOUNDARY(ng) % ad_vbar_north,            &
# endif
# ifdef SOLVE3D
#  ifdef WEST_M3OBC
     &                         BOUNDARY(ng) % ad_u_west,                &
     &                         BOUNDARY(ng) % ad_v_west,                &
#  endif
#  ifdef EAST_M3OBC
     &                         BOUNDARY(ng) % ad_u_east,                &
     &                         BOUNDARY(ng) % ad_v_east,                &
#  endif
#  ifdef SOUTH_M3OBC
     &                         BOUNDARY(ng) % ad_u_south,               &
     &                         BOUNDARY(ng) % ad_v_south,               &
#  endif
#  ifdef NORTH_M3OBC
     &                         BOUNDARY(ng) % ad_u_north,               &
     &                         BOUNDARY(ng) % ad_v_north,               &
#  endif
#  ifdef WEST_TOBC
     &                         BOUNDARY(ng) % ad_t_west,                &
#  endif
#  ifdef EAST_TOBC
     &                         BOUNDARY(ng) % ad_t_east,                &
#  endif
#  ifdef SOUTH_TOBC
     &                         BOUNDARY(ng) % ad_t_south,               &
#  endif
#  ifdef NORTH_TOBC
     &                         BOUNDARY(ng) % ad_t_north,               &
#  endif
# endif
# ifdef SOLVE3D
     &                         BOUNDARY(ng) % ad_t_obc,                 &
     &                         BOUNDARY(ng) % ad_u_obc,                 &
     &                         BOUNDARY(ng) % ad_v_obc,                 &
# endif
     &                         BOUNDARY(ng) % ad_ubar_obc,              &
     &                         BOUNDARY(ng) % ad_vbar_obc,              &
     &                         BOUNDARY(ng) % ad_zeta_obc)
# ifdef PROFILE
      CALL wclock_off (ng, iADM, 7)
# endif

      RETURN
      END SUBROUTINE ad_obc_adjust
!
!***********************************************************************
      SUBROUTINE ad_obc_adjust_tile (ng, tile,                          &
     &                               LBi, UBi, LBj, UBj, LBij, UBij,    &
     &                               IminS, ImaxS, JminS, JmaxS,        &
     &                               Linp,                              &
# ifdef WEST_FSOBC
     &                               ad_zeta_west,                      &
# endif
# ifdef EAST_FSOBC
     &                               ad_zeta_east,                      &
# endif
# ifdef SOUTH_FSOBC
     &                               ad_zeta_south,                     &
# endif
# ifdef NORTH_FSOBC
     &                               ad_zeta_north,                     &
# endif
# ifdef WEST_M2OBC
     &                               ad_ubar_west, ad_vbar_west,        &
# endif
# ifdef EAST_M2OBC
     &                               ad_ubar_east, ad_vbar_east,        &
# endif
# ifdef SOUTH_M2OBC
     &                               ad_ubar_south, ad_vbar_south,      &
# endif
# ifdef NORTH_M2OBC
     &                               ad_ubar_north, ad_vbar_north,      &
# endif
# ifdef SOLVE3D
#  ifdef WEST_M3OBC
     &                               ad_u_west, ad_v_west,              &
#  endif
#  ifdef EAST_M3OBC
     &                               ad_u_east, ad_v_east,              &
#  endif
#  ifdef SOUTH_M3OBC
     &                               ad_u_south, ad_v_south,            &
#  endif
#  ifdef NORTH_M3OBC
     &                               ad_u_north, ad_v_north,            &
#  endif
#  ifdef WEST_TOBC
     &                               ad_t_west,                         &
#  endif
#  ifdef EAST_TOBC
     &                               ad_t_east,                         &
#  endif
#  ifdef SOUTH_TOBC
     &                               ad_t_south,                        &
#  endif
#  ifdef NORTH_TOBC
     &                               ad_t_north,                        &
#  endif
# endif
# ifdef SOLVE3D
     &                               ad_t_obc, ad_u_obc, ad_v_obc,      &
# endif
     &                               ad_ubar_obc, ad_vbar_obc,          &
     &                               ad_zeta_obc)
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
#  ifdef WEST_FSOBC
      real(r8), intent(in) :: ad_zeta_west(LBj:)
#  endif
#  ifdef EAST_FSOBC
      real(r8), intent(in) :: ad_zeta_east(LBj:)
#  endif
#  ifdef SOUTH_FSOBC
      real(r8), intent(in) :: ad_zeta_south(LBi:)
#  endif
#  ifdef NORTH_FSOBC
      real(r8), intent(in) :: ad_zeta_north(LBi:)
#  endif
#  ifdef WEST_M2OBC
      real(r8), intent(in) :: ad_ubar_west(LBj:)
      real(r8), intent(in) :: ad_vbar_west(LBj:)
#  endif
#  ifdef EAST_M2OBC
      real(r8), intent(in) :: ad_ubar_east(LBj:)
      real(r8), intent(in) :: ad_vbar_east(LBj:)
#  endif
#  ifdef SOUTH_M2OBC
      real(r8), intent(in) :: ad_ubar_south(LBi:)
      real(r8), intent(in) :: ad_vbar_south(LBi:)
#  endif
#  ifdef NORTH_M2OBC
      real(r8), intent(in) :: ad_ubar_north(LBi:)
      real(r8), intent(in) :: ad_vbar_north(LBi:)
#  endif
#  ifdef SOLVE3D
#   ifdef WEST_M3OBC
      real(r8), intent(in) :: ad_u_west(LBj:,:)
      real(r8), intent(in) :: ad_v_west(LBj:,:)
#   endif
#   ifdef EAST_M3OBC
      real(r8), intent(in) :: ad_u_east(LBj:,:)
      real(r8), intent(in) :: ad_v_east(LBj:,:)
#   endif
#   ifdef SOUTH_M3OBC
      real(r8), intent(in) :: ad_u_south(LBi:,:)
      real(r8), intent(in) :: ad_v_south(LBi:,:)
#   endif
#   ifdef NORTH_M3OBC
      real(r8), intent(in) :: ad_u_north(LBi:,:)
      real(r8), intent(in) :: ad_v_north(LBi:,:)
#   endif
#   ifdef WEST_TOBC
      real(r8), intent(in) :: ad_t_west(LBj:,:,:)
#   endif
#   ifdef EAST_TOBC
      real(r8), intent(in) :: ad_t_east(LBj:,:,:)
#   endif
#   ifdef SOUTH_TOBC
      real(r8), intent(in) :: ad_t_south(LBi:,:,:)
#   endif
#   ifdef NORTH_TOBC
      real(r8), intent(in) :: ad_t_north(LBi:,:,:)
#   endif
#  endif
#  ifdef SOLVE3D
      real(r8), intent(inout) :: ad_t_obc(LBij:,:,:,:,:,:)
      real(r8), intent(inout) :: ad_u_obc(LBij:,:,:,:,:)
      real(r8), intent(inout) :: ad_v_obc(LBij:,:,:,:,:)
#  endif
      real(r8), intent(inout) :: ad_ubar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: ad_vbar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: ad_zeta_obc(LBij:,:,:,:)

# else

#  ifdef WEST_FSOBC
      real(r8), intent(in) :: ad_zeta_west(LBj:UBj)
#  endif
#  ifdef EAST_FSOBC
      real(r8), intent(in) :: ad_zeta_east(LBj:UBj)
#  endif
#  ifdef SOUTH_FSOBC
      real(r8), intent(in) :: ad_zeta_south(LBi:UBi)
#  endif
#  ifdef NORTH_FSOBC
      real(r8), intent(in) :: ad_zeta_north(LBi:UBi)
#  endif
#  ifdef WEST_M2OBC
      real(r8), intent(in) :: ad_ubar_west(LBj:UBj)
      real(r8), intent(in) :: ad_vbar_west(LBj:UBj)
#  endif
#  ifdef EAST_M2OBC
      real(r8), intent(in) :: ad_ubar_east(LBj:UBj)
      real(r8), intent(in) :: ad_vbar_east(LBj:UBj)
#  endif
#  ifdef SOUTH_M2OBC
      real(r8), intent(in) :: ad_ubar_south(LBi:UBi)
      real(r8), intent(in) :: ad_vbar_south(LBi:UBi)
#  endif
#  ifdef NORTH_M2OBC
      real(r8), intent(in) :: ad_ubar_north(LBi:UBi)
      real(r8), intent(in) :: ad_vbar_north(LBi:UBi)
#  endif
#  ifdef SOLVE3D
#   ifdef WEST_M3OBC
      real(r8), intent(in) :: ad_u_west(LBj:UBj,N(ng))
      real(r8), intent(in) :: ad_v_west(LBj:UBj,N(ng))
#   endif
#   ifdef EAST_M3OBC
      real(r8), intent(in) :: ad_u_east(LBj:UBj,N(ng))
      real(r8), intent(in) :: ad_v_east(LBj:UBj,N(ng))
#   endif
#   ifdef SOUTH_M3OBC
      real(r8), intent(in) :: ad_u_south(LBi:UBi,N(ng))
      real(r8), intent(in) :: ad_v_south(LBi:UBi,N(ng))
#   endif
#   ifdef NORTH_M3OBC
      real(r8), intent(in) :: ad_u_north(LBi:UBi,N(ng))
      real(r8), intent(in) :: ad_v_north(LBi:UBi,N(ng))
#   endif
#   ifdef WEST_TOBC
      real(r8), intent(in) :: ad_t_west(LBj:UBj,N(ng),NT(ng))
#   endif
#   ifdef EAST_TOBC
      real(r8), intent(in) :: ad_t_east(LBj:UBj,N(ng),NT(ng))
#   endif
#   ifdef SOUTH_TOBC
      real(r8), intent(in) :: ad_t_south(LBi:UBi,N(ng),NT(ng))
#   endif
#   ifdef NORTH_TOBC
      real(r8), intent(in) :: ad_t_north(LBi:UBi,N(ng),NT(ng))
#   endif
#  endif
#  ifdef SOLVE3D
      real(r8), intent(inout) :: ad_t_obc(LBij:UBij,N(ng),4,            &
     &                                    Nbrec(ng),2,NT(ng))
      real(r8), intent(inout) :: ad_u_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
      real(r8), intent(inout) :: ad_v_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
#  endif
      real(r8), intent(inout) :: ad_ubar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: ad_vbar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: ad_zeta_obc(LBij:UBij,4,Nbrec(ng),2)
# endif
!
!  Local variable declarations.
!
      integer :: i, ib, it1, it2, j
# ifdef SOLVE3D
      integer :: it, k
# endif
      real(r8) :: fac, fac1, fac2

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Adjoint of adjust open boundary fields with 4DVAR increments.
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
          ib=iwest
          DO j=Jstr,Jend
!>          tl_zeta_west(j)=fac1*tl_zeta_obc(j,ib,it1,Linp)+            &
!>   &                      fac2*tl_zeta_obc(j,ib,it2,Linp)
!>
            ad_zeta_obc(j,ib,it1,Linp)=ad_zeta_obc(j,ib,it1,Linp)+      &
     &                                 fac1*ad_zeta_west(j)
            ad_zeta_obc(j,ib,it2,Linp)=ad_zeta_obc(j,ib,it2,Linp)+      &
     &                                 fac2*ad_zeta_west(j)
          END DO
        END IF
# endif
# ifdef EAST_FSOBC
        IF ((Lobc(ieast,isFsur,ng)).and.EASTERN_EDGE) THEN
          ib=ieast
          DO j=Jstr,Jend
!>          tl_zeta_east(j)=fac1*tl_zeta_obc(j,ib,it1,Linp)+            &
!>   &                      fac2*tl_zeta_obc(j,ib,it2,Linp)
!>
            ad_zeta_obc(j,ib,it1,Linp)=ad_zeta_obc(j,ib,it1,Linp)+      &
     &                                 fac1*ad_zeta_east(j)
            ad_zeta_obc(j,ib,it2,Linp)=ad_zeta_obc(j,ib,it2,Linp)+      &
     &                                 fac2*ad_zeta_east(j)
          END DO
        END IF
# endif
# ifdef SOUTH_FSOBC
        IF ((Lobc(isouth,isFsur,ng)).and.SOUTHERN_EDGE) THEN
          ib=isouth
          DO i=Istr,Iend
!>          tl_zeta_south(i)=fac1*tl_zeta_obc(i,ib,it1,Linp)+           &
!>   &                       fac2*tl_zeta_obc(i,ib,it2,Linp)
!>
            ad_zeta_obc(i,ib,it1,Linp)=ad_zeta_obc(i,ib,it1,Linp)+      &
     &                                 fac1*ad_zeta_south(i)
            ad_zeta_obc(i,ib,it2,Linp)=ad_zeta_obc(i,ib,it2,Linp)+      &
     &                                 fac2*ad_zeta_south(i)
          END DO
        END IF
# endif
# ifdef NORTH_FSOBC
        IF ((Lobc(inorth,isFsur,ng)).and.NORTHERN_EDGE) THEN
          ib=inorth
          DO i=Istr,Iend
!>          tl_zeta_north(i)=fac1*tl_zeta_obc(i,ib,it1,Linp)+           &
!>   &                       fac2*tl_zeta_obc(i,ib,it2,Linp)
!>
            ad_zeta_obc(i,ib,it1,Linp)=ad_zeta_obc(i,ib,it1,Linp)+      &
     &                                 fac1*ad_zeta_north(i)
            ad_zeta_obc(i,ib,it2,Linp)=ad_zeta_obc(i,ib,it2,Linp)+      &
     &                                 fac2*ad_zeta_north(i)
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
          ib=iwest
          DO j=Jstr,Jend
!>          tl_ubar_west(j)=fac1*tl_ubar_obc(j,ib,it1,Linp)+            &
!>   &                      fac2*tl_ubar_obc(j,ib,it2,Linp)
!>
            ad_ubar_obc(j,ib,it1,Linp)=ad_ubar_obc(j,ib,it1,Linp)+      &
     &                                 fac1*ad_ubar_west(j)
            ad_ubar_obc(j,ib,it2,Linp)=ad_ubar_obc(j,ib,it2,Linp)+      &
     &                                 fac2*ad_ubar_west(j)
          END DO
        END IF
# endif
# ifdef EAST_M2OBC
        IF ((Lobc(ieast,isUbar,ng)).and.EASTERN_EDGE) THEN
          ib=ieast
          DO j=Jstr,Jend
!>          tl_ubar_east(j)=fac1*tl_ubar_obc(j,ib,it1,Linp)+            &
!>   &                      fac2*tl_ubar_obc(j,ib,it2,Linp)
!>
            ad_ubar_obc(j,ib,it1,Linp)=ad_ubar_obc(j,ib,it1,Linp)+      &
     &                                 fac1*ad_ubar_east(j)
            ad_ubar_obc(j,ib,it2,Linp)=ad_ubar_obc(j,ib,it2,Linp)+      &
     &                                 fac2*ad_ubar_east(j)
          END DO
        END IF
# endif
# ifdef SOUTH_M2OBC
        IF ((Lobc(isouth,isUbar,ng)).and.SOUTHERN_EDGE) THEN
          ib=isouth
          DO i=IstrU,Iend
!>          tl_ubar_south(i)=fac1*tl_ubar_obc(i,ib,it1,Linp)+           &
!>   &                       fac2*tl_ubar_obc(i,ib,it2,Linp)
!>
            ad_ubar_obc(i,ib,it1,Linp)=ad_ubar_obc(i,ib,it1,Linp)+      &
     &                                 fac1*ad_ubar_south(i)
            ad_ubar_obc(i,ib,it2,Linp)=ad_ubar_obc(i,ib,it2,Linp)+      &
     &                                 fac2*ad_ubar_south(i)
          END DO
        END IF
# endif
# ifdef NORTH_M2OBC
        IF ((Lobc(inorth,isUbar,ng)).and.NORTHERN_EDGE) THEN
          ib=inorth
          DO i=IstrU,Iend
!>          tl_ubar_north(i)=fac1*tl_ubar_obc(i,ib,it1,Linp)+           &
!>   &                       fac2*tl_ubar_obc(i,ib,it2,Linp)
!>
            ad_ubar_obc(i,ib,it1,Linp)=ad_ubar_obc(i,ib,it1,Linp)+      &
     &                                 fac1*ad_ubar_north(i)
            ad_ubar_obc(i,ib,it2,Linp)=ad_ubar_obc(i,ib,it2,Linp)+      &
     &                                 fac2*ad_ubar_north(i)
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
          ib=iwest
          DO j=JstrV,Jend
!>          tl_vbar_west(j)=fac1*tl_vbar_obc(j,ib,it1,Linp)+            &
!>   &                      fac2*tl_vbar_obc(j,ib,it2,Linp)
!>
            ad_vbar_obc(j,ib,it1,Linp)=ad_vbar_obc(j,ib,it1,Linp)+      &
     &                                 fac1*ad_vbar_west(j)
            ad_vbar_obc(j,ib,it2,Linp)=ad_vbar_obc(j,ib,it2,Linp)+      &
     &                                 fac2*ad_vbar_west(j)
          END DO
        END IF
# endif
# ifdef EAST_M2OBC
        IF ((Lobc(ieast,isVbar,ng)).and.EASTERN_EDGE) THEN
          ib=ieast
          DO j=JstrV,Jend
!>          tl_vbar_east(j)=fac1*tl_vbar_obc(j,ib,it1,Linp)+            &
!>   &                      fac2*tl_vbar_obc(j,ib,it2,Linp)
!>
            ad_vbar_obc(j,ib,it1,Linp)=ad_vbar_obc(j,ib,it1,Linp)+      &
     &                                 fac1*ad_vbar_east(j)
            ad_vbar_obc(j,ib,it2,Linp)=ad_vbar_obc(j,ib,it2,Linp)+      &
     &                                 fac2*ad_vbar_east(j)
          END DO
        END IF
# endif
# ifdef SOUTH_M2OBC
        IF ((Lobc(isouth,isVbar,ng)).and.SOUTHERN_EDGE) THEN
          ib=isouth
          DO i=Istr,Iend
!>          tl_vbar_south(i)=fac1*tl_vbar_obc(i,ib,it1,Linp)+           &
!>   &                       fac2*tl_vbar_obc(i,ib,it2,Linp)
!>
            ad_vbar_obc(i,ib,it1,Linp)=ad_vbar_obc(i,ib,it1,Linp)+      &
     &                                 fac1*ad_vbar_south(i)
            ad_vbar_obc(i,ib,it2,Linp)=ad_vbar_obc(i,ib,it2,Linp)+      &
     &                                 fac2*ad_vbar_south(i)
          END DO
        END IF
# endif
# ifdef NORTH_M2OBC
        IF ((Lobc(inorth,isVbar,ng)).and.NORTHERN_EDGE) THEN
          ib=inorth
          DO i=Istr,Iend
!>          tl_vbar_north(i)=fac1*tl_vbar_obc(i,ib,it1,Linp)+           &
!>   &                       fac2*tl_vbar_obc(i,ib,it2,Linp)
!>
            ad_vbar_obc(i,ib,it1,Linp)=ad_vbar_obc(i,ib,it1,Linp)+      &
     &                                 fac1*ad_vbar_north(i)
            ad_vbar_obc(i,ib,it2,Linp)=ad_vbar_obc(i,ib,it2,Linp)+      &
     &                                 fac2*ad_vbar_north(i)
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
          ib=iwest
          DO k=1,N(ng)
            DO j=Jstr,Jend
!>            tl_u_west(j,k)=fac1*tl_u_obc(j,k,ib,it1,Linp)+            &
!>   &                       fac2*tl_u_obc(j,k,ib,it2,Linp)
!>
              ad_u_obc(j,k,ib,it1,Linp)=ad_u_obc(j,k,ib,it1,Linp)+      &
     &                                  fac1*ad_u_west(j,k)
              ad_u_obc(j,k,ib,it2,Linp)=ad_u_obc(j,k,ib,it2,Linp)+      &
     &                                  fac2*ad_u_west(j,k)
            END DO
          END DO
        END IF
#  endif
#  ifdef EAST_M3OBC
        IF ((Lobc(ieast,isUvel,ng)).and.EASTERN_EDGE) THEN
          ib=ieast
          DO k=1,N(ng)
            DO j=Jstr,Jend
!>            tl_u_east(j,k)=fac1*tl_u_obc(j,k,ib,it1,Linp)+            &
!>   &                       fac2*tl_u_obc(j,k,ib,it2,Linp)
!>
              ad_u_obc(j,k,ib,it1,Linp)=ad_u_obc(j,k,ib,it1,Linp)+      &
     &                                  fac1*ad_u_east(j,k)
              ad_u_obc(j,k,ib,it2,Linp)=ad_u_obc(j,k,ib,it2,Linp)+      &
     &                                  fac2*ad_u_east(j,k)
            END DO
          END DO
        END IF
#  endif
#  ifdef SOUTH_M3OBC
        IF ((Lobc(isouth,isUvel,ng)).and.SOUTHERN_EDGE) THEN
          ib=isouth
          DO k=1,N(ng)
            DO i=IstrU,Iend
!>            tl_u_south(i,k)=fac1*tl_u_obc(i,k,ib,it1,Linp)+           &
!>   &                        fac2*tl_u_obc(i,k,ib,it2,Linp)
!>
              ad_u_obc(i,k,ib,it1,Linp)=ad_u_obc(i,k,ib,it1,Linp)+      &
     &                                  fac1*ad_u_south(i,k)
              ad_u_obc(i,k,ib,it2,Linp)=ad_u_obc(i,k,ib,it2,Linp)+      &
     &                                  fac2*ad_u_south(i,k)
            END DO
          END DO
        END IF
#  endif
#  ifdef NORTH_M3OBC
        IF ((Lobc(inorth,isUvel,ng)).and.NORTHERN_EDGE) THEN
          ib=inorth
          DO k=1,N(ng)
            DO i=IstrU,Iend
!>            tl_u_north(i,k)=fac1*tl_u_obc(i,k,ib,it1,Linp)+           &
!>   &                        fac2*tl_u_obc(i,k,ib,it2,Linp)
!>
              ad_u_obc(i,k,ib,it1,Linp)=ad_u_obc(i,k,ib,it1,Linp)+      &
     &                                  fac1*ad_u_north(i,k)
              ad_u_obc(i,k,ib,it2,Linp)=ad_u_obc(i,k,ib,it2,Linp)+      &
     &                                  fac2*ad_u_north(i,k)
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
          ib=iwest
          DO k=1,N(ng)
            DO j=JstrV,Jend
!>            tl_v_west(j,k)=fac1*tl_v_obc(j,k,ib,it1,Linp)+            &
!>   &                       fac2*tl_v_obc(j,k,ib,it2,Linp)
!>
              ad_v_obc(j,k,ib,it1,Linp)=ad_v_obc(j,k,ib,it1,Linp)+      &
     &                                  fac1*ad_v_west(j,k)
              ad_v_obc(j,k,ib,it2,Linp)=ad_v_obc(j,k,ib,it2,Linp)+      &
     &                                  fac2*ad_v_west(j,k)
            END DO
          END DO
        END IF
#  endif
#  ifdef EAST_M3OBC
        IF ((Lobc(ieast,isVvel,ng)).and.EASTERN_EDGE) THEN
          ib=ieast
          DO k=1,N(ng)
            DO j=JstrV,Jend
!>            tl_v_east(j,k)=fac1*tl_v_obc(j,k,ib,it1,Linp)+            &
!>   &                       fac2*tl_v_obc(j,k,ib,it2,Linp)
!>
              ad_v_obc(j,k,ib,it1,Linp)=ad_v_obc(j,k,ib,it1,Linp)+      &
     &                                  fac1*ad_v_east(j,k)
              ad_v_obc(j,k,ib,it2,Linp)=ad_v_obc(j,k,ib,it2,Linp)+      &
     &                                  fac2*ad_v_east(j,k)
            END DO
          END DO
        END IF
#  endif
#  ifdef SOUTH_M3OBC
        IF ((Lobc(isouth,isVvel,ng)).and.SOUTHERN_EDGE) THEN
          ib=isouth
          DO k=1,N(ng)
            DO i=Istr,Iend
!>            tl_v_south(i,k)=fac1*tl_v_obc(j,k,ib,it1,Linp)+           &
!>   &                        fac2*tl_v_obc(j,k,ib,it2,Linp)
!>
              ad_v_obc(i,k,ib,it1,Linp)=ad_v_obc(i,k,ib,it1,Linp)+      &
     &                                  fac1*ad_v_south(i,k)
              ad_v_obc(i,k,ib,it2,Linp)=ad_v_obc(i,k,ib,it2,Linp)+      &
     &                                  fac2*ad_v_south(i,k)
            END DO
          END DO
        END IF
#  endif
#  ifdef NORTH_M3OBC
        IF ((Lobc(inorth,isVvel,ng)).and.NORTHERN_EDGE) THEN
          ib=inorth
          DO k=1,N(ng)
            DO i=Istr,Iend
!>            tl_v_north(i,k)=fac1*tl_v_obc(j,k,ib,it1,Linp)+           &
!>   &                        fac2*tl_v_obc(j,k,ib,it2,Linp)
!>
              ad_v_obc(i,k,ib,it1,Linp)=ad_v_obc(i,k,ib,it1,Linp)+      &
     &                                  fac1*ad_v_north(i,k)
              ad_v_obc(i,k,ib,it2,Linp)=ad_v_obc(i,k,ib,it2,Linp)+      &
     &                                  fac2*ad_v_north(i,k)
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
            ib=iwest
            DO k=1,N(ng)
              DO j=Jstr,Jend
!>              tl_t_west(j,k,it)=fac1*tl_t_obc(j,k,ib,it1,Linp,it)+    &
!>   &                            fac2*tl_t_obc(j,k,ib,it2,Linp,it)
!>
                ad_t_obc(j,k,ib,it1,Linp,it)=                           &
     &                                   ad_t_obc(j,k,ib,it1,Linp,it)+  &
     &                                   fac1*ad_t_west(j,k,it)
                ad_t_obc(j,k,ib,it2,Linp,it)=                           &
     &                                   ad_t_obc(j,k,ib,it2,Linp,it)+  &
     &                                   fac2*ad_t_west(j,k,it)
              END DO
            END DO
          END IF
#  endif
#  ifdef EAST_TOBC
          IF ((Lobc(ieast,isTvar(it),ng)).and.EASTERN_EDGE) THEN
            ib=ieast
            DO k=1,N(ng)
              DO j=Jstr,Jend
!>              tl_t_east(j,k,it)=fac1*tl_t_obc(j,k,ib,it1,Linp,it)+    &
!>   &                            fac2*tl_t_obc(j,k,ib,it2,Linp,it)
!>
                ad_t_obc(j,k,ib,it1,Linp,it)=                           &
     &                                   ad_t_obc(j,k,ib,it1,Linp,it)+  &
     &                                   fac1*ad_t_east(j,k,it)
                ad_t_obc(j,k,ib,it2,Linp,it)=                           &
     &                                   ad_t_obc(j,k,ib,it2,Linp,it)+  &
     &                                   fac2*ad_t_east(j,k,it)
              END DO
            END DO
          END IF
#  endif
#  ifdef SOUTH_TOBC
          IF ((Lobc(isouth,isTvar(it),ng)).and.SOUTHERN_EDGE) THEN
            ib=isouth
            DO k=1,N(ng)
              DO i=Istr,Iend
!>              tl_t_south(i,k,it)=fac1*tl_t_obc(i,k,ib,it1,Linp,it)+   &
!>   &                             fac2*tl_t_obc(i,k,ib,it2,Linp,it)
!>
                ad_t_obc(i,k,ib,it1,Linp,it)=                           &
     &                                   ad_t_obc(i,k,ib,it1,Linp,it)+  &
     &                                   fac1*ad_t_south(i,k,it)
                ad_t_obc(i,k,ib,it2,Linp,it)=                           &
     &                                   ad_t_obc(j,k,ib,it2,Linp,it)+  &
     &                                   fac2*ad_t_south(i,k,it)
              END DO
            END DO
          END IF
#  endif
#  ifdef NORTH_TOBC
          IF ((Lobc(inorth,isTvar(it),ng)).and.NORTHERN_EDGE) THEN
            ib=inorth
            DO k=1,N(ng)
              DO i=Istr,Iend
!>              tl_t_north(i,k,it)=fac1*tl_t_obc(i,k,ib,it1,Linp,it)+   &
!>   &                             fac2*tl_t_obc(i,k,ib,it2,Linp,it)
!>
                ad_t_obc(i,k,ib,it1,Linp,it)=                           &
     &                                   ad_t_obc(i,k,ib,it1,Linp,it)+  &
     &                                   fac1*ad_t_north(i,k,it)
                ad_t_obc(i,k,ib,it2,Linp,it)=                           &
     &                                   ad_t_obc(j,k,ib,it2,Linp,it)+  &
     &                                   fac2*ad_t_north(i,k,it)
              END DO
            END DO
          END IF
#  endif
        END IF
      END DO
# endif

      RETURN
      END SUBROUTINE ad_obc_adjust_tile
#endif
      END MODULE ad_obc_adjust_mod
