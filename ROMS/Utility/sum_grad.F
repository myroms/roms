#include "cppdefs.h"
      MODULE sum_grad_mod
#if defined IS4DVAR
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes the sum of the background cost function       !
!  gradients in v-space                                                !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC sum_grad

      CONTAINS
!
!***********************************************************************
      SUBROUTINE sum_grad (ng, tile, Linp, Lout)
!***********************************************************************
!
      USE mod_param
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      USE mod_forces
# endif
      USE mod_ocean
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, Linp, Lout
!
!  Local variable declarations.
!
# include "tile.h"
!
      CALL sum_grad_tile (ng, tile,                                     &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     Linp, Lout,                                  &
# ifdef ADJUST_WSTRESS
    &                      FORCES(ng) % tl_ustr,                        &
    &                      FORCES(ng) % tl_vstr,                        &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
    &                      FORCES(ng) % tl_tflux,                       &
#  endif
     &                     OCEAN(ng) % tl_t,                            &
     &                     OCEAN(ng) % tl_u,                            &
     &                     OCEAN(ng) % tl_v,                            &
# else
     &                     OCEAN(ng) % tl_ubar,                         &
     &                     OCEAN(ng) % tl_vbar,                         &
# endif
     &                     OCEAN(ng) % tl_zeta)
      RETURN
      END SUBROUTINE sum_grad
!
!***********************************************************************
      SUBROUTINE sum_grad_tile (ng, tile,                               &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           Linp, Lout,                            &
# ifdef ADJUST_WSTRESS
    &                            tl_ustr, tl_vstr,                      &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
    &                            tl_tflux,                              &
#  endif
     &                           tl_t, tl_u, tl_v,                      &
# else
     &                           tl_ubar, tl_vbar,                      &
# endif
     &                           tl_zeta)
!***********************************************************************
!
      USE mod_param
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      USE mod_scalars
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Linp, Lout
!
# ifdef ASSUMED_SHAPE
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: tl_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_vstr(LBi:,LBj:,:,:)
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: tl_tflux(LBi:,LBj:,:,:,:)
#   endif
      real(r8), intent(inout) :: tl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_v(LBi:,LBj:,:,:)
#  else
      real(r8), intent(inout) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_vbar(LBi:,LBj:,:)
#  endif
      real(r8), intent(inout) :: tl_zeta(LBi:,LBj:,:)
# else
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: tl_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: tl_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
#  endif
#  ifdef SOLVE3D 
#   ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: tl_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
#   endif
      real(r8), intent(inout) :: tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: tl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: tl_v(LBi:UBi,LBj:UBj,N(ng),2)
#  else
      real(r8), intent(inout) :: tl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: tl_vbar(LBi:UBi,LBj:UBj,3)
#  endif
      real(r8), intent(inout) :: tl_zeta(LBi:UBi,LBj:UBj,3)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k
# ifdef SOLVE3D
      integer :: itrc
# endif

# include "set_bounds.h"
!
!  Free-surface gradient norm.
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          tl_zeta(i,j,Lout)=tl_zeta(i,j,Linp)+                          &
     &                      tl_zeta(i,j,Lout)
        END DO
      END DO

# ifdef ADJUST_WSTRESS
!
!  Surface momentum stress gradient norm.
!
      DO k=1,Nfrec(ng)
        DO j=JstrR,JendR
          DO i=Istr,IendR
            tl_ustr(i,j,k,Lout)=tl_ustr(i,j,k,Linp)+                    &
     &                          tl_ustr(i,j,k,Lout)
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            tl_vstr(i,j,k,Lout)=tl_vstr(i,j,k,Linp)+                    &
     &                          tl_vstr(i,j,k,Lout)
          END DO
        END DO
      END DO
# endif
# if !defined SOLVE3D
!
!  2D momentum gradient norm.
!
      DO j=JstrR,JendR
        DO i=Istr,IendR
          tl_ubar(i,j,Lout)=tl_ubar(i,j,Linp)+                          &
     &                      tl_ubar(i,j,Lout)
        END DO
      END DO
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          tl_vbar(i,j,Lout)=tl_vbar(i,j,Linp)+                          &
     &                      tl_vbar(i,j,Lout)
        END DO
      END DO
# else
!
!  3D momentum gradient norm.
!
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=Istr,IendR
            tl_u(i,j,k,Lout)=tl_u(i,j,k,Linp)+                          &
     &                       tl_u(i,j,k,Lout)
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            tl_v(i,j,k,Lout)=tl_v(i,j,k,Linp)+                          &
     &                       tl_v(i,j,k,Lout)
          END DO
        END DO
      END DO
!
!  Tracers gradient norm.
!
      DO itrc=1,NT(ng)
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              tl_t(i,j,k,Lout,itrc)=tl_t(i,j,k,Linp,itrc)+              &
     &                              tl_t(i,j,k,Lout,itrc)
            END DO
          END DO
        END DO
      END DO

#  ifdef ADJUST_STFLUX
!
!  Surface tracers flux gradient norm.
!
      DO itrc=1,NT(ng)
        DO k=1,Nfrec(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              tl_tflux(i,j,k,Lout,itrc)=tl_tflux(i,j,k,Linp,itrc)+      &
     &                                  tl_tflux(i,j,k,Lout,itrc)
            END DO
          END DO
        END DO
      END DO
#  endif
# endif

      RETURN
      END SUBROUTINE sum_grad_tile
#endif
      END MODULE sum_grad_mod
