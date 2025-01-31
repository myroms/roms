#include "cppdefs.h"
      MODULE frc_iau_mod

#if defined RBL4DVAR && defined RPCG
!
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2019 The ROMS Group       Andrew M. Moore        !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine is used to compute the IAU forcing for the NLM.        !
!                                                                      !
!  The method is described in:                                         !
!  Bloom et al., 1996: Data Assimilation Using Incremental Analysis    !
!                      Update. Mon. Wea. Rev., 124, 1256-1271.         !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_forces
      USE mod_fourdvar
      USE mod_ocean
      USE mod_scalars
      USE mod_stepping


!
      implicit none
!
      PRIVATE
      PUBLIC :: frc_iau
      PUBLIC :: frc_iau_ini
!
      CONTAINS
!
      SUBROUTINE frc_iau_ini (ng, tile, irec)
!
!=======================================================================
!                                                                      !
!  This subroutine computes the incremental analysis update (IAU)      !
!  that is applied as a forcing term in the non-linear model during    !
!  the first timeIAU seconds of the assimilation window.               !
!                                                                      !
!  The adjoint arrays ad_var(irec) are used as temporary storage for   !
!  the IAU.                                                            !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng        Nested grid number.                                    !
!     tile      Domain partition.                                      !
!     irec      ad_var record containing the fields on input.          !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, irec
!
!  Local variable declarations.
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iADM, 7, __LINE__, __FILE__)
# endif
      CALL frc_iau_ini_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       irec,                                      &
     &                       OCEAN(ng) % ad_zeta,                       &
# ifdef SOLVE3D
     &                       OCEAN(ng) % ad_u,                          &
     &                       OCEAN(ng) % ad_v,                          &
     &                       OCEAN(ng) % ad_t,                          &
# else
     &                       OCEAN(ng) % ad_ubar,                       &
     &                       OCEAN(ng) % ad_vbar,                       &
# endif
     &                       OCEAN(ng) % f_zeta,                        &
# ifdef SOLVE3D
     &                       OCEAN(ng) % f_u,                           &
     &                       OCEAN(ng) % f_v,                           &
     &                       OCEAN(ng) % f_t)
# else
     &                       OCEAN(ng) % f_ubar,                        &
     &                       OCEAN(ng) % f_vbar)
# endif
# ifdef PROFILE
      CALL wclock_off (ng, iADM, 7, __LINE__, __FILE__)
# endif
!
      RETURN
      END SUBROUTINE frc_iau_ini
!
!***********************************************************************
      SUBROUTINE frc_iau_ini_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             IminS, ImaxS, JminS, JmaxS,          &
     &                             irec,                                &
     &                             ad_zeta,                             &
# ifdef SOLVE3D
     &                             ad_u, ad_v, ad_t,                    &
# else
     &                             ad_ubar, ad_vbar,                    &
# endif
     &                             f_zeta,                              &
# ifdef SOLVE3D
     &                             f_u, f_v, f_t)
# else
     &                             f_ubar, f_vbar)
# endif
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, irec
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: ad_zeta(LBi:,LBj:,:)
      real(r8), intent(inout) :: f_zeta(LBi:,LBj:)
#  ifdef SOLVE3D
      real(r8), intent(inout) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_v(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: f_u(LBi:,LBj:,:)
      real(r8), intent(inout) :: f_v(LBi:,LBj:,:)
      real(r8), intent(inout) :: f_t(LBi:,LBj:,:,:)
#  else
      real(r8), intent(inout) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_vbar(LBi:,LBj:,:)
      real(r8), intent(inout) :: f_ubar(LBi:,LBj:)
      real(r8), intent(inout) :: f_vbar(LBi:,LBj:)
#  endif
# else
      real(r8), intent(inout) :: ad_zeta(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: f_zeta(LBi:UBi,LBj:UBj)
#  ifdef SOLVE3D
      real(r8), intent(inout) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_t(LBi:UBi,LBj:UBj,N(ng),2,NT(ng))
      real(r8), intent(inout) :: f_u(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: f_v(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: f_t(LBi:UBi,LBj:UBj,N(ng),NT(ng))
#  else
      real(r8), intent(inout) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: f_ubar(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: f_vbar(LBi:UBi,LBj:UBj)
#  endif
# endif
!
!  Local variable declarations.
!
      integer :: i, it1, it2, j, k
# ifdef SOLVE3D
      integer :: itrc
# endif
      real(r8) :: fac1

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute IAU forcing terms used in the nonlinear model.
!  Use ad_var arrays as temporary storage for incremental analysis
!  update forcing - needed when weak constraing is also activated
!  and used in in the routine frc_iau.
!-----------------------------------------------------------------------
!
!  Set uniform weights in time.
!
      fac1=REAL(timeIAU(ng)/dt(ng),r8)
      fac1=1.0_r8/fac1
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
         f_zeta(i,j)=fac1*ad_zeta(i,j,irec)
         ad_zeta(i,j,irec)=f_zeta(i,j)
        END DO
      END DO

# ifndef SOLVE3D
!
!  Compute 2D-momentum iau forcing terms.
!
      DO j=JstrR,JendR
        DO i=Istr,IendR
           f_ubar(i,j)=fac1*ad_ubar(i,j,irec)
           ad_ubar(i,j,irec)=f_ubar(i,j)
        END DO
      END DO
      DO j=Jstr,JendR
        DO i=IstrR,IendR
           f_vbar(i,j)=fac1*ad_vbar(i,j,irec)
           ad_vbar(i,j,irec)=f_vbar(i,j)
        END DO
      END DO
# endif
# ifdef SOLVE3D
!
!  Compute 3D-momentum iau forcing terms.
!
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=Istr,IendR
           f_u(i,j,k)=fac1*ad_u(i,j,k,irec)
           ad_u(i,j,k,irec)=f_u(i,j,k)
          END DO
        END DO
      END DO
      DO k=1,N(ng)
        DO j=Jstr,JendR
          DO i=IstrR,IendR
           f_v(i,j,k)=fac1*ad_v(i,j,k,irec)
           ad_v(i,j,k,irec)=f_v(i,j,k)
          END DO
        END DO
      END DO
!
!  Compute tracer iau forcing terms.
!
      DO itrc=1,NT(ng)
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
             f_t(i,j,k,itrc)=fac1*ad_t(i,j,k,irec,itrc)
             ad_t(i,j,k,irec,itrc)=f_t(i,j,k,itrc)
            END DO
          END DO
        END DO
      END DO
# endif
!
      RETURN
      END SUBROUTINE frc_iau_ini_tile
!
      SUBROUTINE frc_iau (ng, tile, irec)
!
!=======================================================================
!                                                                      !
!  This subroutine computes the combined weak constraint forcing       !
!  and the incremental analysis update (IAU) when both weak constraint !
!  data assimilation and the IAU are activated.                        !
!  The input adjoint arrays ad_var(irec) hold the IAU forcing.         !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng        Nested grid number.                                    !
!     tile      Domain partition.                                      !
!     irec      ad_var record containing the fields on input.          !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, irec
!
!  Local variable declarations.
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iADM, 7, __LINE__, __FILE__)
# endif
      CALL frc_iau_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   irec,                                          &
     &                   OCEAN(ng) % ad_zeta,                           &
# ifdef SOLVE3D
     &                   OCEAN(ng) % ad_u,                              &
     &                   OCEAN(ng) % ad_v,                              &
     &                   OCEAN(ng) % ad_t,                              &
# else
     &                   OCEAN(ng) % ad_ubar,                           &
     &                   OCEAN(ng) % ad_vbar,                           &
# endif
     &                   OCEAN(ng) % f_zeta,                            &
# ifdef SOLVE3D
     &                   OCEAN(ng) % f_u,                               &
     &                   OCEAN(ng) % f_v,                               &
     &                   OCEAN(ng) % f_t)
# else
     &                   OCEAN(ng) % f_ubar,                            &
     &                   OCEAN(ng) % f_vbar)
# endif
# ifdef PROFILE
      CALL wclock_off (ng, iADM, 7, __LINE__, __FILE__)
# endif
!
      RETURN
      END SUBROUTINE frc_iau
!
!***********************************************************************
      SUBROUTINE frc_iau_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         irec,                                    &
     &                         ad_zeta,                                 &
# ifdef SOLVE3D
     &                         ad_u, ad_v, ad_t,                        &
# else
     &                         ad_ubar, ad_vbar,                        &
# endif
     &                         f_zeta,                                  &
# ifdef SOLVE3D
     &                         f_u, f_v, f_t)
# else
     &                         f_ubar, f_vbar)
# endif
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, irec
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: ad_zeta(LBi:,LBj:,:)
      real(r8), intent(inout) :: f_zeta(LBi:,LBj:)
#  ifdef SOLVE3D
      real(r8), intent(in) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: ad_v(LBi:,LBj:,:,:)
      real(r8), intent(in) :: ad_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: f_u(LBi:,LBj:,:)
      real(r8), intent(inout) :: f_v(LBi:,LBj:,:)
      real(r8), intent(inout) :: f_t(LBi:,LBj:,:,:)
#  else
      real(r8), intent(in) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: ad_vbar(LBi:,LBj:,:)
      real(r8), intent(inout) :: f_ubar(LBi:,LBj:)
      real(r8), intent(inout) :: f_vbar(LBi:,LBj:)
#  endif
# else
      real(r8), intent(in) :: ad_zeta(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: f_zeta(LBi:UBi,LBj:UBj)
#  ifdef SOLVE3D
      real(r8), intent(in) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: ad_t(LBi:UBi,LBj:UBj,N(ng),2,NT(ng))
      real(r8), intent(inout) :: f_u(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: f_v(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: f_t(LBi:UBi,LBj:UBj,N(ng),NT(ng))
#  else
      real(r8), intent(in) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: ad_vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: f_ubar(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: f_vbar(LBi:UBi,LBj:UBj)
#  endif
# endif
!
!  Local variable declarations.
!
      integer :: i, it1, it2, j, k
# ifdef SOLVE3D
      integer :: itrc
# endif

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Update the weak constraint forcing arrays with the IAU forcing
!  which is stored in the ad_var(irec) arrays.
!-----------------------------------------------------------------------
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
         f_zeta(i,j)=f_zeta(i,j)+ad_zeta(i,j,irec)
        END DO
      END DO

# ifndef SOLVE3D
!
!  Add 2D-momentum iau forcing terms.
!
      DO j=JstrR,JendR
        DO i=Istr,IendR
           f_ubar(i,j)=f_ubar(i,j)+ad_ubar(i,j,irec)
        END DO
      END DO
      DO j=Jstr,JendR
        DO i=IstrR,IendR
           f_vbar(i,j)=f_vbar(i,j)+ad_vbar(i,j,irec)
        END DO
      END DO
# endif
# ifdef SOLVE3D
!
!  Add 3D-momentum iau forcing terms.
!
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=Istr,IendR
           f_u(i,j,k)=f_u(i,j,k)+ad_u(i,j,k,irec)
          END DO
        END DO
      END DO
      DO k=1,N(ng)
        DO j=Jstr,JendR
          DO i=IstrR,IendR
           f_v(i,j,k)=f_v(i,j,k)+ad_v(i,j,k,irec)
          END DO
        END DO
      END DO
!
!  Add tracer iau forcing terms.
!
      DO itrc=1,NT(ng)
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
             f_t(i,j,k,itrc)=f_t(i,j,k,itrc)+ad_t(i,j,k,irec,itrc)
            END DO
          END DO
        END DO
      END DO
# endif
!
      RETURN
      END SUBROUTINE frc_iau_tile
#endif
      END MODULE frc_iau_mod
