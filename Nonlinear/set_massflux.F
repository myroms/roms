#include "cppdefs.h"
      MODULE set_massflux_mod

#ifdef SOLVE3D
# ifdef EW_PERIODIC
#  define IU_RANGE Istr,Iend
#  define IV_RANGE Istr,Iend
# else
#  define IU_RANGE Istr,IendR
#  define IV_RANGE IstrR,IendR
# endif

# ifdef NS_PERIODIC
#  define JU_RANGE Jstr,Jend
#  define JV_RANGE Jstr,Jend
# else
#  define JU_RANGE JstrR,JendR
#  define JV_RANGE Jstr,JendR
# endif
!
!========================================== Alexander F. Shchepetkin ===
!  Copyright (c) 2005 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine computes horizontal mass fluxes, Hz*u/n and Hz*v/m.    !
!                                                                      !
!=======================================================================
!
!
      implicit none

      PRIVATE
      PUBLIC  :: set_massflux
# ifdef ADJOINT
      PUBLIC  :: reset_massflux
# endif

      CONTAINS
!
!***********************************************************************
      SUBROUTINE set_massflux (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_ocean
      USE mod_stepping
!
      integer, intent(in) :: ng, tile

# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 12)
# endif
      CALL set_massflux_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        nrhs(ng),                                 &
     &                        OCEAN(ng) % u,                            &
     &                        OCEAN(ng) % v,                            &
     &                        GRID(ng) % Hz,                            &
     &                        GRID(ng) % om_v,                          &
     &                        GRID(ng) % on_u,                          &
     &                        GRID(ng) % Huon,                          &
     &                        GRID(ng) % Hvom)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 12)
# endif
      RETURN
      END SUBROUTINE set_massflux
!
!***********************************************************************
      SUBROUTINE set_massflux_tile (ng, Istr, Iend, Jstr, Jend,         &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              nrhs,                               &
     &                              u, v, Hz, om_v, on_u,               &
     &                              Huon, Hvom)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars

# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
!
      USE exchange_3d_mod
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: nrhs

# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)

      real(r8), intent(out) :: Huon(LBi:,LBj:,:)
      real(r8), intent(out) :: Hvom(LBi:,LBj:,:)
# else
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: om_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)

      real(r8), intent(out) :: Huon(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(out) :: Hvom(LBi:UBi,LBj:UBj,N(ng))
# endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j, k

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute horizontal mass fluxes, Hz*u/n and Hz*v/m.
!-----------------------------------------------------------------------
!
      DO k=1,N(ng)
        DO j=JU_RANGE
          DO i=IU_RANGE
            Huon(i,j,k)=0.5_r8*(Hz(i,j,k)+Hz(i-1,j,k))*u(i,j,k,nrhs)*   &
     &                  on_u(i,j)
          END DO
        END DO
        DO j=JV_RANGE
          DO i=IV_RANGE
            Hvom(i,j,k)=0.5_r8*(Hz(i,j,k)+Hz(i,j-1,k))*v(i,j,k,nrhs)*   &
     &                  om_v(i,j)
          END DO
        END DO
      END DO
# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
!
!  Exchange boundary information.
!
      CALL exchange_u3d_tile (ng, iNLM, Istr, Iend, Jstr, Jend,         &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        NghostPoints,                             &
     &                        Huon)
      CALL exchange_v3d_tile (ng, iNLM, Istr, Iend, Jstr, Jend,         &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        NghostPoints,                             &
     &                        Hvom)
# endif
      RETURN
      END SUBROUTINE set_massflux_tile

# ifdef ADJOINT
!
!***********************************************************************
      SUBROUTINE reset_massflux (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_coupling
      USE mod_grid
      USE mod_ocean
      USE mod_stepping
!
      integer, intent(in) :: ng, tile

#  include "tile.h"
!
#  ifdef PROFILE
      CALL wclock_on (ng, iADM, 12)
#  endif
      CALL reset_massflux_tile (ng, Istr, Iend, Jstr, Jend,             &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          nnew(ng),                               &
     &                          COUPLING(ng) % DU_avg2,                 &
     &                          COUPLING(ng) % DV_avg2,                 &
     &                          OCEAN(ng) % u,                          &
     &                          OCEAN(ng) % v,                          &
     &                          GRID(ng) % Hz,                          &
     &                          GRID(ng) % om_v,                        &
     &                          GRID(ng) % on_u,                        &
     &                          GRID(ng) % Huon,                        &
     &                          GRID(ng) % Hvom)
#  ifdef PROFILE
      CALL wclock_off (ng, iADM, 12)
#  endif
      RETURN
      END SUBROUTINE reset_massflux
!
!***********************************************************************
      SUBROUTINE reset_massflux_tile (ng, Istr, Iend, Jstr, Jend,       &
     &                                LBi, UBi, LBj, UBj,               &
     &                                nnew,                             &
     &                                DU_avg2, DV_avg2,                 &
     &                                u, v, Hz, om_v, on_u,             &
     &                                Huon, Hvom)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars

#  if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
!
      USE exchange_3d_mod
#  endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: nnew

#  ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: DU_avg2(LBi:,LBj:)
      real(r8), intent(in) :: DV_avg2(LBi:,LBj:)
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)

      real(r8), intent(inout) :: Huon(LBi:,LBj:,:)
      real(r8), intent(inout) :: Hvom(LBi:,LBj:,:)
#  else
      real(r8), intent(in) :: DU_avg2(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: DV_avg2(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: om_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)

      real(r8), intent(inout) :: Huon(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: Hvom(LBi:UBi,LBj:UBj,N(ng))
#  endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j, k

      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N(ng)) :: DC
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N(ng)) :: FC


#  include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute intermediate values of mass fluxes Huon and Hvom used by the
!  adjoint model.  The original values can be reinstated by calling
!  "set_massflux" after "ad_omega".
!-----------------------------------------------------------------------

#  ifdef EW_PERIODIC
#   define IU_RANGE Istr,Iend
#   define IV_RANGE Istr,Iend
#  else
#   define IU_RANGE Istr,IendR
#   define IV_RANGE IstrR,IendR
#  endif
#  ifdef NS_PERIODIC
#   define JU_RANGE Jstr,Jend
#   define JV_RANGE Jstr,Jend
#  else
#   define JU_RANGE JstrR,JendR
#   define JV_RANGE Jstr,JendR
#  endif
!
!  Compute mass flux, Hz*u/n.
!
      DO j=JU_RANGE
        DO i=IU_RANGE
          DC(i,0)=0.0_r8
          FC(i,0)=0.0_r8
        END DO
        DO k=1,N(ng)
          DO i=IU_RANGE
            DC(i,k)=0.5_r8*(Hz(i,j,k)+Hz(i-1,j,k))*on_u(i,j)
            DC(i,0)=DC(i,0)+DC(i,k)
          END DO
        END DO
        DO k=N(ng),1,-1
          DO i=IU_RANGE
            Huon(i,j,k)=0.5_r8*(Huon(i,j,k)+u(i,j,k,nnew)*DC(i,k))
            FC(i,0)=FC(i,0)+Huon(i,j,k)
          END DO
        END DO
!
!  Replace with correct vertical mean, DU_avg2.
!
        DO i=IU_RANGE
          DC(i,0)=1.0_r8/DC(i,0)
          FC(i,0)=DC(i,0)*(FC(i,0)-DU_avg2(i,j))
        END DO
        DO k=1,N(ng)
          DO i=IU_RANGE
            Huon(i,j,k)=Huon(i,j,k)-DC(i,k)*FC(i,0)
          END DO
        END DO
!
!  Compute mass flux, Hz*v/m.
!
        IF (j.ge.Jstr) THEN
          DO i=IV_RANGE
            DC(i,0)=0.0_r8
            FC(i,0)=0.0_r8
          END DO
          DO k=1,N(ng)
            DO i=IV_RANGE
              DC(i,k)=0.5_r8*(Hz(i,j,k)+Hz(i,j-1,k))*om_v(i,j)
              DC(i,0)=DC(i,0)+DC(i,k)
            END DO
          END DO
          DO k=N(ng),1,-1
            DO i=IV_RANGE
              Hvom(i,j,k)=0.5_r8*(Hvom(i,j,k)+v(i,j,k,nnew)*DC(i,k))
              FC(i,0)=FC(i,0)+Hvom(i,j,k)
            END DO
          END DO
!
!  Replace with correct vertical mean, DV_avg2.
!
          DO i=IV_RANGE
            DC(i,0)=1.0_r8/DC(i,0)
            FC(i,0)=DC(i,0)*(FC(i,0)-DV_avg2(i,j))
          END DO
          DO k=1,N(ng)
            DO i=IV_RANGE
              Hvom(i,j,k)=Hvom(i,j,k)-DC(i,k)*FC(i,0)
            END DO
          END DO
        ENDIF
      END DO
#  undef IU_RANGE
#  undef JU_RANGE
#  undef IV_RANGE
#  undef JV_RANGE
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
!
!  Exchange boundary information.
!
      CALL exchange_u3d_tile (ng, iADM, Istr, Iend, Jstr, Jend,         &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        NghostPoints,                             &
     &                        Huon)
      CALL exchange_v3d_tile (ng, iADM, Istr, Iend, Jstr, Jend,         &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        NghostPoints,                             &
     &                        Hvom)
#  endif
      RETURN
      END SUBROUTINE reset_massflux_tile
# endif
#endif
      END MODULE set_massflux_mod
