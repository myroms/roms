#include "cppdefs.h"
      MODULE sediment_mod
#ifdef SEDIMENT
!
!==================================================== John C. Warner ===
!  Copyright (c) 2002 ROMS/TOMS Group       Alexander F. Shchepetkin   !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This  routine computes the sediment sources and sinks and adds      !
!  then the global sediment tracer fields. Currently, it includes      !
!  the following:                                                      !
!                                                                      !
!  * Vertical settling of sediment in the water column.                !
!  * Erosive and depositional flux interactions of sediment            !
!    between water column and the bed.                                 !
!  * Transport of multiple grain sizes.                                !
!  * Bed layer stratigraphy.
!  * Seawater/sediment vertical level distribution:                    !
!                                                                      !
!         W-level  RHO-level                                           !
!                                                                      !
!            N     _________                                           !
!                 |         |                                          !
!                 |    N    |                                          !
!          N-1    |_________|  S                                       !
!                 |         |  E                                       !
!                 |   N-1   |  A                                       !
!            2    |_________|  W                                       !
!                 |         |  A                                       !
!                 |    2    |  T                                       !
!            1    |_________|  E                                       !
!                 |         |  R                                       !
!                 |    1    |                                          !
!            0    |_________|_____ bathymetry                          !
!                 |/////////|                                          !
!                 |    1    |                                          !
!            1    |_________|  S                                       !
!                 |         |  E                                       !
!                 |    2    |  D                                       !
!            2    |_________|  I                                       !
!                 |         |  M                                       !
!                 |  Nbed-1 |  E                                       !
!        Nbed-1   |_________|  N                                       !
!                 |         |  T                                       !
!                 |  Nbed   |                                          !
!         Nbed    |_________|                                          !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: sediment

      CONTAINS
!
!***********************************************************************
      SUBROUTINE sediment (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_forces
      USE mod_grid
      USE mod_ocean
      USE mod_stepping
#ifdef BBL_MODEL
      USE mod_bbl
#endif
!
      integer, intent(in) :: ng, tile

# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, 16)
# endif
      CALL sediment_tile (ng, Istr, Iend, Jstr, Jend,                   &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    nnew(ng),                                     &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
# ifdef MASKING
     &                    GRID(ng) % rmask,                             &
# endif
     &                    GRID(ng) % Hz,                                &
     &                    GRID(ng) % z_w,                               &
# ifdef BBL_MODEL
     &                    BBL(ng) % bustrcwmax,                         &
     &                    BBL(ng) % bvstrcwmax,                         &
# else
     &                    FORCES(ng) % bustr,                           &
     &                    FORCES(ng) % bvstr,                           &
# endif
     &                    OCEAN(ng) % t,                                &
     &                    OCEAN(ng) % bed,                              &
     &                    OCEAN(ng) % bed_frac,                         &
     &                    OCEAN(ng) % bed_mass,                         &
     &                    OCEAN(ng) % bottom)
# ifdef PROFILE
      CALL wclock_off (ng, 16)
# endif
      RETURN
      END SUBROUTINE sediment
!
!***********************************************************************
      SUBROUTINE sediment_tile (ng, Istr, Iend, Jstr, Jend,             &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          nnew,                                   &
     &                          pm, pn,                                 &
# ifdef MASKING
     &                          rmask,                                  &
# endif
     &                          Hz, z_w,                                &
# ifdef BBL_MODEL
     &                          bustrcwmax, bvstrcwmax,                 &
# else
     &                          bustr, bvstr,                           &
# endif
     &                          t, bed, bed_frac, bed_mass, bottom)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_sediment
!
      USE bc_3d_mod, ONLY : bc_r3d_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: nnew

# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
#  ifdef BBL_MODEL
      real(r8), intent(in) :: bustrcwmax(LBi:,LBj:)
      real(r8), intent(in) :: bvstrcwmax(LBi:,LBj:)
#  else
      real(r8), intent(in) :: bustr(LBi:,LBj:)
      real(r8), intent(in) :: bvstr(LBi:,LBj:)
#  endif
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: bed(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: bed_frac(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: bed_mass(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: bottom(LBi:,LBj:,:)
# else
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
#  ifdef BBL_MODEL
      real(r8), intent(in) :: bustrcwmax(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: bvstrcwmax(LBi:UBi,LBj:UBj)
#  else
      real(r8), intent(in) :: bustr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: bvstr(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: bed(LBi:UBi,LBj:UBj,Nbed,MBEDP)
      real(r8), intent(inout) :: bed_frac(LBi:UBi,LBj:UBj,Nbed,NST)
      real(r8), intent(inout) :: bed_mass(LBi:UBi,LBj:UBj,Nbed,NST)
      real(r8), intent(inout) :: bottom(LBi:UBi,LBj:UBj,MBOTP)
# endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, indx, ised, j, k, ks

      real(r8) :: cff, cff1, cff2, cff3, cffL, cffR, dltL, dltR
      real(r8) :: cu, ero_flux, tau_water, cff4

      integer, dimension(PRIVATE_1D_SCRATCH_ARRAY,N(ng)) :: ksource

      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N(ng)) :: FC

      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,N(ng)) :: Hz_inv
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,N(ng)) :: Hz_inv2
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,N(ng)) :: Hz_inv3
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,N(ng)) :: qc
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,N(ng)) :: qR
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,N(ng)) :: qL
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,N(ng)) :: WR
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,N(ng)) :: WL

      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,NST) :: dep_mass

      integer :: Ksed, ised2
      real(r8) :: thck_avail, thck_to_add
      real(r8), parameter :: eps = 1.0E-10_r8

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Update bed layer stratigraphy.
!-----------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO i=Istr,Iend 
!
!  Calculate bottom stress.
!
# ifdef BBL_MODEL
          tau_water=0.5_r8*SQRT((bustrcwmax(i,j)+bustrcwmax(i+1,j))*    &
     &                          (bustrcwmax(i,j)+bustrcwmax(i+1,j))+    &
     &                          (bvstrcwmax(i,j)+bvstrcwmax(i,j+1))*    &
     &                          (bvstrcwmax(i,j)+bvstrcwmax(i,j+1)))
# else
          tau_water=0.5_r8*SQRT((bustr(i,j)+bustr(i+1,j))*              &
     &                          (bustr(i,j)+bustr(i+1,j))+              &
     &                          (bvstr(i,j)+bvstr(i,j+1))*              &
     &                          (bvstr(i,j)+bvstr(i,j+1)))
# endif
!
!  Calculate active layer thickness, bottom(i,j,iactv).
!
          bottom(i,j,iactv)=MAX(0.0_r8,                                 &
     &                          0.007_r8*                               &
     &                          (tau_water-bottom(i,j,itauc))*rho0)+    &
     &                          6.0_r8*bottom(i,j,isd50)
!
!  Ensure top bed layer thickness is greater or equal than the active
!  layer thickness. If need to add sediment to top layer, then entrain
!  from lower levels. Create new layers at bottom to maintain Nbed.
!
          IF (bottom(i,j,iactv).gt.bed(i,j,1,ithck)) THEN
            IF (Nbed.eq.1) THEN
              bottom(i,j,iactv)=bed(i,j,1,ithck)
            ELSE
              thck_to_add=bottom(i,j,iactv)-bed(i,j,1,ithck)
              thck_avail=0.0_r8
              Ksed=1                                         ! initialize
              DO k=2,Nbed
                IF (thck_avail.lt.thck_to_add) THEN
                  thck_avail=thck_avail+bed(i,j,k,ithck)
                  Ksed=k
                ENDIF
              END DO
!
!  Catch here if there was not enough bed material.
!
              IF (thck_avail.lt.thck_to_add) THEN
                bottom(i,j,iactv)=bed(i,j,1,ithck)+thck_avail
                thck_to_add=thck_avail
              END IF
!
!  Upate bed mass of top layer and fractional layer.
!
              cff2=(thck_avail-thck_to_add)/                            &
     &             MAX(bed(i,j,Ksed,ithck),eps)
              DO ised=1,NST
                cff1=0.0_r8 
                DO k=1,Ksed
                  cff1=cff1+bed_mass(i,j,k,ised)
                END DO
                bed_mass(i,j,1   ,ised)=cff1-                           &
     &                                  bed_mass(i,j,Ksed,ised)*cff2
                bed_mass(i,j,Ksed,ised)=bed_mass(i,j,Ksed,ised)*cff2
              END DO
!
!  Update thickness of fractional layer Ksed.
!
              bed(i,j,Ksed,ithck)=thck_avail-thck_to_add
!
!  Upate bed fraction of top layer.
!
              cff3=0.0_r8
              DO ised=1,NST
                cff3=cff3+bed_mass(i,j,1,ised)
              END DO
              DO ised=1,NST
                bed_frac(i,j,1,ised)=bed_mass(i,j,1,ised)/              &
     &                               MAX(cff3,eps)
              END DO
!
!  Update bed thickness of top layer.
!
              bed(i,j,1,ithck)=bottom(i,j,iactv)
!
!  Pull all layers closer to the surface.
!
              DO k=Ksed,Nbed
                ks=Ksed-2
                bed(i,j,k-ks,ithck)=bed(i,j,k,ithck)
                bed(i,j,k-ks,iporo)=bed(i,j,k,iporo)
                bed(i,j,k-ks,iaged)=bed(i,j,k,iaged)
                DO ised=1,NST
                  bed_frac(i,j,k-ks,ised)=bed_frac(i,j,k,ised)
                  bed_mass(i,j,k-ks,ised)=bed_mass(i,j,k,ised)
                ENDDO
              END DO
!
!  Add new layers onto the bottom. Split what was in the bottom layer to
!  fill these new empty cells. ("ks" is the number of new layers).
!
              ks=Ksed-2
              cff=1.0_r8/REAL(ks+1,r8)
              DO k=Nbed,Nbed-ks,-1
                bed(i,j,k,ithck)=bed(i,j,Nbed-ks,ithck)*cff
!!              bed(i,j,k,iporo)=bed(i,j,Nbed-ks,iporo)
                bed(i,j,k,iaged)=bed(i,j,Nbed-ks,iaged)
                DO ised=1,NST
                  bed_frac(i,j,k,ised)=bed_frac(i,j,Nbed-ks,ised)
                  bed_mass(i,j,k,ised)=MAX(eps,                         &
     &                                     bed_mass(i,j,Nbed-ks,ised)*  &
     &                                     cff)
                ENDDO
              END DO
            END IF                             ! Nbed > 1
          END IF                               ! increase top bed layer
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Add sediment Source/Sink terms.
!-----------------------------------------------------------------------
!
!  Compute inverse thicknessed to avoid repeated divisions.
!
      J_LOOP : DO j=Jstr,Jend
        DO k=1,N(ng)
          DO i=Istr,Iend
            Hz_inv(i,k)=1.0_r8/Hz(i,j,k)
          END DO
        END DO
        DO k=1,N(ng)-1
          DO i=Istr,Iend
            Hz_inv2(i,k)=1.0_r8/(Hz(i,j,k)+Hz(i,j,k+1))
          END DO
        END DO
        DO k=2,N(ng)-1
          DO i=Istr,Iend
            Hz_inv3(i,k)=1.0_r8/(Hz(i,j,k-1)+Hz(i,j,k)+Hz(i,j,k+1))
          END DO
        END DO
!
!  Copy concentration of suspended sediment into scratch array "qc"
!  (q-central, restrict it to be positive) which is hereafter
!  interpreted as a set of grid-box averaged values for sediment
!  concentration.
!
        SED_LOOP: DO ised=1,NST
          indx=idsed(ised)
          DO k=1,N(ng)
            DO i=Istr,Iend
              qc(i,k)=t(i,j,k,nnew,indx)*Hz_inv(i,k)
            END DO
          END DO
!
!-----------------------------------------------------------------------
!  Vertical sinking of suspended sediment.
!-----------------------------------------------------------------------
!
!  Reconstruct vertical profile of suspended sediment "qc" in terms
!  of a set of parabolic segments within each grid box. Then, compute
!  semi-Lagrangian flux due to sinking.
!
          DO k=N(ng)-1,1,-1
            DO i=Istr,Iend
              FC(i,k)=(qc(i,k+1)-qc(i,k))*Hz_inv2(i,k)
            END DO
          END DO
          DO k=2,N(ng)-1
            DO i=Istr,Iend
              dltR=Hz(i,j,k)*FC(i,k)
              dltL=Hz(i,j,k)*FC(i,k-1)
              cff=Hz(i,j,k-1)+2.0_r8*Hz(i,j,k)+Hz(i,j,k+1)
              cffR=cff*FC(i,k)
              cffL=cff*FC(i,k-1)
!
!  Apply PPM monotonicity constraint to prevent oscillations within the
!  grid box.
!
              IF ((dltR*dltL).le.0.0_r8) THEN
                dltR=0.0_r8
                dltL=0.0_r8
              ELSE IF (ABS(dltR).gt.ABS(cffL)) THEN
                dltR=cffL
              ELSE IF (ABS(dltL).gt.ABS(cffR)) THEN
                dltL=cffR
              END IF
!
!  Compute right and left side values (qR,qL) of parabolic segments
!  within grid box Hz(k); (WR,WL) are measures of quadratic variations. 
!
!  NOTE: Although each parabolic segment is monotonic within its grid
!        box, monotonicity of the whole profile is not guaranteed,
!        because qL(k+1)-qR(k) may still have different sign than
!        qc(k+1)-qc(k).  This possibility is excluded, after qL and qR
!        are reconciled using WENO procedure.
!
              cff=(dltR-dltL)*Hz_inv3(i,k)
              dltR=dltR-cff*Hz(i,j,k+1)
              dltL=dltL+cff*Hz(i,j,k-1)
              qR(i,k)=qc(i,k)+dltR
              qL(i,k)=qc(i,k)-dltL
              WR(i,k)=(2.0_r8*dltR-dltL)**2
              WL(i,k)=(dltR-2.0_r8*dltL)**2
            END DO
          END DO
          cff=1.0E-14_r8
          DO k=2,N(ng)-2
            DO i=Istr,Iend
              dltL=MAX(cff,WL(i,k  ))
              dltR=MAX(cff,WR(i,k+1))
              qR(i,k)=(dltR*qR(i,k)+dltL*qL(i,k+1))/(dltR+dltL)
              qL(i,k+1)=qR(i,k)
            END DO
          END DO
          DO i=Istr,Iend
            FC(i,N(ng))=0.0_r8              ! no-flux boundary condition
# if defined LINEAR_CONTINUATION
            qL(i,N(ng))=qR(i,N(ng)-1)
            qR(i,N(ng))=2.0_r8*qc(i,N(ng))-qL(i,N(ng))
# elif defined NEUMANN
            qL(i,N(ng))=qR(i,N(ng)-1)
            qR(i,N(ng))=1.5*qc(i,N(ng))-0.5_r8*qL(i,N(ng))
# else
            qR(i,N(ng))=qc(i,N(ng))         ! default strictly monotonic
            qL(i,N(ng))=qc(i,N(ng))         ! conditions
            qR(i,N(ng)-1)=qc(i,N(ng))
# endif
# if defined LINEAR_CONTINUATION 
            qR(i,1)=qL(i,2)
            qL(i,1)=2.0_r8*qc(i,1)-qR(i,1)
# elif defined NEUMANN
            qR(i,1)=qL(i,2)
            qL(i,1)=1.5_r8*qc(i,1)-0.5_r8*qR(i,1)
# else  
            qL(i,2)=qc(i,1)                 ! bottom grid boxes are
            qR(i,1)=qc(i,1)                 ! re-assumed to be
            qL(i,1)=qc(i,1)                 ! piecewise constant.
# endif
          END DO
!
!  Apply monotonicity constraint again, since the reconciled interfacial
!  values may cause a non-monotonic behavior of the parabolic segments
!  inside the grid box.
!
          DO k=1,N(ng)
            DO i=Istr,Iend
              dltR=qR(i,k)-qc(i,k)
              dltL=qc(i,k)-qL(i,k)
              cffR=2.0_r8*dltR
              cffL=2.0_r8*dltL
              IF ((dltR*dltL).lt.0.0_r8) THEN
                dltR=0.0_r8
                dltL=0.0_r8
              ELSE IF (ABS(dltR).gt.ABS(cffL)) THEN
                dltR=cffL
              ELSE IF (ABS(dltL).gt.ABS(cffR)) THEN
                dltL=cffR
              END IF
              qR(i,k)=qc(i,k)+dltR
              qL(i,k)=qc(i,k)-dltL
            END DO
          END DO
!
!  After this moment reconstruction is considered complete. The next
!  stage is to compute vertical advective fluxes, FC. It is expected
!  that sinking may occurs relatively fast, the algorithm is designed
!  to be free of CFL criterion, which is achieved by allowing
!  integration bounds for semi-Lagrangian advective flux to use as
!  many grid boxes in upstream direction as necessary.
!
!  In the two code segments below, WL is the z-coordinate of the
!  departure point for grid box interface z_w with the same indices;
!  FC is the finite volume flux; ksource(:,k) is index of vertical
!  grid box which contains the departure point (restricted by N(ng)). 
!  During the search: also add in content of whole grid boxes
!  participating in FC.
!
          cff=dt(ng)*ABS(Wsed(ised,ng))
          DO k=1,N(ng)
            DO i=Istr,Iend
              FC(i,k-1)=0.0_r8
              WL(i,k)=z_w(i,j,k-1)+cff
              WR(i,k)=Hz(i,j,k)*qc(i,k)
              ksource(i,k)=k
            END DO
          END DO
          DO k=1,N(ng)
            DO ks=k,N(ng)-1
              DO i=Istr,Iend
                IF (WL(i,k).gt.z_w(i,j,ks)) THEN
                  ksource(i,k)=ks+1
                  FC(i,k-1)=FC(i,k-1)+WR(i,ks)
                END IF
              END DO
            END DO
          END DO
!
!  Finalize computation of flux: add fractional part.
!
          DO k=1,N(ng)
            DO i=Istr,Iend
              ks=ksource(i,k)
              cu=MIN(1.0_r8,(WL(i,k)-z_w(i,j,ks-1))*Hz_inv(i,ks))
              FC(i,k-1)=FC(i,k-1)+                                      &
     &                  Hz(i,j,ks)*cu*                                  &
     &                  (qL(i,ks)+                                      &
     &                   cu*(0.5_r8*(qR(i,ks)-qL(i,ks))-                &
     &                       (1.5_r8-cu)*                               &
     &                       (qR(i,ks)+qL(i,ks)-2.0_r8*qc(i,ks))))
            END DO
          END DO
          DO k=1,N(ng)
            DO i=Istr,Iend
              qc(i,k)=qc(i,k)+(FC(i,k)-FC(i,k-1))*Hz_inv(i,k)
            END DO
          END DO
!
!-----------------------------------------------------------------------
!  Sediment deposition and resuspension near the bottom.
!-----------------------------------------------------------------------
!
!  The deposition and resuspension of sediment on the bottom "bed"
!  is due to precepitation flux FC(:,0), already computed, and the
!  resuspension (erosion, hence called ero_flux). The resuspension is
!  applied to the bottom-most grid box value qc(:,1) so the total mass
!  is conserved. Restrict "ero_flux" so that "bed" cannot go negative
!  after both fluxes are applied.
!
          cff=1.0_r8/tau_ce(ised,ng)
          DO i=Istr,Iend
            dep_mass(i,ised)=0.0_r8
# ifdef BBL_MODEL
          tau_water=0.5_r8*SQRT((bustrcwmax(i,j)+bustrcwmax(i+1,j))*    &
     &                          (bustrcwmax(i,j)+bustrcwmax(i+1,j))+    &
     &                          (bvstrcwmax(i,j)+bvstrcwmax(i,j+1))*    &
     &                          (bvstrcwmax(i,j)+bvstrcwmax(i,j+1)))
# else
            tau_water=0.5_r8*SQRT((bustr(i,j)+bustr(i+1,j))*            &
     &                            (bustr(i,j)+bustr(i+1,j))+            &
     &                            (bvstr(i,j)+bvstr(i,j+1))*            &
     &                            (bvstr(i,j)+bvstr(i,j+1)))
# endif
!
!  Compute erosion, ero_flux (kg/m2).
!
            cff1=(1.0_r8-bed(i,j,1,iporo))*bed_frac(i,j,1,ised)
            cff2=dt(ng)*Erate(ised,ng)*cff1
            cff3=Srho(ised,ng)*cff1
            ero_flux=MIN(MAX(0.0_r8,cff2*(cff*tau_water-1.0_r8)),       &
     &                   cff3*bottom(i,j,iactv)+FC(i,0))
!
!  Check if depositional.
!
            IF ((ero_flux-FC(i,0)).lt.0.0_r8) THEN
!
!  If first time step of deposit, then store deposit material in
!  temporary array, dep_mass.
!
              IF ((time(ng).gt.(bed(i,j,1,iaged)+1.1_r8*dt(ng))).and.   &
                  (bed(i,j,1,ithck).gt.0.005_r8))THEN
                dep_mass(i,ised)=-(ero_flux-FC(i,0))+eps
              ELSE
!
!  If it is not first time step of deposit, update bed thickness of
!  top layer.
!  
                bed(i,j,1,ithck)=bed(i,j,1,ithck)-                      &
     &                           (ero_flux-FC(i,0))/                    &
     &                           (Srho(ised,ng)*                        &
     &                            (1.0_r8-bed(i,j,1,iporo)))
# ifdef MASKING
                bed(i,j,1,ithck)=bed(i,j,1,ithck)*rmask(i,j)
# endif 
!
!  Upate bed mass of top layer.
!
                bed_mass(i,j,1,ised)=bed_mass(i,j,1,ised)-              &
     &                               (ero_flux-FC(i,0))
                bed(i,j,1,iaged)=time(ng)
              END IF
!
!  Else,  if erosional.
!
            ELSE
!
!  Update bed thickness of top layer.
!  
              bed(i,j,1,ithck)=bed(i,j,1,ithck)-                        &
     &                         (ero_flux-FC(i,0))/                      &
     &                         (Srho(ised,ng)*                          &
     &                          (1.0_r8-bed(i,j,1,iporo)))
# ifdef MASKING
              bed(i,j,1,ithck)=bed(i,j,1,ithck)*rmask(i,j)
# endif 
!
!  Upate bed mass of top layer.
!
              bed_mass(i,j,1,ised)=bed_mass(i,j,1,ised)-                &
     &                             (ero_flux-FC(i,0))
            END IF
!
!  Compute new sedidiment concentrations.
!
            qc(i,1)=qc(i,1)+ero_flux*Hz_inv(i,1)
          END DO
!
!-----------------------------------------------------------------------
!  Update global tracer variables (m Tunits).
!-----------------------------------------------------------------------
!
          DO k=1,N(ng)
            DO i=Istr,Iend
              t(i,j,k,nnew,indx)=Hz(i,j,k)*qc(i,k)
# ifdef MASKING
              t(i,j,k,nnew,indx)=t(i,j,k,nnew,indx)*rmask(i,j)
# endif
            END DO
          END DO
        END DO SED_LOOP
!
!  Upate bed fraction of top layer.
!
        DO i=Istr,Iend
          cff3=0.0_r8
          DO ised=1,NST
             cff3=cff3+bed_mass(i,j,1,ised)
          END DO
          DO ised=1,NST
            bed_frac(i,j,1,ised)=bed_mass(i,j,1,ised)/MAX(cff3,eps)
          END DO
        END DO
!
!  If first time step of deposit, create new layer and combine bottom
!  two bed layers.
!
          DO i=Istr,Iend
            cff=0.0_r8
!
!  Determine if deposition ocurred here.
!
            DO ised=1,NST
               cff=cff+dep_mass(i,ised)
            END DO
            IF (cff.gt.0.0_r8) THEN
              IF (NBED.gt.1) THEN
!
!  Combine bottom layers.
!
                bed(i,j,Nbed,ithck)=bed(i,j,Nbed-1,ithck)+              &
                                    bed(i,j,Nbed  ,ithck)
                bed(i,j,Nbed,iporo)=0.5_r8*(bed(i,j,Nbed-1,iporo)+      &
     &                                      bed(i,j,Nbed,iporo))
                bed(i,j,Nbed,iaged)=0.5_r8*(bed(i,j,Nbed-1,iaged)+      &
     &                                      bed(i,j,Nbed,iaged))
                cff3=0.0_r8
                DO ised=1,NST
                  bed_mass(i,j,Nbed,ised)=bed_mass(i,j,Nbed-1,ised)+    &
     &                                    bed_mass(i,j,Nbed  ,ised)
                  cff3=cff3+bed_mass(i,j,Nbed,ised)
                ENDDO
                cff3=1.0_r8/MAX(cff3,eps)
                DO ised=1,NST
                  bed_frac(i,j,Nbed,ised)=bed_mass(i,j,Nbed,ised)*cff3
                END DO
!
!  Push layers down.
!
                DO k=Nbed-1,2,-1
                  bed(i,j,k,ithck)=bed(i,j,k-1,ithck)
                  bed(i,j,k,iporo)=bed(i,j,k-1,iporo)
                  bed(i,j,k,iaged)=bed(i,j,k-1,iaged)
                  DO ised =1,NST
                    bed_frac(i,j,k,ised)=bed_frac(i,j,k-1,ised)
                    bed_mass(i,j,k,ised)=bed_mass(i,j,k-1,ised)
                  END DO
                END DO
!
!  Set new top layer parameters.
!
                bed(i,j,1,ithck)=0.0_r8
                DO ised=1,NST
                  bed_mass(i,j,1,ised)=0.0_r8
                END DO
              END IF !NBED=1
              cff3=0.0_r8
              DO ised=1,NST
                bed_mass(i,j,1,ised)=bed_mass(i,j,1,ised)+              &
     &                               dep_mass(i,ised)
                cff3=cff3+bed_mass(i,j,1,ised)
                bed(i,j,1,ithck)=bed(i,j,1,ithck)+                      &
     &                           dep_mass(i,ised)/                      &
     &                           (Srho(ised,ng)*                        &
     &                            (1.0_r8-bed(i,j,1,iporo)))
              END DO
              bed(i,j,1,iaged)=time(ng)
              DO ised=1,NST
                bed_frac(i,j,1,ised)=bed_mass(i,j,1,ised)/              &
     &                               MAX(cff3,eps)
!
!  Reset deposition mass.
!
                dep_mass(i,ised)=0.0_r8
              END DO
            END IF
!
!! Re-determine top layer to be at least active layer thickness.
!
!  Ensure top bed layer thickness is greater or equal than active layer
!  thickness. If need to add sed to top layer, then entrain from lower
!  levels. Create new layers at bottom to maintain Nbed.
!
          IF (bottom(i,j,iactv).gt.bed(i,j,1,ithck)) THEN
            IF (Nbed.eq.1) THEN
              bottom(i,j,iactv)=bed(i,j,1,ithck)
            ELSE
              thck_to_add=bottom(i,j,iactv)-bed(i,j,1,ithck)
              thck_avail=0.0_r8
              Ksed=1                                        ! initialize
              DO k=2,Nbed
                IF (thck_avail.lt.thck_to_add) THEN
                  thck_avail=thck_avail+bed(i,j,k,ithck)
                  Ksed=k
                ENDIF
              END DO
!
!  Catch here if there was not enough bed material.
!
              IF (thck_avail.lt.thck_to_add) THEN
                bottom(i,j,iactv)=bed(i,j,1,ithck)+thck_avail
                thck_to_add=thck_avail
              END IF
!
!  Update bed mass of top layer and fractional layer.
!
              cff2=(thck_avail-thck_to_add)/                            &
     &              MAX(bed(i,j,Ksed,ithck),eps)
              DO ised=1,NST
                cff1=0.0_r8
                DO k=1,Ksed
                  cff1=cff1+bed_mass(i,j,k,ised)
                END DO
                bed_mass(i,j,1   ,ised)=cff1-                           &
     &                                  bed_mass(i,j,Ksed,ised)*cff2
                bed_mass(i,j,Ksed,ised)=bed_mass(i,j,Ksed,ised)*cff2
              END DO
!
!  Update thickness of fractional layer ksource_sed.
!
              bed(i,j,Ksed,ithck)=thck_avail-thck_to_add
!
!  Upate bed fraction of top layer.
!
              cff3=0.0_r8
              DO ised=1,NST
                cff3=cff3+bed_mass(i,j,1,ised)
              END DO
              DO ised=1,NST
                bed_frac(i,j,1,ised)=bed_mass(i,j,1,ised)/              &
     &                               MAX(cff3,eps)
              END DO
!
!  Upate bed thickness of top layer.
!
              bed(i,j,1,ithck)=bottom(i,j,iactv)
!
!  Pull all layers closer to the surface.
!
              DO k=Ksed,Nbed
                ks=Ksed-2
                bed(i,j,k-ks,ithck)=bed(i,j,k,ithck)
                bed(i,j,k-ks,iporo)=bed(i,j,k,iporo)
                bed(i,j,k-ks,iaged)=bed(i,j,k,iaged)
                DO ised=1,NST
                  bed_frac(i,j,k-ks,ised)=bed_frac(i,j,k,ised)
                  bed_mass(i,j,k-ks,ised)=bed_mass(i,j,k,ised)
                ENDDO
              END DO
!
!  Add new layers onto the bottom. Split what was in the bottom layer to
!  fill these new empty cells. ("ks" is the number of new layers).
!
              ks=Ksed-2
              cff=1.0_r8/REAL(ks+1,r8)
              DO k=Nbed,Nbed-ks,-1
                bed(i,j,k,ithck)=bed(i,j,Nbed-ks,ithck)*cff
                bed(i,j,k,iaged)=bed(i,j,Nbed-ks,iaged)
                DO ised=1,NST
                  bed_frac(i,j,k,ised)=bed_frac(i,j,Nbed-ks,ised)
                  bed_mass(i,j,k,ised)=MAX(eps,                         &
     &                                     bed_mass(i,j,Nbed-ks,ised)*  &
     &                                     cff)
                ENDDO
              END DO
            END IF  ! Nbed > 1
          END IF  ! increase top bed layer
!
!  Update mean surface properties.
!
          cff1=1.0_r8
          cff2=1.0_r8
          cff3=1.0_r8
          cff4=1.0_r8
          DO ised=1,NST
            cff1=cff1*tau_ce(ised,ng)**bed_frac(i,j,1,ised)
            cff2=cff2*Sd50(ised,ng)**bed_frac(i,j,1,ised)
            cff3=cff3*wsed(ised,ng)**bed_frac(i,j,1,ised)
            cff4=cff4*Srho(ised,ng)**bed_frac(i,j,1,ised)
          END DO
          bottom(i,j,itauc)=cff1
          bottom(i,j,isd50)=cff2
          bottom(i,j,iwsed)=cff3
          bottom(i,j,idens)=cff4
!
        END DO
      END DO J_LOOP
!
!-----------------------------------------------------------------------
!  Apply periodic or gradient boundary conditions to property arrays.
!-----------------------------------------------------------------------
!
      DO ised=1,NST
        CALL bc_r3d_tile (ng, Istr, Iend, Jstr, Jend,                   &
     &                    LBi, UBi, LBj, UBj, 1, Nbed,                  &
     &                    bed_frac(:,:,:,ised))
        CALL bc_r3d_tile (ng, Istr, Iend, Jstr, Jend,                   &
     &                    LBi, UBi, LBj, UBj, 1, Nbed,                  &
     &                    bed_mass(:,:,:,ised))
      END DO
      DO i=1,MBEDP
        CALL bc_r3d_tile (ng, Istr, Iend, Jstr, Jend,                   &
     &                    LBi, UBi, LBj, UBj, 1, Nbed,                  &
     &                    bed(:,:,:,i))
      END DO
      CALL bc_r3d_tile (ng, Istr, Iend, Jstr, Jend,                     &
     &                  LBi, UBi, LBj, UBj, 1, MBOTP, bottom)

      RETURN
      END SUBROUTINE sediment_tile
#endif
      END MODULE sediment_mod
