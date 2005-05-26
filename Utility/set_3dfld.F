#include "cppdefs.h"
      MODULE set_3dfld_mod
#ifdef SOLVE3D
!
!=======================================================================
!  Copyright (c) 2005 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine time-interpolates requested 3D field from snapshots    !
!  of input data.                                                      !
!                                                                      !
!=======================================================================
!
      implicit none

      CONTAINS
!
!***********************************************************************
      SUBROUTINE set_3dfld_tile (ng, model, Istr, Iend, Jstr, Jend,     &
     &                           ifield, LBi, UBi, LBj, UBj, LBk, UBk,  &
     &                           Finp, Fout, update)
!***********************************************************************
!
      USE mod_param
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars

# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
!
      USE exchange_3d_mod
# endif
!
!  Imported variable declarations.
!
      logical, intent(out) :: update

      integer, intent(in) :: ng, model, Istr, Iend, Jstr, Jend, ifield
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk

# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: Finp(LBi:,LBj:,:,:)
      real(r8), intent(out) :: Fout(LBi:,LBj:,:)
# else
      real(r8), intent(in) :: Finp(LBi:UBi,LBj:UBj,LBk:UBk,2)
      real(r8), intent(out) :: Fout(LBi:UBi,LBj:UBj,LBk:UBk)
# endif
!
!  Local variable declarations.
!
      logical :: Lgrided, Lonerec

      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: Tindex, gtype, i, it1, it2, j, k

      real(r8) :: Fval, Tend, fac, fac1, fac2

# include "set_bounds.h"
!
!--------------------------------------------------------------------
!  Set-up requested field for current tile.
!--------------------------------------------------------------------
!
!  Get requested field information from global storage.
!
      Lgrided=Linfo(1,ifield,ng)
      Lonerec=Linfo(3,ifield,ng)
      gtype  =Iinfo(1,ifield,ng)
      Tindex =Iinfo(5,ifield,ng)
      update=.TRUE.
!
!  Set linear-interpolation factors.
!
      it1=3-Tindex
      it2=Tindex
# ifdef CRAY
      fac1=ANINT(Tintrp(it2,ifield,ng)-time(ng))
      fac2=ANINT(time(ng)-Tintrp(it1,ifield,ng))
# else
      fac1=ANINT(Tintrp(it2,ifield,ng)-time(ng),r8)
      fac2=ANINT(time(ng)-Tintrp(it1,ifield,ng),r8)
# endif
!
!  Load time-invariant data. Time interpolation is not necessary.
!
      IF (Lonerec) THEN
        IF (Lgrided) THEN
          DO k=LBk,UBk
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                Fout(i,j,k)=Finp(i,j,k,Tindex)
              END DO
            END DO
          END DO
        ELSE
          Fval=Fpoint(Tindex,ifield,ng)
          DO k=LBk,UBk
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                Fout(i,j,k)=Fval
              END DO
            END DO
          END DO
        END IF
!
!  Time-interpolate from grided or point data.
!
      ELSE IF (((fac1*fac2).ge.0.0_r8).and.                             &
     &        ((fac1+fac2).gt.0.0_r8)) THEN
        fac=1.0_r8/(fac1+fac2)
        fac1=fac*fac1
        fac2=fac*fac2
        IF (Lgrided) THEN
          DO k=LBk,UBk
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                Fout(i,j,k)=fac1*Finp(i,j,k,it1)+fac2*Finp(i,j,k,it2)
              END DO
            END DO
          END DO
        ELSE
          Fval=fac1*Fpoint(it1,ifield,ng)+fac2*Fpoint(it2,ifield,ng)
          DO k=LBk,UBk
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                Fout(i,j,k)=Fval
              END DO
            END DO
          END DO
        END IF
!
!  Activate synchronization flag if a new time record needs to be
!  read in at the next time step.
!
        IF ((time(ng)+dt(ng)).gt.Tintrp(it2,ifield,ng)) THEN
          IF (SOUTH_WEST_TEST) synchro_flag(ng)=.TRUE.
        END IF
!
!  Unable to set-up requested field.  Activate error flag to quit.
!
      ELSE
        IF (SOUTH_WEST_TEST) THEN
          Tend=Finfo(2,ifield,ng)
          WRITE (stdout,10) TRIM(Vname(5,ifield)), tdays(ng), Tend,     &
     &                      Tintrp(it1,ifield,ng)*sec2day,              &
     &                      Tintrp(it2,ifield,ng)*sec2day
  10      FORMAT (/,' SET_3DFLD  - current model time',                 &
     &            ' exceeds ending value for variable: ',a,             &
     &            /,14x,'TDAYS   = ',f15.4,2x,'TEND    = ',f15.4,       &
     &            /,14x,'TINTRP1 = ',f15.4,2x,'TINTRP2 = ',f15.4)
          exit_flag=2
          update=.FALSE.
        END IF
      END IF
# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
!
!  Exchange boundary data.
!
      IF (update) THEN
        IF (gtype.eq.r3dvar) THEN
          CALL exchange_r3d_tile (ng, model, Istr, Iend, Jstr, Jend,    &
     &                            LBi, UBi, LBj, UBj, LBk, UBk,         &
     &                            NghostPoints,                         &
     &                            Fout)
        ELSE IF (gtype.eq.u3dvar) THEN
          CALL exchange_u3d_tile (ng, model, Istr, Iend, Jstr, Jend,    &
     &                            LBi, UBi, LBj, UBj, LBk, UBk,         &
     &                            NghostPoints,                         &
     &                            Fout)
        ELSE IF (gtype.eq.v3dvar) THEN
          CALL exchange_v3d_tile (ng, model, Istr, Iend, Jstr, Jend,    &
     &                            LBi, UBi, LBj, UBj, LBk, UBk,         &
     &                            NghostPoints,                         &
     &                            Fout)
        ELSE IF (gtype.eq.w3dvar) THEN
          CALL exchange_w3d_tile (ng, model, Istr, Iend, Jstr, Jend,    &
     &                            LBi, UBi, LBj, UBj, LBk, UBk,         &
     &                            NghostPoints,                         &
     &                            Fout)
        END IF
      END IF
# endif
      RETURN
      END SUBROUTINE set_3dfld_tile
#endif
      END MODULE set_3dfld_mod
