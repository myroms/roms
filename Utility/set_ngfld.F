#include "cppdefs.h"
      SUBROUTINE set_ngfld (ng, ifield, LBi, UBi, UBj, Istr, Iend,      &
     &                      Jrec, Finp, Fout, update)
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine time-interpolates requested non-grided field from      !
!  time snapshots of input data.                                       !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number (integer).                         !
!     ifield     Field ID.                                             !
!     LBi        Finp/Fout 1st dimension lower-bound value.            !
!     UBi        Finp/Fout 1st dimension upper-bound value.            !
!     UBj        Finp/Fout 2nd dimension upper-bound value, if any.    !
!                  Otherwise, a value of one is expected.              !
!     Istr       Starting location to process in the 1st dimension.    !
!     Iend       Ending location to process in the 1st dimension.      !
!     Jrec       Number of records to process in the 2nd dimenision,   !
!                  if any, Otherwise, a value of one is expected.      !
!     Finp       Latest two-snapshopts of field to interpolate.        !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Fout       Interpolated field.                                   !
!     update     Switch indicating successful interpolation.           !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      logical, intent(out) :: update

      integer, intent(in) :: ng, ifield
      integer, intent(in) :: LBi, UBi, UBj, Istr, Iend, Jrec

      real(r8), intent(in) :: Finp(LBi:UBi,UBj,2)

      real(r8), intent(out) :: Fout(LBi:UBi,UBj)
!
!  Local variable declarations.
!
      logical :: Lonerec

      integer :: Tindex, i, it1, it2, j

      real(r8) :: Tend, fac, fac1, fac2
!
!----------------------------------------------------------------------
!  Set up requested field from data snapshots.
!----------------------------------------------------------------------
!
!  Get requested field information from global storage.
!
      Lonerec=Linfo(3,ifield,ng)
      Tindex=Iinfo(5,ifield,ng)
      update=.true.
!
!  Set linear-interpolation factors.
!
      it1=3-Tindex
      it2=Tindex
      fac1=Tintrp(it2,ifield,ng)-time(ng)
      fac2=time(ng)-Tintrp(it1,ifield,ng)
!
!  Load time-invariant data. Time interpolation is not necessary.
!
      IF (Lonerec) THEN
        DO j=1,Jrec
          DO i=Istr,Iend
            Fout(i,j)=Finp(i,j,Tindex)
          END DO
        END DO
!
!  Time-interpolate.
!
      ELSE IF (((fac1*fac2).ge.0.0_r8).and.(fac1+fac2).gt.0.0_r8) THEN
        fac=1.0_r8/(fac1+fac2)
        fac1=fac*fac1
        fac2=fac*fac2
        DO j=1,Jrec
          DO i=Istr,Iend
            Fout(i,j)=fac1*Finp(i,j,it1)+fac2*Finp(i,j,it2)
          END DO
        END DO
!
!  Activate synchronization flag if a new time record needs to be
!  read in at the next time step.
!
        IF ((time(ng)+dt(ng)).gt.Tintrp(it2,ifield,ng)) THEN
          synchro_flag(ng)=.true.
        END IF
!
!  Unable to interpolate field.  Activate error flag to quit.
!
      ELSE
        Tend=Finfo(2,ifield,ng)
        write(stdout,10) TRIM(Vname(5,ifield)), tdays(ng), Tend,        &
     &                   Tintrp(it1,ifield,ng)*sec2day,                 &
     &                   Tintrp(it2,ifield,ng)*sec2day
  10    format(/,' SET_NGFLD  - current model time',                    &
     &         ' exceeds ending value for variable: ',a,                &
     &         /,14x,'TDAYS   = ',f15.4,2x,'TEND    = ',f15.4,          &
     &         /,14x,'TINTRP1 = ',f15.4,2x,'TINTRP2 = ',f15.4)
        exit_flag=2
        update=.false.
      END IF
      RETURN
      END SUBROUTINE set_ngfld
