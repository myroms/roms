#include "cppdefs.h"
      SUBROUTINE get_cycle (ng, ncid, tfldid, ntime, smday, cycle,      &
     &                      clength, tindex, sfday, tstr, tend, tscale)
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine determines relevant parameters for time cycling        !
!  of data from a input NetCDF file.                                   !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng       Nested grid number (integer).                           !
!     ncid     NetCDF file ID.                                         !
!     tfldid   NetCDF field time variable ID.                          !
!     ntime    Size of time dimension.                                 !
!     smday    Starting model day.                                     !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     cycle    Logical switch indicating cycling of input fields.      !
!     clength  Length of field time cycle.                             !
!     tindex   Starting field time index to read.                      !
!     sfday    Starting field day to process.                          !
!     tstr     Data starting (first record) day.                       !
!     tend     Data Ending (last record) day.                          !
!     tscale   Scale to convert time coordinate to day units.          !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_iounits
      USE mod_netcdf
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      logical, intent(out) :: cycle

      integer, intent(in) :: ng, ncid, tfldid, ntime

      integer, intent(out) :: tindex

      real(r8), intent(in) :: smday

      real(r8), intent(out) :: sfday, tend, clength, tscale
!
!  Local variable declarations.
!
      integer :: i, nvatts, nvdims, status, vartype

      integer, dimension(5) :: vdims

      real(r8) :: mday, tnext, tstart, tstr

      character (len=15) :: attnam, varnam
      character (len=40) :: tunits
!
!-----------------------------------------------------------------------
!  Find time cycling parameters, if any.
!-----------------------------------------------------------------------
!
!  Initialize.
!
      cycle=.false.
      tindex=0
      clength=0.0_r8
      sfday=0.0_r8
      tscale=1.0_r8
!
!  Check if time cycling attribute is present and then read in time
!  cycle length.  Check time coordinate units and determine time
!  scale.  The internal processing of all fields requires time in
!  day units.  Check if more than one time record is available.
!
      IF (ntime.gt.1) THEN
        status=nf_inq_var(ncid,tfldid,varnam,vartype,nvdims,vdims,      &
     &                    nvatts)
        IF (status.eq.nf_noerr) THEN
          DO i=1,nvatts
            status=nf_inq_attname(ncid,tfldid,i,attnam)
            IF (status.eq.nf_noerr) THEN
              IF (TRIM(attnam).eq.'cycle_length') THEN
                cycle=.true.
                status=nf_get_att_TYPE(ncid,tfldid,TRIM(attnam),        &
     &                                 clength)
                IF (status.ne.nf_noerr) THEN
                  WRITE (stdout,10) TRIM(attnam)
                  exit_flag=2
                  RETURN
                END IF
              ELSE IF (TRIM(attnam).eq.'units') THEN
                status=nf_get_att_text(ncid,tfldid,TRIM(attnam),        &
     &                                 tunits)
                IF (status.eq.nf_noerr) THEN
                  IF (tunits(1:6).eq.'second') THEN
                    tscale=sec2day
                  END IF
                ELSE
                  WRITE (stdout,10) TRIM(attnam)
                  exit_flag=2
                  RETURN
                END IF
              END IF
            ELSE
              WRITE (stdout,20) TRIM(varnam)
              exit_flag=2
              RETURN
            END IF
          END DO
        ELSE
          WRITE (stdout,30)
          exit_flag=2
          RETURN
        END IF
      END IF
!
!  Read in starting and ending times for input field.
!
      status=nf_get_var1_TYPE(ncid,tfldid,1,tstr)
      IF (status.ne.nf_noerr) THEN
        WRITE (stdout,40) TRIM(varnam)
        exit_flag=2
        RETURN
      END IF
      status=nf_get_var1_TYPE(ncid,tfldid,ntime,tend)
      IF (status.ne.nf_noerr) THEN
        WRITE (stdout,40) TRIM(varnam)
        exit_flag=2
        RETURN
      END IF
!
!  Search for starting field time and time-index.
!
      tstr=tstr*tscale
      tend=tend*tscale
      IF (cycle) THEN
        mday=MOD(smday,clength)
      ELSE
        mday=smday
      END IF
      IF ((mday.lt.tstr).or.(mday.ge.tend)) THEN
        tindex=ntime
        sfday=tend
      ELSE
        tstart=tstr
        DO i=2,ntime
          status=nf_get_var1_TYPE(ncid,tfldid,i,tnext)
          tnext=tnext*tscale
          IF (status.ne.nf_noerr) THEN
            WRITE (stdout,40) TRIM(varnam)
            exit_flag=2
            RETURN
          END IF
          IF ((tstart.le.mday).and.(mday.lt.tnext)) THEN
            tindex=i-1
            sfday=tstart
            goto 50
          END IF
          tstart=tnext
        END DO
  50    continue
      END IF
!
!  If not cycling, stop execution if there is not field data
!  available for current model time.
!
      IF (.not.cycle.and.(ntime.gt.1)) THEN
        IF (smday.lt.tstr) THEN
          WRITE (stdout,60) TRIM(varnam), tstr, smday
          exit_flag=2
          RETURN
        END IF
      END IF
!
  10  format(/,' GET_CYCLE - error while reading attribute: ',a)
  20  format(/,' GET_CYCLE - error while inquiring attributes for',     &
     &        ' variable: ',a)
  30  format(/,' GET_CYCLE - cannot inquire about time variable from',  &
     &        ' input NetCDF file.')
  40  format(/,' GET_CYCLE - error while reading variable: ',a)
  60  format(/,' GET_CYCLE - starting time for variable: ',a,           &
     &       /,13x,'is greater than current model time. ',              &
     &       /,13x,'TSTR = ',f15.4,2x,'TDAYS = ',f15.4)
      RETURN
      END SUBROUTINE get_cycle
