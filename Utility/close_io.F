#include "cppdefs.h"
      SUBROUTINE close_io
!
!=======================================================================
!  Copyright (c) 2005 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
! This subroutine flushes and closes all IO files.                     !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
!
      implicit none
!
!  Local variable declarations.
!
      integer :: MyError, ng, status
!
!-----------------------------------------------------------------------
!  Close output NetCDF files. Set file indices to closed state.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        IF (ncRSTid(ng).ne.-1) THEN
          status=nf_close(ncRSTid(ng))
          ncRSTid(ng)=-1
        END IF
        IF (ncHISid(ng).ne.-1) THEN
          status=nf_close(ncHISid(ng))
          ncHISid(ng)=-1
        END IF
#ifdef AVERAGES
        IF (ncAVGid(ng).ne.-1) THEN
          status=nf_close(ncAVGid(ng))
          ncAVGid=-1
        END IF
#endif
#ifdef FLOATS
        IF (ncFLTid(ng).ne.-1) THEN
          status=nf_close(ncFLTid(ng))
          ncFLTid(ng)=-1
        END IF
#endif
#if defined FORWARD_READ || defined FORWARD_WRITE
        IF (ncFWDid(ng).ne.-1) THEN
          status=nf_close(ncFWDid(ng))
          ncFWDid(ng)=-1
        END IF
#endif
        IF (LcycleRST(ng)) THEN
          IF (NrecRST(ng).gt.1) THEN
            NrecRST(ng)=2
          ELSE
            NrecRST(ng)=1
          END IF
        END IF
!
!  Report number of time records written.
!
        IF (Master) WRITE (stdout,10) ng, NrecHIS(ng), ng, NrecRST(ng)
#ifdef AVERAGES
        IF (Master) WRITE (stdout,20) ng, NrecAVG(ng)
#endif
#ifdef STATIONS
        IF (Master) WRITE (stdout,30) ng, NrecSTA(ng)
#endif
      END DO
!
!-----------------------------------------------------------------------
!  If applicable, report internal exit errors.
!-----------------------------------------------------------------------
!
      IF (Master.and.blowup.ne.0) THEN
        WRITE (stdout,50)
      END IF
      IF (exit_flag.eq.NoError) THEN
        CALL get_date (date_str)
        IF (Master) WRITE (stdout,40) TRIM(date_str)
      ELSE IF (exit_flag.eq.1) THEN
        IF (Master) WRITE (stdout,50)
      ELSE IF (exit_flag.eq.2) THEN
        IF (Master) WRITE (stdout,60) nf_strerror(ioerror)
      ELSE IF (exit_flag.eq.3) THEN
        IF (Master) WRITE (stdout,70) nf_strerror(ioerror)
      ELSE IF (exit_flag.eq.4) THEN
        IF (Master) WRITE (stdout,80) nf_strerror(ioerror)
      ELSE IF (exit_flag.eq.5) THEN
        IF (Master) WRITE (stdout,90)
      ELSE IF (exit_flag.eq.6) THEN
        IF (Master) WRITE (stdout,100)
      END IF
!
 10   FORMAT (/,' ROMS/TOMS - number of time records written in',       &
     &        ' history  file, Grid ',i2.2,', ',i8.8,/,                 &
     &        13x,'number of time records written in restart ',         &
     &        ' file, Grid ',i2.2,', ',i8.8)
 20   FORMAT (13x,'number of time records written in averages ',        &
     &        'file, Grid ',i2.2,', ',i8.8)
 30   FORMAT (13x,'number of time records written in station ',         &
     &        ' file, Grid ',i2.2,', ',i8.8)
 40   FORMAT (/,' ROMS/TOMS: DONE... ',a)
 50   FORMAT (/,' MAIN: Abnormal termination: BLOWUP.')
 60   FORMAT (/,' ERROR: Abnormal termination: NetCDF INPUT.',/,        &
     &          ' REASON: ',a)
 70   FORMAT (/,' ERROR: Abnormal termination: NetCDF OUTPUT.',/,       &
     &          ' REASON: ',a)
 80   FORMAT (/,' ERROR: OPENCDF - Can not open NetCDF file.',/,        &
     &          ' REASON: ',a)
 90   FORMAT (/,' ERROR: CHECKDEFS - Illegal model configuration.')
100   FORMAT (/,' ERROR: PARAM - Illegal domain partition.')

      RETURN
      END SUBROUTINE close_io
