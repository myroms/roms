#include "cppdefs.h"
      SUBROUTINE close_io
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
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
        IF (ncrstid(ng).ne.-1) THEN
          status=nf_close(ncrstid(ng))
          ncrstid(ng)=-1
        END IF
        IF (nchisid(ng).ne.-1) THEN
          status=nf_close(nchisid(ng))
          nchisid(ng)=-1
        END IF
#ifdef AVERAGES
        IF (ncavgid(ng).ne.-1) THEN
          status=nf_close(ncavgid(ng))
          ncavgid=-1
        END IF
#endif
#ifdef FLOATS
        IF (ncfltid(ng).ne.-1) THEN
          status=nf_close(ncfltid(ng))
          ncfltid(ng)=-1
        END IF
#endif
        IF (lcycle(ng)) THEN
          IF (nrecrst(ng).gt.1) THEN
            nrecrst(ng)=2
          ELSE
            nrecrst(ng)=1
          END IF
        END IF
!
!  Report number of time records written.
!
        IF (Master) WRITE (stdout,10) ng, nrechis(ng), ng, nrecrst(ng)
#ifdef AVERAGES
        IF (Master) WRITE (stdout,20) ng, nrecavg(ng)
#endif
#ifdef STATIONS
        IF (Master) WRITE (stdout,30) ng, nrecsta(ng)
#endif
      END DO
!
!-----------------------------------------------------------------------
!  If applicable, report internal exit errors.
!-----------------------------------------------------------------------
!
      IF (exit_flag.eq.0) THEN
        CALL get_date (date_str)
        IF (Master) WRITE (stdout,40) TRIM(date_str)
      ELSE IF (exit_flag.eq.1) THEN
        IF (Master) WRITE (stdout,50)
      ELSE IF (exit_flag.eq.2) THEN
        IF (Master) WRITE (stdout,60)
      ELSE IF (exit_flag.eq.3) THEN
        IF (Master) WRITE (stdout,70)
      ELSE IF (exit_flag.eq.4) THEN
        IF (Master) WRITE (stdout,80)
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
 60   FORMAT (/,' ERROR: Abnormal termination: NetCDF INPUT.')
 70   FORMAT (/,' ERROR: Abnormal termination: NetCDF OUTPUT.')
 80   FORMAT (/,' ERROR: OPENCDF - Can not open NetCDF file.')
 90   FORMAT (/,' ERROR: CHECKDEFS - Illegal model configuration.')
100   FORMAT (/,' ERROR: PARAM - Illegal domain partition.')

      RETURN
      END SUBROUTINE close_io
