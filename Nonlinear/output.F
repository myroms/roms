#include "cppdefs.h"
      SUBROUTINE output (ng)
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine manages model output. It creates output NetCDF         !
!  files and writes out data into NetCDF files.  If requested,         !
!  it can create several history and/or time-averaged files to         !
!  avoid generating too large files during a single model run.         !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
#ifdef FLOATS
      USE mod_floats
#endif
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars

      implicit none

      integer, intent(in) :: ng

      integer :: ifile, lstr, rem, status

      logical :: ldef
#ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on output data time wall clock.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, 8)
#endif
!
!-----------------------------------------------------------------------
!  If appropriate, process history NetCDF file.
!-----------------------------------------------------------------------
!
!  Create output history NetCDF file or prepare existing file to
!  append new data to it.  Also,  notice that it is possible to
!  create several files during a single model run.
!
      IF (LdefHIS(ng)) THEN
        IF (ndefhis(ng).gt.0) THEN
          IF (idefhis(ng).lt.0) idefhis(ng)=ntstart-1
          IF ((iic(ng)-1).eq.idefhis(ng)) THEN
            idefhis(ng)=idefhis(ng)+ndefhis(ng)
            nrechis(ng)=0
            ifile=(iic(ng)-1)/ndefhis(ng)+1
            IF (Master) THEN
              lstr=len_trim(hisbase(ng))
              WRITE (hisname(ng),10) hisbase(ng)(1:lstr-3),ifile
            END IF
            IF (nchisid(ng).ne.-1) status=nf_close(nchisid(ng))
            IF ((iic(ng).eq.ntstart).and.                               &
     &          (MOD(iic(ng)-1,ndefhis(ng)).gt.0)) THEN
              ldef=ldefout(ng)
            ELSE
              ldef=.true.
            END IF
            CALL def_his (ng, ldef)
            if (exit_flag.ne.0) RETURN
            LwrtHIS(ng)=.true.
          END IF
        ELSE
          IF (iic(ng).eq.ntstart) THEN
            CALL def_his (ng, ldefout(ng))
            if (exit_flag.ne.0) RETURN
            LwrtHIS(ng)=.true.
            LdefHIS(ng)=.false.
          END IF
        END IF
      END IF
!
!  Write out data into history NetCDF file.  Avoid writing initial
!  conditions in perturbation mode computations.
!
      IF (LwrtHIS(ng)) THEN
        IF (ERend.gt.1) THEN
          IF ((iic(ng).gt.ntstart).and.                                 &
     &        (MOD(iic(ng)-1,nhis(ng)).eq.0)) THEN
            CALL wrt_his (ng)
            if (exit_flag.ne.0) RETURN
          END IF
        ELSE
          IF (MOD(iic(ng)-1,nhis(ng)).eq.0) THEN
            CALL wrt_his (ng)
            if (exit_flag.ne.0) RETURN
          END IF
        END IF
      END IF
#ifdef AVERAGES
!
!-----------------------------------------------------------------------
!  If appropriate, process time-average NetCDF file.
!-----------------------------------------------------------------------
!
!  Create output time-averagd NetCDF file or prepare existing file
!  to append new data to it.  Also,  notice that it is possible to
!  create several files during a single model run.
!
      IF (LdefAVG(ng)) THEN
        IF (ndefavg(ng).gt.0) THEN
          IF (idefavg(ng).lt.0) idefavg(ng)=ntstart+navg(ng)-1
          IF ((iic(ng)-1).eq.idefavg(ng)) THEN
            idefavg(ng)=idefavg(ng)+ndefavg(ng)
            ifile=(iic(ng)-navg(ng)-1)/ndefavg(ng)+1
            IF (Master) THEN
              lstr=len_trim(avgbase(ng))
              WRITE (avgname(ng),10) avgbase(ng)(1:lstr-3),ifile
            END IF
            IF (ncavgid(ng).ne.-1) status=nf_close(ncavgid(ng))
            IF ((iic(ng).eq.ntstart+navg(ng)).and.                            &
     &           (MOD(iic(ng)-navg(ng)-1,ndefavg(ng)).ne.0)) THEN
              ldef=ldefout(ng)
            ELSE
              ldef=.true.
            END IF
            CALL def_avg (ng, ldef)
            if (exit_flag.ne.0) RETURN
          END IF
        ELSE
          IF (iic(ng).eq.ntstart) THEN
            CALL def_avg (ng, ldefout(ng))
            if (exit_flag.ne.0) RETURN
            LdefAVG(ng)=.false.
          END IF
        END IF
      END IF
!
!  Write out data into time-averaged NetCDF file.
!
      IF ((iic(ng).gt.ntstart).and.                                     &
     &    (MOD(iic(ng)-1,navg(ng)).eq.0)) THEN
        CALL wrt_avg (ng)
        if (exit_flag.ne.0) RETURN
      END IF
#endif
#ifdef STATIONS
!
!-----------------------------------------------------------------------
!  If appropriate, process stations NetCDF file.
!-----------------------------------------------------------------------
!
      IF (Nstation(ng).gt.0) THEN
!
!  Create output station NetCDF file or prepare existing file to
!  append new data to it.
!
        IF (LdefSTA(ng).and.(iic(ng).eq.ntstart)) THEN
          CALL def_station (ng, ldefout(ng))
          if (exit_flag.ne.0) RETURN
          LdefSTA(ng)=.false.
        END IF
!
!  Write out data into stations NetCDF file.
!
        IF (MOD(iic(ng)-1,nsta(ng)).eq.0) THEN
          CALL wrt_station (ng)
          if (exit_flag.ne.0) RETURN
        END IF
      END IF
#endif
#ifdef FLOATS
!
!-----------------------------------------------------------------------
!  If appropriate, process floats NetCDF file.
!-----------------------------------------------------------------------
!
      IF (Nfloats(ng).gt.0) THEN
!
!  Create output floats NetCDF file or prepare existing file to
!  append new data to it.
!
        IF (LdefFLT(ng)) THEN
          CALL def_floats (ng)
          if (exit_flag.ne.0) RETURN
          LdefFLT(ng)=.false.
        END IF
!
!  Write out data into floats NetCDF file.
!
        IF ((MOD(iic(ng)-1,nflt(ng)).eq.0).and.                         &
     &      ((frrec(ng).eq.0).or.(iic(ng).ne.ntstart))) THEN
          CALL wrt_floats (ng)
          if (exit_flag.ne.0) RETURN
        END IF
      END IF
#endif
!
!-----------------------------------------------------------------------
!  If appropriate, process restart NetCDF file.
!-----------------------------------------------------------------------
!
!  Create output restart NetCDF file or prepare existing file to
!  append new data to it.
!
      IF (LdefRST(ng)) THEN
        CALL def_rst (ng)
        if (exit_flag.ne.0) RETURN
        LdefRST(ng)=.false.
      END IF
!
!  Write out data into restart NetCDF file.
!
      IF ((iic(ng).gt.ntstart).and.                                     &
     &    (MOD(iic(ng)-1,nrst(ng)).eq.0)) THEN
        CALL wrt_rst (ng)
        if (exit_flag.ne.0) RETURN
      END IF
#ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off output data time wall clock.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, 8)
#endif
!
  10  FORMAT (a,'_',i3.3,'.nc')

      RETURN
      END SUBROUTINE output
