#include "cppdefs.h"
#ifdef ADJOINT
      SUBROUTINE ad_output (ng)
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine manages adjoint model output. It creates output NetCDF !
!  files and writes out data into NetCDF files.  If requested,  it can !
!  create several adjoint history files  to avoid generating too large !
!  files during a single model run.                                    !
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
!  Imported variables declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      logical :: ldef, wrtHIS

      integer :: ifile, lstr, status, subs, tile, thread

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on output data time wall clock.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, iADM, 8)
# endif
!
!-----------------------------------------------------------------------
!  If appropriate, process adjoint history NetCDF file.
!-----------------------------------------------------------------------
!
!  Turn off checking for analytical header files.
!
      IF (Lanafile) THEN
        Lanafile=.FALSE.
      END IF
!
!  Create output adjoint NetCDF file or prepare existing file to
!  append new data to it.  Also,  notice that it is possible to
!  create several files during a single model run.
!
      IF (LdefADJ(ng)) THEN
        IF (ndefADJ(ng).gt.0) THEN
          IF (idefADJ(ng).lt.0) idefADJ(ng)=ntstart(ng)-1
          IF ((iic(ng)-1).eq.idefADJ(ng)) THEN
            idefADJ(ng)=idefADJ(ng)-ndefADJ(ng)
            NrecADJ(ng)=0
            ifile=(iic(ng)-1)/ndefADJ(ng)+1
            IF (Master) THEN
              lstr=LEN_TRIM(ADJbase(ng))
              WRITE (ADJname(ng),10) ADJbase(ng)(1:lstr-3),ifile
  10          FORMAT (a,'_',i4.4,'.nc')
            END IF
            IF (ncADJid(ng).ne.-1) status=nf90_close(ncADJid(ng))
            IF ((iic(ng).eq.ntstart(ng)).and.                           &
     &          (MOD(iic(ng)-1,ndefADJ(ng)).gt.0)) THEN
              ldef=ldefout(ng)
            ELSE
              ldef=.TRUE.
            END IF
            CALL ad_def_his (ng, ldef)
            IF (exit_flag.ne.NoError) RETURN
            LwrtADJ(ng)=.TRUE.
!!          LdefADJ(ng)=.FALSE.
          END IF
        ELSE
          IF (iic(ng).eq.ntstart(ng)) THEN
            CALL ad_def_his (ng, ldefout(ng))
            IF (exit_flag.ne.NoError) RETURN
            LwrtADJ(ng)=.TRUE.
            LdefADJ(ng)=.FALSE.
          END IF
        END IF
      END IF
!
!  Write out data into adjoint NetCDF file.
!
      IF (LwrtADJ(ng)) THEN
        IF (LwrtPER(ng)) THEN
          IF ((iic(ng).ne.ntstart(ng)).and.                             &
     &        (MOD(iic(ng)-1,nADJ(ng)).eq.0)) THEN
            CALL ad_wrt_his (ng)
            IF (exit_flag.ne.NoError) RETURN
          END IF
        ELSE
          IF (nADJ(ng).eq.ntimes(ng)) THEN        ! avoid ntstart record
            wrtHIS=(iic(ng).ne.ntstart(ng)).and.                        &
     &             (MOD(iic(ng)-1,nADJ(ng)).eq.0)
          ELSE                                    ! otherwise
            wrtHIS=(MOD(iic(ng)-1,nADJ(ng)).eq.0)
          END IF
          IF (wrtHIS) THEN
            CALL ad_wrt_his (ng)
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END IF
      END IF

# ifdef AVERAGES
!
!-----------------------------------------------------------------------
!  If appropriate, process time-averaged NetCDF file.
!-----------------------------------------------------------------------
!
!  Create output time-averaged NetCDF file or prepare existing file
!  to append new data to it. Also, notice that it is possible to
!  create several files during a single model run.
!
      IF (LdefAVG(ng)) THEN
        IF (ndefAVG(ng).gt.0) THEN
          IF (idefAVG(ng).lt.0) idefAVG(ng)=ntstart(ng)-1
          IF ((iic(ng)-1).eq.idefAVG(ng)) THEN
            idefAVG(ng)=idefAVG(ng)-ndefAVG(ng)+1
            NrecAVG(ng)=0
            ifile=(iic(ng)-1)/ndefAVG(ng)+1
            IF (Master) THEN
              lstr=LEN_TRIM(AVGbase(ng))
              WRITE (AVGname(ng),20) AVGbase(ng)(1:lstr-3),ifile
  20          FORMAT (a,'_',i4.4,'.nc')
            END IF
            IF (ncAVGid(ng).ne.-1) status=nf90_close(ncAVGid(ng))
            IF ((iic(ng).eq.ntstart(ng)).and.                           &
     &          (MOD(iic(ng)-1,ndefAVG(ng)).gt.0)) THEN
              ldef=ldefout(ng)
            ELSE
              ldef=.TRUE.
            END IF
            CALL def_avg (ng, ldef)
            IF (exit_flag.ne.NoError) RETURN
            LwrtAVG(ng)=.TRUE.
!!          LdefAVG(ng)=.FALSE.
          END IF
        ELSE
          IF (iic(ng).eq.ntstart(ng)) THEN
            CALL def_avg (ng, ldefout(ng))
            IF (exit_flag.ne.NoError) RETURN
            LwrtAVG(ng)=.TRUE.
            LdefAVG(ng)=.FALSE.
          END IF
        END IF
      END IF
!
!  Write out data into time-averaged NetCDF file.
!
      IF (LwrtAVG(ng)) THEN
        IF ((iic(ng).ne.ntstart(ng)).and.                               &
     &      (MOD(iic(ng),nAVG(ng)).eq.1)) THEN
          CALL wrt_avg (ng)
          IF (exit_flag.ne.NoError) RETURN
        END IF
      END IF
# endif
# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off output data time wall clock.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, iADM, 8)
# endif
      RETURN
      END SUBROUTINE ad_output
#else
      SUBROUTINE ad_output
      RETURN
      END SUBROUTINE ad_output
#endif
