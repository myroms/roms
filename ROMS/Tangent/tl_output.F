#include "cppdefs.h"
#ifdef TANGENT
      SUBROUTINE tl_output (ng)
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine manages tangent linear model output. It creates output !
!  NetCDF files and writes out data into NetCDF files.  If  requested, !
!  it can create several tangent history files to avoid generating too !
!  large files during a single model run.                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
# ifdef FOUR_DVAR
      USE mod_fourdvar
# endif
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      logical :: ldef

      integer :: ifile, lstr, status

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on output data time wall clock.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, iTLM, 8)
# endif
!
!-----------------------------------------------------------------------
!  If appropriate, process tangent linear history NetCDF file.
!-----------------------------------------------------------------------
!
!  Turn off checking for analytical header files.
!
      IF (Lanafile) THEN
        Lanafile=.FALSE.
      END IF
!
!  Create output tangent NetCDF file or prepare existing file to
!  append new data to it.  Also,  notice that it is possible to
!  create several files during a single model run.
!
      IF (LdefTLM(ng)) THEN
        IF (ndefTLM(ng).gt.0) THEN
          IF (idefTLM(ng).lt.0) THEN
            idefTLM(ng)=((ntstart(ng)-1)/ndefTLM(ng))*ndefTLM(ng)
            IF (idefTLM(ng).lt.iic(ng)-1) THEN
              idefTLM(ng)=idefTLM(ng)+ndefTLM(ng)
            END IF
          END IF
          IF ((iic(ng)-1).eq.idefTLM(ng)) THEN
            idefTLM(ng)=idefTLM(ng)+ndefTLM(ng)
            IF (nTLM(ng).ne.ndefTLM(ng).and.iic(ng).eq.ntstart(ng)) THEN
              idefTLM(ng)=idefTLM(ng)+nTLM(ng)
            END IF
            NrecTLM(ng)=0
            ifile=(iic(ng)-1)/ndefTLM(ng)+1
            IF (Master) THEN
              lstr=LEN_TRIM(TLMbase(ng))
              WRITE (TLMname(ng),10) TLMbase(ng)(1:lstr-3),ifile
  10          FORMAT (a,'_',i4.4,'.nc')
            END IF
            IF (ncTLMid(ng).ne.-1) status=nf90_close(ncTLMid(ng))
            IF ((iic(ng).eq.ntstart(ng)).and.                               &
     &          (MOD(iic(ng)-1,ndefTLM(ng)).gt.0)) THEN
              ldef=ldefout(ng)
            ELSE
              ldef=.TRUE.
            END IF
            CALL tl_def_his (ng, ldef)
            IF (exit_flag.ne.NoError) RETURN
            LwrtTLM(ng)=.TRUE.
!!          LdefTLM(ng)=.FALSE.
          END IF
        ELSE
          IF (iic(ng).eq.ntstart(ng)) THEN
            CALL tl_def_his (ng, ldefout(ng))
            IF (exit_flag.ne.NoError) RETURN
            LwrtTLM(ng)=.TRUE.
            LdefTLM(ng)=.FALSE.
          END IF
        END IF
      END IF
!
!  Write out data into tangent NetCDF file.  Avoid writing initial
!  conditions in perturbation mode computations.
!
      IF (LwrtTLM(ng)) THEN
        IF (LwrtPER(ng)) THEN
          IF ((iic(ng).gt.ntstart(ng)).and.                             &
     &        (MOD(iic(ng)-1,nTLM(ng)).eq.0)) THEN
            CALL tl_wrt_his (ng)
            IF (exit_flag.ne.NoError) RETURN
          END IF
        ELSE
          IF ((MOD(iic(ng)-1,nTLM(ng)).eq.0).and.                       &
     &        ((nrrec(ng).eq.0).or.(iic(ng).ne.ntstart(ng)))) THEN
            CALL tl_wrt_his (ng)
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END IF
      END IF
# if defined FOUR_DVAR
#  ifdef OBSERVATIONS
!
!-----------------------------------------------------------------------
!  If appropriate, process and write model state at observation
!  locations. Compute misfit (model-observations) cost function.
!-----------------------------------------------------------------------
!
      IF (((time(ng)-0.5_r8*dt(ng)).le.ObsTime(ng)).and.                &
     &    (ObsTime(ng).lt.(time(ng)+0.5_r8*dt(ng)))) THEN
        ProcessObs(ng)=.TRUE.
        CALL obs_read (ng, iTLM, .FALSE.)
        CALL obs_write (ng, iTLM)
        CALL obs_cost (ng, iTLM)
      ELSE
        ProcessObs(ng)=.FALSE.
      END IF
#  endif
# endif
# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off output data time wall clock.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, iTLM, 8)
# endif
      RETURN
      END SUBROUTINE tl_output
#else
      SUBROUTINE tl_output
      RETURN
      END SUBROUTINE tl_output
#endif
