#include "cppdefs.h"
#ifdef FLOATS
      SUBROUTINE wrt_floats (ng)
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2007 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine writes simulated drifter trajectories into floats   !
!  NetCDF file.                                                        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_floats
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
      USE mod_stepping
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      integer :: itrc, l, status

      integer, dimension(2) :: start, total

      real(r8), dimension(Nfloats(ng)) :: Tout
!
!-----------------------------------------------------------------------
!  Write out station data at RHO-points.
!-----------------------------------------------------------------------
!
      if (exit_flag.ne.NoError) RETURN
!
!  Set time record index.
!
      tFLTindx(ng)=tFLTindx(ng)+1
      NrecFLT(ng)=NrecFLT(ng)+1
!
      start(1)=1
      total(1)=Nfloats(ng)
      start(2)=tFLTindx(ng)
      total(2)=1
!
!  Write out model time (s).
!
      IF (Master) THEN
        status=nf_put_var1_TYPE(ncFLTid(ng), fltVid(idtime,ng),         &
     &                          tFLTindx(ng), time(ng))
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,10) TRIM(Vname(1,idtime)), tFLTindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out floats X-grid locations.
!
      IF (Master) THEN
        DO l=1,Nfloats(ng)
          IF (FLT(ng)%bounded(l)) THEN
            Tout(l)=FLT(ng)%track(ixgrd,nf(ng),l)
          ELSE
            Tout(l)=spval
          END IF
        END DO
        status=nf_put_vara_TYPE(ncFLTid(ng), fltVid(idXgrd,ng),         &
     &                          start, total, Tout)
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,10) 'Xgrid', tFLTindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out floats Y-grid locations.
!
      IF (Master) THEN
        DO l=1,Nfloats(ng)
          IF (FLT(ng)%bounded(l)) THEN
            Tout(l)=FLT(ng)%track(iygrd,nf(ng),l)
          ELSE
            Tout(l)=spval
          END IF
        END DO
        status=nf_put_vara_TYPE(ncFLTid(ng), fltVid(idYgrd,ng),         &
     &                          start, total, Tout)
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,10) 'Ygrid', tFLTindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# ifdef SOLVE3D
!
!  Write out floats Z-grid locations.
!
      IF (Master) THEN
        DO l=1,Nfloats(ng)
          IF (FLT(ng)%bounded(l)) THEN
            Tout(l)=FLT(ng)%track(izgrd,nf(ng),l)
          ELSE
            Tout(l)=spval
          END IF
        END DO
        status=nf_put_vara_TYPE(ncFLTid(ng), fltVid(idZgrd,ng),         &
     &                          start, total, Tout)
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,10) 'Zgrid', tFLTindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
!
!  Write out floats (lon,lat) or (x,y) locations.
!
      IF (Master) THEN
        DO l=1,Nfloats(ng)
          Tout(l)=FLT(ng)%track(iflon,nf(ng),l)
        END DO
        status=nf_put_vara_TYPE(ncFLTid(ng), fltVid(idglon,ng),         &
     &                          start, total, Tout)
        IF (status.ne.nf_noerr) THEN
          IF (spherical) THEN
            WRITE (stdout,10) 'lon', tFLTindx(ng)
          ELSE
            WRITE (stdout,10) 'x', tFLTindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        DO l=1,Nfloats(ng)
          Tout(l)=FLT(ng)%track(iflat,nf(ng),l)
        END DO
        status=nf_put_vara_TYPE(ncFLTid(ng), fltVid(idglat,ng),         &
     &                          start, total, Tout)
        IF (status.ne.nf_noerr) THEN
          IF (spherical) THEN
            WRITE (stdout,10) 'lat', tFLTindx(ng)
          ELSE
            WRITE (stdout,10) 'y', tFLTindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# ifdef SOLVE3D
!
!  Write out floats depths.
!
      IF (Master) THEN
        DO l=1,Nfloats(ng)
          IF (FLT(ng)%bounded(l)) THEN
            Tout(l)=FLT(ng)%track(idpth,nf(ng),l)
          ELSE
            Tout(l)=spval
          END IF
        END DO
        status=nf_put_vara_TYPE(ncFLTid(ng), fltVid(iddpth,ng),         &
     &                          start, total, Tout)
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,10) 'depth', tFLTindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out density anomaly.
!
      IF (Master) THEN
        DO l=1,Nfloats(ng)
          IF (FLT(ng)%bounded(l)) THEN
            Tout(l)=FLT(ng)%track(ifden,nf(ng),l)
          ELSE
            Tout(l)=spval
          END IF
        END DO
        status=nf_put_vara_TYPE(ncFLTid(ng), fltVid(idDano,ng),         &
     &                          start, total, Tout)
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,10) TRIM(Vname(1,idDano)), tFLTindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out tracer type variables.
!
      IF (Master) THEN
        DO itrc=1,NT(ng)
          DO l=1,Nfloats(ng)
            IF (FLT(ng)%bounded(l)) THEN
              Tout(l)=FLT(ng)%track(itrc+10,nf(ng),l)
            ELSE
              Tout(l)=spval
            END IF
          END DO
          status=nf_put_vara_TYPE(ncFLTid(ng), fltTid(itrc,ng),         &
     &                            start, total, Tout)
          IF (status.ne.nf_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTvar(itrc))), tFLTindx(ng)
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END DO
      END IF
# endif
!
!-----------------------------------------------------------------------
!  Synchronize floats NetCDF file to disk.
!-----------------------------------------------------------------------
!
      IF (Master) THEN
        status=nf_sync(ncFLTid(ng))
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,20)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
  10  FORMAT (/,' WRT_FLOATS - error while writing variable: ',a,/,     &
     &        14x,'into floats NetCDF file for time record: ',i4)
  20  FORMAT (/,' WRT_FLOATS - unable to synchronize stations',         &
     &        1x,'NetCDF file to disk.')
#else
      SUBROUTINE wrt_floats
#endif
      RETURN
      END SUBROUTINE wrt_floats
