#include "cppdefs.h"
#ifdef FLOATS
      SUBROUTINE wrt_floats (ng)
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
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
      if (exit_flag.ne.0) RETURN
!
!  Set time record index.
!
      tfltindx(ng)=tfltindx(ng)+1
      nrecflt(ng)=nrecflt(ng)+1
!
      start(1)=1
      total(1)=Nfloats(ng)
      start(2)=tfltindx(ng)
      total(2)=1
!
!  Write out model time (s).
!
      IF (Master) THEN
        status=nf_put_var1_TYPE(ncfltid(ng), fltVid(idtime,ng),         &
     &                          tfltindx(ng), time(ng))
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,10) TRIM(Vname(1,idtime)), tfltindx(ng)
          exit_flag=3
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
        status=nf_put_vara_TYPE(ncfltid(ng), fltVid(idXgrd,ng),         &
     &                          start, total, Tout)
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,10) 'Xgrid', tfltindx(ng)
          exit_flag=3
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
        status=nf_put_vara_TYPE(ncfltid(ng), fltVid(idYgrd,ng),         &
     &                          start, total, Tout)
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,10) 'Ygrid', tfltindx(ng)
          exit_flag=3
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
        status=nf_put_vara_TYPE(ncfltid(ng), fltVid(idZgrd,ng),         &
     &                          start, total, Tout)
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,10) 'Zgrid', tfltindx(ng)
          exit_flag=3
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
        status=nf_put_vara_TYPE(ncfltid(ng), fltVid(idglon,ng),         &
     &                          start, total, Tout)
        IF (status.ne.nf_noerr) THEN
          IF (spherical) THEN
            WRITE (stdout,10) 'lon', tfltindx(ng)
          ELSE
            WRITE (stdout,10) 'x', tfltindx(ng)
          END IF
          exit_flag=3
          RETURN
        END IF
        DO l=1,Nfloats(ng)
          Tout(l)=FLT(ng)%track(iflat,nf(ng),l)
        END DO
        status=nf_put_vara_TYPE(ncfltid(ng), fltVid(idglat,ng),         &
     &                          start, total, Tout)
        IF (status.ne.nf_noerr) THEN
          IF (spherical) THEN
            WRITE (stdout,10) 'lat', tfltindx(ng)
          ELSE
            WRITE (stdout,10) 'y', tfltindx(ng)
          END IF
          exit_flag=3
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
        status=nf_put_vara_TYPE(ncfltid(ng), fltVid(iddpth,ng),         &
     &                          start, total, Tout)
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,10) 'depth', tfltindx(ng)
          exit_flag=3
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
        status=nf_put_vara_TYPE(ncfltid(ng), fltVid(idDano,ng),         &
     &                          start, total, Tout)
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,10) TRIM(Vname(1,idDano)), tfltindx(ng)
          exit_flag=3
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
          status=nf_put_vara_TYPE(ncfltid(ng), fltTid(itrc,ng),         &
     &                            start, total, Tout)
          IF (status.ne.nf_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTvar(itrc))), tfltindx(ng)
            exit_flag=3
            RETURN
          END IF
        END DO
      END IF
# endif
!
!  Synchronize float NetCDF file to disk.
!
      IF (Master) THEN
        status=nf_sync(ncfltid(ng))
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,20)
          exit_flag=3
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
