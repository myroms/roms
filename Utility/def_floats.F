#include "cppdefs.h"
#ifdef FLOATS
      SUBROUTINE def_floats (ng)
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine creates FLOATS NetCDF file, it defines dimensions,     !
!  attributes, and variables.                                          !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_grid
      USE mod_floats
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
# ifdef SEDIMENT
      USE mod_sediment
# endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      integer, parameter :: Natt = 20

      logical :: got_var(-6:NV)

      integer :: fltdim, i, ifield, itrc, j, l, nrec
      integer :: recdim, status

      integer :: DimIDs(28), fgrd(2), start(2), total(2)

      integer :: def_info, def_var, wrt_info

# ifdef DISTRIBUTE
      integer :: Npts

      real(r8), parameter :: Fspv = 0.0_r8
      real(r8), dimension(Nfloats(ng)*NFV(ng)) :: Fwrk
# endif

      real(r8) :: Aval(4), Tinp(Nfloats(ng))

      character (len=80) :: Vinfo(Natt)
      character (len=80) :: dimnam, fname, ncname
!
!=======================================================================
!  Create a new floats data file.
!=======================================================================
!
!  Deactivate file creation, if a restart run.
!
      IF (frrec(ng).ne.0) LdefFLT=.false.
!
!  Create floats NetCDF file.
!
      IF (exit_flag.ne.0) RETURN
      IF (LdefFLT(ng).and.Master) THEN
        ncname=fltname(ng)
        status=nf_create(TRIM(ncname),nf_clobber,ncfltid(ng))
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,10) TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
      END IF
#ifdef DISTRIBUTE
      CALL mp_bcasti (ng, ncfltid(ng), 1)
#endif
!
!-----------------------------------------------------------------------
!  Define the dimensions of staggered fields.
!-----------------------------------------------------------------------
!
      IF (LdefFLT(ng).and.Master) THEN
# ifdef SOLVE3D
        status=nf_def_dim(ncfltid(ng),'s_rho'   ,N(ng),                 &
     &                    DimIDs( 9))
        status=nf_def_dim(ncfltid(ng),'s_w'     ,N(ng)+1,               &
     &                    DimIDs(10))
        status=nf_def_dim(ncfltid(ng),'tracer'  ,NT(ng),                &
     &                    DimIDs(11))
#  ifdef SEDIMENT
        status=nf_def_dim(ncfltid(ng),'Nbed'    ,Nbed,                  &
     &                    DimIDs(16))
#  endif
#  ifdef ECOSIM
        status=nf_def_dim(ncfltid(ng),'Nphy'    ,Nphy,                  &
     &                    DimIDs(25))
        status=nf_def_dim(ncfltid(ng),'Nbac'    ,Nbac,                  &
     &                    DimIDs(26))
        status=nf_def_dim(ncfltid(ng),'Ndom'    ,Ndom,                  &
     &                    DimIDs(27))
        status=nf_def_dim(ncfltid(ng),'Nfec'    ,Nfec,                  &
     &                    DimIDs(28))
#  endif
# endif
        status=nf_def_dim(ncfltid(ng),'drifter' ,Nfloats(ng),           &
     &                    DimIDs(15))
        status=nf_def_dim(ncfltid(ng),'boundary',4,                     &
     &                    DimIDs(14))
        status=nf_def_dim(ncfltid(ng),TRIM(ADJUSTL(Vname(5,idtime))),   &
     &                    nf_unlimited,DimIDs(12))
        recdim=DimIDs(12)
        fltdim=DimIDs(15)
!
!  Define dimension vectors for point variables.
!
        fgrd(1)=DimIDs(15)
        fgrd(2)=DimIDs(12)
!
!  Initialize unlimited time record dimension.
!
        tfltindx(ng)=0
!
!  Initialize local information variable arrays.
!
        DO i=1,Natt
          DO j=1,80
            Vinfo(i)(j:j)=' '
          END DO
        END DO
        DO i=1,4
          Aval(i)=0.0_r8
        END DO
!
!-----------------------------------------------------------------------
!  Define time-recordless information variables.
!-----------------------------------------------------------------------
!
        status=def_info(ng,ncfltid(ng),ncname,DimIDs)
        if (exit_flag.ne.0) RETURN
!
!-----------------------------------------------------------------------
!  Define variables and their attributes.
!-----------------------------------------------------------------------
!
!  Define model time.
!
        Vinfo( 1)=Vname(1,idtime)
        Vinfo( 2)=Vname(2,idtime)
        IF (INT(time_ref).eq.-2) THEN
          Vinfo( 3)='seconds since 1968-05-23 00:00:00 GMT'
          Vinfo( 4)='modified Julian day number'
          Vinfo(18)='add_offset'
          Aval(1)=jul_off
        ELSE IF (INT(time_ref).eq.-1) THEN
          Vinfo( 3)='seconds since 0000-01-01 00:00:00'
          Vinfo( 4)='360 days in every year, 30 days in each month'
        ELSE IF (INT(time_ref).eq.0) THEN
          Vinfo( 3)='seconds since 0000-01-01 00:00:00'
          Vinfo( 4)='365.25 days in every year'
        ELSE IF (time_ref.gt.0.0_r8) THEN
          Vinfo( 3)='seconds since '//r_text
        END IF
        Vinfo(14)=Vname(4,idtime)
        status=def_var(ncfltid(ng),fltVid(idtime,ng),NF_TYPE,           &
     &                 1,recdim,Aval,Vinfo,ncname)
!
!  Define floats X-grid locations.
!
        Vinfo( 1)='Xgrid'
        Vinfo( 2)='x-grid floats locations'
        Vinfo( 3)='nondimensional'
        Vinfo( 5)='valid_min'
        Vinfo( 6)='valid_max'
        Aval(2)=0.0_r8
        Aval(3)=REAL(Lm(ng)+1,r8)
        Vinfo(14)='Xgrid, scalar, series'
        Vinfo(16)=Vname(1,idtime)
        Vinfo(17)='missing_value'
        Aval(4)=spval
        status=def_var(ncfltid(ng),fltVid(idXgrd,ng),NF_FOUT,           &
     &                 2,fgrd,Aval,Vinfo,ncname)
!
!  Define floats Y-grid locations.
!
        Vinfo( 1)='Ygrid'
        Vinfo( 2)='Y-grid floats locations'
        Vinfo( 3)='nondimensional'
        Vinfo( 5)='valid_min'
        Vinfo( 6)='valid_max'
        Aval(2)=0.0_r8
        Aval(3)=REAL(Mm(ng)+1,r8)
        Vinfo(14)='Ygrid, scalar, series'
        Vinfo(16)=Vname(1,idtime)
        Vinfo(17)='missing_value'
        Aval(4)=spval
        status=def_var(ncfltid(ng),fltVid(idYgrd,ng),NF_FOUT,           &
     &                 2,fgrd,Aval,Vinfo,ncname)
# ifdef SOLVE3D
!
!  Define floats Z-grid locations.
!
        Vinfo( 1)='Zgrid'
        Vinfo( 2)='Z-grid floats locations'
        Vinfo( 3)='nondimensional'
        Vinfo( 5)='valid_min'
        Vinfo( 6)='valid_max'
        Aval(2)=0.0_r8
        Aval(3)=REAL(N(ng),r8)
        Vinfo(14)='Zgrid, scalar, series'
        Vinfo(16)=Vname(1,idtime)
        Vinfo(17)='missing_value'
        Aval(4)=spval
        status=def_var(ncfltid(ng),fltVid(idZgrd,ng),NF_FOUT,           &
     &                 2,fgrd,Aval,Vinfo,ncname)
# endif
!
!  Define floats (lon,lat) or (x,y) locations.
!
        IF (spherical) THEN
          Vinfo( 1)='lon'
          Vinfo( 2)='longitude of floats trajectories'
          Vinfo( 3)='degree_east'
          Vinfo( 5)='valid_min'
          Vinfo( 6)='valid_max'
          Vinfo(14)='lon, scalar, series'
          Vinfo(16)=Vname(1,idtime)
          Vinfo(17)='missing_value'
          Aval(2)=-180.0_r8
          Aval(3)=180.0_r8
          Aval(4)=spval
          status=def_var(ncfltid(ng),fltVid(idglon,ng),NF_FOUT,         &
     &                   2,fgrd,Aval,Vinfo,ncname)
          Vinfo( 1)='lat'
          Vinfo( 2)='latitude of floats trajectories'
          Vinfo( 3)='degree_north'
          Vinfo( 5)='valid_min'
          Vinfo( 6)='valid_max'
          Vinfo(14)='lat, scalar, series'
          Vinfo(16)=Vname(1,idtime)
          Vinfo(17)='missing_value'
          Aval(2)=-90.0_r8
          Aval(3)=90.0_r8
          Aval(4)=spval
          status=def_var(ncfltid(ng),fltVid(idglat,ng),NF_FOUT,         &
     &                   2,fgrd,Aval,Vinfo,ncname)
       ELSE
          Vinfo( 1)='x'
          Vinfo( 2)='x-location of floats trajectories'
          Vinfo( 3)='meter'
          Vinfo(14)='x, scalar, series'
          Vinfo(16)=Vname(1,idtime)
          Vinfo(17)='missing_value'
          Aval(4)=spval
          status=def_var(ncfltid(ng),fltVid(idglon,ng),NF_FOUT,         &
     &                   2,fgrd,Aval,Vinfo,ncname)
          Vinfo( 1)='y'
          Vinfo( 2)='y-location of floats trajectories'
          Vinfo( 3)='meter'
          Vinfo(14)='y, scalar, series'
          Vinfo(16)=Vname(1,idtime)
          Vinfo(17)='missing_value'
          Aval(4)=spval
          status=def_var(ncfltid(ng),fltVid(idglat,ng),NF_FOUT,         &
     &                   2,fgrd,Aval,Vinfo,ncname)
        END IF
# ifdef SOLVE3D
!
!  Define floats depths.
!
        Vinfo( 1)='depth'
        Vinfo( 2)='depth of floats trajectories'
        Vinfo( 3)='meter'
        Vinfo(14)='depth, scalar, series'
        Vinfo(16)=Vname(1,idtime)
        Vinfo(17)='missing_value'
        Aval(4)=spval
        status=def_var(ncfltid(ng),fltVid(iddpth,ng),NF_FOUT,           &
     &                 2,fgrd,Aval,Vinfo,ncname)
!
!  Define density anomaly.
!
        Vinfo( 1)=Vname(1,idDano)
        Vinfo( 2)=Vname(2,idDano)
        Vinfo( 3)=Vname(3,idDano)
        Vinfo(14)=Vname(4,idDano)
        Vinfo(16)=Vname(1,idtime)
        Vinfo(17)='missing_value'
        Aval(4)=spval
        status=def_var(ncfltid(ng),fltVid(idDano,ng),NF_FOUT,           &
     &                 2,fgrd,Aval,Vinfo,ncname)
!
!  Define tracer type variables.
!
        DO itrc=1,NT(ng)
          Vinfo( 1)=Vname(1,idTvar(itrc))
          Vinfo( 2)=Vname(2,idTvar(itrc))
          Vinfo( 3)=Vname(3,idTvar(itrc))
          Vinfo(14)=Vname(4,idTvar(itrc))
          Vinfo(16)=Vname(1,idtime)
          Vinfo(17)='missing_value'
          Aval(4)=spval
#  ifdef SEDIMENT
          DO i=1,NST
            IF (itrc.eq.idsed(i)) THEN
              WRITE (Vinfo(19),20) 1000.0_r8*Sd50(i,ng)
            END IF
          END DO
#  endif
          status=def_var(ncfltid(ng),fltTid(itrc,ng),NF_FOUT,           &
     &                   2,fgrd,Aval,Vinfo,ncname)
        END DO
# endif
!
!  Initialize unlimited time record dimension.
!
        tfltindx(ng)=0
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
        status=nf_enddef(ncfltid(ng))
      END IF
!
!-----------------------------------------------------------------------
!  Write out time-recordless, information variables.
!-----------------------------------------------------------------------
!
      IF (LdefFLT(ng)) THEN
        status=wrt_info(ng,ncfltid(ng),Master,ncname)
        if (exit_flag.ne.0) RETURN
      END IF
!
!=======================================================================
!  Open an existing floats file, check its contents, and prepare for
!  appending data.
!=======================================================================
!
      IF (.not.LdefFLT(ng).and.Master) THEN
!
!  Inquire about the contents of floats NetCDF file:  Inquire about
!  the dimensions and variables.  Check for consistency.
!
        ncname=fltname(ng)
        ifield=0
        CALL opencdf (ng,1,ncname,fname,N(ng),ifield,nrec)
        if (exit_flag.ne.0) RETURN
!
!  Open floats file for read/write.
!
        status=nf_open(TRIM(ncname),nf_write,ncfltid(ng))
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,30) TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
!
!  Initialize logical switches.
!
        DO i=1,NV
          got_var(i)=.false.
        END DO
!
!  Inquire size of the "drifter" dimension.
!
        status=nf_inq_dimid(ncfltid(ng),'drifter',fltdim)
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,40) 'drifter', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        status=nf_inq_dim(ncfltid(ng),fltdim,dimnam,Nfloats(ng))
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,50) 'drifter', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
!
!  Scan variable list from input NetCDF and activate switches for
!  float variables. Get variable IDs.
!
        DO i=1,nvars
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idtime))) THEN
            got_var(idtime)=.true.
            status=nf_inq_varid(ncfltid(ng),TRIM(Vname(1,idtime)),      &
     &                          fltVid(idtime,ng))
          END IF
          IF (TRIM(varnam(i)).eq.'Xgrid') THEN
            got_var(idXgrd)=.true.
            status=nf_inq_varid(ncfltid(ng),'Xgrid',fltVid(idXgrd,ng))
          END IF
          IF (TRIM(varnam(i)).eq.'Ygrid') THEN
            got_var(idYgrd)=.true.
            status=nf_inq_varid(ncfltid(ng),'Ygrid',fltVid(idYgrd,ng))
          END IF
# ifdef SOLVE3D
          IF (TRIM(varnam(i)).eq.'Zgrid') THEN
            got_var(idZgrd)=.true.
            status=nf_inq_varid(ncfltid(ng),'Zgrid',fltVid(idZgrd,ng))
          END IF
# endif
          IF (spherical) THEN
            IF (TRIM(varnam(i)).eq.'lon') THEN
              got_var(idglon)=.true.
              status=nf_inq_varid(ncfltid(ng),'lon',fltVid(idglon,ng))
            END IF
            IF (TRIM(varnam(i)).eq.'lat') THEN
              got_var(idglat)=.true.
              status=nf_inq_varid(ncfltid(ng),'lat',fltVid(idglat,ng))
            END IF
          ELSE
            IF (TRIM(varnam(i)).eq.'x') THEN
              got_var(idglon)=.true.
              status=nf_inq_varid(ncfltid(ng),'x',fltVid(idglon,ng))
            END IF
            IF (TRIM(varnam(i)).eq.'y') THEN
              got_var(idglat)=.true.
              status=nf_inq_varid(ncfltid(ng),'y',fltVid(idglat,ng))
            END IF
          END IF
# ifdef SOLVE3D
          IF (TRIM(varnam(i)).eq.'depth') THEN
            got_var(iddpth)=.true.
            status=nf_inq_varid(ncfltid(ng),'depth',fltVid(iddpth,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idDano))) THEN
            got_var(idDano)=.true.
            status=nf_inq_varid(ncfltid(ng),TRIM(Vname(1,idDano)),      &
     &                          fltVid(idDano,ng))
          END IF
          DO itrc=1,NT(ng)
            IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idTvar(itrc)))) THEN
             got_var(idTvar(itrc))=.true.
             status=nf_inq_varid(ncfltid(ng),                           &
     &                           TRIM(Vname(1,idTvar(itrc))),           &
     &                           fltTid(itrc,ng))
            END IF
          END DO
# endif
        END DO
!
!  Check if floats variables are available in input NetCDF file.
!
        IF (.not.got_var(idtime)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idtime)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idXgrd)) THEN
          WRITE (stdout,60) 'Xgrid', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idYgrd)) THEN
          WRITE (stdout,60) 'Ygrid', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# ifdef SOLVE3D
        IF (.not.got_var(idZgrd)) THEN
          WRITE (stdout,60) 'Zgrid', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
        IF (.not.got_var(idglon)) THEN
          IF (spherical) THEN
            WRITE (stdout,60) 'lon', TRIM(ncname)
          ELSE
            WRITE (stdout,60) 'x', TRIM(ncname)
          END IF
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idglat)) THEN
          IF (spherical) THEN
            WRITE (stdout,60) 'lat', TRIM(ncname)
          ELSE
            WRITE (stdout,60) 'y', TRIM(ncname)
          END IF
          exit_flag=3
          RETURN
        END IF
# ifdef SOLVE3D
        IF (.not.got_var(iddpth)) THEN
          WRITE (stdout,60) 'depth', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idDano)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idDano)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        DO itrc=1,NT(ng)
          IF (.not.got_var(idTvar(itrc))) THEN
            WRITE (stdout,60) TRIM(Vname(1,idTvar(itrc))), TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
# endif
!
!-----------------------------------------------------------------------
!  Initialize floats positions to the appropriate values.
!-----------------------------------------------------------------------
!
!  Set-up floats time record.
!
        IF (frrec(ng).lt.0) THEN
          tfltindx(ng)=tsize
        ELSE
         tfltindx(ng)=ABS(frrec(ng))
!!       tfltindx(ng)=1+(ntstart-1)/nflt(ng)
        END IF
        nrecflt=tfltindx(ng)
!
!  Read in floats nondimentional horizontal positions.
!
        start(1)=1
        total(1)=Nfloats(ng)
        start(2)=tfltindx(ng)
        total(2)=1
        status=nf_get_vara_TYPE(ncfltid(ng),fltVid(idXgrd,ng),          &
     &                          start,total,Tinp)
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,70) 'Xgrid', tfltindx(ng), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        DO l=1,Nfloats(ng)
          IF ((Tinp(l).gt.REAL(Lm(ng)+1,r8)-0.5_r8).or.                 &
     &        (Tinp(l).lt.0.5_r8)) THEN
            FLT(ng)%bounded(l)=.false.
          ELSE
            FLT(ng)%bounded(l)=.true.
            DO i=0,NFT
              FLT(ng)%track(ixgrd,i,l)=Tinp(l)
              FLT(ng)%track(ixrhs,i,l)=0.0_r8
            END DO
          END IF
        END DO
        status=nf_get_vara_TYPE(ncfltid(ng),fltVid(idYgrd,ng),          &
     &                          start,total,Tinp)
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,70) 'Ygrid', tfltindx(ng), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        DO l=1,Nfloats(ng)
          IF ((Tinp(l).gt.REAL(Mm(ng)+1,r8)-0.5_r8).or.                 &
     &        (Tinp(l).lt.0.5_r8)) THEN
            FLT(ng)%bounded(l)=.false.
          ELSE
            FLT(ng)%bounded(l)=.true.
            DO i=0,NFT
              FLT(ng)%track(iygrd,i,l)=Tinp(l)
              FLT(ng)%track(iyrhs,i,l)=0.0_r8
            END DO
          END IF
        END DO
# ifdef SOLVE3D
        status=nf_get_vara_TYPE(ncfltid(ng),fltVid(idZgrd,ng),          &
     &                          start,total,Tinp)
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,70) 'Zgrid', tfltindx(ng), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        DO l=1,Nfloats(ng)
          IF ((Tinp(l).gt.REAL(N(ng),r8)).or.(Tinp(l).lt.0.0_r8)) THEN
            FLT(ng)%bounded(l)=.false.
          ELSE
            FLT(ng)%bounded(l)=.true.
            DO i=0,NFT
              FLT(ng)%track(izgrd,i,l)=Tinp(l)
              FLT(ng)%track(izrhs,i,l)=0.0_r8
            END DO
          END IF
        END DO
# endif
      END IF

# ifdef DISTRIBUTE
!
!  Broadcast restart float information to all nodes.
!
      IF (.not.LdefFLT(ng)) THEN

        Npts=NFV(ng)*(NFT+1)*Nfloats(ng)

        CALL mp_collect (ng, Npts, Fspv, FLT(ng)%track)
!
!  Collect the bounded status switch.
!
        Fwrk=Fspv
        DO l=1,Nfloats(ng)
          IF (FLT(ng)%bounded(l)) THEN
            Fwrk(l)=1.0_r8
          END IF
        END DO
        CALL mp_collect (ng, Nfloats(ng), Fspv, Fwrk)
        DO l=1,Nfloats(ng)
          IF (Fwrk(l).ne.Fspv) THEN
            FLT(ng)%bounded(l)=.true.
          ELSE
            FLT(ng)%bounded(l)=.false.
          END IF
        END DO
      END IF
# endif
!
  10  FORMAT (/,' DEF_FLOATS - unable to create floats NetCDF',         &
     &          ' file: ',a)
  20  FORMAT (1pe11.4,1x,'millimeter')
  30  FORMAT (/,' DEF_FLOATS - unable to open floats NetCDF file: ',a)
  40  FORMAT (/,' DEF_FLOATS - error while inquiring dimension ID ',    &
     &          ' for: ',a,2x,/,13x,'in floats NetCDF file: ',a)
  50  FORMAT (/,' DEF_FLOATS - error while inquiring size of',          &
     &        ' dimension: ',a,2x,/,13x,'in floats NetCDF file: ',a)
  60  FORMAT (/,' DEF_FLOATS - unable to find variable: ',a,2x,         &
     &          ' in floats NetCDF file: ',a)
  70  FORMAT (/,' DEF_FLOATS - error while reading variable: ',a,2x,    &
     &          ' at time record = ', i6.6,/,13x,                       &
     &          ' in floats NetCDF file: ',a)
      RETURN
      END SUBROUTINE def_floats
#else
      SUBROUTINE def_floats
      RETURN
      END SUBROUTINE def_floats
#endif
