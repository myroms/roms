#include "cppdefs.h"
      SUBROUTINE opencdf (ng, nfiles, fname, ncname, nvlev, ifield,     &
     &                    nrec)
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine opens an existing NetCDF file and inquires about       !
!  it contents, and checks for consistency with model dimensions.      !
!                                                                      !
!  Input:                                                              !
!                                                                      !
!    ng        Nested grid number.                                     !
!    nfiles    Number of input NetCDF files.                           !
!    fname     Input NetCDF file name(s).                              !
!    ncname    NetCDF file name containing requested variable.         !
!    nvlev     Number of vertical levels to check for consistency.     !
!    ifield    Index of field to inquire for the size of its time      !
!                dimension, if any.                                    !
!                                                                      !
!  Output:                                                             !
!                                                                      !
!     nrec     Number of time records available for requested field.   !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_scalars
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, nfiles, nvlev, ifield

      integer, intent(out) :: nrec

      character (len=*), intent(in) :: fname(nfiles)

      character (len=80), intent(out) :: ncname
!
!  Local variable declarations.
!
      logical :: SearchVar, timeatt

      integer :: attype, dimid, dimsiz, i, ifile, j
      integer :: ltvar, ltype, ncid, ndims, ngatts
      integer :: nvatts, recdim, status, varid, vartype

      character (len=20) :: dimnam, attnam
      character (len=45) :: text
!
!-----------------------------------------------------------------------
!  Open input NetCDF file(s). If several input NetCDF files (nfiles>1),
!  scan files until the requested variable is found.
!-----------------------------------------------------------------------
!
      SearchVar=.true.
      DO ifile=1,nfiles
        IF (SearchVar) THEN
          ncname=fname(ifile)
          status=nf_open (TRIM(ncname),nf_nowrite,ncid)
          IF (status.ne.nf_noerr) THEN
            WRITE (stdout,10) TRIM(ncname)
            exit_flag=4
            RETURN
          END IF
!
!-----------------------------------------------------------------------
!  Inquire and get global "type" attribute.
!-----------------------------------------------------------------------
!
          status=nf_inq_att (ncid,nf_global,'type',attype,ltype)
          IF (status.eq.nf_noerr) THEN
            status=nf_get_att_text (ncid,nf_global,'type',type)
            IF (status.ne.nf_noerr) THEN
              WRITE (stdout,20) 'type (global)', TRIM(ncname)
              exit_flag=4
              RETURN
            END IF
          ELSE
            WRITE (stdout,30) 'type (global)', TRIM(ncname)
          END IF
!
!-----------------------------------------------------------------------
!  Inquire about the dimensions and variables.
!-----------------------------------------------------------------------
!
          recdim=-1
          status=nf_inq(ncid,ndims,nvars,ngatts,recdim)
          IF (status.eq.nf_noerr) THEN
            IF (nvars.gt.MV) THEN
              WRITE (stdout,40) MV, nvars
              exit_flag=4
              RETURN
            END IF
!
!  Inquire about dimensions. Check dimensions for consistency.
!
            nrec=0
            ltvar=0
            IF (ifield.gt.0) THEN
              ltvar=LEN_TRIM(Vname(5,ifield))
              Tname(ifield)=TRIM(Vname(5,ifield))
            END IF
            DO i=1,ndims
              dimid=i
              status=nf_inq_dim(ncid,dimid,dimnam,dimsiz)
              IF (status.ne.nf_noerr) THEN
                WRITE (stdout,50) dimid, TRIM(ncname)
                exit_flag=4
                RETURN
              END IF
              IF ((TRIM(dimnam).eq.'xi_rho').or.                        &
     &            (TRIM(dimnam).eq.'xi_v')) THEN
                IF (dimsiz.ne.Lm(ng)+2) THEN
                  WRITE (stdout,60) TRIM(dimnam), dimsiz, Lm(ng)+2
                  exit_flag=4
                  RETURN
                END IF
              ELSE IF ((TRIM(dimnam).eq.'xi_u').or.                     &
     &                 (TRIM(dimnam).eq.'xi_psi')) THEN
                IF (dimsiz.ne.Lm(ng)+1) THEN
                  WRITE (stdout,60) TRIM(dimnam), dimsiz, Lm(ng)+1
                  exit_flag=4
                  RETURN
                END IF
              ELSE IF ((TRIM(dimnam).eq.'eta_rho').or.                  &
     &                 (TRIM(dimnam).eq.'eta_u')) THEN
                IF (dimsiz.ne.Mm(ng)+2) THEN
                  WRITE (stdout,60) TRIM(dimnam), dimsiz, Mm(ng)+2
                  exit_flag=4
                  RETURN
                END IF
              ELSE IF ((TRIM(dimnam).eq.'eta_v').or.                    &
     &                 (TRIM(dimnam).eq.'eta_psi')) THEN
                IF (dimsiz.ne.Mm(ng)+1) THEN
                  WRITE (stdout,60) TRIM(dimnam), dimsiz, Mm(ng)+1
                  exit_flag=4
                  RETURN
                END IF
#ifdef SOLVE3D
              ELSE IF (TRIM(dimnam).eq.'s_rho') THEN
                IF (dimsiz.ne.nvlev) THEN
                  WRITE (stdout,60) TRIM(dimnam), dimsiz, nvlev
                  exit_flag=4
                  RETURN
                END IF
#endif
              ELSE IF ((ltvar.gt.0).and.                                &
     &                 (TRIM(dimnam).eq.TRIM(Tname(ifield)))) THEN
                nrec=dimsiz
              END IF
            END DO
!
!  Inquire about variables.  If "nrec" is zero, indicating that the
!  requested variable does not have the specified time dimension;
!  then, check the size of the dimension of its associated time
!  variable.  The associated time variable is reset to the value
!  specified in the "time" attribute, if any.  Check if only water
!  points are available.
!
            timeatt=.false.
            DO i=1,nvars
              varid=i
              vflag(i)=1
              status=nf_inq_var(ncid,varid,varnam(i),vartype,nvdims(i), &
     &                          vdims(1,i),nvatts)
              IF (status.ne.nf_noerr) THEN
                WRITE (stdout,70) varid, TRIM(ncname)
                exit_flag=4
                RETURN
              END IF
              IF (nvdims(i).gt.0) THEN
                DO j=1,nvatts
                  status=nf_inq_attname(ncid,varid,j,attnam)
                  IF (status.eq.nf_noerr) THEN
                    IF (TRIM(attnam).eq.'water_points') THEN
                      vflag(i)=-1
                    END IF
                  ELSE
                    WRITE (stdout,100) TRIM(varnam(i))
                    exit_flag=4
                    RETURN
                  END IF
                END DO
              END IF
              IF ((LEN_TRIM(Vname(1,ifield)).gt.0).and.                 &
     &            (TRIM(varnam(i)).eq.TRIM(Vname(1,ifield)))) THEN
                SearchVar=.false.
                DO j=1,nvatts
                  status=nf_inq_attname(ncid,varid,j,attnam)
                  IF (status.eq.nf_noerr) THEN
                    IF (TRIM(attnam).eq.'time') THEN
                      status=nf_inq_attlen(ncid,varid,TRIM(attnam),     &
     &                                     ltvar)
                      IF (status.ne.nf_noerr) THEN
                        WRITE (stdout,80) TRIM(attnam)
                        exit_flag=4
                        RETURN
                      END IF
                      status=nf_get_att_text(ncid,varid,TRIM(attnam),   &
     &                                       text(1:ltvar))
                      IF (status.ne.nf_noerr) THEN
                        WRITE (stdout,90) TRIM(attnam)
                        exit_flag=4
                        RETURN
                      END IF
                      Tname(ifield)=text(1:ltvar)
                      timeatt=.true.
                    END IF
                  ELSE
                    WRITE (stdout,100) TRIM(varnam(i))
                    exit_flag=4
                    RETURN
                  END IF
                END DO
              END IF
            END DO
            IF (timeatt) THEN
              ltvar=LEN_TRIM(Tname(ifield))
              IF ((ifield.gt.0).and.(nrec.eq.0).and.(ltvar.gt.0)) THEN
                DO i=1,nvars
                  varid=i
                  status=nf_inq_var (ncid,varid,varnam(i),vartype,      &
     &                               nvdims(i),vdims(1,i),nvatts)
                  IF (status.ne.nf_noerr) THEN
                    WRITE (stdout,70) varid, TRIM(ncname)
                    exit_flag=4
                    RETURN
                  END IF
                  IF (TRIM(varnam(i)).eq.TRIM(Tname(ifield))) THEN
                    DO j=1,nvdims(i)
                      status=nf_inq_dim (ncid,vdims(j,i),dimnam,dimsiz)
                      IF (status.ne.nf_noerr) THEN
                        WRITE (stdout,110) TRIM(Tname(ifield))
                        exit_flag=4
                        RETURN
                      END IF
                      IF (INDEX(TRIM(dimnam),'time').ne.0) THEN
                        nrec=dimsiz
                      END IF
                    END DO
                  END IF
                END DO
              END IF
            END IF
          ELSE
            WRITE (stdout,120) TRIM(ncname)
            exit_flag=4
            RETURN
          END IF
!
!-----------------------------------------------------------------------
!  Inquire size of unlimited time record dimension.
!-----------------------------------------------------------------------
!
          tsize=0
          IF (recdim.gt.0) THEN
            status=nf_inq_dim (ncid,recdim,dimnam,tsize)
            IF (status.ne.nf_noerr) THEN
              WRITE (stdout,130) 'time', TRIM(ncname)
              exit_flag=4
              RETURN
            END IF
          END IF
!
!-----------------------------------------------------------------------
!  Close input NetCDF file.
!-----------------------------------------------------------------------
!
          status=nf_close (ncid)
!
!-----------------------------------------------------------------------
!  If appropriate, scan next NetCDF file.
!-----------------------------------------------------------------------
!
        END IF
      END DO
!
  10  FORMAT (/,' OPENCDF - unable to open input NetCDF file: ',a)
  20  FORMAT (/,' OPENCDF - error while reading attribute: ',a,2x,      &
     &          ' in input NetCDF file: ',a)
  30  FORMAT (/,' OPENCDF - cannot find attribute: ',a,2x,              &
     &          ' in input NetCDF file: ',a)
  40  FORMAT (/,' OPENCDF - too small dimension parameter, MV = ',      &
     &       2i5,/,11x,'change file  ncparam.h  and recompile.')
  50  FORMAT (/,' OPENCDF - error while reading dimension ID: ',i3,2x,  &
     &          ' in input NetCDF file: ',a)
  60  FORMAT (/,' OPENCDF - inconsistent size of dimension: ',a,2x,     &
     &       2i5)
  70  FORMAT (/,' OPENCDF - error while inquiring inFORMAT ion for ',   &
     &          ' variable ID: ',i3,2x,' in input NetCDF file: ',a)
  80  FORMAT (/,' OPENCDF - error while inquiring length of',           &
     &          ' attribute: ',a)
  90  FORMAT (/,' OPENCDF - error while reading attribute: ',a)
 100  FORMAT (/,' OPENCDF - error while inquiring attributes for',      &
     &          ' variable: ',a)
 110  FORMAT (/,' OPENCDF - unable to inquire size of time dimesion',   &
     &          ' in variable: ',a)
 120  FORMAT (/,' OPENCDF - unable to inquire about contents of',       &
     &          ' input NetCDF file: ',a)
 130  FORMAT (/,' OPENCDF - error inquiring dimension: ',a,2x,          &
     &          ' in input NetCDF file: ',a)
      RETURN
      END SUBROUTINE opencdf

      SUBROUTINE openids (nfiles, fname, fid, ncname, ncid)
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine determines association between pathnames and NetCDF    !
!  file descriptors. It is used to get the  ID  of an opened NetCDF    !
!  from file name.  This is done to avoid  opening  too  many files    !
!  simultaneously for the same dataset.                                !
!                                                                      !
!  Input:                                                              !
!                                                                      !
!    nfiles    Number of input NetCDF files.                           !
!    fname     NetCDF file name(s).                                    !
!    fid       NetCDF ID(s) associated with fname.                     !
!    ncname    NetCDF file name to check in association table.         !
!    ncid      NetCDF ID to check in association table.                !
!                                                                      !
!  Output:                                                             !
!                                                                      !
!    fid       If updated, it stores the ID associated with fname.     !
!    ncid      If updated, indicates that ncname has been opened       !
!                previously and this is its associated ID.             !
!                                                                      !
!=======================================================================
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: nfiles

      character (len=*), intent(in) :: fname(nfiles)
      character (len=*), intent(in) :: ncname

      integer, intent(inout) :: fid(nfiles)
      integer, intent(inout) :: ncid
!
!  Local variabel declarations.
!
      integer :: i
!
!-----------------------------------------------------------------------
!  Determine association between NetCDF file name and ID descriptor.
!-----------------------------------------------------------------------
!
!  NetCDF file ID descriptors are set to closed state (-1) during
!  initialization.
!
      DO i=1,nfiles
        IF (TRIM(ncname).eq.TRIM(fname(i))) THEN
!
!  NCFILE was just opened. Store its associated ID for future reference.
!
          IF (ncid.ne.-1) fid(i)=ncid
!
!  NCFILE is already open. Get its ID from reference table.
!
          IF ((fid(i).ne.-1).and.(ncid.eq.-1)) ncid=fid(i)
        END IF
      END DO

      RETURN
      END SUBROUTINE openids
