#include "cppdefs.h"
      SUBROUTINE get_ngfld (ng, ifield, ncid, nfiles, fname, update,    &
     &                      LBi, UBi, UBj, UBk, Istr, Iend, Jrec,       &
     &                      Fout)
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine reads in requested non-grided field from specified     !
!  NetCDF file.  A non-grided field has different dimensions  than     !
!  model spatial dimensions.                                           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     ifield     Field ID.                                             !
!     ncid       NetCDF file ID.                                       !
!     nfiles     Number of input NetCDF files.                         !
!     ncfile     NetCDF file name(s).                                  !
!     LBi        "Fout" 1st dimension lower-bound value.               !
!     UBi        "Fout" 1st dimension upper-bound value.               !
!     UBj        "Fout" 2nd dimension upper-bound value, if any.       !
!                  Otherwise, a value of one is expected.              !
!     UBk        "Fout" time dimension upper-bound value, if any.      !
!                  Otherwise, a value of one is expected.              !
!     Istr       Starting location of read data in the 1st dimension.  !
!     Iend       Ending location of read data in the 1st dimension;    !
!                  Number of records read is: Iend-Istr+1.             !
!     Jrec       Number of records read in the 2st dimension.          !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Fout       Read field.                                           !
!     update     Switch indicating reading of the requested field      !
!                  the current time step.                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars

      implicit none
!
!  Imported variable declarations.
!
      logical, intent(out) :: update

      integer, intent(in) :: ng, ifield, ncid, nfiles
      integer, intent(in) :: LBi, UBi, UBj, UBk, Istr, Jrec
      integer, intent(inout) :: Iend

      character (len=*), intent(in) :: fname(nfiles)

      real(r8), intent(inout) :: Fout(LBi:UBi,UBj,UBk)
!
!  Local variable declarations.
!
      logical :: Lgrided, Liocycle, Lonerec, got_var, got_time

      integer :: Nrec, Tid, Tindex, Trec, Vid, Vtype
      integer :: i, ic, j, lfvar, ltvar, npts, nvatts, status, vartype

      integer, dimension(3) :: start, total

      real(r8) :: Aval, Clength, Fmax, Fmin, Tdelta, Tend, Tmono, Tscale
      real(r8) :: Tstart, Tstr, Tval, scale

      real(r8), dimension((UBi-LBi+1)*UBj) :: A

      character (len=45) :: attnam
      character (len=80) :: ncfile
!
!-----------------------------------------------------------------------
!  On first call, inquire about the contents of input NetCDF file.
!-----------------------------------------------------------------------
!
      IF (exit_flag.ne.0) RETURN
!
      IF (iic(ng).eq.0) THEN
!
!  Intialize local variables.
!
        Vid=-1
        Tid=-1        
        Liocycle=.false.
        Lgrided=.false.
        Lonerec=.false.
        got_var=.false.
        got_time=.false.
        Vtype=Iinfo(1,ifield,ng)
        lfvar=LEN_TRIM(Vname(1,ifield))
        ltvar=LEN_TRIM(Tname(ifield))
!
!  Inquire about the dimensions and variables. Check for consistency.
!
        IF (InpThread) THEN
          CALL opencdf (ng, nfiles, fname, ncfile, N(ng), ifield, Nrec)
          IF (exit_flag.ne.0) RETURN
          IF ((UBk.eq.1).and.(Nrec.gt.UBi)) THEN
            WRITE (stdout,10) TRIM(Vname(1,ifield)), UBi, Nrec
            exit_flag=4
            RETURN
          END IF
        END IF
        IF (Iend.eq.0) THEN
          Iend=Nrec
#ifdef DISTRIBUTE
          CALL mp_bcasti (ng, Iend, 1)
#endif
        END IF
!
!  Scan variable list from input NetCDF and check for requested
!  variables.
!
        IF (InpThread) THEN
          DO i=1,nvars
            IF (lfvar.gt.0) THEN
              IF (TRIM(varnam(i)).eq.TRIM(Vname(1,ifield))) THEN
                Vid=i
                got_var=.true.
              END IF
            END IF
            IF (ltvar.gt.0) THEN
              IF (TRIM(varnam(i)).eq.TRIM(Tname(ifield))) THEN
                Tid=i
                got_time=.true.
              END IF
            END IF
          END DO
          Linfo(1,ifield,ng)=Lgrided
          Iinfo(2,ifield,ng)=Vid
          Iinfo(3,ifield,ng)=Tid
          Iinfo(4,ifield,ng)=Nrec
!
!  Terminate execution requested variables are not found.
!
          IF (.not.got_var) THEN
            WRITE (stdout,20) TRIM(Vname(1,ifield)), TRIM(ncfile)
            exit_flag=2
            RETURN
          END IF
          IF (.not.got_time) THEN
            WRITE (stdout,20) TRIM(Tname(ifield)), TRIM(ncfile)
            exit_flag=2
            RETURN
          END IF
!
!  If appropriate, open input NetCDF file for reading.  If processing
!  model forcing (multiple files allowed), check if file for requested
!  field has been already opened and get/save its ID from/to the
!  association table. 
!
          CALL openids (nFfiles(ng), frcname(1,ng), frcids(1,ng),       &
     &                  ncfile, ncid)
          IF (ncid.eq.-1) THEN
            status=nf_open(TRIM(ncfile), nf_nowrite, ncid)
            IF (status.ne.nf_noerr) THEN
              WRITE (stdout,20) TRIM(ncfile)
              exit_flag=2
              RETURN
            END IF
            CALL openids (nFfiles(ng), frcname(1,ng), frcids(1,ng),     &
     &                    ncfile, ncid)
          END IF
!
!  If "scale_factor" attribute is present for a variable, the data are
!  to be multiplied by this factor.  Check if only water points are
!  available.
!
          IF (got_var) THEN
            status=nf_inq_var(ncid, Vid, varnam(Vid), vartype,          &
     &                        nvdims(Vid), vdims(1,Vid), nvatts)
            IF (status.eq.nf_noerr) THEN
              DO i=1,nvatts
                status=nf_inq_attname(ncid, Vid, i, attnam)
                IF (status.eq.nf_noerr) THEN
                  IF (TRIM(attnam).eq.'scale_factor') THEN
                    status=nf_get_att_TYPE(ncid, Vid, TRIM(attnam),     &
     &                                     scale)
                    IF (status.eq.nf_noerr) THEN
                      Fscale(ifield,ng)=Fscale(ifield,ng)*scale
                    ELSE
                      WRITE (stdout,30) TRIM(attnam)
                      exit_flag=2
                      RETURN
                    END IF
                  ELSE IF (TRIM(attnam).eq.'water_points') THEN
                    Iinfo(1,ifield,ng)=-ABS(Iinfo(1,ifield,ng))
                    Vtype=Iinfo(1,ifield,ng)
                  END IF
                ELSE
                  WRITE (stdout,40) TRIM(varnam(Vid))
                  exit_flag=2
                  RETURN
                END IF
              END DO
            ELSE
              WRITE (stdout,50)
              exit_flag=2
              RETURN
            END IF
          END IF
!
!  Determine initial time record to read and cycling switch.
!
          IF (UBk.eq.1) THEN
            Tstart=1
            Tend=1
            Trec=1
            Clength=0.0_r8
            Tscale=1.0_r8
          ELSE
            CALL get_cycle (ng, ncid, Tid, Nrec, tdays(ng), Liocycle,   &
     &                      Clength, Trec, Tstart, Tstr, Tend, Tscale)
            IF (exit_flag.ne.0) RETURN
          END IF
          Linfo(2,ifield,ng)=Liocycle
          Finfo(1,ifield,ng)=Tstart
          Finfo(2,ifield,ng)=Tend
          Finfo(3,ifield,ng)=Clength
          Finfo(4,ifield,ng)=Tscale
!
!  The strategy here is to create a local, monotonically increasing
!  time variable so the interpolation between snapshots is trivial
!  when cycling forcing fields. Subtract one to time record counter
!  "Trec" to avoid doing special case at initialization.
!
          IF (Iend.eq.1) THEN
            Tindex=UBk
          ELSE
            Tindex=1
          END IF
          IF (Liocycle) THEN
            IF (Trec.eq.Nrec) THEN
              IF (tdays(ng).lt.Tend) THEN
                Tmono=Tstart-Clength
              ELSE
                Tmono=tdays(ng)+(Tstart-Clength)
                IF (Tstart.eq.Tend) THEN
                  Tmono=Tmono+(Tstr-MOD(tdays(ng)+Tstr,Clength))
                ELSE
                  Tmono=Tmono+(Tstart-MOD(tdays(ng)+Tstart,Clength))
                END IF
              END IF
              Tmono=Tmono*day2sec
            ELSE
              IF (tdays(ng).gt.Clength) THEN
                Tmono=time(ng)-MOD(tdays(ng)-Tstart,Clength)*day2sec
              ELSE
                Tmono=Tstart*day2sec
              END IF
            END IF
          ELSE
            Tmono=time(ng)
          END IF
          Trec=Trec-1
          Iinfo(5,ifield,ng)=Tindex
          Iinfo(6,ifield,ng)=Trec
          Finfo(5,ifield,ng)=Tmono
!
!  Set switch for one time record dataset. In this case, the input
!  data is always the same and time interpolation is not performed.
!
          IF (Nrec.eq.1) Lonerec=.true.
          Linfo(3,ifield,ng)=Lonerec
        END IF
#ifdef DISTRIBUTE
!
!  Broadcast variable information to all nodes in the group.
!
        CALL mp_bcastl (ng, Linfo(1,ifield,ng), 3)
        CALL mp_bcasti (ng, Iinfo(1,ifield,ng), 6)
        CALL mp_bcastf (ng, Finfo(1,ifield,ng), 5)
#endif
        Tindex=Iinfo(5,ifield,ng)
        Vtime(Tindex,ifield,ng)=Finfo(1,ifield,ng)
      END IF         
!
!-----------------------------------------------------------------------
!  Get requested field information from global storage.
!-----------------------------------------------------------------------
!
      Liocycle=Linfo(2,ifield,ng)
      Lonerec =Linfo(3,ifield,ng)
      Vtype   =Iinfo(1,ifield,ng)
      Vid     =Iinfo(2,ifield,ng)
      Tid     =Iinfo(3,ifield,ng)
      Nrec    =Iinfo(4,ifield,ng)
      Tindex  =Iinfo(5,ifield,ng)
      Trec    =Iinfo(6,ifield,ng)
      Clength =Finfo(3,ifield,ng)
      Tscale  =Finfo(4,ifield,ng)
      Tmono   =Finfo(5,ifield,ng)
!
!-----------------------------------------------------------------------
!  If appropriate, read in new data.
!-----------------------------------------------------------------------
!
      update=.false.
      IF ((Tmono.lt.time(ng)).or.(iic(ng).eq.0).or.                     &
     &    (iic(ng).eq.ntstart)) THEN
        IF (Liocycle) THEN
          Trec=MOD(Trec,Nrec)+1
        ELSE
          Trec=Trec+1
        END IF
        Iinfo(6,ifield,ng)=Trec
        IF (Trec.le.Nrec) THEN
!
!  Set rolling index for two-time record storage of input data.  If
!  "UBk" is unity, input data is stored in recordless array by the
!  calling program.
!
         IF (UBk.eq.1) THEN
           Tindex=1
         ELSE
           Tindex=3-Tindex
         END IF
         Iinfo(5,ifield,ng)=Tindex
!
!  Read in time coordinate.
!
          IF (InpThread.and.((Tid.ge.0).and.(Tid.ne.Vid))) THEN
            status=nf_get_var1_TYPE(ncid, Tid, Trec, Tval)
            IF (InpThread.and.(status.ne.nf_noerr)) THEN
              WRITE (stdout,40) TRIM(Tname(ifield)), Trec
              exit_flag=2
              RETURN
            END IF
          END IF
          IF ((Tid.ge.0).and.(Tid.ne.Vid)) THEN
#ifdef DISTRIBUTE
            CALL mp_bcastf (ng, Tval, 1)
#endif
            Vtime(Tindex,ifield,ng)=Tval*Tscale
          END IF
!
!  Read in non-grided data.
!
          IF (Vid.ge.0) THEN
            Fmin=0.0_r8
            Fmax=0.0_r8
            IF (Jrec.gt.1) THEN
              start(1)=1
              total(1)=(Iend-Istr)+1
              start(2)=1
              total(2)=Jrec
              start(3)=Trec
              total(3)=1
              npts=total(1)*total(2)*total(3)
            ELSE
              start(1)=1
              total(1)=(Iend-Istr)+1
              start(2)=Trec
              total(2)=1
              npts=total(1)*total(2)
            END IF
            IF (InpThread) THEN
              status=nf_get_vara_TYPE(ncid, Vid, start, total, A)
              IF (status.ne.nf_noerr) THEN
                WRITE (stdout,40) TRIM(Vname(1,ifield)), Trec
                exit_flag=2
                RETURN
              END IF
            END IF
#ifdef DISTRIBUTE
            CALL mp_bcastf (ng, A, npts)
#endif
            Fmin=A(1)*Fscale(ifield,ng)
            Fmax=A(1)*Fscale(ifield,ng)
            ic=0
            DO j=1,Jrec
              DO i=Istr,Iend
                ic=ic+1
                Aval=A(ic)*Fscale(ifield,ng)
                Fmin=MIN(Fmin,Aval)
                Fmax=MAX(Fmax,Aval)
                Fout(i,j,Tindex)=Aval
              END DO
            END DO
            Finfo(6,ifield,ng)=Fmin
            Finfo(7,ifield,ng)=Fmax
            IF (Master) THEN
              IF (UBk.eq.1) THEN
                WRITE (stdout,50) TRIM(Vname(2,ifield)), Fmin, Fmax
              ELSE
                WRITE (stdout,60) TRIM(Vname(2,ifield)), Tval*Tscale,   &
     &                            Fmin, Fmax
              END IF
            END IF
            update=.true.
          END IF
        END IF
!
!  Increment the local time variable "Tmono" by the interval between
!  snapshots. If the interval is negative, indicating cycling, add in
!  a cycle length.  Load time value (sec) into "Tintrp" which used
!  during interpolation between snapshots.
!
        IF (.not.Lonerec) THEN
          Tdelta=Vtime(Tindex,ifield,ng)-Vtime(3-Tindex,ifield,ng)
          IF (Liocycle.and.(Tdelta.lt.0.0_r8)) THEN
            Tdelta=Tdelta+Clength
          END IF
          Tmono=Tmono+Tdelta*day2sec
          Finfo(5,ifield,ng)=Tmono
          Tintrp(Tindex,ifield,ng)=Tmono
        END IF
      END IF
!
  10  FORMAT (/,' GET_NGFLD   - too small dimension for variable ',a,   &
     &        ': ',2i5)
  20  FORMAT (/,' GET_NGFLD   - unable to find requested variable: ',a, &
     &       /,15x,'in input NetCDF file: ',a)
  30  FORMAT (/,' GET_NGFLD   - unable to open input NetCDF file: ',a)
  40  FORMAT (/,' GET_NGFLD   - error while reading variable: ',a,2x,   &
     &          ' at TIME index = ',i4)
  50  FORMAT (' GET_NGFLD   - ',a,/,17x,'(Min = ',1pe15.8,' Max = ',    &
     &        1pe15.8,')')
  60  FORMAT (' GET_NGFLD   - ',a,',',t64,'t = ',f12.4,/,17x,'(Min = ', &
     &        1pe15.8,' Max = ',1pe15.8,')')
      RETURN
      END SUBROUTINE get_ngfld
