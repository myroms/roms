#include "cppdefs.h"
      SUBROUTINE get_state (ng, model, msg, ncname, IniRec, Tindex)
!
!=======================================================================
!  Copyright (c) 2005  ROMS/TOMS Adjoint Group                         !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine reads in state variables nonlinear initial conditions  !
!  from requested NetCDF file.                                         !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     msg        Message index for Mstate.                             !
!     ncname     Nonlinear initial conditions NetCDF file name.        !
!     IniRec     Nonlinear initial conditions time record to read.     !
!     Tindex     State variable time index to intialize.               !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_ocean
      USE mod_scalars
#if defined SEDIMENT || defined BBL_MODEL
      USE mod_sediment
#endif
      USE mod_strings

#if (defined EW_PERIODIC || defined NS_PERIODIC) && !defined DISTRIBUTE
!
      USE exchange_2d_mod
# ifdef SOLVE3D
      USE exchange_3d_mod
# endif
#endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, msg, Tindex

      integer, intent(inout) :: IniRec

      character (len=*), intent(in) :: ncname
!
!  Local variable declarations.
!
      logical, dimension(NV) :: get_var, have_var

      integer :: LBi, UBi, LBj, UBj
      integer :: InpRec, gtype, i, ifield, itrc, lstr, lend
      integer :: ncINPid, nvatts, nrec, status, varid, vartype

      integer :: nf_fread

      real(r8), parameter :: Fscl = 1.0_r8

      real(r8) :: Fmax, Fmin, Tmax, time_scale

      character (len=15) :: attnam, tvarnam
      character (len=40) :: tunits
      character (len=80) :: fname
!
!-----------------------------------------------------------------------
!  Inquire about the contents of input NetCDF file:  Inquire about
!  the dimensions and variables.  Check for consistency.
!-----------------------------------------------------------------------
!
      IF (exit_flag.ne.NoError) RETURN
      ifield=0
      CALL opencdf (ng, 1, ncname, fname, N(ng), ifield, nrec)
      IF (exit_flag.ne.NoError) RETURN
!
!  Initialize logical switches.
!
        DO i=1,NV
          get_var(i)=.FALSE.
          have_var(i)=.FALSE.
        END DO
!
!  Determine state variables to read from input NetCDF file.
!
#ifndef ANA_INITIAL
        get_var(idFsur)=.TRUE.
        get_var(idUbar)=.TRUE.
        get_var(idVbar)=.TRUE.
# ifdef SOLVE3D
        get_var(idUvel)=.TRUE.
        get_var(idVvel)=.TRUE.
        DO itrc=1,NAT
          get_var(idTvar(itrc))=.TRUE.
        END DO
# endif
#endif
#ifdef SOLVE3D
# if defined T_PASSIVE && !defined ANA_PASSIVE
        DO itrc=1,NPT
          get_var(inert(itrc))=.TRUE.
        END DO
# endif
# ifdef SEDIMENT
#  ifndef ANA_SEDIMENT
        DO itrc=1,NST
          get_var(idTvar(idsed(itrc)))=.TRUE.
          get_var(idfrac(itrc))=.TRUE.
          get_var(idBmas(itrc))=.TRUE.
        END DO
        DO itrc=1,MBEDP
          get_var(idSbed(itrc))=.TRUE.
        END DO
#   ifdef BEDLOAD
        DO itrc=1,NST
          IF (IniRec.ne.0) THEN
            get_var(idUbld(itrc))=.TRUE.
            get_var(idVbld(itrc))=.TRUE.
          END IF
        END DO
#   endif
#  endif
# endif 
# if defined SEDIMENT || defined BBL_MODEL
        DO itrc=1,MBOTP
          get_var(idBott(itrc))=.TRUE.
        END DO
# endif
# if defined BIOLOGY && !defined ANA_BIOLOGY
        DO itrc=1,NBT
          get_var(idbio(itrc))=.TRUE.
        END DO
# endif
#endif
!
!  Scan variable list from input NetCDF and activate switches for
!  model initialization variables.
!
      DO i=1,nvars
        IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idtime))) THEN
          tvarnam=TRIM(varnam(i))
          have_var(idtime)=.TRUE.
        END IF
        IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idFsur))) THEN
          have_var(idFsur)=.TRUE.
        END IF
        IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idUbar))) THEN
          have_var(idUbar)=.TRUE.
        END IF
        IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVbar))) THEN
          have_var(idVbar)=.TRUE.
        END IF
#ifdef SOLVE3D
        IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idUvel))) THEN
          have_var(idUvel)=.TRUE.
        END IF
        IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVvel))) THEN
          have_var(idVvel)=.TRUE.
        END IF
        DO itrc=1,NT(ng)
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idTvar(itrc)))) THEN
            have_var(idTvar(itrc))=.TRUE.
          END IF
        END DO
# ifdef SEDIMENT
        DO itrc=1,NST
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idfrac(itrc)))) THEN
            have_var(idfrac(itrc))=.TRUE.
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idBmas(itrc)))) THEN
            have_var(idBmas(itrc))=.TRUE.
          END IF
        END DO
        DO itrc=1,MBEDP
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idSbed(itrc)))) THEN
            have_var(idSbed(itrc))=.TRUE.
          END IF
        END DO
#  ifdef BEDLOAD
        DO itrc=1,NST
          IF (TRIM(varnam(itrc)).eq.TRIM(Vname(1,idUbld(itrc)))) THEN
            have_var(idUbld(itrc))=.TRUE.
          END IF
          IF (TRIM(varnam(itrc)).eq.TRIM(Vname(1,idVbld(itrc)))) THEN
            have_var(idVbld(itrc))=.TRUE.
          END IF
        END DO
#  endif
# endif
# if defined SEDIMENT || defined BBL_MODEL
        DO itrc=1,MBOTP
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idBott(itrc)))) THEN
            have_var(idBott(itrc))=.TRUE.
          END IF
        END DO
# endif
#endif
      END DO
!
!  Check if model initialization variables are available in input
!  NetCDF file.
!
      IF (.not.have_var(idtime)) THEN
        WRITE (stdout,10) TRIM(tvarnam), TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.have_var(idFsur).and.get_var(idFsur)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idFsur)), TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.have_var(idUbar).and.get_var(idUbar)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idUbar)), TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.have_var(idVbar).and.get_var(idVbar)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idVbar)), TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
#ifdef SOLVE3D
      IF (.not.have_var(idUvel).and.get_var(idUvel)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idUvel)), TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.have_var(idVvel).and.get_var(idVvel)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idVvel)), TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      DO itrc=1,NT(ng)
        IF (.not.have_var(idTvar(itrc)).and.get_var(idTvar(itrc))) THEN
          WRITE (stdout,10) TRIM(Vname(1,idTvar(itrc))), TRIM(ncname)
          exit_flag=2
          RETURN
        END IF
      END DO
# ifdef SEDIMENT
      DO itrc=1,NST
        IF (.not.have_var(idfrac(itrc)).and.get_var(idfrac(itrc))) THEN
          WRITE (stdout,10) TRIM(Vname(1,idfrac(itrc))), TRIM(ncname)
          exit_flag=2
          RETURN
        END IF
        IF (.not.have_var(idBmas(itrc)).and.get_var(idBmas(itrc))) THEN
          WRITE (stdout,10) TRIM(Vname(1,idBmas(itrc))), TRIM(ncname)
          exit_flag=2
          RETURN
        END IF
      END DO
      DO itrc=1,MBEDP
        IF (.not.have_var(idSbed(itrc)).and.get_var(idSbed(itrc))) THEN
          WRITE (stdout,10) TRIM(Vname(1,idSbed(itrc))), TRIM(ncname)
          exit_flag=2
          RETURN
        END IF
      END DO
#  ifdef BEDLOAD
      DO itrc=1,NST
        IF (.not.have_var(idUbld(itrc)).and.get_var(idUbld(itrc))) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUbld(itrc))), TRIM(ncname)
          exit_flag=2
          RETURN
        END IF
        IF (.not.have_var(idVbld(itrc)).and.get_var(idVbld(itrc))) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVbld(itrc))), TRIM(ncname)
          exit_flag=2
          RETURN
        END IF
      END DO
#  endif
# endif
# if defined SEDIMENT || defined BBL_MODEL
!
! Check only for the first four required properties, but we looked for
! MBOTP properties above.
!
      DO itrc=1,4
        IF (.not.have_var(idBott(itrc)).and.get_var(idBott(itrc))) THEN
          WRITE (stdout,10) TRIM(Vname(1,idBott(itrc))), TRIM(ncname)
          exit_flag=2
          RETURN
        END IF
      END DO
# endif
#endif
!
!-----------------------------------------------------------------------
!  Read model state from input NetCDF file.
!-----------------------------------------------------------------------
!
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)
!
!  Open input NetCDF file.
!
      status=nf_open(TRIM(ncname), nf_nowrite, ncINPid)
      IF (status.ne.nf_noerr) THEN
        WRITE (stdout,20) TRIM(ncname)
        exit_flag=2
        ioerror=status
        RETURN
      END IF
!
!  If using the latest time record from input NetCDF file as the
!  initialization record, find time record index to read and read
!  in model time.
!
      status=nf_inq_varid(ncINPid, TRIM(tvarnam), varid)
      IF (LastRec) THEN
        Tmax=-1.0_r8
        DO i=1,tsize
          status=nf_get_var1_TYPE(ncINPid,varid,i,time(ng))
          IF (status.ne.nf_noerr) THEN
            WRITE (stdout,30) TRIM(tvarnam), i, TRIM(ncname)
            exit_flag=2
            ioerror=status
            RETURN
          END IF
          IF (time(ng).gt.Tmax) THEN
            Tmax=time(ng)
            IniRec=i
          END IF
        END DO
        time(ng)=Tmax
        InpRec=IniRec
      ELSE
        IF ((IniRec.ne.0).and.(IniRec.gt.tsize)) THEN
          WRITE (stdout,40) IniRec, TRIM(ncname), tsize
          exit_flag=2
          RETURN
        END IF
        IF (IniRec.ne.0) THEN
          InpRec=IniRec
        ELSE
          InpRec=1
        END IF
        status=nf_get_var1_TYPE(ncINPid, varid, InpRec, time(ng))
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,30) TRIM(tvarnam), InpRec, TRIM(ncname)
          exit_flag=2
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Inquire time units. Set local time scale.
!
      time_scale=0.0_r8
      status=nf_inq_var(ncINPid, varid, tvarnam, vartype, nvdims,       &
     &                  vdims(1,varid), nvatts)
      IF (status.eq.nf_noerr) THEN
        DO i=1,nvatts
          status=nf_inq_attname(ncINPid, varid, i, attnam)
          IF (status.eq.nf_noerr) THEN
            IF (TRIM(attnam).eq.'units') THEN
              status=nf_get_att_text(ncINPid, varid, TRIM(attnam),      &
     &                               tunits)
              IF (status.eq.nf_noerr) THEN
                IF (tunits(1:3).eq.'day') THEN
                  time_scale=day2sec
                ELSE IF (tunits(1:6).eq.'second') THEN
                 time_scale=1.0_r8
                END IF
              ELSE
                WRITE (stdout,50) TRIM(attnam)
                exit_flag=2
                ioerror=status
                RETURN
              END IF
            END IF
          ELSE
            WRITE (stdout,60) TRIM(tvarnam)
            exit_flag=2
            RETURN
          END IF
        END DO
      ELSE
        WRITE (stdout,70) TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
!
!  Set starting time index and time clock in days.
!
      time(ng)=time(ng)*time_scale
      tdays(ng)=time(ng)*sec2day
      ntstart=NINT((time(ng)-dstart*day2sec)/dt(ng))+1
      IF (ntstart.lt.1) ntstart=1
      ntfirst=ntstart
!
!  Over-write "IniRec" to the actual initial record processed.
!
      IniRec=InpRec
!
!-----------------------------------------------------------------------
!  Read in state initial variables.
!-----------------------------------------------------------------------
!
      lstr=SCAN(ncname,'/',BACK=.TRUE.)+1
      lend=LEN_TRIM(ncname)
      IF (InpThread) THEN
        IF (ERend.gt.ERstr) THEN
          WRITE (stdout,80) TRIM(Mstate(msg)), tdays(ng), Nrun,         &
     &                      ncname(lstr:lend), InpRec, Tindex
        ELSE
          WRITE (stdout,90) TRIM(Mstate(msg)), tdays(ng),               &
     &                      ncname(lstr:lend), InpRec, Tindex
        END IF
      END IF
!
!  Read in free-surface (m).
!
      IF (get_var(idFsur)) THEN
        status=nf_inq_varid(ncINPid, TRIM(Vname(1,idFsur)), varid)
        gtype=vflag(varid)*r2dvar
        status=nf_fread(ng, model, ncINPid, varid, InpRec, gtype,       &
     &                  LBi, UBi, LBj, UBj, 1, 1,                       &
     &                  Fscl, Fmin, Fmax,                               &
#ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
#endif
     &                  OCEAN(ng) % zeta(LBi,LBj,Tindex))
        IF (InpThread) THEN
          IF (status.ne.nf_noerr) THEN
            WRITE (stdout,30) TRIM(Vname(1,idFsur)), InpRec,            &
     &                        TRIM(ncname)
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            WRITE (stdout,100) TRIM(Vname(2,idFsur)), Fmin, Fmax
          END IF
        END IF
      END IF
!
!  Read in 2D momentum component (m/s) in the XI-direction.
!
      IF (get_var(idUbar)) THEN
        status=nf_inq_varid(ncINPid, TRIM(Vname(1,idUbar)), varid)
        gtype=vflag(varid)*u2dvar
        status=nf_fread(ng, model, ncINPid, varid, InpRec, gtype,       &
     &                  LBi, UBi, LBj, UBj, 1, 1,                       &
     &                  Fscl, Fmin, Fmax,                               &
#ifdef MASKING
     &                  GRID(ng) % umask(LBi,LBj),                      &
#endif
     &                  OCEAN(ng) % ubar(LBi,LBj,Tindex))
        IF (InpThread) THEN
          IF (status.ne.nf_noerr) THEN
            WRITE (stdout,30) TRIM(Vname(1,idUbar)), InpRec,            &
     &                        TRIM(ncname)
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            WRITE (stdout,100) TRIM(Vname(2,idUbar)), Fmin, Fmax
          END IF
        END IF
      END IF
!
!  Read in 2D momentum component (m/s) in the ETA-direction.
!
      IF (get_var(idVbar)) THEN
        status=nf_inq_varid(ncINPid, TRIM(Vname(1,idVbar)), varid)
        gtype=vflag(varid)*v2dvar
        status=nf_fread(ng, model, ncINPid, varid, InpRec, gtype,       &
     &                  LBi, UBi, LBj, UBj, 1, 1,                       &
     &                  Fscl, Fmin, Fmax,                               &
#ifdef MASKING
     &                  GRID(ng) % vmask(LBi,LBj),                      &
#endif
     &                  OCEAN(ng) % vbar(LBi,LBj,Tindex))
        IF (InpThread) THEN
          IF (status.ne.nf_noerr) THEN
            WRITE (stdout,30) TRIM(Vname(1,idVbar)), InpRec,            &
     &                        TRIM(ncname)
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            WRITE (stdout,100) TRIM(Vname(2,idVbar)), Fmin, Fmax
          END IF
        END IF
      END IF
#ifdef SOLVE3D
!
!  Read in 3D momentum component (m/s) in the XI-direction.
!
      IF (get_var(idUvel)) THEN
        status=nf_inq_varid(ncINPid, TRIM(Vname(1,idUvel)), varid)
        gtype=vflag(varid)*u3dvar
        status=nf_fread(ng, model, ncINPid, varid, InpRec, gtype,       &
     &                  LBi, UBi, LBj, UBj, 1, N(ng),                   &
     &                  Fscl, Fmin, Fmax,                               &
# ifdef MASKING
     &                  GRID(ng) % umask(LBi,LBj),                      &
# endif
     &                  OCEAN(ng) % u(LBi,LBj,1,Tindex))
        IF (InpThread) THEN
          IF (status.ne.nf_noerr) THEN
            WRITE (stdout,30) TRIM(Vname(1,idUvel)), InpRec,            &
     &                        TRIM(ncname)
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            WRITE (stdout,100) TRIM(Vname(2,idUvel)), Fmin, Fmax
          END IF
        END IF
      END IF
!
!  Read in 3D momentum component (m/s) in the ETA-direction.
!

      IF (get_var(idVvel)) THEN
        status=nf_inq_varid(ncINPid, TRIM(Vname(1,idVvel)), varid)
        gtype=vflag(varid)*v3dvar
        status=nf_fread(ng, model, ncINPid, varid, InpRec, gtype,       &
     &                  LBi, UBi, LBj, UBj, 1, N(ng),                   &
     &                  Fscl, Fmin, Fmax,                               &
# ifdef MASKING
     &                  GRID(ng) % vmask(LBi,LBj),                      &
# endif
     &                  OCEAN(ng) % v(LBi,LBj,1,Tindex))
        IF (InpThread) THEN
          IF (status.ne.nf_noerr) THEN
            WRITE (stdout,30) TRIM(Vname(1,idVvel)), InpRec,            &
     &                        TRIM(ncname)
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            WRITE (stdout,100) TRIM(Vname(2,idVvel)), Fmin, Fmax
          END IF
        END IF
      END IF
!
!  Read in tracer type variables.
!
      DO itrc=1,NT(ng)
        IF (get_var(idTvar(itrc))) THEN
          status=nf_inq_varid(ncINPid, TRIM(Vname(1,idTvar(itrc))),     &
     &                        varid)
          gtype=vflag(varid)*r3dvar
          status=nf_fread(ng, model, ncINPid, varid, InpRec, gtype,     &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    Fscl, Fmin, Fmax,                             &
# ifdef MASKING
     &                    GRID(ng) % rmask(LBi,LBj),                    &
# endif
     &                    OCEAN(ng) % t(LBi,LBj,1,Tindex,itrc))
          IF (InpThread) THEN
            IF (status.ne.nf_noerr) THEN
              WRITE (stdout,30) TRIM(Vname(1,idTvar(itrc))), InpRec,    &
     &                          TRIM(ncname)
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              WRITE (stdout,100) TRIM(Vname(2,idTvar(itrc))), Fmin, Fmax
            END IF
          END IF
        END IF
      END DO

# ifdef SEDIMENT
!
!  Read in sediment fraction of each size class in each bed layer.
!
      DO i=1,NST
        IF (get_var(idfrac(i))) THEN
          status=nf_inq_varid(ncINPid, TRIM(Vname(1,idfrac(i))), varid)
          gtype=vflag(varid)*b3dvar
          status=nf_fread(ng, model, ncINPid, varid, InpRec, gtype,     &
     &                    LBi, UBi, LBj, UBj, 1, Nbed,                  &
     &                    Fscl, Fmin, Fmax,                             &
#  ifdef MASKING
     &                    GRID(ng) % rmask(LBi,LBj),                    &
#  endif
     &                    OCEAN(ng) % bed_frac(LBi,LBj,1,i))
          IF (InpThread) THEN
            IF (status.ne.nf_noerr) THEN
              WRITE (stdout,30) TRIM(Vname(1,idfrac(i))), InpRec,       &
     &                          TRIM(ncname)
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              WRITE (stdout,100) TRIM(Vname(2,idfrac(i))), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Read in sediment mass of each size class in each bed layer.
!
        IF (get_var(idBmas(i))) THEN
          status=nf_inq_varid(ncINPid, TRIM(Vname(1,idBmas(i))), varid)
          gtype=vflag(varid)*b3dvar
          status=nf_fread(ng, model, ncINPid, varid, InpRec, gtype,     &
     &                    LBi, UBi, LBj, UBj, 1, Nbed,                  &
     &                    Fscl, Fmin, Fmax,                             &
#  ifdef MASKING
     &                    GRID(ng) % rmask(LBi,LBj),                    &
#  endif
     &                    OCEAN(ng) % bed_mass(LBi,LBj,1,Tindex,i))
          IF (InpThread) THEN
            IF (status.ne.nf_noerr) THEN
              WRITE (stdout,30) TRIM(Vname(1,idBmas(i))), InpRec,       &
     &                          TRIM(ncname)
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              WRITE (stdout,100) TRIM(Vname(2,idBmas(i))), Fmin, Fmax
            END IF
          END IF
        END IF
      END DO
!
!  Read in sediment properties in each bed layer.
!
      DO i=1,MBEDP
        IF (get_var(idSbed(i))) THEN
          status=nf_inq_varid(ncINPid, TRIM(Vname(1,idSbed(i))), varid)
          gtype=vflag(varid)*b3dvar
          status=nf_fread(ng, model, ncINPid, varid, InpRec, gtype,     &
     &                    LBi, UBi, LBj, UBj, 1, Nbed,                  &
     &                    Fscl, Fmin, Fmax,                             &
#  ifdef MASKING
     &                    GRID(ng) % rmask(LBi,LBj),                    &
#  endif
     &                    OCEAN(ng) % bed(LBi,LBj,1,i))
          IF (InpThread) THEN
            IF (status.ne.nf_noerr) THEN
              WRITE (stdout,30) TRIM(Vname(1,idSbed(i))), InpRec,       &
     &                          TRIM(ncname)
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              WRITE (stdout,100) TRIM(Vname(2,idSbed(i))), Fmin, Fmax
            END IF
          END IF
        END IF
      END DO

#  ifdef BEDLOAD
!
!  Read in sediment fraction of bed load.
!
      DO i=1,NST
        IF (get_var(idUbld(i))) THEN
          status=nf_inq_varid(ncINPid, TRIM(Vname(1,idUbld(i))), varid)
          gtype=vflag(varid)*u2dvar
          status=nf_fread(ng, model, ncINPid, varid, InpRec, gtype,     &
     &                    LBi, UBi, LBj, UBj, 1, 1,                     &
     &                    Fscl, Fmin, Fmax,                             &
#   ifdef MASKING
     &                    GRID(ng) % umask(LBi,LBj),                    &
#   endif
     &                    OCEAN(ng) % bedldu(LBi,LBj,i))
          IF (InpThread) THEN
            IF (status.ne.nf_noerr) THEN
              WRITE (stdout,30) TRIM(Vname(1,idUbld(i))), InpRec,       &
     &                          TRIM(ncname)
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              WRITE (stdout,100) TRIM(Vname(2,idUbld(i))), Fmin, Fmax
            END IF
          END IF
        END IF
!
        IF (get_var(idVbld(i))) THEN
          status=nf_inq_varid(ncINPid, TRIM(Vname(1,idVbld(i))), varid)
          gtype=vflag(varid)*v2dvar
          status=nf_fread(ng, model, ncINPid, varid, InpRec, gtype,     &
     &                    LBi, UBi, LBj, UBj, 1, 1,                     &
     &                    Fscl, Fmin, Fmax,                             &
#   ifdef MASKING
     &                    GRID(ng) % vmask(LBi,LBj),                    &
#   endif
     &                    OCEAN(ng) % bedldv(LBi,LBj,i))
          IF (InpThread) THEN
            IF (status.ne.nf_noerr) THEN
              WRITE (stdout,30) TRIM(Vname(1,idVbld(i))), InpRec,       &
     &                          TRIM(ncname)
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              WRITE (stdout,100) TRIM(Vname(2,idVbld(i))), Fmin, Fmax
            END IF
          END IF
        END IF
      END DO
#  endif
# endif

# if defined SEDIMENT || defined BBL_MODEL
!
!  Read in sediment properties in exposed bed layer.
!
      DO i=1,MBOTP
        IF (get_var(idBott(i)).and.have_var(idBott(i))) THEN
          status=nf_inq_varid(ncINPid, TRIM(Vname(1,idBott(i))), varid)
          gtype=vflag(varid)*r2dvar
          status=nf_fread(ng, model, ncINPid, varid, InpRec, gtype,     &
     &                    LBi, UBi, LBj, UBj, 1, 1,                     &
     &                    Fscl, Fmin, Fmax,                             &
#  ifdef MASKING
     &                    GRID(ng) % rmask(LBi,LBj),                    &
#  endif
     &                    OCEAN(ng) % bottom(LBi,LBj,i))
          IF (InpThread) THEN
            IF (status.ne.nf_noerr) THEN
              WRITE (stdout,30) TRIM(Vname(1,idBott(i))), InpRec,       &
     &                          TRIM(ncname)
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              WRITE (stdout,100) TRIM(Vname(2,idBott(i))), Fmin, Fmax
            END IF
          END IF
        END IF
      END DO
# endif

#endif

#if (defined EW_PERIODIC || defined NS_PERIODIC) && !defined DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Set periodic boundary conditions.
!-----------------------------------------------------------------------
!
      IF (get_var(idFsur)) THEN
        CALL exchange_r2d_tile (ng, model, 1, Lm(ng), 1, Mm(ng),        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          OCEAN(ng) % zeta(:,:,Tindex))
      END IF
      IF (get_var(idUbar)) THEN
        CALL exchange_u2d_tile (ng, model, 1, Lm(ng), 1, Mm(ng),        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          OCEAN(ng) % ubar(:,:,Tindex))
      END IF
      IF (get_var(idVbar)) THEN
        CALL exchange_v2d_tile (ng, model, 1, Lm(ng), 1, Mm(ng),        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          OCEAN(ng) % vbar(:,:,Tindex))
      END IF
# ifdef SOLVE3D
      IF (get_var(idUvel)) THEN
        CALL exchange_u3d_tile (ng, model, 1, Lm(ng), 1, Mm(ng),        &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          OCEAN(ng) % u(:,:,:,Tindex))
      END IF
      IF (get_var(idVvel)) THEN
        CALL exchange_v3d_tile (ng, model, 1, Lm(ng), 1, Mm(ng),        &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          OCEAN(ng) % v(:,:,:,Tindex))
      END IF
      DO itrc=1,NT(ng)
        IF (get_var(idTvar(itrc))) THEN
          CALL exchange_r3d_tile (ng, model, 1, Lm(ng), 1, Mm(ng),      &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            NghostPoints,                         &
     &                            OCEAN(ng) % t(:,:,:,Tindex,itrc))
        END IF
      END DO

#  ifdef SEDIMENT
      DO i=1,NST
        IF (get_var(idfrac(i))) THEN
          CALL exchange_r3d_tile (ng, model, 1, Lm(ng), 1, Mm(ng),      &
     &                            LBi, UBi, LBj, UBj, 1, Nbed,          &
     &                            NghostPoints,                         &
     &                            OCEAN(ng) % bed_frac(:,:,:,i))
        END IF
        IF (get_var(idBmas(i))) THEN
          CALL exchange_r3d_tile (ng, model, 1, Lm(ng), 1, Mm(ng),      &
     &                            LBi, UBi, LBj, UBj, 1, Nbed,          &
     &                            NghostPoints,                         &
     &                            OCEAN(ng) % bed_mass(:,:,:,Tindex,i))
        END IF
#   ifdef BEDLOAD
        IF (get_var(idUbld(i))) THEN
          CALL exchange_u2d_tile (ng, iNLM, Istr, Iend, Jstr, Jend,     &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            NghostPoints,                         &
     &                            OCEAN(ng) % bedldu(:,:,ised))
        END IF
        IF (get_var(idVbld(i))) THEN
          CALL exchange_v2d_tile (ng, iNLM, Istr, Iend, Jstr, Jend,     &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            NghostPoints,                         &
     &                            OCEAN(ng) % bedldv(:,:,ised))
        END IF
#   endif
      END DO
      DO i=1,MBEDP
        IF (get_var(idSbed(i))) THEN
          CALL exchange_r3d_tile (ng, model, 1, Lm(ng), 1, Mm(ng),      &
     &                            LBi, UBi, LBj, UBj, 1, Nbed,          &
     &                            NghostPoints,                         &
     &                            OCEAN(ng) % bed(:,:,:,i))
        END IF
      END DO
#  endif

#  if defined SEDIMENT || defined BBL_MODEL
      DO i=1,MBOTP
        IF (get_var(idBott(i)).and.have_var(idBott(i))) THEN
          CALL exchange_r2d_tile (ng, model, 1, Lm(ng), 1, Mm(ng),      &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            NghostPoints,                         &
     &                            OCEAN(ng) % bottom(:,:,i))
        END IF
      END DO
#  endif
# endif
#endif
!
!-----------------------------------------------------------------------
!  Close input NetCDF file.
!-----------------------------------------------------------------------
!
      status=nf_close(ncINPid)
!
  10  FORMAT (/,' GET_STATE - unable to find model variable: ',a,       &
     &        /,15x,'in input NetCDF file: ',a)
  20  FORMAT (/,' GET_STATE - unable to open input NetCDF file: ',a)
  30  FORMAT (/,' GET_STATE - error while reading variable: ',a,2x,     &
     &        'at time record = ',i3,/,15x,'in input NetCDF file: ',a)
  40  FORMAT (/,' GET_STATE - requested input time record = ',i3,/,     &
     &        15x,'not found in input NetCDF: ',a,/,                    &
     &        15x,'number of available records = ',i3)
  50  FORMAT (/,' GET_STATE - error while reading attribute: ',a)
  60  FORMAT (/,' GET_STATE - error while inquiring attributes for',    &
     &          ' variable: ',a)
  70  FORMAT (/,' GET_STATE - cannot inquire about time variable in',   &
     &          ' input NetCDF file: ',a)
  80  FORMAT (/,6x,'GET_STATE - ',a,t64,'t = ',f12.4,                   &
     &        /,19x,'(Iter=',i4.4,', File: ',a,', Rec=',i4.4,           &
     &        ', Index=',i1,')')
  90  FORMAT (/,6x,'GET_STATE - ',a,t64,'t = ',f12.4,                   &
     &        /,19x,'(File: ',a,', Rec=',i4.4,', Index=',i1,')')
 100  FORMAT (16x,'- ',a,/,19x,'(Min = ',1p,e15.8,                      &
     &        ' Max = ',1p,e15.8,')')
      RETURN
      END SUBROUTINE get_state
