#include "cppdefs.h"
      SUBROUTINE get_initial (ng)
!
!=======================================================================
!  Copyright (c) 2002  ROMS/TOMS Group                                 !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This subroutine reads in primitive variables initial conditions     !
!  from an initial, restart, or history NetCDF file.                   !
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
#ifdef SEDIMENT
      USE mod_sediment
#endif

#if defined EW_PERIODIC || defined NS_PERIODIC || !defined DISTRIBUTE
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
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      logical, dimension(NV) :: got_var
#ifdef SOLVE3D
      logical, dimension(NT(ng)) :: got_trc

      integer :: itrc
#endif

      integer :: LBi, UBi, LBj, UBj
      integer :: gtype, i, ifield, ncinpid
      integer :: nvatts, nrec, status, tindx, varid, vartype

      integer :: nf_fread

      real(r8), parameter :: Fscl = 1.0_r8

      real(r8) :: Fmax, Fmin, Tmax, time_scale

      character (len=15) attnam, tvarnam
      character (len=40) tunits
      character (len=80) fname, ncname
!
!-----------------------------------------------------------------------
!  Inquire about the contents of input NetCDF file:  Inquire about
!  the dimensions and variables.  Check for consistency.
!-----------------------------------------------------------------------
!
      if (exit_flag.ne.0) RETURN
      ncname=ininame(ng)
      ifield=0
      CALL opencdf (ng, 1, TRIM(ncname), fname, N(ng), ifield, nrec)
      if (exit_flag.ne.0) RETURN
!
!  Initialize logical switches.
!
        DO i=1,NV
          got_var(i)=.false.
        END DO
#ifdef SOLVE3D
        DO itrc=1,NT(ng)
          got_trc(itrc)=.false.
        END DO
#endif
!
!  Scan variable list from input NetCDF and activate switches for
!  model initialization variables.
!

      DO i=1,nvars
        IF ((TRIM(varnam(i)).eq.TRIM(Vname(1,idtime))).or.              &
     &      (TRIM(varnam(i)).eq.'scrum_time')) THEN
          tvarnam=TRIM(varnam(i))
          got_var(idtime)=.true.
        END IF
        IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idFsur))) THEN
          got_var(idFsur)=.true.
        END IF
        IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idUbar))) THEN
          got_var(idUbar)=.true.
        END IF
        IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVbar))) THEN
          got_var(idVbar)=.true.
        END IF
#ifdef SOLVE3D
        IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idUvel))) THEN
          got_var(idUvel)=.true.
        END IF
        IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVvel))) THEN
          got_var(idVvel)=.true.
        END IF
        DO itrc=1,NT(ng)
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idTvar(itrc)))) THEN
            got_trc(itrc)=.true.
          END IF
        END DO
#endif
      END DO
!
!  Check if model initialization variables are available in input
!  NetCDF file.
!
      IF (.not.got_var(idtime)) THEN
        WRITE (stdout,10) TRIM(tvarnam), TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.got_var(idFsur)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idFsur)), TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.got_var(idUbar)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idUbar)), TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.got_var(idVbar)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idVbar)), TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
#ifdef SOLVE3D
      IF (.not.got_var(idUvel)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idUvel)), TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.got_var(idVvel)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idVvel)), TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
# ifdef ANA_BIOLOGY
      DO itrc=1,NAT
# else
      DO itrc=1,NT(ng)
# endif
        IF (.not.got_trc(itrc)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idTvar(itrc))), TRIM(ncname)
          exit_flag=2
          RETURN
        END IF
      END DO
#endif
!
!-----------------------------------------------------------------------
!  Read initial conditions from input NetCDF file.
!-----------------------------------------------------------------------
!
      LBi=lbound(GRID(ng)%h,DIM=1)
      UBi=ubound(GRID(ng)%h,DIM=1)
      LBj=lbound(GRID(ng)%h,DIM=2)
      UBj=ubound(GRID(ng)%h,DIM=2)
!
!  Open input NetCDF file.
!
      status=nf_open(TRIM(ncname), nf_nowrite, ncinpid)
      IF (status.ne.nf_noerr) THEN
        WRITE (stdout,20) TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
!
!  If using the latest time record from initial NetCDF file, set
!  initialization record. Set time record index to read.  Read in
!  model time.
!
      status=nf_inq_varid(ncinpid, TRIM(tvarnam), varid)
      IF (LastRec) THEN
        Tmax=-1.0_r8
        DO i=1,tsize
          status=nf_get_var1_TYPE(ncinpid,varid,i,time(ng))
          IF (status.ne.nf_noerr) THEN
            WRITE (stdout,30) TRIM(tvarnam), i, TRIM(ncname)
            exit_flag=2
            RETURN
          END IF
          IF (time(ng).gt.Tmax) THEN
            Tmax=time(ng)
            nrrec=i
          END IF
        END DO
        time(ng)=Tmax
        tindx=nrrec
      ELSE
        IF ((nrrec.ne.0).and.(nrrec.gt.tsize)) THEN
          WRITE (stdout,40) nrrec, TRIM(ncname), tsize
          exit_flag=2
          RETURN
        END IF
        IF (nrrec.ne.0) THEN
          tindx=nrrec
        ELSE
          tindx=1
        END IF
        status=nf_get_var1_TYPE(ncinpid, varid, tindx, time(ng))
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,30) TRIM(tvarnam), tindx, TRIM(ncname)
          exit_flag=2
          RETURN
        END IF
      END IF
!
!  Inquire time units. Set local time scale.
!
      time_scale=0.0_r8
      status=nf_inq_var(ncinpid, varid, tvarnam, vartype, nvdims,       &
     &                  vdims(1,varid), nvatts)
      IF (status.eq.nf_noerr) THEN
        DO i=1,nvatts
          status=nf_inq_attname(ncinpid, varid, i, attnam)
          IF (status.eq.nf_noerr) THEN
            IF (TRIM(attnam).eq.'units') THEN
              status=nf_get_att_text(ncinpid, varid, TRIM(attnam),      &
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
!  Read in free-surface (m).
!
      status=nf_inq_varid(ncinpid, TRIM(Vname(1,idFsur)), varid)
      gtype=vflag(varid)*r2dvar
      status=nf_fread(ng, ncinpid, varid, tindx, gtype,                 &
     &                LBi, UBi, LBj, UBj, 1, 1,                         &
     &                Fscl, Fmin, Fmax,                                 &
#ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#endif
     &                OCEAN(ng) % zeta(LBi,LBj,1))
      IF (InpThread.and.(status.ne.nf_noerr)) THEN
        WRITE (stdout,40) TRIM(Vname(1,idFsur)), tindx, TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
!
!  Read in 2D momentum component (m/s) in the XI-direction.
!
      status=nf_inq_varid(ncinpid, TRIM(Vname(1,idUbar)), varid)
      gtype=vflag(varid)*u2dvar
      status=nf_fread(ng, ncinpid, varid, tindx, gtype,                 &
     &                LBi, UBi, LBj, UBj, 1, 1,                         &
     &                Fscl, Fmin, Fmax,                                 &
#ifdef MASKING
     &                GRID(ng) % umask(LBi,LBj),                        &
#endif
     &                OCEAN(ng) % ubar(LBi,LBj,1))
      IF (InpThread.and.(status.ne.nf_noerr)) THEN
        WRITE (stdout,40) TRIM(Vname(1,idUbar)), tindx, TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
!
!  Read in 2D momentum component (m/s) in the ETA-direction.
!
      status=nf_inq_varid(ncinpid, TRIM(Vname(1,idVbar)), varid)
      gtype=vflag(varid)*v2dvar
      status=nf_fread(ng, ncinpid, varid, tindx, gtype,                 &
     &                LBi, UBi, LBj, UBj, 1, 1,                         &
     &                Fscl, Fmin, Fmax,                                 &
#ifdef MASKING
     &                GRID(ng) % vmask(LBi,LBj),                        &
#endif
     &                OCEAN(ng) % vbar(LBi,LBj,1))
      IF (InpThread.and.(status.ne.nf_noerr)) THEN
        WRITE (stdout,40) TRIM(Vname(1,idVbar)), tindx, TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
#ifdef SOLVE3D
!
!  Read in 3D momentum component (m/s) in the XI-direction.
!
      status=nf_inq_varid(ncinpid, TRIM(Vname(1,idUvel)), varid)
      gtype=vflag(varid)*u3dvar
      status=nf_fread(ng, ncinpid, varid, tindx, gtype,                 &
     &                LBi, UBi, LBj, UBj, 1, N(ng),                     &
     &                Fscl, Fmin, Fmax,                                 &
#ifdef MASKING
     &                GRID(ng) % umask(LBi,LBj),                        &
#endif
     &                OCEAN(ng) % u(LBi,LBj,1,1))
      IF (InpThread.and.(status.ne.nf_noerr)) THEN
        WRITE (stdout,40) TRIM(Vname(1,idUvel)), tindx, TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
!
!  Read in 3D momentum component (m/s) in the ETA-direction.
!
      status=nf_inq_varid(ncinpid, TRIM(Vname(1,idVvel)), varid)
      gtype=vflag(varid)*v3dvar
      status=nf_fread(ng, ncinpid, varid, tindx, gtype,                 &
     &                LBi, UBi, LBj, UBj, 1, N(ng),                     &
     &                Fscl, Fmin, Fmax,                                 &
#ifdef MASKING
     &                GRID(ng) % vmask(LBi,LBj),                        &
#endif
     &                OCEAN(ng) % v(LBi,LBj,1,1))
      IF (InpThread.and.(status.ne.nf_noerr)) THEN
        WRITE (stdout,40) TRIM(Vname(1,idVvel)), tindx, TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
!
!  Read in tracer type variables.
!
# ifdef ANA_BIOLOGY
      DO itrc=1,NAT
# else
      DO itrc=1,NT(ng)
# endif
        status=nf_inq_varid(ncinpid, TRIM(Vname(1,idTvar(itrc))), varid)
        gtype=vflag(varid)*r3dvar
        status=nf_fread(ng, ncinpid, varid, tindx, gtype,               &
     &                  LBi, UBi, LBj, UBj, 1, N(ng),                   &
     &                  Fscl, Fmin, Fmax,                               &
#ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
#endif
     &                  OCEAN(ng) % t(LBi,LBj,1,1,itrc))
        IF (InpThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,40) TRIM(Vname(1,idTvar(itrc))), tindx,         &
     &                      TRIM(ncname)
          exit_flag=2
          RETURN
        END IF
      END DO
# ifdef SEDIMENT
      DO i=1,NST
        status=nf_inq_varid(ncinpid, TRIM(Vname(1,idfrac(i))), varid)
        gtype=vflag(varid)*b3dvar
        status=nf_fread(ng, ncinpid, varid, tindx, gtype,               &
     &                  LBi, UBi, LBj, UBj, 1, Nbed,                    &
     &                  Fscl, Fmin, Fmax,                               &
#ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
#endif
     &                  OCEAN(ng) % bed_frac(LBi,LBj,1,i))
        IF (InpThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,40) TRIM(Vname(1,idfrac(i))), tindx,            &
     &                      TRIM(ncname)
          exit_flag=2
          RETURN
        END IF
      END DO
      DO i=1,MBEDP
        status=nf_inq_varid(ncinpid, TRIM(Vname(1,idSbed(i))), varid)
        gtype=vflag(varid)*b3dvar
        status=nf_fread(ng, ncinpid, varid, tindx, gtype,               &
     &                  LBi, UBi, LBj, UBj, 1, Nbed,                    &
     &                  Fscl, Fmin, Fmax,                               &
#ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
#endif
     &                  OCEAN(ng) % bed(LBi,LBj,1,i))
        IF (InpThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,40) TRIM(Vname(1,idSbed(i))), tindx,            &
     &                      TRIM(ncname)
          exit_flag=2
          RETURN
        END IF
      END DO
# endif
#endif
#if defined EW_PERIODIC || defined NS_PERIODIC || !defined DISTRIBUTE
!
!  Set periodic boundary conditions.
!
      CALL exchange_r2d_tile (ng, 1, Lm(ng), 1, Mm(ng),                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        OCEAN(ng) % zeta(:,:,1))
      CALL exchange_u2d_tile (ng, 1, Lm(ng), 1, Mm(ng),                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        OCEAN(ng) % ubar(:,:,1))
      CALL exchange_v2d_tile (ng, 1, Lm(ng), 1, Mm(ng),                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        OCEAN(ng) % vbar(:,:,1))
# ifdef SOLVE3D
      CALL exchange_u3d_tile (ng, 1, Lm(ng), 1, Mm(ng),                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        OCEAN(ng) % u(:,:,:,1))
      CALL exchange_v3d_tile (ng, 1, Lm(ng), 1, Mm(ng),                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        OCEAN(ng) % v(:,:,:,1))
#  ifdef ANA_BIOL0GY
      DO itrc=1,NAT
#  else
      DO itrc=1,NT(ng)
#  endif
        CALL exchange_r3d_tile (ng, 1, Lm(ng), 1, Mm(ng),               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          OCEAN(ng) % t(:,:,:,1,itrc))
      END DO
#  ifdef SEDIMENT
      DO i=1,NST
        CALL exchange_r3d_tile (ng, 1, Lm(ng), 1, Mm(ng),               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          OCEAN(ng) % bed_frac(:,:,:,i))
      END DO
      DO i=1,MBEDP
        CALL exchange_r3d_tile (ng, 1, Lm(ng), 1, Mm(ng),               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          OCEAN(ng) % bed(:,:,:,i))
      END DO
#  endif
# endif
#endif
!
!  Close input NetCDF file.
!
      status=nf_close(ncinpid)
      IF (Master) WRITE (stdout,80) tdays(ng)
!
  10  format(/,' GET_INITIAL - unable to find model variable: ',a,      &
     &       /,15x,'in input NetCDF file: ',a)
  20  format(/,' GET_INITIAL - unable to open input NetCDF file: ',a)
  30  format(/,' GET_INITIAL - error while reading variable: ',a,2x,    &
     &       'at time record = ',i3,/,15x,'in input NetCDF file: ',a)
  40  format(/,' GET_INITIAL - requested restart time record = ',i3,/,  &
     &       15x,'not found in input NetCDF: ',a,/,                     &
     &       15x,'number of available records = ',i3)
  50  format(/,' GET_INITIAL - error while reading attribute: ',a)
  60  format(/,' GET_INITIAL - error while inquiring attributes for',   &
     &       ' variable: ',a)
  70  format(/,' GET_INITIAL - cannot inquire about time variable in',  &
     &       ' input NetCDF file: ',a)
  80  format(/,' GET_INITIAL - Processed initial conditions,',          &
     &       t64,'t = ',f12.4)
      RETURN
      END SUBROUTINE get_initial
