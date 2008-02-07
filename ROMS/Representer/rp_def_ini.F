#include "cppdefs.h"
#if defined TL_IOMS && defined FOUR_DVAR
      SUBROUTINE rp_def_ini (ng)
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!======================================================================= 
!                                                                      !
!  This routine updates representer model initial NetCDF file used     !
!  for weak constraint, 4DVAR data assimilation.  It defines new       !
!  dimensions and variables.                                           !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
# ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcasti, mp_bcastl
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
      integer, parameter :: Natt = 24

      logical :: Ldefine, got_var(NV)

      integer :: i, j, nrec, nvd, nvd4
      integer :: frecdim, recdim, status, varid
# ifdef SOLVE3D
      integer :: itrc
# endif

      integer :: DimIDs(29), t3dfrc(4), u3dfrc(4), v3dfrc(4)
      integer :: Vsize(4)

      integer :: def_var

      real(r8) :: Aval(6)

      character (len=80) :: Vinfo(Natt)
      character (len=80) :: fname, ncname

# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
!
!=======================================================================
!  Open existing representer model initial conditions and define new
!  new variables.
!=======================================================================
!
      IF (exit_flag.ne.NoError) RETURN
      ncname=IRPname(ng)
!
!  Open NetCDF file.
!
      IF (OutThread) THEN
        status=nf90_open(TRIM(ncname),nf90_write,ncIRPid(ng))
      END IF  
#  ifdef DISTRIBUTE
      CALL mp_bcasti (ng, iRPM, status, 1)
#  endif
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(ncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Inquire about the contents of existing NetCDF file:  Inquire
!  about the dimensions and variables.  Check for consistency.
!
      IF (outThread) THEN
        CALL opencdf (ng, 1, ncname, fname, N(ng), 0, nrec, nvd, Vsize)
        IF (exit_flag.ne.NoError) RETURN
!
!  Check if surface forcing variables have been already defined.
!
        DO i=1,NV
          got_var(i)=.FALSE.
        END DO
!
        DO i=1,nvars
# ifdef ADJUST_WSTRESS
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idUsms))) THEN
            got_var(idUsms)=.TRUE.
            status=nf90_inq_varid(ncIRPid(ng),TRIM(Vname(1,idUsms)),    &
     &                            irpVid(idUsms,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVsms))) THEN
            got_var(idVsms)=.TRUE.
            status=nf90_inq_varid(ncIRPid(ng),TRIM(Vname(1,idVsms)),    &
     &                            irpVid(idVsms,ng))
          END IF
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
          DO itrc=1,NT(ng)
            IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idTsur(itrc)))) THEN
              got_var(idTsur(itrc))=.TRUE.
              status=nf90_inq_varid(ncIRPid(ng),                        &
     &                              TRIM(Vname(1,idTsur(itrc))),        &
     &                              irpVid(idTsur(itrc),ng))
            END IF
          END DO
# endif
        END DO
        Ldefine=.FALSE.
# ifdef ADJUST_WSTRESS
        IF (.not.got_var(idUsms)) Ldefine=.TRUE.
        IF (.not.got_var(idVsms)) Ldefine=.TRUE.
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
          DO itrc=1,NT(ng)
            IF (.not.got_var(idTsur(itrc))) Ldefine=.TRUE.
          END DO
# endif
      END IF
# ifdef DISTRIBUTE
      CALL mp_bcastl (ng, iNLM, Ldefine, 1)
# endif
!
!  Put existing file into define mode so new variables can be added.
!
      IF (Ldefine) THEN
        IF (OutThread) THEN
          status=nf90_redef(ncIRPid(ng))
        END IF
#  ifdef DISTRIBUTE
        CALL mp_bcasti (ng, iRPM, status, 1)
#  endif
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Define the dimensions of staggered fields.
!-----------------------------------------------------------------------
!
      IF (Ldefine.and.OutThread) THEN
        status=nf90_inq_dimid(ncIRPid(ng),'xi_rho' ,DimIDs( 1))
        status=nf90_inq_dimid(ncIRPid(ng),'xi_u'   ,DimIDs( 2))
        status=nf90_inq_dimid(ncIRPid(ng),'xi_v'   ,DimIDs( 3))
        status=nf90_inq_dimid(ncIRPid(ng),'eta_rho',DimIDs( 5))
        status=nf90_inq_dimid(ncIRPid(ng),'eta_u'  ,DimIDs( 6))
        status=nf90_inq_dimid(ncIRPid(ng),'eta_v'  ,DimIDs( 7))
#  ifdef SOLVE3D
        status=nf90_inq_dimid(ncIRPid(ng),'s_rho'  ,DimIDs( 9))
        status=nf90_inq_dimid(ncIRPid(ng),'s_w'    ,DimIDs(10))
#  endif
        status=nf90_inquire(ncIRPid(ng), unlimitedDimID = recdim)
        DimIDs(12)=recdim
        status=nf90_def_dim(ncIRPid(ng),'Nadjust',Nfrec(ng),frecdim)
!
!  Set number of dimensions for output variables.
!
        nvd4=4
!
!  Define dimension vectors for staggered tracer type variables.
!
        t3dfrc(1)=DimIDs( 1)
        t3dfrc(2)=DimIDs( 5)
        t3dfrc(3)=frecdim
        t3dfrc(4)=recdim
!
!  Define dimension vectors for staggered u-momemtum type variables.
!
        u3dfrc(1)=DimIDs( 2)
        u3dfrc(2)=DimIDs( 6)
        u3dfrc(3)=frecdim
        u3dfrc(4)=recdim
!
!  Define dimension vectors for staggered v-momemtum type variables.
!
        v3dfrc(1)=DimIDs( 3)
        v3dfrc(2)=DimIDs( 7)
        v3dfrc(3)=frecdim
        v3dfrc(4)=recdim
!
!  Initialize local information variable arrays.
!
        DO i=1,Natt
          DO j=1,80
            Vinfo(i)=' '
          END DO
        END DO
        DO i=1,6
          Aval(i)=0.0_r8
        END DO
!
!-----------------------------------------------------------------------
!  Define surface forcing variables.
!-----------------------------------------------------------------------

#  ifdef ADJUST_WSTRESS
!
!  Define surface U-momentum stress.  Notice that the stress has its
!  own fixed time-dimension (of size Nfrec) to allow 4DVAR adjustments
!  at other times in addition to initialization time.
!
        IF (.not.got_var(idUsms)) THEN
          Vinfo( 1)=Vname(1,idUsms)
          Vinfo( 2)=Vname(2,idUsms)
          Vinfo( 3)=Vname(3,idUsms)
          Vinfo(14)=Vname(4,idUsms)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(u2dvar,r8)
          Vinfo(24)='_FillValue'
          Aval(6)=0.0_r8
          status=def_var(ncIRPid(ng),irpVid(idUsms,ng),NF_FOUT,         &
     &                   nvd4,u3dfrc,Aval,Vinfo,ncname)
        END IF
!
!  Define surface V-momentum stress.
!
        IF (.not.got_var(idVsms)) THEN
          Vinfo( 1)=Vname(1,idVsms)
          Vinfo( 2)=Vname(2,idVsms)
          Vinfo( 3)=Vname(3,idVsms)
          Vinfo(14)=Vname(4,idVsms)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(v2dvar,r8)
          Vinfo(24)='_FillValue'
          Aval(6)=0.0_r8
          status=def_var(ncIRPid(ng),irpVid(idVsms,ng),NF_FOUT,         &
     &                   nvd4,v3dfrc,Aval,Vinfo,ncname)
        END IF
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
!
!  define surface tracer flux. Notice that additional tracer fluxes are
!  written at their own fixed time-dimension (of size Nfrec) to allow
!  4DVAR adjustments at other times in addition to initial time.
!
        DO itrc=1,NT(ng)
          IF (.not.got_var(idTsur(itrc))) THEN
            Vinfo( 1)=Vname(1,idTsur(itrc))
            Vinfo( 2)=TRIM(Vname(2,idTsur(itrc)))
            Vinfo( 3)=Vname(3,idTsur(itrc))
            IF (itrc.eq.itemp) THEN
              Vinfo(11)='upward flux, cooling'
              Vinfo(12)='downward flux, heating'
            ELSE IF (itrc.eq.isalt) THEN
              Vinfo(11)='upward flux, freshening (net precipitation)'
              Vinfo(12)='downward flux, salting (net evaporation)'
            END IF
            Vinfo(14)=Vname(4,idTsur(itrc))
            Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#   endif
            Vinfo(22)='coordinates'
            Aval(5)=REAL(r2dvar,r8)
            Vinfo(24)='_FillValue'
            Aval(6)=0.0_r8
            status=def_var(ncIRPid(ng),irpVid(idTsur(itrc),ng),         &
     &                     NF_FOUT,nvd4,t3dfrc,Aval,Vinfo,ncname)
          END IF
        END DO
#  endif
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
        status=nf90_enddef(ncIRPid(ng))
      END IF
# endif
!
!=======================================================================
!  Open an existing tangent linear initial file, check its contents, and
!  prepare for appending data.
!=======================================================================
!
      IF (.not.LdefIRP(ng).and.OutThread) THEN

# if !(defined ADJUST_STFLUX || defined ADJUST_WSTRESS)
!
!  Inquire about the contents of initialization NetCDF file:  Inquire
!  about the dimensions and variables.  Check for consistency.
!
        ncname=IRPname(ng)
        CALL opencdf (ng, 1, ncname, fname, N(ng), 0, nrec, nvd, Vsize)
        IF (exit_flag.ne.NoError) RETURN
!
!  Open initialization file for read/write.
!
        status=nf90_open(TRIM(ncname),nf90_write,ncIRPid(ng))
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,30) TRIM(ncname)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
!  Initialize logical switches.
!
        DO i=1,NV
          got_var(i)=.FALSE.
        END DO
# endif
!
!  Scan variable list from input NetCDF and activate switches for
!  initialization variables. Get variable IDs.
!
        DO i=1,nvars
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idtime))) THEN
            got_var(idtime)=.TRUE.
            status=nf90_inq_varid(ncIRPid(ng),TRIM(Vname(1,idtime)),    &
     &                            irpVid(idtime,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idFsur))) THEN
            got_var(idFsur)=.TRUE.
            status=nf90_inq_varid(ncIRPid(ng),TRIM(Vname(1,idFsur)),    &
     &                            irpVid(idFsur,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idUbar))) THEN
            got_var(idUbar)=.TRUE.
            status=nf90_inq_varid(ncIRPid(ng),TRIM(Vname(1,idUbar)),    &
     &                            irpVid(idUbar,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVbar))) THEN
            got_var(idVbar)=.TRUE.
            status=nf90_inq_varid(ncIRPid(ng),TRIM(Vname(1,idVbar)),    &
     &                            irpVid(idVbar,ng))
          END IF
# ifdef SOLVE3D
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idUvel))) THEN
            got_var(idUvel)=.TRUE.
            status=nf90_inq_varid(ncIRPid(ng),TRIM(Vname(1,idUvel)),    &
     &                            irpVid(idUvel,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVvel))) THEN
            got_var(idVvel)=.TRUE.
            status=nf90_inq_varid(ncIRPid(ng),TRIM(Vname(1,idVvel)),    &
     &                            irpVid(idVvel,ng))
          END IF
          DO itrc=1,NT(ng)
            IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idTvar(itrc)))) THEN
             got_var(idTvar(itrc))=.TRUE.
             status=nf90_inq_varid(ncIRPid(ng),                         &
     &                             TRIM(Vname(1,idTvar(itrc))),         &
     &                             irpTid(itrc,ng))
            END IF
          END DO
#  if defined BVF_MIXING || defined LMD_MIXING || \
      defined GLS_MIXING || defined MY25_MIXING
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVvis))) THEN
            got_var(idVvis)=.TRUE.
            status=nf90_inq_varid(ncIRPid(ng),TRIM(Vname(1,idVvis)),    &
     &                            irpVid(idVvis,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idTdif))) THEN
            got_var(idTdif)=.TRUE.
            status=nf90_inq_varid(ncIRPid(ng),TRIM(Vname(1,idTdif)),    &
     &                            irpVid(idTdif,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idSdif))) THEN
            got_var(idSdif)=.TRUE.
            status=nf90_inq_varid(ncIRPid(ng),TRIM(Vname(1,idSdif)),    &
     &                            irpVid(idSdif,ng))
          END IF
#  endif
# endif
        END DO
!
!  Check if initialization variables are available in input NetCDF
!  file.
!
        IF (.not.got_var(idtime)) THEN
          WRITE (stdout,30) TRIM(Vname(1,idtime)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idFsur)) THEN
          WRITE (stdout,30) TRIM(Vname(1,idFsur)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbar)) THEN
          WRITE (stdout,30) TRIM(Vname(1,idUbar)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbar)) THEN
          WRITE (stdout,30) TRIM(Vname(1,idVbar)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# ifdef SOLVE3D
        IF (.not.got_var(idUvel)) THEN
          WRITE (stdout,30) TRIM(Vname(1,idUvel)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVvel)) THEN
          WRITE (stdout,30) TRIM(Vname(1,idVvel)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        DO itrc=1,NT(ng)
          IF (.not.got_var(idTvar(itrc))) THEN
            WRITE (stdout,30) TRIM(Vname(1,idTvar(itrc))), TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
# endif
      END IF
# ifdef DISTRIBUTE
!
!  Broadcast NetCDF variable indices.
!
      CALL mp_bcasti (ng, iRPM, irpVid(:,ng), NV)
# endif
!
  10  FORMAT (/,' RP_DEF_INI - unable to open initial NetCDF file: ',a)
  20  FORMAT (/,' RP_DEF_INI - unable to put in define mode initial',   &
     &        ' NetCDF file: ',a)
  30  FORMAT (/,' RP_DEF_INI - unable to find variable: ',a,2x,         &
     &        ' in initial NetCDF file: ',a)
      RETURN
      END SUBROUTINE rp_def_ini
#else
      SUBROUTINE rp_def_ini
      RETURN
      END SUBROUTINE rp_def_ini
#endif
