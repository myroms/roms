#include "cppdefs.h"
#ifdef FOUR_DVAR
      SUBROUTINE def_ini (ng)
!
!=======================================================================
!  Copyright (c) 2005 ROMS/TOMS Adjoint Group                          !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine creates a initialization NetCDF file used for          !
!  variational data assimilation.  It defines its dimensions,          !
!  attributes, and variables.                                          !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
# if defined SEDIMENT || defined BBL_MODEL
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
      integer, parameter :: Natt = 23

      logical :: Ldefine, got_var(NV)

      integer :: i, j, ifield, nrec, nvd3, nvd4
      integer :: recdim, status, varid

      integer :: CostDim(2), DimIDs(28), t2dgrd(3), u2dgrd(3), v2dgrd(3)

      integer :: def_info, def_var, wrt_info

# ifdef SOLVE3D
      integer :: itrc

      integer :: t3dgrd(4), u3dgrd(4), v3dgrd(4)
# endif

      real(r8) :: Aval(5)

      character (len=80) :: Vinfo(Natt)
      character (len=80) :: fname, ncname
!
!=======================================================================
!  Create a new initial file.
!=======================================================================
!
      IF (exit_flag.ne.NoError) RETURN
      ncname=INIname(ng)
      IF (LdefINI(ng).and.OutThread) THEN
        status=nf_create(TRIM(ncname),nf_clobber,ncINIid(ng))
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,10) TRIM(ncname)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# ifdef DISTRIBUTE
      CALL mp_bcasti (ng, iNLM, ncINIid(ng), 1)
# endif
!
!-----------------------------------------------------------------------
!  Define the dimensions of staggered fields.
!-----------------------------------------------------------------------
!
      IF (LdefINI(ng).and.OutThread) THEN
        status=nf_def_dim(ncINIid(ng),'xi_rho'  ,Lm(ng)+2,              &
     &                    DimIDs( 1))
        status=nf_def_dim(ncINIid(ng),'xi_u'    ,Lm(ng)+1,              &
     &                    DimIDs( 2))
        status=nf_def_dim(ncINIid(ng),'xi_v'    ,Lm(ng)+2,              &
     &                    DimIDs( 3))
        status=nf_def_dim(ncINIid(ng),'xi_psi'  ,Lm(ng)+1,              &
     &                    DimIDs( 4))
        status=nf_def_dim(ncINIid(ng),'eta_rho' ,Mm(ng)+2,              &
     &                    DimIDs( 5))
        status=nf_def_dim(ncINIid(ng),'eta_u'   ,Mm(ng)+2,              &
     &                    DimIDs( 6))
        status=nf_def_dim(ncINIid(ng),'eta_v'   ,Mm(ng)+1,              &
     &                    DimIDs( 7))
        status=nf_def_dim(ncINIid(ng),'eta_psi' ,Mm(ng)+1,              &
     &                    DimIDs( 8))
# if defined WRITE_WATER && defined MASKING
        status=nf_def_dim(ncINIid(ng),'xy_rho'  ,Nxyr(ng),              &
     &                    DimIDs(17))
        status=nf_def_dim(ncINIid(ng),'xy_u'    ,Nxyu(ng),              &
     &                    DimIDs(18))
        status=nf_def_dim(ncINIid(ng),'xy_v'    ,Nxyv(ng),              &
     &                    DimIDs(19))
# endif
# ifdef SOLVE3D
#  if defined WRITE_WATER && defined MASKING
        status=nf_def_dim(ncINIid(ng),'xyz_rho' ,Nxyr(ng)*N(ng),        &
     &                    DimIDs(20))
        status=nf_def_dim(ncINIid(ng),'xyz_u'   ,Nxyu(ng)*N(ng),        &
     &                    DimIDs(21))
        status=nf_def_dim(ncINIid(ng),'xyz_v'   ,Nxyv(ng)*N(ng),        &
     &                    DimIDs(22))
        status=nf_def_dim(ncINIid(ng),'xyz_w'   ,Nxyr(ng)*(N(ng)+1),    &
     &                    DimIDs(23))
#  endif
        status=nf_def_dim(ncINIid(ng),'s_rho'   ,N(ng),                 &
     &                    DimIDs( 9))
        status=nf_def_dim(ncINIid(ng),'s_w'     ,N(ng)+1,               &
     &                    DimIDs(10))
        status=nf_def_dim(ncINIid(ng),'tracer'  ,NT(ng),                &
     &                    DimIDs(11))
#  ifdef SEDIMENT
        status=nf_def_dim(ncINIid(ng),'Nbed'    ,Nbed,                  &
     &                    DimIDs(16))
#   if defined WRITE_WATER && defined MASKING
        status=nf_def_dim(ncINIid(ng),'xybed'   ,Nxyr(ng)*Nbed,         &
     &                    DimIDs(24))
#   endif
#  endif
#  ifdef ECOSIM
        status=nf_def_dim(ncINIid(ng),'Nphy'    ,Nphy,                  &
     &                    DimIDs(25))
        status=nf_def_dim(ncINIid(ng),'Nbac'    ,Nbac,                  &
     &                    DimIDs(26))
        status=nf_def_dim(ncINIid(ng),'Ndom'    ,Ndom,                  &
     &                    DimIDs(27))
        status=nf_def_dim(ncINIid(ng),'Nfec'    ,Nfec,                  &
     &                    DimIDs(28))
#  endif
# endif
        status=nf_def_dim(ncINIid(ng),'boundary',4,                     &
     &                    DimIDs(14))
        status=nf_def_dim(ncINIid(ng),TRIM(ADJUSTL(Vname(5,idtime))),   &
     &                    nf_unlimited,DimIDs(12))
        recdim=DimIDs(12)
!
!  Set dimension for misfit cost function. Since more 2 unlimited
!  dimensions are not allow in NetCDF currently, set the "iteration"
!  dimension to a fix value.
!
        status=nf_def_dim(ncINIid(ng),'state_var',NstateVar(ng)+1,      &
     &                    CostDim(1))
        status=nf_def_dim(ncINIid(ng),'iteration',1000,                 &
     &                    CostDim(2))
!
!  Set number of dimensions for output variables.
!
# if defined WRITE_WATER && defined MASKING
        nvd3=2
        nvd4=2
# else
        nvd3=3
        nvd4=4
# endif
!
!  Define dimension vectors for staggered tracer type variables.
!
# if defined WRITE_WATER && defined MASKING
        t2dgrd(1)=DimIDs(17)
        t2dgrd(2)=DimIDs(12)
#  ifdef SOLVE3D
        t3dgrd(1)=DimIDs(20)
        t3dgrd(2)=DimIDs(12)
#  endif
# else
        t2dgrd(1)=DimIDs( 1)
        t2dgrd(2)=DimIDs( 5)
        t2dgrd(3)=DimIDs(12)
#  ifdef SOLVE3D
        t3dgrd(1)=DimIDs( 1)
        t3dgrd(2)=DimIDs( 5)
        t3dgrd(3)=DimIDs( 9)
        t3dgrd(4)=DimIDs(12)
#  endif
# endif
!
!  Define dimension vectors for staggered u-momemtum type variables.
!
# if defined WRITE_WATER && defined MASKING
        u2dgrd(1)=DimIDs(18)
        u2dgrd(2)=DimIDs(12)
#  ifdef SOLVE3D
        u3dgrd(1)=DimIDs(21)
        u3dgrd(2)=DimIDs(12)
#  endif
# else
        u2dgrd(1)=DimIDs( 2)
        u2dgrd(2)=DimIDs( 6)
        u2dgrd(3)=DimIDs(12)
#  ifdef SOLVE3D
        u3dgrd(1)=DimIDs( 2)
        u3dgrd(2)=DimIDs( 6)
        u3dgrd(3)=DimIDs( 9)
        u3dgrd(4)=DimIDs(12)
#  endif
# endif
!
!  Define dimension vectors for staggered v-momemtum type variables.
!
# if defined WRITE_WATER && defined MASKING
        v2dgrd(1)=DimIDs(19)
        v2dgrd(2)=DimIDs(12)
#  ifdef SOLVE3D
        v3dgrd(1)=DimIDs(22)
        v3dgrd(2)=DimIDs(12)
#  endif
# else
        v2dgrd(1)=DimIDs( 3)
        v2dgrd(2)=DimIDs( 7)
        v2dgrd(3)=DimIDs(12)
#  ifdef SOLVE3D
        v3dgrd(1)=DimIDs( 3)
        v3dgrd(2)=DimIDs( 7)
        v3dgrd(3)=DimIDs( 9)
        v3dgrd(4)=DimIDs(12)
#  endif
# endif
!
!  Initialize unlimited time record dimension.
!
        tINIindx(ng)=0
!
!  Initialize local information variable arrays.
!
        DO i=1,Natt
          DO j=1,80
            Vinfo(i)=' '
          END DO
        END DO
        DO i=1,5
          Aval(i)=0.0_r8
        END DO
!
!-----------------------------------------------------------------------
!  Define time-recordless information variables.
!-----------------------------------------------------------------------
!
        status=def_info(ng,ncINIid(ng),ncname,DimIDs)
        IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Define time-varying variables.
!-----------------------------------------------------------------------
!
!  Define model time.
!
        Vinfo( 1)=Vname(1,idtime)
        Vinfo( 2)=Vname(2,idtime)
        IF (INT(time_ref).eq.-2) THEN
          Vinfo( 3)='seconds since 1968-05-23 00:00:00 GMT'
          Vinfo( 4)='gregorian'
        ELSE IF (INT(time_ref).eq.-1) THEN
          Vinfo( 3)='seconds since 0001-01-01 00:00:00'
          Vinfo( 4)='360_day'
        ELSE IF (INT(time_ref).eq.0) THEN
          Vinfo( 3)='seconds since 0001-01-01 00:00:00'
          Vinfo( 4)='365.25_day'
        ELSE IF (time_ref.gt.0.0_r8) THEN
          WRITE (Vinfo( 3),'(a,1x,a)') 'seconds since', TRIM(r_text)
          Vinfo( 4)='standard'
        END IF
        Vinfo(14)=Vname(4,idtime)
        status=def_var(ncINIid(ng),iniVid(idtime,ng),NF_TYPE,           &
     &                 1,recdim,Aval,Vinfo,ncname)
!
!  Define model-observations misfit cost function.
!
        Vinfo( 1)='cost_function'
        Vinfo( 2)='model-observation misfit cost function'
        Vinfo(17)='missing_value'
        Aval(1)=0.0_r8
        status=def_var(ncINIid(ng),varid,nf_double,                     &
     &                 2,CostDim,Aval,Vinfo,ncname)
!
!  Define free-surface.
!
        Vinfo( 1)=Vname(1,idFsur)
        Vinfo( 2)=Vname(2,idFsur)
        Vinfo( 3)=Vname(3,idFsur)
        Vinfo(14)=Vname(4,idFsur)
        Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
# endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idFsur,ng),r8)
        status=def_var(ncINIid(ng),iniVid(idFsur,ng),NF_FOUT,           &
     &                 nvd3,t2dgrd,Aval,Vinfo,ncname)
!
!  Define 2D momentum in the XI-direction.
!
        Vinfo( 1)=Vname(1,idUbar)
        Vinfo( 2)=Vname(2,idUbar)
        Vinfo( 3)=Vname(3,idUbar)
        Vinfo(14)=Vname(4,idUbar)
        Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_u'
# endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idUbar,ng),r8)
        status=def_var(ncINIid(ng),iniVid(idUbar,ng),NF_FOUT,           &
     &                 nvd3,u2dgrd,Aval,Vinfo,ncname)
!
!  Define 2D momentum in the ETA-direction.
!
        Vinfo( 1)=Vname(1,idVbar)
        Vinfo( 2)=Vname(2,idVbar)
        Vinfo( 3)=Vname(3,idVbar)
        Vinfo(14)=Vname(4,idVbar)
        Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_v'
# endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idVbar,ng),r8)
        status=def_var(ncINIid(ng),iniVid(idVbar,ng),NF_FOUT,           &
     &                 nvd3,v2dgrd,Aval,Vinfo,ncname)
# ifdef SOLVE3D
!
!  Define 3D momentum component in the XI-direction.
!
        Vinfo( 1)=Vname(1,idUvel)
        Vinfo( 2)=Vname(2,idUvel)
        Vinfo( 3)=Vname(3,idUvel)
        Vinfo(14)=Vname(4,idUvel)
        Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_u'
#  endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idUvel,ng),r8)
        status=def_var(ncINIid(ng),iniVid(idUvel,ng),NF_FOUT,           &
     &                 nvd4,u3dgrd,Aval,Vinfo,ncname)
!
!  Define 3D momentum component in the ETA-direction.
!
        Vinfo( 1)=Vname(1,idVvel)
        Vinfo( 2)=Vname(2,idVvel)
        Vinfo( 3)=Vname(3,idVvel)
        Vinfo(14)=Vname(4,idVvel)
        Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_v'
#  endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idVvel,ng),r8)
        status=def_var(ncINIid(ng),iniVid(idVvel,ng),NF_FOUT,           &
     &                 nvd4,v3dgrd,Aval,Vinfo,ncname)
!
!  Define tracer type variables.
!
        DO itrc=1,NT(ng)
          Vinfo( 1)=Vname(1,idTvar(itrc))
          Vinfo( 2)=Vname(2,idTvar(itrc))
          Vinfo( 3)=Vname(3,idTvar(itrc))
          Vinfo(14)=Vname(4,idTvar(itrc))
          Vinfo(16)=Vname(1,idtime)
#  ifdef SEDIMENT
            DO i=1,NST
              IF (itrc.eq.idsed(i)) THEN
                WRITE (Vinfo(19),20) 1000.0_r8*Sd50(i,ng)
              END IF
            END DO
#  endif
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r3dvar,r8)
          status=def_var(ncINIid(ng),iniTid(itrc,ng),NF_FOUT,           &
     &                   nvd4,t3dgrd,Aval,Vinfo,ncname)
        END DO

# endif
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
        status=nf_enddef(ncINIid(ng))
      END IF
!
!-----------------------------------------------------------------------
!  Write out time-recordless, information variables.  Deactive file
!  creation switch.
!-----------------------------------------------------------------------
!
      IF (LdefINI(ng)) THEN
        status=wrt_info(ng, iNLM, ncINIid(ng), OutThread, ncname)
        IF (exit_flag.ne.NoError) RETURN
        LdefINI(ng)=.FALSE.
      END IF
!
!=======================================================================
!  Open an existing initialization file, check its contents, and
!  prepare for appending data.
!=======================================================================
!
      IF (.not.LdefINI(ng).and.OutThread) THEN
!
!  Inquire about the contents of initialization NetCDF file:  Inquire
!  about the dimensions and variables.  Check for consistency.
!
        ncname=INIname(ng)
        ifield=0
        CALL opencdf (ng, 1, ncname, fname, N(ng), ifield, nrec)
        IF (exit_flag.ne.NoError) RETURN
!
!  Open initialization file for read/write.
!
        status=nf_open(TRIM(ncname),nf_write,ncINIid(ng))
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,40) TRIM(ncname)
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
!
!  Scan variable list from input NetCDF and activate switches for
!  initialization variables. Get variable IDs.
!
        DO i=1,nvars
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idtime))) THEN
            got_var(idtime)=.TRUE.
            status=nf_inq_varid(ncINIid(ng),TRIM(Vname(1,idtime)),      &
     &                          iniVid(idtime,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idFsur))) THEN
            got_var(idFsur)=.TRUE.
            status=nf_inq_varid(ncINIid(ng),TRIM(Vname(1,idFsur)),      &
     &                          iniVid(idFsur,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idUbar))) THEN
            got_var(idUbar)=.TRUE.
            status=nf_inq_varid(ncINIid(ng),TRIM(Vname(1,idUbar)),      &
     &                          iniVid(idUbar,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVbar))) THEN
            got_var(idVbar)=.TRUE.
            status=nf_inq_varid(ncINIid(ng),TRIM(Vname(1,idVbar)),      &
     &                          iniVid(idVbar,ng))
          END IF
# ifdef SOLVE3D
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idUvel))) THEN
            got_var(idUvel)=.TRUE.
            status=nf_inq_varid(ncINIid(ng),TRIM(Vname(1,idUvel)),      &
     &                          iniVid(idUvel,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVvel))) THEN
            got_var(idVvel)=.TRUE.
            status=nf_inq_varid(ncINIid(ng),TRIM(Vname(1,idVvel)),      &
     &                          iniVid(idVvel,ng))
          END IF
          DO itrc=1,NT(ng)
            IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idTvar(itrc)))) THEN
             got_var(idTvar(itrc))=.TRUE.
             status=nf_inq_varid(ncINIid(ng),                           &
     &                           TRIM(Vname(1,idTvar(itrc))),           &
     &                           iniTid(itrc,ng))
            END IF
          END DO
# endif
        END DO
!
!  Check if initialization variables are available in input NetCDF
!  file.
!
        IF (.not.got_var(idtime)) THEN
          WRITE (stdout,40) TRIM(Vname(1,idtime)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idFsur)) THEN
          WRITE (stdout,40) TRIM(Vname(1,idFsur)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbar)) THEN
          WRITE (stdout,40) TRIM(Vname(1,idUbar)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbar)) THEN
          WRITE (stdout,40) TRIM(Vname(1,idVbar)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# ifdef SOLVE3D
        IF (.not.got_var(idUvel)) THEN
          WRITE (stdout,40) TRIM(Vname(1,idUvel)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVvel)) THEN
          WRITE (stdout,40) TRIM(Vname(1,idVvel)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        DO itrc=1,NT(ng)
          IF (.not.got_var(idTvar(itrc))) THEN
            WRITE (stdout,40) TRIM(Vname(1,idTvar(itrc))), TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
# endif
!
!  Set unlimited time record dimension to the appropriate value.
!
        tINIindx(ng)=tsize
      END IF
!
  10  FORMAT (/,' DEF_INI - unable to create initial NetCDF file: ',a)
  20  FORMAT (1pe11.4,1x,'millimeter')
  30  FORMAT (/,' DEF_INI - unable to open initial NetCDF file: ',a)
  40  FORMAT (/,' DEF_INI - unable to find variable: ',a,2x,            &
     &        ' in initial NetCDF file: ',a)
      RETURN
      END SUBROUTINE def_ini
#else
      SUBROUTINE def_ini
      RETURN
      END SUBROUTINE def_ini
#endif
