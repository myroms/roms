#include "cppdefs.h"
      SUBROUTINE def_rst (ng)
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine creates restart NetCDF file, it defines its            !
!  dimensions, attributes, and variables.                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
#ifdef SEDIMENT
      USE mod_sediment
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
      integer, parameter :: Natt = 20

      logical :: Ldefine, got_var(NV)

      integer :: i, j, ifield, nrec, nvd3, nvd4
      integer :: recdim, status

      integer :: DimIDs(24), t2dgrd(3), u2dgrd(3), v2dgrd(3)

      integer :: def_info, def_var, wrt_info

#ifdef SOLVE3D
      logical :: got_trc(NT(ng))

      integer :: itrc

      integer :: b3dgrd(4), t3dgrd(4), u3dgrd(4), v3dgrd(4), w3dgrd(4)
#endif

      real(r4) :: Aval(4)

      character (len=65) :: Vinfo(Natt)
      character (len=80) :: fname, ncname
!
!=======================================================================
!  Create a new restart NetCDF file.
!=======================================================================
!
!  Activate creation of restart NetCDF file.  Create a new restart
!  file if during a restart run, the restart filename "rstname"
!  is different than the initial filename "ininame".
!
      IF (exit_flag.ne.0) RETURN
      ncname=rstname(ng)
      Ldefine=.false.
      IF (((nrrec.eq.0).and.(iic(ng).eq.ntstart)).or.                   &
     &    ((nrrec.ne.0).and.(TRIM(ncname).ne.TRIM(ininame(ng))))) THEN
        Ldefine=.true.
      END IF
!
      IF (Ldefine.and.OutThread) THEN
        status=nf_create(TRIM(ncname),nf_clobber,ncrstid(ng))
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,10) TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
      END IF
#ifdef DISTRIBUTE
      CALL mp_bcasti (ng, ncrstid(ng), 1)
#endif
!
!-----------------------------------------------------------------------
!  Define the dimensions of staggered fields.
!-----------------------------------------------------------------------
!
      IF (Ldefine.and.OutThread) THEN
        status=nf_def_dim(ncrstid(ng),'xi_rho'  ,Lm(ng)+2,              &
     &                    DimIDs( 1))
        status=nf_def_dim(ncrstid(ng),'xi_u'    ,Lm(ng)+1,              &
     &                    DimIDs( 2))
        status=nf_def_dim(ncrstid(ng),'xi_v'    ,Lm(ng)+2,              &
     &                    DimIDs( 3))
        status=nf_def_dim(ncrstid(ng),'xi_psi'  ,Lm(ng)+1,              &
     &                    DimIDs( 4))
        status=nf_def_dim(ncrstid(ng),'eta_rho' ,Mm(ng)+2,              &
     &                    DimIDs( 5))
        status=nf_def_dim(ncrstid(ng),'eta_u'   ,Mm(ng)+2,              &
     &                    DimIDs( 6))
        status=nf_def_dim(ncrstid(ng),'eta_v'   ,Mm(ng)+1,              &
     &                    DimIDs( 7))
        status=nf_def_dim(ncrstid(ng),'eta_psi' ,Mm(ng)+1,              &
     &                    DimIDs( 8))
#if defined WRITE_WATER && defined MASKING
        status=nf_def_dim(ncrstid(ng),'xy_rho'  ,Nxyr(ng),              &
     &                    DimIDs(17))
        status=nf_def_dim(ncrstid(ng),'xy_u'    ,Nxyu(ng),              &
     &                    DimIDs(18))
        status=nf_def_dim(ncrstid(ng),'xy_v'    ,Nxyv(ng),              &
     &                    DimIDs(19))
#endif
#ifdef SOLVE3D
# if defined WRITE_WATER && defined MASKING
        status=nf_def_dim(ncrstid(ng),'xyz_rho' ,Nxyr(ng)*N(ng),        &
     &                    DimIDs(20))
        status=nf_def_dim(ncrstid(ng),'xyz_u'   ,Nxyu(ng)*N(ng),        &
     &                    DimIDs(21))
        status=nf_def_dim(ncrstid(ng),'xyz_v'   ,Nxyv(ng)*N(ng),        &
     &                    DimIDs(22))
        status=nf_def_dim(ncrstid(ng),'xyz_w'   ,Nxyr(ng)*(N(ng)+1),    &
     &                    DimIDs(23))
# endif
        status=nf_def_dim(ncrstid(ng),'s_rho'   ,N(ng),                 &
     &                    DimIDs( 9))
        status=nf_def_dim(ncrstid(ng),'s_w'     ,N(ng)+1,               &
     &                    DimIDs(10))
        status=nf_def_dim(ncrstid(ng),'tracer'  ,NT(ng),                &
     &                    DimIDs(11))
# ifdef SEDIMENT
        status=nf_def_dim(ncrstid(ng),'Nbed'    ,Nbed,                  &
     &                    DimIDs(16))
#  if defined WRITE_WATER && defined MASKING
        status=nf_def_dim(ncrstid(ng),'xybed'   ,Nxyr(ng)*Nbed,         &
     &                    DimIDs(24))
#  endif
# endif
#endif
        status=nf_def_dim(ncrstid(ng),'boundary',4,                     &
     &                    DimIDs(14))
        status=nf_def_dim(ncrstid(ng),TRIM(ADJUSTL(Vname(5,idtime))),   &
     &                    nf_unlimited,DimIDs(12))
        recdim=DimIDs(12)
!
!  Set number of dimensions for output variables.
!
#if defined WRITE_WATER && defined MASKING
        nvd3=2
        nvd4=2
#else
        nvd3=3
        nvd4=4
#endif
!
!  Define dimension vectors for staggered tracer type variables.
!
#if defined WRITE_WATER && defined MASKING
        t2dgrd(1)=DimIDs(17)
        t2dgrd(2)=DimIDs(12)
# ifdef SOLVE3D
        t3dgrd(1)=DimIDs(20)
        t3dgrd(2)=DimIDs(12)
# endif
#else
        t2dgrd(1)=DimIDs( 1)
        t2dgrd(2)=DimIDs( 5)
        t2dgrd(3)=DimIDs(12)
# ifdef SOLVE3D
        t3dgrd(1)=DimIDs( 1)
        t3dgrd(2)=DimIDs( 5)
        t3dgrd(3)=DimIDs( 9)
        t3dgrd(4)=DimIDs(12)
# endif
#endif
!
!  Define dimension vectors for staggered u-momemtum type variables.
!
#if defined WRITE_WATER && defined MASKING
        u2dgrd(1)=DimIDs(18)
        u2dgrd(2)=DimIDs(12)
# ifdef SOLVE3D
        u3dgrd(1)=DimIDs(21)
        u3dgrd(2)=DimIDs(12)
# endif
#else
        u2dgrd(1)=DimIDs( 2)
        u2dgrd(2)=DimIDs( 6)
        u2dgrd(3)=DimIDs(12)
# ifdef SOLVE3D
        u3dgrd(1)=DimIDs( 2)
        u3dgrd(2)=DimIDs( 6)
        u3dgrd(3)=DimIDs( 9)
        u3dgrd(4)=DimIDs(12)
# endif
#endif
!
!  Define dimension vectors for staggered v-momemtum type variables.
!
#if defined WRITE_WATER && defined MASKING
        v2dgrd(1)=DimIDs(19)
        v2dgrd(2)=DimIDs(12)
# ifdef SOLVE3D
        v3dgrd(1)=DimIDs(22)
        v3dgrd(2)=DimIDs(12)
# endif
#else
        v2dgrd(1)=DimIDs( 3)
        v2dgrd(2)=DimIDs( 7)
        v2dgrd(3)=DimIDs(12)
# ifdef SOLVE3D
        v3dgrd(1)=DimIDs( 3)
        v3dgrd(2)=DimIDs( 7)
        v3dgrd(3)=DimIDs( 9)
        v3dgrd(4)=DimIDs(12)
# endif
#endif
#ifdef SOLVE3D
!
!  Define dimension vector for staggered w-momemtum type variables.
!
# if defined WRITE_WATER && defined MASKING
        w3dgrd(1)=DimIDs(23)
        w3dgrd(2)=DimIDs(12)
# else
        w3dgrd(1)=DimIDs( 1)
        w3dgrd(2)=DimIDs( 5)
        w3dgrd(3)=DimIDs(10)
        w3dgrd(4)=DimIDs(12)
# endif
!
!  Define dimension vector for sediment bed layer type variables.
!
# if defined WRITE_WATER && defined MASKING
        b3dgrd(1)=DimIDs(24)
        b3dgrd(2)=DimIDs(12)
# else
        b3dgrd(1)=DimIDs( 1)
        b3dgrd(2)=DimIDs( 5)
        b3dgrd(3)=DimIDs(16)
        b3dgrd(4)=DimIDs(12)
# endif
#endif
!
!  Initialize unlimited time record dimension.
!
        trstindx(ng)=0
!
!  Initialize local information variable arrays.
!
        DO i=1,Natt
          DO j=1,65
            Vinfo(i)=' '
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
        status=def_info(ng,ncrstid(ng),ncname,DimIDs)
        IF (exit_flag.ne.0) RETURN
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
          Vinfo( 3)='seconds since 0001-01-01 00:00:00'
          Vinfo( 4)='360 days in every year, 30 days in each month'
        ELSE IF (INT(time_ref).eq.0) THEN
          Vinfo( 3)='seconds since 0001-01-01 00:00:00'
          Vinfo( 4)='365.25 days in every year'
        ELSE IF (time_ref.gt.0.0_r8) THEN
          Vinfo( 3)='seconds since '//r_text
        END IF
        Vinfo(14)=Vname(4,idtime)
        status=def_var(ncrstid(ng),rstVid(idtime,ng),NF_TYPE,           &
     &                 1,recdim,Aval,Vinfo,ncname)
!
!  Define free-surface.
!
        Vinfo( 1)=Vname(1,idFsur)
        Vinfo( 2)=Vname(2,idFsur)
        Vinfo( 3)=Vname(3,idFsur)
        Vinfo(14)=Vname(4,idFsur)
        Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
#endif
        status=def_var(ncrstid(ng),rstVid(idFsur,ng),NF_FRST,           &
     &                 nvd3,t2dgrd,Aval,Vinfo,ncname)
!
!  Define 2D momentum in the XI-direction.
!
        Vinfo( 1)=Vname(1,idUbar)
        Vinfo( 2)=Vname(2,idUbar)
        Vinfo( 3)=Vname(3,idUbar)
        Vinfo(14)=Vname(4,idUbar)
        Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_u'
#endif
        status=def_var(ncrstid(ng),rstVid(idUbar,ng),NF_FRST,           &
     &                 nvd3,u2dgrd,Aval,Vinfo,ncname)
!
!  Define 2D momentum in the ETA-direction.
!
        Vinfo( 1)=Vname(1,idVbar)
        Vinfo( 2)=Vname(2,idVbar)
        Vinfo( 3)=Vname(3,idVbar)
        Vinfo(14)=Vname(4,idVbar)
        Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_v'
#endif
        status=def_var(ncrstid(ng),rstVid(idVbar,ng),NF_FRST,           &
     &                 nvd3,v2dgrd,Aval,Vinfo,ncname)
#ifdef SOLVE3D
!
!  Define 3D momentum component in the XI-direction.
!
        Vinfo( 1)=Vname(1,idUvel)
        Vinfo( 2)=Vname(2,idUvel)
        Vinfo( 3)=Vname(3,idUvel)
        Vinfo(14)=Vname(4,idUvel)
        Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_u'
# endif
        status=def_var(ncrstid(ng),rstVid(idUvel,ng),NF_FRST,           &
     &                 nvd4,u3dgrd,Aval,Vinfo,ncname)
!
!  Define 3D momentum component in the ETA-direction.
!
        Vinfo( 1)=Vname(1,idVvel)
        Vinfo( 2)=Vname(2,idVvel)
        Vinfo( 3)=Vname(3,idVvel)
        Vinfo(14)=Vname(4,idVvel)
        Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_v'
# endif
        status=def_var(ncrstid(ng),rstVid(idVvel,ng),NF_FRST,           &
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
# ifdef SEDIMENT
            DO i=1,NST
              IF (itrc.eq.idsed(i)) THEN
                WRITE (Vinfo(19),20) 1000.0_r8*Sd50(i,ng)
              END IF
            END DO
# endif
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          status=def_var(ncrstid(ng),rstTid(itrc,ng),NF_FRST,           &
     &                   nvd4,t3dgrd,Aval,Vinfo,ncname)
        END DO
!
!  Define density anomaly.
!
        Vinfo( 1)=Vname(1,idDano)
        Vinfo( 2)=Vname(2,idDano)
        Vinfo( 3)=Vname(3,idDano)
        Vinfo(14)=Vname(4,idDano)
        Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
# endif
        status=def_var(ncrstid(ng),rstVid(idDano,ng),NF_FRST,           &
     &                 nvd4,t3dgrd,Aval,Vinfo,ncname)
# ifdef LMD_SKPP
!
!  Define depth of surface boundary layer.
!
        Vinfo( 1)=Vname(1,idHsbl)
        Vinfo( 2)=Vname(2,idHsbl)
        Vinfo( 3)=Vname(3,idHsbl)
        Vinfo(14)=Vname(4,idHsbl)
        Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
#  endif
        status=def_var(ncrstid(ng),rstVid(idHsbl,ng),NF_FRST,           &
     &                 nvd3,t2dgrd,Aval,Vinfo,ncname)
# endif
# ifdef LMD_BKPP
!
!  Define depth of bottom boundary layer.
!
        Vinfo( 1)=Vname(1,idHbbl)
        Vinfo( 2)=Vname(2,idHbbl)
        Vinfo( 3)=Vname(3,idHbbl)
        Vinfo(14)=Vname(4,idHbbl)
        Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
#  endif
        status=def_var(ncrstid(ng),rstVid(idHbbl,ng),NF_FRST,           &
     &                 nvd3,t2dgrd,Aval,Vinfo,ncname)
# endif
# if defined BVF_MIXING || defined LMD_MIXING || \
     defined GLS_MIXING || defined MY25_MIXING
!
!  Define vertical viscosity coefficient.
!
        Vinfo( 1)=Vname(1,idVvis)
        Vinfo( 2)=Vname(2,idVvis)
        Vinfo( 3)=Vname(3,idVvis)
        Vinfo(14)=Vname(4,idVvis)
        Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
#  endif
        status=def_var(ncrstid(ng),rstVid(idVvis,ng),NF_FRST,           &
     &                 nvd4,w3dgrd,Aval,Vinfo,ncname)
!
!  Define vertical diffusion coefficient for potential temperature.
!
        Vinfo( 1)=Vname(1,idTdif)
        Vinfo( 2)=Vname(2,idTdif)
        Vinfo( 3)=Vname(3,idTdif)
        Vinfo(14)=Vname(4,idTdif)
        Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
#  endif
        status=def_var(ncrstid(ng),rstVid(idTdif,ng),NF_FRST,           &
     &                 nvd4,w3dgrd,Aval,Vinfo,ncname)
#  ifdef SALINITY
!
!  Define vertical diffusion coefficient for salinity.
!
        Vinfo( 1)=Vname(1,idSdif)
        Vinfo( 2)=Vname(2,idSdif)
        Vinfo( 3)=Vname(3,idSdif)
        Vinfo(14)=Vname(4,idSdif)
        Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
#   endif
        status=def_var(ncrstid(ng),rstVid(idSdif,ng),NF_FRST,           &
     &                 nvd4,w3dgrd,Aval,Vinfo,ncname)
#  endif
# endif
# ifdef SEDIMENT
!
!  Define bed layer sediment fraction for each size class.
!
        DO i=1,NST
          Vinfo( 1)=Vname(1,idfrac(i))
          Vinfo( 2)=Vname(2,idfrac(i))
          Vinfo( 3)=Vname(3,idfrac(i))
          Vinfo(14)=Vname(4,idfrac(i))
          Vinfo(16)=Vname(1,idtime)
          WRITE (Vinfo(19),20) 1000.0_r8*Sd50(i,ng)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          status=def_var(ncrstid(ng),rstVid(idfrac(i),ng),NF_FRST,      &
     &                   nvd4,b3dgrd,Aval,Vinfo,ncname)
        END DO
!
!  Define bed layer sediment properties.
!
        DO i=1,MBEDP
          Vinfo( 1)=Vname(1,idSbed(i))
          Vinfo( 2)=Vname(2,idSbed(i))
          Vinfo( 3)=Vname(3,idSbed(i))
          Vinfo(14)=Vname(4,idSbed(i))
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          status=def_var(ncrstid(ng),rstVid(idSbed(i),ng),NF_FRST,      &
     &                   nvd4,b3dgrd,Aval,Vinfo,ncname)
        END DO
# endif
#endif
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
        status=nf_enddef(ncrstid(ng))
      END IF
!
!-----------------------------------------------------------------------
!  Write out time-recordless, information variables.
!-----------------------------------------------------------------------
!
      IF (Ldefine) THEN
        status=wrt_info(ng,ncrstid(ng),OutThread,ncname)
        IF (exit_flag.ne.0) RETURN
      END IF
!
!=======================================================================
!  Open an existing restart file, check its contents, and prepare for
!  appending data.
!=======================================================================
!
      IF (.not.Ldefine.and.OutThread) THEN
!
!  Inquire about the contents of restart NetCDF file:  Inquire about
!  the dimensions and variables.  Check for consistency.
!
        ncname=rstname(ng)
        ifield=0
        CALL opencdf (ng,1,ncname,fname,N(ng),ifield,nrec)
        if (exit_flag.ne.0) RETURN
!
!  Open restart file for read/write.
!
        status=nf_open(TRIM(ncname),nf_write,ncrstid(ng))
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,40) TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
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
!  restart variables. Get variable IDs.
!
        DO i=1,nvars
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idtime))) THEN
            got_var(idtime)=.true.
            status=nf_inq_varid(ncrstid(ng),TRIM(Vname(1,idtime)),      &
     &                          rstVid(idtime,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idFsur))) THEN
            got_var(idFsur)=.true.
            status=nf_inq_varid(ncrstid(ng),TRIM(Vname(1,idFsur)),      &
     &                          rstVid(idFsur,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idUbar))) THEN
            got_var(idUbar)=.true.
            status=nf_inq_varid(ncrstid(ng),TRIM(Vname(1,idUbar)),      &
     &                          rstVid(idUbar,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVbar))) THEN
            got_var(idVbar)=.true.
            status=nf_inq_varid(ncrstid(ng),TRIM(Vname(1,idVbar)),      &
     &                          rstVid(idVbar,ng))
          END IF
#ifdef SOLVE3D
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idUvel))) THEN
            got_var(idUvel)=.true.
            status=nf_inq_varid(ncrstid(ng),TRIM(Vname(1,idUvel)),      &
     &                          rstVid(idUvel,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVvel))) THEN
            got_var(idVvel)=.true.
            status=nf_inq_varid(ncrstid(ng),TRIM(Vname(1,idVvel)),      &
     &                          rstVid(idVvel,ng))
          END IF
          DO itrc=1,NT(ng)
            IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idTvar(itrc)))) THEN
             got_trc(itrc)=.true.
             status=nf_inq_varid(ncrstid(ng),                           &
     &                           TRIM(Vname(1,idTvar(itrc))),           &
     &                           rstTid(itrc,ng))
            END IF
          END DO
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idDano))) THEN
            got_var(idDano)=.true.
            status=nf_inq_varid(ncrstid(ng),TRIM(Vname(1,idDano)),      &
     &                          rstVid(idDano,ng))
          END IF
# ifdef LMD_SKPP
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idHsbl))) THEN
            got_var(idHsbl)=.true.
            status=nf_inq_varid(ncrstid(ng),TRIM(Vname(1,idHsbl)),      &
     &                          rstVid(idHsbl,ng))
          END IF
# endif
# ifdef LMD_BKPP
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idHbbl))) THEN
            got_var(idHbbl)=.true.
            status=nf_inq_varid(ncrstid(ng),TRIM(Vname(1,idHbbl)),      &
     &                          rstVid(idHbbl,ng))
          END IF
# endif
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVvis))) THEN
            got_var(idVvis)=.true.
            status=nf_inq_varid(ncrstid(ng),TRIM(Vname(1,idVvis)),      &
     &                          rstVid(idVvis,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idTdif))) THEN
            got_var(idTdif)=.true.
            status=nf_inq_varid(ncrstid(ng),TRIM(Vname(1,idTdif)),      &
     &                          rstVid(idTdif,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idSdif))) THEN
            got_var(idSdif)=.true.
            status=nf_inq_varid(ncrstid(ng),TRIM(Vname(1,idSdif)),      &
     &                          rstVid(idSdif,ng))
          END IF
# ifdef SEDIMENT
          DO itrc=1,NST
            IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idfrac(itrc)))) THEN
              got_var(idfrac(itrc))=.true.
              status=nf_inq_varid(ncrstid(ng),                          &
     &                            TRIM(Vname(1,idfrac(itrc))),          &
     &                            rstVid(idfrac(itrc),ng))
            END IF
          END DO
          DO itrc=1,MBEDP
            IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idSbed(itrc)))) THEN
              got_var(idSbed(itrc))=.true.
              status=nf_inq_varid(ncrstid(ng),                          &
     &                            TRIM(Vname(1,idSbed(itrc))),          &
     &                            rstVid(idSbed(itrc),ng))
            END IF
          END DO
# endif
#endif
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
#ifdef SOLVE3D
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
          IF (.not.got_trc(itrc)) THEN
            WRITE (stdout,40) TRIM(Vname(1,idTvar(itrc))), TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
        IF (.not.got_var(idDano)) THEN
          WRITE (stdout,40) TRIM(Vname(1,idDano)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# ifdef SEDIMENT
        DO i=1,NST
          IF (.not.got_var(idfrac(i))) THEN
            WRITE (stdout,40) TRIM(Vname(1,idfrac(i))), TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
        DO i=1,MBEDP
          IF (.not.got_var(idSbed(i))) THEN
            WRITE (stdout,40) TRIM(Vname(1,idSbed(i))), TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
# endif
#endif
!
!  Set unlimited time record dimension to current value.
!
        IF (lcycle(ng)) THEN
          trstindx(ng)=0
        ELSE
          trstindx(ng)=tsize
        END IF
      END IF
!
  10  format(/,' DEF_RST - unable to create restart NetCDF file: ',a)
  20  format(1pe11.4,1x,'millimeter')
  30  format(/,' DEF_RST - unable to open restart NetCDF file: ',a)
  40  format(/,' DEF_RST - unable to find variable: ',a,2x,             &
     &         ' in restart NetCDF file: ',a)
      RETURN
      END SUBROUTINE def_rst
