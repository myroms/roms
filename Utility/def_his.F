#include "cppdefs.h"
      SUBROUTINE def_his (ng,ldef)
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine creates history NetCDF file, it defines its            !
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

      logical, intent(in) :: ldef
!
!  Local variable declarations.
!
      logical :: got_var(NV)

      integer, parameter :: Natt = 20

      integer :: i, j, ifield, nrec, nvd3, nvd4
      integer :: recdim, status

      integer :: DimIDs(24), t2dgrd(3), u2dgrd(3), v2dgrd(3)

      integer :: def_info, def_var, wrt_info

#ifdef SOLVE3D
      logical :: got_trc(NT(ng))

      integer :: itrc

# ifdef SEDIMENT
      integer :: b3dgrd(4)
# endif
      integer :: t3dgrd(4), u3dgrd(4), v3dgrd(4), w3dgrd(4)
#endif

      real(r8) :: Aval(4)

      character (len=65) :: Vinfo(Natt)
      character (len=80) :: fname, ncname
!
!-----------------------------------------------------------------------
!  Set and report file name.
!-----------------------------------------------------------------------
!
      IF (exit_flag.ne.0) RETURN
      ncname=hisname(ng)
!
      IF (Master) THEN
        IF (ldef) THEN
          WRITE (stdout,10) TRIM(ncname)
        ELSE
          WRITE (stdout,20) TRIM(ncname)
        END IF
      END IF
!
!=======================================================================
!  Create a new history file.
!=======================================================================
!
      IF (ldef.and.OutThread) THEN
        status=nf_create(TRIM(ncname),nf_clobber,nchisid(ng))
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,30) TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
      END IF
#ifdef DISTRIBUTE
      CALL mp_bcasti (ng, nchisid(ng), 1)
#endif
!
!-----------------------------------------------------------------------
!  Define the dimensions of staggered fields.
!-----------------------------------------------------------------------
!
      IF (ldef.and.OutThread) THEN
        status=nf_def_dim(nchisid(ng),'xi_rho'  ,Lm(ng)+2,              &
     &                    DimIDs( 1))
        status=nf_def_dim(nchisid(ng),'xi_u'    ,Lm(ng)+1,              &
     &                    DimIDs( 2))
        status=nf_def_dim(nchisid(ng),'xi_v'    ,Lm(ng)+2,              &
     &                    DimIDs( 3))
        status=nf_def_dim(nchisid(ng),'xi_psi'  ,Lm(ng)+1,              &
     &                    DimIDs( 4))
        status=nf_def_dim(nchisid(ng),'eta_rho' ,Mm(ng)+2,              &
     &                    DimIDs( 5))
        status=nf_def_dim(nchisid(ng),'eta_u'   ,Mm(ng)+2,              &
     &                    DimIDs( 6))
        status=nf_def_dim(nchisid(ng),'eta_v'   ,Mm(ng)+1,              &
     &                    DimIDs( 7))
        status=nf_def_dim(nchisid(ng),'eta_psi' ,Mm(ng)+1,              &
     &                    DimIDs( 8))
#if defined WRITE_WATER && defined MASKING
        status=nf_def_dim(nchisid(ng),'xy_rho'  ,Nxyr(ng),              &
     &                    DimIDs(17))
        status=nf_def_dim(nchisid(ng),'xy_u'    ,Nxyu(ng),              &
     &                    DimIDs(18))
        status=nf_def_dim(nchisid(ng),'xy_v'    ,Nxyv(ng),              &
     &                    DimIDs(19))
#endif
#ifdef SOLVE3D
# if defined WRITE_WATER && defined MASKING
        status=nf_def_dim(nchisid(ng),'xyz_rho' ,Nxyr(ng)*N(ng),        &
     &                    DimIDs(20))
        status=nf_def_dim(nchisid(ng),'xyz_u'   ,Nxyu(ng)*N(ng),        &
     &                    DimIDs(21))
        status=nf_def_dim(nchisid(ng),'xyz_v'   ,Nxyv(ng)*N(ng),        &
     &                    DimIDs(22))
        status=nf_def_dim(nchisid(ng),'xyz_w'   ,Nxyr(ng)*(N(ng)+1),    &
     &                    DimIDs(23))
# endif
        status=nf_def_dim(nchisid(ng),'N'       ,N(ng),                 &
     &                    DimIDs( 9))
        status=nf_def_dim(nchisid(ng),'s_rho'   ,N(ng),                 &
     &                    DimIDs( 9))
        status=nf_def_dim(nchisid(ng),'s_w'     ,N(ng)+1,               &
     &                    DimIDs(10))
        status=nf_def_dim(nchisid(ng),'tracer'  ,NT(ng),                &
     &                    DimIDs(11))
# ifdef SEDIMENT
        status=nf_def_dim(nchisid(ng),'Nbed'    ,Nbed,                  &
     &                    DimIDs(16))
#  if defined WRITE_WATER && defined MASKING
        status=nf_def_dim(nchisid(ng),'xybed'   ,Nxyr(ng)*Nbed,         &
     &                    DimIDs(24))
#  endif
# endif
#endif
        status=nf_def_dim(nchisid(ng),'boundary',4,                     &
     &                    DimIDs(14))
        status=nf_def_dim(nchisid(ng),TRIM(ADJUSTL(Vname(5,idtime))),   &
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
# ifdef SEDIMENT
!
!  Define dimension vector for sediment bed layer type variables.
!
#  if defined WRITE_WATER && defined MASKING
        b3dgrd(1)=DimIDs(24)
        b3dgrd(2)=DimIDs(12)
#  else
        b3dgrd(1)=DimIDs( 1)
        b3dgrd(2)=DimIDs( 5)
        b3dgrd(3)=DimIDs(16)
        b3dgrd(4)=DimIDs(12)
#  endif
# endif
#endif
!
!  Initialize unlimited time record dimension.
!
        thisindx(ng)=0
!
!  Initialize local information variable arrays.
!
        DO i=1,Natt
          DO j=1,65
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
        status=def_info(ng,nchisid(ng),ncname,DimIDs)
        IF (exit_flag.ne.0) RETURN
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
        status=def_var(nchisid(ng),hisVid(idtime,ng),NF_TYPE,           &
     &                 1,recdim,Aval,Vinfo,ncname)
!
!  Define free-surface.
!
        IF (Hout(idFsur,ng)) THEN
          Vinfo( 1)=Vname(1,idFsur)
          Vinfo( 2)=Vname(2,idFsur)
          Vinfo( 3)=Vname(3,idFsur)
          Vinfo(14)=Vname(4,idFsur)
          Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#endif
          status=def_var(nchisid(ng),hisVid(idFsur,ng),NF_FOUT,         &
     &                   nvd3,t2dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define 2D momentum in the XI-direction.
!
        IF (Hout(idUbar,ng)) THEN
          Vinfo( 1)=Vname(1,idUbar)
          Vinfo( 2)=Vname(2,idUbar)
          Vinfo( 3)=Vname(3,idUbar)
          Vinfo(14)=Vname(4,idUbar)
          Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#endif
          status=def_var(nchisid(ng),hisVid(idUbar,ng),NF_FOUT,         &
     &                   nvd3,u2dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define 2D momentum in the ETA-direction.
!
        IF (Hout(idVbar,ng)) THEN
          Vinfo( 1)=Vname(1,idVbar)
          Vinfo( 2)=Vname(2,idVbar)
          Vinfo( 3)=Vname(3,idVbar)
          Vinfo(14)=Vname(4,idVbar)
          Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#endif
          status=def_var(nchisid(ng),hisVid(idVbar,ng),NF_FOUT,         &
     &                   nvd3,v2dgrd,Aval,Vinfo,ncname)
        END IF
#ifdef SOLVE3D
!
!  Define 3D momentum component in the XI-direction.
!
        IF (Hout(idUvel,ng)) THEN
          Vinfo( 1)=Vname(1,idUvel)
          Vinfo( 2)=Vname(2,idUvel)
          Vinfo( 3)=Vname(3,idUvel)
          Vinfo(14)=Vname(4,idUvel)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
# endif
          status=def_var(nchisid(ng),hisVid(idUvel,ng),NF_FOUT,         &
     &                   nvd4,u3dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define 3D momentum component in the ETA-direction.
!
        IF (Hout(idVvel,ng)) THEN
          Vinfo( 1)=Vname(1,idVvel)
          Vinfo( 2)=Vname(2,idVvel)
          Vinfo( 3)=Vname(3,idVvel)
          Vinfo(14)=Vname(4,idVvel)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
# endif
          status=def_var(nchisid(ng),hisVid(idVvel,ng),NF_FOUT,         &
     &                   nvd4,v3dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define 3D momentum component in the S-direction.
!
        IF (Hout(idWvel,ng)) THEN
          Vinfo( 1)=Vname(1,idWvel)
          Vinfo( 2)=Vname(2,idWvel)
          Vinfo( 3)=Vname(3,idWvel)
          Vinfo(14)=Vname(4,idWvel)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          status=def_var(nchisid(ng),hisVid(idWvel,ng),NF_FOUT,         &
     &                   nvd4,w3dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define S-coordinate vertical "omega" momentum component.
!
        IF (Hout(idOvel,ng)) THEN
          Vinfo( 1)=Vname(1,idOvel)
          Vinfo( 2)=Vname(2,idOvel)
          Vinfo( 3)=Vname(3,idOvel)
          Vinfo(14)=Vname(4,idOvel)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          status=def_var(nchisid(ng),hisVid(idOvel,ng),NF_FOUT,         &
     &                   nvd4,w3dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define tracer type variables.
!
        DO itrc=1,NT(ng)
          IF (Hout(idTvar(itrc),ng)) THEN
            Vinfo( 1)=Vname(1,idTvar(itrc))
            Vinfo( 2)=Vname(2,idTvar(itrc))
            Vinfo( 3)=Vname(3,idTvar(itrc))
            Vinfo(14)=Vname(4,idTvar(itrc))
            Vinfo(16)=Vname(1,idtime)
# ifdef SEDIMENT
            DO i=1,NST
              IF (itrc.eq.idsed(i)) THEN
                WRITE (Vinfo(19),40) 1000.0_r8*Sd50(i,ng)
              END IF
            END DO
# endif
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
            status=def_var(nchisid(ng),hisTid(itrc,ng),NF_FOUT,         &
     &                     nvd4,t3dgrd,Aval,Vinfo,ncname)
          END IF
        END DO
!
!  Define density anomaly.
!
        IF (Hout(idDano,ng)) THEN
          Vinfo( 1)=Vname(1,idDano)
          Vinfo( 2)=Vname(2,idDano)
          Vinfo( 3)=Vname(3,idDano)
          Vinfo(14)=Vname(4,idDano)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          status=def_var(nchisid(ng),hisVid(idDano,ng),NF_FOUT,         &
     &                   nvd4,t3dgrd,Aval,Vinfo,ncname)
        END IF
# ifdef LMD_SKPP
!
!  Define depth of surface boundary layer.
!
        IF (Hout(idHsbl,ng)) THEN
          Vinfo( 1)=Vname(1,idHsbl)
          Vinfo( 2)=Vname(2,idHsbl)
          Vinfo( 3)=Vname(3,idHsbl)
          Vinfo(14)=Vname(4,idHsbl)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          status=def_var(nchisid(ng),hisVid(idHsbl,ng),NF_FOUT,         &
     &                   nvd3,t2dgrd,Aval,Vinfo,ncname)
        END IF
# endif
# ifdef LMD_BKPP
!
!  Define depth of bottom boundary layer.
!
        IF (Hout(idHbbl,ng)) THEN
          Vinfo( 1)=Vname(1,idHbbl)
          Vinfo( 2)=Vname(2,idHbbl)
          Vinfo( 3)=Vname(3,idHbbl)
          Vinfo(14)=Vname(4,idHbbl)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          status=def_var(nchisid(ng),hisVid(idHbbl,ng),NF_FOUT,         &
     &                   nvd3,t2dgrd,Aval,Vinfo,ncname)
        END IF
# endif
!
!  Define vertical viscosity coefficient.
!
        IF (Hout(idVvis,ng)) THEN
          Vinfo( 1)=Vname(1,idVvis)
          Vinfo( 2)=Vname(2,idVvis)
          Vinfo( 3)=Vname(3,idVvis)
          Vinfo(14)=Vname(4,idVvis)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          status=def_var(nchisid(ng),hisVid(idVvis,ng),NF_FOUT,         &
     &                   nvd4,w3dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define vertical diffusion coefficient for potential temperature.
!
        IF (Hout(idTdif,ng)) THEN
          Vinfo( 1)=Vname(1,idTdif)
          Vinfo( 2)=Vname(2,idTdif)
          Vinfo( 3)=Vname(3,idTdif)
          Vinfo(14)=Vname(4,idTdif)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          status=def_var(nchisid(ng),hisVid(idTdif,ng),NF_FOUT,         &
     &                   nvd4,w3dgrd,Aval,Vinfo,ncname)
        END IF
# ifdef SALINITY
!
!  Define vertical diffusion coefficient for salinity.
!
        IF (Hout(idSdif,ng)) THEN
          Vinfo( 1)=Vname(1,idSdif)
          Vinfo( 2)=Vname(2,idSdif)
          Vinfo( 3)=Vname(3,idSdif)
          Vinfo(14)=Vname(4,idSdif)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          status=def_var(nchisid(ng),hisVid(idSdif,ng),NF_FOUT,         &
     &                   nvd4,w3dgrd,Aval,Vinfo,ncname)
        END IF
# endif
# if defined GLS_MIXING || defined MY25_MIXING
!
!  Define turbulent kinetic energy.
!
        IF (Hout(idMtke,ng)) THEN
          Vinfo( 1)=Vname(1,idMtke)
          Vinfo( 2)=Vname(2,idMtke)
          Vinfo( 3)=Vname(3,idMtke)
          Vinfo(14)=Vname(4,idMtke)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          status=def_var(nchisid(ng),hisVid(idMtke,ng),NF_FOUT,         &
     &                   nvd4,w3dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define turbulent kinetic energy time length scale.
!
        IF (Hout(idMtls,ng)) THEN
          Vinfo( 1)=Vname(1,idMtls)
          Vinfo( 2)=Vname(2,idMtls)
          Vinfo( 3)=Vname(3,idMtls)
          Vinfo(14)=Vname(4,idMtls)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          status=def_var(nchisid(ng),hisVid(idMtls,ng),NF_FOUT,         &
     &                   nvd4,w3dgrd,Aval,Vinfo,ncname)
        END IF
# endif
!
!  Define surface net heat flux.
!
        IF (Hout(idTsur(itemp),ng)) THEN
          Vinfo( 1)=Vname(1,idTsur(itemp))
          Vinfo( 2)=Vname(2,idTsur(itemp))
          Vinfo( 3)=Vname(3,idTsur(itemp))
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idTsur(itemp))
          Vinfo(16)=Vname(1,idtime)
          status=def_var(nchisid(ng),hisVid(idTsur(itemp),ng),          &
     &                   NF_FOUT,nvd3,t2dgrd,Aval,Vinfo,ncname)
        END IF
# ifdef BULK_FLUXES
!
!  Define latent heat flux.
!
        IF (Hout(idLhea,ng)) THEN
          Vinfo( 1)=Vname(1,idLhea)
          Vinfo( 2)=Vname(2,idLhea)
          Vinfo( 3)=Vname(3,idLhea)
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idLhea)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          status=def_var(nchisid(ng),hisVid(idLhea,ng),NF_FOUT,         &
     &                   nvd3,t2dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define sensible heat flux.
!
        IF (Hout(idShea,ng)) THEN
          Vinfo( 1)=Vname(1,idShea)
          Vinfo( 2)=Vname(2,idShea)
          Vinfo( 3)=Vname(3,idShea)
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idShea)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          status=def_var(nchisid(ng),hisVid(idShea,ng),NF_FOUT,         &
     &                   nvd3,t2dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define longwave radiation flux.
!
        IF (Hout(idLrad,ng)) THEN
          Vinfo( 1)=Vname(1,idLrad)
          Vinfo( 2)=Vname(2,idLrad)
          Vinfo( 3)=Vname(3,idLrad)
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idLrad)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          status=def_var(nchisid(ng),hisVid(idLrad,ng),NF_FOUT,         &
     &                   nvd3,t2dgrd,Aval,Vinfo,ncname)
        END IF
# endif
# ifdef SHORTWAVE
!
!  Define shortwave radiation flux.
!
        IF (Hout(idSrad,ng)) THEN
          Vinfo( 1)=Vname(1,idSrad)
          Vinfo( 2)=Vname(2,idSrad)
          Vinfo( 3)=Vname(3,idSrad)
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idSrad)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          status=def_var(nchisid(ng),hisVid(idSrad,ng),NF_FOUT,         &
     &                   nvd3,t2dgrd,Aval,Vinfo,ncname)
        END IF
# endif
#endif
!
!  Define surface U-momentum stress.
!
        IF (Hout(idUsms,ng)) THEN
          Vinfo( 1)=Vname(1,idUsms)
          Vinfo( 2)=Vname(2,idUsms)
          Vinfo( 3)=Vname(3,idUsms)
          Vinfo(14)=Vname(4,idUsms)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
# endif
          status=def_var(nchisid(ng),hisVid(idUsms,ng),NF_FOUT,         &
     &                   nvd3,u2dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define surface V-momentum stress.
!
        IF (Hout(idVsms,ng)) THEN
          Vinfo( 1)=Vname(1,idVsms)
          Vinfo( 2)=Vname(2,idVsms)
          Vinfo( 3)=Vname(3,idVsms)
          Vinfo(14)=Vname(4,idVsms)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
# endif
          status=def_var(nchisid(ng),hisVid(idVsms,ng),NF_FOUT,         &
     &                   nvd3,v2dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define bottom U-momentum stress.
!
        IF (Hout(idUbms,ng)) THEN
          Vinfo( 1)=Vname(1,idUbms)
          Vinfo( 2)=Vname(2,idUbms)
          Vinfo( 3)=Vname(3,idUbms)
          Vinfo(14)=Vname(4,idUbms)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
# endif
          status=def_var(nchisid(ng),hisVid(idUbms,ng),NF_FOUT,         &
     &                   nvd3,u2dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define bottom V-momentum stress.
!
        IF (Hout(idVbms,ng)) THEN
          Vinfo( 1)=Vname(1,idVbms)
          Vinfo( 2)=Vname(2,idVbms)
          Vinfo( 3)=Vname(3,idVbms)
          Vinfo(14)=Vname(4,idVbms)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
# endif
          status=def_var(nchisid(ng),hisVid(idVbms,ng),NF_FOUT,         &
     &                   nvd3,v2dgrd,Aval,Vinfo,ncname)
        END IF
#ifdef BBL_MODEL
!
!  Define wind-induced, bottom U-wave stress.
!
        IF (Hout(idUbws,ng)) THEN
          Vinfo( 1)=Vname(1,idUbws)
          Vinfo( 2)=Vname(2,idUbws)
          Vinfo( 3)=Vname(3,idUbws)
          Vinfo(14)=Vname(4,idUbws)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
# endif
          status=def_var(nchisid(ng),hisVid(idUbws,ng),NF_FOUT,         &
     &                   nvd3,u2dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define bottom wind-induced, bottom V-wave stress.
!
        IF (Hout(idVbws,ng)) THEN
          Vinfo( 1)=Vname(1,idVbws)
          Vinfo( 2)=Vname(2,idVbws)
          Vinfo( 3)=Vname(3,idVbws)
          Vinfo(14)=Vname(4,idVbws)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
# endif
          status=def_var(nchisid(ng),hisVid(idVbws,ng),NF_FOUT,         &
     &                   nvd3,v2dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define wind-induced, bed wave excursion amplitude.
!
        IF (Hout(idAbed,ng)) THEN
          Vinfo( 1)=Vname(1,idAbed)
          Vinfo( 2)=Vname(2,idAbed)
          Vinfo( 3)=Vname(3,idAbed)
          Vinfo(14)=Vname(4,idAbed)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          status=def_var(nchisid(ng),hisVid(idAbed,ng),NF_FOUT,         &
     &                   nvd3,t2dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define wind-induced, bed wave orbital U-velocity.
!
        IF (Hout(idUbed,ng)) THEN
          Vinfo( 1)=Vname(1,idUbed)
          Vinfo( 2)=Vname(2,idUbed)
          Vinfo( 3)=Vname(3,idUbed)
          Vinfo(14)=Vname(4,idUbed)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
# endif
          status=def_var(nchisid(ng),hisVid(idUbed,ng),NF_FOUT,         &
     &                   nvd3,u2dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define wind-induced, bed wave orbital V-velocity.
!
        IF (Hout(idVbed,ng)) THEN
          Vinfo( 1)=Vname(1,idVbed)
          Vinfo( 2)=Vname(2,idVbed)
          Vinfo( 3)=Vname(3,idVbed)
          Vinfo(14)=Vname(4,idVbed)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
# endif
          status=def_var(nchisid(ng),hisVid(idVbed,ng),NF_FOUT,         &
     &                   nvd3,v2dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define bottom U-momentum above bed.
!
        IF (Hout(idUbot,ng)) THEN
          Vinfo( 1)=Vname(1,idUbot)
          Vinfo( 2)=Vname(2,idUbot)
          Vinfo( 3)=Vname(3,idUbot)
          Vinfo(14)=Vname(4,idUbot)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
# endif
          status=def_var(nchisid(ng),hisVid(idUbot,ng),NF_FOUT,         &
     &                   nvd3,u2dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define bottom V-momentum above bed.
!
        IF (Hout(idVbot,ng)) THEN
          Vinfo( 1)=Vname(1,idVbot)
          Vinfo( 2)=Vname(2,idVbot)
          Vinfo( 3)=Vname(3,idVbot)
          Vinfo(14)=Vname(4,idVbot)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
# endif
          status=def_var(nchisid(ng),hisVid(idVbot,ng),NF_FOUT,         &
     &                   nvd3,v2dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define bed ripple height.
!
        IF (Hout(idHrip,ng)) THEN
          Vinfo( 1)=Vname(1,idHrip)
          Vinfo( 2)=Vname(2,idHrip)
          Vinfo( 3)=Vname(3,idHrip)
          Vinfo(14)=Vname(4,idHrip)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          status=def_var(nchisid(ng),hisVid(idHrip,ng),NF_FOUT,         &
     &                   nvd3,t2dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define bed ripple length.
!
        IF (Hout(idLrip,ng)) THEN
          Vinfo( 1)=Vname(1,idLrip)
          Vinfo( 2)=Vname(2,idLrip)
          Vinfo( 3)=Vname(3,idLrip)
          Vinfo(14)=Vname(4,idLrip)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          status=def_var(nchisid(ng),hisVid(idLrip,ng),NF_FOUT,         &
     &                   nvd3,t2dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define hydraulic bottom roughness.
!
        IF (Hout(idZnot,ng)) THEN
          Vinfo( 1)=Vname(1,idZnot)
          Vinfo( 2)=Vname(2,idZnot)
          Vinfo( 3)=Vname(3,idZnot)
          Vinfo(14)=Vname(4,idZnot)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          status=def_var(nchisid(ng),hisVid(idZnot,ng),NF_FOUT,         &
     &                   nvd3,t2dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define apparent hydraulic bottom roughness.
!
        IF (Hout(idZapp,ng)) THEN
          Vinfo( 1)=Vname(1,idZapp)
          Vinfo( 2)=Vname(2,idZapp)
          Vinfo( 3)=Vname(3,idZapp)
          Vinfo(14)=Vname(4,idZapp)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          status=def_var(nchisid(ng),hisVid(idZapp,ng),NF_FOUT,         &
     &                   nvd3,t2dgrd,Aval,Vinfo,ncname)
        END IF
#endif
#ifdef SEDIMENT
!
!  Define bed layer sediment fraction for each size class.
!
        DO i=1,NST
          IF (Hout(idfrac(i),ng)) THEN
            Vinfo( 1)=Vname(1,idfrac(i))
            Vinfo( 2)=Vname(2,idfrac(i))
            Vinfo( 3)=Vname(3,idfrac(i))
            Vinfo(14)=Vname(4,idfrac(i))
            Vinfo(16)=Vname(1,idtime)
            WRITE (Vinfo(19),40) 1000.0_r8*Sd50(i,ng)
# if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
# endif
            status=def_var(nchisid(ng),hisVid(idfrac(i),ng),NF_FOUT,    &
     &                     nvd4,b3dgrd,Aval,Vinfo,ncname)
          END IF
        END DO
!
!  Define bed layer sediment properties.
!
        DO i=1,MBEDP
          IF (Hout(idSbed(i),ng)) THEN
            Vinfo( 1)=Vname(1,idSbed(i))
            Vinfo( 2)=Vname(2,idSbed(i))
            Vinfo( 3)=Vname(3,idSbed(i))
            Vinfo(14)=Vname(4,idSbed(i))
            Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
# endif
            status=def_var(nchisid(ng),hisVid(idSbed(i),ng),NF_FOUT,    &
     &                     nvd4,b3dgrd,Aval,Vinfo,ncname)
          END IF
        END DO
!
!  Define bottom sediment properties.
!
        DO i=1,MBOTP
          IF (Hout(idBott(i),ng)) THEN
            Vinfo( 1)=Vname(1,idBott(i))
            Vinfo( 2)=Vname(2,idBott(i))
            Vinfo( 3)=Vname(3,idBott(i))
            Vinfo(14)=Vname(4,idBott(i))
            Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
            status=def_var(nchisid(ng),hisVid(idBott(i),ng),NF_FOUT,    &
     &                     nvd3,t2dgrd,Aval,Vinfo,ncname)
          END IF
        END DO
#endif
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
        status=nf_enddef(nchisid(ng))
      END IF
!
!-----------------------------------------------------------------------
!  Write out time-recordless, information variables.
!-----------------------------------------------------------------------
!
      IF (ldef) THEN
        status=wrt_info(ng,nchisid(ng),OutThread,ncname)
        IF (exit_flag.ne.0) RETURN
      END IF
!
!=======================================================================
!  Open an existing history file, check its contents, and prepare for
!  appending data.
!=======================================================================
!
      IF (.not.ldef.and.OutThread) THEN
!
!  Inquire about the contents of history NetCDF file:  Inquire about
!  the dimensions and variables.  Check for consistency.
!
        ifield=0
        CALL opencdf (ng,1,ncname,fname,N(ng),ifield,nrec)
        IF (exit_flag.ne.0) RETURN
!
!  Open history file for read/write.
!
        status=nf_open(TRIM(ncname),nf_write,nchisid(ng))
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,50) TRIM(ncname)
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
!  history variables. Get variable IDs.
!
        DO i=1,nvars
          IF (varnam(i).eq.Vname(1,idtime)) THEN
            got_var(idtime)=.true.
            status=nf_inq_varid(nchisid(ng),TRIM(Vname(1,idtime)),      &
     &                          hisVid(idtime,ng))
          END IF
          IF (varnam(i).eq.Vname(1,idFsur)) THEN
            got_var(idFsur)=.true.
            status=nf_inq_varid(nchisid(ng),TRIM(Vname(1,idFsur)),      &
     &                          hisVid(idFsur,ng))
          END IF
          IF (varnam(i).eq.Vname(1,idUbar)) THEN
            got_var(idUbar)=.true.
            status=nf_inq_varid(nchisid(ng),TRIM(Vname(1,idUbar)),      &
     &                          hisVid(idUbar,ng))
          END IF
          IF (varnam(i).eq.Vname(1,idVbar)) THEN
            got_var(idVbar)=.true.
            status=nf_inq_varid(nchisid(ng),TRIM(Vname(1,idVbar)),      &
     &                          hisVid(idVbar,ng))
          END IF
#ifdef SOLVE3D
          IF (varnam(i).eq.Vname(1,idUvel)) THEN
            got_var(idUvel)=.true.
            status=nf_inq_varid(nchisid(ng),TRIM(Vname(1,idUvel)),      &
     &                          hisVid(idUvel,ng))
          END IF
          IF (varnam(i).eq.Vname(1,idVvel)) THEN
            got_var(idVvel)=.true.
            status=nf_inq_varid(nchisid(ng),TRIM(Vname(1,idVvel)),      &
     &                          hisVid(idVvel,ng))
          END IF
          IF (varnam(i).eq.Vname(1,idWvel)) THEN
            got_var(idWvel)=.true.
            status=nf_inq_varid(nchisid(ng),TRIM(Vname(1,idWvel)),      &
     &                          hisVid(idWvel,ng))
          END IF
          IF (varnam(i).eq.Vname(1,idOvel)) THEN
            got_var(idOvel)=.true.
            status=nf_inq_varid(nchisid(ng),TRIM(Vname(1,idOvel)),      &
     &                          hisVid(idOvel,ng))
          END IF
          DO itrc=1,NT(ng)
            IF (varnam(i).eq.Vname(1,idTvar(itrc))) THEN
              got_trc(itrc)=.true.
              status=nf_inq_varid(nchisid(ng),                          &
     &                            TRIM(Vname(1,idTvar(itrc))),          &
     &                            hisTid(itrc,ng))
            END IF
          END DO
          IF (varnam(i).eq.Vname(1,idDano)) THEN
            got_var(idDano)=.true.
            status=nf_inq_varid(nchisid(ng),TRIM(Vname(1,idDano)),      &
     &                          hisVid(idDano,ng))
          END IF
# ifdef LMD_SKPP
          IF (varnam(i).eq.Vname(1,idHsbl)) THEN
            got_var(idHsbl)=.true.
            status=nf_inq_varid(nchisid(ng),TRIM(Vname(1,idHsbl)),      &
     &                          hisVid(idHsbl,ng))
          END IF
# endif
# ifdef LMD_BKPP
          IF (varnam(i).eq.Vname(1,idHbbl)) THEN
            got_var(idHbbl)=.true.
            status=nf_inq_varid(nchisid(ng),TRIM(Vname(1,idHbbl)),      &
     &                          hisVid(idHbbl,ng))
          END IF
# endif
          IF (varnam(i).eq.Vname(1,idVvis)) THEN
            got_var(idVvis)=.true.
            status=nf_inq_varid(nchisid(ng),TRIM(Vname(1,idVvis)),      &
     &                          hisVid(idVvis,ng))
          END IF
          IF (varnam(i).eq.Vname(1,idTdif)) THEN
            got_var(idTdif)=.true.
            status=nf_inq_varid(nchisid(ng),TRIM(Vname(1,idTdif)),      &
     &                          hisVid(idTdif,ng))
          END IF
          IF (varnam(i).eq.Vname(1,idSdif)) THEN
            got_var(idSdif)=.true.
            status=nf_inq_varid(nchisid(ng),TRIM(Vname(1,idSdif)),      &
     &                          hisVid(idSdif,ng))
          END IF
# if defined GLS_MIXING || defined MY25_MIXING
          IF (varnam(i).eq.Vname(1,idMtke)) THEN
            got_var(idMtke)=.true.
            status=nf_inq_varid(nchisid(ng),TRIM(Vname(1,idMtke)),      &
     &                          hisVid(idMtke,ng))
          END IF
          IF (varnam(i).eq.Vname(1,idMtls)) THEN
            got_var(idMtls)=.true.
            status=nf_inq_varid(nchisid(ng),TRIM(Vname(1,idMtls)),      &
     &                          hisVid(idMtls,ng))
          END IF
# endif
          IF (varnam(i).eq.Vname(1,idTsur(itemp))) THEN
            got_var(idTsur(itemp))=.true.
            status=nf_inq_varid(nchisid(ng),                            &
     &                          TRIM(Vname(1,idTsur(itemp))),           &
     &                          hisVid(idTsur(itemp),ng))
          END IF
# ifdef BULK_FLUXES
          IF (varnam(i).eq.Vname(1,idLhea)) THEN
            got_var(idLhea)=.true.
            status=nf_inq_varid(nchisid(ng),TRIM(Vname(1,idLhea)),      &
     &                          hisVid(idLhea,ng))
          END IF
          IF (varnam(i).eq.Vname(1,idShea)) THEN
            got_var(idShea)=.true.
            status=nf_inq_varid(nchisid(ng),TRIM(Vname(1,idShea)),      &
     &                          hisVid(idShea,ng))
          END IF
          IF (varnam(i).eq.Vname(1,idLrad)) THEN
            got_var(idLrad)=.true.
            status=nf_inq_varid(nchisid(ng),TRIM(Vname(1,idLrad)),      &
     &                          hisVid(idLrad,ng))
          END IF
# endif
# ifdef SHORTWAVE
          IF (varnam(i).eq.Vname(1,idSrad)) THEN
            got_var(idSrad)=.true.
            status=nf_inq_varid(nchisid(ng),TRIM(Vname(1,idSrad)),      &
     &                          hisVid(idSrad,ng))
          END IF
# endif
# ifdef SEDIMENT
          DO itrc=1,NST
            IF (varnam(i).eq.Vname(1,idfrac(itrc))) THEN
              got_var(idfrac(itrc))=.true.
              status=nf_inq_varid(nchisid(ng),                          &
     &                            TRIM(Vname(1,idfrac(itrc))),          &
     &                            hisVid(idfrac(itrc),ng))
            END IF
          END DO
          DO itrc=1,MBEDP
            IF (varnam(i).eq.Vname(1,idSbed(itrc))) THEN
              got_var(idSbed(itrc))=.true.
              status=nf_inq_varid(nchisid(ng),                          &
     &                            TRIM(Vname(1,idSbed(itrc))),          &
     &                            hisVid(idSbed(itrc),ng))
            END IF
          END DO
          DO itrc=1,MBOTP
            IF (varnam(i).eq.Vname(1,idBott(itrc))) THEN
              got_var(idBott(itrc))=.true.
              status=nf_inq_varid(nchisid(ng),                          &
     &                            TRIM(Vname(1,idBott(itrc))),          &
     &                            hisVid(idBott(itrc),ng))
            END IF
          END DO
# endif
#endif
          IF (varnam(i).eq.Vname(1,idUsms)) THEN
            got_var(idUsms)=.true.
            status=nf_inq_varid(nchisid(ng),TRIM(Vname(1,idUsms)),      &
     &                          hisVid(idUsms,ng))
          END IF
          IF (varnam(i).eq.Vname(1,idVsms)) THEN
            got_var(idVsms)=.true.
            status=nf_inq_varid(nchisid(ng),TRIM(Vname(1,idVsms)),      &
     &                          hisVid(idVsms,ng))
          END IF
          IF (varnam(i).eq.Vname(1,idUbms)) THEN
            got_var(idUbms)=.true.
            status=nf_inq_varid(nchisid(ng),TRIM(Vname(1,idUbms)),      &
     &                          hisVid(idUbms,ng))
          END IF
          IF (varnam(i).eq.Vname(1,idVbms)) THEN
            got_var(idVbms)=.true.
            status=nf_inq_varid(nchisid(ng),TRIM(Vname(1,idVbms)),      &
     &                          hisVid(idVbms,ng))
          END IF
#ifdef BBL_MODEL
          IF (varnam(i).eq.Vname(1,idUbws)) THEN
            got_var(idUbws)=.true.
            status=nf_inq_varid(nchisid(ng),TRIM(Vname(1,idUbws)),      &
     &                          hisVid(idUbws,ng))
          END IF
          IF (varnam(i).eq.Vname(1,idVbws)) THEN
            got_var(idVbws)=.true.
            status=nf_inq_varid(nchisid(ng),TRIM(Vname(1,idVbws)),      &
     &                          hisVid(idVbws,ng))
          END IF
          IF (varnam(i).eq.Vname(1,idAbed)) THEN
            got_var(idAbed)=.true.
            status=nf_inq_varid(nchisid(ng),TRIM(Vname(1,idAbed)),      &
     &                          hisVid(idAbed,ng))
          END IF
          IF (varnam(i).eq.Vname(1,idUbed)) THEN
            got_var(idUbed)=.true.
            status=nf_inq_varid(nchisid(ng),TRIM(Vname(1,idUbed)),      &
     &                          hisVid(idUbed,ng))
          END IF
          IF (varnam(i).eq.Vname(1,idVbed)) THEN
            got_var(idVbed)=.true.
            status=nf_inq_varid(nchisid(ng),TRIM(Vname(1,idVbed)),      &
     &                          hisVid(idVbed,ng))
          END IF
          IF (varnam(i).eq.Vname(1,idUbot)) THEN
            got_var(idUbot)=.true.
            status=nf_inq_varid(nchisid(ng),TRIM(Vname(1,idUbot)),      &
     &                          hisVid(idUbot,ng))
          END IF
          IF (varnam(i).eq.Vname(1,idVbot)) THEN
            got_var(idVbot)=.true.
            status=nf_inq_varid(nchisid(ng),TRIM(Vname(1,idVbot)),      &
     &                          hisVid(idVbot,ng))
          END IF
          IF (varnam(i).eq.Vname(1,idHrip)) THEN
            got_var(idHrip)=.true.
            status=nf_inq_varid(nchisid(ng),TRIM(Vname(1,idHrip)),      &
     &                          hisVid(idHrip,ng))
          END IF
          IF (varnam(i).eq.Vname(1,idLrip)) THEN
            got_var(idLrip)=.true.
            status=nf_inq_varid(nchisid(ng),TRIM(Vname(1,idLrip)),      &
     &                          hisVid(idLrip,ng))
          END IF
          IF (varnam(i).eq.Vname(1,idZnot)) THEN
            got_var(idZnot)=.true.
            status=nf_inq_varid(nchisid(ng),TRIM(Vname(1,idZnot)),      &
     &                          hisVid(idZnot,ng))
          END IF
          IF (varnam(i).eq.Vname(1,idZapp)) THEN
            got_var(idZapp)=.true.
            status=nf_inq_varid(nchisid(ng),TRIM(Vname(1,idZapp)),      &
     &                          hisVid(idZapp,ng))
          END IF
#endif
        END DO
!
!  Check if history variables are available in input NetCDF file.
!
        IF (.not.got_var(idtime)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idtime)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idFsur).and.Hout(idFsur,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idFsur)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbar).and.Hout(idUbar,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idUbar)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbar).and.Hout(idVbar,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idVbar)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#ifdef SOLVE3D
        IF (.not.got_var(idUvel).and.Hout(idUvel,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idUvel)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVvel).and.Hout(idVvel,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idVvel)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWvel).and.Hout(idWvel,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idWvel)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idOvel).and.Hout(idOvel,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idOvel)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        DO itrc=1,NT(ng)
          IF (.not.got_trc(itrc).and.Hout(idTvar(itrc),ng)) THEN
            WRITE (stdout,60) TRIM(Vname(1,idTvar(itrc))), TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
        IF (.not.got_var(idDano).and.Hout(idDano,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idDano)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# ifdef LMD_SKPP
        IF (.not.got_var(idHsbl).and.Hout(idHsbl,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idHsbl)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# ifdef LMD_BKPP
        IF (.not.got_var(idHbbl).and.Hout(idHbbl,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idHbbl)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
        IF (.not.got_var(idVvis).and.Hout(idVvis,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idVvis)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTdif).and.Hout(idTdif,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idTdif)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# ifdef SALINITY
        IF (.not.got_var(idSdif).and.Hout(idSdif,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idSdif)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# if defined GLS_MIXING || defined MY25_MIXING
        IF (.not.got_var(idMtke).and.Hout(idMtke,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idMtke)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idMtls).and.Hout(idMtls,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idMtls)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
        IF (.not.got_var(idTsur(itemp)).and.Hout(idTsur(itemp),ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idTsur(itemp))), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# ifdef BULK_FLUXES
        IF (.not.got_var(idLhea).and.Hout(idLhea,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idLhea)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idShea).and.Hout(idShea,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idShea)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idLrad).and.Hout(idLrad,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idLrad)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# ifdef SHORTWAVE
        IF (.not.got_var(idSrad).and.Hout(idSrad,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idSrad)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# ifdef SEDIMENT
        DO itrc=1,NST
          IF (.not.got_var(idfrac(itrc)).and.Hout(idfrac(itrc),ng)) THEN
            WRITE (stdout,60) TRIM(Vname(1,idfrac(itrc))), TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
        DO itrc=1,MBEDP
          IF (.not.got_var(idSbed(itrc)).and.Hout(idSbed(itrc),ng)) THEN
            WRITE (stdout,60) TRIM(Vname(1,idSbed(itrc))), TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
        DO itrc=1,MBOTP
          IF (.not.got_var(idBott(itrc)).and.Hout(idBott(itrc),ng)) THEN
            WRITE (stdout,60) TRIM(Vname(1,idBott(itrc))), TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
# endif
#endif
        IF (.not.got_var(idUsms).and.Hout(idUsms,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idUsms)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVsms).and.Hout(idVsms,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idVsms)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbms).and.Hout(idUbms,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idUbms)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbms).and.Hout(idVbms,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idVbms)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#ifdef BBL_MODEL
        IF (.not.got_var(idUbws).and.Hout(idUbws,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idUbws)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbws).and.Hout(idVbws,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idVbws)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idAbed).and.Hout(idAbed,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idAbed)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbed).and.Hout(idUbed,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idUbed)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbed).and.Hout(idVbed,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idVbed)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbot).and.Hout(idUbot,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idUbot)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbot).and.Hout(idVbot,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idVbot)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idHrip).and.Hout(idHrip,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idHrip)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idLrip).and.Hout(idLrip,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idLrip)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idZnot).and.Hout(idZnot,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idZnot)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idZapp).and.Hout(idZapp,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idZapp)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#endif
!
!  Set unlimited time record dimension to the appropriate value.
!
!!      thisindx(ng)=tsize
        IF (ndefhis(ng).gt.0) THEN
          IF ((ntstart-1).lt.ndefhis(ng)) THEN
            thisindx(ng)=1+((ntstart-1)-                                &
     &                       ndefhis(ng)*((ntstart-1)/ndefhis(ng)))/    &
     &                     nhis(ng)
          ELSE
            thisindx(ng)=((ntstart-1)-                                  &
     &                     ndefhis(ng)*((ntstart-1)/ndefhis(ng)))/      &
     &                   nhis(ng)
          END IF
        ELSE
          thisindx(ng)=(ntstart-1)/nhis(ng)
        END IF
        thisindx(ng)=MIN(thisindx(ng),tsize)
      END IF
!
  10  FORMAT (6x,'DEF_HIS - creating history file: ',a)
  20  FORMAT (6x,'DEF_HIS - inquiring history file: ',a)
  30  FORMAT (/,' DEF_HIS - unable to create history NetCDF file: ',a)
  40  FORMAT (1pe11.4,1x,'millimeter')
  50  FORMAT (/,' DEF_HIS - unable to open history NetCDF file: ',a)
  60  FORMAT (/,' DEF_HIS - unable to find variable: ',a,2x,            &
     &         ' in history NetCDF file: ',a)
      RETURN
      END SUBROUTINE def_his
