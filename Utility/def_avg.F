#include "cppdefs.h"
#ifdef AVERAGES
      SUBROUTINE def_avg (ng,ldef)
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine creates averages NetCDF file, it defines its           !
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
# ifdef SEDIMENT
      USE mod_sediment
# endif
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
      integer, parameter :: Natt = 20

      logical :: got_var(NV)

      integer :: i, itrc, j, nrec, nvd3, nvd4
      integer :: recdim, status

      integer :: DimIDs(28), t2dgrd(3), u2dgrd(3), v2dgrd(3)

      integer :: def_info, def_var, wrt_info

# ifdef SOLVE3D
#  ifdef SEDIMENT
      integer :: b3dgrd(4)
#  endif
      integer :: t3dgrd(4), u3dgrd(4), v3dgrd(4), w3dgrd(4)
# endif

      real(r8) :: Aval(4)

      character (len=80) :: Vinfo(Natt)
      character (len=80) :: fname, ncname
!
!-----------------------------------------------------------------------
!  Set and report file name.
!-----------------------------------------------------------------------
!
      IF (exit_flag.ne.0) RETURN
      ncname=avgname(ng)
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
!  Create a new averages NetCDF file.
!=======================================================================
!
      IF (ldef.and.OutThread) THEN
        status=nf_create(TRIM(ncname),nf_clobber,ncavgid(ng))
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,30) TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
      END IF
#ifdef DISTRIBUTE
      CALL mp_bcasti (ng, ncavgid(ng), 1)
#endif
!
!-----------------------------------------------------------------------
!  Define the dimensions of staggered fields.
!-----------------------------------------------------------------------
!
      IF (ldef.and.OutThread) THEN
        status=nf_def_dim(ncavgid(ng),'xi_rho'  ,Lm(ng)+2,              &
     &                    DimIDs( 1))
        status=nf_def_dim(ncavgid(ng),'xi_u'    ,Lm(ng)+1,              &
     &                    DimIDs( 2))
        status=nf_def_dim(ncavgid(ng),'xi_v'    ,Lm(ng)+2,              &
     &                    DimIDs( 3))
        status=nf_def_dim(ncavgid(ng),'xi_psi'  ,Lm(ng)+1,              &
     &                    DimIDs( 4))
        status=nf_def_dim(ncavgid(ng),'eta_rho' ,Mm(ng)+2,              &
     &                    DimIDs( 5))
        status=nf_def_dim(ncavgid(ng),'eta_u'   ,Mm(ng)+2,              &
     &                    DimIDs( 6))
        status=nf_def_dim(ncavgid(ng),'eta_v'   ,Mm(ng)+1,              &
     &                    DimIDs( 7))
        status=nf_def_dim(ncavgid(ng),'eta_psi' ,Mm(ng)+1,              &
     &                    DimIDs( 8))
# if defined WRITE_WATER && defined MASKING
        status=nf_def_dim(ncavgid(ng),'xy_rho'  ,Nxyr(ng),              &
     &                    DimIDs(17))
        status=nf_def_dim(ncavgid(ng),'xy_u'    ,Nxyu(ng),              &
     &                    DimIDs(18))
        status=nf_def_dim(ncavgid(ng),'xy_v'    ,Nxyv(ng),              &
     &                    DimIDs(19))
# endif
# ifdef SOLVE3D
#  if defined WRITE_WATER && defined MASKING
        status=nf_def_dim(ncavgid(ng),'xyz_rho' ,Nxyr(ng)*N(ng),        &
     &                    DimIDs(20))
        status=nf_def_dim(ncavgid(ng),'xyz_u'   ,Nxyu(ng)*N(ng),        &
     &                    DimIDs(21))
        status=nf_def_dim(ncavgid(ng),'xyz_v'   ,Nxyv(ng)*N(ng),        &
     &                    DimIDs(22))
        status=nf_def_dim(ncavgid(ng),'xyz_w'   ,Nxyr(ng)*(N(ng)+1),    &
     &                    DimIDs(23))
#  endif
        status=nf_def_dim(ncavgid(ng),'s_rho'   ,N(ng),                 &
     &                    DimIDs( 9))
        status=nf_def_dim(ncavgid(ng),'s_w'     ,N(ng)+1,               &
     &                    DimIDs(10))
        status=nf_def_dim(ncavgid(ng),'tracer'  ,NT(ng),                &
     &                    DimIDs(11))
#  ifdef SEDIMENT
        status=nf_def_dim(ncavgid(ng),'Nbed'    ,Nbed,                  &
     &                    DimIDs(16))
#   if defined WRITE_WATER && defined MASKING
        status=nf_def_dim(ncavgid(ng),'xybed'   ,Nxyr(ng)*Nbed,         &
     &                    DimIDs(24))
#   endif
#  endif
#  ifdef ECOSIM
        status=nf_def_dim(ncavgid(ng),'Nphy'    ,Nphy,                  &
     &                    DimIDs(25))
        status=nf_def_dim(ncavgid(ng),'Nbac'    ,Nbac,                  &
     &                    DimIDs(26))
        status=nf_def_dim(ncavgid(ng),'Ndom'    ,Ndom,                  &
     &                    DimIDs(27))
        status=nf_def_dim(ncavgid(ng),'Nfec'    ,Nfec,                  &
     &                    DimIDs(28))
#  endif
# endif
        status=nf_def_dim(ncavgid(ng),'boundary',4,                     &
     &                    DimIDs(14))
        status=nf_def_dim(ncavgid(ng),TRIM(ADJUSTL(Vname(5,idtime))),   &
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
# ifdef SOLVE3D
!
!  Define dimension vector for staggered w-momemtum type variables.
!
#  if defined WRITE_WATER && defined MASKING
        w3dgrd(1)=DimIDs(23)
        w3dgrd(2)=DimIDs(12)
#  else
        w3dgrd(1)=DimIDs( 1)
        w3dgrd(2)=DimIDs( 5)
        w3dgrd(3)=DimIDs(10)
        w3dgrd(4)=DimIDs(12)
#  endif
#  ifdef SEDIMENT
!
!  Define dimension vector for sediment bed layer type variables.
!
#   if defined WRITE_WATER && defined MASKING
        b3dgrd(1)=DimIDs(24)
        b3dgrd(2)=DimIDs(12)
#   else
        b3dgrd(1)=DimIDs( 1)
        b3dgrd(2)=DimIDs( 5)
        b3dgrd(3)=DimIDs(16)
        b3dgrd(4)=DimIDs(12)
#   endif
#  endif
# endif
!
!  Initialize unlimited time record dimension.
!
        tavgindx(ng)=0
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
        status=def_info(ng,ncavgid(ng),ncname,DimIDs)
        IF (exit_flag.ne.0) RETURN
!
!-----------------------------------------------------------------------
!  Define variables and their attributes.
!-----------------------------------------------------------------------
!
!  Define model time.
!
        Vinfo( 1)=Vname(1,idtime)
        Vinfo( 2)='averaged '//TRIM(Vname(2,idtime))
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
        status=def_var(ncavgid(ng),avgVid(idtime,ng),NF_TYPE,           &
     &                 1,recdim,Aval,Vinfo,ncname)
!
!  Define free-surface.
!
        IF (Hout(idFsur,ng)) THEN
          Vinfo( 1)=Vname(1,idFsur)
          Vinfo( 2)='averaged '//TRIM(Vname(2,idFsur))
          Vinfo( 3)=Vname(3,idFsur)
          Vinfo(14)=Vname(4,idFsur)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          status=def_var(ncavgid(ng),avgVid(idFsur,ng),NF_FOUT,         &
     &                   nvd3,t2dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define 2D momentum in the XI-direction.
!
        IF (Hout(idUbar,ng)) THEN
          Vinfo( 1)=Vname(1,idUbar)
          Vinfo( 2)='averaged '//TRIM(Vname(2,idUbar))
          Vinfo( 3)=Vname(3,idUbar)
          Vinfo(14)=Vname(4,idUbar)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
# endif
          status=def_var(ncavgid(ng),avgVid(idUbar,ng),NF_FOUT,         &
     &                   nvd3,u2dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define 2D momentum in the ETA-direction.
!
        IF (Hout(idVbar,ng)) THEN
          Vinfo( 1)=Vname(1,idVbar)
          Vinfo( 2)='averaged '//TRIM(Vname(2,idVbar))
          Vinfo( 3)=Vname(3,idVbar)
          Vinfo(14)=Vname(4,idVbar)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
# endif
          status=def_var(ncavgid(ng),avgVid(idVbar,ng),NF_FOUT,         &
     &                   nvd3,v2dgrd,Aval,Vinfo,ncname)
        END IF
# ifdef SOLVE3D
!
!  Define 3D momentum component in the XI-direction.
!
        IF (Hout(idUvel,ng)) THEN
          Vinfo( 1)=Vname(1,idUvel)
          Vinfo( 2)='averaged '//TRIM(Vname(2,idUvel))
          Vinfo( 3)=Vname(3,idUvel)
          Vinfo(14)=Vname(4,idUvel)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#  endif
          status=def_var(ncavgid(ng),avgVid(idUvel,ng),NF_FOUT,         &
     &                   nvd4,u3dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define 3D momentum component in the ETA-direction.
!
        IF (Hout(idVvel,ng)) THEN
          Vinfo( 1)=Vname(1,idVvel)
          Vinfo( 2)='averaged '//TRIM(Vname(2,idVvel))
          Vinfo( 3)=Vname(3,idVvel)
          Vinfo(14)=Vname(4,idVvel)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#  endif
          status=def_var(ncavgid(ng),avgVid(idVvel,ng),NF_FOUT,         &
     &                   nvd4,v3dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define S-coordinate vertical "omega" momentum component.
!
        IF (Hout(idOvel,ng)) THEN
          Vinfo( 1)=Vname(1,idOvel)
          Vinfo( 2)='averaged '//TRIM(Vname(2,idOvel))
          Vinfo( 3)=Vname(3,idOvel)
          Vinfo(14)=Vname(4,idOvel)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          status=def_var(ncavgid(ng),avgVid(idOvel,ng),NF_FOUT,         &
     &                   nvd4,w3dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define tracer type variables.
!
        DO itrc=1,NT(ng)
          IF (Hout(idTvar(itrc),ng)) THEN
            Vinfo( 1)=Vname(1,idTvar(itrc))
            Vinfo( 2)='averaged '//TRIM(Vname(2,idTvar(itrc)))
            Vinfo( 3)=Vname(3,idTvar(itrc))
            Vinfo(14)=Vname(4,idTvar(itrc))
            Vinfo(16)=Vname(1,idtime)
#  ifdef SEDIMENT
            DO i=1,NST
              IF (itrc.eq.idsed(i)) THEN
                WRITE (Vinfo(19),40) 1000.0_r8*Sd50(i,ng)
              END IF
            END DO
#  endif
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
            status=def_var(ncavgid(ng),avgTid(itrc,ng),NF_FOUT,         &
     &                     nvd4,t3dgrd,Aval,Vinfo,ncname)
          END IF
        END DO
!
!  Define density anomaly.
!
        IF (Hout(idDano,ng)) THEN
          Vinfo( 1)=Vname(1,idDano)
          Vinfo( 2)='averaged '//TRIM(Vname(2,idDano))
          Vinfo( 3)=Vname(3,idDano)
          Vinfo(14)=Vname(4,idDano)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          status=def_var(ncavgid(ng),avgVid(idDano,ng),NF_FOUT,         &
     &                   nvd4,t3dgrd,Aval,Vinfo,ncname)
        END IF
#  ifdef LMD_SKPP
!
!  Define depth of surface boundary layer.
!
        IF (Hout(idHsbl,ng)) THEN
          Vinfo( 1)=Vname(1,idHsbl)
          Vinfo( 2)='averaged '//TRIM(Vname(2,idHsbl))
          Vinfo( 3)=Vname(3,idHsbl)
          Vinfo(14)=Vname(4,idHsbl)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          status=def_var(ncavgid(ng),avgVid(idHsbl,ng),NF_FOUT,         &
     &                   nvd3,t2dgrd,Aval,Vinfo,ncname)
        END IF
#  endif
#  ifdef LMD_BKPP
!
!  Define depth of bottom boundary layer.
!
        IF (Hout(idHbbl,ng)) THEN
          Vinfo( 1)=Vname(1,idHbbl)
          Vinfo( 2)='averaged '//TRIM(Vname(2,idHbbl))
          Vinfo( 3)=Vname(3,idHbbl)
          Vinfo(14)=Vname(4,idHbbl)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          status=def_var(ncavgid(ng),avgVid(idHbbl,ng),NF_FOUT,         &
     &                   nvd3,t2dgrd,Aval,Vinfo,ncname)
        END IF
#  endif
#  ifdef AVERAGES_QUADRATIC
!
!  Define quadratic <uu> term.
!
        Vinfo( 1)=Vname(1,idUUav)
        Vinfo( 2)=TRIM(Vname(2,idUUav))
        Vinfo( 3)=Vname(3,idUUav)
        Vinfo(14)=Vname(4,idUUav)
        Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_u'
#   endif
        status=def_var(ncavgid(ng),avgVid(idUUav,ng),NF_FOUT,           &
     &                 nvd4,u3dgrd,Aval,Vinfo,ncname)
!
!  Define quadratic <uv> term.
!
        Vinfo( 1)=Vname(1,idUVav)
        Vinfo( 2)=TRIM(Vname(2,idUVav))
        Vinfo( 3)=Vname(3,idUVav)
        Vinfo(14)=Vname(4,idUVav)
        Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
#   endif
        status=def_var(ncavgid(ng),avgVid(idUVav,ng),NF_FOUT,           &
     &                 nvd4,t3dgrd,Aval,Vinfo,ncname)
!
!  Define quadratic <vv> term.
!
        Vinfo( 1)=Vname(1,idVVav)
        Vinfo( 2)=TRIM(Vname(2,idVVav))
        Vinfo( 3)=Vname(3,idVVav)
        Vinfo(14)=Vname(4,idVVav)
        Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_v'
#   endif
        status=def_var(ncavgid(ng),avgVid(idVVav,ng),NF_FOUT,           &
     &                 nvd4,v3dgrd,Aval,Vinfo,ncname)
!
!  Define quadratic <ut> and <vt> terms.
!
        DO itrc=1,NAT
          Vinfo( 1)=Vname(1,idUTav(itrc))
          Vinfo( 2)=TRIM(Vname(2,idUTav(itrc)))
          Vinfo( 3)=Vname(3,idUTav(itrc))
          Vinfo(14)=Vname(4,idUTav(itrc))
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#  endif
          status=def_var(ncavgid(ng),avgVid(idUTav(itrc),ng),NF_FOUT,   &
     &                   nvd4,u3dgrd,Aval,Vinfo,ncname)
!
          Vinfo( 1)=Vname(1,idVTav(itrc))
          Vinfo( 2)=TRIM(Vname(2,idVTav(itrc)))
          Vinfo( 3)=Vname(3,idVTav(itrc))
          Vinfo(14)=Vname(4,idVTav(itrc))
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#  endif
          status=def_var(ncavgid(ng),avgVid(idVTav(itrc),ng),NF_FOUT,   &
     &                   nvd4,v3dgrd,Aval,Vinfo,ncname)
        END DO
#  endif
#  ifdef AVERAGES_AKV
!
!  Define vertical viscosity coefficient.
!
        Vinfo( 1)=Vname(1,idVvis)
        Vinfo( 2)='averaged '//TRIM(Vname(2,idVvis))
        Vinfo( 3)=Vname(3,idVvis)
        Vinfo(14)=Vname(4,idVvis)
        Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
#   endif
        status=def_var(ncavgid(ng),avgVid(idVvis,ng),NF_FOUT,           &
     &                 nvd4,w3dgrd,Aval,Vinfo,ncname)
#  endif
#  ifdef AVERAGES_AKT
!
!  Define vertical diffusion coefficient for potential temperature.
!
        Vinfo( 1)=Vname(1,idTdif)
        Vinfo( 2)='averaged '//TRIM(Vname(2,idTdif))
        Vinfo( 3)=Vname(3,idTdif)
        Vinfo(14)=Vname(4,idTdif)
        Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
#   endif
        status=def_var(ncavgid(ng),avgVid(idTdif,ng),NF_FOUT,           &
     &                 nvd4,w3dgrd,Aval,Vinfo,ncname)
#  endif
#  ifdef AVERAGES_AKS
!
!  Define vertical diffusion coefficient for salinity.
!
        Vinfo( 1)=Vname(1,idSdif)
        Vinfo( 2)='averaged '//TRIM(Vname(2,idSdif))
        Vinfo( 3)=Vname(3,idSdif)
        Vinfo(14)=Vname(4,idSdif)
        Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
        status=def_var(ncavgid(ng),avgVid(idSdif,ng),NF_FOUT,           &
     &                 nvd4,w3dgrd,Aval,Vinfo,ncname)
#  endif
#  ifdef AVERAGES_FLUXES
!
!  Define surface net heat flux.
!
        IF (Hout(idTsur(itemp),ng)) THEN
          Vinfo( 1)=Vname(1,idTsur(itemp))
          Vinfo( 2)='averaged '//TRIM(Vname(2,idTsur(itemp)))
          Vinfo( 3)=Vname(3,idTsur(itemp))
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idTsur(itemp))
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          status=def_var(ncavgid(ng),avgVid(idTsur(itemp),ng),NF_FOUT,  &
     &                   nvd3,t2dgrd,Aval,Vinfo,ncname)
        END IF
#   ifdef BULK_FLUXES
!
!  Define latent heat flux.
!
        IF (Hout(idLhea,ng)) THEN
          Vinfo( 1)=Vname(1,idLhea)
          Vinfo( 2)='averaged '//TRIM(Vname(2,idLhea))
          Vinfo( 3)=Vname(3,idLhea)
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idLhea)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#    endif
          status=def_var(ncavgid(ng),avgVid(idLhea,ng),NF_FOUT,         &
     &                   nvd3,t2dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define sensible heat flux.
!
        IF (Hout(idShea,ng)) THEN
          Vinfo( 1)=Vname(1,idShea)
          Vinfo( 2)='averaged '//TRIM(Vname(2,idShea))
          Vinfo( 3)=Vname(3,idShea)
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idShea)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#    endif
          status=def_var(ncavgid(ng),avgVid(idShea,ng),NF_FOUT,         &
     &                   nvd3,t2dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define longwave radiation flux.
!
        IF (Hout(idLrad,ng)) THEN
          Vinfo( 1)=Vname(1,idLrad)
          Vinfo( 2)='averaged '//TRIM(Vname(2,idLrad))
          Vinfo( 3)=Vname(3,idLrad)
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idLrad)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#    endif
          status=def_var(ncavgid(ng),avgVid(idLrad,ng),NF_FOUT,         &
     &                   nvd3,t2dgrd,Aval,Vinfo,ncname)
        END IF
#   endif
#   ifdef SHORTWAVE
!
!  Define shortwave radiation flux.
!
        IF (Hout(idSrad,ng)) THEN
          Vinfo( 1)=Vname(1,idSrad)
          Vinfo( 2)='averaged '//TRIM(Vname(2,idSrad))
          Vinfo( 3)=Vname(3,idSrad)
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idSrad)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#    endif
          status=def_var(ncavgid(ng),avgVid(idSrad,ng),NF_FOUT,         &
     &                   nvd3,t2dgrd,Aval,Vinfo,ncname)
        END IF
#   endif
#  endif
# endif
# ifdef AVERAGES_FLUXES
!
!  Define surface u-momentum stress.
!
        IF (Hout(idUsms,ng)) THEN
          Vinfo( 1)=Vname(1,idUsms)
          Vinfo( 2)='averaged '//TRIM(Vname(2,idUsms))
          Vinfo( 3)=Vname(3,idUsms)
          Vinfo(14)=Vname(4,idUsms)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#  endif
          status=def_var(ncavgid(ng),avgVid(idUsms,ng),NF_FOUT,         &
     &                   nvd3,u2dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define surface v-momentum stress.
!
        IF (Hout(idVsms,ng)) THEN
          Vinfo( 1)=Vname(1,idVsms)
          Vinfo( 2)='averaged '//TRIM(Vname(2,idVsms))
          Vinfo( 3)=Vname(3,idVsms)
          Vinfo(14)=Vname(4,idVsms)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#  endif
          status=def_var(ncavgid(ng),avgVid(idVsms,ng),NF_FOUT,         &
     &                   nvd3,v2dgrd,Aval,Vinfo,ncname)
        END IF
# endif
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
        status=nf_enddef(ncavgid(ng))
      END IF
!
!-----------------------------------------------------------------------
!  Write out time-recordless, information variables.
!-----------------------------------------------------------------------
!
      IF (ldef) THEN
        status=wrt_info(ng,ncavgid(ng),OutThread,ncname)
        IF (exit_flag.ne.0) RETURN
      END IF
!
!=======================================================================
!  Open an existing averages file, check its contents, and prepare
!  for appending data.
!=======================================================================
!
      IF (.not.ldef.and.OutThread) THEN
!
!  Inquire about the contents of averages NetCDF file:  Inquire about
!  the dimensions and variables.  Check for consistency.
!
        CALL opencdf (ng,1,ncname,fname,N(ng),0,nrec)
        IF (exit_flag.ne.0) RETURN
!
!  Open averages file for read/write.
!
        status=nf_open(TRIM(ncname), nf_write, ncavgid(ng))
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
!
!  Scan variable list from input NetCDF and activate switches for
!  average variables. Get variable IDs.
!
        DO i=1,nvars
          IF (varnam(i).eq.Vname(1,idtime)) THEN
            got_var(idtime)=.true.
            status=nf_inq_varid(ncavgid(ng),TRIM(Vname(1,idtime)),      &
     &                          avgVid(idtime,ng))
          END IF
          IF (varnam(i).eq.Vname(1,idFsur)) THEN
            got_var(idFsur)=.true.
            status=nf_inq_varid(ncavgid(ng),TRIM(Vname(1,idFsur)),      &
     &                          avgVid(idFsur,ng))
          END IF
          IF (varnam(i).eq.Vname(1,idUbar)) THEN
            got_var(idUbar)=.true.
            status=nf_inq_varid(ncavgid(ng),TRIM(Vname(1,idUbar)),      &
     &                          avgVid(idUbar,ng))
          END IF
          IF (varnam(i).eq.Vname(1,idVbar)) THEN
            got_var(idVbar)=.true.
            status=nf_inq_varid(ncavgid(ng),TRIM(Vname(1,idVbar)),      &
     &                          avgVid(idVbar,ng))
          END IF
# ifdef SOLVE3D
          IF (varnam(i).eq.Vname(1,idUvel)) THEN
            got_var(idUvel)=.true.
            status=nf_inq_varid(ncavgid(ng),TRIM(Vname(1,idUvel)),      &
     &                          avgVid(idUvel,ng))
          END IF
          IF (varnam(i).eq.Vname(1,idVvel)) THEN
            got_var(idVvel)=.true.
            status=nf_inq_varid(ncavgid(ng),TRIM(Vname(1,idVvel)),      &
     &                          avgVid(idVvel,ng))
          END IF
          IF (varnam(i).eq.Vname(1,idOvel)) THEN
            got_var(idOvel)=.true.
            status=nf_inq_varid(ncavgid(ng),TRIM(Vname(1,idOvel)),      &
     &                          avgVid(idOvel,ng))
          END IF
          DO itrc=1,NT(ng)
            IF (varnam(i).eq.Vname(1,idTvar(itrc))) THEN
             got_var(idTvar(itrc))=.true.
             status=nf_inq_varid(ncavgid(ng),                           &
     &                           TRIM(Vname(1,idTvar(itrc))),           &
     &                           avgTid(itrc,ng))
            END IF
          END DO
          IF (varnam(i).eq.Vname(1,idDano)) THEN
            got_var(idDano)=.true.
            status=nf_inq_varid(ncavgid(ng),TRIM(Vname(1,idDano)),      &
     &                          avgVid(idDano,ng))
          END IF
#  ifdef LMD_SKPP
          IF (varnam(i).eq.Vname(1,idHsbl)) THEN
            got_var(idHsbl)=.true.
            status=nf_inq_varid(ncavgid(ng),TRIM(Vname(1,idHsbl)),      &
     &                          avgVid(idHsbl,ng))
          END IF
#  endif
#  ifdef LMD_BKPP
          IF (varnam(i).eq.Vname(1,idHbbl)) THEN
            got_var(idHbbl)=.true.
            status=nf_inq_varid(ncavgid(ng),TRIM(Vname(1,idHbbl)),      &
     &                          avgVid(idHbbl,ng))
          END IF
#  endif
#  ifdef AVERAGES_QUADRATIC
          IF (varnam(i).eq.Vname(1,idUUav)) THEN
            got_var(idUUav)=.true.
            status=nf_inq_varid(ncavgid(ng),TRIM(Vname(1,idUUav)),      &
     &                          avgVid(idUUav,ng))
          END IF
          IF (varnam(i).eq.Vname(1,idUVav)) THEN
            got_var(idUVav)=.true.
            status=nf_inq_varid(ncavgid(ng),TRIM(Vname(1,idUVav)),      &
     &                          avgVid(idUVav,ng))
          END IF
          IF (varnam(i).eq.Vname(1,idVVav)) THEN
            got_var(idVVav)=.true.
            status=nf_inq_varid(ncavgid(ng),TRIM(Vname(1,idVVav)),      &
     &                          avgVid(idVVav,ng))
          END IF
          DO itrc=1,NAT
            IF (varnam(i).eq.Vname(1,idUTav(itrc))) THEN
             got_var(idUTav(itrc))=.true.
             status=nf_inq_varid(ncavgid(ng),                           &
     &                           TRIM(Vname(1,idUTav(itrc))),           &
     &                           avgVid(idUTav(itrc),ng))
            END IF
            IF (varnam(i).eq.Vname(1,idVTav(itrc))) THEN
             got_var(idVTav(itrc))=.true.
             status=nf_inq_varid(ncavgid(ng),                           &
     &                           TRIM(Vname(1,idVTav(itrc))),           &
     &                           avgVid(idVTav(itrc),ng))
            END IF
          END DO
#  endif
#  ifdef AVERAGES_AKV
          IF (varnam(i).eq.Vname(1,idVvis)) THEN
            got_var(idVvis)=.true.
            status=nf_inq_varid(ncavgid(ng),TRIM(Vname(1,idVvis)),      &
     &                          avgVid(idVvis,ng))
          END IF
#  endif
#  ifdef AVERAGES_AKT
          IF (varnam(i).eq.Vname(1,idTdif)) THEN
            got_var(idTdif)=.true.
           status=nf_inq_varid(ncavgid(ng),TRIM(Vname(1,idTdif)),       &
     &                         avgVid(idTdif,ng))
          END IF
#  endif
#  ifdef AVERAGES_AKS
          IF (varnam(i).eq.Vname(1,idSdif)) THEN
            got_var(idSdif)=.true.
            status=nf_inq_varid(ncavgid(ng),TRIM(Vname(1,idSdif)),      &
     &                          avgVid(idSdif,ng))
          END IF
#  endif
#  ifdef AVERAGES_FLUXES
          IF (varnam(i).eq.Vname(1,idTsur(itemp))) THEN
            got_var(idTsur(itemp))=.true.
            status=nf_inq_varid(ncavgid(ng),                            &
     &                          TRIM(Vname(1,idTsur(itemp))),           &
     &                          avgVid(idTsur(itemp),ng))
          END IF
#   ifdef BULK_FLUXES
          IF (varnam(i).eq.Vname(1,idLhea)) THEN
            got_var(idLhea)=.true.
            status=nf_inq_varid(ncavgid(ng),TRIM(Vname(1,idLhea)),      &
     &                          avgVid(idLhea,ng))
          END IF
          IF (varnam(i).eq.Vname(1,idShea)) THEN
            got_var(idShea)=.true.
            status=nf_inq_varid(ncavgid(ng),TRIM(Vname(1,idShea)),      &
     &                          avgVid(idShea,ng))
          END IF
          IF (varnam(i).eq.Vname(1,idLrad)) THEN
            got_var(idLrad)=.true.
            status=nf_inq_varid(ncavgid(ng),TRIM(Vname(1,idLrad)),      &
     &                          avgVid(idLrad,ng))
          END IF
#   endif
#   ifdef SHORTWAVE
          IF (varnam(i).eq.Vname(1,idSrad)) THEN
            got_var(idSrad)=.true.
            status=nf_inq_varid(ncavgid(ng),TRIM(Vname(1,idSrad)),      &
     &                          avgVid(idSrad,ng))
          END IF
#   endif
#  endif
# endif
# ifdef AVERAGES_FLUXES
          IF (varnam(i).eq.Vname(1,idUsms)) THEN
            got_var(idUsms)=.true.
            status=nf_inq_varid(ncavgid(ng),TRIM(Vname(1,idUsms)),      &
     &                          avgVid(idUsms,ng))
          END IF
          IF (varnam(i).eq.Vname(1,idVsms)) THEN
            got_var(idVsms)=.true.
            status=nf_inq_varid(ncavgid(ng),TRIM(Vname(1,idVsms)),      &
     &                          avgVid(idVsms,ng))
          END IF
# endif
        END DO
!
!  Check if averages variables are available in input NetCDF file.
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
# ifdef SOLVE3D
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
        IF (.not.got_var(idOvel).and.Hout(idOvel,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idOvel)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        DO itrc=1,NT(ng)
          IF (.not.got_var(idTvar(itrc)).and.Hout(idTvar(itrc),ng)) THEN
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
#  ifdef LMD_SKPP
        IF (.not.got_var(idHsbl).and.Hout(idHsbl,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idHsbl)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  ifdef LMD_BKPP
        IF (.not.got_var(idHbbl).and.Hout(idHbbl,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idHbbl)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  ifdef AVERAGES_QUADRATIC
        IF (.not.got_var(idUUav)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idUUav)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUVav)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idUVav)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVVav)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idVVav)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        DO itrc=1,NAT
          IF (.not.got_var(idUTav(itrc))) THEN
            WRITE (stdout,60) TRIM(Vname(1,idUTav(itrc))), TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
          IF (.not.got_var(idVTav(itrc))) THEN
            WRITE (stdout,60) TRIM(Vname(1,idVTav(itrc))), TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
#  endif
#  ifdef AVERAGES_AKV
        IF (.not.got_var(idVvis)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idVvis)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  ifdef AVERAGES_AKT
        IF (.not.got_var(idTdif)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idTdif)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  ifdef AVERAGES_AKS
        IF (.not.got_var(idSdif)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idSdif)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  ifdef AVERAGES_FLUXES
        IF (.not.got_var(idTsur(itemp)).and.Hout(idTsur(itemp),ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idTsur(itemp))), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   ifdef BULK_FLUXES
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
#   endif
#   ifdef SHORTWAVE
        IF (.not.got_var(idSrad).and.Hout(idSrad,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idSrad)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   endif
#  endif
# endif
# ifdef AVERAGES_FLUXES
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
# endif
!
!  Set unlimited time record dimension to the appropriate value.
!
        IF (nrst(ng).eq.navg(ng)) THEN
          IF (ndefavg(ng).gt.0) THEN
            tavgindx(ng)=((ntstart-1)-                                  &
     &                    ndefavg(ng)*((ntstart-1)/ndefavg(ng)))/       &
     &                   navg(ng)
          ELSE
            tavgindx(ng)=(ntstart-1)/navg(ng)
          END IF
        ELSE
          tavgindx(ng)=tsize
        END IF
      END IF
!
!  Set initial average time.
!
      IF (ntsavg(ng).eq.1) THEN
        avgtime(ng)=time(ng)-0.5_r8*REAL(navg(ng),r8)*dt(ng)
      ELSE
        avgtime(ng)=time(ng)+REAL(ntsavg(ng),r8)*dt(ng)-                &
     &              0.5_r8*REAL(navg(ng),r8)*dt(ng)
      END IF
!
  10  FORMAT (6x,'DEF_AVG   - creating averages file: ',a)
  20  FORMAT (6x,'DEF_AVG   - inquiring averages file: ',a)
  30  FORMAT (/,' DEF_AVG - unable to create averages NetCDF file: ',a)
  40  FORMAT (1pe11.4,1x,'millimeter')
  50  FORMAT (/,' DEF_AVG - unable to open averages NetCDF file: ',a)
  60  FORMAT (/,' DEF_AVG - unable to find variable: ',a,2x,            &
     &         ' in averages NetCDF file: ',a)
      RETURN
      END SUBROUTINE def_avg
#else
      SUBROUTINE def_avg
      RETURN
      END SUBROUTINE def_avg
#endif
