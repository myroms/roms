#include "cppdefs.h"
#ifdef ADJOINT
      SUBROUTINE ad_def_his (ng,ldef)
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2007 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine creates adjoint history NetCDF file, it defines its    !
!  dimensions, attributes, and variables.                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
# ifdef FOUR_DVAR
      USE mod_fourdvar
# endif
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
# if defined SEDIMENT_NOT_YET || defined BBL_MODEL_NOT_YET
      USE mod_sediment
# endif
# ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcasti
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
      logical :: got_var(NV)

      integer, parameter :: Natt = 24

      integer :: i, j, itrc, nrec, nvd, nvd3, nvd4
      integer :: recdim, status, varid

      integer :: DimIDs(29), t2dgrd(3), u2dgrd(3), v2dgrd(3)
      integer :: Vsize(4)

      integer :: def_info, def_var, wrt_info

# ifdef SOLVE3D
      integer :: t3dgrd(4), u3dgrd(4), v3dgrd(4), w3dgrd(4)
# endif

      real(r8) :: Aval(6)

      character (len=80) :: Vinfo(Natt)
      character (len=80) :: fname, ncname
!
!-----------------------------------------------------------------------
!  Set and report file name.
!-----------------------------------------------------------------------
!
      IF (exit_flag.ne.NoError) RETURN
      ncname=ADJname(ng)
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
!  Create a new adjoint history file.
!=======================================================================
!
      IF (ldef.and.OutThread) THEN
        status=nf90_create(TRIM(ncname),nf90_clobber,ncADJid(ng))
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,30) TRIM(ncname)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# ifdef DISTRIBUTE
      IF (ldef) THEN
        CALL mp_bcasti (ng, iADM, ncADJid(ng), 1)
      END IF
# endif
!
!-----------------------------------------------------------------------
!  Define the dimensions of staggered fields.
!-----------------------------------------------------------------------
!
      IF (ldef.and.OutThread) THEN
        status=nf90_def_dim(ncADJid(ng),'xi_rho'  ,Lm(ng)+2,            &
     &                      DimIDs( 1))
        status=nf90_def_dim(ncADJid(ng),'xi_u'    ,Lm(ng)+1,            &
     &                      DimIDs( 2))
        status=nf90_def_dim(ncADJid(ng),'xi_v'    ,Lm(ng)+2,            &
     &                      DimIDs( 3))
        status=nf90_def_dim(ncADJid(ng),'xi_psi'  ,Lm(ng)+1,            &
     &                      DimIDs( 4))
        status=nf90_def_dim(ncADJid(ng),'eta_rho' ,Mm(ng)+2,            &
     &                      DimIDs( 5))
        status=nf90_def_dim(ncADJid(ng),'eta_u'   ,Mm(ng)+2,            &
     &                      DimIDs( 6))
        status=nf90_def_dim(ncADJid(ng),'eta_v'   ,Mm(ng)+1,            &
     &                      DimIDs( 7))
        status=nf90_def_dim(ncADJid(ng),'eta_psi' ,Mm(ng)+1,            &
     &                      DimIDs( 8))
# if defined WRITE_WATER && defined MASKING
        status=nf90_def_dim(ncADJid(ng),'xy_rho'  ,Nxyr(ng),            &
     &                      DimIDs(17))
        status=nf90_def_dim(ncADJid(ng),'xy_u'    ,Nxyu(ng),            &
     &                      DimIDs(18))
        status=nf90_def_dim(ncADJid(ng),'xy_v'    ,Nxyv(ng),            &
     &                      DimIDs(19))
# endif
# ifdef SOLVE3D
#  if defined WRITE_WATER && defined MASKING
        status=nf90_def_dim(ncADJid(ng),'xyz_rho' ,Nxyr(ng)*N(ng),      &
     &                      DimIDs(20))
        status=nf90_def_dim(ncADJid(ng),'xyz_u'   ,Nxyu(ng)*N(ng),      &
     &                      DimIDs(21))
        status=nf90_def_dim(ncADJid(ng),'xyz_v'   ,Nxyv(ng)*N(ng),      &
     &                      DimIDs(22))
        status=nf90_def_dim(ncADJid(ng),'xyz_w'   ,Nxyr(ng)*(N(ng)+1),  &
     &                      DimIDs(23))
#  endif
        status=nf90_def_dim(ncADJid(ng),'N'       ,N(ng),               &
     &                      DimIDs( 9))
        status=nf90_def_dim(ncADJid(ng),'s_rho'   ,N(ng),               &
     &                      DimIDs( 9))
        status=nf90_def_dim(ncADJid(ng),'s_w'     ,N(ng)+1,             &
     &                      DimIDs(10))
        status=nf90_def_dim(ncADJid(ng),'tracer'  ,NT(ng),              &
     &                      DimIDs(11))
#  ifdef SEDIMENT_NOT_YET
        status=nf90_def_dim(ncADJid(ng),'Nbed'    ,Nbed,                &
     &                      DimIDs(16))
#   if defined WRITE_WATER && defined MASKING
        status=nf90_def_dim(ncADJid(ng),'xybed'   ,Nxyr(ng)*Nbed,       &
     &                      DimIDs(24))
#   endif
#  endif
#  ifdef ECOSIM
        status=nf90_def_dim(ncADJid(ng),'Nphy'    ,Nphy,                &
     &                      DimIDs(25))
        status=nf90_def_dim(ncADJid(ng),'Nbac'    ,Nbac,                &
     &                      DimIDs(26))
        status=nf90_def_dim(ncADJid(ng),'Ndom'    ,Ndom,                &
     &                      DimIDs(27))
        status=nf90_def_dim(ncADJid(ng),'Nfec'    ,Nfec,                &
     &                      DimIDs(28))
#  endif
# endif
        status=nf90_def_dim(ncADJid(ng),'boundary',4,                   &
     &                      DimIDs(14))
# ifdef FOUR_DVAR
        status=nf90_def_dim(ncADJid(ng),'Nstate  ',NstateVar(ng),       &
     &                      DimIDs(29))
# endif
        status=nf90_def_dim(ncADJid(ng),TRIM(ADJUSTL(Vname(5,idtime))), &
     &                      nf90_unlimited,DimIDs(12))
        recdim=DimIDs(12)
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
# endif
!
!  Initialize unlimited time record dimension.
!
        tADJindx(ng)=0
!
!  Initialize local information variable arrays.
!
        DO i=1,Natt
          DO j=1,80
            Vinfo(i)(j:j)=' '
          END DO
        END DO
        DO i=1,6
          Aval(i)=0.0_r8
        END DO
!
!-----------------------------------------------------------------------
!  Define time-recordless information variables.
!-----------------------------------------------------------------------
!
        status=def_info(ng,ncADJid(ng),ncname,DimIDs)
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
        status=def_var(ncADJid(ng),adjVid(idtime,ng),NF_TYPE,           &
     &                 1,recdim,Aval,Vinfo,ncname)

# ifdef PROPAGATOR
!
!  Define Ritz eigenvalues and Ritz eigenvectors Euclidean norm.
!
        Vinfo( 1)='Ritz_rvalue'
        Vinfo( 2)='Ritz eigenvalues'
        status=def_var(ncADJid(ng),varid,NF_TYPE,                       &
     &                 1,recdim,Aval,Vinfo,ncname)

#  if defined AFT_EIGENMODES
        Vinfo( 1)='Ritz_ivalue'
        Vinfo( 2)='imaginary Ritz eigenvalues'
        status=def_var(ncADJid(ng),varid,NF_TYPE,                       &
     &                 1,recdim,Aval,Vinfo,ncname)
#  endif

        Vinfo( 1)='Ritz_norm'
        Vinfo( 2)='Ritz eigenvectors Euclidean norm'
        status=def_var(ncADJid(ng),varid,NF_TYPE,                       &
     &                 1,recdim,Aval,Vinfo,ncname)
# endif
!
!  Define bathymetry.
!
        Vinfo( 1)=Vname(1,idbath)
        WRITE (Vinfo( 2),'(a,1x,a)') 'adjoint', TRIM(Vname(2,idbath))
        Vinfo( 3)='meter-1'
        Vinfo(14)=Vname(4,idbath)
        Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
# endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idbath,ng),r8)
        status=def_var(ncADJid(ng),adjVid(idbath,ng),NF_FOUT,           &
     &                 nvd3,t2dgrd,Aval,Vinfo,ncname)
!
!  Define free-surface.
!
        IF (Hout(idFsur,ng)) THEN
          Vinfo( 1)=Vname(1,idFsur)
          WRITE (Vinfo( 2),'(a,1x,a)') 'adjoint', TRIM(Vname(2,idFsur))
          Vinfo( 3)='meter-1'
          Vinfo(14)=Vname(4,idFsur)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idFsur,ng),r8)
          status=def_var(ncADJid(ng),adjVid(idFsur,ng),NF_FOUT,         &
     &                   nvd3,t2dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define 2D momentum in the XI-direction.
!
        IF (Hout(idUbar,ng)) THEN
          Vinfo( 1)=Vname(1,idUbar)
          WRITE (Vinfo( 2),'(a,1x,a)') 'adjoint', TRIM(Vname(2,idUbar))
          Vinfo( 3)='second meter-1'
          Vinfo(14)=Vname(4,idUbar)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUbar,ng),r8)
          status=def_var(ncADJid(ng),adjVid(idUbar,ng),NF_FOUT,         &
     &                   nvd3,u2dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define 2D momentum in the ETA-direction.
!
        IF (Hout(idVbar,ng)) THEN
          Vinfo( 1)=Vname(1,idVbar)
          WRITE (Vinfo( 2),'(a,1x,a)') 'adjoint', TRIM(Vname(2,idVbar))
          Vinfo( 3)='second meter-1'
          Vinfo(14)=Vname(4,idVbar)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVbar,ng),r8)
          status=def_var(ncADJid(ng),adjVid(idVbar,ng),NF_FOUT,         &
     &                   nvd3,v2dgrd,Aval,Vinfo,ncname)
        END IF
# ifdef SOLVE3D
!
!  Define 3D momentum component in the XI-direction.
!
        IF (Hout(idUvel,ng)) THEN
          Vinfo( 1)=Vname(1,idUvel)
          WRITE (Vinfo( 2),'(a,1x,a)') 'adjoint', TRIM(Vname(2,idUvel))
          Vinfo( 3)='second meter-1'
          Vinfo(14)=Vname(4,idUvel)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUvel,ng),r8)
          status=def_var(ncADJid(ng),adjVid(idUvel,ng),NF_FOUT,         &
     &                   nvd4,u3dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define 3D momentum component in the ETA-direction.
!
        IF (Hout(idVvel,ng)) THEN
          Vinfo( 1)=Vname(1,idVvel)
          WRITE (Vinfo( 2),'(a,1x,a)') 'adjoint', TRIM(Vname(2,idVvel))
          Vinfo( 3)='second meter-1'
          Vinfo(14)=Vname(4,idVvel)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVvel,ng),r8)
          status=def_var(ncADJid(ng),adjVid(idVvel,ng),NF_FOUT,         &
     &                   nvd4,v3dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define tracer type variables.
!
        DO itrc=1,NT(ng)
          IF (Hout(idTvar(itrc),ng)) THEN
            Vinfo( 1)=Vname(1,idTvar(itrc))
            WRITE (Vinfo( 2),'(a,1x,a)') 'adjoint',                     &
     &                                   TRIM(Vname(2,idTvar(itrc)))
            Vinfo( 3)=Vname(3,idTvar(itrc))
            Vinfo(14)=Vname(4,idTvar(itrc))
            Vinfo(16)=Vname(1,idtime)
#  ifdef SEDIMENT_NOT_YET
            DO i=1,NST
              IF (itrc.eq.idsed(i)) THEN
                WRITE (Vinfo(19),40) 1000.0_r8*Sd50(i,ng)
              END IF
            END DO
#  endif
#  if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#  endif
            Vinfo(22)='coordinates'
            Aval(5)=REAL(r3dvar,r8)
            status=def_var(ncADJid(ng),adjTid(itrc,ng),NF_FOUT,         &
     &                     nvd4,t3dgrd,Aval,Vinfo,ncname)
          END IF
        END DO
!
!  Define density anomaly.
!
        IF (Hout(idDano,ng)) THEN
          Vinfo( 1)=Vname(1,idDano)
          WRITE (Vinfo( 2),'(a,1x,a)') 'adjoint', TRIM(Vname(2,idDano))
          Vinfo( 3)=Vname(3,idDano)
          Vinfo(14)=Vname(4,idDano)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idDano,ng),r8)
          status=def_var(ncADJid(ng),adjVid(idDano,ng),NF_FOUT,         &
     &                   nvd4,t3dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define vertical viscosity coefficient.
!
        IF (Hout(idVvis,ng)) THEN
          Vinfo( 1)=Vname(1,idVvis)
          WRITE (Vinfo( 2),'(a,1x,a)') 'adjoint', TRIM(Vname(2,idVvis))
          Vinfo( 2)=Vname(2,idVvis)
          Vinfo( 3)=Vname(3,idVvis)
          Vinfo(14)=Vname(4,idVvis)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVvis,ng),r8)
          status=def_var(ncADJid(ng),adjVid(idVvis,ng),NF_FOUT,         &
     &                   nvd4,w3dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define vertical diffusion coefficient for potential temperature.
!
        IF (Hout(idTdif,ng)) THEN
          Vinfo( 1)=Vname(1,idTdif)
          WRITE (Vinfo( 2),'(a,1x,a)') 'adjoint', TRIM(Vname(2,idTdif))
          Vinfo( 2)=Vname(2,idTdif)
          Vinfo( 3)=Vname(3,idTdif)
          Vinfo(14)=Vname(4,idTdif)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTdif,ng),r8)
          status=def_var(ncADJid(ng),adjVid(idTdif,ng),NF_FOUT,         &
     &                   nvd4,w3dgrd,Aval,Vinfo,ncname)
        END IF
#  ifdef SALINITY
!
!  Define vertical diffusion coefficient for salinity.
!
        IF (Hout(idSdif,ng)) THEN
          Vinfo( 1)=Vname(1,idSdif)
          WRITE (Vinfo( 2),'(a,1x,a)') 'adjoint', TRIM(Vname(2,idSdif))
          Vinfo( 2)=Vname(2,idSdif)
          Vinfo( 3)=Vname(3,idSdif)
          Vinfo(14)=Vname(4,idSdif)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idSdif,ng),r8)
          status=def_var(ncADJid(ng),adjVid(idSdif,ng),NF_FOUT,         &
     &                   nvd4,w3dgrd,Aval,Vinfo,ncname)
        END IF
#  endif
!
!  Define surface net heat flux.
!
        IF (Hout(idTsur(itemp),ng)) THEN
          Vinfo( 1)=Vname(1,idTsur(itemp))
          WRITE (Vinfo( 2),'(a,1x,a)') 'adjoint',                       &
     &                                 TRIM(Vname(2,idTsur(itemp)))
          Vinfo( 3)=Vname(3,idTsur(itemp))
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idTsur(itemp))
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ncADJid(ng),adjVid(idTsur(itemp),ng),          &
     &                   NF_FOUT,nvd3,t2dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define surface net salt flux.
!
        IF (Hout(idTsur(isalt),ng)) THEN
          Vinfo( 1)=Vname(1,idTsur(isalt))
          WRITE (Vinfo( 2),'(a,1x,a)') 'adjoint',                       &
     &                                 TRIM(Vname(2,idTsur(isalt)))
          Vinfo( 3)=Vname(3,idTsur(isalt))
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idTsur(itemp))
          Vinfo(16)=Vname(1,idtime)
          status=def_var(ncADJid(ng),adjVid(idTsur(isalt),ng),          &
     &                   NF_FOUT,nvd3,t2dgrd,Aval,Vinfo,ncname)
        END IF
# endif
!
!  Define surface U-momentum stress.
!
        IF (Hout(idUsms,ng)) THEN
          Vinfo( 1)=Vname(1,idUsms)
          WRITE (Vinfo( 2),'(a,1x,a)') 'adjoint', TRIM(Vname(2,idUsms))
          Vinfo( 3)=Vname(3,idUsms)
          Vinfo(14)=Vname(4,idUsms)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
# endif
          status=def_var(ncADJid(ng),adjVid(idUsms,ng),NF_FOUT,         &
     &                   nvd3,u2dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define surface V-momentum stress.
!
        IF (Hout(idVsms,ng)) THEN
          Vinfo( 1)=Vname(1,idVsms)
          WRITE (Vinfo( 2),'(a,1x,a)') 'adjoint', TRIM(Vname(2,idVsms))
          Vinfo( 2)=Vname(2,idVsms)
          Vinfo( 3)=Vname(3,idVsms)
          Vinfo(14)=Vname(4,idVsms)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
# endif
          status=def_var(ncADJid(ng),adjVid(idVsms,ng),NF_FOUT,         &
     &                   nvd3,v2dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define bottom U-momentum stress.
!
        IF (Hout(idUbms,ng)) THEN
          Vinfo( 1)=Vname(1,idUbms)
          WRITE (Vinfo( 2),'(a,1x,a)') 'adjoint', TRIM(Vname(2,idUbms))
          Vinfo( 3)=Vname(3,idUbms)
          Vinfo(14)=Vname(4,idUbms)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
# endif
          status=def_var(ncADJid(ng),adjVid(idUbms,ng),NF_FOUT,         &
     &                   nvd3,u2dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define bottom V-momentum stress.
!
        IF (Hout(idVbms,ng)) THEN
          Vinfo( 1)=Vname(1,idVbms)
          WRITE (Vinfo( 2),'(a,1x,a)') 'adjoint', TRIM(Vname(2,idVbms))
          Vinfo( 3)=Vname(3,idVbms)
          Vinfo(14)=Vname(4,idVbms)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
# endif
          status=def_var(ncADJid(ng),adjVid(idVbms,ng),NF_FOUT,         &
     &                   nvd3,v2dgrd,Aval,Vinfo,ncname)
        END IF
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
        status=nf90_enddef(ncADJid(ng))
      END IF
!
!-----------------------------------------------------------------------
!  Write out time-recordless, information variables.
!-----------------------------------------------------------------------
!
      IF (ldef) THEN
        status=wrt_info(ng, iADM, ncADJid(ng), OutThread, ncname)
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!=======================================================================
!  Open an existing adjoint file, check its contents, and prepare for
!  appending data.
!=======================================================================
!
      IF (.not.ldef.and.OutThread) THEN
!
!  Inquire about the contents of adjoint NetCDF file:  Inquire about
!  the dimensions and variables.  Check for consistency.
!
        CALL opencdf (ng, 1, ncname, fname, N(ng), 0, nrec, nvd, Vsize)
        IF (exit_flag.ne.NoError) RETURN
!
!  Open history file for read/write.
!
        status=nf90_open(TRIM(ncname),nf90_write,ncADJid(ng))
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,50) TRIM(ncname)
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
!  adjoint variables. Get variable IDs.
!
        DO i=1,nvars
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idtime))) THEN
            got_var(idtime)=.TRUE.
            status=nf90_inq_varid(ncADJid(ng),TRIM(Vname(1,idtime)),    &
     &                            adjVid(idtime,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idFsur))) THEN
            got_var(idFsur)=.TRUE.
            status=nf90_inq_varid(ncADJid(ng),TRIM(Vname(1,idFsur)),    &
     &                            adjVid(idFsur,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idUbar))) THEN
            got_var(idUbar)=.TRUE.
            status=nf90_inq_varid(ncADJid(ng),TRIM(Vname(1,idUbar)),    &
     &                            adjVid(idUbar,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVbar))) THEN
            got_var(idVbar)=.TRUE.
            status=nf90_inq_varid(ncADJid(ng),TRIM(Vname(1,idVbar)),    &
     &                            adjVid(idVbar,ng))
          END IF
# ifdef SOLVE3D
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idUvel))) THEN
            got_var(idUvel)=.TRUE.
            status=nf90_inq_varid(ncADJid(ng),TRIM(Vname(1,idUvel)),    &
     &                            adjVid(idUvel,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVvel))) THEN
            got_var(idVvel)=.TRUE.
            status=nf90_inq_varid(ncADJid(ng),TRIM(Vname(1,idVvel)),    &
     &                            adjVid(idVvel,ng))
          END IF
          DO itrc=1,NT(ng)
            IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idTvar(itrc)))) THEN
              got_var(idTvar(itrc))=.TRUE.
              status=nf90_inq_varid(ncADJid(ng),                        &
     &                              TRIM(Vname(1,idTvar(itrc))),        &
     &                              adjTid(itrc,ng))
            END IF
          END DO
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idDano))) THEN
            got_var(idDano)=.TRUE.
            status=nf90_inq_varid(ncADJid(ng),TRIM(Vname(1,idDano)),    &
     &                            adjVid(idDano,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idTsur(itemp)))) THEN
            got_var(idTsur(itemp))=.TRUE.
            status=nf90_inq_varid(ncADJid(ng),                          &
     &                            TRIM(Vname(1,idTsur(itemp))),         &
     &                            adjVid(idTsur(itemp),ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVvis))) THEN
            got_var(idVvis)=.TRUE.
            status=nf90_inq_varid(ncADJid(ng),TRIM(Vname(1,idVvis)),    &
     &                            adjVid(idVvis,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idTdif))) THEN
            got_var(idTdif)=.TRUE.
            status=nf90_inq_varid(ncADJid(ng),TRIM(Vname(1,idTdif)),    &
     &                            adjVid(idTdif,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idSdif))) THEN
            got_var(idSdif)=.TRUE.
            status=nf90_inq_varid(ncADJid(ng),TRIM(Vname(1,idSdif)),    &
     &                            adjVid(idSdif,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idTsur(isalt)))) THEN
            got_var(idTsur(isalt))=.TRUE.
            status=nf90_inq_varid(ncADJid(ng),                          &
     &                            TRIM(Vname(1,idTsur(isalt))),         &
     &                            adjVid(idTsur(isalt),ng))
          END IF
# endif
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idUsms))) THEN
            got_var(idUsms)=.TRUE.
            status=nf90_inq_varid(ncADJid(ng),TRIM(Vname(1,idUsms)),    &
     &                            adjVid(idUsms,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVsms))) THEN
            got_var(idVsms)=.TRUE.
            status=nf90_inq_varid(ncADJid(ng),TRIM(Vname(1,idVsms)),    &
     &                            adjVid(idVsms,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idUbms))) THEN
            got_var(idUbms)=.TRUE.
            status=nf90_inq_varid(ncADJid(ng),TRIM(Vname(1,idUbms)),    &
     &                            adjVid(idUbms,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVbms))) THEN
            got_var(idVbms)=.TRUE.
            status=nf90_inq_varid(ncADJid(ng),TRIM(Vname(1,idVbms)),    &
     &                            adjVid(idVbms,ng))
          END IF
        END DO
!
!  Check if adjoint variables are available in input NetCDF file.
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
        IF (.not.got_var(idTsur(itemp)).and.Hout(idTsur(itemp),ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idTsur(itemp))), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTsur(isalt)).and.Hout(idTsur(isalt),ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idTsur(isalt))), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
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
!
!  Set unlimited time record dimension to the appropriate value.
!
!!      tADJindx(ng)=tsize
        IF (ndefADJ(ng).gt.0) THEN
          IF ((ntstart(ng)-1).lt.ndefADJ(ng)) THEN
            tADJindx(ng)=1+((ntstart(ng)-1)-                            &
     &                       ndefADJ(ng)*((ntstart(ng)-1)/ndefADJ(ng)))/&
     &                     nADJ(ng)
          ELSE
            tADJindx(ng)=((ntstart(ng)-1)-                              &
     &                     ndefADJ(ng)*((ntstart(ng)-1)/ndefADJ(ng)))/  &
     &                   nADJ(ng)
          END IF
        ELSE
          tADJindx(ng)=(ntstart(ng)-1)/nADJ(ng)
        END IF
        tADJindx(ng)=MIN(tADJindx(ng),tsize)
      END IF
!
  10  FORMAT (3x,'AD_DEF_HIS   - creating adjoint  file: ',a)
  20  FORMAT (3x,'AD_DEF_HIS   - inquiring adjoint file: ',a)
  30  FORMAT (/,' AD_DEF_HIS - unable to create adjoint NetCDF file: ', &
     &        a)
  40  FORMAT (1pe11.4,1x,'millimeter')
  50  FORMAT (/,' AD_DEF_HIS - unable to open adjoint NetCDF file: ',a)
  60  FORMAT (/,' AD_DEF_HIS - unable to find variable: ',a,2x,         &
     &        ' in adjoint NetCDF file: ',a)
      RETURN
      END SUBROUTINE ad_def_his
#else
      SUBROUTINE ad_def_his
      RETURN
      END SUBROUTINE ad_def_his
#endif
