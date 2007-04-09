#include "cppdefs.h"
#if defined TANGENT || defined TL_IOMS
      SUBROUTINE tl_def_his (ng,ldef)
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2007 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!======================================================================= 
!                                                                      !
!  This routine creates tangent linear history NetCDF file, it defines !
!  its dimensions, attributes, and variables.                          !
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
      ncname=TLMname(ng)
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
!  Create a new tangent linear history file.
!=======================================================================
!
      IF (ldef.and.OutThread) THEN
        status=nf_create(TRIM(ncname),nf_clobber,ncTLMid(ng))
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,30) TRIM(ncname)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# ifdef DISTRIBUTE
      IF (ldef) THEN
        CALL mp_bcasti (ng, iTLM, ncTLMid(ng), 1)
      END IF
# endif
!
!-----------------------------------------------------------------------
!  Define the dimensions of staggered fields.
!-----------------------------------------------------------------------
!
      IF (ldef.and.OutThread) THEN
        status=nf_def_dim(ncTLMid(ng),'xi_rho'  ,Lm(ng)+2,              &
     &                    DimIDs( 1))
        status=nf_def_dim(ncTLMid(ng),'xi_u'    ,Lm(ng)+1,              &
     &                    DimIDs( 2))
        status=nf_def_dim(ncTLMid(ng),'xi_v'    ,Lm(ng)+2,              &
     &                    DimIDs( 3))
        status=nf_def_dim(ncTLMid(ng),'xi_psi'  ,Lm(ng)+1,              &
     &                    DimIDs( 4))
        status=nf_def_dim(ncTLMid(ng),'eta_rho' ,Mm(ng)+2,              &
     &                    DimIDs( 5))
        status=nf_def_dim(ncTLMid(ng),'eta_u'   ,Mm(ng)+2,              &
     &                    DimIDs( 6))
        status=nf_def_dim(ncTLMid(ng),'eta_v'   ,Mm(ng)+1,              &
     &                    DimIDs( 7))
        status=nf_def_dim(ncTLMid(ng),'eta_psi' ,Mm(ng)+1,              &
     &                    DimIDs( 8))
# if defined WRITE_WATER && defined MASKING
        status=nf_def_dim(ncTLMid(ng),'xy_rho'  ,Nxyr(ng),              &
     &                    DimIDs(17))
        status=nf_def_dim(ncTLMid(ng),'xy_u'    ,Nxyu(ng),              &
     &                    DimIDs(18))
        status=nf_def_dim(ncTLMid(ng),'xy_v'    ,Nxyv(ng),              &
     &                    DimIDs(19))
# endif
# ifdef SOLVE3D
#  if defined WRITE_WATER && defined MASKING
        status=nf_def_dim(ncTLMid(ng),'xyz_rho' ,Nxyr(ng)*N(ng),        &
     &                    DimIDs(20))
        status=nf_def_dim(ncTLMid(ng),'xyz_u'   ,Nxyu(ng)*N(ng),        &
     &                    DimIDs(21))
        status=nf_def_dim(ncTLMid(ng),'xyz_v'   ,Nxyv(ng)*N(ng),        &
     &                    DimIDs(22))
        status=nf_def_dim(ncTLMid(ng),'xyz_w'   ,Nxyr(ng)*(N(ng)+1),    &
     &                    DimIDs(23))
#  endif
        status=nf_def_dim(ncTLMid(ng),'N'       ,N(ng),                 &
     &                    DimIDs( 9))
        status=nf_def_dim(ncTLMid(ng),'s_rho'   ,N(ng),                 &
     &                    DimIDs( 9))
        status=nf_def_dim(ncTLMid(ng),'s_w'     ,N(ng)+1,               &
     &                    DimIDs(10))
        status=nf_def_dim(ncTLMid(ng),'tracer'  ,NT(ng),                &
     &                    DimIDs(11))
#  ifdef SEDIMENT_NOT_YET
        status=nf_def_dim(ncTLMid(ng),'Nbed'    ,Nbed,                  &
     &                    DimIDs(16))
#   if defined WRITE_WATER && defined MASKING
        status=nf_def_dim(ncTLMid(ng),'xybed'   ,Nxyr(ng)*Nbed,         &
     &                    DimIDs(24))
#   endif
#  endif
#  ifdef ECOSIM
        status=nf_def_dim(ncTLMid(ng),'Nphy'    ,Nphy,                  &
     &                    DimIDs(25))
        status=nf_def_dim(ncTLMid(ng),'Nbac'    ,Nbac,                  &
     &                    DimIDs(26))
        status=nf_def_dim(ncTLMid(ng),'Ndom'    ,Ndom,                  &
     &                    DimIDs(27))
        status=nf_def_dim(ncTLMid(ng),'Nfec'    ,Nfec,                  &
     &                    DimIDs(28))
#  endif
# endif
        status=nf_def_dim(ncTLMid(ng),'boundary',4,                     &
     &                    DimIDs(14))
# ifdef FOUR_DVAR
        status=nf_def_dim(ncTLMid(ng),'Nstate  ',NstateVar(ng),         &
     &                    DimIDs(29))
# endif
        status=nf_def_dim(ncTLMid(ng),TRIM(ADJUSTL(Vname(5,idtime))),   &
     &                    nf_unlimited,DimIDs(12))
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
        tTLMindx(ng)=0
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
        status=def_info(ng,ncTLMid(ng),ncname,DimIDs)
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
          Vinfo( 3)='seconds since 0000-01-01 00:00:00'
          Vinfo( 4)='360_day'
        ELSE IF (INT(time_ref).eq.0) THEN
          Vinfo( 3)='seconds since 0000-01-01 00:00:00'
          Vinfo( 4)='365.25_day'
        ELSE IF (time_ref.gt.0.0_r8) THEN
          WRITE (Vinfo( 3),'(a,1x,a)') 'seconds since', TRIM(r_text)
          Vinfo( 4)='standard'
        END IF
        Vinfo(14)=Vname(4,idtime)
        status=def_var(ncTLMid(ng),tlmVid(idtime,ng),NF_TYPE,           &
     &                 1,recdim,Aval,Vinfo,ncname)

# ifdef PROPAGATOR
!
!  Define Ritz eigenvalues and Ritz eigenvectors Euclidean norm.
!
        Vinfo( 1)='Ritz_rvalue'
        Vinfo( 2)='real Ritz eigenvalues'
        status=def_var(ncTLMid(ng),varid,NF_TYPE,                       &
     &                 1,recdim,Aval,Vinfo,ncname)

#  if defined FT_EIGENMODES
        Vinfo( 1)='Ritz_ivalue'
        Vinfo( 2)='imaginary Ritz eigenvalues'
        status=def_var(ncTLMid(ng),varid,NF_TYPE,                       &
     &                 1,recdim,Aval,Vinfo,ncname)
#  endif

        Vinfo( 1)='Ritz_norm'
        Vinfo( 2)='Ritz eigenvectors Euclidean norm'
        status=def_var(ncTLMid(ng),varid,NF_TYPE,                       &
     &                 1,recdim,Aval,Vinfo,ncname)
# endif
!
!  Define free-surface.
!
        IF (Hout(idFsur,ng)) THEN
          Vinfo( 1)=Vname(1,idFsur)
          WRITE (Vinfo( 2),40) TRIM(Vname(2,idFsur))
          Vinfo( 3)=Vname(3,idFsur)
          Vinfo(14)=Vname(4,idFsur)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idFsur,ng),r8)
          status=def_var(ncTLMid(ng),tlmVid(idFsur,ng),NF_FOUT,         &
     &                   nvd3,t2dgrd,Aval,Vinfo,ncname)
# if defined FORWARD_WRITE && defined FORWARD_RHS
          Vinfo( 1)=Vname(1,idRzet)
          WRITE (Vinfo( 2),40) TRIM(Vname(2,idRzet))
          Vinfo( 3)=Vname(3,idRzet)
          Vinfo(14)=Vname(4,idRzet)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r2dvar,r8)
          status=def_var(ncTLMid(ng),tlmVid(idRzet,ng),NF_FOUT,         &
     &                   nvd3,t2dgrd,Aval,Vinfo,ncname)
# endif
        END IF
!
!  Define 2D momentum in the XI-direction.
!
        IF (Hout(idUbar,ng)) THEN
          Vinfo( 1)=Vname(1,idUbar)
          WRITE (Vinfo( 2),40) TRIM(Vname(2,idUbar))
          Vinfo( 3)=Vname(3,idUbar)
          Vinfo(14)=Vname(4,idUbar)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUbar,ng),r8)
          status=def_var(ncTLMid(ng),tlmVid(idUbar,ng),NF_FOUT,         &
     &                   nvd3,u2dgrd,Aval,Vinfo,ncname)
# ifdef FORWARD_WRITE
#  ifdef FORWARD_RHS
          Vinfo( 1)=Vname(1,idRu2d)
          WRITE (Vinfo( 2),40) TRIM(Vname(2,idRu2d))
          Vinfo( 3)=Vname(3,idRu2d)
          Vinfo(14)=Vname(4,idRu2d)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(u2dvar,r8)

          status=def_var(ncTLMid(ng),tlmVid(idRu2d,ng),NF_FOUT,         &
     &                   nvd3,u2dgrd,Aval,Vinfo,ncname)
#  endif
#  ifdef SOLVE3D
#   ifdef FORWARD_RHS
          Vinfo( 1)=Vname(1,idRuct)
          WRITE (Vinfo( 2),40) TRIM(Vname(2,idRuct))
          Vinfo( 3)=Vname(3,idRuct)
          Vinfo(14)=Vname(4,idRuct)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#    endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(u2dvar,r8)
          status=def_var(ncTLMid(ng),tlmVid(idRuct,ng),NF_FOUT,         &
     &                   nvd3,u2dgrd,Aval,Vinfo,ncname)
#   endif
          Vinfo( 1)=Vname(1,idUfx1)
          WRITE (Vinfo( 2),40) TRIM(Vname(2,idUfx1))
          Vinfo( 3)=Vname(3,idUfx1)
          Vinfo(14)=Vname(4,idUfx1)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(u2dvar,r8)
          status=def_var(ncTLMid(ng),tlmVid(idUfx1,ng),NF_FOUT,         &
     &                   nvd3,u2dgrd,Aval,Vinfo,ncname)
          Vinfo( 1)=Vname(1,idUfx2)
          WRITE (Vinfo( 2),40) TRIM(Vname(2,idUfx2))
          Vinfo( 3)=Vname(3,idUfx2)
          Vinfo(14)=Vname(4,idUfx2)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(u2dvar,r8)
          status=def_var(ncTLMid(ng),tlmVid(idUfx2,ng),NF_FOUT,         &
     &                   nvd3,u2dgrd,Aval,Vinfo,ncname)
#  endif
# endif
        END IF
!
!  Define 2D momentum in the ETA-direction.
!
        IF (Hout(idVbar,ng)) THEN
          Vinfo( 1)=Vname(1,idVbar)
          WRITE (Vinfo( 2),40) TRIM(Vname(2,idVbar))
          Vinfo( 3)=Vname(3,idVbar)
          Vinfo(14)=Vname(4,idVbar)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVbar,ng),r8)
          status=def_var(ncTLMid(ng),tlmVid(idVbar,ng),NF_FOUT,         &
     &                   nvd3,v2dgrd,Aval,Vinfo,ncname)
# ifdef FORWARD_WRITE
#  ifdef FORWARD_RHS
          Vinfo( 1)=Vname(1,idRv2d)
          WRITE (Vinfo( 2),40) TRIM(Vname(2,idRv2d))
          Vinfo( 3)=Vname(3,idRv2d)
          Vinfo(14)=Vname(4,idRv2d)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(v2dvar,r8)
          status=def_var(ncTLMid(ng),tlmVid(idRv2d,ng),NF_FOUT,         &
     &                   nvd3,v2dgrd,Aval,Vinfo,ncname)
#  endif
#  ifdef SOLVE3D
#   ifdef FORWARD_RHS
          Vinfo( 1)=Vname(1,idRvct)
          WRITE (Vinfo( 2),40) TRIM(Vname(2,idRvct))
          Vinfo( 3)=Vname(3,idRvct)
          Vinfo(14)=Vname(4,idRvct)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#    endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(v2dvar,r8)
          status=def_var(ncTLMid(ng),tlmVid(idRvct,ng),NF_FOUT,         &
     &                   nvd3,v2dgrd,Aval,Vinfo,ncname)
#   endif
          Vinfo( 1)=Vname(1,idVfx1)
          WRITE (Vinfo( 2),40) TRIM(Vname(2,idVfx1))
          Vinfo( 3)=Vname(3,idVfx1)
          Vinfo(14)=Vname(4,idVfx1)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(v2dvar,r8)
          status=def_var(ncTLMid(ng),tlmVid(idVfx1,ng),NF_FOUT,         &
     &                   nvd3,v2dgrd,Aval,Vinfo,ncname)
          Vinfo( 1)=Vname(1,idVfx2)
          WRITE (Vinfo( 2),40) TRIM(Vname(2,idVfx2))
          Vinfo( 3)=Vname(3,idVfx2)
          Vinfo(14)=Vname(4,idVfx2)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(v2dvar,r8)
          status=def_var(ncTLMid(ng),tlmVid(idVfx2,ng),NF_FOUT,         &
     &                   nvd3,v2dgrd,Aval,Vinfo,ncname)
#  endif
# endif
        END IF
# ifdef SOLVE3D
!
!  Define 3D momentum component in the XI-direction.
!
        IF (Hout(idUvel,ng)) THEN
          Vinfo( 1)=Vname(1,idUvel)
          WRITE (Vinfo( 2),40) TRIM(Vname(2,idUvel))
          Vinfo( 3)=Vname(3,idUvel)
          Vinfo(14)=Vname(4,idUvel)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUvel,ng),r8)
          status=def_var(ncTLMid(ng),tlmVid(idUvel,ng),NF_FOUT,         &
     &                   nvd4,u3dgrd,Aval,Vinfo,ncname)
#  if defined FORWARD_WRITE && defined FORWARD_RHS
          Vinfo( 1)=Vname(1,idRu3d)
          WRITE (Vinfo( 2),40) TRIM(Vname(2,idRu3d))
          Vinfo( 3)=Vname(3,idRu3d)
          Vinfo(14)=Vname(4,idRu3d)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(u3dvar,r8)
          status=def_var(ncTLMid(ng),tlmVid(idRu3d,ng),NF_FOUT,         &
     &                   nvd4,u3dgrd,Aval,Vinfo,ncname)
#  endif
        END IF
!
!  Define 3D momentum component in the ETA-direction.
!
        IF (Hout(idVvel,ng)) THEN
          Vinfo( 1)=Vname(1,idVvel)
          WRITE (Vinfo( 2),40) TRIM(Vname(2,idVvel))
          Vinfo( 3)=Vname(3,idVvel)
          Vinfo(14)=Vname(4,idVvel)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVvel,ng),r8)
          status=def_var(ncTLMid(ng),tlmVid(idVvel,ng),NF_FOUT,         &
     &                   nvd4,v3dgrd,Aval,Vinfo,ncname)
#  if defined FORWARD_WRITE && defined FORWARD_RHS
          Vinfo( 1)=Vname(1,idRv3d)
          WRITE (Vinfo( 2),40) TRIM(Vname(2,idRv3d))
          Vinfo( 3)=Vname(3,idRv3d)
          Vinfo(14)=Vname(4,idRv3d)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(v3dvar,r8)
          status=def_var(ncTLMid(ng),tlmVid(idRv3d,ng),NF_FOUT,         &
     &                   nvd4,v3dgrd,Aval,Vinfo,ncname)
#  endif
        END IF
!
!  Define tracer type variables.
!
        DO itrc=1,NT(ng)
          IF (Hout(idTvar(itrc),ng)) THEN
            Vinfo( 1)=Vname(1,idTvar(itrc))
            WRITE (Vinfo( 2),40) TRIM(Vname(2,idTvar(itrc)))
            Vinfo( 3)=Vname(3,idTvar(itrc))
            Vinfo(14)=Vname(4,idTvar(itrc))
            Vinfo(16)=Vname(1,idtime)
#  ifdef SEDIMENT_NOT_YET
            DO i=1,NST
              IF (itrc.eq.idsed(i)) THEN
                WRITE (Vinfo(19),50) 1000.0_r8*Sd50(i,ng)
              END IF
            END DO
#  endif
#  if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#  endif
            Vinfo(22)='coordinates'
            Aval(5)=REAL(r3dvar,r8)
            status=def_var(ncTLMid(ng),tlmTid(itrc,ng),NF_FOUT,         &
     &                     nvd4,t3dgrd,Aval,Vinfo,ncname)
          END IF
        END DO
!
!  Define density anomaly.
!
        IF (Hout(idDano,ng)) THEN
          Vinfo( 1)=Vname(1,idDano)
          WRITE (Vinfo( 2),40) TRIM(Vname(2,idDano))
          Vinfo( 3)=Vname(3,idDano)
          Vinfo(14)=Vname(4,idDano)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idDano,ng),r8)
          status=def_var(ncTLMid(ng),tlmVid(idDano,ng),NF_FOUT,         &
     &                   nvd4,t3dgrd,Aval,Vinfo,ncname)
        END IF

#  if defined FORWARD_MIXING && \
     (defined BVF_MIXING     || defined GLS_MIXING  || \
      defined LMD_MIXING     || defined MY25_MIXING)
!
!  Define vertical viscosity coefficient.
!
        IF (Hout(idVvis,ng)) THEN
          Vinfo( 1)=Vname(1,idVvis)
          WRITE (Vinfo( 2),40) TRIM(Vname(2,idVvis))
          Vinfo( 3)=Vname(3,idVvis)
          Vinfo(14)=Vname(4,idVvis)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVvis,ng),r8)
          status=def_var(ncTLMid(ng),tlmVid(idVvis,ng),NF_FOUT,         &
     &                   nvd4,w3dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define vertical diffusion coefficient for potential temperature.
!
        IF (Hout(idTdif,ng)) THEN
          Vinfo( 1)=Vname(1,idTdif)
          WRITE (Vinfo( 2),40) TRIM(Vname(2,idTdif))
          Vinfo( 3)=Vname(3,idTdif)
          Vinfo(14)=Vname(4,idTdif)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTdif,ng),r8)
          status=def_var(ncTLMid(ng),tlmVid(idTdif,ng),NF_FOUT,         &
     &                   nvd4,w3dgrd,Aval,Vinfo,ncname)
        END IF
#   ifdef SALINITY
!
!  Define vertical diffusion coefficient for salinity.
!
        IF (Hout(idSdif,ng)) THEN
          Vinfo( 1)=Vname(1,idSdif)
          WRITE (Vinfo( 2),40) TRIM(Vname(2,idSdif))
          Vinfo( 3)=Vname(3,idSdif)
          Vinfo(14)=Vname(4,idSdif)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#    endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idSdif,ng),r8)
          status=def_var(ncTLMid(ng),tlmVid(idSdif,ng),NF_FOUT,         &
     &                   nvd4,w3dgrd,Aval,Vinfo,ncname)
        END IF
#   endif
#   if defined GLS_MIXING_NOT_YET || defined MY25_MIXING_NOT_YET
!
!  Define turbulent kinetic energy.
!
        IF (Hout(idMtke,ng)) THEN
          Vinfo( 1)=Vname(1,idMtke)
          WRITE (Vinfo( 2),40) TRIM(Vname(2,idMtke))
          Vinfo( 3)=Vname(3,idMtke)
          Vinfo(14)=Vname(4,idMtke)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#    endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idMtke,ng),r8)
          status=def_var(ncTLMid(ng),tlmVid(idMtke,ng),NF_FOUT,         &
     &                   nvd4,w3dgrd,Aval,Vinfo,ncname)

          Vinfo( 1)=Vname(1,idVmKK)
          WRITE (Vinfo( 2),40) TRIM(Vname(2,idVmKK))
          Vinfo( 3)=Vname(3,idVmKK)
          Vinfo(14)=Vname(4,idVmKK)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#    endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVmKK,ng),r8)
          status=def_var(ncTLMid(ng),tlmVid(idVmKK,ng),NF_FOUT,         &
     &                   nvd4,w3dgrd,Aval,Vinfo,ncname)

        END IF
!
!  Define turbulent kinetic energy time length scale.
!
        IF (Hout(idMtls,ng)) THEN
          Vinfo( 1)=Vname(1,idMtls)
          WRITE (Vinfo( 2),40) TRIM(Vname(2,idMtls))
          Vinfo( 3)=Vname(3,idMtls)
          Vinfo(14)=Vname(4,idMtls)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#    endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idMtls,ng),r8)
          status=def_var(ncTLMid(ng),tlmVid(idMtls,ng),NF_FOUT,         &
     &                   nvd4,w3dgrd,Aval,Vinfo,ncname)

          Vinfo( 1)=Vname(1,idVmLS)
          WRITE (Vinfo( 2),40) TRIM(Vname(2,idVmLS))
          Vinfo( 3)=Vname(3,idVmLS)
          Vinfo(14)=Vname(4,idVmLS)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#    endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVmLS,ng),r8)
          status=def_var(ncTLMid(ng),tlmVid(idVmLS,ng),NF_FOUT,         &
     &                   nvd4,w3dgrd,Aval,Vinfo,ncname)
#    ifdef GLS_MIXING_NOT_YET
          Vinfo( 1)=Vname(1,idVmKP)
          WRITE (Vinfo( 2),40) TRIM(Vname(2,idVmKP))
          Vinfo( 3)=Vname(3,idVmKP)
          Vinfo(14)=Vname(4,idVmKP)
          Vinfo(16)=Vname(1,idtime)
#     if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#     endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVmKP,ng),r8)
          status=def_var(ncTLMid(ng),tlmVid(idVmKP,ng),NF_FOUT,         &
     &                   nvd4,w3dgrd,Aval,Vinfo,ncname)
#    endif
        END IF
#   endif
#  endif
#  ifdef ADJUST_STFLUX
!
!  Define surface net heat flux.
!
        DO itrc=1,NT(ng)
          IF (Hout(idTsur(itrc),ng)) THEN
            Vinfo( 1)=Vname(1,idTsur(itemp))
            WRITE (Vinfo( 2),40) TRIM(Vname(2,idTsur(itrc)))
            Vinfo( 3)=Vname(3,idTsur(itemp))
            Vinfo(14)=Vname(4,idTsur(itemp))
            Vinfo(16)=Vname(1,idtime)
            status=def_var(ncTLMid(ng),tlmVid(idTsur(itrc),ng),         &
     &                   NF_FOUT,nvd3,t2dgrd,Aval,Vinfo,ncname)
          END IF
        END DO
#  endif
# endif
# ifdef ADJUST_WSTRESS
!
!  Define surface U-momentum stress.
!
        IF (Hout(idUsms,ng)) THEN
          Vinfo( 1)=Vname(1,idUsms)
          WRITE (Vinfo( 2),40) TRIM(Vname(2,idUsms))
          Vinfo( 3)=Vname(3,idUsms)
          Vinfo(14)=Vname(4,idUsms)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#  endif
          status=def_var(ncTLMid(ng),tlmVid(idUsms,ng),NF_FOUT,         &
     &                   nvd3,u2dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define surface V-momentum stress.
!
        IF (Hout(idVsms,ng)) THEN
          Vinfo( 1)=Vname(1,idVsms)
          WRITE (Vinfo( 2),40) TRIM(Vname(2,idVsms))
          Vinfo( 2)=Vname(2,idVsms)
          Vinfo( 3)=Vname(3,idVsms)
          Vinfo(14)=Vname(4,idVsms)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#  endif
          status=def_var(ncTLMid(ng),tlmVid(idVsms,ng),NF_FOUT,         &
     &                   nvd3,v2dgrd,Aval,Vinfo,ncname)
        END IF
# endif
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
        status=nf_enddef(ncTLMid(ng))
      END IF
!
!-----------------------------------------------------------------------
!  Write out time-recordless, information variables.
!-----------------------------------------------------------------------
!
      IF (ldef) THEN
        status=wrt_info(ng, iTLM, ncTLMid(ng), OutThread, ncname)
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!=======================================================================
!  Open an existing tangent file, check its contents, and prepare for
!  appending data.
!=======================================================================
!
      IF (.not.ldef.and.OutThread) THEN
!
!  Inquire about the contents of tangent NetCDF file:  Inquire about
!  the dimensions and variables.  Check for consistency.
!
        CALL opencdf (ng, 1, ncname, fname, N(ng), 0, nrec, nvd, Vsize)
        IF (exit_flag.ne.NoError) RETURN
!
!  Open tangent linear history file for read/write.
!
        status=nf_open(TRIM(ncname),nf_write,ncTLMid(ng))
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,60) TRIM(ncname)
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
!  tangent variables. Get variable IDs.
!
        DO i=1,nvars
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idtime))) THEN
            got_var(idtime)=.TRUE.
            status=nf_inq_varid(ncTLMid(ng),TRIM(Vname(1,idtime)),      &
     &                          tlmVid(idtime,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idFsur))) THEN
            got_var(idFsur)=.TRUE.
            status=nf_inq_varid(ncTLMid(ng),TRIM(Vname(1,idFsur)),      &
     &                          tlmVid(idFsur,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idUbar))) THEN
            got_var(idUbar)=.TRUE.
            status=nf_inq_varid(ncTLMid(ng),TRIM(Vname(1,idUbar)),      &
     &                          tlmVid(idUbar,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVbar))) THEN
            got_var(idVbar)=.TRUE.
            status=nf_inq_varid(ncTLMid(ng),TRIM(Vname(1,idVbar)),      &
     &                          tlmVid(idVbar,ng))
          END IF
# ifdef FORWARD_WRITE
#  ifdef FORWARD_RHS
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idRzet))) THEN
            got_var(idRzet)=.TRUE.
            status=nf_inq_varid(ncTLMid(ng),TRIM(Vname(1,idRzet)),      &
     &                          tlmVid(idRzet,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idRu2d))) THEN
            got_var(idRu2d)=.TRUE.
            status=nf_inq_varid(ncTLMid(ng),TRIM(Vname(1,idRu2d)),      &
     &                          tlmVid(idRu2d,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idRv2d))) THEN
            got_var(idRv2d)=.TRUE.
            status=nf_inq_varid(ncTLMid(ng),TRIM(Vname(1,idRv2d)),      &
     &                          tlmVid(idRv2d,ng))
          END IF
#  endif
#  ifdef SOLVE3D
#   ifdef FORWARD_RHS
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idRuct))) THEN
            got_var(idRuct)=.TRUE.
            status=nf_inq_varid(ncTLMid(ng),TRIM(Vname(1,idRuct)),      &
     &                          tlmVid(idRuct,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idRvct))) THEN
            got_var(idRvct)=.TRUE.
            status=nf_inq_varid(ncTLMid(ng),TRIM(Vname(1,idRvct)),      &
     &                          tlmVid(idRvct,ng))
          END IF
#   endif
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idUfx1))) THEN
            got_var(idUfx1)=.TRUE.
            status=nf_inq_varid(ncTLMid(ng),TRIM(Vname(1,idUfx1)),      &
     &                          tlmVid(idUfx1,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idUfx1))) THEN
            got_var(idUfx1)=.TRUE.
            status=nf_inq_varid(ncTLMid(ng),TRIM(Vname(1,idUfx1)),      &
     &                          tlmVid(idUfx1,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idUfx2))) THEN
            got_var(idUfx2)=.TRUE.
            status=nf_inq_varid(ncTLMid(ng),TRIM(Vname(1,idUfx2)),      &
     &                          tlmVid(idUfx2,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVfx1))) THEN
            got_var(idVfx1)=.TRUE.
            status=nf_inq_varid(ncTLMid(ng),TRIM(Vname(1,idVfx1)),      &
     &                          tlmVid(idVfx1,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVfx2))) THEN
            got_var(idVfx2)=.TRUE.
            status=nf_inq_varid(ncTLMid(ng),TRIM(Vname(1,idVfx2)),      &
     &                          tlmVid(idVfx2,ng))
          END IF
#   ifdef FORWARD_RHS
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idRu3d))) THEN
            got_var(idRu3d)=.TRUE.
            status=nf_inq_varid(ncTLMid(ng),TRIM(Vname(1,idRu3d)),      &
     &                          tlmVid(idRu3d,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idRv3d))) THEN
            got_var(idRv3d)=.TRUE.
            status=nf_inq_varid(ncTLMid(ng),TRIM(Vname(1,idRv3d)),      &
     &                          tlmVid(idRv3d,ng))
          END IF
#   endif
#  endif
# endif
# ifdef SOLVE3D
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idUvel))) THEN
            got_var(idUvel)=.TRUE.
            status=nf_inq_varid(ncTLMid(ng),TRIM(Vname(1,idUvel)),      &
     &                          tlmVid(idUvel,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVvel))) THEN
            got_var(idVvel)=.TRUE.
            status=nf_inq_varid(ncTLMid(ng),TRIM(Vname(1,idVvel)),      &
     &                          tlmVid(idVvel,ng))
          END IF
          DO itrc=1,NT(ng)
            IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idTvar(itrc)))) THEN
              got_var(idTvar(itrc))=.TRUE.
              status=nf_inq_varid(ncTLMid(ng),                          &
     &                            TRIM(Vname(1,idTvar(itrc))),          &
     &                            tlmTid(itrc,ng))
            END IF
          END DO
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idDano))) THEN
            got_var(idDano)=.TRUE.
            status=nf_inq_varid(ncTLMid(ng),TRIM(Vname(1,idDano)),      &
     &                          tlmVid(idDano,ng))
          END IF
#  if defined FORWARD_MIXING     && \
     (defined GLS_MIXING_NOT_YET || defined MY25_MIXING_NOT_YET)
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVvis))) THEN
            got_var(idVvis)=.TRUE.
            status=nf_inq_varid(ncTLMid(ng),TRIM(Vname(1,idVvis)),      &
     &                          tlmVid(idVvis,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idTdif))) THEN
            got_var(idTdif)=.TRUE.
            status=nf_inq_varid(ncTLMid(ng),TRIM(Vname(1,idTdif)),      &
     &                          tlmVid(idTdif,ng))
          END IF
#   ifdef SALINITY
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idSdif))) THEN
            got_var(idSdif)=.TRUE.
            status=nf_inq_varid(ncTLMid(ng),TRIM(Vname(1,idSdif)),      &
     &                          tlmVid(idSdif,ng))
          END IF
#   endif
#   if defined GLS_MIXING_NOT_YET || defined MY25_MIXING_NOT_YET
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idMtke))) THEN
            got_var(idMtke)=.TRUE.
            status=nf_inq_varid(ncTLMid(ng),TRIM(Vname(1,idMtke)),      &
     &                          tlmVid(idMtke,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVmKK))) THEN
            got_var(idVmKK)=.TRUE.
            status=nf_inq_varid(ncTLMid(ng),TRIM(Vname(1,idVmKK)),      &
     &                          tlmVid(idVmKK,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idMtls))) THEN
            got_var(idMtls)=.TRUE.
            status=nf_inq_varid(ncTLMid(ng),TRIM(Vname(1,idMtls)),      &
     &                          tlmVid(idMtls,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVmLS))) THEN
            got_var(idVmLS)=.TRUE.
            status=nf_inq_varid(ncTLMid(ng),TRIM(Vname(1,idVmLS)),      &
     &                          tlmVid(idVmLS,ng))
          END IF
#    ifdef GLS_MIXING_NOT_YET
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVmKP))) THEN
            got_var(idVmKP)=.TRUE.
            status=nf_inq_varid(ncTLMid(ng),TRIM(Vname(1,idVmKP)),      &
     &                          tlmVid(idVmKP,ng))
          END IF
#    endif
#   endif
#  endif
#  ifdef ADJUST_STFLUX
          DO itrc=1,NT(ng)
            IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idTsur(itrc)))) THEN
              got_var(idTsur(itrc))=.TRUE.
              status=nf_inq_varid(ncTLMid(ng),                          &
     &                            TRIM(Vname(1,idTsur(itrc))),          &
     &                            tlmVid(idTsur(itrc),ng))
            END IF
          END DO
#  endif
# endif
# ifdef ADJUST_WSTRESS
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idUsms))) THEN
            got_var(idUsms)=.TRUE.
            status=nf_inq_varid(ncTLMid(ng),TRIM(Vname(1,idUsms)),      &
     &                          tlmVid(idUsms,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVsms))) THEN
            got_var(idVsms)=.TRUE.
            status=nf_inq_varid(ncTLMid(ng),TRIM(Vname(1,idVsms)),      &
     &                          tlmVid(idVsms,ng))
          END IF
# endif
        END DO
!
!  Check if tangent variables are available in input NetCDF file.
!
        IF (.not.got_var(idtime)) THEN
          WRITE (stdout,70) TRIM(Vname(1,idtime)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idFsur).and.Hout(idFsur,ng)) THEN
          WRITE (stdout,70) TRIM(Vname(1,idFsur)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbar).and.Hout(idUbar,ng)) THEN
          WRITE (stdout,70) TRIM(Vname(1,idUbar)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbar).and.Hout(idVbar,ng)) THEN
          WRITE (stdout,70) TRIM(Vname(1,idVbar)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# ifdef FORWARD_WRITE
        IF (.not.got_var(idRzet)) THEN
          WRITE (stdout,70) TRIM(Vname(1,idRzet)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idRu2d)) THEN
          WRITE (stdout,70) TRIM(Vname(1,idRu2d)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idRv2d)) THEN
          WRITE (stdout,70) TRIM(Vname(1,idRv2d)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  ifdef SOLVE3D
        IF (.not.got_var(idRuct)) THEN
          WRITE (stdout,70) TRIM(Vname(1,idRuct)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idRvct)) THEN
          WRITE (stdout,70) TRIM(Vname(1,idRvct)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUfx1)) THEN
          WRITE (stdout,70) TRIM(Vname(1,idUfx1)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUfx2)) THEN
          WRITE (stdout,70) TRIM(Vname(1,idUfx2)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVfx1)) THEN
          WRITE (stdout,70) TRIM(Vname(1,idVfx1)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVfx2)) THEN
          WRITE (stdout,70) TRIM(Vname(1,idVfx2)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idRu3d)) THEN
          WRITE (stdout,70) TRIM(Vname(1,idRu3d)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idRv3d)) THEN
          WRITE (stdout,30) TRIM(Vname(1,idRv3d)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
# endif
# ifdef SOLVE3D
        IF (.not.got_var(idUvel).and.Hout(idUvel,ng)) THEN
          WRITE (stdout,70) TRIM(Vname(1,idUvel)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVvel).and.Hout(idVvel,ng)) THEN
          WRITE (stdout,70) TRIM(Vname(1,idVvel)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        DO itrc=1,NT(ng)
          IF (.not.got_var(idTvar(itrc)).and.Hout(idTvar(itrc),ng)) THEN
            WRITE (stdout,70) TRIM(Vname(1,idTvar(itrc))), TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
        IF (.not.got_var(idDano).and.Hout(idDano,ng)) THEN
          WRITE (stdout,70) TRIM(Vname(1,idDano)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  if defined FORWARD_MIXING     && \
     (defined GLS_MIXING_NOT_YET || defined MY25_MIXING_NOT_YET)
        IF (.not.got_var(idVvis).and.Hout(idVvis,ng)) THEN
          WRITE (stdout,70) TRIM(Vname(1,idVvis)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTdif).and.Hout(idTdif,ng)) THEN
          WRITE (stdout,70) TRIM(Vname(1,idTdif)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   ifdef SALINITY
        IF (.not.got_var(idSdif).and.Hout(idSdif,ng)) THEN
          WRITE (stdout,70) TRIM(Vname(1,idSdif)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   endif
#   if defined GLS_MIXING_NOT_YET || defined MY25_MIXING_NOT_YET
        IF (.not.got_var(idMtke).and.Hout(idMtke,ng)) THEN
          WRITE (stdout,70) TRIM(Vname(1,idMtke)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVmKK).and.Hout(idVmKK,ng)) THEN
          WRITE (stdout,70) TRIM(Vname(1,idVmKK)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idMtls).and.Hout(idMtls,ng)) THEN
          WRITE (stdout,70) TRIM(Vname(1,idMtls)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVmLS).and.Hout(idVmLS,ng)) THEN
          WRITE (stdout,70) TRIM(Vname(1,idVmLS)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#    ifdef GSL_MIXING
        IF (.not.got_var(idVmKP).and.Hout(idVmKP,ng)) THEN
          WRITE (stdout,70) TRIM(Vname(1,idVmKP)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#    endif
#   endif
#  endif
#  ifdef ADJUST_STFLUX
        DO itrc=1,NT(ng)
          IF (.not.got_var(idTsur(itrc)).and.Hout(idTsur(itrc),ng)) THEN
            WRITE (stdout,70) TRIM(Vname(1,idTsur(itrc))), TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
#  endif
# endif
# ifdef ADJUST_WSTRESS
        IF (.not.got_var(idUsms).and.Hout(idUsms,ng)) THEN
          WRITE (stdout,70) TRIM(Vname(1,idUsms)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVsms).and.Hout(idVsms,ng)) THEN
          WRITE (stdout,70) TRIM(Vname(1,idVsms)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
!
!  Set unlimited time record dimension to the appropriate value.
!
!!      tTLMindx(ng)=tsize
        IF (ndefTLM(ng).gt.0) THEN
          IF ((ntstart(ng)-1).lt.ndefTLM(ng)) THEN
            tTLMindx(ng)=1+((ntstart(ng)-1)-                            &
     &                       ndefTLM(ng)*((ntstart(ng)-1)/ndefTLM(ng)))/&
     &                     nTLM(ng)
          ELSE
            tTLMindx(ng)=((ntstart(ng)-1)-                              &
     &                     ndefTLM(ng)*((ntstart(ng)-1)/ndefTLM(ng)))/  &
     &                   nTLM(ng)
          END IF
        ELSE
          tTLMindx(ng)=(ntstart(ng)-1)/nTLM(ng)
        END IF
        tTLMindx(ng)=MIN(tTLMindx(ng),tsize)
      END IF
!
  10  FORMAT (3x,'TL_DEF_HIS   - creating tangent  file: ',a)
  20  FORMAT (3x,'TL_DEF_HIS   - inquiring tangent file: ',a)
  30  FORMAT (/,' TL_DEF_HIS - unable to create tangent NetCDF file: ', &
     &        a)
  40  FORMAT ('tangent linear',1x,a)
  50  FORMAT (1pe11.4,1x,'millimeter')
  60  FORMAT (/,' TL_DEF_HIS - unable to open tangent NetCDF file: ',a)
  70  FORMAT (/,' TL_DEF_HIS - unable to find variable: ',a,2x,         &
     &        ' in tangent NetCDF file: ',a)
      RETURN
      END SUBROUTINE tl_def_his
#else
      SUBROUTINE tl_def_his
      RETURN
      END SUBROUTINE tl_def_his
#endif
