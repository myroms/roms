#include "cppdefs.h"
      MODULE ice_output_mod
#if defined ICE_MODEL && defined SOLVE3D
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module defines/writes the ice model prognostic variables into  !
!  output NetCDF files.                                                !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_grid
      USE mod_ice
      USE mod_ncparam
      USE mod_stepping
!
      USE def_var_mod,     ONLY : def_var
      USE nf_fwrite2d_mod, ONLY : nf_fwrite2d
      USE uv_rotate_mod,   ONLY : uv_rotate2d
      USE strings_mod,     ONLY : FoundError
!
      implicit none
!
      PUBLIC :: ice_def_nf90
# if defined PIO_LIB && defined DISTRIBUTE
      PUBLIC :: ice_def_pio
# endif
      PUBLIC :: ice_wrt_nf90
# if defined PIO_LIB && defined DISTRIBUTE
      PUBLIC :: ice_wrt_pio
# endif
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE ice_def_nf90 (ng, model, ldef, VarOut, S,              &
     &                         t2dgrd, u2dgrd, v2dgrd)
!***********************************************************************
!
      USE mod_netcdf
!
!  Imported variable declarations.
!
      logical, intent(in) :: ldef, VarOut(NV,Ngrids)
!
      integer, intent(in) :: ng, model
      integer, intent(in), optional :: t2dgrd(:), u2dgrd(:), v2dgrd(:)
!
      TYPE(T_IO), intent(inout) :: S(Ngrids)
!
!  Local variable declarations.
!
      logical :: LdefVar
!
      integer, parameter :: Natt = 25

      integer :: i, ifield, j, nf, nvd3, status

      integer :: icegrd(3)
!
      real(r8) :: Aval(6)
!
      character (len=120) :: Vinfo(Natt)
      character (len=256) :: ncname
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", ice_def_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Define prognostic ice variables into specified output NetCDF.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ncname=S(ng)%name
!
!  Set number of dimensions for output variables.
!
# if defined WRITE_WATER && defined MASKING
        nvd3=2
# else
        nvd3=3
# endif
!
!  Initialize local information variable arrays.
!
      DO i=1,Natt
        DO j=1,LEN(Vinfo(1))
          Vinfo(i)(j:j)=' '
        END DO
      END DO
      DO i=1,6
        Aval(i)=0.0_r8
      END DO
!
!-----------------------------------------------------------------------
!  Define ice model state prognostic variables.
!-----------------------------------------------------------------------
!
      DEFINE : IF (ldef) THEN

        DO nf=1,nIceS
          IF (iSice(nf).gt.0) THEN
            ifield=iSice(nf)
            IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
              LdefVar=.TRUE.
            ELSE
              LdefVar=VarOut(ifield,ng)
            END IF    
            IF (Ldefvar) THEN
              Vinfo( 1)=Vname(1,ifield)
              Vinfo( 2)=Vname(2,ifield)
              Vinfo( 3)=Vname(3,ifield)
              Vinfo(14)=Vname(4,ifield)
              Vinfo(16)=Vname(1,idtime)
              Vinfo(21)=Vname(6,ifield)
              Vinfo(22)='coordinates'
              Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
!
              SELECT CASE (nf)
                CASE (isUice)
                  icegrd(1:3)=u2dgrd(1:3)
# if defined WRITE_WATER && defined MASKING
                  Vinfo(20)='mask_u'
# endif
                CASE (isVice)
                  icegrd(1:3)=v2dgrd(1:3)
# if defined WRITE_WATER && defined MASKING
                  Vinfo(20)='mask_v'
# endif
                CASE DEFAULT
                  icegrd(1:3)=t2dgrd(1:3)
# if defined WRITE_WATER && defined MASKING
                  Vinfo(20)='mask_rho'
# endif
              END SELECT
!
              status=def_var(ng, model, S(ng)%ncid,                     &
     &                       S(ng)%Vid(ifield),                         &
     &                       NF_FOUT, nvd3, icegrd, Aval, Vinfo, ncname)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END IF
          END IF
        END DO
!
!  Define ice EASTward velocity at RHO-points.
!
        IF (VarOut(idUiER,ng)) THEN
          Vinfo( 1)=Vname(1,idUiER)
          Vinfo( 2)=Vname(2,idUiER)
          Vinfo( 3)=Vname(3,idUiER)
          Vinfo(14)=Vname(4,idUiER)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(21)=Vname(6,idUiER)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUiER,ng),r8)
          status=def_var(ng, model, S(ng)%ncid,                         &
     &                   S(ng)%Vid(idUiER),                             &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define ice NORTHward velocity at RHO-points.
!
        IF (VarOut(idViNR,ng)) THEN
          Vinfo( 1)=Vname(1,idViNR)
          Vinfo( 2)=Vname(2,idViNR)
          Vinfo( 3)=Vname(3,idViNR)
          Vinfo(14)=Vname(4,idViNR)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(21)=Vname(6,idViNR)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idViNR,ng),r8)
          status=def_var(ng, model, S(ng)%ncid,                         &
     &                   S(ng)%Vid(idViNR),                             &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!-----------------------------------------------------------------------
!  Define ice model internal variables.
!-----------------------------------------------------------------------
!
        DO nf=1,nIceF
          IF (iFice(nf).gt.0) THEN
            ifield=iFice(nf)
            IF (VarOut(ifield,ng)) THEN
              Vinfo( 1)=Vname(1,ifield)
              Vinfo( 2)=Vname(2,ifield)
              Vinfo( 3)=Vname(3,ifield)
              Vinfo(14)=Vname(4,ifield)
              Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
              Vinfo(20)='mask_rho'
# endif
              Vinfo(21)=Vname(6,ifield)
              Vinfo(22)='coordinates'
              Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
!
              SELECT CASE (nf)
                CASE (icWdiv)
                  Vinfo(11)='increase ice thickness'
                  Vinfo(12)='decrease ice concentration'
                CASE (icW_ai)
                  Vinfo(11)='freezing'
                  Vinfo(12)='melting'
                CASE (icW_ao, icW_io)
                  Vinfo(11)='melting'
                  Vinfo(12)='freezing'
              END SELECT
!
              status=def_var(ng, model, S(ng)%ncid,                     &
     &                       S(ng)%Vid(ifield),                         &
     &                       NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END IF
          END IF
        END DO
      END IF DEFINE
!
!-----------------------------------------------------------------------
!  Otherwise, check existing output file and prepare for appending
!  data.
!-----------------------------------------------------------------------
!
      QUERY : IF (.not.ldef) THEN
!
!  Initialize locallogical switches.
!
        DO i=1,NV
          got_var(i)=.FALSE.
        END DO
!
!  Scan variable list from input NetCDF and activate switches for
!  sea-ice variables. Get variable IDs.
!
        DO i=1,n_var
          DO nf=1,nIceS
            IF (iSice(nf).gt.0) THEN
              ifield=iSice(nf)
              IF (TRIM(var_name(i)).eq.TRIM(Vname(1,ifield))) THEN
                got_var(ifield)=.TRUE.
                S(ng)%Vid(ifield)=var_id(i)
                CYCLE
              END IF
            END IF
          END DO          
!
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUiER))) THEN
            got_var(idUiER)=.TRUE.
            S(ng)%Vid(idUiER)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idViER))) THEN
            got_var(idViER)=.TRUE.
            S(ng)%Vid(idViER)=var_id(i)
          END IF
!
          DO nf=1,nIceF
            IF (iFice(nf).gt.0) THEN
              ifield=iFice(nf)
              IF (TRIM(var_name(i)).eq.TRIM(Vname(1,ifield))) THEN
                got_var(ifield)=.TRUE.
                S(ng)%Vid(ifield)=var_id(i)
                CYCLE
              END IF
            END IF
          END DO
        END DO
!
!  Check if output variables are available in input NetCDF file.
!
        DO nf=1,nIceS
          IF (iSice(nf).gt.0) THEN
            ifield=iSice(nf)
            IF (.not.got_var(ifield).and.VarOut(ifield,ng)) THEN
              IF (Master) WRITE (stdout,10) TRIM(Vname(1,ifield)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              RETURN
            END IF
          END IF
        END DO
!
        IF (.not.got_var(idUiER).and.VarOut(idUiER,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idUiER)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idViER).and.VarOut(idViER,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idViER)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
!
        DO nf=1,nIceF
          IF (iFice(nf).gt.0) THEN
            ifield=iFice(nf)
            IF (.not.got_var(ifield).and.VarOut(ifield,ng)) THEN
              IF (Master) WRITE (stdout,10) TRIM(Vname(1,ifield)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              RETURN
            END IF
          END IF
        END DO
      END IF QUERY
!
  10  FORMAT (/,' ICE_DEF_NF90 - unable to find variable: ',a,2x,       &
     &        ' in output NetCDF file: ',a)
      RETURN
      END SUBROUTINE ice_def_nf90
!
!*********************************************************************** 
     SUBROUTINE ice_wrt_nf90 (ng, model, tile,                          &
    &                         LBi, UBi, LBj, UBj,                       &
    &                         VarOut, S)
!*********************************************************************** 
!
      USE mod_netcdf
!
!  Imported variable declarations.
!
      logical, intent(in) :: VarOut(NV,Ngrids)
!
      integer, intent(in) :: ng, model, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
      TYPE(T_IO), intent(inout) :: S(Ngrids)
!
!  Local variable declarations.
!
      logical :: LwrtVar
!
      integer :: ifield, ifld
      integer :: gfactor, gtype, status
!
      real(dp) :: scale
!
      real(r8), pointer :: iceField(:,:)
      real(r8), pointer :: iceMask(:,:)
!
      real(r8), allocatable :: Ur2d(:,:)
      real(r8), allocatable :: Vr2d(:,:)
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", ice_wrt_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out ice prognostic variables into specified output NetCDF file.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Set grid type factor to write full (gfactor=1) fields or water
!  points (gfactor=-1) fields only.
!
# if defined WRITE_WATER && defined MASKING
      gfactor=-1
# else
      gfactor=1
# endif
!
!-----------------------------------------------------------------------
!  write out oce ice model state prognostic variables.
!-----------------------------------------------------------------------
!
      DO ifld=1,nIceS
        IF (iSice(ifld).gt.0) THEN
          ifield=iSice(ifld)
          IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
            LwrtVar=.TRUE.
          ELSE
            LwrtVar=VarOut(ifield,ng)
          END IF
!
          IF (LwrtVar) THEN
            IF ((model.eq.iNLM).and.                                    &
                ((S(ng)%ncid.eq.HIS(ng)%ncid).or.                       &
     &           (S(ng)%ncid.eq.QCK(ng)%ncid).or.                       &
     &           (S(ng)%ncid.eq.RST(ng)%ncid))) THEN
              iceField => ICE(ng) % Si(LBi:UBi,LBj:UBj,IUOUT,ifield)
# ifdef AVERAGES
            ELSE IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
              iceField => ICE_SAVG(ifield,ng) % var(LBi:UBi,LBj:UBj)
# endif
            END IF
!
            SELECT CASE (ifld)
              CASE (isUice)
                gtype=gfactor*u2dvar
# ifdef MASKING
                iceMask => GRID(ng) % umask_full
# endif
              CASE (isVice)
                gtype=gfactor*v2dvar
# ifdef MASKING
                iceMask => GRID(ng) % vmask_full
# endif
              CASE DEFAULT
                gtype=gfactor*r2dvar
# ifdef MASKING
                iceMask => GRID(ng) % rmask_full
# endif
            END SELECT
            scale=1.0_dp
!
            status=nf_fwrite2d(ng, model, S(ng)%ncid, ifield,           &
     &                         S(ng)%Vid(ifield),                       &
     &                         S(ng)%Rindex, gtype,                     &
     &                         LBi, UBi, LBj, UBj, scale,               &
# ifdef MASKING
     &                         iceMask,                                 &
# endif
     &                         iceField)
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) WRITE (stdout,10) TRIM(Vname(1,ifield)),      &
     &                                      S(ng)%Rindex,               &
     &                                      TRIM(S(ng)%name)
              exit_flag=3
              ioerror=status
              RETURN
            END IF
          END IF
        END IF
      END DO
!
!  Write ice EASTward and NORTHward velocity component (m/s) at
!  RHO-points.
!
      IF (VarOut(idUiER,ng).and.VarOut(idViNR,ng)) THEN
        IF (.not.allocated(Ur2d)) THEN
          allocate (Ur2d(LBi:UBi,LBj:UBj))
          Ur2d(LBi:UBi,LBj:UBj)=0.0_r8
        END IF
        IF (.not.allocated(Vr2d)) THEN
          allocate (Vr2d(LBi:UBi,LBj:UBj))
          Vr2d(LBi:UBi,LBj:UBj)=0.0_r8
        END IF
        CALL uv_rotate2d (ng, tile, .FALSE., .TRUE.,                    &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    GRID(ng) % CosAngler,                         &
     &                    GRID(ng) % SinAngler,                         &
# ifdef MASKING
     &                    GRID(ng) % rmask_full,                        &
# endif
     &                    ICE(ng) % Si(:,:,IUOUT,isUice),               &
     &                    ICE(ng) % Si(:,:,IUOUT,isVice),               &
     &                    Ur2d, Vr2d)
!
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, S(ng)%ncid, idUiER,               &
     &                     S(ng)%Vid(idUiER),                           &
     &                     S(ng)%Rindex, gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
# endif
     &                     Ur2d)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idUiER)),          &
     &                                  S(ng)%Rindex, TRIM(S(ng)%name)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
        status=nf_fwrite2d(ng, model, S(ng)%ncid, idViNR,               &
     &                     S(ng)%Vid(idViNR),                           &
     &                     S(ng)%Rindex, gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
# endif
     &                     Vr2d)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idViNR)),          &
     &                                  S(ng)%Rindex, TRIM(S(ng)%name)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Ur2d)
        deallocate (Vr2d)
      END IF
!
!-----------------------------------------------------------------------
!  write out ice model internal variables.
!-----------------------------------------------------------------------
!
      DO ifld=1,nIceS
        IF (iSice(ifld).gt.0) THEN
          ifield=iSice(ifld)
          IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
            LwrtVar=.FALSE.
          ELSE
            LwrtVar=VarOut(ifield,ng)
          END IF
!
          IF (LwrtVar) THEN
            IF ((model.eq.iNLM).and.                                    &
                ((S(ng)%ncid.eq.HIS(ng)%ncid).or.                       &
     &           (S(ng)%ncid.eq.QCK(ng)%ncid).or.                       &
     &           (S(ng)%ncid.eq.RST(ng)%ncid))) THEN
              iceField => ICE(ng) % Fi(LBi:UBi,LBj:UBj,ifield)
# ifdef AVERAGES
            ELSE IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
              iceField => ICE_FAVG(ifield,ng) % var(LBi:UBi,LBj:UBj)
# endif
            END IF
!
            scale=1.0_dp
            gtype=gfactor*r2dvar
            status=nf_fwrite2d(ng, model, S(ng)%ncid, ifield,           &
     &                         S(ng)%Vid(ifield),                       &
     &                         S(ng)%Rindex, gtype,                     &
     &                         LBi, UBi, LBj, UBj, scale,               &
# ifdef MASKING
     &                         GRID(ng) % rmask_full,                   &
# endif
     &                         iceField)
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) WRITE (stdout,10) TRIM(Vname(1,ifield)),      &
     &                                      S(ng)%Rindex,               &
     &                                      TRIM(S(ng)%name)
              exit_flag=3
              ioerror=status
              RETURN
            END IF
          END IF
        END IF
      END DO
!
  10  FORMAT (/," ICE_WRT - error while writing variable '",a,          &
     &        "', time record = ",i0,/,11x,'into file: ',a)
!
      RETURN
      END SUBROUTINE ice_wrt_nf90

# if defined PIO_LIB && defined DISTRIBUTE
!
!***********************************************************************
      SUBROUTINE ice_def_pio (ng, model, ldef, VarOut, S,               &
     &                        t2dgrd, u2dgrd, v2dgrd)
!***********************************************************************
!
      USE mod_pio_netcdf
!
!  Imported variable declarations.
!
      logical, intent(in) :: ldef, VarOut(NV,Ngrids)
!
      integer, intent(in) :: ng, model
      integer, intent(in), optional :: t2dgrd(:), u2dgrd(:), v2dgrd(:)
!
      TYPE(T_IO), intent(inout) :: S(Ngrids)
!
!  Local variable declarations.
!
      logical :: LdefVar
!
      integer, parameter :: Natt = 25

      integer :: i, ifield, j, nf, nvd3, status

      integer :: icegrd(3)
!
      real(r8) :: Aval(6)
!
      character (len=120) :: Vinfo(Natt)
      character (len=256) :: ncname
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", ice_def_pio"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Define ice model state prognostic variables.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ncname=S(ng)%name
!!
      DEFINE : IF (ldef) THEN
!
!  Set number of dimensions for output variables.
!
# if defined WRITE_WATER && defined MASKING
        nvd3=2
# else
        nvd3=3
# endif
!
!  Initialize local information variable arrays.
!
        DO i=1,Natt
          DO j=1,LEN(Vinfo(1))
            Vinfo(i)(j:j)=' '
          END DO
        END DO
        DO i=1,6
          Aval(i)=0.0_r8
        END DO
!
!  Define ice model state prognostic variables.
!
        DO nf=1,nIceS
          IF (iSice(nf).gt.0) THEN
            ifield=iSice(nf)
            IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
              LdefVar=.TRUE.
            ELSE
              LdefVar=VarOut(ifield,ng)
            END IF    
            IF (Ldefvar) THEN
              Vinfo( 1)=Vname(1,ifield)
              Vinfo( 2)=Vname(2,ifield)
              Vinfo( 3)=Vname(3,ifield)
              Vinfo(14)=Vname(4,ifield)
              Vinfo(16)=Vname(1,idtime)
              Vinfo(21)=Vname(6,ifield)
              Vinfo(22)='coordinates'
              Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
              S(ng)%pioVar(ifield)%dkind=PIO_FOUT
!
              SELECT CASE (nf)
                CASE (isUice)
                  icegrd(1:3)=u2dgrd(1:3)
# if defined WRITE_WATER && defined MASKING
                  Vinfo(20)='mask_u'
# endif
                  S(ng)%pioVar(ifield)%gtype=u2dvar
                CASE (isVice)
                  icegrd(1:3)=v2dgrd(1:3)
# if defined WRITE_WATER && defined MASKING
                  Vinfo(20)='mask_v'
# endif
                  S(ng)%pioVar(ifield)%gtype=v2dvar
                CASE DEFAULT
                  icegrd(1:3)=t2dgrd(1:3)
# if defined WRITE_WATER && defined MASKING
                  Vinfo(20)='mask_rho'
# endif
                  S(ng)%pioVar(ifield)%gtype=r2dvar
              END SELECT
!
              status=def_var(ng, model, S(ng)%pioFile,                  &
     &                       S(ng)%pioVar(idUice)%vd,                   &
     &                       PIO_FOUT, nvd3, icegrd,                    &
     &                       Aval, Vinfo, ncname)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END IF
          END IF
        END DO
!
!  Define ice EASTward velocity at RHO-points.
!
        IF (VarOut(idUiER,ng)) THEN
          Vinfo( 1)=Vname(1,idUiER)
          Vinfo( 2)=Vname(2,idUiER)
          Vinfo( 3)=Vname(3,idUiER)
          Vinfo(14)=Vname(4,idUiER)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(21)=Vname(6,idUiER)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUiER,ng),r8)
          S(ng)%pioVar(idUiER)%dkind=PIO_FOUT
          S(ng)%pioVar(idUiER)%gtype=r2dvar
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idUiER)%vd,                       &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define ice NORTHward velocity at RHO-points.
!
        IF (VarOut(idViNR,ng)) THEN
          Vinfo( 1)=Vname(1,idViNR)
          Vinfo( 2)=Vname(2,idViNR)
          Vinfo( 3)=Vname(3,idViNR)
          Vinfo(14)=Vname(4,idViNR)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(21)=Vname(6,idViNR)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idViNR,ng),r8)
          S(ng)%pioVar(idViNR)%dkind=PIO_FOUT
          S(ng)%pioVar(idViNR)%gtype=r2dvar
!
          status=def_var(ng, model, S(ng)%pioFile,                      &
     &                   S(ng)%pioVar(idViNR)%vd,                       &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!-----------------------------------------------------------------------
!  Define ice model internal variables.
!-----------------------------------------------------------------------
!
        DO nf=1,nIceF
          IF (iFice(nf).gt.0) THEN
            ifield=iFice(nf)
            IF (VarOut(ifield,ng)) THEN
              Vinfo( 1)=Vname(1,ifield)
              Vinfo( 2)=Vname(2,ifield)
              Vinfo( 3)=Vname(3,ifield)
              Vinfo(14)=Vname(4,ifield)
              Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
              Vinfo(20)='mask_rho'
# endif
              Vinfo(21)=Vname(6,ifield)
              Vinfo(22)='coordinates'
              Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
              S(ng)%pioVar(ifield)%dkind=PIO_FOUT
              S(ng)%pioVar(ifield)%gtype=r2dvar
!
              SELECT CASE (nf)
                CASE (icWdiv)
                  Vinfo(11)='increase ice thickness'
                  Vinfo(12)='decrease ice concentration'
                CASE (icW_ai)
                  Vinfo(11)='freezing'
                  Vinfo(12)='melting'
                CASE (icW_ao, icW_io)
                  Vinfo(11)='melting'
                  Vinfo(12)='freezing'
              END SELECT
!
              status=def_var(ng, model, S(ng)%pioFile,                  &
     &                       S(ng)%pioVar(idAice)%vd,                   &
     &                       PIO_FOUT, nvd3, t2dgrd,                    &
     &                       Aval, Vinfo, ncname)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END IF
          END IF
        END DO
      END IF DEFINE
!
!-----------------------------------------------------------------------
!  Otherwise, check existing output file and prepare for appending
!  data.
!-----------------------------------------------------------------------
!
      QUERY : IF (.not.ldef) THEN
!
!  Initialize locallogical switches.
!
        DO i=1,NV
          got_var(i)=.FALSE.
        END DO
!
!  Scan variable list from input NetCDF and activate switches for
!  sea-ice variables. Get variable IDs.
!
        DO i=1,n_var
          DO nf=1,nIceS
            IF (iSice(nf).gt.0) THEN
              ifield=iSice(nf)
              IF (TRIM(var_name(i)).eq.TRIM(Vname(1,ifield))) THEN
                got_var(ifield)=.TRUE.
                S(ng)%pioVar(ifield)%vd=var_desc(i)
                S(ng)%pioVar(ifield)%dkind=PIO_FOUT
!
                SELECT CASE (nf)
                  CASE (isUice)
                    S(ng)%pioVar(ifield)%gtype=u2dvar
                  CASE (isVice)
                    S(ng)%pioVar(ifield)%gtype=v2dvar
                  CASE DEFAULT
                    S(ng)%pioVar(ifield)%gtype=r2dvar
                END SELECT
!
                CYCLE
              END IF
            END IF
          END DO          
!
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUiER))) THEN
            got_var(idUiER)=.TRUE.
            S(ng)%pioVar(idUiER)%vd=var_desc(i)
            S(ng)%pioVar(idUiER)%dkind=PIO_FOUT
            S(ng)%pioVar(idUiER)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idViER))) THEN
            got_var(idViER)=.TRUE.
            S(ng)%pioVar(idViER)%vd=var_desc(i)
            S(ng)%pioVar(idViER)%dkind=PIO_FOUT
            S(ng)%pioVar(idViER)%gtype=r2dvar
          END IF
!
          DO nf=1,nIceF
            IF (iFice(nf).gt.0) THEN
              ifield=iFice(nf)
              IF (TRIM(var_name(i)).eq.TRIM(Vname(1,ifield))) THEN
                got_var(ifield)=.TRUE.
                S(ng)%pioVar(ifield)%vd=var_desc(i)
                S(ng)%pioVar(ifield)%dkind=PIO_FOUT
                S(ng)%pioVar(ifield)%gtype=r2dvar
                CYCLE
              END IF
            END IF
          END DO
        END DO
!
!  Check if output variables are available in input NetCDF file.
!
        DO nf=1,nIceS
          IF (iSice(nf).gt.0) THEN
            ifield=iSice(nf)
            IF (.not.got_var(ifield).and.VarOut(ifield,ng)) THEN
              IF (Master) WRITE (stdout,10) TRIM(Vname(1,ifield)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              RETURN
            END IF
          END IF
        END DO
!
        IF (.not.got_var(idUiER).and.VarOut(idUiER,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idUiER)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idViER).and.VarOut(idViER,ng)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idViER)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
!
        DO nf=1,nIceF
          IF (iFice(nf).gt.0) THEN
            ifield=iFice(nf)
            IF (.not.got_var(ifield).and.VarOut(ifield,ng)) THEN
              IF (Master) WRITE (stdout,10) TRIM(Vname(1,ifield)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              RETURN
            END IF
          END IF
        END DO
      END IF QUERY
!
  10  FORMAT (/,' ICE_DEF_PIO - unable to find variable: ',a,2x,        &
     &        ' in output NetCDF file: ',a)
!
      RETURN
      END SUBROUTINE ice_def_pio
!
!*********************************************************************** 
     SUBROUTINE ice_wrt_pio (ng, model, tile,                           &
    &                        LBi, UBi, LBj, UBj,                        &
    &                        VarOut, S)
!*********************************************************************** 
!
      USE mod_pio_netcdf
!
!  Imported variable declarations.
!
      logical, intent(in) :: VarOut(NV,Ngrids)
!
      integer, intent(in) :: ng, model, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
      TYPE(T_IO), intent(inout) :: S(Ngrids)
!
!  Local variable declarations.
!
      logical :: LwrtVar
!
      integer :: ifield, ifld
      integer :: gfactor, gtype, status
!
      real(dp) :: scale
!
      real(r8), pointer :: iceField(:,:)
      real(r8), pointer :: iceMask(:,:)
!
      real(r8), allocatable :: Ur2d(:,:)
      real(r8), allocatable :: Vr2d(:,:)
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", ice_wrt_pio"
!
      TYPE (IO_desc_t), pointer :: ioDesc
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out ice prognostic variables into specified output NetCDF file.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Set grid type factor to write full (gfactor=1) fields or water
!  points (gfactor=-1) fields only.
!
# if defined WRITE_WATER && defined MASKING
      gfactor=-1
# else
      gfactor=1
# endif
!
!-----------------------------------------------------------------------
!  write out oce ice model state prognostic variables.
!-----------------------------------------------------------------------
!
      DO ifld=1,nIceS
        IF (iSice(ifld).gt.0) THEN
          ifield=iSice(ifld)
          IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
            LwrtVar=.TRUE.
          ELSE
            LwrtVar=VarOut(ifield,ng)
          END IF
!
          IF (LwrtVar) THEN
            IF ((model.eq.iNLM).and.                                    &
                ((S(ng)%ncid.eq.HIS(ng)%ncid).or.                       &
     &           (S(ng)%ncid.eq.QCK(ng)%ncid).or.                       &
     &           (S(ng)%ncid.eq.RST(ng)%ncid))) THEN
              iceField => ICE(ng) % Si(LBi:UBi,LBj:UBj,IUOUT,ifield)
# ifdef AVERAGES
            ELSE IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
              iceField => ICE_SAVG(ifield,ng) % var(LBi:UBi,LBj:UBj)
# endif
            END IF
!
            SELECT CASE (ifld)
              CASE (isUice)
                IF (S(ng)%pioVar(ifield)%dkind.eq.PIO_double) THEN
                  ioDesc => ioDesc_dp_u2dvar(ng)
                ELSE
                  ioDesc => ioDesc_sp_u2dvar(ng)
                END IF
# ifdef MASKING
                iceMask => GRID(ng) % umask_full
# endif
              CASE (isVice)
                IF (S(ng)%pioVar(ifield)%dkind.eq.PIO_double) THEN
                  ioDesc => ioDesc_dp_v2dvar(ng)
                ELSE
                  ioDesc => ioDesc_sp_v2dvar(ng)
                END IF
# ifdef MASKING
                iceMask => GRID(ng) % vmask_full
# endif
              CASE DEFAULT
                IF (S(ng)%pioVar(ifield)%dkind.eq.PIO_double) THEN
                  ioDesc => ioDesc_dp_r2dvar(ng)
                ELSE
                  ioDesc => ioDesc_sp_r2dvar(ng)
                END IF
# ifdef MASKING
                iceMask => GRID(ng) % rmask_full
# endif
            END SELECT
            scale=1.0_dp
!
            status=nf_fwrite2d(ng, model, S(ng)%pioFile, ifield,        &
     &                         S(ng)%pioVar(ifield), S(ng)%Rindex,      &
     &                         ioDesc,                                  &
     &                         LBi, UBi, LBj, UBj, scale,               &
# ifdef MASKING
     &                         iceMask,                                 &
# endif
     &                         iceField)
            IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
              IF (Master) WRITE (stdout,10) TRIM(Vname(1,ifield)),      &
     &                                      S(ng)%Rindex,               &
     &                                      TRIM(S(ng)%name)
              exit_flag=3
              ioerror=status
              RETURN
            END IF
          END IF
        END IF
      END DO
!
!  Write ice EASTward and NORTHward velocity component (m/s) at
!  RHO-points.
!
      IF (VarOut(idUiER,ng).and.VarOut(idViNR,ng)) THEN
        IF (.not.allocated(Ur2d)) THEN
          allocate (Ur2d(LBi:UBi,LBj:UBj))
          Ur2d(LBi:UBi,LBj:UBj)=0.0_r8
        END IF
        IF (.not.allocated(Vr2d)) THEN
          allocate (Vr2d(LBi:UBi,LBj:UBj))
          Vr2d(LBi:UBi,LBj:UBj)=0.0_r8
        END IF
        CALL uv_rotate2d (ng, tile, .FALSE., .TRUE.,                    &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    GRID(ng) % CosAngler,                         &
     &                    GRID(ng) % SinAngler,                         &
# ifdef MASKING
     &                    GRID(ng) % rmask_full,                        &
# endif
     &                    ICE(ng) % Si(:,:,IUOUT,isUice),               &
     &                    ICE(ng) % Si(:,:,IUOUT,isVice),               &
     &                    Ur2d, Vr2d)
!
        scale=1.0_dp
        IF (S(ng)%pioVar(idUiER)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, S(ng)%pioFile, idUiER,            &
     &                     S(ng)%pioVAR(idUiER), S(ng)%Rindex,          &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
# endif
     &                     Ur2d)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idUiER)),          &
     &                                  S(ng)%Rindex, TRIM(S(ng)%name)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
        IF (S(ng)%pioVar(idViNR)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, S(ng)%pioFile, idViNR,            &
     &                     S(ng)%pioVar(idViNR), S(ng)%Rindex,          &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
# endif
     &                     Vr2d)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idViNR)),          &
     &                                  S(ng)%Rindex, TRIM(S(ng)%name)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Ur2d)
        deallocate (Vr2d)
      END IF
!
!-----------------------------------------------------------------------
!  write out ice model internal variables.
!-----------------------------------------------------------------------
!
      DO ifld=1,nIceS
        IF (iSice(ifld).gt.0) THEN
          ifield=iSice(ifld)
          IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
            LwrtVar=.FALSE.
          ELSE
            LwrtVar=VarOut(ifield,ng)
          END IF
!
          IF (LwrtVar) THEN
            IF ((model.eq.iNLM).and.                                    &
                ((S(ng)%ncid.eq.HIS(ng)%ncid).or.                       &
     &           (S(ng)%ncid.eq.QCK(ng)%ncid).or.                       &
     &           (S(ng)%ncid.eq.RST(ng)%ncid))) THEN
              iceField => ICE(ng) % Fi(LBi:UBi,LBj:UBj,ifield)
# ifdef AVERAGES
            ELSE IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
              iceField => ICE_FAVG(ifield,ng) % var(LBi:UBi,LBj:UBj)
# endif
            END IF
!
            scale=1.0_dp
            IF (S(ng)%pioVar(ifield)%dkind.eq.PIO_double) THEN
              ioDesc => ioDesc_dp_r2dvar(ng)
            ELSE
              ioDesc => ioDesc_sp_r2dvar(ng)
            END IF
            status=nf_fwrite2d(ng, model, S(ng)%pioFile, ifield,        &
     &                         S(ng)%pioVar(ifield), S(ng)%Rindex,      &
     &                         ioDesc,                                  &
     &                         LBi, UBi, LBj, UBj, scale,               &
# ifdef MASKING
     &                         GRID(ng) % rmask_full,                   &
# endif
     &                         iceField)
            IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
              IF (Master) WRITE (stdout,10) TRIM(Vname(1,ifield)),      &
     &                                      S(ng)%Rindex,               &
     &                                      TRIM(S(ng)%name)
              exit_flag=3
              ioerror=status
              RETURN
            END IF
          END IF
        END IF
      END DO
!
  10  FORMAT (/," ICE_WRT - error while writing variable '",a,          &
     &        "', time record = ",i0,/,11x,'into file: ',a)
!
      RETURN
      END SUBROUTINE ice_wrt_pio
# endif
!
#endif
      END MODULE ice_output_mod
