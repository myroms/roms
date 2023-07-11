#include "cppdefs.h"
      MODULE wrt_station_mod
#ifdef STATIONS
!
!git $Id$
!svn $Id: wrt_station.F 1178 2023-07-11 17:50:57Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module writes STATION data into output file using the standard !
!  NetCDF library or the Parallel-IO (PIO) library.                    !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
# ifdef BBL_MODEL
      USE mod_bbl
# endif
      USE mod_forces
      USE mod_grid
      USE mod_iounits
      USE mod_mixing
      USE mod_ncparam
      USE mod_ocean
      USE mod_scalars
# if defined SEDIMENT || defined BBL_MODEL
      USE mod_sedbed
      USE mod_sediment
# endif
      USE mod_stepping
!
      USE extract_sta_mod, ONLY : extract_sta2d
# ifdef SOLVE3D
      USE extract_sta_mod, ONLY : extract_sta3d
# endif
      USE uv_rotate_mod,   ONLY : uv_rotate2d
# ifdef SOLVE3D
      USE uv_rotate_mod,   ONLY : uv_rotate3d
# endif
      USE strings_mod,     ONLY : FoundError
# if defined WEC || defined WEC_VF
      USE wec_output_mod,  ONLY : wec_wrt_station_nf90
#  if defined PIO_LIB && defined DISTRIBUTE
      USE wec_output_mod,  ONLY : wec_wrt_station_pio
#  endif
# endif
!
      implicit none
!
      PUBLIC  :: wrt_station
      PRIVATE :: wrt_station_nf90
# if defined PIO_LIB && defined DISTRIBUTE
      PRIVATE :: wrt_station_pio
# endif
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE wrt_station (ng, tile)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
      integer :: LBi, UBi, LBj, UBj
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
!-----------------------------------------------------------------------
!  Write out history fields according to IO type.
!-----------------------------------------------------------------------
!
      LBi=BOUNDS(ng)%LBi(tile)
      UBi=BOUNDS(ng)%UBi(tile)
      LBj=BOUNDS(ng)%LBj(tile)
      UBj=BOUNDS(ng)%UBj(tile)
!
      SELECT CASE (STA(ng)%IOtype)
        CASE (io_nf90)
          CALL wrt_station_nf90 (ng, iNLM, tile,                        &
     &                           LBi, UBi, LBj, UBj)

# if defined PIO_LIB && defined DISTRIBUTE
        CASE (io_pio)
          CALL wrt_station_pio (ng, iNLM, tile,                         &
     &                          LBi, UBi, LBj, UBj)
# endif
        CASE DEFAULT
          IF (Master) WRITE (stdout,10) STA(ng)%IOtype
          exit_flag=3
      END SELECT
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
  10  FORMAT (' WRT_STATION - Illegal output file type, io_type = ',i0, &
     &        /,15x,'Check KeyWord ''OUT_LIB'' in ''roms.in''.')
!
      RETURN
      END SUBROUTINE wrt_station
!
!***********************************************************************
      SUBROUTINE wrt_station_nf90 (ng, model, tile,                     &
     &                             LBi, UBi, LBj, UBj)
!***********************************************************************
!
      USE mod_netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
!  Local variable declarations.
!
      logical :: Cgrid
!
      integer :: NposB, NposR, NposW
      integer :: Fcount, i, ifield, k, np, status
!
      real(dp) :: scale

      real(r8), dimension(Nstation(ng)) :: Xpos, Ypos, Zpos, psta
# ifdef SOLVE3D
#  ifdef SEDIMENT
      real(r8), dimension(Nstation(ng)*Nbed) :: XposB, YposB, ZposB
      real(r8), dimension(Nstation(ng)*Nbed) :: bsta
#  endif
      real(r8), dimension(Nstation(ng)*(N(ng))) :: XposR, YposR, ZposR
      real(r8), dimension(Nstation(ng)*(N(ng)+1)) :: XposW, YposW, ZposW
      real(r8), dimension(Nstation(ng)*(N(ng)+1)) :: rsta
# endif

      real(r8), allocatable :: Ur2d(:,:)
      real(r8), allocatable :: Vr2d(:,:)
# ifdef SOLVE3D
      real(r8), allocatable :: Ur3d(:,:,:)
      real(r8), allocatable :: Vr3d(:,:,:)
# endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", wrt_station_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out station data at RHO-points.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Set time record index.
!
      STA(ng)%Rindex=STA(ng)%Rindex+1
      Fcount=STA(ng)%Fcount
      STA(ng)%Nrec(Fcount)=STA(ng)%Nrec(Fcount)+1
!
!  Set switch to extract station data at native C-grid position (TRUE)
!  or at RHO-points (FALSE).
!
# ifdef STATIONS_CGRID
      Cgrid=.TRUE.
# else
      Cgrid=.FALSE.
# endif
!
!  Set positions for generic extraction routine.
!
      NposB=Nstation(ng)*Nbed
      NposR=Nstation(ng)*N(ng)
      NposW=Nstation(ng)*(N(ng)+1)
      DO i=1,Nstation(ng)
        Xpos(i)=SCALARS(ng)%SposX(i)
        Ypos(i)=SCALARS(ng)%SposY(i)
        Zpos(i)=1.0_r8
# ifdef SOLVE3D
        DO k=1,N(ng)
          np=k+(i-1)*N(ng)
          XposR(np)=SCALARS(ng)%SposX(i)
          YposR(np)=SCALARS(ng)%SposY(i)
          ZposR(np)=REAL(k,r8)
        END DO
        DO k=0,N(ng)
          np=k+1+(i-1)*(N(ng)+1)
          XposW(np)=SCALARS(ng)%SposX(i)
          YposW(np)=SCALARS(ng)%SposY(i)
          ZposW(np)=REAL(k,r8)
        END DO
#  ifdef SEDIMENT
        DO k=1,Nbed
          np=k+(i-1)*Nbed
          XposB(np)=SCALARS(ng)%SposX(i)
          YposB(np)=SCALARS(ng)%SposY(i)
          ZposB(np)=REAL(k,r8)
        END DO
#  endif
# endif
      END DO
!
!  Write out model time (s).
!
      CALL netcdf_put_fvar (ng, model, STA(ng)%name,                    &
     &                      TRIM(Vname(1,idtime)), time(ng:),           &
     &                      (/STA(ng)%Rindex/), (/1/),                  &
     &                      ncid = STA(ng)%ncid,                        &
     &                      varid = STA(ng)%Vid(idtime))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out free-surface (m).
!
      IF (Sout(idFsur,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idFsur, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, OCEAN(ng)%zeta(:,:,KOUT),            &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idFsur)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idFsur))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF

# if defined SEDIMENT && defined SED_MORPH
!
!  Write out time-varying bathymetry.
!
      IF (Sout(idbath,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idbath, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, GRID(ng)%h,                          &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idbath)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idbath))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
# endif
!
!  Write out 2D momentum component (m/s) in the XI-direction.
!
      IF (Sout(idUbar,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idUbar, u2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, OCEAN(ng)%ubar(:,:,KOUT),            &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idUbar)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idUbar))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out 2D momentum component (m/s) in the ETA-direction.
!
      IF (Sout(idVbar,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idVbar, v2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, OCEAN(ng)%vbar(:,:,KOUT),            &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idVbar)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idVbar))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out 2D Eastward and Northward momentum components (m/s) at
!  RHO-points
!
      IF (Sout(idu2dE,ng).and.Sout(idv2dN,ng)) THEN
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
     &                    OCEAN(ng) % ubar(:,:,KOUT),                   &
     &                    OCEAN(ng) % vbar(:,:,KOUT),                   &
     &                    Ur2d, Vr2d)
!
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idu2dE, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, Ur2d,                                &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idu2dE)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idu2dE))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
        CALL extract_sta2d (ng, model, Cgrid, idv2dN, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, Vr2d,                                &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idv2dN)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idv2dN))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        deallocate (Ur2d)
        deallocate (Vr2d)
      END IF

# ifdef SOLVE3D
!
!  Write out 3D momentum component (m/s) in the XI-direction.
!
      IF (Sout(idUvel,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, model, Cgrid, idUvel, u3dvar,           &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      scale, OCEAN(ng)%u(:,:,:,NOUT),             &
     &                      NposR, XposR, YposR, ZposR, rsta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idUvel)), rsta,              &
     &                        (/1,1,STA(ng)%Rindex/),                   &
     &                        (/N(ng),Nstation(ng),1/),                 &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idUvel))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out 3D momentum component (m/s) in the ETA-direction.
!
      IF (Sout(idVvel,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, model, Cgrid, idVvel, v3dvar,           &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      scale, OCEAN(ng)%v(:,:,:,NOUT),             &
     &                      NposR, XposR, YposR, ZposR, rsta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idVvel)), rsta,              &
     &                        (/1,1,STA(ng)%Rindex/),                   &
     &                        (/N(ng),Nstation(ng),1/),                 &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idVvel))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out 3D Eastward and Northward momentum components (m/s) at
!  RHO-points.
!
      IF (Sout(idu3dE,ng).and.Sout(idv3dN,ng)) THEN
        IF (.not.allocated(Ur3d)) THEN
          allocate (Ur3d(LBi:UBi,LBj:UBj,N(ng)))
          Ur3d(LBi:UBi,LBj:UBj,1:N(ng))=0.0_r8
        END IF
        IF (.not.allocated(Vr3d)) THEN
          allocate (Vr3d(LBi:UBi,LBj:UBj,N(ng)))
          Vr3d(LBi:UBi,LBj:UBj,1:N(ng))=0.0_r8
        END IF
        CALL uv_rotate3d (ng, tile, .FALSE., .TRUE.,                    &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    GRID(ng) % CosAngler,                         &
     &                    GRID(ng) % SinAngler,                         &
#  ifdef MASKING
     &                    GRID(ng) % rmask_full,                        &
#  endif
     &                    OCEAN(ng) % u(:,:,:,NOUT),                    &
     &                    OCEAN(ng) % v(:,:,:,NOUT),                    &
     &                    Ur3d, Vr3d)
!
        scale=1.0_dp
        CALL extract_sta3d (ng, model, Cgrid, idu3dE, r3dvar,           &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      scale, Ur3d,                                &
     &                      NposR, XposR, YposR, ZposR, rsta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idu3dE)), rsta,              &
     &                        (/1,1,STA(ng)%Rindex/),                   &
     &                        (/N(ng),Nstation(ng),1/),                 &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idu3dE))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
        CALL extract_sta3d (ng, model, Cgrid, idv3dN, r3dvar,           &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      scale, Vr3d,                                &
     &                      NposR, XposR, YposR, ZposR, rsta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idv3dN)), rsta,              &
     &                        (/1,1,STA(ng)%Rindex/),                   &
     &                        (/N(ng),Nstation(ng),1/),                 &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idv3dN))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        deallocate (Ur3d)
        deallocate (Vr3d)
      END IF
!
!  Write out "true" vertical velocity (m/s).
!
      IF (Sout(idWvel,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, model, Cgrid, idWvel, w3dvar,           &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      scale, OCEAN(ng)%wvel,                      &
     &                      NposW, XposW, YposW, ZposW, rsta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idWvel)), rsta,              &
     &                        (/1,1,STA(ng)%Rindex/),                   &
     &                        (/N(ng)+1,Nstation(ng),1/),               &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idWvel))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write S-coordinate "omega" vertical velocity (m3/s).
!
      IF (Sout(idOvel,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, model, Cgrid, idOvel, w3dvar,           &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      scale, OCEAN(ng)%W,                         &
     &                      NposW, XposW, YposW, ZposW, rsta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idOvel)), rsta,              &
     &                        (/1,1,STA(ng)%Rindex/),                   &
     &                        (/N(ng)+1,Nstation(ng),1/),               &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idOvel))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out tracer type variables.
!
      DO i=1,NT(ng)
        ifield=idTvar(i)
        IF (Sout(ifield,ng)) THEN
          scale=1.0_dp
          CALL extract_sta3d (ng, model, Cgrid, ifield, r3dvar,         &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        scale, OCEAN(ng)%t(:,:,:,NOUT,i),         &
     &                        NposR, XposR, YposR, ZposR, rsta)
          CALL netcdf_put_fvar (ng, model, STA(ng)%name,                &
     &                          TRIM(Vname(1,idTvar(i))), rsta,         &
     &                          (/1,1,STA(ng)%Rindex/),                 &
     &                          (/N(ng),Nstation(ng),1/),               &
     &                          ncid = STA(ng)%ncid,                    &
     &                          varid = STA(ng)%Tid(i))
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
      END DO
!
!  Write out density anomaly.
!
      IF (Sout(idDano,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, model, Cgrid, idDano, r3dvar,           &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      scale, OCEAN(ng)%rho,                       &
     &                      NposR, XposR, YposR, ZposR, rsta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idDano)), rsta,              &
     &                        (/1,1,STA(ng)%Rindex/),                   &
     &                        (/N(ng),Nstation(ng),1/),                 &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idDano))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF

#  ifdef LMD_SKPP
!
!  Write out depth of surface boundary layer.
!
      IF (Sout(idHsbl,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idHsbl, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, MIXING(ng)%hsbl,                     &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idHsbl)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idHsbl))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#  endif
#  ifdef LMD_BKPP
!
!  Write out depth of bottom boundary layer.
!
      IF (Sout(idHbbl,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idHbbl, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, MIXING(ng)%hbbl,                     &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idHbbl)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idHbbl))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#  endif
!
!  Write out vertical viscosity coefficient.
!
      IF (Sout(idVvis,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, model, Cgrid, idVvis, w3dvar,           &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      scale, MIXING(ng)%Akv,                      &
     &                      NposW, XposW, YposW, ZposW, rsta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idVvis)), rsta,              &
     &                        (/1,1,STA(ng)%Rindex/),                   &
     &                        (/N(ng)+1,Nstation(ng),1/),               &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idVvis))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out vertical diffusion coefficient for potential temperature.
!
      IF (Sout(idTdif,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, model, Cgrid, idTdif, w3dvar,           &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      scale, MIXING(ng)%Akt(:,:,:,itemp),         &
     &                      NposW, XposW, YposW, ZposW, rsta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idTdif)), rsta,              &
     &                        (/1,1,STA(ng)%Rindex/),                   &
     &                        (/N(ng)+1,Nstation(ng),1/),               &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idTdif))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF

#  ifdef SALINITY
!
!  Write out vertical diffusion coefficient for salinity.
!
      IF (Sout(idSdif,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, model, Cgrid, idSdif, w3dvar,           &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      scale, MIXING(ng)%Akt(:,:,:,isalt),         &
     &                      NposW, XposW, YposW, ZposW, rsta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idSdif)), rsta,              &
     &                        (/1,1,STA(ng)%Rindex/),                   &
     &                        (/N(ng)+1,Nstation(ng),1/),               &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idSdif))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#  endif
#  if defined GLS_MIXING || defined MY25_MIXING
!
!  Write out turbulent kinetic energy.
!
      IF (Sout(idMtke,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, model, Cgrid, idMtke, w3dvar,           &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      scale, MIXING(ng)%tke(:,:,:,NOUT),          &
     &                      NposW, XposW, YposW, ZposW, rsta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idMtke)), rsta,              &
     &                        (/1,1,STA(ng)%Rindex/),                   &
     &                        (/N(ng)+1,Nstation(ng),1/),               &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idMtke))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out turbulent kinetic energy times length scale.
!
      IF (Sout(idMtls,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, model, Cgrid, idMtls, w3dvar,           &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      scale, MIXING(ng)%gls(:,:,:,NOUT),          &
     &                      NposW, XposW, YposW, ZposW, rsta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idMtls)), rsta,              &
     &                        (/1,1,STA(ng)%Rindex/),                   &
     &                        (/N(ng)+1,Nstation(ng),1/),               &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idMtls))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#  endif
#  if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
!
!  Write out surface air pressure.
!
      IF (Sout(idPair,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idPair, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%Pair,                     &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idPair)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idPair))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#  endif
#  if defined BULK_FLUXES || defined ECOSIM
!
!  Write out surface winds.
!
      IF (Sout(idUair,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idUair, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%Uwind,                    &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idUair)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idUair))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
      IF (Sout(idVair,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idVair, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%Vwind,                    &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idVair)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idVair))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out 2D Eastward and Northward surface winds (m/s) at
!  RHO-points
!
      IF (Sout(idUaiE,ng).and.Sout(idVaiN,ng)) THEN
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
     &                    FORCES(ng) % Uwind,                           &
     &                    FORCES(ng) % Vwind,                           &
     &                    Ur2d, Vr2d)
!
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idUaiE, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, Ur2d,                                &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idUaiE)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idUairE))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
        CALL extract_sta2d (ng, model, Cgrid, idVaiN, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, Vr2d,                                &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idVairN)), psta,             &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idVaiN))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        deallocate (Ur2d)
        deallocate (Vr2d)
      END IF
#  endif
!
!  Write out surface net heat flux.
!
      IF (Sout(idTsur(itemp),ng)) THEN
        ifield=idTsur(itemp)
        scale=rho0*Cp
        CALL extract_sta2d (ng, model, Cgrid, ifield, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%stflx(:,:,itemp),         &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,ifield)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(ifield))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF

#  ifdef SALINITY
!
!  Write out surface salt flux.
!
      IF (Sout(idTsur(isalt),ng)) THEN
        ifield=idTsur(isalt)
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, ifield, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%stflx(:,:,isalt),         &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,ifield)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(ifield))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#  endif

#  ifdef BULK_FLUXES
!
!  Write out latent heat flux.
!
      IF (Sout(idLhea,ng)) THEN
        scale=rho0*Cp
        CALL extract_sta2d (ng, model, Cgrid, idLhea, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%lhflx,                    &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idLhea)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idLhea))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out sensible heat flux.
!
      IF (Sout(idShea,ng)) THEN
        scale=rho0*Cp
        CALL extract_sta2d (ng, model, Cgrid, idShea, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%shflx,                    &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idShea)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idShea))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out longwave radiation flux.
!
      IF (Sout(idLrad,ng)) THEN
        scale=rho0*Cp
        CALL extract_sta2d (ng, model, Cgrid, idLrad, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%lrflx,                    &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idLrad)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idLrad))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF

#  endif
#  ifdef SHORTWAVE
!
!  Write out shortwave radiation flux.
!
      IF (Sout(idSrad,ng)) THEN
        scale=rho0*Cp
        CALL extract_sta2d (ng, model, Cgrid, idSrad, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%srflx,                    &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idSrad)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idSrad))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#  endif
#  if defined EMINUSP && defined BULK_FLUXES
!
!  Write out E-P (m/s).
!
      IF (Sout(idEmPf,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idEmPf, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%stflux(:,:,isalt),        &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idEmPf)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idEmPf))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out evaporation rate (kg/m2/s).
!
      IF (Sout(idevap,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idevap, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%evap,                     &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idevap)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idevap))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out precipitation rate (kg/m2/s).
!
      IF (Sout(idrain,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idrain, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%rain,                     &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idrain)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idrain))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#  endif
# endif
!
!  Write out surface U-momentum stress.
!
      IF (Sout(idUsms,ng)) THEN
        scale=rho0
        CALL extract_sta2d (ng, model, Cgrid, idUsms, u2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%sustr,                    &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idUsms)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idUsms))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out surface V-momentum stress.
!
      IF (Sout(idVsms,ng)) THEN
        scale=rho0
        CALL extract_sta2d (ng, model, Cgrid, idVsms, v2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%svstr,                    &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idVsms)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idVsms))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out bottom U-momentum stress.
!
      IF (Sout(idUbms,ng)) THEN
        scale=-rho0
        CALL extract_sta2d (ng, model, Cgrid, idUbms, u2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%bustr,                    &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idUbms)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idUbms))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out bottom V-momentum stress.
!
      IF (Sout(idVbms,ng)) THEN
        scale=-rho0
        CALL extract_sta2d (ng, model, Cgrid, idVbms, v2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%bvstr,                    &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idVbms)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idVbms))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF

# ifdef SOLVE3D
#  ifdef BBL_MODEL
!
!  Write out current-induced, bottom U-stress.
!
      IF (Sout(idUbrs,ng)) THEN
        scale=-rho0
        CALL extract_sta2d (ng, model, Cgrid, idUbrs, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, BBL(ng)%bustrc,                      &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idUbrs)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idUbrs))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out current-induced, bottom V-stress.
!
      IF (Sout(idVbrs,ng)) THEN
        scale=-rho0
        CALL extract_sta2d (ng, model, Cgrid, idVbrs, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, BBL(ng)%bvstrc,                      &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idVbrs)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idVbrs))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out wind-induced, bottom U-stress.
!
      IF (Sout(idUbws,ng)) THEN
        scale=rho0
        CALL extract_sta2d (ng, model, Cgrid, idUbws, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, BBL(ng)%bustrw,                      &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idUbws)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idUbws))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out wind-induced, bottom V-wave stress.
!
      IF (Sout(idVbws,ng)) THEN
        scale=rho0
        CALL extract_sta2d (ng, model, Cgrid, idVbws, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, BBL(ng)%bvstrw,                      &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idVbws)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idVbws))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out maximum wind and current, bottom U-stress.
!
      IF (Sout(idUbcs,ng)) THEN
        scale=rho0
        CALL extract_sta2d (ng, model, Cgrid, idUbcs, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, BBL(ng)%bustrcwmax,                  &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idUbcs)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idUbcs))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out maximum wind and current, bottom V-stress.
!
      IF (Sout(idVbcs,ng)) THEN
        scale=rho0
        CALL extract_sta2d (ng, model, Cgrid, idVbcs, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, BBL(ng)%bvstrcwmax,                  &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idVbcs)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idVbcs))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out wind-induced, bed wave orbital U-velocity.
!
      IF (Sout(idUbot,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idUbot, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, BBL(ng)%Ubot,                        &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idUbot)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idUbot))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out wind-induced, bed wave orbital V-velocity.
!
      IF (Sout(idVbot,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idVbot, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, BBL(ng)%Vbot,                        &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idVbot)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idVbot))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out bottom U-velocity above bed.
!
      IF (Sout(idUbur,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idUbur, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, BBL(ng)%Ur,                          &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idUbur)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idUbur))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out bottom V-velocity above bed.
!
      IF (Sout(idVbvr,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idVbvr, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, BBL(ng)%Vr,                          &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idVbvr)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idVbvr))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#  endif

#  ifdef SEDIMENT
!
!  Write out sediment fraction of each size class in each bed layer.
!
      DO i=1,NST
        IF (Sout(idfrac(i),ng)) THEN
          scale=1.0_dp
          CALL extract_sta3d (ng, model, Cgrid, idfrac(i), b3dvar,      &
     &                        LBi, UBi, LBj, UBj, 1, Nbed,              &
     &                        scale, SEDBED(ng)%bed_frac(:,:,:,i),      &
     &                        NposB, XposB, YposB, ZposB, bsta)
          CALL netcdf_put_fvar (ng, model, STA(ng)%name,                &
     &                          TRIM(Vname(1,idfrac(i))), rsta,         &
     &                          (/1,1,STA(ng)%Rindex/),                 &
     &                          (/Nbed,Nstation(ng),1/),                &
     &                          ncid = STA(ng)%ncid,                    &
     &                          varid = STA(ng)%Vid(idfrac(i)))
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Write out sediment mass of each size class in each bed layer.
!
        IF (Sout(idBmas(i),ng)) THEN
          scale=1.0_dp
          CALL extract_sta3d (ng, model, Cgrid, idBmas(i), b3dvar,      &
     &                        LBi, UBi, LBj, UBj, 1, Nbed,              &
     &                        scale,                                    &
     &                        SEDBED(ng)%bed_mass(:,:,:,NOUT,i),        &
     &                        NposB, XposB, YposB, ZposB, bsta)
          CALL netcdf_put_fvar (ng, model, STA(ng)%name,                &
     &                          TRIM(Vname(1,idBmas(i))), rsta,         &
     &                          (/1,1,STA(ng)%Rindex/),                 &
     &                          (/Nbed,Nstation(ng),1/),                &
     &                          ncid = STA(ng)%ncid,                    &
     &                          varid = STA(ng)%Vid(idBmas(i)))
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
      END DO
!
!  Write out sediment properties in each bed layer.
!
      DO i=1,MBEDP
        IF (Sout(idSbed(i),ng)) THEN
          scale=1.0_dp
          CALL extract_sta3d (ng, model, Cgrid, idSbed(i), b3dvar,      &
     &                        LBi, UBi, LBj, UBj, 1, Nbed,              &
     &                        scale, SEDBED(ng)%bed(:,:,:,i),           &
     &                        NposB, XposB, YposB, ZposB, bsta)
          CALL netcdf_put_fvar (ng, model, STA(ng)%name,                &
     &                          TRIM(Vname(1,idSbed(i))), rsta,         &
     &                          (/1,1,STA(ng)%Rindex/),                 &
     &                          (/Nbed,Nstation(ng),1/),                &
     &                          ncid = STA(ng)%ncid,                    &
     &                          varid = STA(ng)%Vid(idSbed(i)))
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
      END DO
#  endif

#  if defined SEDIMENT || defined BBL_MODEL
!
!  Write out exposed sediment layer properties.
!
      DO i=1,MBEDP
        IF (Sout(idBott(i),ng)) THEN
          scale=1.0_dp
          CALL extract_sta2d (ng, model, Cgrid, idBott(i), r2dvar,      &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        scale, SEDBED(ng)%bottom(:,:,i),          &
     &                        Nstation(ng), Xpos, Ypos, psta)
          CALL netcdf_put_fvar (ng, model, STA(ng)%name,                &
     &                          TRIM(Vname(1,idBott(i))), rsta,         &
     &                          (/1,STA(ng)%Rindex/),                   &
     &                          (/Nstation(ng),1/),                     &
     &                          ncid = STA(ng)%ncid,                    &
     &                          varid = STA(ng)%Vid(idBott(i)))
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
      END DO
#  endif
# endif
# ifdef WAVES_HEIGHT
!
!  Write out wind-induced significant wave height.
!
      IF (Sout(idWamp,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idWamp, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Hwave,                  &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idWamp)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idWamp))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
# endif
# ifdef WAVES_LENGTH
!
!  Write out wind-induced mean wavelenght.
!
      IF (Sout(idWlen,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idWlen, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Lwave,                  &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idWlen)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idWlen))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
# endif
# ifdef WAVES_LENGTHP
!
!  Write out wind-induced peak wave wavelength.
!
      IF (Sout(idWlep,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idWlep, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Lwavep,                 &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idWlep)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idWlep))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
# endif
# ifdef WAVES_DIR
!
!  Write out wind-induced mean wave direction.
!
      IF (Sout(idWdir,ng)) THEN
        scale=rad2deg
        CALL extract_sta2d (ng, model, Cgrid, idWdir, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Dwave,                  &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idWdir)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idWdir))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
# endif
# ifdef WAVES_DIRP
!
!  Write out wind-induced peak wave direction.
!
      IF (Sout(idWdip,ng)) THEN
        scale=rad2deg
        CALL extract_sta2d (ng, model, Cgrid, idWdip, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Dwavep,                 &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idWdip)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idWdip))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
# endif
# ifdef WAVES_TOP_PERIOD
!
!  Write out wind-induced surface wave period.
!
      IF (Sout(idWptp,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idWptp, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Pwave_top,              &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idWptp)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idWptp))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
# endif
# ifdef WAVES_BOT_PERIOD
!
!  Write out wind-induced bottom wave period.
!
      IF (Sout(idWpbt,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idWpbt, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Pwave_bot,              &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idWpbt)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idWpbt))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
# endif
#if defined BBL_MODEL
!
!  Write out wind-induced wave bottom orbital velocity.
!
      IF (Sout(idWorb,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idWorb, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Uwave_rms,              &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idWorb)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idWorb))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
# endif
# if defined WAVES_DSPR
!
!  Write out waves directional spreading.
!
      IF (Sout(idWvds,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idWvds, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Wave_ds,                &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idWvds)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idWvds))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out wave spectrum peakedness.
!
      IF (Sout(idWvqp,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idWvqp, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Wave_qp,                &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idWvqp)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idWvqp))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
# endif
# ifdef WET_DRY
!
!  Write out wet/dry mask at RHO-points.
!
      IF (Sout(idRwet,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idRwet, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, GRID(ng)%rmask_wet,                  &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idRwet)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idRwet))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out wet/dry mask at U-points.
!
      IF (Sout(idUwet,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idUwet, u2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, GRID(ng)%umask_wet,                  &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idUwet)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idUwet))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out wet/dry mask at V-points.
!
      IF (Sout(idVwet,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idVwet, v2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, GRID(ng)%vmask_wet,                  &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idVwet)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idVwet))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
# endif
# if defined WEC || defined WEC_VF
!
!-----------------------------------------------------------------------
!  Write out Waves Effect on Currents station variables.
!-----------------------------------------------------------------------
!
      CALL wec_wrt_station_nf90 (ng, model, tile,                       &
    &                            LBi, UBi, LBj, UBj,                    &
    &                            Sout, STA)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
!
!-----------------------------------------------------------------------
!  Synchronize stations NetCDF file to disk.
!-----------------------------------------------------------------------
!
      CALL netcdf_sync (ng, model, STA(ng)%name, STA(ng)%ncid)

      RETURN
      END SUBROUTINE wrt_station_nf90

# if defined PIO_LIB && defined DISTRIBUTE
!
!***********************************************************************
      SUBROUTINE wrt_station_pio (ng, model, tile,                      &
     &                            LBi, UBi, LBj, UBj)
!***********************************************************************
!
      USE mod_pio_netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
!  Local variable declarations.
!
      logical :: Cgrid
!
      integer :: NposB, NposR, NposW
      integer :: Fcount, i, ifield, k, np, status
!
      real(dp) :: scale

      real(r8), dimension(Nstation(ng)) :: Xpos, Ypos, Zpos, psta
#  ifdef SOLVE3D
#   ifdef SEDIMENT
      real(r8), dimension(Nstation(ng)*Nbed) :: XposB, YposB, ZposB
      real(r8), dimension(Nstation(ng)*Nbed) :: bsta
#   endif
      real(r8), dimension(Nstation(ng)*(N(ng))) :: XposR, YposR, ZposR
      real(r8), dimension(Nstation(ng)*(N(ng)+1)) :: XposW, YposW, ZposW
      real(r8), dimension(Nstation(ng)*(N(ng)+1)) :: rsta
#  endif

      real(r8), allocatable :: Ur2d(:,:)
      real(r8), allocatable :: Vr2d(:,:)
#  ifdef SOLVE3D
      real(r8), allocatable :: Ur3d(:,:,:)
      real(r8), allocatable :: Vr3d(:,:,:)
#  endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", wrt_station_pio"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out station data at RHO-points.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Set time record index.
!
      STA(ng)%Rindex=STA(ng)%Rindex+1
      Fcount=STA(ng)%Fcount
      STA(ng)%Nrec(Fcount)=STA(ng)%Nrec(Fcount)+1
!
!  Set switch to extract station data at native C-grid position (TRUE)
!  or at RHO-points (FALSE).
!
#  ifdef STATIONS_CGRID
      Cgrid=.TRUE.
#  else
      Cgrid=.FALSE.
#  endif
!
!  Set positions for generic extraction routine.
!
      NposB=Nstation(ng)*Nbed
      NposR=Nstation(ng)*N(ng)
      NposW=Nstation(ng)*(N(ng)+1)
      DO i=1,Nstation(ng)
        Xpos(i)=SCALARS(ng)%SposX(i)
        Ypos(i)=SCALARS(ng)%SposY(i)
        Zpos(i)=1.0_r8
#  ifdef SOLVE3D
        DO k=1,N(ng)
          np=k+(i-1)*N(ng)
          XposR(np)=SCALARS(ng)%SposX(i)
          YposR(np)=SCALARS(ng)%SposY(i)
          ZposR(np)=REAL(k,r8)
        END DO
        DO k=0,N(ng)
          np=k+1+(i-1)*(N(ng)+1)
          XposW(np)=SCALARS(ng)%SposX(i)
          YposW(np)=SCALARS(ng)%SposY(i)
          ZposW(np)=REAL(k,r8)
        END DO
#   ifdef SEDIMENT
        DO k=1,Nbed
          np=k+(i-1)*Nbed
          XposB(np)=SCALARS(ng)%SposX(i)
          YposB(np)=SCALARS(ng)%SposY(i)
          ZposB(np)=REAL(k,r8)
        END DO
#   endif
#  endif
      END DO
!
!  Write out model time (s).
!
      CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,                &
     &                          TRIM(Vname(1,idtime)), time(ng:),       &
     &                          (/STA(ng)%Rindex/), (/1/),              &
     &                          pioFile = STA(ng)%pioFile,              &
     &                          pioVar = STA(ng)%pioVar(idtime)%vd)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out free-surface (m).
!
      IF (Sout(idFsur,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idFsur, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, OCEAN(ng)%zeta(:,:,KOUT),            &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idFsur)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idFsur)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF

#  if defined SEDIMENT && defined SED_MORPH
!
!  Write out time-varying bathymetry.
!
      IF (Sout(idbath,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idbath, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, GRID(ng)%h,                          &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idbath)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idbath)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#  endif
!
!  Write out 2D momentum component (m/s) in the XI-direction.
!
      IF (Sout(idUbar,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idUbar, u2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, OCEAN(ng)%ubar(:,:,KOUT),            &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idUbar)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idUbar)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out 2D momentum component (m/s) in the ETA-direction.
!
      IF (Sout(idVbar,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idVbar, v2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, OCEAN(ng)%vbar(:,:,KOUT),            &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idVbar)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idVbar)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out 2D Eastward and Northward momentum components (m/s) at
!  RHO-points
!
      IF (Sout(idu2dE,ng).and.Sout(idv2dN,ng)) THEN
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
#  ifdef MASKING
     &                    GRID(ng) % rmask_full,                        &
#  endif
     &                    OCEAN(ng) % ubar(:,:,KOUT),                   &
     &                    OCEAN(ng) % vbar(:,:,KOUT),                   &
     &                    Ur2d, Vr2d)
!
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idu2dE, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, Ur2d,                                &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idu2dE)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idu2dE)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
        CALL extract_sta2d (ng, model, Cgrid, idv2dN, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, Vr2d,                                &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idv2dN)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idv2dN)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        deallocate (Ur2d)
        deallocate (Vr2d)
      END IF

#  ifdef SOLVE3D
!
!  Write out 3D momentum component (m/s) in the XI-direction.
!
      IF (Sout(idUvel,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, model, Cgrid, idUvel, u3dvar,           &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      scale, OCEAN(ng)%u(:,:,:,NOUT),             &
     &                      NposR, XposR, YposR, ZposR, rsta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idUvel)), rsta,          &
     &                            (/1,1,STA(ng)%Rindex/),               &
     &                            (/N(ng),Nstation(ng),1/),             &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idUvel)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out 3D momentum component (m/s) in the ETA-direction.
!
      IF (Sout(idVvel,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, model, Cgrid, idVvel, v3dvar,           &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      scale, OCEAN(ng)%v(:,:,:,NOUT),             &
     &                      NposR, XposR, YposR, ZposR, rsta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idVvel)), rsta,          &
     &                            (/1,1,STA(ng)%Rindex/),               &
     &                            (/N(ng),Nstation(ng),1/),             &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idVvel)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out 3D Eastward and Northward momentum components (m/s) at
!  RHO-points.
!
      IF (Sout(idu3dE,ng).and.Sout(idv3dN,ng)) THEN
        IF (.not.allocated(Ur3d)) THEN
          allocate (Ur3d(LBi:UBi,LBj:UBj,N(ng)))
          Ur3d(LBi:UBi,LBj:UBj,1:N(ng))=0.0_r8
        END IF
        IF (.not.allocated(Vr3d)) THEN
          allocate (Vr3d(LBi:UBi,LBj:UBj,N(ng)))
          Vr3d(LBi:UBi,LBj:UBj,1:N(ng))=0.0_r8
        END IF
        CALL uv_rotate3d (ng, tile, .FALSE., .TRUE.,                    &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    GRID(ng) % CosAngler,                         &
     &                    GRID(ng) % SinAngler,                         &
#   ifdef MASKING
     &                    GRID(ng) % rmask_full,                        &
#   endif
     &                    OCEAN(ng) % u(:,:,:,NOUT),                    &
     &                    OCEAN(ng) % v(:,:,:,NOUT),                    &
     &                    Ur3d, Vr3d)
!
        scale=1.0_dp
        CALL extract_sta3d (ng, model, Cgrid, idu3dE, r3dvar,           &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      scale, Ur3d,                                &
     &                      NposR, XposR, YposR, ZposR, rsta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idu3dE)), rsta,          &
     &                            (/1,1,STA(ng)%Rindex/),               &
     &                            (/N(ng),Nstation(ng),1/),             &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idu3dE)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
        CALL extract_sta3d (ng, model, Cgrid, idv3dN, r3dvar,           &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      scale, Vr3d,                                &
     &                      NposR, XposR, YposR, ZposR, rsta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idv3dN)), rsta,          &
     &                            (/1,1,STA(ng)%Rindex/),               &
     &                            (/N(ng),Nstation(ng),1/),             &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idv3dN)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        deallocate (Ur3d)
        deallocate (Vr3d)
      END IF
!
!  Write out "true" vertical velocity (m/s).
!
      IF (Sout(idWvel,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, model, Cgrid, idWvel, w3dvar,           &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      scale, OCEAN(ng)%wvel,                      &
     &                      NposW, XposW, YposW, ZposW, rsta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idWvel)), rsta,          &
     &                            (/1,1,STA(ng)%Rindex/),               &
     &                            (/N(ng)+1,Nstation(ng),1/),           &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idWvel)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write S-coordinate "omega" vertical velocity (m3/s).
!
      IF (Sout(idOvel,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, model, Cgrid, idOvel, w3dvar,           &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      scale, OCEAN(ng)%W,                         &
     &                      NposW, XposW, YposW, ZposW, rsta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idOvel)), rsta,          &
     &                            (/1,1,STA(ng)%Rindex/),               &
     &                            (/N(ng)+1,Nstation(ng),1/),           &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idOvel)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out tracer type variables.
!
      DO i=1,NT(ng)
        ifield=idTvar(i)
        IF (Sout(ifield,ng)) THEN
          scale=1.0_dp
          CALL extract_sta3d (ng, model, Cgrid, ifield, r3dvar,         &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        scale, OCEAN(ng)%t(:,:,:,NOUT,i),         &
     &                        NposR, XposR, YposR, ZposR, rsta)
          CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,            &
     &                              TRIM(Vname(1,idTvar(i))), rsta,     &
     &                              (/1,1,STA(ng)%Rindex/),             &
     &                              (/N(ng),Nstation(ng),1/),           &
     &                              pioFile = STA(ng)%pioFile,          &
     &                              pioVar = STA(ng)%pioTrc(i)%vd)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
      END DO
!
!  Write out density anomaly.
!
      IF (Sout(idDano,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, model, Cgrid, idDano, r3dvar,           &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      scale, OCEAN(ng)%rho,                       &
     &                      NposR, XposR, YposR, ZposR, rsta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idDano)), rsta,          &
     &                            (/1,1,STA(ng)%Rindex/),               &
     &                            (/N(ng),Nstation(ng),1/),             &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idDano)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF

#   ifdef LMD_SKPP
!
!  Write out depth of surface boundary layer.
!
      IF (Sout(idHsbl,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idHsbl, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, MIXING(ng)%hsbl,                     &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idHsbl)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idHsbl)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#   endif
#   ifdef LMD_BKPP
!
!  Write out depth of bottom boundary layer.
!
      IF (Sout(idHbbl,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idHbbl, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, MIXING(ng)%hbbl,                     &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idHbbl)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idHbbl)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#   endif
!
!  Write out vertical viscosity coefficient.
!
      IF (Sout(idVvis,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, model, Cgrid, idVvis, w3dvar,           &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      scale, MIXING(ng)%Akv,                      &
     &                      NposW, XposW, YposW, ZposW, rsta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idVvis)), rsta,          &
     &                            (/1,1,STA(ng)%Rindex/),               &
     &                            (/N(ng)+1,Nstation(ng),1/),           &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idVvis)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out vertical diffusion coefficient for potential temperature.
!
      IF (Sout(idTdif,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, model, Cgrid, idTdif, w3dvar,           &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      scale, MIXING(ng)%Akt(:,:,:,itemp),         &
     &                      NposW, XposW, YposW, ZposW, rsta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idTdif)), rsta,          &
     &                            (/1,1,STA(ng)%Rindex/),               &
     &                            (/N(ng)+1,Nstation(ng),1/),           &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idTdif)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF

#   ifdef SALINITY
!
!  Write out vertical diffusion coefficient for salinity.
!
      IF (Sout(idSdif,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, model, Cgrid, idSdif, w3dvar,           &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      scale, MIXING(ng)%Akt(:,:,:,isalt),         &
     &                      NposW, XposW, YposW, ZposW, rsta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idSdif)), rsta,          &
     &                            (/1,1,STA(ng)%Rindex/),               &
     &                            (/N(ng)+1,Nstation(ng),1/),           &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idSdif)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#   endif
#   if defined GLS_MIXING || defined MY25_MIXING
!
!  Write out turbulent kinetic energy.
!
      IF (Sout(idMtke,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, model, Cgrid, idMtke, w3dvar,           &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      scale, MIXING(ng)%tke(:,:,:,NOUT),          &
     &                      NposW, XposW, YposW, ZposW, rsta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idMtke)), rsta,          &
     &                            (/1,1,STA(ng)%Rindex/),               &
     &                            (/N(ng)+1,Nstation(ng),1/),           &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idMtke)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out turbulent kinetic energy times length scale.
!
      IF (Sout(idMtls,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, model, Cgrid, idMtls, w3dvar,           &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      scale, MIXING(ng)%gls(:,:,:,NOUT),          &
     &                      NposW, XposW, YposW, ZposW, rsta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idMtls)), rsta,          &
     &                            (/1,1,STA(ng)%Rindex/),               &
     &                            (/N(ng)+1,Nstation(ng),1/),           &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idMtls)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#   endif
#   if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
!
!  Write out surface air pressure.
!
      IF (Sout(idPair,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idPair, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%Pair,                     &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idPair)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idPair)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#   endif
#   if defined BULK_FLUXES || defined ECOSIM
!
!  Write out surface winds.
!
      IF (Sout(idUair,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idUair, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%Uwind,                    &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idUair)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idUair)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
      IF (Sout(idVair,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idVair, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%Vwind,                    &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idVair)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idVair)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out 2D Eastward and Northward surface winds (m/s) at
!  RHO-points
!
      IF (Sout(idUaiE,ng).and.Sout(idVaiN,ng)) THEN
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
     &                    FORCES(ng) % Uwind,                           &
     &                    FORCES(ng) % Vwind,                           &
     &                    Ur2d, Vr2d)
!
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idUaiE, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, Ur2d,                                &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idUaiE)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idUaiE)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
        CALL extract_sta2d (ng, model, Cgrid, idVaiN, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, Vr2d,                                &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idVaiN)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idVaiN)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        deallocate (Ur2d)
        deallocate (Vr2d)
      END IF
#   endif
!
!  Write out surface net heat flux.
!
      IF (Sout(idTsur(itemp),ng)) THEN
        ifield=idTsur(itemp)
        scale=rho0*Cp
        CALL extract_sta2d (ng, model, Cgrid, ifield, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%stflx(:,:,itemp),         &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,ifield)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(ifield)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF

#   ifdef SALINITY
!
!  Write out surface salt flux.
!
      IF (Sout(idTsur(isalt),ng)) THEN
        ifield=idTsur(isalt)
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, ifield, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%stflx(:,:,isalt),         &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,ifield)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(ifield)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#   endif

#   ifdef BULK_FLUXES
!
!  Write out latent heat flux.
!
      IF (Sout(idLhea,ng)) THEN
        scale=rho0*Cp
        CALL extract_sta2d (ng, model, Cgrid, idLhea, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%lhflx,                    &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idLhea)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idLhea)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out sensible heat flux.
!
      IF (Sout(idShea,ng)) THEN
        scale=rho0*Cp
        CALL extract_sta2d (ng, model, Cgrid, idShea, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%shflx,                    &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idShea)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idShea)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out longwave radiation flux.
!
      IF (Sout(idLrad,ng)) THEN
        scale=rho0*Cp
        CALL extract_sta2d (ng, model, Cgrid, idLrad, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%lrflx,                    &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idLrad)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idLrad)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF

#   endif
#   ifdef SHORTWAVE
!
!  Write out shortwave radiation flux.
!
      IF (Sout(idSrad,ng)) THEN
        scale=rho0*Cp
        CALL extract_sta2d (ng, model, Cgrid, idSrad, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%srflx,                    &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idSrad)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idSrad)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#   endif
#   if defined EMINUSP && defined BULK_FLUXES
!
!  Write out E-P (m/s).
!
      IF (Sout(idEmPf,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idEmPf, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%stflux(:,:,isalt),        &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idEmPf)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idEmPf)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out evaportaion rate (kg/m2/s).
!
      IF (Sout(idevap,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idevap, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%evap,                     &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idevap)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idevap)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out precipitation rate (kg/m2/s).
!
      IF (Sout(idrain,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idrain, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%rain,                     &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idrain)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idrain)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#   endif
#  endif
!
!  Write out surface U-momentum stress.
!
      IF (Sout(idUsms,ng)) THEN
        scale=rho0
        CALL extract_sta2d (ng, model, Cgrid, idUsms, u2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%sustr,                    &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idUsms)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idUsms)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out surface V-momentum stress.
!
      IF (Sout(idVsms,ng)) THEN
        scale=rho0
        CALL extract_sta2d (ng, model, Cgrid, idVsms, v2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%svstr,                    &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idVsms)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idVsms)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out bottom U-momentum stress.
!
      IF (Sout(idUbms,ng)) THEN
        scale=-rho0
        CALL extract_sta2d (ng, model, Cgrid, idUbms, u2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%bustr,                    &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idUbms)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idUbms)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out bottom V-momentum stress.
!
      IF (Sout(idVbms,ng)) THEN
        scale=-rho0
        CALL extract_sta2d (ng, model, Cgrid, idVbms, v2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%bvstr,                    &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idVbms)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idVbms)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF

#  ifdef SOLVE3D
#   ifdef BBL_MODEL
!
!  Write out current-induced, bottom U-stress.
!
      IF (Sout(idUbrs,ng)) THEN
        scale=-rho0
        CALL extract_sta2d (ng, model, Cgrid, idUbrs, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, BBL(ng)%bustrc,                      &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idUbrs)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idUbrs)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out current-induced, bottom V-stress.
!
      IF (Sout(idVbrs,ng)) THEN
        scale=-rho0
        CALL extract_sta2d (ng, model, Cgrid, idVbrs, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, BBL(ng)%bvstrc,                      &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idVbrs)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idVbrs)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out wind-induced, bottom U-stress.
!
      IF (Sout(idUbws,ng)) THEN
        scale=rho0
        CALL extract_sta2d (ng, model, Cgrid, idUbws, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, BBL(ng)%bustrw,                      &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idUbws)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idUbws)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out wind-induced, bottom V-wave stress.
!
      IF (Sout(idVbws,ng)) THEN
        scale=rho0
        CALL extract_sta2d (ng, model, Cgrid, idVbws, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, BBL(ng)%bvstrw,                      &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idVbws)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idVbws)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out maximum wind and current, bottom U-stress.
!
      IF (Sout(idUbcs,ng)) THEN
        scale=rho0
        CALL extract_sta2d (ng, model, Cgrid, idUbcs, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, BBL(ng)%bustrcwmax,                  &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idUbcs)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idUbcs)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out maximum wind and current, bottom V-stress.
!
      IF (Sout(idVbcs,ng)) THEN
        scale=rho0
        CALL extract_sta2d (ng, model, Cgrid, idVbcs, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, BBL(ng)%bvstrcwmax,                  &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idVbcs)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idVbcs)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out wind-induced, bed wave orbital U-velocity.
!
      IF (Sout(idUbot,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idUbot, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, BBL(ng)%Ubot,                        &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idUbot)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idUbot)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out wind-induced, bed wave orbital V-velocity.
!
      IF (Sout(idVbot,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idVbot, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, BBL(ng)%Vbot,                        &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idVbot)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idVbot)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out bottom U-velocity above bed.
!
      IF (Sout(idUbur,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idUbur, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, BBL(ng)%Ur,                          &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idUbur)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idUbur)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out bottom V-velocity above bed.
!
      IF (Sout(idVbvr,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idVbvr, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, BBL(ng)%Vr,                          &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idVbvr)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idVbvr)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#   endif

#   ifdef SEDIMENT
!
!  Write out sediment fraction of each size class in each bed layer.
!
      DO i=1,NST
        IF (Sout(idfrac(i),ng)) THEN
          scale=1.0_dp
          CALL extract_sta3d (ng, model, Cgrid, idfrac(i), b3dvar,      &
     &                        LBi, UBi, LBj, UBj, 1, Nbed,              &
     &                        scale, SEDBED(ng)%bed_frac(:,:,:,i),      &
     &                        NposB, XposB, YposB, ZposB, bsta)
          CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,            &
     &                              TRIM(Vname(1,idfrac(i))), rsta,     &
     &                              (/1,1,STA(ng)%Rindex/),             &
     &                              (/Nbed,Nstation(ng),1/),            &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idfrac(i))%vd)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Write out sediment mass of each size class in each bed layer.
!
        IF (Sout(idBmas(i),ng)) THEN
          scale=1.0_dp
          CALL extract_sta3d (ng, model, Cgrid, idBmas(i), b3dvar,      &
     &                        LBi, UBi, LBj, UBj, 1, Nbed,              &
     &                        scale,                                    &
     &                        SEDBED(ng)%bed_mass(:,:,:,NOUT,i),        &
     &                        NposB, XposB, YposB, ZposB, bsta)
          CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,            &
     &                              TRIM(Vname(1,idBmas(i))), rsta,     &
     &                              (/1,1,STA(ng)%Rindex/),             &
     &                              (/Nbed,Nstation(ng),1/),            &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idBmas(i))%vd)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
      END DO
!
!  Write out sediment properties in each bed layer.
!
      DO i=1,MBEDP
        IF (Sout(idSbed(i),ng)) THEN
          scale=1.0_dp
          CALL extract_sta3d (ng, model, Cgrid, idSbed(i), b3dvar,      &
     &                        LBi, UBi, LBj, UBj, 1, Nbed,              &
     &                        scale, SEDBED(ng)%bed(:,:,:,i),           &
     &                        NposB, XposB, YposB, ZposB, bsta)
          CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,            &
     &                              TRIM(Vname(1,idSbed(i))), rsta,     &
     &                              (/1,1,STA(ng)%Rindex/),             &
     &                              (/Nbed,Nstation(ng),1/),            &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idSbed(i))%vd)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
      END DO
#   endif

#   if defined SEDIMENT || defined BBL_MODEL
!
!  Write out exposed sediment layer properties.
!
      DO i=1,MBEDP
        IF (Sout(idBott(i),ng)) THEN
          scale=1.0_dp
          CALL extract_sta2d (ng, model, Cgrid, idBott(i), r2dvar,      &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        scale, SEDBED(ng)%bottom(:,:,i),          &
     &                        Nstation(ng), Xpos, Ypos, psta)
          CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,            &
     &                              TRIM(Vname(1,idBott(i))), rsta,     &
     &                              (/1,STA(ng)%Rindex/),               &
     &                              (/Nstation(ng),1/),                 &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idBott(i))%vd)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
      END DO
#   endif
#  endif
#  ifdef WAVES_HEIGHT
!
!  Write out wind-induced significant wave height.
!
      IF (Sout(idWamp,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idWamp, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Hwave,                  &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idWamp)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idWamp)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#  endif
#  ifdef WAVES_LENGTH
!
!  Write out wind-induced mean wavelenght.
!
      IF (Sout(idWlen,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idWlen, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Lwave,                  &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idWlen)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idWlen)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#  endif
#  ifdef WAVES_LENGTHP
!
!  Write out wind-induced peak wave wavelength.
!
      IF (Sout(idWlep,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idWlep, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Lwavep,                 &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idWlep)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idWlep)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#  endif
#  ifdef WAVES_DIR
!
!  Write out wind-induced mean wave direction.
!
      IF (Sout(idWdir,ng)) THEN
        scale=rad2deg
        CALL extract_sta2d (ng, model, Cgrid, idWdir, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Dwave,                  &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idWdir)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idWdir)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#  endif
#  ifdef WAVES_DIRP
!
!  Write out wind-induced peak wave direction.
!
      IF (Sout(idWdip,ng)) THEN
        scale=rad2deg
        CALL extract_sta2d (ng, model, Cgrid, idWdip, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Dwavep,                 &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idWdip)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idWdip)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#  endif
#  ifdef WAVES_TOP_PERIOD
!
!  Write out wind-induced surface wave period.
!
      IF (Sout(idWptp,ng)) THEN
        scale=rad2deg
        CALL extract_sta2d (ng, model, Cgrid, idWptp, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Pwave_top,              &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idWptp)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idWptp)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#  endif
#  ifdef WAVES_BOT_PERIOD
!
!  Write out wind-induced bottom wave period.
!
      IF (Sout(idWpbt,ng)) THEN
        scale=rad2deg
        CALL extract_sta2d (ng, model, Cgrid, idWpbt, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Pwave_bot,              &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idWpbt)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idWpbt)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#  endif
#  ifdef BBL_MODEL
!
!  Write out wind-induced wave bottom orbital velocity.
!
      IF (Sout(idWorb,ng)) THEN
        scale=rad2deg
        CALL extract_sta2d (ng, model, Cgrid, idWorb, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Ub_swan,                &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idWorb)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idWorb)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#  endif
#  if defined WAVES_DSPR
!
!  Write out waves directional spreading.
!
      IF (Sout(idWvds,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idWvds, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Wave_ds,                &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idWvds)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idWvds)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out wave spectrum peakedness.
!
      IF (Sout(idWvqp,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idWvqp, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Wave_qp,                &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idWvqp)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idWvqp)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#  endif
#  ifdef WET_DRY
!
!  Write out wet/dry mask at RHO-points.
!
      IF (Sout(idRwet,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idRwet, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, GRID(ng)%rmask_wet,                  &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idRwet)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idRwet)%vd)
      END IF
!
!  Write out wet/dry mask at U-points.
!
      IF (Sout(idUwet,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idUwet, u2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, GRID(ng)%umask_wet,                  &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idUwet)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idUwet)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out wet/dry mask at V-points.
!
      IF (Sout(idVwet,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idVwet, v2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, GRID(ng)%vmask_wet,                  &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idVwet)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idVwet)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#  endif
#  if defined WEC || defined WEC_VF
!
!-----------------------------------------------------------------------
!  Write out Waves Effect on Currents station variables.
!-----------------------------------------------------------------------
!
      CALL wec_wrt_station_pio (ng, model, tile,                        &
    &                           LBi, UBi, LBj, UBj,                     &
    &                           Sout, STA)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
!
!-----------------------------------------------------------------------
!  Synchronize stations NetCDF file to disk.
!-----------------------------------------------------------------------
!
      CALL pio_netcdf_sync (ng, model, STA(ng)%name, STA(ng)%pioFile)

      RETURN
      END SUBROUTINE wrt_station_pio
# endif
#endif
      END MODULE wrt_station_mod

