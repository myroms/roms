#include "cppdefs.h"
#ifdef STATIONS
      SUBROUTINE wrt_station (ng)
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine writes out data into stations NetCDF file.             !
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
      USE mod_netcdf
      USE mod_ocean
      USE mod_scalars
# if defined SEDIMENT || defined BBL_MODEL
      USE mod_sediment
# endif
      USE mod_stepping
!
# ifdef SOLVE3D
      USE extract_sta_mod, ONLY : extract_sta2d, extract_sta3d
# else
      USE extract_sta_mod, ONLY : extract_sta2d
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
      logical :: Cgrid

      integer :: NposB, NposR, NposW, LBi, UBi, LBj, UBj
      integer :: i, ifield, k, np, status, tile

      real(r8) :: scale

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
!
      SourceFile='wrt_station.F'
!
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)
!
!-----------------------------------------------------------------------
!  Write out station data at RHO-points.
!-----------------------------------------------------------------------
!
      IF (exit_flag.ne.NoError) RETURN
!
!  Set time record index.
!
      tSTAindx(ng)=tSTAindx(ng)+1
      NrecSTA(ng)=NrecSTA(ng)+1
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
      CALL netcdf_put_fvar (ng, iNLM, STAname(ng),                      &
     &                      TRIM(Vname(1,idtime)), time(ng:),           &
     &                      (/tSTAindx(ng)/), (/1/),                    &
     &                      ncid = ncSTAid(ng),                         &
     &                      varid = staVid(idtime,ng))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out free-surface (m).
!
      IF (Sout(idFsur,ng)) THEN
        scale=1.0_r8
        CALL extract_sta2d (ng, iNLM, Cgrid, idFsur, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, OCEAN(ng)%zeta(:,:,KOUT),            &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STAname(ng),                    &
     &                        TRIM(Vname(1,idFsur)), psta,              &
     &                        (/1,tSTAindx(ng)/), (/Nstation(ng),1/),   &
     &                        ncid = ncSTAid(ng),                       &
     &                        varid = staVid(idFsur,ng))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out 2D momentum component (m/s) in the XI-direction.
!
      IF (Sout(idUbar,ng)) THEN
        scale=1.0_r8
        CALL extract_sta2d (ng, iNLM, Cgrid, idUbar, u2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, OCEAN(ng)%ubar(:,:,KOUT),            &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STAname(ng),                    &
     &                        TRIM(Vname(1,idUbar)), psta,              &
     &                        (/1,tSTAindx(ng)/), (/Nstation(ng),1/),   &
     &                        ncid = ncSTAid(ng),                       &
     &                        varid = staVid(idUbar,ng))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out 2D momentum component (m/s) in the ETA-direction.
!
      IF (Sout(idVbar,ng)) THEN
        scale=1.0_r8
        CALL extract_sta2d (ng, iNLM, Cgrid, idVbar, v2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, OCEAN(ng)%vbar(:,:,KOUT),            &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STAname(ng),                    &
     &                        TRIM(Vname(1,idVbar)), psta,              &
     &                        (/1,tSTAindx(ng)/), (/Nstation(ng),1/),   &
     &                        ncid = ncSTAid(ng),                       &
     &                        varid = staVid(idVbar,ng))
        IF (exit_flag.ne.NoError) RETURN
      END IF

# ifdef SOLVE3D
!
!  Write out 3D momentum component (m/s) in the XI-direction.
!
      IF (Sout(idUvel,ng)) THEN
        scale=1.0_r8
        CALL extract_sta3d (ng, iNLM, Cgrid, idUvel, u3dvar,            &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      scale, OCEAN(ng)%u(:,:,:,NOUT),             &
     &                      NposR, XposR, YposR, ZposR, rsta)
        CALL netcdf_put_fvar (ng, iNLM, STAname(ng),                    &
     &                        TRIM(Vname(1,idUvel)), rsta,              &
     &                        (/1,1,tSTAindx(ng)/),                     &
     &                        (/N(ng),Nstation(ng),1/),                 &
     &                        ncid = ncSTAid(ng),                       &
     &                        varid = staVid(idUvel,ng))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out 3D momentum component (m/s) in the ETA-direction.
!
      IF (Sout(idVvel,ng)) THEN
        scale=1.0_r8
        CALL extract_sta3d (ng, iNLM, Cgrid, idVvel, v3dvar,            &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      scale, OCEAN(ng)%v(:,:,:,NOUT),             &
     &                      NposR, XposR, YposR, ZposR, rsta)
        CALL netcdf_put_fvar (ng, iNLM, STAname(ng),                    &
     &                        TRIM(Vname(1,idVvel)), rsta,              &
     &                        (/1,1,tSTAindx(ng)/),                     &
     &                        (/N(ng),Nstation(ng),1/),                 &
     &                        ncid = ncSTAid(ng),                       &
     &                        varid = staVid(idVvel,ng))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out vertical velocity (m/s).
!
      IF (Sout(idWvel,ng)) THEN
        scale=1.0_r8
        CALL extract_sta3d (ng, iNLM, Cgrid, idWvel, w3dvar,            &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      scale, OCEAN(ng)%wvel,                      &
     &                      NposW, XposW, YposW, ZposW, rsta)
        CALL netcdf_put_fvar (ng, iNLM, STAname(ng),                    &
     &                        TRIM(Vname(1,idWvel)), rsta,              &
     &                        (/1,1,tSTAindx(ng)/),                     &
     &                        (/N(ng)+1,Nstation(ng),1/),               &
     &                        ncid = ncSTAid(ng),                       &
     &                        varid = staVid(idWvel,ng))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write S-coordinate "omega" vertical velocity (m3/s).
!
      IF (Sout(idOvel,ng)) THEN
        scale=1.0_r8
        CALL extract_sta3d (ng, iNLM, Cgrid, idOvel, w3dvar,            &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      scale, OCEAN(ng)%W,                         &
     &                      NposW, XposW, YposW, ZposW, rsta)
        CALL netcdf_put_fvar (ng, iNLM, STAname(ng),                    &
     &                        TRIM(Vname(1,idOvel)), rsta,              &
     &                        (/1,1,tSTAindx(ng)/),                     &
     &                        (/N(ng)+1,Nstation(ng),1/),               &
     &                        ncid = ncSTAid(ng),                       &
     &                        varid = staVid(idOvel,ng))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out tracer type variables.
!
      DO i=1,NT(ng)
        ifield=idTvar(i)
        IF (Sout(ifield,ng)) THEN
          scale=1.0_r8
          CALL extract_sta3d (ng, iNLM, Cgrid, ifield, r3dvar,          &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        scale, OCEAN(ng)%t(:,:,:,NOUT,i),         &
     &                        NposR, XposR, YposR, ZposR, rsta)
          CALL netcdf_put_fvar (ng, iNLM, STAname(ng),                  &
     &                          TRIM(Vname(1,idTvar(i))), rsta,         &
     &                          (/1,1,tSTAindx(ng)/),                   &
     &                          (/N(ng),Nstation(ng),1/),               &
     &                          ncid = ncSTAid(ng),                     &
     &                          varid = staTid(i,ng))
          IF (exit_flag.ne.NoError) RETURN
        END IF
      END DO
!
!  Write out density anomaly.
!
      IF (Sout(idDano,ng)) THEN
        scale=1.0_r8
        CALL extract_sta3d (ng, iNLM, Cgrid, idDano, r3dvar,            &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      scale, OCEAN(ng)%rho,                       &
     &                      NposR, XposR, YposR, ZposR, rsta)
        CALL netcdf_put_fvar (ng, iNLM, STAname(ng),                    &
     &                        TRIM(Vname(1,idDano)), rsta,              &
     &                        (/1,1,tSTAindx(ng)/),                     &
     &                        (/N(ng),Nstation(ng),1/),                 &
     &                        ncid = ncSTAid(ng),                       &
     &                        varid = staVid(idDano,ng))
        IF (exit_flag.ne.NoError) RETURN
      END IF

#  ifdef LMD_SKPP
!
!  Write out depth of surface boundary layer.
!
      IF (Sout(idHsbl,ng)) THEN
        scale=1.0_r8
        CALL extract_sta2d (ng, iNLM, Cgrid, idHsbl, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, MIXING(ng)%hsbl,                     &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STAname(ng),                    &
     &                        TRIM(Vname(1,idHsbl)), psta,              &
     &                        (/1,tSTAindx(ng)/), (/Nstation(ng),1/),   &
     &                        ncid = ncSTAid(ng),                       &
     &                        varid = staVid(idHsbl,ng))
        IF (exit_flag.ne.NoError) RETURN
      END IF
#  endif
#  ifdef LMD_BKPP
!
!  Write out depth of bottom boundary layer.
!
      IF (Sout(idHbbl,ng)) THEN
        scale=1.0_r8
        CALL extract_sta2d (ng, iNLM, Cgrid, idHbbl, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, MIXING(ng)%hbbl,                     &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STAname(ng),                    &
     &                        TRIM(Vname(1,idHbbl)), psta,              &
     &                        (/1,tSTAindx(ng)/), (/Nstation(ng),1/),   &
     &                        ncid = ncSTAid(ng),                       &
     &                        varid = staVid(idHbbl,ng))
        IF (exit_flag.ne.NoError) RETURN
      END IF
#  endif
!
!  Write out vertical viscosity coefficient.
!
      IF (Sout(idVvis,ng)) THEN
        scale=1.0_r8
        CALL extract_sta3d (ng, iNLM, Cgrid, idVvis, w3dvar,            &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      scale, MIXING(ng)%Akv,                      &
     &                      NposW, XposW, YposW, ZposW, rsta)
        CALL netcdf_put_fvar (ng, iNLM, STAname(ng),                    &
     &                        TRIM(Vname(1,idVvis)), rsta,              &
     &                        (/1,1,tSTAindx(ng)/),                     &
     &                        (/N(ng)+1,Nstation(ng),1/),               &
     &                        ncid = ncSTAid(ng),                       &
     &                        varid = staVid(idVvis,ng))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out vertical diffusion coefficient for potential temperature.
!
      IF (Sout(idTdif,ng)) THEN
        scale=1.0_r8
        CALL extract_sta3d (ng, iNLM, Cgrid, idTdif, w3dvar,            &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      scale, MIXING(ng)%Akt(:,:,:,itemp),         &
     &                      NposW, XposW, YposW, ZposW, rsta)
        CALL netcdf_put_fvar (ng, iNLM, STAname(ng),                    &
     &                        TRIM(Vname(1,idTdif)), rsta,              &
     &                        (/1,1,tSTAindx(ng)/),                     &
     &                        (/N(ng)+1,Nstation(ng),1/),               &
     &                        ncid = ncSTAid(ng),                       &
     &                        varid = staVid(idTdif,ng))
        IF (exit_flag.ne.NoError) RETURN
      END IF

#  ifdef SALINITY
!
!  Write out vertical diffusion coefficient for salinity.
!
      IF (Sout(idSdif,ng)) THEN
        scale=1.0_r8
        CALL extract_sta3d (ng, iNLM, Cgrid, idSdif, w3dvar,            &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      scale, MIXING(ng)%Akt(:,:,:,isalt),         &
     &                      NposW, XposW, YposW, ZposW, rsta)
        CALL netcdf_put_fvar (ng, iNLM, STAname(ng),                    &
     &                        TRIM(Vname(1,idSdif)), rsta,              &
     &                        (/1,1,tSTAindx(ng)/),                     &
     &                        (/N(ng)+1,Nstation(ng),1/),               &
     &                        ncid = ncSTAid(ng),                       &
     &                        varid = staVid(idSdif,ng))
        IF (exit_flag.ne.NoError) RETURN
      END IF
#  endif
#  if defined MY25_MIXING || defined GLS_MIXING
!
!  Write out turbulent kinetic energy.
!
      IF (Sout(idMtke,ng)) THEN
        scale=1.0_r8
        CALL extract_sta3d (ng, iNLM, Cgrid, idMtke, w3dvar,            &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      scale, MIXING(ng)%tke(:,:,:,NOUT),          &
     &                      NposW, XposW, YposW, ZposW, rsta)
        CALL netcdf_put_fvar (ng, iNLM, STAname(ng),                    &
     &                        TRIM(Vname(1,idMtke)), rsta,              &
     &                        (/1,1,tSTAindx(ng)/),                     &
     &                        (/N(ng)+1,Nstation(ng),1/),               &
     &                        ncid = ncSTAid(ng),                       &
     &                        varid = staVid(idMtke,ng))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out turbulent kinetic energy times length scale.
!
      IF (Sout(idMtls,ng)) THEN
        scale=1.0_r8
        CALL extract_sta3d (ng, iNLM, Cgrid, idMtls, w3dvar,            &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      scale, MIXING(ng)%gls(:,:,:,NOUT),          &
     &                      NposW, XposW, YposW, ZposW, rsta)
        CALL netcdf_put_fvar (ng, iNLM, STAname(ng),                    &
     &                        TRIM(Vname(1,idMtls)), rsta,              &
     &                        (/1,1,tSTAindx(ng)/),                     &
     &                        (/N(ng)+1,Nstation(ng),1/),               &
     &                        ncid = ncSTAid(ng),                       &
     &                        varid = staVid(idMtls,ng))
        IF (exit_flag.ne.NoError) RETURN
      END IF
#  endif
!
!  Write out surface net heat flux.
!
      IF (Sout(idTsur(itemp),ng)) THEN
        ifield=idTsur(itemp)
        scale=rho0*Cp
        CALL extract_sta2d (ng, iNLM, Cgrid, ifield, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%stflx(:,:,itemp),         &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STAname(ng),                    &
     &                        TRIM(Vname(1,ifield)), psta,              &
     &                        (/1,tSTAindx(ng)/), (/Nstation(ng),1/),   &
     &                        ncid = ncSTAid(ng),                       &
     &                        varid = staVid(ifield,ng))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out surface salt flux.
!
      IF (Sout(idTsur(isalt),ng)) THEN
        ifield=idTsur(isalt)
        scale=1.0_r8
        CALL extract_sta2d (ng, iNLM, Cgrid, ifield, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%stflx(:,:,isalt),         &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STAname(ng),                    &
     &                        TRIM(Vname(1,ifield)), psta,              &
     &                        (/1,tSTAindx(ng)/), (/Nstation(ng),1/),   &
     &                        ncid = ncSTAid(ng),                       &
     &                        varid = staVid(ifield,ng))
        IF (exit_flag.ne.NoError) RETURN
      END IF

#  ifdef BULK_FLUXES
!
!  Write out latent heat flux.
!
      IF (Sout(idLhea,ng)) THEN
        scale=rho0*Cp
        CALL extract_sta2d (ng, iNLM, Cgrid, idLhea, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%lhflx,                    &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STAname(ng),                    &
     &                        TRIM(Vname(1,idLhea)), psta,              &
     &                        (/1,tSTAindx(ng)/), (/Nstation(ng),1/),   &
     &                        ncid = ncSTAid(ng),                       &
     &                        varid = staVid(idLhea,ng))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out sensible heat flux.
!
      IF (Sout(idShea,ng)) THEN
        scale=rho0*Cp
        CALL extract_sta2d (ng, iNLM, Cgrid, idShea, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%shflx,                    &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STAname(ng),                    &
     &                        TRIM(Vname(1,idShea)), psta,              &
     &                        (/1,tSTAindx(ng)/), (/Nstation(ng),1/),   &
     &                        ncid = ncSTAid(ng),                       &
     &                        varid = staVid(idShea,ng))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out longwave radiation flux.
!
      IF (Sout(idLrad,ng)) THEN
        scale=rho0*Cp
        CALL extract_sta2d (ng, iNLM, Cgrid, idLrad, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%lrflx,                    &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STAname(ng),                    &
     &                        TRIM(Vname(1,idLrad)), psta,              &
     &                        (/1,tSTAindx(ng)/), (/Nstation(ng),1/),   &
     &                        ncid = ncSTAid(ng),                       &
     &                        varid = staVid(idLrad,ng))
        IF (exit_flag.ne.NoError) RETURN
      END IF

#  endif
#  ifdef SHORTWAVE
!
!  Write out shortwave radiation flux.
!
      IF (Sout(idSrad,ng)) THEN
        scale=rho0*Cp
        CALL extract_sta2d (ng, iNLM, Cgrid, idSrad, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%srflx,                    &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STAname(ng),                    &
     &                        TRIM(Vname(1,idSrad)), psta,              &
     &                        (/1,tSTAindx(ng)/), (/Nstation(ng),1/),   &
     &                        ncid = ncSTAid(ng),                       &
     &                        varid = staVid(idSrad,ng))
        IF (exit_flag.ne.NoError) RETURN
      END IF
#  endif
#  if defined EMINUSP && defined BULK_FLUXES
!
!  Write out E-P (m/s).
!
      IF (Sout(idEmPf,ng)) THEN
        scale=1.0_r8
        CALL extract_sta2d (ng, iNLM, Cgrid, idEmPf, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%EminusP,                  &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STAname(ng),                    &
     &                        TRIM(Vname(1,idEmPf)), psta,              &
     &                        (/1,tSTAindx(ng)/), (/Nstation(ng),1/),   &
     &                        ncid = ncSTAid(ng),                       &
     &                        varid = staVid(idEmPf,ng))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out evaportaion rate (kg/m2/s).
!
      IF (Sout(idevap,ng)) THEN
        scale=1.0_r8
        CALL extract_sta2d (ng, iNLM, Cgrid, idevap, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%evap,                     &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STAname(ng),                    &
     &                        TRIM(Vname(1,idevap)), psta,              &
     &                        (/1,tSTAindx(ng)/), (/Nstation(ng),1/),   &
     &                        ncid = ncSTAid(ng),                       &
     &                        varid = staVid(idevap,ng))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out precipitation rate (kg/m2/s).
!
      IF (Sout(idrain,ng)) THEN
        scale=1.0_r8
        CALL extract_sta2d (ng, iNLM, Cgrid, idrain, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%rain,                     &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STAname(ng),                    &
     &                        TRIM(Vname(1,idrain)), psta,              &
     &                        (/1,tSTAindx(ng)/), (/Nstation(ng),1/),   &
     &                        ncid = ncSTAid(ng),                       &
     &                        varid = staVid(idrain,ng))
        IF (exit_flag.ne.NoError) RETURN
      END IF
#  endif
# endif
!
!  Write out surface U-momentum stress.
!
      IF (Sout(idUsms,ng)) THEN
        scale=rho0
        CALL extract_sta2d (ng, iNLM, Cgrid, idUsms, u2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%sustr,                    &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STAname(ng),                    &
     &                        TRIM(Vname(1,idUsms)), psta,              &
     &                        (/1,tSTAindx(ng)/), (/Nstation(ng),1/),   &
     &                        ncid = ncSTAid(ng),                       &
     &                        varid = staVid(idUsms,ng))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out surface V-momentum stress.
!
      IF (Sout(idVsms,ng)) THEN
        scale=rho0
        CALL extract_sta2d (ng, iNLM, Cgrid, idVsms, v2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%svstr,                    &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STAname(ng),                    &
     &                        TRIM(Vname(1,idVsms)), psta,              &
     &                        (/1,tSTAindx(ng)/), (/Nstation(ng),1/),   &
     &                        ncid = ncSTAid(ng),                       &
     &                        varid = staVid(idVsms,ng))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out bottom U-momentum stress.
!
      IF (Sout(idUbms,ng)) THEN
        scale=-rho0
        CALL extract_sta2d (ng, iNLM, Cgrid, idUbms, u2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%bustr,                    &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STAname(ng),                    &
     &                        TRIM(Vname(1,idUbms)), psta,              &
     &                        (/1,tSTAindx(ng)/), (/Nstation(ng),1/),   &
     &                        ncid = ncSTAid(ng),                       &
     &                        varid = staVid(idUbms,ng))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out bottom V-momentum stress.
!
      IF (Sout(idVbms,ng)) THEN
        scale=-rho0
        CALL extract_sta2d (ng, iNLM, Cgrid, idVbms, v2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%bvstr,                    &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STAname(ng),                    &
     &                        TRIM(Vname(1,idVbms)), psta,              &
     &                        (/1,tSTAindx(ng)/), (/Nstation(ng),1/),   &
     &                        ncid = ncSTAid(ng),                       &
     &                        varid = staVid(idVbms,ng))
        IF (exit_flag.ne.NoError) RETURN
      END IF

#ifdef SOLVE3D
# ifdef BBL_MODEL
!
!  Write out current-induced, bottom U-stress.
!
      IF (Sout(idUbrs,ng)) THEN
        scale=-rho0
        CALL extract_sta2d (ng, iNLM, Cgrid, idUbrs, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, BBL(ng)%bustrc,                      &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STAname(ng),                    &
     &                        TRIM(Vname(1,idUbrs)), psta,              &
     &                        (/1,tSTAindx(ng)/), (/Nstation(ng),1/),   &
     &                        ncid = ncSTAid(ng),                       &
     &                        varid = staVid(idUbrs,ng))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out current-induced, bottom V-stress.
!
      IF (Sout(idVbrs,ng)) THEN
        scale=-rho0
        CALL extract_sta2d (ng, iNLM, Cgrid, idVbrs, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, BBL(ng)%bvstrc,                      &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STAname(ng),                    &
     &                        TRIM(Vname(1,idVbrs)), psta,              &
     &                        (/1,tSTAindx(ng)/), (/Nstation(ng),1/),   &
     &                        ncid = ncSTAid(ng),                       &
     &                        varid = staVid(idVbrs,ng))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out wind-induced, bottom U-stress.
!
      IF (Sout(idUbws,ng)) THEN
        scale=rho0
        CALL extract_sta2d (ng, iNLM, Cgrid, idUbws, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, BBL(ng)%bustrw,                      &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STAname(ng),                    &
     &                        TRIM(Vname(1,idUbws)), psta,              &
     &                        (/1,tSTAindx(ng)/), (/Nstation(ng),1/),   &
     &                        ncid = ncSTAid(ng),                       &
     &                        varid = staVid(idUbws,ng))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out wind-induced, bottom V-wave stress.
!
      IF (Sout(idVbws,ng)) THEN
        scale=rho0
        CALL extract_sta2d (ng, iNLM, Cgrid, idVbws, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, BBL(ng)%bvstrw,                      &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STAname(ng),                    &
     &                        TRIM(Vname(1,idVbws)), psta,              &
     &                        (/1,tSTAindx(ng)/), (/Nstation(ng),1/),   &
     &                        ncid = ncSTAid(ng),                       &
     &                        varid = staVid(idVbws,ng))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out maximum wind and current, bottom U-stress.
!
      IF (Sout(idUbcs,ng)) THEN
        scale=rho0
        CALL extract_sta2d (ng, iNLM, Cgrid, idUbcs, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, BBL(ng)%bustrcwmax,                  &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STAname(ng),                    &
     &                        TRIM(Vname(1,idUbcs)), psta,              &
     &                        (/1,tSTAindx(ng)/), (/Nstation(ng),1/),   &
     &                        ncid = ncSTAid(ng),                       &
     &                        varid = staVid(idUbcs,ng))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out maximum wind and current, bottom V-stress.
!
      IF (Sout(idVbcs,ng)) THEN
        scale=rho0
        CALL extract_sta2d (ng, iNLM, Cgrid, idVbcs, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, BBL(ng)%bvstrcwmax,                  &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STAname(ng),                    &
     &                        TRIM(Vname(1,idVbcs)), psta,              &
     &                        (/1,tSTAindx(ng)/), (/Nstation(ng),1/),   &
     &                        ncid = ncSTAid(ng),                       &
     &                        varid = staVid(idVbcs,ng))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out wind-induced, bed wave orbital U-velocity.
!
      IF (Sout(idUbot,ng)) THEN
        scale=1.0_r8
        CALL extract_sta2d (ng, iNLM, Cgrid, idUbot, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, BBL(ng)%Ubot,                        &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STAname(ng),                    &
     &                        TRIM(Vname(1,idUbot)), psta,              &
     &                        (/1,tSTAindx(ng)/), (/Nstation(ng),1/),   &
     &                        ncid = ncSTAid(ng),                       &
     &                        varid = staVid(idUbot,ng))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out wind-induced, bed wave orbital V-velocity.
!
      IF (Sout(idVbot,ng)) THEN
        scale=1.0_r8
        CALL extract_sta2d (ng, iNLM, Cgrid, idVbot, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, BBL(ng)%Vbot,                        &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STAname(ng),                    &
     &                        TRIM(Vname(1,idVbot)), psta,              &
     &                        (/1,tSTAindx(ng)/), (/Nstation(ng),1/),   &
     &                        ncid = ncSTAid(ng),                       &
     &                        varid = staVid(idVbot,ng))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out bottom U-velocity above bed.
!
      IF (Sout(idUbur,ng)) THEN
        scale=1.0_r8
        CALL extract_sta2d (ng, iNLM, Cgrid, idUbur, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, BBL(ng)%Ur,                          &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STAname(ng),                    &
     &                        TRIM(Vname(1,idUbur)), psta,              &
     &                        (/1,tSTAindx(ng)/), (/Nstation(ng),1/),   &
     &                        ncid = ncSTAid(ng),                       &
     &                        varid = staVid(idUbur,ng))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Write out bottom V-velocity above bed.
!
      IF (Sout(idVbvr,ng)) THEN
        scale=1.0_r8
        CALL extract_sta2d (ng, iNLM, Cgrid, idVbvr, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, BBL(ng)%Vr,                          &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STAname(ng),                    &
     &                        TRIM(Vname(1,idVbvr)), psta,              &
     &                        (/1,tSTAindx(ng)/), (/Nstation(ng),1/),   &
     &                        ncid = ncSTAid(ng),                       &
     &                        varid = staVid(idVbvr,ng))
        IF (exit_flag.ne.NoError) RETURN
      END IF
# endif

# ifdef SEDIMENT
!
!  Write out sediment fraction of each size class in each bed layer.
!
      DO i=1,NST
        IF (Sout(idfrac(i),ng)) THEN
          scale=1.0_r8
          CALL extract_sta3d (ng, iNLM, Cgrid, idfrac(i), b3dvar,       &
     &                        LBi, UBi, LBj, UBj, 1, Nbed,              &
     &                        scale, OCEAN(ng)%bed_frac(:,:,:,i),       &
     &                        NposB, XposB, YposB, ZposB, bsta)
          CALL netcdf_put_fvar (ng, iNLM, STAname(ng),                  &
     &                          TRIM(Vname(1,idfrac(i))), rsta,         &
     &                          (/1,1,tSTAindx(ng)/),                   &
     &                          (/Nbed,Nstation(ng),1/),                &
     &                          ncid = ncSTAid(ng),                     &
     &                          varid = staVid(idfrac(i),ng))
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Write out sediment mass of each size class in each bed layer.
!
        IF (Sout(idBmas(i),ng)) THEN
          scale=1.0_r8
          CALL extract_sta3d (ng, iNLM, Cgrid, idBmas(i), b3dvar,       &
     &                        LBi, UBi, LBj, UBj, 1, Nbed,              &
     &                        scale,                                    &
     &                        OCEAN(ng)%bed_mass(:,:,:,NOUT,i),         &
     &                        NposB, XposB, YposB, ZposB, bsta)
          CALL netcdf_put_fvar (ng, iNLM, STAname(ng),                  &
     &                          TRIM(Vname(1,idBmas(i))), rsta,         &
     &                          (/1,1,tSTAindx(ng)/),                   &
     &                          (/Nbed,Nstation(ng),1/),                &
     &                          ncid = ncSTAid(ng),                     &
     &                          varid = staVid(idBmas(i),ng))
          IF (exit_flag.ne.NoError) RETURN
        END IF
      END DO
!
!  Write out sediment properties in each bed layer.
!
      DO i=1,MBEDP
        IF (Sout(idSbed(i),ng)) THEN
          scale=1.0_r8
          CALL extract_sta3d (ng, iNLM, Cgrid, idSbed(i), b3dvar,       &
     &                        LBi, UBi, LBj, UBj, 1, Nbed,              &
     &                        scale, OCEAN(ng)%bed(:,:,:,i),            &
     &                        NposB, XposB, YposB, ZposB, bsta)
          CALL netcdf_put_fvar (ng, iNLM, STAname(ng),                  &
     &                          TRIM(Vname(1,idSbed(i))), rsta,         &
     &                          (/1,1,tSTAindx(ng)/),                   &
     &                          (/Nbed,Nstation(ng),1/),                &
     &                          ncid = ncSTAid(ng),                     &
     &                          varid = staVid(idSbed(i),ng))
          IF (exit_flag.ne.NoError) RETURN
        END IF
      END DO
# endif

# if defined SEDIMENT || defined BBL_MODEL
!
!  Write out exposed sediment layer properties.
!
      DO i=1,MBEDP
        IF (Sout(idBott(i),ng)) THEN
          scale=1.0_r8
          CALL extract_sta2d (ng, iNLM, Cgrid, idBott(i), r2dvar,       &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        scale, OCEAN(ng)%bottom(:,:,i),           &
     &                        Nstation(ng), Xpos, Ypos, psta)
          CALL netcdf_put_fvar (ng, iNLM, STAname(ng),                  &
     &                          TRIM(Vname(1,idBott(i))), rsta,         &
     &                          (/1,tSTAindx(ng)/),                     &
     &                          (/Nstation(ng),1/),                     &
     &                          ncid = ncSTAid(ng),                     &
     &                          varid = staVid(idBott(i),ng))
          IF (exit_flag.ne.NoError) RETURN
        END IF
      END DO
# endif
#endif
!
!-----------------------------------------------------------------------
!  Synchronize stations NetCDF file to disk.
!-----------------------------------------------------------------------
!
      CALL netcdf_sync (ng, iNLM, STAname(ng), ncSTAid(ng))

#else
      SUBROUTINE wrt_station
#endif
      RETURN
      END SUBROUTINE wrt_station
