#include "cppdefs.h"
#ifdef STATIONS
      SUBROUTINE wrt_station (ng)
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine writes out data into stations NetCDF file.             !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
#ifdef BBL_MODEL
      USE mod_bbl
#endif
      USE mod_forces
      USE mod_grid
      USE mod_iounits
      USE mod_mixing
      USE mod_ncparam
      USE mod_netcdf
      USE mod_ocean
      USE mod_scalars
# ifdef SEDIMENT
      USE mod_sediment
# endif
      USE mod_stepping
!
# ifdef SOLVE3D
      USE extract_mod, ONLY : extract2d, extract3d
# else
      USE extract_mod, ONLY : extract2d
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

      integer, dimension(3) :: start, total

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

      LBi=lbound(GRID(ng)%h,DIM=1)
      UBi=ubound(GRID(ng)%h,DIM=1)
      LBj=lbound(GRID(ng)%h,DIM=2)
      UBj=ubound(GRID(ng)%h,DIM=2)
!
!-----------------------------------------------------------------------
!  Write out station data at RHO-points.
!-----------------------------------------------------------------------
!
      IF (exit_flag.ne.0) RETURN
!
!  Set time record index.
!
      tstaindx(ng)=tstaindx(ng)+1
      nrecsta(ng)=nrecsta(ng)+1
!
!  Specify parameters for writing profiles to NetCDF files. These
!  are overridden for some variables below
!
      start(1)=1
      total(1)=N(ng)
      start(2)=1
      total(2)=Nstation(ng)
      start(3)=tstaindx(ng)
      total(3)=1
!
!  Set switch to extract station data at native C-grid position (TRUE)
!  or at RHO-points (FALSE).
!
# ifdef STATIONS_CGRID
      Cgrid=.true.
# else
      Cgrid=.false.
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
      IF (Master) THEN
        status=nf_put_var1_TYPE(ncstaid(ng), staVid(idtime,ng),         &
     &                          tstaindx(ng), time(ng))
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,10) TRIM(Vname(1,idtime)), tstaindx(ng)
          exit_flag=3
          RETURN
        END IF
      END IF
!
!  Write out free-surface (m).
!
      IF (Sout(idFsur,ng)) THEN
        scale=1.0_r8
        CALL extract2d (ng, Cgrid, idFsur, r2dvar,                      &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  scale, OCEAN(ng)%zeta(:,:,KOUT),                &
     &                  Nstation(ng), Xpos, Ypos, psta)
        IF (Master) THEN
          status=nf_put_vara_TYPE(ncstaid(ng), staVid(idFsur,ng),       &
     &                            start(2), total(2), psta)
          IF (status.ne.nf_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idFsur)), tstaindx(ng)
            exit_flag=3
            RETURN
          END IF
        END IF
      END IF
!
!  Write out 2D momentum component (m/s) in the XI-direction.
!
      IF (Sout(idUbar,ng)) THEN
        scale=1.0_r8
        CALL extract2d (ng, Cgrid, idUbar, u2dvar,                      &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  scale, OCEAN(ng)%ubar(:,:,KOUT),                &
     &                  Nstation(ng), Xpos, Ypos, psta)
        IF (Master) THEN
          status=nf_put_vara_TYPE(ncstaid(ng), staVid(idUbar,ng),       &
     &                            start(2), total(2), psta)
          IF (status.ne.nf_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbar)), tstaindx(ng)
            exit_flag=3
            RETURN
          END IF
        END IF
      END IF
!
!  Write out 2D momentum component (m/s) in the ETA-direction.
!
      IF (Sout(idVbar,ng)) THEN
        scale=1.0_r8
        CALL extract2d (ng, Cgrid, idVbar, v2dvar,                      &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  scale, OCEAN(ng)%vbar(:,:,KOUT),                &
     &                  Nstation(ng), Xpos, Ypos, psta)
        IF (Master) THEN
          status=nf_put_vara_TYPE(ncstaid(ng), staVid(idVbar,ng),       &
     &                            start(2), total(2), psta)
          IF (status.ne.nf_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbar)), tstaindx(ng)
            exit_flag=3
            RETURN
          END IF
        END IF
      END IF
# ifdef SOLVE3D
!
!  Write out 3D momentum component (m/s) in the XI-direction.
!
      IF (Sout(idUvel,ng)) THEN
        scale=1.0_r8
        CALL extract3d (ng, Cgrid, idUvel, u3dvar,                      &
     &                  LBi, UBi, LBj, UBj, 1, N(ng),                   &
     &                  scale, OCEAN(ng)%u(:,:,:,NOUT),                 &
     &                  NposR, XposR, YposR, ZposR, rsta)
        IF (Master) THEN
          status=nf_put_vara_TYPE(ncstaid(ng), staVid(idUvel,ng),       &
     &                            start, total, rsta)
          IF (status.ne.nf_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUvel)), tstaindx(ng)
            exit_flag=3
            RETURN
          END IF
        END IF
      END IF
!
!  Write out 3D momentum component (m/s) in the ETA-direction.
!
      IF (Sout(idVvel,ng)) THEN
        scale=1.0_r8
        CALL extract3d (ng, Cgrid, idVvel, v3dvar,                      &
     &                  LBi, UBi, LBj, UBj, 1, N(ng),                   &
     &                  scale, OCEAN(ng)%v(:,:,:,NOUT),                 &
     &                  NposR, XposR, YposR, ZposR, rsta)
        IF (Master) THEN
          status=nf_put_vara_TYPE(ncstaid(ng), staVid(idVvel,ng),       &
     &                            start, total, rsta)
          IF (status.ne.nf_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVvel)), tstaindx(ng)
            exit_flag=3
            RETURN
          END IF
        END IF
      END IF
!
!  Write out vertical velocity (m/s).
!
      IF (Sout(idWvel,ng)) THEN
        scale=1.0_r8
        CALL extract3d (ng, Cgrid, idWvel, w3dvar,                      &
     &                  LBi, UBi, LBj, UBj, 0, N(ng),                   &
     &                  scale, OCEAN(ng)%wvel,                          &
     &                  NposW, XposW, YposW, ZposW, rsta)
        IF (Master) THEN
          status=nf_put_vara_TYPE(ncstaid(ng), staVid(idWvel,ng),       &
     &                            start, total+(/1,0,0/), rsta)
          IF (status.ne.nf_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWvel)), tstaindx(ng)
            exit_flag=3
            RETURN
          END IF
        END IF
      END IF
!
!  Write S-coordinate "omega" vertical velocity (m3/s).
!
      IF (Sout(idOvel,ng)) THEN
        scale=1.0_r8
        CALL extract3d (ng, Cgrid, idOvel, w3dvar,                      &
     &                  LBi, UBi, LBj, UBj, 0, N(ng),                   &
     &                  scale, OCEAN(ng)%W,                             &
     &                  NposW, XposW, YposW, ZposW, rsta)
        IF (Master) THEN
          status=nf_put_vara_TYPE(ncstaid(ng), staVid(idOvel,ng),       &
     &                            start, total+(/1,0,0/), rsta)
          IF (status.ne.nf_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idOvel)), tstaindx(ng)
            exit_flag=3
            RETURN
          END IF
        END IF
      END IF
!
!  Write out tracer type variables.
!
      DO i=1,NT(ng)
        ifield=idTvar(i)
        IF (Sout(ifield,ng)) THEN
          scale=1.0_r8
          CALL extract3d (ng, Cgrid, ifield, r3dvar,                    &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    scale, OCEAN(ng)%t(:,:,:,NOUT,i),             &
     &                    NposR, XposR, YposR, ZposR, rsta)
          IF (Master) THEN
            status=nf_put_vara_TYPE(ncstaid(ng), staTid(i,ng),          &
     &                              start, total, rsta)
            IF (status.ne.nf_noerr) THEN
              WRITE (stdout,10) TRIM(Vname(1,idTvar(i))), tstaindx(ng)
              exit_flag=3
              RETURN
            END IF
          END IF
        END IF
      END DO
!
!  Write out density anomaly.
!
      IF (Sout(idDano,ng)) THEN
        scale=1.0_r8
        CALL extract3d (ng, Cgrid, idDano, r3dvar,                      &
     &                  LBi, UBi, LBj, UBj, 1, N(ng),                   &
     &                  scale, OCEAN(ng)%rho,                           &
     &                  NposR, XposR, YposR, ZposR, rsta)
        IF (Master) THEN
          status=nf_put_vara_TYPE(ncstaid(ng), staVid(idDano,ng),       &
     &                            start, total, rsta)
          IF (status.ne.nf_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idDano)), tstaindx(ng)
            exit_flag=3
            RETURN
          END IF
        END IF
      END IF
#  ifdef LMD_SKPP
!
!  Write out depth of surface boundary layer.
!
      IF (Sout(idHsbl,ng)) THEN
        scale=1.0_r8
        CALL extract2d (ng, Cgrid, idHsbl, r2dvar,                      &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  scale, MIXING(ng)%hsbl,                         &
     &                  Nstation(ng), Xpos, Ypos, psta)
        IF (Master) THEN
          status=nf_put_vara_TYPE(ncstaid(ng), staVid(idHsbl,ng),       &
     &                            start(2), total(2), psta)
          IF (status.ne.nf_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idHsbl)), tstaindx(ng)
            exit_flag=3
            RETURN
          END IF
        END IF
      END IF
#  endif
#  ifdef LMD_BKPP
!
!  Write out depth of bottom boundary layer.
!
      IF (Sout(idHbbl,ng)) THEN
        scale=1.0_r8
        CALL extract2d (ng, Cgrid, idHbbl, r2dvar,                      &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  scale, MIXING(ng)%hbbl,                         &
     &                  Nstation(ng), Xpos, Ypos, psta)
        IF (Master) THEN
          status=nf_put_vara_TYPE(ncstaid(ng), staVid(idHbbl,ng),       &
     &                            start(2), total(2), psta)
          IF (status.ne.nf_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idHbbl)), tstaindx(ng)
            exit_flag=3
            RETURN
          END IF
        END IF
      END IF
#  endif
!
!  Write out vertical viscosity coefficient.
!
      IF (Sout(idVvis,ng)) THEN
        scale=1.0_r8
        CALL extract3d (ng, Cgrid, idVvis, w3dvar,                      &
     &                  LBi, UBi, LBj, UBj, 0, N(ng),                   &
     &                  scale, MIXING(ng)%Akv,                          &
     &                  NposW, XposW, YposW, ZposW, rsta)
        IF (Master) THEN
          status=nf_put_vara_TYPE(ncstaid(ng), staVid(idVvis,ng),       &
     &                            start, total+(/1,0,0/), rsta)
          IF (status.ne.nf_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVvis)), tstaindx(ng)
            exit_flag=3
            RETURN
          END IF
        END IF
      END IF
!
!  Write out vertical diffusion coefficient for potential temperature.
!
      IF (Sout(idTdif,ng)) THEN
        scale=1.0_r8
        CALL extract3d (ng, Cgrid, idTdif, w3dvar,                      &
     &                  LBi, UBi, LBj, UBj, 0, N(ng),                   &
     &                  scale, MIXING(ng)%Akt(:,:,:,itemp),             &
     &                  NposW, XposW, YposW, ZposW, rsta)
        IF (Master) THEN
          status=nf_put_vara_TYPE(ncstaid(ng), staVid(idTdif,ng),       &
     &                            start, total+(/1,0,0/), rsta)
          IF (status.ne.nf_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTdif)), tstaindx(ng)
            exit_flag=3
            RETURN
          END IF
        END IF
      END IF
#  ifdef SALINITY
!
!  Write out vertical diffusion coefficient for salinity.
!
      IF (Sout(idSdif,ng)) THEN
        scale=1.0_r8
        CALL extract3d (ng, Cgrid, idSdif, w3dvar,                      &
     &                  LBi, UBi, LBj, UBj, 0, N(ng),                   &
     &                  scale, MIXING(ng)%Akt(:,:,:,isalt),             &
     &                  NposW, XposW, YposW, ZposW, rsta)
        IF (Master) THEN
          status=nf_put_vara_TYPE(ncstaid(ng), staVid(idSdif,ng),       &
     &                            start, total+(/1,0,0/), rsta)
          IF (status.ne.nf_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSdif)), tstaindx(ng)
            exit_flag=3
            RETURN
          END IF
        END IF
      END IF
#  endif
#  if defined MY25_MIXING || defined GLS_MIXING
!
!  Write out turbulent kinetic energy.
!
      IF (Sout(idMtke,ng)) THEN
        scale=1.0_r8
        CALL extract3d (ng, Cgrid, idMtke, w3dvar,                      &
     &                  LBi, UBi, LBj, UBj, 0, N(ng),                   &
     &                  scale, MIXING(ng)%tke(:,:,:,NOUT),              &
     &                  NposW, XposW, YposW, ZposW, rsta)
        IF (Master) THEN
          status=nf_put_vara_TYPE(ncstaid(ng), staVid(idMtke,ng),       &
     &                            start, total+(/1,0,0/), rsta)
          IF (status.ne.nf_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idMtke)), tstaindx(ng)
            exit_flag=3
            RETURN
          END IF
        END IF
      END IF
!
!  Write out turbulent kinetic energy times length scale.
!
      IF (Sout(idMtls,ng)) THEN
        scale=1.0_r8
        CALL extract3d (ng, Cgrid, idMtls, w3dvar,                      &
     &                  LBi, UBi, LBj, UBj, 0, N(ng),                   &
     &                  scale, MIXING(ng)%gls(:,:,:,NOUT),              &
     &                  NposW, XposW, YposW, ZposW, rsta)
        IF (Master) THEN
          status=nf_put_vara_TYPE(ncstaid(ng), staVid(idMtls,ng),       &
     &                            start, total+(/1,0,0/), rsta)
          IF (status.ne.nf_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idMtls)), tstaindx(ng)
            exit_flag=3
            RETURN
          END IF
        END IF
      END IF
#  endif
!
!  Write out surface net heat flux.
!
      IF (Sout(idTsur(itemp),ng)) THEN
        ifield=idTsur(itemp)
        scale=rho0*Cp
        CALL extract2d (ng, Cgrid, ifield, r2dvar,                      &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  scale, FORCES(ng)%stflx(:,:,itemp),             &
     &                  Nstation(ng), Xpos, Ypos, psta)
        IF (Master) THEN
          status=nf_put_vara_TYPE(ncstaid(ng), staVid(ifield,ng),       &
     &                            start(2), total(2), psta)
          IF (status.ne.nf_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,ifield)), tstaindx(ng)
            exit_flag=3
            RETURN
          END IF
        END IF
      END IF
#  ifdef BULK_FLUXES
!
!  Write out latent heat flux.
!
      IF (Sout(idLhea,ng)) THEN
        scale=rho0*Cp
        CALL extract2d (ng, Cgrid, idLhea, r2dvar,                      &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  scale, FORCES(ng)%lhflx,                        &
     &                  Nstation(ng), Xpos, Ypos, psta)
        IF (Master) THEN
          status=nf_put_vara_TYPE(ncstaid(ng), staVid(idLhea,ng),       &
     &                            start(2), total(2), psta)
          IF (status.ne.nf_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idLhea)), tstaindx(ng)
            exit_flag=3
            RETURN
          END IF
        END IF
      END IF
!
!  Write out sensible heat flux.
!
      IF (Sout(idShea,ng)) THEN
        scale=rho0*Cp
        CALL extract2d (ng, Cgrid, idShea, r2dvar,                      &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  scale, FORCES(ng)%shflx,                        &
     &                  Nstation(ng), Xpos, Ypos, psta)
        IF (Master) THEN
          status=nf_put_vara_TYPE(ncstaid(ng), staVid(idShea,ng),       &
     &                            start(2), total(2), psta)
          IF (status.ne.nf_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idShea)), tstaindx(ng)
            exit_flag=3
            RETURN
          END IF
        END IF
      END IF
!
!  Write out longwave radiation flux.
!
      IF (Sout(idLrad,ng)) THEN
        scale=rho0*Cp
        CALL extract2d (ng, Cgrid, idLrad, r2dvar,                      &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  scale, FORCES(ng)%lrflx,                        &
     &                  Nstation(ng), Xpos, Ypos, psta)
        IF (Master) THEN
          status=nf_put_vara_TYPE(ncstaid(ng), staVid(idLrad,ng),       &
     &                            start(2), total(2), psta)
          IF (status.ne.nf_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idLrad)), tstaindx(ng)
            exit_flag=3
            RETURN
          END IF
        END IF
      END IF
#  endif
#  ifdef SHORTWAVE
!
!  Write out shortwave radiation flux.
!
      IF (Sout(idSrad,ng)) THEN
        scale=rho0*Cp
        CALL extract2d (ng, Cgrid, idSrad, r2dvar,                      &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  scale, FORCES(ng)%srflx,                        &
     &                  Nstation(ng), Xpos, Ypos, psta)
        IF (Master) THEN
          status=nf_put_vara_TYPE(ncstaid(ng), staVid(idSrad,ng),       &
     &                            start(2), total(2), psta)
          IF (status.ne.nf_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSrad)), tstaindx(ng)
            exit_flag=3
            RETURN
          END IF
        END IF
      END IF
#  endif
# endif
!
!  Write out surface U-momentum stress.
!
      IF (Sout(idUsms,ng)) THEN
        scale=rho0
        CALL extract2d (ng, Cgrid, idUsms, u2dvar,                      &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  scale, FORCES(ng)%sustr,                        &
     &                  Nstation(ng), Xpos, Ypos, psta)
        IF (Master) THEN
          status=nf_put_vara_TYPE(ncstaid(ng), staVid(idUsms,ng),       &
     &                            start(2), total(2), psta)
          IF (status.ne.nf_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUsms)), tstaindx(ng)
            exit_flag=3
            RETURN
          END IF
        END IF
      END IF
!
!  Write out surface V-momentum stress.
!
      IF (Sout(idVsms,ng)) THEN
        scale=rho0
        CALL extract2d (ng, Cgrid, idVsms, v2dvar,                      &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  scale, FORCES(ng)%svstr,                        &
     &                  Nstation(ng), Xpos, Ypos, psta)
        IF (Master) THEN
          status=nf_put_vara_TYPE(ncstaid(ng), staVid(idVsms,ng),       &
     &                            start(2), total(2), psta)
          IF (status.ne.nf_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVsms)), tstaindx(ng)
            exit_flag=3
            RETURN
          END IF
        END IF
      END IF
!
!  Write out bottom U-momentum stress.
!
      IF (Sout(idUbms,ng)) THEN
        scale=-rho0
        CALL extract2d (ng, Cgrid, idUbms, u2dvar,                      &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  scale, FORCES(ng)%bustr,                        &
     &                  Nstation(ng), Xpos, Ypos, psta)
        IF (Master) THEN
          status=nf_put_vara_TYPE(ncstaid(ng), staVid(idUbms,ng),       &
     &                            start(2), total(2), psta)
          IF (status.ne.nf_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbms)), tstaindx(ng)
            exit_flag=3
            RETURN
          END IF
        END IF
      END IF
!
!  Write out bottom V-momentum stress.
!
      IF (Sout(idVbms,ng)) THEN
        scale=-rho0
        CALL extract2d (ng, Cgrid, idVbms, v2dvar,                      &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  scale, FORCES(ng)%bvstr,                        &
     &                  Nstation(ng), Xpos, Ypos, psta)
        IF (Master) THEN
          status=nf_put_vara_TYPE(ncstaid(ng), staVid(idVbms,ng),       &
     &                            start(2), total(2), psta)
          IF (status.ne.nf_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbms)), tstaindx(ng)
            exit_flag=3
            RETURN
          END IF
        END IF
      END IF
#ifdef SOLVE3D
# ifdef BBL_MODEL
!
!  Write out wind-induced, bottom U-wave stress.
!
      IF (Sout(idUbws,ng)) THEN
        scale=rho0
        CALL extract2d (ng, Cgrid, idUbws, u2dvar,                      &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  scale, BBL(ng)%bustrw,                          &
     &                  Nstation(ng), Xpos, Ypos, psta)
        IF (Master) THEN
          status=nf_put_vara_TYPE(ncstaid(ng), staVid(idUbws,ng),       &
     &                            start(2), total(2), psta)
          IF (status.ne.nf_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbws)), tstaindx(ng)
            exit_flag=3
            RETURN
          END IF
        END IF
      END IF
!
!  Write out wind-induced, bottom V-wave stress.
!
      IF (Sout(idVbws,ng)) THEN
        scale=rho0
        CALL extract2d (ng, Cgrid, idVbws, v2dvar,                      &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  scale, BBL(ng)%bvstrw,                          &
     &                  Nstation(ng), Xpos, Ypos, psta)
        IF (Master) THEN
          status=nf_put_vara_TYPE(ncstaid(ng), staVid(idVbws,ng),       &
     &                            start(2), total(2), psta)
          IF (status.ne.nf_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbws)), tstaindx(ng)
            exit_flag=3
            RETURN
          END IF
        END IF
      END IF
!
!  Write out maximum wind and current, bottom U-wave stress.
!
      IF (Sout(idUbcs,ng)) THEN
        scale=rho0
        CALL extract2d (ng, Cgrid, idUbcs, u2dvar,                      &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  scale, BBL(ng)%bustrcwmax,                      &
     &                  Nstation(ng), Xpos, Ypos, psta)
        IF (Master) THEN
          status=nf_put_vara_TYPE(ncstaid(ng), staVid(idUbcs,ng),       &
     &                            start(2), total(2), psta)
          IF (status.ne.nf_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbcs)), tstaindx(ng)
            exit_flag=3
            RETURN
          END IF
        END IF
      END IF
!
!  Write out maximum wind and current, bottom V-wave stress.
!
      IF (Sout(idVbcs,ng)) THEN
        scale=rho0
        CALL extract2d (ng, Cgrid, idVbcs, v2dvar,                      &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  scale, BBL(ng)%bvstrcwmax,                          &
     &                  Nstation(ng), Xpos, Ypos, psta)
        IF (Master) THEN
          status=nf_put_vara_TYPE(ncstaid(ng), staVid(idVbcs,ng),       &
     &                            start(2), total(2), psta)
          IF (status.ne.nf_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbcs)), tstaindx(ng)
            exit_flag=3
            RETURN
          END IF
        END IF
      END IF
!
!  Write out wind-induced, bed wave orbital U-velocity.
!
      IF (Sout(idUbed,ng)) THEN
        scale=1.0_r8
        CALL extract2d (ng, Cgrid, idUbed, u2dvar,                      &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  scale, BBL(ng)%Ubed,                            &
     &                  Nstation(ng), Xpos, Ypos, psta)
        IF (Master) THEN
          status=nf_put_vara_TYPE(ncstaid(ng), staVid(idUbed,ng),       &
     &                            start(2), total(2), psta)
          IF (status.ne.nf_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbed)), tstaindx(ng)
            exit_flag=3
            RETURN
          END IF
        END IF
      END IF
!
!  Write out wind-induced, bed wave orbital V-velocity.
!
      IF (Sout(idVbed,ng)) THEN
        scale=1.0_r8
        CALL extract2d (ng, Cgrid, idVbed, v2dvar,                      &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  scale, BBL(ng)%Vbed,                            &
     &                  Nstation(ng), Xpos, Ypos, psta)
        IF (Master) THEN
          status=nf_put_vara_TYPE(ncstaid(ng), staVid(idVbed,ng),       &
     &                            start(2), total(2), psta)
          IF (status.ne.nf_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbed)), tstaindx(ng)
            exit_flag=3
            RETURN
          END IF
        END IF
      END IF
!
!  Write out bottom U-momentum above bed.
!
      IF (Sout(idUbot,ng)) THEN
        scale=1.0_r8
        CALL extract2d (ng, Cgrid, idUbot, u2dvar,                      &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  scale, BBL(ng)%Ubot,                            &
     &                  Nstation(ng), Xpos, Ypos, psta)
        IF (Master) THEN
          status=nf_put_vara_TYPE(ncstaid(ng), staVid(idUbot,ng),       &
     &                            start(2), total(2), psta)
          IF (status.ne.nf_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbot)), tstaindx(ng)
            exit_flag=3
            RETURN
          END IF
        END IF
      END IF
!
!  Write out bottom V-momentum above bed.
!
      IF (Sout(idVbot,ng)) THEN
        scale=1.0_r8
        CALL extract2d (ng, Cgrid, idVbot, v2dvar,                      &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  scale, BBL(ng)%Vbot,                            &
     &                  Nstation(ng), Xpos, Ypos, psta)
        IF (Master) THEN
          status=nf_put_vara_TYPE(ncstaid(ng), staVid(idVbot,ng),       &
     &                            start(2), total(2), psta)
          IF (status.ne.nf_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbot)), tstaindx(ng)
            exit_flag=3
            RETURN
          END IF
        END IF
      END IF
# endif
# ifdef SEDIMENT
!
!  Write out sediment fraction of each size class in each bed layer.
!
      total(1)=Nbed                       ! Warning, total(1) redefined
      DO i=1,NST      
        IF (Sout(idfrac(i),ng)) THEN
          scale=1.0_r8
          CALL extract3d (ng, Cgrid, idfrac(i), b3dvar,                 &
     &                    LBi, UBi, LBj, UBj, 1, Nbed,                  &
     &                    scale, OCEAN(ng)%bed_frac(:,:,:,i),           &
     &                    NposB, XposB, YposB, ZposB, bsta)
          IF (Master) THEN
            status=nf_put_vara_TYPE(ncstaid(ng), staVid(idfrac(i),ng),  &
     &                              start, total, bsta)
            IF (status.ne.nf_noerr) THEN
              WRITE (stdout,10) TRIM(Vname(1,idfrac(i))), tstaindx(ng)
              exit_flag=3
              RETURN
            END IF
          END IF
        END IF
!
!  Write out sediment mass of each size class in each bed layer.
!
        IF (Sout(idBmas(i),ng)) THEN
          scale=1.0_r8
          CALL extract3d (ng, Cgrid, idBmas(i), b3dvar,                 &
     &                    LBi, UBi, LBj, UBj, 1, Nbed,                  &
     &                    scale, OCEAN(ng)%bed_mass(:,:,:,i),           &
     &                    NposB, XposB, YposB, ZposB, bsta)
          IF (Master) THEN
            status=nf_put_vara_TYPE(ncstaid(ng), staVid(idBmas(i),ng),  &
     &                              start, total, bsta)
            IF (status.ne.nf_noerr) THEN
              WRITE (stdout,10) TRIM(Vname(1,idBmas(i))), tstaindx(ng)
              exit_flag=3
              RETURN
            END IF
          END IF
        END IF
      END DO
!
!  Write out sediment properties in each bed layer.
!
      DO i=1,MBEDP      
        IF (Sout(idSbed(i),ng)) THEN
          scale=1.0_r8
          CALL extract3d (ng, Cgrid, idSbed(i), b3dvar,                 &
     &                    LBi, UBi, LBj, UBj, 1, Nbed,                  &
     &                    scale, OCEAN(ng)%bed(:,:,:,i),                &
     &                    NposB, XposB, YposB, ZposB, bsta)
          IF (Master) THEN
            status=nf_put_vara_TYPE(ncstaid(ng), staVid(idSbed(i),ng),  &
     &                              start, total, bsta)
            IF (status.ne.nf_noerr) THEN
              WRITE (stdout,10) TRIM(Vname(1,idSbed(i))), tstaindx(ng)
              exit_flag=3
              RETURN
            END IF
          END IF
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
          CALL extract2d (ng, Cgrid, idBott(i), r2dvar,                 &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    scale, OCEAN(ng)%bottom(:,;,i),               &
     &                    Nstation(ng), Xpos, Ypos, psta)
          IF (Master) THEN
            status=nf_put_vara_TYPE(ncstaid(ng), staVid(idBott(i),ng),  &
     &                              start(2), total(2), psta)
            IF (status.ne.nf_noerr) THEN
              WRITE (stdout,10) TRIM(Vname(1,idBott(i))), tstaindx(ng)
              exit_flag=3
              RETURN
            END IF
          END IF
        END IF
      END DO
# endif
#endif
!
!  Synchronize restart NetCDF file to disk.
!
      IF (Master) THEN
        status=nf_sync(ncstaid(ng))
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,20)
          exit_flag=3
          RETURN
        END IF
      END IF
!
  10  FORMAT (/,' WRT_STATION - error while writing variable: ',a,/,    &
     &        15x,'into stations NetCDF file for time record: ',i4)
  20  FORMAT (/,' WRT_STATION - unable to synchronize stations',        &
     &        1x,'NetCDF file to disk.')
#else
      SUBROUTINE wrt_station
#endif
      RETURN
      END SUBROUTINE wrt_station
