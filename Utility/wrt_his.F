#include "cppdefs.h"
      SUBROUTINE wrt_his (ng)
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine writes requested model fields at requested levels      !
!  into history NetCDF file.                                           !
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
#ifdef SEDIMENT
      USE mod_sediment
#endif
      USE mod_stepping
!
#ifdef SOLVE3D
      USE omega_mod, ONLY : scale_omega
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
      integer :: LBi, UBi, LBj, UBj
      integer :: gfactor, gtype, status
#ifdef SOLVE3D
      integer :: i, itrc, j, k, tile
#endif

      integer :: nf_fwrite

      real(r8) :: scale

#ifdef SOLVE3D
      real(r8), allocatable :: wrk(:,:,:)
#endif
!
      LBi=lbound(GRID(ng)%h,DIM=1)
      UBi=ubound(GRID(ng)%h,DIM=1)
      LBj=lbound(GRID(ng)%h,DIM=2)
      UBj=ubound(GRID(ng)%h,DIM=2)
!
!-----------------------------------------------------------------------
!  Write out history fields.
!-----------------------------------------------------------------------
!
      if (exit_flag.ne.0) RETURN
!
!  Set grid type factor to write full (gfactor=1) fields or water
!  points (gfactor=-1) fields only.
!
#if defined WRITE_WATER && defined MASKING
        gfactor=-1
#else
        gfactor=1
#endif
!
!  Set time record index.
!
      thisindx(ng)=thisindx(ng)+1
      nrechis(ng)=nrechis(ng)+1
!
!  Write out model time (s).
!
      IF (OutThread) THEN
        status=nf_put_var1_TYPE(nchisid(ng), hisVid(idtime,ng),         &
     &                          thisindx(ng),time(ng))
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,10) TRIM(Vname(1,idtime)), thisindx(ng)
          exit_flag=3
          RETURN
        END IF
      END IF
!
!  Write out free-surface (m)
!
      IF (Hout(idFsur,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite(ng, nchisid(ng), hisVid(idFsur,ng),            &
     &                   thisindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
#ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#endif
     &                   OCEAN(ng) % zeta(LBi,LBj,KOUT))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idFsur)), thisindx(ng)
          exit_flag=3
          RETURN
        END IF
      END IF
!
!  Write out 2D momentum component (m/s) in the XI-direction.
!
      IF (Hout(idUbar,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u2dvar
        status=nf_fwrite(ng, nchisid(ng), hisVid(idUbar,ng),            &
     &                   thisindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
#ifdef MASKING
     &                   GRID(ng) % umask(LBi,LBj),                     &
#endif
     &                   OCEAN(ng) % ubar(LBi,LBj,KOUT))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUbar)), thisindx(ng)
          exit_flag=3
          RETURN
        END IF
      END IF
!
!  Write out 2D momentum component (m/s) in the ETA-direction.
!
      IF (Hout(idVbar,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v2dvar
        status=nf_fwrite(ng, nchisid(ng), hisVid(idVbar,ng),            &
     &                   thisindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
#ifdef MASKING
     &                   GRID(ng) % vmask(LBi,LBj),                     &
#endif
     &                   OCEAN(ng) % vbar(LBi,LBj,KOUT))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVbar)), thisindx(ng)
          exit_flag=3
          RETURN
        END IF
      END IF
#ifdef SOLVE3D
!
!  Write out 3D momentum component (m/s) in the XI-direction.
!
      IF (Hout(idUvel,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u3dvar
        status=nf_fwrite(ng, nchisid(ng), hisVid(idUvel,ng),            &
     &                   thisindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
# ifdef MASKING
     &                   GRID(ng) % umask(LBi,LBj),                     &
# endif
     &                   OCEAN(ng) % u(LBi,LBj,1,NOUT))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUvel)), thisindx(ng)
          exit_flag=3
          RETURN
        END IF
      END IF
!
!  Write out 3D momentum component (m/s) in the ETA-direction.
!
      IF (Hout(idVvel,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v3dvar
        status=nf_fwrite(ng, nchisid(ng), hisVid(idVvel,ng),            &
     &                   thisindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
# ifdef MASKING
     &                   GRID(ng) % vmask(LBi,LBj),                     &
# endif
     &                   OCEAN(ng) % v(LBi,LBj,1,NOUT))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVvel)), thisindx(ng)
          exit_flag=3
          RETURN
        END IF
      END IF
!
!  Write out S-coordinate omega vertical velocity (m/s).
!
      IF (Hout(idOvel,ng)) THEN
        IF (.not.allocated(wrk)) THEN
          allocate (wrk(LBi:UBi,LBj:UBj,0:N(ng)))
          wrk(LBi:UBi,LBj:UBj,0:N(ng))=0.0_r8
        END IF
        scale=1.0_r8
        gtype=gfactor*w3dvar
        DO tile=0,NtileX(ng)*NtileE(ng)-1
          CALL scale_omega (ng, TILE, LBi, UBi, LBj, UBj, 0, N(ng),     &
     &                      GRID(ng) % pm,                              &
     &                      GRID(ng) % pn,                              &
     &                      OCEAN(ng) % W,                              &
     &                      wrk)
        END DO
        status=nf_fwrite(ng, nchisid(ng), hisVid(idOvel,ng),            &
     &                   thisindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 0, N(ng), scale,           &
# ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
# endif
     &                   wrk(LBi,LBj,0))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idOvel)), thisindx(ng)
          exit_flag=3
          RETURN
        END IF
        deallocate (wrk)
      END IF
!
!  Write out vertical velocity (m/s).
!
      IF (Hout(idWvel,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*w3dvar
        status=nf_fwrite(ng, nchisid(ng), hisVid(idWvel,ng),            &
     &                   thisindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 0, N(ng), scale,           &
# ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
# endif
     &                   OCEAN(ng) % wvel(LBi,LBj,0))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idWvel)), thisindx(ng)
          exit_flag=3
          RETURN
        END IF
      END IF
!
!  Write out tracer type variables.
!
      DO itrc=1,NT(ng)
        IF (Hout(idTvar(itrc),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*r3dvar
          status=nf_fwrite(ng, nchisid(ng), hisTid(itrc,ng),            &
     &                     thisindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
# ifdef MASKING
     &                     GRID(ng) % rmask(LBi,LBj),                   &
# endif
     &                     OCEAN(ng) % t(LBi,LBj,1,NOUT,itrc))
          IF (OutThread.and.(status.ne.nf_noerr)) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTvar(itrc))),              &
     &                        thisindx(ng)
            exit_flag=3
            RETURN
          END IF
        END IF
      END DO
!
!  Write out density anomaly.
!
      IF (Hout(idDano,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r3dvar
        status=nf_fwrite(ng, nchisid(ng), hisVid(idDano,ng),            &
     &                   thisindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
# ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
# endif
     &                   OCEAN(ng) % rho(LBi,LBj,1))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idDano)), thisindx(ng)
          exit_flag=3
          RETURN
        END IF
      END IF
# ifdef LMD_SKPP
!
!  Write out depth surface boundary layer.
!
      IF (Hout(idHsbl,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite(ng, nchisid(ng), hisVid(idHsbl,ng),            &
     &                   thisindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
#  ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#  endif
     &                   MIXING(ng) % hsbl(LBi,LBj))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idHsbl)), thisindx(ng)
          exit_flag=3
          RETURN
        END IF
      END IF
# endif
# ifdef LMD_BKPP
!
!  Write out depth surface boundary layer.
!
      IF (Hout(idHbbl,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite(ng, nchisid(ng), hisVid(idHbbl,ng),            &
     &                   thisindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
#  ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#  endif
     &                   MIXING(ng) % hbbl(LBi,LBj))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idHbbl)), thisindx(ng)
          exit_flag=3
          RETURN
        END IF
      END IF
# endif
!
!  Write out vertical viscosity coefficient.
!
      IF (Hout(idVvis,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*w3dvar
        status=nf_fwrite(ng, nchisid(ng), hisVid(idVvis,ng),            &
     &                   thisindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 0, N(ng), scale,           &
# ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
# endif
     &                   MIXING(ng) % Akv(LBi,LBj,0))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVvis)), thisindx(ng)
          exit_flag=3
          RETURN
        END IF
      END IF
!
!  Write out vertical diffusion coefficient for potential temperature.
!
      IF (Hout(idTdif,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*w3dvar
        status=nf_fwrite(ng, nchisid(ng), hisVid(idTdif,ng),            &
     &                   thisindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 0, N(ng), scale,           &
# ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
# endif
     &                   MIXING(ng) % Akt(LBi,LBj,0,itemp))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idTdif)), thisindx(ng)
          exit_flag=3
          RETURN
        END IF
      END IF
# ifdef SALINITY
!
!  Write out vertical diffusion coefficient for salinity.
!
      IF (Hout(idSdif,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*w3dvar
        status=nf_fwrite(ng, nchisid(ng), hisVid(idSdif,ng),            &
     &                   thisindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 0, N(ng), scale,           &
#  ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#  endif
     &                   MIXING(ng) % Akt(LBi,LBj,0,isalt))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idSdif)), thisindx(ng)
          exit_flag=3
          RETURN
        END IF
      END IF
# endif
# if defined MY25_MIXING || defined GLS_MIXING
!
!  Write out turbulent kinetic energy.
!
      IF (Hout(idMtke,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*w3dvar
        status=nf_fwrite(ng, nchisid(ng), hisVid(idMtke,ng),            &
     &                   thisindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 0, N(ng), scale,           &
#  ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#  endif
     &                   MIXING(ng) % tke(LBi,LBj,0,NOUT))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idMtke)), thisindx(ng)
          exit_flag=3
          RETURN
        END IF
      END IF
!
!  Write out turbulent length scale field.
!
      IF (Hout(idMtls,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*w3dvar
        status=nf_fwrite(ng, nchisid(ng), hisVid(idMtls,ng),            &
     &                   thisindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 0, N(ng), scale,           &
#  ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#  endif
     &                   MIXING(ng) % gls(LBi,LBj,0,NOUT))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idMtls)), thisindx(ng)
          exit_flag=3
          RETURN
        END IF
      END IF
# endif
!
!  Write out surface net heat flux.
!
      IF (Hout(idTsur(itemp),ng)) THEN
        scale=rho0*Cp
        gtype=gfactor*r2dvar
        status=nf_fwrite(ng, nchisid(ng), hisVid(idTsur(itemp),ng),     &
     &                   thisindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
# ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
# endif
     &                   FORCES(ng) % stflx(LBi,LBj,itemp))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idTsur(itemp))),               &
     &                      thisindx(ng)
          exit_flag=3
          RETURN
        END IF
      END IF
# ifdef BULK_FLUXES
!
!  Write out latent heat flux.
!
      IF (Hout(idLhea,ng)) THEN
        scale=rho0*Cp
        gtype=gfactor*r2dvar
        status=nf_fwrite(ng, nchisid(ng), hisVid(idLhea,ng),            &
     &                   thisindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
#  ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#  endif
     &                   FORCES(ng) % lhflx(LBi,LBj))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idLhea)), thisindx(ng)
          exit_flag=3
          RETURN
        END IF
      END IF
!
!  Write out sensible heat flux.
!
      IF (Hout(idShea,ng)) THEN
        scale=rho0*Cp
        gtype=gfactor*r2dvar
        status=nf_fwrite(ng, nchisid(ng), hisVid(idShea,ng),            &
     &                   thisindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
#  ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#  endif
     &                   FORCES(ng) % shflx(LBi,LBj))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idShea)), thisindx(ng)
          exit_flag=3
          RETURN
        END IF
      END IF
!
!  Write out longwave radiation flux.
!
      IF (Hout(idLrad,ng)) THEN
        scale=rho0*Cp
        gtype=gfactor*r2dvar
        status=nf_fwrite(ng, nchisid(ng), hisVid(idLrad,ng),            &
     &                   thisindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
#  ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#  endif
     &                   FORCES(ng) % lrflx(LBi,LBj))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idLrad)), thisindx(ng)
          exit_flag=3
          RETURN
        END IF
      END IF
# endif
# ifdef SHORTWAVE
!
!  Write out shortwave radiation flux.
!
      IF (Hout(idSrad,ng)) THEN
        scale=rho0*Cp
        gtype=gfactor*r2dvar
        status=nf_fwrite(ng, nchisid(ng), hisVid(idSrad,ng),            &
     &                   thisindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
#  ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#  endif
     &                   FORCES(ng) % srflx(LBi,LBj))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idSrad)), thisindx(ng)
          exit_flag=3
          RETURN
        END IF
      END IF
# endif
#endif
!
!  Write out surface U-momentum stress.
!
      IF (Hout(idUsms,ng)) THEN
        scale=rho0
        gtype=gfactor*u2dvar
        status=nf_fwrite(ng, nchisid(ng), hisVid(idUsms,ng),            &
     &                   thisindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
#ifdef MASKING
     &                   GRID(ng) % umask(LBi,LBj),                     &
#endif
     &                   FORCES(ng) % sustr(LBi,LBj))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUsms)), thisindx(ng)
          exit_flag=3
          RETURN
        END IF
      END IF
!
!  Write out surface V-momentum stress.
!
      IF (Hout(idVsms,ng)) THEN
        scale=rho0
        gtype=gfactor*v2dvar
        status=nf_fwrite(ng, nchisid(ng), hisVid(idVsms,ng),            &
     &                   thisindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
#ifdef MASKING
     &                   GRID(ng) % vmask(LBi,LBj),                     &
#endif
     &                   FORCES(ng) % svstr(LBi,LBj))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVsms)), thisindx(ng)
          exit_flag=3
          RETURN
        END IF
      END IF
!
!  Write out bottom U-momentum stress.
!
      IF (Hout(idUbms,ng)) THEN
        scale=-rho0
        gtype=gfactor*u2dvar
        status=nf_fwrite(ng, nchisid(ng), hisVid(idUbms,ng),            &
     &                   thisindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
#ifdef MASKING
     &                   GRID(ng) % umask(LBi,LBj),                     &
#endif
     &                   FORCES(ng) % bustr(LBi,LBj))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUbms)), thisindx(ng)
          exit_flag=3
          RETURN
        END IF
      END IF
!
!  Write out bottom V-momentum stress.
!
      IF (Hout(idVbms,ng)) THEN
        scale=-rho0
        gtype=gfactor*v2dvar
        status=nf_fwrite(ng, nchisid(ng), hisVid(idVbms,ng),            &
     &                   thisindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
#ifdef MASKING
     &                   GRID(ng) % vmask(LBi,LBj),                     &
#endif
     &                   FORCES(ng) % bvstr(LBi,LBj))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVbms)), thisindx(ng)
          exit_flag=3
          RETURN
        END IF
      END IF
#ifdef BBL_MODEL
!
!  Write out wind-induced, bottom U-wave stress.
!
      IF (Hout(idUbws,ng)) THEN
        scale=rho0
        gtype=gfactor*u2dvar
        status=nf_fwrite(ng, nchisid(ng), hisVid(idUbws,ng),            &
     &                   thisindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
# ifdef MASKING
     &                   GRID(ng) % umask(LBi,LBj),                     &
# endif
     &                   BBL(ng) % bustrw(LBi,LBj))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUbws)), thisindx(ng)
          exit_flag=3
          RETURN
        END IF
      END IF
!
!  Write out wind-induced, bottom V-wave stress.
!
      IF (Hout(idVbws,ng)) THEN
        scale=rho0
        gtype=gfactor*v2dvar
        status=nf_fwrite(ng, nchisid(ng), hisVid(idVbws,ng),            &
     &                   thisindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
# ifdef MASKING
     &                   GRID(ng) % vmask(LBi,LBj),                     &
# endif
     &                   BBL(ng) % bvstrw(LBi,LBj))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVbws)), thisindx(ng)
          exit_flag=3
          RETURN
        END IF
      END IF
!
!  Write out wind-induced, bed wave excursion amplitude.
!
      IF (Hout(idAbed,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite(ng, nchisid(ng), hisVid(idAbed,ng),            &
     &                   thisindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
# ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
# endif
     &                   BBL(ng) % Abed(LBi,LBj))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idAbed)), thisindx(ng)
          exit_flag=3
          RETURN
        END IF
      END IF
!
!  Write out wind-induced, bed wave orbital U-velocity.
!
      IF (Hout(idUbed,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u2dvar
        status=nf_fwrite(ng, nchisid(ng), hisVid(idUbed,ng),            &
     &                   thisindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
# ifdef MASKING
     &                   GRID(ng) % umask(LBi,LBj),                     &
# endif
     &                   BBL(ng) % Ubed(LBi,LBj))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUbed)), thisindx(ng)
          exit_flag=3
          RETURN
        END IF
      END IF
!
!  Write out wind-induced, bed wave orbital V-velocity.
!
      IF (Hout(idVbed,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v2dvar
        status=nf_fwrite(ng, nchisid(ng), hisVid(idVbed,ng),            &
     &                   thisindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
# ifdef MASKING
     &                   GRID(ng) % vmask(LBi,LBj),                     &
# endif
     &                   BBL(ng) % Vbed(LBi,LBj))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVbed)), thisindx(ng)
          exit_flag=3
          RETURN
        END IF
      END IF
!
!  Write out bottom U-momentum above bed.
!
      IF (Hout(idUbot,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u2dvar
        status=nf_fwrite(ng, nchisid(ng), hisVid(idUbot,ng),            &
     &                   thisindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
# ifdef MASKING
     &                   GRID(ng) % umask(LBi,LBj),                     &
# endif
     &                   BBL(ng) % Ubot(LBi,LBj))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUbot)), thisindx(ng)
          exit_flag=3
          RETURN
        END IF
      END IF
!
!  Write out bottom V-momentum above bed.
!
      IF (Hout(idVbot,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v2dvar
        status=nf_fwrite(ng, nchisid(ng), hisVid(idVbot,ng),            &
     &                   thisindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
# ifdef MASKING
     &                   GRID(ng) % vmask(LBi,LBj),                     &
# endif
     &                   BBL(ng) % Vbot(LBi,LBj))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVbot)), thisindx(ng)
          exit_flag=3
          RETURN
        END IF
      END IF
!
!  Write out bed ripple height.
!
      IF (Hout(idHrip,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite(ng, nchisid(ng), hisVid(idHrip,ng),            &
     &                   thisindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
# ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
# endif
     &                   BBL(ng) % Hripple(LBi,LBj))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idHrip)), thisindx(ng)
          exit_flag=3
          RETURN
        END IF
      END IF
!
!  Write out bed ripple length.
!
      IF (Hout(idLrip,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite(ng, nchisid(ng), hisVid(idLrip,ng),            &
     &                   thisindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
# ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
# endif
     &                   BBL(ng) % Lripple(LBi,LBj))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idLrip)), thisindx(ng)
          exit_flag=3
          RETURN
        END IF
      END IF
!
!  Write out hydraulic bottom roughness.
!
      IF (Hout(idZnot,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite(ng, nchisid(ng), hisVid(idZnot,ng),            &
     &                   thisindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
# ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
# endif
     &                   BBL(ng) % Zbnot(LBi,LBj))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idZnot)), thisindx(ng)
          exit_flag=3
          RETURN
        END IF
      END IF
!
!  Write out apparent hydraulic bottom roughness.
!
      IF (Hout(idZapp,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite(ng, nchisid(ng), hisVid(idZapp,ng),            &
     &                   thisindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
# ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
# endif
     &                   BBL(ng) % Zbnotc(LBi,LBj))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idZapp)), thisindx(ng)
          exit_flag=3
          RETURN
        END IF
      END IF
#endif
#ifdef SEDIMENT
!
!  Write out bed layer sediment fraction for each size class.
!
      DO i=1,NST
        IF (Hout(idfrac(i),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*b3dvar
          status=nf_fwrite(ng, nchisid(ng), hisVid(idfrac(i),ng),       &
     &                     thisindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, Nbed, scale,          &
# ifdef MASKING
     &                     GRID(ng) % rmask(LBi,LBj),                   &
# endif
     &                     OCEAN(ng) % bed_frac(LBi,LBj,1,i))
          IF (OutThread.and.(status.ne.nf_noerr)) THEN
            WRITE (stdout,10) TRIM(Vname(1,idfrac(i))), thisindx(ng)
            exit_flag=3
            RETURN
          END IF
        END IF
      END DO
!
!  Write out bed layer sediment properties.
!
      DO i=1,MBEDP
        IF (Hout(idSbed(i),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*b3dvar
          status=nf_fwrite(ng, nchisid(ng), hisVid(idSbed(i),ng),       &
     &                     thisindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, Nbed, scale,          &
# ifdef MASKING
     &                     GRID(ng) % rmask(LBi,LBj),                   &
# endif
     &                     OCEAN(ng) % bed(LBi,LBj,1,i))
          IF (OutThread.and.(status.ne.nf_noerr)) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbed(i))), thisindx(ng)
            exit_flag=3
            RETURN
          END IF
        END IF
      END DO
!
!  Write out bottom sediment properties.
!
      DO i=1,MBOTP
        IF (Hout(idBott(i),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*r2dvar
          status=nf_fwrite(ng, nchisid(ng), hisVid(idBott(i),ng),       &
     &                     thisindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, 1, scale,             &
# ifdef MASKING
     &                     GRID(ng) % rmask(LBi,LBj),                   &
# endif
     &                     OCEAN(ng) % bottom(LBi,LBj,i))
          IF (OutThread.and.(status.ne.nf_noerr)) THEN
            WRITE (stdout,10) TRIM(Vname(1,idBott(i))), thisindx(ng)
            exit_flag=3
            RETURN
          END IF
        END IF
      END DO
#endif
!
!  Synchronize history NetCDF file to disk to allow other processes
!  to access data immediately after it is written.
!
      IF (OutThread) THEN
        status=nf_sync(nchisid(ng))
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,20)
          exit_flag=3
          RETURN
        END IF
        WRITE (stdout,30) thisindx(ng)
      END IF
!
  10  FORMAT (/,' WRT_HIS - error while writing variable: ',a,/,11x,    &
     &        'into history NetCDF file for time record: ',i4)
  20  FORMAT (/,' WRT_HIS - unable to synchronize history NetCDF to ',  &
     &        'disk.')
  30  FORMAT (6x,'WRT_HIS - wrote history  fields into time record = ', &
     &        i10.10)
      RETURN
      END SUBROUTINE wrt_his
