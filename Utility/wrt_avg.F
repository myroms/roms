#include "cppdefs.h"
#ifdef AVERAGES
      SUBROUTINE wrt_avg (ng)
!
!=======================================================================
!  Copyright (c) 2005 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This subroutine writes model time-averaged fields into averages     !
!  NetCDF file.                                                        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_average
      USE mod_forces
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
# if defined SEDIMENT || defined BBL_MODEL
      USE mod_sediment
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
      integer :: LBi, UBi, LBj, UBj
      integer :: gfactor, gtype, itrc, status

      integer :: nf_fwrite

      real(r8) :: scale
!
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)
!
!-----------------------------------------------------------------------
!  Write out time-averaged fields when appropriate.
!-----------------------------------------------------------------------
!
      if (exit_flag.ne.NoError) RETURN
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
!  Set time record index.
!
      tAVGindx(ng)=tAVGindx(ng)+1
      NrecAVG(ng)=NrecAVG(ng)+1
!
!  Write out averaged time.
!
      IF (OutThread) THEN
        status=nf_put_var1_TYPE(ncAVGid(ng), avgVid(idtime,ng),         &
     &                          tAVGindx(ng), AVGtime(ng))
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,10) TRIM(Vname(1,idtime)), tAVGindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out free-surface (m).
!
      IF (Hout(idFsur,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite(ng, iNLM, ncAVGid(ng), avgVid(idFsur,ng),      &
     &                   tAVGindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
# ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
# endif
     &                   AVERAGE(ng) % avgzeta(LBi,LBj))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idFsur)), tAVGindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D momentum component (m/s) in the XI-direction.
!
      IF (Hout(idUbar,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u2dvar
        status=nf_fwrite(ng, iNLM, ncAVGid(ng), avgVid(idUbar,ng),      &
     &                   tAVGindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
# ifdef MASKING
     &                   GRID(ng) % umask(LBi,LBj),                     &
# endif
     &                   AVERAGE(ng) % avgu2d(LBi,LBj))

        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUbar)), tAVGindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D momentum component (m/s) in the ETA-direction.
!
      IF (Hout(idVbar,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v2dvar
        status=nf_fwrite(ng, iNLM, ncAVGid(ng), avgVid(idVbar,ng),      &
     &                   tAVGindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
# ifdef MASKING
     &                   GRID(ng) % vmask(LBi,LBj),                     &
# endif
     &                   AVERAGE(ng) % avgv2d(LBi,LBj))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVbar)), tAVGindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# ifdef SOLVE3D
!
!  Write out 3D momentum component (m/s) in the XI-direction.
!
      IF (Hout(idUvel,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u3dvar
        status=nf_fwrite(ng, iNLM, ncAVGid(ng), avgVid(idUvel,ng),      &
     &                   tAVGindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#  ifdef MASKING
     &                   GRID(ng) % umask(LBi,LBj),                     &
#  endif
     &                   AVERAGE(ng) % avgu3d(LBi,LBj,1))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUvel)), tAVGindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D momentum component (m/s) in the ETA-direction.
!
      IF (Hout(idVvel,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v3dvar
        status=nf_fwrite(ng, iNLM, ncAVGid(ng), avgVid(idVvel,ng),      &
     &                   tAVGindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#  ifdef MASKING
     &                   GRID(ng) % vmask(LBi,LBj),                     &
#  endif
     &                   AVERAGE(ng) % avgv3d(LBi,LBj,1))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVvel)), tAVGindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out S-coordinate omega vertical velocity (m/s).
!
      IF (Hout(idOvel,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*w3dvar
        status=nf_fwrite(ng, iNLM, ncAVGid(ng), avgVid(idOvel,ng),      &
     &                   tAVGindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 0, N(ng), scale,           &
#  ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#  endif
     &                   AVERAGE(ng) % avgw3d(LBi,LBj,0))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idOvel)), tAVGindx(ng)
          exit_flag=3
          ioerror=status
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
          status=nf_fwrite(ng, iNLM, ncAVGid(ng), avgTid(itrc,ng),      &
     &                     tAVGindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask(LBi,LBj),                   &
#  endif
     &                     AVERAGE(ng) % avgt(LBi,LBj,1,itrc))
          IF (OutThread.and.(status.ne.nf_noerr)) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTvar(itrc))),              &
     &                        tAVGindx(ng)
            exit_flag=3
            ioerror=status
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
        status=nf_fwrite(ng, iNLM, ncAVGid(ng), avgVid(idDano,ng),      &
     &                   tAVGindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#  ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#  endif
     &                   AVERAGE(ng) % avgrho(LBi,LBj,1))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idDano)), tAVGindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  ifdef LMD_SKPP
!
!  Write out depth of surface boundary layer.
!
      IF (Hout(idHsbl,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite(ng, iNLM, ncAVGid(ng), avgVid(idHsbl,ng),      &
     &                   tAVGindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
#   ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#   endif
     &                   AVERAGE(ng) % avghsbl(LBi,LBj))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idHsbl)), tAVGindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
#  ifdef LMD_BKPP
!
!  Write out depth of bottom boundary layer.
!
      IF (Hout(idHbbl,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite(ng, iNLM, ncAVGid(ng), avgVid(idHbbl,ng),      &
     &                   tAVGindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
#   ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#   endif
     &                   AVERAGE(ng) % avghbbl(LBi,LBj))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idHbbl)), tAVGindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
# endif
# ifdef AVERAGES_QUADRATIC
!
!  Write out quadratic <zeta*zeta> term.
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf_fwrite(ng, iNLM, ncAVGid(ng), avgVid(idZZav,ng),        &
     &                 tAVGindx(ng), gtype,                           &
     &                 LBi, UBi, LBj, UBj, 1, 1, scale,               &
#  ifdef MASKING
     &                 GRID(ng) % rmask(LBi,LBj),                     &
#  endif
     &                 AVERAGE(ng) % avgZZ(LBi,LBj))
      IF (OutThread.and.(status.ne.nf_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idZZav)), tAVGindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out quadratic <ubar*ubar> term.
!
      scale=1.0_r8
      gtype=gfactor*u2dvar
      status=nf_fwrite(ng, iNLM, ncAVGid(ng), avgVid(idU2av,ng),        &
     &                 tAVGindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, 1, 1, scale,                 &
#  ifdef MASKING
     &                 GRID(ng) % umask(LBi,LBj),                       &
#  endif
     &                 AVERAGE(ng) % avgU2(LBi,LBj))
      IF (OutThread.and.(status.ne.nf_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idU2av)), tAVGindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out quadratic <vbar*vbar> term.
!
      scale=1.0_r8
      gtype=gfactor*v2dvar
      status=nf_fwrite(ng, iNLM, ncAVGid(ng), avgVid(idV2av,ng),        &
     &                 tAVGindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, 1, 1, scale,                 &
#  ifdef MASKING
     &                 GRID(ng) % vmask(LBi,LBj),                       &
#  endif
     &                 AVERAGE(ng) % avgV2(LBi,LBj))
      IF (OutThread.and.(status.ne.nf_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idV2av)), tAVGindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
# endif
# ifdef SOLVE3D
#  ifdef AVERAGES_QUADRATIC
!
!  Write out quadratic <u*u> term.
!
      scale=1.0_r8
      gtype=gfactor*u3dvar
      status=nf_fwrite(ng, iNLM, ncAVGid(ng), avgVid(idUUav,ng),        &
     &                 tAVGindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, 1, N(ng), scale,             &
#   ifdef MASKING
     &                 GRID(ng) % umask(LBi,LBj),                       &
#   endif
     &                 AVERAGE(ng) % avgUU(LBi,LBj,1))
      IF (OutThread.and.(status.ne.nf_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idUUav)), tAVGindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out quadratic <u*v> term.
!
      scale=1.0_r8
      gtype=gfactor*r3dvar
      status=nf_fwrite(ng, iNLM, ncAVGid(ng), avgVid(idUVav,ng),        &
     &                 tAVGindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, 1, N(ng), scale,             &
#   ifdef MASKING
     &                 GRID(ng) % rmask(LBi,LBj),                       &
#   endif
     &                 AVERAGE(ng) % avgUV(LBi,LBj,1))
      IF (OutThread.and.(status.ne.nf_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idUVav)), tAVGindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out quadratic <v*v> term.
!
      scale=1.0_r8
      gtype=gfactor*v3dvar
      status=nf_fwrite(ng, iNLM, ncAVGid(ng), avgVid(idVVav,ng),        &
     &                 tAVGindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, 1, N(ng), scale,             &
#   ifdef MASKING
     &                 GRID(ng) % vmask(LBi,LBj),                       &
#   endif
     &                 AVERAGE(ng) % avgVV(LBi,LBj,1))
      IF (OutThread.and.(status.ne.nf_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idVVav)), tAVGindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out quadratic <t*t> term.
!
      DO itrc=1,NAT
        scale=1.0_r8
        gtype=gfactor*r3dvar
        status=nf_fwrite(ng, iNLM, ncAVGid(ng), avgVid(idTTav(itrc),ng),&
     &                   tAVGindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#   ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#   endif
     &                   AVERAGE(ng) % avgTT(LBi,LBj,1,itrc))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idTTav(itrc))),                &
     &                      tAVGindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END DO
!
!  Write out quadratic <u*t> and <v*t> terms.
!
      DO itrc=1,NAT
        scale=1.0_r8
        gtype=gfactor*u3dvar
        status=nf_fwrite(ng, iNLM, ncAVGid(ng), avgVid(idUTav(itrc),ng),&
     &                   tAVGindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#   ifdef MASKING
     &                   GRID(ng) % umask(LBi,LBj),                     &
#   endif
     &                   AVERAGE(ng) % avgUT(LBi,LBj,1,itrc))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUTav(itrc))),                &
     &                      tAVGindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
        scale=1.0_r8
        gtype=gfactor*v3dvar
        status=nf_fwrite(ng, iNLM, ncAVGid(ng), avgVid(idVTav(itrc),ng),&
     &                   tAVGindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#   ifdef MASKING
     &                   GRID(ng) % vmask(LBi,LBj),                     &
#   endif
     &                   AVERAGE(ng) % avgVT(LBi,LBj,1,itrc))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVTav(itrc))),                &
     &                      tAVGindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END DO
#  endif
#  ifdef AVERAGES_AKV
!
!  Write out vertical viscosity coefficient.
!
      scale=1.0_r8
      gtype=gfactor*w3dvar
      status=nf_fwrite(ng, iNLM, ncAVGid(ng), avgVid(idVvis,ng),        &
     &                 tAVGindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, 0, N(ng), scale,             &
#   ifdef MASKING
     &                 GRID(ng) % rmask(LBi,LBj),                       &
#   endif
     &                 AVERAGE(ng) % avgAKv(LBi,LBj,0))
      IF (OutThread.and.(status.ne.nf_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idVvis)), tAVGindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
#  endif
#  ifdef AVERAGES_AKT
!
!  Write out vertical diffusion coefficient for potential temperature.
!
      scale=1.0_r8
      gtype=gfactor*w3dvar
      status=nf_fwrite(ng, iNLM, ncAVGid(ng), avgVid(idTdif,ng),        &
     &                 tAVGindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, 0, N(ng), scale,             &
#   ifdef MASKING
     &                 GRID(ng) % rmask(LBi,LBj),                       &
#   endif
     &                 AVERAGE(ng) % avgAKt(LBi,LBj,0))
      IF (OutThread.and.(status.ne.nf_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idTdif)), tAVGindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
#  endif
#  ifdef AVERAGES_AKS
!
!  Write out vertical diffusion coefficient for salinity.
!
      scale=1.0_r8
      gtype=gfactor*w3dvar
      status=nf_fwrite(ng, iNLM, ncAVGid(ng), avgVid(idSdif,ng),        &
     &                 tAVGindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, 0, N(ng), scale,             &
#   ifdef MASKING
     &                 GRID(ng) % rmask(LBi,LBj),                       &
#   endif
     &                 AVERAGE(ng) % avgAKs(LBi,LBj,0))
      IF (OutThread.and.(status.ne.nf_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idSdif)), tAVGindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
#  endif
#  ifdef AVERAGES_FLUXES
!
!  Write out surface net heat flux.
!
      IF (Hout(idTsur(itemp),ng)) THEN
        scale=rho0*Cp
        gtype=gfactor*r2dvar
        status=nf_fwrite(ng, iNLM, ncAVGid(ng),                         &
     &                   avgVid(idTsur(itemp),ng),                      &
     &                   tAVGindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
#   ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#   endif
     &                   AVERAGE(ng) % avgstf(LBi,LBj))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idTsur(itemp))), tAVGindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface salt flux  (PSU m/s = kg salt/m2/s).
!
      IF (Hout(idTsur(isalt),ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite(ng, iNLM, ncAVGid(ng),                         &
     &                   avgVid(idTsur(isalt),ng),                      &
     &                   tAVGindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
#   ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#   endif
     &                   AVERAGE(ng) % avgswf(LBi,LBj))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idTsur(isalt))), tAVGindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#   ifdef BULK_FLUXES
!
!  Write out latent heat flux.
!
      IF (Hout(idLhea,ng)) THEN
        scale=rho0*Cp
        gtype=gfactor*r2dvar
        status=nf_fwrite(ng, iNLM, ncAVGid(ng), avgVid(idLhea,ng),      &
     &                   tAVGindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
#    ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#    endif
     &                   AVERAGE(ng) % avglhf(LBi,LBj))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idLhea)), tAVGindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out sensible heat flux.
!
      IF (Hout(idShea,ng)) THEN
        scale=rho0*Cp
        gtype=gfactor*r2dvar
        status=nf_fwrite(ng, iNLM, ncAVGid(ng), avgVid(idShea,ng),      &
     &                   tAVGindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
#    ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#    endif
     &                   AVERAGE(ng) % avgshf(LBi,LBj))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idShea)), tAVGindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out longwave radiation flux.
!
      IF (Hout(idLrad,ng)) THEN
        scale=rho0*Cp
        gtype=gfactor*r2dvar
        status=nf_fwrite(ng, iNLM, ncAVGid(ng), avgVid(idLrad,ng),      &
     &                   tAVGindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
#    ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#    endif
     &                   AVERAGE(ng) % avglrf(LBi,LBj))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idLrad)), tAVGindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#    ifdef EMINUSP
!
!  Write out evaportaion rate (kg/m2/s).
!
      IF (Hout(idevap,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite(ng, iNLM, ncAVGid(ng), avgVid(idevap,ng),      &
     &                   tAVGindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
#     ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#     endif
     &                   AVERAGE(ng) % avgevap(LBi,LBj))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idevap)), tAVGindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out precipitation rate (kg/m2/s).
!
      IF (Hout(idrain,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite(ng, iNLM, ncAVGid(ng), avgVid(idrain,ng),      &
     &                   tAVGindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
#     ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#     endif
     &                   AVERAGE(ng) % avgrain(LBi,LBj))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idrain)), tAVGindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#    endif
#   endif
#   ifdef SHORTWAVE
!
!  Write out shortwave radiation flux.
!
      IF (Hout(idSrad,ng)) THEN
        scale=rho0*Cp
        gtype=gfactor*r2dvar
        status=nf_fwrite(ng, iNLM, ncAVGid(ng), avgVid(idSrad,ng),      &
     &                   tAVGindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
#    ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#    endif
     &                   AVERAGE(ng) % avgsrf(LBi,LBj))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idSrad)), tAVGindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#   endif
#  endif
# endif
# ifdef AVERAGES_FLUXES
!
!  Write out surface u-momentun stress.
!
      IF (Hout(idUsms,ng)) THEN
        scale=rho0
        gtype=gfactor*u2dvar
        status=nf_fwrite(ng, iNLM, ncAVGid(ng), avgVid(idUsms,ng),      &
     &                   tAVGindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
#  ifdef MASKING
     &                   GRID(ng) % umask(LBi,LBj),                     &
#  endif
     &                   AVERAGE(ng) % avgsus(LBi,LBj))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUsms)), tAVGindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface v-momentun stress.
!
      IF (Hout(idVsms,ng)) THEN
        scale=rho0
        gtype=gfactor*v2dvar
        status=nf_fwrite(ng, iNLM, ncAVGid(ng), avgVid(idVsms,ng),      &
     &                   tAVGindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
#  ifdef MASKING
     &                   GRID(ng) % vmask(LBi,LBj),                     &
#  endif
     &                   AVERAGE(ng) % avgsvs(LBi,LBj))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVsms)), tAVGindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
# if defined SEDIMENT && defined BEDLOAD
!
!  Write out bed load transport in U-direction.
!
      DO itrc=1,NST
        IF (Hout(idUbld(itrc),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*u2dvar
          status=nf_fwrite(ng, iNLM, ncAVGid(ng),                       &
     &                     avgVid(idUbld(itrc),ng),                     &
     &                     tAVGindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, 1,  scale,            &
#  ifdef MASKING
     &                     GRID(ng) % umask(LBi,LBj),                   &
#  endif
     &                     AVERAGE(ng) % avgbedldu(LBi,LBj,itrc))
          IF (OutThread.and.(status.ne.nf_noerr)) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbld(itrc))), tAVGindx(ng)
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
!
!  Write out bed load transport in V-direction.
!
        IF (Hout(idVbld(itrc),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*v2dvar
          status=nf_fwrite(ng, iNLM, ncAVGid(ng),                       &
     &                     avgVid(idVbld(itrc),ng),                     &
     &                     tAVGindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, 1,  scale,            &
#  ifdef MASKING
     &                     GRID(ng) % vmask(LBi,LBj),                   &
#  endif
     &                     AVERAGE(ng) % avgbedldv(LBi,LBj,itrc))
          IF (OutThread.and.(status.ne.nf_noerr)) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbld(itrc))), tAVGindx(ng)
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
# endif
!
!  Synchronize time-average NetCDF file to disk to allow other processes
!  to access data immediately after it is written.
!
      IF (OutThread) THEN
        status=nf_sync(ncAVGid(ng))
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,20)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        WRITE (stdout,30) tAVGindx(ng)
      END IF
!
  10  FORMAT (/,' WRT_AVG - error while writing variable: ',a,/,11x,    &
     &        'into averages NetCDF file for time record: ',i4)
  20  FORMAT (/,' WRT_AVG - unable to synchronize averages NetCDF to ', &
     &        'disk.')
  30  FORMAT (6x,'WRT_AVG   - wrote averaged fields into time ',        &
     &        'record =',t72,i7.7)
#else
      SUBROUTINE wrt_avg
#endif
      RETURN
      END SUBROUTINE wrt_avg
      
