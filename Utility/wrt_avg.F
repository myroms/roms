#include "cppdefs.h"
#ifdef AVERAGES
      SUBROUTINE wrt_avg (ng)
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
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
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
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
      LBi=lbound(GRID(ng)%h,DIM=1)
      UBi=ubound(GRID(ng)%h,DIM=1)
      LBj=lbound(GRID(ng)%h,DIM=2)
      UBj=ubound(GRID(ng)%h,DIM=2)
!
!-----------------------------------------------------------------------
!  Write out time-averaged fields when appropriate.
!-----------------------------------------------------------------------
!
      if (exit_flag.ne.0) RETURN
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
      tavgindx(ng)=tavgindx(ng)+1
      nrecavg(ng)=nrecavg(ng)+1
!
!  Write out averaged time.
!
      IF (OutThread) THEN
        status=nf_put_var1_TYPE(ncavgid(ng), avgVid(idtime,ng),         &
     &                          tavgindx(ng), avgtime(ng))
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,10) TRIM(Vname(1,idtime)), tavgindx(ng)
          exit_flag=3
          RETURN
        END IF
      END IF
!
! Write out free-surface (m).
!
      IF (Hout(idFsur,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite(ng, ncavgid(ng), avgVid(idFsur,ng),            &
     &                   tavgindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
# ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
# endif
     &                   AVERAGE(ng) % avgzeta(LBi,LBj))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idFsur)), tavgindx(ng)
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
        status=nf_fwrite(ng, ncavgid(ng), avgVid(idUbar,ng),            &
     &                   tavgindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
# ifdef MASKING
     &                   GRID(ng) % umask(LBi,LBj),                     &
# endif
     &                   AVERAGE(ng) % avgu2d(LBi,LBj))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUbar)), tavgindx(ng)
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
        status=nf_fwrite(ng, ncavgid(ng), avgVid(idVbar,ng),            &
     &                   tavgindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
# ifdef MASKING
     &                   GRID(ng) % vmask(LBi,LBj),                     &
# endif
     &                   AVERAGE(ng) % avgv2d(LBi,LBj))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVbar)), tavgindx(ng)
          exit_flag=3
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
        status=nf_fwrite(ng, ncavgid(ng), avgVid(idUvel,ng),            &
     &                   tavgindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#  ifdef MASKING
     &                   GRID(ng) % umask(LBi,LBj),                     &
#  endif
     &                   AVERAGE(ng) % avgu3d(LBi,LBj,1))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUvel)), tavgindx(ng)
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
        status=nf_fwrite(ng, ncavgid(ng), avgVid(idVvel,ng),            &
     &                   tavgindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#  ifdef MASKING
     &                   GRID(ng) % vmask(LBi,LBj),                     &
#  endif
     &                   AVERAGE(ng) % avgv3d(LBi,LBj,1))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVvel)), tavgindx(ng)
          exit_flag=3
          RETURN
        END IF
      END IF
!
!  Write out S-coordinate omega vertical velocity (m/s).
!
      IF (Hout(idOvel,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r3dvar
        status=nf_fwrite(ng, ncavgid(ng), avgVid(idOvel,ng),            &
     &                   tavgindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 0, N(ng), scale,           &
#  ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#  endif
     &                   AVERAGE(ng) % avgw3d(LBi,LBj,0))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idOvel)), tavgindx(ng)
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
          status=nf_fwrite(ng, ncavgid(ng), avgTid(itrc,ng),            &
     &                     tavgindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask(LBi,LBj),                   &
#  endif
     &                     AVERAGE(ng) % avgt(LBi,LBj,1,itrc))
          IF (OutThread.and.(status.ne.nf_noerr)) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTvar(itrc))),              &
     &                        tavgindx(ng)
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
        status=nf_fwrite(ng, ncavgid(ng), avgVid(idDano,ng),            &
     &                   tavgindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#  ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#  endif
     &                   AVERAGE(ng) % avgrho(LBi,LBj,1))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idDano)), tavgindx(ng)
          exit_flag=3
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
        status=nf_fwrite(ng, ncavgid(ng), avgVid(idHsbl,ng),            &
     &                   tavgindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
#   ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#   endif
     &                   AVERAGE(ng) % avghsbl(LBi,LBj))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idHsbl)), tavgindx(ng)
          exit_flag=3
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
        status=nf_fwrite(ng, ncavgid(ng), avgVid(idHbbl,ng),            &
     &                   tavgindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
#   ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#   endif
     &                   AVERAGE(ng) % avghbbl(LBi,LBj))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idHbbl)), tavgindx(ng)
          exit_flag=3
          RETURN
        END IF
      END IF
#  endif
#  ifdef AVERAGES_AKV
!
!  Write out vertical viscosity coefficient.
!
      scale=1.0_r8
      gtype=gfactor*w3dvar
      status=nf_fwrite(ng, ncavgid(ng), avgVid(idVvis,ng),              &
     &                 tavgindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, 0, N(ng), scale,             &
#   ifdef MASKING
     &                 GRID(ng) % rmask(LBi,LBj),                       &
#   endif
     &                 AVERAGE(ng) % avgAKv(LBi,LBj,0))
      IF (OutThread.and.(status.ne.nf_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idVvis)), tavgindx(ng)
        exit_flag=3
        RETURN
      END IF
#  endif
#  ifdef AVERAGES_AKT
!
!  Write out vertical diffusion coefficient for potential temperature.
!
      scale=1.0_r8
      gtype=gfactor*w3dvar
      status=nf_fwrite(ng, ncavgid(ng), avgVid(idTdif,ng),              &
     &                 tavgindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, 0, N(ng), scale,             &
#   ifdef MASKING
     &                 GRID(ng) % rmask(LBi,LBj),                       &
#   endif
     &                 AVERAGE(ng) % avgAKt(LBi,LBj,0))
      IF (OutThread.and.(status.ne.nf_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idTdif)), tavgindx(ng)
        exit_flag=3
        RETURN
      END IF
#  endif
#  ifdef AVERAGES_AKS
!
!  Write out vertical diffusion coefficient for salinity.
!
      scale=1.0_r8
      gtype=gfactor*w3dvar
      status=nf_fwrite(ng, ncavgid(ng), avgVid(idSdif,ng),              &
     &                 tavgindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, 0, N(ng), scale,             &
#   ifdef MASKING
     &                 GRID(ng) % rmask(LBi,LBj),                       &
#   endif
     &                 AVERAGE(ng) % avgAKs(LBi,LBj,0))
      IF (OutThread.and.(status.ne.nf_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idSdif)), tavgindx(ng)
        exit_flag=3
        RETURN
      END IF
#  endif
# endif
!
!  Synchronize time-average NetCDF file to disk to allow other processes
!  to access data immediately after it is written.
!
      IF (OutThread) THEN
        status=nf_sync(ncavgid(ng))
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,20)
          exit_flag=3
          RETURN
        END IF
        WRITE (stdout,30) tavgindx(ng)
      END IF
!
  10  FORMAT (/,' WRT_AVG - error while writing variable: ',a,/,11x,    &
     &        'into averages NetCDF file for time record: ',i4)
  20  FORMAT (/,' WRT_AVG - unable to synchronize averages NetCDF to ', &
     &        'disk.')
  30  FORMAT (6x,'WRT_AVG - wrote averaged fields into time record = ', &
     &        i10.10)
#else
      SUBROUTINE wrt_avg
#endif
      RETURN
      END SUBROUTINE wrt_avg
      
