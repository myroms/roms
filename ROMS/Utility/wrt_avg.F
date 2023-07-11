#include "cppdefs.h"
      MODULE wrt_avg_mod
#if defined AVERAGES    || \
   (defined AD_AVERAGES && defined ADJOINT) || \
   (defined RP_AVERAGES && defined TL_IOMS) || \
   (defined TL_AVERAGES && defined TANGENT)
!
!git $Id$
!svn $Id: wrt_avg.F 1178 2023-07-11 17:50:57Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module writes model time-averaged fields into output file      !
!  using the standard NetCDF library or the Parallel-IO (PIO) library. !
!                                                                      !
!  Notice that only momentum is affected by the full time-averaged     !
!  masks.  If applicable, these mask contains information about        !
!  river runoff and time-dependent wetting and drying variations.      !
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
      USE mod_sedbed
      USE mod_sediment
# endif
# if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
      USE mod_tides
# endif
!
#ifdef ICE_MODEL
      USE ice_output_mod,  ONLY : ice_wrt_nf90
# if defined PIO_LIB && defined DISTRIBUTE
      USE ice_output_mod,  ONLY : ice_wrt_pio
# endif
#endif
      USE nf_fwrite2d_mod, ONLY : nf_fwrite2d
# ifdef SOLVE3D
      USE nf_fwrite3d_mod, ONLY : nf_fwrite3d
# endif
      USE strings_mod,     ONLY : FoundError
# if defined WEC || defined WEC_VF
      USE wec_output_mod,  ONLY : wec_wrt_nf90
# if defined PIO_LIB && defined DISTRIBUTE
      USE wec_output_mod,  ONLY : wec_wrt_pio
# endif
# endif
!
      implicit none
!
      PUBLIC  :: wrt_avg
      PRIVATE :: wrt_avg_nf90
# if defined PIO_LIB && defined DISTRIBUTE
      PRIVATE :: wrt_avg_pio
# endif
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE wrt_avg (ng, tile)
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
!  Write out time-averaged fields according to IO type.
!-----------------------------------------------------------------------
!
      LBi=BOUNDS(ng)%LBi(tile)
      UBi=BOUNDS(ng)%UBi(tile)
      LBj=BOUNDS(ng)%LBj(tile)
      UBj=BOUNDS(ng)%UBj(tile)
!
      SELECT CASE (AVG(ng)%IOtype)
        CASE (io_nf90)
          CALL wrt_avg_nf90 (ng, iNLM, tile,                            &
     &                       LBi, UBi, LBj, UBj)

# if defined PIO_LIB && defined DISTRIBUTE
        CASE (io_pio)
          CALL wrt_avg_pio (ng, iNLM, tile,                             &
     &                      LBi, UBi, LBj, UBj)
# endif
        CASE DEFAULT
          IF (Master) WRITE (stdout,10) AVG(ng)%IOtype
          exit_flag=3
      END SELECT
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
  10  FORMAT (' WRT_AVG - Illegal output file type, io_type = ',i0,     &
     &        /,15x,'Check KeyWord ''OUT_LIB'' in ''roms.in''.')
!
      RETURN
      END SUBROUTINE wrt_avg
!
!***********************************************************************
      SUBROUTINE wrt_avg_nf90 (ng, model, tile,                         &
     &                         LBi, UBi, LBj, UBj)
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
      integer :: Fcount, gfactor, gtype, i, itrc, status
!
      real(dp) :: scale
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", wrt_avg_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out time-averaged fields when appropriate.
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
!  Set time record index.
!
      AVG(ng)%Rindex=AVG(ng)%Rindex+1
      Fcount=AVG(ng)%load
      AVG(ng)%Nrec(Fcount)=AVG(ng)%Nrec(Fcount)+1
!
!  Write out averaged time.
!
      CALL netcdf_put_fvar (ng, model, AVG(ng)%name,                    &
     &                      TRIM(Vname(1,idtime)), AVGtime(ng:),        &
     &                      (/AVG(ng)%Rindex/), (/1/),                  &
     &                      ncid = AVG(ng)%ncid,                        &
     &                      varid = AVG(ng)%Vid(idtime))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out free-surface (m).
!
      IF (Aout(idFsur,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idFsur),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     AVERAGE(ng) % avgzeta)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idFsur)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

# if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
!
!  Write out detided free-surface (m).
!
      IF (Aout(idFsuD,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idFsuD),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     TIDES(ng) % zeta_detided)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idFsuD)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
!
!  Write out 2D momentum component (m/s) in the XI-direction.
!
      IF (Aout(idUbar,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idUbar),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % umask_avg,                        &
# endif
     &                     AVERAGE(ng) % avgu2d)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbar)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

# if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
!
!  Write out detided 2D momentum component (m/s) in the XI-direction.
!
      IF (Aout(idu2dD,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idu2dD),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % umask_avg,                        &
#  endif
     &                     TIDES(ng) % ubar_detided)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idu2dD)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
!
!  Write out 2D momentum component (m/s) in the ETA-direction.
!
      IF (Aout(idVbar,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idVbar),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % vmask_avg,                        &
# endif
     &                     AVERAGE(ng) % avgv2d)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbar)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

# if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
!
!  Write out detided 2D momentum component (m/s) in the ETA-direction.
!
      IF (Aout(idv2dD,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idv2dD),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % vmask_avg,                        &
#  endif
     &                     TIDES(ng) % vbar_detided)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idv2dD)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
!
!  Write out 2D Eastward momentum component (m/s) at RHO-points.
!
      IF (Aout(idu2dE,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idu2dE),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask_avg,                        &
# endif
     &                     AVERAGE(ng) % avgu2dE)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idu2dE)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D Northward momentum component (m/s) at RHO-points.
!
      IF (Aout(idv2dN,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idv2dN),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask_avg,                        &
# endif
     &                     AVERAGE(ng) % avgv2dN)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idv2dN)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

# if defined FORWARD_WRITE && defined SOLVE3D
!
!  Write out time-averaged mass fluxes for 3D momentum coupling.
!
      IF (Aout(idUfx1,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idUfx1),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % umask_avg,                        &
#  endif
     &                     AVERAGE(ng) % avgDU_avg1)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUfx1)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
      IF (Aout(idUfx2,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idUfx2),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % umask_avg,                        &
#  endif
     &                     AVERAGE(ng) % avgDU_avg2)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUfx2)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
      IF (Aout(idVfx1,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idVfx1),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % vmask_avg,                        &
#  endif
     &                     AVERAGE(ng) % avgDV_avg1)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVfx1)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
      IF (Aout(idVfx2,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idVfx2),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % umask_avg,                        &
#  endif
     &                     AVERAGE(ng) % avgDV_avg2)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVfx2)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif

# ifdef SOLVE3D
!
!  Write out 3D momentum component (m/s) in the XI-direction.
!
      IF (Aout(idUvel,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*u3dvar
        status=nf_fwrite3d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idUvel),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % umask_avg,                        &
#  endif
     &                     AVERAGE(ng) % avgu3d)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUvel)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

#  if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
!
!  Write out detided 3D momentum component (m/s) in the XI-direction.
!
      IF (Aout(idu3dD,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*u3dvar
        status=nf_fwrite3d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idu3dD),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % umask_avg,                        &
#   endif
     &                     TIDES(ng) % u_detided)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idu3dD)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
!
!  Write out 3D momentum component (m/s) in the ETA-direction.
!
      IF (Aout(idVvel,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*v3dvar
        status=nf_fwrite3d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idVvel),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % vmask_avg,                        &
#  endif
     &                     AVERAGE(ng) % avgv3d)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVvel)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

#  if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
!
!  Write out detided 3D momentum component (m/s) in the ETA-direction.
!
      IF (Aout(idv3dD,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*v3dvar
        status=nf_fwrite3d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idv3dD),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % vmask_avg,                        &
#   endif
     &                     TIDES(ng) % v_detided)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idv3dD)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
!
!  Write out 3D Eastward momentum component (m/s) at RHO-points.
!
      IF (Aout(idu3dE,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idu3dE),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask_avg,                        &
#  endif
     &                     AVERAGE(ng) % avgu3dE)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idu3dE)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D Northward momentum component (m/s) at RHO-points.
!
      IF (Aout(idv3dN,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idv3dN),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask_avg,                        &
#  endif
     &                     AVERAGE(ng) % avgv3dN)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idv3dN)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out S-coordinate omega vertical velocity (m/s).
!
      IF (Aout(idOvel,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idOvel),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     AVERAGE(ng) % avgw3d)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idOvel)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out "true" vertical velocity (m/s).
!
      IF (Aout(idWvel,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idWvel),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     AVERAGE(ng) % avgwvel)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idOvel)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out tracer type variables.
!
      DO itrc=1,NT(ng)
        IF (Aout(idTvar(itrc),ng)) THEN
          scale=1.0_dp
          gtype=gfactor*r3dvar
          status=nf_fwrite3d(ng, model, AVG(ng)%ncid,                   &
     &                       AVG(ng)%Tid(itrc),                         &
     &                       AVG(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       AVERAGE(ng) % avgt(:,:,:,itrc))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idTvar(itrc))),            &
     &                          AVG(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO

#  if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
!
!  Write out detided active tracer variables.
!
      DO itrc=1,NAT
        IF (Aout(idTrcD(itrc),ng)) THEN
          scale=1.0_dp
          gtype=gfactor*r3dvar
          status=nf_fwrite3d(ng, model, AVG(ng)%ncid,                   &
     &                       AVG(ng)%Vid(idTrcD(itrc)),                 &
     &                       AVG(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       TIDES(ng) % t_detided(:,:,:,itrc))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idTrcD(itrc))),            &
     &                          AVG(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
#  endif
!
!  Write out density anomaly.
!
      IF (Aout(idDano,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idDano),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     AVERAGE(ng) % avgrho)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idDano)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

#  ifdef LMD_SKPP
!
!  Write out depth of surface boundary layer.
!
      IF (Aout(idHsbl,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idHsbl),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     AVERAGE(ng) % avghsbl)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idHsbl)), AVG(ng)%Rindex
          END IF
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
      IF (Aout(idHbbl,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idHbbl),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     AVERAGE(ng) % avghbbl)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idHbbl)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
# endif
!
!  Write out 2D potential vorticity.
!
      IF (Aout(id2dPV,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*p2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(id2dPV),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % pmask,                            &
# endif
     &                     AVERAGE(ng) % avgpvor2d)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,id2dPV)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D relative vorticity.
!
      IF (Aout(id2dRV,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*p2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(id2dRV),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % pmask,                            &
# endif
     &                     AVERAGE(ng) % avgrvor2d)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,id2dRV)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

# ifdef SOLVE3D
!
!  Write out 3D potential vorticity.
!
      IF (Aout(id3dPV,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*p3dvar
        status=nf_fwrite3d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(id3dPV),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % pmask,                            &
#  endif
     &                     AVERAGE(ng) % avgpvor3d)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,id3dPV)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D relative vorticity.
!
      IF (Aout(id3dRV,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*p3dvar
        status=nf_fwrite3d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(id3dRV),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % pmask,                            &
#  endif
     &                     AVERAGE(ng) % avgrvor3d)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,id3dRV)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
!
!  Write out quadratic <zeta*zeta> term.
!
      IF (Aout(idZZav,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idZZav),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     AVERAGE(ng) % avgZZ)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idZZav)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out quadratic <ubar*ubar> term.
!
      IF (Aout(idU2av,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idU2av),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % umask,                            &
# endif
     &                     AVERAGE(ng) % avgU2)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idU2av)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out quadratic <vbar*vbar> term.
!
      IF (Aout(idV2av,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idV2av),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % vmask,                            &
# endif
     &                     AVERAGE(ng) % avgV2)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idV2av)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

# ifdef SOLVE3D
!
!  Write out u-volume flux.
!
      IF (Aout(idHUav,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*u3dvar
        status=nf_fwrite3d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idHUav),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % umask,                            &
#  endif
     &                     AVERAGE(ng) % avgHuon)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idHUav)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out v-volume flux.
!
      IF (Aout(idHVav,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*v3dvar
        status=nf_fwrite3d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idHVav),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % vmask,                            &
#  endif
     &                     AVERAGE(ng) % avgHvom)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idHVav)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out quadratic <u*u> term.
!
      IF (Aout(idUUav,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*u3dvar
        status=nf_fwrite3d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idUUav),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % umask,                            &
#  endif
     &                     AVERAGE(ng) % avgUU)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUUav)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out quadratic <u*v> term.
!
      IF (Aout(idUVav,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idUVav),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     AVERAGE(ng) % avgUV)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUVav)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out quadratic <v*v> term.
!
      IF (Aout(idVVav,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*v3dvar
        status=nf_fwrite3d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idVVav),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % vmask,                            &
#  endif
     &                     AVERAGE(ng) % avgVV)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVVav)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out quadratic <t*t> term.
!
      DO i=1,NT(ng)
        IF (Aout(idTTav(i),ng)) THEN
          scale=1.0_dp
          gtype=gfactor*r3dvar
          status=nf_fwrite3d(ng, model, AVG(ng)%ncid,                   &
     &                       AVG(ng)%Vid(idTTav(i)),                    &
     &                       AVG(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       AVERAGE(ng) % avgTT(:,:,:,i))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idTTav(i))), AVG(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
!
!  Write out active tracer volume fluxes.
!
      DO i=1,NT(ng)
        IF (Aout(iHUTav(i),ng)) THEN
          scale=1.0_dp
          gtype=gfactor*u3dvar
          status=nf_fwrite3d(ng, model, AVG(ng)%ncid,                   &
     &                       AVG(ng)%Vid(iHUTav(i)),                    &
     &                       AVG(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#  ifdef MASKING
     &                       GRID(ng) % umask,                          &
#  endif
     &                       AVERAGE(ng) % avgHuonT(:,:,:,i))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,iHUTav(i))), AVG(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
!
        IF (Aout(iHVTav(i),ng)) THEN
          scale=1.0_dp
          gtype=gfactor*v3dvar
          status=nf_fwrite3d(ng, model, AVG(ng)%ncid,                   &
     &                       AVG(ng)%Vid(iHVTav(i)),                    &
     &                       AVG(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#  ifdef MASKING
     &                       GRID(ng) % vmask,                          &
#  endif
     &                       AVERAGE(ng) % avgHvomT(:,:,:,i))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,iHVTav(i))), AVG(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
!
!  Write out quadratic <u*t> and <v*t> terms.
!
      DO i=1,NT(ng)
        IF (Aout(idUTav(i),ng)) THEN
          scale=1.0_dp
          gtype=gfactor*u3dvar
          status=nf_fwrite3d(ng, model, AVG(ng)%ncid,                   &
     &                       AVG(ng)%Vid(idUTav(i)),                    &
     &                       AVG(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#  ifdef MASKING
     &                       GRID(ng) % umask,                          &
#  endif
     &                       AVERAGE(ng) % avgUT(:,:,:,i))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idUTav(i))), AVG(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
!
        IF (Aout(idVTav(i),ng)) THEN
          scale=1.0_dp
          gtype=gfactor*v3dvar
          status=nf_fwrite3d(ng, model, AVG(ng)%ncid,                   &
     &                       AVG(ng)%Vid(idVTav(i)),                    &
     &                       AVG(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#  ifdef MASKING
     &                       GRID(ng) % vmask,                          &
#  endif
     &                       AVERAGE(ng) % avgVT(:,:,:,i))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idVTav(i))), AVG(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO

#  if defined LMD_MIXING || defined MY25_MIXING || defined GLS_MIXING
!
!  Write out vertical viscosity coefficient.
!
      IF (Aout(idVvis,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idVvis),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     AVERAGE(ng) % avgAKv,                        &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVvis)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out vertical diffusion coefficient for potential temperature.
!
      IF (Aout(idTdif,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idTdif),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     AVERAGE(ng) % avgAKt,                        &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTdif)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

#   ifdef SALINITY
!
!  Write out vertical diffusion coefficient for salinity.
!
      IF (Aout(idSdif,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idSdif),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#    ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#    endif
     &                     AVERAGE(ng) % avgAKs,                        &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSdif)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#   endif
#  endif
#  if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
!
!  Write out surface air pressure.
!
      IF (Aout(idPair,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idPair),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     AVERAGE(ng) % avgPair)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idPair)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
#  if defined BULK_FLUXES
!
!  Write out surface air temperature.
!
      IF (Aout(idTair,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idTair),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     AVERAGE(ng) % avgTair)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTair)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
#  if defined BULK_FLUXES || defined ECOSIM
!
!  Write out surface winds.
!
      IF (Aout(idUair,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idUair),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     AVERAGE(ng) % avgUwind)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUair)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
      IF (Aout(idVair,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idVair),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     AVERAGE(ng) % avgVwind)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVair)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out Eastward/Northward surface wind (m/s) at RHO-points.
!
      IF (Aout(idUaiE,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idUaiE),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask_avg,                        &
#   endif
     &                     AVERAGE(ng) % avguwindE)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUaiE)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
      IF (Aout(idVaiN,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idVaiN),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask_avg,                        &
#   endif
     &                     AVERAGE(ng) % avgvwindN)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVaiN)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
!
!  Write out surface net heat flux.
!
      IF (Aout(idTsur(itemp),ng)) THEN
#  ifdef ADJOINT
#   if defined AD_SENSITIVITY   || defined I4DVAR_ANA_SENSITIVITY || \
       defined OPT_OBSERVATIONS
        scale=1.0_dp/(rho0*Cp)
#   else
        scale=1.0_dp
#   endif
#  else
        scale=rho0*Cp
#  endif
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idTsur(itemp)),                  &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     AVERAGE(ng) % avgstf)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTsur(itemp))),             &
     &                        AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

#  ifdef SALINITY
!
!  Write out surface salt flux  (PSU m/s = kg salt/m2/s).
!
      IF (Aout(idTsur(isalt),ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idTsur(isalt)),                  &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     AVERAGE(ng) % avgswf)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTsur(isalt))),             &
     &                        AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif

#  if defined BULK_FLUXES || defined FRC_COUPLING
!
!  Write out latent heat flux.
!
      IF (Aout(idLhea,ng)) THEN
#   ifdef ADJOINT
        scale=1.0_dp
#   else
        scale=rho0*Cp
#   endif
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idLhea),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     AVERAGE(ng) % avglhf)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idLhea)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out sensible heat flux.
!
      IF (Aout(idShea,ng)) THEN
#   ifdef ADJOINT
        scale=1.0_dp
#   else
        scale=rho0*Cp
#   endif
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idShea),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     AVERAGE(ng) % avgshf)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idShea)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out net longwave radiation flux.
!
      IF (Aout(idLrad,ng)) THEN
#   ifdef ADJOINT
        scale=1.0_dp
#   else
        scale=rho0*Cp
#   endif
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idLrad),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     AVERAGE(ng) % avglrf)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idLrad)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif

#  ifdef BULK_FLUXES
#   ifdef EMINUSP
!
!  Write out evaportaion rate (kg/m2/s).
!
      IF (Aout(idevap,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idevap),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#    endif
     &                     AVERAGE(ng) % avgevap)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idevap)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out precipitation rate (kg/m2/s).
!
      IF (Aout(idrain,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idrain),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#    endif
     &                     AVERAGE(ng) % avgrain)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idrain)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#   endif
#  endif
#  ifdef SHORTWAVE
!
!  Write out net shortwave radiation flux.
!
      IF (Aout(idSrad,ng)) THEN
#    ifdef ADJOINT
        scale=1.0_dp
#    else
        scale=rho0*Cp
#    endif
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idSrad),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#    endif
     &                     AVERAGE(ng) % avgsrf)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSrad)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
# endif
!
!  Write out surface u-momentum stress.
!
      IF (Aout(idUsms,ng)) THEN
# ifdef ADJOINT
#  if defined AD_SENSITIVITY   || defined I4DVAR_ANA_SENSITIVITY || \
      defined OPT_OBSERVATIONS
        scale=1.0_dp/rho0
#  else
        scale=1.0_dp
#  endif
# else
        scale=rho0
# endif
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idUsms),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % umask,                            &
# endif
     &                     AVERAGE(ng) % avgsus)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUsms)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface v-momentum stress.
!
      IF (Aout(idVsms,ng)) THEN
# ifdef ADJOINT
#  if defined AD_SENSITIVITY   || defined I4DVAR_ANA_SENSITIVITY || \
      defined OPT_OBSERVATIONS
        scale=1.0_dp/rho0
#  else
        scale=1.0_dp
#  endif
# else
        scale=rho0
# endif
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idVsms),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % vmask,                            &
# endif
     &                     AVERAGE(ng) % avgsvs)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVsms)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom u-momentum stress.
!
      IF (Aout(idUbms,ng)) THEN
# ifdef ADJOINT
#  if defined AD_SENSITIVITY   || defined I4DVAR_ANA_SENSITIVITY || \
      defined OPT_OBSERVATIONS
        scale=1.0_dp/rho0
#  else
        scale=1.0_dp
#  endif
# else
        scale=rho0
# endif
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idUbms),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % umask,                            &
# endif
     &                     AVERAGE(ng) % avgbus)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbms)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom v-momentum stress.
!
      IF (Aout(idVbms,ng)) THEN
# ifdef ADJOINT
#  if defined AD_SENSITIVITY   || defined I4DVAR_ANA_SENSITIVITY || \
      defined OPT_OBSERVATIONS
        scale=1.0_dp/rho0
#  else
        scale=1.0_dp
#  endif
# else
        scale=rho0
# endif
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idVbms),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % vmask,                            &
# endif
     &                     AVERAGE(ng) % avgbvs)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbms)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

# ifdef BBL_MODEL
!
!  Write out bottom momentum U-stress at RHO-!points.
!
      IF (Aout(idUbrs,ng)) THEN
# ifdef ADJOINT
#  if defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
      defined OPT_OBSERVATIONS
        scale=1.0_dp/rho0
#  else
        scale=1.0_dp
#  endif
# else
        scale=rho0
# endif
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idUbrs),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     AVERAGE(ng) % avgUbrs)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbrs)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom momentum V-stress at RHO-points.
!
      IF (Aout(idVbrs,ng)) THEN
# ifdef ADJOINT
#  if defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
      defined OPT_OBSERVATIONS
        scale=1.0_dp/rho0
#  else
        scale=1.0_dp
#  endif
# else
        scale=rho0
# endif
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idVbrs),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     AVERAGE(ng) % avgVbrs)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbrs)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom wave U-stress at RHO-points.
!
      IF (Aout(idUbws,ng)) THEN
# ifdef ADJOINT
#  if defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
      defined OPT_OBSERVATIONS
        scale=1.0_dp/rho0
#  else
        scale=1.0_dp
#  endif
# else
        scale=rho0
# endif
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idUbws),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     AVERAGE(ng) % avgUbws)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbws)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom wave V-stress at RHO-points.
!
      IF (Aout(idVbws,ng)) THEN
# ifdef ADJOINT
#  if defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
      defined OPT_OBSERVATIONS
        scale=1.0_dp/rho0
#  else
        scale=1.0_dp
#  endif
# else
        scale=rho0
# endif
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idVbws),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     AVERAGE(ng) % avgVbws)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbws)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom max wave current U-stress at RHO-points.
!
      IF (Aout(idUbcs,ng)) THEN
# ifdef ADJOINT
#  if defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
      defined OPT_OBSERVATIONS
        scale=1.0_dp/rho0
#  else
        scale=1.0_dp
#  endif
# else
        scale=rho0
# endif
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idUbcs),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     AVERAGE(ng) % avgUbcs)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbcs)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out max bottom wave current V-stress at RHO-points.
!
      IF (Aout(idVbcs,ng)) THEN
# ifdef ADJOINT
#  if defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
      defined OPT_OBSERVATIONS
        scale=1.0_dp/rho0
#  else
        scale=1.0_dp
#  endif
# else
        scale=rho0
# endif
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idVbcs),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     AVERAGE(ng) % avgVbcs)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbcs)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out max bottom wave current stress magnitude at RHO-points.
!
      IF (Aout(idUVwc,ng)) THEN
# ifdef ADJOINT
#  if defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
      defined OPT_OBSERVATIONS
        scale=1.0_dp/rho0
#  else
        scale=1.0_dp
#  endif
# else
        scale=rho0
# endif
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idUVwc),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     AVERAGE(ng) % avgUVwc)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUVwc)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom wave orbital U-velocity at RHO-points.
!
      IF (Aout(idUbot,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idUbot),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     AVERAGE(ng) % avgUbot)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbot)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom wave orbital V-velocity at RHO-points.
!
      IF (Aout(idVbot,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idVbot),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     AVERAGE(ng) % avgVbot)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbot)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom U-current above bed at RHO-points.
!
      IF (Aout(idUbur,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idUbur),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     AVERAGE(ng) % avgUbur)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbur)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom V-current above the bed at RHO-points.
!
      IF (Aout(idVbvr,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idVbvr),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     AVERAGE(ng) % avgVbvr)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbvr)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#endif

# if defined SEDIMENT && defined BEDLOAD
!
!  Write out bed load transport in U-direction.
!
      DO itrc=1,NST
        IF (Aout(idUbld(itrc),ng)) THEN
          scale=1.0_dp
          gtype=gfactor*u2dvar
          status=nf_fwrite2d(ng, model, AVG(ng)%ncid,                   &
     &                       AVG(ng)%Vid(idUbld(itrc)),                 &
     &                       AVG(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % umask,                          &
#  endif
     &                       SEDBED(ng) % avgbedldu(:,:,itrc))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idUbld(itrc))),            &
     &                          AVG(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
!
!  Write out bed load transport in V-direction.
!
        IF (Aout(idVbld(itrc),ng)) THEN
          scale=1.0_dp
          gtype=gfactor*v2dvar
          status=nf_fwrite2d(ng, model, AVG(ng)%ncid,                   &
     &                       AVG(ng)%Vid(idVbld(itrc)),                 &
     &                       AVG(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % vmask,                          &
#  endif
     &                       SEDBED(ng) % avgbedldv(:,:,itrc))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idVbld(itrc))),            &
     &                          AVG(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
# endif
# ifdef WAVES_HEIGHT
!
!  Write out wind-induced significant wave height.
!
      IF (Aout(idWamp,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idWamp),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     AVERAGE(ng) % avgWamp)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWamp)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out wind-induced significant wave height squared.
!
      IF (Aout(idWam2,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idWam2),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     AVERAGE(ng) % avgWam2)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWam2)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
# ifdef WAVES_LENGTH
!
!  Write out wind-induced mean wavelength.
!
      IF (Aout(idWlen,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idWlen),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     AVERAGE(ng) % avgWlen)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWlen)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
# ifdef WAVES_LENGTHP
!
!  Write out wind-induced peak wave wavelength.
!
      IF (Aout(idWlep,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idWlep),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     AVERAGE(ng) % avgWlep)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWlep)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
# ifdef WAVES_DIR
!
!  Write out wind-induced mean wave direction.
!
      IF (Aout(idWdir,ng)) THEN
        scale=rad2deg
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idWdir),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     AVERAGE(ng) % avgWdir)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWdir)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
# ifdef WAVES_DIRP
!
!  Write out wind-induced peak wave direction.
!
      IF (Aout(idWdip,ng)) THEN
        scale=rad2deg
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idWdip),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     AVERAGE(ng) % avgWdip)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWdip)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
# ifdef WAVES_TOP_PERIOD
!
!  Write out wind-induced surface wave period.
!
      IF (Aout(idWptp,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idWptp),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     AVERAGE(ng) % avgWptp)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWptp)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
# ifdef WAVES_BOT_PERIOD
!
!  Write out wind-induced bottom wave period.

      IF (Aout(idWpbt,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idWpbt),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     AVERAGE(ng) % avgWpbt)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWpbt)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
# if defined BBL_MODEL    || defined BEDLOAD_VANDERA || \
     defined WAV_COUPLING
!
!  Write out wind-induced bottom orbital velocity.
!
      IF (Aout(idWorb,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idWorb),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     AVERAGE(ng) % avgWorb)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWorb)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
# ifdef UV_KIRBY
!
!  Write out U-velocity from Kirby and Chen.
!
      IF (Aout(idUwav,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idUwav),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     AVERAGE(ng) % avgUwav)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUwav)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out V-velocity from Kirby and Chen.
!
      IF (Aout(idVwav,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idVwav),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     AVERAGE(ng) % avgVwav)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVwav)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
# ifdef ICE_MODEL
!
!-----------------------------------------------------------------------
!  Write out sea-ice model variables.
!-----------------------------------------------------------------------
!
        CALL ice_wrt_nf90 (ng, model, tile,                             &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     Aout, AVG)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
# if defined WEC || defined WEC_VF
!
!-----------------------------------------------------------------------
!  Write out Waves Effect on Currents variables.
!-----------------------------------------------------------------------
!
      CALL wec_wrt_nf90 (ng, model, tile,                               &
    &                    LBi, UBi, LBj, UBj,                            &
    &                    Aout, AVG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
!
!-----------------------------------------------------------------------
!  Synchronize time-average NetCDF file to disk to allow other processes
!  to access data immediately after it is written.
!-----------------------------------------------------------------------
!
      CALL netcdf_sync (ng, model, AVG(ng)%name, AVG(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#ifdef NESTING
      IF (Master) WRITE (stdout,20) AVG(ng)%Rindex, ng
#else
      IF (Master) WRITE (stdout,20) AVG(ng)%Rindex
#endif
!
  10  FORMAT (/,' WRT_AVG_NF90 - error while writing variable: ',a,     &
     &        /,16x,'into averages NetCDF file for time record: ',i0)
  20  FORMAT (2x,'WRT_AVG_NF90     - wrote averaged',t40,'fields',t59,  &
# ifdef NESTING
     &        'in record = ',i0,t92,i2.2)
# else
     &        'in record = ',i0)
# endif
!
      RETURN
      END SUBROUTINE wrt_avg_nf90

# if defined PIO_LIB && defined DISTRIBUTE
!
!***********************************************************************
      SUBROUTINE wrt_avg_pio (ng, model, tile,                          &
     &                        LBi, UBi, LBj, UBj)
!***********************************************************************
!
      USE mod_pio_netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
!  Local variable declarations.
!
      integer :: Fcount, i, itrc, status
!
      real(dp) :: scale
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", wrt_avg_pio"
!
      TYPE (IO_desc_t), pointer :: ioDesc
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out time-averaged fields when appropriate.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Set time record index.
!
      AVG(ng)%Rindex=AVG(ng)%Rindex+1
      Fcount=AVG(ng)%load
      AVG(ng)%Nrec(Fcount)=AVG(ng)%Nrec(Fcount)+1
!
!  Write out averaged time.
!
      CALL pio_netcdf_put_fvar (ng, model, AVG(ng)%name,                &
     &                          TRIM(Vname(1,idtime)), AVGtime(ng:),    &
     &                          (/AVG(ng)%Rindex/), (/1/),              &
     &                          pioFile = AVG(ng)%pioFile,              &
     &                          pioVar = AVG(ng)%pioVar(idtime)%vd)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out free-surface (m).
!
      IF (Aout(idFsur,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idFsur)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idFsur),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     AVERAGE(ng) % avgzeta)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idFsur)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

#  if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
!
!  Write out detided free-surface (m).
!
      IF (Aout(idFsuD,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idFsuD)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idFsuD),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     TIDES(ng) % zeta_detided)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idFsuD)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
!
!  Write out 2D momentum component (m/s) in the XI-direction.
!
      IF (Aout(idUbar,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idUbar)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_u2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_u2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idUbar),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % umask_avg,                        &
#  endif
     &                     AVERAGE(ng) % avgu2d)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbar)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

#  if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
!
!  Write out detided 2D momentum component (m/s) in the XI-direction.
!
      IF (Aout(idu2dD,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idu2dD)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_u2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_u2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idu2dD),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % umask_avg,                        &
#   endif
     &                     TIDES(ng) % ubar_detided)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idu2dD)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
!
!  Write out 2D momentum component (m/s) in the ETA-direction.
!
      IF (Aout(idVbar,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idVbar)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_v2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_v2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idVbar),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % vmask_avg,                        &
#  endif
     &                     AVERAGE(ng) % avgv2d)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbar)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

#  if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
!
!  Write out detided 2D momentum component (m/s) in the ETA-direction.
!
      IF (Aout(idv2dD,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idv2dD)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_v2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_v2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idv2dD),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % vmask_avg,                        &
#   endif
     &                     TIDES(ng) % vbar_detided)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idv2dD)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
!
!  Write out 2D Eastward momentum component (m/s) at RHO-points.
!
      IF (Aout(idu2dE,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idu2dE)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idu2dE),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_avg,                        &
#  endif
     &                     AVERAGE(ng) % avgu2dE)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idu2dE)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D Northward momentum component (m/s) at RHO-points.
!
      IF (Aout(idv2dN,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idv2dN)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idv2dN),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_avg,                        &
#  endif
     &                     AVERAGE(ng) % avgv2dN)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idv2dN)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

#  if defined FORWARD_WRITE && defined SOLVE3D
!
!  Write out time-averaged mass fluxes for 3D momentum coupling.
!
      IF (Aout(idUfx1,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idUfx1)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_u2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_u2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idUfx1),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % umask_avg,                        &
#   endif
     &                     AVERAGE(ng) % avgDU_avg1)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUfx1)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
      IF (Aout(idUfx2,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idUfx2)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_u2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_u2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idUfx2),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % umask_avg,                        &
#   endif
     &                     AVERAGE(ng) % avgDU_avg2)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUfx2)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
      IF (Aout(idVfx1,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idVfx1)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_v2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_v2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idVfx1),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % vmask_avg,                        &
#   endif
     &                     AVERAGE(ng) % avgDV_avg1)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVfx1)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
      IF (Aout(idVfx2,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idVfx2)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_v2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_v2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idVfx2),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % umask_avg,                        &
#   endif
     &                     AVERAGE(ng) % avgDV_avg2)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVfx2)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif

#  ifdef SOLVE3D
!
!  Write out 3D momentum component (m/s) in the XI-direction.
!
      IF (Aout(idUvel,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idUvel)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_u3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_u3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idUvel),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % umask_avg,                        &
#   endif
     &                     AVERAGE(ng) % avgu3d)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUvel)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

#   if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
!
!  Write out detided 3D momentum component (m/s) in the XI-direction.
!
      IF (Aout(idu3dD,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idu3dD)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_u3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_u3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idu3dD),                      &
     &                     ioDesc,                                      &
     &                     AVG(ng)%Rindex,                              &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#    ifdef MASKING
     &                     GRID(ng) % umask_avg,                        &
#    endif
     &                     TIDES(ng) % u_detided)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idu3dD)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#   endif
!
!  Write out 3D momentum component (m/s) in the ETA-direction.
!
      IF (Aout(idVvel,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idVvel)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_v3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_v3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idVvel),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % vmask_avg,                        &
#   endif
     &                     AVERAGE(ng) % avgv3d)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVvel)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

#   if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
!
!  Write out detided 3D momentum component (m/s) in the ETA-direction.
!
      IF (Aout(idv3dD,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idv3dD)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_v3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_v3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idv3dD),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#    ifdef MASKING
     &                     GRID(ng) % vmask_avg,                        &
#    endif
     &                     TIDES(ng) % v_detided)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idv3dD)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#   endif
!
!  Write out 3D Eastward momentum component (m/s) at RHO-points.
!
      IF (Aout(idu3dE,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idu3dE)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idu3dE),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % rmask_avg,                        &
#   endif
     &                     AVERAGE(ng) % avgu3dE)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idu3dE)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D Northward momentum component (m/s) at RHO-points.
!
      IF (Aout(idv3dN,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idv3dN)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idv3dN),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % rmask_avg,                        &
#   endif
     &                     AVERAGE(ng) % avgv3dN)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idv3dN)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out S-coordinate omega vertical velocity (m/s).
!
      IF (Aout(idOvel,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idOvel)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_w3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_w3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idOvel),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     AVERAGE(ng) % avgw3d)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idOvel)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out "true" vertical velocity (m/s).
!
      IF (Aout(idWvel,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idWvel)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_w3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_w3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idWvel),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     AVERAGE(ng) % avgwvel)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idOvel)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out tracer type variables.
!
      DO itrc=1,NT(ng)
        IF (Aout(idTvar(itrc),ng)) THEN
          scale=1.0_dp
          IF (AVG(ng)%pioTrc(itrc)%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_r3dvar(ng)
          ELSE
            ioDesc => ioDesc_sp_r3dvar(ng)
          END IF
          status=nf_fwrite3d(ng, model, AVG(ng)%pioFile,                &
     &                       AVG(ng)%pioTrc(itrc),                      &
     &                       AVG(ng)%Rindex,                            &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       AVERAGE(ng) % avgt(:,:,:,itrc))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idTvar(itrc))),            &
     &                          AVG(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO

#   if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
!
!  Write out detided active tracer variables.
!
      DO itrc=1,NAT
        IF (Aout(idTrcD(itrc),ng)) THEN
          scale=1.0_dp
          IF (AVG(ng)%pioVar(idTrcD(itrc))%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_r3dvar(ng)
          ELSE
            ioDesc => ioDesc_sp_r3dvar(ng)
          END IF
          status=nf_fwrite3d(ng, model, AVG(ng)%pioFile,                &
     &                       AVG(ng)%pioVar(idTrcD(itrc)),              &
     &                       AVG(ng)%Rindex,                            &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#    ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#    endif
     &                       TIDES(ng) % t_detided(:,:,:,itrc))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idTrcD(itrc))),            &
     &                          AVG(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
#   endif
!
!  Write out density anomaly.
!
      IF (Aout(idDano,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idDano)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idDano),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     AVERAGE(ng) % avgrho)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idDano)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

#   ifdef LMD_SKPP
!
!  Write out depth of surface boundary layer.
!
      IF (Aout(idHsbl,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idHsbl)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idHsbl),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#    endif
     &                     AVERAGE(ng) % avghsbl)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idHsbl)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#   endif
#   ifdef LMD_BKPP
!
!  Write out depth of bottom boundary layer.
!
      IF (Aout(idHbbl,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idHbbl)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idHbbl),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#    endif
     &                     AVERAGE(ng) % avghbbl)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idHbbl)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#   endif
#  endif
!
!  Write out 2D potential vorticity.
!
      IF (Aout(id2dPV,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(id2dPV)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_p2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_p2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(id2dPV),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % pmask,                            &
#  endif
     &                     AVERAGE(ng) % avgpvor2d)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,id2dPV)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D relative vorticity.
!
      IF (Aout(id2dRV,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(id2dRV)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_p2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_p2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(id2dRV),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % pmask,                            &
#  endif
     &                     AVERAGE(ng) % avgrvor2d)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,id2dRV)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

#  ifdef SOLVE3D
!
!  Write out 3D potential vorticity.
!
      IF (Aout(id3dPV,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(id3dPV)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_p3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_p3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(id3dPV),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % pmask,                            &
#   endif
     &                     AVERAGE(ng) % avgpvor3d)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,id3dPV)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D relative vorticity.
!
      IF (Aout(id3dRV,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(id3dRV)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_p3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_p3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(id3dRV),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % pmask,                            &
#   endif
     &                     AVERAGE(ng) % avgrvor3d)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,id3dRV)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
!
!  Write out quadratic <zeta*zeta> term.
!
      IF (Aout(idZZav,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idZZav)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idZZav),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     AVERAGE(ng) % avgZZ)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idZZav)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out quadratic <ubar*ubar> term.
!
      IF (Aout(idU2av,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idU2av)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_u2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_u2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idU2av),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % umask,                            &
#  endif
     &                     AVERAGE(ng) % avgU2)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idU2av)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out quadratic <vbar*vbar> term.
!
      IF (Aout(idV2av,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idV2av)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_v2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_v2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idV2av),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % vmask,                            &
#  endif
     &                     AVERAGE(ng) % avgV2)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idV2av)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

#  ifdef SOLVE3D
!
!  Write out u-volume flux.
!
      IF (Aout(idHUav,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idHUav)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_u3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_u3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idHUav),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % umask,                            &
#   endif
     &                     AVERAGE(ng) % avgHuon)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idHUav)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out v-volume flux.
!
      IF (Aout(idHVav,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idHVav)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_v3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_v3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idHVav),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % vmask,                            &
#   endif
     &                     AVERAGE(ng) % avgHvom)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idHVav)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out quadratic <u*u> term.
!
      IF (Aout(idUUav,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idUUav)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_u3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_u3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idUUav),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % umask,                            &
#   endif
     &                     AVERAGE(ng) % avgUU)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUUav)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out quadratic <u*v> term.
!
      IF (Aout(idUVav,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idUVav)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idUVav),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     AVERAGE(ng) % avgUV)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUVav)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out quadratic <v*v> term.
!
      IF (Aout(idVVav,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idVVav)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_v3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_v3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idVVav),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % vmask,                            &
#   endif
     &                     AVERAGE(ng) % avgVV)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVVav)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out quadratic <t*t> term.
!
      DO i=1,NT(ng)
        IF (Aout(idTTav(i),ng)) THEN
          scale=1.0_dp
          IF (AVG(ng)%pioVar(idTTav(i))%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_r3dvar(ng)
          ELSE
            ioDesc => ioDesc_sp_r3dvar(ng)
          END IF
          status=nf_fwrite3d(ng, model, AVG(ng)%pioFile,                &
     &                       AVG(ng)%pioVar(idTTav(i)),                 &
     &                       AVG(ng)%Rindex,                            &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       AVERAGE(ng) % avgTT(:,:,:,i))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idTTav(i))), AVG(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
!
!  Write out active tracer volume fluxes.
!
      DO i=1,NT(ng)
        IF (Aout(iHUTav(i),ng)) THEN
          scale=1.0_dp
          IF (AVG(ng)%pioVar(iHUTav(i))%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_u3dvar(ng)
          ELSE
            ioDesc => ioDesc_sp_u3dvar(ng)
          END IF
          status=nf_fwrite3d(ng, model, AVG(ng)%pioFile,                &
     &                       AVG(ng)%pioVar(iHUTav(i)),                 &
     &                       AVG(ng)%Rindex,                            &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#   ifdef MASKING
     &                       GRID(ng) % umask,                          &
#   endif
     &                       AVERAGE(ng) % avgHuonT(:,:,:,i))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,iHUTav(i))), AVG(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
!
        IF (Aout(iHVTav(i),ng)) THEN
          scale=1.0_dp
          IF (AVG(ng)%pioVar(iHVTav(i))%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_v3dvar(ng)
          ELSE
            ioDesc => ioDesc_sp_v3dvar(ng)
          END IF
          status=nf_fwrite3d(ng, model, AVG(ng)%pioFile,                &
     &                       AVG(ng)%pioVar(iHVTav(i)),                 &
     &                       AVG(ng)%Rindex,                            &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#   ifdef MASKING
     &                       GRID(ng) % vmask,                          &
#   endif
     &                       AVERAGE(ng) % avgHvomT(:,:,:,i))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,iHVTav(i))), AVG(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
!
!  Write out quadratic <u*t> and <v*t> terms.
!
      DO i=1,NT(ng)
        IF (Aout(idUTav(i),ng)) THEN
          scale=1.0_dp
          IF (AVG(ng)%pioVar(idUTav(i))%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_u3dvar(ng)
          ELSE
            ioDesc => ioDesc_sp_u3dvar(ng)
          END IF
          status=nf_fwrite3d(ng, model, AVG(ng)%pioFile,                &
     &                       AVG(ng)%pioVar(idUTav(i)),                 &
     &                       AVG(ng)%Rindex,                            &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#   ifdef MASKING
     &                       GRID(ng) % umask,                          &
#   endif
     &                       AVERAGE(ng) % avgUT(:,:,:,i))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idUTav(i))), AVG(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
!
        IF (Aout(idVTav(i),ng)) THEN
          scale=1.0_dp
          IF (AVG(ng)%pioVar(idVTav(i))%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_v3dvar(ng)
          ELSE
            ioDesc => ioDesc_sp_v3dvar(ng)
          END IF
          status=nf_fwrite3d(ng, model, AVG(ng)%pioFile,                &
     &                       AVG(ng)%pioVar(idVTav(i)),                 &
     &                       AVG(ng)%Rindex,                            &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#   ifdef MASKING
     &                       GRID(ng) % vmask,                          &
#   endif
     &                       AVERAGE(ng) % avgVT(:,:,:,i))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idVTav(i))), AVG(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO

#   if defined LMD_MIXING || defined MY25_MIXING || defined GLS_MIXING
!
!  Write out vertical viscosity coefficient.
!
      IF (Aout(idVvis,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idVvis)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_w3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_w3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idVvis),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#    ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#    endif
     &                     AVERAGE(ng) % avgAKv,                        &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVvis)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out vertical diffusion coefficient for potential temperature.
!
      IF (Aout(idTdif,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idTdif)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_w3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_w3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idTdif),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#    ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#    endif
     &                     AVERAGE(ng) % avgAKt,                        &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTdif)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

#    ifdef SALINITY
!
!  Write out vertical diffusion coefficient for salinity.
!
      IF (Aout(idSdif,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idSdif)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_w3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_w3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idSdif),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#     ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#     endif
     &                     AVERAGE(ng) % avgAKs,                        &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSdif)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#    endif
#   endif
#   if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
!
!  Write out surface air pressure.
!
      IF (Aout(idPair,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idPair)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idPair),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#    endif
     &                     AVERAGE(ng) % avgPair)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idPair)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#   endif
#   if defined BULK_FLUXES
!
!  Write out surface air temperature.
!
      IF (Aout(idTair,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idTair)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idTair),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#    endif
     &                     AVERAGE(ng) % avgTair)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTair)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#   endif
#   if defined BULK_FLUXES || defined ECOSIM
!
!  Write out surface winds.
!
      IF (Aout(idUair,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idUair)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idUair),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#    endif
     &                     AVERAGE(ng) % avgUwind)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUair)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
      IF (Aout(idVair,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idVair)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idVair),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#    endif
     &                     AVERAGE(ng) % avgVwind)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVair)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out Eastward/Northward surface wind (m/s) at RHO-points.
!
      IF (Aout(idUaiE,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idUaiE)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idUaiE),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % rmask_avg,                        &
#    endif
     &                     AVERAGE(ng) % avguwindE)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUaiE)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
      IF (Aout(idVaiN,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idVaiN)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idVaiN),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % rmask_avg,                        &
#    endif
     &                     AVERAGE(ng) % avgvwindN)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVaiN)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#   endif
!
!  Write out surface net heat flux.
!
      IF (Aout(idTsur(itemp),ng)) THEN
#   ifdef ADJOINT
#    if defined AD_SENSITIVITY   || defined I4DVAR_ANA_SENSITIVITY || \
       defined OPT_OBSERVATIONS
        scale=1.0_dp/(rho0*Cp)
#    else
        scale=1.0_dp
#    endif
#   else
        scale=rho0*Cp
#   endif
        IF (AVG(ng)%pioVar(idTsur(itemp))%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idTsur(itemp)),               &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     AVERAGE(ng) % avgstf)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTsur(itemp))),             &
     &                        AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

#   ifdef SALINITY
!
!  Write out surface salt flux  (PSU m/s = kg salt/m2/s).
!
      IF (Aout(idTsur(isalt),ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idTsur(isalt))%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idTsur(isalt)),               &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#    endif
     &                     AVERAGE(ng) % avgswf)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTsur(isalt))),             &
     &                        AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#   endif

#   if defined BULK_FLUXES || defined FRC_COUPLING
!
!  Write out latent heat flux.
!
      IF (Aout(idLhea,ng)) THEN
#    ifdef ADJOINT
        scale=1.0_dp
#    else
        scale=rho0*Cp
#    endif
        IF (AVG(ng)%pioVar(idLhea)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idLhea),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#    endif
     &                     AVERAGE(ng) % avglhf)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idLhea)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out sensible heat flux.
!
      IF (Aout(idShea,ng)) THEN
#    ifdef ADJOINT
        scale=1.0_dp
#    else
        scale=rho0*Cp
#    endif
        IF (AVG(ng)%pioVar(idShea)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idShea),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#    endif
     &                     AVERAGE(ng) % avgshf)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idShea)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out net longwave radiation flux.
!
      IF (Aout(idLrad,ng)) THEN
#    ifdef ADJOINT
        scale=1.0_dp
#    else
        scale=rho0*Cp
#    endif
        IF (AVG(ng)%pioVar(idLrad)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idLrad),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#    endif
     &                     AVERAGE(ng) % avglrf)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idLrad)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#   endif

#   ifdef BULK_FLUXES
#    ifdef EMINUSP
!
!  Write out evaportaion rate (kg/m2/s).
!
      IF (Aout(idevap,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idevap)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idevap),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#     ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#     endif
     &                     AVERAGE(ng) % avgevap)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idevap)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out precipitation rate (kg/m2/s).
!
      IF (Aout(idrain,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idrain)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idrain),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#     ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#     endif
     &                     AVERAGE(ng) % avgrain)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idrain)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#    endif
#   endif
#   ifdef SHORTWAVE
!
!  Write out net shortwave radiation flux.
!
      IF (Aout(idSrad,ng)) THEN
#     ifdef ADJOINT
        scale=1.0_dp
#     else
        scale=rho0*Cp
#     endif
        IF (AVG(ng)%pioVar(idSrad)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idSrad),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#     ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#     endif
     &                     AVERAGE(ng) % avgsrf)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSrad)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#   endif
#  endif
!
!  Write out surface u-momentum stress.
!
      IF (Aout(idUsms,ng)) THEN
#  ifdef ADJOINT
#   if defined AD_SENSITIVITY   || defined I4DVAR_ANA_SENSITIVITY || \
      defined OPT_OBSERVATIONS
        scale=1.0_dp/rho0
#   else
        scale=1.0_dp
#   endif
#  else
        scale=rho0
#  endif
        IF (AVG(ng)%pioVar(idUsms)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_u2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_u2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idUsms),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % umask,                            &
#  endif
     &                     AVERAGE(ng) % avgsus)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUsms)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface v-momentum stress.
!
      IF (Aout(idVsms,ng)) THEN
#  ifdef ADJOINT
#   if defined AD_SENSITIVITY   || defined I4DVAR_ANA_SENSITIVITY || \
      defined OPT_OBSERVATIONS
        scale=1.0_dp/rho0
#   else
        scale=1.0_dp
#   endif
#  else
        scale=rho0
#  endif
        IF (AVG(ng)%pioVar(idVsms)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_v2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_v2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idVsms),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % vmask,                            &
#  endif
     &                     AVERAGE(ng) % avgsvs)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVsms)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom u-momentum stress.
!
      IF (Aout(idUbms,ng)) THEN
#  ifdef ADJOINT
#   if defined AD_SENSITIVITY   || defined I4DVAR_ANA_SENSITIVITY || \
      defined OPT_OBSERVATIONS
        scale=1.0_dp/rho0
#   else
        scale=1.0_dp
#   endif
#  else
        scale=rho0
#  endif
        IF (AVG(ng)%pioVar(idUbms)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_u2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_u2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idUbms),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % umask,                            &
#  endif
     &                     AVERAGE(ng) % avgbus)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbms)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom v-momentum stress.
!
      IF (Aout(idVbms,ng)) THEN
#  ifdef ADJOINT
#   if defined AD_SENSITIVITY   || defined I4DVAR_ANA_SENSITIVITY || \
      defined OPT_OBSERVATIONS
        scale=1.0_dp/rho0
#   else
        scale=1.0_dp
#   endif
#  else
        scale=rho0
#  endif
        IF (AVG(ng)%pioVar(idVbms)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_v2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_v2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idVbms),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % vmask,                            &
#  endif
     &                     AVERAGE(ng) % avgbvs)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbms)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

#  ifdef BBL_MODEL
!
!  Write out bottom momentum U-stress at RHO-!points.
!
      IF (Aout(idUbrs,ng)) THEN
#  ifdef ADJOINT
#   if defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
      defined OPT_OBSERVATIONS
        scale=1.0_dp/rho0
#   else
        scale=1.0_dp
#   endif
#  else
        scale=rho0
#  endif
        IF (AVG(ng)%pioVar(idUbrs)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idUbrs),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     AVERAGE(ng) % avgUbrs)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbrs)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom momentum V-stress at RHO-points.
!
      IF (Aout(idVbrs,ng)) THEN
#  ifdef ADJOINT
#   if defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
      defined OPT_OBSERVATIONS
        scale=1.0_dp/rho0
#   else
        scale=1.0_dp
#   endif
#  else
        scale=rho0
#  endif
        IF (AVG(ng)%pioVar(idVbrs)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idVbrs),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     AVERAGE(ng) % avgVbrs)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbrs)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom wave U-stress at RHO-points.
!
      IF (Aout(idUbws,ng)) THEN
#  ifdef ADJOINT
#   if defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
      defined OPT_OBSERVATIONS
        scale=1.0_dp/rho0
#   else
        scale=1.0_dp
#   endif
#  else
        scale=rho0
#  endif
        IF (AVG(ng)%pioVar(idUbws)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idUbws),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     AVERAGE(ng) % avgUbws)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbws)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom wave V-stress at RHO-points.
!
      IF (Aout(idVbws,ng)) THEN
#  ifdef ADJOINT
#   if defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
      defined OPT_OBSERVATIONS
        scale=1.0_dp/rho0
#   else
        scale=1.0_dp
#   endif
#  else
        scale=rho0
#  endif
        IF (AVG(ng)%pioVar(idVbws)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idVbws),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     AVERAGE(ng) % avgVbws)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbws)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom max wave current U-stress at RHO-points.
!
      IF (Aout(idUbcs,ng)) THEN
#  ifdef ADJOINT
#   if defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
      defined OPT_OBSERVATIONS
        scale=1.0_dp/rho0
#   else
        scale=1.0_dp
#   endif
#  else
        scale=rho0
#  endif
        IF (AVG(ng)%pioVar(idUbcs)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idUbcs),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     AVERAGE(ng) % avgUbcs)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbcs)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out max bottom wave current V-stress at RHO-points.
!
      IF (Aout(idVbcs,ng)) THEN
#  ifdef ADJOINT
#   if defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
      defined OPT_OBSERVATIONS
        scale=1.0_dp/rho0
#   else
        scale=1.0_dp
#   endif
#  else
        scale=rho0
#  endif
        IF (AVG(ng)%pioVar(idVbcs)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idVbcs),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     AVERAGE(ng) % avgVbcs)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbcs)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out max bottom wave current stress magnitude at RHO-points.
!
      IF (Aout(idUVwc,ng)) THEN
#  ifdef ADJOINT
#   if defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
      defined OPT_OBSERVATIONS
        scale=1.0_dp/rho0
#   else
        scale=1.0_dp
#   endif
#  else
        scale=rho0
#  endif
        IF (AVG(ng)%pioVar(idUVwc)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idUVwc),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     AVERAGE(ng) % avgUVwc)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUVwc)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom wave orbital U-velocity at RHO-points.
!
      IF (Aout(idUbot,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idUbot)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idUbot),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     AVERAGE(ng) % avgUbot)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbot)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom wave orbital V-velocity at RHO-points.
!
      IF (Aout(idVbot,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idVbot)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idVbot),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     AVERAGE(ng) % avgVbot)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbot)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom U-current above bed at RHO-points.
!
      IF (Aout(idUbur,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idUbur)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idUbur),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     AVERAGE(ng) % avgUbur)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbur)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom V-current above the bed at RHO-points.
!
      IF (Aout(idVbvr,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idVbvr)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idVbvr),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     AVERAGE(ng) % avgVbvr)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbvr)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif

#  if defined SEDIMENT && defined BEDLOAD
!
!  Write out bed load transport in U-direction.
!
      DO itrc=1,NST
        IF (Aout(idUbld(itrc),ng)) THEN
          scale=1.0_dp
          IF (AVG(ng)%pioVar(idUbld(itrc))%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_u2dvar(ng)
          ELSE
            ioDesc => ioDesc_sp_u2dvar(ng)
          END IF
          status=nf_fwrite2d(ng, model, AVG(ng)%pioFile,                &
     &                       AVG(ng)%pioVar(idUbld(itrc)),              &
     &                       AVG(ng)%Rindex,                            &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % umask,                          &
#   endif
     &                       SEDBED(ng) % avgbedldu(:,:,itrc))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idUbld(itrc))),            &
     &                          AVG(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
!
!  Write out bed load transport in V-direction.
!
        IF (Aout(idVbld(itrc),ng)) THEN
          scale=1.0_dp
          IF (AVG(ng)%pioVar(idVbld(itrc))%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_v2dvar(ng)
          ELSE
            ioDesc => ioDesc_sp_v2dvar(ng)
          END IF
          status=nf_fwrite2d(ng, model, AVG(ng)%pioFile,                &
     &                       AVG(ng)%pioVar(idVbld(itrc)),              &
     &                       AVG(ng)%Rindex,                            &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % vmask,                          &
#   endif
     &                       SEDBED(ng) % avgbedldv(:,:,itrc))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idVbld(itrc))),            &
     &                          AVG(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
#  endif
#  ifdef WAVES_HEIGHT
!
!  Write out wind-induced significant wave height.
!
      IF (Aout(idWamp,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idWamp)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idWamp),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     AVERAGE(ng) % avgWamp)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWamp)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out wind-induced significant wave height squared.
!
      IF (Aout(idWam2,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idWam2)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idWam2),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     AVERAGE(ng) % avgWam2)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWam2)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
#  ifdef WAVES_LENGTH
!
!  Write out wind-induced mean wavelength.
!
      IF (Aout(idWlen,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idWlen)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idWlen),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     AVERAGE(ng) % avgWlen)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWlen)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
#  ifdef WAVES_LENGTHP
!
!  Write out wind-induced peak wave wavelength.
!
      IF (Aout(idWlep,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idWlep),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     AVERAGE(ng) % avgWlep)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWlep)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
#  ifdef WAVES_DIR
!
!  Write out wind-induced mean wave direction.
!
      IF (Aout(idWdir,ng)) THEN
        scale=rad2deg
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idWdir),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     AVERAGE(ng) % avgWdir)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWdir)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
#  ifdef WAVES_DIRP
!
!  Write out wind-induced peak wave direction.
!
      IF (Aout(idWdip,ng)) THEN
        scale=rad2deg
        IF (AVG(ng)%pioVar(idWdip)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idWdip),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     AVERAGE(ng) % avgWdip)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWdip)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
#  ifdef WAVES_TOP_PERIOD
!
!  Write out wind-induced surface wave period.
!
      IF (Aout(idWptp,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idWptp),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     AVERAGE(ng) % avgWptp)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWptp)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
#  ifdef WAVES_BOT_PERIOD
!
!  Write out wind-induced bottom wave period.

      IF (Aout(idWpbt,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idWpbt)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idWpbt),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     AVERAGE(ng) % avgWpbt)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWpbt)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
#  if defined BBL_MODEL    || defined BEDLOAD_VANDERA || \
     defined WAV_COUPLING
!
!  Write out wind-induced bottom orbital velocity.
!
      IF (Aout(idWorb,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idWorb)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idWorb),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     AVERAGE(ng) % avgWorb)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWorb)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
#  ifdef UV_KIRBY
!
!  Write out U-velocity from Kirby and Chen.
!
      IF (Aout(idUwav,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idUwav)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idUwav),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     AVERAGE(ng) % avgUwav)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUwav)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out V-velocity from Kirby and Chen.
!
      IF (Aout(idVwav,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idVwav)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idVwav),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     AVERAGE(ng) % avgVwav)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVwav)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
#  ifdef ICE_MODEL
!
!-----------------------------------------------------------------------
!  Write out sea-ice model variables.
!-----------------------------------------------------------------------
!
        CALL ice_wrt_pio (ng, model, tile,                              &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    Aout, AVG)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
#  if defined WEC || defined WEC_VF
!
!-----------------------------------------------------------------------
!  Write out Waves Effect on Current variables.
!-----------------------------------------------------------------------
!
      CALL wec_wrt_pio (ng, model, tile,                                &
    &                   LBi, UBi, LBj, UBj,                             &
    &                   Aout, AVG)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
!
!-----------------------------------------------------------------------
!  Synchronize time-average NetCDF file to disk to allow other processes
!  to access data immediately after it is written.
!-----------------------------------------------------------------------
!
      CALL pio_netcdf_sync (ng, model, AVG(ng)%name, AVG(ng)%pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  ifdef NESTING
      IF (Master) WRITE (stdout,20) AVG(ng)%Rindex, ng
#  else
      IF (Master) WRITE (stdout,20) AVG(ng)%Rindex
#  endif
!
  10  FORMAT (/,' WRT_AVG_PIO - error while writing variable: ',a,      &
     &        /,16x,'into averages NetCDF file for time record: ',i0)
  20  FORMAT (2x,'WRT_AVG_PIO      - wrote averaged',t40,'fields',t59,  &
#  ifdef NESTING
     &        'in record = ',i0,t92,i2.2)
#  else
     &        'in record = ',i0)
#  endif
!
      RETURN
      END SUBROUTINE wrt_avg_pio
# endif


#endif
      END MODULE wrt_avg_mod
