#include "cppdefs.h"
#if defined IS4DVAR && defined LANCZOS
      SUBROUTINE wrt_hessian (ng, kout, nout)
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine writes Hessian eigenvectors into output NetCDF file.   !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_forces
      USE mod_grid
      USE mod_iounits
      USE mod_mixing
      USE mod_ncparam
      USE mod_netcdf
      USE mod_ocean
      USE mod_scalars
      USE mod_stepping
!
# ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti
# endif
      USE nf_fwrite2d_mod, ONLY : nf_fwrite2d
# ifdef SOLVE3D
      USE nf_fwrite3d_mod, ONLY : nf_fwrite3d
# endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, kout, nout
!
!  Local variable declarations.
!
      integer :: LBi, UBi, LBj, UBj
      integer :: i, j, gfactor, gtype, status
# ifdef SOLVE3D
      integer :: itrc, k
# endif
      real(r8) :: scale
!
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)
!
!-----------------------------------------------------------------------
!  Write out Hessian eigenvectors.
!-----------------------------------------------------------------------
!
      IF (exit_flag.ne.NoError) RETURN
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
      tHSSindx(ng)=tHSSindx(ng)+1
      NrecHSS(ng)=NrecHSS(ng)+1
!
!  Write out model time (s).
!
      CALL netcdf_put_fvar (ng, iADM, HSSname(ng),                      &
     &                      TRIM(Vname(1,idtime)), time(ng:),           &
     &                      (/tHSSindx(ng)/), (/1/),                    &
     &                      ncid = ncHSSid(ng),                         &
     &                      varid = hssVid(idtime,ng))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out free-surface.
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf_fwrite2d(ng, iADM, ncHSSid(ng), hssVid(idFsur,ng),      &
     &                   tHSSindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
# ifdef MASKING
     &                   GRID(ng) % rmask,                              &
# endif
     &                   OCEAN(ng)% ad_zeta(:,:,kout))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idFsur)), tHSSindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out 2D momentum component in the XI-direction.
!
      scale=1.0_r8
      gtype=gfactor*u2dvar
      status=nf_fwrite2d(ng, iADM, ncHSSid(ng), hssVid(idUbar,ng),      &
     &                   tHSSindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
# ifdef MASKING
     &                   GRID(ng) % umask,                              &
# endif
     &                   OCEAN(ng) % ad_ubar(:,:,kout))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUbar)), tHSSindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out 2D momentum component in the ETA-direction.
!
      scale=1.0_r8
      gtype=gfactor*v2dvar
      status=nf_fwrite2d(ng, iADM, ncHSSid(ng), hssVid(idVbar,ng),      &
     &                   tHSSindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
# ifdef MASKING
      &                  GRID(ng) % vmask,                              &
# endif
     &                   OCEAN(ng) % ad_vbar(:,:,kout))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVbar)), tHSSindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
# ifdef SOLVE3D
!
!  Write out 3D momentum component in the XI-direction.
!
      scale=1.0_r8
      gtype=gfactor*u3dvar
      status=nf_fwrite3d(ng, iADM, ncHSSid(ng), hssVid(idUvel,ng),      &
     &                   tHSSindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#  ifdef MASKING
     &                   GRID(ng) % umask,                              &
#  endif
     &                   OCEAN(ng) % ad_u(:,:,:,nout))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUvel)), tHSSindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out 3D momentum component in the ETA-direction.
!
      scale=1.0_r8
      gtype=gfactor*v3dvar
      status=nf_fwrite3d(ng, iADM, ncHSSid(ng), hssVid(idVvel,ng),      &
     &                   tHSSindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#  ifdef MASKING
     &                   GRID(ng) % vmask,                              &
#  endif
     &                   OCEAN(ng) % ad_v(:,:,:,nout))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVvel)), tHSSindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out tracer type variables.
!
      DO itrc=1,NT(ng)
        scale=1.0_r8
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, iADM, ncHSSid(ng), hssTid(itrc,ng),      &
     &                     tHSSindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     OCEAN(ng) % ad_t(:,:,:,nout,itrc))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTvar(itrc))),              &
     &                        tHSSindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END DO
#  ifdef ADJUST_STFLUX
!
!  Write out surface net tracers fluxes. Notice that fluxes have their
!  own fixed time-dimension (of size Nfrec) to allow 4DVAR adjustments
!  at other times in addition to initialization time.
!
      DO itrc=1,NT(ng)
        scale=1.0_r8
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, iADM, ncHSSid(ng),                       &
     &                     hssVid(idTsur(itrc),ng),                     &
     &                     tHSSindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,     &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     FORCES(ng) % ad_tflux(:,:,:,nout,itrc))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTsur(itrc))),              &
     &                        tHSSindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END DO
#  endif
# endif
# ifdef ADJUST_WSTRESS
!
!  Write out surface U-momentum stress.  Notice that the stress has its
!  own fixed time-dimension (of size Nfrec) to allow 4DVAR adjustments
!  at other times in addition to initialization time.
!
      scale=1.0_r8
      gtype=gfactor*u3dvar
      status=nf_fwrite3d(ng, iADM, ncHSSid(ng), hssVid(idUsms,ng),      &
     &                   tHSSindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,       &
#  ifdef MASKING
     &                   GRID(ng) % umask,                              &
#  endif
     &                   FORCES(ng) % ad_ustr(:,:,:,kout))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUsms)), tHSSindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out surface V-momentum stress.
!
      scale=1.0_r8
      gtype=gfactor*v3dvar
      status=nf_fwrite3d(ng, iADM, ncHSSid(ng), hssVid(idVsms,ng),      &
     &                   tHSSindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#  ifdef MASKING
     &                   GRID(ng) % vmask,                              &
#  endif
     &                   FORCES(ng) % ad_vstr(:,:,:,kout))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVsms)), tHSSindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
# endif
!
!-----------------------------------------------------------------------
!  Synchronize Hessian eigenvectors NetCDF file to disk to allow other
!  processes to access data immediately after it is written.
!-----------------------------------------------------------------------
!
      IF (OutThread) THEN
        status=nf90_sync(ncHSSid(ng))
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,20) TRIM(HSSname(ng))
          exit_flag=3
          ioerror=status
        END IF
      END IF
# ifdef DISTRIBUTE
      CALL mp_bcasti (ng, iADM, exit_flag, 1)
# endif
# ifdef SOLVE3D
      IF (Master) WRITE (stdout,30) kout, nout, tHSSindx(ng)
# else
      IF (Master) WRITE (stdout,30) kout, tHSSindx(ng)
# endif
!
  10  FORMAT (/,' WRT_HESSIAN - error while writing variable: ',a,/,    &
     &        15x,'into Hessian NetCDF file for time record: ',i4)
  20  FORMAT (/,' WRT_HESSIAN - unable to synchronize to disk Hessian', &
     &        'file:',/,15x,a)
# ifdef SOLVE3D
  30  FORMAT (3x,'WRT_HESSIAN  - wrote Hessian  fields (Index=', i1,    &
     &        ',',i1,') into time record = ',i7.7)
# else
  30  FORMAT (3x,'WRT_HESSIAN  - wrote Hessian  fields (Index=', i1,    &
     &        ') into time record = ',i7.7)
# endif
      RETURN
      END SUBROUTINE wrt_hessian
#else
      SUBROUTINE wrt_hessian
      RETURN
      END SUBROUTINE wrt_hessian
#endif
