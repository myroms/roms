#include "cppdefs.h"
      SUBROUTINE wrt_rst (ng)
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine writes fields into restart NetCDF file.                !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
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
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      integer :: LBi, UBi, LBj, UBj
      integer :: gfactor, gtype, i, itrc, status

      integer :: nf_fwrite

      real(r8) :: scale
!
      LBi=lbound(GRID(ng)%h,DIM=1)
      UBi=ubound(GRID(ng)%h,DIM=1)
      LBj=lbound(GRID(ng)%h,DIM=2)
      UBj=ubound(GRID(ng)%h,DIM=2)
!
!-----------------------------------------------------------------------
!  Write out restart fields.
!-----------------------------------------------------------------------
!
      IF (exit_flag.ne.0) RETURN
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
      trstindx(ng)=trstindx(ng)+1
      nrecrst(ng)=nrecrst(ng)+1
!
!  If requested, set time index to recycle time records in restart
!  file.
!
      IF (lcycle(ng)) THEN
        trstindx(ng)=MOD(trstindx(ng)-1,2)+1
      END IF
!
!  Write out model time (s).
!
      IF (OutThread) THEN
        status=nf_put_var1_TYPE(ncrstid(ng), rstVid(idtime,ng),         &
     &                          trstindx(ng), time(ng))
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,10) TRIM(Vname(1,idtime)), trstindx(ng)
          exit_flag=3
          RETURN
        END IF
      END IF
!
!  Write out free-surface (m)
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf_fwrite(ng, ncrstid(ng), rstVid(idFsur,ng),              &
     &                 trstindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, 1, 1, scale,                 &
#ifdef MASKING
     &                 GRID(ng) % rmask(LBi,LBj),                       &
#endif
     &                 OCEAN(ng) % zeta(LBi,LBj,KOUT))
      IF (OutThread.and.(status.ne.nf_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idFsur)), trstindx(ng)
        exit_flag=3
        RETURN
      END IF
!
!  Write out 2D momentum component (m/s) in the XI-direction.
!
      scale=1.0_r8
      gtype=gfactor*u2dvar
      status=nf_fwrite(ng, ncrstid(ng), rstVid(idUbar,ng),              &
     &                 trstindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, 1, 1, scale,                 &
#ifdef MASKING
     &                 GRID(ng) % umask(LBi,LBj),                       &
#endif
     &                 OCEAN(ng) % ubar(LBi,LBj,KOUT))
      IF (OutThread.and.(status.ne.nf_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idUbar)), trstindx(ng)
        exit_flag=3
        RETURN
      END IF
!
!  Write out 2D momentum component (m/s) in the ETA-direction.
!
      scale=1.0_r8
      gtype=gfactor*v2dvar
      status=nf_fwrite(ng, ncrstid(ng), rstVid(idVbar,ng),              &
     &                 trstindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, 1, 1, scale,                 &
#ifdef MASKING
     &                 GRID(ng) % vmask(LBi,LBj),                       &
#endif
     &                 OCEAN(ng) % vbar(LBi,LBj,KOUT))
      IF (OutThread.and.(status.ne.nf_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idVbar)), trstindx(ng)
        exit_flag=3
        RETURN
      END IF
#ifdef SOLVE3D
!
!  Write out 3D momentum component (m/s) in the XI-direction.
!
      scale=1.0_r8
      gtype=gfactor*u3dvar
      status=nf_fwrite(ng, ncrstid(ng), rstVid(idUvel,ng),              &
     &                 trstindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, 1, N(ng), scale,             &
# ifdef MASKING
     &                 GRID(ng) % umask(LBi,LBj),                       &
# endif
     &                 OCEAN(ng) % u(LBi,LBj,1,NOUT))
      IF (OutThread.and.(status.ne.nf_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idUvel)), trstindx(ng)
        exit_flag=3
        RETURN
      END IF
!
!  Write out momentum component (m/s) in the ETA-direction.
!
      scale=1.0_r8
      gtype=gfactor*v3dvar
      status=nf_fwrite(ng, ncrstid(ng), rstVid(idVvel,ng),              &
     &                 trstindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, 1, N(ng), scale,             &
# ifdef MASKING
     &                 GRID(ng) % vmask(LBi,LBj),                       &
# endif
     &                 OCEAN(ng) % v(LBi,LBj,1,NOUT))
      IF (OutThread.and.(status.ne.nf_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idVvel)), trstindx(ng)
        exit_flag=3
        RETURN
      END IF
!
!  Write out tracer type variables.
!
      DO itrc=1,NT(ng)
        scale=1.0_r8
        gtype=gfactor*r3dvar
        status=nf_fwrite(ng, ncrstid(ng), rstTid(itrc,ng),              &
     &                   trstindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
# ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
# endif
     &                   OCEAN(ng) % t(LBi,LBj,1,NOUT,itrc))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idTvar(itrc))), trstindx(ng)
          exit_flag=3
          RETURN
        END IF
      END DO
!
!  Write out density anomaly.
!
      scale=1.0_r8
      gtype=gfactor*r3dvar
      status=nf_fwrite(ng, ncrstid(ng), rstVid(idDano,ng),              &
     &                 trstindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, 1, N(ng), scale,             &
# ifdef MASKING
     &                 GRID(ng) % rmask(LBi,LBj),                       &
# endif
     &                 OCEAN(ng) % rho(LBi,LBj,1))
      IF (OutThread.and.(status.ne.nf_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idDano)), trstindx(ng)
        exit_flag=3
        RETURN
      END IF
# ifdef LMD_SKPP
!
!  Write out depth of surface boundary layer.
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf_fwrite(ng, ncrstid(ng), rstVid(idHsbl,ng),              &
     &                 trstindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, 1, 1, scale,                 &
#  ifdef MASKING
     &                 GRID(ng) % rmask(LBi,LBj),                       &
#  endif
     &                 MIXING(ng) % hsbl(LBi,LBj))
      IF (OutThread.and.(status.ne.nf_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idHsbl)), trstindx(ng)
        exit_flag=3
        RETURN
      END IF
# endif
# ifdef LMD_BKPP
!
!  Write out depth of bottom boundary layer.
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf_fwrite(ng, ncrstid(ng), rstVid(idHbbl,ng),              &
     &                 trstindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, 1, 1, scale,                 &
#  ifdef MASKING
     &                 GRID(ng) % rmask(LBi,LBj),                       &
#  endif
     &                 MIXING(ng) % hbbl(LBi,LBj))
      IF (OutThread.and.(status.ne.nf_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idHbbl)), trstindx(ng)
        exit_flag=3
        RETURN
      END IF
# endif
# if defined BVF_MIXING  || defined GLS_MIXING || \
     defined MY25_MIXING || defined LMD_MIXING
!
!  Write out vertical viscosity coefficient.
!
      scale=1.0_r8
      gtype=gfactor*r3dvar
      status=nf_fwrite(ng, ncrstid(ng), rstVid(idVvis,ng),              &
     &                 trstindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, 0, N(ng), scale,             &
#  ifdef MASKING
     &                 GRID(ng) % rmask(LBi,LBj),                       &
#  endif
     &                 MIXING(ng) % Akv(LBi,LBj,0))
      IF (OutThread.and.(status.ne.nf_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idVvis)), trstindx(ng)
        exit_flag=3
        RETURN
      END IF
!
!  Write out vertical diffusion coefficient for potential temperature.
!
      scale=1.0_r8
      gtype=gfactor*r3dvar
      status=nf_fwrite(ng, ncrstid(ng), rstVid(idTdif,ng),              &
     &                 trstindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, 0, N(ng), scale,             &
#  ifdef MASKING
     &                 GRID(ng) % rmask(LBi,LBj),                       &
#  endif
     &                 MIXING(ng) % Akt(LBi,LBj,0,itemp))
      IF (OutThread.and.(status.ne.nf_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idTdif)), trstindx(ng)
        exit_flag=3
        RETURN
      END IF
#  ifdef SALINITY
!
!  Write out vertical diffusion coefficient for salinity.
!
      scale=1.0_r8
      gtype=gfactor*r3dvar
      status=nf_fwrite(ng, ncrstid(ng), rstVid(idSdif,ng),              &
     &                 trstindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, 0, N(ng), scale,             &
#   ifdef MASKING
     &                 GRID(ng) % rmask(LBi,LBj),                       &
#   endif
     &                 MIXING(ng) % Akt(LBi,LBj,0,isalt))
      IF (OutThread.and.(status.ne.nf_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idSdif)), trstindx(ng)
        exit_flag=3
        RETURN
      END IF
#  endif
# endif
# ifdef SEDIMENT
!
!  Write out sediment fraction of each size class in each bed layer.
!
      DO i=1,NST
        scale=1.0_r8
        gtype=gfactor*b3dvar
        status=nf_fwrite(ng, ncrstid(ng), rstVid(idfrac(i),ng),         &
     &                   trstindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, Nbed, scale,            &
#  ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#  endif
     &                   OCEAN(ng) % bed_frac(LBi,LBj,1,i))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idfrac(i))), trstindx(ng)
          exit_flag=3
          RETURN
        END IF
      END DO
!
!  Write out sediment mass of each size class in each bed layer.
!
      DO i=1,NST
        scale=1.0_r8
        gtype=gfactor*b3dvar
        status=nf_fwrite(ng, ncrstid(ng), rstVid(idBmas(i),ng),         &
     &                   trstindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, Nbed, scale,            &
#  ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#  endif
     &                   OCEAN(ng) % bed_mass(LBi,LBj,1,i))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idBmas(i))), trstindx(ng)
          exit_flag=3
          RETURN
        END IF
      END DO
!
!  Write out sediment properties in each bed layer.
!
      DO i=1,MBEDP
        scale=1.0_r8
        gtype=gfactor*b3dvar
        status=nf_fwrite(ng, ncrstid(ng), rstVid(idSbed(i),ng),         &
     &                   trstindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, Nbed, scale,            &
#  ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#  endif
     &                   OCEAN(ng) % bed(LBi,LBj,1,i))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idSbed(i))), trstindx(ng)
          exit_flag=3
          RETURN
        END IF
      END DO
# endif
# if defined SEDIMENT || defined BBL_MODEL
!
!  Write out exposed sediment layer properties. Notice that only the
!  first four properties (mean grain diameter, mean grain density,
!  mean settling velocity, and mean critical erosion stress) are
!  written.
!
      DO i=1,4
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite(ng, ncrstid(ng), rstVid(idBott(i),ng),         &
     &                   trstindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, Nbed, scale,            &
#  ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#  endif
     &                   OCEAN(ng) % bottom(LBi,LBj,i))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idBott(i))), trstindx(ng)
          exit_flag=3
          RETURN
        END IF
      END DO
# endif
#endif
!
!  Synchronize history NetCDF file to disk.
!
      IF (OutThread) THEN
        status=nf_sync(ncrstid(ng))
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,20)
          exit_flag=3
          RETURN
        END IF
        WRITE (stdout,30) trstindx(ng)
      END IF
!
  10  FORMAT (/,' WRT_RST - error while writing variable: ',a,/,11x,    &
     &        'into restart NetCDF file for time record: ',i4)
  20  FORMAT (/,' WRT_RST - unable to synchronize restart NetCDF to ',  &
     &        'disk.')
  30  FORMAT (6x,'WRT_RST   - wrote re-start fields into time ',        &
     &        'record:',t63,i10.10)
      RETURN
      END SUBROUTINE wrt_rst
