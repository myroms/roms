#include "cppdefs.h"
      SUBROUTINE wrt_rst (ng)
!
!=======================================================================
!  Copyright (c) 2005 ROMS/TOMS Group                                  !
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
#if defined SEDIMENT || defined BBL_MODEL
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
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)
!
!-----------------------------------------------------------------------
!  Write out restart fields.
!-----------------------------------------------------------------------
!
      IF (exit_flag.ne.NoError) RETURN
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
      tRSTindx(ng)=tRSTindx(ng)+1
      NrecRST(ng)=NrecRST(ng)+1
!
!  If requested, set time index to recycle time records in restart
!  file.
!
      IF (LcycleRST(ng)) THEN
        tRSTindx(ng)=MOD(tRSTindx(ng)-1,2)+1
      END IF
!
!  Write out model time (s).
!
      IF (OutThread) THEN
        status=nf_put_var1_TYPE(ncRSTid(ng), rstVid(idtime,ng),         &
     &                          tRSTindx(ng), time(ng))
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,10) TRIM(Vname(1,idtime)), tRSTindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out free-surface (m)
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf_fwrite(ng, iNLM, ncRSTid(ng), rstVid(idFsur,ng),        &
     &                 tRSTindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, 1, 1, scale,                 &
#ifdef MASKING
     &                 GRID(ng) % rmask(LBi,LBj),                       &
#endif
     &                 OCEAN(ng) % zeta(LBi,LBj,KOUT))
      IF (OutThread.and.(status.ne.nf_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idFsur)), tRSTindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out 2D momentum component (m/s) in the XI-direction.
!
      scale=1.0_r8
      gtype=gfactor*u2dvar
      status=nf_fwrite(ng, iNLM, ncRSTid(ng), rstVid(idUbar,ng),        &
     &                 tRSTindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, 1, 1, scale,                 &
#ifdef MASKING
     &                 GRID(ng) % umask(LBi,LBj),                       &
#endif
     &                 OCEAN(ng) % ubar(LBi,LBj,KOUT))
      IF (OutThread.and.(status.ne.nf_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idUbar)), tRSTindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out 2D momentum component (m/s) in the ETA-direction.
!
      scale=1.0_r8
      gtype=gfactor*v2dvar
      status=nf_fwrite(ng, iNLM, ncRSTid(ng), rstVid(idVbar,ng),        &
     &                 tRSTindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, 1, 1, scale,                 &
#ifdef MASKING
     &                 GRID(ng) % vmask(LBi,LBj),                       &
#endif
     &                 OCEAN(ng) % vbar(LBi,LBj,KOUT))
      IF (OutThread.and.(status.ne.nf_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idVbar)), tRSTindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
#ifdef SOLVE3D
!
!  Write out 3D momentum component (m/s) in the XI-direction.
!
      scale=1.0_r8
      gtype=gfactor*u3dvar
      status=nf_fwrite(ng, iNLM, ncRSTid(ng), rstVid(idUvel,ng),        &
     &                 tRSTindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, 1, N(ng), scale,             &
# ifdef MASKING
     &                 GRID(ng) % umask(LBi,LBj),                       &
# endif
     &                 OCEAN(ng) % u(LBi,LBj,1,NOUT))
      IF (OutThread.and.(status.ne.nf_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idUvel)), tRSTindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out momentum component (m/s) in the ETA-direction.
!
      scale=1.0_r8
      gtype=gfactor*v3dvar
      status=nf_fwrite(ng, iNLM, ncRSTid(ng), rstVid(idVvel,ng),        &
     &                 tRSTindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, 1, N(ng), scale,             &
# ifdef MASKING
     &                 GRID(ng) % vmask(LBi,LBj),                       &
# endif
     &                 OCEAN(ng) % v(LBi,LBj,1,NOUT))
      IF (OutThread.and.(status.ne.nf_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idVvel)), tRSTindx(ng)
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
        status=nf_fwrite(ng, iNLM, ncRSTid(ng), rstTid(itrc,ng),        &
     &                   tRSTindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
# ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
# endif
     &                   OCEAN(ng) % t(LBi,LBj,1,NOUT,itrc))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idTvar(itrc))), tRSTindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END DO
!
!  Write out density anomaly.
!
      scale=1.0_r8
      gtype=gfactor*r3dvar
      status=nf_fwrite(ng, iNLM, ncRSTid(ng), rstVid(idDano,ng),        &
     &                 tRSTindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, 1, N(ng), scale,             &
# ifdef MASKING
     &                 GRID(ng) % rmask(LBi,LBj),                       &
# endif
     &                 OCEAN(ng) % rho(LBi,LBj,1))
      IF (OutThread.and.(status.ne.nf_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idDano)), tRSTindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
# ifdef LMD_SKPP
!
!  Write out depth of surface boundary layer.
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf_fwrite(ng, iNLM, ncRSTid(ng), rstVid(idHsbl,ng),        &
     &                 tRSTindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, 1, 1, scale,                 &
#  ifdef MASKING
     &                 GRID(ng) % rmask(LBi,LBj),                       &
#  endif
     &                 MIXING(ng) % hsbl(LBi,LBj))
      IF (OutThread.and.(status.ne.nf_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idHsbl)), tRSTindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
# endif
# ifdef LMD_BKPP
!
!  Write out depth of bottom boundary layer.
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf_fwrite(ng, iNLM, ncRSTid(ng), rstVid(idHbbl,ng),        &
     &                 tRSTindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, 1, 1, scale,                 &
#  ifdef MASKING
     &                 GRID(ng) % rmask(LBi,LBj),                       &
#  endif
     &                 MIXING(ng) % hbbl(LBi,LBj))
      IF (OutThread.and.(status.ne.nf_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idHbbl)), tRSTindx(ng)
        exit_flag=3
        ioerror=status
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
      status=nf_fwrite(ng, iNLM, ncRSTid(ng), rstVid(idVvis,ng),        &
     &                 tRSTindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, 0, N(ng), scale,             &
#  ifdef MASKING
     &                 GRID(ng) % rmask(LBi,LBj),                       &
#  endif
     &                 MIXING(ng) % Akv(LBi,LBj,0))
      IF (OutThread.and.(status.ne.nf_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idVvis)), tRSTindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out vertical diffusion coefficient for potential temperature.
!
      scale=1.0_r8
      gtype=gfactor*r3dvar
      status=nf_fwrite(ng, iNLM, ncRSTid(ng), rstVid(idTdif,ng),        &
     &                 tRSTindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, 0, N(ng), scale,             &
#  ifdef MASKING
     &                 GRID(ng) % rmask(LBi,LBj),                       &
#  endif
     &                 MIXING(ng) % Akt(LBi,LBj,0,itemp))
      IF (OutThread.and.(status.ne.nf_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idTdif)), tRSTindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
#  ifdef SALINITY
!
!  Write out vertical diffusion coefficient for salinity.
!
      scale=1.0_r8
      gtype=gfactor*r3dvar
      status=nf_fwrite(ng, iNLM, ncRSTid(ng), rstVid(idSdif,ng),        &
     &                 tRSTindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, 0, N(ng), scale,             &
#   ifdef MASKING
     &                 GRID(ng) % rmask(LBi,LBj),                       &
#   endif
     &                 MIXING(ng) % Akt(LBi,LBj,0,isalt))
      IF (OutThread.and.(status.ne.nf_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idSdif)), tRSTindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
#  endif
# endif
# ifdef SEDIMENT
#  ifdef BEDLOAD
!
!  Write out bed load transport in U-direction.
!
      DO i=1,NST
          scale=1.0_r8
          gtype=gfactor*u2dvar
          status=nf_fwrite(ng, iNLM, ncRSTid(ng), rstVid(idUbld(i),ng), &
     &                     tRSTindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, 1,  scale,            &
#   ifdef MASKING
     &                     GRID(ng) % umask(LBi,LBj),                   &
#   endif
     &                     OCEAN(ng) % bedldu(LBi,LBj,i))
          IF (OutThread.and.(status.ne.nf_noerr)) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbld(i))), tRSTindx(ng)
            exit_flag=3
            ioerror=status
            RETURN
          END IF
      END DO
!
!  Write out bed load transport in V-direction.
!
      DO i=1,NST
          scale=1.0_r8
          gtype=gfactor*v2dvar
          status=nf_fwrite(ng, iNLM, ncRSTid(ng), rstVid(idVbld(i),ng), &
     &                     tRSTindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, 1,  scale,            &
#   ifdef MASKING
     &                     GRID(ng) % vmask(LBi,LBj),                   &
#   endif
     &                     OCEAN(ng) % bedldv(LBi,LBj,i))
          IF (OutThread.and.(status.ne.nf_noerr)) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbld(i))), tRSTindx(ng)
            exit_flag=3
            ioerror=status
            RETURN
          END IF
      END DO
#  endif
!
!  Write out sediment fraction of each size class in each bed layer.
!
      DO i=1,NST
        scale=1.0_r8
        gtype=gfactor*b3dvar
        status=nf_fwrite(ng, iNLM, ncRSTid(ng), rstVid(idfrac(i),ng),   &
     &                   tRSTindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, Nbed, scale,            &
#  ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#  endif
     &                   OCEAN(ng) % bed_frac(LBi,LBj,1,i))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idfrac(i))), tRSTindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END DO
!
!  Write out sediment mass of each size class in each bed layer.
!
      DO i=1,NST
        scale=1.0_r8
        gtype=gfactor*b3dvar
        status=nf_fwrite(ng, iNLM, ncRSTid(ng), rstVid(idBmas(i),ng),   &
     &                   tRSTindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, Nbed, scale,            &
#  ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#  endif
     &                   OCEAN(ng) % bed_mass(LBi,LBj,1,NOUT,i))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idBmas(i))), tRSTindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END DO
!
!  Write out sediment properties in each bed layer.
!
      DO i=1,MBEDP
        scale=1.0_r8
        gtype=gfactor*b3dvar
        status=nf_fwrite(ng, iNLM, ncRSTid(ng), rstVid(idSbed(i),ng),   &
     &                   tRSTindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, Nbed, scale,            &
#  ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#  endif
     &                   OCEAN(ng) % bed(LBi,LBj,1,i))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idSbed(i))), tRSTindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END DO
# endif
# if defined SEDIMENT || defined BBL_MODEL
!
!  Write out exposed sediment layer properties. Notice that only the
!  first four properties (mean grain diameter, mean grain density,
!  mean settling velocity, mean critical erosion stress, 
!  ripple length and ripple height) are written.
!
      DO i=1,6
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite(ng, iNLM, ncRSTid(ng), rstVid(idBott(i),ng),   &
     &                   tRSTindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
#  ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#  endif
     &                   OCEAN(ng) % bottom(LBi,LBj,i))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idBott(i))), tRSTindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END DO
# endif
#endif
!
!  Synchronize history NetCDF file to disk.
!
      IF (OutThread) THEN
        status=nf_sync(ncRSTid(ng))
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,20)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#ifdef SOLVE3D
        WRITE (stdout,30) KOUT, NOUT, tRSTindx
#else
        WRITE (stdout,30) KOUT, tRSTindx
#endif
      END IF
!
  10  FORMAT (/,' WRT_RST - error while writing variable: ',a,/,11x,    &
     &        'into restart NetCDF file for time record: ',i4)
  20  FORMAT (/,' WRT_RST - unable to synchronize restart NetCDF to ',  &
     &        'disk.')
#ifdef SOLVE3D
  30  FORMAT (6x,'WRT_RST   - wrote re-start fields (Index=', i1,       &
     &        ',',i1,') into time record = ',i7.7)
#else
  30  FORMAT (6x,'WRT_RST   - wrote re-start fields (Index=', i1,       &
     &        ') into time record = ',i7.7)
#endif
      RETURN
      END SUBROUTINE wrt_rst
