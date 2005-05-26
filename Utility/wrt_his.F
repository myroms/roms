#include "cppdefs.h"
      SUBROUTINE wrt_his (ng)
!
!=======================================================================
!  Copyright (c) 2005 ROMS/TOMS Group                                  !
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
#ifdef SOLVE3D
      USE mod_coupling
#endif
      USE mod_forces
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
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)
!
!-----------------------------------------------------------------------
!  Write out history fields.
!-----------------------------------------------------------------------
!
      if (exit_flag.ne.NoError) RETURN
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
      tHISindx(ng)=tHISindx(ng)+1
      NrecHIS(ng)=NrecHIS(ng)+1
!
!  Write out model time (s).
!
      IF (OutThread) THEN
        status=nf_put_var1_TYPE(ncHISid(ng), hisVid(idtime,ng),         &
     &                          tHISindx(ng),time(ng))
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,10) TRIM(Vname(1,idtime)), tHISindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out free-surface (m)
!
      IF (Hout(idFsur,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite(ng, iNLM, ncHISid(ng), hisVid(idFsur,ng),      &
     &                   tHISindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
#ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#endif
     &                   OCEAN(ng) % zeta(LBi,LBj,KOUT))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idFsur)), tHISindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#if defined FORWARD_WRITE && defined FORWARD_RHS
        status=nf_fwrite(ng, iNLM, ncHISid(ng), hisVid(idRzet,ng),      &
     &                   tHISindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
# ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
# endif
     &                   OCEAN(ng) % rzeta(LBi,LBj,KOUT))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idRzet)), tHISindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#endif
      END IF
!
!  Write out 2D momentum component (m/s) in the XI-direction.
!
      IF (Hout(idUbar,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u2dvar
        status=nf_fwrite(ng, iNLM, ncHISid(ng), hisVid(idUbar,ng),      &
     &                   tHISindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
#ifdef MASKING
     &                   GRID(ng) % umask(LBi,LBj),                     &
#endif
     &                   OCEAN(ng) % ubar(LBi,LBj,KOUT))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUbar)), tHISindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#ifdef FORWARD_WRITE
# ifdef FORWARD_RHS
        status=nf_fwrite(ng, iNLM, ncHISid(ng), hisVid(idRu2d,ng),      &
     &                   tHISindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
#  ifdef MASKING
     &                   GRID(ng) % umask(LBi,LBj),                     &
#  endif
     &                   OCEAN(ng) % rubar(LBi,LBj,KOUT))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idRu2d)), tHISindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
# endif
# ifdef SOLVE3D
#  ifdef FORWARD_RHS
        status=nf_fwrite(ng, iNLM, ncHISid(ng), hisVid(idRuct,ng),      &
     &                   tHISindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
#   ifdef MASKING
     &                   GRID(ng) % umask(LBi,LBj),                     &
#   endif
     &                   COUPLING(ng) % rufrc(LBi,LBj))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idRuct)), tHISindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#  endif
        status=nf_fwrite(ng, iNLM, ncHISid(ng), hisVid(idUfx1,ng),      &
     &                   tHISindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
#  ifdef MASKING
     &                   GRID(ng) % umask(LBi,LBj),                     &
#  endif
     &                   COUPLING(ng) % DU_avg1(LBi,LBj))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUfx1)), tHISindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        status=nf_fwrite(ng, iNLM, ncHISid(ng), hisVid(idUfx2,ng),      &
     &                   tHISindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
#  ifdef MASKING
     &                   GRID(ng) % umask(LBi,LBj),                     &
#  endif
     &                   COUPLING(ng) % DU_avg2(LBi,LBj))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUfx2)), tHISindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
# endif
#endif
      END IF
!
!  Write out 2D momentum component (m/s) in the ETA-direction.
!
      IF (Hout(idVbar,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v2dvar
        status=nf_fwrite(ng, iNLM, ncHISid(ng), hisVid(idVbar,ng),      &
     &                   tHISindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
#ifdef MASKING
     &                   GRID(ng) % vmask(LBi,LBj),                     &
#endif
     &                   OCEAN(ng) % vbar(LBi,LBj,KOUT))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVbar)), tHISindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#ifdef FORWARD_WRITE
# ifdef FORWARD_RHS
        status=nf_fwrite(ng, iNLM, ncHISid(ng), hisVid(idRv2d,ng),      &
     &                   tHISindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
#  ifdef MASKING
     &                   GRID(ng) % vmask(LBi,LBj),                     &
#  endif
     &                   OCEAN(ng) % rvbar(LBi,LBj,KOUT))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idRv2d)), tHISindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
# endif
# ifdef SOLVE3D
#  ifdef FORWARD_RHS
        status=nf_fwrite(ng, iNLM, ncHISid(ng), hisVid(idRvct,ng),      &
     &                   tHISindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
#   ifdef MASKING
     &                   GRID(ng) % vmask(LBi,LBj),                     &
#   endif
     &                   COUPLING(ng) % rvfrc(LBi,LBj))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idRvct)), tHISindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#  endif
        status=nf_fwrite(ng, iNLM, ncHISid(ng), hisVid(idVfx1,ng),      &
     &                   tHISindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
#  ifdef MASKING
     &                   GRID(ng) % vmask(LBi,LBj),                     &
#  endif
     &                   COUPLING(ng) % DV_avg1(LBi,LBj))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVfx1)), tHISindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        status=nf_fwrite(ng, iNLM, ncHISid(ng), hisVid(idVfx2,ng),      &
     &                   tHISindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
#  ifdef MASKING
     &                   GRID(ng) % vmask(LBi,LBj),                     &
#  endif
     &                   COUPLING(ng) % DV_avg2(LBi,LBj))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVfx2)), tHISindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
# endif
#endif
      END IF
#ifdef SOLVE3D
!
!  Write out 3D momentum component (m/s) in the XI-direction.
!
      IF (Hout(idUvel,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u3dvar
        status=nf_fwrite(ng, iNLM, ncHISid(ng), hisVid(idUvel,ng),      &
     &                   tHISindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
# ifdef MASKING
     &                   GRID(ng) % umask(LBi,LBj),                     &
# endif
     &                   OCEAN(ng) % u(LBi,LBj,1,NOUT))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUvel)), tHISindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
# if defined FORWARD_WRITE && defined FORWARD_RHS
        status=nf_fwrite(ng, iNLM, ncHISid(ng), hisVid(idRu3d,ng),      &
     &                   tHISindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#  ifdef MASKING
     &                   GRID(ng) % umask(LBi,LBj),                     &
#  endif
     &                   OCEAN(ng) % ru(LBi,LBj,1,NOUT))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idRu3d)), tHISindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
# endif
      END IF
!
!  Write out 3D momentum component (m/s) in the ETA-direction.
!
      IF (Hout(idVvel,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v3dvar
        status=nf_fwrite(ng, iNLM, ncHISid(ng), hisVid(idVvel,ng),      &
     &                   tHISindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
# ifdef MASKING
     &                   GRID(ng) % vmask(LBi,LBj),                     &
# endif
     &                   OCEAN(ng) % v(LBi,LBj,1,NOUT))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVvel)), tHISindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
# if defined FORWARD_WRITE && defined FORWARD_RHS
        status=nf_fwrite(ng, iNLM, ncHISid(ng), hisVid(idRv3d,ng),      &
     &                   tHISindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#  ifdef MASKING
     &                   GRID(ng) % vmask(LBi,LBj),                     &
#  endif
     &                   OCEAN(ng) % rv(LBi,LBj,1,NOUT))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idRv3d)), tHISindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
# endif
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
        status=nf_fwrite(ng, iNLM, ncHISid(ng), hisVid(idOvel,ng),      &
     &                   tHISindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 0, N(ng), scale,           &
# ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
# endif
     &                   wrk(LBi,LBj,0))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idOvel)), tHISindx(ng)
          exit_flag=3
          ioerror=status
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
        status=nf_fwrite(ng, iNLM, ncHISid(ng), hisVid(idWvel,ng),      &
     &                   tHISindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 0, N(ng), scale,           &
# ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
# endif
     &                   OCEAN(ng) % wvel(LBi,LBj,0))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idWvel)), tHISindx(ng)
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
          status=nf_fwrite(ng, iNLM, ncHISid(ng), hisTid(itrc,ng),      &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
# ifdef MASKING
     &                     GRID(ng) % rmask(LBi,LBj),                   &
# endif
     &                     OCEAN(ng) % t(LBi,LBj,1,NOUT,itrc))
          IF (OutThread.and.(status.ne.nf_noerr)) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTvar(itrc))),              &
     &                        tHISindx(ng)
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
        status=nf_fwrite(ng, iNLM, ncHISid(ng), hisVid(idDano,ng),      &
     &                   tHISindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
# ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
# endif
     &                   OCEAN(ng) % rho(LBi,LBj,1))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idDano)), tHISindx(ng)
          exit_flag=3
          ioerror=status
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
        status=nf_fwrite(ng, iNLM, ncHISid(ng), hisVid(idHsbl,ng),      &
     &                   tHISindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
#  ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#  endif
     &                   MIXING(ng) % hsbl(LBi,LBj))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idHsbl)), tHISindx(ng)
          exit_flag=3
          ioerror=status
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
        status=nf_fwrite(ng, iNLM, ncHISid(ng), hisVid(idHbbl,ng),      &
     &                   tHISindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
#  ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#  endif
     &                   MIXING(ng) % hbbl(LBi,LBj))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idHbbl)), tHISindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
# if defined FORWARD_WRITE && defined LMD_NONLOCAL
!
!  Write out KPP nonlocal transport.
!
      DO i=1,NAT
        scale=1.0_r8
        gtype=gfactor*w3dvar
        status=nf_fwrite(ng, iNLM, ncHISid(ng), hisVid(idGhat(i),ng),   &
     &                   tHISindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 0, N(ng), scale,           &
#  ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#  endif
     &                   MIXING(ng) % ghats(LBi,LBj,0))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idGhat(i))), tHISindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END DO
# endif
!
!  Write out vertical viscosity coefficient.
!
      IF (Hout(idVvis,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*w3dvar
        status=nf_fwrite(ng, iNLM, ncHISid(ng), hisVid(idVvis,ng),      &
     &                   tHISindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 0, N(ng), scale,           &
# ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
# endif
     &                   MIXING(ng) % Akv(LBi,LBj,0))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVvis)), tHISindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out vertical diffusion coefficient for potential temperature.
!
      IF (Hout(idTdif,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*w3dvar
        status=nf_fwrite(ng, iNLM, ncHISid(ng), hisVid(idTdif,ng),      &
     &                   tHISindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 0, N(ng), scale,           &
# ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
# endif
     &                   MIXING(ng) % Akt(LBi,LBj,0,itemp))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idTdif)), tHISindx(ng)
          exit_flag=3
          ioerror=status
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
        status=nf_fwrite(ng, iNLM, ncHISid(ng), hisVid(idSdif,ng),      &
     &                   tHISindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 0, N(ng), scale,           &
#  ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#  endif
     &                   MIXING(ng) % Akt(LBi,LBj,0,isalt))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idSdif)), tHISindx(ng)
          exit_flag=3
          ioerror=status
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
        status=nf_fwrite(ng, iNLM, ncHISid(ng), hisVid(idMtke,ng),      &
     &                   tHISindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 0, N(ng), scale,           &
#  ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#  endif
     &                   MIXING(ng) % tke(LBi,LBj,0,NOUT))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idMtke)), tHISindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out turbulent length scale field.
!
      IF (Hout(idMtls,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*w3dvar
        status=nf_fwrite(ng, iNLM, ncHISid(ng), hisVid(idMtls,ng),      &
     &                   tHISindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 0, N(ng), scale,           &
#  ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#  endif
     &                   MIXING(ng) % gls(LBi,LBj,0,NOUT))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idMtls)), tHISindx(ng)
          exit_flag=3
          ioerror=status
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
        status=nf_fwrite(ng, iNLM, ncHISid(ng),                         &
     &                   hisVid(idTsur(itemp),ng),                      &
     &                   tHISindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
# ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
# endif
     &                   FORCES(ng) % stflx(LBi,LBj,itemp))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idTsur(itemp))),               &
     &                      tHISindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface net salt flux.
!
      IF (Hout(idTsur(isalt),ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite(ng, iNLM, ncHISid(ng),                         &
     &                   hisVid(idTsur(isalt),ng),                      &
     &                   tHISindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
# ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
# endif
     &                   FORCES(ng) % stflx(LBi,LBj,isalt))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idTsur(isalt))),               &
     &                      tHISindx(ng)
          exit_flag=3
          ioerror=status
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
        status=nf_fwrite(ng, iNLM, ncHISid(ng), hisVid(idLhea,ng),      &
     &                   tHISindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
#  ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#  endif
     &                   FORCES(ng) % lhflx(LBi,LBj))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idLhea)), tHISindx(ng)
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
        status=nf_fwrite(ng, iNLM, ncHISid(ng), hisVid(idShea,ng),      &
     &                   tHISindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
#  ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#  endif
     &                   FORCES(ng) % shflx(LBi,LBj))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idShea)), tHISindx(ng)
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
        status=nf_fwrite(ng, iNLM, ncHISid(ng), hisVid(idLrad,ng),      &
     &                   tHISindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
#  ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#  endif
     &                   FORCES(ng) % lrflx(LBi,LBj))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idLrad)), tHISindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  ifdef EMINUSP
!
!  Write out evaporation rate (kg/m2/s).
!
      IF (Hout(idevap,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite(ng, iNLM, ncHISid(ng), hisVid(idevap,ng),      &
     &                   tHISindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
#   ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#   endif
     &                   FORCES(ng) % evap(LBi,LBj))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idevap)), tHISindx(ng)
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
        status=nf_fwrite(ng, iNLM, ncHISid(ng), hisVid(idrain,ng),      &
     &                   tHISindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
#   ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#   endif
     &                   FORCES(ng) % rain(LBi,LBj))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idrain)), tHISindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
# endif
# ifdef SHORTWAVE
!
!  Write out shortwave radiation flux.
!
      IF (Hout(idSrad,ng)) THEN
        scale=rho0*Cp
        gtype=gfactor*r2dvar
        status=nf_fwrite(ng, iNLM, ncHISid(ng), hisVid(idSrad,ng),      &
     &                   tHISindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
#  ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#  endif
     &                   FORCES(ng) % srflx(LBi,LBj))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idSrad)), tHISindx(ng)
          exit_flag=3
          ioerror=status
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
        status=nf_fwrite(ng, iNLM, ncHISid(ng), hisVid(idUsms,ng),      &
     &                   tHISindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
#ifdef MASKING
     &                   GRID(ng) % umask(LBi,LBj),                     &
#endif
     &                   FORCES(ng) % sustr(LBi,LBj))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUsms)), tHISindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface V-momentum stress.
!
      IF (Hout(idVsms,ng)) THEN
        scale=rho0
        gtype=gfactor*v2dvar
        status=nf_fwrite(ng, iNLM, ncHISid(ng), hisVid(idVsms,ng),      &
     &                   tHISindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
#ifdef MASKING
     &                   GRID(ng) % vmask(LBi,LBj),                     &
#endif
     &                   FORCES(ng) % svstr(LBi,LBj))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVsms)), tHISindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom U-momentum stress.
!
      IF (Hout(idUbms,ng)) THEN
        scale=-rho0
        gtype=gfactor*u2dvar
        status=nf_fwrite(ng, iNLM, ncHISid(ng), hisVid(idUbms,ng),      &
     &                   tHISindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
#ifdef MASKING
     &                   GRID(ng) % umask(LBi,LBj),                     &
#endif
     &                   FORCES(ng) % bustr(LBi,LBj))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUbms)), tHISindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom V-momentum stress.
!
      IF (Hout(idVbms,ng)) THEN
        scale=-rho0
        gtype=gfactor*v2dvar
        status=nf_fwrite(ng, iNLM, ncHISid(ng), hisVid(idVbms,ng),      &
     &                   tHISindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
#ifdef MASKING
     &                   GRID(ng) % vmask(LBi,LBj),                     &
#endif
     &                   FORCES(ng) % bvstr(LBi,LBj))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVbms)), tHISindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#ifdef SOLVE3D
# ifdef BBL_MODEL
!
!  Write out current-induced, bottom U-stress at RHO-points.
!
      IF (Hout(idUbrs,ng)) THEN
        scale=-rho0
        gtype=gfactor*r2dvar
        status=nf_fwrite(ng, iNLM, ncHISid(ng), hisVid(idUbrs,ng),      &
     &                   tHISindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
#  ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#  endif
     &                   BBL(ng) % bustrc(LBi,LBj))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUbrs)), tHISindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out current-induced, bottom V-stress at RHO-points.
!
      IF (Hout(idVbrs,ng)) THEN
        scale=-rho0
        gtype=gfactor*r2dvar
        status=nf_fwrite(ng, iNLM, ncHISid(ng), hisVid(idVbrs,ng),      &
     &                   tHISindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
#  ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#  endif
     &                   BBL(ng) % bvstrc(LBi,LBj))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVbrs)), tHISindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out wind-induced, bottom U-stress at RHO-points.
!
      IF (Hout(idUbws,ng)) THEN
        scale=rho0
        gtype=gfactor*r2dvar
        status=nf_fwrite(ng, iNLM, ncHISid(ng), hisVid(idUbws,ng),      &
     &                   tHISindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
#  ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#  endif
     &                   BBL(ng) % bustrw(LBi,LBj))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUbws)), tHISindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out wind-induced, bottom V-stress at RHO-points.
!
      IF (Hout(idVbws,ng)) THEN
        scale=rho0
        gtype=gfactor*r2dvar
        status=nf_fwrite(ng, iNLM, ncHISid(ng), hisVid(idVbws,ng),      &
     &                   tHISindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
#  ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#  endif
     &                   BBL(ng) % bvstrw(LBi,LBj))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVbws)), tHISindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out maximum wind and current, bottom U-stress at RHO-points.
!
      IF (Hout(idUbcs,ng)) THEN
        scale=rho0
        gtype=gfactor*r2dvar
        status=nf_fwrite(ng, iNLM, ncHISid(ng), hisVid(idUbcs,ng),      &
     &                   tHISindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
#  ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#  endif
     &                   BBL(ng) % bustrcwmax(LBi,LBj))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUbcs)), tHISindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out maximum wind and current, bottom V-stress at RHO-points.
!
      IF (Hout(idVbcs,ng)) THEN
        scale=rho0
        gtype=gfactor*r2dvar
        status=nf_fwrite(ng, iNLM, ncHISid(ng), hisVid(idVbcs,ng),      &
     &                   tHISindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
#  ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#  endif
     &                   BBL(ng) % bvstrcwmax(LBi,LBj))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVbcs)), tHISindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out wind-induced, bed wave orbital U-velocity at RHO-points.
!
      IF (Hout(idUbot,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite(ng, iNLM, ncHISid(ng), hisVid(idUbot,ng),      &
     &                   tHISindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
#  ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#  endif
     &                   BBL(ng) % Ubot(LBi,LBj))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUbot)), tHISindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out wind-induced, bed wave orbital V-velocity at RHO-points
!
      IF (Hout(idVbot,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite(ng, iNLM, ncHISid(ng), hisVid(idVbot,ng),      &
     &                   tHISindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
#  ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#  endif
     &                   BBL(ng) % Vbot(LBi,LBj))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVbot)), tHISindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom U-velocity above bed at RHO-points.
!
      IF (Hout(idUbur,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite(ng, iNLM, ncHISid(ng), hisVid(idUbur,ng),      &
     &                   tHISindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
#  ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#  endif
     &                   BBL(ng) % Ur(LBi,LBj))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUbur)), tHISindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom V-velocity above bed at RHO-points.
!
      IF (Hout(idVbvr,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite(ng, iNLM, ncHISid(ng), hisVid(idVbvr,ng),      &
     &                   tHISindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, 1, scale,               &
#  ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#  endif
     &                   BBL(ng) % Vr(LBi,LBj))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVbvr)), tHISindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
# ifdef SEDIMENT
#  ifdef BEDLOAD
!
!  Write out bed load transport in U-direction.
!
      DO i=1,NST
        IF (Hout(idUbld(i),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*u2dvar
          status=nf_fwrite(ng, iNLM, ncHISid(ng), hisVid(idUbld(i),ng), &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, 1,  scale,            &
#   ifdef MASKING
     &                     GRID(ng) % umask(LBi,LBj),                   &
#   endif
     &                     OCEAN(ng) % bedldu(LBi,LBj,i))
          IF (OutThread.and.(status.ne.nf_noerr)) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbld(i))), tHISindx(ng)
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
!
!  Write out bed load transport in V-direction.
!
        IF (Hout(idVbld(i),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*v2dvar
          status=nf_fwrite(ng, iNLM, ncHISid(ng), hisVid(idVbld(i),ng), &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, 1,  scale,            &
#   ifdef MASKING
     &                     GRID(ng) % vmask(LBi,LBj),                   &
#   endif
     &                     OCEAN(ng) % bedldv(LBi,LBj,i))
          IF (OutThread.and.(status.ne.nf_noerr)) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbld(i))), tHISindx(ng)
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
#  endif
!
!  Write out sediment fraction of each size class in each bed layer.
!
      DO i=1,NST
        IF (Hout(idfrac(i),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*b3dvar
          status=nf_fwrite(ng, iNLM, ncHISid(ng), hisVid(idfrac(i),ng), &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, Nbed, scale,          &
#  ifdef MASKING
     &                     GRID(ng) % rmask(LBi,LBj),                   &
#  endif
     &                     OCEAN(ng) % bed_frac(LBi,LBj,1,i))
          IF (OutThread.and.(status.ne.nf_noerr)) THEN
            WRITE (stdout,10) TRIM(Vname(1,idfrac(i))), tHISindx(ng)
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
!
!  Write out sediment mass of each size class in each bed layer.
!
      DO i=1,NST
        IF (Hout(idBmas(i),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*b3dvar
          status=nf_fwrite(ng, iNLM, ncHISid(ng), hisVid(idBmas(i),ng), &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, Nbed, scale,          &
#  ifdef MASKING
     &                     GRID(ng) % rmask(LBi,LBj),                   &
#  endif
     &                     OCEAN(ng) % bed_mass(LBi,LBj,1,NOUT,i))
          IF (OutThread.and.(status.ne.nf_noerr)) THEN
            WRITE (stdout,10) TRIM(Vname(1,idBmas(i))), tHISindx(ng)
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
!
!  Write out sediment properties in each bed layer.
!
      DO i=1,MBEDP
        IF (Hout(idSbed(i),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*b3dvar
          status=nf_fwrite(ng, iNLM, ncHISid(ng), hisVid(idSbed(i),ng), &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, Nbed, scale,          &
#  ifdef MASKING
     &                     GRID(ng) % rmask(LBi,LBj),                   &
#  endif
     &                     OCEAN(ng) % bed(LBi,LBj,1,i))
          IF (OutThread.and.(status.ne.nf_noerr)) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbed(i))), tHISindx(ng)
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
# endif
# if defined SEDIMENT || defined BBL_MODEL
!
!  Write out exposed sediment layer properties.
!
      DO i=1,MBOTP
        IF (Hout(idBott(i),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*r2dvar
          status=nf_fwrite(ng, iNLM, ncHISid(ng), hisVid(idBott(i),ng), &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, 1, scale,             &
#  ifdef MASKING
     &                     GRID(ng) % rmask(LBi,LBj),                   &
#  endif
     &                     OCEAN(ng) % bottom(LBi,LBj,i))
          IF (OutThread.and.(status.ne.nf_noerr)) THEN
            WRITE (stdout,10) TRIM(Vname(1,idBott(i))), tHISindx(ng)
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
# endif
#endif
!
!  Synchronize history NetCDF file to disk to allow other processes
!  to access data immediately after it is written.
!
      IF (OutThread) THEN
        status=nf_sync(ncHISid(ng))
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,20)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#ifdef SOLVE3D
        WRITE (stdout,30) KOUT, NOUT, tHISindx
#else
        WRITE (stdout,30) KOUT, tHISindx
#endif
      END IF
!
  10  FORMAT (/,' WRT_HIS - error while writing variable: ',a,/,11x,    &
     &        'into history NetCDF file for time record: ',i4)
  20  FORMAT (/,' WRT_HIS - unable to synchronize history NetCDF to ',  &
     &        'disk.')
#ifdef SOLVE3D
  30  FORMAT (6x,'WRT_HIS   - wrote history  fields (Index=', i1,       &
     &        ',',i1,') into time record = ',i7.7)
#else
  30  FORMAT (6x,'WRT_HIS   - wrote history  fields (Index=', i1,       &
     &        ') into time record = ',i7.7)
#endif
      RETURN
      END SUBROUTINE wrt_his
