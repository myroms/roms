#include "cppdefs.h"
#ifdef FOUR_DVAR
      SUBROUTINE wrt_ini (ng,Tindex)
!
!=======================================================================
!  Copyright (c) 2005 ROMS/TOMS Adjoint Group                          !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine writes state variables initial conditions into initial !
!  NetCDF file.                                                        !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     Tindex     State variables time index to write.                  !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
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
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Tindex
!
!  Local variable declarations.
!
      integer :: LBi, UBi, LBj, UBj
      integer :: gfactor, gtype, i, itrc, status, varid

      integer :: nf_fwrite

      integer, dimension(2) :: start, total

      real(r8) :: scale
!
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)
!
!-----------------------------------------------------------------------
!  Write out initial conditions.
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
      tINIindx(ng)=tINIindx(ng)+1
      NrecINI(ng)=NrecINI(ng)+1

# ifdef S4DVAR
!
!  If requested, set time index to recycle time record in the nonlinear
!  initial file.
!
      IF (LcycleINI(ng)) THEN
        tINIindx(ng)=MOD(tINIindx(ng)-1,2)+1
      END IF
# endif
!
!  Write out model time (s).
!
      IF (OutThread) THEN
        status=nf_put_var1_TYPE(ncINIid(ng), iniVid(idtime,ng),         &
     &                          tINIindx(ng), time(ng))
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,10) TRIM(Vname(1,idtime)), tINIindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write model-observations misfit cost function value.
!
      IF ((Ipass.eq.1).and.OutThread) THEN
        status=nf_inq_varid(ncINIid(ng),'cost_function',varid)
        IF (status.eq.nf_noerr) THEN
          start(1)=1
          total(1)=NstateVar(ng)+1
          start(2)=Nrun
          total(2)=1
          status=nf_put_vara_TYPE(ncINIid(ng),varid,start,total,        &
     &                            FOURDVAR(ng)%CostFun(0:))
          IF (status.ne.nf_noerr) THEN
            WRITE (stdout,10) 'cost_function', Nrun
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END IF
!
!  Write out free-surface (m)
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf_fwrite(ng, iNLM, ncINIid(ng), iniVid(idFsur,ng),        &
     &                 tINIindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, 1, 1, scale,                 &
# ifdef MASKING
     &                 GRID(ng) % rmask(LBi,LBj),                       &
# endif
     &                 OCEAN(ng) % zeta(LBi,LBj,Tindex))
      IF (OutThread.and.(status.ne.nf_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idFsur)), tINIindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out 2D momentum component (m/s) in the XI-direction.
!
      scale=1.0_r8
      gtype=gfactor*u2dvar
      status=nf_fwrite(ng, iNLM, ncINIid(ng), iniVid(idUbar,ng),        &
     &                 tINIindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, 1, 1, scale,                 &
# ifdef MASKING
     &                 GRID(ng) % umask(LBi,LBj),                       &
# endif
     &                 OCEAN(ng) % ubar(LBi,LBj,Tindex))
      IF (OutThread.and.(status.ne.nf_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idUbar)), tINIindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out 2D momentum component (m/s) in the ETA-direction.
!
      scale=1.0_r8
      gtype=gfactor*v2dvar
      status=nf_fwrite(ng, iNLM, ncINIid(ng), iniVid(idVbar,ng),        &
     &                 tINIindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, 1, 1, scale,                 &
# ifdef MASKING
     &                 GRID(ng) % vmask(LBi,LBj),                       &
# endif
     &                 OCEAN(ng) % vbar(LBi,LBj,Tindex))
      IF (OutThread.and.(status.ne.nf_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idVbar)), tINIindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
# ifdef SOLVE3D
!
!  Write out 3D momentum component (m/s) in the XI-direction.
!
      scale=1.0_r8
      gtype=gfactor*u3dvar
      status=nf_fwrite(ng, iNLM, ncINIid(ng), iniVid(idUvel,ng),        &
     &                 tINIindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, 1, N(ng), scale,             &
#  ifdef MASKING
     &                 GRID(ng) % umask(LBi,LBj),                       &
#  endif
     &                 OCEAN(ng) % u(LBi,LBj,1,Tindex))
      IF (OutThread.and.(status.ne.nf_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idUvel)), tINIindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out 3D momentum component (m/s) in the ETA-direction.
!
      scale=1.0_r8
      gtype=gfactor*v3dvar
      status=nf_fwrite(ng, iNLM, ncINIid(ng), iniVid(idVvel,ng),        &
     &                 tINIindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, 1, N(ng), scale,             &
#  ifdef MASKING
     &                 GRID(ng) % vmask(LBi,LBj),                       &
#  endif
     &                 OCEAN(ng) % v(LBi,LBj,1,Tindex))
      IF (OutThread.and.(status.ne.nf_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idVvel)), tINIindx(ng)
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
        status=nf_fwrite(ng, iNLM, ncINIid(ng), iniTid(itrc,ng),        &
     &                   tINIindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#  ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
#  endif
     &                   OCEAN(ng) % t(LBi,LBj,1,Tindex,itrc))
        IF (OutThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idTvar(itrc))), tINIindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END DO
# endif
!
!-------------------------------------------------------------------------
!  Synchronize initial NetCDF file to disk to allow other processes
!  to access data immediately after it is written.
!-------------------------------------------------------------------------
!
      IF (OutThread) THEN
        status=nf_sync(ncINIid(ng))
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,20)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
# ifdef SOLVE3D
        WRITE (stdout,30) Tindex, Tindex, tINIindx(ng)
# else
        WRITE (stdout,30) Tindex, tINIindx(ng)
# endif
      END IF
!
  10  FORMAT (/,' WRT_INI - error while writing variable: ',a,/,11x,    &
     &        'into initial NetCDF file for time record: ',i4)
  20  FORMAT (/,' WRT_INI - unable to synchronize initial NetCDF to ',  &
     &        'disk.')
# ifdef SOLVE3D
  30  FORMAT (6x,'WRT_INI   - wrote initial  fields (Index=', i1,       &
     &        ',',i1,') into time record = ',i7.7)
# else
  30  FORMAT (6x,'WRT_INI   - wrote initial  fields (Index=', i1,       &
     &        ') into time record = ',i7.7)
# endif
      RETURN
      END SUBROUTINE wrt_ini
#else
      SUBROUTINE wrt_ini
      RETURN
      END SUBROUTINE wrt_ini
#endif
