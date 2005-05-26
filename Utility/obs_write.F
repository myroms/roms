#include "cppdefs.h"
#ifdef FOUR_DVAR
      SUBROUTINE obs_write (ng, model)
!
!=======================================================================
!  Copyright (c) 2005 ROMS/TOMS adjoint group                          !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This subroutine interpolates nonlinear and/or tangent linear model  !
!  state at observation locations when appropriate.  Then,  it writes  !
!  model values and  interpolation  weights  matrix into observations  !
!  NetCDF file.                                                        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_ocean
      USE mod_scalars
      USE mod_stepping
!  
      USE extract_obs_mod, ONLY : extract_obs2d
# ifdef SOLVE3D
      USE extract_obs_mod, ONLY : extract_obs3d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
!
!  Local variable declarations.
!
      logical :: load_Zobs = .FALSE.

      integer, dimension(2) ::  start, total

      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: LBi, UBi, LBj, UBj
      integer :: Istr, Iend, Jstr, Jend
# ifdef DISTRIBUTE
      integer :: Itile, Jtile
# endif
      integer :: Irec, iobs, itrc, iweight, status

      real(r8), parameter :: IniVal = 0.0_r8
!
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)
!
!-----------------------------------------------------------------------
!  Interpolate model state at observation locations.
!-----------------------------------------------------------------------
!
      IF (ProcessObs(ng)) THEN
# if defined S4DVAR
        Irec=Ipass
# else
        Irec=1
# endif
!
!  Get physical non-overlapping horizontal tile bounds.
!
# ifdef DISTRIBUTE
        CALL get_tile (ng, TILE, Itile, Jtile, Istr, Iend, Jstr, Jend)
# else
        Istr=1
        Iend=Lm(ng)
        Jstr=1
        Jend=Mm(ng)
# endif
# include "set_bounds.h"
!
!  Initialize interpolation weights matrix.  Some entries are not
!  computed in the extraction routine.  Set values to zero to avoid
!  problems when writing non initialized values.
!
        DO iobs=1,Mobs
          DO iweight=1,Nweights
            Hmat(iweight,iobs)=IniVal
          END DO
        END DO
        IF (wrtNLmod(ng)) THEN
          DO iobs=1,Mobs
            NLmodVal(iobs,Irec)=IniVal
          END DO
        END IF
# ifdef TLM_OBS
        IF (wrtTLmod(ng)) THEN
          DO iobs=1,Mobs
            TLmodVal(iobs,Ipass)=IniVal
          END DO
        END IF
# endif
!
!  Free-surface observations.
!
        IF (wrtNLmod(ng).and.                                           &
     &      (FOURDVAR(ng)%ObsCount(isFsur).gt.0)) THEN
          CALL extract_obs2d (ng, IstrR, IendR, JstrR, JendR,           &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        rILB(ng), rIUB(ng),                       &
     &                        rJLB(ng), rJUB(ng),                       &
     &                        isFsur, Mobs, Nobs(ng), Nweights,         &
     &                        ObsType, Xobs, Yobs,                      &
     &                        OCEAN(ng)%zeta(:,:,KOUT),                 &
# ifdef MASKING
     &                        GRID(ng)%rmask,                           &
# endif
     &                        NLmodVal(:,Irec), Hmat)
        END IF
# ifdef TLM_OBS
        IF (wrtTLmod(ng).and.                                           &
     &      (FOURDVAR(ng)%ObsCount(isFsur).gt.0)) THEN
          CALL extract_obs2d (ng, IstrR, IendR, JstrR, JendR,           &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        rILB(ng), rIUB(ng),                       &
     &                        rJLB(ng), rJUB(ng),                       &
     &                        isFsur, Mobs, Nobs(ng), Nweights,         &
     &                        ObsType, Xobs, Yobs,                      &
     &                        OCEAN(ng)%tl_zeta(:,:,KOUT),              &
#  ifdef MASKING
     &                        GRID(ng)%rmask,                           &
#  endif
     &                        TLmodVal(:,Ipass), Hmat)
        END IF
# endif
!
!  Vertically integrated u-velocity observations.
!
        IF (wrtNLmod(ng).and.                                           &
     &      (FOURDVAR(ng)%ObsCount(isUbar).gt.0)) THEN
          CALL extract_obs2d (ng, Istr, IendR, JstrR, JendR,            &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        uILB(ng), uIUB(ng),                       &
     &                        uJLB(ng), uJUB(ng),                       &
     &                        isUbar, Mobs, Nobs(ng), Nweights,         &
     &                        ObsType, Xobs, Yobs,                      &
     &                        OCEAN(ng)%ubar(:,:,KOUT),                 &
# ifdef MASKING
     &                        GRID(ng)%umask,                           &
# endif
     &                        NLmodVal(:,Irec), Hmat)
        END IF
# ifdef TLM_OBS
        IF (wrtTLmod(ng).and.                                           &
     &      (FOURDVAR(ng)%ObsCount(isUbar).gt.0)) THEN
          CALL extract_obs2d (ng, Istr, IendR, JstrR, JendR,            &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        uILB(ng), uIUB(ng),                       &
     &                        uJLB(ng), uJUB(ng),                       &
     &                        isUbar, Mobs, Nobs(ng), Nweights,         &
     &                        ObsType, Xobs, Yobs,                      &
     &                        OCEAN(ng)%tl_ubar(:,:,KOUT),              &
#  ifdef MASKING
     &                        GRID(ng)%umask,                           &
#  endif
     &                        TLmodVal(:,Ipass), Hmat)
        END IF
# endif
!
!  Vertically integrated v-velocity observations.
!
        IF (wrtNLmod(ng).and.                                           &
     &      (FOURDVAR(ng)%ObsCount(isVbar).gt.0)) THEN
          CALL extract_obs2d (ng, IstrR, IendR, Jstr, JendR,            &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        vILB(ng), vIUB(ng),                       &
     &                        vJLB(ng), vJUB(ng),                       &
     &                        isVbar, Mobs, Nobs(ng), Nweights,         &
     &                        ObsType, Xobs, Yobs,                      &
     &                        OCEAN(ng)%vbar(:,:,KOUT),                 &
# ifdef MASKING
     &                        GRID(ng)%vmask,                           &
# endif
     &                        NLmodVal(:,Irec), Hmat)
        END IF
# ifdef TLM_OBS
        IF (wrtTLmod(ng).and.                                           &
     &      (FOURDVAR(ng)%ObsCount(isVbar).gt.0)) THEN
          CALL extract_obs2d (ng, IstrR, IendR, Jstr, JendR,            &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        vILB(ng), vIUB(ng),                       &
     &                        vJLB(ng), vJUB(ng),                       &
     &                        isVbar, Mobs, Nobs(ng), Nweights,         &
     &                        ObsType, Xobs, Yobs,                      &
     &                        OCEAN(ng)%tl_vbar(:,:,KOUT),              &
#  ifdef MASKING
     &                        GRID(ng)%vmask,                           &
#  endif
     &                        TLmodVal(:,Ipass), Hmat)
        END IF
# endif

# ifdef SOLVE3D
!
!  3D u-velocity observations.
!
        IF (wrtNLmod(ng).and.                                           &
     &      (FOURDVAR(ng)%ObsCount(isUvel).gt.0)) THEN
          load_Zobs=.TRUE.
          CALL extract_obs3d (ng, Istr, IendR, JstrR, JendR,            &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        uILB(ng), uIUB(ng),                       &
     &                        uJLB(ng), uJUB(ng),                       &
     &                        isUvel, Mobs, Nobs(ng), Nweights,         &
     &                        ObsType, Xobs, Yobs, Zobs,                &
     &                        OCEAN(ng)%u(:,:,:,NOUT),                  &
     &                        GRID(ng)%z_r,                             &
#  ifdef MASKING
     &                        GRID(ng)%umask,                           &
#  endif
     &                        NLmodVal(:,Irec), Hmat)
        END IF
#  ifdef TLM_OBS
        IF (wrtTLmod(ng).and.                                           &
     &      (FOURDVAR(ng)%ObsCount(isUvel).gt.0)) THEN
          load_Zobs=.TRUE.
          CALL extract_obs3d (ng, Istr, IendR, JstrR, JendR,            &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        uILB(ng), uIUB(ng),                       &
     &                        uJLB(ng), uJUB(ng),                       &
     &                        isUvel, Mobs, Nobs(ng), Nweights,         &
     &                        ObsType, Xobs, Yobs, Zobs,                &
     &                        OCEAN(ng)%tl_u(:,:,:,NOUT),               &
     &                        GRID(ng)%z_r,                             &
#   ifdef MASKING
     &                        GRID(ng)%umask,                           &
#   endif
     &                        TLmodVal(:,Ipass), Hmat)
        END IF
#  endif
!
!  3D v-velocity observations.
!
        IF (wrtNLmod(ng).and.                                           &
     &      (FOURDVAR(ng)%ObsCount(isVvel).gt.0)) THEN
          load_Zobs=.TRUE.
          CALL extract_obs3d (ng, IstrR, IendR, Jstr, JendR,            &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        vILB(ng), vIUB(ng),                       &
     &                        vJLB(ng), vJUB(ng),                       &
     &                        isVvel, Mobs, Nobs(ng), Nweights,         &
     &                        ObsType, Xobs, Yobs, Zobs,                &
     &                        OCEAN(ng)%v(:,:,:,NOUT),                  &
     &                        GRID(ng)%z_r,                             &
#  ifdef MASKING
     &                        GRID(ng)%vmask,                           &
#  endif
     &                        NLmodVal(:,Irec), Hmat)
        END IF
#  ifdef TLM_OBS
        IF (wrtTLmod(ng).and.                                           &
     &      (FOURDVAR(ng)%ObsCount(isVvel).gt.0)) THEN
          load_Zobs=.TRUE.
          CALL extract_obs3d (ng, IstrR, IendR, Jstr, JendR,            &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        vILB(ng), vIUB(ng),                       &
     &                        vJLB(ng), vJUB(ng),                       &
     &                        isVvel, Mobs, Nobs(ng), Nweights,         &
     &                        ObsType, Xobs, Yobs, Zobs,                &
     &                        OCEAN(ng)%tl_v(:,:,:,NOUT),               &
     &                        GRID(ng)%z_r,                             &
#   ifdef MASKING
     &                        GRID(ng)%vmask,                           &
#   endif
     &                        TLmodVal(:,Ipass), Hmat)
        END IF
#  endif
!
!  Tracer type observations.
!
        DO itrc=1,NT(ng)
          IF (wrtNLmod(ng).and.                                         &
     &        (FOURDVAR(ng)%ObsCount(isTvar(itrc)).gt.0)) THEN
            load_Zobs=.TRUE.
            CALL extract_obs3d (ng, IstrR, IendR, JstrR, JendR,         &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          rILB(ng), rIUB(ng),                     &
     &                          rJLB(ng), rJUB(ng),                     &
     &                          isTvar(itrc), Mobs, Nobs(ng), Nweights, &
     &                          ObsType, Xobs, Yobs, Zobs,              &
     &                          OCEAN(ng)%t(:,:,:,NOUT,itrc),           &
     &                          GRID(ng)%z_r,                           &
#  ifdef MASKING
     &                          GRID(ng)%rmask,                         &
#  endif
     &                          NLmodVal(:,Irec), Hmat)
          END IF
#  ifdef TLM_OBS
          IF (wrtTLmod(ng).and.                                         &
     &        (FOURDVAR(ng)%ObsCount(isTvar(itrc)).gt.0)) THEN
            load_Zobs=.TRUE.
            CALL extract_obs3d (ng, IstrR, IendR, JstrR, JendR,         &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          rILB(ng), rIUB(ng),                     &
     &                          rJLB(ng), rJUB(ng),                     &
     &                          isTvar(itrc), Mobs, Nobs(ng), Nweights, &
     &                          ObsType, Xobs, Yobs, Zobs,              &
     &                          OCEAN(ng)%tl_t(:,:,:,NOUT,itrc),        &
     &                          GRID(ng)%z_r,                           &
#   ifdef MASKING
     &                          GRID(ng)%rmask,                         &
#   endif
     &                          TLmodVal(:,Ipass), Hmat)
          END IF
#  endif
        END DO
# endif

# ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Collect all extracted data.
!-----------------------------------------------------------------------
!
        IF (wrtNLmod(ng)) THEN 
          CALL mp_collect (ng, model, Nobs(ng), IniVal,                 &
     &                     NLmodVal(1,Irec))
        END IF
#  ifdef TLM_OBS
        IF (wrtTLmod(ng)) THEN 
          CALL mp_collect (ng, model, Nobs(ng), IniVal,                 &
     &                     TLmodVal(1,Ipass))
#  endif
        END IF
!!      IF (load_Zobs) THEN                ! HGA: I need figure out this
!!        CALL mp_collect (ng, model, Nobs(ng) IniVal, Zobs) 
!!      END IF
        CALL mp_collect (ng, model, Mobs*Nweights, IniVal, Hmat)
# endif

!
!-----------------------------------------------------------------------
!  Write out data into observations NetCDF file.
!-----------------------------------------------------------------------
!
!  Nonlinear model state at observation locations.
!
        IF (OutThread.and.wrtNLmod(ng)) THEN 
          start(1)=Irec
          total(1)=1
          start(2)=NstrObs(ng)
          total(2)=Nobs(ng)
          status=nf_put_vara_TYPE(ncOBSid(ng),obsVid(idNLmo,ng),        &
     &                            start,total,NLmodVal(1,Irec))
          IF (status.ne.nf_noerr) then
            WRITE (stdout,10) TRIM(Vname(1,idNLmo)), TRIM(OBSname(ng))
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF

# ifdef TLM_OBS
!
!  Tangent linear model state at observation locations.
!
        IF (OutThread.and.wrtTLmod(ng)) THEN 
          start(1)=Ipass
          total(1)=1
          start(2)=NstrObs(ng)
          total(2)=Nobs(ng)
          status=nf_put_vara_TYPE(ncOBSid(ng),obsVid(idTLmo,ng),        &
     &                            start,total,TLmodVal(1,Ipass))
          IF (status.ne.nf_noerr) then
            WRITE (stdout,10) TRIM(Vname(1,idTLmo)), TRIM(OBSname(ng))
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
# endif
!
!  Interpolation weights matrix.
!
        IF (OutThread.and.wrtNLmod(ng)) THEN 
          start(1)=1
          total(1)=Nweights
          start(2)=NstrObs(ng)
          total(2)=Nobs(ng)
          status=nf_put_vara_TYPE(ncOBSid(ng),obsVid(idModH,ng),        &
     &                            start,total,Hmat)
          IF (status.ne.nf_noerr) then
            WRITE (stdout,10) TRIM(Vname(1,idModH)), TRIM(OBSname(ng))
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
# ifdef SOLVE3D
!
!  Write Z-location of observation in grid coordinates, if applicable.
!  This values are needed elsewhere when using the interpolation
!  weight matrix.  Recall that the depth of observations can be in
!  meters or grid coordinates.  Recall that since the model levels
!  evolve in time, the fractional level coordinate is unknow during
!  the processing of the observations.
!  
        IF (OutThread.and.load_Zobs) THEN
          start(1)=NstrObs(ng)
          total(1)=Nobs(ng)
          status=nf_put_vara_TYPE(ncOBSid(ng),obsVid(idObsZ,ng),        &
     &                            start,total,Zobs)
          IF (status.ne.nf_noerr) then
            WRITE (stdout,10) TRIM(Vname(1,idObsZ)), TRIM(OBSname(ng))
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
# endif
        IF (Master) THEN
          WRITE (stdout,30) NstrObs(ng), NendObs(ng)
        END IF
!
!-----------------------------------------------------------------------
!  Synchronize observations NetCDF file to disk.
!-----------------------------------------------------------------------
!
        IF (OutThread) THEN
          status=nf_sync(ncOBSid(ng))
          IF (status.ne.nf_noerr) THEN
            WRITE (stdout,20)
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END IF
!
  10  FORMAT (/,' OBS_WRITE - error while writing variable: ',a,/,      &
     &        11x,'into observation NetCDF: ',a)
  20  FORMAT (/,' OBS_WRITE - unable to synchronize stations',          &
     &        1x,'NetCDF file to disk.')
  30  FORMAT (' Wrote observation fields from records = ',              &
     &        i7.7,' - ',i7.7)

      RETURN
      END SUBROUTINE obs_write
#else
      SUBROUTINE obs_write
      RETURN
      END SUBROUTINE obs_write
#endif
