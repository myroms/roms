#include "cppdefs.h"
#ifndef ANA_GRID
      SUBROUTINE get_grid (ng, model)
!
!=======================================================================
!  Copyright (c) 2005 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This subroutine reads grid information from GRID NetCDF file.       !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
!
# if (defined EW_PERIODIC || defined NS_PERIODIC) && !defined DISTRIBUTE
      USE exchange_2d_mod
!
# endif

      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
!
!  Local variable declarations.
!
      logical :: gotang = .FALSE.
      logical :: gotel = .FALSE.
      logical :: gotf = .FALSE.
      logical :: gotdm = .FALSE.
      logical :: gotdn = .FALSE.
      logical :: goth = .FALSE.
      logical :: gotlatp = .FALSE.
      logical :: gotlatr = .FALSE.
      logical :: gotlatu = .FALSE.
      logical :: gotlatv = .FALSE.
      logical :: gotlonp = .FALSE.
      logical :: gotlonr = .FALSE.
      logical :: gotlonu = .FALSE.
      logical :: gotlonv = .FALSE.
# ifdef MASKING
      logical :: gotmskp = .FALSE.
      logical :: gotmskr = .FALSE.
      logical :: gotmsku = .FALSE.
      logical :: gotmskv = .FALSE.
# endif
      logical :: gotpm = .FALSE.
      logical :: gotpn = .FALSE.
      logical :: gotsph = .FALSE.
      logical :: gotxl = .FALSE.
      logical :: gotxp = .FALSE.
      logical :: gotxr = .FALSE.
      logical :: gotxu = .FALSE.
      logical :: gotxv = .FALSE.
      logical :: gotyp = .FALSE.
      logical :: gotyr = .FALSE.
      logical :: gotyu = .FALSE.
      logical :: gotyv = .FALSE.
# ifdef ICESHELF
      logical :: gotzice = .FALSE.

      integer :: grdziceid
# endif
      integer :: LBi, UBi, LBj, UBj
# ifdef DISTRIBUTE
      integer :: Itile, Jtile, Istr, Iend, Jstr, Jend
# endif
      integer :: grdangid, grdelid, grdfid, grddmid, grddnid
      integer :: grdhid, grdpmid, grdpnid
      integer :: grdlatpid, grdlatrid, grdlatuid, grdlatvid
      integer :: grdlonpid, grdlonrid, grdlonuid, grdlonvid
# ifdef MASKING
      integer :: mskpid, mskrid, mskuid, mskvid
# endif
      integer :: grdsphid, grdxlid
      integer :: grdxpid, grdxrid, grdxuid, grdxvid
      integer :: grdypid, grdyrid, grdyuid, grdyvid
      integer :: i, ifield, nrec, status

      integer :: nf_fread

      real(r8), parameter :: Fscl = 1.0_r8

      real(r8) :: Fmax, Fmin

      character (len=1 ) :: char1
      character (len=80) :: fname, ncname
!
!-----------------------------------------------------------------------
!  Inquire about the contents of grid NetCDF file:  Inquire about
!  the dimensions and variables.  Check for consistency.
!-----------------------------------------------------------------------
!
      IF (exit_flag.ne.NoError) RETURN
      ncname=GRDname(ng)
      ifield=0
      CALL opencdf (ng, 1, ncname, fname, N(ng), ifield, nrec)
      IF (exit_flag.ne.NoError) RETURN
!
!  Scan variable list from input NetCDF and check for grid variables.
!
      DO i=1,nvars
        IF (TRIM(varnam(i)).eq.'xl') THEN
          grdxlid=i
          gotxl=.TRUE.
        ELSE IF (TRIM(varnam(i)).eq.'el') THEN
          grdelid=i
          gotel=.TRUE.
        ELSE IF (TRIM(varnam(i)).eq.'spherical') THEN
          grdsphid=i
          gotsph=.TRUE.
# ifdef MASKING
        ELSE IF (TRIM(varnam(i)).eq.'mask_rho') THEN
          mskrid=i
          gotmskr=.TRUE.
        ELSE IF (TRIM(varnam(i)).eq.'mask_u') THEN
          mskuid=i
          gotmsku=.TRUE.
        ELSE IF (TRIM(varnam(i)).eq.'mask_v') THEN
          mskvid=i
          gotmskv=.TRUE.
        ELSE IF (TRIM(varnam(i)).eq.'mask_psi') THEN
          mskpid=i
          gotmskp=.TRUE.
# endif
        ELSE IF (TRIM(varnam(i)).eq.'h') THEN
          grdhid=i
          goth=.TRUE.
# ifdef ICESHELF
        ELSE IF (TRIM(varnam(i)).eq.'zice') THEN
          grdziceid=i
          gotzice=.TRUE.
# endif
        ELSE IF (TRIM(varnam(i)).eq.'f') THEN
          grdfid=i
          gotf=.TRUE.
        ELSE IF (TRIM(varnam(i)).eq.'pm') THEN
          grdpmid=i
          gotpm=.TRUE.
        ELSE IF (TRIM(varnam(i)).eq.'pn') THEN
          grdpnid=i
          gotpn=.TRUE.
        ELSE IF (TRIM(varnam(i)).eq.'dndx') THEN
          grddnid=i
          gotdn=.TRUE.
        ELSE IF (TRIM(varnam(i)).eq.'dmde') THEN
          grddmid=i
          gotdm=.TRUE.
        ELSE IF (TRIM(varnam(i)).eq.'x_psi') THEN
          grdxpid=i
          gotxp=.TRUE.
        ELSE IF (TRIM(varnam(i)).eq.'x_rho') THEN
          grdxrid=i
          gotxr=.TRUE.
        ELSE IF (TRIM(varnam(i)).eq.'x_u') THEN
          grdxuid=i
          gotxu=.TRUE.
        ELSE IF (TRIM(varnam(i)).eq.'x_v') THEN
          grdxvid=i
          gotxv=.TRUE.
        ELSE IF (TRIM(varnam(i)).eq.'y_psi') THEN
          grdypid=i
          gotyp=.TRUE.
        ELSE IF (TRIM(varnam(i)).eq.'y_rho') THEN
          grdyrid=i
          gotyr=.TRUE.
        ELSE IF (TRIM(varnam(i)).eq.'y_u') THEN
          grdyuid=i
          gotyu=.TRUE.
        ELSE IF (TRIM(varnam(i)).eq.'y_v') THEN
          grdyvid=i
          gotyv=.TRUE.
        ELSE IF (TRIM(varnam(i)).eq.'lon_psi') THEN
          grdlonpid=i
          gotlonp=.TRUE.
        ELSE IF (TRIM(varnam(i)).eq.'lon_rho') THEN
          grdlonrid=i
          gotlonr=.TRUE.
        ELSE IF (TRIM(varnam(i)).eq.'lon_u') THEN
          grdlonuid=i
          gotlonu=.TRUE.
        ELSE IF (TRIM(varnam(i)).eq.'lon_v') THEN
          grdlonvid=i
          gotlonv=.TRUE.
        ELSE IF (TRIM(varnam(i)).eq.'lat_psi') THEN
          grdlatpid=i
          gotlatp=.TRUE.
        ELSE IF (TRIM(varnam(i)).eq.'lat_rho') THEN
          grdlatrid=i
          gotlatr=.TRUE.
        ELSE IF (TRIM(varnam(i)).eq.'lat_u') THEN
          grdlatuid=i
          gotlatu=.TRUE.
        ELSE IF (TRIM(varnam(i)).eq.'lat_v') THEN
          grdlatvid=i
          gotlatv=.TRUE.
        ELSE IF (TRIM(varnam(i)).eq.'angle') THEN
          grdangid=i
          gotang=.TRUE.
        END IF
      END DO
!
!  Terminate execution if essential grid variables are not found.
!
      IF (.not.gotxl) THEN
        WRITE (stdout,10) 'xl', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.gotel) THEN
        WRITE (stdout,10) 'el', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.gotsph) THEN
        WRITE (stdout,10) 'spherical', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
# ifdef MASKING
      IF (.not.gotmskr) THEN
        WRITE (stdout,10) 'mask_rho', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.gotmsku) THEN
        WRITE (stdout,10) 'mask_u', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.gotmskv) THEN
        WRITE (stdout,10) 'mask_v', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.gotmskp) THEN
        WRITE (stdout,10) 'mask_psi', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
# endif
      IF (.not.goth) THEN
        WRITE (stdout,10) 'h', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
# ifdef ICESHELF
      IF (.not.gotzice) THEN
        WRITE (stdout,10) 'zice', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
# endif
      IF (.not.gotf) THEN
        WRITE (stdout,10) 'f', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.gotpm) THEN
        WRITE (stdout,10) 'pm', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.gotpn) THEN
        WRITE (stdout,10) 'pn', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
# if (defined CURVGRID && defined UV_ADV)
      IF (.not.gotdn) THEN
        WRITE (stdout,10) 'dndx', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.gotdm) THEN
        WRITE (stdout,10) 'dmde', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
# endif
# ifdef CURVGRID
      IF (.not.gotang) THEN
        WRITE (stdout,10) 'angle', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
# endif
!
!  Open grid NetCDF file for reading.
!
      IF (ncGRDid(ng).eq.-1) THEN
        status=nf_open(TRIM(ncname), nf_nowrite, ncGRDid(ng))
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,20) TRIM(ncname)
          exit_flag=2
          ioerror=status
          RETURN
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Read in grid parameters.
!-----------------------------------------------------------------------
!
!  Read in basin lengths.
!
      status=nf_get_var1_TYPE(ncGRDid(ng), grdxlid, 1, xl(ng))
      IF (status.ne.nf_noerr) THEN
        WRITE (stdout,30) 'xl', TRIM(ncname)
        exit_flag=2
        ioerror=status
        RETURN
      END IF
      status=nf_get_var1_TYPE(ncGRDid(ng), grdelid, 1, el(ng))
      IF (status.ne.nf_noerr) THEN
        WRITE (stdout,30) 'el', TRIM(ncname)
        exit_flag=2
        ioerror=status
        RETURN
      END IF
!
!  Read in logical switch for spherical grid configuration.
!
      status=nf_get_var1_text(ncGRDid(ng), grdsphid, 1, char1)
      IF (status.ne.nf_noerr) THEN
        WRITE (stdout,30) 'spherical', TRIM(ncname)
        exit_flag=2
        ioerror=status
        RETURN
      END IF
      IF ((char1.eq.'t').or.(char1.eq.'T')) THEN
        spherical=.TRUE.
      ELSE
        spherical=.FALSE.
      END IF
!
!-----------------------------------------------------------------------
!  Read in grid arrays.
!-----------------------------------------------------------------------
!
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)
# ifdef DISTRIBUTE
      CALL get_tile (ng, MyRank, Itile, Jtile, Istr, Iend, Jstr, Jend)
# endif

# ifdef MASKING
!
!  Read in mask at RHO-points.
!
      status=nf_fread(ng, model, ncGRDid(ng), mskrid, 0, r2dvar,        &
     &                LBi, UBi, LBj, UBj, 1, 1,                         &
     &                Fscl, Fmin, Fmax,                                 &
     &                GRID(ng) % rmask(LBi,LBj),                        &
     &                GRID(ng) % rmask(LBi,LBj))
      IF (InpThread.and.(status.ne.nf_noerr)) THEN
        WRITE (stdout,30) 'mask_rho', TRIM(ncname)
        exit_flag=2
        ioerror=status
        RETURN
      END IF
#  if (defined EW_PERIODIC || defined NS_PERIODIC) && !defined DISTRIBUTE
      CALL exchange_r2d_tile (ng, model, 1, Lm(ng), 1, Mm(ng),          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        GRID(ng) % rmask)
#  endif
!
!  Read in mask at U-points.
!
      status=nf_fread(ng, model, ncGRDid(ng), mskuid, 0, u2dvar,        &
     &                LBi, UBi, LBj, UBj, 1, 1,                         &
     &                Fscl, Fmin, Fmax,                                 &
     &                GRID(ng) % umask(LBi,LBj),                        &
     &                GRID(ng) % umask(LBi,LBj))
      IF (InpThread.and.(status.ne.nf_noerr)) THEN
        WRITE (stdout,30) 'mask_u', TRIM(ncname)
        exit_flag=2
        ioerror=status
        RETURN
      END IF
#  if (defined EW_PERIODIC || defined NS_PERIODIC) && !defined DISTRIBUTE
      CALL exchange_u2d_tile (ng, model, 1, Lm(ng), 1, Mm(ng),          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        GRID(ng) % umask)
#  endif
!
!  Read in mask at V-points.
!
      status=nf_fread(ng, model, ncGRDid(ng), mskvid, 0, v2dvar,        &
     &                LBi, UBi, LBj, UBj, 1, 1,                         &
     &                Fscl, Fmin, Fmax,                                 &
     &                GRID(ng) % vmask(LBi,LBj),                        &
     &                GRID(ng) % vmask(LBi,LBj))
      IF (InpThread.and.(status.ne.nf_noerr)) THEN
        WRITE (stdout,30) 'mask_v', TRIM(ncname)
        exit_flag=2
        ioerror=status
        RETURN
      END IF
#  if (defined EW_PERIODIC || defined NS_PERIODIC) && !defined DISTRIBUTE
      CALL exchange_v2d_tile (ng, model, 1, Lm(ng), 1, Mm(ng),          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        GRID(ng) % vmask)
#  endif
!
!  Read in mask at PSI-points.
!
      status=nf_fread(ng, model, ncGRDid(ng), mskpid, 0, p2dvar,        &
     &                LBi, UBi, LBj, UBj, 1, 1,                         &
     &                Fscl, Fmin, Fmax,                                 &
     &                GRID(ng) % pmask(LBi,LBj),                        &
     &                GRID(ng) % pmask(LBi,LBj))
      IF (InpThread.and.(status.ne.nf_noerr)) THEN
        WRITE (stdout,30) 'mask_psi', TRIM(ncname)
        exit_flag=2
        ioerror=status
        RETURN
      END IF
#  if (defined EW_PERIODIC || defined NS_PERIODIC) && !defined DISTRIBUTE
      CALL exchange_p2d_tile (ng, model, 1, Lm(ng), 1, Mm(ng),          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        GRID(ng) % pmask)
#  endif
# endif
!
!  Read in bathymetry.
!
      status=nf_fread(ng, model, ncGRDid(ng), grdhid, 0, r2dvar,        &
     &                LBi, UBi, LBj, UBj, 1, 1,                         &
     &                Fscl, hmin(ng), hmax(ng),                         &
# ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
# endif
     &                GRID(ng) % h(LBi,LBj))
      IF (InpThread.and.(status.ne.nf_noerr)) THEN
        WRITE (stdout,30) 'h', TRIM(ncname)
        exit_flag=2
        ioerror=status
        RETURN
      END IF
# if (defined EW_PERIODIC || defined NS_PERIODIC) && !defined DISTRIBUTE
      CALL exchange_r2d_tile (ng, model, 1, Lm(ng), 1, Mm(ng),          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        GRID(ng) % h)
# endif
# ifdef ICESHELF
!
!  Read in ice shelf thicknesses.
!
      status=nf_fread(ng, model, ncGRDid(ng), grdziceid, 0, r2dvar,     &
     &                LBi, UBi, LBj, UBj, 1, 1,                         &
     &                Fscl, Fmin, Fmax,                                 &
#  ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#  endif
     &                GRID(ng) % zice(LBi,LBj))
      IF (InpThread.and.(status.ne.nf_noerr)) THEN
        WRITE (stdout,30) 'zice', TRIM(ncname)
        exit_flag=2
        ioerror=status
        RETURN
      END IF
#  if (defined EW_PERIODIC || defined NS_PERIODIC) && !defined DISTRIBUTE
      CALL exchange_r2d_tile (ng, model, 1, Lm(ng), 1, Mm(ng),          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        GRID(ng) % zice)
#  endif
# endif
!
!  Read in Coriolis parameter.
!
      status=nf_fread(ng, model, ncGRDid(ng), grdfid, 0, r2dvar,        &
     &                LBi, UBi, LBj, UBj, 1, 1,                         &
     &                Fscl, Fmin, Fmax,                                 &
# ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
# endif
     &                GRID(ng) % f(LBi,LBj))
      IF (InpThread.and.(status.ne.nf_noerr)) THEN
        WRITE (stdout,30) 'f', TRIM(ncname)
        exit_flag=2
        ioerror=status
        RETURN
      END IF
# if (defined EW_PERIODIC || defined NS_PERIODIC) && !defined DISTRIBUTE
      CALL exchange_r2d_tile (ng, model, 1, Lm(ng), 1, Mm(ng),          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        GRID(ng) % f)
# endif
!
!  Read in coordinate transfomation metrics (m,n) associated with the
!  differential distances in XI and ETA.
!
      status=nf_fread(ng, model, ncGRDid(ng), grdpmid, 0, r2dvar,       &
     &                LBi, UBi, LBj, UBj, 1, 1,                         &
     &                Fscl, Fmin, Fmax,                                 &
# ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
# endif
     &                GRID(ng) % pm(LBi,LBj))
      IF (InpThread.and.(status.ne.nf_noerr)) THEN
        WRITE (stdout,30) 'pm', TRIM(ncname)
        exit_flag=2
        ioerror=status
        RETURN
      END IF
      status=nf_fread(ng, model, ncGRDid(ng), grdpnid, 0, r2dvar,       &
     &                LBi, UBi, LBj, UBj, 1, 1,                         &
     &                Fscl, Fmin, Fmax,                                 &
# ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
# endif
     &                GRID(ng) % pn(LBi,LBj))
      IF (InpThread.and.(status.ne.nf_noerr)) THEN
        WRITE (stdout,30) 'pn', TRIM(ncname)
        exit_flag=2
        ioerror=status
        RETURN
      END IF
# if (defined EW_PERIODIC || defined NS_PERIODIC) && !defined DISTRIBUTE
      CALL exchange_r2d_tile (ng, model, 1, Lm(ng), 1, Mm(ng),          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        GRID(ng) % pm)
      CALL exchange_r2d_tile (ng, model, 1, Lm(ng), 1, Mm(ng),          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        GRID(ng) % pn)
# endif
# if (defined CURVGRID && defined UV_ADV)
!
!  Read in derivatives of inverse metrics factors: rotational factors.
!
      status=nf_fread(ng, model, ncGRDid(ng), grddmid, 0, r2dvar,       &
     &                LBi, UBi, LBj, UBj, 1, 1,                         &
     &                Fscl, Fmin, Fmax,                                 &
# ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
# endif
     &                GRID(ng) % dmde(LBi,LBj))
      IF (InpThread.and.(status.ne.nf_noerr)) THEN
        WRITE (stdout,30) 'dmde', TRIM(ncname)
        exit_flag=2
        ioerror=status
        RETURN
      END IF
      status=nf_fread(ng, model, ncGRDid(ng), grddnid, 0, r2dvar,       &
     &                LBi, UBi, LBj, UBj, 1, 1,                         &
     &                Fscl, Fmin, Fmax,                                 &
# ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
# endif
     &                GRID(ng) % dndx(LBi,LBj))
      IF (InpThread.and.(status.ne.nf_noerr)) THEN
        WRITE (stdout,30) 'dndx', TRIM(ncname)
        exit_flag=2
        ioerror=status
        RETURN
      END IF
#  if (defined EW_PERIODIC || defined NS_PERIODIC) && !defined DISTRIBUTE
      CALL exchange_r2d_tile (ng, model, 1, Lm(ng), 1, Mm(ng),          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        GRID(ng) % dmde)
      CALL exchange_r2d_tile (ng, model, 1, Lm(ng), 1, Mm(ng),          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        GRID(ng) % dndx)
#  endif
# endif
!
!  Read in (x,y) coordinates at PSI-points.
!
      IF (gotxp) THEN
        status=nf_fread(ng, model, ncGRDid(ng), grdxpid, 0, p2dvar,     &
     &                  LBi, UBi, LBj, UBj, 1, 1,                       &
     &                  Fscl, Fmin, Fmax,                               &
# ifdef MASKING
     &                  GRID(ng) % pmask(LBi,LBj),                      &
# endif
     &                  GRID(ng) % xp(LBi,LBj))
        IF (InpThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,30) 'x_p', TRIM(ncname)
          exit_flag=2
          ioerror=status
          RETURN
        END IF
# ifdef DISTRIBUTE
        CALL mp_exchange (ng, model, Istr, Iend, Jstr, Jend,            &
     &                    LBi, UBi, LBj, UBj, 1, 1,                     &
     &                    NghostPoints, .FALSE., .FALSE.,               &
     &                    GRID(ng) % xp(LBi,LBj)) 
# endif
      END IF
      IF (gotyp) THEN
        status=nf_fread(ng, model,  ncGRDid(ng), grdypid, 0, p2dvar,    &
     &                  LBi, UBi, LBj, UBj, 1, 1,                       &
     &                  Fscl, Fmin, Fmax,                               &
# ifdef MASKING
     &                  GRID(ng) % pmask(LBi,LBj),                      &
# endif
     &                  GRID(ng) % yp(LBi,LBj))
        IF (InpThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,30) 'y_p', TRIM(ncname)
          exit_flag=2
          ioerror=status
          RETURN
        END IF
# ifdef DISTRIBUTE
        CALL mp_exchange (ng, model, Istr, Iend, Jstr, Jend,            &
     &                    LBi, UBi, LBj, UBj, 1, 1,                     &
     &                    NghostPoints, .FALSE., .FALSE.,               &
     &                    GRID(ng) % yp(LBi,LBj)) 
# endif
      END IF
!
!  Read in (x,y) coordinates at RHO-points.
!
      IF (gotxr) THEN
        status=nf_fread(ng, model, ncGRDid(ng), grdxrid, 0, r2dvar,     &
     &                  LBi, UBi, LBj, UBj, 1, 1,                       &
     &                  Fscl, Fmin, Fmax,                               &
# ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
# endif
     &                  GRID(ng) % xr(LBi,LBj))
        IF (InpThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,30) 'x_rho', TRIM(ncname)
          exit_flag=2
          ioerror=status
          RETURN
        END IF
# ifdef DISTRIBUTE
        CALL mp_exchange (ng, model, Istr, Iend, Jstr, Jend,            &
     &                    LBi, UBi, LBj, UBj, 1, 1,                     &
     &                    NghostPoints, .FALSE., .FALSE.,               &
     &                    GRID(ng) % xr(LBi,LBj)) 
# endif
      END IF
      IF (gotyr) THEN
        status=nf_fread(ng, model, ncGRDid(ng), grdyrid, 0, r2dvar,     &
     &                  LBi, UBi, LBj, UBj, 1, 1,                       &
     &                  Fscl, Fmin, Fmax,                               &
# ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
# endif
     &                  GRID(ng) % yr(LBi,LBj))
        IF (InpThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,30) 'y_rho', TRIM(ncname)
          exit_flag=2
          ioerror=status
          RETURN
        END IF
# ifdef DISTRIBUTE
        CALL mp_exchange (ng, model, Istr, Iend, Jstr, Jend,            &
     &                    LBi, UBi, LBj, UBj, 1, 1,                     &
     &                    NghostPoints, .FALSE., .FALSE.,               &
     &                    GRID(ng) % yr(LBi,LBj)) 
# endif
      END IF
!
!  Read in (x,y) coordinates at U-points.
!
      IF (gotxu) THEN
        status=nf_fread(ng, model, ncGRDid(ng), grdxuid, 0, u2dvar,     &
     &                  LBi, UBi, LBj, UBj, 1, 1,                       &
     &                  Fscl, Fmin, Fmax,                               &
# ifdef MASKING
     &                  GRID(ng) % umask(LBi,LBj),                      &
# endif
     &                  GRID(ng) % xu(LBi,LBj))
        IF (InpThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,30) 'x_u', TRIM(ncname)
          exit_flag=2
          ioerror=status
          RETURN
        END IF
# ifdef DISTRIBUTE
        CALL mp_exchange (ng, model, Istr, Iend, Jstr, Jend,            &
     &                    LBi, UBi, LBj, UBj, 1, 1,                     &
     &                    NghostPoints, .FALSE., .FALSE.,               &
     &                    GRID(ng) % xu(LBi,LBj)) 
# endif
      END IF
      IF (gotyu) THEN
        status=nf_fread(ng, model, ncGRDid(ng), grdyuid, 0, u2dvar,     &
     &                  LBi, UBi, LBj, UBj, 1, 1,                       &
     &                  Fscl, Fmin, Fmax,                               &
# ifdef MASKING
     &                  GRID(ng) % umask(LBi,LBj),                      &
# endif
     &                  GRID(ng) % yu(LBi,LBj))
        IF (InpThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,30) 'y_u', TRIM(ncname)
          exit_flag=2
          ioerror=status
          RETURN
        END IF
# ifdef DISTRIBUTE
        CALL mp_exchange (ng, model, Istr, Iend, Jstr, Jend,            &
     &                    LBi, UBi, LBj, UBj, 1, 1,                     &
     &                    NghostPoints, .FALSE., .FALSE.,               &
     &                    GRID(ng) % yu(LBi,LBj)) 
# endif
      END IF
!
!  Read in (x,y) coordinates at V-points.
!
      IF (gotxv) THEN
        status=nf_fread(ng, model, ncGRDid(ng), grdxvid, 0, v2dvar,     &
     &                  LBi, UBi, LBj, UBj, 1, 1,                       &
     &                  Fscl, Fmin, Fmax,                               &
# ifdef MASKING
     &                  GRID(ng) % vmask(LBi,LBj),                      &
# endif
     &                  GRID(ng) % xv(LBi,LBj))
        IF (InpThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,30) 'x_v', TRIM(ncname)
          exit_flag=2
          ioerror=status
          RETURN
        END IF
# ifdef DISTRIBUTE
        CALL mp_exchange (ng, model, Istr, Iend, Jstr, Jend,            &
     &                    LBi, UBi, LBj, UBj, 1, 1,                     &
     &                    NghostPoints, .FALSE., .FALSE.,               &
     &                    GRID(ng) % xv(LBi,LBj)) 
# endif
      END IF
      IF (gotyv) THEN
        status=nf_fread(ng, model, ncGRDid(ng), grdyvid, 0, v2dvar,     &
     &                  LBi, UBi, LBj, UBj, 1, 1,                       &
     &                  Fscl, Fmin, Fmax,                               &
# ifdef MASKING
     &                  GRID(ng) % vmask(LBi,LBj),                      &
# endif
     &                  GRID(ng) % yv(LBi,LBj))
        IF (InpThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,30) 'y_v', TRIM(ncname)
          exit_flag=2
          ioerror=status
          RETURN
        END IF
# ifdef DISTRIBUTE
        CALL mp_exchange (ng, model, Istr, Iend, Jstr, Jend,            &
     &                    LBi, UBi, LBj, UBj, 1, 1,                     &
     &                    NghostPoints, .FALSE., .FALSE.,               &
     &                    GRID(ng) % yv(LBi,LBj)) 
# endif
      END IF
!
!  Read in (lon,lat) coordinates at PSI-points.
!
      IF (spherical) THEN
        IF (gotlonp) THEN
          status=nf_fread(ng, model, ncGRDid(ng), grdlonpid, 0, p2dvar, &
     &                    LBi, UBi, LBj, UBj, 1, 1,                     &
     &                    Fscl, Fmin, Fmax,                             &
# ifdef MASKING
     &                    GRID(ng) % pmask(LBi,LBj),                    &
# endif
     &                    GRID(ng) % lonp(LBi,LBj))
          IF (InpThread.and.(status.ne.nf_noerr)) THEN
            WRITE (stdout,30) 'lon_psi', TRIM(ncname)
            exit_flag=2
            ioerror=status
            RETURN
          END IF
# ifdef DISTRIBUTE
          CALL mp_exchange (ng, model, Istr, Iend, Jstr, Jend,          &
     &                      LBi, UBi, LBj, UBj, 1, 1,                   &
     &                      NghostPoints, .FALSE., .FALSE.,             &
     &                      GRID(ng) % lonp(LBi,LBj)) 
# endif
        END IF
        IF (gotlatp) THEN
          status=nf_fread(ng, model, ncGRDid(ng), grdlatpid, 0, p2dvar, &
     &                    LBi, UBi, LBj, UBj, 1, 1,                     &
     &                    Fscl, Fmin, Fmax,                             &
# ifdef MASKING
     &                    GRID(ng) % pmask(LBi,LBj),                    &
# endif
     &                    GRID(ng) % latp(LBi,LBj))
          IF (InpThread.and.(status.ne.nf_noerr)) THEN
            WRITE (stdout,30) 'lat_psi', TRIM(ncname)
            exit_flag=2
            ioerror=status
            RETURN
          END IF
# ifdef DISTRIBUTE
          CALL mp_exchange (ng, model, Istr, Iend, Jstr, Jend,          &
     &                      LBi, UBi, LBj, UBj, 1, 1,                   &
     &                      NghostPoints, .FALSE., .FALSE.,             &
     &                      GRID(ng) % latp(LBi,LBj)) 
# endif
        END IF
      END IF
!
!  Read in (lon,lat) coordinates at RHO-points.
!
      IF (spherical) THEN
        IF (gotlonr) THEN
          status=nf_fread(ng, model, ncGRDid(ng), grdlonrid, 0, r2dvar, &
     &                    LBi, UBi, LBj, UBj, 1, 1,                     &
     &                    Fscl, Fmin, Fmax,                             &
# ifdef MASKING
     &                    GRID(ng) % rmask(LBi,LBj),                    &
# endif
     &                    GRID(ng) % lonr(LBi,LBj))
          IF (InpThread.and.(status.ne.nf_noerr)) THEN
            WRITE (stdout,30) 'lon_rho', TRIM(ncname)
            exit_flag=2
            ioerror=status
            RETURN
          END IF
# ifdef DISTRIBUTE
          CALL mp_exchange (ng, model, Istr, Iend, Jstr, Jend,          &
     &                      LBi, UBi, LBj, UBj, 1, 1,                   &
     &                      NghostPoints, .FALSE., .FALSE.,             &
     &                      GRID(ng) % lonr(LBi,LBj)) 
# endif
        END IF
        IF (gotlatr) THEN
          status=nf_fread(ng, model, ncGRDid(ng), grdlatrid, 0, r2dvar, &
     &                    LBi, UBi, LBj, UBj, 1, 1,                     &
     &                    Fscl, Fmin, Fmax,                             &
# ifdef MASKING
     &                    GRID(ng) % rmask(LBi,LBj),                    &
# endif
     &                    GRID(ng) % latr(LBi,LBj))
          IF (InpThread.and.(status.ne.nf_noerr)) THEN
            WRITE (stdout,30) 'lat_rho', TRIM(ncname)
            exit_flag=2
            ioerror=status
            RETURN
          END IF
# ifdef DISTRIBUTE
          CALL mp_exchange (ng, model, Istr, Iend, Jstr, Jend,          &
     &                      LBi, UBi, LBj, UBj, 1, 1,                   &
     &                      NghostPoints, .FALSE., .FALSE.,             &
     &                      GRID(ng) % latr(LBi,LBj)) 
# endif
        END IF
      END IF
!
!  Read in (lon,lat) coordinates at U-points.
!
      IF (spherical) THEN
        IF (gotlonu) THEN
          status=nf_fread(ng, model, ncGRDid(ng), grdlonuid, 0, u2dvar, &
     &                    LBi, UBi, LBj, UBj, 1, 1,                     &
     &                    Fscl, Fmin, Fmax,                             &
# ifdef MASKING
     &                    GRID(ng) % umask(LBi,LBj),                    &
# endif
     &                    GRID(ng) % lonu(LBi,LBj))
          IF (InpThread.and.(status.ne.nf_noerr)) THEN
            WRITE (stdout,30) 'lon_u', TRIM(ncname)
            exit_flag=2
            ioerror=status
            RETURN
          END IF
# ifdef DISTRIBUTE
          CALL mp_exchange (ng, model, Istr, Iend, Jstr, Jend,          &
     &                      LBi, UBi, LBj, UBj, 1, 1,                   &
     &                      NghostPoints, .FALSE., .FALSE.,             &
     &                      GRID(ng) % lonu(LBi,LBj)) 
# endif
        END IF
        IF (gotlatu) THEN
          status=nf_fread(ng, model, ncGRDid(ng), grdlatuid, 0, u2dvar, &
     &                    LBi, UBi, LBj, UBj, 1, 1,                     &
     &                    Fscl, Fmin, Fmax,                             &
# ifdef MASKING
     &                    GRID(ng) % umask(LBi,LBj),                    &
# endif
     &                    GRID(ng) % latu(LBi,LBj))
          IF (InpThread.and.(status.ne.nf_noerr)) THEN
            WRITE (stdout,30) 'lat_u', TRIM(ncname)
            exit_flag=2
            ioerror=status
            RETURN
          END IF
# ifdef DISTRIBUTE
          CALL mp_exchange (ng, model, Istr, Iend, Jstr, Jend,          &
     &                      LBi, UBi, LBj, UBj, 1, 1,                   &
     &                      NghostPoints, .FALSE., .FALSE.,             &
     &                      GRID(ng) % latu(LBi,LBj)) 
# endif
        END IF
      END IF
!
!  Read in (lon,lat) coordinates at V-points.
!
      IF (spherical) THEN
        IF (gotlonv) THEN
          status=nf_fread(ng, model, ncGRDid(ng), grdlonvid, 0, v2dvar, &
     &                    LBi, UBi, LBj, UBj, 1, 1,                     &
     &                    Fscl, Fmin, Fmax,                             &
# ifdef MASKING
     &                    GRID(ng) % vmask(LBi,LBj),                    &
# endif
     &                    GRID(ng) % lonv(LBi,LBj))
          IF (InpThread.and.(status.ne.nf_noerr)) THEN
            WRITE (stdout,30) 'lon_v', TRIM(ncname)
            exit_flag=2
            ioerror=status
            RETURN
          END IF
# ifdef DISTRIBUTE
          CALL mp_exchange (ng, model, Istr, Iend, Jstr, Jend,          &
     &                      LBi, UBi, LBj, UBj, 1, 1,                   &
     &                      NghostPoints, .FALSE., .FALSE.,             &
     &                      GRID(ng) % lonv(LBi,LBj)) 
# endif
        END IF
        IF (gotlatv) THEN
          status=nf_fread(ng, model, ncGRDid(ng), grdlatvid, 0, v2dvar, &
     &                    LBi, UBi, LBj, UBj, 1, 1,                     &
     &                    Fscl, Fmin, Fmax,                             &
# ifdef MASKING
     &                    GRID(ng) % vmask(LBi,LBj),                    &
# endif
     &                    GRID(ng) % latv(LBi,LBj))
          IF (InpThread.and.(status.ne.nf_noerr)) THEN
            WRITE (stdout,30) 'lat_v', TRIM(ncname)
            exit_flag=2
            ioerror=status
            RETURN
          END IF
# ifdef DISTRIBUTE
          CALL mp_exchange (ng, model, Istr, Iend, Jstr, Jend,          &
     &                      LBi, UBi, LBj, UBj, 1, 1,                   &
     &                      NghostPoints, .FALSE., .FALSE.,             &
     &                      GRID(ng) % latv(LBi,LBj)) 
# endif
        END IF
      END IF
# ifdef CURVGRID
!
!  Read in angle (radians) between XI-axis and EAST at RHO-points.
!
      IF (gotang) THEN
        status=nf_fread(ng, model, ncGRDid(ng), grdangid, 0, r2dvar,    &
     &                  LBi, UBi, LBj, UBj, 1, 1,                       &
     &                  Fscl, Fmin, Fmax,                               &
#  ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
#  endif
     &                  GRID(ng) % angler(LBi,LBj))
        IF (InpThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,30) 'angle', TRIM(ncname)
          exit_flag=2
          ioerror=status
          RETURN
        END IF
#  if (defined EW_PERIODIC || defined NS_PERIODIC) && !defined DISTRIBUTE
        CALL exchange_r2d_tile (ng, model, 1, Lm(ng), 1, Mm(ng),        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          GRID(ng) % angler)
#  endif
      END IF
# endif
!
! Close GRID NetCDF file.
!
      status=nf_close(ncGRDid(ng))
      ncGRDid(ng)=-1
!
  10  FORMAT (/,' GET_GRID - unable to find grid variable: ',a,         &
     &        /,12x,'in grid NetCDF file: ',a)
  20  FORMAT (/,' GET_GRID - unable to open grid NetCDF file: ',a)
  30  FORMAT (/,' GET_GRID - error while reading variable: ',a,         &
     &        /,12x,'in grid NetCDF file: ',a)
#else
      SUBROUTINE get_grid
#endif
      RETURN
      END SUBROUTINE get_grid
