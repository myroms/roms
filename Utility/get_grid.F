#include "cppdefs.h"
#ifndef ANA_GRID
      SUBROUTINE get_grid (ng)
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
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
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      logical :: gotang = .false.
      logical :: gotel = .false.
      logical :: gotf = .false.
      logical :: gotdm = .false.
      logical :: gotdn = .false.
      logical :: goth = .false.
      logical :: gotlatr = .false.
      logical :: gotlonr = .false.
# ifdef MASKING
      logical :: gotmskp = .false.
      logical :: gotmskr = .false.
      logical :: gotmsku = .false.
      logical :: gotmskv = .false.
# endif
      logical :: gotpm = .false.
      logical :: gotpn = .false.
      logical :: gotsph = .false.
      logical :: gotxl = .false.
      logical :: gotxr = .false.
      logical :: gotyr = .false.
# ifdef ICESHELF
      logical :: gotzice = .false.

      integer :: grdziceid
# endif
      integer :: LBi, UBi, LBj, UBj
# ifdef DISTRIBUTE
      integer, parameter :: Nghost = GHOST_POINTS
      integer :: Itile, Jtile, Istr, Iend, Jstr, Jend
# endif
      integer :: grdangid, grdelid, grdfid, grddmid, grddnid
      integer :: grdhid, grdlatrid, grdlonrid, grdpmid, grdpnid
# ifdef MASKING
      integer :: mskpid, mskrid, mskuid, mskvid
# endif
      integer :: grdsphid, grdxlid, grdxrid, grdyrid
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
      IF (exit_flag.ne.0) RETURN
      ncname=grdname(ng)
      ifield=0
      CALL opencdf (ng, 1, ncname, fname, N(ng), ifield, nrec)
      IF (exit_flag.ne.0) RETURN
!
!  Scan variable list from input NetCDF and check for grid variables.
!
      DO i=1,nvars
        IF (TRIM(varnam(i)).eq.'xl') THEN
          grdxlid=i
          gotxl=.true.
        ELSE IF (TRIM(varnam(i)).eq.'el') THEN
          grdelid=i
          gotel=.true.
        ELSE IF (TRIM(varnam(i)).eq.'spherical') THEN
          grdsphid=i
          gotsph=.true.
# ifdef MASKING
        ELSE IF (TRIM(varnam(i)).eq.'mask_rho') THEN
          mskrid=i
          gotmskr=.true.
        ELSE IF (TRIM(varnam(i)).eq.'mask_u') THEN
          mskuid=i
          gotmsku=.true.
        ELSE IF (TRIM(varnam(i)).eq.'mask_v') THEN
          mskvid=i
          gotmskv=.true.
        ELSE IF (TRIM(varnam(i)).eq.'mask_psi') THEN
          mskpid=i
          gotmskp=.true.
# endif
        ELSE IF (TRIM(varnam(i)).eq.'h') THEN
          grdhid=i
          goth=.true.
# ifdef ICESHELF
        ELSE IF (TRIM(varnam(i)).eq.'zice') THEN
          grdziceid=i
          gotzice=.true.
# endif
        ELSE IF (TRIM(varnam(i)).eq.'f') THEN
          grdfid=i
          gotf=.true.
        ELSE IF (TRIM(varnam(i)).eq.'pm') THEN
          grdpmid=i
          gotpm=.true.
        ELSE IF (TRIM(varnam(i)).eq.'pn') THEN
          grdpnid=i
          gotpn=.true.
        ELSE IF (TRIM(varnam(i)).eq.'dndx') THEN
          grddnid=i
          gotdn=.true.
        ELSE IF (TRIM(varnam(i)).eq.'dmde') THEN
          grddmid=i
          gotdm=.true.
        ELSE IF (TRIM(varnam(i)).eq.'x_rho') THEN
          grdxrid=i
          gotxr=.true.
        ELSE IF (TRIM(varnam(i)).eq.'y_rho') THEN
          grdyrid=i
          gotyr=.true.
        ELSE IF (TRIM(varnam(i)).eq.'lon_rho') THEN
          grdlonrid=i
          gotlonr=.true.
        ELSE IF (TRIM(varnam(i)).eq.'lat_rho') THEN
          grdlatrid=i
          gotlatr=.true.
        ELSE IF (TRIM(varnam(i)).eq.'angle') THEN
          grdangid=i
          gotang=.true.
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
      IF (ncgrdid(ng).eq.-1) THEN
        status=nf_open(TRIM(ncname), nf_nowrite, ncgrdid(ng))
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,20) TRIM(ncname)
          exit_flag=2
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
      status=nf_get_var1_TYPE(ncgrdid(ng), grdxlid, 1, xl(ng))
      IF (status.ne.nf_noerr) THEN
        WRITE (stdout,30) 'xl', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      status=nf_get_var1_TYPE(ncgrdid(ng), grdelid, 1, el(ng))
      IF (status.ne.nf_noerr) THEN
        WRITE (stdout,30) 'el', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
!
!  Read in logical switch for spherical grid configuration.
!
      status=nf_get_var1_text(ncgrdid(ng), grdsphid, 1, char1)
      IF (status.ne.nf_noerr) THEN
        WRITE (stdout,30) 'spherical', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF ((char1.eq.'t').or.(char1.eq.'T')) THEN
        spherical=.true.
      ELSE
        spherical=.false.
      END IF
!
!-----------------------------------------------------------------------
!  Read in grid arrays.
!-----------------------------------------------------------------------
!
      LBi=lbound(GRID(ng)%h,DIM=1)
      UBi=ubound(GRID(ng)%h,DIM=1)
      LBj=lbound(GRID(ng)%h,DIM=2)
      UBj=ubound(GRID(ng)%h,DIM=2)
# ifdef DISTRIBUTE
      CALL get_tile (ng, MyRank, Itile, Jtile, Istr, Iend, Jstr, Jend)
# endif

# ifdef MASKING
!
!  Read in mask at RHO-points.
!
      status=nf_fread(ng, ncgrdid(ng), mskrid, 0, r2dvar,               &
     &                LBi, UBi, LBj, UBj, 1, 1,                         &
     &                Fscl, Fmin, Fmax,                                 &
     &                GRID(ng) % rmask(LBi,LBj),                        &
     &                GRID(ng) % rmask(LBi,LBj))
      IF (InpThread.and.(status.ne.nf_noerr)) THEN
        WRITE (stdout,30) 'mask_rho', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
#  if (defined EW_PERIODIC || defined NS_PERIODIC) && !defined DISTRIBUTE
      CALL exchange_r2d_tile (ng, 1, Lm(ng), 1, Mm(ng),                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        GRID(ng) % rmask)
#  endif
!
!  Read in mask at U-points.
!
      status=nf_fread(ng, ncgrdid(ng), mskuid, 0, u2dvar,               &
     &                LBi, UBi, LBj, UBj, 1, 1,                         &
     &                Fscl, Fmin, Fmax,                                 &
     &                GRID(ng) % umask(LBi,LBj),                        &
     &                GRID(ng) % umask(LBi,LBj))
      IF (InpThread.and.(status.ne.nf_noerr)) THEN
        WRITE (stdout,30) 'mask_u', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
#  if (defined EW_PERIODIC || defined NS_PERIODIC) && !defined DISTRIBUTE
      CALL exchange_u2d_tile (ng, 1, Lm(ng), 1, Mm(ng),                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        GRID(ng) % umask)
#  endif
!
!  Read in mask at V-points.
!
      status=nf_fread(ng, ncgrdid(ng), mskvid, 0, v2dvar,               &
     &                LBi, UBi, LBj, UBj, 1, 1,                         &
     &                Fscl, Fmin, Fmax,                                 &
     &                GRID(ng) % vmask(LBi,LBj),                        &
     &                GRID(ng) % vmask(LBi,LBj))
      IF (InpThread.and.(status.ne.nf_noerr)) THEN
        WRITE (stdout,30) 'mask_v', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
#  if (defined EW_PERIODIC || defined NS_PERIODIC) && !defined DISTRIBUTE
      CALL exchange_v2d_tile (ng, 1, Lm(ng), 1, Mm(ng),                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        GRID(ng) % vmask)
#  endif
!
!  Read in mask at PSI-points.
!
      status=nf_fread(ng, ncgrdid(ng), mskpid, 0, p2dvar,               &
     &                LBi, UBi, LBj, UBj, 1, 1,                         &
     &                Fscl, Fmin, Fmax,                                 &
     &                GRID(ng) % pmask(LBi,LBj),                        &
     &                GRID(ng) % pmask(LBi,LBj))
      IF (InpThread.and.(status.ne.nf_noerr)) THEN
        WRITE (stdout,30) 'mask_psi', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
#  if (defined EW_PERIODIC || defined NS_PERIODIC) && !defined DISTRIBUTE
      CALL exchange_p2d_tile (ng, 1, Lm(ng), 1, Mm(ng),                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        GRID(ng) % pmask)
#  endif
# endif
!
!  Read in bathymetry.
!
      status=nf_fread(ng, ncgrdid(ng), grdhid, 0, r2dvar,               &
     &                LBi, UBi, LBj, UBj, 1, 1,                         &
     &                Fscl, hmin(ng), hmax(ng),                         &
# ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
# endif
     &                GRID(ng) % h(LBi,LBj))
      IF (InpThread.and.(status.ne.nf_noerr)) THEN
        WRITE (stdout,30) 'h', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
# if (defined EW_PERIODIC || defined NS_PERIODIC) && !defined DISTRIBUTE
      CALL exchange_r2d_tile (ng, 1, Lm(ng), 1, Mm(ng),                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        GRID(ng) % h)
# endif
# ifdef ICESHELF
!
!  Read in ice shelf thicknesses.
!
      status=nf_fread(ng, ncgrdid(ng), grdziceid, 0, r2dvar,            &
     &                LBi, UBi, LBj, UBj, 1, 1,                         &
     &                Fscl, Fmin, Fmax,                                 &
#  ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
#  endif
     &                GRID(ng) % zice(LBi,LBj))
      IF (InpThread.and.(status.ne.nf_noerr)) THEN
        WRITE (stdout,30) 'zice', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
#  if (defined EW_PERIODIC || defined NS_PERIODIC) && !defined DISTRIBUTE
      CALL exchange_r2d_tile (ng, 1, Lm(ng), 1, Mm(ng),                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        GRID(ng) % zice)
#  endif
# endif
!
!  Read in Coriolis parameter.
!
      status=nf_fread(ng, ncgrdid(ng), grdfid, 0, r2dvar,               &
     &                LBi, UBi, LBj, UBj, 1, 1,                         &
     &                Fscl, Fmin, Fmax,                                 &
# ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
# endif
     &                GRID(ng) % f(LBi,LBj))
      IF (InpThread.and.(status.ne.nf_noerr)) THEN
        WRITE (stdout,30) 'f', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
# if (defined EW_PERIODIC || defined NS_PERIODIC) && !defined DISTRIBUTE
      CALL exchange_r2d_tile (ng, 1, Lm(ng), 1, Mm(ng),                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        GRID(ng) % f)
# endif
!
!  Read in coordinate transfomation metrics (m,n) associated with the
!  differential distances in XI and ETA.
!
      status=nf_fread(ng, ncgrdid(ng), grdpmid, 0, r2dvar,              &
     &                LBi, UBi, LBj, UBj, 1, 1,                         &
     &                Fscl, Fmin, Fmax,                                 &
# ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
# endif
     &                GRID(ng) % pm(LBi,LBj))
      IF (InpThread.and.(status.ne.nf_noerr)) THEN
        WRITE (stdout,30) 'pm', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      status=nf_fread(ng, ncgrdid(ng), grdpnid, 0, r2dvar,              &
     &                LBi, UBi, LBj, UBj, 1, 1,                         &
     &                Fscl, Fmin, Fmax,                                 &
# ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
# endif
     &                GRID(ng) % pn(LBi,LBj))
      IF (InpThread.and.(status.ne.nf_noerr)) THEN
        WRITE (stdout,30) 'pn', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
# if (defined EW_PERIODIC || defined NS_PERIODIC) && !defined DISTRIBUTE
      CALL exchange_r2d_tile (ng, 1, Lm(ng), 1, Mm(ng),                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        GRID(ng) % pm)
      CALL exchange_r2d_tile (ng, 1, Lm(ng), 1, Mm(ng),                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        GRID(ng) % pn)
# endif
# if (defined CURVGRID && defined UV_ADV)
!
!  Read in derivatives of inverse metrics factors: rotational factors.
!
      status=nf_fread(ng, ncgrdid(ng), grddmid, 0, r2dvar,              &
     &                LBi, UBi, LBj, UBj, 1, 1,                         &
     &                Fscl, Fmin, Fmax,                                 &
# ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
# endif
     &                GRID(ng) % dmde(LBi,LBj))
      IF (InpThread.and.(status.ne.nf_noerr)) THEN
        WRITE (stdout,30) 'dmde', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      status=nf_fread(ng, ncgrdid(ng), grddnid, 0, r2dvar,              &
     &                LBi, UBi, LBj, UBj, 1, 1,                         &
     &                Fscl, Fmin, Fmax,                                 &
# ifdef MASKING
     &                GRID(ng) % rmask(LBi,LBj),                        &
# endif
     &                GRID(ng) % dndx(LBi,LBj))
      IF (InpThread.and.(status.ne.nf_noerr)) THEN
        WRITE (stdout,30) 'dndx', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
#  if (defined EW_PERIODIC || defined NS_PERIODIC) && !defined DISTRIBUTE
      CALL exchange_r2d_tile (ng, 1, Lm(ng), 1, Mm(ng),                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        GRID(ng) % dmde)
      CALL exchange_r2d_tile (ng, 1, Lm(ng), 1, Mm(ng),                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        GRID(ng) % dndx)
#  endif
# endif
!
!  Read in (x,y) coordinates at RHO-points.
!
      IF (gotxr) THEN
        status=nf_fread(ng, ncgrdid(ng), grdxrid, 0, r2dvar,            &
     &                  LBi, UBi, LBj, UBj, 1, 1,                       &
     &                  Fscl, Fmin, Fmax,                               &
# ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
# endif
     &                  GRID(ng) % xr(LBi,LBj))
        IF (InpThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,30) 'x_rho', TRIM(ncname)
          exit_flag=2
          RETURN
        END IF
# ifdef DISTRIBUTE
        CALL mp_exchange (ng, Istr, Iend, Jstr, Jend,                   &
     &                    LBi, UBi, LBj, UBj, 1, 1,                     &
     &                    Nghost, .false., .false.,                     &
     &                    GRID(ng) % xr(LBi,LBj)) 
# endif
      END IF
      IF (gotyr) THEN
        status=nf_fread(ng, ncgrdid(ng), grdyrid, 0, r2dvar,            &
     &                  LBi, UBi, LBj, UBj, 1, 1,                       &
     &                  Fscl, Fmin, Fmax,                               &
# ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
# endif
     &                  GRID(ng) % yr(LBi,LBj))
        IF (InpThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,30) 'y_rho', TRIM(ncname)
          exit_flag=2
          RETURN
        END IF
# ifdef DISTRIBUTE
        CALL mp_exchange (ng, Istr, Iend, Jstr, Jend,                   &
     &                    LBi, UBi, LBj, UBj, 1, 1,                     &
     &                    Nghost, .false., .false.,                     &
     &                    GRID(ng) % yr(LBi,LBj)) 
# endif
      END IF
!
!  Read in (lon,lat) coordinates at RHO-points.
!
      IF (spherical) THEN
        IF (gotlonr) THEN
          status=nf_fread(ng, ncgrdid(ng), grdlonrid, 0, r2dvar,        &
     &                    LBi, UBi, LBj, UBj, 1, 1,                     &
     &                    Fscl, Fmin, Fmax,                             &
# ifdef MASKING
     &                    GRID(ng) % rmask(LBi,LBj),                    &
# endif
     &                    GRID(ng) % lonr(LBi,LBj))
          IF (InpThread.and.(status.ne.nf_noerr)) THEN
            WRITE (stdout,30) 'lon_rho', TRIM(ncname)
            exit_flag=2
            RETURN
          END IF
# ifdef DISTRIBUTE
          CALL mp_exchange (ng, Istr, Iend, Jstr, Jend,                 &
     &                      LBi, UBi, LBj, UBj, 1, 1,                   &
     &                      Nghost, .false., .false.,                   &
     &                      GRID(ng) % lonr(LBi,LBj)) 
# endif
        END IF
        IF (gotlatr) THEN
          status=nf_fread(ng, ncgrdid(ng), grdlatrid, 0, r2dvar,        &
     &                    LBi, UBi, LBj, UBj, 1, 1,                     &
     &                    Fscl, Fmin, Fmax,                             &
# ifdef MASKING
     &                    GRID(ng) % rmask(LBi,LBj),                    &
# endif
     &                    GRID(ng) % latr(LBi,LBj))
          IF (InpThread.and.(status.ne.nf_noerr)) THEN
            WRITE (stdout,30) 'lat_rho', TRIM(ncname)
            exit_flag=2
            RETURN
          END IF
# ifdef DISTRIBUTE
          CALL mp_exchange (ng, Istr, Iend, Jstr, Jend,                 &
     &                      LBi, UBi, LBj, UBj, 1, 1,                   &
     &                      Nghost, .false., .false.,                   &
     &                      GRID(ng) % latr(LBi,LBj)) 
# endif
        END IF
      END IF
# ifdef CURVGRID
!
!  Read in angle (radians) between XI-axis and EAST at RHO-points.
!
      IF (gotang) THEN
        status=nf_fread(ng, ncgrdid(ng), grdangid, 0, r2dvar,           &
     &                  LBi, UBi, LBj, UBj, 1, 1,                       &
     &                  Fscl, Fmin, Fmax,                               &
#  ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
#  endif
     &                  GRID(ng) % angler(LBi,LBj))
        IF (InpThread.and.(status.ne.nf_noerr)) THEN
          WRITE (stdout,30) 'angle', TRIM(ncname)
          exit_flag=2
          RETURN
        END IF
#  if (defined EW_PERIODIC || defined NS_PERIODIC) && !defined DISTRIBUTE
        CALL exchange_r2d_tile (ng, 1, Lm(ng), 1, Mm(ng),               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          GRID(ng) % angler)
#  endif
      END IF
# endif
!
! Close GRID NetCDF file.
!
      status=nf_close(ncgrdid(ng))
      ncgrdid(ng)=-1
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
