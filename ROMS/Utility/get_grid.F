#include "cppdefs.h"
#ifndef ANA_GRID
      SUBROUTINE get_grid (ng, model)
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2007 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
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
# if (defined EW_PERIODIC || defined NS_PERIODIC)
      USE exchange_2d_mod
# endif
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
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
# if defined AD_SENSITIVITY || defined OPT_OBSERVATIONS || \
     defined SO_SEMI
      logical :: gotRscope = .FALSE.
      logical :: gotUscope = .FALSE.
      logical :: gotVscope = .FALSE.
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
#  ifdef EW_PERIODIC
      logical :: EWperiodic=.TRUE.
#  else
      logical :: EWperiodic=.FALSE.
#  endif
#  ifdef NS_PERIODIC
      logical :: NSperiodic=.TRUE.
#  else
      logical :: NSperiodic=.FALSE.
#  endif
      integer :: Itile, Jtile
# endif
      integer :: Istr, Iend, Jstr, Jend
      integer :: grdangid, grdelid, grdfid, grddmid, grddnid
      integer :: grdhid, grdpmid, grdpnid
      integer :: grdlatpid, grdlatrid, grdlatuid, grdlatvid
      integer :: grdlonpid, grdlonrid, grdlonuid, grdlonvid
# if defined AD_SENSITIVITY || defined OPT_OBSERVATIONS || \
     defined SO_SEMI
      integer :: Rscopeid, Uscopeid, Vscopeid
# endif
# ifdef MASKING
      integer :: mskpid, mskrid, mskuid, mskvid
# endif
      integer :: grdsphid, grdxlid
      integer :: grdxpid, grdxrid, grdxuid, grdxvid
      integer :: grdypid, grdyrid, grdyuid, grdyvid
      integer :: gtype, i, nrec, nvd, status

      integer :: Vsize(4)

      integer :: nf_fread2d

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
      CALL opencdf (ng, 1, ncname, fname, N(ng), 0, nrec, nvd, Vsize)
      IF (exit_flag.ne.NoError) RETURN
!
!  Set Vsize to zero to deativate interpolation of input data to model
!  grid in "nf_fread2d".
!
      DO i=1,4
        Vsize(i)=0
      END DO
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
# if defined AD_SENSITIVITY || defined OPT_OBSERVATIONS || \
     defined SO_SEMI
        ELSE IF (TRIM(varnam(i)).eq.'scope_rho') THEN
          Rscopeid=i
          gotRscope=.TRUE.
        ELSE IF (TRIM(varnam(i)).eq.'scope_u') THEN
          Uscopeid=i
          gotUscope=.TRUE.
        ELSE IF (TRIM(varnam(i)).eq.'scope_v') THEN
          Vscopeid=i
          gotVscope=.TRUE.
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
# if defined AD_SENSITIVITY || defined OPT_OBSERVATIONS || \
     defined SO_SEMI
      IF (.not.gotRscope) THEN
        WRITE (stdout,10) 'scope_rho', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.gotUscope) THEN
        WRITE (stdout,10) 'scope_u', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.gotVscope) THEN
        WRITE (stdout,10) 'scope_v', TRIM(ncname)
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
        status=nf90_open(TRIM(ncname), nf90_nowrite, ncGRDid(ng))
        IF (status.ne.nf90_noerr) THEN
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
      status=nf90_get_var(ncGRDid(ng), grdxlid, xl(ng))
      IF (status.ne.nf90_noerr) THEN
        WRITE (stdout,30) 'xl', TRIM(ncname)
        exit_flag=2
        ioerror=status
        RETURN
      END IF
      status=nf90_get_var(ncGRDid(ng), grdelid, el(ng))
      IF (status.ne.nf90_noerr) THEN
        WRITE (stdout,30) 'el', TRIM(ncname)
        exit_flag=2
        ioerror=status
        RETURN
      END IF
!
!  Read in logical switch for spherical grid configuration.
!
      status=nf90_get_var(ncGRDid(ng), grdsphid, char1)
      IF (status.ne.nf90_noerr) THEN
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
# else
      Istr=1
      Iend=Lm(ng)
      Jstr=1
      Jend=Mm(ng)
# endif

# ifdef MASKING
!
!  Read in Land/Sea masking arrays.
!
      gtype=r2dvar
      status=nf_fread2d(ng, model, ncGRDid(ng), mskrid, 0, gtype,       &
     &                  Vsize, LBi, UBi, LBj, UBj,                      &
     &                  Fscl, Fmin, Fmax,                               &
     &                  GRID(ng) % rmask(LBi,LBj),                      &
     &                  GRID(ng) % rmask(LBi,LBj))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,30) 'mask_rho', TRIM(ncname)
        END IF
        exit_flag=2
        ioerror=status
        RETURN
      END IF

      gtype=u2dvar
      status=nf_fread2d(ng, model, ncGRDid(ng), mskuid, 0, gtype,       &
     &                  Vsize, LBi, UBi, LBj, UBj,                      &
     &                  Fscl, Fmin, Fmax,                               &
     &                  GRID(ng) % umask(LBi,LBj),                      &
     &                  GRID(ng) % umask(LBi,LBj))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,30) 'mask_u', TRIM(ncname)
        END IF
        exit_flag=2
        ioerror=status
        RETURN
      END IF

      gtype=v2dvar
      status=nf_fread2d(ng, model, ncGRDid(ng), mskvid, 0, gtype,       &
     &                  Vsize, LBi, UBi, LBj, UBj,                      &
     &                  Fscl, Fmin, Fmax,                               &
     &                  GRID(ng) % vmask(LBi,LBj),                      &
     &                  GRID(ng) % vmask(LBi,LBj))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,30) 'mask_v', TRIM(ncname)
        END IF
        exit_flag=2
        ioerror=status
        RETURN
      END IF

      gtype=p2dvar
      status=nf_fread2d(ng, model, ncGRDid(ng), mskpid, 0, gtype,       &
     &                  Vsize, LBi, UBi, LBj, UBj,                      &
     &                  Fscl, Fmin, Fmax,                               &
     &                  GRID(ng) % pmask(LBi,LBj),                      &
     &                  GRID(ng) % pmask(LBi,LBj))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,30) 'mask_psi', TRIM(ncname)
        END IF
        exit_flag=2
        ioerror=status
        RETURN
      END IF
#  if defined EW_PERIODIC || defined NS_PERIODIC
!
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        GRID(ng) % rmask)
      CALL exchange_u2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        GRID(ng) % umask)
      CALL exchange_v2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        GRID(ng) % vmask)
      CALL exchange_p2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        GRID(ng) % pmask)
#  endif
#  ifdef DISTRIBUTE
!
      CALL mp_exchange2d (ng, model, 4, Istr, Iend, Jstr, Jend,         &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    GRID(ng) % pmask,                             &
     &                    GRID(ng) % rmask,                             &
     &                    GRID(ng) % umask,                             &
     &                    GRID(ng) % vmask)
#  endif
# endif

# if defined AD_SENSITIVITY || defined OPT_OBSERVATIONS || \
     defined SO_SEMI
!
!  Read in adjoint sensitivity spatial scope masking arrays.
!
      gtype=r2dvar
      status=nf_fread2d(ng, model, ncGRDid(ng), Rscopeid, 0, gtype,     &
     &                  Vsize, LBi, UBi, LBj, UBj,                      &
     &                  Fscl, Fmin, Fmax,                               &
     &                  GRID(ng) % Rscope(LBi,LBj),                     &
     &                  GRID(ng) % Rscope(LBi,LBj))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,30) 'scope_rho', TRIM(ncname)
        END IF
        exit_flag=2
        ioerror=status
        RETURN
      END IF

      gtype=u2dvar
      status=nf_fread2d(ng, model, ncGRDid(ng), Uscopeid, 0, gtype,     &
     &                  Vsize, LBi, UBi, LBj, UBj,                      &
     &                  Fscl, Fmin, Fmax,                               &
     &                  GRID(ng) % Uscope(LBi,LBj),                     &
     &                  GRID(ng) % Uscope(LBi,LBj))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,30) 'scope_u', TRIM(ncname)
        END IF
        exit_flag=2
        ioerror=status
        RETURN
      END IF

      gtype=v2dvar
      status=nf_fread2d(ng, model, ncGRDid(ng), Vscopeid, 0, gtype,     &
     &                  Vsize, LBi, UBi, LBj, UBj,                      &
     &                  Fscl, Fmin, Fmax,                               &
     &                  GRID(ng) % Vscope(LBi,LBj),                     &
     &                  GRID(ng) % Vscope(LBi,LBj))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,30) 'scope_v', TRIM(ncname)
        END IF
        exit_flag=2
        ioerror=status
        RETURN
      END IF
#  if defined EW_PERIODIC || defined NS_PERIODIC
!
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        GRID(ng) % Rscope)
      CALL exchange_u2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        GRID(ng) % Uscope)
      CALL exchange_v2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        GRID(ng) % Vscope)
#  endif
#  ifdef DISTRIBUTE
!
      CALL mp_exchange2d (ng, model, 3, Istr, Iend, Jstr, Jend,         &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    GRID(ng) % Rscope,                            &
     &                    GRID(ng) % Uscope,                            &
     &                    GRID(ng) % Vscope)
#  endif
# endif
!
!  Read in bathymetry.
!
      gtype=r2dvar
      status=nf_fread2d(ng, model, ncGRDid(ng), grdhid, 0, gtype,       &
     &                  Vsize, LBi, UBi, LBj, UBj,                      &
     &                  Fscl, hmin(ng), hmax(ng),                       &
# ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
# endif
     &                  GRID(ng) % h(LBi,LBj))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,30) 'h', TRIM(ncname)
        END IF
        exit_flag=2
        ioerror=status
        RETURN
      END IF
# if defined EW_PERIODIC || defined NS_PERIODIC
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        GRID(ng) % h)
# endif
# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, model, 1, Istr, Iend, Jstr, Jend,         &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    GRID(ng) % h)
# endif
# ifdef ICESHELF
!
!  Read in ice shelf thicknesses.
!
      gtype=r2dvar
      status=nf_fread2d(ng, model, ncGRDid(ng), grdziceid, 0, gtype,    &
     &                  Vsize, LBi, UBi, LBj, UBj,                      &
     &                  Fscl, Fmin, Fmax,                               &
#  ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
#  endif
     &                  GRID(ng) % zice(LBi,LBj))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,30) 'zice', TRIM(ncname)
        END IF
        exit_flag=2
        ioerror=status
        RETURN
      END IF
#  if defined EW_PERIODIC || defined NS_PERIODIC
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        GRID(ng) % zice)
#  endif
#  ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, model, 1, Istr, Iend, Jstr, Jend,         &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    GRID(ng) % zice)
#  endif
# endif
!
!  Read in Coriolis parameter.
!
      gtype=r2dvar
      status=nf_fread2d(ng, model, ncGRDid(ng), grdfid, 0, gtype,       &
     &                  Vsize, LBi, UBi, LBj, UBj,                      &
     &                  Fscl, Fmin, Fmax,                               &
# ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
# endif
     &                  GRID(ng) % f(LBi,LBj))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,30) 'f', TRIM(ncname)
        END IF
        exit_flag=2
        ioerror=status
        RETURN
      END IF
# if defined EW_PERIODIC || defined NS_PERIODIC
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        GRID(ng) % f)
# endif
# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, model, 1, Istr, Iend, Jstr, Jend,         &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    GRID(ng) % f)
# endif

!
!  Read in coordinate transfomation metrics (m,n) associated with the
!  differential distances in XI and ETA.
!
      gtype=r2dvar
      status=nf_fread2d(ng, model, ncGRDid(ng), grdpmid, 0, gtype,      &
     &                  Vsize, LBi, UBi, LBj, UBj,                      &
     &                  Fscl, Fmin, Fmax,                               &
# ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
# endif
     &                  GRID(ng) % pm(LBi,LBj))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,30) 'pm', TRIM(ncname)
        END IF
        exit_flag=2
        ioerror=status
        RETURN
      END IF

      gtype=r2dvar
      status=nf_fread2d(ng, model, ncGRDid(ng), grdpnid, 0, gtype,      &
     &                  Vsize, LBi, UBi, LBj, UBj,                      &
     &                  Fscl, Fmin, Fmax,                               &
# ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
# endif
     &                  GRID(ng) % pn(LBi,LBj))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,30) 'pn', TRIM(ncname)
        END IF
        exit_flag=2
        ioerror=status
        RETURN
      END IF
# if defined EW_PERIODIC || defined NS_PERIODIC
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        GRID(ng) % pm)
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        GRID(ng) % pn)
# endif
# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, model, 2, Istr, Iend, Jstr, Jend,         &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn)
# endif

# if (defined CURVGRID && defined UV_ADV)
!
!  Read in derivatives of inverse metrics factors: rotational factors.
!
      gtype=r2dvar
      status=nf_fread2d(ng, model, ncGRDid(ng), grddmid, 0, gtype,      &
     &                  Vsize, LBi, UBi, LBj, UBj,                      &
     &                  Fscl, Fmin, Fmax,                               &
# ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
# endif
     &                  GRID(ng) % dmde(LBi,LBj))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,30) 'dmde', TRIM(ncname)
        END IF
        exit_flag=2
        ioerror=status
        RETURN
      END IF

      gtype=r2dvar
      status=nf_fread2d(ng, model, ncGRDid(ng), grddnid, 0, gtype,      &
     &                  Vsize, LBi, UBi, LBj, UBj,                      &
     &                  Fscl, Fmin, Fmax,                               &
# ifdef MASKING
     &                  GRID(ng) % rmask(LBi,LBj),                      &
# endif
     &                  GRID(ng) % dndx(LBi,LBj))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,30) 'dndx', TRIM(ncname)
        END IF
        exit_flag=2
        ioerror=status
        RETURN
      END IF
#  if defined EW_PERIODIC || defined NS_PERIODIC
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        GRID(ng) % dmde)
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        GRID(ng) % dndx)
#  endif
#  ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, model, 2, Istr, Iend, Jstr, Jend,         &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    GRID(ng) % dmde,                              &
     &                    GRID(ng) % dndx)
#  endif
# endif
!
!  Read in (x,y) coordinates at PSI-points.
!
      IF (gotxp) THEN
        gtype=p2dvar
        status=nf_fread2d(ng, model, ncGRDid(ng), grdxpid, 0, gtype,    &
     &                    Vsize, LBi, UBi, LBj, UBj,                    &
     &                    Fscl, Fmin, Fmax,                             &
# ifdef MASKING
     &                    GRID(ng) % pmask(LBi,LBj),                    &
# endif
     &                    GRID(ng) % xp(LBi,LBj))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,30) 'x_p', TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        END IF
# ifdef DISTRIBUTE
        CALL mp_exchange2d (ng, model, 1, Istr, Iend, Jstr, Jend,       &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints, .FALSE., .FALSE.,             &
     &                      GRID(ng) % xp)
# endif
      END IF
      IF (gotyp) THEN
        gtype=p2dvar
        status=nf_fread2d(ng, model,  ncGRDid(ng), grdypid, 0, gtype,   &
     &                    Vsize, LBi, UBi, LBj, UBj,                    &
     &                    Fscl, Fmin, Fmax,                             &
# ifdef MASKING
     &                    GRID(ng) % pmask(LBi,LBj),                    &
# endif
     &                    GRID(ng) % yp(LBi,LBj))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,30) 'y_p', TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        END IF
# ifdef DISTRIBUTE
        CALL mp_exchange2d (ng, model, 1, Istr, Iend, Jstr, Jend,       &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints, .FALSE., .FALSE.,             &
     &                      GRID(ng) % yp)
# endif
      END IF
!
!  Read in (x,y) coordinates at RHO-points.
!
      IF (gotxr) THEN
        gtype=r2dvar
        status=nf_fread2d(ng, model, ncGRDid(ng), grdxrid, 0, gtype,    &
     &                    Vsize, LBi, UBi, LBj, UBj,                    &
     &                    Fscl, Fmin, Fmax,                             &
# ifdef MASKING
     &                    GRID(ng) % rmask(LBi,LBj),                    &
# endif
     &                    GRID(ng) % xr(LBi,LBj))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,30) 'x_rho', TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        END IF
# ifdef DISTRIBUTE
        CALL mp_exchange2d (ng, model, 1, Istr, Iend, Jstr, Jend,       &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints, .FALSE., .FALSE.,             &
     &                      GRID(ng) % xr)
# endif
      END IF
      IF (gotyr) THEN
        gtype=r2dvar
        status=nf_fread2d(ng, model, ncGRDid(ng), grdyrid, 0, gtype,    &
     &                    Vsize, LBi, UBi, LBj, UBj,                    &
     &                    Fscl, Fmin, Fmax,                             &
# ifdef MASKING
     &                    GRID(ng) % rmask(LBi,LBj),                    &
# endif
     &                    GRID(ng) % yr(LBi,LBj))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,30) 'y_rho', TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        END IF
# ifdef DISTRIBUTE
        CALL mp_exchange2d (ng, model, 1, Istr, Iend, Jstr, Jend,       &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints, .FALSE., .FALSE.,             &
     &                      GRID(ng) % yr)
# endif
      END IF
!
!  Read in (x,y) coordinates at U-points.
!
      IF (gotxu) THEN
        gtype=u2dvar
        status=nf_fread2d(ng, model, ncGRDid(ng), grdxuid, 0, gtype,    &
     &                    Vsize, LBi, UBi, LBj, UBj,                    &
     &                    Fscl, Fmin, Fmax,                             &
# ifdef MASKING
     &                    GRID(ng) % umask(LBi,LBj),                    &
# endif
     &                    GRID(ng) % xu(LBi,LBj))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,30) 'x_u', TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        END IF
# ifdef DISTRIBUTE
        CALL mp_exchange2d (ng, model, 1, Istr, Iend, Jstr, Jend,       &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints, .FALSE., .FALSE.,             &
     &                      GRID(ng) % xu)
# endif
      END IF
      IF (gotyu) THEN
        gtype=u2dvar
        status=nf_fread2d(ng, model, ncGRDid(ng), grdyuid, 0, gtype,    &
     &                    Vsize, LBi, UBi, LBj, UBj,                    &
     &                    Fscl, Fmin, Fmax,                             &
# ifdef MASKING
     &                    GRID(ng) % umask(LBi,LBj),                    &
# endif
     &                    GRID(ng) % yu(LBi,LBj))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,30) 'y_u', TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        END IF
# ifdef DISTRIBUTE
        CALL mp_exchange2d (ng, model, 1, Istr, Iend, Jstr, Jend,       &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints, .FALSE., .FALSE.,             &
     &                      GRID(ng) % yu)
# endif
      END IF
!
!  Read in (x,y) coordinates at V-points.
!
      IF (gotxv) THEN
        gtype=v2dvar
        status=nf_fread2d(ng, model, ncGRDid(ng), grdxvid, 0, gtype,    &
     &                    Vsize, LBi, UBi, LBj, UBj,                    &
     &                    Fscl, Fmin, Fmax,                             &
# ifdef MASKING
     &                    GRID(ng) % vmask(LBi,LBj),                    &
# endif
     &                    GRID(ng) % xv(LBi,LBj))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,30) 'x_v', TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        END IF
# ifdef DISTRIBUTE
        CALL mp_exchange2d (ng, model, 1, Istr, Iend, Jstr, Jend,       &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints, .FALSE., .FALSE.,             &
     &                      GRID(ng) % xv)
# endif
      END IF
      IF (gotyv) THEN
        gtype=v2dvar
        status=nf_fread2d(ng, model, ncGRDid(ng), grdyvid, 0, gtype,    &
     &                    Vsize, LBi, UBi, LBj, UBj,                    &
     &                    Fscl, Fmin, Fmax,                             &
# ifdef MASKING
     &                    GRID(ng) % vmask(LBi,LBj),                    &
# endif
     &                    GRID(ng) % yv(LBi,LBj))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,30) 'y_v', TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        END IF
# ifdef DISTRIBUTE
        CALL mp_exchange2d (ng, model, 1, Istr, Iend, Jstr, Jend,       &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints, .FALSE., .FALSE.,             &
     &                      GRID(ng) % yv)
# endif
      END IF
!
!  Read in (lon,lat) coordinates at PSI-points.
!
      IF (spherical) THEN
        IF (gotlonp) THEN
          gtype=p2dvar
          status=nf_fread2d(ng, model, ncGRDid(ng), grdlonpid, 0,       &
     &                      gtype, Vsize, LBi, UBi, LBj, UBj,           &
     &                      Fscl, Fmin, Fmax,                           &
# ifdef MASKING
     &                      GRID(ng) % pmask(LBi,LBj),                  &
# endif
     &                      GRID(ng) % lonp(LBi,LBj))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) 'lon_psi', TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          END IF
# ifdef DISTRIBUTE
          CALL mp_exchange2d (ng, model, 1, Istr, Iend, Jstr, Jend,     &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints, .FALSE., .FALSE.,           &
     &                        GRID(ng) % lonp)
# endif
        END IF
        IF (gotlatp) THEN
          gtype=p2dvar
          status=nf_fread2d(ng, model, ncGRDid(ng), grdlatpid, 0,       &
     &                      gtype, Vsize, LBi, UBi, LBj, UBj,           &
     &                      Fscl, Fmin, Fmax,                           &
# ifdef MASKING
     &                      GRID(ng) % pmask(LBi,LBj),                  &
# endif
     &                      GRID(ng) % latp(LBi,LBj))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) 'lat_psi', TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          END IF
# ifdef DISTRIBUTE
          CALL mp_exchange2d (ng, model, 1, Istr, Iend, Jstr, Jend,     &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints, .FALSE., .FALSE.,           &
     &                        GRID(ng) % latp)
# endif
        END IF
      END IF
!
!  Read in (lon,lat) coordinates at RHO-points.
!
      IF (spherical) THEN
        IF (gotlonr) THEN
          gtype=r2dvar
          status=nf_fread2d(ng, model, ncGRDid(ng), grdlonrid, 0,       &
     &                      gtype, Vsize, LBi, UBi, LBj, UBj,           &
     &                      Fscl, LonMin(ng), LonMax(ng),               &
# ifdef MASKING
     &                      GRID(ng) % rmask(LBi,LBj),                  &
# endif
     &                      GRID(ng) % lonr(LBi,LBj))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) 'lon_rho', TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          END IF
# ifdef DISTRIBUTE
          CALL mp_exchange2d (ng, model, 1, Istr, Iend, Jstr, Jend,     &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints, .FALSE., .FALSE.,           &
     &                        GRID(ng) % lonr)
# endif
        END IF
        IF (gotlatr) THEN
          gtype=r2dvar
          status=nf_fread2d(ng, model, ncGRDid(ng), grdlatrid, 0,       &
     &                      gtype, Vsize, LBi, UBi, LBj, UBj,           &
     &                      Fscl, LatMin(ng), LatMax(ng),               &
# ifdef MASKING
     &                      GRID(ng) % rmask(LBi,LBj),                  &
# endif
     &                      GRID(ng) % latr(LBi,LBj))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) 'lat_rho', TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          END IF
# ifdef DISTRIBUTE
          CALL mp_exchange2d (ng, model, 1, Istr, Iend, Jstr, Jend,     &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints, .FALSE., .FALSE.,           &
     &                        GRID(ng) % latr)
# endif
        END IF
      END IF
!
!  Read in (lon,lat) coordinates at U-points.
!
      IF (spherical) THEN
        IF (gotlonu) THEN
          gtype=u2dvar
          status=nf_fread2d(ng, model, ncGRDid(ng), grdlonuid, 0,       &
     &                      gtype, Vsize, LBi, UBi, LBj, UBj,           &
     &                      Fscl, Fmin, Fmax,                           &
# ifdef MASKING
     &                      GRID(ng) % umask(LBi,LBj),                  &
# endif
     &                      GRID(ng) % lonu(LBi,LBj))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) 'lon_u', TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          END IF
# ifdef DISTRIBUTE
          CALL mp_exchange2d (ng, model, 1, Istr, Iend, Jstr, Jend,     &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints, .FALSE., .FALSE.,           &
     &                        GRID(ng) % lonu)
# endif
        END IF
        IF (gotlatu) THEN
          gtype=u2dvar
          status=nf_fread2d(ng, model, ncGRDid(ng), grdlatuid, 0,       &
     &                      gtype, Vsize, LBi, UBi, LBj, UBj,           &
     &                      Fscl, Fmin, Fmax,                           &
# ifdef MASKING
     &                      GRID(ng) % umask(LBi,LBj),                  &
# endif
     &                      GRID(ng) % latu(LBi,LBj))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) 'lat_u', TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          END IF
# ifdef DISTRIBUTE
          CALL mp_exchange2d (ng, model, 1, Istr, Iend, Jstr, Jend,     &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints, .FALSE., .FALSE.,           &
     &                        GRID(ng) % latu)
# endif
        END IF
      END IF
!
!  Read in (lon,lat) coordinates at V-points.
!
      IF (spherical) THEN
        IF (gotlonv) THEN
          gtype=v2dvar
          status=nf_fread2d(ng, model, ncGRDid(ng), grdlonvid, 0,       &
     &                      gtype, Vsize, LBi, UBi, LBj, UBj,           &
     &                      Fscl, Fmin, Fmax,                           &
# ifdef MASKING
     &                      GRID(ng) % vmask(LBi,LBj),                  &
# endif
     &                      GRID(ng) % lonv(LBi,LBj))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) 'lon_v', TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          END IF
# ifdef DISTRIBUTE
          CALL mp_exchange2d (ng, model, 1, Istr, Iend, Jstr, Jend,     &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints, .FALSE., .FALSE.,           &
     &                        GRID(ng) % lonv)
# endif
        END IF
        IF (gotlatv) THEN
          gtype=v2dvar
          status=nf_fread2d(ng, model, ncGRDid(ng), grdlatvid, 0,       &
     &                      gtype, Vsize, LBi, UBi, LBj, UBj,           &
     &                      Fscl, Fmin, Fmax,                           &
# ifdef MASKING
     &                      GRID(ng) % vmask(LBi,LBj),                  &
# endif
     &                      GRID(ng) % latv(LBi,LBj))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) 'lat_v', TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          END IF
# ifdef DISTRIBUTE
          CALL mp_exchange2d (ng, model, 1, Istr, Iend, Jstr, Jend,     &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints, .FALSE., .FALSE.,           &
     &                        GRID(ng) % latv)
# endif
        END IF
      END IF
!
!  Read in angle (radians) between XI-axis and EAST at RHO-points.
!
      IF (gotang) THEN
        gtype=r2dvar
        status=nf_fread2d(ng, model, ncGRDid(ng), grdangid, 0,          &
     &                    gtype, Vsize, LBi, UBi, LBj, UBj,             &
     &                    Fscl, Fmin, Fmax,                             &
# ifdef MASKING
     &                    GRID(ng) % rmask(LBi,LBj),                    &
# endif
     &                    GRID(ng) % angler(LBi,LBj))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,30) 'angle', TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        END IF
# if defined EW_PERIODIC || defined NS_PERIODIC
        CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,             &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          GRID(ng) % angler)
# endif
# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, model, 1, Istr, Iend, Jstr, Jend,         &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    GRID(ng) % angler)
# endif
      END IF
!
! Close GRID NetCDF file.
!
      status=nf90_close(ncGRDid(ng))
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
