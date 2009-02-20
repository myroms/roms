#include "cppdefs.h"
      MODULE nf_fread3d_bry_mod
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This function writes out a generic floating point 3D boundary array !
!  into an output NetCDF file.                                         !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF output file name (string)                    !
!     ncvname      NetCDF variable name (string)                       !
!     ncid         NetCDF file ID (integer)                            !
!     ncvarid      NetCDF variable ID (integer)                        !
!     tindex       NetCDF time record index to write (integer)         !
!     gtype        Grid type (integer)                                 !
!     LBij         IJ-dimension Lower bound (integer)                  !
!     UBij         IJ-dimension Upper bound (integer)                  !
!     LBk          K-dimension lower bound (integer)                   !
!     UBk          K-dimension upper bound (integer)                   !
!     Nrec         Number of boundary records (integer)                !
!     Ascl         Factor to scale field before writing (real)         !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Amin         Field minimum value (real)                          !
!     Amax         Field maximum value (real)                          !
!     A            3D boundary field to read in (real array)           !
!     status       Error flag (integer)                                !
!                                                                      !
!=======================================================================
!
      implicit none

      CONTAINS

!
!***********************************************************************
      FUNCTION nf_fread3d_bry (ng, model, ncname, ncid,                 &
     &                         ncvname, ncvarid,                        &
     &                         tindex, gtype,                           &
     &                         LBij, UBij, LBk, UBk, Nrec,              &
     &                         Ascl, Amin, Amax,                        &
     &                         A)  RESULT(status)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars

#if !defined PARALLEL_IO && defined DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcastf
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, ncid, ncvarid, tindex, gtype
      integer, intent(in) :: LBij, UBij, LBk, UBk, Nrec

      real(r8), intent(in)  :: Ascl
      real(r8), intent(out) :: Amin
      real(r8), intent(out) :: Amax

      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: ncvname

#ifdef ASSUMED_SHAPE
      real(r8), intent(out) :: A(LBij:,:,:,:)
#else
      real(r8), intent(out) :: A(LBij:UBij,LBk:UBk,4,Nrec)
#endif
!
!  Local variable declarations.
!
      logical, dimension(4) :: bounded

      integer :: bc, ghost, i, ib, ic, ir, j, k, kc, rc, tile
      integer :: IorJ, IJKlen, Imin, Imax, Jmin, Jmax, Klen, Npts
      integer :: Cgrid, Istr, Iend, Jstr, Jend
      integer, dimension(5) :: start, total

      integer :: status

#if !defined PARALLEL_IO && defined DISTRIBUTE
      real(r8), dimension(3) :: buffer
#endif
      real(r8), dimension((UBij-LBij+1)*(UBk-LBk+1)*4*Nrec) :: wrk
!
!-----------------------------------------------------------------------
!  Set starting and ending indices to process.
!-----------------------------------------------------------------------
!
!  Set first and last grid point according to staggered C-grid
!  classification.
!
!  Notice that (Imin,Jmin) and (Imax,Jmax) are the corner of the
!  computational tile. If ghost=0, ghost points are not processed.
!  They will be processed elsewhere by the appropriate call to any
!  of the routines in "mp_exchange.F".  If ghost=1, the ghost points
!  are read.
!
      IF (model.eq.iADM) THEN
        ghost=0                      ! non-overlapping, no ghost points
      ELSE
        ghost=1                      ! overlapping, read ghost points
      END IF

      SELECT CASE (gtype)
        CASE (p2dvar, p3dvar)
          Cgrid=1
        CASE (r2dvar, r3dvar)
          Cgrid=2
        CASE (u2dvar, u3dvar)
          Cgrid=3
        CASE (v2dvar, v3dvar)
          Cgrid=4
        CASE DEFAULT
          Cgrid=2
      END SELECT

#ifdef DISTRIBUTE
      tile=MyRank
      Imin=BOUNDS(ng)%Imin(Cgrid,ghost,tile)
      Imax=BOUNDS(ng)%Imax(Cgrid,ghost,tile)
      Jmin=BOUNDS(ng)%Jmin(Cgrid,ghost,tile)
      Jmax=BOUNDS(ng)%Jmax(Cgrid,ghost,tile)
#else
      tile=-1
      Imin=LBij
      Imax=UBij
      Jmin=LBij
      Jmax=UBij
#endif
      IorJ=IOBOUNDS(ng)%IorJ
      Klen=UBk-LBk+1
      IJKlen=IorJ*Klen
      Npts=IJKlen*4*Nrec
!
!  Get tile bounds.
!
      Istr=BOUNDS(ng)%Istr (tile)
      Iend=BOUNDS(ng)%Iend (tile)
      Jstr=BOUNDS(ng)%Jstr (tile)
      Jend=BOUNDS(ng)%Jend (tile)
!
!  Set switch to process boundary data by their associated tiles.
!
      bounded(iwest )=WESTERN_EDGE
      bounded(ieast )=EASTERN_EDGE
      bounded(isouth)=SOUTHERN_EDGE
      bounded(inorth)=NORTHERN_EDGE
!
!  Set NetCDF dimension counters for processing requested field.
!
      start(1)=1
      total(1)=IorJ
      start(2)=1
      total(2)=Klen
      start(3)=1
      total(3)=4
      start(4)=1
      total(4)=Nrec
      start(5)=tindex
      total(5)=1
!
!-----------------------------------------------------------------------
!  Read in requested data and scale it.
!-----------------------------------------------------------------------
!
      status=nf90_noerr
      IF (InpThread) THEN
        status=nf90_get_var(ncid, ncvarid, wrk, start, total)
        IF (status.eq.nf90_noerr) THEN
          Amin=spval
          Amax=-spval
          DO i=1,Npts
            IF (wrk(i).ge.spval_check) THEN
              wrk(i)=0.0_r8
            ELSE
              wrk(i)=Ascl*wrk(i)
              Amin=MIN(Amin,wrk(i))
              Amax=MAX(Amax,wrk(i))
            END IF
          END DO
        END IF
      END IF

#if !defined PARALLEL_IO && defined DISTRIBUTE
      buffer(1)=REAL(status,r8)
      buffer(2)=Amin
      buffer(3)=Amax
      CALL mp_bcastf (ng, model, buffer)
      status=INT(buffer(1))
      Amin=buffer(2)
      Amax=buffer(3)
#endif

      IF (status.ne.nf90_noerr) THEN
        exit_flag=2
        ioerror=status
        RETURN
      END IF

#if !defined PARALLEL_IO && defined DISTRIBUTE
!
!  Broadcast data to all spawned nodes.
!
      CALL mp_bcastf (ng, model, wrk)
#endif
!
!-----------------------------------------------------------------------
!  Unpack read data.
!-----------------------------------------------------------------------
!
      A=0.0_r8

      DO ir=1,Nrec
        rc=(ir-1)*IJKlen*4
        DO ib=1,4
          IF (bounded(ib)) THEN
            bc=(ib-1)*IJKlen+rc
            IF ((ib.eq.iwest).or.(ib.eq.ieast)) THEN
              DO k=LBk,UBk
                kc=(k-LBk)*IorJ+bc
                DO j=Jmin,Jmax
                  ic=1+(j-LBij)+kc
                  A(j,k,ib,ir)=wrk(ic)
                END DO
              END DO
            ELSE IF ((ib.eq.isouth).or.(ib.eq.inorth)) THEN
              DO k=LBk,UBk
                kc=(k-LBk)*IorJ+bc
                DO i=Imin,Imax
                  ic=1+(i-LBij)+kc
                  A(i,k,ib,ir)=wrk(ic)
                END DO
              END DO
            END IF
          END IF
        END DO
      END DO
  
      RETURN
      END FUNCTION nf_fread3d_bry

      END MODULE nf_fread3d_bry_mod
