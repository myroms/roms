#include "cppdefs.h"
      MODULE nf_fread3d_mod
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This function reads in a generic floating point 3D array from an    !
!  input NetCDF file.                                                  !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number (integer)                          !
!     model      Calling model identifier (integer)                    !
!     ncname     NetCDF file name (string)                             !
!     ncid       NetCDF file ID (integer)                              !
!     ncvname    NetCDF variable name (string)                         !
!     ncvarid    NetCDF variable ID (integer)                          !
!     tindex     NetCDF time record index to read (integer)            !
!     gtype      C-grid type (integer)                                 !
!     Vsize      Variable dimensions in NetCDF file (integer 1D array) !
!     LBi        I-dimension Lower bound (integer)                     !
!     UBi        I-dimension Upper bound (integer)                     !
!     LBj        J-dimension Lower bound (integer)                     !
!     UBj        J-dimension Upper bound (integer)                     !
!     LBk        K-dimension Lower bound (integer)                     !
!     UBk        K-dimension Upper bound (integer)                     !
!     Ascl       Factor to scale field after reading (real).           !
!     Amask      Land/Sea mask, if any (real 3D array)                 !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Amin       Field minimum value (real)                            !
!     Amax       Field maximum value (real)                            !
!     A          Field to read in (real 3D array)                      !
!     nf_fread3d Error flag (integer)                                  !
!                                                                      !
!=======================================================================
!
      implicit none

      CONTAINS

#if defined PARALLEL_IO && defined DISTRIBUTE
!
!***********************************************************************
      FUNCTION nf_fread3d (ng, model, ncname, ncid,                     &
     &                     ncvname, ncvarid,                            &
     &                     tindex, gtype, Vsize,                        &
     &                     LBi, UBi, LBj, UBj, LBk, UBk,                &
     &                     Ascl, Amin, Amax,                            &
# ifdef MASKING
     &                     Amask,                                       &
# endif
     &                     A)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
!
      USE distribute_mod, ONLY : mp_bcasti, mp_reduce, mp_barrier
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, ncid, ncvarid, tindex, gtype
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
      integer, intent(in) :: Vsize(4)

      real(r8), intent(in)  :: Ascl
      real(r8), intent(out) :: Amin
      real(r8), intent(out) :: Amax

      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: ncvname

# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:,LBj:)
#  endif
      real(r8), intent(out) :: A(LBi:,LBj:,LBk:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(out) :: A(LBi:UBi,LBj:UBj,LBk:UBk)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k, ic, Npts, NWpts
      integer :: Imin, Imax, Isize, Jmin, Jmax, Jsize
      integer :: Ioff, Joff, Koff
      integer :: Ilen, Jlen, Klen, IJlen
      integer :: Cgrid, MyType, ghost, status, wtype

      integer, dimension(4) :: start, total

      integer :: nf_fread3d

      real(r8) :: BIG = 1.0E+37_r8
      real(r8) :: Aval

      real(r8), dimension(2) :: buffer

      real(r8), dimension(TileSize(ng)*(UBk-LBk+1)) :: wrk

      character (len=3), dimension(2) :: op_handle
!
!-----------------------------------------------------------------------
!  Set starting and ending indices to process.
!-----------------------------------------------------------------------
!
!  Set first and last grid point according to staggered C-grid
!  classification. Set the offsets for variables with starting 
!  zero-index.  Recall the NetCDF does not support a zero-index.
!
!  Notice that ghost=0, so (Imin,Jmin) and (Imax,Jmax) are the
!  corner of the computational tile and the ghost points are
!  not processed. They will be processed elsewhere by the 
!  appropriate call to any of the routines in "mp_exchange.F".
!
      MyType=gtype

      ghost=0                        ! non-overlapping, no ghost points

      SELECT CASE (ABS(MyType))
        CASE (p2dvar, p3dvar)
          Cgrid=1
          Ioff=0
          Joff=0
          Isize=IOBOUNDS(ng)%IUB_psi-IOBOUNDS(ng)%ILB_psi+1
          Jsize=IOBOUNDS(ng)%JUB_psi-IOBOUNDS(ng)%JLB_psi+1
        CASE (r2dvar, r3dvar)
          Cgrid=2
          Ioff=1
          Joff=1
          Isize=IOBOUNDS(ng)%IUB_rho-IOBOUNDS(ng)%ILB_rho+1
          Jsize=IOBOUNDS(ng)%JUB_rho-IOBOUNDS(ng)%JLB_rho+1
        CASE (u2dvar, u3dvar)
          Cgrid=3
          Ioff=0
          Joff=1
          Isize=IOBOUNDS(ng)%IUB_u-IOBOUNDS(ng)%ILB_u+1
          Jsize=IOBOUNDS(ng)%JUB_u-IOBOUNDS(ng)%JLB_u+1
        CASE (v2dvar, v3dvar)
          Cgrid=4
          Ioff=1
          Joff=0
          Isize=IOBOUNDS(ng)%IUB_u-IOBOUNDS(ng)%ILB_v+1
          Jsize=IOBOUNDS(ng)%JUB_u-IOBOUNDS(ng)%JLB_v+1
        CASE DEFAULT
          Cgrid=1
          Ioff=1
          Joff=1
          Isize=IOBOUNDS(ng)%IUB_rho-IOBOUNDS(ng)%ILB_rho+1
          Jsize=IOBOUNDS(ng)%JUB_rho-IOBOUNDS(ng)%JLB_rho+1
      END SELECT

      Imin=BOUNDS(ng)%Imin(Cgrid,ghost,MyRank)
      Imax=BOUNDS(ng)%Imax(Cgrid,ghost,MyRank)
      Jmin=BOUNDS(ng)%Jmin(Cgrid,ghost,MyRank)
      Jmax=BOUNDS(ng)%Jmax(Cgrid,ghost,MyRank)

      Ilen=Imax-Imin+1
      Jlen=Jmax-Jmin+1
      Klen=UBk-LBk+1
      IJlen=Ilen*Jlen

      IF (LBk.eq.0) THEN
        Koff=1
      ELSE
        Koff=0
      END IF

# if defined READ_WATER && defined MASKING
!
!  If processing water points only, set number of points and type
!  switch.
!
      SELECT CASE (ABS(MyType))
        CASE (p3dvar)
          Npts=IOBOUNDS(ng)%xy_psi
          wtype=p2dvar
        CASE (r3dvar)
          Npts=IOBOUNDS(ng)%xy_rho
          wtype=r2dvar
        CASE (u3dvar)
          Npts=IOBOUNDS(ng)%xy_u
          wtype=u2dvar
        CASE (v3dvar)
          Npts=IOBOUNDS(ng)%xy_v
          wtype=v2dvar
        CASE DEFAULT
          Npts=IOBOUNDS(ng)%xy_rho
          wtype=r2dvar
      END SELECT
      NWpts=(Lm(ng)+2)*(Mm(ng)+2)
      Npts=Npts*Klen
# endif
!
!  Set NetCDF dimension counters for processing requested field.
!
      IF (MyType.gt.0) THEN
        start(1)=Imin+Ioff
        total(1)=Ilen
        start(2)=Jmin+Joff
        total(2)=Jlen
        start(3)=LBk+Koff
        total(3)=Klen
        start(4)=tindex
        total(4)=1
        Npts=IJlen*Klen
# if defined READ_WATER && defined MASKING
      ELSE
        start(1)=1
        total(1)=Npts
        start(2)=1
        total(2)=tindex
# endif
      END IF
!
!  Initialize local array to avoid denormalized numbers. This
!  facilitates processing and debugging.
!
      wrk=0.0_r8
!
!-----------------------------------------------------------------------
!  Parallel I/O: Read in tile data from requested field and scale it.
!-----------------------------------------------------------------------
!
      Amin=BIG
      Amax=-BIG

      status=nf90_get_var(ncid, ncvarid, wrk, start, total)
      nf_fread3d=status

      IF (status.eq.nf90_noerr) THEN
        DO i=1,Npts
          IF (wrk(i).ge.spval_check) THEN
            wrk(i)=0.0_r8               ! land/sea maked with _FillValue
          ELSE
            wrk(i)=Ascl*wrk(i)
            Amin=MIN(Amin,wrk(i))
            Amax=MAX(Amax,wrk(i))
          END IF
        END DO
      ELSE
        RETURN
      END IF          
!
!-----------------------------------------------------------------------
!  Parallel I/O: Unpack read tile data block and load into model array.
!-----------------------------------------------------------------------
!
      IF (MyType.gt.0) THEN
        ic=0
        DO k=LBk,UBk
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
              A(i,j,k)=wrk(ic)
            END DO
          END DO
        END DO
# if defined MASKING || defined READ_WATER
      ELSE
        ic=0
        DO k=LBk,UBk
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              IF (Amask(i,j).gt.0.0_r8) THEN
                ic=ic+1
                A(i,j,k)=wrk(ic)
              ELSE
                A(i,j,k)=0.0_r8
              END IF
            END DO
          END DO
        END DO
# endif
      END IF
!
!-----------------------------------------------------------------------
!  Set global minimun and maximun values: global reduction.
!-----------------------------------------------------------------------
!
      buffer(1)=Amin
      op_handle(1)='MIN'
      buffer(2)=Amax
      op_handle(2)='MAX'
      CALL mp_reduce (ng, model, 2, buffer, op_handle)
      Amin=buffer(1)
      Amax=buffer(2)

      RETURN
      END FUNCTION nf_fread3d

#else

!
!***********************************************************************
      FUNCTION nf_fread3d (ng, model, ncname, ncid,                     &
     &                     ncvname, ncvarid,                            &
     &                     tindex, gtype, Vsize,                        &
     &                     LBi, UBi, LBj, UBj, LBk, UBk,                &
     &                     Ascl, Amin, Amax,                            &
# ifdef MASKING
     &                     Amask,                                       &
# endif
     &                     A)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars

# ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcasti
#  ifdef INLINE_2DIO
      USE distribute_mod, ONLY : mp_scatter2d
#  else
      USE distribute_mod, ONLY : mp_scatter3d
#  endif
# endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, ncid, ncvarid, tindex, gtype
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
      integer, intent(in) :: Vsize(4)

      real(r8), intent(in)  :: Ascl
      real(r8), intent(out) :: Amin
      real(r8), intent(out) :: Amax

      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: ncvname

# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:,LBj:)
#  endif
      real(r8), intent(out) :: A(LBi:,LBj:,LBk:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(out) :: A(LBi:UBi,LBj:UBj,LBk:UBk)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k, ic, Npts, NWpts, status, wtype
      integer :: Imin, Imax, Jmin, Jmax, Koff
      integer :: Ilen, Jlen, Klen, IJlen, MyType
# ifdef DISTRIBUTE
      integer :: Nghost
# endif
      integer, dimension(4) :: start, total

      integer :: nf_fread3d

      real(r8) :: BIG = 1.0E+37_r8
      real(r8) :: Aval

# if defined INLINE_2DIO && defined DISTRIBUTE
      real(r8), dimension(2+(Lm(ng)+2)*(Mm(ng)+2)) :: wrk
# else
      real(r8), dimension(2+(Lm(ng)+2)*(Mm(ng)+2)*(UBk-LBk+1)) :: wrk
# endif
!
!-----------------------------------------------------------------------
!  Set starting and ending indices to process.
!-----------------------------------------------------------------------
!
!  Set first and last grid point according to staggered C-grid
!  classification. Set loops offsets.
!
      MyType=gtype

      SELECT CASE (ABS(MyType))
        CASE (p2dvar, p3dvar)
          Imin=IOBOUNDS(ng)%ILB_psi
          Imax=IOBOUNDS(ng)%IUB_psi
          Jmin=IOBOUNDS(ng)%JLB_psi
          Jmax=IOBOUNDS(ng)%JUB_psi
        CASE (r2dvar, r3dvar)
          Imin=IOBOUNDS(ng)%ILB_rho
          Imax=IOBOUNDS(ng)%IUB_rho
          Jmin=IOBOUNDS(ng)%JLB_rho
          Jmax=IOBOUNDS(ng)%JUB_rho
        CASE (u2dvar, u3dvar)
          Imin=IOBOUNDS(ng)%ILB_u
          Imax=IOBOUNDS(ng)%IUB_u
          Jmin=IOBOUNDS(ng)%JLB_u
          Jmax=IOBOUNDS(ng)%JUB_u
        CASE (v2dvar, v3dvar)
          Imin=IOBOUNDS(ng)%ILB_v
          Imax=IOBOUNDS(ng)%IUB_v
          Jmin=IOBOUNDS(ng)%JLB_v
          Jmax=IOBOUNDS(ng)%JUB_v
        CASE DEFAULT
          Imin=IOBOUNDS(ng)%ILB_rho
          Imax=IOBOUNDS(ng)%IUB_rho
          Jmin=IOBOUNDS(ng)%JLB_rho
          Jmax=IOBOUNDS(ng)%JUB_rho
      END SELECT

      Ilen=Imax-Imin+1
      Jlen=Jmax-Jmin+1
      Klen=UBk-LBk+1
      IJlen=Ilen*Jlen

      IF (LBk.eq.0) THEN
        Koff=0
      ELSE
        Koff=1
      END IF

# ifdef DISTRIBUTE
!
!  Notice that Nghost=0, so (Imin,Jmin) and (Imax,Jmax) are the
!  corner of the computational tile and the ghost points are
!  not processed. They will be processed elsewhere by the 
!  appropriate call to any of the routines in "mp_exchange.F".
!
      Nghost=0
# endif
# if defined READ_WATER && defined MASKING
!
!  If processing water points only, set number of points and type
!  switch.
!
      SELECT CASE (ABS(MyType))
        CASE (p3dvar)
          Npts=IOBOUNDS(ng)%xy_psi
          wtype=p2dvar
        CASE (r3dvar)
          Npts=IOBOUNDS(ng)%xy_rho
          wtype=r2dvar
        CASE (u3dvar)
          Npts=IOBOUNDS(ng)%xy_u
          wtype=u2dvar
        CASE (v3dvar)
          Npts=IOBOUNDS(ng)%xy_v
          wtype=v2dvar
        CASE DEFAULT
          Npts=IOBOUNDS(ng)%xy_rho
          wtype=r2dvar
      END SELECT
      NWpts=(Lm(ng)+2)*(Mm(ng)+2)
#  if !(defined INLINE_2DIO && defined DISTRIBUTE)
      Npts=Npts*Klen
#  endif
# endif
!
!  Set NetCDF dimension counters for processing requested field.
!
      IF (MyType.gt.0) THEN
        start(1)=1
        total(1)=Ilen
        start(2)=1
        total(2)=Jlen
        start(3)=1
        total(3)=Klen
        start(4)=tindex
        total(4)=1
        Npts=IJlen
# if !(defined INLINE_2DIO && defined DISTRIBUTE)
        Npts=Npts*Klen
# endif
# if defined READ_WATER && defined MASKING
      ELSE
        start(1)=1
        total(1)=Npts
        start(2)=1
        total(2)=tindex
# endif
      END IF
!
!  Initialize local array to avoid denormalized numbers. This
!  facilitates processing and debugging.
!
      wrk=0.0_r8
!
!-----------------------------------------------------------------------
!  Serial I/O: Read in requested field and scale it.
!-----------------------------------------------------------------------
!
      Amin=BIG
      Amax=-BIG

# if defined INLINE_2DIO && defined DISTRIBUTE
!
!  If appropriate, process 3D data level by level to reduce memory
!  requirements.
!
      DO k=LBk,UBk
        start(3)=k-Koff+1
        total(3)=1
# endif
        status=nf90_noerr
        IF (InpThread) THEN
          status=nf90_get_var(ncid, ncvarid, wrk, start, total)
          IF (status.eq.nf90_noerr) THEN
            DO i=1,Npts
              IF (wrk(i).ge.spval_check) THEN
                wrk(i)=0.0_r8           ! land/sea maked with _FillValue
              ELSE
                wrk(i)=Ascl*wrk(i)
                Amin=MIN(Amin,wrk(i))
                Amax=MAX(Amax,wrk(i))
              END IF
            END DO
          END IF
        END IF          
# ifdef DISTRIBUTE
        CALL mp_bcasti (ng, model, status)
# endif
        IF (status.ne.nf90_noerr) THEN
          nf_fread3d=status
          RETURN
        END IF
!
!-----------------------------------------------------------------------
!  Serial I/O: Unpack read field.
!-----------------------------------------------------------------------
!
# ifdef DISTRIBUTE
#  ifdef INLINE_2DIO
        CALL mp_scatter2d (ng, model, LBi, UBi, LBj, UBj,               &
     &                     Nghost, MyType, Amin, Amax,                  &
#   if defined READ_WATER && defined MASKING
     &                     NWpts, SCALARS(ng)%IJwater(:,wtype),         &
#   endif
     &                     Npts, wrk, A(:,:,k))
      END DO
#  else
        CALL mp_scatter3d (ng, model, LBi, UBi, LBj, UBj, LBk, UBk,     &
     &                     Nghost, MyType, Amin, Amax,                  &
#   if defined READ_WATER && defined MASKING
     &                     NWpts, SCALARS(ng)%IJwater(:,wtype),         &
#   endif
     &                     Npts, wrk, A)
#  endif
# else        
      IF (MyType.gt.0) THEN
        ic=0
        DO k=LBk,UBk
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
              A(i,j,k)=wrk(ic)
            END DO
          END DO
        END DO
#  if defined MASKING || defined READ_WATER
      ELSE
        ic=0
        DO k=LBk,UBk
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              IF (Amask(i,j).gt.0.0_r8) THEN
                ic=ic+1
                A(i,j,k)=wrk(ic)
              ELSE
                A(i,j,k)=0.0_r8
              END IF
            END DO
          END DO
        END DO
#  endif
      END IF
# endif

      nf_fread3d=status

      RETURN
      END FUNCTION nf_fread3d
#endif
      END MODULE nf_fread3d_mod
