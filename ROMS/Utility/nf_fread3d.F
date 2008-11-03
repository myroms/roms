#include "cppdefs.h"
#undef DEBUG
      FUNCTION nf_fread3d (ng, model, ncid, ncvarid, tindex, gtype,     &
     &                     Vsize, LBi, UBi, LBj, UBj, LBk, UBk, Ascl,   &
     &                     Amin, Amax,                                  &
#ifdef MASKING
     &                     Amask,                                       &
#endif
     &                     A)
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This function reads in a generic floating point 3D array from an    !
!  input NetCDF file.                                                  !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     ncid       NetCDF file ID.                                       !
!     ncvarid    NetCDF variable ID.                                   !
!     tindex     NetCDF time record index to read (integer).           !
!     gtype      Grid type.                                            !
!     Vsize      Variable dimensions in NetCDF file.                   !
!     LBi        I-dimension Lower bound.                              !
!     UBi        I-dimension Upper bound.                              !
!     LBj        J-dimension Lower bound.                              !
!     UBj        J-dimension Upper bound.                              !
!     LBk        K-dimension Lower bound.                              !
!     UBk        K-dimension Upper bound.                              !
!     Ascl       Factor to scale field after reading (real).           !
!     Amask      Land/Sea mask, if any (real).                         !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Amin       Field minimum value (real).                           !
!     Amax       Field maximum value (real).                           !
!     A          Field to read in (real).                              !
#ifdef DISTRIBUTE
!     gtype      Grid type.                                            !
#endif
!     nf_fread3d Error flag (integer).                                 !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars

#ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcasti, mp_scatter
# ifdef DEBUG
      USE distribute_mod, ONLY : mp_dump
# endif
#endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, ncid, ncvarid, tindex
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
      integer, intent(in) :: Vsize(4)

#ifndef DISTRIBUTE
      integer, intent(in) :: gtype
#else
      integer, intent(inout) :: gtype
#endif

      real(r8), intent(in)  :: Ascl
      real(r8), intent(out) :: Amin
      real(r8), intent(out) :: Amax

#ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:UBi,LBj:UBj)
#endif
      real(r8), intent(out) :: A(LBi:UBi,LBj:UBj,LBk:UBk)
!
!  Local variable declarations.
!
      integer :: i, j, k, ic, Npts, NWpts, status, wtype
      integer :: Imin, Imax, Jmin, Jmax, Kmin, Kmax, Koff
      integer :: Ilen, Jlen, Klen, IJlen, MyType, Nghost

      integer, dimension(4) :: start, total

      integer :: nf_fread3d

      real(r8) :: BIG = 1.0E+37_r8
      real(r8) :: Aval

#if defined INLINE_2DIO && defined DISTRIBUTE
      real(r8), dimension(2+(Lm(ng)+2)*(Mm(ng)+2)) :: wrk
#else
      real(r8), dimension(2+(Lm(ng)+2)*(Mm(ng)+2)*(UBk-LBk+1)) :: wrk
#endif

#ifdef DEBUG
      character (len=20) :: Aname
#endif
!
!-----------------------------------------------------------------------
!  Set starting and ending indices to process.
!-----------------------------------------------------------------------

#ifdef DISTRIBUTE
!
!  In some parts of the code the IO is serial and the value of the
!  grid type is only known by the input thread so broadcast its value
!  to all nodes.
!
      CALL mp_bcasti (ng, model, gtype, 1)
#endif
!
!  Set first and last grid point according to staggered C-grid
!  classification. Set loops offsets.  Notice that Nghost is set
!  to zero when processing an adjoint solution.  Therefore, the
!  ghost points will be not assigned in mp_scatter. This is the
!  correct adjoint solution reading.  The ghost points are then
!  assumed to be zero.
!
      MyType=gtype

      SELECT CASE (ABS(MyType))
        CASE (p3dvar)
          Imin=IOBOUNDS(ng)%ILB_psi
          Imax=IOBOUNDS(ng)%IUB_psi
          Jmin=IOBOUNDS(ng)%JLB_psi
          Jmax=IOBOUNDS(ng)%JUB_psi
        CASE (r3dvar)
          Imin=IOBOUNDS(ng)%ILB_rho
          Imax=IOBOUNDS(ng)%IUB_rho
          Jmin=IOBOUNDS(ng)%JLB_rho
          Jmax=IOBOUNDS(ng)%JUB_rho
        CASE (u3dvar)
          Imin=IOBOUNDS(ng)%ILB_u
          Imax=IOBOUNDS(ng)%IUB_u
          Jmin=IOBOUNDS(ng)%JLB_u
          Jmax=IOBOUNDS(ng)%JUB_u
        CASE (v3dvar)
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

      IF (model.eq.iADM) THEN
        Nghost=0
      ELSE  
        Nghost=GHOST_POINTS
      END IF

#if defined READ_WATER && defined MASKING
!
!  If processing water points only, set number of points and type
!  switch.
!
      IF (ABS(MyType).eq.p3dvar) THEN
        Npts=IOBOUNDS(ng)%xy_psi
        wtype=p2dvar
      ELSE IF (ABS(MyType).eq.u3dvar) THEN
        Npts=IOBOUNDS(ng)%xy_u
        wtype=u2dvar
      ELSE IF (ABS(MyType).eq.v3dvar) THEN
        Npts=IOBOUNDS(ng)%xy_v
        wtype=v2dvar
      ELSE
        Npts=IOBOUNDS(ng)%xy_rho
        wtype=r2dvar
      END IF
      NWpts=(Lm(ng)+2)*(Mm(ng)+2)
# if !(defined INLINE_2DIO && defined DISTRIBUTE)
      Npts=Npts*Klen
# endif
#endif
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
#if !(defined INLINE_2DIO && defined DISTRIBUTE)
        Npts=Npts*Klen
#endif
#if defined READ_WATER && defined MASKING
      ELSE
        start(1)=1
        total(1)=Npts
        start(2)=1
        total(2)=tindex
#endif
      END IF
!
!-----------------------------------------------------------------------
!  Read in requested field and scale it.
!-----------------------------------------------------------------------
!
      Amin=BIG
      Amax=-BIG

#if defined INLINE_2DIO && defined DISTRIBUTE
!
!  If appropriate, process 3D data level by level to reduce memory
!  requirements.
!
      Kmin=1
      Kmax=1
      DO k=LBk,UBk
        start(3)=k-Koff+1
        total(3)=1
#else
        Kmin=LBk
        Kmax=UBk
#endif
        status=nf90_noerr
        IF (InpThread) THEN
          status=nf90_get_var(ncid, ncvarid, wrk, start, total)
          IF (status.eq.nf90_noerr) THEN
            DO i=1,Npts
              wrk(i)=Ascl*wrk(i)
              IF (wrk(i).ge.spval) THEN
                wrk(i)=0.0_r8           ! land/sea maked with _FillValue
              ELSE
                Amin=MIN(Amin,wrk(i))
                Amax=MAX(Amax,wrk(i))
              END IF
            END DO
          END IF
        END IF          
#ifdef DISTRIBUTE
        CALL mp_bcasti (ng, model, status, 1)
#endif
        IF (status.ne.nf90_noerr) THEN
          nf_fread3d=status
          RETURN
        END IF
!
!-----------------------------------------------------------------------
!  Unpack read field.
!-----------------------------------------------------------------------
!
#ifdef DISTRIBUTE
        CALL mp_scatter (ng, model, LBi, UBi, LBj, UBj, Kmin, Kmax,     &
     &                   Nghost, MyType, Amin, Amax,                    &
# if defined READ_WATER && defined MASKING
     &                   NWpts, SCALARS(ng)%IJwater(:,wtype),           &
# endif
# ifdef INLINE_2DIO
     &                   Npts, wrk, A(LBi,LBj,k))
      END DO
# else
     &                   Npts, wrk, A)
# endif
#else        
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
#endif
#if defined DISTRIBUTE && defined DEBUG
!
!-----------------------------------------------------------------------
!  If debugging, write distributed data into formatted files.
!-----------------------------------------------------------------------
!
      status=nf90_inquire_variable(ncid, ncvarid, Aname)
      CALL mp_dump (ng, MyRank, MyType, LBi, UBi, LBj, UBj, LBk, UBk,   &
     &              A, Aname)
#endif
      nf_fread3d=status

      RETURN
      END FUNCTION nf_fread3d
