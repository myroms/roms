#include "cppdefs.h"
#undef DEBUG
      FUNCTION nf_fread (ng, ncid, ncvarid, tindex, gtype,              &
     &                   LBi, UBi, LBj, UBj, LBk, UBk, Ascl,            &
     &                   Amin, Amax,                                    &
#ifdef MASKING
     &                   Amask,                                         &
#endif
     &                   A)
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This function reads in a generic floating point array from an       !
!  input NetCDF file.                                                  !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number (integer).                         !
!     ncid       NetCDF file ID (integer).                             !
!     ncvarid    NetCDF variable ID (integer).                         !
!     tindex     NetCDF time record index to read (integer).           !
!     gtype      Grid type.                                            !
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
!     nf_fread   Error flag (integer).                                 !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, ncid, ncvarid, tindex, gtype
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk

      real(r8), intent(in) :: Ascl
      real(r8), intent(out) :: Amin
      real(r8), intent(out) :: Amax

#ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:UBi,LBj:UBj)
#endif
      real(r8), intent(out) :: A(LBi:UBi,LBj:UBj,LBk:UBk)
!
!  Local variable declarations.
!
      logical :: is2D
      integer :: i, j, k, ic, jc, kc, Npts, NWpts, status, wtype
      integer :: Imin, Imax, Jmin, Jmax, Kmin, Kmax, Ioff, Joff, Koff
      integer :: Ilen, Jlen, Klen, IJlen, MyType, Nghost

      integer, dimension(4) :: start, total

      integer :: nf_fread

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
!
!  Set first and last grid point according to staggered C-grid
!  classification. Set loops offsets.
!
      Imin=0
      Imax=Lm(ng)+1
      Jmin=0
      Jmax=Mm(ng)+1
      MyType=ABS(gtype)
      Nghost=GHOST_POINTS
      IF ((MyType.eq.p2dvar).or.(MyType.eq.u2dvar).or.                  &
     &    (MyType.eq.p3dvar).or.(MyType.eq.u3dvar)) Imin=1
      IF ((MyType.eq.p2dvar).or.(MyType.eq.v2dvar).or.                  &
     &    (MyType.eq.p3dvar).or.(MyType.eq.v3dvar)) Jmin=1
      IF (Imin.eq.0) THEN
        Ioff=1
      ELSE
        Ioff=0
      END IF
      IF (Jmin.eq.0) THEN
        Joff=0
      ELSE
        Joff=1
      END IF
      IF (LBk.eq.0) THEN
        Koff=0
      ELSE
        Koff=1
      END IF
      Ilen=Imax-Imin+1
      Jlen=Jmax-Jmin+1
      Klen=UBk-LBk+1
      IJlen=Ilen*Jlen
!
!  Determine if processing 2D or 3D fields.
!
      is2D=.false.
      IF (MyType.le.4) is2D=.true.
#if defined READ_WATER && defined MASKING
!
!  If processing water points only, set number of points and type
!  switch.
!
      IF ((MyType.eq.p2dvar).or.(MyType.eq.p3dvar)) THEN
        Npts=Nxyp(ng)
        wtype=p2dvar
      ELSE IF ((MyType.eq.u2dvar).or.(MyType.eq.u3dvar)) THEN
        Npts=Nxyu(ng)
        wtype=u2dvar
      ELSE IF ((MyType.eq.v2dvar).or.(MyType.eq.v3dvar)) THEN
        Npts=Nxyv(ng)
        wtype=v2dvar
      ELSE
        Npts=Nxyr(ng)
        wtype=r2dvar
      END IF
      NWpts=(Lm(ng)+2)*(Mm(ng)+2)
# if !defined INLINE_2DIO && defined DISTRIBUTE
      IF (.not.is2D) Npts=Npts*Klen
# endif
#endif
!
!  Set NetCDF dimension counters for processing requested field.
!
      IF (gtype.gt.0) THEN
        start(1)=1
        total(1)=Ilen
        start(2)=1
        total(2)=Jlen
        Npts=IJlen
        IF (is2D) THEN
          start(3)=tindex
          total(3)=1
        ELSE
          start(3)=1
          total(3)=Klen
          start(4)=tindex
          total(4)=1
# if !defined INLINE_2DIO && defined DISTRIBUTE
          Npts=Npts*Klen
# endif
        END IF
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
#if defined INLINE_2DIO && defined DISTRIBUTE

!  If appropriate, process 3D data level by level to reduce memory
!  requirements.
!
      Kmin=1
      Kmax=1
      DO k=LBk,UBk
        IF (.not.is2D) THEN
          start(3)=k-Koff+1
          total(3)=1
        END IF
#else
        Kmin=LBk
        Kmax=UBk
#endif
        nf_fread=nf_noerr
        IF (InpThread) THEN
          status=nf_get_vara_TYPE(ncid, ncvarid, start, total, wrk)
          nf_fread=status
          IF (status.ne.nf_noerr) RETURN
          Amin=wrk(1)*Ascl
          Amax=wrk(1)*Ascl
          DO i=1,Npts
            wrk(i)=Ascl*wrk(i)
            Amin=MIN(Amin,wrk(i))
            Amax=MAX(Amax,wrk(i))
          END DO
        END IF          
!
!-----------------------------------------------------------------------
!  Unpack read field.
!-----------------------------------------------------------------------
!
#ifdef DISTRIBUTE
        CALL mp_scatter (ng, LBi, UBi, LBj, UBj, Kmin, Kmax,            &
     &                   Nghost, gtype, Amin, Amax,                     &
# if defined READ_WATER && defined MASKING
     &                   NWpts, SCALARS(ng)%IJwater(1,wtype),           &
# endif
# ifdef INLINE_2DIO
     &                   Npts, wrk, A(LBi,LBj,k))
      END DO
# else
     &                   Npts, wrk, A)
# endif
#else
      IF (gtype.gt.0) THEN
        DO k=LBk,UBk
          kc=(k-Koff)*IJlen
          DO j=Jmin,Jmax
            jc=(j-Joff)*Ilen+kc
            DO i=Imin,Imax
              ic=i+Ioff+jc
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
#ifdef DEBUG
      status=nf_inq_varname(ncid, ncvarid, Aname)
      CALL mp_dump (ng, MyRank, gtype, LBi, UBi, LBj, UBj, LBk, UBk, A, &
     &              Aname)
#endif
      RETURN
      END FUNCTION nf_fread
