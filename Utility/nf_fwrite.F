#include "cppdefs.h"
      FUNCTION nf_fwrite (ng, ncid, ncvarid, tindex, gtype,             &
     &                    LBi, UBi, LBj, UBj, LBk, UBk, Ascl,           &
#ifdef MASKING
     &                    Amask,                                        &
#endif
     &                    A)
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This function writes out a generic floating point array into an     !
!  output NetCDF file.                                                 !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number (integer).                         !
!     ncid       NetCDF file ID (integer).                             !
!     ncvarid    NetCDF variable ID (integer).                         !
!     tindex     NetCDF time record index to write (integer).          !
!     gtype      Grid type. If negative, only write water points.      !
!     LBi        I-dimension Lower bound.                              !
!     UBi        I-dimension Upper bound.                              !
!     LBj        J-dimension Lower bound.                              !
!     UBj        J-dimension Upper bound.                              !
!     LBk        K-dimension Lower bound.                              !
!     UBk        K-dimension Upper bound.                              !
!     Amask      land/Sea mask, if any (real).                         !
!     Ascl       Factor to scale field before writing (real).          !
!     A          Field to write out (real).                            !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     nf_fwrite  Error flag (integer).                                 !
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

#ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:UBi,LBj:UBj)
#endif
      real(r8), intent(in) :: A(LBi:UBi,LBj:UBj,LBk:UBk)
!
!  Local variable declarations.
!
      logical :: is2D
      integer :: i, j, k, ic, jc, kc, Npts
      integer :: Imin, Imax, Jmin, Jmax, Kmin, Kmax, Ioff, Joff, Koff
      integer :: Ilen, Jlen, Klen, IJlen, MyType

      integer, dimension(4) :: start, total

      integer :: nf_fwrite

      real(r8) :: Amin, Amax

#if defined INLINE_2DIO && defined DISTRIBUTE
      real(r8), dimension((Lm(ng)+2)*(Mm(ng)+2)) :: Aout
#else
      real(r8), dimension((Lm(ng)+2)*(Mm(ng)+2)*(UBk-LBk+1)) :: Aout
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

#ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  If distributed-memory set-up, collect tile data from all spawned
!  nodes and store it into a global scratch 1D array, packed in column-
!  major order. If applicable, remove land points.
!-----------------------------------------------------------------------
!
# ifdef INLINE_2DIO

!  If appropriate, process 3D data level by level to reduce memory
!  requirements.
!
      Kmin=1
      Kmax=1
      DO k=LBk,UBk
# else
        Kmin=LBk
        Kmax=UBk
# endif
        CALL mp_gather (ng, LBi, UBi, LBj, UBj, Kmin, Kmax, gtype,      &
     &                  Ascl,                                           &
# ifdef MASKING
     &                  Amask,                                          &
# endif
# ifdef INLINE_2DIO
     &                  A(LBi,LBj,k), Npts, Aout)
# else
     &                  A, Npts, Aout)
# endif
!
!-----------------------------------------------------------------------
!  Write output buffer into NetCDF file.
!-----------------------------------------------------------------------
!
        nf_fwrite=nf_noerr
        IF (OutThread) THEN
          IF (gtype.gt.0) THEN
            start(1)=1
            total(1)=Ilen
            start(2)=1
            total(2)=Jlen
            IF (is2D) THEN
              start(3)=tindex
              total(3)=1
            ELSE
# ifdef INLINE_2DIO
              start(3)=k-Koff+1
              total(3)=1
# else
              start(3)=1
              total(3)=Klen
# endif
              start(4)=tindex
              total(4)=1
            END IF
# ifdef MASKING
          ELSE
            start(1)=1
            total(1)=Npts
            start(2)=tindex
            total(2)=1
# endif
          END IF
          nf_fwrite=nf_put_vara_TYPE(ncid, ncvarid, start, total, Aout)
        END IF
# ifdef INLINE_2DIO
      END DO
# endif
#else
!
!-----------------------------------------------------------------------
!  If serial or shared-memory applications and serial output, pack data
!  into a global 1D array in column-major order. If applicable, remove
!  land points.
!-----------------------------------------------------------------------
!
      IF (gtype.gt.0) THEN
        DO k=LBk,UBk
          kc=(k-Koff)*IJlen
          DO j=Jmin,Jmax
            jc=(j-Joff)*Ilen+kc
            DO i=Imin,Imax
              ic=i+Ioff+jc
              Aout(ic)=A(i,j,k)*Ascl
            END DO
          END DO
        END DO
# ifdef MASKING
      ELSE
        Npts=0
        DO k=LBk,UBk
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              IF (Amask(i,j).gt.0.0_r8) THEN
                Npts=Npts+1
                Aout(Npts)=A(i,j,k)*Ascl
              END IF
            END DO
          END DO
        END DO
# endif
      END IF
!
!-----------------------------------------------------------------------
!  Write output buffer into NetCDF file.
!-----------------------------------------------------------------------
!
      nf_fwrite=nf_noerr
      IF (OutThread) THEN
        IF (gtype.gt.0) THEN
          start(1)=1
          total(1)=Ilen
          start(2)=1
          total(2)=Jlen
          IF (is2D) THEN
            start(3)=tindex
            total(3)=1
          ELSE
            start(3)=1
            total(3)=Klen
            start(4)=tindex
            total(4)=1
          END IF
# ifdef MASKING
        ELSE
          start(1)=1
          total(1)=Npts
          start(2)=tindex
          total(2)=1
# endif
        END IF
        nf_fwrite=nf_put_vara_TYPE(ncid, ncvarid, start, total, Aout)
      END IF
#endif
      RETURN
      END FUNCTION nf_fwrite
