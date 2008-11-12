#include "cppdefs.h"
      MODULE nf_fwrite2d_mod
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This function writes out a generic floating point 2D array into an  !
!  output NetCDF file.                                                 !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number.                                 !
!     model        Calling model identifier.                           !
!     ncid         NetCDF file ID.                                     !
!     ncvarid      NetCDF variable ID.                                 !
!     tindex       NetCDF time record index to write.                  !
!     gtype        Grid type. If negative, only write water points.    !
!     LBi          I-dimension Lower bound.                            !
!     UBi          I-dimension Upper bound.                            !
!     LBj          J-dimension Lower bound.                            !
!     UBj          J-dimension Upper bound.                            !
!     Amask        land/Sea mask, if any (real).                       !
!     Ascl         Factor to scale field before writing (real).        !
!     A            Field to write out (real).                          !
!     SetFillVal   Logical switch to set fill value in land areas      !
!                    (optional).                                       !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     nf_fwrite  Error flag (integer).                                 !
!                                                                      !
!=======================================================================
!
      implicit none

      CONTAINS
!
!***********************************************************************
      FUNCTION nf_fwrite2d (ng, model, ncid, ncvarid, tindex, gtype,    &
     &                      LBi, UBi, LBj, UBj, Ascl,                   &
#ifdef MASKING
     &                      Amask,                                      &
#endif
     &                      A, SetFillVal)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars

#ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcasti, mp_gather
#endif
!
!  Imported variable declarations.
!
      logical, intent(in), optional :: SetFillVal

      integer, intent(in) :: ng, model, ncid, ncvarid, tindex, gtype
      integer, intent(in) :: LBi, UBi, LBj, UBj

      real(r8), intent(in) :: Ascl

#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:,LBj:)
# endif
      real(r8), intent(in) :: A(LBi:,LBj:)
#else
# ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in) :: A(LBi:UBi,LBj:UBj)
#endif
!
!  Local variable declarations.
!
#ifdef MASKING
      logical :: LandFill
#endif
      integer :: i, j, ic, Npts
      integer :: Imin, Imax, Jmin, Jmax
      integer :: Ilen, Jlen, IJlen, MyType, status

      integer, dimension(3) :: start, total

      integer :: nf_fwrite2d

      real(r8), dimension((Lm(ng)+2)*(Mm(ng)+2)) :: Aout
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
        CASE (p2dvar)
          Imin=IOBOUNDS(ng)%ILB_psi
          Imax=IOBOUNDS(ng)%IUB_psi
          Jmin=IOBOUNDS(ng)%JLB_psi
          Jmax=IOBOUNDS(ng)%JUB_psi
        CASE (r2dvar)
          Imin=IOBOUNDS(ng)%ILB_rho
          Imax=IOBOUNDS(ng)%IUB_rho
          Jmin=IOBOUNDS(ng)%JLB_rho
          Jmax=IOBOUNDS(ng)%JUB_rho
        CASE (u2dvar)
          Imin=IOBOUNDS(ng)%ILB_u
          Imax=IOBOUNDS(ng)%IUB_u
          Jmin=IOBOUNDS(ng)%JLB_u
          Jmax=IOBOUNDS(ng)%JUB_u
        CASE (v2dvar)
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
      IJlen=Ilen*Jlen

#ifdef MASKING
!
!  Set switch to replace land areas with fill value, spval.
!
      IF (PRESENT(SetFillVal)) THEN
        LandFill=SetFillVal
      ELSE
        LandFill=tindex.gt.0
      END IF
#endif
!
!  Initialize local array to avoid denormalized numbers. This
!  facilitates processing and debugging.
!
      Aout=0.0_r8

#ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  If distributed-memory set-up, collect tile data from all spawned
!  nodes and store it into a global scratch 1D array, packed in column-
!  major order.
# ifdef MASKING
#  ifdef WRiTE_WATER 
!  Remove land points and pack water points into 1D-array.
#  else
!  Overwrite masked points with special value.
#  endif
# endif
!-----------------------------------------------------------------------
!
      CALL mp_gather (ng, model, LBi, UBi, LBj, UBj, 1, 1,              &
     &                tindex, gtype, Ascl,                              &
# ifdef MASKING
     &                Amask,                                            &
# endif
     &                A, Npts, Aout, SetFillVal)
#else
!
!-----------------------------------------------------------------------
!  If serial or shared-memory applications and serial output, pack data
!  into a global 1D array in column-major order.
# ifdef MASKING
#  ifdef WRiTE_WATER 
!  Remove land points and pack water points into 1D-array.
#  else
!  Overwrite masked points with special value.
#  endif
# endif
!-----------------------------------------------------------------------
!
      IF (gtype.gt.0) THEN
        ic=0
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            ic=ic+1
            Aout(ic)=A(i,j)*Ascl
# ifdef MASKING
            IF ((Amask(i,j).eq.0.0_r8).and.LandFill) THEN
              Aout(ic)=spval
            END IF
# endif
          END DO
        END DO
        Npts=IJlen
# ifdef MASKING
      ELSE
        Npts=0
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            IF (Amask(i,j).gt.0.0_r8) THEN
              Npts=Npts+1
              Aout(Npts)=A(i,j)*Ascl
            END IF
          END DO
        END DO
# endif
      END IF
#endif
!
!-----------------------------------------------------------------------
!  Write output buffer into NetCDF file.
!-----------------------------------------------------------------------
!
      nf_fwrite2d=nf90_noerr
      IF (OutThread) THEN
        IF (gtype.gt.0) THEN
          start(1)=1
          total(1)=Ilen
          start(2)=1
          total(2)=Jlen
          start(3)=tindex
          total(3)=1
#ifdef MASKING
        ELSE
          start(1)=1
          total(1)=Npts
          start(2)=tindex
          total(2)=1
#endif
        END IF
        status=nf90_put_var(ncid, ncvarid, Aout, start, total)
      END IF
#ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Broadcast IO error flag to all nodes.
!-----------------------------------------------------------------------
!
      CALL mp_bcasti (ng, model, status, 1)
#endif
      nf_fwrite2d=status
      RETURN
      END FUNCTION nf_fwrite2d

      END MODULE nf_fwrite2d_mod
