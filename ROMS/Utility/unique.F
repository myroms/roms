#include "cppdefs.h"     
      MODULE unique_mod
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2025 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  This MODULE contains several routines to compute the unique values  !
!  of a vector. It has similar capabilities to the Matlab "unique"     !
!  function, such that:                                                !
!                                                                      !
!    Aout = Ainp(Iout)                                                 !
!                                                                      !
!  where "Aout" is a vector of the unique values and "Iout" is an      !
!  index  vector of first occurence of each repeated value.            !
!                                                                      !
!  It is assumes that the input vector "Ainp" is sorted in ascending   !
!  values. It is primarily used for processing the time coordinate of  !
!  data assimilation observations.                                     !
!                                                                      !
!=======================================================================
!
      USE mod_kinds
!
      implicit none
!
      INTERFACE unique       
        MODULE PROCEDURE unique_1di8
        MODULE PROCEDURE unique_1dr8
      END INTERFACE unique
!
      PUBLIC  :: unique
      PRIVATE
!
      CONTAINS
!
!-----------------------------------------------------------------------
!
!  It Computes unique values for 32-bit integer vector.
!
      SUBROUTINE unique_1di8 (Ninp, Ainp, Nout, Aout, Iout)
!
      integer,                    intent(in ) :: Ninp
      integer (i8b),              intent(in ) :: Ainp(Ninp)
      integer,                    intent(out) :: Nout
      integer,       allocatable, intent(out) :: Iout(:)
      integer (i8b), allocatable, intent(out) :: Aout(:)
!
      logical                                 :: mask(Ninp)
      integer                                 :: i
!
!  Determine unique values of 1D vector Ainp.
!
      mask=.TRUE.
      DO i=Ninp,2,-1
        mask(i)=.not.(Ainp(i-1).eq.Ainp(i))
      END DO
!
!  Allocate an index 
!
      Nout=SIZE(PACK([(i,i=1,Ninp)], mask))
      allocate ( Iout(Nout) )
      Iout=PACK([(i, i=1,Ninp)], mask)
!
!  Copy the unique values.
!
      Aout=Ainp(Iout)
!
      RETURN
      END SUBROUTINE unique_1di8
!
!-----------------------------------------------------------------------
!
!  It Computes unique values for double precision vector.
!
      SUBROUTINE unique_1dr8 (Ninp, Ainp, Nout, Aout, Iout)
!
      integer,                   intent(in ) :: Ninp
      real (r8),                 intent(in ) :: Ainp(Ninp)
      integer,                   intent(out) :: Nout
      integer,      allocatable, intent(out) :: Iout(:)
      real (r8),    allocatable, intent(out) :: Aout(:)
!
      logical                                :: mask(Ninp)
      integer                                :: i
!
!  Determine unique values of 1D vector Ainp.
!
      mask=.TRUE.
      DO i=Ninp,2,-1
        mask(i)=.not.(Ainp(i-1).eq.Ainp(i))
      END DO
!
!  Allocate an index vector of first occurence of each repeated value.
!
      Nout=SIZE(PACK([(i, i=1,Ninp)], mask))
      allocate ( Iout(Nout) )
      Iout=PACK([(i, i=1,Ninp)], mask)
!
!  Copy the unique values.
!
      Aout=Ainp(Iout)
!
      RETURN
      END SUBROUTINE unique_1dr8
!
      END MODULE unique_mod
