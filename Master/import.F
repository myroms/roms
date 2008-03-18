#include "cppdefs.h"
      MODULE import_mod
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module contains several routines to interpolate imported       !
!  to particular model grid during coupling.                           !   
!                                                                      !
!=======================================================================
!
      USE mod_kinds

      implicit none

      PUBLIC :: ROMS_import2d

      CONTAINS
!
!***********************************************************************
      SUBROUTINE ROMS_import2d (ng, id, gtype, scale, add_offset,       &
     &                          InpField,                               &
     &                          Imin, Imax, Jmin, Jmax,                 &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          OutField,                               &
     &                          status)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, id, gtype
      integer, intent(in) :: Imin, Imax, Jmin, Jmax
      integer, intent(in) :: LBi, UBi, LBj, UBj

      integer, intent(out) :: status

      real(r8), intent(in) :: scale, add_offset

#ifdef ASSUMED_SHAPE
      real(r8), intent(in) ::  InpField(:,:)

      real(r8), intent(out) :: OutField(LBi:,LBj:)
#else
      real(r8), intent(in) ::  InpField(:,:)

      real(r8), intent(out) :: OutField(LBi:UBi,LBj:UBj)
#endif
!
!  Local variable declarations.
!
      integer :: i, j
!
!-----------------------------------------------------------------------
!  Interpolate imported 2D field to ROMS grid.
!-----------------------------------------------------------------------
!
      END SUBROUTINE ROMS_import2d

      END MODULE import_mod
