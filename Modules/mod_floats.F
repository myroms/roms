#include "cppdefs.h"
      MODULE mod_floats
#ifdef FLOATS
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  Findex     Indices of spherical coordinates entries in initial      !
!               location arrays, if any.                               ! 
!  Flon       Initial longitude locations, if any.                     !
!  Flat       Initial latitude locations, if any.                      !
!  Tinfo      Float trajectory initial information.                    !
!  bounded    Float bounded status switch.                             !
!  track      Multivariate float trajectory data at several time       !
!               time levels.                                           !
!                                                                      !
!=======================================================================
!
        USE mod_param
!
        implicit none

        TYPE T_FLT

          logical, pointer  :: bounded(:)

          integer, pointer :: Findex(:)

          real(r8), pointer :: Flon(:)
          real(r8), pointer :: Flat(:)
          real(r8), pointer :: Tinfo(:,:)
          real(r8), pointer :: track(:,:,:)

        END TYPE T_FLT

        TYPE (T_FLT), allocatable :: FLT(:)

      CONTAINS

      SUBROUTINE allocate_floats (ng)
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine allocates and initialize all variables in the module   !
!  for all nested grids.                                               !
!                                                                      !
!=======================================================================
!
      USE mod_scalars
!
!  Local variable declarations.
!
      integer, intent(in) :: ng
      real(r8), parameter :: IniVal = 0.0_r8
!
!-----------------------------------------------------------------------
!  Allocate and initialize module variables.
!-----------------------------------------------------------------------
!
      IF (ng.eq.1) allocate ( FLT(Ngrids) )
!
      allocate ( FLT(ng) % bounded(Nfloats(ng)) )
      FLT(ng) % bounded = .false.

      allocate ( FLT(ng) % Findex(0:Nfloats(ng)) )
      FLT(ng) % Findex = 0

      allocate ( FLT(ng) % Flon(Nfloats(ng)) )
      FLT(ng) % Flon = IniVal

      allocate ( FLT(ng) % Flat(Nfloats(ng)) )
      FLT(ng) % Flat = IniVal

      allocate ( FLT(ng) % Tinfo(0:izrhs,Nfloats(ng)) )
      FLT(ng) % Tinfo = IniVal

      allocate ( FLT(ng) % track(NFV(ng),0:NFT,Nfloats(ng)) )
      FLT(ng) % track = IniVal

      RETURN
      END SUBROUTINE allocate_floats
#endif
      END MODULE mod_floats
