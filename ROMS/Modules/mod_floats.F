#include "cppdefs.h"
      MODULE mod_floats
#ifdef FLOATS
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2011 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Findex     Indices of spherical coordinates entries in initial      !
!               location arrays, if any.                               !
!  Flon       Initial longitude locations, if any.                     !
!  Flat       Initial latitude locations, if any.                      !
!  Ftype      Float trajectory type:                                   !
!               Ftype(:) = 1,  neutral density 3D Lagrangian           !
!               Ftype(:) = 2,  isobaric (constant depth) float.        !
!  Tinfo      Float trajectory initial information.                    !
!  bounded    Float bounded status switch.                             !
# if defined SOLVE3D && defined FLOAT_VWALK
!  rwalk      Normally distributed random deviates used in vertical    !
!               random walk.                                           !
# endif
# if defined SOLVE3D && defined FLOAT_STICKY
!  stuck      Reflection switch. Floats that hit the surface are       !
!               reflected and floats that hitthe bottom get stick      !
# endif
!  track      Multivariate float trajectory data at several time       !
!               time levels.                                           !
!                                                                      !
!=======================================================================
!
        USE mod_param
!
        implicit none

        TYPE T_DRIFTER

          logical, pointer  :: bounded(:)

# if defined SOLVE3D && defined FLOAT_STICKY
          logical, pointer  :: stuck(:)
# endif

          integer, pointer :: Findex(:)
          integer, pointer :: Ftype(:)

          real(r8), pointer :: Flon(:)
          real(r8), pointer :: Flat(:)
          real(r8), pointer :: Fz0(:)
          real(r8), pointer :: Tinfo(:,:)
# if defined SOLVE3D && defined FLOAT_VWALK
          real(r8), pointer :: rwalk(:)
# endif
          real(r8), pointer :: track(:,:,:)

        END TYPE T_DRIFTER

        TYPE (T_DRIFTER), allocatable :: DRIFTER(:)

      CONTAINS

      SUBROUTINE allocate_floats (ng)
!
!=======================================================================
!                                                                      !
!  This routine allocates and initialize all variables in the module   !
!  for all nested grids.                                               !
!                                                                      !
!=======================================================================
!
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      integer :: i, iflt

      real(r8), parameter :: IniVal = 0.0_r8
!
!-----------------------------------------------------------------------
!  Allocate module variables.
!-----------------------------------------------------------------------
!
      IF (ng.eq.1) allocate ( DRIFTER(Ngrids) )
!
      allocate ( DRIFTER(ng) % bounded(Nfloats(ng)) )

# if defined SOLVE3D && defined FLOAT_STICKY
      allocate ( DRIFTER(ng) % stuck(Nfloats(ng)) )
# endif

      allocate ( DRIFTER(ng) % Findex(0:Nfloats(ng)) )

      allocate ( DRIFTER(ng) % Ftype(Nfloats(ng)) )

      allocate ( DRIFTER(ng) % Flon(Nfloats(ng)) )

      allocate ( DRIFTER(ng) % Flat(Nfloats(ng)) )

      allocate ( DRIFTER(ng) % Fz0(Nfloats(ng)) )

      allocate ( DRIFTER(ng) % Tinfo(0:izrhs,Nfloats(ng)) )

# if defined SOLVE3D && defined FLOAT_VWALK
      allocate ( DRIFTER(ng) % rwalk(Nfloats(ng)) )
# endif

      allocate ( DRIFTER(ng) % track(NFV(ng),0:NFT,Nfloats(ng)) )
!
!-----------------------------------------------------------------------
!  Initialize module variables.
!-----------------------------------------------------------------------
!
      DRIFTER(ng) % Findex(0) = 0
      DO iflt=1,Nfloats(ng)
        DRIFTER(ng) % bounded(iflt) = .FALSE.
# if defined SOLVE3D && defined FLOAT_STICKY
        DRIFTER(ng) % stuck(iflt) = .FALSE.
# endif
        DRIFTER(ng) % Findex(iflt) = 0
        DRIFTER(ng) % Ftype(iflt) = 0
        DRIFTER(ng) % Flon(iflt) = IniVal
        DRIFTER(ng) % Flat(iflt) = IniVal
        DRIFTER(ng) % Fz0(iflt) = 0
# if defined SOLVE3D && defined FLOAT_VWALK
        DRIFTER(ng) % rwalk = IniVal
# endif
        DO i=0,izrhs
          DRIFTER(ng) % Tinfo(i,iflt) = IniVal
        END DO
        DO i=1,NFV(ng)
          DRIFTER(ng) % track(i,0,iflt) = IniVal
          DRIFTER(ng) % track(i,1,iflt) = IniVal
          DRIFTER(ng) % track(i,2,iflt) = IniVal
          DRIFTER(ng) % track(i,3,iflt) = IniVal
          DRIFTER(ng) % track(i,4,iflt) = IniVal
        END DO
      END DO

      RETURN
      END SUBROUTINE allocate_floats
#endif
      END MODULE mod_floats
