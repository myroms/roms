#include "cppdefs.h"
      MODULE mod_sources
#if defined UV_PSOURCE || defined TS_PSOURCE
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  Dsrc       Direction of point Sources/Sinks:                        !
!               Dsrc(:) = 0,  Along XI-direction.                      !
!               Dsrc(:) > 0,  Along ETA-direction.                     !
!  Fsrc       Point Source/Sinks identification flag:                  !
!               Fsrc(:) = 0,  All Tracer source/sink are off.          !
!               Fsrc(:) = 1,  Only temperature is on.                  !
!               Fsrc(:) = 2,  Only salinity is on.                     !
!               Fsrc(:) = 3,  Both temperature and salinity are on.    !
!               Fsrc(:) = 4,  Both nitrate and salinity are on.        !
!  Isrc       I-grid location of point Sources/Sinks,                  !
!               1 =< Isrc =< Lm(ng).                                   !
!  Jsrc       J-grid location of point Sources/Sinks,                  !
!               1 =< Jsrc =< Mm(ng).                                   !
!  Lsrc       Logical switch indicating which tracer field to apply    !
!               the point Sources/Sinks.                               !
!  Qbar       Vertically integrated mass transport (m3/s) of point     !
!               Sources/Sinks at U- or V-points:                       !
!               Qbar -> positive, if the mass transport is in the      !
!                       positive U- or V-direction.                    !
!               Qbar -> negative, if the mass transport is in the      !
!                       negative U- or V-direction.                    !
!  QbarG      Latest two-time snapshots of vertically integrated       !
!               mass transport (m3/s) of point Sources/Sinks.          !
!  Qshape     Nondimensional shape function to distribute mass         !
!               mass point Sources/Sinks vertically.                   !
!  Qsrc       Mass transport profile (m3/s) of point Sources/Sinks.    !
!  Tsrc       Tracer (tracer units) point Sources/Sinks.               !
!  TsrcG      Latest two-time snapshots of tracer (tracer units)       !
!               point Sources/Sinks.                                   !
!                                                                      !
!=======================================================================
!
        USE mod_kinds

        implicit none

        TYPE T_SOURCES

          logical, pointer :: Lsrc(:,:)

          integer, pointer :: Isrc(:)
          integer, pointer :: Jsrc(:)

          real(r8), pointer :: Dsrc(:)
          real(r8), pointer :: Fsrc(:)
          real(r8), pointer :: Qbar(:)
          real(r8), pointer :: Qshape(:,:)
          real(r8), pointer :: Qsrc(:,:)
          real(r8), pointer :: Tsrc(:,:,:)
          real(r8), pointer :: Xsrc(:)
          real(r8), pointer :: Ysrc(:)

# ifndef ANA_PSOURCE
          real(r8), pointer :: QbarG(:,:)
          real(r8), pointer :: TsrcG(:,:,:,:)
# endif

        END TYPE T_SOURCES

        TYPE (T_SOURCES), allocatable :: SOURCES(:)

      CONTAINS

      SUBROUTINE allocate_sources (ng)
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine allocates and initializes all variables in the module  !
!  for all nested grids.                                               !
!                                                                      !
!=======================================================================
!
      USE mod_param
!
!  Local variable declarations.
!
      integer :: ng
 
      real(r8), parameter :: IniVal = 0.0_r8
!
!-----------------------------------------------------------------------
!  Allocate and initialize module variables.
!-----------------------------------------------------------------------
!
      IF (ng.eq.1) allocate ( SOURCES(Ngrids) )
!
      allocate ( SOURCES(ng) % Isrc(Msrc) )
      SOURCES(ng) % Isrc = 0

      allocate ( SOURCES(ng) % Jsrc(Msrc) )
      SOURCES(ng) % Jsrc = 0

      allocate ( SOURCES(ng) % Lsrc(Msrc,NT(ng)) )
      SOURCES(ng) % Lsrc = .false.

      allocate ( SOURCES(ng) % Dsrc(Msrc) )
      SOURCES(ng) % Dsrc = IniVal

      allocate ( SOURCES(ng) % Fsrc(Msrc) )
      SOURCES(ng) % Fsrc = IniVal

      allocate ( SOURCES(ng) % Qbar(Msrc) )
      SOURCES(ng) % Qbar = IniVal

      allocate ( SOURCES(ng) % Qshape(Msrc,N(ng)) )
      SOURCES(ng) % Qshape = IniVal

      allocate ( SOURCES(ng) % Qsrc(Msrc,N(ng)) )
      SOURCES(ng) % Qsrc = IniVal

      allocate ( SOURCES(ng) % Tsrc(Msrc,N(ng),NT(ng)) )
      SOURCES(ng) % Tsrc = IniVal

      allocate ( SOURCES(ng) % Xsrc(Msrc) )
      SOURCES(ng) % Xsrc = IniVal

      allocate ( SOURCES(ng) % Ysrc(Msrc) )
      SOURCES(ng) % Ysrc = IniVal

# ifndef ANA_PSOURCE
      allocate ( SOURCES(ng) % QbarG(Msrc,2) )
      SOURCES(ng) % QbarG = IniVal

      allocate ( SOURCES(ng) % TsrcG(Msrc,N(ng),2,NT(ng)) )
      SOURCES(ng) % TsrcG = IniVal
# endif
      RETURN
      END SUBROUTINE allocate_sources
#endif
      END MODULE mod_sources
