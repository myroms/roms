#include "cppdefs.h"
      MODULE mod_sources
#if defined UV_PSOURCE || defined TS_PSOURCE || defined Q_PSOURCE
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2007 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
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
!               Fsrc(:) = ... And other combinations.                  !
!                             (We need a more robust logic here)       !
!  Isrc       I-grid location of point Sources/Sinks,                  !
!               1 =< Isrc =< Lm(ng).                                   !
!  Jsrc       J-grid location of point Sources/Sinks,                  !
!               1 =< Jsrc =< Mm(ng).                                   !
!  Lsrc       Logical switch for each source point source/sink datum   !
!               indicating which tracer(s) is (are) activated          !
!               according to the "Fsrc" indentification flag.          !
!  Ltracer    Logical switch indicating which tracer field need to be  !
!               processed for Sources/Sinks terms.                     !
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
          logical, pointer :: Ltracer(:)

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
!                                                                      !
!  This routine allocates and initializes all variables in the module  !
!  for all nested grids.                                               !
!                                                                      !
!=======================================================================
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer :: ng
!
!  Local variable declarations.
! 
      integer :: is, itrc, k

      real(r8), parameter :: IniVal = 0.0_r8
!
!-----------------------------------------------------------------------
!  Allocate module variables.
!-----------------------------------------------------------------------
!
      IF (ng.eq.1) allocate ( SOURCES(Ngrids) )
!
      allocate ( SOURCES(ng) % Isrc(Msrc) )

      allocate ( SOURCES(ng) % Jsrc(Msrc) )

      allocate ( SOURCES(ng) % Lsrc(Msrc,NT(ng)) )

      allocate ( SOURCES(ng) % Ltracer(NT(ng)) )

      allocate ( SOURCES(ng) % Dsrc(Msrc) )

      allocate ( SOURCES(ng) % Fsrc(Msrc) )

      allocate ( SOURCES(ng) % Qbar(Msrc) )

      allocate ( SOURCES(ng) % Qshape(Msrc,N(ng)) )

      allocate ( SOURCES(ng) % Qsrc(Msrc,N(ng)) )

      allocate ( SOURCES(ng) % Tsrc(Msrc,N(ng),NT(ng)) )

      allocate ( SOURCES(ng) % Xsrc(Msrc) )

      allocate ( SOURCES(ng) % Ysrc(Msrc) )

# ifndef ANA_PSOURCE
      allocate ( SOURCES(ng) % QbarG(Msrc,2) )

      allocate ( SOURCES(ng) % TsrcG(Msrc,N(ng),2,NT(ng)) )
# endif
!
!-----------------------------------------------------------------------
!  Initialize module variables.
!-----------------------------------------------------------------------
!
      DO is=1,Msrc
        SOURCES(ng) % Isrc(is) = 0
        SOURCES(ng) % Jsrc(is) = 0
        SOURCES(ng) % Dsrc(is) = IniVal
        SOURCES(ng) % Fsrc(is) = IniVal
        SOURCES(ng) % Xsrc(is) = IniVal
        SOURCES(ng) % Ysrc(is) = IniVal
        SOURCES(ng) % Qbar(is) = IniVal
# ifndef ANA_PSOURCE
        SOURCES(ng) % QbarG(is,1) = IniVal
        SOURCES(ng) % QbarG(is,2) = IniVal
# endif
      END DO
      DO k=1,N(ng)
        DO is=1,Msrc
          SOURCES(ng) % Qshape(is,k) = IniVal
          SOURCES(ng) % Qsrc(is,k) = IniVal
        END DO
      END DO 
      DO itrc=1,NT(ng)
        SOURCES(ng) % Ltracer(itrc) = .FALSE.
        DO is=1,Msrc
          SOURCES(ng) % Lsrc(is,itrc) = .FALSE.
        END DO
        DO k=1,N(ng)
          DO is=1,Msrc
            SOURCES(ng) % Tsrc(is,k,itrc) = IniVal
# ifndef ANA_PSOURCE
            SOURCES(ng) % TsrcG(is,k,1,itrc) = IniVal
            SOURCES(ng) % TsrcG(is,k,2,itrc) = IniVal
# endif
          END DO
        END DO
      END DO

      RETURN
      END SUBROUTINE allocate_sources
#endif
      END MODULE mod_sources
