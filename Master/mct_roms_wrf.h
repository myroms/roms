#include "cppdefs.h"
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2007 The ROMS/TOMS Group       Daniel Schaffer   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module is for coupling ROMS to WRF using the Model Coupling    !
!  Toolkit (MCT; developed at the Argonne National Laboratory)  and    !
!  the WRF I/O API.                                                    !
!                                                                      !
!  Dan Schaffer, FSL, Daniel.S.Schaffer@noaa.gov                       !
!                                                                      !
!=======================================================================
!
      implicit none

      integer, save :: ATM_TO_OCN_T_HANDLE
      integer, save :: OCN_TO_ATM_T_HANDLE

      CONTAINS

      SUBROUTINE initialize_atmos_coupling (ng, tile)
!
!=======================================================================
!                                                                      !
!  Initialize atmosphere and ocean coupling stream.  This is the       !
!  training phase use to constuct MCT parallel interpolators and       !
!  stablish communication patterns.                                    !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_forces
      USE mod_kinds
!
      include 'wrf_io_flags.h'
!
!  Imported variable definitions.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.  Currently, WRF I/O API supports
!  single precision only.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: MyStatus, RealType

      integer :: ng = 1             ! Assume coupling the top-level nest
      integer :: DomainDesc = 0

      integer :: DomainStr(2), DomainEnd(2)
      integer :: MemoryStr(2), MemoryEnd(2)
      integer :: PatchStr(2), PatchEnd(2)

      real(r4) :: FIELD(1,1)

      character (len=80) :: DateStr
      character (len=80) :: MemoryOrder
      character (len=80) :: Stagger
      character (len=80) :: DimNames
      character (len=80) :: VarName

#include "tile.h"
#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set Model Coupling Toolkit (MCT) parameters.
!-----------------------------------------------------------------------
!
      Stagger=""
      DimNames=""
      DateStr=""
      MemoryOrder="XY"

      DomainStr(1)=1
      DomainEnd(1)=Lm(ng)
      DomainStr(2)=1
      DomainEnd(2)=Mm(ng)

      MemoryStr(1)=LBi
      MemoryEnd(1)=UBi
      MemoryStr(2)=LBj
      MemoryEnd(2)=UBj

      PatchStr(1)=Istr
      PatchEnd(1)=Iend
      PatchStr(2)=Jstr
      PatchEnd(2)=Jend

      RealType=WRF_REAL          ! single precision
!!    RealType=WRF_DOUBLE        ! double precision, not yet supported
!
!-----------------------------------------------------------------------
!  Begin training phase: open a coupling stream in which the calling
!  read data from component model.
!-----------------------------------------------------------------------
!
!  Atmosphere to ocean model.
!
      CALL ext_mct_open_for_read_begin ("wrf",                          &
     &                                  MPI_COMM_WORLD,                 &
     &                                  OCN_COMM_WORLD,                 &
     &      "SPARSE_MATRIX_BASE_NAME=wrf_to_roms, COMPONENT_NAME=roms", &
     &                                  ATM_TO_OCN_T_HANDLE,            &
     &                                  MyStatus)
      IF (MyStatus.ne.0) THEN
        CALL finalize_atmos_coupling                                    &
     &                      ("Coupling stream open for read failed")
      END IF
!
!  Ocean to atmosphere model.
!
      CALL ext_mct_open_for_read_begin ("wrf",                          &
     &                                  MPI_COMM_WORLD,                 &
     &                                  OCN_COMM_WORLD,                 &
     &      "SPARSE_MATRIX_BASE_NAME=roms_to_wrf, COMPONENT_NAME=roms", &
     &                                  OCN_TO_ATM_T_HANDLE,            &
     &                                  MyStatus)
      IF (MyStatus.ne.0) THEN
        CALL finalize_atmos_coupling                                    &
     &                      ("Coupling stream open for read failed")
      END IF
!
!-----------------------------------------------------------------------
!  In training phase: construct and cache away coupling-data transfer
!  communication patterns.
!-----------------------------------------------------------------------
!
!  Surface wind stress in the XI-direction.
!
      WRITE (VarName, fmt='(a7)') "USTRESS"
      CALL ext_mct_read_field (ATM_TO_OCN_T_HANDLE,                     &
     &                         DateStr,                                 &
     &                         trim(VarName),                           &
     &                         FIELD,                                   &
     &                         RealType,                                &
     &                         MPI_COMM_WORLD,                          &
     &                         OCN_COMM_WORLD,                          &
     &                         DomainDesc,                              &
     &                         MemoryOrder,                             &
     &                         Stagger,                                 &
     &                         DimNames,                                &
     &                         DomainStr, DomainEnd,                    &
     &                         MemoryStr, MemoryEnd,                    &
     &                         PatchStr, PatchEnd,                      &
     &                         MyStatus)
      IF (MyStatus.ne.0) THEN
        CALL finalize_atmos_coupling ("Coupling training read failed")
      END IF
!
!  Surface wind stress in the ETA-direction.
!
      WRITE (VarName, fmt='(a7)') "VSTRESS"
      CALL ext_mct_read_field (ATM_TO_OCN_T_HANDLE,                     &
     &                         DateStr,                                 &
     &                         TRIM(VarName),                           &
     &                         FIELD,                                   &
     &                         RealType,                                &
     &                         MPI_COMM_WORLD,                          &
     &                         OCN_COMM_WORLD,                          &
     &                         DomainDesc,                              &
     &                         MemoryOrder,                             &
     &                         Stagger,                                 &
     &                         DimNames,                                &
     &                         DomainStr, DomainEnd,                    &
     &                         MemoryStr, MemoryEnd,                    &
     &                         PatchStr, PatchEnd,                      &
     &                         MyStatus)
      IF (MyStatus.ne.0) THEN
        CALL finalize_atmos_coupling ("Coupling training read failed")
      END IF
!
!  Sea surface temperature.
!
      WRITE (VarName, fmt='(a3)') "SST"
      CALL ext_mct_write_field (OCN_TO_ATM_T_HANDLE,                    &
     &                          DateStr,                                &
     &                          TRIM(VarName),                          &
     &                          FIELD,                                  &
     &                          RealType,                               &
     &                          MPI_COMM_WORLD,                         &
     &                          OCN_COMM_WORLD,                         &
     &                          DomainDesc,                             &
     &                          MemoryOrder,                            &
     &                          Stagger,                                &
     &                          DimNames,                               &
     &                          DomainStr, DomainEnd,                   &
     &                          MemoryStr, MemoryEnd,                   &
     &                          PatchStr, PatchEnd,                     &
     &                          MyStatus)
      IF (MyStatus.ne.0) THEN
        CALL finalize_atmos_coupling ("Coupling training write failed")
      END IF
!
!-----------------------------------------------------------------------
!  End of training phase: the coupling stream is referred to by the
!  data handle. 
!-----------------------------------------------------------------------
!
!  Read: atmosphere to ocean model.
!
      CALL ext_mct_open_for_read_commit (ATM_TO_OCN_T_HANDLE,           &
     &                                   MyStatus)
      IF (MyStatus.ne.0) THEN
        CALL finalize_atmos_coupling ("Coupling read commit failed")
      END IF
!
!  Write: ocean to atmosphere model.
!
      CALL ext_mct_open_for_write_commit (OCN_TO_ATM_T_HANDLE,          &
     &                                    MyStatus)
      IF (MyStatus.ne.0) THEN
        CALL finalize_atmos_coupling ("Coupling write commit failed")
      END IF

      END SUBROUTINE initialize_atmos_coupling

      SUBROUTINE atmos_coupling (ng, tile)
!
!=======================================================================
!                                                                      !
!  This subroutine reads and writes the coupling data streams between  !
!  atmosphere and ocean models. Currently, the following data streams  !
!  are processed:                                                      !
!                                                                      !
!     * Kinematic surface wind stress components (m2/s2).              !
!     * Sea surface temperature (Celsius).                             !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_forces
      USE mod_grid
      USE mod_ocean
      USE mod_stepping
!
      implicit none

      integer, intent(in) :: ng, tile

# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 36)
# endif
      CALL atmos_coupling_tile (ng, Istr, Iend, Jstr, Jend,             &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          nrhs(ng),                               &
     &                          GRID(ng) % angler,                      &
     &                          OCEAN(ng) % t,                          &
     &                          FORCES(ng) % sustr,                     &
     &                          FORCES(ng) % svstr)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 36)
# endif
      RETURN
      END SUBROUTINE atmos_coupling
!
!***********************************************************************
      SUBROUTINE atmos_coupling_tile (ng, Istr, Iend, Jstr, Jend,       &
     &                                LBi, UBi, LBj, UBj,               &
     &                                nrhs,                             &
     &                                angler,                           &
     &                                t, sustr, svstr)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_parallel
!
#if defined EW_PERIODIC || defined NS_PERIODIC
      USE exchange_2d_mod, ONLY : exchange_u2d_tile, exchange_v2d_tile
#endif
#ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
#endif
!
      implicit none
!
      include 'wrf_io_flags.h'
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: nrhs
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: angler(LBi:,LBj:)
      real(r8), intent(in) :: t(LBi:,LBj:,:,:,:)

      real(r8), intent(out) :: sustr(LBi:,LBj:)
      real(r8), intent(out) :: svstr(LBi:,LBj:)
# else
      real(r8), intent(in) :: angler(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))

      real(r8), intent(out) :: sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: svstr(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.  Currently, WRF I/O API supports
!  single precision only.
!
#ifdef DISTRIBUTE
# ifdef EW_PERIODIC
      logical :: EWperiodic=.TRUE.
# else
      logical :: EWperiodic=.FALSE.
# endif
# ifdef NS_PERIODIC
      logical :: NSperiodic=.TRUE.
# else
      logical :: NSperiodic=.FALSE.
# endif
#endif
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: MyStatus, RealType, i, j

      integer :: DomainDesc = 0
      integer :: DomainStr(2), DomainEnd(2)
      integer :: MemoryStr(2), MemoryEnd(2)
      integer :: PatchStr(2), PatchEnd(2)

      real(r4) :: var_sum, var_sum_g
      real(r8) :: cff

      real(r4), dimension(LBi:UBi,LBj:UBj) :: tmpX
      real(r4), dimension(LBi:UBi,LBj:UBj) :: tmpY

      character (len=80) :: DateStr
      character (len=80) :: MemoryOrder
      character (len=80) :: Stagger
      character (len=80) :: DimNames
      character (len=80) :: VarName

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set Model Coupling Toolkit (MCT) parameters.
!-----------------------------------------------------------------------
!
      Stagger=""
      DimNames=""
      DateStr=""
      MemoryOrder="XY"

      DomainStr(1)=1
      DomainEnd(1)=Lm(ng)
      DomainStr(2)=1
      DomainEnd(2)=Mm(ng)

      MemoryStr(1)=LBi
      MemoryEnd(1)=UBi
      MemoryStr(2)=LBj
      MemoryEnd(2)=UBj

      PatchStr(1)=Istr
      PatchEnd(1)=Iend
      PatchStr(2)=Jstr
      PatchEnd(2)=Jend

      RealType=WRF_REAL          ! single precision
!!    RealType=WRF_DOUBLE        ! double precision, not yet supported
!
!-----------------------------------------------------------------------
!  Receive kinematic wind stress (m2/s2) from atmosphere model.
!-----------------------------------------------------------------------
!
      IF (Master) THEN
        PRINT *, 'ROMS/TOMS, receiving wind stresses'
      END IF
!
!  Get wind-stress component in the XI-direction.
!
      WRITE (VarName, fmt='(a7)') "USTRESS"
      CALL ext_mct_read_field (ATM_TO_OCN_T_HANDLE,                     &
     &                         DateStr,                                 &
     &                         TRIM(VarName),                           &
     &                         tmpX,                                    &
     &                         RealType,                                &
     &                         MPI_COMM_WORLD,                          &
     &                         OCN_COMM_WORLD,                          &
     &                         DomainDesc,                              &
     &                         MemoryOrder,                             &
     &                         Stagger,                                 &
     &                         DimNames,                                &
     &                         DomainStr, DomainEnd,                    &
     &                         MemoryStr, MemoryEnd,                    &
     &                         PatchStr, PatchEnd,                      &
     &                         MyStatus)
      IF (MyStatus.ne.0) THEN
        CALL finalize_atmos_coupling ("Coupling read failed")
      END IF
!
!  Compute global average.
!
      var_sum=SUM(tmpX(PatchStr(1):PatchEnd(1),                         &
     &                 PatchStr(2):PatchEnd(2)))
      CALL mpi_allreduce (var_sum, var_sum_g, 1, MPI_REAL, MPI_SUM,     &
     &                    OCN_COMM_WORLD, MyStatus)
      IF (MyStatus.ne.MPI_SUCCESS) THEN
        CALL finalize_atmos_coupling ("Coupling global sum failed")
      END IF
      IF (Master) THEN
        cff=(DomainEnd(1)-DomainStr(1)+1)*                              &
     &      (DomainEnd(2)-DomainStr(2)+1)
        PRINT *, 'received U wind-stress from WRF, average: ',          &
     &           var_sum_g/cff
      END IF
!
!  Get wind-stress component in the ETA-direction.
!
      WRITE (VarName, FMT='(a7)') "VSTRESS"
      CALL ext_mct_read_field (ATM_TO_OCN_T_HANDLE,                     &
     &                         DateStr,                                 &
     &                         TRIM(VarName),                           &
     &                         tmpY,                                    &
     &                         RealType,                                &
     &                         MPI_COMM_WORLD,                          &
     &                         OCN_COMM_WORLD,                          &
     &                         DomainDesc,                              &
     &                         MemoryOrder,                             &
     &                         Stagger,                                 &
     &                         DimNames,                                &
     &                         DomainStr, DomainEnd,                    &
     &                         MemoryStr, MemoryEnd,                    &
     &                         PatchStr, PatchEnd,                      &
     &                         MyStatus)
      IF (MyStatus.ne.0) THEN
        CALL finalize_atmos_coupling ("Coupling read failed")
      END IF
!
!  Compute global average.
!
      var_sum=SUM(tmpY(PatchStr(1):PatchEnd(1),                         &
     &                 PatchStr(2):PatchEnd(2)))
      CALL mpi_allreduce (var_sum, var_sum_g, 1, MPI_REAL, MPI_SUM,     &
     &                    OCN_COMM_WORLD, MyStatus)
      IF (MyStatus.ne.MPI_SUCCESS) THEN
        CALL finalize_atmos_coupling ("Coupling global sum failed")
      END IF
      IF (Master) THEN
        cff=(DomainEnd(1)-DomainStr(1)+1)*                              &
     &      (DomainEnd(2)-DomainStr(2)+1)
        PRINT *, 'Received V wind-stress from WRF, average: ',          &
     &           var_sum_g/cff
      END IF
!
!  Scale to kinematic stress.
!
      cff=1.0_r8/rho0
      DO j=Jstr,Jend
        DO i=Istr,Iend
          sustr(i,j)=cff*tmpX(i,j)
          svstr(i,j)=cff*tmpY(i,j)
        END DO
      END DO

#if defined EW_PERIODIC || defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Apply periodic boundary conditions.
!-----------------------------------------------------------------------
!
      CALL exchange_u2d_tile (ng, iNLM, Istr, Iend, Jstr, Jend,         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        sustr)
      CALL exchange_v2d_tile (ng, iNLM, Istr, Iend, Jstr, Jend,         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        svstr)
#endif
#ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Exchange tile boundaries.
!-----------------------------------------------------------------------
!
      CALL mp_exchange2d (ng, iNLM, 2, Istr, Iend, Jstr, Jend,          &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    sustr, svstr)
#endif
!
!-----------------------------------------------------------------------
!  Send sea surface temperature to atmosphere model.
!-----------------------------------------------------------------------
!
!  Load sea surface temperature into temporary array.  
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          tmpX(i,j)=t(i,j,N(ng),nrhs,itemp)
        END DO
      END DO
!
!  Compute global average.
!
      var_sum=SUM(tmpX(PatchStr(1):PatchEnd(1),                         &
                       PatchStr(2):PatchEnd(2)))
      CALL mpi_allreduce (var_sum, var_sum_g, 1, MPI_REAL, MPI_SUM,     &
     &                    OCN_COMM_WORLD, MyStatus)
      IF (MyStatus.ne.MPI_SUCCESS) THEN
        CALL finalize_atmos_coupling ("Coupling global sum failed")
      END IF
!
!  Send sea surface temperature to atmospheric model.
!
      IF (Master) THEN
        cff=(DomainEnd(1)-DomainStr(1)+1)*                              &
     &      (DomainEnd(2)-DomainStr(2)+1) 
        PRINT *, 'Sending SST to WRF, average: ', var_sum_g/cff
      END IF
      WRITE (VarName, fmt='(a3)') "SST"
      CALL ext_mct_write_field (OCN_TO_ATM_T_HANDLE,                    &
     &                          DateStr,                                &
     &                          TRIM(VarName),                          &
     &                          tmpX,                                   &
     &                          RealType,                               &
     &                          MPI_COMM_WORLD,                         &
     &                          OCN_COMM_WORLD,                         &
     &                          DomainDesc,                             &
     &                          MemoryOrder,                            &
     &                          Stagger,                                &
     &                          DimNames,                               &
     &                          DomainStr, DomainEnd,                   &
     &                          MemoryStr, MemoryEnd,                   &
     &                          PatchStr, PatchEnd,                     &
     &                          MyStatus)
      IF (MyStatus.ne.0) THEN
        CALL finalize_atmos_coupling ("Coupling write failed")
      END IF

      RETURN
      END SUBROUTINE atmos_coupling_tile

      SUBROUTINE finalize_atmos_coupling (string)
!
!=======================================================================
!                                                                    ===
!  This routines terminates execution during coupling error.         ===
!                                                                    ===
!=======================================================================
!
!  Imported variable declarations.
!
      character (len=*), intent(in) :: string
!
!  Local variable declarations.
!
      integer :: MyStatus
!
!-----------------------------------------------------------------------
!  Terminate MPI execution environment.
!-----------------------------------------------------------------------
!
      PRINT *, string
      CALL mpi_finalize (MyStatus)

      STOP
      END SUBROUTINE finalize_atmos_coupling
#endif

      END MODULE atm_coupler_mod
