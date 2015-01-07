/*
** svn $Id$
***************************************************** John C. Warner ***
** Copyright (c) 2002-2015 The ROMS/TOMS Group      Hernan G. Arango  **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.txt                                             **
************************************************************************
**                                                                    **
** These routines are use couple ROMS/TOMS to WRF atmosphere model    **
** using the Model Coupling Toolkit (MCT).                            **
**                                                                    **
************************************************************************
*/

      SUBROUTINE initialize_ocn2atm_coupling (ng, tile)
!
!=======================================================================
!                                                                      !
!  Initialize ocean and atmosphere models coupling stream. This is     !
!  the training phase used to constuct  MCT parallel interpolators     !
!  and stablish communication patterns.                                !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_coupler
      USE mod_forces
      USE mod_kinds
      USE mod_scalars
!
!  Imported variable definitions.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
      integer :: Istr, Iend, Jstr, Jend
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: Asize, Jsize, MyError
      integer :: j, jc, nprocs

      integer, allocatable :: length(:)
      integer, allocatable :: start(:)
!
!-----------------------------------------------------------------------
!  Compute lower and upper bounds over a particular domain partition or
!  tile for RHO-, U-, and V-variables. Notice that "set_bounds.h" is
!  not used here because of implementation of periodicity in other
!  models.
!-----------------------------------------------------------------------
!
      Istr=BOUNDS(ng)%Istr(tile)
      Iend=BOUNDS(ng)%Iend(tile)
      Jstr=BOUNDS(ng)%Jstr(tile)
      Jend=BOUNDS(ng)%Jend(tile)
!
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
        IstrR=BOUNDS(ng)%Istr(tile)-1
      ELSE
        IstrR=BOUNDS(ng)%Istr(tile)
      END IF
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
        IendR=BOUNDS(ng)%Iend(tile)+1
      ELSE
        IendR=BOUNDS(ng)%Iend(tile)
      END IF
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
        JstrR=BOUNDS(ng)%Jstr(tile)-1
      ELSE
        JstrR=BOUNDS(ng)%Jstr(tile)
      END IF
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
        JendR=BOUNDS(ng)%Jend(tile)+1
      ELSE
        JendR=BOUNDS(ng)%Jend(tile)
      END IF
!
!-----------------------------------------------------------------------
!  Begin initialization phase.
!-----------------------------------------------------------------------
!
!  Get communicator local rank and size.
!
      CALL mpi_comm_rank (OCN_COMM_WORLD, MyRank, MyError)
      CALL mpi_comm_size (OCN_COMM_WORLD, nprocs, MyError)

#ifndef WAVES_OCEAN
!
!  Initialize MCT coupled model registry.
!
      CALL MCTWorld_init (Nmodels, MPI_COMM_WORLD, OCN_COMM_WORLD,      &
     &                    OCNid)
!
!  Determine start and lengths for domain decomposition.
!
      Jsize=JendR-JstrR+1
      IF (.not.allocated(start)) THEN
        allocate ( start(Jsize) )
      END IF
      IF (.not.allocated(length)) THEN
        allocate ( length(Jsize) )
      END IF
      jc=0
      DO j=JstrR,JendR
        jc=jc+1
        start (jc)=j*(Lm(ng)+2)+IstrR+1
        length(jc)=(IendR-IstrR+1)
      END DO
      CALL GlobalSegMap_init (GSMapROMS, start, length, 0,              &
     &                        OCN_COMM_WORLD, OCNid)
!
!  Initialize attribute vector holding the export data code strings of
!  the atmosphere model. The Asize is the number of grid point on this
!  processor.
!
      Asize=GlobalSegMap_lsize(GSMapROMS, OCN_COMM_WORLD)
      CALL AttrVect_init (atm2ocn_AV, rList=TRIM(ExportList(Iatmos)),   &
     &                    lsize=Asize)
!
!  Initialize attribute vector holding the export data code string of
!  the ocean model.
!
      CALL AttrVect_init (ocn2atm_AV, rList=TRIM(ExportList(Iocean)),   &
     &                    lsize=Asize)
      CALL AttrVect_zero (ocn2atm_AV)
!
!  Initialize a router to the atmosphere model component.
!
      CALL Router_init (ATMid, GSMapROMS, OCN_COMM_WORLD, ROMStoWRF)
!
!  Deallocate working arrays.
!
      IF (allocated(start)) THEN
        deallocate (start)
      END IF
      IF (allocated(length)) THEN
        deallocate (length)
      END IF

      RETURN
      END SUBROUTINE initialize_ocn2atm_coupling

      SUBROUTINE ocn2atm_coupling (ng, tile)
!
!=======================================================================
!                                                                      !
!  This subroutine acquires the coupling data streams between ocean    !
!  and atmosphere models. Currently, the following data streams are    !
!  coded:                                                              !
!                                                                      !
!     (...) WRF  units                                                 !
!     [...] ROMS units                                                 !
!                                                                      !
!  Fields imported WRF model:                                          !
!                                                                      !
!     * Surface atmospheric pressure (Pa), [mb]                        !
!     * Surface air relative humidity (percent), [fraction]            !
!     * Surface (2 m) air temperature (Celsius), [Celsius]             !
!     * Surface (10 m) U-wind speed (m/s), [m/s]                       !
!     * Surface (10 m) V-wind speed (m/s), [m/s]                       !
!     * Cloud fraction (percent/100), [percent/100]                    !
!     * Precipitation (m/s), [kg/m2/s]                                 !
!     * Shortwave radiation (Watts/m2), [Celsius m/s]                  !
!     * Long wave raditaion (Watts/m2), [Celsius m/s]                  !
!     * Latent heat flux (Watts/m2), [Celsius m/s]                     !
!     * Sensible heat flux (Watts/m2), [Celsius m/s]                   !
!     * Net surface heat flux (Watts/2), [Celsius m/s]                 !
!     * Surface U-wind stress (Pa), [m2/s2]                            !
!     * Surface V-wind stress (Pa), [m2/s2]                            !
!                                                                      !
!  Fields exported to WRF model:                                       !
!                                                                      !
!     * Sea surface potential temperature (Celsius), [Celsius]         !
!                                                                      !
!=======================================================================
!
      USE mod_param
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
#include "tile.h"
!
#ifdef PROFILE
      CALL wclock_on (ng, iNLM, 48)
#endif
      CALL ocn2atm_coupling_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj)
#ifdef PROFILE
      CALL wclock_off (ng, iNLM, 48)
#endif

      RETURN
      END SUBROUTINE ocn2atm_coupling
!
!***********************************************************************
      SUBROUTINE ocn2atm_coupling_tile (ng, tile,                       &
     &                                  LBi, UBi, LBj, UBj)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_coupler
      USE mod_forces
      USE mod_ocean
      USE mod_scalars
      USE mod_stepping
      USE mod_iounits
!
      USE distribute_mod, ONLY : mp_reduce
      USE ROMS_import_mod, ONLY : ROMS_import2d
      USE ROMS_export_mod, ONLY : ROMS_export2d
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
!  Local variable declarations.
!
      integer :: Istr, Iend, Jstr, Jend
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: Asize, Iimport, Iexport, MyError
      integer :: gtype, i, id, ifield, j, status

      real(r8) :: add_offset, scale
      real(r8) :: RecvTime, SendTime, buffer(2), wtime(2)

      real(r8) :: my_wtime

      real(r8), pointer :: A(:)

      character (len=3 ), dimension(2) :: op_handle
      character (len=40) :: code
!
!-----------------------------------------------------------------------
!  Compute lower and upper bounds over a particular domain partition or
!  tile for RHO-, U-, and V-variables. Notice that "set_bounds.h" is
!  not used here because of implementation of periodicity in other
!  models.
!-----------------------------------------------------------------------
!
      Istr=BOUNDS(ng)%Istr(tile)
      Iend=BOUNDS(ng)%Iend(tile)
      Jstr=BOUNDS(ng)%Jstr(tile)
      Jend=BOUNDS(ng)%Jend(tile)
!
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
        IstrR=BOUNDS(ng)%Istr(tile)-1
      ELSE
        IstrR=BOUNDS(ng)%Istr(tile)
      END IF
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
        IendR=BOUNDS(ng)%Iend(tile)+1
      ELSE
        IendR=BOUNDS(ng)%Iend(tile)
      END IF
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
        JstrR=BOUNDS(ng)%Jstr(tile)-1
      ELSE
        JstrR=BOUNDS(ng)%Jstr(tile)
      END IF
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
        JendR=BOUNDS(ng)%Jend(tile)+1
      ELSE
        JendR=BOUNDS(ng)%Jend(tile)
      END IF
!
!-----------------------------------------------------------------------
!  Allocate communications array.
!-----------------------------------------------------------------------
!
      Asize=GlobalSegMap_lsize (GSMapROMS, OCN_COMM_WORLD)
      allocate ( A(Asize) )
      A=0.0_r8
!
!  Initialize coupling wait time clocks.
!
      RecvTime=0.0_r8
      SendTime=0.0_r8
!
!-----------------------------------------------------------------------
!  Import fields from atmosphere model (WRF) to ocean model (ROMS).
!  Currently, both atmosphere and ocean model grids are the same.
!  We need to revisit this logic to allow interpolation.
!-----------------------------------------------------------------------
!
!  Schedule receiving fields from atmosphere model.
!
      CALL mpi_comm_rank (OCN_COMM_WORLD, MyRank, MyError)
      buffer(1)=my_wtime(wtime)
      CALL MCT_Recv (atm2ocn_AV, ROMStoWRF, MyError)
      RecvTime=RecvTime+my_wtime(wtime)-buffer(1)
      IF (MyError.ne.0) THEN
        IF (Master) THEN
          WRITE (stdout,10) 'atmosphere model, MyError = ', MyError
        END IF
        exit_flag=2
        RETURN
      END IF
!
!  Receive fields from atmosphere model.
!
      Iimport=0
      DO ifield=1,Nimport(Iocean)
        id=ImportID(Iocean)%val(ifield)
        code=ADJUSTL(Fields(id)%code)
        gtype=Fields(id)%GridType
        scale=Fields(id)%scale
        add_offset=Fields(id)%AddOffset

        SELECT CASE (TRIM(code))

#if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS

          CASE ('Pair')                   ! surface air pressure

            CALL AttrVect_exportRAttr (atm2ocn_AV, TRIM(code), A, Asize)
            Iimport=Iimport+1
            scale=0.01_r8                 ! Pa to mb
            add_offset=0.0_r8
            CALL ROMS_import2d (ng, tile,                               &
     &                          id, gtype, scale, add_offset,           &
     &                          Asize, A,                               &
     &                          IstrR, IendR, JstrR, JendR,             &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fields(id)%ImpMin, Fields(id)%ImpMax,   &
     &                          FORCES(ng)%Pair,                        &
     &                          status)
#endif
#if defined BULK_FLUXES || defined ECOSIM || \
   (defined SHORTWAVE && defined ANA_SRFLUX)

          CASE ('Hair')                   ! surface air humidity

            CALL AttrVect_exportRAttr (atm2ocn_AV, TRIM(code), A, Asize)
            Iimport=Iimport+1
            scale=0.01_r8                 ! percent to fraction
            add_offset=0.0_r8
            CALL ROMS_import2d (ng, tile,                               &
     &                          id, gtype, scale, add_offset,           &
     &                          Asize, A,                               &
     &                          IstrR, IendR, JstrR, JendR,             &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fields(id)%ImpMin, Fields(id)%ImpMax,   &
     &                          FORCES(ng)%Hair,                        &
     &                          status)

          CASE ('Tair')                   ! surface (2m) air temperature

            CALL AttrVect_exportRAttr (atm2ocn_AV, TRIM(code), A, Asize)
            Iimport=Iimport+1
            scale=1.0_r8                  ! Celsius
            add_offset=0.0_r8
            CALL ROMS_import2d (ng, tile,                               &
     &                          id, gtype, scale, add_offset,           &
     &                          Asize, A,                               &
     &                          IstrR, IendR, JstrR, JendR,             &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fields(id)%ImpMin, Fields(id)%ImpMax,   &
     &                          FORCES(ng)%Tair,                        &
     &                          status)
#endif
#if defined BULK_FLUXES || defined ECOSIM

          CASE ('Uwind')                  ! U-wind (10m) component

            CALL AttrVect_exportRAttr (atm2ocn_AV, TRIM(code), A, Asize)
            Iimport=Iimport+1
            scale=1.0_r8                  ! m/s
            add_offset=0.0_r8
            CALL ROMS_import2d (ng, tile,                               &
     &                          id, gtype, scale, add_offset,           &
     &                          Asize, A,                               &
     &                          IstrR, IendR, JstrR, JendR,             &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fields(id)%ImpMin, Fields(id)%ImpMax,   &
     &                          FORCES(ng)%Uwind,                       &
     &                          status)

          CASE ('Vwind')                  ! V-wind (10m) component

            CALL AttrVect_exportRAttr (atm2ocn_AV, TRIM(code), A, Asize)
            Iimport=Iimport+1
            scale=1.0_r8                  ! m/s
            add_offset=0.0_r8
            CALL ROMS_import2d (ng, tile,                               &
     &                          id, gtype, scale, add_offset,           &
     &                          Asize, A,                               &
     &                          IstrR, IendR, JstrR, JendR,             &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fields(id)%ImpMin, Fields(id)%ImpMax,   &
     &                          FORCES(ng)%Vwind,                       &
     &                          status)
#endif
#ifdef CLOUDS

          CASE ('cloud')                  ! cloud fraction

            CALL AttrVect_exportRAttr (atm2ocn_AV, TRIM(code), A, Asize)
            Iimport=Iimport+1
            scale=1.0_r8                  ! percent/100, so 0 to 1
            add_offset=0.0_r8
            CALL ROMS_import2d (ng, tile,                               &
     &                          id, gtype, scale, add_offset,           &
     &                          Asize, A,                               &
     &                          IstrR, IendR, JstrR, JendR,             &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fields(id)%ImpMin, Fields(id)%ImpMax,   &
     &                          FORCES(ng)%cloud,                       &
     &                          status)
#endif
#ifdef BULK_FLUXES

          CASE ('rain')                   ! precipitation

            CALL AttrVect_exportRAttr (atm2ocn_AV, TRIM(code), A, Asize)
            Iimport=Iimport+1
            scale=rho0                    ! kg/m2/s
            add_offset=0.0_r8
            CALL ROMS_import2d (ng, tile,                               &
     &                          id, gtype, scale, add_offset,           &
     &                          Asize, A,                               &
     &                          IstrR, IendR, JstrR, JendR,             &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fields(id)%ImpMin, Fields(id)%ImpMax,   &
     &                          FORCES(ng)%rain,                        &
     &                          status)

          CASE ('LWrad')                  ! longwave radiation

            CALL AttrVect_exportRAttr (atm2ocn_AV, TRIM(code), A, Asize)
            Iimport=Iimport+1
            scale=-1.0_r8/(rho0*Cp)       ! Watts/m2 to Celsius m/s
            add_offset=0.0_r8
            CALL ROMS_import2d (ng, tile,                               &
     &                          id, gtype, scale, add_offset,           &
     &                          Asize, A,                               &
     &                          IstrR, IendR, JstrR, JendR,             &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fields(id)%ImpMin, Fields(id)%ImpMax,   &
     &                          FORCES(ng)%lrflx,                       &
     &                          status)

          CASE ('Lheat')                  ! latent heat flux

            CALL AttrVect_exportRAttr (atm2ocn_AV, TRIM(code), A, Asize)
            Iimport=Iimport+1
            scale=1.0_r8                  ! Watts/m2
            add_offset=0.0_r8
            CALL ROMS_import2d (ng, tile,                               &
     &                          id, gtype, scale, add_offset,           &
     &                          Asize, A,                               &
     &                          IstrR, IendR, JstrR, JendR,             &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fields(id)%ImpMin, Fields(id)%ImpMax,   &
     &                          FORCES(ng)%lhflx,                       &
     &                          status)

          CASE ('Sheat')                  ! sensible heat flux

            CALL AttrVect_exportRAttr (atm2ocn_AV, TRIM(code), A, Asize)
            Iimport=Iimport+1
            scale=1.0_r8                  ! Watts/m2
            add_offset=0.0_r8
            CALL ROMS_import2d (ng, tile,                               &
     &                          id, gtype, scale, add_offset,           &
     &                          Asize, A,                               &
     &                          IstrR, IendR, JstrR, JendR,             &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fields(id)%ImpMin, Fields(id)%ImpMax,   &
     &                          FORCES(ng)%shflx,                       &
     &                          status)
#endif
#ifdef SHORTWAVE

          CASE ('SWrad')                  ! shortwave radiation

            CALL AttrVect_exportRAttr (atm2ocn_AV, TRIM(code), A, Asize)
            Iimport=Iimport+1
            scale=-1.0_r8/(rho0*Cp)       ! Watts/m2 to Celsius m/s
            add_offset=0.0_r8
            CALL ROMS_import2d (ng, tile,                               &
     &                          id, gtype, scale, add_offset,           &
     &                          Asize, A,                               &
     &                          IstrR, IendR, JstrR, JendR,             &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fields(id)%ImpMin, Fields(id)%ImpMax,   &
     &                          FORCES(ng)%srflx,                       &
     &                          status)
#endif
#ifndef BULK_FLUXES

          CASE ('heat')                   ! surface net heat flux

            CALL AttrVect_exportRAttr (atm2ocn_AV, TRIM(code), A, Asize)
            Iimport=Iimport+1
            scale=-1.0_r8/(rho0*Cp)       ! Watts/m2 to Celsius m/s
            add_offset=0.0_r8
            CALL ROMS_import2d (ng, MyRank,                             &
     &                          id, gtype, scale, add_offset,           &
     &                          Asize, A                                &
     &                          IstrR, IendR, JstrR, JendR,             &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fields(id)%ImpMin, Fields(id)%ImpMax,   &
     &                          FORCES(ng)%stflx(:,:,itemp),            &
     &                          status)

          CASE ('Ustr')                   ! surface U-wind stress

            CALL AttrVect_exportRAttr (atm2ocn_AV, TRIM(code), A, Asize)
            Iimport=Iimport+1
            scale=1.0_r8/rho0             ! Pa to m2/s2
            add_offset=0.0_r8
            CALL ROMS_import2d (ng, tile,                               &
     &                          id, gtype, scale, add_offset,           &
     &                          Asize, A,                               &
     &                          IstrR, IendR, JstrR, JendR,             &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fields(id)%ImpMin, Fields(id)%ImpMax,   &
     &                          FORCES(ng)%sustr,                       &
     &                          status)

          CASE ('Vstr')                   ! surface V-wind stress

            CALL AttrVect_exportRAttr (atm2ocn_AV, TRIM(code), A, Asize)
            Iimport=Iimport+1
            scale=1.0_r8/rho0             ! Pa to m2/s2
            add_offset=0.0_r8
            CALL ROMS_import2d (ng, tile,                               &
     &                          id, gtype, scale, add_offset,           &
     &                          Asize, A,                               &
     &                          IstrR, IendR, JstrR, JendR,             &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fields(id)%ImpMin, Fields(id)%ImpMax,   &
     &                          FORCES(ng)%svstr,                       &
     &                          status)
#endif
        END SELECT
      END DO
!
!-----------------------------------------------------------------------
!  Export fields from ocean (ROMS) to atmosphere (WRF) model.
!-----------------------------------------------------------------------
!
!  Schedule sending fields to the atmosphere model.
!
      Iexport=0
      DO ifield=1,Nexport(Iocean)
        id=ExportID(Iocean)%val(ifield)
        code=ADJUSTL(Fields(id)%code)
        gtype=Fields(id)%GridType
        scale=Fields(id)%scale
        add_offset=Fields(id)%AddOffset

        SELECT CASE (TRIM(code))

          CASE ('SST')                    ! sea surface temperature

            CALL ROMS_export2d (ng, tile,                               &
     &                          id, gtype, scale, add_offset,           &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          OCEAN(ng)%t(:,:,N(ng),NOUT,itemp),      &
     &                          Fields(id)%ExpMin, Fields(id)%ExpMax,   &
     &                          Asize, A,                               &
     &                          status)
            CALL AttrVect_importRAttr (ocn2atm_AV, TRIM(code), A, Asize)
            Iexport=Iexport+1

        END SELECT
      END DO
!
!  Send ocean fields to atmosphere model.
!
      IF (Iexport.gt.0) THEN
        buffer(2)=my_wtime(wtime)
        CALL MCT_Send (ocn2atm_AV, ROMStoWRF, MyError)
        SendTime=SendTime+my_wtime(wtime)-buffer(2)
        IF (MyError.ne.0) THEN
          IF (Master) THEN
            WRITE (stdout,20) 'atmosphere model, MyError = ', MyError
          END IF
          exit_flag=2
          RETURN
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Report.
!-----------------------------------------------------------------------
!
      IF (Nthreads(Iocean).gt.1) THEN
        buffer(1)=RecvTime
        buffer(2)=SendTime
        op_handle(1)='SUM'
        op_handle(2)='SUM'
        CALL mp_reduce (ng, iNLM, 2, buffer, op_handle)
        RecvTime=buffer(1)
        SendTime=buffer(2)
      END IF
      IF (Master.and.((Iimport.gt.0).or.(Iexport.gt.0))) THEN
        WRITE (stdout,30) Iimport, Iexport, time_code(ng),              &
     &                    RecvTime, SendTime
        IF (Lreport) THEN
          DO ifield=1,Nimport(Iocean)
            id=ImportID(Iocean)%val(ifield)
            WRITE (stdout,40) 'ROMS Import: ',TRIM(fields(id)%name),    &
     &                        Fields(id)%ImpMin, Fields(id)%ImpMax
          END DO
          DO ifield=1,Nexport(Iocean)
            id=ExportID(Iocean)%val(ifield)
            WRITE (stdout,40) 'ROMS Export: ',TRIM(fields(id)%name),    &
     &                        Fields(id)%ExpMin, Fields(id)%ExpMax
          END DO
        END IF
      END IF
!
!  Deallocate communication arrays.
!
      deallocate (A)
!
 10   FORMAT (' OCN2ATM_COUPLING - error while receiving fields from ', &
     &        a, i4)
 20   FORMAT (' OCN2ATM_COUPLING - error while sending fields to ',     &
     &        a, i4)
 30   FORMAT (6x,'OCN2ATM   - (', i2.2, ') imported and (', i2.2,       &
     &        ') exported fields,', t62, 't = ', a,/, 16x,              &
     &        '- ROMS coupling exchages wait clock (s):',/, 19x,        &
     &        '(Recv= ', 1p,e14.8,0p, ' Send= ', 1p,e14.8,0p,')')
 40   FORMAT (16x,'- ',a,a,                                             &
     &        /,19x,'(Min= ',1p,e15.8,0p,' Max= ',1p,e15.8,0p,')')

      RETURN
      END SUBROUTINE ocn2atm_coupling_tile

      SUBROUTINE finalize_ocn2atm_coupling
!
!========================================================================
!                                                                       !
!  This routine finalizes ocean and atmosphere models coupling data     !
!  streams.                                                             !
!                                                                       !
!========================================================================
!
!  Local variable declarations.
!
      integer :: MyError
!
!-----------------------------------------------------------------------
!  Deallocate MCT environment.
!-----------------------------------------------------------------------
!
      CALL Router_clean (ROMStoWRF, MyError)
      CALL AttrVect_clean (ocn2atm_AV, MyError)
      CALL GlobalSegMap_clean (GSMapROMS, MyError)

      RETURN

      END SUBROUTINE finalize_ocn2atm_coupling
