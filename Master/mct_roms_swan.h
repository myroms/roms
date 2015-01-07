/*
** svn $Id$
***************************************************** John C. Warner ***
** Copyright (c) 2002-2015 The ROMS/TOMS Group      Hernan G. Arango  **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.txt                                             **
************************************************************************
**                                                                    **
** These routines are use couple ROMS/TOMS to SWAN wave model using   **
** the Model Coupling Toolkit (MCT).                                  **
**                                                                    **
************************************************************************
*/

      SUBROUTINE initialize_ocn2wav_coupling (ng, tile)
!
!=======================================================================
!                                                                      !
!  Initialize ocean and wave models coupling stream.  This is the      !
!  training phase used to constuct MCT parallel interpolators and      !
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
!  the wave model. The Asize is the number of grid point on this
!  processor.
!
      Asize=GlobalSegMap_lsize(GSMapROMS, OCN_COMM_WORLD)
      CALL AttrVect_init (wav2ocn_AV, rList=TRIM(ExportList(Iwaves)),   &
     &                    lsize=Asize)
!
!  Initialize attribute vector holding the export data code string of
!  the ocean model.
!
      CALL AttrVect_init (ocn2wav_AV, rList=TRIM(ExportList(Iocean)),   &
     &                    lsize=Asize)
      CALL AttrVect_zero (ocn2wav_AV)
!
!  Initialize a router to the wave model component.
!
      CALL Router_init (WAVid, GSMapROMS, OCN_COMM_WORLD, ROMStoSWAN)
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
      END SUBROUTINE initialize_ocn2wav_coupling

      SUBROUTINE ocn2wav_coupling (ng, tile)
!
!=======================================================================
!                                                                      !
!  This routine acquires the coupling data streams between waves       !
!  and ocean models.  Currently,  the following data streams are       !
!  coded:                                                              !
!                                                                      !
!     (...) SWAN units                                                 !
!     [...] ROMS units                                                 !
!                                                                      !
!  Fields imported from SWAN model:                                    !
!                                                                      !
!     * Wave direction (degrees), [radians]                            !
!     * Significant wave height (m), [m]                               !
!     * Average wave length (m), [m]                                   !
!     * Surface wave relative peak period (s), [s]                     !
!     * Bottom wave period (s), [s]                                    !
!     * Percent of breakig waves (nondimensional), [nondimensional]    !
!     * Wave energy dissipation (W/m2), [m3/s3]                        !
!     * Wave bottom orbital velocity (m/s), [m/s]                      !
!                                                                      !
!  Fields exported to SWAN model:                                      !
!                                                                      !
!     * Bathymetry, bottom elevation (m), [m]                          !
!     * Free-surface, water surface elevation (m), [m]                 !
!     * Depth integrated u-momentum (m/s), [m/s]                       !
!     * Depth integrated v-momentum (m/s), [m/s]                       !
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
      CALL ocn2wav_coupling_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj)
#ifdef PROFILE
      CALL wclock_off (ng, iNLM, 48)
#endif

      RETURN
      END SUBROUTINE ocn2wav_coupling
!
!***********************************************************************
      SUBROUTINE ocn2wav_coupling_tile (ng, tile,                       &
     &                                  LBi, UBi, LBj, UBj)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_coupler
      USE mod_forces
      USE mod_grid
      USE mod_ocean
      USE mod_scalars
      USE mod_stepping
      USE mod_iounits
#ifdef BBL_MODEL
      USE mod_sedbed
#endif
      USE mod_sediment
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
      integer :: gtype, i, id, ifield, ij, j, status

      real(r8), parameter ::  Lwave_max=500.0_r8

      real(r8) :: add_offset, scale
      real(r8) :: RecvTime, SendTime, buffer(2), wtime(2)

      real(r8) :: my_wtime

      real(r8), dimension(LBi:UBi,LBj:UBj) :: Awrk

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
!  Import fields from wave model (SWAN) to ocean model (ROMS).
!  Currently, both waves and ocean model grids are the same.
!  We need to revisit this logic to allow interpolation.
!-----------------------------------------------------------------------
!
!  Schedule receiving fields from wave model.
!
      CALL mpi_comm_rank (OCN_COMM_WORLD, MyRank, MyError)
      buffer(1)=my_wtime(wtime)
      CALL MCT_Recv (wav2ocn_AV, ROMStoSWAN, MyError)
      RecvTime=RecvTime+my_wtime(wtime)-buffer(1)
      IF (MyError.ne.0) THEN
        IF (Master) THEN
          WRITE (stdout,10) 'wave model, MyError = ', MyError
        END IF
        exit_flag=2
        RETURN
      END IF
!
!  Receive fields from wave model.
!
      Iimport=0
      DO ifield=1,Nimport(Iocean)
        id=ImportID(Iocean)%val(ifield)
        code=ADJUSTL(Fields(id)%code)
        gtype=Fields(id)%GridType
        scale=Fields(id)%scale
        add_offset=Fields(id)%AddOffset

        SELECT CASE (TRIM(code))

          CASE ('Wdir')                   ! wave direction

            CALL AttrVect_exportRAttr (wav2ocn_AV, TRIM(code), A, Asize)
            Iimport=Iimport+1
            DO i=1,Asize
              A(i)=MAX(0.0_r8,A(i))
            END DO
            scale=deg2rad                 ! degress to radians
            add_offset=0.0_r8
            CALL ROMS_import2d (ng, tile,                               &
     &                          id, gtype, scale, add_offset,           &
     &                          Asize, A,                               &
     &                          IstrR, IendR, JstrR, JendR,             &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fields(id)%ImpMin, Fields(id)%ImpMax,   &
     &                          FORCES(ng)%Dwave,                       &
     &                          status)

          CASE ('Wamp')                   ! significant wave hight

            CALL AttrVect_exportRAttr (wav2ocn_AV, TRIM(code), A, Asize)
            Iimport=Iimport+1
            DO i=1,Asize
              A(i)=MAX(0.0_r8,A(i))
            END DO
            scale=1.0_r8                  ! m
            add_offset=0.0_r8
            CALL ROMS_import2d (ng, tile,                               &
     &                          id, gtype, scale, add_offset,           &
     &                          Asize, A,                               &
     &                          IstrR, IendR, JstrR, JendR,             &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fields(id)%ImpMin, Fields(id)%ImpMax,   &
     &                          FORCES(ng)%Hwave,                       &
     &                          status)

          CASE ('Wlen')                   ! wave length

            CALL AttrVect_exportRAttr (wav2ocn_AV, TRIM(code), A, Asize)
            Iimport=Iimport+1
            DO i=1,Asize
              A(i)=MAX(0.0_r8,A(i))
              IF (A(i).eq.Infinity) THEN
                A(i)=Lwave_max
              END IF
            END DO
            scale=1.0_r8                  ! m
            add_offset=0.0_r8
            CALL ROMS_import2d (ng, tile,                               &
     &                          id, gtype, scale, add_offset,           &
     &                          Asize, A,                               &
     &                          IstrR, IendR, JstrR, JendR,             &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fields(id)%ImpMin, Fields(id)%ImpMax,   &
     &                          FORCES(ng)%Lwave,                       &
     &                          status)

          CASE ('Wptop')                  ! peak surface wave period

            CALL AttrVect_exportRAttr (wav2ocn_AV, TRIM(code), A, Asize)
            Iimport=Iimport+1
            DO i=1,Asize
              A(i)=MAX(0.0_r8,A(i))
            END DO
            scale=1.0_r8
            add_offset=0.0_r8
            CALL ROMS_import2d (ng, tile,                               &
     &                          id, gtype, scale, add_offset,           &
     &                          Asize, A,                               &
     &                          IstrR, IendR, JstrR, JendR,             &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fields(id)%ImpMin, Fields(id)%ImpMax,   &
     &                          FORCES(ng)%Pwave_top,                   &
     &                          status)

          CASE ('Wpbot')                  ! mean bottom wave period

            CALL AttrVect_exportRAttr (wav2ocn_AV, TRIM(code), A, Asize)
            Iimport=Iimport+1
            DO i=1,Asize
              A(i)=MAX(0.0_r8,A(i))
            END DO
            scale=1.0_r8
            add_offset=0.0_r8
            CALL ROMS_import2d (ng, tile,                               &
     &                          id, gtype, scale, add_offset,           &
     &                          Asize, A,                               &
     &                          IstrR, IendR, JstrR, JendR,             &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fields(id)%ImpMin, Fields(id)%ImpMax,   &
     &                          FORCES(ng)%Pwave_bot,                   &
     &                          status)

          CASE ('Wdiss')                  ! wave dissipation

            CALL AttrVect_exportRAttr (wav2ocn_AV, TRIM(code), A, Asize)
            Iimport=Iimport+1
            DO i=1,Asize
              A(i)=MAX(0.0_r8,A(i))
            END DO
            scale=1.0_r8/rho0
            add_offset=0.0_r8
            CALL ROMS_import2d (ng, tile,                               &
     &                          id, gtype, scale, add_offset,           &
     &                          Asize, A,                               &
     &                          IstrR, IendR, JstrR, JendR,             &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fields(id)%ImpMin, Fields(id)%ImpMax,   &
     &                          FORCES(ng)%Wave_dissip,                 &
     &                          status)

          CASE ('Wubot')                  ! bottom orbital velocity

            CALL AttrVect_exportRAttr (wav2ocn_AV, TRIM(code), A, Asize)
            Iimport=Iimport+1
            DO i=1,Asize
              A(i)=MAX(0.0_r8,A(i))
            END DO
            scale=1.0_r8                  ! m/s
            add_offset=0.0_r8
            CALL ROMS_import2d (ng, tile,                               &
     &                          id, gtype, scale, add_offset,           &
     &                          Asize, A,                               &
     &                          IstrR, IendR, JstrR, JendR,             &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fields(id)%ImpMin, Fields(id)%ImpMax,   &
     &                          FORCES(ng)%Ub_swan,                     &
     &                          status)

#ifdef SVENDSEN_ROLLER

          CASE ('Wbrk')                   ! percent wave breaking

            CALL AttrVect_exportRAttr (wav2ocn_AV, TRIM(code), A, Asize)
            Iimport=Iimport+1
            DO i=1,Asize
              A(i)=MAX(0.0_r8,A(i))
            END DO
            scale=1.0_r8
            add_offset=0.0_r8
            CALL ROMS_import2d (ng, tile,                               &
     &                          id, gtype, scale, add_offset,           &
     &                          Asize, A,                               &
     &                          IstrR, IendR, JstrR, JendR,             &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fields(id)%ImpMin, Fields(id)%ImpMax,   &
     &                          FORCES(ng)%Wave_break,                  &
     &                          status)
#endif
        END SELECT
      END DO
!
!-----------------------------------------------------------------------
!  Export fields from ocean (ROMS) to wave (SWAN) model.
!-----------------------------------------------------------------------
!
!  Schedule sending fields to the wave model.
!
      Iexport=0
      DO ifield=1,Nexport(Iocean)
        id=ExportID(Iocean)%val(ifield)
        code=ADJUSTL(Fields(id)%code)
        gtype=Fields(id)%GridType
        scale=Fields(id)%scale
        add_offset=Fields(id)%AddOffset

        SELECT CASE (TRIM(code))

          CASE ('bath')                   ! bathymetry (depth)

            CALL ROMS_export2d (ng, tile,                               &
     &                          id, gtype, scale, add_offset,           &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          GRID(ng)%h,                             &
     &                          Fields(id)%ExpMin, Fields(id)%ExpMax,   &
     &                          Asize, A,                               &
     &                          status)
            CALL AttrVect_importRAttr (ocn2wav_AV, TRIM(code), A, Asize)
            Iexport=Iexport+1

          CASE ('SSH')                    ! free-surface (water level)

            CALL ROMS_export2d (ng, tile,                               &
     &                          id, gtype, scale, add_offset,           &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          OCEAN(ng)%zeta(:,:,KOUT),               &
     &                          Fields(id)%ExpMin, Fields(id)%ExpMax,   &
     &                          Asize, A,                               &
     &                          status)
            CALL AttrVect_importRAttr (ocn2wav_AV, TRIM(code), A, Asize)
            Iexport=Iexport+1

          CASE ('Ubar')                   ! 2D U-momentum

            CALL ROMS_export2d (ng, tile,                               &
     &                          id, gtype, scale, add_offset,           &
     &                          LBi, UBi, LBj, UBj,                     &
#ifdef SOLVE3D
     &                          OCEAN(ng)%u(:,:,N(ng),NOUT),            &
#else
     &                          OCEAN(ng)%ubar(:,:,KOUT),               &
#endif
     &                          Fields(id)%ExpMin, Fields(id)%ExpMax,   &
     &                          Asize, A,                               &
     &                          status)
            CALL AttrVect_importRAttr (ocn2wav_AV, TRIM(code), A, Asize)
            Iexport=Iexport+1

          CASE ('Vbar')                   ! 2D V-momentum

            CALL ROMS_export2d (ng, tile,                               &
     &                          id, gtype, scale, add_offset,           &
     &                          LBi, UBi, LBj, UBj,                     &
#ifdef SOLVE3D
     &                          OCEAN(ng)%v(:,:,N(ng),NOUT),            &
#else
     &                          OCEAN(ng)%vbar(:,:,KOUT),               &
#endif
     &                          Fields(id)%ExpMin, Fields(id)%ExpMax,   &
     &                          Asize, A,                               &
     &                          status)
            CALL AttrVect_importRAttr (ocn2wav_AV, TRIM(code), A, Asize)
            Iexport=Iexport+1

          CASE ('ZO')                   ! bottom roughness

            DO j=JstrR,JendR
              DO i=IstrR,IendR
#ifdef BBL_MODEL
                Awrk(i,j)=MAX(0.0001_r8,                                &
     &                        SEDBED(ng)%bottom(i,j,izNik)*30.0_r8)
#else
                Awrk(i,j)=MAX(0.0001_r8,rdrg2(ng))
#endif
              END DO
            END DO
            CALL ROMS_export2d (ng, tile,                               &
     &                          id, gtype, scale, add_offset,           &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Awrk,                                   &
     &                          Fields(id)%ExpMin, Fields(id)%ExpMax,   &
     &                          Asize, A,                               &
     &                          status)
            CALL AttrVect_importRAttr (ocn2wav_AV, TRIM(code), A, Asize)
            Iexport=Iexport+1

        END SELECT
      END DO
!
!  Send ocean fields to wave model.
!
      IF (Iexport.gt.0) THEN
        buffer(2)=my_wtime(wtime)
        CALL MCT_Send (ocn2wav_AV, ROMStoSWAN, MyError)
        SendTime=SendTime+my_wtime(wtime)-buffer(2)
        IF (MyError.ne.0) THEN
          IF (Master) THEN
            WRITE (stdout,20) 'wave model, MyError = ', MyError
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
 10   FORMAT (' OCN2WAV_COUPLING - error while receiving fields from ', &
     &        a, i4)
 20   FORMAT (' OCN2WAV_COUPLING - error while sending fields to: ',    &
     &        a, i4)
 30   FORMAT (6x,'OCN2WAV   - (', i2.2, ') imported and (', i2.2,       &
     &        ') exported fields,', t62, 't = ', a,/, 16x,              &
     &        '- ROMS coupling exchanges wait clock (s):',/, 19x,       &
     &        '(Recv= ', 1p,e14.8,0p, ' Send= ', 1p,e14.8,0p,')')
 40   FORMAT (16x,'- ',a,a,                                             &
     &        /,19x,'(Min= ',1p,e15.8,0p,' Max= ',1p,e15.8,0p,')')

      RETURN
      END SUBROUTINE ocn2wav_coupling_tile

      SUBROUTINE finalize_ocn2wav_coupling
!
!========================================================================
!                                                                       !
!  This routine finalizes ocean and wave models coupling data streams.  !
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
      CALL Router_clean (ROMStoSWAN, MyError)
      CALL AttrVect_clean (ocn2wav_AV, MyError)
      CALL GlobalSegMap_clean (GSMapROMS, MyError)

      RETURN

      END SUBROUTINE finalize_ocn2wav_coupling
