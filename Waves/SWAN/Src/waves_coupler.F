#include "cppdefs.h"
      MODULE waves_coupler_mod
#ifdef WAVES_OCEAN
!svn $Id$
!==================================================== John C. Warner ===
!  Copyright (c) 2002-2007 The ROMS/TOMS Group      Hernan G. Arango   !
!   Licensed under a MIT/X style license                               !
!   See License_ROMS.txt                                               !
!=======================================================================
!                                                                      !
!  This module is used to communicate and exchange data between SWAN   !
!  other coupled model(s) using the Model Coupling Toolkit (MCT).      !
!                                                                      !
!=======================================================================
!
!  Componenet model registry.
!
      USE m_MCTWorld, ONLY : MCTWorld_init => init
      USE m_MCTWorld, ONLY : MCTWorld_clean => clean
!
!  Domain decompositin descriptor datatype and assocoiated methods.
!
      USE m_GlobalSegMap, ONLY : GlobalSegMap
      USE m_GlobalSegMap, ONLY : GlobalSegMap_init => init
      USE m_GlobalSegMap, ONLY : GlobalSegMap_lsize => lsize
      USE m_GlobalSegMap, ONLY : GlobalSegMap_clean => clean
      USE m_GlobalSegMap, ONLY : GlobalSegMap_Ordpnts => OrderedPoints
!
!  Field storage data types and associated methods.
!
      USE m_AttrVect, ONLY : AttrVect
      USE m_AttrVect, ONLY : AttrVect_init => init
      USE m_AttrVect, ONLY : AttrVect_zero => zero
      USE m_AttrVect, ONLY : AttrVect_clean => clean
      USE m_AttrVect, ONLY : AttrVect_indxR => indexRA
      USE m_AttrVect, ONLY : AttrVect_importRAttr => importRAttr
      USE m_AttrVect, ONLY : AttrVect_exportRAttr => exportRAttr
!
!  Intercomponent communitcations scheduler.
!
      USE m_Router, ONLY : Router
      USE m_Router, ONLY : Router_init => init
      USE m_Router, ONLY : Router_clean => clean
!
!  Intercomponent transfer.
!
      USE m_Transfer, ONLY : MCT_Send => send
      USE m_Transfer, ONLY : MCT_Recv => recv
!
      implicit none
!
      PRIVATE

      PUBLIC :: initialize_ocean_coupling
      PUBLIC :: ocean_coupling
      PUBLIC :: finalize_ocean_coupling

      include 'mpif.h'
!
!  Declarations.
!
      TYPE(GlobalSegMap) :: GSMapSWAN         ! GloabalSegMap variables
 
      TYPE(AttrVect) :: ToOceanAV             ! AttrVect variables
      TYPE(AttrVect) :: FromOceanAV
 
      type(Router) :: RoutSWAN                ! Router variables

      CONTAINS

      SUBROUTINE initialize_ocean_coupling
!
!=======================================================================
!                                                                      !
!  Initialize waves and ocean models coupling stream.  This is the     !
!  training phase use to constuct  MCT  parallel interpolators and     !
!  stablish communication patterns.                                    !
!                                                                      !
!=======================================================================
!
      USE OCPCOMM4                                                        40.41
      USE SWCOMM3                                                         40.41
      USE M_GENARR                                                        40.31
      USE M_PARALL                                                        40.31
!
      include 'mpif.h'
!
!  Local variable declarations.  
!
      integer :: MyColor, MyCOMM, MyError, MyRank, MyKey, MyValue
      integer :: npoints, gsmsize, nprocs, localsize
      integer :: j, Isize, Jsize

      integer, pointer :: start(:), length(:)

!-----------------------------------------------------------------------
!  Begin initialization phase.
!-----------------------------------------------------------------------
!
! Get local rank and size.
!
      CALL mpi_comm_rank (WAV_COMM_WORLD, MyRank, MyError)
      CALL mpi_comm_size (WAV_COMM_WORLD, nprocs, MyError)
!
! Initialize MCTworld.
!
      CALL MCTWorld_init (ncomps,MPI_COMM_WORLD,WAV_COMM_WORLD,WavId)
!
! Initialize a Global Segment Map for non-haloed transfer of data out of
! SWAN. Determine non-haloed start and length arrays for this processor.
!
      IF (nprocs.eq.1) THEN
        Isize=MXCGL
        Jsize=MYCGL
      ELSE
        IF (MXCGL.gt.MYCGL) THEN
          Isize=MXC-IHALOX*IBLKAD(1)
          Jsize=MYC
        ELSE
          Isize=MXC
          Jsize=MYC-IHALOY*IBLKAD(1)
        END IF
      END IF
!
      allocate ( start(Jsize) )
      allocate ( length(Jsize) )
!
      DO j=1,Jsize
        length(j)=ISIZE
        IF (MXCGL.gt.MYCGL) THEN
          IF (MyRank.eq.0) THEN
            start(j)=MXF+(j-1)*MXCGL
          ELSE
            start(j)=MXF+(j-1)*MXCGL+IHALOX
          ENDIF
        ELSE
          IF (MyRank.eq.0) THEN
            start(j)=MYF+(j-1)*MXCGL
          ELSE
            start(j)=(MYF+IHALOY-1)*MXCGL+1+(j-1)*MXCGL
          ENDIF
        ENDIF
      END DO
!
      gsmsize=Isize*Jsize
!  
      CALL GlobalSegMap_init(GSMapSWAN,start,length,0,                  &
     &                       WAV_COMM_WORLD,WavId)
!
! Initialize Attribute Vector ToOceanAv to hold the data sent to ROMS.
!
      CALL AttrVect_init (ToOceanAV,                                    &
     &     rlist="DISSIP:HSIGN:RTP:SETUP:TMBOT:UBOT:DIR:WLEN:TM01:QB",  &
     &     lsize=gsmsize)
      CALL AttrVect_zero (ToOceanAV)
!
! Initialize Attribute Vector FromOceanAV that will have ROMS data in
! it.
!
      CALL AttrVect_init (FromOceanAV,rList="DEPTH:WLEV:VELX:VELY",     &
     &                    lsize=gsmsize)
      CALL AttrVect_zero (FromOceanAV)
!
! Initialize a router to the Waves component.
!
      CALL Router_init (OcnId,GSMapSWAN,WAV_COMM_WORLD,RoutSWAN)
!
      deallocate (start)
      deallocate (length)

      RETURN
      END SUBROUTINE initialize_ocean_coupling

      SUBROUTINE ocean_coupling (MIP, NVOQP, VOQR, VOQ, IRQ, IVTYPE,    &
     &                           COMPDA)
!
!=======================================================================
!                                                                      !
!  This subroutine reads and writes the coupling data streams between  !
!  ocean and wave models. Currently, the following data streams are    !
!  processed:                                                          !
!                                                                      !
!  Fields sent to the OCEAN Model:                                     !
!                                                                      !
!     * Dwave      Wave direction.                                     !
!     * Hwave      Wave height.                                        !
!     * Lwave      Wave length.                                        !
!     * Pwave_bot  Wave bottom period.                                 !
!     * Pwave_top  Wave surface period.                                !
!     * Wave_break Percent of breakig waves.                           !
!     * Wave_dissip Wave energy dissipation.                           !
!                                                                      !
!  Fields acquired from the OCEAN Model:                               !
!                                                                      !
!     * ubar       Depth integrated xi-direction velocity.             !
!     * vbar       Depth integrated eta-direction velocity.            !
!     * zeta       Water surface elevation.                            !
!     * h          Bottom elevation.                                   !
!                                                                      !
!=======================================================================
!
      USE SWCOMM3                                                         40.41
      USE SWCOMM4                                                         40.41
      USE OUTP_DATA                                                       40.13
      USE M_PARALL                                                        40.51
      USE M_GENARR
      USE SWCOMM1                                                         40.41
!
      implicit none
!
!  Imported variable declarations.
!
      integer :: MIP, IRQ, nvoqp
      integer :: VOQR(NMOVAR), IVTYPE, IP, IX, IY

      real :: COMPDA(MCGRD,MCMVAR)
      real :: VOQ(MIP,NVOQP)
!
!  Local variable declarations.
!
      integer :: MyStatus, i, j, gsmsize, ierr, MyRank
      integer :: MyError, MySize, indx, Istr, Iend, Jstr, Jend
      integer :: Isize, Jsize, INDXG, NPROCS, OFFSET
      integer :: NUMTRANSFER, NNEIGH, HALOSIZE, NUMSENT, INB
      integer :: WHICHWAY, GDEST, GSRC, TAGT, TAGB, TAGR, TAGL
      integer :: TREQUEST,BREQUEST,RREQUEST,LREQUEST,MSIZE

      integer, dimension(MPI_STATUS_SIZE,4) :: status
      integer, pointer :: points(:)

      real :: cff
      real, pointer :: avdata(:), TEMPMCT(:,:)
      real, pointer :: GRECVT(:), GRECVB(:), GRECVR(:), GRECVL(:)
      real, pointer :: GSENDT(:), GSENDB(:), GSENDR(:), GSENDL(:)
!
!-----------------------------------------------------------------------
!  Send wave fields to ROMS.
!-----------------------------------------------------------------------
!
      CALL MPI_COMM_RANK (WAV_COMM_WORLD, MyRank, MyError)
      CALL MPI_COMM_SIZE (WAV_COMM_WORLD, NPROCS, MyError)
!
!  Size is the number of grid point on this processor.
!
      gsmsize=GlobalSegMap_lsize(GSMapSWAN,WAV_COMM_WORLD)
!
!  Load wave data into Attribute Vector Array ToOceanAV.
!
      allocate ( avdata(gsmsize),stat=ierr )
      allocate ( points(gsmsize),stat=ierr )
      avdata=0.0
      points=0
!
!  Ask for points in this tile.
!
      CALL GlobalSegMap_Ordpnts (GSMapSWAN,MyRank,points)
!
!  Load output data into temporary array.
!
      DO IP=1,gsmsize
        avdata(IP)=VOQ(points(IP),VOQR(IVTYPE))
      END DO
!
!  Load output date into MCT attribute.
!
      IF (IVTYPE.eq.7) THEN
        CALL AttrVect_importRAttr (ToOceanAV,"DISSIP",avdata)
      ELSE IF (IVTYPE.eq.10) THEN
        CALL AttrVect_importRAttr (ToOceanAV,"HSIGN",avdata)
      ELSE IF (IVTYPE.eq.12) THEN
        CALL AttrVect_importRAttr (ToOceanAV,"RTP",avdata)
      ELSE IF (IVTYPE.eq.11) THEN
        CALL AttrVect_importRAttr (ToOceanAV,"TM01",avdata)
      ELSE IF (IVTYPE.eq.39) THEN
        CALL AttrVect_importRAttr (ToOceanAV,"SETUP",avdata)
      ELSE IF (IVTYPE.eq.50) THEN
        CALL AttrVect_importRAttr (ToOceanAV,"TMBOT",avdata)
      ELSE IF (IVTYPE.eq.6) THEN
        CALL AttrVect_importRAttr (ToOceanAV,"UBOT",avdata)
      ELSE IF (IVTYPE.eq.13) THEN
        CALL AttrVect_importRAttr (ToOceanAV,"DIR",avdata)
      ELSE IF (IVTYPE.eq.17) THEN
        CALL AttrVect_importRAttr (ToOceanAV,"WLEN",avdata)
      ELSE IF (IVTYPE.eq.8) THEN
        CALL AttrVect_importRAttr (ToOceanAV,"QB",avdata)
      END IF
!
!-----------------------------------------------------------------------
!  Create a restart file.
!-----------------------------------------------------------------------
!
      IF (IRQ.eq.NREOQ) THEN
        CALL BACKUP (AC2, SPCSIG, SPCDIR, KGRPNT, XCGRID, YCGRID)
!
!-----------------------------------------------------------------------
!  Send wave parameters to ROMS.
!-----------------------------------------------------------------------
!
        CALL MCT_SEND (ToOceanAV,RoutSWAN,MyError)
        IF (MYRANK.EQ.0) THEN
          WRITE (*,*) '== SWAN sent wave fields and Myerror= ',MyError
        ENDIF
        IF (MyError.ne.0) THEN
          WRITE (*,*) 'coupling send fail swancplr, MyStatus= ', MyError
          CALL finalize_ocean_coupling ("Coupling failed swancplr")
        END IF
!
!-----------------------------------------------------------------------
!  Receive from ROMS: Depth, Water Level, VELX, and VELY.
!-----------------------------------------------------------------------
!
        NUMTRANSFER=4
        CALL MCT_Recv (FromOceanAV, RoutSWAN, MyError)
        IF (MYRANK.EQ.0) THEN
          WRITE (*,*) '== SWAN recvd ocean fields and Myerror= ',       &
     &                MyError
        END IF
        IF (MyError.ne.0) THEN
          WRITE (*,*) 'coupling fail swancplr, MyStatus= ', MyError
          CALL finalize_ocean_coupling ("Coupling failed swancplr")
        END IF
!
! Pass the non-halo data from MCT into tempmct array.
!
        NNEIGH = IBLKAD(1)
        IF (nprocs.eq.1) THEN
          Istr=1
          Iend=MXC
          Jstr=1
          Jend=MYC
        ELSE
          IF (MXCGL.GT.MYCGL) THEN
            IF (MyRank.eq.0) THEN
              Istr=1
            ELSE
              Istr=IHALOX+1
            ENDIF
            Isize=MXC-IHALOX*IBLKAD(1)
            Iend=Istr+Isize-1
            Jstr=1
            Jend=MYC
            HALOSIZE=IHALOX*MYC
          ELSE
            IF (MyRank.eq.0) THEN
              Jstr=1
            ELSE
              Jstr=IHALOY+1
            ENDIF
            Jsize=MYC-IHALOY*IBLKAD(1)
            Jend=Jstr+Jsize-1
            Istr=1
            Iend=MXC
            HALOSIZE=IHALOY*MXC
          END IF
        END IF
!
        allocate ( TEMPMCT(MXC*MYC,NUMTRANSFER),stat=ierr )
        TEMPMCT=0.0
!
!  Bottom elevation.
!
        CALL AttrVect_exportRAttr (FromOceanAV,"DEPTH",avdata,gsmsize)
        IP=0
        DO IY=Jstr,Jend
          DO IX=Istr,Iend
            IP=IP+1
            INDXG=(IY-1)*MXC+IX
            TEMPMCT(INDXG,1)=avdata(IP)
          END DO
        END DO
!
!  Water surface elevation.
!
        CALL AttrVect_exportRAttr (FromOceanAV,"WLEV",avdata,gsmsize)
        IP=0
        DO IY=Jstr,Jend
          DO IX=Istr,Iend
            IP=IP+1
            INDXG=(IY-1)*MXC+IX
            TEMPMCT(INDXG,2)=avdata(IP)
          END DO
        END DO
!
!  Depth-integrated u-velocity.
!
        CALL AttrVect_exportRAttr(FromOceanAV,"VELX",avdata,gsmsize)
        IP=0
        DO IY=Jstr,Jend
          DO IX=Istr,Iend
            IP=IP+1
            INDXG=(IY-1)*MXC+IX
            TEMPMCT(INDXG,3)=avdata(IP)
          END DO
        END DO
!
!  Depth-integrated v-velocity.
!
        CALL AttrVect_exportRAttr(FromOceanAV,"VELY",avdata,gsmsize)
        IP=0
        DO IY=Jstr,Jend
          DO IX=Istr,Iend
            IP=IP+1
            INDXG=(IY-1)*MXC+IX
            TEMPMCT(INDXG,4)=avdata(IP)
          END DO
        END DO

        IF (NPROCS.GT.1) THEN
!
!  Pack and send halo regions to be exchanged with adjacent tiles.
!  IBLKAD contains the tile data.
!  WHICHWAY: [top, bot, right, left] = [1 2 3 4]
!
          MSIZE=HALOSIZE*NUMTRANSFER
          IF (MXCGL.GT.MYCGL) THEN
            allocate ( GSENDR(MSIZE),stat=ierr )
            allocate ( GSENDL(MSIZE),stat=ierr )
            allocate ( GRECVR(MSIZE),stat=ierr )
            allocate ( GRECVL(MSIZE),stat=ierr )
            GSENDR=0.0
            GSENDL=0.0
            GRECVR=0.0
            GRECVL=0.0
          ELSE
            allocate ( GSENDT(MSIZE),stat=ierr )
            allocate ( GSENDB(MSIZE),stat=ierr )
            allocate ( GRECVT(MSIZE),stat=ierr )
            allocate ( GRECVB(MSIZE),stat=ierr )
            GSENDT=0.0
            GSENDB=0.0
            GRECVT=0.0
            GRECVB=0.0
          ENDIF
          TAGT=1
          TAGB=2
          TAGR=3
          TAGL=4
          DO INB=1,NNEIGH
            OFFSET=0
            WHICHWAY=IBLKAD(3*INB)
            DO NUMSENT=1,NUMTRANSFER
              IP=OFFSET
              IF (WHICHWAY.EQ.1) THEN
                DO IY=MYC-IHALOX-2,MYC-3
                  DO IX=1,MXC
                    IP=IP+1
                    INDXG=(IY-1)*MXC+IX
                    GSENDT(IP)=TEMPMCT(INDXG,NUMSENT)
                  END DO
                END DO
              ELSE IF (WHICHWAY.EQ.2) THEN
                DO IY=IHALOY+1,IHALOY+3
                  DO IX=1,MXC
                    IP=IP+1
                    INDXG=(IY-1)*MXC+IX
                    GSENDB(IP)=TEMPMCT(INDXG,NUMSENT)
                  END DO
                END DO
              ELSE IF (WHICHWAY.EQ.3) THEN
                DO IY=1,MYC
                  DO IX=MXC-IHALOX-2,MXC-3
                    IP=IP+1
                    INDXG=(IY-1)*MXC+IX
                    GSENDR(IP)=TEMPMCT(INDXG,NUMSENT)
                  END DO
                END DO
              ELSE IF (WHICHWAY.EQ.4) THEN
                DO IY=1,MYC
                  DO IX=IHALOX+1,IHALOX+3
                    IP=IP+1
                    INDXG=(IY-1)*MXC+IX
                    GSENDL(IP)=TEMPMCT(INDXG,NUMSENT)
                  END DO
                END DO
              ENDIF
              OFFSET=OFFSET+HALOSIZE
            END DO
          END DO
          DO INB=1,NNEIGH
            GSRC=IBLKAD(3*INB-1)-1
            WHICHWAY=IBLKAD(3*INB)
            IF (WHICHWAY.EQ.1) THEN
              CALL mpi_irecv (GRECVT,MSIZE,SWREAL,                      &
     &                        GSRC,TAGB,WAV_COMM_WORLD,TREQUEST,MyError)
            ELSE IF (WHICHWAY.EQ.2) THEN
              CALL mpi_irecv (GRECVB,MSIZE,SWREAL,                      &
     &                        GSRC,TAGT,WAV_COMM_WORLD,BREQUEST,MyError)
            ELSE IF (WHICHWAY.EQ.3) THEN
              CALL mpi_irecv (GRECVR,MSIZE,SWREAL,                      &
     &                        GSRC,TAGL,WAV_COMM_WORLD,RREQUEST,MyError)
            ELSE IF (WHICHWAY.EQ.4) THEN
              CALL mpi_irecv (GRECVL,MSIZE,SWREAL,                      &
     &                        GSRC,TAGR,WAV_COMM_WORLD,LREQUEST,MyError)
            END IF
          END DO
          DO INB=1,NNEIGH
            GDEST=IBLKAD(3*INB-1)-1
            WHICHWAY=IBLKAD(3*INB)
            IF (WHICHWAY.EQ.1) THEN
              CALL mpi_send (GSENDT,MSIZE,SWREAL,                       &
     &                       GDEST,TAGT,WAV_COMM_WORLD,MyError)
            ELSE IF (WHICHWAY.EQ.2) THEN
              CALL mpi_send (GSENDB,MSIZE,SWREAL,                       &
     &                       GDEST,TAGB,WAV_COMM_WORLD,MyError)
            ELSE IF (WHICHWAY.EQ.4) THEN
              CALL mpi_send (GSENDL,MSIZE,SWREAL,                       &
     &                       GDEST,TAGL,WAV_COMM_WORLD,MyError)
            ELSE IF (WHICHWAY.EQ.3) THEN
              CALL mpi_send (GSENDR,MSIZE,SWREAL,                       &
     &                       GDEST,TAGR,WAV_COMM_WORLD,MyError)
            END IF
          END DO
!
! Receive and unpack halo regions exchanged with adjacent tiles.
! [top, bot, right, left] = [1 2 3 4]
!
          DO INB=1,NNEIGH
            WHICHWAY=IBLKAD(3*INB)
            IF (WHICHWAY.EQ.1) THEN
              CALL mpi_wait (TREQUEST,status(1,1),MyError)
            ELSE IF (WHICHWAY.EQ.2) THEN
              CALL mpi_wait (BREQUEST,status(1,2),MyError)
            ELSE IF (WHICHWAY.EQ.3) THEN
              CALL mpi_wait (RREQUEST,status(1,3),MyError)
            ELSE IF (WHICHWAY.EQ.4) THEN
              CALL mpi_wait (LREQUEST,status(1,4),MyError)
            ENDIF
          END DO
!
          DO INB=1,NNEIGH
            OFFSET=0
            WHICHWAY=IBLKAD(3*INB)
            IF (WHICHWAY.EQ.1) THEN
              DO NUMSENT=1,NUMTRANSFER
                IP=OFFSET
                DO IY=MYC-2,MYC
                  DO IX=1,MXC
                    IP=IP+1
                    INDXG=(IY-1)*MXC+IX
                    TEMPMCT(INDXG,NUMSENT)=GRECVT(IP)
                  END DO
                END DO
                OFFSET=OFFSET+HALOSIZE
              END DO
            ELSE IF (WHICHWAY.EQ.2) THEN
              DO NUMSENT=1,NUMTRANSFER
                IP=OFFSET
                DO IY=1,IHALOY
                  DO IX=1,MXC
                    IP=IP+1
                    INDXG=(IY-1)*MXC+IX
                    TEMPMCT(INDXG,NUMSENT)=GRECVB(IP)
                  END DO
                END DO
                OFFSET=OFFSET+HALOSIZE
              END DO
            ELSE IF (WHICHWAY.EQ.3) THEN
              DO NUMSENT=1,NUMTRANSFER
                IP=OFFSET
                DO IY=1,MYC
                  DO IX=MXC-2,MXC
                    IP=IP+1
                    INDXG=(IY-1)*MXC+IX
                    TEMPMCT(INDXG,NUMSENT)=GRECVR(IP)
                  END DO
                END DO
                OFFSET=OFFSET+HALOSIZE
              END DO
            ELSE IF (WHICHWAY.EQ.4) THEN
              DO NUMSENT=1,NUMTRANSFER
                IP=OFFSET
                DO IY=1,MYC
                  DO IX=1,IHALOX
                    IP=IP+1
                    INDXG=(IY-1)*MXC+IX
                    TEMPMCT(INDXG,NUMSENT)=GRECVL(IP)
                  END DO
                END DO
                OFFSET=OFFSET+HALOSIZE
              END DO
            ENDIF
          END DO
          IF (MXCGL.GT.MYCGL) THEN
            deallocate (GRECVR,GRECVL,GSENDR,GSENDL)
          ELSE
            deallocate (GRECVT,GRECVB,GSENDT,GSENDB)
          END IF
        ENDIF 
!
! Finally insert the full (MXC*MYC) TEMPMCT array into the SWAN
! array for DEPTH and computational array COMPDA. Only insert
! active (wet points) using array KGRPNT.
!
!
!  Insert depth into SWAN array.
!
        IP=0
        DO IY = MYF,MYL
          DO IX = MXF,MXL
            IP=IP+1
            INDX=KGRPNT(IX-MXF+1,IY-MYF+1)
            IF (INDX.GT.1) THEN
              DEPTH(IX,IY)=TEMPMCT(IP,1)
            END IF
          END DO
        END DO
!
! Move values at 'present' time level 2 to 'old' time level 1.
! MCGRD = MXC*MYC+1-#masked cells.
! MXC = # cells x-dir in this tile including halox.
! MYC = # cells y-dir in this tile including haloy.
! COMPDA has only active wet points + 1.
!
        DO INDX = 2, MCGRD
          COMPDA(INDX,JWLV1)=COMPDA(INDX,JWLV2)
          COMPDA(INDX,JVX1) =COMPDA(INDX,JVX2)
          COMPDA(INDX,JVY1) =COMPDA(INDX,JVY2)
        ENDDO
!
! Insert water level, velx, and vely into SWAN arrays.
!
        IP=0
        DO IY=1,MYC
          DO IX=1,MXC
            IP=IP+1
            INDX = KGRPNT(IX,IY)
            IF (INDX.GT.1) THEN
              COMPDA(INDX,JWLV2)=REAL(TEMPMCT(IP,2))
              COMPDA(INDX,JVX2)=REAL(TEMPMCT(IP,3))
              COMPDA(INDX,JVY2)=REAL(TEMPMCT(IP,4))
            END IF
          END DO
        END DO
!
        deallocate (TEMPMCT)
      END IF
      deallocate (avdata, points)

      RETURN
      END SUBROUTINE ocean_coupling

      SUBROUTINE finalize_ocean_coupling (string)
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
      CALL Router_clean (RoutSWAN)
      CALL AttrVect_clean (ToOceanAV)
      CALL AttrVect_clean (FromOceanAV)
      CALL GlobalSegMap_clean (GSMapSWAN)
      CALL MCTWorld_clean ()
      CALL mpi_finalize (MyStatus)

      STOP
      END SUBROUTINE finalize_ocean_coupling
#endif
      END MODULE waves_coupler_mod
