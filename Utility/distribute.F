#include "cppdefs.h"
#ifdef DISTRIBUTE
# define ALLGATHER
# undef  ALLREDUCE

      SUBROUTINE mp_barrier (ng)
!
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!=======================================================================
!                                                                      !
!  This routine bloks the caller until all group members have called   !
!  it.                                                                 !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      integer :: MyError
!
!-----------------------------------------------------------------------
!  Synchronize all distribute-memory nodes in the group.
!-----------------------------------------------------------------------
!
# ifdef MPI
      CALL mpi_barrier (OCN_COMM_WORLD, MyError)
# endif

      RETURN
      END SUBROUTINE mp_barrier

      SUBROUTINE mp_bcastf (ng, A, Asize)
!
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!=======================================================================
!                                                                      !
!  This routine broadcasts a floating-point variable to all processors !
!  in the group. It is called by all the members in the group.         !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     A          Variable to broadcast.                                !
!     Asize      Number of entries to broadcast.                       !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A          Broadcasted variable.                                 !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Asize

      real(r8), intent(inout) :: A(Asize)
!
!  Local variable declarations
!
      integer :: MyError
# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, 39)
# endif
!
!-----------------------------------------------------------------------
!  Broadcast requested variable.
!-----------------------------------------------------------------------
!
# ifdef MPI
      CALL mpi_bcast (A, Asize, MP_FLOAT, MyMaster, OCN_COMM_WORLD,     &
     &                MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        WRITE (stdout,10) 'BCAST', MyRank, MyError
 10     FORMAT (/,' MP_BCASTF - error during ',a,' call, Node = ',i3.3, &
     &            ' Error = ',i3)
        exit_flag=2
      END IF
# endif
# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, 39)
# endif

      RETURN
      END SUBROUTINE mp_bcastf

      SUBROUTINE mp_bcasti (ng, A, Asize)
!
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!=======================================================================
!                                                                      !
!  This routine broadcasts an integer variable to all processors in    !
!  the group. It is called by all the members in the group.            !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     A          Variable to broadcast.                                !
!     Asize      Number of entries to broadcast.                       !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A          Broadcasted variable.                                 !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Asize

      integer, intent(inout) :: A(Asize)
!
!  Local variable declarations
!
      integer :: MyError

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, 39)
# endif
!
!-----------------------------------------------------------------------
!  Broadcast requested variable.
!-----------------------------------------------------------------------
!
# ifdef MPI
      CALL mpi_bcast (A, Asize, MPI_INTEGER, MyMaster, OCN_COMM_WORLD,  &
     &                MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        WRITE (stdout,10) 'BCAST', MyRank, MyError
 10     FORMAT (/,' MP_BCASTI - error during ',a,' call, Node = ',i3.3, &
     &            ' Error = ',i3)
        exit_flag=2
      END IF
# endif
# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, 39)
# endif

      RETURN
      END SUBROUTINE mp_bcasti

      SUBROUTINE mp_bcastl (ng, A, Asize)
!
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!=======================================================================
!                                                                      !
!  This routine broadcasts a logical variable to all processors in     !
!  the group. It is called by all the members in the group.            !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     A          Variable to broadcast.                                !
!     Asize      Number of entries to broadcast.                       !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A          Broadcasted variable.                                 !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Asize

      logical, intent(inout) :: A(Asize)
!
!  Local variable declarations
!
      integer :: MyError

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, 39)
# endif
!
!-----------------------------------------------------------------------
!  Broadcast requested variable.
!-----------------------------------------------------------------------
!
# ifdef MPI
      CALL mpi_bcast (A, Asize, MPI_LOGICAL, MyMaster, OCN_COMM_WORLD,  &
     &                MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        WRITE (stdout,10) 'BCAST', MyRank, MyError
 10     FORMAT (/,' MP_BCASTL - error during ',a,' call, Node = ',i3.3, &
     &            ' Error = ',i3)
        exit_flag=2
      END IF
# endif
# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, 39)
# endif

      RETURN
      END SUBROUTINE mp_bcastl

      SUBROUTINE mp_bcasts (ng, A, Asize)
!
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!=======================================================================
!                                                                      !
!  This routine broadcasts a string variable to all processors in the  !
!  group. It is called by all the members in the group.                !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     A          Variable to broadcast.                                !
!     Asize      Number of entries to broadcast.                       !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A          Broadcasted variable.                                 !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Asize

      character (len=*), intent(inout) :: A
!
!  Local variable declarations
!
      integer :: MyError

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, 39)
# endif
!
!-----------------------------------------------------------------------
!  Broadcast requested variable.
!-----------------------------------------------------------------------
!
# ifdef MPI
      CALL mpi_bcast (A, Asize, MPI_BYTE, MyMaster, OCN_COMM_WORLD,     &
     &                MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        WRITE (stdout,10) 'BCAST', MyRank, MyError
 10     FORMAT (/,' MP_BCASTS - error during ',a,' call, Node = ',i3.3, &
     &            ' Error = ',i3)
        exit_flag=2
      END IF
# endif
# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, 39)
# endif

      RETURN
      END SUBROUTINE mp_bcasts

      SUBROUTINE mp_boundary (ng, Imin, Imax, LBi, UBi, LBk, UBk,       &
     &                        update, A)
!
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!=======================================================================
!                                                                      !
!  This routine exchanges boundary arrays between tiles.               !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     Imin       Starting tile index.                                  !
!     Imax       Ending   tile index.                                  !
!     Jstr       Starting tile index in the J-direction.               !
!     Jend       Ending   tile index in the J-direction.               !
!     LBi        I-dimension Lower bound.                              !
!     UBi        I-dimension Upper bound.                              !
!     LBk        K-dimension Lower bound, if any. Otherwise, a value   !
!                  of one is expected.                                 !
!     LBk        K-dimension Upper bound, if any. Otherwise, a value   !
!                  of one is expected.                                 !
!     UBk        K-dimension Upper bound.                              !
!     update     Switch activated by the node that updated the         !
!                  boundary data.                                      !
!     A          Boundary array (1D or 2D) to process.                 !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A          Updated boundary array (1D or 2D).                    !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
!
      implicit none
!
!  Imported variable declarations.
!
      logical, intent(in) :: update

      integer, intent(in) :: ng, Imin, Imax, LBi, UBi, LBk, UBk

      real(r8), intent(inout) :: A(LBi:UBi,LBk:UBk)
!
!  Local variable declarations.
!
      integer :: Ilen, Ioff, MyError, Nnodes, Npts, i, ik, k, kc, rank

      real(r8), dimension((UBi-LBi+1)*(UBk-LBk+1)) :: Asend
      real(r8), dimension((UBi-LBi+1)*(UBk-LBk+1),                      &
     &                    0:NtileI(ng)*NtileJ(ng)-1) :: Arecv

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, 43)
# endif
!
!-----------------------------------------------------------------------
!  Pack boundary data.  Zero-out boundary array except points updated
!  by the appropriate node, so sum reduction can be perfomed during
!  unpacking.
!-----------------------------------------------------------------------
!
      Ilen=UBi-LBi+1
      Ioff=1-LBi
      Npts=Ilen*(UBk-LBk+1)
      Asend(1:Npts)=0.0_r8
      IF (update) THEN
        DO k=LBk,UBk
          kc=(k-LBk)*Ilen
          DO i=Imin,Imax
            ik=i+Ioff+kc
            Asend(ik)=A(i,k)
          END DO
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Collect data from all nodes.
!-----------------------------------------------------------------------
!
# ifdef MPI
      CALL mpi_allgather (Asend, Npts, MP_FLOAT, Arecv, Npts, MP_FLOAT, &
     &                    OCN_COMM_WORLD, MyError)
# endif
!
!-----------------------------------------------------------------------
!  Unpack data: reduction sum.
!-----------------------------------------------------------------------
!
      Nnodes=NtileI(ng)*NtileJ(ng)-1
      ik=0
      DO k=LBk,UBk
        DO i=LBi,UBi
          A(i,k)=0.0_r8
          ik=ik+1
          DO rank=0,Nnodes
            A(i,k)=A(i,k)+Arecv(ik,rank)
          END DO
        END DO
      END DO
# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, 43)
# endif

      RETURN
      END SUBROUTINE mp_boundary

      SUBROUTINE mp_collect (ng, Npts, Aspv, A)
!
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!=======================================================================
!                                                                      !
!  This routine collects requested buffer from all members in the      !
!  group. Then, it packs distributed data by removing the special      !
!  values. This routine is used when extracting station data from      !
!  tiled arrays.                                                       !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     Npts       Number of extracted data points.                      !
!     Aspv       Special value indicating no data.  This implies that  !
!                  desired data is tile unbouded.                      !
!     A          Extracted data.                                       !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A          Collected data.                                       !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Npts

      real(r8), intent(in) :: Aspv

      real(r8), intent(inout) :: A(Npts)
!
!  Local variable declarations.
!
      integer :: MyError, Nnodes, i, rank

      real(r8), dimension(Npts,0:NtileI(ng)*NtileJ(ng)-1) :: Arecv

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, 44)
# endif
!
!-----------------------------------------------------------------------
!  Collect data from all nodes.
!-----------------------------------------------------------------------
!
      CALL mpi_allgather (A, Npts, MP_FLOAT, Arecv, Npts, MP_FLOAT,     &
     &                    OCN_COMM_WORLD, MyError)
!
!  Pack data according to special values: sum or ignore.
!
      Nnodes=NtileI(ng)*NtileJ(ng)-1
      IF (Aspv.eq.0.0_r8) THEN
        DO i=1,Npts
          A(i)=0.0_r8
          DO rank=0,Nnodes
            A(i)=A(i)+Arecv(i,rank)
          END DO
        END DO
      ELSE
        DO i=1,Npts
          DO rank=0,Nnodes
            IF (Arecv(i,rank).ne.Aspv) THEN
              A(i)=Arecv(i,rank)
            END IF
          END DO
        END DO
      END IF
# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, 44)
# endif

      RETURN
      END SUBROUTINE mp_collect

      SUBROUTINE mp_exchange (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj, LBk, UBk,             &
     &                        Nghost, EWperiodic, NSperiodic, A)
!
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================ Enrique Curchitser ===
!                                                                      !
!  This routine updates the tile overlap halo of a  2D or 3D real      !
!  array. It exchanges the specified number of "ghost-points". In      !
!  order to minimize the number send and recive calls,  the ghost-     !
!  points are included in the buffers. Therefore, the order of the     !
!  pack/send/recive/unpack is crucial.                                 !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     Istr       Starting tile index in the I-direction.               !
!     Iend       Ending   tile index in the I-direction.               !
!     Jstr       Starting tile index in the J-direction.               !
!     Jend       Ending   tile index in the J-direction.               !
!     LBi        I-dimension Lower bound.                              !
!     UBi        I-dimension Upper bound.                              !
!     LBj        J-dimension Lower bound.                              !
!     UBj        J-dimension Upper bound.                              !
!     LBk        K-dimension Lower bound, if any. Otherwise, a value   !
!                  of one is expected.                                 !
!     LBk        K-dimension Upper bound, if any. Otherwise, a value   !
!                  of one is expected.                                 !
!     UBk        K-dimension Upper bound.                              !
!     Nghost     Number of ghost-points in the halo region.            !
!     EWperiodic Switch indicating EW periodicity exchanges.           !
!     NSperiodic Switch indicating NS periodicity exchanges.           !
!     A          Tile array (2D or 3D) to process.                     !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A          Updated tile array (2D or 3D).                        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
!
      implicit none
!
!  Imported variable declarations.
!
      logical, intent(in) :: EWperiodic, NSperiodic

      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk, Nghost

      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj,LBk:UBk)
!
!  Local variable declarations.
!
      logical :: Wexchange, Sexchange, Eexchange, Nexchange

      integer :: i, ikS, ikN, Imin, Imax, Ilen, IKlen, IKsizeS, IKsizeN
      integer :: j, jkW, jkE, Jmin, Jmax, Jlen, JKlen, JKsizeW, JKsizeE
      integer :: k, kc,  m, mc, Klen, MyRankI, MyRankJ, rank
      integer :: Wtile, GsendW, GrecvW, Wtag, Werror, Wrequest
      integer :: Stile, GsendS, GrecvS, Stag, Serror, Srequest
      integer :: Etile, GsendE, GrecvE, Etag, Eerror, Erequest
      integer :: Ntile, GsendN, GrecvN, Ntag, Nerror, Nrequest
      integer :: EWsize
      integer :: NSsize
      integer :: Null_Value

      integer, dimension(-1:NtileI(ng),-1:NtileJ(ng)) :: table
# ifdef MPI
      integer, dimension(MPI_STATUS_SIZE,4) :: status
# endif

# if defined EW_PERIODIC || defined NS_PERIODIC
      integer, parameter :: pp = 1
# else
      integer, parameter :: pp = 0
# endif

      real(r8), dimension((Nghost+pp)*(UBj-LBj+1)*(UBk-LBk+1)) :: AsendW
      real(r8), dimension((Nghost+pp)*(UBj-LBj+1)*(UBk-LBk+1)) :: AsendE
      real(r8), dimension((Nghost+pp)*(UBi-LBi+1)*(UBk-LBk+1)) :: AsendS
      real(r8), dimension((Nghost+pp)*(UBi-LBi+1)*(UBk-LBk+1)) :: AsendN

      real(r8), dimension((Nghost+pp)*(UBj-LBj+1)*(UBk-LBk+1)) :: ArecvW
      real(r8), dimension((Nghost+pp)*(UBj-LBj+1)*(UBk-LBk+1)) :: ArecvE
      real(r8), dimension((Nghost+pp)*(UBi-LBi+1)*(UBk-LBk+1)) :: ArecvS
      real(r8), dimension((Nghost+pp)*(UBi-LBi+1)*(UBk-LBk+1)) :: ArecvN

      character*80 :: Estring

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, 38)
# endif
!
!-----------------------------------------------------------------------
!  Send and recive tile boundary data.
!-----------------------------------------------------------------------
!
!  Set tile partition table for looking up adjacent processes. Notice
!  that a null value is used in places that data transmition is not
!  required.
!
# if defined MPI
      Null_Value=MPI_PROC_NULL
# else
      Null_Value=-1
# endif
      DO j=-1,NtileJ(ng)
        DO i=-1,NtileI(ng)
          table(i,j)=Null_Value
        END DO
      END DO
      rank=0
      DO j=0,NtileJ(ng)-1
        DO i=0,NtileI(ng)-1
          table(i,j)=rank
          IF (MyRank.eq.rank) THEN
            MyRankI=i
            MyRankJ=j
          END IF
          rank=rank+1            
        END DO
      END DO
!
!  Set boundary communication tags and exchange switches.
!
      Wtag=1
      Stag=2
      Etag=3
      Ntag=4
      Wexchange=.false.
      Sexchange=.false.
      Eexchange=.false.
      Nexchange=.false.
!
!-----------------------------------------------------------------------
!  Process Western and Eastern tile boundary data.
!-----------------------------------------------------------------------
!
!  Determine range and length of the distributed tile boundary segments.
!
      Jmin=LBj
      Jmax=UBj
      Jlen=Jmax-Jmin+1
      Klen=UBk-LBk+1
      JKlen=Jlen*Klen
      EWsize=(Nghost+pp)*(UBj-LBj+1)*(UBk-LBk+1)
!
!  Determine the rank of Western and Eastern tiles.  Then, determine
!  the number of ghost-points to send and receive in the West- and
!  East-directions.  The number of ghost-points changes when periodic
!  boundary condition are activated.  The periodicity is as follows:
!
!                      Lm-2  Lm-1  Lm   Lm+1  Lm+2       
!                      -2    -1     0    1     2
!
      IF (EWperiodic) THEN
        IF ((table(MyRankI-1,MyRankJ).eq.Null_Value).and.               &
     &      (NtileI(ng).gt.1)) THEN
          Wtile=table(NtileI(ng)-1,MyRankJ)
          Etile=table(MyRankI+1,MyRankJ)
          GsendW=Nghost
          GsendE=Nghost
          GrecvW=Nghost+1
          GrecvE=Nghost
        ELSE IF ((table(MyRankI+1,MyRankJ).eq.Null_Value).and.          &
     &           (NtileI(ng).gt.1)) THEN
          Wtile=table(MyRankI-1,MyRankJ)
          Etile=table(0,MyRankJ)
          GsendW=Nghost
          GsendE=Nghost+1
          GrecvW=Nghost
          GrecvE=Nghost
        ELSE        
          Wtile=table(MyRankI-1,MyRankJ)
          Etile=table(MyRankI+1,MyRankJ)
          GsendW=Nghost
          GsendE=Nghost
          GrecvW=Nghost
          GrecvE=Nghost
        END IF
      ELSE
        Wtile=table(MyRankI-1,MyRankJ)
        Etile=table(MyRankI+1,MyRankJ)
        GsendW=Nghost
        GsendE=Nghost
        GrecvW=Nghost
        GrecvE=Nghost
      END IF
      IF (Wtile.ne.Null_Value) Wexchange=.true.
      IF (Etile.ne.Null_Value) Eexchange=.true.
!
!  Pack tile boundary data including ghost-points.
!
      JKsizeW=0
      DO m=1,GsendW
        mc=(m-1)*JKlen
        DO k=LBk,UBk
          kc=(k-LBk)*Jlen
          DO j=Jmin,Jmax
            JKsizeW=JKsizeW+1
            jkW=1+(j-Jmin)+kc+mc
            AsendW(jkW)=A(Istr+m-1,j,k)
          END DO
        END DO
      END DO
      JKsizeE=0
      DO m=1,GsendE
        mc=(m-1)*JKlen
        DO k=LBk,UBk
          kc=(k-LBk)*Jlen
          DO j=Jmin,Jmax
            JKsizeE=JKsizeE+1
            jkE=1+(j-Jmin)+kc+mc
            AsendE(jkE)=A(Iend-GsendE+m,j,k)
          END DO
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Send and receive Western and Eastern segments.
!-----------------------------------------------------------------------
!
# if defined MPI
      IF (Wexchange) THEN
        CALL mpi_irecv (ArecvW, EWsize, MP_FLOAT, Wtile, Etag,          &
     &                  OCN_COMM_WORLD, Wrequest, Werror)
      END IF
      IF (Eexchange) THEN
        CALL mpi_irecv (ArecvE, EWsize, MP_FLOAT, Etile, Wtag,          &
     &                  OCN_COMM_WORLD, Erequest, Eerror)
      END IF
      IF (Wexchange) THEN
        CALL mpi_send  (AsendW, JKsizeW, MP_FLOAT, Wtile, Wtag,         &
     &                  OCN_COMM_WORLD, Werror)
      END IF
      IF (Eexchange) THEN
        CALL mpi_send  (AsendE, JKsizeE, MP_FLOAT, Etile, Etag,         &
     &                  OCN_COMM_WORLD, Eerror)
      END IF
      IF (Eexchange) THEN
        CALL mpi_wait (Erequest, status(1,3), Eerror)
        IF (Eerror.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (Eerror, Estring, 80, Eerror)
          print *, Estring
        END IF
      END IF
      IF (Wexchange) THEN
        CALL mpi_wait (Wrequest, status(1,1), Werror)
        IF (Werror.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (Werror, Estring, 80, Werror)
          print *, Estring
        END IF
      END IF
# endif
!
!-----------------------------------------------------------------------
!  Unpack Eastern and Western segments.
!-----------------------------------------------------------------------
!
      IF (Eexchange) THEN
        DO m=1,GrecvE
          mc=(m-1)*JKlen
          DO k=LBk,UBk
            kc=(k-LBk)*Jlen
            DO j=Jmin,Jmax
              jkE=1+(j-Jmin)+kc+mc
              A(Iend+m,j,k)=ArecvE(jkE)
            END DO
          ENDDO
        END DO
      END IF
      IF (Wexchange) THEN
        DO m=GrecvW,1,-1
          mc=(GrecvW-m)*JKlen
          DO k=LBk,UBk
            kc=(k-LBk)*Jlen
            DO j=Jmin,Jmax
              jkW=1+(j-Jmin)+kc+mc
              A(Istr-m,j,k)=ArecvW(jkW)
            END DO
          END DO
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Process Southern and Northern tile boundary data.
!-----------------------------------------------------------------------
!
!  Determine range and length of the distributed tile boundary segments.
!
      Imin=LBi
      Imax=UBi
      Ilen=Imax-Imin+1
      Klen=UBk-LBk+1
      IKlen=Ilen*Klen
      NSsize=(Nghost+pp)*(UBi-LBi+1)*(UBk-LBk+1)
!
!  Determine the rank of Southern and Northern tiles.  Then, determine
!  the number of ghost-points to send and receive in the South- and
!  North-directions.  The number of ghost-points changes when periodic
!  boundary condition are activated.  The periodicity is as follows:
!
!                      Mm-2  Mm-1  Mm   Mm+1  Mm+2       
!                      -2    -1     0    1     2
!
      IF (NSperiodic) THEN
        IF ((table(MyRankI,MyRankJ-1).eq.Null_Value).and.               &
     &      (NtileJ(ng).gt.1)) THEN
          Stile=table(MyRankI,NtileJ(ng)-1)
          Ntile=table(MyRankI,MyRankJ+1)
          GsendS=Nghost
          GsendN=Nghost
          GrecvS=Nghost+1
          GrecvN=Nghost
        ELSE IF ((table(MyRankI,MyRankJ+1).eq.Null_Value).and.          &
     &           (NtileJ(ng).gt.1)) then
          Stile=table(MyRankI,MyRankJ-1)
          Ntile=table(MyRankI,0)
          GsendS=Nghost
          GsendN=Nghost+1
          GrecvS=Nghost
          GrecvN=Nghost
        ELSE
          Stile=table(MyRankI,MyRankJ-1)
          Ntile=table(MyRankI,MyRankJ+1)
          GsendS=Nghost
          GsendN=Nghost
          GrecvS=Nghost
          GrecvN=Nghost
        END IF
      ELSE
        Stile=table(MyRankI,MyRankJ-1)
        Ntile=table(MyRankI,MyRankJ+1)
        GsendS=Nghost
        GsendN=Nghost
        GrecvS=Nghost
        GrecvN=Nghost
      END IF
      IF (Stile.ne.Null_Value) Sexchange=.true.
      IF (Ntile.ne.Null_Value) Nexchange=.true.
!
!  Pack tile boundary data including ghost-points.
!
      IKsizeS=0
      DO m=1,GsendS
        mc=(m-1)*IKlen
        DO k=LBk,UBk
          kc=(k-LBk)*Ilen
          DO i=Imin,Imax
            IKsizeS=IKsizeS+1
            ikS=1+(i-Imin)+kc+mc
            AsendS(ikS)=A(i,Jstr+m-1,k)
          END DO
        END DO
      END DO
      IKsizeN=0
      DO m=1,GsendN
        mc=(m-1)*IKlen
        DO k=LBk,UBk
          kc=(k-LBk)*Ilen
          DO i=Imin,Imax
            IKsizeN=IKsizeN+1
            ikN=1+(i-Imin)+kc+mc
            AsendN(ikN)=A(i,Jend-GsendN+m,k)
          END DO
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Send and receive Southern and Northern segments.
!-----------------------------------------------------------------------
!
# if defined MPI
      IF (Sexchange) THEN
        CALL mpi_irecv (ArecvS, NSsize, MP_FLOAT, Stile, Ntag,          &
     &                  OCN_COMM_WORLD, Srequest, Serror)
      END IF
      IF (Nexchange) THEN
        CALL mpi_irecv (ArecvN, NSsize, MP_FLOAT, Ntile, Stag,          &
     &                  OCN_COMM_WORLD, Nrequest, Nerror)
      END IF
      IF (Sexchange) THEN
        CALL mpi_send  (AsendS, IKsizeS, MP_FLOAT, Stile, Stag,         &
     &                  OCN_COMM_WORLD, Serror)
      END IF
      IF (Nexchange) THEN
        CALL mpi_send  (AsendN, IKsizeN, MP_FLOAT, Ntile, Ntag,         &
     &                  OCN_COMM_WORLD, Nerror)
      END IF
      IF (Nexchange) THEN
        CALL mpi_wait (Nrequest, status(1,4), Nerror)
        IF (Nerror.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (Nerror, Estring, 80, Nerror)
          print *, Estring
        END IF
      END IF
      IF (Sexchange) THEN
        CALL mpi_wait (Srequest, status(1,2), Serror)
        IF (Serror.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (Serror, Estring, 80, Serror)
          print *, Estring
        END IF
      END IF
# endif
!
!-----------------------------------------------------------------------
!  Unpack Northern and Southern segments.
!-----------------------------------------------------------------------
!
      IF (Nexchange) THEN
        DO m=1,GrecvN
          mc=(m-1)*IKlen
          DO k=LBk,UBk
            kc=(k-LBk)*Ilen
            DO i=Imin,Imax
              ikN=1+(i-Imin)+kc+mc
              A(i,Jend+m,k)=ArecvN(ikN)
            END DO
          END DO
        END DO
      END IF
      IF (Sexchange) THEN
        DO m=GrecvS,1,-1
          mc=(GrecvS-m)*IKlen     
          DO k=LBk,UBk
            kc=(k-LBk)*Ilen
            DO i=Imin,Imax
              ikS=1+(i-Imin)+kc+mc
              A(i,Jstr-m,k)=ArecvS(ikS)
            END DO
          END DO
        END DO
      END IF
# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, 38)
# endif

      RETURN
      END SUBROUTINE mp_exchange

      SUBROUTINE mp_gather (ng, LBi, UBi, LBj, UBj, LBk, UBk, gtype,    &
     &                      Ascl,                                       &
# ifdef MASKING
     &                      Amask,                                      &
# endif
     &                      A, Npts, Aout)
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine collects requested real tile array (2D or 3D) from     !
!  each spawned MPI node and stores it into one dimensional global     !
!  array. This routine is used by the  Master Node  to collect and     !
!  pack output data.                                                   !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     LBi        I-dimension Lower bound.                              !
!     UBi        I-dimension Upper bound.                              !
!     LBj        J-dimension Lower bound.                              !
!     UBj        J-dimension Upper bound.                              !
!     LBk        K-dimension Lower bound, if any. Otherwise, a value   !
!                  of one is expected.                                 !
!     LBk        K-dimension Upper bound, if any. Otherwise, a value   !
!                  of one is expected.                                 !
!     UBk        K-dimension Upper bound.                              !
!     gtype      C-grid type. If negative and Land-Sea is available,   !
!                  only water-points processed.                        !
!     Ascl       Factor to scale field before writing.                 !
!     Amask      Land/Sea mask, if any.                                !
!     A          Tile array (2D or 3D) to process.                     !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Npts       Number of points processed in Aout.                   !
!     Aout       Collected data from each node packed into 1D array    !
!                  in column-major order. That is, in the same way     !
!                  that Fortran multi-dimensional arrays are stored    !
!                  in memory.                                          !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_ncparam
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, LBi, UBi, LBj, UBj, LBk, UBk, gtype
      integer, intent(out) :: Npts

      real(r8), intent(in) :: Ascl

# ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in)  :: A(LBi:UBi,LBj:UBj,LBk:UBk)
      real(r8), intent(out) :: Aout((Lm(ng)+2)*(Mm(ng)+2)*(UBk-LBk+1))
!
!  Local variable declarations.
!
      integer :: Itile, Jtile, Nghost
      integer :: Io, Ie, Jo, Je, Ioff, Joff, Koff
      integer :: Ilen, Jlen, Klen, IJlen
      integer :: MyError, MyType, Srequest
      integer :: i, ic, j, jc, k, kc, np, rank

      integer, dimension(0:NtileI(ng)*NtileJ(ng)-1) :: Imin, Imax
      integer, dimension(0:NtileI(ng)*NtileJ(ng)-1) :: Jmin, Jmax
      integer, dimension(0:NtileI(ng)*NtileJ(ng)-1) :: MySize
      integer, dimension(0:NtileI(ng)*NtileJ(ng)-1) :: Rrequest

      integer, dimension(MPI_STATUS_SIZE) :: Rstatus
      integer, dimension(MPI_STATUS_SIZE) :: Sstatus

      real(r8), parameter :: spval = 1.0E+35_r8

      real(r8), dimension((UBi-LBi+1)*(UBj-LBj+1)*(UBk-LBk+1)) :: Asend
      real(r8), dimension((UBi-LBi+1)*(UBj-LBj+1)*(UBk-LBk+1),          &
     &                    NtileI(ng)*NtileJ(ng)-1) :: Arecv

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, 41)
# endif
!
!-----------------------------------------------------------------------
!  Set horizontal starting and ending indices for parallel domain
!  partitions in the XI- and ETA-directions.
!-----------------------------------------------------------------------
!
!  Set first and last grid point according to staggered C-grid
!  classification.
!
      Io=0
      Ie=Lm(ng)+1
      Jo=0
      Je=Mm(ng)+1
      MyType=ABS(gtype)
      IF ((MyType.eq.p2dvar).or.(MyType.eq.u2dvar).or.                  &
     &    (MyType.eq.p3dvar).or.(MyType.eq.u3dvar)) Io=1
      IF ((MyType.eq.p2dvar).or.(MyType.eq.v2dvar).or.                  &
     &    (MyType.eq.p3dvar).or.(MyType.eq.v3dvar)) Jo=1
      IF (Io.eq.0) THEN
        Ioff=1
      ELSE
        Ioff=0
      END IF
      IF (Jo.eq.0) THEN
        Joff=0
      ELSE
        Joff=1
      END IF
      IF (LBk.eq.0) THEN
        Koff=0
      ELSE
        Koff=1
      END IF
      Ilen=Ie-Io+1
      Jlen=Je-Jo+1
      Klen=UBk-LBk+1
      IJlen=Ilen*Jlen
      Npts=IJlen*Klen
!
!  Set physical, non-overlaping (no ghost-points) ranges according to
!  tile rank. Compute size of distributed buffers.
!
      Nghost=0
      DO rank=0,NtileI(ng)*NtileJ(ng)-1
        CALL get_bounds (ng, rank, gtype, Nghost, Itile, Jtile,         &
     &                   Imin(rank), Imax(rank),                        &
     &                   Jmin(rank), Jmax(rank))
        MySize(rank)=(Imax(rank)-Imin(rank)+1)*                         &
     &               (Jmax(rank)-Jmin(rank)+1)*(UBk-LBk+1)
      END DO
!
!-----------------------------------------------------------------------
!  Collect requested array data. 
!-----------------------------------------------------------------------
!
!  Pack and scale input data.
!
      np=0
      DO k=LBk,UBk
        DO j=Jmin(MyRank),Jmax(MyRank)
          DO i=Imin(MyRank),Imax(MyRank)
            np=np+1
            Asend(np)=A(i,j,k)*Ascl
          END DO
        END DO
      END DO
# ifdef MASKING
!
!  If processing water-points only, flag land-points.
!
      IF (gtype.lt.0) THEN
        np=0
        DO k=LBk,UBk
          DO j=Jmin(MyRank),Jmax(MyRank)
            DO i=Imin(MyRank),Imax(MyRank)
              np=np+1
              IF (Amask(i,j).eq.0.0_r8) THEN
                Asend(np)=spval
              END IF
            END DO
          END DO
        END DO
      END IF
# endif
!
!  If master processor, unpack the send buffer since there is not
!  need to distribute.
!
      IF (MyRank.eq.MyMaster) THEN
        np=0
        DO k=LBk,UBk
          kc=(k-Koff)*IJlen
          DO j=Jmin(MyRank),Jmax(MyRank)
            jc=(j-Joff)*Ilen+kc
            DO i=Imin(MyRank),Imax(MyRank)
              np=np+1
              ic=i+Ioff+jc
              Aout(ic)=Asend(np)
            END DO
          END DO
        END DO
      END IF
!
!  Send, receive, and unpacked data.
!
      IF (MyRank.eq.MyMaster) THEN
        DO rank=1,NtileI(ng)*NtileJ(ng)-1
          CALL mpi_irecv (Arecv(1,rank), MySize(rank), MP_FLOAT, rank,  &
     &                    rank+5, OCN_COMM_WORLD, Rrequest(rank),       &
     &                    MyError)
        END DO
        DO rank=1,NtileI(ng)*NtileJ(ng)-1
          CALL mpi_wait (Rrequest(rank), Rstatus, MyError)
          np=0
          DO k=LBk,UBk
            kc=(k-Koff)*IJlen
            DO j=Jmin(rank),Jmax(rank)
              jc=(j-Joff)*Ilen+kc
              DO i=Imin(rank),Imax(rank)
                np=np+1
                ic=i+Ioff+jc
                Aout(ic)=Arecv(np,rank)
              END DO
            END DO
          END DO
        END DO
      ELSE
        CALL mpi_isend (Asend, MySize(MyRank), MP_FLOAT, MyMaster,      &
     &                  MyRank+5, OCN_COMM_WORLD, Srequest, MyError)
        CALL mpi_wait (Srequest, Sstatus, MyError)
      END IF

# ifdef MASKING
!
! If pocessing only water-points, remove land points and repack.
!
      IF ((MyRank.eq.MyMaster).and.(gtype.lt.0)) THEN
        ic=0
        np=IJlen*Klen
        DO i=1,np
          IF (Aout(i).lt.spval) THEN
            ic=ic+1
            Aout(ic)=Aout(i)
          END IF
        END DO
        Npts=ic
      END IF
# endif
# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, 41)
# endif

      RETURN
      END SUBROUTINE mp_gather

      SUBROUTINE mp_reduce (ng, Asize, A, op_handle)
!
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!=======================================================================
!                                                                      !
!  This routine collects and reduces requested variables from all      !
!  nodes in the group.  Then,  it broadcasts reduced variables to      !
!  all nodes in the group.                                             !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     Asize      Number of scalar variables to reduce.                 !
!     A          Vector of scalar variables to reduce.                 !
!     op_handle  Reduction operation handle (string).  The following   !
!                  reduction operations are supported:                 !
!                  'MIN', 'MAX', 'SUM'                                 !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A          Vector of reduced scalar variables.                   !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Asize

      character (len=*), intent(in) :: op_handle(Asize)

      real(r8), intent(inout) :: A(Asize)
!
!  Local variable declarations.
!
      integer :: handle, i, rank
      integer :: Rerror, Srequest, Serror

      integer, dimension(0:NtileI(ng)*NtileJ(ng)-1) :: Rrequest

      integer, dimension(MPI_STATUS_SIZE) :: Rstatus
      integer, dimension(MPI_STATUS_SIZE) :: Sstatus

      real(r8), dimension(Asize,0:NtileI(ng)*NtileJ(ng)-1) :: Arecv
      real(r8), dimension(Asize) :: Areduce
      real(r8), dimension(Asize) :: Asend

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, 40)
# endif
!
!-----------------------------------------------------------------------
!  Collect and reduce requested scalar variables.
!-----------------------------------------------------------------------
!
!  Pack data to reduce.
!
      DO i=1,Asize
        Asend(i)=A(i)
      END DO
!
!  Collect and reduce.
!
# if defined ALLREDUCE
      DO i=1,Asize
        IF (op_handle(i)(1:3).eq.'MIN') THEN
          handle=MPI_MIN
        ELSE IF (op_handle(i)(1:3).eq.'MAX') THEN
          handle=MPI_MAX
        ELSE IF (op_handle(i)(1:3).eq.'SUM') THEN
          handle=MPI_SUM
        END IF
        CALL mpi_allreduce (Asend(i), Areduce(i), 1, MP_FLOAT, handle,  &
     &                      OCN_COMM_WORLD, Rerror)
      END DO
# elif defined ALLGATHER
      CALL mpi_allgather (Asend, Asize, MP_FLOAT,                       &
     &                    Arecv, Asize, MP_FLOAT,                       &
     &                    OCN_COMM_WORLD, Rerror)
      DO i=1,Asize
        Areduce(i)=Arecv(i,0)
        DO rank=1,NtileI(ng)*NtileJ(ng)-1
          IF (op_handle(i)(1:3).eq.'MIN') THEN
            Areduce(i)=MIN(Areduce(i),Arecv(i,rank))
          ELSE IF (op_handle(i)(1:3).eq.'MAX') THEN
            Areduce(i)=MAX(Areduce(i),Arecv(i,rank))
          ELSE IF (op_handle(i)(1:3).eq.'SUM') THEN
            Areduce(i)=Areduce(i)+Arecv(i,rank)
          END IF
        END DO
      END DO
# else
      IF (MyRank.eq.MyMaster) THEN
        DO rank=1,NtileI(ng)*NtileJ(ng)-1
          CALL mpi_irecv (Arecv(1,rank), Asize, MP_FLOAT, rank,         &
     &                    rank+500, OCN_COMM_WORLD, Rrequest(rank),     &
     &                    Rerror)
        END DO
        DO i=1,Asize
          Areduce(i)=Asend(i)
        END DO
        DO rank=1,NtileI(ng)*NtileJ(ng)-1
          CALL mpi_wait (Rrequest(rank), Rstatus, Rerror)
          DO i=1,Asize
            IF (op_handle(i)(1:3).eq.'MIN') THEN
              Areduce(i)=MIN(Areduce(i),Arecv(i,rank))
            ELSE IF (op_handle(i)(1:3).eq.'MAX') THEN
              Areduce(i)=MAX(Areduce(i),Arecv(i,rank))
            ELSE IF (op_handle(i)(1:3).eq.'SUM') THEN
              Areduce(i)=Areduce(i)+Arecv(i,rank)
            END IF
          END DO
        END DO
      ELSE
        CALL mpi_isend (Asend, Asize, MP_FLOAT, MyMaster, MyRank+500,   &
     &                  OCN_COMM_WORLD, Srequest, Serror)
        CALL mpi_wait (Srequest, Sstatus, Serror)
      END IF
!
!  Broadcast reduced variables from process to all processes in the
!  group.
!
      CALL mpi_bcast (Areduce, Asize, MP_FLOAT, MyMaster,               &
     &                OCN_COMM_WORLD, Serror)
# endif
!
!  Unpack.
!
      DO i=1,Asize
        A(i)=Areduce(i)
      END DO
# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, 40)
# endif

      RETURN
      END SUBROUTINE mp_reduce

      SUBROUTINE mp_scatter (ng, LBi, UBi, LBj, UBj, LBk, UBk,          &
     &                       Nghost, gtype, Amin, Amax,                 &
# if defined READ_WATER && defined MASKING
     &                       NWpts, IJ_water,                           &
# endif
     &                       Npts, A, Aout)
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine broadcasts input global data, packed as 1D real array, !
!  to each spawned MPI node.  The output tile array (2D or 3D) has its !
!  ghost-points updated in the halo region.  It is used by the master  !
!  node to scatter input global data to each tiled node.               !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     LBi        I-dimension Lower bound.                              !
!     UBi        I-dimension Upper bound.                              !
!     LBj        J-dimension Lower bound.                              !
!     UBj        J-dimension Upper bound.                              !
!     LBk        K-dimension Lower bound, if any. Otherwise, a value   !
!                  of one is expected.                                 !
!     LBk        K-dimension Upper bound, if any. Otherwise, a value   !
!                  of one is expected.                                 !
!     UBk        K-dimension Upper bound.                              !
!     Nghost     Number of ghost-points in the halo region.            !
!     gtype      C-grid type. If negative and Land-Sea mask is         !
!                  available, only water-points are processed.         !
!     Amin       Input array minimum value.                            !
!     Amax       Input array maximum value.                            !
!     NWpts      Number of water points.                               !
!     IJ_water   IJ-indices for water points.                          !
!     Npts       Number of points to processes in A.                   !
!     A          Input global data from each node packed into 1D array !
!                  in column-major order. That is, in the same way     !
!                  that Fortran multi-dimensional arrays are stored    !
!                  in memory.                                          !
!     Npts       Number of points to processes in A.                   !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Aout       Tile array (2D or 3D) with updated ghost-points in    !
!                  the halo region.                                    !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits 
      USE mod_ncparam
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, LBi, UBi, LBj, UBj, LBk, UBk
      integer, intent(in) :: Nghost, gtype, Npts

# if defined READ_WATER && defined MASKING
      integer, intent(in) :: NWpts
      integer, intent(in) :: IJ_water(NWpts)
# endif
      real(r8), intent(inout) :: Amin, Amax
      real(r8), intent(inout) :: A(Npts+2)
      real(r8), intent(out) :: Aout(LBi:UBi,LBj:UBj,LBk:UBk)
!
!  Local variable declarations.
!
      integer :: Io, Ie, Jo, Je, Ioff, Joff, Koff
      integer :: Ilen, Jlen, Klen, IJlen
      integer :: MyError, MySize, MyType, rank
      integer :: i, ic, ij, j, jc, k, kc, mc, nc

      integer, dimension(0:NtileI(ng)*NtileJ(ng)-1) :: Imin, Imax
      integer, dimension(0:NtileI(ng)*NtileJ(ng)-1) :: Jmin, Jmax
      integer, dimension(0:NtileI(ng)*NtileJ(ng)-1) :: Itile, Jtile

      real(r8), dimension((Lm(ng)+2)*(Mm(ng)+2)*(UBk-LBk+1)) :: Arecv

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, 42)
# endif
!
!-----------------------------------------------------------------------
!  Set horizontal starting and ending indices for parallel domain
!  partitions in the XI- and ETA-directions.
!-----------------------------------------------------------------------
!
!  Set first and last grid point according to staggered C-grid
!  classification. The, set 1D counter offsets.
!
      Io=0
      Ie=Lm(ng)+1
      Jo=0
      Je=Mm(ng)+1
      MyType=ABS(gtype)
      IF ((MyType.eq.p2dvar).or.(MyType.eq.u2dvar).or.                  &
     &    (MyType.eq.p3dvar).or.(MyType.eq.u3dvar)) Io=1
      IF ((MyType.eq.p2dvar).or.(MyType.eq.v2dvar).or.                  &
     &    (MyType.eq.p3dvar).or.(MyType.eq.v3dvar)) Jo=1
      IF (Io.eq.0) THEN
        Ioff=1
      ELSE
        Ioff=0
      END IF
      IF (Jo.eq.0) THEN
        Joff=0
      ELSE
        Joff=1
      END IF
      IF (LBk.eq.0) THEN
        Koff=0
      ELSE
        Koff=1
      END IF
      Ilen=Ie-Io+1
      Jlen=Je-Jo+1
      Klen=UBk-LBk+1
      IJlen=Ilen*Jlen
!
!  Set physical, non-overlaping (no ghost-points) ranges according to
!  tile rank.
!
      DO rank=0,NtileI(ng)*NtileJ(ng)-1
        CALL get_bounds (ng, rank, gtype, Nghost,                       &
     &                   Itile(rank), Jtile(rank),                      &
     &                   Imin(rank), Imax(rank),                        &
     &                   Jmin(rank), Jmax(rank))
      END DO
!
!  Size of broadcast buffer.
!
      IF (gtype.gt.0) THEN
        MySize=IJlen*Klen
      ELSE
        MySize=Npts
      END IF
!
!-----------------------------------------------------------------------
!  Scatter requested array data. 
!-----------------------------------------------------------------------
!
!  If master processor, append minimum and maximum values to the end of
!  the buffer. 
!
      IF (MyRank.eq.MyMaster) Then
        A(MySize+1)=Amin
        A(MySize+2)=Amax
      END IF
      MySize=MySize+2
!
!  Broadcast data to all processes in the group, itself included.
!
      CALL mpi_bcast (A, MySize, MP_FLOAT, MyMaster, OCN_COMM_WORLD,    &
     &                MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        WRITE (stdout,10) 'BCAST', MyRank, MyError
 10     FORMAT (/,' MP_SCATTER - error during ',a,' call, Node = ',     &
     &          i3.3, ' Error = ',i3)
        exit_flag=2
        RETURN
      END IF
!
!  If water points only, fill land points.
!
      IF (gtype.gt.0) THEN
        DO nc=1,MySize-2
          Arecv(nc)=A(nc)
        END DO
# if defined READ_WATER && defined MASKING
      ELSE
        nc=0
        DO k=LBk,UBk
          kc=(k-Koff)*IJlen
          ij=0
          mc=0
          DO j=Jo,Je
            jc=(j-Joff)*Ilen+kc
            DO i=Io,Ie
              ij=ij+1
              ic=i+Ioff+jc
              IF (IJ_water(mc+1).eq.ij) THEN
                mc=mc+1
                nc=nc+1
                Arecv(ic)=A(nc)
              ELSE
                Arecv(ic)=0.0_r8
              ENDIF
            END DO
          END DO
        END DO
# endif
      END IF
!
!  Unpack data buffer.
!
      DO k=LBk,UBk
        kc=(k-Koff)*IJlen
        DO j=Jmin(MyRank),Jmax(MyRank)
          jc=(j-Joff)*Ilen+kc
          DO i=Imin(MyRank),Imax(MyRank)
            ic=i+Ioff+jc
            Aout(i,j,k)=Arecv(ic)
          END DO
        END DO
      END DO
      Amin=A(MySize-1)
      Amax=A(MySize)
      
# if defined EW_PERIODIC || defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Apply periodic boundary conditions.
!-----------------------------------------------------------------------
!
#  ifdef EW_PERIODIC
      IF (Itile(MyRank).eq.(NtileI(ng)-1)) THEN
        DO k=LBk,UBk
          kc=(k-Koff)*IJlen
          DO j=Jmin(MyRank),Jmax(MyRank)
            ic=Ioff+(j-Joff)*Ilen+kc
            Aout(Lm(ng)+1,j,k)=Arecv(1+ic)
            Aout(Lm(ng)+2,j,k)=Arecv(2+ic)
          END DO
        END DO
      END IF
      IF (Itile(MyRank).eq.0) THEN
        DO k=LBk,UBk
          kc=(k-Koff)*IJlen
          DO j=Jmin(MyRank),Jmax(MyRank)
            ic=Ioff+(j-Joff)*Ilen+kc
            Aout(-2,j,k)=Arecv(Lm(ng)-2+ic)
            Aout(-1,j,k)=Arecv(Lm(ng)-1+ic)
            Aout( 0,j,k)=Arecv(Lm(ng)  +ic)
          END DO
        END DO
      END IF
#  endif
#  ifdef NS_PERIODIC
      IF (Jtile(MyRank).eq.(NtileJ(ng)-1)) THEN
        DO k=LBk,UBk
          kc=(k-Koff)*IJlen
          DO i=Imin(MyRank),Imax(MyRank)
            Aout(i,Mm(ng)+1,k)=Arecv(i+Ioff+(1-Joff)*Ilen+kc)
            Aout(i,Mm(ng)+2,k)=Arecv(i+Ioff+(2-Joff)*Ilen+kc)
          END DO
        END DO
      END IF
      IF (Jtile(MyRank).eq.0) THEN
        DO k=LBk,UBk
          kc=(k-Koff)*IJlen
          DO i=Imin(MyRank),Imax(MyRank)
            Aout(i,-2,k)=Arecv(i+Ioff+(Mm(ng)-2-Joff)*Ilen+kc)
            Aout(i,-1,k)=Arecv(i+Ioff+(Mm(ng)-1-Joff)*Ilen+kc)
            Aout(i, 0,k)=Arecv(i+Ioff+(Mm(ng)  -Joff)*Ilen+kc)
          END DO
        END DO
      END IF
#  endif
#  if defined EW_PERIODIC && defined NS_PERIODIC
      IF ((Itile(MyRank).eq.(NtileI(ng)-1)).and.                        &
     &    (Jtile(MyRank).eq.(NtileJ(ng)-1))) THEN
        DO k=LBk,UBk
          kc=(k-Koff)*IJlen
          DO j=1,2
            jc=Ioff+(j-Joff)*Ilen+kc
            Aout(Lm(ng)+1,Mm(ng)+j,k)=Arecv(1+jc)
            Aout(Lm(ng)+2,Mm(ng)+j,k)=Arecv(2+jc)
          END DO
        END DO
      END IF
      IF ((Itile(MyRank).eq.0).and.                                     &
     &    (Jtile(MyRank).eq.(NtileJ(ng)-1))) THEN
        DO k=LBk,UBk
          kc=(k-Koff)*IJlen
          DO j=1,2
            jc=Ioff+(j-Joff)*Ilen+kc
            Aout(-2,Mm(ng)+j,k)=Arecv(Lm(ng)-2+jc)
            Aout(-1,Mm(ng)+j,k)=Arecv(Lm(ng)-1+jc)
            Aout( 0,Mm(ng)+j,k)=Arecv(Lm(ng)  +jc)
          END DO
        END DO
      END IF
      IF ((Itile(MyRank).eq.(NtileI(ng)-1)).and.                        &
     &    (Jtile(MyRank).eq.0)) THEN
        DO k=LBk,UBk
          kc=(k-Koff)*IJlen
          DO j=-2,0
            jc=Ioff+(Mm(ng)+j-Joff)*Ilen+kc
            Aout(Lm(ng)+1,j,k)=Arecv(1+jc)
            Aout(Lm(ng)+2,j,k)=Arecv(2+jc)
          END DO
        END DO
      END IF
      IF ((Itile(MyRank).eq.0).and.                                     &
     &    (Jtile(MyRank).eq.0)) THEN
        DO k=LBk,UBk
          kc=(k-Koff)*IJlen
          DO j=-2,0
            jc=Ioff+(Mm(ng)+j-Joff)*Ilen+kc
            Aout(-2,j,k)=Arecv(Lm(ng)-2+jc)
            Aout(-1,j,k)=Arecv(Lm(ng)-1+jc)
            Aout( 0,j,k)=Arecv(Lm(ng)  +jc)
          END DO
        END DO
      END IF
#  endif
# endif
# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, 42)
# endif

      RETURN
      END SUBROUTINE mp_scatter

      SUBROUTINE mp_dump (ng, tile, gtype,                              &
     &                    ILB, IUB, JLB, JUB, KLB, KUB, A, name)
!
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!=======================================================================
!                                                                      !
!  This routine is used to debug distributed-memory communications.    !
!  It writes field into an ASCII file for further post-processing.     !
!                                                                      !
!=======================================================================
!

      USE mod_param
      USE mod_parallel
      USE mod_ncparam

      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, gtype
      integer, intent(in) :: ILB, IUB, JLB, JUB, KLB, KUB

      real(r8), intent(in) :: A(ILB:IUB,JLB:JUB,KLB:KUB)

      character (len=*) :: name
!
!  Local variable declarations.
!
      common /counter/ nc
      integer :: nc

      logical, save :: first = .true.

      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: Imin, Imax, Ioff, Jmin, Jmax, Joff
      integer :: unit

#include "tile.h"
#include "set_bounds.h"
!
!------------------------------------------------------------------------
!  Write out requested field.
!------------------------------------------------------------------------
!
      IF (first) THEN
        nc=0
        first=.false.
      END IF
      nc=nc+1
      IF (Master) THEN
        WRITE (10,'(a,i3.3,a,a)') 'file ', nc, ': ', TRIM(name)
        CALL my_flush (10)
      END IF
!
!  Write out field including ghost-points.
!
      Imin=0
      Imax=Lm(ng)+1
#ifdef EW_PERIODIC
      Ioff=3
#else
      Ioff=1
#endif
      Jmin=0
      Jmax=Mm(ng)+1
#ifdef NS_PERIODIC
      Joff=3
#else
      Joff=1
#endif
      IF ((gtype.eq.p2dvar).or.(gtype.eq.p3dvar).or.                    &
     &    (gtype.eq.u2dvar).or.(gtype.eq.u3dvar)) THEN
        Imin=1
      END IF
      IF ((gtype.eq.p2dvar).or.(gtype.eq.p3dvar).or.                    &
     &    (gtype.eq.v2dvar).or.(gtype.eq.v3dvar)) THEN
        Jmin=1
      END IF

      unit=(MyRank+1)*1000+nc
      WRITE (unit,*) ILB, IUB, JLB, JUB, KLB, KUB,                      &
     &               Ioff, Joff, Imin, Imax, Jmin, Jmax,                &
     &               A(ILB:IUB,JLB:JUB,KLB:KUB)
      CALL my_flush (unit)
!
!  Write out non-overlapping field.
!
      Imin=IstrR
      Imax=IendR
#ifdef EW_PERIODIC
      Ioff=2
#else
      Ioff=1
#endif
      Jmin=JstrR
      Jmax=JendR
#ifdef NS_PERIODIC
      Joff=2
#else
      Joff=1
#endif
      IF ((gtype.eq.p2dvar).or.(gtype.eq.p3dvar).or.                    &
     &    (gtype.eq.u2dvar).or.(gtype.eq.u3dvar)) THEN
        Imin=Istr
        Ioff=Ioff-1
      END IF
      IF ((gtype.eq.p2dvar).or.(gtype.eq.p3dvar).or.                    &
     &    (gtype.eq.v2dvar).or.(gtype.eq.v3dvar)) THEN
        Jmin=Jstr
        Joff=Joff-1
      END IF

      unit=(MyRank+1)*10000+nc
      WRITE (unit,*) Imin, Imax, Jmin, Jmax, KLB, KUB,                  &
     &               Ioff, Joff, Imin, Imax, Jmin, Jmax,                &
     &               A(Imin:Imax,Jmin:Jmax,KLB:KUB)
      CALL my_flush (unit)

      RETURN
      END SUBROUTINE mp_dump
#else
      SUBROUTINE distribute
      END SUBROUTINE distribute
#endif
