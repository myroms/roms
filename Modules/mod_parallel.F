#include "cppdefs.h"
      MODULE mod_parallel
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This module contains all variables used for parallelization         !
!                                                                      !
!=======================================================================
!
        USE mod_param

        implicit none

#ifdef MPI
        include 'mpif.h'
#endif
!
!  Switch to identify master processor. In serial and shared-memory
!  applications it is always true.
!
        logical :: Master
!
!  Switch to identify which thread is processing input/output files.
!  In distributed-memory applications, this thread can be the master
!  thread or all threads in case of parallel output. In serial and
!  shared-memory applications it is always true.
!
        logical :: InpThread
        logical :: OutThread
!
!  Number of shared-memory parallel threads.  In distributed memory
!  configurations, its value must be equal to one.
!
        integer :: numthreads = 1
!
!  Number distributed memory nodes.
!
        integer :: numnodes = 0

#ifdef AIR_OCEAN
!
!  Parallel nodes assined to atmosphere and ocean models.
!
        integer :: peATM_frst          ! first atmosphere parallel node
        integer :: peATM_last          ! last  atmosphere parallel node 
        integer :: peOCN_frst          ! first ocean parallel node
        integer :: peOCN_last          ! last  ocean parallel node
#endif
!
!  Parallel threads/nodes counters used in critical parallel regions.
!
        integer :: tile_count = 0
        integer :: block_count = 0
        integer :: thread_count = 0
!
!  Profiling variables as function of parallel thread:
!
!    proc          Parallel process ID.
!    Cstr          Starting time for program region.
!    Cend          Ending time for program region.
!    Csum          Accumulated time for progam region.
!
        common /process/ proc
        integer  :: proc(0:1,Ngrids)
!$OMP THREADPRIVATE (/process/)

        common /wallclock/ Cstr, Cend
        real(r8) :: Cstr(0:44,Ngrids)
        real(r8) :: Cend(0:44,Ngrids)
!$OMP THREADPRIVATE (/wallclock/)

        real(r8) :: Csum(0:44,Ngrids)

#ifdef DISTRIBUTE
# ifdef MPI
!
!  Ocean model MPI group communicator handle.
!
        integer :: OCN_COMM_WORLD
# endif
!
!  Distributed-memory master process and rank of the local process.
!
        integer, parameter :: MyMaster = 0
        integer :: MyRank
!
!  Type of message-passage floating point bindings.
!
# ifdef DOUBLE_PRECISION
#  ifdef MPI
!!      integer, parameter :: MP_FLOAT = MPI_DOUBLE_PRECISION
        integer, parameter :: MP_FLOAT = MPI_REAL8
#  endif
# else
#  ifdef MPI
        integer, parameter :: MP_FLOAT = MPI_REAL
!!      integer, parameter :: MP_FLOAT = MPI_REAL4
#  endif
# endif
#endif
        CONTAINS

        SUBROUTINE initialize_parallel
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine initializes and spawn distribute-memory nodes.         !
!                                                                      !
!=======================================================================
!
          USE mod_param
          USE mod_iounits
          USE mod_scalars
!
!  Local variable declarations.
!
          integer :: MyError, ng, nsize
#ifndef DISTRIBUTE
          integer :: my_numthreads
!
!-----------------------------------------------------------------------
!  Initialize shared-memory (OpenMP) or serial configuration.
!-----------------------------------------------------------------------
!
!  Inquire number of threads in parallel region.
!
!$OMP PARALLEL SHARED(numthreads)
          numthreads=my_numthreads()
!$OMP END PARALLEL
          Master=.true.
          InpThread=.true.
          OutThread=.true.
#endif
#ifdef DISTRIBUTE
# ifdef MPI
!
!-----------------------------------------------------------------------
!  Initialize distributed-memory (MPI) configuration.
!-----------------------------------------------------------------------
!
!  Get the number of processes in the group associated with the world
!  communicator.
!
          numthreads=1
          CALL mpi_comm_size (OCN_COMM_WORLD, numnodes, MyError)
          IF (MyError.ne.0) THEN
            WRITE (stdout,10)
  10        FORMAT (/,' ROMS/TOMS - Unable to inquire number of',       &
     &              ' processors in the group.')
            exit_flag=6
            RETURN
          END IF
!
!  Identify master thread and output thread.
!
          Master=.false.
          InpThread=.false.
          OutThread=.false.
          IF (MyRank.eq.MyMaster) THEN
            Master=.true.
            InpThread=.true.
            OutThread=.true.
          END IF
# endif
#endif
          RETURN
        END SUBROUTINE initialize_parallel
      END MODULE mod_parallel
