#include "macros.h"

program concurrent_coupled

! TBD.  Get rid of this when the MCT implementation becomes
! implicit.
! use MCT_WRF_IO_MODULE

  use ASSERT_MODULE, only : Assert

  use MPH_module, only : MPH_COMPONENTS, &   ! MPH multi-component setup.
    MPH_DEBUG, &
    MPH_COMM_JOIN, &
    MPH_COMP_ID, &
    MPH_COMP_NAME, &
    PE_in_component, &
    MPH_LOCAL_WORLD

  use mod_parallel, only : Master
  USE mod_ocean, only : initialize_ocean => initialize
  USE mod_ocean, only : run_ocean => run
  USE mod_ocean, only : finalize_ocean => finalize

#ifdef PROFILE
  USE wall_clocks_module, only : CLOCKS
  USE wall_clocks_module, only : CLOCK_ON
  USE wall_clocks_module, only : CLOCK_OFF
  USE wall_clocks_module, only : CONSTRUCT_CLOCKS => CONSTRUCTOR
  USE wall_clocks_module, only : CLOCK_TIME
  USE wall_clocks_module, only : T_CLOCK_CLASS => T_CLASS
  USE wall_clocks_module, only : GET_MAXIMUMS
#endif

  implicit none
  include 'mpif.h'
  integer :: STATUS
  integer :: EXE_COMM
  integer :: LOCAL_COMM
  integer :: MyRank

  character(LEN=80) :: ID_STR
  character(LEN=80) :: NAME_STR
  character(LEN=132) :: SYS_STR

  real, parameter :: SCALE = 1.0

  integer :: RETURN_VALUE

  integer :: NUM_WRF_STEPS
  integer :: NUM_ROMS_STEPS

#ifdef PROFILE
  type (T_CLOCK_CLASS) :: CLOCK_MAXIMUMS
#endif

  namelist /COUPLE/ NUM_WRF_STEPS, NUM_ROMS_STEPS

  call MPI_INIT(STATUS)

  CALL wrf_termio_dup

  call MPH_debug (0)

  EXE_COMM = MPH_COMPONENTS(name1="wrf", name2="roms")

  open(10, file='couple.nl', form='formatted')
  read(10, COUPLE)
  close(10)

#ifdef PROFILE
  call CONSTRUCT_CLOCKS(CLOCKS)
#endif

  if (PE_in_component("wrf", LOCAL_COMM)) then
    call MPI_COMM_RANK(LOCAL_COMM, MyRank, STATUS)
    call WRF_INIT(LOCAL_COMM)

#ifdef PROFILE
    call MPI_BARRIER(MPI_COMM_WORLD, STATUS)
    call CLOCK_ON(111, CLOCKS)
#endif
    call WRF_RUN(NUM_WRF_STEPS, RETURN_VALUE)

#ifdef PROFILE
    call CLOCK_OFF(111, CLOCKS)
#endif

    ASSERT(RETURN_VALUE == 1)
#ifdef PROFILE
    call GET_MAXIMUMS(CLOCKS, LOCAL_COMM, CLOCK_MAXIMUMS)
    if (MyRank == 0) then
      print *, 'WRF  main loop    TIMINGS  ', CLOCK_TIME(111, CLOCK_MAXIMUMS)
      print *, '  WRF  read       TIMINGS  ', CLOCK_TIME(100, CLOCK_MAXIMUMS)
      print *, '  WRF  write      TIMINGS  ', CLOCK_TIME(101, CLOCK_MAXIMUMS)
      print *, '  WRF  Skew       TIMINGS  ', CLOCK_TIME(120, CLOCK_MAXIMUMS)
      print *, '  WRF  recv/wait  TIMINGS  ', CLOCK_TIME(8, CLOCK_MAXIMUMS)
      call flush(6)
    end if
#endif
    call wrf_finalize

  else if (PE_in_component("roms", LOCAL_COMM)) then
    call initialize_ocean(LOCAL_COMM)

#ifdef PROFILE
    call MPI_BARRIER(MPI_COMM_WORLD, STATUS)
    call CLOCK_ON(112, CLOCKS)
#endif
    call run_ocean(NUM_ROMS_STEPS, RETURN_VALUE)

#ifdef PROFILE
    call CLOCK_OFF(112, CLOCKS)
#endif

    ASSERT(RETURN_VALUE == 1)
#ifdef PROFILE
    call GET_MAXIMUMS(CLOCKS, LOCAL_COMM, CLOCK_MAXIMUMS)
    if (Master) then
      print *, 'ROMS main loop    TIMINGS  ', CLOCK_TIME(112, CLOCK_MAXIMUMS)
      print *, '  ROMS read       TIMINGS  ', CLOCK_TIME(102, CLOCK_MAXIMUMS)
      print *, '  ROMS write      TIMINGS  ', CLOCK_TIME(103, CLOCK_MAXIMUMS)
      print *, '  ROMS Skew       TIMINGS  ', CLOCK_TIME(120, CLOCK_MAXIMUMS)
      print *, '  ROMS recv/wait  TIMINGS  ', CLOCK_TIME(8, CLOCK_MAXIMUMS)
      call flush(6)
    end if
#endif
    call finalize_ocean

    call MPI_FINALIZE(STATUS)

  else
    ASSERT(.false.)
  end if

  call EXT_MCT_IOEXIT ( &
    STATUS )
  ASSERT(STATUS == 0)

end program
