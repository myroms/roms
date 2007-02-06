#include "macros.h"

program sequential_coupled

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
  USE wall_clocks_module, only : GET_MAXIMUMS
  USE wall_clocks_module, only : T_CLOCK_CLASS => T_CLASS
#endif

  implicit none
  integer :: STATUS
  integer :: EXE_COMM
  integer :: MyRank

  character(LEN=80) :: ID_STR
  character(LEN=80) :: NAME_STR
  character(LEN=132) :: SYS_STR
  
  real, parameter :: SCALE = 1.0

  integer :: RETURN_VALUE

  include 'mpif.h'

#ifdef PROFILE
  type (T_CLOCK_CLASS) :: CLOCK_MAXIMUMS
#endif

  integer :: NUM_WRF_STEPS
  integer :: NUM_ROMS_STEPS

  namelist /COUPLE/ NUM_WRF_STEPS, NUM_ROMS_STEPS

  call MPI_INIT(STATUS)

  CALL wrf_termio_dup

  call MPH_debug (0)

  EXE_COMM = MPH_COMPONENTS(name1="wrf", name2="roms")

  open(10, file='couple.nl', form='formatted')
  read(10, COUPLE)
  close(10)

  SYS_STR = "COMPONENT1=wrf, COMPONENT2=roms"

  call EXT_MCT_IOINIT ( SYS_STR, STATUS )

#ifdef PROFILE
  call CONSTRUCT_CLOCKS(CLOCKS)
#endif

  call WRF_INIT(MPI_COMM_WORLD)
  call initialize_ocean(MPI_COMM_WORLD)
  call MPI_COMM_RANK(MPI_COMM_WORLD, MyRank, STATUS)

  RETURN_VALUE = 0

#ifdef PROFILE
  call CLOCK_ON(110, CLOCKS)
#endif

  do while (RETURN_VALUE == 0)

#ifdef PROFILE
  call CLOCK_ON(111, CLOCKS)
#endif

    call WRF_RUN(NUM_WRF_STEPS, RETURN_VALUE)

#ifdef PROFILE
  call CLOCK_OFF(111, CLOCKS)

  call GET_MAXIMUMS(CLOCKS, MPI_COMM_WORLD, CLOCK_MAXIMUMS)
#endif

#ifdef PROFILE
  call CLOCK_ON(112, CLOCKS)
#endif

    call run_ocean(NUM_ROMS_STEPS, RETURN_VALUE)

#ifdef PROFILE
  call CLOCK_OFF(112, CLOCKS)

  call GET_MAXIMUMS(CLOCKS, MPI_COMM_WORLD, CLOCK_MAXIMUMS)
#endif

  end do

#ifdef PROFILE
  call CLOCK_OFF(110, CLOCKS)
  call GET_MAXIMUMS(CLOCKS, MPI_COMM_WORLD, CLOCK_MAXIMUMS)
  if (MyRank == 0) then
    print *, 'Main model loop time ', CLOCK_TIME(110, CLOCK_MAXIMUMS)
    print *, '  WRF  model loop time ', CLOCK_TIME(111, CLOCK_MAXIMUMS)
    print *, '  ROMS model loop time ', CLOCK_TIME(112, CLOCK_MAXIMUMS)
    print *, 'WRF  Read  data time ', CLOCK_TIME(100, CLOCK_MAXIMUMS)
    print *, 'WRF  Write data time ', CLOCK_TIME(101, CLOCK_MAXIMUMS)
    print *, 'ROMS Read  data time ', CLOCK_TIME(102, CLOCK_MAXIMUMS)
    print *, 'ROMS Write data time ', CLOCK_TIME(103, CLOCK_MAXIMUMS)
  end if
#endif

  call EXT_MCT_IOEXIT ( &
    STATUS )
  ASSERT(STATUS == 0)

  call wrf_finalize

  call finalize_ocean

end program

