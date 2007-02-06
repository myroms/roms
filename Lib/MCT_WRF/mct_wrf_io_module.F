!BOP -------------------------------------------------------------------
!
! !MODULE: MCT_WRF_IO_MODULE
!
! !INTERFACE:

module MCT_WRF_IO_MODULE
!
! !DESCRIPTION:

! This documentation provides information specific to the MCT
! coupler reference implementation of the WRF I/O API.  Read
! the WRF Software Design and Implementation document section
! on the WRF I/O API prior to reading this information.
! .bp
! Re-gridding is implemented using the sparse matrix multiplication
! capability provided by MCT.  Currently, the matrix coefficients
! are read from a netCDF file.  This file contains global indices for the
! source and destination array associated weights (see
! the MCT reference guide for further information).
! The user provides the base name of this file to the OPEN routines
! described below.  Since the coefficients vary depending on 
! grid staggering, the full name of the file is constructed by
! combining the base name with the staggering characters passed
! to the read/write routines.
! .bp
! The re-gridding process is implemented by sending the appropriate data
! to the receiving component and then doing all the interpolation
! on the receiver.  In some cases, it would be more optimal to
! have the sender do calculations (to reduce communications bandwidth).
! This optimization is not currently provided.
! .bp
! The behavior of the implementation varies depending on whether
! it is called by model components that are executed sequentially
! or concurrently.  If they are executed concurrently, then every processor
! initializes this package with its component name since
! MCT_WORLD_INIT requires this information (actually it requires
! an integer ID associated with this name).  Note that in this case, every
! processor is associated with only one component.  From this information,
! the package
! determines names and component ids of all components that exchange boundary
! conditions with each other.  These are stored in an object of 
! class T_COMPONENT_TABLE defined in this module.  If the components execute
! sequentially then the component names are not unique to each processor.
! So instead the package is initialized with list of all components.  Here,
! a list of component ids is required by MCT_WORLD_INIT.
! .bp
! Each time an EXT_MCT_OPEN_FOR_WRITE_BEGIN 
! is called, the user is required to specify the component names of both
! the sender and receiver since MCT needs this information to
! setup up a router between the two components.  Technically, in the 
! concurrent case, only 
! one is required
! since the sender is "MY_COMP_NAME" which was obtained when the package
! was initialized.  However, for symmetry, the package will always require both.
! The same logic applies to EXT_MCT_OPEN_FOR_READ_BEGIN.
! .bp
! In the API routines, some of the arguments passed in are unused.
! See the notes in the subroutine documentation.

! !USES:

  use MCT_WRF_IO_UTILITIES_MODULE, only : &
    Assert

  use MCT_WRF_IO_UTILITIES_MODULE, only : MAX_STR

  use MCT_WRF_IO_PIPE_MODULE, only : &
    T_EXCHANGE_PIPE_CLASS => T_CLASS

  implicit none

! !PUBLIC MEMBER FUNCTIONS:
!
  public :: WRF_IO_API_OBJECT  ! An object containing all needed coupling info
  public :: CONSTRUCTOR        ! Constructor for the WRF_IO_API_OBJECT
  public :: EXECUTION_MODE     ! A class logical : Sequential or Concurrent
  public :: EXT_MCT_OPEN_BEGIN ! A routine to open a coupling data stream
  public :: FIRST_OPEN         ! True if first time read or write open occurred
  public :: RECV_DATA_REAL_2D  ! Routine to receive data from another component
  public :: SEND_DATA_REAL_2D  ! Routine to send data to another component

! !REVISION HISTORY:
! 05/28/03 - Dan Schaffer - Initial version

!EOP ___________________________________________________________________

#include "macros.h"

! Class variables
  character(LEN=MAX_STR), save          :: EXECUTION_MODE

! Object variables

  type T_COMPONENT_TABLE
    integer :: COMP_ID
    character(LEN=MAX_STR) :: COMP_NAME
  end type

  type T_CLASS
    integer :: COMM_GLOBAL
    integer :: COMM_LOCAL
    type (T_COMPONENT_TABLE), pointer :: COMPONENT_TABLE(:)
    integer :: NUM_COMPONENTS
    type (T_EXCHANGE_PIPE_CLASS), pointer :: EXCHANGE_PIPES(:)
    logical :: INIT_COMPLETE ! True if instance initialization is complete
    logical :: SEQUENTIAL ! True if the coupling is been executed by
                          ! sequentially executing model components

! For concurrent models where there each processor is associated with 
! only one component.
    character(LEN=MAX_STR) :: COMP_NAME
  end type

  type (T_CLASS), save :: WRF_IO_API_OBJECT

  interface CONSTRUCTOR
    module procedure CONSTRUCT_SEQUENTIAL
    module procedure CONSTRUCT_CONCURRENT
  end interface

contains

!BOP -------------------------------------------------------------------
!
! !IROUTINE: EXT_MCT_OPEN_BEGIN
!
! !INTERFACE:

subroutine EXT_MCT_OPEN_BEGIN ( &
  HIS_COMP_NAME, &
  COMM_GLOBAL, &
  COMM_LOCAL, &
  SYS_DEP_INFO, &
  DATA_HANDLE, &
  ISTATUS &
  )

!
! !USES:
  use MCT_WRF_IO_UTILITIES_MODULE, only : &
    WRF_IO_ERROR, &
    GET_KEY_VALUE

  implicit none
  
! !INPUT PARAMETERS:
!
  character(LEN=*), intent(in   ) :: HIS_COMP_NAME
  integer         , intent(in   ) :: COMM_GLOBAL
  integer         , intent(in   ) :: COMM_LOCAL
  character(LEN=*), intent(in   ) :: SYS_DEP_INFO

! !OUTPUT PARAMETERS:
  integer         , intent(  out) :: DATA_HANDLE

! !INPUT/OUTPUT PARAMETERS:
  integer         , intent(inout) :: ISTATUS

! !DESCRIPTION:  
! This routine opens
! a data stream for a model component that
! is sending/receiving re-gridded boundary conditions to/from 
! another model component.
! .bp
! The name of the component doing the open is parsed from SYS_DEP_INFO.
! The base name of the file containing the coefficients of the
! sparse matrix is also parsed from SYS_DEP_INFO.  
! If it has not yet been done, this routine calls OPEN_SETUP to do 
! some one-time setup.
! Finally, this routine uses ADD_PIPE to construct an entry 
! in the EXCHANGE_PIPES table.
!
! !DEFINED PARAMETERS:
! none

! !LOCAL VARIABLES:
  integer :: HIS_COMP_ID
  integer :: MY_COMP_ID
  character(len=MAX_STR) :: SPARSE_MATRIX_BASE_NAME
  character(len=MAX_STR) :: MY_COMP_NAME

! !SEE ALSO:
!  none

! !REVISION HISTORY:
! 05/28/03 - Dan Schaffer - Initial version
!EOP ___________________________________________________________________

  call GET_KEY_VALUE ( &
    SYS_DEP_INFO, &
    "SPARSE_MATRIX_BASE_NAME", &
    SPARSE_MATRIX_BASE_NAME, &
    ISTATUS, &
    IGNORE_CASE=.false. )

  if ( WRF_IO_ERROR  &
    ( ISTATUS /= 0, "Sparse Matrix Base filename key-value not found", &
    ISTATUS )) return

  call GET_KEY_VALUE ( &
    SYS_DEP_INFO, &
    "COMPONENT_NAME", &
    MY_COMP_NAME, &
    ISTATUS, &
    IGNORE_CASE=.false. )

  if ( WRF_IO_ERROR ( ISTATUS /= 0, "Need to specify name of originating component", ISTATUS )) return

  if (FIRST_OPEN(WRF_IO_API_OBJECT)) then

    call OPEN_SETUP ( &
      WRF_IO_API_OBJECT, &
      COMM_GLOBAL, &
      COMM_LOCAL )
  end if

  MY_COMP_ID = NAME_TO_ID ( &
    WRF_IO_API_OBJECT, &
    MY_COMP_NAME )

  if ( WRF_IO_ERROR ( MY_COMP_ID < 0, "Invalid Component", ISTATUS )) return

  HIS_COMP_ID = NAME_TO_ID ( &
    WRF_IO_API_OBJECT, &
    HIS_COMP_NAME )

  if ( WRF_IO_ERROR ( HIS_COMP_ID < 0, "Invalid Component", ISTATUS )) return

  call ADD_PIPE ( &
    WRF_IO_API_OBJECT, &
    SPARSE_MATRIX_BASE_NAME, &
    MY_COMP_ID, &
    HIS_COMP_ID, &
    DATA_HANDLE )

end subroutine EXT_MCT_OPEN_BEGIN

!BOP -------------------------------------------------------------------
!
! !IROUTINE: CONSTRUCT_CONCURRENT         
!
! !INTERFACE:

  subroutine CONSTRUCT_CONCURRENT ( &
    THIS, &
    MY_COMP_NAME)
!
! !USES:
!
  implicit none

! !INPUT PARAMETERS: 
!
  character(len=*), intent(in) :: MY_COMP_NAME

! !OUTPUT PARAMETERS: 
!
  type (T_CLASS), intent(out) :: THIS

! !DESCRIPTION:  Constructor for the MCT_WRF_IO class in the
! CONCURRENT case.  This
! is not a class in the true sense but is left in this form 
! in case we want to make it into a class.
!
! !DEFINED PARAMETERS:
! none

! !LOCAL VARIABLES: 
! none

! !SEE ALSO:
!  none

! !REVISION HISTORY:
! 05/31/03 - Dan Schaffer - Initial version
!EOP ___________________________________________________________________

  THIS%COMP_NAME = MY_COMP_NAME
  THIS%SEQUENTIAL = .false.
  THIS%INIT_COMPLETE = .false.

  nullify(THIS%COMPONENT_TABLE)
  nullify(THIS%EXCHANGE_PIPES)

end subroutine CONSTRUCT_CONCURRENT        

!BOP -------------------------------------------------------------------
!
! !IROUTINE: CONSTRUCT_SEQUENTIAL         
!
! !INTERFACE:

  subroutine CONSTRUCT_SEQUENTIAL ( &
    THIS, &
    COMP_NAMES )
!
! !USES:
!
  implicit none

! !INPUT PARAMETERS: 
!
  character(len=*), intent(in) :: COMP_NAMES(:)

! !OUTPUT PARAMETERS: 
!
  type (T_CLASS), intent(out) :: THIS

! !DESCRIPTION:  Constructor for the MCT_WRF_IO class for the sequential
! case.  This is not a class in the true sense but is left in this form 
! in case we want to make it into a class.
!
! !DEFINED PARAMETERS:
! none

! !LOCAL VARIABLES: 
  integer :: ISTATUS
  integer :: I

! !SEE ALSO:
!  none

! !REVISION HISTORY:
! 05/31/03 - Dan Schaffer - Initial version
!EOP ___________________________________________________________________

  THIS%NUM_COMPONENTS = size(COMP_NAMES)
  allocate(THIS%COMPONENT_TABLE(THIS%NUM_COMPONENTS), stat=ISTATUS)
  ASSERT(ISTATUS == 0)

  do I = 1, THIS%NUM_COMPONENTS
    THIS%COMPONENT_TABLE(I)%COMP_ID = I
    THIS%COMPONENT_TABLE(I)%COMP_NAME = COMP_NAMES(I)
  end do

  THIS%INIT_COMPLETE = .false.
  THIS%SEQUENTIAL = .true.

  nullify(THIS%EXCHANGE_PIPES)

end subroutine CONSTRUCT_SEQUENTIAL        

!BOP -------------------------------------------------------------------
!
! !IROUTINE: FIRST_OPEN         
!
! !INTERFACE:

logical function FIRST_OPEN ( &
  THIS)
!
! !USES:
!
  implicit none

! !INPUT PARAMETERS: 
!
  type (T_CLASS)   , intent(in)  :: THIS

! !DESCRIPTION:  Determine if this is the first time an open for read
! or write has occurred.
!
! !DEFINED PARAMETERS:
! none

! !LOCAL VARIABLES: 
! none

! !SEE ALSO:
!  none

! !REVISION HISTORY:
! 05/31/03 - Dan Schaffer - Initial version
!EOP ___________________________________________________________________

  FIRST_OPEN = .not. THIS%INIT_COMPLETE

end function FIRST_OPEN        

!BOP -------------------------------------------------------------------
!
! !IROUTINE: OPEN_SETUP         
!
! !INTERFACE:

subroutine OPEN_SETUP ( &
  THIS, &
  COMM_GLOBAL, &
  COMM_LOCAL )
!
! !USES:
!
  use MCT_WRF_IO_UTILITIES_MODULE, only : VALUE_IN_ARRAY
  use MCT_WRF_IO_UTILITIES_MODULE, only : NUM_UNIQUE

  use m_MCTWorld, only : MCT_WORLD_INIT => init

  use MCT_WRF_IO_INTERPOLATOR_MODULE, only : &
    INIT_INTERPOLATOR_CLASS => INITIALIZE

  implicit none

! !INPUT PARAMETERS: 
!
  integer, intent(in) :: COMM_GLOBAL
  integer, intent(in) :: COMM_LOCAL

! !OUTPUT PARAMETERS: 
!
  type (T_CLASS), intent(inout) :: THIS

! !DESCRIPTION:  Do all the package initialization that has to
! wait until the first open occurs.  This includes 
! initializing MCT and creating the list of components in
! the concurrent case.
!
! !DEFINED PARAMETERS:
! none

! !LOCAL VARIABLES: 
  integer :: NPROCS
  integer :: NPROC
  integer :: COMP_INDEX
  integer :: ERROR
  integer :: ISTATUS
  integer :: I
  character(len=MAX_STR), allocatable :: ALL_COMP_NAMES(:)
  include 'mpif.h'

! !SEE ALSO:
!  none

! !REVISION HISTORY:
! 05/31/03 - Dan Schaffer - Initial version
!EOP ___________________________________________________________________

  THIS%COMM_GLOBAL = COMM_GLOBAL
  THIS%COMM_LOCAL = COMM_LOCAL
  THIS%INIT_COMPLETE = .true.

  if (THIS%SEQUENTIAL) then
    call MCT_WORLD_INIT ( &
      THIS%NUM_COMPONENTS, &
      THIS%COMM_GLOBAL, &
      THIS%COMM_LOCAL, &
      MYIDS = THIS%COMPONENT_TABLE(:)%COMP_ID )
  else

  ! Concurrent.  Create the list of components.
    call MPI_COMM_SIZE ( &
      THIS%COMM_GLOBAL, &
      NPROCS, &
      ERROR )
  
    ! Okay here
    allocate(ALL_COMP_NAMES(NPROCS), stat=ISTATUS)
    ASSERT(ISTATUS == 0)
  
    call MPI_ALLGATHER ( &
      THIS%COMP_NAME, &
      MAX_STR, &
      MPI_CHARACTER, &
      ALL_COMP_NAMES, &
      MAX_STR, &
      MPI_CHARACTER, &
      THIS%COMM_GLOBAL, &
      ERROR )
  
    THIS%NUM_COMPONENTS = NUM_UNIQUE ( ALL_COMP_NAMES )
  
    allocate(THIS%COMPONENT_TABLE(THIS%NUM_COMPONENTS), stat=ISTATUS)
    ASSERT(ISTATUS == 0)
  
    COMP_INDEX = 0
  
    do NPROC = 1, NPROCS
      if (.not. VALUE_IN_ARRAY ( &
        ALL_COMP_NAMES(NPROC), &
        THIS%COMPONENT_TABLE(1:COMP_INDEX)%COMP_NAME)) then
  
        COMP_INDEX = COMP_INDEX + 1
        THIS%COMPONENT_TABLE(COMP_INDEX)%COMP_ID = COMP_INDEX
        THIS%COMPONENT_TABLE(COMP_INDEX)%COMP_NAME = ALL_COMP_NAMES(NPROC)
      end if
    end do
  
    deallocate(ALL_COMP_NAMES, stat=ISTATUS)
    ASSERT(ISTATUS == 0)

    do COMP_INDEX = 1, THIS%NUM_COMPONENTS

      if (trim(THIS%COMPONENT_TABLE(COMP_INDEX)%COMP_NAME) == &
        trim(THIS%COMP_NAME)) then

        call MCT_WORLD_INIT ( &
          THIS%NUM_COMPONENTS, &
          THIS%COMM_GLOBAL, &
          THIS%COMM_LOCAL, &
          COMP_INDEX )
      end if
    end do
  
  end if

  call INIT_INTERPOLATOR_CLASS(THIS%NUM_COMPONENTS, THIS%SEQUENTIAL)

end subroutine OPEN_SETUP        

!BOP -------------------------------------------------------------------
!
! !IROUTINE: ADD_PIPE         
!
! !INTERFACE:
subroutine ADD_PIPE ( &
  THIS, &
  SPARSE_MATRIX_BASE_NAME, &
  MY_COMP_ID, &
  HIS_COMP_ID, &
  DATA_HANDLE )

!
! !USES:
!
  use MCT_WRF_IO_PIPE_MODULE, only : &
    CONSTRUCT_EXCHANGE_PIPE => CONSTRUCTOR

  implicit none

! !INPUT PARAMETERS: 
!
  character(len=*), intent(in) :: SPARSE_MATRIX_BASE_NAME
  integer, intent(in) :: MY_COMP_ID
  integer, intent(in) :: HIS_COMP_ID

! !OUTPUT PARAMETERS: 
!
  integer, intent(out) :: DATA_HANDLE

! !INPUT/OUTPUT PARAMETERS: 
!
  type (T_CLASS), intent(out) :: THIS

! !DESCRIPTION:  Construct and add an exchange pipe to THIS
!
! !DEFINED PARAMETERS:
! none

! !LOCAL VARIABLES: 
  type (T_EXCHANGE_PIPE_CLASS) :: NEW_PIPE
  integer :: ISTATUS

! !SEE ALSO:
!  none

! !REVISION HISTORY:
! 05/31/03 - Dan Schaffer - Initial version
!EOP ___________________________________________________________________

  call CONSTRUCT_EXCHANGE_PIPE ( &
    NEW_PIPE, &
    THIS%COMM_LOCAL, &
    SPARSE_MATRIX_BASE_NAME, &
    MY_COMP_ID, &
    HIS_COMP_ID, &
    ISTATUS )

  ASSERT(ISTATUS == 0)

  call EXPAND_EXCHANGE_PIPES ( &
    THIS%EXCHANGE_PIPES, &
    ITEMS= (/NEW_PIPE/) )

  DATA_HANDLE = size(THIS%EXCHANGE_PIPES)

end subroutine ADD_PIPE        

!BOP -------------------------------------------------------------------
!
! !IROUTINE: EXPAND_EXCHANGE_PIPES         
!
! !INTERFACE:

subroutine EXPAND_EXCHANGE_PIPES ( IN, ITEMS, AFTER )

#define IN_TYPE type(T_EXCHANGE_PIPE_CLASS)

#include "arrays.xpnd.code"

#undef IN_TYPE

end subroutine EXPAND_EXCHANGE_PIPES

!BOP -------------------------------------------------------------------
!
! !IROUTINE: NAME_TO_ID         
!
! !INTERFACE:

integer function NAME_TO_ID ( &
  THIS, &
  COMP_NAME )
!
! !USES:
!
  implicit none

! !INPUT PARAMETERS: 
!
  type (T_CLASS)   , intent(in)  :: THIS

! !DESCRIPTION:  Convert a component name to a component id using
! the COMPONENT_TABLE in THIS.
!
! !DEFINED PARAMETERS:
! none

! !LOCAL VARIABLES: 
  character(len=*), intent(in) :: COMP_NAME

  integer :: I

! !SEE ALSO:
!  none

! !REVISION HISTORY:
! 05/31/03 - Dan Schaffer - Initial version
!EOP ___________________________________________________________________

  do I = 1, size(THIS%COMPONENT_TABLE)
    if (THIS%COMPONENT_TABLE(I)%COMP_NAME == COMP_NAME) then
      NAME_TO_ID = I
      return
    end if
  end do

  NAME_TO_ID = -1

end function NAME_TO_ID        

!BOP -------------------------------------------------------------------
!
! !IROUTINE: SEND_DATA_REAL_2D         
!
! !INTERFACE:

subroutine SEND_DATA_REAL_2D ( &
  THIS, &
  DATA_HANDLE, &
  STAGGER, &
  PATCH_START, &
  PATCH_END, &
  MEMORY_START, &
  MEMORY_END, &
  VAR_NAME, &
  FIELD, &
  ISTATUS)

#define IN_TYPE real
#include "send_data.code"
#undef IN_TYPE

end subroutine SEND_DATA_REAL_2D        

!BOP -------------------------------------------------------------------
!
! !IROUTINE: RECV_DATA_REAL_2D         
!
! !INTERFACE:

subroutine RECV_DATA_REAL_2D ( &
  THIS, &
  DATA_HANDLE, &
  STAGGER, &
  PATCH_START, &
  PATCH_END, &
  MEMORY_START, &
  MEMORY_END, &
  VAR_NAME, &
  FIELD, &
  ISTATUS)

#define IN_TYPE real
#include "recv_data.code"
#undef IN_TYPE

end subroutine RECV_DATA_REAL_2D        

end module MCT_WRF_IO_MODULE
