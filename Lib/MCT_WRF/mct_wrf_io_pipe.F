!BOP -------------------------------------------------------------------
!
! !MODULE: MCT_WRF_IO_PIPE_MODULE - Class for a data exchange pipe
!          
! !INTERFACE:

 module MCT_WRF_IO_PIPE_MODULE
!
! !DESCRIPTION:
!   Class that provides the capability to open a data exchange pipe
! between two model components.

! !USES:
  use MCT_WRF_IO_UTILITIES_MODULE, only : &
    Assert

  use MCT_WRF_IO_UTILITIES_MODULE, only : MAX_STR

  use MCT_WRF_IO_INTERPOLATOR_MODULE, only : &
    T_INTERPOLATOR_CLASS => T_CLASS

  implicit none

#include "macros.h"

  private	! except

! !PUBLIC TYPES:

  public :: T_CLASS 

!

! !PUBLIC MEMBER FUNCTIONS:
!
  public :: CONSTRUCTOR        ! Construct the object
  public :: DESTRUCTOR         ! Destruct the object
  public :: COMMIT             ! Commit the pipe for reading or writing
  public :: COMMITTED          ! True if the pipe has been committed
  public :: EXPAND             ! Add a variable to the pipe
  public :: SEND               ! Send data into a pipe
  public :: RECV               ! Receive data from a pipe

! !DEFINED PARAMETERS:
! none

! !SEE ALSO:
!  none

! !REVISION HISTORY:
! 05/21/03 - Dan Schaffer - Initial version

!EOP ___________________________________________________________________

  type T_CLASS
    integer :: MY_ID           ! Id of my component
    integer :: HIS_ID          ! Id of component with which exchange is done

    character(LEN=MAX_STR) :: &
       SPARSE_MATRIX_BASE_NAME ! Base name of file containing
                               ! weights and indices for interpolation
    integer :: COMM_LOCAL      ! Communicator my process belongs to

    type (T_INTERPOLATOR_CLASS), pointer :: &
      INTERPOLATORS(:)         ! Interpolators associated with THIS

    logical :: COMMITTED       ! False initially.  
    			       ! True after a read/write commit
  end type

  interface SEND
    module procedure SEND_REAL_2D
  end interface

  interface RECV
    module procedure RECV_REAL_2D
  end interface

  contains

!BOP -------------------------------------------------------------------
!
! !IROUTINE: CONSTRUCTOR         - Constructor
!
! !INTERFACE:

 subroutine CONSTRUCTOR ( &
   THIS, &
   COMM_LOCAL, &
   SPARSE_MATRIX_BASE_NAME, &
   MY_ID, &
   HIS_ID, &
   ISTATUS)
!
! !USES:
!
  implicit none

! !INPUT PARAMETERS: 
!
  integer          , intent(in)  :: COMM_LOCAL
  character(LEN=*) , intent(in)  :: SPARSE_MATRIX_BASE_NAME
  integer          , intent(in)  :: MY_ID
  integer          , intent(in)  :: HIS_ID

! !OUTPUT PARAMETERS: 
!
  type (T_CLASS)   , intent(out) :: THIS

! !INPUT/OUTPUT PARAMETERS: 
!
  integer          , intent(inout)  :: ISTATUS

! !DESCRIPTION:  Construct a pipe for exchanging data between two components.
! Store away the information passed in.  Nullify the list of 
! interpolators associated with this pipe.
! 
!
! !DEFINED PARAMETERS:
! none

! !LOCAL VARIABLES: 
! none

! !SEE ALSO:
! none

! !REVISION HISTORY:
! 05/21/03 - Dan Schaffer - Initial version
!EOP ___________________________________________________________________

  nullify(THIS%INTERPOLATORS)
  THIS%COMMITTED = .false.
  THIS%SPARSE_MATRIX_BASE_NAME = SPARSE_MATRIX_BASE_NAME
  THIS%MY_ID = MY_ID
  THIS%HIS_ID = HIS_ID
  THIS%COMM_LOCAL = COMM_LOCAL
 end subroutine CONSTRUCTOR        

!BOP -------------------------------------------------------------------
!
! !IROUTINE: DESTRUCTOR         - Destructor
!
! !INTERFACE:

 subroutine DESTRUCTOR ( &
   THIS)
!
! !USES:
!
  implicit none

! !INPUT/OUTPUT PARAMETERS: 
!
  type (T_CLASS)   , intent(inout) :: THIS

! !DESCRIPTION:  Destruct a pipe.
! 
!
! !DEFINED PARAMETERS:
! none

! !LOCAL VARIABLES: 
! none

! !SEE ALSO:
! none

! !REVISION HISTORY:
! 05/21/03 - Dan Schaffer - Initial version
!EOP ___________________________________________________________________

  print *, 'MCT_WRF_IO_PIPE_MODULE DESTRUCTOR routine stubbed for now'
 end subroutine DESTRUCTOR        

!BOP -------------------------------------------------------------------
!
! !IROUTINE: COMMIT         
!
! !INTERFACE:

 subroutine COMMIT (THIS)
!
! !USES:
!
  implicit none

! !INPUT PARAMETERS: 
!
! none

! !OUTPUT PARAMETERS: 
!
  type (T_CLASS)   , intent(out) :: THIS

! !DESCRIPTION:  Mark a pipe as committed for reading or writing.
! This will indicate to a caller that if a read or write occurs
! it is time to send/receive data instead of continue adding the variable
! to the pipe.
!
! !DEFINED PARAMETERS:
! none

! !LOCAL VARIABLES: 
! none

! !SEE ALSO:
!  none

! !REVISION HISTORY:
! 05/23/03 - Dan Schaffer - Initial version
!EOP ___________________________________________________________________

  THIS%COMMITTED = .true.
 end subroutine COMMIT        

!BOP -------------------------------------------------------------------
!
! !IROUTINE: COMMIT         
!
! !INTERFACE:

 function COMMITTED (THIS)
!
! !USES:
!
  implicit none

! !INPUT PARAMETERS: 
!
  type (T_CLASS)   , intent(in) :: THIS

! !OUTPUT PARAMETERS: 
!
  logical          :: COMMITTED

! !DESCRIPTION:  True if a pipe has been committed for
! reading or writing.
!
! !DEFINED PARAMETERS:
! none

! !LOCAL VARIABLES: 
! none

! !SEE ALSO:
!  none

! !REVISION HISTORY:
! 05/23/03 - Dan Schaffer - Initial version
!EOP ___________________________________________________________________

  COMMITTED = THIS%COMMITTED
 end function COMMITTED        

!BOP -------------------------------------------------------------------
!
! !IROUTINE: EXPAND
!
! !INTERFACE:

 subroutine EXPAND ( &
   THIS, &
   VAR_NAME, &
   VAR_TYPE, &
   STAGGER, &
   DOMAIN_START, &
   DOMAIN_END, &
   PATCH_START, &
   PATCH_END, &
   RECEIVER, &
   ISTATUS &
   )
!
! !USES:

  use MCT_WRF_IO_INTERPOLATOR_MODULE, only : &
    CONSTRUCT_INTERPOLATOR => CONSTRUCTOR, &
    EXTEND_SRC_INTERPOLATOR => EXTEND_SRC, &
    EXTEND_DST_INTERPOLATOR => EXTEND_DST
!
  implicit none

! !INPUT PARAMETERS: 
!
  character(LEN=*) , intent(in)  :: VAR_NAME
  integer          , intent(in)  :: VAR_TYPE
  character(LEN=*) , intent(in)  :: STAGGER
  integer          , intent(in)  :: DOMAIN_START(:)
  integer          , intent(in)  :: DOMAIN_END(:)
  integer          , intent(in)  :: PATCH_START(:)
  integer          , intent(in)  :: PATCH_END(:)
  logical          , intent(in)  :: RECEIVER

! !INPUT/OUTPUT PARAMETERS: 
!
  type (T_CLASS)   , intent(inout) :: THIS
  integer          , intent(inout) :: ISTATUS

! !DESCRIPTION:  Search for the interpolator with given ID.  If not 
! found, construct the interpolator.  If found, extend it with the 
! given VAR_NAME.  For now we are assuming that patches are 2-D.
!
! !DEFINED PARAMETERS:
! none

! !LOCAL VARIABLES: 
  integer :: INTERP_NUM
  integer :: STARTS (PATCH_END(2) - PATCH_START(2) + 1)
  integer :: LENGTHS(PATCH_END(2) - PATCH_START(2) + 1)

  type (T_INTERPOLATOR_CLASS) :: NEW_INTERPOLATOR

  character(LEN=MAX_STR) :: ID

! !SEE ALSO:
!  none

! !REVISION HISTORY:
! 05/23/03 - Dan Schaffer - Initial version
!EOP ___________________________________________________________________

  call INTERPOLATOR_ID(THIS, STAGGER, ID)

  if (FIND_INTERPOLATOR(THIS%INTERPOLATORS, ID, INTERP_NUM)) then
    if (RECEIVER) then
      call EXTEND_DST_INTERPOLATOR ( &
        THIS%INTERPOLATORS(INTERP_NUM), &
        VAR_NAME, &
        VAR_TYPE, &
        ISTATUS )
    else
      call EXTEND_SRC_INTERPOLATOR ( &
        THIS%INTERPOLATORS(INTERP_NUM), &
        VAR_NAME, &
        VAR_TYPE, &
        ISTATUS )
    end if
  else
    call COMPUTE_VECTORS ( &
      DOMAIN_START, &
      DOMAIN_END, &
      PATCH_START, &
      PATCH_END, &
      STARTS, &
      LENGTHS, &
      ISTATUS )

    call CONSTRUCT_INTERPOLATOR ( &
      NEW_INTERPOLATOR, &
      VAR_NAME, &
      VAR_TYPE, &
      STARTS, &
      LENGTHS, &
      THIS%COMM_LOCAL, &
      THIS%MY_ID, &
      THIS%HIS_ID, &
      ID, &
      RECEIVER, &
      ISTATUS )

    call EXPAND_INTERPOLATORS ( &
      THIS%INTERPOLATORS, &
      ITEMS = (/ NEW_INTERPOLATOR /) )

  end if
    
 end subroutine EXPAND

!BOP -------------------------------------------------------------------
!
! !IROUTINE: SEND_REAL_2D         
!
! !INTERFACE:

 subroutine SEND_REAL_2D ( &
   THIS, &
   STAGGER, &
   VAR_NAME, &
   DATA, &
   ISTATUS )

#define IN_TYPE real
#include "pipe_send.code"
#undef IN_TYPE

 end subroutine SEND_REAL_2D        

!BOP -------------------------------------------------------------------
!
! !IROUTINE: RECV_REAL_2D         
!
! !INTERFACE:

 subroutine RECV_REAL_2D ( &
   THIS, &
   STAGGER, &
   VAR_NAME, &
   DATA, &
   ISTATUS )

#define IN_TYPE real
#include "pipe_recv.code"
#undef IN_TYPE

 end subroutine RECV_REAL_2D        

!BOP -------------------------------------------------------------------
!
! !IROUTINE: FIND_INTERPOLATOR         
!
! !INTERFACE:

 logical function FIND_INTERPOLATOR (INTERPOLATOR_LIST, ID, INTERP_NUM)
!
! !USES:
  use MCT_WRF_IO_INTERPOLATOR_MODULE, only : &
    INTERPOLATOR_MATCHES => MATCHES
!
  implicit none

! !INPUT PARAMETERS: 
!
  type (T_INTERPOLATOR_CLASS), pointer     :: INTERPOLATOR_LIST(:)
  character(LEN=*), intent(in)  :: ID

! !OUTPUT PARAMETERS: 
!
  integer          , intent(out) :: INTERP_NUM

! !DESCRIPTION:  Look for an interpolator with the given ID in
! INTERPOLATOR_LIST.  If found, fill INTERP_NUM.  Otherwise
! just return false.
!
! !DEFINED PARAMETERS:
! none

! !LOCAL VARIABLES: 
  integer :: I

! !SEE ALSO:
!  none

! !REVISION HISTORY:
! 05/23/03 - Dan Schaffer - Initial version
!EOP ___________________________________________________________________

  if (.not. associated(INTERPOLATOR_LIST)) then
    FIND_INTERPOLATOR = .false.
    return
  else
    do I = 1, size(INTERPOLATOR_LIST)
      if (INTERPOLATOR_MATCHES(INTERPOLATOR_LIST(I), ID)) then
        FIND_INTERPOLATOR = .true.
        INTERP_NUM = I
        return
      end if
    end do
    FIND_INTERPOLATOR = .false.
    INTERP_NUM = -1
  end if
 end function FIND_INTERPOLATOR        

!BOP -------------------------------------------------------------------
!
! !IROUTINE: COMPUTE_VECTORS         
!
! !INTERFACE:

 subroutine COMPUTE_VECTORS ( &
   DOMAIN_START,  &
   DOMAIN_END,  &
   PATCH_START,  &
   PATCH_END,  &
   STARTS,  &
   LENGTHS,  &
   ISTATUS )
!
! !USES:
!
  implicit none

! !INPUT PARAMETERS: 
!
  integer          , intent(in)  :: DOMAIN_START(:)
  integer          , intent(in)  :: DOMAIN_END(:)
  integer          , intent(in)  :: PATCH_START(:)
  integer          , intent(in)  :: PATCH_END(:)

! !OUTPUT PARAMETERS: 
!
  integer          , intent(out)  :: STARTS(:)
  integer          , intent(out)  :: LENGTHS(:)
  integer          , intent(out)  :: ISTATUS

! !DESCRIPTION:  Figure out the starts and lengths of the vectors
! that need to be described to MCT.
!
! !DEFINED PARAMETERS:
! none

! !LOCAL VARIABLES: 
  integer :: I
  integer :: IM_WORLD

! !SEE ALSO:
!  none

! !REVISION HISTORY:
! 05/23/03 - Dan Schaffer - Initial version
!EOP ___________________________________________________________________

   LENGTHS = PATCH_END(1) - PATCH_START(1) + 1
   IM_WORLD = DOMAIN_END(1) - DOMAIN_START(1) + 1
   STARTS(1) = (PATCH_START(2) - 1) * IM_WORLD + PATCH_START(1)

   do I = 2, size(STARTS)
     STARTS(I) = STARTS(I-1) + IM_WORLD
   end do
 end subroutine COMPUTE_VECTORS        

!BOP -------------------------------------------------------------------
!
! !IROUTINE: EXPAND_INTERPOLATORS         
!
! !INTERFACE:

subroutine EXPAND_INTERPOLATORS ( &
  IN, &
  ITEMS, &
  AFTER)

#define IN_TYPE type(T_INTERPOLATOR_CLASS)
#include "arrays.xpnd.code"
#undef IN_TYPE

end subroutine EXPAND_INTERPOLATORS        

!BOP -------------------------------------------------------------------
!
! !IROUTINE: COMPUTE_INTERPOLATOR_NAME         
!
! !INTERFACE:

subroutine INTERPOLATOR_ID ( &
  THIS, &
  STAGGER, &
  ID)
!
! !USES:
!
  implicit none

! !INPUT PARAMETERS: 
!
  type (T_CLASS)   , intent(in)  :: THIS
  character(len=*) , intent(in)  :: STAGGER

! !OUTPUT PARAMETERS: 
!
  character(len=*), intent(out) :: ID

! !DESCRIPTION:  Set INTERPOLATOR_ID to the concatentation of
! THIS%SPARSE_MATRIX_BASE_NAME and STAGGER.
!
! !DEFINED PARAMETERS:
! none

! !LOCAL VARIABLES: 
! none

! !SEE ALSO:
!  none

! !REVISION HISTORY:
! 05/01/03 - Dan Schaffer - Initial version
!EOP ___________________________________________________________________

  if (trim(STAGGER) == "") then
    ID = trim(THIS%SPARSE_MATRIX_BASE_NAME)
  else
    ID = trim(THIS%SPARSE_MATRIX_BASE_NAME) // "." // trim(STAGGER)
  end if

end subroutine INTERPOLATOR_ID        

 end module MCT_WRF_IO_PIPE_MODULE
