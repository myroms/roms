#include "macros.h"
!BOP -------------------------------------------------------------------
!
! !ROUTINE: EXT_MCT_IOINIT
!
! !INTERFACE:

subroutine EXT_MCT_IOINIT ( &
  SYS_DEP_INFO, &
  ISTATUS &
  )
!
! !USES:

  use MCT_WRF_IO_MODULE, only : WRF_IO_API_OBJECT
  use MCT_WRF_IO_MODULE, only : CONSTRUCTOR
  use MCT_WRF_IO_MODULE, only : EXECUTION_MODE

  use MCT_WRF_IO_UTILITIES_MODULE, only : WRF_IO_ERROR
  use MCT_WRF_IO_UTILITIES_MODULE, only : GET_KEY_VALUE
  use MCT_WRF_IO_UTILITIES_MODULE, only : MAX_STR
  use MCT_WRF_IO_UTILITIES_MODULE, only : TO_LOWER
  use MCT_WRF_IO_UTILITIES_MODULE, only : EXPAND
  use MCT_WRF_IO_UTILITIES_MODULE, only : Assert

  implicit none

! !INPUT PARAMETERS:
!
  character(LEN=*), intent(in   ) :: SYS_DEP_INFO

! !OUTPUT PARAMETERS:
!
  integer         , intent(  out) :: ISTATUS

! !DESCRIPTION:  Initialize the MCT coupler implementation of the WRF
! I/O API.  Parse SYS_DEP_INFO to determine if the coupler is being
! used by a sequential or concurrent calling sequence. 
! Parse the list of all component names from SYS_DEP_INFO.
! A constructor for the WRF_IO_API class is then called.

!
! !DEFINED PARAMETERS:
! none

! !LOCAL VARIABLES:
  character(LEN=MAX_STR) :: COMP_NAME

  integer, parameter :: MAX_COMPONENTS = 9
  character(LEN=MAX_STR) :: COMP_NAME_KEY
  character(LEN=MAX_STR), pointer :: COMP_NAMES(:)
  logical, save :: FIRST_TIME = .true.


! !SEE ALSO:

! !REVISION HISTORY:
! 05/28/03 - Dan Schaffer - Initial version
!EOP ___________________________________________________________________

  if (FIRST_TIME) then
    FIRST_TIME = .false.
  else
    return
  end if

  call GET_KEY_VALUE ( &
    SYS_DEP_INFO, &
    "COMPONENT1", &
    COMP_NAME, &
    ISTATUS, &
    IGNORE_CASE=.true.)

  ! Sequential case
  if (ISTATUS == 0) then

    EXECUTION_MODE = 'SEQUENTIAL'

    allocate(COMP_NAMES(0), stat=ISTATUS)

    call EXPAND ( &
       COMP_NAMES, &
       ITEMS = (/ COMP_NAME /))

    do while (ISTATUS == 0)
  
      write(COMP_NAME_KEY, '("COMPONENT", I1.1)') size(COMP_NAMES) + 1
      call GET_KEY_VALUE ( &
        SYS_DEP_INFO, &
        trim(COMP_NAME_KEY), &
        COMP_NAME, &
        ISTATUS, &
        IGNORE_CASE=.true.)

      if (ISTATUS ==0) then

        call EXPAND ( &
          COMP_NAMES, &
	  ITEMS = (/ COMP_NAME /))
  
      end if
  
    end do

    call CONSTRUCTOR ( &
      WRF_IO_API_OBJECT, &
      COMP_NAMES )

    deallocate(COMP_NAMES, stat=ISTATUS)
    ASSERT(ISTATUS == 0)

  else

    EXECUTION_MODE = 'CONCURRENT'

    call GET_KEY_VALUE ( &
      SYS_DEP_INFO, &
      "COMPONENT_NAME", &
      COMP_NAME, &
      ISTATUS, &
      IGNORE_CASE=.true.)

    call CONSTRUCTOR ( &
      WRF_IO_API_OBJECT, &
      COMP_NAME)

  end if

end subroutine EXT_MCT_IOINIT

!BOP -------------------------------------------------------------------
!
! !ROUTINE: EXT_MCT_OPEN_FOR_READ_BEGIN
!
! !INTERFACE:

subroutine EXT_MCT_OPEN_FOR_READ_BEGIN ( &
  HIS_COMP_NAME, &
  COMM_GLOBAL, &
  COMM_LOCAL, &
  SYS_DEP_INFO, &
  DATA_HANDLE, &
  ISTATUS &
  )

!
! !USES:
  use MCT_WRF_IO_MODULE, only : EXT_MCT_OPEN_BEGIN

  use MCT_WRF_IO_UTILITIES_MODULE, only : &
    WRF_IO_ERROR

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
! is exchanging re-gridded boundary conditions with another model component.
! See EXT_MCT_OPEN_BEGIN for the details.
!
! !DEFINED PARAMETERS:
! none

! !LOCAL VARIABLES:
! none

! !SEE ALSO:
!  none

! !REVISION HISTORY:
! 05/28/03 - Dan Schaffer - Initial version
!EOP ___________________________________________________________________

call EXT_MCT_OPEN_BEGIN ( &
  HIS_COMP_NAME, &
  COMM_GLOBAL, &
  COMM_LOCAL, &
  SYS_DEP_INFO, &
  DATA_HANDLE, &
  ISTATUS &
  )

end subroutine EXT_MCT_OPEN_FOR_READ_BEGIN

!BOP -------------------------------------------------------------------
!
! !ROUTINE: EXT_MCT_OPEN_FOR_READ_COMMIT
!
! !INTERFACE:

subroutine EXT_MCT_OPEN_FOR_READ_COMMIT ( &
  DATA_HANDLE, &
  ISTATUS )

!
! !USES:

  use MCT_WRF_IO_MODULE, only : WRF_IO_API_OBJECT

  use MCT_WRF_IO_PIPE_MODULE, only : &
    COMMIT_PIPE => COMMIT

  implicit none

! !INPUT PARAMETERS:
!
  integer,          intent(in   ) :: DATA_HANDLE
  integer,          intent(inout) :: ISTATUS

! !DESCRIPTION:  
! A call to this routine indicates that the "training" phase
! is complete.  Now the package knows all the variables
! for which communication will be aggregated for each exchange pipe.

!
! !DEFINED PARAMETERS:
! none

! !LOCAL VARIABLES:
! none

! !SEE ALSO:
!  none

! !REVISION HISTORY:
! 05/28/03 - Dan Schaffer - Initial version
!EOP ___________________________________________________________________

  ISTATUS = 0
  call COMMIT_PIPE(WRF_IO_API_OBJECT%EXCHANGE_PIPES(DATA_HANDLE))

end subroutine EXT_MCT_OPEN_FOR_READ_COMMIT

!BOP -------------------------------------------------------------------
!
! !ROUTINE: EXT_MCT_OPEN_FOR_READ         
!
! !INTERFACE:

 subroutine EXT_MCT_OPEN_FOR_READ ( & 
  HIS_COMP_NAME, &
  COMM_GLOBAL, &
  COMM_LOCAL, &
  SYS_DEP_INFO, &
  DATA_HANDLE, &
  ISTATUS &
  )

!
! !USES:
!
  use MCT_WRF_IO_UTILITIES_MODULE, only : &
    WRF_IO_ERROR

  implicit none

! !INPUT PARAMETERS: 
!
  character(LEN=*), intent(in   ) :: HIS_COMP_NAME
  integer         , intent(in   ) :: COMM_GLOBAL
  integer         , intent(in   ) :: COMM_LOCAL
  character(LEN=*), intent(in   ) :: SYS_DEP_INFO

! !OUTPUT PARAMETERS: 
!
  integer         , intent(  out) :: DATA_HANDLE
  integer,          intent(inout) :: ISTATUS

! !DESCRIPTION:  
! Although not currently implemented, this routine would eventually
! provide a means to open a data stream for immediate reception
! of re-gridded data.  Under
! this scenario, there would be no "training" as discussed in 
! the API documentation.  For this implementation, that means
! there would be no aggregation of inter-component communication; 
! each boundary array would be received by itself.
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

  if ( WRF_IO_ERROR  &
    ( .true., "EXT_MCT_OPEN_FOR_READ routine stubbed for now", &
    ISTATUS )) return

 end subroutine EXT_MCT_OPEN_FOR_READ        

!BOP -------------------------------------------------------------------
!
! !ROUTINE: EXT_MCT_READ_FIELD
!
! !INTERFACE:

subroutine EXT_MCT_READ_FIELD ( &
  DATA_HANDLE, &
  DATE_STR, &
  VAR_NAME, &
  FIELD, &
  FIELD_TYPE, &
  COMM_GLOBAL, &
  COMM_LOCAL, &
  DOMAIN_DESC, &
  MEMORY_ORDER, &
  STAGGER, &
  DIM_NAMES, &
  DOMAIN_START, &
  DOMAIN_END, &
  MEMORY_START, &
  MEMORY_END, &
  PATCH_START, &
  PATCH_END, &
  ISTATUS )

!
! !USES:
  use MCT_WRF_IO_MODULE, only : WRF_IO_API_OBJECT
  use MCT_WRF_IO_MODULE, only : FIRST_OPEN
  use MCT_WRF_IO_MODULE, only : RECV_DATA_REAL_2D

  use MCT_WRF_IO_UTILITIES_MODULE, only : Assert
  use MCT_WRF_IO_UTILITIES_MODULE, only : &
    WRF_IO_ERROR

  use MCT_WRF_IO_PIPE_MODULE, only : &
    PIPE_COMMITTED => COMMITTED, &
    EXPAND_PIPE => EXPAND

  implicit none

! !INPUT PARAMETERS:
!
  integer,          intent(in   ) :: DATA_HANDLE
  character(LEN=*), intent(in   ) :: DATE_STR        ! Not used
  character(LEN=*), intent(in   ) :: VAR_NAME

  integer,          intent(in   ) :: FIELD_TYPE      
  integer,          intent(in   ) :: COMM_GLOBAL     
  integer,          intent(in   ) :: COMM_LOCAL     
  integer,          intent(in   ) :: DOMAIN_DESC     ! Not used
  character(LEN=*), intent(in   ) :: MEMORY_ORDER    ! Not used
  character(LEN=*), intent(in   ) :: STAGGER
  character(LEN=*), intent(in   ) :: DIM_NAMES       ! Not used
  integer,          intent(in   ) :: DOMAIN_START(len_trim(MEMORY_ORDER))
  integer,          intent(in   ) :: DOMAIN_END  (len_trim(MEMORY_ORDER))
  integer,          intent(in   ) :: MEMORY_START(len_trim(MEMORY_ORDER))
  integer,          intent(in   ) :: MEMORY_END  (len_trim(MEMORY_ORDER))
  integer,          intent(in   ) :: PATCH_START (len_trim(MEMORY_ORDER))
  integer,          intent(in   ) :: PATCH_END   (len_trim(MEMORY_ORDER))

! !OUTPUT PARAMETERS:
  real,               intent(  out) :: FIELD( &
    MEMORY_START(1):MEMORY_END(1),MEMORY_START(2):MEMORY_END(2))

! !INPUT/OUTPUT PARAMETERS:
  integer,          intent(inout) :: ISTATUS
  
! !DESCRIPTION:  
! This routine does one of several things depending on the conditions
! upon which it is called.  If it is called prior to
! a call to EXT_MCT_OPEN_FOR_READ_BEGIN then there is an error condition.
! If it is called after a call to EXT_MCT_OPEN_FOR_READ_BEGIN but
! before a call to EXT_MCT_OPEN_FOR_READ_COMMIT then
! the pipe between this component and the sending component is expanded
! to include the variable whose name is VAR_NAME.  If 
! EXT_MCT_OPEN_FOR_READ_COMMIT has been called then RECV_DATA_REAL_2D
! handles the communication, re-gridding and unpacking of the data.
! 

!
! !DEFINED PARAMETERS:
! none

! !LOCAL VARIABLES:
  include 'wrf_io_flags.h'

! !SEE ALSO:
!  none

! !REVISION HISTORY:
! 05/28/03 - Dan Schaffer - Initial version
!EOP ___________________________________________________________________

  if (WRF_IO_ERROR ( &
    FIRST_OPEN(WRF_IO_API_OBJECT), &
    "Read occurred before open", ISTATUS)) return

  if (WRF_IO_ERROR ( &
    FIELD_TYPE /= WRF_REAL, &
    "Implementation currently only supports field type WRF_REAL", &
    ISTATUS)) return

  if (PIPE_COMMITTED(WRF_IO_API_OBJECT%EXCHANGE_PIPES(DATA_HANDLE))) then

    select case (FIELD_TYPE)

      case (WRF_REAL)

        call RECV_DATA_REAL_2D ( &
          WRF_IO_API_OBJECT, &
          DATA_HANDLE, &
          STAGGER, &
          PATCH_START, &
          PATCH_END, &
          MEMORY_START, &
          MEMORY_END, &
          VAR_NAME, &
          FIELD, &
          ISTATUS )

      case default
        ASSERT(.false.)

    end select

  else
    call EXPAND_PIPE ( &
      WRF_IO_API_OBJECT%EXCHANGE_PIPES(DATA_HANDLE), &
      VAR_NAME, &
      FIELD_TYPE, &
      STAGGER, &
      DOMAIN_START, &
      DOMAIN_END, &
      PATCH_START, &
      PATCH_END, &
      .true., &      ! A receiver
      ISTATUS )
  end if
  
end subroutine EXT_MCT_READ_FIELD

!BOP -------------------------------------------------------------------
!
! !ROUTINE: EXT_MCT_OPEN_FOR_WRITE_BEGIN
!
! !INTERFACE:

subroutine EXT_MCT_OPEN_FOR_WRITE_BEGIN ( &
  HIS_COMP_NAME, &
  COMM_GLOBAL, &
  COMM_LOCAL, &
  SYS_DEP_INFO, &
  DATA_HANDLE, &
  ISTATUS &
  )

!
! !USES:

  use MCT_WRF_IO_MODULE, only : EXT_MCT_OPEN_BEGIN

  use MCT_WRF_IO_UTILITIES_MODULE, only : &
    WRF_IO_ERROR

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
! Same as EXT_MCT_OPEN_FOR_READ_BEGIN except for a sender of
! re-gridded information.

!
! !DEFINED PARAMETERS:
! none

! !LOCAL VARIABLES:
! none

! !SEE ALSO:
!  none

! !REVISION HISTORY:
! 05/28/03 - Dan Schaffer - Initial version
!EOP ___________________________________________________________________

call EXT_MCT_OPEN_BEGIN ( &
  HIS_COMP_NAME, &
  COMM_GLOBAL, &
  COMM_LOCAL, &
  SYS_DEP_INFO, &
  DATA_HANDLE, &
  ISTATUS &
  )

end subroutine EXT_MCT_OPEN_FOR_WRITE_BEGIN

!BOP -------------------------------------------------------------------
!
! !ROUTINE: EXT_MCT_OPEN_FOR_WRITE_COMMIT
!
! !INTERFACE:

subroutine EXT_MCT_OPEN_FOR_WRITE_COMMIT ( &
  DATA_HANDLE, &
  ISTATUS )

!
! !USES:
!  none

  use MCT_WRF_IO_MODULE, only : WRF_IO_API_OBJECT

  use MCT_WRF_IO_PIPE_MODULE, only : COMMIT_PIPE => COMMIT

  implicit none

! !INPUT PARAMETERS:
!
  integer,          intent(in   ) :: DATA_HANDLE
  integer,          intent(inout) :: ISTATUS

! !DESCRIPTION:  
! Analogous to EXT_MCT_OPEN_FOR_READ_COMMIT

!
! !DEFINED PARAMETERS:
! none

! !LOCAL VARIABLES:
! none

! !SEE ALSO:
!  none

! !REVISION HISTORY:
! 05/28/03 - Dan Schaffer - Initial version
!EOP ___________________________________________________________________
  ISTATUS = 0
  call COMMIT_PIPE(WRF_IO_API_OBJECT%EXCHANGE_PIPES(DATA_HANDLE))

end subroutine EXT_MCT_OPEN_FOR_WRITE_COMMIT

!BOP -------------------------------------------------------------------
!
! !ROUTINE: EXT_MCT_OPEN_FOR_WRITE         
!
! !INTERFACE:

 subroutine EXT_MCT_OPEN_FOR_WRITE ( & 
  HIS_COMP_NAME, &
  COMM_GLOBAL, &
  COMM_LOCAL, &
  SYS_DEP_INFO, &
  DATA_HANDLE, &
  ISTATUS &
  )

!
! !USES:
!
  use MCT_WRF_IO_UTILITIES_MODULE, only : &
    WRF_IO_ERROR

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
! Although not currently implemented, this routine would eventually
! provide a means to open a data stream for immediate send
! of re-gridded data.  Under
! this scenario, there would be no "training" as discussed in 
! the API documentation.  For this implementation, that means
! there would be no aggregation of inter-component communication; 
! each boundary array would be sent by itself.
!
! !DEFINED PARAMETERS:
! none

! !LOCAL VARIABLES: 
! none

! !SEE ALSO:
!  none

! !REVISION HISTORY:
! 06/02/03 - Dan Schaffer - Initial version
!EOP ___________________________________________________________________

  if ( WRF_IO_ERROR  &
    ( .true., "EXT_MCT_OPEN_FOR_WRITE routine stubbed for now", &
    ISTATUS )) return

 end subroutine EXT_MCT_OPEN_FOR_WRITE        

!BOP -------------------------------------------------------------------
!
! !ROUTINE: EXT_MCT_WRITE_FIELD
!
! !INTERFACE:

subroutine EXT_MCT_WRITE_FIELD ( &
  DATA_HANDLE, &
  DATE_STR, &
  VAR_NAME, &
  FIELD, &
  FIELD_TYPE, &
  COMM_GLOBAL, &
  COMM_LOCAL, &
  DOMAIN_DESC, &
  MEMORY_ORDER, &
  STAGGER, &
  DIM_NAMES, &
  DOMAIN_START, &
  DOMAIN_END, &
  MEMORY_START, &
  MEMORY_END, &
  PATCH_START, &
  PATCH_END, &
  ISTATUS )

!
! !USES:
  use MCT_WRF_IO_MODULE, only : WRF_IO_API_OBJECT
  use MCT_WRF_IO_MODULE, only : SEND_DATA_REAL_2D
  use MCT_WRF_IO_MODULE, only : FIRST_OPEN

  use MCT_WRF_IO_UTILITIES_MODULE, only : Assert
  use MCT_WRF_IO_UTILITIES_MODULE, only : &
    WRF_IO_ERROR

  use MCT_WRF_IO_PIPE_MODULE, only : &
    PIPE_COMMITTED => COMMITTED, &
    EXPAND_PIPE => EXPAND

  implicit none

! !INPUT PARAMETERS:
!
  integer,          intent(in   ) :: DATA_HANDLE
  character(LEN=*), intent(in   ) :: DATE_STR        ! Not used
  character(LEN=*), intent(in   ) :: VAR_NAME

  integer,          intent(in   ) :: FIELD_TYPE     
  integer,          intent(in   ) :: COMM_GLOBAL   
  integer,          intent(in   ) :: COMM_LOCAL   
  integer,          intent(in   ) :: DOMAIN_DESC     ! Not used
  character(LEN=*), intent(in   ) :: MEMORY_ORDER    ! Not used
  character(LEN=*), intent(in   ) :: STAGGER
  character(LEN=*), intent(in   ) :: DIM_NAMES       ! Not used
  integer,          intent(in   ) :: DOMAIN_START(len_trim(MEMORY_ORDER))
  integer,          intent(in   ) :: DOMAIN_END  (len_trim(MEMORY_ORDER))
  integer,          intent(in   ) :: MEMORY_START(len_trim(MEMORY_ORDER))
  integer,          intent(in   ) :: MEMORY_END  (len_trim(MEMORY_ORDER))
  integer,          intent(in   ) :: PATCH_START (len_trim(MEMORY_ORDER))
  integer,          intent(in   ) :: PATCH_END   (len_trim(MEMORY_ORDER))
  real,             intent(in   ) :: FIELD( &
    MEMORY_START(1):MEMORY_END(1),MEMORY_START(2):MEMORY_END(2))      

! !INPUT/OUTPUT PARAMETERS:
  integer,          intent(inout) :: ISTATUS

! !DESCRIPTION:  
! This routine does one of several things depending on the conditions
! upon which it is called.  If it is called prior to
! a call to EXT_MCT_OPEN_FOR_WRITE_BEGIN then there is an error condition.
! If it is called after a call to EXT_MCT_OPEN_FOR_WRITE_BEGIN but before
! a call to EXT_MCT_OPEN_FOR_WRITE_COMMIT then
! the pipe between this component and the sending component is expanded
! to include the variable whose name is VAR_NAME.  If 
! EXT_MCT_OPEN_FOR_WRITE_COMMIT has been called then SEND_DATA_REAL_2D
! is called to send the data into the pipe.

!
! !DEFINED PARAMETERS:
! none

! !LOCAL VARIABLES:
  include 'wrf_io_flags.h'

! !SEE ALSO:
!  none

! !REVISION HISTORY:
! 05/28/03 - Dan Schaffer - Initial version
!EOP ___________________________________________________________________

  if (WRF_IO_ERROR ( &
    FIRST_OPEN(WRF_IO_API_OBJECT), &
    "Write occurred before open", ISTATUS)) return

  if (WRF_IO_ERROR ( &
    FIELD_TYPE /= WRF_REAL, &
    "Implementation currently only supports field type WRF_REAL", &
    ISTATUS)) return

  if (PIPE_COMMITTED(WRF_IO_API_OBJECT%EXCHANGE_PIPES(DATA_HANDLE))) then

    select case (FIELD_TYPE)

      case (WRF_REAL)

        call SEND_DATA_REAL_2D ( &
          WRF_IO_API_OBJECT, &
          DATA_HANDLE, &
          STAGGER, &
          PATCH_START, &
          PATCH_END, &
          MEMORY_START, &
          MEMORY_END, &
          VAR_NAME, &
          FIELD, &
          ISTATUS )

      case default
        ASSERT(.false.)

    end select

  else
    call EXPAND_PIPE ( &
      WRF_IO_API_OBJECT%EXCHANGE_PIPES(DATA_HANDLE), &
      VAR_NAME, &
      FIELD_TYPE, &
      STAGGER, &
      DOMAIN_START, &
      DOMAIN_END, &
      PATCH_START, &
      PATCH_END, &
      .false., &      ! Not a receiver
      ISTATUS )
  end if
  
end subroutine EXT_MCT_WRITE_FIELD

!BOP -------------------------------------------------------------------
!
! !ROUTINE: EXT_MCT_INQUIRY         
!
! !INTERFACE:

 subroutine EXT_MCT_INQUIRY ( &
   INQUIRY, &
   RESULT, &
   ISTATUS)
!
! !USES:
!
  use MCT_WRF_IO_MODULE, only : EXECUTION_MODE

  use MCT_WRF_IO_UTILITIES_MODULE, only : &
    WRF_IO_ERROR

  implicit none

! !INPUT PARAMETERS: 
!
  character(len=*) , intent(in)  :: INQUIRY

! !OUTPUT PARAMETERS: 
!
  character(len=*) , intent(out) :: RESULT
  integer          , intent(out) :: ISTATUS

! !DESCRIPTION:  Provide the user with information about this 
! implementation of the WRF I/O API.
!
! !DEFINED PARAMETERS:
! none

! !LOCAL VARIABLES: 
! none

! !SEE ALSO:
!  none

! !REVISION HISTORY:
! 06/02/03 - Dan Schaffer - Initial version
!EOP ___________________________________________________________________

   select case (trim(INQUIRY))

     case ("RANDOM_READ")
       RESULT = 'NO'

     case ("RANDOM_WRITE")
       RESULT = 'NO'

     case ("SERIAL_READ")
       RESULT = 'YES'

     case ("SERIAL_WRITE")
       RESULT = 'YES'

     case ("OPEN_READ")
       RESULT = 'NO'

     case ("OPEN_WRITE")
       RESULT = 'NO'

     case ("OPEN_COMMIT_READ")
       RESULT = 'YES'

     case ("OPEN_COMMIT_WRITE")
       RESULT = 'YES'

     case ("MEDIUM")
       RESULT = 'COUPLED'

     case ("PARALLEL_IO")
       RESULT = 'YES'

     case ("SELF_DESCRIBING")
       RESULT = 'NO'

     case ("SUPPORT_METADATA")
       RESULT = 'NO'

     case ("SUPPORT_3D_FIELDS")
       RESULT = 'NO'

     case ("EXECUTION_MODE")
       RESULT = EXECUTION_MODE

     case default
       if (WRF_IO_ERROR(.true., "Invalid Inquiry", ISTATUS)) return
   end select
 end subroutine EXT_MCT_INQUIRY        

!BOP -------------------------------------------------------------------
!
! !ROUTINE: EXT_MCT_IOEXIT
!
! !INTERFACE:

subroutine EXT_MCT_IOEXIT ( &
  ISTATUS &
  )
!
! !USES:
! none

  implicit none

! !INPUT PARAMETERS:
!
! none

! !OUTPUT PARAMETERS:
!
  integer         , intent(  out) :: ISTATUS

! !DESCRIPTION:  Finalize the MCT coupler implementation of the WRF
! I/O API.

!
! !DEFINED PARAMETERS:
! none

! !LOCAL VARIABLES:
! none

! !SEE ALSO:

! !REVISION HISTORY:
! 06/12/03 - Dan Schaffer - Initial version
!EOP ___________________________________________________________________

  ISTATUS = 0

end subroutine EXT_MCT_IOEXIT

