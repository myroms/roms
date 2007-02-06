!BOP -------------------------------------------------------------------
!
! !MODULE: MCT_WRF_IO_INTERPOLATOR_MODULE - Interpolator class for WRF I/O API
!          
! !INTERFACE:

 module MCT_WRF_IO_INTERPOLATOR_MODULE
!
! !DESCRIPTION:
!   A class that provides the capability to interpolate data stored
! in MCT attribute vectors between two groups (components).
! 
! Here is a typical scenario that will invoke this class and objects.
! We will see how the decompositions are constructed and then used
! to create routers.  
! .bp
! Suppose we have model 1 that reads data from model 2 and writes data
! to model 2.  Model 2 writes data to model 1 and reads data from model 1.
! .bp
! Training phase.  Note that the tags referenced below are all identical
! to each other in the concurrent case (they are not needed).
! .bp
! A. Model 1 read field
! .bp
! The receiver, model 1, uses the sparse matrix to construct a receiver
! decomposition of the sender's (model 2) grid.  The component id of this
! decomposition is 1.  It is sent to model 2 with tag1.
! .bp
! B. Model 1 write field
! .bp
! The sender, model 1, constructs decomposition of its own grid.  The 
! component id is 1.  This decomposition is sent to model 2 with tag2.
! .bp
! C. Model 2 write field
! .bp
! The sender, model 2, constructs a decomposition of its own grid.  The
! component id is 2.  This decomposition is sent to model 1 with tag3.
! .bp
! D. Model 2 read field
! .bp
! The receiver, model 2, uses the sparse matrix to construct a receiver
! decomposition of the sender's (model 1) grid.  The component id of
! this decomposition is 2.  It is sent to model 1 with tag4.
! .bp
! Commit phase for the SEQUENTIAL case where model 1 runs before model 2.
! The key point is that the first model 1 read_field is skipped.
! .bp
! E. Model 1 write field
! .bp
! Model 1 receives the decomposition sent by model 2 in step D (tag4).
! This decomposition is the second of two arguments to the MCT_ROUTER_INIT
! routine.  The first is its own decomposition of its own grid constructed
! in step B.
! .bp
! F. Model 2 write field
! .bp
! Model 2 receives the decomposition sent by model 1 in step A (tag1). 
! This decomposition is the second argument to MCT_ROUTER_INIT.  The
! first is its own decomposition of its own grid determined in step C.
! .bp
! G. Model 2 read field
! .bp
! Model 2 receives the decomposition sent by model 1 in step B (tag2).
! This is the second argument to MCT_ROUTER_INIT.  The first is its
! decomposition of the sender's grid constructed in step D.
! .bp
! H. Model 1 read field
! .bp
! Model 1 receives the decomposition sent by model 2 in step C (tag3).
! This is the second argument to MCT_ROUTER_INIT.  The first is its 
! decomposition of the sender's grid constructed in step A.
! .bp
! Commit phase for the CONCURRENT case. The two main differences are:
! .bp
! 1. Model 1 read field is called before Model 1 write field.
! .bp
! 2. The tags are all set to a single value since they are not needed.
! .bp
! A second set of tags are used for the actual interpolation.
! Again, for concurrent, they are
! not used.  For sequential, each component keeps track of which
! tag it is at.  A new COMM_TAG is assigned to each INTERPOLATOR
! constructed.  The COMM_TAGS match between the components.
! So if an interpolator is constructed for component 1 where it
! writes to component 2 and a matching interpolator is constructed
! where component 2 reads from component 1 then both will have
! COMM_TAG tag1.

! !USES:

  use MCT_WRF_IO_UTILITIES_MODULE, only : ERROR_CHECK
  use MCT_WRF_IO_UTILITIES_MODULE, only : MAX_STR

  use m_AttrVect, only : MCT_ATTR_VECTOR => AttrVect

  use M_GlobalSegmap, only : GlobalSegMap

  use m_Router, only : MCT_ROUTER => Router

  use m_SparseMatrix, only : MCT_SPARSE_MATRIX => SparseMatrix

  use MCT_WRF_IO_UTILITIES_MODULE, only : Assert

  implicit none

  private	! except

! !PUBLIC TYPES:
!
  public :: T_CLASS

! !PUBLIC MEMBER FUNCTIONS:
!
  public :: INITIALIZE  ! Initializes class variables
  public :: CONSTRUCTOR
  public :: DESTRUCTOR
  public :: EXTEND_SRC
  public :: EXTEND_DST
  public :: INTERPOLATE_SRC
  public :: INTERPOLATE_DST
  public :: MATCHES

! !DEFINED PARAMETERS:
  integer, parameter :: UNDEFINED_TAG = -999
  integer, parameter :: BASE_DECOMP_TAG = 100
  integer, parameter :: BASE_COMM_TAG = 1000
  integer, parameter :: TO = 1
  integer, parameter :: FROM = 2

! !SEE ALSO:
!  none

! !REVISION HISTORY:
! 05/21/03 - Dan Schaffer - Initial version

!EOP ___________________________________________________________________

!  CLASS VARIABLES:

  ! Counter of decomp tags
  integer, save :: CURRENT_DECOMP_TAG

  ! Array over sending components, receiving components, and
  ! directions (TO or FROM) of tags
  ! used to send/receive decompositions for setting up the routers.
  integer, save, pointer :: DECOMP_TAGS(:,:,:)

  ! Array over components of tags used to send/receive the data
  integer, save, pointer :: COMM_TAGS(:)

  ! Indicates if objects are for sequential or concurrent coupling
  logical, save :: SEQUENTIAL

  type T_CLASS
!  INSTANCE VARIABLES:

    ! The source attribute vector decomposed on the source grid
    type (MCT_ATTR_VECTOR) :: SRC_ATTR_VECT_SRC_GRID

    ! The destination attribute vector decomposed on the source grid
    type (MCT_ATTR_VECTOR) :: DST_ATTR_VECT_SRC_GRID

    ! The destination attribute vector decomposed on the destination grid
    type (MCT_ATTR_VECTOR) :: DST_ATTR_VECT_DST_GRID

    ! The MCT object containing communication patterns required to
    ! re-grid an attribute vector
    type (MCT_ROUTER) :: ROUTER

    ! The sparse matrix used to construct the re-gridding operation
    type (MCT_SPARSE_MATRIX) :: SPARSE_MATRIX

    character(LEN=MAX_STR) :: ID  ! Used to locate this interpolator from a list

    integer :: CURRENT_VARIABLE ! Keeps track of what variables have been
                                ! Packed into or unpacked from the attribute 
				! vector

    ! The source grid decomposed using the source component decomposition
    type (GlobalSegMap) :: SRC_DECOMP_OF_SRC_GRID

    ! The destination grid decomposed using the source component decomposition
    type (GlobalSegMap) :: DST_DECOMP_OF_SRC_GRID

    ! The destination grid decomposed using the destination component 
    ! decomposition
    type (GlobalSegMap) :: DST_DECOMP_OF_DST_GRID

    integer :: MY_COMM

    logical :: CONSTRUCTION_FINISHED

    integer, pointer :: SEND_DECOMP_REQUEST(:)
    integer, pointer :: RECV_DECOMP_REQUEST(:)

    integer :: MY_COMP_ID  ! Id of this component

    integer :: HIS_COMP_ID ! Id of the component with which this component is
                           ! communicating.

    integer :: TAG ! Tag of the interpolator that is used to match up
                   ! messages with the other component

  end type

  interface INTERPOLATE_SRC
    module procedure INTERPOLATE_SRC_REAL_2D
  end interface

  interface INTERPOLATE_DST
    module procedure INTERPOLATE_DST_REAL_2D
  end interface

#include "macros.h"

  include 'mpif.h'

  contains

!BOP -------------------------------------------------------------------
!
! !IROUTINE: INITIALIZE         
!
! !INTERFACE:

subroutine INITIALIZE ( &
  NUM_COMPONENTS, &
  SEQ)
!
! !USES:
!
  implicit none

! !INPUT PARAMETERS: 
!
  integer          , intent(in)  :: NUM_COMPONENTS
  logical          , intent(in)  :: SEQ

! !OUTPUT PARAMETERS: 
!
! none

! !DESCRIPTION:  Initializes class variables
!
! !DEFINED PARAMETERS:
! none

! !LOCAL VARIABLES: 
  integer :: ISTATUS

! !SEE ALSO:
!  none

! !REVISION HISTORY:
! 08/28/03 - Dan Schaffer - Initial version
!EOP ___________________________________________________________________

#ifdef DEBUG
  print *, 'in initialize'
#endif

  allocate(DECOMP_TAGS(NUM_COMPONENTS, NUM_COMPONENTS, 2), stat=ISTATUS)
  ASSERT(ISTATUS == 0)

  DECOMP_TAGS = UNDEFINED_TAG

  CURRENT_DECOMP_TAG = BASE_DECOMP_TAG

  allocate(COMM_TAGS(NUM_COMPONENTS), stat=ISTATUS)
  ASSERT(ISTATUS == 0)

  COMM_TAGS = BASE_COMM_TAG

  SEQUENTIAL = SEQ

end subroutine INITIALIZE        

!BOP -------------------------------------------------------------------
!
! !IROUTINE: CONSTRUCTOR      - Construct the object
!
! !INTERFACE:

 subroutine CONSTRUCTOR ( &
   THIS, &
   VAR_NAME, &
   VAR_TYPE, &
   STARTS, &
   LENGTHS, &
   COMM, &
   MY_ID, &
   HIS_ID, &
   SPARSE_MATRIX_FNAME, &
   RECEIVER, &
   ISTATUS)
!
! !USES:
!
  use MCT_WRF_IO_UTILITIES_MODULE, only : &
    MY_PE

  use M_GlobalSegMap, only : &
    MCT_DECOMP_INIT => init, &
    MCT_DECOMP_DESTROY => clean, &
    MCT_DECOMP_SIZE => LSIZE

  use m_SparseMatrix, only : &
    SparseMatrix, &
    MCT_SM_INIT => INIT, &
    MCT_SM_DESTROY => CLEAN, &
    MCT_SM_IMPORT_ROWS => importGlobalRowIndices, &
    MCT_SM_IMPORT_COLS => importGlobalColumnIndices, &
    MCT_SM_IMPORT_ELTS => importMatrixElements

  use m_SparseMatrixComms, only : &
    MCT_DISTRIBUTE_SM => ScatterByRow

  use m_GlobalSegMapComms, only : MCT_ISEND_DECOMP => ISEND

  use m_SparseMatrixToMaps, only : &
    MCT_SM_TO_DECOMP => SparseMatrixToXGlobalSegMap

  use m_GlobalToLocal, only : &
    MCT_SM_GLOBAL_TO_LOCAL => GlobalToLocalMatrix

  use m_AttrVect, only : &
    MCT_INIT_ATTR_VECTOR => INIT

  implicit none

! !INPUT PARAMETERS: 
!
  character(LEN=*),  intent(in)  :: VAR_NAME    ! Name of the 1st variable
                                                ! to be included in the 
						! interpolator
  integer,           intent(in)  :: VAR_TYPE    ! Type/kind of var to interp
  integer,           intent(in)  :: STARTS(:)   ! Starts of my PE's vectors
  integer,           intent(in)  :: LENGTHS(:)  ! Lengths of my PE's vectors
  integer,           intent(in)  :: COMM        ! Communicator of my group
  integer,           intent(in)  :: MY_ID       ! Identifier of my component
  integer,           intent(in)  :: HIS_ID      ! Identifier of his component
  character(LEN=*),  intent(in) &
                         :: SPARSE_MATRIX_FNAME ! Name of file containing
    						! sparse matrix weights/indices
  logical,           intent(in)  :: RECEIVER    ! True if my proc is 
  						! receiving data.

! !INPUT/OUTPUT PARAMETERS: 
  integer,           intent(inout) :: ISTATUS

! !OUTPUT PARAMETERS: 
!
  type (T_CLASS)   , intent(out) :: THIS

! !DESCRIPTION:  This routine constructs an interpolator.
! For senders, a sender decomposition of the sender grid is constructed.
! For receivers, receiver decompositions of the receiver and sender
! grids are computed.
! An exchange of these decompositions between the receiver and sender
! is initiated.
! Also, attribute vectors on the sender grid and receiver grids are created.
! In addition, a sparse matrix is constructed.
!
! !DEFINED PARAMETERS:
! none

! !LOCAL VARIABLES: 
  integer :: SRC_GRID_SIZE
  integer :: DST_GRID_SIZE
  integer :: NUM_SM_ELEMENTS
  integer, pointer :: SM_ROWS(:)
  integer, pointer :: SM_COLS(:)
  real(kind=8), pointer :: SM_WEIGHTS(:)
  real(kind=4), pointer :: SM_WEIGHTS_R4(:)
  integer :: N
  integer :: NCID
  integer :: DIMID
  integer :: VARID

  type (SparseMatrix) :: ROOT_SM

  include 'wrf_io_flags.h'
  include 'netcdf.inc'

! !SEE ALSO:
! none

! !REVISION HISTORY:
! 05/21/03 - Dan Schaffer - Initial version
!EOP ___________________________________________________________________

  THIS%CONSTRUCTION_FINISHED = .false.
  THIS%MY_COMM = COMM
  THIS%CURRENT_VARIABLE = 0
  THIS%ID = SPARSE_MATRIX_FNAME
  THIS%MY_COMP_ID = MY_ID
  THIS%HIS_COMP_ID = HIS_ID
  THIS%TAG = COMM_TAGS(THIS%MY_COMP_ID)

#ifdef DEBUG
  print *, 'got to interpolator constructor new tag is ', &
    this%my_comp_id, this%tag
#endif

  if (SEQUENTIAL) then
    COMM_TAGS(THIS%MY_COMP_ID) = COMM_TAGS(THIS%MY_COMP_ID) + 4
  end if

  if (.not. RECEIVER) then

    call MCT_DECOMP_INIT ( &
      THIS%SRC_DECOMP_OF_SRC_GRID, &
      STARTS, &
      LENGTHS, &
      0, &
      COMM, &
      MY_ID )

    select case (VAR_TYPE)
      case (WRF_REAL)
        call MCT_INIT_ATTR_VECTOR ( &
	  THIS%SRC_ATTR_VECT_SRC_GRID, &
	  RLIST=VAR_NAME, &
	  LSIZE=MCT_DECOMP_SIZE(THIS%SRC_DECOMP_OF_SRC_GRID, COMM))

      case default
        ASSERT(.false.)
    end select

    if (MY_PE(COMM) == 0) then

      if (SEQUENTIAL) then
        if (DECOMP_TAGS(THIS%MY_COMP_ID, THIS%HIS_COMP_ID, TO) &
          == UNDEFINED_TAG) then
  
          DECOMP_TAGS(THIS%MY_COMP_ID, THIS%HIS_COMP_ID, TO) = &
	    CURRENT_DECOMP_TAG
        end if
      else
        DECOMP_TAGS(THIS%MY_COMP_ID, THIS%HIS_COMP_ID, TO) = BASE_DECOMP_TAG
      end if

#ifdef DEBUG
      print *, 'got to interpolate constructor 1, ids, TAG, COMP_ID', &
	THIS%MY_COMP_ID, THIS%HIS_COMP_ID, &
        DECOMP_TAGS(THIS%MY_COMP_ID, THIS%HIS_COMP_ID, TO), &
	THIS%SRC_DECOMP_OF_SRC_GRID%COMP_ID
#endif


      call MCT_ISEND_DECOMP ( &
        THIS%SRC_DECOMP_OF_SRC_GRID, &
        THIS%HIS_COMP_ID, &
	DECOMP_TAGS(THIS%MY_COMP_ID, THIS%HIS_COMP_ID, TO), &
        THIS%SEND_DECOMP_REQUEST, &
	ISTATUS)

      ASSERT(ISTATUS == 0)

      if (SEQUENTIAL) then
        CURRENT_DECOMP_TAG = CURRENT_DECOMP_TAG + 4
      end if

      deallocate(THIS%SEND_DECOMP_REQUEST, stat=ISTATUS)

      ASSERT(ISTATUS == 0)
    end if


  else ! receiver

    if (MY_PE(COMM) == 0) then

      ISTATUS = NF_OPEN(SPARSE_MATRIX_FNAME, NF_NOWRITE, NCID )

      call ERROR_CHECK(ISTATUS == 0, &
        'Sparse matrix file bad or not found: ' // SPARSE_MATRIX_FNAME)


      ISTATUS = NF_INQ_DIMID(NCID, "num_links", DIMID)

      call ERROR_CHECK(ISTATUS == 0, &
        'Invalid sparse matrix file num_links dim: ' // SPARSE_MATRIX_FNAME)

      ISTATUS = NF_INQ_DIMLEN(NCID, DIMID, NUM_SM_ELEMENTS)

      call ERROR_CHECK(ISTATUS == 0, &
        'Invalid sparse matrix file num_links dim: ' // SPARSE_MATRIX_FNAME)


      ISTATUS = NF_INQ_DIMID(NCID, "src_grid_size", DIMID)

      call ERROR_CHECK(ISTATUS == 0, &
        'Invalid sparse matrix file src_grid_size: ' // SPARSE_MATRIX_FNAME)

      ISTATUS = NF_INQ_DIMLEN(NCID, DIMID, SRC_GRID_SIZE)

      call ERROR_CHECK(ISTATUS == 0, &
        'Invalid sparse matrix file src_grid_size: ' // SPARSE_MATRIX_FNAME)


      ISTATUS = NF_INQ_DIMID(NCID, "dst_grid_size", DIMID)

      call ERROR_CHECK(ISTATUS == 0, &
        'Invalid sparse matrix file dst_grid_size: ' // SPARSE_MATRIX_FNAME)

      ISTATUS = NF_INQ_DIMLEN(NCID, DIMID, DST_GRID_SIZE)

      call ERROR_CHECK(ISTATUS == 0, &
        'Invalid sparse matrix file dst_grid_size: ' // SPARSE_MATRIX_FNAME)


      allocate(SM_ROWS(NUM_SM_ELEMENTS), stat=ISTATUS)
      ASSERT(ISTATUS==0)

      allocate(SM_COLS(NUM_SM_ELEMENTS), stat=ISTATUS)
      ASSERT(ISTATUS==0)

      allocate(SM_WEIGHTS(NUM_SM_ELEMENTS), stat=ISTATUS)
      ASSERT(ISTATUS==0)

      allocate(SM_WEIGHTS_R4(NUM_SM_ELEMENTS), stat=ISTATUS)
      ASSERT(ISTATUS==0)


      ISTATUS = NF_INQ_VARID(NCID, "src_address", VARID)

      call ERROR_CHECK(ISTATUS == 0, &
        'Invalid sparse matrix file src_address: ' // SPARSE_MATRIX_FNAME)

      ISTATUS = NF_GET_VAR_INT(NCID, VARID, SM_COLS)

      call ERROR_CHECK(ISTATUS == 0, &
        'Invalid sparse matrix file src_address: ' // SPARSE_MATRIX_FNAME)


      ISTATUS = NF_INQ_VARID(NCID, "dst_address", VARID)

      call ERROR_CHECK(ISTATUS == 0, &
        'Invalid sparse matrix file dst_address: ' // SPARSE_MATRIX_FNAME)

      ISTATUS = NF_GET_VAR_INT(NCID, VARID, SM_ROWS)

      call ERROR_CHECK(ISTATUS == 0, &
        'Invalid sparse matrix file dst_address: ' // SPARSE_MATRIX_FNAME)


      ISTATUS = NF_INQ_VARID(NCID, "remap_matrix", VARID)

      call ERROR_CHECK(ISTATUS == 0, &
        'Invalid sparse matrix file remap_matrix: ' // SPARSE_MATRIX_FNAME)

      ISTATUS = NF_GET_VARS_DOUBLE(NCID, VARID, &
        (/1,1/), &  ! Start
        (/1, NUM_SM_ELEMENTS/), & ! Count
	(/1,1/), & ! stride
	SM_WEIGHTS)

      call ERROR_CHECK(ISTATUS == 0, &
        'Invalid sparse matrix file remap_matrix: ' // SPARSE_MATRIX_FNAME)

      SM_WEIGHTS_R4 = SM_WEIGHTS

      call MCT_SM_INIT(ROOT_SM, DST_GRID_SIZE, SRC_GRID_SIZE, NUM_SM_ELEMENTS)
      call MCT_SM_IMPORT_ROWS(ROOT_SM, SM_ROWS, size(SM_ROWS))
      call MCT_SM_IMPORT_COLS(ROOT_SM, SM_COLS, size(SM_COLS))
      call MCT_SM_IMPORT_ELTS(ROOT_SM, SM_WEIGHTS_R4, size(SM_WEIGHTS_R4))

      deallocate(SM_ROWS, stat=ISTATUS)
      ASSERT(ISTATUS==0)

      deallocate(SM_COLS, stat=ISTATUS)
      ASSERT(ISTATUS==0)

      deallocate(SM_WEIGHTS, stat=ISTATUS)
      ASSERT(ISTATUS==0)

      deallocate(SM_WEIGHTS_R4, stat=ISTATUS)
      ASSERT(ISTATUS==0)

    end if ! MY_PE(COMM) == 0

    call MCT_DECOMP_INIT ( &
      THIS%DST_DECOMP_OF_DST_GRID, &
      STARTS, &
      LENGTHS, &
      0, &
      COMM, &
      MY_ID )

    call MCT_DISTRIBUTE_SM ( &
      THIS%DST_DECOMP_OF_DST_GRID, &
      ROOT_SM, &
      THIS%SPARSE_MATRIX, &
      0, &
      COMM, &
      ISTATUS )

    ASSERT(ISTATUS == 0)
    call MCT_SM_TO_DECOMP ( &
      THIS%SPARSE_MATRIX, &
      THIS%DST_DECOMP_OF_SRC_GRID, &
      0, &
      COMM, &
      MY_ID )

#ifdef DEBUG
    print *, 'Got to constructor receiver 1, DST_DECOMP_OF_SRC_GRID id ', &
      THIS%DST_DECOMP_OF_SRC_GRID%comp_id, MY_ID
#endif

    ! Compute the process local rows
    call MCT_SM_GLOBAL_TO_LOCAL ( &
      THIS%SPARSE_MATRIX, &
      THIS%DST_DECOMP_OF_DST_GRID, &
      'row', &
      COMM )

    ! Compute the process local columns
    call MCT_SM_GLOBAL_TO_LOCAL ( &
      THIS%SPARSE_MATRIX, &
      THIS%DST_DECOMP_OF_SRC_GRID, &
      'column', &
      COMM )

    select case (VAR_TYPE)
      case (WRF_REAL)
        call MCT_INIT_ATTR_VECTOR ( &
	  THIS%DST_ATTR_VECT_SRC_GRID, &
	  RLIST=VAR_NAME, &
	  LSIZE=MCT_DECOMP_SIZE(THIS%DST_DECOMP_OF_SRC_GRID, COMM))

        call MCT_INIT_ATTR_VECTOR ( &
	  THIS%DST_ATTR_VECT_DST_GRID, &
	  RLIST=VAR_NAME, &
	  LSIZE=MCT_DECOMP_SIZE(THIS%DST_DECOMP_OF_DST_GRID, COMM))

      case default
        ASSERT(.false.)
    end select

    if (MY_PE(COMM) == 0) then

      if (SEQUENTIAL) then
        if (DECOMP_TAGS(THIS%MY_COMP_ID, THIS%HIS_COMP_ID, FROM) &
            == UNDEFINED_TAG) then

          DECOMP_TAGS(THIS%MY_COMP_ID, THIS%HIS_COMP_ID, FROM) &
	     = CURRENT_DECOMP_TAG
        end if
      else
        DECOMP_TAGS(THIS%MY_COMP_ID, THIS%HIS_COMP_ID, FROM) &
	   = BASE_DECOMP_TAG
      end if

#ifdef DEBUG
      print *, 'got to interpolate constructor 3, ids, TAG, COMP_ID', &
	THIS%MY_COMP_ID, THIS%HIS_COMP_ID, &
        DECOMP_TAGS(THIS%MY_COMP_ID, THIS%HIS_COMP_ID, FROM), &
	THIS%DST_DECOMP_OF_SRC_GRID%COMP_ID
#endif

      call MCT_ISEND_DECOMP ( &
        THIS%DST_DECOMP_OF_SRC_GRID, &
        THIS%HIS_COMP_ID, &
	DECOMP_TAGS(THIS%MY_COMP_ID, THIS%HIS_COMP_ID, FROM), &
        THIS%SEND_DECOMP_REQUEST, &
	ISTATUS)

      ASSERT(ISTATUS == 0)

      if (SEQUENTIAL) then
        CURRENT_DECOMP_TAG = CURRENT_DECOMP_TAG + 4
      end if

      deallocate(THIS%SEND_DECOMP_REQUEST, stat=ISTATUS)

      ASSERT(ISTATUS == 0)

    end if

  end if

end subroutine CONSTRUCTOR      

!BOP -------------------------------------------------------------------
!
! !IROUTINE: DESTRUCTOR         
!
! !INTERFACE:

 subroutine DESTRUCTOR        (THIS)
!
! !USES:
!
  implicit none

! !INPUT PARAMETERS: 
!
! None

! !OUTPUT PARAMETERS: 
!
  type (T_CLASS)   , intent(out) :: THIS

! !DESCRIPTION:  Destroys the WRF Coupler I/O API Interpolator
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

  print *, 'MCT_WRF_IO_INTERPOLATOR_MODULE DESTRUCTOR routine stubbed for now'
 end subroutine DESTRUCTOR        

!BOP -------------------------------------------------------------------
!
! !IROUTINE: EXTEND_SRC         
!
! !INTERFACE:

 subroutine EXTEND_SRC (THIS, VAR_NAME, VAR_TYPE, ISTATUS)
!
! !USES:
!
  use m_AttrVect, only : &
    MCT_EXTEND_ATTR_VECT_REAL => appendRAttr

  implicit none

! !INPUT PARAMETERS: 
!
  character(LEN=*) , intent(in)  :: VAR_NAME    ! Variable name
  integer          , intent(in)  :: VAR_TYPE    ! Type/kind of var to interp

! !INPUT/OUTPUT PARAMETERS: 
!
  type (T_CLASS), intent(inout) :: THIS         ! Object
  integer       , intent(inout) :: ISTATUS

! !DESCRIPTION:  Add a variable to the SRC attribute vector associated
! with THIS.
!
! !DEFINED PARAMETERS:
! none

! !LOCAL VARIABLES: 
  include 'wrf_io_flags.h'

! !SEE ALSO:
!  none

! !REVISION HISTORY:
! 05/21/03 - Dan Schaffer - Initial version
!EOP ___________________________________________________________________

  select case (VAR_TYPE)
  
    case (WRF_REAL)

      call MCT_EXTEND_ATTR_VECT_REAL ( &
        THIS%SRC_ATTR_VECT_SRC_GRID, &
        VAR_NAME, &
        ISTATUS )

      case default
        ASSERT(.false.)
    end select

  ASSERT(ISTATUS == 0)

end subroutine EXTEND_SRC        

!BOP -------------------------------------------------------------------
!
! !IROUTINE: EXTEND_DST         
!
! !INTERFACE:

subroutine EXTEND_DST (THIS, VAR_NAME, VAR_TYPE, ISTATUS)
!
! !USES:
!
  use m_AttrVect, only : &
    MCT_EXTEND_ATTR_VECT_REAL => appendRAttr

  implicit none

! !INPUT PARAMETERS: 
!
  character(LEN=*) , intent(in)  :: VAR_NAME    ! Variable name
  integer          , intent(in)  :: VAR_TYPE    ! Type/kind of var to interp

! !INPUT/OUTPUT PARAMETERS: 
!
  type (T_CLASS), intent(inout) :: THIS         ! Object
  integer       , intent(inout) :: ISTATUS


! !DESCRIPTION:  Add a variable to the DST attribute vectors associated
! with THIS.
!
! !DEFINED PARAMETERS:
! none

! !LOCAL VARIABLES: 
  include 'wrf_io_flags.h'

! !SEE ALSO:
!  none

! !REVISION HISTORY:
! 05/21/03 - Dan Schaffer - Initial version
!EOP ___________________________________________________________________

  select case (VAR_TYPE)
  
    case (WRF_REAL)

      call MCT_EXTEND_ATTR_VECT_REAL ( &
        THIS%DST_ATTR_VECT_SRC_GRID, &
        VAR_NAME, &
        ISTATUS )
      ASSERT(ISTATUS == 0)

      call MCT_EXTEND_ATTR_VECT_REAL ( &
        THIS%DST_ATTR_VECT_DST_GRID, &
        VAR_NAME, &
        ISTATUS )
      ASSERT(ISTATUS == 0)

    case default
      ASSERT(.false.)
  end select

end subroutine EXTEND_DST        

!BOP -------------------------------------------------------------------
!
! !IROUTINE: INTERPOLATE_SRC_REAL_2D  - Interpolate the given src DATA
!
! !INTERFACE:

subroutine INTERPOLATE_SRC_REAL_2D (THIS, VAR_NAME, MSG_ID, DATA, ISTATUS)
!
! !USES:
!

  use m_String, only : &
    String, &
    MPEU_STRING_TO_CHAR => toChar

  use MCT_WRF_IO_UTILITIES_MODULE, only : WRF_IO_ERROR

  use MCT_WRF_IO_UTILITIES_MODULE, only : MY_PE

  use m_AttrVect, only : &
    MCT_PUT_ATTR_VECT_DATA => importRAttr, &
    MCT_ATTR_VECT_LENGTH => nRAttr, &
    MCT_ATTR_NAME => getRList
  

  use m_GlobalSegMapComms, only : MCT_RECV_DECOMP => RECV

  use m_GlobalSegMapComms, only : MCT_BCAST_DECOMP => BCAST

  use m_Router, only : MCT_ROUTER_INIT => INIT

  use m_Transfer, only : MCT_ATTR_VECT_ISEND => ISEND
  use m_Transfer, only : MCT_ATTR_VECT_WAITSEND => WAITSEND

  implicit none

! !INPUT PARAMETERS: 
!
  character(LEN=*) , intent(in)  :: VAR_NAME    ! Variable name
  integer          , intent(in)  :: MSG_ID      ! Message identifier
  real             , intent(in)  :: DATA(:,:)

! !INPUT/OUTPUT PARAMETERS: 
!
  type (T_CLASS)   , intent(inout) :: THIS
  integer          , intent(inout) :: ISTATUS

! !DESCRIPTION:  If the router for this interpolation has not yet
! been fully constructed, receive the other component's decomposition and
! initialize the router.
! .bp
! Pack the given data into the appropriate part of
! the SRC attribute vector associated with this.  If this is the
! last variable to be packed then Route the attribute vector.
!
! !DEFINED PARAMETERS:
! none

! !LOCAL VARIABLES: 
  type (String) :: NAME

  real, pointer :: RVECT(:)
  integer :: I
  integer :: J
  integer :: II

! !SEE ALSO:
!  none

! !REVISION HISTORY:
! 05/21/03 - Dan Schaffer - Initial version
!EOP ___________________________________________________________________

  if (.not. THIS%CONSTRUCTION_FINISHED) then

    if (MY_PE(THIS%MY_COMM) == 0) then

      if (SEQUENTIAL) then
        ASSERT(DECOMP_TAGS(THIS%HIS_COMP_ID, THIS%MY_COMP_ID, FROM) /= UNDEFINED_TAG)

#ifdef DEBUG
        print *, 'got to interp src 1', &
          DECOMP_TAGS(THIS%HIS_COMP_ID, THIS%MY_COMP_ID, FROM)
        call flush(6)
#endif

      else
        DECOMP_TAGS(THIS%HIS_COMP_ID, THIS%MY_COMP_ID, FROM) = BASE_DECOMP_TAG
      end if

      call MCT_RECV_DECOMP ( &
        THIS%DST_DECOMP_OF_SRC_GRID, &
        THIS%HIS_COMP_ID, &
	DECOMP_TAGS(THIS%HIS_COMP_ID, THIS%MY_COMP_ID, FROM), &
	ISTATUS)

      ASSERT(ISTATUS == 0)

    end if

    call MCT_BCAST_DECOMP ( &
      THIS%DST_DECOMP_OF_SRC_GRID, &
      0, &
      THIS%MY_COMM, &
      ISTATUS )

#ifdef DEBUG
    print *, 'in interp src, done broadcasting decomp'
#endif

    ASSERT(ISTATUS == 0)

    call MCT_ROUTER_INIT ( &
      THIS%SRC_DECOMP_OF_SRC_GRID, &
      THIS%DST_DECOMP_OF_SRC_GRID, &
      THIS%MY_COMM, &
      THIS%ROUTER )

#ifdef DEBUG
    print *, 'in interp src, router comp1id ', this%router%comp1id
    print *, 'in interp src, router comp2id ', this%router%comp2id
    print *, 'in interp src, router nprocs ', this%router%nprocs
    print *, 'in interp src, router pe_list ', this%router%pe_list
    print *, 'in interp src, router num_segs ', this%router%num_segs
#endif

    THIS%CONSTRUCTION_FINISHED = .true.

  end if

  THIS%CURRENT_VARIABLE = THIS%CURRENT_VARIABLE + 1

  call MCT_ATTR_NAME ( &
    NAME, &
    THIS%CURRENT_VARIABLE, &
    THIS%SRC_ATTR_VECT_SRC_GRID )

  if (WRF_IO_ERROR ( &
    trim(VAR_NAME) /= trim(MPEU_STRING_TO_CHAR(NAME)), &
    "Unexpected variable name found", ISTATUS)) return
    
  allocate(RVECT(size(DATA)), stat=ISTATUS)
  ASSERT(ISTATUS==0)
  II = 0
  do J = 1, size(DATA,2)
    do I = 1, size(DATA,1)
      II = II + 1
      RVECT(II) = DATA(I,J)
    end do
  end do

  call MCT_PUT_ATTR_VECT_DATA ( &
    THIS%SRC_ATTR_VECT_SRC_GRID, &
    VAR_NAME, &
    RVECT)

  deallocate(RVECT, stat=ISTATUS)
  ASSERT(ISTATUS == 0)

  if (THIS%CURRENT_VARIABLE == &
      MCT_ATTR_VECT_LENGTH(THIS%SRC_ATTR_VECT_SRC_GRID)) then

#ifdef DEBUG
    print *, THIS%MY_COMP_ID, 'got to interp src 11 tag ', THIS%TAG
    call flush(6)
#endif

    call MCT_ATTR_VECT_ISEND ( &
      THIS%SRC_ATTR_VECT_SRC_GRID, &
      THIS%ROUTER, &
      TAG=THIS%TAG)

! This appeared to be necessary to make progress in a timely manner
! in the concurrent case. Otherwise, the receiving model would sit 
! for awhile when run on MPICH on ijet
!
! However, I can't do this for the sequential case because
! if a processor sends to itself, it certainly can't have
! a block on that send which is what waitsend does!!!
    if (.not. SEQUENTIAL) then
      call MCT_ATTR_VECT_WAITSEND ( &
        THIS%ROUTER )
    end if

#ifdef DEBUG
    print *, THIS%MY_COMP_ID, 'got to interp src 12 tag ', THIS%TAG
    call flush(6)
#endif

    THIS%CURRENT_VARIABLE = 0

  end if

end subroutine INTERPOLATE_SRC_REAL_2D        

!BOP -------------------------------------------------------------------
!
! !IROUTINE: INTERPOLATE_DST_REAL_2D  - Interpolate to produce the dst DATA
!
! !INTERFACE:

subroutine INTERPOLATE_DST_REAL_2D     (THIS, VAR_NAME, MSG_ID, DATA, ISTATUS)
!
! !USES:
!
  use MCT_WRF_IO_UTILITIES_MODULE, only :  WRF_IO_ERROR

  use MCT_WRF_IO_UTILITIES_MODULE, only :  MY_PE

  use m_String, only : &
    String, &
    MPEU_STRING_TO_CHAR => toChar

  use m_AttrVect, only : &
    MCT_GET_ATTR_VECT_DATA => exportRAttr, &
    MCT_ATTR_VECT_LENGTH => nRAttr, &
    MCT_ATTR_NAME => getRList, &
    MCT_ATTR_VECT_SIZE => lsize
  
  use m_GlobalSegMapComms, only : MCT_RECV_DECOMP => RECV

  use m_GlobalSegMapComms, only : MCT_BCAST_DECOMP => BCAST

  use m_Router, only : MCT_ROUTER_INIT => INIT

  use m_Transfer, only : MCT_ATTR_VECT_RECV => RECV

  use m_MatAttrVectMul, only : &
    MCT_SPARSE_MATRIX_MULTIPLY => sMatAvMult

  implicit none

! !INPUT PARAMETERS: 
!
  character(LEN=*) , intent(in)  :: VAR_NAME    ! Variable name
  integer          , intent(in)  :: MSG_ID      ! Message identifier

! !OUTPUT PARAMETERS: 
  real             , intent(out)  :: DATA(:,:)

! !INPUT/OUTPUT PARAMETERS: 
!
  type (T_CLASS)   , intent(inout) :: THIS
  integer          , intent(inout) :: ISTATUS

! !DESCRIPTION:  
! If the router for this interpolation has not yet
! been fully constructed, receive the other component's decomposion and
! initialize the router.
! .bp
! If this is the first DST variable to be interpolated,
! use the Router to receive the data and then execute
! the sparse matrix multiply to do the interpolation .  For all variables, 
! unpack data corresponding to VAR_NAME from the resulting attribute vector 
! into DATA.
!
! !DEFINED PARAMETERS:
! none

! !LOCAL VARIABLES: 
  type (String) :: NAME

  real, pointer :: RVECT(:)
  integer :: LSIZE
  integer :: I
  integer :: J
  integer :: II


! !SEE ALSO:
!  none

! !REVISION HISTORY:
! 05/21/03 - Dan Schaffer - Initial version
!EOP ___________________________________________________________________

  if (.not. THIS%CONSTRUCTION_FINISHED) then

    if (MY_PE(THIS%MY_COMM) == 0) then

      if (SEQUENTIAL) then
        ASSERT(DECOMP_TAGS(THIS%HIS_COMP_ID, THIS%MY_COMP_ID, TO) /= UNDEFINED_TAG)

#ifdef DEBUG
        print *, 'got to interp dst 1', &
          DECOMP_TAGS(THIS%HIS_COMP_ID, THIS%MY_COMP_ID, TO)
#endif

      else
        DECOMP_TAGS(THIS%HIS_COMP_ID, THIS%MY_COMP_ID, TO) = BASE_DECOMP_TAG
      end if

      call MCT_RECV_DECOMP ( &
        THIS%SRC_DECOMP_OF_SRC_GRID, &
        THIS%HIS_COMP_ID, &
	DECOMP_TAGS(THIS%HIS_COMP_ID, THIS%MY_COMP_ID, TO), &
	ISTATUS)

      ASSERT(ISTATUS == 0)

    end if

    call MCT_BCAST_DECOMP ( &
      THIS%SRC_DECOMP_OF_SRC_GRID, &
      0, &
      THIS%MY_COMM, &
      ISTATUS )

#ifdef DEBUG
    print *, 'in interp dst, done broadcasting decomp'
#endif

    ASSERT(ISTATUS == 0)


    call MCT_ROUTER_INIT ( &
      THIS%DST_DECOMP_OF_SRC_GRID, &
      THIS%SRC_DECOMP_OF_SRC_GRID, &
      THIS%MY_COMM, &
      THIS%ROUTER )

#ifdef DEBUG
    print *, 'in interp dst, router comp1id ', this%router%comp1id
    print *, 'in interp dst, router comp2id ', this%router%comp2id
    print *, 'in interp dst, router nprocs ', this%router%nprocs
    print *, 'in interp dst, router pe_list ', this%router%pe_list
    print *, 'in interp dst, router num_segs ', this%router%num_segs
#endif

    THIS%CONSTRUCTION_FINISHED = .true.

  end if

  THIS%CURRENT_VARIABLE = THIS%CURRENT_VARIABLE + 1

  call MCT_ATTR_NAME ( &
    NAME, &
    THIS%CURRENT_VARIABLE, &
    THIS%DST_ATTR_VECT_DST_GRID )

  if (WRF_IO_ERROR ( &
    trim(VAR_NAME) /= trim(MPEU_STRING_TO_CHAR(NAME)), &
    "Unexpected variable name found", ISTATUS)) return
    
  if (THIS%CURRENT_VARIABLE == 1) then

#ifdef DEBUG
    print *, THIS%MY_COMP_ID, 'receiving with tag ', THIS%TAG
    call flush(6)
#endif

    call MCT_ATTR_VECT_RECV ( &
      THIS%DST_ATTR_VECT_SRC_GRID, &
      THIS%ROUTER, &
      TAG=THIS%TAG)

#ifdef DEBUG
    print *, THIS%MY_COMP_ID, 'done receiving with tag ', THIS%TAG
    call flush(6)
#endif

#ifdef DEBUG
    print *, 'got to interp dst 3'
    call flush(6)
#endif

    call MCT_SPARSE_MATRIX_MULTIPLY ( &
      THIS%DST_ATTR_VECT_SRC_GRID, &
      THIS%SPARSE_MATRIX, &
      THIS%DST_ATTR_VECT_DST_GRID )

#ifdef DEBUG
      print *, 'got to interp dst 4'
      call flush(6)
#endif

  end if

  allocate(RVECT(size(DATA)),stat=ISTATUS)
  ASSERT(ISTATUS==0)
  call MCT_GET_ATTR_VECT_DATA ( &
    THIS%DST_ATTR_VECT_DST_GRID, &
    VAR_NAME, &
    RVECT, &
    LSIZE)

  II = 0
  do J = 1, size(DATA,2)
    do I = 1, size(DATA,1)
      II = II + 1
      DATA(I,J) = RVECT(II)
    end do
  end do

  deallocate(RVECT, stat=ISTATUS)
  ASSERT(ISTATUS==0)
  if (THIS%CURRENT_VARIABLE == &
      MCT_ATTR_VECT_LENGTH(THIS%DST_ATTR_VECT_DST_GRID)) then

    THIS%CURRENT_VARIABLE = 0

  end if
end subroutine INTERPOLATE_DST_REAL_2D        

!BOP -------------------------------------------------------------------
!
! !IROUTINE: MATCHES         
!
! !INTERFACE:

logical function MATCHES (THIS, ID)
!
! !USES:
!
  implicit none

! !INPUT PARAMETERS: 
!
  type (T_CLASS)   , intent(in)  :: THIS
  character(LEN=*) , intent(in) :: ID

! !OUTPUT PARAMETERS: 
!
! None

! !DESCRIPTION:  Returns true if the identifier of THIS is ID
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

  MATCHES = trim(ID) == trim(THIS%ID)
end function MATCHES        

end module MCT_WRF_IO_INTERPOLATOR_MODULE
