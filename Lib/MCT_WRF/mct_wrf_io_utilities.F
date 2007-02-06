!BOP -------------------------------------------------------------------
!
! !MODULE: MCT_WRF_IO_UTILITIES_MODULE - MCT WRF I/O API package
!          
! !INTERFACE:

 module MCT_WRF_IO_UTILITIES_MODULE
!
! !DESCRIPTION:
!   Various utilities for MCT WRF I/O API

! !USES:

  implicit none

#include "macros.h"

  private	! except

! !PUBLIC TYPES:
!
! !PUBLIC MEMBER FUNCTIONS:
!
  public :: MY_PE          ! My process number in the given communicator
  public :: WRF_IO_ERROR   ! Check and process error messages
  public :: GET_KEY_VALUE  ! Retrieve a value from a key-value pair
                           ! in a string.
  public :: NUM_UNIQUE     ! Determine the number of unique values in an array
  public :: VALUE_IN_ARRAY ! Determine if a value is contained in an array
  public :: Assert         ! Check if the provided logical expression is true.
  public :: TO_LOWER       ! Convert to lower case
  public :: EXPAND         ! Expand an array of various object types
  public :: MAX_STR        ! Maximum size of any string.
  public :: ERROR_CHECK    ! Check an error condition

! !DEFINED PARAMETERS:

! !SEE ALSO:
!  none

! !REVISION HISTORY:
! 05/30/03 - Dan Schaffer - Initial version

!EOP ___________________________________________________________________


  interface GET_KEY_VALUE
    module procedure GET_KEY_VALUE_STRING
  end interface

  interface NUM_UNIQUE
    module procedure NUM_UNIQUE_I4
    module procedure NUM_UNIQUE_CHAR
  end interface

  interface VALUE_IN_ARRAY
    module procedure VALUE_IN_ARRAY_I4
    module procedure VALUE_IN_ARRAY_CHAR
  end interface

  interface Assert
    module procedure ASSERT_STRING
    module procedure ASSERT_LINE_FILE
  end interface

  interface EXPAND
    module procedure EXPAND_STRING
    module procedure EXPAND_R4
    module procedure EXPAND_R8
    module procedure EXPAND_I4
    module procedure EXPAND_I8
    module procedure EXPAND_C4
    module procedure EXPAND_C8
    module procedure EXPAND_L4
    module procedure EXPAND_L8
  end interface

! !DEFINED PARAMETERS:
  integer, parameter :: MAX_STR = 128

!NAMED CONSTANTS USED TO DEFINE VARIABLE TYPES
  INTEGER, PARAMETER :: L1K = 1        !1 byte integer
  INTEGER, PARAMETER :: L2K = 2        !2 byte integer
  INTEGER, PARAMETER :: L4K = 4        !4 byte integer

  INTEGER, PARAMETER :: I1K = SELECTED_INT_KIND(2) !1 byte integer
  INTEGER, PARAMETER :: I2K = SELECTED_INT_KIND(4) !2 byte integer
  INTEGER, PARAMETER :: I4K = SELECTED_INT_KIND(9) !4 byte integer
  INTEGER, PARAMETER :: R4K = KIND(1.0)            !single precision
  INTEGER, PARAMETER :: &
    R8K = SELECTED_REAL_KIND(2*PRECISION(1.0_R4K)) !double precision

contains

!BOP -------------------------------------------------------------------
!
! !IROUTINE: MY_PE         - My process number in COMM
!
! !INTERFACE:

 integer function MY_PE ( &
   COMM)
!
! !USES:
!
  implicit none

! !INPUT PARAMETERS: 
!
  integer          , intent(in)  :: COMM

! !OUTPUT PARAMETERS: 
!
! !DESCRIPTION:  Return rank of my process in COMM
!
! !DEFINED PARAMETERS:
! none

! !LOCAL VARIABLES: 
  integer                        :: ERROR

! !SEE ALSO:
! none

! !REVISION HISTORY:
! 05/30/03 - Dan Schaffer - Initial version
!EOP ___________________________________________________________________

   call MPI_COMM_RANK(COMM, MY_PE, ERROR)
   ASSERT(ERROR == 0)
 end function MY_PE        


!BOP -------------------------------------------------------------------
!
! !IROUTINE: WRF_IO_ERROR         
!
! !INTERFACE:

logical function WRF_IO_ERROR ( &
  CONDITION, &
  MSG, &
  ISTATUS)
!
! !USES:
!
  implicit none

! !INPUT PARAMETERS: 
!
  logical          , intent(in)  :: CONDITION
  character(len=*) , intent(in)  :: MSG

! !OUTPUT PARAMETERS: 
!
  integer          , intent(out) :: ISTATUS

! !DESCRIPTION:  
! Check to see if CONDITION is true.  If so, print MSG,
! store WRF_MCT_IO_ERROR in ISTATUS and return .true.
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

  include 'wrf_io_flags.h'

  if (CONDITION) then
    print *, MSG
    ISTATUS = WRF_MCT_IO_ERROR
    WRF_IO_ERROR = .true.
  else
    ISTATUS = 0
    WRF_IO_ERROR = .false.
  end if
end function WRF_IO_ERROR        

!BOP -------------------------------------------------------------------
!
! !IROUTINE: GET_KEY_VALUE_STRING
!
! !INTERFACE:
subroutine GET_KEY_VALUE_STRING ( &
  STRING, &
  KEY, &
  VALUE, &
  ISTATUS, &
  IGNORE_CASE )

#define STRING_TYPE
#define TYPE character(len=MAX_STR)

#include "get_key_value.code"
#undef STRING_TYPE

#undef TYPE

end subroutine GET_KEY_VALUE_STRING

!BOP -------------------------------------------------------------------
!
! !IROUTINE: NUM_UNIQUE_I4         
!
! !INTERFACE:

function NUM_UNIQUE_I4 ( &
  ARR) result(RES)
#define IN_TYPE integer(kind=I4K)
#define INIT_VALUE -9999
#include "num_unique.code"
end function NUM_UNIQUE_I4        

!BOP -------------------------------------------------------------------
!
! !IROUTINE: NUM_UNIQUE_CHAR         
!
! !INTERFACE:

function NUM_UNIQUE_CHAR ( &
  ARR) result(RES)
#define IN_TYPE character(len=MAX_STR)
#define INIT_VALUE ""
#include "num_unique.code"
end function NUM_UNIQUE_CHAR        

!BOP -------------------------------------------------------------------
!
! !IROUTINE: VALUE_IN_ARRAY_I4         
!
! !INTERFACE:

function VALUE_IN_ARRAY_I4 ( &
  VALUE, &
  ARRAY) result(RES)

#define IN_TYPE integer(kind=I4K)
#include "value_in_arr.code"

end function VALUE_IN_ARRAY_I4

!BOP -------------------------------------------------------------------
!
! !IROUTINE: VALUE_IN_ARRAY_CHAR         
!
! !INTERFACE:

function VALUE_IN_ARRAY_CHAR ( &
  VALUE, &
  ARRAY) result(RES)

#define IN_TYPE character(len=MAX_STR)
#include "value_in_arr.code"

end function VALUE_IN_ARRAY_CHAR

!BOP -------------------------------------------------------------------
!
! !IROUTINE: ASSERT_STRING         
!
! !INTERFACE:

subroutine ASSERT_STRING ( &
  LOGICAL_EXPRESSION, &
  ERROR_STRING)
!
! !USES:
!
  implicit none

! !INPUT PARAMETERS: 
!
  logical         , intent(in) :: LOGICAL_EXPRESSION
  character(LEN=*), intent(in) :: ERROR_STRING

! !OUTPUT PARAMETERS: 
!
! none

! !DESCRIPTION:  If the given logical expression is false then
! print the given error string and stop.
!
! !DEFINED PARAMETERS:
! none

! !LOCAL VARIABLES: 
! none

! !SEE ALSO:
!  none

! !REVISION HISTORY:
! 12/29/03 - Dan Schaffer - Initial version
!EOP ___________________________________________________________________
  if (.not. LOGICAL_EXPRESSION ) then
    print *, 'ASSERTION ERROR: ', ERROR_STRING
    stop
  end if

end subroutine ASSERT_STRING        

!BOP -------------------------------------------------------------------
!
! !IROUTINE: ASSERT_LINE_FILE         
!
! !INTERFACE:

subroutine ASSERT_LINE_FILE ( &
  LOGICAL_EXPRESSION, &
  line, &
  file)
!
! !USES:
!
  implicit none

! !INPUT PARAMETERS: 
!
  logical         , intent(in) :: LOGICAL_EXPRESSION
  integer         , intent(in) :: line
  character(LEN=*), intent(in) :: file

! !OUTPUT PARAMETERS: 
!
! none

! !DESCRIPTION:  If the given logical expression is false then
! print out an error message including the given line number and
! file name and stop
!
! !DEFINED PARAMETERS:
! none

! !LOCAL VARIABLES: 
! none

! !SEE ALSO:
!  none

! !REVISION HISTORY:
! 12/29/03 - Dan Schaffer - Initial version
!EOP ___________________________________________________________________

  if ( .not. LOGICAL_EXPRESSION ) then
   print *, "ASSERTION Error in line ", line, " of file ", file
   stop
  end if

end subroutine ASSERT_LINE_FILE        

!BOP -------------------------------------------------------------------
!
! !IROUTINE: TO_LOWER         
!
! !INTERFACE:

function TO_LOWER ( &
  STRING)
!
! !USES:
!
  implicit none

! !INPUT PARAMETERS: 
!
character(len=*), intent(in) :: STRING


! !OUTPUT PARAMETERS: 
!
character(len=len(STRING)) :: TO_LOWER

! !DESCRIPTION:  Convert letters in the given string to lower case.
!
! !DEFINED PARAMETERS:
! none

! !LOCAL VARIABLES: 
  integer :: I

! !SEE ALSO:
!  none

! !REVISION HISTORY:
! 08/04/03 - Dan Schaffer - Initial version
!EOP ___________________________________________________________________

  do I = 1, len(STRING)
    if (STRING(I:I) >= 'A' .and. STRING(I:I) <= 'Z') then
      TO_LOWER(I:I) = ACHAR(IACHAR(STRING(I:I)) + IACHAR('a') - IACHAR('A'))
    else
      TO_LOWER(I:I) = STRING(I:I)
    end if
  end do

end function TO_LOWER        

!BOP -------------------------------------------------------------------
!
! !IROUTINE: EXPAND_STRING
!
! !INTERFACE:
! 
subroutine EXPAND_STRING ( IN, ITEMS, AFTER )

#define IN_TYPE character(len=MAX_STR)

#include "arrays.xpnd.code"

#undef IN_TYPE

end subroutine EXPAND_STRING

!BOP -------------------------------------------------------------------
!
! !IROUTINE: EXPAND_R4
!
! !INTERFACE:
! 
subroutine EXPAND_R4 ( IN, ITEMS, AFTER )

#define IN_TYPE real(kind=4)

#include "arrays.xpnd.code"

#undef IN_TYPE

end subroutine EXPAND_R4

!BOP -------------------------------------------------------------------
!
! !IROUTINE: EXPAND_R8
!
! !INTERFACE:
! 
subroutine EXPAND_R8 ( IN, ITEMS, AFTER )

#define IN_TYPE real(kind=8)

#include "arrays.xpnd.code"

#undef IN_TYPE

end subroutine EXPAND_R8

!BOP -------------------------------------------------------------------
!
! !IROUTINE: EXPAND_I4
!
! !INTERFACE:
! 
subroutine EXPAND_I4 ( IN, ITEMS, AFTER )

#define IN_TYPE integer(kind=4)

#include "arrays.xpnd.code"

#undef IN_TYPE

end subroutine EXPAND_I4

!BOP -------------------------------------------------------------------
!
! !IROUTINE: EXPAND_I8
!
! !INTERFACE:
! 
subroutine EXPAND_I8 ( IN, ITEMS, AFTER )

#define IN_TYPE integer(kind=8)

#include "arrays.xpnd.code"

#undef IN_TYPE

end subroutine EXPAND_I8

!BOP -------------------------------------------------------------------
!
! !IROUTINE: EXPAND_C4
!
! !INTERFACE:
! 
subroutine EXPAND_C4 ( IN, ITEMS, AFTER )

#define IN_TYPE complex(kind=4)

#include "arrays.xpnd.code"

#undef IN_TYPE

end subroutine EXPAND_C4

!BOP -------------------------------------------------------------------
!
! !IROUTINE: EXPAND_C8
!
! !INTERFACE:
! 
subroutine EXPAND_C8 ( IN, ITEMS, AFTER )

#define IN_TYPE complex(kind=8)

#include "arrays.xpnd.code"

#undef IN_TYPE

end subroutine EXPAND_C8

!BOP -------------------------------------------------------------------
!
! !IROUTINE: EXPAND_L4
!
! !INTERFACE:
! 
subroutine EXPAND_L4 ( IN, ITEMS, AFTER )

#define IN_TYPE logical(kind=4)

#include "arrays.xpnd.code"

#undef IN_TYPE

end subroutine EXPAND_L4

!BOP -------------------------------------------------------------------
!
! !IROUTINE: EXPAND_L8
!
! !INTERFACE:
! 
subroutine EXPAND_L8 ( IN, ITEMS, AFTER )

#define IN_TYPE logical(kind=8)

#include "arrays.xpnd.code"

#undef IN_TYPE

end subroutine EXPAND_L8

!BOP -------------------------------------------------------------------
!
! !IROUTINE: ABORT         - Abort a program
!
! !INTERFACE:

subroutine ABORT ( &
  MSG)
!
! !USES:
!
  implicit none

! !INPUT PARAMETERS: 
!
  character(len=*) , intent(in)  :: MSG

! !OUTPUT PARAMETERS: 
!
! none

! !DESCRIPTION:  Print MSG and stop
!
! !DEFINED PARAMETERS:
! none

! !SEE ALSO:
! none

! !REVISION HISTORY:
! 06/27/03 - Dan Schaffer - Initial version
!EOP ___________________________________________________________________

  print *, MSG
  stop
end subroutine ABORT        

!BOP -------------------------------------------------------------------
!
! !IROUTINE: ERROR_CHECK         - Check an error condition
!
! !INTERFACE:

subroutine ERROR_CHECK ( &
  COND, &
  MSG)
!
! !USES:
!
  implicit none

! !INPUT PARAMETERS: 
!
  logical          , intent(in)  :: COND
  character(len=*) , intent(in)  :: MSG

! !OUTPUT PARAMETERS: 
!
! none

! !DESCRIPTION:  If COND is false then ABORT with MSG
!
! !DEFINED PARAMETERS:
! none

! !LOCAL VARIABLES: 
! none

! !SEE ALSO:
!  none

! !REVISION HISTORY:
! 06/27/03 - Dan Schaffer - Initial version
!EOP ___________________________________________________________________

  if (.not. COND) then
    call ABORT(MSG)
  end if
end subroutine ERROR_CHECK        

      SUBROUTINE CHR_DELETE(ch,str)
!     Delete all occurances of character CH in STR.
!
!     INPUT:
!        ch: a character                  [A1]
!       str: character variable           [A*]
!
!     OUTPUT:
!       str: modified character variable  [A*]
!
!      LIBRARY: [STRPAK]
!       STATUS: [STABLE]
!     CATEGORY: [LEVEL 1]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 02-JAN-1995
!     REVISIONS HISTORY:
!
!       08-NOV-1997 RD STEWART
!          Revised to meet standardized FLIB coding conventions
!       21-DEC-1998 RD STEWART
!          Changed F77 style comments to F90.  Routine renamed
!          from DELCH to CHR_ERASE.
!
!     COMMENTS: The length of STR on return, minus trailing blanks,
!               will always be less than or equal to the original
!               length minus trailing blanks.
!
      CHARACTER(LEN=1) ch,buf
      CHARACTER(LEN=*) str
      INTEGER ILEN,I,J
      INTENT(in) :: ch
      INTENT(inout) :: str

      ILEN=LEN(str)
      j=0
      DO i=1,ILEN
        IF (str(i:i).NE.ch) THEN
!         Keep this character
          j=j+1
          buf=str(i:i)
          str(j:j) = buf
        ENDIF
      ENDDO
      DO i=j+1,ILEN
        str(i:i)=' '
      ENDDO

      RETURN
      END subroutine CHR_DELETE

      SUBROUTINE STR_FAM(str1,str2, j1,j2,iflg)
!     Attempt to locate a case insensitive match for STR1 in STR2.  If
!     the search is successfully (iflg=0), the location of STR1 in
!     STR2 is returned as J1,J2, i.e., STR1 occurs at STR2(J1:J2).
!     Trailing blanks in STR1 are ignored.
!
!     INPUT:
!       STR1: a chracter variable                              [A*]
!       STR2: a chracter variable                              [A*]
!
!     OUTPUT:
!         J1: location of the 1st letter of STR1 in STR2       [I4]
!         J2: location of the last letter of STR1 in STR2      [I4]
!         IFLG: status flag                                    [I4]
!               -1: sub-string STR1 not found in STR2
!                0: found STR1 in STR2
!
!      LIBRARY: [STRPAK]
!       STATUS: [BETA]
!     CATEGORY: [LEVEL 2]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 08-NOV-1997
!     REVISIONS HISTORY:
!       30-DEC-1998 RD STEWART
!          Changed F77 style comments to F90.
!
!     COMMENTS:
!
      CHARACTER(LEN=*) STR1,STR2
      INTEGER J1,J2,IFLG
      INTEGER i,j,ilen,jlen
      LOGICAL mat,again,more,found
      INTENT(in) :: str1,str2
      INTENT(out) :: j1,j2,iflg

      !Initialize output parameters.
      iflg=0
      j1=1
      j2=1

!     Check for bad input.
      ilen=len_trim(str1)  !length minus trailing blanks
      jlen=LEN(str2)
      IF (jlen.EQ.0) THEN
        !STR2 is a null string
        iflg=-1
        RETURN
       ELSEIF (ilen.EQ.0) THEN
         ilen=LEN(str1)
         IF (ilen.EQ.0) THEN
           !STR1 is a null string
           iflg=-1
           RETURN
         ENDIF
       ENDIF

!     Initial condition
      j1=1
      j2=1
      iflg=-1   !no sub-string found

!     Begin search for sub-string
      j=0
      again=(jlen.GE.ilen)   !STR2 must be longer than STR1 to find a match.
      DO WHILE (again)
        j=j+1
        found=.false.
        mat=CHR_IMAT(str1(1:1),str2(j:j))
        IF (mat) THEN
!         ith letter in STR1 matches the jth letter
!         in STR2. Set position of 1st letter and
!         continue until all of STR1 is found, a character
!         in STR1 doesn't match one in STR2, or the
!         end of STR1 or STR2 is encountered.
          j1=j
          i=1
          found=(ilen.EQ.1)
          more=((.not.found).AND.((j.LT.jlen).AND.(i.LT.ilen)))
          DO WHILE (more)
            i=i+1
            j=j+1
            mat=CHR_IMAT(str1(i:i),str2(j:j))
            found=(mat.AND.(i.EQ.ilen))
            more=((mat).AND.((i.LT.ilen).AND.(j.LT.jlen)))
          ENDDO
          IF (found) THEN
            j2=j
            iflg=0
          ELSE
!           Reset index to STR1
            j=j1
          ENDIF
        ENDIF
        again=((.not.found).AND.(j.LT.jlen))
      ENDDO
      IF (iflg.NE.0) THEN
        j1=1
        j2=1
      ENDIF

      RETURN
      END subroutine STR_FAM

      SUBROUTINE STR_FXM(str1,str2, j1,j2,iflg)
!     Attempt to locate an exact match for STR1 in STR2.  If the
!     search is successfully (iflg=0), the location of STR1 in
!     STR2 is returned as J1,J2, i.e., STR1 occurs at STR2(J1:J2).
!     Trailing blanks in STR1 are ignored.
!
!     INPUT:
!       STR1: a chracter variable                              [A*]
!       STR2: a chracter variable                              [A*]
!
!     OUTPUT:
!         J1: location of the 1st letter of STR1 in STR2       [I4]
!         J2: location of the last letter of STR1 in STR2      [I4]
!         IFLG: status flag                                    [I4]
!               -1: sub-string STR1 not found in STR2
!                0: found STR1 in STR2
!
!      LIBRARY: [STRPAK]
!       STATUS: [BETA]
!     CATEGORY: [LEVEL 2]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 08-NOV-1997
!     REVISIONS HISTORY:
!       30-DEC-1998 RD STEWART
!          Changed F77 style comments to F90.
!
!     COMMENTS:
!
      CHARACTER(LEN=*) STR1,STR2
      INTEGER J1,J2,IFLG
      INTEGER i,j,ilen,jlen
      LOGICAL mat,again,more,found
      INTENT(in) :: str1,str2
      INTENT(out) :: j1,j2,iflg

      !Initialize output parameters.
      iflg=0
      j1=1
      j2=1

!     Check for bad input.
      ilen=len_trim(str1)  !length minus trailing blanks
      jlen=LEN(str2)
      IF (jlen.EQ.0) THEN
        !STR2 is a null string
        iflg=-1
        RETURN
       ELSEIF (ilen.EQ.0) THEN
         ilen=LEN(str1)
         IF (ilen.EQ.0) THEN
           !STR1 is a null string
           iflg=-1
           RETURN
         ENDIF
       ENDIF

!     Initial condition
      iflg=-1   !no sub-string found

!     Begin search for sub-string
      j=0
      again=(jlen.GE.ilen)   !STR2 must be longer than STR1 to find a match.
      DO WHILE (again)
        j=j+1
        found=.false.
        mat=(str1(1:1).EQ.str2(j:j))
        IF (mat) THEN
!         ith letter in STR1 matches the jth letter
!         in STR2. Set position of 1st letter and
!         continue until all of STR1 is found, a character
!         in STR1 doesn't match one in STR2, or the
!         end of STR1 or STR2 is encountered.
          j1=j
          i=1
          found=(ilen.EQ.1)
          more=((.not.found).AND.((j.LT.jlen).AND.(i.LT.ilen)))
          DO WHILE (more)
            i=i+1
            j=j+1
            mat=(str1(i:i).EQ.str2(j:j))
            found=(mat.AND.(i.EQ.ilen))
            more=((mat).AND.((i.LT.ilen).AND.(j.LT.jlen)))
          ENDDO
          IF (found) THEN
            j2=j
            iflg=0
          ELSE
!           Reset index to STR1
            j=j1
          ENDIF
        ENDIF
        again=((.not.found).AND.(j.LT.jlen))
      ENDDO
      IF (iflg.NE.0) THEN
        j1=1
        j2=1
      ENDIF

      RETURN
      END SUBROUTINE STR_FXM

FUNCTION CHR_IMAT(ch1,ch2)
!     CHR_IMAT is returned as TRUE if CH1 is an upper-
!     or lowercase match for CH2.
!
!     INPUT:
!        ch1: a character                      [A1]
!        ch2: a character                      [A1]
!
!     OUTPUT:
!       CHR_IMAT: logical variable indicating  [L4]
!                 whether CH1 is a case
!                 insensitive match to CH2
!
!      LIBRARY: [STRPAK]
!       STATUS: [STABLE]
!     CATEGORY: [LEVEL 1]
!
!     AUTHOR: Robert D. Stewart
!             5127 West 3rd Ave.
!             Kennewick, WA 99336
!             trebor@mad.scientist.com
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 07-JUL-1997
!     REVISIONS HISTORY:
!
!       08-NOV-1997 RD STEWART
!          Revised to meet standardized FLIB coding conventions
!       31-DEC-1998 RD STEWART
!          Changed comments to Fortran 90 style and renamed
!          routine from ISMAT to CHR_IMAT.  Removed (trival)
!          option for case sensitive match between CH1 and CH2.
!          Made algorithm more computationally efficient.
!
!     COMMENTS:
!
      CHARACTER(LEN=1) ch1,ch2
      LOGICAL Ok
      INTEGER i1,i2
      INTENT(in) :: ch1,ch2
      LOGICAL CHR_IMAT

      ok=(ch1.EQ.ch2)
      IF (.not.ok) THEN
!       Check for case insensitive match. NOTE: the ASCII
!       code for an uppercase letter can be found by adding
!       32 to the ASCII code for the corresponding
!       lowercase letter.
        i1 = ICHAR(ch1)
        i2 = ICHAR(ch2)
        IF ((i1.GE.65).AND.(i1.LE.90)) THEN
          !ch1 is uppercase
          ok=((i1+32).EQ.i2)
        ELSEIF ((i1.GE.97).AND.(i1.LE.122)) THEN
          !ch1 is lowercase
          ok=((i1-32).EQ.i2)
        ENDIF
      ENDIF
      CHR_IMAT = ok

      RETURN
      END FUNCTION CHR_IMAT

 end module MCT_WRF_IO_UTILITIES_MODULE
