module Parser_mod
!@sum procedures to read parameters from the rundeck into the database
!@auth I. Aleinov and T. Clune
!@ver 1.1     
  implicit none

  public :: Parser_type ! data type
  public :: parse
  public :: stripComment
  public :: skipHeader
  public :: isEndData
  public :: setCommentCharacters
  public :: setBeginData
  public :: setEndData
  public :: splitTokens
  public :: writeFormatted

  public :: MAX_LEN_LINE
  public :: MAX_LEN_TOKEN
  integer, parameter :: MAX_COMMENT_CHARACTERS = 3
  integer, parameter :: MAX_TOKEN_SEPARATORS   = 2
  integer, parameter :: MAX_LEN_TOKEN = 32
  integer, parameter :: MAX_LEN_LINE = 256
  character(len=*), parameter :: ENTIRE_LINE = '(a256)' ! MAX_LEN_LINE
  character(len=*), parameter :: ENTIRE_TOKEN = '(a33)' ! MAX_LEN_TOKEN

  type Parser_type
    character(len=MAX_COMMENT_CHARACTERS) :: commentCharacters = '!#' ! legacy default
    character(len=MAX_TOKEN_SEPARATORS) :: tokenSeparators = ' =,'     ! legacy default
    character(len=MAX_LEN_LINE) :: beginData = '&&PARAMETERS'         ! legacy default
    character(len=MAX_LEN_LINE) :: endData = '&&END_PARAMETERS'     ! legacy default
  end type Parser_type

  type (Parser_type), save :: globalParser

  interface writeFormatted
    module procedure writeFormatted_parser
  end interface

contains

  function strip_comment( str ) result(newStr)
    ! remove comment at the end of the line. comments symbols: !#
    character*(*), intent(in) :: str
    character(len=len(str)) :: newStr
    
    newStr = stripComment(globalParser, str)

  end function strip_comment

  subroutine skip_junk( str )
    character*(*) str

    do while ( len_trim( str ) > 0 .and. scan( str, ' =,' ) == 1 )
      str = str(2:)
    enddo
    return
  end subroutine skip_junk

  subroutine sread_int( str, value )
    character*(*) str
    integer value
    integer n

    read ( str, * ) value

    ! remove chars till the next [ =,]
    n = scan( str, ' =,' )
    str = str(n+1:)

    call skip_junk( str )

    return
  end subroutine sread_int

  subroutine sread_real( str, value )
    character*(*) str
    real*8 value
    integer n

    read ( str, * ) value

    ! remove chars till the next [ =,]
    n = scan( str, ' =,' )
    str = str(n+1:)

    call skip_junk( str )

  end subroutine sread_real

  subroutine sread_char( str, value )
    character*(*) str
    character*(*) value
    !character*256 tstr
    integer n, n1

    ! replace '=' with space if not quoted
    n1 = scan( str, '''' )
    n  = scan( str, '=' )
    if ( n>0 .and. ( n1==0 .or. n<n1 ) ) str(n:n) = ' '

    read ( str, * ) value

    if ( scan( str, '''' ) == 1 ) then  ! quoted string
      str = str(2:)
      n = scan( str, '''' )
      str = str(n+1:)
    else  ! remove chars till the next [ =,]
      n = scan( str, ' =,' ) 
      str = str(n+1:)
    endif

    call skip_junk( str )

  end subroutine sread_char

  subroutine parse_params( kunit )
  !@sub Parse the information in the file <kunit>
  !@+ and load the parameters into the database. It overwrites
  !@+ existing parameters.
  !@+
  !@var integer, intent(in) :: kunit
  !@+     Unit number of open rundeck file (it should contain a block
  !@+     starting with &&PARAMETERS and ending with &&END_PARAMETERS)

    use Dictionary_mod
    integer, parameter :: MAXDIM=128
    integer, intent(in) :: kunit
    character*256 bufs
    character*32 name
    character*1 type
    integer np
    integer ivars(MAXDIM)
    real*8 rvars(MAXDIM)
    character(len=MAX_CHAR_LEN+1) :: cvars(MAXDIM)

    ! skip unrelated stuff
    do
      read( kunit, '(a256)', err=666, end=667 ) bufs
      if ( len_trim(bufs) < 1 ) cycle
      read( bufs, * ) name
      if ( name == '&&PARAMETERS' ) exit
    enddo

    do
      read( kunit, '(a256)', err=666, end=666 ) bufs

      if ( len_trim(bufs) < 1 ) cycle
      if ( len_trim(bufs) > 255 ) &
           call stop_model("parse_params: rundeck line too long",255)

      bufs = strip_comment( bufs )
      call skip_junk( bufs )

      if ( len_trim(bufs) < 1 ) cycle

      !read the name of the variable
      call sread_char( bufs, name )

      if ( name == '&&END_PARAMETERS' ) exit  ! end of list 

      if ( len_trim(bufs) < 1 ) then
        print *,'PARSER: no values were given to param: ', name
        call stop_model('PARSER error',255)
      endif

      ! now check the type of variables
      if ( scan( bufs, '''' ) > 0 ) then
        type = 'c'
      else if ( scan( bufs, '.eEdD' ) > 0 ) then
        type = 'r'
      else
        type = 'i'
      endif

      select case ( type )
      case ('i')
        np = 0
        do while ( len_trim(bufs) > 0 )
          np = np+1
          if (np > MAXDIM) call stop_model("parse_params: increase MAXDIM",255)
          call sread_int( bufs, ivars(np) )
        end do
        call set_param( name, ivars, np, 'or' )
      case ('r')
        np = 0
        do while ( len_trim(bufs) > 0 )
          np = np+1
          if (np > MAXDIM) call stop_model("parse_params: increase MAXDIM",255)
          call sread_real( bufs, rvars(np) )
        end do
        call set_param( name, rvars, np, 'or' )
      case ('c')
        np = 0
        do while ( len_trim(bufs) > 0 )
          np = np+1
          if (np > MAXDIM) call stop_model("parse_params: increase MAXDIM",255)
          call sread_char( bufs, cvars(np) )
        end do
        call set_param( name, cvars, np, 'or' )
      end select

    enddo

    return
666 print *, 'PARSER: Error reading params'
    call stop_model( 'PARSER: Error reading params', 255 )
667 print *, 'PARSER: No &&PARAMETERS or &&END_PARAMETERS found'
    call stop_model( &
    &     'PARSER: No &&PARAMETERS or &&END_PARAMETERS found',255)
  end subroutine parse_params

  subroutine parseLine(this, line, name, value)
!@sum Attempts to create a key value pair by parsing a line of text
!@+ into separate tokens (via splitTokens).  Throws an exception if
!@+ less than two tokens are found.   (Should probably be modified
!@+ to also throw exception if first token is not an exceptable
!@+ key.)
    use Attributes_mod
    use AttributeDictionary_mod
    use AbstractAttribute_mod
    use StringUtilities_mod
    type (Parser_type), intent(in) :: this
    character(len=*), intent(in) :: line
    character(len=*), intent(out) :: name
    class (AbstractAttribute), pointer :: value

    character(len=MAX_LEN_TOKEN), pointer :: tokens(:)

    tokens => splitTokens(this, line)
    if (size(tokens) < 2) then
      call stop_model('Parser_mod: syntax error in input unit.', 14)
      return
    end if

    name = tokens(1)
    if (isInteger(tokens(2))) then
      if (size(tokens) == 2) then
!!$        value = newAttribute(stringToInteger(tokens(2))) ! scalar
         allocate(value,source=newAttribute(stringToInteger(tokens(2)))) ! scalar
      else
!!$        value = newAttribute(stringToInteger(tokens(2:)))
        allocate(value,source=newAttribute(stringToInteger(tokens(2:))))
      end if
    else if (isReal64(tokens(2))) then
      if (size(tokens) == 2) then
!!$        value = newAttribute(stringToRealDP(tokens(2))) ! scalar
        allocate(value,source=newAttribute(stringToRealDP(tokens(2))))
      else
!!$        value = newAttribute(stringToRealDP(tokens(2:)))
        allocate(value,source=newAttribute(stringToRealDP(tokens(2:))))
      end if
    else if (isLogical(tokens(2))) then
      if (size(tokens) == 2) then
!!$        value = newAttribute(stringToLogical(tokens(2))) ! scalar
        allocate(value, source=newAttribute(stringToLogical(tokens(2)))) ! scalar
      else
!!$        value = newAttribute(stringToLogical(tokens(2:)))
        allocate(value, source=newAttribute(stringToLogical(tokens(2:))))
      end if
    else ! is string
      if (size(tokens) == 2) then
!!$        value = newAttribute(tokens(2)) ! scalar
        allocate(value, source=newAttribute(tokens(2))) ! scalar
      else
!!$        value = newAttribute(tokens(2:))
        allocate(value, source=newAttribute(tokens(2:)))
      end if
    end if

    deallocate(tokens)

  end subroutine parseLine

  function parse(this, unit, status) result(dictionary)
!@sum Parses text from input unit to populate a Dictionary object.
!@+ Skips header section at top.
    use AttributeDictionary_mod
    use AbstractAttribute_mod
    use AttributeHashMap_mod, only: MAX_LEN => MAX_LEN_KEY

    type (Parser_type), intent(in) :: this
    integer, intent(in) :: unit
    integer, intent(out) :: status
    type (AttributeDictionary) :: dictionary

    character(len=MAX_LEN_LINE) :: line
    character(len=MAX_LEN) :: attributeName
    class (AbstractAttribute), pointer :: attributeValue

    status = 0 ! unless ...
    dictionary = newAttributeDictionary() ! TODO - fix memory leak here

    call skipHeader(this, unit, status)
    if (status /= 0) return
    
    do
      read(unit,fmt=ENTIRE_LINE,iostat=status) line
      if (status /= 0) exit
      if (isEndData(this, line)) exit
      
      line = stripComment(this, line)
      if (len_trim(line) == 0) cycle ! skip

      call parseLine(this, line, attributeName, attributeValue)
      call dictionary%insert(trim(attributeName), attributeValue)
      deallocate(attributeValue)

    end do

  end function parse
  
  subroutine setCommentCharacters(this, commentCharacters)
!@sum Set the characters which should be interpreted as comment
!@+ when processing input.  Defaults are conventional Fortran "!#",
!@+ but this routine can be used to override.
    type (Parser_type), intent(inout) :: this
    character(len=*), intent(in) :: commentCharacters

    this%commentCharacters = commentCharacters

  end subroutine setCommentCharacters

  subroutine setTokenSeparators(this, tokenSeparators)
!@sum Override default characters used to separate tokens.
    type (Parser_type), intent(inout) :: this
    character(len=*), intent(in) :: tokenSeparators
    this%tokenSeparators = tokenSeparators
  end subroutine setTokenSeparators

  function stripComment(this, str) result(newStr)
!@sum Eliminate trailing comments in line of text. 
!TO DO - ignore comments in strings
    type (Parser_type), intent(in) :: this
    character*(*), intent(in) :: str
    character(len=len(str)) :: newStr

    integer :: n

    n = scan(str, trim(this%commentCharacters))
    select case (n)
    case (0)
      newStr = trim(str)
    case (1:)
      newStr = str(:n-1)
    end select

  end function stripComment

  subroutine setBeginData(this, beginData)
!@sum Override default string signfying beginning of data,
!@+ i.e. end of header.
    type (Parser_type), intent(inout) :: this
    character(len=*), intent(in) :: beginData
    this%beginData = beginData
  end subroutine setBeginData

  logical function isEndData(this, string)
!@sum Return true if string matches the semaphore for the end of data.
    type (Parser_type), intent(in) :: this
    character(len=*), intent(in) :: string
    
    isEndData = (trim(this%endData) == trim(string))

  end function isEndData

  subroutine setEndData(this, endData)
!@sum Override default string signfying end of data.
    type (Parser_type), intent(inout) :: this
    character(len=*), intent(in) :: endData
    this%endData = endData
  end subroutine setEndData

  subroutine skipHeader(this, unit, status)
!@sum Read from unit until semaphore for beginning of data
!@+ is located.  Returns nonzero status if semaphore is not
!@+ found before EOF.
    type (Parser_type), intent(in) :: this
    integer, intent(in) :: unit
    integer, intent(out) :: status

    character(len=MAX_LEN_LINE) :: line
    
    do
      read(unit,fmt=ENTIRE_LINE,iostat=status) line
      if (status < 0) then
        return
      end if

      if (trim(line) == this%beginData) exit
    end do

  end subroutine skipHeader

  logical function isInteger(string)
!@sum Returns true if string can be converted to an integer.
!@+ Should possibly relocate to GenericType_mod
    character(len=*), intent(in) :: string
    integer :: integerValue
    integer :: status
    
    read(string,'(i20)',iostat=status) integerValue
    isInteger = (status == 0)

  end function isInteger

  logical function isReal64(string)
!@sum Returns true if string can be converted to a double.
!@+ Should possibly relocate to GenericType_mod
    character(len=*), intent(in) :: string
    real*8 :: real64Value
    integer :: status
    
    read(string,'(g37.30)',iostat=status) real64Value
    isReal64 = (status == 0)

  end function isReal64

  logical function isLogical(string)
!@sum Returns true if string can be converted to a logical.
!@+ Allows for Fortran defaults AND other strings that are 
!@+ convenient equivalents.  (See implementation below.)
!@+ Should possibly relocate to GenericType_mod
    use StringUtilities_mod, only: toLowerCase
    character(len=*), intent(in) :: string
    
    select case (trim(toLowerCase(string)))
    case ('t','f','true','false','.true.','.false.')
      isLogical = .true.
    case default
      isLogical = .false.
    end select

  end function isLogical

  function splitTokens(this, string) result(tokens)
!@ Attempt to split a string into a set of independent tokens
!@+ based upon separaters specified in the Parser object.
!@+ A crude attempt is made to ensure that the 1st separator is '='.
!@+ A more comprehensive treatement of this should be made.
    type (Parser_type), intent(in) :: this
    character(len=*), intent(in) :: string
    character(len=MAX_LEN_TOKEN), pointer :: tokens(:)

    character(len=len(string)) :: buffer
    integer :: numTokens
    integer :: i, idxStart
    integer :: idxNextSeparator

    numTokens = countTokens(this, string)
    allocate(tokens(numTokens))
    
    buffer = adjustl(string)
    i = 0
    do  ! pull off tokens until
      if (len_trim(buffer) == 0) exit ! done
      i = i + 1
      idxStart = skipEmbeddedSeparators(buffer)
      idxNextSeparator = scan(buffer(idxStart:), trim(this%tokenSeparators))

      if (idxNextSeparator > 0) then
        tokens(i) = trim(buffer(:idxStart+idxNextSeparator-2))
        buffer = adjustl(buffer(idxStart+idxNextSeparator:))
      else
        tokens(i) = trim(buffer) ! take the rest
        exit
      end if
    end do

    if (numTokens >= 1) then
      if (scan(trim(tokens(1)), ' ') /= 0) then
        call stop_model('Parser_mod: Illegal syntax.  "=" not first separator.', 14)
      end if
    end if

  end function splitTokens

  integer function countTokens(this, string)
!@sum Sweep through string to count tokens.  This can then
!@+ be used to allocate the token array and a 2nd sweep
!@+ will fill.
!TO DO - logic with splitTokens() is duplicated. 
    type (Parser_type), intent(in) :: this
    character(len=*), intent(in) :: string

    character(len=len(string)) :: buffer
    integer :: numTokens
    integer :: idxNextSeparator, idxStart

    buffer = adjustl(string)
    numTokens = 0
    do
      if (len_trim(buffer) == 0) exit ! done
      numTokens = numTokens + 1
      ! skip strings which might have embedded separator characters
      idxStart = skipEmbeddedSeparators(buffer)
      idxNextSeparator = scan(buffer(idxStart:), this%tokenSeparators)
      if (idxNextSeparator > 0) then
        buffer = adjustl(buffer(idxStart + idxNextSeparator:))
      else
        exit
      end if
    end do

    countTokens = numTokens

  end function countTokens

  integer function skipEmbeddedSeparators(buffer) result(idxStart)
    character(len=*), intent(in) :: buffer
!@sum Location of end of string, which might have an embedded separator.
    if (scan(buffer, '"') == 1) then
      idxStart = scan(buffer(2:),'"') + 1 + 1 ! include 1st char
    else if (scan(buffer, "'") == 1) then
      idxStart = scan(buffer(2:),"'") + 1 + 1 ! include 1st char
    else
      idxStart = 1
    end if
  end function skipEmbeddedSeparators
    
  subroutine writeFormatted_parser(this, unit, dictionary)
!@sum Write dictionary as a text file.  Inverse of
!@+ parse().
    use AttributeDictionary_mod
    use AttributeHashMap_mod
    use AbstractAttribute_mod
    type (Parser_type), intent(in) :: this
    integer, intent(in) :: unit
    class (AttributeDictionary), intent(in) :: dictionary

    character(len=MAX_LEN_LINE) :: line
    type (AttributeHashMapIterator) :: iter
    class (AbstractAttribute), pointer :: t

    write(unit,'(a)') trim(this%beginData)

    iter = dictionary%begin()
    do while (iter /= dictionary%last())
      t => iter%value()
      line = trim(iter%key()) // ' = ' // trim(t%toString())
      write(unit,'(a)') trim(line)
      call iter%next()
    end do
    write(unit,'(a)') trim(this%endData)
    
 end subroutine writeFormatted_parser

end module PARSER_MOD
