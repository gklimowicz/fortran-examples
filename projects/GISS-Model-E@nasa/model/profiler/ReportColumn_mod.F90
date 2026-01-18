module ReportColumn_mod
   use Timer_mod
   use TimeFormatUtilities_mod
   implicit none
   private

   public :: ReportColumn_type
   public :: newColumn
   public :: getFieldWidth
   public :: getHeader
   public :: getUnits
   public :: getField
   public :: setScale
   public :: setPrecision
   
   integer, parameter :: MAX_FIELD_WIDTH = 20
   type ReportColumn_type
      private
      real(kind=r64) :: scale = 1.
      integer :: columnType = 0
      integer :: fieldWidth = 0
      integer :: numDigits  = 0
      integer :: precision  = 2
      character(len=MAX_FIELD_WIDTH) :: header = ' '
      character(len=MAX_FIELD_WIDTH) :: units  = ' '
   end type ReportColumn_type

   integer, parameter, public :: NAME_COLUMN           = 1
   integer, parameter, public :: INCLUSIVE_TIME_COLUMN = 2
   integer, parameter, public :: EXCLUSIVE_TIME_COLUMN = 3
   integer, parameter, public :: SPACER_COLUMN         = 4
   integer, parameter, public :: TRIP_COUNTS_COLUMN    = 5
   integer, parameter, public :: MAXIMUM_TIME_COLUMN   = 6
   integer, parameter, public :: MINIMUM_TIME_COLUMN   = 7
   integer, parameter, public :: AVERAGE_TIME_COLUMN   = 8
   integer, parameter, public :: MIN_PROCESS_COLUMN    = 9
   integer, parameter, public :: MAX_PROCESS_COLUMN    = 10


   interface getHeader
      module procedure getHeader_column
   end interface

   interface getUnits
      module procedure getUnits_column
   end interface

contains

   function newColumn(columnType, fieldWidth) result(column)
      integer, intent(in) :: columnType
      integer, intent(in) :: fieldWidth
      type (ReportColumn_type) :: column

      column%columnType = columnType
      column%fieldWidth = fieldWidth

      select case (columnType)
      case (NAME_COLUMN)
         column%header = 'Name'
         column%units  = ''
         column%numDigits = 0 ! unused
      case (INCLUSIVE_TIME_COLUMN)
         column%header = 'Inclusive'
         column%units  = 'seconds'
         column%numDigits = 2
      case (EXCLUSIVE_TIME_COLUMN)
         column%header = 'Exclusive'
         column%units  = 'seconds'
         column%numDigits = 2
      case (TRIP_COUNTS_COLUMN)
         column%header = 'Trips'
         column%units  = ''
         column%numDigits = 0
      case (MAXIMUM_TIME_COLUMN)
         column%header = 'Maximum'
         column%units  = 'seconds'
         column%numDigits = 2
      case (MINIMUM_TIME_COLUMN)
         column%header = 'Minimum'
         column%units  = 'seconds'
         column%numDigits = 2
      case (AVERAGE_TIME_COLUMN)
         column%header = 'Average'
         column%units  = 'seconds'
         column%numDigits = 2
      case (MIN_PROCESS_COLUMN)
         column%header = 'MIN PROCESS'
         column%units  = ''
         column%numDigits = 0
      case (MAX_PROCESS_COLUMN)
         column%header = 'MAX PROCESS'
         column%units  = ''
         column%numDigits = 0
      case default
         column%header = 'unknown'
      end select

   end function newColumn

   integer function getFieldWidth(this) result(fieldWidth)
      type (ReportColumn_type), intent(in) :: this

      fieldWidth = this%fieldWidth
   end function getFieldWidth
   
   function getHeader_column(this) result(header)
      type (ReportColumn_type), intent(in) :: this
      character(len=this%fieldWidth) :: header
      header = center(this%header(1:this%fieldWidth))
   end function getHeader_column

   function getUnits_column(this) result(units)
      type (ReportColumn_type), intent(in) :: this
      character(len=this%fieldWidth) :: units
      units = center(this%units(1:this%fieldWidth))
   end function getUnits_column

   function center(string) result(centeredString)
      character(len=*), intent(in) :: string
      character(len=len(string)) centeredString

      integer :: shift
      integer :: n

      n = len(string)
      shift = (n - len_trim(string))/2
      centeredString = repeat(" ", shift) // string(1:n-shift)

   end function center

   function getField(this, name, timer) result(field)
      type (ReportColumn_type), intent(in) :: this
      character(len=*), intent(in) :: name
      type (Timer_type), intent(in) :: timer
      character(len=this%fieldWidth) :: field

      integer :: numDigitsLHS
      integer :: iQuantity
      real(kind=r64) :: quantity
      character(len=10) :: fmt

      if (this%columnType == NAME_COLUMN) then
         field = name
      elseif (any(this%columnType == (/ TRIP_COUNTS_COLUMN, MIN_PROCESS_COLUMN, MAX_PROCESS_COLUMN /))) then
         write(fmt,'("(I",I2.2,")")') this%fieldWidth
         select case (this%columnType)
         case (TRIP_COUNTS_COLUMN)
            iQuantity = getNumTrips(timer)
         case (MIN_PROCESS_COLUMN)
            iQuantity = getMinProcess(timer)
         case (MAX_PROCESS_COLUMN)
            iQuantity = getMaxProcess(timer)
         end select
         write(field,fmt) iQuantity
      else
         select case (this%columnType)
         case (INCLUSIVE_TIME_COLUMN)
            quantity = getInclusiveTime(timer)
         case (EXCLUSIVE_TIME_COLUMN)
            quantity = getExclusiveTime(timer)
         case (MAXIMUM_TIME_COLUMN)
            quantity = getMaximumTime(timer)
         case (MINIMUM_TIME_COLUMN)
            quantity = getMinimumTime(timer)
         case (AVERAGE_TIME_COLUMN)
            quantity = getAverageTripTime(timer)
         case default
            quantity = -1
         end select
         numDigitsLHS = this%fieldWidth - this% precision - 1
         field = formatSeconds(quantity * this%scale, numDigitsLHS, this%precision)
      end if
   end function getField

   subroutine setScale(this, scale, units)
      type (ReportColumn_type), intent(inOut) :: this
      real(kind=r64), intent(in) :: scale
      character(len=*), intent(in) :: units
      this%units = units
      this%scale = scale
   end subroutine setScale

   subroutine setPrecision(this, precision)
      type (ReportColumn_type), intent(inOut) :: this
      integer, intent(in) :: precision
      this%precision = precision
   end subroutine setPrecision

end module ReportColumn_mod
