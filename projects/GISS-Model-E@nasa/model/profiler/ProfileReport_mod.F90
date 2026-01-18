module ProfileReport_mod
   use ReportColumn_mod
   use TimerList_mod
   implicit none
   private

   public :: ProfileReport_type
   public :: newReport
   public :: addColumn
   public :: getHeader
   public :: getUnits
   public :: getLine
   public :: generateReport
   public :: delete


   public :: generateParallelReport

   public :: MAX_RECORD_LENGTH

   type ProfileReport_type
      private
      integer :: recordWidth = 0
      type (ReportColumn_type), pointer :: columns(:) => null()
   end type ProfileReport_type

   interface generateReport
      module procedure generateReportAsStringDefault
      module procedure generateReportAsString
   end interface

   interface generateParallelReport
      module procedure generateParallelReportDefault
      module procedure generateParallelReportFromList
   end interface

   interface getHeader
      module procedure getHeader_report
   end interface

   interface getUnits
      module procedure getUnits_report
   end interface

   interface getLine
      module procedure getLine_report
   end interface

   interface delete
      module procedure delete_report
   end interface

   integer, parameter :: MAX_RECORD_LENGTH = 132
contains

   function newReport() result(report)
      type (ProfileReport_type) :: report
      report%recordWidth = 0
      allocate(report%columns(0))
   end function newReport

   subroutine addColumn(this, column)
      type (ProfileReport_type) :: this
      type (ReportColumn_type), intent(in) :: column
      integer :: n

      call extend(this%columns)
      n = getNumColumnns(this)
      this%columns(n) = column
      this%recordWidth = this%recordWidth + getFieldWidth(column) + 1

   contains

      subroutine extend(columns)
         type (ReportColumn_type), pointer :: columns(:)
         type (ReportColumn_type), allocatable :: tmp(:)
         integer :: n

         n = size(columns)
         allocate(tmp(n))
         tmp = columns
         deallocate(columns)
         allocate(columns(n+1))
         columns(:n) = tmp
         deallocate(tmp)
      end subroutine extend
   end subroutine addColumn

   function getHeader_report(this) result(record)
      type (ProfileReport_type) :: this
      character(len=this%recordWidth) :: record

      integer :: i, j0, j1

      record = ''
      j0 = 1
      do i = 1, getNumColumnns(this)
         j1 = j0 + getFieldWidth(this%columns(i)) - 1
         record(j0:j1) = getHeader(this%columns(i))
         record(j1+1:j1+1) = ' '
         j0 = j1 + 2
      end do

   end function getHeader_report
      
   function getUnits_report(this) result(record)
      type (ProfileReport_type) :: this
      character(len=this%recordWidth) :: record

      integer :: i, j0, j1

      record = ''
      j0 = 1
      do i = 1, getNumColumnns(this)
         j1 = j0 + getFieldWidth(this%columns(i)) - 1
         record(j0:j1) = getUnits(this%columns(i))
         record(j1+1:j1+1) = ' '
         j0 = j1 + 2
      end do

   end function getUnits_report

   function getLine_report(this, list, lineNumber) result(record)
      use Timer_mod
      type (ProfileReport_type) :: this
      type (TimerList_type), intent(in) :: list
      integer, intent(in) :: lineNumber
      character(len=this%recordWidth) :: record

      integer :: i, j0, j1

      record = ''
      j0 = 1
      do i = 1, getNumColumnns(this)
         j1 = j0 + getFieldWidth(this%columns(i)) - 1
         record(j0:j1) = getField(this%columns(i), getName(list,lineNumber), getTimer(list,lineNumber))
         record(j1+1:j1+1) = ' '
         j0 = j1 + 2
      end do
   end function getLine_report

   integer function getNumColumnns(this) result(numColumns)
      type (ProfileReport_type), intent(in) :: this
      numColumns = size(this%columns)
   end function getNumColumnns

   function generateReportAsStringDefault(this) result(report)
      type (ProfileReport_type), intent(in) :: this
      character(len=MAX_RECORD_LENGTH), pointer ::  report(:)
      report => generateReport(this, getDefaultList())
   end function generateReportAsStringDefault

   function generateParallelReportDefault(this, communicator) result(report)
      type (ProfileReport_type), intent(in) :: this
      integer, intent(in) :: communicator
      character(len=MAX_RECORD_LENGTH), pointer ::  report(:)
      report => generateParallelReport(this, getDefaultList(), communicator)
   end function generateParallelReportDefault

   function generateParallelReportFromList(this, timerList, communicator) result(report)
      type (ProfileReport_type), intent(in) :: this
      type (TimerList_type), intent(in) :: timerList
      integer, intent(in) :: communicator
      character(len=MAX_RECORD_LENGTH), pointer ::  report(:)

      type (TimerList_type) :: globalList

      globalList = gather(timerList, communicator)
      report => generateReport(this, globalList)

   end function generateParallelReportFromList

   function generateReportAsString(this, timerList) result(report)
      type (ProfileReport_type), intent(in) :: this
      type (TimerList_type), intent(in) :: timerList
      character(len=MAX_RECORD_LENGTH), pointer ::  report(:)

      integer :: i, n

      n = getNumTimers(timerList)
      allocate(report(5+n))
      report(1)(:) = repeat("-", this%recordWidth)
      report(2)(:) = getHeader(this)
      report(3)(:) = getUnits(this)
      report(4)(:) = repeat("-", this%recordWidth)
      do i = 1, n
         report(4+i)(:) = getLine(this, timerList, i)
      end do
      report(5+n)(:) = repeat("-", this%recordWidth)

   end function generateReportAsString

   subroutine writeReportToFile(this, timerList, ioUnit)
      type (ProfileReport_type), intent(in) :: this
      type (TimerList_type), intent(in) :: timerList
      integer, intent(in) :: ioUnit
      character(len=this%recordWidth), pointer :: report(:)
      integer :: i

      report => generateReport(this, timerList)
      do i = 1, size(report)
         write(ioUnit,'(a)') trim(report(i))
      end do
      deallocate(report)

   end subroutine writeReportToFile

   subroutine delete_report(this)
      type (ProfileReport_type), intent(inOut) :: this
      deallocate(this%columns)
   end subroutine delete_report

end module ProfileReport_mod
