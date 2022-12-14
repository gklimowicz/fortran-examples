#include "rundeck_opts.h"

module TracerSurfaceSource_mod
  use TracerSource_mod
  use timestream_mod, only : timestream
  implicit none
  private

  public :: TracerSurfaceSource
  public :: initSurfaceSource
  public :: readEmissionHeader
  public :: readMonthly
  public :: readSurfaceSource
  public :: bracketDay
  public :: bracketYear

  public :: PARSE_SUCCESS
  public :: PARSE_MISSING_FIELD
  public :: PARSE_UNEXPECTED_VALUE

  type, extends(TracerSource) :: TracerSurfaceSource
    character(len=30) :: sourceName ! holds source name, read from file header, e.g. to be
    ! placed into sname arrays.
    character(len=30) :: sourceLname ! holds long source name
!@var EMstream interface for reading and time-interpolating emissions file if it is netcdf.
!@+   See usage notes in timestream_mod.
    type (timestream) :: EMstream
    character(len=10) :: tracerName ! tracer name read from emis file header (should match trname)
    character(len=1) :: frequency ! (annual 'a' or monthly 'm') read from emis file header
    character(len=1) :: resolution ! horiz. resolution (C? M? F?) read from emis file header
    character(len=9) :: years
    integer :: yearStart ! starting year/decade for a transient emissions file
    integer :: yearEnd   ! ending year/decade for a transient emissions file
    integer :: yearStep  ! interval between records in a transient emissions file
    logical :: firstTrip = .true.

    integer :: monthA = -1 ! first month for the current interpolation
    real*8, allocatable :: month1cache(:,:) ! used for interpolating from file source
    real*8, allocatable :: month2cache(:,:)
    logical :: saveCache = .true. ! set to false to use less memory (but more frequent reads)
  end type TracerSurfaceSource

  integer, parameter :: PARSE_SUCCESS = 0
  integer, parameter :: PARSE_MISSING_FIELD = 1
  integer, parameter :: PARSE_UNEXPECTED_VALUE = 2

contains

  subroutine initSurfaceSource(this, tracerName, fileName, sectorNames, checkname)
    use SystemTools, only : stLinkStatus,stFileList
    use TracerSource_mod, only: N_MAX_SECT
    use Dictionary_mod, only : sync_param
    USE FILEMANAGER, only: openunit,closeunit,is_fbsa
    use pario, only : par_open,par_close,read_attr
    USE DOMAIN_DECOMP_ATM, only: GRID
    use TimeConstants_mod, only: HOURS_PER_DAY
    use SpecialIO_mod, only: write_parallel,read_parallel
    type (TracerSurfaceSource), intent(inout) :: this
    character(len=*), intent(in) :: tracerName
    character(len=*), intent(in) :: fileName
    character(len=300) :: out_line
    character*10, intent(in):: sectorNames(:)
    logical, intent(in) :: checkName
    logical :: diurnalFileExists = .false.

    integer :: nsect, nn, i, j, iu, fid
    integer :: linkstatus, nfiles, ifile, ios, jyr
    integer, parameter :: max_fname_len=128
    character(len=max_fname_len), allocatable :: flist(:)
    character(len=max_fname_len) :: thisline
    character(len=max_fname_len+8) :: fileToRead
    character(len=4) :: c4
    character*32 :: pname
    character*35 :: fname
    character*124 :: tr_sectors_are
    integer :: numTrSectors
    character(len=80) :: name ! sector
    real*8 :: sumDiurnal
    real*8, parameter :: diurnalSumTolerance=1.d-4
    character*80 :: targetVariable

    if(is_fbsa(fileName)) then ! binary file. Use old method.
      call openunit(fileName,iu,.true.)
      call readEmissionHeader(this, tracerName, iu, checkname)
      call closeunit(iu)
    else
      fileToRead=fileName ! default (e.g. if file is not a directory, or
                          ! the directory search doesn't find a good file)
      call stLinkStatus(trim(fileName), linkstatus)
      if(linkstatus==2) then  ! this is a directory. todo: no hard-coded retcodes
        ! The file is a directory. Do something similar to subroutine check_metadata
        ! in timestream_mod to determine any useable YYYY.nc file in the directory.
        allocate(flist(1000)) ! 1000 files maximum
        call stFileList(trim(fileName),flist,nfiles)
        do ifile=1,nfiles
          thisline = adjustl(flist(ifile))
          if(len_trim(thisline).ne.7) cycle
          if(thisline(5:7).ne.'.nc') cycle
          c4 = thisline(1:4)
          read(c4,*,iostat=ios) jyr
          if(ios.ne.0) cycle
          if(jyr.lt.0) cycle
          ! acceptable file. define it and exit:
          fileToRead=trim(fileName)//'/'//trim(thisline)
          exit
        end do
        deallocate(flist)
      end if ! directory search

      ! continue reading netCDF file:
      this%tracerName = tracerName
      this%sourceName = 'notfound'
      fid = par_open(grid,trim(fileToRead),'read')
      ! First try to read the variable attribute to get source name:
      call read_attr(grid,fid,this%tracerName,'source',i,this%sourceName)
      ! If that fails, look for the source attribute of the variable
      ! that varname attribute points to (like init_stream would):
      if(trim(this%sourceName).eq.'notfound') then
        targetVariable=this%tracerName
        call read_attr(grid,fid,'global',trim(this%tracerName)//'name',&
        & i,targetVariable)
        call read_attr(grid,fid,trim(targetVariable),'source',&
        & i,this%sourceName)
      endif
      ! If that fails, look for a global source attribute:
      if(trim(this%sourceName).eq.'notfound') then
        call read_attr(grid,fid,'global','source',i,this%sourceName)
      endif
      ! If even that fails, stop the model:
      call par_close(grid,fid)
      if(trim(this%sourceName).eq.'notfound') then
        call stop_model('source name not found in file '//trim(fileName),255)
      endif
    endif

    ! append ' source' to the long name, and '_src' to the short name
    this%sourceLname = trim(this%sourceName)//' source'
    this%sourceName = trim(this%sourceName)//'_src'

    ! -- begin sector stuff --
    tr_sectors_are = ' '
    pname=trim(trim(fileName)//'_sect')
    call sync_param(pname,tr_sectors_are)
    numTrSectors = 0

    i=1
    do while(i < len(tr_sectors_are))
      j=index(tr_sectors_are(i:len(tr_sectors_are))," ")
      if (j > 1) then
        numTrSectors = numTrSectors + 1
        i=i+j
      else
        i=i+1
      end if
    enddo
    if(numTrSectors > n_max_sect)  &
         &     call stop_model("num_tr_sectors problem",255)
    this%num_tr_sectors = numTrSectors

    if(numTrSectors > 0) then
      read(tr_sectors_are,*) this%tr_sect_name(1:numTrSectors)

      do nsect=1, numTrSectors
        name = trim(this%tr_sect_name(nsect))
        this%tr_sect_index(nsect) = 0
        loop_nn: do nn=1, size(sectorNames)
          if(trim(name) == trim(sectorNames(nn))) then
            this%tr_sect_index(nsect) = nn
            exit loop_nn
          endif
        enddo loop_nn
      enddo
    endif

    ! -- begin diurnal stuff -- 
    fname=trim('diurnal_'//trim(fileName))
    ! governed by file existance:
    inquire(file=trim(fname), exist=diurnalFileExists)
    if(diurnalFileExists)then
       this%applyDiurnalCycle=.true.
       write(out_line,*)'Applying diurnal cycle to file '//fileName
       call write_parallel(trim(out_line))
       call openunit(fname,iu,.false.,.true.)
       call read_parallel(this%diurnalCycle,iu)
       ! check that the diurnal cycle's sum is close to the number of hours
       ! in a day (meaning it's hourly average would be a factor of 1.):
       sumDiurnal=SUM(this%diurnalCycle)
       if(  sumDiurnal > HOURS_PER_DAY + diurnalSumTolerance  &
     & .or. sumDiurnal < HOURS_PER_DAY - diurnalSumTolerance) then
         write(out_line,*) &
     &   trim(fname),' sum is ',sumDiurnal,' not',HOURS_PER_DAY
         call write_parallel(trim(out_line))
         call stop_model('Problem with emissions diurnal cycle.',255)
       end if
       call closeunit(iu)
    end if

  end subroutine initSurfaceSource

  subroutine parseHeader(this, name, str,error)
!@sum parse_header gets the informantion from emissions file's
!@+  header and reports back the meta-data. 
!@auth Greg Faluvegi
!@auth Tom Clune (significant refactoring)

    type (TracerSurfaceSource), intent(inout) :: this
    character(len=*), intent(in) :: name
    integer :: error
    character*(*) str
    character(len=16) :: yearStepStr

    integer :: rc
    character(len=9) :: years

    rc = parseField(str, 'freq', this%frequency, ['a','m','d'])
    call checkRc(rc, 2, 3, error)

    rc = parseField(str, 'name', this%tracerName, [name])
    call checkRc(rc, 4, 5, error)

    rc = parseField(str, 'source', this%sourceName)
    call checkRc(rc, 6, 6, error)

#ifndef CUBED_SPHERE
    rc = parseField(str, 'res', this%resolution, ['M','F','C'])
#else
    rc = parseField(str, 'res', this%resolution)
#endif
    call checkRc(rc, 7, 8, error)

    rc = parseField(str, 'y', years)
    if (rc == PARSE_MISSING_FIELD) then
      this%yearStart = 0
      this%yearEnd = 0
    else
      read(years(1:4),'(I4)') this%yearStart
      read(years(6:9),'(I4)') this%yearEnd

      if(this%yearStart /= this%yearEnd)then
        if(this%yearStart < 0 .or. this%yearStart > 3000) error=9
        if(this%yearEnd   < 0 .or. this%yearEnd   > 3000) error=9
      endif
    endif

    if(this%frequency == 'd')then
      ! daily files, while non-transient, must have years defined:
      if(this%yearStart==0 .or. this%yearEnd==0)error=12
      this%yearStep=1  ! default= 1-day steps
    else
      this%yearStep=10 ! default=decades for backwards compatability
    endif

    rc = parseField(str, 'del', yearStepStr)
    if (rc /= PARSE_MISSING_FIELD) read(yearStepStr,*) this%yearStep
    if(this%yearStep == 0) then
      error=10
    else
      ! check for integer number of slices:
      if(this%frequency == 'd')then
        if (this%yearStep /= 1) error=13
      else
        if (mod((this%yearEnd-this%yearStart), this%yearStep) /=0. ) &
           & error=11
      endif
    endif

  end subroutine parseHeader

  integer function parseField(string, field, value, expected) result(rc)
!@sum Given a string from the header, parse it to extract a value from the given field, and
!@+ optionally check to ensure that it is among the expected results.
    use StringUtilities_mod, only: toLowerCase
    character(len=*), intent(inout) :: string
    character(len=*), intent(in) :: field
    character(len=*), intent(out) :: value
    character(len=*), intent(in), optional :: expected(:)

    integer :: n1, n2

    n1 = scan(string,'=')
    if (n1 == 0) then
      rc = PARSE_MISSING_FIELD
      return
    end if
    n2 = scan(string,' ')
    rc = PARSE_SUCCESS

    if (string(1:n1-1) /= trim(field)) rc = PARSE_MISSING_FIELD

    read(string(n1+1:n2-1),*) value
    string = string(n2+1:)
    if (present(expected)) then
      if (all(toLowerCase(trim(value)) /= toLowerCase(expected(:)))) &
           &          rc = PARSE_UNEXPECTED_VALUE
    end if

  end function parseField

  subroutine checkRc(rc, errMissing, errUnexpected, error)
!@sum A bit of a kludge to allow different error codes for each field in the header.
    integer, intent(in) :: rc
    integer, intent(in) :: errMissing
    integer, intent(in) :: errUnexpected
    integer, intent(inout) :: error

    select case (rc)
    case (PARSE_SUCCESS)
      error = 0
    case (PARSE_MISSING_FIELD)
      error = errMissing
    case (PARSE_UNEXPECTED_VALUE)
      error = errUnexpected
    end select
  end subroutine checkRc

  subroutine readEmissionHeader(this, name, iu,checkname)
!@sum readEmissionHeader reads the emissions file's header and
!@+   reports back the meta-data. 
!@auth Greg Faluvegi

    use SpecialIO_mod, only: write_parallel
    USE FileManager, only: openunit,closeunit

    implicit none

    type (TracerSurfaceSource), intent(inout) :: this
    character(len=*), intent(in) :: name
    integer, intent(in) :: iu
    integer :: error
    character*80 :: message,header
    character(len=300) :: out_line
    logical :: checkname

    error=0
    read(iu)header
    if(len_trim(header) < 1)error=1

    call parseHeader(this, name, trim(header), error)

    select case(error)
    case default ! nothing
    case(1) ; message='readEmissionHeader: missing header'
    case(2) ; message='readEmissionHeader: problem with freq'
    case(3) ; message='readEmissionHeader: a,m,d are choices for freq'
    case(4) ; message='readEmissionHeader: problem with tracer name'
    case(5) ; message='readEmissionHeader: tracer name mismatch'
    case(6) ; message='readEmissionHeader: problem with source'
    case(7) ; message='readEmissionHeader: problem with res'
    case(8) ; message='readEmissionHeader: M and F are choices for res'
    case(9) ; message='readEmissionHeader: transient years seem wrong'
    case(10); message='readEmissionHeader: yearStep(e.g. kstep) is zero'
    case(11); message='readEmissionHeader: trans yrs step/years suspect'
    case(12); message='readEmissionHeader: daily emis years undefined'
    case(13); message='readEmissionHeader: 1-day steps expected for now'
    end select
    if(error > 0) then
      if(error == 5 .and. .not. checkname)then
        continue
      else
        write(out_line,*) trim(header)
        call write_parallel(trim(out_line))
        write(out_line,*) trim(message)
        call write_parallel(trim(out_line))
        call stop_model('problem reading emissions',255)
      endif
    endif

  end subroutine readEmissionHeader

  !TODO - pass in CLOCK object rather than year, day.
  subroutine readMonthly(this, iu, data, xyear, xday, grid)
!@sum Read in monthly sources and interpolate to current day
!@+   Input: this, the object containing the source metadata; 
!@+          iu, the fileUnit;
!@+          xyear, xday - year and day of year
!@+   Output: interpolated data array + two monthly data arrays
!@+   "this" object caches the data to reduce I/O
!@auth Greg Faluvegi, Jean Lerner and others
!@auth Tom Clune - refactored into this module and added caching

    use SpecialIO_mod, only: readt_parallel
    use SpecialIO_mod, only: write_parallel
    use MpiSupport_mod, only: am_i_root

    use domain_decomp_atm, only : DIST_GRID, getDomainBounds

    implicit none

    type (TracerSurfaceSource), intent(inout) :: this
    type (DIST_GRID), intent(in) :: grid
    real*8, DIMENSION(GRID%I_STRT_HALO:GRID%I_STOP_HALO, &
         &                  GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: &
         &     data
    integer ::  iu
    integer, intent(in) :: xyear, xday

    integer :: J_0, J_1, I_0, I_1

    call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1)
    call getDomainBounds(grid, I_STRT=I_0, I_STOP=I_1)

    if (isTransient(this)) then
      call readMonthlyTransient(xday, xyear)
    else 
      call readMonthlyNontransient(xday)
    endif

    return

  contains

    subroutine readMonthlyNontransient(dayOfYear)
      integer, intent(in) :: dayOfYear
      integer :: monthA, monthB
      character(len=300) :: out_line
      real*8 :: frac

      !**** Interpolate two months of data to current day
      call bracketDay(dayOfYear, monthA, monthB, frac)

      if (monthA /= this%monthA) then ! new month
        this%monthA = monthA
        call updateCache(this, I_0, I_1, J_0, J_1)
        call accumulate(this%month1cache, monthA, 1.d0)
        call accumulate(this%month2cache, monthB, 1.d0)
      endif

      data(I_0:I_1,J_0:J_1) = &
           &   this%month1Cache(:,:)*frac +  this%month2Cache(:,:)*(1-frac)

      write(out_line,*) trim(this%tracerName),' ', trim(this%sourceName), &
           &     ' interp to now ', frac

      if (.not. this%saveCache) call releaseCache(this)
      call write_parallel(trim(out_line))

    end subroutine readMonthlyNontransient

    subroutine readMonthlyTransient(dayOfYear, year)
      integer, intent(in) :: dayOfYear
      integer, intent(in) :: year
      integer :: monthA, monthB
      integer :: yearA, yearB
      character(len=300) :: out_line

      integer :: offset1, offset2

      real*8 :: alpha ! interpolation parameter for (decadal) periods
      real*8 :: frac  ! interpolation parameter for months

      call bracketYear(this, dayOfYear, year, offset1, offset2, alpha, yearA, yearB)

      call bracketDay(dayOfYear, monthA, monthB, frac)

      data = 0
      if (monthA /= this%monthA) then ! new month
        this%monthA = monthA
        call updateCache(this, I_0, I_1, J_0, J_1)

        call accumulate(this%month1cache, monthA+offset1, 1-alpha)
        call accumulate(this%month1cache, monthA+offset2, alpha)

        call accumulate(this%month2cache, monthB+offset1, 1-alpha)
        call accumulate(this%month2cache, monthB+offset2, alpha)
      end if

      ! note that the cache's are not halo'd
      data(I_0:I_1,J_0:J_1) = &
           &     this%month1cache(:,:)*frac + this%month2cache(:,:)*(1-frac)

      if (.not. this%saveCache) call releaseCache(this)

      write(out_line,'(a,1X,a,a4,F9.4,a21,I4,a13,I4,a22,F9.4)') &
           &     trim(this%tracerName),trim(this%sourceName),' at ', &
           &     100.d0*alpha, &
           &     '% of period this day ',yearA,' to this day ',yearB, &
           &     ' and monthly fraction=',frac
      call write_parallel(trim(out_line))

    end subroutine readMonthlyTransient

    subroutine accumulate(data, record, weight)
      use SpecialIO_mod, only: rewind_parallel
      use FileManager, only : nameUnit
      real*8, intent(inout) :: data(:,:)
      integer, intent(in) :: record
      real*8 :: weight
      real*8 :: tmp(grid%i_strt_halo:grid%i_stop_halo, &
           &     grid%j_strt_halo:grid%j_stop_halo)

      call rewind_parallel(iu)
      call skipHeader(iu)
      call readt_parallel(grid,iu,nameunit(iu), tmp, record)
      data(:,:) = data(:,:) + weight * tmp(I_0:I_1,J_0:J_1)

    end subroutine accumulate

  end subroutine readMonthly

  !
  ! Bracket the two months for the given day of year
  !
  subroutine bracketDay(dayOfYear, monthA, monthB, frac)
    use TimeConstants_mod, only: INT_MONTHS_PER_YEAR
    use JulianCalendar_mod, only: MID_JULIAN_DAY_IN_MONTH
    integer, intent(in) :: dayOfYear
    integer, intent(out) :: monthA
    integer, intent(out) :: monthB
    real*8, intent(out) :: frac

    if (dayOfYear <= MID_JULIAN_DAY_IN_MONTH(1)) then ! dayOfYear in Jan 1-15, first month is Dec
      monthA = 0
    else                      ! dayOfYear is after Jan 16
      do monthA = 1, INT_MONTHS_PER_YEAR
        if (dayOfYear <= MID_JULIAN_DAY_IN_MONTH(monthA+1)) exit
      end do
    end if
    monthB = 1 + monthA

    frac = float(MID_JULIAN_DAY_IN_MONTH(monthB)-dayOfYear) / &
         &     (MID_JULIAN_DAY_IN_MONTH(monthB)-MID_JULIAN_DAY_IN_MONTH(monthB-1))

    ! adjust to 12 month cycle
    monthA = 1 + modulo(monthA - 1, INT_MONTHS_PER_YEAR)
    monthB = 1 + modulo(monthB - 1,  INT_MONTHS_PER_YEAR)

  end subroutine bracketDay

  !
  ! Bracket year by years contained in file
  !
  subroutine bracketYear(this, dayOfYear, year, offset1, offset2, alpha, yearA, yearB)
    use TimeConstants_mod, only: INT_MONTHS_PER_YEAR
    use TimeConstants_mod, only: INT_DAYS_PER_YEAR
    type (TracerSurfaceSource), intent(in) :: this
    integer, intent(in) :: dayOfYear
    integer, intent(in) :: year
    integer, intent(out) :: offset1
    integer, intent(out) :: offset2
    real*8, intent(out) :: alpha
    integer, intent(out) :: yearA
    integer, intent(out) :: yearB

    integer, parameter :: HALF_YEAR = 1 + (INT_DAYS_PER_YEAR-1)/2
    integer :: period1, period2, numPeriods

    numPeriods = 1 + (this%yearEnd - this%yearStart) / this%yearStep

    if (isBefore(dayOfYear,year,HALF_YEAR,this%yearStart)) then
      ! use start year
      yearA = this%yearStart
      yearB = this%yearStart
      period1 = 1
      period2 = 1
      alpha = 0.0d0
    elseif (isBefore(HALF_YEAR,this%yearEnd,dayOfYear,year)) then
      ! use final year
      yearA = this%yearEnd
      yearB = this%yearEnd
      period1 = numPeriods
      period2 = numPeriods
      alpha = 1.0d0
    else ! bracket year 
      do period1 = 1, numPeriods - 1
        yearA = this%yearStart + (period1-1)*this%yearStep
        yearB = this%yearStart + (period1)*this%yearStep
        if (isBefore(dayOfYear,year,HALF_YEAR,yearB)) exit
      end do
      period2 = 1 + period1
      alpha = real(year - yearA) / this%yearStep
    end if
    offset1 = (period1-1) * INT_MONTHS_PER_YEAR
    offset2 = (period2-1) * INT_MONTHS_PER_YEAR
  end subroutine bracketYear

  logical function isBefore(day1, year1, day2, year2)
    integer, intent(in) :: day1, year1
    integer, intent(in) :: day2, year2
    isBefore = (year1 < year2) .or. (year1 == year2 .and. day1 < day2)
  end function isBefore

  subroutine skipHeader(iu)
    use domain_decomp_atm, only : AM_I_ROOT
    integer, intent(in) :: iu
    if (AM_I_ROOT()) read(iu)
  end subroutine skipHeader

  logical function isTransient(this)
    type (TracerSurfaceSource), intent(in) :: this
    isTransient = (this%yearStart /= this%yearEnd)
  end function isTransient

  subroutine updateCache(this, i0, i1, j0, j1)
    type (TracerSurfaceSource), intent(inout) :: this
    integer, intent(in) :: i0, i1, j0, j1
    if (.not. allocated(this%month1cache)) then
      allocate(this%month1cache(i0:i1,j0:j1))
      allocate(this%month2cache(i0:i1,j0:j1))
    end if
    this%month1cache = 0
    this%month2cache = 0
  end subroutine updateCache

  subroutine releaseCache(this)
    type (TracerSurfaceSource), intent(inout) :: this

    if (allocated(this%month1cache)) then
      deallocate(this%month1cache)
      deallocate(this%month2cache)
    end if

  end subroutine releaseCache

  subroutine readSurfaceSource(tracerName, this, fname, checkname, sfc_src, xyear, xday, isChemTracer)
    USE DOMAIN_DECOMP_ATM, only: GRID,  readt_parallel, write_parallel
    use Domain_decomp_atm, only: getDomainBounds
    USE FILEMANAGER, only: openunit,closeunit, nameunit,is_fbsa
    use TimeConstants_mod, only: EARTH_DAYS_PER_YEAR
    use timestream_mod, only : init_stream,read_stream
    use dictionary_mod, only : get_param
    character(len=*), intent(in) :: tracerName
    type (TracerSurfaceSource), intent(inout) :: this
    character(*), intent(in) :: fname
    logical, intent(in) :: checkname
    real*8, intent(inout) :: sfc_src(grid%i_strt_halo:,grid%j_strt_halo:)
    integer, intent(in) :: xyear, xday
    logical, intent(in) :: isChemTracer

    integer :: iu,k,ipos,kx,iposDay,kstep=10
    character(len=300) :: out_line
    real*8 :: alpha
    real*8, dimension(GRID%I_STRT_HALO:GRID%I_STOP_HALO, &
         &                  GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: &
         & sfc_a,sfc_b

    INTEGER :: J_1, J_0, I_0, I_1
    integer :: cyclic_yr,master_yr,nc_emis_use_ppm_interp

    if(.not.is_fbsa(fname)) then

      if(this%firstTrip) then
        this%firstTrip = .false.
        call get_param('master_yr',master_yr)
        if (isChemTracer) then
          call get_param('o3_yr',cyclic_yr,default=master_yr)
          select case (tracerName)
          case ('NOx')
            call get_param('NOx_yr',cyclic_yr,default=cyclic_yr)
          case ('CO')
            call get_param('CO_yr',cyclic_yr,default=cyclic_yr)
          case ('Alkenes', 'Paraffin')
            call get_param('VOC_yr',cyclic_yr,default=cyclic_yr)
          end select
        else
          call get_param('aer_int_yr',cyclic_yr,default=master_yr)
          select case (tracerName)
          case ('SO2', 'SO4', 'M_ACC_SU', 'M_AKK_SU', 'ASO4__01')
            call get_param('SO2_int_yr',cyclic_yr,default=cyclic_yr)
          case ('NH3')
            call get_param('NH3_int_yr',cyclic_yr,default=cyclic_yr)
          case ('BCII', 'BCB', 'M_BC1_BC', 'M_BOC_BC', 'AECOB_01')
            call get_param('BC_int_yr',cyclic_yr,default=cyclic_yr)
          case ('OCII', 'OCB', 'M_OCC_OC', 'M_BOC_OC', 'AOCOB_01',&
                'vbsAm2', 'vbsAm1', 'vbsAz', 'vbsAp1', 'vbsAp2',&
                'vbsAp3', 'vbsAp4', 'vbsAp5', 'vbsAp6')
            call get_param('OC_int_yr',cyclic_yr,default=cyclic_yr)
          end select
        end if
        cyclic_yr=ABS(cyclic_yr)
        call get_param('nc_emis_use_ppm_interp',nc_emis_use_ppm_interp,&
          & default=1)
        if (nc_emis_use_ppm_interp==1) then
          call init_stream(grid,this%EMstream,trim(fname), &
             trim(this%tracername),0d0,1d30,'ppm',xyear,xday, &
             cyclic = (cyclic_yr > 0) )
        else
          call init_stream(grid,this%EMstream,trim(fname), &
             trim(this%tracername),0d0,1d30,'linm2m',xyear,xday, &
             cyclic = (cyclic_yr > 0) )
        endif
      endif
      call read_stream(grid,this%EMstream,xyear,xday,sfc_src)

    else

    CALL getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1)
    CALL getDomainBounds(grid, I_STRT=I_0, I_STOP=I_1)

    ! open file and read its header:
    call openunit(fname,iu,.true.)
    call skipHeader(iu)

    ! now read the data: (should be in kg/m2/s please)

    ! -------------- non-transient emissions ----------------------------!
    if (this%yearStart==this%yearEnd.or.this%frequency=='d') then 

      select case (this%frequency)
      case ('a')        ! annual file, only read first time
        if (this%firstTrip) then
          call readt_parallel(grid,iu,fname,sfc_src(:,:),1)
          write(out_line,*)trim(this%tracerName), &
               &        ' ',trim(this%sourceName),&
               &             ' ann source read.'
          call write_parallel(trim(out_line))
          this%firstTrip = .false.
        endif
      case ('m')        ! monthly file, interpolate to now
        call readMonthly(this,iu,sfc_src(:,:), xyear,xday, grid)
      case ('d')        ! daily file, don't interpolate
        if(xyear < this%yearStart .or. xyear > this%yearEnd)then
          write(out_line,*)'Year ',xyear,' out of range of daily ', &
     &    'tracer source file years: ',this%yearStart,this%yearEnd,'.'
          call write_parallel(trim(out_line))
          call stop_model('xyear bad for daily emis reading',255)
        endif
        iposDay=NINT(EARTH_DAYS_PER_YEAR*(xyear-this%yearStart)+xday)
        call readt_parallel(grid,iu,fname,sfc_src(:,:),iposDay)
      end select

      ! --------------- transient emissions -------------------------------!
    else                
      select case (this%frequency)
      case ('a')        ! annual file, only read first time + new steps
        kstep=this%yearStep
        ipos=1
        alpha=0.d0 ! before start year, use start year value
        kx=this%yearStart ! just for printing
        if(xyear>this%yearEnd .or.&
             &      (xyear==this%yearEnd.and.xday>=183))then
          alpha=1.d0 ! after end year, use end year value     
          ipos=(this%yearEnd-this%yearStart)/kstep
          kx=this%yearEnd-kstep
        endif
        do k=this%yearStart,this%yearEnd-kstep,kstep
          if(xyear>k .or. (xyear==k.and.xday>=183)) then
            if(xyear<k+kstep.or.(xyear==k+kstep.and.xday<183))then
              ipos=1+(k-this%yearStart)/kstep ! (integer artithmatic)
              alpha=(EARTH_DAYS_PER_YEAR*(0.5+real(xyear-1-k))+xday) / &
                   &      (EARTH_DAYS_PER_YEAR*real(kstep))
!              alpha = real(365*(xyear-k) + xday-183,kind=8) / real(365*kstep,kind=8)
              kx=k
              exit
            endif
          endif
        enddo

        call readt_parallel(grid,iu,fname,sfc_a(:,:),ipos)
        call readt_parallel(grid,iu,fname,sfc_b(:,:),1)

        sfc_src(I_0:I_1,J_0:J_1)=sfc_a(I_0:I_1,J_0:J_1)* &
             &      (1.d0-alpha)+sfc_b(I_0:I_1,J_0:J_1)*alpha

        write(out_line,'(a,1X,a,a4,F9.4,a16,I4,a8,I4)') &
             &      trim(this%tracerName),&
             &           trim(this%sourceName),' at ',&
             &      100.d0*alpha,'% of period mid ',kx,' to mid ',kx+kstep
        call write_parallel(trim(out_line))

      case ('m')        ! monthly file, interpolate to now
        call readMonthly(this,iu,sfc_src(:,:),xyear,xday, grid)
      case ('d')
        call stop_model('Transient, daily tracer src not allowed.',255)
      end select

    endif
    call closeunit(iu)

    endif ! giss format or not

  end subroutine readSurfaceSource

end module TracerSurfaceSource_mod
