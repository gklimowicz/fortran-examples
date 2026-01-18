      module util
!@sum Utility routines for Ent_standalone and giss_LSM_standalone runs
!@+   Utilities for managing time variables.

      implicit none
      private
      save

      public YEARSEC, TimeDiff

      contains

      integer function YEARSEC(year) result(seconds_in_year)
!@sum YEARSEC  Seconds in a year for leap and non-leap years.
      use ent_const, only : sday
          integer :: year
          !integer, parameter :: sday = 86400 !second in a day

          if (IsLeapYear(year)) then
            seconds_in_year = 366*24*sday
          else
            seconds_in_year = 365*24*sday
          end if
      end function YEARSEC

      LOGICAL FUNCTION IsLeapYear(year) Result(Leap)
!@sum IsLeapYear  Return true if leap year.
!@+   Is leap year if YEAR is a number that is
!@+   exactly divisible by 4, except for century years which
!@+   must also be divisible by 400.
!
!     INPUT:
!         year: 4-digit number                  [I4]
!
!     OUTPUT:
!        TP_IsLeapYear: logical value           [L}
!
!      LIBRARY: [TIMPAK]
!       STATUS: [BETA]
!
!     AUTHOR: Robert D. Stewart
!             Pacific Northwest National Laboratory (PNNL)
!             PO Box 999, MSIN K3-55
!             Richland, WA 99352-0999  USA
!             trebor@purdue.edu
!             http://www.pnl.gov/berc/staff/rds.html
!
!     CREATION DATE: 19-APR-1999
!     REVISIONS HISTORY:
!
!     COMMENTS:
!
      integer,intent(in) :: year
!      logical :: Leap

      IF (MOD(year,4).EQ.0) THEN
        !Possible leap year
        IF (MOD(year,100).EQ.0) THEN
          !Test for special case
          IF (MOD(year,400).EQ.0) THEN
            Leap=.true.
          ELSE
            Leap=.false.
          ENDIF
        ELSE
          Leap=.true.
        ENDIF
      ELSE
        Leap=.false.
      ENDIF

      END FUNCTION IsLeapYear

      
      !************************************************************************
      integer function JulianDay(time) Result(jday)
!@sum JulianDay  Julian day of a leap or non-leap year given time data struct.
      use ent_types, only : timestruct
      type(timestruct),intent(in) :: time
      integer :: m

      integer, parameter :: daysinmonth(12) 
     &     = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
      integer, parameter :: daysinmonthleap(12) 
     &     = (/31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)      

      jday = 0.0
      if (IsLeapYear(time%year)) then
        if (time%month.gt.1) then
          do m=1,(time%month-1)
            jday = jday + daysinmonthleap(m)
          end do
        endif
      else 
        if (time%month.gt.1) then
          do m=1,(time%month-1)
            jday = jday + daysinmonth(m)
          end do
        endif
      endif
      jday = jday + time%day
      end function JulianDay

      !************************************************************************

      integer function TimeDiff(time, prevtime) Result(dtsec)
!@sum TimeDiff Time difference between time and prevtime in seconds
      use ent_types, only : timestruct

      type(timestruct), intent(in) :: time, prevtime
      integer :: jday, jdayprev, hoursec, hoursecprev
      integer, parameter :: sday = 86400  !seconds in a day

      jday = JulianDay(time)
      jdayprev = JulianDay(prevtime)
      hoursec = time%hour*3600 + time%minute*60.0 + time%seconds
      hoursecprev = prevtime%hour*3600 + prevtime%minute*60.0
     &     + prevtime%seconds
      
      if (prevtime%year.le.time%year) then
        dtsec = (jday-jdayprev)*sday + hoursec - hoursecprev
      else
        dtsec = (jdayprev*sday + YEARSEC(prevtime%year))-jday*sday 
     &       + hoursec - hoursecprev
      end if
      end function TimeDiff

      !************************************************************************
!      integer function get_year(YYYYMMDDHHMMSS) Result(year)
!          integer :: YYYYMMDDHHMMSS
!          year =  FLOOR(YYYYMMDDHHMMSS/10000000000)
!      end function get_year
      !************************************************************************
!      integer function get_month(YYYYMMDDHHMMSS) Result(month)
!          integer :: YYYYMMDDHHMMSS
!          month = FLOOR(FRACTION(YYYYMMDDHHMMSS/10000000000)*100)
!      end function get_month
      !************************************************************************
!      integer function get_day(YYYYMMDDHHMMSS) Result(day)
!          integer :: YYYYMMDDHHMMSS
!          day = FLOOR(FRACTION(YYYYMMDDHHMMSS/1000000000000)*100)
!      end function get_day
      !************************************************************************
!      real function get_hourfrac(YYYYMMDDHHMMSS) Result(hourfrac)
!          integer :: YYYYMMDDHHMMSS
!          real :: minutes, sec
!          hourfrac = FRACTION(YYYYMMDDHHMMSS/100000000)*100
!          minutes = FLOOR(FRACTION(hourfrac)*100)
!          sec = FRACTION(hourfrac*10000)*100
!          hourfrac = FLOOR(hourfrac)
!          hourfrac = hourfrac + minutes/60 + seconds/3600
!      end function get_hourfrac
      !************************************************************************

      END MODULE util
