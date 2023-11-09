#include "rundeck_opts.h"

      MODULE MODEL_COM
!@sum  MODEL_COM Main model variables, independent of resolution
!@auth Original Development Team
      use ModelClock_mod
      use TimeConstants_mod, only: INT_MONTHS_PER_YEAR,INT_DAYS_PER_YEAR
      use AbstractOrbit_mod, only: AbstractOrbit
      use AbstractCalendar_mod, only: AbstractCalendar
      IMPLICIT NONE
      SAVE

      CHARACTER*132 XLABEL !@var XLABEL=runID+brief description of run
      INTEGER :: LRUNID    !@var Run name stored in XLABEL(1:LRUNID)
      INTEGER :: LHEAD=15  !@var length of crucial beg of module_headers

!**** Model control parameters:
!@dbparam KOCEAN: if 0 => specified, if 1 => predicted ocean
      integer :: KOCEAN = 1

!**** Default simulation year. If set to zero, transient run.
!@dbparam master_yr year of simulation. This value will define aero_yr,
!@+       aer_int_yr, albsn_yr, crops_yr, ghg_yr, o3_yr, s0_yr, volc_yr,
!@+       variable_orb_par,orb_par_year_bp,and PI_run, unless any of these
!@+       are specifically defined.
      INTEGER ::  master_yr = 1951

!**** Diagnostic control parameters
!@dbparam NMONAV number of months in a diagnostic accuml. period
!@dbparam NIPRNT number of instantaneous initial printouts
      integer :: NMONAV=1, NIPRNT=1

C**** (Simplified) Calendar Related Terms
!@param JDperY,JMperY    number of days,months per year
!@var   JDendOfM(0:12)   last Julian day in month
!@var   JDmidOfM(0:13)   middle Julian day in month
      !integer, PARAMETER :: JDPERY = 365, JMPERY = 12  !Obsolete
      integer :: JDendOfM(0:INT_MONTHS_PER_YEAR) = (
     *     /0,31,59,90,120,151,181,212,243,273,304,334,365/)
      integer :: JDmidOfM(0:INT_MONTHS_PER_YEAR+1) = (
     *     /-15,16,45,75,106,136,167,197,228,259,289,320,350,381/)

!@var AMON,AMONTH(0:12)  (3-4 letter) names for current,all months
!@var AMON0  (3-4 letter) name of first month of the current acc-period
      CHARACTER*4 :: AMON='none',AMON0='none', AMONTH(0:12) = (/'IC  ',
     *  'JAN ','FEB ','MAR ','APR ','MAY ','JUNE',
     *  'JULY','AUG ','SEP ','OCT ','NOV ','DEC '/)

!@var NDAY and IYEAR1 relate CALENDAR TIME and INTERNAL TIME Itime :
!@var NDAY number of Internal Time Units per day (1 ITU = DTsrc sec)
!@nlparam IYEAR1  year 1 of internal clock (Itime=0 to 365*NDAY)
      INTEGER :: NDAY,IYEAR1=-1   !@var relate internal to calendar time

      class (AbstractOrbit), allocatable :: orbit
      class (AbstractCalendar), allocatable :: calendar

!@var modelEclock encapsulates current time with reference to a calendar
      type (ModelClock), public :: modelEClock
!@var modelEclockI encapsulates start time of model run
      type (ModelClock), public :: modelEClockI

!@var ITIME current time in ITUs (1 ITU = DTsrc sec, currently 1 hour)
      INTEGER :: Itime
!@var ItimeI,ItimeE   time at start,end of run
!@var Itime0          time at start of current accumulation period
!@var JMON0,JDATE0,JYEAR0,JHOUR0 date-info about Itime0 (beg.of acc.per)
      INTEGER :: ItimeI,ItimeE,   Itime0,JMON0,JDATE0,JYEAR0,JHOUR0

!@nlparam HOURI,DATEI,MONTHI,YEARI   start of model run
!@nlparam HOURE,DATEE,MONTHE,YEARE   end of model run
      INTEGER :: HOURI=0 , DATEI=1, MONTHI=1, YEARI=-1
     *          ,HOURE=0 , DATEE=1, MONTHE=1, YEARE=-1

!@dbparam DTSRC source time step (s)   = 1 ITU
      REAL*8 :: DTsrc = 3600.

!@dbparam KCOPY: if 1 => acc, if 3 => +od are saved
!@dbparam KRSF:  .rsf is written at beginning of every KRSF month
!@dbparam Ndisk:  DT_saversf    =  Ndisk *DTsrc fort.1/fort.2 saves
!@dbparam Nssw:   DT_checkSsw   =  Nssw  *DTsrc
      INTEGER :: KCOPY = 1, KRSF = 120, NDisk = 24, Nssw = 1

!**** Accounting variables
!@dbparam IRAND last seed used by rand.number generator
!@var KDISK next rsf (fort.)1 or 2 to be written to
      INTEGER :: IRAND=123456789, KDISK=1
!@param rsf_file_name names of restart files
      CHARACTER(6), PARAMETER :: rsf_file_name(4)=
     &         (/'fort.1','fort.2','fort.3','fort.4'/)
!@var MDYN,MCNDS,MRAD,MSURF,MDIAG,MELSE timing-indices
      INTEGER  MDYN,MCNDS,MRAD,MSURF,MDIAG,MELSE

!@param NSAMPL number of diagnostic sampling schemes
      INTEGER, PARAMETER :: NSAMPL = 12
!@var IDACC(NSAMPL) counters for diagn. accumulations
      INTEGER, DIMENSION(NSAMPL) :: IDACC


!**** IO read/write flags used by the io_xyz routines
!@param IOWRITE Flag used for writing normal restart files
!@param IOWRITE_SINGLE Flag used for saving diags in single precision
!@param IOWRITE_MON Flag used for saving restart part only (no diags)
!@param IOREAD Flag used for reading in (composite) restart files
!@param IOREADNT Flag used for reading in restart files (w/o tracers)
!@param IRSFIC Flag used for reading in restart part to start NEW run
!@param IRSFICNT Flag used for reading restart (w/o tracers) for NEW run
!@param IRSFICNO Flag used for reading restart (w/o ocean) for NEW run
!@param IRERUN Flag used for reading in restart part to extend OLD run
      INTEGER, PARAMETER :: ioread=1,ioread_acc=2,
     *     irerun=3,irsfic=4,irsficnt=5,ioreadnt=6,irsficno=7,
     *     ioread_nodiag=8,
     *     iowrite=-1,iowrite_single=-2,iowrite_mon=-3

!@nlparam QCHECK TRUE for running diagnostic checks
      LOGICAL :: QCHECK = .FALSE.

!@var stop_on TRUE stops the model (set with "kill -15 PID)
      LOGICAL :: stop_on = .FALSE.

! these do not belong here but are needed until IWRITE/KDIAG etc.
! are removed from the INPUTZ namelist (harder than it sounds)
      INTEGER :: IWRITE_sv,JWRITE_sv,ITWRITE_sv
      INTEGER, DIMENSION(13) :: KDIAG_sv

      type ModelE_Clock_type
        integer :: iTime
      end type ModelE_Clock_type

      contains

      ! Use rundeck parameters to determine which orbit to use
      function makeOrbit() result(orbit)
      use AbstractOrbit_mod
      use Dictionary_mod
      use DOMAIN_DECOMP_1d, only: am_i_root
      implicit none
      class (AbstractOrbit), allocatable :: orbit

      character(len=80) :: planetName

      planetName = 'Earth'       ! default
      call sync_param('planetName', planetName)

      select case (planetName)
      case ('Earth','earth','EARTH')

        if (AM_I_ROOT()) print*,'Using standard Earth orbit'
        allocate(orbit, source=makeEarthOrbit())

      case default

        allocate(orbit, source=makePlanetOrbit(planetName))
         
      end select

      ! Send orbit description to stdout
      if (am_i_root()) call orbit%print()

      end function makeOrbit

      function makeEarthOrbit() result(orbit)
      use AbstractOrbit_mod
      use Earth365DayOrbit_mod
      use ParameterizedEarthOrbit_mod
      use DOMAIN_DECOMP_1d, only: am_i_root
      use Dictionary_mod
      use Constant, only: planetParams
      implicit none
      class (AbstractOrbit), allocatable :: orbit

      integer :: variable_orb_par
      integer :: orb_par_year_bp
      real*8 :: eccen
      real*8 :: obliq
      real*8 :: omegt

      real*8 :: pYear

      if (is_set_param("variable_orb_par")) then
        call get_param( "variable_orb_par", variable_orb_par )
      else
        if (master_yr == 0) then
          variable_orb_par=1
        else
          variable_orb_par=0
        endif
      endif

      if (is_set_param("orb_par_year_bp")) then
        call get_param( "orb_par_year_bp", orb_par_year_bp )
      else
        if (master_yr == 0) then
          orb_par_year_bp=0
        else
          orb_par_year_bp=1950-master_yr
        endif
      endif

      select case (variable_orb_par)
      case (1) 
        allocate(orbit, 
     &        source=newParameterizedEarthOrbit(orb_par_year_bp))


      case (0)  ! orbital parameters fixed from year orb_par_year_bp
        pyear=1950.-orb_par_year_bp ! here "present" means "1950"
        allocate(orbit, source=Earth365DayOrbit(pYear))
        if (am_i_root()) then
          write(6,*) 'Fixed orbital parameters from year',pyear,' CE:'
        end if
      case (-1) ! orbital parameters fixed, directly set
        eccen = planetParams%getEccentricity()
        obliq = planetParams%getObliquity()
        omegt = planetParams%getLongitudeAtPeriapsis()
        allocate(orbit, source=Earth365DayOrbit(eccen, obliq, omegt))
        if (am_i_root()) then
          write(6,*) 'Orbital Parameters Specified:'
        end if
      case default  ! set from defaults (defined in CONSTANT module)
        allocate(orbit, source=Earth365DayOrbit())
      end select

      if (am_i_root()) then
        eccen = orbit%getEccentricity()
        obliq = orbit%getObliquity()
        omegt = orbit%getLongitudeAtPeriapsis()
      end if
      
      end function makeEarthOrbit

      function makePlanetOrbit(planetName) result(orbit)
      use PlanetaryOrbit_mod, only: PlanetaryOrbit
      use BaseTime_mod
      use TimeInterval_mod
      use Rational_mod
      use DOMAIN_DECOMP_1d, only: am_i_root
      use Dictionary_mod
      use Constant, only: planetParams
      type (PlanetaryOrbit) :: orbit
      character(len=*), intent(in) :: planetName

      real*8 :: eccentricity
      real*8 :: obliquity ! in degrees
      real*8 :: longitudeAtPeriapsis ! in degrees
      real*8 :: orbitalPeriod ! in seconds
      real*8 :: rotationPeriod ! in seconds
      real*8 :: meanDistance ! in A.U.
      real*8 :: s
      type (TimeInterval) :: secondsPerDay
      type (TimeInterval) :: secondsPerYear
      type (TimeInterval) :: period
      integer :: daysPerYear

      associate (p => planetParams)

      orbit = PlanetaryOrbit(p)

      secondsPerYear = TimeInterval(
     & Rational(p%getSiderealOrbitalPeriod(), tolerance=1.d-6))
      s = 1.d0 / 
     &     (1/p%getSiderealRotationPeriod() - 
     &     1/p%getSiderealOrbitalPeriod())
      secondsPerDay = TimeInterval( Rational(s, tolerance=1.d-6) )
      daysPerYear = max(1,nint(secondsPerYear / secondsPerDay))
      end associate

      if (AM_I_ROOT()) then
         write(*,*) 'Planet :: ' // trim(planetName)
         write(*,*)'Using planetary calendar:', s
         period = orbit%getSiderealRotationPeriod()
         write(*,*)'siderealRotationPeriod: ', real(period)
         period = orbit%getSiderealOrbitalPeriod()
         write(*,*)'siderealOrbitalPeriod: ', real(period)
         write(*,*)'meanDistance: ', orbit%getMeanDistance()
         write(*,*) '  Precession (degs from ve):',
     &        orbit%getLongitudeAtPeriapsis()
         write(*,*)'   Days per year: ', daysPerYear
         write(*,*)'   Seconds per day: ', real(secondsPerDay)
         write(*,*) '  Eccentricity:', orbit%getEccentricity()
         write(*,*) '  Obliquity (degs):',orbit%getObliquity()
      end if

      end function makePlanetOrbit

!TODO move to ModelClock class
      logical function isBeginningAccumPeriod(clock)
      use CalendarMonth_mod
      type (ModelClock) :: clock
      integer :: months
      type (CalendarMonth) :: cMonth

      integer :: month, day, year
      month = clock%getMonth()
      day = clock%getDayOfYear()
      year = clock%getYear()
      months=(year-Jyear0)*INT_MONTHS_PER_YEAR + month - JMON0

      cMonth = calendar%getCalendarMonth(month-1, year)

      isBeginningAccumPeriod = 
     &     clock%isBeginningOfDay() .and. 
     &     months.ge.NMONAV .and. 
     &     day.eq.1+cmonth%lastDayInMonth

      end function isBeginningAccumPeriod

      END MODULE MODEL_COM

      MODULE MDIAG_COM
!@sum  MDIAG_COM information common to all diagnostics
!@auth Original Development Team
      implicit none

      integer, parameter ::
     &     sname_strlen=30,units_strlen=30,lname_strlen=80

C**** Accumulating_period information
      INTEGER, DIMENSION(12) :: MONACC  !@var MONACC(1)=#Januaries, etc
      CHARACTER*12 :: ACC_PERIOD='PARTIAL'    !@var string MONyyr1-yyr2

!@param ia_cpl idacc-index currently associated with DTsrc, placed
!@+            here for visibility to non-atmospheric components
      integer, parameter :: ia_cpl=1 ! currently has to be == 1

!@dbparam make_timeaxis whether scaled monthly output files should contain
!@+       a time axis.  This option has not yet been introduced for all
!@+       diagnostic categories, and is not meaningful for some.
      logical, public :: make_timeaxis=.false.

      END MODULE MDIAG_COM

      subroutine reset_mdiag
!@sum reset info common to all diagnostics
      use model_com, only : idacc,modelEClock,
     &     itime,itime0,nday,iyear1,jyear0,jmon0,jdate0,jhour0,amon0
      implicit none
      integer jd0
      idacc(1:12)=0
      idacc(12)=1

      call modelEclock%get(year=jyear0, month=jmon0, dayOfYear=jd0,
     & date=jdate0, hour=jhour0, amn=amon0)
      itime0=itime

      return
      end subroutine reset_mdiag

      SUBROUTINE aPERIOD (JMON1,JYR1,months,years,moff,  aDATE,LDATE)
!@sum  aPERIOD finds a 7 or 12-character name for an accumulation period
!@+   if the earliest month is NOT the beginning of the 2-6 month period
!@+   the name will reflect that fact ONLY for 2 or 3-month periods
!@auth Reto A. Ruedy
      USE MODEL_COM, only : AMONTH
      use TimeConstants_mod, only: INT_MONTHS_PER_YEAR
      implicit none
!@var JMON1,JYR1 month,year of beginning of period 1
      INTEGER JMON1,JYR1
!@var JMONM,JMONL middle,last month of period
      INTEGER JMONM,JMONL
!@var months,years length of 1 period,number of periods
      INTEGER months,years
!@var moff = # of months from beginning of period to JMON1 if months<12
      integer moff
!@var yr1,yr2 (end)year of 1st and last period
      INTEGER yr1,yr2
!@var aDATE date string: MONyyr1(-yyr2)
      character*12 aDATE
      character(len=4) :: i5toc4 ! function in shared/Utilities.F90
!@var LDATE length of date string (7 or 12)
      INTEGER LDATE

      LDATE = 7                  ! if years=1
      if(years.gt.1) LDATE = 12

      aDATE(1:12)=' '
      aDATE(1:3)=AMONTH(JMON1)        ! letters 1-3 of month IF months=1
      yr1=JYR1
      JMONL=JMON1+months-1
      if(JMONL.GT.INT_MONTHS_PER_YEAR) then
         yr1=yr1+1
         JMONL=JMONL-INT_MONTHS_PER_YEAR
      end if
      if (moff.gt.0.and.months.le.3) then  ! earliest month is NOT month
        JMONL = 1 + mod(10+jmon1,INT_MONTHS_PER_YEAR)  ! 1 of the 2-3 month pd
        yr1=JYR1
        if (jmon1.gt.1) yr1=yr1+1
      end if
      yr2=yr1+years-1
      write(aDATE(4:7),'(a4)') i5toc4(yr1)
      if(years.gt.1) write(aDATE(8:12),'(a1,a4)') '-',i5toc4(yr2)

      if(months.gt.INT_MONTHS_PER_YEAR) aDATE(1:1)='x'       ! should not happen
      if(months.le.1 .or. months.gt.INT_MONTHS_PER_YEAR) return

!**** 1<months<13: adjust characters 1-3 of aDATE (=beg) if necessary:
!**** beg=F?L where F/L=letter 1 of First/Last month for 2-11 mo.periods
!****    =F+L                                        for 2 month periods
!****    =FML where M=letter 1 of Middle month       for 3 month periods
!****    =FnL where n=length of period if n>3         4-11 month periods
      aDATE(3:3)=AMONTH(JMONL)(1:1)            ! we know: months>1
      IF (months.eq.2) then
        aDATE(2:2)='+'
        return
      end if
      if (months.eq.3) then
        JMONM = JMONL-1
        if (moff.eq.1) jmonm = jmon1+1
        if (jmonm.gt.INT_MONTHS_PER_YEAR) then
          jmonm = jmonm - INT_MONTHS_PER_YEAR
        end if  
        if (jmonm.le.0 ) jmonm = jmonm + INT_MONTHS_PER_YEAR
        aDATE(2:2)=AMONTH(JMONM)(1:1)
        return
      end if
      if (moff.gt.0) then  ! can't tell non-consec. from consec. periods
        jmon1 = jmon1-moff
        if (jmon1.le.0) jmon1 = jmon1 + INT_MONTHS_PER_YEAR
        JMONL=JMON1+months-1
        if (jmonl.gt.INT_MONTHS_PER_YEAR) then
          jmonl = jmonl - INT_MONTHS_PER_YEAR
        end if
        aDATE(1:1)=AMONTH(JMON1)(1:1)
        aDATE(3:3)=AMONTH(JMONL)(1:1)
      end if
      IF (months.ge.4.and.months.le.9) write (aDATE(2:2),'(I1)') months
      IF (months.eq.10) aDATE(2:2)='X'         ! roman 10
      IF (months.eq.11) aDATE(2:2)='B'         ! hex   11
      IF (months.eq.6) THEN                    !    exceptions:
         IF (JMON1.eq. 5) aDATE(1:3)='NHW'     ! NH warm season May-Oct
         IF (JMON1.eq.11) aDATE(1:3)='NHC'     ! NH cold season Nov-Apr
      END IF
      IF (months.eq.7) THEN                    !    to avoid ambiguity:
         IF (JMON1.eq. 1) aDATE(1:3)='J7L'     ! Jan-Jul J7J->J7L
         IF (JMON1.eq. 7) aDATE(1:3)='L7J'     ! Jul-Jan J7J->L7J
      END IF
      IF (months.eq.INT_MONTHS_PER_YEAR) THEN
C****    beg=ANn where the period ends with month n if n<10 (except 4)
         aDATE(1:3)='ANN'                      ! regular annual mean
         IF (JMONL.le. 9) WRITE(aDATE(3:3),'(I1)') JMONL
         IF (JMONL.eq. 4) aDATE(1:3)='W+C'     ! NH warm+cold seasons
         IF (JMONL.eq.10) aDATE(1:3)='C+W'     ! NH cold+warm seasons
         IF (JMONL.eq.11) aDATE(1:3)='ANM'     ! meteor. annual mean
      END IF
      return
      end SUBROUTINE aPERIOD

      MODULE TIMINGS
!@sum  TIMINGS contains variables for keeping track of computing time
!@auth Gavin Schmidt
      IMPLICIT NONE
      SAVE
!@param NTIMEMAX maximum number of possible time accumulators
      INTEGER, PARAMETER :: NTIMEMAX=12
!@var NTIMEACC actual number of time accumulators
      INTEGER :: NTIMEACC = 0
!@var TIMING array that holds timing info
      REAL*8, DIMENSION(0:NTIMEMAX) :: TIMING
!@var TIMESTR array that holds timing info description
      CHARACTER*12, DIMENSION(NTIMEMAX) :: TIMESTR

      END MODULE TIMINGS

      SUBROUTINE SET_TIMER(STR,MINDEX)
!@sum  SET_TIMER sets an index of TIMING for a particular description
!@auth Gavin Schmidt
      USE TIMINGS
      IMPLICIT NONE
!@var STR string that describes timing accumulator
      CHARACTER*12, INTENT(IN) :: STR
!@var MINDEX index for that accumulator
      INTEGER, INTENT(OUT) :: MINDEX
      INTEGER N

C**** Check whether index has been set
      DO N=1,NTIMEACC
        IF (STR.EQ.TIMESTR(N)) THEN
          MINDEX=N
          RETURN
        END IF
      END DO
C**** Otherwise increase number of indexes
      NTIMEACC = NTIMEACC + 1
      IF (NTIMEACC.gt.NTIMEMAX) call stop_model(
     &     "Too many timing indices: increase NTIMEMAX",255)
      MINDEX = NTIMEACC
      TIMESTR(MINDEX) = STR
C****
      RETURN
      END SUBROUTINE SET_TIMER

      SUBROUTINE TIMER (NOW,MSUM)
!@sum  TIMER keeps track of elapsed CPU time in hundredths of seconds
!@auth Gary Russell
      USE TIMINGS
      USE GETTIME_MOD
      IMPLICIT NONE
      REAL*8, INTENT(OUT) :: NOW     !@var NOW current CPU time (seconds)
      INTEGER, INTENT(INOUT) :: MSUM !@var MSUM index for running total
      REAL*8 :: INC                  !@var INC time since last call
      REAL*8, SAVE :: LAST = 0       !@var LAST  last CPU time
      REAL*8 :: CMAX                 !@var CMAX max.count before 0-reset

      CALL GETTIME(NOW, CMAX)
      INC  = NOW - LAST
      if(inc<0) inc=inc+cmax ! offset system_clock reset
      TIMING(MSUM)  = TIMING(MSUM) + INC
      LAST = NOW
      RETURN
      END SUBROUTINE TIMER

      SUBROUTINE TIMEOUT (BEGIN,MIN,MOUT)
!@sum  TIMEOUT redistributes timing info between counters
!@auth Gary Russell
      USE TIMINGS
      USE GETTIME_MOD
      IMPLICIT NONE
!@var MBEGIN CPU time start of section (.01 s)
      REAL*8, INTENT(IN) :: BEGIN
      INTEGER, INTENT(INOUT) :: MIN  !@var MIN index to be added to
      INTEGER, INTENT(INOUT) :: MOUT !@var MOUT index to be taken from
      REAL*8 :: INC                  !@var INC time since MBEGIN
      REAL*8 :: NOW                  !@var NOW current CPU time (s)
      REAL*8 :: CMAX                 !@var CMAX max.count before 0-reset

      CALL GETTIME(NOW, CMAX)
      INC  = NOW - BEGIN
      if(inc<0) inc=inc+cmax ! offset system_clock reset
      TIMING(MIN)  = TIMING(MIN)  + INC
      TIMING(MOUT) = TIMING(MOUT) - INC
      RETURN
      END SUBROUTINE TIMEOUT

      subroutine getdte(It,Nday,Iyr0,Jyr,Jmn,Jd,Jdate,Jhour,amn)
!@sum  getdte gets julian calendar info from internal timing info
!@auth Gavin Schmidt
      use TimeConstants_mod, only : HOURS_PER_DAY, INT_DAYS_PER_YEAR
      USE MODEL_COM, only : amonth, calendar
      use CalendarMonth_mod
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: It,Nday,Iyr0
      INTEGER, INTENT(OUT) :: Jyr,Jmn,Jd,Jdate,Jhour
      CHARACTER*4, INTENT(OUT) :: amn
      type (CalendarMonth) :: cMonth

      Jyr=Iyr0+It/(Nday*INT_DAYS_PER_YEAR)
      Jd=1+It/Nday-(Jyr-Iyr0)*INT_DAYS_PER_YEAR
      Jmn=1
      cMonth = calendar%getCalendarMonth(Jmn, Jyr)

      do while (Jd.GT.cMonth%lastDayinMonth)
        Jmn=Jmn+1
        cMonth = calendar%getCalendarMonth(Jmn, Jyr)
      end do

      cMonth = calendar%getCalendarMonth(Jmn-1, Jyr)
      Jdate=Jd-cMonth%lastDayinMonth
      Jhour=mod(It*HOURS_PER_DAY/Nday,HOURS_PER_DAY)
      amn=amonth(Jmn)

      return
      end subroutine getdte

      SUBROUTINE DAILY_cal(end_of_day)
!@sum  DAILY performs daily tasks at end-of-day and maybe at (re)starts
!@auth Original Development Team
!@calls getdte
      use model_com, only: modelEclock
      USE MODEL_COM, only : itime,iyear1,nday,aMON
      use ModelClock_mod
      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: end_of_day   !!!!! NOT USED ?????
      integer :: year, month, day, hour, date

C****
C**** CALCULATE THE DAILY CALENDAR
C****
      call modelEclock%get(year=year, month=month, dayOfYear=day, 
     *     date=date, hour=hour, amn=amon)

      RETURN
      END SUBROUTINE DAILY_cal

#ifdef USE_ESMF
!-------------------------------------------------------------------------------
      subroutine init_esmf_clock_for_modelE(interval, clock)
!-------------------------------------------------------------------------------
      use TimeConstants_mod, only: HOURS_PER_DAY
      use MODEL_COM, only : itimei,itimee,nday,iyear1
      use ESMF
      implicit none
      integer :: interval ! timestep in seconds
      type (ESMF_clock)              :: clock

      type (ESMF_time) :: startTime
      type (ESMF_time) :: stopTime
      type (ESMF_timeinterval) :: timeStep
      type (ESMF_calendar) :: gregorianCalendar

      integer :: rc
      integer :: jday
      integer :: YEARI,MONTHI,DATEI,HOURI,MINTI
      integer :: YEARE,MONTHE,DATEE,HOURE,MINTE
      CHARACTER*4 :: cmon

      call getdte(itimei,nday,iyear1,YEARI,MONTHI,jday,DATEI,HOURI,cmon)
      MINTI = nint(mod( 
     &     mod(Itimei*HOURS_PER_DAY/Nday,HOURS_PER_DAY) * 60d0, 60d0))
      call getdte(itimee,nday,iyear1,YEARE,MONTHE,jday,DATEE,HOURE,cmon)
      MINTE = nint(mod( 
     &     mod(Itimee*HOURS_PER_DAY/Nday,HOURS_PER_DAY) * 60d0, 60d0))

    ! initialize calendar to be Gregorian type
      gregorianCalendar = esmf_calendarcreate(ESMF_CALKIND_GREGORIAN,
     &     name="GregorianCalendar", rc=rc)
      call stop_if_error(rc,'creating calendar')

      call ESMF_CalendarSetDefault(ESMF_CALKIND_GREGORIAN, rc=rc)
      call stop_if_error(rc,'creating calendar')

    ! initialize start time
      call ESMF_timeset(startTime,
     &     YY=YEARI,
     &     MM=MONTHI,
     &     DD=DATEI,
     &      H=HOURI,
     &      M=MINTI,
     &      S=0,
     &     calendar=gregorianCalendar, rc=rc)
      call stop_if_error(rc,'setting initial time')
!      write(*,*)'Time Set Start: ',STARTTIME

    ! initialize stop time
      call ESMF_timeset(stopTime,
     &     YY=YEARE,
     &     MM=MONTHE,
     &     DD=DATEE,
     &      H=HOURE,
     &      M=MINTE,
     &      S=0,
     &     calendar=gregorianCalendar, rc=rc)
      call stop_if_error(rc,'setting final time')
!      write(*,*)'Time Set End: ',ENDTIME

    ! initialize time interval
      call ESMF_timeintervalset(timeStep, S=INT(interval), rc=rc)
      call stop_if_error(rc,'setting time interval')

    ! initialize the clock with the above values
      clock = esmf_clockcreate(timeStep, startTime, stoptime=stopTime,
     &     name="ApplClock",rc=rc)
      call stop_if_error(rc,'creating clock')

      call ESMF_ClockSet ( clock, CurrTime=startTime, rc=rc)
      call stop_if_error(rc,'setting clock')

      contains

      subroutine stop_if_error(retcode,errmsg)
      integer :: retcode
      character(len=*) :: errmsg
      If (retcode /= ESMF_SUCCESS) call stop_model(
     &     'init_esmf_clock_for_modelE: '//trim(errmsg), 255)
      end subroutine stop_if_error

      end subroutine init_esmf_clock_for_modelE
#endif
