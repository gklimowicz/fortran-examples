#include "rundeck_opts.h"

! Soures of this doc:
! https://simplex.giss.nasa.gov/gcm/doc/nlparams.txt
!
! ISTART controls how the model first picks up the initial conditions to
! start the run. We have a number of options depending on how much
! information is already available.
! 
!     Cold Starts
!     ============
!
!     ISTART=2
!         Observed start. This sets atmospheric values to observations,
!         based on a particular format of AIC file. As for ISTART=1,
!         input files are required for ground and ocean variables.
!
!
!     Restart from checkpoint files
!     =============================
!
!     Checkpoint files, called ``fort.1.nc`` and ``fort.2.nc``, are written
!     frequently.  In case ModelE is terminated prematurely (i.e. before
!     the end time specified in the rundeck / I file), they may be used
!     to seamlessly continue the run where it left off.
!
!     Premature termination may happen for a number of reasons:
!        a) The run is stopped gracefully by the user.
!        b) The run exceeds its time limit on the supercomputer and is
!           stopped forcefully
!        c) The supercomputer crashes, and all jobs are stopped forcefully.
!        d) ModelE has a bug that causes it to crash.
!
!     Checkpoint files may be corrupted or otherwise unreadable, if
!     ModelE terminates while writing them.  For that reason, ModelE
!     alternates between writing to the names ``fort.1.nc`` and
!     ``fort.2.nc``.
!
!     ISTART=10   (DEFAULT if not specified in I file)
!         This is used internally to pick up from the checkpoint
!         file (the later of fort.1 and fort.2). Does not ever need to
!         be set in the rundeck.
! 
!     ISTART=11
!     ISTART=12
!         This is used internally to pick up from an instantaneous rsf
!         file (fort.1.nc for ISTART=11 or fort.2.nc for ISTART=12).
!
!     ISTART=13
!         Restarts from the OLDEST checkpoint file; the reverse of
!         ISTART=10
!
!     ISTART=14
!         Restarts from (presumably symlinked file) named fort.4.nc
!         Should link to fort.1.nc or fort.2.nc written on a previous
!         run.
!
!     Restart from .rsf files
!     =======================
!
!     .rsf files are written at the beginning of every KRSF months.
!     They contain everything 
!     in the checkpoint files EXCEPT diagnostic accumulation status
!
!     ISTART=9
!         Continuation of an old run that was stopped at the beginning
!         of a diagnostic accumulation period.  Use this to restart from
!         .rsf files:
!            a) Set AIC=myrestartfile.rsf
!            b) Set ISTART=9
!
!     Perturbation Experiments
!     ======================== 
!
!     ISTART=8
!         This is a restart from a model configuration identical to the
!         run now starting. This is for perturbation experiments, etc.
!
!         Start of a new run - parameters from rundeck+defaults
!         In particular: Itime is set to ItimeI, radiation and all
!         diagnostic accumulations are performed in the first hour,
!         IRAND is set to its default (unless reset in the rundeck).
!
!         Note: Since itime_tr0 defaults to Itime, tracers will be
!             reinitialized. To have them keep their setttings from the
!             rsf file, set itime_tr0 to < ItimeI (e.g. 0) for all
!             tracers in the rundeck parameters. If you use tracers that
!             depend on (Itime-itime_tr0), you need to set itime_tr0 to
!             the starting time of the rsf file for continuity.
! 
!
!     Obsolete ISTART Values
!     ======================
!
!     ISTART=1 (OBSOLETE)
!         Default start. This sets atmospheric variables to constants
!         and requires input files for ground values (a GIC file), and
!         ocean values (OIC) if required.  ISTART=1 may still work - if
!         I remember correctly, it was used a long time ago for
!         benchmarking when we were asked to submit a version that
!         needed no input files. It may still be useful for simpler
!         versions of the model, maybe for a different planet or
!         simplified earth (e.g. no topography, all desert, ...).
! 
!     ISTART=3-7 (OBSOLETE)
!     ----------
!         These were reserved for starting up a more complex model from
!         the state obtained by spinning up a simpler model, e.g. a
!         coupled model from an atmospheric model, a tracer run from a
!         run without tracers, etc. I'm not sure whether those options
!         are still needed or can be achieved without using the ISTART
!         parameter. They were kind of place holders to deal with
!         changes in the model restart file.

!     ISTART=3 (OBSOLETE)
!         Not used.
! 
!     ISTART=4 (OBSOLETE)
!         A restart from an rsf file from a previous run, but the ocean
!         is reinitialised. Needs an initial OIC file (for fully coupled
!         models).
! 
!     ISTART=5 (OBSOLETE)
!         A restart from an rsf file from a previous run, but no
!         tracers. This is only useful for tracer runs that need to be
!         initialised with a particular model state.
! 
!     ISTART=6 (OBSOLETE)
!         A restart from an rsf file from a previous run that might not
!         have had the same land-ocean mask. This makes sure to reset
!         snow values, pbl values and ocean values accordingly.
! 
!     ISTART=7 (OBSOLETE)
!         A restart from an rsf file from a previous run with the same
!         land-ocean mask. This still makes sure to set snow values and
!         ocean values. This is used mainly for converted model II'
!         data.
! 
!     ISTART<0 (OBSOLETE) This option is used by the post-processing
!         program to run the model to generate nice diagnostics. This
!         should never need to be set manually.  ISTART<0 may still work
!         if the model is run with "old I/O" but is not needed with "new
!         I/O". It was meant as a device to bridge the transition period
!         from old to new I/O.

      subroutine GISS_modelE(qcRestart,coldRestart,iFile,max_wall_time)
!@sum  MAIN GISS modelE main time-stepping routine
!@auth Original Development Team
!@ver  2009/05/11 (Based originally on B399)
!@var max_wall_time Time (in seconds) this ModelE is to run for.
!@+      Once it notices it has exceeded this allottment, it will wind
!@+      down and exit.
      USE FILEMANAGER, only : openunit,closeunit
      USE TIMINGS, only : ntimemax,ntimeacc,timing,timestr
      USE Dictionary_mod
      Use Parser_mod
      USE MODEL_COM, only: modelEclock
     &     , ItimeI, Itime, Ndisk
     &     , Jyear0, JMON0, Iyear1, ItimeE, Itime0
     &     , NIPRNT, XLABEL, LRUNID, MELSE, Nssw, stop_on
     &     , iowrite_single, isBeginningAccumPeriod
     &     , KCOPY,KRSF, NMONAV, IRAND, iowrite_mon, MDIAG, NDAY
     &     , rsf_file_name, iowrite, KDISK, dtSRC, MSURF
     &     , calendar
      USE DOMAIN_DECOMP_1D, only: AM_I_ROOT,broadcast,sumxpe
      USE RANDOM
      USE GETTIME_MOD
      USE MDIAG_COM, only : monacc,acc_period
#ifdef USE_MPP
      USE fms_mod,         only : fms_init, fms_end
#endif
#ifdef USE_FVCORE
      USE FV_INTERFACE_MOD, only: fvstate
      USE FV_INTERFACE_MOD, only: Checkpoint,Compute_Tendencies
#endif
      use TimeConstants_mod, only: SECONDS_PER_MINUTE,
     &                             INT_MONTHS_PER_YEAR
      use TimerPackage_mod, only: startTimer => start
      use TimerPackage_mod, only: stopTimer => stop
      use SystemTimers_mod
      use Timer_mod, only : getWTime   ! Tells time in seconds
      use seaice_com, only : si_ocn,iceocn ! temporary until precip_si,
      use fluxes, only : atmocn,atmice     ! precip_oc calls are moved
      use CalendarMonth_mod, only: LEN_MONTH_ABBREVIATION
      use Constant, only: initializeConstants
#ifdef CACHED_SUBDD
      USE SUBDD_MOD, only : write_monthly_files,write_daily_files,
     &     days_per_file,write_one_file
#endif
      implicit none
C**** Command line options
      logical, intent(in) :: qcRestart
      logical, intent(in) :: coldRestart
      character(len=*), intent(in) :: iFile
      integer :: max_wall_time

      INTEGER K,M,MSTART,MNOW,months,ioerr,Ldate,istart
      INTEGER :: MDUM = 0

      character(len=80) :: filenm

      REAL*8, DIMENSION(NTIMEMAX) :: PERCENT
      REAL*8, DIMENSION(0:NTIMEMAX) ::TIMING_glob = 0.
      REAL*8 start,now, DTIME,TOTALT

      CHARACTER aDATE*14, i5toc4*4 ! function in shared/Utilities.F90
      CHARACTER*8 :: string_go='___GO___'      ! green light
      CHARACTER*8 :: str
      integer :: iflag=1
      external sig_stop_model
#ifdef USE_GDB_FOR_FPE_BACKTRACE
      external sig_exception
#endif
      logical :: start9

      integer :: iu_IFILE
      real*8 :: tloopbegin, tloopend
      integer :: hour, month, day, date, year
      character(len=LEN_MONTH_ABBREVIATION) :: amon
      real*8 :: wtime0, wtime1   ! Start and end wall times (rank 0 only)

#ifdef CACHED_SUBDD
      character(len=8) :: yyyymmdd
#endif

#ifdef USE_SYSUSAGE
      do i_su=0,max_su
        call sysusage(i_su,0)
      enddo
#endif

C****
C**** Reading rundeck (I-file) options
C****
      call openunit(trim(ifile),iu_IFILE,.false.,.true.)
      call parse_params(iu_IFILE)
      call closeunit(iu_IFILE)

      call initializeModelE()

      ! Only the root node pays attention to allotted wall time
      if (AM_I_ROOT()) then
        wtime0 = getWTime()
        wtime1 = wtime0 + max_wall_time
      end if

C****
C**** INITIALIZATIONS
C****
         CALL TIMER (NOW,MDUM)

C**** Read input/ic files
      CALL INPUT (istart,ifile,coldRestart)

C**** Set run_status to "run in progress"
      call write_run_status("Run in progress...",1)

      IF (AM_I_ROOT()) Then
         open(3,file='flagGoStop',form='FORMATTED',status='REPLACE')
         write (3,'(A8)') string_go
         close (3)
      END IF
      call sys_signal( 15, sig_stop_model )  ! works only on single CPU
#ifdef USE_GDB_FOR_FPE_BACKTRACE
      call sys_signal( 8, sig_exception )
#endif
      START=NOW
      DO M=1,NTIMEACC
        START= START-TIMING(M)
      END DO

      call modelEclock%get(hour=hour, date=date, year=year,amn=amon)

      if (AM_I_ROOT())
     *   WRITE (6,'(A,11X,A4,I5,A5,I3,A4,I3,6X,A,I4,I10)')
     *   '0NASA/GISS Climate Model (re)started',
     *   'Year', year, aMon, date, ', Hr', hour,
     *   'Internal clock: DTsrc-steps since 1/1/',Iyear1,ITIME

         CALL TIMER (NOW,MELSE)

      call sys_flush(6)

C****
C**** MAIN LOOP
C****
      call gettime(tloopbegin)
      start9 = (istart == 9)

      main_loop: DO WHILE (Itime.lt.ItimeE)
        call startTimer('Main Loop')

      if (Ndisk > 0) then
        if (mod(Itime-ItimeI,Ndisk).eq.0 .or. start9) then
         start9 = .false.
         call checkpointModelE()
         call timer(NOW,MELSE)
        END IF
      end if

      if (modelEclock%isBeginningOfDay()) then
        call startNewDay()
      end if

#ifdef CACHED_SUBDD
      call set_subdd_period()
#endif

      call atm_phase1
      atmocn%updated=.true.

C****
C**** SURFACE INTERACTION AND GROUND CALCULATION
C****
C**** NOTE THAT FLUXES ARE APPLIED IN TOP-DOWN ORDER SO THAT THE
C**** FLUXES FROM ONE MODULE CAN BE SUBSEQUENTLY APPLIED TO THAT BELOW
C****
C**** APPLY PRECIPITATION TO SEA/LAKE/LAND ICE
      call startTimer('Surface')
      CALL PRECIP_SI(si_ocn,iceocn,atmice)  ! move to ocean_driver
      CALL PRECIP_OC(atmocn,iceocn)         ! move to ocean_driver

C**** CALCULATE SURFACE FLUXES (and, for now, this procedure
C**** also drives "surface" components that are on the atm grid)
      CALL SURFACE
      call stopTimer('Surface')
         CALL TIMER (NOW,MSURF)

      call ocean_driver

! phase 2 changes surf pressure which affects the ocean
      call atm_phase2

#ifdef CACHED_SUBDD
      if(write_one_file .and. itime+1.eq.itimee) then ! run finished
        filenm = 'allsteps.subdd'//XLABEL(1:LRUNID)
      elseif(write_daily_files .and.
     &     mod(itime+1,days_per_file*nday).eq.0) then
        write(yyyymmdd,'(i4,i2.2,i2.2)') year,month,date
        filenm=yyyymmdd//'.subdd'//XLABEL(1:LRUNID)
      else
        filenm = ''
      endif
      if(filenm.ne.'') call write_subdd_accfile (filenm)
#endif

C****
C**** UPDATE Internal MODEL TIME AND CALL DAILY IF REQUIRED
C****
      call modelEclock%nextTick()
      call modelEclock%get(year=year, month=month, dayOfYear=day,
     &     date=date, hour=hour, amn=amon)
      Itime=Itime+1                       ! DTsrc-steps since 1/1/Iyear1

      if (modelEclock%isBeginningOfDay()) THEN ! NEW DAY
        months=(year-Jyear0)*INT_MONTHS_PER_YEAR + month - JMON0
        call startTimer('Daily')
        call dailyUpdates
        call TIMER (NOW,MELSE)
        call stopTimer('Daily')
      end if                                  !  NEW DAY

#ifdef USE_FVCORE
! Since dailyUpdates currently adjusts surf pressure,
! moving this call to the atm driver will change results.
! todo 2: fold this into the fv run procedure.
       Call Compute_Tendencies(fvstate)
#endif

C****
C**** CALL DIAGNOSTIC ROUTINES
C****
        call startTimer('Diagnostics')

C**** PRINT CURRENT DIAGNOSTICS (INCLUDING THE INITIAL CONDITIONS)
      IF (NIPRNT.GT.0) THEN
        acc_period='PARTIAL      '
        filenm='PARTIAL.acc'//XLABEL(1:LRUNID)
        call io_rsf (filenm,Itime,iowrite_single,ioerr)
        call print_diags(1)
        NIPRNT=NIPRNT-1
        call set_param( "NIPRNT", NIPRNT, 'o' )
      END IF

C**** THINGS TO DO BEFORE ZEROING OUT THE ACCUMULATING ARRAYS
C**** (after the end of a diagn. accumulation period)
      if (isBeginningAccumPeriod(modelEClock)) then

C**** PRINT DIAGNOSTIC TIME AVERAGED QUANTITIES
        call aPERIOD (JMON0,JYEAR0,months,1,0, aDATE(1:12),Ldate)
        acc_period=aDATE(1:12)
        WRITE (aDATE(8:14),'(A3,a4)') aMON(1:3),i5toc4(year)
        call print_diags(0)
C**** SAVE ONE OR BOTH PARTS OF THE FINAL RESTART DATA SET
        IF (KCOPY.GT.0) THEN
C**** KCOPY > 0 : SAVE THE DIAGNOSTIC ACCUM ARRAYS IN SINGLE PRECISION
          monacc = 0
          do k=JMON0,JMON0+NMONAV-1
            m = k
            if(m.gt.INT_MONTHS_PER_YEAR) m = m - INT_MONTHS_PER_YEAR
            monacc(m) = monacc(m) + 1
          end do
          filenm=aDATE(1:7)//'.acc'//XLABEL(1:LRUNID)
          call io_rsf (filenm,Itime,iowrite_single,ioerr)
#ifdef CACHED_SUBDD
          if(write_monthly_files) then
            filenm=aDATE(1:7)//'.subdd'//XLABEL(1:LRUNID)
            call write_subdd_accfile (filenm)
          endif
#endif
        EndIf  !  (KCOPY > 0)
!**** KRSF > 0 : ALSO SAVE THE RESTART INFORMATION
        If (KRSF > 0)  Then
          If (Modulo(YEAR*INT_MONTHS_PER_YEAR+MONTH-1,KRSF) == 0)  Then
            CALL RFINAL (IRAND)
            call set_param( "IRAND", IRAND, 'o' )
            filenm='1'//aDATE(8:14)//'.rsf'//XLABEL(1:LRUNID)
            call io_rsf(filenm,Itime,iowrite_mon,ioerr)
#if defined( USE_FVCORE )
            call Checkpoint(fvstate, filenm)
#endif
          END IF
        EndIf  !  (KRSF > 0)

C**** PRINT AND ZERO OUT THE TIMING NUMBERS
        CALL TIMER (NOW,MDIAG)
        CALL SUMXPE(TIMING, TIMING_glob, increment=.true.)
        if (am_i_root()) then
          TOTALT=SUM(TIMING_glob(1:NTIMEACC)) ! over all processors
          DO M=1,NTIMEACC
            PERCENT(M) = 100d0*TIMING_glob(M)/(TOTALT+.00001)
          END DO
          TOTALT=SUM(TIMING(1:NTIMEACC)) ! on the root processor
          TOTALT=TOTALT/SECONDS_PER_MINUTE     ! seconds -> minutes
          DTIME = NDAY*TOTALT/(Itime-Itime0) ! minutes/day
          WRITE (6,'(/A,F7.2,A,/(8(A13,F5.1/))//)')
     *         '0TIME',DTIME,'(MINUTES) ',
     *         (TIMESTR(M),PERCENT(M),M=1,NTIMEACC)
        end if
        TIMING = 0
        START= NOW

      END IF  ! beginning of accumulation period

C**** CPU TIME FOR CALLING DIAGNOSTICS
      call stopTimer('Diagnostics')
      CALL TIMER (NOW,MDIAG)

C**** TEST FOR TERMINATION OF RUN
      IF (MOD(Itime,Nssw).eq.0) then
       IF (AM_I_ROOT()) then
        ! iflag values:
        !     0: Stopping as per user requestUser-requested stop
        !     1: Running
        !     2
        ! ------ stop_on ==> iflag =3
        if (iflag == 1) then  ! iflag==1 means it's running
          if (stop_on) then
              iflag = 3
          end if
        end if

        ! ------ User-requested stop ==> iflag=0
        if (iflag == 1) then
          open(3,file='flagGoStop',form='FORMATTED',status='OLD'
     &         ,err=210)
          read (3,'(A8)',end=210) str
          close (3)
 210            continue
          IF (str .ne. string_go) iflag=0
        endif

        ! ------ Timeout ==> iflag=2
        if (iflag == 1) then
          if (getWTime() >= wtime1) iflag = 2
        end if

        call broadcast(iflag)
       else
        call broadcast(iflag)
       end if
      endif
      select case (iflag)
        case (0)
           WRITE (6,'("0Flag to continue run has been turned off.")')
        case (2)
           WRITE (6,'("0Reached maximum wall clock time.")')
        case (3)
           WRITE (6,'("0Got signal 15.")')
      end select

      if (iflag .ne. 1) exit main_loop

      call stopTimer('Main Loop')
      END DO main_loop
C****
C**** END OF MAIN LOOP
C****

      call gettime(tloopend)
      if (AM_I_ROOT())
     *     write(*,*) "Time spent in the main loop in seconds:",
     *     tloopend-tloopbegin

C**** ALWAYS PRINT OUT RSF FILE WHEN EXITING
      CALL RFINAL (IRAND)
      call set_param( "IRAND", IRAND, 'o' )
      call io_rsf(rsf_file_name(KDISK),Itime,iowrite,ioerr)

      call finalize_atm

      if (AM_I_ROOT()) then
      WRITE (6,'(A,I1,45X,A4,I5,A5,I3,A4,I3,A,I8)')
     *  '0Restart file written on fort.',KDISK,'Year',year,
     *     aMON, date,', Hr',hour,'  Internal clock time:',ITIME
      end if

C**** RUN TERMINATED BECAUSE IT REACHED TAUE (OR SS6 WAS TURNED ON)

      call printSysTimers()

      IF (AM_I_ROOT())
     *   WRITE (6,'(/////4(1X,33("****")/)//,A,I8
     *             ///4(1X,33("****")/))')
     *  ' PROGRAM TERMINATED NORMALLY - Internal clock time:',ITIME

      IF (Itime.ge.ItimeE) then
         call reportProfile((itimee-itimei)*dtSRC)
         if (AM_I_ROOT()) call print_unused_param(6)
         CALL stop_model (
     &     'Terminated normally (reached maximum time)',13)
      END IF

      select case(iflag)
        case (0)
          CALL stop_model('Run stopped with sswE',12)  ! voluntary stop
        case(2)
          call stop_model('Reached maximum wall clock time.',12)
        case(3)
          call stop_model('Got signal 15',12)
      end select

#ifdef USE_MPP
      call fms_end( )
#endif

      contains

      subroutine initializeModelE()
      USE DOMAIN_DECOMP_1D, ONLY : init_app, am_i_root
      use Model_com, only: orbit, calendar, makeOrbit
      use Dictionary_mod
      USE MODEL_COM, only : master_yr
      use AbstractOrbit_mod, only: AbstractOrbit
      implicit none

      call initializeSysTimers()

#ifdef USE_MPP
      call fms_init( )
#endif
      call initializeConstants()
      call init_app()
      call initializeDefaultTimers()

      if (is_set_param("master_yr")) then
        call get_param( "master_yr", master_yr )
      else
        call stop_model('Please define master_yr in the rundeck.',255)
      endif

      allocate(orbit, source=makeOrbit())
      call orbit%setVerbose(am_I_root())

      allocate(calendar, source=orbit%makeCalendar())
      call calendar%setVerbose(am_I_root())

      if (am_i_root()) call calendar%print(2000)

      call alloc_drv_atm()
      call alloc_drv_ocean()

      end subroutine initializeModelE

      subroutine startNewDay()
      use model_com, only: modelEclock, calendar
      use CalendarMonth_mod
C**** INITIALIZE SOME DIAG. ARRAYS AT THE BEGINNING OF SPECIFIED DAYS
      logical :: newmonth
      integer :: month, day_of_month, year
      integer :: day_of_year
      type (CalendarMonth) :: cMonth

      year = modelEclock%getYear()
      month = modelEclock%getMonth()
      day_of_month = modelEclock%getDate()
      day_of_year = modelEclock%getDayOfYear()

        if (am_i_root()) then
          print '(A,I9,A,I0.4,A1,I0.2,A1,I0.2)',
     &       '---------- Main Loop, itime=',itime,
     &       ' day=',year,'-',month,'-',day_of_month
        end if


      cMonth = calendar%getCalendarMonth(month=month-1,year=year)
      newmonth = (day_of_year == 1+ cMonth%lastDayInMonth)
      call daily_DIAG(newmonth) ! atmosphere
      if (isBeginningAccumPeriod(modelEClock)) then
C**** THINGS THAT GET DONE AT THE BEGINNING OF EVERY ACC.PERIOD
        call reset_ADIAG(0)
        call reset_ODIAG(0)
#ifndef STANDALONE_OCEAN
        call reset_glaacc
#endif
      endif
      end subroutine startNewDay

!TODO fv, fv_fname, and fv_dfname are  not yet passed as arguments
!TODO exist except when building an FV version
      subroutine checkpointModelE
!@sum Every Ndisk Time Steps (DTsrc), starting with the first one,
!@+ write restart information alternately onto 2 disk files
      use MODEL_COM, only: rsf_file_name,kdisk,irand
      use MODEL_COM, only: itime
#ifdef USE_FVCORE
      USE FV_INTERFACE_MOD, only: Checkpoint,fvstate
#endif

      integer :: hour, date
      character(len=LEN_MONTH_ABBREVIATION) :: amon

      call modelEclock%get(hour=hour, date=date, amn=amon)

      CALL rfinal(IRAND)
      call set_param( "IRAND", IRAND, 'o' )
      call io_rsf(rsf_file_name(KDISK),Itime,iowrite,ioerr)
#if defined( USE_FVCORE )
      call checkpoint(fvstate, rsf_file_name(KDISK))
#endif
      if (AM_I_ROOT())
     *     WRITE (6,'(A,I1,45X,A4,I5,A5,I3,A4,I3,A,I8)')
     *     '0Restart file written on fort.',KDISK,'Year',
     *     year,aMon,date,', Hr',hour,'  Internal clock time:',ITIME
      kdisk=3-kdisk  ! Swap next fort.X.nc file

      end subroutine checkpointModelE

      subroutine initializeDefaultTimers()
      use TimerPackage_mod, only: initialize
      use TimerPackage_mod, only: addTimer

      call initialize()
      call addTimer('Main Loop')
      call addTimer(' Atm. Dynamics')
      call addTimer(' MELT_SI()')
      call addTimer(' CONDSE()')
      call addTimer(' RADIA()')
      call addTimer(' Surface')
      call addTimer('  Precip')
      call addTimer('  SURFACE()')
      call addTimer(' DYNSI()')
      call addTimer(' UNDERICE()')
      call addTimer(' GROUND_SI()')
      call addTimer(' GROUND_LI()')
      call addTimer(' GROUND_LK()')
      call addTimer(' RIVERF()')
      call addTimer(' OCEANS')
      call addTimer(' Daily')
      call addTimer(' Diagnostics')

      end subroutine initializeDefaultTimers

      subroutine reportProfile(elapsedTimeInSeconds)
      use TimerPackage_mod
      use DOMAIN_DECOMP_1D, only: AM_I_ROOT
      use TimeConstants_mod, only: INT_MINUTES_PER_DAY

#ifdef USE_MPI
      include 'mpif.h'
#endif
      real*8, intent(in) :: elapsedTimeInSeconds
      character(len=MAX_RECORD_LENGTH), pointer :: lines(:)
      integer :: i

      type (ProfileReport_type) :: report
      type (ReportColumn_type) :: column
      real*8 :: totalTime

      call finalize()
      report = newReport()
      call addColumn(report, newColumn(NAME_COLUMN, fieldWidth=20))

      column = newColumn(INCLUSIVE_TIME_COLUMN, fieldWidth=7)
      totalTime = getInclusiveTime(getTimer('main'))
      call setPrecision(column, 2)
      call setScale(column, 100/totalTime, ' %')
      call addColumn(report, column)

      column = newColumn(INCLUSIVE_TIME_COLUMN, fieldWidth=11)
      call setPrecision(column, 5)
      call setScale(column, INT_MINUTES_PER_DAY/elapsedTimeInSeconds,
     &     'min/day')
      call addColumn(report, column)

      column = newColumn(EXCLUSIVE_TIME_COLUMN, fieldWidth=11)
      call setPrecision(column, 5)
      call setScale(column, INT_MINUTES_PER_DAY/elapsedTimeInSeconds,
     &     'min/day')
      call addColumn(report, column)

      call addColumn(report, newColumn(TRIP_COUNTS_COLUMN,fieldWidth=6))

      column = newColumn(MAXIMUM_TIME_COLUMN, fieldWidth=10)
      call setPrecision(column, 6)
      call addColumn(report, column)

      call addColumn(report, newColumn(MAX_PROCESS_COLUMN,fieldWidth=4))

      column = newColumn(MINIMUM_TIME_COLUMN, fieldWidth=10)
      call setPrecision(column, 6)
      call addColumn(report, column)
      call addColumn(report, newColumn(MIN_PROCESS_COLUMN,fieldWidth=4))

      column = newColumn(AVERAGE_TIME_COLUMN, fieldWidth=10)
      call setPrecision(column, 6)
      call addColumn(report, column)

#ifdef USE_MPI
      lines => generateParallelReport(report, MPI_COMM_WORLD)
#else
      lines => generateReport(report)
#endif

      if (AM_I_ROOT()) then
         do i = 1, size(lines)
            write(*,'(a)') trim(lines(i))
         end do
      end if
      deallocate(lines)
      call delete(report)

      end subroutine reportProfile

      end subroutine GISS_modelE

      subroutine dailyUpdates
      use fluxes, only : atmocn
      implicit none

      call daily_CAL(.true.)    ! end_of_day
      call daily_OCEAN(.true.,atmocn)  ! end_of_day
      call daily_ATM(.true.)

      return
      end subroutine dailyUpdates

      subroutine sig_stop_model
      USE MODEL_COM, only : stop_on
      implicit none
      stop_on = .true.
      print *,"got signal 15"
      call sys_flush(6)
      end subroutine sig_stop_model

#ifdef USE_GDB_FOR_FPE_BACKTRACE
      subroutine sig_exception
#ifdef COMPILER_Intel8
      USE IFPORT
#endif
      implicit none
      character*80 :: str
      integer :: pid, retcode
#ifdef COMPILER_Intel8
      character*6 :: gdb = 'gdb-ia'
#else
      character*6 :: gdb = 'gdb   '
#endif

      write(6,*) "got signal 8"
      write(6,*) "floating point exception"
      call sys_flush(6)
      
      pid = getpid()
      write(str,'(a6,a,i16)')
     &     gdb, " -ex='set confirm off' -ex bt -ex quit -p ", pid
      !write(0,*) "command: ", str

      retcode = system(str)

      call stop_model("floating point exception",255)
      end subroutine sig_exception
#endif

      subroutine init_Model
!@sum This program reads most of parameters from the database (DB)
!@+   get_param( "A", X ) reads parameter A into variable X
!@+   if "A" is not in the database, it will generate an error
!@+   message and stop
!@+   sync_param( "B", Y ) reads parameter B into variable Y
!@+   if "B" is not in the database, then Y is unchanged and its
!@+   value is saved in the database as "B" (here sync = synchronize)
      USE MODEL_COM, only : NIPRNT,master_yr
     *     ,NMONAV,Ndisk,Nssw,KCOPY,KRSF,KOCEAN,IRAND,ItimeI
      USE DOMAIN_DECOMP_1D, only: AM_I_ROOT
      USE Dictionary_mod
      USE MDIAG_COM, only : make_timeaxis
      implicit none
      integer :: dummy_int

C**** Rundeck parameters:
      call sync_param( "NMONAV", NMONAV )
      call sync_param( "NIPRNT", NIPRNT )
      call sync_param( "Ndisk", Ndisk )
      call sync_param( "Nssw", Nssw )
      call sync_param( "KCOPY", KCOPY )
      call sync_param( "KRSF",   KRSF )
      call sync_param( "KOCEAN", KOCEAN )
      call sync_param( "IRAND", IRAND )
      if (is_set_param("master_yr")) then
        call get_param( "master_yr", master_yr )
      else
        call stop_model('Please define master_yr in the rundeck.',255)
      endif
      dummy_int = 0
      call sync_param("make_timeaxis",dummy_int)
      make_timeaxis = dummy_int==1
      RETURN
C****
      end subroutine init_Model


      SUBROUTINE INPUT (istart,ifile,coldRestart)
C****
C**** THIS SUBROUTINE SETS THE PARAMETERS IN THE C ARRAY, READS IN THE
C**** INITIAL CONDITIONS, AND CALCULATES THE DISTANCE PROJECTION ARRAYS
C****
      use TimeInterval_mod
      USE FILEMANAGER, only : openunit,closeunit
      USE TIMINGS, only : timing,ntimeacc
      USE Dictionary_mod
      USE MODEL_COM, only :
     *      xlabel,lrunid,nmonav,qcheck,irand
     *     ,nday,dtsrc,kdisk,jmon0,jyear0
     *     ,iyear1,itime,itimei,itimee
     *     ,idacc,modelEclock, modelEclockI
     *     ,aMONTH,aMON0
     *     ,ioread,ioread_acc,irerun,irsfic
     *     ,melse,Itime0,Jdate0
     *     ,Jhour0,rsf_file_name
     *     ,HOURI,DATEI,MONTHI,YEARI ,HOURE,DATEE,MONTHE,YEARE
     &     ,iwrite_sv,jwrite_sv,itwrite_sv,kdiag_sv
      USE RANDOM
      USE DOMAIN_DECOMP_1D, only : AM_I_ROOT
#ifndef STANDALONE_OCEAN
#if (defined TRACERS_ON) || (defined TRACERS_OCEAN)
      USE RESOLUTION, only : LM ! atm reference for init_tracer hack
#endif
#endif

      use TimeConstants_mod, only: INT_HOURS_PER_DAY
      use ModelClock_mod, only: ModelClock
      use Time_mod, only: Time, newTime
      use MODEL_COM, only: calendar, orbit
      use CalendarMonth_mod, only: LEN_MONTH_ABBREVIATION
      use TimeInterval_mod
      use Rational_mod, only: nint

      IMPLICIT NONE
!@var istart  postprocessing(-1)/start(1-8)/restart(>8)  option
      integer, intent(out) :: istart
      character(*), intent(in) :: ifile
      logical, intent(in) :: coldRestart
!@dbparam init_topog_related : set = 1 if IC and topography are incompatible
      integer :: init_topog_related = 0
!@dbparam do_IC_fixups : set = 1 if surface IC are to be checked/corrected
      integer :: do_IC_fixups = 0
!@var iu_AIC,iu_IFILE unit numbers for input files
      INTEGER iu_AIC,iu_IFILE
      INTEGER I,J,L,K,LID1,LID2,NOFF,ioerr
!@nlparam IHRI,TIMEE,IHOURE   end of model run
!@var  IHRI,IHOURE start and end of run in hours (from 1/1/IYEAR1 hr 0)
!@nlparam IRANDI  random number seed to perturb init.state (if>0)
      INTEGER :: IHRI=-1,TIMEE=-1,IHOURE=-1,IRANDI=0
      INTEGER IhrX
      INTEGER KDISK_restart   ! Name of fort.X.nc file from which we restarted
      LOGICAL :: is_coldstart
      CHARACTER NLREC*80,RLABEL*132

      INTEGER :: IWRITE=0,JWRITE=0,ITWRITE=23
      INTEGER, DIMENSION(13) :: KDIAG

      NAMELIST/INPUTZ/ ISTART,IRANDI
     *     ,IWRITE,JWRITE,ITWRITE,QCHECK,KDIAG
     *     ,IHOURE, TIMEE,HOURE,DATEE,MONTHE,YEARE,IYEAR1
C****    List of parameters that are disregarded at restarts
     *     ,        HOURI,DATEI,MONTHI,YEARI
      NAMELIST/INPUTZ_cold/ ISTART,IRANDI
     *     ,IWRITE,JWRITE,ITWRITE,QCHECK,KDIAG
     *     ,IHOURE, TIMEE,HOURE,DATEE,MONTHE,YEARE,IYEAR1
C****    List of parameters that are disregarded at restarts
     *     ,        HOURI,DATEI,MONTHI,YEARI
      integer istart_fixup
      character*132 :: bufs
      integer, parameter :: MAXLEN_RUNID = 32

      type (Time) :: modelETimeI, tmpTime, modelETime0, modelETimeE
      type (Time) :: modelETime
      integer :: hour, month, day, date, year

      character(len=80) :: tmpStr
      character(len=LEN_MONTH_ABBREVIATION) :: amon
      type (TimeInterval) :: dtSrcUsed
      type (TimeInterval) :: secsPerDay

C****
C**** Default setting for ISTART : restart from latest save-file (10)
C****
      ISTART=10

C**** All diagn. are enabled unless KDIAG is changed in the rundeck
      KDIAG(1:12)=0
      KDIAG(13)=9

C**** Set global default timing descriptions
C**** Other speciality descriptions can be added/used locally
      NTIMEACC = 0

C****
C**** Print Header and Label (2 lines) from rundeck
C****
      call openunit(trim(ifile),iu_IFILE,.false.,.true.)
      if (AM_I_ROOT()) WRITE (6,'(A,40X,A/)') '0','GISS CLIMATE MODEL'
      READ(iu_IFILE,'(A80)') XLABEL(1:80),NLREC
      NOFF=0
      IF (XLABEL(73:80).EQ.'        ') NOFF=8   ! for 72-column rundecks
      XLABEL(81-NOFF:132)=NLREC(1:52+NOFF)
      if (AM_I_ROOT()) WRITE (6,'(A,A/)') '0',XLABEL
      RLABEL = XLABEL !@var RLABEL rundeck-label

C****
C**** Print preprocessing options (if any are defined)
C****
      IF(AM_I_ROOT()) call print_and_check_PPopts

C****
C**** Read parameters from the rundeck to database and namelist
C****
      !call parse_params(iu_IFILE)
      ! skip "&&PARAMETERS" section
      do
        read( iu_IFILE, *, err=910, end=910 ) bufs
        if ( bufs == '&&END_PARAMETERS' ) exit
      enddo

      READ (iu_IFILE,NML=INPUTZ,ERR=900)
      if (coldRestart) READ (iu_IFILE,NML=INPUTZ_cold,ERR=900)
      call closeunit(iu_IFILE)

      IWRITE_sv = IWRITE
      JWRITE_sv = JWRITE
      ITWRITE_sv = ITWRITE
      KDIAG_sv = KDIAG

      if (istart.le.0) then
        call stop_model('pdE not supported',255)
      end if

C**** Get those parameters which are needed in this subroutine
      call get_param( "DTsrc", DTsrc )
      if(is_set_param("IRAND"))  call get_param( "IRAND", IRAND )

      if (istart.lt.9) then
C***********************************************************************
C****                                                               ****
C****                  INITIAL STARTS - ISTART: 2, 8                ****
C****                                                               ****
C****   Current settings: 2 - from observed data                    ****
C****                     8 - from current model M-file - no resets ****
C****                                                               ****
C***********************************************************************
C****
C**** Set quantities that are derived from the namelist parameters
C****
!@var NDAY=(1 day)/DTsrc : even integer; adjust DTsrc to be commensurate
        NDAY = 2*nint(calendar%getSecondsPerDay()/(DTsrc*2))
        dtSrcUsed = TimeInterval(calendar%getSecondsPerDay() / NDAY)
        DTsrc = real(dtSrcUsed)
        call set_param( "DTsrc", DTsrc, 'o')

C**** Get Start Time; at least YearI HAS to be specified in the rundeck
        IF (YearI.lt.0) then
          IF (AM_I_ROOT())
     *      WRITE(6,*) 'Please choose a proper start year yearI, not',
     *                  yearI
          call stop_model('INPUT: yearI not provided',255)
        END IF
        IF (Iyear1.lt.0) Iyear1 = yearI
        tmpTime = newTime(calendar)
        modelETime0 = newTime(calendar)

        call tmpTime%setByDate(yearI, monthI, dateI, hourI)
        call modelEtime0%setByDate(iyear1, month=1, date=1, hour=0)

        IhrI = nint((tmpTime -modelEtime0)/calendar%getSecondsPerHour())
        ITimeI = nint((tmpTime - modelEtime0)/ dtSrcUsed)
        Itime=ItimeI
        IF (IhrI.lt.0) then
          IF (AM_I_ROOT())
     *      WRITE(6,*) 'Improper start time OR Iyear1=',Iyear1,
     *      ' > yearI;',' yearI,monthI,dateI,hourI=',
     *      yearI,monthI,dateI,hourI
          call stop_model(
     &      'INPUT: Improper start date or base year Iyear1',255)
        END IF

        IF (ISTART.EQ.2) THEN
C****
C**** Cold Start: ISTART=2
C****
          XLABEL(1:80)='Observed atmospheric data from NMC tape'

C**** Set flag to initialise topography-related variables
          init_topog_related = 1

        ELSE IF (ISTART==8) THEN
C****
C**** Data from current type of RESTART FILE
C****
! no need to read SRHR,TRHR,FSF,TSFREZ,diag.arrays
          call io_rsf("AIC",IhrX,irsfic,ioerr)

          tmpTime = modelEtime0
          call tmpTime%add(calendar%getSecondsPerHour()*Ihrx)

          modelEtimeI = modelEtime0
          call modelEtimeI%add(calendar%getSecondsPerHour()*IhrI)

C**** Check consistency of starting time
          IF( ((modelEtimeI%getDayOfYear()/=tmpTime%getDayOfYear()) .or.
     &      (modelEtimeI%getHour() /= tmpTime%getHour())) ) then
            WRITE (6,*) ' Difference in hours between ',
     &       'Starting date and Data date:',
     &       modelEtimeI%getDayOfYear(),'d:',modelEtimeI%getHour(),'h ',
     &       tmpTime%getDayOfYear(),'d:',tmpTime%getHour(),'h '
            WRITE (6,*) 'Please change HOURI,DATEI,MONTHI'
             call stop_model('INPUT: start date inconsistent with data',
     &       255)
          ENDIF
        END IF

C**** Set flags to initialise some variables related to topography
        call sync_param( "init_topog_related", init_topog_related )

        IF (init_topog_related == 1) then
          do_IC_fixups = 1      ! new default, not necessarily final
        ENDIF

        IF (AM_I_ROOT())
     *    WRITE(6,'(A,i3,1x,a4,i5,a3,i3,3x,a,i2/" ",a)')
     *    '0Model started on',datei,aMONTH(monthi),yeari,' Hr',houri,
     *    'ISTART =',ISTART,XLABEL(1:80) ! report input file label
        XLABEL = RLABEL       ! switch to rundeck label

      else                    ! initial versus restart
C***********************************************************************
C****                                                               ****
C****                  RESTARTS: ISTART > 8                         ****
C****                                                               ****
C****   Current settings: 9 - from own model M-file                 ****
C****                    10 - from later of fort.1 or fort.2        ****
C****                    11 - from fort.1                           ****
C****                    12 - from fort.2                           ****
C****               13 & up - from earlier of fort.1 or fort.2      ****
C****                                                               ****
C***********************************************************************
C****
C**** DATA FROM end-of-month RESTART FILE     ISTART=9
C**** mainly used for REPEATS and delayed EXTENSIONS
        IF(ISTART==9) THEN      !  diag.arrays are not read in
          call io_rsf("AIC",Itime,irerun,ioerr)
          WRITE (6,'(A,I2,A,I11,A,A/)') '0Model restarted; ISTART=',
     *      ISTART,', TIME=',Itime,' ',XLABEL(1:80) ! sho input file label
          XLABEL = RLABEL       ! switch to rundeck label
          TIMING = 0
        ELSE
C****
C**** RESTART ON DATA SETS 1 OR 2, ISTART=10 or more
C****
C**** CHOOSE DATA SET TO RESTART ON
          IF(ISTART==11 .OR. ISTART==12) THEN
            ! ISTART=11: Use fort.1.nc
            ! ISTART=12: Use fort.2.nc
            KDISK=ISTART-10
          ELSEIF(ISTART==10 .OR. ISTART==13) THEN
            call find_later_rsf(kdisk)
            IF (ISTART==13) KDISK=3-KDISK ! Use earlier fort file, not later
          ENDIF
          if (istart == 14) then
              kdisk = 4    ! Start from fort.4.nc file
          end if
          call io_rsf(rsf_file_name(KDISK),Itime,ioread,ioerr)
          KDISK_restart = KDISK   ! Fort file we started from
          if (AM_I_ROOT())
     *      WRITE (6,'(A,I2,A,I11,A,A/)') '0RESTART DISK READ, UNIT',
     *      KDISK,', Time=',Itime,' ',XLABEL(1:80)

          ! Set up the first checkpoint file to write
          if (istart == 14) then
            ! Set from I file
            call get_param('KDISK', kdisk)
          else if (istart == 10) then
            ! Keep KDISK after reading from the later restart file, so that
            ! the same file is overwritten first; in case of trouble,
            ! the earlier restart file will still be available
          else if (istart.gt.10) then
            ! If user specified a fort.X.nc file,
            ! switch to the checkpoint file we did NOT start from
            KDISK=3-KDISK
          end if

        ENDIF

      endif                   ! initial versus restart

C***********************************************************************
C****                                                              *****
C****       INITIAL- AND RESTARTS: Final Initialization steps      *****
C****                                                              *****
C***********************************************************************

C**** initialize Lrunid (length of the identifying part of XLABEL)
C****
      lid1 = INDEX(XLABEL,'(') -1
      if (lid1.lt.1) lid1=MAXLEN_RUNID+1
      lid2 = INDEX(XLABEL,' ') -1
      if (lid2.lt.1) lid2=MAXLEN_RUNID+1
      LRUNID = min(lid1,lid2)
      IF (LRUNID.gt.MAXLEN_RUNID) call stop_model
     *     ('INPUT: Rundeck name too long. Shorten to 32 char or less'
     *     ,255)

C**** Update ItimeE only if YearE or IhourE is specified in the rundeck
C****
      modelETime0 = newTime(calendar)
      call modelEtime0%setByDate(iyear1, month=1, date=1, hour=0)

      ! dtSrcUsed is of type TimeInterval to guarantee exact arithmetic
      ! use real(...) to convert for convenience in other calculations.
      dtSrcUsed = TimeInterval(calendar%getSecondsPerDay() / NDAY)
      DTsrc = real(dtSrcUsed)

      modelETimeE = newTime(calendar)
      if (timee .lt. 0) then
        timee = houre*nday/INT_HOURS_PER_DAY
        call modelEtimeE%setByDate(yearE, monthE, dateE, houre)
      else
        call modelEtimeE%setByDate(yearE, monthE, dateE, 0)
        call modelEtimeE%add(DTsrcUsed * timee)
      end if
      ITimeE = nint((modelEtimeE - modelEtime0) / dtSrcUsed)


C**** Check consistency of DTsrc with NDAY
      if (is_set_param("DTsrc") .and.
     &     nint(calendar%getSecondsPerDay()/DTsrc) .ne. NDAY) then
        if (AM_I_ROOT()) then
          secsPerDay = calendar%getSecondsPerDay()
          write(6,*) 'DTsrc=',DTsrc,' has to stay at/be set to',
     &               real(secsPerDay / NDAY)
        end if
        call stop_model('INPUT: DTsrc inappropriately set',255)
      end if
      call set_param( "DTsrc", DTsrc, 'o' )   ! copy DTsrc into DB

C**** NMONAV has to be 1(default),2,3,4,6,12, i.e. a factor of 12
      if(is_set_param("NMONAV")) call get_param( "NMONAV", NMONAV )
      if (NMONAV.lt.1 .or. MOD(12,NMONAV).ne.0) then
        write (6,*) 'NMONAV has to be 1,2,3,4,6 or 12, not',NMONAV
        call stop_model('INPUT: nmonav inappropriately set',255)
      end if
      if (AM_I_ROOT())
     *     write (6,*) 'Diag. acc. period:',NMONAV,' month(s)'

C**** Updating Parameters: If any of them changed beyond this line
C**** use set_param(.., .., 'o') to update them in the database (DB)

C**** Get the rest of parameters from DB or put defaults to DB

      call init_Model

C**** Set date information

      modelETime0 = newTime(calendar)
      call modelEtime0%setByDate(iyear1, month=1, date=1, hour=0)
      call modelEtime0%add(dtSrcUsed * itime0)

      jyear0 = modelEtime0%getYear()
      jmon0 = modelEtime0%getMonth()
      jdate0 = modelEtime0%getDate()
      jhour0 = modelEtime0%getHour()
      amon0 = modelEtime0%getAbbreviation()

      modelETime = newTime(calendar)
      call modelEtime%setByDate(Iyear1, 1, 1, 0)
      modelEclockI = ModelClock(modelETime, dtSrcUsed, itimeI)
      call modelETime%add( dtSrcUsed * itime )

      modelEclock = ModelClock(modelEtime, dtSrcUsed, itime)

      ! These next two lines are not necessary - but act as
      ! a (poor) test that the alternate constructor for clocks
      ! is working.
      tmpStr = modelEclock%toString()
      modelEclock = ModelClock(tmpStr, calendar, dtSrcUsed)

      ! In the case of parameterized orbits, the year must now be set.
      ! Unfortunately, year is not available when orbit and calendar are
      ! established.
      year = modelEclock%getYear()
      call orbit%setYear(real(year,kind=8))

      CALL DAILY_cal(.false.)                  ! not end_of_day

#ifndef STANDALONE_OCEAN
#if (defined TRACERS_ON) || (defined TRACERS_OCEAN)
C**** Initialise tracer parameters and diagnostics
C**** MUST be before other init routines
      call laterInitTracerMetadata()
      call InitTracerMetadataAtmOcnCpler()
#endif
#endif

!!! hack: may be prevented if post-processing option is eliminated
      istart_fixup = istart
      if (istart==8 .and. do_IC_fixups==1) istart_fixup = 9

      is_coldstart = (istart<9 .and. init_topog_related == 1)
! long version:
!      is_coldstart = istart==2 .or. (istart==8 .and. init_topog_related == 1)

      call INPUT_ocean (istart,istart_fixup,
     &     do_IC_fixups,is_coldstart)

      call INPUT_atm(istart,istart_fixup,
     &     do_IC_fixups,is_coldstart,
     &     KDISK_restart, IRANDI)

      if (istart.le.9) then
        call reset_adiag(0)
        call reset_odiag(0)
#ifndef STANDALONE_OCEAN
        call reset_glaacc
#endif
      endif

      CALL SET_TIMER("       OTHER",MELSE)

      if (AM_I_ROOT()) then
         WRITE (6,INPUTZ)
         call print_param( 6 )
         WRITE (6,'(A7,12I6)') "IDACC=",(IDACC(I),I=1,12)
      end if

#ifdef DEFER_ACC_READ
      if(istart.ge.10) then
        ! Reading of diagnostic accmulation arrays deferred until
        ! full metadata is known.
        ! Todo: defer reading of most other arrays as well.
        call io_rsf(rsf_file_name(kdisk_restart),itime,ioread_acc,ioerr)
      endif
#endif

#ifdef CACHED_SUBDD
      ! Initialize subdaily diagnostics
      call parse_subdd
      call reset_cached_subdd
      if(istart.ge.10) then
        call read_subdd_rsf(trim(rsf_file_name(kdisk_restart))//'.nc')
      endif
#endif

C****
      RETURN
C****
C**** TERMINATE BECAUSE OF IMPROPER PICK-UP
C****
 900   write (6,*) 'Error in NAMELIST parameters'
      call stop_model('Error in NAMELIST parameters',255)
 910   write (6,*) 'Error readin I-file'
      call stop_model('Error reading I-file',255)


      END SUBROUTINE INPUT

      subroutine print_and_check_PPopts
!@sum prints preprocessor options in english and checks some
!@+  interdependencies. (moved from subroutine INPUT).
!@+  Called by root thread only.
      use runtimecontrols_mod, only: tracers_gasexch_ocean,
     &   tracers_oceanbiology, tracers_gasexch_ocean_cfc,
     &   tracers_gasexch_ocean_co2
      implicit none

#if (defined TRACERS_ON) || (defined TRACERS_OCEAN)
      write(6,*) 'This program includes tracer code'
#endif
#ifdef TRACERS_WATER
      write(6,*) '...and water tracer code'
#ifndef TRACERS_ON
      call stop_model(
     &' Water tracers need TRACERS_ON as well as TRACERS_WATER',255)
#endif
#endif
#ifdef TRACERS_OCEAN
      write(6,*) '...and ocean tracer code'
#endif
#ifdef TRACERS_SPECIAL_O18
      write(6,*) '...and water isotope code'
#ifdef TRACERS_WISO_O17
      write(6,*) '...with the H217O isotope tracer'
#endif
#ifndef TRACERS_WATER
      call stop_model('Water isotope tracers need TRACERS_WATER '//
     *'as well as TRACERS_SPECIAL_O18',255)
#endif
#endif
#ifdef TRACERS_SPECIAL_Lerner
      write(6,*) '...and Jean/David tracers and chemistry'
#endif
      if (tracers_gasexch_ocean) then
        write(6,*) '          '
        write(6,*) '...and Natassa Romanou air-sea GAS EXCHANGE'
        if (tracers_oceanbiology) then
          write(6,*) '          '
          write(6,*)'...and Natassa Romanou/Watson Gregg ocean biology '
        endif
        if (tracers_gasexch_ocean_cfc)
     &           write(6,*) '****CFC flux across air/sea interface****'
        if (tracers_gasexch_ocean_co2)
     &           write(6,*) '****CO2 flux across air/sea interface****'
      endif
#ifdef TRACERS_SPECIAL_Shindell
      write(6,*) '...and Drew Shindell tracers and chemistry'
#endif
#ifdef TRACERS_TERP
#ifdef TRACERS_SPECIAL_Shindell
      write(6,*) '...and Terpenes tracer'
#else
      call stop_model('Terpenes tracer needs tropo chemistry',255)
#endif
#endif  /* TRACERS_TERP */
#ifdef CALCULATE_FLAMMABILITY
      write(6,*) '...and calculating sfc veg flammability'
#endif
#if(defined CALCULATE_LIGHTNING)||(defined TRACERS_SPECIAL_Shindell)
      write(6,*) '...and calculating lightning flash rate'
#endif
#ifdef DYNAMIC_BIOMASS_BURNING
      write(6,*) '...and dynamic biomass burning srcs by flammability'
#endif
#ifdef TRACERS_AEROSOLS_Koch
      write(6,*) '...and Dorothy Koch aerosols'
#ifdef TRACERS_AEROSOLS_VBS
      write(6,*) '...and VBS organics'
#endif
#endif
#ifdef TRACERS_AEROSOLS_SEASALT
      write(6,*) '...and sea salt aerosols'
#endif
#ifdef TRACERS_AMP
      write(6,*) '...and aerosol microphysics'
#endif
#ifdef TRACERS_TOMAS
      write(6,*) '...and TOMAS aerosol microphysics'
#endif
#ifdef TRACERS_AEROSOLS_SOA
#ifdef TRACERS_SPECIAL_Shindell
      write(6,*) '...and secondary organic aerosols'
#else
      call stop_model('SOA version needs tropo chemistry',255)
#endif
#endif  /* TRACERS_AEROSOLS_SOA */
#ifdef SOA_DIAGS
#ifdef TRACERS_AEROSOLS_SOA
      write(6,*) '...and additional SOA diagnostics'
#else
      call stop_model('SOA_DIAGS needs TRACERS_AEROSOLS_SOA',255)
#endif  /* TRACERS_AEROSOLS_SOA */
#endif  /* SOA_DIAGS */
#ifdef TRACERS_AEROSOLS_OCEAN
      write(6,*) '...and oceanic organic aerosol sources'
#endif  /* TRACERS_AEROSOLS_OCEAN */
#ifdef TRACERS_DRYDEP
      write(6,*) '...and tracer dry deposition'
#endif
#ifdef EDGAR_HYDE_SOURCES
      write(6,*) '...and EDGAR HYDE sources instead of GISS'
#endif
#ifdef SHINDELL_STRAT_EXTRA
      write(6,*) '...and Drew Shindell extra strat tracers'
#endif
#ifdef INTERACTIVE_WETLANDS_CH4
      write(6,*) '...and interactive CH4 wetlands emissions'
#endif
#ifdef NUDGE_ON
      write(6,*) '...and nudging of meteorology'
#ifdef MERRA_NUDGING
      write(6,*) '   WITH MERRA WINDS.'
#endif
#endif
#ifdef ACCMIP_LIKE_DIAGS
      write(6,*) '...and ACCMIP set of diagnostics'
#ifndef SHINDELL_STRAT_EXTRA
      call stop_model
     & ('SHINDELL_STRAT_EXTRA should be on for ACCMIP_LIKE_DIAGS',255)
#endif
#endif /* ACCMIP_LIKE_DIAGS */

      return
      end subroutine print_and_check_PPopts

