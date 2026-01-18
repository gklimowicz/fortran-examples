      program qc
!@sum  qc query c-array and other parameters
!@auth G. Schmidt/R. Ruedy
      USE CONSTANT
      USE MODEL_COM
      USE TIMINGS
      USE Dictionary_mod
      USE FILEMANAGER
      use TimeConstants_mod, only: HOURS_PER_DAY, DAYS_PER_YEAR
      IMPLICIT NONE
      CHARACTER*80 FILEIN
      INTEGER N,NARGS,K,iargc,KFILE,I,days_togo,itm,iu_RSF
      INTEGER :: ioerr=0, KSTART=1, ItimeMax=0
      REAL*8 TOT,yrs_togo,FAC,FACT,xfac,hour
      LOGICAL :: QCALL = .FALSE., QCMIN=.FALSE., QCRESTART=.FALSE.
!@var QCRESTART if TRUE compute max Itime and do printout for "runpm"

      NARGS = IARGC()
      IF(NARGS.LE.0)  GO TO 800
C**** check for arguments
 10   CALL GETARG(KSTART,FILEIN)
      IF (FILEIN(1:1).eq."-") THEN
        SELECT CASE (FILEIN(2:2))
        CASE ("a","A")
          QCALL=.TRUE.
        CASE ("t","T")
          QCMIN=.TRUE.
          xfac = 1.
          if (filein(3:3).ne.' ') read(filein(3:80),*) xfac
        CASE ("r","R")
          QCRESTART=.TRUE.
        END SELECT
        KSTART=KSTART+1
        GOTO 10
      END IF

      DO K=KSTART,NARGS
      if (qcall .and. k.gt.kstart) write (6,*)
      CALL GETARG (K,FILEIN)
      !OPEN (10,FILE=FILEIN,FORM='UNFORMATTED',STATUS='OLD',err=850)
      ioerr=0
      call openunit(FILEIN,iu_RSF,.true.,.true.)
      call io_label(iu_RSF,Itime,itm,ioread,ioerr)
      call closeunit(iu_RSF)
      if (ioerr.eq.1) go to 860
      !CLOSE (10)

      if ( QCRESTART ) then ! find max Itime; skip the rest of the loop
        ItimeMax = max ( ItimeMax, Itime )
        cycle
      endif

      call getdte(Itime,Nday,Iyear1,Jyear,Jmon,Jday,Jdate,Jhour,amon)
      hour=mod(itime,nday)*HOURS_PER_DAY/nday

      WRITE (6,900) ITIME,JMON,JDATE,JYEAR,HOUR,XLABEL(1:50)
      TOT=0
      DO N=1,NTIMEACC
        TOT = TOT + TIMING(N)
      END DO
      IF (Itime-Itime0.gt.0) THEN
        FACT = NDAY/(60.*(Itime-Itime0))
      ELSE
        FACT = 0.
      END IF
      IF (QCMIN) THEN ! output in minutes
        FAC = FACT*xfac
      ELSE            ! output in percentages
        if(TOT.gt.0) FAC = 100./TOT
      END IF
      IF (TOT.gt.0) THEN
        WRITE (6,906) FACT*TOT,(TIMESTR(N),FAC*TIMING(N),N=1,3)
        WRITE (6,907) (TIMESTR(N),FAC*TIMING(N),N=4,NTIMEACC)
      END IF
      IF (QCALL) THEN
        call print_param(6)
c       write (6,*) "IDACC = ",(IDACC(I),I=1,12)
        call getdte(ItimeI,Nday,Iyear1,Jyear,Jmon,Jday,Jdate,Jhour,amon)
        hour=mod(itime,nday)*HOURS_PER_DAY/nday
        WRITE (6,900) ITIMEI,JMON,JDATE,JYEAR,HOUR,' = start of run'
        call getdte(ItimeE,Nday,Iyear1,Jyear,Jmon,Jday,Jdate,Jhour,amon)
        hour=mod(itime,nday)*HOURS_PER_DAY/nday
        WRITE (6,900) ITIMEE,JMON,JDATE,JYEAR,HOUR,' =   end of run'
        if(itimee.ge.itime) then
          days_togo = (Itimee-itime+nday-1)/nday
          yrs_togo  = (Itimee-itime)/(DAYS_PER_YEAR*nday)
          write(XLABEL(29:50),'(I10,a12)') days_togo,'  days to go'
          if (days_togo.eq.1) XLABEL(44:44) = ' '
         call getdte(Itime,Nday,Iyear1,Jyear,Jmon,Jday,Jdate,Jhour,amon)
         hour=mod(itime,nday)*HOURS_PER_DAY/nday
          if (yrs_togo.ge.2.)
     *    write(XLABEL(29:50),'(f10.1,a12)') yrs_togo,' years to go'
        WRITE (6,900) ITIME,JMON,JDATE,JYEAR,HOUR,XLABEL(1:50)
      end if
      END IF
      END DO
      if ( QCRESTART ) then ! print data needed for restart
        call getdte(
     &       ItimeMax,Nday,Iyear1,Jyear,Jmon,Jday,Jdate,Jhour,amon)
        write(6,"('QCRESTART_DATA: ',I10,1X,I2,'-',I2.2,'-',I4.4)")
c        write(6,*)
     &       ItimeMax*nint(HOURS_PER_DAY)/Nday, Jmon, Jdate, Jyear
      endif
      Stop
  800 continue
      write(6,*) " Usage: qc [-a] [-t] [-r] FILE1 [FILE2...]"
      write(6,*) " Where FILE1,2.. are modelE files ",
     *     "(rsf/acc/fort.[12])"
      write(6,*) "  -a  output all header information, otherwise"
      write(6,*) "      only timing info is printed"
      write(6,*) "  -tx output in minutes/x for each sub-section,"
      write(6,*) "      otherwise percentage time is printed"
      write(6,*) "  -r  find maximal Itime and do printout in a format,"
      write(6,*) "      convenient for restart script (runpm)"
      STOP
  850 continue
      write(6,*) "Cannot open file ",FILEIN
      STOP
  860 continue
      write(6,*) "Error reading file ", FILEIN
      STOP
 900  FORMAT (I10,1X,I2,'/',I2.2,'/',I4.4,' hr',f4.1,1X,A)
 906  FORMAT (' TIME',F7.2,' (MINUTES)',3(1X,A12,F5.1))
 907  FORMAT (10(22X,3(1X,A12,F5.1) / ))
      end

