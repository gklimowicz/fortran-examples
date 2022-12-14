      PROGRAM mkdeep
!@sum  mkdeep creates the mixed layer temperature climatology needed of
!@+    deep ocean diffusion runs
!@auth G. Russell/L. Nazarenko/G. Schmidt
      IMPLICIT NONE
      CHARACTER TITLE*80,FNAME*60,RUNID*20
      CHARACTER*3,dimension(12) :: amonth=(/'JAN','FEB','MAR',
     *   'APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC'/)
      INTEGER IM,JM,M,MONACC(12),K,MTG3,IYR,IYRE,IYRI,NFILES,IARGC
     *     ,lf,itime
!!    REAL*8 TG3M(IM,JM,12),TG3(IM,JM)
      real*8, allocatable :: TG3M(:,:,:),TG3(:,:)

      write(*,*) 'enter grid dimensions im,jm:'
      read(*,*) im,jm
      allocate (TG3M(im,jm,12),TG3(im,jm))

C**** ZERO OUT ACCUMULATING ARRAY
      TG3M = 0.

      NFILES=IARGC()
      IF (NFILES.le.0) THEN
        PRINT*,"mkdeep: make climatology to initialise deep ocean runs"
        PRINT*,"Usage: mkdeep.exe set_of_oda_files"
        PRINT*,"  after linking the oda_files to the current directory"
        STOP
      END IF

      IYRE=0
      IYRI=9999
      monacc=0

C**** loop over input files
      DO K=1,NFILES
        CALL GETARG(K,FNAME)
C****   get year and month from file name FNAME = monYYYY.odaRUNID
        READ(FNAME(4:7),'(I4)') IYR
        IYRE=MAX(IYRE,IYR)
        IYRI=MIN(IYRI,IYR)
        MTG3=0
        DO M=1,12
          IF (FNAME(1:3).eq.AMONTH(M)(1:3)) THEN
            monacc(M)=monacc(M)+1
            MTG3=M
          END IF
        END DO
        if (MTG3.eq.0) then
          PRINT*,trim(FNAME)," non-standard, month=",FNAME(1:3)
          stop
        end if

        open(10,file=FNAME,form='unformatted')
        read(10,err=910) itime,tg3
        close(10)

C**** ACCUMULATE TG3M SEPARATELY FOR EACH MONTH
        TG3M(:,:,MTG3)=TG3M(:,:,MTG3)+TG3(:,:)

      END DO
      RUNID=FNAME(12:len_trim(FNAME))
C**** DIVIDE EACH MONTH BY NUMBER OF ACCUMULATIONS
      DO M=1,12
        TG3M(:,:,M)=TG3M(:,:,M)/monacc(M)
      END DO

      TITLE="mean daily sums of TG3 for each month, TG3M"
      FNAME="TG3M."//trim(RUNID)
      lf=len_trim(FNAME)
      WRITE(FNAME(lf+1:lf+10),'(a1,I4.0,a1,I4.0)') '.',IYRI,'-',IYRE

      open(10,file=FNAME,form='unformatted')
      WRITE(10) TITLE,TG3M
      close(10)

      STOP
  910 WRITE(6,*) 'READ ERROR ON FILE ',FNAME
      STOP
      END
