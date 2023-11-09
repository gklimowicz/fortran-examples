*
* $Id: mnstin.F,v 1.1.1.1 1996/03/07 14:31:32 mclareni Exp $
*
* $Log: mnstin.F,v $
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
      SUBROUTINE MNSTIN(CRDBUF,IERR)
*
* $Id: d506dp.inc,v 1.1.1.1 1996/03/07 14:31:32 mclareni Exp $
*
* $Log: d506dp.inc,v $
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506dp.inc
*
C ************ DOUBLE PRECISION VERSION *************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
CC Called from MNREAD.
CC Implements the SET INPUT command to change input units.
CC If command is: 'SET INPUT'   'SET INPUT 0'   or  '*EOF',
CC                 or 'SET INPUT , ,  ',
CC                reverts to previous input unit number,if any.
CC
CC      If it is: 'SET INPUT n'  or  'SET INPUT n filename',
CC                changes to new input file, added to stack
CC
CC      IERR = 0: reading terminated normally
CC             2: end-of-data on primary input file
CC             3: unrecoverable read error
CC             4: unable to process request
CC
*
* $Id: d506cm.inc,v 1.1.1.1 1996/03/07 14:31:32 mclareni Exp $
*
* $Log: d506cm.inc,v $
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506cm.inc
*
      PARAMETER (MNE=100 , MNI=50)
      PARAMETER (MNIHL=MNI*(MNI+1)/2)
      CHARACTER*10 CPNAM
      COMMON
     1/MN7NAM/ CPNAM(MNE)
     2/MN7EXT/ U(MNE)     ,ALIM(MNE)  ,BLIM(MNE)
     3/MN7ERR/ ERP(MNI)   ,ERN(MNI)   ,WERR(MNI)  ,GLOBCC(MNI)
     4/MN7INX/ NVARL(MNE) ,NIOFEX(MNE),NEXOFI(MNI)
     5/MN7INT/ X(MNI)     ,XT(MNI)    ,DIRIN(MNI)
     6/MN7FX2/ XS(MNI)    ,XTS(MNI)   ,DIRINS(MNI)
     7/MN7DER/ GRD(MNI)   ,G2(MNI)    ,GSTEP(MNI) ,GIN(MNE) ,DGRD(MNI)
     8/MN7FX3/ GRDS(MNI)  ,G2S(MNI)   ,GSTEPS(MNI)
     9/MN7FX1/ IPFIX(MNI) ,NPFIX
     A/MN7VAR/ VHMAT(MNIHL)
     B/MN7VAT/ VTHMAT(MNIHL)
     C/MN7SIM/ P(MNI,MNI+1),PSTAR(MNI),PSTST(MNI) ,PBAR(MNI),PRHO(MNI)
C
      PARAMETER (MAXDBG=10, MAXSTK=10, MAXCWD=20, MAXP=30, MAXCPT=101)
      PARAMETER (ZERO=0.0,  ONE=1.0,   HALF=0.5)
      COMMON
     D/MN7NPR/ MAXINT ,NPAR   ,MAXEXT ,NU
     E/MN7IOU/ ISYSRD ,ISYSWR ,ISYSSA ,NPAGWD ,NPAGLN ,NEWPAG
     E/MN7IO2/ ISTKRD(MAXSTK) ,NSTKRD ,ISTKWR(MAXSTK) ,NSTKWR
     F/MN7TIT/ CFROM  ,CSTATU ,CTITL  ,CWORD  ,CUNDEF ,CVRSN ,COVMES
     G/MN7FLG/ ISW(7) ,IDBG(0:MAXDBG) ,NBLOCK ,ICOMND
     H/MN7MIN/ AMIN   ,UP     ,EDM    ,FVAL3  ,EPSI   ,APSI  ,DCOVAR
     I/MN7CNV/ NFCN   ,NFCNMX ,NFCNLC ,NFCNFR ,ITAUR,ISTRAT,NWRMES(2)
     J/MN7ARG/ WORD7(MAXP)
     K/MN7LOG/ LWARN  ,LREPOR ,LIMSET ,LNOLIM ,LNEWMN ,LPHEAD
     L/MN7CNS/ EPSMAC ,EPSMA2 ,VLIMLO ,VLIMHI ,UNDEFI ,BIGEDM,UPDFLT
     M/MN7RPT/ XPT(MAXCPT)    ,YPT(MAXCPT)
     N/MN7CPT/ CHPT(MAXCPT)
     o/MN7XCR/ XMIDCR ,YMIDCR ,XDIRCR ,YDIRCR ,KE1CR  ,KE2CR
      CHARACTER CTITL*50, CWORD*(MAXCWD), CUNDEF*10, CFROM*8,
     +          CVRSN*6,  COVMES(0:3)*22, CSTATU*10, CHPT*1
      LOGICAL   LWARN, LREPOR, LIMSET, LNOLIM, LNEWMN, LPHEAD
      CHARACTER CRDBUF*(*),CUNIT*10,CFNAME*64,CGNAME*64,CANSWR*1
      CHARACTER CMODE*16
      LOGICAL LOPEN,LREWIN,NONAME,LNAME,MNUNPT
      NONAME = .TRUE.
      IERR = 0
      IF (INDEX(CRDBUF,'*EOF') .EQ. 1) GO TO 190
      IF (INDEX(CRDBUF,'*eof') .EQ. 1) GO TO 190
      LEND = LEN(CRDBUF)
C                               look for end of SET INPUT command
        DO 20 IC= 8,LEND
        IF (CRDBUF(IC:IC) .EQ. ' ') GO TO 25
        IF (CRDBUF(IC:IC) .EQ. ',') GO TO 53
   20   CONTINUE
      GO TO 200
   25 CONTINUE
C         look for end of separator between command and first argument
      ICOL = IC+1
         DO 50 IC= ICOL,LEND
         IF (CRDBUF(IC:IC) .EQ. ' ') GO TO 50
         IF (CRDBUF(IC:IC) .EQ. ',') GO TO 53
         GO TO 55
   50 CONTINUE
      GO TO 200
   53 IC = IC + 1
   55 IC1 = IC
C                      see if "REWIND" was requested in command
      LREWIN = .FALSE.
      IF (INDEX(CRDBUF(1:IC1),'REW') .GT. 5)  LREWIN=.TRUE.
      IF (INDEX(CRDBUF(1:IC1),'rew') .GT. 5)  LREWIN=.TRUE.
C                      first argument begins in or after col IC1
      DO 75 IC= IC1,LEND
      IF (CRDBUF(IC:IC) .EQ. ' ') GO TO 75
      IF (CRDBUF(IC:IC) .EQ. ',') GO TO 200
      GO TO 80
   75 CONTINUE
      GO TO 200
   80 IC1 = IC
C                        first argument really begins in col IC1
      DO 100 IC= IC1+1,LEND
      IF (CRDBUF(IC:IC) .EQ. ' ') GO TO 108
      IF (CRDBUF(IC:IC) .EQ. ',') GO TO 108
  100 CONTINUE
      IC = LEND + 1
  108 IC2 = IC-1
C                            end of first argument is in col IC2
  110 CONTINUE
      CUNIT = CRDBUF(IC1:IC2)
      WRITE (ISYSWR,'(A,A)') ' UNIT NO. :',CUNIT
      READ (CUNIT,'(BN,F10.0)',ERR=500) FUNIT
      IUNIT = FUNIT
      IF (IUNIT .EQ. 0)  GO TO 200
C                             skip blanks and commas, find file name
      DO 120 IC= IC2+1,LEND
      IF (CRDBUF(IC:IC) .EQ. ' ') GO TO 120
      IF (CRDBUF(IC:IC) .EQ. ',') GO TO 120
      GO TO 130
  120 CONTINUE
      GO TO 131
  130 CONTINUE
      CFNAME = CRDBUF(IC:LEND)
      NONAME = .FALSE.
      WRITE (ISYSWR, '(A,A)') ' FILE NAME IS:',CFNAME
C              ask if file exists, if not ask for name and open it
  131 CONTINUE
      INQUIRE(UNIT=IUNIT,OPENED=LOPEN,NAMED=LNAME,NAME=CGNAME)
      IF (LOPEN) THEN
         IF (NONAME) THEN
             GO TO 136
         ELSE
             IF (.NOT.LNAME) CGNAME='unknown'
             WRITE (ISYSWR,132) IUNIT,CGNAME,CFNAME
  132        FORMAT (' UNIT',I3,' ALREADY OPENED WITH NAME:',A/
     +                  '                 NEW NAME IGNORED:',A)
         ENDIF
      ELSE
C                new file, open it
         WRITE (ISYSWR,135) IUNIT
  135    FORMAT (' UNIT',I3,' IS NOT OPENED.')
         IF (NONAME) THEN
            WRITE (ISYSWR,'(A)') ' NO FILE NAME GIVEN IN COMMAND.'
            IF (ISW(6) .LT. 1)  GO TO 800
            WRITE (ISYSWR,'(A)') ' PLEASE GIVE FILE NAME:'
            READ (ISYSRD,'(A)') CFNAME
         ENDIF
         OPEN (UNIT=IUNIT,FILE=CFNAME,STATUS='OLD',ERR=600)
         WRITE (ISYSWR,'(A)') ' FILE OPENED SUCCESSFULLY.'
      ENDIF
C                                     . .   file is correctly opened
  136 IF (LREWIN) GO TO 150
      IF (ISW(6) .LT. 1)  GO TO 300
      WRITE (ISYSWR,137)  IUNIT
  137 FORMAT (' SHOULD UNIT',I3,' BE REWOUND?' )
      READ  (ISYSRD,'(A)')  CANSWR
      IF (CANSWR.NE.'Y' .AND. CANSWR.NE.'y') GO TO 300
  150 REWIND IUNIT
      GO TO 300
C                      *EOF
  190 CONTINUE
      IF (NSTKRD .EQ. 0)  THEN
         IERR = 2
         GO TO 900
         ENDIF
C                      revert to previous input file
  200 CONTINUE
      IF (NSTKRD .EQ. 0)  THEN
          WRITE (ISYSWR, '(A,A)') ' COMMAND IGNORED:',CRDBUF
          WRITE (ISYSWR, '(A)') ' ALREADY READING FROM PRIMARY INPUT'
      ELSE
        ISYSRD = ISTKRD(NSTKRD)
        NSTKRD = NSTKRD - 1
        IF (NSTKRD .EQ. 0)  ISW(6) = IABS(ISW(6))
        IF (ISW(5) .GE. 0)  THEN
          INQUIRE(UNIT=ISYSRD,NAMED=LNAME,NAME=CFNAME)
          CMODE = 'BATCH MODE      '
          IF (ISW(6) .EQ. 1)  CMODE = 'INTERACTIVE MODE'
          IF (.NOT.LNAME) CFNAME='unknown'
          IF (MNUNPT(CFNAME))  CFNAME='unprintable'
          WRITE (ISYSWR,290) CMODE,ISYSRD,CFNAME
  290     FORMAT (' INPUT WILL NOW BE READ IN ',A,' FROM UNIT NO.',I3/
     +    ' FILENAME: ',A)
        ENDIF
      ENDIF
      GO TO 900
C                      switch to new input file, add to stack
  300 CONTINUE
      IF (NSTKRD .GE. MAXSTK)  THEN
          WRITE (ISYSWR, '(A)') ' INPUT FILE STACK SIZE EXCEEDED.'
          GO TO 800
          ENDIF
      NSTKRD = NSTKRD + 1
      ISTKRD(NSTKRD) = ISYSRD
      ISYSRD = IUNIT
C                   ISW(6) = 0 for batch, =1 for interactive, and
C                      =-1 for originally interactive temporarily batch
      IF (ISW(6) .EQ. 1)  ISW(6) = -1
      GO TO 900
C                      format error
  500 CONTINUE
      WRITE (ISYSWR,'(A,A)') ' CANNOT READ FOLLOWING AS INTEGER:',CUNIT
      GO TO 800
  600 CONTINUE
      WRITE (ISYSWR, 601) CFNAME
  601 FORMAT (' SYSTEM IS UNABLE TO OPEN FILE:',A)
C                      serious error
  800 CONTINUE
      IERR = 3
  900 CONTINUE
      RETURN
      END
