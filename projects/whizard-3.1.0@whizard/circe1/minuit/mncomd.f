*
* $Id: mncomd.F,v 1.1.1.1 1996/03/07 14:31:29 mclareni Exp $
*
* $Log: mncomd.F,v $
* Revision 1.1.1.1  1996/03/07 14:31:29  mclareni
* Minuit
*
*
      SUBROUTINE MNCOMD(FCN,CRDBIN,ICONDN,FUTIL)
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
CC        Called by user.  'Reads' a command string and executes.
CC     Equivalent to MNEXCM except that the command is given as a
CC          character string.
CC
CC     ICONDN = 0: command executed normally
CC              1: command is blank, ignored
CC              2: command line unreadable, ignored
CC              3: unknown command, ignored
CC              4: abnormal termination (e.g., MIGRAD not converged)
CC              5: command is a request to read PARAMETER definitions
CC              6: 'SET INPUT' command
CC              7: 'SET TITLE' command
CC              8: 'SET COVAR' command
CC              9: reserved
CC             10: END command
CC             11: EXIT or STOP command
CC             12: RETURN command
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
      DIMENSION PLIST(MAXP)
      CHARACTER COMAND*(MAXCWD)
      CHARACTER CLOWER*26, CUPPER*26
      LOGICAL LEADER
C
      EXTERNAL FCN,FUTIL
      CHARACTER*(*) CRDBIN
      CHARACTER*100 CRDBUF
      DATA CLOWER/'abcdefghijklmnopqrstuvwxyz'/
      DATA CUPPER/'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
C
      LENBUF = LEN(CRDBIN)
      CRDBUF = CRDBIN
      ICONDN = 0
C     record not case-sensitive, get upper case, strip leading blanks
      LEADER = .TRUE.
      IPOS = 1
         DO 110 I= 1, MIN(MAXCWD,LENBUF)
         IF (CRDBUF(I:I) .EQ. '''') GO TO 111
         IF (CRDBUF(I:I) .EQ. ' ')  THEN
           IF (LEADER) IPOS = IPOS + 1
           GO TO 110
         ENDIF
         LEADER = .FALSE.
           DO 108 IC= 1, 26
           IF (CRDBUF(I:I) .EQ. CLOWER(IC:IC)) CRDBUF(I:I)=CUPPER(IC:IC)
  108      CONTINUE
  110    CONTINUE
  111 CONTINUE
C                     blank or null command
      IF (IPOS .GT. LENBUF)  THEN
         WRITE (ISYSWR,'(A)') ' BLANK COMMAND IGNORED.'
         ICONDN = 1
         GO TO 900
      ENDIF
C                                           . .   preemptive commands
C               if command is 'PARAMETER'
      IF (CRDBUF(IPOS:IPOS+2) .EQ. 'PAR')    THEN
         ICONDN = 5
         LPHEAD = .TRUE.
         GO TO 900
         ENDIF
C               if command is 'SET INPUT'
      IF (CRDBUF(IPOS:IPOS+6) .EQ. 'SET INP')  THEN
         ICONDN = 6
         LPHEAD = .TRUE.
         GO TO 900
         ENDIF
C              if command is 'SET TITLE'
      IF (CRDBUF(IPOS:IPOS+6) .EQ. 'SET TIT')  THEN
         ICONDN = 7
         LPHEAD = .TRUE.
         GO TO 900
         ENDIF
C               if command is 'SET COVARIANCE'
      IF (CRDBUF(IPOS:IPOS+6) .EQ. 'SET COV')   THEN
         ICONDN = 8
         LPHEAD = .TRUE.
         GO TO 900
         ENDIF
C               crack the command . . . . . . . . . . . . . . . .
      CALL MNCRCK(CRDBUF(IPOS:LENBUF),MAXCWD,COMAND,LNC,
     +                            MAXP,  PLIST, LLIST, IERR,ISYSWR)
      IF (IERR .GT. 0) THEN
            WRITE (ISYSWR,'(A)') ' COMMAND CANNOT BE INTERPRETED'
            ICONDN = 2
            GO TO 900
      ENDIF
C
      CALL MNEXCM(FCN,COMAND(1:LNC),PLIST,LLIST,IERR,FUTIL)
      ICONDN = IERR
  900 RETURN
      END
