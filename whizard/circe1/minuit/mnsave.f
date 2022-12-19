*
* $Id: mnsave.F,v 1.1.1.1 1996/03/07 14:31:31 mclareni Exp $
*
* $Log: mnsave.F,v $
* Revision 1.1.1.1  1996/03/07 14:31:31  mclareni
* Minuit
*
*
      SUBROUTINE MNSAVE
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
CC       Writes current parameter values and step sizes onto file ISYSSA
CC          in format which can be reread by Minuit for restarting.
CC       The covariance matrix is also output if it exists.
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
      LOGICAL LOPEN,LNAME
      CHARACTER CGNAME*64, CFNAME*64, CANSWR*1
C
      INQUIRE(UNIT=ISYSSA,OPENED=LOPEN,NAMED=LNAME,NAME=CGNAME)
      IF (LOPEN) THEN
         IF (.NOT.LNAME) CGNAME='UNNAMED FILE'
         WRITE (ISYSWR,32) ISYSSA,CGNAME
   32    FORMAT (' CURRENT VALUES WILL BE SAVED ON UNIT',I3,': ',A/)
      ELSE
C                new file, open it
         WRITE (ISYSWR,35) ISYSSA
   35    FORMAT (' UNIT',I3,' IS NOT OPENED.')
         IF (ISW(6) .EQ. 1) THEN
            WRITE (ISYSWR,'(A)') ' PLEASE GIVE FILE NAME:'
            READ (ISYSRD,'(A)') CFNAME
            OPEN (UNIT=ISYSSA,FILE=CFNAME,STATUS='NEW',ERR=600)
            CGNAME = CFNAME
         ELSE
            GO TO 650
         ENDIF
      ENDIF
C                               file is now correctly opened
      IF (ISW(6) .EQ. 1)  THEN
         WRITE (ISYSWR,37)  ISYSSA
   37    FORMAT (' SHOULD UNIT',I3,' BE REWOUND BEFORE WRITING TO IT?' )
         READ  (ISYSRD,'(A)')  CANSWR
         IF (CANSWR.EQ.'Y' .OR. CANSWR.EQ.'y') REWIND ISYSSA
      ENDIF
C                               and rewound if requested
      WRITE (ISYSSA,'(10HSET TITLE )',ERR=700)
      WRITE (ISYSSA,'(A)')  CTITL
      WRITE (ISYSSA,'(10HPARAMETERS)')
      NLINES = 3
C                                write out parameter values
      DO 200 I= 1, NU
      IF (NVARL(I) .LT. 0)  GO TO 200
      NLINES = NLINES + 1
      IINT = NIOFEX(I)
      IF (NVARL(I) .GT. 1)  GO TO 100
C         parameter without limits
      WRITE (ISYSSA,1001)  I,CPNAM(I),U(I),WERR(IINT)
      GO TO 200
C         parameter with limits
  100 CONTINUE
      WRITE (ISYSSA,1001) I,CPNAM(I),U(I),WERR(IINT),ALIM(I),BLIM(I)
 1001 FORMAT (1X,I5,1H',A10,1H',4E13.5)
  200 CONTINUE
      WRITE (ISYSSA,'(A)')  ' '
      NLINES = NLINES + 1
C                                  write out covariance matrix, if any
      IF (ISW(2) .LT. 1)  GO TO 750
      WRITE (ISYSSA,1003,ERR=700)  NPAR
 1003 FORMAT ('SET COVARIANCE',I6)
      NPAR2 = NPAR*(NPAR+1)/2
      WRITE (ISYSSA,1004) (VHMAT(I),I=1,NPAR2)
 1004 FORMAT (BN,7E11.4,3X)
      NCOVAR = NPAR2/7 + 1
      IF (MOD(NPAR2,7) .GT. 0)  NCOVAR = NCOVAR + 1
      NLINES = NLINES + NCOVAR
      WRITE (ISYSWR, 501) NLINES,ISYSSA,CGNAME(1:45)
  501 FORMAT (1X,I5,' RECORDS WRITTEN TO UNIT',I4,':',A)
      IF (NCOVAR .GT. 0) WRITE (ISYSWR, 502) NCOVAR
  502 FORMAT (' INCLUDING',I5,' RECORDS FOR THE COVARIANCE MATRIX.'/)
      GO TO 900
C                                           some error conditions
  600 WRITE (ISYSWR,'(A,I4)') ' I/O ERROR: UNABLE TO OPEN UNIT',ISYSSA
      GO TO 900
  650 WRITE (ISYSWR,'(A,I4,A)') ' UNIT',ISYSSA,' IS NOT OPENED.'
      GO TO 900
  700 WRITE (ISYSWR,'(A,I4)') ' ERROR: UNABLE TO WRITE TO UNIT',ISYSSA
      GO TO 900
  750 WRITE (ISYSWR,'(A)') ' THERE IS NO COVARIANCE MATRIX TO SAVE.'
C
  900 RETURN
      END
