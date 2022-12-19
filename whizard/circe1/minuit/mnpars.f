*
* $Id: mnpars.F,v 1.1.1.1 1996/03/07 14:31:31 mclareni Exp $
*
* $Log: mnpars.F,v $
* Revision 1.1.1.1  1996/03/07 14:31:31  mclareni
* Minuit
*
*
      SUBROUTINE MNPARS(CRDBUF,ICONDN)
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
CC        Called from MNREAD and user-callable
CC    Implements one parameter definition, that is:
CC       parses the string CRDBUF and calls MNPARM
C
C output conditions:
C        ICONDN = 0    all OK
C        ICONDN = 1    error, attempt to define parameter is ignored
C        ICONDN = 2    end of parameter definitions
C
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
C
      DIMENSION PLIST(MAXP)
      CHARACTER CNAMK*10, CRDBUF*(*) , CELMNT*20 , COMAND*(MAXCWD)
C
      LENBUF = LEN(CRDBUF)
C                     find out whether fixed or free-field format
      KAPO1 = INDEX(CRDBUF,'''')
      IF (KAPO1 .EQ. 0)  GO TO 150
      KAPO2 = INDEX(CRDBUF(KAPO1+1:),'''')
      IF (KAPO2 .EQ. 0)  GO TO 150
C          new (free-field) format
      KAPO2 = KAPO2 + KAPO1
C                             skip leading blanks if any
         DO 115 ISTART=1, KAPO1-1
         IF (CRDBUF(ISTART:ISTART) .NE. ' ')  GO TO 120
  115    CONTINUE
         GO TO 210
  120 CONTINUE
C                               parameter number integer
      CELMNT = CRDBUF(ISTART:KAPO1-1)
      READ (CELMNT,'(BN,F20.0)',ERR=180) FK
      K = FK
      IF (K .LE. 0)  GO TO 210
      CNAMK = 'PARAM '//CELMNT
      IF (KAPO2-KAPO1 .GT. 1) CNAMK = CRDBUF(KAPO1+1:KAPO2-1)
C  special handling if comma or blanks and a comma follow 'name'
        DO 135 ICY= KAPO2+1,LENBUF
        IF (CRDBUF(ICY:ICY) .EQ. ',') GO TO 139
        IF (CRDBUF(ICY:ICY) .NE. ' ') GO TO 140
  135 CONTINUE
        UK = 0.
        WK = 0.
        A  = 0.
        B = 0.
      GO TO 170
  139 CONTINUE
      ICY = ICY+1
  140 CONTINUE
      IBEGIN = ICY
      CALL MNCRCK(CRDBUF(IBEGIN:),MAXCWD,COMAND,LNC,
     +                             MAXP,PLIST,LLIST, IERR,ISYSWR)
      IF (IERR .GT. 0)  GO TO 180
      UK = PLIST(1)
      WK = 0.
      IF (LLIST .GE. 2)  WK = PLIST(2)
      A = 0.
      IF (LLIST .GE. 3)  A = PLIST(3)
      B = 0.
      IF (LLIST .GE. 4)  B = PLIST(4)
      GO TO 170
C          old (fixed-field) format
  150 CONTINUE
      READ (CRDBUF, 158,ERR=180)  XK,CNAMK,UK,WK,A,B
  158 FORMAT (BN,F10.0, A10, 4F10.0)
      K = XK
      IF (K .EQ. 0)  GO TO 210
C          parameter format cracked, implement parameter definition
  170 CALL MNPARM(K,CNAMK,UK,WK,A,B,IERR)
      ICONDN = IERR
      RETURN
C          format or other error
  180 CONTINUE
      ICONDN = 1
      RETURN
C        end of data
  210 CONTINUE
      ICONDN = 2
      RETURN
      END
