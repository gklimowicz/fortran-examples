*
* $Id: mnerrs.F,v 1.1.1.1 1996/03/07 14:31:29 mclareni Exp $
*
* $Log: mnerrs.F,v $
* Revision 1.1.1.1  1996/03/07 14:31:29  mclareni
* Minuit
*
*
      SUBROUTINE MNERRS(NUMBER,EPLUS,EMINUS,EPARAB,GCC)
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
CC    Called by user, utility routine to get MINOS errors
CC    If NUMBER is positive, then it is external parameter number,
CC                  if negative, it is -internal number.
CC    values returned by MNERRS:
CC       EPLUS, EMINUS are MINOS errors of parameter NUMBER,
CC       EPARAB is 'parabolic' error (from error matrix).
CC                 (Errors not calculated are set = 0.)
CC       GCC is global correlation coefficient from error matrix
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
      IEX = NUMBER
      IF (NUMBER .LT. 0)  THEN
         IIN = -NUMBER
         IF (IIN .GT. NPAR)  GO TO 900
         IEX = NEXOFI(IIN)
      ENDIF
      IF (IEX .GT. NU .OR. IEX .LE. 0)  GO TO 900
      IIN = NIOFEX(IEX)
      IF (IIN .LE. 0)  GO TO 900
C             IEX is external number, IIN is internal number
      EPLUS = ERP(IIN)
        IF (EPLUS.EQ.UNDEFI)  EPLUS=0.
      EMINUS= ERN(IIN)
        IF (EMINUS.EQ.UNDEFI) EMINUS=0.
      CALL MNDXDI(X(IIN),IIN,DXDI)
      NDIAG = IIN*(IIN+1)/2
      EPARAB = ABS(DXDI*SQRT(ABS(UP*VHMAT(NDIAG))))
C              global correlation coefficient
      GCC = 0.
      IF (ISW(2) .LT. 2)  GO TO 990
      GCC = GLOBCC(IIN)
      GO TO 990
C                  ERROR.  parameter number not valid
  900 EPLUS = 0.
      EMINUS = 0.
      EPARAB = 0.
      GCC = 0.
  990 RETURN
      END
