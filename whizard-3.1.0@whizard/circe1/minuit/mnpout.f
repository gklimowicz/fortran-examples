*
* $Id: mnpout.F,v 1.1.1.1 1996/03/07 14:31:31 mclareni Exp $
*
* $Log: mnpout.F,v $
* Revision 1.1.1.1  1996/03/07 14:31:31  mclareni
* Minuit
*
*
      SUBROUTINE MNPOUT(IUEXT,CHNAM,VAL,ERR,XLOLIM,XUPLIM,IUINT)
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
CC     User-called
CC   Provides the user with information concerning the current status
CC          of parameter number IUEXT. Namely, it returns:
CC        CHNAM: the name of the parameter
CC        VAL: the current (external) value of the parameter
CC        ERR: the current estimate of the parameter uncertainty
CC        XLOLIM: the lower bound (or zero if no limits)
CC        XUPLIM: the upper bound (or zero if no limits)
CC        IUINT: the internal parameter number (or zero if not variable,
CC           or negative if undefined).
CC  Note also:  If IUEXT is negative, then it is -internal parameter
CC           number, and IUINT is returned as the EXTERNAL number.
CC     Except for IUINT, this is exactly the inverse of MNPARM
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
      CHARACTER*(*) CHNAM
      XLOLIM = 0.
      XUPLIM = 0.
      ERR = 0.
      IF (IUEXT .EQ. 0)  GO TO 100
      IF (IUEXT .LT. 0)  THEN
C                   internal parameter number specified
         IINT = -IUEXT
         IF (IINT .GT. NPAR) GO TO 100
         IEXT = NEXOFI(IINT)
         IUINT = IEXT
      ELSE
C                    external parameter number specified
         IEXT = IUEXT
         IF (IEXT .EQ. 0)   GO TO 100
         IF (IEXT .GT. NU)  GO TO 100
         IINT = NIOFEX(IEXT)
         IUINT = IINT
      ENDIF
C                     in both cases
         NVL = NVARL(IEXT)
         IF (NVL .LT. 0) GO TO 100
      CHNAM = CPNAM(IEXT)
      VAL = U(IEXT)
      IF (IINT .GT. 0)  ERR = WERR(IINT)
      IF (NVL .EQ. 4) THEN
         XLOLIM = ALIM(IEXT)
         XUPLIM = BLIM(IEXT)
      ENDIF
      RETURN
C                parameter is undefined
  100 IUINT = -1
      CHNAM = 'undefined'
      VAL = 0.
      RETURN
      END
