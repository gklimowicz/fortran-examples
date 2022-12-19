*
* $Id: mnhes1.F,v 1.1.1.1 1996/03/07 14:31:30 mclareni Exp $
*
* $Log: mnhes1.F,v $
* Revision 1.1.1.1  1996/03/07 14:31:30  mclareni
* Minuit
*
*
      SUBROUTINE MNHES1(FCN,FUTIL)
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
CC      Called from MNHESS and MNGRAD
CC      Calculate first derivatives (GRD) and uncertainties (DGRD)
CC         and appropriate step sizes GSTEP
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
      EXTERNAL FCN,FUTIL
      LOGICAL LDEBUG
      CHARACTER CBF1*22
      LDEBUG = (IDBG(5) .GE. 1)
      IF (ISTRAT .LE. 0) NCYC = 1
      IF (ISTRAT .EQ. 1) NCYC = 2
      IF (ISTRAT .GT. 1) NCYC = 6
      IDRV = 1
      NPARX = NPAR
      DFMIN = 4.*EPSMA2*(ABS(AMIN)+UP)
C                                     main loop over parameters
      DO 100 I= 1, NPAR
      XTF = X(I)
      DMIN = 4.*EPSMA2*ABS(XTF)
      EPSPRI = EPSMA2 + ABS(GRD(I)*EPSMA2)
      OPTSTP = SQRT(DFMIN/(ABS(G2(I))+EPSPRI))
      D = 0.2 * ABS(GSTEP(I))
      IF (D .GT. OPTSTP)  D = OPTSTP
      IF (D .LT. DMIN)  D = DMIN
      CHGOLD = 10000.
C                                       iterate reducing step size
      DO 50 ICYC= 1, NCYC
      X(I) = XTF + D
      CALL MNINEX(X)
      CALL FCN(NPARX,GIN,FS1,U,4,FUTIL)
      NFCN = NFCN + 1
      X(I) = XTF - D
      CALL MNINEX(X)
      CALL FCN(NPARX,GIN,FS2,U,4,FUTIL)
      NFCN = NFCN + 1
      X(I) = XTF
C                                       check if step sizes appropriate
      SAG = 0.5*(FS1+FS2-2.0*AMIN)
      GRDOLD = GRD(I)
      GRDNEW = (FS1-FS2)/(2.0*D)
      DGMIN = EPSMAC*(ABS(FS1)+ABS(FS2))/D
      IF (LDEBUG) WRITE (ISYSWR,11) I,IDRV,GSTEP(I),D,G2(I),GRDNEW,SAG
   11 FORMAT (I4,I2,6G12.5)
      IF (GRDNEW .EQ. ZERO)  GO TO 60
      CHANGE = ABS((GRDOLD-GRDNEW)/GRDNEW)
      IF (CHANGE.GT.CHGOLD .AND. ICYC.GT.1)  GO TO 60
      CHGOLD = CHANGE
      GRD(I) = GRDNEW
      GSTEP(I) = SIGN(D,GSTEP(I))
C                  decrease step until first derivative changes by <5%
      IF (CHANGE .LT. 0.05) GO TO 60
      IF (ABS(GRDOLD-GRDNEW) .LT. DGMIN)  GO TO 60
      IF (D .LT. DMIN)  THEN
         CALL MNWARN('D','MNHES1','Step size too small for 1st drv.')
         GO TO 60
      ENDIF
      D = 0.2*D
   50 CONTINUE
C                                       loop satisfied = too many iter
      WRITE (CBF1,'(2G11.3)') GRDOLD,GRDNEW
      CALL MNWARN('D','MNHES1','Too many iterations on D1.'//CBF1)
   60 CONTINUE
      DGRD(I) = MAX(DGMIN,ABS(GRDOLD-GRDNEW))
  100 CONTINUE
C                                        end of first deriv. loop
      CALL MNINEX(X)
      RETURN
      END
