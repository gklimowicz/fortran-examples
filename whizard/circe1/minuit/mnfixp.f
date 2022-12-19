*
* $Id: mnfixp.F,v 1.1.1.1 1996/03/07 14:31:29 mclareni Exp $
*
* $Log: mnfixp.F,v $
* Revision 1.1.1.1  1996/03/07 14:31:29  mclareni
* Minuit
*
*
      SUBROUTINE MNFIXP(IINT,IERR)
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
CC        removes parameter IINT from the internal (variable) parameter
CC        list, and arranges the rest of the list to fill the hole.
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
      DIMENSION YY(MNI)
C                           first see if it can be done
      IERR = 0
      IF (IINT.GT.NPAR .OR. IINT.LE.0)  THEN
         IERR = 1
         WRITE (ISYSWR,'(A,I4)')
     +       ' MINUIT ERROR.  ARGUMENT TO MNFIXP=',IINT
         GO TO 300
      ENDIF
      IEXT = NEXOFI(IINT)
      IF (NPFIX .GE. MNI) THEN
         IERR = 1
         WRITE (ISYSWR,'(A,I4,A,I4)') ' MINUIT CANNOT FIX PARAMETER',
     +   IEXT,' MAXIMUM NUMBER THAT CAN BE FIXED IS',MNI
         GO TO 300
      ENDIF
C                           reduce number of variable parameters by one
      NIOFEX(IEXT) = 0
      NOLD = NPAR
      NPAR = NPAR - 1
C                       save values in case parameter is later restored
      NPFIX = NPFIX + 1
      IPFIX(NPFIX) = IEXT
      LC = IINT
      XS(NPFIX) = X(LC)
      XTS(NPFIX) = XT(LC)
      DIRINS(NPFIX) = WERR(LC)
      GRDS(NPFIX) = GRD(LC)
      G2S(NPFIX) = G2(LC)
      GSTEPS(NPFIX) = GSTEP(LC)
C                        shift values for other parameters to fill hole
      DO 100  IK= IEXT+1, NU
         IF  (NIOFEX(IK) .GT. 0)  THEN
         LC = NIOFEX(IK) - 1
         NIOFEX(IK) = LC
         NEXOFI(LC) = IK
         X(LC)     = X(LC+1)
         XT(LC)    = XT(LC+1)
         DIRIN(LC) = DIRIN(LC+1)
         WERR(LC)  = WERR(LC+1)
         GRD(LC)   = GRD(LC+1)
         G2(LC)    = G2(LC+1)
         GSTEP(LC) = GSTEP(LC+1)
         ENDIF
  100 CONTINUE
      IF (ISW(2) .LE. 0)  GO TO 300
C                    remove one row and one column from variance matrix
      IF (NPAR .LE. 0)  GO TO 300
      DO 260 I= 1, NOLD
      M = MAX(I,IINT)
      N = MIN(I,IINT)
      NDEX = M*(M-1)/2 + N
  260 YY(I)=VHMAT(NDEX)
      YYOVER = 1.0/YY(IINT)
      KNEW = 0
      KOLD = 0
      DO 294 I= 1, NOLD
      DO 292 J= 1, I
      KOLD = KOLD + 1
      IF (J.EQ.IINT .OR. I.EQ.IINT)  GO TO 292
      KNEW = KNEW + 1
      VHMAT(KNEW) = VHMAT(KOLD) - YY(J)*YY(I)*YYOVER
  292 CONTINUE
  294 CONTINUE
  300 RETURN
      END
