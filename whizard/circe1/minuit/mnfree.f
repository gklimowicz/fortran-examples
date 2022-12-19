*
* $Id: mnfree.F,v 1.1.1.1 1996/03/07 14:31:29 mclareni Exp $
*
* $Log: mnfree.F,v $
* Revision 1.1.1.1  1996/03/07 14:31:29  mclareni
* Minuit
*
*
      SUBROUTINE MNFREE(K)
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
CC        Restores one or more fixed parameter(s) to variable status
CC        by inserting it into the internal parameter list at the
CC        appropriate place.
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
C--       K = 0 means restore all parameters
C--       K = 1 means restore the last parameter fixed
C--       K = -I means restore external parameter I (if possible)
C--       IQ = fix-location where internal parameters were stored
C--       IR = external number of parameter being restored
C--       IS = internal number of parameter being restored
      IF (K .GT. 1)  WRITE (ISYSWR,510)
      IF (NPFIX .LT. 1)  WRITE (ISYSWR,500)
      IF (K.EQ.1 .OR. K.EQ.0)  GO TO 40
C                   release parameter with specified external number
      KA = IABS(K)
      IF (NIOFEX(KA) .EQ. 0)  GO TO 15
      WRITE (ISYSWR,540)
  540 FORMAT (' IGNORED.  PARAMETER SPECIFIED IS ALREADY VARIABLE.')
      RETURN
   15 IF (NPFIX .LT. 1)  GO TO 21
      DO 20 IK= 1, NPFIX
      IF (IPFIX(IK) .EQ. KA)  GO TO 24
   20 CONTINUE
   21 WRITE (ISYSWR,530) KA
  530 FORMAT (' PARAMETER',I4,' NOT FIXED.  CANNOT BE RELEASED.')
      RETURN
   24 IF (IK .EQ. NPFIX)  GO TO 40
C                   move specified parameter to end of list
      IPSAV = KA
      XV = XS(IK)
      XTV = XTS(IK)
      DIRINV = DIRINS(IK)
      GRDV = GRDS(IK)
      G2V = G2S(IK)
      GSTEPV = GSTEPS(IK)
         DO 30 I= IK+1,NPFIX
         IPFIX(I-1) = IPFIX(I)
         XS(I-1) = XS(I)
         XTS(I-1) = XTS(I)
         DIRINS(I-1) = DIRINS(I)
         GRDS(I-1) = GRDS(I)
         G2S(I-1) = G2S(I)
         GSTEPS(I-1) = GSTEPS(I)
   30    CONTINUE
      IPFIX(NPFIX) = IPSAV
      XS(NPFIX) = XV
      XTS(NPFIX) = XTV
      DIRINS(NPFIX) = DIRINV
      GRDS(NPFIX) = GRDV
      G2S(NPFIX) = G2V
      GSTEPS(NPFIX) = GSTEPV
C                restore last parameter in fixed list  -- IPFIX(NPFIX)
   40 CONTINUE
      IF (NPFIX .LT. 1)  GO TO 300
      IR = IPFIX(NPFIX)
      IS = 0
      DO 100 IK= NU, IR, -1
        IF (NIOFEX(IK) .GT. 0) THEN
         LC = NIOFEX(IK) + 1
         IS = LC - 1
         NIOFEX(IK) = LC
         NEXOFI(LC) = IK
         X(LC)     = X(LC-1)
         XT(LC)    = XT(LC-1)
         DIRIN(LC) = DIRIN(LC-1)
         WERR(LC)  = WERR(LC-1)
         GRD(LC)   = GRD(LC-1)
         G2(LC)    = G2(LC-1)
         GSTEP(LC) = GSTEP(LC-1)
        ENDIF
  100 CONTINUE
      NPAR = NPAR + 1
      IF (IS .EQ. 0)   IS = NPAR
      NIOFEX(IR) = IS
      NEXOFI(IS) = IR
      IQ = NPFIX
      X(IS) = XS(IQ)
      XT(IS) = XTS(IQ)
      DIRIN(IS) = DIRINS(IQ)
      WERR(IS)  = DIRINS(IQ)
      GRD(IS) = GRDS(IQ)
      G2(IS) = G2S(IQ)
      GSTEP(IS) = GSTEPS(IQ)
      NPFIX = NPFIX - 1
      ISW(2) = 0
      DCOVAR = 1.
      IF (ISW(5)-ITAUR .GE. 1)  WRITE(ISYSWR,520) IR,CPNAM(IR)
      IF (K.EQ.0)  GO TO 40
  300 CONTINUE
C         if different from internal, external values are taken
      CALL MNEXIN(X)
  400 RETURN
  500 FORMAT (' CALL TO MNFREE IGNORED.  THERE ARE NO FIXED PA',
     + 'RAMETERS'/)
  510 FORMAT (' CALL TO MNFREE IGNORED.  ARGUMENT GREATER THAN ONE'/)
  520 FORMAT (20X, 9HPARAMETER,I4,2H, ,A10,' RESTORED TO VARIABLE.')
      END
