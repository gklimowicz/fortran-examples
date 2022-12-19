*
* $Id: mnemat.F,v 1.1.1.1 1996/03/07 14:31:29 mclareni Exp $
*
* $Log: mnemat.F,v $
* Revision 1.1.1.1  1996/03/07 14:31:29  mclareni
* Minuit
*
*
      SUBROUTINE MNEMAT(EMAT,NDIM)
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
      DIMENSION EMAT(NDIM,NDIM)
CC        Calculates the external error matrix from the internal
CC        to be called by user, who must dimension EMAT at (NDIM,NDIM)
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
      IF (ISW(2) .LT. 1)  RETURN
      IF (ISW(5) .GE. 2)  WRITE (ISYSWR,'(/A,I4,A,I3,A,G10.3)')
     +    ' EXTERNAL ERROR MATRIX.    NDIM=',NDIM,'    NPAR=',NPAR,
     +    '    ERR DEF=',UP
C                    size of matrix to be printed
      NPARD = NPAR
      IF (NDIM .LT. NPAR)  THEN
        NPARD = NDIM
        IF (ISW(5) .GE. 0) WRITE (ISYSWR,'(A,A)') ' USER-DIMENSIONED ',
     +      ' ARRAY EMAT NOT BIG ENOUGH. REDUCED MATRIX CALCULATED.'
      ENDIF
C                 NPERLN is the number of elements that fit on one line
      NPERLN = (NPAGWD-5)/10
      NPERLN = MIN(NPERLN,13)
      IF (ISW(5).GE. 1 .AND. NPARD.GT.NPERLN)  WRITE (ISYSWR,'(A)')
     +     ' ELEMENTS ABOVE DIAGONAL ARE NOT PRINTED.'
C                 I counts the rows of the matrix
      DO 110 I= 1, NPARD
         CALL MNDXDI(X(I),I,DXDI)
         KGA = I*(I-1)/2
         DO 100 J= 1, I
            CALL MNDXDI(X(J),J,DXDJ)
            KGB = KGA + J
            EMAT(I,J) = DXDI * VHMAT(KGB) * DXDJ * UP
            EMAT(J,I) = EMAT(I,J)
  100    CONTINUE
  110 CONTINUE
C                    IZ is number of columns to be printed in row I
      IF (ISW(5) .GE. 2)  THEN
      DO 160 I= 1, NPARD
         IZ = NPARD
         IF (NPARD .GE. NPERLN)  IZ = I
         DO 150 K= 1, IZ, NPERLN
           K2 = K + NPERLN - 1
           IF (K2 .GT. IZ)  K2=IZ
           WRITE (ISYSWR,'(1X,13E10.3)')  (EMAT(I,KK),KK=K,K2)
  150    CONTINUE
  160 CONTINUE
      ENDIF
      RETURN
      END
