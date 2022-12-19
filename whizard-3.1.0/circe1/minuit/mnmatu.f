*
* $Id: mnmatu.F,v 1.1.1.1 1996/03/07 14:31:30 mclareni Exp $
*
* $Log: mnmatu.F,v $
* Revision 1.1.1.1  1996/03/07 14:31:30  mclareni
* Minuit
*
*
      SUBROUTINE MNMATU(KODE)
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
CC        prints the covariance matrix v when KODE=1.
CC        always prints the global correlations, and
CC        calculates and prints the individual correlation coefficients
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
      DIMENSION VLINE(MNI)
      ISW2 = ISW(2)
      IF (ISW2 .LT. 1)  THEN
          WRITE (ISYSWR,'(1X,A)')  COVMES(ISW2)
          GO TO 500
      ENDIF
      IF (NPAR .EQ. 0)  THEN
          WRITE (ISYSWR,'('' MNMATU: NPAR=0'')')
          GO TO 500
          ENDIF
C                                       . . . . .external error matrix
      IF (KODE .EQ. 1)  THEN
         ISW5 = ISW(5)
         ISW(5) = 2
         CALL MNEMAT(P,MAXINT)
           IF (ISW2.LT.3)  WRITE (ISYSWR,'(1X,A)')  COVMES(ISW2)
         ISW(5) = ISW5
      ENDIF
C                                       . . . . . correlation coeffs. .
      IF (NPAR .LE. 1)   GO TO 500
      CALL MNWERR
C     NCOEF is number of coeff. that fit on one line, not to exceed 20
      NCOEF = (NPAGWD-19)/6
      NCOEF = MIN(NCOEF,20)
      NPARM = MIN(NPAR,NCOEF)
      WRITE (ISYSWR, 150) (NEXOFI(ID),ID=1,NPARM)
  150 FORMAT (/36H PARAMETER  CORRELATION COEFFICIENTS  /
     +         18H       NO.  GLOBAL   ,20I6)
      DO 200 I= 1, NPAR
         IX = NEXOFI(I)
         NDI = I*(I+1)/2
           DO 170 J= 1, NPAR
           M = MAX(I,J)
           N = MIN(I,J)
           NDEX = M*(M-1)/2 + N
           NDJ = J*(J+1)/2
  170      VLINE(J) = VHMAT(NDEX)/SQRT(ABS(VHMAT(NDI)*VHMAT(NDJ)))
         NPARM = MIN(NPAR,NCOEF)
         WRITE (ISYSWR,171)   IX, GLOBCC(I), (VLINE(IT),IT=1,NPARM)
  171    FORMAT (6X,I3,2X,F7.5,1X,20F6.3)
         IF (I.LE.NPARM) GO TO 200
            DO 190 ISO= 1, 10
            NSOFAR = NPARM
            NPARM = MIN(NPAR,NSOFAR+NCOEF)
            WRITE (ISYSWR,181)  (VLINE(IT),IT=NSOFAR+1,NPARM)
  181       FORMAT (19X,20F6.3)
            IF (I .LE. NPARM) GO TO 192
  190       CONTINUE
  192    CONTINUE
  200 CONTINUE
      IF (ISW2.LT.3)  WRITE (ISYSWR,'(1X,A)')  COVMES(ISW2)
  500 RETURN
      END
