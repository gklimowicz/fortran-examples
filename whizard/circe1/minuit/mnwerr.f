*
* $Id: mnwerr.F,v 1.1.1.1 1996/03/07 14:31:32 mclareni Exp $
*
* $Log: mnwerr.F,v $
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
      SUBROUTINE MNWERR
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
CC          Calculates the WERR, external parameter errors,
CC      and the global correlation coefficients, to be called
CC      whenever a new covariance matrix is available.
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
C                         calculate external error if v exists
      IF (ISW(2) .GE. 1) THEN
      DO 100 L= 1, NPAR
        NDEX = L*(L+1)/2
        DX = SQRT(ABS(VHMAT(NDEX)*UP))
        I = NEXOFI(L)
        IF (NVARL(I) .GT. 1)  THEN
          AL = ALIM(I)
          BA = BLIM(I) - AL
          DU1 = AL + 0.5 *(SIN(X(L)+DX) +1.0) * BA - U(I)
          DU2 = AL + 0.5 *(SIN(X(L)-DX) +1.0) * BA - U(I)
          IF (DX .GT. 1.0)  DU1 = BA
          DX = 0.5 * (ABS(DU1) + ABS(DU2))
        ENDIF
        WERR(L) = DX
  100 CONTINUE
      ENDIF
C                          global correlation coefficients
      IF (ISW(2) .GE. 1) THEN
         DO 130 I= 1, NPAR
            GLOBCC(I) = 0.
            K1 = I*(I-1)/2
            DO 130 J= 1, I
               K = K1 + J
               P(I,J) = VHMAT(K)
  130          P(J,I) = P(I,J)
         CALL MNVERT(P,MAXINT,MAXINT,NPAR,IERR)
         IF (IERR .EQ. 0)   THEN
            DO 150 IIN= 1, NPAR
               NDIAG = IIN*(IIN+1)/2
               DENOM = P(IIN,IIN)*VHMAT(NDIAG)
               IF (DENOM.LE.ONE .AND. DENOM.GE.ZERO)  THEN
                   GLOBCC(IIN) = 0.
               ELSE
                   GLOBCC(IIN) = SQRT(1.0-1.0/DENOM)
               ENDIF
  150       CONTINUE
         ENDIF
      ENDIF
      RETURN
      END
