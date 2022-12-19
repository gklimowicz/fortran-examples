*
* $Id: mnseek.F,v 1.1.1.1 1996/03/07 14:31:31 mclareni Exp $
*
* $Log: mnseek.F,v $
* Revision 1.1.1.1  1996/03/07 14:31:31  mclareni
* Minuit
*
*
      SUBROUTINE MNSEEK(FCN,FUTIL)
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
CC   Performs a rough (but global) minimization by monte carlo search.
CC        Each time a new minimum is found, the search area is shifted
CC        to be centered at the best value.  Random points are chosen
CC        uniformly over a hypercube determined by current step sizes.
CC   The Metropolis algorithm accepts a worse point with probability
CC      exp(-d/UP), where d is the degradation.  Improved points
CC      are of course always accepted.  Actual steps are random
CC      multiples of the nominal steps (DIRIN).
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
      EXTERNAL FCN,FUTIL
      PARAMETER (TWOPI=2.0*3.141593)
      DIMENSION  XBEST(MNI), XMID(MNI)
      MXFAIL = WORD7(1)
      IF (MXFAIL .LE. 0)  MXFAIL=100+20*NPAR
      MXSTEP = 10*MXFAIL
      IF (AMIN .EQ. UNDEFI)  CALL MNAMIN(FCN,FUTIL)
      ALPHA = WORD7(2)
      IF (ALPHA .LE. ZERO)  ALPHA=3.
      IF (ISW(5) .GE. 1)  WRITE (ISYSWR, 3) MXFAIL,MXSTEP,ALPHA
    3 FORMAT (' MNSEEK: MONTE CARLO MINIMIZATION USING METROPOLIS',
     + ' ALGORITHM'/' TO STOP AFTER',I6,' SUCCESSIVE FAILURES, OR',
     + I7,' STEPS'/' MAXIMUM STEP SIZE IS',F9.3,' ERROR BARS.')
      CSTATU= 'INITIAL  '
      IF (ISW(5) .GE. 2)  CALL MNPRIN(2,AMIN)
      CSTATU = 'UNCHANGED '
      IFAIL = 0
      RNUM = ZERO
      RNUM1 = ZERO
      RNUM2 = ZERO
      NPARX = NPAR
      FLAST = AMIN
C              set up step sizes, starting values
      DO 10 IPAR =  1, NPAR
      IEXT = NEXOFI(IPAR)
      DIRIN(IPAR) = 2.0*ALPHA*WERR(IPAR)
      IF (NVARL(IEXT) .GT. 1)  THEN
C              parameter with limits
         CALL MNDXDI(X(IPAR),IPAR,DXDI)
         IF (DXDI .EQ. ZERO)  DXDI=1.
         DIRIN(IPAR) = 2.0*ALPHA*WERR(IPAR)/DXDI
         IF (ABS(DIRIN(IPAR)).GT.TWOPI)  DIRIN(IPAR)=TWOPI
         ENDIF
      XMID(IPAR) = X(IPAR)
   10 XBEST(IPAR) = X(IPAR)
C                              search loop
      DO 500 ISTEP= 1, MXSTEP
      IF (IFAIL .GE. MXFAIL)  GO TO 600
        DO 100 IPAR= 1, NPAR
        CALL MNRN15(RNUM1,ISEED)
        CALL MNRN15(RNUM2,ISEED)
  100   X(IPAR) = XMID(IPAR) + 0.5*(RNUM1+RNUM2-1.)*DIRIN(IPAR)
      CALL MNINEX(X)
      CALL FCN(NPARX,GIN,FTRY,U,4,FUTIL)
      NFCN = NFCN + 1
      IF (FTRY .LT. FLAST)  THEN
         IF (FTRY .LT. AMIN)  THEN
            CSTATU = 'IMPROVEMNT'
            AMIN = FTRY
            DO 200 IB= 1, NPAR
  200       XBEST(IB) = X(IB)
            IFAIL = 0
            IF (ISW(5) .GE. 2) CALL MNPRIN(2,AMIN)
            ENDIF
         GO TO 300
      ELSE
         IFAIL = IFAIL + 1
C                   Metropolis algorithm
         BAR = (AMIN-FTRY)/UP
         CALL MNRN15(RNUM,ISEED)
         IF (BAR .LT. LOG(RNUM))  GO TO 500
      ENDIF
C                    Accept new point, move there
  300 CONTINUE
      DO 350 J= 1, NPAR
      XMID(J) = X(J)
  350 CONTINUE
      FLAST = FTRY
  500 CONTINUE
C                               end search loop
  600 CONTINUE
      IF (ISW(5) .GT. 1) WRITE (ISYSWR,601) IFAIL
  601 FORMAT(' MNSEEK:',I5,' SUCCESSIVE UNSUCCESSFUL TRIALS.')
      DO 700 IB= 1, NPAR
  700 X(IB) = XBEST(IB)
      CALL MNINEX(X)
      IF (ISW(5) .GE. 1)  CALL MNPRIN(2,AMIN)
      IF (ISW(5) .EQ. 0)  CALL MNPRIN(0,AMIN)
      RETURN
      END
