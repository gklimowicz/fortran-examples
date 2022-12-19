*
* $Id: mngrad.F,v 1.1.1.1 1996/03/07 14:31:29 mclareni Exp $
*
* $Log: mngrad.F,v $
* Revision 1.1.1.1  1996/03/07 14:31:29  mclareni
* Minuit
*
*
      SUBROUTINE MNGRAD(FCN,FUTIL)
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
CC       Called from MNSET
CC       Interprets the SET GRAD command, which informs MINUIT whether
CC       the first derivatives of FCN will be calculated by the user
CC       inside FCN.  It can check the user's derivative calculation
CC       by comparing it with a finite difference approximation.
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
C
      EXTERNAL FCN,FUTIL
      CHARACTER*4 CGOOD,CBAD,CNONE,CWD
      LOGICAL LNONE
      DIMENSION GF(MNI)
      PARAMETER (CGOOD='GOOD',CBAD=' BAD',CNONE='NONE')
C
      ISW(3) = 1
      NPARX = NPAR
      IF (WORD7(1) .GT. ZERO)  GO TO 2000
C                  get user-calculated first derivatives from FCN
      DO 30 I= 1, NU
   30 GIN(I) = UNDEFI
      CALL MNINEX(X)
      CALL FCN(NPARX,GIN,FZERO,U,2,FUTIL)
      NFCN = NFCN + 1
      CALL MNDERI(FCN,FUTIL)
      DO 40 I= 1, NPAR
   40 GF(I) = GRD(I)
C                    get MINUIT-calculated first derivatives
      ISW(3) = 0
      ISTSAV = ISTRAT
      ISTRAT = 2
      CALL MNHES1(FCN,FUTIL)
      ISTRAT = ISTSAV
      WRITE (ISYSWR,51)
   51 FORMAT(/' CHECK OF GRADIENT CALCULATION IN FCN'/12X,'PARAMETER',
     + 6X,9HG(IN FCN) ,3X,9HG(MINUIT) ,2X,'DG(MINUIT)',3X,9HAGREEMENT)
      ISW(3) = 1
      LNONE = .FALSE.
      DO 100 LC = 1, NPAR
      I = NEXOFI(LC)
      CWD = CGOOD
      ERR = DGRD(LC)
      IF (ABS(GF(LC)-GRD(LC)) .GT. ERR)  CWD = CBAD
      IF (GIN(I) .EQ. UNDEFI)  THEN
          CWD = CNONE
          LNONE = .TRUE.
          GF(LC) = 0.
          ENDIF
      IF (CWD .NE. CGOOD)  ISW(3) = 0
      WRITE (ISYSWR,99) I,CPNAM(I),GF(LC),GRD(LC),ERR,CWD
   99 FORMAT (7X,I5,2X ,A10,3E12.4,4X ,A4)
  100 CONTINUE
      IF (LNONE) WRITE (ISYSWR,'(A)')
     +  '  AGREEMENT=NONE  MEANS FCN DID NOT CALCULATE THE DERIVATIVE'
      IF (ISW(3) .EQ. 0)  WRITE (ISYSWR,1003)
 1003 FORMAT(/' MINUIT DOES NOT ACCEPT DERIVATIVE CALCULATIONS BY FCN'/
     + ' TO FORCE ACCEPTANCE, ENTER "SET GRAD    1"'/)
C
 2000 CONTINUE
      RETURN
      END
