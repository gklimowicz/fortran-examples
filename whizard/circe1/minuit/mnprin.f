*
* $Id: mnprin.F,v 1.1.1.1 1996/03/07 14:31:31 mclareni Exp $
*
* $Log: mnprin.F,v $
* Revision 1.1.1.1  1996/03/07 14:31:31  mclareni
* Minuit
*
*
      SUBROUTINE MNPRIN  (INKODE,FVAL)
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
CC        Prints the values of the parameters at the time of the call.
CC        also prints other relevant information such as function value,
CC        estimated distance to minimum, parameter errors, step sizes.
CC
C         According to the value of IKODE, the printout is:
C    IKODE=INKODE= 0    only info about function value
C                  1    parameter values, errors, limits
C                  2    values, errors, step sizes, internal values
C                  3    values, errors, step sizes, first derivs.
C                  4    values, parabolic errors, MINOS errors
C    when INKODE=5, MNPRIN chooses IKODE=1,2, or 3, according to ISW(2)
C
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
      CHARACTER*14 COLHDU(6),COLHDL(6), CX2,CX3,CGETX
      CHARACTER*11 CNAMBF, CBLANK
      CHARACTER  CHEDM*10, CHEVAL*15
      PARAMETER (CGETX='PLEASE GET X..')
      DATA CBLANK/'          '/
C
      IF (NU .EQ. 0)  THEN
       WRITE (ISYSWR,'(A)') ' THERE ARE CURRENTLY NO PARAMETERS DEFINED'
       GO TO 700
      ENDIF
C                  get value of IKODE based in INKODE, ISW(2)
      IKODE = INKODE
      IF (INKODE .EQ. 5) THEN
         IKODE = ISW(2)+1
         IF (IKODE .GT. 3)  IKODE=3
      ENDIF
C                  set 'default' column headings
      DO 5 K= 1, 6
      COLHDU(K) = 'UNDEFINED'
    5 COLHDL(K) = 'COLUMN HEAD'
C              print title if Minos errors, and title exists.
      IF (IKODE.EQ.4 .AND. CTITL.NE.CUNDEF)
     +            WRITE (ISYSWR,'(/A,A)')  ' MINUIT TASK: ',CTITL
C              report function value and status
      IF (FVAL .EQ. UNDEFI) THEN
         CHEVAL = ' unknown       '
      ELSE
         WRITE (CHEVAL,'(G15.7)') FVAL
      ENDIF
         IF (EDM .EQ. BIGEDM) THEN
            CHEDM = ' unknown  '
         ELSE
            WRITE (CHEDM, '(E10.2)') EDM
         ENDIF
      NC = NFCN-NFCNFR
      WRITE (ISYSWR,905)  CHEVAL,CFROM,CSTATU,NC,NFCN
  905 FORMAT (/' FCN=',A,' FROM ',A8,'  STATUS=',A10,I6,' CALLS',
     +         I9,' TOTAL')
      M = ISW(2)
      IF (M.EQ.0 .OR. M.EQ.2 .OR. DCOVAR.EQ.ZERO) THEN
        WRITE (ISYSWR,907) CHEDM,ISTRAT,COVMES(M)
  907   FORMAT (21X,'EDM=',A,'    STRATEGY=',I2,6X,A)
      ELSE
        DCMAX = 1.
        DC = MIN(DCOVAR,DCMAX) * 100.
        WRITE (ISYSWR,908) CHEDM,ISTRAT,DC
  908   FORMAT (21X,'EDM=',A,'  STRATEGY=',I1,'  ERROR MATRIX',
     +     ' UNCERTAINTY=',F5.1,'%')
      ENDIF
C
      IF (IKODE .EQ. 0)  GO TO 700
C               find longest name (for Rene!)
      NTRAIL = 10
      DO 20 I= 1, NU
         IF (NVARL(I) .LT. 0)  GO TO 20
         DO 15 IC= 10,1,-1
            IF (CPNAM(I)(IC:IC) .NE. ' ') GO TO 16
   15    CONTINUE
         IC = 1
   16    LBL = 10-IC
         IF (LBL .LT. NTRAIL)  NTRAIL=LBL
   20 CONTINUE
      NADD = NTRAIL/2 + 1
      IF (IKODE .EQ. 1)  THEN
         COLHDU(1) = '              '
         COLHDL(1) = '      ERROR   '
         COLHDU(2) = '      PHYSICAL'
         COLHDU(3) = ' LIMITS       '
         COLHDL(2) = '    NEGATIVE  '
         COLHDL(3) = '    POSITIVE  '
      ENDIF
      IF (IKODE .EQ. 2)  THEN
         COLHDU(1) = '              '
         COLHDL(1) = '      ERROR   '
         COLHDU(2) = '    INTERNAL  '
         COLHDL(2) = '    STEP SIZE '
         COLHDU(3) = '    INTERNAL  '
         COLHDL(3) = '      VALUE   '
      ENDIF
      IF (IKODE .EQ. 3)  THEN
         COLHDU(1) = '              '
         COLHDL(1) = '      ERROR   '
         COLHDU(2) = '       STEP   '
         COLHDL(2) = '       SIZE   '
         COLHDU(3) = '      FIRST   '
         COLHDL(3) = '   DERIVATIVE '
      ENDIF
      IF (IKODE .EQ. 4)  THEN
         COLHDU(1) = '    PARABOLIC '
         COLHDL(1) = '      ERROR   '
         COLHDU(2) = '        MINOS '
         COLHDU(3) = 'ERRORS        '
         COLHDL(2) = '   NEGATIVE   '
         COLHDL(3) = '   POSITIVE   '
      ENDIF
C
      IF (IKODE .NE. 4)  THEN
         IF (ISW(2) .LT. 3) COLHDU(1)='  APPROXIMATE '
         IF (ISW(2) .LT. 1) COLHDU(1)=' CURRENT GUESS'
      ENDIF
      NCOL = 3
      WRITE (ISYSWR, 910) (COLHDU(KK),KK=1,NCOL)
      WRITE (ISYSWR, 911) (COLHDL(KK),KK=1,NCOL)
  910 FORMAT (/'  EXT PARAMETER ',     13X       ,6A14)
  911 FORMAT ( '  NO.   NAME    ','    VALUE    ',6A14)
C
C                                        . . . loop over parameters . .
      DO 200 I= 1, NU
      IF (NVARL(I) .LT. 0)  GO TO 200
      L = NIOFEX(I)
      CNAMBF = CBLANK(1:NADD)//CPNAM(I)
      IF (L .EQ. 0)  GO TO 55
C              variable parameter.
      X1 = WERR(L)
      CX2 = CGETX
      CX3 = CGETX
      IF (IKODE .EQ. 1) THEN
         IF (NVARL(I) .LE. 1) THEN
            WRITE (ISYSWR, 952)  I,CNAMBF,U(I),X1
            GO TO 200
         ELSE
         X2 = ALIM(I)
         X3 = BLIM(I)
         ENDIF
      ENDIF
      IF (IKODE .EQ. 2) THEN
         X2 = DIRIN(L)
         X3 = X(L)
      ENDIF
      IF (IKODE .EQ. 3) THEN
         X2 = DIRIN(L)
         X3 = GRD(L)
         IF (NVARL(I).GT.1 .AND. ABS(COS(X(L))) .LT. 0.001)
     +      CX3 = '** at limit **'
      ENDIF
      IF (IKODE .EQ. 4) THEN
         X2 = ERN(L)
           IF (X2.EQ.ZERO)   CX2=' '
           IF (X2.EQ.UNDEFI) CX2='   at limit   '
         X3 = ERP(L)
           IF (X3.EQ.ZERO)   CX3=' '
           IF (X3.EQ.UNDEFI) CX3='   at limit   '
      ENDIF
      IF (CX2.EQ.CGETX) WRITE (CX2,'(G14.5)') X2
      IF (CX3.EQ.CGETX) WRITE (CX3,'(G14.5)') X3
      WRITE (ISYSWR,952)   I,CNAMBF,U(I),X1,CX2,CX3
  952 FORMAT (I4,1X,A11,2G14.5,2A)
C               check if parameter is at limit
      IF (NVARL(I) .LE. 1 .OR. IKODE .EQ. 3)  GO TO 200
      IF (ABS(COS(X(L))) .LT. 0.001)  WRITE (ISYSWR,1004)
 1004 FORMAT (1H ,32X,42HWARNING -   - ABOVE PARAMETER IS AT LIMIT.)
      GO TO 200
C
C                                print constant or fixed parameter.
   55 CONTINUE
                          COLHDU(1) = '   constant   '
      IF (NVARL(I).GT.0)  COLHDU(1) = '     fixed    '
      IF (NVARL(I).EQ.4 .AND. IKODE.EQ.1) THEN
        WRITE (ISYSWR,'(I4,1X,A11,G14.5,A,2G14.5)')
     +     I,CNAMBF,U(I),COLHDU(1),ALIM(I),BLIM(I)
      ELSE
        WRITE (ISYSWR,'(I4,1X,A11,G14.5,A)')  I,CNAMBF,U(I),COLHDU(1)
      ENDIF
  200 CONTINUE
C
      IF (UP.NE.UPDFLT)  WRITE (ISYSWR,'(31X,A,G10.3)') 'ERR DEF=',UP
  700 CONTINUE
      RETURN
      END
