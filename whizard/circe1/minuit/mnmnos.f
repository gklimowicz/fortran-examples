*
* $Id: mnmnos.F,v 1.1.1.1 1996/03/07 14:31:30 mclareni Exp $
*
* $Log: mnmnos.F,v $
* Revision 1.1.1.1  1996/03/07 14:31:30  mclareni
* Minuit
*
*
      SUBROUTINE MNMNOS(FCN,FUTIL)
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
CC        Performs a MINOS error analysis on those parameters for
CC        which it is requested on the MINOS command by calling 
CC        MNMNOT for each parameter requested.
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
      IF (NPAR .LE. 0)  GO TO 700
      NGOOD = 0
      NBAD = 0
      NFCNMI = NFCN
C                                      . loop over parameters requested
      DO 570 KNT= 1, NPAR
      IF (INT(WORD7(2)) .EQ. 0) THEN
          ILAX = NEXOFI(KNT)
      ELSE
          IF (KNT .GE. 7)  GO TO 580
          ILAX = INT(WORD7(KNT+1))
          IF (ILAX .EQ. 0)  GO TO 580
          IF (ILAX .GT. 0 .AND. ILAX .LE. NU) THEN
             IF (NIOFEX(ILAX) .GT. 0)  GO TO 565
          ENDIF
          WRITE (ISYSWR,564) ILAX
  564     FORMAT (' PARAMETER NUMBER ',I5,' NOT VARIABLE. IGNORED.')
          GO TO 570
      ENDIF
  565 CONTINUE
C                                         calculate one pair of M E's
      ILAX2 = 0
      CALL MNMNOT(FCN,ILAX,ILAX2,VAL2PL,VAL2MI,FUTIL)
      IF (LNEWMN)  GO TO 650
C                                          update NGOOD and NBAD
      IIN = NIOFEX(ILAX)
      IF (ERP(IIN) .GT. ZERO) THEN
         NGOOD=NGOOD+1
      ELSE
         NBAD=NBAD+1
      ENDIF
      IF (ERN(IIN) .LT. ZERO) THEN
         NGOOD=NGOOD+1
      ELSE
         NBAD=NBAD+1
      ENDIF
  570 CONTINUE
C                                           end of loop . . . . . . .
  580 CONTINUE
C                                        . . . . printout final values .
      CFROM = 'MINOS   '
      NFCNFR = NFCNMI
      CSTATU= 'UNCHANGED '
      IF (NGOOD.EQ.0.AND.NBAD.EQ.0) GO TO 700
      IF (NGOOD.GT.0.AND.NBAD.EQ.0) CSTATU='SUCCESSFUL'
      IF (NGOOD.EQ.0.AND.NBAD.GT.0) CSTATU='FAILURE   '
      IF (NGOOD.GT.0.AND.NBAD.GT.0) CSTATU='PROBLEMS  '
      IF (ISW(5) .GE. 0) CALL MNPRIN(4,AMIN)
      IF (ISW(5) .GE. 2) CALL MNMATU(0)
      GO TO 900
C                                        . . . new minimum found . . . .
  650 CONTINUE
      CFROM = 'MINOS   '
      NFCNFR = NFCNMI
      CSTATU= 'NEW MINIMU'
      IF (ISW(5) .GE. 0) CALL MNPRIN(4,AMIN)
      WRITE (ISYSWR,675)
  675 FORMAT(/50H NEW MINIMUM FOUND.  GO BACK TO MINIMIZATION STEP./1H ,
     +60(1H=)/60X,1HV/60X,1HV/60X,1HV/57X,7HVVVVVVV/58X,5HVVVVV/59X,
     +3HVVV/60X,1HV//)
      GO TO 900
  700 WRITE (ISYSWR,'(A)') ' THERE ARE NO MINOS ERRORS TO CALCULATE.'
  900 RETURN
      END
