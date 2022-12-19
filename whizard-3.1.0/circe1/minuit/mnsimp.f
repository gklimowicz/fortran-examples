*
* $Id: mnsimp.F,v 1.2 1996/03/15 18:02:54 james Exp $
*
* $Log: mnsimp.F,v $
* Revision 1.2  1996/03/15 18:02:54  james
*     Modified Files:
* mnderi.F eliminate possible division by zero
* mnexcm.F suppress print on STOP when print flag=-1
*          set FVAL3 to flag if FCN already called with IFLAG=3
* mninit.F set version 96.03
* mnlims.F remove arguments, not needed
* mnmigr.F VLEN -> LENV in debug print statement
* mnparm.F move call to MNRSET to after NPAR redefined, to zero all
* mnpsdf.F eliminate possible division by zero
* mnscan.F suppress printout when print flag =-1
* mnset.F  remove arguments in call to MNLIMS
* mnsimp.F fix CSTATU so status is PROGRESS only if new minimum
* mnvert.F eliminate possible division by zero
*
* Revision 1.1.1.1  1996/03/07 14:31:31  mclareni
* Minuit
*
*
      SUBROUTINE MNSIMP(FCN,FUTIL)
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
CC        Performs a minimization using the simplex method of Nelder
CC        and Mead (ref. -- Comp. J. 7,308 (1965)).
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
      DIMENSION Y(MNI+1)
      DATA ALPHA,BETA,GAMMA,RHOMIN,RHOMAX / 1.0, 0.5, 2.0, 4.0, 8.0/
      IF (NPAR .LE. 0)  RETURN
      IF (AMIN .EQ. UNDEFI)  CALL MNAMIN(FCN,FUTIL)
      CFROM = 'SIMPLEX '
      NFCNFR = NFCN
      CSTATU= 'UNCHANGED '
      NPFN=NFCN
      NPARP1=NPAR+1
      NPARX = NPAR
      RHO1 = 1.0 + ALPHA
      RHO2 = RHO1 + ALPHA*GAMMA
      WG = 1.0/FLOAT(NPAR)
      IF (ISW(5) .GE. 0) WRITE(ISYSWR,100) EPSI
  100 FORMAT(' START SIMPLEX MINIMIZATION.    CONVERGENCE WHEN EDM .LT.'
     +,E10.2 )
         DO 2 I= 1, NPAR
         DIRIN(I) = WERR(I)
           CALL MNDXDI(X(I),I,DXDI)
           IF (DXDI .NE. ZERO) DIRIN(I)=WERR(I)/DXDI
         DMIN = EPSMA2*ABS(X(I))
         IF (DIRIN(I) .LT. DMIN)  DIRIN(I)=DMIN
    2    CONTINUE
C**       choose the initial simplex using single-parameter searches
    1 CONTINUE
      YNPP1 = AMIN
      JL = NPARP1
      Y(NPARP1) = AMIN
      ABSMIN = AMIN
      DO 10 I= 1, NPAR
      AMING = AMIN
      PBAR(I) = X(I)
      BESTX = X(I)
      KG = 0
      NS = 0
      NF = 0
    4 X(I) = BESTX + DIRIN(I)
      CALL MNINEX(X)
      CALL FCN(NPARX,GIN, F, U, 4, FUTIL)
      NFCN = NFCN + 1
      IF (F .LT. AMING)  GO TO 6
C         failure
      IF (KG .EQ. 1)  GO TO 8
      KG = -1
      NF = NF + 1
      DIRIN(I) = DIRIN(I) * (-0.4)
      IF (NF .LT. 3)  GO TO 4
C         stop after three failures
      BESTX = X(I)
      DIRIN(I) = DIRIN(I) * 3.0
      AMING = F
      GO TO 8
C
C         success
    6 BESTX = X(I)
      DIRIN(I) = DIRIN(I) * 3.0
      AMING = F
      CSTATU= 'PROGRESS  '
      KG = 1
      NS = NS + 1
      IF (NS .LT. 6)  GO TO 4
C
C         3 failures or 6 successes or
C         local minimum found in ith direction
    8 Y(I) = AMING
      IF (AMING .LT. ABSMIN)  JL = I
      IF (AMING .LT. ABSMIN)  ABSMIN = AMING
      X(I) = BESTX
      DO 9 K= 1, NPAR
    9 P(K,I) = X(K)
   10 CONTINUE
      JH = NPARP1
      AMIN=Y(JL)
      CALL MNRAZZ(YNPP1,PBAR,Y,JH,JL)
      DO 20 I= 1, NPAR
   20 X(I) = P(I,JL)
      CALL MNINEX(X)
      IF (ISW(5) .GE. 1)  CALL MNPRIN(5,AMIN)
      EDM = BIGEDM
      SIG2 = EDM
      NCYCL=0
C                                        . . . . .  start main loop
   50 CONTINUE
      IF (SIG2 .LT. EPSI .AND. EDM.LT.EPSI)     GO TO 76
      SIG2 = EDM
      IF ((NFCN-NPFN) .GT. NFCNMX)  GO TO 78
C         calculate new point * by reflection
      DO 60 I= 1, NPAR
      PB = 0.
      DO 59 J= 1, NPARP1
   59 PB = PB + WG * P(I,J)
      PBAR(I) = PB - WG * P(I,JH)
   60 PSTAR(I)=(1.+ALPHA)*PBAR(I)-ALPHA*P(I,JH)
      CALL MNINEX(PSTAR)
      CALL FCN(NPARX,GIN,YSTAR,U,4,FUTIL)
      NFCN=NFCN+1
      IF(YSTAR.GE.AMIN) GO TO 70
C         point * better than jl, calculate new point **
      CSTATU = 'PROGRESS  '
      DO 61 I=1,NPAR
   61 PSTST(I)=GAMMA*PSTAR(I)+(1.-GAMMA)*PBAR(I)
      CALL MNINEX(PSTST)
      CALL FCN(NPARX,GIN,YSTST,U,4,FUTIL)
      NFCN=NFCN+1
C         try a parabola through ph, pstar, pstst.  min = prho
      Y1 = (YSTAR-Y(JH)) * RHO2
      Y2 = (YSTST-Y(JH)) * RHO1
      RHO = 0.5 * (RHO2*Y1 -RHO1*Y2) / (Y1 -Y2)
      IF (RHO .LT. RHOMIN)  GO TO 66
      IF (RHO .GT. RHOMAX)  RHO = RHOMAX
      DO 64 I= 1, NPAR
   64 PRHO(I) = RHO*PBAR(I) + (1.0-RHO)*P(I,JH)
      CALL MNINEX(PRHO)
      CALL FCN(NPARX,GIN,YRHO, U,4,FUTIL)
      NFCN = NFCN + 1
      IF (YRHO .LT. AMIN)     CSTATU = 'PROGRESS  '
      IF (YRHO .LT. Y(JL) .AND. YRHO .LT. YSTST)  GO TO 65
      IF (YSTST .LT. Y(JL))  GO TO 67
      IF (YRHO .GT. Y(JL))  GO TO 66
C         accept minimum point of parabola, PRHO
   65 CALL MNRAZZ (YRHO,PRHO,Y,JH,JL)
      GO TO 68
   66 IF (YSTST .LT. Y(JL))  GO TO 67
      CALL MNRAZZ(YSTAR,PSTAR,Y,JH,JL)
      GO TO 68
   67 CALL MNRAZZ(YSTST,PSTST,Y,JH,JL)
   68 NCYCL=NCYCL+1
      IF (ISW(5) .LT. 2)  GO TO 50
      IF (ISW(5) .GE. 3 .OR. MOD(NCYCL, 10) .EQ. 0) CALL MNPRIN(5,AMIN)
      GO TO 50
C         point * is not as good as jl
   70 IF (YSTAR .GE. Y(JH))  GO TO 73
      JHOLD = JH
      CALL MNRAZZ(YSTAR,PSTAR,Y,JH,JL)
      IF (JHOLD .NE. JH)  GO TO 50
C         calculate new point **
   73 DO 74 I=1,NPAR
   74 PSTST(I)=BETA*P(I,JH)+(1.-BETA)*PBAR(I)
      CALL MNINEX (PSTST)
      CALL FCN(NPARX,GIN,YSTST,U,4,FUTIL)
      NFCN=NFCN+1
      IF(YSTST.GT.Y(JH)) GO TO 1
C     point ** is better than jh
      IF (YSTST .LT. AMIN)     CSTATU = 'PROGRESS  '
      IF (YSTST .LT. AMIN)  GO TO 67
      CALL MNRAZZ(YSTST,PSTST,Y,JH,JL)
      GO TO 50
C                                        . . . . . .  end main loop
   76 IF (ISW(5) .GE. 0)  WRITE(ISYSWR,'(A)')
     +                    ' SIMPLEX MINIMIZATION HAS CONVERGED.'
      ISW(4) = 1
      GO TO 80
   78 IF (ISW(5) .GE. 0)  WRITE(ISYSWR,'(A)')
     +                    ' SIMPLEX TERMINATES WITHOUT CONVERGENCE.'
      CSTATU= 'CALL LIMIT'
      ISW(4) = -1
      ISW(1) = 1
   80 DO 82 I=1,NPAR
      PB = 0.
      DO 81 J=1,NPARP1
   81 PB = PB + WG * P(I,J)
   82 PBAR(I) = PB - WG * P(I,JH)
      CALL MNINEX(PBAR)
      CALL FCN(NPARX,GIN,YPBAR,U,4,FUTIL)
      NFCN=NFCN+1
      IF (YPBAR .LT. AMIN)  CALL MNRAZZ(YPBAR,PBAR,Y,JH,JL)
      CALL MNINEX(X)
      IF (NFCNMX+NPFN-NFCN .LT. 3*NPAR)  GO TO 90
      IF (EDM .GT. 2.0*EPSI)  GO TO 1
   90 IF (ISW(5) .GE. 0)  CALL MNPRIN(5, AMIN)
      RETURN
      END
