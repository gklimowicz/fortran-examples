*
* $Id: mnwarn.F,v 1.1.1.1 1996/03/07 14:31:32 mclareni Exp $
*
* $Log: mnwarn.F,v $
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
      SUBROUTINE MNWARN(COPT,CORG,CMES)
C     If COPT='W', CMES is a WARning message from CORG.
C     If COPT='D', CMES is a DEBug message from CORG.
C         If SET WARnings is in effect (the default), this routine
C             prints the warning message CMES coming from CORG.
C         If SET NOWarnings is in effect, the warning message is
C             stored in a circular buffer of length MAXMES.
C         If called with CORG=CMES='SHO', it prints the messages in
C             the circular buffer, FIFO, and empties the buffer.
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
      CHARACTER COPT*1, CORG*(*), CMES*(*), CTYP*7
      PARAMETER (MAXMES=10)
      CHARACTER     ORIGIN(MAXMES,2)*10, WARMES(MAXMES,2)*60
      COMMON/MN7WRC/ORIGIN,              WARMES
      COMMON/MN7WRI/NFCWAR(MAXMES,2),ICIRC(2)
      CHARACTER ENGLSH*20
C
      IF (CORG(1:3).EQ.'SHO' .AND. CMES(1:3).EQ.'SHO')  GO TO 200
C             Either print warning or put in buffer
      IF (COPT .EQ. 'W')  THEN
        ITYP = 1
        IF (LWARN) THEN
          WRITE (ISYSWR,'(A,A/A,A)') ' MINUIT WARNING IN ',CORG,
     +              ' ============== ',CMES
          RETURN
        ENDIF
      ELSE
        ITYP = 2
        IF (LREPOR) THEN
          WRITE (ISYSWR,'(A,A/A,A)') ' MINUIT DEBUG FOR  ',CORG,
     +              ' ============== ',CMES
          RETURN
        ENDIF
      ENDIF
C                 if appropriate flag is off, fill circular buffer
         IF (NWRMES(ITYP) .EQ. 0)  ICIRC(ITYP) = 0
         NWRMES(ITYP) = NWRMES(ITYP) + 1
         ICIRC(ITYP) = ICIRC(ITYP) + 1
         IF (ICIRC(ITYP) .GT. MAXMES) ICIRC(ITYP) = 1
         IC = ICIRC(ITYP)
         ORIGIN(IC,ITYP) = CORG
         WARMES(IC,ITYP) = CMES
         NFCWAR(IC,ITYP) = NFCN
      RETURN
C
C             'SHO WARnings', ask if any suppressed mess in buffer
  200 CONTINUE
      IF (COPT .EQ. 'W') THEN
        ITYP = 1
        CTYP = 'WARNING'
      ELSE
        ITYP = 2
        CTYP = '*DEBUG*'
      ENDIF
      IF (NWRMES(ITYP) .GT. 0) THEN
         ENGLSH = ' WAS SUPPRESSED.  '
         IF (NWRMES(ITYP) .GT. 1) ENGLSH = 'S WERE SUPPRESSED.'
         WRITE (ISYSWR,'(/1X,I5,A,A,A,A/)') NWRMES(ITYP),
     +    ' MINUIT ',CTYP,' MESSAGE', ENGLSH
         NM = NWRMES(ITYP)
         IC = 0
         IF (NM .GT. MAXMES) THEN
              WRITE (ISYSWR,'(A,I2,A)')  ' ONLY THE MOST RECENT ',
     +          MAXMES,' WILL BE LISTED BELOW.'
              NM = MAXMES
              IC = ICIRC(ITYP)
         ENDIF
         WRITE (ISYSWR,'(A)') '  CALLS  ORIGIN         MESSAGE'
           DO 300 I= 1, NM
           IC = IC + 1
           IF (IC .GT. MAXMES)  IC = 1
           WRITE (ISYSWR,'(1X,I6,1X,A,1X,A)')
     +           NFCWAR(IC,ITYP),ORIGIN(IC,ITYP),WARMES(IC,ITYP)
 300       CONTINUE
         NWRMES(ITYP) = 0
         WRITE (ISYSWR,'(1H )')
      ENDIF
      RETURN
      END
