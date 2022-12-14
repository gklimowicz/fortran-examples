#include "rundeck_opts.h"
      MODULE WORKJK
      USE DIAG_COM, ONLY : JM,LM
      implicit none
      REAL*8, DIMENSION(JM,LM,2):: DPJK
      REAL*8, DIMENSION(2,LM,2) :: DPHEM
      REAL*8, DIMENSION(LM,2) :: DPGLOB

      END MODULE WORKJK


!------------------------------------------------

      BLOCK DATA BDWP
C****
C**** TITLES FOR SUBROUTINE DIAG7
C****
      use mdiag_com, only : lname_strlen,sname_strlen,units_strlen
      COMMON/D7COM/LNAME,SNAME,UNITS
      CHARACTER(len=lname_strlen) :: LNAME(12)
      CHARACTER(len=sname_strlen) :: SNAME(12)
      CHARACTER(len=units_strlen) ::  UNITS(12)
      DATA LNAME/
     1'WAVE POWER FOR U NEAR 850 MB AND EQUATOR  ',
     2'WAVE POWER FOR V NEAR 850 MB AND EQUATOR  ',
     3'WAVE POWER FOR U NEAR 300 MB AND EQUATOR  ',
     4'WAVE POWER FOR V NEAR 300 MB AND EQUATOR  ',
     5'WAVE POWER FOR U NEAR 50 MB AND EQUATOR   ',
     6'WAVE POWER FOR V NEAR 50 MB AND EQUATOR   ',
     7'WAVE POWER FOR PHI AT 922 MB AND 50 DEG N.',
     8'WAVE POWER FOR PHI AT 700 MB AND 50 DEG N.',
     9'WAVE POWER FOR PHI AT 500 MB AND 50 DEG N.',
     A'WAVE POWER FOR PHI AT 300 MB AND 50 DEG N.',
     B'WAVE POWER FOR PHI AT 100 MB AND 50 DEG N.',
     C'WAVE POWER FOR PHI AT 10 MB AND 50 DEG N. '/
!      .........1.........2.........3.........4.........5.........6
      DATA UNITS/
     1'DAY*(m/s)^2   ','DAY*(m/s)^2   ','10 DAY*(m/s)^2',
     4'DAY*(m/s)^2   ','10 DAY*(m/s)^2','DAY*(m/s)^2   ',
     7'10**3 DAY*m^2 ','10**3 DAY*m^2 ','10**3 DAY*m^2 ',
     A'10**3 DAY*m^2 ','10**4 DAY*m^2 ','10**4 DAY*m^2 '/
      DATA SNAME/
     1'WPU850EQU'   ,'WPV850EQU'   ,'WPU300EQU'   ,'WPV300EQU'   ,
     5'WPU50EQU'    ,'WPV50EQU'    ,'WPPHI922_50N','WPPHI700_50N',
     9'WPPHI500_50N','WPPHI300_50N','WPPHI100_50N','WPPHI10_50N' /

      END BLOCK DATA BDWP

!------------------------------------------------

!------------------------------------------------

      MODULE BDIJ
!@sum  stores information for outputting lon-lat diagnostics
      use DIAG_COM
      IMPLICIT NONE
      SAVE

!@var SENTDSE stand.eddy northw. transport of dry static energy * 16
!@var TENTDSE trans.eddy northw. transport of dry static energy * 16
      REAL*8, DIMENSION(IM,JM) :: SENTDSE,TENTDSE

      contains

      function mark (val,ibar,undef)
!@sum  mark selects a character (color) based on value and color bar
!@auth R. Ruedy
      real*8 val,undef
      integer ibar,n
      character*1 mark
      
      if (val .eq. undef) then
        mark=' '
      else
      select case (ibar)
      case (ib_pct)                                ! 0.....100 %
        n = 2.5 + val
        if (val .ge. 20.) n=23
        if (val .le.  0.) n= 1
        mark = cbar(ib_pct)(n:n)
      case (ib_pos)                                ! 0++++++++++
        n = 2.5 + val
c          non-unif scaling: (currently not used)
c          if (n .gt. 13) n = (n+123)/10
        if (n .gt. 38) n=38
        if (n .lt. 1 .or. val .le. 0.) n= 1
        mark = cbar(ib_pos)(n:n)
      case (ib_npp)                                ! ---0+++++++
        n = 11.5 + val
        if (n .gt. 38.) n=38
        if (n .lt.  1 ) n= 1
        mark = cbar(ib_npp)(n:n)
      case (ib_nnp)                                ! -------0+++
        n = 28.5 + val
        if (n .gt. 38.) n=38
        if (n .lt.  1 ) n= 1
        mark = cbar(ib_nnp)(n:n)
      case (ib_hyb)                      ! hybrid: multiple scales
        n = 2.5 + val
        if (n .gt. 28) n=(n+263)/10
        if (n .gt. 35) n=(n+180)/6
        if (n .gt. 37) n=37
        if (val .le.  0.) n=1
        mark = cbar(ib_hyb)(n:n)
      case (ib_ntr)                !tracers       ! ---0+++++++
        if (val.lt.0.) then
          n = 11.5-LOG(-val)/LOG(2.)
          if (n .le.  0) n= 1
          if (n .gt. 11) n=11
        else if (val.eq.0.) then
          n = 11
        else
          n = 11.5+LOG( val)/LOG(2.)
          if (n .lt. 11) n=11
          if (n .ge. 38) n=38
        end if
        mark = cbar(ib_npp)(n:n)                  ! use ib_npp
      end select
      end if

      return
      end function mark

      function ib_of_legnd (leg)
!@sum  ib_of_legnd finds the 'colorbar' for the given legend
!@auth R. Ruedy
      integer ib_of_legnd, leg

      ib_of_legnd = ib_pos
      if (legend(leg)(7:8) .eq. ',Z') ib_of_legnd = ib_nnp
      if (legend(leg)(7:8) .eq. ',9') ib_of_legnd = ib_npp
      if (index(legend(leg)(21:40),'-') .gt. 0) ib_of_legnd = ib_hyb
      if (index(legend(leg),'100 ') .gt. 0) ib_of_legnd = ib_pct
      if (legend(leg)(1:4) .eq. '9=-5') ib_of_legnd = ib_ntr

      return
      end function ib_of_legnd

      END MODULE BDIJ

!--------------------------------------------------------


!====================
!      SUBROUTINE DIAGKS
      module DIAGKS
C****
C**** THIS ROUTINE PRODUCES A SUMMARY OF KEY NUMBERS CALCULATED IN
C**** OTHER DIAGNOSTIC ROUTINES
C****
C**** CONTENTS OF KEYNR
C****
C****  N
C****
C****  1 MONTH
C****  2 TOTAL CLOUD COVER (PERCENT)
C****  3 SNOW AND ICE COVERAGE OF GLOBE (PERCENT)
C****  4 SNOW AND ICE COVERAGE OF NORTHERN HEMISPHERE (PERCENT)
C****  5 SNOW COVER--NORTHERN HEMSIPHERE (PERCENT)
C****  6 ICE COVER--NORTHERN HEMISPHERE (PERCENT)
C****  7 PLANETARY ALBEDO (PERCENT)
C****  8 SOLAR RADIATION ABSORBED BY ATMOSPHERE (WT/M**2)
C****  9 SOLAR RADIATION ABSORBED BY PLANET (WT/M**2)
C**** 10 NET HEAT AT GROUND (WT/M**2)
Cobso    ANGULAR MOMENTUM PER UNIT AREA (10**10 J*SEC/M**2)
Cobso    EVAPORATION (.1 MM/DAY)
C**** 11 PRECIPITATION (.1 MM/DAY)
C**** 12 SENSIBLE HEAT FLUX INTO GROUND (ABS.VALUE)
C**** 13 LATENT HEAT FLUX INTO GROUND (ABS.VALUE)
C**** 14 MEAN GROUND TEMPERATURE (DEGREES K)
C**** 15 MEAN GLOBAL ATMOSPHERIC TEMPERATURE (DEGREES K)
C**** 16 MERID. TEMPERATURE GRADIENT (N.HEMISPHERE)
C**** 17 MERID. TEMPERATURE GRADIENT (S.HEMISPHERE)
C**** 18 MEAN TROPOSPHERIC EKE-NORTHERN HEMISPHERE
C**** 19 MEAN TROPOSPHERIC EKE-SOUTHERN HEMISPHERE
C**** 20 MEAN TROPOSPHERIC ZKE-NORTHERN HEMISPHERE
C**** 21 MEAN TROPOSPHERIC ZKE-SOUTHERN HEMISPHERE
C**** 22 MEAN TROPOSPHERIC EPE-NORTHERN HEMISPHERE
C**** 23 MEAN TROPOSPHERIC ZPE-NORTHERN HEMISPHERE
C**** 24 MEAN EDDY KINETIC ENERGY AT EQUATOR
C**** 25 MAX. MEAN EDDY KINETIC ENERGY IN MID NORTH LATITUDES
C**** 26 MAX. ZONAL WIND (U COMPONENT) IN TROPOSPHERE (NH), M/SEC
C**** 27 LATITUDE CORRESPONDING TO 26
C**** 28 MAX. ZONAL WIND (U COMPONENT) IN TROPOSPHERE (SH), M/SEC
C**** 29 LATITUDE CORRESPONDING TO 28
C**** 30-32: 31 IS LARGEST VALUE OF STREAM FUNCTION, POSITIVE OR
C****    NEGATIVE; 30 AND 32 ARE THE MAGNITUDES OF THE LARGEST VALUES OF
C****    OPPOSITE SIGN TO THE NORTH AND SOUTH RESPECTIVELY
C**** 33-42 REFER TO NORTHERN HEMISPHERE ONLY
C**** 33 MAX.NORTHWARD TRANS. OF DRY STATIC ENERGY BY STANDING EDDIES
C**** 34 MAX.NORTHWARD TRANS. OF DRY STATIC ENERGY BY EDDIES
C**** 35 MAX. TOTAL NORTH. TRANS. OF DRY STATIC ENERGY
C**** 36 MAX.NORTHWARD TRANS. OF STATIC ENERGY BY EDDIES
C**** 37 MAX.TOTAL NORTH. TRANS. OF STATIC ENERGY
C**** 38 LATITUDE CORRESPONDING TO 37
C**** 39 MAX. NORTH. TRANS. OF ANGULAR MOMENTUM BY STANDING EDDIES
C**** 40 MAX. NORTH. TRANS. OF ANGULAR MOMENTUM BY EDDIES
C**** 41 MAX. TOTAL NORTH. TRANS. OF ANGULAR MOMENTUM
C**** 42 LATITUDE CORRESPONDING TO 41
C****
      USE CONSTANT, only : twopi
      use model_com, only: modelEclock
      USE MODEL_COM, only : JHOUR0,
     &     JDATE0,JMON0,AMON,AMON0,JYEAR0,
     &     Itime,ItimeI,Itime0,XLABEL,AMONTH,nday
      USE GEOM, only : DLAT,DXYP,LAT_DG
      USE ATM_COM, only : pmidl00
      USE DIAG_COM, only : jm,lm,keyct,keynr,nkeynr
      USE GC_COM, only : ned,jeq
      USE Dictionary_mod
      IMPLICIT NONE
      PRIVATE
      SAVE

      public KEYDJ,KEYJKT,KEYJKJ,KEYJLS,KEYJKE,KEYJKN,KEYIJ
     &     ,KEYD4,DIAGKN

      !REAL*8, DIMENSION(JM) :: FLAT
      REAL*8, DIMENSION(JM,LM) :: FKEY
      !REAL*8, DIMENSION(JM) :: ASUM
      !REAL*8, DIMENSION(2) :: HSUM
      !INTEGER, DIMENSION(2*NED) :: IK

      INTEGER ::
     &     I,I35,I70,IEND,ISIGN,
     &     J,J60,JMAX,JNDEX,JSTART,
     &     K,KEYMAX,KNDEX,
     &     LL,LMAX,LNLM,LNM,LSLM,LSM

      REAL*8 ::
     &     A,BIG,CEPT,CHECK,
     &     HN,HS,
     &     SAVE,DAYS,TEQ,TNOR,TSOU,
     &     UNLM,UNM,USLM,USM,X60

      contains
C****
C**** ENTRIES CALLED FROM DIAGJ
C****
!      ENTRY KEYDJ (N,FGLOB,FNH)
      subroutine KEYDJ(name,FGLOB,FNH)
      use mdiag_com, only : sname_strlen
      character(len=sname_strlen) :: name
      real*8 FGLOB,FNH

      SELECT CASE ( name )
      CASE ('J_totcld')           ; KEYNR( 2,keyct) = NINT(FGLOB)
      CASE ('J_snow_cover')       ; KEYNR( 5,keyct) = NINT(FNH)
      CASE ('J_ocn_lak_ice_frac') ; KEYNR( 6,keyct) = NINT(FNH)
      CASE ('J_plan_alb')         ; KEYNR( 7,keyct) = NINT(10.*fglob)
      CASE ('J_sw_abs_atm')       ; KEYNR( 8,keyct) = NINT(fglob)
      CASE ('J_net_rad_p0')       ; KEYNR( 9,keyct) = NINT(fglob)
      CASE ('J_nt_ht_z0')         ; KEYNR(10,keyct) = NINT(fglob)
      CASE ('J_prec')             ; KEYNR(11,keyct) = NINT(10.*fglob)
      CASE ('J_snsht_flx')        ; KEYNR(12,keyct) = NINT(-fglob)
      CASE ('J_evht_flx')         ; KEYNR(13,keyct) = NINT(-fglob)
      CASE ('J_tg1')              ; KEYNR(14,keyct) = NINT(.1*fglob)
!!!   CASE ('J_tair')             ; KEYNR(15,keyct) = NINT(.1*fglob)
      end select
      RETURN
      end subroutine KEYDJ
C****
C**** ENTRIES CALLED FROM DIAGJL VIA JLMAP OR FROM DIAGJK VIA JKMAP
C****
!      ENTRY KEYJKT (GSUM,ASUM)
      subroutine KEYJKT (GSUM,ASUM)
      real*8 GSUM
      REAL*8, DIMENSION(JM) :: ASUM
C**** TEMPERATURES
C      JEQ=2.+.5*(JM-1.)
      TEQ=.5*(ASUM(JEQ-1)+ASUM(JEQ))
      X60=TWOPI/(12.*DLAT)
      J60=.5+X60
      A=DXYP(J60+1)*(X60+.5-J60)
      TSOU=ASUM(J60+1)*A
      TNOR=ASUM(JM-J60)*A
      DO 210 J=1,J60
      A=A+DXYP(J)
      TSOU=TSOU+ASUM(J)*DXYP(J)
  210 TNOR=TNOR+ASUM(JM+1-J)*DXYP(J)
      KEYNR(16,KEYCT)=NINT(TEQ-TNOR/A)
      KEYNR(17,KEYCT)=NINT(TEQ-TSOU/A)
      KEYNR(15,KEYCT)=NINT(GSUM)

      RETURN
      end subroutine KEYJKT
C****
!      ENTRY KEYJKJ (L,FLAT)

      subroutine KEYJKJ (L,FLAT)
      integer L
      REAL*8, DIMENSION(JM) :: FLAT
C**** JET STREAMS
      IF (L.LT.LM) GO TO 220
      DO 216 LL=1,LM
      IF (pmidl00(ll).LT.200.) GO TO 218
  216 CONTINUE
  218 LMAX=LL-1
  220 IF (L.GT.LMAX) RETURN
      USLM=-999999.
      DO 222 J=3,JEQ
      IF (FLAT(J).LT.USLM) GO TO 222
      USLM=FLAT(J)
      JMAX=J
  222 CONTINUE
      CEPT=.5*(FLAT(JMAX-1)-FLAT(JMAX+1))/
     *  (FLAT(JMAX-1)-2.*FLAT(JMAX)+FLAT(JMAX+1))
      LSLM=INT((JMAX-1.5+CEPT)*DLAT*360/TWOPI+.5)-90
      UNLM=-999999.
      DO 224 J=JEQ,JM-1
      IF (FLAT(J).LT.UNLM) GO TO 224
      UNLM=FLAT(J)
      JMAX=J
  224 CONTINUE
      CEPT=.5*(FLAT(JMAX-1)-FLAT(JMAX+1))/
     *  (FLAT(JMAX-1)-2.*FLAT(JMAX)+FLAT(JMAX+1))
      LNLM=INT((JMAX-1.5+CEPT)*DLAT*360/TWOPI+.5)-90
      IF (L.LT.LMAX) GO TO 226
      USM=USLM
      LSM=LSLM
      UNM=UNLM
      LNM=LNLM
      RETURN
  226 IF (USLM.LT.USM) GO TO 228
      USM=USLM
      LSM=LSLM
  228 IF (UNLM.LT.UNM) GO TO 230
      UNM=UNLM
      LNM=LNLM
  230 IF (L.NE.1) RETURN
      KEYNR(26,KEYCT)=.1*UNM+.5
      KEYNR(27,KEYCT)=LNM
      KEYNR(28,KEYCT)=.1*USM+.5
      KEYNR(29,KEYCT)=-LSM
      RETURN
      end subroutine KEYJKJ
C****
!      ENTRY KEYJLS (L,FLAT)
      subroutine KEYJLS (L,FLAT)
      integer L
      REAL*8, DIMENSION(JM) :: FLAT
C**** STREAM FUNCTION
      DO 290 J=2,JM
  290 FKEY(J,L)=FLAT(J)
      IF (L.NE.1) RETURN
  300 SAVE=0.
      HS=0.
      HN=0.
      DO 310 K=1,LM
      DO 310 I=2,JM
      CHECK=ABS(FKEY(I,K))
      IF (CHECK.LT.SAVE) GO TO 310
      SAVE=CHECK
      JNDEX=I
      KNDEX=K
  310 CONTINUE
      SAVE=FKEY(JNDEX,KNDEX)
      ISIGN=1
      IF (SAVE.GT.0.0) ISIGN=-1
      IF (JNDEX.LT.4) GO TO 325
      IEND=JNDEX-1
      DO 320 K=1,LM
      DO 320 I=2,IEND
      CHECK=FKEY(I,K)*ISIGN
  320 IF (CHECK.GT.HS)HS=CHECK
  325 CONTINUE
      IF (JNDEX.GT.(JM-2))GO TO 335
      JSTART=JNDEX+1
      DO 330 K=1,LM
      DO 330 I=JSTART,JM
      CHECK=FKEY(I,K)*ISIGN
  330 IF (CHECK.GT.HN)HN=CHECK
  335 CONTINUE
      KEYNR(30,KEYCT)=ABS(HN)+0.5
      KEYNR(31,KEYCT)=NINT(SAVE)
      KEYNR(32,KEYCT)=ABS(HS)+0.5
      RETURN
      end subroutine KEYJLS
C****
!      ENTRY KEYJKE (NT,HSUM,ASUM)
      subroutine KEYJKE (NT,HSUM,ASUM)
      integer NT
      REAL*8, DIMENSION(2) :: HSUM
      REAL*8, DIMENSION(JM) :: ASUM
C**** EDDY AND ZONAL KINETIC ENERGY
      IF (NT.EQ.19) GO TO 450
      KEYNR(18,KEYCT)=NINT(HSUM(2))
      KEYNR(19,KEYCT)=NINT(HSUM(1))
      KEYNR(20,KEYCT)=KEYNR(20,KEYCT)-NINT(HSUM(2))
      KEYNR(21,KEYCT)=KEYNR(21,KEYCT)-NINT(HSUM(1))
      KEYNR(24,KEYCT)=NINT(ASUM(JEQ))

      BIG=-99999.
      I35=2.+(JM-1.)*125./180.
      I70=2.+(JM-1.)*160./180.
      DO 440 I=I35,I70
      IF (ASUM(I).LT.BIG) GO TO 440
      BIG=ASUM(I)
  440 CONTINUE
      KEYNR(25,KEYCT)=NINT(BIG)
      RETURN
  450 KEYNR(20,KEYCT)=KEYNR(20,KEYCT)+NINT(HSUM(2))
      KEYNR(21,KEYCT)=KEYNR(21,KEYCT)+NINT(HSUM(1))
      RETURN
      end subroutine KEYJKE
C****
!      ENTRY KEYJKN (NT,ASUM,SUMFAC)
      subroutine KEYJKN (NT,ASUM,SUMFAC)
      integer NT
      REAL*8, DIMENSION(JM) :: ASUM
      REAL*8 SUMFAC
C**** NORTHWARD TRANSPORTS : KEYNR(33:42)
  500 BIG=-99999.
      DO 510 I=JEQ+1,JM
      IF (ASUM(I).LT.BIG) GO TO 510
      BIG=ASUM(I)
      JNDEX=I
  510 CONTINUE
      BIG=BIG*SUMFAC
      KEYNR(NT,KEYCT)=NINT(BIG)
      if (nt==37 .or. nt==41) KEYNR(NT+1,KEYCT)=NINT(LAT_DG(JNDEX,2))
      RETURN
      end subroutine KEYJKN
C****
C**** ENTRY CALLED FROM DIAGIJ
C****
!      ENTRY KEYIJ(PISG,PISN)
      subroutine KEYIJ(PISG,PISN)
      REAL*8 PISG,PISN
      KEYNR(3,KEYCT)=NINT(PISG)
      KEYNR(4,KEYCT)=NINT(PISN)
      RETURN
      end subroutine KEYIJ
C****
C**** ENTRY CALLED FROM DIAG4
C****
!      ENTRY KEYD4 (IK)
      subroutine KEYD4 (IK)
      INTEGER, DIMENSION(2*NED) :: IK
      KEYNR(22,KEYCT)=(IK(10)+IK(20)+5)/10
      KEYNR(23,KEYCT)=(IK(8)+IK(18)+5)/10
      RETURN
      end subroutine KEYD4
C****
!      ENTRY DIAGKN
      subroutine DIAGKN
      use model_com, only: modelEclock
C**** PRINTS THE TABLE OF KEY NUMBERS
C****
      integer :: year, month, hour, date

      call modelEclock%get(year=year, month=month, hour=hour,
     *     date=date)
      DAYS=(Itime-Itime0)/FLOAT(nday)
      KEYNR(1,KEYCT)=JMON0
      IF (Itime.eq.ItimeI+1) KEYNR(1,KEYCT)=0
      IF (KEYCT.GE.2) THEN
        if (KEYNR(1,KEYCT-1).EQ.JMON0) KEYCT=KEYCT-1
      ENDIF
      WRITE(6,901) XLABEL
      WRITE(6,910) JYEAR0,AMON0,JDATE0,JHOUR0,
     *  YEAR,AMON,DATE,HOUR,ITIME,DAYS
      WRITE(6,902)
      DO 810 K=1,KEYCT
      IF (KEYNR(1,K).EQ.1) WRITE (6,905)
  810 WRITE(6,905) AMONTH(KEYNR(1,K)),(KEYNR(I,K),I=2,42)
      WRITE (6,915)
C**** additional key diag for Nino 3.4
      WRITE(6,906)
      DO K=1,KEYCT
        WRITE(6,907) AMONTH(KEYNR(1,K)),KEYNR(43,K)
      END DO
C****
      KEYCT=KEYCT+1
      KEYMAX=49
      IF (KEYNR(1,1).NE.0) KEYMAX=48
      IF (KEYCT.LE.KEYMAX) RETURN
C**** ROLL UP KEY NUMBERS 1 YEAR AT A TIME
      DO 820 K=1,36
      DO 820 I=1,NKEYNR
  820 KEYNR(I,K)=KEYNR(I,K+KEYMAX-36)
      DO 880 K=37,50
      DO 880 I=1,NKEYNR
  880 KEYNR(I,K)=0
      KEYCT=37
      RETURN
  901 FORMAT('1',A)
  902 FORMAT ('0',7X,'SN+IC NH NH AL AB NT NT PR        T   T-OF-ATM  EK
     *E   ZKE           EKE   JET-STREAMS STREAM-FN NOR-TRAN NOR-TRAN NO
     *RTH-TRANS'/
     *         5X,'CL GL    SN OI BE BY RD HT EC SN LAT OF  GL  GRAD ---
     *-- ----- EPE ZPE ------ NORTH SOUTH --------- DRY-STAT STAT-ENR ZO
     *N MOMENTM'/
     *      5X,'CV OB NH CV CV DO AT P0 Z0 IP HT  HT GD  OB NH SH NH ','
     *SH NH SH  NH  NH EQ  ML VL LT VL LT NH MAX SH SE ED TL ED TL LT SE
     * ED TL LT'/)
  905 FORMAT (1X,A3,4I3,I2,I4,5I3,I4,I3,I4,6I3,2I4,I3,I4,5I3,I4,11I3)
  906 FORMAT (1X,'Nino3.4*100')
  907 FORMAT (1X,A4,10I5)
  910 FORMAT ('0',15X,'From:',I6,A6,I2,',  Hr',I3,
     *  6X,'To:',I6,A6,I2,', Hr',I3,
     *  '  Model-Time:',I9,5X,'Dif:',F7.2,' Days')
  915 FORMAT('0')
      end subroutine DIAGKN

      END MODULE DIAGKS


!==================



      MODULE DIAG_SERIAL
      USE DIAG_COM, ONLY : IM, JM
      USE DOMAIN_DECOMP_ATM, only : GRID
      USE DOMAIN_DECOMP_1D, only : DIST_GRID, AM_I_ROOT
      USE DIAGKS

      PRIVATE
      PUBLIC :: JLMAP
      PUBLIC :: MAPTXT
      PUBLIC :: IJ_avg
      PUBLIC :: scale_ijlmap
      PUBLIC :: focean_glob
      PUBLIC :: diag_gather,diagj,diagjk,diagil,diag7p,diagij,diagcp,
     &          diag5p,diagdd,diagdh,diag4,diagkn,diag_scatter,
     &          diag_isccp

!ESMF: These globalsum routines are private to this module and execute
!      serially in a single processor.
      INTERFACE GLOBALSUM
        MODULE PROCEDURE GLOBALSUM_J
        MODULE PROCEDURE GLOBALSUM_JK
      END INTERFACE

      !REAL*8 :: FLAND_glob(IM,JM)
      !REAL*8 :: FEARTH_glob(IM,JM)
      REAL*8 :: FOCEAN_glob(IM,JM)
      REAL*8 :: FLICE_glob(IM,JM)
      REAL*8 :: ZATMO_glob(IM,JM)

      CONTAINS




      SUBROUTINE GLOBALSUM_J(grd_dum, garr, gsum,
     &                       hsum, istag, iskip, all)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum
      REAL*8,            INTENT(IN) :: garr(grd_dum%jm_world)
      REAL*8,            INTENT(OUT):: gsum
      REAL*8, OPTIONAL,  INTENT(OUT):: hsum(2)
      INTEGER,OPTIONAL,  INTENT(IN) :: istag
      INTEGER,OPTIONAL,  INTENT(IN) :: iskip
      LOGICAL,OPTIONAL,  INTENT(IN) :: all

      INTEGER :: IM, JM, J, ier
      LOGICAL :: istag_, iskip_


      IM   = grd_dum%IM_WORLD
      JM   = grd_dum%JM_WORLD

      istag_ = .false.
      If (Present(istag)) Then
        If (istag == 1) istag_ = .true.
      End If

      iskip_ = .false.
      If (Present(iskip)) Then
        If (iskip == 1) iskip_ = .true.
      End If

      If (istag_) then
        gsum = sum(garr(2:JM),1)
      ElseIf (iskip_) then
        gsum = sum(garr(2:JM-1),1)
      Else
        gsum = sum(garr(1:JM),1)
      EndIf
      If (Present(hsum)) then
        If (istag_) then
          hsum(1)   = Sum( garr(2     :JM/2),1   )
          hsum(2)   = Sum( garr(2+JM/2:JM  ),1   )
          hsum(1)   = hsum(1) + 0.5*garr(1+JM/2)
          hsum(2)   = hsum(2) + 0.5*garr(1+JM/2)
        Else
          hsum(1)   = Sum( garr(1     :JM/2),1   )
          hsum(2)   = Sum( garr(1+JM/2:JM  ),1   )
        EndIf
      EndIf

      END SUBROUTINE GLOBALSUM_J


      SUBROUTINE GLOBALSUM_JK(grd_dum, garr, gsum, hsum, istag, all)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum
      REAL*8,            INTENT(IN) :: garr(              : ,:)
!     REAL*8,            INTENT(IN) :: garr(grd_dum%jm_world,:)
      REAL*8,            INTENT(OUT):: gsum(size(garr,2))
      REAL*8, OPTIONAL,  INTENT(OUT):: hsum(2,size(garr,2))
      INTEGER,OPTIONAL,  INTENT(IN) :: istag
      LOGICAL,OPTIONAL,  INTENT(IN) :: all

      INTEGER :: k
      INTEGER :: ier
      INTEGER :: IM, JM
      LOGICAL :: istag_

      IM   = grd_dum%IM_WORLD
      JM   = grd_dum%JM_WORLD
      istag_ = .false.
      If (Present(istag)) Then
        If (istag == 1) istag_ = .true.
      End If

      If (istag_) then
        gsum = sum(garr(2:JM,:),1)
      Else
        gsum = sum(garr(1:JM,:),1)
      EndIf
      If (Present(hsum)) then
        If (istag_) then
          hsum(1,:)   = Sum( garr(2     :JM/2,:),1   )
          hsum(2,:)   = Sum( garr(2+JM/2:JM  ,:),1   )
          hsum(1,:)   = hsum(1,:) + 0.5*garr(1+JM/2,:)
          hsum(2,:)   = hsum(2,:) + 0.5*garr(1+JM/2,:)
        Else
          hsum(1,:)   = Sum( garr(1     :JM/2,:),1   )
          hsum(2,:)   = Sum( garr(1+JM/2:JM  ,:),1   )
        EndIf
      EndIf

      END SUBROUTINE GLOBALSUM_JK


!------------------------------------------------

      SUBROUTINE DIAGJ
!@sum DIAGJ produces area weighted statistics of zonal budget diags
!@+   based on settings and quantities found in j_defs
!@auth G. Schmidt/R. Ruedy/G. Russell
      USE CONSTANT, only : teeny
      USE DOMAIN_DECOMP_ATM, only : GRID
      use model_com, only: modelEclock
      USE MODEL_COM, only : idacc,jhour0,jdate0,amon,amon0,
     &     jyear0,itime,itime0,nday,xlabel,lrunid
      USE DIAG_COM, only : aj=>aj_out,areg=>areg_out,
     &     ntype_out,ntype,nreg,kaj,terrain,
     &     QDIAG,kdiag,namreg,ia_j,iden_j,iden_reg,
     *     scale_j,stitle_j,lname_j,name_j,units_j,
     *     fmt_j,fmt_reg
     &     ,jm=>jm_budg,dxyp_budg,lat_budg,Qbp
      USE MDIAG_COM, only : sname_strlen,units_strlen,lname_strlen
      USE MDIAG_COM, only : acc_period
      IMPLICIT NONE
      INTEGER, PARAMETER :: INC=1+(JM-1)/24
C**** Arrays needed for full output
      REAL*8, DIMENSION(JM+3,KAJ) :: BUDG
      CHARACTER*16, DIMENSION(KAJ) :: TITLEO
      CHARACTER(len=lname_strlen), DIMENSION(KAJ) :: LNAMEO
      CHARACTER(len=sname_strlen), DIMENSION(KAJ) :: SNAMEO
      CHARACTER(len=units_strlen), DIMENSION(KAJ) :: UNITSO

      REAL*8 :: BYIACC,FGLOB,GSUM,GWT,FHEM(2),HWT(2),DAYS

      REAL*8, DIMENSION(JM) :: FLAT,HSUMJ,HWTJ
      INTEGER, DIMENSION(JM) :: MLAT

      INTEGER :: IACC,J,JR,K,M,N,IT,iu_Ibp,n_out

      CHARACTER*200    :: fmt903
      CHARACTER*200    :: fmt918
      integer :: year, hour, date

      call modelEclock%get(year=year, hour=hour, date=date)

      fmt903 = "('0',131('-')/20X,'G      NH     SH   ',24I4)"
      fmt918 = "('0',16X,23(1X,A4)/17X,23(1X,A4)/1X,131('-'))"

C**** OPEN PLOTTABLE OUTPUT FILE IF DESIRED
      IF (QDIAG)
     &     call open_j(trim(acc_period)//'.j'//XLABEL(1:LRUNID)
     *     ,ntype_out,jm,lat_budg)

      DAYS=(Itime-Itime0)/FLOAT(nday)
C****
C**** LOOP OVER SURFACE TYPES: 1 TO NTYPE
C****
      IF (KDIAG(1).GT.7) GO TO 510

      DO M=1,NTYPE_OUT
        if(.not.Qbp(M)) cycle
        WRITE (6,901) XLABEL
        WRITE (6,902) TERRAIN(M),JYEAR0,AMON0,JDATE0,JHOUR0,
     *       YEAR,AMON,DATE,HOUR,ITIME,DAYS
#if (defined COMPILER_PGI)
        write(6,*) "skipping some info due to PGI bugs :-("
#else
        write(6,fmt=fmt903) NINT(LAT_BUDG(JM:INC:-INC))
#endif
        WRITE (6,905)
        n_out = 0
        DO N=1,kaj
          if(trim(stitle_j(n)).eq.'no output') cycle
          if(IDACC(IA_J(N)).eq.0) then
            ! may happen for ia_j == ia_12hr when starting
            ! at a time not divisible by 12 hrs
            write(6,*) 'skipping '//trim(NAME_J(N))//
     &           ' output due to insufficient accumulation steps'
            cycle
          endif
          n_out = n_out + 1
          FLAT(:)=AJ(:,N,M)*SCALE_J(N)/IDACC(IA_J(N))
          HSUMJ(:)=FLAT(:)*DXYP_BUDG(:)
          IF(IDEN_J(N).GT.0) THEN
            HWTJ(:)=AJ(:,IDEN_J(N),M)/IDACC(IA_J(IDEN_J(N)))
          ELSE
            HWTJ(:)=1d0
          ENDIF
          FLAT(:)=FLAT(:)/(HWTJ(:)+teeny)
          HWTJ(:)=HWTJ(:)*DXYP_BUDG(:)
          DO J=1,JM
            MLAT(J)=NINTlimit(FLAT(J) )
          END DO
          CALL GLOBALSUM(GRID, HSUMJ, GSUM, FHEM)
          CALL GLOBALSUM(GRID, HWTJ,  GWT,  HWT)
          FHEM(:)=FHEM(:)/(HWT(:)+teeny)
          FGLOB=GSUM/(GWT+teeny)
          IF (M.EQ.1) CALL KEYDJ (name_j(N),FGLOB,FHEM(2))
C**** Save BUDG for full output
          BUDG(1:JM,N)=FLAT(1:JM)
          BUDG(JM+1,N)=FHEM(1)
          BUDG(JM+2,N)=FHEM(2)
          BUDG(JM+3,N)=FGLOB
          TITLEO(N)=STITLE_J(N)
          LNAMEO(N)=LNAME_J(N)
          SNAMEO(N)=NAME_J(N)
          UNITSO(N)=UNITS_J(N)
          IF(INDEX(FMT_J(N),'24I').GT.0) THEN ! integer format
            WRITE (6,FMT_J(N)) STITLE_J(N),FGLOB,FHEM(2),FHEM(1),
     *           (MLAT(J),J=JM,INC,-INC)
          ELSE                                ! real format
            WRITE (6,FMT_J(N)) STITLE_J(N),FGLOB,FHEM(2),FHEM(1),
     *           (FLAT(J),J=JM,INC,-INC)
          ENDIF
        END DO ! end loop over quantities
#if (defined COMPILER_PGI)
        write(6,*) "skipping some info due to PGI bugs :-("
#else
        write(6,fmt=(fmt903)) NINT(LAT_BUDG(JM:INC:-INC))
#endif
        WRITE (6,905)
        IF (QDIAG) CALL POUT_J(TITLEO,SNAMEO,LNAMEO,UNITSO,BUDG,n_out,
     *       TERRAIN(M),M)
      END DO ! end loop over types
      if(qdiag) call close_j

      IF (.not.Qbp(ntype_out+1)) RETURN

  510 CONTINUE
C****
C**** PRODUCE REGIONAL STATISTICS
C****
      WRITE (6,901) XLABEL
      WRITE (6,902) '   (REGIONS)    ',
     *        JYEAR0,AMON0,JDATE0,JHOUR0,
     *        YEAR,AMON,DATE,HOUR,ITIME,DAYS
#if (defined COMPILER_PGI)
      write(6,*) "skipping some info due to PGI bugs :-("
#else
         write(6,fmt=fmt918) RESHAPE( (/NAMREG(1,1:23),NAMREG(2,1:23)/),
     *        (/23*2/) )
#endif
c     write(6,fmt=fmt918) NAMREG(1,1:23)

      DO N=1,kaj
        if(trim(stitle_j(n)).eq.'no output') cycle
        BYIACC=1./(IDACC(IA_J(N))+teeny)
        DO JR=1,23
          FLAT(JR)=AREG(JR,N)*SCALE_J(N)*BYIACC
        ENDDO
        IF(IDEN_REG(N).GT.0) THEN ! ratio
          DO JR=1,23
            FLAT(JR)=FLAT(JR)
     &           *IDACC(IA_J(IDEN_REG(N)))/(AREG(JR,IDEN_REG(N))+TEENY)
          ENDDO
        ENDIF
C**** select output format based on field name
        IF(TRIM(FMT_REG(N)).NE.'not computed') THEN
          IF(INDEX(FMT_REG(N),'23I').GT.0) THEN ! integer format
            WRITE (6,FMT_REG(N)) STITLE_J(N),(NINT(FLAT(JR)),JR=1,23)
          ELSE                                  ! real format
            WRITE (6,FMT_REG(N)) STITLE_J(N),(FLAT(JR),JR=1,23)
          ENDIF
        ENDIF
      END DO
#if (defined COMPILER_PGI)
      write(6,*) "skipping some info due to PGI bugs :-("
#else
      write(6,fmt=fmt918) RESHAPE( (/NAMREG(1,1:23),NAMREG(2,1:23)/),
     *                              (/23*2/) )
#endif
      WRITE (6,905)
      RETURN
C****
  901 FORMAT ('1',A)
  902 FORMAT ('0** BUDGETS ',A16,'**  From:',I6,A6,I2,',  Hr',I3,
     *  6X,'To:',I6,A6,I2,', Hr',I3,'  Model-Time:',I9,5X,
     *  'Dif:',F7.2,' Days')
  903 FORMAT ('0',131('-')/20X,'G      NH     SH   ',24I4)
  905 FORMAT (1X,131('-'))
  918 FORMAT ('0',16X,23(1X,A4)/17X,23(1X,A4)/1X,131('-'))
      END SUBROUTINE DIAGJ

      SUBROUTINE JKJL_TITLEX
!@sum  JKJL_TITLEX titles etc for composite jl, jk output
!@auth G. Schmidt/J. Lerner
      use filemanager
      USE DIAG_COM
      USE GC_COM, only : lname_gc,sname_gc,kagcx
      IMPLICIT NONE
      INTEGER :: K,kk,iu_Ijk
      LOGICAL qIjk,Ql(KAJL),Qk(KAGCx)
      character*80 line

      k = kajl

      inquire(file='Ijk',exist=qIjk)
      if(.not.qIjk) then
         call openunit('Ijk',iu_Ijk,.false.,.false.)
         write(iu_Ijk, FMT='(a)') 'list of JL-fields'
         do kk = 1,k
           write(iu_Ijk, '(i3,1x,a)') kk,lname_jl(kk)
         end do
      else if(kdiag(2).gt.0) then
         Ql=.false.
         call openunit('Ijk',iu_Ijk,.false.,.true.)
         read (iu_Ijk,'(a)',end=20) line
   10    read (iu_Ijk,'(a)',end=20) line
         if(line(1:1).eq.'l') go to 20
         read(line,'(i3)') kk
         Ql(kk)=.true.
         go to 10
   20    continue
         do kk=1,KAJL
           if(.not.Ql(kk)) sname_jl(kk)='skip'
         end do
      end if

      k = kagcx
      if(.not.qIjk) then
         write(iu_Ijk, FMT='(a)') 'list of GC-fields'
         do kk = 1,k
           write (iu_Ijk, '(i3,1x,a)') kk,lname_gc(kk)
         end do
         call closeunit(iu_Ijk)
      else if(kdiag(2).gt.0) then
         Qk=.false.
   30    read (iu_Ijk,'(a)',end=40) line
         read(line,'(i3)') kk
         Qk(kk)=.true.
         go to 30
   40    continue
         do kk=1,KAGCx
           if(.not.Qk(kk)) sname_gc(kk)='skip'
         end do
         call closeunit(iu_Ijk)
       end if

      RETURN
      END SUBROUTINE JKJL_TITLEX


      SUBROUTINE DIAGJK
      USE CONSTANT, only :
     &     grav,rgas,kapa,twopi,bygrav,tf,teeny,radius
      USE DOMAIN_DECOMP_ATM, only : GRID
      USE RESOLUTION, only : ls1=>ls1_nominal,pmtop
      USE MODEL_COM, only : xlabel,lrunid,DTsrc,IDACC
      USE ATM_COM, only : lm_req,req_fac_d,req_fac_m
      USE DYNAMICS, only : do_gwdrag
      USE GEOM, only :
     &     AREAG,BYDXYP,COSV,DXV,DXYP,DXYV,DYP,FCOR,WTJ,lat_dg
      USE MDIAG_COM, only : acc_period
      USE GC_COM, only : kep,
     &     agc=>agc_out,scale_gc,ia_gc,units_gc,sname_gc,lname_gc,
     &     jgrid_gc,pow_gc,denom_gc
      USE DIAG_COM, only : im,jm,lm,fim,byim,imh,
     &     kdiag,qdiag,linect,ia_dga,p1000k,
     &     plm,ple_dn,ple,pmb,kgz,
     &     aij,aijl,aijk,ijk_q, ij_zpmb1,
     &     ajl,scale_jl,ia_jl,units_jl,sname_jl,lname_jl,
     &     jgrid_jl,pow_jl,
     &     asjl,
     &     jk_tx,jk_hght,jk_q,jk_rh,jk_cldh2o,jk_cldwtr,jk_cldice,
     &     jl_dudfmdrg,jl_dumtndrg,jl_dushrdrg,
     &     jl_dumcdrgm10,jl_dumcdrgp10,jl_mcdrgpm10,
     &     jl_dumcdrgm20,jl_dumcdrgp20,jl_mcdrgpm20,
     &     jl_dumcdrgm40,jl_dumcdrgp40,jl_mcdrgpm40,
     &     jl_dudtsdif,jl_sumdrg,jl_dudtvdif,jl_dudtsdrg,jl_dtdtsdrg,
     &     jl_dtdyn,jl_mcmflx,jl_mcdflx,jl_rad_cool,
     &     jl_rhe,jl_trbke,jl_sshr,jl_trbhr,
     &     jl_trbdlht,jl_mcldht,jl_mcheat,jl_mcdeep,jl_mcshlw,jl_mcdry,
     &     jl_csizmc,jl_csizss,jl_damdc,jl_dammc,
     &     jl_uepac,jl_vepac,jl_wepac,jl_uwpac,jl_vwpac,jl_wwpac
     &    ,jl_cnumwm,jl_cnumim,jl_cnumws,jl_cnumis
      USE DIAG_COM_RAD
      USE GCDIAG
      USE WORKJK
      IMPLICIT NONE

      REAL*8, DIMENSION(JM+LM) :: ONES
      REAL*8, DIMENSION(JM,LM) :: AX,BX,CX,DX,EX
      REAL*8, DIMENSION(JM,LM) :: DENOM
      REAL*8, DIMENSION(JM,LM_REQ) :: ARQX
      REAL*8, DIMENSION(LM_REQ) :: BYPKS
      REAL*8, DIMENSION(0:IMH) :: AN,BN

      REAL*8, DIMENSION(JM,LM,2) :: DSJK
      REAL*8, DIMENSION(2,LM,2) :: DSHEM
      REAL*8, DIMENSION(LM,2) :: DSGLOB

Cbmp - ADDED
      REAL*8, DIMENSION(JM,LM) :: DPHJK
      REAL*8, DIMENSION(JM,LM) :: PIHJK
      REAL*8, DIMENSION(2,LM)  :: DHtemp
      REAL*8, DIMENSION(LM)    :: DGtemp
Cbmp - ADDED

      REAL*8, DIMENSION(LM) :: PM,PKM
      REAL*8, DIMENSION(JM,2) :: PJ
      REAL*8, DIMENSION(JM,kgz+1,4) :: AMPLTD,PHASE
      INTEGER, PARAMETER, DIMENSION(5) :: MW=(/1,2,3,6,9/)

      REAL*8, DIMENSION(JM,LM) :: ! outputs of EPFLXP
     &     DUDS,DMF,DEF,DMFR,DEFR,ER1,ER2

      INTEGER ::
     &     I,IP1,IX,J,J1,K,K1,KDN,KM,KUP,L,LR,M,N,ND

      REAL*8 ::
     &     BYDP2,BYDPK,BYFSQ,BYIADA,
     &     BYN,BYRCOS,DALPHA,DE4TI,
     &     DP,DPG,DPH,DPTI,DTHETA,EL4QI,ELOFIM,
     &     GBYRSQ,GSQ,PDN,PMK,PUP,
     &     PUTI,PVTI,QDP,SCALES,SCALET,SDDP,SKEI,
     &     SN,SNAMI,SNDEGI,SNELGI,SQM,SQN,SZNDEG,
     &     SZNELG,THETA,TX,UDXN,UDXS,UX

C**** avoid printing out diagnostics that are not yet defined
      if (IDACC(ia_dga).eq.0) RETURN

C**** OPEN PLOTTABLE OUTPUT FILE IF DESIRED
      IF(QDIAG) call open_jl(trim(acc_period)//'.jk'//XLABEL(1:LRUNID)
     *     ,jm,lm,lm_req,lat_dg)

C**** INITIALIZE CERTAIN QUANTITIES
      call JKJL_TITLEX

      KM=LM
      DO L=1,LM
        PKM(L)=PLM(L)**KAPA
        PM(L)=PLE(L)
      END DO
      BYPKS(1:LM_REQ)=1./(REQ_FAC_M(:)*PMTOP)**KAPA

      ONES=1.

      LINECT=65
C      WRITE (6,901)
      write(6,*) ' DEG K/DAY  = 0.01*SDAY*GRAV/SHA (= 8.445) W/(m^2*mb)'
      write(6,*) ' 10**18 JOULES = .864 * 10**30 GM*cm^2/s/DAY'
      BYIADA=1./(IDACC(ia_dga)+teeny)
      DO J=1,JM
        PJ(J,1)=0
        PJ(J,2)=0
        DO K=1,KM
          PJ(J,1)=PJ(J,1)+AGC(J,K,JK_DPA)
          PJ(J,2)=PJ(J,2)+AGC(J,K,JK_DPB)
        END DO
      END DO
C****
C**** INITIALIZE DELTA SIGMA IN PRESSURE COORDINATES
C****
C
C     J1=1 -> standard  Grid
C     J1=2 -> staggered Grid
C
      DO J1=1,2
         if (J1==1) then
           DO K=1,KM
            DO J=1,JM
               DPJK(J,K,J1) = AGC(J,K,J1)
               DSJK(J,K,J1)=AGC(J,K,J1)/(PJ(J,J1)+teeny)
               DPHJK(J,K) = AGC(J,K,J1)*WTJ(J,2,J1)
               PIHJK(J,K) = PJ(J,J1)*WTJ(J,2,J1)
            END DO
           END DO
           CALL GLOBALSUM(GRID, DPHJK, DGtemp, DHtemp)
           DPHEM(:,:,J1) = DHtemp
           DPGLOB(:,J1)  = DGtemp
           CALL GLOBALSUM(GRID, PIHJK, DGtemp, DHtemp)
           DSHEM(:,:,J1) =  DPHEM(:,:,J1)/(DHtemp+teeny)
           DSGLOB(:,J1)  =  DPGLOB(:,J1) /(DGtemp+teeny)
         else
           DO K=1,KM
            DO J=2,JM
               DPJK(J,K,J1) = AGC(J,K,J1)
               DSJK(J,K,J1)=AGC(J,K,J1)/(PJ(J,J1)+teeny)
               DPHJK(J,K) = AGC(J,K,J1)*WTJ(J,2,J1)
               PIHJK(J,K) = PJ(J,J1)*WTJ(J,2,J1)
            END DO
           END DO
           CALL GLOBALSUM(GRID, DPHJK, DGtemp, DHtemp, istag=1)
           DPHEM(:,:,J1) = DHtemp
           DPGLOB(:,J1)  = DGtemp
           CALL GLOBALSUM(GRID, PIHJK, DGtemp, DHtemp, istag=1)
           DSHEM(:,:,J1) =  DPHEM(:,:,J1)/(DHtemp+teeny)
           DSGLOB(:,J1)  =  DPGLOB(:,J1) /(DGtemp+teeny)
         endif
      END DO

cC****
cC**** Calculate a density field on tracer grid, edge pressure
cC****     (not quite ok if K-1 is underground?)
c      REAL*8, DIMENSION(JM,LM) :: RHO
c        PME(L)=PLE_DN(L)
c      DO K=1,KM-1
c      DO J=1,JM
c        IF (K.eq.1) RHO(J,K)=100.*PME(K)/(RGAS*(tf+
c     *        AGC(J,K  ,jk_temp)/(AGC(J,K  ,jk_dpa)+teeny)))
c        IF (K.gt.1) RHO(J,K)=100.*PME(K)/(RGAS*(tf+
c     *        AGC(J,K-1,jk_temp)/(AGC(J,K-1,jk_dpa)+teeny)+
c     *       (AGC(J,K  ,jk_temp)/(AGC(J,K  ,jk_dpa)+teeny)
c     *       -AGC(J,K-1,jk_temp)/(AGC(J,K-1,jk_dpa)+teeny))
c     *       *(PME(K)-PLM(K-1))/(PLM(K)-PLM(K-1))))
c      END DO
c      END DO

C****
C**** PROGNOSTIC QUANTITIES AT CONSTANT PRESSURE
C****
C**** # OF GRIDPOINTS, DELTA P, S.D. OF DELTA P
      n = JK_NPTSAVG
      SCALET = scale_gc(n)/idacc(ia_gc(n))
      CALL JLMAP(LNAME_GC(n),SNAME_GC(n),UNITS_GC(n),POW_GC(n),
     &     PLM,AGC(1,1,n),SCALET,ONES,ONES,LS1-1,1,JGRID_GC(n))
      n = JK_DPB
      SCALET = scale_gc(n)/idacc(ia_gc(n))
      CALL JLMAP(LNAME_GC(n),SNAME_GC(n),UNITS_GC(n),POW_GC(n),
     &     PLM,AGC(1,1,n),SCALET,ONES,ONES,LS1-1,2,JGRID_GC(n))
      DO 98 J=2,JM
      DO 98 K=1,LS1-1
      BYN=1./(AGC(J,K,JK_NPTSAVG)+1.D-10)
      AX(J,K)=0.
      SDDP=(AGC(J,K,JK_DPSQR)-AGC(J,K,JK_DPB)*AGC(J,K,JK_DPB)*BYN)*BYN
   98 IF (SDDP.GT.0.) AX(J,K)=SQRT(SDDP)
      n = jk_stdev_dp
      SCALET = scale_gc(n)
      CALL JLMAP(LNAME_gc(n),SNAME_gc(n),UNITS_GC(n),POW_GC(n),
     &     PLM,AX,SCALET,ONES,ONES,LS1-1,2,JGRID_GC(n))
C**** TEMPERATURE, HEIGHT, SPECIFIC AND RELATIVE HUMIDITY
      call print_generic_jl(jk_tx,sjl_index=1,vsumfac=10d0)
      call print_generic_jl(jk_hght,sjl_index=2)
      call print_generic_jl(jk_q)
      call print_generic_jl(jk_rh)

C**** CLOUD WATER CONTENT (total/liquid/ice)
      call print_generic_jl(jk_cldh2o)
      call print_generic_jl(jk_cldwtr)
      call print_generic_jl(jk_cldice)

C**** U AND V WINDS, STREAM FUNCTION
      n = JK_U
      SCALET = SCALE_GC(n)
      CALL JKMAP(LNAME_GC(n),SNAME_GC(n),UNITS_GC(n),POW_GC(n),
     &    PLM,AGC(1,1,n),SCALET,ONES,ONES,KM,2,JGRID_GC(n))
      n = JK_V
      SCALET = SCALE_GC(n)
      CALL JKMAP(LNAME_GC(n),SNAME_GC(n),UNITS_GC(n),POW_GC(n),
     &    PLM,AGC(1,1,n),SCALET,ONES,ONES,KM,2,JGRID_GC(n))
      n = jk_Vstar
      SCALET = SCALE_GC(n)
      CALL JKMAP(LNAME_GC(n),SNAME_GC(n),UNITS_GC(n),POW_GC(n),
     &    PLM,AGC(1,1,N),SCALET,ONES,ONES,KM,2,JGRID_GC(n))

      n = jk_psi
      SCALET = SCALE_GC(n)/IDACC(IA_GC(n))
      CALL JLMAP(LNAME_gc(n),SNAME_gc(n),UNITS_GC(n),POW_GC(n),
     &     PM,AGC(1,1,N),SCALET,ONES,ONES,KM,2,JGRID_GC(n))
      n = jk_psi_tem
      SCALET = SCALE_GC(n)/IDACC(IA_GC(n))
      CALL JLMAP(LNAME_gc(n),SNAME_gc(n),UNITS_GC(n),POW_GC(n),
     &     PM,AGC(1,1,N),SCALET,ONES,ONES,KM,2,JGRID_GC(n))
C**** RESIDUAL VERTICAL VELOCITY (W*)
      n = jk_Wstar
      SCALET=SCALE_GC(n)/IDACC(IA_GC(n))
      CALL JLMAP(LNAME_GC(n),SNAME_GC(n),UNITS_GC(n),POW_GC(n),
     &  PM,AGC(1,1,N),SCALET,ONES,ONES,KM-1,2,JGRID_GC(n))
C**** VERTICAL WINDS
      n = JK_VVEL
      SCALET = SCALE_GC(n)/IDACC(IA_GC(n))
      CALL JLMAP(LNAME_GC(n),SNAME_GC(n),UNITS_GC(n),POW_GC(n),
     &     PM,AGC(1,2,n),SCALET,ONES,ONES,KM-1,2,JGRID_GC(n))
C****
C**** CALCULATIONS FOR STANDING EDDIES
C****
        AX=0.
        BX=0.
        CX=0.
        EX=0.
      DO 170 J=2,JM
      DO 170 K=1,KM
      DPTI=0.
      PUTI=0.
      PVTI=0.
      DE4TI=0.
      EL4QI=0.
      SKEI=0.
      SNDEGI=0.
      SNELGI=0.
      SNAMI=0.
      I=IM
      DO IP1=1,IM
        IF (AIJK(I,J,K,IJK_DPB).GT.0.) THEN
          DPTI=DPTI+AIJK(I,J,K,IJK_DPB)
          BYDPK=1./(AIJK(I,J,K,IJK_DPB)+teeny)
          PUTI=PUTI+AIJK(I,J,K,IJK_UB)
          PVTI=PVTI+AIJK(I,J,K,IJK_VB)
          DE4TI=DE4TI+AIJK(I,J,K,IJK_DSE)
          QDP = .25*(AIJL(I,J-1,K,IJK_Q)+AIJL(IP1,J-1,K,IJK_Q)+
     &               AIJL(I,J  ,K,IJK_Q)+AIJL(IP1,J  ,K,IJK_Q))
          EL4QI=EL4QI+QDP
          SKEI=SKEI+(AIJK(I,J,K,IJK_UB)*AIJK(I,J,K,IJK_UB)
     *         +AIJK(I,J,K,IJK_VB)*AIJK(I,J,K,IJK_VB))*BYDPK
          SNDEGI=SNDEGI+(AIJK(I,J,K,IJK_DSE)*AIJK(I,J,K,IJK_VB)*BYDPK)
          SNELGI=SNELGI+(QDP*AIJK(I,J,K,IJK_VB)*BYDPK)
          SNAMI=SNAMI+AIJK(I,J,K,IJK_UB)*AIJK(I,J,K,IJK_VB)*BYDPK
        ENDIF
        I=IP1
      ENDDO
      AX(J,K)=SKEI-(PUTI*PUTI+PVTI*PVTI)/(DPTI+teeny)
      SZNDEG=DE4TI*PVTI/(DPTI+teeny)
      SZNELG=EL4QI*PVTI/(DPTI+teeny)
      EX(J,K)=(SNELGI-SZNELG)*DXV(J)
      BX(J,K)=(SNDEGI-SZNDEG)*DXV(J)
      CX(J,K)=SNAMI-PUTI*PVTI/(DPTI+teeny)
  170 CONTINUE
C**** STANDING EDDY, EDDY AND TOTAL KINETIC ENERGY
      n = jk_seke
      SCALET = scale_gc(n)
      CALL JKMAP(LNAME_gc(n),SNAME_gc(n),UNITS_GC(n),POW_GC(n),
     &     PLM,AX,SCALET,ONES,ONES,KM,2,JGRID_gc(n))
      n = jk_eddke
      SCALET = scale_gc(n)
      CALL JKMAP(LNAME_gc(n),SNAME_gc(n),UNITS_GC(n),POW_GC(n),
     &     PLM,AGC(1,1,n),SCALET,ONES,ONES,KM,2,JGRID_gc(n))
      n = jk_totke
      SCALET = SCALE_GC(n)
      CALL JKMAP(LNAME_GC(n),SNAME_GC(n),UNITS_GC(n),POW_GC(n),
     &    PLM,AGC(1,1,n),SCALET,ONES,ONES,KM,2,JGRID_GC(n))
C**** POTENTIAL TEMPERATURE, POTENTIAL VORTICITY
      DO 205 LR=1,LM_REQ
      DO 205 J=1,JM
  205 ARQX(J,LR)=ASJL(J,LR,1)*BYIADA+TF
      N = JK_THETA
      SCALET = SCALE_GC(n)
      SCALES = P1000K
      CALL JKMAP(LNAME_GC(N),SNAME_GC(N),UNITS_GC(N),POW_GC(n),
     &     PLM,AGC(1,1,N),SCALET,ONES,ONES,KM,2,JGRID_GC(N),
     &     ARQX,SCALES,ONES,BYPKS)
      N = JK_POTVORT
      SCALET = SCALE_GC(n)/IDACC(IA_GC(n))
      CALL JLMAP(LNAME_GC(n),SNAME_GC(n),UNITS_GC(n),POW_GC(n),
     &     PLM,AGC(1,1,N),SCALET,ONES,ONES,KM,2,JGRID_GC(n))
C****
C**** NORTHWARD TRANSPORTS AT CONSTANT PRESSURE
C****
C**** NORTHWARD TRANSPORT OF SENSIBLE HEAT BY EDDIES
      n = jk_eddntsh
      SCALET = scale_gc(n)
      CALL JKMAP(LNAME_gc(n),SNAME_gc(n),UNITS_GC(n),POW_GC(n),
     &    PLM,AGC(1,1,N),SCALET,ONES,ONES,KM,2,JGRID_GC(n))
C**** NORTHWARD TRANSPORT OF DRY STATIC ENERGY BY STANDING EDDIES,
C****   EDDIES, AND TOTAL
C**** Individual wave transports commented out. (gas - 05/2001)
      N = jk_nt_dse_se
      SCALET = scale_gc(n)
      CALL JKMAP(LNAME_gc(n),SNAME_gc(n),UNITS_GC(n),POW_GC(n),
     &     PLM,BX,SCALET,ONES,ONES,KM,2,JGRID_gc(n))
      n = jk_nt_dse_e
      SCALET = scale_gc(n)
      CALL JKMAP(LNAME_gc(n),SNAME_gc(n),UNITS_GC(n),POW_GC(n),
     &     PLM,AGC(1,1,N),SCALET,ONES,ONES,KM,2,JGRID_gc(n))
      n = jk_totntdse
      SCALET = scale_gc(n)
      CALL JKMAP(LNAME_gc(n),SNAME_gc(n),UNITS_GC(n),POW_GC(n),
     &     PLM,AGC(1,1,n),SCALET,ONES,ONES,KM,2,JGRID_gc(n))
C**** NORTHWARD TRANSPORT OF LATENT HEAT BY STAND.EDDY, EDDIES AND TOTAL
C**** New way! (Direct diag from QDYNAM, on layers)
      n = jl_nt_lh_e; nd = denom_gc(n)
      SCALET = scale_gc(n)/idacc(ia_gc(n))
      denom = agc(:,:,nd)/idacc(ia_gc(nd))+teeny
      dx(1:jm,1:lm) = agc(1:jm,1:lm,n)/denom
      CALL jlMAP(LNAME_gc(n),SNAME_gc(n),UNITS_gc(n),POW_gc(n),
     &     PLM,DX,SCALET,ONES,ONES,lm,2,JGRID_gc(n))
      n = jl_totntlh; nd = denom_gc(n)
      SCALET = SCALE_gc(n)/idacc(ia_gc(n))
      denom = agc(:,:,nd)/idacc(ia_gc(nd))+teeny
      dx(1:jm,1:lm) = agc(1:jm,1:lm,n)/denom
      CALL jlMAP(LNAME_gc(n),SNAME_gc(n),UNITS_gc(n),POW_gc(n),
     &     PLM,DX,SCALET,ONES,ONES,lm,2,JGRID_gc(n))
C**** NORTHWARD TRANSPORT OF LATENT HEAT BY STAND. EDDY, EDDIES AND TOTA
C**** Old Way!  (estimate from DIAGB using 2nd order advection on CP)
      n = jk_nt_lh_se
      SCALET = scale_gc(n)
      CALL JKMAP(LNAME_gc(n),SNAME_gc(n),UNITS_GC(n),POW_GC(n),
     &     PLM,EX,SCALET,ONES,ONES,KM,2,JGRID_gc(n))
      n = jk_eddntlh
      SCALET = scale_gc(n)
      CALL JKMAP(LNAME_gc(n),SNAME_gc(n),UNITS_GC(n),POW_GC(n),
     &     PLM,AGC(1,1,n),SCALET,ONES,ONES,KM,2,JGRID_gc(n))
      n = jk_totntlh
      SCALET = SCALE_GC(n)
      CALL JKMAP(LNAME_GC(n),SNAME_GC(n),UNITS_GC(n),POW_GC(n),
     &     PLM,AGC(1,1,n),SCALET,ONES,ONES,KM,2,JGRID_GC(n))
C**** NORTHWARD TRANSPORT OF STATIC ENERGY BY EDDIES AND TOTAL
      n = jk_nt_see
      SCALET = scale_gc(n)
      CALL JKMAP(LNAME_gc(n),SNAME_gc(n),UNITS_GC(n),POW_GC(n),
     &     PLM,AGC(1,1,N),SCALET,ONES,ONES,KM,2,JGRID_GC(n))
      n = jk_tot_nt_se
      SCALET = scale_gc(n)
      CALL JKMAP(LNAME_gc(n),SNAME_gc(n),UNITS_GC(n),POW_GC(n),
     &     PLM,AGC(1,1,N),SCALET,ONES,ONES,KM,2,JGRID_GC(n))
C**** NORTHWARD TRANSPORT OF KINETIC ENERGY
      n = jk_totntke
      SCALET = SCALE_GC(n)
      CALL JKMAP(LNAME_GC(n),SNAME_GC(n),UNITS_GC(n),POW_GC(n),
     &     PLM,AGC(1,1,n),SCALET,ONES,ONES,KM,2,JGRID_GC(n))
C**** NOR. TRANS. OF MOM, BY STANDING EDDIES, EDDIES, AND TOTAL ANG. MOM
      n = jk_nt_am_stand_eddy
      SCALET = scale_gc(n)
      CALL JKMAP(LNAME_gc(n),SNAME_gc(n),UNITS_GC(n),POW_GC(n),
     &     PLM,CX,SCALET,ONES,ONES,KM,2,JGRID_GC(n))
      n = jk_eddntmom
      SCALET = scale_gc(n)
      CALL JKMAP(LNAME_gc(n),SNAME_gc(n),UNITS_GC(n),POW_GC(n),
     &     PLM,AGC(1,1,n),SCALET,ONES,ONES,KM,2,JGRID_GC(n))
      n = jk_totntmom
      SCALET = scale_gc(n)
      CALL JKMAP(LNAME_gc(n),SNAME_gc(n),UNITS_GC(n),POW_GC(n),
     &     PLM,AGC(1,1,N),SCALET,ONES,ONES,KM,2,JGRID_GC(n))
C****
C**** DYNAMIC CONVERGENCE OF ENERGY
C****
      n = jk_dyn_conv_dse
      SCALET = scale_gc(n)
      CALL JKMAP(LNAME_gc(n),SNAME_gc(n),UNITS_GC(n),POW_GC(n),
     &     PLM,AGC(1,1,N),SCALET,ONES,ONES,KM,2,JGRID_gc(n))
      n = jk_dyn_conv_eddy_geop
      SCALET = scale_gc(n)
      CALL JKMAP(LNAME_gc(n),SNAME_gc(n),UNITS_GC(n),POW_GC(n),
     &     PLM,AGC(1,1,N),SCALET,ONES,ONES,KM,2,JGRID_gc(n))
C**** BAROCLINIC EKE GENERATION, P-K BY EDDY PRESSURE GRADIENT FORCE
      n = jk_barekegen
      SCALET = SCALE_GC(n)
      CALL JKMAP(LNAME_GC(n),SNAME_GC(n),UNITS_GC(n),POW_GC(n),
     &     PLM,AGC(1,1,n),SCALET,ONES,ONES,KM,2,JGRID_GC(n))
      n = jk_p2kedpgf
      SCALET = SCALE_GC(n)
      CALL JKMAP(LNAME_GC(n),SNAME_GC(n),UNITS_GC(n),POW_GC(n),
     &    PLM,AGC(1,1,n),SCALET,ONES,ONES,KM,2,JGRID_GC(n))
C****
C**** VERTICAL TRANSPORTS
C****
C**** VERTICAL TRANSPORT OF GEOPOTENTIAL ENERGY BY EDDIES
      n = jk_vtgeoeddy
      SCALET = SCALE_GC(n)/IDACC(IA_GC(n))
      CALL JLMAP(LNAME_GC(n),SNAME_GC(n),UNITS_GC(n),POW_GC(n),
     &     PM,AGC(1,1,n),SCALET,ONES,ONES,KM-1,2,JGRID_GC(n))
C**** VERTICAL TRANSPORT OF DRY STATIC ENERGY BY EDDIES AND TOTAL
      n = jk_eddvtdse
      SCALET = SCALE_GC(n)/IDACC(IA_GC(n))
      CALL JLMAP(LNAME_gc(n),SNAME_gc(n),UNITS_GC(n),POW_GC(n),
     &     PM,AGC(1,1,n),SCALET,ONES,ONES,KM-1,2,JGRID_gc(n))
      n = jk_totvtdse
      SCALET = SCALE_GC(n)/IDACC(IA_GC(n))
      CALL JLMAP(LNAME_GC(n),SNAME_GC(n),UNITS_GC(n),POW_GC(n),
     &     PM,AGC(1,1,n),SCALET,ONES,ONES,KM-1,2,JGRID_GC(n))
C**** VERTICAL TRANSPORT OF LATENT HEAT BY EDDIES AND TOTAL
C**** New way!
      n = jl_vt_lh_e
      SCALET = scale_gc(n)/idacc(ia_gc(n))
      CALL jlMAP(LNAME_gc(n),SNAME_gc(n),UNITS_gc(n),POW_gc(n),
     &     PM,AGC(1,1,N),SCALET,ONES,ONES,lm-1,2,JGRID_gc(n))
      n = jl_totvtlh
      SCALET = SCALE_gc(n)/idacc(ia_gc(n))
      CALL jlMAP(LNAME_gc(n),SNAME_gc(n),UNITS_gc(n),POW_gc(n),
     &     PM,Agc(1,1,n),SCALET,ONES,ONES,lm-1,2,JGRID_gc(n))
C**** VERTICAL TRANSPORT OF LATENT HEAT BY EDDIES AND TOTAL
C**** Old way!
      n = jk_eddvtlh
      SCALET = SCALE_GC(n)/IDACC(IA_GC(n))
      CALL JLMAP(LNAME_gc(n),SNAME_gc(n),UNITS_GC(n),POW_GC(n),
     &     PM,AGC(1,1,n),SCALET,ONES,ONES,KM-1,2,JGRID_gc(n))
      n = jk_totvtlh
      SCALET = SCALE_GC(n)/IDACC(IA_GC(n))
      CALL JLMAP(LNAME_GC(n),SNAME_GC(n),UNITS_GC(n),POW_GC(n),
     &     PM,AGC(1,1,n),SCALET,ONES,ONES,KM-1,2,JGRID_GC(n))
C**** VERTICAL TRANSPORT OF STATIC ENERGY BY EDDIES AND TOTAL
      n = jk_vt_se_eddy
      SCALET = SCALE_GC(n)/IDACC(IA_GC(n))
      CALL JLMAP(LNAME_gc(n),SNAME_gc(n),UNITS_GC(n),POW_GC(n),
     &     PM,AGC(1,1,N),SCALET,ONES,ONES,KM-1,2,JGRID_gc(n))
      n = jk_tot_vt_se
      SCALET = SCALE_GC(n)/IDACC(IA_GC(n))
      CALL JLMAP(LNAME_gc(n),SNAME_gc(n),UNITS_GC(n),POW_GC(n),
     &     PM,AGC(1,1,N),SCALET,ONES,ONES,KM-1,2,JGRID_gc(n))
C**** VERTICAL TRANSPORT OF KINETIC ENERGY
      n = jk_totvtke
      SCALET = SCALE_GC(n)/IDACC(IA_GC(n))
      CALL JLMAP(LNAME_GC(n),SNAME_GC(n),UNITS_GC(n),POW_GC(n),
     &     PM,AGC(1,1,n),SCALET,ONES,ONES,KM-1,2,JGRID_GC(n))
C**** VERTICAL TRANSPORT OF ANGULAR MOMENTUM BY LARGE SCALE MOTIONS
      n = jk_vtameddy; nd = denom_gc(n)
      SCALET = SCALE_GC(n)/IDACC(IA_GC(n))
      denom = agc(:,:,nd)/idacc(ia_gc(nd))+teeny
      AX = 0.
      DO K=1,KM-1
      DO J=1,JM
        AX(J,K)=AGC(J,K,n)/DENOM(J,K) !/RHO(J,K)
      END DO
      END DO
      CALL JLMAP(LNAME_GC(n),SNAME_GC(n),UNITS_GC(n),POW_GC(n),
     &     PM,AX,SCALET,ONES,ONES,KM-1,2,JGRID_GC(n))
      n = jk_totvtam
      SCALET = SCALE_GC(n)/IDACC(IA_GC(n))
      DO K=1,KM-1
      DO J=1,JM
        AX(J,K)=AGC(J,K,n)/DENOM(J,K) !/RHO(J,K)
      END DO
      END DO
      CALL JLMAP(LNAME_GC(n),SNAME_GC(n),UNITS_GC(n),POW_GC(n),
     &     PM,AX,SCALET,ONES,ONES,KM-1,2,JGRID_GC(n))
C**** VERTICAL TRANSPORT OF POTENTIAL VORTICITY TOTAL AND BY EDDIES
      n = jk_vtpv
      SCALET = SCALE_GC(n)/IDACC(IA_GC(n))
      CALL JLMAP(LNAME_GC(n),SNAME_GC(n),UNITS_GC(n),POW_GC(n),
     &     PM,AGC(1,1,n),SCALET,ONES,ONES,KM-1,2,JGRID_GC(n))
      n = jk_vtpveddy
      SCALET = SCALE_GC(n)/IDACC(IA_GC(n))
      CALL JLMAP(LNAME_GC(n),SNAME_GC(n),UNITS_GC(n),POW_GC(n),
     &     PM,AGC(1,1,N),SCALET,ONES,ONES,KM-1,2,JGRID_GC(n))
C**** NOR. TRANSPORT OF QUASI-GEOSTROPHIC POT. VORTICITY BY EDDIES
      n = jk_nt_eqgpv
      SCALET = scale_gc(n)
      CALL JKMAP(LNAME_gc(n),SNAME_gc(n),UNITS_GC(n),POW_GC(n),
     &     PLM,AGC(1,1,N),SCALET,ONES,ONES,KM,2,JGRID_gc(n))
C****
C**** Wave Energy (ELIASSEN PALM) FLUX:  NORTHWARD, VERTICAL, DIVERGENCE
C****
c      DO 510 K=1,KM
c      AX(1,K)=0.
c      DO 510 J=2,JM
c      UX=AGC(J,K,JK_U)/(AGC(J,K,JK_DPB)+teeny)
c      IF (ABS(UX).GE.teeny) GO TO 510
c      SN=+1.
c      IF (UX.LT.0.) SN=-1.
c      UX=SN*teeny
c  510 AX(J,K)=(AGC(J,K,JK_EDDNTGEO))/UX!*DXV(J)
c      n = jk_we_flx_nor
c      SCALET = scale_gc(n)
c      CALL JKMAP(LNAME_gc(n),SNAME_gc(n),UNITS_GC(n),POW_GC(n),
c     &     PLM,AX,SCALET,ONES,ONES,KM,2,JGRID_gc(n))
c      DO 520 K=1,KM-1
c      BX(1,K)=0.
c      BX(JM,K)=0.
c      DO 520 J=1,JM
c      IF (J.NE.1.AND.J.NE.JM) GO TO 516  ! corrected 4-25-2000
c      IF (J.EQ.1) UX=.5*(AGC(J+1,K,JK_U)+AGC(J+1,K+1,JK_U))/
c     *     (AGC(J+1,K,JK_DPB)+AGC(J+1,K+1,JK_DPB)+teeny)
c      IF (J.EQ.JM) UX=.5*(AGC(J,K,JK_U)+AGC(J,K+1,JK_U))/
c     *     (AGC(J,K,JK_DPB)+AGC(J,K+1,JK_DPB)+teeny)
c      GO TO 518
c  516 UX=(AGC(J,K,JK_U)+AGC(J+1,K,JK_U)
c     &   +AGC(J,K+1,JK_U)+AGC(J+1,K+1,JK_U))/
c     *   (AGC(J,K,JK_DPB)+AGC(J+1,K,JK_DPB)+
c     &    AGC(J,K+1,JK_DPB)+AGC(J+1,K+1,JK_DPB)+teeny)
c  518 IF (ABS(UX).GE.teeny) GO TO 520
c      SN=+1.
c      IF (UX.LT.0.) SN=-1.
c      UX=SN*teeny
c  520 BX(J,K)=AGC(J,K,JK_VTGEOEDDY)/(UX*RHO(J,K))
c      n = jk_epflx_v
c      SCALET = SCALE_GC(n)/IDACC(IA_GC(n))
c      CALL JLMAP(LNAME_gc(n),SNAME_gc(n),UNITS_GC(n),POW_GC(n),
c     &     PM,BX,SCALET,ONES,ONES,KM-1,2,JGRID_gc(n))
c      DO 530 K=1,KM
c      CX(1,K)=0.
c      CX(JM,K)=0.
c      DO 530 J=2,JM-1
c  530 CX(J,K)=.25*(AX(J+1,K)-AX(J,K))
c      DO 540 K=1,KM-1
c      DO 540 J=1,JM
c      CX(J,K)=CX(J,K)-BX(J,K)
c  540 CX(J,K+1)=CX(J,K+1)+BX(J,K)
c      n = jk_we_flx_div
c      SCALET = scale_gc(n)
c      CALL JKMAP(LNAME_gc(n),SNAME_gc(n),UNITS_GC(n),POW_GC(n),
c     &     PLM,CX,SCALET,BYDXYP,ONES,KM,2,JGRID_gc(n))
C****
C**** D/DY OF Q-G POTENTIAL VORTICITY AND REFRACTION INDICES
C****
C**** PRELIMINARIES:  VERTICAL DERIVATIVES AND N**2
      GSQ=GRAV*GRAV
      GBYRSQ=GRAV*GRAV/(RGAS*RGAS)
      IF (AGC(2,KM,JK_DPA).LT.teeny) GO TO 670  ! diagB not called yet
      DO 600 J=1,JM
      K1=1
  580 IF (AGC(J,K1,JK_DPA).GT.teeny) GO TO 590
      AX(J,K1)=0.
      BX(J,K1)=0.
      DX(J,K1)=0.
      K1=K1+1
      GO TO 580
  590 KDN=K1
      PDN=PM(KDN)+.5*AGC(J,KDN,JK_DPA)/(AGC(J,KDN,JK_NPTSAVG1)+teeny)
      DO 600 K=K1,KM
      DP=AGC(J,K,JK_DPA)
      PMK=PM(K)+.5*AGC(J,K,JK_DPA)/(AGC(J,K,JK_NPTSAVG1)+teeny)
      KUP=K+1
      IF (K.EQ.KM) KUP=KM
      PUP=PM(KUP)+.5*AGC(J,KUP,JK_DPA)/(AGC(J,KUP,JK_NPTSAVG1)+teeny)
      DALPHA=(AGC(J,KUP,JK_TEMP)/(AGC(J,KUP,JK_DPA)+teeny)+TF)/PUP-
     *  (AGC(J,KDN,JK_TEMP)/(AGC(J,KDN,JK_DPA)+teeny)+TF)/PDN
      DTHETA=AGC(J,KUP,JK_THETA)/(AGC(J,KUP,JK_DPA)+teeny)-
     *  AGC(J,KDN,JK_THETA)/(AGC(J,KDN,JK_DPA)+teeny)
      THETA=AGC(J,K,JK_THETA)/(AGC(J,K,JK_DPA)+teeny)
      TX=AGC(J,K,JK_TEMP)/(AGC(J,K,JK_DPA)+teeny)+TF
      IF (ABS(DTHETA).GE.teeny) GO TO 595
      SN=+1.
      IF (DTHETA.LT.0.) SN=-1.
      DTHETA=SN*teeny
  595 DX(J,K)=DP*FCOR(J)*PMK*THETA*(PUP-PDN)/(TX*DTHETA*DXYP(J))
      AX(J,K)=DALPHA/(PUP-PDN-teeny)
C**** CALCULATE N**2 AT PRESSURE LATITUDES
      BX(J,K)=-DP*GSQ*PMK*DTHETA/(RGAS*TX*THETA*(PUP-PDN-teeny))
      KDN=K
  600 PDN=PMK
C**** CALCULATE  Q12 = (D(UDX) + F*DA)/DA
      DO 620 K=1,KM
      UDXS=0.
      DO 610 J=1,JM-1
      UDXN=AGC(J+1,K,JK_U)/(AGC(J+1,K,JK_DPB)+teeny)*DXV(J+1)
      CX(J,K)=(UDXS-UDXN+FCOR(J))/DXYP(J)
  610 UDXS=UDXN
      CX(JM,K)=(UDXS+FCOR(JM))/DXYP(JM)
C**** FIND DQ/DY = (Q12(J)-Q12(J-1)+Q3(J)-Q3(J-1))/DY
      DO 620 J=JM,2,-1
      DP=AGC(J,K,JK_DPB)
      AX(J,K)=DP*(CX(J,K)-CX(J-1,K) + (AX(J,K)-AX(J-1,K))*
     *  (DX(J,K)+DX(J-1,K))/
     &     (AGC(J,K,JK_DPA)+AGC(J-1,K,JK_DPA)+teeny))/DYP(3)
  620 CONTINUE
      n = jk_del_qgpv
      SCALET = scale_gc(n)
      CALL JKMAP(LNAME_gc(n),SNAME_gc(n),UNITS_GC(n),POW_GC(n),
     &     PLM,AX,SCALET,ONES,ONES,KM,2,JGRID_gc(n))
C**** TERMS FOR THE REFRACTION INDEX EXPRESSION
      DO 640 J=2,JM
      BYFSQ=2.*DXYV(J)*DXYV(J)/(FCOR(J-1)*FCOR(J-1)+FCOR(J)*FCOR(J))
      DO 640 K=1,KM
      BYDP2=1./(AGC(J-1,K,JK_DPA)+AGC(J,K,JK_DPA)+teeny)
      TX=BYDP2*(AGC(J-1,K,JK_TEMP)+AGC(J,K,JK_TEMP))+TF
      DX(J,K)=GBYRSQ/(TX*TX)
      SQN=BYDP2*(BX(J-1,K)+BX(J,K))
      CX(J,K)=SQN*BYFSQ
      UX=AGC(J,K,JK_U)
      IF (ABS(UX).GE.teeny) GO TO 635
      SN=+1.
      IF (UX.LT.0.) SN=-1.
      UX=SN*teeny
  635 AX(J,K)=AX(J,K)/UX
  640 CONTINUE
C**** COMBINE TERMS, PRINT OUT REFRACTION INDICES
      SCALET = 1.D8
      IX = jk_refr_ind_wave1-1
      DO 660 M=1,5
      SQM=MW(M)*MW(M)
      DO 650 J=2,JM
      BYRCOS=1./(RADIUS*RADIUS*COSV(J)*COSV(J))
      DO 650 K=1,KM
      DP=AGC(J,K,JK_DPB)
  650 BX(J,K)=DP*(CX(J,K)*(AX(J,K)-SQM*BYRCOS)-.25*DX(J,K))
  660 CALL JKMAP(LNAME_gc(M+IX),SNAME_gc(M+IX),UNITS_GC(M+IX),
     &     POW_GC(M+IX),
     &    PLM,BX,SCALET,ONES,ONES,KM,2,JGRID_gc(1+IX))
  670 CONTINUE
C**** SKIP REMAINING MAPS IF DATA NOT AVAILABLE
      IF (AGC(1,1,JK_EPFLXNCP).NE.0.) GO TO 799
C****
C**** CHANGE OF THE MEAN FIELDS OF WIND AND TEMPERATURE
C****
C**** WIND: RATE OF CHANGE, ADVECTION, EDDY CONVERGENCE
      IF (IDACC(ia_dga).GT.1) then
        n = JK_TOTDUDT
        SCALET = SCALE_GC(n)/IDACC(IA_GC(n))
        CALL JLMAP(LNAME_GC(n),SNAME_GC(n),UNITS_GC(n),POW_GC(n),
     &       PLM,AGC(1,1,n),SCALET,ONES,ONES,KM,2,JGRID_GC(n))
      ENDIF
C**** Depending on whether EP fluxes have been specially calculated
C**** output full or approximate version
      IF (KEP.gt.0) THEN
        if (kdiag(2).ne.7)
     &       CALL EPFLXP(.true.,DUDS,DMF,DEF,DMFR,DEFR,ER1,ER2)
      ELSE ! these are not very good
        n = jk_dudtmadv
        SCALET = SCALE_GC(n)/IDACC(IA_GC(n))
        CALL JLMAP(LNAME_GC(n),SNAME_GC(n),UNITS_GC(n),POW_GC(n),
     &       PLM,AGC(1,1,n),SCALET,ONES,ONES,KM,2,JGRID_GC(n))
        n = jk_dudt_econv; nd = denom_gc(n)
        SCALET = scale_gc(n)/IDACC(IA_GC(n))
        denom = agc(:,:,nd)/idacc(ia_gc(nd))+teeny
        ax = agc(:,:,n)/denom
        CALL JLMAP(LNAME_gc(n),SNAME_gc(n),UNITS_GC(n),POW_GC(n),
     &       PLM,AX,SCALET,ONES,ONES,KM,2,jgrid_gc(n))
C**** WIND: TRANSFORMED ADVECTION, LAGRANGIAN CONVERGENCE (DEL.F)
        n = jk_dudttem
        SCALET = SCALE_GC(n)/IDACC(IA_GC(n))
        CALL JLMAP(LNAME_GC(n),SNAME_GC(n),UNITS_GC(n),POW_GC(n),
     &       PLM,AGC(1,1,n),SCALET,ONES,ONES,KM,2,JGRID_GC(n))
        n = jk_dudt_epdiv
        SCALET = SCALE_GC(n)/IDACC(IA_GC(n))
        CALL JLMAP(LNAME_gc(n),SNAME_gc(n),UNITS_GC(n),POW_GC(n),
     &       PLM,AGC(1,1,N),SCALET,ONES,ONES,KM,2,jgrid_gc(n))
      END IF

C**** WIND: DU/DT BY STRAT. DRAG -  MTN, DEFORM., SHEAR ...
      if (DO_GWDRAG) then
      SCALET = scale_jl(jl_dudfmdrg)/idacc(ia_jl(jl_dudfmdrg))
      n = jl_dumtndrg
      CALL JLMAP(LNAME_JL(n),SNAME_JL(n),UNITS_JL(n),POW_JL(n),
     *     PLM,AJL(1,1,n),SCALET,ONES,ONES,LM,2,JGRID_JL(n))
      n = jl_dudfmdrg
      CALL JLMAP(LNAME_JL(n),SNAME_JL(n),UNITS_JL(n),POW_JL(n),
     *     PLM,AJL(1,1,n),SCALET,ONES,ONES,LM,2,JGRID_JL(n))
      n = jl_dushrdrg
      CALL JLMAP(LNAME_JL(n),SNAME_JL(n),UNITS_JL(n),POW_JL(n),
     *     PLM,AJL(1,1,n),SCALET,ONES,ONES,LM,2,JGRID_JL(n))
      n = jl_mcdrgpm10
      CALL JLMAP(LNAME_jl(n),SNAME_jl(n),UNITS_JL(n),POW_JL(n),
     &     PLM,AJL(1,1,N),SCALET,ONES,ONES,LM,2,JGRID_jl(n))
      n = jl_mcdrgpm40
      CALL JLMAP(LNAME_jl(n),SNAME_jl(n),UNITS_JL(n),POW_JL(n),
     &     PLM,AJL(1,1,N),SCALET,ONES,ONES,LM,2,JGRID_jl(n))
      n = jl_mcdrgpm20
      CALL JLMAP(LNAME_jl(n),SNAME_jl(n),UNITS_JL(n),POW_JL(n),
     &     PLM,AJL(1,1,N),SCALET,ONES,ONES,LM,2,JGRID_jl(n))
C**** DU/DT BY STRAT. DRAG - TOTAL
      n = jl_sumdrg
      CALL JLMAP(LNAME_JL(n),SNAME_JL(n),UNITS_JL(n),POW_JL(n),
     *     PLM,AJL(1,1,N),SCALET,ONES,ONES,LM,2,JGRID_JL(n))
C**** CHANGE IN EAST WIND BY STRATOSPHERIC DIFFUSION
      n = jl_dudtsdif
      CALL JLMAP(LNAME_jl(n),SNAME_jl(n),UNITS_JL(n),POW_JL(n),
     *     PLM,AJL(1,1,n),SCALET,ONES,ONES,LM,2,JGRID_JL(n))
C**** CHANGE IN EAST WIND BY VERTICAL DIFFUSION
      n = jl_dudtvdif
      CALL JLMAP(LNAME_jl(n),SNAME_jl(n),UNITS_JL(n),POW_JL(n),
     *     PLM,AJL(1,1,n),SCALET,ONES,ONES,LM,2,JGRID_JL(n))
      end if  ! do_GWDRAG

C**** DU/DT BY SDRAG
      n = jl_dudtsdrg
      SCALET = scale_jl(n)/idacc(ia_jl(n))
      CALL JLMAP(LNAME_JL(n),SNAME_JL(n),UNITS_JL(n),POW_JL(n),
     &     PLM,AJL(1,1,n),SCALET,ONES,ONES,LM,2,JGRID_JL(n))
C**** TEMPERATURE: RATE OF CHANGE, ADVECTION, EDDY CONVERGENCE
      IF (IDACC(ia_dga).GT.1) then
        n = JK_TOTDTDT
        SCALET = SCALE_GC(n)/IDACC(IA_GC(n))
        CALL JLMAP(LNAME_GC(n),SNAME_GC(n),UNITS_GC(n),POW_GC(n),
     &       PLM,AGC(1,1,n),SCALET,ONES,ONES,KM,2,JGRID_GC(n))
      end if
      n = JK_DTDTMADV
      SCALET = SCALE_GC(n)/IDACC(IA_GC(n))
      CALL JLMAP(LNAME_GC(n),SNAME_GC(n),UNITS_GC(n),POW_GC(n),
     &     PLM,AGC(1,1,n),SCALET,ONES,ONES,KM,2,JGRID_GC(n))
      n = jk_dtempdt_econv; nd = denom_gc(n)
      SCALET = scale_gc(n)/idacc(ia_gc(n))
      denom = agc(:,:,nd)/idacc(ia_gc(nd))+teeny
      cx = agc(:,:,n)/denom
      CALL JLMAP(LNAME_gc(n),SNAME_gc(n),UNITS_GC(n),POW_GC(n),
     &      PLM,CX,SCALET,ONES,ONES,KM,2,JGRID_GC(n))
C**** TEMPERATURE: TRANSFORMED ADVECTION
      n = JK_DTDTTEM
      SCALET = SCALE_GC(n)/IDACC(IA_GC(n))
      CALL JLMAP(LNAME_GC(n),SNAME_GC(n),UNITS_GC(n),POW_GC(n),
     &     PLM,AGC(1,1,n),SCALET,ONES,ONES,KM,2,JGRID_GC(n))
C**** CHANGE IN TEMPERATURE BY STRATOSPHERIC DRAG
c      call print_generic_jl(jl_dtdtsdrg) ! PKM factor prevents use
      n = jl_dtdtsdrg
      SCALET = scale_jl(n)/idacc(ia_jl(n))
      CALL JLMAP(LNAME_JL(n),SNAME_JL(n),UNITS_JL(n),POW_JL(n),
     &     PLM,AJL(1,1,n),SCALET,ONES,PKM,KM,2,JGRID_JL(n))
C**** CHANGE IN TEMPERATURE BY DYNAMICS
      call print_generic_jl(jl_dtdyn)

  799 CONTINUE

C****
C**** Transplanted from DIAGJL
C****
      LINECT=65
C**** MASS FLUX MOIST CONVECTION
      call print_generic_jl(jl_mcmflx) ! should only go to lm-1
C**** DOWNDRAFT FLUX MOIST CONVECTION
      call print_generic_jl(jl_mcdflx) ! should only go to lm-1

C****
C**** RADIATION, CONDENSATION AND CONVECTION
C****
C**** SOLAR AND THERMAL RADIATION HEATING
      call print_generic_jl(jl_srhr,sjl_index=3)
      call print_generic_jl(jl_trcr,sjl_index=4)
      call print_generic_jl(jl_rad_cool,sjl_index=5)

C**** TOTAL, SUPER SATURATION, CONVECTIVE CLOUD COVER, EFFECTIVE RH
      call print_generic_jl(jl_totcld)
      call print_generic_jl(jl_sscld)
      call print_generic_jl(jl_mccld)
      call print_generic_jl(jl_rhe)

C**** WATER CLOUD COVER AND ICE CLOUD COVER
      call print_generic_jl(jl_wcld)
      call print_generic_jl(jl_icld)

C**** WATER AND ICE CLOUD  optical depth
      call print_generic_jl(jl_wcod)
      call print_generic_jl(jl_icod)

C**** Water and ice cloud particle sizes (weighted by opt depths)
      call print_generic_jl(jl_wcsiz)
      call print_generic_jl(jl_icsiz)

C**** TURBULENT KINETIC ENERGY
      call print_generic_jl(jl_trbke)

C**** HEATING BY LARGE SCALE COND., MOIST CONVECTION AND TURBULENCE
      call print_generic_jl(jl_sshr)
      call print_generic_jl(jl_trbhr)
      call print_generic_jl(jl_trbdlht)
      call print_generic_jl(jl_mcldht)
      call print_generic_jl(jl_mcheat)
      call print_generic_jl(jl_mcdeep)
      call print_generic_jl(jl_mcshlw)
      call print_generic_jl(jl_mcdry)

C**** Weighted average cloud sizes
      call print_generic_jl(jl_csizmc)
      call print_generic_jl(jl_csizss)

#ifdef CLD_AER_CDNC
C**** Weighted average warm cloud droplet number
      call print_generic_jl(jl_cnumwm)
      call print_generic_jl(jl_cnumws)

C**** Weighted average cold cloud droplet number
      call print_generic_jl(jl_cnumim)
      call print_generic_jl(jl_cnumis)

#endif
C**** this output is not required (very similar to jl_sscld etc.)
c       n=JL_CLDMC
c       SCALET = scale_jl(n)/idacc(ia_jl(n))
c       CALL JLMAP(LNAME_JL(n),SNAME_JL(n),UNITS_JL(n),POW_JL(n),
c      &     PLM,AJL(1,1,n),SCALET,ONES,ONES,LM,2,JGRID_JL(n))
c       n=JL_CLDSS
c       SCALET = scale_jl(n)/idacc(ia_jl(n))
c       CALL JLMAP(LNAME_JL(n),SNAME_JL(n),UNITS_JL(n),POW_JL(n),
c      &     PLM,AJL(1,1,n),SCALET,ONES,ONES,LM,2,JGRID_JL(n))

C****
C**** ENERGY
C****
C**** AVAILABLE POTENTIAL ENERGY
      n = JL_APE; nd = denom_gc(n)
      SCALET = scale_gc(n)/idacc(ia_gc(n))
      denom = agc(:,:,nd)/idacc(ia_gc(nd))+teeny
      ax(1:jm,1:lm) = agc(1:jm,1:lm,n)/denom
      CALL JLMAP(LNAME_GC(n),SNAME_GC(n),UNITS_GC(n),POW_GC(n),
     &      PLM,AX,SCALET,ONES,ONES,LM,2,JGRID_GC(n))
C****
C**** NORTHWARD TRANSPORTS
C****
C**** NOR. TRANSPORT OF QUASI-GEOSTROPHIC POT. VORTICITY BY EDDIES
c      DO 50 J=2,JM
c      DXCOSV(J)=DXV(J)*COSV(J)
c   50 CONTINUE
c      DO 366 L=1,LM
c      CX(1,L)=0.
c      CX(2,L)=DXCOSV(2)*(AGC(2,L,JL_TOTNTMOM)-AGC(2,L,JL_ZMFNTMOM))+
c     &     .25*FIM*FCOR(2)*COSP(2)*(AGC(2,L,JL_47)+AGC(3,L,JL_47))
c      DO 364 J=3,JM-1
c      DAM4=DXCOSV(J)*(AGC(J,L,JL_TOTNTMOM)-AGC(J,L,JL_ZMFNTMOM))
c      CX(J,L)=DAM4+.25*FIM*FCOR(J)*COSP(J)*
c     &     (AGC(J,L,JL_47)+AGC(J-1,L,JL_47))
c      CX(J-1,L)=CX(J-1,L)-DAM4
c  364 CONTINUE
c      CX(JM-1,L)=CX(JM-1,L)-DXCOSV(JM)*(AGC(JM,L,JL_TOTNTMOM)-
c     &     AGC(JM,L,JL_ZMFNTMOM))
c      CX(JM,L)=0.
c  366 CONTINUE
C****
C**** VERTICAL TRANSPORTS
C****
C**** VERTICAL TRANSPORT OF ANGULAR MOMENTUM BY SMALL SCALE MOTIONS
      call print_generic_jl(jl_damdc)
      call print_generic_jl(jl_dammc)

C****
C**** MERIDIONAL LUNES
C****
C**** U, V AND W VELOCITY FOR EAST PACIFIC
      call print_generic_jl(jl_uepac)
      call print_generic_jl(jl_vepac)
      call print_generic_jl(jl_wepac)

C**** U, V AND W VELOCITY FOR WEST PACIFIC
      call print_generic_jl(jl_uwpac)
      call print_generic_jl(jl_vwpac)
      call print_generic_jl(jl_wwpac)

C****
C**** ELIASSEN-PALM FLUX : NORTHWARD, VERTICAL, DIVERGENCE
C****
      n = JL_EPFLXN; nd = denom_gc(n)
      SCALET = scale_gc(n)/idacc(ia_gc(n))
      denom = agc(:,:,nd)/idacc(ia_gc(nd))+teeny
      ax(1:jm,1:lm) = agc(1:jm,1:lm,n)/denom
      CALL JLMAP(LNAME_GC(n),SNAME_GC(n),UNITS_GC(n),POW_GC(n),
     &     PLM,AX,SCALET,ONES,ONES,LM,2,JGRID_GC(n))
      n = JL_EPFLXV; nd = denom_gc(n)
      SCALET = scale_gc(n)/idacc(ia_gc(n))
      denom = agc(:,:,nd)/idacc(ia_gc(nd))+teeny
C**** scale with density for m^2/s^2 unit. Note that RHO is really a JK.
      AX = 0.
      DO L=1,LM-1
      DO J=1,JM
        AX(J,L)=AGC(J,L,n)/DENOM(J,L) !RHO(J,L)
      END DO
      END DO
      CALL JLMAP(LNAME_GC(n),SNAME_GC(n),UNITS_GC(n),POW_GC(n),
     &     PLE,AX,SCALET,ONES,ONES,LM-1,2,JGRID_GC(n))
!      n = JL_EPFLXVm2ps2  ! in m2/s2
!      SCALET=.125*RADIUS
!      CALL JLMAP(LNAME_JL(n),SNAME_JL(n),UNITS_JL(n),POW_JL(n),
!     &     PLE,AGC(1,1,JL_EPFLXV),SCALET,COSBYPDA,BYDSIG,LM-1,1,
!     &     JGRID_JL(n))
      n = jl_epflx_div; nd = denom_gc(n)
      denom = agc(:,:,nd)/idacc(ia_gc(nd))+teeny
      ax(1:jm,1:lm) = agc(1:jm,1:lm,n)/denom
      SCALET = scale_gc(n)/idacc(ia_gc(n))
      CALL JLMAP(LNAME_gc(n),SNAME_gc(n),UNITS_GC(n),POW_GC(n),
     &     PLM,AX,SCALET,ONES,ONES,LM,2,jgrid_gc(n))

C****
C**** FOURIER ANALYSIS OF GEOPOTENTIAL HEIGHTS FOR WAVE NUMBERS 1 TO 4,
C****   AMPLITUDE AND PHASE
C****
            LINECT=63
      ELOFIM=.5*TWOPI-TWOPI/FIM

      DO K=1,KGZ
      DO N=1,4
      AMPLTD(1,K,N)=0.
      AMPLTD(JM,K,N)=0.
      PHASE(1,K,N)=0.
      PHASE(JM,K,N)=0.
      ENDDO
      DO J=2,JM-1
      CALL FFT (AIJ(1,J,IJ_ZPMB1-1+K),AN,BN)
      DO N=1,4
      AMPLTD(J,K,N)=SQRT(AN(N)*AN(N)+BN(N)*BN(N))
      PHASE(J,K,N)=(ATAN2(BN(N),AN(N))-TWOPI)/N+ELOFIM
      IF (PHASE(J,K,N).LE.-.5*TWOPI) PHASE(J,K,N)=PHASE(J,K,N)+TWOPI
      PHASE(J,K,N)=-PHASE(J,K,N)
      ENDDO
      ENDDO
      ENDDO
      SCALET = BYIADA
      IX = jl_phi_amp_wave1-1
      DO N=1,4
      CALL JLMAP(LNAME_gc(N+ix),SNAME_gc(N+ix),UNITS_gc(N+ix),
     &        POW_gc(N+ix),
     &      PMB,AMPLTD(1,1,N),SCALET,ONES,ONES,kgz,2,jgrid_gc(N+ix))
      ENDDO
      SCALET = 360./TWOPI
      IX = jl_phi_phase_wave1-1
      DO N=1,4
      CALL JLMAP(LNAME_gc(N+ix),SNAME_gc(N+ix),UNITS_gc(N+ix),
     &        POW_gc(N+ix),
     &      PMB,PHASE(1,1,N),SCALET,ONES,ONES,kgz,2,jgrid_gc(N+ix))
      ENDDO

      if(qdiag) call close_jl

      RETURN
  901 FORMAT (
     *  ' DEG K/DAY  = 0.01*SDAY*GRAV/SHA (= 8.445) W/(m^2*mb)'/
     *  ' 10**18 JOULES = .864 * 10**30 GM*cm^2/s/DAY')
      END SUBROUTINE DIAGJK


      SUBROUTINE JKMAP(LNAME,SNAME,UNITS,POW10P,
     &     PM,AX,SCALET,SCALEJ,SCALEK,KMAX,JWT,J1,
     *  ARQX,SCALER,SCALJR,SCALLR)
      USE CONSTANT, only : teeny
      USE DOMAIN_DECOMP_ATM, only : GRID
      use model_com, only: modelEclock
      USE MODEL_COM, only :
     &     JDATE0,JMON0,AMON0,AMON,JYEAR0,XLABEL
      USE WORKJK
      USE GEOM, only :
     &     LAT_DG,WTJ
      USE MDIAG_COM, only : acc_period
     &     ,sname_strlen,units_strlen,lname_strlen
      USE DIAG_COM, only : jm,lm,
     &     QDIAG,lm_req,inc=>incj,linect
      IMPLICIT NONE

!@var units string containing output field units
      CHARACTER(LEN=units_strlen) :: UNITS,UNITS_WITH_SCALE
!@var lname string describing output field
      CHARACTER(LEN=lname_strlen) :: LNAME
!@var sname string referencing output field
      CHARACTER(LEN=sname_strlen) :: SNAME
!@var title string, formed as concatentation of lname//units
      CHARACTER(LEN=64) :: TITLE

      INTEGER, DIMENSION(JM) :: MLAT
      REAL*8, DIMENSION(JM) :: FLAT,ASUM
      REAL*8, DIMENSION(2) :: AHEM,AHEML

      INTEGER :: J1,JWT,KMAX
      REAL*8 :: SCALET,SCALER,PRTFAC
      INTEGER :: POW10P
      REAL*8, DIMENSION(JM,LM) :: AX
      REAL*8, DIMENSION(JM,LM_REQ) :: ARQX
      REAL*8, DIMENSION(JM) :: SCALEJ,SCALJR
      REAL*8, DIMENSION(LM) :: SCALEK
      REAL*8, DIMENSION(LM_REQ) :: SCALLR
      REAL*8, DIMENSION(LM+LM_REQ) :: PM

      REAL*8, DIMENSION(JM,LM) :: CX

      CHARACTER*4 DASH,WORD(4)
      DATA DASH/'----'/,WORD/'SUM','MEAN',' ','.1*'/

      INTEGER :: IWORD,J,JHEMI,K,L ,ksx,klmax
      REAL*8 :: AGLOB,FGLOB,FLATJ,G1,H1,H2,SUMFAC

      REAL*8, DIMENSION(JM+3,LM+LM_REQ+1) :: XJL ! for binary output
      CHARACTER XLB*16,CLAT*16,CPRES*16,CBLANK*16,TITLEO*80,TPOW*8
      DATA CLAT/'LATITUDE'/,CPRES/'PRESSURE (MB)'/,CBLANK/' '/

      optional :: ARQX,SCALER,SCALJR,SCALLR

      integer :: year, date

      call modelEclock%get(year=year, date=date)

      if ( present(ARQX) ) goto 777

      if(sname.eq.'skip') return
C form title string
      units_with_scale = units
      PRTFAC = 10.**(-pow10p)
      title = trim(lname)//' ('//trim(units)//')'
      if(pow10p.ne.0) then
         WRITE (tpow, '(i3)') pow10p  ! checked by BMP OK for parallel output
         tpow='10**'//trim(adjustl(tpow))
         units_with_scale=trim(tpow)//' '//trim(units_with_scale)
         title = trim(lname)//' ('//trim(units_with_scale)//')'
      endif
C****
C**** PRODUCE A LATITUDE BY LAYER TABLE OF THE ARRAY A
C****
   10 LINECT=LINECT+KMAX+7
      IF (LINECT.LE.60) GO TO 20
      WRITE (6,907)
     *  XLABEL(1:105),JDATE0,AMON0,JYEAR0,DATE,AMON,YEAR
      LINECT=KMAX+8
   20 WRITE (6,901) TITLE,(DASH,J=J1,JM,INC)
      WRITE (6,904) WORD(JWT),(NINT(LAT_DG(J,J1)),J=JM,J1,-INC)
      WRITE (6,905) (DASH,J=J1,JM,INC)
         DO 40 L=1,LM+LM_REQ+1
         DO 40 J=1,JM+3
   40    XJL(J,L) = -1.D30
         KSX = 0            ! KSX = LAYERS GENERATED AT ENTRY
      CX = 0.0
  100 DO 110 J=J1,JM
      DO 110 K=1,KMAX
  110 CX(J,K)=AX(J,K)*SCALET*SCALEJ(J)*SCALEK(K)
         KLMAX = KMAX+KSX
C**** HORIZONTAL SUMS AND TABLE ENTRIES
      AHEM(:) = 0.
      DO K=KMAX,1,-1
         If (J1==1) then  ! Standard Grid
            DO J=1,JM
               FLAT(J)  = CX(J,K)/(DPJK(J,K,J1)+teeny)
               XJL(J,K) = FLAT(J)*PRTFAC
               FLAT(J)  = FLAT(J)*PRTFAC
               IF (DPJK(J,K,J1).EQ.0.) XJL(J,K) = -1.D30
               MLAT(J)=NINT(MIN(1d5,MAX(-1d5,FLAT(J)))) ! prevent too large int?
            END DO
            CALL GLOBALSUM(GRID, CX(:,K)*WTJ(:,JWT,J1)*PRTFAC,
     *                     AGLOB, AHEML)
         Else             ! Staggered Grid
            DO J=2,JM
               FLAT(J)  = CX(J,K)/(DPJK(J,K,J1)+teeny)
               XJL(J,K) = FLAT(J)*PRTFAC
               FLAT(J)  = FLAT(J)*PRTFAC
               IF (DPJK(J,K,J1).EQ.0.) XJL(J,K) = -1.D30
               MLAT(J)=NINT(MIN(1d5,MAX(-1d5,FLAT(J)))) ! prevent too large int?
            END DO
            CALL GLOBALSUM(GRID, CX(:,K)*WTJ(:,JWT,J1)*PRTFAC,
     *                     AGLOB, AHEML, istag=1)
         EndIf
         AHEM(:) = AHEM(:) + AHEML(:)
         H1=AHEML(1)/(DPHEM(1,K,J1)+teeny)
         H2=AHEML(2)/(DPHEM(2,K,J1)+teeny)
         G1=AGLOB/(DPGLOB(K,J1)+teeny)
         XJL(JM+3,K)=H1   ! SOUTHERN HEM
         XJL(JM+2,K)=H2   ! NORTHERN HEM
         XJL(JM+1,K)=G1   ! GLOBAL
         WRITE (6,902) PM(K),G1,H2,H1,(MLAT(J),J=JM,J1,-INC)
         CALL KEYNRL (SNAME,K,FLAT)
      END DO

C**** VERTICAL SUMS
      WRITE (6,905) (DASH,J=J1,JM,INC)
      SUMFAC=1.
      IWORD=3
      IF ( SNAME.EQ.'temp' .OR. ! make sumfac an argument to avoid this
     &     SNAME.EQ.'v' .OR.
     &     SNAME.EQ.'tot_nt_dse' .OR.
     &     SNAME.EQ.'tot_nt_se' .OR.
     &     SNAME.EQ.'tot_nt_am') THEN
         SUMFAC=10.
         IWORD=4
      ENDIF
      DO 180 J=J1,JM
      ASUM(J)=0.
      DO 170 K=1,KMAX
  170 ASUM(J)=ASUM(J)+CX(J,K)*PRTFAC
      ASUM(J)=ASUM(J)/SUM(DPJK(J,:,J1))
         XJL(J,LM+LM_REQ+1)=ASUM(J)
  180 MLAT(J)=NINT(ASUM(J)*SUMFAC)

      aglob = 0.
      ahem(:) = ahem(:)*sumfac
      do jhemi=1,2
         aglob = aglob + ahem(jhemi)
         ahem(jhemi) = ahem(jhemi)/sum(dphem(jhemi,:,j1))
      enddo
      aglob = aglob/sum(dpglob(:,j1))
         XJL(JM+3,LM+LM_REQ+1)=AHEM(1)/SUMFAC   ! SOUTHERN HEM
         XJL(JM+2,LM+LM_REQ+1)=AHEM(2)/SUMFAC   ! NORTHERN HEM
         XJL(JM+1,LM+LM_REQ+1)=AGLOB/SUMFAC     ! GLOBAL
         XLB=' '//acc_period(1:3)//' '//acc_period(4:12)//'  '
         TITLEO=TITLE//XLB
         IF(QDIAG) CALL POUT_JL(TITLEO,LNAME,SNAME,UNITS_WITH_SCALE,
     *        J1,KLMAX,XJL,PM,CLAT,CPRES)
      WRITE (6,903) WORD(IWORD),AGLOB,AHEM(2),AHEM(1),
     *  (MLAT(J),J=JM,J1,-INC)
         CALL KEYVSUMS(SNAME,AGLOB,AHEM,ASUM,SUMFAC)
      RETURN
C****
!      ENTRY JKMAPS(LNAME,SNAME,UNITS,POW10P,
!     &     PM,AX,SCALET,SCALEJ,SCALEK,KMAX,JWT,J1,
!     *  ARQX,SCALER,SCALJR,SCALLR)
 777  continue
      if(sname.eq.'skip') return


C form title string
      units_with_scale = units
      title = trim(lname)//' ('//trim(units)//')'
      PRTFAC = 10.**(-pow10p)
      if(pow10p.ne.0) then
         write(tpow,'(i3)') pow10p
         tpow='10**'//trim(adjustl(tpow))
         units_with_scale=trim(tpow)//' '//trim(units_with_scale)
         title = trim(lname)//' ('//trim(units_with_scale)//')'
      endif
         KSX = 3
         DO 205 L=1,LM+LM_REQ+1
         DO 205 J=1,JM+3
  205    XJL(J,L) = -1.D30
      LINECT=LINECT+KMAX+10
      IF (LINECT.LE.60) GO TO 230
      WRITE (6,907)
     *  XLABEL(1:105),JDATE0,AMON0,JYEAR0,DATE,AMON,YEAR
      LINECT=KMAX+11
  230 CONTINUE
C**** PRODUCE UPPER STRATOSPHERE NUMBERS FIRST
      WRITE (6,901) TITLE,(DASH,J=J1,JM,INC)
      WRITE (6,904) WORD(JWT),(NINT(LAT_DG(J,J1)),J=JM,J1,-INC)
      WRITE (6,905) (DASH,J=J1,JM,INC)

      DO L=LM_REQ,1,-1
         If (J1==1) then   ! Standard Grid
            DO J=1,JM
               FLATJ=ARQX(J,L)*SCALER*SCALJR(J)*SCALLR(L)
               XJL(J,L+KMAX) = FLATJ
c              FLATJ=FLATJ*PRTFAC
C               MLAT(J)=NINT(FLATJ)
               MLAT(J)=NINT(MIN(1d5,MAX(-1d5,FLATJ)))
               FLAT(J) = FLATJ*WTJ(J,JWT,J1)
            END DO
            CALL GLOBALSUM(GRID, FLAT(:), FGLOB, AHEM)
            FGLOB=FGLOB/JWT
         Else              ! Staggered Grid
            DO J=2,JM
               FLATJ=ARQX(J,L)*SCALER*SCALJR(J)*SCALLR(L)
               XJL(J,L+KMAX) = FLATJ
c              FLATJ=FLATJ*PRTFAC
C               MLAT(J)=NINT(FLATJ)
               MLAT(J)=NINT(MIN(1d5,MAX(-1d5,FLATJ)))
               FLAT(J) = FLATJ*WTJ(J,JWT,J1)
            END DO
            CALL GLOBALSUM(GRID, FLAT(:), FGLOB, AHEM, istag=1)
            FGLOB=FGLOB/JWT
         EndIf
         XJL(JM+3,L+KMAX)=AHEM(1)   ! SOUTHERN HEM
         XJL(JM+2,L+KMAX)=AHEM(2)   ! NORTHERN HEM
         XJL(JM+1,L+KMAX)=FGLOB     ! GLOBAL
         WRITE (6,902) PM(L+LM),FGLOB,AHEM(2),AHEM(1),
     *                (MLAT(J),J=JM,J1,-INC)
      END DO

      GO TO 100
  901 FORMAT ('0',30X,A64,'  CP'/1X,32('-'),24A4)
  902 FORMAT (1X,F7.3,3F8.1,1X,24I4)
  903 FORMAT (A6,2X,3F8.1,1X,24I4)
  904 FORMAT (' P(MB)   ',A4,' G      NH      SH  ',24I4)
  905 FORMAT (1X,32('-'),24A4)
  907 FORMAT ('1',A,I3,1X,A3,I5,' - ',I3,1X,A3,I5)
      END SUBROUTINE JKMAP


      SUBROUTINE JLMAP(LNAME,SNAME,UNITS,POW10P,
     &     PL,AX,SCALET,SCALEJ,SCALEL,LMAX,JWT,J1,
     *  ARQX,SCALER,SCALJR,SCALLR)
C****
C**** THIS ROUTINE PRODUCES LAYER BY LATITUDE TABLES ON THE LINE
C**** PRINTER.  THE INTERIOR NUMBERS OF THE TABLE ARE CALCULATED AS
C****               AX * SCALET * SCALEJ * SCALEL.
C**** WHEN JWT=1, THE INSIDE NUMBERS ARE NOT AREA WEIGHTED AND THE
C****    HEMISPHERIC AND GLOBAL NUMBERS ARE SUMMATIONS.
C**** WHEN JWT=2, ALL NUMBERS ARE PER UNIT AREA.
C**** J1 INDICATES PRIMARY OR SECONDARY GRID.
C**** THE BOTTOM LINE IS CALCULATED AS THE SUMMATION OF DSIG TIMES THE
C**** NUMBERS ABOVE (POSSIBLY MULTIPLIED BY A FACTOR OF 10)
C****
      USE DOMAIN_DECOMP_ATM, only : GRID
      use model_com, only: modelEclock
      USE MODEL_COM, only :
     &     JDATE0,AMON,AMON0,JYEAR0,XLABEL
      USE DYNAMICS, only : DSIG,SIGE
      USE GEOM, only :
     &     LAT_DG,WTJ
      USE MDIAG_COM, only : acc_period
     &     ,sname_strlen,units_strlen,lname_strlen
      Use DIAG_COM,  Only: JM,JMby2,LM,LM_REQ,KGZ,QDIAG,INC=>INCJ,LINECT
      IMPLICIT NONE

!@var units string containing output field units
      CHARACTER(LEN=units_strlen) :: UNITS,UNITS_WITH_SCALE
!@var lname string describing output field
      CHARACTER(LEN=lname_strlen) :: LNAME
!@var sname string referencing output field
      CHARACTER(LEN=sname_strlen) :: SNAME
!@var title string, formed as concatentation of lname//units
      CHARACTER(LEN=64) :: TITLE

      REAL*8, DIMENSION(JM) :: FLAT,ASUM
      REAL*8, DIMENSION(2) :: FHEM,HSUM


      INTEGER :: J1,JWT,LMAX
      REAL*8 :: SCALET,SCALER,PRTFAC
      INTEGER :: POW10P
      REAL*8, DIMENSION(JM,LMAX) :: AX
      REAL*8, DIMENSION(JM,LM_REQ) :: ARQX
      REAL*8, DIMENSION(JM) :: SCALEJ,SCALJR
      REAL*8, DIMENSION(:) :: SCALEL
      REAL*8, DIMENSION(:) :: SCALLR
      REAL*8, DIMENSION(:) :: PL

      CHARACTER*4 DASH,WORD(4)
      DATA DASH/'----'/,WORD/'SUM','MEAN',' ','.1*'/

      Logical :: QFINAL  !  produce final line in table
      INTEGER :: IWORD,J,JH,K,L  ,ksx,klmax
      REAL*8 :: FGLOB,GSUM,SDSIG,SUMFAC

      REAL*8, DIMENSION(JM+3,LM+LM_REQ+1+KGZ) :: XJL ! for binary output
      CHARACTER XLB*16,CLAT*16,CPRES*16,CBLANK*16,TITLEO*80,TPOW*8
      DATA CLAT/'LATITUDE'/,CPRES/'PRESSURE (MB)'/,CBLANK/' '/
      optional :: ARQX,SCALER,SCALJR,SCALLR

      integer :: year, date

      call modelEclock%get(year=year, date=date)

      if ( present(ARQX) ) goto 777

      if(sname.eq.'skip') return

      QFINAL = .not. (sname(1:7).eq.'phi_amp' .or.     
     &                sname(1:7).eq.'phi_pha' .or.
     &                sname.eq.'wcod' .or. sname.eq.'icod')

C form title string
      units_with_scale = units
      PRTFAC = 10.**(-pow10p)
      title = trim(lname)//' ('//trim(units)//')'
      if(pow10p.ne.0) then
         write(tpow,'(i3)') pow10p
         tpow='10**'//trim(adjustl(tpow))
         units_with_scale=trim(tpow)//' '//trim(units_with_scale)
         title = trim(lname)//' ('//trim(units_with_scale)//')'
      endif
C****
C**** PRODUCE A LATITUDE BY LAYER TABLE OF THE ARRAY A
C****
   10 LINECT=LINECT+LMAX+7
      IF (LINECT.LE.60) GO TO 20
      WRITE (6,907) XLABEL(1:105),JDATE0,AMON0,JYEAR0,DATE,AMON,YEAR
      LINECT=LMAX+8
   20 WRITE (6,901) TITLE,(DASH,J=J1,JM,INC)
      WRITE (6,904) WORD(JWT),(NINT(LAT_DG(J,J1)),J=JM,J1,-INC)
      WRITE (6,905) (DASH,J=J1,JM,INC)
      XJL(:,:) = -1d30
         KSX = 0            ! KSX = LAYERS GENERATED AT ENTRY
  100 If (QFINAL)  SDSIG=1.-SIGE(LMAX+1)
         KLMAX = LMAX+KSX
      DO 110 J=1,JM
  110 ASUM(J)=0.
      HSUM(1)=0.
      HSUM(2)=0.
      GSUM=0.
      SUMFAC=1.
      IWORD=3
      if(sname.eq.'dudt_mtndrg') then ! make sumfac an argument ... ???
         SUMFAC=10.                   ! ... to avoid this if-block  ???
         IWORD=4
      endif

      FLAT = 0.0
      DO L=LMAX,1,-1
         If (J1==1) then      ! Standard Grid
            DO J=1,JM
               FLAT(J)=AX(J,L)*SCALET*SCALEJ(J)*SCALEL(L)
               XJL(J,L) = FLAT(J)   *PRTFAC
               FLAT(J)=FLAT(J)*PRTFAC
               If (QFINAL)  ASUM(J)=ASUM(J)+FLAT(J)*DSIG(L)/SDSIG
            END DO
            CALL GLOBALSUM(GRID, FLAT(:)*WTJ(:,JWT,J1),
     *                           FGLOB, FHEM)
            FGLOB=FGLOB/JWT
         Else                 ! Staggered Grid
            DO J=2,JM
               FLAT(J)=AX(J,L)*SCALET*SCALEJ(J)*SCALEL(L)
               XJL(J,L) = FLAT(J)   *PRTFAC
               FLAT(J)=FLAT(J)*PRTFAC
               If (QFINAL)  ASUM(J)=ASUM(J)+FLAT(J)*DSIG(L)/SDSIG
            END DO
            CALL GLOBALSUM(GRID, FLAT(:)*WTJ(:,JWT,J1),
     *                           FGLOB, FHEM, istag=1)
            FGLOB=FGLOB/JWT
         EndIf
         XJL(JM+3,L)=FHEM(1)   ! SOUTHERN HEM
         XJL(JM+2,L)=FHEM(2)   ! NORTHERN HEM
         XJL(JM+1,L)=FGLOB     ! GLOBAL
      WRITE (6,902) PL(L),FGLOB,FHEM(2),FHEM(1),
     &        (NINT(MIN(1d5,MAX(-1d5,FLAT(J)))),J=JM,J1,-INC)
         CALL KEYNRL (SNAME,L,FLAT)
         If (QFINAL)  HSUM(1)=HSUM(1)+FHEM(1)*SUMFAC*DSIG(L)/SDSIG
         If (QFINAL)  HSUM(2)=HSUM(2)+FHEM(2)*SUMFAC*DSIG(L)/SDSIG
         If (QFINAL)  GSUM=GSUM+FGLOB*SUMFAC*DSIG(L)/SDSIG
      END DO

      WRITE (6,905) (DASH,J=J1,JM,INC)
      If (.not. QFINAL)  Return
cBMP      ASUM(jmby2+1)=ASUM(jmby2+1)/J1
         DO 180 J=J1,JM
  180    XJL(J   ,LM+LM_REQ+1)=ASUM(J)
         XJL(JM+3,LM+LM_REQ+1)=HSUM(1)/SUMFAC   ! SOUTHERN HEM
         XJL(JM+2,LM+LM_REQ+1)=HSUM(2)/SUMFAC   ! NORTHERN HEM
         XJL(JM+1,LM+LM_REQ+1)=GSUM/SUMFAC      ! GLOBAL
         XLB=' '//acc_period(1:3)//' '//acc_period(4:12)//'  '
         TITLEO=TITLE//XLB
         IF(QDIAG) CALL POUT_JL(TITLEO,LNAME,SNAME,UNITS_WITH_SCALE,
     *        J1,KLMAX,XJL,PL,CLAT,CPRES)
      WRITE (6,903) WORD(IWORD),GSUM,HSUM(2),HSUM(1),
     *  (NINT(MIN(1d5,MAX(-1d5,ASUM(J)*SUMFAC))),J=JM,J1,-INC)
      RETURN
C****
!      ENTRY JLMAPS(LNAME,SNAME,UNITS,POW10P,
!     &     PL,AX,SCALET,SCALEJ,SCALEL,LMAX,JWT,J1,
!     *  ARQX,SCALER,SCALJR,SCALLR)
 777  continue

      if(sname.eq.'skip') return
C form title string
      units_with_scale = units
      title = trim(lname)//' ('//trim(units)//')'
      PRTFAC = 10.**(-pow10p)
      if(pow10p.ne.0) then
         write(tpow,'(i3)') pow10p
         tpow='10**'//trim(adjustl(tpow))
         units_with_scale=trim(tpow)//' '//trim(units_with_scale)
         title = trim(lname)//' ('//trim(units_with_scale)//')'
      endif
         KSX = 3
         DO 205 L=1,LM+LM_REQ
         DO 205 J=1,JM
  205    XJL(J,L) = -1.D30
      LINECT=LINECT+LMAX+10
      IF (LINECT.LE.60) GO TO 200
      WRITE (6,907) XLABEL(1:105),JDATE0,AMON0,JYEAR0,DATE,AMON,YEAR
      LINECT=LMAX+11
  200 CONTINUE
C**** PRODUCE UPPER STRATOSPHERE NUMBERS FIRST
      WRITE (6,901) TITLE,(DASH,J=J1,JM,INC)
      WRITE (6,904) WORD(JWT),(NINT(LAT_DG(J,J1)),J=JM,J1,-INC)
      WRITE (6,905) (DASH,J=J1,JM,INC)
      DO L=LM_REQ,1,-1
         If (J1==1) then     ! Standard Grid
            DO J=1,JM
               FLAT(J)=ARQX(J,L)*SCALER*SCALJR(J)*SCALLR(L)
               XJL(J,L+LMAX) = FLAT(J)
c              FLAT(J)=FLAT(J)*PRTFAC
            END DO
            CALL GLOBALSUM(GRID, FLAT(:)*WTJ(:,JWT,J1),
     *                           FGLOB, FHEM)
            FGLOB=FGLOB/JWT
         Else                ! Staggered Grid
            DO J=2,JM
               FLAT(J)=ARQX(J,L)*SCALER*SCALJR(J)*SCALLR(L)
               XJL(J,L+LMAX) = FLAT(J)
c              FLAT(J)=FLAT(J)*PRTFAC
            END DO
            CALL GLOBALSUM(GRID, FLAT(:)*WTJ(:,JWT,J1),
     *                           FGLOB, FHEM, istag=1)
            FGLOB=FGLOB/JWT
         EndIf
         XJL(JM+3,L+LMAX)=FHEM(1)   ! SOUTHERN HEM
         XJL(JM+2,L+LMAX)=FHEM(2)   ! NORTHERN HEM
         XJL(JM+1,L+LMAX)=FGLOB     ! GLOBAL
  230 WRITE (6,902) PL(L+LM),FGLOB,FHEM(2),FHEM(1),
     *  (NINT(MIN(1d5,MAX(-1d5,FLAT(J)))),J=JM,J1,-INC)
      END DO
      GO TO 100
  901 FORMAT ('0',30X,A64/2X,32('-'),24A4)
  902 FORMAT (1X,F8.3,3F8.1,1X,24I4)
  903 FORMAT (1X,A6,2X,3F8.1,1X,24I4)
  904 FORMAT ('  P(MB)   ',A4,' G      NH      SH  ',24I4)
  905 FORMAT (2X,32('-'),24A4)
  907 FORMAT ('1',A,I3,1X,A3,I5,' - ',I3,1X,A3,I5)
      END SUBROUTINE JLMAP

      subroutine print_generic_jl(jl_index,sjl_index,vsumfac)
      use constant, only : teeny
      use atm_com, only : lm_req
      use model_com, only: modelEclock
      use model_com, only : idacc,
     &     jdate0,amon0,amon,jyear0,xlabel
      use dynamics, only : dsig
      use geom, only : lat_dg,dxyp
      use mdiag_com, only : acc_period
     &     ,sname_strlen,units_strlen,lname_strlen
      use diag_com, only : jm,lm,qdiag,inc=>incj,linect,
     &     ajl,denom_jl,ia_jl,scale_jl,sname_jl,lname_jl,units_jl,
     &     pow_jl,lgrid_jl,ctr_ml,edg_ml,ctr_cp,edg_cp,plm,ple,
     &     asjl,scale_sjl,ia_sjl
      implicit none
      integer, intent(in) :: jl_index
      integer, intent(in), optional :: sjl_index
      real*8, intent(in), optional :: vsumfac
      real*8, parameter :: missing=-1d30
      real*8 :: scalet,prtfac
      real*8, dimension(lm) :: wtl
      real*8, dimension(lm+lm_req) :: pm
c      real*8, dimension(:,:), allocatable :: anum,aden,xjl
      real*8, dimension(jm+3,lm+lm_req+1) :: anum,aden,xjl
      real*8, dimension(jm+3) :: flat
      real*8 :: fsh,fnh,fglob,vsfac
      integer :: j,l,lmax,jl_denom,lgrid

      character(len=4) :: dash='----'
      character(len=16) :: clat='LATITUDE'
      character(len=16) :: cpres='PRESSURE (MB)'
      character xlb*16,titleo*80,tpow*8
      character(len=units_strlen) :: units,units_with_scale
      character(len=lname_strlen) :: lname
      character(len=sname_strlen) :: sname
      character(len=64) :: title
      character(len=4) :: vsword

      integer :: year, date

      call modelEclock%get(year=year, date=date)

      if(present(vsumfac)) then
        vsfac = vsumfac
        if(vsfac.eq.10. ) vsword=' .1*'
        if(vsfac.eq.100.) vsword='.01*'
      else
        vsfac = 1d0
        vsword = '    '
      endif
c
c check whether radiation equilibrium layers are to be output
c
      if(present(sjl_index)) then
        lmax = lm+lm_req
      else
        lmax = lm
      endif

c      allocate(anum(jm+3,lmax+1),aden(jm+3,lmax+1),xjl(jm+3,lmax+1))

c
c create weighting for vertical averages if none has been declared.
c dsig is not the correct weight for qtys at layer edges,
c but we simply follow tradition here.
c
      lgrid = lgrid_jl(jl_index)
      jl_denom = denom_jl(jl_index)
      if(jl_denom.eq.0) then
        wtl(:) = dsig(:)
      else
        wtl(:) = 1.
      endif

c
c compose the numerator
c
      anum = 0d0
      prtfac = 10.**(-pow_jl(jl_index))
      scalet = prtfac*scale_jl(jl_index)/idacc(ia_jl(jl_index))
      do l=1,lm
      do j=1,jm
        anum(j,l) = ajl(j,l,jl_index)*scalet*dxyp(j)*wtl(l)
      enddo
      enddo
      if(present(sjl_index)) then
        scalet = scale_sjl(sjl_index)/idacc(ia_sjl(sjl_index))
        do l=1,lm_req
        do j=1,jm
          anum(j,lm+l) = asjl(j,l,sjl_index)*scalet*dxyp(j)
        enddo
        enddo
      endif
      anum(1:jm,lm+lm_req+1) = sum(anum(1:jm,1:lm),dim=2)
      anum(jm+3,:) = sum(anum(1:jm/2,:),dim=1)
      anum(jm+2,:) = sum(anum(jm/2+1:jm,:),dim=1)
      anum(jm+1,:) = anum(jm+2,:)+anum(jm+3,:)

c
c compose the denominator
c
      aden = 0d0
      do l=1,lm
        aden(1:jm,l) = dxyp(1:jm)*wtl(l)
      enddo
      do l=lm+1,lmax
        aden(1:jm,l) = dxyp(1:jm)
      enddo
      if(jl_denom.gt.0) then
        scalet = 1d0/(idacc(ia_jl(jl_denom))+teeny)
        do l=1,lm
        do j=1,jm
          aden(j,l) = aden(j,l)*ajl(j,l,jl_denom)*scalet
        enddo
        enddo
      endif
      aden(1:jm,lm+lm_req+1) = sum(aden(1:jm,1:lm),dim=2)
      aden(jm+3,:) = sum(aden(1:jm/2,:),dim=1)
      aden(jm+2,:) = sum(aden(jm/2+1:jm,:),dim=1)
      aden(jm+1,:) = aden(jm+2,:)+aden(jm+3,:)

c
c compose the ratio
c
      where(aden.ne.0.)
        xjl = anum/aden
      elsewhere
        xjl = missing
      end where

c
c copy metadata into local strings
c
      sname = sname_jl(jl_index)
      lname = lname_jl(jl_index)
      units = units_jl(jl_index)

c
c choose a vertical coordiate scale
c
      if(lgrid.eq.ctr_ml .or. lgrid.eq.ctr_cp) then
        pm(1:lmax) = plm(1:lmax)
      else
        pm(1:lm) = ple(1:lm)
      endif

c
c form title string
c
      units_with_scale = units
      title = trim(lname)//' ('//trim(units)//')'
      if(pow_jl(jl_index).ne.0) then
         write (tpow, '(i3)') pow_jl(jl_index)
         tpow='10**'//trim(adjustl(tpow))
         units_with_scale=trim(tpow)//' '//trim(units_with_scale)
         title = trim(lname)//' ('//trim(units_with_scale)//')'
      endif
      if(lgrid.eq.ctr_cp .or. lgrid.eq.edg_cp) title(61:64)='(CP)'

      linect = linect + lmax + 7
      if(linect.gt.60) then
        WRITE(6,907) XLABEL(1:105)
     &       ,JDATE0,AMON0,JYEAR0,DATE,AMON,YEAR
        linect = lmax+8
      endif

c
c print table
c
      WRITE (6,901) TITLE,(DASH,J=1,JM,INC)
      WRITE (6,904) 'MEAN',(NINT(LAT_DG(J,1)),J=JM,1,-INC)
      WRITE (6,905) (DASH,J=1,JM,INC)
      do l=lmax,1,-1
        where(xjl(:,l).ne.missing)
          flat(:) = xjl(:,l)
        elsewhere
          flat(:) = 0.
        end where
        fglob=flat(jm+1)
        fnh  =flat(jm+2)
        fsh  =flat(jm+3)
        write(6,902) pm(l),fglob,fnh,fsh,(nint(flat(j)),j=jm,inc,-inc)
      enddo
      WRITE(6,905) (DASH,J=1,JM,INC)
      l = lm+lm_req+1
      where(xjl(:,l).ne.missing)
        flat(:) = xjl(:,l)
      elsewhere
        flat(:) = 0.
      end where
      call keyvsums(sname,flat(jm+1),flat(jm+3:jm+2:-1),flat,1d0)
      flat = flat*vsfac
      fglob=flat(jm+1)
      fnh  =flat(jm+2)
      fsh  =flat(jm+3)
      write(6,903) vsword,fglob,fnh,fsh,(nint(flat(j)),j=jm,inc,-inc)

c
c write to binary file if requested
c
      if(qdiag) then
        XLB=' '//acc_period(1:3)//' '//acc_period(4:12)//'  '
        TITLEO=TITLE//XLB
        CALL POUT_JL(TITLEO,LNAME,SNAME,UNITS_WITH_SCALE,
     *        1,LMAX,XJL,PM,CLAT,CPRES)
      endif

c      deallocate(anum,aden,xjl)

      return
  901 FORMAT ('0',30X,A64/2X,32('-'),24A4)
  902 FORMAT (1X,F8.3,3F8.1,1X,24I4)
  903 FORMAT (1X,A6,2X,3F8.1,1X,24I4)
  904 FORMAT ('  P(MB)   ',A4,' G      NH      SH  ',24I4)
  905 FORMAT (2X,32('-'),24A4)
  907 FORMAT ('1',A,I3,1X,A3,I5,' - ',I3,1X,A3,I5)
      END SUBROUTINE PRINT_GENERIC_JL

      SUBROUTINE DIAGIL
!@sum  DIAGIL prints out longitude/height diagnostics
!@auth Original Development Team
      USE MODEL_COM, only : idacc,xlabel,lrunid,dtsrc
      use dynamics, only : bydsig
      use mdiag_com, only : acc_period
     &     ,sname_strlen,units_strlen,lname_strlen
      USE DIAG_COM, only : im,lm,aijl,lm_req,qdiag
     &     ,ia_src,ia_rad,ia_dga,plm,ple,linect
     &     ,IJL_U,IJL_V,IJK_TX,IJL_W,IJK_RH,IJL_MC
     &     ,ia_ijl,denom_ijl
      USE DIAG_COM_RAD, only : IJL_RC
      USE GC_COM, only : j5s,j5n,j5suv,j5nuv,j50n,j70n
      USE CONSTANT, only : grav,rgas,by3,sha,bygrav,teeny
      USE GEOM, only : dxyp
      IMPLICIT NONE
      CHARACTER(len=sname_strlen) :: sname
      CHARACTER(len=units_strlen) :: unit
      CHARACTER(len=lname_strlen) :: lname
      REAL*8, DIMENSION(LM) :: ONES
      REAL*8, DIMENSION(IM,LM) :: XIL
      INTEGER, PARAMETER :: KAILX=15
      character(len=sname_strlen), dimension(kailx) :: name_il
      character(len=units_strlen), dimension(kailx) :: units_il
      character(len=lname_strlen), dimension(kailx) :: lname_il
      real*8, dimension(kailx) :: scale_il
      integer, dimension(kailx) :: ia_il,j1_il,j2_il,qty_il
      real*8 :: bydj,bydjuv
      integer :: k,kd,j
c
      do k=1,kailx
         write(name_il(k),'(a3,i3.3)') 'AIL',k
         lname_il(k) = 'unused'
         units_il(k) = 'unused'
         scale_il(k) = 1.
         ia_il(k)    = 0
      enddo

C**** some scaling numbers for the equatorial diags.
      bydj   = 1./REAL(j5n-j5s+1,KIND=8)
      bydjuv = 1./REAL(j5nuv-j5suv+1,KIND=8)
C****
      k=0
c
      k = k + 1
      qty_il(k) = IJL_U
      name_il(k) = 'u_equator'
      lname_il(k) = 'ZONAL WIND (U COMPONENT) AROUND +/- 5 DEG'
      units_il(k) = 'm/s'
      scale_il(k) = bydjuv
      ia_il(k)    = ia_dga
      j1_il(k) = j5suv
      j2_il(k) = j5nuv
      k = k + 1
      qty_il(k) = IJL_V
      name_il(k) = 'v_equator'
      lname_il(k) = 'MERIDIONAL WIND (V COMPONENT) AROUND +/- 5 DEG'
      units_il(k) = 'm/s'
      scale_il(k) = bydjuv
      ia_il(k)    = ia_dga
      j1_il(k) = j5suv
      j2_il(k) = j5nuv
      k = k + 1
      qty_il(k) = IJL_W
      name_il(k) = 'vvel_equator'
      lname_il(k) = 'VERTICAL VELOCITY AROUND +/- 5 DEG'
      units_il(k) = '10**-4 m/s'
      scale_il(k) = -1d4*RGAS*BYGRAV*bydj
      ia_il(k)    = ia_dga
      j1_il(k) = j5s
      j2_il(k) = j5n
      k = k + 1
      qty_il(k) = IJK_TX
      name_il(k) = 'temp_equator'
      lname_il(k) = 'TEMPERATURE AROUND +/- 5 DEG  (Const. Pres.)'
      units_il(k) = 'C'
      scale_il(k) = 1. !bydj
      ia_il(k)    = ia_dga
      j1_il(k) = j5s
      j2_il(k) = j5n
      k = k + 1
      qty_il(k) = IJK_RH
      name_il(k) = 'rh_equator'
      lname_il(k) = 'RELATIVE HUMIDITY AROUND +/- 5 DEG (Const. Pres.)'
      units_il(k) = '%'
      scale_il(k) = 1d2 !*bydj
      ia_il(k)    = ia_dga
      j1_il(k) = j5s
      j2_il(k) = j5n
      k = k + 1
      qty_il(k) = IJL_MC
      name_il(k) = 'mcheat_equator'
      lname_il(k) = 'MOIST CONVECTIVE HEATING AROUND +/- 5 DEG'
      units_il(k) =  '.01 W/(m^2*mb)'           !'10**13 WATTS/DSIG'
      scale_il(k) =  100.*100.*BYGRAV*SHA/DTsrc !100d-13*SHA/(GRAV*DTsrc)
      ia_il(k)    = ia_src
      j1_il(k) = j5s
      j2_il(k) = j5n
      k = k + 1
      qty_il(k) = IJL_RC
      name_il(k) = 'rad_cool_equator'
      lname_il(k) = 'TOTAL RADIATIVE COOLING AROUND +/- 5 DEG'
      units_il(k) =  '.01 W/(m^2*mb)' !'10**13 WATTS/DSIG'
      scale_il(k) = -100.             !-1d-13
      ia_il(k)    = ia_rad
      j1_il(k) = j5s
      j2_il(k) = j5n
      k = k + 1
      qty_il(k) = IJL_W
      name_il(k) = 'vvel_50N'
      lname_il(k) = 'VERTICAL VELOCITY AT 50 N'
      units_il(k) = '10**-4 m/s'
      scale_il(k) = -1d4*RGAS*BYGRAV
      ia_il(k)    = ia_dga
      j1_il(k) = j50n
      j2_il(k) = j50n
      k = k + 1
      qty_il(k) = IJK_TX
      name_il(k) = 'temp_50N'
      lname_il(k) = 'TEMPERATURE AT 50 N (Const. Pres.)'
      units_il(k) = 'C'
      scale_il(k) = 1.
      ia_il(k)    = ia_dga
      j1_il(k) = j50n
      j2_il(k) = j50n
      k = k + 1
      qty_il(k) = IJL_RC
      name_il(k) = 'rad_cool_50N'
      lname_il(k) = 'TOTAL RADIATIVE COOLING AT 50 N'
      units_il(k) = '.01 W/(m^2*mb)'  !'10**13 WATTS/UNIT SIGMA'
      scale_il(k) = -100.             ! 1d-13
      ia_il(k)    = ia_rad
      j1_il(k) = j50n
      j2_il(k) = j50n
      k = k + 1
      qty_il(k) = IJL_U
      name_il(k) = 'u_50N'
      lname_il(k) = 'ZONAL WIND AT 50 N'
      units_il(k) = 'm/s'
      scale_il(k) = 0.5
      ia_il(k)    = ia_dga
      j1_il(k) = j50n
      j2_il(k) = j50n+1
      k = k + 1
      qty_il(k) = IJL_W
      name_il(k) = 'vvel_70N'
      lname_il(k) = 'VERTICAL VELOCITY AT 70 N'
      units_il(k) = '10**-4 m/s'
      scale_il(k) = -1d4*RGAS*BYGRAV
      ia_il(k)    = ia_dga
      j1_il(k) = j70n
      j2_il(k) = j70n
      k = k + 1
      qty_il(k) = IJK_TX
      name_il(k) = 'temp_70N'
      lname_il(k) = 'TEMPERATURE AT 70 N (Const. Pres.)'
      units_il(k) = 'C'
      scale_il(k) = 1.
      ia_il(k)    = ia_dga
      j1_il(k) = j70n
      j2_il(k) = j70n
      k = k + 1
      qty_il(k) = IJL_RC
      name_il(k) = 'rad_cool_70N'
      lname_il(k) = 'TOTAL RADIATIVE COOLING AT 70 N'
      units_il(k) = '.01 W/(m^2*mb)' !'10**13 WATTS/UNIT SIGMA'
      scale_il(k) = -100.            !-1d-13
      ia_il(k)    = ia_rad
      j1_il(k) = j70n
      j2_il(k) = j70n
      k = k + 1
      qty_il(k) = IJL_U
      name_il(k) = 'u_70N'
      lname_il(k) = 'ZONAL WIND AT 70 N'
      units_il(k) = 'm/s'
      scale_il(k) = 0.5
      ia_il(k)    = ia_dga
      j1_il(k) = j70n
      j2_il(k) = j70n+1
c

C**** OPEN PLOTTABLE OUTPUT FILE IF DESIRED
      IF(QDIAG) call open_il(trim(acc_period)//'.il'//XLABEL(1:LRUNID)
     *     ,im,lm,lm_req)

C**** INITIALIZE CERTAIN QUANTITIES
      ONES(1:LM)=1.

      linect = 65

      DO K=1,KAILX
        sname=name_il(k)
        lname=lname_il(k)
        unit=units_il(k)
        if (lname.ne.'unused') then
        XIL=SUM(AIJL(:,j1_il(k):j2_il(k),:,qty_il(k)),dim=2)
     &         *SCALE_IL(K)/IDACC(IA_IL(K))
        kd = denom_ijl(qty_il(k))
        if(kd.gt.0) then ! weighted by air mass
          XIL = XIL*IDACC(IA_IJL(KD))/
     &         (SUM(AIJL(:,j1_il(k):j2_il(k),:,kd),dim=2)+teeny)
        endif
        SELECT CASE (sname)
! Centered in L; primary grid; hor. mean; vert. sum
        CASE ('u_equator','v_equator','u_70N','u_50N')
          CALL ILMAP(sname,lname,unit,PLM,XIL,ONES,LM,2,1)
! Vertical edges; primary grid; hor. mean; vert. sum
        CASE ('vvel_equator','vvel_50N','vvel_70N')
          CALL ILMAP(sname,lname,unit,PLE,XIL,ONES,LM-1,2,1)
! Centered in L; primary grid; hor. mean; vert. sum
        CASE ('temp_equator','rh_equator','temp_50N','temp_70N')
          CALL ILMAP(sname,lname,unit,PLM,XIL,ONES,LM,2,1)
! Centered in L; primary grid; hor. mean; vert. sum
        CASE ('mcheat_equator')
          CALL ILMAP(sname,lname,unit,PLM,XIL,ONES,LM,2,1)
! Centered in L; primary grid; hor. mean; vert. mean
        CASE ('rad_cool_equator') ! also 'rad_cool_50N','rad_cool_70N'
          CALL ILMAP(sname,lname,unit,PLM,XIL,ONES,LM,2,1)
        END SELECT
        end if
      END DO
      if(qdiag) call close_il
      RETURN
      END SUBROUTINE DIAGIL


      SUBROUTINE ILMAP (sname,lname,unit,PL,AX,SCALEL,LMAX,JWT
     *     ,ISHIFT)
      USE CONSTANT, only : twopi
      use model_com, only: modelEclock
      USE MODEL_COM, only : jdate0,amon,amon0
     *     ,jyear0,xlabel
      use dynamics, only : dsig,sige
      USE GEOM, only : dlon,lon_dg
      use mdiag_com, only : acc_period
     &     ,sname_strlen,units_strlen,lname_strlen
      USE DIAG_COM, only : im,jm,lm,qdiag,inc=>inci,linect
      IMPLICIT NONE
      CHARACTER XLB*80,CWORD*8
      character(len=sname_strlen), intent(in) :: sname
      character(len=units_strlen), intent(in) :: unit
      character(len=lname_strlen), intent(in) :: lname
      CHARACTER*64 :: TITLE
      CHARACTER*4, PARAMETER :: DASH='----'
      CHARACTER*4, DIMENSION(2), PARAMETER :: WORD=(/'SUM ','MEAN'/)
      CHARACTER*16, PARAMETER :: CBLANK=' ', CLAT='LONGITUDE',
     *     CPRES='PRESSURE (MB)'
      REAL*8, DIMENSION(LM), INTENT(IN) :: PL,SCALEL
      REAL*8, DIMENSION(IM,LM), INTENT(IN) :: AX
      REAL*8 :: XIL(IM,LM),ZONAL(LM) ! used for post-proc
      REAL*8 :: FGLOB,FLON,GSUM,SDSIG
      REAL*8, DIMENSION(IM) :: ASUM
      INTEGER, DIMENSION(IM) :: MLON
      INTEGER, INTENT(IN) :: JWT,ISHIFT
      INTEGER :: I,L,LMAX
      integer :: year, date

      call modelEclock%get(year=year, date=date)

C****
C**** PRODUCE A LONGITUDE BY LAYER TABLE OF THE ARRAY A
C****
!@var ISHIFT: When=2, print longitude indices off center (U-grid)
      LINECT=LINECT+LMAX+7
      IF (LINECT.GT.60) THEN
        WRITE (6,907) XLABEL(1:105),JDATE0,AMON0,JYEAR0,DATE,AMON,YEAR
        LINECT=LMAX+8
      END IF
      SDSIG=1.-SIGE(LMAX+1)

      TITLE=trim(lname)//" ("//trim(unit)//")"
      WRITE (6,901) TITLE,(DASH,I=1,IM,INC)
      IF (ISHIFT.EQ.1) WRITE (6,904) WORD(JWT),(I,I=1,IM,INC)
      IF (ISHIFT.EQ.2) WRITE (6,906) WORD(JWT),(I,I=1,IM,INC)
      WRITE (6,905) (DASH,I=1,IM,INC)

      ASUM(:)=0. ; GSUM=0.
      DO L=LMAX,1,-1
        FGLOB=0.
        DO I=1,IM
          FLON=AX(I,L)*SCALEL(L)
          XIL(I,L)=FLON
          MLON(I)=NINT(MIN(1d5,MAX(-1d5,FLON)))
          ASUM(I)=ASUM(I)+FLON*DSIG(L)/SDSIG
          FGLOB=FGLOB+FLON
        END DO
        FGLOB=FGLOB/IM
        IF (JWT.EQ.1) FGLOB=FGLOB*TWOPI/DLON
        ZONAL(L)=FGLOB
        WRITE (6,902) PL(L),FGLOB,(MLON(I),I=1,IM,INC)
        GSUM=GSUM+FGLOB*DSIG(L)/SDSIG
      END DO
      MLON(:)=NINT(ASUM(:))

      WRITE (6,905) (DASH,I=1,IM,INC)
      WRITE (6,903) GSUM,(MLON(I),I=1,IM,INC)
C**** Output for post-processing
      CWORD=WORD(JWT)           ! pads out to 8 characters
      XLB=TITLE
      XLB(65:80)=' '//acc_period(1:3)//' '//acc_period(4:12)//'  '
      IF(QDIAG) CALL POUT_IL(XLB,sname,lname,unit,1,ISHIFT,LMAX,XIL
     *     ,PL,CLAT,CPRES,ASUM,GSUM,ZONAL)
      RETURN
C****
  901 FORMAT ('0',30X,A64,2('-'),/1X,16('-'),36A3)
  902 FORMAT (F8.3,F8.1,1X,36I3)
  903 FORMAT (F16.1,1X,36I3)
  904 FORMAT (' P(MB)',6X,A4,1X,36I3)  ! U-grid (i.e., centers)
  905 FORMAT (1X,16('-'),36A3)
  906 FORMAT (' P(MB)',6X,A4,36I3)     ! V-grid (i.e., edges)
  907 FORMAT ('1',A,I3,1X,A3,I5,' - ',I3,1X,A3,I5)
C****
      END SUBROUTINE ILMAP

      SUBROUTINE DIAG7P
C****
C**** THIS ENTRY PRINTS THE TABLES
C****
      use model_com, only: modelEclock
      USE MODEL_COM, only :
     &     IDACC,JDATE0,AMON,AMON0,JYEAR0,XLABEL,lrunid
      USE DIAG_COM, only : qdiag,im,ia_12hr,ia_inst
      USE gc_COM, only :
     &     nwav_dag,wave,Max12HR_sequ,Min12HR_sequ
      USE MDIAG_COM, only : acc_period,
     &     sname_strlen,units_strlen,lname_strlen
      IMPLICIT NONE

      INTEGER, DIMENSION(44) :: IPOWER
      REAL*8, DIMENSION(120) :: POWER
      REAL*8, DIMENSION(43) :: XPOWER
      REAL*8, DIMENSION(13) :: FPE
!     Arrays for pdE
      REAL*8, DIMENSION(43+1,NWAV_DAG+1) :: FPOWER
      REAL*8, DIMENSION(41,2) :: period_e
      REAL*8, DIMENSION(nwav_dag) :: xnwav
!@var COMP_WAVE complex form of WAVE. Correct arg. to subr. MEM
      COMPLEX*16, DIMENSION(Max12HR_sequ) :: COMP_WAVE
      CHARACTER XLB*14,CLAT*16,CPRES*16,CBLANK*16,TITLEO*80
      DATA CLAT/'PERIOD EASTWARD'/,CPRES/'N'/,CBLANK/' '/

      INTEGER, PARAMETER :: MMAX=12,NUAMAX=120,NUBMAX=15

      COMMON/D7COM/LNAME,SNAME,UNITS
      CHARACTER TITLE(12)*66
      CHARACTER(len=lname_strlen) :: LNAME(12)
      CHARACTER(len=sname_strlen) :: SNAME(12)
      CHARACTER(len=units_strlen) ::  UNITS(12)

      REAL*8, DIMENSION(12) :: SCALET
      DATA SCALET/1.,1., .1,1., .1,1., 4*1.D-3,1.D-4,1.D-5/

      REAL*8 :: BYIA12,PNU,POWX,VAR

      INTEGER ::
     &     IDACC9,K,KPAGE,KQ,KTABLE,
     &     M,MMAXP1,N,NMAX,NS,NUA,NX
      INTEGER :: ic

      integer :: year, date

      call modelEclock%get(year=year, date=date)

      NMAX=NWAV_DAG
      IDACC9=IDACC(ia_12hr)
      IF (IDACC9.LE.MMAX) RETURN
C**** PATCH NEEDED IF SEVERAL RESTART FILES WERE ACCUMULATED
      IF (IDACC(ia_inst).LE.1) GO TO 320
      IDACC9=Min12HR_sequ           ! in case a February was included
      BYIA12=1./IDACC(ia_inst)
      WAVE(:,:,:,:)=WAVE(:,:,:,:)*BYIA12
  320 CONTINUE
      IF (IDACC9.GT.Max12HR_sequ) IDACC9=Max12HR_sequ

C**** OPEN PLOTTABLE OUTPUT FILE IF DESIRED
C**** NOTE: there are 41 periods.  Fpower is padded for POUT
C**** Save inverse period (day) so coordinate is monotonic
      IF(QDIAG) then
        do k=1,41
          period_e(41+1-k,1) = (k-25)  !/60.
          period_e(41+1-k,2) = (k-17)  !/60.
        end do
          call open_wp(trim(acc_period)//'.wp'//XLABEL(1:LRUNID)
     *     ,41,NMAX,0,period_e)
      do n=1,nmax
        xnwav(n) = n
      end do
      XLB=' '//acc_period(1:3)//' '//acc_period(4:12)
      fpower = -1.E30
      end if
C****
C**** OUTPUT WAVE POWER AT THE EQUATOR
C****
      MMAXP1=MMAX+1
      DO 400 KPAGE=1,2
      WRITE (6,907) XLABEL(1:105),JDATE0,AMON0,JYEAR0,DATE,AMON,YEAR
      DO 390 KTABLE=1,3
      KQ=3*(KPAGE-1)+KTABLE
      TITLE(KQ)=TRIM(LNAME(KQ))//" ("//TRIM(UNITS(KQ))//") "
      WRITE (6,901) TITLE(KQ)
      DO 380 NX=1,NMAX
      N=NMAX+1-NX
      do ic=1,Max12HR_sequ
        comp_wave(ic)=cmplx( WAVE(1,ic,N,KQ) , WAVE(2,ic,N,KQ) )
      end do
      CALL MEM (COMP_WAVE,IDACC9,MMAX,NUAMAX,NUBMAX,POWER,FPE,
     *  VAR,PNU)
      POWX=.5*POWER(1)
      DO 330 NUA=2,27
  330 POWX=POWX+POWER(NUA)
      XPOWER(1)=SCALET(KQ)*POWX/26.5
      POWX=0.
      DO 340 NUA=28,34
  340 POWX=POWX+POWER(NUA)
      XPOWER(2)=SCALET(KQ)*POWX/7.
      XPOWER(3)=SCALET(KQ)*(POWER(35)+POWER(36)+POWER(37)+POWER(38))/4.
      XPOWER(4)=SCALET(KQ)*(POWER(39)+POWER(40))/2.
      DO 350 NUA=41,76
  350 XPOWER(NUA-36)=SCALET(KQ)*POWER(NUA)
      POWX=.5*POWER(1)
      DO 360 NUA=77,120
  360 POWX=POWX+POWER(NUA)
      XPOWER(41)=SCALET(KQ)*POWX/44.5
      XPOWER(42)=10.*SCALET(KQ)*VAR
      XPOWER(43)=1000.*SCALET(KQ)*(VAR-PNU)
      DO 370 NS=1,43
      IPOWER(NS)=XPOWER(NS)+.5
  370 CONTINUE
        DO NS=1,41
          FPOWER(41+1-NS,N)=XPOWER(NS)
        END DO
          FPOWER(42,N)=XPOWER(42)
          FPOWER(43,N)=XPOWER(43)
  380 WRITE (6,902) N,(IPOWER(NS),NS=1,43)
      WRITE (6,903) (FPE(M),M=1,MMAXP1)
      TITLEO=TITLE(KQ)//'*60day'//XLB
      IF(QDIAG) CALL POUT_WP(TITLEO,LNAME(KQ),SNAME(KQ),UNITS(KQ),
     *     1,NMAX,FPOWER,xnwav,CLAT,CPRES)
  390 CONTINUE
  400 CONTINUE
C****
C**** OUTPUT WAVE POWER AT 50 DEG NORTH
C****
      DO 500 KPAGE=3,4
      WRITE (6,907) XLABEL(1:105),JDATE0,AMON0,JYEAR0,DATE,AMON,YEAR
      DO 490 KTABLE=1,3
      KQ=3*(KPAGE-1)+KTABLE
      TITLE(KQ)=TRIM(LNAME(KQ))//" ("//TRIM(UNITS(KQ))//") "
  410 WRITE (6,911) TITLE(KQ)
      DO 480 NX=1,NMAX
      N=NMAX+1-NX
      do ic=1,Max12HR_sequ
        comp_wave(ic)=cmplx( WAVE(1,ic,N,KQ) , WAVE(2,ic,N,KQ) )
      end do
      CALL MEM (COMP_WAVE,IDACC9,MMAX,NUAMAX,NUBMAX,POWER,FPE,
     *  VAR,PNU)
      DO 420 M=1,MMAXP1
  420 FPE(M)=1000.*SCALET(KQ)*FPE(M)
      POWX=.5*POWER(1)
      DO 430 NUA=2,45
  430 POWX=POWX+POWER(NUA)
      XPOWER(1)=SCALET(KQ)*POWX/44.5
      DO 440 NUA=46,81
  440 XPOWER(NUA-44)=SCALET(KQ)*POWER(NUA)
      XPOWER(38)=SCALET(KQ)*(POWER(82)+POWER(83))/2.
      XPOWER(39)=SCALET(KQ)*(POWER(84)+POWER(85)+POWER(86)+POWER(87))/4.
      POWX=0.
      DO 450 NUA=88,94
  450 POWX=POWX+POWER(NUA)
      XPOWER(40)=SCALET(KQ)*POWX/7.
      POWX=.5*POWER(1)
      DO 460 NUA=95,120
  460 POWX=POWX+POWER(NUA)
      XPOWER(41)=SCALET(KQ)*POWX/26.5
      XPOWER(42)=10.*SCALET(KQ)*VAR
      XPOWER(43)=1000.*SCALET(KQ)*(VAR-PNU)
      DO 470 NS=1,43
      IPOWER(NS)=XPOWER(NS)+.5
  470 CONTINUE
        DO NS=1,41
          FPOWER(41+1-NS,N)=XPOWER(NS)
        END DO
          FPOWER(42,N)=XPOWER(42)
          FPOWER(43,N)=XPOWER(43)
  480 WRITE (6,902) N,(IPOWER(NS),NS=1,43)
      WRITE (6,903) (FPE(M),M=1,MMAXP1)
      TITLEO=TITLE(KQ)//'*60day'//XLB
      IF(QDIAG) CALL POUT_WP(TITLEO,LNAME(KQ),SNAME(KQ),UNITS(KQ),
     *     2,NMAX,FPOWER,xnwav,CLAT,CPRES)
  490 CONTINUE
  500 CONTINUE
      if(qdiag) call close_wp
      RETURN
C****
  901 FORMAT ('0',30X,A64,8X,'*1/60 (1/DAY)'/'   PERIOD EASTWARD--',
     *   35('---')/' N    -2      *-3   -3.3      -4       -5    -6   -7
     *.5  -10-12-15-20-30-60    60 30 20 15 12 10    7.5    6     5
     *   4*   VAR ERR'/'   --',40('---'))
  902 FORMAT (I2,41I3,I4,I4)
  903 FORMAT ('   --',40('---')/(1X,13F10.4))
  907 FORMAT ('1',A,I3,1X,A3,I5,' - ',I3,1X,A3,I5)
  911 FORMAT ('0',30X,A64,8X,'*1/60 (1/DAY)'/'   PERIOD EASTWARD--',
     *  35('---')/               ' N   *-4       -5    -6   -7.5  -10-12
     *-15-20-30-60    60 30 20 15 12 10    7.5    6     5        4
     * 3.3    3*       2    VAR ERR'/'   --',40('---'))
      END SUBROUTINE DIAG7P


      SUBROUTINE MEM (SERIES,ITM,MMAX,NUAMAX,NUBMAX,POWER,FPE,VAR,PNU)
      USE CONSTANT, only : pi
      IMPLICIT NONE
      DIMENSION C(1800),S(1800),B1(62),B2(62),A(12),AA(11),P(13)
      DIMENSION SERIES(*),POWER(*),FPE(*)
      REAL*8 ARG,PP,POWERX,P,C,S,POWER,FPE
      COMPLEX*16 CI,CSUM,A,AA,B1,B2,ANOM,ADEN
      COMPLEX*16 SERIES
      REAL*8 :: PNU,VAR
      INTEGER ::
     &     I,ITM,L,M,MMAX,MMAXP1,NU,NUA,
     &     NUAMAX,NUB,NUBMAX,NUMAX,NUTM
      CI=CMPLX(0.D0,1.D0)
      MMAXP1=MMAX+1
C**COSINE AND SINE FUNCTION
      NUMAX=NUAMAX*NUBMAX
      DO 20 NU=1,NUMAX
      ARG=2.0*PI*FLOAT(NU)/FLOAT(NUMAX)
      C(NU)=DCOS(ARG)
   20 S(NU)=DSIN(ARG)
   50 PP=0.0
      DO 60 I=1,ITM
   60 PP=PP+SERIES(I)*CONJG(SERIES(I))
      P(1)=PP/FLOAT(ITM)
      VAR=P(1)
      M=1
      B1(1)=SERIES(1)
      B2(ITM-1)=SERIES(ITM)
      DO 70 I=2,ITM-1
      B1(I)=SERIES(I)
   70 B2(I-1)=SERIES(I)
      GO TO 80
  100 DO 110 I=1,M
  110 AA(I)=A(I)
      M=M+1
      DO 120 I=1,ITM-M
      B1(I)=B1(I)-CONJG(AA(M-1))*B2(I)
  120 B2(I)=B2(I+1)-AA(M-1)*B1(I+1)
   80 ANOM=CMPLX(0.D0,0.D0)
      ADEN=CMPLX(0.D0,0.D0)
      DO 90 I=1,ITM-M
      ANOM=ANOM+CONJG(B1(I))*B2(I)
   90 ADEN=ADEN+B1(I)*CONJG(B1(I))+B2(I)*CONJG(B2(I))
      A(M)=(ANOM+ANOM)/ADEN
      P(M+1)=P(M)*(1.0-CONJG(A(M))*A(M))
      IF (M.EQ.1) GO TO 100
  130 CONTINUE
      DO 140 I=1,M-1
  140 A(I)=AA(I)-A(M)*CONJG(AA(M-I))
      IF (M.LT.MMAX) GO TO 100
C**FINAL PREDICTION ERROR
      DO 150 M=1,MMAXP1
  150 FPE(M)=P(M)*FLOAT(ITM+M-1)/FLOAT(ITM-M+1)
      DO 180 NUA=1,NUAMAX
      POWERX=0.
C**FREQUENCY BAND AVERAGE
      DO 170 NUB=1,NUBMAX
      NU=NUB+NUA*NUBMAX+(NUMAX-3*NUBMAX-1)/2
      CSUM=1.
      DO 160 M=1,MMAX
      NUTM=MOD(NU*M-1,NUMAX)+1
  160 CSUM=CSUM-A(M)*(C(NUTM)-CI*S(NUTM))
  170 POWERX=POWERX+P(MMAXP1)/(CSUM*CONJG(CSUM))
      POWER(NUA)=.5*POWERX/FLOAT(NUBMAX)
  180 CONTINUE
      PNU=0.0
      DO 210 L=1,NUAMAX
  210 PNU=PNU+POWER(L)
      PNU=PNU/(.5*NUAMAX)
      RETURN
      END SUBROUTINE MEM

      subroutine IJ_MAPk (k,smap,smapj,gm,igrid,jgrid,irange,
     &     name,lname,units)
!@sum IJ_MAPk returns the map data and related terms for the k-th field
      USE CONSTANT, only : grav,rgas,twopi,sha,kapa,bygrav,tf,undef
     *     ,teeny
      use model_com, only: modelEclock
      USE MODEL_COM, only : DTsrc,IDACC,
     &     JHOUR0,JDATE0,AMON,AMON0,JYEAR0,
     &     NDAY,Itime,Itime0,XLABEL,LRUNID
      USE DIAG_COM
      USE BDIJ
      USE MDIAG_COM, only : acc_period,
     &     sname_strlen,units_strlen,lname_strlen
      IMPLICIT NONE

      REAL*8, DIMENSION(IM,JM) :: anum,adenom,smap
      REAL*8, DIMENSION(JM) :: smapj
      integer, intent(in) :: k
      integer i,j,l,k1,k2,iwt,igrid,jgrid,irange,n1,n2
      character(len=sname_strlen) name
      character(len=units_strlen) units
      character(len=lname_strlen) lname
      real*8 :: gm,nh,sh, off, byiacc, scalek, an2Zan1
!@var  isumz,isumg = 1 or 2 if zon,glob sums or means are appropriate
      integer isumz,isumg
      integer :: year

      year = modelEclock%getYear()

      isumz = 2 ; isumg = 2  !  default: in most cases MEANS are needed
      if (k.eq.ij_dsev) isumz = 1

c**** Find & scale the numerators and find the appropriate denominators
c****
      name = name_ij(k) ; lname = lname_ij(k) ; units = units_ij(k)
      iwt = iw_all ; jgrid = jgrid_ij(k) ; irange = ir_ij(k)
      igrid = igrid_ij(k)

c**** the standard cases: aij(.,.,k) or aij(.,.,k)/aij(.,.,k1)
      byiacc = 1./(idacc(ia_ij(k))+teeny)
      do j=1,jm
      do i=1,im
        anum(i,j) = aij(i,j,k)*(scale_ij(k)*byiacc)
      end do
      end do

c**** ratios (the denominators)
      if (denom_ij(k) .gt. 0) then
        adenom(:,:) = aij(:,:,denom_ij(k))/
     &       (idacc( ia_ij(denom_ij(k)) )+teeny)
      else
        adenom = 1.             ! default
      end if

c**** remaining special cases for compound quantities

c**** precomputed fields: northward tranports by eddies
      if (k.eq.ij_ntdsese) then                   ! standing eddies
        byiacc=1./(idacc(ia_ij(ij_dsev))+teeny)   ; irange = ir_m95_265
        anum=SENTDSE*(byiacc*scale_ij(ij_dsev))  ;  igrid = 2; jgrid = 2
        isumz = 1 ; isumg = 2

      else if (k.eq.ij_ntdsete) then                  ! transient eddies
        byiacc=1./(idacc(ia_ij(ij_dsev))+teeny)   ; irange = ir_m1_3
        anum=TENTDSE*(byiacc*scale_ij(ij_dsev))  ;  igrid = 2; jgrid = 2
        isumz = 1 ; isumg = 2

      end if

c**** Find final field and zonal, global, and hemispheric means
      call ij_avg (anum,adenom,wt_ij(1,1,iwt),jgrid,isumz,isumg,  ! in
     *             smap,smapj,gm,nh,sh)                    ! out

c**** fill in some key numbers
      if (k .eq. IJ_RSIT) call keyij(gm,nh)
      if (k .eq. IJ_TGO) call keyij2(anum)

      return

      end subroutine ij_mapk


      subroutine ij_avg (anum,aden,wtij,jgrid,isumz,isumg,
     *                   smap,smapj,gm,nh,sh)
!@sum ij_avg finds num/den and various averages from num and den
!@auth R.Ruedy
      USE DOMAIN_DECOMP_ATM, only : GRID
      USE CONSTANT, only :  undef
      USE GEOM, only : wtj
      USE DIAG_COM, only : im,jm,fim
      IMPLICIT NONE

      real*8, dimension(im,jm) :: anum,aden,wtij,smap
      real*8, dimension(jm) :: smapj
      real*8, dimension(2) :: znumh,zdenh
      real*8  gm,nh,sh, sumj,wt

cBMP - added
      real*8, dimension(jm) :: sumjA, wtA
cBMP - added

      integer k,i,j,jgrid,isumz,isumg


c**** find smap,smapj  from the numerators and denominators
      smap = undef ; smapj = undef
      znumh = 0. ; zdenh = 0.

      If (Jgrid==1) then     ! Standard Grid
         DO J=1,JM
            sumj = 0. ; wt = 0.
            DO i=1,im
               sumj = sumj + anum(i,j)*wtij(i,j)
               wt   = wt   + aden(i,j)*wtij(i,j)
               if (aden(i,j)*wtij(i,j).ne.0.)
     *               smap(i,j)=anum(i,j)/aden(i,j)
            END DO
            if (isumz.eq.1) wt = 1.
            if (wt .gt. 0.) smapj(j) = sumj/wt
            sumjA(j) = sumj*wtj(j,isumg,jgrid)
            wtA(j)   = wt*wtj(j,isumg,jgrid)
         END DO
         CALL GLOBALSUM(GRID, sumjA(:), gm, znumh)
         CALL GLOBALSUM(GRID,   wtA(:), gm, zdenh)
      Else                   ! Staggered Grid
         DO J=2,JM
            sumj = 0. ; wt = 0.
            DO i=1,im
               sumj = sumj + anum(i,j)*wtij(i,j)
               wt   = wt   + aden(i,j)*wtij(i,j)
               if (aden(i,j)*wtij(i,j).ne.0.)
     *              smap(i,j)=anum(i,j)/aden(i,j)
            END DO
            if (isumz.eq.1) wt = 1.
            if (wt .gt. 0.) smapj(j) = sumj/wt
            sumjA(j) = sumj*wtj(j,isumg,jgrid)
            wtA(j)   = wt*wtj(j,isumg,jgrid)
         END DO
         CALL GLOBALSUM(GRID, sumjA(:), gm, znumh, istag=1)
         CALL GLOBALSUM(GRID,   wtA(:), gm, zdenh, istag=1)
      EndIf

c**** find hemispheric and global means
      nh = undef ; sh = undef ; gm = undef
      if (zdenh(1).gt.0.) sh = znumh(1)/zdenh(1)
      if (zdenh(2).gt.0.) nh = znumh(2)/zdenh(2)
      if (zdenh(1)+zdenh(2).gt.0.) gm = (znumh(1)+znumh(2))/
     /                                   (zdenh(1)+zdenh(2))
      if (isumg.eq.1) then
        sh = znumh(1) ; nh = znumh(2) ; gm = znumh(1)+znumh(2)
      end if

      return
      end subroutine ij_avg


      SUBROUTINE DIAGIJ
!@sum  DIAGIJ produces lat-lon fields as maplets (6/page) or full-page
!@+    digital maps, and binary (netcdf etc) files (if qdiag=true)
!@auth Gary Russell,Reto Ruedy
      USE CONSTANT, only : sha,teeny
      USE DOMAIN_DECOMP_ATM, only : GRID
      use model_com, only: modelEclock
      USE MODEL_COM, only :
     &     JHOUR0,JDATE0,AMON,AMON0,JYEAR0,
     &     NDAY,Itime,Itime0,XLABEL,LRUNID,idacc
      USE DYNAMICS, only : ido_gwdrag
      USE RAD_COM, only : cloud_rad_forc
      USE LAKES_COM, only : flake
      USE GEOM, only : DXV
      !USE VEG_COM, only : vdata
      USE DIAG_COM
      USE DIAG_COM_RAD
      USE GCDIAG
      USE GC_COM, only : AGC
      USE BDIJ
      USE MDIAG_COM, only : acc_period,
     &     sname_strlen,units_strlen,lname_strlen
      IMPLICIT NONE

!@var Qk: if Qk(k)=.true. field k still has to be processed
      logical, dimension (kaij) :: Qk
!@var Iord: Index array, fields are processed in order Iord(k), k=1,2,..
!@+     only important for fields 1->nmaplets+nmaps (appear in printout)
!@+     Iord(k)=0 indicates that a blank space replaces a maplet
      INTEGER Iord(kaij+10),nmaplets,nmaps ! 10 extra blank maplets
      INTEGER kmaplets
      REAL*8, DIMENSION(IM,JM) :: SMAP
      REAL*8, DIMENSION(JM) :: SMAPJ
      CHARACTER(len=sname_strlen) :: name
      CHARACTER(len=units_strlen) :: units
      CHARACTER(len=lname_strlen) :: lname
      CHARACTER xlb*32,title*48
!@var LINE virtual half page (with room for overstrikes)
      CHARACTER*133 LINE(53)
      logical qIij
      INTEGER ::   I,J,K,L,M,N,kcolmn,nlines,igrid,jgrid,irange,
     &     iu_Iij,koff

      REAL*8 ::
     &     DAYS,ZNDE16,DPTI,PVTI,gm,
     &     DE4TI,BYDPK,SZNDEG
      integer :: year, hour, date

      call modelEclock%get(year=year, hour=hour, date=date)

C**** OPEN PLOTTABLE OUTPUT FILE IF DESIRED
      IF(QDIAG) call open_ij(trim(acc_period)//'.ij'//XLABEL(1:LRUNID)
     *     ,im,jm)

C**** INITIALIZE CERTAIN QUANTITIES
C**** standard printout
      kmaplets = 57
      nmaplets = kmaplets + 6*isccp_diags + 2*cloud_rad_forc +        
     *           iDO_GWDRAG                  
      nmaps = 2
      iord(1:kmaplets) = (/
     *  ij_topo,    ij_fland,   ij_rsoi,     ! pg  1  row 1
     *  ij_rsnw,    ij_snow,    ij_rsit,     !        row 2
     *  ij_prec,    ij_evap,    ij_shdt,     ! pg  2  row 1
     *  ij_beta,    ij_rune,    ij_tg1,      !        row 2
     *  ij_ws,      ij_jet ,    ij_dtdp,     ! pg  3  row 1
     *  ij_wsdir,   ij_jetdir,  ij_sstabx,   !        row 2
     *  ij_netrdp,  ij_srnfp0,  ij_btmpw,    ! pg  4  row 1
     *  ij_srtr,    ij_srincg,  ij_clr_srincg, !      row 2
     *  ij_albp,    ij_albv,    ij_trnfp0,   ! pg  5  row 1
     *  ij_albg,    ij_albgv,   ij_neth,     !        row 2
     *  ij_dsev,    ij_ntdsese, ij_ntdsete,  ! pg  6  row 1
     *  ij_gwtr,    ij_wmsum,   ij_qm,       !        row 2
     *  ij_cldcv,   ij_dcnvfrq, ij_scnvfrq,  ! pg  7  row 1
     *  ij_pmccld,  ij_pdcld,   ij_pscld,    !        row 2
     *  ij_wtrcld,  ij_optdw,   ij_cldtppr,  ! pg  8  row 1
     *  ij_icecld,  ij_optdi,   ij_cldtpt,   !        row 2
     *  ij_cldcv1,  ij_cldt1p,  ij_cldt1t,   ! pg  9  row 1
     *  ij_pcldl,   ij_pcldm,   ij_pcldh,    !        row 2
     *  ij_pblht,   ij_gusti,   ij_mccon /)  ! pg  10 row 1

C**** include ISCCP diags if requested
      if (isccp_diags.eq.1) then
        iord(kmaplets+1:kmaplets+6) = (/ij_taui,ij_ctpi,ij_lcldi,
     *                                  ij_mcldi,ij_hcldi,ij_tcldi/)
        kmaplets=kmaplets+6
      else
        lname_ij(ij_lcldi)='unused'
        lname_ij(ij_mcldi)='unused'
        lname_ij(ij_hcldi)='unused'
        lname_ij(ij_tcldi)='unused'
        lname_ij(ij_taui) ='unused'
        lname_ij(ij_ctpi) ='unused'
      end if

C**** include CRF diags if requested
      if (cloud_rad_forc.eq.1) then
        iord(kmaplets+1:kmaplets+2) = (/ij_swcrf,ij_lwcrf/)
        kmaplets=kmaplets+2
      else
        lname_ij(ij_swcrf)='unused'
        lname_ij(ij_lwcrf)='unused'
      end if

C**** Fill in maplet indices for gravity wave diagnostics
      do k=1,iDO_GWDRAG
        iord(k+kmaplets) = ij_gw1+k-1  !i.e. first entry is ij_gw1
      end do

C**** Add the full-page maps (nmaps)
      iord(nmaplets+1:nmaplets+nmaps) = (/ij_slp,ij_ts/)
c**** always skip unused fields
      Qk = .true.
      do k=1,kaij
        if (index(lname_ij(k),'unused').gt.0) Qk(k) = .false.
      end do

      inquire (file='Iij',exist=qIij)
      if (.not.qIij) kdiag(3)=0
      call set_ijout (nmaplets,nmaps,Iord,Qk,iu_Iij)
      xlb=acc_period(1:3)//' '//acc_period(4:12)//' '//XLABEL(1:LRUNID)
C****
      DAYS=(Itime-Itime0)/FLOAT(nday)
C**** Collect the appropriate weight-arrays in WT_IJ
      do J=1,JM
      do i=1,im
        wt_ij(i,j,iw_all) = 1.
cgsfc        wt_ij(i,j,2) = focean(i,j)
        wt_ij(i,j,iw_ocn) = focean_glob(i,j)
cgsfc        wt_ij(i,j,3) = flake(i,j)
        wt_ij(i,j,iw_lake) = aij(i,j,ij_lk)/(idacc(ia_ij(ij_lk))+teeny)
cgsfc        wt_ij(i,j,4) = flice(i,j)
        wt_ij(i,j,iw_lice) = flice_glob(i,j)
cgsfc        wt_ij(i,j,5) = fearth(i,j)
        wt_ij(i,j,iw_soil) = 1.d0 - wt_ij(i,j,iw_ocn)
     &       - wt_ij(i,j,iw_lake) -  wt_ij(i,j,iw_lice)
cgsfc        wt_ij(i,j,6) = fearth(i,j)*(vdata(i,j,1)+vdata(i,j,10))
        wt_ij(i,j,iw_bare) = wt_ij(i,j,iw_soil)
     &       *(1.d0 - aij(i,j,ij_fveg)/(idacc(ia_ij(ij_fveg))+teeny))
cgsfc        wt_ij(i,j,7) = fearth(i,j)*(1.-(vdata(i,j,1)+vdata(i,j,10)))
        wt_ij(i,j,iw_veg) = wt_ij(i,j,iw_soil)
     &       *aij(i,j,ij_fveg)/(idacc(ia_ij(ij_fveg))+teeny)
        wt_ij(i,j,iw_land) = wt_ij(i,j,iw_soil) + wt_ij(i,j,iw_lice)
      end do
      end do
C**** CACULATE STANDING AND TRANSIENT EDDY NORTHWARD TRANSPORT OF DSE
      SENTDSE = 0
      TENTDSE = 0
      DO J=2,JM
      DO K=1,LM
        DPTI=0.
        PVTI=0.
        DE4TI=0.
        DO I=1,IM
          IF (AIJK(I,J,K,IJK_DPB).GT.0.) THEN
            DPTI=DPTI+AIJK(I,J,K,IJK_DPB)
            BYDPK=1./(AIJK(I,J,K,IJK_DPB)+teeny)
            PVTI=PVTI+AIJK(I,J,K,IJK_VB)
            DE4TI=DE4TI+AIJK(I,J,K,IJK_DSE)
            SENTDSE(I,J)=SENTDSE(I,J)
     *        +(AIJK(I,J,K,IJK_DSE)*AIJK(I,J,K,IJK_VB)*BYDPK)
          END IF
        END DO
        SZNDEG=DE4TI*PVTI/(DPTI+teeny)
        DO I=1,IM
          SENTDSE(I,J)=SENTDSE(I,J)-SZNDEG*byim
        END DO
      END DO
      END DO
      DO J=2,JM
        ZNDE16=0.
        DO L=1,LM
c          ZNDE16=ZNDE16+(SHA*AGC(J,L,JK_ZMFNTSH)+AGC(J,L,JK_ZMFNTGEO))
          ZNDE16=ZNDE16+AGC(J,L,JK_TOTNTDSE)-
     &         SHA*AGC(J,L,JK_EDDNTSH)-AGC(J,L,JK_EDDNTGEO)
        END DO
        ZNDE16=ZNDE16*DXV(J)*byim
        DO I=1,IM
          SENTDSE(I,J)=SENTDSE(I,J)*DXV(J)
          TENTDSE(I,J)=AIJ(I,J,IJ_DSEV)-ZNDE16-SENTDSE(I,J)
        END DO
      END DO

C**** Fill in the undefined pole box duplicates
      DO N=1,KAIJ
      IF (JGRID_ij(N).EQ.2) CYCLE
      DO I=1,IM
         AIJ(I,1,N)=AIJ(1,1,N)
      END DO
      DO I=1,IM
         AIJ(I,JM,N)=AIJ(1,JM,N)
      END DO
      END DO

C**** Print out 6-map pages
      do n=1,nmaplets
        if (mod(n-1,6) .eq. 0) then
c**** print header lines
          WRITE (6,901) XLABEL
          WRITE (6,902) JYEAR0,AMON0,JDATE0,JHOUR0,
     *      YEAR,AMON,DATE,HOUR,ITIME,DAYS
        end if
        kcolmn = 1 + mod(n-1,3)
        if (kcolmn .eq. 1) line=' '
c**** Find, then display the appropriate array
        k = Iord(n)
        if (k .gt. 0 .and. Qk(k)) then
          call ij_mapk (k,smap,smapj,gm,igrid,jgrid,irange,name,lname,
     &          units)
          title=trim(lname)//' ('//trim(units)//')'
          call maptxt(smap,smapj,gm,irange,title,line,kcolmn,nlines)
          if(qdiag) call pout_ij(title//xlb,name,lname,units,
     *                            smap,smapj,gm,igrid,jgrid)
          Qk(k) = .false.
        end if
c**** copy virtual half-page to paper if appropriate
        if (kcolmn.eq.3 .or. n.eq.nmaplets) then
          do k=1,nlines
            write (6,'(a133)') line(k)
          end do
        end if
      end do

C**** Print out full-page digital maps
      do n=nmaplets+1,nmaplets+nmaps
        k = Iord(n)
      if (k.le.0 .or. .not.Qk(k)) cycle
        call ij_mapk (k,smap,smapj,gm,igrid,jgrid,irange,name,lname,
     &     units)
        title=trim(lname)//' ('//trim(units)//')'
        call ijmap (title//xlb,smap,smapj,jgrid)
        if(qdiag) call pout_ij(title//xlb,name,lname,units,smap,smapj,
     *                          gm,igrid,jgrid)
        Qk(k) = .false.
      end do

C**** produce binary files of remaining fields if appropriate
      do k=1,kaij
        if (Qk(k)) then
          call ij_mapk (k,smap,smapj,gm,igrid,jgrid,irange,name,lname,
     &          units)
          title=trim(lname)//' ('//trim(units)//')'
          if (qdiag) call pout_ij(title//xlb,name,lname,units,smap,smapj
     *         ,gm,igrid,jgrid)
        end if
      end do
      if (qdiag) then
        call close_ij
        if (kdiag(3).lt.8) CALL IJKMAP (iu_Iij)
        if (kdiag(3).lt.8) CALL IJLMAP (iu_Iij)
      end if

      RETURN
C****
  901 FORMAT ('1',A)
  902 FORMAT ('0',15X,'From:',I6,A6,I2,',  Hr',I3,
     *  6X,'To:',I6,A6,I2,', Hr',I3,'  Model-Time:',I9,5X,
     *  'Dif:',F7.2,' Days')
      END SUBROUTINE DIAGIJ

      subroutine maptxt (smap,smapj,gm,irange,title,line,kcol,nlines)
!@sum  maptxt prints a maplet onto 1/3 of a virtual half-page line(1-51)
!@auth R.Ruedy
      use constant, only : undef
      use diag_com, only : im,jm,inci,incj
      use bdij

      IMPLICIT NONE

      real*8, dimension(im,jm) :: smap
      real*8, dimension(jm) :: smapj
      real*8  gm,val,zmax,zmin
      integer irange,kcol,k,k1,ifm,nz1,nz2,n1,n2,ibar,i,j
     *  ,nlines
      character(len=133), dimension(53) :: line
      character(len=40) :: title

c**** find first and last column for zonal means nz1,nz2 and maps n1,n2
      nz1 = 2 + (kcol-1)*(9+im/inci) ; nz2 = nz1 + 4
      n1  = nz2 + 2  ;  n2 = n1 + im/inci - 1
c**** pick color bar and format for zonal mean according to the range
      ibar = ib_of_legnd(irange)
      zmax = 999.9 ; zmin = -99.9
      if (kcol.gt.1) then
        zmax = 99999.9 ; zmin = -9999.9
        nz1 = nz1 - 2
      end if
      ifm = 2
      do j=1,jm
        if (smapj(j).eq.undef) cycle
        if (smapj(j).gt.zmax .or. smapj(j).lt.zmin) ifm = 1
      end do

c**** title on line 1
      line(1)(n1-4:n2) = title(1:40) ; line(1)(1:1) = '0'  ! line feed
c**** use 2 lines for each shown latitude because of overstrike
      k=0
      do j=jm,1,-incj
        k = k+2
c**** zonal mean
        val=0.
        if (smapj(j) .ne. undef) val = smapj(j)
        if (kcol.eq.1) then
          if (ifm.eq.1) write(line(k)(nz1:nz2),'(i5)') nint(min(1d5,max(
     *         -1d5,val)))
          if (ifm.eq.2) write(line(k)(nz1:nz2),'(f5.1)') val
        else
          if (ifm.eq.1) write(line(k)(nz1:nz2),'(i7)') nint(min(1d5,max(
     *         -1d5,val)))
          if (ifm.eq.2) write(line(k)(nz1:nz2),'(f7.1)') val
        end if
c**** mark selected longitudes for each selected latitude
        k1 = n1-1
        do i=1,im,inci
          k1 = k1+1
          val = undef
          if (smap(i,j).ne.undef) val = smap(i,j)*fac_legnd(irange)
c**** for angles, change range from -180->180 to 0-360 (fac=.1)
          if (irange.eq.ir_angl .and. val.lt.-.5) val = val+36.
          if (irange.eq.ir_angl .and. val.le..1) val = .1
          line(k)(k1:k1) = mark(val,ibar,undef)
          if (wt_ij(i,j,iw_land).gt..5) line(k+1)(k1:k1)=line(k)(k1:k1)
          line(k+1)(1:1) = '+'             !  overstrike
        end do
      end do
c**** below map, show global mean and mark each quarter with a '+'
      k = k+2
      if (gm .ge. zmin .and. gm .le. zmax) then
        if (kcol.eq.1) write(line(k)(nz1:nz2),'(f5.1)') gm
        if (kcol.gt.1) write(line(k)(nz1:nz2),'(f7.1)') gm
      else
        if (gm.eq.undef) then
          if (kcol.eq.1) write(line(k)(nz1:nz2),'(a)') "Undef"
          if (kcol.gt.1) write(line(k)(nz1:nz2),'(a)') "  Undef"
        else
          if (kcol.eq.1) write(line(k)(nz1:nz2),'(i5)') nint(min(1d5
     *         ,max(-1d5,gm)))
          if (kcol.gt.1) write(line(k)(nz1:nz2),'(i7)') nint(min(1d5
     *         ,max(-1d5,gm)))
        end if
      end if
      do k1 = n1,n2,im/(4*inci)
        line(k)(k1:k1) = '+'
      end do
c**** last line: legend (there is a little more space in col 1 and 2)
      if(kcol.lt.3) n2 = n1 + 39
      line(k+1)(n1:n2) = legend(irange)
      nlines = k+1

      return
      end subroutine maptxt


      SUBROUTINE IJMAP (title,smap,smapj,jgrid)
C**** Print out full-page digital maps
      USE CONSTANT, only :  undef
      use model_com, only: modelEclock
      USE MODEL_COM, only :
     &     NDAY,JHOUR0,JDATE0,AMON,AMON0,
     &     JYEAR0,Itime,Itime0,XLABEL,lrunid
      USE GEOM, only :
     &     LAT_DG,LON_DG
      use diag_com, only : im,jm,inc=>inci,wt_ij,iw_land
      IMPLICIT NONE

      CHARACTER*48 TITLE

      CHARACTER(LEN=3), DIMENSION(IM) :: LINE
      CHARACTER(LEN=9) :: AVG

      REAL*8, DIMENSION(IM,JM) :: SMAP
      REAL*8, DIMENSION(JM) :: SMAPJ
      REAL*8 :: DAYS
      INTEGER :: I,J,jgrid
      integer :: year, hour, date

      call modelEclock%get(year=year, hour=hour, date=date)

C**** WRITE HEADER LINES
      DAYS=(Itime-Itime0)/FLOAT(nday)
      WRITE(6,901)XLABEL
      WRITE (6,902) JYEAR0,AMON0,JDATE0,JHOUR0,
     *  YEAR,AMON,DATE,HOUR,ITIME,DAYS
      WRITE(6,900) TITLE(1:48)
      DO I=1,IM
        WRITE(LINE(I),'(I3)') I
      end do
      AVG='     MEAN'
      WRITE (6,910) (LINE(I),I=1,IM,INC),AVG
      WRITE(6,940)
      WRITE(6,940)

C**** PRINT MAP
      DO J=JM,jgrid,-1
        do i=1,im
          IF (SMAP(I,J).LT.999.5.AND.SMAP(I,J).GE.-99.5) then
            line(i) = '   '
            write (line(i),'(i3)') nint(SMAP(I,J))
          else if (SMAP(I,J).eq.undef) then
            line(i) = '   '
          else
            line(i) = ' **'
          end if
        end do
        WRITE(AVG,'(a1,F8.2)') ' ',SMAPJ(J)
        if (SMAPJ(J).eq.undef) AVG='         '
        WRITE (6,920) NINT(LAT_DG(J,jgrid)),J,(LINE(I),I=1,IM,INC),AVG
        DO I=1,IM
          IF (wt_ij(i,j,iw_land).lt..5) LINE(I)='   '
        end do
        WRITE (6,925) (LINE(I),I=1,IM,INC)
        WRITE (6,925) (LINE(I),I=1,IM,INC)
        IF (JM.LE.24) WRITE (6,940)
      END DO
      WRITE (6,930) (LON_DG(I,jgrid),I=1,IM,INC*2)
      RETURN
C****
  900 FORMAT('0',45X,A48)
  901 FORMAT ('1',A)
  902 FORMAT ('0',15X,'From:',I6,A6,I2,',  Hr',I3,
     *  6X,'To:',I6,A6,I2,', Hr',I3,'  Model-Time:',I9,5X,
     *  'Dif:',F7.2,' Days')
  910 FORMAT('0LAT  J/I  ',36A3,A9)
  920 FORMAT(2I4,3X,36A3,A9)
  925 FORMAT('+',10X,36A3)
  930 FORMAT('0  LONG',2X,19F6.1)
  940 FORMAT(' ')
      END SUBROUTINE IJMAP


      subroutine set_ijout (nmaplets,nmaps,Iord,Qk,iu_Iij)
!@sum set_ijout either lists or sets the fields to be processed
!@auth Reto A. Ruedy
      USE DIAG_COM
      USE BDIJ
      use filemanager

      IMPLICIT NONE
      character*80 line
      logical Qk(kaij),Qktmp(kaij)
      INTEGER Iord(kaij+10),nmaplets,nmaps,iu_Iij,k,
     *   n,kmap(3)

c**** Just list what's available - then do same for ijk-fields
      if (kdiag(3) .eq. 0) then
        Qktmp = Qk
        call openunit('Iij',iu_Iij,.false.,.false.)
        write (iu_Iij,'(a)') 'List of fields shown as maplets'
        do n=1,nmaplets
          k = Iord(n)
          Qktmp(k) = .false.
          if (k.le.0) then
             write (iu_Iij,'(i3,1x,a)') k, '  blank maplet'
          else
             write (iu_Iij,'(i3,1x,a)') k,lname_ij(k)
          end if
        end do
        write (iu_Iij,'(a)') 'List of fields shown as 1-pg maps'
        do n=nmaplets+1,nmaplets+nmaps
          k = Iord(n)
          Qktmp(k) = .false.
          if (k.le.0) then
             cycle
          else
             write (iu_Iij,'(i3,1x,a)') k,lname_ij(k)
          end if
        end do
        write (iu_Iij,'(a)') 'List of other fields in binary output'
        do k=1,kaij
          if (.not.Qktmp(k)) cycle
          write (iu_Iij,'(i3,1x,a)') k,lname_ij(k)
        end do
        kdiag(3)=9
        CALL IJKMAP (iu_Iij)
        kdiag(3)=0
        return
      end if

c**** Redefine nmaplets,nmaps,Iord,Qk if  kdiag(3) > 0
      call openunit('Iij',iu_Iij,.false.,.true.)

      nmaplets = 0 ; nmaps = 0 ; Iord = 0 ; Qk = .false.

      kmap = 0 ; n=0 ; k=0
   10 read (iu_Iij,'(a)',end=20) line
      if (line(1:1) .eq. 'l') go to 20
      if (line(1:1) .eq. 'L') then
        n=n+1
        go to 10
      end if
      k = k+1
      read(line,'(i3)') Iord(k)
      if (Iord(k).gt.0) Qk(Iord(k)) = .true.
      kmap(n) = kmap(n) + 1
      go to 10

   20 nmaplets = kmap(1) ; nmaps = kmap(2)
      if (.not.qdiag .or. kdiag(3).eq.1) call closeunit(iu_Iij)
      return
      end subroutine set_ijout

      SUBROUTINE DIAGCP
!@sum  DIAGCP produces tables of the conservation diagnostics
!@auth Gary Russell/Gavin Schmidt
      USE DOMAIN_DECOMP_ATM, only : GRID
      use model_com, only: modelEclock
      USE MODEL_COM, only :
     &     idacc,jhour0,jdate0,amon,amon0,
     &     jyear0,nday,itime,itime0,xlabel,lrunid
      USE GEOM, only :
     &     areag,WTJ
      USE DIAG_COM, only :  fim,qdiag,
     &     consrv,kcon,scale_con,title_con,nsum_con,ia_con,kcmx,
     *     inc=>incj,ia_inst,name_consrv,lname_consrv,
     *     units_consrv,jm=>jm_budg,dxyp_budg,lat_budg
      USE GC_COM, only : jeq
      USE DIAG_ZONAL, only : xwon
      USE MDIAG_COM, only : acc_period,
     &     sname_strlen,units_strlen,lname_strlen
      IMPLICIT NONE

      REAL*8, DIMENSION(JM,KCMX) :: CSJ
      INTEGER, DIMENSION(JM) :: MAREA
      REAL*8, DIMENSION(KCON) :: FGLOB
      REAL*8, DIMENSION(2,KCON) :: FHEM
      INTEGER, DIMENSION(JM,KCON) :: MLAT
      REAL*8, DIMENSION(JM+3,KCON) :: CNSLAT
      CHARACTER*4, PARAMETER :: HEMIS(2) = (/' SH ',' NH '/),
     *     DASH = ('----')

      INTEGER :: j,jhemi,jnh,jp1,jpm,jv1,jvm,jx,n,n1
      REAL*8 :: aglob,ahem,days
C**** Arrays needed for full output and pdE
      CHARACTER*38, DIMENSION(KCON) :: TITLEO
      integer :: year, hour, date

      call modelEclock%get(year=year, hour=hour, date=date)

C**** OPEN PLOTTABLE OUTPUT FILE IF DESIRED
      IF (QDIAG)
     *     call open_jc(trim(acc_period)//'.jc'//XLABEL(1:LRUNID),
     *     jm,lat_budg)

C**** CALCULATE SCALING FACTORS
      IF (IDACC(ia_inst).LT.1) IDACC(ia_inst)=1
C**** CALCULATE SUMMED QUANTITIES
C**** LOOP BACKWARDS SO THAT INITIALISATION IS DONE BEFORE SUMMATION!
      DO J=1,JM
        DO N=KCMX,1,-1
          IF (NSUM_CON(N).eq.0) THEN
            CONSRV(J,N)=0.
          ELSEIF(IDACC(IA_CON(N)).eq.0) then
            ! may happen for ia_con == ia_12hr when starting
            ! at a time not divisible by 12 hrs
            CONSRV(J,N)=0.
          ELSEIF (NSUM_CON(N).gt.0) THEN
            CONSRV(J,NSUM_CON(N))=CONSRV(J,NSUM_CON(N))+CONSRV(J,N)
     *           *SCALE_CON(N)*IDACC(ia_inst)/(IDACC(IA_CON(N))+1d-20)
          END IF
        END DO
      END DO

C**** CALCULATE ALL CONSERVED QUANTITIES ON TRACER GRID
      cnslat = 0.
      csj = 0.
      N1=1
      DO N=N1,KCMX
        if(IDACC(IA_CON(N)).eq.0) then
          ! may happen for ia_con == ia_12hr when starting
          ! at a time not divisible by 12 hrs
          cycle
        endif
        DO J=1,JM
          CSJ(J,N)    = CONSRV(J,N)*SCALE_CON(N)/
     &                           (IDACC(IA_CON(N))+1d-20)
          CNSLAT(J,N) = CSJ(J,N)
          CSJ(J,N)    = CSJ(J,N)*DXYP_BUDG(J)
        END DO
      END DO
      CALL GLOBALSUM(GRID, CSJ(:,N1:KCMX),
     &                     FGLOB(N1:KCMX), FHEM(:,N1:KCMX))
      FGLOB(N1:KCMX)=FGLOB(N1:KCMX)/AREAG
      FHEM(1,N1:KCMX)=FHEM(1,N1:KCMX)/(.5*AREAG)
      FHEM(2,N1:KCMX)=FHEM(2,N1:KCMX)/(.5*AREAG)

      AGLOB=1.D-10*AREAG*XWON
      AHEM=1.D-10*(.5*AREAG)*XWON
C**** LOOP OVER HEMISPHERES
      DAYS=(Itime-Itime0)/FLOAT(nday)
      DO N=1,KCMX
        DO J=1,JM
          MLAT(J,N)=NINT(CNSLAT(J,N))
        END DO
        CNSLAT(JM+1,N)=FHEM(1,N)
        CNSLAT(JM+2,N)=FHEM(2,N)
        CNSLAT(JM+3,N)=FGLOB(N)
          titleo(n)=title_con(n)
      END DO
      DO JHEMI=2,1,-1
        WRITE (6,901) XLABEL
        WRITE (6,902) JYEAR0,AMON0,JDATE0,JHOUR0,
     *       YEAR,AMON,DATE,HOUR,ITIME,DAYS
        JP1=1+(JHEMI-1)*(JEQ-1)
        JPM=JHEMI*(JEQ-1)

C**** WRITE TABLES
        DO J=JP1,JPM
          MAREA(J)=1.D-10*XWON*DXYP_BUDG(J)+.5
        END DO
        DO N=1,KCMX
          IF(N.EQ.1 .OR. N.EQ.26) THEN ! AM/KE get their own section
            WRITE(6,907)
            WRITE(6,903) (DASH,J=JP1,JPM,INC)
            WRITE(6,904) HEMIS(JHEMI),
     &           (NINT(LAT_BUDG(JX)),JX=JPM,JP1,-INC)
            WRITE(6,903) (DASH,J=JP1,JPM,INC)
          ENDIF
          WRITE (6,905) TITLE_CON(N),FGLOB(N),FHEM(JHEMI,N),
     *         (MLAT(JX,N),JX=JPM,JP1,-INC)
          IF(N.EQ.25 .OR. N.EQ.KCMX) THEN ! AM/KE get their own section
            WRITE (6,906) AGLOB,AHEM,(MAREA(JX),JX=JPM,JP1,-INC)
          ENDIF
        END DO

      END DO
      IF (QDIAG) CALL POUT_JC(titleo,name_consrv,lname_consrv,
     *     units_consrv,cnslat,kcmx)

      if(qdiag) call close_JC
      RETURN
C****
  901 FORMAT ('1',A)
  902 FORMAT ('0Conservation Quantities       From:',
     *  I6,A6,I2,',  Hr',I3,  6X,  'To:',I6,A6,I2,', Hr',I3,
     *  '  Model-Time:',I9,5X,'Dif:',F7.2,' Days')
  903 FORMAT (1X,25('--'),13(A4,'--'))
  904 FORMAT (35X,'GLOBAL',A7,2X,13I6)
  905 FORMAT (A32,2F9.2,1X,13I6)
  906 FORMAT ('0AREA (10**10 m^2)',14X,2F9.1,1X,13I6)
  907 FORMAT ('0')
      END SUBROUTINE DIAGCP

      SUBROUTINE DIAG5P
!@sum  DIAG5P PRINTS THE SPECTRAL ANALYSIS TABLES
!@auth Gary Russell
      USE CONSTANT, only : grav,rgas,teeny
      use model_com, only: modelEclock
      USE MODEL_COM, only :
     &     IDACC,JHOUR0,JDATE0,
     &     AMON,AMON0,JYEAR0,XLABEL
      USE GEOM, only : DXYV
      USE DIAG_COM, only : im,jm,lm,fim,
     &     aijk,ia_d5s,ia_filt,ia_12hr,ia_d5f,ia_d5d,ia_dga
     *     ,ia_inst,kdiag
      USE GC_COM, only : speca,atpe,agc,kspeca,ktpe,nhemi,nspher,klayer,
     &     istrat,jeq
      USE DIAG_ZONAL, only : xwon
      USE GCDIAG
      USE DYNAMICS, only : DT
      IMPLICIT NONE

      REAL*8, DIMENSION(IM) :: X
      REAL*8, DIMENSION(KSPECA) :: SCALET,F0,FNSUM
      INTEGER, DIMENSION(KSPECA) :: MN

      REAL*8, DIMENSION(KTPE,NHEMI) :: FATPE

      INTEGER, PARAMETER :: IZERO=0

      INTEGER, DIMENSION(KTPE), PARAMETER ::
     &     MAPEOF=(/3,8,10,11,13,15,17,20/)

      CHARACTER*8 :: LATITD(4) = (/
     *     'SOUTHERN','NORTHERN',' EQUATOR','45 NORTH'/)
      CHARACTER*16 :: SPHERE(4)=
     *     (/'TROPOSPHERE     ','LOW STRATOSPHERE',
     *       'MID STRATOSPHERE','UPP STRATOSPHERE'/)
      REAL*8,DIMENSION(4) :: SCALEK=(/1.,1.,10.,10./)

      INTEGER ::
     &     I,IUNITJ,IUNITW,J,J45N,
     &     K,KPAGE,KROW,KSPHER,L,
     &     M,MAPE,MTPE,N,NM,NM1, LATF,LATL

      REAL*8 :: FACTOR,FNM
      integer :: year, hour, date

      call modelEclock%get(year=year, hour=hour, date=date)

      NM=1+IM/2
      IF (IDACC(ia_inst).LT.1) IDACC(ia_inst)=1
      J45N=2.+.75*(JM-1.)
C****
C**** STANDING KINETIC ENERGY
C****
      DO 710 K=1,NSPHER
      DO 710 N=1,NM
  710 SPECA(N,1,K)=0.
      DO 770 L=1,LM
        KSPHER=KLAYER(L)
      DO 770 J=2,JM
      IF (AGC(J,L,JK_DPB).LE.teeny) GO TO 770
      FACTOR=FIM*DXYV(J)/AGC(J,L,JK_DPB)
      DO 769 K=0,1
      DO 720 I=1,IM
  720 X(I)=AIJK(I,J,L,IJK_UB+K)
      CALL FFTE (X,X)
      IF (J.EQ.JEQ) GO TO 750
      DO 730 N=1,NM
  730 SPECA(N,1,KSPHER)=SPECA(N,1,KSPHER)+X(N)*FACTOR
      IF (J.NE.J45N) GO TO 769
      DO 740 N=1,NM
  740 SPECA(N,1,KSPHER+2)=SPECA(N,1,KSPHER+2)+X(N)*FACTOR
      GO TO 769
  750 DO 760 N=1,NM
      SPECA(N,1,KSPHER+2)=SPECA(N,1,KSPHER+2)+X(N)*FACTOR
      SPECA(N,1,KSPHER)=SPECA(N,1,KSPHER)+.5*X(N)*FACTOR
  760 SPECA(N,1,KSPHER+1)=SPECA(N,1,KSPHER+1)+.5*X(N)*FACTOR
      KSPHER=KSPHER+K
  769 CONTINUE
  770 CONTINUE
C****
  600 SCALET(1)=100.D-17/(GRAV*IDACC(ia_dga)+teeny)
      SCALET(19)=100.D-17/(GRAV*IDACC(ia_inst))
      SCALET(20)=SCALET(19)*RGAS
      SCALET(2)=SCALET(19)*IDACC(ia_inst)/(IDACC(ia_d5d)+teeny)
      SCALET(3)=SCALET(2)*RGAS
      SCALET(4)=100.D-12/(GRAV*DT*IDACC(ia_d5f)+teeny)
      SCALET(5)=SCALET(4)
      SCALET(6)=SCALET(4)
      SCALET(7)=100.D-12/(GRAV*DT*(IDACC(ia_d5d)+teeny))
      SCALET(8)=SCALET(7)*RGAS
      SCALET(9)=100.D-12/(GRAV*DT*(IDACC(ia_d5s)+teeny))
      SCALET(10)=SCALET(9)*RGAS
      SCALET(11)=SCALET(10)
      SCALET(12)=SCALET(9)
      SCALET(13)=SCALET(10)
      SCALET(14)=100.D-12/(GRAV*DT*(IDACC(ia_filt)+teeny))
      SCALET(15)=SCALET(14)*RGAS
      if(IDACC(ia_12hr).eq.0) then
        ! may happen when starting at a time not divisible by 12 hrs
        SCALET(16)=0.
      else
        SCALET(16)=100.D-12/(GRAV*DT*(.5*IDACC(ia_12hr)+teeny))
      endif
      SCALET(17)=SCALET(16)*RGAS
      SCALET(18)=100.D-17/(GRAV*IDACC(ia_dga)+teeny)
      DO 605 K=1,KSPECA
  605 SCALET(K)=XWON*SCALET(K)
      IUNITJ=17
      IUNITW=12
      LATF=1
      LATL=4
      IF (KDIAG(5).GT.0) LATL=4-KDIAG(5)  ! just zones 1->(4-kd5)
      IF (KDIAG(5).LT.0.AND.KDIAG(5).GT.-5) LATF=-KDIAG(5)
      IF (KDIAG(5).LT.0) LATL=LATF        ! just 1 zone
      DO 690 KPAGE=LATF,LATL   ! one for each lat.zone SH/NH/EQ/45N
C**** WRITE HEADINGS
      WRITE (6,901) XLABEL
      WRITE (6,902) JYEAR0,AMON0,JDATE0,JHOUR0,YEAR,AMON,DATE,
     *  IUNITJ,IUNITW
      DO 670 KROW=1,2+ISTRAT !one for each level (trp/lstr/mstr/ustr)
      IF (JM.GE.25.AND.KROW.EQ.2) WRITE (6,901)
      WRITE (6,903) LATITD(KPAGE),SPHERE(KROW)
      KSPHER=4*(KROW-1)+KPAGE
C**** WRITE KINETIC AND AVAILABLE POTENTIAL ENERGY BY WAVE NUMBER
      DO 610 M=1,KSPECA
      F0(M)=SPECA(1,M,KSPHER)*SCALET(M)*SCALEK(KROW)
      MN(M)=NINT(F0(M))
  610 FNSUM(M)=0.
      WRITE (6,904) MN
      DO 630 N=2,NM
      KSPHER=4*(KROW-1)+KPAGE
      DO 620 M=1,KSPECA
      FNM=SPECA(N,M,KSPHER)*SCALET(M)*SCALEK(KROW)
      MN(M)=NINT(FNM)
  620 FNSUM(M)=FNSUM(M)+FNM
      NM1=N-1
  630 WRITE (6,905) NM1,MN
      DO 640 M=1,KSPECA
  640 MN(M)=NINT(FNSUM(M))
      WRITE (6,906) MN
      DO 650 M=1,KSPECA
  650 MN(M)=NINT(FNSUM(M)+F0(M))
      WRITE (6,907) MN
  670 CONTINUE
      IF (KPAGE.GE.3) GO TO 690
C**** WRITE TOTAL POTENTIAL ENERGY
      DO 680 MTPE=1,KTPE
      MAPE=MAPEOF(MTPE)
         FATPE(MTPE,KPAGE)=ATPE(MTPE,KPAGE)*SCALET(MAPE)/RGAS
  680 MN(MTPE)=NINT(FATPE(MTPE,KPAGE))
      WRITE (6,909) (MN(MTPE),MTPE=1,8)
      IF (KPAGE.NE.2) GO TO 690
      DO 685 M=1,KSPECA
  685 SCALET(M)=SCALET(M)*10.
      IUNITJ=16
      IUNITW=11
  690 CONTINUE
      RETURN
C****
  901 FORMAT ('1',A)
  902 FORMAT ('0**  Spectral Analysis **      From:',
     *  I6,A6,I2,',  Hr',I3,  6X,  'To:',I6,A6,I2,', Hr',I3,
     *  '       UNITS 10**',I2,' JOULES AND 10**',I2,' WATTS')
  903 FORMAT ('0',50X,A8,1X,A16/
     *  13X,'MEAN',19X,'DYNAMICS',25X,'SOURCES',16X,'FILTER',8X,
     *     'DAILY',4X,'PR SURF',5X,'LAST'/
     *'   N    SKE   KE   APE    KADV  KCOR   P-K  KDYN  PDYN   ',
     *     'KCNDS PCNDS   PRAD KSURF PSURF   KFIL  PFIL   KGMP  PGMP',
     *     '    KE',6X,'KE   APE')
  904 FORMAT ( '0  0',I7,I5,I6,I8,4I6,I8,I6,I7,2I6,I7,I6,I7,2I6,I8,I6/)
  905 FORMAT (     I4,I7,I5,I6,I8,4I6,I8,I6,I7,2I6,I7,I6,I7,2I6,I8,I6)
  906 FORMAT (' EDDY',I6,I5,I6,I8,4I6,I8,I6,I7,2I6,I7,I6,I7,2I6,I8,I6)
  907 FORMAT ('0TOTL',I6,I5,I6,I8,4I6,I8,I6,I7,2I6,I7,I6,I7,2I6,I8,I6)
  909 FORMAT (/'0TPE',I18,I32,I14,I7,I12,2I13,I20)
      END SUBROUTINE DIAG5P

      SUBROUTINE DIAGDD
!@sum  DIAGDD prints out diurnal cycle diagnostics
!@auth G. Russell
      use model_com, only: modelEclock
      use TimeConstants_mod, only: HOURS_PER_DAY
      USE MODEL_COM, only :
     &     idacc,JDATE0,AMON,AMON0,JYEAR0,XLABEL,LRUNID,NDAY
      USE DIAG_COM, only :   kdiag,qdiag,units_dd,ndiupt,
     &     adiurn,ijdd,namdd,ndiuvar,hr_in_day,scale_dd,lname_dd,name_dd
     *     ,denom_dd,ia_12hr
      USE MDIAG_COM, only : acc_period
      IMPLICIT NONE

      REAL*8, DIMENSION(HR_IN_DAY+1) :: XHOUR
      INTEGER, DIMENSION(HR_IN_DAY+1) :: MHOUR
      REAL*8 :: AVE,AVED,AVEN,BYIDAC
      INTEGER :: I,IH,IREGF,IREGL,IS,K,KP,KQ,KR,NDAYS,KF,KNDIU,KR1,KR2
      CHARACTER*16, DIMENSION(NDIUVAR) :: UNITSO,LNAMEO,SNAMEO
      REAL*8, DIMENSION(HR_IN_DAY+1,NDIUVAR) :: FHOUR
      CHARACTER :: CPOUT*2
      integer :: year, date

      call modelEclock%get(year=year, date=date)

C****
      NDAYS=IDACC(ia_12hr)/2
      IF (NDAYS.LE.0) RETURN
      BYIDAC=HOURS_PER_DAY/(NDAY*NDAYS)
C****
      IREGF=1
      IREGL=NDIUPT-KDIAG(6)       ! kd6=KDIAG(6)>0: skip last kd6 points
      IF (KDIAG(6).LT.0.AND.KDIAG(6).GE.-NDIUPT) IREGF=-KDIAG(6)
      IF (KDIAG(6).LT.0) IREGL=IREGF       ! kd6<0: show only point -kd6
C**** for netcdf limits, loop in steps of 2000
      KNDIU=0
      DO KQ=1,NDIUVAR
        IF (LNAME_DD(KQ) == "unused") CYCLE
        KNDIU=KNDIU+1
      END DO
      DO KF=1,1+(KNDIU*(IREGL-IREGF+1)-1)/2000
C**** OPEN PLOTTABLE OUTPUT FILE IF DESIRED
      KR1=IREGF+(KF-1)*INT(2000/KNDIU)
      KR2=MIN(IREGL,IREGF+KF*INT(2000/KNDIU)-1)
      IF (QDIAG) THEN
        CPOUT=""
        IF (KNDIU*(IREGL-IREGF+1)/2000 > 1) THEN ! more than one file
          IF (KF <= 9) THEN
            WRITE(CPOUT(1:1),'(I1)') KF
          ELSE
            WRITE(CPOUT(1:2),'(I2)') KF
          END IF
        END IF
        call open_diurn (trim(acc_period)//'.diurn'//trim(cpout)
     *      //XLABEL(1:LRUNID),hr_in_day,KNDIU,KR1,KR2)
      END IF
C**** LOOP OVER EACH BLOCK OF DIAGS
      DO KR=KR1,KR2
        WRITE (6,901) XLABEL(1:105),JDATE0,AMON0,JYEAR0,DATE,AMON,YEAR
        WRITE (6,903) NAMDD(KR),IJDD(1,KR),IJDD(2,KR),(I,I=1,HR_IN_DAY)
C**** KP packs the quantities for postprocessing (skipping unused)
        KP = 0
        DO KQ=1,NDIUVAR
          IF (MOD(KQ-1,5).eq.0) WRITE(6,*)
          IF (LNAME_DD(KQ).eq."unused") CYCLE
          KP = KP+1
          IF(DENOM_DD(KQ).EQ.0) THEN
C**** NORMAL QUANTITIES
            AVE=0.
            DO IH=1,HR_IN_DAY
              AVE=AVE+ADIURN(KQ,KR,IH)
              XHOUR(IH)=ADIURN(KQ,KR,IH)*SCALE_DD(KQ)*BYIDAC
            END DO
            XHOUR(HR_IN_DAY+1)=AVE/FLOAT(HR_IN_DAY)*SCALE_DD(KQ)*BYIDAC
          ELSE
C**** RATIO OF TWO QUANTITIES
            AVEN=0.
            AVED=0.
            DO IH=1,HR_IN_DAY
              AVEN=AVEN+ADIURN(KQ,KR,IH)
              AVED=AVED+ADIURN(DENOM_DD(KQ),KR,IH)
              XHOUR(IH)=ADIURN(KQ,KR,IH)*SCALE_DD(KQ)/
     *             (ADIURN(DENOM_DD(KQ),KR,IH)+1D-20)
            END DO
            XHOUR(HR_IN_DAY+1)=AVEN*SCALE_DD(KQ)/(AVED+1D-20)
          ENDIF
          DO IS=1,HR_IN_DAY+1
            FHOUR(IS,KP)=XHOUR(IS)
            MHOUR(IS)=NINT(XHOUR(IS))
          END DO
          WRITE (6,904) LNAME_DD(KQ),MHOUR
          SNAMEO(KP)=NAME_DD(KQ)(1:16)
          LNAMEO(KP)=LNAME_DD(KQ)(1:16)
          UNITSO(KP)=UNITS_DD(KQ)(1:16)
        END DO
        IF (QDIAG) CALL POUT_DIURN(SNAMEO,LNAMEO,UNITSO,FHOUR,
     *       NAMDD(KR),IJDD(1,KR),IJDD(2,KR),HR_IN_DAY,KP)
      END DO
      IF (QDIAG) call close_diurn
      END DO

      RETURN
C****
  901 FORMAT ('1',A,I3,1X,A3,I5,' - ',I3,1X,A3,I5)
  903 FORMAT ('0',A4,I2,',',I2,' ',I2,23I5,'  AVE')
  904 FORMAT (A8,25I5)
      END SUBROUTINE DIAGDD

      SUBROUTINE DIAGDH
!@sum  DIAGDH prints out hourly diurnal cycle diagnostics
!@+       It uses the same quantities as DIAGHH and shares some arrays
!@+       When radiation is not called every hour this will not average
!@+       exactly to same numbers as in DIAGDD.
!@auth J. Lerner
#ifdef USE_HDIURN
      use TimeConstants_mod, only: HOURS_PER_DAY
      USE MODEL_COM, only :   modelEclock, NDAY, calendar,
     &     idacc,JDATE0,AMON,AMON0,JYEAR0,XLABEL,LRUNID
      USE DIAG_COM, only :   kdiag,qdiag,units_dd,hr_in_month
     *     ,hdiurn,ijdd,namdd,ndiuvar,hr_in_day,scale_dd,lname_dd
     *     ,name_dd,denom_dd,ia_12hr,NDIUPT
      USE MDIAG_COM, only : acc_period
      use CalendarMonth_mod
      IMPLICIT NONE
      REAL*8, DIMENSION(HR_IN_MONTH) :: XHOUR
      INTEGER, DIMENSION(HR_IN_MONTH) :: MHOUR
      INTEGER :: I,IH,IH0,IREGF,IREGL,IS,JD,jdayofm,K,KP,KQ,KR,NDAYS,KF,
     &     KNDIU,KR1,KR2
      CHARACTER*16, DIMENSION(NDIUVAR) :: UNITSO,LNAMEO,SNAMEO
      REAL*8, DIMENSION(HR_IN_MONTH,NDIUVAR) :: FHOUR
      CHARACTER :: CPOUT*2
      integer :: year, month, date
      type (CalendarMonth) :: cMonth

      call modelEclock%get(year=year, month=month, date=date)

C****
      NDAYS=IDACC(ia_12hr)/2
      IF (NDAYS.LE.0) RETURN
C****
C**** KP packs the quantities for postprocessing (skipping unused)

      cMonth = calendar%getCalendarMonth(month, year)
      jdayofM = cMonth%daysInMonth
      IREGF=1
      IREGL=NDIUPT-KDIAG(13)      ! kd13=KDIAG(13)>0: skip last kd13 pts
      IF (KDIAG(13).LT.0.AND.KDIAG(13).GE.-NDIUPT) IREGF=-KDIAG(13)
      IF (KDIAG(13).LT.0) IREGL=IREGF       ! kd13<0: show only pt -kd13
C**** for netcdf limits, loop in steps of 2000
      KNDIU=0
      DO KQ=1,NDIUVAR
        IF (LNAME_DD(KQ) == "unused") CYCLE
        KNDIU=KNDIU+1
      END DO
      DO KF=1,1+(KNDIU*(IREGL-IREGF+1)-1)/2000
      KR1=IREGF+(KF-1)*INT(2000/KNDIU)
      KR2=MIN(IREGL,IREGF+KF*INT(2000/KNDIU)-1)
      IF (QDIAG) THEN
        CPOUT=""
        IF (KNDIU*(IREGL-IREGF+1)/2000 > 1) THEN ! more than one file
          IF (KF <= 9) THEN
            WRITE(CPOUT(1:1),'(I1)') KF
          ELSE
            WRITE(CPOUT(1:2),'(I2)') KF
          END IF
        END IF
C**** OPEN PLOTTABLE OUTPUT FILE IF DESIRED
        call open_hdiurn (trim(acc_period)//'.hdiurn'//trim(cpout)
     &       //XLABEL(1:LRUNID),hr_in_month,KNDIU,KR1,KR2)
      END IF
C**** LOOP OVER EACH BLOCK OF DIAGS
      DO KR=KR1,KR2
        WRITE (6,901) XLABEL(1:105),JDATE0,AMON0,JYEAR0,DATE,AMON,YEAR
        WRITE (6,903)NAMDD(KR),IJDD(1,KR),IJDD(2,KR),(I,I=1,HR_IN_DAY)
C**** KP packs the quantities for postprocessing (skipping unused)
        KP = 0
        DO KQ=1,NDIUVAR
          IF (MOD(KQ-1,5).eq.0) WRITE(6,*)
          IF (LNAME_DD(KQ).eq."unused") CYCLE
          KP = KP+1
          IF(DENOM_DD(KQ).EQ.0) THEN
C**** NORMAL QUANTITIES
            DO IH=1,HR_IN_MONTH
              XHOUR(IH)=HDIURN(KQ,KR,IH)*SCALE_DD(KQ)*
     &                  (HOURS_PER_DAY/NDAY)
            END DO
          ELSE
C**** RATIO OF TWO QUANTITIES
            DO IH=1,HR_IN_MONTH
              XHOUR(IH)=HDIURN(KQ,KR,IH)*SCALE_DD(KQ)/
     *             (HDIURN(DENOM_DD(KQ),KR,IH)+1D-20)
            END DO
          ENDIF
          DO IS=1,HR_IN_MONTH
            FHOUR(IS,KP)=XHOUR(IS)
            MHOUR(IS)=NINT(XHOUR(IS))
          END DO
          ih0 = 1
          do jd = 1,jdayofm
            WRITE (6,904) LNAME_DD(KQ),(MHOUR(i),i=ih0,ih0+23),jd
            ih0 = ih0+24
          end do
          SNAMEO(KP)=NAME_DD(KQ)(1:16)
          LNAMEO(KP)=LNAME_DD(KQ)(1:16)
          UNITSO(KP)=UNITS_DD(KQ)(1:16)
        END DO
        IF (QDIAG) CALL POUT_HDIURN(SNAMEO,LNAMEO,UNITSO,FHOUR,
     *     NAMDD(KR),IJDD(1,KR),IJDD(2,KR),HR_IN_MONTH,KP)
      END DO
      IF (QDIAG) call close_hdiurn
      END DO
#endif
      RETURN
C****
  901 FORMAT ('1',A,I3,1X,A3,I5,' - ',I3,1X,A3,I5)
  903 FORMAT ('0',A4,I2,',',I2,' ',I2,23I5,'  Day')
  904 FORMAT (A8,24I5,I5)
      END SUBROUTINE DIAGDH


      SUBROUTINE DIAG4
!@sum  DIAG4 prints out a time history of the energy diagnostics
!@auth G. Russell
      USE CONSTANT, only :
     &     grav,rgas,bygrav
      use model_com, only: modelEclock
      use TimeConstants_mod, only: HOURS_PER_DAY
      USE MODEL_COM, only :
     &     IDACC,JHOUR0,JDATE0,AMON,AMON0,
     &     JYEAR0,NDAY,Itime0,XLABEL
      USE DIAG_COM, only : im,jm,lm,fim,
     &     ia_inst,ia_d4a,nda4
      USE GC_COM, only : istrat,energy,ned,nehist,hist_days
      USE DIAG_ZONAL, only : xwon
      IMPLICIT NONE

      REAL*8, DIMENSION(2) :: FAC
      REAL*8, DIMENSION(NED) :: SCALET
      REAL*8, DIMENSION(2*NED) :: SUME
      INTEGER, DIMENSION(2*NED) :: IK
      REAL*8, DIMENSION(NEHIST,HIST_DAYS+1) :: EHIST

      INTEGER ::
     &     I,IDACC5,ItimeX,IDAYX,IDAYXM,K,K0,KS,KN,KSPHER
      REAL*8 :: TOFDYX
      integer :: year, hour, date

      call modelEclock%get(year=year, hour=hour, date=date)

      IDACC5=IDACC(ia_d4a)
      IF (IDACC5.LE.0) RETURN
      IF (IDACC(ia_inst).LT.1) IDACC(ia_inst)=1
      SCALET(1)=100.D-18*BYGRAV
      SCALET(2)=SCALET(1)
      SCALET(3)=SCALET(1)
      SCALET(4)=SCALET(1)
c     SCALET(5)=.5*SCALET(1)
      SCALET(5)=SCALET(1)
      SCALET(6)=SCALET(5)
      SCALET(7)=SCALET(1)*RGAS
      SCALET(8)=SCALET(7)
      SCALET(9)=SCALET(7)
      SCALET(10)=SCALET(7)
      DO K=1,NED
        SCALET(K)=XWON*SCALET(K)/IDACC(ia_inst)
      END DO
C****
      DO K0=1,MIN(1+ISTRAT,2)
        WRITE (6,901) XLABEL
        IF (K0.eq.1) THEN
          FAC(1) = 1.
          FAC(2) = 10.  ! a factor of 10 for LOW STRAT
          WRITE (6,902) JYEAR0,AMON0,JDATE0,JHOUR0,YEAR,AMON,DATE
     *         ,HOUR
          WRITE (6,903)
        ELSE
          FAC(1) = 10.  ! 10 goes from 10^18 to 10^17
          FAC(2) = 100. ! another factor of 10 for HIGH STRAT
          WRITE (6,906) JYEAR0,AMON0,JDATE0,JHOUR0,YEAR,AMON,DATE
     *         ,HOUR
          WRITE (6,907)
        END IF
        SUME(:)=0.
        DO I=1,IDACC5
          ItimeX=Itime0+I*NDA4-1
          IDAYX=1+ItimeX/NDAY
          IDAYXM=MOD(IDAYX,100000)
          TOFDYX=MOD(ItimeX,NDAY)*HOURS_PER_DAY/NDAY
          DO KSPHER=1,2
            DO K=1,NED
              KS=K+(KSPHER-1)*NED
              KN=KS+(K0-1)*2*NED
              IF (KN.le.NEHIST) THEN
                EHIST(KN,I)=ENERGY(KN,I)*SCALET(K)*FAC(KSPHER)
                IK(KS)=EHIST(KN,I)+.5
                SUME(KS)=SUME(KS)+ENERGY(KN,I)
              ELSE
                IK(KS)=-999
              END IF
            END DO
          END DO
          WRITE (6,904) IDAYXM,TOFDYX,IK
        END DO
        DO KSPHER=1,2
          DO K=1,NED
            KS=K+(KSPHER-1)*NED
            KN=KS+(K0-1)*2*NED
            IF (KN.le.NEHIST) THEN
              EHIST(KN,HIST_DAYS+1)=SUME(KS)*SCALET(K)*FAC(KSPHER)
     *             /IDACC5
              IK(KS)=EHIST(KN,HIST_DAYS+1)+.5
            ELSE
              IK(KS)=-999
            END IF
          END DO
        END DO
        WRITE (6,905) IK
        IF (K0.eq.1) CALL KEYD4 (IK)
      END DO
      RETURN
C****
  901 FORMAT ('1',A)
  902 FORMAT ('0** ENERGY HISTORY **      From:',I6,A6,I2,',  Hr',I3,
     *  6X,'To:',I6,A6,I2,', Hr',I3,
     *  '    UNITS OF 10**18 JOULES')
  903 FORMAT ('0',15X,21('-'),' TROPOSPHERE ',22('-'),5X,21('-'),
     *  '  LOW STRAT. * 10 ',17('-')/8X,2(11X,'ZKE',8X,'EKE',7X,
     *     'SEKE',9X,
     * 'ZPE',10X,'EPE')/3X,'DAY  HOUR     SH   NH    SH   NH     1    2
     *    SH    NH     SH    NH      SH   NH    SH   NH    SH   NH     S
     *H    NH     SH    NH'/1X,132('='))
  904 FORMAT (I6,F6.1,1X,3(I6,I5),2(I7,I6),2X,3(I6,I5),2(I7,I6))
  905 FORMAT (1X,132('=')/8X,'MEAN ',3(I6,I5),2(I7,I6),2X,3(I6,I5),
     *  2(I7,I6))
  906 FORMAT ('0** ENERGY HISTORY **      From:',I6,A6,I2,',  Hr',I3,
     *  6X,'To:',I6,A6,I2,', Hr',I3,
     *  '    UNITS OF 10**17 JOULES')
  907 FORMAT ('0',15X,19('-'),' MID STRATOSPHERE ',19('-'),5X,18('-'),
     *  ' HIGH STRAT. * 10  ',19('-')/8X,2(11X,'ZKE',8X,'EKE',7X,
     *   'NHKE',9X,
     *  'ZPE',10X,'EPE')/3X,'DAY  HOUR     SH   NH    SH   NH    1    2
     *    SH    NH     SH    NH      SH   NH    SH   NH    1    2      S
     *H    NH     SH    NH'/1X,132('='))
  920 FORMAT (1X)
      END SUBROUTINE DIAG4


      subroutine KEYVSUMS (QUANT,GSUM,HSUM,ASUM,SUMFAC)
      use diag_com, only : jm
      use mdiag_com, only : sname_strlen
      implicit none
!@var quant string designating the quantity for which to save keynrs
      CHARACTER(LEN=sname_strlen) :: QUANT
      REAL*8, DIMENSION(JM) :: ASUM
      REAL*8, DIMENSION(2) :: HSUM
      REAL*8 :: GSUM,SUMFAC
      if(quant.eq.'tx') CALL KEYJKT (GSUM,ASUM)
      if(quant.eq.'eddy_ke') CALL KEYJKE (18,HSUM,ASUM)
      if(quant.eq.'tot_ke') CALL KEYJKE (19,HSUM,ASUM)
      if(quant.eq.'nt_dse_stand_eddy') CALL KEYJKN (33,ASUM,SUMFAC)
      if(quant.eq.'nt_dse_eddy') CALL KEYJKN (34,ASUM,SUMFAC)
      if(quant.eq.'tot_nt_dse') CALL KEYJKN (35,ASUM,SUMFAC)
!!!   if(quant.eq.'nt_lh_e') CALL KEYJKN (??,ASUM,SUMFAC)
!!!   if(quant.eq.'tot_nt_lh') CALL KEYJKN (??,ASUM,SUMFAC)
      if(quant.eq.'nt_se_eddy') CALL KEYJKN (36,ASUM,SUMFAC)
      if(quant.eq.'tot_nt_se') CALL KEYJKN (37,ASUM,SUMFAC)
!!!   if(quant.eq.'tot_nt_ke') CALL KEYJKN (??,ASUM,SUMFAC)
      if(quant.eq.'nt_u_stand_eddy') CALL KEYJKN (39,ASUM,SUMFAC)
      if(quant.eq.'nt_u_eddy') CALL KEYJKN (40,ASUM,SUMFAC)
      if(quant.eq.'tot_nt_u') CALL KEYJKN (41,ASUM,SUMFAC)
      RETURN
      end subroutine keyvsums


      subroutine keynrl(quant,l,flat)
      use diag_com, only : jm
      use mdiag_com, only : sname_strlen
      implicit none
      integer :: l
      REAL*8, DIMENSION(JM) :: FLAT
!@var quant string designating the quantity for which to save keynrs
      CHARACTER(LEN=sname_strlen) :: QUANT
      if(quant.eq.'u') CALL KEYJKJ (L,FLAT)
      if(quant.eq.'psi_cp') CALL KEYJLS (L,FLAT)
      return
      end subroutine keynrl

      SUBROUTINE IJKMAP (iu_Iij)
!@sum  IJKMAP output 3-D constant pressure output fields
!@auth G. Schmidt
C**** Note that since many IJK diags are weighted w.r.t pressure, all
C**** diagnostics must be divided by the accumulated pressure
C**** All titles/names etc. implicitly assume that this will be done.
      USE CONSTANT, only : grav,sha,undef
      USE MODEL_COM, only : XLABEL,LRUNID,idacc
      USE ATM_COM, only : pmidl00
      USE DIAG_COM, only : im,jm,lm
     &     ,kdiag,jgrid_ijk,aijk,denom_ijk
     *     ,scale_ijk,off_ijk,name_ijk,lname_ijk,units_ijk,kaijk
     *     ,ia_dga
      use mdiag_com, only : acc_period
      use gcdiag, only : ijk_dpb
      use filemanager
      IMPLICIT NONE

      CHARACTER XLB*24,TITLEX*56
      CHARACTER*80 TITLEL(LM)
      REAL*8 SMAP(IM,JM,LM),SMAPJK(JM,LM),SMAPK(LM)
      REAL*8 flat,dp
      CHARACTER*8 CPRESS(LM)
      INTEGER i,j,l,kxlb,ni,k,iu_Iij
      logical, dimension (kaijk) :: Qk

C****
C**** INITIALIZE CERTAIN QUANTITIES
C****

      Qk = .true.
      do k=1,kaijk
        if (lname_ijk(k).eq.'unused') Qk(k) = .false.
      end do
      if (kdiag(3).eq.9) then
         write (iu_Iij,'(a)') 'list of 3-d fields'
         do k=1,kaijk
           if (lname_ijk(k).ne.'unused')
     *        write (iu_Iij,'(i3,1x,a)') k,lname_ijk(k)
         end do
         call closeunit(iu_Iij)
         return
      else if (kdiag(3).gt.1) then
         Qk = .false.
   10    read (iu_Iij,'(i3)',end=20) k
         Qk(k) = .true.
         go to 10
   20    continue
         call closeunit(iu_Iij)
      end if

C**** OPEN PLOTTABLE OUTPUT FILE
      call open_ijk(trim(acc_period)//'.ijk'//XLABEL(1:LRUNID),im,jm,lm)
      KXLB = INDEX(XLABEL(1:11),'(')-1
      IF(KXLB.le.0) KXLB = 10
      XLB = ' '
      XLB(1:13)=acc_period(1:3)//' '//acc_period(4:12)
      XLB(15:14+KXLB) = XLABEL(1:KXLB)
C****
C**** Complete 3D-field titles
C****
      DO L=1,LM
        WRITE(CPRESS(L),'(F8.3)') pmidl00(l)
      END DO

C**** Select fields
      DO K=1,Kaijk
        if (.not.Qk(k).or.k.eq.ijk_dpb) cycle
        SMAP(:,:,:) = UNDEF
        SMAPJK(:,:) = UNDEF
        SMAPK(:)    = UNDEF
        IF(denom_ijk(k).eq.0) THEN ! no weighting, assuming qty on a-grid
          TITLEX = lname_ijk(k)(1:17)//"   at        mb ("//
     *         trim(units_ijk(k))//")"
          DO L=1,LM
            DO J=1,JM
              NI = 0
              FLAT = 0.
              DO I=1,IM
                SMAP(I,J,L) = SCALE_IJK(k)*AIJK(I,J,L,K)/IDACC(ia_dga)
                FLAT = FLAT+SMAP(I,J,L)
                NI = NI+1
              END DO
              IF (NI.GT.0) SMAPJK(J,L) = FLAT/NI
            END DO
            WRITE(TITLEX(23:30),'(A)') CPRESS(L)
            TITLEL(L) = TITLEX//XLB
          END DO
        ELSE    ! for now, assuming weight is b-grid dp
          TITLEX = lname_ijk(k)(1:17)//"   at        mb ("//
     *         trim(units_ijk(k))//", UV grid)"
          DO L=1,LM
            DO J=2,JM
              NI = 0
              FLAT = 0.
              DO I=1,IM
                DP=AIJK(I,J,L,IJK_DPB)
                IF(DP.GT.0.) THEN
                  SMAP(I,J,L)=SCALE_IJK(K)*AIJK(I,J,L,K)/DP+OFF_IJK(K)
                  FLAT = FLAT+SMAP(I,J,L)
                  NI = NI+1
                END IF
              END DO
              IF (NI.GT.0) SMAPJK(J,L) = FLAT/NI
            END DO
            WRITE(TITLEX(23:30),'(A)') CPRESS(L)
            TITLEL(L) = TITLEX//XLB
          END DO
        END IF
        CALL POUT_IJK(TITLEL,name_ijk(k),lname_ijk(k),units_ijk(k)
     *       ,SMAP,SMAPJK,SMAPK,jgrid_ijk(k))
      END DO
C****
      call close_ijk

      return
      END SUBROUTINE IJKMAP

      SUBROUTINE IJLMAP (iu_Iij)
C****
C**** ijl output
C****
      USE CONSTANT, only : undef
      USE MODEL_COM, only : XLABEL,LRUNID,idacc
      USE DIAG_COM, only : im,jm,lm
     &     ,aijl,kaijl,ijkgridc
     &     ,ia_ijl,denom_ijl,name_ijl,lname_ijl,units_ijl,scale_ijl,
     &     lgrid_ijl,jgrid_ijl,
     &     ctr_ml,edg_ml,ctr_cp
      USE MDIAG_COM, only : acc_period
      USE ATM_COM, only : pmidl00
      use filemanager
      IMPLICIT NONE
      INTEGER iu_Iij

      CHARACTER XLB*24,TITLEX*56
      CHARACTER*80 TITLEL(LM)
      REAL*8 SMAP(IM,JM,LM),SMAPJL(JM,LM),SMAPL(LM)
      REAL*8 XDEN(IM,JM,LM),XDENJL(JM,LM)
      CHARACTER*11 CPRESS(LM)
      CHARACTER*3 CLEV(LM)
      INTEGER i,j,l,kxlb,k,kd,id
      KXLB = INDEX(XLABEL(1:11),'(')-1
      IF(KXLB.le.0) KXLB = 10
      XLB = ' '
      XLB(1:13)=acc_period(1:3)//' '//acc_period(4:12)
      XLB(15:14+KXLB) = XLABEL(1:KXLB)
      DO L=1,LM
        WRITE(CPRESS(L),'(F8.3,A3)') pmidl00(l),' mb'
        WRITE(CLEV(L),'(A1,I2.2)') 'L',l
      END DO

      call open_ijl(trim(acc_period)//'.ijl'//XLABEL(1:LRUNID),im,jm,lm,
     &     kaijl,name_ijl,lname_ijl,units_ijl,lgrid_ijl)

c
c loop over quantities
c
C     Fill in the undefined pole box duplicates
      do k=1,kaijl
        if (jgrid_ijl(k) == 2) CYCLE
        do L=1,lm
           aijl(2:im, 1,L,k) = aijl(1, 1,L,k)
           aijl(2:im,jm,L,k) = aijl(1,jm,L,k)
        end do
      end do

      do k=1,kaijl
        if(trim(lname_ijl(k)).eq.'no output') cycle
c
c scale this quantity and compute its zonal means
c
        kd=denom_ijl(k)
        id=2
        if (kd.gt.0) id=3
        call scale_ijlmap(id,aijl(:,:,:,k),aijl(:,:,:,kd),scale_ijl(k)
     *       ,kd,idacc(ia_ijl(k)),idacc(ia_ijl(kd)),smap,smapjl)
C****
C**** Complete 3D-field titles
C****
        titlex = trim(lname_ijl(k))//' ('//trim(units_ijl(k))//')'
        if(lgrid_ijl(k).eq.ctr_ml.or.lgrid_ijl(k).eq.edg_ml) then
          do l=1,lm
            titlel(l)=trim(titlex)//' '//clev(l)
            titlel(l)(57:80)=xlb
          enddo
        elseif(lgrid_ijl(k).eq.ctr_cp) then
          do l=1,lm
            titlel(l)=trim(titlex)//' '//cpress(l)
            titlel(l)(57:80)=xlb
          enddo
        endif
c write field
        call pout_ijl(titlel,name_ijl(k),lname_ijl(k),units_ijl(k)
     &       ,smap,smapjl,smapl,ijkgridc)
      enddo

      call close_ijl

      return
      end subroutine ijlmap

      subroutine scale_ijlmap(nmap,aijl1,aijl2,scale,kd,id1,id2,smap,
     *     smapjl)
!@sum scale and calculate zonal means for 3D diag
      use constant, only : undef
      use diag_com, only : im,jm,lm
      use geom, only : dxyp
      implicit none
      integer, intent(in) :: nmap  ! type of diag
      real*8, intent(in) :: aijl1(im,jm,lm),aijl2(im,jm,lm)
      real*8, intent(in) :: scale
      integer, intent(in) :: id1,id2,kd
      real*8, intent(out) :: smap(im,jm,lm)
      real*8, intent(out) :: smapjl(jm,lm)
      real*8 :: xden(im,jm,lm),xdenjl(jm,lm)
      integer j,l

      if (nmap.eq.1) then  ! area weighting
        do l=1,lm
          do j=1,jm
            smap(:,j,l) = scale*aijl1(:,j,l)/(id1*dxyp(j))
          end do
        end do
        smapjl = sum(smap,dim=1)/im
      elseif (nmap.eq.2) then ! simple cases - no division by area
        smap = scale*aijl1(:,:,:)/id1
        smapjl = sum(smap,dim=1)/im
      elseif (nmap.eq.3) then ! ratio
        smap = scale*aijl1(:,:,:)/id1
        smapjl = sum(smap,dim=1)
        xden = aijl2(:,:,:)/id2
        xdenjl = sum(xden,dim=1)
        where(xden.ne.0.)
          smap = smap/xden
        elsewhere
          smap = undef
        end where
        where(xdenjl.ne.0.)
          smapjl = smapjl/xdenjl
        elsewhere
          smapjl = undef
        end where
      else
        write(6,*) "Incorrect nmap type in scale_ijlmap",nmap
        call stop_model('scale_ijlmap: undefined nmap type',255)
      end if

      return
      end subroutine scale_ijlmap

      function NINTlimit( x )
      real*8 x
      integer NINTlimit
      real*8 y
      y = min (  2147483647.d0, x )
      y = max ( -2147483647.d0, y )
      NINTlimit = NINT( y )
      return
      end function NINTlimit

      subroutine diag_isccp
!@sum diag_isccp prints out binary and prt output for isccp histograms
!@auth Gavin Schmidt
      USE MODEL_COM, only : xlabel,lrunid,idacc
      USE GEOM, only : dxyp,lat_dg
      USE DIAG_COM, only : im,jm,fim,aisccp,ntau,npres,nisccp
     *     ,qdiag,ia_src,isccp_press,isccp_taum,aij,ij_tcldi,ij_scldi
     &     ,isccp_late,wisccp
      USE MDIAG_COM, only : acc_period,
     &     sname_strlen,units_strlen,lname_strlen
      IMPLICIT NONE

      CHARACTER*80 :: TITLE(nisccp) = (/
     *     "ISCCP CLOUD FREQUENCY (NTAU,NPRES) % 60S-30S",
     *     "ISCCP CLOUD FREQUENCY (NTAU,NPRES) % 30S-15S",
     *     "ISCCP CLOUD FREQUENCY (NTAU,NPRES) % 15S-15N",
     *     "ISCCP CLOUD FREQUENCY (NTAU,NPRES) % 15N-30N",
     *     "ISCCP CLOUD FREQUENCY (NTAU,NPRES) % 30N-60N" /)
      REAL*8 AX(ntau-1,npres,nisccp)
      INTEGER N,ITAU,IPRESS,J,I

      character(len=sname_strlen) :: sname
      character(len=units_strlen) :: units
      character(len=lname_strlen) :: lname

C**** write out scaled results
      do n=nisccp,1,-1   ! north to south
        title(n)(61:72) = acc_period
        write(6,100) title(n)
        AX(1:ntau-1,:,n) = 100.*AISCCP(2:ntau,:,n)/wisccp(n)
        do ipress=1,npres
          write(6,101) isccp_press(ipress),(AX(itau,ipress,n),itau=1
     *         ,ntau-1)
        end do
        write(6,*)
      end do
C**** write the binary file
      if (qdiag) then
        call open_isccp(trim(acc_period)//'.isccp'//
     *       XLABEL(1:LRUNID),ntau-1,npres,nisccp)
         sname='pcld'
         lname='cloud cover histogram'
         units='%'
         call pout_isccp(title,sname,lname,units,ax,isccp_taum
     *        ,real(isccp_press,kind=8))
         call close_isccp
      endif
      RETURN

 100  FORMAT (1X,A80/1X,72('-')/3X,
     *     'PRESS\TAU    0.  1.3  3.6  9.4  23   60   > ')
 101  FORMAT (5X,I3,7X,6F5.1)

      end subroutine diag_isccp

      subroutine keyij2(sst)
!@sum output basic ENSO diagnostics to PRT file
!@auth Gavin Schmidt
      USE GEOM, only : adxyp,lat_dg,lon_dg
      USE DIAG_COM, only : keynr,keyct
      IMPLICIT NONE
      REAL*8, DIMENSION(IM,JM) :: SST
      REAL*8 nino34,wt
      INTEGER i,j

      nino34=0.
      wt=0.
      DO J=1,JM
        DO I=1,IM
          if (lat_dg(j,1).ge.-6.   .and. lat_dg(j,1).le.6 .and.
     *        lon_dg(i,1).ge.-170. .and. lon_dg(i,1).le.-120.) then
            nino34=nino34+adxyp(i,j)*sst(i,j)
            wt=wt+adxyp(i,j)
          end if
        END DO
      END DO
      print*,"nino",nino34/wt,nino34,wt
      keynr(43,KEYCT)=NINT(100.0*nino34/wt)

      end subroutine keyij2

      SUBROUTINE DIAG_GATHER
      USE ATM_COM, only : ZATMO
      USE FLUXES, only : FOCEAN, FLICE
cddd      USE LAKES_COM, only : FLAKE
cddd      USE GHY_COM, only : FEARTH
cddd      use ent_com, only : entcells
cddd      use ent_mod, only : ent_get_exports
      USE DIAG_COM, only : IM, AIJ,  AIJ_loc, AJ,   AJ_loc,
     *     AIJK, AIJK_loc,
     *     ASJL, ASJL_loc, AJL,  AJL_loc , CONSRV, CONSRV_loc, TSFREZ,
     *     TSFREZ_loc, WT_IJ
      USE GC_COM, only : AGC, AGC_loc
      USE DOMAIN_DECOMP_ATM, ONLY : GRID
      USE DOMAIN_DECOMP_1D, ONLY : am_i_root,pack_data
      USE CONSTANT, only : NaN
      IMPLICIT NONE
cddd      INTEGER :: J_0, J_1, J_0H, J_1H
cddd      REAL*8, ALLOCATABLE :: tmp(:,:)
cddd      REAL*8, ALLOCATABLE :: fract_vege(:,:)
cddd      INTEGER i,j

      call alloc_ijdiag_glob

      call Gather_Diagnostics()

#ifdef TRACERS_ON
      call gather_trdiag
#endif

! Now the external arrays
      !!CALL PACK_DATA(GRID, fland, fland_glob)
      !!CALL PACK_DATA(GRID, fearth, fearth_glob)
      CALL PACK_DATA(GRID, focean, focean_glob)
      CALL PACK_DATA(GRID, flice, flice_glob)
      CALL PACK_DATA(GRID, zatmo, zatmo_glob)

cddd      CALL GET(GRID, J_STRT=J_0, J_STOP=J_1,
cddd     &     J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)
cddd      ALLOCATE(tmp(IM, J_0H:J_1H))
cddd
cddd      wt_ij(:,:,1) = 1.
      wt_ij(:,:,:) = 1.  ! NaN
      !!CALL PACK_DATA(GRID, focean, wt_ij(:,:,2))
      !!CALL PACK_DATA(GRID, flake,  wt_ij(:,:,3))  ! not correct
      !!CALL PACK_DATA(GRID, flice,  wt_ij(:,:,4))
      !!CALL PACK_DATA(GRID, fearth, wt_ij(:,:,5))  ! not correct

cddd      ALLOCATE(fract_vege(IM, J_0H:J_1H))
cddd      call ent_get_exports( entcells(1:IM,J_0:J_1),
cddd     &           fraction_of_vegetated_soil=fract_vege(1:IM,J_0:J_1) )
cddd      tmp(:,J_0:J_1) = fearth(:,J_0:J_1) * (1.d0-fract_vege(:,J_0:J_1))
cddd      CALL PACK_DATA(GRID, tmp, wt_ij(:,:,6))
cddd      tmp(:,J_0:J_1) = fearth(:,J_0:J_1) * fract_vege(:,J_0:J_1)
cddd      CALL PACK_DATA(GRID, tmp, wt_ij(:,:,7))
cddd      DEALLOCATE(fract_vege)
cddd      DEALLOCATE(tmp)

      call gather_odiags () ; call gather_icdiags ()

      END SUBROUTINE DIAG_GATHER

      SUBROUTINE DIAG_SCATTER
      USE DIAG_COM, only : AIJ, AIJ_loc, AJ,  AJ_loc,
     *     AIJK, AIJK_loc, ASJL, ASJL_loc,
     *     AJL,  AJL_loc, TSFREZ, TSFREZ_loc
      USE GC_COM, only : AGC, AGC_loc
      USE DOMAIN_DECOMP_ATM, ONLY : GRID
      USE DOMAIN_DECOMP_1D, ONLY : UNPACK_DATA, UNPACK_DATAj
      USE DOMAIN_DECOMP_1D, ONLY : am_i_root
      IMPLICIT NONE

#ifndef CUBED_SPHERE
      CALL UNPACK_DATAj(GRID, AJ,  AJ_loc)
      CALL UNPACK_DATAj(GRID, ASJL, ASJL_loc)
      CALL UNPACK_DATAj(GRID, AJL,  AJL_loc)
#endif
      CALL UNPACK_DATAj(GRID, AGC, AGC_loc)
      CALL UNPACK_DATA (GRID, AIJ, AIJ_loc)
      CALL UNPACK_DATA (GRID, AIJK, AIJK_loc)
      CALL UNPACK_DATA (GRID, TSFREZ,  TSFREZ_loc)

#ifdef TRACERS_ON
      call scatter_trdiag
#endif
      call dealloc_ijdiag_glob

      END SUBROUTINE DIAG_SCATTER

      END MODULE DIAG_SERIAL

      subroutine print_diags(partial)
!@sum print_diag prints out binary and ascii diag output.
!@auth  Original Development Team
      USE DOMAIN_DECOMP_1D, only : AM_I_ROOT
      USE DIAG_SERIAL
      USE MODEL_COM, only : itime,itimeI
      USE DIAG_COM, only : kdiag,keynr,keyct,isccp_diags
      IMPLICIT NONE
!@var partial : accum period is complete (if =0) or partial (if =1)
      INTEGER, INTENT(IN) :: partial
 
#ifdef SCM
      return
#endif

      call calc_derived_aij
      call calc_derived_aijk
      if(isccp_diags.eq.1) call diag_isccp_prep
      IF (KDIAG(12).LT.9) CALL diag_OCEAN_prep

      CALL DIAG_GATHER

      IF (AM_I_ROOT()) THEN

      IF (KDIAG(1).LT.9) CALL DIAGJ_PREP
      IF (KDIAG(1).LT.9) CALL DIAGJ
      IF (KDIAG(2).LT.9) CALL DIAGJL_PREP
      IF (KDIAG(2).LT.9) CALL DIAGGC_PREP
      IF (KDIAG(2).LT.9) CALL DIAGJK
      IF (KDIAG(10).LT.9) CALL DIAGIL
      IF (KDIAG(7).LT.9) CALL DIAG7P
      IF (KDIAG(3).LT.9) CALL DIAGIJ
      IF (KDIAG(9).LT.9) CALL DIAGCP
      IF (KDIAG(5).LT.9) CALL DIAG5P
      IF (partial.eq.0 .and. KDIAG(6).LT.9) CALL DIAGDD  ! full period
      IF (KDIAG(13).LT.9) CALL DIAGDH
      IF (KDIAG(4).LT.9) CALL DIAG4

      END IF

      IF (KDIAG(11).LT.9) CALL diag_RIVER

      IF (AM_I_ROOT()) THEN

      IF (KDIAG(12).LT.9) CALL diag_OCEAN
      IF (KDIAG(12).LT.9) CALL diag_ICEDYN
      IF (isccp_diags.eq.1) CALL diag_ISCCP
      IF (partial.eq.0 .or. Itime.LE.ItimeI+1) THEN  ! full period or IC
        CALL DIAGKN
      ELSE                      ! RESET THE UNUSED KEYNUMBERS TO ZERO
        KEYNR(1:42,KEYCT)=0
      END IF
#ifdef TRACERS_ON
      IF (KDIAG(8).LT.9) then
        CALL DIAGJLT
        CALL DIAGIJT
        CALL DIAGIJLT
      end if
#endif
#if (defined TRACERS_ON) || (defined TRACERS_OCEAN)
      IF (KDIAG(8).LT.9) CALL DIAGTCP
#endif
      END IF ! AM_I_ROOT

      CALL DIAG_SCATTER

      return
      end subroutine print_diags
