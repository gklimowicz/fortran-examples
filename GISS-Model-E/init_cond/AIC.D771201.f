C**** AIC.RES_M20A.D771201   DUMMY AIC.D771201
C****
C**** AIC.D771201.S       Create NCAR Initial Conditions    1999/03/26
C****
C**** Compile: fco0 AIC72X46.CRE
C**** Link:   f77 -64 -mips4 AIC72X46.o /u/exec/HNTRPS.o -o AIC72X46.exe
C****
C**** AIC72X46.S reads an NCAR file of observed atmospheric
C**** variables for the specified time, interpolates the data
C**** horizontally and vertically, and produces a datafile of
C**** the prognostic variables necessary to start the GCM.
C****
C**** A(1,J,K) is centered on the Greenwich meridion.
C****
C**** NCAR input file           I.C. file created         Ocean frac
C**** ---------------           -----------------         ----------
C**** NCARIC.R2H2H.D7712010     AIC.C72X46L9.D771201      Z72X46N
C**** NCARIC.R2H2H.D7712010     AIC.C72X46L12.D771201     Z72X46N
C****
      use resolution, only: im,jm,lm,ls1,psf,plbot
      CHARACTER*8  FILEO
      CHARACTER*80 FILEIN,FILEZ,TITLEI, TITLEO,TITLEH(5),TITLEV(6)
      CHARACTER*1 :: GRID='B'
C**** DLAT_DG=width of a non-polar latitude zone in degrees
      PARAMETER (ltm=ls1-1)
      REAL*4 POUT(lm)
      REAL*4 A(144,73,57),AH(IM,JM,57),AHV(IM,JM,4*lm+3),
     *     WTA(144,73), SIGE(0:LTM),PE(0:lm),
     *     RHJ(JM,7:19),UJ(JM,13:19),TJ(JM,13:19)
      INTEGER :: JDOFM(12) =
     *      (/0,31,59,90,120,151,181,212,243,273,304,334/)
      EQUIVALENCE (A,AHV)
      COMMON /FIXDCB/ FOCEAN(IM,JM), ZATMO(IM,JM)
C**** Observed input data : vertical grid
      REAL :: PLEV(19)=(/1000.,850.,700.,500.,400.,300.,
     *            250.,200.,150.,100., 70., 50.,
     *             30., 10.,  3.,  1.,  .3, .05,  0./)
C**** Titles and input/output files
      DATA TITLEH/'U wind after horizontal interpolation',
     *            'V wind after horizontal interpolation',
     *            'Temperature after horizontal interpolation',
     *            'Geopotential Height after horizontal interpolation',
     *            'Relative Humidity after horizontal interpolation'/
      DATA TITLEV/'U wind after vertical interpolation',
     *            'V wind after vertical interpolation',
     *            'Temperature after vertical interpolation',
     *            'Specific Humidity after vertical interpolation',
     *            'Surface Temperature after vertical interpolation',
     *            'Sea Level Pressure after vertical interpolation'/
C****
C**** output pressure levels
      do L=1,lm+1
        pe(L-1)=PLbot(L)
      end do

C**** model zone widths in degrees (defaults)
      DLAT_DG=180./JM                     ! even spacing (default)
      IF (JM.eq.46) DLAT_DG=180./(JM-1)   ! 1/2 box at pole for5
cc    IF (JM.eq.24) DLAT_DG=180./(JM-1)   ! 1/2 box at pole,orig 8x10
      IF (JM.eq.24) DLAT_DG=180./(JM-1.5) ! 1/4 box at pole,real 8x10

C**** Override grid/lat.spacing using command line argument (e.g. C3.00)
      titlei = ' '
      if(iargc() > 0) then
        call getarg(1,titlei(1:50))
        if(iargc() > 1 )  call getarg(2,titlei(52:80))
        i1=1
        if(index('bBcC',titlei(1:1)) > 0) then
          i1=2
          if(index('bBcC',titlei(1:1)) > 2) GRID='C'
        end if
        if(titlei(i1:80).ne.' ') read(titlei(i1:80),*) dlat_dg
      end if

      WTA = 1.
      KDEBUG = 1
      FILEZ  = 'TOPO'           !  topography file
      FILEIN = 'NCARIC.144x73.D7712010'     !  from NCARtape
      FILEO  = 'OUT'//GRID  !  new atm. I.C.
      write(FILEO(5:8),'(f4.2)') dlat_dg
      if(FILEO(5:5)==' ') FILEO(5:5)='0'

C**** Date for observed data
      NYR    = 1977   ! year
      NMM    = 12     ! month
      NDD    = 1      ! day
      NZZ    = 0      ! time of day (hrs)
      JDAY   = JDOFM(NMM) + NDD   ! Julian date

      write (TITLEO(56:80),'(a4,i4,a1,i2.2,a1,i2.2,a7,i2,a2)')
     *       'NMC ',NYR,'/',nmm,'/',ndd,', hour ',nzz,'  '

C**** Calculate model pressure levels POUT
      PSURF  = PE(0) !  global mean surface pressure
      SIGE(0) = 1.
      DO 10 L=1,LTM
      SIGE(L) = (PE(L)-PE(LTM)) / (PSURF-PE(LTM))
   10 POUT(L) =      (PE(LTM) +
     *                 (PSURF-PE(LTM))*.5*(SIGE(L-1)+SIGE(L)))
      DO 20 L=LTM+1,lm
   20 POUT(L) =     (.5*(PE(L-1)+PE(L)))

C**** Read in ocean fraction FOCEAN and atmospheric topography ZATMO
      OPEN  (11,FILE=FILEZ,FORM='UNFORMATTED',STATUS='OLD')
      READ  (11) TITLEI,FOCEAN
      WRITE (6,*) 'FOCEAN read from unit 11: ',TITLEI
      READ  (11)
      READ  (11)
      READ  (11)
      READ  (11) TITLEI,ZATMO
      WRITE (6,*) 'ZATMO read from unit 11: ',TITLEI
      CLOSE (11)

C**** Read in NCAR atmospheric variables at 2.5 x 2.5 resolution
C****
      OPEN (1,FILE=FILEIN,FORM='UNFORMATTED',STATUS='OLD')
      DO 110 K=1,57
      READ (1) TITLEI, (A(I,1,K),I=1,144*73)
  110 WRITE (6,*) 'Read from unit 1: ',TITLEI
      CLOSE (1)
C****
C**** Perform horizontal interpolation from NCAR 2.5 x 2.5 resolution
C****
      IF(GRID.eq.'C')  then
C**** Interpolate winds for C grid
        CALL HNTRP0 (144,73,71.5,72., IM,JM,.5,180./dlat_dg, 0.)
        DO K=1,12
          CALL HNTRPP (WTA,A(1,1,K),AH(1,1,K))
        end do
        CALL HNTRP0 (144,73,71.5,72., IM,JM-1,0.,180./dlat_dg, 0.)
        DO K=13,24
          CALL HNTRP  (WTA,A(1,1,K),AH(1,1,K))
        end do
      else
C**** Interpolate winds for B grid
        CALL HNTRP0 (144,73,71.5,72., IM,JM-1,.5,180./dlat_dg, 0.)
        DO K=1,24
          CALL HNTRP  (WTA,A(1,1,K),AH(1,1,K))
        end do
      end if
C**** Interpolate quantities defined on the A grid
      CALL HNTRP0 (144,73,71.5,72., IM,JM,0.,180./dlat_dg, 0.)
      DO 270 K=25,57
  270 CALL HNTRPP (WTA,A(1,1,K),AH(1,1,K))
C****
C**** Fill in stratospheric zonal wind and temperature and
C**** relative humidity from climatology, set v-winds to 0
C****
      CALL STRDAT (RHJ,TJ,UJ,JM,JDAY,GRID,DLAT_DG,PLEV(7))
C****
C**** Produce line printer maps of horizontally interpolated data
C****
C     CALL MAP  (IM,JM,NINT(PLEV( 4)),TITLEH(1),AH(1,1, 4))
C     CALL MAP  (IM,JM,NINT(PLEV(12)),TITLEH(1),AH(1,1,12))
C     CALL MAP  (IM,JM,NINT(PLEV( 4)),TITLEH(2),AH(1,1,16))
Cx    CALL MAP1 (IM,JM,NINT(PLEV( 4)),TITLEH(3),AH(1,1,28),
Cx   *           WTA,1.,273.15,26)
Cx    CALL MAP1 (IM,JM,NINT(PLEV(13)),TITLEH(3),AH(1,1,41),
Cx   *           WTA,1.,273.15,26)
Cx    CALL MAP1 (IM,JM,NINT(PLEV( 4)),TITLEH(4),AH(1,1,46),
Cx   *           WTA,.1,500.,26)
Cx    CALL MAP  (IM,JM,NINT(PLEV( 4)),TITLEH(5),AH(1,1,60))
C****
C**** Perform vertical interpolation
C****
      CALL VERT (IM,JM,GRID,DLAT_DG, PLEV,AH,RHJ,TJ,UJ,ZATMO, JDAY,
     *           lm,LTM,SIGE,AHV, PE)
C****
C**** Write Initial Conditions to output datafile
C****
      OPEN (2,FILE=FILEO,FORM='UNFORMATTED')
C**** Write out surface pressure
      TITLEO(1:55) = 'SURFACE PRESSURE (mb)'
      WRITE (2) TITLEO, (AHV(I,1,1),I=1,IM*JM)
C**** Write out U velocity component
      DO 420 L=1,lm
      WRITE (TITLEO(1:55),942) 'U VELOCITY COMPONENT (m/s) at',POUT(L)
  420 WRITE (2) TITLEO, (AHV(I,1,1+L),I=1,IM*JM)
C**** Write out V velocity component
      DO 430 L=1,lm
      WRITE (TITLEO(1:55),942) 'V VELOCITY COMPONENT (m/s) at',POUT(L)
  430 WRITE (2) TITLEO, (AHV(I,1,1+lm+L),I=1,IM*JM)
C**** Write out temperature
      DO 440 L=1,lm
      WRITE (TITLEO(1:55),942) 'TEMPERATURE (K) at',POUT(L)
  440 WRITE (2) TITLEO, (AHV(I,1,1+lm*2+L),I=1,IM*JM)
C**** Write out specific humidity
      DO 450 L=1,lm
      WRITE (TITLEO(1:55),942) 'SPECIFIC HUMIDITY at',POUT(L)
  450 WRITE (2) TITLEO, (AHV(I,1,1+lm*3+L),I=1,IM*JM)
C**** Write out surface air temperature
      TITLEO(1:55) = 'SURFACE AIR TEMPERATURE (K)'
      WRITE (2) TITLEO, (AHV(I,1,1+lm*4+1),I=1,IM*JM)
C**** Write out sea level pressure
      TITLEO(1:55) = 'SEA LEVEL PRESSURE (mb)'
      WRITE (2) TITLEO, (AHV(I,1,1+lm*4+2),I=1,IM*JM)
      CLOSE (2)
C****
C**** Produce line printer maps of data
C****
C     IPL = NINT(974.*SIG(lm/2)+PTOP)
C     CALL MAP  (IM,JM,IPL,TITLEV(1),AHV(1,1,1+lm/2))
C     CALL MAP  (IM,JM,IPL,TITLEV(2),AHV(1,1,1+lm/2+lm))
C     CALL MAP1 (IM,JM,IPL,TITLEV(3),AHV(1,1,1+lm/2+2*lm),WTA,1.,
C    *           273.15,1)
C     CALL MAP1 (IM,JM,IPL,TITLEV(4),AHV(1,1,1+lm/2+3*lm),WTA,1.e5,
C                0.,1)
C     CALL MAP1 (IM,JM,0,TITLEV(5),AHV(1,1,2+4*lm),WTA,1.,273.15,1)
C     CALL MAP1 (IM,JM,-1000,TITLEV(6),AHV(1,1,3+4*lm),WTA,1.,1000.,1)
  942 FORMAT (A,F10.5,' (mb)')
      END

      SUBROUTINE VERT (IM,JM,GRID,DLATDG,PLEV,AIN,RHJ,TJ,UJ,ZATMO, JDAY,
     *                 LMA,LTM,SIGE,AOUT, PE)
C****
C**** VERT performs vertical interpolation of prognostic quantities
C**** from pressure surfaces to sigma surfaces.  The input arrays
C**** have been horizontally interpolated to the approprate grid.
C****
C**** Input: IM,JM = horizontal dimensions based on resolution
C****         GRID = grid arrangement for location of winds
C****         PLEV = NCAR pressure levels
C****   AIN( 1:12) = U = west-east wind component (m/s)
C****   AIN(13:24) = V = south-north wind component (m/s)
C****   AIN(25:36) = T = air temperature (K)
C****   AIN(37:48) = H = geopotential heights (m)
C****   AIN(49:54) = RH = relative humidity (%) (to 300mb only)
C****      AIN(55) = TSURF = surface air temperature (K)
C****      AIN(56) = PTROP = tropopause pressure (mb)
C****      AIN(57) = TTROP = tropopause temperature (K)
C****          RHJ = relative humidity (%)    250->0mb zonal climatology
C****       UTJ(1) = west-east wind  (m/s)     30->0mb zonal climatology
C****       UTJ(2) = air temperature (K)       30->0mb zonal climatology
C****        ZATMO = atmospheric topography (m)
C****          LMA = number of vertical layers
C****          LTM = last layer with sigma coordinates
C****         SIGE = sigma values for layer edges below LTM
C****           PE = constant pressure layer edges (mb) for LTM->LM
C****
C**** Output: AOUT = PS = surface pressure (mb)
C****            U(LMA) = west-east wind component (m/s)
C****            V(LMA) = south-north wind component (m/s)
C****            T(LMA) = air temperature (K)
C****            Q(LMA) = specific humidity (1)
C****             TSURF = surface temperature (K)
C****             PSLEV = sea level pressure (mb)
C****
C****             PE = pressure at upper layer edge (mb) 0->LTM-1
      REAL*4 PLEV(19),AIN(IM,JM,57),RHJ(JM,7:19),TJ(JM,13:19),
     *       ZATMO(IM,JM),RLAT(360),DXYU(360),XA(99),XB(99),
     *   UJ(JM,13:19),SIGE(0:LTM),AOUT(IM,JM,4*LMA+3)
      REAL*4 U(0:20),V(0:20),T(0:20),R(0:20),P(0:20),PE(0:LMA)
      REAL*8 BETA
      CHARACTER GRID*1
      QSAT(TM,PR) = .622*EXP(21.65604-5417.983/TM)/PR
      if(LMA.gt.99) stop 'increase dim. of xa,xb in vert to LMA'
C****
C**** Define function statements to locate prognostic variables in
C**** the output array, AOUT
C****
      LPS  =         1
      LU1  =         2
      LV1  = 1*LMA + 2
      LT0  = 2*LMA + 1
      LQ0  = 3*LMA + 1
      LTS  = 4*LMA + 2
      LPSL = 4*LMA + 3
C**** Define constant parameters
      KNCAR = 12
      KMAX  = 19
      GRAV  = 9.81
      RGAS  = 287.
      BETA0 = .0065
      GBYR  = GRAV/RGAS
      GBYRB = GBYR/BETA0
C**** Calculate spherical geometry
      TWOPI= 6.2831853
      DLON = TWOPI/IM
ccc   DLATDG = NINT(360./(JM-1))/2.
      DLAT = TWOPI*DLATDG/360.
      FJEQ = .5*(1+JM)
      SINS = -1.
      DO 20 J=1,JM-1
      SINN = SIN(DLAT*(J+.5-FJEQ))
      RLAT(J) = DLATDG*(J-FJEQ)
      DXYU(J) = DLON*(SINN-SINS)
   20 SINS = SINN
      RLAT(1) = -90.
      RLAT(JM) = 90.
      DXYU(1) = 2.*DXYU(1)
      DXYU(JM)= DXYU(1)
C****
C**** Create one dimensional arrays for pressure, temperature and
C**** humidity including surface values, pressure levels above the
C**** surface, and tropopause values.
C****
      DO 340 J=1,JM
      DO 340 I=1,IM
C**** Determine lowest input pressure level above the surface
      HS = ZATMO(I,J)
      DO 110 K=1,KNCAR
      K1=K
      IF(AIN(I,J,K+36).gt.HS)  GO TO 120
  110 continue
      WRITE (6,911) I,J,HS,AIN(I,J,48)
      STOP 'AIC: VERT: 110'
  120 CONTINUE
C**** Determine highest input pressure level below the tropopause
      PTROP = AIN(I,J,56)
      DO 130 K=KNCAR,1,-1
      KTM=K
      IF(PLEV(K).gt.PTROP)  GO TO 140
  130 continue
      WRITE (6,913) I,J,PTROP,PLEV(1)
      STOP 'AIC: VERT: 130'
  140 CONTINUE
C**** Load stratospheric values into local arrays
      DO 150 K=KTM+1,KMAX
      P(K+2-K1) = PLEV(K)
      IF(K.LE.KNCAR) T(K+2-K1) = AIN(I,J,K+24)
      IF(K.GT.KNCAR) T(K+2-K1) = TJ(J,K)
C**** To put Q=0 in stratosphere, set R... = 0. below (CA)
  150 R(K+2-K1) = 0. ! RHJ(J,K)
CA150 R(K+2-K1) = 0.
C**** Load tropopause values into local arrays
      P(KTM+2-K1) = PTROP
      T(KTM+2-K1) = AIN(I,J,57)
      CALL TSBOOK(JDAY,PTROP,RLAT(J),TTRPPC,QTROP)
      R(KTM+2-K1) = 0. ! 100.*QTROP/QSAT(TTRPPC,PTROP)
CALT  R(KTM+2-K1) = 0.
C**** Load tropospheric values into local arrays
      DO 160 K=K1,KTM
      P(K+1-K1) = PLEV(K)
      T(K+1-K1) = AIN(I,J,K+24)
      IF(K.le.6)  then
        R(K+1-K1) = AIN(I,J,K+48)
      else
        R(K+1-K1) = R(KTM+2-K1)+(AIN(I,J,54)-R(KTM+2-K1))*
     *              (PLEV(K)-PTROP)/(PLEV(6)-PTROP)
      endif
  160 continue
C**** Calculate surface pressure
      IF(K1.le.1)  GO TO 170
      BETA = (T(1)-AIN(I,J,K1+23)) / (AIN(I,J,K1+36)-AIN(I,J,K1+35))
      IF(ABS(BETA).lt.1.E-5)  GO TO 170
      P(0) = PLEV(K1)*(1.-BETA*(AIN(I,J,K1+36)-HS)/T(1))**(-GBYR/BETA)
      GO TO 180
  170 P(0) = PLEV(K1)*EXP(GBYR*(AIN(I,J,K1+36)-HS)/T(1))
  180 IF(PTROP.le..9*P(0))  GO TO 190
      WRITE (6,*) ' Tropopause pressure too high, PTROP,PS=',PTROP,P(0)
      STOP 180
C**** Calculate surface air temperature and surface relative humidity
  190 IF(K1.gt.1)  GO TO 200
      T(0) = AIN(I,J,25)
      R(0) = AIN(I,J,49)
      GO TO 300
  200 T(0) = T(1) + (AIN(I,J,K1+23)-T(1))*(P(0)-P(1))/(PLEV(K1-1)-P(1))
Calt* T(0) = AIN(I,J,55)  ! currently not used - not consistent with Qs
      IF(K1.le.7)  then
        R(0) = R(1)+(AIN(I,J,K1+47)-R(1))*(P(0)-P(1))/(PLEV(K1-1)-P(1))
      else
        R(0) = R(1)+(RHJ(J,K1-1)-R(1))*(P(0)-P(1))/(PLEV(K1-1)-P(1))
      endif
C****
C**** Calculate output values defined on the A grid
C****
C**** Surface pressure, surface air temperature, sea level pressure
C****
  300 AOUT(I,J,LPS)  = P(0)
      AOUT(I,J,LTS)  = T(0)
      AOUT(I,J,LPSL) = P(0)*(1.+BETA0*ZATMO(I,J)/T(0))**GBYRB
C****
C**** Integrate temperature and humidity over sigma layers
C****
      DO 310 L=0,LTM-1
  310 PE(L) = SIGE(L)*(P(0)-PE(LTM)) + PE(LTM)
      CALL VNTRP1 (KMAX+2-K1,P,T, LMA,PE,XA)
      CALL VNTRP1 (KMAX+2-K1,P,R, LMA,PE,XB)
C**** Convert relative to specific humidity, save T and RH in AOUT
      DO 320 L=1,LMA
      PL = .5*(PE(L-1)+PE(L))
      AOUT(I,J,LT0+L) = XA(L)
  320 AOUT(I,J,LQ0+L) = max(3.e-6,.01*XB(L)*QSAT(XA(L),PL))
  340 continue
C****
C**** Winds are defined on the B grid
C****
      DO 400 K=1,KMAX
  400 P(K)=PLEV(K)
      IF(GRID.ne.'B')  GO TO 500
      I=IM
      DO 480 J=1,JM-1
      DO 480 IP1=1,IM
      P(0) = .5*((AOUT(I,J  ,LPS)+AOUT(IP1,J  ,LPS))*DXYU(J  ) +
     +                (AOUT(I,J+1,LPS)+AOUT(IP1,J+1,LPS))*DXYU(J+1)) /
     /            (DXYU(J)+DXYU(J+1))
      U(0)=0.
      V(0)=0.
      DO 410 K=1,KMAX
      IF(K.LE.KNCAR) THEN
         U(K)=AIN(I,J,K)
         V(K)=AIN(I,J,K+KNCAR)
      ELSE
         U(K)=UJ(J,K)
         V(K)=0.
      END IF
  410 CONTINUE
      DO 420 L=0,LTM-1
  420 PE(L) = SIGE(L)*(P(0)-PE(LTM)) + PE(LTM)
      CALL VNTRP1 (KMAX,P,U, LMA,PE,XA)
      CALL VNTRP1 (KMAX,P,V, LMA,PE,XB)
      DO 430 L=1,LMA
      AOUT(I,J+1,LU1-1+L) = XA(L)
      AOUT(I,J+1,LV1-1+L) = XB(L)
CNEWG AOUT(I,J,LU1-1+L) = XA(L)
CNEWG AOUT(I,J,LV1-1+L) = XB(L)
  430 CONTINUE
  480 I=IP1
C**** Zero out undefined winds
      DO 490 L=LU1,LT0
      DO 490 I=1,IM
      AOUT(I,1,L) = 0.
CNEWG AOUT(I,JM,L) = 0.
  490 CONTINUE
      RETURN
C****
C**** Winds are defined on the C grid
C****
C**** U component
  500 I=IM
      DO 550 J=1,JM
      DO 550 IP1=1,IM
      P(0) = .5*(AOUT(I,J,LPS)+AOUT(IP1,J,LPS))
      U(0)=0.
      DO 510 K=1,KNCAR
  510 U(K)=AIN(I,J,K)
      DO 515 K=KNCAR+1,KMAX
  515 U(K)=UJ(J,K)
      DO 520 L=0,LTM-1
  520 PE(L) = SIGE(L)*(P(0)-PE(LTM)) + PE(LTM)
      CALL VNTRP1 (KMAX,P,U, LMA,PE,XA)
      DO 530 L=1,LMA
  530 AOUT(I,J,LU1-1+L) = XA(L)
  550 I=IP1
C**** V component
      DO 580 J=1,JM-1
      DO 580 I=1,IM
      P(0) = (AOUT(I,J,LPS)*DXYU(J)+AOUT(I,J+1,LPS)*DXYU(J+1)) /
     /            (DXYU(J)+DXYU(J+1))
      V(0)=0.
      DO 560 K=1,KNCAR
  560 V(K)=AIN(I,J,K+KNCAR)
      DO 565 K=KNCAR+1,KMAX
  565 V(K)=0.
      DO 570 L=0,LTM-1
  570 PE(L) = SIGE(L)*(P(0)-PE(LTM)) + PE(LTM)
      CALL VNTRP1 (KMAX,P,V, LMA,PE,XA)
      DO 580 L=1,LMA
  580 AOUT(I,J,LV1-1+L) = XA(L)
C**** Zero out undefined winds
      DO 590 L=LV1,LT0
      DO 590 I=1,IM
  590 AOUT(I,JM,L) = 0.
      RETURN
C****
  911 FORMAT ('0Surface topography exceeds highest geopotential ',
     *  'height in subroutine VERT.'/' IM,JM,HS,AIN(56)=',2I4,2E16.6)
  913 FORMAT ('0Tropopause pressure exceeds first pressure level.'/
     *  ' IM,JM,AIN(64),PLEV(1)=',2I4,2E16.6)
      END

      SUBROUTINE VNTRP1 (KM,P,AIN,  LMA,PE,AOUT)
C**** Vertically interpolates a 1-D array
C**** Input:       KM = number of input pressure levels
C****            P(K) = input pressure levels (mb)
C****          AIN(K) = input quantity at level P(K)
C****             LMA = number of vertical layers of output grid
C****           PE(L) = output pressure levels (mb) (edges of layers)
C**** Output: AOUT(L) = output quantity: mean between PE(L-1) & PE(L)
C****
      REAL*4 P(0:KM),AIN(0:KM),    PE(0:LMA),AOUT(LMA)
C****
      PDN = PE(0)
      ADN = AIN(0)
      K=1
C**** Ignore input levels below ground level pe(0)=p(0)
      IF(P(1).GT.PE(0)) THEN
         DO K1=2,KM
         K=K1
         IF(P(K).LT.PE(0)) THEN  ! interpolate to ground level
           ADN=AIN(K)+(AIN(K-1)-AIN(K))*(PDN-P(K))/(P(K-1)-P(K))
           GO TO 300
         END IF
         END DO
         STOP 'VNTRP1 - error - should not get here'
      END IF
C**** Integrate - connecting input data by straight lines
  300 DO 330 L=1,LMA
      ASUM = 0.
      PSUM = 0.
      PUP = PE(L)
  310 IF(P(K).le.PUP)  GO TO 320
      PSUM = PSUM + (PDN-P(K))
      ASUM = ASUM + (PDN-P(K))*(ADN+AIN(K))/2.
      PDN  = P(K)
      ADN  = AIN(K)
      K=K+1
      IF(K.LE.KM) GO TO 310
      stop 'VNTRP1 - should not happen'
C****
  320 AUP  = AIN(K) + (ADN-AIN(K))*(PUP-P(K))/(PDN-P(K))
      PSUM = PSUM + (PDN-PUP)
      ASUM = ASUM + (PDN-PUP)*(ADN+AUP)/2.
      AOUT(L) = ASUM/PSUM
      PDN = PUP
  330 ADN = AUP
C****
      RETURN
      END

      SUBROUTINE STRDAT (RHJ,TJ,UJ,JM,JDAY,GRID,DLATDG,PLEV)
C****
C**** STRDAT uses observed stratospheric data to fill in an initial
C**** conditions array above 50 mb, and above 300mb for rel.humidity
C****
C**** Input:  IM,JM   = longitude and latitude dimensions
C****         JDAY    = Julian date of request
C****         GRID    = grid arrangement for velocity components
C****         PLEV    = pressure levels starting at 250mb
C**** Output:  RHJ    = relative humidity (%)  250mb->0mb
C****           TJ    = temperature (C)         30mb->0mb
C****           UJ    = zonal wind (m/s)        30mb->0mb
C****
      REAL*4 RHJ(JM,13),TJ(JM,7),UJ(JM,7),PLEV(13)
      CHARACTER*1 GRID
      QSAT(TM,PR) = .622*EXP(21.65604-5417.983/TM)/PR
C****
      FJEQ = .5*(1+JM)
      DO 210 LS=1,6
      PRES = PLEV(LS)
      DO 210 J=1,JM
C**** Fill in upper-tropospheric rel.humidity from observations
      XLAT= DLATDG*(J-FJEQ)
      CALL TSBOOK (JDAY,PRES,XLAT,TBK,SBK)
  210 RHJ(J,LS)=100.*SBK/QSAT(TBK,PRES)
C****
      DO 250 LS=7,13
      PRES = PLEV(LS)
      UJ(JM,LS-6)=0.
      DO 250 J=1,JM
C**** Fill in stratospheric temperatures and humidity from climatology
      XLAT=DLATDG*(J-FJEQ)
      CALL TSBOOK (JDAY,PRES,XLAT,TBK,SBK)
      RHJ(J,LS)=100.*SBK/QSAT(TBK,PRES)
      TJ(J,LS-6)=TBK
C**** Fill in stratospheric zonal wind from observations
      IF(GRID.ne.'B')  GO TO 220
      IF(J.eq.JM)  GO TO 250
      XLAT = XLAT + .5*DLATDG
  220 CALL UBOOK (JDAY,PRES,XLAT,UBK)
      UJ(J,LS-6)=UBK
  250 continue
      RETURN
      END

      SUBROUTINE UBOOK (JDAY,P,XLAT,UU)
C****
C**** UBOOK returns observed zonal velocities in the stratosphere
C**** Input: JDAY = Julian day of the year
C****        P    = pressure (mb)
C****        XLAT = latitude (degrees)
C**** Output: UU  = zonal velocity (m/s)
C****
      PARAMETER (JM=19,LM=8)
      REAL*4 U(JM,LM)
      COMMON /NU19X8/ NU(JM,LM)
      save COSDAY,U
      DATA PI/3.14159265/, JDAY0/-1/
C****
      IF(JDAY.eq.JDAY0)  GO TO 200
      JDAY0 = JDAY
      COSDAY = COS(2.*PI*(JDAY-16)/365.)
      DO 110 L=1,LM
      DO 110 J=1,JM
      JHX=JM+1-J
  110 U(J,L) = .5*(NU(J,L)+NU(JHX,L)) + .5*(NU(J,L)-NU(JHX,L))*COSDAY
C****
  200 XJ = (XLAT+90.)/10.+1.
      J = XJ
      IF(J.ge.JM)  J = JM-1
      YL = 4.+2.*ALOG10(P+1.e-6)
      IF(YL.LT.1.)  YL = 1.
      L = YL
      IF(L.ge.LM)  L = LM-1
      UU = (U(J,L  )*(J+1-XJ) + U(J+1,L  )*(XJ-J))*(L+1-YL)
     *   + (U(J,L+1)*(J+1-XJ) + U(J+1,L+1)*(XJ-J))*(YL-L)
      RETURN
      END

      SUBROUTINE TSBOOK (JDAY,P,XLAT,T,Q)
C****
C**** TSBOOK returns observed temperature and humidity in the
C**** stratosphere from the McClatchey atlas
C**** Input: JDAY = Julian day of the year
C****        P    = pressure (mb)
C****        XLAT = latitude (degrees)
C**** Output: T = temperature (K)
C****         Q = specific humidity (1)
C****
      DATA PI/3.141592653/
C****
      PHASEX = PI*(JDAY-16-91.25)/182.5
      SEASON = .5*SIN(PHASEX)
      IF(XLAT.lt.0.)  SEASON = -SEASON
      WINTER = .5-SEASON
      SUMMER = .5+SEASON
      DEGL45 = 0.
      ABSLAT = ABS(XLAT)
      IF(ABSLAT.gt.15.)  DEGL45 = (ABSLAT-15.)/30.
      IF(ABSLAT.gt.45.)  DEGL45 = (60.-ABSLAT)/15.
      IF(DEGL45.lt.-1.)  DEGL45 = -1.
      DEGLXX = 1.-DEGL45
      NATMW=3
      NATMS=2
      CALL PHDATM (P,NATMW,TW,QW)
      CALL PHDATM (P,NATMS,TS,QS)
      T = TW*WINTER + TS*SUMMER
      Q = QW*WINTER + QS*SUMMER
      NATMW=5
      NATMS=4
      IF(ABSLAT.lt.45.)  NATMS=1
      IF(ABSLAT.lt.45.)  NATMW=1
      CALL PHDATM (P,NATMW,TW,QW)
      CALL PHDATM (P,NATMS,TS,QS)
      IF(ABSLAT.gt.60.)  GO TO 100
      T = T*DEGL45 + DEGLXX*(TW*WINTER+TS*SUMMER)
      Q = Q*DEGL45 + DEGLXX*(QW*WINTER+QS*SUMMER)
      GO TO 110
 100  continue
      WWSS = TW*WINTER + TS*SUMMER
      T    = EXP(DEGL45*ALOG(T)+DEGLXX*ALOG(WWSS))
      WWSS = QW*WINTER + QS*SUMMER
      Q    = EXP(DEGL45*ALOG(Q)+DEGLXX*ALOG(WWSS))
 110  continue
      RETURN
      END

      SUBROUTINE PHDATM (P,N,T,Q)
C****
C**** PHDATM returns temperature and humidity for a given pressure
C**** and location using the McClatchy (1972) atmospheric data
C****
C**** Input:  P = pressure (mb)
C****         N = 1  GIVES ATMOSPHERE DATA FOR  TROPICAL LATITUDES
C****             2  GIVES ATMOSPHERE DATA FOR  MIDLATITUDE SUMMER
C****             3  GIVES ATMOSPHERE DATA FOR  MIDLATITUDE WINTER
C****             4  GIVES ATMOSPHERE DATA FOR  SUBARCTIC SUMMER
C****             5  GIVES ATMOSPHERE DATA FOR  SUBARCTIC WINTER
C**** Output: T = TEMPERATURE (K)
C****         Q = SPECIFIC HUMIDITY (GRAMS WATER VAPOR)/(GRAMS AIR)
C****
      COMMON /CDAT/ PRES(33,5),DENS(33,5),TEMP(33,5),WVAP(33,5)
C****
      IF(P.ge.PRES(1,N))  GO TO 30
      DO 10 K1=2,33
      K=K1
      IF(P.ge.PRES(K,N))  GO TO 20
   10 continue
C**** Input pressure is less than minimum pressure level
      T = TEMP(33,N)
      Q = WVAP(33,N)/DENS(33,N)
      RETURN
C**** Input pressure is between pressure levels K-1 and K
   20 FRAC = (P-PRES(K,N)) / (PRES(K-1,N)-PRES(K,N))
      T = FRAC*TEMP(K-1,N) + (1.-FRAC)*TEMP(K,N)
      Q = FRAC*WVAP(K-1,N)/DENS(K-1,N) + (1.-FRAC)*WVAP(K,N)/DENS(K,N)
      RETURN
C**** Input pressure exceeds maximum pressure level
   30 T = TEMP(1,N)
      Q = WVAP(1,N)/DENS(1,N)
      RETURN
      END

      BLOCK DATA
C*****
C***** McClatchey (1972) atmospheric data
C*****
      COMMON /CDAT/ PRS1(33),PRS2(33),PRS3(33),PRS4(33),PRS5(33),
     1              DNS1(33),DNS2(33),DNS3(33),DNS4(33),DNS5(33),
     2              TMP1(33),TMP2(33),TMP3(33),TMP4(33),TMP5(33),
     3              WVP1(33),WVP2(33),WVP3(33),WVP4(33),WVP5(33)
      COMMON /NU19X8/ NU(19,8)
C****
C*1** TROPICAL LATITUDES     MCCLATCHEY (1972) ATMOSPHERE DATA VS HEIGHT
C****
      DATA PRS1/  1.013E 3, 9.040E 2, 8.050E 2, 7.150E 2, 6.330E 2,
     *  5.590E 2, 4.920E 2, 4.320E 2, 3.780E 2, 3.290E 2, 2.860E 2,
     *  2.470E 2, 2.130E 2, 1.820E 2, 1.560E 2, 1.320E 2, 1.110E 2,
     *  9.370E 1, 7.890E 1, 6.660E 1, 5.650E 1, 4.800E 1, 4.090E 1,
     *  3.500E 1, 3.000E 1, 2.570E 1, 1.220E 1, 6.000E 0, 3.050E 0,
     *  1.590E 0, 8.540E-1, 5.790E-2, 3.000E-4/
      DATA DNS1/  1.167E 3, 1.064E 3, 9.689E 2, 8.756E 2, 7.951E 2,
     *  7.199E 2, 6.501E 2, 5.855E 2, 5.258E 2, 4.708E 2, 4.202E 2,
     *  3.740E 2, 3.316E 2, 2.929E 2, 2.578E 2, 2.260E 2, 1.972E 2,
     *  1.676E 2, 1.382E 2, 1.145E 2, 9.515E 1, 7.938E 1, 6.645E 1,
     *  5.618E 1, 4.763E 1, 4.045E 1, 1.831E 1, 8.600E 0, 4.181E 0,
     *  2.097E 0, 1.101E 0, 9.210E-2, 5.000E-4/
      DATA TMP1/    300.0,    294.0,    288.0,    284.0,    277.0,
     *    270.0,    264.0,    257.0,    250.0,    244.0,    237.0,
     *    230.0,    224.0,    217.0,    210.0,    204.0,    197.0,
     *    195.0,    199.0,    203.0,    207.0,    211.0,    215.0,
     *    217.0,    219.0,    221.0,    232.0,    243.0,    254.0,
     *    265.0,    270.0,    219.0,    210.0/
      DATA WVP1/    1.9E 1,   1.3E 1,   9.3E 0,   4.7E 0,   2.2E 0,
     *    1.5E 0,   8.5E-1,   4.7E-1,   2.5E-1,   1.2E-1,   5.0E-2,
     *    1.7E-2,   6.0E-3,   1.8E-3,   1.0E-3,   7.6E-4,   6.4E-4,
     *    5.6E-4,   5.0E-4,   4.9E-4,   4.5E-4,   5.1E-4,   5.1E-4,
     *    5.4E-4,   6.0E-4,   6.7E-4,   3.6E-4,   1.1E-4,   4.3E-5,
     *    1.9E-5,   6.3E-6,   1.4E-7,   1.0E-9/
C****
C*2** MIDLATITUDE SUMMER     MCCLATCHEY (1972) ATMOSPHERE DATA VS HEIGHT
C****
      DATA PRS2/  1.013E 3, 9.020E 2, 8.020E 2, 7.100E 2, 6.280E 2,
     *  5.540E 2, 4.870E 2, 4.260E 2, 3.720E 2, 3.240E 2, 2.810E 2,
     *  2.430E 2, 2.090E 2, 1.790E 2, 1.530E 2, 1.300E 2, 1.110E 2,
     *  9.500E 1, 8.120E 1, 6.950E 1, 5.950E 1, 5.100E 1, 4.370E 1,
     *  3.760E 1, 3.220E 1, 2.770E 1, 1.320E 1, 6.520E 0, 3.330E 0,
     *  1.760E 0, 9.510E-1, 6.710E-2, 3.000E-4/
      DATA DNS2/  1.191E 3, 1.080E 3, 9.757E 2, 8.846E 2, 7.998E 2,
     *  7.211E 2, 6.487E 2, 5.830E 2, 5.225E 2, 4.669E 2, 4.159E 2,
     *  3.693E 2, 3.269E 2, 2.882E 2, 2.464E 2, 2.104E 2, 1.797E 2,
     *  1.535E 2, 1.305E 2, 1.110E 2, 9.453E 1, 8.056E 1, 6.872E 1,
     *  5.867E 1, 5.014E 1, 4.288E 1, 1.322E 1, 6.519E 0, 3.330E 0,
     *  1.757E 0, 9.512E-1, 6.706E-2, 5.000E-4/
      DATA TMP2/    294.0,    290.0,    285.0,    279.0,    273.0,
     *    267.0,    261.0,    255.0,    248.0,    242.0,    235.0,
     *    229.0,    222.0,    216.0,    216.0,    216.0,    216.0,
     *    216.0,    216.0,    217.0,    218.0,    219.0,    220.0,
     *    222.0,    223.0,    224.0,    234.0,    245.0,    258.0,
     *    270.0,    276.0,    218.0,    210.0/
      DATA WVP2/    1.4E 1,   9.3E 0,   5.9E 0,   3.3E 0,   1.9E 0,
     *    1.0E 0,   6.1E-1,   3.7E-1,   2.1E-1,   1.2E-1,   6.4E-2,
     *    2.2E-2,   6.0E-3,   1.8E-3,   1.0E-3,   7.6E-4,   6.4E-4,
     *    5.6E-4,   5.0E-4,   4.9E-4,   4.5E-4,   5.1E-4,   5.1E-4,
     *    5.4E-4,   6.0E-4,   6.7E-4,   3.6E-4,   1.1E-4,   4.3E-5,
     *    1.9E-5,   6.3E-6,   1.4E-7,   1.0E-9/
C****
C*3** MIDLATITUDE WINTER     MCCLATCHEY (1972) ATMOSPHERE DATA VS HEIGHT
C****
      DATA PRS3/  1.018E 3, 8.973E 2, 7.897E 2, 6.938E 2, 6.081E 2,
     *  5.313E 2, 4.627E 2, 4.016E 2, 3.473E 2, 2.992E 2, 2.568E 2,
     *  2.199E 2, 1.882E 2, 1.610E 2, 1.378E 2, 1.178E 2, 1.007E 2,
     *  8.610E 1, 7.350E 1, 6.280E 1, 5.370E 1, 4.580E 1, 3.910E 1,
     *  3.340E 1, 2.860E 1, 2.430E 1, 1.110E 1, 5.180E 0, 2.530E 0,
     *  1.290E 0, 6.820E-1, 4.670E-2, 3.000E-4/
      DATA DNS3/  1.301E 3, 1.162E 3, 1.037E 3, 9.230E 2, 8.282E 2,
     *  7.411E 2, 6.614E 2, 5.886E 2, 5.222E 2, 4.619E 2, 4.072E 2,
     *  3.496E 2, 2.999E 2, 2.572E 2, 2.206E 2, 1.890E 2, 1.620E 2,
     *  1.388E 2, 1.188E 2, 1.017E 2, 8.690E 1, 7.421E 1, 6.338E 1,
     *  5.415E 1, 4.624E 1, 3.950E 1, 1.783E 1, 7.924E 0, 3.625E 0,
     *  1.741E 0, 8.954E-1, 7.051E-2, 5.000E-4/
      DATA TMP3/    272.2,    268.7,    265.2,    261.7,    255.7,
     *    249.7,    243.7,    237.7,    231.7,    225.7,    219.7,
     *    219.2,    218.7,    218.2,    217.7,    217.2,    216.7,
     *    216.2,    215.7,    215.2,    215.2,    215.2,    215.2,
     *    215.2,    215.2,    215.2,    217.4,    227.8,    243.2,
     *    258.5,    265.7,    230.7,    210.2/
      DATA WVP3/    3.5E 0,   2.5E 0,   1.8E 0,   1.2E 0,   6.6E-1,
     *    3.8E-1,   2.1E-1,   8.5E-2,   3.5E-2,   1.6E-2,   7.5E-3,
     *    6.9E-3,   6.0E-3,   1.8E-3,   1.0E-3,   7.6E-4,   6.4E-4,
     *    5.6E-4,   5.0E-4,   4.9E-4,   4.5E-4,   5.1E-4,   5.1E-4,
     *    5.4E-4,   6.0E-4,   6.7E-4,   3.6E-4,   1.1E-4,   4.3E-5,
     *    1.9E-5,   6.3E-6,   1.4E-7,   1.0E-9/
C****
C*4** SUBARCTIC SUMMER       MCCLATCHEY (1972) ATMOSPHERE DATA VS HEIGHT
C****
      DATA PRS4/  1.010E 3, 8.960E 2, 7.929E 2, 7.000E 2, 6.160E 2,
     *  5.410E 2, 4.730E 2, 4.130E 2, 3.590E 2, 3.107E 2, 2.677E 2,
     *  2.300E 2, 1.977E 2, 1.700E 2, 1.460E 2, 1.250E 2, 1.080E 2,
     *  9.280E 1, 7.980E 1, 6.860E 1, 5.890E 1, 5.070E 1, 4.360E 1,
     *  3.750E 1, 3.227E 1, 2.780E 1, 1.340E 1, 6.610E 0, 3.400E 0,
     *  1.810E 0, 9.870E-1, 7.070E-2, 3.000E-4/
      DATA DNS4/  1.220E 3, 1.110E 3, 9.971E 2, 8.985E 2, 8.077E 2,
     *  7.244E 2, 6.519E 2, 5.849E 2, 5.231E 2, 4.663E 2, 4.142E 2,
     *  3.559E 2, 3.059E 2, 2.630E 2, 2.260E 2, 1.943E 2, 1.671E 2,
     *  1.436E 2, 1.235E 2, 1.062E 2, 9.128E 1, 7.849E 1, 6.750E 1,
     *  5.805E 1, 4.963E 1, 4.247E 1, 1.338E 1, 6.614E 0, 3.404E 0,
     *  1.817E 0, 9.868E-1, 7.071E-2, 5.000E-4/
      DATA TMP4/    287.0,    282.0,    276.0,    271.0,    266.0,
     *    260.0,    253.0,    246.0,    239.0,    232.0,    225.0,
     *    225.0,    225.0,    225.0,    225.0,    225.0,    225.0,
     *    225.0,    225.0,    225.0,    225.0,    225.0,    225.0,
     *    225.0,    226.0,    228.0,    235.0,    247.0,    262.0,
     *    274.0,    277.0,    216.0,    210.0/
      DATA WVP4/    9.1E 0,   6.0E 0,   4.2E 0,   2.7E 0,   1.7E 0,
     *    1.0E 0,   5.4E-1,   2.9E-1,   1.3E-2,   4.2E-2,   1.5E-2,
     *    9.4E-3,   6.0E-3,   1.8E-3,   1.0E-3,   7.6E-4,   6.4E-4,
     *    5.6E-4,   5.0E-4,   4.9E-4,   4.5E-4,   5.1E-4,   5.1E-4,
     *    5.4E-4,   6.0E-4,   6.7E-4,   3.6E-4,   1.1E-4,   4.3E-5,
     *    1.9E-5,   6.3E-6,   1.4E-7,   1.0E-9/
C****
C*5** SUBARCTIC WINTER       MCCLATCHEY (1972) ATMOSPHERE DATA VS HEIGHT
C****
      DATA PRS5/  1.013E 3, 8.878E 2, 7.775E 2, 6.798E 2, 5.932E 2,
     *  5.158E 2, 4.467E 2, 3.853E 2, 3.308E 2, 2.829E 2, 2.418E 2,
     *  2.067E 2, 1.766E 2, 1.510E 2, 1.291E 2, 1.103E 2, 9.431E 1,
     *  8.058E 1, 6.882E 1, 5.875E 1, 5.014E 1, 4.277E 1, 3.647E 1,
     *  3.109E 1, 2.649E 1, 2.256E 1, 1.020E 1, 4.701E 0, 2.243E 0,
     *  1.113E 0, 5.719E-1, 4.016E-2, 3.000E-4/
      DATA DNS5/  1.372E 3, 1.193E 3, 1.058E 3, 9.366E 2, 8.339E 2,
     *  7.457E 2, 6.646E 2, 5.904E 2, 5.226E 2, 4.538E 2, 3.879E 2,
     *  3.315E 2, 2.834E 2, 2.422E 2, 2.071E 2, 1.770E 2, 1.517E 2,
     *  1.300E 2, 1.113E 2, 9.529E 1, 8.155E 1, 6.976E 1, 5.966E 1,
     *  5.100E 1, 4.358E 1, 3.722E 1, 1.645E 1, 7.368E 0, 3.330E 0,
     *  1.569E 0, 7.682E-1, 5.695E-2, 5.000E-4/
      DATA TMP5/    257.1,    259.1,    255.9,    252.7,    247.7,
     *    240.9,    234.1,    227.3,    220.6,    217.2,    217.2,
     *    217.2,    217.2,    217.2,    217.2,    217.2,    216.6,
     *    216.0,    215.4,    214.8,    214.1,    213.6,    213.0,
     *    212.4,    211.8,    211.2,    216.0,    222.2,    234.7,
     *    247.0,    259.3,    245.7,    210.0/
      DATA WVP5/    1.2E 0,   1.2E 0,   9.4E-1,   6.8E-1,   4.1E-1,
     *    2.0E-1,   9.8E-2,   5.4E-2,   1.1E-2,   8.4E-3,   5.5E-3,
     *    3.8E-3,   2.6E-3,   1.8E-3,   1.0E-3,   7.6E-4,   6.4E-4,
     *    5.6E-4,   5.0E-4,   4.9E-4,   4.5E-4,   5.1E-4,   5.1E-4,
     *    5.4E-4,   6.0E-4,   6.7E-4,   3.6E-4,   1.1E-4,   4.3E-5,
     *    1.9E-5,   6.3E-6,   1.4E-7,   1.0E-9/
C****
C**** Zonal winds for January,    every 10 degrees in latitude and
C**** pressures .03, .1, .3, 1, 3, 10, 32, 100 (mb).  Source:
C**** Barnett and Corney (1985) from 20-70 and Groves (1971) in
C**** the tropics (check with David Rind).
C****
      DATA NU/ -5,-15,-30,-50,-58,-54,-44,-28,-12,
     *      3, 20, 43, 62, 73, 72, 62, 29, 18,  8,
     *        -10,-25,-44,-60,-62,-55,-47,-35,-20,
     *     -3, 15, 40, 66, 80, 80, 55, 25,  9,  4,
     *         -9,-18,-32,-45,-53,-53,-51,-42,-26,
     *    -10,  5, 30, 55, 74, 71, 40, 10,  4, -2,
     *         -7,-10,-20,-30,-37,-43,-43,-41,-35,
     *    -25,-10, 15, 40, 52, 45, 21,  6, -1, -4,
     *         -5, -2,-13,-21,-26,-30,-35,-37,-34,
     *    -28,-25,  0, 27, 41, 35, 25,  6,  0, -4,
     *         -4, -4, -6,-10,-15,-21,-28,-33,-33,
     *    -28,-20, -7, 10, 17, 25, 25, 15, 11,  5,
     *         -2, -0,  1,  1, -3,-12,-22,-27,-27,
     *    -23,-13, -3,  4, 10, 20, 23, 20, 12,  5,
     *         -1,  2,  4,  6,  8,  5, -4,-10,-14,
     *     -6,  5, 20, 24, 22, 15, 16, 12,  5,  2/
      END
