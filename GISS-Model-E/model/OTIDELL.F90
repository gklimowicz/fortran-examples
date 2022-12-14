!**** OTIDELL.F90     Fortran source code for Ocean Tides on Lat-Lon grid     2015/02/12

      Module OTIDE
!****
!**** In the Earth's rectangular geocentric coordinate system (X,Y,Z),
!**** X-axis passes through the International Date Line and the Equator,
!**** Y-axis passes through 90 W and the Equator, and
!**** Z-axis passes through the North Pole.
!****
!**** In the local topographic coordinate system (U,V,W),
!**** U-axis points eastward,
!**** V-axis points northward, and
!**** W-axis points upward (away from center of Earth)
!****
      Real*8,Parameter :: &
         EDAYzY = 365, &      !  number of days in tropical year (should be 365.2425 for Gregorian calendar)
         VE2000 = 79.5, &     !  days from 2000/01/01/00 to Vernal Equinos of year 2000
         GUNIV = 6.672d-11, & !  universal gravity coef. (m^3/kg*s^2)
         ASEMI = 1.4960d11, & !  semi-major axis of Earth's orbit (m)
         MMOON = 7.35028d22,& !  mass of Moon (kg)
         MSUN  = 1.98941d30   !  mass of Sun (kg)
      Real*8 :: &
         XMC,YMC,ZMC, &           !  location of Moon in (X,Y,Z) coordinates (m)
         XSC,YSC,ZSC, &           !  location of Sun in (X,Y,Z) coordinates (m)
         XMCzD3,YMCzD3,ZMCzD3, &  !  = (XMC,YMC,ZMC) / (distance to Moon)^3
         XSCzD3,YSCzD3,ZSCzD3     !  = (XSC,YSC,ZSC) / (distance to Sun)^3
      Real*8,Allocatable,Dimension(:,:) :: &
         WTIDE, &           !  vertical perturbation of gravity (m/s^2)
         XUIJ,YUIJ,ZUIJ, &  !  (X,Y,Z) coordinates of C-grid UIJ on Earth's surface
         XVIJ,YVIJ,ZVIJ, &  !  (X,Y,Z) coordinates of C-grid VIJ on Earth's surface
         XAIJ,YAIJ,ZAIJ, &  !  (X,Y,Z) coordinates of center AIJ on Earth's surface
         XUC,YUC,ZUC,    &  !  Unit vector of U in (X,Y,Z) coordinates at C-grid U cells
         XVC,YVC,ZVC,    &  !  Unit vector of V in (X,Y,Z) coordinates at C-grid V cells
         XUD,YUD,ZUD,    &  !  Unit vector of U in (X,Y,Z) coordinates at D-grid U cells
         XVD,YVD,ZVD,    &  !  Unit vector of V in (X,Y,Z) coordinates at D-grid V cells
         XWA,YWA,ZWA        !  Unit vector of W in (X,Y,Z) coordinates at center A cells
      Real*8,Allocatable,Dimension(:) :: &
         XSTN,YSTN,ZSTN, &  !  (X,Y,Z) coordinates of center of strait
         XUST,YUST,ZUST     !  Unit vector in (X,Y,Z) coordinates of strait directions
      EndModule OTIDE


      Subroutine OTIDE0
!**** OTIDE0 calculates grid cell locations on Earth's surface
      Use OTIDE
      Use CONSTANT,   Only: RADIUS
      Use OCEANR_DIM, Only: OGRID
      Use OCEANRES,   Only: IMO,JMO
      Use OCEAN,      Only: SINI=>SINIC,COSI=>COSIC,SINU,COSU, DLON,DLAT,FJEQ, &
                            SINP=>SINPO,COSP=>COSPO,SINV=>SINVO,COSV=>COSVO
      Use STRAITS,    Only: NMST, IST,JST, XST,YST
      Use DOMAIN_DECOMP_1D, Only: AM_I_ROOT
      Implicit None
      Integer :: I,J,N, J1,JN,JNP,J1H,JNH
      Real*8  :: RLON1,RLON2,RLAT1,RLAT2, dX,dY,dZ,D

!**** Extract domain decomposition
      J1 = OGRID%J_STRT  ;  JN = OGRID%J_STOP  ;  JNP = Min(JN,JMO-1)  ;  J1H = Max(J1-1,1)  ;  JNH = Min(JN+1,JMO)

!**** Allocate space
      Allocate (WTIDE(IMO,J1H:JNH), &
                 XUIJ(IMO,J1H:JNH),YUIJ(IMO,J1H:JNH),ZUIJ(IMO,J1H:JNH), &
                 XVIJ(IMO,J1H:JNH),YVIJ(IMO,J1H:JNH),ZVIJ(IMO,J1H:JNH), &
                 XAIJ(IMO,J1H:JNH),YAIJ(IMO,J1H:JNH),ZAIJ(IMO,J1H:JNH), &
                  XUC(IMO,J1H:JNH), YUC(IMO,J1H:JNH), ZUC(IMO,J1H:JNH), &
                  XVC(IMO,J1H:JNH), YVC(IMO,J1H:JNH), ZVC(IMO,J1H:JNH), &
                  XUD(IMO,J1H:JNH), YUD(IMO,J1H:JNH), ZUD(IMO,J1H:JNH), &
                  XVD(IMO,J1H:JNH), YVD(IMO,J1H:JNH), ZVD(IMO,J1H:JNH), &
                  XWA(IMO,J1H:JNH), YWA(IMO,J1H:JNH), ZWA(IMO,J1H:JNH), &
                  XSTN(NMST),YSTN(NMST),ZSTN(NMST), &
                  XUST(NMST),YUST(NMST),ZUST(NMST))

!**** Geocentric coordinates on Earth's surface
      Do 10 J=J1,JN  ;  Do 10 I=1,IMO
      XUIJ(I,J) = RADIUS*COSP(J)*COSU(I)
      YUIJ(I,J) = RADIUS*COSP(J)*SINU(I)
      ZUIJ(I,J) = RADIUS*SINP(J)
      XAIJ(I,J) = RADIUS*COSP(J)*COSI(I)
      YAIJ(I,J) = RADIUS*COSP(J)*SINI(I)
   10 ZAIJ(I,J) = RADIUS*SINP(J)
      Do 20 J=J1,JNP  ;  Do 20 I=1,IMO
      XVIJ(I,J) = RADIUS*COSV(J)*COSI(I)
      YVIJ(I,J) = RADIUS*COSV(J)*SINI(I)
   20 ZVIJ(I,J) = RADIUS*SINV(J)

!**** Compute velocity unit vectors in geocentric coordinates
      Do 30 J=J1,JN  ;  Do 30 I=1,IMO
      XUC(I,J) = - SINU(I)  ;  XVD(I,J) = - SINP(J)*COSU(I)  ;  XWA(I,J) = COSP(J)*COSI(I)
      YUC(I,J) =   COSU(I)  ;  YVD(I,J) = - SINP(J)*COSU(I)  ;  YWA(I,J) = COSP(J)*SINI(I)
      ZUC(I,J) = 0          ;  ZVD(I,J) =   COSP(J)          ;  ZWA(I,J) = SINP(J)
   30 Continue
      Do 40 J=J1,JNP  ;  Do 40 I=1,IMO
      XUD(I,J) = - SINI(I)  ;  XVC(I,J) = - SINV(J)*COSI(I)
      YUD(I,J) =   COSI(I)  ;  YVC(I,J) = - SINV(J)*COSI(I)
      ZUD(I,J) = 0          ;  ZVC(I,J) =   COSV(J)
   40 Continue

!**** Ocean straits on Earth's surface
      If (.not. AM_I_ROOT())  Return
      Do 50 N=1,NMST
      RLON1 = DLON*(IST(N,1)-.5   + .5*XST(N,1))
      RLON2 = DLON*(IST(N,2)-.5   + .5*XST(N,2))
      RLAT1 = DLAT*(JST(N,1)-FJEQ + .5*YST(N,1))
      RLAT2 = DLAT*(JST(N,2)-FJEQ + .5*YST(N,2))
      XSTN(N) = RADIUS*Cos(.5*(RLAT1+RLAT2))*Cos(.5*(RLON1+RLON2))
      YSTN(N) = RADIUS*Cos(.5*(RLAT1+RLAT2))*Sin(.5*(RLON1+RLON2))
      ZSTN(N) = RADIUS*Sin(.5*(RLAT1+RLAT2))
      dX = Cos(RLAT2)*Cos(RLON2) - Cos(RLAT1)*Cos(RLON1)
      dY = Cos(RLAT2)*Sin(RLON2) - Cos(RLAT1)*Sin(RLON1)
      dZ = Sin(RLAT2)            - Sin(RLAT1)
      D  = Sqrt (dX**2 + dY**2 + dZ**2)
      XUST(N) = dX / D  !  horizontal unit vector along the strait in
      YUST(N) = dY / D  !  (X,Y,Z) coordinates
   50 ZUST(N) = dZ / D
      Return
      EndSubroutine OTIDE0


      Subroutine OTIDEW (TIME)
!**** OTIDEW computes vertical perturbations of gravity due to Lunar and Solar tides.
      Use OTIDE
      Use CONSTANT,   Only: TWOPI,RADIUS
      Use RAD_COM,    Only: ECCEN=>ECCN,OBLIQD=>OBLIQ,OMEGVPD=>OMEGT
      Use OCEANR_DIM, Only: OGRID
      Use OCEANRES,   Only: IMO
      Use OCEAN,      Only: SINI=>SINIC,COSI=>COSIC,SINU,COSU, DLON,DLAT,FJEQ, &
                            SINP=>SINPO,COSP=>COSPO,SINV=>SINVO,COSV=>COSVO, LMOM=>LMM
      Use TimeConstants_MOD, Only: SECONDS_PER_DAY
      Implicit  None
      Real*8,Intent(In) :: TIME
!**** Local variables
      Integer :: I,J, J1,JN
      Real*8  :: DAY, &
                 MOOLON,MOOLAT,MOODIS, DMCmEC, &
                 SUNLON,SUNLAT,SUNDIS, DSCmEC, SIND,COSD,EQTIME, &
                 XGM,YGM,ZGM, &  !  gravity perturbations by Moon in (X,Y,Z) coordinates
                 XGS,YGS,ZGS, &  !  gravity perturbations by Sun  in (X,Y,Z) coordinates
                 zDMCmIJe3,   &  !  cube of reciprocal of distance from Moon center to cell (I,J)
                 zDSCmIJe3       !  cube of reciprocal of distance from Sun  center to cell (I,J)

!**** Extract domain decomposition
      J1 = OGRID%J_STRT  ;  JN = OGRID%J_STOP

!**** Locate Moon and Sun for given day, Model longitude starts from IDL
      DAY = TIME / SECONDS_PER_DAY
      Call MOON  (DAY, MOOLON,MOOLAT,MOODIS)
      Call ORBIT (OBLIQD,ECCEN,OMEGVPD, VE2000,EDAYzY,DAY, SUNDIS,SIND,COSD,SUNLON,SUNLAT,EQTIME)
      MOOLON = MOOLON + TWOPI*.5
      SUNLON = SUNLON + TWOPI*.5

!**** Compute distances between centers (m)
      DMCmEC = RADIUS*MOODIS  !  distance between Moon and Earth (m)
      DSCmEC = ASEMI *SUNDIS  !  distance between Sun and Earth (m)

!**** Compute location of Moon and Sun in (X,Y,Z) coordinates
      XMC = DMCmEC*Cos(MOOLAT)*Cos(MOOLON)
      YMC = DMCmEC*Cos(MOOLAT)*Sin(MOOLON)
      ZMC = DMCmEC*Sin(MOOLAT)
      XSC = DSCmEC*Cos(SUNLAT)*Cos(SUNLON)
      YSC = DSCmEC*Cos(SUNLAT)*Sin(SUNLON)
      ZSC = DSCmEC*Sin(SUNLAT)
      XMCzD3 = XMC / DMCmEC**3
      YMCzD3 = YMC / DMCmEC**3
      ZMCzD3 = ZMC / DMCmEC**3
      XSCzD3 = XSC / DSCmEC**3
      YSCzD3 = YSC / DSCmEC**3
      ZSCzD3 = ZSC / DSCmEC**3

!**** Loop over Model A-grid cells
      Do 10 J=J1,JN  ;  Do 10 I=1,IMO
      If (LMOM(I,J) <= 0)  GoTo 10
!**** Compute spatial perturbation of gravity due to Moon and Sun in (X,Y,Z) coordinates
      zDMCmIJe3 = ((XMC-XAIJ(I,J))**2 + (YMC-YAIJ(I,J))**2 + (ZMC-ZAIJ(I,J))**2) ** (-1.5)
      zDSCmIJe3 = ((XSC-XAIJ(I,J))**2 + (YSC-YAIJ(I,J))**2 + (ZSC-ZAIJ(I,J))**2) ** (-1.5)
      XGM = GUNIV*MMOON*((XMC-XAIJ(I,J))*zDMCmIJe3 - XMCzD3)
      YGM = GUNIV*MMOON*((YMC-YAIJ(I,J))*zDMCmIJe3 - YMCzD3)
      ZGM = GUNIV*MMOON*((ZMC-ZAIJ(I,J))*zDMCmIJe3 - ZMCzD3)
      XGS = GUNIV*MSUN *((XSC-XAIJ(I,J))*zDSCmIJe3 - XSCzD3) 
      YGS = GUNIV*MSUN *((YSC-YAIJ(I,J))*zDSCmIJe3 - YSCzD3)
      ZGS = GUNIV*MSUN *((ZSC-ZAIJ(I,J))*zDSCmIJe3 - ZSCzD3)
!**** Project (X,Y,Z) gravity perturbations onto W unit vectors on A-grid 
      WTIDE(I,J) = XWA(I,J)*(XGM+XGS) + YWA(I,J)*(YGM+YGS) + ZWA(I,J)*(ZGM+ZGS)
   10 Continue
      Return
      EndSubroutine OTIDEW


      Subroutine OTIDEV (QEVEN,DT1, UC,VC,UD,VD)
!**** OTIDEV accelerates horizontal ocean currents due to Lunar and Solar tides.
!**** Horizontal mass flux through straits is also accelerated.
      Use OTIDE
      Use OCEANR_DIM, Only: OGRID
      Use OCEANRES,   Only: IMO,JMO,LMO
      Use OCEAN,      Only: LMOU=>LMU,LMOV=>LMV
      Use STRAITS,    Only: NMST,LMST,MMST,MUST,DIST
      Implicit None
      Logical,Intent(In) :: QEVEN
      Real*8, Intent(In) :: DT1
      Real*8, Intent(InOut),Dimension(IMO,OGRID%J_STRT_HALO:OGRID%J_STOP_HALO,LMO) :: UC,VC,UD,VD
!**** Local variables
      Integer :: I,J,L,N, J1,JN,JNP
      Real*8  :: UACC,VACC,STACC, &  ! acceleration (m/s^2)
                 XGM,YGM,ZGM, &  !  gravity perturbations by Moon in (X,Y,Z) coordinates
                 XGS,YGS,ZGS, &  !  gravity perturbations by Sun  in (X,Y,Z) coordinates
                 zDMCmIJe3,   &  !  cube of reciprocal of distance from Moon center to cell (I,J)
                 zDSCmIJe3       !  cube of reciprocal of distance from Sun  center to cell (I,J)
 
!**** Extract domain decomposition
      J1 = OGRID%J_STRT  ;  JN = OGRID%J_STOP  ;  JNP = Min(JN,JMO-1)

!**** Loop over UC-grid locations
      Do 10 J=J1,JN  ;  Do 10 I=1,IMO
      If (LMOU(I,J) <= 0)  GoTo 10
!**** Compute spatial perturbation of gravity due to Moon and Sun in (X,Y,Z) coordinates
      zDMCmIJe3 = ((XMC-XUIJ(I,J))**2 + (YMC-YUIJ(I,J))**2 + (ZMC-ZUIJ(I,J))**2) ** (-1.5)
      zDSCmIJe3 = ((XSC-XUIJ(I,J))**2 + (YSC-YUIJ(I,J))**2 + (ZSC-ZUIJ(I,J))**2) ** (-1.5)
      XGM = GUNIV*MMOON*((XMC-XUIJ(I,J))*zDMCmIJe3 - XMCzD3)
      YGM = GUNIV*MMOON*((YMC-YUIJ(I,J))*zDMCmIJe3 - YMCzD3)
      ZGM = GUNIV*MMOON*((ZMC-ZUIJ(I,J))*zDMCmIJe3 - ZMCzD3)
      XGS = GUNIV*MSUN *((XSC-XUIJ(I,J))*zDSCmIJe3 - XSCzD3)
      YGS = GUNIV*MSUN *((YSC-YUIJ(I,J))*zDSCmIJe3 - YSCzD3)
      ZGS = GUNIV*MSUN *((ZSC-ZUIJ(I,J))*zDSCmIJe3 - ZSCzD3)
!**** Project (X,Y,Z) gravity pertubations onto U unit vectors on UC grid
!**** Project (X,Y,Z) gravity pertubations onto V unit vectors on VD grid = UC grid
      UACC = XUC(I,J)*(XGM+XGS) + YUC(I,J)*(YGM+YGS)  !  ZUC(I,J) = 0
      VACC = XVD(I,J)*(XGM+XGS) + YVD(I,J)*(YGM+YGS) + ZVD(I,J)*(ZGM+ZGS)
      UC(I,J,1:LMOU(I,J)) = UC(I,J,1:LMOU(I,J)) + UACC*DT1
      VD(I,J,1:LMOU(I,J)) = VD(I,J,1:LMOU(I,J)) + VACC*DT1
   10 Continue

!**** Loop over VC-grid locations
      Do 20 J=J1,JNP  ;  Do 20 I=1,IMO
      If (LMOV(I,J) <= 0)  GoTo 20
!**** Calculate spatial perturbation of gravity due to Moon and Sun in (X,Y,Z) coordinates
      zDMCmIJe3 = ((XMC-XVIJ(I,J))**2 + (YMC-YVIJ(I,J))**2 + (ZMC-ZVIJ(I,J))**2) ** (-1.5)
      zDSCmIJe3 = ((XSC-XVIJ(I,J))**2 + (YSC-YVIJ(I,J))**2 + (ZSC-ZVIJ(I,J))**2) ** (-1.5)
      XGM = GUNIV*MMOON*((XMC-XVIJ(I,J))*zDMCmIJe3 - XMCzD3)
      YGM = GUNIV*MMOON*((YMC-YVIJ(I,J))*zDMCmIJe3 - YMCzD3)
      ZGM = GUNIV*MMOON*((ZMC-ZVIJ(I,J))*zDMCmIJe3 - ZMCzD3)
      XGS = GUNIV*MSUN *((XSC-XVIJ(I,J))*zDSCmIJe3 - XSCzD3)
      YGS = GUNIV*MSUN *((YSC-YVIJ(I,J))*zDSCmIJe3 - YSCzD3)
      ZGS = GUNIV*MSUN *((ZSC-ZVIJ(I,J))*zDSCmIJe3 - ZSCzD3)
!**** Project (X,Y,Z) gravity pertubations onto V unit vectors on VC grid
!**** Project (X,Y,Z) gravity pertubations onto U unit vectors on UD grid = VC grid
      VACC = XVC(I,J)*(XGM+XGS) + YVC(I,J)*(YGM+YGS) + ZVC(I,J)*(ZGM+ZGS)
      UACC = XUD(I,J)*(XGM+XGS) + YUD(I,J)*(YGM+YGS)  !  ZUD(I,J) = 0
      VC(I,J,1:LMOV(I,J)) = VC(I,J,1:LMOV(I,J)) + VACC*DT1
      UD(I,J,1:LMOV(I,J)) = UD(I,J,1:LMOV(I,J)) + UACC*DT1
   20 Continue

      If (.not. QEVEN)  Return
!**** Loop over ocean straits, accelerate mass flux
      Do 30 N=1,NMST
!**** Calculate spatial perturbation of gravity due to Moon and Sun in (X,Y,Z) coordinates
      zDMCmIJe3 = ((XMC-XSTN(N))**2 + (YMC-YSTN(N))**2 + (ZMC-ZSTN(N))**2) ** (-1.5)
      zDSCmIJe3 = ((XSC-XSTN(N))**2 + (YSC-YSTN(N))**2 + (ZSC-ZSTN(N))**2) ** (-1.5)
      XGM = GUNIV*MMOON*((XMC-XSTN(N))*zDMCmIJe3 - XMCzD3)
      YGM = GUNIV*MMOON*((YMC-YSTN(N))*zDMCmIJe3 - YMCzD3)
      ZGM = GUNIV*MMOON*((ZMC-ZSTN(N))*zDMCmIJe3 - ZMCzD3)
      XGS = GUNIV*MSUN *((XSC-XSTN(N))*zDSCmIJe3 - XSCzD3)
      YGS = GUNIV*MSUN *((YSC-YSTN(N))*zDSCmIJe3 - YSCzD3)
      ZGS = GUNIV*MSUN *((ZSC-ZSTN(N))*zDSCmIJe3 - ZSCzD3)
!**** Project (X,Y,Z) gravity pertubations onto strait unit vactors
!**** coordinates onto strait vector of local topographic coordinates
      STACC = XUST(N)*(XGM+XGS) + YUST(N)*(YGM+YGS) + ZUST(N)*(ZGM+ZGS)
      Do 30 L=1,LMST(N)
   30 MUST(L,N) = MUST(L,N) + STACC*DT1*MMST(L,N)/DIST(N)
      Return
      EndSubroutine OTIDEV


      Subroutine MOON (DAY, MOOLON,MOOLAT,MOODIS)
!****
!**** Input: DAY = days measured since 2000 January 1, hour 0
!****
!**** Output: MOOLON = Earth's LONgitude directly beneath the MOOn
!****         MOOLAT = Earth's LATitude directly beneath the MOOn
!****         MOODIS = DIStance from MOOn to Earth in Earth radii units
!****
!**** Constants: OBLIQ = dihedral angle between Earth's equatorial
!****                    plane and orbital plane of year 2000
!****           TAofVE = true anomaly of vernal equinox of year 2000
!****           MAofVE = mean anomaly of vernal equinox of year 2000
!****           EDAYzY = tropical year = Earth days per year = 365.2425
!****           VE200D = days from 2000 January 1, hour 0 till vernal
!****                    equinox of year 2000 = 31 + 29 + 19 + 7.5/24
!****
!**** Intermediate quantities:
!****       RA = right ascension of the Moon =
!****          = longitude of Moon in geocentric equatorial coordinates
!****   VEQLON = longitude of Greenwich Meridion in geocentric
!****            equatorial coordinates at vernal equinox of year 2000
!****   ROTATE = change in longitude in geocentric nonrotating
!****            equatorial coordinates from a point's location on
!****            vernal equinox to is current location where the point
!****            is fixed on rotating Earth
!****
!**** Angles, longitude and latitude are measured in radians unless
!**** the variable is terminated by a capital D in which case it is
!**** measured in degrees
!****
      Use CONSTANT, Only: TWOPI
      Implicit None
      Real*8,Intent(In)  :: DAY
      Real*8,Intent(Out) :: MOOLON,MOOLAT,MOODIS
!**** Local variables
      Real*8,Parameter :: PIz180=TWOPI/360, OBLIQ=.409101d0, TAofVE=1.345728d0, MAofVE=1.313255d0, &
                          EDAYzY=365.2425d0, VE200D=79.3125
      Integer :: I, RADEG,RAMIN, MLDEG,MLMIN
      Real*8  :: TJC, ECLOND,ECLATD,HOPARD, ECLON,ECLAT,HOPAR, &
                 XEC,YEC,ZEC, XEQ,YEQ,ZEQ, RA,VEQLON,ROTATE, &
                 ECLONN(3,0:6),ECLATN(3,4),HOPARN(3,0:4)  !  in degrees
      Data ECLONN /218.32d0,   0.   ,  481267.883d0,&
                     6.29d0, 134.9d0,  477198.85d0, &
                    -1.27d0, 259.2d0, -413335.38d0, &
                      .66d0, 235.7d0,  890534.23d0, &
                      .21d0, 269.9d0,  954397.70d0, &
                     -.19d0, 357.5d0,   35999.05d0, &
                     -.11d0, 186.6d0,  966404.05d0/
      Data ECLATN /  5.13d0,  93.3d0,  483202.03d0, &
                      .28d0, 228.2d0,  960400.87d0, &
                     -.28d0, 318.3d0,    6003.18d0, &
                     -.17d0, 217.6d0, -407332.20d0/
      Data HOPARN / .9508d0,   0.   ,       0.    , &
                    .0518d0, 134.9d0,  477198.85d0, &
                    .0095d0, 259.2d0, -413335.38d0, &
                    .0078d0, 235.7d0,  890534.23d0, &
                    .0028d0, 269.9d0,  954397.70d0/
!**** Calculate time in Julian centuries
      TJC = (DAY-.5) / (100*365.25)
!****
!**** Geocentric Ecliptic coordinates:
!**** EcLat = angle between an object and its perpendicular projection
!****         onto the Ecliptic plane with Earth center as angle vertex
!**** EcLon = angle from vernal equinox vector (ray from Earth to Sun
!****         at vernal equinox) to projection onto the ecliptic plane
!****         of vector of current Earth to object
!****
!**** Calculate ecliptic longitude ECLOND in degrees and radians ECLON
      ECLOND = ECLONN(1,0) + ECLONN(3,0)*TJC
      Do 10 I=1,6
   10 ECLOND = ECLOND + ECLONN(1,I)*Sin(PIz180*(ECLONN(2,I)+ECLONN(3,I)*TJC))
      ECLOND = Mod (ECLOND, 360d0)
      ECLON  = ECLOND*PIz180
!**** Calculate ecliptic latitude ECLATD in degrees and radians ECLAT
      ECLATD = 0
      Do 20 I=1,4
   20 ECLATD = ECLATD + ECLATN(1,I)*Sin(PIz180*(ECLATN(2,I)+ECLATN(3,I)*TJC))
      ECLAT  = ECLATD*PIz180
!**** Calculate horizontal parallax in degrees HOPARD and radians HOPAR
      HOPARD = HOPARN(1,0)
      Do 30 I=1,4
   30 HOPARD = HOPARD + HOPARN(1,I)*Cos(PIz180*(HOPARN(2,I)+HOPARN(3,I)*TJC))
      HOPAR  = HOPARD*PIz180
!**** Calculate MOODIS in Earth radii
      MOODIS = 1 / Sin(HOPAR)
!**** Calculate geocentric ecliptic rectangular coordinates
      XEC = Cos(ECLAT)*Cos(ECLON)
      YEC = Cos(ECLAT)*Sin(ECLON)
      ZEC = Sin(ECLAT)
!****
!**** Geocentric nonrotating Equatorial coordinates:
!**** Declination = angle between an object and its perpendicular
!****               projection onto the Equatorial plane with Earth
!****               center as angle vertex
!**** RightAscention = angle from the vernal equinox vector (ray from
!****                  Earth to Sun at vernal equinox) to projection
!****                  onto equatorial plane of vector of current Earth
!****                  to object
!****
!**** Calculate geocentric equatorial rectangular coordinates
      XEQ = XEC
      YEQ = YEC*Cos(OBLIQ) - ZEC*Sin(OBLIQ)
      ZEQ = ZEC*Cos(OBLIQ) + YEC*Sin(OBLIQ)
!**** Calculate right ascension and MOOLON
      RA     = ATan2 (YEQ,XEQ)
      VEQLON = TWOPI*VE200D - TWOPI/2 + MAofVE - TAofVE  !  modulo 2*Pi
      ROTATE = TWOPI*(DAY-VE200D)*(EDAYzY+1)/EDAYzY
      MOOLON = Modulo (RA-ROTATE-VEQLON, TWOPI)  ;  If (MOOLON > TWOPI/2) MOOLON = MOOLON - TWOPI
!**** Calculate declination = MOOLAT in radians
      MOOLAT = ASin (ZEQ)
      Return
      EndSubroutine MOON
