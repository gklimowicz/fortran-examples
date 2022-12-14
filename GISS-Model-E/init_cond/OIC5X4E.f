C**** OIC5X4E.CRE    CREate Ocean Initial Conditions    2008/03/04
C****
C**** Fortran Compile and Go:  FCG OIC5X4E.CRE
C****
C**** To change to different resolution:
C**** Set: IM,JM,LMO,LMO_MIN (both in MAIN and in subroutine FIXED)
C**** Set: ZOE(0:LMO) =      (in MAIN)
C**** Set: FILEI? and FILOUT (in MAIN)
C**** Set: "Fill in missing WOA98 data by hand"  (before  Do 490 J=1,JM)
C****
C**** Calculate WOA98 temperature and salinity data from monthly MC
C**** DataFiles to a December 1, hour 0 data.
C**** Interpolate December 1, hour 0 data from 1x1 to 2Hx2 resolution.
C**** Convert WOA98 ocean data from temperature to potential specific
C**** enthalpy.
C**** Calculate ocean bottom topography for present resolution.
C**** Fill in missing data for G and S and check all cells are filled.
C**** Integrated G and S to model layers; calculate vertical gradients.
C****
C**** Input:
C**** $OBS/WOA98/TO1X1.MC = Ocean Temperature (C) at 1x1 and 33 levels
C**** $OBS/WOA98/SO1X1.MC = Ocean Salinity (psu) at 1x1 and 33 levels
C**** $IFDIR/Z72X46N_gas.1_nocasp_21k_ice5g = 5x4 Z file
C**** $IFDIR/OFTABLE_NEW            = TABLE of Ocean Functions
C**** $IFDIR/GIC.E046DM20A.1DEC1955 = Sea-Ice Cover (%) & Mass (kg/m^2)
C****
C**** Output: OIC5X4_21k_ice5g = Ocean IC for December 1
C****
      Implicit  Real*8 (A-H,M,O-Z)
      Integer*4,Parameter :: IM=72,JM=46,LMO=13,LMO_MIN=2  !  modelE
      Real*8,   Parameter :: GRAV=9.80665d0, RHOW=1000, RHOO=1035,
     *                       zRT3=.577350269d0, zRT12=.288675135d0,
C    *  ZOE(0:LMO) = (/ 0,  12,  30,  56,  92,  140, 202, 280, 376,
C    *                     492, 626, 776, 940, 1116,1302,1496,1696,
C    *                    1900,2106,2312,2518, 2724,2930,3136,3342,
C    *                    3548,3754,3960,4166, 4372,4578,4784,4990 /),
     *  dZO1  = 12,
     *  ZORAT = 1.5,
     *  ZOE(0:LMO) = (/ (dZO1*(ZORAT**L-1)/(ZORAT-1),L=0,LMO) /),
     *  dZO(LMO) = (/ (ZOE(L)-ZOE(L-1),L=1,LMO) /)
      Integer*4 LMOM
      Real*8    MFO(LMO,LMO), MO(IM,JM,LMO)
      Common /FIXDCB/ FOCEAN(IM,JM),ZOCEAN(IM,JM), LMOM(IM,JM)
C****
      Integer*4,Parameter :: IM1=360,JM1=180,KM1=33  !  WOA98 resolution
      Real*8,   Parameter :: DATMIS = -999999,  !  MISsing DATa value
     *  Z1(KM1) = (/ 0,  10,  20,  30,   50,  75, 100, 125,
     *             150, 200, 250, 300,  400, 500, 600, 700,
     *             800, 900,1000,1100, 1200,1300,1400,1500,
     *            1750,2000,2500,3000, 3500,4000,4500,5000, 5500 /)
C**** Variables WOA98 horizontal resolution and levels, Oct to Jan
      Real*4 TWOA(IM1,JM1,KM1,11:12), !  tempreature (C)
     *       SWOA(IM1,JM1,KM1,11:12)  !  salinity (psu)
C**** Variables Dec 1, WOA98 horizontal resolution and levels
      Real*8 WTA(IM1,JM1),     !  weight array for interpolation
     *    TWOAD1(IM1,JM1,KM1), !  temperature (C)
     *    SWOAD1(IM1,JM1,KM1), !  salinity (psu)
C**** Variables Dec 1, Model horizontal resolution, WOA98 levels
     *    TD1(IM,JM,KM1),  !  temperature (C)
     *    SD1(IM,JM,KM1),  !  salinity (psu)
     *    PD1(IM,JM,KM1),  !  pressure (Pa)
     *    VD1(IM,JM,KM1),  !  specific volume (m^3/kg)
     *    GD1(IM,JM,KM1)   !  potential specific enthalpy (J/kg)
      Integer*4 KM(IM,JM)  !  number of WOA98 levels for each I,J
C**** Variables for single vertical column, WOA98 levels
      Real*8     GK(KM1),  !  potential specific enthalpy (J/kg)
     *           SK(KM1),  !  salinity (psu)
C**** Variables for Model layers, Model horizontal resolution
     *    G0(LMO,IM,JM),   !  mean potential specific enthalpy (J/kg)
     *    GZ(LMO,IM,JM),   !  vertical gradient of GL (J/kg)
     *    S0(LMO,IM,JM),   !  mean salinity (psu)
     *    SZ(LMO,IM,JM)    !  vertical gradient of salinity (psu)
C**** Two-dimensional variables on Model horizontal resolution
      Real*8  RSI(IM,JM),  !  sea ice cover (1)
     *      HSI(4,IM,JM),  !  sea ice heat content (J/m^2)
     *       SNOW(IM,JM),  !  snow mass over sea ice (kg/m^2)
     *       MSI2(IM,JM),  !  sea ice mass of layer 2 (kg/m^2)
     *        MSI(IM,JM),  !  sea ice mass of whole cell (kg/m^2)
     *       ZLOD1(IM,JM), !  liquid ocean height (m) on December 1
     *       PSLD1(IM,JM), !  sea level pressure (Pa) on December 1
C**** Ocean functions of Pressure (Pa), Temperature (C), and Salinity
     *       PHPTS,        !  potential specific enthalpy (J/kg)
     *      VOLPTS,VOLGSP  !  specific vloume (m^3/kg)
C****
      Character*80 TITLE, OBS,IFDIR
      Character*96 FILEIT,FILEIS,FILEIZ,FILEIO,FILEII, FILOUT
      Integer*4 I,J,K,L,M, KMG,KMS, LM,
     *       ITMIN,JTMIN,KTMIN, ISMIN,JSMIN,KSMIN, IGMIN,JGMIN,KGMIN,
     *       ITMAX,JTMAX,KTMAX, ISMAX,JSMAX,KSMAX, IGMAX,JGMAX,KGMAX
      Real*8 TMIN,TMAX, SMIN,SMAX, GMIN,GMAX,
     *       GZOSI,GAREA, PE,ZE, PUP,PDN, VUP,VDN, dMCOL
      integer iargc
C****
c      Call GetEnv ('OBS'  ,OBS)
c      Call GetEnv ('IFDIR',IFDIR)
      OBS    = '/discover/nobackup/projects/giss/OBS'
      IFDIR  = '/discover/nobackup/projects/giss/prod_input_files'
      FILEIT = Trim(OBS)   // '/WOA98/TO1X1.MC'
      FILEIS = Trim(OBS)   // '/WOA98/SO1X1.MC'
      FILEIZ = Trim(IFDIR) // '/Z72X46N_gas.1_nocasp_21k_ice5g'
      FILEIO = Trim(IFDIR) // '/OFTABLE_NEW'
      FILEII = Trim(IFDIR) // '/GIC.E046D3M20A.1DEC1955'
      FILOUT = Trim(IFDIR) // '/OIC5X4_21k_ice5g'
C****
C**** Read WOA98 temperature and salinity data MC DataFiles
C****
C**** Temperature (C)
      Open (1, File=FILEIT, Form='Unformatted', Status='Old')
      Do 10 M=1,10   !  Skip data for January through October
      Do 10 K=1,KM1
   10 Read (1)
      Do 20 M=11,12  !  Read data for November and December
      Do 20 K=1,KM1
   20 Read (1) TITLE,TWOA(:,:,K,M)
      Close(1)
C**** Salinity (psu)
      Open (1, File=FILEIS, Form='Unformatted', Status='Old')
      Do 30 M=1,10   !  Skip data for January through October
      Do 30 K=1,KM1
   30 Read (1)
      Do 40 M=11,12  !  Read data for November and December
      Do 40 K=1,KM1
   40 Read (1) TITLE,SWOA(:,:,K,M)
      Close(1)
      Write (0,*) 'WOA98 data read in'
C****
C**** Calculate WOA98 temperature and salinity data from monthly MC
C**** DataFiles to a December 1, hour 0 data (average Nov and Dec data)
C****
      TMIN = 1d6  ;  TMAX = -1d6  ;  SMIN = 1d6  ;  SMAX = -1d6
      Do 120 K=1,KM1
      Do 120 J=1,JM1
      Do 120 I=1,IM1
C**** Temperature (C)
      TWOAD1(I,J,K) = DATMIS
      If (TWOA(I,J,K,11)==DATMIS .or. TWOA(I,J,K,12)==DATMIS)  GoTo 110
      TWOAD1(I,J,K) = (TWOA(I,J,K,11) + TWOA(I,J,K,12)) / 2
      If (TMIN > TWOAD1(I,J,K))  Then
          TMIN = TWOAD1(I,J,K) ; ITMIN=I ; JTMIN=J ; KTMIN=K ; EndIf
      If (TMAX < TWOAD1(I,J,K))  Then
          TMAX = TWOAD1(I,J,K) ; ITMAX=I ; JTMAX=J ; KTMAX=K ; EndIf
C**** Salinity (psu)
  110 SWOAD1(I,J,K) = DATMIS
      If (SWOA(I,J,K,11)==DATMIS .or. SWOA(I,J,K,12)==DATMIS)  GoTo 120
      SWOAD1(I,J,K) = (SWOA(I,J,K,11) + SWOA(I,J,K,12)) / 2
      If (SMIN > SWOAD1(I,J,K))  Then
          SMIN = SWOAD1(I,J,K) ; ISMIN=I ; JSMIN=J ; KSMIN=K ; EndIf
      If (SMAX < SWOAD1(I,J,K))  Then
          SMAX = SWOAD1(I,J,K) ; ISMAX=I ; JSMAX=J ; KSMAX=K ; EndIf
  120 Continue
      Write (0,*) 'I,J,K,TMIN=',ITMIN,JTMIN,KTMIN,TMIN
      Write (0,*) 'I,J,K,TMAX=',ITMAX,JTMAX,KTMAX,TMAX
      Write (0,*) 'I,J,K,SMIN=',ISMIN,JSMIN,KSMIN,SMIN
      Write (0,*) 'I,J,K,SMAX=',ISMAX,JSMAX,KSMAX,SMAX
C****
C**** Interpolate December 1, hour 0 data from 1x1 to 5x4 resolution
C****
      Call HNTR80 (IM1,JM1,0d0,60d0, IM,JM,0d0,240d0, DATMIS)
      Do 230 K=1,KM1                        !  240 minutes = 4 degrees
C**** Temperature (C)
      Do 210 J=1,JM1
      Do 210 I=1,IM1
      WTA(I,J) = 1
  210 If (TWOAD1(I,J,K) == DATMIS)  WTA(I,J) = 0
      Call HNTR8P (WTA, TWOAD1(1,1,K), TD1(1,1,K))
C**** Salinity (psu)
      Do 220 J=1,JM1
      Do 220 I=1,IM1
      WTA(I,J) = 1
  220 If (SWOAD1(I,J,K) == DATMIS)  WTA(I,J) = 0
  230 Call HNTR8P (WTA, SWOAD1(1,1,K), SD1(1,1,K))
C****
C**** Convert ocean temperature to potential specific enthalpy
C****
      GMIN = 1d6  ;  GMAX = -1d6
      Do 330 J=1,JM
      Do 330 I=1,IM
C**** Calculate potential specific enthalpy for the first level
      K = 1
      If (TD1(I,J,K) <= DATMIS .or. SD1(I,J,K) <= DATMIS)  GoTo 320
      PD1(I,J,K) = 0
      VD1(I,J,K) = VOLPTS (PD1(I,J,K), TD1(I,J,K), SD1(I,J,K)*1d-3)
      GD1(I,J,K) =  PHPTS (PD1(I,J,K), TD1(I,J,K), SD1(I,J,K)*1d-3)
C**** Calculate potential specific enthalpy for subsequent levels
C**** Use 3 iterations of PD1(I,J,K) and VD1(I,J,K) to converge
  310 K = K + 1  ;  If (K > KM1)  GoTo 330
      If (TD1(I,J,K) <= DATMIS .or. SD1(I,J,K) <= DATMIS)  GoTo 320
      PD1(I,J,K) = PD1(I,J,K-1) + GRAV*(Z1(K)-Z1(K-1)) / VD1(I,J,K-1)
      VD1(I,J,K) = VOLPTS (PD1(I,J,K), TD1(I,J,K), SD1(I,J,K)*1d-3)
      PD1(I,J,K) = PD1(I,J,K-1) +
     +             GRAV*(Z1(K)-Z1(K-1)) * 2/(VD1(I,J,K-1)+VD1(I,J,K))
      VD1(I,J,K) = VOLPTS (PD1(I,J,K), TD1(I,J,K), SD1(I,J,K)*1d-3)
      PD1(I,J,K) = PD1(I,J,K-1) +
     +             GRAV*(Z1(K)-Z1(K-1)) * 2/(VD1(I,J,K-1)+VD1(I,J,K))
      VD1(I,J,K) = VOLPTS (PD1(I,J,K), TD1(I,J,K), SD1(I,J,K)*1d-3)
      GD1(I,J,K) =  PHPTS (PD1(I,J,K), TD1(I,J,K), SD1(I,J,K)*1d-3)
      If (GMIN > GD1(I,J,K))  Then
          GMIN = GD1(I,J,K) ; IGMIN=I ; JGMIN=J ; KGMIN=K ; EndIf
      If (GMAX < GD1(I,J,K))  Then
          GMAX = GD1(I,J,K) ; IGMAX=I ; JGMAX=J ; KGMAX=K ; EndIf
      GoTo 310
C**** Missing data is encountered for level K
  320 GD1(I,J,K:KM1) = DATMIS
  330 Continue
      Write (0,*) 'I,J,K,GMIN=',IGMIN,JGMIN,KGMIN,GMIN
      Write (0,*) 'I,J,K,GMAX=',IGMAX,JGMAX,KGMAX,GMAX
C****
C**** Calculate ocean bottom topography for present resolution.
C**** Read in ocean table for VGSP (specific volume)
C****
      Open (11, File=FILEIZ, Form='Unformatted', Status='Old')
      Open (22, File=FILEIO, Form='Unformatted', Status='Old')
      Call FIXED (ZOE,dZO, MFO)
C****
C**** Fill in missing data for G and S and check all cells are filled
C****
C**** Determine KM, the number of valid WOA98 levels for each I,J cell
      Do 420 J=1,JM
      Do 420 I=1,IM
      Do 410 K=1,KM1
  410 If (GD1(I,J,K) <= DATMIS .or. SD1(I,J,K) <= DATMIS)  GoTo 420
C     K = KM1+1
  420 KM(I,J) = K-1
C**** Check whether each ocean cell has enough WOA98 levels to
C**** determine all model layers.  Copy data from adjacent longitudes.
      Do 450 J=2,JM-1
      Do 450 I=1,IM
      If (ZOE(LMOM(I,J)) <= Z1(KM(I,J)))  GoTo 450
      Is1 = I-1  ;  If(I==1)  Is1 = IM
      Ia1 = I+1  ;  If(I==IM) Ia1 = 1
      KMC = KM(I,J)  ;  KMW = KM(Is1,J)  ;  KME = KM(Ia1,J)
      If (KMC >= Max(KMW,KME))  Goto 450
      If (KMC <  Min(KMW,KME))  Then
        KMMIN =  Min(KMW,KME)
        GD1(I,J,KMC+1:KMMIN) = .5*(GD1(Is1,J,KMC+1:KMMIN) +
     *                             GD1(Ia1,J,KMC+1:KMMIN))
        SD1(I,J,KMC+1:KMMIN) = .5*(SD1(Is1,J,KMC+1:KMMIN) +
     *                             SD1(Ia1,J,KMC+1:KMMIN))
        Write (6,943) 'West and East',I,J,KMC+1,KMMIN
        KM(I,J) = KMMIN  ;  KMC = KMMIN  ;  EndIf
      If (KME-KMW) 430,450,440
  430 GD1(I,J,KMC+1:KMW) = GD1(Is1,J,KMC+1:KMW)
      SD1(I,J,KMC+1:KMW) = SD1(Is1,J,KMC+1:KMW)
      Write (6,943) 'West',I,J,KMC+1,KMW
      KM(I,J) = KMW
      GoTo 450
  440 GD1(I,J,KMC+1:KME) = GD1(Ia1,J,KMC+1:KME)
      SD1(I,J,KMC+1:KME) = SD1(Ia1,J,KMC+1:KME)
      Write (6,943) 'East',I,J,KMC+1,KME
      KM(I,J) = KME
  450 Continue
C**** Copy data from adjacent latitudes if necessary
      Do 480 J=2,JM-1
      Do 480 I=1,IM
      If (ZOE(LMOM(I,J)) <= Z1(KM(I,J)))  GoTo 480
      KMC = KM(I,J)  ;  KMS = KM(I,J-1)  ;  KMN = KM(I,J+1)
      If (KMC >= Max(KMS,KMN))  GoTo 480
      If (KMC <  Min(KMS,KMN))  Then
        KMMIN =  Min(KMS,KMN)
        GD1(I,J,KMC+1:KMMIN) = .5*(GD1(I,J-1,KMC+1:KMMIN) +
     *                             GD1(I,J+1,KMC+1:KMMIN))
        SD1(I,J,KMC+1:KMMIN) = .5*(SD1(I,J-1,KMC+1:KMMIN) +
     *                             SD1(I,J+1,KMC+1:KMMIN))
        Write (6,943) 'South and North',I,J,KMC+1,KMMIN
        KM(I,J) = KMMIN  ;  KMC = KMMIN  ;  EndIf
      If (KMN-KMS) 460,480,470
  460 GD1(I,J,KMC+1:KMS) = GD1(I,J-1,KMC+1:KMS)
      SD1(I,J,KMC+1:KMS) = SD1(I,J-1,KMC+1:KMS)
      Write (6,943) 'South',I,J,KMC+1,KMS
      KM(I,J) = KMS
      GoTo 480
  470 GD1(I,J,KMC+1:KMN) = GD1(I,J+1,KMC+1:KMN)
      SD1(I,J,KMC+1:KMN) = SD1(I,J+1,KMC+1:KMN)
      Write (6,943) 'North',I,J,KMC+1,KMN
      KM(I,J) = KMN
  480 Continue
C**** Fill in missing WOA98 data by hand
C     KM(119,16) = KM(122,16)  ;  GD1(119,16,32 ) = GD1(122,16,32 )
C                                 SD1(119,16,32 ) = SD1(122,16,32 )
C     KM(120,16) = KM(122,16)  ;  GD1(120,16,32 ) = GD1(122,16,32 )
C                                 SD1(120,16,32 ) = SD1(122,16,32 )
C     KM( 75,65) = KM( 78,65)  ;  GD1( 75,65,28:) = GD1( 78,65,28:)
C                                 SD1( 75,65,28:) = SD1( 78,65,28:)
C     KM(102,80) = KM(102,81)  ;  GD1(102,80,  :) = GD1(101,82,  :)
C                                 SD1(102,80,  :) = SD1(101,82,  :)
C     KM( 38,83) = KM( 40,83)  ;  GD1( 38,83,12:) = GD1( 40,83,12:)
C                                 SD1( 38,83,12:) = SD1( 40,83,12:)
C     KM( 39,83) = KM( 40,83)  ;  GD1( 39,83,12:) = GD1( 40,83,12:)
C                                 SD1( 39,83,12:) = SD1( 40,83,12:)
C**** Check whether all necessary data has been filled
      KERROR = 0
      Do 490 J=1,JM
      Do 490 I=1,IM
      LM = LMOM(I,J)
      If (LM == 0)  GoTo 490
      If (Z1(KM(I,J)) > ZOE(LM-1))  GoTo 490
      Write (0,949) I,J,KM(I,J),LM,Z1(KM(I,J)),ZOE(LM-1)
      Write (6,949) I,J,KM(I,J),LM,Z1(KM(I,J)),ZOE(LM-1)
      KERROR = KERROR + 1
  490 Continue
      If (KERROR == 0)  GoTo 500
      Write (0,*) 'Fill WOA98 data by hand above  Do 490 J=1,JM'
      Write (6,*) 'Fill WOA98 data by hand above  Do 490 J=1,JM'
      Stop 490
C****
C**** Integrate G and S to model layers; calculate vertical gradients
C****
  500 G0(:,:,:) = DATMIS  ;  GZ(:,:,:) = DATMIS
      S0(:,:,:) = DATMIS  ;  SZ(:,:,:) = DATMIS
      Do 510 J=1,JM
      Do 510 I=1,IM
      If (FOCEAN(I,J) == 0)  GoTo 510
      KMIJ = KM(I,J)
      GK(1:KMIJ) = GD1(I,J,1:KMIJ)
      SK(1:KMIJ) = SD1(I,J,1:KMIJ)
      Call VLKtoLZ (KMIJ,LMOM(I,J), Z1,ZOE, GK, G0(1,I,J),GZ(1,I,J))
      Call VLKtoLZ (KMIJ,LMOM(I,J), Z1,ZOE, SK, S0(1,I,J),SZ(1,I,J))
  510 Continue
      Write (0,*) 'Vertical interpolation completed'
C****
C**** Calculate MSI from GIC.144X90.DEC01.1
C**** Read in ZOSI and PSL from prior Model simulation
C****
C**** Calculate MSI
      Open (1, File=FILEII, Form='Unformatted', Status='Old')
      Read (1)
      Read (1) TITLE,RSI,HSI,SNOW,MSI2
      Close(1)
C     MSI(:,:) = .01d0*RSI(:,:) * (MSI1(:,:) + MSI2(:,:))
      MSI(:,:) =       RSI(:,:) * (SNOW(:,:) + MSI2(:,:) + .1d0*916.6d0)
      ZLOD1(:,:) =            - MSI(:,:)/RHOW
      PSLD1(:,:) = 101325
      Write (0,*) 'Sea ice data read in'
C****
C**** Iteratively solve for MO so that integrated Z matches ZOSI
C****
      MO(:,:,:) = 0
      Do 730 J=1,JM
      Do 730 I=1,IM
      If (FOCEAN(I,J) == 0)  GoTo 730
      LM = LMOM(I,J)
      MO(I,J,1:LM) = RHOO*dZO(1:LM)  !  initial guess
C**** Add heights from each layer from ZSOLID
  710 PE = PSLD1(I,J) - 101325
      ZE = - ZOCEAN(I,J)
      DO 720 L=1,LM
      PUP = PE + MO(I,J,L)*GRAV*(.5-zRT12)
      PDN = PE + MO(I,J,L)*GRAV*(.5+zRT12)
      VUP = VOLGSP (G0(L,I,J)-zRT3*GZ(L,I,J),
     *             (S0(L,I,J)-zRT3*SZ(L,I,J))*1d-3, PUP)
      VDN = VOLGSP (G0(L,I,J)+zRT3*GZ(L,I,J),
     *             (S0(L,I,J)+zRT3*SZ(L,I,J))*1d-3, PDN)
      PE  = PE + MO(I,J,L)*GRAV
  720 ZE  = ZE + MO(I,J,L)*(VUP+VDN)*.5
C**** Modify MO and go back for another iteration
      dMCOL = RHOO*(ZE-ZLOD1(I,J))  !  excess mass in column
      MO(I,J,1:LM) = MO(I,J,1:LM) - dMCOL*MFO(1:LM,LM)
      If (Abs(ZE-ZLOD1(I,J)) > 1d-6)  GoTo 710
  730 Continue
      Write (0,*) 'Ocean mass calculated'
C****
C**** Write OIC output DataFile
C****
      Open (20, File=FILOUT, Form='Unformatted')
      TITLE = 'OIC 5x4 L13 21k ice5g from WOA98 12/01  ' //
     *        'TITLE,MO,G0,GZ,S0,SZ'
      Write (20) TITLE,
     *  (((Sngl(MO(I,J,L)),I=1,IM),J=1,JM),L=1,LMO),
     *  (((Sngl(G0(L,I,J)),I=1,IM),J=1,JM),L=1,LMO),
     *  (((Sngl(GZ(L,I,J)),I=1,IM),J=1,JM),L=1,LMO),
     *  (((Sngl(S0(L,I,J)*1d-3),I=1,IM),J=1,JM),L=1,LMO),
     *  (((Sngl(SZ(L,I,J)*1d-3),I=1,IM),J=1,JM),L=1,LMO)
C**** Liquid Ocean Mass (kg/m^2)
C     TITLE = 'LIQUID OCEAN MASS (kg/m^2) of Layer  1 at    6 m  ' //
C    *        'WOA98    1900:1997/12/01      '
C     Do 810 L=1,LMO
C     Write (TITLE(37:38),981) L
C     Write (TITLE(43:46),982) Nint(.5*(ZOE(L)+ZOE(L-1)))
C     Write (20) TITLE,Sngl(MO(:,:,L))
C 810 Write (0,983) Trim(TITLE),MO(1,30,L)
C**** Potential Specific Enthalpy (J/kg)
C     TITLE(1:50) = 'POTENTIAL SPECIFIC ENTHALPY (J/kg) at    6 m'
C     Do 820 L=1,LMO
C     Write (TITLE(39:42),982) Nint(.5*(ZOE(L)+ZOE(L-1)))
C     Write (20) TITLE,Sngl(G0(L,:,:))
C 820 Write (0,983) Trim(TITLE),G0(L,1,30)
C**** Vertical Gradient of Potential Specific Enthalpy (J/kg)
C     TITLE(1:50) = 'VERT GRAD of POT SPEC ENTHALPY (J/kg) at    6 m'
C     Do 830 L=1,LMO
C     Write (TITLE(42:45),982) Nint(.5*(ZOE(L)+ZOE(L-1)))
C     Write (20) TITLE,Sngl(GZ(L,:,:))
C 830 Write (0,983) Trim(TITLE),GZ(L,1,30)
C**** Salinity (psu)
C     TITLE(1:50) = 'SALINITY (psu) at    6 m'
C     Do 840 L=1,LMO
C     Write (TITLE(19:22),982) Nint(.5*(ZOE(L)+ZOE(L-1)))
C     Write (20) TITLE,Sngl(S0(L,:,:))
C 840 Write (0,983) Trim(TITLE),S0(L,1,30)
C**** Vertical Gradient of Salinity (psu)
C     TITLE(1:50) = 'VERTICAL GRADIENT of SALINITY (psu) at    6 m'
C     Do 850 L=1,LMO
C     Write (TITLE(40:43),982) Nint(.5*(ZOE(L)+ZOE(L-1)))
C     Write (20) TITLE,Sngl(SZ(L,:,:))
C 850 Write (0,983) Trim(TITLE),SZ(L,1,30)
      Close (20)

      do i=1,im
        do j=1,jm
          do l=1,lmo
            if((mo(i,j,l).gt.0).and.(s0(l,i,j).lt.0) ) print*
     *           ,"ijl s0 g0 mo",i,j,l,s0(l,i,j),g0(l,i,j),mo(i,j,l)
          enddo
        enddo
      enddo

      GoTo 999
C****
  943 Format (' WOA98 data is filled from ',A,'.  I,J,KMC+1:KMNEW =',
     *        4I4)
  949 Format (' WOA98 data is missing.  I,J,KM,LM,Z1(KM),ZOE(LM-1) =',
     *        4I4,2F8.1)
  981 Format (I2)
  982 Format (I4)
  983 Format (1X,A,F20.10)
      

 999  End

C**** C540C.FOR   Atmosphere-Ocean Model Common Blocks    2007/12/13
C****
      Subroutine FIXED (ZOE,dZO, MFO)
C****
C**** Read in fixed arrays and calculate arrays derived from them
C****
C     Use C540, Only: IM,JM,LMA,LMO, ZOE,dZO,MFO,
C    *    FOCEAN,FLAKE,FGRND,FGICE,FWATER,FLAND, ZATMO,
C    *    ZOCEAN,ZLAKE,ZSOLID,zFLAKE, LMAM,LMAU,LMAV, LMOM,LMOU,LMOV,
C    *    VGSP,TGSP,CGS,HGSP
      Implicit  Real*8 (A-H,M,O-Z)
      Integer*4,Parameter :: IM=72,JM=46,LMO=13,LMO_MIN=2
      Real*8    ZOE(0:LMO),dZO(LMO), MFO(LMO,LMO)
      Character TITLE*80
      Common /FIXDCB/ FOCEAN(IM,JM),ZOCEAN(IM,JM), LMOM(IM,JM)
      Common /OFUNCB/ VGSP(-2:40,0:40,0:39)
C**** Read in FIXDCB: FOCEAN, FLAKE, FGRND, FGICE, ZATMO, ZOCEAN, ZLAKE
      Call READR4 (11,IM*JM,FOCEAN,FOCEAN)
      Read (11)  !  skip FLAKE
      Read (11)  !  skip FGRND
      Read (11)  !  skip FGICE
      Read (11)  !  skip ZATMO
      Call READR4 (11,IM*JM,ZOCEAN,ZOCEAN)
      print*, "ANL!!",focean(1,4),zocean(1,4)
      Close (11)
C****
C**** Calculate arrays for ocean layering
C****
C**** Calculate LMOM from ZOE(L) and ZOCEAN (stored in ZSOLID)
      Do 330 J=1,JM
      Do 330 I=1,IM
      L=0
      If (FOCEAN(I,J) == 0)  GoTo 320
      Do 310 L=LMO_MIN,LMO-1
  310 If (ZOCEAN(I,J) <= .5*(ZOE(L)+ZOE(L+1)))  GoTo 320
C     L=LMO
  320 LMOM(I,J) = L
  330 ZOCEAN(I,J) = ZOE(L)
C**** Calculate MFO
      Do 370 LM=1,LMO
      MFO(:LM-1,LM) = Nint(1048576*dZO(:LM-1)/ZOE(LM)) / 1048576d0
      MFO(LM   ,LM) = 1 - Sum(MFO(:LM-1,LM))
  370 MFO(LM+1:,LM) = 0
C**** Read in tables of ocean functions
      Read  (22, Err=838) TITLE,VGSP
      Write (6,*) 'VGSP read from unit 22: ',TITLE
      Close (22)
  400 Return
C****
  838 Write (0,*) 'VGSP not read in from unit 22, should be OFUNTABLE'
      Return
      End

      Subroutine READR4 (IUNIT,NM,DATAR4,DATAR8)
C****
C**** READR4 reads a record from unit IUNIT containing TITLE,DATAR4,
C**** converts the Real*4 array DATAR4 to the Real*8 array DATAR8, and
C**** writes a line to unit 6 containing the TITLE just read.
C****
      Real*4 DATAR4(NM)
      Real*8 DATAR8(NM)
      Character*80 TITLE
C****
      Read (IUNIT) TITLE,DATAR4
      Do 10 I=NM,1,-1
   10 DATAR8(I) = DATAR4(I)
      Write (6,901) IUNIT,TITLE
      Return
C****
  901 Format (' Read on unit',I3,': ',A80)
      End

C**** HNTR8.FOR   Horizontal Interpolation Program Real*8   2007/11/13
C****
      Subroutine HNTR80 (IMA,JMA,OFFIA,DLATA,
     *                   IMB,JMB,OFFIB,DLATB, DATMIS)
C****
C**** HNTR80 fills in the common block HNTRCB with coordinate
C**** parameters that will be used by subsequent calls to HNTR8.
C**** The 5 Real input values are expected to be Real*8.
C****
C**** Input: IMA = number of cells in east-west direction of grid A
C****        JMA = number of cells in north-south direction of grid A
C****      OFFIA = number of cells of grid A in east-west direction
C****              from IDL (180) to western edge of cell IA=1
C****      DLATA = minutes of latitude for non-polar cells on grid A
C****        IMB = number of cells in east-west direction of grid B
C****        JMB = number of cells in north-south direction of grid B
C****      OFFIB = number of cells of grid B in east-west direction
C****              from IDL (180) to western edge of cell IB=1
C****      DLATB = minutes of latitude for non-polar cells on grid B
C****     DATMIS = missing data value inserted in output array B when
C****              cell (IB,JB) has integrated value 0 of WTA
C****
C**** Output: common block /HNTRCB/
C**** SINA(JA) = sine of latitude of northern edge of cell JA on grid A
C**** SINB(JB) = sine of latitude of northern edge of cell JB on grid B
C**** FMIN(IB) = fraction of cell IMIN(IB) on grid A west of cell IB
C**** FMAX(IB) = fraction of cell IMAX(IB) on grid A east of cell IB
C**** GMIN(JB) = fraction of cell JMIN(JB) on grid A south of cell JB
C**** GMAX(JB) = fraction of cell JMAX(JB) on grid A north of cell JB
C**** IMIN(IB) = western most cell of grid A that intersects cell IB
C**** IMAX(IB) = eastern most cell of grid A that intersects cell IB
C**** JMIN(JB) = southern most cell of grid A that intersects cell JB
C**** JMAX(JB) = northern most cell of grid A that intersects cell JB
C****
      Implicit Real*8 (A-H,O-Z)
      Parameter (TWOPI=6.283185307179586477d0)
      Real*8 OFFIA,DLATA, OFFIB,DLATB, DATMIS,DATMCB
      Common /HNTRCB/ SINA(0:5401),SINB(0:5401),
     *       FMIN(10800),FMAX(10800),GMIN(5401),GMAX(5401),
     *       IMIN(10800),IMAX(10800),JMIN(5401),JMAX(5401),
     *       DATMCB, INA,JNA, INB,JNB
C****
      INA = IMA  ;  JNA = JMA
      INB = IMB  ;  JNB = JMB
      DATMCB = DATMIS
      If (IMA<1 .or. IMA>10800 .or. JMA<1 .or. JMA>5401 .or.
     *    IMB<1 .or. IMB>10800 .or. JMB<1 .or. JMB>5401)  GoTo 400
C****
C**** Partitions in east-west (I) direction
C**** Domain, around the globe, is scaled to fit from 0 to IMA*IMB
C****
      DIA = IMB  !  width of single A grid cell in scaled domain
      DIB = IMA  !  width of single B grid cell in scaled domain
      IA  = 1
      RIA = (IA+OFFIA - IMA)*IMB  !  scaled longitude of eastern edge
      IB  = IMB
      Do 150 IBp1=1,IMB
      RIB = (IBp1-1+OFFIB)*IMA    !  scaled longitude of eastern edge
  110 If (RIA-RIB)  120,130,140
  120 IA  = IA  + 1
      RIA = RIA + DIA
      GoTo 110
C**** Eastern edges of cells IA of grid A and IB of grid B coincide
  130 IMAX(IB) = IA
      FMAX(IB) = 0
      IA  = IA  + 1
      RIA = RIA + DIA
      IMIN(IBp1) = IA
      FMIN(IBp1) = 0
      GoTo 150
C**** Cell IA of grid A contains western edge of cell IB of grid B
  140 IMAX(IB) = IA
      FMAX(IB) = (RIA-RIB)/DIA
      IMIN(IBp1) = IA
      FMIN(IBp1) = 1-FMAX(IB)
  150 IB = IBp1
      IMAX(IMB) = IMAX(IMB) + IMA
C****
C**** Partitions in the north-south (J) direction
C**** Domain is measured in minutes (1/60-th of a degree)
C****
      FJEQA = .5*(1+JMA)
      Do 210 JA=1,JMA-1
      RJA = (JA+.5-FJEQA)*DLATA  !  latitude in minutes of northern edge
  210 SINA(JA) = Sin (RJA*TWOPI/(360*60))
      SINA(0)  = -1
      SINA(JMA)=  1
C****
      FJEQB = .5*(1+JMB)
      Do 220 JB=1,JMB-1
      RJB = (JB+.5-FJEQB)*DLATB  !  latitude in minutes of northern edge
  220 SINB(JB) = Sin (RJB*TWOPI/(360*60))
      SINB(0)  = -1
      SINB(JMB)=  1
C****
      JMIN(1) = 1
      GMIN(1) = 0
      JA = 1
      Do 350 JB=1,JMB-1
  310 If (SINA(JA)-SINB(JB))  320,330,340
  320 JA = JA + 1
      GoTo 310
C**** Northern edges of cells JA of grid A and JB of grid B coincide
  330 JMAX(JB) = JA
      GMAX(JB) = 0
      JA = JA + 1
      JMIN(JB+1) = JA
      GMIN(JB+1) = 0
      GoTo 350
C**** Cell JA of grid A contains northern edge of cell JB of grid B
  340 JMAX(JB) = JA
      GMAX(JB) = SINA(JA) - SINB(JB)
      JMIN(JB+1) = JA
      GMIN(JB+1) = SINB(JB) - SINA(JA-1)
  350 Continue
      JMAX(JMB) = JMA
      GMAX(JMB) = 0
      Return
C****
C**** Invalid parameters or dimensions out of range
C****
  400 Write (0,940) IMA,JMA,OFFIA,DLATA, IMB,JMB,OFFIB,DLATB, DATMIS
      Stop 400
  940 Format ('0Arguments received by HNTRP0 in order:'/
     *   2I12,' = IMA,JMA = array dimensions for A grid'/
     *  E24.8,' = OFFIA   = fractional number of grid cells from',
     *                    ' IDL to western edge of grid cell I=1'/
     *  E24.8,' = DLATA   = minutes of latitude for interior grid cell'/
     *   2I12,' = IMB,JMB = array dimensions for B grid'/
     *  E24.8,' = OFFIB   = fractional number of grid cells from',
     *                    ' IDL to western edge of grid cell I=1'/
     *  E24.8,' = DLATB   = minute of latitude for interior grid cell'/
     *  E24.8,' = DATMIS  = missing data value to be put in B array',
     *                    ' when integrated WTA = 0'/
     *  '0These arguments are invalid or out of range.')
      End

      Subroutine HNTR8 (WTA,A,B)
C****
C**** HNTR8 performs a horizontal interpolation of per unit area or per
C**** unit mass quantities defined on grid A, calculating the quantity
C**** on grid B.  B grid values that cannot be calculated because the
C**** covering A grid boxes have WTA = 0, are set to the value DATMIS.
C**** The area weighted integral of the quantity is conserved.
C**** The 3 Real input values are expected to be Real*8.
C****
C**** Input: WTA = weighting array for values on the A grid
C****          A = per unit area or per unit mass quantity
C**** Output:  B = horizontally interpolated quantity on B grid
C****
      Implicit Real*8 (A-H,O-Z)
      Real*8 WTA(*), A(*), B(*), DATMIS
      Common /HNTRCB/ SINA(0:5401),SINB(0:5401),
     *       FMIN(10800),FMAX(10800),GMIN(5401),GMAX(5401),
     *       IMIN(10800),IMAX(10800),JMIN(5401),JMAX(5401),
     *       DATMIS, IMA,JMA, IMB,JMB
C****
C**** Interpolate the A grid onto the B grid
C****
      Do 20 JB=1,JMB
      JAMIN = JMIN(JB)
      JAMAX = JMAX(JB)
      Do 20 IB=1,IMB
      IJB  = IB + IMB*(JB-1)
      WEIGHT= 0
      VALUE = 0
      IAMIN = IMIN(IB)
      IAMAX = IMAX(IB)
      Do 10 JA=JAMIN,JAMAX
      G = SINA(JA)-SINA(JA-1)
      If (JA==JAMIN)  G = G - GMIN(JB)
      If (JA==JAMAX)  G = G - GMAX(JB)
      Do 10 IAREV=IAMIN,IAMAX
      IA  = 1 + Mod(IAREV-1,IMA)
      IJA = IA + IMA*(JA-1)
      F   = 1
      If (IAREV==IAMIN)  F = F - FMIN(IB)
      If (IAREV==IAMAX)  F = F - FMAX(IB)
      WEIGHT = WEIGHT + F*G*WTA(IJA)
   10 VALUE  = VALUE  + F*G*WTA(IJA)*A(IJA)
      B(IJB) = DATMIS
      If (WEIGHT.ne.0)  B(IJB) = VALUE/WEIGHT
   20 Continue
      Return
      End

      Subroutine HNTR8P (WTA,A,B)
C****
C**** HNTR8P is similar to HNTR8 but polar values are replaced by
C**** their longitudinal mean.
C**** The 3 Real input values are expected to be Real*8.
C****
      Implicit Real*8 (A-H,O-Z)
      Real*8 WTA(*), A(*), B(*)
      Common /HNTRCB/ SINA(0:5401),SINB(0:5401),
     *       FMIN(10800),FMAX(10800),GMIN(5401),GMAX(5401),
     *       IMIN(10800),IMAX(10800),JMIN(5401),JMAX(5401),
     *       DATMIS, IMA,JMA, IMB,JMB
C****
      Call HNTR8 (WTA,A,B)
C****
C**** Replace individual values near the poles by longitudinal mean
C****
      Do 40 JB=1,JMB,JMB-1
      BMEAN  = DATMIS
      WEIGHT = 0
      VALUE  = 0
      Do 10 IB=1,IMB
      IJB  = IB + IMB*(JB-1)
      If (B(IJB) == DATMIS)  GoTo 20
      WEIGHT = WEIGHT + 1
      VALUE  = VALUE  + B(IJB)
   10 Continue
      If (WEIGHT.ne.0)  BMEAN = VALUE/WEIGHT
   20 Do 30 IB=1,IMB
      IJB  = IB + IMB*(JB-1)
   30 B(IJB) = BMEAN
   40 Continue
      Return
      End

      Subroutine VLKtoLZ (KM,LM, MK,ME, RK, RL,RZ)        !  2008/02/05
C****
C**** VLKtoLZ assumes a continuous piecewise linear tracer distribution,
C**** defined by input tracer concentrations RK at KM specific points.
C**** MK in the downward vertical mass coordinate.
C**** R(M) = {RK(K-1)*[MK(K)-M] + RK(K)*[M-MK(K-1)]} / [MK(K)-MK(K-1)]
C****               when MK(K-1) < M < MK(K).
C**** R(M) = RK(1)  when M < MK(1).
C**** R(M) = is undefined when MK(KM) < M.
C****
C**** VLKtoLZ integrates this tracer distribution over the LM output
C**** layers defined by their layer edges ME, calculating the tracer
C**** mass RM of each layer and the vertical gradient RZ.
C**** RNEW(M) = RL(L) + RZ(L)*[M-MC(L)]/dM(L) when ME(L-1) < M < ME(L)
C**** where MC(L) = .5*[ME(L-1)+ME(L)] and dM(L) = ME(L)-ME(L-1)
C**** Mean concentration of output layers is RL(L) = RM(L)/dM(L).
C****
C**** If ME(L-1) < MK(KM) < ME(L), then RL(L) and RZ(L) are calculated
C**** from the input profile up to MK(KM); RL(L+1:LM) and RZ(L+1:LM)
C**** for deeper layers are undefined, set to DATMIS.
C****
C**** Input:  KM = number of input edges
C****         LM = number of output cells
C****         MK = mass coordinates of input points (kg/m^2)
C****         ME = mass coordinates of output layer edges (kg/m^2)
C****         RK = tracer concentration at input points
C****
C**** Output: RL = mean tracer concentration of each output layer
C****         RZ = vertical gradient of tracer mass of each output layer
C****
C**** Internal: RM = integrated tracer mass of output layers (kg/m^2)
C****           RQ = integrated tracer mass times mass (kg^2/m^4)
C****
      Implicit Real*8 (A-H,M-Z)
      Parameter (DATMIS = -999999)
      Real*8 MK(KM),ME(0:LM), RK(KM), RL(LM),RZ(LM), RM(1024),RQ(1024)
C     If (LM > 1024)  Stop 'LM exceeds internal dimentions in VLKtoLZ'
C****
      RM(1:LM) = 0
      RQ(1:LM) = 0
      K = 1
      L = 1
      MC = .5*(ME(L)+ME(L-1))
C****
C**** Integrate layers with M < MK(1)
C****
      If (ME(0) < MK(1))  GoTo 20
C**** MK(1) <= ME(0), determine K such that MK(K-1) <= ME(0) < MK(K)
   10 If (K == KM)  GoTo 200  ;  K = K+1
      If (MK(K) <= ME(0))  GoTo 10
      GoTo 130  !  MK(K-1) <= ME(0) < MK(K)
C**** ME(0) < MK(1), determine output cell containing MK(1)
   20 If (MK(1) < ME(L))  GoTo 30
C**** ME(L-1) < ME(L) < MK(1), integrate RM from ME(L-1) to ME(L)
      RM(L) = RK(1)*(ME(L)-ME(L-1))
      RQ(L) = 0
      If (L == LM)  GoTo 300  ;  L = L+1  ;  MC = .5*(ME(L)+ME(L-1))
      GoTo 20
C**** ME(L-1) < MK(1) < ME(L), integrate RM from ME(L-1) to MK(1)
   30 RM(L) = RK(1)*(MK(1)-ME(L-1))
      RQ(L) = RK(1)*(MK(1)-ME(L-1))*(.5*(MK(1)+ME(L-1))-MC)
      If (K == KM)  GoTo 220  ;  K = K+1
C****
C**** Integrate layers with MK(1) < M < MK(KM)
C****
  100 If (ME(L) < MK(K))  GoTo 120
C**** ME(L-1) < MK(K-1) < MK(K) < ME(L), integrate from MK(K-1) to MK(K)
      RM(L) = RM(L) + (RK(K)-RK(K-1))*(MK(K)+MK(K-1))/2 +
     +                RK(K-1)*MK(K)-RK(K)*MK(K-1)
      RQ(L) = RQ(L) +
     +  (RK(K)-RK(K-1))*(MK(K)*MK(K)+MK(K)*MK(K-1)+MK(K-1)*MK(K-1))/3 +
     +  (RK(K-1)*(MK(K)+MC)-RK(K)*(MK(K-1)+MC))*(MK(K)+MK(K-1))/2 +
     +  (RK(K)*MK(K-1)-RK(K-1)*MK(K))*MC
      If (K == KM)  GoTo 220  ;  K = K+1
      GoTo 100
C**** ME(L-1) < MK(K-1) < ME(L) < MK(K), integrate from MK(K-1) to ME(L)
  120 RM(L) = RM(L) + ((RK(K)-RK(K-1))*(ME(L)+MK(K-1))/2 +
     +                 (RK(K-1)*MK(K)-RK(K)*MK(K-1))) * (ME(L)-MK(K-1))
     /              / (MK(K)-MK(K-1))
      RQ(L) = RQ(L) +
     +  ((RK(K)-RK(K-1))*(ME(L)*ME(L)+ME(L)*MK(K-1)+MK(K-1)*MK(K-1))/3 +
     +   (RK(K-1)*(MK(K)+MC)-RK(K)*(MK(K-1)+MC))*(ME(L)+MK(K-1))/2 +
     +   (RK(K)*MK(K-1)-RK(K-1)*MK(K))*MC) * (ME(L)-MK(K-1)) /
     /  (MK(K)-MK(K-1))
      If (L == LM)  GoTo 300  ;  L = L+1  ;  MC = .5*(ME(L)+ME(L-1))
  130 If (MK(K) < ME(L))  GoTo 160
C**** MK(K-1) < ME(L-1) < ME(L) < MK(K), integrate from ME(L-1) to ME(L)
  140 RM(L) = ((RK(K)-RK(K-1))*(ME(L)+ME(L-1))/2 +
     +         (RK(K-1)*MK(K)-RK(K)*MK(K-1))) * (ME(L)-ME(L-1)) /
     /        (MK(K)-MK(K-1))
      RQ(L) =
     +  ((RK(K)-RK(K-1))*(ME(L)*ME(L)+ME(L)*ME(L-1)+ME(L-1)*ME(L-1))/3 +
     +   (RK(K-1)*(MK(K)+MC)-RK(K)*(MK(K-1)+MC))*(ME(L)+ME(L-1))/2 +
     +   (RK(K)*MK(K-1)-RK(K-1)*MK(K))*MC) * (ME(L)-ME(L-1)) /
     /  (MK(K)-MK(K-1))
      If (L == LM)  GoTo 300  ;  L = L+1  ;  MC = .5*(ME(L)+ME(L-1))
      If (ME(L) < MK(K))  GoTo 140
C**** MK(K-1) < ME(L-1) < MK(K) < ME(L), integrate from ME(L-1) to MK(K)
  160 RM(L) = RM(L) + ((RK(K)-RK(K-1))*(MK(K)+ME(L-1))/2 +
     +                 (RK(K-1)*MK(K)-RK(K)*MK(K-1))) * (MK(K)-ME(L-1))
     /              / (MK(K)-MK(K-1))
      RQ(L) = RQ(L) +
     +  ((RK(K)-RK(K-1))*(MK(K)*MK(K)+MK(K)*ME(L-1)+ME(L-1)*ME(L-1))/3 +
     +   (RK(K-1)*(MK(K)+MC)-RK(K)*(MK(K-1)+MC))*(MK(K)+ME(L-1))/2 +
     +   (RK(K)*MK(K-1)-RK(K-1)*MK(K))*MC) * (MK(K)-ME(L-1)) /
     /  (MK(K)-MK(K-1))
      If (K == KM)  GoTo 220  ;  K = K+1
      GoTo 100
C****
C**** Calculate RL and RZ from RM and RQ when MK(KM) < ME(LM)
C****
C**** MK(KM) <= ME(0)
  200 RL(:) = DATMIS
      RZ(:) = DATMIS
      Return
C**** ME(L-1) < MK(KM) < ME(L)
  220 Do 230 LL=1,L-1
      RL(LL) =   RM(LL) / (ME(LL)-ME(LL-1))
  230 RZ(LL) = 6*RQ(LL) / (ME(LL)-ME(LL-1))**2
      RL(L)  =   RM(L)  / (MK(KM)-ME(L-1))
      RZ(L)  = 6*(RQ(L) + .5*(ME(L)-MK(KM))*RM(L)) / (MK(KM)-ME(L-1))**2
C**** Vertical gradient is extrapolated half way to .5*[MK(KM)+ME(L)]
      RZ(L)  = RZ(L) * (.5*(MK(KM)+ME(L))-ME(L-1)) / (MK(KM)-ME(L-1))
      RL(L+1:LM) = DATMIS
      RZ(L+1:LM) = DATMIS
      Return
C****
C**** Calculate RL and RZ from RM and RQ when ME(LM) < MK(KM)
C****
  300 Do 310 L=1,LM
      RL(L) =   RM(L) / (ME(L)-ME(L-1))
  310 RZ(L) = 6*RQ(L) / (ME(L)-ME(L-1))**2
      Return
      End

C**** OFunction.SUB   Useful functions of ocean parameters   2007/11/20
C****
      Function VOLPTS (PIN,T,SIN)
C****
C**** VOLPTS calculates the specific volume of sea water as a
C**** function of pressure, temperature and salinity.
C**** The reference for this function is:
C**** "Tenth report of the joint panel on oceanographic tables and
C**** standards", Sidney, British Columbia, Canada, 1-5 September
C**** 1980, sponsored by Unesco, ICES, SCOR, IAPSO.
C**** Also see:
C**** N.P. Fofonoff, 1985.  Physical Properties of Seawater: A New
C**** Salinity Scale and Equation of State for Seawater.  Journal
C**** of Geophysical Research, volume 90, pp 3332-3342.
C****
C**** Input: PIN (Pa) = pressure above normal atmospheric pressure,
C****                   from 0 to 1.E8
C****          T (C)  = temperature, from -2 to 40
C****        SIN (1)  = salinity (kg NaCl/kg sea water), from 0 to .042
C****
C**** Output: VOLPTS (m^3/kg) = specific volume of sea water
C****
      Implicit Real*8 (A-Z)
      Data A0,A1,A2,A3,A4,A5 /999.842594, 6.793952D-2,
     *  -9.095290D-3, 1.001685D-4, -1.120083D-6, 6.536332D-9/
      Data B0,B1,B2,B3,B4 /8.24493D-1, -4.0899D-3, 7.6438D-5,
     *  -8.2467D-7, 5.3875D-9/
      Data C0,C1,C2 /-5.72466D-3, 1.0227D-4, -1.6546D-6/
      Data D0 /4.8314D-4/
      Data E0,E1,E2,E3,E4 /19652.21, 148.4206, -2.327105,
     *  1.360477D-2, -5.155288D-5/
      Data F0,F1,F2,F3 /54.6746, -.603459, 1.09987D-2, -6.1670D-5/
      Data G0,G1,G2 /7.944D-2, 1.6483D-2, -5.3009D-4/
      Data H0,H1,H2,H3 /3.239908, 1.43713D-3, 1.16092D-4, -5.77905D-7/
      Data I0,I1,I2 /2.2838D-3, -1.0981D-5, -1.6078D-6/
      Data J0 /1.91075D-4/
      Data K0,K1,K2 /8.50935D-5, -6.12293D-6, 5.2787D-8/
      Data M0,M1,M2 /-9.9348D-7, 2.0816D-8, 9.1697D-10/
C****
      P  = PIN*1.D-5
      S  = SIN*1.D3
      S32= S*Sqrt(S)
      KW = E0+(E1+(E2+(E3+E4*T)*T)*T)*T
      AW = H0+(H1+(H2+H3*T)*T)*T
      BW = K0+(K1+K2*T)*T
      KO = KW + (F0+(F1+(F2+F3*T)*T)*T)*S + (G0+(G1+G2*T)*T)*S32
      A  = AW + (I0+(I1+I2*T)*T)*S + J0*S32
      B  = BW + (M0+(M1+M2*T)*T)*S
      K  = KO + A*P + B*P**2
      DENSTW = A0+(A1+(A2+(A3+(A4+A5*T)*T)*T)*T)*T
      DENST0 = DENSTW + (B0+(B1+(B2+(B3+B4*T)*T)*T)*T)*S
     *                + (C0+(C1+C2*T)*T)*S32
     *                +  D0*S**2
      VOLPTS = (1.-P/K)/DENST0
      Return
      End

      Function SHCPTS (PIN,T,SIN)
C****
C**** SHCPTS calculates the specific heat capacity of sea water as
C**** a function of pressure, temperature and salinity.
C**** The reference for this function is:
C**** N.P. Fofonoff and R.C. Millard Jr., 1983.  Algorithms for
C**** Computation of Fundamental Properties of Seawater.  UNESCO
C**** Technical Papers in Marine Science, volume 44.
C**** Also see:
C**** N.P. Fofonoff, 1985.  Physical Properties of Seawater: A New
C**** Salinity Scale and Equation of State for Seawater.  Journal
C**** of Geophysical Research, volume 90, pp 3332-3342.
C****
C**** Input: PIN (Pa) = pressure above normal atmospheric pressure,
C****                   from 0 to 1.E8
C****          T (C)  = temperature, from 0 to 35
C****        SIN (1)  = salinity (kg NaCl/kg sea water), from 0 to .04
C****
C**** Output: SHCPTS (J/kg*C) = specific heat capacity of sea water
C****                           with standard deviation error of
C****                           .636 (J/C*kg)
C****
      Implicit Real*8 (A-Z)
      Data A000/ 4217.4    /, A001/-7.643575  /, A002/  .1770383 /,
     *     A010/-3.720283  /, A011/  .1072763 /, A012/-4.07718D-3/,
     *     A020/  .1412855 /, A021/-1.38385D-3/, A022/ 5.148D-5  /,
     *     A030/-2.654387D-3/,
     *     A040/ 2.093236D-5/,
     *     A100/-4.9592D-1 /, A101/ 4.9247D-3 /, A102/-1.2331D-4 /,
     *     A110/ 1.45747D-2/, A111/-1.28315D-4/, A112/-1.517D-6  /,
     *     A120/-3.13885D-4/, A121/ 9.802D-7  /, A122/ 3.122D-8  /,
     *     A130/ 2.0357D-6 /, A131/ 2.5941D-8 /,
     *     A140/ 1.7168D-8 /, A141/-2.9179D-10/,
     *     A200/ 2.4931D-4 /, A201/-2.9558D-6 /, A202/ 9.971D-8  /,
     *     A210/-1.08645D-5/, A211/ 1.17054D-7/,
     *     A220/ 2.87533D-7/, A221/-2.3905D-9 /,
     *     A230/-4.0027D-9 /, A231/ 1.8448D-11/,
     *     A240/ 2.2956D-11/,
     *     A300/-5.422D-8  /, A301/ 5.540D-10 /,
     *     A310/ 2.6380D-9 /, A311/-1.7682D-11/, A312/-1.4300D-12/,
     *     A320/-6.5637D-11/, A321/ 3.513D-13 /,
     *     A330/ 6.136D-13 /
C****
      P  = PIN*1.D-5
      S  = SIN*1.D3
      S32= S*Sqrt(S)
      SHCPTS = A000+(A010+(A020+(A030+A040*T)*T)*T)*T
     *  +   S*(A001+(A011+ A021                 *T)*T)
     *  + S32*(A002+(A012+ A022                 *T)*T)
     *  +     (A100+(A110+(A120+(A130+A140*T)*T)*T)*T
     *  +   S*(A101+(A111+(A121+(A131+A141*T)*T)*T)*T)
     *  + S32*(A102+(A112+ A122                 *T)*T)
     *  +     (A200+(A210+(A220+(A230+A240*T)*T)*T)*T
     *  +   S*(A201+(A211+(A221+ A231        *T)*T)*T)
     *  + S32* A202
     *  +     (A300+(A310+(A320+ A330        *T)*T)*T
     *  +   S*(A301+(A311+ A321                 *T)*T)
     *  + S32*       A312                          *T )*P)*P)*P
      Return
      End

      Function ATGPTS (PIN,T,SIN)
C****
C**** ATGPTS calculates the adiabatic lapse rate of sea water as
C**** a function of pressure, temperature and salinity.
C**** The reference for this function is:
C**** N.P. Fofonoff and R.C. Millard Jr., 1983.  Algorithms for
C**** Computation of Fundamental Properties of Seawater.  UNESCO
C**** Technical Papers in Marine Science, volume 44.
C****
C**** Input: PIN (Pa) = pressure above normal atmospheric pressure
C****          T (C)  = temperature
C****        SIN (1)  = salinity (kg NaCl/kg sea water)
C****
C**** Output: ATGPTS (C/Pa) = adiabatic lapse tate of sea water, at
C****                         S = .035, error < .006 (C) when used
C****                         to calculate potential temperature
C****
      Implicit Real*8 (A-Z)
      Data A000/ 3.5803D-5 /, A010/ 8.5258D-6 /, A020/-6.8360D-8 /,
     *                        A030/ 6.6228D-10/,
     *     A001/ 1.8932D-6 /, A011/-4.2393D-8 /,
     *     A100/ 1.8741D-8 /, A110/-6.7795D-10/, A120/ 8.7330D-12/,
     *                        A130/-5.4481D-14/,
     *     A101/-1.1351D-10/, A111/ 2.7759D-12/,
     *     A200/-4.6206D-13/, A210/ 1.8676D-14/, A220/-2.1687D-16/
C****
      P = PIN*1.D-4
      S = SIN*1.D3 - 35.
      ATGPTS = A000+(A010+(A020+A030*T)*T)*T
     *    + S*(A001+ A011*T)
     *    +   (A100+(A110+(A120+A130*T)*T)*T
     *    + S*(A101+ A111*T)
     *    +   (A200+(A210+ A220        *T)*T)*P)*P
      ATGPTS = ATGPTS*1.D-4
      Return
      End

      Function PTPTS (P,T,S)
C****
C**** PTPTS calculates the potential temperature of sea water as
C**** a function of pressure, temperature and salinity.
C**** At pressures above 0, PTPTS solves the differential equation:
C**** dT/dP = ATG = TK/SHC * dSVOL/dT .
C**** At pressure = 0, PTPTS = T.
C****
C**** Input: P (Pa) = pressure above normal atmospheric pressure
C****        T (C)  = temperature
C****        S (1)  = salinity (kg NaCl/kg sea water)
C****
C**** Output: PTPTS (C) = potential temperature of sea water,
C****                     with maximum error of .004 (C) ?
C****
      Implicit Real*8 (A-H,O-Z)
      NM = 1 + Abs(P)/2.D6
      DP = P/NM
      T0 = T
      DO 10 N=NM,1,-1
      P0 = N*DP
      T1 = T0 - .5*DP*ATGPTS(P0-.25*DP,T0,S)
   10 T0 = T0 -    DP*ATGPTS(P0-.50*DP,T1,S)
      PTPTS = T0
      Return
      End

      Function DELHTS (T,SIN)
C****
C**** DELHTS calculates the change of specific enthalpy of sea
C**** water as salinity changes from 0 to an input value, as a
C**** function of temperature at atmospheric pressure.
C**** The reference for this function is:
C**** Frank J. Millero and Wing H. Leung, 1976.  The Thermodynamics
C**** of Seawater at One Atmosphere.  American Journal of Science,
C**** volume 276.
C****
C**** Input: T (C) = temperature, from -2 to 40
C****        S (1) = salinity (kg NaCl/kg sea water), from 0 to .04
C****
C**** Output: DELHTS (J/kg) = change of specific heat of sea water
C****
      Implicit Real*8 (A-Z)
      Data A01/ 3.4086D-3/, A03/ 7.9350D-4/, A02/-4.7989D-4/,
     *     A11/-6.3798D-5/, A13/ 1.0760D-4/, A12/ 6.3787D-6/,
     *     A21/ 1.3877D-6/, A23/-6.3923D-7/, A22/-1.1647D-7/,
     *     A31/-1.0512D-8/, A33/ 8.60D-9  /, A32/ 5.717D-10/
C****
      S = SIN*1.D3
      DELHTS = (A01+(A11+(A21+A31*T)*T)*T
     *       + (A03+(A13+(A23+A33*T)*T)*T)*Sqrt(S)
     *       + (A02+(A12+(A22+A32*T)*T)*T)*S)*S*1.D3
      Return
      End

      Function HETPTS (P,T,S)
C****
C**** HETPTS calculates the specific enthalpy of sea water as a
C**** function of pressure, temperature and salinity.
C**** At pressures above 0, HETPTS solves the differential equation:
C**** dH/dP = V .
C**** At pressure = 0 and salinities above 0,
C**** H(0,T,S) = H(0,T,0) + DELHTS(T,S)
C**** At pressure = 0 and salinity = 0, HETPTS solves the differential
C**** equation:  dH/dT = C .
C****
C**** Input: P (Pa) = pressure above normal atmospheric pressure
C****        T (C)  = temperature
C****        S (1)  = salinity (kg NaCl/kg sea water)
C****
C**** Output: HETPTS (J/kg) = specific enthalpy of sea water
C****
      Implicit Real*8 (A-H,O-Z)
C****
C**** Calculate H(P,T,S) - H(0,T,S) by integrating  dH/dP = V  at
C**** constant entropy and salinity
C****
      H0 = 0.
      T0 = T
      If (P==0)  GoTo 20
      NM = 1 + Abs(P)/2.D6
      DP = P/NM
C**** For first step, integrate T down 1/2 DP, and H down full DP
      P1 = DP*(NM-.5)
      TX = T0 - .25*DP*ATGPTS(P1+.375*DP,T0,S)
      T1 = T0 - .50*DP*ATGPTS(P1+.250*DP,TX,S)
      H0 = H0 +     DP*VOLPTS(P1,T1,S)
C**** For subsequent steps, integrate T and H down full DP
      Do 10 N=NM-1,1,-1
      P1 = DP*(N-.5)
      T0 = T1 - .5*DP*ATGPTS(P1+.75*DP,T1,S)
      T1 = T1 -    DP*ATGPTS(P1+.50*DP,T0,S)
   10 H0 = H0 +    DP*VOLPTS(P1,T1,S)
C**** For last step, integrate T down 1/2 DP
      TX = T1 - .25*DP*ATGPTS(.375*DP,T1,S)
      T0 = T1 - .50*DP*ATGPTS(.250*DP,TX,S)
C****
C**** Calculate H(P,T,S) - H(0,T,0) = H(P,T,S) - H(0,T,S) + DELHTS(T,S)
C****
   20 H0 = H0 + DELHTS(T0,S)
C****
C**** Calculate H(0,T,0) - H(0,0,0) by integrating  dH/dT = C  at
C**** constant pressure and salinity
C****
      NM = 1 + Abs(T0)
      DT = T0/NM
      Do 30 N=0,NM-1
      T1 = DT*(N+.5)
   30 H0 = H0 + DT*SHCPTS(0.D0,T1,0.D0)
C****
      HETPTS = H0
      Return
      End

      Function PHPTS (P,T,S)
C****
C**** PHPTS calculates the potential specific enthalpy of sea water
C**** as a function of pressure, temperature and salinity.
C****
C**** Input: P (Pa) = pressure above normal atmospheric pressure
C****        T (C)  = temperature
C****        S (1)  = salinity (kg NaCl/kg sea water)
C****
C**** Output: PHPTS (J/kg) = potential specific enthalpy of sea water
C****
      Implicit Real*8 (A-H,O-Z)
      A = PTPTS(P,T,S)
      PHPTS = HETPTS(0.D0,A,S)
      Return
      End

C**** OFUNTABLE.SUB   TABLE FUNctions of Ocean parameters   2007/11/20
C****
      Function VOLGSP (G,S,P)
C****
C**** VOLGSP returns a linearly interpolated specific volume from
C**** an input table that depends on potential specific enthalpy,
C**** salinity, and pressure.
C****
C**** Input: G (J/kg) = potential specific enthalpy,
C****                   from -8000 to 160000
C****           S (1) = salinity (kg NaCl/kg sea water), from 0 to .04
C****          P (Pa) = pressure above normal atmospheric pressure,
C****                   from 0 to 78.E6
C****
C**** Output: VOLGSP (m**3/kg) = specific volume of sea water
C****
      Implicit Real*8 (A-H,O-Z)
      Common /OFUNCB/ V(-2:40,0:40,0:39)
C****
      GG = G/4000.
      SS = S*1000.
      PP = P/2.D6
      IG = Int(GG+2.) - 2
      If (IG < -2)  IG = -2
      If (IG > 39)  IG = 39
      JS = SS
C     If (JS <  0)  JS =  0
      If (JS > 39)  JS = 39
      KP = PP
      If (KP <  0)  KP =  0
      If (KP > 38)  KP = 38
C****
      VOLGSP = (KP-PP+1)*((JS-SS+1)*((IG-GG+1)*V(IG  ,JS  ,KP  )
     *                             + (GG-IG  )*V(IG+1,JS  ,KP  ))
     *                  + (SS-JS  )*((IG-GG+1)*V(IG  ,JS+1,KP  )
     *                             + (GG-IG  )*V(IG+1,JS+1,KP  )))
     *       + (PP-KP  )*((JS-SS+1)*((IG-GG+1)*V(IG  ,JS  ,KP+1)
     *                             + (GG-IG  )*V(IG+1,JS  ,KP+1))
     *                  + (SS-JS  )*((IG-GG+1)*V(IG  ,JS+1,KP+1)
     *                             + (GG-IG  )*V(IG+1,JS+1,KP+1)))
      Return
      End

