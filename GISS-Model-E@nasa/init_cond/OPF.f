C**** OPF.E.CRE   CREate Ocean Polar Filter for Model E   2008/07/17
C****
C**** Fortan Compile and Go:  FCG OPF.E.CRE
C****
C**** To change to different resolution:
C**** Set: IM,JM,LMO,J1O,LMO_MIN (in MAIN and in SGEOM)
C**** Set: NBASM and INDM        (in MAIN)
C**** Set: dZO(LMO)              (in MAIN)
C**** Set: FILEIZ and FILOUT     (in MAIN)
C****
      Program MAIN
      Implicit  Real*8 (A-H,O-Z)
      Integer*4,Parameter ::
     *  IM   = 144, !  number of grid cells in longitude
     *  JM    = 90, !  number of grid cells in latitude
     *  LMO   = 32, !  maximum number of ocean layers
     *  J1O   =  4, !  most southern latitude cell where ocean exists
     *  LMO_MIN= 2, !  minimum number of ocean layers over ocean
     *  NBASM = 10, !  maximum number of ocean basins at any latitude
c     *  INDM  = 941701,    !  number of Real*4 coefficients in REDUCO
     *  INDM  = 1476847,    !  number of Real*4 coefficients in REDUCO
     *  JLATD = 180/(JM-1) !  LATitudinal spacing in Degrees
      Real*8,Parameter ::
c     *  dZO(LMO) = (/ 12d0, 18d0, 27d0, 40.5d0, 60.75d0, 91.125d0,
c     *                136.6875d0, 205.03125d0, 307.546875d0,
c     *                461.3203125d0, 691.98046875d0,
c     *                1037.970703125d0, 1556.9560546875d0 /)  !  13 Lay
     *  dZO(LMO) = (/ 12, 18, 26, 36,  48, 62, 78, 96,
     *               116,134,150,164, 176,186,184,200,
     *               204,206,206,206, 206,206,206,206,
     *               206,206,206,206, 206,206,206,206 /)  !  32 Layers
      Character*80,Parameter ::
     *  FILEIZ = 'Z144X90N_nocasp.1',  !  Z input file
     *  FILOUT = 'OPF.E2HX2.L32'       !  OPF file created
c     *  FILOUT = 'OPF.E2HX2.L13'       !  OPF file created
      Character*80 ::
c     *  TITLE1 = 'Ocean Polar Filter, 2Hx2, L13, Z file: ',  !  output
     *  TITLE1 = 'Ocean Polar Filter, 2Hx2, L32, Z file: ',  !  output
     *  TITLE2 = 'Ocean Polar Filter Arrays for Z file: '    !  TITLEs
C****
C**** The size of the output file depends on NBASM and INDM.  Start
C**** with large values for these two parameters.  The printer output
C**** file of this program will show what values are needed.  You may
C**** rerun this program with those values for NBASM and INDM.
C****
C**** The ocean dynamics time step is chosen as the largest convenient
C**** time step such that gravity waves of length 2*dY do not cause
C**** the numerical solution of the momentum equation to diverge.
C**** Shorter gravity waves that are resolved at high latitudes in the
C**** zonal direction could cause the numerical solution to diverge.
C**** To prevent this, certain fields are spectrally analyzed and
C**** coefficients of waves shorter than 2*dY are multiplied by the
C**** wave length divided by 2*dY.
C****
C**** For an ocean basin that is IWID grid cells wide, the fields of
C**** interest are defined on IWID-1 grid cell edges.  The fields are
C**** spectrally analyzed on 2*IWID periodic values which include 0's
C**** on each coast and IWID-1 reflected values of opposite sign.
C**** The cosine spectral coefficients are all zero, but the sine
C**** coefficients for wave numbers 1 to IWID-1 are generally nonzero.
C**** The shortest wave (wave number IWID-1) in a basin that is IWID
C**** cells wide has length 2*IWID*dX/(IWID-1).
C****
C**** Since filtered field values are fixed linear combinations of
C**** the IWID-1 unfiltered values, the matrix coefficents of linear
C**** combinations are calculated once in this program, but are used
C**** repetively in the ocean polar filter subroutine.  The matrix
C**** coefficients depend upon IWID and latitude (dX).
C****
C**** The polar filter is not applied at the poles (J = 1 or JM) nor
C**** between 40S and 40N.  The polar filter is applied from J1O to
C**** JSMAX and from JNMIN to JM-1.
C****
      Real*8,Parameter :: TWOPI = 6.283185307179586477d0,
     *  FJEQ = .5*(1+JM)  !  J coordinate of EQuator
      Integer*4,Parameter ::
     *  JSMAX = FJEQ-40./JLATD, !  most northern latitude in SH (40S)
                                !  where polar filter is applied
     *  JNMIN = JM+1-JSMAX      !  most southern latitude in NH (40N)
                                !  where polar filter is applied
     * ,JOMAX = JM-1-(JNMIN-JSMAX-1)  !  = JM-1 - unfiltered J's
      Integer*2 NBAS(LMO,2:JM-1), !  number of ocean basins
     *    IMIN(NBASM,LMO,2:JM-1), !  western most cell in basin
     *    IMAX(NBASM,LMO,2:JM-1), !  eastern most cell in basin
     *    IWID(NBASM,LMO,2:JM-1)  !  number of cells in basin
     * ,IMINm1(NBASM,LMO,J1O:JOMAX), !  western most cell in basin - 1
     *  IWIDm1(NBASM,LMO,J1O:JOMAX)  !  number of cells in basin - 1
      Integer*4 LMOM(IM,JM), !  number of ocean layers for each column
     *       IWMIN(2:JSMAX), !  minimum ocean width for polar filter
     *    INDEX(IM,2:JSMAX)  !  index to concatenated reduction maticies
      Real*4 FOCEAN(IM,JM),  !  horizontal ocean fraction
     *       ZOCEAN(IM,JM),  !  water depth (m)
     *        REDUCO(INDM)   !  concatenation of reduction matricies
      Real*8 ZOE(0:LMO),     !  approximate ocean depths at layer edges
     *       REDUC(IM,IM)    !  single polar filter reduction matrix
      Character*80  TITLE
      Character*256 IFDIR,FILEIN
      Common /GEOMCB/ DXYP(JM),DXP(JM),DYP(JM),DXV(JM)
C**** Diagnostic output array CBAS assumes JSMAX <= 60
      Character*1,Parameter ::
     *  CBAS(2:60) = (/ '2','3','4','5', '6','7','8','9','o',
     *              '1','2','3','4','5', '6','7','8','9','o',
     *              '1','2','3','4','5', '6','7','8','9','o',
     *              '1','2','3','4','5', '6','7','8','9','o',
     *              '1','2','3','4','5', '6','7','8','9','o',
     *              '1','2','3','4','5', '6','7','8','9','o' /)
      Character*1  QBAS(IM,2:JSMAX)  !  existence of basin, wide enough?
C****
      Call SGEOM
      TITLE1 = Trim(TITLE1) // ' ' // FILEIZ
      TITLE2 = Trim(TITLE2) // ' ' // FILEIZ
      Call GetEnv ('IFDIR',IFDIR)
      FILEIN = Trim(IFDIR) // '/' // FILEIZ
      QBAS(:,:) = ' '
      ZOE(0) = 0
      Do 10 L=1,LMO
   10 ZOE(L) = ZOE(L-1) + dZO(L)
      Write (6,901) ZOE
C****
C**** Read in FOCEAN and ZOCEAN
C****
      Open  (11, File=FILEIN, Form='Unformatted', Status='Old', Err=801)
      Read  (11) TITLE,FOCEAN
      Write (6,*) 'FOCEAN read from FILEIZ: ',TITLE
      Read  (11)  !  skip FLAKE
      Read  (11)  !  skip FGRND
      Read  (11)  !  skip FGICE
      Read  (11)  !  skip ZATMO
      Read  (11) TITLE,ZOCEAN
      Write (6,*) 'ZOCEAN read from FILEIZ: ',TITLE
      Close (11)
C****
C**** Calculate LMOM
C****
      Do 40 J=1,JM
      Do 40 I=1,IM
      LMOM(I,J) = 0
      If (FOCEAN(I,J) <= 0)  GoTo 40
      Do 20 L=LMO_MIN,LMO-1
   20 If (ZOCEAN(I,J) <= .5*(ZOE(L)+ZOE(L+1)))  GoTo 30
C     L = LMO
   30 LMOM(I,J) = L
   40 Continue
C****
C**** Calculate minimum number of cells in ocean basin for polar filter
C**** to be applied.  The length of the shortest wave in this minimum
C**** basin is 2*IWMIN*DXP(J)/(IWMIN-1) which must be less than 2*DYP(3)
C****
      IWMIN(2:JSMAX) = Ceiling (DYP(3) / (DYP(3)-DXP(2:JSMAX)))
      Write (6,912) (J,J=2,JSMAX)
      Write (6,913) IWMIN(2:JSMAX)
C****
C**** Determine ocean basins for each layer and latitude
C****   NBAS(L,J) = number of ocean basins at this layer and latitude
C**** IMIN(N,L,J) = western most ocean cell for N-th ocean basin
C**** IMAX(N,L,J) = eastern most ocean cell for N-th ocean basin
C**** IWID(N,L,J) = number of ocean cells for N-th ocean basin
C**** QBAS(WID,J) = existence of ocean basin with given WIDth
C****
      NBAS(:,:)   = 0
      IMIN(:,:,:) = 0
      IWID(:,:,:) = 0
      NBASN = 0
      Do 310 J=2,JM-1
      If (J > JSMAX .and. J < JNMIN)  GoTo 310
      JA=J  ;  If(J > JM/2) JA=JM+1-J  !  Absolute latitude, reflect NH
      Do 300 L=1,LMO                                     !  cells to SH
      N = 0
C**** Locate the first land cell
      IL1=0
      If (LMOM(IM,J) < L)  GoTo 220
      Do 210 IL1=1,IM-1
      If (LMOM(IL1,J) < L)  GoTo 220
  210 Continue
C**** All grid cells are ocean cells
      N = 1
      IMIN(N,L,J) = 1
      IMAX(N,L,J) = IM
      IWID(N,L,J) = IM
      QBAS(IM,JA)  = 'Y'
      GoTo 290
C**** IL1 is a land cell, find first subsequent ocean cell
  220 Do 230 IW1=IL1+1,IM
  230 If (LMOM(IW1,J) >= L)  GoTo 240
C**** There are no more ocean cells
      GoTo 290
C**** IW1 is the western most cell of an ocean basin,
C**** find the eastern most cell of this basin
  240 Do 250 IL1=IW1+1,IM
  250 If (LMOM(IL1,J) < L)  GoTo 270
C**** Cell (IM,J) is an ocean cell,
C**** find the eastern most cell of this basin, search past IM
      Do 260 IL1=1,IM
  260 If (LMOM(IL1,J) < L)  GoTo 280
      Stop 'OPF: 260 Inconsistancy'
C**** IL1-1 is the eastern most cell of the present ocean basin
  270 IWIDE = IL1 - IW1
      If (IWIDE < IWMIN(JA)) Then
        QBAS(IWIDE,JA) = 'N'  ;  GoTo 220  ;  EndIf  !  basin too narrow
      N = N+1
      IMIN(N,L,J) = IW1
      IMAX(N,L,J) = IL1-1
      IWID(N,L,J) = IWIDE
      QBAS(IWIDE,JA) = 'Y' !  Yes, basin is sufficiently wide
      GoTo 220
C**** IL1-1 is the eastern most cell of the present ocean basin,
C**** but this is the last ocean basin because IL1 crossed IM
  280 IWIDE = IM+IL1 - IW1
      If (IWIDE < IWMIN(JA)) Then
        QBAS(IWIDE,JA) = 'N'  ;  GoTo 290  ;  EndIf  !  basin too narrow
      N = N+1
      IMIN(N,L,J) = IW1
      IMAX(N,L,J) = IL1-1
      IWID(N,L,J) = IWIDE
      QBAS(IWIDE,JA) = 'Y' !  Yes, basin is sufficiently wide
  290 NBAS(L,J) = N
      Write (6,929) J,L,N,(IMIN(K,L,J),IMAX(K,L,J),K=1,N)
      If (N > NBASN)  NBASN = N
      If (NBASN > NBASM)  GoTo 829
  300 Continue
  310 Continue
      Write (6,*)
      Write (6,*) 'NBASN needed =',NBASN
      Write (6,*) 'NBASM input  =',NBASM
      Write (0,*) 'NBASN needed =',NBASN
      Write (0,*) 'NBASM input  =',NBASM
      Write (6,931) (CBAS(JA),JA=2,JSMAX)
      Do 320 I=1,IM
  320 Write (6,932) I,(QBAS(I,JA),JA=2,JSMAX),'+'
C****
C**** Loop over absolute latitudes and ocean basin widths.
C**** Determine starting locations in the concatenated reduction matrix,
C**** calculate entries in each reduction matrix, and load them into
C**** the concatenated array.
C**** INDEX(IWID,JA) = starting location of matrix in REDUNO
C**** REDUCO(INDEX)  = concatenation of polar filter reduction matricies
C****
      INDX=0
      Do 480 JA=2,JSMAX
      Do 480 IWIDE=2,IM-1
      If (QBAS(IWIDE,JA).ne.'Y')  GoTo 480
C****
C**** Calculate the local reduction matrix for given IWIDE and JA.
C**** A(I) is defined on grid cell edges of which there are IWIDE-1.
C**** A(I) = A(I) - sum REDUC(K,I)*A(K)  for  K = 1,IWIDE-1
C****
      KM = 2*IWIDE
C**** Zero out local reduction factor matrix
      Do 410 I=1,IWIDE-1
      Do 410 K=1,IWIDE-1
  410 REDUC(K,I) = 0
C**** Loop over wave numbers, N, from IWIDE-1 down to 1.
C**** Calculate the reduction for each wave number and the grid cell
C**** contribution of it from all other grid cells.
      Do 440 N=IWIDE-1,1,-1
      REDUCN = 1 - DXP(JA)*IWIDE/(DYP(3)*N)                 !  ocean
CATMO REDUCN = 1 - DXP(JA)/(DYP(3)*Sin(TWOPI*N*.25/IWIDE))  !  atmos
      If (REDUCN <= 0)  GoTo 440
      Do 420 I=1,IWIDE-1
      Do 420 K=1,IWIDE-1
  420 REDUC(K,I) = REDUC(K,I) +
     +  Sin(TWOPI*N*K/KM)*(Sin(TWOPI*N*I/KM)*REDUCN*4/KM)
C     Do 430 I=IWIDE/2+1,IWIDE
C     Do 430 K=1,IWIDE-1
C 430 REDUC(K,I) = REDUC(IWIDE-K,IWIDE-I)
  440 Continue
      If (IWIDE > 6)  GoTo 460
      Write (6,944) JA,IWIDE
      Do 450 I=1,IWIDE-1
  450 Write (6,945) (REDUC(K,I),K=1,IWIDE-1)
  460 Continue
C**** Load the reduction contribution matrix into the large
C**** coefficient array
      INDEX(IWIDE-1,JA) = INDX
      Do 470 I=1,IWIDE-1
      Do 470 K=1,IWIDE-1
      INDX=INDX+1
      If (INDX  > INDM)  GoTo 847
  470 REDUCO(INDX) = REDUC(K,I)
  480 Continue
      Write (6,*)
      Write (6,*) 'INDM needed =',INDX
      Write (6,*) 'INDM input  =',INDM
      Write (0,*) 'INDM needed =',INDX
      Write (0,*) 'INDM input  =',INDM
      Write (0,*) 'Rerun OPF??.CRE program substituting into code' //
     *            ' needed values for NBASM and INDM'
C****
C**** Remove unfiltered latitudes between 40S and 40N from the arrays
C**** IMIN and IWID and subtract 1 to produce standard format Model E
C****
      Do 510 J=J1O,JSMAX
      Do 510 L=1,LMO
      Do 510 N=1,NBAS(L,J)
      IMINm1(N,L,J) = IMIN(N,L,J) - 1
      If (IWID(N,L,J) < IM)
     *  Then  ;  IWIDm1(N,L,J) = IWID(N,L,J) - 1
        Else  ;  IWIDm1(N,L,J) = IM  ;  EndIf
  510 Continue
      JOFF = JNMIN - JSMAX - 1  !  number of discarded latitude bands
      Do 520 J=JNMIN,JM-1
      Do 520 L=1,LMO
      Do 520 N=1,NBAS(L,J)
      IMINm1(N,L,J-JOFF) = IMIN(N,L,J) - 1
      If (IWID(N,L,J) < IM) 
     *  Then  ;  IWIDm1(N,L,J-JOFF) = IWID(N,L,J) - 1
        Else  ;  IWIDm1(N,L,J-JOFF) = IM  ;  EndIf
  520 Continue
C****
C**** Write locating parameters and concatenated reduction matricies
C****
      Open  (2, File=FILOUT, Form='Unformatted')
      Write (2) TITLE1, NBASN,INDX
C     Write (2) TITLE2, NBAS, IMIN(:NBASN,:,:), IWID(:NBASN,:,:),
C    *          INDEX, REDUCO(:INDX)
      Write (2) TITLE2, NBAS(:,J1O:JSMAX),NBAS(:,JNMIN:JM-1), 
     *          IMINm1(:NBASN,:,:), IWIDm1(:NBASN,:,:),
     *          INDEX, REDUCO(:INDX)
      Close (2)
      Write (6,*)
      Write (6,*) 'Output 2: ',Trim(TITLE1)
      Write (6,*) '          ',Trim(TITLE2)
      GoTo 999
C****
  801 Write (0,*) 'Error opening FILEIZ: ',FILEIN
      Write (0,*) 'Was environment variable IFDIR set in your .profile'
      Stop 801
  829 Write (0,*) 'NBASN exceeds parameter NBASM =',NBASM
      Stop 829
  847 Write (0,*) 'INDX exceeds parameter INDM =',INDM
      Stop 847
C****
  901 Format (/' ZOE =',8F8.1 / (6X,8F8.1))
  912 Format (/'     JA = ',32I4)
  913 Format ( ' WIDmin = ',32I4 /)
  929 Format (' J,L,NBAS=',3I3,12(I5,I3))
  931 Format (/' Existence of Ocean Basin' / 4X,'JA=',120A1)
  932 Format (I4,'  +',120A1)
  944 Format ('0J,IWID=',2I3,5X,'Polar Filter Reduction Matrix')
  945 Format (1X,13F10.4)
  999 End

      Subroutine SGEOM
C****
C**** SGEOM calculates the spherical geometry for the C grid
C****
      Implicit Real*8 (A-H,O-Z)
      Integer*4,Parameter ::
     *  IM    = 144,       !  number of grid cells in longitude
     *  JM    =  90,       !  number of grid cells in latitude
     *  JLATD = 180/(JM-1) !  LATitudinal spacing in Degrees
      Real*8,Parameter :: TWOPI=6.283185307179586477d0,
     *  RADIUS = 6371000,
     *  FJEQ   = .5*(1+JM),
     *  DLON   = TWOPI/IM,
     *  DLAT   = TWOPI*JLATD/360
      Common /GEOMCB/ DXYP(JM),DXP(JM),DYP(JM),DXV(JM)
C****
C**** Calculates the spherical geometry for the C grid
C****
C**** Geometric parameters centered at secondary latitudes
      Do 10 J=1,JM-1
      SINS = Sin (DLAT*(J  -FJEQ))
      SINN = Sin (DLAT*(J+1-FJEQ))
   10 DXV(J)  = DLON*RADIUS*(SINN-SINS)/DLAT
C**** Geometric parameters centered at primary latitudes
      Do 20 J=2,JM-1
      SINS    = Sin (DLAT*(J-.5-FJEQ))
      SINN    = Sin (DLAT*(J+.5-FJEQ))
      DXYP(J) = DLON*RADIUS*RADIUS*(SINN-SINS)
      DXP(J)  = .5*(DXV(J-1)+DXV(J))
   20 DYP(J)  = DLAT*RADIUS
      DXYP(JM)= DLON*RADIUS*RADIUS*(1.-SINN)
      DXYP(1) = DXYP(JM)
      DXP(1)  = .5*DXV(1)
      DXP(JM) = .5*DXV(JM-1)
      DYP(1)  = .5*DLAT*RADIUS
      DYP(JM) = .5*DLAT*RADIUS
      Return
      End
