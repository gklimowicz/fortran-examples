C**** EXCHANGE.F    Cube - Sphere intersection points    2008/09/10
C****
      Module LLvsCS  !  Longitude-Latitude grid vs Cube-Sphere grid
C****
C**** Sphere of radius 1 is embedded inside a cube [-1:1,-1:1,-1:1].
C**** O (lOngitude) and A (lAtitude) are fixed spherical coordinates.
C**** X,Y,Z are fixed rectangular coordinates; the X axis passes
C**** through O=0 and A=0, and the Z axis intersects the poles.
C****
C**** Faces of the cube are ordered 1 to 6 as follows:
C**** K=1: X=+1, (I,J) measures (+Y,+Z), face is tangent at A=0, O=0.
C**** K=2: Y=+1, (I,J) measures (-X,+Z), face is tangent at A=0, O=90.
C**** K=3: Z=+1, (I,J) measures (-X,-Y), face is tangent at A=90.
C**** K=4: X=-1, (I,J) measures (-Z,-Y), face is tangent at A=0, O=180.
C**** K=5: Y=-1, (I,J) measures (-Z,+X), face is tangent at A=0, O=-90.
C**** K=6: Z=-1, (I,J) measures (+Y,+X), face is tangent at A=90.
C****
C**** G = index of longitude cell may be edge (0:GM) or center (1:GM)
C**** H = index of latitude cell may be edge (0:HM) or center (1:HM)
C**** I = index on cube or sphere may be edge (0:IM) or center (1:IM)
C**** J = index on cube or sphere may be edge (0:IM) or center (1:IM)
C**** K = face number on cube or projected onto sphere (1:6)
C****
      Integer*4,Parameter ::
c     *    GM = 72, !  number of cells in longitude, multiple of 8
c     *    HM = 46, !  number of cells in latitude, must be even
c     *    GM = 144, !  number of cells in longitude, multiple of 8
c     *    GM = 288, !  number of cells in longitude, multiple of 8
c     *    HM = 180, !  number of cells in latitude, must be even
     *    GM = 360, !  number of cells in longitude, multiple of 8
     *    HM = 180, !  number of cells in latitude, must be even
     *    IM = 90, !  number of linear cells on cube face, must be even
c     *  NMX1 = 32, !  maximum number of Ns for line intersections
     *  NMX1 = 96, !  maximum number of Ns for line intersections
c     *  NMX2 = 256,   !  maximum number of Ns for polygons in cell
     *  NMX2 = 1024,   !  maximum number of Ns for polygons in cell
     *  NMAX = 300000  !  maximum total number of polygons on sphere
      Real*16,Parameter :: PRECIS=1q-20
      Real*16 :: TWOPIQ,  !  2*pi (quartic precision)
C**** Geometric parameters for the Longitude-Latitude grid
     *  OofG(0:GM),  !  longitude edges of Lon-Lat grid cells (radians)
     *  AofH(0:HM),  !  latitude edges of Lon-Lat grid cells (radians)
     *  SinOofG(0:GM), CosOofG(0:GM), TanOofG(0:GM),
     *  CscOofG(0:GM), SecOofG(0:GM), CotOofG(0:GM),
     *  SinAofH(0:HM), TanAofH(0:HM), CotAofH(0:HM),
     *     AREALL(HM), !  spherical area of Lon-Lat grid cells
C**** Geometric parameters for the Cube-Sphere grid
     *     XofI(0:IM), !  edge of Cube-Sphere cells on cube (-1:1)
     *  AREACS(IM,IM)  !  spherical area of Cube-Sphere grid cells
C**** Parameters for intersections of constant H lines (latitude)
      Integer*4 ::
     *        NMH(GM,HM), !  # of H segments from (G-1,H) to (G,H)
     *    IH(NMX1,GM,HM), !  I of cube-sphere cell south of segment
     *    JH(NMX1,GM,HM), !  J of cube-sphere cell south of segment
     *    KH(NMX1,GM,HM)  !  K of cube-sphere cell south of segment
      Real*16 ::
     *  OH(0:NMX1,GM,HM)  !  longitude of segemnt ends
C**** Parameters for intersections of constant J lines
      Integer*4 ::
     *        NMJ(IM,0:IM,6), !  # of J segments from (I-1,J) to (I,J)
     *    GJ(NMX1,IM,0:IM,6), !  G of lon-lat cell containing cell
     *    HJ(NMX1,IM,0:IM,6)  !  H of lon-lat cell containing cell
      Real*16 ::
     *  OJ(0:NMX1,IM,0:IM,6), !  longitude of segment ends
     *  AJ(0:NMX1,IM,0:IM,6), !  latitude of segment ends
     *  QJ(0:NMX1,IM,0:IM,6)  !  cube coordinate: QJ=[-1,1] as I=[0,IM]
C**** EGA = exchange grid area = intersection polygon of LL and CS grids
C**** EGAs grouped by location inside same Longitude-Latitude cell
      Integer*4 ::
     *       NMofGH(GM,HM), !  # of EGAs in lon-lat cell (G,H)
     *  IofNGH(NMX2,GM,HM), !  I of cube-sphere cell that contains EGA
     *  JofNGH(NMX2,GM,HM), !  J of cube-sphere cell that contains EGA
     *  KofNGH(NMX2,GM,HM)  !  K of cube-sphere cell that contains EGA
      Real*16 ::
     *  SofNGH(NMX2,GM,HM)  !  sperical area of each EGA
C**** EGAs grouped by location inside same Cube-Sphere cell
      Integer*4 ::
     *       NMofIJK(IM,IM,6), !  # of EGAs in cell (I,J,K)
     *  GofNIJK(NMX2,IM,IM,6), !  G of lon-lat cell that contains EGA
     *  HofNIJK(NMX2,IM,IM,6)  !  H of lon-lat cell that contains EGA
      Real*16 ::
     *  SofNIJK(NMX2,IM,IM,6)  !  spherical area of each EGA
C**** EGAs ungrouped
      Integer*4 ::
     *       NMofN, !  number of EGAs in cell (I,J,K)
     *  GofN(NMAX), !  G of Lon-Lat cell that contains N-th EGA
     *  HofN(NMAX), !  H of Lon-Lat cell that contains N-th EGA
     *  IofN(NMAX), !  I of Cube-Sphere cell that contains N-th EGA
     *  JofN(NMAX), !  J of Cube-Sphere cell that contains N-th EGA
     *  KofN(NMAX)  !  K of Cube-Sphere cell that contains N-th EGA
      Real*8 ::
     *  SofN(NMAX)  !  spherical area of N-th EGA
      EndModule LLvsCS
c*

      subroutine accurate_xgrid()

!@sum Longitude-Latitude grid cells intersect Cube-Sphere grid cells
!@+   on a sphere, their common areas being formed by polygons.
!@+   The set of such polygons is called the Exchange grid.
!@+   This program calculates the area of each Exchange grid cell.
!@+
!@+   Output of this program is the Exchange grid areas (Sof..) and
!@+   their associated numbered grid cells on the Longitude-Latitude
!@+   grid (G,H) or the Cube-Sphere grid (I,J,K).
!@+   Longitudinal numbering of grid cells is rather arbitrary and can
!@+   be modified by rotating the G's or the faces of the cube.
!@+   In the comments embedded in this program, it is assumed that
!@+   Face 1 of the cube projects onto longitudes between 45 W and 45 E
!@+   and that the western edges of cells with G = 1 coincide with the
!@+   International Date Line.  If Face 1 projects onto longitudes
!@+   between 135 E and 135 W and the western edges of cells with G = 1
!@+   coincide with the Greenwich Meridian, then the output arrays
!@+   (Sof.., G, H, I, J, K) are still correct and do not need to be
!@+   modified (although the embedded comments are not correct).

!@auth Gary Russell

      Use LLvsCS
      Implicit None
      Integer,External :: HofA,IofX
      Real*16,External :: AREA
      integer :: fid,vid,intshift
      include 'netcdf.inc'
      integer, allocatable, dimension(:,:) :: ijcub,ijlatlon
      integer, allocatable, dimension(:) :: itile
      real*8, allocatable, dimension(:) :: xarea
      real*8 :: areacheck
      character*200 :: ofi
      Integer :: G,H, I,J,K, N,M, GR,HR,G2,G3,G4, IR,JR, STATUS,IN
      Real*16 :: O,A, X,Y,Z, XSQ,YSQ, S,SGH,SIJK,SGLOB
C****
      TWOPIQ = 8 * ATan(1q0)
      Call GEOMLL
      Call GEOMCS
C****
C**** Define end points of H lines, compute OH, IH, JH and KH
C****
      Do 130 G=GM/2+1,5*GM/8  !  .1 to 45 E
      Do 120 H=HM/2,HM-1      !  0 to 89.9 N
       NMH(G,H) = 1
      OH(0,G,H) = OofG(G-1)
      OH(1,G,H) = OofG(G)
      Z = TanAofH(H) * SecOofG(G)  !  as measured from face 1
      If (Z > 1)  GoTo 110
C**** Eastern end of H line is inside face 1
      Y = TanOofG(G)
      IH(1,G,H) = IofX (Y)
      JH(1,G,H) = IofX (Z)
      KH(1,G,H) = 1
      GoTo 120
C**** Eastern end of H line is inside face 3
  110 X = CotAofH(H) * CosOofG(G)
      Y = CotAofH(H) * SinOofG(G)
      IH(1,G,H) = IofX (-X)
      JH(1,G,H) = IofX (-Y)
      KH(1,G,H) = 3
  120 Continue
C**** North pole H lines reside inside face 3
       NMH(G,HM) = 1
      OH(0,G,HM) = OofG(G-1)
      OH(1,G,HM) = OofG(G)
      IH(1,G,HM) = IM/2
      JH(1,G,HM) = IM/2
  130 KH(1,G,HM) = 3
C****
C**** Calculate additional segment intersections of H lines
C****
      Do 240 H=HM/2,HM-1  !  0 to 89.9 N
      A = AofH(H)
C**** I line intersections inside face 1
      Do 210 I=IM/2+1,IM-1  !  .1 to 44.9 E
      Y = XofI(I)
      Z = Sqrt(1 + Y**2) * TanAofH(H)
      If (Z > 1)  GoTo 210
      O = ATan2 (Y,1q0)
      G = (O/TWOPIQ + .5)*GM + 1
      J = IofX (Z)
      Call HLINE (G,H, I,J,1, O)
  210 Continue
C**** J line intersections inside face 1
      Do 220 J=IM/2+1,IM
      Z = XofI(J)
      If (TanAofH(H) > Z .or. Z*Sqrt(.5q0) > TanAofH(H))  GoTo 220
      O = ACos (TanAofH(H) / Z)  !  0 < O < 45 E
      Y = Tan (O)
      G = (O/TWOPIQ + .5)*GM + 1
      I = IofX (Y)
      Call HLINE (G,H, I,J,1, O)
  220 Continue
C**** I line intersections inside face 3
      Do 230 I=1,IM/2-1
      X = - XofI(I)
      YSQ = CotAofH(H)**2 - X**2
      If (YSQ > X**2-PRECIS)    GoTo 230  !  skip point if O > 45 E
      If (YSQ < PRECIS)         GoTo 230  !  skip X,Y,Z point (1,0,1)
      Y = Sqrt (YSQ)
      O = ACos (X * TanAofH(H))
      G = (O/TWOPIQ + .5)*GM + 1
      J = IofX (-Y)
      Call HLINE (G,H, I,J,3, O)
  230 Continue
C**** J line intersections inside face 3
      Do 240 J=1,IM/2-1
      Y = - XofI(J)
      XSQ = CotAofH(H)**2 - Y**2
      If (XSQ < Y**2+PRECIS)  GoTo 240  !  skip point if O > 45 E
      If (XSQ > 1)            GoTo 240  !  skipped point is in face 1
      X = Sqrt (XSQ)
      O = ASin (Y * TanAofH(H))
      G = (O/TWOPIQ + .5)*GM + 1
      I = IofX (-X)
      Call HLINE (G,H, I,J+1,3, O)
  240 Continue
C****
C**** Define initial J line intersections
C****
C**** Initialize J line intersection in face 1, quadrant 1 (0:45E)
      Do 310 J=IM/2+1,IM
      Z = XofI(J)
      Y = 0         !  I=IM/2
      O = 0         !  I=IM/2
      A = ATan (Z)  !  I=IM/2
      Do 310 I=IM/2+1,IM
      OJ(0,I,J,1) = O  !  from previous I
      AJ(0,I,J,1) = A  !  from previous I
      QJ(0,I,J,1) = Y  !  from previous I
      Y = XofI(I)
      O = ATan2 (Y,1q0)
      A = ATan (Z / Sqrt(1 + Y**2))
       NMJ(I,J,1) = 1
      GJ(1,I,J,1) = (O/TWOPIQ + .5-PRECIS)*GM + 1
      HJ(1,I,J,1) = HofA (A)
      OJ(1,I,J,1) = O
      AJ(1,I,J,1) = A
  310 QJ(1,I,J,1) = Y
C**** Initialize J line intersection in face 3, quadrant 1 (0:90E)
      Do 330 J=0,IM/2-1
      Y = - XofI(J)
      X = 1                          !  I=0
      O = ATan2 (Y,1q0)              !  I=0
      A = ATan (1 / Sqrt(1 + Y**2))  !  I=0
      Do 330 I=1,IM/2
      OJ(0,I,J,3) = O    !  from previous I
      AJ(0,I,J,3) = A    !  from previous I
      QJ(0,I,J,3) = - X  !  from previous I
      X = - XofI(I)
      O = ATan2 (Y,X)
      A = ATan (1 / Sqrt(X**2 + Y**2))
       NMJ(I,J,3) = 1
      GJ(1,I,J,3) = (O/TWOPIQ + .5-PRECIS)*GM + 1
      HJ(1,I,J,3) = HofA (A)
      OJ(1,I,J,3) = O
      AJ(1,I,J,3) = A
  330 QJ(1,I,J,3) = - X
C****
C**** Calculate additional J line intersections
C****
C**** Additional J line intersection in face 1, quadrant 1 (0:45E)
      Do 420 J=IM/2+1,IM
      Z = XofI(J)
      Do 410 G=GM/2+1,5*GM/8-1  !  G line intersetions
      Y = TanOofG(G)
      O = OofG(G)
      A = ATan (Z * CosOofG(G))
      I = IofX (Y)
      H = HofA (A)
  410 Call JLINE (I,J,1, G,H, O,A, Y)
      Do 420 H=HM/2+1,HM-1  !  H line intersections
      A = AofH(H)
      If (TanAofH(H) >= Z-PRECIS)  GoTo 420             !  A >= 45 N
      If (Z*Sqrt(.5q0)+PRECIS >= TanAofH(H))  GoTo 420  !  O >= 45 E
      O = ACos (TanAofH(H) / Z)  !  0 < O < 45 E
      Y = Tan (O)
      I = IofX (Y)
      G = (O/TWOPIQ + .5)*GM + 1
      Call JLINE (I,J,1, G,H+1, O,A, Y)
  420 Continue
C**** Additional J line intersection in face 3, quadrant 1 (0:90E)
      Do 440 J=0,IM/2-1
      Y = - XofI(J)
      Do 430 G=GM/2+1,3*GM/4-1
      If (G == 5*GM/8)  GoTo 430
      X = Y * CotOofG(G)
      If (X > 1)  GoTo 430
      O = OofG(G)
      A = ATan (SinOofG(G) / Y)
      I = IofX (-X)
      H = HofA (A)
      Call JLINE (I,J,3, G,H, O,A, -X)
  430 Continue
      Do 440 H=2*HM/3+1,HM-1
      A = AofH(H)
      XSQ = CotAofH(H)**2 - Y**2
      If (XSQ < PRECIS .or. XSQ > 1)  GoTo 440
      X = Sqrt (XSQ)
      O = ASin (Y * TanAofH(H))
      I = IofX (-X)
      G = (O/TWOPIQ + .5)*GM + 1
      Call JLINE (I,J,3, G,H, O,A, -X)
  440 Continue
C****
C**** Calculate Exchange grid areas
C****
      NMofGH(  GM/2+1:5*GM/8,HM/2+1:HM) = 0
      SofNGH(:,GM/2+1:5*GM/8,HM/2+1:HM) = 0
      Do 520 H=HM/2+1,HM      !  0 to 90 N
      Do 520 G=GM/2+1,5*GM/8  !  0 to 45 E
C**** Subtract areas that are south of H (latitude) lines
      Do 510 N=1,NMH(G,H-1)
      S = (OH(N,G,H-1) - OH(N-1,G,H-1)) * (SinAofH(H-1) + 1)
      I =  IH(N,G,H-1)
      J =  JH(N,G,H-1)  ;  If(H-1 == HM/2) J = J+1
      K =  KH(N,G,H-1)
      Call EXGRID (G,H, I,J,K, M)
  510 SofNGH(M,G,H) = SofNGH(M,G,H) - S
C**** Add areas that are north of H (latitude) lines
      Do 520 N=1,NMH(G,H)
      S = (OH(N,G,H) - OH(N-1,G,H)) * (SinAofH(H) + 1)
      I = IH(N,G,H)
      J = JH(N,G,H)
      K = KH(N,G,H)
      Call EXGRID (G,H, I,J,K, M)
  520 SofNGH(M,G,H) = SofNGH(M,G,H) + S
C**** Add and subtract areas from constant J lines for K=1
      Do 530 J=IM/2+1,IM
      Do 530 I=IM/2+1,IM  !  0 to 45 E
      Do 530 N=1,NMJ(I,J,1)
      S = AREA (OJ(N-1,I,J,1),OJ(N,I,J,1), AJ(N-1,I,J,1),AJ(N,I,J,1))
      G = GJ(N,I,J,1)
      H = HJ(N,I,J,1)
      Call EXGRID (G,H, I,J,1, M)
      SofNGH(M,G,H) = SofNGH(M,G,H) + S
      If (J==IM)  GoTo 530
      Call EXGRID (G,H, I,J+1,1, M)
      SofNGH(M,G,H) = SofNGH(M,G,H) - S
  530 Continue
C**** Add and subtract areas from constant J lines for K=3
      Do 540 J=0,IM/2-1
      Do 540 I=1,J  !  0 to 45 E
      Do 540 N=1,NMJ(I,J,3)
      S = AREA (OJ(N-1,I,J,3),OJ(N,I,J,3), AJ(N-1,I,J,3),AJ(N,I,J,3))
      G = GJ(N,I,J,3)
      H = HJ(N,I,J,3)
      Call EXGRID (G,H, I,J+1,3, M)
      SofNGH(M,G,H) = SofNGH(M,G,H) - S
      If (J==0)  GoTo 540
      Call EXGRID (G,H, I,J,3, M)
      SofNGH(M,G,H) = SofNGH(M,G,H) + S
  540 Continue
C**** Add and subtract areas from constant I lines for K=3
C**** Reflect J line about 45 E to obtain I lines
      Do 550 JR=0,IM/2-1
      I = JR
      Do 550 IR=JR+1,IM/2  !  45 to 90 E
      J = IR
      Do 550 N=1,NMJ(IR,JR,3)
      S = AREA (OJ(N-1,IR,JR,3),OJ(N,IR,JR,3),
     *          AJ(N-1,IR,JR,3),AJ(N,IR,JR,3))
      G = 5*GM/4 - GJ(N,IR,JR,3) + 1
      H = HJ(N,IR,JR,3)
      Call EXGRID (G,H, I+1,J,3, M)
      SofNGH(M,G,H) = SofNGH(M,G,H) - S
      If (I==0)  GoTo 550
      Call EXGRID (G,H, I,J,3, M)
      SofNGH(M,G,H) = SofNGH(M,G,H) + S
  550 Continue
C****
C**** Fill in remainder of GHEX common block by reflection and rotation
C****
C**** Reflect data about 0 longitude, (0:45E,0:90N) to (0:45W,0:90N)
      Do 610 H=HM/2,HM         !  0 to 90 N
      Do 610 GR=GM/2+1,5*GM/8  !  0 to 45 E
      G = GM - GR + 1          !  0 to 45 W
      NMofGH(G,H) = NMofGH(GR,H)
      Do 610 N=1,NMofGH(GR,H)
      SofNGH(N,G,H) = SofNGH(N,GR,H)
      If (KofNGH(N,GR,H) == 1)
     *  Then  ;  IofNGH(N,G,H) = IM - IofNGH(N,GR,H) + 1
                 JofNGH(N,G,H) = JofNGH(N,GR,H)
                 KofNGH(N,G,H) = 1
        Else  ;  IofNGH(N,G,H) = IofNGH(N,GR,H)
                 JofNGH(N,G,H) = IM - JofNGH(N,GR,H) + 1
                 KofNGH(N,G,H) = 3  ;  EndIf
  610 Continue
C**** Reflect data about equator, (45W:45E,0:90N) to (45W:45E,0:90S)
      Do 620 HR=HM/2+1,HM       !  0 to 90 N
      H = HM - HR + 1           !  0 to 90 S
      Do 620 G=3*GM/8+1,5*GM/8  !  45 W to 90 E
      NMofGH(G,H) = NMofGH(G,HR)
      Do 620 N=1,NMofGH(G,HR)
      SofNGH(N,G,H) = SofNGH(N,G,HR)
      If (KofNGH(N,G,HR) == 1)
     *  Then  ;  IofNGH(N,G,H) = IofNGH(N,G,HR)
                 JofNGH(N,G,H) = IM - JofNGH(N,G,HR) + 1
                 KofNGH(N,G,H) = 1
        Else  ;  IofNGH(N,G,H) = IM - JofNGH(N,G,HR) + 1
                 JofNGH(N,G,H) = IM - IofNGH(N,G,HR) + 1
                 KofNGH(N,G,H) = 6  ;  EndIf
  620 Continue
C**** Rotate data 90 degrees in longitude
      Do 630 H=1,HM             !  90 S to 90 N
      Do 630 G=3*GM/8+1,5*GM/8  !  45 W to 45 E
      G2 = G + GM/4
      G3 = G + GM/2  ;  If(G3 > GM) G3 = G3 - GM
      G4 = G - GM/4
      NMofGH(G2,H) = NMofGH(G,H)
      NMofGH(G3,H) = NMofGH(G,H)
      NMofGH(G4,H) = NMofGH(G,H)
      Do 630 N=1,NMofGH(G,H)
      SofNGH(N,G2,H) = SofNGH(N,G,H)
      SofNGH(N,G3,H) = SofNGH(N,G,H)
      SofNGH(N,G4,H) = SofNGH(N,G,H)
      SelectCase (KofNGH(N,G,H))
        Case (1)  ;  IofNGH(N,G2,H) = IofNGH(N,G,H)
                     JofNGH(N,G2,H) = JofNGH(N,G,H)
                     KofNGH(N,G2,H) = 2
                     IofNGH(N,G3,H) = IM - JofNGH(N,G,H) + 1
                     JofNGH(N,G3,H) = IofNGH(N,G,H)
                     KofNGH(N,G3,H) = 4
                     IofNGH(N,G4,H) = IM - JofNGH(N,G,H) + 1
                     JofNGH(N,G4,H) = IofNGH(N,G,H)
                     KofNGH(N,G4,H) = 5
        Case (3)  ;  IofNGH(N,G2,H) = IM - JofNGH(N,G,H) + 1
                     JofNGH(N,G2,H) = IofNGH(N,G,H)
                     KofNGH(N,G2,H) = 3
                     IofNGH(N,G3,H) = IM - IofNGH(N,G,H) + 1
                     JofNGH(N,G3,H) = IM - JofNGH(N,G,H) + 1
                     KofNGH(N,G3,H) = 3
                     IofNGH(N,G4,H) = JofNGH(N,G,H)
                     JofNGH(N,G4,H) = IM - IofNGH(N,G,H) + 1
                     KofNGH(N,G4,H) = 3
        Case (6)  ;  IofNGH(N,G2,H) = JofNGH(N,G,H)
                     JofNGH(N,G2,H) = IM - IofNGH(N,G,H) + 1
                     KofNGH(N,G2,H) = 6
                     IofNGH(N,G3,H) = IM - IofNGH(N,G,H) + 1
                     JofNGH(N,G3,H) = IM - JofNGH(N,G,H) + 1
                     KofNGH(N,G3,H) = 6
                     IofNGH(N,G4,H) = IM - JofNGH(N,G,H) + 1
                     JofNGH(N,G4,H) = IofNGH(N,G,H)
                     KofNGH(N,G4,H) = 6  ;  EndSelect
  630 Continue
C**** Check areas of Lon-Lat cells
      Do 640 H=1,HM
      Do 640 G=1,GM
      SGH = Sum (SofNGH(:NMofGH(G,H),G,H))
      If (Abs(SGH-AREALL(H)) < 1q-24)  GoTo 640
      Write (0,*) 'SGH   =',SGH,NMofGH(G,H),G,H
      Write (0,*) 'AREALL=',AREALL(H)
      Stop 640
  640 Continue
      Write (0,*) 'SofNGH is correct.'
C****
C**** Define variables on C-S grid from those on L-L grid
C****
      NMofIJK(  :,:,:) = 0
      SofNIJK(:,:,:,:) = 0
      Do 710 H=1,HM
      Do 710 G=1,GM
      Do 710 N=1,NMofGH(G,H)
      S = SofNGH(N,G,H)
      I = IofNGH(N,G,H)
      J = JofNGH(N,G,H)
      K = KofNGH(N,G,H)
      M = NMofIJK(I,J,K) + 1
      If (M > NMX2) Stop '710: NMofIJK > NMX2'
      NMofIJK(  I,J,K) = M
      GofNIJK(M,I,J,K) = G
      HofNIJK(M,I,J,K) = H
  710 SofNIJK(M,I,J,K) = S
C**** Check areas of Cube-Sphere cells
      Do 720 K=1,6
      Do 720 J=1,6
      Do 720 I=1,6
      SIJK = Sum (SofNIJK(:NMofIJK(I,J,K),I,J,K))
      If (Abs(SIJK-AREACS(I,J)) < 1q-24)  GoTo 720
      Write (0,*) 'SIJK  =',SIJK,NMofIJK(I,J,K),I,J,K
      Write (0,*) 'AREACS=',AREACS(I,J)
      Stop 720
  720 Continue
      Write (0,*) 'SofNIJK is correct.'
C****
C**** Define variables on the Exchange grid from those on L-L grid
C****
      M = 0
      SofN(:) = 0
      Do 810 H=1,HM
      Do 810 G=1,GM
      Do 810 N=1,NMofGH(G,H)
      M = M+1
      If (M > NMAX) Stop '810: NMofN > NMAX'
      GofN(M) = G
      HofN(M) = H
      IofN(M) = IofNGH(N,G,H)
      JofN(M) = JofNGH(N,G,H)
      KofN(M) = KofNGH(N,G,H)
  810 SofN(M) = SofNGH(N,G,H)
      NMofN   = M
      Write (0,*) 'NMofN =',NMofN
C**** Check global area of Excahnge grid cells
      SGLOB = Sum (SofN(1:NMofN))
      If (Abs(SGLOB-2*TWOPIQ) < 1q-12)  GoTo 820
      Write (0,*) 'SGLOB =',SGLOB
      Write (0,*) '4*pi  =',2*TWOPIQ
      Stop 820
  820 Write (0,*) 'SofN is correct.'

C***  Output data in netcdf file
      ofi='exexch.nc'
      
      status = nf_open(trim(ofi),nf_write,fid)
      if (status .ne. NF_NOERR) write(*,*) NF_STRERROR(status)

      allocate(itile(NMofN))
      allocate(ijcub(2,NMofN))
      allocate(ijlatlon(2,NMofN))
      allocate(xarea(NMofN))

      areacheck=0.0
      intshift=GM*170/360
      if (GM*170/360-intshift .ne. 0.0) then
         write(*,*) "STOP"
         stop
      else
         do in=1,NMofN
            itile(in)=KofN(in)
            
            ijcub(1,in)=IofN(in)
            ijcub(2,in)=JofN(in)
            
            ijlatlon(1,in)=GofN(in)+GM*170/360
            if (ijlatlon(1,in) .gt. GM) then
               ijlatlon(1,in)=ijlatlon(1,in)-GM
            endif
            ijlatlon(2,in)=HofN(in)
            
            xarea(in)=SofN(in)
c            write(*,*) "xarea=",SofN(in)
            
            areacheck=areacheck+SofN(in)
         enddo
      endif
      write(*,*) "areacheck=",areacheck/(2.*TWOPIQ)

      status = nf_inq_varid(fid,'tile1',vid)
      write(*,*) NF_STRERROR(status)
      status = nf_put_var_int(fid,vid,itile)

      status = nf_inq_varid(fid,'tile1_cell',vid)
      status = nf_put_var_int(fid,vid,ijcub)
      write(*,*) NF_STRERROR(status)

      status = nf_inq_varid(fid,'tile2_cell',vid)
      status = nf_put_var_int(fid,vid,ijlatlon)
      write(*,*) NF_STRERROR(status)

      status = nf_inq_varid(fid,'xgrid_area',vid)
      status = nf_put_var_double(fid,vid,xarea)
      write(*,*) NF_STRERROR(status)


      status = nf_close(fid)
      End


      Subroutine GEOMLL
C****
C**** Calculate GEOMetry of Longitude-Latitude grid cells
C****
      Use LLvsCS, Only: GM,HM, TWOPIQ, OofG,AofH,
     *    SinOofG,CosOofG,TanOofG, CscOofG,SecOofG,CotOofG,
     *    SinAofH,TanAofH,CotAofH, AREALL
      Implicit None
      Integer*4 :: G,H, HDLATD
C**** Functions of Longitude
      Do 10 G=0,GM
      OofG(G)    = - TWOPIQ/2 + G*TWOPIQ/GM
      SinOofG(G) = Sin (OofG(G))
      CosOofG(G) = Cos (OofG(G))
      TanOofG(G) = Tan (OofG(G))
      If (SinOofG(G) /= 0)  CscOofG(G) = 1 / SinOofG(G)
      If (CosOofG(G) /= 0)  SecOofG(G) = 1 / CosOofG(G)
   10 If (TanOofG(G) /= 0)  CotOofG(G) = 1 / TanOofG(G)
C**** Functions of Latitude
c**** up to 1x1 resolution. for .5 degrees use HDLATD = Nint(360./HM) / 2
      HDLATD = Nint (180./HM)  !  latitude spacing in degrees
      Do 20 H=1,HM-1
   20 AofH(H) = HDLATD*(H-HM/2)*TWOPIQ/360  !  unequal spacing at poles
      AofH(0) = - TWOPIQ / 4
      AofH(HM) =  TWOPIQ / 4
      Do 30 H=0,HM
      SinAofH(H) = Sin (AofH(H))
      TanAofH(H) = Tan (AofH(H))
   30 If (TanAofH(H) /= 0)  CotAofH(H) = 1 / TanAofH(H)
C**** Longitude-Latitude grid cell area
      Do 40 H=1,HM
   40 AREALL(H) = (SinAofH(H) - SinAofH(H-1)) * TWOPIQ/GM
      Return
      End


      Subroutine GEOMCS
C****
C**** Calculate GEOMetry of Cube-Sphere grid cells
C****
      Use LLvsCS, Only: IM, XofI,AREACS
      Implicit None
      Integer*4 :: I,J
      Real*16 ::
     *   ACOR, !  latitude of cube point (1,1,1) projected onto sphere
     *  dEDGE, !  equal spacing on sphere of projected cube edge cells
     *  X,Y,Z, !  fixed rectangular coordinates on cube surface
     *    A,B, !  projection from X,Y to surface of sphere
     *  W(0:IM,0:IM)  !  spherical angle projected from first quadrant
                      !  right angle on cube at vertex I,J
C**** Compute XofI(I)
      ACOR  = ATan(1/Sqrt(2q0))  !  lat of (1,1,1) projected onto sphere
      dEDGE = 2*ACOR / IM        !  spacing of cube edge on sphere
      Do 10 I=0,IM
   10 XofI(I) = Tan (- ACOR + I*dEDGE) * Sqrt(2q0)  !  range is [-1,1]
      XofI(0) = -1  ;  XofI(IM/2) = 0  ;  XofI(IM) = 1
C**** Compute angle projected from first quadrant right angle at I,J
      Do 20 J=0,IM
      Do 20 I=0,IM
      X = XofI(I)
      Y = XofI(J)
      Z = 1
      A = X / (Sqrt(X**2 + Y**2 + Z**2))
      B = Y / (Sqrt(X**2 + Y**2 + Z**2))
   20 W(I,J) = ACos (- A * B / Sqrt((1-A**2)*(1-B**2)))  !  S.W.Russell
C**** Area formula is sum of 4 angles - 2*pi =
C**** = W(I-1,J-1) + [pi - W(I,J-1)] + [pi - W(I-1,J)] + W(I,J) - 2*pi
      Do 30 I=1,IM
      Do 30 J=1,IM
   30 AREACS(I,J) = W(I-1,J-1) - W(I-1,J) - W(I,J-1) + W(I,J)
      Return
      End


      Function IofX (X)
C****
C**** IofX returns the value of I such that XofI(I-1) < X <= XofI(I)
C****
      Use LLvsCS, Only: IM, PRECIS, XofI
      Implicit None
      Integer*4 :: IofX, I
      Real*16   :: X
C****
      Do 10 I=0,IM
   10 If (X <= XofI(I)+PRECIS)  GoTo 20
   20 IofX = I
      Return
      End


      Function HofA (A)
C****
C**** HofA returns the value of H such that AofH(H-1) < A <= AofH(H)
C****
      Use LLvsCS, Only: HM, PRECIS, AofH
      Implicit None
      Integer*4 :: HofA, H
      Real*16   :: A
C****
      Do 10 H=0,HM
   10 If (A <= AofH(H)+PRECIS)  GoTo 20
   20 HofA = H
      Return
      End


      Subroutine HLINE (G,H, I,J,K, O)
C****
C**** Adds another intersection point to constant H (latitude) line
C****
      Use LLvsCS, Only: GM,HM,NMX1, PRECIS, NMH,IH,JH,KH,OH
      Implicit None
      Integer*4 :: G,H, I,J,K, NNEW,N
      Real*16   :: O
C****
      Do 10 NNEW=1,NMH(G,H)
   10 If (O < OH(NNEW,G,H)+PRECIS)  GoTo 20
      Stop 'HLINE 1'
   20 If (NMH(G,H)==NMX1)  Stop 'HLINE 20: NMH > NMX1'
      NMH(G,H) = NMH(G,H) + 1
      Do 30 N=NMH(G,H),NNEW+1,-1
      OH(N,G,H) = OH(N-1,G,H)
      IH(N,G,H) = IH(N-1,G,H)
      JH(N,G,H) = JH(N-1,G,H)
   30 KH(N,G,H) = KH(N-1,G,H)
      OH(NNEW,G,H) = O
      IH(NNEW,G,H) = I
      JH(NNEW,G,H) = J
      KH(NNEW,G,H) = K
      Return
      End


      Subroutine JLINE (I,J,K, G,H, O,A, Q)
C****
C**** Adds another intersection point to constant J lines
C****
      Use LLvsCS, Only: IM,NMX1, PRECIS, NMJ,GJ,HJ,OJ,AJ,QJ
      Implicit None
      Integer*4 :: I,J,K, G,H, NNEW,N
      Real*16   :: O,A, Q
C****
      Do 10 NNEW=1,NMJ(I,J,K)
   10 If (Q < QJ(NNEW,I,J,K)+PRECIS)  GoTo 20
      Stop 'JLINE 10'
   20 If (NMJ(I,J,K)==NMX1)  Stop 'JLINE 20: NMJ > NMX1'
      NMJ(I,J,K) = NMJ(I,J,K) + 1
      Do 30 N=NMJ(I,J,K),NNEW+1,-1
      GJ(N,I,J,K) = GJ(N-1,I,J,K)
      HJ(N,I,J,K) = HJ(N-1,I,J,K)
      OJ(N,I,J,K) = OJ(N-1,I,J,K)
      AJ(N,I,J,K) = AJ(N-1,I,J,K)
   30 QJ(N,I,J,K) = QJ(N-1,I,J,K)
      GJ(NNEW,I,J,K) = G
      HJ(NNEW,I,J,K) = H
      OJ(NNEW,I,J,K) = O
      AJ(NNEW,I,J,K) = A
      QJ(NNEW,I,J,K) = Q
      Return
      End


      Function AREA (O1,O2, A1,A2)
C****
C**** AREA calculates the area of the spherical triangle whose
C**** verticies are the south pole, (O1,A1) and (O2,A2).  S.W.Russell
C****
      Use LLvsCS, Only: TWOPIQ
      Implicit None
      Real*16 :: AREA, O1,O2, A1,A2,
     *  COSdO, SINA1,COSA1,SINA2,COSA2, COSB,COSG, ALPHA,BETA,GAMMA
C****
      COSdO = Cos (O2 - O1)
      SINA1 = Sin (A1)  ;  COSA1 = Cos (A1)
      SINA2 = Sin (A2)  ;  COSA2 = Cos (A2)
      COSB  = (COSdO*SINA1*COSA2 - COSA1*SINA2) /
     /        Sqrt(1 - (COSdO*COSA1*COSA2 + SINA1*SINA2)**2)
      COSG  = (COSdO*SINA2*COSA1 - COSA2*SINA1) /
     /        Sqrt(1 - (COSdO*COSA1*COSA2 + SINA1*SINA2)**2)
C****
      ALPHA = O2 - O1      !  south pole angle
      BETA  = ACos (COSB)  !  angle at (O1,A1)
      GAMMA = ACos (COSG)  !  angle at (O2,A2)
      AREA  = ALPHA + BETA + GAMMA - TWOPIQ/2
      Return
      End


      Subroutine EXGRID (G,H, I,J,K, M)
C****
C**** EXGRID organizes areas of the exchange grid
C**** Output: M = M-th exchange grid area for grid cell (G,H)
C****
      Use LLvsCS, Only: GM,HM,IM,NMX2,
     *                  NMofGH,IofNGH,JofNGH,KofNGH, SofNGH
      Implicit None
      Integer*4 :: G,H, I,J,K, M, N
C**** Locate M from 1 to NMofGH(G,H)
      Do 10 N=1,NMofGH(G,H)
   10 If (I==IofNGH(N,G,H) .and. J==JofNGH(N,G,H). and.
     *    K==KofNGH(N,G,H))  GoTo 20
C     N = NMofGH(G,H) + 1
      If (N > NMX2)  Stop 'EXGRID 10: NMofGH > NMX2'
      NMofGH(G,H)   = N
      IofNGH(N,G,H) = I
      JofNGH(N,G,H) = J
      KofNGH(N,G,H) = K
   20 M = N
      Return
      End
c*
