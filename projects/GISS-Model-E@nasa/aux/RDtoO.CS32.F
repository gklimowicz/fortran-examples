C**** RDtoO.CS32.F   Create River Direction from DtoO   2009/03/18
C****
      Integer*4,Parameter :: IM=32
      Real*4,   Parameter :: TWOPI=6.2831853, XMEANDER=1.1
      Integer*4 IRBID(IM,IM,6), ID(IM,IM,6),JD(IM,IM,6),KD(IM,IM,6)
      Real*4    FOCEAN(0:IM+1,0:IM+1,6),ZATMO(0:IM+1,0:IM+1,6),
     *            RBID(0:IM+1,0:IM+1,6), DtoO(0:IM+1,0:IM+1,6)
      Real*8    DtoOmin, CtoC
      Character*80  IFDIR, ZFILE,STNFILE,RDFILE, TITLE
C**** From Cube-Sphere (I,J,K) to Geocentric coordinates (A,B,C)
C**** Geocentric coordinates: (1,0,0) = (0E,0N), (0,0,1) = NorthPole
      Real*4 A(0:IM+1,0:IM+1,6),B(0:IM+1,0:IM+1,6),C(0:IM+1,0:IM+1,6)
C****
      ZFILE   = 'Z_CS32'
      STNFILE = 'STN_CS32'
      RDFILE  = 'RDtoO.CS32'
C**** Calculate haloed geocentric coordinates
      Call GEOMCS (IM,A,B,C)
C****
C**** Read in FOCEAN and ZATMO from Z file
C****
c      Call GetEnv ('IFDIR',IFDIR)
c      ZFILE = Trim(IFDIR) // '/' // ZFILE
      Open (11, File=ZFILE, Form='Unformatted', Status='Old')
      Read (11) TITLE,FOCEAN(1:IM,1:IM,:)
      Read (11)
      Read (11)
      Read (11)
      Read (11) TITLE,ZATMO(1:IM,1:IM,:)
      Close(11)
C**** Define halo cells at face edges for FOCEAN
      CORNER = 1
      Call HALOODD  (IM,FOCEAN,1,2,3,5,6,CORNER)
      Call HALOEVEN (IM,FOCEAN,2,4,3,1,6,CORNER)
      Call HALOODD  (IM,FOCEAN,3,4,5,1,2,CORNER)
      Call HALOEVEN (IM,FOCEAN,4,6,5,3,2,CORNER)
      Call HALOODD  (IM,FOCEAN,5,6,1,3,4,CORNER)
      Call HALOEVEN (IM,FOCEAN,6,2,1,5,4,CORNER)
C**** Define halo cells at face edges for ZATMO
      CORNER = 1e10
      Call HALOODD  (IM,ZATMO,3,4,5,1,2,CORNER)  !  includes Greenland
      Call HALOEVEN (IM,ZATMO,6,2,1,5,4,CORNER)  !  includes Antarctica
C****
C**** Read in River-Basin-ID and Distance-to-Ocean from STN-30p.bin
C****
      Open (1, File=STNFILE, Form='Unformatted', Status='Old')
      Read (1)
      Read (1)
      Read (1) TITLE,IRBID
      Read (1) TITLE,DtoO(1:IM,1:IM,:)
      Close(1)
C**** Define halo cells at face edges for RBID
      RBID(1:IM,1:IM,:) = IRBID(:,:,:)
      CORNER = 0
      Call HALOODD  (IM,RBID,1,2,3,5,6,CORNER)
      Call HALOEVEN (IM,RBID,2,4,3,1,6,CORNER)
      Call HALOODD  (IM,RBID,3,4,5,1,2,CORNER)
      Call HALOEVEN (IM,RBID,4,6,5,3,2,CORNER)
      Call HALOODD  (IM,RBID,5,6,1,3,4,CORNER)
      Call HALOEVEN (IM,RBID,6,2,1,5,4,CORNER)
C**** Define halo cells at face edges for DtoO
      CORNER = 1e10
      Call HALOODD  (IM,DtoO,1,2,3,5,6,CORNER)
      Call HALOEVEN (IM,DtoO,2,4,3,1,6,CORNER)
      Call HALOODD  (IM,DtoO,3,4,5,1,2,CORNER)
      Call HALOEVEN (IM,DtoO,4,6,5,3,2,CORNER)
      Call HALOODD  (IM,DtoO,5,6,1,3,4,CORNER)
      Call HALOEVEN (IM,DtoO,6,2,1,5,4,CORNER)
C****
C**** Determine downstream river cell from Distance-to-Ocean
C****
      Do 190 K=1,6
      Do 190 J=1,IM
      Do 190 I=1,IM
      ID(I,J,K) = I  ;  JD(I,J,K) = J  ;  KD(I,J,K) = K
      DtoOmin = 1d20
      If (FOCEAN(I,J,K) > 0)  GoTo 190
C**** Check direction 1
      If (RBID(I,J,K) == RBID(I+1,J+1,K) .and.
     *    DtoO(I,J,K) >  DtoO(I+1,J+1,K) .and.
     *    DtoOmin     >  DtoO(I+1,J+1,K))  Then
          DtoOmin     =  DtoO(I+1,J+1,K)
          ID(I,J,K) = I+1  ;  JD(I,J,K) = J+1  ;  GoTo 120  ;  EndIf
      CtoC = (I+.5 - J-.5)**2 - (I+.5 + J+.5 - IM-1)**2 + 1d10
      If (DtoOmin > CtoC .and. FOCEAN(I+1,J+1,K) > 0)  Then
          DtoOmin = CtoC
          ID(I,J,K) = I+1  ;  JD(I,J,K) = J+1  ;  EndIf
C**** Check direction 2
  120 If (RBID(I,J,K) == RBID(I,J+1,K) .and.
     *    DtoO(I,J,K) >  DtoO(I,J+1,K) .and.
     *    DtoOmin     >  DtoO(I,J+1,K))  Then
          DtoOmin     =  DtoO(I,J+1,K)
          ID(I,J,K) = I  ;  JD(I,J,K) = J+1  ;  GoTo 130  ;  EndIf
      CtoC = (2*J+1 - IM-1)**2 - (2*I - IM-1)**2 + 1d5
      If (DtoOmin > CtoC .and. FOCEAN(I,J+1,K) > 0)  Then
          DtoOmin = CtoC
          ID(I,J,K) = I  ;  JD(I,J,K) = J+1  ;  EndIf
C**** Check direction 3
  130 If (RBID(I,J,K) == RBID(I-1,J+1,K) .and.
     *    DtoO(I,J,K) >  DtoO(I-1,J+1,K) .and.
     *    DtoOmin     >  DtoO(I-1,J+1,K))  Then
          DtoOmin     =  DtoO(I-1,J+1,K)
          ID(I,J,K) = I-1  ;  JD(I,J,K) = J+1  ;  GoTo 140  ;  EndIf
      CtoC = (I-.5 + J+.5 - IM-1)**2 - (I-.5 - J-.5)**2 + 1d10
      If (DtoOmin > CtoC .and. FOCEAN(I-1,J+1,K) > 0)  Then
          DtoOmin = CtoC
          ID(I,J,K) = I-1  ;  JD(I,J,K) = J+1  ;  EndIf
C**** Check direction 4
  140 If (RBID(I,J,K) == RBID(I-1,J,K) .and.
     *    DtoO(I,J,K) >  DtoO(I-1,J,K) .and.
     *    DtoOmin     >  DtoO(I-1,J,K))  Then
          DtoOmin     =  DtoO(I-1,J,K)
          ID(I,J,K) = I-1  ;  JD(I,J,K) = J  ;  GoTo 150  ;  EndIf
      CtoC = (2*I-1 - IM-1)**2 - (2*J - IM-1)**2 + 1d5
      If (DtoOmin > CtoC .and. FOCEAN(I-1,J,K) > 0)  Then
          DtoOmin = CtoC
          ID(I,J,K) = I-1  ;  JD(I,J,K) = J  ;  EndIf
C**** Check direction 5
  150 If (RBID(I,J,K) == RBID(I-1,J-1,K) .and.
     *    DtoO(I,J,K) >  DtoO(I-1,J-1,K) .and.
     *    DtoOmin     >  DtoO(I-1,J-1,K))  Then
          DtoOmin     =  DtoO(I-1,J-1,K)
          ID(I,J,K) = I-1  ;  JD(I,J,K) = J-1  ;  GoTo 160  ;  EndIf
      CtoC = (I-.5 - J+.5)**2 - (I-.5 + J-.5 - IM-1)**2 + 1d10
      If (DtoOmin > CtoC .and. FOCEAN(I-1,J-1,K) > 0)  Then
          DtoOmin = CtoC
          ID(I,J,K) = I-1  ;  JD(I,J,K) = J-1  ;  EndIf
C**** Check direction 6
  160 If (RBID(I,J,K) == RBID(I,J-1,K) .and.
     *    DtoO(I,J,K) >  DtoO(I,J-1,K) .and.
     *    DtoOmin     >  DtoO(I,J-1,K))  Then
          DtoOmin     =  DtoO(I,J-1,K)
          ID(I,J,K) = I  ;  JD(I,J,K) = J-1  ;  GoTo 170  ;  EndIf
      CtoC = (2*J-1 - IM-1)**2 - (2*I - IM-1)**2 + 1d5
      If (DtoOmin > CtoC .and. FOCEAN(I,J-1,K) > 0)  Then
          DtoOmin = CtoC
          ID(I,J,K) = I  ;  JD(I,J,K) = J-1  ;  EndIf
C**** Check direction 7
  170 If (RBID(I,J,K) == RBID(I+1,J-1,K) .and.
     *    DtoO(I,J,K) >  DtoO(I+1,J-1,K) .and.
     *    DtoOmin     >  DtoO(I+1,J-1,K))  Then
          DtoOmin     =  DtoO(I+1,J-1,K)
          ID(I,J,K) = I+1  ;  JD(I,J,K) = J-1  ;  GoTo 180  ;  EndIf
      CtoC = (I+.5 + J-.5 - IM-1)**2 - (I+.5 - J+.5)**2 + 1d10
      If (DtoOmin > CtoC .and. FOCEAN(I+1,J-1,K) > 0)  Then
          DtoOmin = CtoC
          ID(I,J,K) = I+1  ;  JD(I,J,K) = J-1  ;  EndIf
C**** Check direction 8
  180 If (RBID(I,J,K) == RBID(I+1,J,K) .and.
     *    DtoO(I,J,K) >  DtoO(I+1,J,K) .and.
     *    DtoOmin     >  DtoO(I+1,J,K))  Then
          DtoOmin     =  DtoO(I+1,J,K)
          ID(I,J,K) = I+1  ;  JD(I,J,K) = J  ;  GoTo 190  ;  EndIf
      CtoC = (2*I+1 - IM-1)**2 - (2*J - IM-1)**2 + 1d5
      If (DtoOmin > CtoC .and. FOCEAN(I+1,J,K) > 0)  Then
          DtoOmin = CtoC
          ID(I,J,K) = I+1  ;  JD(I,J,K) = J  ;  EndIf
  190 Continue
C****
C**** Determine downstream river cell for Antarctica from ZATMO
C****
      zRT2 = 1 / Sqrt(2.)
      K=6
      Do 210 J=1,IM*4/5  !  Antarctica only, skip Tasmania
      Do 210 I=1,IM*4/5  !  Antarctica only, skip tip of South America
      ID(I,J,K) = I  ;  JD(I,J,K) = J  ;  KD(I,J,K) = K  ;  GRAD = -1
      If (FOCEAN(I,J,K) > 0)  GoTo 210
C**** Check direction 1
      DIST =  Sqrt ((A(I,J,K) - A(I+1,J+1,K))**2 +
     +              (B(I,J,K) - B(I+1,J+1,K))**2 +
     +              (C(I,J,K) - C(I+1,J+1,K))**2)
      GRADnew = (ZATMO(I,J,K) - ZATMO(I+1,J+1,K)) / DIST
      If (GRADnew >= GRAD)  Then
        ID(I,J,K) = I+1  ;  JD(I,J,K) = J+1  ;  GRAD = GRADnew  ;  EndIf
C**** Check direction 2
      DIST =  Sqrt ((A(I,J,K) - A(I,J+1,K))**2 +
     +              (B(I,J,K) - B(I,J+1,K))**2 +
     +              (C(I,J,K) - C(I,J+1,K))**2)
      GRADnew = (ZATMO(I,J,K) - ZATMO(I,J+1,K)) / DIST
      If (GRADnew >= GRAD)  Then
        ID(I,J,K) = I  ;  JD(I,J,K) = J+1  ;  GRAD = GRADnew  ;  EndIf
C**** Check direction 3
      DIST =  Sqrt ((A(I,J,K) - A(I-1,J+1,K))**2 +
     +              (B(I,J,K) - B(I-1,J+1,K))**2 +
     +              (C(I,J,K) - C(I-1,J+1,K))**2)
      GRADnew = (ZATMO(I,J,K) - ZATMO(I-1,J+1,K)) / DIST
      If (GRADnew >= GRAD)  Then
        ID(I,J,K) = I-1  ;  JD(I,J,K) = J+1  ;  GRAD = GRADnew  ;  EndIf
C**** Check direction 4
      DIST =  Sqrt ((A(I,J,K) - A(I-1,J,K))**2 +
     +              (B(I,J,K) - B(I-1,J,K))**2 +
     +              (C(I,J,K) - C(I-1,J,K))**2)
      GRADnew = (ZATMO(I,J,K) - ZATMO(I-1,J,K)) / DIST
      If (GRADnew >= GRAD)  Then
        ID(I,J,K) = I-1  ;  JD(I,J,K) = J  ;  GRAD = GRADnew  ;  EndIf
C**** Check direction 5
      DIST =  Sqrt ((A(I,J,K) - A(I-1,J-1,K))**2 +
     +              (B(I,J,K) - B(I-1,J-1,K))**2 +
     +              (C(I,J,K) - C(I-1,J-1,K))**2)
      GRADnew = (ZATMO(I,J,K) - ZATMO(I-1,J-1,K)) / DIST
      If (GRADnew >= GRAD)  Then
        ID(I,J,K) = I-1  ;  JD(I,J,K) = J-1  ;  GRAD = GRADnew  ;  EndIf
C**** Check direction 6
      DIST =  Sqrt ((A(I,J,K) - A(I,J-1,K))**2 +
     +              (B(I,J,K) - B(I,J-1,K))**2 +
     +              (C(I,J,K) - C(I,J-1,K))**2)
      GRADnew = (ZATMO(I,J,K) - ZATMO(I,J-1,K)) / DIST
      If (GRADnew >= GRAD)  Then
        ID(I,J,K) = I  ;  JD(I,J,K) = J-1  ;  GRAD = GRADnew  ;  EndIf
C**** Check direction 7
      DIST =  Sqrt ((A(I,J,K) - A(I+1,J-1,K))**2 +
     +              (B(I,J,K) - B(I+1,J-1,K))**2 +
     +              (C(I,J,K) - C(I+1,J-1,K))**2)
      GRADnew = (ZATMO(I,J,K) - ZATMO(I+1,J-1,K)) / DIST
      If (GRADnew >= GRAD)  Then
        ID(I,J,K) = I+1  ;  JD(I,J,K) = J-1  ;  GRAD = GRADnew  ;  EndIf
C**** Check direction 8
      DIST =  Sqrt ((A(I,J,K) - A(I+1,J,K))**2 +
     +              (B(I,J,K) - B(I+1,J,K))**2 +
     +              (C(I,J,K) - C(I+1,J,K))**2)
      GRADnew = (ZATMO(I,J,K) - ZATMO(I+1,J,K)) / DIST
      If (GRADnew >= GRAD)  Then
        ID(I,J,K) = I+1  ;  JD(I,J,K) = J  ;  GRAD = GRADnew  ;  EndIf
  210 Continue
C****
C**** Determine downstream river cell for Greenland from ZATMO
C****
      zRT2 = 1 / Sqrt(2.)
      K=3
      Do 310 J=Nint(IM*.30)+1  ,Nint(IM*.5)
      Do 310 I=Nint(IM*.97)+1-J,Nint(IM*1.15)
      ID(I,J,K) = I  ;  JD(I,J,K) = J  ;  KD(I,J,K) = K  ;  GRAD = -1
      If (FOCEAN(I,J,K) > 0)  GoTo 310
C**** Check direction 1
      DIST =  Sqrt ((A(I,J,K) - A(I+1,J+1,K))**2 +
     +              (B(I,J,K) - B(I+1,J+1,K))**2 +
     +              (C(I,J,K) - C(I+1,J+1,K))**2)
      GRADnew = (ZATMO(I,J,K) - ZATMO(I+1,J+1,K)) / DIST
      If (GRADnew >= GRAD)  Then
        ID(I,J,K) = I+1  ;  JD(I,J,K) = J+1  ;  GRAD = GRADnew  ;  EndIf
C**** Check direction 2
      DIST =  Sqrt ((A(I,J,K) - A(I,J+1,K))**2 +
     +              (B(I,J,K) - B(I,J+1,K))**2 +
     +              (C(I,J,K) - C(I,J+1,K))**2)
      GRADnew = (ZATMO(I,J,K) - ZATMO(I,J+1,K)) / DIST
      If (GRADnew >= GRAD)  Then
        ID(I,J,K) = I  ;  JD(I,J,K) = J+1  ;  GRAD = GRADnew  ;  EndIf
C**** Check direction 3
      DIST =  Sqrt ((A(I,J,K) - A(I-1,J+1,K))**2 +
     +              (B(I,J,K) - B(I-1,J+1,K))**2 +
     +              (C(I,J,K) - C(I-1,J+1,K))**2)
      GRADnew = (ZATMO(I,J,K) - ZATMO(I-1,J+1,K)) / DIST
      If (GRADnew >= GRAD)  Then
        ID(I,J,K) = I-1  ;  JD(I,J,K) = J+1  ;  GRAD = GRADnew  ;  EndIf
C**** Check direction 4
      DIST =  Sqrt ((A(I,J,K) - A(I-1,J,K))**2 +
     +              (B(I,J,K) - B(I-1,J,K))**2 +
     +              (C(I,J,K) - C(I-1,J,K))**2)
      GRADnew = (ZATMO(I,J,K) - ZATMO(I-1,J,K)) / DIST
      If (GRADnew >= GRAD)  Then
        ID(I,J,K) = I-1  ;  JD(I,J,K) = J  ;  GRAD = GRADnew  ;  EndIf
C**** Check direction 5
      DIST =  Sqrt ((A(I,J,K) - A(I-1,J-1,K))**2 +
     +              (B(I,J,K) - B(I-1,J-1,K))**2 +
     +              (C(I,J,K) - C(I-1,J-1,K))**2)
      GRADnew = (ZATMO(I,J,K) - ZATMO(I-1,J-1,K)) / DIST
      If (GRADnew >= GRAD)  Then
        ID(I,J,K) = I-1  ;  JD(I,J,K) = J-1  ;  GRAD = GRADnew  ;  EndIf
C**** Check direction 6
      DIST =  Sqrt ((A(I,J,K) - A(I,J-1,K))**2 +
     +              (B(I,J,K) - B(I,J-1,K))**2 +
     +              (C(I,J,K) - C(I,J-1,K))**2)
      GRADnew = (ZATMO(I,J,K) - ZATMO(I,J-1,K)) / DIST
      If (GRADnew >= GRAD)  Then
        ID(I,J,K) = I  ;  JD(I,J,K) = J-1  ;  GRAD = GRADnew  ;  EndIf
C**** Check direction 7
      DIST =  Sqrt ((A(I,J,K) - A(I+1,J-1,K))**2 +
     +              (B(I,J,K) - B(I+1,J-1,K))**2 +
     +              (C(I,J,K) - C(I+1,J-1,K))**2)
      GRADnew = (ZATMO(I,J,K) - ZATMO(I+1,J-1,K)) / DIST
      If (GRADnew >= GRAD)  Then
        ID(I,J,K) = I+1  ;  JD(I,J,K) = J-1  ;  GRAD = GRADnew  ;  EndIf
C**** Check direction 8
      DIST =  Sqrt ((A(I,J,K) - A(I+1,J,K))**2 +
     +              (B(I,J,K) - B(I+1,J,K))**2 +
     +              (C(I,J,K) - C(I+1,J,K))**2)
      GRADnew = (ZATMO(I,J,K) - ZATMO(I+1,J,K)) / DIST
      If (GRADnew >= GRAD)  Then
        ID(I,J,K) = I+1  ;  JD(I,J,K) = J  ;  GRAD = GRADnew  ;  EndIf
  310 Continue
C****
C**** Determine proper coordinates of river flow that leaves a face
C****
      Call RDODD  (IM,ID,JD,KD,1,2,3,5,6)
      Call RDEVEN (IM,ID,JD,KD,2,4,3,1,6)
      Call RDODD  (IM,ID,JD,KD,3,4,5,1,2)
      Call RDEVEN (IM,ID,JD,KD,4,6,5,3,2)
      Call RDODD  (IM,ID,JD,KD,5,6,1,3,4)
      Call RDEVEN (IM,ID,JD,KD,6,2,1,5,4)
C****
C**** Hand redefine downstream coordinates
C****
      ID(24, 2,3) = 25  ;  JD(24, 2,3) =  2  ;  KD(24, 2,3) = 3
      ID(19, 3,4) = 18  ;  JD(19, 3,4) =  3  ;  KD(19, 3,4) = 4
      ID( 5,21,5) =  5  ;  JD( 5,21,5) = 22  ;  KD( 5,21,5) = 5
      ID( 5,22,5) =  5  ;  JD( 5,22,5) = 23  ;  KD( 5,22,5) = 5
      ID( 6,22,5) =  7  ;  JD( 6,22,5) = 22  ;  KD( 6,22,5) = 5
      ID( 8,23,5) =  9  ;  JD( 8,23,5) = 23  ;  KD( 8,23,5) = 5
      ID( 3,27,5) =  2  ;  JD( 3,27,5) = 28  ;  KD( 3,27,5) = 5
C****
C**** Write out binary river direction file
C****
      Open (54, File=RDFILE, Form='Unformatted')
      TITLE = 'River Direction File from Distance-to-Ocean    IM=32'
      Write (54) TITLE,ID,JD,KD
      Close (54)
      Write (0,*) 'File written: ',Trim(RDFILE)
      Write (0,*) Trim(TITLE)
      End

      Subroutine GEOMCS (IM,A,B,C)
C****
C**** Calculate geocentric coordinates (A,B,C) from the center of each
C**** cube-sphere cell
C****
      Real*8,Parameter :: TWOPI = 6.283185307179586477d0,
     *                   RADIUS = 6371  !  (km)
      Real*4 A(0:IM+1,0:IM+1,6),B(0:IM+1,0:IM+1,6),C(0:IM+1,0:IM+1,6)
      Real*4 XofI(7200)
C**** Compute XofI(I)
      SQRT2 = Sqrt (2.)
      ACOR  = ATan (1 / SQRT2)  !  latitude of (1,1,1) projected onto s
      dEDGE = 2*ACOR/IM         !  spacing of cube edge
      Do 10 I=1,IM
   10 XofI(I) = Tan (- ACOR + (I-.5)*dEDGE) * SQRT2  !  range is (-1,1)
C**** Loop over cells on a Cube-Sphere face
      Do 20 J=1,IM
      Do 20 I=1,IM
C**** Compute A,B,C for face 1
      X = 1
      Y = XofI(I)
      Z = XofI(J)
      A(I,J,1) = RADIUS * X / Sqrt(X**2 + Y**2 + Z**2)
      B(I,J,1) = RADIUS * Y / Sqrt(X**2 + Y**2 + Z**2)
      C(I,J,1) = RADIUS * Z / Sqrt(X**2 + Y**2 + Z**2)
C**** Compute A,B,C for face 2
      X = - XofI(I)
      Y = 1
      Z = XofI(J)
      A(I,J,2) = RADIUS * X / Sqrt(X**2 + Y**2 + Z**2)
      B(I,J,2) = RADIUS * Y / Sqrt(X**2 + Y**2 + Z**2)
      C(I,J,2) = RADIUS * Z / Sqrt(X**2 + Y**2 + Z**2)
C**** Compute A,B,C for face 3
      X = - XofI(I)
      Y = - XofI(J)
      Z = 1
      A(I,J,3) = RADIUS * X / Sqrt(X**2 + Y**2 + Z**2)
      B(I,J,3) = RADIUS * Y / Sqrt(X**2 + Y**2 + Z**2)
      C(I,J,3) = RADIUS * Z / Sqrt(X**2 + Y**2 + Z**2)
C**** Compute A,B,C for face 4
      X = - 1
      Y = - XofI(J)
      Z = - XofI(I)
      A(I,J,4) = RADIUS * X / Sqrt(X**2 + Y**2 + Z**2)
      B(I,J,4) = RADIUS * Y / Sqrt(X**2 + Y**2 + Z**2)
      C(I,J,4) = RADIUS * Z / Sqrt(X**2 + Y**2 + Z**2)
C**** Compute A,B,C for face 5
      X = XofI(J)
      Y = - 1
      Z = - XofI(I)
      A(I,J,5) = RADIUS * X / Sqrt(X**2 + Y**2 + Z**2)
      B(I,J,5) = RADIUS * Y / Sqrt(X**2 + Y**2 + Z**2)
      C(I,J,5) = RADIUS * Z / Sqrt(X**2 + Y**2 + Z**2)
C**** Compute A,B,C for face 6
      X = XofI(J)
      Y = XofI(I)
      Z = - 1
      A(I,J,6) = RADIUS * X / Sqrt(X**2 + Y**2 + Z**2)
      B(I,J,6) = RADIUS * Y / Sqrt(X**2 + Y**2 + Z**2)
   20 C(I,J,6) = RADIUS * Z / Sqrt(X**2 + Y**2 + Z**2)
C**** Define halo cells at face edges
      CORNER = 1e10
      Call HALOODD  (IM,A,1,2,3,5,6,CORNER)
      Call HALOODD  (IM,B,1,2,3,5,6,CORNER)
      Call HALOODD  (IM,C,1,2,3,5,6,CORNER)
      Call HALOEVEN (IM,A,2,4,3,1,6,CORNER)
      Call HALOEVEN (IM,B,2,4,3,1,6,CORNER)
      Call HALOEVEN (IM,C,2,4,3,1,6,CORNER)
      Call HALOODD  (IM,A,3,4,5,1,2,CORNER)
      Call HALOODD  (IM,B,3,4,5,1,2,CORNER)
      Call HALOODD  (IM,C,3,4,5,1,2,CORNER)
      Call HALOEVEN (IM,A,4,6,5,3,2,CORNER)
      Call HALOEVEN (IM,B,4,6,5,3,2,CORNER)
      Call HALOEVEN (IM,C,4,6,5,3,2,CORNER)
      Call HALOODD  (IM,A,5,6,1,3,4,CORNER)
      Call HALOODD  (IM,B,5,6,1,3,4,CORNER)
      Call HALOODD  (IM,C,5,6,1,3,4,CORNER)
      Call HALOEVEN (IM,A,6,2,1,5,4,CORNER)
      Call HALOEVEN (IM,B,6,2,1,5,4,CORNER)
      Call HALOEVEN (IM,C,6,2,1,5,4,CORNER)
C****
      Return
      EndSubroutine GEOMCS

      Subroutine HALOODD (IM,Z,K,KR,KT,KL,KB,CORNER)
C****
C**** Define values of for halo cells at edges of odd faces
C****
      Real*4 Z(0:IM+1,0:IM+1,6)
C**** Right edge
      Do 10 J=1,IM
   10 Z(IM+1,J,K) = Z(1,J,KR)
C**** Top edge
      Do 20 I=1,IM
   20 Z(I,IM+1,K) = Z(1,IM+1-I,KT)
C**** Left edge
      Do 30 J=1,IM
   30 Z(0,J,K) = Z(IM+1-J,IM,KL)
C**** Bottom edge
      Do 40 I=1,IM
   40 Z(I,0,K) = Z(I,IM,KB)
C**** Undefined corner values
      Z(0,0   ,K) = CORNER  ;  Z(IM+1,0   ,K) = CORNER
      Z(0,IM+1,K) = CORNER  ;  Z(IM+1,IM+1,K) = CORNER
      Return
      End

      Subroutine HALOEVEN (IM,Z,K,KR,KT,KL,KB,CORNER)
C****
C**** Define values of for halo cells at edges of even faces
C****
      Real*4 Z(0:IM+1,0:IM+1,6)
C**** Right edge
      Do 10 J=1,IM
   10 Z(IM+1,J,K) = Z(IM+1-J,1,KR)
C**** Top edge
      Do 20 I=1,IM
   20 Z(I,IM+1,K) = Z(I,1,KT)
C**** Left edge
      Do 30 J=1,IM
   30 Z(0,J,K) = Z(IM,J,KL)
C**** Bottom edge
      Do 40 I=1,IM
   40 Z(I,0,K) = Z(IM,IM+1-I,KB)
C**** Undefined corner values
      Z(0,0   ,K) = CORNER  ;  Z(IM+1,0   ,K) = CORNER
      Z(0,IM+1,K) = CORNER  ;  Z(IM+1,IM+1,K) = CORNER
      Return
      End

      Subroutine RDODD (IM,ID,JD,KD,K,KR,KT,KL,KB)
C****
C**** Replace downstream halo values with proper coordinates
C**** for odd face
C****
      Integer*4 ID(IM,IM,6),JD(IM,IM,6),KD(IM,IM,6)
C**** Check whether downstream cell leaves face K at right edge
      Do 10 J=1,IM
      If (ID(IM,J,K) == IM+1)  Then
          ID(IM,J,K) = 1
          KD(IM,J,K) = KR  ;  EndIf
   10 Continue
C**** Check whether downstream cell leaves face K at top edge
      Do 20 I=1,IM
      If (JD(I,IM,K) == IM+1)  Then
          JD(I,IM,K) = IM+1 - ID(I,IM,K)
          ID(I,IM,K) = 1
          KD(I,IM,K) = KT  ;  EndIf
   20 Continue
C**** Check whether downstream cell leaves face K at left edge
      Do 30 J=1,IM
      If (ID(1,J,K) == 0)  Then
          ID(1,J,K) = IM+1 - JD(1,J,K)
          JD(1,J,K) = IM
          KD(1,J,K) = KL  ;  EndIf
   30 Continue
C**** Check whether downstream cell leaves face K at bottom edge
      Do 40 I=1,IM
      If (JD(I,1,K) == 0)  Then
          JD(I,1,K) = IM
          KD(I,1,K) = KB  ;  EndIf
   40 Continue
      Return
      End

      Subroutine RDEVEN (IM,ID,JD,KD,K,KR,KT,KL,KB)
C****
C**** Replace downstream halo values with proper coordinates
C**** for odd face
C****
      Integer*4 ID(IM,IM,6),JD(IM,IM,6),KD(IM,IM,6)
C**** Check whether downstream cell leaves face K at right edge
      Do 10 J=1,IM
      If (ID(IM,J,K) == IM+1)  Then
          ID(IM,J,K) = IM+1 - JD(IM,J,K)
          JD(IM,J,K) = 1
          KD(IM,J,K) = KR  ;  EndIf
   10 Continue
C**** Check whether downstream cell leaves face K at top edge
      Do 20 I=1,IM
      If (JD(I,IM,K) == IM+1)  Then
          JD(I,IM,K) = 1
          KD(I,IM,K) = KT  ;  EndIf
   20 Continue
C**** Check whether downstream cell leaves face K at left edge
      Do 30 J=1,IM
      If (ID(1,J,K) == 0)  Then
          ID(1,J,K) = IM
          KD(1,J,K) = KL  ;  EndIf
   30 Continue
C**** Check whether downstream cell leaves face K at bottom edge
      Do 40 I=1,IM
      If (JD(I,1,K) == 0)  Then
          JD(I,1,K) = IM+1 - ID(I,1,K)
          ID(I,1,K) = IM
          KD(I,1,K) = KB  ;  EndIf
   40 Continue
      Return
      End
