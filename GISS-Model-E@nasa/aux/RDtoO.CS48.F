C**** RDtoO.CS48.F   Create River Direction from DtoO   2009/03/13
C****
      Integer*4,Parameter :: IM=48
      Real*4,   Parameter :: TWOPI=6.2831853, XMEANDER=1.1
      Integer*4 IRBID(IM,IM,6), ID(IM,IM,6),JD(IM,IM,6),KD(IM,IM,6)
      Real*4    FOCEAN(0:IM+1,0:IM+1,6),ZATMO(0:IM+1,0:IM+1,6),
     *            RBID(0:IM+1,0:IM+1,6), DtoO(0:IM+1,0:IM+1,6)
      Character*80  IFDIR, ZFILE,STNFILE,RDFILE, TITLE
C**** From Cube-Sphere (I,J,K) to Geocentric coordinates (A,B,C)
C**** Geocentric coordinates: (1,0,0) = (0E,0N), (0,0,1) = NorthPole
      Real*4 A(0:IM+1,0:IM+1,6),B(0:IM+1,0:IM+1,6),C(0:IM+1,0:IM+1,6)
C****
      ZFILE   = 'Z_CS48'
      STNFILE = 'STN_CS48'
      RDFILE  = 'RDtoO.CS48'
C**** Calculate haloed geocentric coordinates
      Call GEOMCS (IM,A,B,C)
C****
C**** Read in FOCEAN and ZATMO from Z file
C****
      Call GetEnv ('IFDIR',IFDIR)
      ZFILE = Trim(IFDIR) // '/' // ZFILE
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
      Do 10 K=1,6
      Do 10 J=1,IM
      Do 10 I=1,IM
      ID(I,J,K) = I  ;  JD(I,J,K) = J  ;  KD(I,J,K) = K
      DtoOmin = 1e10
      If (FOCEAN(I,J,K) > 0)  GoTo 10
C**** Check direction 1
      If (DtoO(I+1,J+1,K) < DtoOmin .and.
     *    DtoO(I+1,J+1,K) < DtoO(I,J,K) .and.
     *   (RBID(I+1,J+1,K) ==RBID(I,J,K) .or. FOCEAN(I+1,J+1,K) > 0))
     *  Then  ;  DtoOmin = DtoO(I+1,J+1,K)
                 ID(I,J,K) = I+1  ;  JD(I,J,K) = J+1  ;  EndIf
C**** Check direction 2
      If (DtoO(I,J+1,K) < DtoOmin .and.
     *    DtoO(I,J+1,K) < DtoO(I,J,K) .and.
     *   (RBID(I,J+1,K) ==RBID(I,J,K) .or. FOCEAN(I,J+1,K) > 0))
     *  Then  ;  DtoOmin = DtoO(I,J+1,K)
                 ID(I,J,K) = I  ;  JD(I,J,K) = J+1  ;  EndIf
C**** Check direction 3
      If (DtoO(I-1,J+1,K) < DtoOmin .and.
     *    DtoO(I-1,J+1,K) < DtoO(I,J,K) .and.
     *   (RBID(I-1,J+1,K) ==RBID(I,J,K) .or. FOCEAN(I-1,J+1,K) > 0))
     *  Then  ;  DtoOmin = DtoO(I-1,J+1,K)
                 ID(I,J,K) = I-1  ;  JD(I,J,K) = J+1  ;  EndIf
C**** Check direction 4
      If (DtoO(I-1,J,K) < DtoOmin .and.
     *    DtoO(I-1,J,K) < DtoO(I,J,K) .and.
     *   (RBID(I-1,J,K) ==RBID(I,J,K) .or. FOCEAN(I-1,J,K) > 0))
     *  Then  ;  DtoOmin = DtoO(I-1,J,K)
                 ID(I,J,K) = I-1  ;  JD(I,J,K) = J  ;  EndIf
C**** Check direction 5
      If (DtoO(I-1,J-1,K) < DtoOmin .and.
     *    DtoO(I-1,J-1,K) < DtoO(I,J,K) .and.
     *   (RBID(I-1,J-1,K) ==RBID(I,J,K) .or. FOCEAN(I-1,J-1,K) > 0))
     *  Then  ;  DtoOmin = DtoO(I-1,J-1,K)
                 ID(I,J,K) = I-1  ;  JD(I,J,K) = J-1  ;  EndIf
C**** Check direction 6
      If (DtoO(I,J-1,K) < DtoOmin .and.
     *    DtoO(I,J-1,K) < DtoO(I,J,K) .and.
     *   (RBID(I,J-1,K) ==RBID(I,J,K) .or. FOCEAN(I,J-1,K) > 0))
     *  Then  ;  DtoOmin = DtoO(I,J-1,K)
                 ID(I,J,K) = I  ;  JD(I,J,K) = J-1  ;  EndIf
C**** Check direction 7
      If (DtoO(I+1,J-1,K) < DtoOmin .and.
     *    DtoO(I+1,J-1,K) < DtoO(I,J,K) .and.
     *   (RBID(I+1,J-1,K) ==RBID(I,J,K) .or. FOCEAN(I+1,J-1,K) > 0))
     *  Then  ;  DtoOmin = DtoO(I+1,J-1,K)
                 ID(I,J,K) = I+1  ;  JD(I,J,K) = J-1  ;  EndIf
C**** Check direction 8
      If (DtoO(I+1,J,K) < DtoOmin .and.
     *    DtoO(I+1,J,K) < DtoO(I,J,K) .and.
     *   (RBID(I+1,J,K) ==RBID(I,J,K) .or. FOCEAN(I+1,J,K) > 0))
     *  Then  ;  DtoOmin = DtoO(I+1,J,K)
                 ID(I,J,K) = I+1  ;  JD(I,J,K) = J  ;  EndIf
   10 Continue
C****
C**** Determine downstream river cell for Antarctica from ZATMO
C****
      zRT2 = 1 / Sqrt(2.)
      K=6
      Do 20 J=1,IM*4/5  !  Antarctica only, skip Tasmania
      Do 20 I=1,IM*4/5  !  Antarctica only, skip tip of South America
      ID(I,J,K) = I  ;  JD(I,J,K) = J  ;  KD(I,J,K) = K  ;  GRAD = -1
      If (FOCEAN(I,J,K) > 0)  GoTo 20
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
   20 Continue
C****
C**** Determine downstream river cell for Greenland from ZATMO
C****
      zRT2 = 1 / Sqrt(2.)
      K=3
      Do 30 J=Nint(IM*.30)+1  ,Nint(IM*.5)
      Do 30 I=Nint(IM*.97)+1-J,Nint(IM*1.15)
      ID(I,J,K) = I  ;  JD(I,J,K) = J  ;  KD(I,J,K) = K  ;  GRAD = -1
      If (FOCEAN(I,J,K) > 0)  GoTo 30
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
   30 Continue
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
      ID(46, 2,2) = 47  ;  JD(46, 2,2) =  1  ;  KD(46, 2,2) = 2
      ID(45, 3,2) = 46  ;  JD(45, 3,2) =  2  ;  KD(45, 3,2) = 2
      ID(46, 3,2) = 46  ;  JD(46, 3,2) =  2  ;  KD(46, 3,2) = 2
      ID(46,26,2) = 46  ;  JD(46,26,2) = 27  ;  KD(46,26,2) = 2
      ID(24,38,2) = 25  ;  JD(24,38,2) = 38  ;  KD(24,38,2) = 2
      ID(23,42,2) = 23  ;  JD(23,42,2) = 41  ;  KD(23,42,2) = 2
      ID(16,11,3) = 16  ;  JD(16,11,3) = 12  ;  KD(16,11,3) = 3
      ID(16,12,3) = 16  ;  JD(16,12,3) = 13  ;  KD(16,12,3) = 3
      ID(15,13,3) = 15  ;  JD(15,13,3) = 14  ;  KD(15,13,3) = 3
      ID(16,13,3) = 16  ;  JD(16,13,3) = 14  ;  KD(16,13,3) = 3
      ID(15,14,3) = 15  ;  JD(15,14,3) = 15  ;  KD(15,14,3) = 3
      ID(19,14,3) = 18  ;  JD(19,14,3) = 14  ;  KD(19,14,3) = 3
      ID(15,15,3) = 15  ;  JD(15,15,3) = 16  ;  KD(15,15,3) = 3
      ID( 7,36,3) =  8  ;  JD( 7,36,3) = 36  ;  KD( 7,36,3) = 3
      ID( 6,37,3) =  7  ;  JD( 6,37,3) = 37  ;  KD( 6,37,3) = 3
      ID( 7,37,3) =  8  ;  JD( 7,37,3) = 37  ;  KD( 7,37,3) = 3
      ID( 5,38,3) =  6  ;  JD( 5,38,3) = 38  ;  KD( 5,38,3) = 3
      ID( 6,38,3) =  6  ;  JD( 6,38,3) = 37  ;  KD( 6,38,3) = 3
      ID( 5,39,3) =  6  ;  JD( 5,39,3) = 39  ;  KD( 5,39,3) = 3
      ID( 5,40,3) =  5  ;  JD( 5,40,3) = 39  ;  KD( 5,40,3) = 3
      ID( 5,41,3) =  5  ;  JD( 5,41,3) = 40  ;  KD( 5,41,3) = 3
      ID( 5,42,3) =  5  ;  JD( 5,42,3) = 41  ;  KD( 5,42,3) = 3
      ID( 5,43,3) =  5  ;  JD( 5,43,3) = 42  ;  KD( 5,43,3) = 3
      ID( 5,44,3) =  5  ;  JD( 5,44,3) = 43  ;  KD( 5,44,3) = 3
      ID(26, 2,4) = 25  ;  JD(26, 2,4) =  3  ;  KD(26, 2,4) = 4
      ID(25, 3,4) = 24  ;  JD(25, 3,4) =  4  ;  KD(25, 3,4) = 4
      ID(28, 4,4) = 27  ;  JD(28, 4,4) =  4  ;  KD(28, 4,4) = 4
      ID(27, 5,4) = 26  ;  JD(27, 5,4) =  5  ;  KD(27, 5,4) = 4
      ID(31,45,4) = 31  ;  JD(31,45,4) = 44  ;  KD(31,45,4) = 4
      ID(14, 1,5) = 14  ;  JD(14, 1,5) = 48  ;  KD(14, 1,5) = 4
      ID(11,17,5) = 12  ;  JD(11,17,5) = 17  ;  KD(11,17,5) = 5
      ID(12,32,5) = 13  ;  JD(12,32,5) = 32  ;  KD(12,32,5) = 5
      ID(13,32,5) = 14  ;  JD(13,32,5) = 32  ;  KD(13,32,5) = 5
      ID(17,36,5) = 17  ;  JD(17,36,5) = 37  ;  KD(17,36,5) = 5
      ID(15,37,5) = 16  ;  JD(15,37,5) = 37  ;  KD(15,37,5) = 5
      ID(16,37,5) = 17  ;  JD(16,37,5) = 37  ;  KD(16,37,5) = 5
      ID(17,37,5) = 18  ;  JD(17,37,5) = 37  ;  KD(17,37,5) = 5
      ID(17,38,5) = 17  ;  JD(17,38,5) = 37  ;  KD(17,38,5) = 5
      ID( 3,43,5) =  3  ;  JD( 3,43,5) = 44  ;  KD( 3,43,5) = 5
C****
C**** Write out binary river direction file
C****
      Open (54, File=RDFILE, Form='Unformatted')
      TITLE = 'River Direction File from Distance-to-Ocean    IM=48'
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
