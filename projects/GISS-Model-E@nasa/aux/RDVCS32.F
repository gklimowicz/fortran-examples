C**** RDVCS32.F                                      2009/03/11
C****
C**** Generates a PostScript file to display the river flow for
C**** Cube-Sphere 32x32 horizontal resolution
C****
C**** Usage: RDVCS32 RiverDirectionAsciiFile
C****
      Implicit  None
      Integer*4,External  :: IArgC
      Integer*4,Parameter :: IM=32, NHATCH=8  !  cells between hatching
      Integer*4 I,J,K, IU,JU,KU, ID,JD,KD, NARGS, COUNT,
     &          IDOWNS(IM,IM,6), IMOUTH(0:IM+1,0:IM+1,6),
     &          JDOWNS(IM,IM,6), JMOUTH(0:IM+1,0:IM+1,6),
     &          KDOWNS(IM,IM,6), KMOUTH(0:IM+1,0:IM+1,6)
      Real*4,Parameter :: TWOPI = 6.2831853
      Real*4    FOCEAN(IM,IM,6), yb, xl, x, y, dBox, temp, MapHeight,
     &          MapWidth, dX,dY,THETA
      Character*80  ZFILE,PSFILE, IFDIR,FILEIN, TITLE
      Character*128 BUFFER
C****
      ZFILE  = 'Z_CS32'
      PSFILE = 'RDVCS32.PS'
      NARGS  = IArgC()
      If (NARGS == 0)  GoTo 800  !  usage message
C**** Read in FOCEAN data
c      Call GetEnv ('IFDIR',IFDIR)
c      FILEIN = Trim(IFDIR) // '/' // ZFILE
      FILEIN = ZFILE
      Open (11,Status='Old',File=FILEIN,Form='Unformatted',Err=801)
      Read (11) TITLE,FOCEAN
      Close(11)
C**** Read in River Direction file
      Call GetArg (1,FILEIN)
      Open (54,Status='Old',File=FILEIN,Form='Unformatted',Err=801)
      Read (54) TITLE,IDOWNS,JDOWNS,KDOWNS
      Close(54)
      Write (6,*) 'river dorection file: ',TITLE
C****
C**** Find river mouths
C****
      Do 130 KU=1,6
      Do 120 JU=1,IM
      Do 120 IU=1,IM
      IMOUTH(IU,JU,KU) = 0
      JMOUTH(IU,JU,KU) = 0
      KMOUTH(IU,JU,KU) = 0
      If (FOCEAN(IU,JU,KU) > 0)  GoTo 120
      COUNT = 0
      I = IU  ;  J = JU  ;  K = KU
  110 ID = IDOWNS(I,J,K)  ;  IMOUTH(IU,JU,KU) = ID
      JD = JDOWNS(I,J,K)  ;  JMOUTH(IU,JU,KU) = JD
      KD = KDOWNS(I,J,K)  ;  KMOUTH(IU,JU,KU) = KD
      If (FOCEAN(ID,JD,KD) > 0)  GoTo 120
      If (ID==I .and. JD==J .and. KD==K)  GoTo 120
      COUNT = COUNT + 1
      If (COUNT > 2*IM)  Then
        Write (6,911) IU,JU,KU
        GoTo 120  ;  EndIf
      I = ID  ;  J = JD  ;  K = KD
      GoTo 110
  120 Continue
  130 Continue
C**** Fill in halo edges for river mouths
      Call RMHODD  (IM,IMOUTH,JMOUTH,KMOUTH,1,2,3,5,6)
      Call RMHEVEN (IM,IMOUTH,JMOUTH,KMOUTH,2,4,3,1,6)
      Call RMHODD  (IM,IMOUTH,JMOUTH,KMOUTH,3,4,5,1,2)
      Call RMHEVEN (IM,IMOUTH,JMOUTH,KMOUTH,4,6,5,3,2)
      Call RMHODD  (IM,IMOUTH,JMOUTH,KMOUTH,5,6,1,3,4)
      Call RMHEVEN (IM,IMOUTH,JMOUTH,KMOUTH,6,2,1,5,4)
C****
C**** Redefine downstram cells at face edges to point to halo cells
C****
      Call RDHODD  (IM,IDOWNS,JDOWNS,KDOWNS,1,2,3,5,6)
      Call RDHEVEN (IM,IDOWNS,JDOWNS,KDOWNS,2,4,3,1,6)
      Call RDHODD  (IM,IDOWNS,JDOWNS,KDOWNS,3,4,5,1,2)
      Call RDHEVEN (IM,IDOWNS,JDOWNS,KDOWNS,4,6,5,3,2)
      Call RDHODD  (IM,IDOWNS,JDOWNS,KDOWNS,5,6,1,3,4)
      Call RDHEVEN (IM,IDOWNS,JDOWNS,KDOWNS,6,2,1,5,4)
C****
C**** Write PS file
C****
C**** Determine PostScript page parameters
      dBOX = 540.0/IM
      MapWidth  = IM*dBox
      MapHeight = IM*dBox
      xl = (612.0 - MapWidth)/2.0
      yb = (792.0 - MapHeight)/2.0
      temp = 0.5*dBox
C**** Open PostScript file
      Open (2, File=PSFILE, Form='Formatted')
C**** PostScript header data: set up rotated page
      Write (2,'("%!PS-Adobe-3.0    PostScript file from RDVCS32")')
      Write (2,'("%%Orientation: Portrait")')
      Write (2,'("/rl {rlineto} def" / "/s {stroke} def")')
      Write (2,'("/m {moveto} def" / "/l {lineto} def")')
      Write (2,'("/b {m  ",F8.2," 0 rl  0 ",F8.2," rl" / F8.2," 0 rl ",
     &          " closepath fill" / "} bind def")') dBox, dBox, -dBox
C**** Loop over cube faces
      Do 400 KU=1,6
      Write (TITLE(78:80),920) KU
      Write (2,'("gsave")')
      Write (2,'("/a { gsave newpath translate rotate " /
     &          F8.2," dup scale" /
     &          "-0.571 0.143 m  0 0.143 l" / "0 0.429 l 0.429 0 l" /
     &          "0 -0.429 l  0 -0.143 l" / "-0.571 -0.143 l" /
     &          "closepath fill grestore" / "} bind def")') temp
C**** Shade in grid cell with any ocean
      Write (2,'("0.95 setgray")')
      Do 210 J = 1, IM
      y = yb + (J - 1)*dBox
      Do 210 I = 1, IM
      If (FOCEAN(I,J,KU) > .5)  Then
         x = xl + (i - 1)*dBox
         Write (2,'(F8.2,F8.2," b")') x, y  ;  EndIf
  210 Continue
C**** Draw arrows
      Write (2,'("0 setgray")')
      Do 220 J = 1, IM
      Y = yb + (j - 0.5)*dBox
      Do 220 I = 1, IM
      If (FOCEAN(I,J,KU) > 0)  GoTo 220
      dX = IDOWNS(I,J,KU) - I
      dY = JDOWNS(I,J,KU) - J
      If (dX==0 .and. dY==0)  GoTo 220
      X = xl + (i - 0.5)*dBox
      THETA = ATan2(dY,dX) * 360 / TWOPI
      Write (2,'(F8.2,F8.2,F8.2," a")') theta,x,y
  220 Continue
C****
C**** Draw drainage basin separators
C****
      Write (2,'("1 setlinecap   1 setlinewidth")')
C**** Check separators between cells (I,J) and (I+1,J)
      Do 240 J=1,IM
      Y = YB + (J-.5)*dBox
      Do 240 I=0,IM
      X = XL + I*dBox
      If (IMOUTH(I,J,KU) == IMOUTH(I+1,J,KU) .and.
     *    JMOUTH(I,J,KU) == JMOUTH(I+1,J,KU) .and.
     *    KMOUTH(I,J,KU) == KMOUTH(I+1,J,KU))  GoTo 240
C**** Check lower half separator between cells (I,J) and (I+1,J)
      If (IMOUTH(I,J-1,KU) == IMOUTH(I+1,J,KU) .and.
     *    JMOUTH(I,J-1,KU) == JMOUTH(I+1,J,KU) .and.
     *    KMOUTH(I,J-1,KU) == KMOUTH(I+1,J,KU))  Then
          Write (2,923) X, Y, -dBOX/2, -dBOX/2
          GoTo 230  ;  EndIf
      If (IMOUTH(I,J,KU) == IMOUTH(I+1,J-1,KU) .and.
     *    JMOUTH(I,J,KU) == JMOUTH(I+1,J-1,KU) .and.
     *    KMOUTH(I,J,KU) == KMOUTH(I+1,J-1,KU))  Then
          Write (2,923) X, Y, dBOX/2, -dBOX/2
          GoTo 230  ;  EndIf
      Write (2,923) X, Y, 0, -dBOX/2
C**** Check upper half separator between cells (I,J) and (I+1,J)
  230 If (IMOUTH(I,J,KU) == IMOUTH(I+1,J+1,KU) .and.
     *    JMOUTH(I,J,KU) == JMOUTH(I+1,J+1,KU) .and.
     *    KMOUTH(I,J,KU) == KMOUTH(I+1,J+1,KU))  Then
          Write (2,923) X, Y, dBOX/2, dBOX/2
          GoTo 240  ;  EndIf
      If (IMOUTH(I,J+1,KU) == IMOUTH(I+1,J,KU) .and.
     *    JMOUTH(I,J+1,KU) == JMOUTH(I+1,J,KU) .and.
     *    KMOUTH(I,J+1,KU) == KMOUTH(I+1,J,KU))  Then
          Write (2,923) X, Y, -dBOX/2, dBOX/2
          GoTo 240  ;  EndIf
      Write (2,923) X, Y, 0, dBOX/2
  240 Continue
C**** Check separators between cells (I,J) and (I,J+1)
      Do 260 J=0,IM
      Y = YB + J*dBox
      Do 260 I=1,IM
      X = XL + (I-.5)*dBox
      If (IMOUTH(I,J,KU) == IMOUTH(I,J+1,KU) .and.
     *    JMOUTH(I,J,KU) == JMOUTH(I,J+1,KU) .and.
     *    KMOUTH(I,J,KU) == KMOUTH(I,J+1,KU))  GoTo 260
C**** Check left half separator between cells (I,J) and (I,J+1)
      If (IMOUTH(I-1,J,KU) == IMOUTH(I,J+1,KU) .and.
     *    JMOUTH(I-1,J,KU) == JMOUTH(I,J+1,KU) .and.
     *    KMOUTH(I-1,J,KU) == KMOUTH(I,J+1,KU))  GoTo 250
      If (IMOUTH(I,J,KU) == IMOUTH(I-1,J+1,KU) .and.
     *    JMOUTH(I,J,KU) == JMOUTH(I-1,J+1,KU) .and.
     *    KMOUTH(I,J,KU) == KMOUTH(I-1,J+1,KU))  GoTo 250
      Write (2,923) X, Y, -dBOX/2, 0
C**** Check right half separator between cells (I,J) and (I,J+1)
  250 If (IMOUTH(I,J,KU) == IMOUTH(I+1,J+1,KU) .and.
     *    JMOUTH(I,J,KU) == JMOUTH(I+1,J+1,KU) .and.
     *    KMOUTH(I,J,KU) == KMOUTH(I+1,J+1,KU))  GoTo 260
      If (IMOUTH(I+1,J,KU) == IMOUTH(I,J+1,KU) .and.
     *    JMOUTH(I+1,J,KU) == JMOUTH(I,J+1,KU) .and.
     *    KMOUTH(I+1,J,KU) == KMOUTH(I,J+1,KU))  GoTo 260
      Write (2,923) X, Y, dBOX/2, 0
  260 Continue
C****
C**** Draw cross hatch lines that show border
C****
      Write (2,'(".1 setlinewidth")')
C**** Bottom to top hatch lines
      Do 310 I=0,IM,NHATCH
      X = XL + I*dBOX
      Write (2,'(F8.2,F8.2," m")') X, YB
      Write (2,'("0 ",F8.2," rl")') MapHeight
  310 Write (2,'("stroke")')
C**** Left to right hatch lines
      Do 320 J=0,IM,NHATCH
      Y = YB + J*dBOX
      Write (2,'(F8.2,F8.2," m")') XL, Y
      Write (2,'(F8.2," 0 rl")') MapWidth
  320 Write (2,'("stroke")')
C**** Draw title
      Write (2,'("save")')
      Write (2,'("/Helvetica findfont 14 scalefont setfont")')
      Write (2,'("(",A,")" / "dup stringwidth pop 648 exch sub 2 div" /
     &          "684 m show")') TITLE
C**** Draw date
C     Write (2,'(F8.2," 36 m (2002/10/03) show")') xl
C**** Draw NASA/GISS logo
      Write (2,'("currentfont 0.7 scalefont setfont" / "(NASA/GISS)",
     &      " dup stringwidth pop ",F8.2," exch sub ",F8.2," m show")')
     &      xl+MapWidth, yb-12.6
      Write (2,'("restore")')
      Write (2,'("grestore")')
      Write (2,'("showpage")')
  400 Continue
      Close (2,status='keep')
C**** Send file to printer
c     Write (BUFFER,'("lp -d ",A," ",A)') Trim(PrinterName),Trim(PSfile)
c     Call System (BUFFER)
      Write (6,*) 'File written: ',Trim(PSFILE)
      GoTo 999
C****
  800 Write (6,*) 'Usage: RDVCS32 RiverDirectionFile    2009/03/11'
      GoTo 999
  801 Write (0,*) 'Error opening input file: ',FILEIN
      Stop 901
C 811 Write (0,*) 'Cell (IU,JU,KU) flows into an infinite loop'
C     Write (0,*) 'IU,JU,KU =',IU,JU,KU
C     Stop 811
C****
  911 Format ('Cell (',I2.2,',',I2.2,',',I2.2,
     *        ') flows into an infinite loop')
  920 Format ('K=',I1)
  923 Format (F8.2,F8.2,' m ',2F8.2,' rl  s')
  999 End

      Subroutine RMHODD (IM,IMOUTH,JMOUTH,KMOUTH,K,KR,KT,KL,KB)
C****
C**** Fill in river-mouth grid coordinates at halo edges
C**** In RMHODD, face K+1 is adjacent along (IM,J,K).
C****
      Integer*4 IMOUTH(0:IM+1,0:IM+1,6),JMOUTH(0:IM+1,0:IM+1,6),
     *          KMOUTH(0:IM+1,0:IM+1,6)
C**** Right edge
      Do 10 J=1,IM
      IMOUTH(IM+1,J,K) = IMOUTH(1,J,KR)
      JMOUTH(IM+1,J,K) = JMOUTH(1,J,KR)
   10 KMOUTH(IM+1,J,K) = KMOUTH(1,J,KR)
C**** Top edge
      Do 20 I=1,IM
      IMOUTH(I,IM+1,K) = IMOUTH(1,IM+1-I,KT)
      JMOUTH(I,IM+1,K) = JMOUTH(1,IM+1-I,KT)
   20 KMOUTH(I,IM+1,K) = KMOUTH(1,IM+1-I,KT)
C**** Left edge
      Do 30 J=1,IM
      IMOUTH(0,J,K) = IMOUTH(IM+1-J,IM,KL)
      JMOUTH(0,J,K) = JMOUTH(IM+1-J,IM,KL)
   30 KMOUTH(0,J,K) = KMOUTH(IM+1-J,IM,KL)
C**** Bottom edge
      Do 40 I=1,IM
      IMOUTH(I,0,K) = IMOUTH(I,IM,KB)
      JMOUTH(I,0,K) = JMOUTH(I,IM,KB)
   40 KMOUTH(I,0,K) = KMOUTH(I,IM,KB)
      Return
      End

      Subroutine RMHEVEN (IM,IMOUTH,JMOUTH,KMOUTH,K,KR,KT,KL,KB)
C****
C**** Fill in river-mouth grid coordinates at halo edges
C**** In RMHEVEN, face K+1 is adjacent along (I,IM,K).
C****
      Integer*4 IMOUTH(0:IM+1,0:IM+1,6),JMOUTH(0:IM+1,0:IM+1,6),
     *          KMOUTH(0:IM+1,0:IM+1,6)
C**** Right edge
      Do 10 J=1,IM
      IMOUTH(IM+1,J,K) = IMOUTH(IM+1-J,1,KR)
      JMOUTH(IM+1,J,K) = JMOUTH(IM+1-J,1,KR)
   10 KMOUTH(IM+1,J,K) = KMOUTH(IM+1-J,1,KR)
C**** Top edge
      Do 20 I=1,IM
      IMOUTH(I,IM+1,K) = IMOUTH(I,1,KT)
      JMOUTH(I,IM+1,K) = JMOUTH(I,1,KT)
   20 KMOUTH(I,IM+1,K) = KMOUTH(I,1,KT)
C**** Left edge
      Do 30 J=1,IM
      IMOUTH(0,J,K) = IMOUTH(IM,J,KL)
      JMOUTH(0,J,K) = JMOUTH(IM,J,KL)
   30 KMOUTH(0,J,K) = KMOUTH(IM,J,KL)
C**** Bottom edge
      Do 40 I=1,IM
      IMOUTH(I,0,K) = IMOUTH(IM,IM+1-I,KR)
      JMOUTH(I,0,K) = JMOUTH(IM,IM+1-I,KR)
   40 KMOUTH(I,0,K) = KMOUTH(IM,IM+1-I,KR)
      Return
      End

      Subroutine RDHODD (IM,IDOWNS,JDOWNS,KDOWNS,K,KR,KT,KL,KB)
C****
C**** Replace downstream grid coordinate with halo values.
C**** In RDHODD, face K+1 is adjacent along (IM,J,K).
C****
      Integer*4 IDOWNS(IM,IM,6),JDOWNS(IM,IM,6),KDOWNS(IM,IM,6)
C**** Check whether downstream cell leaves face K at right edge
      Do 10 J=1,IM
      If (KDOWNS(IM,J,K) == KR)  Then
          IDOWNS(IM,J,K) = IM + 1
          KDOWNS(IM,J,K) = K  ;  EndIf
   10 Continue
C**** Check whether downstream cell leaves face K at top edge
      Do 20 I=1,IM
      If (KDOWNS(I,IM,K) == KT)  Then
          IDOWNS(I,IM,K) = IM + 1 - JDOWNS(I,IM,K)
          JDOWNS(I,IM,K) = IM + 1
          KDOWNS(I,IM,K) = K  ;  EndIf
   20 Continue
C**** Check whether downstream cell leaves face K at left edge
      Do 30 J=1,IM
      If (KDOWNS(1,J,K) == KL)  Then
          JDOWNS(1,J,K) = IM + 1 - IDOWNS(1,J,K)
          IDOWNS(1,J,K) = 0
          KDOWNS(1,J,K) = K  ;  EndIf
   30 Continue
C**** Check whether downstream cell leaves face K at bottom edge
      Do 40 I=1,IM
      If (KDOWNS(I,1,K) == KB)  Then
          JDOWNS(I,1,K) = 0
          KDOWNS(I,1,K) = K  ;  EndIf
   40 Continue
      Return
      End

      Subroutine RDHEVEN (IM,IDOWNS,JDOWNS,KDOWNS,K,KR,KT,KL,KB)
C****
C**** Replace downstream grid coordinate with halo values.
C**** In RDHEVEN, face K+1 is adjacent along (I,IM,K).
C****
      Integer*4 IDOWNS(IM,IM,6),JDOWNS(IM,IM,6),KDOWNS(IM,IM,6)
C**** Check whether downstream cell leaves face K at right edge
      Do 10 J=1,IM
      If (KDOWNS(IM,J,K) == KR)  Then
          JDOWNS(IM,J,K) = IM + 1 - IDOWNS(IM,J,K)
          IDOWNS(IM,J,K) = IM + 1
          KDOWNS(IM,J,K) = K  ;  EndIf
   10 Continue
C**** Check whether downstream cell leaves face K at top edge
      Do 20 I=1,IM
      If (KDOWNS(I,IM,K) == KT)  Then
          JDOWNS(I,IM,K) = IM + 1
          KDOWNS(I,IM,K) = K  ;  EndIf
   20 Continue
C**** Check whether downstream cell leaves face K at left edge
      Do 30 J=1,IM
      If (KDOWNS(1,J,K) == KL)  Then
          IDOWNS(1,J,K) = 0
          KDOWNS(1,J,K) = K  ;  EndIf
   30 Continue
C**** Check whether downstream cell leaves face K at bottom edge
      Do 40 I=1,IM
      If (KDOWNS(I,1,K) == KB)  Then
          IDOWNS(I,1,K) = IM + 1 - JDOWNS(I,1,K)
          JDOWNS(I,1,K) = 0
          KDOWNS(I,1,K) = K  ;  EndIf
   40 Continue
      Return
      End
