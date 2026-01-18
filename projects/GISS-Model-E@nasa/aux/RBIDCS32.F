C**** RBIDCS32.F                                   2009/03/10
C****
C**** Generates a PostScript file to display the river basins
C**** for the Cube-Sphere 32x32 horizontal resolution
C****
      Implicit  None
      Integer*4,External  :: IArgC
      Integer*4,Parameter :: IM=32
      Integer*4 IDRB(0:IM+1,0:IM+1,6), I,J,KU
      Real*4,Parameter :: TWOPI = 6.2831853
      Real*4    FOCEAN(IM,IM,6), X,Y, XL,YB, dBox, temp, MapHeight,
     &          MapWidth
      Character*80  ZFILE,STNFILE,PSFILE, TITLE
      Character*128 BUFFER
C****
      ZFILE   = 'Z_CS32'
      STNFILE = 'STN_CS32'
      PSFILE  = 'RBIDCS32.PS'
C**** Read in FOCEAN data
      Call GetEnv ('IFDIR',TITLE)
      ZFILE = Trim(TITLE) // '/' // ZFILE
      Open (11, Status='Old', File=ZFILE, Form='Unformatted')
      Read (11) TITLE,FOCEAN
      Close(11)
C**** Read in River Basin ID file
      Open (1, Status='Old', File=STNFILE, Form='Unformatted')
      Read (1)
      Read (1)
      Read (1) TITLE,IDRB(1:IM,1:IM,:)
      Close(1)
      Write (0,*) 'River Basin ID file read: ',TITLE
C**** Fill in halo edges for river mouths
      Call IDHODD  (IM,IDRB,1,2,3,5,6)
      Call IDHEVEN (IM,IDRB,2,4,3,1,6)
      Call IDHODD  (IM,IDRB,3,4,5,1,2)
      Call IDHEVEN (IM,IDRB,4,6,5,3,2)
      Call IDHODD  (IM,IDRB,5,6,1,3,4)
      Call IDHEVEN (IM,IDRB,6,2,1,5,4)
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
      Write (2,'("%!PS-Adobe-3.0    PostScript file from RBIDCS32")')
      Write (2,'("%%Orientation: Portrait")')
      Write (2,'("/rl {rlineto} def" / "/s {stroke} def")')
      Write (2,'("/m {moveto} def" / "/l {lineto} def")')
      Write (2,'("/b {m  ",F8.2," 0 rl  0 ",F8.2," rl" / F8.2," 0 rl ",
     &          " closepath fill" / "} bind def")') dBox, dBox, -dBox
C**** Loop over cube faces
      TITLE = 'River Basins from ' // STNFILE
      Do 400 KU=1,6
      Write (TITLE(78:80),920) KU
      Write (2,'("gsave")')
      Write (2,'("/a { gsave newpath translate rotate " /
     &          F8.2," dup scale" /
     &          "-0.571 0.143 m  0 0.143 l" / "0 0.429 l 0.429 0 l" /
     &          "0 -0.429 l  0 -0.143 l" / "-0.571 -0.143 l" /
     &          "closepath fill grestore" / "} bind def")') temp
C**** Shade in grid cell with any ocean
      Write (2,'("0.9 setgray")')
      Do 210 J = 1, IM
      y = yb + (J - 1)*dBox
      Do 210 I = 1, IM
      If (FOCEAN(I,J,KU) > .5)  Then
         x = xl + (i - 1)*dBox
         Write (2,'(F8.2,F8.2," b")') x, y  ;  EndIf
  210 Continue
C****
C**** Draw drainage basin separators
C****
      Write (2,'("0 setgray")')
      Write (2,'("1 setlinecap" / "1 setlinewidth")')
C**** Check separators between cells (I,J) and (I+1,J)
      Do 240 J=1,IM
      Y = YB + (J-.5)*dBox
      Do 240 I=0,IM
      X = XL + I*dBox
      If (IDRB(I,J,KU) == IDRB(I+1,J,KU))  GoTo 240
C**** Check lower half separator between cells (I,J) and (I+1,J)
      If (IDRB(I,J-1,KU) == IDRB(I+1,J,KU))  Then
          Write (2,923) X, Y, -dBOX/2, -dBOX/2
          GoTo 230  ;  EndIf
      If (IDRB(I,J,KU) == IDRB(I+1,J-1,KU))  Then
          Write (2,923) X, Y, dBOX/2, -dBOX/2
          GoTo 230  ;  EndIf
      Write (2,923) X, Y, 0, -dBOX/2
C**** Check upper half separator between cells (I,J) and (I+1,J)
  230 If (IDRB(I,J,KU) == IDRB(I+1,J+1,KU))  Then
          Write (2,923) X, Y, dBOX/2, dBOX/2
          GoTo 240  ;  EndIf
      If (IDRB(I,J+1,KU) == IDRB(I+1,J,KU))  Then
          Write (2,923) X, Y, -dBOX/2, dBOX/2
          GoTo 240  ;  EndIf
      Write (2,923) X, Y, 0, dBOX/2
  240 Continue
C**** Check separators between cells (I,J) and (I,J+1)
      Do 260 J=0,IM
      Y = YB + J*dBox
      Do 260 I=1,IM
      X = XL + (I-.5)*dBox
      If (IDRB(I,J,KU) == IDRB(I,J+1,KU))  GoTo 260
C**** Check left half separator between cells (I,J) and (I,J+1)
      If (IDRB(I-1,J,KU) == IDRB(I  ,J+1,KU))  GoTo 250
      If (IDRB(I  ,J,KU) == IDRB(I-1,J+1,KU))  GoTo 250
      Write (2,923) X, Y, -dBOX/2, 0
C**** Check right half separator between cells (I,J) and (I,J+1)
  250 If (IDRB(I  ,J,KU) == IDRB(I+1,J+1,KU))  GoTo 260
      If (IDRB(I+1,J,KU) == IDRB(I  ,J+1,KU))  GoTo 260
      Write (2,923) X, Y, dBOX/2, 0
  260 Continue
C**** Draw border
      Write (2,'(".5 setlinewidth")')
      Write (2,'(F8.2,F8.2," m")') xl, yb
      Write (2,'("0 ",F8.2," rl" / F8.2," 0 rl" / " 0 ",F8.2," rl")')
     &        MapHeight, MapWidth, -MapHeight
      Write (2,'("closepath" / "stroke")')
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
      Write (0,*) 'File written: ',Trim(PSFILE)
      Write (0,*) 'TITLE: ',Trim(TITLE)
C**** Send file to printer
c     Write (BUFFER,'("lp -d ",A," ",A)') Trim(PrinterName),Trim(PSfile)
c     Call System (BUFFER)
      GoTo 999
C****
  920 Format ('K=',I1)
  923 Format (F8.2,F8.2,' m ',2F8.2,' rl  s')
  999 End

      Subroutine IDHODD (IM,IDRB,K,KR,KT,KL,KB)
C****
C**** Fill in river-basin ID grid coordinates at halo edges
C**** In IDHODD, face K+1 is adjacent along (IM,J,K).
C****
      Integer*4 IDRB(0:IM+1,0:IM+1,6)
C**** Right edge
      Do 10 J=1,IM
   10 IDRB(IM+1,J,K) = IDRB(1,J,KR)
C**** Top edge
      Do 20 I=1,IM
   20 IDRB(I,IM+1,K) = IDRB(1,IM+1-I,KT)
C**** Left edge
      Do 30 J=1,IM
   30 IDRB(0,J,K) = IDRB(IM+1-J,IM,KL)
C**** Bottom edge
      Do 40 I=1,IM
   40 IDRB(I,0,K) = IDRB(I,IM,KB)
      Return
      End

      Subroutine IDHEVEN (IM,IDRB,K,KR,KT,KL,KB)
C****
C**** Fill in river-basin ID grid coordinates at halo edges
C**** In IDHEVEN, face K+1 is adjacent along (I,IM,K).
C****
      Integer*4 IDRB(0:IM+1,0:IM+1,6)
C**** Right edge
      Do 10 J=1,IM
   10 IDRB(IM+1,J,K) = IDRB(IM+1-J,1,KR)
C**** Top edge
      Do 20 I=1,IM
   20 IDRB(I,IM+1,K) = IDRB(I,1,KT)
C**** Left edge
      Do 30 J=1,IM
   30 IDRB(0,J,K) = IDRB(IM,J,KL)
C**** Bottom edge
      Do 40 I=1,IM
   40 IDRB(I,0,K) = IDRB(IM,IM+1-I,KR)
      Return
      End
