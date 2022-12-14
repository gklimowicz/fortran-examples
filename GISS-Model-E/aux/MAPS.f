C****  MAPS.GCM     REAL*4 Line Printer Maps for Model III     3/17/93
C****
      SUBROUTINE MAP (IM,JM,IHOUR,TITLE,A)
C****
C**** MAP calls the appropriate resolution line printer map
C****
      PARAMETER (IMJMAX=144*90)
      REAL*4 WEIGHT(IMJMAX)
C****
      IF(IM.EQ.36)  GO TO 110
      IF(IM.EQ.72)  GO TO 120
      IF(IM.EQ.144) GO TO 130
      IF(IM.EQ.12)  GO TO 140
      WRITE (6,901) IM,JM
      STOP 5
C****
      ENTRY CMAP (JM,IHOUR,TITLE,A)
  110 IF(36*JM.GT.IMJMAX)  WRITE (6,902) 36,JM
      DO 111 I=1,36*JM
  111 WEIGHT(I) = 1.
      CALL CMAP1 (JM,IHOUR,TITLE,A,WEIGHT,1.,0.,26)
      RETURN
C****
      ENTRY MMAP (JM,IHOUR,TITLE,A)
  120 IF(72*JM.GT.IMJMAX)  WRITE (6,902) 72,JM
      DO 121 I=1,72*JM
  121 WEIGHT(I) = 1.
      CALL MMAP1 (JM,IHOUR,TITLE,A,WEIGHT,1.,0.,26)
      RETURN
C****
      ENTRY FMAP (JM,IHOUR,TITLE,A)
  130 IF(144*JM.GT.IMJMAX)  WRITE (6,902) 144,JM
      DO 131 I=1,144*JM
  131 WEIGHT(I) = 1.
      CALL FMAP1 (JM,IHOUR,TITLE,A,WEIGHT,1.,0.,26)
      RETURN
C****
      ENTRY WMAP (JM,IHOUR,TITLE,A)
  140 IF(12*JM.GT.IMJMAX)  WRITE (6,902) 12*JM
      DO 141 I=1,12*JM
  141 WEIGHT(I) = 1.
      CALL WMAP1 (JM,IHOUR,TITLE,A,WEIGHT,1.,0.,26)
      RETURN
C****
C**** THIS ENTRY USED TO CHANGE SCALE, DIF OR KDSN
C****
      ENTRY MAP0 (IM,JM,IHOUR,TITLE,A,SCALE,DIF,KDSN)
      IF(IM.EQ.36)  GO TO 210
      IF(IM.EQ.72)  GO TO 220
      IF(IM.EQ.144) GO TO 230
      IF(IM.EQ.12)  GO TO 240
      WRITE (6,901) IM,JM
      STOP 40
C****
      ENTRY CMAP0 (JM,IHOUR,TITLE,A,SCALE,DIF,KDSN)
  210 IF(36*JM.GT.IMJMAX)  WRITE (6,902) 36,JM
      DO 211 I=1,36*JM
  211 WEIGHT(I) = 1.
      CALL  CMAP1 (JM,IHOUR,TITLE,A,WEIGHT,SCALE,DIF,KDSN)
      RETURN
C****
      ENTRY MMAP0 (JM,IHOUR,TITLE,A,SCALE,DIF,KDSN)
  220 IF(72*JM.GT.IMJMAX)  WRITE (6,902) 72,JM
      DO 221 I=1,72*JM
  221 WEIGHT(I) = 1.
      CALL  MMAP1 (JM,IHOUR,TITLE,A,WEIGHT,SCALE,DIF,KDSN)
      RETURN
C****
      ENTRY FMAP0 (JM,IHOUR,TITLE,A,SCALE,DIF,KDSN)
  230 IF(144*JM.GT.IMJMAX)  WRITE (6,902) 144,JM
      DO 231 I=1,144*JM
  231 WEIGHT(I) = 1.
      CALL  FMAP1 (JM,IHOUR,TITLE,A,WEIGHT,SCALE,DIF,KDSN)
      RETURN
C****
      ENTRY WMAP0 (JM,IHOUR,TITLE,A,SCALE,DIF,KDSN)
  240 IF(12*JM.GT.IMJMAX)  WRITE (6,902) 12,JM
      DO 241 I=1,12*JM
  241 WEIGHT(I) = 1.
      CALL  WMAP1 (JM,IHOUR,TITLE,A,WEIGHT,SCALE,DIF,KDSN)
      RETURN
C****
  901 FORMAT ('0Inappropriate IM and JM calling MAP routines:',2I5)
  902 FORMAT ('0MAP routines received IM times JM that exceeds',
     *  ' the internal dimension of WT.  Input IM,JM =',2I5)
      END

      SUBROUTINE MAP1 (IM,JM,IHOUR,TITLE,A,WEIGHT,SCALE,DIF,KDSN)
C****
C**** MAP1 is used to change the weighting for latiitudinal means
C****
      IF(IM.EQ.36)  GO TO 10
      IF(IM.EQ.72)  GO TO 20
      IF(IM.EQ.144) GO TO 30
      IF(IM.EQ.12)  GO TO 40
      WRITE (6,901) IM,JM
      STOP 60
C****
   10 CALL CMAP1 (JM,IHOUR,TITLE,A,WEIGHT,SCALE,DIF,KDSN)
      RETURN
   20 CALL MMAP1 (JM,IHOUR,TITLE,A,WEIGHT,SCALE,DIF,KDSN)
      RETURN
   30 CALL FMAP1 (JM,IHOUR,TITLE,A,WEIGHT,SCALE,DIF,KDSN)
      RETURN
   40 CALL WMAP1 (JM,IHOUR,TITLE,A,WEIGHT,SCALE,DIF,KDSN)
      RETURN
C****
  901 FORMAT ('0Inappropriate IM and JM calling MAP routines:',2I5)
      END

      SUBROUTINE CMAP1 (JM,IHOUR,TITLE,A,WEIGHT,SCALE,DIF,KDSN)
C****
C**** CMAP1 prints on the line printer the contents of the array A
C**** which is dimensioned 36 by JM.  If JM = 24 and KDSN > 0, the
C**** output indicates the continents for the medium resolution.
C****
      use TimeConstants_mod, only: INT_HOURS_PER_DAY, INT_DAYS_PER_YEAR
      PARAMETER (IM=36)
      REAL*4 A(36,1),WEIGHT(36,1)
      CHARACTER*80 TITLE,TITLEI
      COMMON /MAPCOM/ ALAT(24),HSUM(2),GSUM,APOINT(36,24)
      COMMON /CMAPFX/ RLAT(24),DXYP(24),HWT(2),PLAND(36,24)
      CHARACTER*1 CHAR(35)
      CHARACTER*3 LAND(36),LINE(36),BLANK,PERIOD,MINUS,PLUS
      CHARACTER*5 XCHAR
C****
      DATA JMLAST/0/, KDLAST/0/
      DATA SKIP /-1.e20/
      DATA CHAR/'1','2','3','4','5','6','7','8','9','A','B','C',
     *          'D','E','F','G','H','I','J','K','L','M','N','O',
     *          'P','Q','R','S','T','U','V','W','X','Y','Z'/
      DATA BLANK/'   '/, PERIOD/' ..'/, MINUS/' --'/, PLUS/' ++'/
      save jmlast,kdlast,skip,char,BLANK,PERIOD,MINUS,PLUS,FJEQ
C****
      IF(JM.EQ.JMLAST)  GO TO 200
      JMLAST=JM
C****
C**** CALCULATE SPHERICAL GEOMETRY
C****
      TWOPI  = 6.2831853
      RADIUS= 1.
      DLON  = TWOPI/IM
      DLAT  = NINT(180./(JM-1)) * TWOPI/360.
      FJEQ  = .5*(1+JM)
      SINS  = -1.
      DO 120 J=1,JM
      SINN = SIN(DLAT*(J+.5-FJEQ))
      IF(J.EQ.JM)  SINN = 1.
      DXYP(J) = DLON*RADIUS*RADIUS*(SINN-SINS)
      RLAT(J) = DLAT*(J-FJEQ)*360/TWOPI
  120 SINS = SINN
      DXYP(1) = DXYP( 1)*IM
      DXYP(JM)= DXYP(JM)*IM
      RLAT(1) = -90.
      RLAT(JM) = 90.
C****
C**** Read in PLAND from disk if JM = 24 and KDSN > 0
C****
  200 IF(KDSN.LT.0)  GO TO 700
      IF(KDSN.EQ.KDLAST)  GO TO 300
      KDLAST=KDSN
      IF(JM.NE.24 .OR. KDSN.EQ.0)  GO TO 210
      READ  (KDSN) TITLEI,PLAND
      REWIND KDSN
  210 DO 220 I=1,IM
  220 LAND(I) = BLANK
C****
C**** WRITE TITLE AND LONGITUDES AT TOP OF PAGE
C****
  300 JYEAR = IHOUR/(INT_HOURS_PER_DAY*INT_DAYS_PER_YEAR)
      JDAY  = IHOUR/INT_HOURS_PER_DAY-INT_DAYS_PER_YEAR*JYEAR
      JHOUR = IHOUR-INT_HOURS_PER_DAY*(INT_DAYS_PER_YEAR*JYEAR+JDAY)
      JYEAR = JYEAR+1958
      IF(TITLE(2:2).NE.'&')  WRITE (6,930) IHOUR,TITLE,JYEAR,JDAY,JHOUR
      WRITE (6,931) (LON,LON=-180,150,30),(I,I=1,IM)
C****
C**** OUTSIDE J LOOP
C****
      DO 520 JHEMI=2,1,-1
      HWT(JHEMI)  = 0.
      HSUM(JHEMI) = 0.
      JMAX=JM-(2-JHEMI)*(FJEQ-.5)
      JMIN= 1+(JHEMI-1)*(FJEQ-.5)
      DO 520 J=JMAX,JMIN,-1
C**** Produce periods over land only if JM = 24 and KDSN > 0
      IF(JM.NE.24 .OR. KDSN.GT.0)  GO TO 330
      DO 320 I=1,IM
      LAND(I) = BLANK
      IF(PLAND(I,J).GE..5)  LAND(I) = PERIOD
  320 CONTINUE
  330 WRITE (6,933) LAND
C****
C**** OUTSIDE I LOOP
C****
      SUM = 0
      WT  = 0.
      DO 450 I=1,IM
      WT = WT + WEIGHT(I,J)
C**** If WEIGHT=0 or A=SKIP, blank out the grid point
      IF(WEIGHT(I,J).NE.0. .AND. A(I,J).NE.SKIP)  GO TO 410
      LINE(I) = BLANK
      GO TO 450
C****
  410 APOINT(I,J) = A(I,J)*SCALE-DIF
      SUM = SUM + APOINT(I,J)*WEIGHT(I,J)
C****
      IF(APOINT(I,J).GT.-359.5)  GO TO 420
      LINE(I) = MINUS
      GO TO 450
C****
  420 IF(APOINT(I,J).LT.3599.5)  GO TO 430
      LINE(I) = PLUS
      GO TO 450
C****
  430 WRITE (XCHAR,'(F5.0)') APOINT(I,J)
      LINE(I) = XCHAR(2:4)
C****
      IF(APOINT(I,J).GT.-99.5)  GO TO 440
      READ (XCHAR,'(I3)') KTENS
      LINE(I)(1:1) = MINUS(3:3)
      LINE(I)(2:2) = CHAR(-KTENS)
      GO TO 450
C****
  440 IF(APOINT(I,J).LT.999.5)  GO TO 450
      READ (XCHAR,'(I2)') KHUNDS
      LINE(I)(1:1) = CHAR(KHUNDS)
  450 CONTINUE
C****
C**** END OF OUTSIDE I LOOP, PRODUCE NUMBERS FOR ONE LATITUDE
C****
      IF(J.GT.1 .AND. J.LT.JM)  GO TO 510
      WT  = WEIGHT(1,J)
      SUM = 0.
      IF(WEIGHT(1,J).NE.0. .AND. A(1,J).NE.SKIP)
     *  SUM = (A(1,J)*SCALE-DIF)*WEIGHT(1,J)
  510 ALAT(J) = 0.
      IF(WT.GT.0.)  ALAT(J) = SUM/WT
      WRITE (6,951) RLAT(J),J,LINE,J,ALAT(J)
      HWT(JHEMI)  = HWT(JHEMI)  + WT*DXYP(J)
  520 HSUM(JHEMI) = HSUM(JHEMI) + SUM*DXYP(J)
C****
C**** END OF OUTSIDE J LOOP,PRODUCE BOTTOM LONGITUDES AND GLOBAL MEAN
C****
      GSUM = HSUM(1) + HSUM(2)
      GWT  = HWT(1)  + HWT(2)
      IF(HWT(1).NE.0.)  HSUM(1) = HSUM(1)/HWT(1)
      IF(HWT(2).NE.0.)  HSUM(2) = HSUM(2)/HWT(2)
      IF(GWT   .NE.0.)  GSUM    = GSUM   /GWT
      WRITE (6,961) (I,I=1,IM),HSUM(2),HSUM(1),
     *              (LON,LON=-180,150,30),GSUM
      IF(JM.NE.24)  WRITE (6,962)
      RETURN
C****
C**** CALCULATE GLOBAL AND LATITUDINAL AVERAGES BUT DO NOT PRODUCE
C**** A LINE PRINTER MAP
C****
  700 IF(JM.NE.24)  RETURN
      DO 720 JHEMI=1,2
      HWT(JHEMI)  = 0.
      HSUM(JHEMI) = 0.
      JMAX=JM-(2-JHEMI)*(FJEQ-.5)
      JMIN= 1+(JHEMI-1)*(FJEQ-.5)
      DO 720 J=JMIN,JMAX
      IMAX=IM
      IF(J.EQ.1 .OR. J.EQ.JM)  IMAX=1
      WT  = 0.
      SUM = 0.
      DO 710 I=1,IMAX
      WT  = WT  + WEIGHT(I,J)
      APOINT(I,J) = SKIP
      IF(A(I,J).EQ.SKIP)  GO TO 710
      APOINT(I,J) = A(I,J)*SCALE-DIF
      SUM = SUM + APOINT(I,J)*WEIGHT(I,J)
  710 CONTINUE
      ALAT(J) = 0.
      IF(WT.NE.0.)  ALAT(J) = SUM/WT
      HWT(JHEMI)  = HWT(JHEMI)  + WT*DXYP(J)
  720 HSUM(JHEMI) = HSUM(JHEMI) + SUM*DXYP(J)
      GSUM = HSUM(1) + HSUM(2)
      GWT  = HWT(1)  + HWT(2)
      IF(HWT(1).NE.0.)  HSUM(1) = HSUM(1)/HWT(1)
      IF(HWT(2).NE.0.)  HSUM(2) = HSUM(2)/HWT(2)
      IF(GWT   .NE.0.)  GSUM    = GSUM   /GWT
      RETURN
C****
  930 FORMAT ('1  IHOUR =',I8,A86,I10,'  DAY',I4,', HOUR',I3//)
  931 FORMAT (7X,12I9/'0  LAT       ',36I3,'        MEAN'/)
  933 FORMAT (13X,36A3)
  951 FORMAT (F6.1,I5,2X,36A3,I5,F7.1)
  961 FORMAT (/'0  LAT       ',36I3,'   NH',F7.1/124X,'SH',F7.1/
     *  7X,12I9,9X,'GL',F7.1)
  962 FORMAT (' If JM is not equal to 24, the hemispheric and',
     *  ' global averages may be incorrect.')
      END

      SUBROUTINE MMAP1 (JM,IHOUR,TITLE,A,WEIGHT,SCALE,DIF,KDSN)
C****
C**** MMAP1 prints on the line printer the contents of the array A
C**** which is dimensioned 72 by JM.  If JM = 46 and KDSN > 0, the
C**** output indicates the continents for the 4x5 resolution.
C****
      use TimeConstants_mod, only: INT_HOURS_PER_DAY, INT_DAYS_PER_YEAR
      PARAMETER (IM=72)
      REAL*4 A(IM,1),WEIGHT(IM,1)
      CHARACTER*80 TITLE,TITLEI
      COMMON /FIXDCB/ FOCEAN(IM,46)
      COMMON /MAPCOM/ ALAT(46),HSUM(2),GSUM,APOINT(IM,46)
      COMMON /MMAPFX/ RLAT(46),DXYP(46),WT(46),SUM(46),HWT(2)
      CHARACTER*1 CHAR(35)
      CHARACTER*3 LINE(IM),BLANK,MINUS,PLUS
      CHARACTER*4 RCOL(2)
      CHARACTER*5 XCHAR
      LOGICAL*4 QLAND
C****
      DATA JMLAST/0/
      DATA SKIP /-1.e20/
      DATA CHAR/'1','2','3','4','5','6','7','8','9','A','B','C',
     *          'D','E','F','G','H','I','J','K','L','M','N','O',
     *          'P','Q','R','S','T','U','V','W','X','Y','Z'/
      DATA BLANK/'   '/, MINUS/' --'/, PLUS/' ++'/
      DATA RCOL/'LAT ','MEAN'/
      save jmlast,kdlast,skip,char,BLANK,PERIOD,MINUS,PLUS,FJEQ,rcol
C****
      QLAND = JM.EQ.46 .AND. KDSN.GT.0
      IF(JM.EQ.JMLAST)  GO TO 200
      JMLAST = JM
C****
C**** CALCULATE SPHERICAL GEOMETRY
C****
      TWOPI = 6.2831853
      RADIUS= 1.
      DLON  = TWOPI/IM
      DLAT  = NINT(180./(JM-1)) * TWOPI/360.
      FJEQ  = .5*(1+JM)
      SINS  = -1.
      DO 120 J=1,JM
      SINN = SIN(DLAT*(J+.5-FJEQ))
      IF(J.EQ.JM) SINN = 1.
      DXYP(J) = DLON*RADIUS*RADIUS*(SINN-SINS)
      RLAT(J) = DLAT*(J-FJEQ)*360/TWOPI
  120 SINS = SINN
      DXYP(1) = DXYP( 1)*IM
      DXYP(JM)= DXYP(JM)*IM
      RLAT(1) = -90.
      RLAT(JM) = 90.
C****
C**** Determine date from IHOUR
C****
  200 IF(KDSN.LT.0)  GO TO 700
      JYEAR = IHOUR/(INT_HOURS_PER_DAY*INT_DAYS_PER_YEAR)
      JDAY  = IHOUR/INT_HOURS_PER_DAY-INT_DAYS_PER_YEAR*JYEAR
      JHOUR = IHOUR-INT_HOURS_PER_DAY*(INT_DAYS_PER_YEAR*JYEAR+JDAY)
      JYEAR = JYEAR+1958
      DO 230 J=1,JM
      WT(J)  = 0.
  230 SUM(J) = 0.
      DO 240 JHEMI=1,2
      HWT(JHEMI)  = 0.
  240 HSUM(JHEMI) = 0.
C****
C**** OUTSIDE LOOP OVER HEMISPHERES (WESTERN, EASTERN)
C****
      DO 610 IHEMI=1,2
      IMIN=1+(IHEMI-1)*IM/2
      IMAX=IHEMI*IM/2
C**** Write TITLE and longitudes at top of page
      IF(TITLE(2:2).NE.'&')  WRITE (6,930) IHOUR,TITLE,JYEAR,JDAY,JHOUR
      LONMIN = 180*(IHEMI-2)
      LONMAX = LONMIN+180
      WRITE (6,931) (LON,LON=LONMIN,LONMAX,15),
     *              (I,I=IMIN,IMAX),RCOL(IHEMI)
C****
C**** OUTSIDE J LOOP
C****
      DO 570 JHEMI=2,1,-1
      JMAX=JM-(2-JHEMI)*(FJEQ-.5)
      JMIN= 1+(JHEMI-1)*(FJEQ-.5)
      DO 570 J=JMAX,JMIN,-1
C****
C**** OUTSIDE I LOOP
C****
      DO 450 I=IMIN,IMAX
      WT(J) = WT(J) + WEIGHT(I,J)
C**** If WEIGHT=0 or A=SKIP, blank out the grid point
c     if(j.eq.0)write(0,*) j,jmin,jmax,jhemi,jm,fjeq ! leave in, opt err
      IF(WEIGHT(I,J).NE.0. .AND. A(I,J).NE.SKIP)  GO TO 410
      LINE(I) = BLANK
      GO TO 450
C****
  410 APOINT(I,J) = A(I,J)*SCALE-DIF
      SUM(J) = SUM(J) + APOINT(I,J)*WEIGHT(I,J)
C****
      IF(APOINT(I,J).GT.-359.5)  GO TO 420
      LINE(I) = MINUS
      GO TO 450
C****
  420 IF(APOINT(I,J).LT.3599.5)  GO TO 430
      LINE(I) = PLUS
      GO TO 450
C****
  430 WRITE (XCHAR,'(F5.0)')  APOINT(I,J)
      LINE(I) = XCHAR(2:4)
C****
      IF(APOINT(I,J).GT.-99.5)  GO TO 440
      READ (XCHAR,'(I3)') KTENS
      LINE(I)(1:1) = MINUS(3:3)
      LINE(I)(2:2) = CHAR(-KTENS)
      GO TO 450
C****
  440 IF(APOINT(I,J).LT.999.5)  GO TO 450
      READ (XCHAR,'(I2)') KHUNDS
      LINE(I)(1:1) = CHAR(KHUNDS)
  450 CONTINUE
C****
C**** END OF OUTSIDE I LOOP, PRODUCE NUMBERS FOR ONE LATITUDE
C****
      IF(.NOT.QLAND)  GO TO 520
      DO 510 I=IMIN,IMAX
      IF(FOCEAN(I,J).GE..5)  GO TO 510
      IF(LINE(I)(1:1).EQ.' ')  LINE(I)(1:1) = ':'
      IF(LINE(I)(2:2).EQ.' ')  LINE(I)(2:2) = ':'
  510 CONTINUE
  520 IF(IHEMI.EQ.2)  GO TO 530
      WRITE (6,951) RLAT(J),J,(LINE(I),I=IMIN,IMAX),J,RLAT(J)
      GO TO 570
  530 IF(J.NE.1 .AND. J.NE.JM)  GO TO 540
      WT(J)  = WEIGHT(1,J)
      SUM(J) = 0.
      IF(WEIGHT(1,J).NE.0. .AND. A(1,J).NE.SKIP)
     *  SUM(J) = (A(1,J)*SCALE-DIF)*WEIGHT(1,J)
  540 ALAT(J) = 0.
      IF(WT(J).GT.0.)  ALAT(J) = SUM(J)/WT(J)
      WRITE (6,951) RLAT(J),J,(LINE(I),I=IMIN,IMAX),J,ALAT(J)
      HWT(JHEMI)  = HWT(JHEMI)  + WT(J)*DXYP(J)
      HSUM(JHEMI) = HSUM(JHEMI) + SUM(J)*DXYP(J)
  570 CONTINUE
C****
C**** END OF OUTSIDE J LOOP,PRODUCE BOTTOM LONGITUDES AND GLOBAL MEAN
C****
      IF(IHEMI.EQ.1)
     *  WRITE (6,960) (I,I=IMIN,IMAX),(LON,LON=LONMIN,LONMAX,15)
  610 CONTINUE
      GSUM = HSUM(1) + HSUM(2)
      GWT  = HWT(1)  + HWT(2)
      IF(HWT(1).NE.0.)  HSUM(1) = HSUM(1)/HWT(1)
      IF(HWT(2).NE.0.)  HSUM(2) = HSUM(2)/HWT(2)
      IF(GWT   .NE.0.)  GSUM    = GSUM   /GWT
      LONMAX=LONMAX-15
      WRITE (6,961) (I,I=IMIN,IMAX),HSUM(2),HSUM(1),
     *              (LON,LON=LONMIN,LONMAX,15),GSUM
      IF(JM.NE.46)  WRITE (6,962)
      RETURN
C****
C**** CALCULATE GLOBAL AND LATITUDINAL AVERAGES BUT DO NOT PRODUCE
C**** A LINE PRINTER MAP
C****
  700 IF(JM.NE.46)  RETURN
      DO 720 JHEMI=1,2
      HWT(JHEMI)  = 0.
      HSUM(JHEMI) = 0.
      JMAX=JM-(2-JHEMI)*(FJEQ-.5)
      JMIN= 1+(JHEMI-1)*(FJEQ-.5)
      DO 720 J=JMIN,JMAX
      IMAX=IM
      IF(J.EQ.1 .OR. J.EQ.JM)  IMAX=1
      WTJ  = 0.
      SUMJ = 0.
      DO 710 I=1,IMAX
      WTJ  = WTJ  + WEIGHT(I,J)
      APOINT(I,J) = SKIP
      IF(A(I,J).EQ.SKIP)  GO TO 710
      APOINT(I,J) = A(I,J)*SCALE-DIF
      SUMJ = SUMJ + APOINT(I,J)*WEIGHT(I,J)
  710 CONTINUE
      ALAT(J) = 0.
      IF(WTJ.NE.0.)  ALAT(J) = SUMJ/WTJ
      HWT(JHEMI)  = HWT(JHEMI)  + WTJ*DXYP(J)
  720 HSUM(JHEMI) = HSUM(JHEMI) + SUMJ*DXYP(J)
      GSUM = HSUM(1) + HSUM(2)
      GWT  = HWT(1)  + HWT(2)
      IF(HWT(1).NE.0.)  HSUM(1) = HSUM(1)/HWT(1)
      IF(HWT(2).NE.0.)  HSUM(2) = HSUM(2)/HWT(2)
      IF(GWT   .NE.0.)  GSUM    = GSUM   /GWT
      RETURN
C****
  930 FORMAT ('1  IHOUR =',I8,A86,I10,'  DAY',I4,', HOUR',I3//)
  931 FORMAT (6X,13I9/'0  LAT       ',36I3,8X,A4/)
  933 FORMAT ('+',12X,36A3)
  951 FORMAT (F6.1,I5,2X,36A3,I5,F7.1)
  960 FORMAT (/'0  LAT       ',36I3//6X,13I9)
  961 FORMAT (/'0  LAT       ',36I3,'   NH',F7.1/124X,'SH',F7.1/
     *  6X,12I9,10X,'GL',F7.1)
  962 FORMAT (' If JM is not equal to 46, the hemispheric and',
     *  ' global averages may be incorrect.')
      END

      SUBROUTINE FMAP1 (JM,IHOUR,TITLE,A,WEIGHT,SCALE,DIF,KDSN)
C****
C**** FMAP1 prints on the line printer the contents of the array A
C**** which is dimensioned 144 by JM.  If JM = 90 and KDSN > 0, the
C**** output indicates the continents for the new fine resolution.
C****
      use TimeConstants_mod, only: INT_HOURS_PER_DAY, INT_DAYS_PER_YEAR
      PARAMETER (IM=144)
      REAL*4 A(IM,1),WEIGHT(IM,1)
      CHARACTER*80 TITLE,TITLEI
      COMMON /FIXDCB/ FOCEAN(IM,90)
      COMMON /MAPCOM/ ALAT(90),HSUM(2),GSUM,APOINT(IM,90)
      COMMON /FMAPFX/ RLAT(180),DXYP(180),WT(180),SUM(180),HWT(2)
      CHARACTER*1 CHAR(35)
      CHARACTER*3 LINE(IM),BLANK,MINUS,PLUS
      CHARACTER*4 RCOL(4)
      CHARACTER*5 XCHAR
      LOGICAL*4 QLAND
C****
      DATA JMLAST/0/
      DATA SKIP /-1.e20/
      DATA CHAR/'1','2','3','4','5','6','7','8','9','A','B','C',
     *          'D','E','F','G','H','I','J','K','L','M','N','O',
     *          'P','Q','R','S','T','U','V','W','X','Y','Z'/
      DATA BLANK/'   '/, MINUS/' --'/, PLUS/' ++'/
      DATA RCOL/3*'LAT ','MEAN'/
      save jmlast,kdlast,skip,char,BLANK,PERIOD,MINUS,PLUS,FJEQ,rcol
C****
      QLAND = JM.EQ.90 .AND. KDSN.GT.0
      IF(JM.EQ.JMLAST)  GO TO 200
      JMLAST = JM
C****
C**** CALCULATE SPHERICAL GEOMETRY
C****
      TWOPI  = 6.283185307
      RADIUS = 1.
      DLON = TWOPI/IM
      DLAT = NINT(180./(JM-1)) * TWOPI/360.
      FJEQ = .5*(1+JM)
      SINS = -1.
      DO 120 J=1,JM
      SINN = SIN(DLAT*(J+.5-FJEQ))
      IF(J.EQ.JM) SINN = 1.
      DXYP(J) = DLON*RADIUS*RADIUS*(SINN-SINS)
      RLAT(J) = DLAT*(J-FJEQ)*360/TWOPI
  120 SINS=SINN
      RLAT(1) = -90.
      RLAT(JM) = 90.
C****
C**** Determine date from IHOUR
C****
  200 IF(KDSN.LT.0)  GO TO 700
      JYEAR = IHOUR/(INT_HOURS_PER_DAY*INT_DAYS_PER_YEAR)
      JDAY  = IHOUR/INT_HOURS_PER_DAY-INT_DAYS_PER_YEAR*JYEAR
      JHOUR = IHOUR-INT_HOURS_PER_DAY*(INT_DAYS_PER_YEAR*JYEAR+JDAY)
      JYEAR = JYEAR+1958
      DO 230 J=1,JM
      WT(J)  = 0.
  230 SUM(J) = 0.
      DO 240 JHEMI=1,2
      HWT(JHEMI)  = 0.
  240 HSUM(JHEMI) = 0.
C****
C**** Outside loop over longitudinal quarter spheres and
C**** latitudinal hemispheres
C****
      DO 610 IHEMI=1,4
      IMIN=1+(IHEMI-1)*IM/4
      IMAX=IHEMI*IM/4
      DO 610 JHEMI=2,1,-1
      JMAX=JM-(2-JHEMI)*(FJEQ-.5)
      JMIN= 1+(JHEMI-1)*(FJEQ-.5)
C**** Write TITLE and longitudes at top of page
      IF(TITLE(2:2).NE.'&')  WRITE (6,930) IHOUR,TITLE,JYEAR,JDAY,JHOUR
      LONMIN = 90*(IHEMI-3)
      LONMAX = LONMIN+90
      WRITE (6,931) (LON,LON=LONMIN,LONMAX,10),
     *              (I,I=IMIN,IMAX),RCOL(IHEMI)
C****
C**** J loop within a latitudinal hemisphere
C****
      DO 570 J=JMAX,JMIN,-1
C****
C**** I loop within a longitudinal quarter sphere
C****
      DO 450 I=IMIN,IMAX
      WT(J) = WT(J) + WEIGHT(I,J)
C**** If WEIGHT=0 or A=SKIP, blank out the grid point
      IF(WEIGHT(I,J).NE.0. .AND. A(I,J).NE.SKIP)  GO TO 410
      LINE(I) = BLANK
      GO TO 450
C****
  410 APOINT(I,J) = A(I,J)*SCALE-DIF
      SUM(J) = SUM(J) + APOINT(I,J)*WEIGHT(I,J)
C****
      IF(APOINT(I,J).GT.-359.5)  GO TO 420
      LINE(I) = MINUS
      GO TO 450
C****
  420 IF(APOINT(I,J).LT.3599.5)  GO TO 430
      LINE(I) = PLUS
      GO TO 450
C****
  430 WRITE (XCHAR,'(F5.0)')  APOINT(I,J)
      LINE(I) = XCHAR(2:4)
C****
      IF(APOINT(I,J).GT.-99.5)  GO TO 440
      READ (XCHAR,'(I3)') KTENS
      LINE(I)(1:1) = MINUS(3:3)
      LINE(I)(2:2) = CHAR(-KTENS)
      GO TO 450
C****
  440 IF(APOINT(I,J).LT.999.5)  GO TO 450
      READ (XCHAR,'(I2)') KHUNDS
      LINE(I)(1:1) = CHAR(KHUNDS)
  450 CONTINUE
C****
C**** End of I loop, produce numbers for one latitude
C****
      IF(.NOT.QLAND)  GO TO 520
      DO 510 I=IMIN,IMAX
      IF(FOCEAN(I,J).GE..5)  GO TO 510
      IF(LINE(I)(1:1).EQ.' ')  LINE(I)(1:1) = ':'
      IF(LINE(I)(2:2).EQ.' ')  LINE(I)(2:2) = ':'
  510 CONTINUE
  520 IF(IHEMI.EQ.4)  GO TO 530
      WRITE (6,951) RLAT(J),J,(LINE(I),I=IMIN,IMAX),J,RLAT(J)
      GO TO 570
  530 IF(J.NE.1 .AND. J.NE.JM)  GO TO 540
      WT(J)  = WEIGHT(1,J)
      SUM(J) = 0.
      IF(WEIGHT(1,J).NE.0. .AND. A(1,J).NE.SKIP)
     *  SUM(J) = (A(1,J)*SCALE-DIF)*WEIGHT(1,J)
  540 ALAT(J) = 0.
      IF(WT(J).GT.0.)  ALAT(J) = SUM(J)/WT(J)
      WRITE (6,951) RLAT(J),J,(LINE(I),I=IMIN,IMAX),J,ALAT(J)
      HWT(JHEMI)  = HWT(JHEMI)  + WT(J)*DXYP(J)
      HSUM(JHEMI) = HSUM(JHEMI) + SUM(J)*DXYP(J)
  570 CONTINUE
C****
C**** End of J loop, produce bottom longitudes and global mean
C****
      IF(IHEMI.LT.4 .OR. JHEMI.EQ.2)
     *  WRITE (6,960) (I,I=IMIN,IMAX),(LON,LON=LONMIN,LONMAX,10)
  610 CONTINUE
      GSUM = HSUM(1) + HSUM(2)
      GWT  =  HWT(1) +  HWT(2)
      IF(HWT(1).NE.0)  HSUM(1) = HSUM(1)/HWT(1)
      IF(HWT(2).NE.0)  HSUM(2) = HSUM(2)/HWT(2)
      IF(GWT   .NE.0)  GSUM    = GSUM   /GWT
      LONMAX=LONMAX-10
      WRITE (6,961) (I,I=IMIN,IMAX),HSUM(2),HSUM(1),
     *              (LON,LON=LONMIN,LONMAX,10),GSUM
      RETURN
C****
C**** CALCULATE GLOBAL AND LATITUDINAL AVERAGES BUT DO NOT PRODUCE
C**** A LINE PRINTER MAP
C****
  700 IF(JM.NE.90)  RETURN
      DO 720 JHEMI=1,2
      HWT(JHEMI)  = 0.
      HSUM(JHEMI) = 0.
      JMAX=JM-(2-JHEMI)*(FJEQ-.5)
      JMIN= 1+(JHEMI-1)*(FJEQ-.5)
      DO 720 J=JMIN,JMAX
      IMAX=IM
      IF(J.EQ.1 .OR. J.EQ.JM)  IMAX=1
      WTJ  = 0.
      SUMJ = 0.
      DO 710 I=1,IMAX
      WTJ  = WTJ  + WEIGHT(I,J)
      APOINT(I,J) = SKIP
      IF(A(I,J).EQ.SKIP)  GO TO 710
      APOINT(I,J) = A(I,J)*SCALE-DIF
      SUMJ = SUMJ + APOINT(I,J)*WEIGHT(I,J)
  710 CONTINUE
      ALAT(J) = 0.
      IF(WTJ.NE.0.)  ALAT(J) = SUMJ/WTJ
      HWT(JHEMI)  = HWT(JHEMI)  + WTJ*DXYP(J)
  720 HSUM(JHEMI) = HSUM(JHEMI) + SUMJ*DXYP(J)
      GSUM = HSUM(1) + HSUM(2)
      GWT  = HWT(1)  + HWT(2)
      IF(HWT(1).NE.0.)  HSUM(1) = HSUM(1)/HWT(1)
      IF(HWT(2).NE.0.)  HSUM(2) = HSUM(2)/HWT(2)
      IF(GWT   .NE.0.)  GSUM    = GSUM   /GWT
      RETURN
C****
  930 FORMAT ('1  IHOUR =',I8,A86,I10,'  DAY',I4,', HOUR',I3//)
  931 FORMAT (3X,10I12/'0  LAT       ',36I3,8X,A4/)
  951 FORMAT (F6.1,I5,2X,36A3,I5,F7.1)
  960 FORMAT (/'0  LAT       ',36I3//3X,10I12)
  961 FORMAT (/'0  LAT       ',36I3,'   NH',F7.1/124X,'SH',F7.1/
     *  3X,9I12,13X,'GL',F7.1)
      END

      SUBROUTINE WMAP1 (JM,IHOUR,TITLE,A,WEIGHT,SCALE,DIF,KDSN)
C****
C**** WMAP1 prints on the line printer the contents of the array A
C**** which is dimensioned 12 by JM.  If JM = 24 and KDSN > 0, the
C**** output indicates the continents for the wonder resolution.
C****
      use TimeConstants_mod, only: INT_HOURS_PER_DAY, INT_DAYS_PER_YEAR
      PARAMETER (IM=12)
      REAL*4 A(12,1),WEIGHT(12,1)
      CHARACTER*80 TITLE,TITLEI
      COMMON /MAPCOM/ ALAT(24),ASHEMI,ANHEMI,AGLOB,APOINT(12,24)
      COMMON /WMAPFX/ RLAT(24),DXYP(24),PLAND(12,24)
      CHARACTER*1 CHAR(35)
      CHARACTER*4 LAND(12),LINE(12),BLANK,PERIOD,MINUS,PLUS
      CHARACTER*6 XCHAR
C****
      DATA JMLAST/0/, KDLAST/0/
      DATA SKIP /-1.e20/
      DATA CHAR/'1','2','3','4','5','6','7','8','9','A','B','C',
     *          'D','E','F','G','H','I','J','K','L','M','N','O',
     *          'P','Q','R','S','T','U','V','W','X','Y','Z'/
      DATA BLANK/'    '/, PERIOD/'  ..'/, MINUS/' ---'/, PLUS/' +++'/
      save jmlast,kdlast,skip,char,BLANK,PERIOD,MINUS,PLUS
C****
      IF(JM.EQ.JMLAST)  GO TO 200
      JMLAST = JM
C****
C**** CALCULATE SPHERICAL GEOMETRY
C****
      TWOPI  = 6.2831853
      RADIUS = 1.
      DLON = TWOPI/IM
      DLAT = NINT(180./(JM-1)) * TWOPI/360.
      FJEQ = .5*(1+JM)
      SINS = -1.
      DO 120 J=1,JM
      SINN = SIN(DLAT*(J+.5-FJEQ))
      IF(J.EQ.JM) SINN = 1.
      DXYP(J) = DLON*RADIUS*RADIUS*(SINN-SINS)
      RLAT(J) = DLAT*(J-FJEQ)*360/TWOPI
  120 SINS=SINN
      RLAT(1) = -90.
      RLAT(JM) = 90.
C****
C**** Read in PLAND from disk if JM = 24 and KDSN > 0
C****
  200 IF(KDSN.LT.0)  GO TO 700
      JYEAR = IHOUR/(INT_HOURS_PER_DAY*INT_DAYS_PER_YEAR)
      JDAY  = IHOUR/INT_HOURS_PER_DAY-INT_DAYS_PER_YEAR*JYEAR
      JHOUR = IHOUR-INT_HOURS_PER_DAY*(INT_DAYS_PER_YEAR*JYEAR+JDAY)
      JYEAR = JYEAR+1958
      IF(KDSN.EQ.KDLAST)  GO TO 300
      KDLAST=KDSN
      IF(JM.NE.24 .OR. KDSN.EQ.0)  GO TO 210
      READ  (KDSN) TITLEI,PLAND
      REWIND KDSN
  210 DO 220 I=1,IM
  220 LAND(I) = BLANK
C****
C**** WRITE TITLE AND LONGITUDES AT TOP OF PAGE
C****
  300 IF(TITLE(2:2).NE.'&')  WRITE (6,930) IHOUR,TITLE,JYEAR,JDAY,JHOUR
      WRITE (6,931) (LON,LON=-60,40,20),(I,I=1,IM)
C****
C**** OUTSIDE J LOOP
C****
      GWT    = 0.
      GACCUM = 0.
      HWT    = 0.
      HACCUM = 0.
      DO 600 J=JM,1,-1
      IF(J.NE.JM/2)  GO TO 310
      ANHEMI = 0.
      IF(HWT.GT.0.)  ANHEMI = HACCUM/HWT
      HWT    = 0.
      HACCUM = 0.
C**** Produce periods over land only if JM = 24 and KDSN > 0
  310 IF(JM.NE.24 .OR. KDSN.EQ.0)  GO TO 330
      DO 320 I=1,IM
      LAND(I) = BLANK
      IF(PLAND(I,J).GE..5)  LAND(I) = PERIOD
  320 CONTINUE
  330 WRITE (6,933) LAND
C****
C**** OUTSIDE I LOOP
C****
      WT    = 0.
      ACCUM = 0
      DO 500 I=1,IM
      WT    = WT + WEIGHT(I,J)
C**** If WEIGHT=0 or A=SKIP, blank out the grid point
      IF(WEIGHT(I,J).NE.0. .AND. A(I,J).NE.SKIP)  GO TO 410
      LINE(I) = BLANK
      GO TO 500
C****
  410 APOINT(I,J) = A(I,J)*SCALE-DIF
      ACCUM = ACCUM + APOINT(I,J)*WEIGHT(I,J)
C****
      IF(APOINT(I,J).GT.-3599.5)  GO TO 420
      LINE(I) = MINUS
      GO TO 500
C****
  420 IF(APOINT(I,J).LT.35999.5)  GO TO 430
      LINE(I) = PLUS
      GO TO 500
C****
  430 WRITE (XCHAR,'(F6.0)') APOINT(I,J)
      LINE(I) = XCHAR(2:5)
C****
      IF(APOINT(I,J).GT.-999.5)  GO TO 440
      READ (XCHAR,'(I3)') KHUNDS
      LINE(I)(1:1) = MINUS(3:3)
      LINE(I)(2:2) = CHAR(-KHUNDS)
      GO TO 500
C****
  440 IF(APOINT(I,J).LT.9999.5)  GO TO 500
      READ (XCHAR,'(I2)') KTHOUS
      LINE(I)(1:1) = CHAR(KTHOUS)
  500 CONTINUE
C****
C**** END OF OUTSIDE I LOOP, PRODUCE NUMBERS FOR ONE LATITUDE
C****
      IF(J.NE.1 .AND. J.NE.JM)  GO TO 510
      WT    = WEIGHT(1,J)
      ACCUM = 0.
      IF(A(1,J).NE.SKIP .AND. WEIGHT(1,J).NE.0.)
     *  ACCUM = (A(1,J)*SCALE-DIF)*WEIGHT(1,J)
  510 ALAT(J) = 0.
      IF(WT.GT.0.)  ALAT(J) = ACCUM/WT
      WRITE (6,951) RLAT(J),J,LINE,J,ALAT(J)
      HWT    = HWT    + WT*DXYP(J)
      HACCUM = HACCUM + ACCUM*DXYP(J)
      GWT    = GWT    + WT*DXYP(J)
  600 GACCUM = GACCUM + ACCUM*DXYP(J)
C****
C**** END OF OUTSIDE J LOOP, PRODUCE BOTTOM LONGITUDES AND GLOBAL MEAN
C****
      ASHEMI = 0.
      IF(HWT.GT.0.)  ASHEMI = HACCUM/HWT
      AGLOB = 0.
      IF(GWT.GT.0.)  AGLOB = GACCUM/GWT
      WRITE (6,961) (I,I=1,IM),ANHEMI,ASHEMI,
     *              (LON,LON=-60,40,20),AGLOB
      RETURN
C****
C**** CALCULATE GLOBAL AND LATITUDINAL AVERAGES BUT DO NOT PRODUCE A
C**** LINE PRINTER MAP
C****
  700 GWT    = 0.
      GACCUM = 0.
      HWT    = 0.
      HACCUM = 0.
      DO 720 J=1,JM
      IMAX=IM
      IF(J.EQ.1 .OR. J.EQ.JM)  IMAX=1
      IF(J.NE.JM/2+1)  GO TO 705
      ASHEMI = 0.
      IF(HWT.GT.0.)  ASHEMI = HACCUM/HWT
      HWT    = 0.
      HACCUM = 0.
  705 CONTINUE
      WT    = 0.
      ACCUM = 0.
      DO 710 I=1,IMAX
      APOINT(I,J) = SKIP
      WT    = WT    + WEIGHT(I,J)
      IF(A(I,J).EQ.SKIP)  GO TO 710
      APOINT(I,J) = A(I,J)*SCALE-DIF
      ACCUM = ACCUM + APOINT(I,J)*WEIGHT(I,J)
  710 CONTINUE
      ALAT(J) = 0.
      IF(WT.NE.0.)  ALAT(J) = ACCUM/WT
      HWT    = HWT    + WT*DXYP(J)
      HACCUM = HACCUM + ACCUM*DXYP(J)
      GWT    = GWT    + WT*DXYP(J)
  720 GACCUM = GACCUM + ACCUM*DXYP(J)
      ANHEMI = 0.
      IF(HWT.GT.0.)  ANHEMI = HACCUM/HWT
      AGLOB = 0.
      IF(GWT.NE.0.)  AGLOB = GACCUM/GWT
      RETURN
C****
  930 FORMAT ('1  IHOUR =',I8,A86,I10,'  DAY',I4,', HOUR',I3//)
  931 FORMAT (39X,6I8/'0',32X,'LAT       ',12I4,'        MEAN'/)
  933 FORMAT (43X,12A4)
  951 FORMAT (F36.1,I5,2X,12A4,I5,F7.1)
  961 FORMAT (/'0',32X,'LAT       ',12I4,'   NH',F7.1/94X,'SH',F7.1/
     *  39X,6I8,7X,'GL',F7.1)
      END
