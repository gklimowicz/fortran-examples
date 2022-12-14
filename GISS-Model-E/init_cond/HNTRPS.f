C**** HNTRPS.FOR   Horizontal Interpolation Program REAL*8   2004/03/17
C****
      SUBROUTINE HNTRP0 (INA,JNA,OFFIA,DIVJA,
     *                   INB,JNB,OFFIB,DIVJB, SKIB)
C****
C**** HNTRP performs a horizontal interpolation of per unit area or per
C**** unit mass quantities defined on grid A, calculating the quantity
C**** on grid B.  B grid values that cannot be calculated because the
C**** covering A grid boxes have WTA = 0, are set to the value of SKIP.
C**** The area weighted integral of the quantity is conserved.
C****
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (TWOPI=6.283185307179586477d0)
C     REAL*4 WTA(INA,JNA),A(INA,JNA),B(INB,JNB),
      REAL*4 WTA(*),      A(*),      B(*),
     *       OFFIA,DIVJA, OFFIB,DIVJB, SKIB,SKIP
      COMMON /HNTRCB/ SINA(0:361),SINB(0:361),
     *       FMIN(720),FMAX(720),GMIN(361),GMAX(361),
     *       IMIN(720),IMAX(720),JMIN(362),JMAX(362)
      LOGICAL*4 QMPOLE
      DATA IMA,JMA,IMB,JMB/4*0/, SKIP/0/
C****
      IMA = INA
      JMA = JNA
      IMB = INB
      JMB = JNB
      SKIP = SKIB
      IF(IMA.lt.1 .or. IMA.gt.720 .or. JMA.lt.1 .or. JMA.gt.361 .or.
     *   IMB.lt.1 .or. IMB.gt.720 .or. JMB.lt.1 .or. JMB.gt.361)
     *   GO TO 400
C****
C**** Partitions in the I direction
C**** RIA = longitude in degrees of right edge of grid box IA on grid A
C**** RIB = longitude in degrees of right edge of grid box IB of grid B
C**** IMIN(IB) = box on grid A containing left edge of box IB on B
C**** IMAX(IB) = box on grid A containing right edge of box IB on B
C**** FMIN(IB) = fraction of box IMIN(IB) on A that is left of box IB
C**** FMAX(IB) = fraction of box IMAX(IB) on A that is right of box IB
C****
      DIA = 360d0/IMA
      DIB = 360d0/IMB
      IA  = 1
      RIA = (IA+OFFIA)*DIA - 360
      IB  = IMB
      DO 150 IBP1=1,IMB
      RIB = (IBP1-1+OFFIB)*DIB
  110 IF(RIA-RIB)  120,130,140
  120 IA  = IA  + 1
      RIA = RIA + DIA
      GO TO 110
C**** Right edge of A box IA and right edge of B box IB coincide
  130 IMAX(IB) = IA
      FMAX(IB) = 0.
      IA  = IA  + 1
      RIA = RIA + DIA
      IMIN(IBP1) = IA
      FMIN(IBP1) = 0.
      GO TO 150
C**** A box IA contains right edge of B box IB
  140 IMAX(IB) = IA
      FMAX(IB) = (RIA-RIB)/DIA
      IMIN(IBP1) = IA
      FMIN(IBP1) = 1.-FMAX(IB)
  150 IB = IBP1
      IMAX(IMB) = IMAX(IMB) + IMA
C       WRITE (0,*) 'IMIN=',(IMIN(I),I=1,IMB)
C       WRITE (0,*) 'IMAX=',(IMAX(I),I=1,IMB)
C       WRITE (0,*) 'FMIN=',(FMIN(I),I=1,IMB)
C       WRITE (0,*) 'FMAX=',(FMAX(I),I=1,IMB)
C****
C**** Partitions in the J direction
C****
C**** RJA = latitude in radians at top edge of box JA on grid A
C**** SINA(JA) = sine of latitude of top edge of box JA on grid A
      OFFJA = (DIVJA-JMA)/2.
      DJA   = .5*TWOPI/DIVJA
      DO 210 JA=1,JMA-1
      RJA = (JA+OFFJA)*DJA - .25*TWOPI
  210 SINA(JA) = DSIN(RJA)
      SINA(0)  = -1.
      SINA(JMA)=  1.
C**** RJB = latitude in radians at top edge of box JB on grid B
C**** SINB(JB) = sine of latitude of top edge of box JB on grid B
      OFFJB = (DIVJB-JMB)/2.
      DJB   = .5*TWOPI/DIVJB
      DO 220 JB=1,JMB-1
      RJB = (JB+OFFJB)*DJB - .25*TWOPI
  220 SINB(JB) = DSIN(RJB)
      SINB(0)  = -1.
      SINB(JMB)=  1.
C****
C**** JMIN(JB) = index of box of A that contains bottom edge of box JB
C**** JMAX(JB) = index of box of A that contains top edge of box JB
C**** GMIN(JB) = fraction of box JMIN(JB) on A grid that is below box JB
C**** GMAX(JB) = fraction of box JMAX(JB) on A grid that is above box JB
C****
      JMIN(1) = 1
      GMIN(1) = 0.
      JA = 1
      DO 350 JB=1,JMB-1
  310 IF(SINA(JA)-SINB(JB))  320,330,340
  320 JA = JA + 1
      GO TO 310
C**** Top edge of A box JA and top edge of B box JB coincide
  330 JMAX(JB) = JA
      GMAX(JB) = 0.
      JA = JA + 1
      JMIN(JB+1) = JA
      GMIN(JB+1) = 0.
      GO TO 350
C**** A box JA contains top edge of B box JB
  340 JMAX(JB) = JA
      GMAX(JB) = SINA(JA)-SINB(JB)
      JMIN(JB+1) = JA
      GMIN(JB+1) = SINB(JB)-SINA(JA-1)
  350 CONTINUE
      JMAX(JMB) = JMA
      GMAX(JMB) = 0.
C       WRITE (0,*) 'JMIN=',(JMIN(J),J=1,JMB)
C       WRITE (0,*) 'JMAX=',(JMAX(J),J=1,JMB)
C       WRITE (0,*) 'GMIN=',(GMIN(J),J=1,JMB)
C       WRITE (0,*) 'GMAX=',(GMAX(J),J=1,JMB)
      RETURN
C****
C**** Invalid parameters or dimensions out of range
C****
  400 WRITE (0,940) IMA,JMA,OFFIA,DIVJA, IMB,JMB,OFFIB,DIVJB, SKIP
      STOP 400
  940 FORMAT ('0Arguments received by HNTRP0 in order:'/
     *   2I12,' = IMA,JMA = array dimensions for A grid'/
     *  E24.8,' = OFFIA   = fractional number of grid boxes from',
     *                    ' IDL to left edge of grid box I=1'/
     *  E24.8,' = DIVJA   = number of whole grid boxes from SP to NP'/
     *   2I12,' = IMB,JMB = array dimensions for B grid'/
     *  E24.8,' = OFFIB   = fractional number of grid boxes from',
     *                    ' IDL to left edge of grid box I=1'/
     *  E24.8,' = DIVJB   = number of whole grid boxes from SP to NP'/
     *  E24.8,' = SKIP    = value to be put in B array when B',
     *  ' grid box is subset of A grid boxes with WTA = 0'/
     *  '0These arguments are invalid or out of range.')
C****

      ENTRY HNTRP (WTA,A,B)
C****
C**** HNTRP performs the horizontal interpolation
C**** Input: WTA = weighting array for values on the A grid
C****          A = per unit area or per unit mass quantity
C**** Output:  B = horizontally interpolated quantity on B grid
C****
      QMPOLE = .FALSE.
      GO TO 500

      ENTRY HNTRPP (WTA,A,B)
C****
C**** HNTRPP is similar to HNTRP but polar values are replaced by
C**** their longitudinal mean
C****
      QMPOLE = .TRUE.
C****
C**** Interpolate the A grid onto the B grid
C****
  500 DO 520 JB=1,JMB
      JAMIN = JMIN(JB)
      JAMAX = JMAX(JB)
      DO 520 IB=1,IMB
      IJB  = IB + IMB*(JB-1)
      WEIGHT= 0.
      VALUE = 0.
      IAMIN = IMIN(IB)
      IAMAX = IMAX(IB)
      DO 510 JA=JAMIN,JAMAX
      G = SINA(JA)-SINA(JA-1)
      IF(JA.eq.JAMIN)  G = G - GMIN(JB)
      IF(JA.eq.JAMAX)  G = G - GMAX(JB)
      DO 510 IAREV=IAMIN,IAMAX
      IA  = 1+MOD(IAREV-1,IMA)
      IJA = IA + IMA*(JA-1)
      F   = 1.
      IF(IAREV.eq.IAMIN)  F = F - FMIN(IB)
      IF(IAREV.eq.IAMAX)  F = F - FMAX(IB)
      WEIGHT = WEIGHT + F*G*WTA(IJA)
  510 VALUE  = VALUE  + F*G*WTA(IJA)*A(IJA)
      B(IJB) = SKIP
      IF(WEIGHT.ne.0.)  B(IJB) = VALUE/WEIGHT
  520 continue
C****
C**** Replace individual values near the poles by longitudinal mean
C****
      IF(.NOT.QMPOLE)  RETURN
      DO 630 JB=1,JMB,JMB-1
      WEIGHT = 0.
      VALUE  = 0.
      DO 610 IB=1,IMB
      IJB  = IB + IMB*(JB-1)
      IF(B(IJB).eq.SKIP)  GO TO 610
      WEIGHT = WEIGHT + 1.
      VALUE  = VALUE  + B(IJB)
  610 continue
      BMEAN = SKIP
      IF(WEIGHT.ne.0.)  BMEAN = VALUE/WEIGHT
      DO 620 IB=1,IMB
      IJB  = IB + IMB*(JB-1)
  620 B(IJB) = BMEAN
  630 continue
      RETURN
      END
