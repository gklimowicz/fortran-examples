!	To evaluate the cubic spline interpolant at a specified point
!
!	XB : (input) point at which interpolation is required
!	N : (input) Number of points in the table
!	X : (input) Array of length N, containing the abscissas
!	F : (input) Array of length N, containing the function values at X(I)
!	C : (input) Array of length 3*N containing the spline coefficients
!		which should have been calculated using SPLINE
!	DFB : (output) First derivative of spline at x=XB
!	DDFB : (output) Second derivative of spline at x=XB
!	IER : (output) error parameter, IER=0 if execution is successful
!		IER=24 implies XB is outside the range of table on higher side
!		IER=25 implies XB is outside the range of table on lower side
!		IER=201 implies N<2
!	SPLEVL will be the interpolated value at x=XB
!
!	Required routines : None

      FUNCTION SPLEVL(XB,N,X,F,C,DFB,DDFB,IER)
      IMPLICIT REAL*16(A-H,O,P,R-Z)
      IMPLICIT LOGICAL(Q)
      DIMENSION X(N),F(N),C(3,N)
      SAVE
      DATA LOW/0/

      SPLEVL=0.0
      IF(N.LE.1) THEN
        IER=201
        RETURN
      ENDIF

!	QASCND is true if table is in ascending order
      QASCND=X(N).GT.X(1)
      IER=0

      IF(LOW.LT.1.OR.LOW.GE.N) THEN
!	If the previous value of LOW is inadmissible, set the range to (1,N)
        LOW=1
        IGH=N
      ELSE
        IGH=LOW+1
      ENDIF

1000  IF((XB.LT.X(LOW).AND.XB.LT.X(IGH)).OR.
     1    (XB.GT.X(LOW).AND.XB.GT.X(IGH))) THEN
!	Extend the range
        IF(XB.GT.X(LOW).EQV.QASCND) THEN
!	Extend the range on higher side
          IF(IGH.GE.N) THEN
            IER=24
            LOW=N-1
          ELSE
            NIGH=MIN(N,IGH+2*(IGH-LOW))
            LOW=IGH
            IGH=NIGH
            GO TO 1000
          ENDIF
        ELSE
!	Extend the range on lower side
          IF(LOW.LE.1) THEN
            IER=25
          ELSE
            NIGH=LOW
            LOW=MAX(1,LOW-2*(IGH-LOW))
            IGH=NIGH
            GO TO 1000
          ENDIF
        ENDIF
      ELSE

!	Once the point is bracketed between two tabular points locate it by bisection
1500    IF(IGH-LOW.GT.1.AND.XB.NE.X(LOW)) THEN
          MID=(LOW+IGH)/2
          IF(XB.LE.X(MID).EQV.XB.LE.X(LOW)) THEN
            LOW=MID
          ELSE
            IGH=MID
          ENDIF
          GO TO 1500
        ENDIF
      ENDIF

!	Calculate the interpolant and its derivatives
      DX=XB-X(LOW)
      SPLEVL=((C(3,LOW)*DX+C(2,LOW))*DX+C(1,LOW))*DX+F(LOW)
      DFB=(3.*C(3,LOW)*DX+2.*C(2,LOW))*DX+C(1,LOW)
      DDFB=6.*C(3,LOW)*DX+2.*C(2,LOW)
      END
