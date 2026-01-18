!	To draw a smooth curve passing through a set of data points
!	    using cubic spline interpolation
!
!	NTAB : (input) Number of points in the table
!	X : (input) Array of length NTAB containing X values
!	F : (input) Array of length NTAB containing function values at X(I)
!	C : (output) Array of length 3*NTAB which will contain the spline coefficients
!	NP : (input) Number of points at which interpolation is to be calculated
!	XP : (output) Real array of length NP containing the x values at
!	           NP uniformly spaced points for use in plotting
!	FP : (output) Real array of length NP containing interpolated
!		function values at XP(I)
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=202 implies NP.LE.1
!		other values may be set by SPLINE
!
!	Arrays XP and FP can be used to draw a smooth curve through the
!	tabulated points.
!
!	Required routines : SPLINE, SPLEVL

      SUBROUTINE SMOOTH(NTAB,X,F,C,NP,XP,FP,IER)
!      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(NTAB),F(NTAB),C(3,NTAB),XP(NP),FP(NP)

      CALL SPLINE(X,F,NTAB,C,IER)
      IF(IER.GT.100) RETURN

      IF(NP.LE.1) THEN
        IER=202
        RETURN
      ENDIF

      DX=(X(NTAB)-X(1))/(NP-1)
      DO 1000 I=1,NP
        XP(I)=X(1)+DX*(I-1)
        FP(I)=SPLEVL(XP(I),NTAB,X,F,C,DFB,DDFB,IER)
1000  CONTINUE
      END
