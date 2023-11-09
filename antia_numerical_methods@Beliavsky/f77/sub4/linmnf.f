!	To perform a line search for minimum of several variables as required
!	by direction set method
!	This routine should not be used for any other purpose
!
!	X0 : (input/output) Starting value for the line search. After execution
!		it should contain the distance to minimiser along the line
!	X1 : (input/output) Initial estimate for the minimum along the line.
!		 This value will be modified by the subroutine.
!	F0 : (input/output) The function value at X1, this value must be supplied
!	REPS : (input) Required relative accuracy 
!	AEPS : (input) Required absolute accuracy,
!		This criterion is only used to terminate line search under
!		certain conditions and not generally applicable to the line
!		search. These values should be same as what is used in subroutine NMINF
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=56 implies that subroutine failed to find acceptable point
!			in this case the starting point itself is accepted
!	F : (input) Name of the subroutine to calculate the function value
!	V : (input/output) Real array of length 2N. First N element specify
!		the direction in which minimisation is required. The next
!		N elements will contain the coordinates of the minimiser
!		found by LINMNF.
!	XI : (input) Real array of length N containing the coordinates of
!		starting point for line search
!	N : (input) Number of variables in the function to be minimised
!	NUM : (output) Integer variable to keep count of the number of
!		function evaluations used so far
!
!	SUBROUTINE F(N,X,FX) to calculate the required function, must be supplied
!		by the user. Here N is the number of variables, FX is the
!		function value at X. X is a real array of length N.
!
!	Required routines :  FLN, F
!
      SUBROUTINE LINMNF(X0,X1,F0,REPS,AEPS,IER,F,V,XI,N,NUM)
!      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL F
      PARAMETER(NIT=25,T1=9.D0,GR=1.618034D0,GC=0.381966D0)
      DIMENSION V(N),XI(2*N)

      IER=0
      F1=FLN(F,X1,V,XI,N,NUM)
      IF(F1.LT.F0) THEN
!	Choose the next point further down on the same side
        X2=X0+GR*(X1-X0)
      ELSE
!	choose the next point on opposite side
        X2=X0-GR*(X1-X0)
      ENDIF
      F2=FLN(F,X2,V,XI,N,NUM)

      DO 2000 I=1,NIT
!	Parabolic interpolation
        R=(X2-X1)*(F2-F0)
        T=(X2-X0)*(F2-F1)
        P=(X2-X0)*T-(X2-X1)*R
        T=2.*(T-R)
        IF(T.GT.0) P=-P
        T=ABS(T)
        XP=X0
        IF(ABS(P).LT.ABS(T1*T*(X2-X1))) THEN
!	Try the interpolated value
          XP=X2+P/T
          FP=FLN(F,XP,V,XI,N,NUM)
          IF(FP.LT.F0) THEN
!	accept the point
            X0=XP
            F0=FP
            RETURN
          ENDIF
        ENDIF

        IF(F1.LT.F0) THEN
          X0=X1
          F0=F1
          RETURN
        ELSE IF(F2.LT.F0) THEN
          X0=X2
          F0=F2
          RETURN
        ELSE IF(XP.NE.X0) THEN
!	subdivide the interval
          IF(XP-X0.GT.0.0.EQV.X1-X0.GT.0.0) THEN
            IF(ABS(XP-X0).LT.0.5*ABS(X1-X0)) THEN
              X1=XP
              F1=FP
            ELSE
!	use golden section
              X1=X0+GC*(X1-X0)
              F1=FLN(F,X1,V,XI,N,NUM)
            ENDIF
          ELSE
            IF(ABS(XP-X0).LT.0.5*ABS(X2-X0)) THEN
              X2=XP
              F2=FP
            ELSE
!	use golden section
              X2=X0+GC*(X2-X0)
              F2=FLN(F,X2,V,XI,N,NUM)
            ENDIF
          ENDIF

!	If interpolated point is not acceptable use golden section
        ELSE IF(ABS(X2-X0).GT.ABS(X1-X0)) THEN
          X2=X0+GC*(X2-X0)
          F2=FLN(F,X2,V,XI,N,NUM)
        ELSE
          X1=X0+GC*(X1-X0)
          F1=FLN(F,X1,V,XI,N,NUM)
        ENDIF

!	If the change in function value is too small, then quit
        IF(MIN(F1-F0,F2-F0).LT.MAX(REPS*ABS(F0),AEPS)) RETURN
2000  CONTINUE

!	fails to find an acceptable point
      IER=56
      END
