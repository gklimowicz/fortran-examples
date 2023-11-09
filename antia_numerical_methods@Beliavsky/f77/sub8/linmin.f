!	To perform a line search for minimum of several variables as required
!	by Quasi-Newton methods
!	This routine should not be used for any other purpose
!
!	X1 : (input/output) Starting value for the line search. After execution
!		it should contain the distance to minimiser along the line
!	X2 : (input/output) Initial estimate for the minimum along the line.
!		 This value will be modified by the subroutine.
!	F1 : (input/output) The function value at X1, this value must be supplied
!	DF1 : (output) The first derivative along the search direction at X1
!	REPS : (input) Required relative accuracy 
!	AEPS : (input) Required absolute accuracy,
!		These criterion is only used to terminate line search under
!		certain conditions and not generally applicable to the line
!		search. These values should be same as what is used in subroutine BFGS
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=54 implies that subroutine failed to find acceptable point
!			but the function value at last point is less than
!			that at beginning.
!		IER=55 implies that subroutine failed to find acceptable point
!			even though the interval has been reduced to required accuracy
!			but the function value at last point is less than
!			that at beginning.
!		IER=527 implies that iteration failed to find any point where
!			the function value is less than the starting value
!		IER=528 implies that iteration failed to find any point where
!			the function value is less than the starting value
!			even though the interval has been reduced to required accuracy
!	F : (input) Name of the subroutine to calculate the function value
!		and its derivatives
!	V : (input/output) Real array of length 3N. First N element specify
!		the direction in which minimisation is required. The next
!		N elements will contain the coordinates of the minimiser
!		found by LINMIN. The last 3N elements will contain the gradient
!		vector at the minimiser.
!	XI : (input) Real array of length N containing the coordinates of
!		starting point for line search
!	N : (input) Number of variables in the function to minimised
!	NUM : (output) Integer variable to keep count of the number of
!		function evaluations used so far
!
!	SUBROUTINE F(N,X,FX,G) to calculate the required function, must be supplied
!		by the user. Here N is the number of variables, FX is the
!		function value at X and G is the gradient vector. X and G
!		are real arrays of length N.
!
!	Required routines :  FLNM, F
!
      SUBROUTINE LINMIN(X1,X2,F1,DF1,REPS,AEPS,IER,F,V,XI,N,NUM)
      IMPLICIT LOGICAL(Q)
      IMPLICIT REAL*8(A-H,O,P,R-Z)
      EXTERNAL F
      PARAMETER(NIT=15,RHO=.01D0,SIGMA=.1D0,T1=9.,T2=.9D0,T3=.5D0)
      DIMENSION V(3*N),XI(N)

      IER=0
!	Select the bracketing phase
      QB=.FALSE.
      F2=FLNM(F,X2,DF2,V,XI,N,NUM)
      DX1=X2-X1
      F0=F1
      DF0=DF1
      X0=X1

      DO 2000 I=1,NIT
        FC=F0+DF0*RHO*(X2-X0)
        IF(ABS(DF2).LE.-SIGMA*DF0.AND.F2.LE.FC) THEN
!	Found an acceptable point
          X1=X2
          F1=F2
          DF1=DF2
          RETURN
        ENDIF
!	Test for bracketing
        IF(.NOT.QB) THEN
          IF(F2.GT.FC.OR.F2.GT.F1.OR.DF2.GE.0) QB=.TRUE.
        ENDIF

!	Hermite cubic interpolation
        DF12=(F2-F1)/DX1
        R=2.*DF2+DF1-3.*DF12
        R1=(3.*DF12-DF2-DF1)**2-DF1*DF2
        DX=0.0
        IF(R1.GT.0.0) THEN
          R1=SIGN(SQRT(R1),R)
          R=R+R1
          DX=-DF2*DX1/R
        ELSE
!	try parabolic interpolation
          R2=2.*(DF12-DF2)
          IF(R2.NE.0.0) DX=DX1*DF2/R2
        ENDIF

        IF(QB) THEN
!	Minimum is bracketed and hence improve on the bracket
          IF(DX.LT.-T2*DX1) DX=-T2*DX1
          IF(DX.GT.-T3*DX1) DX=-T3*DX1
          XA=X2+DX
          FA=FLNM(F,XA,DFA,V,XI,N,NUM)
          FC=F0+DF0*RHO*(XA-X0)

          IF(ABS(DFA).LE.-SIGMA*DF0.AND.FA.LE.FC) THEN
!	The new point is acceptable
            X1=XA
            F1=FA
            DF1=DFA
            RETURN
          ELSE IF(FA.GT.FC.OR.FA.GE.F1.OR.DFA.GT.0.0) THEN
            X2=XA
            F2=FA
            DF2=DFA
          ELSE
            X1=XA
            F1=FA
            DF1=DFA
          ENDIF

          DX1=X2-X1
          IF(ABS(DX1).LT.MAX(REPS*ABS(X2),AEPS)) THEN
!	If the interval is too small, then quit
            IER=528
            IF(F2.LE.F0) THEN
!	Accept the last point in any case
              X1=X2
              F1=F2
              DF1=DF2
              IER=55
            ENDIF
            RETURN
          ENDIF
        ELSE
!	Minimum hasn't been bracketed, choose the point further down.
          IF(DX.LT.X2-X1.AND.DX.GT.0) DX=X2-X1
          IF(DX.GT.T1*(X2-X1).OR.DX.LE.0.0) DX=T1*(X2-X1)
          X1=X2
          X2=X2+DX
          DX1=DX
          F1=F2
          DF1=DF2
          F2=FLNM(F,X2,DF2,V,XI,N,NUM)
        ENDIF
2000  CONTINUE

!	Iteration has failed to find an acceptable point
      IF(F2.LE.F0) THEN
!	accept this point if function value is smaller
        F1=F2
        X1=X2
        DF1=DF2
        IER=54
        RETURN
      ELSE IF(F1.LE.F0.AND.X1.NE.X0) THEN
        IER=54
        RETURN
      ENDIF

!	No acceptable point found
      IER=527
      END
