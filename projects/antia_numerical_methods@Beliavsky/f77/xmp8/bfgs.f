!     PROGRAM TO MINIMISE A FUNCTION OF SEVERAL VARIABLES USING 
!     QUASI-NEWTON METHOD WITH BFGS FORMULA.

      PROGRAM MINMIS
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL FCN
      DIMENSION X0(20),G(20),H(2,200),WK(60)

!     EXAMPLE 8.4: ROSENBROCK'S FUNCTION

51    FORMAT('   STARTING VALUES =',1P4D14.6/(2X,5D14.6))
52    FORMAT('  IER =',I4,5X,'MINIMISER:',1P4D14.6/(2X,5D14.6))
53    FORMAT('  NO. OF FUNCTION EVALUATIONS =',I6,5X,'MINIMUM ='
     1       ,1PD14.6/'    GRADIENT VECTOR =',4D14.6/(2X,5D14.6))
54    FORMAT(2X,1P5D14.6)
55    FORMAT(5X,'INVERSE OF HESSIAN MATRIX')

      REPS=1.D-11
      AEPS=1.D-13
      NVAR=2

100   PRINT *,'TYPE STARTING VALUES X0(1) TO X0(NVAR)'
      PRINT *,'                 (QUITS WHEN X0(1)<-100)'
      READ *,(X0(I),I=1,NVAR)
      IF(X0(1).LT.-100) STOP
      WRITE(6,51) (X0(I),I=1,NVAR)
      CALL BFGS(NVAR,X0,F,G,H,NUM,REPS,AEPS,IER,FCN,WK)
      WRITE(6,52) IER,(X0(I),I=1,NVAR)
      WRITE(6,53) NUM,F,(G(I),I=1,NVAR)
      WRITE(6,55)
      DO 200 I=1,NVAR
200   WRITE(6,54) (H(I,J),J=1,NVAR)
      GO TO 100
      END

!     -------------------------------------------

!	To minimise a function of several variables using quasi-Newton method
!
!	N : (input) Number of variables
!	X : (input/output) Real array of length N containing the initial
!		guess for the minimum.
!		After execution it should contain the coordinates of minimiser
!	F : (output) The function value at X
!	G : (output) Real array of length N containing the gradient vector at X
!	H : (output) Real array of length N*N containing the estimated
!		Hessian matrix at X. The first dimension of H should be N.
!	NUM : (output) Number of function evaluations used by the subroutine
!	REPS : (input) Required relative accuracy
!	AEPS : (input) Required absolute accuracy, all components of the
!		Minimiser will be calculated with accuracy MAX(AEPS, REPS*ABS(X(I)))
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=53 implies that Hessian matrix at final point is probably singular
!			The iteration may have converged to a saddle point
!		IER=503 implies that N < 1, in which case no calculations are done
!		IER=526 implies that iteration failed to converge to specified accuracy
!		Other values may be set by LINMIN
!	FCN : (input) Name of the subroutine to calculate the function value
!		and its derivatives
!	WK : Real array of length 3N used as scratch space
!
!	SUBROUTINE FCN(N,X,F,G) to calculate the required function, must be supplied
!		by the user. Here N is the number of variables, F is the
!		function value at X and G is the gradient vector. X and G
!		are real arrays of length N. F and G must be calculated by FCN.
!
!	Required routines : LINMIN, FLNM, FCN
!
      SUBROUTINE BFGS(N,X,F,G,H,NUM,REPS,AEPS,IER,FCN,WK)
      IMPLICIT LOGICAL(Q)
      IMPLICIT REAL*8(A-H,O,P,R-Z)
      EXTERNAL FCN
      PARAMETER(NIT=200)
      DIMENSION X(N),WK(3*N),H(N,N),G(N)

      IER=0
      IF(N.LT.1) THEN
        IER=503
        RETURN
      ENDIF

      DO 2000 I=1,N
!	Initialise the Hessian matrix to unit matrix
        DO 1800 J=1,N
1800    H(J,I)=0.0
        H(I,I)=1.0
2000  CONTINUE
!	If some variable needs to be kept fixed the corresponding
!	diagonal elements, H(I,I) should be set to zero.

      CALL FCN(N,X,F,G)
      DF=ABS(F)
      IF(DF.EQ.0.0) DF=1
      N2=2*N
      NUM=1
      H2=1.

!	The iteration loop
      DO 4000 IT=1,NIT
        DF1=0.0
        H1=H2
        H2=0.0
!	Calculating the search direction WK =S^(k)
        DO 2400 I=1,N
          WK(I)=0.0
          DO 2200 J=1,N
            H2=H2+ABS(H(I,J))
2200      WK(I)=WK(I)-H(I,J)*G(J)
          DF1=DF1+WK(I)*G(I)
2400    CONTINUE

        IF(DF1.EQ.0.0) THEN
!	If gradient vanishes, then quit
!	If Hessian matrix appears to be singular, set the error flag
          IF(ABS(H2/H1).GT.1.3D0) IER=53
          RETURN
        ENDIF

!	Initial guess for line search
        X1=0
        X2=MIN(1.D0,-2.*DF/DF1)
        F1=F
        IF(X2.LE.0.0) X2=1
        CALL LINMIN(X1,X2,F1,DF1,REPS,AEPS,IER1,FCN,WK,X,N,NUM)
        IF(IER1.GT.0) IER=IER1
!	If line search fails, then quit
        IF(IER1.GT.100) RETURN

!	The convergence test
        QC=.TRUE.
        DO 2800 I=1,N
          X(I)=WK(N+I)
          WK(N+I)=X1*WK(I)
          IF(ABS(WK(N+I)).GT.MAX(REPS*ABS(X(I)),AEPS)) QC=.FALSE.
          GI=WK(N2+I)
          WK(N2+I)=WK(N2+I)-G(I)
          G(I)=GI
2800    CONTINUE
!	It is possible to apply convergence check on Function value using DF
!	instead of X(I)
        DF=F-F1
        F=F1
        IF(QC) THEN
          IF(ABS(H2/H1).GT.1.3D0) IER=53
          RETURN
        ENDIF

!	Update the matrix using BFGS formula
        DO 3200 I=1,N
          WK(I)=0.0
          DO 3200 J=1,N
3200    WK(I)=WK(I)+H(I,J)*WK(N2+J)
        GHG=0.0
        DG=0.0
        DO 3400 I=1,N
          DG=DG+WK(N2+I)*WK(N+I)
3400    GHG=GHG+WK(I)*WK(N2+I)
        R1=(1.+GHG/DG)/DG
        DO 3600 J=1,N
          DO 3600 I=1,N
3600   H(I,J)=H(I,J)+R1*WK(N+I)*WK(N+J)-(WK(N+I)*WK(J)+WK(I)*WK(N+J))/DG

4000  CONTINUE

!	Iteration fails to converge
      IER=526
      END

!     --------------------------------------------

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

!     --------------------------------------------------

!	Function routine to calculate the function value and its derivative
!	as required for line search
!
!	FCN : (input) Name of subroutine to calculate the required function
!	X : (input) Parameter along the line to specify the point where
!		function evaluation is required
!	DF : (output) First derivative of function along the line at X
!	V : (input/output) Real array of length 3N, first N elements specify the
!		direction of line search. Next N elements will contain the
!		coordinates of the point at which function is evaluated,
!		while the last N elements contain the gradient vector at the point
!	X0 : (input) Real array of length N, containing the coordinates
!		of starting point for line search
!	N : (input) Number of variables in the function to be minimised
!	NUM : (input/output) Integer variable to keep count of function evaluations
!
!	SUBROUTINE FCN(N,X,FX,G) to calculate the required function, must be supplied
!		by the user. Here N is the number of variables, FX is the
!		function value at X and G is the gradient vector. X and G
!		are real arrays of length N.
!
!	Required routines : FCN

      FUNCTION FLNM(FCN,X,DF,V,X0,N,NUM)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION V(3*N),X0(N)

      NUM=NUM+1
      N2=2*N
!	The coordinates of the required point
      DO 1000 I=1,N
1000  V(N+I)=X0(I)+V(I)*X

      CALL FCN(N,V(N+1),FLNM,V(N2+1))
!	The first derivative along the search direction
      DF=0.0
      DO 2000 I=1,N
2000  DF=DF+V(I)*V(N2+I)
      END

!     ----------------------------------------------------

      SUBROUTINE FCN(N,X,F,G)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(*),G(*)

!     ROSENBROCK'S FUNCTION AND ITS GRADIENT

      F=100.*(X(2)-X(1)**2)**2+(1-X(1))**2+1.
      G(1)=-400.*(X(2)-X(1)**2)*X(1)-2.*(1.-X(1))
      G(2)=200.*(X(2)-X(1)**2)
      END
