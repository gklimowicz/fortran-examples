!     PROGRAM FOR MINIMISING A FUNCTION OF N VARIABLES USING
!     DIRECTION SET METHOD

      PROGRAM MINI
      IMPLICIT REAL*16(A-H,O-Z)
      EXTERNAL FCN
      DIMENSION X0(20),WK(900)

!     EXAMPLE 8.5: ROSENBROCK'S FUNCTION

51    FORMAT('   STARTING VALUES :',1P4D14.6/(2X,5D14.6))
52    FORMAT('   IER =',I4,5X,'MINIMISER =',1P4D14.6/(2X,5D14.6))
53    FORMAT('   NO. OF FUNCTION EVALUATIONS =',I7,5X,'MINIMUM =',
     1        1PD14.6)

      REPS=1.Q-33
      AEPS=1.Q-35
      NVAR=2

100   PRINT *,'TYPE STARTING VALUES X0(1) ... X0(NVAR)'
      PRINT *,'                 (QUITS WHEN X0(1)<-100)'
      READ *,(X0(I),I=1,NVAR)
      IF(X0(1).LT.-100) STOP

      WRITE(6,51) (X0(I),I=1,NVAR)
      CALL NMINF(NVAR,X0,F,NUM,REPS,AEPS,IER,FCN,WK)
      WRITE(6,52) IER,(X0(I),I=1,NVAR)
      WRITE(6,53) NUM,F
      GO TO 100
      END

!     -------------------------------------------

!	Function routine to calculate the function value as required for
!	line search without derivatives
!
!	FCN : (input) Name of subroutine to calculate the required function
!	X : (input) Parameter along the line to specifying the point where
!		function evaluation is required
!	V : (input/output) Real array of length 2N, first N elements specify the
!		direction of line search. After execution next N elements will contain
!		the coordinates of the point at which function is evaluated.
!	X0 : (input) Real array of length N, containing the coordinates
!		of the starting point for line search
!	N : (input) Number of variables in the function to be minimised
!	NUM : (input/output) Integer variable to keep count of function evaluations
!
!	SUBROUTINE FCN(N,X,FX) to calculate the required function, must be supplied
!		by the user. Here N is the number of variables, FX is the
!		function value at X. X is a real array of length N.
!
!	Required routines : FCN

      FUNCTION FLN(FCN,X,V,X0,N,NUM)
      IMPLICIT REAL*16(A-H,O-Z)
      DIMENSION V(N),X0(2*N)

      NUM=NUM+1
!	coordinates of the required points
      DO 1000 I=1,N
1000  X0(N+I)=X0(I)+V(I)*X
      CALL FCN(N,X0(N+1),FLN)
      END

!     --------------------------------------------

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
      IMPLICIT REAL*16(A-H,O-Z)
      EXTERNAL F
      PARAMETER(NIT=25,T1=9.Q0,GR=1.618034Q0,GC=0.381966Q0)
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

!     --------------------------------------------------

!	To minimise a function of several variables using direction set method
!
!	N : (input) Number of variables
!	X : (input/output) Real array of length N containing the initial
!		guess for the minimum.
!		After execution it should contain the coordinates of minimiser
!	F : (output) The function value at X
!	NUM : (output) Number of function evaluations used by the subroutine
!	REPS : (input) Required relative accuracy
!	AEPS : (input) Required absolute accuracy, iteration will stop when
!		change in function value is less than MAX(AEPS, REPS*ABS(F))
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=504 implies that N < 2, in which case no calculations are done
!		IER=529 implies that iteration failed to converge to specified accuracy
!		Other values may be set by LINMNF
!	FCN : (input) Name of the subroutine to calculate the function value
!	WK : Real array of length N*(2N+2) used as scratch space
!
!	SUBROUTINE FCN(N,X,F) to calculate the required function, must be supplied
!		by the user. Here N is the number of variables, F is the
!		function value at X. X is a real array of length N.
!
!	Required routines : LINMNF, FLN, SVD, FCN
!
      SUBROUTINE NMINF(N,X,F,NUM,REPS,AEPS,IER,FCN,WK)
      IMPLICIT LOGICAL(Q)
      IMPLICIT REAL*16(A-H,O,P,R-Z)
      EXTERNAL FCN
      PARAMETER(NIT=200)
      DIMENSION X(N),WK(N,2*(N+1))

      IER=0
      IF(N.LE.1) THEN
        IER=504
        RETURN
      ENDIF

!	Initialise the direction set matrix to identity matrix
      DO 2000 I=1,N
        WK(I,N+3)=0.0
        WK(I,N+4)=0.0
        DO 1800 J=1,N
1800    WK(J,I)=0.0
        WK(I,I)=1.0
2000  CONTINUE

      CALL FCN(N,X,F)
      NUM=1

!	The main iteration loop
      DO 6000 IT=1,NIT
        FI=F
        DO 4000 K=1,N

!	The starting point for line search
          DO 2400 I=1,N
2400      WK(I,N+1)=X(I)
          F0=F
          DFMAX=0.0
          KMAX=1

          DO 2800 I=1,N
            X1=0.0
!	Use previous value as initial approximation to minimum
            X2=WK(I,N+3)
            FUN=F0
            IF(X2.EQ.0.0) X2=1
            CALL LINMNF(X1,X2,F0,REPS,AEPS,IER1,FCN,WK(1,I),WK(1,N+1),
     1                  N,NUM)
            IF(IER1.GT.0) IER=IER1

            WK(I,N+3)=X1
!	Estimate of second derivative along the line
            IF(X1.NE.0.0) WK(I,N+4)=ABS((FUN-F0)/X1**2)
!	The new starting point
            DO 2600 J=1,N
2600        WK(J,N+1)=WK(J,N+1)+X1*WK(J,I)
            IF(FUN-F0.GE.DFMAX.AND.I.LE.N-K+1) THEN
              DFMAX=FUN-F0
              KMAX=I
            ENDIF
2800      CONTINUE

!	Remove the KMAX th direction
          DO 3000 I=KMAX,N-1
            WK(I,N+3)=WK(I+1,N+3)
            WK(I,N+4)=WK(I+1,N+4)
            DO 3000 J=1,N
3000      WK(J,I)=WK(J,I+1)
!	Add a new direction
          DO 3200 I=1,N
3200      WK(I,N)=WK(I,N+1)-X(I)

          X1=0.0
          X2=1
          FUN=F
          WK(N,N+4)=0.0
!	Starting point for the final line search in the loop
          DO 3400 I=1,N
3400      WK(I,N+1)=X(I)
          CALL LINMNF(X1,X2,F,REPS,AEPS,IER1,FCN,WK(1,N),WK(1,N+1),
     1                N,NUM)
          IF(IER1.GT.0) IER=IER1

          WK(N,N+3)=X1
          IF(X1.NE.0.0) WK(N,N+4)=ABS((FUN-F)/X1**2)
          DO 3600 J=1,N
3600      X(J)=X(J)+X1*WK(J,N)
4000    CONTINUE

        IF(ABS(F-FI).LT.MAX(REPS*ABS(F),AEPS)) RETURN

!	The matrix V for SVD
        DO 4500 J=1,N
          IF(WK(J,N+4).GT.0.0) THEN
            DO 4400 I=1,N
4400        WK(I,J)=WK(I,J)/SQRT(WK(J,N+4))
          ENDIF
4500    CONTINUE
        M=N
        LA=N
        LV=N
        CALL SVD(N,M,WK,WK(1,N+3),WK(1,N+1),LA,LV,WK(1,N+2),IER1)
        IF(IER1.GT.0) THEN
          IER=IER1
          RETURN
        ENDIF

        DO 4600 I=1,N
          WK(I,N+3)=SQRT(FUN-F)*WK(I,N+1)
          WK(I,N+4)=0.0
          IF(WK(I,N+1).NE.0.0) WK(I,N+4)=1./WK(I,N+1)**2
4600    CONTINUE

6000  CONTINUE
      IER=529
      END

!     ----------------------------------------------------

!	To calculate the Singular Value Decomposition of a matrix A=U D V-transpose
!
!	N : (input) Number of variables
!	M : (input) Number of equations
!	A : (input/output) Matrix of coefficients of size LA*N
!		After execution it will contain the matrix U
!	V : (output) The matrix V of size LV*N
!	SIGMA : (output) Array of length N, containing the singular values
!	LA : (input) Actual value of first dimension of A in the calling program
!	LV : (input) Actual value of first dimension of V in the calling program
!	E : Scratch array of length N
!	IER : (output) Error parameter, IER=0 if execution is successful
!		IER=12 QR iteration failed to converge to required accuracy
!		IER=105 implies N.LE.0, N.GT.LV, M.LE.0, M.GT.LA, N.GT.M
!
!	Required routines : None

      SUBROUTINE SVD(N,M,A,V,SIGMA,LA,LV,E,IER)
      IMPLICIT REAL*16(A-H,O-Z)
      PARAMETER(ITMAX=30,REPS=1.Q-33)
!	For REAL*4 use REPS=6.E-8
!      PARAMETER(ITMAX=30,REPS=6.E-8)
      DIMENSION A(LA,N),V(LV,N),SIGMA(N),E(N)

      IF(N.GT.M.OR.N.LE.0.OR.M.LE.0.OR.M.GT.LA.OR.N.GT.LV) THEN
        IER=105
        RETURN
      ENDIF

      IER=0
!	Reduction to Bidiagonal form using Householder transformations
      G=0
      RMAX=0

      DO 3000 I=1,N
!	Off-diagonal elements of bidiagonal form
        E(I)=G
        S=0
        DO 1200 J=I,M
1200    S=S+A(J,I)**2

        IF(S.LE.0.0) THEN
!	transformation not required
          G=0
        ELSE
          F=A(I,I)
          G=SQRT(S)
          IF(F.GE.0.0) G=-G
          H=F*G-S
          A(I,I)=F-G

          DO 1800 J=I+1,N
            S=0
            DO 1400 K=I,M
1400        S=S+A(K,I)*A(K,J)
            F=S/H
            DO 1600 K=I,M
1600        A(K,J)=A(K,J)+F*A(K,I)
1800      CONTINUE
        ENDIF

!	Diagonal elements of bidiagonal form
        SIGMA(I)=G
        S=0
        DO 2000 J=I+1,N
2000    S=S+A(I,J)**2

        IF(S.LE.0.0) THEN
!	Transformation not required
          G=0
        ELSE
          F=A(I,I+1)
          G=SQRT(S)
          IF(F.GE.0.0) G=-G
          H=F*G-S
          A(I,I+1)=F-G
          DO 2200 J=I+1,N
!	Temporary storage of intermediate results
2200      E(J)=A(I,J)/H

          DO 2800 J=I+1,M
            S=0
            DO 2400 K=I+1,N
2400        S=S+A(J,K)*A(I,K)
            DO 2600 K=I+1,N
2600        A(J,K)=A(J,K)+S*E(K)
2800      CONTINUE
        ENDIF
        R1=ABS(SIGMA(I))+ABS(E(I))
        IF(R1.GT.RMAX) RMAX=R1
3000  CONTINUE

!	Accumulation of right hand transformation in array V
      DO 4000 I=N,1,-1
        IF(G.NE.0.0) THEN
          H=A(I,I+1)*G
          DO 3200 J=I+1,N
3200      V(J,I)=A(I,J)/H

          DO 3800 J=I+1,N
            S=0
            DO 3400 K=I+1,N
3400        S=S+A(I,K)*V(K,J)
            DO 3600 K=I+1,N
3600        V(K,J)=V(K,J)+S*V(K,I)
3800      CONTINUE
        ENDIF

        DO 3900 J=I+1,N
          V(I,J)=0.0
          V(J,I)=0.0
3900    CONTINUE
        V(I,I)=1
        G=E(I)
4000  CONTINUE

!	Accumulation of left hand transformation overwritten on matrix A
      DO 5000 I=N,1,-1
        G=SIGMA(I)
        DO 4200 J=I+1,N
4200    A(I,J)=0
        IF(G.NE.0.0) THEN
          H=A(I,I)*G

          DO 4700 J=I+1,N
            S=0
            DO 4400 K=I+1,M
4400        S=S+A(K,I)*A(K,J)
            F=S/H
            DO 4600 K=I,M
4600        A(K,J)=A(K,J)+F*A(K,I)
4700      CONTINUE

          DO 4800 J=I,M
4800      A(J,I)=A(J,I)/G
        ELSE
          DO 4900 J=I,M
4900      A(J,I)=0.0
        ENDIF
        A(I,I)=A(I,I)+1
5000  CONTINUE

!	Diagonalisation of the bidiagonal form
      AEPS=REPS*RMAX
!	Loop over the singular values
      DO 8000 K=N,1,-1
!	The QR transformation
        DO 7500 ITR=1,ITMAX

!	Test for splitting
          DO 5200 L=K,1,-1
            IF(ABS(E(L)).LT.AEPS) GO TO 6000
            IF(ABS(SIGMA(L-1)).LT.AEPS) GO TO 5400
5200      CONTINUE

!	cancellation of E(L) if L>1
5400      C=0.0
          S=1.0
          DO 5800 I=L,K
            F=S*E(I)
            E(I)=C*E(I)
            IF(ABS(F).LT.AEPS) GO TO 6000
            G=SIGMA(I)
            SIGMA(I)=SQRT(F*F+G*G)
            C=G/SIGMA(I)
            S=-F/SIGMA(I)

            DO 5600 J=1,M
              R1=A(J,L-1)
              R2=A(J,I)
              A(J,L-1)=R1*C+R2*S
              A(J,I)=C*R2-S*R1
5600        CONTINUE
5800      CONTINUE

6000      Z=SIGMA(K)
          IF(L.EQ.K) THEN
!	QR iteration has converged
            IF(Z.LT.0.0) THEN
              SIGMA(K)=-Z
              DO 6200 J=1,N
6200          V(J,K)=-V(J,K)
            ENDIF
            GO TO 8000
          ENDIF

          IF(ITR.EQ.ITMAX) THEN
            IER=12
            GO TO 7500
          ENDIF

!	calculating shift from bottom 2x2 minor
          X=SIGMA(L)
          Y=SIGMA(K-1)
          G=E(K-1)
          H=E(K)
          F=((Y-Z)*(Y+Z)+(G-H)*(G+H))/(2.*H*Y)
          G=SQRT(1.+F*F)
          IF(F.LT.0.0) G=-G
          F=((X-Z)*(X+Z)+H*(Y/(F+G)-H))/X

!	next QR transformation
          C=1.0
          S=1.0
!	Given's rotation
          DO 7000 I=L+1,K
            G=E(I)
            Y=SIGMA(I)
            H=S*G
            G=C*G
            E(I-1)=SQRT(F*F+H*H)
            C=F/E(I-1)
            S=H/E(I-1)
            F=C*X+S*G
            G=C*G-S*X
            H=S*Y
            Y=C*Y

            DO 6400 J=1,N
              X=V(J,I-1)
              Z=V(J,I)
              V(J,I-1)=C*X+S*Z
              V(J,I)=C*Z-S*X
6400        CONTINUE

            SIGMA(I-1)=SQRT(F*F+H*H)
            IF(SIGMA(I-1).NE.0.0) THEN
              C=F/SIGMA(I-1)
              S=H/SIGMA(I-1)
            ENDIF
            F=C*G+S*Y
            X=C*Y-S*G
            DO 6600 J=1,M
              Y=A(J,I-1)
              Z=A(J,I)
              A(J,I-1)=C*Y+S*Z
              A(J,I)=C*Z-S*Y
6600        CONTINUE
7000      CONTINUE

          E(L)=0
          E(K)=F
          SIGMA(K)=X
7500    CONTINUE
8000  CONTINUE
      END

!     ------------------------------------------

      SUBROUTINE FCN(N,X,F)
      IMPLICIT REAL*16(A-H,O-Z)
      DIMENSION X(*)

!     ROSENBROCK'S FUNCTION

      F=100.*(X(2)-X(1)**2)**2+(1-X(1))**2
      END
