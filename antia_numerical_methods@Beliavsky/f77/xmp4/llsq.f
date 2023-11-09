!     PROGRAM TO OBTAIN A LEAST SQUARES APPROXIMATION USING
!     SINGULAR VALUE DECOMPOSITION IN 3 DIMENSIONS
 
      PROGRAM LFIT
!      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL PHI
      DIMENSION G(10000,50),X(3,9000),B(10000),FX(20,20,20),V(50,50)
      DIMENSION YF(20,20,20),WK(100),EF(20,20,20),D(50),X0(100)
      DIMENSION COV(100,100)
 
!     FUNCTION WITH RANDOM ERROR TO GENERATE INPUT DATA
 
      F(Y,Y1,Y2)=((231*Y*Y1-315)*Y*Y2+105)*Y1*Y2-5+RANGAU(SEED)*1.D-3
 
51    FORMAT('   IER =',I4,5X,'NO. OF DATA PTS =',I4,5X,
     1       'NUMBER OF BASIS FUNCTIONS =',I3/5X,'CHI SQUARE =',1PE12.4)
52    FORMAT('   COEF =',1P5E14.6/(4X,5E14.6))
 
      REPS=1.E-7
      K=3
      IX=3
      IU=10000
      IV=50
      SEED=11
!	No. of data points is fixed by dimension
      N=20
 
!    For fit in 3 dimensions the number of basis functions should be
!    M0**3 where M0 is number of basis functions along each dimension.
100   PRINT *, 'TYPE  M=NUMBER OF BASIS FUNCTIONS'
      PRINT *,'        (QUITS WHEN M.LE.0)'
      READ *,M
      IF(M.LE.0) STOP
 
!     GENERATING THE DATA POINTS WITH RANDOM ERROR
 
      H=1.D0/(N-1.)
      DO 1000 I=1,N
        X0(I)=(I-1)*H
1000  CONTINUE
      DO 2000 I=1,N
        DO 2000 J=1,N
          DO 2000 J1=1,N
            FX(I,J,J1)=F(X0(I),X0(J),X0(J1))
            EF(I,J,J1)=1.E-3
            X(1,I+N*(J-1)+N*N*(J1-1))=X0(I)
            X(2,I+N*(J-1)+N*N*(J1-1))=X0(J)
            X(3,I+N*(J-1)+N*N*(J1-1))=X0(J1)
2000  CONTINUE
      N1=N*N*N
 
      CALL LLSQ(N1,M,K,X,IX,FX,EF,B,G,V,IU,IV,D,YF,WK,PHI,
     1          REPS,CS,COV,IER)
      WRITE(6,51) IER,N1,M,CS
      WRITE(6,52) (B(I),I=1,M)
      GO TO 100
      END
 
!     -----------------------------------------------------
 
!	Linear least squares fit in K dimensions
!
!	N : (input) Number of data points to be fitted
!	M : (input) Number of basis functions to be used
!	K : (input) Number of dimensions
!	X : (input) Real array of length IX*N containing the coordinates
!		of points at which function values is available
!		X(I,J) is the Ith coordinate of Jth point
!		The points may have arbitrary distribution in K-space
!	IX : (input) First dimension of X in the calling program
!		IX .GE. K 
!	F : (input) Real array of length N containing the function values
!		F(I) should be function value at X(1,I),...,X(K,I)
!	EF : (input) Real array of length N containing the estimated error
!		in F(I). 
!	A : (output) Real array of length N containing the fitted coefficients
!		Note that although the number of coefficients is M, the
!		rest of array is used as scratch space
!	U : (output) Real array of length IU*M containing the matrix U of SVD
!		of the design matrix
!	V : (output) Real array of length IV*M containing the matrix V of SVD
!		of the design matrix
!	IU : (input) First dimension of U in the calling program (IU.GE.N)
!	IV : (input) First dimension of V in the calling program (IV.GE.M)
!	SIGMA : (output) Real array of length M containing the singular values
!		of the design matrix
!	Y : (output) Real array of length N containing the values of fitted
!		function at each of the tabular points
!	WK : Real array of length M used as scratch space
!	PHI : (input) Name of subroutine to calculate the basis functions
!		at any given point
!	REPS : (input) Required accuracy for solution of equations using SVD
!		singular values less than REPS times maximum will be set to zero
!	CHISQ : (output) The value of Chi square at minimum
!       COV : (output) Array of length IV*M containing the covariance
!               matrix of the fitted parameters. COV(I,I) will be the
!               variance in A(I).
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=606 implies that M>N, M.LE.0, N.LE.0, or K>IX
!		IER=607 implies that EF(I).LE.0 for some I
!		No calculations are done in both these cases
!		Other values of IER may be set by SVD
!
!	The SUBROUTINE PHI(M,X,Y) must be supplied by the user to calculate
!	the required basis functions. M is the number of basis functions,
!	X is a real array of length K containing the coordinates of point
!	at which the basis function needs to be calculated. Y is a real
!	array of length M containing the computed basis functions at X
!
!       THE ARGUMENTS OF THIS FUNCTION HAVE CHANGED FROM THE EARLIER VERSION.
!	NOW THERE IS AN ADDITIONAL ARGUMENT COV TO CALCULATE THE COVARIANCE
!	MATRIX.
!
!	Required routines : SVD, SVDEVL, PHI
!	
      SUBROUTINE LLSQ(N,M,K,X,IX,F,EF,A,U,V,IU,IV,SIGMA,Y,WK,PHI,
     1      REPS,CHISQ,COV,IER)
!      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(IX,N),F(N),EF(N),A(N),Y(N),WK(M)
      DIMENSION U(IU,M),V(IV,M),SIGMA(M),COV(IV,M)
 
      IF(M.GT.N.OR.M.LE.0.OR.N.LE.0.OR.K.GT.IX) THEN
        IER=606
        RETURN
      ENDIF
 
!	Setting up the design matrix and the RHS
      IER=0
      DO 2000 I=1,N
        IF(EF(I).LE.0) THEN
          IER=607
          RETURN
        ENDIF
        A(I)=F(I)/EF(I)
        CALL PHI(M,X(1,I),WK)
        DO 2000 J=1,M
          U(I,J)=WK(J)/EF(I)
2000  CONTINUE
 
      CALL SVD(M,N,U,V,SIGMA,IU,IV,WK,IER)
      IF(IER.GT.100) RETURN
 
!	Calculate the least squares solution
      CALL SVDEVL(M,N,U,V,SIGMA,IU,IV,A,WK,REPS)
 
!	Computing the \chi^2 from fitted coefficients
      CHISQ=0.0
      DO 3000 I=1,N
        CALL PHI(M,X(1,I),WK)
        S1=0.0
        DO 2500 J=1,M
          S1=S1+A(J)*WK(J)
2500    CONTINUE
        Y(I)=S1
        CHISQ=CHISQ+((F(I)-Y(I))/EF(I))**2
3000  CONTINUE

!       Computing the covariance matrix
      SIGMAX=0
      DO I=1,M
        IF(SIGMA(I).GT.SIGMAX) SIGMAX=SIGMA(I)
      ENDDO
      DO I=1,M
      DO J=1,I
        COV(J,I)=0.0
        DO IK=1,M
        IF(SIGMA(IK).GT.REPS*SIGMAX) COV(J,I)=COV(J,I)
     1            +V(J,IK)*V(I,IK)/SIGMA(IK)**2
        ENDDO
        COV(I,J)=COV(J,I)
      ENDDO
      ENDDO
 
      END
 
!     ---------------------------------------------------------
 
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
!      IMPLICIT REAL*8(A-H,O-Z)
!      PARAMETER(ITMAX=30,REPS=1.D-16)
!	For REAL*4 use REPS=6.E-8
      PARAMETER(ITMAX=30,REPS=6.E-8)
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
 
!     --------------------------------------------------------------
 
!	To evaluate the solution of a system of linear equations using SVD
!
!	N : (input) Number of variables
!	M : (input) Number of equations
!	U : (input) array of size LU*N containing the left-hand transformation
!	V : (input) array of size LV*N containing the right-hand transformation
!	SIGMA : (input) array of size N containing the singular values
!	LU : (input) First dimension of array U in the calling program
!	LV : (input) First dimension of array V in the calling program
!	B : (input/output) Array of length M containing the RHS
!		after execution it will contain the solution
!	WK : Scratch array of length N
!	REPS : (input) Relative accuracy.
!               All singular values < REPS*(Max of singular values)
!		will be reduced to zero
!
!	Required routines : None

      SUBROUTINE SVDEVL(N,M,U,V,SIGMA,LU,LV,B,WK,REPS)
!      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION U(LU,N),V(LV,N),SIGMA(N),B(*),WK(N)

!	Finding the largest singular value
      SMAX=0.0
      DO 2000 I=1,N
        IF(SIGMA(I).GT.SMAX) SMAX=SIGMA(I)
2000  CONTINUE

      AEPS=SMAX*REPS
      DO 3000 I=1,N
        S=0.0
!	Only SIGMA(I) > AEPS contribute to the solution
        IF(SIGMA(I).GT.AEPS) THEN
          DO 2400 J=1,M
2400      S=S+U(J,I)*B(J)
          S=S/SIGMA(I)
        ENDIF
        WK(I)=S
3000  CONTINUE

      DO 4000 I=1,N
        S=0.0
        DO 3400 J=1,N
3400    S=S+V(I,J)*WK(J)
        B(I)=S
4000  CONTINUE
      END
 
!     ---------------------------------------------------
 
!	To generate random numbers with Gaussian probability distribution
!	It generates random numbers with zero mean and variance of 1.
!	
!	SEED : (input/output) real seed, it should be positive and
!		less than AM. It is updated by the routine and should
!		not be modified between two calls, unless a fresh
!		sequence is required
!
!       THE ARGUMENT OF THIS FUNCTION HAS CHANGED AS COMPARED TO EARLIER
!       VERSION AS THE SEED IS NOW REAL INSTEAD OF INTEGER.
!
!	Required routines : None

      FUNCTION RANGAU(SEED)
!      IMPLICIT REAL*8(A-H,O-Z)
!	Retain the following declaration even for REAL*4 version
!	otherwise AN and A will be rounded
      REAL*8 AM,A,AC,AN
      PARAMETER(AM=2147483648D0,A=45875D0,AC=453816693D0)
      PARAMETER(PI=3.14159265358979324D0,AN=2147483647D0)

      R1=MOD(A*SEED+AC,AM)
      IF(SEED.EQ.0.0) SEED=0.1D0
      RANGAU=SQRT(2.D0*LOG(AN/SEED))*COS(2.0*PI*R1/AN)
      SEED=R1
      END
 
!     ---------------------------------------------------
 
!	  To calculate the set of basis functions in 3 dimensions
!     These are the monomials in 3 D in natural order

      SUBROUTINE PHI(N,X,F)
!      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION F(N),X(3)
 
!	The number N is assumed to be a cube of integer
      M1=(N*1.D0)**0.3333D0+0.1
      M2=N/M1**2
      DO 1000 I=1,M1
        DO 1000 I1=1,M1
          DO 1000 J=1,M2
            F1=1.0
            F2=1.0
            F3=1.0
            IF(I.GT.1) F1=X(1)**(I-1)
            IF(I1.GT.1) F2=X(2)**(I1-1)
            IF(J.GT.1) F3=X(3)**(J-1)
            F(I+(J-1)*M1**2+(I1-1)*M1)=F1*F2*F3
1000  CONTINUE
      END
 
