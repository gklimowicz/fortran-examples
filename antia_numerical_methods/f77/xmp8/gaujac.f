!     PROGRAM TO CALCULATE WEIGHTS AND ABSCISSAS OF GAUSS-JACOBI QUADRATURE FORMULAS
!     WITH WEIGHT FUNCTION (1-X)**ALP*(1+X)**BETA

      PROGRAM GAUSSL
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(65),W(65),COF(3,65),WK(9000),RIN(0:65)
 
51    FORMAT(5X,'IER =',I4,'    N =',I4,5X,'ALP =',1PD14.6,5X,'BETA =',
     1       D14.6/9X,'ABSCISSAS',10X,'WEIGHTS'/(D22.14,5X,D22.14))
52    FORMAT(/'   MAXIMUM ERROR USING ',I4,' POINTS =',1PD12.4)
 

100   PRINT *,'TYPE  NT=NO. OF POINTS IN QUADRATURE FORMULA, ALP, BETA'
      PRINT *,'           (QUITS WHEN NT.LE.0)'
      READ *, NT,ALP,BETA
      IF(NT.LE.0) STOP

      CALL GAUJAC(NT,ALP,BETA,W,X,IER,WK)
      WRITE(6,51) IER,NT,ALP,BETA,(X(I),W(I),I=1,NT)

!     USE THE WEIGHTS AND ABSCISSAS TO EVALUATE THE INTEGRAL OF X**N FOR TESTING
!     Only N=0,1,2 are used, so the test is not complete
 
      ERR=0.0
      RIN(0)=2.D0**(ALP+BETA+1)*GAMMA(ALP+1.0)*GAMMA(BETA+1.0)/
     1       GAMMA(ALP+BETA+2.0)
      RIN(1)=RIN(0)*(2.0*(1+BETA)/(2+ALP+BETA)-1.0)
      RIN(2)=RIN(0)*(1.0+4*(2+BETA)*(1+BETA)/((3+ALP+BETA)*(2+ALP+BETA))
     1       -4*(1+BETA)/(2+ALP+BETA))
      DO 2000 J=0,2
        RI=0.0
        DO 1000 I=1,NT
          RI=RI+W(I)*X(I)**J
1000    CONTINUE
        REX=RIN(J)
        ERR=MAX(ERR,ABS(RI-REX))
2000  CONTINUE
 
      WRITE(6,52) NT,ERR
      GO TO 100
      END
 
!     --------------------------------------------------
 
!	To calculate Gamma function for any real value of XG
!	Use GAMMAL for calculating the logarithm of Gamma function
!	which may be useful for large arguments or when argument is
!	close to a negative integer.
!
!	Required routines : None
 
      FUNCTION GAMMA(XG)
      IMPLICIT REAL*8(A-H,O-Z)
!	PIS is SQRT(2*PI)
      PARAMETER(PI=3.14159265358979323846D0,PIS=2.5066282746310005024D0)
      DIMENSION A(2),B(3),A1(6),B1(7)

!	The coefficients for rational function approximations
      DATA A/1.767971449569122937D+00,  2.909421117928672645D-01/
      DATA B/8.333333333333231537D-02,  1.445531763554246280D-01,
     1       2.012779361583001035D-02/
      DATA A1/3.905731686764559737D+03,  2.204952264401381785D+03,
     1       -1.932467485468849660D+03,  4.643360871045442213D+02,
     1       -4.818088806916028754D+01,  1.896853765546068169D+00/
      DATA B1/3.918055655523400310D+03, -1.088116266563809683D+02,
     1        8.203258626193993149D+02, -9.289402000761705906D+01,
     1        6.521113026294866877D+01, -6.090618615608719044D+00,
     1        1.475909104740280784D+00/
 
      X=ABS(XG)
      IF(X.GT.1000.0) THEN
!	Use asymptotic approximation (Stirling formula)
        GX=(1+1.D0/(12*X)+1./(288*X*X)
     1         -139/(51840*X**3)-571./(2488320D0*X**4))
        GAMMA=X**(X-0.5)*EXP(-X)*PIS*GX

      ELSE IF(X.GT.8.0) THEN
!	Use rational function approximation for Log(Gamma) 
        Y=1./X**2
        RMK=((B(3)*Y+B(2))*Y+B(1))/((A(2)*Y+A(1))*Y+1)
        GAMMA=X**(X-0.5)*EXP(-X)*PIS*EXP(RMK/X)

      ELSE IF(X.GE.2.0) THEN
!	Use rational function approximation for (Gamma) over [2,3]
!	after translating the range if necessary
        F1=1.0
        X1=X
2500    IF(X1.LE.3) GO TO 3000
        F1=F1*(X1-1)
        X1=X1-1
        GO TO 2500
3000    IF(X1.EQ.3) THEN
          GAMMA=F1*2
        ELSE IF(X1.EQ.2) THEN
          GAMMA=F1
        ENDIF

        FN=(((((B1(7)*X1+B1(6))*X1+B1(5))*X1+B1(4))*X1+B1(3))*X1+
     1          B1(2))*X1+B1(1)
        FD=(((((A1(6)*X1+A1(5))*X1+A1(4))*X1+A1(3))*X1+A1(2))*X1+
     1          A1(1))*X1+1
        GAMMA=F1*FN/FD
      ELSE IF(X.GT.0.0) THEN
!	Use rational function approximation for (Gamma) over [2,3]
!	after translating the range if necessary
        F1=1./X
        X1=X+1
        IF(X.LT.1) THEN
          F1=F1/X1
          X1=X1+1
        ENDIF
        IF(X1.EQ.2) GAMMA=F1

        FN=(((((B1(7)*X1+B1(6))*X1+B1(5))*X1+B1(4))*X1+B1(3))*X1+
     1          B1(2))*X1+B1(1)
        FD=(((((A1(6)*X1+A1(5))*X1+A1(4))*X1+A1(3))*X1+A1(2))*X1+
     1          A1(1))*X1+1
        GAMMA=F1*FN/FD
 
      ENDIF

      IF(XG.GT.0.0) RETURN
      IX=X
      IF(X.GT.IX) THEN
        GAMMA=PI/(XG*SIN(PI*X)*GAMMA)
      ELSE
        GAMMA=(-1)**IX/0.0
      ENDIF
 
      END
 
!     --------------------------------------------------
 
!     To calculate weights and abscissas of a Gauss-Jacobi quadrature formula
!
!     N : (input) Number of points in the required quadrature formula
!     ALP : (input) Exponent of 1-x in weight function
!     BETA : (input) Exponent of 1+x in weight function
!		w(x)=(1-x)**ALP*(1+x)**BETA
!     W : (output) Array of length N, which will contain the weights
!     A : (output) Array of length N containing the abscissas
!     IER : (output) Error parameter, IER=0 implies successful execution
!		IER=302 implies N.LE.0
!		IER=313 implies ALP.LE.-1 or BETA.LE.-1
!		IER=321 implies that some coefficient becomes imaginary
!			during calculations.
!		In both these cases calculations are abandoned.
!		Other values may be set by TQL2
!     WK : Real array of length N*(N+2)+3*(N+1) used as scratch space
!
!     Required routines : GAUSRC, TQL2, GAMMA
 
      SUBROUTINE GAUJAC(N,ALP,BETA,W,A,IER,WK)
      IMPLICIT REAL*8(A-H,O,P,R-Z)
      DIMENSION A(N),W(N),WK(3,*)
 
      IF(ALP.LE.-1.0.OR.BETA.LE.-1.0) THEN
        IER=313
        RETURN
      ENDIF

      RI0=2.D0**(ALP+BETA+1)*GAMMA(ALP+1.0)*GAMMA(BETA+1.0)/
     1       GAMMA(ALP+BETA+2.0)

!     Coefficients of recurrence relation
      WK(1,1)=1.D0+(ALP+BETA)/2.0
      WK(2,1)=(ALP-BETA)/2.0
!	WK(3,1) is not required and would give 0/0 form in some cases
!	Hence do not calculate it.

      DO 2500 I=2,N+1
        A1=2.0*I*(I+ALP+BETA)*(2*I-2+ALP+BETA)
        WK(1,I)=(2*I-2+ALP+BETA)*(ALP+BETA+2*I-1)*(ALP+BETA+2*I)/A1
        WK(2,I)=(2*I-1.0+ALP+BETA)*(ALP*ALP-BETA*BETA)/A1
        WK(3,I)=2.0*(I-1.D0+ALP)*(I-1.0+BETA)*(ALP+BETA+2*I)/A1
2500  CONTINUE
 
      CALL GAUSRC(N,W,A,WK,RI0,IER,WK(1,N+2))
      END
 
!     --------------------------------------------------
 
!     To calculate weights and abscissas of a quadrature formula with
!     specified weight function when recurrence relation for orthogonal
!     polynomial is known.
!
!     N : (input) Number of points in the required quadrature formula
!     W : (output) Array of length N, which will contain the weights
!     AB : (output) Array of length N containing the abscissas
!     COF : (input) Array of length 3*N containing the coefficients of
!		the recurrence relation for orthogonal polynomials
!		P_i(x)=(COF(1,i)*x+COF(2,i))P_{i-1}(x) - COF(3,i)*P_{i-2}(x)
!     RI0 : (input) The integral of weight function over the required
!		interval.
!     IER : (output) Error parameter, IER=0 implies successful execution
!		IER=302 implies N.LE.0
!		IER=321 implies that some coefficient becomes imaginary
!			during calculations.
!		In both these cases calculations are abandoned.
!		Other values may be set by TQL2
!     WK : Real array of length N*(N+2) used as scratch space
!
!     Required routines : TQL2
 
      SUBROUTINE GAUSRC(N,W,AB,COF,RI0,IER,WK)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (REPS=1.D-15)
!	For REAL*4 use REPS=6.E-8
!      PARAMETER (REPS=6.E-8)
      DIMENSION WK(N,N+2),COF(3,N),W(N),AB(N)
 
      IER=302
      IF(N.LE.0) RETURN
      LJ=N
 
!     Calculate the coefficients of symmetric tridiagonal matrix
      DO 2000 I=1,N
        WK(I,N+1)=-COF(2,I)/COF(1,I)
        IF(I.LT.N) THEN
          R1=COF(3,I+1)/(COF(1,I)*COF(1,I+1))
          IF(R1.GE.0.0) THEN
            WK(I+1,N+2)=SQRT(R1)
          ELSE
            IER=321
            RETURN
          ENDIF
        ENDIF
        DO 1800 J=1,N
          WK(J,I)=0.0
1800    CONTINUE
        WK(I,I)=1.0
2000  CONTINUE
 
!     Find eigenvalues and eigenvectors of the tridiagonal matrix
      CALL TQL2(WK,N,LJ,WK(1,N+1),WK(1,N+2),REPS,IER)
      IF(IER.GT.0) RETURN
 
!     Calculate the abscissas and weights
      DO 3000 I=1,N
        AB(I)=WK(I,N+1)
        W(I)=WK(1,I)**2*RI0
3000  CONTINUE
 
      END

!	------------------------------------------------------------------

!	To find eigenvalues and eigenvectors of Z T Z^T using QL algorithm
!	where T is a symmetric tridiagonal matrix and Z is an orthogonal matrix.
!	If Z is the transformation matrix to reduce original matrix to
!	tridiagonal matrix, it will calculate the eigenvectors of original matrix
!
!	Z : (input/output) Real array of length IZ*N which should contain
!		the transformation matrix required to reduce original real
!		symmetric matrix to tridiagonal form. To find eigenvectors
!		of a symmetric tridiagonal matrix, set Z to unit matrix
!		After execution Z will contain the eigenvector of the original
!		matrix Z T Z^T. Z(i,j) should contain the ith component of
!		jth eigenvector
!	N : (input) Order of matrix
!	IZ : (input) The first dimension of array Z as declared in the
!		calling program. (IZ.GE.N)
!	D : (input/output) Real array of length N, containing the diagonal
!		elements of the tridiagonal matrix, D(i)=T(i,i).
!		After execution it will contain the eigenvalues of the matrix
!	E : (input/output) Real array of length N containing the off-diagonal
!		elements of the tridiagonal matrix, E(i)=T(i,i+1)=T(i+1,i)
!		It is used as scratch space and its contents will be destroyed
!		during execution.
!	REPS : (input) Required tolerance, it should be of order of machine
!		accuracy
!	IER : (output) Error parameter; IER=0 implies successful execution
!		IER=108 implies that N.LE.1 or N>IZ, in which case no
!			calculations are performed
!		IER=143 implies that the QL algorithm failed to converge
!			for some eigenvalue, the calculations are abandoned
!
!	Required routines : None
!	
      SUBROUTINE TQL2(Z,N,IZ,D,E,REPS,IER)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NIT=30)
      DIMENSION Z(IZ,N),D(N),E(N)

      IF(N.LE.0.OR.N.GT.IZ) THEN
        IER=108
        RETURN
      ENDIF
      IER=0
      DO 2000 I=2,N
2000  E(I-1)=E(I)
      E(N)=0
      B=0
      F=0

      DO 4000 L=1,N
        H=REPS*(ABS(D(L))+ABS(E(L)))
        IF(B.LT.H) B=H
!	Look for small off-diagonal elements
        DO 2200 M=L,N
          IF(ABS(E(M)).LE.B) GO TO 2300
2200    CONTINUE
        M=N
2300    IF(M.EQ.L) THEN
!	one eigenvalue is isolated
          D(L)=D(L)+F
          GO TO 4000
        ENDIF

!	Loop for QL transformation 
        DO 3600 IT=1,NIT
!	Find shift
          G=D(L)
          P=(D(L+1)-G)/(2*E(L))
          R=SQRT(P*P+1)
          IF(P.LT.0.0) R=-R
          D(L)=E(L)/(P+R)
          H=G-D(L)
          DO 2600 I=L+1,N
2600      D(I)=D(I)-H
          F=F+H

!	The QL transformation
          P=D(M)
          C=1
          S=0
!	Given's rotations
          DO 3400 I=M-1,L,-1
            G=C*E(I)
            H=C*P
            IF(ABS(P).GE.ABS(E(I))) THEN
              C=E(I)/P
              R=SQRT(C*C+1)
              E(I+1)=S*P*R
              S=C/R
              C=1/R
            ELSE
              C=P/E(I)
              R=SQRT(C*C+1)
              E(I+1)=S*E(I)*R
              S=1./R
              C=C/R
            ENDIF
            P=C*D(I)-S*G
            D(I+1)=H+S*(C*G+S*D(I))

!	Transforming the eigenvectors
            DO 3000 K=1,N
              H=Z(K,I+1)
              Z(K,I+1)=S*Z(K,I)+C*H
              Z(K,I)=C*Z(K,I)-S*H
3000        CONTINUE
3400      CONTINUE

          E(L)=S*P
          D(L)=C*P
          IF(ABS(E(L)).LE.B) THEN
!	One eigenvalue is isolated
            D(L)=D(L)+F
            GO TO 4000
          ENDIF
3600    CONTINUE
!	QL iteration fails to converge
        IER=143
        RETURN
4000  CONTINUE

!	Sort eigenvalues in ascending order by straight selection
      DO 5000 I=1,N-1
        K=I
        P=D(I)
        DO 4200 J=I+1,N
          IF(D(J).LT.P) THEN
            K=J
            P=D(J)
          ENDIF
4200    CONTINUE
        IF(K.NE.I) THEN
!	exchange the eigenvalues and eigenvectors
          D(K)=D(I)
          D(I)=P
          DO 4400 J=1,N
            P=Z(J,I)
            Z(J,I)=Z(J,K)
            Z(J,K)=P
4400      CONTINUE
        ENDIF
5000  CONTINUE
      END
