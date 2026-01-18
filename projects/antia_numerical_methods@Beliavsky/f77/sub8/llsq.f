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
      IMPLICIT REAL*8(A-H,O-Z)
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
