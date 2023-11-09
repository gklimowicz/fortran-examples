!	To calculate linear least squares fit to B-spline basis functions in
!	1 dimension
!
!	N : (input) Number of data points to be fitted
!	X : (input) Real array of length N containing the coordinates
!		of points at which function values is available
!	F : (input) Real array of length N containing the function values
!		F(I) should be function value at X(I)
!	EF : (input) Real array of length N containing the estimated error
!		in F(I). 
!	K : (input) Order of B-splines required, K=4 gives cubic B-splines
!	A : (output) Real array of length LA*(NO+K-2) containing the matrix
!		U of SVD of the design matrix
!	LA : (input) First dimension of A in the calling program (LA.GE.2N)
!	V : (output) Real array of length IV*(NO+K-2) containing the matrix
!		V of SVD of the design matrix
!	IV : (input) First dimension of V, COV in the calling program (IV.GE.NO+K-2)
!	SIGMA : (output) Real array of length NO+K-2 containing the singular
!		values of the design matrix
!	C : (output) Real array of length 2N containing the fitted coefficients
!		Note that although the number of coefficients is NO+K-2, the
!		rest of array is used as scratch space
!	XF : (input) Real array of size NO, containing
!		the knots used for defining B-spline basis functions.
!		The knots must be distinct and in ascending order.
!	NO : (input) Number of knots for B-splines, the number of basis
!		functions would be NO+K-2
!	Y : (output) Real array of length N containing the values of fitted
!		function at each of the tabular points
!	IFLG : (input/output) Integer specifying the type of calculation required
!		IFLG=0 The matrix will be calculated and solved for coefficients
!			the fitted values Y and CHISQ are also calculated
!		IFLG=1 The matrix will be calculated and SVD
!			is obtained, but coefficients are not calculated
!		IFLG=2 The SVD of matrix is assumed to be available in
!			arrays A, V, SIGMA and coefficients C are calculated
!		IFLG=3 The SVD of matrix is assumed to be available in arrays
!			A, V, SIGMA and coefficients C are calculated and in
!			addition fitted values Y and CHISQ are also calculated
!	WK : Real array of length 4*NO+5*K+2 used as scratch space
!	REPS : (input) Required accuracy for solution of equations using SVD
!		singular values less than REPS times maximum will be set to zero
!	RLM : (input) Parameter lambda for smoothing. If RLM.LE.0 no smoothing
!		is applied
!	IDE : (input) Order of derivative to be used for smoothing
!		This is used only when RLM>0. IDE=1 for first derivative
!		and IDE=2 for second derivative smoothing
!	CHISQ : (output) The value of Chi square at minimum
!       COV : (output) Array of length IV*M containing the covariance
!               matrix of the fitted parameters. COV(I,I) will be the
!               variance in C(I).
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=608 implies that NO+K-2>N or K<2
!		IER=609 implies that RLM>0 and IDE is not acceptable
!		IER=610 implies that EF(I).LE.0 for some I
!		No calculations are done in all these cases
!		Other values of IER may be set by SVD or BSPLIN
!
!       THE ARGUMENTS OF THIS FUNCTION HAVE CHANGED FROM THE EARLIER VERSION.
!	NOW THERE IS AN ADDITIONAL ARGUMENT COV TO CALCULATE THE COVARIANCE
!	MATRIX.
!
!	Required routines : BSPLIN, BSPEVL, SVD, SVDEVL
!
      SUBROUTINE BSPFIT(N,X,F,EF,K,A,LA,V,IV,SIGMA,C,XF,NO,Y,IFLG,WK,
     1       REPS,RLM,IDE,CHISQ,COV,IER)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(N),A(LA,NO+K-2),WK(4*NO+5*K+2),C(N),F(N),XF(NO),EF(N),
     1   Y(N),SIGMA(NO+K-2),V(IV,NO+K-2),COV(IV,NO+K-2)
 
      IF(N.LT.NO+K-2.OR.K.LT.2) THEN
        IER=608
        RETURN
      ENDIF
 
      IF(RLM.GT.0.0.AND.(IDE.LT.1.OR.IDE.GT.2)) THEN
        IER=609
        RETURN
      ENDIF
!	N1 is the number of equations to be solved
      N1=N
      IF(RLM.GT.0.0) N1=2*N
 
      IF(IFLG.LE.1) THEN
!	Set up the matrix equation and obtain SVD
!	M is the number of coefficients to be determined
        M=NO+K-2
        NDB=M+2
        NDERIV=0
        IF(RLM.GT.0.0) THEN
          NDERIV=IDE
          NB=NDB-1
          IF(IDE.EQ.2) NB=2*NDB-1
        ENDIF
 
!	Set up the matrix for equations
        DO 2500 I=1,N
          XB=X(I)
          IF(EF(I).LE.0) THEN
            IER=610
            RETURN
          ENDIF
          CALL BSPLIN(XF,NO,K,XB,NDERIV,WK,WK(NDB),WK(NDB*2),LEFT,IER,
     1                WK(3*NDB))
          IF(IER.GT.100) RETURN
          DO 2200 J=1,M
            A(I,J)=WK(J)/EF(I)
            IF(RLM.GT.0.0) A(I+N,J)=RLM*WK(NB+J)
2200      CONTINUE
2500    CONTINUE
        CALL SVD(M,N1,A,V,SIGMA,LA,IV,WK,IER)
        IF(IER.GT.100) RETURN
 
        IF(IFLG.EQ.1) THEN
          IFLG=2
          RETURN
        ENDIF
 
      ENDIF
 
!	Setup the RHS and solve the equations
      DO 3000 I=1,N
        C(I)=F(I)/EF(I)
        IF(RLM.GT.0.0) C(I+N)=0.0
3000  CONTINUE
 
      CALL SVDEVL(M,N1,A,V,SIGMA,LA,IV,C,WK,REPS)
 
      IF(IFLG.EQ.2) RETURN
      IFLG=2
 
!	Calculate the \chi^2
      CHISQ=0.0
      NDERIV=0
      DO 4000 I=1,N
        Y(I)=BSPEVL(NO,XF,K,NDERIV,C,X(I),DF,DDF,WK,IER)
        CHISQ=CHISQ+((F(I)-Y(I))/EF(I))**2
4000  CONTINUE

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
