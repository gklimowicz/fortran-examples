!	Least squares polynomial fit using orthogonal polynomials in 1 dimension
!	Modified version of POLFIT to fit multiple sets of function values
!	This routine is meant to be used for fit in multiple dimensions
!
!	N : (input) Number of data points
!	M : (input) Required degree of polynomial
!	NUM : (input) Number of different RHS (function values) to be fitted
!		Each set must be defined over the same abscissas and
!		with same weights.
!	X : (input) Real array of length N containing the abscissas
!	F : (input) Real array of length N*NUM containing the function values
!		for each RHS, F(I,J) is the value at X(I) for Jth set
!		The first dimension of F is assumed to be exactly equal
!		to N to minimise storage requirement.
!	SIG : (input) Real array of length N containing the estimated
!	        error in the function values, SIG(I) is the error in F(I,J).
!		The errors are assumed to be the same for all NUM data sets.
!	A : (output) Real array of length (M+1)*NUM containing the coefficients
!		for the fit for each RHS. The first dimension of A
!		is assumed to be M+1 to minimise storage requirements
!		A(I,J) is the Ith coefficient for Jth RHS.
!	ALP, BETA : (output) Real arrays of length M+1, containing the coefficients
!		required for defining the orthogonal polynomials
!	GAM : (output) Real array of length M+1, containing the quantities
!		\gamma_i for the orthogonal polynomials
!	WK : Real array of length 2N used as scratch space
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=601 implies that N<M+1 or M<0 or N < 1
!		IER=621 implies that GAM(I) vanishes at some I and calculations
!			are abandoned
!	
!	The fitted polynomial can be calculated at any value of x using POLEVL
!
!	Required routines : None
 
      SUBROUTINE POLFIT1(N,M,NUM,X,F,SIG,A,ALP,BETA,GAM,WK,IER)
!      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(N),F(N,NUM),SIG(N),A(M+1,NUM),ALP(M+1),BETA(M+1),
     1      WK(N,2),GAM(M+1)
 
      IF(M.GE.N.OR.M.LT.0.OR.N.LT.1) THEN
        IER=601
        RETURN
      ENDIF
 
!	initialisation
      I1=1
      I2=2
      GAM0=0.0
      DO 2000 I=1,N
        WK(I,I1)=0.0
        WK(I,I2)=1.0
        GAM0=GAM0+1/SIG(i)**2
2000  CONTINUE
      GAM(1)=GAM0
      BETA(1)=0
      IER=0
 
!	Loop over the degree
      DO 5000 J=0,M
        IF(GAM(J+1).LE.0) THEN
          IER=621
          RETURN
        ENDIF
 
        DO 3200 J1=1,NUM
          S=0.0
          DO 3000 I=1,N
3000      S=S+F(I,J1)*WK(I,I2)/SIG(I)**2
          A(J+1,J1)=S/GAM(J+1)
3200    CONTINUE
        IF(J.EQ.M) RETURN

        S=0
        DO 3400 I=1,N
3400    S=S+X(I)*(WK(I,I2)/SIG(I))**2
        ALP(J+1)=S/GAM(J+1)
        GAM0=0.0
        DO 3600 I=1,N
          WK(I,I1)=(X(I)-ALP(J+1))*WK(I,I2)-BETA(J+1)*WK(I,I1)
          GAM0=GAM0+(WK(I,I1)/SIG(I))**2
3600    CONTINUE
        BETA(J+2)=GAM0/GAM(J+1)
        GAM(J+2)=GAM0
        IT=I1
        I1=I2
        I2=IT
5000  CONTINUE
 
      END
