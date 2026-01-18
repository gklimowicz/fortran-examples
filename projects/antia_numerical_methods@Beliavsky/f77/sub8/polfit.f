!	Least squares polynomial fit using orthogonal polynomials in 1 dimension
!
!	N : (input) Number of data points
!	M : (input) Required degree of polynomial
!	X : (input) Real array of length N containing the abscissas
!	F : (input) Real array of length N containing the function values
!		F(I) is the value at X(I)
!	SIG : (input) Real array of length N containing the estimated
!	        error in the function values, SIG(I) is the error in F(I)
!	A : (output) Real array of length M+1 containing the coefficients
!		for the fit
!	ALP, BETA : (output) Real arrays of length M+1, containing the coefficients
!		required for defining the orthogonal polynomials
!	Y : (output) Real array of length N containing the fitted values at X(I)
!	H : (output) Real array of length M+1 containing the \chi^2 values
!		for residuals using polynomial fit of degrees 0,..,M
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
!       THE ARGUMENTS OF THIS FUNCTION HAS CHANGED AS COMPARED TO EARLIER
!	VERSION AS THE 5TH ARGUMENT IS NOW ERROR INSTEAD OF WEIGHT
!
!	Required routines : None

      SUBROUTINE POLFIT(N,M,X,F,SIG,A,ALP,BETA,Y,H,GAM,WK,IER)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(N),F(N),SIG(N),A(M+1),ALP(M+1),BETA(M+1),Y(N),WK(N,2),
     1          H(M+1),GAM(M+1)

      IF(M.GE.N.OR.M.LT.0.OR.N.LT.1) THEN
        IER=601
        RETURN
      ENDIF

!	Initialisation
      I1=1
      I2=2
      GAM0=0.0
      H0=0.0
      DO 2000 I=1,N
        WK(I,I1)=0.0
        WK(I,I2)=1.0
        Y(I)=0.0
        GAM0=GAM0+1/SIG(I)**2
        H0=H0+(F(I)/SIG(I))**2
2000  CONTINUE
      GAM(1)=GAM0
      BETA(1)=0
      IER=0

!	Loop over the degree of polynomial
      DO 5000 J=0,M
        S=0.0
        DO 3000 I=1,N
3000    S=S+F(I)*WK(I,I2)/SIG(I)**2
        IF(GAM(J+1).LE.0) THEN
          IER=621
          RETURN
        ENDIF

!	The coefficient a_j
        A(J+1)=S/GAM(J+1)
        H0=H0-A(J+1)**2*GAM(J+1)
        H(J+1)=H0
        DO 3200 I=1,N
3200    Y(I)=Y(I)+A(J+1)*WK(I,I2)
        IF(J.EQ.M) RETURN
        S=0
        DO 3400 I=1,N
3400    S=S+X(I)*(WK(I,I2)/SIG(I))**2
!	The coefficient \alpha_{j+1}
        ALP(J+1)=S/GAM(J+1)
        GAM0=0.0
        DO 3600 I=1,N
          WK(I,I1)=(X(I)-ALP(J+1))*WK(I,I2)-BETA(J+1)*WK(I,I1)
          GAM0=GAM0+(WK(I,I1)/SIG(I))**2
3600    CONTINUE
!	The coefficient \beta_{j+1}
        BETA(J+2)=GAM0/GAM(J+1)
!	The coefficient \gamma_{j+1}
        GAM(J+2)=GAM0
!	Interchange indices I1, I2 so that only last two columns of WK are stored
        IT=I1
        I1=I2
        I2=IT
5000  CONTINUE

      END
