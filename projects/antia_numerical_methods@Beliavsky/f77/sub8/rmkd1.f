!	To evaluate a rational function and its derivative at any point,
!	the constant term in denominator is assumed to be unity.
!
!	M : (input) Degree of numerator
!	K : (input) Degree of denominator
!	A : (input) Real array of length M+1 containing the coefficients
!		of the polynomial in numerator. A(1) is the constant term
!		while A(M+1) is the coefficient of X**M.
!	B : (input) Real array of length K containing the coefficients of
!		the polynomial in denominator. B(1) is the coefficient of X
!		while B(K) is the coefficient of X**K. The constant term
!		is assumed to be unity and is not supplied.
!	X : (input) The value of X at which the rational function needs
!		to be evaluated
!	DF : (output) The first derivative of the rational function at X
!	The value of rational function will be returned as RMKD1
!
!	Required routines : None
 
      FUNCTION RMKD1(M,K,A,B,X,DF)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(M+1),B(K+1)
 
      RMK=A(M+1)
      DN=0.0
      DO 2000 I=M,1,-1
        DN=DN*X+RMK
        RMK=RMK*X+A(I)
2000  CONTINUE
 
      DEN=B(K)
      DEN1=0.0
      DO 2500 I=K-1,1,-1
        DEN1=DEN1*X+DEN
        DEN=DEN*X+B(I)
2500  CONTINUE
      DEN1=DEN1*X+DEN
      DEN=DEN*X+1.0
      RMKD1=RMK/DEN
      DF=DN/DEN-RMKD1*DEN1/DEN
      END
