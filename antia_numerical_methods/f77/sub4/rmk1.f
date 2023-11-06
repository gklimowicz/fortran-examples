!	To evaluate a rational function at any point, the constant term
!	in denominator is assumed to be unity.
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
!	The value of rational function will be returned as RMK1
!
!	Required routines : None
 
      FUNCTION RMK1(M,K,A,B,X)
!      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(M+1),B(K)
 
      RMK1=A(M+1)
      DO 2000 I=M,1,-1
        RMK1=RMK1*X+A(I)
2000  CONTINUE
 
      DEN=B(K)
      DO 2500 I=K-1,1,-1
        DEN=DEN*X+B(I)
2500  CONTINUE
      DEN=DEN*X+1.0
      RMK1=RMK1/DEN
      END
