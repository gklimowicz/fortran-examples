!     TO EVALUATE A RATIONAL FUNCTION AND ITS DERIVATIVE
!     CONSTANT TERM IN DENOMINATOR IS ASSUMED TO BE 1

      PROGRAM RATNAL
      IMPLICIT REAL*16(A-H,O-Z)
      DIMENSION A(100),B(100)
 
51    FORMAT('   M =',I4,4X,'K =',I3,5X,'X =',1PD14.6,
     1       3X,'RMK(X) =',D14.6)
52    FORMAT('   RMK(X) =',1PD14.6,5X,9HRMK'(X) =,D14.6)
53    FORMAT('   COEFFICIENTS IN NUMERATOR :',1P3D14.6/(5D14.6))
54    FORMAT('   COEFFICIENTS IN DENOMINATOR :',1P3D14.6/(5D14.6))

      PRINT *,'TYPE M, K = DEGREE OF POLYNOMIAL IN NUMERATOR AND'
     1        ,' DENOMINATOR'
      READ *,M,K
	  PRINT *,'TYPE COEFFICIENTS IN NUMERATOR, STARTING FROM HIGHEST'
     1      ,' DEGREE TERM'
      READ *,(A(I),I=M+1,1,-1)
      WRITE(6,53) (A(I),I=M+1,1,-1)
	  PRINT *,'TYPE COEFFICIENTS IN DENOMINATOR, STARTING FROM HIGHEST'
     1      ,' DEGREE TERM'
      READ *,(B(I),I=K,1,-1)
      WRITE(6,54) (B(I),I=K,1,-1),1.0
 
100   PRINT *,'TYPE X     (QUITS WHEN X<-100)'
      READ *,X
      IF(X.LT.-100) STOP
      F1=RMK1(M,K,A,B,X)
      WRITE(6,51) M,K,X,F1
      F2=RMKD1(M,K,A,B,X,DF)
      WRITE(6,52) F2,DF
      GO TO 100
      END
 
!     ----------------------------------------------------------
 
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
      IMPLICIT REAL*16(A-H,O-Z)
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
 
!     ----------------------------------------------------------
 
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
      IMPLICIT REAL*16(A-H,O-Z)
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
