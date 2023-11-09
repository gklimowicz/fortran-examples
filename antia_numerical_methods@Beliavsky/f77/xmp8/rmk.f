!     TO EVALUATE A RATIONAL FUNCTION AND ITS DERIVATIVE

      PROGRAM RATNAL
      IMPLICIT REAL*8(A-H,O-Z)
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
      READ *,(B(I),I=K+1,1,-1)
      WRITE(6,54) (B(I),I=K+1,1,-1)
 
100   PRINT *,'TYPE X     (QUITS WHEN X<-100)'
      READ *,X
      IF(X.LT.-100) STOP
      F1=RMK(M,K,A,B,X)
      WRITE(6,51) M,K,X,F1
      F2=RMKD(M,K,A,B,X,DF)
      WRITE(6,52) F2,DF
      GO TO 100
      END
 
!     ------------------------------------------------------------
 
!	To evaluate a rational function at any point
!
!	M : (input) Degree of numerator
!	K : (input) Degree of denominator
!	A : (input) Real array of length M+1 containing the coefficients
!		of the polynomial in numerator. A(1) is the constant term
!		while A(M+1) is the coefficient of X**M.
!	B : (input) Real array of length K+1 containing the coefficients
!		of the polynomial in denominator. B(1) is the constant term
!		while B(K+1) is the coefficient of X**K.
!	X : (input) The value of X at which the rational function needs
!		to be evaluated
!	The value of rational function will be returned as RMK
!
!	Required routines : None
 
      FUNCTION RMK(M,K,A,B,X)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(M+1),B(K+1)
 
      RMK=A(M+1)
      DO 2000 I=M,1,-1
        RMK=RMK*X+A(I)
2000  CONTINUE
 
      DEN=B(K+1)
      DO 2500 I=K,1,-1
        DEN=DEN*X+B(I)
2500  CONTINUE
      RMK=RMK/DEN
      END
 
!     ------------------------------------------------------------
 
!	To evaluate a rational function and its derivative at any point
!
!	M : (input) Degree of numerator
!	K : (input) Degree of denominator
!	A : (input) Real array of length M+1 containing the coefficients
!		of the polynomial in numerator. A(1) is the constant term
!		while A(M+1) is the coefficient of X**M.
!	B : (input) Real array of length K+1 containing the coefficients
!		of the polynomial in denominator. B(1) is the constant term
!		while B(K+1) is the coefficient of X**K.
!	X : (input) The value of X at which the rational function needs
!		to be evaluated
!	DF : (output) The first derivative of the rational function at X
!	The value of rational function will be returned as RMKD
!
!	Required routines : None
 
      FUNCTION RMKD(M,K,A,B,X,DF)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(M+1),B(K+1)
 
      RMK=A(M+1)
      DN=0.0
      DO 2000 I=M,1,-1
        DN=DN*X+RMK
        RMK=RMK*X+A(I)
2000  CONTINUE
 
      DEN=B(K+1)
      DEN1=0.0
      DO 2500 I=K,1,-1
        DEN1=DEN1*X+DEN
        DEN=DEN*X+B(I)
2500  CONTINUE
      RMKD=RMK/DEN
      DF=DN/DEN-RMKD*DEN1/DEN
      END
