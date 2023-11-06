!	To calculate coefficients of rational function Chebyshev approximation
!	using the known coefficients of Chebyshev expansion
!
!	M : (input) Degree of polynomial in the numerator
!	K : (input) Degree of polynomial in the denominator
!	A : (output) Real array of length M+K+1 containing the coefficients
!		of rational function approximation. A(I) is the coefficient
!		of T_I(x) in the denominator, the constant term being 1.
!		A(K+J+1) is the coefficient of T_J(x) in the numerator
!	C : (input) Real array of length M+2K+1 containing the coefficients
!		of Chebyshev expansion. C(I+1) is the coefficient of T_I(x).
!		Coefficient of T_0(x) is C(1)/2
!	IER : (output) error parameter, IER=0 implies successful execution
!		IER=612 implies that M<0 or K<0, no calculations are done
!		Other values of IER may be set by GAUELM
!	WK : Real array of length K*K used as scratch space
!	IWK : Integer array of length K used as scratch space
!
!	Required routines : GAUELM
!
      SUBROUTINE CHEBAP(M,K,A,C,IER,WK,IWK)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(M+K+1),C(M+2*K+1),IWK(K),WK(K,K)

      IF(M.LT.0.OR.K.LT.0) THEN
        IER=612
        RETURN
      ENDIF

      IER=0
      IF(K.GT.0) THEN
!	Setting up the coefficients of matrix of linear equations
!	to calculate the coefficients in the denominator
        DO 2200 I=1,K
          DO 2000 J=1,K
            WK(J,I)=0.5*(C(M+J+I+1)+C(IABS(M+J-I)+1))
2000      CONTINUE
!	The right-hand side vector for linear equations
          A(I)=-C(M+I+1)
2200    CONTINUE

        NUM=1
        LJ=K
        IFLG=0
!	Solve the system of linear equations
        CALL GAUELM(K,NUM,WK,A,DET,IWK,LJ,IER,IFLG)
        IF(IER.GT.0) RETURN
      ENDIF

!	Calculating the coefficients in the numerator
      DO 3000 I=1,M+1
        AI=C(I)
        DO 2600 J=1,K
2600    AI=AI+0.5*A(J)*(C(I+J)+C(IABS(I-1-J)+1))
        A(K+I)=AI
3000  CONTINUE
      A(K+1)=0.5*A(K+1)
      END
