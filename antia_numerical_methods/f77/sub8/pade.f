!	To calculate coefficients of Pade approximation using the known
!	coefficients of Maclaurin series
!
!	M : (input) Degree of polynomial in the numerator
!	K : (input) Degree of polynomial in the denominator
!	A : (output) Real array of length M+K+1 containing the coefficients
!		of Pade approximation. A(I) is the coefficient of x**I in
!		the denominator, the constant term being 1. A(K+J+1) is
!		the coefficient of x**J in the numerator
!	C : (input) Real array of length M+K+1 containing the coefficients
!		of Maclaurin series. C(I+1) is the coefficient of x**I
!	IER : (output) error parameter, IER=0 implies successful execution
!		IER=612 implies that M<0 or K<0, no calculations are done
!		Other values of IER may be set by GAUELM
!	WK : Real array of length K*K used as scratch space
!	IWK : Integer array of length K used as scratch space
!
!	Required routines : GAUELM
!
      SUBROUTINE PADE(M,K,A,C,IER,WK,IWK)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(M+K+1),C(M+K+1),IWK(K),WK(K,K)

      IF(M.LT.0.OR.K.LT.0) THEN
        IER=612
        RETURN
      ENDIF

      IF(K.GT.0) THEN
!	Setting up the coefficients matrix of linear equations to
!	calculate the coefficients in the denominator
        DO 2200 I=1,K
          DO 2000 J=1,K
            WK(J,I)=0.0
            IF(M+J-I+1.GT.0) WK(J,I)=C(M+J-I+1)
2000      CONTINUE
!	The right-hand side vector
          A(I)=-C(M+I+1)
2200    CONTINUE

        NUM=1
        LJ=K
        IFLG=0
!	Solve the system of linear equations
        CALL GAUELM(K,NUM,WK,A,DET,IWK,LJ,IER,IFLG)
        IF(IER.GT.0) RETURN
      ENDIF

!	Calculate the coefficients in the numerator
      A(K+1)=C(1)
      IF(M.GT.0) THEN
        DO 3000 I=2,M+1
          AI=C(I)
          IF(K.GT.0) THEN
            DO 2600 J=1,MIN(K,I-1)
2600        AI=AI+C(I-J)*A(J)
          ENDIF
          A(K+I)=AI
3000    CONTINUE
      ENDIF

      END
