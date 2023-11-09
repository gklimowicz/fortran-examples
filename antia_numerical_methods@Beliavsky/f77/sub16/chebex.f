!	To calculate the coefficients of Chebyshev expansion of a given
!	function, which can be evaluated at any required point
!	This subroutine uses orthogonality relation over discrete points
!	to compute the coefficients.
!
!	N : (input) Number of coefficients required, this number should
!		be much larger than the actual number required. The
!		accuracy of computed coefficients increases with N.
!	C : (output) Real array of length N containing the required
!		coefficients. The constant term is C(1)/2, while C(I+1)
!		is the coefficient of T_I(x).
!	FUN : (input) Name of function routine to compute the given function
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=613 implies N<10, in which case no calculations are done
!	
!	There is no test for accuracy and this has to be ascertained by
!	doing another calculation with larger N (say 2N).
!
!	FUNCTION FUN(X) must be supplied by the user
!
!	Required routines : FUN
 
      SUBROUTINE CHEBEX(N,C,FUN,IER)
      IMPLICIT REAL*16(A-H,O-Z)
      PARAMETER(PI=3.14159265358979323846264338327950288Q0)
      DIMENSION C(N)
 
      IF(N.LT.10) THEN
        IER=613
        RETURN
      ENDIF
      IER=0
 
      DO 2000 I=1,N
        C(I)=0.0
2000  CONTINUE
 
!     evaluate the sums over n points
      DO 3000 I=1,N
        X=COS((2*I-1)*PI/(2*N))
        FX=FUN(X)
        C(1)=C(1)+FX
        C(2)=C(2)+FX*X
        T0=1.0
        T1=X
        DO 2500 J=3,N
          T2=2*X*T1-T0
          C(J)=C(J)+FX*T2
          T0=T1
          T1=T2
2500    CONTINUE
3000  CONTINUE
 
      DO 3500 I=1,N
        C(I)=2.0*C(I)/N
3500  CONTINUE
      END
