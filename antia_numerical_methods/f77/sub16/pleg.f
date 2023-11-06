!	To calculate Legendre polynomial P_l(X)
!
!	L : (input) Order of polynomial (L.GE.0)
!	X : (input) Argument at which the value of polynomial is required
!	P : (output) Real array of length L+1, which will contain the
!		calculated values of polynomials. P(j+1) will contain P_j(X)
!
!	Required routines : None
 
      SUBROUTINE PLEG(L,X,P)
      IMPLICIT REAL*16(A-H,O-Z)
      DIMENSION P(L+1)
 
      IF(L.LT.0) RETURN
      P(1)=1
      P(2)=X
      DO 2000 N=2,L
        P(N+1)=((2*N-1)*X*P(N)-(N-1)*P(N-1))/N
2000  CONTINUE
      END
