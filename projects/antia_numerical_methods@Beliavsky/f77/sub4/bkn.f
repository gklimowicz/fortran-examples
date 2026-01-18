!	To calculate the modified Bessel function of second kind of
!		integral order for real argument
!		For X.LE.0 the function is not defined and in this case
!		no calculations are done and no warning is issued.
!
!	N : (input) Order of Bessel function required, N must be positive
!	X : (input) Argument at which the value is required
!	BK : (output) Real array of length at least ABS(N)+1
!		which will contain the value of Bessel function of order
!		0,1,...,N. BK(I+1) will contain Bessel function of order I
!
!	Required routines : BK0, BK1
!
 
      SUBROUTINE BKN(N,X,BK)
!      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION BK(N+1)
 
      BK(1)=BK0(X)
      BK(2)=BK1(X)
      IF(X.LE.0.0) RETURN
!	Use the recurrence relation in Forward direction
      DO 2000 I=2,ABS(N)
        BK(I+1)=2*(I-1)*BK(I)/X+BK(I-1)
2000  CONTINUE
      END
