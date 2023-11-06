!	To calculate the Bessel function of second kind of integral order
!		for real argument
!
!	N : (input) Order of Bessel function required, N may be negative
!		or positive. For N=0,1 use BY0 and BY1 respectively.
!	X : (input) Argument at which the value is required
!	BY : (output) Real array of length at least ABS(N)+1
!		which will contain the value of Bessel function of order
!		0,1,...,ABS(N). BY(I+1) will contain Bessel function of order I
!		or -I (if N<0).
!
!	Required routines : BY0, BY1
!
 
      SUBROUTINE BYN(N,X,BY)
      IMPLICIT REAL*16(A-H,O-Z)
      DIMENSION BY(N+1)
 
      BY(1)=BY0(X)
      BY(2)=BY1(X)
      IF(X.LE.0.0) RETURN
!	Use the recurrence relation in Forward direction
      DO 2000 I=2,ABS(N)
        BY(I+1)=2*(I-1)*BY(I)/X-BY(I-1)
2000  CONTINUE
 
      IF(N.LT.0) THEN
        DO 2500 I=2,ABS(N)+1,2
          BY(I)=-BY(I)
2500    CONTINUE
      ENDIF
      END
