!	To calculate the error in rational function approximation
!	For use with subroutine REMES. It is called by BRENTM to find
!	extrema of error curve. It may not be used for any other purpose.
!
!	X : (input) the value at which error is to be calculated
!	FM = (FUN(X) - FUND(X)*R_mk(X))*SI  is the calculated error
!	The parameters for rational function are passed through common block
!	M,K are degree of numerator and denominator
!	A is a real array containing the coefficient of rational function
!		approximation
!	SI is -1 for maximum and +1 for minimum. The error is multiplied
!		by SI to make it a minimum.
!		For initial scan SI>10 and function is not evaluated.
!
!	Functions FUN(x) and FUND(x) must be provided by the user to
!		seek approximation of form FUN(X) = FUND(X)*RMK(X)
!
!	Required routines : FUN, FUND
!
      FUNCTION FM(X)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NMAX=50)
      COMMON/ZZFN/A(NMAX),SI,M,K

      NK=M+K+2
      FM=0.0
      IF(SI.LT.10) FM=FUN(X)

!	Calculate the numerator using nested multiplication
      FN=A(K+M+1)
      DO 2200 J=1,M
2200  FN=FN*X+A(NK-J-1)

      IF(K.GT.0) THEN
!	Calculate the denominator using nested multiplication
        FD=A(K)
        DO 2300 J=1,K-1
2300    FD=FD*X+A(K-J)
        FD=FD*X+1
      ELSE
        FD=1.
      ENDIF

      FM=(FM-FUND(X)*FN/FD)*SIGN(1.D0,SI)
      END
