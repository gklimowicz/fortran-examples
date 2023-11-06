!	To integrate a function over finite interval using Simpson's rule
!
!	RI : (output) Calculated value of the integral
!	XL : (input) The lower limit
!	XU : (input) The upper limit
!	REPS : (input) The required relative accuracy
!	AEPS : (input) The required absolute accuracy
!		The estimated error should be less than MAX(AEPS,REPS*ABS(RI))
!	DIF : (output) estimated (absolute) error achieved by the subroutine
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=30 implies specified accuracy was not achieved
!			DIF will contain the estimated accuracy
!	N : (output) Number of function evaluations used by subroutine
!	FUN : (input) Name of the function routine to calculate the integrand
!		FUNCTION FUN(X) must be supplied by the user.
!
!	Required routines : FUN

      SUBROUTINE SIMSON(RI,XL,XU,REPS,AEPS,DIF,IER,N,FUN)
      IMPLICIT REAL*16(A-H,O-Z)
      PARAMETER(NMIN=3,NMAX=13)

      FEND=FUN(XL)+FUN(XU)
      EVEN=0.0
      ODD=0.0
      IER=0
      RI=0.0
      DIF=0.0

      N=2
!	starting with 2+1 points, subdivide the intervals into 2 until convergence
      H=(XU-XL)
      IF(H.EQ.0.0) RETURN

      DO 3000 I=1,NMAX
        H=H/2.
        EVEN=EVEN+ODD
        ODD=0.0
        X1=XL+H
        N2=N/2
        H2=2.*H

        DO 1000 J=1,N2
          X=X1+H2*(J-1)
          ODD=ODD+FUN(X)
1000    CONTINUE
!	Estimate for the integral
        R1=(FEND+4.*ODD+2.*EVEN)*H/3.

        DIF=R1-RI
        RI=R1
!	To avoid spurious convergence in first few trials skip the convergence test
        IF(I.LE.NMIN) GO TO 3000
        IF(ABS(DIF).LT.MAX(REPS*ABS(R1),AEPS)) RETURN
3000  N=N*2

      N=N/2
      IER=30
      END
