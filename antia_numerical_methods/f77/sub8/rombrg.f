!	To integrate a function over finite interval using Romberg integration
!
!	RI : (output) Calculated value of the integral
!	A : (input) The lower limit
!	B : (input) The upper limit
!	GI : (input/output) Real array of length NMAX (=13), containing
!		the expected values of exponents \gamma_i in error expansion
!		If GI(I).LE.0 it will be set to 2I, the correct value for
!		a smooth function
!	REPS : (input) The required relative accuracy
!	AEPS : (input) The required absolute accuracy
!		The estimated error should be less than MAX(AEPS,REPS*ABS(RI))
!	DIF : (output) estimated (absolute) error achieved by the subroutine
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=30 implies specified accuracy was not achieved
!			DIF will contain the estimated accuracy
!		IER=33 implies that N>NPT (=100) in which case it is set to 2
!	N : (input/output) On input it should contain the number of function evaluations
!		to be used for first estimate. If N<2 or N>NPT it is set to 2.
!		After execution it will contain the number of function
!		evaluations actually used by subroutine
!	FUN : (input) Name of the function routine to calculate the integrand
!		FUNCTION FUN(X) must be supplied by the user.
!
!	Required routines : FUN

      SUBROUTINE ROMBRG(RI,A,B,GI,REPS,AEPS,DIF,IER,N,FUN)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NMIN=3,NMAX=13,NPT=100)
      DIMENSION GI(NMAX),T(NMAX,NMAX)

      DO 1000 I=1,NMAX
        IF(GI(I).LE.0.0) GI(I)=2*I
1000  CONTINUE

      IER=0
      IF(N.LE.1) N=2
      IF(N.GT.NPT) THEN
        IER=33
        N=2
      ENDIF
!	Contribution from the end points
      S1=0.5*(FUN(A)+FUN(B))
!	First time use all points
      ND=1
      DIF=0.0

      DO 4000 I=1,NMAX
        H=(B-A)/(N-1)

!	Add new points to the sum
        DO 2200 J=2,N-1,ND
          Y=A+(J-1)*H
2200    S1=S1+FUN(Y)
!	The trapezoidal rule approximation
        T(I,1)=S1*H

!	The Richardson's extrapolation
        DO 2400 J=1,I-1
          FJ=2.**GI(J)
          T(I,J+1)=T(I,J)+(T(I,J)-T(I-1,J))/(FJ-1)
          DIF1=ABS(T(I,J)-T(I-1,J))
!	Find the minimum difference between the last two rows of T-table
          IF(DIF1.LT.DIF.OR.J.EQ.1) THEN
            DIF=DIF1
            RI=T(I,J)
          ENDIF
2400    CONTINUE

!	On second and subsequent pass add only new points to the sum
        ND=2
        IF(I.LE.NMIN) GO TO 4000
        IF(DIF.LT.MAX(REPS*ABS(RI),AEPS)) RETURN

4000  N=2*N-1

!	Routine fails to converge
      IER=30
      N=(N+1)/2
      END
