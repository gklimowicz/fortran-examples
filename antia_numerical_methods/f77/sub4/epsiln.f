!	To integrate a function over finite interval using Epsilon algorithm
!
!	RI : (output) Calculated value of the integral
!	A : (input) The lower limit
!	B : (input) The upper limit
!	REPS : (input) The required relative accuracy
!	AEPS : (input) The required absolute accuracy
!		The estimated error should be less than MAX(AEPS,REPS*ABS(RI))
!	DIF : (output) estimated (absolute) error achieved by the subroutine
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=30 implies specified accuracy was not achieved
!			DIF will contain the estimated accuracy
!		IER=33 implies that N>NPT (=100) in which case it is set to 2
!		IER=34 implies that at some stage denominator vanished while
!			calculating epsilon table. This value is ignored.
!		IER=35 implies that roundoff error appears to be dominating
!	N : (input/output) On input it should contain the number of function
!		evaluations to be used for first estimate. If N<2 or N>NPT it
!		is set to 2. After execution it will contain the number of
!		function evaluations actually used by subroutine
!	FUN : (input) Name of the function routine to calculate the integrand
!		FUNCTION FUN(X) must be supplied by the user.
!
!	Required routines : FUN

      SUBROUTINE EPSILN(RI,A,B,REPS,AEPS,DIF,IER,N,FUN)
!      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NMIN=4,NMAX=13,NPT=100)
      DIMENSION T(NMAX,NMAX)

      IER=0
      IF(N.LE.1) N=2
      IF(N.GT.NPT) THEN
        IER=33
        N=2
      ENDIF

!	Sum of the end points for trapezoidal rule
      S1=0.5*(FUN(A)+FUN(B))
      ND=1
      RI=0.0
      T(1,1)=0.0

      DO 4000 I=1,NMAX-1
        H=(B-A)/(N-1)

        DO 2200 J=2,N-1,ND
          Y=A+(J-1)*H
2200    S1=S1+FUN(Y)
!	The trapezoidal rule approximation
        T(I,2)=S1*H
        T(I+1,1)=0.0
        RI1=RI
        IF(I.GE.2) THEN
          DIF=ABS(T(I,2)-T(I-1,2))
          RI=T(I,2)
        ENDIF

!	Construct the Epsilon table
        DO 2400 J=3,I+1
          DEN=T(I-J+3,J-1)-T(I-J+2,J-1)

!	If denominator is zero set the error flag
          IF(DEN.NE.0.0) THEN
            T(I-J+2,J)=T(I-J+3,J-2)+1./DEN
          ELSE
            IER=34
            T(I-J+2,J)=T(I-J+3,J-2)
          ENDIF

2400    CONTINUE

        IF(I.GT.4) THEN
!	DIF is the minimum difference between two rows of epsilon table
          DO 2600 J=4,I-1,2
            DIF1=ABS(T(I-J+2,J)-T(I-J+1,J))
            IF(DIF1.LT.DIF) THEN
              DIF=DIF1
              RI=T(I-J+2,J)
            ENDIF
2600      CONTINUE
        ENDIF

        ND=2
        IF(I.LE.NMIN) GO TO 4000
        IF(I.GT.6.AND.DIF.GT.DIF0) THEN
!	Roundoff error appears to be dominating, retain the previous value of RI
          IER=35
          RI=RI1
          RETURN
        ENDIF
        DIF0=DIF
        IF(DIF.LT.MAX(REPS*ABS(RI),AEPS)) RETURN

4000  N=2*N-1

!	Integral fails to converge
      IER=30
      N=(N+1)/2
      END
