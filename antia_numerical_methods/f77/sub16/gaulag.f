!	To integrate a function over semi-infinite interval using a
!	combination of Gauss-Legendre and Gauss-Laguerre formulas
!
!	RINT : (output) Calculated value of the integral
!	A : (input) The lower limit
!	A1 : (input/output) The point at which integral has to be broken
!		A1 will be adjusted by the subroutine.
!	REPS : (input) The required relative accuracy
!	AEPS : (input) The required absolute accuracy
!		The estimated error should be less than MAX(AEPS,REPS*ABS(RINT))
!	DIF : (output) estimated (absolute) error achieved by the subroutine
!	F : (input) Name of the function routine to calculate the integrand
!	F2 : (input) Name of the function routine to calculate F(X)*EXP(X)
!	NP : (output) Number of function evaluations used by the subroutine
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=30 implies specified accuracy was not achieved by GAUSS
!		IER=37 implies specified accuracy was not achieved by LAGURE
!		IER=38 implies specified accuracy was not achieved by both
!		GAUSS and LAGURE
!		In all cases DIF will contain the estimated accuracy
!
!		FUNCTION F(X) and F2(X) must be supplied by the user.
!
!	Required routines : GAUSS, LAGURE, F, F2
!
      SUBROUTINE GAULAG(RINT,A,A1,REPS,AEPS,DIF,F,F2,NP,IER)
      IMPLICIT REAL*16(A-H,O-Z)
      PARAMETER(AMAX=50.)
      EXTERNAL F,F2

      IER1=0
      R1=0.0
      R2=0.0
      DIF=0.0
      NP=0
      NPT=0
      IF(A1.LT.A) A1=A

!	To calculate integral over [A1,Infinity)
2200  CALL LAGURE(R1,A1,AEPS,REPS,DIF1,F2,NPT1,IER)
      NP=NP+NPT1
      IF(IER.EQ.0) GO TO 2500
!	If LAGURE fails then increase A1
      T1=A1
      A1=MAX(A1*2.,A1+2.)
      IF(A1.LT.AMAX) GO TO 2200
!	To avoid the possibility of getting in infinite loop A1 is not
!	increased beyond AMAX
      IER1=37
      IER=0
      A1=T1

!	To calculate integral over [A,A1]
2500  IF(A1-A.GT.AEPS) CALL GAUSS(R2,A,A1,16,REPS,AEPS,DIF,IER,NPT,F)
      IER=IER+IER1
      IF(IER.GT.IER1.AND.IER1.GT.0) IER=38
      RINT=R1+R2
      DIF=ABS(DIF)+ABS(DIF1)
      NP=NP+NPT
      END
