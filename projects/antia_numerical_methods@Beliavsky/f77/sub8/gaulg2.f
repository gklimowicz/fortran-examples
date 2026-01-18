!     To integrate a function with logarithmic singularity over (0,A]
!     using a combination of Gaussian formulas
!
!     RINT : (output) Calculated value of the integral
!     A : (input) The upper limit
!     A1 : (input/output) The point at which integral has to be broken
!		A1 will be adjusted by the subroutine.
!     REPS : (input) The required relative accuracy
!     AEPS : (input) The required absolute accuracy
!		The estimated error should be less than MAX(AEPS,REPS*ABS(RINT))
!     DIF : (output) estimated (absolute) error achieved by the subroutine
!     F : (input) Name of the function routine to calculate the integrand
!     F2 : (input) Name of the function routine to calculate F(X)/LOG(X)
!     NP : (output) Number of function evaluations used by the subroutine
!     IER : (output) Error parameter, IER=0 implies successful execution
!		IER=31 implies specified accuracy was not achieved by GAUSS over [A1,A]
!		IER=32 implies specified accuracy was not achieved by GAULOG
!		IER=34 implies specified accuracy was not achieved by GAUSS over (0,A1]
!		In case of multiple failures second digit of IER will be sum
!			of these values.
!		In all cases DIF will contain the estimated accuracy
!
!     	FUNCTION F(X) and F2(X) must be supplied by the user.
!
!     Required routines : GAUSS, GAULOG, F, F2
 
      SUBROUTINE GAULG2(RINT,A,A1,REPS,AEPS,DIF,F,F2,NP,IER)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(AMN=1.D-2)
      EXTERNAL F,F2
 
      IER1=0
      IER2=0
      R1=0.0
      R2=0.0
      R3=0.0
      A2=0.0
      DIF=0.0
      DIF2=0.0
      NP=0
      NPT=0
      NPT2=0
      IF(A1.GT.A) A1=A
      IF(A1.LE.0.0) A1=A
 
!     Evaluate the integral over (0,A1]
2200  CALL GAULOG(R1,A1,AEPS,REPS,DIF1,F2,NPT1,IER)
      NP=NP+NPT1
      IF(IER.EQ.0) GO TO 2500
!     If GAULOG fails decrease A1
      T1=A1
      A1=A1/2.
      IF(A1.GT.AMN) GO TO 2200
!     To prevent infinite loop do not reduce A1 below AMN
      IER1=2
      IER=0
      A1=T1
 
!     Evaluate the integral over [A1,A]
2500  IF(A-A1.GT.AEPS) CALL GAUSS(R2,A1,A,16,REPS,AEPS,DIF,IER,NPT,F)
      IF(IER.GT.0) IER=1
 
!     Evaluate the regular part over [0,A1]
      IF(A1.NE.1) CALL GAUSS(R3,A2,A1,16,REPS,AEPS,DIF2,IER2,NPT2,F2)
      IF(IER2.GT.0) IER2=4
      IER=IER+IER1+IER2
      IF(IER.GT.0) IER=IER+30
      RINT=R1+R2-R3*LOG(A1)
      DIF=ABS(DIF)+ABS(DIF1)+ABS(DIF2)
      NP=NP+NPT+NPT2
      END
