!	To calculate the incomplete Beta function I_x(a,b) using
!	the integral, called by BETAP
!
!	A,B : (input) Arguments for the complete Beta function
!	X : (input) Upper limit of integration defining the incomplete
!                       Beta function
!
!	Required routines : GAMMAL, ADPINT, KRONRD, FBETA

      FUNCTION BETAI(A,B,X)
      IMPLICIT REAL*16(A-H,O-Z)
      EXTERNAL FBETA
      COMMON/BETA/AA,BB

      AA=A
      BB=B
      XL=0
      XU=X
      REPS=1.Q-14
      AEPS=1.Q-280
      CALL ADPINT(RINT,XL,XU,REPS,AEPS,DIF,FBETA,IER,NPT,NMAX)
      B1=LOG(RINT)+GAMMAL(A+B)-GAMMAL(A)-GAMMAL(B)
      BETAI=EXP(B1)
      IF(BETAI.GT.1.0) BETAI=1.0
      END
