!     PROGRAM FOR INTEGRATING F(X)/SQRT((X-A)(B-X)) OVER (A,B)

      PROGRAM INTEG
      IMPLICIT REAL*16(A-H,O-Z)
      EXTERNAL F
!		Pass upper limit B through common block to F(x)
      COMMON B

!     EXERCISE 6.23 (I_4)

51    FORMAT('   IER =',I4,5X,'A =',1PD14.6,5X,'B=',D14.6/5X,
     1  'NO. OF FUNCTION EVALUATIONS =',I7/5X,'INTEGRAL =',D14.6,5X,
     2       'ESTIMATED ERROR =',3D14.6)

      REPS=1.Q-33
      AEPS=1.Q-39

100   PRINT *,'TYPE  A=LOWER LIMIT,  B=UPPER LIMIT   QUITS WHEN(A.EQ.B)'
      READ *,A,B
      IF(A.EQ.B) STOP
      CALL GAUCBY(RI,A,B,REPS,AEPS,DIF,IER,N,F)
      WRITE(6,51) IER,A,B,N,RI,DIF
      GO TO 100
      END
 
!     ---------------------------------------------

!	To integrate a function over finite interval using Gauss-Chebyshev formulas
!	Calculates the integral of FUN(X)/SQRT((X-A)*(B-X))
!
!	RINT : (output) Calculated value of the integral
!	A : (input) The lower limit
!	B : (input) The upper limit
!	REPS : (input) The required relative accuracy
!	AEPS : (input) The required absolute accuracy
!		The estimated error should be less than MAX(AEPS,REPS*ABS(RINT))
!	DIF : (output) estimated (absolute) error achieved by the subroutine
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=30 implies specified accuracy was not achieved
!			DIF will contain the estimated accuracy
!	NPT : (output) Number of function evaluations used by the subroutine
!	FUN : (input) Name of the function routine to calculate the
!		integrand multiplied by SQRT((X-A)*(B-X))
!
!		FUNCTION FUN(X) must be supplied by the user.
!
!	Required routines : FUN

      SUBROUTINE GAUCBY(RINT,A,B,REPS,AEPS,DIF,IER,NPT,FUN)
      IMPLICIT REAL*16(A-H,O-Z)
      PARAMETER(NMAX=13,PI=3.14159265358979323846264338327950288Q0)
 
 
      N=1
      IER=0
      RINT=0.0
      NPT=0
      A0=(B+A)/2.
      DA=(B-A)/2.
 
      DO 5000 I=1,NMAX
        N=N*2
        DX=PI/(2*N)
 
!	Apply N-point formula after transforming the range to (-1,1)
        R1=0.0
        DO 3000 J=1,N
          A1=(2*J-1)*DX
          R1=R1+FUN(A0+DA*COS(A1))
3000    CONTINUE
        R1=R1*DX*2.
 
        DIF=R1-RINT
        RINT=R1
        NPT=NPT+N
        IF(I.GT.3.AND.ABS(DIF).LT.MAX(AEPS,ABS(RINT)*REPS)) RETURN
5000  CONTINUE

!	Integral fails to converge 
      IER=30
      END

!     ---------------------------------------------

      FUNCTION F(X)
      IMPLICIT REAL*16(A-H,O-Z)
      COMMON B

!     THE INTEGRAND (WITHOUT THE WEIGHT FUNCTION) FOR SUBROUTINE GAUCBY

      F=1.Q0/SQRT((X+1)*(B+X))
      END
