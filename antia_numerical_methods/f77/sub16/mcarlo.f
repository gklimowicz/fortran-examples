!	Multiple integration over a hyper-rectangle in n-dimensions
!	using Monte-Carlo technique
!
!	A : (input) Real array of length N containing the lower limit
!		along each dimension
!	B : (input) Real array of length N containing the upper limit
!		along each dimension
!	N : (input) The number of dimensions
!	NPT : (input) Maximum number of function evaluations to be used
!	F : (input) Name of the function routine to calculate the integrand
!		FUNCTION(N,X) should calculate the integrand, where N is the
!		number of dimensions and X is a real array of length N containing
!		the coordinates of the point where integrand is to be calculated
!	RI : (output) The calculated value of the integral
!	REPS : (input) The required relative accuracy
!	AEPS : (input) The required absolute accuracy
!		The estimated error should be less than MAX(AEPS,REPS*ABS(RI))
!	ERR : (output) estimated (absolute) error achieved by the subroutine
!		It is 2.576 times the estimated standard deviation
!		divided by SQRT(NP).
!	NP : (output) Number of function evaluations used by subroutine
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=39 implies specified accuracy was not achieved in
!			which case ERR will contain the estimated accuracy
!		IER=311 implies N<1 or N>50 and no calculations are done
!
!	FUNCTION F(N,X) should be supplied by the user
!	
!	Required routines : RANF, F
!
      SUBROUTINE MCARLO(A,B,N,NPT,F,RI,REPS,AEPS,ERR,NP,IER)
      IMPLICIT REAL*16(A-H,O-Z)
      PARAMETER(NMAX=50,NCHK=100)
!	To override the internal RANF if it exists
      EXTERNAL RANF
      DIMENSION A(N),B(N),H(NMAX),XA(NMAX)

      IER=311
      RI=0.0
      NP=0
      IF(N.GT.NMAX.OR.N.LT.1) RETURN

      IER=0
      HH=1.0
      DO 1000 I=1,N
        H(I)=B(I)-A(I)
1000  HH=HH*H(I)

      RI1=0.0
      VAR1=0.0
!	Seed for random number generator, should be changed if another routine is used
      ISEED=-12345
      NPT1=NCHK

      DO 2000 I=1,NPT
!	Generating the abscissas
        DO 1500 J=1,N
1500    XA(J)=A(J)+H(J)*RANF(ISEED)
        F1=F(N,XA)
        RI1=RI1+F1
        VAR1=VAR1+F1*F1

        IF(I.EQ.NPT1) THEN
!	Compute intermediate sums to check for convergence
          RI=RI1/I
          VAR=VAR1/I-RI*RI
          IF(VAR.LT.0.0) VAR=0.0
          ERR=2.576Q0*HH*SQRT(VAR/NPT1)
          RI=RI*HH
          NP=I
          IF(ERR.LT.MAX(AEPS/HH,ABS(RI)*REPS)) RETURN
          NPT1=2*NPT1
        ENDIF
2000  CONTINUE

!	Integral fails to converge
      RI=RI1/NPT
      VAR=VAR1/NPT-RI*RI
      ERR=2.576Q0*HH*SQRT(VAR/NPT)
      RI=RI*HH
      IF(ERR.GT.MAX(AEPS,ABS(RI)*REPS)) IER=39
      NP=NPT
      END
