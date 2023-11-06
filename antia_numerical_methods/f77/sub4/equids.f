!	Multiple integration over a hyper-rectangle in n-dimensions
!	using equidistributed sequences
!	Because of large roundoff error in summation it may not be
!	advisable to use this routine in single precision.
!
!	A : (input) Real array of length N containing the lower limit
!		along each dimension
!	B : (input) Real array of length N containing the upper limit
!		along each dimension
!	N : (input) The number of dimensions
!	NPT : (input) Maximum number of function evaluations to be used
!	F : (input) Name of the function routine to calculate the integrand
!		FUNCTION F(N,X) should calculate the integrand, where N is the
!		number of dimensions and X is a real array of length N containing
!		the coordinates of the point where integrand is to be calculated
!	S1 : (output) The calculated value of the integral
!	S2 : (output) Another approximation to the value of the integral
!		For smooth functions S2 is expected to be better approximation
!	REPS : (input) The required relative accuracy
!	AEPS : (input) The required absolute accuracy
!		The estimated error should be less than MAX(AEPS,REPS*ABS(S2))
!	DIF : (output) estimated (absolute) error achieved by the subroutine
!	NP : (output) Number of function evaluations used by subroutine
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=39 implies specified accuracy was not achieved in
!			which case DIF will contain the estimated accuracy
!		IER=312 implies N<1 or N>21 and no calculations are done
!
!	FUNCTION F(N,X) must be supplied by the user
!	
!	Required routines :  F
!
      SUBROUTINE EQUIDS(A,B,N,NPT,F,S1,S2,REPS,AEPS,DIF,NP,IER)
!      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NCHK=100,NMAX=21)
      DIMENSION A(N),B(N),H(NMAX),XA(NMAX),WT(NMAX),AT(NMAX)
!	The first 21 prime numbers
      DATA AT/2.,3.,5.,7.,11.,13.,17.,19.,23.,29.,31.,37.,41.,
     1        43.,47.,53.,59.,61.,67.,71.,73./

      IER=312
      S1=F(N,A)
      S2=S1
      NP=1
      IF(N.GT.NMAX.OR.N.LT.1) RETURN

      IER=0
      HH=1.0
      DO 1000 I=1,N
        H(I)=B(I)-A(I)
        HH=HH*H(I)
!	The irrational numbers for generating equidistributed sequences
        WT(I)=SQRT(AT(I))
1000  CONTINUE

      RI=0.0
      RI1=0.0
      DIF=0.0
      NPT1=NCHK
      DO 2000 I=1,NPT
!	Generate the abscissas using equidistributed sequences
        DO 1500 J=1,N
          A1=I*WT(J)
          A1=2.*ABS(A1-INT(A1+0.5))*H(J)
          XA(J)=A(J)+A1
1500    CONTINUE
!	Accumulate the sum
        S1=S1+2.*(F(N,XA)-RI)
        S2=S2+S1

        IF(MOD(I,NCHK).EQ.0) THEN
!	To control the roundoff error form partial sums
          SS1=RI+S1/(2*I+1)
          DIFF=S2/((I+1.)**2)
          S2=0.0
          S1=S1-(2*I+1)*DIFF
!	The new approximation to the average value of function
          RI=RI+DIFF

          IF(I.EQ.NPT1) THEN
!	Check for convergence
            DIF1=DIF
            DIF=ABS(RI-RI1)
            RI1=RI
            IF(DIF+DIF1.LT.MAX(AEPS/HH,ABS(RI)*REPS).AND.I.GT.5*NCHK)
     1         THEN
              S1=SS1*HH
              S2=RI*HH
              DIF=(DIF+DIF1)*HH
              NP=I+1
              RETURN
            ENDIF

            NPT1=NPT1*2
          ENDIF
        ENDIF
2000  CONTINUE

!	Integral fails to converge
      IER=39
      S1=(RI+S1/(2*NPT+1))*HH
      S2=(RI+S2/((NPT+1.)**2))*HH
      DIF=(DIF+DIF1)*HH
      NP=NPT+1
      END
