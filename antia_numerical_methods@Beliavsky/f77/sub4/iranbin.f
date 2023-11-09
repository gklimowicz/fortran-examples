!       To generate random numbers with Binomial distribution 
!       For large values of mean it approximates the distribution by a
!       normal distribution with the same mean and standard deviation.
!
!       SEED : (input/output) real seed, it should be negative during
!              the first call and should not be modified between two
!              calls with the same N, P. If a new value of N or P is used
!              it should be reset to a negative value.
!       N   :  (input) Number of trials in the binomial distribution
!       P   :  (input) probability of the event
!       C   :  (output) Real array of length N used to store the
!              cumulative probability table for use in subsequent
!              calculations. This array should not be modified by the
!              user. This is used only if the mean N*P<RMAX and there
!              is no underflow while calculating (1-P)**N
!
!       Required routines : None
!
!       If P is close to 1 the approximation by Gaussian may not be very
!       good and in that case increasing RMAX may help although it will
!       take longer time to generate the numbers.

      FUNCTION IRANBIN(SEED,N,P,C)
!      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION C(*)
      REAL*8 AM,A,AC,AN
      PARAMETER(AM=2147483648D0,A=45875,AC=453816693D0,AN=2147483647D0)
      PARAMETER(PI=3.14159265358979324D0,RMAX=100)

      RMU=N*P
      SIG=SQRT(RMU*(1-P))
      P0=(1-P)**N

!       For large mu of if P0 gives underflow use normal distribution
      IF(RMU.GT.RMAX.OR.P0.EQ.0.0) THEN
!       Initialise the seed during first call
        IF(SEED.LT.0.0) SEED=MOD(ABS(SEED),AM)
        R1=MOD(A*SEED+AC,AM)
        IF(SEED.EQ.0.0) SEED=0.1
        RANGAU=SQRT(2.0*LOG(AN/SEED))*COS(2.0*PI*R1/AN)
        SEED=R1
        IRANBIN=RMU+SIG*RANGAU+0.5
        IF(IRANBIN.LT.0) IRANBIN=0
        IF(IRANBIN.GT.N) IRANBIN=N
      ELSE
!       Initialise the array during the first call
        IF(SEED.LT.0) THEN
          C(1)=P0
          DO I=2,N+1
            P0=P0*P*(N-I+2)/((1-P)*(I-1))
            C(I)=C(I-1)+P0
          ENDDO
          SEED=ABS(SEED)
        ENDIF

!       generate the random number
        SEED=MOD(SEED*A+AC,AM)
        R=SEED/AN
        DO I=1,N+1
          IF(R.LT.C(I)) EXIT
        ENDDO
        IF(I.GT.N+1) I=N+1
        IRANBIN=I-1
      ENDIF
      END
