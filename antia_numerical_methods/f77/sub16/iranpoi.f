!       To generate random numbers with Poisson distribution 
!       For large values of mean it approximates the distribution by a
!       normal distribution with the same mean and standard deviation.
!
!       SEED : (input/output) real seed, it should be negative during
!              the first call and should not be modified between two
!              calls with the same RMU. If a new value of RMU is used
!              it should be reset to a negative value.
!       RMU :  (input) Mean value of the Poisson distribution
!       P   :  (output) Real array of length NMAX used to store the
!              cumulative probability table for use in subsequent
!              calculations. This array should not be modified by the
!              user. This is used only if RMU is smaller than some
!              critical value.
!
!       Required routines : None
!

      FUNCTION IRANPOI(SEED,RMU,P)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION P(*)
      PARAMETER(AM=2147483648Q0,A=45875,AC=453816693Q0,AN=2147483647Q0)
      PARAMETER(PI=3.14159265358979323846264338327950288Q0,NMAX=800)

      SIG=SQRT(RMU)
      N=MAX(20Q0,RMU+7*SIG)

!       For large mu use normal distribution
      IF(N.GT.NMAX) THEN
!       Initialise the seed during first call
        IF(SEED.LT.0.0) SEED=MOD(ABS(SEED),AM)
        R1=MOD(A*SEED+AC,AM)
        IF(SEED.EQ.0.0) SEED=0.1Q0
        RANGAU=SQRT(2.Q0*LOG(AN/SEED))*COS(2.0*PI*R1/AN)
        SEED=R1
        IRANPOI=RMU+SIG*RANGAU+0.5
        IF(IRANPOI.LT.0) IRANPOI=0
      ELSE
!       Initialise the array during the first call
        IF(SEED.LT.0) THEN
          P0=EXP(-RMU)
          P(1)=P0
          DO I=2,N
            P0=P0*RMU/(I-1)
            P(I)=P(I-1)+P0
          ENDDO
          SEED=ABS(SEED)
        ENDIF

!       generate the random number
        SEED=MOD(SEED*A+AC,AM)
        R=SEED/AN
        DO I=1,N
          IF(R.LT.P(I)) EXIT
        ENDDO
        IF(I.GT.N) I=N
        IRANPOI=I-1
      ENDIF
      END

