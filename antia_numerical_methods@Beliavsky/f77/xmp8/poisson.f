!       To simulate a Poisson distribution
      PROGRAM RANPOI
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION C(2000),IP(2000)

51    FORMAT(I5,1P5D13.5)
52    FORMAT('    I    PROB(I) ')
53    FORMAT(/,'  Mean =',1PD13.5,'   Chi-square = ',D13.5,
     1  '  Degree of freedom =',I4)
54    FORMAT(' Probability of finding a value of CHI or larger =',
     1               1PD13.5)

      RN=1.D4
      NUM=RN
100   PRINT *,' TYPE RMU = Mean of the Poisson distribution (>0,<700)'
      READ *,RMU
!       Larger value of RMU may give underflow
      IF(RMU.GT.700.OR.RMU.LE.0) STOP
      SEED=-1
      DO I=1,2000
        IP(I)=0.0
      ENDDO
!       Generate NUM random numbers with binomial distribution and
!       place them in the appropriate bin
      DO I=1,NUM
        II=IRANPOI(SEED,RMU,C)
        IP(II+1)=IP(II+1)+1
      ENDDO

!       Apply the Chi-square test to check the distribution
!       p0 is the expected probability
      P0=EXP(-RMU)
      CHI=0.0
      WRITE(6,52)
      NU=0
      P1=0
!       Use only the bins that have non-zero counts
      DO I=1,2000
        IF(IP(I).GT.0) THEN
          WRITE(6,51)I-1,IP(I)/RN,P0
          CHI=CHI+(P0*RN-IP(I))**2/(P0*RN)
          NU=NU+1
        ELSE
          P1=P1+P0
        ENDIF
        P0=P0*RMU/I
      ENDDO
!       Add the contribution from bins that had no event
!       NU should actually be higher as all bins with no event may not
!       be consecutive. Thus the Chi-square test may fail.
      CHI=CHI+P1*RN
      WRITE(6,53) RMU,CHI,NU
      PROB=1-GAMMAP(NU/2D0,CHI/2)
      WRITE(6,54) PROB
      GO TO 100
      END

!       ------------------------------------------------------------

!	To calculate Gamma function for any real value of XG
!	Use GAMMAL for calculating the logarithm of Gamma function
!	which may be useful for large arguments or when argument is
!	close to a negative integer.
!
!	Required routines : None
 
      FUNCTION GAMMA(XG)
      IMPLICIT REAL*8(A-H,O-Z)
!	PIS is SQRT(2*PI)
      PARAMETER(PI=3.14159265358979323846D0,PIS=2.5066282746310005024D0)
      DIMENSION A(2),B(3),A1(6),B1(7)

!	The coefficients for rational function approximations
      DATA A/1.767971449569122937D+00,  2.909421117928672645D-01/
      DATA B/8.333333333333231537D-02,  1.445531763554246280D-01,
     1       2.012779361583001035D-02/
      DATA A1/3.905731686764559737D+03,  2.204952264401381785D+03,
     1       -1.932467485468849660D+03,  4.643360871045442213D+02,
     1       -4.818088806916028754D+01,  1.896853765546068169D+00/
      DATA B1/3.918055655523400310D+03, -1.088116266563809683D+02,
     1        8.203258626193993149D+02, -9.289402000761705906D+01,
     1        6.521113026294866877D+01, -6.090618615608719044D+00,
     1        1.475909104740280784D+00/
 
      X=ABS(XG)
      IF(X.GT.1000.0) THEN
!	Use asymptotic approximation (Stirling formula)
        GX=(1+1.D0/(12*X)+1./(288*X*X)
     1         -139/(51840*X**3)-571./(2488320D0*X**4))
        GAMMA=X**(X-0.5)*EXP(-X)*PIS*GX

      ELSE IF(X.GT.8.0) THEN
!	Use rational function approximation for Log(Gamma) 
        Y=1./X**2
        RMK=((B(3)*Y+B(2))*Y+B(1))/((A(2)*Y+A(1))*Y+1)
        GAMMA=X**(X-0.5)*EXP(-X)*PIS*EXP(RMK/X)

      ELSE IF(X.GE.2.0) THEN
!	Use rational function approximation for (Gamma) over [2,3]
!	after translating the range if necessary
        F1=1.0
        X1=X
2500    IF(X1.LE.3) GO TO 3000
        F1=F1*(X1-1)
        X1=X1-1
        GO TO 2500
3000    IF(X1.EQ.3) THEN
          GAMMA=F1*2
        ELSE IF(X1.EQ.2) THEN
          GAMMA=F1
        ENDIF

        FN=(((((B1(7)*X1+B1(6))*X1+B1(5))*X1+B1(4))*X1+B1(3))*X1+
     1          B1(2))*X1+B1(1)
        FD=(((((A1(6)*X1+A1(5))*X1+A1(4))*X1+A1(3))*X1+A1(2))*X1+
     1          A1(1))*X1+1
        GAMMA=F1*FN/FD
      ELSE IF(X.GT.0.0) THEN
!	Use rational function approximation for (Gamma) over [2,3]
!	after translating the range if necessary
        F1=1./X
        X1=X+1
        IF(X.LT.1) THEN
          F1=F1/X1
          X1=X1+1
        ENDIF
        IF(X1.EQ.2) GAMMA=F1

        FN=(((((B1(7)*X1+B1(6))*X1+B1(5))*X1+B1(4))*X1+B1(3))*X1+
     1          B1(2))*X1+B1(1)
        FD=(((((A1(6)*X1+A1(5))*X1+A1(4))*X1+A1(3))*X1+A1(2))*X1+
     1          A1(1))*X1+1
        GAMMA=F1*FN/FD
 
      ENDIF

      IF(XG.GT.0.0) RETURN
      IX=X
      IF(X.GT.IX) THEN
        GAMMA=PI/(XG*SIN(PI*X)*GAMMA)
      ELSE
        GAMMA=(-1)**IX/0.0
      ENDIF
 
      END

!       ------------------------------------------------------------

!	To calculate the Logarithm of Gamma function for a real argument
!	For negative values it give ln(abs(Gamma(x)))
!
!	Required routines : None
 
      FUNCTION GAMMAL(X)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(4),B(4),A1(6),B1(6)
!	PI2L=LOG(2*PI)/2
      PARAMETER (PI2L=0.918938533204672741D0,PI=3.141592653589793238D0)

!	The coefficients of rational function approximations
      DATA A/ 1.013142782275024216D-2,  7.645657825398191944D-1,
     1        3.381172379819461227D-4,  1.595363637547538209D-2/
      DATA B/0.0, 8.333333333338911768D-2, 8.442856404442060242D-4,
     1            6.093603832366013704D-2/

      DATA A1/ 4.681163846241230144D0,  3.208225429683256526D0,
     1         5.145525793448859216D-1, 1.581117883959157936D-2,
     2        -6.398416804905407512D-5, 5.264566254181773919D-7/
      DATA B1/ 2.938038561191284576D0,  1.489364948862436743D0,
     1        -5.466291543917642961D0,  1.972497734170110410D-1,
     2         7.830146473241555157D-1, 5.756753067834747499D-2/
 
      T=ABS(X)
      IT=T
      IF(T.GE.10.0) THEN
!	Use asymptotic approximation for T>10
        Y=1/T
        FN=((B(4)*Y+B(3))*Y+B(2))*Y+B(1)
        FD=(((A(4)*Y+A(3))*Y+A(2))*Y+A(1))*Y+1
        GAMMAL=FN/FD-T+PI2L+(T-0.5)*LOG(T)
      ELSE
 
!	Use approximation for [4,5]
        FAC=0.0
        T1=T
        IF(T.LT.4) THEN
          DO I=IT,3
            FAC=FAC-LOG(T1)
            T1=T1+1
          ENDDO
        ELSE IF(T.GT.5) THEN
          DO I=IT,5,-1
            T1=T1-1
            FAC=FAC+LOG(T1)
          ENDDO
        ENDIF
        FN=((((B1(6)*T1+B1(5))*T1+B1(4))*T1+B1(3))*T1+B1(2))*T1+B1(1)
        FD=(((((A1(6)*T1+A1(5))*T1+A1(4))*T1+A1(3))*T1+A1(2))*T1+
     1          A1(1))*T1+1
        GAMMAL=FN/FD+FAC
      ENDIF

      IF(X.LE.0) THEN
        IF(T.GT.IT) THEN
          GAMMAL=LOG(PI)-LOG(ABS(SIN(PI*T)*T))-GAMMAL
        ELSE
          GAMMAL=(-1)**IT/0.0
        ENDIF
      ENDIF
 
      END

!       ------------------------------------------------------------

!	To calculate the incomplete Gamma function P(a,x) for
!	positive real arguments
!	It returns a value of -1 for a.le.0 or x<0
!
!	A : (input) Argument for the complete Gamma function
!	X : (input) Upper limit of integration defining the incomplete
!                       Gamma function
!
!	Required routines : GAMMA, GAMMAL

      FUNCTION GAMMAP(A,X)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION C(500),B(500)

      IF(A.LE.0.OR.X.LT.0) THEN
        GAMMAP=-1
        RETURN
      ENDIF

      IF(X.LT.3) THEN
!	Use power series
        S1=1
        T=A
        DO I=1,500
          T=-T*X/I
          S1=S1+T/(A+I)
          IF(ABS(T/(A+I)).LT.1.D-14) EXIT
        ENDDO
        IF(A.LT.140) THEN
          GAMMAP=X**A*S1/GAMMA(A+1)
        ELSE
          GAMMAP=0.0
        ENDIF
      ELSE IF(A.LT.1.2D0*X) THEN
!	Use continued fraction
        C(1)=1
        B(1)=X
        C(2)=1
        B(2)=X+1-A
        C1=C(2)/B(2)
        DO I=2,498,2
          C(I+1)=X*C(I)+(I/2)*C(I-1)
          B(I+1)=X*B(I)+(I/2)*B(I-1)
          C(I+2)=C(I+1)+(I/2+1-A)*C(I)
          B(I+2)=B(I+1)+(I/2+1-A)*B(I)
          IF(ABS(B(I+2)).GT.1.D200) THEN
            C(I+1)=C(I+1)/1.D200
            B(I+1)=B(I+1)/1.D200
            C(I+2)=C(I+2)/1.D200
            B(I+2)=B(I+2)/1.D200
          ENDIF
          IF(ABS(B(I+2)).LT.1.D-200) THEN
            C(I+1)=C(I+1)*1.D200
            B(I+1)=B(I+1)*1.D200
            C(I+2)=C(I+2)*1.D200
            B(I+2)=B(I+2)*1.D200
          ENDIF
          C2=C(I+2)/B(I+2)
          IF(ABS(C2-C1).LT.1.D-12) EXIT
          C1=C2
        ENDDO
        G1=-X+A*LOG(X)+LOG(C2)-GAMMAL(A)
        GAMMAP=1-EXP(G1)
      ELSE
!	Use the power series for a>x
        S2=1
        T=1
        DO I=1,500
          T=T*X/(A+I)
          S2=S2+T
          IF(ABS(T).LT.1.D-14) EXIT
        ENDDO
        G1=-X+A*LOG(X)+LOG(S2)-GAMMAL(A+1)
        GAMMAP=EXP(G1)
      ENDIF
      END

!       ------------------------------------------------------------

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
!       Increasing NMAX may give underflow while calculating EXP(-RMU)

      FUNCTION IRANPOI(SEED,RMU,P)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION P(*)
      PARAMETER(AM=2147483648D0,A=45875,AC=453816693D0,AN=2147483647D0)
      PARAMETER(PI=3.14159265358979324D0,NMAX=800)

      SIG=SQRT(RMU)
      N=MAX(20D0,RMU+7*SIG)

!       For large mu use normal distribution
      IF(N.GT.NMAX) THEN
!       Initialise the seed during first call
        IF(SEED.LT.0.0) SEED=MOD(ABS(SEED),AM)
        R1=MOD(A*SEED+AC,AM)
        IF(SEED.EQ.0.0) SEED=0.1D0
        RANGAU=SQRT(2.D0*LOG(AN/SEED))*COS(2.0*PI*R1/AN)
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

