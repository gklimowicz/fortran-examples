!       Calculate the correlation coefficient and the probability of
!       finding a value as high as that for uncorrelated dats sets

      PROGRAM XYCOR
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(1000),Y(1000)

51    FORMAT(/,'  N =',I5,'  Correlation Coefficient =',F8.4,/,
     1  'Probability of finding this value or higher for uncorrelated',
     2  ' sequence =',1PD12.3)

      SEED=9
100   PRINT *,'TYPE N = NO. OF DATA POINTS (QUITS WHEN N.LE.0)'
      READ *,N
      IF(N.LE.0.OR.N.GT.1000) STOP
!       Generate the data
      DO I=1,N
        R=RAN1(SEED)
        X(I)=I*1D-2+1.D-3*(R-0.5)
        R=RAN1(SEED)
        Y(I)=I*0.001D-2+1.D-2*(R-0.5)
      ENDDO

!       Calculate the correlation coefficient between X and Y
      XMEAN=0.0
      YMEAN=0.0
      DO I=1,N
        XMEAN=XMEAN+X(I)
        YMEAN=YMEAN+Y(I)
      ENDDO
      XMEAN=XMEAN/N
      YMEAN=YMEAN/N
      SXX=0.0
      SYY=0.0
      SXY=0.0
      DO I=1,N
        SXX=SXX+(X(I)-XMEAN)**2
        SYY=SYY+(Y(I)-YMEAN)**2
        SXY=SXY+(X(I)-XMEAN)*(Y(I)-YMEAN)
      ENDDO
      RCOR=SXY/SQRT(SXX*SYY)
      PROB=PCOR(N-2,ABS(RCOR))
      WRITE(6,51) N,RCOR,PROB
      GO TO 100
      END

!       -----------------------------------------------------

!	To calculate the Error function for a real argument
!
!	Required routines : None
 
      FUNCTION ERF(X)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(7),B(8),A1(7),B1(8)

!	The coefficients of rational function approximations
      DATA A/ 4.891837297874833514D-01,  1.110175596090841322D-01,
     1        1.526977188817169289D-02,  1.388143322498740953D-03,
     1        8.446154421878829637D-05,  3.239915935631695227D-06,
     1        6.200069065781009292D-08/
      DATA B/ 1.128379167095512575D+00,  1.758583405424390318D-01,
     1        5.411290829613803886D-02,  3.805760982281134630D-03,
     1        4.314620532020313106D-04,  1.423201530735308681D-05,
     1        6.188698151904667110D-07,  2.293915472550064153D-09/
      DATA A1/ 3.040844316651787557D+01,  3.358927547920980388D+02,
     1         1.703170048554673596D+03,  4.133195606956137156D+03,
     1         4.555611776312694034D+03,  1.935778559646575488D+03,
     1         2.051089518116697193D+02/
      DATA B1/ 5.641895835477562749D-01,  1.687403209467950089D+01,
     1         1.813522721872712655D+02,  8.779664433851135986D+02,
     1         1.965115658619443782D+03,  1.865558781108286245D+03,
     1         5.828865035607128761D+02,  2.558559157228883880D+01/
 
      IF(ABS(X).LT.2.0) THEN
!	Use approximation for Erf(x)/x
        Y=X*X
        FN=((((((B(8)*Y+B(7))*Y+B(6))*Y+B(5))*Y+B(4))*Y+B(3))*Y+B(2))*Y
     1      +B(1)
        FD=((((((A(7)*Y+A(6))*Y+A(5))*Y+A(4))*Y+A(3))*Y+A(2))*Y
     1      +A(1))*Y+1
        ERF=X*FN/FD
      ELSE
 
!	Use approximation for x*Exp(x^2)*Erfc(x)
        X1=1./X**2
        FN=((((((B1(8)*X1+B1(7))*X1+B1(6))*X1+B1(5))*X1+B1(4))*X1+
     1      B1(3))*X1+B1(2))*X1+B1(1)
        FD=((((((A1(7)*X1+A1(6))*X1+A1(5))*X1+A1(4))*X1+A1(3))*X1+
     1      A1(2))*X1+A1(1))*X1+1
        ERF=1.-EXP(-X*X)*FN/(FD*ABS(X))
        IF(X.LT.0.0) ERF=-ERF
      ENDIF
 
      END

!       -----------------------------------------------------

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

!       -----------------------------------------------------

!	To calculate the probability that two uncorrelated sequences
!	of length (n+2) will give a correlation coefficient exceeding |x|
!	For large even values of n the series can give large roundoff
!	in which case the distribution is approximated by Normal distribution
!
!       N  : (input) Length of data sets should be N+2
!       XX : (input) The function will calculate the probability for
!               correlation exceeding |x| for uncorrelated sequences
!
!	Required routines : GAMMAL, ERF

      FUNCTION PCOR(N,XX)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(PI=3.14159265358979323846D0,PIS=1.7724538509055144D0)

      X=ABS(XX)
      PCOR=0.0
      IF(X.GT.1.0) RETURN
      IF(2*(N/2).EQ.N) THEN
!	if n is even
        L=(N-2)/2
        SS=X
        P=X
        PMAX=P
        DO K=1,l
          P=-(L-K+1)*P*X*X/K
          SS=SS+P/(2*K+1.D0)
          IF(ABS(P).GT.PMAX) PMAX=ABS(P)
        ENDDO
        PL=GAMMAL((N+1.D0)/2)-GAMMAL(N/2.D0)
        PCOR=2.*EXP(PL)*SS/PIS
        IF(PMAX.GT.1.D5.OR.PCOR.GT.1.0) PCOR=ERF(X*SQRT(N/2.D0))
      ELSE
!	if n is odd
        L=(N-3)/2
        SS=SQRT(1.-X*X)
        P=SS
        IF(N.EQ.1) SS=0.0
        DO K=1,L
          P=P*(2*K)*(1-X*X)/(2*K+1.D0)
          SS=SS+P
        ENDDO
        PCOR=(ASIN(X)+X*SS)*2/PI
      ENDIF
      PCOR=1-PCOR
      END

!       -----------------------------------------------------

!	To generate uniformly distributed random numbers in interval (0,1)
!
!	SEED : (input/output) is a real value used as the seed
!		It should be positive during initial call and
!		should not be modified between different calls.
!
!	Required routines : None

      FUNCTION RAN1(SEED)
      IMPLICIT REAL*8(A-H,O-Z)
!	Retain the following declaration even for REAL*4 version
!	otherwise AM and AC will be rounded
      REAL*8 AM,A,AC,AN
      PARAMETER(AM=2147483648D0,A=45875,AC=453816693D0,AN=2147483647D0)

      SEED=MOD(SEED*A+AC,AM)
      RAN1=SEED/AN
      END
