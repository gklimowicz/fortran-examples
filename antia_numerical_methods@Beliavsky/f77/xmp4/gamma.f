!     TO EVALUATE GAMMA FUNCTION AND ITS NATURAL LOGARITHM

      PROGRAM GAM
!      IMPLICIT REAL*8(A-H,O-Z)
 
51    FORMAT('   X =',1PE14.6,5X,'GAMMA(X) =',E14.6,5X,
     1       'LOG(|GAMMA(X)|) =',E14.6)

 
100   PRINT *,'TYPE X    (QUITS WHEN X<-1000)'
      READ *,X
      IF(X.LT.-1000) STOP
      F1=GAMMA(X)
      F2=GAMMAL(X)
      WRITE(6,51) X,F1,F2
      GO TO 100
      END
 
!     ----------------------------------------------------
 
!	To calculate Gamma function for any real value of XG
!	Use GAMMAL for calculating the logarithm of Gamma function
!	which may be useful for large arguments or when argument is
!	close to a negative integer.
!
!	Required routines : None
 
      FUNCTION GAMMA(XG)
!      IMPLICIT REAL*8(A-H,O-Z)
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
 
!     ----------------------------------------------------
 
!	To calculate the Logarithm of Gamma function for a real argument
!	For negative values it give ln(abs(Gamma(x)))
!
!	Required routines : None
 
      FUNCTION GAMMAL(X)
!      IMPLICIT REAL*8(A-H,O-Z)
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
