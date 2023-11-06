!	To calculate Gamma function for any real value of XG
!	Use GAMMAL for calculating the logarithm of Gamma function
!	which may be useful for large arguments or when argument is
!	close to a negative integer.
!
!	Required routines : None
 
      FUNCTION GAMMA(XG)
      IMPLICIT REAL*16(A-H,O-Z)
!	PIS is SQRT(2*PI)
      PARAMETER(PI=3.14159265358979323846Q0,PIS=2.5066282746310005024Q0)
      DIMENSION A(2),B(3),A1(6),B1(7)

!	The coefficients for rational function approximations
      DATA A/1.767971449569122937Q+00,  2.909421117928672645Q-01/
      DATA B/8.333333333333231537Q-02,  1.445531763554246280Q-01,
     1       2.012779361583001035Q-02/
      DATA A1/3.905731686764559737Q+03,  2.204952264401381785Q+03,
     1       -1.932467485468849660Q+03,  4.643360871045442213Q+02,
     1       -4.818088806916028754Q+01,  1.896853765546068169Q+00/
      DATA B1/3.918055655523400310Q+03, -1.088116266563809683Q+02,
     1        8.203258626193993149Q+02, -9.289402000761705906Q+01,
     1        6.521113026294866877Q+01, -6.090618615608719044Q+00,
     1        1.475909104740280784Q+00/
 
      X=ABS(XG)
      IF(X.GT.1000.0) THEN
!	Use asymptotic approximation (Stirling formula)
        GX=(1+1.Q0/(12*X)+1./(288*X*X)
     1         -139/(51840*X**3)-571./(2488320Q0*X**4))
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
