!	To calculate the Bessel function of first and second kind of order one
!	for real positive argument
!	For XB.LE.0 the function of second kind is not defined and subroutine
!	will return zero value  without any error message or flag.
!
!	XB : (input) Argument at which the functions need to be evaluated
!	BJ1 : (output) Calculated value of Bessel function of first kind
!	BY1 : (output) Calculated value of Bessel function of second kind
!
!	Required routines : None
 
      SUBROUTINE BJY1(XB,BJ1,BY1)
      IMPLICIT REAL*16(A-H,O-Z)
!	EUGAM is the Euler's constant
      PARAMETER(EUGAM=0.5772156649015328606Q0,PI=3.14159265358979324Q0)
      DIMENSION A(7),B(8),A1(5),B1(5),A2(5),B2(5),A0(7),B0(8)
 
!	Coefficients of rational function approximations
      DATA A/ 1.156878510067059849Q-02,  6.749406787503859073Q-05,
     1        2.614560837317535451Q-07,  7.408815126464007290Q-10,
     1        1.577637796406197189Q-12,  2.432945305413635549Q-15,
     1        2.257446839754248784Q-18/
      DATA B/ 5.000000000000000074Q-01, -5.671560744966475746Q-02,
     1        1.914864631812345532Q-03, -2.821407888958585592Q-05,
     1        2.103168789640803591Q-07, -8.322474383730280556Q-10,
     1        1.678871778708754849Q-12, -1.372424374400306547Q-15/

      DATA A0/ 1.186963690270342970Q-02,  7.123839029323002736Q-05,
     1         2.848196759538669928Q-07,  8.365303089083305885Q-10,
     1         1.857096246589275584Q-12,  3.012506935004947491Q-15,
     1         2.996951174746838817Q-18/
      DATA B0/ 2.500000000000000299Q-01, -7.515759077432437273Q-02,
     1         3.430771992327672576Q-03, -6.022315614557372919Q-05,
     1         5.067136874996839630Q-07, -2.197514674456554803Q-09,
     1         4.768619679411702014Q-12, -4.139491442515065355Q-15/

      DATA A1/ 8.659888261699365129Q+01,
     1         1.932665751369749084Q+03,  1.172714583536277145Q+04,
     1         1.256737699073784218Q+04, -6.147124347503755010Q+02/
      DATA B1/ 1.000000000000000011Q+00,
     1         8.671607011699346720Q+01,  1.942669862370300601Q+03,
     1         1.194181952104744095Q+04,  1.371467864585746530Q+04/

      DATA A2/  1.021472573795463627Q+02,  2.807865400111916226Q+03,
     1          2.280402060738415865Q+04,  4.121116954504273053Q+04,
     1          3.501974669280301705Q+03/
      DATA B2/  3.749999999999999461Q-01,  3.820268245483084309Q+01,
     1          1.042753017477090289Q+03,  8.289951986135169400Q+03,
     1          1.371889615877945967Q+04/
 
      X=ABS(XB)
      IF(X.LT.8.0) THEN
!	Use rational function approximations
        Y=X*X
        FN=((((((B(8)*Y+B(7))*Y+B(6))*Y+B(5))*Y+B(4))*Y+B(3))*Y+B(2))*Y
     1          +B(1)
        FD=((((((A(7)*Y+A(6))*Y+A(5))*Y+A(4))*Y+A(3))*Y+A(2))*Y+A(1))*Y
     1          +1
        BJ1=X*FN/FD
 
        FN=((((((B0(8)*Y+B0(7))*Y+B0(6))*Y+B0(5))*Y+B0(4))*Y+B0(3))*Y
     1          +B0(2))*Y +B0(1)
        FD=((((((A0(7)*Y+A0(6))*Y+A0(5))*Y+A0(4))*Y+A0(3))*Y+A0(2))*Y
     1          +A0(1))*Y +1
        IF(X.GT.0.0) BY1=2.*(BJ1*(LOG(X/2)+EUGAM)-1./X-X*FN/FD)/PI
 
      ELSE
!	Use rational function approximations for P_1 and Q_1
        X1=1./X**2
        FN=(((B1(5)*X1+B1(4))*X1+B1(3))*X1+B1(2))*X1+B1(1)
        FD=((((A1(5)*X1+A1(4))*X1+A1(3))*X1+A1(2))*X1+A1(1))*X1+1
        P1=FN/FD
 
        FN=(((B2(5)*X1+B2(4))*X1+B2(3))*X1+B2(2))*X1+B2(1)
        FD=((((A2(5)*X1+A2(4))*X1+A2(3))*X1+A2(2))*X1+A2(1))*X1+1
        Q1=FN/(FD*ABS(X))

        X1=ABS(X)
        X2=X1-0.75Q0*PI
        BY1=SQRT(2./(PI*X1))*(P1*SIN(X2)+Q1*COS(X2))
        BJ1=SQRT(2./(PI*X1))*(P1*COS(X2)-Q1*SIN(X2))
      ENDIF

      IF(XB.LE.0.0) BY1=0.0
 
      END