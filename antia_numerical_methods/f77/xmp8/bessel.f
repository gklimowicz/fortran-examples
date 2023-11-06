!     TO EVALUATE BESSEL FUNCTION OF VARIOUS TYPES OF INTEGRAL ORDER
!     It evaluates the Bessel functions of first and second kind,
!     modified Bessel functions of first and second kind and the
!     spherical Bessel function.
!     All these may not be defined for all values of X and N.
 
      PROGRAM BESSEL
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL BJ0,BJ1,BJN,BI0,BI1,BIN,BK0,BK1,BKN,BY0,BY1,BYN
      DIMENSION B(5000)
 
51    FORMAT('   X =',1PD14.6,5X,'N =',I4,5X,'J0(X) =',D14.6,5X,
     1       'J1(X) =',D14.6)
52    FORMAT('   BJN:   J0(X) =',1PD14.6,5X,'J1(X) =',D14.6,5X,
     1       'JN(X) =',D14.6)
53    FORMAT('   BJY0:   J0(X) =',1PD14.6,5X,'Y0(X) =',D14.6)
54    FORMAT('   BJY1:   J1(X) =',1PD14.6,5X,'Y1(X) =',D14.6)
55    FORMAT(5X,'Y0(X) =',1PD14.6,5X,'Y1(X) =',D14.6)
56    FORMAT('   BYN:   Y0(X) =',1PD14.6,5X,'Y1(X) =',D14.6,5X,
     1       'YN(X) =',D14.6)
57    FORMAT('   SPHBJN:   j0(X) =',1PD14.6,5X,'j1(X) =',D14.6,5X,
     1       'jN(X) =',D14.6)
58    FORMAT(5X,'I0(X) =',1PD14.6,5X,'I1(X) =',D14.6)
59    FORMAT('   BIN:   I0(X) =',1PD14.6,5X,'I1(X) =',D14.6,5X,
     1       'IN(X) =',D14.6)
60    FORMAT(5X,'K0(X) =',1PD14.6,5X,'K1(X) =',D14.6)
61    FORMAT('   BKN:   K0(X) =',1PD14.6,5X,'K1(X) =',D14.6,5X,
     1       'KN(X) =',D14.6)
 
!     For each type the Bessel functions of order 0,1 and N are printed
!     The functions of order 0, 1 are evaluated separately also
 
100   PRINT *,'TYPE X, N    (QUITS WHEN X<-1000)'
      READ *,X,N
      IF(X.LT.-1000) STOP
      NA=ABS(N)+1
 
!     Bessel function of first kind
      FJ0=BJ0(X)
      FJ1=BJ1(X)
      WRITE(6,51) X,N,FJ0,FJ1
      CALL BJN(N,X,B)
      WRITE(6,52) B(1),B(2),B(NA)
 
!     Bessel function of first and second kind
      CALL BJY0(X,GJ0,GY0)
      WRITE(6,53) GJ0,GY0
      CALL BJY1(X,GJ1,GY1)
      WRITE(6,54) GJ1,GY1
 
!     Bessel function of second kind
      FY0=BY0(X)
      FY1=BY1(X)
      WRITE(6,55) FY0,FY1
      CALL BYN(N,X,B)
      WRITE(6,56) B(1),B(2),B(NA)
 
!     Spherical Bessel function
      CALL SPHBJN(N,X,B)
      WRITE(6,57) B(1),B(2),B(NA)
 
!     Modified Bessel function of first kind
      FI0=BI0(X)
      FI1=BI1(X)
      WRITE(6,58) FI0,FI1
      CALL BIN(N,X,B)
      WRITE(6,59) B(1),B(2),B(NA)
 
!     Modified Bessel function of second kind
      FK0=BK0(X)
      FK1=BK1(X)
      WRITE(6,60) FK0,FK1
      CALL BKN(N,X,B)
      WRITE(6,61) B(1),B(2),B(NA)
 
      GO TO 100
      END
 
!     ----------------------------------------------------
 
!	To calculate the Bessel function of order zero for real argument
!
!	Required routines : None
 
      FUNCTION BJ0(X)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(PI=3.14159265358979324D0)
      DIMENSION A(8),B(8),A1(5),B1(5),A2(5),B2(5)
 
!	Coefficients of rational function approximations
      DATA A/ 1.293686560051304152D-02,  8.573459862295151747D-05,
     1        3.854769244046149702D-07,  1.308534328117880493D-09,
     1        3.512360907188715842D-12,  7.512575042421009221D-15,
     1        1.229302278444845702D-17,  1.311883486088925264D-20/
      DATA B/ 9.999999999999999878D-01, -2.370631343994868513D-01,
     1        1.247651819849453565D-02, -2.529374255010058573D-04,
     1        2.411267406461247155D-06, -1.159484705672466498D-08,
     1        2.730546745501229851D-11, -2.517936655103065990D-14/

      DATA A1/  8.911849018950665793D+01,
     1          2.078818787053760203D+03,  1.366258799766718466D+04,
     1          1.800383785973922830D+04,  4.923440494847201509D+02/
      DATA B1/  9.999999999999999908D-01,
     1          8.904817768950681616D+01,  2.072664795311476688D+03,
     1          1.352584337655551999D+04,  1.723138433448795889D+04/

      DATA A2/  1.046366195300779895D+02,  2.980727727381642723D+03,
     1          2.570418404044668245D+04,  5.291014161889741749D+04,
     1          9.497228391199055149D+03/
      DATA B2/ -1.249999999999999534D-01, -1.300633525376058087D+01,
     1         -3.651542590150084093D+02, -3.016744074771875522D+03,
     1         -5.251679479249748063D+03/
 
      IF(ABS(X).LT.8.0) THEN
!	Use rational function approximation
        Y=X*X
        FN=((((((B(8)*Y+B(7))*Y+B(6))*Y+B(5))*Y+B(4))*Y+B(3))*Y+B(2))*Y
     1      +B(1)
        FD=(((((((A(8)*Y+A(7))*Y+A(6))*Y+A(5))*Y+A(4))*Y+A(3))*Y+A(2))*Y
     1      +A(1))*Y+1
        BJ0=FN/FD
 
      ELSE
!	Use rational function approximations for P_0 and Q_0
        X1=1./X**2
        FN=(((B1(5)*X1+B1(4))*X1+B1(3))*X1+B1(2))*X1+B1(1)
        FD=((((A1(5)*X1+A1(4))*X1+A1(3))*X1+A1(2))*X1+A1(1))*X1+1
        P0=FN/FD
 
        FN=(((B2(5)*X1+B2(4))*X1+B2(3))*X1+B2(2))*X1+B2(1)
        FD=((((A2(5)*X1+A2(4))*X1+A2(3))*X1+A2(2))*X1+A2(1))*X1+1
        Q0=FN/(FD*ABS(X))

        X1=ABS(X)
        BJ0=SQRT(2./(PI*X1))*(P0*COS(X1-PI/4)-Q0*SIN(X1-PI/4))
      ENDIF
 
      END
 
!     ----------------------------------------------------
 
!	To calculate the Bessel function of order one for real argument
!
!	Required routines : None
 
      FUNCTION BJ1(X)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(PI=3.14159265358979324D0)
      DIMENSION A(7),B(8),A1(5),B1(5),A2(5),B2(5)
 
!	Coefficients of rational function approximations
      DATA A/ 1.156878510067059849D-02,  6.749406787503859073D-05,
     1        2.614560837317535451D-07,  7.408815126464007290D-10,
     1        1.577637796406197189D-12,  2.432945305413635549D-15,
     1        2.257446839754248784D-18/
      DATA B/ 5.000000000000000074D-01, -5.671560744966475746D-02,
     1        1.914864631812345532D-03, -2.821407888958585592D-05,
     1        2.103168789640803591D-07, -8.322474383730280556D-10,
     1        1.678871778708754849D-12, -1.372424374400306547D-15/
      DATA A1/  8.659888261699365129D+01,
     1          1.932665751369749084D+03,  1.172714583536277145D+04,
     1          1.256737699073784218D+04, -6.147124347503755010D+02/
      DATA B1/  1.000000000000000011D+00,
     1          8.671607011699346720D+01,  1.942669862370300601D+03,
     1          1.194181952104744095D+04,  1.371467864585746530D+04/
      DATA A2/  1.021472573795463627D+02,  2.807865400111916226D+03,
     1          2.280402060738415865D+04,  4.121116954504273053D+04,
     1          3.501974669280301705D+03/
      DATA B2/  3.749999999999999461D-01,  3.820268245483084309D+01,
     1          1.042753017477090289D+03,  8.289951986135169400D+03,
     1          1.371889615877945967D+04/
 
      IF(ABS(X).LT.8.0) THEN
!	Use rational function approximation
        Y=X*X
        FN=((((((B(8)*Y+B(7))*Y+B(6))*Y+B(5))*Y+B(4))*Y+B(3))*Y+B(2))*Y
     1      +B(1)
        FD=((((((A(7)*Y+A(6))*Y+A(5))*Y+A(4))*Y+A(3))*Y+A(2))*Y+A(1))*Y
     1      +1
        BJ1=X*FN/FD
 
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
        X2=X1-0.75D0*PI
        BJ1=SQRT(2./(PI*X1))*(P1*COS(X2)-Q1*SIN(X2))
        IF(X.LT.0.0) BJ1=-BJ1
      ENDIF
 
      END
 
!     ----------------------------------------------------
 
!	To calculate the Bessel function of integral order for real argument
!
!	N : (input) Order of Bessel function required, N may be negative
!		or positive. For N=0,1 use BJ0 and BJ1 respectively.
!	XB : (input) Argument at which the value is required
!	BJ : (output) Real array of length at least 
!		ABS(N)+16+MAX(25,5*SQRT(N))
!		which will contain the value of Bessel function of order
!		0,1,...,ABS(N). BJ(I+1) will contain Bessel function of order I
!		or -I (if N<0)
!		Remaining elements of array are used as scratch space
!
!	Required routines : BJ0, BJ1
!
      SUBROUTINE BJN(N,XB,BJ)
      IMPLICIT REAL*8(A-H,O-Z)
!	REPS should be less than the machine accuracy
      PARAMETER(REPS=1.D-17)
      DIMENSION BJ(*)
 
51    FORMAT('  BJN FAILED AT N =',I6,'  X =',1PD13.5,'  S =',D13.5)
 
      X=ABS(XB)
      NA=ABS(N)
      IF(XB.EQ.0) THEN
        BJ(1)=1.0
        DO 1800 I=2,NA+1
          BJ(I)=0.0
1800    CONTINUE
        RETURN

      ELSE IF(NA.LT.X.OR.NA.LT.2) THEN
!	Use the recurrence relation in the forward direction 
        BJ(1)=BJ0(X)
        BJ(2)=BJ1(X)
        DO 2000 I=2,NA
          BJ(I+1)=2*(I-1)*BJ(I)/X-BJ(I-1)
2000    CONTINUE
      ELSE IF(X.LE.4.0) THEN

!	Use series expansion to calculate  BJ(NA), BJ(NA-1)
        XA=X*X/4.0
        T0=X/2.0
        DO 2400 I=2,NA
          T0=T0*X/(2.0*I)
          IF(I.GE.NA-1) THEN
            T=T0
            S=T0
            DO 2200 J=1,50
              T=-T*XA/(J*(J+I))
              S=S+T
              IF(ABS(T).LT.ABS(S*REPS)) GO TO 2300
2200        CONTINUE
2300        BJ(I+1)=S
          ENDIF
2400    CONTINUE
 
        DO 2600 I=NA-1,1,-1
          BJ(I)=2.*I*BJ(I+1)/X-BJ(I+2)
2600    CONTINUE

      ELSE
!	Use the recurrence relation in the backward direction 
        N1=NA+MAX(25.D0,5.*SQRT(1.D0*NA))
        IF(X.LT.NA/6.0) N1=NA+MAX(25.D0,5.*SQRT(1.D0*NA))/LOG(NA*0.5/X)
        IF(NA.LT.X+15) N1=N1+15
        IF(N1-NA.LT.15) N1=15+NA
        BJ(N1+1)=0.0
        BJ(N1)=1.0
        DO 3200 I=N1-1,1,-1
          BJ(I)=2.*I*BJ(I+1)/X-BJ(I+2)
3200    CONTINUE

        S=BJ(1)
        DO 3400 I=3,N1,2
          S=S+2*BJ(I)
3400    CONTINUE

        DO 3600 I=1,NA+1
          BJ(I)=BJ(I)/S
3600    CONTINUE
!	If ABS(BJ(NA+2))<1/REPS, then the required accuracy may not be achieved
!	hence printout an error message
        IF(ABS(BJ(NA+2))*REPS.LT.1.0) PRINT 51,N,X,BJ(NA+2)
      ENDIF
 
!       IF(XB*N.LT.0.0) THEN    ! use this if the following doesn't work
      IF(N.LT.0.XOR.XB.LT.0.0) THEN
        DO 3800 I=2,NA+1,2
          BJ(I)=-BJ(I)
3800    CONTINUE
      ENDIF
      END
 
!     ----------------------------------------------------
 
!	To calculate the Bessel function of first and second kind of order zero
!	for real positive argument
!	For XB.LE.0 the function of second kind is not defined and subroutine
!	will return zero value  without any error message or flag.
!
!	XB : (input) Argument at which the functions need to be evaluated
!	BJ0 : (output) Calculated value of Bessel function of first kind
!	BY0 : (output) Calculated value of Bessel function of second kind
!
!	Required routines : None
 
      SUBROUTINE BJY0(XB,BJ0,BY0)
      IMPLICIT REAL*8(A-H,O-Z)
!	EUGAM is the Euler's constant
      PARAMETER(EUGAM=0.5772156649015328606D0,PI=3.14159265358979324D0)
      DIMENSION A(8),B(8),A1(5),B1(5),A2(5),B2(5),A0(7),B0(8)
 
!	Coefficients of rational function approximations
      DATA A/ 1.293686560051304152D-02,  8.573459862295151747D-05,
     1        3.854769244046149702D-07,  1.308534328117880493D-09,
     1        3.512360907188715842D-12,  7.512575042421009221D-15,
     1        1.229302278444845702D-17,  1.311883486088925264D-20/
      DATA B/ 9.999999999999999878D-01, -2.370631343994868513D-01,
     1        1.247651819849453565D-02, -2.529374255010058573D-04,
     1        2.411267406461247155D-06, -1.159484705672466498D-08,
     1        2.730546745501229851D-11, -2.517936655103065990D-14/

      DATA A0/ 1.089079731266387424D-02,  5.954632605213292419D-05,
     1         2.150109922530480401D-07,  5.641082188778387960D-10,
     1         1.102341761343675716D-12,  1.539990321465010920D-15,
     1         1.263081729204829936D-18/
      DATA B0/ 2.500000000000000006D-01, -2.071480067183403591D-02,
     1         5.553511120900719150D-04, -6.804373640943337406D-06,
     1         4.346149688717712144D-08, -1.505635199331021665D-10,
     1         2.703193275976574669D-13, -1.993047807317608951D-16/

      DATA A1/ 8.911849018950665793D+01,
     1         2.078818787053760203D+03,  1.366258799766718466D+04,
     1         1.800383785973922830D+04,  4.923440494847201509D+02/
      DATA B1/ 9.999999999999999908D-01,
     1         8.904817768950681616D+01,  2.072664795311476688D+03,
     1         1.352584337655551999D+04,  1.723138433448795889D+04/

      DATA A2/  1.046366195300779895D+02,  2.980727727381642723D+03,
     1          2.570418404044668245D+04,  5.291014161889741749D+04,
     1          9.497228391199055149D+03/
      DATA B2/ -1.249999999999999534D-01, -1.300633525376058087D+01,
     1         -3.651542590150084093D+02, -3.016744074771875522D+03,
     1         -5.251679479249748063D+03/
 
      X=ABS(XB)
      IF(X.LT.8.0) THEN
!	Use rational function approximations for 
        Y=X*X
        FN=((((((B(8)*Y+B(7))*Y+B(6))*Y+B(5))*Y+B(4))*Y+B(3))*Y+B(2))*Y
     1      +B(1)
        FD=(((((((A(8)*Y+A(7))*Y+A(6))*Y+A(5))*Y+A(4))*Y+A(3))*Y+A(2))*Y
     1      +A(1))*Y+1
        BJ0=FN/FD
 
        FN=((((((B0(8)*Y+B0(7))*Y+B0(6))*Y+B0(5))*Y+B0(4))*Y+B0(3))*Y
     1      +B0(2))*Y +B0(1)
        FD=((((((A0(7)*Y+A0(6))*Y+A0(5))*Y+A0(4))*Y+A0(3))*Y+A0(2))*Y
     1      +A0(1))*Y +1
        IF(X.GT.0.0) BY0=2.*(BJ0*(LOG(X/2)+EUGAM)+Y*FN/FD)/PI
 
      ELSE
!	Use rational function approximations for P_0 and Q_0
        X1=1./X**2
        FN=(((B1(5)*X1+B1(4))*X1+B1(3))*X1+B1(2))*X1+B1(1)
        FD=((((A1(5)*X1+A1(4))*X1+A1(3))*X1+A1(2))*X1+A1(1))*X1+1
        P0=FN/FD
 
        FN=(((B2(5)*X1+B2(4))*X1+B2(3))*X1+B2(2))*X1+B2(1)
        FD=((((A2(5)*X1+A2(4))*X1+A2(3))*X1+A2(2))*X1+A2(1))*X1+1
        Q0=FN/(FD*ABS(X))

        X1=ABS(X)
        BY0=SQRT(2./(PI*X1))*(P0*SIN(X1-PI/4)+Q0*COS(X1-PI/4))
        BJ0=SQRT(2./(PI*X1))*(P0*COS(X1-PI/4)-Q0*SIN(X1-PI/4))
      ENDIF

      IF(XB.LE.0) BY0=0.0
 
      END
 
!     ----------------------------------------------------
 
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
      IMPLICIT REAL*8(A-H,O-Z)
!	EUGAM is the Euler's constant
      PARAMETER(EUGAM=0.5772156649015328606D0,PI=3.14159265358979324D0)
      DIMENSION A(7),B(8),A1(5),B1(5),A2(5),B2(5),A0(7),B0(8)
 
!	Coefficients of rational function approximations
      DATA A/ 1.156878510067059849D-02,  6.749406787503859073D-05,
     1        2.614560837317535451D-07,  7.408815126464007290D-10,
     1        1.577637796406197189D-12,  2.432945305413635549D-15,
     1        2.257446839754248784D-18/
      DATA B/ 5.000000000000000074D-01, -5.671560744966475746D-02,
     1        1.914864631812345532D-03, -2.821407888958585592D-05,
     1        2.103168789640803591D-07, -8.322474383730280556D-10,
     1        1.678871778708754849D-12, -1.372424374400306547D-15/

      DATA A0/ 1.186963690270342970D-02,  7.123839029323002736D-05,
     1         2.848196759538669928D-07,  8.365303089083305885D-10,
     1         1.857096246589275584D-12,  3.012506935004947491D-15,
     1         2.996951174746838817D-18/
      DATA B0/ 2.500000000000000299D-01, -7.515759077432437273D-02,
     1         3.430771992327672576D-03, -6.022315614557372919D-05,
     1         5.067136874996839630D-07, -2.197514674456554803D-09,
     1         4.768619679411702014D-12, -4.139491442515065355D-15/

      DATA A1/ 8.659888261699365129D+01,
     1         1.932665751369749084D+03,  1.172714583536277145D+04,
     1         1.256737699073784218D+04, -6.147124347503755010D+02/
      DATA B1/ 1.000000000000000011D+00,
     1         8.671607011699346720D+01,  1.942669862370300601D+03,
     1         1.194181952104744095D+04,  1.371467864585746530D+04/

      DATA A2/  1.021472573795463627D+02,  2.807865400111916226D+03,
     1          2.280402060738415865D+04,  4.121116954504273053D+04,
     1          3.501974669280301705D+03/
      DATA B2/  3.749999999999999461D-01,  3.820268245483084309D+01,
     1          1.042753017477090289D+03,  8.289951986135169400D+03,
     1          1.371889615877945967D+04/
 
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
        X2=X1-0.75D0*PI
        BY1=SQRT(2./(PI*X1))*(P1*SIN(X2)+Q1*COS(X2))
        BJ1=SQRT(2./(PI*X1))*(P1*COS(X2)-Q1*SIN(X2))
      ENDIF

      IF(XB.LE.0.0) BY1=0.0
 
      END
 
!     ----------------------------------------------------
 
!	To calculate the modified Bessel function of first kind of order
!	zero for real argument
!
!	Required routines : None

      FUNCTION BI0(X)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(8),B(8),A1(9),B1(10)
 
!	Coefficients of rational function approximations
      DATA A/ -1.212778758454003780D-02,  7.395384394332586467D-05,
     1        -2.981475020389602859D-07,  8.749000589372218583D-10,
     1        -1.925066925538111917D-12,  3.116611043626875576D-15,
     1        -3.403826199284281585D-18,  1.919794284150553073D-21/
      DATA B/  9.999999999999999887D-01,  2.378722124154600447D-01,
     1         1.266700694780821586D-02,  2.627214102531807755D-04,
     1         2.599754169250266946D-06,  1.322628469861000483D-08,
     1         3.380151377715495065D-11,  3.502398414704805956D-14/
 
      DATA A1/ -4.368454479954000936D+01,  8.662328585154270157D+02,
     1         -9.865641284565305256D+03,  6.940667740497879456D+04,
     1         -3.022161789164256804D+05,  7.741268251525742446D+05,
     1         -1.016413121409283393D+06,  5.111938862294702143D+05,
     1         -5.159799972015923803D+04/
      DATA B1/  3.989422804014326451D-01, -1.737774413557676433D+01,
     1          3.434265110737608756D+02, -3.893820477243345847D+03,
     1          2.722034716144037839D+04, -1.173589448296665233D+05,
     1          2.954547271385017668D+05, -3.737076509558829360D+05,
     1          1.685151343593986274D+05, -9.954986502601715062D+03/
 
      XA=ABS(X)
      IF(XA.LT.8.0) THEN
!	Use rational function approximation
        Y=X*X
        FN=((((((B(8)*Y+B(7))*Y+B(6))*Y+B(5))*Y+B(4))*Y+B(3))*Y+B(2))*Y
     1      +B(1)
        FD=(((((((A(8)*Y+A(7))*Y+A(6))*Y+A(5))*Y+A(4))*Y+A(3))*Y+A(2))*Y
     1      +A(1))*Y+1
        BI0=FN/FD
 
      ELSE
!	Use rational function approximation to the asymptotic form
        Y=1./XA
        FN=((((((((B1(10)*Y+B1(9))*Y+B1(8))*Y+B1(7))*Y+B1(6))*Y+B1(5))*Y
     1      +B1(4))*Y+B1(3))*Y+B1(2))*Y+B1(1)
        FD=((((((((A1(9)*Y+A1(8))*Y+A1(7))*Y+A1(6))*Y+A1(5))*Y+A1(4))*Y+
     1      A1(3))*Y+A1(2))*Y+A1(1))*Y+1
        BI0=EXP(XA)*FN/(FD*SQRT(XA))
      ENDIF
 
      END
 
!     ----------------------------------------------------
 
!	To calculate the modified Bessel function of first kind of order
!	one for real argument
!
!	Required routines : None
 
      FUNCTION BI1(X)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(7),B(8),A1(9),B1(10)
 
!	Coefficients of rational function approximations
      DATA A/ -1.070229423222827071D-02,  5.648033546364466297D-05,
     1        -1.918252346392549379D-07,  4.558004979539079070D-10,
     1        -7.621816709677626962D-13,  8.342527256774311367D-16,
     1        -4.621507477175875649D-19/
      DATA B/  4.999999999999999914D-01,  5.714885288388592026D-02,
     1         1.963513444884162262D-03,  2.981702267497806226D-05,
     1         2.328548099504502597D-07,  9.862614413072000150D-10,
     1         2.192079167003372805D-12,  2.056894175269540721D-15/
 
      DATA A1/ -4.354982033508663805D+01,  8.607755018195304592D+02,
     1         -9.764617274576599078D+03,  6.833415275089609019D+04,
     1         -2.952159059017153952D+05,  7.462070230806506945D+05,
     1         -9.526843724015979184D+05,  4.423280833889137276D+05,
     1         -2.993689102971858935D+04/
      DATA B1/  3.989422804014327116D-01, -1.752346799070278231D+01,
     1          3.498681897994668814D+02, -4.022298493285640672D+03,
     1          2.868368405169782754D+04, -1.275734365963396269D+05,
     1          3.390186096946305841D+05, -4.802175709866444585D+05,
     1          2.931298028513980846D+05, -5.264175913253416995D+04/
 
      XA=ABS(X)
      IF(XA.LT.8.0) THEN
!	Use rational function approximation
        Y=X*X
        FN=((((((B(8)*Y+B(7))*Y+B(6))*Y+B(5))*Y+B(4))*Y+B(3))*Y+B(2))*Y
     1       +B(1)
        FD=((((((A(7)*Y+A(6))*Y+A(5))*Y+A(4))*Y+A(3))*Y+A(2))*Y+A(1))*Y
     1       +1
        BI1=X*FN/FD
 
      ELSE
!	Use rational function approximation to the asymptotic form
        Y=1./XA
        FN=((((((((B1(10)*Y+B1(9))*Y+B1(8))*Y+B1(7))*Y+B1(6))*Y+B1(5))*Y
     1       + B1(4))*Y+B1(3))*Y+B1(2))*Y+B1(1)
        FD=((((((((A1(9)*Y+A1(8))*Y+A1(7))*Y+A1(6))*Y+A1(5))*Y+A1(4))*Y+
     1       A1(3))*Y+A1(2))*Y+A1(1))*Y+1
        BI1=EXP(XA)*FN/(FD*SQRT(XA))
        IF(X.LT.0.0) BI1=-BI1
      ENDIF
 
      END
 
!     ----------------------------------------------------
 
!	To calculate the modified Bessel function of first kind of
!		integral order for real argument
!
!	N : (input) Order of Bessel function required, N must be positive
!	XB : (input) Argument at which the value is required
!	BI : (output) Real array of length at least N+16+MAX(25, 5*SQRT(N))
!		which will contain the value of Bessel function of order
!		0,1,...,N. BI(I+1) will contain Bessel function of order I
!		Remaining elements of array are used as scratch space
!
!	Required routines : BI0, BI1
!
 
      SUBROUTINE BIN(N,XB,BI)
      IMPLICIT REAL*8(A-H,O-Z)
!	REPS should be less than the machine accuracy
      PARAMETER(REPS=1.D-17)
      DIMENSION BI(*)
 
51    FORMAT('  BIN FAILED AT N =',I6,'  X =',1PD13.5,'  S =',D13.5)
 
      X=ABS(XB)
      NA=ABS(N)
      IF(XB.EQ.0.0) THEN
        BI(1)=1.0
        DO 1500 I=2,NA+1
          BI(I)=0.0
1500    CONTINUE
        RETURN
      ENDIF

      IF(NA.LT.X-10.OR.NA.LT.2) THEN
!	Use the recurrence relation in forward direction
 
        BI(1)=BI0(X)
        BI(2)=BI1(X)
        DO 2000 I=2,NA
          BI(I+1)=-2*(I-1)*BI(I)/X+BI(I-1)
2000    CONTINUE

      ELSE IF(X.LE.4.0) THEN

!	Use series expansion to calculate BI(NA), BI(NA-1)
        XA=X*X/4.0
        T0=X/2.0
        DO 2400 I=2,NA
          T0=T0*X/(2.0*I)
          IF(I.GE.NA-1) THEN
            T=T0
            S=T0
            DO 2200 J=1,50
              T=T*XA/(J*(J+I))
              S=S+T
              IF(ABS(T).LT.ABS(S*REPS)) GO TO 2300
2200        CONTINUE
2300        BI(I+1)=S
          ENDIF
2400    CONTINUE

        DO 2600 I=NA-1,1,-1
          BI(I)=2.*I*BI(I+1)/X+BI(I+2)
2600    CONTINUE
      ELSE
 
!	Use the recurrence relation in backward direction
        N1=NA+MAX(25.D0,5.*SQRT(1.D0*NA))
        IF(X.LT.NA/6.0) N1=NA+MAX(25.D0,5.*SQRT(1.D0*NA))/LOG(NA*0.5/X)
        IF(NA.LT.X+15) N1=N1+15
        IF(N1-NA.LT.15) N1=15+NA
        BI(N1+1)=0.0
        BI(N1)=1.0
        DO 3200 I=N1-1,1,-1
          BI(I)=2.*I*BI(I+1)/X+BI(I+2)
3200    CONTINUE

        S=BI(1)/BI0(X)
        DO 3600 I=1,NA+1
          BI(I)=BI(I)/S
3600    CONTINUE
!	If ABS(BI(NA+2))<1/REPS, then the required accuracy may not be achieved
!	hence printout an error message
        IF(ABS(BI(NA+2))*REPS.LT.1.0) PRINT 51,N,X,BI(NA+2)
      ENDIF
 
      IF(XB.LT.0.0) THEN
        DO 3800 I=2,NA+1,2
          BI(I)=-BI(I)
3800    CONTINUE
      ENDIF
      END
 
!     ----------------------------------------------------
 
!	To calculate the modified Bessel function of second kind of order
!	zero for real argument
!	For X<0 the function is not defined and the routine simply returns
!	a value of zero without any warning or error flag.
!
!	Required routines : None
 
      FUNCTION BK0(X)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(EUGAM=0.5772156649015328606D0)
      DIMENSION A(7),B(8),A1(5),B1(6),A0(8),B0(8)
 
!	Coefficients of rational function approximations
      DATA A/ -1.011123277211023982D-02,  5.027323349741487160D-05,
     1        -1.604132766768703653D-07,  3.571438627971404230D-10,
     1        -5.582724661115911120D-13,  5.702487045740630357D-16,
     1        -2.945917638250090849D-19/
      DATA B/  2.499999999999999993D-01,  2.090969180697244468D-02,
     1         5.713038828706075545D-04,  7.220998182565402322D-06,
     1         4.832471102493292948D-08,  1.789925692526897035D-10,
     1         3.530871144986696274D-13,  2.972558105712627660D-16/
 
      DATA A0/ -1.212778758454003780D-02,  7.395384394332586467D-05,
     1         -2.981475020389602859D-07,  8.749000589372218583D-10,
     1         -1.925066925538111917D-12,  3.116611043626875576D-15,
     1         -3.403826199284281585D-18,  1.919794284150553073D-21/
      DATA B0/  9.999999999999999887D-01,  2.378722124154600447D-01,
     1          1.266700694780821586D-02,  2.627214102531807755D-04,
     1          2.599754169250266946D-06,  1.322628469861000483D-08,
     1          3.380151377715495065D-11,  3.502398414704805956D-14/
 
      DATA A1/  1.134095488162070337D+01,  4.133208417436337182D+01,
     1          5.704583485964969346D+01,  2.690522251864472800D+01,
     1          2.767771031868579433D+00/
      DATA B1/  1.253314137315500251D+00,  1.405711481662802554D+01,
     1          5.011348668524983318D+01,  6.592870149979143648D+01,
     1          2.752549950796738039D+01,  1.796256223905248475D+00/
 
      BK0=0.0
      IF(X.LE.0.0) RETURN
      IF(X.LT.8.0) THEN
!	Use rational function approximation
        Y=X*X
        FN=((((((B0(8)*Y+B0(7))*Y+B0(6))*Y+B0(5))*Y+B0(4))*Y+B0(3))*Y
     1       +B0(2))*Y+B0(1)
        FD=(((((((A0(8)*Y+A0(7))*Y+A0(6))*Y+A0(5))*Y+A0(4))*Y+A0(3))*Y
     1       +A0(2))*Y+A0(1))*Y+1
        BI0=FN/FD
 
        FN=((((((B(8)*Y+B(7))*Y+B(6))*Y+B(5))*Y+B(4))*Y+B(3))*Y+B(2))*Y
     1       +B(1)
        FD=((((((A(7)*Y+A(6))*Y+A(5))*Y+A(4))*Y+A(3))*Y+A(2))*Y+A(1))*Y
     1       +1
        BK0=Y*FN/FD-BI0*(EUGAM+LOG(X/2))
 
      ELSE
!	Use rational function approximation to the asymptotic form
        X1=1./X
        FN=((((B1(6)*X1+B1(5))*X1+B1(4))*X1+B1(3))*X1+B1(2))*X1+B1(1)
        FD=((((A1(5)*X1+A1(4))*X1+A1(3))*X1+A1(2))*X1+A1(1))*X1+1
        BK0=EXP(-X)*FN/(FD*SQRT(X))
      ENDIF
 
      END
 
!     ----------------------------------------------------
 
!	To calculate the modified Bessel function of second kind of order
!	one for real argument
!	For X<0 the function is not defined and the routine simply returns
!	a value of zero without any warning or error flag.
!
!	Required routines : None
 
      FUNCTION BK1(X)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(EUGAM=0.5772156649015328606D0)
      DIMENSION A(7),B(8),A1(5),B1(6),A0(7),B0(8)
 
!	Coefficients of rational function approximations
      DATA A/ -1.070229423222827071D-02,  5.648033546364466297D-05,
     1        -1.918252346392549379D-07,  4.558004979539079070D-10,
     1        -7.621816709677626962D-13,  8.342527256774311367D-16,
     1        -4.621507477175875649D-19/
      DATA B/  4.999999999999999914D-01,  5.714885288388592026D-02,
     1         1.963513444884162262D-03,  2.981702267497806226D-05,
     1         2.328548099504502597D-07,  9.862614413072000150D-10,
     1         2.192079167003372805D-12,  2.056894175269540721D-15/
 
      DATA A0/ -1.097271232519001047D-02,  5.948570919145243152D-05,
     1         -2.079707688602524007D-07,  5.097946369393536825D-10,
     1         -8.813373771643053620D-13,  9.993295069392185250D-16,
     1         -5.743938665570093767D-19/
      DATA B0/  2.499999999999999663D-01,  7.538182191870271423D-02,
     1          3.497906054669938529D-03,  6.321709197334740028D-05,
     1          5.569209160486120478D-07,  2.585145364373725340D-09,
     1          6.185681870407510042D-12,  6.175435537988524458D-15/
 
      DATA A1/  9.546087640477785968D+00,  2.842458429353285707D+01,
     1          3.039341689594721892D+01,  9.857685145186666128D+00,
     1          4.376813258161116587D-01/
      DATA B1/  1.253314137315500259D+00,  1.243423939735687111D+01,
     1          3.996465306538935273D+01,  5.017830258820480688D+01,
     1          2.351074680068346432D+01,  2.993271886379231231D+00/
 
      BK1=0.0
      IF(X.LE.0.0) RETURN
      IF(X.LT.8.0) THEN
!	Use rational function approximation
        Y=X*X
        FN=((((((B(8)*Y+B(7))*Y+B(6))*Y+B(5))*Y+B(4))*Y+B(3))*Y+B(2))*Y
     1      +B(1)
        FD=((((((A(7)*Y+A(6))*Y+A(5))*Y+A(4))*Y+A(3))*Y+A(2))*Y+A(1))*Y
     1      +1
 
        BI1=X*FN/FD
        FN=((((((B0(8)*Y+B0(7))*Y+B0(6))*Y+B0(5))*Y+B0(4))*Y+B0(3))*Y
     1      +B0(2))*Y+B0(1)
        FD=((((((A0(7)*Y+A0(6))*Y+A0(5))*Y+A0(4))*Y+A0(3))*Y+A0(2))*Y
     1      +A0(1))*Y+1
        BK1=-X*FN/FD+BI1*(EUGAM+LOG(X/2))+1./X
 
      ELSE
!	Use rational function approximation to the asymptotic form
        X1=1./X
        FN=((((B1(6)*X1+B1(5))*X1+B1(4))*X1+B1(3))*X1+B1(2))*X1+B1(1)
        FD=((((A1(5)*X1+A1(4))*X1+A1(3))*X1+A1(2))*X1+A1(1))*X1+1
        BK1=EXP(-X)*FN/(FD*SQRT(X))
      ENDIF
 
      END
 
!     ----------------------------------------------------
 
!	To calculate the modified Bessel function of second kind of
!		integral order for real argument
!		For X.LE.0 the function is not defined and in this case
!		no calculations are done and no warning is issued.
!
!	N : (input) Order of Bessel function required, N must be positive
!	X : (input) Argument at which the value is required
!	BK : (output) Real array of length at least ABS(N)+1
!		which will contain the value of Bessel function of order
!		0,1,...,N. BK(I+1) will contain Bessel function of order I
!
!	Required routines : BK0, BK1
!
 
      SUBROUTINE BKN(N,X,BK)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION BK(N+1)
 
      BK(1)=BK0(X)
      BK(2)=BK1(X)
      IF(X.LE.0.0) RETURN
!	Use the recurrence relation in Forward direction
      DO 2000 I=2,ABS(N)
        BK(I+1)=2*(I-1)*BK(I)/X+BK(I-1)
2000  CONTINUE
      END
 
!     ----------------------------------------------------
 
!	To calculate the Bessel function of second kind of order zero
!	for real positive argument
!	For X.LE.0 the function is not defined and subroutine will return
!	a zero value without any error message or flag.
!
!	Required routines : None
 
      FUNCTION BY0(X)
      IMPLICIT REAL*8(A-H,O-Z)
!	EUGAM is the Euler's constant
      PARAMETER(EUGAM=0.5772156649015328606D0,PI=3.14159265358979324D0)
      DIMENSION A(8),B(8),A1(5),B1(5),A2(5),B2(5),A0(7),B0(8)
 
!	Coefficients of rational function approximations
      DATA A/ 1.293686560051304152D-02,  8.573459862295151747D-05,
     1        3.854769244046149702D-07,  1.308534328117880493D-09,
     1        3.512360907188715842D-12,  7.512575042421009221D-15,
     1        1.229302278444845702D-17,  1.311883486088925264D-20/
      DATA B/ 9.999999999999999878D-01, -2.370631343994868513D-01,
     1        1.247651819849453565D-02, -2.529374255010058573D-04,
     1        2.411267406461247155D-06, -1.159484705672466498D-08,
     1        2.730546745501229851D-11, -2.517936655103065990D-14/

      DATA A0/ 1.089079731266387424D-02,  5.954632605213292419D-05,
     1         2.150109922530480401D-07,  5.641082188778387960D-10,
     1         1.102341761343675716D-12,  1.539990321465010920D-15,
     1         1.263081729204829936D-18/
      DATA B0/ 2.500000000000000006D-01, -2.071480067183403591D-02,
     1         5.553511120900719150D-04, -6.804373640943337406D-06,
     1         4.346149688717712144D-08, -1.505635199331021665D-10,
     1         2.703193275976574669D-13, -1.993047807317608951D-16/

      DATA A1/ 8.911849018950665793D+01,
     1         2.078818787053760203D+03,  1.366258799766718466D+04,
     1         1.800383785973922830D+04,  4.923440494847201509D+02/
      DATA B1/ 9.999999999999999908D-01,
     1         8.904817768950681616D+01,  2.072664795311476688D+03,
     1         1.352584337655551999D+04,  1.723138433448795889D+04/

      DATA A2/  1.046366195300779895D+02,  2.980727727381642723D+03,
     1          2.570418404044668245D+04,  5.291014161889741749D+04,
     1          9.497228391199055149D+03/
      DATA B2/ -1.249999999999999534D-01, -1.300633525376058087D+01,
     1         -3.651542590150084093D+02, -3.016744074771875522D+03,
     1         -5.251679479249748063D+03/
 
      BY0=0.0
      IF(X.LE.0) RETURN

      IF(X.LT.8.0) THEN
!	Use rational function approximations for 
        Y=X*X
        FN=((((((B(8)*Y+B(7))*Y+B(6))*Y+B(5))*Y+B(4))*Y+B(3))*Y+B(2))*Y
     1      +B(1)
        FD=(((((((A(8)*Y+A(7))*Y+A(6))*Y+A(5))*Y+A(4))*Y+A(3))*Y+A(2))*Y
     1      +A(1))*Y+1
        BJ0=FN/FD
 
        FN=((((((B0(8)*Y+B0(7))*Y+B0(6))*Y+B0(5))*Y+B0(4))*Y+B0(3))*Y
     1      +B0(2))*Y +B0(1)
        FD=((((((A0(7)*Y+A0(6))*Y+A0(5))*Y+A0(4))*Y+A0(3))*Y+A0(2))*Y
     1      +A0(1))*Y +1
        BY0=2.*(BJ0*(LOG(X/2)+EUGAM)+Y*FN/FD)/PI
 
      ELSE
!	Use rational function approximations for P_0 and Q_0
        X1=1./X**2
        FN=(((B1(5)*X1+B1(4))*X1+B1(3))*X1+B1(2))*X1+B1(1)
        FD=((((A1(5)*X1+A1(4))*X1+A1(3))*X1+A1(2))*X1+A1(1))*X1+1
        P0=FN/FD
 
        FN=(((B2(5)*X1+B2(4))*X1+B2(3))*X1+B2(2))*X1+B2(1)
        FD=((((A2(5)*X1+A2(4))*X1+A2(3))*X1+A2(2))*X1+A2(1))*X1+1
        Q0=FN/(FD*ABS(X))

        X1=ABS(X)
        BY0=SQRT(2./(PI*X1))*(P0*SIN(X1-PI/4)+Q0*COS(X1-PI/4))
!        BJ0=SQRT(2./(PI*X1))*(P0*COS(X1-PI/4)-Q0*SIN(X1-PI/4))
      ENDIF
 
      END
 
!     ----------------------------------------------------
 
!	To calculate the Bessel function of second kind of order one
!	for real positive argument
!	For X.LE.0 the function is not defined and subroutine will return
!	a zero value without any error message or flag.
!
!	Required routines : None
 
      FUNCTION BY1(X)
      IMPLICIT REAL*8(A-H,O-Z)
!	EUGAM is the Euler's constant
      PARAMETER(EUGAM=0.5772156649015328606D0,PI=3.14159265358979324D0)
      DIMENSION A(7),B(8),A1(5),B1(5),A2(5),B2(5),A0(7),B0(8)
 
!	Coefficients of rational function approximations
      DATA A/ 1.156878510067059849D-02,  6.749406787503859073D-05,
     1        2.614560837317535451D-07,  7.408815126464007290D-10,
     1        1.577637796406197189D-12,  2.432945305413635549D-15,
     1        2.257446839754248784D-18/
      DATA B/ 5.000000000000000074D-01, -5.671560744966475746D-02,
     1        1.914864631812345532D-03, -2.821407888958585592D-05,
     1        2.103168789640803591D-07, -8.322474383730280556D-10,
     1        1.678871778708754849D-12, -1.372424374400306547D-15/

      DATA A0/ 1.186963690270342970D-02,  7.123839029323002736D-05,
     1         2.848196759538669928D-07,  8.365303089083305885D-10,
     1         1.857096246589275584D-12,  3.012506935004947491D-15,
     1         2.996951174746838817D-18/
      DATA B0/ 2.500000000000000299D-01, -7.515759077432437273D-02,
     1         3.430771992327672576D-03, -6.022315614557372919D-05,
     1         5.067136874996839630D-07, -2.197514674456554803D-09,
     1         4.768619679411702014D-12, -4.139491442515065355D-15/

      DATA A1/ 8.659888261699365129D+01,
     1         1.932665751369749084D+03,  1.172714583536277145D+04,
     1         1.256737699073784218D+04, -6.147124347503755010D+02/
      DATA B1/ 1.000000000000000011D+00,
     1         8.671607011699346720D+01,  1.942669862370300601D+03,
     1         1.194181952104744095D+04,  1.371467864585746530D+04/

      DATA A2/  1.021472573795463627D+02,  2.807865400111916226D+03,
     1          2.280402060738415865D+04,  4.121116954504273053D+04,
     1          3.501974669280301705D+03/
      DATA B2/  3.749999999999999461D-01,  3.820268245483084309D+01,
     1          1.042753017477090289D+03,  8.289951986135169400D+03,
     1          1.371889615877945967D+04/
 
      BY1=0.0 
      IF(X.LE.0) RETURN
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
        BY1=2.*(BJ1*(LOG(X/2)+EUGAM)-1./X-X*FN/FD)/PI
 
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
        X2=X1-0.75D0*PI
        BY1=SQRT(2./(PI*X1))*(P1*SIN(X2)+Q1*COS(X2))
!        BJ1=SQRT(2./(PI*X1))*(P1*COS(X2)-Q1*SIN(X2))
      ENDIF
 
      END
 
!     ----------------------------------------------------
 
!	To calculate the Bessel function of second kind of integral order
!		for real argument
!
!	N : (input) Order of Bessel function required, N may be negative
!		or positive. For N=0,1 use BY0 and BY1 respectively.
!	X : (input) Argument at which the value is required
!	BY : (output) Real array of length at least ABS(N)+1
!		which will contain the value of Bessel function of order
!		0,1,...,ABS(N). BY(I+1) will contain Bessel function of order I
!		or -I (if N<0).
!
!	Required routines : BY0, BY1
!
 
      SUBROUTINE BYN(N,X,BY)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION BY(N+1)
 
      BY(1)=BY0(X)
      BY(2)=BY1(X)
      IF(X.LE.0.0) RETURN
!	Use the recurrence relation in Forward direction
      DO 2000 I=2,ABS(N)
        BY(I+1)=2*(I-1)*BY(I)/X-BY(I-1)
2000  CONTINUE
 
      IF(N.LT.0) THEN
        DO 2500 I=2,ABS(N)+1,2
          BY(I)=-BY(I)
2500    CONTINUE
      ENDIF
      END
 
!     ----------------------------------------------------
 
!	To calculate the spherical Bessel function of integral order
!	(j_n(x)=Sqrt(PI/(2x))*J_{n+1/2}(x)) for a real argument
!
!	N : (input) Order of Bessel function required, N may be negative
!		or positive. 
!	XB : (input) Argument at which the value is required
!	BJ : (output) Real array of length at least 
!		ABS(N)+16+MAX(25,5*SQRT(N))
!		which will contain the value of Bessel function of order
!		0,1,...,ABS(N). BJ(I+1) will contain Bessel function of order I
!		or -I (if N<0)
!		Remaining elements of array are used as scratch space
!
!	Required routines : None
 
      SUBROUTINE SPHBJN(N,XB,BJ)
      IMPLICIT REAL*8(A-H,O-Z)
!	REPS should be less than the machine accuracy
      PARAMETER(REPS=1.D-17)
      DIMENSION BJ(*)
 
51    FORMAT(' SPHBJN FAILED AT N =',I6,'  X =',1PD13.5,'  S =',D13.5)
 
      X=ABS(XB)
      NA=ABS(N)
      IF(XB.EQ.0) THEN
        BJ(1)=1.0
        DO 1800 I=2,NA+1
          BJ(I)=0.0
1800    CONTINUE
        RETURN
      ENDIF
      IF(N.GT.0) THEN
        IF(N.LT.X.OR.N.LT.2) THEN
!	Use the recurrence relation in the forward direction 
 
          BJ(1)=SIN(X)/X
          BJ(2)=SIN(X)/X**2-COS(X)/X
          DO 2000 I=2,N
            BJ(I+1)=(2*I-1)*BJ(I)/X-BJ(I-1)
2000      CONTINUE

        ELSE IF(X.LE.4.0) THEN

!	Use series expansion to calculate  BJ(NA), BJ(NA-1)
          XA=X*X/4.0
          T0=1.0
          DO 2400 I=1,NA
            T0=T0*X/(2.0*I+1)
            IF(I.GE.NA-1) THEN
              T=T0
              S=T0
              DO 2200 J=1,50
                T=-T*XA/(J*(J+I+0.5D0))
                S=S+T
                IF(ABS(T).LT.ABS(S*REPS)) GO TO 2300
2200          CONTINUE
2300          BJ(I+1)=S
            ENDIF
2400      CONTINUE
 
          DO 2600 I=NA-1,1,-1
            BJ(I)=(2.*I+1)*BJ(I+1)/X-BJ(I+2)
2600      CONTINUE
        ELSE
 
!	Use the recurrence relation in the backward direction 
          N1=NA+MAX(25.D0,5.*SQRT(1.D0*NA))
          IF(X.LT.NA/6.0) N1=NA+MAX(25.D0,5.*SQRT(1D0*NA))/LOG(NA*0.5/X)
          IF(NA.LT.X+15) N1=N1+15
          IF(N1-NA.LT.15) N1=15+NA
          BJ(N1+1)=0.0
          BJ(N1)=1.0
          DO 3200 I=N1-1,1,-1
            BJ(I)=(2.*I+1)*BJ(I+1)/X-BJ(I+2)
3200      CONTINUE

          S=BJ(1)*X/SIN(X)
          DO 3400 I=1,NA+1
            BJ(I)=BJ(I)/S
3400      CONTINUE

!	If ABS(BJ(NA+2))<1/REPS, then the required accuracy may not be achieved
          IF(ABS(BJ(NA+2))*REPS.LT.1.0) PRINT 51,N,X,BJ(NA+2)
        ENDIF
 
      ELSE
 
!	For negative N use the recurrence relation in the forward direction 
        BJ(1)=SIN(X)/X
        BJ(2)=COS(X)/X
        DO 3600 I=2,ABS(N)
          BJ(I+1)=(-2*I+3)*BJ(I)/X-BJ(I-1)
3600    CONTINUE
      ENDIF
 
      IF(XB.LT.0.0) THEN
        DO 4000 I=2,ABS(N)+1,2
          BJ(I)=-BJ(I)
4000    CONTINUE
      ENDIF
      END
 
