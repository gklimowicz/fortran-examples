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
      IMPLICIT REAL*16(A-H,O-Z)
!	EUGAM is the Euler's constant
      PARAMETER(EUGAM=0.5772156649015328606Q0,PI=3.14159265358979324Q0)
      DIMENSION A(8),B(8),A1(5),B1(5),A2(5),B2(5),A0(7),B0(8)
 
!	Coefficients of rational function approximations
      DATA A/ 1.293686560051304152Q-02,  8.573459862295151747Q-05,
     1        3.854769244046149702Q-07,  1.308534328117880493Q-09,
     1        3.512360907188715842Q-12,  7.512575042421009221Q-15,
     1        1.229302278444845702Q-17,  1.311883486088925264Q-20/
      DATA B/ 9.999999999999999878Q-01, -2.370631343994868513Q-01,
     1        1.247651819849453565Q-02, -2.529374255010058573Q-04,
     1        2.411267406461247155Q-06, -1.159484705672466498Q-08,
     1        2.730546745501229851Q-11, -2.517936655103065990Q-14/

      DATA A0/ 1.089079731266387424Q-02,  5.954632605213292419Q-05,
     1         2.150109922530480401Q-07,  5.641082188778387960Q-10,
     1         1.102341761343675716Q-12,  1.539990321465010920Q-15,
     1         1.263081729204829936Q-18/
      DATA B0/ 2.500000000000000006Q-01, -2.071480067183403591Q-02,
     1         5.553511120900719150Q-04, -6.804373640943337406Q-06,
     1         4.346149688717712144Q-08, -1.505635199331021665Q-10,
     1         2.703193275976574669Q-13, -1.993047807317608951Q-16/

      DATA A1/ 8.911849018950665793Q+01,
     1         2.078818787053760203Q+03,  1.366258799766718466Q+04,
     1         1.800383785973922830Q+04,  4.923440494847201509Q+02/
      DATA B1/ 9.999999999999999908Q-01,
     1         8.904817768950681616Q+01,  2.072664795311476688Q+03,
     1         1.352584337655551999Q+04,  1.723138433448795889Q+04/

      DATA A2/  1.046366195300779895Q+02,  2.980727727381642723Q+03,
     1          2.570418404044668245Q+04,  5.291014161889741749Q+04,
     1          9.497228391199055149Q+03/
      DATA B2/ -1.249999999999999534Q-01, -1.300633525376058087Q+01,
     1         -3.651542590150084093Q+02, -3.016744074771875522Q+03,
     1         -5.251679479249748063Q+03/
 
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
