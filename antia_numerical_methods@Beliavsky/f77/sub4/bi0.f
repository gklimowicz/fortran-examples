!	To calculate the modified Bessel function of first kind of order
!	zero for real argument
!
!	Required routines : None

      FUNCTION BI0(X)
!      IMPLICIT REAL*8(A-H,O-Z)
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