!	To calculate the modified Bessel function of first kind of order
!	zero for real argument
!
!	Required routines : None

      FUNCTION BI0(X)
      IMPLICIT REAL*16(A-H,O-Z)
      DIMENSION A(8),B(8),A1(9),B1(10)
 
!	Coefficients of rational function approximations
      DATA A/ -1.212778758454003780Q-02,  7.395384394332586467Q-05,
     1        -2.981475020389602859Q-07,  8.749000589372218583Q-10,
     1        -1.925066925538111917Q-12,  3.116611043626875576Q-15,
     1        -3.403826199284281585Q-18,  1.919794284150553073Q-21/
      DATA B/  9.999999999999999887Q-01,  2.378722124154600447Q-01,
     1         1.266700694780821586Q-02,  2.627214102531807755Q-04,
     1         2.599754169250266946Q-06,  1.322628469861000483Q-08,
     1         3.380151377715495065Q-11,  3.502398414704805956Q-14/
 
      DATA A1/ -4.368454479954000936Q+01,  8.662328585154270157Q+02,
     1         -9.865641284565305256Q+03,  6.940667740497879456Q+04,
     1         -3.022161789164256804Q+05,  7.741268251525742446Q+05,
     1         -1.016413121409283393Q+06,  5.111938862294702143Q+05,
     1         -5.159799972015923803Q+04/
      DATA B1/  3.989422804014326451Q-01, -1.737774413557676433Q+01,
     1          3.434265110737608756Q+02, -3.893820477243345847Q+03,
     1          2.722034716144037839Q+04, -1.173589448296665233Q+05,
     1          2.954547271385017668Q+05, -3.737076509558829360Q+05,
     1          1.685151343593986274Q+05, -9.954986502601715062Q+03/
 
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
