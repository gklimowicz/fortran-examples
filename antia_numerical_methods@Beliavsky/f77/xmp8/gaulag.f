!     PROGRAM FOR INTEGRATION OVER INFINITE REGION USING GAUSSIAN FORMULAS

      PROGRAM INTEG
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL F,F2

!     EXAMPLE 6.10

51    FORMAT('   IER =',I4,5X,'A =',1PD14.6,5X,'A1 =',D14.6/
     1       5X,'NO. OF FUNCTION',
     2       ' EVALUATIONS =',I7/5X,'INTEGRAL =',D14.6,5X,
     3       'ESTIMATED ERROR =',D14.6)
52    FORMAT('   INITIAL VALUE OF A1 =',1PD14.6)

      REPS=1.D-12
      AEPS=1.D-13

100   PRINT *,'TYPE  A=LOWER LIMIT,  A1=POINT AT WHICH INTEGRAL IS TO'
     1       , ' BE SPLIT'
      PRINT *,'                    (QUITS WHEN A<-100)'
      READ *,A,A1
      IF(A.LT.-100) STOP
      WRITE(6,52) A1
      CALL GAULAG(RI,A,A1,REPS,AEPS,DIF,F,F2,N,IER)
      WRITE(6,51) IER,A,A1,N,RI,DIF
      GO TO 100
      END
 
!     ---------------------------------------------------
 
!	To integrate a function over semi-infinite interval using a
!	combination of Gauss-Legendre and Gauss-Laguerre formulas
!
!	RINT : (output) Calculated value of the integral
!	A : (input) The lower limit
!	A1 : (input/output) The point at which integral has to be broken
!		A1 will be adjusted by the subroutine.
!	REPS : (input) The required relative accuracy
!	AEPS : (input) The required absolute accuracy
!		The estimated error should be less than MAX(AEPS,REPS*ABS(RINT))
!	DIF : (output) estimated (absolute) error achieved by the subroutine
!	F : (input) Name of the function routine to calculate the integrand
!	F2 : (input) Name of the function routine to calculate F(X)*EXP(X)
!	NP : (output) Number of function evaluations used by the subroutine
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=30 implies specified accuracy was not achieved by GAUSS
!		IER=37 implies specified accuracy was not achieved by LAGURE
!		IER=38 implies specified accuracy was not achieved by both
!		GAUSS and LAGURE
!		In all cases DIF will contain the estimated accuracy
!
!		FUNCTION F(X) and F2(X) must be supplied by the user.
!
!	Required routines : GAUSS, LAGURE, F, F2
!
      SUBROUTINE GAULAG(RINT,A,A1,REPS,AEPS,DIF,F,F2,NP,IER)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(AMAX=50.)
      EXTERNAL F,F2

      IER1=0
      R1=0.0
      R2=0.0
      DIF=0.0
      NP=0
      NPT=0
      IF(A1.LT.A) A1=A

!	To calculate integral over [A1,Infinity)
2200  CALL LAGURE(R1,A1,AEPS,REPS,DIF1,F2,NPT1,IER)
      NP=NP+NPT1
      IF(IER.EQ.0) GO TO 2500
!	If LAGURE fails then increase A1
      T1=A1
      A1=MAX(A1*2.,A1+2.)
      IF(A1.LT.AMAX) GO TO 2200
!	To avoid the possibility of getting in infinite loop A1 is not
!	increased beyond AMAX
      IER1=37
      IER=0
      A1=T1

!	To calculate integral over [A,A1]
2500  IF(A1-A.GT.AEPS) CALL GAUSS(R2,A,A1,16,REPS,AEPS,DIF,IER,NPT,F)
      IER=IER+IER1
      IF(IER.GT.IER1.AND.IER1.GT.0) IER=38
      RINT=R1+R2
      DIF=ABS(DIF)+ABS(DIF1)
      NP=NP+NPT
      END
 
!     -----------------------------------------------
 
!	To integrate a function over finite interval using Gauss-Legendre formulas
!
!	RINT : (output) Calculated value of the integral
!	A : (input) The lower limit
!	B : (input) The upper limit
!	NP : (input/output) The Gauss-Legendre formula to be used, (NP=2,4,8,16,32)
!		For other values of NP it will be set to 8.
!		Subroutine will use composite NP-point formula
!	REPS : (input) The required relative accuracy
!	AEPS : (input) The required absolute accuracy
!		The estimated error should be less than MAX(AEPS,REPS*ABS(RINT))
!	DIF : (output) estimated (absolute) error achieved by the subroutine
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=30 implies specified accuracy was not achieved
!			DIF will contain the estimated accuracy
!		IER=36 implies that NP was not 2,4,8,16 or 32. In which case
!			it is set to 8.
!	NPT : (output) Number of function evaluations used by the subroutine
!	FUN : (input) Name of the function routine to calculate the integrand
!		FUNCTION FUN(X) must be supplied by the user.
!
!	Required routines : FUN

      SUBROUTINE GAUSS(RINT,A,B,NP,REPS,AEPS,DIF,IER,NPT,FUN)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NMAX=9)
      DIMENSION W(31),X(31)

!	Weights and abscissas for Gauss-Legendre quadrature.
!	For N-point formula W(K)=W(N-K+1) and X(K)=-X(N-K+1)
!		For K=1,2,...,N/2. Hence only half points are tabulated.
!	For 2-point W(1); 4-point W(2), W(3); 8-point W(4),...,W(7);
!	16-point W(8),...,W(15); 32-point W(16),...,W(31) are the
!	weights corresponding to abscissas X(I).

      DATA W/1.0D0,
     1       0.34785484513745385737D0, 0.65214515486254614263D0,
     2       0.10122853629037625915D0, 0.22238103445337447054D0,
     3       0.31370664587788728734D0, 0.36268378337836198297D0,
     4       0.02715245941175409485D0, 0.06225352393864789286D0,
     5       0.09515851168249278481D0, 0.12462897125553387205D0,
     6       0.14959598881657673208D0, 0.16915651939500253819D0,
     7       0.18260341504492358887D0, 0.18945061045506849629D0,
     8       0.00701861000947009660D0, 0.01627439473090567061D0,
     9       0.02539206530926205945D0, 0.03427386291302143310D0,
     1       0.04283589802222668066D0, 0.05099805926237617620D0,
     2       0.05868409347853554714D0, 0.06582222277636184684D0,
     3       0.07234579410884850623D0, 0.07819389578707030647D0,
     4       0.08331192422694675522D0, 0.08765209300440381114D0,
     4       0.09117387869576388471D0, 0.09384439908080456564D0,
     5       0.09563872007927485942D0, 0.09654008851472780057D0/

      DATA X/0.57735026918962576451D0,
     1       0.86113631159405257522D0, 0.33998104358485626480D0,
     2       0.96028985649753623168D0, 0.79666647741362673959D0,
     3       0.52553240991632898582D0, 0.18343464249564980494D0,
     4       0.98940093499164993260D0, 0.94457502307323257608D0,
     5       0.86563120238783174388D0, 0.75540440835500303390D0,
     6       0.61787624440264374845D0, 0.45801677765722738634D0,
     7       0.28160355077925891323D0, 0.09501250983763744019D0,
     8       0.99726386184948156354D0, 0.98561151154526833540D0,
     9       0.96476225558750643077D0, 0.93490607593773968917D0,
     1       0.89632115576605212397D0, 0.84936761373256997013D0,
     2       0.79448379596794240696D0, 0.73218211874028968039D0,
     3       0.66304426693021520098D0, 0.58771575724076232904D0,
     4       0.50689990893222939002D0, 0.42135127613063534536D0,
     5       0.33186860228212764978D0, 0.23928736225213707454D0,
     6       0.14447196158279649349D0, 0.04830766568773831623D0/

      N=1
      DX=B-A
      IER=0
      RINT=0.0
      NPT=0

      NO=-1
      IF(NP.EQ.2) NO=1
      IF(NP.EQ.4) NO=2
      IF(NP.EQ.8) NO=4
      IF(NP.EQ.16) NO=8
      IF(NP.EQ.32) NO=16
      IF(NO.LT.0) THEN
!	If NP-point formula is not available use NP=8
        NP=8
        NO=4
        IER=36
      ENDIF
!	X(NO),...,X(NP2) are the abscissas for the formula
      NP2=NO+NP/2-1

!	Subdivide the interval until convergence
      DO 5000 I=1,NMAX

        R1=0.0
        DO 3000 J=1,N
          A1=A+(J-1)*DX
          AT=DX/2.
          BT=A1+DX/2.
!	To reduce roundoff errors sum over each subinterval is evaluated separately
          S1=0.0
          DO 2000 K=NO,NP2
2000      S1=S1+W(K)*(FUN(AT*X(K)+BT)+FUN(BT-AT*X(K)))
          R1=R1+S1
3000    CONTINUE
        R1=R1*DX/2.

!	convergence check
        DIF=R1-RINT
        RINT=R1
        NPT=NPT+N*NP
        IF(I.GT.1.AND.ABS(DIF).LT.MAX(AEPS,ABS(RINT)*REPS)) RETURN
        DX=DX/2.
        N=N*2
5000  CONTINUE

!	Integral fails to converge
      IER=30
      END
 
!     -------------------------------------------------------

!	To integrate a function over semi-infinite interval using Gauss-Laguerre formulas
!
!	RINT : (output) Calculated value of the integral
!	A : (input) The lower limit
!	AEPS : (input) The required absolute accuracy
!	REPS : (input) The required relative accuracy
!		The estimated error should be less than MAX(AEPS,REPS*ABS(RINT))
!	DIF : (output) estimated (absolute) error achieved by the subroutine
!	F : (input) Name of the function routine to calculate the
!		integrand (multiplied by EXP(X))
!	NPT : (output) Number of function evaluations used by the subroutine
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=30 implies specified accuracy was not achieved
!			DIF will contain the estimated accuracy
!
!	Function F(X) must be supplied by the user
!
!	Required routines : F

      SUBROUTINE LAGURE(RINT,A,AEPS,REPS,DIF,F,NPT,IER)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION W(62),X(62)

!	Weights and abscissas for Gauss-Laguerre quadrature
!	W(N-1),...,W(2N-2), are the weights for N-point rule and
!	X(N-1),...,X(2N-2), the corresponding abscissas
!	Weights and abscissas are available for N=2,4,8,16,32

      DATA (X(I),I=1,30)/   0.585786437627D0,     3.414213562373D0,
     *.322547689619392312D0,1.74576110115834658D0,4.53662029692112798D0,
     *9.39507091230113313D0,.170279632305101000D0,.903701776799379912D0,
     *2.25108662986613069D0,4.26670017028765879D0,7.04590540239346570D0,
     *10.7585160101809952D0,15.7406786412780046D0,22.8631317368892641D0,
     *.087649410478927840D0,.462696328915080832D0,1.14105777483122686D0,
     *2.12928364509838062D0,3.43708663389320665D0,5.07801861454976791D0,
     *7.07033853504823413D0,9.43831433639193878D0,12.2142233688661587D0,
     *15.4415273687816171D0,19.1801568567531349D0,23.5159056939919085D0,
     *28.5787297428821404D0,34.5833987022866258D0,41.9404526476883326D0,
     *51.7011603395433184D0/

!	Some of these weights will give underflow in REAL*4 arithmetic
!	If necessary replace the last two numbers by zero to avoid warning
      DATA (X(I),I=31,62)/0.0444893658332670184D0,.234526109519618537D0,
     *.576884629301886426D0,1.07244875381781763D0,1.72240877644464544D0,
     *2.52833670642579488D0,3.49221327302199449D0,4.61645676974976739D0,
     *5.90395850417424395D0,7.35812673318624111D0,8.98294092421259610D0,
     *10.7830186325399721D0,12.7636979867427251D0,14.9311397555225573D0,
     *17.2924543367153148D0,19.8558609403360547D0,22.6308890131967745D0,
     *25.6286360224592478D0,28.8621018163234747D0,32.3466291539647370D0,
     *36.1004948057519738D0,40.1457197715394415D0,44.5092079957549380D0,
     *49.2243949873086392D0,54.3337213333969073D0,59.8925091621340182D0,
     *65.9753772879350528D0,72.6876280906627086D0,80.1874469779135231D0,
     *88.7353404178923987D0,98.8295428682839726D0,111.751398097937695D0/
      DATA (W(I),I=1,30)/    0.853553390593D0,     0.146446609407D0,
     *.603154104341633602D0,.357418692437799687D0,.38887908515005384D-1,
     *.53929470556132745D-3,.369188589341637530D0,.418786780814342956D0,
     *.175794986637171806D0,.33343492261215651D-1,.27945362352256725D-2,
     *.90765087733582131D-4,.84857467162725315D-6,.10480011748715104D-8,
     *.206151714957800994D0,.331057854950884166D0,.265795777644214153D0,
     *.136296934296377540D0,.47328928694125219D-1,.11299900080339453D-1,
     *.18490709435263109D-2,.20427191530827846D-3,.14844586873981299D-4,
     *.68283193308711996D-6,.18810248410796732D-7,.28623502429738816D-9,
     *.2127079033224103D-11,.6297967002517868D-14,.5050473700035513D-17,
     *.4161462370372855D-21/
      DATA (W(I),I=31,62)/0.109218341952384971D0,0.210443107938813234D0,
     *.235213229669848005D0,.195903335972881043D0,.129983786286071761D0,
     *.70578623865717442D-1,.31760912509175070D-1,.11918214834838557D-1,
     *.37388162946115248D-2,.98080330661495513D-3,.21486491880136419D-3,
     *.39203419679879472D-4,.59345416128686329D-5,.74164045786675522D-6,
     *.76045678791207815D-7,.63506022266258067D-8,.42813829710409289D-9,
     *.2305899491891336D-10,.9799379288727094D-12,.3237801657729266D-13,
     *.8171823443420719D-15,.1542133833393823D-16,.2119792290163619D-18,
     *.2054429673788045D-20,.1346982586637395D-22,.5661294130397359D-25,
     *.1418560545463037D-27,.1913375494454224D-30,.1192248760098222D-33,
     *.2671511219240137D-37,.1338616942106256D-41,.4510536193898974D-47/


      IER=0
      EXPA=EXP(-A)
!	The 2-point formula
      R1=(F(X(1)+A)*W(1)+F(X(2)+A)*W(2))*EXPA
      NPT=2
      N=2

!	Use higher order formula until convergence
      DO 2000 J=2,5
        N=N*2
        R2=0.0
        DO 1000 I=N-1,2*N-2
1000    R2=R2+F(X(I)+A)*W(I)
        R2=R2*EXPA

        NPT=NPT+N
        DIF=R2-R1
        RINT=R2
        IF(ABS(DIF).LT.MAX(AEPS,REPS*ABS(RINT))) RETURN
        R1=R2
2000  CONTINUE

!	Integral fails to converge
      IER=30
      RETURN
      END

!     ---------------------------------------------
 
      FUNCTION F(X)
      IMPLICIT REAL*8(A-H,O-Z)

!     THE INTEGRAND FOR SUBROUTINE GAUSS

      F=X/(1.+EXP(X))
      END

!     -------------------------------
 
      FUNCTION F2(X)
      IMPLICIT REAL*8(A-H,O-Z)

!     THE INTEGRAND*EXP(X) FOR SUBROUTINE LAGURE

      F2=X/(1+EXP(-X))
      END
