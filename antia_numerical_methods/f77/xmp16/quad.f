!     PROGRAM FOR INTEGRATION OVER A FINITE INTERVAL
!     USING SIMPSON'S RULE OR ROMBRG INTEGRATION OR
!     EPSILON ALGORITHM OR GAUSS-LEGENDRE FORMULAS OR
!     ADAPTIVE INTEGRATION

      PROGRAM QUAD
      IMPLICIT REAL*16(A-H,O-Z)
      DIMENSION GI(20)
      EXTERNAL FUN

!     EXAMPLE 6.1: INTEGRATE SQRT(X) OVER [0,1]

51    FORMAT('   EXPONENTS IN ERROR EXPANSION:'/13F6.2)
52    FORMAT(I8,' POINT COMPOSITE GAUSS-LEGENDRE FORMULA')
53    FORMAT('   A =',1PD14.6,5X,'B =',D14.6,5X,'IER =',I4,5X,'IT=',I2/
     1       '    NO. OF FUNCTION EVALUATIONS =',I7,5X,'INTEGRAL =',
     2       D14.6/5X,'ESTIMATED ERROR =',2D14.6)

      REPS=1.Q-33
      AEPS=1.Q-34

100   PRINT *,'TYPE   A=LOWER LIMIT,  B=UPPER LIMIT,',
     1        '    (QUITS WHEN A.EQ.B)'
      READ *,A,B
      IF(A.EQ.B) STOP
      PRINT *,'TYPE IT=1/2/3/4/5  FOR SIMSON/ROMBRG/EPSILN/GAUSS/ADPINT'
      READ *,IT
      IF(IT.LT.1.OR.IT.GT.5) STOP
      REX=(B**1.5Q0-A**1.5Q0)/1.5Q0
      IF(IT.EQ.1) THEN
        CALL SIMSON(RI,A,B,REPS,AEPS,DIF,IER,NPT,FUN)
      ELSE IF(IT.EQ.2) THEN
        PRINT *,'TYPE (GI(I),I=1,13)  THE EXPONENTS IN ERROR EXPANSION'
        NPT=0
        READ *,(GI(I),I=1,13)
        CALL ROMBRG(RI,A,B,GI,REPS,AEPS,DIF,IER,NPT,FUN)
        WRITE(6,51) (GI(I),I=1,13)
      ELSE IF(IT.EQ.3) THEN
        NPT=0
        CALL EPSILN(RI,A,B,REPS,AEPS,DIF,IER,NPT,FUN)
      ELSE IF(IT.EQ.4) THEN
        PRINT *, 'TYPE NP=2/4/8/16/32   GAUSSIAN FORMULA TO BE USED'
        READ *,NP
        CALL GAUSS(RI,A,B,NP,REPS,AEPS,DIF,IER,NPT,FUN)
        WRITE(6,52) NP
      ELSE 
        NMAX=16385
        CALL ADPINT(RI,A,B,REPS,AEPS,DIF,FUN,IER,NPT,NMAX)
      ENDIF
      WRITE(6,53) A,B,IER,IT,NPT,RI,DIF
      GO TO 100
      END
 
!     ---------------------------------------------------------

!	To integrate a function over finite interval using adaptive control
!	of step size
!
!	RINT : (output) Calculated value of the integral
!	XL : (input) The lower limit
!	XU : (input) The upper limit
!	REPS : (input) The required relative accuracy
!	AEPS : (input) The required absolute accuracy
!		The estimated error should be less than MAX(AEPS,REPS*ABS(RINT))
!	DIF : (output) estimated (absolute) error achieved by the subroutine
!	F : (input) Name of the function routine to calculate the integrand
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=31 implies specified accuracy was not achieved on
!			at least one subinterval
!		IER=32 implies that this failure occurred more than IFMAX (=5) times
!		IER=325 implies that subroutine failed to attain required
!			accuracy using NMAX function evaluations
!		In all cases DIF will contain the estimated accuracy
!	NPT : (output) Number of function evaluations used by subroutine
!	NMAX : (input/output) Maximum number of function evaluations to be tried
!		If NMAX.LE.0 it is set to MAXPT (=100000)
!
!		FUNCTION F(X) must be supplied by the user.
!
!	Required routines : KRONRD (or GAUS16), F
!
!	The weights in KRONRD are accurate only to REAL*16 precision and hence if
!	higher precision is required it will be preferable to use GAUS16, although
!	it is less efficient.

      SUBROUTINE ADPINT(RINT,XL,XU,REPS,AEPS,DIF,F,IER,NPT,NMAX)
      IMPLICIT REAL*16(A-H,O,P,R-Z)
      IMPLICIT LOGICAL(Q)
      PARAMETER(IPMAX=200,IFMAX=5,MAXPT=100000)
      EXTERNAL F
      DIMENSION XU1(IPMAX)

      IER=0
      IFAIL=0
      RINT=0.0
      DIF=0.0
      IF(XL.EQ.XU) RETURN
      IF(NMAX.LE.0) NMAX=MAXPT
      AEPSL=AEPS
      NPT=0
      RL=XL
      RU=XU
      IU=0

!	To evaluate the integral over [RL,RU]
1000  CALL KRONRD(FINT,RL,RU,DIF0,NP,F)
!1000  CALL GAUS16(FINT,RL,RU,DIF0,NP,F)
      NPT=NPT+NP
      RM=0.5*(RL+RU)
!	Q=.TRUE. if the interval cannot be divided further
      Q=IU.GE.IPMAX.OR.RM.EQ.RL.OR.RM.EQ.RU

      IF(DIF0.LT.MAX(ABS(FINT)*REPS,AEPSL).OR.Q) THEN
!	Accept the value of FINT if adequate convergence or if the interval
!	cannot be subdivided further
        RINT=RINT+FINT
        DIF=DIF+DIF0
        IF(Q.AND.DIF0.GT.MAX(ABS(RINT)*REPS,AEPSL)) THEN
!	Integration fails to converge on this subinterval. Go to the next subinterval
          IER=31
          IFAIL=IFAIL+1
          IF(IFAIL.GT.IFMAX) THEN
!	If failure is frequent then adjust the convergence criterion.
            IER=32
            AEPSL=DIF*0.5
          ENDIF
        ENDIF

!	If all subintervals are exhausted then return
        IF(IU.LE.0) RETURN

!	otherwise try next subinterval
        RL=RU
        RU=XU1(IU)
        IU=IU-1
      ELSE

!	Subdivide the current interval and try again
        IU=IU+1
        XU1(IU)=RU
        RU=RM
      ENDIF

      IF(NPT.LT.NMAX) GO TO 1000
!	If the number of function evaluations has exceeded the limit then
!	try a last call to estimate the integral over the remaining interval
      IER=325
      RU=XU
      CALL KRONRD(FINT,RL,RU,DIF0,NP,F)
!      CALL GAUS16(FINT,RL,RU,DIF0,NP,F)
      NPT=NPT+NP
      RINT=RINT+FINT
      DIF=DIF+DIF0
      END

!     ----------------------------------------------------------

!	To integrate a function over finite interval using Epsilon algorithm
!
!	RI : (output) Calculated value of the integral
!	A : (input) The lower limit
!	B : (input) The upper limit
!	REPS : (input) The required relative accuracy
!	AEPS : (input) The required absolute accuracy
!		The estimated error should be less than MAX(AEPS,REPS*ABS(RI))
!	DIF : (output) estimated (absolute) error achieved by the subroutine
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=30 implies specified accuracy was not achieved
!			DIF will contain the estimated accuracy
!		IER=33 implies that N>NPT (=100) in which case it is set to 2
!		IER=34 implies that at some stage denominator vanished while
!			calculating epsilon table. This value is ignored.
!		IER=35 implies that roundoff error appears to be dominating
!	N : (input/output) On input it should contain the number of function
!		evaluations to be used for first estimate. If N<2 or N>NPT it
!		is set to 2. After execution it will contain the number of
!		function evaluations actually used by subroutine
!	FUN : (input) Name of the function routine to calculate the integrand
!		FUNCTION FUN(X) must be supplied by the user.
!
!	Required routines : FUN

      SUBROUTINE EPSILN(RI,A,B,REPS,AEPS,DIF,IER,N,FUN)
      IMPLICIT REAL*16(A-H,O-Z)
      PARAMETER(NMIN=4,NMAX=13,NPT=100)
      DIMENSION T(NMAX,NMAX)

      IER=0
      IF(N.LE.1) N=2
      IF(N.GT.NPT) THEN
        IER=33
        N=2
      ENDIF

!	Sum of the end points for trapezoidal rule
      S1=0.5*(FUN(A)+FUN(B))
      ND=1
      RI=0.0
      T(1,1)=0.0

      DO 4000 I=1,NMAX-1
        H=(B-A)/(N-1)

        DO 2200 J=2,N-1,ND
          Y=A+(J-1)*H
2200    S1=S1+FUN(Y)
!	The trapezoidal rule approximation
        T(I,2)=S1*H
        T(I+1,1)=0.0
        RI1=RI
        IF(I.GE.2) THEN
          DIF=ABS(T(I,2)-T(I-1,2))
          RI=T(I,2)
        ENDIF

!	Construct the Epsilon table
        DO 2400 J=3,I+1
          DEN=T(I-J+3,J-1)-T(I-J+2,J-1)

!	If denominator is zero set the error flag
          IF(DEN.NE.0.0) THEN
            T(I-J+2,J)=T(I-J+3,J-2)+1./DEN
          ELSE
            IER=34
            T(I-J+2,J)=T(I-J+3,J-2)
          ENDIF

2400    CONTINUE

        IF(I.GT.4) THEN
!	DIF is the minimum difference between two rows of epsilon table
          DO 2600 J=4,I-1,2
            DIF1=ABS(T(I-J+2,J)-T(I-J+1,J))
            IF(DIF1.LT.DIF) THEN
              DIF=DIF1
              RI=T(I-J+2,J)
            ENDIF
2600      CONTINUE
        ENDIF

        ND=2
        IF(I.LE.NMIN) GO TO 4000
        IF(I.GT.6.AND.DIF.GT.DIF0) THEN
!	Roundoff error appears to be dominating, retain the previous value of RI
          IER=35
          RI=RI1
          RETURN
        ENDIF
        DIF0=DIF
        IF(DIF.LT.MAX(REPS*ABS(RI),AEPS)) RETURN

4000  N=2*N-1

!	Integral fails to converge
      IER=30
      N=(N+1)/2
      END

!     ----------------------------------------------------------

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
      IMPLICIT REAL*16(A-H,O-Z)
      PARAMETER(NMAX=9)
      DIMENSION W(31),X(31)

!	Weights and abscissas for Gauss-Legendre quadrature.
!	For N-point formula W(K)=W(N-K+1) and X(K)=-X(N-K+1)
!		For K=1,2,...,N/2. Hence only half points are tabulated.
!	For 2-point W(1); 4-point W(2), W(3); 8-point W(4),...,W(7);
!	16-point W(8),...,W(15); 32-point W(16),...,W(31) are the
!	weights corresponding to abscissas X(I).

      DATA W  /  1.0D+00,
     *           3.4785484513745385737306394922199943Q-1,
     *           6.5214515486254614262693605077800183Q-1,
     *           1.0122853629037625915253135430996220Q-1,
     *           2.2238103445337447054435599442624126Q-1,
     *           3.1370664587788728733796220198660110Q-1,
     *           3.6268378337836198296515044927719567Q-1,
     *           2.7152459411754094851780572456019197Q-2,
     *           6.2253523938647892862843836994377013Q-2,
     *           9.5158511682492784809925107602247148Q-2,
     *           1.2462897125553387205247628219201544Q-1,
     *           1.4959598881657673208150173054747880Q-1,
     *           1.6915651939500253818931207903036093Q-1,
     *           1.8260341504492358886676366796921865Q-1,
     *           1.8945061045506849628539672320828440Q-1,
     *           7.0186100094700966004070637388536822Q-3,
     *           1.6274394730905670605170562206387365Q-2,
     *           2.5392065309262059455752589789223263Q-2,
     *           3.4273862913021433102687732252373180Q-2,
     *           4.2835898022226680656878646606125681Q-2,
     *           5.0998059262376176196163244689520896Q-2,
     *           5.8684093478535547145283637300170933Q-2,
     *           6.5822222776361846837650063706937930Q-2,
     *           7.2345794108848506225399356478487401Q-2,
     *           7.8193895787070306471740918828308613Q-2,
     *           8.3311924226946755222199074604350642Q-2,
     *           8.7652093004403811142771462751801258Q-2,
     *           9.1173878695763884712868577111635117Q-2,
     *           9.3844399080804565639180237668115127Q-2,
     *           9.5638720079274859419082002204134859Q-2,
     *           9.6540088514727800566764830063574027Q-2/

      DATA X  /  5.7735026918962576450914878050195760Q-1,
     *           8.6113631159405257522394648889280965Q-1,
     *           3.3998104358485626480266575910324466Q-1,
     *           9.6028985649753623168356086856947230Q-1,
     *           7.9666647741362673959155393647583006Q-1,
     *           5.2553240991632898581773904918924604Q-1,
     *           1.8343464249564980493947614236018303Q-1,
     *           9.8940093499164993259615417345033306Q-1,
     *           9.4457502307323257607798841553460820Q-1,
     *           8.6563120238783174388046789771239397Q-1,
     *           7.5540440835500303389510119484744268Q-1,
     *           6.1787624440264374844667176404879140Q-1,
     *           4.5801677765722738634241944298357810Q-1,
     *           2.8160355077925891323046050146049710Q-1,
     *           9.5012509837637440185319335424958114Q-2,
     *           9.9726386184948156354498112866504099Q-1,
     *           9.8561151154526833540017504463090231Q-1,
     *           9.6476225558750643077381192811827550Q-1,
     *           9.3490607593773968917091913483540993Q-1,
     *           8.9632115576605212396530724371921239Q-1,
     *           8.4936761373256997013369300496774309Q-1,
     *           7.9448379596794240696309729897042952Q-1,
     *           7.3218211874028968038742666509126756Q-1,
     *           6.6304426693021520097511516866323809Q-1,
     *           5.8771575724076232904074547640182703Q-1,
     *           5.0689990893222939002374747437782170Q-1,
     *           4.2135127613063534536411943617242740Q-1,
     *           3.3186860228212764977991680573018860Q-1,
     *           2.3928736225213707454460320916550261Q-1,
     *           1.4447196158279649348518637359881043Q-1,
     *           4.8307665687738316234812570440502563Q-2/


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

!     ----------------------------------------------------------

!	To integrate a function over a finite interval using Gauss-Kronrod formula
!	For use with ADPINT
!	Since the weights and abscissas are not accurate to REAL*16
!	accuracy, this routine will not achive the maximum accuracy
!	permissible by arithmetic. Use GAUS16 instead of KRONRD
!	if high accuracy is required.
!
!	RI : (output) Calculated value of the integral
!	A : (input) The lower limit
!	B : (input) The upper limit
!	DIF : (output) estimated (absolute) error achieved by the subroutine
!	N : (output) Number of function evaluations used by subroutine
!	F : (input) Name of the function routine to calculate the integrand
!
!	FUNCTION F(X) must be supplied by the user
!
!	Required routines : F

      SUBROUTINE KRONRD(RI,A,B,DIF,N,F)
      IMPLICIT REAL*16(A-H,O-Z)
      DIMENSION  W7(4),A7(4),WK7(4),WK15(4),AK15(4)

!	W7 and A7 are the weights and abscissas for the 7-point Gauss formula
!	WK7 are the weights for these points in Kronrod formula
!	WK15 and AK15 are the weights and abscissas for the remaining points
!	in Kronrod formula.
!	Because of symmetry only half the points are given.

      DATA W7  /0.1294849661688696932706114326790820183Q0,
     *          0.2797053914892766679014677714237795825Q0,
     *          0.3818300505051189449503697754889751339Q0,
     *          0.4179591836734693877551020408163265306Q0/
      DATA A7  /0.9491079123427585245261896840478512624Q0,
     *          0.7415311855993944398638647732807884070Q0,
     *          0.4058451513773971669066064120769614633Q0, 0.0/
      DATA WK7 /0.0630920926299785532907006631892042866Q0,
     *          0.1406532597155259187451895905102379204Q0,
     *          0.1903505780647854099132564024210136828Q0,
     *          0.2094821410847278280129991748917142637Q0/
      DATA WK15/0.0229353220105292249637320080589695920Q0,
     *          0.1047900103222501838398763225415180174Q0,
     *          0.1690047266392679028265834265985502841Q0,
     *          0.2044329400752988924141619992346490847Q0/
      DATA AK15/0.9914553711208126392068546975263285166Q0,
     *          0.8648644233597690727897127886409262012Q0,
     *          0.5860872354676911302941448382587295984Q0,
     *          0.2077849550078984676006894037732449134Q0/

      AT=(B-A)/2.
      BT=(B+A)/2.
      FBT=F(BT)
      R1=W7(4)*FBT
      RI=WK7(4)*FBT
      DO 2000 K=1,3
        F1=F(AT*A7(K)+BT)
        F2=F(BT-AT*A7(K))
!	7-point Gauss-Legendre formula
        R1=R1+W7(K)*(F1+F2)
!	15-point Kronrod formula
        RI=RI+WK7(K)*(F1+F2)
2000  CONTINUE

      DO 2500 K=1,4
2500  RI=RI+WK15(K)*(F(AT*AK15(K)+BT)+F(BT-AT*AK15(K)))

      RI=RI*AT
      R1=R1*AT
      DIF=ABS(RI-R1)
      N=15
      END

!     ----------------------------------------------------------

!	To integrate a function over finite interval using Romberg integration
!
!	RI : (output) Calculated value of the integral
!	A : (input) The lower limit
!	B : (input) The upper limit
!	GI : (input/output) Real array of length NMAX (=13), containing
!		the expected values of exponents \gamma_i in error expansion
!		If GI(I).LE.0 it will be set to 2I, the correct value for
!		a smooth function
!	REPS : (input) The required relative accuracy
!	AEPS : (input) The required absolute accuracy
!		The estimated error should be less than MAX(AEPS,REPS*ABS(RI))
!	DIF : (output) estimated (absolute) error achieved by the subroutine
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=30 implies specified accuracy was not achieved
!			DIF will contain the estimated accuracy
!		IER=33 implies that N>NPT (=100) in which case it is set to 2
!	N : (input/output) On input it should contain the number of function evaluations
!		to be used for first estimate. If N<2 or N>NPT it is set to 2.
!		After execution it will contain the number of function
!		evaluations actually used by subroutine
!	FUN : (input) Name of the function routine to calculate the integrand
!		FUNCTION FUN(X) must be supplied by the user.
!
!	Required routines : FUN

      SUBROUTINE ROMBRG(RI,A,B,GI,REPS,AEPS,DIF,IER,N,FUN)
      IMPLICIT REAL*16(A-H,O-Z)
      PARAMETER(NMIN=3,NMAX=13,NPT=100)
      DIMENSION GI(NMAX),T(NMAX,NMAX)

      DO 1000 I=1,NMAX
        IF(GI(I).LE.0.0) GI(I)=2*I
1000  CONTINUE

      IER=0
      IF(N.LE.1) N=2
      IF(N.GT.NPT) THEN
        IER=33
        N=2
      ENDIF
!	Contribution from the end points
      S1=0.5*(FUN(A)+FUN(B))
!	First time use all points
      ND=1
      DIF=0.0

      DO 4000 I=1,NMAX
        H=(B-A)/(N-1)

!	Add new points to the sum
        DO 2200 J=2,N-1,ND
          Y=A+(J-1)*H
2200    S1=S1+FUN(Y)
!	The trapezoidal rule approximation
        T(I,1)=S1*H

!	The Richardson's extrapolation
        DO 2400 J=1,I-1
          FJ=2.**GI(J)
          T(I,J+1)=T(I,J)+(T(I,J)-T(I-1,J))/(FJ-1)
          DIF1=ABS(T(I,J)-T(I-1,J))
!	Find the minimum difference between the last two rows of T-table
          IF(DIF1.LT.DIF.OR.J.EQ.1) THEN
            DIF=DIF1
            RI=T(I,J)
          ENDIF
2400    CONTINUE

!	On second and subsequent pass add only new points to the sum
        ND=2
        IF(I.LE.NMIN) GO TO 4000
        IF(DIF.LT.MAX(REPS*ABS(RI),AEPS)) RETURN

4000  N=2*N-1

!	Routine fails to converge
      IER=30
      N=(N+1)/2
      END

!     ----------------------------------------------------------

!	To integrate a function over finite interval using Simpson's rule
!
!	RI : (output) Calculated value of the integral
!	XL : (input) The lower limit
!	XU : (input) The upper limit
!	REPS : (input) The required relative accuracy
!	AEPS : (input) The required absolute accuracy
!		The estimated error should be less than MAX(AEPS,REPS*ABS(RI))
!	DIF : (output) estimated (absolute) error achieved by the subroutine
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=30 implies specified accuracy was not achieved
!			DIF will contain the estimated accuracy
!	N : (output) Number of function evaluations used by subroutine
!	FUN : (input) Name of the function routine to calculate the integrand
!		FUNCTION FUN(X) must be supplied by the user.
!
!	Required routines : FUN

      SUBROUTINE SIMSON(RI,XL,XU,REPS,AEPS,DIF,IER,N,FUN)
      IMPLICIT REAL*16(A-H,O-Z)
      PARAMETER(NMIN=3,NMAX=13)

      FEND=FUN(XL)+FUN(XU)
      EVEN=0.0
      ODD=0.0
      IER=0
      RI=0.0
      DIF=0.0

      N=2
!	starting with 2+1 points, subdivide the intervals into 2 until convergence
      H=(XU-XL)
      IF(H.EQ.0.0) RETURN

      DO 3000 I=1,NMAX
        H=H/2.
        EVEN=EVEN+ODD
        ODD=0.0
        X1=XL+H
        N2=N/2
        H2=2.*H

        DO 1000 J=1,N2
          X=X1+H2*(J-1)
          ODD=ODD+FUN(X)
1000    CONTINUE
!	Estimate for the integral
        R1=(FEND+4.*ODD+2.*EVEN)*H/3.

        DIF=R1-RI
        RI=R1
!	To avoid spurious convergence in first few trials skip the convergence test
        IF(I.LE.NMIN) GO TO 3000
        IF(ABS(DIF).LT.MAX(REPS*ABS(R1),AEPS)) RETURN
3000  N=N*2

      N=N/2
      IER=30
      END

!     ----------------------------------------------------------
 
      FUNCTION FUN(X)
      IMPLICIT REAL*16(A-H,O-Z)

!     THE INTEGRAND

      FUN=SQRT(X)
      END
