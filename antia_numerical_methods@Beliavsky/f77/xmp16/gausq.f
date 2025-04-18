!     PROGRAM FOR INTEGRATING F(X)/SQRT(X) OVER (0,A] 

      PROGRAM INTEG
      IMPLICIT REAL*16(A-H,O-Z)
      EXTERNAL F,F2

!     EXERCISE 6.14

51    FORMAT('   IER =',I4,5X,'A =',1PD14.6,5X,'A1 =',D14.6,5X,
     1       'NO. OF FUNCTION EVALUATIONS =',I7/5X,'INTEGRAL =',
     2       D14.6,5X,'ESTIMATED ERROR =',3D14.6)

      REPS=1.Q-33
      AEPS=1.Q-39

100   PRINT *,'TYPE  A=UPPER LIMIT,  A1=POINT AT WHICH INTEGRAL IS TO'
     1       , ' BE SPLIT'
      PRINT *,'                     (QUITS WHEN A<0)'
      READ *,A,A1
      IF(A.LT.0.0) STOP
      CALL GAUSQ2(RI,A,A1,REPS,AEPS,DIF,F,F2,N,IER)
      WRITE(6,51) IER,A,A1,N,RI,DIF
      GO TO 100
      END
 
!     ---------------------------------------------------
 
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
 
!     ---------------------------------------------

!	To integrate a function over (0,A] using a combination of
!	Gaussian formulas
!
!	RINT : (output) Calculated value of the integral
!	A : (input) The upper limit
!	A1 : (input/output) The point at which integral has to be broken
!		A1 will be adjusted by the subroutine.
!	REPS : (input) The required relative accuracy
!	AEPS : (input) The required absolute accuracy
!		The estimated error should be less than MAX(AEPS,REPS*ABS(RINT))
!	DIF : (output) estimated (absolute) error achieved by the subroutine
!	F : (input) Name of the function routine to calculate the integrand
!	F2 : (input) Name of the function routine to calculate F(X)*SQRT(X)
!	NP : (output) Number of function evaluations used by the subroutine
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=30 implies specified accuracy was not achieved by GAUSS
!		IER=37 implies specified accuracy was not achieved by GAUSQ
!		IER=38 implies specified accuracy was not achieved by both
!		GAUSS and GAUSQ
!		In all cases DIF will contain the estimated accuracy
!
!	FUNCTION F(X) and F2(X) must be supplied by the user.
!
!	Required routines : GAUSS, GAUSQ, F, F2
 
      SUBROUTINE GAUSQ2(RINT,A,A1,REPS,AEPS,DIF,F,F2,NP,IER)
      IMPLICIT REAL*16(A-H,O-Z)
      PARAMETER(AMN=1.Q-2)
      EXTERNAL F,F2
 
      IER1=0
      R1=0.0
      R2=0.0
      DIF=0.0
      NP=0
      NPT=0
      IF(A1.GT.A) A1=A
      IF(A1.LE.0.0) A1=A
 
!	Evaluate the integral over (0,A1]
2200  CALL GAUSQ(R1,A1,AEPS,REPS,DIF1,IER,NPT1,F2)
      NP=NP+NPT1
      IF(IER.EQ.0) GO TO 2500
!	If GAUSQ fails decrease A1
      T1=A1
      A1=A1/2.
      IF(A1.GT.AMN) GO TO 2200
!	To prevent infinite loop do not reduce A1 below AMN
      IER1=37
      IER=0
      A1=T1
 
!	Evaluate the integral over [A1,A]
2500  IF(A-A1.GT.AEPS) CALL GAUSS(R2,A1,A,16,REPS,AEPS,DIF,IER,NPT,F)
      IER=IER+IER1
      IF(IER.GT.IER1.AND.IER1.GT.0) IER=38
      RINT=R1+R2
      DIF=ABS(DIF)+ABS(DIF1)
      NP=NP+NPT
      END

!     --------------------------------------------------------

!	To integrate F(X)/SQRT(X) over (0,A] using Gaussian formulas with
!	1/SQRT(X) weight function
!
!	RINT : (output) Calculated value of the integral
!	A : (input) The upper limit
!	AEPS : (input) The required absolute accuracy
!	REPS : (input) The required relative accuracy
!		The estimated error should be less than MAX(AEPS,REPS*ABS(RINT))
!	DIF : (output) estimated (absolute) error achieved by the subroutine
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=30 implies specified accuracy was not achieved
!			DIF will contain the estimated accuracy
!	NPT : (output) Number of function evaluations used by the subroutine
!	FUN : (input) Name of the function routine to calculate the
!		integrand (multiplied by SQRT(X))
!
!	Function FUN(X) must be supplied by the user
!
!	Required routines : FUN

      SUBROUTINE GAUSQ(RINT,A,AEPS,REPS,DIF,IER,NPT,FUN)
      IMPLICIT REAL*16(A-H,O-Z)
      DIMENSION W(31),X(31)
 
!     WEIGHTS AND ABSCISSAS FOR GAUSS-LEGENDRE QUADRATURE
!     FOR N-POINT FORMULA  W(K) = W(N-K+1) AND X(K) = -X(N-K+1)
!     HENCE ONLY THE VALUES FOR K = 1 ... (N+1)/2 ARE TABULATED
!     FOR 2-POINTS (W(1))
!     4-POINTS (W(2),W(3)) , 8-POINTS (W(4) ... W(7)),
!     16-POINTS (W(8) ... W(15)) , 32-POINTS (W(16) ... W(31))
!     ARE THE WEIGHTS CORRESPONDING TO ABSCISSAS X(I)
!     FOR WEIGHT FUNCTION OF 1/SQRT(X) THE WEIGHTS ARE 2*W(I)
!     AND ABSCISSAS ARE X(I)**2
 
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
      IER=0
!	The 1-point formula
      R1=2.Q0*W(1)*FUN(A*X(1)**2)
      R1=R1*SQRT(A)
      NPT=0
 
!	Try higher order formulas until convergence 
      DO 3000 I=2,5
        N=N*2
 
        S1=0.0
        DO 2000 K=N,2*N-1
2000    S1=S1+2.Q0*W(K)*FUN(A*(X(K))**2)
        S1=S1*SQRT(A)
        NPT=NPT+N
        DIF=S1-R1
        RINT=S1
        IF(ABS(DIF).LT.MAX(AEPS,REPS*ABS(RINT))) RETURN
        R1=S1
3000  CONTINUE
 
!	Integral fails to converge
      IER=30
      END

!     ---------------------------------------------
 
      FUNCTION F(X)
      IMPLICIT REAL*16(A-H,O-Z)

!     THE INTEGRAND FOR SUBROUTINE GAUSS

      F=SQRT(X)*COS(X)
      END

!     -------------------------------
 
      FUNCTION F2(X)
      IMPLICIT REAL*16(A-H,O-Z)

!     THE INTEGRAND*SQRT(X) FOR USE WITH SUBROUTINE GAUSQ

      F2=X*COS(X)
      END

