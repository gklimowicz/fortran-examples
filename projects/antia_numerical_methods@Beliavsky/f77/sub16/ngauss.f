!	Multiple integration over a hyper-rectangle in n-dimensions
!	using product Gauss-Legendre formulas with given number of points
!
!	A : (input) Real array of length N containing the lower limit
!		along each dimension
!	B : (input) Real array of length N containing the upper limit
!		along each dimension
!	N : (input) The number of dimensions, N>0 and N<NMAX (=21)
!	M : (input) Integer array of length N specifying the formula
!		to be used along each dimension. M(J)-point formula will
!		be used along Jth dimension, M(J) should be 2,4,8,16 or 32
!		otherwise IER is set to 306 and no calculations are done
!	IND : (input) Integer array of length N specifying the number
!		of subintervals to be used along each dimension. IND(J)>0
!		otherwise IER is set to 306 and no calculations are done
!	F : (input) Name of the function routine to calculate the integrand
!		FUNCTION(N,X) should calculate the integrand, where N is the
!		number of dimensions and X is a real array of length N containing
!		the coordinates of the point where integrand is to be calculated
!	RI : (output) The calculated value of the integral
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=305 implies N<1 or N.GE.NMAX, in which case no calculations are done
!		IER=306 implies M(J) is not 2,4,8,16 or 32 or IND(J)<1 for some J
!			in which case no calculations are done
!		IER=307 implies that number of points exceeded MAXPT and
!			no calculations are done
!	NUM : (output) Number of function evaluations used by subroutine
!	MAXPT : (input/output) Maximum number of function evaluations permitted
!		If MAXPT <1 it is set to a default value of MAXPTS (=1100000)
!	
!	Required routines : F
!
      SUBROUTINE NGAUSS(A,B,N,M,IND,F,RI,IER,NUM,MAXPT)
      IMPLICIT REAL*16(A-H,O-Z)
      PARAMETER(NMAX=21,MAXPTS=1100000)
      DIMENSION IP(NMAX),NPT(NMAX),H(NMAX),XA(NMAX),WT(NMAX)
      DIMENSION A(N),B(N),IND(N),M(N),W(31),X(31)

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

      IER=305
      RI=0.0
      NUM=0
      IF(N.GE.NMAX.OR.N.LT.1) RETURN
      IF(MAXPT.LE.0) MAXPT=MAXPTS
      IER=307

!	calculate the number of function evaluations required
      NUM=M(1)*IND(1)
      DO 200 I=2,N
200   NUM=NUM*M(I)*IND(I)
      IF(NUM.GT.MAXPT) RETURN

!	Initialisation
      IER=0
      DO 1000 I=1,N
        IP(I)=0
        NPT(I)=M(I)*IND(I)-1
        IF(M(I).NE.2.AND.M(I).NE.4.AND.M(I).NE.8.AND.M(I).NE.16.
     1      AND.M(I).NE.32) IER=306
        IF(IND(I).LT.1) IER=306
1000  H(I)=(B(I)-A(I))/(2*IND(I))
      IF(IER.NE.0) RETURN
      DO 1200 I=N+1,NMAX
        H(I)=1.0
        WT(I)=1.0
1200  CONTINUE

!	Loop for sum over N dimensions
      K=N

3000  DO 3100 I=K,1,-1
        M2=M(I)/2
!	The abscissas are X(NO),...,X(NO+M2-1)
        NO=M2
        H1=H(I)
        J1=IP(I)/M(I)
        J2=IP(I)-J1*M(I)
!	Use the (J2+1)th point in (J1+1)th subinterval
        X1=A(I)+(2*J1+1)*H1
        IF(J2.LT.M2) THEN
!	For the first M2 abscissas
          XA(I)=X1+H1*X(NO+J2)
          WT(I)=W(NO+J2)*WT(I+1)
        ELSE IF(J2-M2.LT.M2) THEN
!	For the next M2 abscissas
          XA(I)=X1-H1*X(NO+J2-M2)
          WT(I)=W(NO+J2-M2)*WT(I+1)
        ELSE
!	For Gaussian formula with odd number of points use the abscissa at x=0
          XA(I)=X1
          WT(I)=W(NO+M2)*WT(I+1)
        ENDIF
3100  CONTINUE

!	Add the new point to running sum
      RI=RI+WT(1)*F(N,XA)
      K=1
3200  IF(IP(K).GE.NPT(K)) GO TO 3400
!	try next point along Kth dimension
      IP(K)=IP(K)+1
      GO TO 3000

!	If Kth dimension is exhausted go to next one
3400  IP(K)=0
      K=K+1
      IF(K.LE.N) GO TO 3200

!	If all points are exhausted compute the value of integral
      DO 4000 I=1,N
4000  RI=RI*H(I)
      END
