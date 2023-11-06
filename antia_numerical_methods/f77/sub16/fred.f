!	To solve linear Fredholm equation using quadrature method
!
!	M : (input) Number of abscissas to be used in quadrature formula
!	A : (input) Lower limit of the integral
!	B : (input) Upper limit of the integral
!	WT : (input/output) Real array of length M containing the
!		weights used in quadrature formula.
!		If IQ is negative the weights must be supplied by the user
!		otherwise they are calculated by the routine
!	X : (input/output) Real array of length M containing the
!		abscissas used in quadrature formula.
!		If IQ is negative the abscissas must be supplied by the user
!		otherwise they are calculated by the routine
!	F : (output) Real array of length M containing the calculated solution
!		F(I) is the value of solution at X(I).
!	FC : (output) Real array of length M containing the calculated solution
!		after applying the deferred correction. This will be
!		relevant only if IQ=1 and IT=1,2. In other cases a dummy
!		array of any length may be supplied.
!	FG : (input) Name of the function routine used to calculate the right
!		hand side g(X). For IT=3 g(x) is not required, but in that
!		case this function is used to calculate an initial guess
!		for the eigenfunctions. In most cases the inverse iteration
!		converges from essentially arbitrary initial guess and it
!		is enough to set FG(X) to any nonzero value.
!	FKER : (input) Name of the function routine used to calculate the
!		kernel K(x,t)
!	EI : (input/output) Initial guess for the eigenvalue. After execution
!		it will contain the computed eigenvalue. Used only for IT=3
!	WK : Real array used as scratch space. Its length should be M*M for
!		IT=1,2, while for IT=3 it should be 2M*M+M.
!	IWK : Integer array of length M used as a scratch space
!	IQ : (input) Integer variable to specify the quadrature formula to be used.
!		If IQ=1, then trapezoidal rule is used and deferred correction
!			is calculated using Gregory's formula
!		If IQ=2, the Simpson's 1/3 rule is used
!		If IQ=4,8,16,32 a composite rule using IQ point
!			Gauss-Legendre formula is used.
!		In all these cases the weights and abscissas are calculated
!		If IQ is negative then it is assumed that weights and abscissas
!		are supplied in arrays WT and X.
!		Other values of IQ will cause an error return.
!	IT : (input) Integer variable to specify the type of integral equation
!		If IT=1 Fredholm equation of the first kind is solved
!		If IT=2 Fredholm equation of the second kind is solved
!		If IT=3 Fredholm equation of the third kind (eigenvalue
!			problem) is solved
!	REPS : (input) Required relative accuracy in calculating eigenvalue
!		and eigenvectors. It is not used for IT=1,2. It only specifies
!		the accuracy to which inverse iteration converges. It does not
!		control the truncation error.
!	IER : (output) The error parameter, IER=0 implies successful execution
!		IER=-11 implies that the calculations are actually performed
!			using a smaller number of abscissas than M
!		IER=706 implies M<3, IT>3, IT.LE.0 or M is not sufficient
!			to apply the required quadrature formula. No calculations
!			are done.
!		IER=707 implies that IQ is not acceptable and no calculations
!			are done
!		Other values of IER may be set by GAUELM and INVIT
!
!	FUNCTION FG(X) and FUNCTION FKER(X,T) must be supplied by the user.
!		FG is the right hand side function g(x) and FKER(X,T) is 
!		the kernel. The integral equation is specified in the form
!		given by Eq.(12.1)
!
!	Required routines : GAUELM, INVIT, FG, FKER
!	
      SUBROUTINE FRED(M,A,B,WT,X,F,FC,FG,FKER,EI,WK,IWK,IQ,IT,REPS,IER)
      IMPLICIT REAL*16(A-H,O-Z)
      DIMENSION WT(M),F(M),FC(*),X(M),WK(M,*),IWK(M)
      DIMENSION W(31),XA(31)

!	Weights and abscissas for Gauss-Legendre quadrature.
!	For N-point formula W(K)=W(N-K+1) and XA(K)=-XA(N-K+1)
!		For K=1,2,...,N/2. Hence only half points are tabulated.
!	For 2-point W(1); 4-point W(2), W(3); 8-point W(4),...,W(7);
!	16-point W(8),...,W(15); 32-point W(16),...,W(31) are the
!	weights corresponding to abscissas XA(I).

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

      DATA XA /  5.7735026918962576450914878050195760Q-1,
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

      IER=0
      IF(M.LT.3.OR.IT.GT.3.OR.IT.LE.0) THEN
        IER=706
        RETURN
      ENDIF

!	M should not be changed since it is used in dimension statement
      N=M
      IF(IQ.EQ.1) THEN
!	Use the trapezoidal rule
        H=(B-A)/(N-1)
        WT(1)=0.5*H
        WT(N)=WT(1)
        X(1)=A
        X(N)=B
        DO 2000 I=2,N-1
          X(I)=A+(I-1)*H
2000    WT(I)=H

      ELSE IF(IQ.EQ.2) THEN
!	Use the Simpson's 1/3 rule, if N is even, then reduce it by 1
        N=2*((N-1)/2)+1
        H=(B-A)/(N-1)
        WT(1)=H/3.
        X(1)=A
        DO 2100 I=2,N-1,2
          X(I)=A+(I-1)*H
          X(I+1)=A+I*H
          WT(I)=4.*H/3.
2100    WT(I+1)=2.*H/3.
        WT(N)=WT(1)
        X(N)=B

      ELSE IF(IQ.GE.0) THEN
!	Try Gauss-Legendre formulas
        NO=-1
        IF(IQ.EQ.4) NO=1
        IF(IQ.EQ.8) NO=3
        IF(IQ.EQ.16) NO=7
        IF(IQ.EQ.32) NO=15
        IF(NO.LT.0) THEN
!	If IQ is not acceptable then quit
          IER=707
          RETURN
        ENDIF
        N=(N/IQ)*IQ
        IF(N.LT.IQ) THEN
!	If the number of points is not sufficient, then quit
          IER=706
          RETURN
        ENDIF

!	Setup the weights and abscissas for Gauss-Legendre formula
        H=IQ*(B-A)/N
        DO 2300 I=1,N,IQ
          A1=A+(I-1)*H/IQ+H/2.
          DO 2300 I1=1,IQ/2
            WT(I+I1-1)=W(NO+I1)*H/2.
            WT(I+IQ-I1)=WT(I+I1-1)
            X(I+I1-1)=A1-XA(NO+I1)*H/2.
            X(I+IQ-I1)=A1+XA(NO+I1)*H/2.
2300    CONTINUE
      ENDIF
      IF(M.NE.N) IER=-11

!	Setting up the equation matrix and the right hand side
      DO 3000 I=1,N
        F(I)=FG(X(I))
        DO 2800 J=1,N
          WK(J,I)=WT(I)*FKER(X(J),X(I))
2800    CONTINUE
        IF(IT.EQ.2) WK(I,I)=WK(I,I)-1.
3000  CONTINUE

      NUM=1
      LJ=M
      IF(IT.LE.2) THEN
!	Solve the system of linear equations
        IFLG=0
        CALL GAUELM(N,NUM,WK,F,DET,IWK,LJ,IER,IFLG)
        IF(IER.GT.0) RETURN
      ELSE
!	Solve the eigenvalue problem
        IFLG=0
        P=EI
        NIT=0
        CALL INVIT(WK,N,LJ,P,F,IFLG,EI,RC,REPS,WK(1,N+1),IWK,NIT,IER)
      ENDIF

      IF(IQ.EQ.1.AND.IT.NE.3) THEN
!	Apply the deferred correction
        DO 3200 I=1,N
          U2=FKER(X(I),X(2))*F(2)
          U2M=FKER(X(I),X(N-1))*F(N-1)
          T1=FKER(X(I),X(1))*F(1)-U2+FKER(X(I),X(N))*F(N)-U2M
          T2=T1+FKER(X(I),X(3))*F(3)-U2+FKER(X(I),X(N-2))*F(N-2)-U2M
          FC(I)=H*T1/12.+H*T2/24.
3200    CONTINUE
        CALL GAUELM(N,NUM,WK,FC,DET,IWK,LJ,IER,IFLG)

        DO 4000 I=1,N
4000    FC(I)=F(I)+FC(I)
      ENDIF
      END
