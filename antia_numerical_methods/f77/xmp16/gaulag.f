!     PROGRAM FOR INTEGRATION OVER INFINITE REGION USING GAUSSIAN FORMULAS

      PROGRAM INTEG
      IMPLICIT REAL*16(A-H,O-Z)
      EXTERNAL F,F2

!     EXAMPLE 6.10

51    FORMAT('   IER =',I4,5X,'A =',1PD14.6,5X,'A1 =',D14.6/
     1       5X,'NO. OF FUNCTION',
     2       ' EVALUATIONS =',I7/5X,'INTEGRAL =',D14.6,5X,
     3       'ESTIMATED ERROR =',D14.6)
52    FORMAT('   INITIAL VALUE OF A1 =',1PD14.6)

      REPS=1.Q-33
      AEPS=1.Q-36

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
      IMPLICIT REAL*16(A-H,O-Z)
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
      IMPLICIT REAL*16(A-H,O-Z)
      DIMENSION W(62),X(62)

!	Weights and abscissas for Gauss-Laguerre quadrature
!	W(N-1),...,W(2N-2), are the weights for N-point rule and
!	X(N-1),...,X(2N-2), the corresponding abscissas
!	Weights and abscissas are available for N=2,4,8,16,32


      DATA X  /  5.8578643762690495119831127579030182Q-01,
     *           3.4142135623730950488016887242096982Q+00,
     *           3.2254768961939231180036145910436758Q-01,
     *           1.7457611011583465756868167125179473Q+00,
     *           4.5366202969211279832792853849571364Q+00,
     *           9.3950709123011331292335364434205469Q+00,
     *           1.7027963230510099978886185660829757Q-01,
     *           9.0370177679937991218602022355509030Q-01,
     *           2.2510866298661306893071183669686348Q+00,
     *           4.2667001702876587936494218269006407Q+00,
     *           7.0459054023934656972793254821193568Q+00,
     *           1.0758516010180995224059956788032032Q+01,
     *           1.5740678641278004578028761158402829Q+01,
     *           2.2863131736889264105700534297413110Q+01,
     *           8.7649410478927840360198097340102636Q-02,
     *           4.6269632891508083188083826066425360Q-01,
     *           1.1410577748312268568779450181051781Q+00,
     *           2.1292836450983806163261590706581807Q+00,
     *           3.4370866338932066452351070167518546Q+00,
     *           5.0780186145497679129230583081419303Q+00,
     *           7.0703385350482341303959894708013374Q+00,
     *           9.4383143363919387839472467291128199Q+00,
     *           1.2214223368866158736939124608824205Q+01,
     *           1.5441527368781617076764774162213519Q+01,
     *           1.9180156856753134854663140949739646Q+01,
     *           2.3515905693991908531823187275156244Q+01,
     *           2.8578729742882140367520613709877714Q+01,
     *           3.4583398702286625814527687177752595Q+01,
     *           4.1940452647688332635472233025186606Q+01,
     *           5.1701160339543318364342697119673885Q+01,
     *           4.4489365833267018418850194524358936Q-02,
     *           2.3452610951961853745290956130193197Q-01,
     *           5.7688462930188642649155256937756379Q-01,
     *           1.0724487538178176330409139771776169Q+00,
     *           1.7224087764446454411309329279747156Q+00,
     *           2.5283367064257948811241999055589149Q+00,
     *           3.4922132730219944896088033997715356Q+00,
     *           4.6164567697497673877620522961721438Q+00,
     *           5.9039585041742439465615214915812638Q+00,
     *           7.3581267331862411132219897371949988Q+00,
     *           8.9829409242125961033782475267656929Q+00,
     *           1.0783018632539972067501949138058345Q+01,
     *           1.2763697986742725114969033082205887Q+01,
     *           1.4931139755522557319796964687308652Q+01,
     *           1.7292454336715314789235718383619840Q+01,
     *           1.9855860940336054739789944584146077Q+01,
     *           2.2630889013196774488675779339387325Q+01,
     *           2.5628636022459247767476176876849794Q+01,
     *           2.8862101816323474744342640711515499Q+01,
     *           3.2346629153964737003232165423674708Q+01,
     *           3.6100494805751973804017118947901278Q+01,
     *           4.0145719771539441536209343928925638Q+01,
     *           4.4509207995754937975906604377521372Q+01,
     *           4.9224394987308639176722221806566310Q+01,
     *           5.4333721333396907332867181551169358Q+01,
     *           5.9892509162134018196130475324669364Q+01,
     *           6.5975377287935052796563076119348581Q+01,
     *           7.2687628090662708638675349087774606Q+01,
     *           8.0187446977913523067491638568690214Q+01,
     *           8.8735340417892398689355449524296604Q+01,
     *           9.8829542868283972559184478409520473Q+01,
     *           1.1175139809793769521366471653944923Q+02/


      DATA W  /  8.5355339059327376220042218105242454Q-01,
     *           1.4644660940672623779957781894757548Q-01,
     *           6.0315410434163360163596602381807856Q-01,
     *           3.5741869243779968664149201745809158Q-01,
     *           3.8887908515005384272438168156209922Q-02,
     *           5.3929470556132745010379056762059309Q-04,
     *           3.6918858934163752992058283937570413Q-01,
     *           4.1878678081434295607697858133333403Q-01,
     *           1.7579498663717180569965986677677358Q-01,
     *           3.3343492261215651522132534934406461Q-02,
     *           2.7945362352256725249389241479286418Q-03,
     *           9.0765087733582131042385014933573514Q-05,
     *           8.4857467162725315448680183089320506Q-07,
     *           1.0480011748715103816150885355156921Q-09,
     *           2.0615171495780099433427363674130636Q-01,
     *           3.3105785495088416599298309870968019Q-01,
     *           2.6579577764421415259950202065032923Q-01,
     *           1.3629693429637753997554751352639752Q-01,
     *           4.7328928694125218978062339278099457Q-02,
     *           1.1299900080339453231249045970117299Q-02,
     *           1.8490709435263108642917678325206456Q-03,
     *           2.0427191530827846012601833842122623Q-04,
     *           1.4844586873981298771351506755069150Q-05,
     *           6.8283193308711995643955959032657740Q-07,
     *           1.8810248410796732138815992041799925Q-08,
     *           2.8623502429738816196306262915599456Q-10,
     *           2.1270790332241029673903361097755404Q-12,
     *           6.2979670025178677871744621455152583Q-15,
     *           5.0504737000355128204021323330302525Q-18,
     *           4.1614623703728551904264835611619532Q-22,
     *           1.0921834195238497113613133794342890Q-01,
     *           2.1044310793881323293606207139170982Q-01,
     *           2.3521322966984800539494110667639529Q-01,
     *           1.9590333597288104341324790118177027Q-01,
     *           1.2998378628607176060721682232593595Q-01,
     *           7.0578623865717441560164332043306989Q-02,
     *           3.1760912509175070305825521162864987Q-02,
     *           1.1918214834838557056544650575558017Q-02,
     *           3.7388162946115247896612284779588774Q-03,
     *           9.8080330661495513223063061130806851Q-04,
     *           2.1486491880136418802319948368572206Q-04,
     *           3.9203419679879472043269568278231443Q-05,
     *           5.9345416128686328783558289377302396Q-06,
     *           7.4164045786675522190708220213025518Q-07,
     *           7.6045678791207814811192654594345483Q-08,
     *           6.3506022266258067424277710855230020Q-09,
     *           4.2813829710409288788136058209794941Q-10,
     *           2.3058994918913360792733680961757377Q-11,
     *           9.7993792887270940633345522599474347Q-13,
     *           3.2378016577292664623104264614174741Q-14,
     *           8.1718234434207194332018605917701098Q-16,
     *           1.5421338333938233721785594912914511Q-17,
     *           2.1197922901636186120409347437315661Q-19,
     *           2.0544296737880454266557098760188216Q-21,
     *           1.3469825866373951558051934047774799Q-23,
     *           5.6612941303973593711263443238216036Q-26,
     *           1.4185605454630369059514293389207924Q-28,
     *           1.9133754944542243093712782968307795Q-31,
     *           1.1922487600982223565416453283118294Q-34,
     *           2.6715112192401369859986789395807877Q-38,
     *           1.3386169421062562827190570142264186Q-42,
     *           4.5105361938989742322234283013237457Q-48/

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
      IMPLICIT REAL*16(A-H,O-Z)

!     THE INTEGRAND FOR SUBROUTINE GAUSS

      F=X/(1.+EXP(X))
      END

!     -------------------------------
 
      FUNCTION F2(X)
      IMPLICIT REAL*16(A-H,O-Z)

!     THE INTEGRAND*EXP(X) FOR SUBROUTINE LAGURE

      F2=X/(1+EXP(-X))
      END
