!     PROGRAM FOR INTEGRATION OF F(X)*EXP(-X*X) OVER infinite interval

      PROGRAM INTEG
      IMPLICIT REAL*16(A-H,O-Z)
      EXTERNAL F

!     EXERCISE 6.24 (I_1)

51    FORMAT('   IER =',I4,5X,
     1  'NO. OF FUNCTION EVALUATIONS =',I7/5X,'INTEGRAL =',1PD14.6,5X,
     2       'ESTIMATED ERROR =',D14.6/5X,'EXACT VALUE =',D14.6)

      REPS=1.Q-33
      AEPS=1.Q-39
      PI=2.*ACOS(0.Q0)
      REX=SQRT(PI)*EXP(-0.25Q0)

      CALL HERMIT(RI,AEPS,REPS,DIF,F,N,IER)
      WRITE(6,51) IER,N,RI,DIF,REX
      END
 
!     ---------------------------------------------

!     To integrate a function over infinite interval using Gauss-Hermite formulas
!
!     RINT : (output) Calculated value of the integral
!     AEPS : (input) The required absolute accuracy
!     REPS : (input) The required relative accuracy
!     	The estimated error should be less than MAX(AEPS,REPS*ABS(RINT))
!     DIF : (output) estimated (absolute) error achieved by the subroutine
!     F : (input) Name of the function routine to calculate the
!     	integrand (multiplied by EXP(X**2))
!     NPT : (output) Number of function evaluations used by the subroutine
!     IER : (output) Error parameter, IER=0 implies successful execution
!     	IER=30 implies specified accuracy was not achieved
!     		DIF will contain the estimated accuracy
!
!     Function F(X) must be supplied by the user
!
!	Required routines : F

      SUBROUTINE HERMIT(RINT,AEPS,REPS,DIF,F,NPT,IER)
      IMPLICIT REAL*16(A-H,O-Z)
      DIMENSION W(31),X(31)
 
!     Weights and abscissas for Gauss-Hermite quadrature
!     Because of symmetry only positive abscissas are listed
!     W(N/2),...,W(N-1), are the weights for N-point rule and
!     X(N/2),...,X(N-1), the corresponding abscissas
!     Weights and abscissas are available for N=2,4,8,16,32
 

      DATA X  /  7.0710678118654752440084436210484918Q-01,
     *           5.2464762327529031788406025383474093Q-01,
     *           1.6506801238857845558833411111207452Q+00,
     *           3.8118699020732211685471888558369149Q-01,
     *           1.1571937124467801947207657790631012Q+00,
     *           1.9816567566958429258546306397693112Q+00,
     *           2.9306374202572440192235027052435995Q+00,
     *           2.7348104613815245215828040196501441Q-01,
     *           8.2295144914465589258245449673394312Q-01,
     *           1.3802585391988807963720896696945812Q+00,
     *           1.9517879909162539774346554149598867Q+00,
     *           2.5462021578474813621593287054458959Q+00,
     *           3.1769991619799560268139945592637001Q+00,
     *           3.8694479048601226987194240980148130Q+00,
     *           4.6887389393058183646884986487456130Q+00,
     *           1.9484074156939932670874128953220071Q-01,
     *           5.8497876543593244846695754401145483Q-01,
     *           9.7650046358968283848470487198177257Q-01,
     *           1.3703764109528718381617056486445837Q+00,
     *           1.7676541094632016046276732585265015Q+00,
     *           2.1694991836061121733057055950150543Q+00,
     *           2.5772495377323174540309293011413068Q+00,
     *           2.9924908250023742062854940760637751Q+00,
     *           3.4171674928185707358739272956352736Q+00,
     *           3.8537554854714446438878729210936244Q+00,
     *           4.3055479533511984452634865319312600Q+00,
     *           4.7771645035025963930357940568894015Q+00,
     *           5.2755509865158801278190604813963664Q+00,
     *           5.8122259495159138327659661536595846Q+00,
     *           6.4094981492696604121737637415270924Q+00,
     *           7.1258139098307275727952076034223051Q+00/

      DATA W  /  8.8622692545275801364908374167057276Q-01,
     *           8.0491409000551283650604918448068342Q-01,
     *           8.1312835447245177143034557189888271Q-02,
     *           6.6114701255824129103041597449588209Q-01,
     *           2.0780232581489187954325862028570085Q-01,
     *           1.7077983007413475456203056436445699Q-02,
     *           1.9960407221136761920609045254409736Q-04,
     *           5.0792947901661374191351734179052181Q-01,
     *           2.8064745852853367536946333537965308Q-01,
     *           8.3810041398985829415420734900115654Q-02,
     *           1.2880311535509973683464299931172100Q-02,
     *           9.3228400862418052991427730553690579Q-04,
     *           2.7118600925378815120189143224359391Q-05,
     *           2.3209808448652106533874942318486174Q-07,
     *           2.6548074740111822447092636605021063Q-10,
     *           3.7523835259280239286681838890711916Q-01,
     *           2.7745814230252989813769891854161759Q-01,
     *           1.5126973407664248257514711461402814Q-01,
     *           6.0458130955912614186585760783289777Q-02,
     *           1.7553428831573430303437844611009365Q-02,
     *           3.6548903266544280791256571224128517Q-03,
     *           5.3626836552797204597023810150086303Q-04,
     *           5.4165840618199825580019393926719719Q-05,
     *           3.6505851295623760573703241874583511Q-06,
     *           1.5741677925455940292686925792946215Q-07,
     *           4.0988321647708966182350410137972719Q-09,
     *           5.9332914633966386145115682155834451Q-11,
     *           4.2150102113264475729694452118301675Q-13,
     *           1.1973440170928486658286818995090002Q-15,
     *           9.2317365365182922334944200720747836Q-19,
     *           7.3106764273841623932742784550630736Q-23/
 
 
      IER=0
!     The 2-point formula
      R1=W(1)*(F(X(1))+F(-X(1)))
      NPT=2
      N=2
 
!     Use higher order formula until convergence
      DO 2000 J=2,5
        N=N*2
        R2=0.0
        DO 1000 I=N/2,N-1
1000    R2=R2+(F(X(I))+F(-X(I)))*W(I)
 
        NPT=NPT+N
        DIF=R2-R1
        RINT=R2
        IF(ABS(DIF).LT.MAX(AEPS,REPS*ABS(RINT))) RETURN
        R1=R2
2000  CONTINUE
 
!     Integral fails to converge
      IER=30
      RETURN
      END

!     ---------------------------------------------

      FUNCTION F(X)
      IMPLICIT REAL*16(A-H,O-Z)

!     THE INTEGRAND*EXP(X*X) FOR SUBROUTINE HERMIT

      F=COS(X)
      END
