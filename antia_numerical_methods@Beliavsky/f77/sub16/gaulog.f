!     To integrate a function with logarithmic singularity using Gaussian formulas
!	Since the weights and abscissas are not accurate to REAL*16
!	accuracy, this routine will not achive the maximum accuracy
!	permissible by arithmetic. It may be possible to achieve
!	higher accuracy using ADPINT with GAUS16.
!
!     RINT : (output) Calculated value of the integral
!     A : (input) The upper limit
!     AEPS : (input) The required absolute accuracy
!     REPS : (input) The required relative accuracy
!		The estimated error should be less than MAX(AEPS,REPS*ABS(RINT))
!     DIF : (output) estimated (absolute) error achieved by the subroutine
!     F : (input) Name of the function routine to calculate the
!		integrand (divided by LOG(A/X))
!     NPT : (output) Number of function evaluations used by the subroutine
!     IER : (output) Error parameter, IER=0 implies successful execution
!     	IER=30 implies specified accuracy was not achieved
!     		DIF will contain the estimated accuracy
!
!     Function F(X) must be supplied by the user
!     Note that subroutine calculates integral of F(X)*LOG(A/X)
!
!	Required routines : F

      SUBROUTINE GAULOG(RINT,A,AEPS,REPS,DIF,F,NPT,IER)
      IMPLICIT REAL*16(A-H,O-Z)
      DIMENSION W(30),X(30)
 
!     Weights and abscissas for Gaussian formula with logarithmic singularity
!     W(N-1),...,W(2N-2), are the weights for N-point rule and
!     X(N-1),...,X(2N-2), the corresponding abscissas
!     Weights and abscissas are available for N=2,4,8,16
 
      DATA (X(I),I=1,30)/1.1200880616697618295720548894768Q-01,
     *        6.0227690811873810275708022533804Q-01,
     *        4.1448480199383220803321310156357Q-02,
     *        2.4527491432060225193967575952330Q-01,
     *        5.5616545356027583718018435437603Q-01,
     *        8.4898239453298517464784918808468Q-01,
     *        1.3320244160892465012252672653536Q-02,
     *        7.9750429013894938409827729838524Q-02,
     *        1.9787102932618805379447616088820Q-01,
     *        3.5415399435190941967146360530582Q-01,
     *        5.2945857523491727770614970171970Q-01,
     *        7.0181452993909996383715267162352Q-01,
     *        8.4937932044110667604830920304715Q-01,
     *        9.5332645005635978876737967876114Q-01,
     *        3.8978344871159090954375938463874Q-03,
     *        2.3028945616873200451676266521095Q-02,
     *        5.8280398306240319722764201402098Q-02,
     *        1.0867836509105388173796912262945Q-01,
     *        1.7260945490984372439747684888913Q-01,
     *        2.4793705447057823633620061636912Q-01,
     *        3.3209454912991687050322528213811Q-01,
     *        4.2218391058194830849764583089701Q-01,
     *        5.1508247338146232501017958196979Q-01,
     *        6.0755612044772847474264194887061Q-01,
     *        6.9637565322821385232564092866513Q-01,
     *        7.7843256587326524309414044521296Q-01,
     *        8.5085026971539096880619590573013Q-01,
     *        9.1108685722227183475290609868700Q-01,
     *        9.5702557170354212258034366596252Q-01,
     *        9.8704780024798446604416946162482Q-01/

      DATA (W(I),I=1,30)/7.1853931903038444066551020089099Q-01,
     *        2.8146068096961555933448979910901Q-01,
     *        3.8346406814513512485004652234303Q-01,
     *        3.8687531777476262733600823455435Q-01,
     *        1.9043512695014241536136001454740Q-01,
     *        8.4898239453298517464784918808468Q-01,
     *        1.6441660472800288683147256953944Q-01,
     *        2.3752561002330602050134856290643Q-01,
     *        2.2684198443191912636878040291254Q-01,
     *        1.7575407900607024498805621134754Q-01,
     *        1.1292403024675905185500044134637Q-01,
     *        5.7872210717782072398527966800967Q-02,
     *        2.0979073742132978043461523909241Q-02,
     *        3.6864071040276190133523212374786Q-03,
     *        6.0791710043591145085177022050780Q-02,
     *        1.0291567751758202281665163290891Q-01,
     *        1.2235566204600909185180938370661Q-01,
     *        1.2756924693701593233951253210084Q-01,
     *        1.2301357460007090831897038729127Q-01,
     *        1.1184724485548575522933168696643Q-01,
     *        9.6596385152124398489528192055775Q-02,
     *        7.9356664351473205727177277568989Q-02,
     *        6.1850494581965271972248730674684Q-02,
     *        4.5435246507726723814277153028712Q-02,
     *        3.1098974751581848285136209225991Q-02,
     *        1.9459765927360870294836740237693Q-02,
     *        1.0776254963205542126752129705866Q-02,
     *        4.9725428900876496102471599749318Q-03,
     *        1.6782011100511972493950116557883Q-03,
     *        2.8235376466843678894875084672700Q-04/
 
      IER=0
!     The 2-point formula
      R1=(F(A*X(1))*W(1)+F(A*X(2))*W(2))*A
      NPT=2
      N=2
 
!     Use higher order formula until convergence
      DO 2000 J=2,4
        N=N*2
        R2=0.0
        DO 1000 I=N-1,2*N-2
1000    R2=R2+F(X(I)*A)*W(I)
        R2=R2*A
 
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
