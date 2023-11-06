!     PROGRAM FOR INTEGRATION OF F(X)*EXP(-X*X) OVER infinite interval

      PROGRAM INTEG
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL F

!     EXERCISE 6.24 (I_1)

51    FORMAT('   IER =',I4,5X,
     1  'NO. OF FUNCTION EVALUATIONS =',I7/5X,'INTEGRAL =',1PD14.6,5X,
     2       'ESTIMATED ERROR =',D14.6/5X,'EXACT VALUE =',D14.6)

      REPS=1.D-14
      AEPS=1.D-19
      PI=2.*ACOS(0.D0)
      REX=SQRT(PI)*EXP(-0.25D0)

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
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION W(31),X(31)
 
!     Weights and abscissas for Gauss-Hermite quadrature
!     Because of symmetry only positive abscissas are listed
!     W(N/2),...,W(N-1), are the weights for N-point rule and
!     X(N/2),...,X(N-1), the corresponding abscissas
!     Weights and abscissas are available for N=2,4,8,16,32
 
      DATA (X(I),I=1,31)/   0.707106781186548D0,  0.524647623275290D0,
     *1.650680123885785D0,  0.3811869902073221D0, 1.1571937124467802D0,
     *1.98165675669584293D0,2.93063742025724402D0,0.27348104613815245D0,
     *0.82295144914465589D0,1.3802585391988808D0, 1.9517879909162540D0,
     *2.54620215784748136D0,3.17699916197995603D0,3.8694479048601227D0,
     *4.68873893930581836D0,0.19484074156939933D0,0.58497876543593245D0,
     *0.97650046358968284D0,1.37037641095287184D0,1.76765410946320160D0,
     *2.16949918360611217D0,2.57724953773231745D0,2.99249082500237421D0,
     *3.41716749281857074D0,3.85375548547144464D0,4.30554795335119845D0,
     *4.77716450350259639D0,5.27555098651588013D0,5.81222594951591383D0,
     *6.409498149269660412D0,7.1258139098307275728D0/

      DATA (W(I),I=1,31)/   0.886226925452758D0,  0.804914090005513D0,
     *8.13128354472452D-2,  0.6611470125582413D0, 0.207802325814892D0,
     *1.707798300741347D-2, 1.9960407221136762D-4,0.50792947901661374D0,
     *0.2806474585285342D0, 8.3810041398985349D-2,1.2880311535509972D-2,
     *9.3228400862418052D-4,2.7118600925378815D-5,2.3209808448652106D-7,
     *2.654807474011182D-10,0.37523835259280239D0,0.27745814230252993D0,
     *0.15126973407664245D0,6.0458130955912614D-2,1.7553428831573430D-2,
     *3.6548903266544280D-3,5.3626836552797204D-4,5.4165840618199826D-5,
     *3.6505851295623761D-6,1.5741677925455940D-7,4.0988321647708966D-9,
     *5.933291463396639D-11,4.215010211326448D-13,1.197344017092849D-15,
     *9.2317365365182922D-19,7.3106764273841624D-23/
 
 
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
      IMPLICIT REAL*8(A-H,O-Z)

!     THE INTEGRAND*EXP(X*X) FOR SUBROUTINE HERMIT

      F=COS(X)
      END
