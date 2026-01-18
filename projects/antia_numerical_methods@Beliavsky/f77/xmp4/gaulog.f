!     PROGRAM FOR INTEGRATING F(X)*LOG(X) OVER (0,A] 

      PROGRAM INTEG
!      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL F,F2

!     EXERCISE 6.2 (I_5)

51    FORMAT('   IER =',I4,5X,'A =',1PE14.6,5X,'A1 =',E14.6/
     1      'NO. OF FUNCTION EVALUATIONS =',I7/5X,'INTEGRAL =',E14.6,5X,
     2      'ESTIMATED ERROR =',3E14.6)
52    FORMAT('  INITIAL VALUE OF A1 = ',1PE14.6)

      REPS=1.E-6
      AEPS=1.E-9

100   PRINT *,'TYPE  A=UPPER LIMIT,  A1=POINT AT WHICH INTEGRAL IS TO'
     1       , ' BE SPLIT'
      PRINT *,'                          (QUITS WHEN A<0)'
      READ *,A,A1
      IF(A.LT.0.0) STOP
      WRITE(6,52) A1
      CALL GAULG2(RI,A,A1,REPS,AEPS,DIF,F,F2,N,IER)
      WRITE(6,51) IER,A,A1,N,RI,DIF
      GO TO 100
      END
 
!     ---------------------------------------------------
 
!     To integrate a function with logarithmic singularity over (0,A]
!     using a combination of Gaussian formulas
!
!     RINT : (output) Calculated value of the integral
!     A : (input) The upper limit
!     A1 : (input/output) The point at which integral has to be broken
!		A1 will be adjusted by the subroutine.
!     REPS : (input) The required relative accuracy
!     AEPS : (input) The required absolute accuracy
!		The estimated error should be less than MAX(AEPS,REPS*ABS(RINT))
!     DIF : (output) estimated (absolute) error achieved by the subroutine
!     F : (input) Name of the function routine to calculate the integrand
!     F2 : (input) Name of the function routine to calculate F(X)/LOG(X)
!     NP : (output) Number of function evaluations used by the subroutine
!     IER : (output) Error parameter, IER=0 implies successful execution
!		IER=31 implies specified accuracy was not achieved by GAUSS over [A1,A]
!		IER=32 implies specified accuracy was not achieved by GAULOG
!		IER=34 implies specified accuracy was not achieved by GAUSS over (0,A1]
!		In case of multiple failures second digit of IER will be sum
!			of these values.
!		In all cases DIF will contain the estimated accuracy
!
!     	FUNCTION F(X) and F2(X) must be supplied by the user.
!
!     Required routines : GAUSS, GAULOG, F, F2
 
      SUBROUTINE GAULG2(RINT,A,A1,REPS,AEPS,DIF,F,F2,NP,IER)
!      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(AMN=1.D-2)
      EXTERNAL F,F2
 
      IER1=0
      IER2=0
      R1=0.0
      R2=0.0
      R3=0.0
      A2=0.0
      DIF=0.0
      DIF2=0.0
      NP=0
      NPT=0
      NPT2=0
      IF(A1.GT.A) A1=A
      IF(A1.LE.0.0) A1=A
 
!     Evaluate the integral over (0,A1]
2200  CALL GAULOG(R1,A1,AEPS,REPS,DIF1,F2,NPT1,IER)
      NP=NP+NPT1
      IF(IER.EQ.0) GO TO 2500
!     If GAULOG fails decrease A1
      T1=A1
      A1=A1/2.
      IF(A1.GT.AMN) GO TO 2200
!     To prevent infinite loop do not reduce A1 below AMN
      IER1=2
      IER=0
      A1=T1
 
!     Evaluate the integral over [A1,A]
2500  IF(A-A1.GT.AEPS) CALL GAUSS(R2,A1,A,16,REPS,AEPS,DIF,IER,NPT,F)
      IF(IER.GT.0) IER=1
 
!     Evaluate the regular part over [0,A1]
      IF(A1.NE.1) CALL GAUSS(R3,A2,A1,16,REPS,AEPS,DIF2,IER2,NPT2,F2)
      IF(IER2.GT.0) IER2=4
      IER=IER+IER1+IER2
      IF(IER.GT.0) IER=IER+30
      RINT=R1+R2-R3*LOG(A1)
      DIF=ABS(DIF)+ABS(DIF1)+ABS(DIF2)
      NP=NP+NPT+NPT2
      END
 
!     ---------------------------------------------

!     To integrate a function with logarithmic singularity using Gaussian formulas
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
!      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION W(30),X(30)
 
!     Weights and abscissas for Gaussian formula with logarithmic singularity
!     W(N-1),...,W(2N-2), are the weights for N-point rule and
!     X(N-1),...,X(2N-2), the corresponding abscissas
!     Weights and abscissas are available for N=2,4,8,16
 
      DATA (X(I),I=1,30)/1.120088061669761830D-1,6.022769081187381028D-1
     *,       4.144848019938322080D-2, 2.452749143206022519D-1,
     *        5.561654535602758372D-1, 8.489823945329851746D-1,
     *        1.332024416089246501D-2, 7.975042901389493841D-2,
     *        1.978710293261880538D-1, 3.541539943519094197D-1,
     *        5.294585752349172777D-1, 7.018145299390999638D-1,
     *        8.493793204411066760D-1, 9.533264500563597888D-1,
     *        3.897834487115909095D-3, 2.302894561687320045D-2,
     *        5.828039830624031972D-2, 1.086783650910538817D-1,
     *        1.726094549098437244D-1, 2.479370544705782363D-1,
     *        3.320945491299168705D-1, 4.221839105819483085D-1,
     *        5.150824733814623250D-1, 6.075561204477284747D-1,
     *        6.963756532282138523D-1, 7.784325658732652431D-1,
     *        8.508502697153909688D-1, 9.110868572222718348D-1,
     *        9.570255717035421226D-1, 9.870478002479844660D-1/

      DATA (W(I),I=1,30)/7.185393190303844407D-1,2.814606809696155593D-1
     *,       3.834640681451351249D-1, 3.868753177747626273D-1,
     *        1.904351269501424154D-1, 3.922548712995983245D-2,
     *        1.644166047280028868D-1, 2.375256100233060205D-1,
     *        2.268419844319191264D-1, 1.757540790060702450D-1,
     *        1.129240302467590519D-1, 5.787221071778207240D-2,
     *        2.097907374213297804D-2, 3.686407104027619013D-3,
     *        6.079171004359114509D-2, 1.029156775175820228D-1,
     *        1.223556620460090919D-1, 1.275692469370159323D-1,
     *        1.230135746000709083D-1, 1.118472448554857552D-1,
     *        9.659638515212439849D-2, 7.935666435147320573D-2,
     *        6.185049458196527197D-2, 4.543524650772672381D-2,
     *        3.109897475158184829D-2, 1.945976592736087029D-2,
     *        1.077625496320554213D-2, 4.972542890087649610D-3,
     *        1.678201110051197249D-3, 2.823537646684367889D-4/
 
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

!     --------------------------------------------------------

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
!      IMPLICIT REAL*8(A-H,O-Z)
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

!     ---------------------------------------------
 
      FUNCTION F(X)
!      IMPLICIT REAL*8(A-H,O-Z)

!     THE INTEGRAND FOR SUBROUTINE GAUSS

      F=-LOG(X)*SIN(X)
      END

!     -------------------------------
 
      FUNCTION F2(X)
!      IMPLICIT REAL*8(A-H,O-Z)

!     THE INTEGRAND/LOG(1/X) FOR SUBROUTINE GAULOG

      F2=SIN(X)
      END

