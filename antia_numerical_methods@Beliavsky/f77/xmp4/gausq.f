!     PROGRAM FOR INTEGRATING F(X)/SQRT(X) OVER (0,A] 

      PROGRAM INTEG
!      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL F,F2

!     EXERCISE 6.14

51    FORMAT('   IER =',I4,5X,'A =',1PE14.6,5X,'A1 =',E14.6,5X,
     1       'NO. OF FUNCTION EVALUATIONS =',I7/5X,'INTEGRAL =',
     2       E14.6,5X,'ESTIMATED ERROR =',3E14.6)

      REPS=1.E-7
      AEPS=1.E-9

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
!      IMPLICIT REAL*8(A-H,O-Z)
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
      IER=0
!	The 1-point formula
      R1=2.D0*W(1)*FUN(A*X(1)**2)
      R1=R1*SQRT(A)
      NPT=0
 
!	Try higher order formulas until convergence 
      DO 3000 I=2,5
        N=N*2
 
        S1=0.0
        DO 2000 K=N,2*N-1
2000    S1=S1+2.D0*W(K)*FUN(A*(X(K))**2)
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
!      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(AMN=1.D-2)
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

      F=SQRT(X)*COS(X)
      END

!     -------------------------------
 
      FUNCTION F2(X)
!      IMPLICIT REAL*8(A-H,O-Z)

!     THE INTEGRAND*SQRT(X) FOR USE WITH SUBROUTINE GAUSQ

      F2=X*COS(X)
      END

