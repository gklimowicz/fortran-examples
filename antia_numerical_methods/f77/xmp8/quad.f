!     PROGRAM FOR INTEGRATION OVER A FINITE INTERVAL
!     USING SIMPSON'S RULE OR ROMBRG INTEGRATION OR
!     EPSILON ALGORITHM OR GAUSS-LEGENDRE FORMULAS OR
!     ADAPTIVE INTEGRATION

      PROGRAM QUAD
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION GI(20)
      EXTERNAL FUN

!     EXAMPLE 6.1: INTEGRATE SQRT(X) OVER [0,1]

51    FORMAT('   EXPONENTS IN ERROR EXPANSION:'/13F6.2)
52    FORMAT(I8,' POINT COMPOSITE GAUSS-LEGENDRE FORMULA')
53    FORMAT('   A =',1PD14.6,5X,'B =',D14.6,5X,'IER =',I4,5X,'IT=',I2/
     1       '    NO. OF FUNCTION EVALUATIONS =',I7,5X,'INTEGRAL =',
     2       D14.6/5X,'ESTIMATED ERROR =',2D14.6)

      REPS=1.D-13
      AEPS=1.D-18

100   PRINT *,'TYPE   A=LOWER LIMIT,  B=UPPER LIMIT,',
     1        '    (QUITS WHEN A.EQ.B)'
      READ *,A,B
      IF(A.EQ.B) STOP
      PRINT *,'TYPE IT=1/2/3/4/5  FOR SIMSON/ROMBRG/EPSILN/GAUSS/ADPINT'
      READ *,IT
      IF(IT.LT.1.OR.IT.GT.5) STOP
      REX=(B**1.5D0-A**1.5D0)/1.5D0
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

      SUBROUTINE ADPINT(RINT,XL,XU,REPS,AEPS,DIF,F,IER,NPT,NMAX)
      IMPLICIT REAL*8(A-H,O,P,R-Z)
      IMPLICIT LOGICAL(Q)
      PARAMETER(IPMAX=100,IFMAX=5,MAXPT=100000)
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
      IMPLICIT REAL*8(A-H,O-Z)
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
      IMPLICIT REAL*8(A-H,O-Z)
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

!     ----------------------------------------------------------

!	To integrate a function over a finite interval using Gauss-Kronrod formula
!	For use with ADPINT
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
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION  W7(4),A7(4),WK7(4),WK15(4),AK15(4)

!	W7 and A7 are the weights and abscissas for the 7-point Gauss formula
!	WK7 are the weights for these points in Kronrod formula
!	WK15 and AK15 are the weights and abscissas for the remaining points
!	in Kronrod formula.
!	Because of symmetry only half the points are given.

      DATA W7  /0.12948496616886969327D0, 0.27970539148927666790D0,
     *          0.38183005050511894495D0, 0.41795918367346938775D0/
      DATA A7  /0.94910791234275852452D0, 0.74153118559939443986D0,
     *          0.40584515137739716690D0, 0.0/
      DATA WK7 /0.06309209262997855329D0, 0.14065325971552591874D0,
     *          0.19035057806478540991D0, 0.20948214108472782801D0/
      DATA WK15/0.02293532201052922496D0, 0.10479001032225018383D0,
     *          0.16900472663926790282D0, 0.20443294007529889241D0/
      DATA AK15/0.99145537112081263920D0, 0.86486442335976907278D0,
     *          0.58608723546769113029D0, 0.20778495500789846760D0/

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
      IMPLICIT REAL*8(A-H,O-Z)
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
      IMPLICIT REAL*8(A-H,O-Z)
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
      IMPLICIT REAL*8(A-H,O-Z)

!     THE INTEGRAND

      FUN=SQRT(X)
      END
