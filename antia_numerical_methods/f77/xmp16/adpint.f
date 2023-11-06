!     PROGRAM FOR ADAPTIVE INTEGRATION OVER A FINITE INTERVAL
!	This program uses GAUS16 instead of KRONRD and hence
!	is less efficient as compared to the version used in quad.f

      PROGRAM QUAD
      IMPLICIT REAL*16(A-H,O-Z)
      EXTERNAL FUN

!     EXAMPLE 6.1: INTEGRATE SQRT(X) OVER [0,1]

53    FORMAT('   A =',1PD14.6,5X,'B =',D14.6,5X,'IER =',I4/
     1       '    NO. OF FUNCTION EVALUATIONS =',I7,5X,'INTEGRAL =',
     2       D14.6/5X,'ESTIMATED ERROR =',D14.6)

      REPS=1.Q-24
      AEPS=1.Q-35

100   PRINT *,'TYPE   A=LOWER LIMIT,  B=UPPER LIMIT,',
     1        '    (QUITS WHEN A.EQ.B)'
      READ *,A,B
      IF(A.EQ.B) STOP
      NMAX=16385
      CALL ADPINT(RI,A,B,REPS,AEPS,DIF,FUN,IER,NPT,NMAX)
      WRITE(6,53) A,B,IER,NPT,RI,DIF
      GO TO 100
      END
 

!     ----------------------------------------------------------

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
!
!	The weights in KRONRD are accurate only to REAL*16 precision and hence if
!	higher precision is required it will be preferable to use GAUS16, although
!	it is less efficient.

      SUBROUTINE ADPINT(RINT,XL,XU,REPS,AEPS,DIF,F,IER,NPT,NMAX)
      IMPLICIT REAL*16(A-H,O,P,R-Z)
      IMPLICIT LOGICAL(Q)
      PARAMETER(IPMAX=200,IFMAX=5,MAXPT=100000)
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
!1000  CALL KRONRD(FINT,RL,RU,DIF0,NP,F)
1000  CALL GAUS16(FINT,RL,RU,DIF0,NP,F)
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
!      CALL KRONRD(FINT,RL,RU,DIF0,NP,F)
      CALL GAUS16(FINT,RL,RU,DIF0,NP,F)
      NPT=NPT+NP
      RINT=RINT+FINT
      DIF=DIF+DIF0
      END

!     ----------------------------------------------------------

!	To integrate a function over a finite interval using 16 point
!	Gauss-Legendre formula, for use with ADPINT
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

      SUBROUTINE GAUS16(RI,A,B,DIF,N,F)
      IMPLICIT REAL*16(A-H,O-Z)
      DIMENSION  W8(4),A8(4),W16(8),A16(8)

!	W8 and A8 are the weights and abscissas for the 8-point Gauss formula
!	W16 and A16 are the weights and abscissas for the 16-point Gauss formula
!	Because of symmetry only half the points are given.

      DATA W8 /  1.0122853629037625915253135430996220Q-1,
     *           2.2238103445337447054435599442624126Q-1,
     *           3.1370664587788728733796220198660110Q-1,
     *           3.6268378337836198296515044927719567Q-1/
      DATA A8 /  9.6028985649753623168356086856947230Q-1,
     *           7.9666647741362673959155393647583006Q-1,
     *           5.2553240991632898581773904918924604Q-1,
     *           1.8343464249564980493947614236018303Q-1/

      DATA W16/  2.7152459411754094851780572456019197Q-2,
     *           6.2253523938647892862843836994377013Q-2,
     *           9.5158511682492784809925107602247148Q-2,
     *           1.2462897125553387205247628219201544Q-1,
     *           1.4959598881657673208150173054747880Q-1,
     *           1.6915651939500253818931207903036093Q-1,
     *           1.8260341504492358886676366796921865Q-1,
     *           1.8945061045506849628539672320828440Q-1/
      DATA A16/  9.8940093499164993259615417345033306Q-1,
     *           9.4457502307323257607798841553460820Q-1,
     *           8.6563120238783174388046789771239397Q-1,
     *           7.5540440835500303389510119484744268Q-1,
     *           6.1787624440264374844667176404879140Q-1,
     *           4.5801677765722738634241944298357810Q-1,
     *           2.8160355077925891323046050146049710Q-1,
     *           9.5012509837637440185319335424958114Q-2/

      AT=(B-A)/2.
      BT=(B+A)/2.
      R1=0.0
!	8-point Gauss-Legendre formula
      DO 2000 K=1,4
        R1=R1+W8(K)*(F(AT*A8(K)+BT)+F(BT-AT*A8(K)))
2000  CONTINUE

      RI=0.0
!	16-point Gauss-Legendre formula
      DO 2500 K=1,8
        RI=RI+W16(K)*(F(AT*A16(K)+BT)+F(BT-AT*A16(K)))
2500  CONTINUE

      RI=RI*AT
      R1=R1*AT
      DIF=ABS(RI-R1)
      N=24
      END

!     ----------------------------------------------------------
 
      FUNCTION FUN(X)
      IMPLICIT REAL*16(A-H,O-Z)

!     THE INTEGRAND

      FUN=SQRT(X)
      END
