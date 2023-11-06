!     PROGRAM TO EVALUATE THE CAUCHY PRINCIPAL VALUE BY REWRITING THE INTEGRAL
!     AS EXPLAINED IN SECTION 6.6.10. THE INTEGRAND IS SUPPLIED BY
!     FUNCTION ROUTINE FUN, WHOSE NAME SHOULD NOT BE CHANGED
!     SINCE IT ALSO OCCURS IN FUNCTION FUNP.

      PROGRAM CAUCHYP
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL FUN,FUNP

!     EXERCISE 6.26 : INTEGRATE EXP(X)/X OVER [-1,1]

51    FORMAT('   IER =',i4,'  TOTAL INTEGRAL =',1PD14.6,5X,
     1  'NO. OF FUNCTION EVALUATIONS =',I7/5X,'ESTIMATED ERROR =',D14.6)
52    FORMAT('   A =',1PD14.6,5X,'B =',D14.6,5X,'C =',D14.6)

      REPS=1.D-13
      AEPS=1.D-15

!     THE  POSITION OF SINGULARITY
      C=0.0D0

100   PRINT *,'TYPE   A=LOWER LIMIT,  B=UPPER LIMIT,  C=SINGULARITY'
      PRINT *, '               (QUITS WHEN A.GE.B)'
      READ *,A,B,C
      IF(A.GE.B) STOP
      WRITE(6,52) A,B,C
      CALL CAUCHY(RI,A,B,C,REPS,AEPS,DIF,FUN,FUNP,IER,NPT)
      WRITE(6,51) IER,RI,NPT,DIF
      GO TO 100

      END

!	----------------------------------------------------------

!	To evaluate the Cauchy principal value of an integral over a finite interval
!
!	RI : (output) Calculated value of the integral
!	A : (input) The lower limit
!	B : (input) The upper limit (B > A)
!	C : (input) Location of the singularity (A < C < B)
!	REPS : (input) The required relative accuracy
!	AEPS : (input) The required absolute accuracy
!		The estimated error should be less than MAX(AEPS,REPS*ABS(RI))
!	DIF : (output) estimated (absolute) error achieved by the subroutine
!	F : (input) Name of the function routine to calculate the integrand
!	FUNP : (input) Name of the function routine to calculate F(C+x)+F(C-x)
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=304 implies A>B, A>C or C>B in which case no calculations
!			are done
!		Other values may be set by subroutine ADPINT which is called
!		twice. The returned value of IER is IER1+IER2*2, where IER1
!		and IER2 are the returned values of IER from two calls to ADPINT
!		In this case DIF will contain the estimated accuracy
!	NPT : (output) Number of function evaluations used by subroutine
!
!		FUNCTION F(X) must be supplied by the user.
!		FUNCTION FUNP(X) to calculate F(C+X)+F(C-X) should also be
!		supplied by the user. The value of C is available through
!		common block. Simplest version for FUNP may be
!
!		FUNCTION FUNP(X)
!		IMPLICIT REAL*8(A-H,O-Z)
!		COMMON/CAUFN/C
!		FUNP=F(C+X)+F(C-X)
!		END
!
!		If F(C+X)+F(C-X) can be combined roundoff error may be reduced.
!		There is no provision to pass the name F to FUNP, so it
!		will have to put explicitly.
!
!	Required routines : ADPINT, KRONRD, F, FUNP
 
      SUBROUTINE CAUCHY(RI,A,B,C,REPS,AEPS,DIF,F,FUNP,IER,NPT)
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL F,FUNP
!	To pass the value of C to FUNCTION F or FUNP
      COMMON/CAUFN/CC
 
      IF(A.GT.B.OR.A.GT.C.OR.C.GT.B) THEN
        IER=304
        RETURN
      ENDIF
 
!     FIRST EVALUATE THE SINGULAR PART
 
      CC=C
      R=MIN(C-A,B-C)
      AA=0.0
      NMAX=0
      NPT1=0
      DIF1=0.0
      RI1=0.0
      CALL ADPINT(RI1,AA,R,REPS,AEPS,DIF1,FUNP,IER,NPT1,NMAX)
 
!     EVALUATE THE REMAINING PORTION
 
      IF(C-A.GT.B-C) THEN
        AA=A
        BB=C-R
      ELSE
        AA=C+R
        BB=B
      ENDIF
      NPT2=0
      DIF2=0.0
      RI2=0.0
      IER1=0
      IF(ABS(BB-AA).GT.AEPS)
     1  CALL ADPINT(RI2,AA,BB,REPS,AEPS,DIF2,F,IER1,NPT2,NMAX)
 
      RI=RI1+RI2
!     FUNCTION FUNP REQUIRES TWO EVALUATIONS OF FUN
      NPT=2*NPT1+NPT2
      DIF=ABS(DIF1)+ABS(DIF2)
      IER=IER+2*IER1
      END

!	---------------------------

!		FUNCTION FUNP(X)
!		IMPLICIT REAL*8(A-H,O-Z)
!		COMMON/CAUFN/C
!		FUNP=F(C+X)+F(C-X)
!		END
 
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
 
      FUNCTION FUN(X)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/CAUFN/CC

!	The integrand

      FUN=EXP(X)/X
      END

!     --------------------------------------------------

      FUNCTION FUNP(X)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/CAUFN/CC

!      FUNP=FUN(CC+X)+FUN(CC-X)
!	This function has been simplified to avoid roundoff error
!	close to the singularity
      FUNP=2.*SINH(X)/X
      END
