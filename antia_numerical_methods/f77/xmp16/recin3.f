!     PROGRAM FOR RECURSIVE EVALUATION OF TRIPLE INTEGRAL
!     ADPIN1, ADPIN2 AND KRONR1, KRONR2 ARE COPIES OF ADPINT AND KRONRD
!     TO AVOID RECURSIVE CALLS.

      PROGRAM INTEG
      IMPLICIT REAL*16(A-H,O-Z)
      EXTERNAL F,F1,F2
      COMMON/FN/ALP,XX,YY,NP

!     EXERCISE 6.44 : INTEGRATION OVER ELLIPSOID

51    FORMAT('  IER =',I4,5X,'ALP =',1PD14.6,5X,
     1 'NO. OF FUNCTION EVALUATIONS =',I10/5X,'INTEGRAL =',D14.6,
     2 5X,'ESTIMATED ERROR =',D14.6)

      REPS=1.Q-20
      AEPS=1.Q-21

!     LIMITS FOR INTEGRATION W.R.T. X
      A=-1
      B=1

!     ALP IS THE COORDINATE A IN THE EXERCISE
100   PRINT *,'TYPE A,      (QUITS WHEN A<0)'
      READ *,ALP
      IF(ALP.LT.0.0) STOP

!     THE NUMBER OF FUNCTION EVALUATIONS ARE ACCUMULATED IN THE
!     VARIABLE NP IN THE COMMON BLOCK
      NP=0
      NMAX=6000

!     INTEGRATE W.R.T. X

      CALL ADPINT(RI,A,B,REPS,AEPS,DIF,F,IER,NPT,NMAX)
      WRITE(6,51) IER,ALP,NP,RI,DIF
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

!   -----------------------------------------------------

!	To integrate a function over a finite interval using Gauss-Kronrod formula
!	For use with ADPINT
!	Since the weights and abscissas are not accurate to REAL*16
!	accuracy, this routine will not achive the maximum accuracy
!	permissible by arithmetic. Use GAUS16 instead of KRONRD
!	if high accuracy is required.
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
      IMPLICIT REAL*16(A-H,O-Z)
      DIMENSION  W7(4),A7(4),WK7(4),WK15(4),AK15(4)

!	W7 and A7 are the weights and abscissas for the 7-point Gauss formula
!	WK7 are the weights for these points in Kronrod formula
!	WK15 and AK15 are the weights and abscissas for the remaining points
!	in Kronrod formula.
!	Because of symmetry only half the points are given.

      DATA W7  /0.1294849661688696932706114326790820183Q0,
     *          0.2797053914892766679014677714237795825Q0,
     *          0.3818300505051189449503697754889751339Q0,
     *          0.4179591836734693877551020408163265306Q0/
      DATA A7  /0.9491079123427585245261896840478512624Q0,
     *          0.7415311855993944398638647732807884070Q0,
     *          0.4058451513773971669066064120769614633Q0, 0.0/
      DATA WK7 /0.0630920926299785532907006631892042866Q0,
     *          0.1406532597155259187451895905102379204Q0,
     *          0.1903505780647854099132564024210136828Q0,
     *          0.2094821410847278280129991748917142637Q0/
      DATA WK15/0.0229353220105292249637320080589695920Q0,
     *          0.1047900103222501838398763225415180174Q0,
     *          0.1690047266392679028265834265985502841Q0,
     *          0.2044329400752988924141619992346490847Q0/
      DATA AK15/0.9914553711208126392068546975263285166Q0,
     *          0.8648644233597690727897127886409262012Q0,
     *          0.5860872354676911302941448382587295984Q0,
     *          0.2077849550078984676006894037732449134Q0/

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
!
!	The weights in KRONRD are accurate only to REAL*16 precision and hence if
!	higher precision is required it will be preferable to use GAUS16, although
!	it is less efficient.

      SUBROUTINE ADPIN1(RINT,XL,XU,REPS,AEPS,DIF,F,IER,NPT,NMAX)
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
1000  CALL KRONRD1(FINT,RL,RU,DIF0,NP,F)
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
      CALL KRONRD1(FINT,RL,RU,DIF0,NP,F)
!      CALL GAUS16(FINT,RL,RU,DIF0,NP,F)
      NPT=NPT+NP
      RINT=RINT+FINT
      DIF=DIF+DIF0
      END

!   -----------------------------------------------------

!	To integrate a function over a finite interval using Gauss-Kronrod formula
!	For use with ADPINT
!	Since the weights and abscissas are not accurate to REAL*16
!	accuracy, this routine will not achive the maximum accuracy
!	permissible by arithmetic. Use GAUS16 instead of KRONRD
!	if high accuracy is required.
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

      SUBROUTINE KRONRD1(RI,A,B,DIF,N,F)
      IMPLICIT REAL*16(A-H,O-Z)
      DIMENSION  W7(4),A7(4),WK7(4),WK15(4),AK15(4)

!	W7 and A7 are the weights and abscissas for the 7-point Gauss formula
!	WK7 are the weights for these points in Kronrod formula
!	WK15 and AK15 are the weights and abscissas for the remaining points
!	in Kronrod formula.
!	Because of symmetry only half the points are given.

      DATA W7  /0.1294849661688696932706114326790820183Q0,
     *          0.2797053914892766679014677714237795825Q0,
     *          0.3818300505051189449503697754889751339Q0,
     *          0.4179591836734693877551020408163265306Q0/
      DATA A7  /0.9491079123427585245261896840478512624Q0,
     *          0.7415311855993944398638647732807884070Q0,
     *          0.4058451513773971669066064120769614633Q0, 0.0/
      DATA WK7 /0.0630920926299785532907006631892042866Q0,
     *          0.1406532597155259187451895905102379204Q0,
     *          0.1903505780647854099132564024210136828Q0,
     *          0.2094821410847278280129991748917142637Q0/
      DATA WK15/0.0229353220105292249637320080589695920Q0,
     *          0.1047900103222501838398763225415180174Q0,
     *          0.1690047266392679028265834265985502841Q0,
     *          0.2044329400752988924141619992346490847Q0/
      DATA AK15/0.9914553711208126392068546975263285166Q0,
     *          0.8648644233597690727897127886409262012Q0,
     *          0.5860872354676911302941448382587295984Q0,
     *          0.2077849550078984676006894037732449134Q0/

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

!     ------------------------------------

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

      SUBROUTINE ADPIN2(RINT,XL,XU,REPS,AEPS,DIF,F,IER,NPT,NMAX)
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
1000  CALL KRONRD2(FINT,RL,RU,DIF0,NP,F)
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
      CALL KRONRD2(FINT,RL,RU,DIF0,NP,F)
!      CALL GAUS16(FINT,RL,RU,DIF0,NP,F)
      NPT=NPT+NP
      RINT=RINT+FINT
      DIF=DIF+DIF0
      END

!   -----------------------------------------------------

!	To integrate a function over a finite interval using Gauss-Kronrod formula
!	For use with ADPINT
!	Since the weights and abscissas are not accurate to REAL*16
!	accuracy, this routine will not achive the maximum accuracy
!	permissible by arithmetic. Use GAUS16 instead of KRONRD
!	if high accuracy is required.
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

      SUBROUTINE KRONRD2(RI,A,B,DIF,N,F)
      IMPLICIT REAL*16(A-H,O-Z)
      DIMENSION  W7(4),A7(4),WK7(4),WK15(4),AK15(4)

!	W7 and A7 are the weights and abscissas for the 7-point Gauss formula
!	WK7 are the weights for these points in Kronrod formula
!	WK15 and AK15 are the weights and abscissas for the remaining points
!	in Kronrod formula.
!	Because of symmetry only half the points are given.

      DATA W7  /0.1294849661688696932706114326790820183Q0,
     *          0.2797053914892766679014677714237795825Q0,
     *          0.3818300505051189449503697754889751339Q0,
     *          0.4179591836734693877551020408163265306Q0/
      DATA A7  /0.9491079123427585245261896840478512624Q0,
     *          0.7415311855993944398638647732807884070Q0,
     *          0.4058451513773971669066064120769614633Q0, 0.0/
      DATA WK7 /0.0630920926299785532907006631892042866Q0,
     *          0.1406532597155259187451895905102379204Q0,
     *          0.1903505780647854099132564024210136828Q0,
     *          0.2094821410847278280129991748917142637Q0/
      DATA WK15/0.0229353220105292249637320080589695920Q0,
     *          0.1047900103222501838398763225415180174Q0,
     *          0.1690047266392679028265834265985502841Q0,
     *          0.2044329400752988924141619992346490847Q0/
      DATA AK15/0.9914553711208126392068546975263285166Q0,
     *          0.8648644233597690727897127886409262012Q0,
     *          0.5860872354676911302941448382587295984Q0,
     *          0.2077849550078984676006894037732449134Q0/

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

!     ------------------------------------

      FUNCTION F(X)
      IMPLICIT REAL*16(A-H,O-Z)
      COMMON/FN/ALP,XX,YY,NP
      EXTERNAL F1

!     THE REQUIRED INTEGRAND IS THE INTEGRAL OVER Y AND Z DIRECTIONS

!     LIMITS FOR INTEGRATION IN Y DIRECTION
      A=0
      B=SQRT(2.*(1.-X*X))
      REPS=1.Q-20
      AEPS=1.Q-21
      NMAX=9000

!     STORE THE VALUE OF X IN THE COMMON BLOCK FOR USE BY F1 AND F2
      XX=X

!     INTEGRATE W. R. T. Y
      CALL ADPIN1(F,A,B,REPS,AEPS,DIF,F1,IER,NPT,NMAX)
      IF(IER.GT.0) STOP 11
      END

!     ----------------------------------------------

      FUNCTION F1(Y)
      IMPLICIT REAL*16(A-H,O-Z)
      COMMON/FN/ALP,XX,YY,NP
      EXTERNAL F2

!     THE REQUIRED INTEGRAND IS THE INTEGRAL OVER Z DIRECTION

!     LIMITS FOR INTEGRATION IN Z DIRECTION

      A=0
      B=(1.-XX*XX-Y*Y/2.)
      F1=0.0
      IF(B.LE.0.0) RETURN
      B=SQRT(3.*B)
      REPS=1.Q-20
      AEPS=1.Q-21
      NMAX=20000

!     STORE THE VALUE OF Y IN THE COMMON BLOCK FOR USE BY F2
      YY=Y

!     INTEGRATE W.R.T. Z TO OBTAIN THE FUNCTION F1
      CALL ADPIN2(F1,A,B,REPS,AEPS,DIF,F2,IER,NPT,NMAX)
      if(ier.gt.0) print *,npt,dif,xx,y
      IF(IER.GT.0) STOP 12

!     ACCUMULATE THE NO. OF FUNCTION EVALUATIONS IN THE COMMON BLOCK
      NP=NP+NPT
      END

!     ------------------------------------

      FUNCTION F2(Z)
      IMPLICIT REAL*16(A-H,O-Z)
      COMMON/FN/ALP,XX,YY,NP

!     THE REQUIRED INTEGRAND
!     THE VALUES OF X AND Y ARE TAKEN FROM THE COMMON BLOCK
!     TAKING ACCOUNT OF SYMMETRY ALONG Y AND Z AXES THE INTEGRAND
!     IS MULTIPLIED BY 4 SINCE THE RANGE OF INTEGRATION IS
!     RESTRICTED TO POSITIVE QUADRANT IN Y-Z PLANE

      AF=(XX-ALP)**2+YY**2+Z**2
      F2=0.0
      IF(AF.GT.0.0) F2=4./SQRT(AF)
      END
