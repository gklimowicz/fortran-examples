!     PROGRAM TO SOLVE FREDHOLM EQUATION USING COLLOCATION METHOD

      PROGRAM INTEQ
      IMPLICIT REAL*16(A-H,O,P,R-Z)
      IMPLICIT LOGICAL(Q)
      DIMENSION X(65),WK(65,65),INT(65),F(65)
      EXTERNAL FG,FKER,PHI,PSI

!     EXAMPLE 13.4 : FREDHOLM EQUATION OF THE FIRST KIND

51    FORMAT('   IER =',I4,5X,'N =',I3,5X,'COEF. =',1P3D14.6/
     1       (2X,5D14.6))
52    FORMAT('    X =',1PD14.6,5X,'F(X) =',2D14.6)

      REPS=1.Q-7
      AEPS=1.Q-8
      A=0.0
      B=1.0
      IQ=0
      IT=1

!     For collocation method the number of points is equal to the
!     number of basis functions.
100   PRINT *,'TYPE N=NO. OF PTS     (QUITS WHEN N.LE.0)'
      READ *,N
      IF(N.LE.0) STOP

!     CHOOSE UNIFORMLY SPACED POINTS AS THE COLLOCATION POINTS

      DO 1000 I=1,N
        X(I)=DFLOAT(I)/N
1000  CONTINUE

      CALL FREDCO(N,A,B,F,X,REPS,AEPS,WK,INT,IQ,IT,IER)
      WRITE(6,51) IER,N,(F(J),J=1,N)

!     CALCULATING THE SOLUTION AT SOME SELECTED POINTS

      DO 2000 I=1,5
        XI=(I-1)*0.25Q0
        FI=0.0
        DO 1500 J=1,N
1500    FI=FI+PHI(J,XI)*F(J)
        WRITE(6,52) XI,FI
2000  CONTINUE

      GO TO 100
      END

!     -------------------------------------------------

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

!     --------------------------------------------------

!	To solve linear Fredholm equation of first or second kind using
!	collocation method
!	It can be used to solve Volterra equations by defining the kernel
!	to be zero for t>x
!
!	N : (input) Number of points to be used in the collocation method
!	A : (input) Lower limit of the integral
!	B : (input) Upper limit of the integral
!	F : (output) Real array of length N containing the calculated coefficients
!		of expansion
!		solution = SUM_I F(I)*PHI(I,X)
!	X : (input) Real array of length N containing the points to
!		be used for collocation.
!	REPS : (input) Required relative accuracy to which integrals are to
!		be calculated
!	AEPS : (input) Required absolute accuracy to which integrals are to
!		be calculated
!		REPS and AEPS are passed on to subroutine ADPINT for calculating
!		the integrals when IQ=0. Otherwise these variables are not used.
!	WK : Real array of length N*N used as scratch space.
!	IWK : Integer array of length N used as a scratch space
!	IQ : (input) Integer variable to specify the treatment for integrals
!		PSI(I,X)=Integral[FKER(X,T)*PHI(I,T) dT] over [A,B]
!		If IQ=0, then the integrals are evaluated using subroutine
!			ADPINT, using function routine FUNK to calculate the
!			integrand, which in turn requires, PHI(I,T) and FKER(X,T). 
!		Otherwise the integrals are calculated using a user supplied
!			routine PSI(I,X).
!	IT : (input) Integer variable to specify the type of integral equation
!		If IT=1 Fredholm equation of the first kind is solved
!		If IT=2 Fredholm equation of the second kind is solved
!	IER : (output) The error parameter, IER=0 implies successful execution
!		IER=708 implies N<1, IT>2 or IT.LE.0, No calculations are done.
!		Other values of IER may be set by GAUELM and ADPINT
!
!	FUNCTION FG(X), FUNCTION PHI(I,X) and FUNCTION FKER(X,T) (for IQ=0)
!	or FUNCTION PSI(I,X) (for IQ.NE.0) must be supplied by the user.
!	Names of these function routines are fixed. FG(X) is the right
!	hand side function g(x), FKER(X,T) is the kernel, PHI(I,X) calculates
!	the basis functions phi_i(x), while PSI(I,X) calculates the integrals
!	defined above. The common block ZZFRED is used to pass on the variables
!	to FUNK. XI and II are the values of X and I.
!
!	Required routines : GAUELM, ADPINT, GAUS16 (or KRONRD), FUNK, FG, FKER, PHI, PSI
!	
      SUBROUTINE FREDCO(N,A,B,F,X,REPS,AEPS,WK,IWK,IQ,IT,IER)
      IMPLICIT REAL*16(A-H,O-Z)
      EXTERNAL FUNK
!	To pass arguments to FUNK
      COMMON/ZZFRED/XI,II
      DIMENSION F(N),X(N),WK(N,N),IWK(N)

      IER=0
      IF(N.LT.1.OR.IT.GT.2.OR.IT.LE.0) THEN
        IER=708
        RETURN
      ENDIF
      NMAX=10000

!	Setting up the system of linear equations
      DO 2000 J=1,N
        F(J)=FG(X(J))
        DO 2000 I=1,N
          IF(IQ.EQ.0) THEN
!	Evaluate the integrals numerically, split the range into two 
!	to tackle possible discontinuity at t=X
            XI=X(I)
            II=J
            CALL ADPINT(RI,A,XI,REPS,AEPS,DIF,FUNK,IER,NPT,NMAX)
            IF(IER.GT.100) RETURN
            CALL ADPINT(RI1,XI,B,REPS,AEPS,DIF1,FUNK,IER,NPT1,NMAX)
            WK(I,J)=RI+RI1
            IF(IER.GT.100) RETURN
          ELSE
!	Calculate the integrals PSI(I,X) using user supplied routine
            WK(I,J)=PSI(J,X(I))
          ENDIF
          IF(IT.EQ.2) WK(I,J)=WK(I,J)-PHI(J,X(I))
2000  CONTINUE

!	Solve the resulting system of linear equations
      NUM=1
      LJ=N
      IFLG=0
      CALL GAUELM(N,NUM,WK,F,DET,IWK,LJ,IER,IFLG)
      END

!     ---------------------------------------------------------

!	Solution of a system of linear equations using Gaussian elimination
!	with partial pivoting
!
!	N : (input) Number of equations to be solved
!	NUM : (input) Number of different sets (each with N equations) of
!	        equations to be solved
!	A : (input/output) The matrix of coefficient of size LJ*N
!	        A(I,J) is the coefficient of x_J in Ith equation
!	     	at output it will contain the triangular decomposition
!	X : (input/output) The matrix containing right hand sides (size LJ*NUM)
!	        X(I,J) is the Ith element of Jth right hand side
!	     	at output it will contain the solutions
!	DET : (output) The determinant of the matrix
!	INC : (output) Integer array of length N containing information about
!		interchanges performed during elimination
!	LJ : (input) First dimension of arrays A and X in calling program
!	IER : (output) Error flag, IER=0 signifies successful execution
!		IER=101 implies (N.LE.0 or N.GT.LJ) 
!		IER=121 implies some pivot turned out to be zero and hence
!			matrix must be nearly singular
!	IFLG : (input) Integer parameter to specify the type of computation required
!		If IFLG.LE.0, both elimination and solution are
!			done and IFLG is set to 2
!		If IFLG=1, only elimination is done and IFLG is set to 2
!		If IFLG.GE.2 only solution is calculated, the triangular
!			decomposition should have been calculated earlier
!
!	Required routines : None

      SUBROUTINE GAUELM(N,NUM,A,X,DET,INC,LJ,IER,IFLG)
      IMPLICIT REAL*16(A-H,O-Z)
!	For complex matrices use the following statements instead
!      IMPLICIT REAL*8(R)
!      IMPLICIT COMPLEX*16(A-H,S-Z)

      DIMENSION A(LJ,N),INC(N),X(LJ,NUM)

      IF(N.LE.0.OR.N.GT.LJ) THEN
        IER=101
        RETURN
      ENDIF

      IER=121
      IF(IFLG.LE.1) THEN
!	Perform elimination

        DET=1.0
        DO 2600 K=1,N-1
!	Find the maximum element in the Kth column
          R1=0.0
          KM=K
          DO 2200 L=K,N
            IF(ABS(A(L,K)).GT.R1) THEN
              R1=ABS(A(L,K))
              KM=L
            ENDIF
2200      CONTINUE

          INC(K)=KM
          IF(KM.NE.K) THEN
!	Interchange the rows if needed
            DO 2300 L=K,N
              T1=A(K,L)
              A(K,L)=A(KM,L)
2300        A(KM,L)=T1
            DET=-DET
          ENDIF

          DET=DET*A(K,K)
          IF(A(K,K).EQ.0.0) RETURN
!	To check for singular or nearly singular matrices replace this
!	statement by, where REPS is approximately \hcross*Max(A(I,J))
!         IF(ABS(A(K,K)).LT.REPS) RETURN
          DO 2500 L=K+1,N
            A(L,K)=A(L,K)/A(K,K)
            DO 2500 L1=K+1,N
2500      A(L,L1)=A(L,L1)-A(L,K)*A(K,L1)
2600    CONTINUE
        DET=DET*A(N,N)
        INC(N)=N
!	If pivot is zero then return, IER has been set to 121
        IF(A(N,N).EQ.0.0) RETURN
!	To check for singular or nearly singular matrices replace this
!	statement by, where REPS is approximately \hcross*Max(A(I,J))
!         IF(ABS(A(N,N)).LT.REPS) RETURN

        IER=0
        IF(IFLG.EQ.1) THEN
          IFLG=2
          RETURN
        ENDIF
        IFLG=2
      ENDIF

      IER=0
!	Solution for the NUM different right-hand sides
      DO 5000 J=1,NUM
        DO 3000 K=1,N-1
!	Forward substitution
          IF(K.NE.INC(K)) THEN
            T1=X(K,J)
            X(K,J)=X(INC(K),J)
            X(INC(K),J)=T1
          ENDIF
          DO 3000 L=K+1,N
3000    X(L,J)=X(L,J)-A(L,K)*X(K,J)

!	back-substitution
        X(N,J)=X(N,J)/A(N,N)
        DO 3300 K=N-1,1,-1
          DO 3200 L=N,K+1,-1
3200      X(K,J)=X(K,J)-X(L,J)*A(K,L)
3300    X(K,J)=X(K,J)/A(K,K)
5000  CONTINUE
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

!     ---------------------------------------------------

!	Function routine to calculate the integrand for calculating
!	PSI(I,X) as required by subroutine FREDCO.
!
!	Function FKER(X,T) is the kernel K(x,t) and PHI(I,T) is the Ith
!	basis function, phi_i(t). The argument X and I are passed through
!	the common block.
!
!	FUNCTION FKER(X,T) and FUNCTION PHI(I,T) must be supplied by the user
!
!	Required routines : FKER, PHI

      FUNCTION FUNK(T)
      IMPLICIT REAL*16(A-H,O-Z)
!	To pass parameters from subroutine FREDCO
      COMMON/ZZFRED/X,I

      FUNK=FKER(X,T)*PHI(I,T)
      END

!     ------------------------------------------

      FUNCTION FG(X)
      IMPLICIT REAL*16(A-H,O-Z)
      PARAMETER(PI=3.1415926535Q0)

!     RHS FUNCTION

      FG=(EXP(1.+X)-1.)/(X+1)
      END

!     ---------------------------------------

      FUNCTION FKER(X,T)
      IMPLICIT REAL*16(A-H,O-Z)

!     THE KERNEL

      FKER=EXP(X*T)
      END

!     ------------------------------

      FUNCTION PSI(I,X)
      IMPLICIT REAL*16(A-H,O-Z)

!     DUMMY FUNCTION, SINCE INTEGRALS ARE EVALUATED NUMERICALLY

      PSI=0.0
      END

!     -------------------------------------

      FUNCTION PHI(I,X)
      IMPLICIT REAL*16(A-H,O-Z)

!     THE BASIS FUNCTIONS

      PHI=0.0
      IF(X.EQ.0.0) RETURN
      PHI=X**(I-1)
      END
