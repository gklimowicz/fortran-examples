!     PROGRAM TO OBTAIN RATIONAL FUNCTION MINIMAX APPROXIMATION FOR
!     MATHEMATICAL FUNCTIONS USING REMES ALGORITHM

      PROGRAM MINMAX
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(20),WK(400),IWK(20),F(200),E(25),X(200)

!     EXAMPLE 10.10: RELATIVE MINIMAX APPROXIMATION TO ARC TAN(X)

51    FORMAT('   IER =',I4,5X,'MAX. ERROR =',1PD14.6/
     1       '   COEF. IN NUMERATOR',4D16.8/(2X,5D16.8))
52    FORMAT('   COEF. IN DENOMINATOR',1P4D16.8/(2X,5D16.8))
53    FORMAT('  EXTREMA IN ERROR CURVE',1P4D14.6/(2X,5D14.6))
54    FORMAT('   M =',I3,5X,'K =',I3,5X,'N =',I3,5X,'IFLG =',I3)
55    FORMAT('  INITIAL GUESS FOR COEFFICIENTS :',1P2D14.6/
     1       (2X,1P5D14.6))
56    FORMAT('  INITIAL GUESS FOR EXTREMA OF ERROR CURVE :',1P2D14.6/
     1       (2X,1P5D14.6))

100   PRINT *,'TYPE M=DEGREE OF NUMERATOR,  K=DEGREE OF DENOMINATOR'
      PRINT *,'     (QUITS WHEN (M+K).LE.0)'
      READ *,M,K
      IF(M+K.LE.0) STOP
      PRINT *,'TYPE N=NO. OF PTS FOR INITIAL SEARCH OF EXTREMA'
      READ *,N
      PRINT *,'TYPE XL=LOWER LIMIT,  XU=UPPER LIMIT, IFLG=0/1/2'
      READ *,XL,XU,IFLG
      MM=M
      KK=K
      IF(XL.EQ.XU) STOP

      EPS=1.D-8
      EPSM=1.D-5
      IE=M+K+2
      WRITE(6,54) M,K,N,IFLG
      IF(IFLG.EQ.0) THEN
        PRINT *,'TYPE INITIAL GUESS FOR COEF'
        READ *,(A(I),I=1,M+K+1)
        WRITE(6,55)(A(I),I=1,M+K+1)
      ELSE IF(IFLG.EQ.2) THEN
        PRINT *,'TYPE INITIAL GUESS FOR EXTREMA'
        READ *,(E(I),I=1,M+K+2)
        WRITE(6,56) (E(I),I=1,M+K+2)
      ENDIF

      CALL REMES(M,K,N,XL,XU,A,X,F,E,IE,EMAX,EPS,EPSM,IFLG,IER,WK,IWK)
      WRITE(6,51) IER,EMAX,(A(I),I=K+1,M+K+1)
      WRITE(6,52) (A(I),I=1,K)
      WRITE(6,53) (E(I),I=1,IE)
      GO TO 100
      END

!     -------------------------------------------

!	To minimise a function in one dimension using Brent's method
!
!	A,B,X : (input/output) Triplet which brackets the minimum.
!		After execution these will contain the final triplet
!		with X giving the best estimate for minimiser.
!		X must be between A and B; and F(X)<MIN(F(A),F(B))
!	FX : (output) The function value at X.
!	REPS : (input) Required relative accuracy
!	AEPS : (input) Required absolute accuracy
!		The bracketing interval is subdivided until
!		ABS(B-A) < MAX(AEPS, REPS*ABS(X))
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=51 implies that subroutine failed to reduce the
!			bracketing interval to required level
!		IER=523 implies that initial values of A,X,B do not bracket
!			the minimum
!	F : (input) Name of the function routine to calculate the function
!		which is to be minimised
!
!	FUNCTION F(X) to calculate the required function, must be supplied
!		by the user.
!
!	Required routines : F

      SUBROUTINE BRENTM(A,B,X,FX,REPS,AEPS,IER,F)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NIT=75,GR=0.381966D0)

      IER=0
      IF(A.GT.B) THEN
!	Interchange A and B
        T=A
        A=B
        B=T
      ENDIF
      FA=F(A)
      FB=F(B)
      FX=F(X)
      IF(FA.LT.FX.OR.FB.LT.FX.OR.X.LE.A.OR.X.GE.B) THEN
        IER=523
        RETURN
      ENDIF

      E=0
      V=X
      W=X
      FV=FX
      FW=FX

!	Loop for iteration
      DO 2000 I=1,NIT
        XM=0.5*(A+B)
        EPS2=MAX(REPS*ABS(X),AEPS)
        EPS=0.5*EPS2

!	The convergence test
        IF(ABS(X-XM).LT.EPS2-0.5*(B-A)) THEN
          IER=0
          RETURN
        ENDIF

        P=0
        T=0
        R=0
        IF(ABS(E).GT.EPS) THEN
!	Parabolic interpolation
          R=(X-W)*(FX-FV)
          T=(X-V)*(FX-FW)
          P=(X-V)*T-(X-W)*R
          T=2*(T-R)
          IF(T.GT.0) P=-P
          T=ABS(T)
          R=E
          E=D
        ENDIF

        IF(ABS(P).LT.ABS(.5*T*R).AND.P.GT.T*(A-X).AND.P.LT.T*(B-X)) THEN
!	accept the interpolated point
          D=P/T
          U=X+D
          IF(U-A.LT.EPS2.OR.B-U.LT.EPS2) THEN
!	If it is too close to end points shift it by EPS at least
            D=EPS
            IF(X.GE.XM) D=-EPS
          ENDIF
        ELSE
!	Perform golden section search
          E=B-X
          IF(X.GE.XM) E=A-X
          D=GR*E
        ENDIF
        IF(ABS(D).GE.EPS) THEN
          U=X+D
        ELSE
!	Shift the point by at least EPS
          U=X+SIGN(EPS,D)
        ENDIF
        FU=F(U)

!	Updating the bracketing triplet
        IF(FU.LE.FX) THEN
          IF(U.LT.X) THEN
!	(A, U, X) is the triplet
            B=X
          ELSE
!	(X, U, B) is the triplet
            A=X
          ENDIF
          V=W
          FV=FW
          W=X
          FW=FX
          X=U
          FX=FU
        ELSE
          IF(U.LT.X) THEN
!	(U, X, B) is the triplet
            A=U
          ELSE
!	(A, X, U) is the triplet
            B=U
          ENDIF
          IF(FU.LE.FW.OR.W.EQ.X) THEN
            V=W
            FV=FW
            W=U
            FW=FU
          ELSE IF(FU.LE.FV.OR.V.EQ.X.OR.V.EQ.W) THEN
            V=U
            FV=FU
          ENDIF
        ENDIF

2000  CONTINUE

!	Iteration fails to converge
      IER=51
      END

!     --------------------------------------------------

!	To calculate the error in rational function approximation
!	For use with subroutine REMES. It is called by BRENTM to find
!	extrema of error curve. It may not be used for any other purpose.
!
!	X : (input) the value at which error is to be calculated
!	FM = (FUN(X) - FUND(X)*R_mk(X))*SI  is the calculated error
!	The parameters for rational function are passed through common block
!	M,K are degree of numerator and denominator
!	A is a real array containing the coefficient of rational function
!		approximation
!	SI is -1 for maximum and +1 for minimum. The error is multiplied
!		by SI to make it a minimum.
!		For initial scan SI>10 and function is not evaluated.
!
!	Functions FUN(x) and FUND(x) must be provided by the user to
!		seek approximation of form FUN(X) = FUND(X)*RMK(X)
!
!	Required routines : FUN, FUND
!
      FUNCTION FM(X)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NMAX=50)
      COMMON/ZZFN/A(NMAX),SI,M,K

      NK=M+K+2
      FM=0.0
      IF(SI.LT.10) FM=FUN(X)

!	Calculate the numerator using nested multiplication
      FN=A(K+M+1)
      DO 2200 J=1,M
2200  FN=FN*X+A(NK-J-1)

      IF(K.GT.0) THEN
!	Calculate the denominator using nested multiplication
        FD=A(K)
        DO 2300 J=1,K-1
2300    FD=FD*X+A(K-J)
        FD=FD*X+1
      ELSE
        FD=1.
      ENDIF

      FM=(FM-FUND(X)*FN/FD)*SIGN(1.D0,SI)
      END

!     -------------------------------------------

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
      IMPLICIT REAL*8(A-H,O-Z)
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

!     ----------------------------------

!	To calculate coefficients of Rational function minimax approximation
!	for a given function over a finite interval using Remes algorithm
!
!	M : (input) Required degree of polynomial in the numerator
!	K : (input) Required degree of polynomial in the denominator
!	N : (input/output) Number of points to be used for scanning the
!		extrema of error curves. If N<3(M+K+2) it will be set to
!		this value.
!	XL : (input) Lower limit of the interval over which approximation is required
!	XU : (input) Upper limit of the interval over which approximation is required
!	A : (input/output) Real array of length M+K+1 containing the coefficients
!		of approximation. A(I) is the coefficient of x**I in
!		the denominator, the constant term being 1. A(K+J+1) is
!		the coefficient of x**J in the numerator. If IFLG=0 the
!		initial guess to coefficients must be supplied.
!	X : (output) Real array of length N containing the points at which
!		function value is calculated for scanning the extrema in
!		error curve
!	F : (output) Real array of length N containing the function values at X(I)
!	EX : (input/output) Real array of length M+K+5 containing the
!		extrema of error curve. If IFLG=2, then the initial guess
!		for these extrema must be supplied
!	IE : (input/output) Number of extrema in error curve. This value
!		must be supplied if IFLG=2
!	EMAX : (output) Maximum error in the calculated approximation
!	EPS : (input) Required accuracy, the iteration for calculating the
!		coefficients of approximations is continued until
!		difference in maximum error is less than EPS
!	EPSM : (input) Required accuracy to which extrema in error curve
!		are determined.
!	IFLG : (input) Integer parameter specifying the nature of initial
!		approximation that is supplied. 
!		If IFLG=0 then iteration is started from initial guess for
!			coefficients supplied in array A
!		If IFLG=1 then no initial guess is required. Iteration is
!			started by assuming that the extrema of error curve
!			coincides with those of T_{M+K+1}(x).
!		If IFLG=2 then iteration is started from initial guess for
!			position of extrema in error curve. These must be
!			supplied in array EX and IE should be set to the number
!			of extrema
!	IER : (output) error parameter, IER=0 implies successful execution
!		IER=614 implies that M<0, K<0, XU.LE.XL or M+K+2>NMAX (=50)
!			in this case no calculations are done
!		IER=632 implies that the Remes iteration did not converge
!			to the specified accuracy
!		IER=633 implies that at some stage the error curve does not
!			have the required number of extrema.
!		Other values of IER may be set by GAUELM or BRENTM
!	WK : Real array of length (K+M+2)**2 used as scratch space
!	IWK : Integer array of length K+M+2 used as scratch space
!
!	Functions FUN(X) and FUND(X) must be supplied by the user.
!	Here X, FUN and FUND are real variables and we seek approximation
!	of the form FUN(x) ~ FUND(X)*R_mk(x)
!	To obtain minimax approximation with respect to absolute error
!	set FUND(X)=1 and FUN(X) to required function
!	To obtain minimax approximation with respect to relative error
!	set FUN(X)=1 and FUND(X) to reciprocal of required function.
!
!	Required routines : GAUELM, BRENTM, FM, FUN, FUND
!
      SUBROUTINE REMES(M,K,N,XL,XU,A,X,F,EX,IE,EMAX,EPS,EPSM,IFLG,IER,
     1                 WK,IWK)
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL FM
      PARAMETER(NIT=30,NJT=10,NMAX=50,PI=3.14159265358979324D0)
      COMMON/ZZFN/AA(NMAX),SI,MM,KK
      DIMENSION A(M+K+2),F(N),IWK(K+M+2),WK(K+M+2,*),EX(*),X(N)

      NK=M+K+2
      IF(M.LT.0.OR.K.LT.0.OR.NK.GT.NMAX.OR.XU.LE.XL) THEN
        IER=614
        RETURN
      ENDIF

!	Copy the value in the common block for use by function FM
      MM=M
      KK=K
      IF(N.LT.3*NK) N=3*NK

!	Generating the mesh points for crude scan of extrema
      H=(XU-XL)/(N-1.)
      DO 2000 I=1,N
        X(I)=XL+(I-1)*H
        F(I)=FUN(X(I))
2000  CONTINUE

      IER=0
      NUM=1
      LJ=NK
      REPS=0.0

!	Loop for Remes iteration
      DO 5000 IT=1,NIT
        E1=0.0
!	Copy the coefficients to common block for function FM
        DO 2200 I=1,M+K+1
2200    AA(I)=A(I)

        IF(IT.EQ.1.AND.IFLG.GE.1) THEN
          IF(IFLG.EQ.1) THEN
!	Use the extrema of Chebyshev polynomial T_{M+K+1} as initial guess
            EX(1)=X(1)
            IE=1
            DO 2300 I=2,NK-1
2300        EX(I)=0.5*(XL+XU)+0.5*(XU-XL)*COS((NK-I)*PI/(NK-1.))
            EX(NK)=XU
          ENDIF
          EI=0.0
          EMAX=0.0
          EMIN=0.0
        ELSE

!	Locate the extrema of the error curve
          DO 2400 I=1,N
!	Flag to avoid calculating the function repeatedly
            SI=22.0
            EI=F(I)+FM(X(I))
            IF(I.GT.2) THEN
              IF(E1.GE.E2.EQV.E1.GE.EI) THEN
!	An extrema is bracketed
                SI=1.
!	To convert maxima to minima
                IF(E1.GT.E2.OR.E1.GT.EI) SI=-1.
                AI=X(I-2)
                B=X(I)
                XI=X(I-1)
                CALL BRENTM(AI,B,XI,FX,REPS,EPSM,IER,FM)
                IF(IER.GT.0) RETURN

                IE=IE+1
                EX(IE)=XI
                WK(IE,1)=FX*SI
                IF(ABS(FX).GT.EMAX) EMAX=ABS(FX)
                IF(ABS(FX).LT.EMIN) EMIN=ABS(FX)
!	To ensure that dimensions of EX are not exceeded
                IF(IE.GT.M+K+4) GO TO 2500
              ENDIF
            ELSE IF(I.EQ.1) THEN
!	The end point is always included in the list of extrema
              IE=1
              EX(1)=X(1)
              EMAX=ABS(EI)
              EMIN=EMAX
              WK(I,1)=EI
              E0=EI
            ENDIF
            E2=E1
            E1=EI
2400      CONTINUE

          IE=IE+1
!	The end point is always included in the list of extrema
          EX(IE)=X(N)
          WK(IE,1)=EI
          EMAX=MAX(EMAX,ABS(EI))
          EMIN=MIN(EMIN,ABS(EI))

          IF(IE.LT.NK) THEN
!	If the number of extrema is less than required then quit
            IER=633
            RETURN
          ENDIF
2500      IF(IE.GT.NK) THEN

!	Remove one extrema from the list
            IE=IE-1
            IF(ABS(WK(IE+1,1)).GT.ABS(WK(1,1))) THEN
!	remove the first extrema
              EMIN=ABS(WK(2,1))
              E0=WK(2,1)
              DO 2600 I=1,IE
                EX(I)=EX(I+1)
                IF(ABS(WK(I+1,1)).LT.EMIN) EMIN=ABS(WK(I+1,1))
2600          WK(I,1)=WK(I+1,1)
            ELSE
!	remove the last extrema
              EMIN=ABS(WK(1,1))
              DO 2700 I=2,IE
                IF(ABS(WK(I,1)).LT.EMIN) EMIN=ABS(WK(I,1))
2700          CONTINUE
            ENDIF
!	Repeat until the number of extrema = NK
            GO TO 2500
          ENDIF

!	Convergence check, quit if difference between various extrema is
!	less than 1%
          IF(EMAX-EMIN.LT.1.D-2*EMAX) RETURN
          EI=MIN(0.5*(EMAX+EMIN),2.*EMIN)
          IF(E0.LT.0.0) EI=-EI
        ENDIF

!	Iteration to determine the coefficients
        DO 4000 JT=1,NJT
          AI=1
!	Setting up the equation matrix
          DO 3600 I=1,NK
            FI=FUN(EX(I))
            FD=FUND(EX(I))
            DO 3000 J=1,K
3000        WK(I,J)=-(FI-AI*EI)*EX(I)**J
            WK(I,K+1)=FD
            DO 3200 J=1,M
3200        WK(I,J+K+1)=FD*EX(I)**J
            WK(I,NK)=AI
            A(I)=FI
            AI=-AI
3600      CONTINUE

          IFL=0
          CALL GAUELM(NK,NUM,WK,A,DET,IWK,LJ,IER,IFL)
          IF(IER.GT.0) RETURN
          DIF=ABS(A(NK)-EI)
          EI=A(NK)
!	convergence check
          IF(DIF.LT.EPS.OR.DIF.LT.0.3D0*(EMAX-EMIN)) GO TO 5000
          IF(IT.EQ.1.AND.IFLG.GE.1) EMAX=0.3D0*ABS(A(NK))
4000    CONTINUE

!	Even if iteration for calculating coefficients fails continue further

5000  CONTINUE

!	Remes iteration fails to converge
      IER=632
      END

!     ------------------------------------

      FUNCTION FUND(X)
      IMPLICIT REAL*8(A-H,O-Z)
!     SET FUND=1/F(X) FOR RELATIVE MINIMAX APPROXIMATION

      IF(X.EQ.0.0) THEN
!     TAKE THE LIMITING VALUE
      FUND=1.
      ELSE
      X1=SQRT(X)
      FUND=X1/ATAN(X1)
      ENDIF
      END

!     ------------------------------------

      FUNCTION FUN(X)
      IMPLICIT REAL*8(A-H,O-Z)

!     SET FUN=1 FOR RELATIVE MINIMAX APPROXIMATION 
      FUN=1.0
      END
