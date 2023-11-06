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
