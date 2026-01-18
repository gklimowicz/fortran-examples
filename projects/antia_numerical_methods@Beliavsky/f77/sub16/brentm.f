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
      IMPLICIT REAL*16(A-H,O-Z)
      PARAMETER(NIT=75,GR=0.381966Q0)

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
