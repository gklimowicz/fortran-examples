!	To minimise a function in one dimension using Golden section search
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
!		IER=50 implies that function values at A, X, B are not
!			distinct, most likely due to roundoff errors
!		IER=51 implies that subroutine failed to reduce the
!			bracketing interval to required level
!		IER=522 implies that initial values of A,X,B do not bracket the minimum
!	F : (input) Name of the function routine to calculate the function
!		which is to be minimised
!
!	FUNCTION F(X) to calculate the required function, must be supplied
!		by the user.
!
!	Required routines : F

      SUBROUTINE GOLDEN(A,B,X,FX,REPS,AEPS,IER,F)
      IMPLICIT REAL*16(A-H,O-Z)
      PARAMETER(GR=0.381966Q0,NIT=100)

      FA=F(A)
      FB=F(B)
      FX=F(X)
      IF(FX.GT.MIN(FA,FB).OR.(X-A.GT.0.EQV.X-B.GT.0)) THEN
!	(A, X, B) does not bracket a minimum
        IER=522
        RETURN
      ENDIF
      IER=0

      DO 2000 I=1,NIT
!	Choose the next point in the larger section
        IF(ABS(X-A).GT.ABS(X-B)) THEN
          XP=X+GR*(A-X)
          FP=F(XP)
          IF(FP.LT.FX) THEN
!	(A, XP, X) bracket the minimum
            B=X
            X=XP
            FB=FX
            FX=FP
          ELSE
!	(XP, X, B) bracket the minimum
            A=XP
            FA=FP
          ENDIF
        ELSE
          XP=X+GR*(B-X)
          FP=F(XP)
          IF(FP.LT.FX) THEN
!	(X, XP, B) bracket the minimum
            A=X
            X=XP
            FA=FX
            FX=FP
          ELSE
!	(A, X, XP) bracket the minimum
            B=XP
            FB=FP
          ENDIF
        ENDIF

        IF(ABS(B-A).LT.MAX(ABS(X)*REPS,AEPS)) RETURN
        IF(FX.GE.FA.AND.FX.GE.FB) THEN
!	Roundoff errors may be dominating
          IER=50
          RETURN
        ENDIF
2000  CONTINUE

!	Fails to converge to required accuracy
      IER=51
      END
