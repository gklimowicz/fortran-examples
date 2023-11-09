!	To bracket a minimum in one dimension
!
!	A,B : (input/output) Starting values for the search. After
!		execution they will contain the limits of bracketing triplet
!	X : (output) Real variable which will contain the point X between
!		A and B, such that F(X) < MIN(F(A),F(B))
!	F : (input) Name of the function routine to calculate the function
!		which is to be minimised
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=501 implies that A=B, in which case no calculations are done
!		IER=521 implies that subroutine failed to find a bracketing triplet
!
!	FUNCTION F(X) to calculate the required function, must be supplied
!		by the user.
!
!	Required routines : F

      SUBROUTINE BRACKM(A,B,X,F,IER)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(GR=1.618034D0,RMAX=10.,NTRY=200)

      IF(A.EQ.B) THEN
        IER=501
        RETURN
      ENDIF
      IER=0

      FA=F(A)
      FB=F(B)
      IF(FB.GT.FA) THEN
!	Exchange A and B
        T=A
        A=B
        B=T
        T=FA
        FA=FB
        FB=T
      ENDIF
!	First guess for X
      X=B+GR*(B-A)
      FX=F(X)

      DO 2000 I=1,NTRY
        IF(FB.LT.FX) THEN
!	The minimum is bracketed by the triplet (A,B,X), exchange B and X
          T=B
          B=X
          X=T
          RETURN
        ENDIF

!	Next guess by parabolic interpolation
        R=(B-A)*(FB-FX)
        T=(B-X)*(FB-FA)
        P=(B-X)*T-(B-A)*R
        T=2.*(T-R)
        IF(T.GT.0) P=-P
        T=ABS(T)

        IF(T.GT.0.AND.(P.GT.0.EQV.(X-B).GT.0)) THEN
          IF(ABS(P).LT.ABS(T*(X-B))) THEN
!	Interpolated point is between B and X
            XP=B+P/T
            FP=F(XP)
            IF(FP.LT.FX) THEN
!	Minimum is bracketed by (B,XP,X)
              A=B
              B=X
              X=XP
              RETURN
            ELSE IF(FP.GT.FB) THEN
!	Minimum is bracketed by (A,B,XP)
              X=B
              B=XP
              RETURN
            ENDIF
!	If the minimum is not bracketed, then reject the interpolated point
            XP=X+GR*(X-B)
            FP=F(XP)

          ELSE IF(ABS(P).LT.ABS(T*(X-B)*RMAX)) THEN
!	Interpolated point is between X and the allowed limit
            XP=B+P/T
            FP=F(XP)
            IF(FP.LT.FX) THEN
!	If the minimum is not bracketed, then reject the interpolated point
              XP=X+GR*(X-B)
              FP=F(XP)
            ENDIF
          ELSE
!	Reject the interpolated point if it is beyond the allowed limit
            XP=B+RMAX*(X-B)
            FP=F(XP)
          ENDIF
        ELSE
!	Reject the interpolated point if it is on the other side of B
          XP=X+GR*(X-B)
          FP=F(XP)
        ENDIF

        A=B
        B=X
        X=XP
        FA=FB
        FB=FX
        FX=FP
2000  CONTINUE

!	Fails to bracket the minimum
      IER=521
      END
