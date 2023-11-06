!	Real zero of a given function using Brent's method
!
!	A : (input/output) Lower limit of interval where zero is located
!	B : (input/output) Upper limit of interval where zero is located
!		The limits A, B are updated by the subroutine
!	X : (output) Computed value of the zero
!	REPS : (input) Required relative accuracy
!	AEPS : (input) Required absolute accuracy
!		The estimated error should be less than MAX(AEPS, REPS*ABS(X))
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=427 implies that function has same sign at x=A and B
!			in which case no calculations are done
!		IER=428 implies that iteration failed to converge to specified
!			accuracy
!	F : (input) Name of the function routine to calculate the function
!		FUNCTION F(X) must be supplied by the user.
!
!	Required routines : F

      SUBROUTINE BRENT(A,B,X,REPS,AEPS,IER,F)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NIT=75)

      IER=0
      FA=F(A)
      FB=F(B)
      IF(FA.GT.0.0.EQV.FB.GT.0.0) THEN
!	If the function has the same sign at both the end points, then quit
        IER=427
        RETURN
      ENDIF

      C=A
      FC=FA
      D=B-A
      E=D

!	The iteration loop
      DO 2000 I=1,NIT
        IF(ABS(FC).LT.ABS(FB)) THEN
!	If F(C) < F(B), then interchange B and C
          A=B
          B=C
          C=A
          FA=FB
          FB=FC
          FC=FA
        ENDIF

        EPS=MAX(REPS*ABS(B),AEPS)
        DX=0.5*(C-B)
        IF(ABS(DX).LT.EPS.OR.FB.EQ.0) THEN
!	Iteration has converged
          IER=0
          X=B
          RETURN
        ENDIF

        IF(ABS(E).LT.EPS.OR.ABS(FA).LE.ABS(FB)) THEN
!	If the previous step is too small or if F(A).LE.F(B), perform bisection
          D=DX
          E=DX
        ELSE
          R3=FB/FA
          IF(A.EQ.C) THEN
!	Try linear interpolation
            P=2.*DX*R3
            P1=1.-R3
          ELSE
!	Try inverse quadratic interpolation
            R1=FA/FC
            R2=FB/FC
            P=R3*(2.*DX*R1*(R1-R2)-(B-A)*(R2-1))
            P1=(R1-1.)*(R2-1.)*(R3-1.)
          ENDIF
          IF(P.GT.0) THEN
            P1=-P1
          ELSE
            P=-P
          ENDIF

          E=D
          IF(2.*P.LT.3.*DX*P1-ABS(EPS*P1).AND.P.LT.ABS(0.5*E*P1)) THEN
!	Accept the interpolated value
            D=P/P1
          ELSE
!	otherwise, perform bisection
            D=DX
            E=DX
          ENDIF
        ENDIF

        A=B
        FA=FB
        IF(ABS(D).GT.EPS) THEN
          B=B+D
        ELSE
!	If the change is too small, shift by EPS
          B=B+SIGN(EPS,DX)
        ENDIF
        FB=F(B)
        IF(FB.GT.0.0.EQV.FC.GT.0.0) THEN
!	If F(B) and F(C) have the same sign, then replace C by A
          C=A
          FC=FA
          D=B-A
          E=D
        ENDIF
2000  CONTINUE

!	Iteration fails to converge
      IER=428
      X=B
      END
