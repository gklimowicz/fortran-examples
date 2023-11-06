!     PROGRAM FOR RECURSIVE SOLUTION OF A SYSTEM OF TWO NONLINEAR EQUATIONS
!     IT SCANS FOR ZEROS ALONG THE FIRST VARIABLE BY LOOKING FOR SIGN
!     CHANGES AT SPECIFIED SEARCH STEP. THEN IT USED SECANT ITERATION OR
!     BRENT'S METHOD FOR FINDING THE ROOTS.
!     BRENT1 AND SECAN1 ARE COPIES OF BRENT AND SECANT RESPECTIVELY
!     THIS IS REQUIRED TO AVOID CALLING THE FUNCTIONS RECURSIVELY

      PROGRAM ROOT
      IMPLICIT REAL*16(A-H,O-Z)
      EXTERNAL F
      COMMON/FN/XX,YY,FY,IT

!     EXERCISE 7.49 (I)

51    FORMAT('   SEARCH STEP =',1PD12.4,3X,'FOR ROOTS IN (',D12.4,
     1       ',',D12.4,')    IT=',I2)
52    FORMAT('  IER =',I4,3X,'ROOT IN (',1PD12.4,',',D12.4,') ='
     1      ,2E14.6/5X,'F =',2E12.4)
53    FORMAT('   IER =',I4,3X,'ROOT =',1P2D14.6,3X,'F =',2D12.4)

      REPS=1.Q-30
      AEPS=1.Q-33

100   PRINT *,'TYPE DX=SEARCH STEP, XL=LOWER LIMIT, XU=UPPER LIMIT'
      PRINT *,'           (QUITS WHEN DX(XL-XU).GT.0)'
      READ *,DX,XL,XU
      IF(DX.GT.0.0.EQV.XL-XU.GT.0.0) STOP
      PRINT *,'TYPE IT=1/2   FOR   SECANT/BRENT'
      READ *,IT
      WRITE(6,51) DX,XL,XU,IT

!     LOCATE ZEROS BY LOOKING FOR SIGN CHANGES

      F1=F(XL)
      DO 1000 X1=XL+DX,XU,DX
        FX=F(X1)
        IF(.NOT.(FX.GT.0.0.EQV.F1.GT.0.0)) THEN

!     FIND ACCURATE VALUE OF THE ROOT

          A=X1-DX
          B=X1
          X0=(A+B)*0.5
          IF(IT.EQ.1) THEN
            CALL SECANT(X0,A,B,X,REPS,AEPS,IER,F)
            WRITE(6,52) IER,A,B,X,YY,F(X),FY
          ELSE
            CALL BRENT(A,B,X,REPS,AEPS,IER,F)
            WRITE(6,53) IER,X,YY,F(X),FY
          ENDIF
        ENDIF
        F1=FX
1000  CONTINUE
      GO TO 100
      END

!     -------------------------------------------

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
      IMPLICIT REAL*16(A-H,O-Z)
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

!     -------------------------------------------

!	Real zero of a given function using secant iteration
!
!	X0 : (input) Initial guess for the zero
!	XL : (input) Lower limit of interval where zero is expected
!	XU : (input) Upper limit of interval where zero is expected
!	X : (output) Computed value of the zero
!	REPS : (input) Required relative accuracy
!	AEPS : (input) Required absolute accuracy
!		The estimated error should be less than MAX(AEPS, REPS*ABS(X))
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=40 implies that function value is equal at two points
!			and it is not possible to continue the iteration
!		IER=402 implies XL>X0 or XU<X0, in which case no calculations are done
!		IER=422 implies that iteration goes outside the specified limits
!		IER=423 implies that iteration failed to converge to specified accuracy
!	FUN : (input) Name of the function routine to calculate the function
!		FUNCTION FUN(X) must be supplied by the user.
!
!	Required routines : FUN

      SUBROUTINE SECANT(X0,XL,XU,X,REPS,AEPS,IER,FUN)
      IMPLICIT REAL*16(A-H,O-Z)
      PARAMETER(NIS=75)

      IER=0
      IF(XL.GT.X0.OR.XU.LT.X0) THEN
!	If X0 is outside the specified interval (XL,XU) then quit
        IER=402
        RETURN
      ENDIF

      X=X0
!	Select the increment for the next point X+DX
      DX=(XU-X0)/100.
      IF(ABS(DX).LT.5.*MAX(REPS*ABS(X),AEPS)) DX=(XL-X0)/5.
      IF(ABS(DX).GT.0.1Q0*MAX(ABS(X),100.*AEPS))
     1    DX=SIGN(0.1Q0*MAX(ABS(X),100.*AEPS),DX)
      F1=0.0

      DO 1000 L=1,NIS
        F=FUN(X)
        DX1=DX

        IF(F1-F.EQ.0.0) THEN
          IF(F.EQ.0.0) RETURN
!	If F1=F and F.NE.0, then quit
          IER=40
          RETURN
        ENDIF

!	The secant iteration
        IF(L.GT.1) DX=DX1*F/(F1-F)
        X=X+DX
        F1=F

        IF(ABS(DX).LT.MAX(REPS*ABS(X),AEPS).AND.L.GT.2) RETURN
        IF(X.LT.XL.OR.X.GT.XU) THEN
!	If iteration goes outside the specified limits (XL,XU), then quit
          IER=422
          RETURN
        ENDIF
1000  CONTINUE

!	The iteration fails to converge
      IER=423
      END

!     -------------------------------------------

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

      SUBROUTINE BRENT1(A,B,X,REPS,AEPS,IER,F)
      IMPLICIT REAL*16(A-H,O-Z)
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

!     -------------------------------------------

!	Real zero of a given function using secant iteration
!
!	X0 : (input) Initial guess for the zero
!	XL : (input) Lower limit of interval where zero is expected
!	XU : (input) Upper limit of interval where zero is expected
!	X : (output) Computed value of the zero
!	REPS : (input) Required relative accuracy
!	AEPS : (input) Required absolute accuracy
!		The estimated error should be less than MAX(AEPS, REPS*ABS(X))
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=40 implies that function value is equal at two points
!			and it is not possible to continue the iteration
!		IER=402 implies XL>X0 or XU<X0, in which case no calculations are done
!		IER=422 implies that iteration goes outside the specified limits
!		IER=423 implies that iteration failed to converge to specified accuracy
!	FUN : (input) Name of the function routine to calculate the function
!		FUNCTION FUN(X) must be supplied by the user.
!
!	Required routines : FUN

      SUBROUTINE SECAN1(X0,XL,XU,X,REPS,AEPS,IER,FUN)
      IMPLICIT REAL*16(A-H,O-Z)
      PARAMETER(NIS=75)

      IER=0
      IF(XL.GT.X0.OR.XU.LT.X0) THEN
!	If X0 is outside the specified interval (XL,XU) then quit
        IER=402
        RETURN
      ENDIF

      X=X0
!	Select the increment for the next point X+DX
      DX=(XU-X0)/100.
      IF(ABS(DX).LT.5.*MAX(REPS*ABS(X),AEPS)) DX=(XL-X0)/5.
      IF(ABS(DX).GT.0.1Q0*MAX(ABS(X),100.*AEPS))
     1    DX=SIGN(0.1Q0*MAX(ABS(X),100.*AEPS),DX)
      F1=0.0

      DO 1000 L=1,NIS
        F=FUN(X)
        DX1=DX

        IF(F1-F.EQ.0.0) THEN
          IF(F.EQ.0.0) RETURN
!	If F1=F and F.NE.0, then quit
          IER=40
          RETURN
        ENDIF

!	The secant iteration
        IF(L.GT.1) DX=DX1*F/(F1-F)
        X=X+DX
        F1=F

        IF(ABS(DX).LT.MAX(REPS*ABS(X),AEPS).AND.L.GT.2) RETURN
        IF(X.LT.XL.OR.X.GT.XU) THEN
!	If iteration goes outside the specified limits (XL,XU), then quit
          IER=422
          RETURN
        ENDIF
1000  CONTINUE

!	The iteration fails to converge
      IER=423
      END

!     ----------------------------------------------------

      FUNCTION F(X)
      IMPLICIT REAL*16(A-H,O-Z)
      EXTERNAL FY
      COMMON/FN/XX,YY,FF,IT

      REPS=1.Q-7
      AEPS=1.Q-8

!     LIMITS ON Y

      A=-1.57Q0
      B=1.57Q0
      X0=1

!     STORE THE VALUE OF X IN THE COMMON BLOCK FOR USE BY FY
      XX=X

!     FOR A GIVEN VALUE OF X SOLVE THE SECOND EQ. FY=0 FOR Y

      IF(IT.EQ.1) THEN
        CALL SECAN1(X0,A,B,Y,REPS,AEPS,IER,FY)
      ELSE
        CALL BRENT1(A,B,Y,REPS,AEPS,IER,FY)
      ENDIF
      IF(IER.GT.0) STOP 11

!     STORE THE VALUE OF Y IN THE COMMON BLOCK FOR USE IN THE MAIN PROGRAM
      YY=Y
      FF=FY(Y)

!     THE FIRST EQUATION
      F=SIN(4/(X*X+Y*Y+0.1Q0))
      END

!     ------------------------------------------

      FUNCTION FY(Y)
      IMPLICIT REAL*16(A-H,O-Z)
      COMMON/FN/XX,YY,FF,IT

!     THE SECOND EQUATION
      FY=XX-TAN(Y)
      END
