!     PROGRAM FOR FINDING REAL ROOTS OF A NONLINEAR EQUATION IN A
!     SPECIFIED INTERVAL.  THE ROOTS ARE LOCATED BY LOOKING FOR SIGN
!     CHANGES IN THE FUNCTION VALUE AT THE SPECIFIED SEARCH STEP. 
!     AFTER THAT SECANT ITERATION OR BRENT'S METHOD IS USED TO CALCULATE
!     THE ROOT.

      PROGRAM ROOT
!      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL F

!     EXAMPLE 7.2

51    FORMAT('   SEARCH STEP =',1PE12.4,3X,'FOR ROOTS IN (',E12.4,
     1       ',',E12.4,')    IT=',I2)
52    FORMAT('  IER =',I4,3X,'ROOT IN (',1PE12.4,',',E12.4,') ='
     1       ,E14.6,3X,'F =',E10.2)
53    FORMAT('   IER =',I4,5X,'ROOT =',1PE14.6,5X,'F =',E14.6)

      REPS=1.E-6
      AEPS=1.E-8

!     THE SEARCH STEP SHOULD BE CHOSEN CAREFULLY, SINCE A PAIR OF
!     CLOSE ROOTS MAY BE MISSED IF THE STEP IS LARGE. SINGULARITIES
!     WHERE THE FUNCTION CHANGES SIGN ARE ALSO LOCATED AS ROOTS.

100   PRINT *,'TYPE DX=SEARCH STEP, XL=LOWER LIMIT, XU=UPPER LIMIT'
      PRINT *,'           (QUITS WHEN DX(XL-XU).GT.0)'
      READ *,DX,XL,XU
      IF(DX.GT.0.0.EQV.XL-XU.GT.0.0) STOP
      PRINT *,'TYPE IT=1/2    FOR   SECANT/BRENT'
      READ *,IT
      WRITE(6,51) DX,XL,XU,IT

!     LOCATE THE ROOTS BY LOOKING FOR SIGN CHANGES

      F1=F(XL)
      DO 1000 X1=XL+DX,XU,DX
        FX=F(X1)
!        IF(.NOT.(FX.GT.0.0.EQV.F1.GT.0.0)) THEN
!       To eliminate the singularity use the following condition
        IF(.NOT.(FX.GT.0.0.EQV.F1.GT.0.0).AND.
     1    (TAN(X1)*TAN(X1-DX).GT.0.0)) THEN

!     FIND ACCURATE VALUE OF THE ROOT

          IF(DX.GT.0.0) THEN
            A=X1-DX
            B=X1
          ELSE
            B=X1-DX
            A=X1
          ENDIF
          X0=(A+B)*0.5
          IF(IT.EQ.1) THEN
            CALL SECANT(X0,A,B,X,REPS,AEPS,IER,F)
            WRITE(6,52) IER,A,B,X,F(X)
          ELSE
            CALL BRENT(A,B,X,REPS,AEPS,IER,F)
            WRITE(6,53) IER,X,F(X)
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
!      IMPLICIT REAL*8(A-H,O-Z)
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

!     ----------------------------------------------------

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
!      IMPLICIT REAL*8(A-H,O-Z)
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
      IF(ABS(DX).GT.0.1D0*MAX(ABS(X),100.*AEPS))
     1    DX=SIGN(0.10*MAX(ABS(X),100.*AEPS),DX)
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
!      IMPLICIT REAL*8(A-H,O-Z)

!	THE REQUIRED FUNCTION WHOSE ZERO IS REQUIRED
      F=X-TAN(X)
      END
