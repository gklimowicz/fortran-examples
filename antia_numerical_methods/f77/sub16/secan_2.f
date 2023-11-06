!	Real zero of a given function using secant iteration
!	Function is calculated as FX*2**JF
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
!		FUNCTION FUN(X,JF) must be supplied by the user.
!		The function value should be FUN*2**JF
!
!	Required routines : FUN

      SUBROUTINE SECAN_2(X0,XL,XU,X,REPS,AEPS,IER,FUN)
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
      JF1=0

      DO 1000 L=1,NIS
        F=FUN(X,JF)
        DX1=DX
        F1=F1*2.Q0**(JF1-JF)

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
        JF1=JF

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
