!     Real zero of a given function using secant iteration
!     Function is calculated as FX*2**JF
!     This subroutine uses reverse communication to calculate function
!     values. If IER<0 the function should be evaluated and SECANI should
!     be called back with new function value. Calculation of function
!     value should not change any other variables in the call statement.
!
!     X0 : (input) Initial guess for the zero
!     XL : (input) Lower limit of interval where zero is expected
!     XU : (input) Upper limit of interval where zero is expected
!     X : (output) Value of x at which the function evaluation is required.
!		If IER=0 then it will contain the final value of zero computed
!		by the routine.
!     F : (input) Calculated value of the function at X.
!		If subroutine exits with IER<0, then the calling routine should
!		calculate the function value at X and call SECANI with this value
!		stored in F and JF. Other variables should not be changed.
!     JF : (input) The exponent of function value, the function value
!		should be F*2**JF
!     REPS : (input) Required relative accuracy
!     AEPS : (input) Required absolute accuracy
!     		The estimated error should be less than MAX(AEPS, REPS*ABS(X))
!     IER : (input/output) Error parameter, IER=0 implies successful execution
!		Before the first call IER should be set to zero
!		IER<0 implies that execution is not over and the subroutine needs
!			a new function evaluation at X. After calculating the
!			function value SECANI should be called back.
!     		IER=40 implies that function value is equal at two points
!     			and it is not possible to continue the iteration
!     		IER=402 implies XL>X0 or XU<X0, in which case no calculations are done
!     		IER=422 implies that iteration goes outside the specified limits
!     		IER=423 implies that iteration failed to converge to specified accuracy
!
!	Required routines : None (Function is calculated by the calling program)

      SUBROUTINE SECANI(X0,XL,XU,X,F,JF,REPS,AEPS,IER)
      IMPLICIT REAL*16(A-H,O-Z)
      PARAMETER(NIS=75)
      SAVE
 
      IFL=-IER
!     Jump to proper location and continue execution
      IF(IFL.EQ.1) GO TO 1500

!     Initial call to subroutine, start from beginning
      IF(XL.GT.X0.OR.XU.LT.X0) THEN
!     If X0 is outside the specified interval (XL,XU) then quit
        IER=402
        RETURN
      ENDIF
 
      X=X0
!     Select the increment for the next point X+DX
      DX=(XU-X0)/100.
      IF(ABS(DX).LT.5.*MAX(REPS*ABS(X),AEPS)) DX=(XL-X0)/5.
      IF(ABS(DX).GT.0.1Q0*MAX(ABS(X),100.*AEPS))
     1       DX=SIGN(0.1Q0*MAX(ABS(X),100.*AEPS),DX)
      F1=0.0
      JF1=0
      L=0
 
1000  L=L+1
      IER=-1
!     To evaluate the function at X
      RETURN

1500  DX1=DX
      F1=F1*2.Q0**(JF1-JF)
      IER=0
 
      IF(F1-F.EQ.0.0) THEN
        IF(F.EQ.0.0) RETURN
!     If F1=F and F.NE.0, then quit
        IER=40
        RETURN
      ENDIF
 
!     The secant iteration
      IF(L.GT.1) DX=DX1*F/(F1-F)
      X=X+DX
      F1=F
      JF1=JF
 
      IF(ABS(DX).LT.MAX(REPS*ABS(X),AEPS).AND.L.GT.2) RETURN
      IF(X.LT.XL.OR.X.GT.XU) THEN
!     If iteration goes outside the specified limits (XL,XU), then quit
        IER=422
        RETURN
      ENDIF
      IF(L.LT.NIS) GO TO 1000
 
!     The iteration fails to converge
      IER=423
      END
