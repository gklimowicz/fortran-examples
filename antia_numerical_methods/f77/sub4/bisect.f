!	Real zero of a continuous function using bisection
!
!	X : (output) Computed value of the zero using interpolation
!	XL : (input/output) Lower limit of interval containing the zero
!	XU : (input/output) Upper limit of interval containing the zero
!		These limits will be replaced by refined limits after execution
!	NB : (input) Number of bisections to be performed
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=-1 implies that function value is zero at some trial point
!		IER=401 implies NB.LE.0 in which case no calculations are done
!		IER=421 implies that function has same sign at both limits
!		and it is not possible to apply bisection
!	F : (input) Name of the function routine to calculate the function
!		FUNCTION F(X) must be supplied by the user.
!
!	Required routines : F

      SUBROUTINE BISECT(X,XL,XU,NB,IER,F)
!      IMPLICIT REAL*8(A-H,O-Z)

      IF(NB.LE.0) THEN
        IER=401
        RETURN
      ENDIF

      FL=F(XL)
      FU=F(XU)
      IER=-1
!	If the function value is zero at either point then quit
      IF(FL.EQ.0.0) THEN
        X=XL
        RETURN
      ELSE IF(FU.EQ.0.0) THEN
        X=XU
        RETURN
      ENDIF

      IF(FL.GT.0.0.EQV.FU.GT.0.0) THEN
!	If the function has the same sign at both end points then quit
        IER=421
        RETURN
      ENDIF

!	Loop for bisection
      DO 1000 K=1,NB
        X=(XL+XU)/2.
        FX=F(X)
!	If function is zero then quit
        IF(FX.EQ.0.0) RETURN
        IF(FX.GT.0.0.EQV.FU.GT.0.0) THEN
          XU=X
          FU=FX
        ELSE
          XL=X
          FL=FX
        ENDIF
1000  CONTINUE

!	linear interpolation between XL and XU
      X=(XL*FU-XU*FL)/(FU-FL)
      IER=0
      END
