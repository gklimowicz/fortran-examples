!	Real zero of a given function using Newton-Raphson iteration
!
!	X0 : (input) Initial guess for the zero
!	XL : (input) Lower limit of interval where zero is expected
!	XU : (input) Upper limit of interval where zero is expected
!	X : (output) Computed value of the zero
!	REPS : (input) Required relative accuracy
!	AEPS : (input) Required absolute accuracy
!		The estimated error should be less than MAX(AEPS, REPS*ABS(X))
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=-k implies that the zero is multiple and iteration
!			is modified to account for multiplicity k.
!		IER=126-100*k implies that the derivative is zero and the zero is
!			detected to be multiple. Iteration is terminated.
!		IER=403 implies XL>X0 or XU<X0, in which case no calculations are done
!		IER=424 implies that iteration goes outside the specified limits
!		IER=425 implies that iteration failed to converge to specified accuracy
!		IER=426 implies that the derivative is zero and iteration
!			is terminated
!	FUN : (input) Name of the function routine to calculate the function
!		FUNCTION FUN(X,DX) must be supplied by the user.
!		Here DX is the first derivative of FUN at X.
!
!	Required routines : FUN

      SUBROUTINE NEWRAP(X0,XL,XU,X,REPS,AEPS,IER,FUN)
!      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NIT=75)

      IER=0
      IF(XL.GT.X0.OR.XU.LT.X0) THEN
!	If X0 is outside the specified limits then quit
        IER=403
        RETURN
      ENDIF

      MR=1
      X=X0
      DX=1
      DXR=1
      L0=1

      DO 1000 L=1,NIT
        F=FUN(X,DF)
        DX1=DX
        DXR1=DXR

        IF(DF.EQ.0.0) THEN
!	If the derivative is zero, then quit
          IF(F.EQ.0.0) RETURN
!	If the function is nonzero, then set the error flag
          IER=426
          IF(MR.GT.1) IER=126-100*IER
          RETURN
        ENDIF

!	Newton-Raphson iteration for root with multiplicity MR
        DX=-MR*F/DF
        DXR=DX/DX1
        X=X+DX

        IF(ABS(DX).LT.MAX(REPS*ABS(X),AEPS).AND.L.GT.2) RETURN
        IF(X.LT.XL.OR.X.GT.XU) THEN
!	If the iteration goes outside the specified limits, then quit
          IER=424
          RETURN
        ENDIF

        IF(L-L0.GT.3) THEN
!	If 3 iterations have been done with same multiplicity, then
!	get a new estimate for multiplicity MR
          MR1=MR
          MR2=MR
          IF(DXR.LT.0.99D0) MR1=MR/(1.-DXR)+0.5
          IF(DXR1.LT.0.99D0) MR2=MR/(1.-DXR1)+0.5
!	Accept the new value of MR only if both estimates match
          IF(MR1.EQ.MR2.AND.MR1.NE.MR) THEN
            IER=-MR1
            L0=L
            MR=MR1
            IF(MR.EQ.1) IER=0
          ENDIF
        ENDIF
1000  CONTINUE

!	Iteration fails to converge
      IER=425
      END
