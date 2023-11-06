!	To minimise a function in one dimension using Hermite cubic interpolation
!
!	X1,X2 : (input/output) Two starting values for iteration
!		These will be updated during execution
!		with X2 giving the best estimate for minimiser.
!	F2 : (output) The function value at X2, which should be the minimum
!	D2F : (output) Estimated value of second derivative at X2, to detect
!		if the extremum is a minimum or maximum. If D2F>0 then
!		X2 should be minimum, otherwise it will be maximum
!	REPS : (input) Required relative accuracy
!	AEPS : (input) Required absolute accuracy
!		The Minimiser will be calculated with accuracy MAX(AEPS, REPS*ABS(X2))
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=52 implies that iteration has converged to a maximum
!		IER=502 implies that initial values X1=X2 and no calculations are done
!		IER=524 implies that iteration cannot be continued 
!		IER=525 implies that iteration failed to converge to specified accuracy
!	F : (input) Name of the function routine to calculate the function
!		which is to be minimised
!
!	FUNCTION F(X,DX) to calculate the required function, must be supplied
!		by the user. Here DX is the first derivative of F(X) at X
!
!	Required routines : F

      SUBROUTINE DAVIDM(X1,X2,F2,D2F,REPS,AEPS,IER,F)
      IMPLICIT REAL*16(A-H,O-Z)
      PARAMETER(NIT=75)

      IER=0
      IF(X1.EQ.X2) THEN
        IER=502
        RETURN
      ENDIF

      F1=F(X1,DF1)
      F2=F(X2,DF2)
      DX1=X2-X1

!	Loop for iteration
      DO 2000 I=1,NIT
        DF12=(F2-F1)/DX1
        R=2.*DF2+DF1-3.*DF12
        R1=(3.*DF12-DF2-DF1)**2-DF1*DF2
        IF(R1.GT.0.0) THEN
!	Perform Hermite cubic interpolation
          R1=SIGN(SQRT(R1),R)
          R=R+R1
          DX=-DF2*DX1/R
        ELSE
!	Perform parabolic interpolation
          R2=2.*(DF12-DF2)
          IF(R2.NE.0.0) THEN
            DX=DX1*DF2/R2
          ELSE
!	Both cubic and parabolic interpolation fail, hence quit
            IER=524
            RETURN
          ENDIF
        ENDIF

        X1=X2
        X2=X2+DX
        D2F=(DF2-DF1)/DX1
        DX1=DX
        F1=F2
        DF1=DF2
        F2=F(X2,DF2)

!	Convergence check
        IF(X1.EQ.X2.OR.(ABS(DX).LT.MAX(REPS*ABS(X2),AEPS).AND.I.GT.2))
     1     THEN
!	If it is maximum set error flag
          IF(D2F.LE.0) IER=52
          RETURN
        ENDIF
2000  CONTINUE

!	Iteration fails to converge
      IER=525
      END
