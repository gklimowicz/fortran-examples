!	To calculate the Bessel function of integral order for real argument
!
!	N : (input) Order of Bessel function required, N may be negative
!		or positive. For N=0,1 use BJ0 and BJ1 respectively.
!	XB : (input) Argument at which the value is required
!	BJ : (output) Real array of length at least 
!		ABS(N)+16+MAX(25,5*SQRT(N))
!		which will contain the value of Bessel function of order
!		0,1,...,ABS(N). BJ(I+1) will contain Bessel function of order I
!		or -I (if N<0)
!		Remaining elements of array are used as scratch space
!
!	Required routines : BJ0, BJ1
!
      SUBROUTINE BJN(N,XB,BJ)
      IMPLICIT REAL*16(A-H,O-Z)
!	REPS should be normally less than machine accuracy, but
!	since BJ0 and BJ1 are only calculated to accuracy of 1.D-15
!	it is not set to machine precision.
      PARAMETER(REPS=1.Q-19)
      DIMENSION BJ(*)
 
51    FORMAT('  BJN FAILED AT N =',I6,'  X =',1PD13.5,'  S =',D13.5)
 
      X=ABS(XB)
      NA=ABS(N)
      IF(XB.EQ.0) THEN
        BJ(1)=1.0
        DO 1800 I=2,NA+1
          BJ(I)=0.0
1800    CONTINUE
        RETURN

      ELSE IF(NA.LT.X.OR.NA.LT.2) THEN
!	Use the recurrence relation in the forward direction 
        BJ(1)=BJ0(X)
        BJ(2)=BJ1(X)
        DO 2000 I=2,NA
          BJ(I+1)=2*(I-1)*BJ(I)/X-BJ(I-1)
2000    CONTINUE
      ELSE IF(X.LE.4.0) THEN

!	Use series expansion to calculate  BJ(NA), BJ(NA-1)
        XA=X*X/4.0
        T0=X/2.0
        DO 2400 I=2,NA
          T0=T0*X/(2.0*I)
          IF(I.GE.NA-1) THEN
            T=T0
            S=T0
            DO 2200 J=1,50
              T=-T*XA/(J*(J+I))
              S=S+T
              IF(ABS(T).LT.ABS(S*REPS)) GO TO 2300
2200        CONTINUE
2300        BJ(I+1)=S
          ENDIF
2400    CONTINUE
 
        DO 2600 I=NA-1,1,-1
          BJ(I)=2.*I*BJ(I+1)/X-BJ(I+2)
2600    CONTINUE

      ELSE
!	Use the recurrence relation in the backward direction 
        N1=NA+MAX(25.Q0,5.*SQRT(1.Q0*NA))
        IF(X.LT.NA/6.0) N1=NA+MAX(25.Q0,5.*SQRT(1.Q0*NA))/LOG(NA*0.5/X)
        IF(NA.LT.X+15) N1=N1+15
        IF(N1-NA.LT.15) N1=15+NA
        BJ(N1+1)=0.0
        BJ(N1)=1.0
        DO 3200 I=N1-1,1,-1
          BJ(I)=2.*I*BJ(I+1)/X-BJ(I+2)
3200    CONTINUE

        S=BJ(1)
        DO 3400 I=3,N1,2
          S=S+2*BJ(I)
3400    CONTINUE

        DO 3600 I=1,NA+1
          BJ(I)=BJ(I)/S
3600    CONTINUE
!	If ABS(BJ(NA+2))<1/REPS, then the required accuracy may not be achieved
!	hence printout an error message
        IF(ABS(BJ(NA+2))*REPS.LT.1.0) PRINT 51,N,X,BJ(NA+2)
      ENDIF
 
!       IF(XB*N.LT.0.0) THEN    ! use this if the following doesn't work
      IF(N.LT.0.XOR.XB.LT.0.0) THEN
        DO 3800 I=2,NA+1,2
          BJ(I)=-BJ(I)
3800    CONTINUE
      ENDIF
      END
