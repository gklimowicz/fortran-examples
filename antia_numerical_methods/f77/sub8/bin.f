!	To calculate the modified Bessel function of first kind of
!		integral order for real argument
!
!	N : (input) Order of Bessel function required, N must be positive
!	XB : (input) Argument at which the value is required
!	BI : (output) Real array of length at least N+16+MAX(25, 5*SQRT(N))
!		which will contain the value of Bessel function of order
!		0,1,...,N. BI(I+1) will contain Bessel function of order I
!		Remaining elements of array are used as scratch space
!
!	Required routines : BI0, BI1
!
 
      SUBROUTINE BIN(N,XB,BI)
      IMPLICIT REAL*8(A-H,O-Z)
!	REPS should be less than the machine accuracy
      PARAMETER(REPS=1.D-17)
      DIMENSION BI(*)
 
51    FORMAT('  BIN FAILED AT N =',I6,'  X =',1PD13.5,'  S =',D13.5)
 
      X=ABS(XB)
      NA=ABS(N)
      IF(XB.EQ.0.0) THEN
        BI(1)=1.0
        DO 1500 I=2,NA+1
          BI(I)=0.0
1500    CONTINUE
        RETURN
      ENDIF

      IF(NA.LT.X-10.OR.NA.LT.2) THEN
!	Use the recurrence relation in forward direction
 
        BI(1)=BI0(X)
        BI(2)=BI1(X)
        DO 2000 I=2,NA
          BI(I+1)=-2*(I-1)*BI(I)/X+BI(I-1)
2000    CONTINUE

      ELSE IF(X.LE.4.0) THEN

!	Use series expansion to calculate BI(NA), BI(NA-1)
        XA=X*X/4.0
        T0=X/2.0
        DO 2400 I=2,NA
          T0=T0*X/(2.0*I)
          IF(I.GE.NA-1) THEN
            T=T0
            S=T0
            DO 2200 J=1,50
              T=T*XA/(J*(J+I))
              S=S+T
              IF(ABS(T).LT.ABS(S*REPS)) GO TO 2300
2200        CONTINUE
2300        BI(I+1)=S
          ENDIF
2400    CONTINUE

        DO 2600 I=NA-1,1,-1
          BI(I)=2.*I*BI(I+1)/X+BI(I+2)
2600    CONTINUE
      ELSE
 
!	Use the recurrence relation in backward direction
        N1=NA+MAX(25.D0,5.*SQRT(1.D0*NA))
        IF(X.LT.NA/6.0) N1=NA+MAX(25.D0,5.*SQRT(1.D0*NA))/LOG(NA*0.5/X)
        IF(NA.LT.X+15) N1=N1+15
        IF(N1-NA.LT.15) N1=15+NA
        BI(N1+1)=0.0
        BI(N1)=1.0
        DO 3200 I=N1-1,1,-1
          BI(I)=2.*I*BI(I+1)/X+BI(I+2)
3200    CONTINUE

        S=BI(1)/BI0(X)
        DO 3600 I=1,NA+1
          BI(I)=BI(I)/S
3600    CONTINUE
!	If ABS(BI(NA+2))<1/REPS, then the required accuracy may not be achieved
!	hence printout an error message
        IF(ABS(BI(NA+2))*REPS.LT.1.0) PRINT 51,N,X,BI(NA+2)
      ENDIF
 
      IF(XB.LT.0.0) THEN
        DO 3800 I=2,NA+1,2
          BI(I)=-BI(I)
3800    CONTINUE
      ENDIF
      END
