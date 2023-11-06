!	To calculate the spherical Bessel function of integral order
!	(j_n(x)=Sqrt(PI/(2x))*J_{n+1/2}(x)) for a real argument
!	Because of low exponent range this routine may give overflow in
!	many cases, even when the actual value is within acceptable limits
!
!	N : (input) Order of Bessel function required, N may be negative
!		or positive. 
!	XB : (input) Argument at which the value is required
!	BJ : (output) Real array of length at least 
!		ABS(N)+16+MAX(25,5*SQRT(N))
!		which will contain the value of Bessel function of order
!		0,1,...,ABS(N). BJ(I+1) will contain Bessel function of order I
!		or -I (if N<0)
!		Remaining elements of array are used as scratch space
!
!	Required routines : None
 
      SUBROUTINE SPHBJN(N,XB,BJ)
!      IMPLICIT REAL*8(A-H,O-Z)
!	REPS should be less than the machine accuracy
      PARAMETER(REPS=1.E-8)
      DIMENSION BJ(*)
 
51    FORMAT(' SPHBJN FAILED AT N =',I6,'  X =',1PE13.5,'  S =',E13.5)
 
      X=ABS(XB)
      NA=ABS(N)
      IF(XB.EQ.0) THEN
        BJ(1)=1.0
        DO 1800 I=2,NA+1
          BJ(I)=0.0
1800    CONTINUE
        RETURN
      ENDIF
      IF(N.GT.0) THEN
        IF(N.LT.X.OR.N.LT.2) THEN
!	Use the recurrence relation in the forward direction 
 
          BJ(1)=SIN(X)/X
          BJ(2)=SIN(X)/X**2-COS(X)/X
          DO 2000 I=2,N
            BJ(I+1)=(2*I-1)*BJ(I)/X-BJ(I-1)
2000      CONTINUE

        ELSE IF(X.LE.4.0) THEN

!	Use series expansion to calculate  BJ(NA), BJ(NA-1)
          XA=X*X/4.0
          T0=1.0
          DO 2400 I=1,NA
            T0=T0*X/(2.0*I+1)
            IF(I.GE.NA-1) THEN
              T=T0
              S=T0
              DO 2200 J=1,50
                T=-T*XA/(J*(J+I+0.5D0))
                S=S+T
                IF(ABS(T).LT.ABS(S*REPS)) GO TO 2300
2200          CONTINUE
2300          BJ(I+1)=S
            ENDIF
2400      CONTINUE
 
          DO 2600 I=NA-1,1,-1
            BJ(I)=(2.*I+1)*BJ(I+1)/X-BJ(I+2)
2600      CONTINUE
        ELSE
 
!	Use the recurrence relation in the backward direction 
          N1=NA+MAX(15.0,3.*SQRT(1.0*NA))
          IF(X.LT.NA/3.0) N1=NA+MAX(15.0,3.*SQRT(1.0*NA))/LOG(NA/X)
          IF(NA.LT.X+15) N1=N1+10
          IF(N1-NA.LT.12) N1=12+NA
          BJ(N1+1)=0.0
          BJ(N1)=1.0
          DO 3200 I=N1-1,1,-1
            BJ(I)=(2.*I+1)*BJ(I+1)/X-BJ(I+2)
3200      CONTINUE

          S=BJ(1)*X/SIN(X)
          DO 3400 I=1,NA+1
            BJ(I)=BJ(I)/S
3400      CONTINUE

!	If ABS(BJ(NA+2))<1/REPS, then the required accuracy may not be achieved
          IF(ABS(BJ(NA+2))*REPS.LT.1.0) PRINT 51,N,X,BJ(NA+2)
        ENDIF
 
      ELSE
 
!	For negative N use the recurrence relation in the forward direction 
        BJ(1)=SIN(X)/X
        BJ(2)=COS(X)/X
        DO 3600 I=2,NA
          BJ(I+1)=(-2*I+3)*BJ(I)/X-BJ(I-1)
3600    CONTINUE
      ENDIF
 
      IF(XB.LT.0.0) THEN
        DO 4000 I=2,ABS(N)+1,2
          BJ(I)=-BJ(I)
4000    CONTINUE
      ENDIF
      END
 
