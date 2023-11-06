!     PROGRAM TO INTEGRATE OSCILLATORY FUNCTIONS USING FILON'S FORMULA

      PROGRAM INTEG
!      IMPLICIT REAL*8(A-H,O,P,R-Z)
      IMPLICIT LOGICAL(Q)
      EXTERNAL FUN

!     EXERCISE 6.21 : INTEGRAL I1
 
51    FORMAT('   IER =',I4,5X,'A =',1PE14.6,5X,'B =',E14.6,5X,
     1       'RK =',E14.5/'  NO. OF FUNCTION EVALUATIONS =',I7,5X,
     2       'INTEGRAL =',E14.6/9X,'ESTIMATED ERROR =',E14.6,5X,
     3       'EXACT VALUE =',E14.6)

      REPS=1.E-6
      AEPS=1.E-8
      QSIN=.TRUE.
!      QSIN=.FALSE.
      B=4.*ACOS(0.D0)

100   PRINT *, 'TYPE  A=LOWER LIMIT,  B=UPPER LIMIT,'
      READ *,A,B
      PRINT *, 'RK=COEFFICIENT OF X IN OSCILLATORY FACTOR'
      PRINT *, '             (QUITS WHEN RK.EQ.0)'
      READ *,RK
      IF(RK.EQ.0) STOP
      CALL FILON(RI,A,B,RK,QSIN,REPS,AEPS,DIF,IER,N,FUN)
      REX=-B*RK/(RK*RK-1)
      WRITE(6,51) IER,A,B,RK,N,RI,DIF,REX
      GO TO 100
      END
 
!     ---------------------------------------------------
 
!	To calculate integrals with oscillatory integrand of the form
!	FUN(x)*SIN(RK*x)  or  FUN(x)*COS(RK*x)
!
!	RI : (output) Computed value of the integral
!	XL : (input) The lower limit
!	XU : (input) The upper limit
!	RK : (input) Coefficient of x in the oscillatory term in the integrand
!	QSIN : (input) Logical variable to specify the form of integrand
!		If QSIN=.TRUE. the oscillatory factor is SIN(RK*x)
!		If QSIN=.FALSE. the oscillatory factor is COS(RK*x)
!	REPS : (input) The required relative accuracy
!	AEPS : (input) The required absolute accuracy
!		The estimated error should be less than MAX(AEPS,REPS*ABS(RI))
!	DIF : (output) estimated (absolute) error achieved by the subroutine
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=30 implies specified accuracy was not achieved
!			DIF will contain the estimated accuracy
!	N : (output) Number of function evaluations used by subroutine
!	FUN : (input) Name of the function routine to calculate the integrand
!		excluding the oscillatory factor
!		FUNCTION FUN(X) must be supplied by the user.
!	
!	Required routines : FUN

      SUBROUTINE FILON(RI,XL,XU,RK,QSIN,REPS,AEPS,DIF,IER,N,FUN)
      IMPLICIT LOGICAL(Q)
!      IMPLICIT REAL*8(A-H,O,P,R-Z)
!	THC should be about (100*hcross)**(1/6)
!      PARAMETER(NMIN=5,NMAX=15,THC=0.005D0)
!	For REAL*4 version use this value
      PARAMETER(NMIN=5,NMAX=15,THC=0.15)

      IF(QSIN) THEN
!	For Sin(kx)
        FEND=FUN(XL)*COS(RK*XL)-FUN(XU)*COS(RK*XU)
        EVEN=0.5*(FUN(XL)*SIN(RK*XL)+FUN(XU)*SIN(RK*XU))
      ELSE
!	For Cos(kx)
        FEND=FUN(XU)*SIN(RK*XU)-FUN(XL)*SIN(RK*XL)
        EVEN=0.5*(FUN(XL)*COS(RK*XL)+FUN(XU)*COS(RK*XU))
      ENDIF

      ODD=0.0
      IER=0
      RI=0.0
      DIF=0.0
      N=2
      H=(XU-XL)
      IF(H.EQ.0.0) RETURN

      DO 3000 I=1,NMAX
        H=H/2.
        EVEN=EVEN+ODD
        ODD=0.0
        X1=XL+H
        H2=2.*H

!	Starting with 3 points subdivide the intervals into 2 until convergence
        DO 1000 J=1,N/2
          X=X1+H2*(J-1)
          IF(QSIN) THEN
            ODD=ODD+FUN(X)*SIN(RK*X)
          ELSE
            ODD=ODD+FUN(X)*COS(RK*X)
          ENDIF
1000    CONTINUE

        T=RK*H
        IF(ABS(T).GT.THC) THEN
!	Use normal functions
          ALPHA=(T*T+SIN(T)*(T*COS(T)-2.*SIN(T)))/T**3
          BETA=2.*(T+COS(T)*(T*COS(T)-2.*SIN(T)))/T**3
          GAMMA=4.*(SIN(T)-T*COS(T))/T**3
        ELSE
!	Use Taylor series expansion
          ALPHA=2.*T**3*(1.+T*T*(-1.+T*T/15.)/7.)/45.
          BETA=2.D0/3.+2.*T*T*(1.D0/15.+T*T*(-2.D0/105.+T*T/567.))
          GAMMA=4.D0/3.+T*T*(-2.D0/15.+T*T*(1.D0/210.-T*T/11340.))
        ENDIF
        R1=H*(ALPHA*FEND+BETA*EVEN+GAMMA*ODD)

        DIF=R1-RI
        RI=R1
        IF(I.LE.NMIN) GO TO 3000
        IF(ABS(DIF).LT.MAX(REPS*ABS(R1),AEPS)) RETURN
3000  N=N*2

      N=N/2
      IER=30
      END
 
!     ---------------------------------------------
 
      FUNCTION FUN(X)
!      IMPLICIT REAL*8(A-H,O,P,R-Z)

!     SPECIFY THE INTEGRAND (OMIT THE OSCILLATORY FACTOR COS(KX) OR SIN(KX))

      FUN=X*COS(X)
      END
