!     PROGRAM TO FIND COMPLEX ROOTS OF A NONLINEAR EQUATION
!     USING SECANT OR MULLER'S METHOD 

!	In some cases (e.g., with starting value (7,7)) SECANC
!	gives a spurious convergence to some arbitrary value.
!	This arises because the iteration suddenly tends to very
!	large value of X, where the function value is very large
!	and hence it comes back to the previous value of X and
!	converges. MULLER may not have this problem as it looks at
!	three previous values to calculate next iteration. This problem
!	can be avoided by giving more realistic limits on root. If
!	the subroutine is modified to keep within a rectangular region
!	around the expected root, then this problem will not arise if
!	the rectangle is chosen carefully.

      PROGRAM ZERO
      IMPLICIT COMPLEX*16(C)
      IMPLICIT REAL*8(A,B,D-H,O-Z)
      EXTERNAL CF
      DIMENSION CZERO(20)

!     EXAMPLE 7.6:  COMPLEX ZEROS OF Z+SIN(Z)

51    FORMAT('   IER =',I4,5X,' ZERO USING SECANC =',1P2D14.6)
52    FORMAT('   IER =',I4,5X,' ZERO USING MULLER =',1P2D14.6)
53    FORMAT('  STARTING VALUE =',1P2D14.6)

      REPS=1.D-7
      AEPS=1.D-7
      RMX=1.D2
      NZ=0


100   PRINT *,'TYPE CX0 = COMPLEX STARTING VALUE  (QUITS WHEN CX0=-100)'
      READ *,CX0
      IF(CX0.EQ.-100) STOP
      WRITE(6,53) CX0
      CALL SECANC(CX0,RMX,CX,REPS,AEPS,IER,CF)
      WRITE(6,51) IER,CX

      CX1=CX0*0.999D0
      CX2=CX0*1.001D0
      CALL MULLER(CX1,CX2,CX0,REPS,AEPS,IER,CF,NZ,CZERO,RMX)
      WRITE(6,52) IER,CX0
      GO TO 100
      END

!     -------------------------------------------------------

!	Complex zero of a given function using Muller iteration with deflation
!
!	CX1,CX2,CX3 : (input/output) Complex starting values for iteration
!		These will be updated during execution and CX3 should be
!		the best estimate for the zero.
!	REPS : (input) Required relative accuracy
!	AEPS : (input) Required absolute accuracy
!		The estimated error should be less than MAX(AEPS, REPS*ABS(CX3))
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=42 implies that roundoff errors appear to be dominating
!			and calculations are terminated.
!		IER=43 implies that iteration failed to converge to specified accuracy
!			but a weaker convergence criterion is satisfied
!		IER=404 implies 2 of the 3 starting values are identical and
!			no calculations are done
!		IER=431 implies that iteration goes outside the specified limits
!		IER=432 implies that iteration failed to converge to specified accuracy
!		IER=433 implies that denominator for calculating the iteration
!			vanishes and it is not possible to continue further
!	CF : (input) Name of the function routine to calculate the function value
!	NZ : (input) Number of zeros already known (for deflation)
!	CZERO : (input) Complex array of length NZ containing the known
!		zeros for deflation.
!	RMAX : (input) Maximum magnitude of zeros. Iteration will be terminated
!		when ABS(CX) > RMAX
!	
!	FUNCTION CF(CX) must be supplied by the user. Here CF and CX are
!		both complex variables.
!
!	Required routines : CF

      SUBROUTINE MULLER(CX1,CX2,CX3,REPS,AEPS,IER,CF,NZ,CZERO,RMAX)
      IMPLICIT COMPLEX*16(C)
      IMPLICIT REAL*8(A,B,D-H,O-Z)
      PARAMETER(REPS0=1.D-4,NIT=50)
      DIMENSION CZERO(NZ)

      IF(CX2.EQ.CX1.OR.CX2.EQ.CX3.OR.CX3.EQ.CX1) THEN
!	If two of the three starting values are equal then quit
        IER=404
        RETURN
      ENDIF

      CF1=CF(CX1)
      CF2=CF(CX2)
!	Perform deflation
      DO 2000 J=1,NZ
        CF1=CF1/(CX1-CZERO(J))
2000  CF2=CF2/(CX2-CZERO(J))

!	The divided difference f[CX1,CX2]
      CH1=(CF2-CF1)/(CX2-CX1)
      IER=0
      DX1=ABS(CX3-CX2)

!	Loop for Muller iteration
      DO 5000 I=1,NIT
        CF3=CF(CX3)
!	Perform deflation
        DO 3000 J=1,NZ
3000    CF3=CF3/(CX3-CZERO(J))

        IF(CX3.EQ.CX1.OR.CX3.EQ.CX2) RETURN
        IF(CF3.EQ.0.0) RETURN
        CH2=(CF3-CF2)/(CX3-CX2)
        C2=(CH2-CH1)/(CX3-CX1)
        CG=CH2+(CX3-CX2)*C2
        CI=SQRT(CG*CG-4.*CF3*C2)
        CD=CG+CI
        CD1=CG-CI
        IF(ABS(CD1).GT.ABS(CD)) CD=CD1

        IF(CD.EQ.0.0) THEN
!	If denominator is zero, then quit
          IER=433
          RETURN
        ENDIF

        CDX=-2.*CF3/CD
        CX1=CX2
        CX2=CX3
!	The new approximation to zero
        CX3=CX3+CDX
        CF1=CF2
        CF2=CF3
        CH1=CH2

        DX=ABS(CDX)
        IF(I.GT.2.AND.DX.LT.MAX(REPS*ABS(CX3),AEPS)) RETURN
        IF(I.GT.5.AND.DX.LT.MAX(REPS0*ABS(CX3),AEPS)
     1      .AND.DX.GT.DX1) THEN
!	Roundoff errors appear to be dominating, hence quit
          IER=42
          RETURN
        ENDIF

        DX1=DX
        IF(ABS(CX3).GT.RMAX) THEN
!	Iteration goes outside the specified range
          IER=431
          RETURN
        ENDIF
5000  CONTINUE

!	Iteration fails to converge.
      IER=43
      IF(DX.LT.MAX(REPS0*ABS(CX3),AEPS)) RETURN
      IER=432
      END

!     -------------------------------------------------------------

!	Complex zero of a given function using secant iteration 
!
!	X0 : (input) Initial guess for the zero
!	R : (input) limiting magnitude for the zero
!	X : (output) Computed value of the zero
!	REPS : (input) Required relative accuracy
!	AEPS : (input) Required absolute accuracy
!		The estimated error should be less than MAX(AEPS, REPS*ABS(X))
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=40 implies that function value is equal at two points
!			and it is not possible to continue the iteration
!		IER=402 implies ABS(X0)>R in which case no calculations are done
!		IER=422 implies that iteration goes outside the specified limits
!		IER=423 implies that iteration failed to converge to specified accuracy
!	FUN : (input) Name of the function routine to calculate the function
!		FUNCTION FUN(X) must be supplied by the user.
!
!	Required routines : FUN

      SUBROUTINE SECANC(X0,R,X,REPS,AEPS,IER,FUN)
      IMPLICIT REAL*8(A-C,O-T)
      IMPLICIT COMPLEX*16(D-H,U-Z)
      PARAMETER(NIS=75)

      IER=0
      IF(ABS(X0).GT.R) THEN
!	If X0 is outside the specified interval (XL,XU) then quit
        IER=402
        RETURN
      ENDIF

      X=X0
!	Select the increment for the next point X+DX
      DX=X*0.001D0
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
        IF(ABS(X).GT.R) THEN
!	If iteration goes outside the specified limits (XL,XU), then quit
          IER=422
          RETURN
        ENDIF
1000  CONTINUE

!	The iteration fails to converge
      IER=423
      END

!     ----------------------------------------------

      FUNCTION CF(Z)
      IMPLICIT COMPLEX*16(C,Z)

!	The required function

      CF=Z+SIN(Z)
      END
