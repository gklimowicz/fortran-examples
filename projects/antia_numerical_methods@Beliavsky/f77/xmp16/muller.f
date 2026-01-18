!     PROGRAM TO FIND COMPLEX ROOTS OF A NONLINEAR EQUATION
!     USING MULLER'S METHOD WITH DEFLATION TO REMOVE KNOWN ZEROS

      PROGRAM ZERO
      IMPLICIT COMPLEX*32(C)
      IMPLICIT REAL*16(A,B,D-H,O-Z)
      EXTERNAL CF
      DIMENSION CZERO(20),CX(20)

!     EXAMPLE 7.6:  COMPLEX ZEROS OF Z+SIN(Z)

51    FORMAT('   IER =',I4,5X,I3,' ZEROS CALCULATED')
52    FORMAT(2X,1P2D14.6,2X,2D14.6)
53    FORMAT('  STARTING VALUES'/(2X,1P2D14.6,2X,2D14.6))

      REPS=1.Q-30
      AEPS=1.Q-33
      RMX=1.Q2
      NZ=0

!     THE ZEROS WILL KEEP ACCUMULATING IN THE ARRAY CZERO

100   PRINT *,'TYPE NUM=NO. OF ZEROS TO BE TRIED  (QUITS WHEN NUM.EQ.0)' 
      READ *,NUM
      IF(NUM.EQ.0) STOP
      PRINT *,'TYPE NUM COMPLEX STARTING VALUES'
      READ *,(CX(I),I=1,NUM)
      WRITE(6,53) (CX(I),I=1,NUM)
      CALL ZROOT(NUM,CX,CZERO,NZ,REPS,AEPS,IER,RMX,CF)
      WRITE(6,51) IER,NZ
      WRITE(6,52) (CZERO(I),I=1,NZ)
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
      IMPLICIT COMPLEX*32(C)
      IMPLICIT REAL*16(A,B,D-H,O-Z)
      PARAMETER(REPS0=1.Q-4,NIT=50)
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

!     ---------------------------------------------

!	Complex zeros of a given function using Muller iteration with deflation
!
!	N : (input) Number of zeros to be determined
!	CX : (input) Complex array of length N containing the initial guesses
!		for the zeros
!	CZERO : (input/output) Complex array of length N+NZ containing the
!		computed values of the zeros. The first NZ zeros which are
!		already known must be supplied as input while other zeros will be added
!	NZ : (input/output) Number of known zeros. At input it should
!		contain the number already known which are included in
!		array CZERO. This number will be incremented as more zeros
!		are determined successfully.
!	REPS : (input) Required relative accuracy
!	AEPS : (input) Required absolute accuracy
!		The estimated error should be less than MAX(AEPS, REPS*ABS(X))
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=41 implies that iteration did not converge to satisfactory accuracy for
!			at least one zero, but value is acceptable at lower accuracy.
!		IER=429 implies that iteration failed to find at least one zero.
!			The number of zeros found can be checked from the value of NZ.
!	RMAX : (input) Maximum magnitude of zeros. Iteration will be terminated
!		when ABS(CZ) > RMAX
!	CF : (input) Name of the function routine to calculate the function
!		FUNCTION CF(CZ) must be supplied by the user, where CF and CZ
!		are complex. For use with MULER2 the function is calculated
!		in the form CF*2**JF and FUNCTION CF(CZ,JF) should be provided.
!
!	Required routines : MULLER (or MULER2), CF

      SUBROUTINE ZROOT(N,CX,CZERO,NZ,REPS,AEPS,IER,RMAX,CF)
      IMPLICIT COMPLEX*32(C)
      IMPLICIT REAL*16(A,B,D-H,O-Z)
      EXTERNAL CF
      DIMENSION CZERO(NZ+N),CX(N)

      IER=0
      IF(N.LE.0) RETURN
      IF(NZ.LT.0) NZ=0

      DO 1000 I=1,N
!	The starting values for Muller's iteration
        CX3=CX(I)
        CDX=0.01Q0*CX3
        IF(I.LT.N.AND.CDX.EQ.0.0) CDX=0.1Q0*(CX(I+1)-CX(I))
        IF(I.GT.1.AND.CDX.EQ.0.0) CDX=0.1Q0*(CX(I-1)-CX(I))
        IF(CDX.EQ.0.0) CDX=1.Q-3
!	These values may need to be changed in some cases
        CX2=CX3+CDX
        CX1=CX3-CDX

!	Find the next zero
        CALL MULLER(CX1,CX2,CX3,REPS,AEPS,IER1,CF,NZ,CZERO,RMAX)
!	For MULER2 use the following statements instead
!         IER1=0
! 500     CALL MULER2(CX1,CX2,CX3,REPS,AEPS,IER1,CF0,CX,IX,NZ,CZERO,RMAX)
!         IF(IER1.LT.0) THEN
!           CF0=CF(CX,IX)
!           GO TO 500
!         ENDIF

        IF(IER1.LT.100) THEN
!	The zero is accepted
          NZ=NZ+1
          CZERO(NZ)=CX3
          IF(IER.EQ.0.AND.IER1.NE.0) IER=41
        ELSE
          IER=429
        ENDIF
1000  CONTINUE
      END

!     ----------------------------------------------

      FUNCTION CF(Z)
      IMPLICIT COMPLEX*32(C,Z)

!	The required function for finding zeros
      CF=Z+SIN(Z)
      END
