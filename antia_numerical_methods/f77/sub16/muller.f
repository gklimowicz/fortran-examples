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
