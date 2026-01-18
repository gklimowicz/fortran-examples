!	Complex zero of a given function using Muller iteration with deflation
!	Function is assumed to be calculated as CF*2**IX
!	This subroutine uses reverse communication to calculate function
!	values. If IER<0 the function should be evaluated and MULER2
!	should be called back with new function value. Calculation of
!	function value should not change any other variables in the
!	call statement. First call should be with IER=0
!
!	CX1,CX2,CX3 : (input/output) Complex starting values for iteration
!		These will be updated during execution and CX3 should be
!		the best estimate for the zero.
!	REPS : (input) Required relative accuracy
!	AEPS : (input) Required absolute accuracy
!		The estimated error should be less than MAX(AEPS, REPS*ABS(CX3))
!	IER : (input/output) Error parameter, IER=0 implies successful execution
!		For the first call IER should be set to 0.
!		IER<0 implies that execution is not over and the subroutine
!			needs a new function evaluation at z=CX. After calculating
!			the function value MULER2 should be called back.
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
!	CF : (input) Calculated value of the function at z=CX
!		If subroutine exits with IER<0, then the calling routine
!		should calculate the function value at CX and call MULER2
!		with this value stored in CF and IX. Other variables should
!		not be changed.
!	CX : (output) Value of z at which the function evaluation is needed
!	IX : (input) The exponent of function value, The function value 
!		should be CF*2**IX
!	NZ : (input) Number of zeros already known (for deflation)
!	CZERO : (input) Complex array of length NZ containing the known
!		zeros for deflation.
!	RMAX : (input) Maximum magnitude of zeros. Iteration will be terminated
!		when ABS(CX) > RMAX
!	
!	Required routines : None (Function has to be calculated by calling program)
!
      SUBROUTINE MULER2(CX1,CX2,CX3,REPS,AEPS,IER,CF,CX,IX,NZ,
     1                  CZERO,RMAX)
      IMPLICIT COMPLEX*8(C)
      IMPLICIT REAL*4(A,B,D-H,O-Z)
      PARAMETER(REPS0=1.D-3,NIT=50)
      DIMENSION CZERO(NZ)
      SAVE

      IFL=-IER
      IF(IFL.GT.0.AND.IFL.LT.4) THEN
!	Jump to proper location and continue execution
        GO TO (1000,2000,3000), IFL
      ENDIF
!	Initial call to subroutine, start from beginning
      IF(CX2.EQ.CX1.OR.CX2.EQ.CX3.OR.CX3.EQ.CX1) THEN
!	If 2 of the starting values are equal, then quit
        IER=404
        RETURN
      ENDIF

      IER=-1
      CX=CX1
!	To evaluate the function at CX1
      RETURN
1000  CF1=CF
      IF1=IX
!	perform deflation
      DO 1500 J=1,NZ
1500  CF1=CF1/(CX1-CZERO(J))

      IER=-2
      CX=CX2
!	To calculate the function value at CX2
      RETURN
2000  CF2=CF
      IF2=IX
!	perform deflation
      DO 2500 J=1,NZ
2500  CF2=CF2/(CX2-CZERO(J))

      CF1I=CF1*2.**(IF1-IF2)
      CH1=(CF2-CF1I)/(CX2-CX1)
      I=0
      DX1=ABS(CX3-CX2)
      CX=CX3
      IER=-3
!	To calculate the function value at CX3
      RETURN

!	Loop for the Muller iteration
3000  I=I+1
      CF3=CF
      IF3=IX
      IF(NZ.GT.0) THEN
!	perform deflation
        DO 3500 J=1,NZ
3500    CF3=CF3/(CX3-CZERO(J))
      ENDIF

      IER=0
      IF(CX3.EQ.CX1.OR.CX3.EQ.CX2) RETURN
      IF(CF3.EQ.0.0) RETURN
      CF2A=CF2*2.**(IF2-IF3)
      CH2=(CF3-CF2A)/(CX3-CX2)
      CH1A=CH1*2.**(IF2-IF3)
      C2=(CH2-CH1A)/(CX3-CX1)
      CG=CH2+(CX3-CX2)*C2
      CI=SQRT(CG*CG-4.*CF3*C2)
      CD=CG+CI
      CD1=CG-CI
      IF(ABS(CD1).GT.ABS(CD)) CD=CD1

      IF(CD.EQ.0.0) THEN
!	If the denominator is zero, then quit
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
      IF1=IF2
      IF2=IF3
      CH1=CH2
      DX=ABS(CDX)
      IF(I.GT.2.AND.DX.LT.MAX(REPS*ABS(CX3),AEPS)) RETURN
      IF(I.GT.5.AND.DX.LT.MAX(REPS0*ABS(CX3),AEPS)
     1   .AND.DX.GT.DX1) THEN
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

      IF(I.LT.NIT) THEN
        IER=-3
        CX=CX3
!	To calculate the function value at CX3 and continue the loop
        RETURN
      ENDIF

!	Iteration fails to converge
      IER=43
      IF(DX.LT.MAX(REPS0*ABS(CX3),AEPS)) RETURN
      IER=432
      END
