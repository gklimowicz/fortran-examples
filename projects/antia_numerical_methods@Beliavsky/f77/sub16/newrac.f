!	Complex zero of a given function using Newton-Raphson iteration
!
!	CX : (input/output) Initial guess for the zero, after execution
!		it will contain the computed zero.
!	REPS : (input) Required relative accuracy
!	AEPS : (input) Required absolute accuracy
!		The estimated error should be less than MAX(AEPS, REPS*ABS(CX))
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=425 implies that iteration failed to converge to specified accuracy
!		IER=426 implies that derivative is zero and the iteration
!			has to be abandoned at this stage.
!	CFUN : (input) Name of the function routine to calculate the function
!		FUNCTION CFUN(CX,CDX) must be supplied by the user.
!		Here CDX is the first derivative of CFUN.
!
!	Required routines : CFUN

      SUBROUTINE NEWRAC(CX,REPS,AEPS,IER,CFUN)
      IMPLICIT COMPLEX*32(C)
      IMPLICIT REAL*16(A,B,D-H,O-Z)
      PARAMETER(NIT=75)

      IER=0
      DO 1000 L=1,NIT
        CF=CFUN(CX,CDF)
        IF(CDF.EQ.0.0) THEN
          IF(CF.EQ.0.0) RETURN
          IER=426
          RETURN
        ENDIF

!	Newton-Raphson iteration
        CDX=-CF/CDF
        CX=CX+CDX
        IF(ABS(CDX).LT.MAX(REPS*ABS(CX),AEPS).AND.L.GT.2) RETURN
1000  CONTINUE

      IER=425
      END
