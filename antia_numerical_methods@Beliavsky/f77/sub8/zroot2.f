!	Complex zeros of a given function using Muller iteration with deflation
!	For use with MULER2 which requires function in form CF*2**IDET
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
!		FUNCTION CF(CZ,JF) must be supplied by the user, where CF and CZ
!		are complex. To avoid overflow and underflow the function is
!		calculated in the form CF*2**JF.
!
!	Required routines : MULER2, CF

      SUBROUTINE ZROOT2(N,CX,CZERO,NZ,REPS,AEPS,IER,RMAX,CF)
      IMPLICIT COMPLEX*16(C)
      IMPLICIT REAL*8(A,B,D-H,O-Z)
      EXTERNAL CF
      DIMENSION CZERO(NZ+N),CX(N)

      IF(N.LE.0) RETURN
      IF(NZ.LT.0) NZ=0
      IER=0

      DO 1000 I=1,N
!	The starting values for Muller's iteration
        CX3=CX(I)
        CDX=0.01D0*CX3
        IF(I.LT.N.AND.CDX.EQ.0.0) CDX=0.1D0*(CX(I+1)-CX(I))
        IF(I.GT.1.AND.CDX.EQ.0.0) CDX=0.1D0*(CX(I-1)-CX(I))
        IF(CDX.EQ.0.0) CDX=1.D-3
!	These values may need to be changed in some cases
        CX2=CX3+CDX
        CX1=CX3-CDX

!	Find the next zero
!        CALL MULLER(CX1,CX2,CX3,REPS,AEPS,IER1,CF,NZ,CZERO,RMAX)
!	For MULER use the preceding statements instead
        IER1=0
500     CALL MULER2(CX1,CX2,CX3,REPS,AEPS,IER1,CF0,CX0,IX,NZ,CZERO,RMAX)
        IF(IER1.LT.0) THEN
          CF0=CF(CX0,IX)
          GO TO 500
        ENDIF

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
