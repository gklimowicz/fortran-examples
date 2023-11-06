!	Roots of a polynomial with complex coefficients using Laguerre iteration
!
!	N : (input) The degree of the polynomial
!	COF : (input) Complex array of length N+1 containing the coefficients of
!		the polynomial. COF(I) is the coefficient of X**(I-1)
!	CX : (output) Complex array of length N, containing the computed roots
!	IER : (output) Error parameter, IER=0 for successful execution
!		IER=k*11 implies that iteration for refining the roots failed
!			to converge for k of the roots
!		IER=406 implies that N.LE.0 and no calculations are done
!		IER=408 implies that A(N+1)=0 and no calculations are done
!		IER=430 implies that iteration failed to converge for some root
!	QREFIN : (input) Logical parameter to decide if roots need to be refined
!		If QREFIN=.TRUE. the roots are refined using original polynomial
!		otherwise no refinement is done.
!	CWK : Complex array of length N+1 used as a scratch space
!
!	Required routines : LAGITC
!
      SUBROUTINE POLYC(N,COF,CX,IER,QREFIN,CWK)
      IMPLICIT COMPLEX*32(C)
      IMPLICIT LOGICAL(Q)
      IMPLICIT REAL*16(A,B,D-H,O,P,R-Z)
      DIMENSION COF(N+1),CX(N),CWK(N+1)

      IF(N.LE.0) THEN
        IER=406
        RETURN
      ENDIF
      IF(COF(N+1).EQ.0.0) THEN
        IER=408
        RETURN
      ENDIF

      DO 1000 I=1,N+1
1000  CWK(I)=COF(I)
      NP=N
      CXR=0.0
      IER=0

!	Find the next root
2000  CALL LAGITC(NP,CWK,CXR,IER1)
      IF(IER1.NE.0) THEN
!	If iteration fails to converge, try once more
        CXR=1.123456Q0
        CALL LAGITC(NP,CWK,CXR,IER1)
        IF(IER1.NE.0) THEN
!	If iteration fails again, then quit
          IER=430
          RETURN
        ENDIF
      ENDIF

      CXRT=CXR
      IF(NP.LT.N.AND.QREFIN) THEN
!	Refine the roots with original polynomial
        CALL LAGITC(N,COF,CXRT,IER1)
        IF(IER1.NE.0) THEN
!	If iteration fails to converge then retain the old value
          IER=IER+11
          CXRT=CXR
        ENDIF
      ENDIF

      IF(NP.GT.1) THEN
!	Perform deflation using unrefined root
        CN0=CWK(NP)
        CWK(NP)=CWK(NP+1)
        DO 2400 I=NP-1,1,-1
          CN=CXR*CWK(I+1)+CN0
          CN0=CWK(I)
          CWK(I)=CN
2400    CONTINUE
      ENDIF

      NP=NP-1
      CX(N-NP)=CXRT
!	If any more roots are left, find them
      IF(NP.GT.0) GO TO 2000
      END
