!	Root of a polynomial with real coefficients using Laguerre iteration
!
!	N : (input) The degree of the polynomial
!	A : (input) Real array of length N+1 containing the coefficients of
!		the polynomial. A(1) is the constant term and A(N+1) is the
!		coefficient of X**N
!	CXI : (input/output) Complex variable containing the initial guess,
!		 after execution it will contain the computed root
!	IER : (output) Error parameter, IER=0 for successful execution
!		IER=438 implies that denominator is zero and iteration cannot
!			be continued further
!		IER=439 implies that iteration has failed to converge
!
!	Required routines : None

      SUBROUTINE LAGITR(N,A,CXI,IER)
      IMPLICIT COMPLEX*32(C)
      IMPLICIT LOGICAL(Q)
      IMPLICIT REAL*16(A,B,D-H,O,P,R-Z)
      PARAMETER(ITMAX=50,REPS=1.Q-4,AEPS=1.Q-6)
      DIMENSION A(N+1)

      IER=0
      CDX1=ABS(CXI)+1.0
      QC=.FALSE.
      IC=ITMAX

!	The Laguerre iteration
      DO 2000 I=1,ITMAX
!	Calculate the polynomial and its derivatives
        CF=A(N+1)
        CFP=0.0
        CFPP=0.0
        DO 1000 J=N,1,-1
          CFPP=CXI*CFPP+2.*CFP
          CFP=CXI*CFP+CF
          CF=CXI*CF+A(J)
1000    CONTINUE

        CH=(N-1)*((N-1)*CFP*CFP-N*CF*CFPP)
        CH=SQRT(CH)
        CDEN=CFP+CH
        IF(ABS(CFP-CH).GT.ABS(CDEN)) CDEN=CFP-CH

        IF(CDEN.NE.0.0) THEN
!	Laguerre's iteration
          CDX=-N*CF/CDEN
          IF(ABS(CDX).LT.MAX(REPS*ABS(CXI),AEPS).AND.I.GT.1.AND.
     1        (.NOT.QC)) THEN
            QC=.TRUE.
            IC=I
          ENDIF
          IF(QC.AND.ABS(CDX/CDX1).GT.1.0) RETURN
          IF(I-IC.GT.5.AND.ABS(CDX/CDX1).GT.0.99Q0) RETURN
          CDX1=CDX
          IF(CDX.EQ.0.0) RETURN
          CXI=CXI+CDX
        ELSE
          IF(CF.EQ.0.0) RETURN
!	If denominator vanishes, then quit
          IER=438
          RETURN
        ENDIF
2000  CONTINUE

!	Iteration fails to converge
      IF(.NOT.QC) IER=439
      END
