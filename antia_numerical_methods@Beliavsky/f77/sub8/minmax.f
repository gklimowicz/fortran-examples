!	To calculate coefficients of Rational function minimax approximation
!	for a tabulated function
!
!	M : (input) Required degree of polynomial in the numerator
!	K : (input) Required degree of polynomial in the denominator
!	N : (input) Number of points in the table of values. N> M+K+1
!	A : (input/output) Real array of length M+K+2 containing the coefficients
!		of approximation. A(I) is the coefficient of x**I in
!		the denominator, the constant term being 1. A(K+J+1) is
!		the coefficient of x**J in the numerator. The 
!		initial guess for coefficients must be supplied.
!	X : (input) Real array of length N containing the points at which
!		function value is available
!	F : (input) Real array of length N containing the function values at X(I)
!	EPS : (input) Required accuracy, the iteration is continued until
!		change in EMAX is less than EPS
!	EMAX : (output) Maximum error in the calculated approximation
!	IER : (output) error parameter, IER=0 implies successful execution
!		IER=615 implies that M<0, K<0 or N<M+K+2
!			in this case no calculations are done
!		IER=634 implies that the iteration did not converge
!			to the specified accuracy
!		Other values of IER may be set by SIMPX
!	WK : Real array of length (K+M+5)*(3N+1) used as scratch space
!	IWK : Integer array of length 3N+K+M+3 used as scratch space
!
!	Required routines : SIMPX
!
      SUBROUTINE MINMAX(M,K,N,A,X,F,EPS,EMAX,IER,WK,IWK)
      IMPLICIT REAL*8(A-H,O,P,R-Z)
      IMPLICIT LOGICAL(Q)
      PARAMETER(NIT=30)
      DIMENSION A(M+K+2),F(N),IWK(3*N+M+K+3),WK(M+K+5,3*N+1),X(N)

      NK=M+K+2
      IF(M.LT.0.OR.K.LT.0.OR.N.LT.NK) THEN
        IER=615
        RETURN
      ENDIF

      IER=0
      LJ=M+K+5
      NV=3*N
!	The number of constraints
      M3=NK+1
      QF=.TRUE.
      NC=NV+M3
      N2=2*N

!	Loop for iteration
      DO 5000 IT=1,NIT

!	Finding the maximum error in approximation
        EMAX=0.0
        DO 2400 I=1,N
          XI=X(I)
          FD=0.0
          DO 2100 J=K,1,-1
2100      FD=FD*XI+A(J)
          FD=FD*XI+1
          WK(LJ,I)=FD

          FN=A(M+K+1)
          DO 2200 J=M+K,K+1,-1
2200      FN=FN*XI+A(J)
          EI=F(I)-FN/FD
          IF(ABS(EI).GT.EMAX) EMAX=ABS(EI)
2400    CONTINUE

!	Setting up the tableau for Simplex algorithm
        DO 3600 I=1,N
          XI=X(I)
          FD=WK(LJ,I)
          FI=F(I)

          IWK(I+M3+1)=I
          IWK(I+N+M3+1)=I+N
          IWK(I+N2+M3+1)=I+N2
          WK(1,I+1)=-FI+EMAX-EMAX*FD
          WK(1,I+N+1)=FI+EMAX-EMAX*FD
          WK(1,I+N2+1)=1
          WK(2,I+1)=FD
          WK(2,I+N+1)=FD
          WK(2,I+N2+1)=0.0
          S3=0.0
          DO 2800 J=1,K
            TI=XI**J
            S3=S3+TI
            WK(2+J,I+1)=TI*(EMAX-FI)
            WK(2+J,I+N+1)=TI*(EMAX+FI)
            WK(2+J,I+N2+1)=TI
2800      CONTINUE
          S1=S3*(EMAX-FI)+1
          S2=S3*(EMAX+FI)-1
          WK(3+K,I+1)=1
          WK(3+K,I+N+1)=-1
          WK(3+K,I+N2+1)=0.0
          DO 3000 J=1,M
            TI=XI**J
            S1=S1+TI
            S2=S2-TI
            WK(3+K+J,I+1)=TI
            WK(3+K+J,I+N+1)=-TI
            WK(3+K+J,I+N2+1)=0.0
3000      CONTINUE
          WK(M+K+4,I+1)=-S1
          WK(M+K+4,I+N+1)=-S2
          WK(M+K+4,I+N2+1)=-S3
3600    CONTINUE

        WK(1,1)=0
        WK(2,1)=1
        IWK(2)=3*N+1
        DO 3200 J=3,M+K+4
          IWK(J)=3*N+J-1
3200    WK(J,1)=0.0

        CALL SIMPX(WK,LJ,NC,M3,NV,QF,IWK,IWK(M3+1),IER,EPS)
        IF(IER.GT.0) RETURN

!	The maximum error
        EI=WK(1,1)
!	Obtaining the coefficients from the tableau
        DO 4000 I=M3+2,NC+1
          IF(IWK(I).GT.NV+1) A(IWK(I)-NV-1)=WK(1,I-M3)
4000    CONTINUE
        DO 4200 I=2,M3+1
          IF(IWK(I).GT.NV+1) A(IWK(I)-NV-1)=0.0
4200    CONTINUE
        DO 4400 I=1,M+K+1
4400    A(I)=A(I)-A(M+K+2)
        DIF=EMAX-EI
        EMAX=EI

!	The convergence test
        IF(DIF.LT.EPS.OR.K.EQ.0.OR.IER.GT.0) RETURN
5000  CONTINUE
      IER=634
      END
