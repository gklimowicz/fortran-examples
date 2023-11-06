!     PROGRAM TO FIND RATIONAL FUNCTION MINIMAX APPROXIMATION TO DISCRETE DATA

      PROGRAM MINMX
      IMPLICIT REAL*16(A-H,O-Z)
      DIMENSION B(20),WK(7000),IWK(321),F(100),X(100)

!     EXERCISE 10.57 : GENERATING MINIMAX APPROXIMATION OF ARC TAN(X) 
      FM(Y)=ATAN(Y)

51    FORMAT('   IER =',I4,5X,'MAXIMUM ERROR =',1PD14.6/
     1       '  COEF. IN NUMERATOR :',3D16.8/(2X,5D16.8))
52    FORMAT('  COEF. IN DENOMINATOR :',1P3D16.8/(2X,5D16.8))
53    FORMAT('   M =',I3,5X,'K =',I3,5X,'N =',I3)
54    FORMAT('   INITIAL GUESS FOR COEF. :',1P3D14.6/(2X,5D14.6))

100   PRINT *,'TYPE M=DEGREE OF NUMERATOR,  K=DEGREE OF DENOMINATOR'
      READ *,M,K
      PRINT *,'TYPE N=NO. OF DATA PTS OVER [0,1]   (QUITS WHEN N.LE.0)'
      READ *,N
      IF(N.LE.0) STOP

!     THE ITERATION MAY NOT CONVERGE UNLESS INITIAL APPROXIMATION IS
!     SUFFICIENTLY CLOSE TO THE ACTUAL VALUE

      PRINT *,'TYPE INITIAL GUESS FOR COEF'
      READ *,(B(I),I=1,M+K+1)
      WRITE(6,53) M,K,N
      WRITE(6,54) (B(I),I=1,M+K+1)

!     GENERATING THE INPUT DATA SET 

      H=2.Q0/(N-1.)
      DO 2000 I=1,N
        XI=-1+(I-1)*H
        X(I)=XI
        F(I)=FM(XI)
2000  CONTINUE

      EPS=1.Q-6
      CALL MINMAX(M,K,N,B,X,F,EPS,EMAX,IER,WK,IWK)
      WRITE(6,51) IER,EMAX,(B(I),I=K+1,M+K+1)
      WRITE(6,52) (B(I),I=1,K)
      GO TO 100
      END

!     -------------------------------------------

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
      IMPLICIT REAL*16(A-H,O,P,R-Z)
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

!     --------------------------------------------

!	To solve a linear programming problem in the standard form
!		 using the simplex method
!
!	A : (input/output) Real array of length IA*(N-M+1) containing
!		the tableau of simplex algorithm
!		A(1,I+1)=c_i, the cost coefficients
!		Rows 2 to M+1 contain constraints with A(j,1)=b_j and A(j,i+1)=a_i
!		Row M+2 contains the cost coefficients for auxiliary problem
!		when QF=.FALSE.
!	IA : (input) First dimension of array A as declared in the calling
!		program. IA .GE. M+2
!	N : (input) Number of variables, each is constrained to be .GE.0
!	M : (input) Number of constraints of form a^T X = b_i .GE. 0
!	NV : (input) Number of variables, excluding the artificial variables
!	QF : (input) Logical variable to decide which objective function
!		is to be minimised.
!		QF=.TRUE. if the main objective function specified by
!			the first row of A is to be minimised.
!		QF=.FALSE. if the auxiliary objective function specified by
!			the last row of A is to be minimised.
!	ID : (input/output) integer array of length M+1 which contains
!		information about interchange of variables on LHS
!	IV : (input/output) integer array of length N-M+1 which contains
!		information about interchange of variables on RHS
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=57 implies that the objective function is unbounded from below
!		IER=531 implies that the simplex algorithm failed to find
!			the optimal feasible vector
!	AEPS : (input) Required accuracy, any coefficient <AEPS, may be
!		assumed to be zero
!
!	Required routines : None

      SUBROUTINE SIMPX(A,IA,N,M,NV,QF,ID,IV,IER,AEPS)
      IMPLICIT LOGICAL(Q)
      IMPLICIT REAL*16(A-H,O,P,R-Z)
      DIMENSION A(IA,N-M+1),ID(M+1),IV(N-M+1)
      PARAMETER(NIT=20)

      IF(QF) THEN
!	Minimise the objective function in the first row
        JF=1
        M1=M+1
      ELSE
!	Minimise the objective function in the last row
        JF=M+2
        M1=JF
      ENDIF
      N1=N-M+1
      IER=0

      DO 4000 IT=1,NIT*(N+M)
!	Find the minimum of the reduced cost coefficients
        RMIN=A(JF,2)
        K=2
        DO 2000 J=3,N1
          IF(A(JF,J).LT.RMIN) THEN
            RMIN=A(JF,J)
            K=J
          ENDIF
2000    CONTINUE

        IF(RMIN.GE.0.0) THEN
!	The objective function cannot be reduced further
          IF(QF.OR.A(JF,1).LT.-AEPS) RETURN
!	Check if any artificial variable is on the LHS
          DO 2300 I=2,M+1
            IF(ID(I).GT.NV.AND.ABS(A(I,1)).LE.AEPS) THEN
              A(I,1)=0.0
              L=I
              K=0
              RMAX=0.0
              DO 2200 J=2,N1
                IF(ABS(A(I,J)).GT.RMAX.AND.IV(J).LE.NV) THEN
                  RMAX=ABS(A(I,J))
                  K=J
                ENDIF
2200          CONTINUE
!	To exchange the artificial variable
              IF(K.GT.0) GO TO 2600
            ENDIF
2300      CONTINUE
          RETURN
        ENDIF

!	Finding allowed change in the Kth variable
        RMIN=0.0
        L=0
        DO 2400 J=2,M+1
          IF(A(J,K).GT.AEPS) THEN
            R1=A(J,1)/A(J,K)
            IF(R1.LT.RMIN.OR.L.EQ.0) THEN
              RMIN=R1
              L=J
            ENDIF
          ENDIF
2400    CONTINUE
        IF(L.EQ.0) THEN
!	The objective function is unbounded from below
          IER=57
          RETURN
        ENDIF

!	exchange the variables
2600    L1=ID(L)
        ID(L)=IV(K)
        IV(K)=L1
        DO 3000 J=1,N1
          IF(J.NE.K) THEN
            DO 2800 I=1,M1
              IF(I.NE.L) A(I,J)=A(I,J)-A(I,K)*A(L,J)/A(L,K)
2800        CONTINUE
          ENDIF
3000    CONTINUE

        DO 3200 J=1,N1
          IF(J.NE.K) A(L,J)=A(L,J)/A(L,K)
3200    CONTINUE
        DO 3400 I=1,M1
          IF(I.NE.L) A(I,K)=-A(I,K)/A(L,K)
3400    CONTINUE
        A(L,K)=1./A(L,K)
4000  CONTINUE

!	Iteration fails to converge
      IER=531
      END
