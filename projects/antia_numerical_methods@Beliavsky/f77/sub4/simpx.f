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
!      IMPLICIT REAL*8(A-H,O,P,R-Z)
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
