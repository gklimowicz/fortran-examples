!	To solve a linear programming problem in the standard form arising
!		 in L_1 minimisation problems using the simplex method
!
!	A : (input/output) Real array of length IA*(N-M+1) containing
!		the tableau of simplex algorithm
!		A(1,I+1)=c_i, the cost coefficients
!		Rows 2 to M+1 contain constraints with A(j,1)=b_j and A(j,i+1)=a_i
!	IA : (input) First dimension of array A as declared in the calling
!		program. IA .GE. M+2
!	N : (input) Number of variables, each is constraint to be .GE.0
!	M : (input) Number of constraints of form a^T X = b_i .GE. 0
!	ID : (input/output) integer array of length M+1 which contains
!		information about interchange of variables on LHS
!	IV : (input/output) integer array of length N-M+1 which contains
!		information about interchange of variables on RHS
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=63 implies that the objective function is unbounded from below
!		IER=635 implies that the simplex algorithm failed to find
!			the optimal feasible vector
!	AEPS : (input) Required accuracy, any coefficient <AEPS, may be
!		assumed to be zero
!
!	Required routines : None
!
      SUBROUTINE SIMPL1(A,IA,N,M,ID,IV,IER,AEPS)
!      IMPLICIT REAL*8(A-H,O,P,R-Z)
      IMPLICIT LOGICAL(Q)
      DIMENSION A(IA,N-M+1),ID(M+1),IV(N-M+1)
      PARAMETER(NIT=20)

      JF=1
      M1=M+1
      N1=N-M+1
      IER=0

!	The simplex iteration
      DO 4000 IT=1,NIT*(N+M)
!	Finding the minimum reduced cost coefficient
        RMIN=0.0
        K=0
        DO 2000 J=2,N1
          IF(A(JF,J).LT.RMIN) THEN
            RMIN=A(JF,J)
            K=J
          ELSE IF(IV(J).LE.M.AND.2-A(JF,J).LT.RMIN) THEN
            RMIN=2-A(JF,J)
            K=-J
          ENDIF
2000    CONTINUE
        IF(RMIN.GE.-AEPS) RETURN

!	Finding the pivot element
        K1=ABS(K)
        RMIN=0.0
        L=0
        DO 2400 J=2,M+1
          AJ=A(J,K1)
          IF(K.LT.0) AJ=-AJ
          IF(AJ.GT.AEPS) THEN
            R1=A(J,1)/AJ
            IF(R1.LT.RMIN.OR.L.EQ.0) THEN
              RMIN=R1
              L=J
            ENDIF
          ENDIF
2400    CONTINUE
        IF(L.EQ.0) THEN
!	The objective function is unbounded from below
          IER=63
          RETURN
        ENDIF

        IF(K.LT.0) THEN
          A(JF,K1)=-2+A(JF,K1)
          IV(K1)=-IV(K1)
        ENDIF

!	Exchange the variables
        L1=ID(L)
        ID(L)=IV(K1)
        IV(K1)=L1
        DO 3000 J=1,N1
          IF(J.NE.K1) THEN
            R1=A(L,J)/A(L,K1)
            DO 2800 I=1,M1
              IF(I.NE.L) A(I,J)=A(I,J)-A(I,K1)*R1
2800        CONTINUE
          ENDIF
3000    CONTINUE

        R1=ABS(A(L,K1))
        DO 3200 J=1,N1
          IF(J.NE.K1) A(L,J)=A(L,J)/R1
3200    CONTINUE
        DO 3400 I=1,M1
          IF(I.NE.L) A(I,K1)=-A(I,K1)/A(L,K1)
3400    CONTINUE
        A(L,K1)=1./R1
4000  CONTINUE

!	Iteration fails to converge
      IER=635
      END
