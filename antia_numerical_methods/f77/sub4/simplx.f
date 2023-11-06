!	To solve a linear programming problem using the simplex method
!
!	A : (input/output) Real array of length IA*(N+M2+1) containing
!		the tableau of simplex algorithm
!		A(1,I+1)=c_i, the cost coefficients
!		Rows 2 to M1+1 contain constraint of first (see below) type
!		Rows M1+2 to M1+M2+1 contain constraint of second type
!		Rows M1+M2+2 to M1+M2+M3+1 contain constraint of third type
!		For all constraints A(j,1)=b_j and A(j,i+1)=a_i
!	IA : (input) First dimension of array A as declared in the calling
!		program. IA .GE. M1+M2+M3+2
!	N : (input) Number of variables, each is constrained to be .GE.0
!	M1 : (input) Number of constraints of form a^T X .LE. b_i .GE. 0
!	M2 : (input) Number of constraints of form a^T X .GE. b_i .GE. 0
!	M3 : (input) Number of constraints of form a^T X = b_i .GE. 0
!	X : (output) Real array of length N, containing the optimal feasible
!		vector
!	F : (output) Optimum value of the objective function
!	IWK : Integer array of length N+M1+2M2+M3+1 used as scratch space
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=57 implies that the objective function is unbounded from below
!		IER=58 implies that there is no basic feasible vector
!		IER=505 implies that N,M1,M2 or M3 is negative or 
!			IA < M1+M2+M3+2, no calculations are done
!		IER=506 implies that some coefficient b_i for constraints
!			is negative, no calculations are done
!		IER=531 implies that the simplex algorithm failed to find
!			the optimal feasible vector
!		IER=532 implies that a basic feasible vector to start the
!			calculations could not be found
!	AEPS : (input) Required accuracy, any coefficient <AEPS, may be
!		assumed to be zero
!
!	Required routines : SIMPX
!
      SUBROUTINE SIMPLX(A,IA,N,M1,M2,M3,X,F,IWK,IER,AEPS)
      IMPLICIT LOGICAL(Q)
!      IMPLICIT REAL*8(A-H,O,P,R-Z)
      DIMENSION A(IA,N+M2+1),X(N),IWK(N+M1+2*M2+M3+1)

      IF(MIN(M1,M2,M3,N).LT.0.OR.IA.LT.M1+M2+M3+2) THEN
        IER=505
        RETURN
      ENDIF

      IER=0
!	Testing for negative b_i
      M=M1+M2+M3
      DO 1500 I=2,M+1
        IF(A(I,1).LT.0.0) IER=506
1500  CONTINUE
      IF(IER.GT.0) RETURN

      A(1,1)=0.0
      IF(M2.GT.0) THEN
!	Add extra artificial variables
        DO 2200 J=N+2,N+M2+1
          DO 2000 I=1,M+1
2000      A(I,J)=0.0
2200    A(M1+J-N,J)=-1.
      ENDIF
      N1=N+M+M2
      NM=N+M2+1
      NV=N+M1+M2
      DO 2300 I=1,NM-1
2300  IWK(I+M+1)=I
      DO 2400 I=1,M
2400  IWK(I+1)=NM-1+I

      IF(M2+M3.GT.0) THEN
!	Cost coefficients for the auxiliary objective function
        DO 2800 I=1,NM
          S=0.0
          DO 2600 J=M1+2,M+1
2600      S=S+A(J,I)
2800    A(M+2,I)=-S
        QF=.FALSE.

!	solve the auxiliary problem to get a basic feasible vector
        CALL SIMPX(A,IA,N1,M,NV,QF,IWK,IWK(M+1),IER,AEPS)
        IF(IER.GT.0) RETURN
        IF(A(M+2,1).LT.-AEPS) THEN
!	If the minimum is positive, then no basic feasible vector exists
          IER=58
          RETURN
        ENDIF

!	Remove the artificial variables from the tableau
        DO 4000 I=2,NM
          IF(IWK(I+M).GT.NV.AND.I.LE.N1-M+1) THEN
            IS=1
            DO 3200 J=I+M+1,N1
              IF(IWK(J).LE.NV) GO TO 3300
              IS=IS+1
3200        CONTINUE
3300        N1=N1-IS
            IF(N1-M+1.GE.I) THEN
              DO 3600 J=I,N1-M+1
                IWK(J+M)=IWK(J+M+IS)
                DO 3600 K=1,M+1
3600          A(K,J)=A(K,J+IS)
            ENDIF
          ENDIF
4000    CONTINUE

        IF(N1.NE.NV) THEN
!	If artificial variables are not removed, then quit
          IER=532
          RETURN
        ENDIF

      ENDIF
      QF=.TRUE.
      CALL SIMPX(A,IA,N1,M,NV,QF,IWK,IWK(M+1),IER,AEPS)
!	The minimum value of the objective function
      F=-A(1,1)
      DO 5200 I=2,M+1
        IF(IWK(I).LE.N) X(IWK(I))=A(I,1)
5200  CONTINUE
      DO 5400 I=M+2,N1
        IF(IWK(I).LE.N) X(IWK(I))=0.0
5400  CONTINUE

      END
