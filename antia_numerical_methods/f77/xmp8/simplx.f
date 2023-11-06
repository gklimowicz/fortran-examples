!     PROGRAM TO SOLVE A LINEAR PROGRAMMING PROBLEM USING SIMPLEX ALGORITHM

      PROGRAM LP
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(20,40),X(20),IWK(50)

51    FORMAT('  N =',I5,5X,'M1 =',I5,5X,'M2 =',I5,5X,'M3 =',I5)
52    FORMAT('  IER =',I4,5X,'MINIMUM =',1PD14.6,5X/
     1       '  OPTIMAL FEASIBLE VECTOR :',3D14.6/(2X,5D14.6))
53    FORMAT('   COST COEFFICIENTS :',1P4D14.6/(2X,5D14.6))
54    FORMAT('    CONSTRAINTS')
55    FORMAT(2X,1PD14.6,2X,4D14.6/(2X,5D14.6))

      REPS=1.D-7
      IA=20

100   PRINT *,'TYPE N=NO. OF VARIABLES     (QUITS WHEN N.LE.0)'
      READ *,N
      IF(N.LE.0) STOP
      PRINT *,'TYPE M1,M2,M3=NO. OF CONSTRAINTS OF <, > AND = TYPES'
      READ *,M1,M2,M3
      WRITE(6,51) N,M1,M2,M3
      PRINT *,'TYPE COST COEFFICIENTS'
      READ *,(A(1,J),J=2,N+1)
      WRITE(6,53) (A(1,J),J=2,N+1)

      WRITE(6,54)
      M=M1+M2+M3
      PRINT *,'TYPE COEFFICIENTS OF CONSTRAINTS, WITH RHS FIRST'
      DO 500 I=1,M
        PRINT *,I,'TH CONSTRAINT'
        READ *,(A(I+1,J),J=1,N+1)
        WRITE(6,55) (A(I+1,J),J=1,N+1)
500   CONTINUE
      CALL SIMPLX(A,IA,N,M1,M2,M3,X,F,IWK,IER,REPS)
      WRITE(6,52) IER,F,(X(I),I=1,N)
      GO TO 100
      END

!     ----------------------------------------------

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
      IMPLICIT REAL*8(A-H,O,P,R-Z)
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
      IMPLICIT REAL*8(A-H,O,P,R-Z)
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
