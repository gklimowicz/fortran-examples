!	Solution of a system of linear equations using Crout's algorithm 
!	with partial pivoting
!
!	N : (input) Number of equations to be solved
!	NUM : (input) Number of different sets (each with N equations) of
!	         equations to be solved
!	A : (input/output) The matrix of coefficient of size LJ*N
!	         A(I,J) is the coefficient of x_J in Ith equation
!	     at output it will contain the triangular decomposition
!	X : (input/output) The matrix containing right hand sides (size LJ*NUM)
!	        X(I,J) is the Ith element of Jth right hand side
!	     	at output it will contain the solutions
!	DET, IDET : (output) The determinant of the matrix = DET*2**IDET
!	INC : (output) Integer array of length N containing information about
!		interchanges performed during elimination
!	LJ : (input) First dimension of arrays A and X in calling program
!	IER : (output) Error flag, IER=0 signifies successful execution
!		IER=102 implies (N.LE.0 or N.GT.LJ) 
!		IER=122 implies some pivot turned out to be zero and hence
!			matrix must be nearly singular
!	IFLG : (input) Integer parameter which determines the type of computation
!		required.
!		If IFLG.LE.0, both elimination and solution are calculated
!			and IFLG is set to 2
!		If IFLG=1, only elimination is done and IFLG is set to 2
!		If IFLG.GE.2 only solution is calculated, the triangular
!		    	decomposition should have been calculated earlier
!	WK : scratch array of length N
!
!	Required routines : None

      SUBROUTINE CROUT(N,NUM,A,X,DET,IDET,INC,LJ,IER,IFLG,WK)
      IMPLICIT REAL*8(A-H,O-Z)

      DIMENSION A(LJ,N),INC(N),X(LJ,NUM),WK(N)

      IF(N.LE.0.OR.N.GT.LJ) THEN
        IER=102
        RETURN
      ENDIF
      IER=122

      IF(IFLG.LE.1) THEN
!	Perform LU decomposition
        DO 1500 I=1,N
          R1=0.0
          DO 1200 J=1,N
            IF(ABS(A(I,J)).GT.R1) R1=ABS(A(I,J))
1200      CONTINUE
          WK(I)=R1
!	If any row is zero, then quit
          IF(R1.EQ.0.0) RETURN
1500    CONTINUE

        DET=1.0
        IDET=0
        DO 2600 K=1,N
          R1=0.0
          KM=K
!	Generate the Kth column of L
          DO 2200 L=K,N
            D1=A(L,K)
            DO 2100 L1=1,K-1
2100        D1=D1-A(L,L1)*A(L1,K)
            A(L,K)=D1

!	Finding the pivot
            R2=ABS(D1/WK(L))
            IF(R2.GT.R1) THEN
              R1=R2
              KM=L
            ENDIF
2200      CONTINUE

          INC(K)=KM
!	Interchange the rows if needed
          IF(KM.NE.K) THEN
            DET=-DET
            DO 2300 L=1,N
              T1=A(K,L)
              A(K,L)=A(KM,L)
2300        A(KM,L)=T1
            T1=WK(K)
            WK(K)=WK(KM)
            WK(KM)=T1
          ENDIF

          DET=DET*A(K,K)
!	If the pivot is zero, then quit
          IF(A(K,K).EQ.0.0) RETURN
!	To check for singular or nearly singular matrices replace this statement by
!         IF(ABS(A(K,K)).LT.REPS) RETURN

          IF(DET.NE.0.0) THEN
!	Scale the value of the determinant DET
2350        IF(ABS(DET).GT.32.) THEN
              DET=DET*0.03125D0
              IDET=IDET+5
              GO TO 2350
            ENDIF

2370        IF(ABS(DET).LT.0.03125D0) THEN
              DET=DET*32.
              IDET=IDET-5
              GO TO 2370
            ENDIF
          ENDIF

!	Generate the Kth row of U
          DO 2500 L=K+1,N
            D1=A(K,L)
            DO 2400 L1=1,K-1
2400        D1=D1-A(K,L1)*A(L1,L)
            A(K,L)=D1/A(K,K)
2500      CONTINUE
2600    CONTINUE
        IER=0

        IF(IFLG.EQ.1) THEN
          IFLG=2
          RETURN
        ENDIF
        IFLG=2
      ENDIF

      IER=0
!	Solution for NUM different right-hand sides
      DO 5000 J=1,NUM
!	Forward substitution
        DO 3000 K=1,N
          IF(K.NE.INC(K)) THEN
            T1=X(K,J)
            X(K,J)=X(INC(K),J)
            X(INC(K),J)=T1
          ENDIF

          D1=X(K,J)
          DO 2800 L=1,K-1
2800      D1=D1-A(K,L)*X(L,J)
          X(K,J)=D1/A(K,K)
3000    CONTINUE

!	Back-substitution
        DO 3300 K=N-1,1,-1
          D1=X(K,J)
          DO 3200 L=N,K+1,-1
3200      D1=D1-X(L,J)*A(K,L)
3300    X(K,J)=D1
5000  CONTINUE
      END
