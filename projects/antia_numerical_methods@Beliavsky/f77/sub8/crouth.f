!	Solution of a system of linear equations using Crout's algorithm 
!            with iterative refinement
!
!	N : (input) Number of equations to be solved
!	NUM : (input) Number of different sets (each with N equations) of
!	         equations to be solved
!	A : (input) The matrix of coefficient of size LJ*N
!	         A(I,J) is the coefficient of x_J in Ith equation
!	B : (output) Array of size LJ*N containing triangular decomposition of matrix A
!	X : (input/output) The matrix containing right hand sides (size LJ*NUM)
!	        X(I,J) is the Ith element of Jth right hand side
!	     	at output it will contain the solutions
!	DET, IDET : (output) The determinant of the matrix = DET*2**IDET
!	INC : (output) Integer array of length N containing information about
!		interchanges performed during elimination
!	LJ : (input) First dimension of arrays A and X in calling program
!	REPS : (input) Required relative precision in solution
!	IER : (output) Error flag, IER=0 signifies successful execution
!		IER=11 implies that iterative refinement did not converge
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
!	WK : Scratch array of length (2*N+NUM)
!		at output WK(I) will contain estimated error in Ith solution
!
!	Required routines : CROUT

      SUBROUTINE CROUTH(N,NUM,A,B,X,DET,IDET,INC,LJ,REPS,IER,IFLG,WK)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NITR=10)
!	If the compiler supports higher precision use this statement
!	Otherwise the iterative refinement process may not converge
!	as it is necessary to use higher precision arithmetic while
!	calculating the residuals.
      REAL*16 D1,D2
      DIMENSION A(LJ,N),B(LJ,N),X(LJ,NUM),INC(N),WK(2*N+NUM)

      IF(N.LE.0.OR.N.GT.LJ) THEN
        IER=102
        RETURN
      ENDIF

      IF(IFLG.LE.1) THEN
!	Preserving the matrix for calculating the residuals
        DO 1000 I=1,N
          DO 1000 J=1,N
1000    B(J,I)=A(J,I)
        IFLG1=1
!	Perform LU decomposition using CROUT
        CALL CROUT(N,NUM,B,X,DET,IDET,INC,LJ,IER,IFLG1,WK)
        IF(IER.GT.100) RETURN
      ENDIF

      IF(IFLG.EQ.1) THEN
        IFLG=2
        RETURN
      ENDIF
      IFLG=2
      NUM1=1
      IER=0

!	Solving the systems with NUM different right-hand sides
      DO 5000 J=1,NUM
        DO 2000 I=1,N
          WK(I)=X(I,J)
!	Preserving the RHS for calculating residuals
          WK(I+N)=WK(I)
          X(I,J)=0.0
2000    CONTINUE

        RP1=0
!	The iterative refinement
        DO 3000 IT=1,NITR
          CALL CROUT(N,NUM1,B,WK,DET,IDET,INC,LJ,IER1,IFLG1,WK)
          R1=0.0
          R2=0.0
          DO 2200 I=1,N
            IF(ABS(WK(I)).GT.R1) R1=ABS(WK(I))
            X(I,J)=X(I,J)+WK(I)
            IF(ABS(X(I,J)).GT.R2) R2=ABS(X(I,J))
2200      CONTINUE

!	The error estimate
          IF(IT.EQ.2) WK(J+2*N)=R1
          IF(R2.EQ.0.0) GO TO 5000
          RP=R1/R2
          IF(RP.LT.REPS) GO TO 5000
          IF(RP.GT.RP1.AND.IT.GT.1) THEN
            IER=11
            GO TO 5000
          ENDIF
          RP1=RP

!	Calculating the residue
          DO 2600 I=1,N
            D1=WK(I+N)
            DO 2400 K=1,N
!              D1=D1-A(I,K)*X(K,J)
!	To force double length accumulation of residuals
             D2=A(I,K)
             D1=D1-D2*X(K,J)
2400        CONTINUE
            WK(I)=D1
2600      CONTINUE
3000    CONTINUE
        IER=11

5000  CONTINUE
!	The error estimates
      DO 6000 I=1,NUM
6000  WK(I)=WK(I+2*N)
      END
