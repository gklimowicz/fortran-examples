!     PROGRAM TO SOLVE A SYSTEM OF LINEAR EQUATIONS USING GAUSSIAN ELIMINATION
!     OR CROUT'S DECOMPOSITION OR CROUT'S DECOMPOSITION WITH ITERATIVE REFINEMENT
!     OR CHOLESKY'S DECOMPOSITION FOR SYMMETRIC POSITIVE DEFINITE MATRIX
!     OR SINGULAR VALUE DECOMPOSITION

!     CHOLSK ASSUMES THE MATRIX TO BE POSITIVE DEFINITE AND SYMMETRIC
!     OTHER ROUTINES WILL WORK FOR A GENERAL MATRIX

      PROGRAM LINEAR
!      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(20),WK(60),AM(20,20),INT(20),BM(20,20),SIGMA(20)

!     EXAMPLE 3.5 : HILBERT MATRIX

51    FORMAT('GAUELM :   IER =',I4,5X,'N =',I3,4X,'DET =',1PE14.6)
52    FORMAT('CROUT :   IER =',I4,5X,'N =',I3,4X,'DET =',1PE14.6,5X,
     1       'IDET =',I7)
53    FORMAT('CROUTH :   IER =',I4,4X,'N =',I3,4X,'DET =',1PE14.6,5X,
     1       'IDET =',I7/14X,'ERROR ESTIMATE =',E10.2)
54    FORMAT(' SOLUTION :',(1P5E14.6))
55    FORMAT('  IER =',I4,5X,'N =',I3,4X,'SINGULAR VALUES :',1P2E14.6/
     1       (1X,5E14.6))
56    FORMAT('CHOLSK :   IER =',I4,5X,'N =',I3,4X,'DET =',1PE14.6)

      NUM=1
      LJ=20

100   PRINT *,'TYPE N=NO. OF EQUATIONS       (QUITS WHEN N.LE.0)'
      READ *,N
      IF(N.LE.0) STOP
      PRINT *,'TYPE IT=1/2/3/4/5  FOR  GAUELM/CROUT/CROUTH/CHOLSK/SVD'
      READ *,IT

!     GENERATING THE N X N HILBERT MATRIX AND THE RIGHT HAND SIDE

      DO 1000 I=1,N
        R1=0.0
        DO 800 J=1,N
          AM(I,J)=1.D0/(I+J-1.)
          R1=R1+AM(I,J)
800     CONTINUE
        X(I)=R1
1000  CONTINUE

      IFLG=0
      IF(IT.EQ.1) THEN
        CALL GAUELM(N,NUM,AM,X,DET,INT,LJ,IER,IFLG)
        WRITE(6,51) IER,N,DET
      ELSE IF(IT.EQ.2) THEN
        CALL CROUT(N,NUM,AM,X,DET,IDET,INT,LJ,IER,IFLG,WK)
        WRITE(6,52) IER,N,DET,IDET
      ELSE IF(IT.EQ.3) THEN
        REPS=1.E-7
        CALL CROUTH(N,NUM,AM,BM,X,DET,IDET,INT,LJ,REPS,IER,IFLG,WK)
        WRITE(6,53) IER,N,DET,IDET,WK(1)
      ELSE IF(IT.EQ.4) THEN
        CALL CHOLSK(N,NUM,AM,X,DET,LJ,IER,IFLG)
        WRITE(6,56) IER,N,DET
      ELSE
        M=N
        LU=LJ
        REPS=1.E-7
        CALL SVD(N,M,AM,BM,SIGMA,LU,LJ,WK,IER)
        WRITE(6,55) IER,N,(SIGMA(I),I=1,N)
        CALL SVDEVL(N,M,AM,BM,SIGMA,LU,LJ,X,WK,REPS)
      ENDIF

      WRITE(6,54) (X(I),I=1,N)
      GO TO 100
      END

!     ---------------------------------------------

!     Solution of a system of linear equations with real symmetric
!            positive definite matrix using Cholesky decomposition
!
!     N : (input) Number of equations to be solved
!     NUM : (input) Number of different sets (each with N equations) of
!              equations to be solved
!     A : (input/output) The matrix of coefficient of size ND*N
!              	A(I,J) is the coefficient of x_J in Ith equation
!          	at output it will contain the triangular decomposition
!     X : (input/output) The matrix containing right hand sides (size ND*NUM)
!               X(I,J) is the Ith element of Jth right hand side
!          	at output it will contain the solutions
!     DET : (output) The determinant of the matrix
!     ND : (input) First dimension of arrays A and X in calling program
!     IER : (output) Error flag, IER=0 signifies successful execution
!     		IER=103 implies (N.LE.0 or N.GT.ND)
!     		IER=123 implies some pivot turned out to be zero
!     IFLG : (input) Integer parameter which determines the type of computation
!		required.
!		If IFLG.LE.0, both elimination and solution are calculated
!		    and IFLG is set to 2
!		If IFLG=1, only elimination is done and IFLG is set to 2
!		If IFLG.GE.2 only solution is calculated, the triangular
!		    decomposition should have been calculated earlier
!
!	Required routines : None
 
      SUBROUTINE CHOLSK(N,NUM,A,X,DET,ND,IER,IFLG)
!      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(ND,N),X(ND,NUM)
 
      IF(N.LE.0.OR.N.GT.ND) THEN
        IER=103
        RETURN
      ENDIF
 
      IER=123
      IF(IFLG.LE.1) THEN
!     Perform triangular decomposition
 
        DET=1.0
        DO 2000 K=1,N
          DO 1500 I=1,K-1
            IF(A(I,I).EQ.0.0) RETURN
            SUM=A(K,I)
            DO 1200 J=1,I-1
              SUM=SUM-A(I,J)*A(K,J)
1200        CONTINUE
            A(K,I)=SUM/A(I,I)
1500      CONTINUE
          SUM=A(K,K)
          DO 1800 J=1,K-1
            SUM=SUM-A(K,J)**2
1800      CONTINUE
          IF(SUM.LE.0.0) RETURN
          A(K,K)=SQRT(SUM)
          DET=DET*SUM
2000    CONTINUE
 
        IER=0
        IF(IFLG.EQ.1) THEN
          IFLG=2
          RETURN
        ENDIF
        IFLG=2
      ENDIF
 
 
      IER=0
!     Solution for the NUM different right-hand sides
      DO 5000 J=1,NUM
!     Forward substitution
        X(1,J)=X(1,J)/A(1,1)
        DO 3100 K=2,N
          DO 3000 L=1,K-1
3000      X(K,J)=X(K,J)-A(K,L)*X(L,J)
3100    X(K,J)=X(K,J)/A(K,K)
 
!     back-substitution
        X(N,J)=X(N,J)/A(N,N)
        DO 3300 K=N-1,1,-1
          DO 3200 L=N,K+1,-1
3200      X(K,J)=X(K,J)-X(L,J)*A(L,K)
3300    X(K,J)=X(K,J)/A(K,K)
5000  CONTINUE
      END
 

!     ------------------------------------------------------

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
!      IMPLICIT REAL*8(A-H,O-Z)

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

!     ----------------------------------------------

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
!      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NITR=10)
!	If the compiler supports higher precision use this statement
!	Otherwise the iterative refinement process may not converge
!	as it is necessary to use higher precision arithmetic while
!	calculating the residuals.
      REAL*8 D1,D2
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

!     ----------------------------------------------

!	Solution of a system of linear equations using Gaussian elimination
!	with partial pivoting
!
!	N : (input) Number of equations to be solved
!	NUM : (input) Number of different sets (each with N equations) of
!	        equations to be solved
!	A : (input/output) The matrix of coefficient of size LJ*N
!	        A(I,J) is the coefficient of x_J in Ith equation
!	     	at output it will contain the triangular decomposition
!	X : (input/output) The matrix containing right hand sides (size LJ*NUM)
!	        X(I,J) is the Ith element of Jth right hand side
!	     	at output it will contain the solutions
!	DET : (output) The determinant of the matrix
!	INC : (output) Integer array of length N containing information about
!		interchanges performed during elimination
!	LJ : (input) First dimension of arrays A and X in calling program
!	IER : (output) Error flag, IER=0 signifies successful execution
!		IER=101 implies (N.LE.0 or N.GT.LJ) 
!		IER=121 implies some pivot turned out to be zero and hence
!			matrix must be nearly singular
!	IFLG : (input) Integer parameter to specify the type of computation required
!		If IFLG.LE.0, both elimination and solution are
!			done and IFLG is set to 2
!		If IFLG=1, only elimination is done and IFLG is set to 2
!		If IFLG.GE.2 only solution is calculated, the triangular
!			decomposition should have been calculated earlier
!
!	Required routines : None

      SUBROUTINE GAUELM(N,NUM,A,X,DET,INC,LJ,IER,IFLG)
!      IMPLICIT REAL*8(A-H,O-Z)
!	For complex matrices use the following statements instead
!      IMPLICIT REAL*8(R)
!      IMPLICIT COMPLEX*16(A-H,S-Z)

      DIMENSION A(LJ,N),INC(N),X(LJ,NUM)

      IF(N.LE.0.OR.N.GT.LJ) THEN
        IER=101
        RETURN
      ENDIF

      IER=121
      IF(IFLG.LE.1) THEN
!	Perform elimination

        DET=1.0
        DO 2600 K=1,N-1
!	Find the maximum element in the Kth column
          R1=0.0
          KM=K
          DO 2200 L=K,N
            IF(ABS(A(L,K)).GT.R1) THEN
              R1=ABS(A(L,K))
              KM=L
            ENDIF
2200      CONTINUE

          INC(K)=KM
          IF(KM.NE.K) THEN
!	Interchange the rows if needed
            DO 2300 L=K,N
              T1=A(K,L)
              A(K,L)=A(KM,L)
2300        A(KM,L)=T1
            DET=-DET
          ENDIF

          DET=DET*A(K,K)
          IF(A(K,K).EQ.0.0) RETURN
!	To check for singular or nearly singular matrices replace this
!	statement by, where REPS is approximately \hcross*Max(A(I,J))
!         IF(ABS(A(K,K)).LT.REPS) RETURN
          DO 2500 L=K+1,N
            A(L,K)=A(L,K)/A(K,K)
            DO 2500 L1=K+1,N
2500      A(L,L1)=A(L,L1)-A(L,K)*A(K,L1)
2600    CONTINUE
        DET=DET*A(N,N)
        INC(N)=N
!	If pivot is zero then return, IER has been set to 121
        IF(A(N,N).EQ.0.0) RETURN
!	To check for singular or nearly singular matrices replace this
!	statement by, where REPS is approximately \hcross*Max(A(I,J))
!         IF(ABS(A(N,N)).LT.REPS) RETURN

        IER=0
        IF(IFLG.EQ.1) THEN
          IFLG=2
          RETURN
        ENDIF
        IFLG=2
      ENDIF

      IER=0
!	Solution for the NUM different right-hand sides
      DO 5000 J=1,NUM
        DO 3000 K=1,N-1
!	Forward substitution
          IF(K.NE.INC(K)) THEN
            T1=X(K,J)
            X(K,J)=X(INC(K),J)
            X(INC(K),J)=T1
          ENDIF
          DO 3000 L=K+1,N
3000    X(L,J)=X(L,J)-A(L,K)*X(K,J)

!	back-substitution
        X(N,J)=X(N,J)/A(N,N)
        DO 3300 K=N-1,1,-1
          DO 3200 L=N,K+1,-1
3200      X(K,J)=X(K,J)-X(L,J)*A(K,L)
3300    X(K,J)=X(K,J)/A(K,K)
5000  CONTINUE
      END

!     ----------------------------------------------

!	To calculate the Singular Value Decomposition of a matrix A=U D V-transpose
!
!	N : (input) Number of variables
!	M : (input) Number of equations
!	A : (input/output) Matrix of coefficients of size LA*N
!		After execution it will contain the matrix U
!	V : (output) The matrix V of size LV*N
!	SIGMA : (output) Array of length N, containing the singular values
!	LA : (input) Actual value of first dimension of A in the calling program
!	LV : (input) Actual value of first dimension of V in the calling program
!	E : Scratch array of length N
!	IER : (output) Error parameter, IER=0 if execution is successful
!		IER=12 QR iteration failed to converge to required accuracy
!		IER=105 implies N.LE.0, N.GT.LV, M.LE.0, M.GT.LA, N.GT.M
!
!	Required routines : None

      SUBROUTINE SVD(N,M,A,V,SIGMA,LA,LV,E,IER)
!      IMPLICIT REAL*8(A-H,O-Z)
!      PARAMETER(ITMAX=30,REPS=1.D-16)
!	For REAL*4 use REPS=6.E-8
      PARAMETER(ITMAX=30,REPS=6.E-8)
      DIMENSION A(LA,N),V(LV,N),SIGMA(N),E(N)

      IF(N.GT.M.OR.N.LE.0.OR.M.LE.0.OR.M.GT.LA.OR.N.GT.LV) THEN
        IER=105
        RETURN
      ENDIF

      IER=0
!	Reduction to Bidiagonal form using Householder transformations
      G=0
      RMAX=0

      DO 3000 I=1,N
!	Off-diagonal elements of bidiagonal form
        E(I)=G
        S=0
        DO 1200 J=I,M
1200    S=S+A(J,I)**2

        IF(S.LE.0.0) THEN
!	transformation not required
          G=0
        ELSE
          F=A(I,I)
          G=SQRT(S)
          IF(F.GE.0.0) G=-G
          H=F*G-S
          A(I,I)=F-G

          DO 1800 J=I+1,N
            S=0
            DO 1400 K=I,M
1400        S=S+A(K,I)*A(K,J)
            F=S/H
            DO 1600 K=I,M
1600        A(K,J)=A(K,J)+F*A(K,I)
1800      CONTINUE
        ENDIF

!	Diagonal elements of bidiagonal form
        SIGMA(I)=G
        S=0
        DO 2000 J=I+1,N
2000    S=S+A(I,J)**2

        IF(S.LE.0.0) THEN
!	Transformation not required
          G=0
        ELSE
          F=A(I,I+1)
          G=SQRT(S)
          IF(F.GE.0.0) G=-G
          H=F*G-S
          A(I,I+1)=F-G
          DO 2200 J=I+1,N
!	Temporary storage of intermediate results
2200      E(J)=A(I,J)/H

          DO 2800 J=I+1,M
            S=0
            DO 2400 K=I+1,N
2400        S=S+A(J,K)*A(I,K)
            DO 2600 K=I+1,N
2600        A(J,K)=A(J,K)+S*E(K)
2800      CONTINUE
        ENDIF
        R1=ABS(SIGMA(I))+ABS(E(I))
        IF(R1.GT.RMAX) RMAX=R1
3000  CONTINUE

!	Accumulation of right hand transformation in array V
      DO 4000 I=N,1,-1
        IF(G.NE.0.0) THEN
          H=A(I,I+1)*G
          DO 3200 J=I+1,N
3200      V(J,I)=A(I,J)/H

          DO 3800 J=I+1,N
            S=0
            DO 3400 K=I+1,N
3400        S=S+A(I,K)*V(K,J)
            DO 3600 K=I+1,N
3600        V(K,J)=V(K,J)+S*V(K,I)
3800      CONTINUE
        ENDIF

        DO 3900 J=I+1,N
          V(I,J)=0.0
          V(J,I)=0.0
3900    CONTINUE
        V(I,I)=1
        G=E(I)
4000  CONTINUE

!	Accumulation of left hand transformation overwritten on matrix A
      DO 5000 I=N,1,-1
        G=SIGMA(I)
        DO 4200 J=I+1,N
4200    A(I,J)=0
        IF(G.NE.0.0) THEN
          H=A(I,I)*G

          DO 4700 J=I+1,N
            S=0
            DO 4400 K=I+1,M
4400        S=S+A(K,I)*A(K,J)
            F=S/H
            DO 4600 K=I,M
4600        A(K,J)=A(K,J)+F*A(K,I)
4700      CONTINUE

          DO 4800 J=I,M
4800      A(J,I)=A(J,I)/G
        ELSE
          DO 4900 J=I,M
4900      A(J,I)=0.0
        ENDIF
        A(I,I)=A(I,I)+1
5000  CONTINUE

!	Diagonalisation of the bidiagonal form
      AEPS=REPS*RMAX
!	Loop over the singular values
      DO 8000 K=N,1,-1
!	The QR transformation
        DO 7500 ITR=1,ITMAX

!	Test for splitting
          DO 5200 L=K,1,-1
            IF(ABS(E(L)).LT.AEPS) GO TO 6000
            IF(ABS(SIGMA(L-1)).LT.AEPS) GO TO 5400
5200      CONTINUE

!	cancellation of E(L) if L>1
5400      C=0.0
          S=1.0
          DO 5800 I=L,K
            F=S*E(I)
            E(I)=C*E(I)
            IF(ABS(F).LT.AEPS) GO TO 6000
            G=SIGMA(I)
            SIGMA(I)=SQRT(F*F+G*G)
            C=G/SIGMA(I)
            S=-F/SIGMA(I)

            DO 5600 J=1,M
              R1=A(J,L-1)
              R2=A(J,I)
              A(J,L-1)=R1*C+R2*S
              A(J,I)=C*R2-S*R1
5600        CONTINUE
5800      CONTINUE

6000      Z=SIGMA(K)
          IF(L.EQ.K) THEN
!	QR iteration has converged
            IF(Z.LT.0.0) THEN
              SIGMA(K)=-Z
              DO 6200 J=1,N
6200          V(J,K)=-V(J,K)
            ENDIF
            GO TO 8000
          ENDIF

          IF(ITR.EQ.ITMAX) THEN
            IER=12
            GO TO 7500
          ENDIF

!	calculating shift from bottom 2x2 minor
          X=SIGMA(L)
          Y=SIGMA(K-1)
          G=E(K-1)
          H=E(K)
          F=((Y-Z)*(Y+Z)+(G-H)*(G+H))/(2.*H*Y)
          G=SQRT(1.+F*F)
          IF(F.LT.0.0) G=-G
          F=((X-Z)*(X+Z)+H*(Y/(F+G)-H))/X

!	next QR transformation
          C=1.0
          S=1.0
!	Given's rotation
          DO 7000 I=L+1,K
            G=E(I)
            Y=SIGMA(I)
            H=S*G
            G=C*G
            E(I-1)=SQRT(F*F+H*H)
            C=F/E(I-1)
            S=H/E(I-1)
            F=C*X+S*G
            G=C*G-S*X
            H=S*Y
            Y=C*Y

            DO 6400 J=1,N
              X=V(J,I-1)
              Z=V(J,I)
              V(J,I-1)=C*X+S*Z
              V(J,I)=C*Z-S*X
6400        CONTINUE

            SIGMA(I-1)=SQRT(F*F+H*H)
            IF(SIGMA(I-1).NE.0.0) THEN
              C=F/SIGMA(I-1)
              S=H/SIGMA(I-1)
            ENDIF
            F=C*G+S*Y
            X=C*Y-S*G
            DO 6600 J=1,M
              Y=A(J,I-1)
              Z=A(J,I)
              A(J,I-1)=C*Y+S*Z
              A(J,I)=C*Z-S*Y
6600        CONTINUE
7000      CONTINUE

          E(L)=0
          E(K)=F
          SIGMA(K)=X
7500    CONTINUE
8000  CONTINUE
      END

!     ----------------------------------------------

!	To evaluate the solution of a system of linear equations using SVD
!
!	N : (input) Number of variables
!	M : (input) Number of equations
!	U : (input) array of size LU*N containing the left-hand transformation
!	V : (input) array of size LV*N containing the right-hand transformation
!	SIGMA : (input) array of size N containing the singular values
!	LU : (input) First dimension of array U in the calling program
!	LV : (input) First dimension of array V in the calling program
!	B : (input/output) Array of length M containing the RHS
!		after execution it will contain the solution
!	WK : Scratch array of length N
!	REPS : (input) Relative accuracy.
!               All singular values < REPS*(Max of singular values)
!		will be reduced to zero
!
!	Required routines : None

      SUBROUTINE SVDEVL(N,M,U,V,SIGMA,LU,LV,B,WK,REPS)
!      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION U(LU,N),V(LV,N),SIGMA(N),B(*),WK(N)

!	Finding the largest singular value
      SMAX=0.0
      DO 2000 I=1,N
        IF(SIGMA(I).GT.SMAX) SMAX=SIGMA(I)
2000  CONTINUE

      AEPS=SMAX*REPS
      DO 3000 I=1,N
        S=0.0
!	Only SIGMA(I) > AEPS contribute to the solution
        IF(SIGMA(I).GT.AEPS) THEN
          DO 2400 J=1,M
2400      S=S+U(J,I)*B(J)
          S=S/SIGMA(I)
        ENDIF
        WK(I)=S
3000  CONTINUE

      DO 4000 I=1,N
        S=0.0
        DO 3400 J=1,N
3400    S=S+V(I,J)*WK(J)
        B(I)=S
4000  CONTINUE
      END
