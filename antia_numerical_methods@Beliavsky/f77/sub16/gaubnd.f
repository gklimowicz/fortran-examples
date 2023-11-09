!     Solution of a system of linear equations using Gaussian elimination
!     	for a band matrix
!
!     N : (input) Number of equations to be solved
!     KB : (input) Bandwidth of matrix A(I,J)=0 if ABS(I-J)>KB
!     NUM : (input) Number of different sets (each with N equations) of
!		equations to be solved
!     A : (input/output) The matrix of coefficient of size LJ*(3*KB+1)
!		A(I,J-I+KB+1) is the coefficient of x_j in Ith equation
!		at output it will contain the triangular decomposition
!     X : (input/output) The matrix containing right hand sides (size LJ*NUM)
!		X(I,J) is the Ith element of Jth right hand side
!		at output it will contain the solutions
!     DET, IDET : (output) The determinant of the matrix = DET*2**IDET
!     INC : (output) Integer array of length N containing information about
!		interchanges performed during elimination
!     LJ : (input) First dimension of arrays A and X in calling program
!     IER : (output) Error flag, IER=0 signifies successful execution
!		IER=104 implies (N.LE.0 or N.GT.LJ or KB.GT.N)
!		IER=124 implies some pivot turned out to be zero and hence
!	     		matrix must be nearly singular
!     IFLG : (input) Integer variable used as a flag to specify the type
!		of computation required
!     		If IFLG=-1, both elimination and solution are calculated
!     			without pivoting and IFLG is set to 2
!		If IFLG=0, both elimination and solution are computed
!     			with partial pivoting and IFLG is set to 2
!     		If IFLG=1, only elimination is done with pivoting and IFLG is set to 2
!     		If IFLG.GE.2 only solution is calculated, the triangular
!     			decomposition should have been calculated earlier
!     WK : Real array of length 3*KB+1 used as scratch space
!
!	Required routines : None
 
      SUBROUTINE GAUBND(N,KB,NUM,A,X,DET,IDET,INC,LJ,IER,IFLG,WK)
      IMPLICIT REAL*16(A-H,O-Z)
!     For complex matrices use the following statements instead
!      IMPLICIT REAL*8(R)
!      IMPLICIT COMPLEX*16(A-H,S-Z)
 
      DIMENSION A(LJ,3*KB+1),INC(N),X(LJ,NUM),WK(3*KB+1)
 
      IF(N.LE.0.OR.N.GT.LJ.OR.KB.GT.N) THEN
        IER=104
        RETURN
      ENDIF
 
      KB1=KB+1
      IER=124
      IF(IFLG.LE.1) THEN
!     Perform elimination
        DO 2000 I=1,N
          DO 2000 J=2*KB+2,3*KB+1
            A(I,J)=0.0
2000    CONTINUE
 
        DET=1.0
        IDET=0
        DO 2600 K=1,N-1
!     Find the maximum element in the Kth column
          R1=0.0
          KM=K
          IF(IFLG.GE.0) THEN
            DO 2200 L=K,MIN(N,K+KB)
              IF(ABS(A(L,K-L+KB1)).GT.R1) THEN
                R1=ABS(A(L,K-L+KB1))
                KM=L
              ENDIF
2200        CONTINUE
          ENDIF
 
          INC(K)=KM
          IF(KM.NE.K) THEN
!     Interchange the rows if needed
            DO 2300 L=K,MIN(N,2*KB+K)
              WK(L-K+1)=A(K,L-K+KB1)
2300        CONTINUE
            DO 2400 L=K,MIN(N,2*KB+K)
              A(K,L-K+KB1)=A(KM,L-KM+KB1)
2400        A(KM,L-KM+KB1)=WK(L-K+1)
            DET=-DET
          ENDIF
 
          DET=DET*A(K,KB1)
          IF(A(K,KB1).EQ.0.0) RETURN
!     To check for singular or nearly singular matrices replace this
!     statement by, where REPS is approximately \hcross*Max(A(I,J))
!         IF(ABS(A(K,KB1)).LT.REPS) RETURN
          IF(DET.NE.0.0) THEN
 
!     Scale the value of the determinant DET
2350        IF(ABS(DET).GT.32.) THEN
              DET=DET*0.03125Q0
              IDET=IDET+5
              GO TO 2350
            ENDIF
 
2370        IF(ABS(DET).LT.0.03125Q0) THEN
              DET=DET*32.
              IDET=IDET-5
              GO TO 2370
            ENDIF
          ENDIF
 
          DO 2500 L=K+1,MIN(N,K+KB)
            A(L,K-L+KB1)=A(L,K-L+KB1)/A(K,KB1)
            DO 2500 L1=K+1,MIN(N,2*KB+K)
2500      A(L,L1-L+KB1)=A(L,L1-L+KB1)-A(L,K-L+KB1)*A(K,L1-K+KB1)
2600    CONTINUE
        DET=DET*A(N,KB1)
        INC(N)=N
!     If pivot is zero then return, IER has been set to 124
        IF(A(N,KB1).EQ.0.0) RETURN
!     To check for singular or nearly singular matrices replace this
!     statement by, where REPS is approximately \hcross*Max(A(I,J))
!         IF(ABS(A(N,KB1)).LT.REPS) RETURN
 
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
        DO 3000 K=1,N-1
!     Forward substitution
          IF(K.NE.INC(K)) THEN
            T1=X(K,J)
            X(K,J)=X(INC(K),J)
            X(INC(K),J)=T1
          ENDIF
          DO 3000 L=K+1,MIN(N,K+KB)
3000    X(L,J)=X(L,J)-A(L,K-L+KB1)*X(K,J)
 
!     back-substitution
        X(N,J)=X(N,J)/A(N,KB1)
        DO 3300 K=N-1,1,-1
          DO 3200 L=MIN(N,K+2*KB),K+1,-1
3200      X(K,J)=X(K,J)-X(L,J)*A(K,L-K+KB1)
3300    X(K,J)=X(K,J)/A(K,KB1)
5000  CONTINUE
      END
