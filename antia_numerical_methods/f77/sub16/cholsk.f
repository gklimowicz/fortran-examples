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
      IMPLICIT REAL*16(A-H,O-Z)
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
 
