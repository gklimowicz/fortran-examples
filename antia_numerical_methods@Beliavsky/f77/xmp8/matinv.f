!     PROGRAM TO FIND INVERSE OF A REAL MATRIX

      PROGRAM MATIV
      IMPLICIT REAL*8(A-H,O-Z)
!      IMPLICIT COMPLEX*16(A-H,S-Z)
      DIMENSION X(20,20),AM(20,20),INC(20)

51    FORMAT('   IER =',I4,5X,'N =',I5/
     1       '    THE INVERSE MATRIX IS')
52    FORMAT('   THE MATRIX IS')
53    FORMAT(5(2X,2F8.3))
54    FORMAT(2(2X,1P2D14.6))

      LJ=20

100   PRINT *,'TYPE N=ORDER OF MATRIX   (QUITS WHEN N.LE.0)'
      READ *,N
      IF(N.LE.0) STOP
      WRITE(6,52)
      PRINT *,'TYPE THE MATRIX ROW-WISE'

!     READ THE MATRIX 

      DO 1000 I=1,N
        PRINT *,I,'TH ROW'
        READ *,(AM(I,J),J=1,N)
        WRITE(6,53) (AM(I,J),J=1,N)
1000  CONTINUE

      CALL MATINV(N,LJ,AM,X,INC,IER)
      WRITE(6,51) IER,N
      DO 2000 I=1,N
2000  WRITE(6,54) (X(I,J),J=1,N)

      GO TO 100
      END

!	---------------------------------------------------------

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
      IMPLICIT REAL*8(A-H,O-Z)
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

!     ---------------------------------------------------------

!	To calculate inverse of a square matrix
!
!	N : (input) Order of matrix
!	IA : (input) The first dimension of arrays A and AI as specified
!		in the calling program
!	A : (input) Real array of length IA*N containing the matrix
!	AI : (output) Real array of length IA*N which will contain the
!		calculated inverse of A
!	IWK : (output) Integer array of length N used as scratch space
!	IER : (output) The error parameter, IER=0 implies successful execution
!		Nonzero values may be set by subroutine GAUELM
!	It is possible to use CROUT instead of GAUELM for calculating
!	the triangular decomposition, but in that case an extra real scratch
!	array of size N will be required.
!
!	Required routines : GAUELM
 
      SUBROUTINE MATINV(N,IA,A,AI,IWK,IER)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(IA,N),AI(IA,N),IWK(N)
 
      DO 1000 I=1,N
        DO 800 J=1,N
          AI(J,I)=0.0
800     CONTINUE
        AI(I,I)=1.0
1000  CONTINUE
 
      NUM=N
      IFLG=0
      CALL GAUELM(N,NUM,A,AI,DET,IWK,IA,IER,IFLG)
 
      END
