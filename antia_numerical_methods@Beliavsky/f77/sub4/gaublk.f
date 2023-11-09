!	To solve a system of linear equations arising from finite difference
!	approximation of ordinary differential equations
!
!	N : (input) Number of mesh points used.
!	M : (input) Number of first order differential equations in the system
!	ML : (input) Number of boundary conditions at the first boundary t=T(1)
!	A : (input/output) Real array of length (M+ML)*2M*N containing
!		the matrix of equations. After execution it will contain
!		the triangular decomposition
!	IFLG : (input/output) Integer variable used as a flag to decide
!		the nature of computation.
!		If IFLG=0 the triangular decomposition, determinant and
!			solution of equations is computed. IFLG is set to 2
!		If IFLG=1 only the triangular decomposition and determinant
!			are calculated and IFLG is set to 2
!		If IFLG=2 then it is assumed that triangular decomposition
!			is already done and available in A and only the solution
!			of system of equations is solved
!	DET : (output) Scaled value of determinant of the matrix
!	IDET : (output) exponent of determinant, the value of determinant
!		is DET*2**IDET
!	INC : (input/output) Integer array of length M*N containing the
!		information about interchanges used by Gaussian elimination.
!		It is calculated if IFLG=0 or 1, while for IFLG=2 it must
!		be supplied from previous calculation.
!	X : (input/output) Real array of length M*N containing the right hand
!		side of equation at input. After execution it will be overwritten
!		by the solution if IFLG=0 or 2. For IFLG=1 it is not used.
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=735 implies that the matrix is singular
!
!	Required routines : None

      SUBROUTINE GAUBLK(N,M,ML,A,IFLG,DET,IDET,INC,X,IER)
!      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(M+ML,2*M,N),INC(M,N),X(M,N)

      IF(IFLG.LE.1) THEN
        IDET=0
        DET=1.
!	The number of rows in each block of the matrix
        MR=M+ML
!	The number of columns in each block of the matrix
        MC=2*M
        IER=735

        DO 3000 J=1,N
          IF(J.EQ.N) THEN
            MR=M
            MC=M
          ENDIF

          DO 2000 K=1,MIN(M,MR-1)
            RMAX=ABS(A(K,K,J))
            KMAX=K
!	Find the pivot
            DO 1200 KI=K+1,MR
              R1=ABS(A(KI,K,J))
              IF(R1.GT.RMAX) THEN
                RMAX=R1
                KMAX=KI
              ENDIF
1200        CONTINUE
            INC(K,J)=KMAX

            IF(KMAX.NE.K) THEN
!	exchange rows K and KMAX
              DET=-DET
              DO 1300 KI=K,MC
                AT=A(K,KI,J)
                A(K,KI,J)=A(KMAX,KI,J)
1300          A(KMAX,KI,J)=AT
            ENDIF

            DET=DET*A(K,K,J)
!	If the pivot is zero, then quit
            IF(A(K,K,J).EQ.0.0) RETURN

!	Gaussian elimination
            DO 1500 KI=K+1,MR
              A(KI,K,J)=A(KI,K,J)/A(K,K,J)
              DO 1500 KJ=K+1,MC
                A(KI,KJ,J)=A(KI,KJ,J)-A(KI,K,J)*A(K,KJ,J)
1500        CONTINUE
2000      CONTINUE

          IF(DET.NE.0.0) THEN
!	Scale the determinant if necessary
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

!	Copy the overlapping elements into the next block
          IF(J.LT.N) THEN
            DO 2600 K=1,ML
              DO 2600 KI=1,M
2600        A(K,KI,J+1)=A(K+M,KI+M,J)
          ENDIF
3000    CONTINUE
        INC(M,N)=M
        DET=DET*A(M,M,N)
        IF(A(M,M,N).EQ.0.0) RETURN
        IER=0

        IF(IFLG.EQ.1) THEN
          IFLG=2
          RETURN
        ENDIF
        IFLG=2
      ENDIF

!	Solve the system of linear equations
      IER=0
      MR=M+ML
      DO 3100 J=1,N
        IF(J.EQ.N) MR=M
        DO 3100 K=1,MIN(M,MR-1)
          KK=INC(K,J)
          IF(K.NE.KK) THEN
!	exchange the corresponding elements of RHS
            XT=X(K,J)
            X(K,J)=X(KK,J)
            X(KK,J)=XT
          ENDIF

!	Gaussian elimination
          DO 2800 L=K+1,MR
2800      X(L,J)=X(L,J)-A(L,K,J)*X(K,J)
3100  CONTINUE

!	back-substitution
      MC=M
      DO 3500 J=N,1,-1
        DO 3300 K=M,1,-1
          D1=X(K,J)
            DO 3200 L=MC,K+1,-1
3200        D1=D1-X(L,J)*A(K,L,J)
3300    X(K,J)=D1/A(K,K,J)
        MC=2*M
3500  CONTINUE
      END
