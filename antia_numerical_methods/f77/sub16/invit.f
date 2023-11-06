!	Real eigenvalue and eigenvector of a real matrix using inverse iteration
!
!	A : (input) Real array of length IA*M containing the matrix elements
!	M : (input) Order of the matrix
!	IA : (input) The first dimension of A as declared in the calling program
!		IA.GE.M
!	P : (input/output) Initial value of the shift. This will be modified
!		by the program if IFLG>0
!	U : (input/output) Real array of length M, which should specify the
!		initial approximation to eigenvector. After execution it
!		will contain the calculated eigenvector.
!	IFLG : (input) Integer variable to specify the type of iteration required
!		If IFLG=0 the shift P is kept fixed
!		If IFLG=1 the shift P is varied using Rayleigh quotient
!		If IFLG=2 the shift P is varied using max(V_s+1)
!	EI : (output) Estimated eigenvalue using simple inverse iteration
!	ERC : (output) Estimated eigenvalue using Rayleigh quotient
!	REPS : (input) Required absolute accuracy. Iteration is terminated
!		when all components of eigenvector and the eigenvalue have
!		converged to REPS.
!	WK : Real array of length M*(M+1) used as scratch space
!	IWK : Integer array of length M used as scratch space
!	NIT : (input/output) Number of iterations required. If it is
!		zero or negative NIT is set to NIT0=100
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=106 implies that M.LE.1 or M>IA, in which case no
!			calculations are done
!		IER=141 implies that vector is zero at some stage and
!			calculations are aborted
!		IER=142 implies that inverse iteration has failed to converge
!
!	Required routines : GAUELM
!
      SUBROUTINE INVIT(A,M,IA,P,U,IFLG,EI,ERC,REPS,WK,IWK,NIT,IER)
      IMPLICIT REAL*16(A-H,O-Z)
      DIMENSION A(IA,M),U(M),WK(M,M+1),IWK(M)
      PARAMETER(NIT0=100)

      IF(M.LE.1.OR.M.GT.IA) THEN
        IER=106
        RETURN
      ENDIF

!	Copy the matrix to WK and apply the shift
      DO 1100 I=1,M
        WK(I,M+1)=U(I)
        DO 1000 J=1,M
1000    WK(J,I)=A(J,I)
1100  WK(I,I)=A(I,I)-P

      NUM=1
      LJ=M
      IFL=1
      IF(NIT.LE.0) NIT=NIT0
!	Perform Gaussian elimination on A-pI
      CALL GAUELM(M,NUM,WK,U,DET,IWK,LJ,IER,IFL)
      IF(IER.GT.0) RETURN

      EPI=0.0
!	Loop for inverse iteration
      DO 5000 J=1,NIT
        CALL GAUELM(M,NUM,WK,U,DET,IWK,LJ,IER,IFL)
        IF(IER.GT.0) RETURN

!	Normalising the vector U
        R1=0.0
        KM=0
        DO 3500 K=1,M
          IF(R1.LT.ABS(U(K))) THEN
            R1=ABS(U(K))
            KM=K
          ENDIF
3500    CONTINUE
        UKM=U(KM)
        IF(UKM.EQ.0.0) THEN
!	If the vector is zero, then quit
          IER=141
          RETURN
        ENDIF

!	The eigenvalue
        EI=WK(KM,M+1)/UKM+P
        S1=0.0
        S2=0.0
!	Calculating the Rayleigh quotient
        DO 4000 K=1,M
          S1=S1+U(K)*WK(K,M+1)
!	For complex eigenvalues use the following statement instead of the
!	preceding one
!         S1=S1+CONJG(U(K))*WK(K,M+1)
          S2=S2+ABS(U(K))**2
          U(K)=U(K)/UKM
4000    CONTINUE
        ERC=P+S1/S2

!	Convergence check
        R1=ABS(EI-EPI)
        DO 4500 I=1,M
          R1=MAX(R1,ABS(WK(I,M+1)-U(I)))
          WK(I,M+1)=U(I)
4500    CONTINUE
        IF(ABS(R1).LT.REPS) RETURN
        EPI=EI

        IF(IFLG.GE.1) THEN
!	Update the shift
          P=ERC
          IF(IFLG.EQ.2) P=EI
!	Setting up the new matrix A-pI
          DO 4700 I=1,M
            DO 4600 K=1,M
4600        WK(K,I)=A(K,I)
4700      WK(I,I)=A(I,I)-P
          IFL=1
          CALL GAUELM(M,NUM,WK,U,DET,IWK,LJ,IER,IFL)
          IF(IER.GT.0) RETURN
        ENDIF

5000  CONTINUE
!	Iteration fails to converge
      IER=142
      END
