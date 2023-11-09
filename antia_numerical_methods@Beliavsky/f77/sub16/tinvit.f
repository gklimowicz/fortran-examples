!	To find a specified eigenvalue and eigenvector of a real symmetric
!	tridiagonal matrix using inverse iteration
!
!	E : (input) Real array of length N containing the off-diagonal elements
!		of the tridiagonal matrix, E(i+1)=A(i,i+1)=A(i+1,i)
!	D : (input) Real array of length N containing the diagonal elements
!		of the tridiagonal matrix, D(i)=A(i,i)
!	N : (input) Order of the matrix
!	EL : (input) lower limit on eigenvalue
!	EU : (input) upper limit on eigenvalue
!	EI : (output) The calculated eigenvalue
!	EV : (input/output) Real array of length N containing the calculated
!		eigenvector. If IFLG.NE.0 then EV should contain the previous
!		eigenvector, so that initial guess to EV is orthogonal to it
!	REPS : (input) Relative accuracy to which eigenvalues or eigenvectors
!		are calculated
!	IFLG : (input) Integer variable used as a flag to denote how initial
!		approximation to eigenvector is obtained. If IFLG=0 then
!		initial guess to eigenvector is generated using random numbers.
!		Otherwise, initial guess is chosen to be orthogonal to EV
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=144 implies that inverse iteration failed to converge
!	U : Real array of length N used as scratch space
!	B : Real array of length 4N used as scratch space
!	NUM : (output) Number of iterations required by the subroutine
!
!	FUNCTION RAN1(SEED) is required to generate random numbers.
!	If a different routine is used, the seed should be changed
!	appropriately.
!
!	Required routines : RAN1
!
      SUBROUTINE TINVIT(E,D,N,EL,EU,EI,EV,REPS,IFLG,IER,U,B,NUM)
      IMPLICIT REAL*16(A-H,O-Z)
      PARAMETER(NIT=100)
      EXTERNAL RAN1
      DIMENSION E(N),D(N),EV(N),B(4,N),U(N)
!	The seed for random number generator
      DATA SEED/12345/

      E(1)=0.0
      IF(IFLG.EQ.0) THEN
!	Choose a random initial guess for eigenvector
        DO 1000 I=1,N
          U(I)=RAN1(SEED)
          EV(I)=U(I)
1000    CONTINUE
      ELSE

!	Choose the initial vector to be orthogonal to EV
        S=0.0
        EVS=0.0
        DO 1200 I=1,N
          U(I)=RAN1(SEED)
          EVS=EVS+EV(I)**2
1200    S=S+U(I)*EV(I)
        S=S/EVS
        DO 1400 I=1,N
          EV(I)=U(I)-EV(I)*S
1400    U(I)=EV(I)
      ENDIF

      EI=(EL+EU)/2.
!	Stretch the bounds to allow for roundoff errors
      EL1=EL-REPS*ABS(EL)
      EU1=EU+REPS*ABS(EU)
      IER=0

!	Loop for inverse iteration
      DO 5000 J=1,NIT
        NUM=J

!	Set up the matrix A-pI
        DO 2000 I=1,N
          B(1,I)=E(I)
          B(2,I)=D(I)-EI
          IF(I.LT.N) B(3,I)=E(I+1)
          B(4,I)=0.0
2000    CONTINUE
        B(3,N)=0.0

!	Gaussian elimination with partial pivoting
        DO 3000 K=1,N-1
          IF(ABS(B(2,K)).LT.ABS(B(1,K+1))) THEN
            DO 2200 I=1,3
              T=B(I,K+1)
              B(I,K+1)=B(I+1,K)
              B(I+1,K)=T
2200        CONTINUE
            T=EV(K+1)
            EV(K+1)=EV(K)
            EV(K)=T
          ENDIF

!	If a pivot is zero make it nonzero
          IF(B(2,K).EQ.0.0) B(2,K)=REPS
          R=-B(1,K+1)/B(2,K)
          B(2,K+1)=B(2,K+1)+R*B(3,K)
          B(3,K+1)=B(3,K+1)+R*B(4,K)
          EV(K+1)=EV(K+1)+R*EV(K)
3000    CONTINUE
        IF(B(2,N).EQ.0.0) B(2,N)=REPS

!	Back-substitution
        EV(N)=EV(N)/B(2,N)
        EV(N-1)=(EV(N-1)-EV(N)*B(3,N-1))/B(2,N-1)
        DO 3400 K=N-2,1,-1
          EV(K)=(EV(K)-EV(K+1)*B(3,K)-EV(K+2)*B(4,K))/B(2,K)
3400    CONTINUE

!	Calculating the Rayleigh quotient (DEI)
        R1=0.0
        KM=0
        S1=0.0
        S2=0.0
        DO 3500 K=1,N
          S1=S1+EV(K)*EV(K)
          S2=S2+EV(K)*U(K)
          IF(R1.LT.ABS(EV(K))) THEN
            R1=ABS(EV(K))
            KM=K
          ENDIF
3500    CONTINUE

        DEI=S2/S1
        UKM=EV(KM)
        DVS=0.0
!	Normalising the eigenvector
        DO 4000 K=1,N
          EV(K)=EV(K)/UKM
          DVS=DVS+ABS(EV(K)-U(K))
4000    CONTINUE

!	The convergence check
        REP=REPS*ABS(EI)
        IF(ABS(DEI).LT.REP.OR.R1*REP.GT.1.0.OR.DVS.LT.REPS) RETURN
        DO 4500 I=1,N
          U(I)=EV(I)
4500    CONTINUE
        IF(EI+DEI.GT.EL1.AND.EI+DEI.LT.EU1) THEN
!	The new shift
          EI=EI+DEI
        ELSE
!	If the new shift is outside the specified bounds, then modify it
          HI=MIN(ABS(EI-EL),ABS(EI-EU))
          EI=EI+SIGN(HI*0.5,DEI)
        ENDIF
5000  CONTINUE

!	Iteration fails to converge
      IER=144
      END
