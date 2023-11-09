!     PROGRAM TO SOLVE A SYSTEM OF NONLINEAR EQUATIONS USING BROYDEN'S METHOD
!     IT USES DAVIDENKO'S METHOD TO IMPROVE CHANCES OF CONVERGENCE

      PROGRAM NONLIN
      IMPLICIT REAL*16(A-H,O-Z)
      EXTERNAL FCN
      COMMON/FUN/THETA,X0(20),F0(20)
      DIMENSION X(20),F(20),WK(800),IWK(20)
      DATA (X(I),I=1,8)/-0.25Q0,0.1Q0,-0.25Q0,0.35Q0,-0.25Q0,0.6Q0,
     1   -0.25Q0,0.85Q0/

!     EXAMPLE 7.15


51    FORMAT(/'   F =',1P5D14.6/(6X,5D14.6))
52    FORMAT('   IER =',I4,'    X =',1P4D14.6/(6X,5D14.6))

      NP=8
      REPS=1.Q-25
      AEPS=1.Q-35

      CALL DAVIDN_B(FCN,NP,X,F,REPS,AEPS,IER,WK,IWK)
      WRITE(6,52) IER,(X(I),I=1,NP)
      WRITE(6,51) (F(I),I=1,NP)
      END

!     -------------------------------------------------------

!	Real root of a system of nonlinear equations using Broyden's method
!
!	FCN : (input) Name of the subroutine to calculate the vector function
!	NP : (input) Number of variables and equations to be solved
!	X : (input/output) Real array of length NP containing the starting
!		values. After execution it should contain the computed roots
!	F : (output) Real array of length NP containing the function values
!		at the point specified by array X
!	REPS : (input) Required relative accuracy
!	AEPS : (input) Required absolute accuracy
!		The estimated error should be less than MAX(AEPS, REPS*ABS(X))
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=441 implies that GAUELM failed to solve the required
!			system of linear equations
!		IER=442 implies that iteration did not converge to required
!			accuracy
!	WK : Real array of length MAX(2*NP*NP,NP*(NP+4)) used as scratch space
!	IWK : integer array of length NP used as scratch space
!
!	SUBROUTINE FCN(NP,X,F) must be supplied by the user.
!	Here X and F are real arrays of length NP containing the value of
!	vector function at point X.
!	The system of equations to be solved are F(I)=0
!
!	Required routines : GAUELM, FCN
!
      SUBROUTINE BROYDN(FCN,NP,X,F,REPS,AEPS,IER,WK,IWK)
      IMPLICIT LOGICAL(Q)
      IMPLICIT REAL*16(A-H,O,P,R-Z)
      DIMENSION X(NP),F(NP),WK(NP,*),IWK(NP)

      LJ=NP
      IER=0
      CALL FCN(NP,X,F)

!	Generating initial approximation to the Jacobian
      DO 2000 IP=1,NP
        X1=X(IP)
        DK=0.01Q0*ABS(X1)
        IF(DK.EQ.0.0) DK=MAX(1.Q-6,AEPS*100.)
        X(IP)=X1+DK
        CALL FCN(NP,X,WK)
        DO 1800 J=1,NP
1800    WK(J,NP+IP)=(WK(J,1)-F(J))/DK
        X(IP)=X1
2000  CONTINUE

      NUM=NP
      IFLG=0
      DO 2500 I=1,NP
        DO 2400 J=1,NP
2400    WK(J,I)=0.0
2500  WK(I,I)=1.0

!	To calculate H = B**(-1)
      CALL GAUELM(NP,NUM,WK(1,NP+1),WK,DET,IWK,LJ,IER1,IFLG)
      IF(IER1.GT.0) THEN
!	If GAUELM fails, then quit
        IER=441
        RETURN
      ENDIF

      DO 2800 I=1,NP
        WK(I,NP+1)=F(I)
        WK(I,NP+2)=X(I)
2800  CONTINUE

!	Broyden's iteration
      DO 5000 IPAS=1,200

!	Convergence check
        QCHK=.TRUE.
        DO 3500 J=1,NP
          S1=0.0
          DO 3200 K=1,NP
3200      S1=S1+WK(J,K)*WK(K,NP+1)
          IF(ABS(S1).GT.MAX(REPS*ABS(WK(J,NP+2)),AEPS)) QCHK=.FALSE.
          X(J)=WK(J,NP+2)-S1
          WK(J,NP+3)=-S1
3500    CONTINUE
        IF(QCHK) RETURN

        CALL FCN(NP,X,F)

!	Updating the inverse matrix
        DO 4000 J=1,NP
          WK(J,NP+4)=F(J)-WK(J,NP+1)
          WK(J,NP+1)=F(J)
          WK(J,NP+2)=X(J)
4000    CONTINUE

        SS=0.0
        DO 4300 J=1,NP
          S1=0.0
          S2=0.0
          DO 4100 K=1,NP
            S1=S1+WK(K,NP+3)*WK(K,J)
            S2=S2+WK(J,K)*WK(K,NP+4)
4100      CONTINUE
          F(J)=S1
          X(J)=S2-WK(J,NP+3)
          SS=SS+S1*WK(J,NP+4)
4300    CONTINUE

        DO 4400 J=1,NP
          DO 4400 K=1,NP
            WK(K,J)=WK(K,J)-X(K)*F(J)/SS
4400    CONTINUE
5000  CONTINUE

!	Iteration fails to converge
      IER=442
      END

!     -------------------------------------------------------

!	Real root of a system of nonlinear equations using Davidenko's method
!	For use with subroutine BROYDN
!
!	FCN : (input) Name of the subroutine to calculate the vector function
!	NP : (input) Number of variables and equations to be solved
!	X : (input/output) Real array of length NP containing the starting
!		values. After execution it should contain the computed roots
!	F : (output) Real array of length NP containing the function values
!		at the point specified by array X
!	REPS : (input) Required relative accuracy
!	AEPS : (input) Required absolute accuracy
!		The estimated error should be less than MAX(AEPS, REPS*ABS(X(I)))
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=407 implies NP > NMAX or NP < 1, in which case no
!			calculations are done
!		IER=440 implies that THETA is nonzero even after 100 steps
!		Other values may be set by subroutines NEWTON or BROYDN
!	WK : Scratch real array of length NP*(NP+1) (for use with NEWTON) or
!		length MAX(2*NP*NP,NP*(NP+4)) (for use with BROYDN) 
!	IWK : integer array of length NP used as scratch space
!	THETA : which is read by the subroutine is the parameter introduced
!		into the equations. Normally THETA=1 to start with and
!		initial guess should satisfy the equations with this parameter
!		THETA=0 corresponds to the equations that need to be solved.
!		The value of THETA should be gradually reduced from 1 to 0
!
!	SUBROUTINE FCN(NP,X,F) or SUBROUTINE FCN(NP,X,F,DF)
!	must be supplied by the user. Here X and F are real arrays of length
!	NP with F containing the value of vector function at point X. DF is
!	real array of length NP*NP containing the Jacobian which is required
!	by subroutine NEWTON. If subroutine BROYDN is used then DF is not
!	required. The subroutine FCN should parameterise the equations in
!	terms of THETA as explained in section 7.16. A sample routine is
!	appended at the end of this file.
!
!	Required routines : BROYDN (or NEWTON), GAUELM, FCN
!
      SUBROUTINE DAVIDN_B(FCN,NP,X,F,REPS,AEPS,IER,WK,IWK)
      IMPLICIT REAL*16(A-H,O-Z)
      PARAMETER(NMAX=200)
      EXTERNAL FCN
      DIMENSION X(NP),F(NP),WK(*),IWK(NP)
!	To pass on parameters to FCN this common block must be included
!	in subroutine FCN
      COMMON/DNFUN/THETA,X0(NMAX),F0(NMAX)

51    FORMAT(/5X,'THETA =',F6.2/)
52    FORMAT('   IER =',I4,'    X =',1P4D14.6/(6X,5D14.6))

      IF(NP.GT.NMAX.OR.NP.LT.1) THEN
        IER=407
        RETURN
      ENDIF

      PRINT *,' TYPE IN THE INITIAL GUESS'
      READ *,(X(I),I=1,NP)

!	To initialise the parameterised function
      THETA=0.0
!	For NEWTON use the following form
!      CALL FCN(NP,X,F0,WK)
      CALL FCN(NP,X,F0)
      DO 500 I=1,NP
500   X0(I)=X(I)
      REPS0=MAX(1.Q-3,REPS)
      AEPS0=MAX(1.Q-4,AEPS)

      DO 1000 I=1,100
        PRINT *,' TYPE THETA'
        READ *,THETA
        PRINT 51, THETA
        IF(THETA.EQ.0.0) THEN
          REPS0=REPS
          AEPS0=AEPS
        ENDIF

!	Use either NEWTON or BROYDN
!        CALL NEWTON(FCN,NP,X,F,REPS0,AEPS0,IER,WK,IWK)
       CALL BROYDN(FCN,NP,X,F,REPS0,AEPS0,IER,WK,IWK)
        PRINT 52,IER,(X(J),J=1,NP)
        IF(IER.GT.0) RETURN
        IF(THETA.EQ.0.0) RETURN
1000  CONTINUE

      IER=440
      END

!	-------------------
!
!	For use with NEWTON include DF
!
!	SUBROUTINE FCN(NP,X,F)
!	IMPLICIT REAL*16(A-H,O-Z)
!	PARAMETER (NMAX=200)
!	DIMENSION X(NP),F(NP)
!	COMMON/DNFUN/THETA,X0(NMAX),F0(NMAX)
!
!	call the required function to calculate the vector function
!	There is no provision to pass the name FUN, it has to be put explicitly
!	CALL FUN(NP,X,F)
!
!	APPLY PARAMETRISATION IN TERMS OF THETA
!
!	DO 2222 I=1,NP
!	F(I)=F(I)-THETA*F0(I)
!2222	CONTINUE
!	END

!  -----------------------------------------------------

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
      IMPLICIT REAL*16(A-H,O-Z)
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

!     --------------------------------------------------

      SUBROUTINE FCN(N,X,F)
      IMPLICIT REAL*16(A-H,O-Z)
      DIMENSION X(*),F(*)
      COMMON/DNFUN/THETA,X0(200),F0(200)

!     SPECIFY THE SYSTEM OF NONLINEAR EQUATIONS  F(I)=0, FOR I=1,...,N

      DO 2000 I=1,N
        FI=1.Q0/I**2
        DO 1500 J=1,N,2
          IF(I.EQ.1) FI=FI+X(J)
          IF(I.GT.1) FI=FI+X(J)*X(J+1)**(I-1)
1500    CONTINUE

!     PARAMETERISATION FOR DAVIDENKO'S METHOD

        F(I)=FI-THETA*F0(I)
2000  CONTINUE
      END
