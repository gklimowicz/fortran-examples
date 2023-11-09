!     PROGRAM TO SOLVE FREDHOLM EQUATIONS USING QUADRATURE METHOD

      PROGRAM FREDHM
!      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(65),WT(65),WK(65,135),INT(65),F(65),FC(65)
      EXTERNAL FG,FKER

!     EXAMPLE 13.1 : FREDHOLM EQUATION OF THE SECOND KIND

51    FORMAT('   IER =',I4,5X,'NO. OF PTS =',I5,5X,'IQ =',I4/
     1       16X,1HX,11X,'SOLUTION',5X,'COR. SOL.(IF IQ=1)')
52    FORMAT(I6,2X,1P3E16.6)
53    FORMAT(/,'  VALUES CALCULATED USING QUADRATURE fORMULA')
54    FORMAT('   IER =',I4,5X,'NO. OF PTS =',I5,5X,'IQ =',I4/
     1       16X,1HX,11X,'SOLUTION')

      REPS=1.E-6
      IT=2
      A=0
      B=1

100   PRINT *,'N=NO. OF PTS,   IQ=1/2/4/8/16/32   (QUITS WHEN N.LE.0)'
      PRINT *,'IQ=1 FOR TRAPEZOIDAL RULE,  IQ=2 FOR SIMPSON''S RULE'
      PRINT *,'IQ=4/8/16/32   FOR GAUSSIAN FORMULAS'
      READ *,N,IQ
      IF(N.LE.0) STOP

      CALL FRED(N,A,B,WT,X,F,FC,FG,FKER,E0,WK,INT,IQ,IT,REPS,IER)
!	FC is calculated only for IQ=1
      IF(IQ.EQ.1) THEN
        WRITE(6,51) IER,N,IQ
        DO 1000 I=1,N
          WRITE(6,52) I,X(I),F(I),FC(I)
1000    CONTINUE
      ELSE
        WRITE(6,54) IER,N,IQ
        DO 1200 I=1,N
          WRITE(6,52) I,X(I),F(I)
1200    CONTINUE
      ENDIF

!     FOR FREDHOLM EQUATION OF SECOND KIND USE QUADRATURE FORMULA TO
!     CALCULATE THE SOLUTION AT INTERMEDIATE POINTS

      WRITE(6,53)
      DO 2000 I=1,6
        XI=(I-1)*0.2D0
        FI=-FG(XI)
        DO 1500 J=1,N
1500    FI=FI+WT(J)*FKER(XI,X(J))*F(J)
        WRITE(6,52) I,XI,FI
2000  CONTINUE

      GO TO 100
      END

!     -------------------------------------------------

!	To solve linear Fredholm equation using quadrature method
!
!	M : (input) Number of abscissas to be used in quadrature formula
!	A : (input) Lower limit of the integral
!	B : (input) Upper limit of the integral
!	WT : (input/output) Real array of length M containing the
!		weights used in quadrature formula.
!		If IQ is negative the weights must be supplied by the user
!		otherwise they are calculated by the routine
!	X : (input/output) Real array of length M containing the
!		abscissas used in quadrature formula.
!		If IQ is negative the abscissas must be supplied by the user
!		otherwise they are calculated by the routine
!	F : (output) Real array of length M containing the calculated solution
!		F(I) is the value of solution at X(I).
!	FC : (output) Real array of length M containing the calculated solution
!		after applying the deferred correction. This will be
!		relevant only if IQ=1 and IT=1,2. In other cases a dummy
!		array of any length may be supplied.
!	FG : (input) Name of the function routine used to calculate the right
!		hand side g(X). For IT=3 g(x) is not required, but in that
!		case this function is used to calculate an initial guess
!		for the eigenfunctions. In most cases the inverse iteration
!		converges from essentially arbitrary initial guess and it
!		is enough to set FG(X) to any nonzero value.
!	FKER : (input) Name of the function routine used to calculate the
!		kernel K(x,t)
!	EI : (input/output) Initial guess for the eigenvalue. After execution
!		it will contain the computed eigenvalue. Used only for IT=3
!	WK : Real array used as scratch space. Its length should be M*M for
!		IT=1,2, while for IT=3 it should be 2M*M+M.
!	IWK : Integer array of length M used as a scratch space
!	IQ : (input) Integer variable to specify the quadrature formula to be used.
!		If IQ=1, then trapezoidal rule is used and deferred correction
!			is calculated using Gregory's formula
!		If IQ=2, the Simpson's 1/3 rule is used
!		If IQ=4,8,16,32 a composite rule using IQ point
!			Gauss-Legendre formula is used.
!		In all these cases the weights and abscissas are calculated
!		If IQ is negative then it is assumed that weights and abscissas
!		are supplied in arrays WT and X.
!		Other values of IQ will cause an error return.
!	IT : (input) Integer variable to specify the type of integral equation
!		If IT=1 Fredholm equation of the first kind is solved
!		If IT=2 Fredholm equation of the second kind is solved
!		If IT=3 Fredholm equation of the third kind (eigenvalue
!			problem) is solved
!	REPS : (input) Required relative accuracy in calculating eigenvalue
!		and eigenvectors. It is not used for IT=1,2. It only specifies
!		the accuracy to which inverse iteration converges. It does not
!		control the truncation error.
!	IER : (output) The error parameter, IER=0 implies successful execution
!		IER=-11 implies that the calculations are actually performed
!			using a smaller number of abscissas than M
!		IER=706 implies M<3, IT>3, IT.LE.0 or M is not sufficient
!			to apply the required quadrature formula. No calculations
!			are done.
!		IER=707 implies that IQ is not acceptable and no calculations
!			are done
!		Other values of IER may be set by GAUELM and INVIT
!
!	FUNCTION FG(X) and FUNCTION FKER(X,T) must be supplied by the user.
!		FG is the right hand side function g(x) and FKER(X,T) is 
!		the kernel. The integral equation is specified in the form
!		given by Eq.(12.1)
!
!	Required routines : GAUELM, INVIT, FG, FKER
!	
      SUBROUTINE FRED(M,A,B,WT,X,F,FC,FG,FKER,EI,WK,IWK,IQ,IT,REPS,IER)
!      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION WT(M),F(M),FC(*),X(M),WK(M,*),IWK(M)
      DIMENSION W(31),XA(31)

!	Weights and abscissas for Gauss-Legendre quadrature.
!	For N-point formula W(K)=W(N-K+1) and XA(K)=-XA(N-K+1)
!		For K=1,2,...,N/2. Hence only half points are tabulated.
!	For 2-point W(1); 4-point W(2), W(3); 8-point W(4),...,W(7);
!	16-point W(8),...,W(15); 32-point W(16),...,W(31) are the
!	weights corresponding to abscissas XA(I).

      DATA W/1.0D0,
     1       0.34785484513745385737D0, 0.65214515486254614263D0,
     2       0.10122853629037625915D0, 0.22238103445337447054D0,
     3       0.31370664587788728734D0, 0.36268378337836198297D0,
     4       0.02715245941175409485D0, 0.06225352393864789286D0,
     5       0.09515851168249278481D0, 0.12462897125553387205D0,
     6       0.14959598881657673208D0, 0.16915651939500253819D0,
     7       0.18260341504492358887D0, 0.18945061045506849629D0,
     8       0.00701861000947009660D0, 0.01627439473090567061D0,
     9       0.02539206530926205945D0, 0.03427386291302143310D0,
     1       0.04283589802222668066D0, 0.05099805926237617620D0,
     2       0.05868409347853554714D0, 0.06582222277636184684D0,
     3       0.07234579410884850623D0, 0.07819389578707030647D0,
     4       0.08331192422694675522D0, 0.08765209300440381114D0,
     4       0.09117387869576388471D0, 0.09384439908080456564D0,
     5       0.09563872007927485942D0, 0.09654008851472780057D0/

      DATA XA/0.57735026918962576451D0,
     1        0.86113631159405257522D0, 0.33998104358485626480D0,
     2        0.96028985649753623168D0, 0.79666647741362673959D0,
     3        0.52553240991632898582D0, 0.18343464249564980494D0,
     4        0.98940093499164993260D0, 0.94457502307323257608D0,
     5        0.86563120238783174388D0, 0.75540440835500303390D0,
     6        0.61787624440264374845D0, 0.45801677765722738634D0,
     7        0.28160355077925891323D0, 0.09501250983763744019D0,
     8        0.99726386184948156354D0, 0.98561151154526833540D0,
     9        0.96476225558750643077D0, 0.93490607593773968917D0,
     1        0.89632115576605212397D0, 0.84936761373256997013D0,
     2        0.79448379596794240696D0, 0.73218211874028968039D0,
     3        0.66304426693021520098D0, 0.58771575724076232904D0,
     4        0.50689990893222939002D0, 0.42135127613063534536D0,
     5        0.33186860228212764978D0, 0.23928736225213707454D0,
     6        0.14447196158279649349D0, 0.04830766568773831623D0/

      IER=0
      IF(M.LT.3.OR.IT.GT.3.OR.IT.LE.0) THEN
        IER=706
        RETURN
      ENDIF

!	M should not be changed since it is used in dimension statement
      N=M
      IF(IQ.EQ.1) THEN
!	Use the trapezoidal rule
        H=(B-A)/(N-1)
        WT(1)=0.5*H
        WT(N)=WT(1)
        X(1)=A
        X(N)=B
        DO 2000 I=2,N-1
          X(I)=A+(I-1)*H
2000    WT(I)=H

      ELSE IF(IQ.EQ.2) THEN
!	Use the Simpson's 1/3 rule, if N is even, then reduce it by 1
        N=2*((N-1)/2)+1
        H=(B-A)/(N-1)
        WT(1)=H/3.
        X(1)=A
        DO 2100 I=2,N-1,2
          X(I)=A+(I-1)*H
          X(I+1)=A+I*H
          WT(I)=4.*H/3.
2100    WT(I+1)=2.*H/3.
        WT(N)=WT(1)
        X(N)=B

      ELSE IF(IQ.GE.0) THEN
!	Try Gauss-Legendre formulas
        NO=-1
        IF(IQ.EQ.4) NO=1
        IF(IQ.EQ.8) NO=3
        IF(IQ.EQ.16) NO=7
        IF(IQ.EQ.32) NO=15
        IF(NO.LT.0) THEN
!	If IQ is not acceptable then quit
          IER=707
          RETURN
        ENDIF
        N=(N/IQ)*IQ
        IF(N.LT.IQ) THEN
!	If the number of points is not sufficient, then quit
          IER=706
          RETURN
        ENDIF

!	Setup the weights and abscissas for Gauss-Legendre formula
        H=IQ*(B-A)/N
        DO 2300 I=1,N,IQ
          A1=A+(I-1)*H/IQ+H/2.
          DO 2300 I1=1,IQ/2
            WT(I+I1-1)=W(NO+I1)*H/2.
            WT(I+IQ-I1)=WT(I+I1-1)
            X(I+I1-1)=A1-XA(NO+I1)*H/2.
            X(I+IQ-I1)=A1+XA(NO+I1)*H/2.
2300    CONTINUE
      ENDIF
      IF(M.NE.N) IER=-11

!	Setting up the equation matrix and the right hand side
      DO 3000 I=1,N
        F(I)=FG(X(I))
        DO 2800 J=1,N
          WK(J,I)=WT(I)*FKER(X(J),X(I))
2800    CONTINUE
        IF(IT.EQ.2) WK(I,I)=WK(I,I)-1.
3000  CONTINUE

      NUM=1
      LJ=M
      IF(IT.LE.2) THEN
!	Solve the system of linear equations
        IFLG=0
        CALL GAUELM(N,NUM,WK,F,DET,IWK,LJ,IER,IFLG)
        IF(IER.GT.0) RETURN
      ELSE
!	Solve the eigenvalue problem
        IFLG=0
        P=EI
        NIT=0
        CALL INVIT(WK,N,LJ,P,F,IFLG,EI,RC,REPS,WK(1,N+1),IWK,NIT,IER)
      ENDIF

      IF(IQ.EQ.1.AND.IT.NE.3) THEN
!	Apply the deferred correction
        DO 3200 I=1,N
          U2=FKER(X(I),X(2))*F(2)
          U2M=FKER(X(I),X(N-1))*F(N-1)
          T1=FKER(X(I),X(1))*F(1)-U2+FKER(X(I),X(N))*F(N)-U2M
          T2=T1+FKER(X(I),X(3))*F(3)-U2+FKER(X(I),X(N-2))*F(N-2)-U2M
          FC(I)=H*T1/12.+H*T2/24.
3200    CONTINUE
        CALL GAUELM(N,NUM,WK,FC,DET,IWK,LJ,IER,IFLG)

        DO 4000 I=1,N
4000    FC(I)=F(I)+FC(I)
      ENDIF
      END

!     --------------------------------------------------

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

!     -----------------------------------------------------

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
!      IMPLICIT REAL*8(A-H,O-Z)
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

!     ------------------------------------------

      FUNCTION FG(X)
!      IMPLICIT REAL*8(A-H,O-Z)

!     RHS FUNCTION G(X)

      FG=-X**3
      END

!     ---------------------------------------

      FUNCTION FKER(X,T)
!      IMPLICIT REAL*8(A-H,O-Z)

!     THE KERNEL

      IF(T.GT.X) THEN
        FKER=-X*(1.-T)
      ELSE
        FKER=-T*(1.-X)
      ENDIF
      END
