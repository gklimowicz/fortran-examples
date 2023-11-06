!     PROGRAM TO OBTAIN POLYNOMIAL L1-APPROXIMATION FOR DISCRETE DATA
!     BOTH POLYL1 AND LINL1 ARE USED TO CALCULATE THE SAME APPROXIMATION

      PROGRAM L1POLY
      IMPLICIT REAL*16(A-H,O,P,R-Z)
      DIMENSION B(20),WK(9000),IWK(230),F(200),X(200),g(100,200)

!     EXAMPLE 10.13
!     TO GENERATE DATA SET WITH RANDOM ERROR
      FM(Y)=((231*Y*Y-315)*Y*Y+105)*Y*Y-5+RANGAU(SEED)*1.Q-3

51    FORMAT('   IER =',I4,5X,'DEGREE =',I3,5X,'NO. OF PTS =',I4,
     1       5X,'ERROR SUM =',1PD14.6/'   COEF. :',4D14.6/(2X,5D14.6))
52    FORMAT('   USING LINL1 :   IER =',I4,5X,
     1       5X,'ERROR SUM =',1PD14.6/'   COEF. :',4D14.6/(2X,5D14.6))

      EPS=1.Q-9
      IG=100
      XL=0.0
      SEED=2

100   PRINT *,'TYPE M=DEGREE,  N=NO. OF DATA PTS'
      PRINT *,'         (QUITS WHEN N.EQ.0)'
      READ *,M,N
      IF(N.LE.0) STOP

!     GENERATING INPUT DATA SET
      H=1.Q0/(N-1.)
      DO 2000 I=1,N
        XI=XL+(I-1)*H
        X(I)=XI
        F(I)=FM(XI)

!	THE BASIS FUNCTION FOR USE WITH LINL1
        G(1,I)=1.0
        DO 2000 J=1,M
          G(J+1,I)=XI**J
2000  CONTINUE

      CALL POLYL1(M,N,B,X,F,EPS,ESUM,IER,WK,IWK)
      WRITE(6,51) IER,M,N,ESUM,(B(I),I=1,M+1)

!	NO. OF BASIS FUNCTIONS
      M1=M+1
      CALL LINL1(M1,N,B,F,G,IG,EPS,ESUM,IER,WK,IWK)
      WRITE(6,52) IER,ESUM,(B(I),I=1,M+1)
      GO TO 100
      END

!     -------------------------------------------

!	To calculate coefficients of linear L_1 approximation in terms of
!	specified basis functions for a tabulated function
!
!	M : (input) Number of basis functions in the required approximation
!	N : (input) Number of points in the table of values. N> M
!	A : (output) Real array of length M+1 containing the coefficients
!		of approximation. A(I) is the coefficient of Phi_I(x) in the
!		approximation
!	F : (input) Real array of length N containing the function values
!	G : (input) Real array of length IG*N, containing the values of
!		basis functions at each point in table. G(I,J) is the
!		value of Phi_I(x) at Jth tabular point.
!	IG : (input) First dimension of array G as declared in the calling
!		program, IG .GE. M
!	EPS : (input) Estimate of roundoff error, passed on to SIMPL1
!	ESUM : (output) L_1 norm of the residual at the calculated approximation
!	IER : (output) error parameter, IER=0 implies successful execution
!		IER=616 implies that M<1  or N<M+1
!			in this case no calculations are done
!		Other values of IER may be set by SIMPL1
!	WK : Real array of length (N+2)*(M+3) used as scratch space
!	IWK : Integer array of length N+M+3 used as scratch space
!
!	Required routines : SIMPL1
!
 
      SUBROUTINE LINL1(M,N,A,F,G,IG,EPS,ESUM,IER,WK,IWK)
      IMPLICIT REAL*16(A-H,O,P,R-Z)
      IMPLICIT LOGICAL(Q)
      DIMENSION A(M+2),F(N),IWK(N+M+3),WK(N+2,M+3),X(N),G(IG,N)
 
      IF(M.LE.0.OR.N.LT.M+1) THEN
        IER=616
        RETURN
      ENDIF
 
      IER=0
      LJ=N+2
      NV=M+1
      M3=N
 
!	Setting up the tableau for simplex algorithm
      DO 2100 J=1,M+2
2100  WK(1,J)=0.0
      DO 2300 I=1,N
        FI=F(I)
        SI=SIGN(1.Q0,FI)
        IWK(I+1)=I*SI
        WK(I+1,1)=FI*SI
        WK(1,1)=WK(1,1)-FI*SI
        S=0.0
        DO 2200 J=1,M
          TI=G(J,I)*SI
          S=S+TI
          WK(I+1,J+1)=TI
          WK(1,J+1)=WK(1,J+1)-TI
2200    CONTINUE
        WK(I+1,M+2)=-S
        WK(1,M+2)=WK(1,M+2)+S
2300  CONTINUE
 
      DO 2400 J=1,M+1
        IWK(N+1+J)=N+J
2400  CONTINUE
 
      NC=NV+M3
      CALL SIMPL1(WK,LJ,NC,M3,IWK,IWK(M3+1),IER,EPS)

!	L_1 norm of the residual
      ESUM=-WK(1,1)
!	Finding the coefficients from the tableau
      DO 3200 I=M3+2,NC+1
        IF(IWK(I).GT.N) A(IWK(I)-N)=0.0
3200  CONTINUE
      DO 3400 I=2,M3+1
        IF(IWK(I).GT.N) A(IWK(I)-N)=WK(I,1)
3400  CONTINUE
      DO 3600 I=1,M
3600  A(I)=A(I)-A(M+1)
      END

!     -------------------------------------------------------------

!	To calculate coefficients of polynomial L_1 approximation
!	for a tabulated function
!
!	M : (input) Required degree of polynomial
!	N : (input) Number of points in the table of values. N> M+1
!	A : (output) Real array of length M+2 containing the coefficients
!		of approximation. A(I+1) is the coefficient of x**I
!	X : (input) Real array of length N containing the points at which
!		function value is available
!	F : (input) Real array of length N containing the function values at X(I)
!	EPS : (input) Estimate of roundoff error, passed on to SIMPL1
!	ESUM : (output) L_1 norm of the residual at the calculated approximation
!	IER : (output) error parameter, IER=0 implies successful execution
!		IER=616 implies that M<0 or N<M+2
!			in this case no calculations are done
!		Other values of IER may be set by SIMPL1
!	WK : Real array of length (N+2)*(M+3) used as scratch space
!	IWK : Integer array of length N+M+3 used as scratch space
!
!	Required routines : SIMPL1
!
      SUBROUTINE POLYL1(M,N,A,X,F,EPS,ESUM,IER,WK,IWK)
      IMPLICIT REAL*16(A-H,O,P,R-Z)
      IMPLICIT LOGICAL(Q)
      DIMENSION A(M+2),F(N),IWK(N+M+3),WK(N+2,M+3),X(N)

      IF(M.LT.0.OR.N.LT.M+2) THEN
        IER=616
        RETURN
      ENDIF

      IER=0
      LJ=N+2
      NV=M+2
!	The number of constraints
      M3=N

!	Setting up the tableau for simplex algorithm
      DO 2100 J=1,M+3
2100  WK(1,J)=0.0
      DO 2300 I=1,N
        FI=F(I)
        SI=SIGN(1.Q0,FI)
        IWK(I+1)=I*SI
        WK(I+1,1)=FI*SI
        WK(1,1)=WK(1,1)-FI*SI
        WK(I+1,2)=SI
        WK(1,2)=WK(1,2)-SI
        S=SI
        DO 2200 J=1,M
          TI=X(I)**J*SI
          S=S+TI
          WK(I+1,J+2)=TI
          WK(1,J+2)=WK(1,J+2)-TI
2200    CONTINUE
        WK(I+1,M+3)=-S
        WK(1,M+3)=WK(1,M+3)+S
2300  CONTINUE

      DO 2400 J=1,M+2
        IWK(N+1+J)=N+J
2400  CONTINUE

      NC=NV+M3
      CALL SIMPL1(WK,LJ,NC,M3,IWK,IWK(M3+1),IER,EPS)

!	L_1 norm of the residual
      ESUM=-WK(1,1)
!	Finding the coefficients from the tableau
      DO 3200 I=M3+2,NC+1
        IF(IWK(I).GT.N) A(IWK(I)-N)=0.0
3200  CONTINUE
      DO 3400 I=2,M3+1
        IF(IWK(I).GT.N) A(IWK(I)-N)=WK(I,1)
3400  CONTINUE
      DO 3600 I=1,M+1
3600  A(I)=A(I)-A(M+2)
      END

!     --------------------------------------------------------------

!	To solve a linear programming problem in the standard form arising
!		 in L_1 minimisation problems using the simplex method
!
!	A : (input/output) Real array of length IA*(N-M+1) containing
!		the tableau of simplex algorithm
!		A(1,I+1)=c_i, the cost coefficients
!		Rows 2 to M+1 contain constraints with A(j,1)=b_j and A(j,i+1)=a_i
!	IA : (input) First dimension of array A as declared in the calling
!		program. IA .GE. M+2
!	N : (input) Number of variables, each is constraint to be .GE.0
!	M : (input) Number of constraints of form a^T X = b_i .GE. 0
!	ID : (input/output) integer array of length M+1 which contains
!		information about interchange of variables on LHS
!	IV : (input/output) integer array of length N-M+1 which contains
!		information about interchange of variables on RHS
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=63 implies that the objective function is unbounded from below
!		IER=635 implies that the simplex algorithm failed to find
!			the optimal feasible vector
!	AEPS : (input) Required accuracy, any coefficient <AEPS, may be
!		assumed to be zero
!
!	Required routines : None
!
      SUBROUTINE SIMPL1(A,IA,N,M,ID,IV,IER,AEPS)
      IMPLICIT REAL*16(A-H,O,P,R-Z)
      IMPLICIT LOGICAL(Q)
      DIMENSION A(IA,N-M+1),ID(M+1),IV(N-M+1)
      PARAMETER(NIT=20)

      JF=1
      M1=M+1
      N1=N-M+1
      IER=0

!	The simplex iteration
      DO 4000 IT=1,NIT*(N+M)
!	Finding the minimum reduced cost coefficient
        RMIN=0.0
        K=0
        DO 2000 J=2,N1
          IF(A(JF,J).LT.RMIN) THEN
            RMIN=A(JF,J)
            K=J
          ELSE IF(IV(J).LE.M.AND.2-A(JF,J).LT.RMIN) THEN
            RMIN=2-A(JF,J)
            K=-J
          ENDIF
2000    CONTINUE
        IF(RMIN.GE.-AEPS) RETURN

!	Finding the pivot element
        K1=ABS(K)
        RMIN=0.0
        L=0
        DO 2400 J=2,M+1
          AJ=A(J,K1)
          IF(K.LT.0) AJ=-AJ
          IF(AJ.GT.AEPS) THEN
            R1=A(J,1)/AJ
            IF(R1.LT.RMIN.OR.L.EQ.0) THEN
              RMIN=R1
              L=J
            ENDIF
          ENDIF
2400    CONTINUE
        IF(L.EQ.0) THEN
!	The objective function is unbounded from below
          IER=63
          RETURN
        ENDIF

        IF(K.LT.0) THEN
          A(JF,K1)=-2+A(JF,K1)
          IV(K1)=-IV(K1)
        ENDIF

!	Exchange the variables
        L1=ID(L)
        ID(L)=IV(K1)
        IV(K1)=L1
        DO 3000 J=1,N1
          IF(J.NE.K1) THEN
            R1=A(L,J)/A(L,K1)
            DO 2800 I=1,M1
              IF(I.NE.L) A(I,J)=A(I,J)-A(I,K1)*R1
2800        CONTINUE
          ENDIF
3000    CONTINUE

        R1=ABS(A(L,K1))
        DO 3200 J=1,N1
          IF(J.NE.K1) A(L,J)=A(L,J)/R1
3200    CONTINUE
        DO 3400 I=1,M1
          IF(I.NE.L) A(I,K1)=-A(I,K1)/A(L,K1)
3400    CONTINUE
        A(L,K1)=1./R1
4000  CONTINUE

!	Iteration fails to converge
      IER=635
      END

!     ---------------------------------------------------

!	To generate random numbers with Gaussian probability distribution
!	It generates random numbers with zero mean and variance of 1.
!	
!	SEED : (input/output) real seed, it should be positive and
!		less than AM. It is updated by the routine and should
!		not be modified between two calls, unless a fresh
!		sequence is required
!
!       THE ARGUMENT OF THIS FUNCTION HAS CHANGED AS COMPARED TO EARLIER
!       VERSION AS THE SEED IS NOW REAL*16 INSTEAD OF INTEGER.
!
!	Required routines : None

      FUNCTION RANGAU(SEED)
      IMPLICIT REAL*16(A-H,O-Z)
      PARAMETER(AM=2147483648Q0,A=45875Q0,AC=453816693Q0,
     1                  AN=2147483647Q0)
      PARAMETER(PI=3.14159265358979323846264338327950288Q0)

      R1=1+MOD(A*SEED+AC,AM)
      IF(SEED.EQ.0.0) SEED=0.1Q0
      RANGAU=SQRT(2.Q0*LOG(AN/SEED))*COS(2.0*PI*R1/AN)
      SEED=R1
      END
