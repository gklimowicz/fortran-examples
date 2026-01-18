!     PROGRAM TO SOLVE BOUNDARY VALUE PROBLEMS IN 
!     ORDINARY DIFFERENTIAL EQUATIONS USING FINITE DIFFERENCE METHOD

      PROGRAM BVP
      IMPLICIT REAL*16(A-H,O-Z)
      EXTERNAL EQN,BCS
      DIMENSION PAR(10),X(2,801),WK(12500),IWK(1700),T(801),XC(2,801)

!     EXAMPLE 12.13

51    FORMAT('   IER =',I4,5X,'NO. OF PTS =',I5,5X,'LAMDA =',1PD14.6/
     1       11X,'T',20X,'SOLUTION',18X,'CORRECTED SOLUTION')
52    FORMAT(I4,2X,1PD14.6,2X,2D14.6,2X,2D14.6)

      M=2
      ML=1
      IFLAG=0
      REPS=1.Q-14

!     PASS THE PARAMETER LAMDA IN THE EQUATION VIA THE ARRAY PAR

100   PRINT *,'TYPE N=NO. OF PTS,  LAMDA        (QUITS WHEN N.LE.0)'
      READ *,N,PAR(1)
      IF(N.LE.0) STOP

!     SET UP THE MESH WITH UNIFORM SPACING AND THE INITIAL VALUES

      H=1.Q0/(N-1)
      DO 500 I=1,N
        X(1,I)=1.
        X(2,I)=1.
500   T(I)=(I-1)*H
      CALL FDM(N,M,ML,PAR,X,XC,T,EQN,BCS,IWK,WK,IFLAG,REPS,IER)
      WRITE(6,51)IER,N,PAR(1)
      DO 1000 I=1,N,10
        WRITE(6,52) I,T(I),X(1,I),X(2,I),XC(1,I),XC(2,I)
1000  CONTINUE
      GO TO 100
      END
 
!     ------------------------------------------------------

!	To solve two-point boundary value problem in ordinary differential
!	equations using finite difference method
!
!	N : (input) Number of mesh points to be used.
!	M : (input) Number of first order differential equations in the system
!	ML : (input) Number of boundary conditions at the first boundary t=T(1)
!	PAR : (input) Real array to be passed on to EQN and BCS for calculating
!		the equations and boundary conditions. This array is not used
!		by FDM, but is merely used to pass on any extra parameters
!		that may be required to specify the equations
!	X : (input/output) Real array of length M*N containing the solution.
!		It should contain the initial guess to solution at the time
!		of calling. After execution it will contain the calculated
!		solution. X(i,j) is the ith component of solution at jth mesh point.
!		First dimension of X in calling program must be M
!	XC : (output) Real array of length M*N containing the solution after
!		applying the deferred correction. It is stored in the same
!		format as X.
!	T : (input) Real array of length N containing the mesh points.
!		These points must be in ascending or descending order.
!		For calculating the deferred correction the mesh spacing
!		must be uniform.
!	EQN : (input) Name of subroutine to specify the differential equation
!		to be solved.
!	BCS : (input) Name of subroutine to calculate the boundary conditions
!		at t=T(1) and T(N)
!		Outline of a sample routine for EQN and BCS can be found
!		at the end of this file
!	IWK : Integer array of length M*N used as scratch space
!	WK : Real array of length (M+ML)*2M*(N+1) used as scratch space
!	IFLAG : (input) Integer variable used as a flag to decide the type of
!		computation required.
!		IFLAG=0 implies that equations are nonlinear and deferred
!			correction is to be calculated. In this case X must
!			contain the initial guess and mesh spacing must be uniform
!		IFLAG=1 implies that equations are nonlinear and deferred
!			correction is not required. In this case X must contain
!			the initial guess and mesh spacing could be arbitrary
!		IFLAG=2 implies that equations are linear and deferred correction
!			is required. In this case initial guess is not required,
!			but mesh spacing must be uniform
!		IFLAG=3 implies that equations are linear and deferred correction
!			is not required. In this case initial guess is not
!			required and mesh spacing can be arbitrary
!	REPS : (input) Required accuracy. This is only used to check convergence
!		of Newton's method for solving finite difference equations.
!		The truncation error depends on mesh spacing and is not
!		controlled in this routine.
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=704 implies that N<3, M.LE.ML, or ML.LE.0, in which case
!			no calculations are done
!		IER=734 implies that N<5 and deferred correction is requested.
!			In this case deferred correction is not calculated
!		IER=735 implies that the finite difference matrix is singular
!		IER=736 implies that mesh spacing is not uniform and
!			deferred correction is not calculated
!		IER=737 implies that Newton's iteration for solving the
!			finite difference equations failed to converge.
!
!	SUBROUTINE EQN and BCS must be supplied by the user.
!	SUBROUTINE EQN(J,M,ML,PAR,A,B,Y,F,T) calculates the right hand
!		sides for differential equations By'_i=f_i(t,y,par)
!		J is the serial number of mesh point at which calculation
!		is required. M is the number of first order differential
!		equations in the system, ML the number of boundary conditions
!		at the first point, PAR is a real array which can be used
!		to pass on any required parameters to define the equations.
!		A and B are real arrays of length (M+ML)*M defining the
!		differential equation B Y'=f(T,PAR,Y) and
!		A(I,K)=dF_I/dY(K) is the Jacobian matrix
!		Y is real array of length M specifying the solution at t=T
!		F is a real array of length M containing the right hand
!		sides of differential equations f_I as defined above.
!		F, A and B must be calculated by the subroutine.
!
!	SUBROUTINE BCS(M,ML,PAR,BC,G,T1,TN,Y1,YN) calculates the boundary
!		conditions at both boundaries. M is the number of differential
!		equations, ML is the number of boundary condition at the
!		first mesh point t=T1. PAR is a real array which can be used
!		to pass any required parameters. BC is a real array of length
!		(M+ML)*M which should contain the coefficients of boundary
!		conditions. First ML rows will specify the boundary conditions
!		at t=T1, while remaining rows will specify those at t=TN.
!		G is a real array of length M specifying the boundary conditions
!		G(I)=g_i(T1,PAR,Y1) (I.LE.ML) are the boundary conditions
!		at T1 (g_i=0), while G(I)=g_i(TN,PAR,YN) (I.GT.ML) are the
!		boundary conditions at TN. BC is the Jacobian matrix dg_i/dY(K)
!		Y1 and YN are real arrays of length M specifying the
!		solution at t=T1 and TN.
!
!	Required routines : SETMAT, GAUBLK, EQN, BCS
!
      SUBROUTINE FDM(N,M,ML,PAR,X,XC,T,EQN,BCS,IWK,WK,IFLAG,REPS,IER)
      IMPLICIT REAL*16(A-H,O-Z)
      EXTERNAL EQN,BCS
      PARAMETER(NIT=20,EPS=1.Q-30)
      DIMENSION PAR(*),X(M,N),XC(M,N),IWK(M,N),WK(M+ML,2*M,N+1),T(N)

      IF(N.LT.3.OR.M.LE.ML.OR.ML.LE.0) THEN
        IER=704
        RETURN
      ENDIF

      DO 1100 I=1,NIT
!	Set up the matrix of finite difference equations
        CALL SETMAT(N,M,ML,WK,WK(1,1,N+1),X,XC,T,PAR,EQN,BCS)
        IFLG=0
!	Solve the system of linear equations
        CALL GAUBLK(N,M,ML,WK,IFLG,DET,IDET,IWK,XC,IER)
        IF(IER.GT.0) RETURN

!	Checking for convergence
        RERR=0.0
        DO 1000 J=1,N
          J1=J+1
          IF(J.EQ.N) J1=J-1
          DO 1000 K=1,M
            XJ=X(K,J)+XC(K,J)
            R2=ABS(XJ)+ABS(X(K,J)-X(K,J1))+EPS
            RE=ABS(XC(K,J)/R2)
            IF(RE.GT.RERR) RERR=RE
            X(K,J)=XJ
1000    CONTINUE
        IF(RERR.LT.REPS.OR.IFLAG.GT.1) GO TO 1150
1100  CONTINUE

!	Newton's iteration fails to converge
      IER=737
      RETURN

1150  IF(IFLAG.EQ.1.OR.IFLAG.GT.2) RETURN
      IF(N.LE.4) THEN
        IER=734
        RETURN
      ENDIF

!	Calculate the deferred correction
      HH=T(2)-T(1)
      DO 2000 J=1,N-1
        TJ=0.5*(T(J)+T(J+1))
        H=T(J+1)-T(J)
        IF(ABS(H-HH).GT.1.Q-4*ABS(HH)) THEN
!	Mesh spacing is not uniform, hence quit
          IER=736
          RETURN
        ENDIF

        JJ=J
        CALL EQN(JJ,M,ML,PAR,WK(1,1,N+1),WK(1,M+1,N+1),X(1,J),XC(1,J+1),
     1           TJ)
        DO 1200 I=1,M
          IF(J.EQ.1) THEN
            WK(M+1,I,N+1)=(-2*X(I,1)+7*X(I,2)-9*X(I,3)+5*X(I,4)-X(I,5))
     1                     /24.
            WK(M+1,I+M,N+1)=(43*X(I,1)-112*X(I,2)+102*X(I,3)-40*X(I,4)+
     1                       7*X(I,5))*H/192.
          ELSE IF(J.EQ.N-1) THEN
            WK(M+1,I,N+1)=(2*X(I,N)-7*X(I,N-1)+9*X(I,N-2)-5*X(I,N-3)
     1                     +X(I,N-4))/24.
            WK(M+1,I+M,N+1)=(43*X(I,N)-112*X(I,N-1)+102*X(I,N-2)-
     1                       40*X(I,N-3)+7*X(I,N-4))*H/192.
          ELSE
            WK(M+1,I,N+1)=(-X(I,J-1)+3.*X(I,J)-3.*X(I,J+1)+X(I,J+2))/24.
            WK(M+1,I+M,N+1)=(X(I,J-1)-X(I,J)-X(I,J+1)+X(I,J+2))*H/16.
          ENDIF
1200    CONTINUE

!	Set up the RHS for deferred correction
        DO 1500 I=1,M
          XD=0.0
          DO 1400 K=1,M
           XD=XD+WK(I,K+M,N+1)*WK(M+1,K,N+1)-WK(I,K,N+1)*WK(M+1,K+M,N+1)
1400      CONTINUE
          XC(ML+I,J)=XD
1500    CONTINUE
2000  CONTINUE

      DO 2200 I=1,ML
2200  XC(I,1)=0.0
      DO 2300 I=ML+1,M
2300  XC(I,N)=0.0

!	Calculate deferred correction
      CALL GAUBLK(N,M,ML,WK,IFLG,DET,IDET,IWK,XC,IER)
      DO 2500 J=1,N
        DO 2500 I=1,M
2500  XC(I,J)=XC(I,J)+X(I,J)
      END

!	--------------------------------
!	
!	SUBROUTINE EQN(J,M,ML,PAR,A,B,Y,F,T)
!       IMPLICIT REAL*16(A,B,D-H,O-Z)
!	DIMENSION A(M+ML,M),B(M+ML,M),Y(M),PAR(*),F(M)
!
!	DO 1000 I=1,M
!	F(I)=f_i(T,PAR,Y)
!	  DO 1000 K=1,M
!	    A(K,I)=df_k/dY(i)
!	    B(K,I)=b_{ki}(T,PAR)
! 1000  CONTINUE
!       END
!
!	SUBROUTINE BCS(M,ML,PAR,BC,G,T1,TN,Y1,YN)
!       IMPLICIT REAL*16(A,B,D-H,O-Z)
!	DIMENSION PAR(*),BC(M+ML,M),G(M),Y1(M),YN(M)
!
!	DO 1000 I=1,M
!	  IF(I.LE.ML) THEN
!	    G(I)=g_i(T1,PAR,Y1)
!	  ELSE
!	    G(I)=g_i(TN,PAR,YN)
!         ENDIF
!	  DO 1000 K=1,M
!	    BC(I,K)=dg_I/dY(K)
! 1000  CONTINUE
!       END

!     ----------------------------------------------------------

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
      IMPLICIT REAL*16(A-H,O-Z)
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

!     -------------------------------------------------------

!	To setup the matrix of system of linear equations arising from
!	finite difference approximation of ordinary differential equations
!	This routine is called by FDM or GEVP
!
!	N : (input) Number of mesh points used.
!	M : (input) Number of first order differential equations in the system
!	ML : (input) Number of boundary conditions at the first boundary t=T(1)
!	A : (output) Real array of length (M+ML)*2M*N containing
!		the matrix of equations.
!	BC : (output) Real array of length (M+ML)*(M+1) containing the
!		coefficients of boundary conditions
!	X : (input) Real array of length M*N containing the current approximation
!		to the solution.
!	XC : (output) Real array of length M*N containing the right hand
!		side of finite difference equations calculated by the routine
!	T : (input) Real array of length N containing the mesh points.
!	PAR : (input) Real array containing the parameters to be passed
!		on to subroutines EQN and BCS. This array is not used by
!		the subroutine.
!	EQN : (input) Name of subroutine to calculate the equation matrix
!	BCS : (input) Name of subroutine to calculate the boundary conditions
!
!	SUBROUTINE EQN(J,M,ML,PAR,A,B,Y,F,T) and
!	SUBROUTINE BCS(M,ML,PAR,BC,G,T1,TN,Y1,YN) must be supplied by the user.
!	The form of these routines is described in documentation for FDM or GEVP
!
!	Required routines : EQN, BCS
!
      SUBROUTINE SETMAT(N,M,ML,A,BC,X,XC,T,PAR,EQN,BCS)
      IMPLICIT REAL*16(A-H,O-Z)
      DIMENSION A(M+ML,2*M,N),BC(M+ML,M+1),X(M,N),XC(M,N),PAR(*),T(N)

!	Loop over the mesh points
      DO 1500 K=1,N-1
!	t_{k+1/2}
        TK=0.5*(T(K)+T(K+1))
        H=T(K+1)-T(K)
        DO 800 I=1,M
!	y_{k+1/2}
800     BC(I,1)=0.5*(X(I,K)+X(I,K+1))
        KK=K
!	Calculate the equation matrix at t_{k+1/2}
        CALL EQN(KK,M,ML,PAR,A(1,1,K),A(1,M+1,K),BC,XC(ML+1,K),TK)

!	Setup the RHS of finite difference equations
        DO 1100 J=1,M
          XK=XC(ML+J,K)*H
          DO 1000 I=1,M
1000      XK=XK-A(J,M+I,K)*(X(I,K+1)-X(I,K))
          XC(ML+J,K)=XK
1100    CONTINUE

!	Setup the finite difference matrix
        DO 1500 J=1,M
          DO 1200 I=M,1,-1
            A(ML+I,J,K)=-A(I,J+M,K)-0.5*H*A(I,J,K)
1200      A(ML+I,J+M,K)=A(I,J+M,K)-0.5*H*A(I,J,K)
          DO 1500 I=1,ML
1500  A(I,J+M,K)=0.0

!	The boundary conditions
      CALL BCS(M,ML,PAR,BC,BC(1,M+1),T(1),T(N),X(1,1),X(1,N))

!	Boundary conditions at the first boundary
      DO 3000 I=1,ML
        XC(I,1)=-BC(I,M+1)
        DO 3000 J=1,M
3000  A(I,J,1)=BC(I,J)

!	Boundary conditions at the second boundary
      DO 4000 I=ML+1,M
        XC(I,N)=-BC(I,M+1)
        DO 4000 J=1,M
4000  A(I,J,N)=BC(I,J)
      END

!     -------------------------------------------

      SUBROUTINE EQN(J,M,ML,PAR,A,B,X,F,T)
      IMPLICIT REAL*16(A-H,O-Z)
      DIMENSION A(M+ML,M),B(M+ML,M),X(M),PAR(*),F(*)

!     THE DIFFERENTIAL EQ.  X1'=X2,  X2'=LAMDA*SINH(LAMDA*X1)
!     WITH LAMDA=PAR(1)

      DO 1000 I=1,M
        F(I)=0.0
        DO 1000 K=1,M
          A(K,I)=0.0
1000  B(K,I)=0.0

      B(1,1)=1.
      A(1,2)=1.
      B(2,2)=1.
      A(2,1)=PAR(1)**2*COSH(PAR(1)*X(1))
      
      F(1)=X(2)
      F(2)=PAR(1)*SINH(PAR(1)*X(1))
      END

!     ----------------------------------------------

      SUBROUTINE BCS(M,ML,PAR,BC,G,T1,T2,X1,X2)
      IMPLICIT REAL*16(A-H,O-Z)
      DIMENSION PAR(*),BC(M+ML,M),G(M),X1(*),X2(*)

!     THE BOUNDARY CONDITIONS  X1(T=T1)=0,  X1(T=T2)=1

      BC(1,1)=1.0
      BC(1,2)=0.0
      BC(2,1)=1.0
      BC(2,2)=0.0

      G(1)=X1(1)
      G(2)=X2(1)-1.
      END
