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
!      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL EQN,BCS
      PARAMETER(NIT=20,EPS=1.D-30)
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
        IF(ABS(H-HH).GT.1.D-4*ABS(HH)) THEN
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
!       IMPLICIT REAL*4(A,B,D-H,O-Z)
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
!       IMPLICIT REAL*4(A,B,D-H,O-Z)
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
