!	To solve two-point boundary value problem in ordinary differential
!	equations using expansion method with B-spline basis functions
!
!	NK : (input) Number of knots to be used for calculating B-splines.
!	K : (input) Order of B-splines to be used, K=4 for cubic B-splines
!	M : (input) Number of first order differential equations in the system
!	ML : (input) Number of boundary conditions at the first boundary t=T(1)
!	PAR : (input) Real array to be passed on to EQN and BCS for calculating
!		the equations and boundary conditions. This array is not used
!		by BSPODE, but is merely used to pass on any extra parameters
!		that may be required to specify the equations
!	X : (output) Real array of length M*N containing the solution.
!		X(i,j) is the ith component of solution at jth mesh point
!		TX(j). First dimension of X in calling program must be M
!	A : (input/output) Real array of length (NK+K-2)*M containing the
!		coefficients of expansion. At the time of calling it should
!		contain the initial guess. After execution it will contain
!		the final values. The first dimension of A must be NK+K-2.
!	T : (input) Real array of length NK containing the knots.
!		These points must be in ascending order with T(1) and T(NK)
!		as the two boundaries. 
!	N : (input) The number of mesh points to be used for calculating
!		the coefficients. The solution will be calculated at
!		all these points. N.GE.NK+K-2
!	TX : (input/output) Real array of length N containing the mesh points
!		For IFLAG>1 the mesh points must be supplied, while for
!		other values of IFLAG the routine calculates the mesh point
!		assuming uniform spacing. 
!	EQN : (input) Name of subroutine to specify the differential equation
!		to be solved.
!	BCS : (input) Name of subroutine to calculate the boundary conditions
!		at t=T(1) and T(NK)
!	WK : Real array of length M*(N+1)*(M*(NK+K-2)+7)+(M*(NK+K-2))**2
!		used as scratch space
!	IFLAG : (input) Integer variable used as a flag to decide the type of
!		computation required.
!		IFLAG=0 implies that equations are nonlinear and mesh points
!			are not supplied in array TX. These would be calculated
!			assuming uniform spacing.
!		IFLAG=1 implies that equations are linear and mesh points
!			are not supplied in array TX. These would be calculated
!			assuming uniform spacing.
!		IFLAG=2 implies that equations are nonlinear and mesh points
!			are supplied in array TX.
!		IFLAG=3 implies that equations are linear and mesh points
!			are supplied in array TX.
!	REPS : (input) Required accuracy. This is only used to check convergence
!		of Newton's method for solving the resulting system of equations.
!		The truncation error depends on NK and K and is not
!		controlled in this routine.
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=705 implies that NK<3, M.LE.ML, ML.LE.0, or N.LT.NK+K-2
!			in which case no calculations are done
!		IER=741 implies that Newton's iteration for solving the
!			system of equations failed to converge.
!		Other values may be set by BSPLIN or SVD
!
!	SUBROUTINE EQN and BCS must be supplied by the user.
!	SUBROUTINE EQN(J,M,ML,PAR,A,B,Y,F,T) calculates the right hand
!		sides for differential equations y'_i=f_i(t,y,par)
!		J is the serial number of mesh point at which calculation
!		is required. M is the number of first order differential
!		equations in the system, ML the number of boundary conditions
!		at the first point, PAR is a real array which can be used
!		to pass on any required parameters to define the equations.
!		A and B are real arrays of length M*M defining the
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
!		M*M which should contain the coefficients of boundary
!		conditions. First ML rows will specify the boundary conditions
!		at t=T1, while remaining rows will specify those at t=TN.
!		G is a real array of length M specifying the boundary conditions
!		G(I)=g_i(T1,PAR,Y1) (I.LE.ML) are the boundary conditions
!		at T1 (g_i=0), while G(I)=g_i(TN,PAR,YN) (I.GT.ML) are the
!		boundary conditions at TN. BC is the Jacobian matrix dg_i/dY(K)
!		Y1 and YN are real arrays of length M specifying the
!		solution at t=T1 and TN.
!
!	Required routines : BSPLIN, BSPEVL, SVD, SVDEVL, EQN, BCS
!
 
      SUBROUTINE BSPODE(NK,K,M,ML,PAR,X,A,T,N,TX,EQN,BCS,WK,IFLAG,
     1            REPS,IER)
!      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NIT=50,EPS=1.D-30)
      DIMENSION PAR(*),X(M,N),T(NK),TX(N),WK(M,N+1,*),A(NK+K-2,M)
 
      NB=NK+K-2
      IER=0
      IF(NK.LT.3.OR.M.LE.ML.OR.ML.LE.0.OR.N.LT.NB) THEN
        IER=705
        RETURN
      ENDIF
 
      IF(IFLAG.LE.1) THEN
!	Setup the mesh assuming uniform spacing
        H=(T(NK)-T(1))/(N-1)
        DO 2000 I=1,N
          TX(I)=T(1)+(I-1)*H
2000    CONTINUE
      ENDIF
      NE=N*M+M
      ME=NB*M
 
!	The iteration loop
      NDERIV=1
      DO 5000 IT=1,NIT
 
!	Setup the equation matrix
        DO 3000 I=1,N
          XB=TX(I)
          CALL BSPLIN(T,NK,K,XB,NDERIV,WK(1,1,ME+1),WK(1,1,ME+2),
     1                WK(1,1,ME+3),LEFT,IER,WK(1,1,ME+4))
          IF(IER.GT.100) RETURN
          DO 2500 J=1,M
            X(J,I)=BSPEVL(NK,T,K,NDERIV,A(1,J),XB,DF,DDF,WK(1,1,ME+3)
     1            ,IER)
            IF(IER.GT.100) RETURN
            X(J,I+1)=DF
2500      CONTINUE
          JI=I
          CALL EQN(JI,M,ML,PAR,WK(1,1,ME+3),WK(1,1,ME+4),X(1,I),
     1                WK(1,1,ME+5),XB)
 
          DO 2800 J=1,M
            S=0
            DO 2600 J1=1,M
              S=S+WK(J,J1,ME+4)*X(J1,I+1)
2600        CONTINUE
            WK(J,I,ME+6)=-S+WK(J,1,ME+5)
            DO 2700 J1=1,M*NB
              JK=(J1-1)/NB+1
              JI=J1-(JK-1)*NB
              IF(JI.EQ.0) JI=NB
              WK(J,I,J1)=WK(J,JK,ME+4)*WK(JI,1,ME+2)-WK(J,JK,ME+3)*
     1                   WK(JI,1,ME+1)
2700        CONTINUE
 
2800      CONTINUE
3000    CONTINUE
 
        T1=T(1)
        CALL BSPLIN(T,NK,K,T1,NDERIV,WK(1,1,ME+1),WK(1,1,ME+2),
     1              WK(1,1,ME+3),LEFT,IER,WK(1,1,ME+4))
        IF(IER.GT.100) RETURN
        T2=T(NK)
        CALL BSPLIN(T,NK,K,T2,NDERIV,WK(1,1,ME+2),WK(1,1,ME+3),
     1              WK(1,1,ME+4),LEFT,IER,WK(1,1,ME+5))
        IF(IER.GT.100) RETURN
 
        CALL BCS(M,ML,PAR,WK(1,1,ME+3),WK(1,1,ME+4),T1,T2,X(1,1),X(1,N))
 
        DO 3500 I=1,M
          WK(I,N+1,ME+6)=-WK(I,1,ME+4)
          DO 3200 J1=1,M*NB
            JK=(J1-1)/NB+1
            JI=J1-(JK-1)*NB
            IF(JI.EQ.0) JI=NB
            IF(I.LE.ML) THEN
              WK(I,N+1,J1)=WK(I,JK,ME+3)*WK(JI,1,ME+1)
            ELSE
              WK(I,N+1,J1)=WK(I,JK,ME+3)*WK(JI,1,ME+2)
            ENDIF
3200      CONTINUE
3500    CONTINUE
 
!	Solve the system of equations using SVD
        CALL SVD(ME,NE,WK,WK(1,1,ME+7),WK(1,1,ME+1),NE,ME,WK(1,1,ME+2),
     1          IER)
        IF(IER.GT.100) RETURN
 
        CALL SVDEVL(ME,NE,WK,WK(1,1,ME+7),WK(1,1,ME+1),NE,ME,
     1              WK(1,1,ME+6),WK(1,1,ME+2),REPS)
 
!	Convergence check
        RERR=0.0
        DO 4000 J=1,NB
          J1=J+1
          IF(J.EQ.NB) J1=J-1
          DO 4000 I=1,M
            JK=J+(I-1)*NB
            XJ=A(J,I)+WK(JK,1,ME+6)
            R2=ABS(XJ)+ABS(X(J,I)-X(J1,I))+EPS
            RE=ABS(WK(JK,1,ME+6)/R2)
            IF(RE.GT.RERR) RERR=RE
            A(J,I)=XJ
4000    CONTINUE
        IF(RERR.LT.REPS.OR.IFLAG.EQ.1.OR.IFLAG.EQ.3) RETURN
 
5000  CONTINUE
      IER=741
      RETURN
 
      END

!	--------------------------------
!	
!	SUBROUTINE EQN(J,M,ML,PAR,A,B,Y,F,T)
!       IMPLICIT REAL*4(A,B,D-H,O-Z)
!	DIMENSION A(M,M),B(M,M),Y(M),PAR(*),F(M)
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
!	DIMENSION PAR(*),BC(M,M),G(M),Y1(M),YN(M)
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
