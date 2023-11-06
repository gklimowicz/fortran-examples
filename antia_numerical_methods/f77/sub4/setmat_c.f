!	To setup the matrix of system of linear equations arising from
!	finite difference approximation of ordinary differential equations
!	This routine is called by GEVP_C, the matrix is complex
!
!	N : (input) Number of mesh points used.
!	M : (input) Number of first order differential equations in the system
!	ML : (input) Number of boundary conditions at the first boundary t=T(1)
!	A : (output) Complex array of length (M+ML)*2M*N containing
!		the matrix of equations.
!	BC : (output) Complex array of length (M+ML)*(M+1) containing the
!		coefficients of boundary conditions
!	X : (input) Complex array of length M*N containing the current approximation
!		to the solution.
!	XC : (output) Complex array of length M*N containing the right hand
!		side of finite difference equations calculated by the routine
!	T : (input) Real array of length N containing the mesh points.
!	PAR : (input) Complex array containing the parameters to be passed
!		on to subroutines EQN and BCS. This array is not used by
!		the subroutine.
!	EQN : (input) Name of subroutine to calculate the equation matrix
!	BCS : (input) Name of subroutine to calculate the boundary conditions
!
!	SUBROUTINE EQN(J,M,ML,PAR,A,B,Y,F,T) and
!	SUBROUTINE BCS(M,ML,PAR,BC,G,T1,TN,Y1,YN) must be supplied by the user.
!	The form of these routines is described in documentation for GEVP_C
!
!	Required routines : EQN, BCS
!
      SUBROUTINE SETMAT_C(N,M,ML,A,BC,X,XC,T,PAR,EQN,BCS)
!      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT REAL*4(H,O,P,R,T)
      IMPLICIT COMPLEX*8(A-G,U-Z)
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
