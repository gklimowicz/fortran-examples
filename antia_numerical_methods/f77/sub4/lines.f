!	To solve a system of nonlinear parabolic differential equation using
!	the method of lines
!	The differential equation is assumed to be of the form
!
!	du/dt = f_i( d^2u/dx^2, du/dx, u, x, t),    i=1,2,...,NE
!
!	with Dirichlet boundary conditions
!
!	u(x0,t) = g1(t);  	u(xN,t) = g2(t)
!
!	This subroutine is used through subroutine MSTEP or RKM for solving
!	a system of ordinary differential equations. This routine specifies
!	the equivalent system of ordinary differential equations.
!	The parameters for the partial differential equations are passed
!	through the common block ZZLINE. Elements of this common block must
!	be initialised before calling MSTEP or RKM.
!
!	T : (input) Value of "time" where the derivatives need to be calculated.
!	N : (input) Number of ordinary differential equations to be solved
!		after applying the method of lines, N=NE*(NX-2)
!	U : (input) Real array of length N containing the estimated
!		solution at t.
!	DU : (output) Real array of length N which will contain the derivatives
!		of U at T.
!
!	The arguments in common block are as follows:
!
!	DX : (input) The step length in X to be used for computations.
!		It is assumed to be uniform.
!	X0 : (input) Lower limit on X where the solution is required
!	XN : (input) Upper limit on X where the solution is required.
!		Solution is computed in the interval (X0,XN)
!	X : (input) Real array of length NX-2 containing the mesh points used
!		in X direction. The end points X0 and XN are not included.
!		This array must be set before calling MSTEP or RKM.
!	U0 : (output) Real array of length NE containing the solution
!		at the first boundary X0. It is calculated using the
!		boundary conditions calculated by subroutine BC.
!	UN : (output) Real array of length NE containing the solution
!		at the second boundary XN. It is calculated using the
!		boundary conditions calculated by subroutine BC.
!	UX : Real array of length at least 2*NE used as scratch space
!	NE : (input) Number of differential equations in the parabolic system
!	NX : (input) Number of mesh points in the X direction.
!
!	The dimensions of all arrays in the common block must match those
!	in the calling program and should be sufficient to define the
!	equations. No checks are made and there is no error exit.
!	It is the responsibility of user to ensure consistency.
!
!	This subroutine requires subroutines FCN to calculate the derivatives
!	defining the differential equations and subroutine BC to calculate
!	the boundary conditions. The names of these routines are fixed and
!	cannot be passed on as arguments.
!
!	SUBROUTINE FCN(NE,X,T,U,UX,UXX,DU) and SUBROUTINE BC(NE,T,X0,XN,U0,UN)
!	must be supplied by the user 
!	Subroutine FCN should calculate the derivatives DU = dU/dt at
!	specified T, X, U, UX (=dU/dX), UXX (=d^2U/dX^2)
!	Here NE is the number of parabolic equations in the system, T and X
!	are the values of t and x, U, UX, UXX and DU are real arrays of
!	length NE containing the solution and its derivative as described
!	above.
!
!	Subroutine BC should calculate the solution at boundaries at x=X0 and XN
!	Here NE is the number of parabolic equations in the system,
!	T is the value of t at which boundary conditions are required.
!	X0 and XN are the two boundaries in X where boundary conditions are
!	applied. U0 and UN are real arrays of length NE containing the
!	solution at X0 and XN. U0 and UN must be calculated by subroutine BC.
!
!	Required routines : FCN, BC

      SUBROUTINE LINES(T,N,U,DU)
!      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NP=101,NQ=5)
      COMMON/ZZLINE/DX,X0,XN,X(NP),U0(NQ),UN(NQ),UX(NQ,2),NX,NE
      DIMENSION U(N),DU(N)

!	Calculate the boundary values
      CALL BC(NE,T,X0,XN,U0,UN)
      DO 2000 I=1,N,NE
        J=(I-1)/NE+1

        DO 1500 K=1,NE
          IF(I.EQ.1) THEN
            UM=U0(K)
          ELSE
            UM=U((J-2)*NE+K)
          ENDIF

          IF(I+NE.GT.N) THEN
            UP=UN(K)
          ELSE
            UP=U(J*NE+K)
          ENDIF

          UX(K,1)=(UP-UM)/(2.*DX)
          UX(K,2)=(UP-2.*U(I-1+K)+UM)/DX**2
1500    CONTINUE

        CALL FCN(NE,X(J),T,U(I),UX,UX(1,2),DU(I))
2000  CONTINUE
      END
