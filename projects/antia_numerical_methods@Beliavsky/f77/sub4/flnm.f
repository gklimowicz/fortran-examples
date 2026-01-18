!	Function routine to calculate the function value and its derivative
!	as required for line search
!
!	FCN : (input) Name of subroutine to calculate the required function
!	X : (input) Parameter along the line to specify the point where
!		function evaluation is required
!	DF : (output) First derivative of function along the line at X
!	V : (input/output) Real array of length 3N, first N elements specify the
!		direction of line search. Next N elements will contain the
!		coordinates of the point at which function is evaluated,
!		while the last N elements contain the gradient vector at the point
!	X0 : (input) Real array of length N, containing the coordinates
!		of starting point for line search
!	N : (input) Number of variables in the function to be minimised
!	NUM : (input/output) Integer variable to keep count of function evaluations
!
!	SUBROUTINE FCN(N,X,FX,G) to calculate the required function, must be supplied
!		by the user. Here N is the number of variables, FX is the
!		function value at X and G is the gradient vector. X and G
!		are real arrays of length N.
!
!	Required routines : FCN

      FUNCTION FLNM(FCN,X,DF,V,X0,N,NUM)
!      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION V(3*N),X0(N)

      NUM=NUM+1
      N2=2*N
!	The coordinates of the required point
      DO 1000 I=1,N
1000  V(N+I)=X0(I)+V(I)*X

      CALL FCN(N,V(N+1),FLNM,V(N2+1))
!	The first derivative along the search direction
      DF=0.0
      DO 2000 I=1,N
2000  DF=DF+V(I)*V(N2+I)
      END
