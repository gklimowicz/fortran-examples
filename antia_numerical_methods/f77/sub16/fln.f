!	Function routine to calculate the function value as required for
!	line search without derivatives
!
!	FCN : (input) Name of subroutine to calculate the required function
!	X : (input) Parameter along the line to specifying the point where
!		function evaluation is required
!	V : (input/output) Real array of length 2N, first N elements specify the
!		direction of line search. After execution next N elements will contain
!		the coordinates of the point at which function is evaluated.
!	X0 : (input) Real array of length N, containing the coordinates
!		of the starting point for line search
!	N : (input) Number of variables in the function to be minimised
!	NUM : (input/output) Integer variable to keep count of function evaluations
!
!	SUBROUTINE FCN(N,X,FX) to calculate the required function, must be supplied
!		by the user. Here N is the number of variables, FX is the
!		function value at X. X is a real array of length N.
!
!	Required routines : FCN

      FUNCTION FLN(FCN,X,V,X0,N,NUM)
      IMPLICIT REAL*16(A-H,O-Z)
      DIMENSION V(N),X0(2*N)

      NUM=NUM+1
!	coordinates of the required points
      DO 1000 I=1,N
1000  X0(N+I)=X0(I)+V(I)*X
      CALL FCN(N,X0(N+1),FLN)
      END
