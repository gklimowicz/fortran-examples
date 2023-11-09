!	Subroutine to calculate the Chi square function for a nonlinear
!	least squares fit, for use with subroutine BFGS.
!
!	N : (input) Number of parameters to be fitted
!	A : (input) Real array of length N containing the parameter values
!		at which the function needs to be evaluated
!	F : (output) The function value at specified parameter values
!	G : (output) Real array of length N containing the gradient vector
!		G(I) will contain dF/dA(I)
!
!	The data points are passed through common block ZZFUN, which
!	must be initialised in the calling program before calling BFGS
!
!	FX : (input) Real array of length NP containing the function values
!	X : (input) Real array of length NP containing the abscissas at which
!		function values are available
!	EF : (input) Real array of length NP containing the estimated errors in 
!		function values, for use in deciding the weights for each point
!		Although EF should contain the error estimate, but in many
!		cases it is found that multiplying all values by a suitable
!		constant can improve the convergence of BFGS dramatically
!	FX1 : (output) Real array of length NP containing the fitted value
!		of the function at each tabular point
!	NN : (input) The number of data points in the table to be fitted
!
!	The parameter NP must be equal to the dimension of arrays as
!	declared in the calling program.
!
!	This routine requires subroutine FCN to calculate the required
!	function which has to be fitted. There is no provision to pass on
!	the name of this subroutine and hence it must be changed explicitly
!	to the required name. SUBROUTINE FCN(N,A,X,F,DF) must be supplied
!	by the user. Here N is the number of parameters, A is a real array
!	of length N containing the values of parameters and X is the value
!	of independent variable where the fitting function needs to be evaluated.
!	F is the calculated function value and DF is a real array of length N
!	containing the calculated derivatives.
!
!	Required routines : FCN

      SUBROUTINE NLLSQ(N,A,F,G)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NP=100)
      COMMON/ZZFUN/FX(NP),X(NP),EF(NP),FX1(NP),NN
      DIMENSION A(N),G(N),DF(NP)
 
!     GENERATING THE FUNCTION FOR MINIMISATION
 
      IF(N.GT.NP) STOP 129
      F=0.0
      DO 1000 I=1,N
1000  G(I)=0.0
      DO 3000 I=1,NN
        CALL FCN(N,A,X(I),FA,DF)
        FX1(I)=FA
 
!     SUM OF SQUARED DIFFERENCE AND ITS GRADIENT
 
        F=F+((FX(I)-FA)/EF(I))**2
        DO 2400 J=1,N
          G(J)=G(J)+2*(FA-FX(I))*DF(J)/EF(I)**2
2400    CONTINUE
3000  CONTINUE
      END
