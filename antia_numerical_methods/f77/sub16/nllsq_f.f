!	Subroutine to calculate the Chi square function for a nonlinear
!	least squares fit, for use with subroutine NMINF.
!	Version of NLLSQ for use with NMINF, when derivatives are not available
!
!	N : (input) Number of parameters to be fitted
!	A : (input) Real array of length N containing the parameter values
!		at which the function needs to be evaluated
!	F : (output) The function value at specified parameter values
!
!	The data points are passed through common block ZZFUN, which
!	must be initialised in the calling program before calling NMINF
!
!	FX : (input) Real array of length NP containing the function values
!	X : (input) Real array of length NP containing the abscissas at which
!		function values are available
!	EF : (input) Real array of length NP containing the estimated errors in 
!		function values, for use in deciding the weights for each point
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
!	to the required name. SUBROUTINE FCN(N,A,X,F) must be supplied
!	by the user. Here N is the number of parameters, A is a real array
!	of length N containing the values of parameters and X is the value
!	of independent variable where the fitting function needs to be evaluated.
!	F is the calculated function value.
!
!	Required routines : FCN

      SUBROUTINE NLLSQ_F(N,A,F)
      IMPLICIT REAL*16(A-H,O-Z)
      PARAMETER(NP=100)
      COMMON/ZZFUN/FX(NP),X(NP),EF(NP),FX1(NP),NN
      DIMENSION A(N)
 
!     GENERATING THE FUNCTION FOR MINIMISATION
 
      F=0.0
      DO 3000 I=1,NN
        CALL FCN(N,A,X(I),FA)
        FX1(I)=FA
 
!     SUM OF SQUARED DIFFERENCE
 
        F=F+((FX(I)-FA)/EF(I))**2
3000  CONTINUE
      END
