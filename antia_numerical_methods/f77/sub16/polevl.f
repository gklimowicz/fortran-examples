!	Evaluating the fitted polynomial and its derivatives at any value
!	of x using known coefficients of orthogonal polynomials
!	Should be used to evaluate the polynomial using coefficients calculated
!	by POLFIT.
!
!	M : (input) Degree of polynomial
!	A : (input) Real array of length M+1 containing the coefficients
!		of the fit
!	ALP, BETA : (input) Real arrays of length M+1, containing the coefficients
!		required for defining the orthogonal polynomials
!		A, ALP, BETA could be calculated using POLFIT
!	X : (input) Value of x at which polynomial needs to be evaluated
!	F : (output) Calculated value of the polynomial at X
!	DF : (output) First derivative of F(X) at X
!	DDF : (output) Second derivative of F(X) at X
!	
!	Required routines : None

      SUBROUTINE POLEVL(M,A,ALP,BETA,X,F,DF,DDF)
      IMPLICIT REAL*16(A-H,O-Z)
      DIMENSION A(M+1),ALP(M+1),BETA(M+1)

      F=A(M)+(X-ALP(M))*A(M+1)
      F1=A(M+1)
      DF=A(M+1)
      DF1=0.0
      DDF=0.0
      DDF1=0.0
      IF(M.LE.1) RETURN

!	Clenshaw's recurrence for F, DF and DDF
      DO 1000 J=M-2,0,-1
        DD=2.*DF+(X-ALP(J+1))*DDF-BETA(J+2)*DDF1
        D=F+(X-ALP(J+1))*DF-BETA(J+2)*DF1
        FF=A(J+1)+(X-ALP(J+1))*F-BETA(J+2)*F1
        F1=F
        F=FF
        DF1=DF
        DF=D
        DDF1=DDF
        DDF=DD
1000  CONTINUE
      END
