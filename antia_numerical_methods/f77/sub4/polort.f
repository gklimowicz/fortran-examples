!	Evaluating the orthogonal polynomial basis functions at any value
!	of x using known coefficients
!	Should be used to evaluate the basis using coefficients calculated
!	by POLFIT or POLFIT1.
!
!	M : (input) Degree of polynomial
!	ALP, BETA : (input) Real arrays of length M+1, containing the coefficients
!		required for defining the orthogonal polynomials
!		ALP, BETA could be calculated using POLFIT
!	X : (input) Value of x at which polynomials needs to be evaluated
!	F : (output) Real array of length M+1 containing the value of
!		orthogonal polynomials at X
!	DF : (output) Real array of length M+1 containing first
!		derivative of F at X
!	DDF : (output) Real array of length M+1 containing second
!		derivative of F at X
!	
!	Required routines : None

 
      SUBROUTINE POLORT(M,ALP,BETA,X,F,DF,DDF)
!      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION ALP(M+1),BETA(M+1),F(M+1),DF(M+1),DDF(M+1)

      F(1)=1.0
      F(2)=(X-ALP(1))
      DF(1)=0.0
      DF(2)=1.0
      DDF(1)=0.0
      DDF(2)=0.0
 
!	The recurrence relations
      DO 1000 J=3,M+1
        DDF(J)=2.*DF(J-1)+(X-ALP(J-1))*DDF(J-1)-BETA(J-1)*DDF(J-2)
        DF(J)=F(J-1)+(X-ALP(J-1))*DF(J-1)-BETA(J-1)*DF(J-2)
        F(J)=(X-ALP(J-1))*F(J-1)-BETA(J-1)*F(J-2)
1000  CONTINUE
      END
