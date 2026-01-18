!	To evaluate a polynomial and its derivatives at any point
!
!	N : (input) Degree of the polynomial
!	A : (input) Real array of length N+1 containing the coefficients
!		of the polynomial. A(1) is the constant term and A(N+1) is
!		the coefficient of X**N
!	X : (input) The value of x at which the polynomial is to be evaluated
!	ND : (input) Number of derivatives required. The subroutine always
!		calculates the first derivative irrespective of value of ND
!	PD : (output) Real array of length PD containing the calculated
!		values of the derivatives. PD(I) will contain the Ith derivative
!	The polynomial value is returned through POLD
!
!	Required routines : None

      FUNCTION POLD(N,A,X,ND,PD)
      IMPLICIT REAL*16(A-H,O-Z)
      DIMENSION A(N+1),PD(ND)
 
      POLD=A(N+1)
      DO 2000 I=1,ND
        PD(I)=0.0
2000  CONTINUE
 
      DO 2500 J=N,1,-1
 
        DO 2200 I=ND,2,-1
          PD(I)=PD(I)*X+I*PD(I-1)
2200    CONTINUE
        PD(1)=PD(1)*X+POLD
        POLD=POLD*X+A(J)
2500  CONTINUE
      END
