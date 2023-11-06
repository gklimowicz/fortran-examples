!	To compute integral of a tabulated function using cubic splines
!
!	A : (input) Lower limit of integration
!	B : (input) Upper limit of integration (B > A)
!	SINT : (output) Computed value of integral using cubic splines
!	TINT : (output) Computed value of integral using trapezoidal rule
!	N : (input) Number of tabular points
!	X : (input) Array of length N containing abscissas (in ascending order)
!	F : (input) Array of length N containing the function values
!	C : (input) Array of length 3*N containing the coefficients of splines
!	IER : (output) Error parameter, IER=0 implies successful execution
!		IER=31 implies lower limit A is outside the table
!		IER=32 implies upper limit B is outside the table
!		IER=301 implies A>B or X(1)>X(N) and no calculations are done
!		Other values may be set by SPLEVL
!
!	Required routines : SPLEVL
!
      SUBROUTINE SPLINT(A,B,SINT,TINT,N,X,F,C,IER)
!      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(N),F(N),C(3,N)

      IER=0
      SINT=0.0
      TINT=0.0
      IF(A.EQ.B) RETURN
      IER=301
      IF(A.GT.B.OR.X(1).GE.X(N)) RETURN

      A1=A
!	Evaluate the function at x=A for trapezoidal rule
      F1=SPLEVL(A1,N,X,F,C,DFB,DDFB,IER)
      IF(IER.GT.100) RETURN
      IER=0
      IF(A.LT.X(1).OR.A.GT.X(N)) IER=31
      IF(B.LT.X(1).OR.B.GT.X(N)) IER=32

!	Integrating over the N-1 subintervals
      DO 2000 I=1,N-1
        IF(A1.LT.X(I+1).OR.I.EQ.N-1) THEN
          B1=MIN(B,X(I+1))
          IF(I.EQ.N-1) B1=B
          D1=A1-X(I)
          D2=B1-X(I)
!	Add integral of cubic spline over [A1,B1]
          SINT=SINT+F(I)*(D2-D1)+C(1,I)*(D2*D2-D1*D1)/2.+
     1         C(2,I)*(D2**3-D1**3)/3.+C(3,I)*(D2**4-D1**4)/4.
          F2=F(I+1)
          IF(B1.NE.X(I+1)) F2=SPLEVL(B1,N,X,F,C,DFB,DDFB,IER1)
!	Trapezoidal rule approximation to integral
          TINT=TINT+0.5*(B1-A1)*(F1+F2)
          IF(B.LE.X(I+1)) RETURN
          A1=B1
          F1=F2
        ENDIF
2000  CONTINUE
      END
